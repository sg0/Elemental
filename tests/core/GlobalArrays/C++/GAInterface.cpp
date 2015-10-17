/*
   Copyright (c) 2009-2014, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   Copyright (c) 2014, Sayan Ghosh, University of Houston
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
http://opensource.org/licenses/BSD-2-Clause
*/
/*
 * This test approximates a Hartree-Fock
 * application, all the ranks perform Acc 
 * (or Axpy) on different patches of the 
 * matrix during an epoch, then Flush all, 
 * then a Barrier, then another epoch 
 * where all the ranks perform Get (on
 * their patch)
 * Some of the MPI functions are not defined
 * in El, hence this test mixes MPI routines
 * and MPI from El. This is nasty, but at one
 * point would be made better.
 * This is implemented using MPI-3
 */
#include "El.hpp"
#include <cassert>
using namespace El;

#define ITER 		10
#define DIM 		20
#define AXPY_DIM 	4

//#define DIM 		1024
//#define AXPY_DIM 	256
//#define DIM 		8192
//#define AXPY_DIM 	2048
//#define DIM 		16384
//#define AXPY_DIM 	512

#define FOP_ROOT 	0

#if MPI_VERSION < 3
#  error SORRY, THE TEST ONLY WORKS WITH MPI VERSION > 3
#endif

Int main (Int argc, char *argv[])
{
    Initialize (argc, argv);
    mpi::Comm comm = mpi::COMM_WORLD;
    mpi::Window win;
    const Int commRank = mpi::Rank (comm);
    const Int commSize = mpi::Size (comm);
    double t1, t2, seconds;
    void *win_base;
    long counter = 0, next = 0;

    assert (DIM % AXPY_DIM == 0);

    try
    {
	t1 = MPI_Wtime();
	
	GlobalArrays < double > gaint;
	// GA parameters
	Int ndim = 2;
	Int dims[2] = {DIM, DIM};
	Int ld[2] = {AXPY_DIM, AXPY_DIM};
	Int lo[2] = {-1, -1};
	Int hi[2] = {-2, -2};

	// nb handle
	Int nbhdl = -1;
	double scale = 1.0;
	
	// Initialize Global Array
	gaint.GA_Initialize ();

	// Create a Global Array
	Int g_a = gaint.GA_Create (ndim, dims, "Elemental Global Arrays");
	// Global Array for Read increment
	Int len = 1;
	Int g_cnt = gaint.GA_Create (1, &len, "Counter Global Arrays");
	// access local portions of global array
	int loclo[2], lochi[2]; // local low and hi
	int locld[1];
	// find my distribution
	gaint.NGA_Distribution( g_a, commRank, loclo, lochi );
        // access local patch
	double *ptr;
	gaint.NGA_Access (g_a, loclo, lochi, &ptr, locld);
	const int height = lochi[0] - loclo[0] + 1;
	const int width = lochi[1] - loclo[1] + 1;
	const int ldim = *locld; // only 1st dimension matters
	// update local portion of global matrix
	for (int j = 0; j < width; ++j) {
	    for (int i = 0; i < height; ++i) {
		ptr[i + j*ldim] = 2.0;
	    }
	}
	
	for (Int k = 0; k < ITER; ++k)
	{
	    if (commRank == 0)
		std::cout << "Iteration " << k << std::endl;

	    // local matrix patch that would be sent to the GA
	    Matrix < double >B (AXPY_DIM, AXPY_DIM);
	    Identity (B, AXPY_DIM, AXPY_DIM);
    
	    counter = 0; 
	    next = 0;
	    
	    // AXPY into parts of the DistMatrix
	    int zero = 0;
	    next = gaint.NGA_Read_inc (g_cnt, &zero, long (1));
	    for (Int i = 0; i < DIM; i += AXPY_DIM)
	    {
		if (counter == next)
		{
		    for (Int j = 0; j < DIM; j += AXPY_DIM)
		    {
			// wait for pending transfers before firing accumulates
			gaint.NGA_NbWait (&nbhdl);
			
			// calculate indices of the patch
			lo[0] = i;
			lo[1] = j;
			hi[0] = lo[0] + AXPY_DIM;
			hi[1] = lo[1] + AXPY_DIM;

			hi[0] = hi[0]-1;
			hi[1] = hi[1]-1;

			// nb accumulate
			gaint.NGA_NbAcc (g_a, lo, hi, B.Buffer(), ld, (void *)&scale, &nbhdl);
			//gaint.NGA_Acc (g_a, lo, hi, B.Buffer(), ld, &scale);
#if DEBUG > 2
			std::cout << std::to_string(commRank) + ": AXPY patch: " 
			    + std::to_string(i) + "," + std::to_string(j) 
			    << std::endl;
#endif
		    }
	            next = gaint.NGA_Read_inc (g_cnt, &zero, long (1));
		}
		counter++;
	    }
	    
	    // sync is flushall and barrier
	    gaint.GA_Sync ();

	    counter = 0;
	    next = 0;
	    
	    // Bring my updated patch to me from DistMatrix
	    Matrix < double >C;
	    Zeros (C, AXPY_DIM, AXPY_DIM);
	    
	    next = gaint.NGA_Read_inc (g_cnt, &zero, long (1));
	    for (Int i = 0; i < DIM; i += AXPY_DIM)
	    {
		if (counter == next)
		{
		    for (Int j = 0; j < DIM; j += AXPY_DIM)
		    {
			// calculate indices of the patch
			lo[0] = i;
			lo[1] = j;
			hi[0] = lo[0] + AXPY_DIM;
			hi[1] = lo[1] + AXPY_DIM;

			hi[0] = hi[0]-1;
			hi[1] = hi[1]-1;

			// get patches from GA
			gaint.NGA_Get (g_a, lo, hi, C.Buffer(), ld); 
#if DEBUG > 2
			std::cout << std::to_string(commRank) + ": GET patch: " 
			    + std::to_string(i) + "," + std::to_string(j) 
			    << std::endl;
#endif
		    }
		    next = gaint.NGA_Read_inc (g_cnt, &zero, long (1));
		}
		counter++;
	    }
	    
	    gaint.GA_Sync ();
		    
#if DEBUG > 1
	    for (Int j = 0; j < commSize; j++)
	    {
		if (j == commRank)
		{
		    if (DIM <= 20)
			gaint.GA_Print (g_a);
		}
	    }
	    mpi::Barrier ( comm );
	    for (Int j = 0; j < commSize; j++)
	    {
		if (j == commRank)
		{
		    // Process 0 can now locally print its copy of A
		    if (DIM <= 20)
			Print (C, "Patch of A");
		}
	    }
	    mpi::Barrier ( comm );
#endif
	}
	    
	// Collectively terminate GA
	gaint.GA_Terminate ();    

	t2 = MPI_Wtime();
	seconds = (t2 - t1)/ITER;
	double total_secs;

	MPI_Reduce(&seconds, &total_secs, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (commRank == 0)
	    printf("Time taken for AXPY (secs):%lf \n", total_secs);
    }
    catch (std::exception & e)
    {
	ReportException (e);
    }

    
    mpi::Finalize ();
    return 0;
}
