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
 * matrix during an epoch, 
 * then a Barrier, then another epoch 
 * where all the ranks perform Get (on
 * their patch)
 * Some of the MPI functions are not defined
 * in El, hence this test mixes MPI routines
 * and MPI from El. This is nasty, but at one
 * point would be made better.
 * This also requires MPI-3, as we have used
 * MPI-3 fetch and op to simulate a global
 * counter. We could later use MPI-2 version of
 * this function if necessary.
 */
#include "El.hpp"
#include <cassert>
using namespace El;

#define ITER 		1
//#define DIM 		1000
//#define AXPY_DIM 	100
//#define DIM 		20
//#define AXPY_DIM 	4
#define DIM 		8
#define AXPY_DIM 	2

#define ALPHA		2.0
#define FOP_ROOT 	0

#if MPI_VERSION < 3
#  error SORRY, THE TEST ONLY WORKS WITH MPI VERSION > 3
#endif

long ReadInc (MPI_Win win, MPI_Aint offset, long inc)
{
    long otemp;			
    MPI_Fetch_and_op (&inc, &otemp, MPI_LONG, FOP_ROOT, offset, MPI_SUM,
	    win);
    MPI_Win_flush (FOP_ROOT, win);

    return otemp;
}

int main (int argc, char *argv[])
{
    Initialize (argc, argv);
    mpi::Comm comm = mpi::COMM_WORLD;
    MPI_Win win;
    const Int commRank = mpi::Rank (comm);
    const Int commSize = mpi::Size (comm);
    double t1, t2, seconds;
    void *win_base;
    long counter, next = 0;

    assert (DIM % AXPY_DIM == 0);

    try
    {
	// Initialization
	// Allocate memory and create window for ReadInc
	MPI_Win_allocate (sizeof (long), sizeof (long), MPI_INFO_NULL,
		comm.comm, &win_base, &win);
	memset (win_base, 0, sizeof (long));
	MPI_Win_lock_all (MPI_MODE_NOCHECK, win);

	// Create window
	Grid grid (comm);

	// Create an DIM X DIM distributed matrix over the given grid
	DistMatrix < double, MC, MR > A (DIM, DIM, grid);

	// Set every entry of A to zero
	Zeros (A, DIM, DIM);

	// Print the original A
	if (DIM <= 20)
	    Print (A, "Original distributed A");

	t1 = MPI_Wtime();
	for (Int k = 0; k < ITER; ++k)
	{
	    if (commRank == 0)
		std::cout << "Iteration " << k << std::endl;

	    // Open up a LOCAL_TO_GLOBAL interface to A
	    AxpyInterface < double >interface;
	    interface.Attach (LOCAL_TO_GLOBAL, A);

	    Matrix < double >B (AXPY_DIM, AXPY_DIM);
	    Identity (B, AXPY_DIM, AXPY_DIM);
	    // AXPY into parts of the DistMatrix
	    counter = ReadInc (win, 0, (long) 1);
	    for (int i = 0; i < DIM; i += AXPY_DIM)
	    {
		if (counter == next)
		{
		    for (int j = 0; j < DIM; j += AXPY_DIM)
		    {
			interface.Axpy (ALPHA, B, i, j);
#if DEBUG > 2
			std::cout << std::to_string(commRank) + ": AXPY patch: " 
			    + std::to_string(i) + "," + std::to_string(j) 
			    << std::endl;
#endif
		    }
		    counter = ReadInc (win, 0, (long) 1);
		}
		next++;
	    }

	    interface.Detach ();

#if DEBUG > 1
	    for (int j = 0; j < commSize; j++)
	    {
		if (j == commRank)
		{
		    if (DIM <= 20)
			Print (A, "Updated distributed A");
		}
	    }
#endif
	    // Reattach to A, but in the GLOBAL_TO_LOCAL direction
	    interface.Attach (GLOBAL_TO_LOCAL, A);
	    Matrix < double >C;
	    Zeros (C, AXPY_DIM, AXPY_DIM);
	    // Bring my updated patch to me from DistMatrix
	    for (int i = 0; i < DIM; i += AXPY_DIM)
	    {
		if (counter == next)
		{
		    for (int j = 0; j < DIM; j += AXPY_DIM)
		    {
			interface.Axpy (1.0, C, i, j);
#if DEBUG > 2
			std::cout << std::to_string(commRank) + ": GET patch: " 
			    + std::to_string(i) + "," + std::to_string(j) 
			    << std::endl;
#endif
		    }
		    counter = ReadInc (win, 0, (long) 1);
		}
		next++;
	    }
	    interface.Detach ();
#if DEBUG > 1
	    if (DIM <= 20 && commSize < 16)
	    {
		for (int j = 0; j < commSize; j++)
		{
		    if (j == commRank)
		    {
			Print (A, "Updated distributed A");
		    }
		}
		mpi::Barrier ( comm );
		for (int j = 0; j < commSize; j++)
		{
		    if (j == commRank)
		    {
			// Process 0 can now locally print its copy of A
			Print (C, "Patch of A");
		    }
		}
		mpi::Barrier ( comm );
	    }
	    else
	    {
		if ( commRank == 0 && k == (ITER-1) )
		    std::cout << "Inifinity norm of local matrix after " 
			<< k+1 << " iterations: " << InfinityNorm ( C ) << "\n";
	    }
#endif
	}
	t2 = MPI_Wtime();
	seconds = (t2 - t1); ///ITER;
	double total_secs;

	MPI_Reduce(&seconds, &total_secs, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (commRank == 0)
	    printf("Time taken (secs):%lf \n", total_secs);
	
	// clear window object for FOP
	MPI_Win_unlock_all (win);
	MPI_Win_free (&win);
    }
    catch (std::exception & e)
    {
	ReportException (e);
    }

    mpi::Finalize ();
    return 0;
}
