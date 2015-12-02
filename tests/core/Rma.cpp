/*
   Copyright (c) 2009-2015, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include <cassert>
using namespace El;

void calc_indices (int dims[], int block[], int map[], int size, int rank, int hi[], int lo[])
{
    int right_lo[2] = { -1, -1 };
    int bottom_lo[2] = { -1, -1 };

    const int nprow = block[0];
    const int npcol = block[1];

    assert (size == (nprow * npcol));

    // 2d process grid
    int * grid = new int [nprow * npcol];
    for (int j = 0; j < npcol; j++)
	for (int i = 0; i < nprow; i++)
	    grid[i + j*nprow] = i + j*nprow;

    mpi::Barrier (mpi::COMM_WORLD);

    // calculate hi/lo
    bool isbreak = false;
    for (int i = 0; i < nprow; i++)
    {
	for (int j = nprow; j < (nprow + npcol); j++)
	{
	    int jj = j - nprow;
	    int crank = grid[i * npcol + jj];
	    if (crank == rank)
	    {
		// lo     
		lo[0] = map[i];
		lo[1] = map[j];
		// bottom neighbor
		if (j < ((npcol + nprow) - 1))
		{
		    bottom_lo[0] = map[i];
		    bottom_lo[1] = map[j + 1];
		}
		// right neighbor
		if (i < (nprow - 1))
		{
		    right_lo[0] = map[i + 1];
		    right_lo[1] = map[j];
		}
		// hi           
		// get width -- hi[0] from right neighbor
		if (right_lo[0] == -1 && right_lo[1] == -1)
		    hi[0] = dims[0] - 1;
		else
		{
		    if (nprow == 1)
			hi[0] = right_lo[0];
		    else
			hi[0] = right_lo[0] - 1;
		}
		// get height from bottom neighbor
		if (bottom_lo[0] == -1 && bottom_lo[1] == -1)
		    hi[1] = dims[1] - 1;
		else
		{
		    if (npcol == 1)
			hi[1] = bottom_lo[1];
		    else
			hi[1] = bottom_lo[1] - 1;
		}

		isbreak = true;
		break;
	    }
	}
	if (isbreak)
	    break;
    }
}

    int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    int dims[2] = { 8, 8 };
    int block[2] = { 2, 2 };
    int map[4] = { 0, 5, 0, 5 };

    int hi[2], lo[2];

    try 
    {
	const Int m = dims[1];
	const Int n = dims[0];

	Grid g( comm );

	// hi/lo calculation
	calc_indices (dims, block, map, commSize, commRank, hi, lo);
	DistMatrix<double> A(m, n, true, g);
	//DistMatrix<double> A(m, n, g);
	RmaInterface<double> interface;

	interface.Attach( A );

	const Int height = hi[1] - lo[1] + 1;
	const Int width = hi[0] - lo[0] + 1;
	Matrix<double> X( height, width );	

	for( Int j=0; j<X.Width(); ++j )
	    for( Int i=0; i<X.Height(); ++i )
		X.Set(i,j,1);

	const Int i = lo[1];
	const Int j = lo[0];

	// prints
	for (int i = 0; i < commSize; i++)
	{
	    if (i == commRank)
	    {
		cout << "rank #" << i << ": lo{0,1}: " << lo[0] << ", " << lo[1]
		    << ": hi{0,1}: " << hi[0] << ", " << hi[1] 
		    << " X height/width: " << X.Height() << "/" << X.Width() << endl;
	    }
	    mpi::Barrier( mpi::COMM_WORLD );
	}

	for (int i = 0; i < commSize; i++)
	{
	    if (i == commRank)
	    {
		std::cout << "# " << i;
		Print( X, "X" );
	    }
	    mpi::Barrier( mpi::COMM_WORLD );
	}

	// start communication
	for( Int k=0; k<10; ++k )
	{
	    if( commRank == 0 )
		std::cout << "Iteration " << k << std::endl;

	    interface.Acc( X, i, j );

	    interface.Flush( X );	
	}

	mpi::Barrier( mpi::COMM_WORLD );
	Print( A, "A" );
	mpi::Barrier( mpi::COMM_WORLD );

	// detach
	interface.Detach();
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
