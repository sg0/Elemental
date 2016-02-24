/*
   Copyright (c) 2009-2014, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   Copyright (c) 2014, Jeff Hammond (Intel)
   Copyright (c) 2014, Sayan Ghosh (Washington State University)
   All rights reserved.

Authors:
Sayan Ghosh
Washington State University

This file contains implementation of few functions from Global Arrays
C API (http://hpc.pnl.gov/globalarrays/api/capi.shtml) using Elemental
RMAInterface (which uses MPI-3 one sided API). The data 
distribution mechanism of GA and Elemental are completely different, 
hence some GA functions which assumes contiguous data patches are 
deliberately not implemented. Only 1/2-dimensional arrays are supported.

This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause

This is a subset of Global Arrays toolkit, the intention is to use
as many core Elemental functions as possible in building this. Only for
remote put/gets El::RmaInterface should be used.
*/
#include "El.hpp"

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY) && defined(EL_ENABLE_RMA_GLOBAL_ARRAYS)
namespace El
{
// constructors    
template<typename T>
GlobalArrays< T >::GlobalArrays()
    : ga_initialized( false ),
      ga_handles( 0 ),
      matrices_( 0 ),
      nbhdls_( 0 )
{}

// destructor
template<typename T>
GlobalArrays< T >::~GlobalArrays()
{
    if( std::uncaught_exception() )
    {
	std::ostringstream os;
	os << "Uncaught exception detected during GlobalArrays destructor "
	    "that required a call to GA_Terminate. Instead of allowing for the "
	    "possibility of GA_Terminate throwing another exception and "
	    "resulting in a 'terminate', we instead immediately dump the "
	    "call stack (if not in RELEASE mode) since the program will "
	    "likely hang:" << std::endl;
	std::cerr << os.str();
	DEBUG_ONLY( DumpCallStack() )
    }
    else
	GA_Terminate();
}

// At present, up to 2D GA is supported, if ndim is 1
// then it is assumed the resulting GA will *only*
// be used for Read_increment operation
// and in that case, RMAInterface object is instantiated
// to nullptr, so is the DistMatrix object -- in other
// words, no DistMatrix is created for a 1D GA
// chunk -- minimum blocking size for each dimension
template<typename T>
Int GlobalArrays< T >::GA_Create(Int ndim, Int dims[], const char *array_name, Int chunk[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Create" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array"); 
	
    const Int handle = ga_handles.size();

    // create a GA instance and push
    // it into the ga_handles vector
    ga_handles.push_back( GA() );

    // default grid
    const Grid& grid = DefaultGrid();	
    const Int p = grid.Size();
    const Int rank = grid.Rank();

    // 2-D array, create DistMatrix
    if (ndim == 2)
    {
	ga_handles[handle].ndim = 2;
	ga_handles[handle].length = (dims[0] * dims[1]);
	
	// create distmatrix over default grid, i.e mpi::COMM_WORLD
	// using MC x MR distribution
	// this won't be allocated until RmaInterface->Attach
	// dim[0] = height, dim[1] = width
#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
	ga_handles[handle].DM = new DistMatrix< T, MC, MR  >( dims[1], dims[0], true, grid );
#else
	ga_handles[handle].DM = new DistMatrix< T, MC, MR  >( dims[1], dims[0], grid );
#endif
	DistMatrix< T, MC, MR > &D = *(ga_handles[handle].DM);
	
	//FIXME minimize metadata size
	// resize Int vectors for storing local heights and widths
	ga_handles[handle].ga_local_height.resize( p );
	ga_handles[handle].ga_local_width.resize( p );
	// resize Int vectors for storing local lo/hi
	ga_handles[handle].ga_hi.resize( 2*p );
	ga_handles[handle].ga_lo.resize( 2*p );

	const Int localHeight = D.LocalHeight();
	const Int localWidth = D.LocalWidth();
	
        const Int my_rank = D.DistRank();  
	Int lo[2];
	Int hi[2];
	
	// GA is distributed evenly
	// vertical strips
	if (chunk == nullptr)
	{ 
	    ga_handles[handle].patchHeight = dims[1];

	    Int strips = (dims[0] / p);
	    const Int rem = (dims[0] % p);
	    if ( rank == (p - 1) )
		strips += rem;
	
	    ga_handles[handle].patchWidth = strips;
	}
	else if (chunk[0] == dims[0]) // horizontal strips
	{
	    ga_handles[handle].patchWidth = dims[0];

	    Int strips = (dims[1] / p);
	    const Int rem = (dims[1] % p);
	    if ( rank == (p - 1) )
		strips += rem;

	    ga_handles[handle].patchHeight = strips;
	}
	else if (chunk[1] == dims[1]) // vertical strips
	{
	    ga_handles[handle].patchHeight = dims[1];

	    Int strips = (dims[0] / p);
	    const Int rem = (dims[0] % p);
	    if ( rank == (p - 1) )
		strips += rem;
	
	    ga_handles[handle].patchWidth = strips;
	}
	else
	{
	    ga_handles[handle].patchWidth = chunk[0];
	    ga_handles[handle].patchHeight = chunk[1];

	    Int strips_h = (dims[1] / p);
	    Int strips_w = (dims[0] / p);

	    const Int rem_h = (dims[1] % p);
	    const Int rem_w = (dims[0] % p);

	    if ( rank == (p - 1) )
	    {
		strips_h += rem_h;
		strips_w += rem_w;
	    }

	    if (chunk[1] != strips_h)
		ga_handles[handle].patchHeight = strips_h;

	    if (chunk[0] != strips_w)
		ga_handles[handle].patchWidth = strips_w;
	}

	// lo/hi
	lo[1] = 0; 
	lo[0] = rank * ga_handles[handle].patchHeight;
	hi[0] = lo[0] + ga_handles[handle].patchWidth;
	hi[1] = lo[1] + ga_handles[handle].patchHeight;
	hi[0] -= 1;
	hi[1] -= 1;

	// allocate local matrix for use during access/release
	const Int currentIndex = matrices_.size();
	matrices_.push_back( matrix_params_() );

	matrices_[currentIndex].ga_index_ = handle;
	matrices_[currentIndex].is_accumulate_ = false;

	matrices_[currentIndex].M_ 
	    = new Matrix< T >( ga_handles[handle].patchHeight, 
		    ga_handles[handle].patchWidth );
	
	const Int pos = rank * 2;
	ga_handles[handle].ga_lo[pos]     = lo[0];
	ga_handles[handle].ga_lo[pos + 1] = lo[1];	
	ga_handles[handle].ga_hi[pos]     = hi[0];
	ga_handles[handle].ga_hi[pos + 1] = hi[1];
	// gather every PE's lo/hi
	Int *vlo = ga_handles[handle].ga_lo.data();
	Int *vhi = ga_handles[handle].ga_hi.data();
	mpi::AllGather <Int>( vlo, 2, grid.Comm() );
	mpi::AllGather <Int>( vhi, 2, grid.Comm() );

	// store local heights and widths
	ga_handles[handle].ga_local_height[my_rank] = localHeight;
	ga_handles[handle].ga_local_width[my_rank] = localWidth;
	// NOTE: for MC, MR distribution, DistComm() == VCComm()
	// In GAInterface, this is alright as we only deal
	// with MC, MR distribution
	Int *vheight = ga_handles[handle].ga_local_height.data();
	Int *vwidth = ga_handles[handle].ga_local_width.data();
	mpi::AllGather <Int>( vheight, 1, D.DistComm() );
	mpi::AllGather <Int>( vwidth, 1, D.DistComm() );

	// call rmainterface constructor
	ga_handles[handle].rmaint = new RmaInterface< T >();
	// attach DM to RMAInterface
	ga_handles[handle].rmaint->Attach( D );
    }
    else if (ndim == 1) // fetch-and-op
    {
	ga_handles[handle].ndim = 1;
	// length of GA
	ga_handles[handle].length = *(dims);
	// start access epoch on FOP window
	mpi::WindowAllocate<T>( ga_handles[handle].length, grid.Comm(), 
		ga_handles[handle].fop_win );
	mpi::WindowLock( ga_handles[handle].fop_win );
    }
    else
	LogicError ("Up to 2 dimensions supported presently");

    ga_handles[handle].is_destroyed = false;
    mpi::Barrier( grid.Comm() );

    return handle;
}

// mimic GA irregular distribution
template<typename T>
Int GlobalArrays< T >::GA_Create_irreg(Int ndim, Int dims[], const char *array_name, Int block[], Int map[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Create_irreg" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array"); 
	
    const Int handle = ga_handles.size();

    // create a GA instance and push
    // it into the ga_handles vector
    ga_handles.push_back( GA() );

    // default grid
    const Grid& grid = DefaultGrid();	
    const Int p = grid.Size();

    // 2-D array, create DistMatrix
    if (ndim == 2)
    {
	const Int nprow = block[0];
	const Int npcol = block[1];

	if (p != (nprow * npcol))
	    LogicError( "Incorrect block parameter specified for irregular GA distribution" );
	
	ga_handles[handle].ndim = 2;
	ga_handles[handle].length = (dims[0] * dims[1]);
	
	// create distmatrix over default grid, i.e mpi::COMM_WORLD
	// using MC x MR distribution
	// this won't be allocated until RmaInterface->Attach
	// dim[0] = width, dim[1] = height
#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
	ga_handles[handle].DM = new DistMatrix< T, MC, MR  >( dims[1], dims[0], true, grid );
#else
	ga_handles[handle].DM = new DistMatrix< T, MC, MR  >( dims[1], dims[0], grid );
#endif
	DistMatrix< T, MC, MR > &D = *(ga_handles[handle].DM);
	
	// resize Int vectors for storing local heights and widths
	ga_handles[handle].ga_local_height.resize( p );
	ga_handles[handle].ga_local_width.resize( p );
	// resize Int vectors for storing local lo/hi
	ga_handles[handle].ga_hi.resize( 2*p );
	ga_handles[handle].ga_lo.resize( 2*p );

	const Int localHeight = D.LocalHeight();
	const Int localWidth = D.LocalWidth();
	
        const Int my_rank = D.DistRank(); 

	Int hi[2]        = { -1, -2 };
	Int lo[2]        = { -1, -1 };
	// To calculate dimensions of blocks,
	// for each block, we need the start
	// indices of the bottom neighbor 
	// and the right neighbor of a 
	// particular block
	Int right_lo[2]  = { -1, -1 };
	Int bottom_lo[2] = { -1, -1 };
  
	// FIXME remove this, and use El::Grid instead,
	// not sure if there is a function that returns
	// a process id when supplied with (i, j)
	// 2d process grid
	Int * grid2D = new Int[p];
	for (Int j = 0; j < npcol; j++)
	    for (Int i = 0; i < nprow; i++)
		grid2D[i + j*nprow] = i + j*nprow;

	// rank in the default grid
	Int rank = grid.Rank();

	// calculate hi/lo
	bool isbreak = false;
	for (int i = 0; i < nprow; i++)
	{
	    for (int j = nprow; j < (nprow + npcol); j++)
	    {
		int jj = j - nprow;
		int crank = grid2D[i * npcol + jj];
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

	ga_handles[handle].patchHeight = hi[1] - lo[1] + 1;
	ga_handles[handle].patchWidth  = hi[0] - lo[0] + 1;

	// allocate local matrix for use during access/release
	const Int currentIndex = matrices_.size();
	matrices_.push_back( matrix_params_() );

	matrices_[currentIndex].ga_index_ = handle;
	matrices_[currentIndex].is_accumulate_ = false;

	matrices_[currentIndex].M_ 
	    = new Matrix< T >( ga_handles[handle].patchHeight, 
		    ga_handles[handle].patchWidth );
	
	// lo/hi
	const Int pos = rank * 2;
	ga_handles[handle].ga_lo[pos]     = lo[0];
	ga_handles[handle].ga_lo[pos + 1] = lo[1];	
	ga_handles[handle].ga_hi[pos]     = hi[0];
	ga_handles[handle].ga_hi[pos + 1] = hi[1];
	// gather every PE's lo/hi
	Int *vlo = ga_handles[handle].ga_lo.data();
	Int *vhi = ga_handles[handle].ga_hi.data();
	mpi::AllGather <Int>( vlo, 2, grid.Comm() );
	mpi::AllGather <Int>( vhi, 2, grid.Comm() );		
	
	// store local heights and widths
	ga_handles[handle].ga_local_height[my_rank] = localHeight;
	ga_handles[handle].ga_local_width[my_rank] = localWidth;
	// NOTE: for MC, MR distribution, DistComm() == VCComm()
	// In GAInterface, this is alright as we only deal
	// with MC, MR distribution
	Int *vheight = ga_handles[handle].ga_local_height.data();
	Int *vwidth = ga_handles[handle].ga_local_width.data();
	mpi::AllGather <Int>( vheight, 1, D.DistComm() );
	mpi::AllGather <Int>( vwidth, 1, D.DistComm() );

	// call rmainterface constructor
	ga_handles[handle].rmaint = new RmaInterface< T >();
	// attach DM to RMAInterface
	ga_handles[handle].rmaint->Attach( D );
	
	delete[] grid2D;
    }
    else if (ndim == 1) // fetch-and-op
    {
	ga_handles[handle].ndim = 1;
	// length of GA
	ga_handles[handle].length = *(dims);
	// start access epoch on FOP window
	mpi::WindowAllocate<T>( ga_handles[handle].length, grid.Comm(), 
		ga_handles[handle].fop_win );
	mpi::WindowLock( ga_handles[handle].fop_win );
    }
    else
	LogicError ("Up to 2 dimensions supported presently");
	
    ga_handles[handle].is_destroyed = false;
    mpi::Barrier( grid.Comm() );

    return handle;
}

template<typename T>
Int GlobalArrays< T >::GA_Duplicate(Int g_a, const char *array_name)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Duplicate" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");   
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible, as Global Arrays has been destroyed");
	    )

    const Int handle = ga_handles.size();
	
    const Grid& grid = DefaultGrid();

    // create a GA instance and push
    // it into the ga_handles vector
    ga_handles.push_back( GA() );
    
    const Int ndim = ga_handles[g_a].ndim;
    
    if (ndim == 2)
    {
	ga_handles[handle].ndim = 2;
	ga_handles[handle].length = ga_handles[g_a].length;
	
	DistMatrix< T, MC, MR >& GADM = *(ga_handles[g_a].DM);
	
	// dimensions
	Int dim[2];
	dim[0] = GADM.Height();
	dim[1] = GADM.Width();

	// grid
	const Grid& grid = GADM.Grid();
	const Int p = grid.Size();

	// create distmatrix
#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
	ga_handles[handle].DM = new DistMatrix< T, MC, MR >( dim[0], dim[1], true, grid );
#else
	ga_handles[handle].DM = new DistMatrix< T, MC, MR >( dim[0], dim[1], grid );
#endif
	DistMatrix< T, MC, MR > &D = *(ga_handles[handle].DM);
	
	// ideal patch width and height
	ga_handles[handle].patchHeight = ga_handles[g_a].patchHeight;
	ga_handles[handle].patchWidth = ga_handles[g_a].patchWidth;

	// allocate local matrix for use during access/release
	const Int currentIndex = matrices_.size();
	matrices_.push_back( matrix_params_() );

	matrices_[currentIndex].ga_index_ = handle;
	matrices_[currentIndex].is_accumulate_ = false;

	matrices_[currentIndex].M_ 
	    = new Matrix< T >( ga_handles[handle].patchHeight, 
		    ga_handles[handle].patchWidth );
	
	// resize Int vectors for storing local heights and widths
	ga_handles[handle].ga_local_height.resize( p );
	ga_handles[handle].ga_local_width.resize( p );		
	// resize Int vectors for storing local lo/hi
	ga_handles[handle].ga_hi.resize( 2*p );
	ga_handles[handle].ga_lo.resize( 2*p );

	// locally copy 
	MemCopy<Int> (ga_handles[handle].ga_local_height.data(),
		ga_handles[g_a].ga_local_height.data(),
		p);

	MemCopy<Int> (ga_handles[handle].ga_local_width.data(), 
		ga_handles[g_a].ga_local_width.data(), 
		p);
	
	MemCopy<Int> (ga_handles[handle].ga_lo.data(), 
		ga_handles[g_a].ga_lo.data(), 
		2 * p);

	MemCopy<Int> (ga_handles[handle].ga_hi.data(), 
		ga_handles[g_a].ga_hi.data(),
		2 * p);

	// call rmainterface constructor
	ga_handles[handle].rmaint = new RmaInterface< T >();
	// attach DM to RMAInterface
	ga_handles[handle].rmaint->Attach( D );
    }
    else if (ndim == 1)// fetch-and-op
    {
	ga_handles[handle].ndim = 1;
	
	// length of GA
	ga_handles[handle].length = ga_handles[g_a].length;
	const Int bufferSize = ga_handles[handle].length * sizeof(T);

	// start access epoch on FOP window
	mpi::WindowAllocate<T>( ga_handles[handle].length, grid.Comm(), 
		ga_handles[handle].fop_win );
	mpi::WindowLock( ga_handles[handle].fop_win );    
    }
    else
	LogicError ("Up to 2 dimensions supported presently");
	
    ga_handles[handle].is_destroyed = false;
    mpi::Barrier( grid.VCComm() );
    
    return handle;
}

#if 0
// GA_Add using EntrywiseMap + UpdateSubmatrix
// g_c = alpha * g_a  +  beta * g_b;
template<typename T>
void GlobalArrays< T >::GA_Add(T *alpha, Int g_a, T* beta, Int g_b, Int g_c)
{
    DEBUG_ONLY( 
    	CallStackEntry cse( "GlobalArrays::GA_Add" ) 
    	if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    	if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    	if (g_b < 0 || g_b >= ga_handles.size())
	LogicError ("Invalid GA handle");
    	if (g_c < 0 || g_c >= ga_handles.size())
	LogicError ("Invalid GA handle");
    	if ( ga_handles[g_a].ndim == 1 
	    || ga_handles[g_b].ndim == 1 
	    || ga_handles[g_c].ndim == 1 )
	LogicError ("A 1D GA is not allowed for this operation");
    	if (ga_handles[g_a].is_destroyed
	    || ga_handles[g_b].is_destroyed
	    || ga_handles[g_c].is_destroyed)
	    LogicError ("Operation not possible, as Global Arrays has been destroyed");
    	)

    DistMatrix< T, MC, MR >& GA = *(ga_handles[g_a].DM);
    DistMatrix< T, MC, MR >& GB = *(ga_handles[g_b].DM);
    DistMatrix< T, MC, MR >& GC = *(ga_handles[g_c].DM);
  
    const Int g_a_height = GA.Height();
    const Int g_a_width = GA.Width();
    const Int g_b_height = GB.Height();
    const Int g_b_width = GB.Width();   
    const Int g_c_height = GC.Height();
    const Int g_c_width = GC.Width();
   
    const Grid& grid = GA.Grid();
    
    T a = *alpha;
    T b = *beta;
    T one = T( 1 );
    
    // arrays must have the same shape
    if ((g_a_width != g_b_width) 
	    || (g_a_width != g_c_width) 
	    || (g_b_width != g_c_width))
	LogicError ("Global Arrays of different widths cannot be added. ");
    if ((g_a_height != g_b_height) 
	    || (g_a_height != g_c_height) 
	    || (g_b_height != g_c_height))
	LogicError ("Global Arrays of different heights cannot be added. ");

    // GC = GB
    Copy( GB, GC );
    
    // Entrywise map for b * GB
    if (b != one)
    {
	auto valMap = [b]( T val ) { return (val * b); }; 
        EntrywiseMap( GC, function<T(T)>(valMap) );
    }
    
    // Add a*GA into GC(0,0)
    vector<Int> I(g_b_height), J(g_b_width);
    for( Int i=0; i<g_b_height; ++i ) I[i] = i;
    for( Int j=0; j<g_b_width; ++j ) J[j] = j;

    UpdateSubmatrix( GC, I, J, a, GA );
}
#endif

// g_c = alpha * g_a  +  beta * g_b;
// GA_Add using Gemm
template<typename T>
void GlobalArrays< T >::GA_Add(T *alpha, Int g_a, T* beta, Int g_b, Int g_c)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Add" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (g_b < 0 || g_b >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (g_c < 0 || g_c >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 
		|| ga_handles[g_b].ndim == 1 
		|| ga_handles[g_c].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed
		|| ga_handles[g_b].is_destroyed
		|| ga_handles[g_c].is_destroyed)
	    LogicError ("Operation not possible, as Global Arrays has been destroyed");
	    )

    DistMatrix< T, MC, MR >& GA = *(ga_handles[g_a].DM);
    DistMatrix< T, MC, MR >& GB = *(ga_handles[g_b].DM);
    DistMatrix< T, MC, MR >& GC = *(ga_handles[g_c].DM);
  
    const Int g_a_height = GA.Height();
    const Int g_a_width = GA.Width();
    const Int g_b_height = GB.Height();
    const Int g_b_width = GB.Width();   
    const Int g_c_height = GC.Height();
    const Int g_c_width = GC.Width();

    T a = *((T *)alpha);
    T b = *((T *)beta);
    
    // arrays must have the same shape
    if ((g_a_width != g_b_width) 
	    || (g_a_width != g_c_width) 
	    || (g_b_width != g_c_width))
	LogicError ("Global Arrays of different widths cannot be added. ");
    if ((g_a_height != g_b_height) 
	    || (g_a_height != g_c_height) 
	    || (g_b_height != g_c_height))
	LogicError ("Global Arrays of different heights cannot be added. ");

    DistMatrix< T, MC, MR > Bd( g_a_height, g_a_width, GA.Grid() );
    // Initialize
    Identity( Bd, g_a_height, g_a_width );
    // gc = gb
    Copy( GB, GC );

    // add
    // FIXME dont hardcode
    Int nb = 96;
    GemmAlgorithm alg = GEMM_SUMMA_C;
    
    const Orientation orientA = CharToOrientation( 'N' );
    const Orientation orientB = CharToOrientation( 'N' );
    SetBlocksize( nb );

    // GC = a*GA*Bd + b*GC
    Gemm( orientA, orientB, a, GA, Bd, b, GC, alg);

    Bd.Empty();
}

// copies g_a into g_b, must be of same shape
template<typename T>
void GlobalArrays< T >::GA_Copy(Int g_a, Int g_b)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Copy" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (g_b < 0 || g_b >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 
		|| ga_handles[g_b].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed
		|| ga_handles[g_b].is_destroyed)
	    LogicError ("Operation not possible, as Global Arrays has been destroyed");
	    )

    DistMatrix< T, MC, MR >& GA = *(ga_handles[g_a].DM);
    DistMatrix< T, MC, MR >& GB = *(ga_handles[g_b].DM);

    const Int g_a_height = GA.Height();
    const Int g_a_width = GA.Width();
    const Int g_b_height = GB.Height();
    const Int g_b_width = GB.Width();   

    // arrays must have the same shape
    if (g_a_width != g_b_width)
	LogicError ("Global Arrays of different widths cannot be copied. ");
    if (g_a_height != g_b_height)
	LogicError ("Global Arrays of different heights cannot be copied. ");

    // copy
    Copy( GA, GB );
}

// print GA
// void GA_Print(Int ga)
template<typename T>
void GlobalArrays< T >::GA_Print(Int g_a)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Print" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible, as Global Arrays has been destroyed");
	    )

    DistMatrix< T, MC, MR >& DM = *(ga_handles[g_a].DM);
    
    // TODO should be printed in following format:
    /*
    global array: array_A[1:4,1:2],  handle: -1000
                  1           2  
	      ----------- -----------
	1       2.00000     2.00000
	2       2.00000     2.00000
	3       2.00000     2.00000
	4       2.00000     2.00000
    */

    Print( DM );
    mpi::Barrier( DM.DistComm() );
}

// deallocates ga and frees associated resources
// TODO modify code such that destroyed handles could
// be reused
template<typename T>
void GlobalArrays< T >::GA_Destroy(Int g_a)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Destroy" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (ga_handles[g_a].is_destroyed)
	    return;
	    )

    if (ga_handles[g_a].ndim == 1)
    {
	// clear window object for FOP
	mpi::WindowUnlock( ga_handles[g_a].fop_win );
	mpi::WindowFree ( ga_handles[g_a].fop_win );
    }
    
    if (ga_handles[g_a].ndim == 2)
    {
	// clear vectors that holds local dims
	ga_handles[g_a].ga_local_height.clear();
	ga_handles[g_a].ga_local_width.clear();
	ga_handles[g_a].ga_lo.clear();
	ga_handles[g_a].ga_hi.clear();

	// detach RmaInterface
	// will end access epoch
	ga_handles[g_a].rmaint->Detach();
	delete (ga_handles[g_a].rmaint);

	ga_handles[g_a].DM->Empty();
	delete (ga_handles[g_a].DM);

	// NOTE: Do we need to erase ga entry from global 
	// ga_handles vector?
	// Erasing would mess up g_a handle values, 
	// as GAs could be destroyed at different times
	// so at present vector is cleared only on terminate
	// ga_handles.erase( ga_handles.begin() + g_a );
    }

    // nullify pointers
    ga_handles[g_a].ndim 		= -1;
    ga_handles[g_a].length 		= -1;
    ga_handles[g_a].patchHeight		= -1;
    ga_handles[g_a].patchWidth		= -1;
    ga_handles[g_a].rma_local_pending 	= false;
    ga_handles[g_a].fop_win 		= mpi::WIN_NULL;
    ga_handles[g_a].rmaint 		= nullptr;
    ga_handles[g_a].DM 			= nullptr;

    // delete matrices associated with this GA
    for (typename std::vector< matrix_params_ >::iterator it = matrices_.begin();
	    it != matrices_.end();)
    {
	if ( it->ga_index_ == g_a )
	{
	    it->M_->Empty();
	    delete( it->M_ );
	    it->M_ = nullptr;
	    it = matrices_.erase( it );
	}
	else
	    ++it;
    }

    // delete handles associated with this global array
    for (typename std::vector< nbhdl_t_ >::iterator it = nbhdls_.begin();
	    it != nbhdls_.end();)
    {
	if ( it->nbhandle_ == g_a )
	{
	    it->M_ = nullptr;
	    it = nbhdls_.erase( it );
	}
	else
	    ++it;
    }
    
    ga_handles[g_a].is_destroyed = true;
}

// A := 0.5 * ( A + A' )
template<typename T>
void GlobalArrays< T >::GA_Symmetrize (Int g_a)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Symmetrize" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    DistMatrix< T, MC, MR >& GA = *(ga_handles[g_a].DM);
    
    // create intermediate GAs using g_a parameters
    Int dims[2];
    dims[1] = GA.Width();
    dims[0] = GA.Height();
   
    if (dims[0] != dims[1])
	LogicError ("Cannot symmetrize a non-square matrix");
    
    Int g_at = GA_Create( 2, dims, "transpose array", nullptr );
    Int g_c = GA_Duplicate( g_at, "final array" );

    // transpose
    GA_Transpose( g_a, g_at );

    // add: scale * ( g_a + g_at )
    // g_c = alpha * g_a  +  beta * g_at;
    T scale = T( 0.5 );
    GA_Add( &scale, g_a, &scale, g_at, g_c );
    
    // copy g_c into g_a
    GA_Copy( g_c, g_a );
    
    // destroy the intermediate GAs
    GA_Destroy( g_at );
    GA_Destroy( g_c );
}

// C := alpha*op( A )*op( B ) + beta*C
// where op( X ) = X or X' (transpose)
// A = m x k -- B = k x n -- C = m x n
template<typename T>
void GlobalArrays< T >::GA_Dgemm(char ta, char tb, Int m, Int n, Int k, T alpha, Int g_a, Int g_b, T beta, Int g_c )
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Dgemm" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (g_b < 0 || g_b >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (g_c < 0 || g_c >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 
		|| ga_handles[g_b].ndim == 1 
		|| ga_handles[g_c].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");    
	    if (ga_handles[g_a].is_destroyed
		|| ga_handles[g_b].is_destroyed
		|| ga_handles[g_c].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays have been destroyed");
	    // TODO extra checks pertaining m, n, k
	    )

    DistMatrix< T, MC, MR >& ADM = *(ga_handles[g_a].DM);
    DistMatrix< T, MC, MR >& BDM = *(ga_handles[g_b].DM);
    DistMatrix< T, MC, MR >& CDM = *(ga_handles[g_c].DM);
    
    // set algorithmic block size
    // multiple of 16 (cache-line-size)
    Int nb = 96;
    // keep A stationary
    //GemmAlgorithm alg = GEMM_SUMMA_A;
    //GemmAlgorithm alg = GEMM_SUMMA_B;
    GemmAlgorithm alg = GEMM_SUMMA_C;
    //GemmAlgorithm alg = GEMM_DEFAULT;

    SetBlocksize( nb );

    // Note: We are flipping the matrices because the way it has been
    // distributed, not DistMatrix(height, width), but DM(width, height)
    // in keeping with GA's layout
    const Orientation orientA = CharToOrientation( ((ta == 'T' || ta == 't') ? 'N' : 'T') );
    const Orientation orientB = CharToOrientation( ((tb == 'T' || tb == 't') ? 'N' : 'T') );

    Gemm( orientA, orientB, alpha, ADM, BDM, beta, CDM, alg);
}

// FIXME GA struct should have a grid pointer,
// currently grid information is fetched from
// DM.Grid() or DefaultGrid(), but this is hacky

// fills a global array with a specific value
// collective
template<typename T>
void GlobalArrays< T >::GA_Fill(Int g_a, T* value)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Fill" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");    
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    T a = *value;

    if (ga_handles[g_a].ndim == 1)
    {
	T * fop_base = reinterpret_cast<T *>( mpi::GetWindowBase( ga_handles[g_a].fop_win ) );
	for (Int i = 0; i < ga_handles[g_a].length; i++) fop_base[i] = a;
	// default grid
	const Grid &grid = DefaultGrid();
	mpi::Barrier( grid.Comm() );
    }
    else // fill DM with alpha
    {
	DistMatrix< T, MC, MR >& DM = *(ga_handles[g_a].DM);
	Fill( DM, a );
	mpi::Barrier( DM.DistComm() );
    }
}

// allocate and initialize internal data structures in Global Arrays.
template<typename T>
void GlobalArrays< T >::GA_Initialize()
{
   DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Initialize" ) )

   if (!ga_initialized)
       ga_initialized = true;
}

// single buffer reduce
#define REDUCE(x, n, op) \
    do { \
	const Grid &grid = DefaultGrid(); \
	switch (op) \
	{ \
	    case '+': \
		mpi::AllReduce( x, n, mpi::SUM, grid.Comm() ); \
	        break; \
	    case '*': \
		mpi::AllReduce( x, n, mpi::PROD, grid.Comm() ); \
	        break; \
	    case 'X': \
		mpi::AllReduce( x, n, mpi::MAX, grid.Comm() ); \
	    	break; \
	    case 'N': \
		mpi::AllReduce( x, n, mpi::MIN, grid.Comm() ); \
	    	break; \
	    default: \
		LogicError ("Unsupported global operation specified"); \
	} \
    } while (0)

// NOTE: transfers -- BXFER is local+remote
// blocking in keeping with the GA model 
// LNBXFER is locally nonblocking, needs
// NbWait to enforce local completion
// NBXFER is locally blocking
// (does not apply to Get, which ensures 
// local+remote completion anyway)

// locally nonblocking
#define LNBXFER(type, g_a, M, i, j) \
    do { \
	switch (type) \
	{ \
	    case 'A': \
		ga_handles[g_a].rmaint->Iacc (M, i, j); \
	        break; \
	    case 'P': \
		ga_handles[g_a].rmaint->Iput (M, i, j); \
	        break; \
	    case 'G': \
		ga_handles[g_a].rmaint->Iget (M, i, j); \
	    	break; \
	    default: \
		LogicError ("Unsupported nonblocking transfer type"); \
	    ga_handles[g_a].rma_local_pending = true; \
	} \
    } while (0)

// blocking transfers -- blocks the calling process
// till transfer is initiated (locally complete)
#define NBXFER(type, g_a, M, i, j) \
    do { \
	switch (type) \
	{ \
	    case 'A': \
		ga_handles[g_a].rmaint->Acc (M, i, j); \
	        break; \
	    case 'P': \
		ga_handles[g_a].rmaint->Put (M, i, j); \
	        break; \
	    case 'G': \
		ga_handles[g_a].rmaint->Get (M, i, j); \
	    	break; \
	    default: \
		LogicError ("Unsupported nonblocking transfer type"); \
	} \
    } while (0)

// global commutative operations
template<typename T>
void GlobalArrays< T >::GA_Gop(T x[], Int n, char op)
{
   DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Gop" ) )
   if (!ga_initialized)
       LogicError ("Global Arrays must be initialized before any operations");
   // op
   REDUCE( x, n, op );
}

// dot product
template<typename T>
T GlobalArrays< T >::GA_Dot(Int g_a, Int g_b)
{
   DEBUG_ONLY( 
	   CallStackEntry cse( "GlobalArrays::GA_Dot" ) 
	   if (!ga_initialized)
	   LogicError ("Global Arrays must be initialized before any operations");
	   if (g_a < 0 || g_a >= ga_handles.size())
	   LogicError ("Invalid GA handle");
	   if (g_b < 0 || g_b >= ga_handles.size())
	   LogicError ("Invalid GA handle");
	   if ( ga_handles[g_a].ndim == 1 
	       || ga_handles[g_b].ndim == 1 )
	   LogicError ("A 1D GA is not allowed for this operation"); 
	   if (ga_handles[g_a].is_destroyed
	       || ga_handles[g_b].is_destroyed)
	   LogicError ("Operation not possible as Global Arrays has been destroyed");
	   )

   const DistMatrix< T, MC, MR >& A = *(ga_handles[g_a].DM);
   const DistMatrix< T, MC, MR >& B = *(ga_handles[g_b].DM);

   T result = Dot( A, B );

   return result;
}

// barrier
// TODO dont sync over comm world
// fetch comm from ga handles
template<typename T>
void GlobalArrays< T >::GA_Sync()
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Sync" )
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    )

    // ensure all GA operations are complete	
    for (Int i = 0; i < ga_handles.size(); i++)
    {
	if (ga_handles[i].ndim == 2 && !ga_handles[i].is_destroyed)
	{	
	    ga_handles[i].rmaint->Flush(); // flush all
	    ga_handles[i].rma_local_pending = false;
	}
    }

    const Grid &grid = DefaultGrid();
    mpi::Barrier( grid.Comm() );
}

// delete all active arrays and destroy internal data structures.
template<typename T>
void GlobalArrays< T >::GA_Terminate()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Terminate" ) )
    
    if (!ga_initialized) // already terminated
	return;

    // destroy GAs as applicable
    for (Int i = 0; i < ga_handles.size(); i++)
    	GA_Destroy( i );
    
    ga_initialized = false;
    
    matrices_.clear();
    ga_handles.clear();
    nbhdls_.clear();
    
    const Grid &grid = DefaultGrid();
    mpi::Barrier( grid.Comm() );
}

// B = A'
template<typename T>
void GlobalArrays< T >::GA_Transpose(Int g_a, Int g_b)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::GA_Transpose" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (g_b < 0 || g_b >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 
		|| ga_handles[g_b].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed
		|| ga_handles[g_b].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    DistMatrix< T, MC, MR >& ADM = *(ga_handles[g_a].DM);
    DistMatrix< T, MC, MR >& BDM = *(ga_handles[g_b].DM);
    
    Transpose( ADM, BDM );
}

template<typename T>
void GlobalArrays< T >::NGA_Distribution( Int g_a, Int iproc, Int lo[], Int hi[] )
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Distribution" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    const Int patchHeight = ga_handles[g_a].patchHeight;
    const Int patchWidth = ga_handles[g_a].patchWidth;
    
    // lo, hi
    if (patchHeight == 0 || patchWidth == 0)
    {
	lo[0] = -1; lo[1] = -1;
	hi[0] = -2; hi[1] = -2;
    }
    else
    {
	const Int pos = iproc * 2;
	hi[0] = ga_handles[g_a].ga_hi[pos];
	hi[1] = ga_handles[g_a].ga_hi[pos + 1];
	lo[0] = ga_handles[g_a].ga_lo[pos];
	lo[1] = ga_handles[g_a].ga_lo[pos + 1];
    }
}

// Inquires the shape of a global array
template<typename T>
void GlobalArrays< T >::NGA_Inquire( Int g_a, Int * ndim, Int dims[] )
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Inquire" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    if (ga_handles[g_a].ndim == 1)
    {
	*ndim = 1;
	*dims = ga_handles[g_a].length;
    }
    else
    {
	DistMatrix< T, MC, MR >&Y = *(ga_handles[g_a].DM);
	*ndim = 2; // number of dims is always 2
	dims[1] = Y.Height();
	dims[0] = Y.Width();
    }
}

// accesses data locally allocated for a global array   
// NOTE The returned leading dimension will be width,
// because in GA, ldim is width of a matrix. Unlike,
// Elemental, where ldim is height. Therefore, in NGA_
// communication functions, we need to perform a transpose
// to correctly put values in the output buffer
template<typename T>
void GlobalArrays< T >::NGA_Access(Int g_a, Int lo[], Int hi[], T** ptr, Int ld[])
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Access" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (lo[0] == -1 && hi[0] == -2)
	    LogicError("Invalid coordinate axes");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    // check lo/hi
    const Int localHeight = ga_handles[g_a].patchHeight;
    const Int localWidth = ga_handles[g_a].patchWidth;
    const Int patchHeight = hi[1] - lo[1] + 1;
    const Int patchWidth = hi[0] - lo[0] + 1;

    if (patchHeight < 0 || patchHeight > localHeight)
	LogicError ("Incorrect coordinates in hi[1],lo[1] supplied");
    if (patchWidth < 0 || patchWidth > localWidth)
	LogicError ("Incorrect coordinates in hi[0],lo[0] supplied");

    Int m_index = -99;
    for (Int i = 0; i < matrices_.size(); i++)
    {
	if (matrices_[i].ga_index_ == g_a 
		&& !matrices_[i].is_accumulate_)
	{
	    m_index = i;
	    break;
	}
    }

    Matrix< T >& M = *(matrices_[m_index].M_);
    matrices_[m_index].lo[0] = lo[0];
    matrices_[m_index].lo[1] = lo[1];
    matrices_[m_index].hi[0] = hi[0];
    matrices_[m_index].hi[1] = hi[1];

    const Int i = lo[1];
    const Int j = lo[0];
    
    // get
    NBXFER ('G', g_a, M, i, j);
    
    *ptr = M.Buffer();
    *ld = localHeight;
}

// this is used when (local portion of) GA is accessed for reading
template<typename T>
void GlobalArrays< T >::NGA_Release(Int g_a, Int lo[], Int hi[])
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Release" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    if (lo[0] == -1 && hi[0] == -2) // means no elements are stored locally
	    return;
	    // Note: This is just checking to see whether 
	    // NGA_Release was called after calling NGA_Access
	    Int m_index = -99;
	    for (Int i = 0; i < matrices_.size(); i++)
	    {
		if (matrices_[i].ga_index_ == g_a 
			&& !matrices_[i].is_accumulate_
			&& matrices_[i].hi[0] == hi[0]
			&& matrices_[i].hi[1] == hi[1]
			&& matrices_[i].lo[0] == lo[0]
			&& matrices_[i].lo[1] == lo[1])
		{
		    m_index = i;
		    break;
		}
	    }
            if (m_index < 0)
		LogicError ("Calling NGA_Release without calling NGA_Access before is not allowed");
	    )
}

// this is used when (local portion of) GA is accessed for writing
template<typename T>
void GlobalArrays< T >::NGA_Release_update(Int g_a, Int lo[], Int hi[])
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Release" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    if (lo[0] == -1 && hi[0] == -2) // means no elements are stored locally
	    return;
	    )
    
    Int m_index = -99;
    for (Int i = 0; i < matrices_.size(); i++)
    {
	if (matrices_[i].ga_index_ == g_a 
		&& !matrices_[i].is_accumulate_
		&& matrices_[i].hi[0] == hi[0]
		&& matrices_[i].hi[1] == hi[1]
		&& matrices_[i].lo[0] == lo[0]
		&& matrices_[i].lo[1] == lo[1])
	{
	    m_index = i;
	    break;
	}
    }
    if (m_index < 0)
	LogicError ("Calling NGA_Release without calling NGA_Access before is not allowed");
    
    // sync global array with local info
    const Int i = lo[1];
    const Int j = lo[0];

    Matrix< T >& M = *(matrices_[m_index].M_);

    // Put - locally blocking transfer
    NBXFER('P', g_a, M, i, j);
}

template<typename T>
void  GlobalArrays< T >::NGA_Acc(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], T* alpha)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Acc" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (buf == NULL)
	    LogicError ("Input buffer cannot be NULL");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    T a = *alpha;
    T one = T( 1 );
    
    // calculate local height and width
    const Int width = hi[0] - lo[0] + 1;
    const Int height = hi[1] - lo[1] + 1;
    const Int ldim = *ld; // ldim for GA buffer

    Matrix< T > M;	

    if (a == one)
	M.Attach( height, width, buf, ldim );
    else
    {
        const Int eldim = Max( height, 1 );
        // create a matrix for nonblocking transfer
        const Int currentIndex = matrices_.size();
        matrices_.push_back( matrix_params_() );

        matrices_[currentIndex].ga_index_ = g_a;
        matrices_[currentIndex].is_accumulate_ = true;

        matrices_[currentIndex].M_ = new Matrix< T >( height, width );
        M = *(matrices_[currentIndex].M_);
	
	T * inbuf = M.Buffer();
	
	for (Int j = 0; j < width; j++)
	    MemCopy( &inbuf[j*eldim], &buf[j*ldim], height );

	const Int numEntries = width * height;
	blas::Scal( numEntries, a, inbuf, 1 );
    }

    const Int i = lo[1];
    const Int j = lo[0];

    // Acc - locally blocking transfer
    NBXFER ('A', g_a, M, i, j);
}

template<typename T>
void GlobalArrays< T >::NGA_Get(Int g_a, Int lo[], Int hi[], T* buf, Int ld[])
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Get" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (buf == NULL)
	    LogicError ("Output buffer cannot be NULL");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    // calculate height and width from lo and hi
    const Int width = hi[0] - lo[0] + 1;
    const Int height = hi[1] - lo[1] + 1;
    const Int ldim = *ld; // ldim for GA buffer
   
    const Int i = lo[1];
    const Int j = lo[0];
	
    Matrix< T > A;
    A.Attach( height, width, buf, ldim );
        
    // get
    NBXFER ('G', g_a, A, i, j);
}

// locally nonblocking scaled accumulate
template<typename T>
void GlobalArrays< T >::NGA_NbAcc(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], T* alpha, ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_NbAcc" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (buf == NULL)
	    LogicError ("Input buffer cannot be NULL");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    T a = *alpha;
    T one = T( 1 );

    // calculate local height and width
    const Int width = hi[0] - lo[0] + 1;
    const Int height = hi[1] - lo[1] + 1;
    const Int ldim = *ld; // ldim for GA Buffer	   

    Matrix< T > M;
        
    if (a == one)
	M.Attach( height, width, buf, ldim );
    else
    {
	const Int eldim = Max( height, 1 );
	// create a matrix for nonblocking transfer
	const Int currentIndex = matrices_.size();
	matrices_.push_back( matrix_params_() );      

	matrices_[currentIndex].ga_index_ = g_a;
	matrices_[currentIndex].is_accumulate_ = true;

	matrices_[currentIndex].M_ = new Matrix< T >( height, width );	
	M = *(matrices_[currentIndex].M_);	

	T * inbuf = M.Buffer();

	for (Int j = 0; j < width; j++)
	    MemCopy( &inbuf[j*eldim], &buf[j*ldim], height );

	const Int numEntries = width * height;
	blas::Scal( numEntries, a, inbuf, 1 );
    }

    const Int i = lo[1];
    const Int j = lo[0];
   
    // Acc - (locally) nonblocking transfer
    LNBXFER ('A', g_a, M, i, j);
    
    // nb handle management
    const Int hdlIndex = nbhdls_.size();
    nbhdls_.push_back( nbhdl_t_() );
    nbhdls_[hdlIndex].M_ = &M;
    nbhdls_[hdlIndex].nbhandle_ = g_a;
    nbhdls_[hdlIndex].active_ = true;

    *nbhandle = g_a;
}

template<typename T>
void GlobalArrays< T >::NGA_NbGet(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_NbGet" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (buf == NULL)
	    LogicError ("Output buffer cannot be NULL");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    // calculate height and width from lo and hi
    const Int width = hi[0] - lo[0] + 1;
    const Int height = hi[1] - lo[1] + 1;
    const Int ldim = *ld; // ldim for GA buffer
   
    const Int i = lo[1];
    const Int j = lo[0];
	
    Matrix< T > A;
    A.Attach( height, width, buf, ldim );
        
    // get
    LNBXFER ('G', g_a, A, i, j);

    // nb handle management
    const Int hdlIndex = nbhdls_.size();
    nbhdls_.push_back( nbhdl_t_() );
    nbhdls_[hdlIndex].nbhandle_ = g_a;
    nbhdls_[hdlIndex].M_ = &A;
    nbhdls_[hdlIndex].active_ = true;
    
    *nbhandle = g_a;
}

template<typename T>
void GlobalArrays< T >::NGA_NbPut(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_NbPut" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (buf == NULL)
	    LogicError ("Input buffer cannot be NULL");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    // calculate local height and width
    const Int width = hi[0] - lo[0] + 1;
    const Int height = hi[1] - lo[1] + 1;
    const Int ldim = *ld; // ldim for GA buffer

    // create a matrix
    Matrix< T > A;
    A.Attach( height, width, buf, ldim );

    const Int i = lo[1];
    const Int j = lo[0];
    
    // Put - (locally) nonblocking transfer
    LNBXFER ('P', g_a, A, i, j);
 
    // nb handle management
    const Int hdlIndex = nbhdls_.size();
    nbhdls_.push_back( nbhdl_t_() );
    nbhdls_[hdlIndex].M_ = &A;
    nbhdls_[hdlIndex].nbhandle_ = g_a;   
    nbhdls_[hdlIndex].active_ = true;
    
    *nbhandle = g_a;
}

// The wait operation ensures only *local* completion
template<typename T>
void GlobalArrays< T >::NGA_NbWait(ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_NbWait" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (ga_handles.size() == 0 || *nbhandle >= ga_handles.size())
	    return;
	    )
    
    if (*nbhandle != -1)
    {
	// find nb handle index
	Int hdlIndex = -1;
	for (Int i = 0; i < nbhdls_.size(); i++)
	{
	    if (*nbhandle == nbhdls_[i].nbhandle_
		    && nbhdls_[i].active_ == true)
	    {
		hdlIndex = i;
		break;
	    }
	}
	if (ga_handles[*nbhandle].rma_local_pending)
	{
	    Matrix< T > M = *(nbhdls_[hdlIndex].M_);
	    ga_handles[*nbhandle].rmaint->LocalFlush( M );
	    ga_handles[*nbhandle].rma_local_pending = false;
	    nbhdls_[hdlIndex].active_ = false;
	}

	// release handle
	*nbhandle = -1;
    }
}

template<typename T>
void GlobalArrays< T >::NGA_Put(Int g_a, Int lo[], Int hi[], T* buf, Int ld[])
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Put" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (buf == NULL)
	    LogicError ("Input buffer cannot be NULL");
	    if ( ga_handles[g_a].ndim == 1 )
	    LogicError ("A 1D GA is not allowed for this operation");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    // calculate local height and width
    const Int width = hi[0] - lo[0] + 1;
    const Int height = hi[1] - lo[1] + 1;
    const Int ldim = *ld; // ldim for GA buffer

    // create a matrix by attaching buffer
    Matrix< T > A;
    A.Attach( height, width, buf, ldim );

    const Int i = lo[1];
    const Int j = lo[0];

    // Put -  local blocking transfer
    NBXFER('P', g_a, A, i, j);
}

// ensures remote completion
template<typename T>
T GlobalArrays< T >::NGA_Read_inc(Int g_a, Int subscript[], T inc)
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "GlobalArrays::NGA_Read_inc" ) 
	    if (!ga_initialized)
	    LogicError ("Global Arrays must be initialized before any operations on the global array");
	    if (g_a < 0 || g_a >= ga_handles.size())
	    LogicError ("Invalid GA handle");
	    if (ga_handles[g_a].is_destroyed)
	    LogicError ("Operation not possible as Global Arrays has been destroyed");
	    )

    T prev;

    if (ga_handles[g_a].ndim == 2) // use RMAInterface atomic function
    {
	const Int i = subscript[1];
	const Int j = subscript[0];
	prev = ga_handles[g_a].rmaint->AtomicIncrement( i, j, inc );
    }
    else
	prev = mpi::ReadInc( ga_handles[g_a].fop_win, *(subscript), inc );
    
    return prev;   
}

#define PROTO(T) template class GlobalArrays<T>;
#ifndef PROTO_INT
# define PROTO_INT(T) PROTO(T)
#endif

#ifndef PROTO_REAL 
# define PROTO_REAL(T) PROTO(T)
#endif
#ifndef PROTO_FLOAT
# define PROTO_FLOAT PROTO_REAL(float)
#endif
#ifndef PROTO_DOUBLE
# define PROTO_DOUBLE PROTO_REAL(double)
#endif

#ifndef EL_NO_INT_PROTO
PROTO_INT(Int)
#endif

#ifndef EL_NO_REAL_PROTO
# if !defined(EL_NO_FLOAT_PROTO)
PROTO_FLOAT
# endif
# if !defined(EL_NO_DOUBLE_PROTO)
PROTO_DOUBLE
# endif
#endif
} // namespace El
#endif // EL_ENABLE_RMA_GLOBAL_ARRAYS
