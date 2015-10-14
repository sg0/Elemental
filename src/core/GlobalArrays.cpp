/*
   Copyright (c) 2009-2014, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   Copyright (c) 2014, Jeff Hammond (Intel)
   Copyright (c) 2014, Sayan Ghosh (Washington State University)
   All rights reserved.

Authors:
Sayan Ghosh

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
    : ga_initialized (false), ga_handles (0), 
      fop_win (mpi::WIN_NULL) 
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

template<typename T>
Int GlobalArrays< T >::GA_Create(Int ndim, Int dims[], const char *array_name)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Create" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array"); 
    // FIXME ndim = 1 implicitly assumes fetch-and-op, hence 
    // no DistMatrix is created and entire RMAInterface is bypassed

    Int handle = ga_handles.size();

    // call GA constructor to initialize DM
    GA ga;
    ga.ndim = ndim;
    
    // 2-D array, create DistMatrix
    if (ndim > 1)
    {
	// call rmainterface/dm constructor
	RmaInterface< T > * rmaint = new RmaInterface< T >();
	// create distmatrix over default grid, i.e mpi::COMM_WORLD
	// using MC x MR distribution
	// this won't be allocated until RmaInterface->Attach
	DistMatrix< T > * DM = new DistMatrix< T  >( dims[0], dims[1] );

	const Grid& grid = DM->Grid();

	const Int p = grid.Size();
	const Int my_rank = grid.VCRank();
	// create Int vectors for storing local heights and widths
	std::vector< Int > * hvect = new std::vector< Int >();
	std::vector< Int > * wvect = new std::vector< Int >();

	// copy objects 
	ga.rmaint = rmaint;
	ga.DM = DM;
	ga.ga_local_height = hvect;
	ga.ga_local_width = wvect;
	ga.ga_local_height->resize( p );
	ga.ga_local_width->resize( p );
	ga.pending_rma_op = false;

	// store local heights and widths
	ga.ga_local_height->at( my_rank ) = DM->LocalHeight();
	ga.ga_local_width->at( my_rank ) = DM->LocalWidth();
	// FIXME mpi allgather in el mpi does not have 
	// MPI_IN_PLACE port
	const Int local_height = ga.ga_local_height->at( my_rank );
	const Int local_width = ga.ga_local_width->at( my_rank );
	mpi::AllGather <Int>( &local_height, 1, ga.ga_local_height->data(), 1, grid.VCComm() );
	mpi::AllGather <Int>( &local_width, 1, ga.ga_local_width->data(), 1, grid.VCComm() );

	// push into vector
	ga_handles.push_back( ga );

	// attach DM for RMA ops
	DistMatrix< T > &D = *DM;
	ga_handles[handle].rmaint->Attach( D );
	// zero out allocated DM
	Zeros (D, dims[0], dims[1]);
    }
    else // fetch-and-op
    {
	const Grid& grid = DefaultGrid();
	// length of GA
	ga.length = *(dims);
	const Int bufferSize = ga.length * sizeof(T);
	// start access epoch on FOP window
	mpi::WindowAllocate( bufferSize, grid.VCComm(), fop_win );
	mpi::WindowLock( fop_win );

	ga.pending_rma_op = false;
	// push into vector
	ga_handles.push_back( ga );
    }
	
    const Grid& grid = DefaultGrid();
    mpi::Barrier( grid.VCComm() );

    return handle;
}

// DM might be resized here
template<typename T>
Int GlobalArrays< T >::GA_Duplicate(Int g_a, const char *array_name)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Duplicate" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");   
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");

    Int handle = ga_handles.size();
    const Int ndim = ga_handles[g_a].ndim;
    // call GA constructor to initialize DM
    GA ga;
    ga.ndim = ndim;

    if (ndim > 1)
    {
	DistMatrix< T >& GADM = *(ga_handles[g_a].DM);
	Int dim[2];
	dim[0] = GADM.Height();
	dim[1] = GADM.Width();
	const Grid& grid = GADM.Grid();

	// call rmainterface/dm constructor
	RmaInterface< T > * rmaint = new RmaInterface< T >();
	DistMatrix< T > * DM = new DistMatrix< T >( dim[0], dim[1], grid );
	const Int p = grid.Size();
	const Int my_rank = grid.VCRank();
	// create Int vectors for storing local heights and widths
	std::vector< Int > * hvect = new std::vector< Int >();
	std::vector< Int > * wvect = new std::vector< Int >();

	// copy objects 
	ga.rmaint = rmaint;
	ga.DM = DM;
	ga.ga_local_height = hvect;
	ga.ga_local_width = wvect;
	ga.ga_local_height->resize( p );
	ga.ga_local_width->resize( p );
	ga.pending_rma_op = false;

	// store local heights and widths
	ga.ga_local_height->at( my_rank ) = DM->LocalHeight();
	ga.ga_local_width->at( my_rank ) = DM->LocalWidth();
	// FIXME mpi allgather in el mpi does not have 
	// MPI_IN_PLACE port
	const Int local_height = ga.ga_local_height->at( my_rank );
	const Int local_width = ga.ga_local_width->at( my_rank );
	mpi::AllGather <Int>( &local_height, 1, ga.ga_local_height->data(), 1, grid.VCComm() );
	mpi::AllGather <Int>( &local_width, 1, ga.ga_local_width->data(), 1, grid.VCComm() );

	// push into vector
	ga_handles.push_back( ga );

	// attach DM for RMA ops
	DistMatrix< T > &D = *DM;
	ga_handles[handle].rmaint->Attach( D );  
	// zero out allocated DM
	Zeros (D, dims[0], dims[1]);
    }
    else // fetch-and-op
    {
	const Grid& grid = DefaultGrid();
	// length of GA
	ga.length = ga_handles[g_a].length;
	const Int bufferSize = ga.length * sizeof(T);
	// start access epoch on FOP window
	mpi::WindowAllocate( bufferSize, grid.VCComm(), fop_win );
	mpi::WindowLock( fop_win );

	ga.pending_rma_op = false;
	// push into vector
	ga_handles.push_back( ga );
    }
    
    const Grid& grid = DefaultGrid();
    mpi::Barrier( grid.VCComm() );

    return handle;
}

// g_c = alpha * g_a  +  beta * g_b;
template<typename T>
void GlobalArrays< T >::GA_Add(void *alpha, Int g_a, void* beta, Int g_b, Int g_c)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Add" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_b < 0 || g_b > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_c < 0 || g_c > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 
	    || ga_handles[g_b].ndim == 1 
	    || ga_handles[g_c].ndim == 1 )
	LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix<T>& GA = *(ga_handles[g_a].DM);
    DistMatrix<T>& GB = *(ga_handles[g_b].DM);
    DistMatrix<T>& GC = *(ga_handles[g_c].DM);

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

    DistMatrix< T, MC, MR > Bd;
    Bd.Resize( g_a_height, g_a_width );
    Identity( Bd, g_a_height, g_a_width );

    // add
    // FIXME dont hardcode
    Int nb = 128;
    GemmAlgorithm alg = GEMM_SUMMA_A;
    
    const Orientation orientA = CharToOrientation( 'N' );
    const Orientation orientB = CharToOrientation( 'N' );
    SetBlocksize( nb );

    // GB = a*GA*Bd + b*GB
    Gemm( orientA, orientB, a, GA, Bd, b, GB, alg);
    // GC = GB
    Copy( GB, GC );
}

// copies g_a into g_b, must be of same shape
template<typename T>
void GlobalArrays< T >::GA_Copy(Int g_a, Int g_b)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Allocate" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_b < 0 || g_b > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 
	    || ga_handles[g_b].ndim == 1 )
	LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix<T>& GA = *(ga_handles[g_a].DM);
    DistMatrix<T>& GB = *(ga_handles[g_b].DM);

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
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Print_distribution" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 )
	LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix<T>& DM = *(ga_handles[g_a].DM);
    Print( DM );
}

// deallocates ga and frees associated resources
template<typename T>
void GlobalArrays< T >::GA_Destroy(Int g_a)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Destroy" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");

    if (ga_handles[g_a].ndim == 1)
    {
	// clear window object for FOP
	mpi::WindowUnlock( fop_win );
	mpi::WindowFree ( fop_win );
    }
    else
    {
	// detach RmaInterface
	ga_handles[g_a].rmaint->Detach();
	// set matrix size to 0 and free mem
	ga_handles[g_a].DM->EmptyData();
	// clear vectors that holds local dims
	ga_handles[g_a].ga_local_height->clear();
	ga_handles[g_a].ga_local_width->clear();
	// erase ga entry from global ga_handles vector
	// FIXME erasing would mess up g_a handle values, 
	// as GAs could be destroyed at different times
	// so at present vector is cleared only on terminate
	// ga_handles.erase( ga_handles.begin() + g_a );
    }
}

// A := 0.5 * ( A + A' )
template<typename T>
void GlobalArrays< T >::GA_Symmetrize (Int g_a)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Symmetrize" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 )
	LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix<T>& GA = *(ga_handles[g_a].DM);
    Int width = GA.Width();
    Int height = GA.Height();

    // create intermediate GAs using g_a parameters
    Int dims[2];
    dims[0] = height;
    dims[1] = width;
    Int g_at = GA_Create( 2, dims, "transpose array" );
    Int g_c = GA_Create( 2, dims, "final array" );

    // transpose
    GA_Transpose( g_a, g_at );

    // add: scale * ( g_a + g_at )
    // g_c = alpha * g_a  +  beta * g_at;
    double scale = 0.5;
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
void GlobalArrays< T >::GA_Dgemm(char ta, char tb, Int m, Int n, Int k, double alpha, Int g_a, Int g_b, double beta, Int g_c )
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Dgemm" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_b < 0 || g_b > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_c < 0 || g_c > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 
	    || ga_handles[g_b].ndim == 1 
	    || ga_handles[g_c].ndim == 1 )
	LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix<T>& ADM = *(ga_handles[g_a].DM);
    DistMatrix<T>& BDM = *(ga_handles[g_b].DM);
    DistMatrix<T>& CDM = *(ga_handles[g_c].DM);

    // set algorithmic block size	
    Int nb = 96;
    GemmAlgorithm alg = GEMM_SUMMA_B;
    
    const Orientation orientA = CharToOrientation( ((ta == 'T' || ta == 't') ? 'T' : 'N') );
    const Orientation orientB = CharToOrientation( ((tb == 'T' || tb == 't') ? 'T' : 'N') );
    SetBlocksize( nb );

    T a = static_cast<T>( alpha );
    T b = static_cast<T>( beta );

    Gemm( orientA, orientB, a, ADM, BDM, b, CDM, alg);
}

// fills a global array with a specific value
template<typename T>
void GlobalArrays< T >::GA_Fill(Int g_a, void *value)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Fill" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
   
    T zero = static_cast<T>( 0.0 );
    T * alpha;
    if (value == NULL)
	alpha = (T *)(&zero);
    else
	alpha = (T *)(value);

    if (ga_handles[g_a].ndim == 1)
    {
	T * fop_base = reinterpret_cast<T *>( mpi::GetWindowBase( fop_win ) );
	for (Int i = 0; i < ga_handles[g_a].length; i++) fop_base[i] = *alpha;
    }
    else
    {
	DistMatrix<T>& DM = *(ga_handles[g_a].DM);
	// fill DM with alpha
	Fill( DM, *alpha );
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
		mpi::AllReduce( x, n, mpi::SUM, grid.VCComm() ); \
	        break; \
	    case '*': \
		mpi::AllReduce( x, n, mpi::PROD, grid.VCComm() ); \
	        break; \
	    case 'X': \
		mpi::AllReduce( x, n, mpi::MAX, grid.VCComm() ); \
	    	break; \
	    case 'N': \
		mpi::AllReduce( x, n, mpi::MIN, grid.VCComm() ); \
	    	break; \
	    default: \
		LogicError ("Unsupported global operation specified"); \
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
   DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Dot" ) )
   if (!ga_initialized)
       LogicError ("Global Arrays must be initialized before any operations");
   if (g_a < 0 || g_a > ga_handles.size())
       LogicError ("Invalid GA handle");
   if (g_b < 0 || g_b > ga_handles.size())
       LogicError ("Invalid GA handle");
   if ( ga_handles[g_a].ndim == 1 
	   || ga_handles[g_b].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

   const DistMatrix<T>& A = *(ga_handles[g_a].DM);
   const DistMatrix<T>& B = *(ga_handles[g_b].DM);

   T result = Dot( A, B );

   return result;
}

// barrier
// TODO dont sync over comm world
// fetch comm from ga handles
template<typename T>
void GlobalArrays< T >::GA_Sync()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Sync" ) )
    if (!ga_initialized)
       LogicError ("Global Arrays must be initialized before any operations on the global array");

    // ensure all GA operations are complete	
    for (Int i = 0; i < ga_handles.size(); i++)
    {
	if (ga_handles[i].pending_rma_op && ga_handles[i].ndim > 1)
	    ga_handles[i].rmaint->Flush();
    }

    const Grid &grid = DefaultGrid();
    mpi::Barrier( grid.VCComm() );
}

// delete all active arrays and destroy internal data structures.
template<typename T>
void GlobalArrays< T >::GA_Terminate()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Terminate" ) )
    
    if (!ga_initialized) // already terminated
	return;

    // detach if not yet already
    for (Int i = 0; i < ga_handles.size(); i++)
    {
	GA_Destroy (i);
	ga_handles.erase( ga_handles.begin() + i );
    }
    ga_handles.clear();
    
    ga_initialized = false;
    const Grid &grid = DefaultGrid();
    mpi::Barrier( grid.VCComm() );
}

// B = A'
template<typename T>
void GlobalArrays< T >::GA_Transpose(Int g_a, Int g_b)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Transpose" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_b < 0 || g_b > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 
	   || ga_handles[g_b].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix<T>& ADM = *(ga_handles[g_a].DM);
    DistMatrix<T>& BDM = *(ga_handles[g_b].DM);
    
    Transpose( ADM, BDM );
}

// returns a data range to my PE from iproc
// assuming GA contiguous default distribution
// Note: El data distribution is not contiguous
// TODO
template<typename T>
void GlobalArrays< T >::NGA_Distribution( Int g_a, Int iproc, Int lo[], Int hi[] )
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Distribution" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    const Int local_height = ga_handles[g_a].ga_local_height->at( iproc );
    const Int local_width = ga_handles[g_a].ga_local_width->at( iproc );

    // lo, hi
    if (local_height == 0 || local_width == 0)
    {
	lo[0] = -1; lo[1] = -1;
	hi[0] = -2; hi[1] = -2;
    }
    else
    {
	hi[0] = local_height;
	hi[1] = local_width;
	lo[0] = 0;
	lo[1] = 0;
    }
}

// Inquires the shape of a global array
template<typename T>
void GlobalArrays< T >::NGA_Inquire( Int g_a, Int * ndim, Int dims[] )
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Inquire" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix< T >&Y = *(ga_handles[g_a].DM);
    *ndim = 2; // number of dims is always 2
    dims[0] = Y.Height();
    dims[1] = Y.Width();
}

// accesses data locally allocated for a global array    
template<typename T>
void GlobalArrays< T >::NGA_Access(Int g_a, Int lo[], Int hi[], void** ptr, Int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Access" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");
    if (lo[0] == -1 && hi[0] == -2)
	LogicError("Invalid coordinate axes");

    DistMatrix< T >&Y = *(ga_handles[g_a].DM);
    T ** buffer = reinterpret_cast<T **>( ptr );
    // pointer to local portion of DM
    *buffer = Y.Buffer();
}

// transfers
// locally blocking transfers
#define BXFER(type, g_a, M, i, j) \
    do { \
        ga_handles[g_a].pending_rma_op = true; \
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
		LogicError ("Unsupported blocking transfer type"); \
	} \
    } while (0)

// locally nonblocking transfers
#define NBXFER(type, g_a, M, i, j) \
    do { \
        ga_handles[g_a].pending_rma_op = true; \
	switch (type) \
	{ \
	    case 'A': \
		ga_handles[g_a].rmaint->Iacc (M, i, j); \
	        break; \
	    case 'P': \
		ga_handles[g_a].rmaint->Iput (M, i, j); \
	        break; \
	    case 'G': \
		ga_handles[g_a].rmaint->Get (M, i, j); \
	    	break; \
	    default: \
		LogicError ("Unsupported nonblocking transfer type"); \
	} \
    } while (0)

template<typename T>
void  GlobalArrays< T >::NGA_Acc(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], void* alpha)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Acc" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    T one = static_cast<T> (1.0);
    T * a;

    if (alpha == NULL)
	a = (T *)&one;
    else
	a = (T *)alpha;
    
    T * buffer = reinterpret_cast<T *>( buf );
    
    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    if (*a != one) // then a * buf
    {
	for (Int i = 0; i < width; i++)
	    for (Int j = 0; j < height; j++)
		buffer[i + j*ldim] *= *a;
    }

    // create a matrix
    Matrix< T > A;
    A.Attach( height, width, buffer, ldim );

    const Int i = lo[0];
    const Int j = lo[1];

    // Acc - (locally) blocking transfer
    BXFER ('A', g_a, A, i, j);		
}

template<typename T>
void GlobalArrays< T >::NGA_Get(Int g_a, Int lo[], Int hi[], void* buf, Int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Get" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Output buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    // calculate height and width from lo and hi
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;
   
    const Int i = lo[0];
    const Int j = lo[1];

    T * buffer = reinterpret_cast<T *>( buf );
    // create Matrix<T> for get
    Matrix< T > A (height, width); 
    // get
    BXFER ('G', g_a, A, i, j);
    T * inbuf = A.Buffer();

    for (Int i = 0; i < width; i++)
	for (Int j = 0; j < height; j++)
	    buffer[i + j*ldim] = inbuf[i + j*ldim];
}
 
template<typename T>
void GlobalArrays< T >::NGA_NbAcc(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], void* alpha, ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbAcc" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    T one = static_cast<T> (1.0);
    T * a;

    if (alpha == NULL)
	a = (T *)&one;
    else
	a = (T *)alpha;
    
    T * buffer = reinterpret_cast<T *>( buf );

    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    if (*a != one) // then a * buf
    {
	for (Int i = 0; i < width; i++)
	    for (Int j = 0; j < height; j++)
		buffer[i + j*ldim] *= *a;
    }

    // create a local matrix
    Matrix< T > A;
    A.Attach( height, width, buffer, ldim );

    const Int i = lo[0];
    const Int j = lo[1];
    
    // Acc - (locally) nonblocking transfer
    NBXFER ('A', g_a, A, i, j);
    *nbhandle = g_a;
}

template<typename T>
void GlobalArrays< T >::NGA_NbGet(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], ga_nbhdl_t* nbhandle)
{
    NGA_Get (g_a, lo, hi, buf, ld);
    *nbhandle = -1; // no nb get implementation in el::rmainterface
}

template<typename T>
void GlobalArrays< T >::NGA_NbPut(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbPut" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    T * buffer = reinterpret_cast<T *>( buf ); 

    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    // create a matrix
    Matrix< T > A;
    A.Attach( height, width, buffer, ldim );

    const Int i = lo[0];
    const Int j = lo[1];
    
    // Put - (locally) nonblocking transfer
    NBXFER ('P', g_a, A, i, j);
    *nbhandle = g_a;
}
 
template<typename T>
Int GlobalArrays< T >::NGA_NbTest(ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbTest" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");

    if (*nbhandle != -1)
    {
	bool status = ga_handles[*nbhandle].rmaint->Testall ();
	if (status) // release handle
	{
            ga_handles[*nbhandle].pending_rma_op = false;
	    *nbhandle = -1;
	    return 1;
	}
    }

    return 0;
}
 
template<typename T>
void GlobalArrays< T >::NGA_NbWait(ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbWait" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
   
    if (*nbhandle != -1)
    {
	ga_handles[*nbhandle].rmaint->Waitall ();
        ga_handles[*nbhandle].pending_rma_op = false;
	// release handle
	*nbhandle = -1;
    }
}

template<typename T>
void GlobalArrays< T >::NGA_Put(Int g_a, Int lo[], Int hi[], void* buf, Int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Put" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    T * buffer = reinterpret_cast<T *>( buf );

    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    // create a matrix
    Matrix< T > A;
    A.Attach( height, width, buffer, ldim );

    const Int i = lo[0];
    const Int j = lo[1];

    // Put - (locally) nonblocking transfer
    BXFER('P', g_a, A, i, j);
    ga_handles[g_a].rmaint->Flush (A);
}

// FIXME At present, it is implicitly assumed that a 1-D GA will
// *only* be used for read-increment/fetch-and-op atomic
// operation, and a 2-D GA would not work for this case
template<typename T>
long GlobalArrays< T >::NGA_Read_inc(Int g_a, Int subscript[], long inc)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Read_inc" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
    if (*(subscript) < 0 || *(subscript) > ga_handles[g_a].length)
	LogicError ("Invalid subscript for read increment operation");
    if ( ga_handles[g_a].ndim > 1 )
       LogicError ("A 2-D GA is not allowed for this operation");

    // rank 0 is FOP root
    long prev = mpi::ReadInc( fop_win, *(subscript), inc, 0 );
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
