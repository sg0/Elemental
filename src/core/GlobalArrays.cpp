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
    : ga_initialized (false),
      ga_handles (0)
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
template<typename T>
Int GlobalArrays< T >::GA_Create(Int ndim, Int dims[], const char *array_name)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Create" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array"); 
	
    Int handle = ga_handles.size();
    // call GA constructor
    ga_handles.push_back( GA() );
    
    // 2-D array, create DistMatrix
    if (ndim == 2)
    {
	ga_handles[handle].ndim = 2;
	ga_handles[handle].length = (dims[0] * dims[1]);
	// call rmainterface/dm constructor
	RmaInterface< T > * rmaint = new RmaInterface< T >();
	// create distmatrix over default grid, i.e mpi::COMM_WORLD
	// using MC x MR distribution
	// this won't be allocated until RmaInterface->Attach
	DistMatrix< T > * DM = new DistMatrix< T  >( dims[0], dims[1] );
	DistMatrix< T > &D = *DM;

	const Grid& grid = D.Grid();

	const Int p = grid.Size();
	const Int my_rank = grid.VCRank();
	// create Int vectors for storing local heights and widths
	std::vector< Int > * hvect = new std::vector< Int >();
	std::vector< Int > * wvect = new std::vector< Int >();

	// copy objects 
	ga_handles[handle].DM = DM;
	ga_handles[handle].ga_local_height = hvect;
	ga_handles[handle].ga_local_width = wvect;
	ga_handles[handle].ga_local_height->resize( p );
	ga_handles[handle].ga_local_width->resize( p );
	
	// store local heights and widths
	ga_handles[handle].ga_local_height->at( my_rank ) = DM->LocalHeight();
	ga_handles[handle].ga_local_width->at( my_rank ) = DM->LocalWidth();
	// FIXME mpi allgather in el mpi does not have 
	// MPI_IN_PLACE port
	const Int local_height = ga_handles[handle].ga_local_height->at( my_rank );
	const Int local_width = ga_handles[handle].ga_local_width->at( my_rank );
	mpi::AllGather <Int>( &local_height, 1, ga_handles[handle].ga_local_height->data(), 1, grid.VCComm() );
	mpi::AllGather <Int>( &local_width, 1, ga_handles[handle].ga_local_width->data(), 1, grid.VCComm() );

	// attach DM for RMA ops
	ga_handles[handle].rmaint = rmaint;
	ga_handles[handle].rmaint->Attach( D );
	// zero out allocated DM
	//Zeros (D, dims[0], dims[1]);
    }
    else if (ndim == 1) // fetch-and-op
    {
	const Grid& grid = DefaultGrid();	
	ga_handles[handle].ndim = 1;
	// length of GA
	ga_handles[handle].length = *(dims);
	const Int bufferSize = ga_handles[handle].length * sizeof(T);
	// start access epoch on FOP window
	mpi::WindowAllocate( bufferSize, grid.VCComm(), ga_handles[handle].fop_win );
	mpi::WindowLock( ga_handles[handle].fop_win );
    }
    else
	LogicError ("Up to 2 dimensions supported presently");
	
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
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");

    const Grid& grid = DefaultGrid();
    Int handle = ga_handles.size();

    // call GA constructor
    ga_handles.push_back( GA() );
    const Int ndim = ga_handles[g_a].ndim;

    if (ndim == 2)
    {
	ga_handles[handle].ndim = 2;
	ga_handles[handle].length = ga_handles[g_a].length;
	
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
	ga_handles[handle].DM = DM;
	ga_handles[handle].ga_local_height = hvect;
	ga_handles[handle].ga_local_width = wvect;
	ga_handles[handle].ga_local_height->resize( p );
	ga_handles[handle].ga_local_width->resize( p );
	
	// store local heights and widths
	ga_handles[handle].ga_local_height->at( my_rank ) = DM->LocalHeight();
	ga_handles[handle].ga_local_width->at( my_rank ) = DM->LocalWidth();
	// FIXME mpi allgather in el mpi does not have 
	// MPI_IN_PLACE port
	const Int local_height = ga_handles[handle].ga_local_height->at( my_rank );
	const Int local_width = ga_handles[handle].ga_local_width->at( my_rank );
	mpi::AllGather <Int>( &local_height, 1, ga_handles[handle].ga_local_height->data(), 1, grid.VCComm() );
	mpi::AllGather <Int>( &local_width, 1, ga_handles[handle].ga_local_width->data(), 1, grid.VCComm() );

	// attach DM for RMA ops
	DistMatrix< T > &D = *DM;
	ga_handles[handle].rmaint = rmaint;
	ga_handles[handle].rmaint->Attach( D );  
	// zero out allocated DM
	//Zeros (D, dim[0], dim[1]);
    }
    else if (ndim == 1)// fetch-and-op
    {
	ga_handles[handle].ndim = 1;
	// length of GA
	ga_handles[handle].length = ga_handles[g_a].length;
	const Int bufferSize = ga_handles[handle].length * sizeof(T);
	// start access epoch on FOP window
	mpi::WindowAllocate( bufferSize, grid.VCComm(), ga_handles[handle].fop_win );
	mpi::WindowLock( ga_handles[handle].fop_win );
    }
    else
	LogicError ("Up to 2 dimensions supported presently");
    
    mpi::Barrier( grid.VCComm() );

    return handle;
}

// g_c = alpha * g_a  +  beta * g_b;
template<typename T>
void GlobalArrays< T >::GA_Add(T *alpha, Int g_a, T* beta, Int g_b, Int g_c)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Add" ) )
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

    DistMatrix<T>& GA = *(ga_handles[g_a].DM);
    DistMatrix<T>& GB = *(ga_handles[g_b].DM);
    DistMatrix<T>& GC = *(ga_handles[g_c].DM);

    const Int g_a_height = GA.Height();
    const Int g_a_width = GA.Width();
    const Int g_b_height = GB.Height();
    const Int g_b_width = GB.Width();   
    const Int g_c_height = GC.Height();
    const Int g_c_width = GC.Width();
    
    T a = *alpha;
    T b = *beta;
    
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
    // algorithmic blocksize
    Int nb = 96;
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
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_b < 0 || g_b >= ga_handles.size())
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
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Print" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 )
	LogicError ("A 1D GA is not allowed for this operation");

    DistMatrix<T>& DM = *(ga_handles[g_a].DM);
    Print( DM );
}

// deallocates ga and frees associated resources
// TODO modify code such that destroyed handles could
// be reused
template<typename T>
void GlobalArrays< T >::GA_Destroy(Int g_a)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Destroy" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");

    if (ga_handles[g_a].ndim == 1)
    {
	// clear window object for FOP
	mpi::WindowUnlock( ga_handles[g_a].fop_win );
	mpi::WindowFree ( ga_handles[g_a].fop_win );
    }
    else if (ga_handles[g_a].ndim == 2)
    {
	// detach RmaInterface
	// will end access epoch
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
    else
	LogicError ("GA was already destroyed");

    // nullify pointers
    ga_handles[g_a].ndim = -1;
    ga_handles[g_a].length = -1;
    ga_handles[g_a].pending_rma_op = false;		
    ga_handles[g_a].rmaint = nullptr;
    ga_handles[g_a].DM = nullptr;
    ga_handles[g_a].ga_local_height = nullptr;
    ga_handles[g_a].ga_local_width = nullptr;
}

// A := 0.5 * ( A + A' )
template<typename T>
void GlobalArrays< T >::GA_Symmetrize (Int g_a)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Symmetrize" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
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
    T scale = static_cast<T>( 0.5 );
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

// FIXME GA struct should have a grid pointer,
// currently grid information is fetched from
// DM.Grid() or DefaultGrid(), but this is hacky

// fills a global array with a specific value
// collective
template<typename T>
void GlobalArrays< T >::GA_Fill(Int g_a, T* value)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Fill" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");    
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");

    T a = *value;

    if (ga_handles[g_a].ndim == 1)
    {
	T * fop_base = reinterpret_cast<T *>( mpi::GetWindowBase( ga_handles[g_a].fop_win ) );
	for (Int i = 0; i < ga_handles[g_a].length; i++) fop_base[i] = a;
	// default grid
	const Grid &grid = DefaultGrid();
	mpi::Barrier( grid.VCComm() );
    }
    else // fill DM with alpha
    {
	DistMatrix<T>& DM = *(ga_handles[g_a].DM);
	const Grid &grid = DM.Grid();
	Fill( DM, a );
	mpi::Barrier( grid.VCComm() );
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
   if (g_a < 0 || g_a >= ga_handles.size())
       LogicError ("Invalid GA handle");
   if (g_b < 0 || g_b >= ga_handles.size())
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

    ga_initialized = false;
    ga_handles.clear();
    
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
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if (g_b < 0 || g_b >= ga_handles.size())
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
    if (g_a < 0 || g_a >= ga_handles.size())
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
	hi[0] = local_height - 1;
	hi[1] = local_width - 1;
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
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");

    if (ga_handles[g_a].ndim == 1)
    {
	*ndim = 1;
	*dims = ga_handles[g_a].length;
    }
    else
    {
	DistMatrix< T >&Y = *(ga_handles[g_a].DM);
	*ndim = 2; // number of dims is always 2
	dims[0] = Y.Height();
	dims[1] = Y.Width();
    }
}

// accesses data locally allocated for a global array    
template<typename T>
void GlobalArrays< T >::NGA_Access(Int g_a, Int lo[], Int hi[], T** ptr, Int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Access" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");
    if (lo[0] == -1 && hi[0] == -2)
	LogicError("Invalid coordinate axes");

    DistMatrix< T >&Y = *(ga_handles[g_a].DM);
    *ld = Max( Y.LocalHeight(), 1 );
    // pointer to local portion of DM
    *ptr = Y.Buffer();
}

// NOTE: transfers -- both BXFER and NBXFER are locally
// blocking in keeping with the GA model (does
// not apply to Get, which ensures local+remote
// completion)

// nonblocking transfers

// locally nonblocking
#define LNBXFER(type, g_a, M, i, j) \
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

// blocking transfers -- blocks the calling process
// till transfer is done (remotely completes)
#define BXFER(type, g_a, M, i, j) \
    do { \
        ga_handles[g_a].pending_rma_op = true; \
	switch (type) \
	{ \
	    case 'A': \
		ga_handles[g_a].rmaint->Acc (M, i, j); \
	        ga_handles[g_a].rmaint->Flush (M); \
	        break; \
	    case 'P': \
		ga_handles[g_a].rmaint->Put (M, i, j); \
	        ga_handles[g_a].rmaint->Flush (M); \
	        break; \
	    case 'G': \
		ga_handles[g_a].rmaint->Get (M, i, j); \
	    	break; \
	    default: \
		LogicError ("Unsupported blocking transfer type"); \
	} \
        ga_handles[g_a].pending_rma_op = false; \
    } while (0)

// locally blocking transfer
#define NBXFER(type, g_a, M, i, j) \
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
		LogicError ("Unsupported nonblocking transfer type"); \
	} \
    } while (0)

// NOTE: the input buffer should not be updated, previously, this code
// was performing a * buf[], which changed buf[], which had side effects
template<typename T>
void  GlobalArrays< T >::NGA_Acc(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], T* alpha)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Acc" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    T one = static_cast<T> (1.0);
    T a = *alpha;
    
    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    // create a matrix
    Matrix< T > A;
    A.Attach( height, width, buf, ldim );

    if (a != one)
    {
	for (Int i = 0; i < height; i++)
	    for (Int j = 0; j < width; j++)
		buf[i*ldim + j] *= a;
    }

    const Int i = lo[0];
    const Int j = lo[1];

    // Acc - blocking transfer
    BXFER ('A', g_a, A, i, j);	

    // safe to update local buffer after
    // locally blocking transfer
    if (a != one)
    {
	for (Int i = 0; i < height; i++)
	    for (Int j = 0; j < width; j++)
		buf[i*ldim + j] /= a;
    }
}

template<typename T>
void GlobalArrays< T >::NGA_Get(Int g_a, Int lo[], Int hi[], T* buf, Int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Get" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
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

    // create Matrix<T> for get
    Matrix< T > A; 
    A.Attach( height, width, buf, ldim );

    // get
    BXFER ('G', g_a, A, i, j);
}
 
template<typename T>
void GlobalArrays< T >::NGA_NbAcc(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], T* alpha, ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbAcc" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    T one = static_cast<T> (1.0);
    T a = *alpha;

    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    // create a matrix
    Matrix< T > A (height, width, ldim);
    A.Attach( height, width, buf, ldim );

    if (a != one)
    {
	for (Int i = 0; i < height; i++)
	    for (Int j = 0; j < width; j++)
		buf[i*ldim + j] *= a;
    }

    const Int i = lo[0];
    const Int j = lo[1];
    
    // Acc - (locally) nonblocking transfer
    NBXFER ('A', g_a, A, i, j);
    *nbhandle = g_a;

    // safe to update local buffer after
    // locally blocking transfer
    if (a != one)
    {
	for (Int i = 0; i < height; i++)
	    for (Int j = 0; j < width; j++)
		buf[i*ldim + j] /= a;
    }
}

template<typename T>
void GlobalArrays< T >::NGA_NbGet(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], ga_nbhdl_t* nbhandle)
{
    NGA_Get (g_a, lo, hi, buf, ld);
    *nbhandle = -1; // no nb get implementation in el::rmainterface
}

template<typename T>
void GlobalArrays< T >::NGA_NbPut(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbPut" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    // create a matrix
    Matrix< T > A;
    A.Attach( height, width, buf, ldim );

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
    if (*nbhandle >= ga_handles.size())
	return -1;

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
    if (*nbhandle >= ga_handles.size())
	return;  
    
    if (*nbhandle != -1)
    {
	ga_handles[*nbhandle].rmaint->Waitall ();
        ga_handles[*nbhandle].pending_rma_op = false;
	// release handle
	*nbhandle = -1;
    }
}

template<typename T>
void GlobalArrays< T >::NGA_Put(Int g_a, Int lo[], Int hi[], T* buf, Int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Put" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");
    if (buf == NULL)
	LogicError ("Input buffer cannot be NULL");
    if ( ga_handles[g_a].ndim == 1 )
       LogicError ("A 1D GA is not allowed for this operation");

    // calculate local height and width
    const Int width = hi[1] - lo[1] + 1;
    const Int height = hi[0] - lo[0] + 1;
    const Int ldim = *ld;

    // create a matrix
    Matrix< T > A;
    A.Attach( height, width, buf, ldim );

    const Int i = lo[0];
    const Int j = lo[1];

    // Put - (locally) nonblocking transfer
    BXFER('P', g_a, A, i, j);
}

// ensures remote completion
template<typename T>
T GlobalArrays< T >::NGA_Read_inc(Int g_a, Int subscript[], T inc)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Read_inc" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a >= ga_handles.size())
	LogicError ("Invalid GA handle");

    T prev;

    if (ga_handles[g_a].ndim == 2) // use RMAInterface atomic function
    {
	const Int i = subscript[0];
	const Int j = subscript[1];
	prev = ga_handles[g_a].rmaint->AtomicIncrement( i, j, inc );
    }
    else // rank 0 is default FOP root for 1-D FOP GA, use mpi function directly
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
