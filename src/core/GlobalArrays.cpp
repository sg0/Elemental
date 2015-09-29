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
    : ga_initialized (false), ga_handles (0) 
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

// INFO: type argument will most probably be removed
// from the interface
template<typename T>
Int GlobalArrays< T >::GA_Create(Int type, Int ndim, Int dims[], const char *array_name)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Create" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");   

    Int handle = ga_handles.size();

    Int dim[2] = {1, 1};
    if (ndim == 1)
	dim[0] = dims[0];
    else
    {
	dim[0] = dims[0];
	dim[1] = dims[1];
    }

    // call GA constructor to initialize DM
    GA ga;
    // call rmainterface/dm constructor
    RmaInterface< T > * rmaint = new RmaInterface< T >();
    // create distmatrix over default grid, i.e mpi::COMM_WORLD
    DistMatrix< T > * DM = new DistMatrix< T >( dim[0], dim[1] );
    const Grid& grid = DM->Grid();
    const Int p = grid.Size();
    const Int my_rank = grid.VCRank();
    // create Int vectors for storing local heights and widths
    std:vector< Int > * hvect = new std:vector< Int >( p );
    std:vector< Int > * wvect = new std:vector< Int >( p );
    
    // copy objects 
    ga.rmaint = rmaint;
    ga.DM = DM;
    ga.ga_local_height = hvect;
    ga.ga_local_width = wvect;

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

    DistMatrix< T >& GADM = *(ga_handles[g_a].DM);
    Int dim[2];
    dim[0] = GADM.Height();
    dim[1] = GADM.Width();
    const Grid& grid = GADM.Grid();
    // call GA constructor to initialize DM
    GA ga;
    // call rmainterface/dm constructor
    RmaInterface< T > * rmaint = new RmaInterface< T >();
    DistMatrix< T > * DM = new DistMatrix< T >( dim[0], dim[1], grid );
    const Int p = grid.Size();
    const Int my_rank = grid.VCRank();
    // create Int vectors for storing local heights and widths
    std:vector< Int > * hvect = new std:vector< Int >( p );
    std:vector< Int > * wvect = new std:vector< Int >( p );

    // copy objects 
    ga.rmaint = rmaint;
    ga.DM = DM;
    ga.ga_local_height = hvect;
    ga.ga_local_width = wvect;

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

    return handle;
}

#define FALLBACK_STRIDE 8
#define SCALAR_x_GA(scalar, g_a) \
    do { \
	mpi::Window fop_win; \
	long *fop_win_base; \
	Int local_height, local_width; \
	/* calculate block sizes */ \
	local_height =  ((ga_handles[g_a].DM->Height() / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].DM->Height() / 64)); \
	local_height = (( local_height >= ga_handles[g_a].DM->Height() ) ? (ga_handles[g_a].DM->Height() / 2) : local_height ); \
	local_width =  ((ga_handles[g_a].DM->Width() / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].DM->Width() / 64)); \
	local_width = (( local_width >= ga_handles[g_a].DM->Width() ) ? (ga_handles[g_a].DM->Width() / 2) : local_width ); \
	/* initialize matrix (block, block) to hold a portion of GA */ \
	Matrix < T >A; \
	Zeros( A, local_height, local_width ); \
	/* initialize fetch-and-op window */ \
	mpi::Comm comm = ga_handles[g_a].DM->DistComm(); \
	mpi::WindowAllocate (sizeof (long), comm, fop_win); \
	memset (fop_win_base, 0, sizeof (long)); \
	mpi::WindowLock( fop_win ); \
	long counter = 0, next = 0; \
	next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	/* initialize GA dim and strides */ \
	Int ga_height = ga_handles[g_a].DM->Height(); \
	Int ga_width = ga_handles[g_a].DM->Width(); \
	Int height_stride = local_height; \
	Int width_stride = local_width; \
	for (Int i = 0; i < ga_height; i += height_stride) \
	{ \
	    if (counter == next) \
	    { \
		for (Int j = 0; j < ga_width; j += width_stride) \
		{ \
		    ga_handles[g_a].rmaint->Getx( scalar, A, i, j ); \
		    ga_handles[g_a].rmaint->Iput( A, i, j ); \
		} \
		next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	    } \
	    counter++; \
	} \
	ga_handles[g_a].rmaint->Flush( A ); \
	mpi::Barrier( comm ); \
	mpi::WindowUnlock(fop_win ); \
	mpi::WindowFree(fop_win ); \
    } while (0)

// g_c = alpha * g_a + beta * g_b
// g_c, g_a, g_b are GA struct type,
// not integer handles
// all 3 arrays have same shape
// and distributed over same comm
#define GA_ADD(alpha, beta, g_a, g_b, g_c) \
    do { \
	mpi::Window fop_win; \
	long *fop_win_base; \
	Int local_height, local_width; \
	/* calculate block sizes */ \
	local_height =  ((ga_handles[g_a].DM->Height() / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].DM->Height() / 64)); \
	local_height = (( local_height >= ga_handles[g_a].DM->Height() ) ? (ga_handles[g_a].DM->Height() / 2) : local_height ); \
	local_width =  ((ga_handles[g_a].DM->Width() / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].DM->Width() / 64)); \
	local_width = (( local_width >= ga_handles[g_a].DM->Width() ) ? (ga_handles[g_a].DM->Width() / 2) : local_width ); \
	Matrix < T >A; \
	Matrix < T >B; \
	Zeros( A, local_height, local_width ); \
	Zeros( B, local_height, local_width ); \
	mpi::Comm comm = ga_handles[g_a].DM->DistComm(); \
	mpi::WindowAllocate (sizeof (long), comm, fop_win); \
	memset (fop_win_base, 0, sizeof (long)); \
	mpi::WindowLock( fop_win ); \
	/* initialize GA dim and strides */ \
	Int ga_height = ga_handles[g_a].DM->Height(); \
	Int ga_width = ga_handles[g_a].DM->Width(); \
	Int height_stride = local_height; \
	Int width_stride = local_width; \
	long counter = 0, next = 0; \
	mpi::Barrier( comm ); \
	next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	for (Int i = 0; i < ga_height; i += height_stride) \
	{ \
	    if (counter == next) \
	    { \
		for (Int j = 0; j < ga_width; j += width_stride) \
		{ \
		    ga_handles[g_a].rmaint->Getx( alpha, A, i, j ); \
		    ga_handles[g_b].rmaint->Getx( beta, B, i, j ); \
		    ga_handles[g_c].rmaint->Iput( A, i, j ); \
		    ga_handles[g_c].rmaint->Iacc( B, i, j ); \
		} \
		next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	    } \
	    counter++; \
	} \
	ga_handles[g_a].rmaint->Flush( A ); \
	ga_handles[g_b].rmaint->Flush( B ); \
	ga_handles[g_c].rmaint->Flush( A ); \
	ga_handles[g_c].rmaint->Flush( B ); \
	mpi::Barrier( comm ); \
	mpi::WindowUnlock(fop_win ); \
	mpi::WindowFree(fop_win ); \
    } while (0)

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

    T * a = (T *)alpha;
    T * b = (T *)beta;
    
    // arrays must have the same shape
    if ((ga_handles[g_a].DM->Width() != ga_handles[g_b].DM->Width()) 
	    || (ga_handles[g_a].DM->Width() != ga_handles[g_c].DM->Width()) 
	    || (ga_handles[g_b].DM->Width() != ga_handles[g_c].DM->Width()))
	LogicError ("Global Arrays of different widths cannot be added. ");
    if ((ga_handles[g_a].DM->Height() != ga_handles[g_b].DM->Height()) 
	    || (ga_handles[g_a].DM->Height() != ga_handles[g_c].DM->Height()) 
	    || (ga_handles[g_b].DM->Height() != ga_handles[g_c].DM->Height()))
	LogicError ("Global Arrays of different heights cannot be added. ");

    // add
    GA_ADD (*a, *b, g_a, g_b, g_c);
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

    // arrays must have the same shape
    if ((ga_handles[g_a].DM->Width() != ga_handles[g_b].DM->Width()))
	LogicError ("Global Arrays of different widths cannot be copied. ");
    if ((ga_handles[g_a].DM->Height() != ga_handles[g_b].DM->Height()))
	LogicError ("Global Arrays of different heights cannot be copied. ");

    DistMatrix<T>& ADM = *(ga_handles[g_a].DM);
    DistMatrix<T>& BDM = *(ga_handles[g_b].DM);
    // copy
    Copy( ADM, BDM );
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

    // detach RmaInterface
    ga_handles[g_a].rmaint->Detach();
    ga_handles[g_a].DM->Empty();
    ga_handles[g_a].ga_local_height->clear();
    ga_handles[g_a].ga_local_width->clear();
    // erase
    ga_handles.erase( ga_handles.begin() + g_a );
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

    Int width = ga_handles[g_a].DM->Width();
    Int height = ga_handles[g_a].DM->Height();

    // create intermediate GAs using g_a parameters
    Int dims[2];
    dims[0] = height;
    dims[1] = width;
    Int g_at = GA_Create( 0, 2, dims, "transpose array");
    Int g_c = GA_Create( 0, 2, dims, "final array");

    // transpose
    GA_Transpose( g_a, g_at );

    // add: g_a + g_at
    // g_c = alpha * g_a  +  beta * g_b;
    double *scale;
    *scale = 0.5;
    GA_Add((void *)scale, g_a, (void *)scale, g_at, g_c);

    // copy, g_a = g_c
    GA_Copy (g_a, g_c);
    
    // destroy the intermediate GAs
    GA_Destroy (g_at);
    GA_Destroy (g_c);
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

    DistMatrix<T>& ADM = *(ga_handles[g_a].DM);
    DistMatrix<T>& BDM = *(ga_handles[g_b].DM);
    DistMatrix<T>& CDM = *(ga_handles[g_c].DM);

    // set algorithmic block size	
    Int nb = 96;
    GemmAlgorithm alg = GEMM_SUMMA_A;
    
    const Orientation orientA = CharToOrientation( ((ta == 'T' || ta == 't') ? 'T' : 'N') );
    const Orientation orientB = CharToOrientation( ((ta == 'T' || ta == 't') ? 'T' : 'N') );
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
	alpha = (T *)(value);
    else
	alpha = (T *)(&zero);

    DistMatrix<T>& DM = *(ga_handles[g_a].DM);
    // fill DM with alpha
    Fill( DM, *alpha );
}

// allocate and initialize internal data structures in Global Arrays.
template<typename T>
void GlobalArrays< T >::GA_Initialize()
{
   DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Initialize" ) )

   if (!ga_initialized)
       ga_initialized = true;
}

// barrier
// TODO dont sync over comm world
// fetch comm from ga handles
template<typename T>
void GlobalArrays< T >::GA_Sync()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Sync" ) )

    // ensure all GA operations are complete	
    for (Int i = 0; i < ga_handles.size(); i++)
	    ga_handles[i].rmaint->Waitall ();

    mpi::Barrier (mpi::COMM_WORLD);
}

// delete all active arrays and destroy internal data structures.
template<typename T>
void GlobalArrays< T >::GA_Terminate()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Terminate" ) )
    
    if (!ga_initialized) // already terminated
	return;

    ga_initialized = false;
    // detach if not yet already
    for (Int i = 0; i < ga_handles.size(); i++)
    {
	ga_handles[i].rmaint->Detach();
	ga_handles[i].DM->Empty();
	ga_handles[i].ga_local_height->clear();
	ga_handles[i].ga_local_width->clear();
    }
    ga_handles.clear();

    mpi::Barrier (mpi::COMM_WORLD);
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

// accesses data locally allocated for a global array    
template<typename T>
void GlobalArrays< T >::NGA_Access(Int g_a, Int lo[], Int hi[], void** ptr, Int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Access" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");
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
	} \
    } while (0)

// locally nonblocking transfers
#define NBXFER(type, g_a, M, i, j) \
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
		ga_handles[g_a].rmaint->Get (M, i, j); \
	    	break; \
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

    T zero = static_cast<T> (0.0);
    T one = static_cast<T> (1.0);
    T * a;

    if (alpha == NULL)
	a = (T *)&one;
    else
	a = (T *)alpha;
    T * buffer = (T *) buf;
    
    if (*a == zero) // no-op
	return;

    // calculate local height and width
    Int width = hi[1] - lo[1] + 1;
    Int height = hi[0] - lo[0] + 1;
    Int ldim = ld[0];

    // create a matrix
    Matrix< T > A;
    Zeros( A, height, width );
    T * source = A.Buffer ();

    if (*a != one) // then a * buf
    {
	for (Int i = 0; i < width; i++)
	    for (Int j = 0; j < height; j++)
		buffer[i + j*ldim] *= *a;
    }

    // memcopy to source
    MemCopy( source, buffer, (width * height) );

    const Int i = lo[0];
    const Int j = lo[1];

    // Acc - (locally) blocking transfer
    BXFER ('A', g_a, A, i, j);		
    ga_handles[g_a].rmaint->Flush (A);
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

     // calculate height and width from lo and hi
    Int width = hi[1] - lo[1] + 1;
    Int height = hi[0] - lo[0] + 1;
   
    const Int i = lo[0];
    const Int j = lo[1];

    // create Matrix<T> for get
    Matrix< T > A; 
    Zeros( A, height, width );
    BXFER ('G', g_a, A, i, j);

    T * buffer = A.Buffer();
    T * outbuf = reinterpret_cast<T *>(buf);

    MemCopy( outbuf, buffer, (width * height) );
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

    T zero = static_cast<T> (0.0);
    T one = static_cast<T> (1.0);
    T * a;

    if (alpha == NULL)
	a = (T *)&one;
    else
	a = (T *)alpha;
    T * buffer = reinterpret_cast<T *>( buf );

    if (*a == zero) // no-op
	return;

    // calculate local height and width
    Int width = hi[1] - lo[1] + 1;
    Int height = hi[0] - lo[0] + 1;
    Int ldim = ld[0];

    // create a local matrix
    Matrix< T > A;
    Zeros( A, height, width );
    T * source = A.Buffer ();

    if (*a != one) // then a * buf
    {
	for (Int i = 0; i < width; i++)
	    for (Int j = 0; j < height; j++)
		buffer[i + j*ldim] *= *a;
    }

    // memcopy to source
    MemCopy( source, buffer, (width * height) );

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

    T * buffer = reinterpret_cast<T *>( buf ); 

    // calculate local height and width
    Int width = hi[1] - lo[1] + 1;
    Int height = hi[0] - lo[0] + 1;
    Int ldim = ld[0];

    // create a matrix
    Matrix< T > A;
    Zeros( A, height, width );
    T * source = A.Buffer ();

    // memcopy to source
    MemCopy( source, buffer, (width * height) );

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

    T * buffer = reinterpret_cast<T *>( buf );

    // calculate local height and width
    Int width = hi[1] - lo[1] + 1;
    Int height = hi[0] - lo[0] + 1;

    // create a matrix
    Matrix< T > A;
    Zeros( A, height, width );
    T * source = A.Buffer ();

    // memcopy to source
    MemCopy( source, buffer, (width * height) );

    const Int i = lo[0];
    const Int j = lo[1];

    // Put - (locally) nonblocking transfer
    BXFER('P', g_a, A, i, j);
    ga_handles[g_a].rmaint->Flush (A);
}

template<typename T>
long GlobalArrays< T >::NGA_Read_inc(Int g_a, Int ndim, Int subscript[], long inc)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Read_inc" ) )
    if (!ga_initialized)
	LogicError ("Global Arrays must be initialized before any operations on the global array");
    if (g_a < 0 || g_a > ga_handles.size())
	LogicError ("Invalid GA handle");

    long prev;

    if (ndim == 2)
    {
	Int i = subscript[0];
	Int j = subscript[1];
	prev = ga_handles[g_a].rmaint->CompareAndSwap (i, j, inc);
    }
    else // 1 dimension
    {
	Int j = subscript[0];
	prev = ga_handles[g_a].rmaint->CompareAndSwap (1, j, inc);
    }

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
