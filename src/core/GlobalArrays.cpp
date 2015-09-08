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
// constructor    
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
	    "that required a call to Detach. Instead of allowing for the "
	    "possibility of Detach throwing another exception and "
	    "resulting in a 'terminate', we instead immediately dump the "
	    "call stack (if not in RELEASE mode) since the program will "
	    "likely hang:" << std::endl;
	std::cerr << os.str();
	DEBUG_ONLY( DumpCallStack() )
    }
    else
	GA_Terminate();
}

// create an integer handle and return
template<typename T>
int GlobalArrays< T >::GA_Create_handle()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Create_handle" ) )

    const Int numCreated = ga_handles.size();
    ga_handles.resize( numCreated + 1 );
    ga_handles[numCreated].handle = numCreated;
    ga_handles[numCreated].status = CREATED;

    return numCreated;
}

// set dims for the handle
template<typename T>
void GlobalArrays< T >::GA_Set_data (int g_a, int ndim, int dims[], int type)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Set_data" ) )

    if (ndim > 2)
	LogicError ("Up to 2-D arrays are supported. ");

    ga_handles[g_a].ndims = ndim;
    for (int i = 0; i < ndim; i++)
	ga_handles[g_a].dims[i] = dims[i];
    if (ndim == 1)
	ga_handles[g_a].dims[1] = 1;
    // default grid
    // TODO implement c api functions to manipulate pgroups
    ga_handles[g_a].comm = mpi::COMM_WORLD;
    ga_handles[g_a].status = SET;
}

// last stage of GA creation,
// assumes GA_Set_data and GA_Create_handle
// is already called
template<typename T>
int GlobalArrays< T >::GA_Allocate(int g_a)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Allocate" ) )

    Grid grid (ga_handles[g_a].comm);
    DistMatrix< T > D (grid); 
    D.Resize( ga_handles[g_a].dims[0], ga_handles[g_a].dims[1] ); 
    Zeros( D, ga_handles[g_a].dims[0], ga_handles[g_a].dims[1] );
    ga_handles[g_a].DM = D;
    
    ga_handles[g_a].rmaint.Attach( ga_handles[g_a].DM );
    
    ga_handles[g_a].status = ALLOCATED;
    ga_handles[g_a].pending_transfer = false;
    
    return 1;
}

// creates a new array with the same properties as 
// the given array. new array handle returned
template<typename T>
int GlobalArrays< T >::GA_Duplicate(int g_a, const char *array_name)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Duplicate" ) )

    int g_b = GA_Create_handle ();

    ga_handles[g_b].dims[0] = ga_handles[g_a].dims[0];
    ga_handles[g_b].dims[1] = ga_handles[g_a].dims[1];
    ga_handles[g_b].status = ga_handles[g_a].status;
    ga_handles[g_b].pending_transfer = false;

    ga_handles[g_b].comm = ga_handles[g_a].comm;
    if (ga_handles[g_b].dims[0] > 0 && ga_handles[g_b].dims[1] > 0)
    {
	Grid grid (ga_handles[g_a].comm);
	DistMatrix< T > D (grid); 
	D.Resize( ga_handles[g_a].dims[0], ga_handles[g_a].dims[1] ); 
	Zeros( D, ga_handles[g_a].dims[0], ga_handles[g_a].dims[1] );
	ga_handles[g_b].DM = D;
    }
 
    if (ga_handles[g_a].status == ALLOCATED)
    	ga_handles[g_b].rmaint.Attach (ga_handles[g_b].DM);

    return g_b;
}


#define FALLBACK_STRIDE 8
#define SCALAR_x_GA(scalar, g_a) \
    do { \
	mpi::Window fop_win; \
	long *fop_win_base; \
	int ga_local_dims[2] = {1, 1}; \
	/* calculate block sizes */ \
	for (int i = 0; i < ga_handles[g_a].ndims; i++) \
	{ \
	    ga_local_dims[i] =  ((ga_handles[g_a].dims[i] / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].dims[i] / 64)); \
	    ga_local_dims[i] = (( ga_local_dims[i] >= ga_handles[g_a].dims[i] ) ? (ga_handles[g_a].dims[i] / 2) : ga_local_dims[i] ); \
	} \
	/* initialize matrix (block, block) to hold a portion of GA */ \
	Matrix < T >A (ga_local_dims[0], ga_local_dims[1]); \
	Zeros (A, ga_local_dims[0], ga_local_dims[1]); \
	/* initialize fetch-and-op window */ \
	mpi::WindowAllocate (sizeof (long), ga_handles[g_a].comm.comm, fop_win); \
	memset (fop_win_base, 0, sizeof (long)); \
	mpi::WindowLock( fop_win ); \
	long counter = 0, next = 0; \
	next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	/* initialize GA dim and strides */ \
	Int ga_height = ga_handles[g_a].dims[1]; \
	Int ga_width = ga_handles[g_a].dims[0]; \
	Int height_stride = ga_local_dims[1]; \
	Int width_stride = ga_local_dims[0]; \
	for (int i = 0; i < ga_height; i += height_stride) \
	{ \
	    if (counter == next) \
	    { \
		for (int j = 0; j < ga_width; j += width_stride) \
		{ \
		    ga_handles[g_a].rmaint.Getx( scalar, A, i, j ); \
		    ga_handles[g_a].rmaint.Iput( A, i, j ); \
		} \
		next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	    } \
	    counter++; \
	} \
	ga_handles[g_a].rmaint.Flush( A ); \
	mpi::Barrier( ga_handles[g_a].comm ); \
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
	int ga_local_dims[2] = {1, 1}; \
	for (int i = 0; i < ga_handles[g_a].ndims; i++) \
	{ \
	    ga_local_dims[i] =  ((ga_handles[g_a].dims[i] / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].dims[i] / 64)); \
	    ga_local_dims[i] = (( ga_local_dims[i] >= ga_handles[g_a].dims[i] ) ? (ga_handles[g_a].dims[i] / 2) : ga_local_dims[i] ); \
	} \
	Matrix < T >A (ga_local_dims[0], ga_local_dims[1]); \
	Matrix < T >B (ga_local_dims[0], ga_local_dims[1]); \
	Zeros (A, ga_local_dims[0], ga_local_dims[1]); \
	Zeros (B, ga_local_dims[0], ga_local_dims[1]); \
	mpi::WindowAllocate (sizeof (long), ga_handles[g_a].comm.comm, fop_win); \
	memset (fop_win_base, 0, sizeof (long)); \
	mpi::WindowLock( fop_win ); \
	long counter = 0, next = 0; \
	mpi::Barrier( ga_handles[g_a].comm ); \
	next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	for (int i = 0; i < ga_handles[g_a].dims[1]; i += ga_local_dims[1]) \
	{ \
	    if (counter == next) \
	    { \
		for (int j = 0; j < ga_handles[g_a].dims[0]; j += ga_local_dims[0]) \
		{ \
		    ga_handles[g_a].rmaint.Getx( alpha, A, i, j ); \
		    ga_handles[g_b].rmaint.Getx( beta, B, i, j ); \
		    ga_handles[g_c].rmaint.Iput( A, i, j ); \
		    ga_handles[g_c].rmaint.Iacc( B, i, j ); \
		} \
		next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	    } \
	    counter++; \
	} \
	ga_handles[g_a].rmaint.Flush( A ); \
	ga_handles[g_b].rmaint.Flush( B ); \
	ga_handles[g_c].rmaint.Flush( A ); \
	ga_handles[g_c].rmaint.Flush( B ); \
	mpi::Barrier( ga_handles[g_a].comm ); \
	mpi::WindowUnlock(fop_win ); \
	mpi::WindowFree(fop_win ); \
    } while (0)

// copy g_a into g_b
#define GA_COPY(g_a, g_b) \
    do { \
	mpi::Window fop_win; \
	long *fop_win_base; \
	int ga_local_dims[2] = {1, 1}; \
	for (int i = 0; i < ga_handles[g_a].ndims; i++) \
	{ \
	    ga_local_dims[i] =  ((ga_handles[g_a].dims[i] / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].dims[i] / 64)); \
	    ga_local_dims[i] = (( ga_local_dims[i] >= ga_handles[g_a].dims[i] ) ? (ga_handles[g_a].dims[i] / 2) : ga_local_dims[i] ); \
	} \
	Matrix < T >A (ga_local_dims[0], ga_local_dims[1]); \
	Zeros (A, ga_local_dims[0], ga_local_dims[1]); \
	mpi::WindowAllocate (sizeof (long), ga_handles[g_a].comm.comm, fop_win); \
	memset (fop_win_base, 0, sizeof (long));\
	mpi::WindowLock( fop_win ); \
	long counter = 0, next = 0; \
	mpi::Barrier( ga_handles[g_a].comm ); \
	next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	for (int i = 0; i < ga_handles[g_a].dims[1]; i += ga_local_dims[1]) \
	{ \
	    if (counter == next) \
	    { \
		for (int j = 0; j < ga_handles[g_a].dims[0]; j += ga_local_dims[0]) \
		{ \
		    ga_handles[g_a].rmaint.Get( A, i, j ); \
		    ga_handles[g_b].rmaint.Iput( A, i, j ); \
		} \
		next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	    } \
	    counter++; \
	} \
	ga_handles[g_a].rmaint.Flush( A ); \
	ga_handles[g_b].rmaint.Flush( A ); \
	mpi::Barrier( ga_handles[g_a].comm ); \
	mpi::WindowUnlock(fop_win ); \
	mpi::WindowFree(fop_win ); \
    } while (0)


// g_c = alpha * g_a  +  beta * g_b;
template<typename T>
void GlobalArrays< T >::GA_Add(void *alpha, int g_a, void* beta, int g_b, int g_c)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Add" ) )
    
    T * a = (T *)alpha;
    T * b = (T *)beta;
    
    // arrays must have the same shape
    for (int i = 0; i < ga_handles[g_a].ndims; i++)
    {
	if ((ga_handles[g_a].dims[i] != ga_handles[g_b].dims[i]) 
		|| (ga_handles[g_a].dims[i] != ga_handles[g_c].dims[i]) 
		|| (ga_handles[g_b].dims[i] != ga_handles[g_c].dims[i]))
	    LogicError ("Global Arrays of different shapes cannot be added. ");
    }

    // add
    GA_ADD (*a, *b, g_a, g_b, g_c);
}

// copies g_a into g_b, must be of same shape
template<typename T>
void GlobalArrays< T >::GA_Copy(int g_a, int g_b)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Allocate" ) )

    // arrays must have the same shape
    for (int i = 0; i < ga_handles[g_a].ndims; i++)
    {
	if (ga_handles[g_a].dims[i] != ga_handles[g_b].dims[i])
	    LogicError ("Global Arrays of different shapes cannot be copied. ");
    }

    // copy
    // GA_COPY (g_a, g_b)
    Copy (ga_handles[g_a].DM, ga_handles[g_b].DM);
}

// deallocates ga and frees associated resources
template<typename T>
void GlobalArrays< T >::GA_Destroy(int g_a)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Destroy" ) )
	
    // detach RmaInterface
    ga_handles[g_a].rmaint.Detach();
    // erase
    ga_handles.erase( ga_handles.begin() + g_a );
}

// A := 0.5 * ( A + A' )
template<typename T>
void GlobalArrays< T >::GA_Symmetrize (int g_a)
{
    // create an empty GA for transpose
    int g_at = GA_Duplicate( g_a, "transpose array" ); 
    int g_c = GA_Duplicate( g_a, "final array" ); 

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
void GlobalArrays< T >::GA_Dgemm(char ta, char tb, int m, int n, int k, double alpha, int g_a, int g_b, double beta, int g_c )
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Dgemm" ) )
    
    // set algorithmic block size	
    int nb = 96;
    GemmAlgorithm alg = GEMM_SUMMA_A;
    
    const Orientation orientA = CharToOrientation( ((ta == 'T' || ta == 't') ? 'T' : 'N') );
    const Orientation orientB = CharToOrientation( ((ta == 'T' || ta == 't') ? 'T' : 'N') );
    SetBlocksize( nb );

    T a = static_cast<T>( alpha );
    T b = static_cast<T>( beta );

    Gemm( orientA, orientB, a, ga_handles[g_a].DM, ga_handles[g_b].DM, 
		b, ga_handles[g_c].DM, alg);
}

#define FILL_GA(value, g_a) \
    do { \
	mpi::Window fop_win; \
	long *fop_win_base; \
	int ga_local_dims[2] = {1, 1}; \
	for (int i = 0; i < ga_handles[g_a].ndims; i++) \
	{ \
	    ga_local_dims[i] =  ((ga_handles[g_a].dims[i] / 64) == 0 ? FALLBACK_STRIDE : (ga_handles[g_a].dims[i] / 64)); \
	    ga_local_dims[i] = (( ga_local_dims[i] >= ga_handles[g_a].dims[i] ) ? (ga_handles[g_a].dims[i] / 2) : ga_local_dims[i] ); \
	} \
	/* initialize matrix (block, block) as a portion of GA */ \
	Int local_height = ga_local_dims[1]; \
	Int local_width = ga_local_dims[0]; \
	Matrix < T >A (local_width, local_height); \
	T* buffer = A.Buffer(); \
	const int ldim = A.LDim(); \
	/* fill local matrix with particular value */ \
	for( int j = 0; j < local_height; j++ ) \
	{ \
	    for( int i = 0; i < local_width; i++ ) \
	    { \
		buffer[i + j * ldim] = value; \
	    } \
	} \
	mpi::WindowAllocate (sizeof (long), ga_handles[g_a].comm.comm, fop_win); \
	memset (fop_win_base, 0, sizeof (long));\
	mpi::WindowLock( fop_win ); \
	long counter = 0, next = 0; \
	mpi::Barrier( ga_handles[g_a].comm ); \
	next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	/* update GA with chunks initialized with particular value */ \
	Int ga_height = ga_handles[g_a].dims[1]; \
	Int ga_width = ga_handles[g_a].dims[0]; \
	for (int i = 0; i < ga_height; i += local_height) \
	{ \
	    if (counter == next) \
	    { \
		for (int j = 0; j < ga_width; j += local_width) \
		{ \
		    ga_handles[g_a].rmaint.Iput( A, i, j ); \
		} \
		next = mpi::ReadInc( fop_win, 0, (long) 1 ); \
	    } \
	    counter++; \
	} \
	ga_handles[g_a].rmaint.Flush( A ); \
	mpi::Barrier( ga_handles[g_a].comm ); \
	mpi::WindowUnlock(fop_win ); \
	mpi::WindowFree(fop_win ); \
    } while (0)

// fills a global array with a specific value
template<typename T>
void GlobalArrays< T >::GA_Fill(int g_a, void *value)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Fill" ) )
    
    T * v = (T *)(value);
    if (*v != 0)
	FILL_GA (*v, g_a);
    else
	Zeros( ga_handles[g_a].DM, ga_handles[g_a].dims[0], ga_handles[g_a].dims[1] );
}

// allocate and initialize internal data structures in Global Arrays.
template<typename T>
void GlobalArrays< T >::GA_Initialize()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Initialize" ) )

   if (ga_handles.size() == 0)
	ga_initialized = true;
}

// barrier
// TODO dont sync over comm world
// fetch comm from ga handles
template<typename T>
void GlobalArrays< T >::GA_Sync()
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Sync" ) )
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
	ga_handles[i].rmaint.Detach();
    ga_handles.clear();
}

// B = A'
template<typename T>
void GlobalArrays< T >::GA_Transpose(int g_a, int g_b)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::GA_Transpose" ) )

    Transpose( ga_handles[g_a].DM, ga_handles[g_b].DM );
}

// Inquires for the data range on a specified processor
template<typename T>
void GlobalArrays< T >::NGA_Distribution(int g_a, int iproc, int lo[], int hi[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Distribution" ) )

    int dim[2] = {-1, -1};
    Grid grid (ga_handles[g_a].comm);
    const int my_rank = grid.VCRank();

    // find the width and height of the submatrix held by process iproc
    dim[0] = ga_handles[g_a].DM.LocalWidth();
    dim[1] = ga_handles[g_a].DM.LocalHeight();

    // broadcast iproc's dim to everyone
    if (iproc != my_rank)
	mpi::Broadcast(dim, 2, iproc, ga_handles[g_a].comm);

    if (dim[0] == 0 && dim[1] == 0) // in case iproc does not own a submatrix
    {
	lo[0] = -1; lo[1] = -1;
	hi[0] = -2; hi[1] = -2;
    }
    else
    {
	lo[0] = 0; lo[1] = 0;
	hi[0] = dim[0]; hi[1] = dim[1];
    }
}

// accesses data locally allocated for a global array    
template<typename T>
void GlobalArrays< T >::NGA_Access(int g_a, int lo[], int hi[], void *ptr, int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Access" ) )

    int width, height;
    // calculate height and width from lo and hi if possible
    if (lo[0] != -1 && hi[0] != -2)
    {
	width = hi[0] - lo[0];
	height = hi[1] - lo[1];
    }
    else // full local dimensions
    {
	width = ga_handles[g_a].DM.LocalWidth();
	height = ga_handles[g_a].DM.LocalHeight();
    } 

    T * buffer = (T *)ptr;
    T * Abuf = (T *)ga_handles[g_a].DM.Buffer();
	
    memcpy (buffer, Abuf, (height * width * sizeof(T)));
}

// transfers
// locally blocking transfers
#define BXFER(str, g_a, M, i, j) \
    do { \
	switch (str) \
	{ \
	    case 'A': \
		ga_handles[g_a].rmaint.Acc (M, i, j); \
	        break; \
	    case 'P': \
		ga_handles[g_a].rmaint.Put (M, i, j); \
	        break; \
	    case 'G': \
		ga_handles[g_a].rmaint.Get (M, i, j); \
	    	break; \
	} \
    } while (0)

// locally nonblocking transfers
#define NBXFER(str, g_a, M, i, j) \
    do { \
	switch (str) \
	{ \
	    case 'A': \
		ga_handles[g_a].rmaint.Iacc (M, i, j); \
	        break; \
	    case 'P': \
		ga_handles[g_a].rmaint.Iput (M, i, j); \
	        break; \
	    case 'G': \
		ga_handles[g_a].rmaint.Get (M, i, j); \
	    	break; \
	} \
    } while (0)

template<typename T>
void  GlobalArrays< T >::NGA_Acc(int g_a, int lo[], int hi[], void* buf, int ld[], void* alpha)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Acc" ) )

    T * a = (T *)(alpha);
    
    if (*a == 0.0) // no-op
	return;

    T * buffer = (T *) buf;

    // calculate height and width from lo and hi
    int width = hi[0] - lo[0];
    int height = hi[1] - lo[1];
    int ldim = width;

    if (*a != 1.0) // then a * buf
    {
	for (int i = 0; i < height; i++)
	    for (int j = 0; j < width; j++)
	    buffer[i * width + j] *= *a;
    }

    // unpack buf into Matrix<T>
    Matrix< T > A;
    A.Attach(height, width, buffer, ldim);

    // Acc - (locally) blocking transfer
    BXFER ('A', g_a, A, lo[0], hi[0]);
    ga_handles[g_a].rmaint.Flush (A);
}

template<typename T>
void GlobalArrays< T >::NGA_Get(int g_a, int lo[], int hi[], void* buf, int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Get" ) )

     // calculate height and width from lo and hi
    int width = hi[0] - lo[0];
    int height = hi[1] - lo[1];
    
    // unpack buf into Matrix<T>
    Matrix< T > A (height, width);
    BXFER ('G', g_a, A, lo[0], hi[0]);

    T * buffer = A.Buffer();
    T * inbuf = (T *)buf;

    const int ldim = A.LDim();

    for( int j = 0; j < height; j++ )
    for( int i = 0; i < width; i++ )
	inbuf[i+j*ldim] = buffer[i+j*ldim];
}
 
template<typename T>
void GlobalArrays< T >::NGA_NbAcc(int g_a,int lo[], int hi[],void* buf,int ld[],void* alpha, ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbAcc" ) )

    T * a = (T *)alpha;
    T * buffer = (T *) buf;

    if (*a == 0.0) // no-op
	return;

    // calculate height and width from lo and hi
    int width = hi[0] - lo[0];
    int height = hi[1] - lo[1];
    int ldim = width;

    if (*a != 1.0) // then a * buf
    {
	for (int i = 0; i < height; i++)
	    for (int j = 0; j < width; j++)
	    buffer[i * width + j] *= *a;
    }

    // unpack buf into Matrix<T>
    Matrix< T > A;
    A.Attach(height, width, buffer, ldim);

    // Acc - (locally) nonblocking transfer
    NBXFER ('A', g_a, A, lo[0], hi[0]);

    *nbhandle = ga_handles[g_a].handle;
    if ( !ga_handles[g_a].pending_transfer ) 
	ga_handles[g_a].pending_transfer = true; 
}

template<typename T>
void GlobalArrays< T >::NGA_NbGet(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle)
{
    NGA_Get (g_a, lo, hi, buf, ld);
    *nbhandle = -1; // no nb get implementation in el::rmainterface
}

template<typename T>
void GlobalArrays< T >::NGA_NbPut(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbPut" ) )

    // calculate height and width from lo and hi
    int width = hi[0] - lo[0];
    int height = hi[1] - lo[1];
    const int ldim = width;

    T * buffer = (T *) buf;

    // unpack buf into Matrix<T>
    Matrix< T > A;
    A.Attach(height, width, buffer, ldim);

    // Put - (locally) nonblocking transfer
    NBXFER ('P', g_a, A, lo[0], hi[0]);

    *nbhandle = ga_handles[g_a].handle;
    if ( !ga_handles[g_a].pending_transfer ) 
	ga_handles[g_a].pending_transfer = true; 
}
 
template<typename T>
int GlobalArrays< T >::NGA_NbTest(ga_nbhdl_t* nbhandle)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_NbTest" ) )

    if (*nbhandle != -1)
    {
	bool status = ga_handles[*nbhandle].rmaint.Testall ();
	if (status) // release handle
	{
	    ga_handles[*nbhandle].pending_transfer = false; 
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
    
    if (*nbhandle != -1)
    {
	ga_handles[*nbhandle].rmaint.Waitall ();
	ga_handles[*nbhandle].pending_transfer = false; 
	// release handle
	*nbhandle = -1;
    }
}

template<typename T>
void GlobalArrays< T >::NGA_Put(int g_a, int lo[], int hi[], void* buf, int ld[])
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Put" ) )

    // calculate height and width from lo and hi
    int width = hi[0] - lo[0];
    int height = hi[1] - lo[1];
    const int ldim = width;

    T * buffer = (T *) buf;

    // unpack buf into Matrix<T>
    Matrix< T > A;
    A.Attach(height, width, buffer, ldim);

    // Put - (locally) nonblocking transfer
    BXFER('P', g_a, A, lo[0], hi[0]);
    ga_handles[g_a].rmaint.Flush (A);
}

template<typename T>
long GlobalArrays< T >::NGA_Read_inc(int g_a, int subscript[], long inc)
{
    DEBUG_ONLY( CallStackEntry cse( "GlobalArrays::NGA_Read_inc" ) )

    long prev;

    if (ga_handles[g_a].ndims == 2)
    {
	int i = subscript[0];
	int j = subscript[1];
	prev = ga_handles[g_a].rmaint.CompareAndSwap (i, j, inc);
    }
    else // 1 dimension
    {
	int j = subscript[0];
	prev = ga_handles[g_a].rmaint.CompareAndSwap (1, j, inc);
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
