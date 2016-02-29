/*
   Copyright (c) 2009-2015, Jack Poulson
                      2013, Jeff Hammond
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_MPI_HPP
#define EL_IMPORTS_MPI_HPP
#include <limits>

namespace El {

namespace mpi {

#if defined(EL_HAVE_MPI3_NONBLOCKING_COLLECTIVES) || \
    defined(EL_HAVE_MPIX_NONBLOCKING_COLLECTIVES)
#define EL_HAVE_NONBLOCKING 1
#else
#define EL_HAVE_NONBLOCKING 0
#endif

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
#ifdef EL_HAVE_MPI3_NONBLOCKING_COLLECTIVES
#define EL_NONBLOCKING_COLL(name) MPI_ ## name
#else
#define EL_NONBLOCKING_COLL(name) MPIX_ ## name
#endif
#endif

struct Comm
{
    MPI_Comm comm;
    Comm( MPI_Comm mpiComm=MPI_COMM_NULL ) : comm(mpiComm) { }
};
inline bool operator==( const Comm& a, const Comm& b )
{ return a.comm == b.comm; }
inline bool operator!=( const Comm& a, const Comm& b )
{ return a.comm != b.comm; }

struct Group
{
    MPI_Group group;
    Group( MPI_Group mpiGroup=MPI_GROUP_NULL ) : group(mpiGroup) { }
};
inline bool operator==( const Group& a, const Group& b )
{ return a.group == b.group; }
inline bool operator!=( const Group& a, const Group& b )
{ return a.group != b.group; }

struct Op
{
    MPI_Op op;
    Op( MPI_Op mpiOp=MPI_OP_NULL ) : op(mpiOp) { }
};
inline bool operator==( const Op& a, const Op& b )
{ return a.op == b.op; }
inline bool operator!=( const Op& a, const Op& b )
{ return a.op != b.op; }

// Datatype definitions
// TODO: Convert these to structs/classes
typedef MPI_Aint Aint;
typedef MPI_Datatype Datatype;
typedef MPI_Errhandler ErrorHandler;
typedef MPI_Request Request;
typedef MPI_Status Status;
typedef MPI_Message Message;
typedef MPI_User_function UserFunction;
#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY)
typedef MPI_Win Window;
const Window WIN_NULL = MPI_WIN_NULL;
typedef enum
{
	STRICT_ACC_ORDERING 	= 0,
	PARTIAL_ACC_ORDERING 	= 2,
	NO_ACC_ORDERING 	= 4
} acc_order_t;
#endif

typedef MPI_Info Info;
// Standard constants
const int ANY_SOURCE = MPI_ANY_SOURCE;
const int ANY_TAG = MPI_ANY_TAG;
const int ERR_RANK = MPI_ERR_RANK;
#ifdef EL_HAVE_MPI_QUERY_THREAD
const int THREAD_SINGLE = MPI_THREAD_SINGLE;
const int THREAD_FUNNELED = MPI_THREAD_FUNNELED;
const int THREAD_SERIALIZED = MPI_THREAD_SERIALIZED;
const int THREAD_MULTIPLE = MPI_THREAD_MULTIPLE;
#else
const int THREAD_SINGLE = 0;
const int THREAD_FUNNELED = 1;
const int THREAD_SERIALIZED = 2;
const int THREAD_MULTIPLE = 3;
#endif
const int DATATYPE_NULL = MPI_DATATYPE_NULL;
const int UNDEFINED = MPI_UNDEFINED;
const Group GROUP_NULL = MPI_GROUP_NULL;
const Comm COMM_NULL = MPI_COMM_NULL;
const Comm COMM_SELF = MPI_COMM_SELF;
const Comm COMM_WORLD = MPI_COMM_WORLD;
const ErrorHandler ERRORS_RETURN = MPI_ERRORS_RETURN;
const ErrorHandler ERRORS_ARE_FATAL = MPI_ERRORS_ARE_FATAL;
const Group GROUP_EMPTY = MPI_GROUP_EMPTY;
const Request REQUEST_NULL = MPI_REQUEST_NULL;

const Op MAX = MPI_MAX;
const Op MIN = MPI_MIN;
const Op MAXLOC = MPI_MAXLOC;
const Op MINLOC = MPI_MINLOC;
const Op PROD = MPI_PROD;
const Op SUM = MPI_SUM;
const Op REPLACE = MPI_REPLACE;
const Op LOGICAL_AND = MPI_LAND;
const Op LOGICAL_OR = MPI_LOR;
const Op LOGICAL_XOR = MPI_LXOR;
const Op BINARY_AND = MPI_BAND;
const Op BINARY_OR = MPI_BOR;
const Op BINARY_XOR = MPI_BXOR;

// Added constant(s)
const int MIN_COLL_MSG = 1; // minimum message size for collectives
inline int Pad( int count ) { return std::max(count,MIN_COLL_MSG); }

bool CommSameSizeAsInteger();
bool GroupSameSizeAsInteger();

// Environment routines
void Initialize( int& argc, char**& argv );
int InitializeThread( int& argc, char**& argv, int required );
void Finalize();
bool Initialized();
bool Finalized();
int QueryThread();
void Abort( Comm comm, int errCode );
double Time();
void Create( UserFunction* func, bool commutes, Op& op );
void Free( Op& op );

// Communicator manipulation
int WorldRank();
int Rank( Comm comm );
int Size( Comm comm );
void Create( Comm parentComm, Group subsetGroup, Comm& subsetComm );
void Dup( Comm original, Comm& duplicate );
void Split( Comm comm, int color, int key, Comm& newComm );
void Free( Comm& comm );
bool Congruent( Comm comm1, Comm comm2 );
void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler );

// Cartesian communicator routines
void CartCreate
( Comm comm, int numDims, const int* dimensions, const int* periods, 
  bool reorder, Comm& cartComm );
void CartSub
( Comm comm, const int* remainingDims, Comm& subComm );

// Group manipulation
int Rank( Group group );
int Size( Group group );
void CommGroup( Comm comm, Group& group );
void Dup( Group group, Group& newGroup );
void Union( Group groupA, Group groupB, Group& newGroup );
void Incl( Group group, int n, const int* ranks, Group& subGroup );
void Excl( Group group, int n, const int* ranks, Group& subGroup );
void Difference( Group parent, Group subset, Group& complement );
void Free( Group& group );
int Translate( Group origGroup, int origRank, Group newGroup );
int Translate( Comm  origComm,  int origRank, Group newGroup );
int Translate( Group origGroup, int origRank, Comm  newComm  );
int Translate( Comm  origComm,  int origRank, Comm  newComm  );
void Translate
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks );
void Translate
( Comm origComm,  int size, const int* origRanks, 
  Group newGroup,                 int* newRanks );
void Translate
( Group origGroup, int size, const int* origRanks, 
  Comm newComm,                    int* newRanks );
void Translate
( Comm origComm, int size, const int* origRanks, 
  Comm newComm,                  int* newRanks );

// MPI-3 one-sided
// ===============
#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY)

// Window creation/update/delete
// -----------------------------
void WindowLock( int rank, Window & window );
void WindowLock( Window & window );
void WindowUnlock( int rank, Window & window );
void WindowUnlock( Window & window );
// --------
template<typename R>
void WindowCreate( R * baseptr, Aint entries, Comm comm, Window & window );
template<typename R>
void WindowCreate( Complex<R> * baseptr, Aint entries, Comm comm, Window & window );
template<typename R>
void WindowAllocate( Aint entries, Comm comm, Window & window );   
void * GetWindowBase( Window & window );
// --------
void WindowFree( Window & window );
void SetWindowProp( Window & window, acc_order_t prop );

// TODO come up with a generic
// interface for UDD
// create user-defined data type 
// -----------------------------
template<typename R>
void CreateContigType( int count, Datatype * contig_type );	    
template<typename R>
void CreateVectorType( int count, int blocklength, int stride, 
	Datatype * vector_type );    
void DestroyType( Datatype * subtype );

// Atomic operations
// -----------------
template<typename R>
R ReadInc (Window & win, Aint offset, 
	R inc, int fop_root = 0);
template<typename R>
Complex<R> ReadInc (Window & win, Aint offset, 
	Complex<R> inc, int fop_root = 0);

// One-sided operations
// --------------------
// put
// ---
template<typename R>
void Iput (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window);
template<typename R>
void Iput (const Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window);
template<typename R>
void Iput (const R* source, int target_rank, Aint disp, 
	Datatype & dtype, Window & window);
template<typename R>
void Iput (const Complex<R>* source, int target_rank,
           Aint disp, Datatype & dtype, Window & window);
template<typename R>
void Iput( const R * source, int origin_count, int target_rank, Aint disp, 
	int target_count, Datatype & target_type, Window & window );
template<typename R>
void Iput( const Complex<R>* source, int origin_count, int target_rank, Aint disp, 
	int target_count, Datatype & target_type, Window & window );

// get
// ---
template<typename R>
void Iget (R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window);
template<typename R>
void Iget (Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window);
template<typename R>
void Iget (R* source, int target_rank, Aint disp, 
	Datatype & dtype, Window & window);
template<typename R>
void Iget (Complex<R>* source, int target_rank, Aint disp, 
	Datatype & dtype, Window & window);
template<typename R>
void Iget( R * source, int origin_count, int target_rank, Aint disp, 
	int target_count, Datatype & target_type, Window & window );
template<typename R>
void Iget( Complex<R>* source, int origin_count, int target_rank, Aint disp, 
	int target_count, Datatype & target_type, Window & window );

// acc
// ---
template<typename R>
void Iacc (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Op op, Window & window);
template<typename R>
void Iacc (const Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Op op, Window & window);
// op == SUM
template<typename R>
void Iacc (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window);
template<typename R>
void Iacc (const Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window);
// derived data type
template<typename R>
void Iacc( const R * source, int origin_count, int target_rank, Aint disp, 
	int target_count, Datatype & target_type, Op op, Window & window );
template<typename R>
void Iacc( const Complex<R>* source, int origin_count, int target_rank, Aint disp, 
	int target_count, Datatype & target_type, Op op, Window & window );

// Synchronization
// ---------------
void Flush( int target_rank, Window& window );
void Flush( Window& window );
void FlushLocal( int target_rank, Window& window );
void FlushLocal( Window& window );
void WindowSync( Window& window );

// Collectives
// -----------
// In-place AllGather
// ------------------
template<typename R>
void AllGather
( R* rbuf, int rc, Comm comm );
template<typename R>
void AllGather
( Complex<R>* rbuf, int rc, Comm comm );

#endif // EL_ENABLE_RMA_AXPY

// Utilities
void Barrier( Comm comm );
#if MPI_VERSION>=3 && defined(EL_USE_IBARRIER_FOR_AXPY)
void IBarrier( Comm comm, Request& request );
#endif
void RequestFree( Request& request );
void Wait( Request& request );
void Wait( Request& request, Status& status );
//TODO add another function for getting statuses
void WaitAny (int numRequests, Request * requests, Int * index);
void WaitAll( int numRequests, Request* requests );
void WaitAll( int numRequests, Request* requests, Status* statuses );
bool Test( Request& request );
bool Test( Request& request, Status& status );
bool Testany( int count, Request* requests );
bool Testany( int count, Request* requests, int& indx );
bool Testany( int count, Request* requests, int& indx, Status& status );
bool IProbe( int source, int tag, Comm comm, Status& status );
bool IProbe( int source, Comm comm, Status& status );
bool IProbe( Comm comm, Status& status );
void Probe ( int source, int tag, Comm comm, Status & status );
void Probe ( int source, Comm comm, Status & status );
void Probe ( Comm comm, Status & status );
// matching probe
bool IMprobe( int source, int tag, Comm comm, Status& status, Message& message );

template<typename T>
int GetCount( Status& status );

// Point-to-point communication
// ============================

// Send
// ----
template<typename R>
void TaggedSend( const R* buf, int count, int to, int tag, Comm comm );
template<typename R>
void TaggedSend( const Complex<R>* buf, int count, int to, int tag, Comm comm );
// If the tag is irrelevant
template<typename T>
void Send( const T* buf, int count, int to, Comm comm );
// If the send-count is one
template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm );
// If the send-count is one and the tag is irrelevant
template<typename T>
void Send( T b, int to, Comm comm );

// Non-blocking send
// -----------------
template<typename R>
void TaggedISend
( const R* buf, int count, int to, int tag, Comm comm, Request& request );
template<typename R>
void TaggedISend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request );
// If the tag is irrelevant
template<typename T>
void ISend( const T* buf, int count, int to, Comm comm, Request& request );
// If the send count is one
template<typename T>
void TaggedISend( T b, int to, int tag, Comm comm, Request& request );
// If the send count is one and the tag is irrelevant
template<typename T>
void ISend( T b, int to, Comm comm, Request& request );

// Non-blocking synchronous Send
// -----------------------------
template<typename R>
void TaggedISSend
( const R* buf, int count, int to, int tag, Comm comm, Request& request );
template<typename R>
void TaggedISSend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request );
// If the tag is irrelevant
template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request& request );
// If the send count is one
template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request& request );
// If the send count is one and the tag is irrelevant
template<typename T>
void ISSend( T b, int to, Comm comm, Request& request );

// Recv
// ----
template<typename R>
void TaggedRecv( R* buf, int count, int from, int tag, Comm comm );
template<typename R>
void TaggedRecv( Complex<R>* buf, int count, int from, int tag, Comm comm );
// If the tag is irrelevant
template<typename T>
void Recv( T* buf, int count, int from, Comm comm );
// If the recv count is one
template<typename T>
T TaggedRecv( int from, int tag, Comm comm );
// If the recv count is one and the tag is irrelevant
template<typename T>
T Recv( int from, Comm comm );
// matched recv
template<typename R>
void TaggedMrecv( R* buf, int count, Message & message );
template<typename R>
void TaggedMrecv( Complex<R>* buf, int count, Message & message );
template<typename R>
void TaggedRecvS( R* buf, int count, int from, int tag, Comm comm, Status & status );
template<typename R>
void TaggedRecvS( Complex<R>* buf, int count, int from, int tag, Comm comm, Status & status );

// Non-blocking recv
// -----------------
template<typename R>
void TaggedIRecv
( R* buf, int count, int from, int tag, Comm comm, Request& request );
template<typename R>
void TaggedIRecv
( Complex<R>* buf, int count, int from, int tag, Comm comm, Request& request );
// If the tag is irrelevant
template<typename T>
void IRecv( T* buf, int count, int from, Comm comm, Request& request );
// If the recv count is one
template<typename T>
T TaggedIRecv( int from, int tag, Comm comm, Request& request );
// If the recv count is one and the tag is irrelevant
template<typename T>
T IRecv( int from, Comm comm, Request& request );

// SendRecv
// --------
template<typename R>
void TaggedSendRecv
( const R* sbuf, int sc, int to,   int stag,
        R* rbuf, int rc, int from, int rtag, Comm comm );
template<typename R>
void TaggedSendRecv
( const Complex<R>* sbuf, int sc, int to,   int stag,
        Complex<R>* rbuf, int rc, int from, int rtag, Comm comm );
// If the tags are irrelevant
template<typename T>
void SendRecv
( const T* sbuf, int sc, int to,
        T* rbuf, int rc, int from, Comm comm );
// If the send and recv counts are one
template<typename T>
T TaggedSendRecv( T sb, int to, int stag, int from, int rtag, Comm comm );
// If the send and recv counts are one and the tags don't matter
template<typename T>
T SendRecv( T sb, int to, int from, Comm comm );

// Single-buffer SendRecv
// ----------------------
template<typename R>
void TaggedSendRecv
( R* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template<typename R>
void TaggedSendRecv
( Complex<R>* buf, int count, int to, int stag, int from, int rtag, Comm comm );
// If the tags don't matter
template<typename T>
void SendRecv( T* buf, int count, int to, int from, Comm comm );

// Collective communication
// ========================

// Broadcast
// ---------
template<typename R>
void Broadcast( R* buf, int count, int root, Comm comm );
template<typename R>
void Broadcast( Complex<R>* buf, int count, int root, Comm comm );
// If the message length is one
template<typename T>
void Broadcast( T& b, int root, Comm comm );

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
// Non-blocking broadcast
// ----------------------
template<typename R>
void IBroadcast
( R* buf, int count, int root, Comm comm, Request& request );
template<typename R>
void IBroadcast
( Complex<R>* buf, int count, int root, Comm comm, Request& request );
// If the message length is one
template<typename T>
void IBroadcast( T& b, int root, Comm comm, Request& request );
#endif

// Gather
// ------
template<typename R>
void Gather
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm );
template<typename R>
void  Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm );

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
// Non-blocking gather
// -------------------
template<typename R>
void IGather
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm, Request& request );
template<typename R>
void IGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm, Request& request );
#endif

// Gather with variable recv sizes
// -------------------------------
template<typename R>
void Gather
( const R* sbuf, int sc,
        R* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template<typename R>
void Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );

// AllGather
// ---------
template<typename R>
void AllGather
( const R* sbuf, int sc,
        R* rbuf, int rc, Comm comm );
template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm );

// AllGather with variable recv sizes
// ----------------------------------
template<typename R>
void AllGather
( const R* sbuf, int sc,
        R* rbuf, const int* rcs, const int* rds, Comm comm );
template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm );

// Scatter
// -------
template<typename R>
void Scatter
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm );
template<typename R>
void Scatter
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm );
// In-place option
template<typename R>
void Scatter( R* buf, int sc, int rc, int root, Comm comm );
template<typename R>
void Scatter( Complex<R>* buf, int sc, int rc, int root, Comm comm );

// AllToAll
// --------
template<typename R>
void AllToAll
( const R* sbuf, int sc,
        R* rbuf, int rc, Comm comm );
template<typename R>
void AllToAll
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm );

// AllToAll with non-uniform send/recv sizes
// -----------------------------------------
template<typename R>
void AllToAll
( const R* sbuf, const int* scs, const int* sds,
        R* rbuf, const int* rcs, const int* rds, Comm comm );
template<typename R>
void AllToAll
( const Complex<R>* sbuf, const int* scs, const int* sds,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm );

// Reduce
// ------
template<typename T>
void Reduce
( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm );
template<typename R>
void Reduce
( const Complex<R>* sbuf, Complex<R>* rbuf, int count, Op op, 
  int root, Comm comm );
// Default to mpi::SUM
template<typename T>
void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm );
// With a message-size of one
template<typename T>
T Reduce( T sb, Op op, int root, Comm comm );
// With a message-size of one and default to mpi::SUM
template<typename T>
T Reduce( T sb, int root, Comm comm );

// Single-buffer reduce
// --------------------
template<typename T>
void Reduce( T* buf, int count, Op op, int root, Comm comm );
template<typename R>
void Reduce( Complex<R>* buf, int count, Op op, int root, Comm comm );
// Default to mpi::SUM
template<typename T>
void Reduce( T* buf, int count, int root, Comm comm );

// AllReduce
// ---------
template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Op op, Comm comm );
template<typename R>
void AllReduce
( const Complex<R>* sbuf, Complex<R>* rbuf, int count, Op op, Comm comm );
// Default to mpi::SUM
template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm );
// If the message-length is one
template<typename T>
T AllReduce( T sb, Op op, Comm comm );
// If the message-length is one (and default to mpi::SUM)
template<typename T>
T AllReduce( T sb, Comm comm );

// Single-buffer AllReduce
// -----------------------
template<typename T>
void AllReduce( T* buf, int count, Op op, Comm comm );
template<typename R>
void AllReduce( Complex<R>* buf, int count, Op op, Comm comm );
// Default to mpi::SUM
template<typename T>
void AllReduce( T* buf, int count, Comm comm );

// ReduceScatter
// -------------
template<typename R>
void ReduceScatter
( R* sbuf, R* rbuf, int rc, Op op, Comm comm );
template<typename R>
void ReduceScatter
( Complex<R>* sbuf, Complex<R>* rbuf, int rc, Op op, Comm comm );
// Default to mpi::SUM
template<typename T>
void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm );

// Single-buffer ReduceScatter
// ---------------------------
template<typename R>
void ReduceScatter( R* buf, int rc, Op op, Comm comm );
template<typename R>
void ReduceScatter( Complex<R>* buf, int rc, Op op, Comm comm );
// Default to mpi::SUM
template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm );

// Variable-length ReduceScatter
// -----------------------------
template<typename R>
void ReduceScatter
( const R* sbuf, R* rbuf, const int* rcs, Op op, Comm comm );
template<typename R>
void ReduceScatter
( const Complex<R>* sbuf, Complex<R>* rbuf, const int* rcs, Op op, Comm comm );
// Default to mpi::SUM
template<typename T>
void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Comm comm );

template<typename R>
void MaxLocFunc
( void* in, void* out, int* length, mpi::Datatype* datatype );
template<typename R>
void MaxLocPairFunc
( void* in, void* out, int* length, mpi::Datatype* datatype );

template<typename R>
void MinLocFunc
( void* in, void* out, int* length, mpi::Datatype* datatype );
template<typename R>
void MinLocPairFunc
( void* in, void* out, int* length, mpi::Datatype* datatype );

template<typename R> mpi::Datatype& ValueIntType();
template<> mpi::Datatype& ValueIntType<Int>();
template<> mpi::Datatype& ValueIntType<float>();
template<> mpi::Datatype& ValueIntType<double>();

template<typename R> mpi::Datatype& ValueIntPairType();
template<> mpi::Datatype& ValueIntPairType<Int>();
template<> mpi::Datatype& ValueIntPairType<float>();
template<> mpi::Datatype& ValueIntPairType<double>();

template<typename R> void CreateValueIntType();
template<typename R> void CreateValueIntPairType();
template<typename R> void DestroyValueIntType();
template<typename R> void DestroyValueIntPairType();

template<typename R> mpi::Op MaxLocOp();
template<> mpi::Op MaxLocOp<Int>();
template<> mpi::Op MaxLocOp<float>();
template<> mpi::Op MaxLocOp<double>();

template<typename R> mpi::Op MaxLocPairOp();
template<> mpi::Op MaxLocPairOp<Int>();
template<> mpi::Op MaxLocPairOp<float>();
template<> mpi::Op MaxLocPairOp<double>();

template<typename R> void CreateMaxLocOp();
template<> void CreateMaxLocOp<Int>();
template<> void CreateMaxLocOp<float>();
template<> void CreateMaxLocOp<double>();

template<typename R> void CreateMaxLocPairOp();
template<> void CreateMaxLocPairOp<Int>();
template<> void CreateMaxLocPairOp<float>();
template<> void CreateMaxLocPairOp<double>();

template<typename R> void DestroyMaxLocOp();
template<> void DestroyMaxLocOp<Int>();
template<> void DestroyMaxLocOp<float>();
template<> void DestroyMaxLocOp<double>();

template<typename R> void DestroyMaxLocPairOp();
template<> void DestroyMaxLocPairOp<Int>();
template<> void DestroyMaxLocPairOp<float>();
template<> void DestroyMaxLocPairOp<double>();

template<typename R> mpi::Op MinLocOp();
template<> mpi::Op MinLocOp<Int>();
template<> mpi::Op MinLocOp<float>();
template<> mpi::Op MinLocOp<double>();

template<typename R> mpi::Op MinLocPairOp();
template<> mpi::Op MinLocPairOp<Int>();
template<> mpi::Op MinLocPairOp<float>();
template<> mpi::Op MinLocPairOp<double>();

template<typename R> void CreateMinLocOp();
template<> void CreateMinLocOp<Int>();
template<> void CreateMinLocOp<float>();
template<> void CreateMinLocOp<double>();

template<typename R> void CreateMinLocPairOp();
template<> void CreateMinLocPairOp<Int>();
template<> void CreateMinLocPairOp<float>();
template<> void CreateMinLocPairOp<double>();

template<typename R> void DestroyMinLocOp();
template<> void DestroyMinLocOp<Int>();
template<> void DestroyMinLocOp<float>();
template<> void DestroyMinLocOp<double>();

template<typename R> void DestroyMinLocPairOp();
template<> void DestroyMinLocPairOp<Int>();
template<> void DestroyMinLocPairOp<float>();
template<> void DestroyMinLocPairOp<double>();

template<typename T> Datatype TypeMap();
template<> inline Datatype TypeMap<byte>() { return MPI_UNSIGNED_CHAR; }
template<> inline Datatype TypeMap<int>() { return MPI_INT; }
template<> inline Datatype TypeMap<unsigned>() { return MPI_UNSIGNED; }
template<> inline Datatype TypeMap<long int>() { return MPI_LONG_INT; }
template<> inline Datatype TypeMap<long unsigned>()
{ return MPI_UNSIGNED_LONG; }
template<> inline Datatype TypeMap<long long int>()
{
#ifdef EL_HAVE_MPI_LONG_LONG
    return MPI_LONG_LONG_INT;
#else
    throw std::runtime_error("MPI_LONG_LONG_INT does not exist");
    return 0;
#endif
}
template<>
inline Datatype TypeMap<unsigned long long>()
{
#ifdef EL_HAVE_MPI_LONG_LONG
    return MPI_UNSIGNED_LONG_LONG;
#else
    throw std::runtime_error("MPI_UNSIGNED_LONG_LONG does not exist");
    return 0;
#endif
}
template<> inline Datatype TypeMap<float>() { return MPI_FLOAT; }
template<> inline Datatype TypeMap<double>() { return MPI_DOUBLE; }
template<> inline Datatype TypeMap<Complex<float>>()
{ return MPI_COMPLEX; }
template<> inline Datatype TypeMap<Complex<double>>()
{ return MPI_DOUBLE_COMPLEX; }

template<> inline Datatype TypeMap<ValueInt<Int>>()
{ return ValueIntType<Int>(); }
template<> inline Datatype TypeMap<ValueInt<float>>()
{ return ValueIntType<float>(); }
template<> inline Datatype TypeMap<ValueInt<double>>()
{ return ValueIntType<double>(); }

template<> inline Datatype TypeMap<ValueIntPair<Int>>()
{ return ValueIntPairType<Int>(); }
template<> inline Datatype TypeMap<ValueIntPair<float>>()
{ return ValueIntPairType<float>(); }
template<> inline Datatype TypeMap<ValueIntPair<double>>()
{ return ValueIntPairType<double>(); }

template<typename T>
void SparseAllToAll
( const std::vector<T>& sendBuffer,
  const std::vector<int>& sendCounts, 
  const std::vector<int>& sendOffs,
        std::vector<T>& recvBuffer,
  const std::vector<int>& recvCounts, 
  const std::vector<int>& recvOffs,
        mpi::Comm comm );

void VerifySendsAndRecvs
( const std::vector<int>& sendCounts,
  const std::vector<int>& recvCounts, mpi::Comm comm );

} // mpi
} // elem

#endif // ifndef EL_IMPORTS_MPI_HPP
