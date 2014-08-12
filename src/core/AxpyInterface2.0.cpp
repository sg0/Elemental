/*
This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include <cassert>

// TODO Use DDT for put/get/acc when EL_USE_DERIVED_TYPE is defined
// TODO bring back const interfaces
// TODO localflush
// TODO error checks
namespace El
{
template<typename T>
AxpyInterface2<T>::AxpyInterface2()
    : GlobalArrayPut_(0), GlobalArrayGet_(0),
    putVectors_(0), getVectors_(0), dataRequests_(0), matrixBase_(0),
    dataRequestStatuses_(0), toBeAttachedForPut_(false), 
    toBeAttachedForGet_(false), attached_(false), 
    detached_(false)
{ }

template<typename T>
AxpyInterface2<T>::AxpyInterface2( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::AxpyInterface2"))

    attached_ 		= false;
    detached_ 		= false;
    toBeAttachedForGet_ = true;
    toBeAttachedForPut_ = true;
    GlobalArrayPut_ 	= &Z;
    GlobalArrayGet_ 	= &Z;

    const Int p = Z.Grid ().Size();
    putVectors_.resize( p );
    getVectors_.resize( p );
    dataRequests_.resize (p);
    dataRequestStatuses_.resize (p);    
    requestRequests_.resize (p);
    requestRequestStatuses_.resize (p);    
    matrixBase_.resize (p);
}

template<typename T>
AxpyInterface2<T>::~AxpyInterface2()
{
        if( std::uncaught_exception() )
        {
            std::ostringstream os;
            os << "Uncaught exception detected during AxpyInterface2 destructor "
            "that required a call to Detach. Instead of allowing for the "
            "possibility of Detach throwing another exception and "
            "resulting in a 'terminate', we instead immediately dump the "
            "call stack (if not in RELEASE mode) since the program will "
            "likely hang:" << std::endl;
            std::cerr << os.str();
            DEBUG_ONLY(DumpCallStack())
        }
        else
        {
            Detach();
        }
}

template<typename T>
void AxpyInterface2<T>::Attach( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Attach"))
    // attached_ will be only set in Attach
    // and only unset in Detach
    if (!attached_)
        attached_ = true;
    else
        LogicError("Must detach before reattaching.");

    // if DistMatrix is non-const, all one-sided
    // transfers -- put, get and acc are possible
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
    {
        GlobalArrayPut_ 	= &Z;
        toBeAttachedForPut_ 	= true;
        GlobalArrayGet_ 	= &Z;
        toBeAttachedForGet_ 	= true;
    }
    const Grid& g = Z.Grid();
    const Int p = g.Size ();

    if (putVectors_.size() != p)
    {
	getVectors_.resize( p );
	putVectors_.resize( p );
	dataRequests_.resize (p);
	dataRequestStatuses_.resize (p);
	matrixBase_.resize (p);
	requestRequests_.resize (p);
	requestRequestStatuses_.resize (p);    
    }
}

template<typename T>
Int AxpyInterface2<T>::NextIndex
( Int dataSize, std::deque <std::vector<T>> &dataVectors,
  std::deque <mpi::Request> &requests,
  std::deque <bool> &requestStatus,
  std::deque <T *> &matrixBase,
  T * base_address )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndex"))
    const Int Index = Int(requests.size ());
    
    dataVectors.resize (Index + 1);
    dataVectors[Index].resize (dataSize);
    requests.push_back (mpi::REQUEST_NULL);
    requestStatus.push_back (true);
    // stores Matrix base address by index
    matrixBase.push_back (base_address);
    
    return Index;
}

template<typename T>
void AxpyInterface2<T>::Put( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Put"))

    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
        LogicError("Global matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayPut_;
    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix");

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;

    const Int XLDim = Z.LDim();
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    const Int YLDim = Y.LDim ();

    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        // number of entries in my PE
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;
	
        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            T* XBuffer = Z.Buffer();
     	    const Int index =
	    NextIndex (numEntries, putVectors_[destination],
		    dataRequests_[destination],
		    dataRequestStatuses_[destination],
		    matrixBase_[destination],
		    XBuffer);

	    DEBUG_ONLY (if
			(Int (putVectors_[destination][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)
	    T *sendBuffer = putVectors_[destination][index].data ();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }
	    // put request 
	    mpi::TaggedISSend (sendBuffer, numEntries, destination, 
			DATA_PUT_TAG, g.VCComm (), 
			dataRequests_[destination][index]);
        }
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Get( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Get"))
	// a call to Attach with a non-const DistMatrix must set
	// toBeAttachedForGet_ also, if not then it is assumed that
	// the DistMatrix isn't attached
	if ( !toBeAttachedForGet_ )
	    LogicError ("Cannot perform this operation as matrix is not attached.");
    DistMatrix<T>& X = *GlobalArrayGet_;

    const Int height = X.Height ();
    const Int width = X.Width ();
    if (i + height > X.Height () || j + width > X.Width ())
	LogicError ("Invalid AxpyGlobalToLocal submatrix");

    const Grid & g = X.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int p = g.Size ();

    std::vector<T> recvVector_;

    T* XBuffer = Z.Buffer();
    // Send out the requests to all processes in the grid
    // 0 byte send, no need to cancel or post matching 
    // receive
    for (Int rank = 0; rank < p; ++rank)
    {
	const Int bufferSize = 0;
	const Int index =
	    NextIndex (0, putVectors_[rank],
		    requestRequests_[rank],
		    requestRequestStatuses_[rank],
		    matrixBase_[rank],
		    XBuffer);

	// Copy the request header into the send buffer
	T *sendBuffer = putVectors_[rank][index].data ();
	// Begin the non-blocking send
	mpi::TaggedISSend (sendBuffer, bufferSize, rank, REQUEST_GET_TAG, g.VCComm (),
		requestRequests_[rank][index]);
    }
    // Receive all of the replies
    Int numReplies = 0;
    while (numReplies < p)
    {
	mpi::Status status;

	if (mpi::IProbe
		(mpi::ANY_SOURCE, DATA_GET_TAG, g.VCComm (), status))
	{
	    const Int source = status.MPI_SOURCE;

	    // Ensure that we have a recv buffer
	    const Int count = mpi::GetCount < byte > (status);
	    recvVector_.resize (count);
	    T *recvBuffer = recvVector_.data ();

	    // Receive the data
	    mpi::TaggedRecv
		(recvBuffer, count, source, DATA_GET_TAG, g.VCComm ());

	    // Compute the local heights and offsets
	    const Int myRow = g.Row ();
	    const Int myCol = g.Col ();
	    const Int colAlign = (X.ColAlign () + i) % r;
	    const Int rowAlign = (X.RowAlign () + j) % c;
	    const Int colShift = Shift (myRow, colAlign, r);
	    const Int rowShift = Shift (myCol, rowAlign, c);
	    const Int localHeight = Length (height, colShift, r);
	    const Int localWidth = Length (width, rowShift, c);

	    // Unpack the local matrix
	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *YCol = X.Buffer (0, rowShift + t * c);
		const T *XCol = &recvBuffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[colShift + s * r] = XCol[s];
	    }
	    ++numReplies;
	    recvVector_.resize ( 0 );
	}
    }
}

// accumulate = Update Y(i:i+height-1,j:j+width-1) += X,
// where X is height x width
template<typename T>
void AxpyInterface2<T>::Acc( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Acc"))

    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
        LogicError("Global matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayPut_;
    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix");

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;

    const Int XLDim = Z.LDim();
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    const Int YLDim = Y.LDim ();

    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        // number of entries in my PE
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;
	
        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            T* XBuffer = Z.Buffer();
     	    const Int index =
	    NextIndex (numEntries, putVectors_[destination],
		    dataRequests_[destination],
		    dataRequestStatuses_[destination],
		    matrixBase_[destination],
		    XBuffer);

	    DEBUG_ONLY (if
			(Int (putVectors_[destination][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)
	    T *sendBuffer = putVectors_[destination][index].data ();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }
	    // acc request 
	    mpi::TaggedISSend (sendBuffer, numEntries, destination, 
			DATA_ACC_TAG, g.VCComm (), 
			dataRequests_[destination][index]);
        }
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

// get index of a matrix for a particular process 
template<typename T>
Int AxpyInterface2<T>::GetIndexForMatrix ( Matrix<T>& Z, int rank )
{
    typename std::deque<T *>::iterator dit;
    dit = std::find ( matrixBase_[rank].begin(), 
	    matrixBase_[rank].end(), Z.LockedBuffer ());
    const Int index = (dit - matrixBase_[rank].begin());
    //std::cout << "matrixBase size: " << matrixBase_[rank].size () << "\n";
    assert ( index != matrixBase_[rank].size () );
    
    return index;
}

// progress communication for a particular matrix
// this could be used to progress sends and recvs
template<typename T>
void AxpyInterface2<T>::ProgressMatrix ( Matrix<T>& Z, int rank )
{
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    const Int index = GetIndexForMatrix ( Z, rank );

    if ( !dataRequestStatuses_[rank][index] ) // nothing to do
	return;
    // wait
    mpi::Wait ( dataRequests_[rank][index] );
    dataRequestStatuses_[rank][index] = false;
    getVectors_[rank][index].resize (0);
    putVectors_[rank][index].resize (0);
}

// flush ensures local and remote completion
// this interface assumes a send has been issued
// and will post a matching receive and progress
template<typename T>
void AxpyInterface2<T>::Flush( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Flush"))
	if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
	    LogicError("Must initiate transfer before flushing.");
    
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    mpi::Status status;

    if ( mpi::IProbe (mpi::ANY_SOURCE, mpi::ANY_TAG, g.VCComm (), status) )
    {
	switch (status.MPI_TAG)
	{
	    case DATA_PUT_TAG:
		{
		    const Int count = mpi::GetCount <T> (status);
		    HandleLocalToGlobalData ( Z, i, j, count, status.MPI_SOURCE );
		    break;
		}
	    case DATA_ACC_TAG:
		{
		    const Int count = mpi::GetCount <T> (status);
		    HandleLocalToGlobalAcc ( Z, i, j, count, status.MPI_SOURCE );
		    break;
		}
	    case REQUEST_GET_TAG:
		{
		    const Int count = mpi::GetCount <T> (status);
		    HandleGlobalToLocalData ( Z, i, j, count, status.MPI_SOURCE );
		    break;
		}
	}
    }
}

template<typename T>
void AxpyInterface2<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Flush"))
    Flush ( Z, 0, 0 );	
}

template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalData ( Matrix<T>& Z, Int i, Int j, 
	Int count, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    int height = Z.Height();
    int width = Z.Width();
    std::vector<T> getVector_;

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)
	T* Buffer = Z.Buffer();
    getVector_.resize (count);
    DEBUG_ONLY (if
	    (Int (getVector_.size ()) != count) LogicError ("Not enough space allocated");)

	T *getBuffer = getVector_.data ();
    mpi::TaggedRecv (getBuffer, count, source, DATA_PUT_TAG, g.VCComm ());

    // Update Y
    const T *XBuffer = reinterpret_cast < const T * >(getBuffer);
    const Int colAlign = (Y.ColAlign () + i) % r;
    const Int rowAlign = (Y.RowAlign () + j) % c;
    const Int colShift = Shift (myRow, colAlign, r);
    const Int rowShift = Shift (myCol, rowAlign, c);

    const Int localHeight = Length (height, colShift, r);
    const Int localWidth = Length (width, rowShift, c);
    const Int iLocalOffset = Length (i, Y.ColShift (), r);
    const Int jLocalOffset = Length (j, Y.RowShift (), c);

    for (Int t = 0; t < localWidth; ++t)
    {
	T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
	const T *XCol = &XBuffer[t * localHeight];
	for (Int s = 0; s < localHeight; ++s)
	    YCol[s] = XCol[s];
    }
    getVector_.resize ( 0 );
    // Free the memory
    //ProgressMatrix ( Z, source );
}

// replica of above function except this accumulates
template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalAcc ( Matrix<T>& Z, Int i, Int j, 
	Int count, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    int height = Z.Height();
    int width = Z.Width();
    std::vector<T> getVector_;

    	DEBUG_ONLY (if (count < Int (sizeof (T)))
		LogicError ("Count was too small");)
	T* Buffer = Z.Buffer();
        getVector_.resize (count);
	DEBUG_ONLY (if
		(Int (getVector_.size ()) != count) LogicError ("Not enough space allocated");)
	
	T *getBuffer = getVector_.data ();
	mpi::TaggedRecv (getBuffer, count, source, DATA_ACC_TAG, g.VCComm ());

	// Update Y
	const T *XBuffer = reinterpret_cast < const T * >(getBuffer);
	const Int colAlign = (Y.ColAlign () + i) % r;
	const Int rowAlign = (Y.RowAlign () + j) % c;
	const Int colShift = Shift (myRow, colAlign, r);
	const Int rowShift = Shift (myCol, rowAlign, c);

	const Int localHeight = Length (height, colShift, r);
	const Int localWidth = Length (width, rowShift, c);
	const Int iLocalOffset = Length (i, Y.ColShift (), r);
	const Int jLocalOffset = Length (j, Y.RowShift (), c);

	for (Int t = 0; t < localWidth; ++t)
	{
	    T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
	    const T *XCol = &XBuffer[t * localHeight];
	    for (Int s = 0; s < localHeight; ++s)
		YCol[s] += XCol[s];
	}
        getVector_.resize ( 0 );
	// Free the memory
	//ProgressMatrix ( Z, source );
}

// handle request for data, post a matching issend
template < typename T >
void AxpyInterface2<T>::HandleGlobalToLocalData ( Matrix<T>& Z, Int i, Int j, 
	Int count, Int source )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface::HandleGlobalToLocalData"))

	if( i < 0 || j < 0 )
	    LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForGet_ )
	LogicError("Local matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;
    std::vector<T> putVector_;

    const Int XLDim = Z.LDim();
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    const Int colShift = Shift (myRow, colAlign, r);
    const Int rowShift = Shift (myCol, rowAlign, c);
    const Int localHeight = Length (height, colShift, r);
    const Int localWidth = Length (width, rowShift, c);

    const Int iLocalOffset = Length (i, Y.ColShift (), r);
    const Int jLocalOffset = Length (j, Y.RowShift (), c);

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)

    putVector_.resize (count);
    DEBUG_ONLY (if
	    (Int (putVector_.size ()) != count) LogicError ("Not enough space allocated");)

    T *sendBuffer = putVector_.data ();

    for (Int t = 0; t < localWidth; ++t)
    {
	T *sendCol = &sendBuffer[t * localHeight];
	const T *XCol = Y.LockedBuffer (iLocalOffset, jLocalOffset + t);
	MemCopy (sendCol, XCol, localHeight);
    }

    // Fire off non-blocking send
    mpi::TaggedSend (sendBuffer, count, source, 
	    DATA_GET_TAG, g.VCComm ());
    // clear
    //ProgressMatrix ( Z, source );

    putVector_.resize ( 0 );
}

// detach collectively 
template<typename T>
void AxpyInterface2<T>::Detach()
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Detach"))
	// destructor will call detach again...
	if (detached_)
	    return;
    if( !attached_ )
	LogicError("Must attach before detaching.");

    const Grid& g = ( toBeAttachedForPut_ ?
	    GlobalArrayPut_->Grid() :
	    GlobalArrayGet_->Grid() );

    mpi::Barrier( g.VCComm() );

    attached_ 		= false;
    detached_ 		= true;
    toBeAttachedForPut_ = false;
    toBeAttachedForGet_ = false;

    GlobalArrayPut_ 	= 0;
    GlobalArrayGet_ 	= 0;

    putVectors_.clear();
    getVectors_.clear();

    dataRequests_.clear();
    dataRequestStatuses_.clear();
    requestRequests_.clear();
    requestRequestStatuses_.clear();    
    
    matrixBase_.clear();
}

template class AxpyInterface2<Int>;
template class AxpyInterface2<float>;
template class AxpyInterface2<double>;
template class AxpyInterface2<Complex<float>>;
template class AxpyInterface2<Complex<double>>;

} // namespace El
