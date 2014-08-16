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
namespace El
{
template<typename T>
AxpyInterface2<T>::AxpyInterface2()
    : GlobalArrayPut_(0), GlobalArrayGet_(0),
    sendVectors_(0), recvVectors_(0), replyVectors_(0), coordVectors_(0), 
    requestVectors_(0), sendRequests_(0), recvRequests_(0), 
    replyRequests_(0), requestRequests_(0), sendIJRequests_(0),
    matrixBase_(0), sendRequestStatuses_(0), requestRequestStatuses_(0), 
    replyRequestStatuses_(0), recvRequestStatuses_(0), 
    sendIJRequestStatuses_(0), toBeAttachedForPut_(false), 
    toBeAttachedForGet_(false), attached_(false), detached_(false)
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
    
    requestVectors_.resize( p );
    sendVectors_.resize( p );
    recvVectors_.resize( p );
    replyVectors_.resize( p );
    coordVectors_.resize( p );
    
    sendRequests_.resize (p);
    recvRequests_.resize (p);
    replyRequests_.resize (p);
    requestRequests_.resize (p);
    sendIJRequests_.resize (p);
    
    sendRequestStatuses_.resize (p);    
    recvRequestStatuses_.resize (p);    
    requestRequestStatuses_.resize (p);    
    replyRequestStatuses_.resize (p);    
    sendIJRequestStatuses_.resize (p);    
    
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

    if (sendVectors_.size() != p)
    {
	recvVectors_.resize( p );
	sendVectors_.resize( p );
	replyVectors_.resize( p );
	requestVectors_.resize( p );
	coordVectors_.resize( p );
	
	sendRequests_.resize (p);
	recvRequests_.resize (p);
	requestRequests_.resize (p);
	replyRequests_.resize (p);
	sendIJRequests_.resize (p);    
	
	sendRequestStatuses_.resize (p);
	recvRequestStatuses_.resize (p);
	requestRequestStatuses_.resize (p);    
	replyRequestStatuses_.resize (p);    
	sendIJRequestStatuses_.resize (p);    
	
	matrixBase_.resize (p);
    }
}

template<typename T>
Int AxpyInterface2<T>::NextIndex
( Int rank, Int dataSize, 
  std::deque <std::vector<T>> &dataVectors,
  std::deque <mpi::Request> &requests,
  std::deque <bool> &requestStatus,
  T * base_address )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndex"))
    const Int Index = Int(requests.size ());
    
    dataVectors.resize (Index + 1);
    dataVectors[Index].resize (dataSize);

    requests.push_back (mpi::REQUEST_NULL);
    requestStatus.push_back ( true );
    
    // stores Matrix base address by index
    matrixBase_[rank].push_back (base_address);
    
    return Index;
}

template<typename T>
Int AxpyInterface2<T>::NextIndex
( Int rank, Int i, Int j, Int dataSize, 
  std::deque <std::vector<T>> &dataVectors,
  std::deque <mpi::Request> &requestData,
  std::deque <bool> &requestDataStatus,
  T * base_address) 
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndex"))
    const Int Index = Int(requestData.size ());
   
    dataVectors.resize (Index + 1);
    dataVectors[Index].resize (dataSize);

    coordVectors_[rank].resize (Index + 1);
    coordVectors_[rank][Index].resize (2);
    coordVectors_[rank][Index][0] = i;
    coordVectors_[rank][Index][1] = j;
    
    requestData.push_back (mpi::REQUEST_NULL);
    requestDataStatus.push_back ( true );
    
    sendIJRequests_[rank].push_back (mpi::REQUEST_NULL);
    sendIJRequestStatuses_[rank].push_back ( true );
    
    // stores Matrix base address by index
    matrixBase_[rank].push_back (base_address);
    
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
	    NextIndex (destination, i, j, numEntries, 
		    sendVectors_[destination],
		    sendRequests_[destination],
		    sendRequestStatuses_[destination],
		    XBuffer);

	    DEBUG_ONLY (if
			(Int (sendVectors_[destination][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)

	    T *sendBuffer = sendVectors_[destination][index].data ();
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
			sendRequests_[destination][index]);
	
	    // send coordinates
	    Int *coord = coordVectors_[destination][index].data ();
	    mpi::TaggedISend (coord, 2, destination, COORD_IJ_TAG, g.VCComm (), 
		    sendIJRequests_[destination][index]);
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

    const Int height = Z.Height ();
    const Int width = Z.Width ();
    
    if (i + height > X.Height () || j + width > X.Width ())
	LogicError ("Invalid AxpyGlobalToLocal submatrix");

    const Grid & g = X.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int p = g.Size ();

    std::vector<T> recvVector_;

    T* XBuffer = Z.Buffer();
    // Send out the requests to all processes in the grid
    for (Int rank = 0; rank < p; ++rank)
    {
	// we just use the request objects for progress
	const Int index =
	    NextIndex (rank, i, j, 1, 
		    requestVectors_[rank],
		    requestRequests_[rank],
		    requestRequestStatuses_[rank],
		    XBuffer);
	// send request
	T *requestBuffer = requestVectors_[rank][index].data();
	mpi::TaggedISSend (requestBuffer, 1, rank, 
		REQUEST_GET_TAG, g.VCComm(),
		requestRequests_[rank][index]);
    		
	// send coordinates
	Int *coord = coordVectors_[rank][index].data ();
	mpi::TaggedISend (coord, 2, rank, 
		COORD_IJ_TAG, g.VCComm (), 
		sendIJRequests_[rank][index]);
    }

    // Receive all of the replies
    Int numReplies = 0;
    while (numReplies < p)
    {
	mpi::Status status;
	HandleGlobalToLocalData ( Z );
	if (mpi::IProbe
		(mpi::ANY_SOURCE, DATA_GET_TAG, g.VCComm (), status))
	{
	    const Int source = status.MPI_SOURCE;
	    // Ensure that we have a recv buffer
	    const Int count = mpi::GetCount <T> (status);
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
	    recvVector_.clear();
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

    const Int YLDim = Y.LDim();

    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        // number of entries in my PE
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;
	
        if( numEntries != 0  )
        {
            const Int destination = receivingRow + r*receivingCol;
            T* XBuffer = Z.Buffer();
     	    const Int index =
	    NextIndex (destination, i, j, numEntries, 
		    sendVectors_[destination],
		    sendRequests_[destination],
		    sendRequestStatuses_[destination],
		    XBuffer);

	    DEBUG_ONLY (if
			(Int (sendVectors_[destination][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)
	    	
	    T *sendBuffer = sendVectors_[destination][index].data ();
	    for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

	    // send data 
	    mpi::TaggedISSend (sendBuffer, numEntries, destination, 
	    	DATA_ACC_TAG, g.VCComm(), 
	   	sendRequests_[destination][index]);

	    // send coordinates
	    Int *coord = coordVectors_[destination][index].data();
	    mpi::TaggedISend (coord, 2, destination, 
		    COORD_IJ_TAG, g.VCComm(), 
		    sendIJRequests_[destination][index]);
	}
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

// progress communication for a particular matrix
// progress requests
template<typename T>
bool AxpyInterface2<T>::TestRequests ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    
    Int index;
    typename std::deque<T *>::iterator dit;

    for (int rank = 0; rank < p; ++rank)
    {
	dit = std::find ( matrixBase_[rank].begin(), 
		matrixBase_[rank].end(), Z.LockedBuffer());
	index = (dit - matrixBase_[rank].begin());
	
	if ( index == matrixBase_[rank].size() )
	    continue;
    	if ( requestRequestStatuses_[rank].size() == 0 )
	    continue;
	// test all send requests related to matrix
	const Int numStatuses = requestRequestStatuses_[rank].size();
 	for (int i = 0; i < numStatuses; i++)
	{
	    requestRequestStatuses_[rank][i] = !mpi::Test ( requestRequests_[rank][i] );
	    if ( requestRequestStatuses_[rank][i] )
		return false;
	}
	// okay to deallocate
	requestVectors_[rank][index].clear();
    }
    return true;
}

// progress sends
template<typename T>
bool AxpyInterface2<T>::TestSends ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    
    Int index;
    typename std::deque<T *>::iterator dit;

    for (int rank = 0; rank < p; ++rank)
    {
	dit = std::find ( matrixBase_[rank].begin(), 
		matrixBase_[rank].end(), Z.LockedBuffer());
	index = (dit - matrixBase_[rank].begin());
	
	if ( index == matrixBase_[rank].size() )
	    continue;
    	if ( sendRequestStatuses_[rank].size() == 0 )
	    continue;
	
	// test all sends related to matrix
	const Int numStatuses = sendRequestStatuses_[rank].size();
 	for (int i = 0; i < numStatuses; i++)
	{
	    sendRequestStatuses_[rank][i] = !mpi::Test ( sendRequests_[rank][i] );
	    if ( sendRequestStatuses_[rank][i] )
		return false;
	}

	// if test is true, then it is safe to free buffer
	sendVectors_[rank][index].clear();
    }
    return true;
}

// progress recvs
template<typename T>
bool AxpyInterface2<T>::TestRecvs ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    
    Int index;
    typename std::deque<T *>::iterator dit;

    for (int rank = 0; rank < p; ++rank)
    {
	dit = std::find ( matrixBase_[rank].begin(), 
		matrixBase_[rank].end(), Z.LockedBuffer());
	index = (dit - matrixBase_[rank].begin());
	
	if ( index == matrixBase_[rank].size() )
	    continue;
    	if ( recvRequestStatuses_[rank].size() == 0 )
	    continue;

	// test all sends related to matrix
	const Int numStatuses = recvRequestStatuses_[rank].size();
 	for (int i = 0; i < numStatuses; i++)
	{
	    recvRequestStatuses_[rank][i] = !mpi::Test ( recvRequests_[rank][i] );
	    if ( recvRequestStatuses_[rank][i] )
		return false;
	}
	
	// if test is true, then it is safe to free buffer
	recvVectors_[rank][index].clear();
    }
    return true;
}

// progress replies
template<typename T>
bool AxpyInterface2<T>::TestReplies ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    
    Int index;
    typename std::deque<T *>::iterator dit;

    for (int rank = 0; rank < p; ++rank)
    {
	dit = std::find ( matrixBase_[rank].begin(), 
		matrixBase_[rank].end(), Z.LockedBuffer());
	index = (dit - matrixBase_[rank].begin());
	
	if ( index == matrixBase_[rank].size() )
	    continue;
    	if ( replyRequestStatuses_[rank].size() == 0 )
	    continue;

	// test all sends related to matrix
	const Int numStatuses = replyRequestStatuses_[rank].size();
 	for (int i = 0; i < numStatuses; i++)
	{
	    replyRequestStatuses_[rank][i] = !mpi::Test ( replyRequests_[rank][i] );
	    if ( replyRequestStatuses_[rank][i] )
		return false;
	}

	// if test is true, then it is safe to free buffer
	replyVectors_[rank][index].clear();
    }
    return true;
}

// progress coordinate sends
template<typename T>
bool AxpyInterface2<T>::TestSendsCoord ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    
    Int index;
    typename std::deque<T *>::iterator dit;

    for (int rank = 0; rank < p; ++rank)
    {
	dit = std::find ( matrixBase_[rank].begin(), 
		matrixBase_[rank].end(), Z.LockedBuffer());
	index = (dit - matrixBase_[rank].begin());
	
	if ( index == matrixBase_[rank].size() )
	    continue;
    	if ( sendIJRequestStatuses_[rank].size() == 0 )
	    continue;

	// test all sends related to matrix
	const Int numStatuses = sendIJRequestStatuses_[rank].size();
 	for (int i = 0; i < numStatuses; i++)
	{
	    sendIJRequestStatuses_[rank][i] = !mpi::Test ( sendIJRequests_[rank][i] );
	    if ( sendIJRequestStatuses_[rank][i] )
		return false;
	}
	// if test is true, then it is safe to free buffer
	coordVectors_[rank][index].clear();
    }
    return true;
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

    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    
    mpi::Status status;
	
    bool DONE = false;
    mpi::Request nb_bar_request;
    bool nb_bar_active = false;

    while ( !DONE )
    {
	// similar to HandleXYZ functions in original AxpyInterface
	if ( mpi::IProbe (mpi::ANY_SOURCE, mpi::ANY_TAG, g.VCComm(), status) )
	{
	    switch (status.MPI_TAG)
	    {
		case DATA_PUT_TAG:
		    {
			const Int count = mpi::GetCount <T> (status);
			HandleLocalToGlobalData ( Z, count, status.MPI_SOURCE );
			break;
		    }
		case DATA_ACC_TAG:
		    {
			const Int count = mpi::GetCount <T> (status);
			HandleLocalToGlobalAcc ( Z, count, status.MPI_SOURCE );
			break;
		    }
		case REQUEST_GET_TAG:
		    {
			HandleGlobalToLocalData ( Z );
			break;
		    }
	    }
	}
	if ( nb_bar_active )
	{
	    DONE = mpi::Test ( nb_bar_request );
	}
	else
	{
	    // check if all sends (data or request) are 
	    // complete for a particular matrix
	    if ( TestSends( Z ) 
	    	    && TestRecvs( Z ) 
	    	    && TestRequests( Z ) 
	    	    && TestReplies ( Z ) 
	    	    && TestSendsCoord ( Z ) )
	    {
		mpi::IBarrier ( g.VCComm(), nb_bar_request );
		nb_bar_active = true;
	    }
	}
    }
}

// all communications pertaining to matrix Z
template<typename T>
void AxpyInterface2<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Flush"))
    Flush ( Z, 0, 0 );	
}

template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalData ( Matrix<T>& Z, Int count, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    int height = Z.Height();
    int width = Z.Width();
    // data vector
    std::vector<T> getVector_;
    getVector_.resize (count);

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)
    DEBUG_ONLY (if (Int (getVector_.size ()) != count) 
	    LogicError ("Not enough space allocated");)

    // post receive for coordinates
    Int coord[2];
    mpi::TaggedRecv (coord, 2, source, 
	    COORD_IJ_TAG, g.VCComm());
    Int i = coord[0]; 
    Int j = coord[1];
 
    // post receive for data
    T *getBuffer = getVector_.data();
    mpi::TaggedRecv (getBuffer, count, source, 
	    DATA_PUT_TAG, g.VCComm());
        
    // Update Y
    const T *XBuffer = reinterpret_cast < const T * >(getBuffer);
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;
    const Int colShift = Shift (myRow, colAlign, r);
    const Int rowShift = Shift (myCol, rowAlign, c);

    const Int localHeight = Length (height, colShift, r);
    const Int localWidth = Length (width, rowShift, c);
    const Int iLocalOffset = Length (i, Y.ColShift(), r);
    const Int jLocalOffset = Length (j, Y.RowShift(), c);

    for (Int t = 0; t < localWidth; ++t)
    {
	T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
	const T *XCol = &XBuffer[t * localHeight];
	for (Int s = 0; s < localHeight; ++s)
	    YCol[s] = XCol[s];
    }
    // Free the memory
    getVector_.clear();
}

// replica of above function except this accumulates
template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalAcc ( Matrix<T>& Z, Int count, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    const int height = Z.Height();
    const int width = Z.Width();

    // data buffer
    std::vector<T> getVector_;
    getVector_.resize (count);

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)

    DEBUG_ONLY (if (Int (getVector_.size ()) != count) 
	    LogicError ("Not enough space allocated");)
    
    // post receive for coordinates
    Int coord[2];
    mpi::TaggedRecv (coord, 2, source, 
	    COORD_IJ_TAG, g.VCComm());
    Int i = coord[0]; Int j = coord[1];
    
    // post receive for data
    T *getBuffer = getVector_.data();
    mpi::TaggedRecv (getBuffer, count, source, 
	    DATA_ACC_TAG, g.VCComm());

    // Update Y
    const T *XBuffer = reinterpret_cast < const T * >(getBuffer);
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;
    const Int colShift = Shift (myRow, colAlign, r);
    const Int rowShift = Shift (myCol, rowAlign, c);

    const Int localHeight = Length (height, colShift, r);
    const Int localWidth = Length (width, rowShift, c);
    const Int iLocalOffset = Length (i, Y.ColShift(), r);
    const Int jLocalOffset = Length (j, Y.RowShift(), c);

    for (Int t = 0; t < localWidth; ++t)
    {
	T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
	const T *XCol = &XBuffer[t * localHeight];
	for (Int s = 0; s < localHeight; ++s)
	    YCol[s] += XCol[s];
    }
    // Free the memory
    getVector_.clear();
}

// handle request for data, post a matching issend
template < typename T >
void AxpyInterface2<T>::HandleGlobalToLocalData ( Matrix<T>& Z )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface::HandleGlobalToLocalData"))

    if ( !toBeAttachedForGet_ )
	LogicError("Local matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myRow = g.Row();
    const Int myCol = g.Col();

    mpi::Status status;

    if (mpi::IProbe (mpi::ANY_SOURCE, REQUEST_GET_TAG, g.VCComm (), status))
    {
	const Int source = status.MPI_SOURCE;
	// dummy var for receiving request
	// we don't use this anyway
	T dummy_=0;
	// post receive request for get
	mpi::TaggedRecv (&dummy_, 1, source, 
		REQUEST_GET_TAG, g.VCComm());
    
	// post receive for coordinates
	Int coord[2];
	mpi::TaggedRecv (coord, 2, source, 
		COORD_IJ_TAG, g.VCComm());
    	Int i = coord[0]; Int j = coord[1];
 
	const Int colAlign = (Y.ColAlign() + i) % r;
	const Int rowAlign = (Y.RowAlign() + j) % c;

	const Int XLDim = Z.LDim();
	// local matrix width and height
	const Int height = Z.Height();
	const Int width = Z.Width();

	const Int colShift = Shift (myRow, colAlign, r);
	const Int rowShift = Shift (myCol, rowAlign, c);
	const Int localHeight = Length (height, colShift, r);
	const Int localWidth = Length (width, rowShift, c);

	const Int iLocalOffset = Length (i, Y.ColShift (), r);
	const Int jLocalOffset = Length (j, Y.RowShift (), c);

	const Int numEntries = localHeight * localWidth;

	DEBUG_ONLY (if (numEntries < Int (sizeof (T)))
		LogicError ("Count was too small");)

	T* XBuffer = Z.Buffer();
	const Int index =
	    NextIndex (source, numEntries, replyVectors_[source],
		    replyRequests_[source],
		    replyRequestStatuses_[source],
		    XBuffer);

	DEBUG_ONLY (if
		(Int (replyVectors_[source][index].size ()) !=
		 numEntries) LogicError ("Error in NextIndex");)
	
	T *replyBuffer = replyVectors_[source][index].data ();
	for (Int t = 0; t < localWidth; ++t)
	{
	    T *sendCol = &replyBuffer[t * localHeight];
	    const T *XCol = Y.LockedBuffer (iLocalOffset, jLocalOffset + t);
	    MemCopy (sendCol, XCol, localHeight);
	}

	// Fire off non-blocking send
	mpi::TaggedISSend (replyBuffer, numEntries, source, 
		DATA_GET_TAG, g.VCComm (), 
		replyRequests_[source][index]);
    }
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

    sendVectors_.clear();
    coordVectors_.clear();
    recvVectors_.clear();
    replyVectors_.clear();
    requestVectors_.clear();

    sendIJRequests_.clear();
    sendRequests_.clear();
    recvRequests_.clear();
    requestRequests_.clear();
    replyRequests_.clear();
    
    sendRequestStatuses_.clear();
    sendIJRequestStatuses_.clear();
    recvRequestStatuses_.clear();
    requestRequestStatuses_.clear();    
    replyRequestStatuses_.clear();    
    
    matrixBase_.clear();
}

template class AxpyInterface2<Int>;
template class AxpyInterface2<float>;
template class AxpyInterface2<double>;
template class AxpyInterface2<Complex<float>>;
template class AxpyInterface2<Complex<double>>;

} // namespace El
