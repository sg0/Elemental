/*
This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include <cassert>

#if MPI_VERSION>=3
// TODO Use DDT for put/get/acc when EL_USE_DERIVED_TYPE is defined
// TODO bring back const interfaces
namespace El
{
    template<typename T>
	AxpyInterface2<T>::AxpyInterface2()
	: GlobalArrayPut_(0), GlobalArrayGet_(0),
	sendDataStatuses_(0), sendCoordStatuses_(0),
	recvDataStatuses_(0), recvCoordStatuses_(0),
	sendDataRequests_(0), sendCoordRequests_(0),
	recvDataRequests_(0), recvCoordRequests_(0),
	sendData_(0), recvData_(0),
	sendCoord_(0), recvCoord_(0),
	put_win_(0), acc_win_(0), getrq_win_(0),
	put_win_base_(0), acc_win_base_(0), getrq_win_base_(0),
	toBeAttachedForGet_(false), toBeAttachedForPut_(false),
	attached_(false), detached_(false)
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

    const Grid& g = Z.Grid();
    const Int p = g.Size ();
    
    if ( sendData_.empty() )
    {
	sendDataStatuses_.resize (p);
	sendCoordStatuses_.resize (p);
	recvDataStatuses_.resize (p);
	recvCoordStatuses_.resize (p);
	
	sendDataRequests_.resize (p);
	sendCoordRequests_.resize (p);
	recvDataRequests_.resize (p);
	recvCoordRequests_.resize (p);
    	
	sendData_.resize (p);
	sendCoord_.resize (p);
	recvData_.resize (p);
	recvCoord_.resize (p);
    }

    // count window related
    put_win_base_ = new long;
    mpi::WindowCreate ( put_win_base_, sizeof(long), 
	    g.VCComm(), put_win_ );
    memset (put_win_base_, 0, sizeof (long));
    mpi::WindowLock (put_win_);

    acc_win_base_ = new long;
    mpi::WindowCreate ( acc_win_base_, sizeof(long), 
	    g.VCComm(), acc_win_ );
    memset (acc_win_base_, 0, sizeof (long));
    mpi::WindowLock (acc_win_);

    getrq_win_base_ = new long;
    mpi::WindowCreate ( getrq_win_base_, sizeof(long), 
	    g.VCComm(), getrq_win_ );
    memset (getrq_win_base_, 0, sizeof (long));
    mpi::WindowLock (getrq_win_);
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

    // the matrix base_ is not known until 
    // an update operation (put/get/acc)
    // so it is kept blank

    // if DistMatrix is non-const, all one-sided
    // transfers -- put, get and acc are possible
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
    {
	GlobalArrayPut_ 	= &Z;
	toBeAttachedForPut_ 	= true;
	GlobalArrayGet_ 	= &Z;
	toBeAttachedForGet_ 	= true;

	const Grid& g = Z.Grid();
	const Int p = g.Size ();

	if ( sendData_.empty() )
	{
	    sendDataStatuses_.resize (p);
	    sendCoordStatuses_.resize (p);
	    recvDataStatuses_.resize (p);
	    recvCoordStatuses_.resize (p);

	    sendDataRequests_.resize (p);
	    sendCoordRequests_.resize (p);
	    recvDataRequests_.resize (p);
	    recvCoordRequests_.resize (p);

	    sendData_.resize (p);
	    sendCoord_.resize (p);
	    recvData_.resize (p);
	    recvCoord_.resize (p);
	}
	// count window related
	put_win_base_ = new long;
	mpi::WindowCreate ( put_win_base_, sizeof(long), 
		g.VCComm(), put_win_ );
	memset (put_win_base_, 0, sizeof (long));
	mpi::WindowLock (put_win_);

	acc_win_base_ = new long;
	mpi::WindowCreate ( acc_win_base_, sizeof(long), 
		g.VCComm(), acc_win_ );
	memset (acc_win_base_, 0, sizeof (long));
	mpi::WindowLock (acc_win_);

	getrq_win_base_ = new long;
	mpi::WindowCreate ( getrq_win_base_, sizeof(long), 
		g.VCComm(), getrq_win_ );
	memset (getrq_win_base_, 0, sizeof (long));
	mpi::WindowLock (getrq_win_);
    }
}

template<typename T>
    Int AxpyInterface2<T>::NextIndexData
    (Int dataSize,
     std::deque < std::vector < T >> &data,
     std::deque < mpi::Request > &requests,
     std::deque < bool > &requestStatuses)
  {
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndexData"))
	const Int numCreated = data.size ();
    DEBUG_ONLY (if (numCreated != Int (requests.size ()) ||
		    numCreated !=
		    Int (requestStatuses.size ()))LogicError
		("size mismatch");)
	for (Int i = 0; i < numCreated; ++i)
	{
	    // If this request is still running, test to see if it finished.
	  if (requestStatuses[i])
	    {
	      const bool finished = mpi::Test (requests[i]);
	      requestStatuses[i] = !finished;
	    }

	  if (!requestStatuses[i])
	    {
	      requestStatuses[i] = true;
	      data[i].resize (dataSize);
	      return i;
	    }
	}

    data.resize (numCreated + 1);
    data[numCreated].resize (dataSize);
    requests.push_back (mpi::REQUEST_NULL);
    requestStatuses.push_back (true);
    
    return numCreated;
  }

template<typename T>
    Int AxpyInterface2<T>::NextIndexCoord
    (std::deque < std::array<Int, 3> > &coord,
    std::deque < mpi::Request > &requests,
     std::deque < bool > &requestStatuses)
  {
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndexCoord"))
    
    const Int numCreated = coord.size ();
    DEBUG_ONLY (if (numCreated != Int (requests.size ()) ||
		    numCreated !=
		    Int (requestStatuses.size ()))LogicError
		("size mismatch");)
	
	for (Int i = 0; i < numCreated; ++i)
	{
	    // If this request is still running, test to see if it finished.
	  if (requestStatuses[i])
	    {
	      const bool finished = mpi::Test (requests[i]);
	      requestStatuses[i] = !finished;
	    }

	  if (!requestStatuses[i])
	    {
	      requestStatuses[i] = true;
	      return i;
	    }
	}

    coord.resize (numCreated + 1);
    requests.push_back (mpi::REQUEST_NULL);
    requestStatuses.push_back (true);
    
    return numCreated;
  }

template<typename T>
void AxpyInterface2<T>::Iput( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Iput"))

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
            const T* XBuffer = Z.LockedBuffer();

	    const Int dindex =
	    NextIndexData (numEntries,
		    sendData_[destination],
		    sendDataRequests_[destination],
		    sendDataStatuses_[destination]);

	    DEBUG_ONLY (if
			(Int (sendData_[destination][dindex].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)

	    T *sendBuffer = sendData_[destination][dindex].data ();
            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }
	    // put request 
	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
			DATA_PUT_TAG, g.VCComm (), 
			sendDataRequests_[destination][dindex]);
		    
	    // send coordinates
	    const Int cindex =
	    NextIndexCoord (sendCoord_[destination],
		    sendCoordRequests_[destination],
		    sendCoordStatuses_[destination]);

	    Int *coord = reinterpret_cast<Int *>(sendCoord_[destination][cindex].data ());
	    coord[0] = i; 
	    coord[1] = j;
	    coord[2] = numEntries;

	    mpi::TaggedISend (coord, 3, destination, COORD_PUT_TAG, g.VCComm (), 
		    sendCoordRequests_[destination][cindex]);
	
	    // put count
	    mpi::ReadInc (put_win_, 0, 1, destination);
	}

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Iget( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Iget"))
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

    const T* XBuffer = Z.LockedBuffer();
    // Send out the requests to all processes in the grid
    for (Int rank = 0; rank < p; ++rank)
    {
	// send coordinates
	const Int cindex =
	    NextIndexCoord (sendCoord_[rank],
		    sendCoordRequests_[rank],
		    sendCoordStatuses_[rank]);

	Int *coord = reinterpret_cast<Int *>(sendCoord_[rank][cindex].data ());
        coord[0] = i; 
	coord[1] = j;
        coord[2] = -1;

	mpi::TaggedISend (coord, 3, rank, 
		REQUEST_GET_TAG, g.VCComm (), 
		sendCoordRequests_[rank][cindex]);
	// get request count
	mpi::ReadInc (getrq_win_, 0, 1, rank);
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
void AxpyInterface2<T>::Iacc( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Iacc"))

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
            const T* XBuffer = Z.LockedBuffer();

	    // send data 
     	    const Int dindex =
	    NextIndexData (numEntries,
		    sendData_[destination],
		    sendDataRequests_[destination],
		    sendDataStatuses_[destination]);

	    DEBUG_ONLY (if
			(Int (sendData_[destination][dindex].size ()) !=
			 numEntries) LogicError ("Error in NextIndexData");)
	    
	    T *sendBuffer = sendData_[destination][dindex].data ();
	    for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
	    	DATA_ACC_TAG, g.VCComm(), 
	   	sendDataRequests_[destination][dindex]);
		
	   // send coordinates
	    const Int cindex =
		NextIndexCoord (sendCoord_[destination],
		    sendCoordRequests_[destination],
		    sendCoordStatuses_[destination]);

	    Int *coord = reinterpret_cast<Int *>(sendCoord_[destination][cindex].data());
	    coord[0] = i; 
	    coord[1] = j;
	    coord[2] = numEntries;

	    mpi::TaggedISend (coord, 3, destination, 
		    COORD_ACC_TAG, g.VCComm(), 
		    sendCoordRequests_[destination][cindex]);   

	    // acc count
	    mpi::ReadInc (acc_win_, 0, 1, destination);
	}

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
bool AxpyInterface2<T>::TestRequests( Matrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::TestRequests"))
    
    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int p = g.Size ();

    for (Int i = 0; i < p; ++i)
    {
	// coord recvs
	const Int numrecvCoordRequests = recvCoordRequests_[i].size ();
	for (Int j = 0; j < numrecvCoordRequests; ++j)
	{
	    recvCoordStatuses_[i][j] = 
		!mpi::Test (recvCoordRequests_[i][j]);
	    if (recvCoordStatuses_[i][j])
		return false;
	} 
	
	// coord sends
	const Int numsendCoordRequests = sendCoordRequests_[i].size ();
	for (Int j = 0; j < numsendCoordRequests; ++j)
	{
	    sendCoordStatuses_[i][j] = 
		!mpi::Test (sendCoordRequests_[i][j]);
	    if (sendCoordStatuses_[i][j])
		return false;
	}

	// data recvs
	const Int numrecvDataRequests = recvDataRequests_[i].size ();
	for (Int j = 0; j < numrecvDataRequests; ++j)
	{
	    recvDataStatuses_[i][j] = 
		!mpi::Test (recvDataRequests_[i][j]);
	    if (recvDataStatuses_[i][j])
		return false;
	}

	// data sends
	const Int numsendDataRequests = sendDataRequests_[i].size ();
	for (Int j = 0; j < numsendDataRequests; ++j)
	{
	    sendDataStatuses_[i][j] = 
		!mpi::Test (sendDataRequests_[i][j]);
	    if (sendDataStatuses_[i][j])
		return false;
	}
    }

    return true;
}

// flush ensures local and remote completion
// this interface assumes a send has been issued
// and will post a matching receive and progress
    template<typename T>
void AxpyInterface2<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Flush"))
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
	LogicError("Must initiate transfer before flushing.");

    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int p = g.Size ();
    const Int me = g.VCRank();
    
    // get my put/get/acc recv counts
    const Int put_count = mpi::ReadInc (put_win_, 0, 0, me);
    const Int acc_count = mpi::ReadInc (acc_win_, 0, 0, me);
    const Int getrq_count = mpi::ReadInc (getrq_win_, 0, 0, me);
	    
    TestRequests (Z);

    for (Int count = 0; count < put_count; ++count)
    {
	mpi::Status status;
	if ( mpi::IProbe (mpi::ANY_SOURCE, DATA_PUT_TAG, g.VCComm(), status) )
	    HandleLocalToGlobalData ( Z, status.MPI_SOURCE );
    }

    for (Int count = 0; count < acc_count; ++count)
    {
	mpi::Status status;
	if ( mpi::IProbe (mpi::ANY_SOURCE, DATA_ACC_TAG, g.VCComm(), status) )
	    HandleLocalToGlobalAcc ( Z, status.MPI_SOURCE );
    }

    for (Int count = 0; count < getrq_count; ++count)
    {
	mpi::Status status;
	if ( mpi::IProbe (mpi::ANY_SOURCE, REQUEST_GET_TAG, g.VCComm(), status) )
	    HandleGlobalToLocalData ( Z ); 
    }

}

template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalData ( Matrix<T>& Z, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    int height = Z.Height();
    int width = Z.Width();

    // post receive for coordinates
    Int coord[3];
    mpi::TaggedRecv (coord, 3, source, 
	    COORD_PUT_TAG, g.VCComm());
    Int i = coord[0]; 
    Int j = coord[1];
    Int count = coord[2];
 
    // data vector
    std::vector<T> getVector_;
    getVector_.resize (count);

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)
    DEBUG_ONLY (if (Int (getVector_.size ()) != count) 
	    LogicError ("Not enough space allocated");)

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
void AxpyInterface2<T>::HandleLocalToGlobalAcc ( Matrix<T>& Z, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    const int height = Z.Height();
    const int width = Z.Width();

    // post receive for coordinates
    Int coord[3];
    mpi::TaggedRecv (coord, 3, source, 
	    COORD_ACC_TAG, g.VCComm());
    Int i = coord[0]; 
    Int j = coord[1];
    Int count = coord[2];
 
    // data buffer
    std::vector<T> getVector_;
    getVector_.resize (count);

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)

    DEBUG_ONLY (if (Int (getVector_.size ()) != count) 
	    LogicError ("Not enough space allocated");)   

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

// handle request for data, post a matching isend
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

	// post receive for coordinates
	Int coord[3];
	mpi::TaggedRecv (coord, 3, source, 
		REQUEST_GET_TAG, g.VCComm());
    	Int i = coord[0]; 
	Int j = coord[1];
	// we need the localwidth/height here,
	// used also to calculate numEntries
 
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
	NextIndexData (numEntries,
		    sendData_[source],
		    sendDataRequests_[source],
		    sendDataStatuses_[source]);

	DEBUG_ONLY (if
		(Int (sendData_[source][index].size ()) !=
		 numEntries) LogicError ("Error in NextIndex");)
	
	T *replyBuffer = sendData_[source][index].data ();
	
	for (Int t = 0; t < localWidth; ++t)
	{
	    T *sendCol = &replyBuffer[t * localHeight];
	    const T *XCol = Y.LockedBuffer (iLocalOffset, jLocalOffset + t);
	    MemCopy (sendCol, XCol, localHeight);
	}

	// Fire off non-blocking send
	mpi::TaggedISend (replyBuffer, numEntries, source, 
		DATA_GET_TAG, g.VCComm (), 
		sendDataRequests_[source][index]);
    }
}

// blocking update routines
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

    const Int XLDim = Z.LDim();

    const Int height = Z.Height();
    const Int width = Z.Width();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    std::vector<Int> dataindices_;
    dataindices_.resize (p);

    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();

    const T* XBuffer = Z.LockedBuffer();
 
    // prepost receives for coordinates
    for ( int rank = 0; rank < p; ++rank )
    {
	const Int index =
	    NextIndexCoord (recvCoord_[rank],
		    recvCoordRequests_[rank],
		    recvCoordStatuses_[rank]);
	    
	dataindices_[rank] = index;
	Int *coord_ = recvCoord_[rank][index].data();
	mpi::TaggedIRecv (coord_, 3, rank, COORD_PUT_TAG, g.VCComm(), 
		 recvCoordRequests_[rank][index]);
    } 
    
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    const Int colAlign = (Y.ColAlign() + i) % r;	    
    const Int rowAlign = (Y.RowAlign() + j) % c;

    const Int YLDim = Y.LDim();
    
    // send coordinates and data size
    for( Int step=0; step<p; ++step )
    {
	const Int colShift = Shift( receivingRow, colAlign, r );
	const Int rowShift = Shift( receivingCol, rowAlign, c );

	const Int localHeight = Length( height, colShift, r );
	const Int localWidth = Length( width, rowShift, c );

	const Int numEntries = localHeight * localWidth;
	    
	// target rank	    
	const Int destination = receivingRow + r*receivingCol;
	
	const Int index =
	NextIndexCoord (sendCoord_[destination],
		    sendCoordRequests_[destination],
		    sendCoordStatuses_[destination]);

	int * coord_ = sendCoord_[destination][index].data ();
	coord_[0] = i; 
	coord_[1] = j;
	coord_[2] = numEntries;

	// post receive for coordinates
	mpi::TaggedISend (coord_, 3, destination, 
		COORD_PUT_TAG, g.VCComm(), 
		sendCoordRequests_[destination][index]);

	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }

    // wait for my coordinates xfer to be over
    for (Int i = 0; i < p; ++i)
    {
	// coord receives
	const Int numRecvCoordRequests = recvCoordRequests_[i].size ();
	for (Int j = 0; j < numRecvCoordRequests; ++j)
	{
	    if (recvCoordStatuses_[i][j])
	    {
		mpi::Wait (recvCoordRequests_[i][j]);
		recvCoordStatuses_[i][j] = false;
	    }
	}
	// coord sends
	const Int numSendCoordRequests = sendCoordRequests_[i].size ();
	for (Int j = 0; j < numSendCoordRequests; ++j)
	{
	    if (sendCoordStatuses_[i][j])
	    {
		mpi::Wait (recvCoordRequests_[i][j]);
		sendCoordStatuses_[i][j] = false;
	    }
	}
    }

    // prepost receives for data
    // should be some way to get the index!
    for ( int rank = 0; rank < p; ++rank )
    {
	const int index = dataindices_[rank];
	const int i = recvCoord_[rank][index][0]; 
	const int j = recvCoord_[rank][index][1]; 
	const int numEntries = recvCoord_[rank][index][2];

	// post recv for data	
	if ( numEntries > 0 )
	{
	    const Int index =
		NextIndexData (numEntries,
		    recvData_[rank],
		    recvDataRequests_[rank],
		    recvDataStatuses_[rank]);

	    DEBUG_ONLY (if
		    (Int (recvData_[rank][index].size ()) !=
		     numEntries) LogicError ("Error in NextIndexData");)

            T *recvData = recvData_[rank][index].data ();

	    mpi::TaggedIRecv (recvData, numEntries, rank, 
		    DATA_PUT_TAG, g.VCComm (), 
		    recvDataRequests_[rank][index]);
	}
    } 
 
    // sends for data
    receivingRow = myProcessRow;
    receivingCol = myProcessCol;

    for( Int step=0; step<p; ++step )
    {
	const Int colShift = Shift( receivingRow, colAlign, r );
	const Int rowShift = Shift( receivingCol, rowAlign, c );

	const Int localHeight = Length( height, colShift, r );
	const Int localWidth = Length( width, rowShift, c );

	const Int numEntries = localHeight * localWidth;

	// send data 
	if( numEntries > 0  )
	{
	    // target rank	    
	    const Int destination = receivingRow + r*receivingCol;    

	    const Int index =
		NextIndexData (numEntries,
			sendData_[destination],
			sendDataRequests_[destination],
			sendDataStatuses_[destination]);

	    DEBUG_ONLY (if
		    (Int (sendData_[destination][index].size ()) !=
		     numEntries) LogicError ("Error in NextIndex");)

	    T *sendBuffer = sendData_[destination][index].data ();

	    for( Int t=0; t<localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t*localHeight];
		const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
		for( Int s=0; s<localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift+s*r];
	    }

	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
		    DATA_PUT_TAG, g.VCComm(), 
		    sendDataRequests_[destination][index]);
	}

	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;	    
    }

    // wait for my data xfer to be over
    for (Int i = 0; i < p; ++i)
    {
	// data receives
	const Int numrecvDataRequests = recvDataRequests_[i].size ();
	for (Int j = 0; j < numrecvDataRequests; ++j)
	{
	    if (recvDataStatuses_[i][j])
	    {
		mpi::Wait (recvDataRequests_[i][j]);
		recvDataStatuses_[i][j] = false;
	    }		    
	}
	// data sends
	const Int numsendDataRequests = sendDataRequests_[i].size ();
	for (Int j = 0; j < numsendDataRequests; ++j)
	{
	    if (sendDataStatuses_[i][j])
	    {
		mpi::Wait (sendDataRequests_[i][j]);
		sendDataStatuses_[i][j] = false;
	    }		    
	}
    }   
	
    // accumulate as data xfer is over
    // there must be a way to get index
    for ( int rank = 0; rank < p; ++rank )
    {
	const int index = dataindices_[rank];
	const int i = recvCoord_[rank][index][0]; 
	const int j = recvCoord_[rank][index][1]; 
	const int numEntries = recvCoord_[rank][index][2];
	
	// data recv'd, now accumulate	
	if ( numEntries > 0 )
	{
	    // Update Y
	    const T *Buffer = reinterpret_cast < const T * >(recvData_[rank][index].data());
	    
	    const Int colAlign = (Y.ColAlign () + i) % r;
	    const Int rowAlign = (Y.RowAlign () + j) % c;

	    const Int colShift = Shift (g.Row(), colAlign, r);
	    const Int rowShift = Shift (g.Col(), rowAlign, c);

	    const Int localHeight = Length (height, colShift, r);
	    const Int localWidth = Length (width, rowShift, c);

	    const Int iLocalOffset = Length (i, Y.ColShift (), r);
	    const Int jLocalOffset = Length (j, Y.RowShift (), c);

	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
		const T *XCol = &Buffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[s] = XCol[s];
	    }
	}
    } 

    dataindices_.clear();
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

    const Int XLDim = Z.LDim();	
     
    std::vector<Int> dataindices_;
    dataindices_.resize (p);

    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();

    const T* XBuffer = Z.LockedBuffer();
 
    // prepost receives for coordinates
    for ( Int rank = 0; rank < p; ++rank )
    {
	const Int index =
	    NextIndexCoord (recvCoord_[rank],
		    recvCoordRequests_[rank],
		    recvCoordStatuses_[rank]);
	    
	dataindices_[rank] = index;
	Int *coord_ = recvCoord_[rank][index].data();
	mpi::TaggedIRecv (coord_, 3, rank, COORD_PUT_TAG, g.VCComm(), 
		 recvCoordRequests_[rank][index]);
    } 
    
    // send coordinates
    for( Int rank=0; rank<p; ++rank )
    {
	const Int index =
	NextIndexCoord (sendCoord_[rank],
		    sendCoordRequests_[rank],
		    sendCoordStatuses_[rank]);

	int * coord_ = sendCoord_[rank][index].data ();
	coord_[0] = i; 
	coord_[1] = j;
	coord_[2] = -1;

	// post receive for coordinates
	mpi::TaggedISend (coord_, 3, rank, 
		COORD_PUT_TAG, g.VCComm(), 
		sendCoordRequests_[rank][index]);
    }

    // wait for my coordinates xfer to be over
    for (Int i = 0; i < p; ++i)
    {
	// coord receives
	const Int numRecvCoordRequests = recvCoordRequests_[i].size ();
	for (Int j = 0; j < numRecvCoordRequests; ++j)
	{
	    if (recvCoordStatuses_[i][j])
	    {
		mpi::Wait (recvCoordRequests_[i][j]);
		recvCoordStatuses_[i][j] = false;
	    }
	}
	// coord sends
	const Int numSendCoordRequests = sendCoordRequests_[i].size ();
	for (Int j = 0; j < numSendCoordRequests; ++j)
	{
	    if (sendCoordStatuses_[i][j])
	    {
		mpi::Wait (recvCoordRequests_[i][j]);
		sendCoordStatuses_[i][j] = false;
	    }
	}
    }
 
    // exchange data
    // data sends
    for( Int source=0; source<p; ++source )
    {
	const Int rindex = dataindices_[source];
	const Int i = recvCoord_[source][rindex][0]; 
	const Int j = recvCoord_[source][rindex][1]; 

	const Int myRow = g.Row ();
	const Int myCol = g.Col ();

	const Int colAlign = (X.ColAlign() + i) % r;
	const Int rowAlign = (X.RowAlign() + j) % c;

	const Int XLDim = Z.LDim();
	// local matrix width and height
	const Int height = Z.Height();
	const Int width = Z.Width();

	const Int colShift = Shift (myRow, colAlign, r);
	const Int rowShift = Shift (myCol, rowAlign, c);
	const Int localHeight = Length (height, colShift, r);
	const Int localWidth = Length (width, rowShift, c);

	const Int iLocalOffset = Length (i, X.ColShift (), r);
	const Int jLocalOffset = Length (j, X.RowShift (), c);

	const Int numEntries = localHeight * localWidth;

	DEBUG_ONLY (if (numEntries < Int (sizeof (T)))
		LogicError ("Count was too small");)

	T* XBuffer = Z.Buffer();
	const Int index =
	NextIndexData (numEntries,
		    sendData_[source],
		    sendDataRequests_[source],
		    sendDataStatuses_[source]);

	DEBUG_ONLY (if
		(Int (sendData_[source][index].size ()) !=
		 numEntries) LogicError ("Error in NextIndex");)
	
	T *replyBuffer = sendData_[source][index].data ();
	
	for (Int t = 0; t < localWidth; ++t)
	{
	    T *sendCol = &replyBuffer[t * localHeight];
	    const T *XCol = X.LockedBuffer (iLocalOffset, jLocalOffset + t);
	    MemCopy (sendCol, XCol, localHeight);
	}

	// Fire off non-blocking send
	mpi::TaggedISend (replyBuffer, numEntries, source, 
		DATA_GET_TAG, g.VCComm (), 
		sendDataRequests_[source][index]);
    }
 
    // data receives
    std::vector<T> recvVector_;

    // Receive all of the replies
    Int numReplies = 0;
    while ( numReplies < p )
    {
	mpi::Status status;
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
	}
	recvVector_.clear();
    }
    // wait for my data xfer to be over
    for (Int i = 0; i < p; ++i)
    {
	// data sends
	const Int numsendDataRequests = sendDataRequests_[i].size ();
	for (Int j = 0; j < numsendDataRequests; ++j)
	{
	    if (sendDataStatuses_[i][j])
	    {
		mpi::Wait (sendDataRequests_[i][j]);
		sendDataStatuses_[i][j] = false;
	    }		    
	}
    } 

    dataindices_.clear();

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

    const Int XLDim = Z.LDim();

    const Int height = Z.Height();
    const Int width = Z.Width();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    std::vector<Int> dataindices_;
    dataindices_.resize (p);

    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();

    const T* XBuffer = Z.LockedBuffer();
 
    // prepost receives for coordinates
    for ( int rank = 0; rank < p; ++rank )
    {
	const Int index =
	    NextIndexCoord (recvCoord_[rank],
		    recvCoordRequests_[rank],
		    recvCoordStatuses_[rank]);
	    
	dataindices_[rank] = index;
	Int *coord_ = recvCoord_[rank][index].data();
	mpi::TaggedIRecv (coord_, 3, rank, COORD_PUT_TAG, g.VCComm(), 
		 recvCoordRequests_[rank][index]);
    } 
    
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    const Int colAlign = (Y.ColAlign() + i) % r;	    
    const Int rowAlign = (Y.RowAlign() + j) % c;

    const Int YLDim = Y.LDim();
    
    // send coordinates and data size
    for( Int step=0; step<p; ++step )
    {
	const Int colShift = Shift( receivingRow, colAlign, r );
	const Int rowShift = Shift( receivingCol, rowAlign, c );

	const Int localHeight = Length( height, colShift, r );
	const Int localWidth = Length( width, rowShift, c );

	const Int numEntries = localHeight * localWidth;
	    
	// target rank	    
	const Int destination = receivingRow + r*receivingCol;
	
	const Int index =
	NextIndexCoord (sendCoord_[destination],
		    sendCoordRequests_[destination],
		    sendCoordStatuses_[destination]);

	int * coord_ = sendCoord_[destination][index].data ();
	coord_[0] = i; 
	coord_[1] = j;
	coord_[2] = numEntries;

	// post receive for coordinates
	mpi::TaggedISend (coord_, 3, destination, 
		COORD_PUT_TAG, g.VCComm(), 
		sendCoordRequests_[destination][index]);

	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }

    // wait for my coordinates xfer to be over
    for (Int i = 0; i < p; ++i)
    {
	// coord receives
	const Int numRecvCoordRequests = recvCoordRequests_[i].size ();
	for (Int j = 0; j < numRecvCoordRequests; ++j)
	{
	    if (recvCoordStatuses_[i][j])
	    {
		mpi::Wait (recvCoordRequests_[i][j]);
		recvCoordStatuses_[i][j] = false;
	    }
	}
	// coord sends
	const Int numSendCoordRequests = sendCoordRequests_[i].size ();
	for (Int j = 0; j < numSendCoordRequests; ++j)
	{
	    if (sendCoordStatuses_[i][j])
	    {
		mpi::Wait (recvCoordRequests_[i][j]);
		sendCoordStatuses_[i][j] = false;
	    }
	}
    }

    // prepost receives for data
    // should be some way to get the index!
    for ( int rank = 0; rank < p; ++rank )
    {
	const int index = dataindices_[rank];
	const int i = recvCoord_[rank][index][0]; 
	const int j = recvCoord_[rank][index][1]; 
	const int numEntries = recvCoord_[rank][index][2];

	// post recv for data	
	if ( numEntries > 0 )
	{
	    const Int index =
		NextIndexData (numEntries,
		    recvData_[rank],
		    recvDataRequests_[rank],
		    recvDataStatuses_[rank]);

	    DEBUG_ONLY (if
		    (Int (recvData_[rank][index].size ()) !=
		     numEntries) LogicError ("Error in NextIndexData");)

            T *recvData = recvData_[rank][index].data ();

	    mpi::TaggedIRecv (recvData, numEntries, rank, 
		    DATA_PUT_TAG, g.VCComm (), 
		    recvDataRequests_[rank][index]);
	}
    } 
 
    // sends for data
    receivingRow = myProcessRow;
    receivingCol = myProcessCol;

    for( Int step=0; step<p; ++step )
    {
	const Int colShift = Shift( receivingRow, colAlign, r );
	const Int rowShift = Shift( receivingCol, rowAlign, c );

	const Int localHeight = Length( height, colShift, r );
	const Int localWidth = Length( width, rowShift, c );

	const Int numEntries = localHeight * localWidth;

	// send data 
	if( numEntries > 0  )
	{
	    // target rank	    
	    const Int destination = receivingRow + r*receivingCol;    

	    const Int index =
		NextIndexData (numEntries,
			sendData_[destination],
			sendDataRequests_[destination],
			sendDataStatuses_[destination]);

	    DEBUG_ONLY (if
		    (Int (sendData_[destination][index].size ()) !=
		     numEntries) LogicError ("Error in NextIndex");)

	    T *sendBuffer = sendData_[destination][index].data ();

	    for( Int t=0; t<localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t*localHeight];
		const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
		for( Int s=0; s<localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift+s*r];
	    }

	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
		    DATA_PUT_TAG, g.VCComm(), 
		    sendDataRequests_[destination][index]);
	}

	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;	    
    }

    // wait for my data xfer to be over
    for (Int i = 0; i < p; ++i)
    {
	// data receives
	const Int numrecvDataRequests = recvDataRequests_[i].size ();
	for (Int j = 0; j < numrecvDataRequests; ++j)
	{
	    if (recvDataStatuses_[i][j])
	    {
		mpi::Wait (recvDataRequests_[i][j]);
		recvDataStatuses_[i][j] = false;
	    }		    
	}
	// data sends
	const Int numsendDataRequests = sendDataRequests_[i].size ();
	for (Int j = 0; j < numsendDataRequests; ++j)
	{
	    if (sendDataStatuses_[i][j])
	    {
		mpi::Wait (sendDataRequests_[i][j]);
		sendDataStatuses_[i][j] = false;
	    }		    
	}
    }   
	
    // accumulate as data xfer is over
    // there must be a way to get index
    for ( int rank = 0; rank < p; ++rank )
    {
	const int index = dataindices_[rank];
	const int i = recvCoord_[rank][index][0]; 
	const int j = recvCoord_[rank][index][1]; 
	const int numEntries = recvCoord_[rank][index][2];
	
	// data recv'd, now accumulate	
	if ( numEntries > 0 )
	{
	    // Update Y
	    const T *Buffer = reinterpret_cast < const T * >(recvData_[rank][index].data());
	    
	    const Int colAlign = (Y.ColAlign () + i) % r;
	    const Int rowAlign = (Y.RowAlign () + j) % c;

	    const Int colShift = Shift (g.Row(), colAlign, r);
	    const Int rowShift = Shift (g.Col(), rowAlign, c);

	    const Int localHeight = Length (height, colShift, r);
	    const Int localWidth = Length (width, rowShift, c);

	    const Int iLocalOffset = Length (i, Y.ColShift (), r);
	    const Int jLocalOffset = Length (j, Y.RowShift (), c);

	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
		const T *XCol = &Buffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[s] += XCol[s];
	    }
	}
    } 

    dataindices_.clear();
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

    sendDataStatuses_.clear ();
    sendCoordStatuses_.clear ();
    recvDataStatuses_.clear ();
    recvCoordStatuses_.clear ();

    sendDataRequests_.clear ();
    sendCoordRequests_.clear ();
    recvDataRequests_.clear ();
    recvCoordRequests_.clear ();

    sendData_.clear ();
    sendCoord_.clear ();
    recvData_.clear ();
    recvCoord_.clear ();

    mpi::WindowUnlock (put_win_);
    mpi::WindowFree (put_win_);

    mpi::WindowUnlock (acc_win_);
    mpi::WindowFree (acc_win_);
    
    mpi::WindowUnlock (getrq_win_);
    mpi::WindowFree (getrq_win_);

    delete put_win_base_;
    delete acc_win_base_;
    delete getrq_win_base_;
}

template class AxpyInterface2<Int>;
template class AxpyInterface2<float>;
template class AxpyInterface2<double>;
template class AxpyInterface2<Complex<float>>;
template class AxpyInterface2<Complex<double>>;

} // namespace El
#endif
