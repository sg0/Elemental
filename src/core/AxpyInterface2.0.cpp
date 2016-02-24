/*
   Copyright (c) 2009-2015, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   All rights reserved.

   Authors:
   This interface is mainly due to Martin Schatz, but it was put into its
   current form by Jack Poulson.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El
{
template<typename T>
AxpyInterface2<T>::AxpyInterface2()
    : attached_( false ), detached_( true ), 
      GlobalMat_( 0 ), sendDummy_( 0 ), recvDummy_( 0 )
{ }

template<typename T>
AxpyInterface2<T>::AxpyInterface2( DistMatrix<T>& Z )
    : sendDummy_( 0 ), recvDummy_( 0 )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::AxpyInterface2" ) )

    if (detached_)
    {
	attached_ = true;
	detached_ = false;
	GlobalMat_ = &Z;

	const Int p = Z.Grid().Size();

	// eom
	sentEomTo_.resize( p, false );
	haveEomFrom_.resize( p, false );
	eomSendRequests_.resize( p );

	// send data
	sendingData_.resize( p );
	sendingRequest_.resize( p );
	sendingReply_.resize( p );

	// send requests
	dataSendRequests_.resize( p );
	requestSendRequests_.resize( p );
	replySendRequests_.resize( p );

	// recv data
	dataVectors_.resize( p );
	requestVectors_.resize( p );
	replyVectors_.resize( p );
    }
}

template <typename T> AxpyInterface2 <T>::~AxpyInterface2()
{
    if( attached_ && !detached_ )
    {
        if( std::uncaught_exception() )
        {
            const Grid& g = GlobalMat_->Grid();
            std::ostringstream os;
            os << g.Rank()
               <<
               "Uncaught exception detected during AxpyInterface destructor "
               "that required a call to Detach. Instead of allowing for the "
               "possibility of Detach throwing another exception and "
               "resulting in a 'terminate', we instead immediately dump the "
               "call stack (if not in RELEASE mode) since the program will "
               "likely hang:" << std::endl;
            std::cerr << os.str();
            DEBUG_ONLY( DumpCallStack() )
        }
        else
            Detach();
    }
}

template <typename T>
void AxpyInterface2 <T>::Attach( DistMatrix <T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Attach" ) )

    if (detached_)
    {
	attached_ = true;
	detached_ = false;
	GlobalMat_ = &Z;

	const Int p = Z.Grid().Size();

	// eom
	sentEomTo_.resize( p, false );
	haveEomFrom_.resize( p, false );
	eomSendRequests_.resize( p );

	// send data
	sendingData_.resize( p );
	sendingRequest_.resize( p );
	sendingReply_.resize( p );

	// send requests
	dataSendRequests_.resize( p );
	requestSendRequests_.resize( p );
	replySendRequests_.resize( p );

	// recv data
	dataVectors_.resize( p );
	requestVectors_.resize( p );
	replyVectors_.resize( p );
    }
}

template <typename T> bool AxpyInterface2 <T>::Finished()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Finished" );

    if( !attached ) LogicError( "Not attached" ); )
    
    const Grid& g = GlobalMat_->Grid();

    const Int p = g.Size();
    bool finished = true;

    for( Int rank = 0; rank < p; ++rank )
    {
        if( !sentEomTo_[rank] || !haveEomFrom_[rank] )
        {
            finished = false;
            break;
        }
    }

    return finished;
}

template <typename T> void AxpyInterface2 <T>::HandleEoms()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::HandleEoms" ) )
        
    const Grid& g = GlobalMat_->Grid();
    const Int p = g.Size();
   
    // progress
    UpdateRequestStatuses();

    // Try to progress our EOM sends
    for( Int i = 0; i < p; ++i )
    {
        if( !sentEomTo_[i] )
        {
            bool shouldSendEom = true;
            const Int numSends = sendingData_[i].size();

            for( Int j = 0; j < numSends; ++j )
            {
                if( sendingData_[i][j] )
                {
                    shouldSendEom = false;
                    break;
                }
            }

            const Int numRequests = sendingRequest_[i].size();

            for( Int j = 0; j < numRequests; ++j )
            {
                if( !shouldSendEom || sendingRequest_[i][j] )
                {
                    shouldSendEom = false;
                    break;
                }
            }

            const Int numReplies = sendingReply_[i].size();

            for( Int j = 0; j < numReplies; ++j )
            {
                if( !shouldSendEom || sendingReply_[i][j] )
                {
                    shouldSendEom = false;
                    break;
                }
            }

            if( shouldSendEom )
            {
                mpi::Request& request = eomSendRequests_[i];
                mpi::TaggedISend
                ( &sendDummy_, 1, i, EOM_TAG, g.VCComm(), request );
                sentEomTo_[i] = true;
            }
        }
    }

    mpi::Status status;

    if( mpi::IProbe( mpi::ANY_SOURCE, EOM_TAG, g.VCComm(), status ) )
    {
        const Int source = status.MPI_SOURCE;
        mpi::TaggedRecv( &recvDummy_, 1, source, EOM_TAG, g.VCComm() );
        haveEomFrom_[source] = true;
    }
}

template <typename T> void AxpyInterface2 <T>::HandleLocalToGlobalData( TransferType ttype )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::HandleLocalToGlobalData" ) )

    DistMatrix <T>& Y = *GlobalMat_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRow = g.Row();
    const Int myCol = g.Col();
    mpi::Status status;

    if( mpi::IProbe( mpi::ANY_SOURCE, DATA_TAG, g.VCComm(), status ) )
    {
        // Message exists, so recv and pack
        const Int count = mpi::GetCount <byte> ( status );

        DEBUG_ONLY( if( count < Int( 3 * sizeof( Int ) + sizeof( T ) ) )
                    LogicError( "Count was too small" ); )
        const Int source = status.MPI_SOURCE;

        recvVector_.resize( count );
        byte* recvBuffer = recvVector_.data();
        mpi::TaggedRecv( recvBuffer, count, source, DATA_TAG, g.VCComm() );
        // Extract the header
        byte* head = recvBuffer;
        const Int i = *reinterpret_cast <const Int*>( head );
        head += sizeof( Int );
        const Int j = *reinterpret_cast <const Int*>( head );
        head += sizeof( Int );
        const Int height = *reinterpret_cast <const Int*>( head );
        head += sizeof( Int );
        const Int width = *reinterpret_cast <const Int*>( head );
        head += sizeof( Int );

        DEBUG_ONLY( if( height < 0 || width < 0 )
                    RuntimeError
                    ( "Unpacked heights were negative:\n",
                      "  i=     ", i, std::hex, "(", i, ")\n", std::dec,
                      "  j=     ", j, std::hex, "(", j, ")\n", std::dec,
                      "  height=", height, std::hex, "(", height, ")\n",
                      std::dec, "  width= ", width, std::hex, "(", width,
                      ")\n" );
                    if( i < 0
                            || j <
                            0 ) RuntimeError( "Unpacked offsets were negative:\n",
                                              "  i=     ", i, std::hex, "(", i,
                                              ")\n", std::dec, "  j=     ", j,
                                              std::hex, "(", j, ")\n", std::dec,
                                              "  height=", height, std::hex, "(",
                                              height, ")\n", std::dec, "  width= ",
                                              width, std::hex, "(", width, ")\n" );
                     if( i + height > Y.Height()
                                || j + width >
                                Y.Width() )RuntimeError
                            ( "Unpacked submatrix was out of bounds:\n", "  i=     ",
                              i, std::hex, "(", i, ")\n", std::dec, "  j=     ", j,
                              std::hex, "(", j, ")\n", std::dec, "  height=", height,
                              std::hex, "(", height, ")\n", std::dec, "  width= ",
                              width, std::hex, "(", width, ")\n" ); )
                    
	// Update Y
        const T* XBuffer = reinterpret_cast <const T*>( head );

        const Int colAlign = ( Y.ColAlign() + i ) % r;
        const Int rowAlign = ( Y.RowAlign() + j ) % c;
        const Int colShift = Shift( myRow, colAlign, r );
        const Int rowShift = Shift( myCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int iLocalOffset = Length( i, Y.ColShift(), r );
        const Int jLocalOffset = Length( j, Y.RowShift(), c );

	switch( ttype )
	{
	    case ACC:
		{
		    for( Int t = 0; t < localWidth; ++t )
		    {
			T* YCol = Y.Buffer( iLocalOffset, jLocalOffset + t );
			const T* XCol = &XBuffer[t * localHeight];

			for( Int s = 0; s < localHeight; ++s )
			    YCol[s] += XCol[s];
		    }
		    break;
		}
	    case PUT:
		{
		    for( Int t = 0; t < localWidth; ++t )
		    {
			T* YCol = Y.Buffer( iLocalOffset, jLocalOffset + t );
			const T* XCol = &XBuffer[t * localHeight];

			MemCopy<T>( YCol, XCol, localHeight );
		    }
		    break;
		}
	    default:
		LogicError( "Illegal Transfer type" );
	}

        // Free the memory for the recv buffer
        recvVector_.clear();
    }
}

template <typename T>
void AxpyInterface2 <T>::HandleGlobalToLocalRequest()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::HandleGlobalToLocalRequest" ) )

    DistMatrix <T>& X = *GlobalMat_;
    const Grid& g = X.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRow = g.Row();
    const Int myCol = g.Col();
    mpi::Status status;

    if( mpi::IProbe( mpi::ANY_SOURCE, DATA_REQUEST_TAG, g.VCComm(), status ) )
    {
        // Request exists, so recv
        const Int source = status.MPI_SOURCE;
        const Int recvSize = 4 * sizeof( Int );
        recvVector_.resize( recvSize );
        byte* recvBuffer = recvVector_.data();
        mpi::TaggedRecv
        ( recvBuffer, recvSize, source, DATA_REQUEST_TAG, g.VCComm() );
        // Extract the header
        const byte* recvHead = recvBuffer;
        const Int i = *reinterpret_cast <const Int*>( recvHead );
        recvHead += sizeof( Int );
        const Int j = *reinterpret_cast <const Int*>( recvHead );
        recvHead += sizeof( Int );
        const Int height = *reinterpret_cast <const Int*>( recvHead );
        recvHead += sizeof( Int );
        const Int width = *reinterpret_cast <const Int*>( recvHead );
        recvHead += sizeof( Int );
        const Int colAlign = ( X.ColAlign() + i ) % r;
        const Int rowAlign = ( X.RowAlign() + j ) % c;
        const Int colShift = Shift( myRow, colAlign, r );
        const Int rowShift = Shift( myCol, rowAlign, c );
        const Int iLocalOffset = Length( i, X.ColShift(), r );
        const Int jLocalOffset = Length( j, X.RowShift(), c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;
        const Int bufferSize = 2 * sizeof( Int ) + numEntries * sizeof( T );
        const Int index = ReadyForSend( bufferSize, replyVectors_[source],
                                        replySendRequests_[source],
                                        sendingReply_[source] );
        // Pack the reply header
        byte* sendBuffer = replyVectors_[source][index].data();
        byte* sendHead = sendBuffer;
        *reinterpret_cast <Int*>( sendHead ) = myRow;
        sendHead += sizeof( Int );
        *reinterpret_cast <Int*>( sendHead ) = myCol;
        sendHead += sizeof( Int );
        // Pack the payload
        T* sendData = reinterpret_cast <T*>( sendHead );

        for( Int t = 0; t < localWidth; ++t )
        {
            T* sendCol = &sendData[t * localHeight];
            const T* XCol = X.LockedBuffer( iLocalOffset, jLocalOffset + t );
            MemCopy( sendCol, XCol, localHeight );
        }

        // Fire off non-blocking send
        mpi::TaggedISend
        ( sendBuffer, bufferSize, source, DATA_REPLY_TAG, g.VCComm(),
          replySendRequests_[source][index] );
    }
}

template <typename T>
void AxpyInterface2 <T>::Acc( Matrix <T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Acc" ) )

    if( attached_ )
        LocalToGlobal( Z, i, j );
    else
        LogicError( "Cannot issue Acc before attaching DistMatrix." );
}

template <typename T>
void AxpyInterface2 <T>::Put( Matrix <T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Acc" ) )

    if( attached_ )
        LocalToGlobal( Z, i, j );
    else
	LogicError( "Cannot issue Put before attaching DistMatrix." );
}

template <typename T>
void AxpyInterface2 <T>::Get( Matrix <T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Get" ) )

    if( attached_ )
        GlobalToLocal( Z, i, j );
    else
	LogicError( "Cannot issue Get before attaching DistMatrix." );
}

// Update Y(i:i+height-1,j:j+width-1) += alpha X, where X is height x width
template <typename T>
void AxpyInterface2 <T>::LocalToGlobal( Matrix <T>& X, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::LocalToGlobal" ) )
    DistMatrix <T>& Y = *GlobalMat_;

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( i + X.Height() > Y.Height() || j + X.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int height = X.Height();
    const Int width = X.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    for( Int step = 0; step < p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;

        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r * receivingCol;
            const Int bufferSize =
                3 * sizeof( Int ) + ( numEntries + 1 ) * sizeof( T );
            const Int index =
                ReadyForSend( bufferSize, dataVectors_[destination],
                              dataSendRequests_[destination],
                              sendingData_[destination] );

            DEBUG_ONLY( if
                        ( Int( dataVectors_[destination][index].size() ) !=
                          bufferSize ) LogicError( "Error in ReadyForSend" ); )
            // Pack the header
            byte* sendBuffer = dataVectors_[destination][index].data();

            byte* head = sendBuffer;
            *reinterpret_cast <Int*>( head ) = i;
            head += sizeof( Int );
            *reinterpret_cast <Int*>( head ) = j;
            head += sizeof( Int );
            *reinterpret_cast <Int*>( head ) = height;
            head += sizeof( Int );
            *reinterpret_cast <Int*>( head ) = width;
            head += sizeof( Int );
            // Pack the payload
            T* sendData = reinterpret_cast <T*>( head );
            const T* XBuffer = X.LockedBuffer();
            const Int XLDim = X.LDim();

            for( Int t = 0; t < localWidth; ++t )
            {
                T* thisSendCol = &sendData[t * localHeight];
                const T* thisXCol = &XBuffer[( rowShift + t * c ) * XLDim];

                for( Int s = 0; s < localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift + s * r];
            }

            // Fire off the locally-blocking send
            mpi::TaggedISend
            ( sendBuffer, bufferSize, destination, DATA_TAG, g.VCComm(),
              dataSendRequests_[destination][index] );
	    mpi::Wait( dataSendRequests_[destination][index] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

// Update Y += alpha X(i:i+height-1,j:j+width-1), where X is the dist-matrix
template <typename T>
void AxpyInterface2 <T>::GlobalToLocal( Matrix <T>& Y, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::GlobalToLocal" ) )

    const DistMatrix <T>& X = *GlobalMat_;
    const Int height = Y.Height();
    const Int width = Y.Width();

    if( i + height > X.Height() || j + width > X.Width() )
        LogicError( "Invalid AxpyGlobalToLocal submatrix" );

    const Grid& g = X.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    // Send out the requests to all processes in the grid
    for( Int rank = 0; rank < p; ++rank )
    {
        const Int bufferSize = 4 * sizeof( Int );
        const Int index = ReadyForSend( bufferSize, requestVectors_[rank],
                                        requestSendRequests_[rank],
                                        sendingRequest_[rank] );
        // Copy the request header into the send buffer
        byte* sendBuffer = requestVectors_[rank][index].data();
        byte* head = sendBuffer;
        *reinterpret_cast <Int*>( head ) = i;
        head += sizeof( Int );
        *reinterpret_cast <Int*>( head ) = j;
        head += sizeof( Int );
        *reinterpret_cast <Int*>( head ) = height;
        head += sizeof( Int );
        *reinterpret_cast <Int*>( head ) = width;
        head += sizeof( Int );
        // Begin the non-blocking send
        mpi::TaggedISend
        ( sendBuffer, bufferSize, rank, DATA_REQUEST_TAG, g.VCComm(),
          requestSendRequests_[rank][index] );
    }

    // Receive all of the replies
    Int numReplies = 0;

    while( numReplies < p )
    {
        HandleGlobalToLocalRequest();
        mpi::Status status;

        if( mpi::IProbe
            ( mpi::ANY_SOURCE, DATA_REPLY_TAG, g.VCComm(), status ) )
        {
            const Int source = status.MPI_SOURCE;
            // Ensure that we have a recv buffer
            const Int count = mpi::GetCount <byte> ( status );
            recvVector_.resize( count );
            byte* recvBuffer = recvVector_.data();
            // Receive the data
            mpi::TaggedRecv
            ( recvBuffer, count, source, DATA_REPLY_TAG, g.VCComm() );
            // Unpack the reply header
            const byte* head = recvBuffer;
            const Int row = *reinterpret_cast <const Int*>( head );
            head += sizeof( Int );
            const Int col = *reinterpret_cast <const Int*>( head );
            head += sizeof( Int );
            const T* recvData = reinterpret_cast <const T*>( head );
            // Compute the local heights and offsets
            const Int colAlign = ( X.ColAlign() + i ) % r;
            const Int rowAlign = ( X.RowAlign() + j ) % c;
            const Int colShift = Shift( row, colAlign, r );
            const Int rowShift = Shift( col, rowAlign, c );
            const Int localHeight = Length( height, colShift, r );
            const Int localWidth = Length( width, rowShift, c );

            // Unpack the local matrix
            for( Int t = 0; t < localWidth; ++t )
            {
                T* YCol = Y.Buffer( 0, rowShift + t * c );
                const T* XCol = &recvData[t * localHeight];

                for( Int s = 0; s < localHeight; ++s )
                    YCol[colShift + s * r] = XCol[s];
            }

            ++numReplies;
        }
    }
}

template <typename T>
Int AxpyInterface2 <T>::ReadyForSend
( Int sendSize,
  std::deque <std::vector <byte>>& sendVectors,
  std::deque <mpi::Request>& requests,
  std::deque <bool>& requestStatuses )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::ReadyForSend" ) )
    const Int numCreated = sendVectors.size();

    DEBUG_ONLY( if( numCreated != Int( requests.size() ) ||
                    numCreated !=
                    Int( requestStatuses.size() ) )LogicError
                ( "size mismatch" ); )
    
    for( Int i = 0; i < numCreated; ++i )
    {
	// If this request is still running, test to see if it finished.
	if( requestStatuses[i] )
	{
	    const bool finished = mpi::Test( requests[i] );
	    requestStatuses[i] = !finished;
	}

	if( !requestStatuses[i] )
	{
	    requestStatuses[i] = true;
	    sendVectors[i].resize( sendSize );
	    return i;
	}
    }

    sendVectors.resize( numCreated + 1 );
    sendVectors[numCreated].resize( sendSize );
    requests.push_back( mpi::REQUEST_NULL );
    requestStatuses.push_back( true );
    return numCreated;
}

template <typename T> bool AxpyInterface2 <T>::ReturnRequestStatuses()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::ReturnRequestStatuses" ) )
    const Grid& g = GlobalMat_->Grid();
    const Int p = g.Size();

    for( Int i = 0; i < p; ++i )
    {
        const Int numDataSendRequests = dataSendRequests_[i].size();

        for( Int j = 0; j < numDataSendRequests; ++j )
        {
            if( sendingData_[i][j] )
                sendingData_[i][j] = !mpi::Test( dataSendRequests_[i][j] );

            if( sendingData_[i][j] )
                return false;
        }

        const Int numRequestSendRequests = requestSendRequests_[i].size();

        for( Int j = 0; j < numRequestSendRequests; ++j )
        {
            if( sendingRequest_[i][j] )
                sendingRequest_[i][j] = !mpi::Test( requestSendRequests_[i][j] );

            if( sendingRequest_[i][j] )
                return false;
        }

        const Int numReplySendRequests = replySendRequests_[i].size();

        for( Int j = 0; j < numReplySendRequests; ++j )
        {
            if( sendingReply_[i][j] )
                sendingReply_[i][j] = !mpi::Test( replySendRequests_[i][j] );

            if( sendingReply_[i][j] )
                return false;
        }
    }

    return true;
}

template <typename T> void AxpyInterface2 <T>::UpdateRequestStatuses()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::UpdateRequestStatuses" ) )
	
    const Grid& g = GlobalMat_->Grid();
    const Int p = g.Size();

    for( Int i = 0; i < p; ++i )
    {
        const Int numDataSendRequests = dataSendRequests_[i].size();

        for( Int j = 0; j < numDataSendRequests; ++j )
            if( sendingData_[i][j] )
                sendingData_[i][j] = !mpi::Test( dataSendRequests_[i][j] );

        const Int numRequestSendRequests = requestSendRequests_[i].size();

        for( Int j = 0; j < numRequestSendRequests; ++j )
            if( sendingRequest_[i][j] )
                sendingRequest_[i][j] = !mpi::Test( requestSendRequests_[i][j] );

        const Int numReplySendRequests = replySendRequests_[i].size();

        for( Int j = 0; j < numReplySendRequests; ++j )
            if( sendingReply_[i][j] )
                sendingReply_[i][j] = !mpi::Test( replySendRequests_[i][j] );
    }
}

// Remote completion
template <typename T> void AxpyInterface2 <T>::Flush( TransferType ttype )
{
    switch( ttype )
    {
	case ACC:
	    {
		HandleLocalToGlobalData( ACC );
		break;
	    }
	case PUT:
	    {
		HandleLocalToGlobalData( PUT );
		break;
	    }
	case GET:
	    {
		HandleGlobalToLocalRequest();
		break;
	    }
	default:
	    LogicError( "Illegal transfer type" );
    }
}

template <typename T> void AxpyInterface2 <T>::Detach()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Detach" ) )

    if ( detached_ )
	return;
    if( !attached_  )
        LogicError( "Must attach before detaching." );

    const Grid& g = GlobalMat_->Grid();
 
    while( !Finished() )
    {
	Flush( PUT );
	Flush( ACC );
	Flush( GET );

        HandleEoms();
    }   
    
    attached_ = false;
    detached_ = true;

    sentEomTo_.clear();
    haveEomFrom_.clear();
    eomSendRequests_.clear();
    
    sendingData_.clear();
    sendingRequest_.clear();
    sendingReply_.clear();
    
    dataVectors_.clear();
    recvVector_.clear();
    requestVectors_.clear();
    replyVectors_.clear();
    
    dataSendRequests_.clear();
    requestSendRequests_.clear();
    replySendRequests_.clear();

    mpi::Barrier( g.VCComm() );
}

#define PROTO(T) template class AxpyInterface2<T>;
#include "El/macros/Instantiate.h"

} // namespace El
