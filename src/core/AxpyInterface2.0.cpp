/*
This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include <cassert>

namespace El
{
template<typename T>
AxpyInterface2<T>::AxpyInterface2()
    : GlobalArrayPut_( 0 ), GlobalArrayGet_( 0 ),
      matrices_( 0 ), coords_( 0 ), dataVectors_( 0 ),
      toBeAttachedForGet_( false ), toBeAttachedForPut_( false ),
      attached_( false ), detached_( true )
{ }

template<typename T>
AxpyInterface2<T>::AxpyInterface2( DistMatrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::AxpyInterface2" ) )
    attached_ 		= false;
    detached_ 		= true;
    toBeAttachedForGet_ = false;
    toBeAttachedForPut_ = false;
    GlobalArrayPut_ 	= 0;
    GlobalArrayGet_ 	= 0;
}

template<typename T>
AxpyInterface2<T>::AxpyInterface2( const DistMatrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::AxpyInterface2" ) )
 
    attached_ 		= false;
    detached_ 		= true;
    
    toBeAttachedForGet_ = false;
    toBeAttachedForPut_ = false;
    
    GlobalArrayPut_ 	= 0;
    GlobalArrayGet_ 	= 0;
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
        DEBUG_ONLY( DumpCallStack() )
    }
    else
        Detach();
}

template<typename T>
Int AxpyInterface2<T>::NextIndexData(
    Int target,
    Int dataSize,
    const void* base_address,
    Int* mindex )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::NextIndexData" ) )
    assert( base_address != NULL );
    Int matrixIndex = 0;
    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    const Int numMatrices = matrices_.size();

    // search for matrix base
    for( Int m = 0; m < numMatrices; m++ )
    {
        if( matrices_[m].base_ == base_address )
        {
            matrixIndex = m;
            break;
        }

        // uninitiated, first time
        if( matrices_[m].base_ == NULL )
        {
            matrices_[m].base_ = base_address;
            matrixIndex = m;
            break;
        }

        matrixIndex = m+1;
    }

    // need to create new object
    if( matrixIndex == numMatrices )
    {
        struct matrix_params_ mp;
        mp.data_.resize( p );
        mp.requests_.resize( p );
        mp.statuses_.resize( p );
        mp.base_ = NULL;
        // push back new matrix_params created
        // with default constructor
        matrices_.push_back( mp );
        matrices_[matrixIndex].base_ = base_address;
    }

    // go through the request, data,
    // status objects
    const Int numCreated = matrices_[matrixIndex].data_[target].size();

    DEBUG_ONLY( if( numCreated != Int( matrices_[matrixIndex].requests_[target].size() )
                    ||  numCreated != Int( matrices_[matrixIndex].statuses_[target].size() ) )
                LogicError( "size mismatch" ); )

    for( Int i = 0; i < numCreated; ++i )
    {
	// If this request is still running,
	// test to see if it finished.
	if( matrices_[matrixIndex].statuses_[target][i] )
	{
	    const bool finished = mpi::Test( matrices_[matrixIndex].requests_[target][i] );
	    matrices_[matrixIndex].statuses_[target][i] = !finished;
	}

        if( !matrices_[matrixIndex].statuses_[target][i] )
        {
	    matrices_[matrixIndex].statuses_[target][i] = true;
            matrices_[matrixIndex].data_[target][i].resize( dataSize );
            *mindex = matrixIndex;
            return i;
        }
    }

    matrices_[matrixIndex].data_[target].resize( numCreated + 1 );
    matrices_[matrixIndex].data_[target][numCreated].resize( dataSize );
    matrices_[matrixIndex].requests_[target].push_back( mpi::REQUEST_NULL );
    matrices_[matrixIndex].statuses_[target].push_back( true );
    *mindex = matrixIndex;
    return numCreated;
}

template<typename T>
Int AxpyInterface2<T>::NextIndexCoord(
    Int i, Int j,
    Int target,
    const void* base_address,
    Int* cindex )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::NextIndexCoord" ) )
    assert( base_address != NULL );
    Int coordIndex = 0;
    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    const Int numCoords = coords_.size();

    // search for matrix base
    for( Int m = 0; m < numCoords; m++ )
    {
        if( coords_[m].base_ == base_address )
        {
            coordIndex = m;
            break;
        }

        if( coords_[m].base_ == NULL )
        {
            coords_[m].base_ = base_address;
            coordIndex = m;
            break;
        }

        coordIndex = m+1;
    }

    // need to create new object
    if( coordIndex == numCoords )
    {
        struct coord_params_ cp;
        cp.coord_.resize( p );
        cp.requests_.resize( p );
        cp.statuses_.resize( p );
        cp.base_ = NULL;
        // push back new matrix_params created
        // with default constructor
        coords_.push_back( cp );
        coords_[coordIndex].base_ = base_address;
    }

    // go through the request, data,
    // status objects
    const Int numCreated = coords_[coordIndex].coord_[target].size();

    DEBUG_ONLY( if( numCreated != Int( coords_[coordIndex].requests_[target].size() )
                    ||  numCreated != Int( matrices_[coordIndex].statuses_[target].size() ) )
                LogicError( "size mismatch" ); )
    
    for( Int i = 0; i < numCreated; ++i )
    {
	// If this request is still running,
        // test to see if it finished.
        if( coords_[coordIndex].statuses_[target][i] )
	{
	    const bool finished = mpi::Test( coords_[coordIndex].requests_[target][i] );
	    coords_[coordIndex].statuses_[target][i] = !finished;
        }

        if( !coords_[coordIndex].statuses_[target][i] )
        {
	    coords_[coordIndex].statuses_[target][i] = true;
            coords_[coordIndex].coord_[target][i][0] = i;
            coords_[coordIndex].coord_[target][i][1] = j;
            *cindex = coordIndex;
            return i;
        }
    }

    coords_[coordIndex].coord_[target].resize( numCreated + 1 );
    coords_[coordIndex].coord_[target][numCreated][0] = i;
    coords_[coordIndex].coord_[target][numCreated][1] = j;
    coords_[coordIndex].requests_[target].push_back( mpi::REQUEST_NULL );
    coords_[coordIndex].statuses_[target].push_back( true );
    *cindex = coordIndex;
    return numCreated;
}

template<typename T>
void AxpyInterface2<T>::Attach( DistMatrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Attach" ) )

    // attached_ will be only set in Attach
    // and only unset in Detach
    if( !attached_ && detached_ )
    {
        attached_ = true;
        detached_ = false;
    }
    else
        LogicError( "Must detach before reattaching." );

    const Grid& g = Z.Grid();
    const Int p = g.Size();

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

        if( dataVectors_.empty() )
            dataVectors_.resize( p );

        if( matrices_.empty() )
        {
            struct matrix_params_ mp;
            mp.data_.resize( p );
            mp.requests_.resize( p );
            mp.statuses_.resize( p );
            mp.base_ = NULL;
            // push back new matrix_params created
            // with default constructor
            matrices_.push_back( mp );
        }

        if( coords_.empty() )
        {
            struct coord_params_ cp;
            cp.coord_.resize( p );
            cp.requests_.resize( p );
            cp.statuses_.resize( p );
            cp.base_ = NULL;
            // push back new matrix_params created
            // with default constructor
            coords_.push_back( cp );
        }
    }

    mpi::Barrier( g.VCComm() );
}

template<typename T>
void AxpyInterface2<T>::Attach( const DistMatrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Attach" ) )

    // attached_ will be only set in Attach
    // and only unset in Detach
    if( !attached_ && detached_ )
    {
        attached_ = true;
        detached_ = false;
    }
    else
        LogicError( "Must detach before reattaching." );

    const Grid& g = Z.Grid();
    const Int p = g.Size();

    // the matrix base_ is not known until
    // an update operation (put/get/acc)
    // so it is kept blank
    // if DistMatrix is non-const, all one-sided
    // transfers -- put, get and acc are possible
    if( !toBeAttachedForGet_ )
    {
        GlobalArrayPut_ 	= 0;
        toBeAttachedForPut_ 	= false;
        GlobalArrayGet_ 	= &Z;
        toBeAttachedForGet_ 	= true;

        if( matrices_.empty() )
        {
            struct matrix_params_ mp;
            mp.data_.resize( p );
            mp.requests_.resize( p );
            mp.statuses_.resize( p );
            mp.base_ = NULL;
            // push back new matrix_params created
            // with default constructor
            matrices_.push_back( mp );
        }

        if( coords_.empty() )
        {
            struct coord_params_ cp;
            cp.coord_.resize( p );
            cp.requests_.resize( p );
            cp.statuses_.resize( p );
            cp.base_ = NULL;
            // push back new matrix_params created
            // with default constructor
            coords_.push_back( cp );
        }
    }

    mpi::Barrier( g.VCComm() );
}

// end-to-end blocking put/acc routines
template<typename T>
void AxpyInterface2<T>::Eput( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Eput" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Grid& g = Y.Grid();
    const Int XLDim = Z.LDim();
    const Int height = Z.Height();
    const Int width = Z.Width();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const T* XBuffer = Z.LockedBuffer();
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int YLDim = Y.LDim();
    Int matrix_index, coord_index;

    // data/coord send
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;

        if( numEntries > 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            // data
            const Int dindex =
                NextIndexData( destination,
                               numEntries,
                               Buffer,
                               &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][dindex].size() ) !=
                          numEntries ) LogicError( "Error in NextIndexData" ); )
                T* sendBuffer = matrices_[matrix_index].data_[destination][dindex].data();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            mpi::TaggedISend( sendBuffer, numEntries, destination,
                              DATA_PUT_TAG, g.VCComm(),
                              matrices_[matrix_index].requests_[destination][dindex] );
            // coordinates
            const Int cindex =
                NextIndexCoord( i, j,
                                destination,
                                Buffer,
                                &coord_index );
            Int* coord_ = reinterpret_cast<Int*>( coords_[coord_index].coord_[destination][cindex].data() );
            coord_[0] = i;
            coord_[1] = j;
            coord_[2] = numEntries;
            // post receive for coordinates
            mpi::TaggedISend( coord_, 3, destination,
                              COORD_PUT_TAG, g.VCComm(),
                              coords_[coord_index].requests_[destination][cindex] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }

    // poke
    Test( Z );
    // data/coord receive
    std::vector<T> recvVector_;

    for( Int step=0; step<p; ++step )
    {
        mpi::Status status;

        if( mpi::IProbe( mpi::ANY_SOURCE, DATA_PUT_TAG, g.VCComm(), status ) )
        {
            const Int source = status.MPI_SOURCE;
            // coordinates
            Int coord[3];
            mpi::TaggedRecv( coord, 3, source,
                             COORD_PUT_TAG, g.VCComm() );
            Int i = coord[0];
            Int j = coord[1];
            Int count = coord[2];
            recvVector_.resize( count );
            T* recvBuffer = recvVector_.data();
            // data
            mpi::TaggedRecv( recvBuffer, count, source,
                             DATA_PUT_TAG, g.VCComm() );
            // Update Y
            const T* XBuffer = recvBuffer;
            const Int colAlign = ( Y.ColAlign() + i ) % r;
            const Int rowAlign = ( Y.RowAlign() + j ) % c;
            const Int colShift = Shift( g.Row(), colAlign, r );
            const Int rowShift = Shift( g.Col(), rowAlign, c );
            const Int localHeight = Length( height, colShift, r );
            const Int localWidth = Length( width, rowShift, c );
            const Int iLocalOffset = Length( i, Y.ColShift(), r );
            const Int jLocalOffset = Length( j, Y.RowShift(), c );

            for( Int t = 0; t < localWidth; ++t )
            {
                T* YCol = Y.Buffer( iLocalOffset, jLocalOffset + t );
                const T* XCol = &XBuffer[t * localHeight];
                MemCopy( YCol, XCol, localHeight );
            }
        }
    }

    // wait
    Wait( Z );
    recvVector_.clear();
}

template<typename T>
void AxpyInterface2<T>::Eput( Matrix<T>& Z, Int i, Int j )
{ Eput( const_cast<const Matrix<T>&>( Z ), i, j ); }

// end to end blocking routines
template<typename T>
void AxpyInterface2<T>::Eacc( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Eacc" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Grid& g = Y.Grid();
    const Int XLDim = Z.LDim();
    const Int height = Z.Height();
    const Int width = Z.Width();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const T* XBuffer = Z.LockedBuffer();
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int YLDim = Y.LDim();
    // data/coord receive
    std::vector<T> recvVector_;
    Int matrix_index, coord_index;

    // data/coord send
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;

        if( numEntries > 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            // data
            const Int dindex =
                NextIndexData( destination,
                               numEntries,
                               Buffer,
                               &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][dindex].size() ) !=
                          numEntries ) LogicError( "Error in NextIndexData" ); )
                
	    T* sendBuffer = matrices_[matrix_index].data_[destination][dindex].data();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            mpi::TaggedISend( sendBuffer, numEntries, destination,
                              DATA_ACC_TAG, g.VCComm(),
                              matrices_[matrix_index].requests_[destination][dindex] );
            // coordinates
            const Int cindex =
                NextIndexCoord( i, j,
                                destination,
                                Buffer,
                                &coord_index );
            Int* coord_ = reinterpret_cast<Int*>( coords_[coord_index].coord_[destination][cindex].data() );
            coord_[0] = i;
            coord_[1] = j;
            coord_[2] = numEntries;
            // post receive for coordinates
            mpi::TaggedISend( coord_, 3, destination,
                              COORD_ACC_TAG, g.VCComm(),
                              coords_[coord_index].requests_[destination][cindex] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }

    // test for requests
    Test( Z );

    for( Int step=0; step<p; ++step )
    {
        mpi::Status status;

        if( mpi::IProbe( mpi::ANY_SOURCE, DATA_ACC_TAG, g.VCComm(), status ) )
        {
            const Int source = status.MPI_SOURCE;
            // coordinates
            Int coord[3];
            mpi::TaggedRecv( coord, 3, source,
                             COORD_ACC_TAG, g.VCComm() );
            Int i = coord[0];
            Int j = coord[1];
            Int count = coord[2];
            recvVector_.resize( count );
            T* recvBuffer = recvVector_.data();
            // data
            mpi::TaggedRecv( recvBuffer, count, source,
                             DATA_ACC_TAG, g.VCComm() );
            // Update Y
            const T* XBuffer = recvBuffer;
            const Int colAlign = ( Y.ColAlign() + i ) % r;
            const Int rowAlign = ( Y.RowAlign() + j ) % c;
            const Int colShift = Shift( g.Row(), colAlign, r );
            const Int rowShift = Shift( g.Col(), rowAlign, c );
            const Int localHeight = Length( height, colShift, r );
            const Int localWidth = Length( width, rowShift, c );
            const Int iLocalOffset = Length( i, Y.ColShift(), r );
            const Int jLocalOffset = Length( j, Y.RowShift(), c );

            for( Int t = 0; t < localWidth; ++t )
            {
                T* YCol = Y.Buffer( iLocalOffset, jLocalOffset + t );
                const T* XCol = &XBuffer[t * localHeight];

                for( Int s = 0; s < localHeight; ++s )
                    YCol[s] += XCol[s];
            }
        }
    }

    // wait for requests
    Wait( Z );
    recvVector_.clear();
}

template<typename T>
void AxpyInterface2<T>::Eacc( Matrix<T>& Z, Int i, Int j )
{ Eacc( const_cast<const Matrix<T>&>( Z ), i, j ); }

template<typename T>
void AxpyInterface2<T>::Get( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Get" ) )

    // a call to Attach with a non-const DistMatrix must set
    // toBeAttachedForGet_ also, if not then it is assumed that
    // the DistMatrix isn't attached
    if( !toBeAttachedForGet_ )
        LogicError( "Cannot perform this operation as matrix is not attached." );

    const DistMatrix<T>& X = *GlobalArrayGet_;
    const Int height = Z.Height();
    const Int width = Z.Width();

    if( i + height > X.Height() || j + width > X.Width() )
        LogicError( "Invalid submatrix for Iget" );

    T* XBuffer = Z.Buffer();
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );
    const Grid& g = X.Grid();
    const Int p = g.Size();
    const Int r = g.Height();
    const Int c = g.Width();
    Int coord_index;
    std::vector<T> recvVector_;

    // Send out the requests to all processes in the grid
    for( Int rank = 0; rank < p; ++rank )
    {
        const Int cindex =
            NextIndexCoord( i, j,
                            rank,
                            Buffer,
                            &coord_index );
        Int* coord = reinterpret_cast<Int*>( coords_[coord_index].coord_[rank][cindex].data() );
        coord[0] = i;
        coord[1] = j;
        coord[2] = -1;
        mpi::TaggedISend( coord, 3, rank,
                          REQUEST_GET_TAG, g.VCComm(),
                          coords_[coord_index].requests_[rank][cindex] );
    }

    // Receive all of the replies
    Int numReplies = 0;

    while( numReplies < p )
    {
        mpi::Status status;
        HandleGlobalToLocalData( Z );

        if( mpi::IProbe
            ( mpi::ANY_SOURCE, DATA_GET_TAG, g.VCComm(), status ) )
        {
            const Int source = status.MPI_SOURCE;
            // Ensure that we have a recv buffer
            const Int count = mpi::GetCount <T> ( status );
            recvVector_.resize( count );
            T* recvBuffer = recvVector_.data();
        
	    // Receive the data
            mpi::TaggedRecv
            ( recvBuffer, count, source, DATA_GET_TAG, g.VCComm() );
            
	    // Compute the local heights and offsets
            const Int myRow = g.Row();
            const Int myCol = g.Col();
            const Int colAlign = ( X.ColAlign() + i ) % r;
            const Int rowAlign = ( X.RowAlign() + j ) % c;
            const Int colShift = Shift( myRow, colAlign, r );
            const Int rowShift = Shift( myCol, rowAlign, c );
            const Int localHeight = Length( height, colShift, r );
            const Int localWidth = Length( width, rowShift, c );

            // Unpack the local matrix
            for( Int t = 0; t < localWidth; ++t )
            {
                //T *YCol = X.Buffer (0, rowShift + t * c);
                T* YCol = Z.Buffer( 0, rowShift + t * c );
                const T* XCol = &recvBuffer[t * localHeight];

                for( Int s = 0; s < localHeight; ++s )
                    YCol[colShift + s * r] = XCol[s];
            }

            ++numReplies;
            recvVector_.clear();
        }
    }
}

// nonblocking, no local completion
template<typename T>
void AxpyInterface2<T>::Iput( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Iput" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int XLDim = Z.LDim();
    Int matrix_index, coord_index;
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    const T* XBuffer = Z.LockedBuffer();
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );

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
            const Int dindex =
                NextIndexData( destination,
                               numEntries,
                               Buffer,
                               &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][dindex].size() ) !=
                          numEntries ) LogicError( "Error in NextIndexData" ); )
            T* sendBuffer = matrices_[matrix_index].data_[destination][dindex].data();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            // put request
            mpi::TaggedISend( sendBuffer, numEntries, destination,
                              DATA_PUT_TAG, g.VCComm(),
                              matrices_[matrix_index].requests_[destination][dindex] );
            // send coordinates
            const Int cindex =
                NextIndexCoord( i, j,
                                destination,
                                Buffer,
                                &coord_index );
            Int* coord_ = reinterpret_cast<Int*>( coords_[coord_index].coord_[destination][cindex].data() );
            coord_[0] = i;
            coord_[1] = j;
            coord_[2] = numEntries;
            // post receive for coordinates
            mpi::TaggedISend( coord_, 3, destination,
                              COORD_PUT_TAG, g.VCComm(),
                              coords_[coord_index].requests_[destination][cindex] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Iput( Matrix<T>& Z, Int i, Int j )
{ Iput( const_cast<const Matrix<T>&>( Z ), i, j ); }

template<typename T>
void AxpyInterface2<T>::Iget( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Iget" ) )

    // a call to Attach with a non-const DistMatrix must set
    // toBeAttachedForGet_ also, if not then it is assumed that
    // the DistMatrix isn't attached
    if( !toBeAttachedForGet_ )
        LogicError( "Cannot perform this operation as matrix is not attached." );

    const DistMatrix<T>& X = *GlobalArrayGet_;
    const Int height = Z.Height();
    const Int width = Z.Width();
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );
    Int coord_index;

    if( i + height > X.Height() || j + width > X.Width() )
        LogicError( "Invalid submatrix for Iget" );

    const Grid& g = X.Grid();
    const Int p = g.Size();

    // Send out the requests to all processes in the grid
    for( Int rank = 0; rank < p; ++rank )
    {
        // send coordinates
        const Int cindex =
            NextIndexCoord( i, j,
                            rank,
                            Buffer,
                            &coord_index );
        Int* coord_ = reinterpret_cast<Int*>( coords_[coord_index].coord_[rank][cindex].data() );
        coord_[0] = i;
        coord_[1] = j;
        coord_[2] = -1;
        // post receive for coordinates
        mpi::TaggedISend( coord_, 3, rank,
                          REQUEST_GET_TAG, g.VCComm(),
                          coords_[coord_index].requests_[rank][cindex] );
    }
}

// accumulate = Update Y(i:i+height-1,j:j+width-1) += X,
// where X is height x width
template<typename T>
void AxpyInterface2<T>::Iacc( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Iacc" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;
    Int matrix_index, coord_index;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int XLDim = Z.LDim();
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    const T* XBuffer = Z.LockedBuffer();
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );

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
            const Int dindex =
                NextIndexData( destination,
                               numEntries,
                               Buffer,
                               &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][dindex].size() ) !=
                          numEntries ) LogicError( "Error in NextIndexData" ); )
            T* sendBuffer = matrices_[matrix_index].data_[destination][dindex].data();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            // acc request
            mpi::TaggedISend( sendBuffer, numEntries, destination,
                              DATA_ACC_TAG, g.VCComm(),
                              matrices_[matrix_index].requests_[destination][dindex] );
            // send coordinates
            const Int cindex =
                NextIndexCoord( i, j,
                                destination,
                                Buffer,
                                &coord_index );
            Int* coord_ = reinterpret_cast<Int*>( coords_[coord_index].coord_[destination][cindex].data() );
            coord_[0] = i;
            coord_[1] = j;
            coord_[2] = numEntries;
            mpi::TaggedISend( coord_, 3, destination,
                              COORD_ACC_TAG, g.VCComm(),
                              coords_[coord_index].requests_[destination][cindex] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Iacc( Matrix<T>& Z, Int i, Int j )
{ Iacc( const_cast<const Matrix<T>&>( Z ), i, j ); }

// nonblocking, local completion
template<typename T>
void AxpyInterface2<T>::Put( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Put" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int XLDim = Z.LDim();
    Int matrix_index, coord_index;
    
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();

    // copy local matrix buffer
    const Int my_rank = g.VCRank();
    const Int numCreated = dataVectors_[my_rank].size();
    dataVectors_[my_rank].resize( numCreated + 1 );
    dataVectors_[my_rank][numCreated].resize( width * height );
    
    const void* Buffer = static_cast <void*>( const_cast <T*>( Z.LockedBuffer() ) );
    T* ZBuffer = reinterpret_cast <T*>( dataVectors_[my_rank][numCreated].data() );
    MemCopy( ZBuffer, reinterpret_cast <const T*>( Buffer ),
             height * width );
    T* XBuffer = reinterpret_cast <T*>( ZBuffer );

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
            // data
            const Int dindex =
                NextIndexData( destination,
                               numEntries,
                               Buffer,
                               &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][dindex].size() ) !=
                          numEntries ) LogicError( "Error in NextIndexData" ); )
            T* sendBuffer = matrices_[matrix_index].data_[destination][dindex].data();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            mpi::TaggedISend( sendBuffer, numEntries, destination,
                              DATA_PUT_TAG, g.VCComm(),
                              matrices_[matrix_index].requests_[destination][dindex] );
            // send coordinates
            const Int cindex =
                NextIndexCoord( i, j,
                                destination,
                                Buffer,
                                &coord_index );
            Int* coord_ = reinterpret_cast<Int*>( coords_[coord_index].coord_[destination][cindex].data() );
            coord_[0] = i;
            coord_[1] = j;
            coord_[2] = numEntries;
            mpi::TaggedISend( coord_, 3, destination,
                              COORD_PUT_TAG, g.VCComm(),
                              coords_[coord_index].requests_[destination][cindex] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Put( Matrix<T>& Z, Int i, Int j )
{ Put( const_cast<const Matrix<T>&>( Z ), i, j ); }

// input buffer could be modified upon exit
// from this function
template<typename T>
void AxpyInterface2<T>::Acc( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Acc" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;
    Int matrix_index, coord_index;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int XLDim = Z.LDim();

    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    
    // copy local matrix buffer
    const Int my_rank = g.VCRank();
    const Int numCreated = dataVectors_[my_rank].size();
    dataVectors_[my_rank].resize( numCreated + 1 );
    dataVectors_[my_rank][numCreated].resize( width * height );
    
    const void* Buffer = static_cast <void*>( const_cast <T*>( Z.LockedBuffer() ) );
    T* ZBuffer = reinterpret_cast <T*>( dataVectors_[my_rank][numCreated].data() );
    MemCopy( ZBuffer, reinterpret_cast <const T*>( Buffer ),
             height * width );
    T* XBuffer = reinterpret_cast <T*>( ZBuffer );

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
            // data
            const Int dindex =
                NextIndexData( destination,
                               numEntries,
                               Buffer,
                               &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][dindex].size() ) !=
                          numEntries ) LogicError( "Error in NextIndexData" ); )
            T* sendBuffer = matrices_[matrix_index].data_[destination][dindex].data();

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            // acc request
            mpi::TaggedISend( sendBuffer, numEntries, destination,
                              DATA_ACC_TAG, g.VCComm(),
                              matrices_[matrix_index].requests_[destination][dindex] );
            // send coordinates
            const Int cindex =
                NextIndexCoord( i, j,
                                destination,
                                Buffer,
                                &coord_index );
            Int* coord_ = reinterpret_cast<Int*>( coords_[coord_index].coord_[destination][cindex].data() );
            coord_[0] = i;
            coord_[1] = j;
            coord_[2] = numEntries;
            mpi::TaggedISend( coord_, 3, destination,
                              COORD_ACC_TAG, g.VCComm(),
                              coords_[coord_index].requests_[destination][cindex] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Acc( Matrix<T>& Z, Int i, Int j )
{ Acc( const_cast<const Matrix<T>&>( Z ), i, j ); }

// waitany implementation
// cannot use mpi::Waitany
// as of now because request
// objects are vector of deques
template<typename T>
void AxpyInterface2<T>::WaitAny( const Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::WaitAny" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex, coordIndex;
    const Int numMatrices = matrices_.size();
    const Int numCoords = coords_.size();
    const void* base_address = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );

    // search for matrix base
    for( Int m = 0; m < numMatrices; m++ )
    {
        if( matrices_[m].base_ == base_address )
        {
            matrixIndex = m;
            break;
        }

        matrixIndex = m+1;
    }

    // search for matrix base in coords
    for( Int c = 0; c < numCoords; c++ )
    {
        if( coords_[c].base_ == base_address )
        {
            coordIndex = c;
            break;
        }

        coordIndex = c+1;
    }

    // matrix not found
    if( matrixIndex == numMatrices &&
        coordIndex == numCoords )
        return;

    // data
    for( int rank = 0; rank < p; ++rank )
    {
        if( matrices_[matrixIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numDataStatuses = matrices_[matrixIndex].requests_[rank].size();

        for( int i = 0; i < numDataStatuses; i++ )
        {
            if( !matrices_[matrixIndex].statuses_[rank][i] )
            {
                mpi::Wait( matrices_[matrixIndex].requests_[rank][i] );
                matrices_[matrixIndex].statuses_[rank][i] = true;
                return;
            }
        }
    }

    // coordinates
    for( int rank = 0; rank < p; ++rank )
    {
        if( coords_[coordIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numCoordStatuses = coords_[coordIndex].requests_[rank].size();

        for( int i = 0; i < numCoordStatuses; i++ )
        {
            if( !coords_[coordIndex].statuses_[rank][i] )
            {
                mpi::Wait( coords_[coordIndex].requests_[rank][i] );
                coords_[coordIndex].statuses_[rank][i] = true;
                return;
            }
        }
    }
}

template<typename T>
void AxpyInterface2<T>::WaitAny( Matrix<T>& Z )
{ WaitAny( const_cast<const Matrix<T>&>( Z ) ); }

template<typename T>
void AxpyInterface2<T>::Wait( const Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Wait" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex, coordIndex;
    const Int numMatrices = matrices_.size();
    const Int numCoords = coords_.size();
    const void* base_address = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );

    // search for matrix base
    for( Int m = 0; m < numMatrices; m++ )
    {
        if( matrices_[m].base_ == base_address )
        {
            matrixIndex = m;
            break;
        }

        matrixIndex = m+1;
    }

    // search for matrix base in coords
    for( Int c = 0; c < numCoords; c++ )
    {
        if( coords_[c].base_ == base_address )
        {
            coordIndex = c;
            break;
        }

        coordIndex = c+1;
    }

    // matrix not found
    if( matrixIndex == numMatrices &&
        coordIndex == numCoords )
        return;

    // data
    for( int rank = 0; rank < p; ++rank )
    {
        if( matrices_[matrixIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numDataStatuses = matrices_[matrixIndex].requests_[rank].size();

        for( int i = 0; i < numDataStatuses; i++ )
        {
            mpi::Wait( matrices_[matrixIndex].requests_[rank][i] );
            matrices_[matrixIndex].statuses_[rank][i] = true;
        }
    }

    // coordinates
    for( int rank = 0; rank < p; ++rank )
    {
        if( coords_[coordIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numCoordStatuses = coords_[coordIndex].requests_[rank].size();

        for( int i = 0; i < numCoordStatuses; i++ )
        {
            mpi::Wait( coords_[coordIndex].requests_[rank][i] );
            coords_[coordIndex].statuses_[rank][i] = true;
        }
    }
}

template<typename T>
void AxpyInterface2<T>::Wait( Matrix<T>& Z )
{ Wait( const_cast<const Matrix<T>&>( Z ) ); }

template<typename T>
void AxpyInterface2<T>::Waitall()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Waitall" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex, coordIndex;
    const Int numMatrices = matrices_.size();
    const Int numCoords = coords_.size();

    // data
    for( int matrixIndex = 0; matrixIndex < numMatrices; ++matrixIndex )
    {
        for( int rank = 0; rank < p; ++rank )
        {
            const Int numDataStatuses = matrices_[matrixIndex].requests_[rank].size();

            for( int i = 0; i < numDataStatuses; i++ )
            {
                mpi::Wait( matrices_[matrixIndex].requests_[rank][i] );
                matrices_[matrixIndex].statuses_[rank][i] = true;
            }
        }
    }

    // coordinates
    for( int coordIndex = 0; coordIndex < numCoords; ++coordIndex )
    {
        for( int rank = 0; rank < p; ++rank )
        {
            const Int numCoordStatuses = coords_[coordIndex].requests_[rank].size();

            for( int i = 0; i < numCoordStatuses; i++ )
            {
                mpi::Wait( coords_[coordIndex].requests_[rank][i] );
                coords_[coordIndex].statuses_[rank][i] = true;
            }
        }
    }
}

template<typename T>
bool AxpyInterface2<T>::Test( const Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Test" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex, coordIndex;
    const Int numMatrices = matrices_.size();
    const Int numCoords = coords_.size();
    const void* base_address = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );

    // search for matrix base
    for( Int m = 0; m < numMatrices; m++ )
    {
        if( matrices_[m].base_ == base_address )
        {
            matrixIndex = m;
            break;
        }

        matrixIndex = m+1;
    }

    // search for matrix base in coords
    for( Int c = 0; c < numCoords; c++ )
    {
        if( coords_[c].base_ == base_address )
        {
            coordIndex = c;
            break;
        }

        coordIndex = c+1;
    }

    // matrix not found
    if( matrixIndex == numMatrices &&
        coordIndex == numCoords )
        return true;

    for( int rank = 0; rank < p; ++rank )
    {
        if( matrices_[matrixIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numDataStatuses = matrices_[matrixIndex].requests_[rank].size();

        for( int i = 0; i < numDataStatuses; i++ )
        {
            matrices_[matrixIndex].statuses_[rank][i] =
                !mpi::Test( matrices_[matrixIndex].requests_[rank][i] );

            if( matrices_[matrixIndex].statuses_[rank][i] )
                return false;
        }
    }

    for( int rank = 0; rank < p; ++rank )
    {
        if( coords_[coordIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numCoordStatuses = coords_[coordIndex].requests_[rank].size();

        for( int i = 0; i < numCoordStatuses; i++ )
        {
            coords_[coordIndex].statuses_[rank][i] =
                !mpi::Test( coords_[coordIndex].requests_[rank][i] );

            if( coords_[coordIndex].statuses_[rank][i] )
                return false;
        }
    }

    return true;
}

template<typename T>
bool AxpyInterface2<T>::Test( Matrix<T>& Z )
{ return Test( const_cast<const Matrix<T>&>( Z ) ); }

template<typename T>
bool AxpyInterface2<T>::TestAny( const Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::TestAny" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex, coordIndex;
    const Int numMatrices = matrices_.size();
    const Int numCoords = coords_.size();
    const void* base_address = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );

    // search for matrix base
    for( Int m = 0; m < numMatrices; m++ )
    {
        if( matrices_[m].base_ == base_address )
        {
            matrixIndex = m;
            break;
        }

        matrixIndex = m+1;
    }

    // search for matrix base in coords
    for( Int c = 0; c < numCoords; c++ )
    {
        if( coords_[c].base_ == base_address )
        {
            coordIndex = c;
            break;
        }

        coordIndex = c+1;
    }

    // matrix not found
    if( matrixIndex == numMatrices &&
        coordIndex == numCoords )
        return true;

    for( int rank = 0; rank < p; ++rank )
    {
        if( matrices_[matrixIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numDataStatuses = matrices_[matrixIndex].requests_[rank].size();

        for( int i = 0; i < numDataStatuses; i++ )
        {
            matrices_[matrixIndex].statuses_[rank][i] =
                !mpi::Test( matrices_[matrixIndex].requests_[rank][i] );

            if( matrices_[matrixIndex].statuses_[rank][i] )
                continue;
            else
                return true;
        }
    }

    for( int rank = 0; rank < p; ++rank )
    {
        if( coords_[coordIndex].statuses_[rank].size() == 0 )
            continue;

        const Int numCoordStatuses = coords_[coordIndex].requests_[rank].size();

        for( int i = 0; i < numCoordStatuses; i++ )
        {
            coords_[coordIndex].statuses_[rank][i] =
                !mpi::Test( coords_[coordIndex].requests_[rank][i] );

            if( coords_[coordIndex].statuses_[rank][i] )
                continue;
            else
                return true;
        }
    }

    return false;
}

template<typename T>
bool AxpyInterface2<T>::TestAny( Matrix<T>& Z )
{ return TestAny( const_cast<const Matrix<T>&>( Z ) ); }

template<typename T>
bool AxpyInterface2<T>::Testall()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Testall" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    const Int numMatrices = matrices_.size();
    const Int numCoords = coords_.size();

    // data
    for( int matrixIndex = 0; matrixIndex < numMatrices; ++matrixIndex )
    {
        for( int rank = 0; rank < p; ++rank )
        {
            if( matrices_[matrixIndex].statuses_[rank].size() == 0 )
                continue;

            const Int numDataStatuses = matrices_[matrixIndex].requests_[rank].size();

            for( int i = 0; i < numDataStatuses; i++ )
            {
                matrices_[matrixIndex].statuses_[rank][i] =
                    !mpi::Test( matrices_[matrixIndex].requests_[rank][i] );

                if( matrices_[matrixIndex].statuses_[rank][i] )
                    return false;
            }
        }
    }

    // coordinates
    for( int coordIndex = 0; coordIndex < numCoords; ++coordIndex )
    {
        for( int rank = 0; rank < p; ++rank )
        {
            if( coords_[coordIndex].statuses_[rank].size() == 0 )
                continue;

            const Int numCoordStatuses = coords_[coordIndex].requests_[rank].size();

            for( int i = 0; i < numCoordStatuses; i++ )
            {
                coords_[coordIndex].statuses_[rank][i] =
                    !mpi::Test( coords_[coordIndex].requests_[rank][i] );

                if( coords_[coordIndex].statuses_[rank][i] )
                    return false;
            }
        }
    }

    return true;
}

// This is non-collective flush
// This will ensure local+remote completion
// if Z is  const then only Put/Acc is possible
template<typename T>
void AxpyInterface2<T>::Flush( const Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Flush" ) )

    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer before flushing." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    bool DONE = false;
    mpi::Status status;

    while( !DONE )
    {
        if( mpi::IProbe( mpi::ANY_SOURCE, mpi::ANY_TAG, g.VCComm(), status ) )
        {
            switch( status.MPI_TAG )
            {
		case DATA_PUT_TAG:
		    {
			HandleLocalToGlobalData( Z, status.MPI_SOURCE );
			break;
		    }
		case DATA_ACC_TAG:
		    {
			HandleLocalToGlobalAcc( Z, status.MPI_SOURCE );
			break;
		    }
	    }
        }

        // wait for requests to
        // complete one by one
        WaitAny( Z );
        DONE = Test( Z );
    }
}

template<typename T>
void AxpyInterface2<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Flush" ) )

    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer before flushing." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    bool DONE = false;
    mpi::Status status;

    while( !DONE )
    {
        if( mpi::IProbe( mpi::ANY_SOURCE, mpi::ANY_TAG, g.VCComm(), status ) )
        {
            switch( status.MPI_TAG )
            {
		case DATA_PUT_TAG:
		    {
			HandleLocalToGlobalData( Z, status.MPI_SOURCE );
			break;
		    }
		case DATA_ACC_TAG:
		    {
			HandleLocalToGlobalAcc( Z, status.MPI_SOURCE );
			break;
		    }
		case REQUEST_GET_TAG:
		    {
			HandleGlobalToLocalData( Z );
			break;
		    }
            }
        }

        // wait for requests to
        // complete one by one
        WaitAny( Z );
        DONE = Test( Z );
    }
}

template <typename T>
void AxpyInterface2<T>::HandleLocalToGlobalData( const Matrix<T>& Z, Int source )
{
    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRow = g.Row();
    const Int myCol = g.Col();
    int height = Z.Height();
    int width = Z.Width();
    
    // post receive for coordinates
    Int coord[3];
    mpi::TaggedRecv( coord, 3, source,
                     COORD_PUT_TAG, g.VCComm() );
    Int i = coord[0];
    Int j = coord[1];
    Int count = coord[2];
    // data vector
    std::vector<T> getVector_;
    getVector_.resize( count );

    DEBUG_ONLY( if( count < Int( sizeof( T ) ) )
                LogicError( "Count was too small" ); )
    DEBUG_ONLY( if( Int( getVector_.size() ) != count )
                LogicError( "Not enough space allocated" ); )
            
    // post receive for data
    T* getBuffer = getVector_.data();
    mpi::TaggedRecv( getBuffer, count, source,
                     DATA_PUT_TAG, g.VCComm() );
    
    // Update Y
    const T* XBuffer = const_cast <const T*>( getBuffer );
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int colShift = Shift( myRow, colAlign, r );
    const Int rowShift = Shift( myCol, rowAlign, c );
    const Int localHeight = Length( height, colShift, r );
    const Int localWidth = Length( width, rowShift, c );
    const Int iLocalOffset = Length( i, Y.ColShift(), r );
    const Int jLocalOffset = Length( j, Y.RowShift(), c );

    for( Int t = 0; t < localWidth; ++t )
    {
        T* YCol = Y.Buffer( iLocalOffset, jLocalOffset + t );
        const T* XCol = &XBuffer[t * localHeight];
        MemCopy( YCol, XCol, localHeight );
    }

    // Free the memory
    getVector_.clear();
}

template<typename T>
void AxpyInterface2<T>::HandleLocalToGlobalData( Matrix<T>& Z, Int source )
{ HandleLocalToGlobalData( const_cast<const Matrix<T>&>( Z ), source ); }

// replica of above function except this accumulates
template <typename T>
void AxpyInterface2<T>::HandleLocalToGlobalAcc( const Matrix<T>& Z, Int source )
{
    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRow = g.Row();
    const Int myCol = g.Col();
    const int height = Z.Height();
    const int width = Z.Width();

    // post receive for coordinates
    Int coord[3];
    mpi::TaggedRecv( coord, 3, source,
                     COORD_ACC_TAG, g.VCComm() );
    Int i = coord[0];
    Int j = coord[1];
    Int count = coord[2];
    // data buffer
    std::vector<T> getVector_;
    getVector_.resize( count );

    DEBUG_ONLY( if( count < Int( sizeof( T ) ) )
                LogicError( "Count was too small" ); )
    DEBUG_ONLY( if( Int( getVector_.size() ) != count )
                LogicError( "Not enough space allocated" ); )
            
    // post receive for data
    T* getBuffer = getVector_.data();
    mpi::TaggedRecv( getBuffer, count, source,
                     DATA_ACC_TAG, g.VCComm() );
    // Update Y
    const T* XBuffer = const_cast <const T*>( getBuffer );
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int colShift = Shift( myRow, colAlign, r );
    const Int rowShift = Shift( myCol, rowAlign, c );
    const Int localHeight = Length( height, colShift, r );
    const Int localWidth = Length( width, rowShift, c );
    const Int iLocalOffset = Length( i, Y.ColShift(), r );
    const Int jLocalOffset = Length( j, Y.RowShift(), c );

    for( Int t = 0; t < localWidth; ++t )
    {
        T* YCol = Y.Buffer( iLocalOffset, jLocalOffset + t );
        const T* XCol = &XBuffer[t * localHeight];

        for( Int s = 0; s < localHeight; ++s )
            YCol[s] += XCol[s];
    }

    // Free the memory
    getVector_.clear();
}

template<typename T>
void AxpyInterface2<T>::HandleLocalToGlobalAcc( Matrix<T>& Z, Int source )
{ HandleLocalToGlobalAcc( const_cast<const Matrix<T>&>( Z ), source ); }

// handle request for data, post a matching isend
template <typename T>
void AxpyInterface2<T>::HandleGlobalToLocalData( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface::HandleGlobalToLocalData" ) )

    if( !toBeAttachedForGet_ )
        LogicError( "Local matrix cannot be updated" );

    const DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myRow = g.Row();
    const Int myCol = g.Col();
    Int i, j;
    Int matrix_index;
    std::vector<T> recvVector_;
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );
    const Int XLDim = Z.LDim();
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    for( Int step = 0; step < p; step++ )
    {
        mpi::Status status;

        if( mpi::IProbe( mpi::ANY_SOURCE, REQUEST_GET_TAG, g.VCComm(), status ) )
        {
            const Int source = status.MPI_SOURCE;
            // post receive for coordinates
            Int coord[3];
            mpi::TaggedRecv( coord, 3, source,
                             REQUEST_GET_TAG, g.VCComm() );
            i = coord[0];
            j = coord[1];
            // we need the localwidth/height here,
            // used also to calculate numEntries
            const Int colAlign = ( Y.ColAlign() + i ) % r;
            const Int rowAlign = ( Y.RowAlign() + j ) % c;
            const Int colShift = Shift( myRow, colAlign, r );
            const Int rowShift = Shift( myCol, rowAlign, c );
            const Int localHeight = Length( height, colShift, r );
            const Int localWidth = Length( width, rowShift, c );
            const Int iLocalOffset = Length( i, Y.ColShift(), r );
            const Int jLocalOffset = Length( j, Y.RowShift(), c );
            const Int numEntries = localHeight * localWidth;

            DEBUG_ONLY( if( numEntries < Int( sizeof( T ) ) )
                        LogicError( "Count was too small" ); )
                const Int index =
                    NextIndexData( source,
                                   numEntries,
                                   Buffer,
                                   &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[source][index].size() ) !=
                          numEntries ) LogicError( "Error in NextIndexData" ); )
                T* replyBuffer = matrices_[matrix_index].data_[source][index].data();

            for( Int t = 0; t < localWidth; ++t )
            {
                T* sendCol = &replyBuffer[t * localHeight];
                const T* XCol = Y.LockedBuffer( iLocalOffset, jLocalOffset + t );
                MemCopy( sendCol, XCol, localHeight );
            }

            // Fire off non-blocking send
            mpi::TaggedISend( replyBuffer, numEntries, source,
                              DATA_GET_TAG, g.VCComm(),
                              matrices_[matrix_index].requests_[source][index] );
        }

        // receive data
        if( mpi::IProbe
            ( mpi::ANY_SOURCE, DATA_GET_TAG, g.VCComm(), status ) )
        {
            const Int source = status.MPI_SOURCE;
            // Ensure that we have a recv buffer
            const Int count = mpi::GetCount <T> ( status );
            recvVector_.resize( count );
            T* recvBuffer = recvVector_.data();
            // Receive the data
            mpi::TaggedRecv
            ( recvBuffer, count, source, DATA_GET_TAG, g.VCComm() );
            // Compute the local heights and offsets
            const Int myRow = g.Row();
            const Int myCol = g.Col();
            const Int colAlign = ( Y.ColAlign() + i ) % r;
            const Int rowAlign = ( Y.RowAlign() + j ) % c;
            const Int colShift = Shift( myRow, colAlign, r );
            const Int rowShift = Shift( myCol, rowAlign, c );
            const Int localHeight = Length( height, colShift, r );
            const Int localWidth = Length( width, rowShift, c );

            // Unpack the local matrix
            for( Int t = 0; t < localWidth; ++t )
            {
                T* YCol = Z.Buffer( 0, rowShift + t * c );
                const T* XCol = &recvBuffer[t * localHeight];

                for( Int s = 0; s < localHeight; ++s )
                    YCol[colShift + s * r] = XCol[s];
            }
        }
    }

    recvVector_.clear();
}

// detach collectively
template<typename T>
void AxpyInterface2<T>::Detach()
{
    DEBUG_ONLY( CallStackEntry cse( "AxpyInterface2::Detach" ) )

    // destructor will call detach again...
    if( detached_ )
        return;

    if( !attached_ )
        LogicError( "Must attach before detaching." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();

    mpi::Barrier( g.VCComm() );
    
    attached_ 		= false;
    detached_ 		= true;
    
    toBeAttachedForPut_ = false;
    toBeAttachedForGet_ = false;
    
    GlobalArrayPut_ 	= 0;
    GlobalArrayGet_ 	= 0;

    if( !dataVectors_.empty() )
        dataVectors_.clear();

    matrices_.clear();
    coords_.clear();
}

#define PROTO(T) template class AxpyInterface2<T>;
#include "El/macros/Instantiate.h"

} // namespace El
