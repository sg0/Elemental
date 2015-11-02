/*
   Copyright (c) 2009-2014, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   Copyright (c) 2014, Jeff Hammond (Intel)
   Copyright (c) 2014, Sayan Ghosh (Washington State University)
   All rights reserved.

Authors:
Jeff Hammond adapted the RMA interface from the AXPY one.

This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include <cassert>

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY)  
namespace El
{
// constructor
template<typename T>
RmaInterface<T>::RmaInterface()
    : GlobalArrayPut_( 0 ), GlobalArrayGet_( 0 ),
      matrices_( 0 ), 
      window( mpi::WIN_NULL ),
      putVector_( 0 ), getVector_( 0 ),
      toBeAttachedForPut_( false ), toBeAttachedForGet_( false ),
      attached_( false ), detached_( true )
{}

template<typename T>
RmaInterface<T>::RmaInterface( DistMatrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::RmaInterface" ) )
    attached_ 			= false;
    detached_ 			= true;
    toBeAttachedForGet_ 	= false;
    toBeAttachedForPut_ 	= false;
    GlobalArrayPut_ 		= 0;
    GlobalArrayGet_ 		= 0;
    window 		        = mpi::WIN_NULL;
}

// until attach, I am not setting anything
// which might not be a good thing to do,
// but would modify this eventually
template<typename T>
RmaInterface<T>::RmaInterface( const DistMatrix<T>& X )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::RmaInterface" ) )
    attached_ 			= false;
    detached_ 			= true;
    toBeAttachedForGet_ 	= false;
    toBeAttachedForPut_ 	= false;
    GlobalArrayPut_ 		= 0;
    GlobalArrayGet_ 		= 0;
    window 	    		= mpi::WIN_NULL;
}

template<typename T>
RmaInterface<T>::~RmaInterface()
{
    {
        if( std::uncaught_exception() )
        {
            std::ostringstream os;
            os << "Uncaught exception detected during RmaInterface destructor "
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

template<typename T>
void RmaInterface<T>::Attach( DistMatrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Attach" ) )

    // attached_ will be only set in Attach
    // and only unset in Detach
    if( !attached_ && detached_ )
    {
        attached_ = true;
        detached_ = false;
    }
    else
        LogicError( "Must detach before reattaching." );

    // if DistMatrix is non-const, all one-sided
    // transfers -- put, get and acc are possible
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
    {
        GlobalArrayPut_ 		= &Z;
        toBeAttachedForPut_ 		= true;
        GlobalArrayGet_ 		= &Z;
        toBeAttachedForGet_ 		= true;
        const Grid& g = Z.Grid();
        const Int p = g.Size();

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

        if( putVector_.empty() )
        {
            getVector_.resize( p );
            putVector_.resize( p );
        }
#if defined(EL_USE_WIN_CREATE_FOR_RMA) && \
	!defined(EL_USE_WIN_ALLOC_FOR_RMA)
        // TODO rma related checks
        // creation of window
        const Int numEntries = Z.LocalHeight() * Z.LocalWidth();
        const Int bufferSize = numEntries * sizeof( T );
        void * baseptr = reinterpret_cast<void *>( Z.Buffer() );
        mpi::WindowCreate( baseptr, bufferSize, g.VCComm(), window );
        mpi::WindowLock( window );
#endif
#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
        const Int numEntries = Z.LocalHeight() * Z.LocalWidth();
        const Int bufferSize = numEntries * sizeof( T );
	
	mpi::WindowAllocate( bufferSize, g.VCComm(), window );
       
	T * baseptr = reinterpret_cast< T *>( mpi::GetWindowBase( window ) );
	Z.SetWindowBase( baseptr );

	mpi::WindowLock( window );
#endif
	mpi::Barrier( g.VCComm() );
    }
}

// for gets
template<typename T>
void RmaInterface<T>::Attach( const DistMatrix<T>& X )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Attach" ) )

    if( !attached_ && detached_ )
    {
        attached_ = true;
        detached_ = false;
    }
    else
        LogicError( "Must detach before reattaching." );

    if( !toBeAttachedForGet_ )
    {
        GlobalArrayGet_ 	= &X;
        toBeAttachedForGet_ 	= true;
        GlobalArrayPut_ 	= 0;
        toBeAttachedForPut_ 	= false;
        const Grid& g = X.Grid();
        const Int p = g.Size();

        if( getVector_.size() != p )
            getVector_.resize( p );

#if defined(EL_USE_WIN_CREATE_FOR_RMA) &&\
	!defined(EL_USE_WIN_ALLOC_FOR_RMA)
        //TODO rma related checks
        const Int numEntries = X.LocalHeight() * X.LocalWidth();
        const Int bufferSize = numEntries * sizeof( T );
        void * baseptr = static_cast<void *>( const_cast<T *>( X.LockedBuffer() ) );
        mpi::WindowCreate( baseptr, bufferSize, g.VCComm(), window );
        mpi::WindowLock( window );
#endif
#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)       
	LogicError ("Const DistMatrix cannot be modified, select EL_USE_WIN_CREATE_FOR_RMA");
#endif
	mpi::Barrier( g.VCComm() );
    }
}

// for standard passive rma
template<typename T>
Int RmaInterface<T>::NextIndex
( Int dataSize,
  std::deque <std::vector<T>>& dataVectors )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::NextIndex" ) )
    
    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "No DistMatrix attached" );
    if (dataSize < 0)
	LogicError ("Resize factor cannot be negative");
    
    const Int numCreated = dataVectors.size();
    dataVectors.resize( numCreated + 1 );
    dataVectors[numCreated].resize( dataSize );

    return numCreated;
}

// for request-based passive rma
template<typename T>
Int RmaInterface<T>::NextIndex(
    Int target,
    Int dataSize,
    const void* base_address,
    Int* mindex )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::NextIndex" ) )
    if ( base_address == NULL )
        LogicError( "Base address is NULL" );
    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "No DistMatrix attached" );
    if (dataSize < 0)
	LogicError ("Resize factor cannot be negative");

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

    DEBUG_ONLY( if( numCreated != Int( matrices_[matrixIndex].requests_[target].size() ) ||
                    numCreated != Int( matrices_[matrixIndex].statuses_[target].size() ) )
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

// request based RMA operations
template<typename T>
void RmaInterface<T>::Rput( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Rput" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do rma related checks
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
    const Int dm_height = Y.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    const T* XBuffer = Z.LockedBuffer();
    const void* Buffer = static_cast<void*>( const_cast<T*>( Z.LockedBuffer() ) );
    Int matrix_index;

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
            const Int index =
                NextIndex( destination,
                           numEntries,
                           Buffer,
                           &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][index].size() ) !=
                          numEntries ) LogicError( "Error in NextIndex" ); )
            T* sendBuffer = reinterpret_cast<T*>( matrices_[matrix_index].data_[destination][index].data() );
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );
	    Int iMapped, jMapped;
	    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
	    {
		iMapped = i; 
		jMapped = ( j / c );
	    }
	    else // 2-D process grid 
	    {
		iMapped = ( i / r ); 
		jMapped = ( j / c );
	    }

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }    
	    // put
	    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
	    mpi::Rput( sendBuffer, numEntries, destination, 
		    disp, numEntries, window, matrices_[matrix_index].requests_[destination][index] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void RmaInterface<T>::Rput( Matrix<T>& Z, Int i, Int j )
{
    Rput( const_cast<const Matrix<T>&>( Z ), i, j );
}

// accumulate = Update Y(i:i+height-1,j:j+width-1) += X,
// where X is height x width
template<typename T>
void RmaInterface<T>::Racc( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Racc" ) )

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated." );

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative." );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix." );

    //do rma related checks
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int XLDim = Z.LDim();
    const Int YLDim = Y.LDim();
    // local matrix width and height
    const Int dm_height = Y.Height();
    const Int height = Z.Height();
    const Int width = Z.Width();
    const T* XBuffer = Z.LockedBuffer();
    const void* Buffer = static_cast <void*>( const_cast <T*>( Z.LockedBuffer() ) );
    Int matrix_index;
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

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
            const Int index =
                NextIndex( destination,
                           numEntries,
                           Buffer,
                           &matrix_index );

            DEBUG_ONLY( if
                        ( Int( matrices_[matrix_index].data_[destination][index].size() ) !=
                          numEntries ) LogicError( "Error in NextIndex" ); )
            
	    T* sendBuffer = reinterpret_cast<T*>( matrices_[matrix_index].data_[destination][index].data() );
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );

	    Int iMapped, jMapped;
	    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
	    {
		iMapped = i; 
		jMapped = ( j / c );
	    }
	    else // 2-D process grid 
	    {
		iMapped = ( i / r ); 
		jMapped = ( j / c );
	    }

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }
	    // acc
	    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
	    mpi::Racc( sendBuffer, numEntries, destination, 
		    disp, numEntries, window, matrices_[matrix_index].requests_[destination][index] );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void RmaInterface<T>::Racc( Matrix<T>& Z, Int i, Int j )
{
    Racc( const_cast<const Matrix<T>&>( Z ), i, j );
}

// Locally Blocking
template<typename T>
void RmaInterface<T>::Put( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Put" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do rma related checks
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
    const Int dm_height = Y.Height();
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    const T* XBuffer = Z.LockedBuffer();

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
            const Int index =
                NextIndex( numEntries,
                           putVector_[destination] );
            T* sendBuffer = putVector_[destination][index].data();
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );

	    Int iMapped, jMapped;
	    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
	    {
		iMapped = i; 
		jMapped = ( j / c );
	    }
	    else // 2-D process grid 
	    {
		iMapped = ( i / r ); 
		jMapped = ( j / c );
	    }

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }
	    
	    // put
	    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
	    mpi::Iput( sendBuffer, numEntries, destination, disp, numEntries, window );
            mpi::FlushLocal( destination, window );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void RmaInterface<T>::Put( Matrix<T>& Z, Int i, Int j )
{
    Put( const_cast<const Matrix<T>&>( Z ), i, j );
}

template<typename T>
void RmaInterface<T>::Acc( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Acc" ) )
    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do rma related checks
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
    const Int ZLDim = Z.LDim();
    // local matrix width and height
    const Int dm_height = Y.Height();
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    const T* XBuffer = Z.LockedBuffer();

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
            const Int index =
                NextIndex( numEntries,
                           putVector_[destination] );
            T* sendBuffer = putVector_[destination][index].data();
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );

	    Int iMapped, jMapped;
	    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
	    {
		iMapped = i; 
		jMapped = ( j / c );
	    }
	    else // 2-D process grid 
	    {
		iMapped = ( i / r ); 
		jMapped = ( j / c );
	    }
	    
            for( Int t = 0; t < localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t * localHeight];
                const T* thisXCol = &XBuffer[( rowShift + t * c ) * ZLDim];

                for( Int s = 0; s < localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift + s * r];
            }
	    
	    // acc
	    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
	    mpi::Iacc( sendBuffer, numEntries, destination, disp, numEntries, window );
	    mpi::FlushLocal( destination, window );
        }
            
        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void RmaInterface<T>::Acc( Matrix<T>& Z, Int i, Int j )
{
    Acc( const_cast<const Matrix<T>&>( Z ), i, j );
}

template<typename T>
void RmaInterface<T>::Get( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Get" ) )

    // a call to Attach with a non-const DistMatrix must set
    // toBeAttachedForGet_ also, if not then it is assumed that
    // the DistMatrix isn't attached
    if( !toBeAttachedForGet_ )
        LogicError( "Cannot perform this operation as matrix is not attached." );

    const DistMatrix<T>& X = *GlobalArrayGet_;
    const Grid& g = X.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myRow = g.Row();
    const Int myCol = g.Col();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    // local width and height
    const Int dm_height = X.Height();
    const Int height = Z.Height();
    const Int width = Z.Width();

    if( i + height > X.Height() || j + width > X.Width() )
        LogicError( "Submatrix out of bounds of global matrix" );

    const Int colAlign = ( X.ColAlign() + i ) % r;
    const Int rowAlign = ( X.RowAlign() + j ) % c;
    const Int XLDim = X.LDim();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;

        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            const Int index = NextIndex( numEntries,
                              getVector_[destination] );
            T* getBuffer = getVector_[destination][index].data();
	    const Int remoteHeight = Length( dm_height, receivingRow, X.ColAlign(), r );
	    //const Int remoteHeight = Length( dm_height, receivingRow, colAlign, r );

	    Int iMapped, jMapped;
	    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
	    {
		iMapped = i; 
		jMapped = ( j / c );
	    }
	    else // 2-D process grid 
	    {
		iMapped = ( i / r ); 
		jMapped = ( j / c );
	    }
	    // get
	    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
	    mpi::Iget( getBuffer, numEntries, destination, disp, numEntries, window );
	    // no difference between localflush
	    // and flush for Get
	    mpi::FlushLocal( destination, window );

            // update local matrix
            for( Int t=0; t<localWidth; ++t )
            {
                T* YCol = Z.Buffer( 0,rowShift+t*c );
                const T* XCol = &getBuffer[t * localHeight];

                for( Int s = 0; s < localHeight; ++s )
                    YCol[colShift+s*r] = XCol[s];
            }
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

// non-blocking interface
template<typename T>
void RmaInterface<T>::Iput( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Iput" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do rma related checks
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
    const Int dm_height = Y.Height();
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    const T* XBuffer = Z.LockedBuffer();

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
            const Int index =
                NextIndex(
                    numEntries,
                    putVector_[destination] );
            T* sendBuffer = putVector_[destination][index].data();
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );

	    Int iMapped, jMapped;
	    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
	    {
		iMapped = i; 
		jMapped = ( j / c );
	    }
	    else // 2-D process grid 
	    {
		iMapped = ( i / r ); 
		jMapped = ( j / c );
	    }
            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
	    }
	    // put
	    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
	    mpi::Iput( sendBuffer, numEntries, destination, disp, numEntries, window );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

// accumulate = Update Y(i:i+height-1,j:j+width-1) += X,
// where X is height x width
template<typename T>
void RmaInterface<T>::Iacc( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Iacc" ) )

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated." );

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative." );

    DistMatrix<T>& Y = *GlobalArrayPut_;

    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError( "Submatrix out of bounds of global matrix." );

    //TODO rma related checks
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = ( Y.ColAlign() + i ) % r;
    const Int rowAlign = ( Y.RowAlign() + j ) % c;
    const Int XLDim = Z.LDim();
    const Int YLDim = Y.LDim();
    // local matrix width and height
    const Int dm_height = Y.Height();
    const Int height = Z.Height();
    const Int width = Z.Width();
    const T* XBuffer = Z.LockedBuffer();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

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
            const Int index =
                NextIndex(
                    numEntries,
                    putVector_[destination] );
            T* sendBuffer = putVector_[destination][index].data();
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );

	    Int iMapped, jMapped;
	    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
	    {
		iMapped = i; 
		jMapped = ( j / c );
	    }
	    else // 2-D process grid 
	    {
		iMapped = ( i / r ); 
		jMapped = ( j / c );
	    }

            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                const T* thisXCol = &XBuffer[( rowShift+t*c )*XLDim];

                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
	    }
	    // acc
	    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
	    mpi::Iacc( sendBuffer, numEntries, destination, disp, numEntries, window );
        }

        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void RmaInterface<T>::Iput( Matrix<T>& Z, Int i, Int j )
{
    Iput( const_cast<const Matrix<T>&>( Z ), i, j );
}

template<typename T>
void RmaInterface<T>::Iacc( Matrix<T>& Z, Int i, Int j )
{
    Iacc( const_cast<const Matrix<T>&>( Z ), i, j );
}

// Local completion of all ops upon
// return
template<typename T>
void RmaInterface<T>::LocalFlush()
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::LocalFlush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    mpi::FlushLocal( window );
}

// Local completion (specific to Z) upon
// return
template<typename T>
void RmaInterface<T>::LocalFlush( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::LocalFlush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    // if there are no request based RMA pending
    // for Z, then this functions acts like Flush
    // local all
    if( !anyPendingXfers( Z ) )
        LocalFlush();
    else
        Wait( Z );
}

// there is no use as of now in
// passing Z, as mpi3 flush enforces
// completion of *all* operations on
// process window
template<typename T>
void RmaInterface<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Flush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    mpi::Flush( window );
}

template<typename T>
void RmaInterface<T>::Flush()
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Flush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    mpi::Flush( window );
}

template<typename T>
bool RmaInterface<T>::anyPendingXfers( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::anyPendingXfers" ) )
    // by default, number of matrices
    // == number of processes
    Int matrixIndex;
    const Int numMatrices = matrices_.size();
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

    // matrix not found
    if( matrixIndex == numMatrices )
        return false;

    return true;
}

// waitany implementation
// cannot use mpi::Waitany
// as of now because request
// objects are vector of deques
template<typename T>
void RmaInterface<T>::WaitAny( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::WaitAny" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex;
    const Int numMatrices = matrices_.size();
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

    // matrix not found
    if( matrixIndex == numMatrices )
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
}

template<typename T>
void RmaInterface<T>::Wait( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Wait" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex;
    const Int numMatrices = matrices_.size();
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

    // matrix not found
    if( matrixIndex == numMatrices )
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
}

template<typename T>
void RmaInterface<T>::Waitall()
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Waitall" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex;
    const Int numMatrices = matrices_.size();

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
}

template<typename T>
bool RmaInterface<T>::Test( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Test" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex;
    const Int numMatrices = matrices_.size();
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

    // matrix not found
    if( matrixIndex == numMatrices )
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

    return true;
}

// TODO Use mpi::Testany instead of mpi::Test
// at present request object is vector
// of deques, so cannot convert it to
// an array required by Testany
template<typename T>
bool RmaInterface<T>::TestAny( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::TestAny" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    Int matrixIndex;
    const Int numMatrices = matrices_.size();
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

    // matrix not found
    if( matrixIndex == numMatrices )
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

    return false;
}

template<typename T>
bool RmaInterface<T>::Testall()
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Testall" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
        LogicError( "Must initiate transfer at first." );

    const Grid& g = ( toBeAttachedForPut_ ?
                      GlobalArrayPut_->Grid() :
                      GlobalArrayGet_->Grid() );
    const Int p = g.Size();
    const Int numMatrices = matrices_.size();

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

    return true;
}

// atomic operations
template<typename T>
T RmaInterface<T>::AtomicIncrement( Int i, Int j, T incr )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::AtomicIncrement" ) )

    if( i < 0 || j < 0 )
        LogicError( "Submatrix offsets must be non-negative" );

    if( !toBeAttachedForPut_ )
        LogicError( "Global matrix cannot be updated" );

    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int dm_height = Y.Height();
    // owner of the cell
    const Int owner = Y.Owner (i, j);
    // calculate height of owner to calculate remote 
    // displacement for RMA
    const Int remoteHeight = Length( dm_height, Y.RowOwner( i ), Y.ColAlign(), r );
    Int iMapped, jMapped;
    if (remoteHeight == dm_height) // 1-D process grid, each PE gets a column strip
    {
	iMapped = i; 
	jMapped = ( j / c );
    }
    else // 2-D process grid 
    {
	iMapped = ( i / r ); 
	jMapped = ( j / c );
    }

    // displacement
    mpi::Aint disp = (iMapped + jMapped * remoteHeight) * sizeof( T );
    // fetch from owner and increment using MPI-3 fetch-and-op
    T value = mpi::ReadInc( window, disp, incr, owner );

    return value;
}

template<typename T>
void RmaInterface<T>::Detach()
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Detach" ) )

    // destructor will call detach again...
    if( detached_ )
        return;

    if( !attached_ )
        LogicError( "Must attach before detaching." );

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
    
    putVector_.clear();
    getVector_.clear();
    
    matrices_.clear();
    
    // data window
    mpi::WindowUnlock( window );
    mpi::WindowFree( window );
}

#define PROTO(T) template class RmaInterface<T>;
#include "El/macros/Instantiate.h"

} // namespace El
#endif // EL_ENABLE_RMA_AXPY

