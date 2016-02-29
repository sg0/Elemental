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

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY)  
namespace El
{
// constructor
template<typename T>
RmaInterface<T>::RmaInterface()
    : GlobalArrayPut_( 0 ), GlobalArrayGet_( 0 ),
      window( mpi::WIN_NULL ),
      putVector_( 0 ), getVector_( 0 ),
      toBeAttachedForPut_( false ), toBeAttachedForGet_( false ),
      attached_( false ), detached_( true ),
      anyPendingGets_( false ), pending_gets_( 0 )
{}

template<typename T>
RmaInterface<T>::RmaInterface( DistMatrix<T>& X )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::RmaInterface" ) )
   
    attached_ 			= false;
    detached_ 			= true;
    toBeAttachedForGet_ 	= false;
    toBeAttachedForPut_ 	= false;
    anyPendingGets_		= false;
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
    anyPendingGets_		= false;
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
	anyPendingGets_			= false;

	const Grid& g = Z.Grid();
        const Int p = g.Size();

        if( putVector_.empty() )
        {
            getVector_.resize( p );
            putVector_.resize( p );
        }
        
	// creation of window
#if defined(EL_USE_WIN_CREATE_FOR_RMA) && \
	!defined(EL_USE_WIN_ALLOC_FOR_RMA)
        const Int numEntries = Z.LocalHeight() * Z.LocalWidth();
        T * baseptr = Z.Buffer();
       
	mpi::WindowCreate( baseptr, numEntries, Z.DistComm(), window );
        mpi::WindowLock( window );
#endif
	// allocation of window
#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
	const Int numEntries = Z.LocalHeight() * Z.LocalWidth();
	if( !Z.ForRMA() ) // fall back to win create
	{
	    T * baseptr = Z.Buffer();
	    mpi::WindowCreate( baseptr, numEntries, Z.DistComm(), window );
	}
	else // DM constructor with RMA flag set, hence DM not allocated
	{
	    mpi::WindowAllocate<T>( numEntries, Z.DistComm(), window );
	    T * baseptr = reinterpret_cast<T *>( mpi::GetWindowBase( window ) );
	    Z.SetWindowBase( baseptr );
	}

	mpi::WindowLock( window );
#endif
	mpi::Barrier( Z.DistComm() );
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
	anyPendingGets_		= false;

        const Grid& g = X.Grid();
        const Int p = g.Size();

        if( getVector_.size() != p )
            getVector_.resize( p );

        // window creation
#if defined(EL_USE_WIN_CREATE_FOR_RMA) &&\
	!defined(EL_USE_WIN_ALLOC_FOR_RMA)
        const Int numEntries = X.LocalHeight() * X.LocalWidth();
        T * baseptr = const_cast<T *>( X.LockedBuffer() );
        mpi::WindowCreate( baseptr, numEntries, X.DistComm(), window );
        mpi::WindowLock( window );
#endif
        // window allocation
#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)       
	const Int numEntries = X.LocalHeight() * X.LocalWidth();
	if( !X.ForRMA() ) // fall-back to window creation
	{
	    T * baseptr = const_cast<T *>( X.LockedBuffer() );
	    mpi::WindowCreate( baseptr, numEntries, X.DistComm(), window );
	    mpi::WindowLock( window );
	}
	else
	    LogicError ("Const DistMatrix cannot be modified, select EL_USE_WIN_CREATE_FOR_RMA");
	// FIXME throwing LogicError in this case is unjustified, fix this
	/*
	{
	    mpi::WindowAllocate( bufferSize, X.DistComm(), window );
	    T * baseptr = reinterpret_cast<T *>( mpi::GetWindowBase( window ) );
	    X.SetWindowBase( baseptr );
	    mpi::WindowLock( window );
	}
	*/
#endif
	mpi::Barrier( X.DistComm() );
    }
}

template<typename T>
Int RmaInterface<T>::NextIndex
( Int dataSize,
  std::deque <std::vector<T>>& dataVectors )
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "RmaInterface::NextIndex" )
	    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	    LogicError( "No DistMatrix attached" );
	    if (dataSize < 0)
	    LogicError ("Resize factor cannot be negative");
	    )
    
    const Int numCreated = dataVectors.size();
    dataVectors.resize( numCreated + 1 );
    dataVectors[numCreated].resize( dataSize );

    return numCreated;
}

// Locally Blocking
template<typename T>
void RmaInterface<T>::Put( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Put" ) )

    Iput( Z, i, j );
    LocalFlush( Z );
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

    Iacc( Z, i, j );
    LocalFlush( Z );
}

template<typename T>
void RmaInterface<T>::Acc( Matrix<T>& Z, Int i, Int j )
{
    Acc( const_cast<const Matrix<T>&>( Z ), i, j );
}

template<typename T>
void RmaInterface<T>::Iget( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "RmaInterface::Iget" )
	    if( i < 0 || j < 0 )
	    LogicError( "Submatrix offsets must be non-negative" );
	    if( !toBeAttachedForGet_ )
	    LogicError( "Cannot perform this operation as matrix is not attached." );
	    )

    const DistMatrix<T>& Y = *GlobalArrayGet_;
    anyPendingGets_ = true;
    
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
    const Int dm_width = Y.Width();
    const Int height = Z.Height();
    const Int width = Z.Width();
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    const Int YLDim = Y.LDim();
    T* ZBuffer = Z.Buffer();

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
	    // height of remote PE holding a DM chunk 
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );
            const Int destination = receivingRow + r*receivingCol;
            const Int index =
                NextIndex( numEntries,
                           getVector_[destination] );
	    T * getBuffer = getVector_[destination][index].data();

	    // keep track of the vector for local completion later
	    const Int currentIndex = pending_gets_.size();
	    pending_gets_.push_back( pending_get_() );

	    pending_gets_[currentIndex].base_ = Z.Buffer();
	    pending_gets_[currentIndex].index_ = index;
	    pending_gets_[currentIndex].destination_ = destination;
	    pending_gets_[currentIndex].colShift_ = colShift;
	    pending_gets_[currentIndex].rowShift_ = rowShift;
	    pending_gets_[currentIndex].localHeight_ = localHeight;
	    pending_gets_[currentIndex].localWidth_ = localWidth;

	    Int iMapped, jMapped;
	    bool isbreak = false;
	    for (Int i_ = i; i_ < dm_height && !isbreak; ++i_)
	    {
		for (Int j_ = j; j_ < dm_width && !isbreak; ++j_)
		{
		    if ( Y.Owner( i_, j_ ) == destination )
		    {
			isbreak = true;
			iMapped = i_ / r;
			jMapped = j_ / c;
		    }
		}
	    }

	    // starting displacement
	    const mpi::Aint disp = iMapped + jMapped * remoteHeight;

#if defined(EL_USE_DDT_FOR_RMA) 
	    // create a new UDD for MPI contig type
	    mpi::Datatype dtype_;

	    // create type
	    mpi::CreateVectorType<T>( localWidth, localHeight, 
		    remoteHeight, &dtype_ );    
	 
	    // get
	    mpi::Iget( getBuffer, numEntries, destination, 
		    disp, 1, dtype_, window );
	
	    mpi::DestroyType( &dtype_ );
#else
	    for( Int t = 0; t < localWidth; ++t )
	    {
		T* getCol = &getBuffer[t * localHeight];

		// displacement of column chunks
		mpi::Aint t_disp = disp + t * remoteHeight;
		
		// get
		mpi::Iget( getCol, localHeight, destination, 
			t_disp, localHeight, window );
	    }
#endif
        }
            
        receivingRow = ( receivingRow + 1 ) % r;

        if( receivingRow == 0 )
            receivingCol = ( receivingCol + 1 ) % c;
    }
}

template<typename T>
void RmaInterface<T>::Get( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Get" ) )

    Iget( Z, i, j );
    LocalFlush( Z );
}

// non-blocking interface
template<typename T>
void RmaInterface<T>::Iput( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "RmaInterface::Iput" )
	    if( i < 0 || j < 0 )
	    LogicError( "Submatrix offsets must be non-negative" );
	    if( !toBeAttachedForPut_ )
	    LogicError( "Global matrix cannot be updated" );
	    )

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
    const Int dm_width = Y.Width();
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
	    // find the number of entries to send
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );
            const Int destination = receivingRow + r*receivingCol;
            const Int index =
                NextIndex( numEntries,
                           putVector_[destination] );
            T* sendBuffer = putVector_[destination][index].data();

	    Int iMapped, jMapped;
	    bool isbreak = false;
	    for (Int i_ = i; i_ < dm_height && !isbreak; ++i_)
	    {
		for (Int j_ = j; j_ < dm_width && !isbreak; ++j_)
		{
		    if ( Y.Owner( i_, j_ ) == destination )
		    {
			isbreak = true;
			iMapped = i_ / r;
			jMapped = j_ / c;
		    }
		}
	    }

	    // starting displacement
	    const mpi::Aint disp = iMapped + jMapped * remoteHeight;

#if defined(EL_USE_DDT_FOR_RMA) 
            // create a new UDD for subarray type
	    mpi::Datatype dtype_;
 
	    // create type
	    mpi::CreateVectorType<T>( localWidth, localHeight, 
		    remoteHeight, &dtype_ );
	    
	    for( Int t = 0; t < localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t * localHeight];
		const T* thisXCol = &XBuffer[( rowShift + t * c ) * ZLDim];

		for( Int s = 0; s < localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift + s * r];
	    }
	    
	    // put
	    mpi::Iput( sendBuffer, numEntries, destination, 
	    	    disp, 1, dtype_, window );
	    
	    mpi::DestroyType( &dtype_ );
#else
	    for( Int t = 0; t < localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t * localHeight];
		const T* thisXCol = &XBuffer[( rowShift + t * c ) * ZLDim];

		for( Int s = 0; s < localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift + s * r];

		// displacement of column chunks
		mpi::Aint t_disp = disp + t * remoteHeight;

		// locally nonblocking put
		mpi::Iput( thisSendCol, localHeight, destination, 
			t_disp, localHeight, window );
	    }
#endif
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
    DEBUG_ONLY( 
	    CallStackEntry cse( "RmaInterface::Iacc" )
	    if( i < 0 || j < 0 )
	    LogicError( "Submatrix offsets must be non-negative" );
	    if( !toBeAttachedForPut_ )
	    LogicError( "Global matrix cannot be updated" );
	    )

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
    const Int dm_width = Y.Width();
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
	    // find the number of entries to send
	    const Int remoteHeight = Length( dm_height, receivingRow, Y.ColAlign(), r );
            const Int destination = receivingRow + r*receivingCol;
            const Int index =
                NextIndex( numEntries,
                           putVector_[destination] );
            T* sendBuffer = putVector_[destination][index].data();

	    Int iMapped, jMapped;
	    bool isbreak = false;
	    for (Int i_ = i; i_ < dm_height && !isbreak; ++i_)
	    {
		for (Int j_ = j; j_ < dm_width && !isbreak; ++j_)
		{
		    if ( Y.Owner( i_, j_ ) == destination )
		    {
			isbreak = true;
			iMapped = i_ / r;
			jMapped = j_ / c;
		    }
		}
	    }

	    // starting displacement
	    const mpi::Aint disp = iMapped + jMapped * remoteHeight;

#if defined(EL_USE_DDT_FOR_RMA) 
            // create a new UDD for subarray type
	    mpi::Datatype dtype_;

	    // create type
	    mpi::CreateVectorType<T>( localWidth, localHeight, 
		    remoteHeight, &dtype_ );    
	    
	    for( Int t = 0; t < localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t * localHeight];
		const T* thisXCol = &XBuffer[( rowShift + t * c ) * ZLDim];

		for( Int s = 0; s < localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift + s * r];
	    }
	    
	    // acc
	    mpi::Iacc( sendBuffer, numEntries, destination, 
		    disp, 1, dtype_, mpi::SUM, window );
	    
	    mpi::DestroyType( &dtype_ );
#else
	    for( Int t = 0; t < localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t * localHeight];
		const T* thisXCol = &XBuffer[( rowShift + t * c ) * ZLDim];

		for( Int s = 0; s < localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift + s * r];

		// displacement of column chunks
		mpi::Aint t_disp = disp + t * remoteHeight;

		// locally nonblocking acc
		mpi::Iacc( thisSendCol, localHeight, destination, 
			t_disp, localHeight, mpi::SUM, window );
	    }
#endif
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

// flush only puts/accs
template<typename T>
void RmaInterface<T>::LocalFlush( const Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::LocalFlush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    mpi::FlushLocal( window );
}

// flush only puts/accs
template<typename T>
void RmaInterface<T>::Flush( const Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Flush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    mpi::Flush( window );
}

// flush put/accs, and check for gets
template<typename T>
void RmaInterface<T>::LocalFlush( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::LocalFlush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    // local flush on RMA window
    mpi::FlushLocal( window );
    
    // process pending gets if any
    if (anyPendingGets_ == true)
    {
	// DM params
	const DistMatrix<T>& Y = *GlobalArrayGet_;

	const Grid& g = Y.Grid();
	const Int r = g.Height();
	const Int c = g.Width();

	// find the index of Z in pending_gets_
	for (Int mindex = 0; mindex < pending_gets_.size(); mindex++)
	{
	    // compare base address of matrices
	    if (pending_gets_[mindex].base_ == Z.Buffer()
		    && pending_gets_[mindex].is_active_ == true)
	    {
		// get metadata before initiating flush
		const Int destination = pending_gets_[mindex].destination_;
		const Int index = pending_gets_[mindex].index_;
		T* getBuffer = getVector_[destination][index].data();

		const Int colShift = pending_gets_[mindex].colShift_;
		const Int rowShift = pending_gets_[mindex].rowShift_;
		const Int localHeight = pending_gets_[mindex].localHeight_;
		const Int localWidth = pending_gets_[mindex].localWidth_;

		// update local matrix
		for( Int t = 0; t < localWidth; ++t )
		{
		    T* getCol = &getBuffer[t * localHeight];
		    T* YCol = Z.Buffer( 0,rowShift+t*c );
		    const T* XCol = getCol;

		    for( Int s = 0; s < localHeight; ++s )
			YCol[colShift+s*r] = XCol[s];
		}

		pending_gets_[mindex].is_active_ = false;
	    }
	}

	anyPendingGets_ = false;
    }
}

// Local completion (specific to Z) upon
// return
template<typename T>
void RmaInterface<T>::LocalFlush()
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::LocalFlush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    // local flush on RMA window
    mpi::FlushLocal( window );
}

template<typename T>
void RmaInterface<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Flush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    // local/remote flush on RMA window
    mpi::Flush( window );

    // process pending gets if any
    if (anyPendingGets_ == true)
    {
	// DM params
	const DistMatrix<T>& Y = *GlobalArrayGet_;

	const Grid& g = Y.Grid();
	const Int r = g.Height();
	const Int c = g.Width();

	// find the index of Z in pending_gets_
	for (Int mindex = 0; mindex < pending_gets_.size(); mindex++)
	{
	    // compare base address of matrices
	    if (pending_gets_[mindex].base_ == Z.Buffer()
		    && pending_gets_[mindex].is_active_ == true)
	    {
		// get metadata before initiating flush
		const Int destination = pending_gets_[mindex].destination_;
		const Int index = pending_gets_[mindex].index_;
		T* getBuffer = getVector_[destination][index].data();

		const Int colShift = pending_gets_[mindex].colShift_;
		const Int rowShift = pending_gets_[mindex].rowShift_;
		const Int localHeight = pending_gets_[mindex].localHeight_;
		const Int localWidth = pending_gets_[mindex].localWidth_;

		// update local matrix
		for( Int t = 0; t < localWidth; ++t )
		{
		    T* getCol = &getBuffer[t * localHeight];
		    T* YCol = Z.Buffer( 0,rowShift+t*c );
		    const T* XCol = getCol;

		    for( Int s = 0; s < localHeight; ++s )
			YCol[colShift+s*r] = XCol[s];
		}

		pending_gets_[mindex].is_active_ = false;
	    }
	}
	
	anyPendingGets_ = false;
    }
}

template<typename T>
void RmaInterface<T>::Flush()
{
    DEBUG_ONLY( CallStackEntry cse( "RmaInterface::Flush" ) )

    if( !toBeAttachedForPut_ || !toBeAttachedForGet_ )
	LogicError( "Must initiate transfer before flushing." );

    // local/remote flush on RMA window
    mpi::Flush( window );
}

// atomic operations
template<typename T>
T RmaInterface<T>::AtomicIncrement( Int i, Int j, T incr )
{
    DEBUG_ONLY( 
	    CallStackEntry cse( "RmaInterface::AtomicIncrement" )
	    if( i < 0 || j < 0 )
	    LogicError( "Submatrix offsets must be non-negative" );
	    if( !toBeAttachedForPut_ )
	    LogicError( "Global matrix cannot be updated" );
	    )

    DistMatrix<T>& Y = *GlobalArrayPut_;
    
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int dm_height = Y.Height();
    const Int dm_width = Y.Width();
   
    // owner of the cell
    const Int owner = Y.Owner (i, j);
    
    // calculate height of owner to calculate remote 
    // displacement for RMA
    const Int remoteHeight = Length( dm_height, Y.RowOwner( i ), Y.ColAlign(), r );
    const Int iMapped = i / r;
    const Int jMapped = j / c;

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

    const mpi::Comm& comm = ( toBeAttachedForPut_ ?
                              GlobalArrayPut_->DistComm() :
                              GlobalArrayGet_->DistComm() );

    mpi::Flush( window );
    
    attached_ 		= false;
    detached_ 		= true;
    
    toBeAttachedForPut_ = false;
    toBeAttachedForGet_ = false;
    anyPendingGets_	= false;
    
    GlobalArrayPut_ 	= 0;
    GlobalArrayGet_ 	= 0;
    
    putVector_.clear();
    getVector_.clear();

    pending_gets_.clear();
    
    // data window
    mpi::WindowUnlock( window );
    mpi::WindowFree( window );
    window = mpi::WIN_NULL;
}

#define PROTO(T) template class RmaInterface<T>;
#include "El/macros/Instantiate.h"

} // namespace El
#endif // EL_ENABLE_RMA_AXPY
