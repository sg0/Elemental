/*
   Copyright (c) 2009-2014, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   Copyright (c) 2014, Jeff Hammond (Intel)
   Copyright (c) 2014, Sayan Ghosh (University of Houston)
   All rights reserved.

Authors:
Jeff Hammond adapted the RMA interface from the AXPY one.

This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include <assert.h>

// TODO Complete the const interfaces...
// TODO RMA related checks pending (e.g bounds checking)...
// TODO Consider DDT
// TODO Add disp as a parameter to MPI one sided functions?
// TODO Use DEBUG_ONLY or something that EL provides instead of assert
// TODO Make variable names rma friendly
#if MPI_VERSION>=3
namespace El
{

template<typename T>
RmaInterface<T>::RmaInterface()
    : GlobalArrayPut_(0), GlobalArrayGet_(0),
    putVector_(0), getVector_(0), window (MPI_WIN_NULL),
    toBeAttachedForPut_(false), toBeAttachedForGet_(false),
    attached_(false), detached_(false)
{ }

template<typename T>
RmaInterface<T>::RmaInterface( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::RmaInterface"))

    attached_ 		= false;
    detached_ 		= false;
    toBeAttachedForGet_ = true;
    toBeAttachedForPut_ = true;
    GlobalArrayPut_ 	= &Z;
    GlobalArrayGet_ 	= &Z;
    window 	    	= MPI_WIN_NULL;

    const Int p = Z.Grid ().Size();
    putVector_.resize( p );
    getVector_.resize( p );
}

template<typename T>
RmaInterface<T>::RmaInterface( const DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::RmaInterface"))

    attached_ 		= false;
    detached_ 		= false;
    toBeAttachedForGet_ = true;
    toBeAttachedForPut_ = false;
    GlobalArrayGet_ 	= &X;
    GlobalArrayPut_ 	= 0;
    window 	    	= MPI_WIN_NULL;

    const Int p = X.Grid ().Size ();
    getVector_.resize( p );
    putVector_.resize( p );
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
            DEBUG_ONLY(DumpCallStack())
        }
        else
        {
            Detach();
        }
    }
}

template<typename T>
void RmaInterface<T>::Attach( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Attach"))
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
    
    if (putVector_.size() != p)
    {
	getVector_.resize( p );
	putVector_.resize( p );
    }

    // do rma related checks
    const Int numEntries = Z.LocalHeight () * Z.LocalWidth ();
    const Int bufferSize = numEntries * sizeof(T);
    void* baseptr = reinterpret_cast<void*>(Z.Buffer ());
    assert(baseptr != NULL);

    mpi::WindowCreate (baseptr, bufferSize, g.VCComm (), window);
    mpi::WindowLock (window);
}

template<typename T>
void RmaInterface<T>::Attach( const DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Attach"))
    if (!attached_)
        attached_ = true;
    else
        LogicError("Must detach before reattaching.");

    if( !toBeAttachedForGet_ && !toBeAttachedForPut_)
    {
        GlobalArrayGet_ 	= &X;
        toBeAttachedForGet_ 	= true;
        GlobalArrayPut_ 	= 0;
        toBeAttachedForPut_ 	= false;
    }
    else
        LogicError("Cannot update Global matrix.");

    const Grid& g = X.Grid();
    const Int p = g.Size ();

    if (putVector_.size() != p)
    {
	getVector_.resize( p );
	putVector_.resize( p );
    }

    //do rma related checks
    const Int numEntries = X.LocalHeight () * X.LocalWidth ();
    const Int bufferSize = numEntries * sizeof(T);
    void* baseptr = (void*)X.LockedBuffer ();
    assert (baseptr != NULL);

    mpi::WindowCreate (baseptr, bufferSize, g.VCComm (), window);
    mpi::WindowLock (window);
}

template<typename T>
void RmaInterface<T>::Put( T alpha, Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Put"))

	if( i < 0 || j < 0 )
	    LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
	LogicError("Global matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayPut_;
    //do rma related checks
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

    const Int iLocalOffset = Length( i, Y.ColShift (), r );
    const Int jLocalOffset = Length( j, Y.RowShift (), c );

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
	    const Int bufferSize = numEntries * sizeof(T);

	    putVector_[destination].resize( bufferSize );
	    byte* sendBuffer = putVector_[destination].data();
	    T* sendData = reinterpret_cast<T*>(sendBuffer);
	    const T* XBuffer = Z.LockedBuffer();

	    for( Int t=0; t<localWidth; ++t )
	    {
		T* thisSendCol = &sendData[t*localHeight];
		const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
		for( Int s=0; s<localHeight; ++s )
		    thisSendCol[s] = alpha*thisXCol[colShift+s*r];
		// put
		mpi::Aint disp =  (iLocalOffset + (jLocalOffset+t) * Y.LDim ()) * sizeof(T);
		mpi::Iput_nolocalflush (reinterpret_cast<void*>(sendBuffer + t*localHeight*sizeof(T)), 
			localHeight * sizeof(T), destination, disp, localHeight * sizeof(T), window);	    
	    }
	    // local flush, okay to clear buffers after this
	    mpi::FlushLocal (destination, window);
	    // clear
	    putVector_[destination].resize (0);
	}
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void RmaInterface<T>::Put( T alpha, const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Put"))
}

template<typename T>
void RmaInterface<T>::Get( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Get"))
    // a call to Attach with a non-const DistMatrix must set
    // toBeAttachedForGet_ also, if not then it is assumed that
    // the DistMatrix isn't attached
    if ( !toBeAttachedForGet_ )
        LogicError ("Cannot perform this operation as matrix is not attached.");

    const DistMatrix<T> &X = *GlobalArrayGet_;

    const Grid & g = X.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int p = g.Size ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();

    // local width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    if (i + height > X.Height () || j + width > X.Width ())
        LogicError("Submatrix out of bounds of global matrix");

    const Int colAlign = (X.ColAlign() + i) % r;
    const Int rowAlign = (X.RowAlign() + j) % c;
    
    const Int iLocalOffset = Length (i, X.ColShift (), r);
    const Int jLocalOffset = Length (j, X.RowShift (), c);

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
            const Int bufferSize = numEntries * sizeof(T);

            getVector_[destination].resize (bufferSize);
            byte *getBuffer = getVector_[destination].data ();
	    
	    // get
	    for( Int t=0; t<localWidth; ++t )
	    {
		mpi::Aint disp =  (iLocalOffset + (jLocalOffset+t) * X.LDim ()) * sizeof(T);
		mpi::Iget_nolocalflush (reinterpret_cast<void*>(getBuffer + t*localHeight*sizeof(T)), 
			localHeight * sizeof(T), destination, disp, localHeight * sizeof(T), window);	    
	    }
	    // no difference between localflush 
	    // and flush for Get
	    mpi::FlushLocal (destination, window);
	    T* getData = reinterpret_cast<T*>(getBuffer);
	    // update local matrix 	    
	    for( Int t=0; t<localWidth; ++t )
	    {
		T *YCol = Z.Buffer (0,rowShift+t*c);
		const T *XCol = &getData[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[colShift+s*r] = XCol[s];
	    }
	    // clear
	    getVector_[destination].resize (0);
        }
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void RmaInterface<T>::Get( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Get"))
}

// scaled accumulate = Update Y(i:i+height-1,j:j+width-1) += alpha X, 
// where X is height x width
template<typename T>
void RmaInterface<T>::Acc( T alpha, Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Acc"))

    if ( !toBeAttachedForPut_ )
        LogicError("Global matrix cannot be updated.");
    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative.");

    DistMatrix<T>& Y = *GlobalArrayPut_;

    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix.");

    //do rma related checks

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
	
    const Int iLocalOffset = Length( i, Y.ColShift (), r );
    const Int jLocalOffset = Length( j, Y.RowShift (), c );
    
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
	    const Int bufferSize = numEntries * sizeof(T);

	    putVector_[destination].resize( bufferSize );
	    byte* sendBuffer = putVector_[destination].data();
	    T* sendData = reinterpret_cast<T*>(sendBuffer);
	    const T* XBuffer = Z.LockedBuffer();

	    //src *= scale
	    for( Int t=0; t<localWidth; ++t )
	    {
		T* thisSendCol = &sendData[t*localHeight];
		const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
		for( Int s=0; s<localHeight; ++s )
		    thisSendCol[s] = alpha*thisXCol[colShift+s*r];
	    }
	    // acc
	    for( Int t=0; t<localWidth; ++t )
	    {
		mpi::Aint disp =  (iLocalOffset + (jLocalOffset+t) * Y.LDim ()) * sizeof(T);
		//mpi::Iacc_nolocalflush (reinterpret_cast<void*>(sendBuffer + t*localHeight*sizeof(T)), 
		//	localHeight * sizeof(T), destination, disp, localHeight * sizeof(T), window);	    
		mpi::Iacc (reinterpret_cast<void*>(sendBuffer + t*localHeight*sizeof(T)), 
			localHeight * sizeof(T), destination, disp, localHeight * sizeof(T), window);	    
	    }
	    // local flush, okay to clear buffers after this
	    //mpi::FlushLocal (destination, window);
	    // clear
	    putVector_[destination].resize (0);
	}
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void RmaInterface<T>::Acc( T alpha, const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Acc"))
}

template<typename T>
void RmaInterface<T>::Flush( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Flush"))
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
        LogicError("Must initiate transfer before flushing.");

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do rma related checks
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;

    // local width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    // find destination
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight*localWidth;

        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            mpi::Flush ( destination, window );
        }

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void RmaInterface<T>::Flush( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Flush"))
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
        LogicError("Must initiate transfer before flushing.");

    const DistMatrix<T>& Y = *GlobalArrayGet_;

    //do rma related checks
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;

    // local width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    // find destination
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight*localWidth;

        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            mpi::Flush ( destination, window );
        }

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

// Are these only useful when the user wants to 
// get/put the entire DistMatrix to it's local 
// PE/everyone in world ?
template<typename T>
void RmaInterface<T>::Flush( Matrix<T>& Z )
{

    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Flush"))
    
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
        LogicError("Must initiate transfer before flushing.");

    DistMatrix<T>& Y = *GlobalArrayPut_;
    
    //do rma related checks
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    // i = j = 0 - leftmost coordinates of DistMatrix
    const Int colAlign = Y.ColAlign() % r;
    const Int rowAlign = Y.RowAlign() % c;

    // local width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    // find destination
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight*localWidth;

        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            mpi::Flush ( destination, window );
        }

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void RmaInterface<T>::Flush( const Matrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Flush"))

    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
        LogicError("Must initiate transfer before flushing.");

    // rma checks, see if Z is not NULL, etc
    const DistMatrix<T>& Y = *GlobalArrayGet_;

    //do rma related checks
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    // i = j = 0 - leftmost coordinates of DistMatrix
    const Int colAlign = Y.ColAlign() % r;
    const Int rowAlign = Y.RowAlign() % c;

    // local width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    // find destination
    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight*localWidth;

        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            mpi::Flush ( destination, window );
        }

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void RmaInterface<T>::Detach()
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Detach"))
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

    putVector_.clear();
    getVector_.clear();

    mpi::WindowUnlock (window);
    mpi::WindowFree (window);
}

template class RmaInterface<Int>;
template class RmaInterface<float>;
template class RmaInterface<double>;
template class RmaInterface<Complex<float>>;
template class RmaInterface<Complex<double>>;

} // namespace El
#endif
