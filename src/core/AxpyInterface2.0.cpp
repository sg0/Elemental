/*
This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include <cassert>

// TODO Use DDT for put/get/acc when EL_USE_DERIVED_TYPE is defined
namespace El
{
template<typename T>
AxpyInterface2<T>::AxpyInterface2()
    : GlobalArrayPut_(0), GlobalArrayGet_(0),
    putVectors_(0), getVectors_(0), dataRequests_(0),
    dataRequestStatuses_(0), matrixBase_(0), opKind_(0),
    toBeAttachedForPut_(false), toBeAttachedForGet_(false),
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

    const Int p = Z.Grid ().Size();
    putVectors_.resize( p );
    getVectors_.resize( p );
    dataRequests_.resize (p);
    dataRequestStatuses_.resize (p);
    matrixBase_.resize (p);
    opKind_.resize (p);
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
	opKind_.resize (p);
    }
}

template<typename T>
Int AxpyInterface2<T>::NextIndex
( Int dataSize, std::deque <std::vector<T>> &dataVectors,
  std::deque <mpi::Request> &requests,
  std::deque <bool> &requestStatus,
  std::deque <Int> &opKind,
  Int op,
  std::deque <T *> &matrixBase,
  T * base_address)
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndex"))
    const Int Index = Int(requests.size ());
    
    dataVectors.resize (Index + 1);
    dataVectors[Index].resize (dataSize);
    requests.push_back (mpi::REQUEST_NULL);
    requestStatus.push_back (true);
    opKind.push_back (op);
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
		    opKind_[destination],
		    DATA_PUT_TAG,
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

    const Int XLDim = X.LDim ();

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
            const Int source = receivingRow + r*receivingCol;	    
	    T* XBuffer = Z.Buffer();
	    const Int index =
	    NextIndex (numEntries, getVectors_[source],
		    dataRequests_[source],
		    dataRequestStatuses_[source],
		    opKind_[source],
		    DATA_GET_TAG,
		    matrixBase_[source],
		    XBuffer);
            
	    DEBUG_ONLY (if
			(Int (getVectors_[source][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)
	    T *getBuffer = reinterpret_cast<T *>(getVectors_[source][index].data ());
            // get request
	    mpi::TaggedIRecv (getBuffer, numEntries, source, DATA_GET_TAG, 
		    g.VCComm (), dataRequests_[source][index]);
	}
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
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
		    opKind_[destination],
		    DATA_ACC_TAG,
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

// local accumulate, Z += Get Y(i:i+height-1,j:j+width-1),
// where Z is local matrix height x width
template<typename T>
void AxpyInterface2<T>::LocalAcc( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::LocalAcc"))
    // a call to Attach with a non-const DistMatrix must set
    // toBeAttachedForGet_ also, if not then it is assumed that
    // the DistMatrix isn't attached
    if ( !toBeAttachedForGet_ )
        LogicError ("Cannot perform this operation as matrix is not attached.");
    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");

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

    const Int XLDim = X.LDim ();

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
            const Int source = receivingRow + r*receivingCol;
 	    T *XBuffer = Z.Buffer ();
	    const Int index =
	    NextIndex (numEntries, getVectors_[source],
		    dataRequests_[source],
		    dataRequestStatuses_[source],
		    opKind_[source],
		    DATA_LCC_TAG,
		    matrixBase_[source],
		    XBuffer);

	    DEBUG_ONLY (if
			(Int (getVectors_[source][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)
	    T *getBuffer = reinterpret_cast<T *>(getVectors_[source][index].data ());
            // get request
	    mpi::TaggedRecv (getBuffer, numEntries, source, 
		    DATA_LCC_TAG, g.VCComm ());
            // acc to local matrix
            for( Int t=0; t<localWidth; ++t )
            {
                T *YCol = Z.Buffer (0,rowShift+t*c);
                const T *XCol = &getBuffer[t * localHeight];
                for (Int s = 0; s < localHeight; ++s)
                    YCol[colShift+s*r] += XCol[s];
            }
        }
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

// get index of a matrix 
template<typename T>
Int AxpyInterface2<T>::GetIndexForMatrix ( Matrix<T>& Z, const Int rank )
{
    typename std::deque<T *>::iterator dit;
    dit = std::find ( matrixBase_[rank].begin(), 
	    matrixBase_[rank].end(), Z.LockedBuffer ());
    const Int index = (dit - matrixBase_[rank].begin());
    //std::cout << "matrixBase size: " << matrixBase_[rank].size () << "\n";
    assert ( index != matrixBase_[rank].size () );

    return index;
}

// get operation associated with a matrix 
template<typename T>
Int AxpyInterface2<T>::GetMatrixType ( Matrix<T>& Z, const Int rank )
{
    const Int index = GetIndexForMatrix ( Z, rank );
    return opKind_[rank][index];
}

// progress communication for a particular matrix
// this could be used to progress sends and recvs
template<typename T>
void AxpyInterface2<T>::ProgressMatrix ( Matrix<T>& Z, const Int rank )
{
    const Int index = GetIndexForMatrix ( Z, rank );
    if ( !dataRequestStatuses_[rank][index] ) // nothing to do
	return;
    // wait
    mpi::Wait ( dataRequests_[rank][index] );
    dataRequestStatuses_[rank][index] = false;
    getVectors_[rank][index].resize (0);
    putVectors_[rank][index].resize (0);
}

// local matrix could be updated after local flush
template<typename T>
void AxpyInterface2<T>::LocalFlush( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::LocalFlush"))
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
        LogicError("Must initiate transfer before flushing.");

    DistMatrix<T>& Y = *GlobalArrayPut_;

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
	    // local flush has no meaning for global to local
	    Int type = GetMatrixType ( Z, destination );
	    if ( type == DATA_GET_TAG || type == DATA_LCC_TAG )
	    {
		Flush ( Z, i, j );
		return;
	    }
	    ProgressMatrix ( Z, destination );
        }

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
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
	    // ensure local completion for local to 
	    // global transfers
	    const Int type = GetMatrixType ( Z, destination );
	    switch ( type )
	    {
		case DATA_PUT_TAG:
		    {
			ProgressMatrix ( Z, destination );
			HandleLocalToGlobalData ( Z, i, j );
			break;
		    }
		case DATA_ACC_TAG:
		    {
			ProgressMatrix ( Z, destination );
			HandleLocalToGlobalAcc ( Z, i, j );
			break;
		    }
		case DATA_GET_TAG:
		    {
			HandleGlobalToLocalData ( Z, i, j );
			break;
		    }
		case DATA_LCC_TAG:
		    {
			HandleGlobalToLocalAcc ( Z, i, j );
			break;
		    }
	    }
	}
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Flush"))
    Flush ( Z, 0, 0 );	
}

    template < typename T >
void AxpyInterface2<T>::HandleGlobalToLocalData ( Matrix<T>& Z, Int i, Int j )
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
		    opKind_[destination],
		    DATA_PUT_TAG,
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
	    // Fire off non-blocking send
	    mpi::TaggedSend (sendBuffer, numEntries, destination, 
		    DATA_GET_TAG, g.VCComm ());
	    // clear
	    ProgressMatrix ( Z, destination );
	}
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }
}

// replica of above, except the tag is different
    template < typename T >
void AxpyInterface2<T>::HandleGlobalToLocalAcc ( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface::HandleGlobalToLocalAcc"))

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
		    opKind_[destination],
		    DATA_PUT_TAG,
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
	    // Fire off non-blocking send
	    mpi::TaggedSend (sendBuffer, numEntries, destination, 
		    DATA_LCC_TAG, g.VCComm ());
	    // clear
	    ProgressMatrix ( Z, destination );
	}
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }
}

template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalData ( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface::HandleLocalToGlobalData"))

    if( i < 0 || j < 0 )
	    LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
	LogicError("Global matrix cannot be updated");   

    const DistMatrix<T> &X = *GlobalArrayPut_;

    const Grid & g = X.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int p = g.Size ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    std::vector < std::vector<T> > getVector;
    getVector.resize (p);

    // local width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    if (i + height > X.Height () || j + width > X.Width ())
	LogicError("Submatrix out of bounds of global matrix");

    const Int colAlign = (X.ColAlign() + i) % r;
    const Int rowAlign = (X.RowAlign() + j) % c;
    const Int XLDim = X.LDim ();

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
            T* XBuffer = Z.Buffer();
     	    const Int index =
	    NextIndex (numEntries, putVectors_[destination],
		    dataRequests_[destination],
		    dataRequestStatuses_[destination],
		    opKind_[destination],
		    DATA_GET_TAG,
		    matrixBase_[destination], 
		    XBuffer);
	    DEBUG_ONLY (if
			(Int (getVectors_[destination][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)
	    T *getBuffer = getVectors_[destination][index].data ();

	    mpi::TaggedRecv (getBuffer, numEntries, destination, 
		    DATA_PUT_TAG, g.VCComm ());
	    // update local matrix
	    for( Int t=0; t<localWidth; ++t )
	    {
		T *YCol = Z.Buffer (0,rowShift+t*c);
		const T *XCol = &getBuffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[colShift+s*r] = XCol[s];
	    }
	    // clear
	    ProgressMatrix ( Z, destination );
	}
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }
}
// replica of above function except this accumulates
// we need this function to avoid an Iprobe 
template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalAcc ( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface::HandleLocalToGlobalData"))

    if( i < 0 || j < 0 )
	    LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
	LogicError("Global matrix cannot be updated");   
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

    const Int XLDim = X.LDim ();

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
            T* XBuffer = Z.Buffer();
     	    const Int index =
	    NextIndex (numEntries, putVectors_[destination],
		    dataRequests_[destination],
		    dataRequestStatuses_[destination],
		    opKind_[destination],
		    DATA_GET_TAG,
		    matrixBase_[destination], 
		    XBuffer);
	    DEBUG_ONLY (if
			(Int (getVectors_[destination][index].size ()) !=
			 numEntries) LogicError ("Error in NextIndex");)
	    T *getBuffer = getVectors_[destination][index].data ();

	    mpi::TaggedRecv (getBuffer, numEntries, destination, 
		    DATA_ACC_TAG, g.VCComm ());
	    // update local matrix
	    for( Int t=0; t<localWidth; ++t )
	    {
		T *YCol = Z.Buffer (0,rowShift+t*c);
		const T *XCol = &getBuffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[colShift+s*r] += XCol[s];
	    }
	    // clear
	    ProgressMatrix ( Z, destination );
	}
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }
}

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
    matrixBase_.clear();
    opKind_.clear();
}

template class AxpyInterface2<Int>;
template class AxpyInterface2<float>;
template class AxpyInterface2<double>;
template class AxpyInterface2<Complex<float>>;
template class AxpyInterface2<Complex<double>>;

} // namespace El
