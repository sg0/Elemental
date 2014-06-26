/*
   Copyright (c) 2009-2014, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   Copyright (c) 2014, Jeff Hammond (Intel)
   All rights reserved.

   Authors:
   Jeff Hammond adapted the RMA interface from the AXPY one.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
 
// This is direct copy-paste from
// El two-sided implementation with
// point-to-point replaced by one-sided

// If you're seeing this then at this 
// point I just want to compile

#if MPI_VERSION>=3
namespace El {

    // dont care about const 
    // interfaces now
template<typename T>
RmaInterface<T>::RmaInterface( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::RmaInterface"))
	const Int p = Z.Grid ().Size ();
}

template<typename T>
RmaInterface<T>::RmaInterface( const DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::RmaInterface"))
	const Int p = X.Grid ().Size ();
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
    
}

template<typename T>
void RmaInterface<T>::Attach( const DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Attach"))
}

template<typename T>
void RmaInterface<T>::Put( T alpha, Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Put"))

    DistMatrix<T>& Y = *localToGlobalMat_;
    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix");

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

    // put local matrix cells in
    // correct places in global array
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
            const Int bufferSize = 4*sizeof(Int) + (numEntries+1)*sizeof(T);
            // Pack the header
	    // make variable names rma friendly
            byte* sendBuffer;
            byte* head = sendBuffer;
            *reinterpret_cast<Int*>(head) = i; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = j; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = height; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = width; head += sizeof(Int);
            *reinterpret_cast<T*>(head) = alpha; head += sizeof(T);

            // Pack the payload
	    // consider ddt here
            T* sendData = reinterpret_cast<T*>(head);
            const T* XBuffer = Y.LockedBuffer();
            const Int XLDim = Y.LDim();
            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendData[t*localHeight];
                const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

	    mpi::Iput (sendBuffer, bufferSize, destination, bufferSize, window);
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
	const DistMatrix < T > &X = *globalToLocalMat_;

    const Int height = Z.Height ();
    const Int width = Z.Width ();
    if (i + height > X.Height () || j + width > X.Width ())
	LogicError ("Invalid AxpyGlobalToLocal submatrix");

    const Grid & g = X.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int p = g.Size ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();

    for (Int rank = 0; rank < p; ++rank)
    {
	// this is horrendously wrong, but 
	// just for compiling
	const Int buffersize = height * width * sizeof(T);
	getVector_.resize (buffersize);
	byte *getBuffer = getVector_.data ();
	// how do we know the data size 
	mpi::Iget (getBuffer, buffersize, rank, buffersize, window);
	// Extract the header
	byte *head = getBuffer;
	const Int i = *reinterpret_cast < const Int * >(head);
	head += sizeof (Int);
	const Int j = *reinterpret_cast < const Int * >(head);
	head += sizeof (Int);
	const Int height = *reinterpret_cast < const Int * >(head);
	head += sizeof (Int);
	const Int width = *reinterpret_cast < const Int * >(head);
	head += sizeof (Int);
	const T alpha = *reinterpret_cast < const T * >(head);
	head += sizeof (T);  

	// Update Y
	const T *XBuffer = reinterpret_cast < const T * >(head);
	const Int colAlign = (X.ColAlign () + i) % r;
	const Int rowAlign = (X.RowAlign () + j) % c;
	const Int colShift = Shift (myRow, colAlign, r);
	const Int rowShift = Shift (myCol, rowAlign, c);

	const Int localHeight = Length (height, colShift, r);
	const Int localWidth = Length (width, rowShift, c);
	const Int iLocalOffset = Length (i, X.ColShift (), r);
	const Int jLocalOffset = Length (j, X.RowShift (), c);

	for (Int t = 0; t < localWidth; ++t)
	{
	    T *YCol = Z.Buffer (iLocalOffset, jLocalOffset + t);
	    const T *XCol = &XBuffer[t * localHeight];
	    for (Int s = 0; s < localHeight; ++s)
		YCol[s] += alpha * XCol[s];
	}
    }
}

template<typename T>
void RmaInterface<T>::Get( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Get"))
}

// scaled accumulate
template<typename T>
void RmaInterface<T>::Acc( T alpha, Matrix<T>& Z, mpi::Op &op, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Acc"))

    DistMatrix<T>& Y = *localToGlobalMat_;
    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix");

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

    // put local matrix cells in
    // correct places in global array
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
            const Int bufferSize = 4*sizeof(Int) + (numEntries+1)*sizeof(T);
            // Pack the header
	    // make variable names rma friendly
            byte* sendBuffer;
            byte* head = sendBuffer;
            *reinterpret_cast<Int*>(head) = i; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = j; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = height; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = width; head += sizeof(Int);
            *reinterpret_cast<T*>(head) = alpha; head += sizeof(T);

            // Pack the payload
	    // consider ddt here
            T* sendData = reinterpret_cast<T*>(head);
            const T* XBuffer = Z.LockedBuffer();
            const Int XLDim = Z.LDim();
            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendData[t*localHeight];
                const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

	    mpi::Iacc (sendBuffer, bufferSize, destination, bufferSize, op, window);
	}

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void RmaInterface<T>::Acc( T alpha, const Matrix<T>& Z, mpi::Op &op, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Acc"))
}

template<typename T>
void RmaInterface<T>::Flush( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Flush"))
}

template<typename T>
void RmaInterface<T>::Detach()
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Detach"))
}

template class RmaInterface<Int>;
template class RmaInterface<float>;
template class RmaInterface<double>;
template class RmaInterface<Complex<float>>;
template class RmaInterface<Complex<double>>;

} // namespace El
#endif
