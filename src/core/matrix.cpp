/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

// Constructors and destructors
// ============================

// Create a 0x0 matrix
template<typename T>
Matrix<T>::Matrix( bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ), lfm_(0,0,1,nullptr)
{ }

// Create a matrix with the specified dimensions
// TODO: Intelligent choice of leading dimension?
template<typename T>
Matrix<T>::Matrix( DimPair dims, bool fixed )
: Matrix(Layout(dims,Max(dims.m,1)),fixed) { }
template<typename T>
Matrix<T>::Matrix( Int height, Int width, bool fixed )
: Matrix(Layout(height,width,fixed)) { }

// Create a matrix with the specified dimensions and leading dimension
template<typename T>
Matrix<T>::Matrix( Layout layout, bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Matrix");
        Check( layout );
    )
    memory_.Require( layout.LDim()*layout.Width() );
    lfm_.layout = layout;
    lfm_.buffer = memory_.Buffer();
}
template<typename T>
Matrix<T>::Matrix( Int height, Int width, Int ldim, bool fixed )
: Matrix(Layout(height,width,ldim),fixed) { }

template<typename T>
Matrix<T>::Matrix( const elem::LockedFlatMatrix<T>& lfm, bool fixed )
: viewType_( fixed ? LOCKED_VIEW_FIXED: LOCKED_VIEW ), lfm_(lfm)
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Matrix");
        Check( lfm.layout );
    )
}
// Deprecated
template<typename T>
Matrix<T>::Matrix
( Int height, Int width, const T* buffer, Int ldim, bool fixed )
: Matrix(elem::LockedFlatMatrix<T>(height,width,ldim,buffer)) { }

template<typename T>
Matrix<T>::Matrix( elem::FlatMatrix<T>& fm, bool fixed )
: viewType_( fixed ? VIEW_FIXED: VIEW ), lfm_(fm)
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Matrix");
        Check( fm.layout );
    )
}
// Deprecated
template<typename T>
Matrix<T>::Matrix( Int height, Int width, T* buffer, Int ldim, bool fixed )
: Matrix(elem::FlatMatrix<T>(height,width,ldim,buffer)) { }

template<typename T>
Matrix<T>::Matrix( const Matrix<T>& A )
: viewType_(OWNER)
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Matrix( const Matrix& )"))
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct a Matrix with itself");
}

template<typename T>
Matrix<T>::Matrix( Matrix<T>&& A )
: viewType_(A.viewType_), memory_(std::move(A.memory_))
{
    lfm_.layout = A.Layout();
    lfm_.buffer = nullptr; 
    std::swap( lfm_.buffer, A.lfm_.buffer );
}

template<typename T>
Matrix<T>::~Matrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
const Matrix<T>&
Matrix<T>::operator=( const Matrix<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::operator=");
        if( Locked() )
            LogicError("Cannot assign to a locked view");
        if( viewType_ != OWNER && A.Dimensions() != Dimensions() )
            LogicError("Cannot assign to a view of different dimensions");
    )
    if( viewType_ == OWNER )
        Resize( A.Dimensions() );
    const Int height = Height();
    const Int width = Width();
    const Int ldim = LDim();
    const Int ldimOfA = A.LDim();
    const T* src = A.LockedBuffer();
    T* dst = Buffer();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
        MemCopy( &dst[j*ldim], &src[j*ldimOfA], height );
    return *this;
}

template<typename T>
Matrix<T>&
Matrix<T>::operator=( Matrix<T>&& A )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::operator=( Matrix&& )"))
    if( this == &A )
        LogicError("Tried to move to self");
    viewType_ = A.viewType_;
    lfm_.layout = A.Layout();
    std::swap( lfm_.buffer, A.lfm_.buffer );
    memory_.ShallowSwap( A.memory_ );
    return *this;
}

template<typename T>
void
Matrix<T>::Empty()
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Empty()");
        if ( FixedSize() )
            LogicError("Cannot empty a fixed-size matrix" );
    )
    Empty_();
}

template<typename T>
void
Matrix<T>::Resize( DimPair dims )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Resize");
        Check( dims );
        if( FixedSize() && dims != Dimensions() )
            LogicError("Cannot change the size of this matrix");
        if( Viewing() && (dims.m>Height() || dims.n>Width()) )
            LogicError("Cannot increase the size of this matrix");
    )
    Resize_( dims );
}
template<typename T>
void
Matrix<T>::Resize( Int height, Int width )
{ Resize(DimPair(height,width)); }

template<typename T>
void
Matrix<T>::Resize( Layout layout )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Resize(height,width,ldim)");
        Check( layout );
        if( FixedSize() && layout != Layout() )
            LogicError("Cannot change the layout of this matrix");
        if( Viewing() && 
            (layout.Height() > Height() || 
             layout.Width() > Width() ||
             layout.LDim() != LDim() )
            LogicError("Cannot increase the size of this matrix");
    )
    Resize_( layout );
}
template<typename T>
void
Matrix<T>::Resize( Int height, Int width, Int ldim )
{ Resize(Layout(height,width,ldim)); }

template<typename T>
void
Matrix<T>::Control( elem::FlatMatrix<T>& fm )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Control");
        if( FixedSize() )
            LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Control_( fm );
}
template<typename T>
void
Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim )
{ Control(elem::FlatMatrix<T>(height,width,ldim,buffer)); }

template<typename T>
void
Matrix<T>::Attach( elem::FlatMatrix<T>& fm )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Attach");
        if( FixedSize() )
            LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Attach_( fm );
}
template<typename T>
void
Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim )
{ Attach(elem::FlatMatrix<T>(height,width,ldim,buffer)); }

template<typename T>
void
Matrix<T>::LockedAttach( const elem::LockedFlatMatrix<T>& lfm )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::LockedAttach");
        if( FixedSize() )
            LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    LockedAttach_(lfm);
}
template<typename T>
void
Matrix<T>::LockedAttach( Int height, Int width, const T* buffer, Int ldim )
{ LockedAttach(elem::LockedFlatMatrix<T>(height,width,ldim,buffer)); }

// Basic queries
// =============

template<typename T>
Int Matrix<T>::Height() const { return lfm_.Height(); }
template<typename T>
Int Matrix<T>::Width() const { return lfm_.Width(); }
template<typename T>
Int Matrix<T>::LDim() const { return lfm_.LDim(); }
template<typename T>
Int Matrix<T>::MemorySize() const { return memory_.Size(); }

template<typename T>
Int Matrix<T>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(Height(),Width(),offset); }

template<typename T>
T*
Matrix<T>::Buffer()
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Buffer");
        if( Locked() )
            LogicError("Cannot return non-const buffer of locked Matrix");
    )
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return const_cast<T*>(lfm_.buffer);
}
template<typename T>
T*
Matrix<T>::Buffer( IndPair inds )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Buffer");
        if( Locked() )
            LogicError("Cannot return non-const buffer of locked Matrix");
    )
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return &const_cast<T*>(lfm_.buffer)[inds.i+inds.j*LDim()];
}
template<typename T>
T* Matrix<T>::Buffer( Int i, Int j )
{ return Buffer(IndPair(i,j)); }

template<typename T>
const T*
Matrix<T>::LockedBuffer() const
{ return lfm_.buffer; }
template<typename T>
const T*
Matrix<T>::LockedBuffer( IndPair inds ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::LockedBuffer"))
    return &lfm_.buffer[inds.i+inds.j*LDim()];
}
template<typename T>
const T*
Matrix<T>::LockedBuffer( Int i, Int j ) const
{ return LockedBuffer(IndPair(i,j)); }

template<typename T>
DimPair Matrix<T>::Dimensions() const { return lfm_.Dimensions(); }
template<typename T>
const elem::Layout& Matrix<T>::Layout() const { return lfm_.Layout(); }
template<typename T>
const elem::LockedFlatMatrix<T>& 
Matrix<T>::LockedFlatMatrix() const { return lfm_; }
template<typename T>
elem::FlatMatrix<T>
Matrix<T>::FlatMatrix() 
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::FlatMatrix");
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    elem::FlatMatrix<T> fm;    
    fm.layout = lfm_.layout;
    fm.buffer = const_cast<T*>(lfm_.buffer);
    return fm;
}

template<typename T>
bool Matrix<T>::Owner() const { return IsOwner( viewType_ ); }
template<typename T>
bool Matrix<T>::Viewing() const { return !IsOwner( viewType_ ); }
template<typename T>
bool Matrix<T>::Shrinkable() const { return IsShrinkable( viewType_ ); }
template<typename T>
bool Matrix<T>::FixedSize() const { return !IsShrinkable( viewType_ ); }
template<typename T>
bool Matrix<T>::Locked() const { return IsLocked( viewType_ ); }

// Single entry manipulation
// =========================

template<typename T>
T
Matrix<T>::Get( IndPair inds ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Get");
        Check( inds );
    )
    return Get_( inds );
}
template<typename T>
T Matrix<T>::Get( Int i, Int j ) const
{ return Get(IndPair(i,j)); }

template<typename T>
Base<T>
Matrix<T>::GetRealPart( IndPair inds ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetRealPart");
        Check( inds );
    )
    return elem::RealPart( Get_(inds) );
}
template<typename T>
Base<T>
Matrix<T>::GetRealPart( Int i, Int j ) const
{ return GetRealPart(IndPair(i,j)); }

template<typename T>
Base<T>
Matrix<T>::GetImagPart( IndPair inds ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetImagPart");
        Check( inds );
    )
    return elem::ImagPart( Get_(inds) );
}
template<typename T>
Base<T>
Matrix<T>::GetImagPart( Int i, Int j ) const
{ return GetImagPart(IndPair(i,j)); }

template<typename T>
void
Matrix<T>::Set( IndPair inds, T alpha ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Set");
        Check( inds );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    Set_(inds) = alpha;
}
template<typename T>
void
Matrix<T>::Set( Int i, Int j, T alpha ) 
{ Set(IndPair(i,j),alpha); }

template<typename T>
void 
Matrix<T>::SetRealPart( IndPair inds, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetRealPart");
        Check( inds );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    elem::SetRealPart( Set_(inds), alpha );
}
template<typename T>
void 
Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
{ SetRealPart(IndPair(i,j),alpha); }

template<typename T>
void 
Matrix<T>::SetImagPart( IndPair inds, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetImagPart");
        Check( inds );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    ComplainIfReal();
    elem::SetImagPart( Set_(inds), alpha );
}
template<typename T>
void 
Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
{ SetImagPart(IndPair(i,j),alpha); }

template<typename T>
void
Matrix<T>::Update( IndPair inds, T alpha ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Update");
        Check( inds );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    Set_(inds) += alpha;
}
template<typename T>
void
Matrix<T>::Update( Int i, Int j, T alpha ) 
{ Update(IndPair(i,j),alpha); }

template<typename T>
void 
Matrix<T>::UpdateRealPart( IndPair inds, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateRealPart");
        Check( inds );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    elem::UpdateRealPart( Set_(inds), alpha );
}
template<typename T>
void 
Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
{ UpdateRealPart(IndPair(i,j),alpha); }

template<typename T>
void 
Matrix<T>::UpdateImagPart( IndPair inds, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateImagPart");
        Check( inds );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    ComplainIfReal();
    elem::UpdateImagPart( Set_(inds), alpha );
}
template<typename T>
void 
Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
{ UpdateImagPart(IndPair(i,j),alpha); }

template<typename T>
void
Matrix<T>::MakeReal( IndPair inds )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::MakeReal"))
    Set( inds, GetRealPart(inds) );
}
template<typename T>
void
Matrix<T>::MakeReal( Int i, Int j )
{ MakeReal(IndPair(i,j)); }

template<typename T>
void
Matrix<T>::Conjugate( IndPair inds )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Conjugate"))
    Set( inds, elem::Conj(Get(inds)) );
}
template<typename T>
void
Matrix<T>::Conjugate( Int i, Int j )
{ Conjugate(IndPair(i,j)); }

// Arbitrary submatrix manipulation
// ================================

template<typename T>
void
Matrix<T>::Get
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Matrix<T>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Get"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.Resize( m, n );
    T* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    const T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            bufSub[i+j*ldSub] = buf[rowInd[i]+jSub*ld];
        }
    }
}

template<typename T>
Matrix<T> 
Matrix<T>::Get
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    Matrix<T> ASub;
    Get( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
void
Matrix<T>::GetRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::GetRealPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.Resize( m, n );
    Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    const T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            bufSub[i+j*ldSub] = RealPart(buf[rowInd[i]+jSub*ld]);
        }
    }
}

template<typename T>
void
Matrix<T>::GetImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::GetImagPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.Resize( m, n );
    Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    const T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            bufSub[i+j*ldSub] = ImagPart(buf[rowInd[i]+jSub*ld]);
        }
    }
}

template<typename T>
Matrix<Base<T>> 
Matrix<T>::GetRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    Matrix<Base<T>> ASub;
    GetRealPart( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
Matrix<Base<T>>
Matrix<T>::GetImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    Matrix<Base<T>> ASub;
    GetImagPart( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
void 
Matrix<T>::Set
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Set"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const T* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            buf[rowInd[i]+jSub*ld] = bufSub[i+j*ldSub];
        }
    }
}

template<typename T>
void
Matrix<T>::SetRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::SetRealPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            SetRealPart( buf[rowInd[i]+jSub*ld], bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::SetImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::SetImagPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            SetImagPart( buf[rowInd[i]+jSub*ld], bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::Update
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  T alpha, const Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Update"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const T* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            buf[rowInd[i]+jSub*ld] += alpha*bufSub[i+j*ldSub];
        }
    }
}

template<typename T>
void
Matrix<T>::UpdateRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::UpdateRealPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            UpdateRealPart( buf[rowInd[i]+jSub*ld], alpha*bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::UpdateImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::UpdateImagPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            UpdateImagPart( buf[rowInd[i]+jSub*ld], alpha*bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::MakeReal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::MakeReal"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            MakeReal( buf[rowInd[i]+jSub*ld] );
        }
    }
}

template<typename T>
void
Matrix<T>::Conjugate
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Conjugate"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(Check(IndPair(rowInd[i],colInd[j])))
            Conjugate( buf[rowInd[i]+jSub*ld] );
        }
    }
}

// Diagonal manipulation
// =====================

template<typename T>
void
Matrix<T>::GetDiagonal( Matrix<T>& d, Int offset ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetDiagonal");
        if( d.Locked() )
            LogicError("d must not be a locked view");
    )
    const Int diagLength = DiagonalLength(offset);
    d.Resize( diagLength, 1 );
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        d.Set_(k,0) = Get_(k+iOff,k+jOff);
}

template<typename T>
void
Matrix<T>::GetRealPartOfDiagonal( Matrix<Base<T>>& d, Int offset ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetRealPartOfDiagonal");
        if( d.Locked() )
            LogicError("d must not be a locked view");
    )
    const Int diagLength = DiagonalLength(offset);
    d.Resize( diagLength, 1 );
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        d.Set_(k,0) = elem::RealPart( Get_(k+iOff,k+jOff) );
}

template<typename T>
void
Matrix<T>::GetImagPartOfDiagonal( Matrix<Base<T>>& d, Int offset ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetImagPartOfDiagonal");
        if( d.Locked() )
            LogicError("d must not be a locked view");
    )
    const Int diagLength = DiagonalLength(offset);
    d.Resize( diagLength, 1 );
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        d.Set_( k, 0 ) = elem::ImagPart( Get_(k+iOff,k+jOff) );
}

template<typename T>
Matrix<T>
Matrix<T>::GetDiagonal( Int offset ) const
{ 
    Matrix<T> d;
    GetDiagonal( d, offset );
    return d;
}

template<typename T>
Matrix<Base<T>>
Matrix<T>::GetRealPartOfDiagonal( Int offset ) const
{ 
    Matrix<Base<T>> d;
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
Matrix<Base<T>>
Matrix<T>::GetImagPartOfDiagonal( Int offset ) const
{ 
    Matrix<Base<T>> d;
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
void
Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        Set_(k+iOff,k+jOff) = d.Get_(k,0);
}

template<typename T>
void
Matrix<T>::SetRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetRealPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        elem::SetRealPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

template<typename T>
void
Matrix<T>::SetImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetImagPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    ComplainIfReal();
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        elem::SetImagPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

template<typename T>
void
Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        Set_(k+iOff,k+jOff) += d.Get(k,0);
}

template<typename T>
void
Matrix<T>::UpdateRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateRealPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        elem::UpdateRealPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

template<typename T>
void
Matrix<T>::UpdateImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateImagPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    ComplainIfReal();
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        elem::UpdateImagPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

// Private auxilliary routines
// ===========================

template<typename T>
void
Matrix<T>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T>
void
Matrix<T>::Check( DimPair dims ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Check"))
    if( dims.m < 0 || dims.n < 0 )
        LogicError("Height and width must be non-negative");
}

template<typename T>
void
Matrix<T>::Check( Layout layout ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Check"))
    Check( dims );
    if( layout.LDim() < layout.Height() )
        LogicError("Leading dimension must be no less than height");
    if( layout.LDim() == 0 )
        LogicError("Leading dimension cannot be zero (for BLAS compatibility)");
}

template<typename T>
void
Matrix<T>::Check( IndPair inds ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Check"))
    if( inds.i < 0 || inds.j < 0 )
        LogicError("Indices must be non-negative");
    if( inds.i >= Height() || inds.j >= Width() )
        LogicError
        ("Out of bounds: (",inds.i,",",inds.j,") of ",
         Height()," x ",Width()," Matrix");
}

template<typename T>
void
Matrix<T>::ShallowSwap( Matrix<T>& A )
{
    memory_.ShallowSwap( A.memory_ );
    std::swap( viewType_, A.viewType_ );
    std::swap( lfm_, A.lfm_ );
}

template<typename T>
void
Matrix<T>::Empty_()
{
    memory_.Empty();
    lfm_.dims.m = 0;
    lfm_.dims.n = 0;
    lfm_.ldim = 1;
    lfm_.buffer = nullptr;
    viewType_ = (ViewType)( viewType_ & ~LOCKED_VIEW );
}

template<typename T>
void
Matrix<T>::Resize_( DimPair dims )
{
    bool reallocate = dims.m > LDim() || dims.n > Width();
    lfm_.layout.dims = dims;
    // Only change ldim when necessary. Simply 'shrink' view is possible
    if( reallocate )
    {
        lfm_.layout.ldim = Max( dims.m, 1 );
        memory_.Require( LDim()*Width() );
        lfm_.buffer = memory_.Buffer();
    }
}
template<typename T>
void
Matrix<T>::Resize_( Int height, Int width )
{ Resize_(DimPair(height,width)); }

template<typename T>
void
Matrix<T>::Resize_( Layout layout )
{
    bool reallocate = layout.Height() >  LDim()  || 
                      layout.Width()  >  Width() || 
                      layout.LDim()   != LDim();
    lfm_.layout.dims = layout.Dimensions();
    if( reallocate )
    {
        lfm_.layout.ldim = layout.LDim();
        memory_.Require(LDim()*Width());
        lfm_.buffer = memory_.Buffer();
    }
}
template<typename T>
void
Matrix<T>::Resize_( Int height, Int width, Int ldim )
{ Resize_(Layout(height,width,ldim)); }

template<typename T>
void
Matrix<T>::Control_( elem::FlatMatrix<T>& fm )
{
    memory_.Empty();
    lfm = fm;
    viewType_ = (ViewType)( viewType_ & ~LOCKED_VIEW );
}

template<typename T>
void
Matrix<T>::Attach_( elem::FlatMatrix<T>& fm )
{
    memory_.Empty();
    lfm_ = fm;
    viewType_ = (ViewType)( ( viewType_ & ~LOCKED_OWNER ) | VIEW );
}

template<typename T>
void
Matrix<T>::LockedAttach_( const elem::LockedFlatMatrix<T>& lfm )
{
    memory_.Empty();
    lfm_ = lfm;
    viewType_ = (ViewType)( viewType_ | VIEW );
}

template<typename T>
const T&
Matrix<T>::Get_( IndPair inds ) const
{ return lfm_.buffer[inds.i+inds.j*LDim()]; }

template<typename T>
T&
Matrix<T>::Set_( IndPair inds ) 
{
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return (const_cast<T*>(lfm_.buffer))[inds.i+inds.j*LDim()];
}

template class Matrix<Int>;
#ifndef DISABLE_FLOAT
template class Matrix<float>;
#endif // ifndef DISABLE_FLOAT
template class Matrix<double>;
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template class Matrix<Complex<float>>;
#endif // ifndef DISABLE_FLOAT
template class Matrix<Complex<double>>;
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
