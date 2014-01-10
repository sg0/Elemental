/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_VIEW_IMPL_HPP
#define ELEM_CORE_VIEW_IMPL_HPP

namespace elem {

// View a full matrix
// ==================

template<typename T>
inline void View( Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        if( IsLocked(B.viewType_) )
            LogicError("Cannot grab an unlocked view of a locked matrix");
    )
    A.memory_.Empty();
    A.lfm_ = B.lfm_;
    A.viewType_ = VIEW;
}

template<typename T>
inline Matrix<T> View( Matrix<T>& B )
{
    Matrix<T> A;
    View( A, B );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void View( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("View"))
    A.Empty();
    A.grid_ = B.grid_;
    A.dimPair_ = B.Dimensions();
    A.colAlign_ = B.ColAlign();
    A.rowAlign_ = B.RowAlign();
    HandleDiagPath( A, B );
    A.viewType_ = VIEW;
    A.SetShifts();
    if( A.Participating() )
        View( A.Matrix(), B.Matrix() );
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> View( DistMatrix<T,U,V>& B )
{
    DistMatrix<T,U,V> A(B.Grid());
    View( A, B );
    return A;
}

template<typename T>
inline void LockedView( Matrix<T>& A, const Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    A.memory_.Empty();
    A.lfm_ = B.lfm_;
    A.viewType_ = LOCKED_VIEW;
}

template<typename T>
inline Matrix<T> LockedView( const Matrix<T>& B )
{
    Matrix<T> A;
    LockedView( A, B );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void LockedView( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    A.Empty();
    A.grid_ = B.grid_;
    A.dimPair_ = B.Dimensions();
    A.colAlign_ = B.ColAlign();
    A.rowAlign_ = B.RowAlign();
    HandleDiagPath( A, B );
    A.viewType_ = LOCKED_VIEW;
    if( A.Participating() )
    {
        A.colShift_ = B.ColShift();
        A.rowShift_ = B.RowShift();
        LockedView( A.Matrix(), B.LockedMatrix() );
    }
    else 
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> LockedView( const DistMatrix<T,U,V>& B )
{
    DistMatrix<T,U,V> A(B.Grid());
    LockedView( A, B );
    return A;
}

// View a submatrix of the specified size starting at the specified coordinate
// ===========================================================================

template<typename T>
inline void View
( Matrix<T>& A, Matrix<T>& B, IndPair indPair, DimPair dimPair )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        const Int i = indPair.i;
        const Int j = indPair.j;
        const Int m = dimPair.m;
        const Int n = dimPair.n;
        if( i < 0 || j < 0 )
            LogicError("Indices must be non-negative");
        if( m < 0 || n < 0 )
            LogicError("Height and width must be non-negative");
        if( (i+m) > B.Height() || (j+n) > B.Width() )
            LogicError
            ("Trying to view outside of a Matrix: (",i,",",j,") up to (",
             i+m-1,",",j+n-1,") of ",B.Height()," x ",B.Width()," Matrix");
        if( IsLocked(B.viewType_) )
            LogicError("Cannot grab an unlocked view of a locked matrix");
    )
    A.memory_.Empty();
    A.lfm_.layout.dimPair = dimPair;
    A.lfm_.layout.ldim = B.LDim();
    A.lfm_.buffer = B.Buffer(indPair);
    A.viewType_ = VIEW;
}

template<typename T>
inline Matrix<T> View( Matrix<T>& B, IndPair indPair, DimPair dimPair )
{
    Matrix<T> A;
    View( A, B, indPair, dimPair );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void View
( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B, IndPair indPair, DimPair dimPair )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        B.Check( indPair, dimPair );
    )
    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();
    A.Empty();
    A.grid_ = &B.Grid();
    A.dimPair_ = dimPair;
    A.colAlign_ = (B.ColAlign()+indPair.i) % colStride;
    A.rowAlign_ = (B.RowAlign()+indPair.j) % rowStride;
    HandleDiagPath( A, B );
    A.viewType_ = VIEW;
    A.SetShifts();
    if( A.Participating() )
    {
        const Int localHeightBehind = Length(indPair.i,B.ColShift(),colStride);
        const Int localWidthBehind = Length(indPair.j,B.RowShift(),rowStride);
        const Int localHeight = Length(dimPair.m,A.ColShift(),colStride);
        const Int localWidth = Length(dimPair.n,A.RowShift(),rowStride);
        View
        ( A.Matrix(), B.Matrix(), 
          IndPair(localHeightBehind,localWidthBehind), 
          DimPair(localHeight,localWidth) );
    }
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> View
( DistMatrix<T,U,V>& B, IndPair indPair, DimPair dimPair )
{
    DistMatrix<T,U,V> A(B.Grid());
    View( A, B, indPair, dimPair );
    return A;
}

template<typename T>
inline void LockedView
( Matrix<T>& A, const Matrix<T>& B, IndPair indPair, DimPair dimPair )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView");
        const Int i = indPair.i;
        const Int j = indPair.j;
        const Int m = dimPair.m;
        const Int n = dimPair.n;
        if( i < 0 || j < 0 )
            LogicError("Indices must be non-negative");
        if( m < 0 || n < 0 )
            LogicError("Height and width must be non-negative");
        if( (i+m) > B.Height() || (j+n) > B.Width() )
            LogicError
            ("Trying to view outside of a Matrix: (",i,",",j,") up to (",
             i+m-1,",",j+n-1,") of ",B.Height()," x ",B.Width()," Matrix");
    )
    A.memory_.Empty();
    A.lfm_.layout = B.Layout();
    A.lfm_.buffer = B.LockedBuffer(indPair);
    A.viewType_ = LOCKED_VIEW;
}

template<typename T>
inline Matrix<T> LockedView
( const Matrix<T>& B, IndPair indPair, DimPair dimPair )
{
    Matrix<T> A;
    LockedView( A, B, indPair, dimPair );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void LockedView
( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B,
  IndPair indPair, DimPair dimPair )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView");
        B.Check( indPair, dimPair );
    )
    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();
    const Int colRank = B.ColRank();
    const Int rowRank = B.RowRank();
    A.Empty();
    A.grid_ = &B.Grid();
    A.dimPair_ = dimPair;
    A.colAlign_ = (B.ColAlign()+indPair.i) % colStride;
    A.rowAlign_ = (B.RowAlign()+indPair.j) % rowStride;
    HandleDiagPath( A, B );
    A.viewType_ = LOCKED_VIEW;
    A.SetShifts();
    if( A.Participating() )
    {
        const Int localHeightBehind = Length(indPair.i,B.ColShift(),colStride);
        const Int localWidthBehind = Length(indPair.j,B.RowShift(),rowStride);
        const Int localHeight = Length(dimPair.m,A.ColShift(),colStride);
        const Int localWidth = Length(dimPair.n,A.RowShift(),rowStride);
        LockedView
        ( A.Matrix(), B.LockedMatrix(), 
          IndPair(localHeightBehind,localWidthBehind), 
          DimPair(localHeight,localWidth) );
    }
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> LockedView
( const DistMatrix<T,U,V>& B, IndPair indPair, DimPair dimPair )
{
    DistMatrix<T,U,V> A(B.Grid());
    LockedView( A, B, indPair, dimPair );
    return A;
}

// View a submatrix [iBeg:iEnd),[jBeg:jEnd)
// ========================================

template<typename T>
void View( Matrix<T>& A, Matrix<T>& B, IndPair begPair, IndPair endPair )
{ 
    DEBUG_ONLY(CallStackEntry cse("View"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    View( A, B, begPair, dimPair );
}

template<typename T>
Matrix<T> View( Matrix<T>& B, IndPair begPair, IndPair endPair )
{
    DEBUG_ONLY(CallStackEntry cse("View"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    return View( B, begPair, dimPair ); 
}

template<typename T,Distribution U,Distribution V>
void View
( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B, IndPair begPair, IndPair endPair )
{
    DEBUG_ONLY(CallStackEntry cse("View"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    View( A, B, begPair, dimPair ); 
}

template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> View( DistMatrix<T,U,V>& B, IndPair begPair, IndPair endPair )
{
    DEBUG_ONLY(CallStackEntry cse("View"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    return View( B, begPair, dimPair ); 
} 

template<typename T>
void LockedView
( Matrix<T>& A, const Matrix<T>& B, IndPair begPair, IndPair endPair )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    LockedView( A, B, begPair, dimPair );
}

template<typename T>
Matrix<T> LockedView( const Matrix<T>& B, IndPair begPair, IndPair endPair )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    return LockedView( B, begPair, dimPair );
}

template<typename T,Distribution U,Distribution V>
void LockedViewRange
( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B,
  IndPair begPair, IndPair endPair )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    LockedView( A, B, begPair, dimPair );
}

template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> LockedView
( const DistMatrix<T,U,V>& B, IndPair begPair, IndPair endPair )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    DimPair dimPair;
    dimPair.m = endPair.i-begPair.i;
    dimPair.n = endPair.j-begPair.j;
    return LockedView( B, begPair, dimPair );
}

// View two horizontally-connected matrices
// ========================================

template<typename T>
inline void View1x2( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View1x2");
        if( IsLocked(BL.viewType_) || IsLocked(BR.viewType_) )
            LogicError("Cannot grab an unlocked view of a locked matrix");
        if( BL.Height() != BR.Height() )
            LogicError("1x2 must have consistent height to combine");
        if( BL.LDim() != BR.LDim() )
            LogicError("1x2 must have consistent ldims to combine");
        if( BR.Buffer() != (BL.Buffer()+BL.LDim()*BL.Width()) )
            LogicError("1x2 must have contiguous memory");
    )
    A.memory_.Empty();
    A.lfm_.layout.dimPair.m = BL.Height();
    A.lfm_.layout.dimPair.n = BL.Width() + BR.Width();
    A.lfm_.layout.ldim = BL.LDim();
    A.lfm_.buffer = BL.Buffer();
    A.viewType_ = VIEW;
}

template<typename T>
inline Matrix<T> View1x2( Matrix<T>& BL, Matrix<T>& BR )
{
    Matrix<T> A;
    View1x2( A, BL, BR );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void View1x2
( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View1x2");
        AssertConforming1x2( BL, BR );
        BL.CheckSame( BR.Grid() );
    )
    A.Empty();
    A.grid_ = BL.grid_;
    A.dimPair_.m = BL.Height();
    A.dimPair_.n = BL.Width() + BR.Width();
    A.colAlign_ = BL.ColAlign();
    A.rowAlign_ = BL.RowAlign();
    HandleDiagPath( A, BL );
    A.viewType_ = VIEW;
    A.SetShifts();
    if( A.Participating() )
        View1x2( A.Matrix(), BL.Matrix(), BR.Matrix() );
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> View1x2( DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR )
{
    DistMatrix<T,U,V> A(BL.Grid());
    View1x2( A, BL, BR );
    return A;
}

template<typename T>
inline void LockedView1x2
( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView1x2");
        if( BL.Height() != BR.Height() )
            LogicError("1x2 must have consistent height to combine");
        if( BL.LDim() != BR.LDim() )
            LogicError("1x2 must have consistent ldims to combine");
        if( BR.LockedBuffer() != (BL.LockedBuffer()+BL.LDim()*BL.Width()) )
            LogicError("1x2 must have contiguous memory");
    )
    A.memory_.Empty();
    A.lfm_.layout.dimPair.m = BL.Height();
    A.lfm_.layout.dimPair.n = BL.Width() + BR.Width();
    A.lfm_.layout.ldim = BL.LDim();
    A.lfm_.buffer = BL.LockedBuffer();
    A.viewType_ = LOCKED_VIEW;
}

template<typename T>
inline Matrix<T> LockedView1x2( const Matrix<T>& BL, const Matrix<T>& BR )
{
    Matrix<T> A;
    LockedView1x2( A, BL, BR );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void LockedView1x2
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView1x2");
        AssertConforming1x2( BL, BR );
        BL.CheckSame( BR.Grid() );
    )
    A.Empty();
    A.grid_ = BL.grid_;
    A.dimPair_.m = BL.Height();
    A.dimPair_.n = BL.Width() + BR.Width();
    A.colAlign_ = BL.ColAlign();
    A.rowAlign_ = BL.RowAlign();
    HandleDiagPath( A, BL );
    A.viewType_ = LOCKED_VIEW;
    A.SetShifts();
    if( A.Participating() )
        LockedView1x2( A.Matrix(), BL.LockedMatrix(), BR.LockedMatrix() );
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> LockedView1x2
( const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR )
{
    DistMatrix<T,U,V> A(BL.Grid());
    LockedView1x2( A, BL, BR );
    return A;
}

// View two vertically-connected matrices
// ======================================

template<typename T>
inline void View2x1( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x1");
        if( IsLocked(BT.viewType_) || IsLocked(BB.viewType_) )
            LogicError("Cannot grab an unlocked view of a locked matrix");
        if( BT.Width() != BB.Width() )
            LogicError("2x1 must have consistent width to combine");
        if( BT.LDim() != BB.LDim() )
            LogicError("2x1 must have consistent ldim to combine");
        if( BB.Buffer() != (BT.Buffer() + BT.Height()) )
            LogicError("2x1 must have contiguous memory");
    )
    A.memory_.Empty();
    A.lfm_.layout.dimPair.m = BT.Height() + BB.Height();
    A.lfm_.layout.dimPair.n = BT.Width();
    A.lfm_.layout.ldim = BT.LDim();
    A.lfm_.buffer = BT.Buffer();
    A.viewType_ = VIEW;
}

template<typename T>
inline Matrix<T> View2x1( Matrix<T>& BT, Matrix<T>& BB )
{
    Matrix<T> A;
    View2x1( A, BT, BB );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void View2x1
( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x1");
        AssertConforming2x1( BT, BB );
        BT.CheckSame( BB.Grid() );
    )
    A.Empty();
    A.grid_ = BT.grid_;
    A.dimPair_.m = BT.Height() + BB.Height();
    A.dimPair_.n = BT.Width();
    A.colAlign_ = BT.ColAlign();
    A.rowAlign_ = BT.RowAlign();
    HandleDiagPath( A, BT );
    A.viewType_ = LOCKED_VIEW;
    A.SetShifts();
    if( A.Participating() )
        View2x1( A.Matrix(), BT.Matrix(), BB.Matrix() );
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> View2x1( DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB )
{
    DistMatrix<T,U,V> A(BT.Grid());
    View2x1( A, BT, BB );
    return A;
}

template<typename T>
inline void LockedView2x1
( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x1");
        if( BT.Width() != BB.Width() )
            LogicError("2x1 must have consistent width to combine");
        if( BT.LDim() != BB.LDim() )
            LogicError("2x1 must have consistent ldim to combine");
        if( BB.LockedBuffer() != (BT.LockedBuffer() + BT.Height()) )
            LogicError("2x1 must have contiguous memory");
    )
    A.memory_.Empty();
    A.lfm_.layout.dimPair.m = BT.Height() + BB.Height();
    A.lfm_.layout.dimPair.n = BT.Width();
    A.lfm_.layout.ldim = BT.LDim();
    A.lfm_.buffer = BT.LockedBuffer();
    A.viewType_ = LOCKED_VIEW;
}

template<typename T>
inline Matrix<T> LockedView2x1
( const Matrix<T>& BT, const Matrix<T>& BB )
{
    Matrix<T> A;
    LockedView2x1( A, BT, BB );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void LockedView2x1
(       DistMatrix<T,U,V>& A,
  const DistMatrix<T,U,V>& BT,
  const DistMatrix<T,U,V>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x1");
        AssertConforming2x1( BT, BB );
        BT.CheckSame( BB.Grid() );
    )
    A.Empty();
    A.grid_ = BT.grid_;
    A.dimPair_.m = BT.Height() + BB.Height();
    A.dimPair_.n = BT.Width();
    A.colAlign_ = BT.ColAlign();
    A.rowAlign_ = BT.RowAlign();
    HandleDiagPath( A, BT );
    A.viewType_ = LOCKED_VIEW;
    A.SetShifts();
    if( A.Participating() )
        LockedView2x1( A.Matrix(), BT.LockedMatrix(), BB.LockedMatrix() );
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> LockedView2x1
( const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB )
{
    DistMatrix<T,U,V> A(BT.Grid());
    LockedView2x1( A, BT, BB );
    return A;
}

// View a two-by-two set of connected matrices
// ===========================================

template<typename T>
inline void View2x2
( Matrix<T>& A,
  Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x2");
        if( IsLocked(BTL.viewType_) || IsLocked(BTR.viewType_) ||
            IsLocked(BBL.viewType_) || IsLocked(BBR.viewType_) )
            LogicError("Cannot grab an unlocked view of a locked matrix");
        if( BTL.Width() != BBL.Width()   ||
            BTR.Width() != BBR.Width()   ||
            BTL.Height() != BTR.Height() ||
            BBL.Height() != BBR.Height()   )
            LogicError("2x2 must conform to combine");
        if( BTL.LDim() != BTR.LDim() ||
            BTR.LDim() != BBL.LDim() ||
            BBL.LDim() != BBR.LDim()   )
            LogicError("2x2 must have consistent ldims to combine");
        if( BBL.Buffer() != (BTL.Buffer() + BTL.Height()) ||
            BBR.Buffer() != (BTR.Buffer() + BTR.Height()) ||
            BTR.Buffer() != (BTL.Buffer() + BTL.LDim()*BTL.Width()) )
            LogicError("2x2 must have contiguous memory");
    )
    A.memory_.Empty();
    A.lfm_.layout.dimPair.m = BTL.Height() + BBL.Height();
    A.lfm_.layout.dimPair.n = BTL.Width() + BTR.Width();
    A.lfm_.layout.ldim = BTL.LDim();
    A.lfm_.buffer = BTL.Buffer();
    A.viewType_ = VIEW;
}

template<typename T>
inline Matrix<T> View2x2
( Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR )
{
    Matrix<T> A;
    View2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void View2x2
( DistMatrix<T,U,V>& A,
  DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR,
  DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x2");
        AssertConforming2x2( BTL, BTR, BBL, BBR );
        BTL.CheckSame( BTR.Grid() );
        BTL.CheckSame( BBL.Grid() );
        BTL.CheckSame( BBR.Grid() );
    )
    A.Empty();
    A.grid_ = BTL.grid_;
    A.dimPair_.m = BTL.Height() + BBL.Height();
    A.dimPair_.n = BTL.Width() + BTR.Width();
    A.colAlign_ = BTL.ColAlign();
    A.rowAlign_ = BTL.RowAlign();
    HandleDiagPath( A, BTL );
    A.viewType_ = VIEW;
    A.SetShifts();
    if( A.Participating() )
        View2x2
        ( A.Matrix(), BTL.Matrix(), BTR.Matrix(),
                      BBL.Matrix(), BBR.Matrix() ); 
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> View2x2
( DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR,
  DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR )
{
    DistMatrix<T,U,V> A(BTL.Grid());
    View2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T>
inline void LockedView2x2
(       Matrix<T>& A,
  const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x2");
        if( BTL.Width() != BBL.Width()   ||
            BTR.Width() != BBR.Width()   ||
            BTL.Height() != BTR.Height() ||
            BBL.Height() != BBR.Height()   )
            LogicError("2x2 must conform to combine");
        if( BTL.LDim() != BTR.LDim() ||
            BTR.LDim() != BBL.LDim() ||
            BBL.LDim() != BBR.LDim()   )
            LogicError("2x2 must have consistent ldims to combine");
        if( BBL.LockedBuffer() != (BTL.LockedBuffer()+BTL.Height()) ||
            BBR.LockedBuffer() != (BTR.LockedBuffer()+BTR.Height()) ||
            BTR.LockedBuffer() != (BTL.LockedBuffer()+BTL.LDim()*BTL.Width()) )
            LogicError("2x2 must have contiguous memory");
    )
    A.memory_.Empty();
    A.lfm_.layout.dimPair.m = BTL.Height() + BBL.Height();
    A.lfm_.layout.dimPair.n = BTL.Width() + BTR.Width();
    A.lfm_.layout.ldim = BTL.LDim();
    A.lfm_.buffer = BTL.LockedBuffer();
    A.viewType_ = LOCKED_VIEW;
}

template<typename T>
inline Matrix<T> LockedView2x2
( const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR )
{
    Matrix<T> A;
    LockedView2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void LockedView2x2
(       DistMatrix<T,U,V>& A,
  const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR,
  const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x2");
        AssertConforming2x2( BTL, BTR, BBL, BBR );
        BTL.CheckSame( BTR.Grid() );
        BTL.CheckSame( BBL.Grid() );
        BTL.CheckSame( BBR.Grid() );
    )
    A.Empty();
    A.grid_ = BTL.grid_;
    A.dimPair_.m = BTL.Height() + BBL.Height();
    A.dimPair_.n = BTL.Width() + BTR.Width();
    A.colAlign_ = BTL.ColAlign();
    A.rowAlign_ = BTL.RowAlign();
    HandleDiagPath( A, BTL );
    A.viewType_ = LOCKED_VIEW;
    A.SetShifts();
    if( A.Participating() )
        LockedView2x2
        ( A.Matrix(), BTL.LockedMatrix(), BTR.LockedMatrix(),
                      BBL.LockedMatrix(), BBR.LockedMatrix() );
}

template<typename T,Distribution U,Distribution V>
inline DistMatrix<T,U,V> LockedView2x2
( const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR,
  const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR )
{
    DistMatrix<T,U,V> A(BTL.Grid());
    LockedView2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

// Utilities for handling the extra information needed for [MD,* ] and [* ,MD]
// ===========================================================================

template<typename T,Distribution U,Distribution V>
inline void HandleDiagPath
( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B ) { }

template<typename T>
inline void HandleDiagPath
( DistMatrix<T,MD,STAR>& A, const DistMatrix<T,MD,STAR>& B )
{ A.root_ = B.root_; } 

template<typename T>
inline void HandleDiagPath
( DistMatrix<T,STAR,MD>& A, const DistMatrix<T,STAR,MD>& B )
{ A.root_ = B.root_; } 

} // namespace elem

#endif // ifndef ELEM_CORE_VIEW_IMPL_HPP
