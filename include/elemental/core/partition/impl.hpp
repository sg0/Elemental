/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_PARTITION_IMPL_HPP
#define ELEM_CORE_PARTITION_IMPL_HPP

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

// Partition from the bottom up
// ============================

template<typename T>
inline void
PartitionUp
( M& A, M& AT,
        M& AB, Int heightAB )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUp"))
    PartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionUp
( DM& A, DM& AT,
         DM& AB, Int heightAB )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUp"))
    PartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T>
inline void
LockedPartitionUp
( const M& A, M& AT,
              M& AB, Int heightAB )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUp"))
    LockedPartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, Int heightAB )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUp"))
    LockedPartitionDown( A, AT, AB, A.Height()-heightAB );
}

// Partition from the top down
// ===========================

template<typename T>
inline void
PartitionDown
( M& A, M& AT,
        M& AB, Int heightAT ) 
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    View( AT, A, IndPair(0,       0), DimPair(heightAT,A.Width()) );
    View( AB, A, IndPair(heightAT,0), DimPair(heightAB,A.Width()) );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionDown
( DM& A, DM& AT,
         DM& AB, Int heightAT )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    View( AT, A, IndPair(0,       0), DimPair(heightAT,A.Width()) );
    View( AB, A, IndPair(heightAT,0), DimPair(heightAB,A.Width()) );
}

template<typename T>
inline void
LockedPartitionDown
( const M& A, M& AT,
              M& AB, Int heightAT ) 
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    LockedView( AT, A, IndPair(0,       0), DimPair(heightAT,A.Width()) );
    LockedView( AB, A, IndPair(heightAT,0), DimPair(heightAB,A.Width()) );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, Int heightAT )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    LockedView( AT, A, IndPair(0,       0), DimPair(heightAT,A.Width()) );
    LockedView( AB, A, IndPair(heightAT,0), DimPair(heightAB,A.Width()) );
}

// Partition leftward starting from the right side
// ===============================================

template<typename T>
inline void
PartitionLeft( M& A, M& AL, M& AR, Int widthAR )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionLeft"))
    PartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionLeft( DM& A, DM& AL, DM& AR, Int widthAR )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionLeft"))
    PartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T>
inline void
LockedPartitionLeft( const M& A, M& AL, M& AR, Int widthAR )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionLeft"))
    LockedPartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionLeft( const DM& A, DM& AL, DM& AR, Int widthAR )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionLeft"))
    LockedPartitionRight( A, AL, AR, A.Width()-widthAR );
}

// Partition rightward starting from the left side
// ===============================================

template<typename T>
inline void
PartitionRight( M& A, M& AL, M& AR, Int widthAL )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    View( AL, A, IndPair(0,0      ), DimPair(A.Height(),widthAL) );
    View( AR, A, IndPair(0,widthAL), DimPair(A.Height(),widthAR) );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionRight( DM& A, DM& AL, DM& AR, Int widthAL )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    View( AL, A, IndPair(0,0      ), DimPair(A.Height(),widthAL) );
    View( AR, A, IndPair(0,widthAL), DimPair(A.Height(),widthAR) );
}

template<typename T>
inline void
LockedPartitionRight( const M& A, M& AL, M& AR, Int widthAL )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    LockedView( AL, A, IndPair(0,0      ), DimPair(A.Height(),widthAL) );
    LockedView( AR, A, IndPair(0,widthAL), DimPair(A.Height(),widthAR) );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionRight( const DM& A, DM& AL, DM& AR, Int widthAL )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    LockedView( AL, A, IndPair(0,0      ), DimPair(A.Height(),widthAL) );
    LockedView( AR, A, IndPair(0,widthAL), DimPair(A.Height(),widthAR) );
}

// Partition up the main diagonal
// ==============================

template<typename T>
inline void
PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpDiagonal"))
    PartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpDiagonal"))
    PartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
inline void
LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpDiagonal"))
    LockedPartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpDiagonal"))
    LockedPartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

// Partition up an arbitrary diagonal
// ==================================

template<typename T>
inline void
PartitionUpOffsetDiagonal
( Int offset,
  M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpOffsetDiagonal"))
    PartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionUpOffsetDiagonal
( Int offset,
  DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpOffsetDiagonal"))
    PartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T>
inline void
LockedPartitionUpOffsetDiagonal
( Int offset,
  const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpOffsetDiagonal"))
    LockedPartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionUpOffsetDiagonal
( Int offset,
  const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpOffsetDiagonal"))
    LockedPartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

// Partition down the main diagonal
// ================================

template<typename T>
inline void
PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownDiagonal"))
    PartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownDiagonal"))
    PartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
inline void
LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownDiagonal"))
    LockedPartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownDiagonal"))
    LockedPartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

// Partition down an arbitrary diagonal
// ====================================

template<typename T>
inline void
PartitionDownOffsetDiagonal
( Int offset,
  M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);
    
    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    View( ATL, A, IndPair(0,   0   ), DimPair(mCut,  nCut  ) );
    View( ATR, A, IndPair(0,   nCut), DimPair(mCut,  n-nCut) );
    View( ABL, A, IndPair(mCut,0   ), DimPair(m-mCut,nCut  ) );
    View( ABR, A, IndPair(mCut,nCut), DimPair(m-mCut,n-nCut) );
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionDownOffsetDiagonal
( Int offset,
  DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);

    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    View( ATL, A, IndPair(0,   0   ), DimPair(mCut,  nCut  ) );
    View( ATR, A, IndPair(0,   nCut), DimPair(mCut,  n-nCut) );
    View( ABL, A, IndPair(mCut,0   ), DimPair(m-mCut,nCut  ) );
    View( ABR, A, IndPair(mCut,nCut), DimPair(m-mCut,n-nCut) );
}

template<typename T>
inline void
LockedPartitionDownOffsetDiagonal
( Int offset,
  const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);
    
    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    LockedView( ATL, A, IndPair(0,   0   ), DimPair(mCut,  nCut  ) );
    LockedView( ATR, A, IndPair(0,   nCut), DimPair(mCut,  n-nCut) );
    LockedView( ABL, A, IndPair(mCut,0   ), DimPair(m-mCut,nCut  ) );
    LockedView( ABR, A, IndPair(mCut,nCut), DimPair(m-mCut,n-nCut) );
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionDownOffsetDiagonal
( Int offset,
  const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);
    
    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    LockedView( ATL, A, IndPair(0,   0   ), DimPair(mCut,  nCut  ) );
    LockedView( ATR, A, IndPair(0,   nCut), DimPair(mCut,  n-nCut) );
    LockedView( ABL, A, IndPair(mCut,0   ), DimPair(m-mCut,nCut  ) );
    LockedView( ABR, A, IndPair(mCut,nCut), DimPair(m-mCut,n-nCut) );
}

#undef DM
#undef M

} // namespace elem

#endif // ifndef ELEM_CORE_PARTITION_IMPL_HPP
