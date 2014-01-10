/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP
#define ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP

namespace elem {

template<typename T> 
class AbstractDistMatrix
{
public:
    typedef AbstractDistMatrix<T> type;

    // Constructors and destructors
    // ============================
#ifndef SWIG
    // Move constructor
    AbstractDistMatrix( type&& A );
#endif
    virtual ~AbstractDistMatrix();

    // Assignment and reconfiguration
    // ==============================
#ifndef SWIG
    // Move assignment
    type& operator=( type&& A );
#endif
    void Empty();
    void EmptyData();
    void SetGrid( const elem::Grid& grid );
    void Resize( DimPair dimPair );
    void Resize( Int height, Int width );
    void Resize( Layout layout );
    void Resize( Int height, Int width, Int ldim );
    void MakeConsistent();
    void Align( Int colAlign, Int rowAlign );
    void Align( IndPair aligns );
    void AlignCols( Int colAlign );
    void AlignRows( Int rowAlign );
    void FreeAlignments();
    void SetRoot( Int root );
    // Can be overridden
    // -----------------
    virtual void AlignWith( const elem::DistData& data );
    virtual void AlignColsWith( const elem::DistData& data );
    virtual void AlignRowsWith( const elem::DistData& data );

    // Basic queries
    // =============

    // Global information
    // ------------------
    Int Height() const;
    Int Width() const;
    DimPair Dimensions() const;
    Int DiagonalLength( Int offset=0 ) const;
    const elem::Grid& Grid() const;
    bool Viewing() const;
    bool Locked()  const;
    bool ColConstrained() const;
    bool RowConstrained() const;
    Int ColAlign() const;
    Int RowAlign() const;
    Int DistSize() const;
    Int CrossSize() const;
    Int RedundantSize() const;
    Int Root() const;
    // Rank in DistComm
    Int Owner( Int i, Int j ) const; 
    Int Owner( IndPair indPair ) const; 
    // Rank in ColComm
    Int RowOwner( Int i ) const; 
    // Rank in RowComm
    Int ColOwner( Int j ) const; 

    // Local information
    // -----------------
    Int LocalHeight() const;
    Int LocalWidth() const;
    Int LocalDimensions() const;
    Int LDim() const;
    size_t AllocatedMemory() const;
    elem::Matrix<T>& Matrix();
    T* Buffer( Int iLoc=0, Int jLoc=0 );
    T* Buffer( IndPair indPairLoc );
    const T* LockedBuffer( Int iLoc=0, Int jLoc=0 ) const;
    const T* LockedBuffer( IndPair indPairLoc ) const;
    const elem::Matrix<T>& LockedMatrix() const;
    bool Participating() const;
    Int ColShift() const;
    Int RowShift() const;
    Int ColRank() const;
    Int RowRank() const;
    Int DistRank() const;
    Int CrossRank() const;
    Int RedundantRank() const;
    bool IsLocal( Int i, Int j ) const;
    bool IsLocal( IndPair indPair ) const;
    bool IsLocalRow( Int i ) const; 
    bool IsLocalCol( Int j ) const;
    // Debug throws if row 'i' is not locally owned
    Int LocalRow( Int i ) const; 
    // Debug throws if column 'j' is not locally owned
    Int LocalCol( Int j ) const; 

    // MUST be overridden
    // ------------------
    virtual elem::DistData DistData() const = 0;
    virtual mpi::Comm DistComm() const = 0;
    virtual mpi::Comm CrossComm() const = 0;
    virtual mpi::Comm RedundantComm() const = 0;
    virtual mpi::Comm ColComm() const = 0;
    virtual mpi::Comm RowComm() const = 0;
    virtual Int ColStride() const = 0;
    virtual Int RowStride() const = 0;

    // Single-entry manipulation
    // =========================

    // Global manipulation
    // -------------------
    T Get( Int i, Int j    ) const;
    T Get( IndPair indPair ) const;
    BASE(T) GetRealPart( Int i, Int j    ) const;
    BASE(T) GetImagPart( Int i, Int j    ) const;
    BASE(T) GetRealPart( IndPair indPair ) const;
    BASE(T) GetImagPart( IndPair indPair ) const;
    void Set( Int i, Int j,    T alpha );
    void Set( IndPair indPair, T alpha );
    void SetRealPart( Int i, Int j,    BASE(T) alpha );
    void SetImagPart( Int i, Int j,    BASE(T) alpha );
    void SetRealPart( IndPair indPair, BASE(T) alpha );
    void SetImagPart( IndPair indPair, BASE(T) alpha );
    void Update( Int i, Int j,    T alpha );
    void Update( IndPair indPair, T alpha );
    void UpdateRealPart( Int i, Int j,    BASE(T) alpha );
    void UpdateImagPart( Int i, Int j,    BASE(T) alpha );
    void UpdateRealPart( IndPair indPair, BASE(T) alpha );
    void UpdateImagPart( IndPair indPair, BASE(T) alpha );
    void MakeReal( Int i, Int j    );
    void MakeReal( IndPair indPair );
    void Conjugate( Int i, Int j    );
    void Conjugate( IndPair indPair );

    // Local manipulation
    // ------------------
    T GetLocal( Int iLoc, Int jLoc ) const;
    T GetLocal( IndPair indPairLoc ) const;
    BASE(T) GetLocalRealPart( Int iLoc, Int jLoc ) const;
    BASE(T) GetLocalImagPart( Int iLoc, Int jLoc ) const;
    BASE(T) GetLocalRealPart( IndPair indPairLoc ) const;
    BASE(T) GetLocalImagPart( IndPair indPairLoc ) const;
    void SetLocal( Int iLoc, Int jLoc, T alpha );
    void SetLocal( IndPair indPairLoc, T alpha );
    void SetLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void SetLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void SetLocalRealPart( IndPair indPairLoc, BASE(T) alpha );
    void SetLocalImagPart( IndPair indPairLoc, BASE(T) alpha );
    void UpdateLocal( Int iLoc, Int jLoc, T alpha );
    void UpdateLocal( IndPair indPairLoc, T alpha );
    void UpdateLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void UpdateLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void UpdateLocalRealPart( IndPair indPairLoc, BASE(T) alpha );
    void UpdateLocalImagPart( IndPair indPairLoc, BASE(T) alpha );
    void MakeRealLocal( Int iLoc, Int jLoc );
    void MakeRealLocal( IndPair indPairLoc );
    void ConjugateLocal( Int iLoc, Int jLoc );
    void ConjugateLocal( IndPair indPairLoc );

    // Arbitrary submatrix manipulation
    // ================================

    // Global manipulation
    // -------------------
    template<Distribution U,Distribution V>
    void Get
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      DistMatrix<T,U,V>& ASub ) const;
    DistMatrix<T,STAR,STAR> Get
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    void GetRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      DistMatrix<BASE(T),STAR,STAR>& ASub ) const;
    void GetImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      DistMatrix<BASE(T),STAR,STAR>& ASub ) const;
    DistMatrix<BASE(T),STAR,STAR> GetRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    DistMatrix<BASE(T),STAR,STAR> GetImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    template<Distribution U,Distribution V>
    void Set
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      const DistMatrix<T,U,V>& ASub );
    template<Distribution U,Distribution V>
    void SetRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      const DistMatrix<BASE(T),U,V>& ASub );
    template<Distribution U,Distribution V>
    void SetImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      const DistMatrix<BASE(T),U,V>& ASub );
    template<Distribution U,Distribution V>
    void Update
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      T alpha, const DistMatrix<T,U,V>& ASub );
    template<Distribution U,Distribution V>
    void UpdateRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      BASE(T) alpha, const DistMatrix<BASE(T),U,V>& ASub );
    template<Distribution U,Distribution V>
    void UpdateImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      BASE(T) alpha, const DistMatrix<BASE(T),U,V>& ASub );
    void MakeReal
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd );
    void Conjugate
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd );

    // Local manipulation
    // ------------------
    // TODO

    // Diagonal manipulation
    // =====================
    // TODO: Standardize for all of DistMatrix

protected:
    ViewType viewType_;
    DimPair dims_;
    Memory<T> auxMemory_;
    elem::Matrix<T> matrix_;
    
    BoolPair constraints_;
    IndPair aligns_, shifts_;
    Int root_;
    const elem::Grid* grid_;

    // Build around a particular grid
    AbstractDistMatrix( const elem::Grid& g );

    // Exchange metadata with A
    virtual void ShallowSwap( type& A );

    void SetShifts();
    void SetColShift();
    void SetRowShift();
    void SetGrid();

    // Sanity checks
    void ComplainIfReal() const;
    void AssertNotLocked() const;
    void AssertNotStoring() const;
    void Check( IndPair indPair ) const;
    void Check( IndPair indPair, DimPair dimPair ) const;
    void CheckSame( const elem::Grid& grid ) const;
    void CheckSame( DimPair dimPair ) const;

    void AlignAndResize
    ( Int colAlign, Int rowAlign, Int height, Int width, bool force=false );
    void AlignColsAndResize
    ( Int colAlign, Int height, Int width, bool force=false );
    void AlignRowsAndResize
    ( Int rowAlign, Int height, Int width, bool force=false );

#ifndef SWIG
    template<typename S,Distribution U,Distribution V> 
    friend void View( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& B );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView( DistMatrix<S,U,V>& A, const DistMatrix<S,U,V>& B );
    template<typename S,Distribution U,Distribution V> 
    friend void View
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& B, 
      IndPair indPair, DimPair dimPair );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView
    ( DistMatrix<S,U,V>& A, const DistMatrix<S,U,V>& B,
      IndPair indPair, DimPair dimPair );
    template<typename S,Distribution U,Distribution V> 
    friend void View1x2
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& BL, DistMatrix<S,U,V>& BR );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView1x2
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BL, const DistMatrix<S,U,V>& BR );
    template<typename S,Distribution U,Distribution V> 
    friend void View2x1
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& BT, DistMatrix<S,U,V>& BB );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView2x1
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BT, const DistMatrix<S,U,V>& BB );
    template<typename S,Distribution U,Distribution V> 
    friend void View2x2
    ( DistMatrix<S,U,V>& A,
      DistMatrix<S,U,V>& BTL, DistMatrix<S,U,V>& BTR,
      DistMatrix<S,U,V>& BBL, DistMatrix<S,U,V>& BBR );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView2x2
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BTL, const DistMatrix<S,U,V>& BTR,
      const DistMatrix<S,U,V>& BBL, const DistMatrix<S,U,V>& BBR );

    template<typename S,Distribution U,Distribution V>
    friend class DistMatrix;
#endif // ifndef SWIG
};

DEBUG_ONLY(
    template<typename T>
    void AssertConforming1x2
    ( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR );

    template<typename T>
    void AssertConforming2x1
    ( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB );

    template<typename T>
    void AssertConforming2x2
    ( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR,
      const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR );
)

} // namespace elem

#endif // ifndef ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP
