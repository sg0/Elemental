/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_MATRIX_HPP
#define ELEM_CORE_MATRIX_HPP

namespace elem {

template<typename T> // T must support addition and be MPI compatible
class Matrix
{
public:    
    typedef Matrix<T> type;
    typedef elem::FlatMatrix<T> FM;
    typedef elem::LockedFlatMatrix<T> LFM;

    // Constructors and destructors
    // ============================
    // Create a 0x0 matrix
    Matrix( bool fixed=false );
    // Create a matrix with the specified dimensions
    Matrix( DimPair dims, bool fixed=false );
    Matrix( Int height, Int width, bool fixed=false );
    // Create a matrix with the specified dimensions and leading dimension
    Matrix( const elem::Layout& layout, bool fixed=false );
    Matrix( Int height, Int width, Int ldim, bool fixed=false );
    // Construct around an existing buffer
    Matrix( FM& fm, bool fixed=false );
    Matrix( Int height, Int width, T* buffer, Int ldim, bool fixed=false );
    // Construct around an existing immutable buffer
    Matrix( const LFM& lfm, bool fixed=false );
    Matrix
    ( Int height, Int width, const T* buffer, Int ldim, bool fixed=false );
    // Copy constructor
    Matrix( const Matrix<T>& A );
#ifndef SWIG
    // Move constructor
    Matrix( Matrix<T>&& A );
#endif
    virtual ~Matrix();

    // Assignment and reconfiguration
    // ==============================
    const Matrix<T>& operator=( const Matrix<T>& A );
#ifndef SWIG
    // Move assignment
    Matrix<T>& operator=( Matrix<T>&& A );
#endif
    void Empty();
    void Resize( DimPair dims );
    void Resize( Int height, Int width );
    void Resize( const elem::Layout& layout );
    void Resize( Int height, Int width, Int ldim );
    void Attach( FM& fm );
    void Attach( Int height, Int width, T* buffer, Int ldim );
    void LockedAttach( const LFM& fm );
    void LockedAttach( Int height, Int width, const T* buffer, Int ldim );
    // Use this memory *as if it were not a view*, but do not take control of 
    // its deallocation. If Resize() forces reallocation, this buffer is 
    // released from control but not deleted.
    void Control( FM& fm );
    void Control( Int height, Int width, T* buffer, Int ldim );

    // Basic queries
    // =============
    Int Height() const;
    Int Width() const;
    Int LDim() const;
    Int MemorySize() const;
    Int DiagonalLength( Int offset=0 ) const;
    T* Buffer();
    T* Buffer( Int i, Int j );
    T* Buffer( IndPair inds );
    const T* LockedBuffer() const;
    const T* LockedBuffer( Int i, Int j ) const;
    const T* LockedBuffer( IndPair inds ) const;
    DimPair Dimensions() const;
    const elem::Layout& Layout() const;
    FM FlatMatrix();
    const LFM& LockedFlatMatrix() const;
    bool Owner()      const;
    bool Shrinkable() const;
    bool FixedSize()  const;
    bool Viewing()    const;
    bool Locked()     const;

    // Single-entry manipulation
    // =========================
    T Get( Int i, Int j ) const;
    T Get( IndPair inds ) const;
    BASE(T) GetRealPart( Int i, Int j ) const;
    BASE(T) GetImagPart( Int i, Int j ) const;
    BASE(T) GetRealPart( IndPair inds ) const;
    BASE(T) GetImagPart( IndPair inds ) const;
    void Set( Int i, Int j, T alpha );
    void Set( IndPair inds, T alpha );
    void SetRealPart( Int i, Int j, BASE(T) alpha );
    void SetImagPart( Int i, Int j, BASE(T) alpha );
    void SetRealPart( IndPair inds, BASE(T) alpha );
    void SetImagPart( IndPair inds, BASE(T) alpha );
    void Update( Int i, Int j, T alpha );
    void Update( IndPair inds, T alpha );
    void UpdateRealPart( Int i, Int j, BASE(T) alpha );
    void UpdateImagPart( Int i, Int j, BASE(T) alpha );
    void UpdateRealPart( IndPair inds, BASE(T) alpha );
    void UpdateImagPart( IndPair inds, BASE(T) alpha );
    void MakeReal( Int i, Int j );
    void MakeReal( IndPair inds );
    void Conjugate( Int i, Int j );
    void Conjugate( IndPair inds );

    // Arbitrary submatrix manipulation
    // ================================
    void Get
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      Matrix<T>& ASub ) const;
    Matrix<T> Get
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    void GetRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      Matrix<BASE(T)>& ASub ) const;
    void GetImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      Matrix<BASE(T)>& ASub ) const;
    Matrix<BASE(T)> GetRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    Matrix<BASE(T)> GetImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    void Set
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      const Matrix<T>& ASub );
    void SetRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      const Matrix<BASE(T)>& ASub );
    void SetImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
      const Matrix<BASE(T)>& ASub );
    void Update
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      T alpha, const Matrix<T>& ASub );
    void UpdateRealPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      BASE(T) alpha, const Matrix<BASE(T)>& ASub );
    void UpdateImagPart
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      BASE(T) alpha, const Matrix<BASE(T)>& ASub );
    void MakeReal
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd );
    void Conjugate
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd );

    // Diagonal manipulation
    // =====================
    void GetDiagonal( Matrix<T>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal( Matrix<BASE(T)>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal( Matrix<BASE(T)>& d, Int offset=0 ) const;
    Matrix<T> GetDiagonal( Int offset=0 ) const;
    Matrix<BASE(T)> GetRealPartOfDiagonal( Int offset=0 ) const;
    Matrix<BASE(T)> GetImagPartOfDiagonal( Int offset=0 ) const;
    void SetDiagonal( const Matrix<T>& d, Int offset=0 );
    void SetRealPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );
    void SetImagPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );
    void UpdateDiagonal( const Matrix<T>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );
    void UpdateImagPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );

private:
    ViewType viewType_;
    LFM lfm_;
    Memory<T> memory_;

    // Sanity checks
    void ComplainIfReal() const;
    void Check( DimPair dims ) const;
    void Check( const elem::Layout& layout ) const;
    void Check( IndPair inds ) const;
    void Check( IndPair inds, DimPair dims ) const;

    // Exchange metadata with A
    void ShallowSwap( Matrix<T>& A );

    const T& Get_( IndPair inds ) const;
    T& Set_( IndPair inds );
    void Empty_();
    void Resize_( Int height, Int width );
    void Resize_( DimPair dims );
    void Resize_( Int height, Int width, Int ldim );
    void Resize_( const elem::Layout& layout );
    void Attach_( FM& fm );
    void LockedAttach_( const LFM& lfm );
    void Control_( FM& fm );
    
#ifndef SWIG
    template <typename F> 
    friend class Matrix;
    template <typename F,Distribution U,Distribution V> 
    friend class DistMatrix;
    friend class AbstractDistMatrix<T>;

    friend void View<T>( Matrix<T>& A, Matrix<T>& B );
    friend void View<T>
    ( Matrix<T>& A, Matrix<T>& B, IndPair inds, DimPair dims );
    friend void View1x2<T>( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR );
    friend void View2x1<T>( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB );
    friend void View2x2<T>
    ( Matrix<T>& A, Matrix<T>& BTL, Matrix<T>& BTR,
                    Matrix<T>& BBL, Matrix<T>& BBR );

    friend void LockedView<T>( Matrix<T>& A, const Matrix<T>& B );
    friend void LockedView<T>
    ( Matrix<T>& A, const Matrix<T>& B, IndPair inds, DimPair dims );
    friend void LockedView1x2<T>
    ( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR );
    friend void LockedView2x1<T>
    ( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB );
    friend void LockedView2x2<T>
    ( Matrix<T>& A, const Matrix<T>& BTL, const Matrix<T>& BTR,
                    const Matrix<T>& BBL, const Matrix<T>& BBR );
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef ELEM_CORE_MATRIX_HPP
