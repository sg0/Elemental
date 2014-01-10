/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_MC_MR_DECL_HPP
#define ELEM_CORE_DISTMATRIX_MC_MR_DECL_HPP

namespace elem {

// Partial specialization to A[MC,MR].
//
// The columns of these matrices will be distributed among columns of the
// process grid, and the rows will be distributed among rows of the process
// grid.

template<typename T>
class DistMatrix<T,MC,MR> : public AbstractDistMatrix<T>
{
public:
    typedef DistMatrix<T,MC,MR> DM;
    typedef AbstractDistMatrix<T> ADM;

    // Constructors and destructors
    // ============================
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elem::Grid& g=DefaultGrid() );
    // Create a matrix with the specified dimensions
    DistMatrix( DimPair dims, const elem::Grid& g=DefaultGrid() );
    DistMatrix( Int height, Int width, const elem::Grid& g=DefaultGrid() );
    // Create a matrix with the specified dimensions and alignments
    DistMatrix( DimPair dims, IndPair aligns, const elem::Grid& g );
    DistMatrix
    ( Int height, Int width, Int colAlign, Int rowAlign, 
      const elem::Grid& g );
    // Create a matrix with the specified dimensions, alignments, 
    // and leading dimension
    DistMatrix( DimPair dims, IndPair aligns, Int ldim, const elem::Grid& g );
    DistMatrix
    ( Int height, Int width, Int colAlign, Int rowAlign, 
      Int ldim, const elem::Grid& g );
    // Configure around an existing buffer
    DistMatrix
    ( DimPair dims, IndPair aligns, const T* buffer, Int ldim, 
      const elem::Grid& g );
    DistMatrix
    ( Int height, Int width, Int colAlign, Int rowAlign,
      const T* buffer, Int ldim, const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlign, Int rowAlign,
      T* buffer, Int ldim, const elem::Grid& g );
    DistMatrix
    ( DimPair dims, IndPair aligns, T* buffer, Int ldim, 
      const elem::Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix( const DM& A );
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

#ifndef SWIG
    // Move constructor
    DistMatrix( DM&& A );
#endif

    ~DistMatrix();

    // Assignment and reconfiguration
    // ==============================
    const DM& operator=( const DistMatrix<T,MC,  MR  >& A );
    const DM& operator=( const DistMatrix<T,MC,  STAR>& A );
    const DM& operator=( const DistMatrix<T,STAR,MR  >& A );
    const DM& operator=( const DistMatrix<T,MD,  STAR>& A );
    const DM& operator=( const DistMatrix<T,STAR,MD  >& A );
    const DM& operator=( const DistMatrix<T,MR,  MC  >& A );
    const DM& operator=( const DistMatrix<T,MR,  STAR>& A );
    const DM& operator=( const DistMatrix<T,STAR,MC  >& A );
    const DM& operator=( const DistMatrix<T,VC,  STAR>& A );
    const DM& operator=( const DistMatrix<T,STAR,VC  >& A );
    const DM& operator=( const DistMatrix<T,VR,  STAR>& A );
    const DM& operator=( const DistMatrix<T,STAR,VR  >& A );
    const DM& operator=( const DistMatrix<T,STAR,STAR>& A );
    const DM& operator=( const DistMatrix<T,CIRC,CIRC>& A );
 #ifndef SWIG
    // Move assignment
    DM& operator=( DM&& A );
#endif
    virtual void AlignWith( const elem::DistData& data );
    virtual void AlignColsWith( const elem::DistData& data );
    virtual void AlignRowsWith( const elem::DistData& data );
    void Attach
    ( Int height, Int width, Int colAlign, Int rowAlign,
      T* buffer, Int ldim, const elem::Grid& grid );
    void Attach
    ( DimPair dims, IndPair aligns, T* buffer, Int ldim, 
      const elem::Grid& grid );
    void LockedAttach
    ( Int height, Int width, Int colAlign, Int rowAlign,
      const T* buffer, Int ldim, const elem::Grid& grid );      
    void LockedAttach
    ( DimPair dims, IndPair aligns, const T* buffer, Int ldim, 
      const elem::Grid& grid );      
    void Attach
    ( Matrix<T>& A, Int colAlign, Int rowAlign, const elem::Grid& grid );
    void Attach( Matrix<T>& A, IndPair aligns, const elem::Grid& grid );
    void LockedAttach
    ( const Matrix<T>& A, Int colAlign, Int rowAlign, const elem::Grid& grid );
    void LockedAttach
    ( const Matrix<T>& A, IndPair aligns, const elem::Grid& grid );
    // Equate/Update with scattered summation of A[MC,* ] across process rows
    void SumScatterFrom( const DistMatrix<T,MC,STAR>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,MC,STAR>& A );
    // Equate/Update with scattered summation of A[* ,MR] across process cols
    void SumScatterFrom( const DistMatrix<T,STAR,MR>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,MR>& A );
    // Equate/Update with scattered summation of A[* ,* ] across entire grid
    void SumScatterFrom( const DistMatrix<T,STAR,STAR>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR>& A );
    // Auxiliary routines for avoiding inefficient unpackings
    void AdjointFrom( const DistMatrix<T,STAR,MC>& A );
    void AdjointFrom( const DistMatrix<T,MR,STAR>& A );
    void AdjointSumScatterFrom( const DistMatrix<T,MR,STAR>& A );
    void AdjointSumScatterUpdate( T alpha, const DistMatrix<T,MR,STAR>& A );
    void TransposeFrom
    ( const DistMatrix<T,STAR,MC>& A, bool conjugate=false );
    void TransposeFrom
    ( const DistMatrix<T,MR,STAR>& A, bool conjugate=false );
    void TransposeSumScatterFrom
    ( const DistMatrix<T,MR,STAR>& A, bool conjugate=false );
    void TransposeSumScatterUpdate
    ( T alpha, const DistMatrix<T,MR,STAR>& A, bool conjugate=false );

    // Basic queries
    // =============
    virtual elem::DistData DistData() const;
    virtual mpi::Comm DistComm() const;
    virtual mpi::Comm CrossComm() const;
    virtual mpi::Comm RedundantComm() const;
    virtual mpi::Comm ColComm() const;
    virtual mpi::Comm RowComm() const;
    virtual Int RowStride() const;
    virtual Int ColStride() const;

    // Single-entry manipulation
    // =========================

    // Arbitrary submatrix manipulation
    // ================================

    // Diagonal manipulation
    // =====================
    void GetDiagonal( DistMatrix<T,MD,STAR>& d, Int offset=0 ) const;
    void GetDiagonal( DistMatrix<T,STAR,MD>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal
    ( DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal
    ( DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 ) const;
    DistMatrix<T,MD,STAR> GetDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),MD,STAR> GetRealPartOfDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),MD,STAR> GetImagPartOfDiagonal( Int offset=0 ) const;
    void SetDiagonal( const DistMatrix<T,MD,STAR>& d, Int offset=0 );
    void SetDiagonal( const DistMatrix<T,STAR,MD>& d, Int offset=0 );
    void SetRealPartOfDiagonal
    ( const DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 );
    void SetImagPartOfDiagonal
    ( const DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 );
    void SetRealPartOfDiagonal
    ( const DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 );
    void SetImagPartOfDiagonal
    ( const DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 );

private:
    template<typename S,class Function>
    void GetDiagonalHelper
    ( DistMatrix<S,MD,STAR>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void GetDiagonalHelper
    ( DistMatrix<S,STAR,MD>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const DistMatrix<S,MD,STAR>& d, Int offset, Function func );
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const DistMatrix<S,STAR,MD>& d, Int offset, Function func );

    void CopyFromDifferentGrid( const type& A );
#ifndef SWIG
    template<typename S,Distribution U,Distribution V>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef ELEM_CORE_DISTMATRIX_MC_MR_DECL_HPP
