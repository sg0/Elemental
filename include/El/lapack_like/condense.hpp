/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CONDENSE_HPP
#define EL_CONDENSE_HPP

namespace El {

// Bidiag
// ======

// Return the packed reduction to bidiagonal form
// ----------------------------------------------
template<typename F>
void Bidiag( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ );

template<typename F>
void Bidiag
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<F>& tP, AbstractDistMatrix<F>& tQ );

namespace bidiag {

// Overwrite A with B = Q^H A P and additionally return P and Q
template<typename F>
void Explicit
( Matrix<F>& A, 
  Matrix<F>& P, Matrix<F>& Q );
template<typename F>
void Explicit
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& Q );

// Only return the condensed bidiagonal matrix
template<typename F>
void ExplicitCondensed( Matrix<F>& A ); 
template<typename F>
void ExplicitCondensed( AbstractDistMatrix<F>& A );

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, 
        AbstractDistMatrix<F>& B );

template<typename F>
void ApplyP
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B );
template<typename F>
void ApplyP
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, 
        AbstractDistMatrix<F>& B );

} // namespace bidiag

// HermitianTridiag
// ================

namespace HermitianTridiagApproachNS {
enum HermitianTridiagApproach
{
    HERMITIAN_TRIDIAG_NORMAL, // Keep the current grid
    HERMITIAN_TRIDIAG_SQUARE, // Drop to a square process grid
    HERMITIAN_TRIDIAG_DEFAULT // Square grid algorithm only if already square
};
}
using namespace HermitianTridiagApproachNS;

template<typename F>
struct HermitianTridiagCtrl {
    HermitianTridiagApproach approach;
    GridOrder order;
    SymvCtrl<F> symvCtrl;

    HermitianTridiagCtrl<F>()
    : approach(HERMITIAN_TRIDIAG_SQUARE), order(ROW_MAJOR)
    { }
};

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t );
template<typename F>
void HermitianTridiag
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t,
  const HermitianTridiagCtrl<F>& ctrl=HermitianTridiagCtrl<F>() );

namespace herm_tridiag {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void ExplicitCondensed
( UpperOrLower uplo, AbstractDistMatrix<F>& A,
  const HermitianTridiagCtrl<F>& ctrl=HermitianTridiagCtrl<F>() );

template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, 
        AbstractDistMatrix<F>& B );

} // namespace herm_tridiag

// Hessenberg
// ==========
template<typename F>
void Hessenberg( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t );
template<typename F>
void Hessenberg
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t );

namespace hessenberg {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void ExplicitCondensed( UpperOrLower uplo, AbstractDistMatrix<F>& A );

template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, 
        AbstractDistMatrix<F>& B );

} // namespace hessenberg

} // namespace El

#endif // ifndef EL_CONDENSE_HPP
