/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SVT_PIVOTEDQR_HPP
#define EL_SVT_PIVOTEDQR_HPP

namespace El {
namespace svt {

// Preprocess with numSteps iterations of pivoted QR factorization

template<typename F>
Int PivotedQR( Matrix<F>& A, Base<F> tau, Int numSteps, bool relative )
{
    DEBUG_ONLY(
        CallStackEntry cse("svt::PivotedQR");
        if( numSteps > std::min(A.Height(),A.Width()) )
            LogicError("number of steps is too large");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<F> ACopy( A ), t;
    Matrix<Real> d;
    Matrix<Int> p;
    QRCtrl<Base<F>> qrCtrl;
    qrCtrl.boundRank = true;
    qrCtrl.maxRank = numSteps;
    QR( ACopy, t, d, p, qrCtrl );
    auto ACopyUpper = ACopy( IR(0,numSteps), IR(0,n) );

    Matrix<F> U( ACopyUpper ), V;
    Matrix<Real> s;
    MakeTrapezoidal( UPPER, U );

    SVDCtrl<Real> svdCtrl;
    svdCtrl.thresholded = true;
    svdCtrl.tol = tau;
    svdCtrl.relative = relative;
    SVD( U, s, V, svdCtrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    InversePermuteRows( V, p );
    Matrix<F> RThresh;
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, ACopy, t );
    DiagonalScale( RIGHT, NORMAL, d, ACopy );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

    return ZeroNorm( s );
}

template<typename F>
Int PivotedQR
( AbstractDistMatrix<F>& APre, Base<F> tau, Int numSteps, bool relative )
{
    DEBUG_ONLY(
        CallStackEntry cse("svt::PivotedQR");
        if( numSteps > std::min(APre.Height(),APre.Width()) )
            LogicError("number of steps is too large");
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy( A );
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Real,MD,STAR> d(g);
    DistMatrix<Int,VR,STAR> p(g);
    QRCtrl<Base<F>> qrCtrl;
    qrCtrl.boundRank = true;
    qrCtrl.maxRank = numSteps;
    QR( ACopy, t, d, p, qrCtrl );
    auto ACopyUpper = ACopy( IR(0,numSteps), IR(0,n) );

    DistMatrix<F> U( ACopyUpper ), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    MakeTrapezoidal( UPPER, U );

    SVDCtrl<Real> svdCtrl;
    svdCtrl.thresholded = true;
    svdCtrl.tol = tau;
    svdCtrl.relative = relative;
    SVD( U, s, V, svdCtrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    InversePermuteRows( V, p );
    DistMatrix<F> RThresh(g);
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, ACopy, t );
    DiagonalScale( RIGHT, NORMAL, d, ACopy );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_PIVOTEDQR_HPP
