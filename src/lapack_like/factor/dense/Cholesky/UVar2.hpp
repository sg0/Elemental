/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CHOLESKY_UVAR2_HPP
#define EL_CHOLESKY_UVAR2_HPP

#include "./UVar3.hpp"

// TODO: Reverse variants

namespace El {
namespace cholesky {

template<typename F> 
inline void
UVar2( Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar2");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        Herk( UPPER, ADJOINT, F(-1), A01, F(1), A11 );
        cholesky::UVar3Unb( A11 );
        Gemm( ADJOINT, NORMAL, F(-1), A02, A01, F(1), A12 );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A12 );
    }
}

template<typename F> 
inline void
UVar2( AbstractDistMatrix<F>& APre )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar2");
        if( APre.Height() != APre.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = APre.Grid();

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    DistMatrix<F,MC,  STAR> A01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,MR,  STAR> X11Adj_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X11Adj_MR_MC(g);
    DistMatrix<F,MR,  STAR> X12Adj_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X12Adj_MR_MC(g);
    DistMatrix<F> X11(g), X12(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        A01_MC_STAR.AlignWith( A01 );
        A01_MC_STAR = A01;
        X11Adj_MR_STAR.AlignWith( A01 );
        LocalGemm( ADJOINT, NORMAL, F(1), A01, A01_MC_STAR, X11Adj_MR_STAR );
        X11Adj_MR_MC.AlignWith( A11 );
        Contract( X11Adj_MR_STAR, X11Adj_MR_MC );
        X11.AlignWith( A11 );
        Adjoint( X11Adj_MR_MC, X11 );
        Axpy( F(-1), X11, A11 );

        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        X12Adj_MR_STAR.AlignWith( A02 );
        LocalGemm( ADJOINT, NORMAL, F(1), A02, A01_MC_STAR, X12Adj_MR_STAR );
        X12Adj_MR_MC.AlignWith( A12 );
        Contract( X12Adj_MR_STAR, X12Adj_MR_MC );
        X12.AlignWith( A12 );
        Adjoint( X12Adj_MR_MC, X12 );
        Axpy( F(-1), X12, A12 );

        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_UVAR2_HPP
