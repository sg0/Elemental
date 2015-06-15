/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_TNTT_HPP
#define EL_TRR2K_TNTT_HPP

namespace El {
namespace trr2k {

// E := alpha A' B + beta C' D' + gamma E
template<typename T>
void Trr2kTNTT
( UpperOrLower uplo,
  Orientation orientA,
  Orientation orientC, Orientation orientD,
  T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,  const AbstractDistMatrix<T>& CPre, const AbstractDistMatrix<T>& DPre,
  T gamma,       AbstractDistMatrix<T>& EPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("trr2k::Trr2kTNTT");
        if( EPre.Height() != EPre.Width()  || APre.Height() != CPre.Height() ||
            APre.Width()  != EPre.Height() || CPre.Width()  != EPre.Height() ||
            BPre.Width()  != EPre.Width()  || DPre.Height() != EPre.Width()  ||
            APre.Height() != BPre.Height() || CPre.Height() != DPre.Width() )
            LogicError("Nonconformal Trr2kNNTT");
    )
    const Int n = EPre.Height();
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = EPre.Grid();

    auto APtr = ReadProxy<T,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<T,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadProxy<T,MC,MR>( &CPre );      auto& C = *CPtr;
    auto DPtr = ReadProxy<T,MC,MR>( &DPre );      auto& D = *DPtr;
    auto EPtr = ReadWriteProxy<T,MC,MR>( &EPre ); auto& E = *EPtr;

    DistMatrix<T,STAR,MC  > A1_STAR_MC(g), C1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,VR,  STAR> D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > D1Trans_STAR_MR(g);

    A1_STAR_MC.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1_VR_STAR.AlignWith( E );
    D1Trans_STAR_MR.AlignWith( E );

    const Range<Int> outerInd( 0, n );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ind1,     outerInd );
        auto B1 = B( ind1,     outerInd );
        auto C1 = C( ind1,     outerInd );
        auto D1 = D( outerInd, ind1     );

        A1_STAR_MC = A1;
        C1_STAR_MC = C1;
        Transpose( B1, B1Trans_MR_STAR );
        D1_VR_STAR = D1;
        Transpose( D1_VR_STAR, D1Trans_STAR_MR, (orientD==ADJOINT) );
        LocalTrr2k
        ( uplo, orientA, TRANSPOSE, orientC, NORMAL,
          alpha, A1_STAR_MC, B1Trans_MR_STAR, 
          beta,  C1_STAR_MC, D1Trans_STAR_MR, gamma, E );
    }
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_TNTT_HPP