/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

template<typename T>
void
elemental::basic::internal::HemmLL
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmLL");
#endif
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Width() )
        basic::internal::HemmLLA( alpha, A, B, beta, C );
    else
        basic::internal::HemmLLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::HemmLLA
( T alpha, const DistMatrix<T,MC,MR>& A, 
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmLLA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.Grid();

    DistMatrix<T,MC,MR> 
        BL(g), BR(g),
        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,MR>
        CL(g), CR(g),
        C0(g), C1(g), C2(g);

    DistMatrix<T,MC,STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Adj_STAR_MR(g);
    DistMatrix<T,MC,MR  > Z1(g);
    DistMatrix<T,MC,STAR> Z1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    basic::Scal( beta, C );
    LockedPartitionRight
    ( B, BL, BR, 0 );
    PartitionRight
    ( C, CL, CR, 0 );
    while( CL.Width() < C.Width() )
    {
        LockedRepartitionRight 
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );

        B1_MC_STAR.AlignWith( A );
        B1_VR_STAR.AlignWith( A );
        B1Adj_STAR_MR.AlignWith( A );
        Z1_MC_STAR.AlignWith( A );
        Z1_MR_STAR.AlignWith( A );
        Z1.AlignWith( C1 );
        Z1_MC_STAR.ResizeTo( C1.Height(), C1.Width() );
        Z1_MR_STAR.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        B1_MC_STAR = B1;
        B1_VR_STAR = B1_MC_STAR;
        B1Adj_STAR_MR.AdjointFrom( B1_VR_STAR );
        Z1_MC_STAR.SetToZero();
        Z1_MR_STAR.SetToZero();
        basic::internal::LocalSymmetricAccumulateLL
        ( ADJOINT, 
          alpha, A, B1_MC_STAR, B1Adj_STAR_MR, Z1_MC_STAR, Z1_MR_STAR );

        Z1_MR_MC.SumScatterFrom( Z1_MR_STAR );
        Z1 = Z1_MR_MC;
        Z1.SumScatterUpdate( (T)1, Z1_MC_STAR );
        basic::Axpy( (T)1, Z1, C1 );
        //--------------------------------------------------------------------//
        B1_MC_STAR.FreeAlignments();
        B1_VR_STAR.FreeAlignments();
        B1Adj_STAR_MR.FreeAlignments();
        Z1_MC_STAR.FreeAlignments();
        Z1_MR_STAR.FreeAlignments();
        Z1.FreeAlignments();

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::HemmLLC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmLLC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AColPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  ARowPan(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> 
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);

    DistMatrix<T,MC,MR> 
        CT(g),  C0(g),  CAbove(g),
        CB(g),  C1(g),  CBelow(g),
                C2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> AColPan_MC_STAR(g);
    DistMatrix<T,STAR,MC  > ARowPan_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Adj_MR_STAR(g);

    // Start the algorithm
    basic::Scal( beta, C );
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( CB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        ARowPan.LockedView1x2( A10, A11 );

        AColPan.LockedView2x1
        ( A11,
          A21 );

        CAbove.View2x1
        ( C0,
          C1 );

        CBelow.View2x1
        ( C1,
          C2 );

        AColPan_MC_STAR.AlignWith( CBelow );
        ARowPan_STAR_MC.AlignWith( CAbove );
        B1Adj_MR_STAR.AlignWith( C );
        //--------------------------------------------------------------------//
        AColPan_MC_STAR = AColPan;
        ARowPan_STAR_MC = ARowPan;
        AColPan_MC_STAR.MakeTrapezoidal( LEFT, LOWER );
        ARowPan_STAR_MC.MakeTrapezoidal( RIGHT, LOWER, -1 );

        B1Adj_MR_STAR.AdjointFrom( B1 );

        basic::internal::LocalGemm
        ( NORMAL, ADJOINT, 
          alpha, AColPan_MC_STAR, B1Adj_MR_STAR, (T)1, CBelow );

        basic::internal::LocalGemm
        ( ADJOINT, ADJOINT, 
          alpha, ARowPan_STAR_MC, B1Adj_MR_STAR, (T)1, CAbove );
        //--------------------------------------------------------------------//
        AColPan_MC_STAR.FreeAlignments();
        ARowPan_STAR_MC.FreeAlignments();
        B1Adj_MR_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );

        SlidePartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, T alpha,  
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,MC,  STAR>& B_MC_STAR,
  const DistMatrix<T,STAR,MR  >& BAdjOrTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalSymmetricAccumulateLL");
    if( A.Grid() != B_MC_STAR.Grid() || 
        B_MC_STAR.Grid() != BAdjOrTrans_STAR_MR.Grid() ||
        BAdjOrTrans_STAR_MR.Grid() != Z_MC_STAR.Grid() ||
        Z_MC_STAR.Grid() != Z_MR_STAR.Grid() )
        throw logic_error( "{A,B,Z} must be distributed over the same grid." );
    if( A.Height() != A.Width() || 
        A.Height() != B_MC_STAR.Height() ||
        A.Height() != BAdjOrTrans_STAR_MR.Width() ||
        A.Height() != Z_MC_STAR.Height() ||
        A.Height() != Z_MR_STAR.Height() ||
        B_MC_STAR.Width() != BAdjOrTrans_STAR_MR.Height() ||
        BAdjOrTrans_STAR_MR.Height() != Z_MC_STAR.Width() ||
        Z_MC_STAR.Width() != Z_MR_STAR.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalSymmetricAccumulateLL: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B[MC,* ] ~ " << B_MC_STAR.Height() << " x " 
                               << B_MC_STAR.Width() << "\n"
            << "  B^H[* ,MR] ~ " << BAdjOrTrans_STAR_MR.Height() << " x " 
                               << BAdjOrTrans_STAR_MR.Width() << "\n"
            << "  Z[MC,* ] ~ " << Z_MC_STAR.Height() << " x " 
                               << Z_MC_STAR.Width() << "\n"
            << "  Z[MR,* ] ` " << Z_MR_STAR.Height() << " x " 
                               << Z_MR_STAR.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( B_MC_STAR.ColAlignment() != A.ColAlignment() || 
        BAdjOrTrans_STAR_MR.RowAlignment() != A.RowAlignment() ||
        Z_MC_STAR.ColAlignment() != A.ColAlignment() ||
        Z_MR_STAR.ColAlignment() != A.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,MC,STAR> 
        BT_MC_STAR(g),  B0_MC_STAR(g),
        BB_MC_STAR(g),  B1_MC_STAR(g),
                        B2_MC_STAR(g);

    DistMatrix<T,STAR,MR>
        BLAdjOrTrans_STAR_MR(g), BRAdjOrTrans_STAR_MR(g),
        B0AdjOrTrans_STAR_MR(g), B1AdjOrTrans_STAR_MR(g), 
        B2AdjOrTrans_STAR_MR(g);

    DistMatrix<T,MC,STAR>
        ZT_MC_STAR(g),  Z0_MC_STAR(g),
        ZB_MC_STAR(g),  Z1_MC_STAR(g),
                        Z2_MC_STAR(g);

    DistMatrix<T,MR,STAR>
        ZT_MR_STAR(g),  Z0_MR_STAR(g),
        ZB_MR_STAR(g),  Z1_MR_STAR(g),
                        Z2_MR_STAR(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( B_MC_STAR, BT_MC_STAR,
                 BB_MC_STAR, 0 );
    LockedPartitionRight
    ( BAdjOrTrans_STAR_MR, BLAdjOrTrans_STAR_MR, BRAdjOrTrans_STAR_MR, 0 );
    PartitionDown
    ( Z_MC_STAR, ZT_MC_STAR,
                 ZB_MC_STAR, 0 );
    PartitionDown
    ( Z_MR_STAR, ZT_MR_STAR,
                 ZB_MR_STAR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( BT_MC_STAR,  B0_MC_STAR,
         /**********/ /**********/
                       B1_MC_STAR,
          BB_MC_STAR,  B2_MC_STAR );

        LockedRepartitionRight
        ( BLAdjOrTrans_STAR_MR, /**/ BRAdjOrTrans_STAR_MR,
          B0AdjOrTrans_STAR_MR, /**/ B1AdjOrTrans_STAR_MR, 
                                     B2AdjOrTrans_STAR_MR );

        RepartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
         /**********/ /**********/
                       Z1_MC_STAR,
          ZB_MC_STAR,  Z2_MC_STAR );

        RepartitionDown
        ( ZT_MR_STAR,  Z0_MR_STAR,
         /**********/ /**********/
                       Z1_MR_STAR,
          ZB_MR_STAR,  Z2_MR_STAR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( LEFT, LOWER );
        basic::internal::LocalGemm
        ( NORMAL, orientation, alpha, D11, B1AdjOrTrans_STAR_MR, 
          (T)1, Z1_MC_STAR );
        D11.MakeTrapezoidal( LEFT, LOWER, -1 );

        basic::internal::LocalGemm
        ( orientation, NORMAL, 
          alpha, D11, B1_MC_STAR, (T)1, Z1_MR_STAR );

        basic::internal::LocalGemm
        ( NORMAL, orientation, alpha, A21, B1AdjOrTrans_STAR_MR, 
          (T)1, Z2_MC_STAR );

        basic::internal::LocalGemm
        ( orientation, NORMAL, 
          alpha, A21, B2_MC_STAR, (T)1, Z1_MR_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( BT_MC_STAR,  B0_MC_STAR,
                       B1_MC_STAR,
         /**********/ /**********/
          BB_MC_STAR,  B2_MC_STAR );

        SlideLockedPartitionRight
        ( BLAdjOrTrans_STAR_MR,                       /**/ BRAdjOrTrans_STAR_MR,
          B0AdjOrTrans_STAR_MR, B1AdjOrTrans_STAR_MR, /**/ B2AdjOrTrans_STAR_MR 
        );

        SlidePartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
                       Z1_MC_STAR,
         /**********/ /**********/
          ZB_MC_STAR,  Z2_MC_STAR );

        SlidePartitionDown
        ( ZT_MR_STAR,  Z0_MR_STAR,
                       Z1_MR_STAR,
         /**********/ /**********/
          ZB_MR_STAR,  Z2_MR_STAR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::HemmLL
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, float alpha,  
  const DistMatrix<float,MC,  MR  >& A,
  const DistMatrix<float,MC,  STAR>& B_MC_STAR,
  const DistMatrix<float,STAR,MR  >& BAdjOrTrans_STAR_MR,
        DistMatrix<float,MC,  STAR>& Z_MC_STAR,
        DistMatrix<float,MR,  STAR>& Z_MR_STAR );

template void elemental::basic::internal::HemmLL
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, double alpha,  
  const DistMatrix<double,MC,  MR  >& A,
  const DistMatrix<double,MC,  STAR>& B_MC_STAR,
  const DistMatrix<double,STAR,MR  >& BAdjOrTrans_STAR_MR,
        DistMatrix<double,MC,  STAR>& Z_MC_STAR,
        DistMatrix<double,MR,  STAR>& Z_MR_STAR );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::HemmLL
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, scomplex alpha,  
  const DistMatrix<scomplex,MC,  MR  >& A,
  const DistMatrix<scomplex,MC,  STAR>& B_MC_STAR,
  const DistMatrix<scomplex,STAR,MR  >& BAdjOrTrans_STAR_MR,
        DistMatrix<scomplex,MC,  STAR>& Z_MC_STAR,
        DistMatrix<scomplex,MR,  STAR>& Z_MR_STAR );

template void elemental::basic::internal::HemmLL
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, dcomplex alpha,  
  const DistMatrix<dcomplex,MC,  MR  >& A,
  const DistMatrix<dcomplex,MC,  STAR>& B_MC_STAR,
  const DistMatrix<dcomplex,STAR,MR  >& BAdjOrTrans_STAR_MR,
        DistMatrix<dcomplex,MC,  STAR>& Z_MC_STAR,
        DistMatrix<dcomplex,MR,  STAR>& Z_MR_STAR );
#endif
