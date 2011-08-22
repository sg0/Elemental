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
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;

template<typename F> // represents a real or complex number
void
elemental::advanced::LU
( DistMatrix<F,MC,MR>& A, DistMatrix<int,VC,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::LU");
    if( A.Grid() != p.Grid() )
        throw logic_error( "A and p must be distributed over the same grid." );
    if( A.Height() != p.Height() ) 
        throw logic_error( "A and p must be the same height." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AB(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  
                         A20(g), A21(g), A22(g);

    DistMatrix<int,VC,STAR>
        pT(g),  p0(g), 
        pB(g),  p1(g),
                p2(g);

    // Temporary distributions
    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);
    DistMatrix<int,STAR,STAR> p1_STAR_STAR(g);

    // Pivot composition
    vector<int> image, preimage;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( p, pT,
         pB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( pT,  p0,
         /**/ /**/
               p1,
          pB,  p2 );

        AB.View1x2( ABL, ABR );

        int pivotOffset = A01.Height();
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_MR.AlignWith( A22 );
        A21_MC_STAR.AlignWith( A22 );
        A11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        p1_STAR_STAR.ResizeTo( p1.Height(), 1 );
        //--------------------------------------------------------------------//
        A21_MC_STAR = A21;
        A11_STAR_STAR = A11;
        advanced::internal::PanelLU
        ( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR, pivotOffset );

        advanced::internal::ComposePivots
        ( p1_STAR_STAR, image, preimage, pivotOffset );
        advanced::internal::ApplyRowPivots( AB, image, preimage, pivotOffset );

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR = A12;
        basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, (F)1, A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR = A12_STAR_VR;
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, A21_MC_STAR, A12_STAR_MR, (F)1, A22 );

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        A21 = A21_MC_STAR;
        p1 = p1_STAR_STAR;
        //--------------------------------------------------------------------//
        A12_STAR_VR.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A21_MC_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown
        ( pT,  p0,
               p1,
         /**/ /**/
          pB,  p2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::advanced::LU
( DistMatrix<float,MC,MR>& A, DistMatrix<int,VC,STAR>& p );

template void
elemental::advanced::LU
( DistMatrix<double,MC,MR>& A, DistMatrix<int,VC,STAR>& p );

#ifndef WITHOUT_COMPLEX
template void
elemental::advanced::LU
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<int,VC,STAR>& p );

template void
elemental::advanced::LU
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<int,VC,STAR>& p );
#endif
