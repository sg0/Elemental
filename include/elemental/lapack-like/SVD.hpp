/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {

namespace svd {

template<typename R>
inline void
SimpleSVD
( DistMatrix<R,MC,  MR>& A,
  DistMatrix<R,VR,STAR>& s,
  DistMatrix<R,MC,  MR>& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVD");
#endif
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int subdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& grid = A.Grid();

    // Bidiagonalize A
    Bidiag( A );

    // Grab copies of the diagonal and superdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( grid ), 
                          e_MD_STAR( grid );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    DistMatrix<R,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                            e_STAR_STAR( e_MD_STAR );

    // Initialize U and VTrans to the appropriate identity matrices.
    DistMatrix<R,VC,STAR> U_VC_STAR( grid );
    DistMatrix<R,STAR,VR> VTrans_STAR_VR( grid );
    U_VC_STAR.AlignWith( A );
    VTrans_STAR_VR.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( k, n, VTrans_STAR_VR );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VTrans
    Matrix<R>& ULocal = U_VC_STAR.LocalMatrix();
    Matrix<R>& VTransLocal = VTrans_STAR_VR.LocalMatrix();
    lapack::BidiagQRAlg
    ( uplo, k, VTransLocal.Width(), ULocal.Height(),
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      VTransLocal.Buffer(), VTransLocal.LDim(), 
      ULocal.Buffer(), ULocal.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VTrans into a standard matrix dist.
    DistMatrix<R,MC,MR> B( A );
    if( m >= n )
    {
        DistMatrix<R,MC,MR> AT( grid ),
                            AB( grid );
        DistMatrix<R,VC,STAR> UT_VC_STAR( grid ), 
                              UB_VC_STAR( grid );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        Transpose( VTrans_STAR_VR, V );
    }
    else
    {
        DistMatrix<R,MC,MR> VT( grid ), 
                            VB( grid );
        DistMatrix<R,STAR,VR> VTransL_STAR_VR( grid ), VTransR_STAR_VR( grid );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionRight( VTrans_STAR_VR, VTransL_STAR_VR, VTransR_STAR_VR, m );
        Transpose( VTransL_STAR_VR, VT );
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, 0, B, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, 1, B, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, -1, B, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, 0, B, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SimpleSVD
( DistMatrix<Complex<R>,MC,  MR>& A,
  DistMatrix<R,         VR,STAR>& s,
  DistMatrix<Complex<R>,MC,  MR>& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVD");
#endif
    typedef Complex<R> C;
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int subdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& grid = A.Grid();

    // Bidiagonalize A
    DistMatrix<C,STAR,STAR> tP( grid ), tQ( grid );
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and superdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( grid ),
                          e_MD_STAR( grid );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, subdiagonal );

    // on each process
    DistMatrix<R,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                            e_STAR_STAR( e_MD_STAR );

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<C,VC,STAR> U_VC_STAR( grid );
    DistMatrix<C,STAR,VR> VAdj_STAR_VR( grid );
    U_VC_STAR.AlignWith( A );
    VAdj_STAR_VR.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( k, n, VAdj_STAR_VR );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VAdj
    Matrix<C>& ULocal = U_VC_STAR.LocalMatrix();
    Matrix<C>& VAdjLocal = VAdj_STAR_VR.LocalMatrix();
    lapack::BidiagQRAlg
    ( uplo, k, VAdjLocal.Width(), ULocal.Height(),
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      VAdjLocal.Buffer(), VAdjLocal.LDim(), 
      ULocal.Buffer(), ULocal.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VAdj into a standard matrix dist.
    DistMatrix<C,MC,MR> B( A );
    if( m >= n )
    {
        DistMatrix<C,MC,MR> AT( grid ),
                            AB( grid );
        DistMatrix<C,VC,STAR> UT_VC_STAR( grid ),
                              UB_VC_STAR( grid );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        Adjoint( VAdj_STAR_VR, V );
    }
    else
    {
        DistMatrix<C,MC,MR> VT( grid ), 
                            VB( grid );
        DistMatrix<C,STAR,VR> VAdjL_STAR_VR( grid ), VAdjR_STAR_VR( grid );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionRight( VAdj_STAR_VR, VAdjL_STAR_VR, VAdjR_STAR_VR, m );
        Adjoint( VAdjL_STAR_VR, VT );
        MakeZeros( VB );
    }

    // Backtransform U and VAdj
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 1, B, tP, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, -1, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, B, tP, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
CheckScale
( DistMatrix<F,MC,MR>& A, bool& needRescaling, typename Base<F>::type& scale )
{
    typedef typename Base<F>::type R;

    scale = 1;
    needRescaling = false;
    const R oneNormOfA = Norm( A, ONE_NORM );
    const R safeMin = lapack::MachineSafeMin<R>();
    const R precision = lapack::MachinePrecision<R>();
    const R smallNumber = safeMin/precision;
    const R bigNumber = 1/smallNumber;
    const R rhoMin = Sqrt(smallNumber);
    const R rhoMax = std::min( Sqrt(bigNumber), 1/Sqrt(Sqrt(safeMin)) );

    if( oneNormOfA > 0 && oneNormOfA < rhoMin )
    {
        needRescaling = true;
        scale = rhoMin/oneNormOfA;
    }
    else if( oneNormOfA > rhoMax )
    {
        needRescaling = true;
        scale = rhoMax/oneNormOfA;
    }
}

} // namespace svd

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H.              //
// On exit, A is overwritten with U.                                          //
//----------------------------------------------------------------------------//
template<typename F>
inline void
SVD
( DistMatrix<F,                     MC,  MR>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F,                     MC,  MR>& V )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    typedef typename Base<F>::type R;

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    svd::CheckScale( A, needRescaling, scale );
    if( needRescaling )
        Scal( (F)scale, A );

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::SimpleSVD( A, s, V );
    }
    else
    {
        // Lower bidiagonalization is not yet supported, so we instead play a 
        // trick to get the SVD of A.
        Adjoint( A, V );
        svd::SimpleSVD( V, s, A );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        Scal( 1/scale, s );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem