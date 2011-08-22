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
#include <ctime>
#include "elemental.hpp"
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

void Usage()
{
    cout << "Generates random matrix then solves for its QR factorization.\n\n"
         << "  QR <r> <c> <m> <n> <nb> <correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  m: height of matrix\n"
         << "  n: width of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename R> // represents a real number
void TestCorrectness
( bool printMatrices,
  const DistMatrix<R,MC,MR>& A,
        DistMatrix<R,MC,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int minDim = std::min(m,n);

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<R,MC,MR> Z(m,n,g);
    Z.SetToIdentity();
    advanced::ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, Z );
    advanced::ApplyPackedReflectors( LEFT, LOWER, VERTICAL, FORWARD, 0, A, Z );

    DistMatrix<R,MC,MR> ZUpper(g);
    ZUpper.View( Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<R,MC,MR> X(minDim,minDim,g);
    X.SetToIdentity();

    // Form X := I - Q^H Q
    basic::Axpy( (R)-1, ZUpper, X );

    R oneNormOfError = advanced::Norm( X, ONE_NORM );
    R infNormOfError = advanced::Norm( X, INFINITY_NORM );
    R frobNormOfError = advanced::Norm( X, FROBENIUS_NORM );
    if( g.VCRank() == 0 )
    {
        cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
             << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Testing if A = QR..." << endl;
    }

    // Form Q R
    DistMatrix<R,MC,MR> U( A );
    U.MakeTrapezoidal( LEFT, UPPER );
    advanced::ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, U );

    // Form Q R - A
    basic::Axpy( (R)-1, AOrig, U );
    
    R oneNormOfA = advanced::Norm( AOrig, ONE_NORM );
    R infNormOfA = advanced::Norm( AOrig, INFINITY_NORM );
    R frobNormOfA = advanced::Norm( AOrig, FROBENIUS_NORM );
    oneNormOfError = advanced::Norm( U, ONE_NORM );
    infNormOfError = advanced::Norm( U, INFINITY_NORM );
    frobNormOfError = advanced::Norm( U, FROBENIUS_NORM );
    if( g.VCRank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - QR||_1  = " << oneNormOfError << "\n"
             << "    ||A - QR||_oo = " << infNormOfError << "\n"
             << "    ||A - QR||_F  = " << frobNormOfError << endl;
    }
}

#ifndef WITHOUT_COMPLEX
template<typename R> // represents a real number
void TestCorrectness
( bool printMatrices,
  const DistMatrix<complex<R>,MC,  MR  >& A,
  const DistMatrix<complex<R>,STAR,STAR>& t,
        DistMatrix<complex<R>,MC,  MR  >& AOrig )
{
    typedef complex<R> C;

    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int minDim = std::min(m,n);

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<C,MC,MR> Z(m,n,g);
    Z.SetToIdentity();
    advanced::ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, CONJUGATED, 0, A, t, Z );
    advanced::ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, FORWARD, CONJUGATED, 0, A, t, Z );
    
    DistMatrix<C,MC,MR> ZUpper(g);
    ZUpper.View( Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<C,MC,MR> X(minDim,minDim,g);
    X.SetToIdentity();

    // Form X := I - Q^H Q
    basic::Axpy( (C)-1, ZUpper, X );

    R oneNormOfError = advanced::Norm( X, ONE_NORM );
    R infNormOfError = advanced::Norm( X, INFINITY_NORM );
    R frobNormOfError = advanced::Norm( X, FROBENIUS_NORM );
    if( g.VCRank() == 0 )
    {
        cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
             << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
    }

    if( g.VCRank() == 0 )
        cout << "  Testing if A = QR..." << endl;

    // Form Q R
    DistMatrix<C,MC,MR> U( A );
    U.MakeTrapezoidal( LEFT, UPPER );
    advanced::ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, CONJUGATED, 0, A, t, U );

    // Form Q R - A
    basic::Axpy( (C)-1, AOrig, U );
    
    R oneNormOfA = advanced::Norm( AOrig, ONE_NORM );
    R infNormOfA = advanced::Norm( AOrig, INFINITY_NORM );
    R frobNormOfA = advanced::Norm( AOrig, FROBENIUS_NORM );
    oneNormOfError = advanced::Norm( U, ONE_NORM );
    infNormOfError = advanced::Norm( U, INFINITY_NORM );
    frobNormOfError = advanced::Norm( U, FROBENIUS_NORM );
    if( g.VCRank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - QR||_1  = " << oneNormOfError << "\n"
             << "    ||A - QR||_oo = " << infNormOfError << "\n"
             << "    ||A - QR||_F  = " << frobNormOfError << endl;
    }
}
#endif // WITHOUT_COMPLEX

template<typename F> // represents a real or complex field
void TestQR
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g );

template<>
void TestQR<double>
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> A(g);
    DistMatrix<R,MC,MR> AOrig(g);

    A.ResizeTo( m, n );

    A.SetToRandom();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting QR factorization...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    advanced::QR( A );
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = advanced::internal::QRGFlops<R>( m, n, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, AOrig );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestQR< complex<double> >
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g )
{
    typedef complex<double> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,  MR  > A(g);
    DistMatrix<C,STAR,STAR> t(g);
    DistMatrix<C,MC,  MR  > AOrig(g);

    A.ResizeTo( m, n );

    A.SetToRandom();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting QR factorization...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    advanced::QR( A, t );
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = advanced::internal::QRGFlops<C>( m, n, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, t, AOrig );
}
#endif

int 
main( int argc, char* argv[] )
{
    Init( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 8 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try
    {
        int argNum = 0;
        const int r = atoi(argv[++argNum]);
        const int c = atoi(argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int n = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        const Grid g( comm, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
            cout << "Will test QR" << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestQR<double>
        ( testCorrectness, printMatrices, m, n, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestQR<dcomplex>
        ( testCorrectness, printMatrices, m, n, g );
#endif
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:\n"
             << e.what() << endl;
    }   
    Finalize();
    return 0;
}
