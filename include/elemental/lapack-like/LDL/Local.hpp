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

template<typename F>
inline void
LDLH( Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("LDLH");
#endif
    internal::LDLVar3( ADJOINT, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LDLT( Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("LDLT");
#endif
    internal::LDLVar3( TRANSPOSE, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

namespace internal {

// Unblocked serial LDL _without_ partial pivoting
template<typename F> 
inline void
LDLVar3( Orientation orientation, Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("LDLVar3");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( d.Viewing() && (d.Height() != A.Height() || d.Width() != 1) )
        throw std::logic_error
        ("d must be a column vector the same height as A");
    if( orientation == NORMAL )
        throw std::logic_error("Can only perform LDL^T or LDL^H");
#endif
    const int n = A.Height();
    if( !d.Viewing() )
        d.ResizeTo( n, 1 );
    if( n == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    F* ABuffer = A.Buffer();
    F* dBuffer = d.Buffer();
    const int ldim = A.LDim();
    for( int j=0; j<n; ++j )
    {
        const int a21Height = n - (j+1);

        // Extract and store the diagonal of D
        const F alpha11 = ABuffer[j+j*ldim];
        if( alpha11 == (F)0 )
            throw SingularMatrixException();
        dBuffer[j] = alpha11; 

        // A22 := A22 - a21 (a21/alpha11)^[T/H]
        const F* RESTRICT a21 = &ABuffer[(j+1)+j*ldim];
        if( orientation == ADJOINT )
        {
            for( int k=0; k<a21Height; ++k )
            {
                const F beta = Conj(a21[k]/alpha11);
                F* RESTRICT A22Col = &ABuffer[(j+1)+(j+1+k)*ldim];
                for( int i=k; i<a21Height; ++i )
                    A22Col[i] -= a21[i]*beta;
            }
        }
        else
        {
            for( int k=0; k<a21Height; ++k )
            {
                const F beta = a21[k]/alpha11;
                F* RESTRICT A22Col = &ABuffer[(j+1)+(j+1+k)*ldim];
                for( int i=k; i<a21Height; ++i )
                    A22Col[i] -= a21[i]*beta;
            }
        }
        
        // a21 := a21 / alpha11
        for( int i=0; i<a21Height; ++i )
            a21[i] /= alpha11;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

} // namespace elem
