/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T> 
void Wilkinson( Matrix<T>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("Wilkinson"))
    const Int n = 2*k+1;
    Zeros( A, n, n );
    FillDiagonal( A, T(1), -1 );
    FillDiagonal( A, T(1),  1 );
    
    for( Int j=0; j<=k; ++j )
        A.Set( j, j, T(k-j) );
    for( Int j=k+1; j<n; ++j )
        A.Set( j, j, T(j-k) );
}

template<typename T>
void Wilkinson( AbstractDistMatrix<T>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("Wilkinson"))
    const Int n = 2*k+1;
    Zeros( A, n, n );
    FillDiagonal( A, T(1), -1 );
    FillDiagonal( A, T(1),  1 );
    
    for( Int j=0; j<=k; ++j )
        A.Set( j, j, T(k-j) );
    for( Int j=k+1; j<n; ++j )
        A.Set( j, j, T(j-k) );
}

#define PROTO(T) \
  template void Wilkinson( Matrix<T>& A, Int k ); \
  template void Wilkinson( AbstractDistMatrix<T>& A, Int k );

#include "El/macros/Instantiate.h"

} // namespace El
