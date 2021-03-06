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
void Pei( Matrix<T>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    Ones( P, n, n );
    ShiftDiagonal( P, alpha );
}

template<typename T>
void Pei( AbstractDistMatrix<T>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    Ones( P, n, n );
    ShiftDiagonal( P, alpha );
}

template<typename T>
void Pei( AbstractBlockDistMatrix<T>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    Ones( P, n, n );
    ShiftDiagonal( P, alpha );
}

#define PROTO(T) \
  template void Pei( Matrix<T>& P, Int n, T alpha ); \
  template void Pei( AbstractDistMatrix<T>& P, Int n, T alpha ); \
  template void Pei( AbstractBlockDistMatrix<T>& P, Int n, T alpha );

#include "El/macros/Instantiate.h"

} // namespace El
