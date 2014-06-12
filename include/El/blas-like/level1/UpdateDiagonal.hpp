/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_UPDATEDIAGONAL_HPP
#define EL_UPDATEDIAGONAL_HPP

// This is essentially equivalent to SetDiagonal, but with s/Set/Update/g.

namespace El {

template<typename T,typename S>
inline void
UpdateDiagonal( Matrix<T>& A, S alpha, Int offset=0 )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
    {
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A.Update(i,j,alpha);
    }
}

template<typename T,typename S>
inline void
UpdateDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset=0 )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int i = j-offset;
        if( i >= 0 && i < height && A.IsLocalRow(i) )
        {
            const Int iLoc = A.LocalRow(i);
            A.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

} // namespace El

#endif // ifndef EL_UPDATEDIAGONAL_HPP