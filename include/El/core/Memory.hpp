/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MEMORY_HPP
#define EL_MEMORY_HPP

namespace El {

template<typename G>
class Memory
{
    size_t size_;
    G* buffer_;
#if MPI_VERSION>=3 && defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
    bool rma_;
#endif
public:
    Memory();
    Memory( size_t size );
    ~Memory();

    Memory( Memory<G>&& mem );
    Memory<G>& operator=( Memory<G>&& mem );

    // Exchange metadata with 'mem'
    void ShallowSwap( Memory<G>& mem );

    G* Buffer() const;
    size_t Size()   const;
#if MPI_VERSION>=3 && defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
    // rma get/set
    bool GetRMA() const;
    void SetRMA();
    // assign window baseptr to buffer_
    void Preallocated( size_t size, G * baseptr);
#endif
    G* Require( size_t size );
    void Release();
    void Empty();
};

} // namespace El

#endif // ifndef EL_MEMORY_HPP
