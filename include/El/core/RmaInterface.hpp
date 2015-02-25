/*
   Copyright (c) 2009-2014, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   Copyright (c) 2014, Jeff Hammond (Intel)
   All rights reserved.

   Authors:
   Jeff Hammond adapted the RMA interface from the AXPY one.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RMAINTERFACE_HPP
#define EL_RMAINTERFACE_HPP

namespace El {
#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY)
template<typename T>
class RmaInterface
{
public:
    RmaInterface();
    ~RmaInterface();

    RmaInterface(       DistMatrix<T,MC,MR>& Z );
    RmaInterface( const DistMatrix<T,MC,MR>& Z );

    void Attach(       DistMatrix<T,MC,MR>& Z );
    void Attach( const DistMatrix<T,MC,MR>& Z );

    // Local completion
    void Put(       Matrix<T>& Z, Int i, Int j );
    void Put( const Matrix<T>& Z, Int i, Int j );

    void Get(       Matrix<T>& Z, Int i, Int j );

    void Acc(       Matrix<T>& Z, Int i, Int j );
    void Acc( const Matrix<T>& Z, Int i, Int j );

    // No local completion
    void Iput(       Matrix<T>& Z, Int i, Int j );
    void Iput( const Matrix<T>& Z, Int i, Int j );

    void Iacc(       Matrix<T>& Z, Int i, Int j );
    void Iacc( const Matrix<T>& Z, Int i, Int j );

    // Request based RMA
    void Rput(       Matrix<T>& Z, Int i, Int j );
    void Rput( const Matrix<T>& Z, Int i, Int j );

    void Racc(       Matrix<T>& Z, Int i, Int j );
    void Racc( const Matrix<T>& Z, Int i, Int j );

    // Synchronization routines
    void Flush(            Matrix<T>& Z );
    void Flush(      const Matrix<T>& Z );
    void LocalFlush( const Matrix<T>& Z );
    void LocalFlush(       Matrix<T>& Z );
    void LocalFlush();

    void Detach();

private:

    mpi::Window window;

    // struct for passing data
    // for request based rma
    struct matrix_params_
    {
        const void *base_;
        std::vector<std::deque<std::vector<T>>>
        data_;
        std::vector<std::deque<mpi::Request>>
                                           requests_;
        std::vector<std::deque<bool>>
                                   statuses_;
    };

    std::vector<struct matrix_params_> matrices_;

    // buffers for rma
    std::vector<std::deque<std::vector<T>>>
    getVector_, putVector_;

    DistMatrix<T,MC,MR>* GlobalArrayPut_;
    const DistMatrix<T,MC,MR>* GlobalArrayGet_;

    bool toBeAttachedForPut_, toBeAttachedForGet_,
         attached_, detached_;

    // next index for data
    Int NextIndex (
        Int dataSize,
        std::deque <std::vector<T>>& dataVectors );

    Int NextIndex (
        Int target,
        Int dataSize,
        const void* base_address,
        Int* mindex );

    // only relevant for request-based
    // passive RMA
    bool anyPendingXfers (       Matrix<T>& Z );
    bool anyPendingXfers ( const Matrix<T>& Z );

    bool Testall();
    bool Test(          Matrix<T>& Z );
    bool Test(    const Matrix<T>& Z );
    bool TestAny(       Matrix<T>& Z );
    bool TestAny( const Matrix<T>& Z );

    void Waitall();
    void Wait(          Matrix<T>& Z );
    void Wait(    const Matrix<T>& Z );
    void WaitAny(       Matrix<T>& Z );
    void WaitAny( const Matrix<T>& Z );
};
#endif // EL_ENABLE_RMA_AXPY
} // namespace El
#endif // ifndef EL_RMAINTERFACE_HPP
