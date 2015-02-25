/*
   Copyright (c) 2009-2015, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   All rights reserved.

   Authors:
   This interface is mainly due to Martin Schatz, but it was put into its
   current form by Jack Poulson.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_AXPYINTERFACE_HPP
#define EL_AXPYINTERFACE_HPP

namespace El {

namespace AxpyTypeNS {
enum AxpyType { LOCAL_TO_GLOBAL, GLOBAL_TO_LOCAL };
}
using namespace AxpyTypeNS;

template<typename T>
class AxpyInterface
{
public:
    AxpyInterface();
    ~AxpyInterface();

    AxpyInterface( AxpyType type,       DistMatrix<T,MC,MR>& Z );
    AxpyInterface( AxpyType type, const DistMatrix<T,MC,MR>& Z );

    void Attach( AxpyType type,       DistMatrix<T,MC,MR>& Z );
    void Attach( AxpyType type, const DistMatrix<T,MC,MR>& Z );

    void Axpy( T alpha,       Matrix<T>& Z, Int i, Int j );
    void Axpy( T alpha, const Matrix<T>& Z, Int i, Int j );

    void Detach();

private:
#if MPI_VERSION>=3 && defined(EL_USE_IBARRIER_FOR_AXPY)
    static const Int
    DATA_TAG        =1,
    DATA_REQUEST_TAG=2,
    DATA_REPLY_TAG  =3;
#else
    static const Int
    DATA_TAG        =1,
    EOM_TAG         =2,
    DATA_REQUEST_TAG=3,
    DATA_REPLY_TAG  =4;
#endif

    //request object for polling on Issends
    bool attachedForLocalToGlobal_, attachedForGlobalToLocal_;

    DistMatrix<T,MC,MR>* localToGlobalMat_;
    const DistMatrix<T,MC,MR>* globalToLocalMat_;

#if MPI_VERSION>=3 && defined(EL_USE_IBARRIER_FOR_AXPY)
#else
    std::vector<bool> sentEomTo_, haveEomFrom_;
    std::vector<mpi::Request> eomSendRequests_;
#endif

    std::vector<std::deque<bool>>
                               sendingData_, sendingRequest_, sendingReply_;
    std::vector<std::deque<mpi::Request>>
                                       dataSendRequests_, requestSendRequests_, replySendRequests_;

    std::vector<byte> recvVector_;
    std::vector<std::deque<std::vector<byte>>>
    dataVectors_, requestVectors_, replyVectors_;

    byte sendDummy_, recvDummy_;

    // Progress functions
#if MPI_VERSION>=3 && defined(EL_USE_IBARRIER_FOR_AXPY)
    bool ReturnRequestStatuses();
#else
    // Check if we are done with this attachment's work
    bool Finished();
    void HandleEoms();
    void StartSendingEoms();
    void FinishSendingEoms();
    void UpdateRequestStatuses();
#endif

    Int ReadyForSend
    ( Int sendSize,
      std::deque<std::vector<byte>>& sendVectors,
      std::deque<mpi::Request>& requests,
      std::deque<bool>& requestStatuses );

    void HandleLocalToGlobalData();
    void HandleGlobalToLocalRequest();

    void AxpyLocalToGlobal( T alpha, const Matrix<T>& X, Int i, Int j );
    void AxpyGlobalToLocal( T alpha,       Matrix<T>& Y, Int i, Int j );
};

} // namespace El
#endif // ifndef EL_AXPYINTERFACE_HPP
