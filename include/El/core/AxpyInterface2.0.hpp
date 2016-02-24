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
#ifndef EL_AXPYINTERFACE2_HPP
#define EL_AXPYINTERFACE2_HPP

namespace El {

namespace TransferTypeNS {
enum TransferType { PUT, ACC, GET };
}
using namespace TransferTypeNS;

template<typename T>
class AxpyInterface2
{   
public:
    AxpyInterface2();
    ~AxpyInterface2();
   
    AxpyInterface2( DistMatrix<T,MC,MR>& Z );

    // collective
    void Attach( DistMatrix<T,MC,MR>& Z ); 

    // locally nonblocking 
    void Put( Matrix<T>& X, Int i, Int j );
    void Acc( Matrix<T>& X, Int i, Int j );
    void Get( Matrix<T>& Y, Int i, Int j );

    // remote completion
    void Flush( TransferType ttype );

    // collective
    void Detach();

private:
    /*
    static const Int 
        ACC_TAG        	=1, 
        EOM_ACC_TAG     =2, 
	GET_REQUEST_TAG =3, 
        GET_REPLY_TAG   =4,
        EOM_GET_TAG     =5;
    */
    static const Int 
        DATA_TAG        =1, 
        EOM_TAG         =2, 
        DATA_REQUEST_TAG=3, 
        DATA_REPLY_TAG  =4;
    
    bool attached_, detached_;

    // pointer to DistMatrix
    DistMatrix<T,MC,MR>* GlobalMat_;

    std::vector<bool> sentEomTo_, haveEomFrom_;
    std::vector<mpi::Request> eomSendRequests_;
    
    std::vector<std::deque<bool>> 
        sendingData_, sendingRequest_, sendingReply_;
    std::vector<std::deque<mpi::Request>> 
        dataSendRequests_, requestSendRequests_, replySendRequests_;
    
    std::vector<byte> recvVector_;
    std::vector<std::deque<std::vector<byte>>>
        dataVectors_, requestVectors_, replyVectors_;
    
    byte sendDummy_, recvDummy_;

    // Progress functions
    bool ReturnRequestStatuses();
    // Check if we are done with this attachment's work
    bool Finished();
    void HandleEoms();
    void StartSendingEoms();
    void FinishSendingEoms();
    void UpdateRequestStatuses();

    Int ReadyForSend
    ( Int sendSize,
      std::deque<std::vector<byte>>& sendVectors,
      std::deque<mpi::Request>& requests, 
      std::deque<bool>& requestStatuses );

    void HandleLocalToGlobalData( TransferType ttype ); // put/acc
    void HandleGlobalToLocalRequest(); // get

    void LocalToGlobal( Matrix<T>& X, Int i, Int j );
    void GlobalToLocal( Matrix<T>& Y, Int i, Int j );
};

} // namespace El
#endif // ifndef EL_AXPYINTERFACE2_HPP
