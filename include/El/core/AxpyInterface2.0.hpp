/*
   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_AXPYINTERFACE2_HPP
#define EL_AXPYINTERFACE2_HPP

namespace El {
template<typename T>
class AxpyInterface2
{
public:
    AxpyInterface2();
    ~AxpyInterface2();

    AxpyInterface2(       DistMatrix<T,MC,MR>& Z );
    AxpyInterface2( const DistMatrix<T,MC,MR>& Z );

    void Attach(       DistMatrix<T,MC,MR>& Z );
    void Attach( const DistMatrix<T,MC,MR>& Z );

    void Put( Matrix<T>& Z, Int i, Int j );
    void Put( const Matrix<T>& Z, Int i, Int j );

    void Get(       Matrix<T>& Z, Int i, Int j );

    void Acc(       Matrix<T>& Z, Int i, Int j );
    void Acc( const Matrix<T>& Z, Int i, Int j );

    void Flush( Matrix<T>& Z, Int i, Int j );
    void Flush( const Matrix<T>& Z, Int i, Int j );
    void Flush(       Matrix<T>& Z );
    void Flush( const Matrix<T>& Z );
    
    void LocalFlush( Matrix<T>& Z, Int i, Int j );
    void LocalFlush( const Matrix<T>& Z, Int i, Int j );

    void Detach();

private:
   
    static const Int 
        DATA_PUT_TAG      =1, 
        DATA_GET_TAG      =2,
        DATA_ACC_TAG   	  =3,
        REQUEST_GET_TAG   =4,
	IJ_TAG		  =5;
     
    /* Request objects for send, recv and request op */
    std::vector<std::deque<mpi::Request>> 
	sendRequests_, requestRequests_, 
	recvRequests_, replyRequests_,
	sendIJRequests_;
    
    /* Request statuses for send, recv and request op */
    std::vector<std::deque<bool>>  
	sendRequestStatuses_, 
	recvRequestStatuses_, 
	requestRequestStatuses_,
	replyRequestStatuses_,
	sendIJRequestStatuses_;
    
    /* Stores matrix base addresses */
     std::vector<std::deque<T *>>  
	matrixBase_;

     /* Stores i, j coordinates */
    std::vector<std::deque<std::vector<Int>>>
	coordVectors_;
    
    /* Receive and Send vectors */
    std::vector<std::deque<std::vector<T>>>
        recvVectors_, sendVectors_, 
	replyVectors_, requestVectors_;

    // need to add const here...
    DistMatrix<T,MC,MR>* GlobalArrayPut_;
    DistMatrix<T,MC,MR>* GlobalArrayGet_;
   
    bool toBeAttachedForPut_, toBeAttachedForGet_, 
	 attached_, detached_;
     
    Int NextIndex ( Int dataSize, 
	    std::deque <std::vector<T>> &dataVectors,
	    std::deque <mpi::Request> &requests,
	    std::deque <bool> &requestStatus,
	    std::deque <T *> &matrixBase,
	    T * base_address );
    // note: this is just a placeholder
    // would be replaced soon
    Int NextIndex ( Int dataSize, 
	    std::deque <std::vector<T>> &dataVectors,
	    std::deque <mpi::Request> &requests,
	    std::deque <bool> &requestStatus,
	    std::deque <T *> &matrixBase,
	    const T * base_address );

    Int NextIndex
	( Int i, Int j, Int dataSize, 
	  std::deque <std::vector<Int>> &coordVectors,
	  std::deque <std::vector<T>> &dataVectors,
	  std::deque <mpi::Request> &requestData,
	  std::deque <bool> &requestDataStatus,
	  std::deque <mpi::Request> &requestCoord,
	  std::deque <bool> &requestCoordStatus,
	  std::deque <T *> &matrixBase,
	  T * base_address);

    Int NextIndex
	( Int i, Int j, Int dataSize, 
	  std::deque <std::vector<Int>> &coordVectors,
	  std::deque <std::vector<T>> &dataVectors,
	  std::deque <mpi::Request> &requestData,
	  std::deque <bool> &requestDataStatus,
	  std::deque <mpi::Request> &requestCoord,
	  std::deque <bool> &requestCoordStatus,
	  std::deque <T *> &matrixBase,
	  const T * base_address);

    /* Test */
    // probably we need const interfaces also?
    bool TestRequests      ( Matrix<T>& Z );
    bool TestReplies       ( Matrix<T>& Z );
    bool TestSends         ( Matrix<T>& Z );
    bool TestRecvs    	   ( Matrix<T>& Z );
    bool TestSendsCoord    ( Matrix<T>& Z );

    void HandleGlobalToLocalData( Matrix<T>& Z );
    void HandleLocalToGlobalData( Matrix<T>& Z, Int count, Int source );
    void HandleLocalToGlobalAcc(  Matrix<T>& Z, Int count, Int source );
};
} // namespace El
#endif // ifndef EL_AXPYINTERFACE2_HPP
