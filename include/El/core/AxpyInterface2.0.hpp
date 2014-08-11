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
    void LocalAcc(  Matrix<T>& Z, Int i, Int j );

    void Flush( Matrix<T>& Z, Int i, Int j );
    void Flush( const Matrix<T>& Z, Int i, Int j );
    void Flush(       Matrix<T>& Z );
    void Flush( const Matrix<T>& Z );
    
    void LocalFlush( Matrix<T>& Z, Int i, Int j );
    void LocalFlush( const Matrix<T>& Z, Int i, Int j );

    void Detach();

private:
   
    static const Int 
        DATA_PUT_TAG   =1, 
        DATA_GET_TAG   =2,
        DATA_ACC_TAG   =3,
        DATA_LCC_TAG   =4;
     
    /* Meta */
    std::vector<std::deque<mpi::Request>> 
	dataRequests_;
    std::vector<std::deque<bool>>  
	dataRequestStatuses_;
    std::vector<std::deque<T *>> 
	matrixBase_; 
    std::vector<std::deque< Int >> 
	opKind_; 
    /* Data */
    std::vector<std::deque<std::vector<T>>>
        getVectors_, putVectors_;
 
    DistMatrix<T,MC,MR>* GlobalArrayPut_;
    DistMatrix<T,MC,MR>* GlobalArrayGet_;
    
    bool toBeAttachedForPut_, toBeAttachedForGet_, 
	 attached_, detached_, sends_complete_;
     
   Int GetIndexForMatrix ( Matrix<T>& Z, const Int rank );
   void ProgressMatrix ( Matrix<T>& Z, const Int rank );
   Int GetMatrixType ( Matrix<T>& Z, const Int rank );

   Int NextIndex (Int dataSize, 
	   std::deque<std::vector<T>> &dataVectors,
	   std::deque<mpi::Request> &requests,
	   std::deque<bool> &requestStatus,
	   std::deque<Int> &opKind,
	   Int op,
	   std::deque<T *> &matrixBase,
	   T * base);
    
    void HandleLocalToGlobalData( Matrix<T>& Z, Int i, Int j );
    void HandleGlobalToLocalData( Matrix<T>& Z, Int i, Int j );
    void HandleLocalToGlobalAcc(  Matrix<T>& Z, Int i, Int j );
    void HandleGlobalToLocalAcc(  Matrix<T>& Z, Int i, Int j );
};
} // namespace El
#endif // ifndef EL_AXPYINTERFACE2_HPP
