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
	COORD_IJ_TAG      =5;

    struct matrix_params_
    {
	T *base_;
	std::vector<std::deque<std::vector<T>>>
	    data_;
	std::vector<std::deque<mpi::Request>> 
	    requests_;
	std::vector<std::deque<bool>> 
	    statuses_;
    };
        	
    std::vector<struct matrix_params_> matrices_;

    // need to add const here...
    DistMatrix<T,MC,MR>* GlobalArrayPut_;
    DistMatrix<T,MC,MR>* GlobalArrayGet_;
   
    bool toBeAttachedForPut_, toBeAttachedForGet_, 
	 attached_, detached_;
     
    Int NextIndex ( 
	Int target,
	Int dataSize, 
	T * base_address,
	Int *matrix_index);
 
    Int NextIndex ( 
	Int target,
	Int dataSize, 
	const T * base_address,
	Int *matrix_index);
    
    /* Test */
    bool TestRequests      ( Matrix<T>& Z );
    bool TestRequests      ( const Matrix<T>& Z );

    void HandleGlobalToLocalData( Matrix<T>& Z );
    void HandleLocalToGlobalData( Matrix<T>& Z, Int count, Int source );
    void HandleLocalToGlobalAcc(  Matrix<T>& Z, Int count, Int source );

    void HandleGlobalToLocalData( const Matrix<T>& Z );
    void HandleLocalToGlobalData( const Matrix<T>& Z, Int count, Int source );
    void HandleLocalToGlobalAcc(  const Matrix<T>& Z, Int count, Int source );
};
} // namespace El
#endif // ifndef EL_AXPYINTERFACE2_HPP
