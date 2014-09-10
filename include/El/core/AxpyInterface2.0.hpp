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

    // nonblocking update routines
    // requires flush for completion
    void Iput( Matrix<T>& Z, Int i, Int j );
    void Iput( const Matrix<T>& Z, Int i, Int j );

    void Iget(       Matrix<T>& Z, Int i, Int j );

    void Iacc(       Matrix<T>& Z, Int i, Int j );
    void Iacc( const Matrix<T>& Z, Int i, Int j );

    void Flush(       Matrix<T>& Z );
    void Flush( const Matrix<T>& Z );
    
    // blocking update routines
    void Put( Matrix<T>& Z, Int i, Int j );
    void Put( const Matrix<T>& Z, Int i, Int j );

    void Get(       Matrix<T>& Z, Int i, Int j );

    void Acc(       Matrix<T>& Z, Int i, Int j );
    void Acc( const Matrix<T>& Z, Int i, Int j );

    void Detach();

private:
   
    static const Int 
        DATA_PUT_TAG      =1, 
        DATA_GET_TAG      =2,
        DATA_ACC_TAG   	  =3,
        REQUEST_GET_TAG   =4,
	COORD_IJ_TAG      =5;

    // struct for passing data
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

    // struct for passing coordinates
    struct coord_params_
    {
	T *base_;
	std::vector<std::deque<std::array<Int, 2>>>
	    coord_;
	std::vector<std::deque<mpi::Request>> 
	    requests_;
	std::vector<std::deque<bool>> 
	    statuses_;
    };
        	
    std::vector<struct coord_params_> coords_;

    // TODO need to add const here...
    DistMatrix<T,MC,MR>* GlobalArrayPut_;
    DistMatrix<T,MC,MR>* GlobalArrayGet_;
   
    bool toBeAttachedForPut_, toBeAttachedForGet_, 
	 attached_, detached_;
     
    Int NextIndexMatrix ( 
	Int target,
	Int dataSize, 
	T * base_address,
	Int *matrix_index);
 
    Int NextIndexMatrix ( 
	Int target,
	Int dataSize, 
	const T * base_address,
	Int *matrix_index);
 
    Int NextIndexCoord (
	Int i, Int j,
	Int target,
	T * base_address,
	Int *matrix_index);
 
    Int NextIndexCoord ( 
	Int i, Int j,
	Int target,
	const T * base_address,
	Int *matrix_index);
   
    /* Test */
    bool TestMatrix      ( Matrix<T>& Z );
    bool TestMatrix      ( const Matrix<T>& Z );
    
    bool TestCoord      ( Matrix<T>& Z );
    bool TestCoord      ( const Matrix<T>& Z );

    /* Wait */
    void WaitMatrix      ( Matrix<T>& Z );
    void WaitMatrix      ( const Matrix<T>& Z );

    // these are only used for nonblocking
    // update rountines
    void HandleGlobalToLocalData( Matrix<T>& Z );
    void HandleLocalToGlobalData( Matrix<T>& Z, Int count, Int source );
    void HandleLocalToGlobalAcc(  Matrix<T>& Z, Int count, Int source );

    void HandleGlobalToLocalData( const Matrix<T>& Z );
    void HandleLocalToGlobalData( const Matrix<T>& Z, Int count, Int source );
    void HandleLocalToGlobalAcc(  const Matrix<T>& Z, Int count, Int source );
};
} // namespace El
#endif // ifndef EL_AXPYINTERFACE2_HPP
