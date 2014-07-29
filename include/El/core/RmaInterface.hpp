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
#if MPI_VERSION>=3
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

    void Put( Matrix<T>& Z, Int i, Int j );
    void Put( const Matrix<T>& Z, Int i, Int j );

    void Get(       Matrix<T>& Z, Int i, Int j );

    void Acc(       Matrix<T>& Z, Int i, Int j );
    void Acc( const Matrix<T>& Z, Int i, Int j );
    void LocalAcc(  Matrix<T>& Z, Int i, Int j );

    void Flush(       Matrix<T>& Z, Int i, Int j );
    void Flush( const Matrix<T>& Z, Int i, Int j );
    void Flush(       Matrix<T>& Z );
    void Flush( const Matrix<T>& Z );

    void Detach();

private:
    mpi::Window window;

   std::vector<std::vector<T>>
        getVector_, putVector_;

    DistMatrix<T,MC,MR>* GlobalArrayPut_;
    const DistMatrix<T,MC,MR>* GlobalArrayGet_;
    
    bool toBeAttachedForPut_, toBeAttachedForGet_, 
	 attached_, detached_;
 };
#endif //MPI-3
} // namespace El
#endif // ifndef EL_RMAINTERFACE_HPP
