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
    void Put( Matrix<T>& Z, Int i, Int j );
    void Put( const Matrix<T>& Z, Int i, Int j );

    void Get(            Matrix<T>& Z, Int i, Int j );

    void Acc(       Matrix<T>& Z, Int i, Int j );
    void Acc( const Matrix<T>& Z, Int i, Int j );

    // No local completion 
    void Iput( Matrix<T>& Z, Int i, Int j );
    void Iput( const Matrix<T>& Z, Int i, Int j );

    void Iacc(       Matrix<T>& Z, Int i, Int j );
    void Iacc( const Matrix<T>& Z, Int i, Int j );

    void Iget( Matrix<T>& Z, Int i, Int j );

    // atomic routines

    // element-wise atomic increment, returns
    // previous value held in (i, j)
    T AtomicIncrement( Int i, Int j, T incr );

    // Synchronization routines
    void Flush(            Matrix<T>& Z );
    void Flush(      const Matrix<T>& Z );
    void LocalFlush( const Matrix<T>& Z );
    void LocalFlush(       Matrix<T>& Z );
    void Flush();
    void LocalFlush();

    void Detach();

private:
    // window for data
    mpi::Window window;

    // buffers for rma 
    std::vector<std::deque<std::vector<T>>>
        getVector_, putVector_;

    // UDD vector
    std::vector< mpi::Datatype > 
	uddtypes_;

    // metadata for pending gets
    struct pending_get_
    {
	bool is_active_;  	// is this tuple active?
	T * base_;              // base address of Matrix<T>
	Int index_;
	Int destination_;
	Int colShift_;
	Int rowShift_;
	Int remoteHeight_; 	// height of a DM chunk in remote process
	Int localHeight_;	// height/width of local chunk
	Int localWidth_;
	// initialize
	pending_get_() :
	    is_active_( true ),
	    base_( nullptr ),
	    index_( -1 ),
	    destination_( -1 ),
	    colShift_( -1 ),
	    rowShift_( -1 ),
	    remoteHeight_( -1 ),
	    localHeight_( -1 ),
	    localWidth_( -1 )
	{}
    };
    
    // vector of pending get parameters
    std::vector < pending_get_ > pending_gets_;

    DistMatrix<T,MC,MR>* GlobalArrayPut_;
    const DistMatrix<T,MC,MR>* GlobalArrayGet_;
    
    bool toBeAttachedForPut_, toBeAttachedForGet_, 
	 attached_, detached_, anyPendingGets_;

    // next index for data
    Int NextIndex ( 
	    Int dataSize, 
	    std::deque <std::vector<T>> &dataVectors );
};
#endif // EL_ENABLE_RMA_AXPY
} // namespace El
#endif // ifndef EL_RMAINTERFACE_HPP
