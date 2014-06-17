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
#include "El-lite.hpp"

namespace El {

template<typename T>
RmaInterface<T>::RmaInterface( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::RmaInterface"))
}

template<typename T>
RmaInterface<T>::RmaInterface( const DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::RmaInterface"))
}

template<typename T>
RmaInterface<T>::~RmaInterface()
{
    {
        if( std::uncaught_exception() )
        {
           std::ostringstream os;
           os << "Uncaught exception detected during RmaInterface destructor "
                 "that required a call to Detach. Instead of allowing for the "
                 "possibility of Detach throwing another exception and "
                 "resulting in a 'terminate', we instead immediately dump the "
                 "call stack (if not in RELEASE mode) since the program will "
                 "likely hang:" << std::endl;
           std::cerr << os.str();
           DEBUG_ONLY(DumpCallStack())
        }
        else
        {
            Detach();
        }
    }
}

template<typename T>
void RmaInterface<T>::Attach( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Attach"))
}

template<typename T>
void RmaInterface<T>::Attach( const DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Attach"))
}

template<typename T>
void RmaInterface<T>::Put( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Put"))
}

template<typename T>
void RmaInterface<T>::Put( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Put"))
}

template<typename T>
void RmaInterface<T>::Get( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Get"))
}

template<typename T>
void RmaInterface<T>::Get( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Get"))
}

template<typename T>
void RmaInterface<T>::Acc( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Acc"))
}

template<typename T>
void RmaInterface<T>::Acc( const Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Acc"))
}

template<typename T>
void RmaInterface<T>::Detach()
{
    DEBUG_ONLY(CallStackEntry cse("RmaInterface::Detach"))
}

template class RmaInterface<Int>;
template class RmaInterface<float>;
template class RmaInterface<double>;
template class RmaInterface<Complex<float>>;
template class RmaInterface<Complex<double>>;

} // namespace El
