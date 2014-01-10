/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_TYPES_DECL_HPP
#define ELEM_CORE_TYPES_DECL_HPP

namespace elem {

typedef unsigned char byte;

// If these are changes, you must make sure that they have 
// existing MPI datatypes. This is only sometimes true for 'long long'
#ifdef USE_64BIT_INTS
typedef long long int Int;
typedef long long unsigned Unsigned;
#else
typedef int Int;
typedef unsigned Unsigned;
#endif
 
typedef Complex<float>  scomplex; 
typedef Complex<double> dcomplex;

struct BoolPair
{
    bool x, y;
    BoolPair( const BoolPair& bools ) : x(bools.x), y(bools.y) { }
    BoolPair( Int x_=0, Int y_=0 ) : x(x_), y(y_) { }
};
inline bool operator==( const BoolPair& a, const BoolPair& b )
{ return memcmp(&a,&b,sizeof(BoolPair))==0; }
inline bool operator!=( const BoolPair& a, const BoolPair& b )
{ return !operator==(a,b); }
inline BoolPair Transpose( const BoolPair& a )
{ return BoolPair(a.y,a.x); }

struct IndPair 
{ 
    Int i, j; 
    IndPair( const IndPair& inds ) : i(inds.i), j(inds.j) { }
    IndPair( Int i_=0, Int j_=0 ) : i(i_), j(j_) { }
};
inline bool operator==( const IndPair& a, const IndPair& b )
{ return memcmp(&a,&b,sizeof(IndPair))==0; }
inline bool operator!=( const IndPair& a, const IndPair& b )
{ return !operator==(a,b); }
inline IndPair Transpose( const IndPair& a )
{ return IndPair(a.j,a.i); }

struct DimPair 
{ 
    Int m, n; 

    DimPair( const DimPair& dims ) : m(dims.m), n(dims.n) { }
    DimPair( Int m_=0, Int n_=0 ) : m(m_), n(n_) { }
};
inline bool operator==( const DimPair& a, const DimPair& b )
{ return memcmp(&a,&b,sizeof(DimPair))==0; }
inline bool operator!=( const DimPair& a, const DimPair& b )
{ return !operator==(a,b); }
inline DimPair Transpose( const DimPair& a )
{ return DimPair(a.n,a.m); }

struct Layout
{
    DimPair dims;
    Int ldim;

    Layout( const Layout& layout ) 
    : dims(layout.dims), ldim(layout.ldim) { }
    Layout( const DimPair& dims_, Int ldim_ ) 
    : dims(dims_), ldim(ldim_) { }
    Layout( Int height, Int width, Int ldim_ ) 
    : dims(height,width), ldim(ldim_) { }

    Int Height() const { return dims.m; }
    Int Width() const { return dims.n; }
    Int LDim() const { return ldim; }
    DimPair Dimensions() const { return dims; }
};
inline bool operator==( const Layout& a, const Layout& b )
{ return memcmp(&a,&b,sizeof(Layout))==0; }
inline bool operator!=( const Layout& a, const Layout& b )
{ return !operator==(a,b); }

template<typename T>
struct FlatMatrix
{
    elem::Layout layout;
    T* buffer;

    FlatMatrix( const FlatMatrix<T>& fm ) 
    : layout(fm.layout), buffer(fm.buffer) { }
    FlatMatrix( const elem::Layout& layout_, T* buffer_ ) 
    : layout(layout_), buffer(buffer_) { }
    FlatMatrix( const DimPair& dims, Int ldim, T* buffer_ )
    : layout(dims,ldim), buffer(buffer_) { }
    FlatMatrix( Int height, Int width, Int ldim, T* buffer_ )
    : layout(height,width,ldim), buffer(buffer_) { }

    Int Height() const { return layout.Height(); }
    Int Width() const { return layout.Width(); }
    Int LDim() const { return layout.LDim(); }
    const elem::Layout& Layout() const { return layout; }
    DimPair Dimensions() const { return layout.Dimensions(); }
    T* Buffer() { return buffer; }
    const T* Buffer() const { return buffer; }
};
template<typename T>
inline bool operator==( const FlatMatrix<T>& a, const FlatMatrix<T>& b )
{ return memcmp(&a,&b,sizeof(FlatMatrix<T>))==0; }
template<typename T>
inline bool operator!=( const FlatMatrix<T>& a, const FlatMatrix<T>& b )
{ return !operator==(a,b); }

template<typename T>
struct LockedFlatMatrix
{
    elem::Layout layout;
    const T* buffer;

    LockedFlatMatrix( const FlatMatrix<T>& fm ) 
    : layout(fm.layout), buffer(fm.buffer) { }
    LockedFlatMatrix( const LockedFlatMatrix<T>& fm ) 
    : layout(fm.layout), buffer(fm.buffer) { }
    LockedFlatMatrix( const elem::Layout& layout_, const T* buffer_ ) 
    : layout(layout_), buffer(buffer_) { }
    LockedFlatMatrix( const DimPair& dims, Int ldim, const T* buffer_ )
    : layout(dims,ldim), buffer(buffer_) { }
    LockedFlatMatrix( Int height, Int width, Int ldim, const T* buffer_ )
    : layout(height,width,ldim), buffer(buffer_) { }

    Int Height() const { return layout.Height(); }
    Int Width() const { return layout.Width(); }
    Int LDim() const { return layout.LDim(); }
    const elem::Layout& Layout() const { return layout; }
    DimPair Dimensions() const { return layout.Dimensions(); }
    const T* Buffer() const { return buffer; }
};
template<typename T>
inline bool operator==
( const LockedFlatMatrix<T>& a, const LockedFlatMatrix<T>& b )
{ return memcmp(&a,&b,sizeof(LockedFlatMatrix<T>))==0; }
template<typename T>
inline bool operator!=
( const LockedFlatMatrix<T>& a, const LockedFlatMatrix<T>& b )
{ return !operator==(a,b); }

template<typename Real>
struct ValueInt
{
    Real value;
    Int index;

    static bool Lesser( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value < b.value; }
    static bool Greater( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value > b.value; }
};

template<typename Real>
struct ValueInt<Complex<Real> >
{
    Complex<Real> value;
    Int index;

    static bool Lesser( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return Abs(a.value) < Abs(b.value); }
    static bool Greater( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return Abs(a.value) > Abs(b.value); }
};

template<typename Real>
struct ValueIntPair
{
    Real value;
    Int indices[2];
    
    static bool Lesser( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value < b.value; }
    static bool Greater( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value > b.value; }
};

template<typename Real>
struct ValueIntPair<Complex<Real> >
{
    Complex<Real> value;
    Int indices[2];
    
    static bool Lesser
    ( const ValueIntPair<Real>& a, const ValueIntPair<Real>& b )
    { return Abs(a.value) < Abs(b.value); }
    static bool Greater
    ( const ValueIntPair<Real>& a, const ValueIntPair<Real>& b )
    { return Abs(a.value) > Abs(b.value); }
};

// For the safe computation of products. The result is given by 
//   product = rho * exp(kappa*n)
// where rho lies in (usually on) the unit circle and kappa is real-valued.
template<typename F>
struct SafeProduct
{
    F rho;
    BASE(F) kappa;
    Int n;

    SafeProduct( Int numEntries );
};

// The basic eigenvalue structure of a Hermitian matrix
struct Inertia
{
    Int numPositive, numNegative, numZero;
};

namespace conjugation_wrapper {
enum Conjugation
{
    UNCONJUGATED,
    CONJUGATED
};
}
using namespace conjugation_wrapper;

namespace distribution_wrapper {
enum Distribution
{
    MC,   // Col of a matrix distribution
    MD,   // Diagonal of a matrix distribution
    MR,   // Row of a matrix distribution
    VC,   // Col-major vector distribution
    VR,   // Row-major vector distribution
    STAR, // Give to every process
    CIRC  // Give to a single process
};
std::string DistToString( Distribution distribution );
Distribution StringToDist( std::string s );
}
using namespace distribution_wrapper;

namespace viewtype_wrapper {
enum ViewType
{
    OWNER = 0x0,
    VIEW = 0x1,
    OWNER_FIXED = 0x2,
    VIEW_FIXED = 0x3,
    LOCKED_OWNER = 0x4, // unused
    LOCKED_VIEW = 0x5,
    LOCKED_OWNER_FIXED = 0x6, // unused
    LOCKED_VIEW_FIXED = 0x7
};
static inline bool IsOwner( ViewType v ) 
{ return ( v & VIEW  ) == 0; }
static inline bool IsViewing( ViewType v )
{ return ( v & VIEW  ) != 0; }
static inline bool IsShrinkable( ViewType v )
{ return ( v & OWNER_FIXED ) == 0; }
static inline bool IsFixedSize( ViewType v )
{ return ( v & OWNER_FIXED ) != 0; }
static inline bool IsUnlocked( ViewType v )
{ return ( v & LOCKED_OWNER     ) == 0; }
static inline bool IsLocked( ViewType v )
{ return ( v & LOCKED_OWNER     ) != 0; }
}
using namespace viewtype_wrapper;

namespace forward_or_backward_wrapper {
enum ForwardOrBackward
{
    FORWARD,
    BACKWARD
};
}
using namespace forward_or_backward_wrapper;

namespace grid_order_wrapper {
enum GridOrder
{
    ROW_MAJOR,
    COLUMN_MAJOR
};
}
using namespace grid_order_wrapper;

namespace left_or_right_wrapper {
enum LeftOrRight
{
    LEFT,
    RIGHT
};
char LeftOrRightToChar( LeftOrRight side );
LeftOrRight CharToLeftOrRight( char c );
}
using namespace left_or_right_wrapper;

namespace sort_type_wrapper {
enum SortType
{
    UNSORTED,
    DESCENDING,
    ASCENDING
};
}
using namespace sort_type_wrapper;

namespace norm_type_wrapper {
enum NormType
{
    ONE_NORM,           // Operator one norm
    INFINITY_NORM,      // Operator infinity norm
    ENTRYWISE_ONE_NORM, // One-norm of vectorized matrix
    MAX_NORM,           // Maximum entry-wise magnitude
    NUCLEAR_NORM,       // One-norm of the singular values
    FROBENIUS_NORM,     // Two-norm of the singular values
    TWO_NORM            // Infinity-norm of the singular values
};
}
using namespace norm_type_wrapper;

namespace orientation_wrapper {
enum Orientation
{
    NORMAL,
    TRANSPOSE,
    ADJOINT
};
char OrientationToChar( Orientation orientation );
Orientation CharToOrientation( char c );
}
using namespace orientation_wrapper;

namespace unit_or_non_unit_wrapper {
enum UnitOrNonUnit
{
    NON_UNIT,
    UNIT
};
char UnitOrNonUnitToChar( UnitOrNonUnit diag );
UnitOrNonUnit CharToUnitOrNonUnit( char c );
}
using namespace unit_or_non_unit_wrapper;

namespace upper_or_lower_wrapper {
enum UpperOrLower
{
    LOWER,
    UPPER
};
char UpperOrLowerToChar( UpperOrLower uplo );
UpperOrLower CharToUpperOrLower( char c );
}
using namespace upper_or_lower_wrapper;

namespace vertical_or_horizontal_wrapper {
enum VerticalOrHorizontal
{
    VERTICAL,
    HORIZONTAL
};
}
using namespace vertical_or_horizontal_wrapper;

} // namespace elem

#endif // ifndef ELEM_CORE_TYPES_DECL_HPP
