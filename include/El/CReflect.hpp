/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CREFLECT_C_HPP
#define EL_CREFLECT_C_HPP

#define EL_CATCH \
  catch( std::bad_alloc& e ) \
  { El::ReportException(e); return EL_ALLOC_ERROR; } \
  catch( El::ArgException& e ) \
  { El::ReportException(e); return EL_ARG_ERROR; } \
  catch( std::logic_error& e ) \
  { El::ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { El::ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { El::ReportException(e); return EL_ERROR; }

#define EL_TRY(payload) \
  try { payload; } EL_CATCH \
  return EL_SUCCESS;

#define EL_RC(TYPE,INPUT) reinterpret_cast<TYPE>(INPUT)

namespace El {

template<typename T>
struct CReflectType { typedef T type; };

// ElInt and Int are typedef's
/*
template<> struct CReflectType<ElInt> { typedef Int type; };
template<> struct CReflectType<Int> { typedef ElInt type; };
*/

template<> struct CReflectType<complex_float> 
{ typedef Complex<float> type; };
template<> struct CReflectType<complex_double> 
{ typedef Complex<double> type; };

template<> struct CReflectType<Complex<float>> 
{ typedef complex_float type; };
template<> struct CReflectType<Complex<double>> 
{ typedef complex_double type; };

#define CREFLECT(T) typename CReflectType<T>::type

template<typename T>
inline void DynamicCastCheck( T* A )
{ if( A == nullptr ) RuntimeError("Dynamic cast failed"); }

inline string CReflect( const char* name ) { return string(name); }
// NOTE: This creates a deep copy and the pointer should be deleted later
inline char* CReflect( const string& name ) 
{
    const auto size = name.size();
    char* buffer = new char[size+1];
    memcpy( buffer, name.c_str(), size+1 );
    return buffer;
}

inline Range<Int> CReflect( ElRange_i rangeC )
{ return Range<Int>(rangeC.beg,rangeC.end); }
inline ElRange_i CReflect( Range<Int> range )
{ 
    ElRange_i rangeC; 
    rangeC.beg = range.beg; 
    rangeC.end = range.end; 
    return rangeC; 
}

inline Range<float> CReflect( ElRange_s rangeC )
{ return Range<float>(rangeC.beg,rangeC.end); }
inline ElRange_s CReflect( Range<float> range )
{ 
    ElRange_s rangeC; 
    rangeC.beg = range.beg; 
    rangeC.end = range.end; 
    return rangeC; 
}

inline Range<double> CReflect( ElRange_d rangeC )
{ return Range<double>(rangeC.beg,rangeC.end); }
inline ElRange_d CReflect( Range<double> range )
{ 
    ElRange_d rangeC; 
    rangeC.beg = range.beg; 
    rangeC.end = range.end; 
    return rangeC; 
}

inline Orientation CReflect( ElOrientation orient ) 
{ return static_cast<Orientation>(orient); }
inline ElOrientation CReflect( Orientation orient )
{ return static_cast<ElOrientation>(orient); }

inline LeftOrRight CReflect( ElLeftOrRight side )
{ return static_cast<LeftOrRight>(side); }
inline ElLeftOrRight CReflect( LeftOrRight side )
{ return static_cast<ElLeftOrRight>(side); }

inline UpperOrLower CReflect( ElUpperOrLower uplo )
{ return static_cast<UpperOrLower>(uplo); }
inline ElUpperOrLower CReflect( UpperOrLower uplo )
{ return static_cast<ElUpperOrLower>(uplo); }

inline UnitOrNonUnit CReflect( ElUnitOrNonUnit diag )
{ return static_cast<UnitOrNonUnit>(diag); }
inline ElUnitOrNonUnit CReflect( UnitOrNonUnit diag )
{ return static_cast<ElUnitOrNonUnit>(diag); }

inline VerticalOrHorizontal CReflect( ElVerticalOrHorizontal dir )
{ return static_cast<VerticalOrHorizontal>(dir); }
inline ElVerticalOrHorizontal CReflect( VerticalOrHorizontal dir )
{ return static_cast<ElVerticalOrHorizontal>(dir); }

inline ForwardOrBackward CReflect( ElForwardOrBackward order )
{ return static_cast<ForwardOrBackward>(order); }
inline ElForwardOrBackward CReflect( ForwardOrBackward order )
{ return static_cast<ElForwardOrBackward>(order); }

inline Conjugation CReflect( ElConjugation conjugation )
{ return static_cast<Conjugation>(conjugation); }
inline ElConjugation CReflect( Conjugation conjugation )
{ return static_cast<ElConjugation>(conjugation); }

// Dist
// ----
inline Dist   CReflect( ElDist dist ) { return static_cast<  Dist>(dist); }
inline ElDist CReflect(   Dist dist ) { return static_cast<ElDist>(dist); }

// Grid
// ----
inline   GridOrder     CReflect( ElGridOrderType order )
{ return static_cast<  GridOrder    >(order); }
inline ElGridOrderType CReflect(   GridOrder     order )
{ return static_cast<ElGridOrderType>(order); }

inline Grid* CReflect( ElGrid grid )
{ return EL_RC(Grid*,grid); }
inline ElGrid CReflect( Grid* grid )
{ return (ElGrid)EL_RC(struct ElGrid_sDummy*,grid); }

inline const Grid* CReflect( ElConstGrid grid )
{ return EL_RC(const Grid*,grid); }
inline ElConstGrid CReflect( const Grid* grid )
{ return (ElConstGrid)EL_RC(const struct ElGrid_sDummy*,grid); }

// Complex<T>
// ----------
inline complex_float* CReflect( Complex<float>* buffer )
{ return EL_RC(complex_float*,buffer); }

inline complex_double* CReflect( Complex<double>* buffer )
{ return EL_RC(complex_double*,buffer); }

inline const complex_float* CReflect( const Complex<float>* buffer )
{ return EL_RC(const complex_float*,buffer); }

inline const complex_double* CReflect( const Complex<double>* buffer )
{ return EL_RC(const complex_double*,buffer); }

inline Complex<float>* CReflect( complex_float* buffer )
{ return EL_RC(Complex<float>*,buffer); }

inline Complex<double>* CReflect( complex_double* buffer )
{ return EL_RC(Complex<double>*,buffer); }

inline const Complex<float>* CReflect( const complex_float* buffer )
{ return EL_RC(const Complex<float>*,buffer); }

inline const Complex<double>* CReflect( const complex_double* buffer )
{ return EL_RC(const Complex<double>*,buffer); }

inline Complex<float> CReflect( complex_float alpha )
{ return Complex<float>(alpha.real,alpha.imag); }

inline Complex<double> CReflect( complex_double alpha )
{ return Complex<double>(alpha.real,alpha.imag); }

inline complex_float CReflect( Complex<float> alpha )
{ complex_float beta; beta.real = alpha.real(); beta.imag = alpha.imag();
  return beta; }

inline complex_double CReflect( Complex<double> alpha )
{ complex_double beta; beta.real = alpha.real(); beta.imag = alpha.imag();
  return beta; }

// Analogues for real variables and integers
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
inline Int CReflect( Int alpha ) { return alpha; }
inline Int* CReflect( Int* buffer ) { return buffer; }
inline const Int* CReflect( const Int* buffer ) { return buffer; }
/*
inline ElInt CReflect( Int alpha ) { return alpha; }

inline ElInt* CReflect( Int*   buffer ) { return buffer; }
inline Int*   CReflect( ElInt* buffer ) { return buffer; }

inline const ElInt* CReflect( const Int*   buffer ) { return buffer; }
inline const Int*   CReflect( const ElInt* buffer ) { return buffer; }
*/

inline float CReflect( float alpha) { return alpha; }
inline double CReflect( double alpha ) { return alpha; }

inline float* CReflect( float* buffer ) { return buffer; }
inline double* CReflect( double* buffer ) { return buffer; }

inline const float* CReflect( const float* buffer ) { return buffer; }
inline const double* CReflect( const double* buffer ) { return buffer; }

inline ValueInt<Int> CReflect( ElValueInt_i entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_i CReflect( ValueInt<Int> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<Int>* CReflect( ElValueInt_i* entryC )
{ return EL_RC(ValueInt<Int>*,entryC); }
inline ElValueInt_i* CReflect( ValueInt<Int>* entryC )
{ return EL_RC(ElValueInt_i*,entryC); }

inline ValueInt<float> CReflect( ElValueInt_s entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_s CReflect( ValueInt<float> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<float>* CReflect( ElValueInt_s* entryC )
{ return EL_RC(ValueInt<float>*,entryC); }
inline ElValueInt_s* CReflect( ValueInt<float>* entryC )
{ return EL_RC(ElValueInt_s*,entryC); }

inline ValueInt<double> CReflect( ElValueInt_d entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_d CReflect( ValueInt<double> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<double>* CReflect( ElValueInt_d* entryC )
{ return EL_RC(ValueInt<double>*,entryC); }
inline ElValueInt_d* CReflect( ValueInt<double>* entryC )
{ return EL_RC(ElValueInt_d*,entryC); }

inline ValueInt<Complex<float>> CReflect( ElValueInt_c entryC )
{ return {CReflect(entryC.value),entryC.index}; }
inline ElValueInt_c CReflect( ValueInt<Complex<float>> entry )
{ return {CReflect(entry.value),entry.index}; }

inline ValueInt<Complex<float>>* CReflect( ElValueInt_c* entryC )
{ return EL_RC(ValueInt<Complex<float>>*,entryC); }
inline ElValueInt_c* CReflect( ValueInt<Complex<float>>* entryC )
{ return EL_RC(ElValueInt_c*,entryC); }

inline ValueInt<Complex<double>> CReflect( ElValueInt_z entryC )
{ return {CReflect(entryC.value),entryC.index}; }
inline ElValueInt_z CReflect( ValueInt<Complex<double>> entry )
{ return {CReflect(entry.value),entry.index}; }

inline ValueInt<Complex<double>>* CReflect( ElValueInt_z* entryC )
{ return EL_RC(ValueInt<Complex<double>>*,entryC); }
inline ElValueInt_z* CReflect( ValueInt<Complex<double>>* entryC )
{ return EL_RC(ElValueInt_z*,entryC); }

inline ValueIntPair<Int> CReflect( ElValueIntPair_i entryC )
{ return {entryC.value,{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_i CReflect( ValueIntPair<Int> entry )
{ return {entry.value,{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<float> CReflect( ElValueIntPair_s entryC )
{ return {entryC.value,{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_s CReflect( ValueIntPair<float> entry )
{ return {entry.value,{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<double> CReflect( ElValueIntPair_d entryC )
{ return {entryC.value,{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_d CReflect( ValueIntPair<double> entry )
{ return {entry.value,{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<Complex<float>> CReflect( ElValueIntPair_c entryC ) { return {CReflect(entryC.value),{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_c CReflect( ValueIntPair<Complex<float>> entry )
{ return {CReflect(entry.value),{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<Complex<double>> CReflect( ElValueIntPair_z entryC )
{ return {CReflect(entryC.value),{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_z CReflect( ValueIntPair<Complex<double>> entry )
{ return {CReflect(entry.value),{entry.indices[0],entry.indices[1]}}; }

// Matrix
// ------
inline Matrix<Int>* CReflect( ElMatrix_i A )
{ return EL_RC(Matrix<Int>*,A); }

inline Matrix<float>* CReflect( ElMatrix_s A )
{ return EL_RC(Matrix<float>*,A); }

inline Matrix<double>* CReflect( ElMatrix_d A )
{ return EL_RC(Matrix<double>*,A); }

inline Matrix<Complex<float>>* CReflect( ElMatrix_c A )
{ return EL_RC(Matrix<Complex<float>>*,A); }

inline Matrix<Complex<double>>* CReflect( ElMatrix_z A )
{ return EL_RC(Matrix<Complex<double>>*,A); }

inline const Matrix<Int>* CReflect( ElConstMatrix_i A )
{ return EL_RC(const Matrix<Int>*,A); }

inline const Matrix<float>* CReflect( ElConstMatrix_s A )
{ return EL_RC(const Matrix<float>*,A); }

inline const Matrix<double>* CReflect( ElConstMatrix_d A )
{ return EL_RC(const Matrix<double>*,A); }

inline const Matrix<Complex<float>>* CReflect( ElConstMatrix_c A )
{ return EL_RC(const Matrix<Complex<float>>*,A); }

inline const Matrix<Complex<double>>* CReflect( ElConstMatrix_z A )
{ return EL_RC(const Matrix<Complex<double>>*,A); }

inline ElMatrix_i CReflect( Matrix<Int>* A )
{ return (ElMatrix_i)EL_RC(struct ElMatrix_iDummy*,A); }

inline ElMatrix_s CReflect( Matrix<float>* A )
{ return (ElMatrix_s)EL_RC(struct ElMatrix_sDummy*,A); }

inline ElMatrix_d CReflect( Matrix<double>* A )
{ return (ElMatrix_d)EL_RC(struct ElMatrix_dDummy*,A); }

inline ElMatrix_c CReflect( Matrix<Complex<float>>* A )
{ return (ElMatrix_c)EL_RC(struct ElMatrix_cDummy*,A); }

inline ElMatrix_z CReflect( Matrix<Complex<double>>* A )
{ return (ElMatrix_z)EL_RC(struct ElMatrix_zDummy*,A); }

inline ElConstMatrix_i CReflect( const Matrix<Int>* A )
{ return (ElConstMatrix_i)EL_RC(const struct ElMatrix_iDummy*,A); }

inline ElConstMatrix_s CReflect( const Matrix<float>* A )
{ return (ElConstMatrix_s)EL_RC(const struct ElMatrix_sDummy*,A); }

inline ElConstMatrix_d CReflect( const Matrix<double>* A )
{ return (ElConstMatrix_d)EL_RC(const struct ElMatrix_dDummy*,A); }

inline ElConstMatrix_c CReflect( const Matrix<Complex<float>>* A )
{ return (ElConstMatrix_c)EL_RC(const struct ElMatrix_cDummy*,A); }

inline ElConstMatrix_z CReflect( const Matrix<Complex<double>>* A )
{ return (ElConstMatrix_z)EL_RC(const struct ElMatrix_zDummy*,A); }

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY) && defined(EL_ENABLE_RMA_GLOBAL_ARRAYS)
// GlobalArrays
// ------------
inline GlobalArrays<Int>* CReflect( ElGlobalArrays_i A )
{ return EL_RC(GlobalArrays<Int>*,A); }

inline GlobalArrays<float>* CReflect( ElGlobalArrays_s A )
{ return EL_RC(GlobalArrays<float>*,A); }

inline GlobalArrays<double>* CReflect( ElGlobalArrays_d A )
{ return EL_RC(GlobalArrays<double>*,A); }

inline ElGlobalArrays_i CReflect( GlobalArrays<Int>* A )
{ return (ElGlobalArrays_i)EL_RC(struct ElGlobalArrays_iDummy*,A); }

inline ElGlobalArrays_s CReflect( GlobalArrays<float>* A )
{ return (ElGlobalArrays_s)EL_RC(struct ElGlobalArrays_sDummy*,A); }

inline ElGlobalArrays_d CReflect( GlobalArrays<double>* A )
{ return (ElGlobalArrays_d)EL_RC(struct ElGlobalArrays_dDummy*,A); }
#endif

// AbstractDistMatrix
// ------------------
inline AbstractDistMatrix<Int>* 
CReflect( ElDistMatrix_i A )
{ return EL_RC(AbstractDistMatrix<Int>*,A); }

inline AbstractDistMatrix<float>* 
CReflect( ElDistMatrix_s A )
{ return EL_RC(AbstractDistMatrix<float>*,A); }

inline AbstractDistMatrix<double>* 
CReflect( ElDistMatrix_d A )
{ return EL_RC(AbstractDistMatrix<double>*,A); }

inline AbstractDistMatrix<Complex<float>>* 
CReflect( ElDistMatrix_c A )
{ return EL_RC(AbstractDistMatrix<Complex<float>>*,A); }

inline AbstractDistMatrix<Complex<double>>* 
CReflect( ElDistMatrix_z A )
{ return EL_RC(AbstractDistMatrix<Complex<double>>*,A); }

inline const AbstractDistMatrix<Int>* 
CReflect( ElConstDistMatrix_i A )
{ return EL_RC(const AbstractDistMatrix<Int>*,A); }

inline const AbstractDistMatrix<float>* 
CReflect( ElConstDistMatrix_s A )
{ return EL_RC(const AbstractDistMatrix<float>*,A); }

inline const AbstractDistMatrix<double>* 
CReflect( ElConstDistMatrix_d A )
{ return EL_RC(const AbstractDistMatrix<double>*,A); }

inline const AbstractDistMatrix<Complex<float>>* 
CReflect( ElConstDistMatrix_c A )
{ return EL_RC(const AbstractDistMatrix<Complex<float>>*,A); }

inline const AbstractDistMatrix<Complex<double>>* 
CReflect( ElConstDistMatrix_z A )
{ return EL_RC(const AbstractDistMatrix<Complex<double>>*,A); }

inline ElDistMatrix_i
CReflect( AbstractDistMatrix<Int>* A )
{ return (ElDistMatrix_i)EL_RC(struct ElDistMatrix_iDummy*,A); }

inline ElDistMatrix_s 
CReflect( AbstractDistMatrix<float>* A )
{ return (ElDistMatrix_s)EL_RC(struct ElDistMatrix_sDummy*,A); }

inline ElDistMatrix_d 
CReflect( AbstractDistMatrix<double>* A )
{ return (ElDistMatrix_d)EL_RC(struct ElDistMatrix_dDummy*,A); }

inline ElDistMatrix_c 
CReflect( AbstractDistMatrix<Complex<float>>* A )
{ return (ElDistMatrix_c)EL_RC(struct ElDistMatrix_cDummy*,A); }

inline ElDistMatrix_z 
CReflect( AbstractDistMatrix<Complex<double>>* A )
{ return (ElDistMatrix_z)EL_RC(struct ElDistMatrix_zDummy*,A); }

inline ElConstDistMatrix_i
CReflect( const AbstractDistMatrix<Int>* A )
{ return (ElConstDistMatrix_i)EL_RC(const struct ElDistMatrix_iDummy*,A); }

inline ElConstDistMatrix_s 
CReflect( const AbstractDistMatrix<float>* A )
{ return (ElConstDistMatrix_s)EL_RC(const struct ElDistMatrix_sDummy*,A); }

inline ElConstDistMatrix_d 
CReflect( const AbstractDistMatrix<double>* A )
{ return (ElConstDistMatrix_d)EL_RC(const struct ElDistMatrix_dDummy*,A); }

inline ElConstDistMatrix_c 
CReflect( const AbstractDistMatrix<Complex<float>>* A )
{ return (ElConstDistMatrix_c)EL_RC(const struct ElDistMatrix_cDummy*,A); }

inline ElConstDistMatrix_z 
CReflect( const AbstractDistMatrix<Complex<double>>* A )
{ return (ElConstDistMatrix_z)EL_RC(const struct ElDistMatrix_zDummy*,A); }

/* Graph
   ----- */
inline Graph* CReflect( ElGraph graph )
{ return EL_RC(Graph*,graph); }

inline const Graph* CReflect( ElConstGraph graph )
{ return EL_RC(const Graph*,graph); }

inline ElGraph CReflect( Graph* graph )
{ return EL_RC(ElGraph,graph); }

inline ElConstGraph CReflect( const Graph* graph )
{ return EL_RC(ElConstGraph,graph); }

/* DistGraph
   --------- */
inline DistGraph* CReflect( ElDistGraph graph )
{ return EL_RC(DistGraph*,graph); }

inline const DistGraph* CReflect( ElConstDistGraph graph )
{ return EL_RC(const DistGraph*,graph); }

inline ElDistGraph CReflect( DistGraph* graph )
{ return EL_RC(ElDistGraph,graph); }

inline ElConstDistGraph CReflect( const DistGraph* graph )
{ return EL_RC(ElConstDistGraph,graph); }

/* SparseMatrix
   ------------ */
inline SparseMatrix<Int>* CReflect( ElSparseMatrix_i A )
{ return EL_RC(SparseMatrix<Int>*,A); }

inline SparseMatrix<float>* CReflect( ElSparseMatrix_s A )
{ return EL_RC(SparseMatrix<float>*,A); }

inline SparseMatrix<double>* CReflect( ElSparseMatrix_d A )
{ return EL_RC(SparseMatrix<double>*,A); }

inline SparseMatrix<Complex<float>>* CReflect( ElSparseMatrix_c A )
{ return EL_RC(SparseMatrix<Complex<float>>*,A); }

inline SparseMatrix<Complex<double>>* CReflect( ElSparseMatrix_z A )
{ return EL_RC(SparseMatrix<Complex<double>>*,A); }

inline const SparseMatrix<Int>* CReflect( ElConstSparseMatrix_i A )
{ return EL_RC(const SparseMatrix<Int>*,A); }

inline const SparseMatrix<float>* CReflect( ElConstSparseMatrix_s A )
{ return EL_RC(const SparseMatrix<float>*,A); }

inline const SparseMatrix<double>* CReflect( ElConstSparseMatrix_d A )
{ return EL_RC(const SparseMatrix<double>*,A); }

inline const SparseMatrix<Complex<float>>* CReflect( ElConstSparseMatrix_c A )
{ return EL_RC(const SparseMatrix<Complex<float>>*,A); }

inline const SparseMatrix<Complex<double>>* CReflect( ElConstSparseMatrix_z A )
{ return EL_RC(const SparseMatrix<Complex<double>>*,A); }

inline ElSparseMatrix_i CReflect( SparseMatrix<Int>* A )
{ return (ElSparseMatrix_i)EL_RC(struct ElSparseMatrix_iDummy*,A); }

inline ElSparseMatrix_s CReflect( SparseMatrix<float>* A )
{ return (ElSparseMatrix_s)EL_RC(struct ElSparseMatrix_sDummy*,A); }

inline ElSparseMatrix_d CReflect( SparseMatrix<double>* A )
{ return (ElSparseMatrix_d)EL_RC(struct ElSparseMatrix_dDummy*,A); }

inline ElSparseMatrix_c CReflect( SparseMatrix<Complex<float>>* A )
{ return (ElSparseMatrix_c)EL_RC(struct ElSparseMatrix_cDummy*,A); }

inline ElSparseMatrix_z CReflect( SparseMatrix<Complex<double>>* A )
{ return (ElSparseMatrix_z)EL_RC(struct ElSparseMatrix_zDummy*,A); }

inline ElConstSparseMatrix_i CReflect( const SparseMatrix<Int>* A )
{ return (ElConstSparseMatrix_i)EL_RC(const struct ElSparseMatrix_iDummy*,A); }

inline ElConstSparseMatrix_s CReflect( const SparseMatrix<float>* A )
{ return (ElConstSparseMatrix_s)EL_RC(const struct ElSparseMatrix_sDummy*,A); }

inline ElConstSparseMatrix_d CReflect( const SparseMatrix<double>* A )
{ return (ElConstSparseMatrix_d)EL_RC(const struct ElSparseMatrix_dDummy*,A); }

inline ElConstSparseMatrix_c CReflect( const SparseMatrix<Complex<float>>* A )
{ return (ElConstSparseMatrix_c)EL_RC(const struct ElSparseMatrix_cDummy*,A); }

inline ElConstSparseMatrix_z CReflect( const SparseMatrix<Complex<double>>* A )
{ return (ElConstSparseMatrix_z)EL_RC(const struct ElSparseMatrix_zDummy*,A); }

/* DistSparseMatrix
   ---------------- */
inline DistSparseMatrix<Int>* CReflect( ElDistSparseMatrix_i A )
{ return EL_RC(DistSparseMatrix<Int>*,A); }

inline DistSparseMatrix<float>* CReflect( ElDistSparseMatrix_s A )
{ return EL_RC(DistSparseMatrix<float>*,A); }

inline DistSparseMatrix<double>* CReflect( ElDistSparseMatrix_d A )
{ return EL_RC(DistSparseMatrix<double>*,A); }

inline DistSparseMatrix<Complex<float>>* CReflect( ElDistSparseMatrix_c A )
{ return EL_RC(DistSparseMatrix<Complex<float>>*,A); }

inline DistSparseMatrix<Complex<double>>* CReflect( ElDistSparseMatrix_z A )
{ return EL_RC(DistSparseMatrix<Complex<double>>*,A); }

inline const DistSparseMatrix<Int>* CReflect( ElConstDistSparseMatrix_i A )
{ return EL_RC(const DistSparseMatrix<Int>*,A); }

inline const DistSparseMatrix<float>* CReflect( ElConstDistSparseMatrix_s A )
{ return EL_RC(const DistSparseMatrix<float>*,A); }

inline const DistSparseMatrix<double>* CReflect( ElConstDistSparseMatrix_d A )
{ return EL_RC(const DistSparseMatrix<double>*,A); }

inline const DistSparseMatrix<Complex<float>>* CReflect
( ElConstDistSparseMatrix_c A )
{ return EL_RC(const DistSparseMatrix<Complex<float>>*,A); }

inline const DistSparseMatrix<Complex<double>>* CReflect
( ElConstDistSparseMatrix_z A )
{ return EL_RC(const DistSparseMatrix<Complex<double>>*,A); }

inline ElDistSparseMatrix_i CReflect( DistSparseMatrix<Int>* A )
{ return (ElDistSparseMatrix_i)EL_RC(struct ElDistSparseMatrix_iDummy*,A); }

inline ElDistSparseMatrix_s CReflect( DistSparseMatrix<float>* A )
{ return (ElDistSparseMatrix_s)EL_RC(struct ElDistSparseMatrix_sDummy*,A); }

inline ElDistSparseMatrix_d CReflect( DistSparseMatrix<double>* A )
{ return (ElDistSparseMatrix_d)EL_RC(struct ElDistSparseMatrix_dDummy*,A); }

inline ElDistSparseMatrix_c CReflect( DistSparseMatrix<Complex<float>>* A )
{ return (ElDistSparseMatrix_c)EL_RC(struct ElDistSparseMatrix_cDummy*,A); }

inline ElDistSparseMatrix_z CReflect( DistSparseMatrix<Complex<double>>* A )
{ return (ElDistSparseMatrix_z)EL_RC(struct ElDistSparseMatrix_zDummy*,A); }

inline ElConstDistSparseMatrix_i CReflect( const DistSparseMatrix<Int>* A )
{ return (ElConstDistSparseMatrix_i)
  EL_RC(const struct ElDistSparseMatrix_iDummy*,A); }

inline ElConstDistSparseMatrix_s CReflect
( const DistSparseMatrix<float>* A )
{ return (ElConstDistSparseMatrix_s)
  EL_RC(const struct ElDistSparseMatrix_sDummy*,A); }

inline ElConstDistSparseMatrix_d CReflect
( const DistSparseMatrix<double>* A )
{ return (ElConstDistSparseMatrix_d)
  EL_RC(const struct ElDistSparseMatrix_dDummy*,A); }

inline ElConstDistSparseMatrix_c CReflect
( const DistSparseMatrix<Complex<float>>* A )
{ return (ElConstDistSparseMatrix_c)
  EL_RC(const struct ElDistSparseMatrix_cDummy*,A); }

inline ElConstDistSparseMatrix_z CReflect
( const DistSparseMatrix<Complex<double>>* A )
{ return (ElConstDistSparseMatrix_z)
  EL_RC(const struct ElDistSparseMatrix_zDummy*,A); }

/* DistMultiVec
   ------------ */
inline DistMultiVec<Int>* CReflect( ElDistMultiVec_i A )
{ return EL_RC(DistMultiVec<Int>*,A); }

inline DistMultiVec<float>* CReflect( ElDistMultiVec_s A )
{ return EL_RC(DistMultiVec<float>*,A); }

inline DistMultiVec<double>* CReflect( ElDistMultiVec_d A )
{ return EL_RC(DistMultiVec<double>*,A); }

inline DistMultiVec<Complex<float>>* CReflect( ElDistMultiVec_c A )
{ return EL_RC(DistMultiVec<Complex<float>>*,A); }

inline DistMultiVec<Complex<double>>* CReflect( ElDistMultiVec_z A )
{ return EL_RC(DistMultiVec<Complex<double>>*,A); }

inline const DistMultiVec<Int>* CReflect( ElConstDistMultiVec_i A )
{ return EL_RC(const DistMultiVec<Int>*,A); }

inline const DistMultiVec<float>* CReflect( ElConstDistMultiVec_s A )
{ return EL_RC(const DistMultiVec<float>*,A); }

inline const DistMultiVec<double>* CReflect( ElConstDistMultiVec_d A )
{ return EL_RC(const DistMultiVec<double>*,A); }

inline const DistMultiVec<Complex<float>>* CReflect( ElConstDistMultiVec_c A )
{ return EL_RC(const DistMultiVec<Complex<float>>*,A); }

inline const DistMultiVec<Complex<double>>* CReflect( ElConstDistMultiVec_z A )
{ return EL_RC(const DistMultiVec<Complex<double>>*,A); }

inline ElDistMultiVec_i CReflect( DistMultiVec<Int>* A )
{ return (ElDistMultiVec_i)EL_RC(struct ElDistMultiVec_iDummy*,A); }

inline ElDistMultiVec_s CReflect( DistMultiVec<float>* A )
{ return (ElDistMultiVec_s)EL_RC(struct ElDistMultiVec_sDummy*,A); }

inline ElDistMultiVec_d CReflect( DistMultiVec<double>* A )
{ return (ElDistMultiVec_d)EL_RC(struct ElDistMultiVec_dDummy*,A); }

inline ElDistMultiVec_c CReflect( DistMultiVec<Complex<float>>* A )
{ return (ElDistMultiVec_c)EL_RC(struct ElDistMultiVec_cDummy*,A); }

inline ElDistMultiVec_z CReflect( DistMultiVec<Complex<double>>* A )
{ return (ElDistMultiVec_z)EL_RC(struct ElDistMultiVec_zDummy*,A); }

inline ElConstDistMultiVec_i CReflect( const DistMultiVec<Int>* A )
{ return (ElConstDistMultiVec_i)EL_RC(const struct ElDistMultiVec_iDummy*,A); }

inline ElConstDistMultiVec_s CReflect( const DistMultiVec<float>* A )
{ return (ElConstDistMultiVec_s)EL_RC(const struct ElDistMultiVec_sDummy*,A); }

inline ElConstDistMultiVec_d CReflect( const DistMultiVec<double>* A )
{ return (ElConstDistMultiVec_d)EL_RC(const struct ElDistMultiVec_dDummy*,A); }

inline ElConstDistMultiVec_c CReflect( const DistMultiVec<Complex<float>>* A )
{ return (ElConstDistMultiVec_c)EL_RC(const struct ElDistMultiVec_cDummy*,A); }

inline ElConstDistMultiVec_z CReflect( const DistMultiVec<Complex<double>>* A )
{ return (ElConstDistMultiVec_z)EL_RC(const struct ElDistMultiVec_zDummy*,A); }

inline ElDistData CReflect( const DistData& data )
{
    ElDistData distData;
    distData.colDist = CReflect(data.colDist);
    distData.rowDist = CReflect(data.rowDist);
    distData.colAlign = data.colAlign;
    distData.rowAlign = data.rowAlign;
    distData.root = data.root;
    distData.grid = CReflect(data.grid);
    return distData;
}

inline DistData CReflect( const ElDistData& distData )
{
    DistData data;
    data.colDist = CReflect(distData.colDist);
    data.rowDist = CReflect(distData.rowDist);
    data.colAlign = distData.colAlign;
    data.rowAlign = distData.rowAlign;
    data.root = distData.root;
    data.grid = CReflect(distData.grid);
    return data;
}

inline ElSafeProduct_s CReflect( const SafeProduct<float>& prod )
{ 
    ElSafeProduct_s prodC;    
    prodC.rho = prod.rho;
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}
inline ElSafeProduct_d CReflect( const SafeProduct<double>& prod )
{ 
    ElSafeProduct_d prodC;    
    prodC.rho = prod.rho;
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}
inline ElSafeProduct_c CReflect( const SafeProduct<Complex<float>>& prod )
{ 
    ElSafeProduct_c prodC;    
    prodC.rho = CReflect(prod.rho);
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}
inline ElSafeProduct_z CReflect( const SafeProduct<Complex<double>>& prod )
{ 
    ElSafeProduct_z prodC;    
    prodC.rho = CReflect(prod.rho);
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}

inline SafeProduct<float> CReflect( const ElSafeProduct_s& prodC )
{ 
    SafeProduct<float> prod( prodC.n );
    prod.rho = prodC.rho;
    prod.kappa = prodC.kappa;
    return prod;
}
inline SafeProduct<double> CReflect( const ElSafeProduct_d& prodC )
{ 
    SafeProduct<double> prod( prodC.n );
    prod.rho = prodC.rho;
    prod.kappa = prodC.kappa;
    return prod;
}
inline SafeProduct<Complex<float>> CReflect( const ElSafeProduct_c& prodC )
{ 
    SafeProduct<Complex<float>> prod( prodC.n );
    prod.rho = CReflect(prodC.rho);
    prod.kappa = prodC.kappa;
    return prod;
}
inline SafeProduct<Complex<double>> CReflect( const ElSafeProduct_z& prodC )
{ 
    SafeProduct<Complex<double>> prod( prodC.n );
    prod.rho = CReflect(prodC.rho);
    prod.kappa = prodC.kappa;
    return prod;
}

// Input/Output
// ------------
inline ElFileFormat CReflect( FileFormat format )
{ return static_cast<ElFileFormat>(format); }
inline FileFormat CReflect( ElFileFormat format )
{ return static_cast<FileFormat>(format); }

inline ElColorMap CReflect( ColorMap map )
{ return static_cast<ElColorMap>(map); }
inline ColorMap CReflect( ElColorMap map )
{ return static_cast<ColorMap>(map); }

// BLAS-like
// ---------
inline ElGemmAlgorithm CReflect( GemmAlgorithm alg )
{ return static_cast<ElGemmAlgorithm>(alg); }
inline GemmAlgorithm CReflect( ElGemmAlgorithm alg )
{ return static_cast<GemmAlgorithm>(alg); }

template<typename T>
inline ElSymvCtrl
CReflect( const SymvCtrl<T>& ctrl )
{ 
    ElSymvCtrl ctrlC;
    ctrlC.bsize = ctrl.bsize;
    ctrlC.avoidTrmvBasedLocalSymv = ctrl.avoidTrmvBasedLocalSymv;
    return ctrlC;
}

template<typename T>
inline SymvCtrl<T>
CReflect( const ElSymvCtrl& ctrlC )
{ 
    SymvCtrl<T> ctrl;
    ctrl.bsize = ctrlC.bsize;
    ctrl.avoidTrmvBasedLocalSymv = ctrlC.avoidTrmvBasedLocalSymv;
    return ctrl;
}

// LAPACK-like
// -----------

inline ElSortType CReflect( SortType type )
{ return static_cast<ElSortType>(type); }

inline SortType CReflect( ElSortType type )
{ return static_cast<SortType>(type); }

// Permutations
// ^^^^^^^^^^^^

inline ElPermutationMeta CReflect( const PermutationMeta& meta )
{
    ElPermutationMeta metaC;    

    metaC.align = meta.align;
    metaC.comm = meta.comm.comm;

    const Int commSize = mpi::Size( meta.comm );
    metaC.sendCounts = new int[commSize];
    metaC.sendDispls = new int[commSize];
    metaC.recvCounts = new int[commSize];
    metaC.recvDispls = new int[commSize];
    MemCopy( metaC.sendCounts, meta.sendCounts.data(), commSize );
    MemCopy( metaC.sendDispls, meta.sendDispls.data(), commSize ); 
    MemCopy( metaC.recvCounts, meta.recvCounts.data(), commSize );
    MemCopy( metaC.recvDispls, meta.recvDispls.data(), commSize );

    metaC.numSendIdx = meta.sendIdx.size();
    metaC.numRecvIdx = meta.recvIdx.size();
    metaC.sendIdx   = new int[metaC.numSendIdx];
    metaC.sendRanks = new int[metaC.numSendIdx];
    metaC.recvIdx   = new int[metaC.numRecvIdx];
    metaC.recvRanks = new int[metaC.numRecvIdx];
    MemCopy( metaC.sendIdx,   meta.sendIdx.data(),   metaC.numSendIdx );
    MemCopy( metaC.sendRanks, meta.sendRanks.data(), metaC.numSendIdx );
    MemCopy( metaC.recvIdx,   meta.recvIdx.data(),   metaC.numRecvIdx );
    MemCopy( metaC.recvRanks, meta.recvRanks.data(), metaC.numRecvIdx );

    return metaC;
}

inline PermutationMeta CReflect( const ElPermutationMeta& metaC )
{
    PermutationMeta meta;

    meta.align = metaC.align;
    meta.comm = metaC.comm;

    int commSize;
    MPI_Comm_size( metaC.comm, &commSize );
    meta.sendCounts = 
        vector<int>( metaC.sendCounts, metaC.sendCounts+commSize );
    meta.sendDispls = 
        vector<int>( metaC.sendDispls, metaC.sendDispls+commSize );
    meta.recvCounts =
        vector<int>( metaC.recvCounts, metaC.recvCounts+commSize );
    meta.recvDispls =
        vector<int>( metaC.recvDispls, metaC.recvDispls+commSize );

    meta.sendIdx = 
        vector<int>( metaC.sendIdx, metaC.sendIdx+metaC.numSendIdx );
    meta.sendRanks =
        vector<int>( metaC.sendRanks, metaC.sendRanks+metaC.numSendIdx );
    meta.recvIdx =
        vector<int>( metaC.recvIdx, metaC.recvIdx+metaC.numRecvIdx );
    meta.recvRanks =
        vector<int>( metaC.recvRanks, metaC.recvRanks+metaC.numRecvIdx );

    return meta;
}

// Condensed form
// ^^^^^^^^^^^^^^
inline ElHermitianTridiagApproach 
CReflect( HermitianTridiagApproach approach )
{ return static_cast<ElHermitianTridiagApproach>( approach ); }

inline HermitianTridiagApproach 
CReflect( ElHermitianTridiagApproach approach )
{ return static_cast<HermitianTridiagApproach>( approach ); }

template<typename F>
inline ElHermitianTridiagCtrl
CReflect( const HermitianTridiagCtrl<F>& ctrl )
{ 
    ElHermitianTridiagCtrl ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.order = CReflect(ctrl.order);
    ctrlC.symvCtrl = CReflect(ctrl.symvCtrl);
    return ctrlC;
}

template<typename F>
inline HermitianTridiagCtrl<F>
CReflect( const ElHermitianTridiagCtrl& ctrlC )
{ 
    HermitianTridiagCtrl<F> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.order = CReflect(ctrlC.order);
    ctrl.symvCtrl = CReflect<F>(ctrlC.symvCtrl);
    return ctrl;
}

// Decompositions
// ^^^^^^^^^^^^^^

/* Pencil */
inline ElPencil CReflect( Pencil pencil )
{ return static_cast<ElPencil>(pencil); }

inline Pencil CReflect( ElPencil pencil )
{ return static_cast<Pencil>(pencil); }

/* HermitianSDCCtrl */
inline ElHermitianSDCCtrl_s CReflect( const HermitianSDCCtrl<float>& ctrl )
{
    ElHermitianSDCCtrl_s ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}
inline ElHermitianSDCCtrl_d CReflect( const HermitianSDCCtrl<double>& ctrl )
{
    ElHermitianSDCCtrl_d ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline HermitianSDCCtrl<float> CReflect( const ElHermitianSDCCtrl_s& ctrlC )
{
    HermitianSDCCtrl<float> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}
inline HermitianSDCCtrl<double> CReflect( const ElHermitianSDCCtrl_d& ctrlC )
{
    HermitianSDCCtrl<double> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* HermitianEigSubset */
inline ElHermitianEigSubset_s CReflect
( const HermitianEigSubset<float>& subset )
{
    ElHermitianEigSubset_s subsetC;
    subsetC.indexSubset = subset.indexSubset;
    subsetC.lowerIndex = subset.lowerIndex;
    subsetC.upperIndex = subset.upperIndex;
    subsetC.rangeSubset = subset.rangeSubset;
    subsetC.lowerBound = subset.lowerBound;
    subsetC.upperBound = subset.upperBound;
    return subsetC;
}
inline ElHermitianEigSubset_d CReflect
( const HermitianEigSubset<double>& subset )
{
    ElHermitianEigSubset_d subsetC;
    subsetC.indexSubset = subset.indexSubset;
    subsetC.lowerIndex = subset.lowerIndex;
    subsetC.upperIndex = subset.upperIndex;
    subsetC.rangeSubset = subset.rangeSubset;
    subsetC.lowerBound = subset.lowerBound;
    subsetC.upperBound = subset.upperBound;
    return subsetC;
}

inline HermitianEigSubset<float> CReflect
( const ElHermitianEigSubset_s& subsetC )
{
    HermitianEigSubset<float> subset;
    subset.indexSubset = subsetC.indexSubset;
    subset.lowerIndex = subsetC.lowerIndex;
    subset.upperIndex = subsetC.upperIndex;
    subset.rangeSubset = subsetC.rangeSubset;
    subset.lowerBound = subsetC.lowerBound;
    subset.upperBound = subsetC.upperBound;
    return subset;
}
inline HermitianEigSubset<double> CReflect
( const ElHermitianEigSubset_d& subsetC )
{
    HermitianEigSubset<double> subset;
    subset.indexSubset = subsetC.indexSubset;
    subset.lowerIndex = subsetC.lowerIndex;
    subset.upperIndex = subsetC.upperIndex;
    subset.rangeSubset = subsetC.rangeSubset;
    subset.lowerBound = subsetC.lowerBound;
    subset.upperBound = subsetC.upperBound;
    return subset;
}

/* HermitianEigCtrl */
inline ElHermitianEigCtrl_s CReflect( const HermitianEigCtrl<float>& ctrl )
{
    ElHermitianEigCtrl_s ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline ElHermitianEigCtrl_d CReflect( const HermitianEigCtrl<double>& ctrl )
{
    ElHermitianEigCtrl_d ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_c 
CReflect( const HermitianEigCtrl<Complex<float>>& ctrl )
{
    ElHermitianEigCtrl_c ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_z
CReflect( const HermitianEigCtrl<Complex<double>>& ctrl )
{
    ElHermitianEigCtrl_z ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline HermitianEigCtrl<float> CReflect( const ElHermitianEigCtrl_s& ctrlC )
{
    HermitianEigCtrl<float> ctrl;
    ctrl.tridiagCtrl = CReflect<float>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<double> CReflect( const ElHermitianEigCtrl_d& ctrlC )
{
    HermitianEigCtrl<double> ctrl;
    ctrl.tridiagCtrl = CReflect<double>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<float>> 
CReflect( const ElHermitianEigCtrl_c& ctrlC )
{
    HermitianEigCtrl<Complex<float>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<float>>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<double>> 
CReflect( const ElHermitianEigCtrl_z& ctrlC )
{
    HermitianEigCtrl<Complex<double>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<double>>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}

/* PolarCtrl */
inline ElPolarCtrl CReflect( const PolarCtrl& ctrl )
{
    ElPolarCtrl ctrlC;
    ctrlC.qdwh = ctrl.qdwh;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.maxIts = ctrl.maxIts;
    ctrlC.numIts = ctrl.numIts;
    return ctrlC;
}

inline PolarCtrl CReflect( const ElPolarCtrl& ctrlC )
{
    PolarCtrl ctrl;
    ctrl.qdwh = ctrlC.qdwh;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.maxIts = ctrlC.maxIts;
    ctrl.numIts = ctrlC.numIts;
    return ctrl;
}

/* SVDCtrl */
inline SVDCtrl<float> CReflect( const ElSVDCtrl_s& ctrlC )
{
    SVDCtrl<float> ctrl;
    ctrl.seqQR = ctrlC.seqQR;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.thresholded = ctrlC.thresholded;
    ctrl.relative = ctrlC.relative;
    ctrl.tol = ctrlC.tol;
    return ctrl;
}

inline SVDCtrl<double> CReflect( const ElSVDCtrl_d& ctrlC )
{
    SVDCtrl<double> ctrl;
    ctrl.seqQR = ctrlC.seqQR;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.thresholded = ctrlC.thresholded;
    ctrl.relative = ctrlC.relative;
    ctrl.tol = ctrlC.tol;
    return ctrl;
}

inline ElSVDCtrl_s CReflect( const SVDCtrl<float>& ctrl )
{
    ElSVDCtrl_s ctrlC;
    ctrlC.seqQR = ctrl.seqQR;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.thresholded = ctrl.thresholded;
    ctrlC.relative = ctrl.relative;
    ctrlC.tol = ctrl.tol;
    return ctrlC;
}

inline ElSVDCtrl_d CReflect( const SVDCtrl<double>& ctrl )
{
    ElSVDCtrl_d ctrlC;
    ctrlC.seqQR = ctrl.seqQR;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.thresholded = ctrl.thresholded;
    ctrlC.relative = ctrl.relative;
    ctrlC.tol = ctrl.tol;
    return ctrlC;
}

/* HessQRCtrl */
inline ElHessQRCtrl CReflect( const HessQRCtrl& ctrl )
{
    ElHessQRCtrl ctrlC;
    ctrlC.distAED = ctrl.distAED;
    ctrlC.blockHeight = ctrl.blockHeight;
    ctrlC.blockWidth = ctrl.blockWidth;
    return ctrlC;
}

inline HessQRCtrl CReflect( const ElHessQRCtrl& ctrlC )
{
    HessQRCtrl ctrl;
    ctrl.distAED = ctrlC.distAED;
    ctrl.blockHeight = ctrlC.blockHeight;
    ctrl.blockWidth = ctrlC.blockWidth;
    return ctrl;
}

/* SDCCtrl */
inline ElSDCCtrl_s CReflect( const SDCCtrl<float>& ctrl )
{
    ElSDCCtrl_s ctrlC;    
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.random = ctrl.random;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}
inline ElSDCCtrl_d CReflect( const SDCCtrl<double>& ctrl )
{
    ElSDCCtrl_d ctrlC;    
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.random = ctrl.random;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline SDCCtrl<float> CReflect( const ElSDCCtrl_s& ctrlC )
{
    SDCCtrl<float> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.random = ctrlC.random;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}
inline SDCCtrl<double> CReflect( const ElSDCCtrl_d& ctrlC )
{
    SDCCtrl<double> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.random = ctrlC.random;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* SchurCtrl */
inline ElSchurCtrl_s CReflect( const SchurCtrl<float>& ctrl )
{
    ElSchurCtrl_s ctrlC;
    ctrlC.useSDC = ctrl.useSDC;
    ctrlC.qrCtrl = CReflect( ctrl.qrCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    return ctrlC;
}
inline ElSchurCtrl_d CReflect( const SchurCtrl<double>& ctrl )
{
    ElSchurCtrl_d ctrlC;
    ctrlC.useSDC = ctrl.useSDC;
    ctrlC.qrCtrl = CReflect( ctrl.qrCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    return ctrlC;
}

inline SchurCtrl<float> CReflect( const ElSchurCtrl_s& ctrlC )
{
    SchurCtrl<float> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.qrCtrl = CReflect( ctrlC.qrCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    return ctrl;
}
inline SchurCtrl<double> CReflect( const ElSchurCtrl_d& ctrlC )
{
    SchurCtrl<double> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.qrCtrl = CReflect( ctrlC.qrCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    return ctrl;
}

// Factorizations
// ^^^^^^^^^^^^^^
inline ElLDLPivotType CReflect( LDLPivotType pivotType )
{ return static_cast<ElLDLPivotType>( pivotType ); }

inline LDLPivotType CReflect( ElLDLPivotType pivotType )
{ return static_cast<LDLPivotType>( pivotType ); }

inline ElLDLPivotCtrl_s CReflect( const LDLPivotCtrl<float>& ctrl )
{
    ElLDLPivotCtrl_s ctrlC;
    ctrlC.pivotType = CReflect(ctrl.pivotType);
    ctrlC.gamma = ctrl.gamma;
    return ctrlC;
}
inline ElLDLPivotCtrl_d CReflect( const LDLPivotCtrl<double>& ctrl )
{
    ElLDLPivotCtrl_d ctrlC;
    ctrlC.pivotType = CReflect(ctrl.pivotType);
    ctrlC.gamma = ctrl.gamma;
    return ctrlC;
}

inline LDLPivotCtrl<float> CReflect( ElLDLPivotCtrl_s ctrlC )
{
    LDLPivotCtrl<float> ctrl;
    ctrl.pivotType = CReflect(ctrlC.pivotType);
    ctrl.gamma = ctrlC.gamma;
    return ctrl;
}
inline LDLPivotCtrl<double> CReflect( ElLDLPivotCtrl_d ctrlC )
{
    LDLPivotCtrl<double> ctrl;
    ctrl.pivotType = CReflect(ctrlC.pivotType);
    ctrl.gamma = ctrlC.gamma;
    return ctrl;
}

inline ElLDLPivot CReflect( const LDLPivot& pivot )
{
    ElLDLPivot pivotC;
    pivotC.nb = pivot.nb;
    pivotC.from[0] = pivot.from[0];
    pivotC.from[1] = pivot.from[1];
    return pivotC;
}

inline LDLPivot CReflect( const ElLDLPivot& pivotC )
{
    LDLPivot pivot;
    pivot.nb = pivotC.nb;
    pivot.from[0] = pivotC.from[0];
    pivot.from[1] = pivotC.from[1];
    return pivot;
}

inline ElInertiaType CReflect( const InertiaType& inertia )
{ 
    ElInertiaType inertiaC;
    inertiaC.numPositive = inertia.numPositive;
    inertiaC.numNegative = inertia.numNegative;
    inertiaC.numZero = inertia.numZero;
    return inertiaC;
}

inline InertiaType CReflect( const ElInertiaType& inertiaC )
{ 
    InertiaType inertia;
    inertia.numPositive = inertiaC.numPositive;
    inertia.numNegative = inertiaC.numNegative;
    inertia.numZero = inertiaC.numZero;
    return inertia;
}

inline ElQRCtrl_s CReflect( const QRCtrl<float>& ctrl )
{ 
    ElQRCtrl_s ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    return ctrlC;
}
inline ElQRCtrl_d CReflect( const QRCtrl<double>& ctrl )
{ 
    ElQRCtrl_d ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    return ctrlC;
}

inline QRCtrl<float> CReflect( const ElQRCtrl_s& ctrlC )
{ 
    QRCtrl<float> ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    return ctrl;
}
inline QRCtrl<double> CReflect( const ElQRCtrl_d& ctrlC )
{ 
    QRCtrl<double> ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    return ctrl;
}

// Properties
// ^^^^^^^^^^
inline ElNormType CReflect( NormType type )
{ return static_cast<ElNormType>(type); }
inline NormType CReflect( ElNormType type )
{ return static_cast<NormType>(type); }

inline ElPseudospecNorm CReflect( PseudospecNorm psNorm )
{ return static_cast<ElPseudospecNorm>(psNorm); }
inline PseudospecNorm CReflect( ElPseudospecNorm psNorm )
{ return static_cast<PseudospecNorm>(psNorm); }

inline ElSnapshotCtrl CReflect( const SnapshotCtrl& ctrl )
{
    ElSnapshotCtrl ctrlC;
    ctrlC.realSize = ctrl.realSize; 
    ctrlC.imagSize = ctrl.imagSize;
    ctrlC.imgSaveFreq = ctrl.imgSaveFreq;
    ctrlC.numSaveFreq = ctrl.numSaveFreq;
    ctrlC.imgDispFreq = ctrl.imgDispFreq;
    ctrlC.imgSaveCount = ctrl.imgSaveCount;
    ctrlC.numSaveCount = ctrl.numSaveCount;
    ctrlC.imgDispCount = ctrl.imgDispCount;
    ctrlC.imgBase = CReflect(ctrl.imgBase);
    ctrlC.numBase = CReflect(ctrl.numBase);
    ctrlC.imgFormat = CReflect(ctrl.imgFormat);
    ctrlC.numFormat = CReflect(ctrl.numFormat);
    ctrlC.itCounts = ctrl.itCounts;
    return ctrlC;
}
inline SnapshotCtrl CReflect( const ElSnapshotCtrl& ctrlC )
{
    SnapshotCtrl ctrl;
    ctrl.realSize = ctrlC.realSize; 
    ctrl.imagSize = ctrlC.imagSize;
    ctrl.imgSaveFreq = ctrlC.imgSaveFreq;
    ctrl.numSaveFreq = ctrlC.numSaveFreq;
    ctrl.imgDispFreq = ctrlC.imgDispFreq;
    ctrl.imgSaveCount = ctrlC.imgSaveCount;
    ctrl.numSaveCount = ctrlC.numSaveCount;
    ctrl.imgDispCount = ctrlC.imgDispCount;
    ctrl.imgBase = CReflect(ctrlC.imgBase);
    ctrl.numBase = CReflect(ctrlC.numBase);
    ctrl.imgFormat = CReflect(ctrlC.imgFormat);
    ctrl.numFormat = CReflect(ctrlC.numFormat);
    ctrl.itCounts = ctrlC.itCounts;
    return ctrl;
}

inline ElPseudospecCtrl_s CReflect( const PseudospecCtrl<float>& ctrl )
{
    ElPseudospecCtrl_s ctrlC;
    ctrlC.norm = CReflect(ctrl.norm);
    ctrlC.blockWidth = ctrl.norm;
    ctrlC.schur = ctrl.schur;
    ctrlC.forceComplexSchur = ctrl.forceComplexSchur;
    ctrlC.forceComplexPs = ctrl.forceComplexPs;
    ctrlC.schurCtrl = CReflect(ctrl.schurCtrl);
    ctrlC.maxIts = ctrl.maxIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.deflate = ctrl.deflate;
    ctrlC.arnoldi = ctrl.arnoldi;
    ctrlC.basisSize = ctrl.basisSize;
    ctrlC.reorthog = ctrl.reorthog;
    ctrlC.progress = ctrl.progress;
    ctrlC.snapCtrl = CReflect(ctrl.snapCtrl);
    return ctrlC;
}
inline ElPseudospecCtrl_d CReflect( const PseudospecCtrl<double>& ctrl )
{
    ElPseudospecCtrl_d ctrlC;
    ctrlC.norm = CReflect(ctrl.norm);
    ctrlC.blockWidth = ctrl.norm;
    ctrlC.schur = ctrl.schur;
    ctrlC.forceComplexSchur = ctrl.forceComplexSchur;
    ctrlC.forceComplexPs = ctrl.forceComplexPs;
    ctrlC.schurCtrl = CReflect(ctrl.schurCtrl);
    ctrlC.maxIts = ctrl.maxIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.deflate = ctrl.deflate;
    ctrlC.arnoldi = ctrl.arnoldi;
    ctrlC.basisSize = ctrl.basisSize;
    ctrlC.reorthog = ctrl.reorthog;
    ctrlC.progress = ctrl.progress;
    ctrlC.snapCtrl = CReflect(ctrl.snapCtrl);
    return ctrlC;
}

inline PseudospecCtrl<float> CReflect( const ElPseudospecCtrl_s& ctrlC )
{
    PseudospecCtrl<float> ctrl;
    ctrl.norm = CReflect(ctrlC.norm);
    ctrl.blockWidth = ctrlC.norm;
    ctrl.schur = ctrlC.schur;
    ctrl.forceComplexSchur = ctrlC.forceComplexSchur;
    ctrl.forceComplexPs = ctrlC.forceComplexPs;
    ctrl.schurCtrl = CReflect(ctrlC.schurCtrl);
    ctrl.maxIts = ctrlC.maxIts;
    ctrl.tol = ctrlC.tol;
    ctrl.deflate = ctrlC.deflate;
    ctrl.arnoldi = ctrlC.arnoldi;
    ctrl.basisSize = ctrlC.basisSize;
    ctrl.reorthog = ctrlC.reorthog;
    ctrl.progress = ctrlC.progress;
    ctrl.snapCtrl = CReflect(ctrlC.snapCtrl);
    return ctrl;
}
inline PseudospecCtrl<double> CReflect( const ElPseudospecCtrl_d& ctrlC )
{
    PseudospecCtrl<double> ctrl;
    ctrl.norm = CReflect(ctrlC.norm);
    ctrl.blockWidth = ctrlC.norm;
    ctrl.schur = ctrlC.schur;
    ctrl.forceComplexSchur = ctrlC.forceComplexSchur;
    ctrl.forceComplexPs = ctrlC.forceComplexPs;
    ctrl.schurCtrl = CReflect(ctrlC.schurCtrl);
    ctrl.maxIts = ctrlC.maxIts;
    ctrl.tol = ctrlC.tol;
    ctrl.deflate = ctrlC.deflate;
    ctrl.arnoldi = ctrlC.arnoldi;
    ctrl.basisSize = ctrlC.basisSize;
    ctrl.reorthog = ctrlC.reorthog;
    ctrl.progress = ctrlC.progress;
    ctrl.snapCtrl = CReflect(ctrlC.snapCtrl);
    return ctrl;
}

inline ElSpectralBox_s CReflect( const SpectralBox<float>& box )
{
    ElSpectralBox_s boxC;
    boxC.center = CReflect(box.center);
    boxC.realWidth = box.realWidth;
    boxC.imagWidth = box.imagWidth;
    return boxC;
}

inline ElSpectralBox_d CReflect( const SpectralBox<double>& box )
{
    ElSpectralBox_d boxC;
    boxC.center = CReflect(box.center);
    boxC.realWidth = box.realWidth;
    boxC.imagWidth = box.imagWidth;
    return boxC;
}

inline SpectralBox<float> CReflect( const ElSpectralBox_s boxC )
{
    SpectralBox<float> box;
    box.center = CReflect(boxC.center);
    box.realWidth = CReflect(boxC.realWidth);
    box.imagWidth = CReflect(boxC.imagWidth);
    return box;
}

inline SpectralBox<double> CReflect( const ElSpectralBox_d boxC )
{
    SpectralBox<double> box;
    box.center = CReflect(boxC.center);
    box.realWidth = CReflect(boxC.realWidth);
    box.imagWidth = CReflect(boxC.imagWidth);
    return box;
}

// Solvers
// ^^^^^^^
inline ElTikhonovAlg CReflect( TikhonovAlg alg )
{ return static_cast<ElTikhonovAlg>(alg); }
inline TikhonovAlg CReflect( ElTikhonovAlg alg )
{ return static_cast<TikhonovAlg>(alg); }

inline ElRidgeAlg CReflect( RidgeAlg alg )
{ return static_cast<ElRidgeAlg>(alg); }
inline RidgeAlg CReflect( ElRidgeAlg alg )
{ return static_cast<RidgeAlg>(alg); }

// Optimization
// ------------
inline ElRegularization CReflect( Regularization penalty )
{ return static_cast<ElRegularization>(penalty); }
inline Regularization CReflect( ElRegularization penalty )
{ return static_cast<Regularization>(penalty); }

/* Linear programs
   ^^^^^^^^^^^^^^^ */
inline ElLPApproach CReflect( LPApproach approach )
{ return static_cast<ElLPApproach>(approach); }
inline LPApproach CReflect( ElLPApproach approach )
{ return static_cast<LPApproach>(approach); }

inline ElLPIPFLineSearchCtrl_s CReflect
( const lp::IPFLineSearchCtrl<float>& ctrl )
{
    ElLPIPFLineSearchCtrl_s ctrlC;
    ctrlC.gamma     = ctrl.gamma;
    ctrlC.beta      = ctrl.beta;
    ctrlC.psi       = ctrl.psi;
    ctrlC.stepRatio = ctrl.stepRatio;
    ctrlC.print     = ctrl.print;
    return ctrlC;
}
inline ElLPIPFLineSearchCtrl_d CReflect
( const lp::IPFLineSearchCtrl<double>& ctrl )
{
    ElLPIPFLineSearchCtrl_d ctrlC;
    ctrlC.gamma     = ctrl.gamma;
    ctrlC.beta      = ctrl.beta;
    ctrlC.psi       = ctrl.psi;
    ctrlC.stepRatio = ctrl.stepRatio;
    ctrlC.print     = ctrl.print;
    return ctrlC;
}
inline lp::IPFLineSearchCtrl<float> CReflect
( ElLPIPFLineSearchCtrl_s ctrlC )
{
    lp::IPFLineSearchCtrl<float> ctrl;
    ctrl.gamma     = ctrlC.gamma;
    ctrl.beta      = ctrlC.beta;
    ctrl.psi       = ctrlC.psi;
    ctrl.stepRatio = ctrlC.stepRatio;
    ctrl.print     = ctrlC.print;
    return ctrl;
}
inline lp::IPFLineSearchCtrl<double> CReflect
( ElLPIPFLineSearchCtrl_d ctrlC )
{
    lp::IPFLineSearchCtrl<double> ctrl;
    ctrl.gamma     = ctrlC.gamma;
    ctrl.beta      = ctrlC.beta;
    ctrl.psi       = ctrlC.psi;
    ctrl.stepRatio = ctrlC.stepRatio;
    ctrl.print     = ctrlC.print;
    return ctrl;
}

/* Direct conic form
   """"""""""""""""" */
inline ElLPDirectKKTSystem CReflect( lp::direct::KKTSystem system )
{ return static_cast<ElLPDirectKKTSystem>(system); }
inline lp::direct::KKTSystem CReflect( ElLPDirectKKTSystem system )
{ return static_cast<lp::direct::KKTSystem>(system); }

inline ElLPDirectADMMCtrl_s CReflect( const lp::direct::ADMMCtrl<float>& ctrl )
{
    ElLPDirectADMMCtrl_s ctrlC;
    ctrlC.rho     = ctrl.rho;
    ctrlC.alpha   = ctrl.alpha;
    ctrlC.maxIter = ctrl.maxIter;
    ctrlC.absTol  = ctrl.absTol;
    ctrlC.relTol  = ctrl.relTol;
    ctrlC.inv     = ctrl.inv;
    ctrlC.print   = ctrl.print;
    return ctrlC;
}
inline ElLPDirectADMMCtrl_d CReflect( const lp::direct::ADMMCtrl<double>& ctrl )
{
    ElLPDirectADMMCtrl_d ctrlC;
    ctrlC.rho     = ctrl.rho;
    ctrlC.alpha   = ctrl.alpha;
    ctrlC.maxIter = ctrl.maxIter;
    ctrlC.absTol  = ctrl.absTol;
    ctrlC.relTol  = ctrl.relTol;
    ctrlC.inv     = ctrl.inv;
    ctrlC.print   = ctrl.print;
    return ctrlC;
}
inline lp::direct::ADMMCtrl<float> CReflect( ElLPDirectADMMCtrl_s ctrlC )
{
    lp::direct::ADMMCtrl<float> ctrl;
    ctrl.rho     = ctrlC.rho;
    ctrl.alpha   = ctrlC.alpha;
    ctrl.maxIter = ctrlC.maxIter;
    ctrl.absTol  = ctrlC.absTol;
    ctrl.relTol  = ctrlC.relTol;
    ctrl.inv     = ctrlC.inv;
    ctrl.print   = ctrlC.print;
    return ctrl;
}
inline lp::direct::ADMMCtrl<double> CReflect( ElLPDirectADMMCtrl_d ctrlC )
{
    lp::direct::ADMMCtrl<double> ctrl;
    ctrl.rho     = ctrlC.rho;
    ctrl.alpha   = ctrlC.alpha;
    ctrl.maxIter = ctrlC.maxIter;
    ctrl.absTol  = ctrlC.absTol;
    ctrl.relTol  = ctrlC.relTol;
    ctrl.inv     = ctrlC.inv;
    ctrl.print   = ctrlC.print;
    return ctrl;
}

inline ElLPDirectIPFCtrl_s CReflect( const lp::direct::IPFCtrl<float>& ctrl )
{
    ElLPDirectIPFCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElLPDirectIPFCtrl_d CReflect( const lp::direct::IPFCtrl<double>& ctrl )
{
    ElLPDirectIPFCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline lp::direct::IPFCtrl<float> CReflect( ElLPDirectIPFCtrl_s ctrlC )
{
    lp::direct::IPFCtrl<float> ctrl(false);
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.system            = CReflect(ctrlC.system);
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}
inline lp::direct::IPFCtrl<double> CReflect( ElLPDirectIPFCtrl_d ctrlC )
{
    lp::direct::IPFCtrl<double> ctrl(false);
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.system            = CReflect(ctrlC.system);
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}

inline ElLPDirectMehrotraCtrl_s CReflect
( const lp::direct::MehrotraCtrl<float>& ctrl )
{
    ElLPDirectMehrotraCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElLPDirectMehrotraCtrl_d CReflect
( const lp::direct::MehrotraCtrl<double>& ctrl )
{
    ElLPDirectMehrotraCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline lp::direct::MehrotraCtrl<float> CReflect
( ElLPDirectMehrotraCtrl_s ctrlC )
{
    lp::direct::MehrotraCtrl<float> ctrl(false);
    ctrl.primalInitialized  = ctrlC.primalInitialized;
    ctrl.dualInitialized    = ctrlC.dualInitialized;
    ctrl.tol                = ctrlC.tol;
    ctrl.maxIts             = ctrlC.maxIts;
    ctrl.maxStepRatio       = ctrlC.maxStepRatio;
    ctrl.system             = CReflect(ctrlC.system);
    ctrl.print              = ctrlC.print;
    return ctrl;
}
inline lp::direct::MehrotraCtrl<double> CReflect
( ElLPDirectMehrotraCtrl_d ctrlC )
{
    lp::direct::MehrotraCtrl<double> ctrl(false);
    ctrl.primalInitialized  = ctrlC.primalInitialized;
    ctrl.dualInitialized    = ctrlC.dualInitialized;
    ctrl.tol                = ctrlC.tol;
    ctrl.maxIts             = ctrlC.maxIts;
    ctrl.maxStepRatio       = ctrlC.maxStepRatio;
    ctrl.system             = CReflect(ctrlC.system);
    ctrl.print              = ctrlC.print;
    return ctrl;
}

inline ElLPDirectCtrl_s CReflect( const lp::direct::Ctrl<float>& ctrl )
{
    ElLPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.admmCtrl     = CReflect(ctrl.admmCtrl);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElLPDirectCtrl_d CReflect( const lp::direct::Ctrl<double>& ctrl )
{
    ElLPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.admmCtrl     = CReflect(ctrl.admmCtrl);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline lp::direct::Ctrl<float> CReflect( ElLPDirectCtrl_s ctrlC )
{
    lp::direct::Ctrl<float> ctrl(false);
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.admmCtrl     = CReflect(ctrlC.admmCtrl);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline lp::direct::Ctrl<double> CReflect( ElLPDirectCtrl_d ctrlC )
{
    lp::direct::Ctrl<double> ctrl(false);
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.admmCtrl     = CReflect(ctrlC.admmCtrl);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElLPAffineIPFCtrl_s CReflect( const lp::affine::IPFCtrl<float>& ctrl )
{
    ElLPAffineIPFCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElLPAffineIPFCtrl_d CReflect( const lp::affine::IPFCtrl<double>& ctrl )
{
    ElLPAffineIPFCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline lp::affine::IPFCtrl<float> CReflect( ElLPAffineIPFCtrl_s ctrlC )
{
    lp::affine::IPFCtrl<float> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}
inline lp::affine::IPFCtrl<double> CReflect( ElLPAffineIPFCtrl_d ctrlC )
{
    lp::affine::IPFCtrl<double> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}

inline ElLPAffineMehrotraCtrl_s CReflect
( const lp::affine::MehrotraCtrl<float>& ctrl )
{
    ElLPAffineMehrotraCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElLPAffineMehrotraCtrl_d CReflect
( const lp::affine::MehrotraCtrl<double>& ctrl )
{
    ElLPAffineMehrotraCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline lp::affine::MehrotraCtrl<float> CReflect
( ElLPAffineMehrotraCtrl_s ctrlC )
{
    lp::affine::MehrotraCtrl<float> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.maxStepRatio      = ctrlC.maxStepRatio;
    ctrl.print             = ctrlC.print;
    return ctrl;
}
inline lp::affine::MehrotraCtrl<double> CReflect
( ElLPAffineMehrotraCtrl_d ctrlC )
{
    lp::affine::MehrotraCtrl<double> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.maxStepRatio      = ctrlC.maxStepRatio;
    ctrl.print             = ctrlC.print;
    return ctrl;
}

inline ElLPAffineCtrl_s CReflect( const lp::affine::Ctrl<float>& ctrl )
{
    ElLPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElLPAffineCtrl_d CReflect( const lp::affine::Ctrl<double>& ctrl )
{
    ElLPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline lp::affine::Ctrl<float> CReflect( ElLPAffineCtrl_s ctrlC )
{
    lp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline lp::affine::Ctrl<double> CReflect( ElLPAffineCtrl_d ctrlC )
{
    lp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Quadratic programs
   ^^^^^^^^^^^^^^^^^^ */
inline ElQPApproach CReflect( QPApproach approach )
{ return static_cast<ElQPApproach>(approach); }
inline QPApproach CReflect( ElQPApproach approach )
{ return static_cast<QPApproach>(approach); }

inline ElQPIPFLineSearchCtrl_s CReflect
( const qp::IPFLineSearchCtrl<float>& ctrl )
{
    ElQPIPFLineSearchCtrl_s ctrlC;
    ctrlC.gamma     = ctrl.gamma;
    ctrlC.beta      = ctrl.beta;
    ctrlC.psi       = ctrl.psi;
    ctrlC.stepRatio = ctrl.stepRatio;
    ctrlC.print     = ctrl.print;
    return ctrlC;
}
inline ElQPIPFLineSearchCtrl_d CReflect
( const qp::IPFLineSearchCtrl<double>& ctrl )
{
    ElQPIPFLineSearchCtrl_d ctrlC;
    ctrlC.gamma     = ctrl.gamma;
    ctrlC.beta      = ctrl.beta;
    ctrlC.psi       = ctrl.psi;
    ctrlC.stepRatio = ctrl.stepRatio;
    ctrlC.print     = ctrl.print;
    return ctrlC;
}
inline qp::IPFLineSearchCtrl<float> CReflect
( ElQPIPFLineSearchCtrl_s ctrlC )
{
    qp::IPFLineSearchCtrl<float> ctrl;
    ctrl.gamma     = ctrlC.gamma;
    ctrl.beta      = ctrlC.beta;
    ctrl.psi       = ctrlC.psi;
    ctrl.stepRatio = ctrlC.stepRatio;
    ctrl.print     = ctrlC.print;
    return ctrl;
}
inline qp::IPFLineSearchCtrl<double> CReflect
( ElQPIPFLineSearchCtrl_d ctrlC )
{
    qp::IPFLineSearchCtrl<double> ctrl;
    ctrl.gamma     = ctrlC.gamma;
    ctrl.beta      = ctrlC.beta;
    ctrl.psi       = ctrlC.psi;
    ctrl.stepRatio = ctrlC.stepRatio;
    ctrl.print     = ctrlC.print;
    return ctrl;
}

/* Direct conic form
   """"""""""""""""" */
inline ElQPDirectKKTSystem CReflect( qp::direct::KKTSystem system )
{ return static_cast<ElQPDirectKKTSystem>(system); }
inline qp::direct::KKTSystem CReflect( ElQPDirectKKTSystem system )
{ return static_cast<qp::direct::KKTSystem>(system); }

inline ElQPDirectIPFCtrl_s CReflect( const qp::direct::IPFCtrl<float>& ctrl )
{
    ElQPDirectIPFCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElQPDirectIPFCtrl_d CReflect( const qp::direct::IPFCtrl<double>& ctrl )
{
    ElQPDirectIPFCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline qp::direct::IPFCtrl<float> CReflect( ElQPDirectIPFCtrl_s ctrlC )
{
    qp::direct::IPFCtrl<float> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.system            = CReflect(ctrlC.system);
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}
inline qp::direct::IPFCtrl<double> CReflect( ElQPDirectIPFCtrl_d ctrlC )
{
    qp::direct::IPFCtrl<double> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.system            = CReflect(ctrlC.system);
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}

inline ElQPDirectMehrotraCtrl_s CReflect
( const qp::direct::MehrotraCtrl<float>& ctrl )
{
    ElQPDirectMehrotraCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElQPDirectMehrotraCtrl_d CReflect
( const qp::direct::MehrotraCtrl<double>& ctrl )
{
    ElQPDirectMehrotraCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.system            = CReflect(ctrl.system);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline qp::direct::MehrotraCtrl<float> CReflect
( ElQPDirectMehrotraCtrl_s ctrlC )
{
    qp::direct::MehrotraCtrl<float> ctrl;
    ctrl.primalInitialized  = ctrlC.primalInitialized;
    ctrl.dualInitialized    = ctrlC.dualInitialized;
    ctrl.tol                = ctrlC.tol;
    ctrl.maxIts             = ctrlC.maxIts;
    ctrl.maxStepRatio       = ctrlC.maxStepRatio;
    ctrl.system             = CReflect(ctrlC.system);
    ctrl.print              = ctrlC.print;
    return ctrl;
}
inline qp::direct::MehrotraCtrl<double> CReflect
( ElQPDirectMehrotraCtrl_d ctrlC )
{
    qp::direct::MehrotraCtrl<double> ctrl;
    ctrl.primalInitialized  = ctrlC.primalInitialized;
    ctrl.dualInitialized    = ctrlC.dualInitialized;
    ctrl.tol                = ctrlC.tol;
    ctrl.maxIts             = ctrlC.maxIts;
    ctrl.maxStepRatio       = ctrlC.maxStepRatio;
    ctrl.system             = CReflect(ctrlC.system);
    ctrl.print              = ctrlC.print;
    return ctrl;
}

inline ElQPDirectCtrl_s CReflect( const qp::direct::Ctrl<float>& ctrl )
{
    ElQPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElQPDirectCtrl_d CReflect( const qp::direct::Ctrl<double>& ctrl )
{
    ElQPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline qp::direct::Ctrl<float> CReflect( ElQPDirectCtrl_s ctrlC )
{
    qp::direct::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline qp::direct::Ctrl<double> CReflect( ElQPDirectCtrl_d ctrlC )
{
    qp::direct::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElQPAffineIPFCtrl_s CReflect( const qp::affine::IPFCtrl<float>& ctrl )
{
    ElQPAffineIPFCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElQPAffineIPFCtrl_d CReflect( const qp::affine::IPFCtrl<double>& ctrl )
{
    ElQPAffineIPFCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.centering         = ctrl.centering;
    ctrlC.lineSearchCtrl    = CReflect(ctrl.lineSearchCtrl);
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline qp::affine::IPFCtrl<float> CReflect( ElQPAffineIPFCtrl_s ctrlC )
{
    qp::affine::IPFCtrl<float> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}
inline qp::affine::IPFCtrl<double> CReflect( ElQPAffineIPFCtrl_d ctrlC )
{
    qp::affine::IPFCtrl<double> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.centering         = ctrlC.centering;
    ctrl.lineSearchCtrl    = CReflect(ctrlC.lineSearchCtrl);
    ctrl.print             = ctrlC.print;
    return ctrl;
}

inline ElQPAffineMehrotraCtrl_s CReflect
( const qp::affine::MehrotraCtrl<float>& ctrl )
{
    ElQPAffineMehrotraCtrl_s ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline ElQPAffineMehrotraCtrl_d CReflect
( const qp::affine::MehrotraCtrl<double>& ctrl )
{
    ElQPAffineMehrotraCtrl_d ctrlC;
    ctrlC.primalInitialized = ctrl.primalInitialized;
    ctrlC.dualInitialized   = ctrl.dualInitialized;
    ctrlC.tol               = ctrl.tol;
    ctrlC.maxIts            = ctrl.maxIts;
    ctrlC.maxStepRatio      = ctrl.maxStepRatio;
    ctrlC.print             = ctrl.print;
    return ctrlC;
}
inline qp::affine::MehrotraCtrl<float> CReflect
( ElQPAffineMehrotraCtrl_s ctrlC )
{
    qp::affine::MehrotraCtrl<float> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.maxStepRatio      = ctrlC.maxStepRatio;
    ctrl.print             = ctrlC.print;
    return ctrl;
}
inline qp::affine::MehrotraCtrl<double> CReflect
( ElQPAffineMehrotraCtrl_d ctrlC )
{
    qp::affine::MehrotraCtrl<double> ctrl;
    ctrl.primalInitialized = ctrlC.primalInitialized;
    ctrl.dualInitialized   = ctrlC.dualInitialized;
    ctrl.tol               = ctrlC.tol;
    ctrl.maxIts            = ctrlC.maxIts;
    ctrl.maxStepRatio      = ctrlC.maxStepRatio;
    ctrl.print             = ctrlC.print;
    return ctrl;
}

inline ElQPAffineCtrl_s CReflect( const qp::affine::Ctrl<float>& ctrl )
{
    ElQPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElQPAffineCtrl_d CReflect( const qp::affine::Ctrl<double>& ctrl )
{
    ElQPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.ipfCtrl      = CReflect(ctrl.ipfCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline qp::affine::Ctrl<float> CReflect( ElQPAffineCtrl_s ctrlC )
{
    qp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline qp::affine::Ctrl<double> CReflect( ElQPAffineCtrl_d ctrlC )
{
    qp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.ipfCtrl      = CReflect(ctrlC.ipfCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

} // namespace El

#endif // ifndef EL_CREFLECT_C_HPP
