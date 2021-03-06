/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.
   Sayan Ghosh,
   Washington State University

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY) && defined(EL_ENABLE_RMA_GLOBAL_ARRAYS)
extern "C" {
#define GA_BASE(SIG,SIGBASE,T) \
  /* GlobalArrays<T>::GlobalArrays() */ \
  ElError ElGlobalArraysConstruct_ ## SIG ( ElGlobalArrays_ ## SIG * A ) \
  { EL_TRY( *A = CReflect( new GlobalArrays<T> ) ) } \
  /* GlobalArrays<T>::~GlobalArrays() */ \
  ElError ElGlobalArraysDestruct_ ## SIG ( ElGlobalArrays_ ## SIG A ) \
  { EL_TRY( delete CReflect(A) ) } \
  
#define GA_OPS(SIG,SIGBASE,T) \
  /* int GlobalArrays< T >::GA_Create(int ndim, int dims[], const char *array_name, int chunk[]) */ \
  ElError ElGlobalArraysCreate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt ndim, \
	  ElInt dims[], const char *array_name, ElInt chunk[], ElInt* g_a ) \
  { EL_TRY( *g_a = CReflect(A)->GA_Create (ndim, dims, array_name, chunk) ) } \
  /* Int  GA_Create_irreg(Int ndim, Int dims[], const char *array_name, Int block[], Int map[]); */ \
  ElError ElGlobalArraysCreateIrreg_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt ndim, \
	  ElInt dims[], const char *array_name, ElInt block[], ElInt map[], ElInt* g_a ) \
  { EL_TRY( *g_a = CReflect(A)->GA_Create_irreg (ndim, dims, array_name, block, map) ) } \
  /* void GlobalArrays<T>::GA_Copy(int g_a, int g_b) */ \
  ElError ElGlobalArraysCopy_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt g_b ) \
  { EL_TRY( CReflect(A)->GA_Copy (g_a, g_b) ) } \
  /* void GlobalArrays<T>::GA_Print(int g_a) */ \
  ElError ElGlobalArraysPrint_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a ) \
  { EL_TRY( CReflect(A)->GA_Print (g_a) ) } \
  /* void GlobalArrays<T>::GA_Symmetrize(int g_a) */ \
  ElError ElGlobalArraysSymmetrize_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a ) \
  { EL_TRY( CReflect(A)->GA_Symmetrize (g_a) ) } \
  /* void GlobalArrays<T>::GA_Destroy(int g_a) */ \
  ElError ElGlobalArraysDestroy_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a ) \
  { EL_TRY( CReflect(A)->GA_Destroy (g_a) ) } \
  /* void GlobalArrays<T>::GA_Add(void *alpha, int g_a, void* beta, int g_b, int g_c) */ \
  ElError ElGlobalArraysAdd_ ## SIG ( ElGlobalArrays_ ## SIG A, CREFLECT(T)* alpha, ElInt g_a, CREFLECT(T)* beta, ElInt g_b, ElInt g_c ) \
  { EL_TRY( CReflect(A)->GA_Add(CReflect(alpha), g_a, CReflect(beta), g_b, g_c) ) } \
  /* T GlobalArrays<T>::GA_Dot(int g_a, int g_b); */ \
  EL_EXPORT ElError ElGlobalArraysDot_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt g_b, CREFLECT(T)* dotproduct ) \
  { EL_TRY( *dotproduct = CReflect(A)->GA_Dot(g_a, g_b) ) } \
  /* void GlobalArrays<T>::GA_Dgemm(char ta, char tb, int m, int n, int k, T alpha, int g_a, int g_b, T beta, int g_c ) */ \
  ElError ElGlobalArraysDgemm_ ## SIG ( ElGlobalArrays_ ## SIG A, char ta, char tb, ElInt m, ElInt n, ElInt k, \
	CREFLECT(T) alpha, ElInt g_a, ElInt g_b, CREFLECT(T) beta, ElInt g_c ) \
  { EL_TRY( CReflect(A)->GA_Dgemm(ta, tb, m, n, k, alpha, g_a, g_b, beta, g_c ) ) } \
  /* void GA_Gop(T x[], Int n, char op); */ \
  EL_EXPORT ElError ElGlobalArraysOp_ ## SIG ( ElGlobalArrays_ ## SIG A,  CREFLECT(T)* x, ElInt n, char op ) \
  { EL_TRY( CReflect(A)->GA_Gop(CReflect(x), n, op) ) } \
  /* int GlobalArrays<T>::GA_Duplicate(int g_a, const char* array_name) */ \
  ElError ElGlobalArraysDuplicate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, \
	  const char *array_name, ElInt* g_dup ) \
  { EL_TRY( *g_dup = CReflect(A)->GA_Duplicate(g_a, array_name) ) } \
  /* void GlobalArrays<T>::GA_Fill(int g_a, void *value) */ \
  ElError ElGlobalArraysFill_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, CREFLECT(T)* value ) \
  { EL_TRY( CReflect(A)->GA_Fill(g_a, CReflect(value)) ) } \
  /* void GlobalArrays<T>::GA_Initialize() */ \
  ElError ElGlobalArraysInitialize_ ## SIG ( ElGlobalArrays_ ## SIG A ) \
  { EL_TRY( CReflect(A)->GA_Initialize() ) } \
  /* void GlobalArrays<T>::GA_Sync() */ \
  ElError ElGlobalArraysSync_ ## SIG ( ElGlobalArrays_ ## SIG A ) \
  { EL_TRY( CReflect(A)->GA_Sync() ) } \
  /* void GlobalArrays<T>::GA_Terminate() */ \
  ElError ElGlobalArraysTerminate_ ## SIG ( ElGlobalArrays_ ## SIG A ) \
  { EL_TRY( CReflect(A)->GA_Terminate() ) } \
  /* void GlobalArrays<T>::GA_Transpose(int g_a, int g_b) */ \
  ElError ElGlobalArraysTranspose_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt g_b ) \
  { EL_TRY( CReflect(A)->GA_Transpose(g_a, g_b) ) } \
  /* void GlobalArrays<T>::NGA_Access(int g_a, int lo[], int hi[], T** ptr, int ld[]) */ \
  ElError ElGlobalArraysAccess_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], ElInt hi[], CREFLECT(T)** ptr, ElInt ld[] ) \
  { EL_TRY( CReflect(A)->NGA_Access(g_a, lo, hi, ptr, ld) ) } \
  /* void GlobalArrays<T>::NGA_Release(int g_a, int lo[], int hi[]) */ \
  ElError ElGlobalArraysRelease_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], ElInt hi[] ) \
  { EL_TRY( CReflect(A)->NGA_Release(g_a, lo, hi) ) } \
  /* void GlobalArrays<T>::NGA_Release_update(int g_a, int lo[], int hi[]) */ \
  ElError ElGlobalArraysReleaseUpdate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], ElInt hi[] ) \
  { EL_TRY( CReflect(A)->NGA_Release_update(g_a, lo, hi) ) } \
  /* void GlobalArrays<T>::NGA_Distribution(int g_a, int iproc, int lo[], int hi[]) */ \
  ElError ElGlobalArraysDistribution_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt iproc, ElInt lo[], ElInt hi[] ) \
  { EL_TRY( CReflect(A)->NGA_Distribution(g_a, iproc, lo, hi) ) } \
  /* void NGA_Inquire(Int g_a, Int * ndim, Int dims[]) */ \
  ElError ElGlobalArraysInquire_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt * ndim, ElInt dims[] ) \
  { EL_TRY( CReflect(A)->NGA_Inquire(g_a, ndim, dims) ) } \
  /* void GlobalArrays<T>::NGA_Acc(int g_a, int lo[], int hi[],void* ptr,int ld[],void* alpha) */ \
  ElError ElGlobalArraysAccumulate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], CREFLECT(T)* ptr, ElInt ld[], CREFLECT(T)* alpha ) \
  { EL_TRY( CReflect(A)->NGA_Acc(g_a, lo, hi, CReflect(ptr), ld, CReflect(alpha)) ) } \
  /* void GlobalArrays<T>::NGA_Get(int g_a, int lo[], int hi[], void* buf, int ld[]) */ \
  ElError ElGlobalArraysGet_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], CREFLECT(T)* buf, ElInt ld[] ) \
  { EL_TRY( CReflect(A)->NGA_Get(g_a, lo, hi, CReflect(buf), ld) ) } \
  /* void GlobalArrays<T>::NGA_NbAcc(int g_a,int lo[], int hi[],void* ptr,int ld[],void* alpha, ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBAccumulate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], CREFLECT(T)* ptr, ElInt ld[], CREFLECT(T)* alpha, ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbAcc(g_a, lo, hi, CReflect(ptr), ld, CReflect(alpha), nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_NbGet(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBGet_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], CREFLECT(T)* buf, ElInt ld[], ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbGet(g_a, lo, hi, CReflect(buf), ld, nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_NbPut(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBPut_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], CREFLECT(T)* ptr, ElInt ld[], ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbPut(g_a, lo, hi, CReflect(ptr), ld, nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_NbWait(ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBWait_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbWait(nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_Put(int g_a, int lo[], int hi[], void* buf, int ld[]) */ \
  ElError ElGlobalArraysPut_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	  ElInt hi[], CREFLECT(T)* ptr, ElInt ld[] ) \
  { EL_TRY( CReflect(A)->NGA_Put(g_a, lo, hi, CReflect(ptr), ld) ) } \
  /* long GlobalArrays<T>::NGA_Read_inc(int g_a, int ndim, int subscript[], long inc) */ \
  ElError ElGlobalArraysReadIncrement_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, \
	  ElInt subscript[], CREFLECT(T) inc, CREFLECT(T)* prev ) \
  { EL_TRY( *prev = CReflect(A)->NGA_Read_inc( g_a, subscript, CReflect(inc) ) ) } \

#define C_PROTO(SIG,SIGBASE,T) \
  GA_BASE(SIG,SIGBASE,T) \
  GA_OPS(SIG,SIGBASE,T) \

#ifndef C_PROTO_INT
# define C_PROTO_INT(SIG,T) C_PROTO(SIG,SIG,T)
#endif

#ifndef C_PROTO_REAL 
# define C_PROTO_REAL(SIG,T) C_PROTO(SIG,SIG,T)
#endif
#ifndef C_PROTO_FLOAT
# define C_PROTO_FLOAT C_PROTO_REAL(s,float)
#endif
#ifndef C_PROTO_DOUBLE
# define C_PROTO_DOUBLE C_PROTO_REAL(d,double)
#endif

#ifndef EL_NO_INT_PROTO
C_PROTO_INT(i,Int)
#endif

#ifndef EL_NO_REAL_PROTO
# if !defined(EL_NO_FLOAT_PROTO)
C_PROTO_FLOAT
# endif
# if !defined(EL_NO_DOUBLE_PROTO)
C_PROTO_DOUBLE
# endif
#endif
}
#endif // end of MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY) && defined(EL_ENABLE_RMA_GLOBAL_ARRAYS)     
