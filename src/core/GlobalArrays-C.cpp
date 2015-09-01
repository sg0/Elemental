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

extern "C" {
#define GA_BASE(SIG,SIGBASE,T) \
  /* GlobalArrays<T>::GlobalArrays() */ \
  ElError ElGlobalArraysCreate_ ## SIG ( ElGlobalArrays_ ## SIG * A ) \
  { EL_TRY( *A = CReflect( new GlobalArrays<T> ) ) } \
  /* GlobalArrays<T>::~GlobalArrays() */ \
  ElError ElGlobalArraysDestruct_ ## SIG ( ElGlobalArrays_ ## SIG A ) \
  { EL_TRY( delete CReflect(A) ) } \
  
#define GA_OPS(SIG,SIGBASE,T) \
  /* int GlobalArrays<T>::GA_Create_handle() */ \
  ElError ElGlobalArraysCreateHandle_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt* g_a ) \
  { EL_TRY( *ga = CReflect(A)->GA_Create_handle() ) } \
  /* void GlobalArrays<T>::GA_Set_data (int g_a, int ndim, int dims[], int type) */ \
  ElError ElGlobalArraysSetData_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt ndim, ElInt dims[], ElInt type ) \
  { EL_TRY( CReflect(A)->GA_Set_data (g_a, ndim, dims, type) ) } \
  /* void GlobalArrays<T>::GA_Allocate (int g_a) */ \
  ElError ElGlobalArraysAllocate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a ) \
  { EL_TRY( CReflect(A)->GA_Allocate (g_a) ) } \
  /* void GlobalArrays<T>::GA_Copy(int g_a, int g_b) */ \
  ElError ElGlobalArraysCopy_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt g_b ) \
  { EL_TRY( CReflect(A)->GA_Copy (g_a, g_b) ) } \
  /* void GlobalArrays<T>::GA_Destroy(int g_a) */ \
  ElError ElGlobalArraysDestroy_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a ) \
  { EL_TRY( CReflect(A)->GA_Destroy (g_a) ) } \
  /* void GlobalArrays<T>::GA_Add(void *alpha, int g_a, void* beta, int g_b, int g_c) */ \
  ElError ElGlobalArraysAdd_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt* alpha, ElInt g_a, ElInt* beta, ElInt g_b, ElInt g_c ) \
  { EL_TRY( CReflect(A)->GA_Add(alpha, g_a, beta, g_b, g_c) ) } \
  /* void GlobalArrays<T>::GA_Dgemm(char ta, char tb, int m, int n, int k, double alpha, int g_a, int g_b, double beta, int g_c ) */ \
  ElError ElGlobalArraysDgemm_ ## SIG ( ElGlobalArrays_ ## SIG A, char ta, ElInt m, ElInt n, ElInt k, \
	double alpha, ElInt g_a, ElInt g_b, double beta, ElInt g_c ) \
  { EL_TRY( CReflect(A)->GA_Dgemm(ta, tb, m, n, k, alpha, g_a, g_b, beta, g_c ) ) } \
  /* int GlobalArrays<T>::GA_Duplicate(int g_a, char* array_name) */ \
  ElError ElGlobalArraysDuplicate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, char* array_name, ElInt* g_dup ) \
  { EL_TRY( *g_dup = CReflect(A)->GA_Duplicate(g_a, array_name) ) } \
  /* void GlobalArrays<T>::GA_Fill(int g_a, void *value) */ \
  ElError ElGlobalArraysFill_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt* value ) \
  { EL_TRY( CReflect(A)->GA_Fill(g_a, value) ) } \
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
  /* void GlobalArrays<T>::NGA_Access(int g_a, int lo[], int hi[], void *ptr, int ld[]) */ \
  ElError ElGlobalArraysAccess_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], ElInt hi[], ElInt* ptr, ElInt ld[] ) \
  { EL_TRY( CReflect(A)->NGA_Access(g_a, lo, hi, ptr, ld) ) } \
  /* void GlobalArrays<T>::NGA_Acc(int g_a, int lo[], int hi[],void* buf,int ld[],void* alpha) */ \
  ElError ElGlobalArraysAccumulate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], ElInt* ptr, ElInt ld[], ElInt* alpha ) \
  { EL_TRY( CReflect(A)->NGA_Acc(g_a, lo, hi, ptr, ld, alpha) ) } \
  /* void GlobalArrays<T>::NGA_Get(int g_a, int lo[], int hi[], void* buf, int ld[]) */ \
  ElError ElGlobalArraysGet_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], ElInt* ptr, ElInt ld[] ) \
  { EL_TRY( CReflect(A)->NGA_Acc(g_a, lo, hi, ptr, ld, alpha) ) } \
  /* void GlobalArrays<T>::NGA_NbAcc(int g_a,int lo[], int hi[],void* buf,int ld[],void* alpha, ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBAccumulate_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], ElInt* ptr, ElInt ld[], ElInt* alpha, ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbAcc(g_a, lo, hi, ptr, ld, alpha, nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_NbGet(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBGet_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], ElInt* ptr, ElInt ld[], ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbGet(g_a, lo, hi, buf, ld, nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_NbPut(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBPut_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	ElInt hi[], ElInt* ptr, ElInt ld[], ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbPut(g_a, lo, hi, ptr, ld, nbhandle) ) } \
  /* int GlobalArrays<T>::NGA_NbTest(ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBTest_ ## SIG A ( ElGlobalArrays_ ## SIG A, ElInt* nbhandle, ElInt* status ) \
  { EL_TRY( *status = CReflect(A)->NGA_NbTest(nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_NbWait(ga_nbhdl_t* nbhandle) */ \
  ElError ElGlobalArraysNBWait_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt* nbhandle ) \
  { EL_TRY( CReflect(A)->NGA_NbWait(nbhandle) ) } \
  /* void GlobalArrays<T>::NGA_Put(int g_a, int lo[], int hi[], void* buf, int ld[]) */ \
  ElError ElGlobalArraysPut_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, ElInt lo[], \
	  ElInt hi[], ElInt* ptr, ElInt ld[] ) \
  { EL_TRY( CReflect(A)->NGA_Put(g_a, lo, hi, ptr, ld) ) } \
  /* long GlobalArrays<T>::NGA_Read_inc(int g_a, int subscript[], long inc) */ \
  ElError ElGlobalArraysReadIncrement_ ## SIG ( ElGlobalArrays_ ## SIG A, ElInt g_a, \
	  ElInt subscript[], ElInt inc, ElInt* prev ) \
  { EL_TRY( *prev = CReflect(A)->NGA_Read_inc(g_a, subscript, inc) ) } \

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
}
