#pragma once
#ifndef EL_GLOBALARRAYS_C_H
#define EL_GLOBALARRAYS_C_H

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY) && defined(EL_ENABLE_RMA_GLOBAL_ARRAYS)
#ifdef __cplusplus
extern "C" {
#endif
/* An anonymous struct meant as a placeholder for GlobalArrays<T>
-------------------------------------------------------- */
typedef struct ElGlobalArrays_iDummy* ElGlobalArrays_i;
typedef struct ElGlobalArrays_sDummy* ElGlobalArrays_s;
typedef struct ElGlobalArrays_dDummy* ElGlobalArrays_d;

/* GlobalArrays<T>::GlobalArrays()
   ------------------------------- */
EL_EXPORT ElError ElGlobalArraysConstruct_i( ElGlobalArrays_i* A );
EL_EXPORT ElError ElGlobalArraysConstruct_s( ElGlobalArrays_s* A );
EL_EXPORT ElError ElGlobalArraysConstruct_d( ElGlobalArrays_d* A );

/* GlobalArrays<T>::~GlobalArrays() 
   -------------------------------- */
EL_EXPORT ElError ElGlobalArraysDestruct_i( ElGlobalArrays_i A );
EL_EXPORT ElError ElGlobalArraysDestruct_s( ElGlobalArrays_s A );
EL_EXPORT ElError ElGlobalArraysDestruct_d( ElGlobalArrays_d A );

/* int GlobalArrays< T >::GA_Create(int type, int ndim, int dims[], const char *array_name); 
   ---------------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysCreate_i ( ElGlobalArrays_i A, ElInt type, ElInt ndim, 
	ElInt dims[], const char* array_name, ElInt* g_a );
EL_EXPORT ElError ElGlobalArraysCreate_s ( ElGlobalArrays_s A, ElInt type, ElInt ndim, 
	ElInt dims[], const char* array_name, ElInt* g_a );
EL_EXPORT ElError ElGlobalArraysCreate_d ( ElGlobalArrays_d A, ElInt type, ElInt ndim, 
	ElInt dims[], const char* array_name, ElInt* g_a );

/* void GlobalArrays<T>::GA_Copy(int g_a, int g_b); 
 ------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysCopy_i( ElGlobalArrays_i A, ElInt g_a, ElInt g_b );
EL_EXPORT ElError ElGlobalArraysCopy_s( ElGlobalArrays_s A, ElInt g_a, ElInt g_b );
EL_EXPORT ElError ElGlobalArraysCopy_d( ElGlobalArrays_d A, ElInt g_a, ElInt g_b );

/* void GlobalArrays<T>::GA_Print(int g_a); 
 -----------------------------------------------------*/
EL_EXPORT ElError ElGlobalArraysPrint_i( ElGlobalArrays_i A, ElInt g_a );
EL_EXPORT ElError ElGlobalArraysPrint_s( ElGlobalArrays_s A, ElInt g_a );
EL_EXPORT ElError ElGlobalArraysPrint_d( ElGlobalArrays_d A, ElInt g_a );

/* void GlobalArrays<T>::GA_Symmetrize(int g_a); 
 -----------------------------------------------------*/
EL_EXPORT ElError ElGlobalArraysSymmetrize_i( ElGlobalArrays_i A, ElInt g_a );
EL_EXPORT ElError ElGlobalArraysSymmetrize_s( ElGlobalArrays_s A, ElInt g_a );
EL_EXPORT ElError ElGlobalArraysSymmetrize_d( ElGlobalArrays_d A, ElInt g_a );

/* void GlobalArrays<T>::GA_Destroy(int g_a);
   ------------------------------------------------*/
EL_EXPORT ElError ElGlobalArraysDestroy_i( ElGlobalArrays_i A, ElInt g_a );
EL_EXPORT ElError ElGlobalArraysDestroy_s( ElGlobalArrays_s A, ElInt g_a );
EL_EXPORT ElError ElGlobalArraysDestroy_d( ElGlobalArrays_d A, ElInt g_a );

/* void GlobalArrays<T>::GA_Add(void *alpha, int g_a, void* beta, int g_b, int g_c); 
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysAdd_i( ElGlobalArrays_i A, void* alpha, ElInt g_a, void* beta, ElInt g_b, ElInt g_c );
EL_EXPORT ElError ElGlobalArraysAdd_s( ElGlobalArrays_s A, void* alpha, ElInt g_a, void* beta, ElInt g_b, ElInt g_c );
EL_EXPORT ElError ElGlobalArraysAdd_d( ElGlobalArrays_d A, void* alpha, ElInt g_a, void* beta, ElInt g_b, ElInt g_c );

/* void GlobalArrays<T>::GA_Dgemm(char ta, char tb, int m, int n, int k, double alpha, int g_a, int g_b, double beta, int g_c );
   ----------------------------------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysDgemm_i( ElGlobalArrays_i A, char ta, char tb, ElInt m, ElInt n, ElInt k, 
	double alpha, ElInt g_a, ElInt g_b, double beta, ElInt g_c );
EL_EXPORT ElError ElGlobalArraysDgemm_s( ElGlobalArrays_s A, char ta, char tb, ElInt m, ElInt n, ElInt k, 
	double alpha, ElInt g_a, ElInt g_b, double beta, ElInt g_c );
EL_EXPORT ElError ElGlobalArraysDgemm_d( ElGlobalArrays_d A, char ta, char tb, ElInt m, ElInt n, ElInt k, 
	double alpha, ElInt g_a, ElInt g_b, double beta, ElInt g_c );

/* int GlobalArrays<T>::GA_Duplicate(int g_a, const char *array_name);
 ------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysDuplicate_i( ElGlobalArrays_i A, ElInt g_a, 
	const char* array_name, ElInt* g_dup );
EL_EXPORT ElError ElGlobalArraysDuplicate_s( ElGlobalArrays_s A, ElInt g_a, 
	const char* array_name, ElInt* g_dup );
EL_EXPORT ElError ElGlobalArraysDuplicate_d( ElGlobalArrays_d A, ElInt g_a, 
	const char* array_name, ElInt* g_dup );

/* void GlobalArrays<T>::GA_Fill(int g_a, void *value);
 -------------------------------------- */
EL_EXPORT ElError ElGlobalArraysFill_i( ElGlobalArrays_i A, ElInt g_a, void* value );
EL_EXPORT ElError ElGlobalArraysFill_s( ElGlobalArrays_s A, ElInt g_a, void* value );
EL_EXPORT ElError ElGlobalArraysFill_d( ElGlobalArrays_d A, ElInt g_a, void* value );

/* void GlobalArrays<T>::GA_Initialize();
 ------------------------------------- */
EL_EXPORT ElError ElGlobalArraysInitialize_i( ElGlobalArrays_i A );
EL_EXPORT ElError ElGlobalArraysInitialize_s( ElGlobalArrays_s A );
EL_EXPORT ElError ElGlobalArraysInitialize_d( ElGlobalArrays_d A );

/* void GlobalArrays<T>::GA_Sync();
 --------------------------------------- */
EL_EXPORT ElError ElGlobalArraysSync_i( ElGlobalArrays_i A );
EL_EXPORT ElError ElGlobalArraysSync_s( ElGlobalArrays_s A );
EL_EXPORT ElError ElGlobalArraysSync_d( ElGlobalArrays_d A );

/* void GlobalArrays<T>::GA_Terminate();
 ---------------------------------------- */
EL_EXPORT ElError ElGlobalArraysTerminate_i( ElGlobalArrays_i A );
EL_EXPORT ElError ElGlobalArraysTerminate_s( ElGlobalArrays_s A );
EL_EXPORT ElError ElGlobalArraysTerminate_d( ElGlobalArrays_d A );

/* void GlobalArrays<T>::GA_Transpose(int g_a, int g_b);
 ------------------------------------------ */
EL_EXPORT ElError ElGlobalArraysTranspose_i( ElGlobalArrays_i A, ElInt g_a, ElInt g_b );
EL_EXPORT ElError ElGlobalArraysTranspose_s( ElGlobalArrays_s A, ElInt g_a, ElInt g_b );
EL_EXPORT ElError ElGlobalArraysTranspose_d( ElGlobalArrays_d A, ElInt g_a, ElInt g_b );

/* void GlobalArrays<T>::NGA_Access(int g_a, int lo[], int hi[], void *ptr, int ld[]);
 ---------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysAccess_i( ElGlobalArrays_i A, ElInt g_a, ElInt lo[], ElInt hi[], ElInt* ptr, ElInt ld[] );
EL_EXPORT ElError ElGlobalArraysAccess_s( ElGlobalArrays_s A, ElInt g_a, ElInt lo[], ElInt hi[], float* ptr, ElInt ld[] );
EL_EXPORT ElError ElGlobalArraysAccess_d( ElGlobalArrays_d A, ElInt g_a, ElInt lo[], ElInt hi[], double* ptr, ElInt ld[] );

/* void GlobalArrays<T>::NGA_Distribution(int g_a, int iproc, int lo[], int hi[]);
 ---------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysDistribution_i( ElGlobalArrays_i A, ElInt g_a, ElInt iproc, ElInt lo[], ElInt hi[] );
EL_EXPORT ElError ElGlobalArraysDistribution_s( ElGlobalArrays_s A, ElInt g_a, ElInt iproc, ElInt lo[], ElInt hi[] );
EL_EXPORT ElError ElGlobalArraysDistribution_d( ElGlobalArrays_d A, ElInt g_a, ElInt iproc, ElInt lo[], ElInt hi[] );

/* void GlobalArrays<T>::NGA_Acc(int g_a, int lo[], int hi[],void* buf,int ld[],void* alpha);
 ---------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysAccumulate_i( ElGlobalArrays_i A, ElInt g_a, ElInt lo[], 
	ElInt hi[], ElInt* ptr, ElInt ld[], void* alpha );
EL_EXPORT ElError ElGlobalArraysAccumulate_s( ElGlobalArrays_s A, ElInt g_a, ElInt lo[], 
	ElInt hi[], float* ptr, ElInt ld[], void* alpha );
EL_EXPORT ElError ElGlobalArraysAccumulate_d( ElGlobalArrays_d A, ElInt g_a, ElInt lo[], 
	ElInt hi[], double* ptr, ElInt ld[], void* alpha );

/* void GlobalArrays<T>::NGA_Get(int g_a, int lo[], int hi[], void* buf, int ld[]); 
 ------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysGet_i( ElGlobalArrays_i A, ElInt g_a, ElInt lo[], 
	ElInt hi[], ElInt* buf, ElInt ld[] );
EL_EXPORT ElError ElGlobalArraysGet_s( ElGlobalArrays_s A, ElInt g_a, ElInt lo[], 
	ElInt hi[], float* buf, ElInt ld[] );
EL_EXPORT ElError ElGlobalArraysGet_d( ElGlobalArrays_d A, ElInt g_a, ElInt lo[], 
	ElInt hi[], double* buf, ElInt ld[] );

/* void GlobalArrays<T>::NGA_NbAcc(int g_a,int lo[], int hi[],void* buf,int ld[],void* alpha, ga_nbhdl_t* nbhandle);
 ----------------------------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysNBAccumulate_i( ElGlobalArrays_i A, ElInt g_a, ElInt lo[], 
	ElInt hi[], ElInt* ptr, ElInt ld[], void* alpha, ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBAccumulate_s( ElGlobalArrays_s A, ElInt g_a, ElInt lo[], 
	ElInt hi[], float* ptr, ElInt ld[], void* alpha, ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBAccumulate_d( ElGlobalArrays_d A, ElInt g_a, ElInt lo[], 
	ElInt hi[], double* ptr, ElInt ld[], void* alpha, ElInt* nbhandle );

/* void GlobalArrays<T>::NGA_NbGet(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle);
 --------------------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysNBGet_i( ElGlobalArrays_i A, ElInt g_a, ElInt lo[], 
	ElInt hi[], ElInt* buf, ElInt ld[], ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBGet_s( ElGlobalArrays_s A, ElInt g_a, ElInt lo[], 
	ElInt hi[], float* buf, ElInt ld[], ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBGet_d( ElGlobalArrays_d A, ElInt g_a, ElInt lo[], 
	ElInt hi[], double* buf, ElInt ld[], ElInt* nbhandle );

/* void GlobalArrays<T>::NGA_NbPut(int g_a, int lo[], int hi[], void* ptr, int ld[], ga_nbhdl_t* nbhandle);
  --------------------------------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysNBPut_i( ElGlobalArrays_i A, ElInt g_a, ElInt lo[], 
	ElInt hi[], ElInt* ptr, ElInt ld[], ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBPut_s( ElGlobalArrays_s A, ElInt g_a, ElInt lo[], 
	ElInt hi[], float* ptr, ElInt ld[], ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBPut_d( ElGlobalArrays_d A, ElInt g_a, ElInt lo[], 
	ElInt hi[], double* ptr, ElInt ld[], ElInt* nbhandle );

/* int GlobalArrays<T>::NGA_NbTest(ga_nbhdl_t* nbhandle);
 ------------------------------------------------ */
EL_EXPORT ElError ElGlobalArraysNBTest_i( ElGlobalArrays_i A, ElInt* nbhandle, ElInt* status );
EL_EXPORT ElError ElGlobalArraysNBTest_s( ElGlobalArrays_s A, ElInt* nbhandle, ElInt* status );
EL_EXPORT ElError ElGlobalArraysNBTest_d( ElGlobalArrays_d A, ElInt* nbhandle, ElInt* status );

/* void GlobalArrays<T>::NGA_NbWait(ga_nbhdl_t* nbhandle);
 ------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysNBWait_i( ElGlobalArrays_i A, ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBWait_s( ElGlobalArrays_s A, ElInt* nbhandle );
EL_EXPORT ElError ElGlobalArraysNBWait_d( ElGlobalArrays_d A, ElInt* nbhandle );

/* void GlobalArrays<T>::NGA_Put(int g_a, int lo[], int hi[], void* ptr, int ld[]); 
  ---------------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysPut_i( ElGlobalArrays_i A, ElInt g_a, ElInt lo[], 
	ElInt hi[], ElInt* ptr, ElInt ld[] );
EL_EXPORT ElError ElGlobalArraysPut_s( ElGlobalArrays_s A, ElInt g_a, ElInt lo[], 
	ElInt hi[], float* ptr, ElInt ld[] );
EL_EXPORT ElError ElGlobalArraysPut_d( ElGlobalArrays_d A, ElInt g_a, ElInt lo[], 
	ElInt hi[], double* ptr, ElInt ld[] );

/* long GlobalArrays<T>::NGA_Read_inc(int g_a, int ndim, int subscript[], long inc);
 ----------------------------------------------------------------- */
EL_EXPORT ElError ElGlobalArraysReadIncrement_i( ElGlobalArrays_i A, ElInt g_a, ElInt ndim, 
	ElInt subscript[], ElInt inc, ElInt* prev );
EL_EXPORT ElError ElGlobalArraysReadIncrement_s( ElGlobalArrays_s A, ElInt g_a, ElInt ndim, 
	ElInt subscript[], ElInt inc, ElInt* prev );
EL_EXPORT ElError ElGlobalArraysReadIncrement_d( ElGlobalArrays_d A, ElInt g_a, ElInt ndim, 
	ElInt subscript[], ElInt inc, ElInt* prev );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY) && defined(EL_ENABLE_RMA_GLOBAL_ARRAYS) */
#endif /* ifndef EL_GLOBALARRAYS_C_H */
