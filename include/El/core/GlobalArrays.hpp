#pragma once
#ifndef EL_GLOBALARRAYS_HPP
#define EL_GLOBALARRAYS_HPP

namespace El {

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY) && defined(EL_ENABLE_RMA_GLOBAL_ARRAYS)
template<typename T>
class GlobalArrays 
{
	public:
		GlobalArrays();
		~GlobalArrays();

		typedef Int ga_nbhdl_t;

		// Subset of GA C API
		// TODO follow either NGA_ or GA_ in function names
		Int  GA_Create(Int type, Int ndim, Int dims[], const char *array_name);
		Int  GA_Duplicate(Int g_a, const char *array_name);
		/* Create rmainterface object, and attach DM corresponding to g_a */
		Int  GA_Allocate(Int g_a);
		void GA_Copy(Int g_a, Int g_b); 
		/* At present, this is similar to Print (DistMatrix<T> ...) */
		void GA_Print_distribution(Int g_a);
		/* call Detach on rmainterface, empty DM corresponding to g_a and erase entry from ga_handles */
		void GA_Destroy(Int g_a);
		void GA_Add(void *alpha, Int g_a, void* beta, Int g_b, Int g_c); 
		void GA_Dgemm(char ta, char tb, Int m, Int n, Int k, double alpha, Int g_a, Int g_b, double beta, Int g_c );
		void GA_Fill(Int g_a, void *value);
		void GA_Initialize();
		void GA_Sync();
		/* Check if Detach was called already on rmainterface objects, and erase ga_handles */
		void GA_Terminate();
		void GA_Transpose(Int g_a, Int g_b);
		void NGA_Access(Int g_a, Int lo[], Int hi[], void *ptr, Int ld[]);
		void NGA_Acc(Int g_a, Int lo[], Int hi[],void* buf,Int ld[],void* alpha);
		void NGA_Get(Int g_a, Int lo[], Int hi[], void* buf, Int ld[]); 
		void NGA_NbAcc(Int g_a,Int lo[], Int hi[],void* buf,Int ld[],void* alpha, ga_nbhdl_t* nbhandle);
		void NGA_NbGet(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], ga_nbhdl_t* nbhandle);
		void NGA_NbPut(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], ga_nbhdl_t* nbhandle);
		Int  NGA_NbTest(ga_nbhdl_t* nbhandle);
		void NGA_NbWait(ga_nbhdl_t* nbhandle);
		void NGA_Put(Int g_a, Int lo[], Int hi[], void* buf, Int ld[]); 
		long NGA_Read_inc(Int g_a, Int ndim, Int subscript[], long inc);
		void NGA_Distribution(Int g_a, Int iproc, Int lo[], Int hi[]);
		void GA_Symmetrize(Int g_a);

	private:
		bool ga_initialized;
	
		typedef struct GA_t
		{
		    DistMatrix   < T, MC, MR >* DM;     // distmatrix instance
		    RmaInterface < T >*         rmaint; // rmainterface instance
		} GA;

		// vector of GA handles
		std::vector < GA > ga_handles;
};
#endif // EL_ENABLE_RMA_GLOBAL_ARRAYS
} // namespace El
#endif // ifndef EL_GLOBALARRAYS_HPP
