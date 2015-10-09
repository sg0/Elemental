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
		void GA_Initialize();
		// TODO follow either NGA_ or GA_ in function names
		Int  GA_Create(Int ndim, Int dims[], const char *array_name);
		Int  GA_Duplicate(Int g_a, const char *array_name);
		/* Create rmainterface object, and attach DM corresponding to g_a */
		Int  GA_Allocate(Int g_a);
		void GA_Print(Int g_a);
		/* call Detach on rmainterface, entry is not erased from ga_handles until terminate */
		void GA_Destroy(Int g_a);	
		/* Check if Detach was called already on rmainterface objects, and erase ga_handles */
		void GA_Terminate();	

		void GA_Fill(Int g_a, void *value);
		void GA_Copy(Int g_a, Int g_b); 
		void GA_Sync();
		void GA_Gop(T x[], Int n, char op);
		T GA_Dot( Int g_a, Int g_b );

		void NGA_Distribution(Int g_a, Int iproc, Int lo[], Int hi[]);
		void NGA_Inquire(Int g_a, Int * ndim, Int dims[]);
		void NGA_Access(Int g_a, Int lo[], Int hi[], void **ptr, Int ld[]);
		
		void NGA_Acc(Int g_a, Int lo[], Int hi[],void* buf,Int ld[],void* alpha);
		void NGA_Get(Int g_a, Int lo[], Int hi[], void* buf, Int ld[]); 
		void NGA_NbAcc(Int g_a,Int lo[], Int hi[],void* buf,Int ld[],void* alpha, ga_nbhdl_t* nbhandle);
		void NGA_NbGet(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], ga_nbhdl_t* nbhandle);
		void NGA_NbPut(Int g_a, Int lo[], Int hi[], void* buf, Int ld[], ga_nbhdl_t* nbhandle);
		Int  NGA_NbTest(ga_nbhdl_t* nbhandle);
		void NGA_NbWait(ga_nbhdl_t* nbhandle);
		void NGA_Put(Int g_a, Int lo[], Int hi[], void* buf, Int ld[]); 
		long NGA_Read_inc(Int g_a, Int ndim, Int subscript[], long inc);

		void GA_Symmetrize(Int g_a);
		void GA_Transpose(Int g_a, Int g_b);
		void GA_Add(void *alpha, Int g_a, void* beta, Int g_b, Int g_c); 
		void GA_Dgemm(char ta, char tb, Int m, Int n, Int k, double alpha, Int g_a, Int g_b, double beta, Int g_c );

	private:
		bool ga_initialized;
	
		typedef struct GA_t
		{
		    // distmatrix_mc_mr instance
		    DistMatrix   < T >* DM;     
		    // rmainterface instance
		    RmaInterface < T >* rmaint; 
		    // local height/width of distmatrix
		    std::vector< Int >* ga_local_height;
		    std::vector< Int >* ga_local_width;
		    bool pending_rma_op;
		} GA;
		
		// vector of GA handles
		std::vector < GA > ga_handles;
};
#endif // EL_ENABLE_RMA_GLOBAL_ARRAYS
} // namespace El
#endif // ifndef EL_GLOBALARRAYS_HPP
