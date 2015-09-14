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
		GlobalArrays(DistMatrix< T > & DM);
		GlobalArrays(DistMatrix< T > & DM, Int height, Int width);
		~GlobalArrays();

		typedef Int ga_nbhdl_t;

		// Subset of GA C API
		// TODO follow either NGA_ or GA_ in function names
		int  GA_Create_handle();
		void GA_Set_data (int g_a, int ndim, int dims[], int type);
		int  GA_Allocate(int g_a);
		int  GA_Create(int type, int ndim, int dims[], const char *array_name);
		void GA_Copy(int g_a, int g_b); 
		void GA_Print_distribution(int g_a);
		void GA_Destroy(int g_a);
		void GA_Add(void *alpha, int g_a, void* beta, int g_b, int g_c); 
		void GA_Dgemm(char ta, char tb, int m, int n, int k, double alpha, int g_a, int g_b, double beta, int g_c );
		int  GA_Duplicate(int g_a, const char *array_name);
		void GA_Fill(int g_a, void *value);
		void GA_Initialize();
		void GA_Sync();
		void GA_Terminate();
		void GA_Transpose(int g_a, int g_b);
		void NGA_Access(int g_a, int lo[], int hi[], void *ptr, int ld[]);
		void NGA_Acc(int g_a, int lo[], int hi[],void* buf,int ld[],void* alpha);
		void NGA_Get(int g_a, int lo[], int hi[], void* buf, int ld[]); 
		void NGA_NbAcc(int g_a,int lo[], int hi[],void* buf,int ld[],void* alpha, ga_nbhdl_t* nbhandle);
		void NGA_NbGet(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle);
		void NGA_NbPut(int g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle);
		int  NGA_NbTest(ga_nbhdl_t* nbhandle);
		void NGA_NbWait(ga_nbhdl_t* nbhandle);
		void NGA_Put(int g_a, int lo[], int hi[], void* buf, int ld[]); 
		long NGA_Read_inc(int g_a, int subscript[], long inc);
		void NGA_Distribution(int g_a, int iproc, int lo[], int hi[]);
		void GA_Symmetrize(int g_a);

	private:
		bool ga_initialized;
		bool ga_dm_dim_initialized;
	
		// ga creation status
		typedef enum ga_status_
		{
		    UNDEFINED, // undefined, while init
		    CREATED, // handle created
		    SET, // called ga_set_data
		    ALLOCATED // memory allocated
		} ga_status_t;

		struct GA
		{
			int handle; // integer handle
			int ndims; // number of dimensions
			int dims[2]; // x and y dims of distmatrix
			bool pending_transfer; // whether there is a pending xfer to/from this ga
			ga_status_t status; // whether GA is set, allocated or just handle created
			DistMatrix < T > DM; // distmatrix  
			RmaInterface < T > rmaint; // rmainterface
		};

		// vector of GA handles
		std::vector < struct GA > ga_handles;
};
#endif // EL_ENABLE_RMA_GLOBAL_ARRAYS
} // namespace El
#endif // ifndef EL_GLOBALARRAYS_HPP
