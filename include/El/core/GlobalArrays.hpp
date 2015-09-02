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

		// ga creation status
		typedef enum ga_status_
		{
		    CREATED, // handle created
		    SET, // called ga_set_data
		    ALLOCATED // memory allocated
		} ga_status_t;

		// Subset of GA C API
		int GA_Create_handle();
		void GA_Set_data (int g_a, int ndim, int dims[], int type);
		int  GA_Allocate(int g_a);
		void GA_Copy(int g_a, int g_b); 
		void GA_Destroy(int g_a);
		void GA_Add(void *alpha, int g_a, void* beta, int g_b, int g_c); 
		void GA_Dgemm(char ta, char tb, int m, int n, int k, double alpha, int g_a, int g_b, double beta, int g_c );
		int  GA_Duplicate(int g_a, char* array_name);
		void GA_Fill(int g_a, void *value);
		void GA_Initialize();
		void GA_Sync();
		void GA_Terminate();
		void GA_Transpose(int g_a, int g_b);
		void NGA_Access(int g_a, int lo[], int hi[], void *ptr, int ld[]);
		void NGA_Acc(int g_a, int lo[], int hi[],void* buf,int ld[],void* alpha);
		int  NGA_Allocate(int g_a);
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

		struct GA
		{
		    int handle; // integer handle
		    mpi::Datatype dtype; // data type of ga
		    int ndims; // number of dimensions
		    int dims[2]; // x and y dims of distmatrix
		    bool pending_transfer; // whether there is a pending xfer to/from this ga
		    ga_status_t status; // whether GA is set, allocated or just handle created
		    mpi::Comm comm; // comm for window allocation and everything
		    DistMatrix < T > DM; // distmatrix      
		    RmaInterface < T > rmaint; // rma object
		};
		// handle vector of GAs
		std::vector < struct GA > ga_handles;
	};
#endif // EL_ENABLE_RMA_GLOBAL_ARRAYS
} // namespace El
#endif // ifndef EL_GLOBALARRAYS_HPP
