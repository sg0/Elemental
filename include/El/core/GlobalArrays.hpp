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
		Int  GA_Create(Int ndim, Int dims[], const char *array_name, Int chunk[]);
		Int  GA_Create_irreg(Int ndim, Int dims[], const char *array_name, Int block[], Int map[]);
		Int  GA_Duplicate(Int g_a, const char *array_name);
		/* Create rmainterface object, and attach DM corresponding to g_a */
		Int  GA_Allocate(Int g_a);
		void GA_Print(Int g_a);
		/* call Detach on rmainterface, entry is not erased from ga_handles until terminate */
		void GA_Destroy(Int g_a);	
		/* Check if Detach was called already on rmainterface objects, and erase ga_handles */
		void GA_Terminate();	

		void GA_Fill(Int g_a, T *value);
		void GA_Copy(Int g_a, Int g_b); 
		void GA_Sync();
		void GA_Gop(T x[], Int n, char op);
		T GA_Dot( Int g_a, Int g_b );

		void NGA_Distribution(Int g_a, Int iproc, Int lo[], Int hi[]);
		void NGA_Inquire(Int g_a, Int * ndim, Int dims[]);
		void NGA_Access(Int g_a, Int lo[], Int hi[], T** ptr, Int ld[]);
		void NGA_Release(Int g_a, Int lo[], Int hi[]);
		void NGA_Release_update(Int g_a, Int lo[], Int hi[]);

		void NGA_Acc(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], T* alpha);
		void NGA_Get(Int g_a, Int lo[], Int hi[], T* buf, Int ld[]); 
		void NGA_NbAcc(Int g_a,Int lo[], Int hi[], T* buf, Int ld[], T* alpha, ga_nbhdl_t* nbhandle);
		void NGA_NbGet(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], ga_nbhdl_t* nbhandle);
		void NGA_NbPut(Int g_a, Int lo[], Int hi[], T* buf, Int ld[], ga_nbhdl_t* nbhandle);
		void NGA_NbWait(ga_nbhdl_t* nbhandle);
		void NGA_Put(Int g_a, Int lo[], Int hi[], T* buf, Int ld[]); 
		
		T NGA_Read_inc(Int g_a, Int subscript[], T inc);

		void GA_Symmetrize(Int g_a);
		void GA_Transpose(Int g_a, Int g_b);
		void GA_Add(T *alpha, Int g_a, T* beta, Int g_b, Int g_c); 
		void GA_Dgemm(char ta, char tb, Int m, Int n, Int k, T alpha, Int g_a, Int g_b, T beta, Int g_c );

	private:
		bool ga_initialized;
	
		// find the actual GA handle
		Int FindGAIndex(Int g_a);

		struct GA
		{
		    // ga handle
		    Int ga_handle_; 
		    // distmatrix_mc_mr instance
		    DistMatrix   < T, MC, MR >* DM;     
		    // rmainterface instance
		    RmaInterface < T >* rmaint; 
		    // local height/width of distmatrix
		    std::vector< Int > ga_local_height;
		    std::vector< Int > ga_local_width;
		    // lo/hi coordinates - 2*p
		    std::vector< Int > ga_lo;
		    std::vector< Int > ga_hi;
		    // is there a pending rma op on
		    // this GA
		    bool rma_local_pending;
		    // is the GA destroyed already?
		    bool is_destroyed;
		    // this is only required when a DistMatrix
		    // is not associated with GA, at present
		    // this is the case with fetch-and-op
		    Int length;
		    // number of dimensions
		    Int ndim;
		    // GA patch width
		    Int patchWidth;
		    // GA patch height
		    Int patchHeight;
		    // MPI_Window for fetch-and-op/read-inc
		    // for 1D GA
		    mpi::Window fop_win;
		    // initialize
		    GA() :
			ga_handle_( -1 ),                // GA handle
			DM( nullptr ),                   // DM
		        rmaint( nullptr ),               // RMA interface
	                ga_local_height(),               // local height
		        ga_local_width(),                // local width
			ga_lo(),                      	 // lo
			ga_hi(),                         // hi
		        rma_local_pending( false ),      // rma_local_pending
			is_destroyed( true ),            // is_destroyed
		        length( -1 ),                    // length
		        ndim( -1 ),                      // ndim
			patchWidth( -1 ),                // patchWidth
			patchHeight( -1 ),               // patchHeight
		        fop_win( mpi::WIN_NULL )         // fop_win
		    {}
		};
				
		// vector of GA handles
		std::vector < GA > ga_handles;

		// Local matrix for holding nga_acess related
		// data...also used for locally nonblocking
		// accumulate
		struct matrix_params_
		{
		    Int ga_index_;
		    Matrix< T > * M_;
		    Int lo[2];
		    Int hi[2];
		    bool is_accumulate_;
		    // initialize
		    matrix_params_() :
			ga_index_( -1 ), 		// ga handle the matrix is associated with
			M_( nullptr ),   		// matrix to be allocated
			lo{-1, -1},  			// lo index associated with matrix
			hi{-2, -2},  			// hi index associated with matrix
			is_accumulate_( false )		// is this matrix used in locally nonblocking accumulate
			{}
		};

		// vector of matrix parameters
		std::vector < matrix_params_ > matrices_;

		// local completion handles
		struct nbhdl_t_
		{
		    bool active_;
		    Int nbhandle_;
		    Matrix< T > * M_;
		    // initialize
		    nbhdl_t_() :
			active_( true ),
			nbhandle_( -1 ),
			M_( nullptr )
		    {}
		};
		
		std::vector < nbhdl_t_ > nbhdls_;
};
#endif // EL_ENABLE_RMA_GLOBAL_ARRAYS
} // namespace El
#endif // ifndef EL_GLOBALARRAYS_HPP
