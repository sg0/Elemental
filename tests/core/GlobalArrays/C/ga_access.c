#include<stdio.h>
#include<stdlib.h>

#if defined(USE_ELEMENTAL)
#include "El.h"
ElGlobalArrays_d eldga;
#else
#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#endif

//#define SIZE 12996
#define SIZE 10

int main(int argc, char **argv)
{
    int rank, nprocs;
    int g_A, dims[2], ld, block[2], *map, lo[2], hi[2];
    double val=1e-11, *A=NULL;
#if defined(USE_ELEMENTAL)
    ElInitialize( &argc, &argv );
    ElMPICommRank( MPI_COMM_WORLD, &rank );
    ElMPICommSize( MPI_COMM_WORLD, &nprocs );
    // instantiate el::global array
    ElGlobalArraysConstruct_d( &eldga );
    // initialize global arrays
    ElGlobalArraysInitialize_d( eldga );
#else  
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MA_init(C_DBL, 1000, 1000);

    GA_Initialize();
#endif

    // block/map
    map = (int *)malloc(sizeof(int) * (1 + nprocs));
    if (NULL == map) {
	printf ("memory allocation failed\n");
	MPI_Abort (MPI_COMM_WORLD, 99);
    }

    block[0] = nprocs;
    block[1] = 1;

    for (int i = 0; i < nprocs; i++) {
	map[i] = i;
    }
    map[nprocs] = 0;

    dims[0] = nprocs;
    dims[1] = SIZE;

    // create irregular GA and initialize
#if defined(USE_ELEMENTAL)
    ElGlobalArraysCreateIrreg_d (eldga, 2, dims, "array_A", block, map, &g_A);
    ElGlobalArraysFill_d (eldga, g_A, &val);
#else
    g_A = NGA_Create_irreg(C_DBL, 2, dims, "array_A", block, map);
    GA_Fill(g_A, &val);
#endif
    // access local chunk and update
    lo[0] = rank;
    hi[0] = rank;
    lo[1] = 0;
    hi[1] = SIZE - 1;
#if defined(USE_ELEMENTAL)
    ElGlobalArraysAccess_d (eldga, g_A, lo, hi, &A, &ld);
#else
    NGA_Access (g_A, lo, hi, &A, &ld);
#endif
    if (rank == 0)
	printf ("ldim = %d and val = %e\n", ld, val);
    if (rank == 0)
	printf ("BEFORE: \n");
#if defined(USE_ELEMENTAL)
    ElGlobalArraysPrint_d (eldga, g_A);
#else
    GA_Print (g_A);
#endif
    
    int width = hi[0] - lo[0] + 1;
    int height = hi[1] - lo[1] + 1;

    for (int j = 0; j < width; j++)
	for (int i = 0; i < height; i++)
	    A[i + j * ld] *= val;

#if defined(USE_ELEMENTAL)
    ElGlobalArraysRelease_d (eldga, g_A, lo, hi);
#else
    NGA_Release(g_A, lo, hi);
#endif
    if (rank == 0)
	printf ("AFTER: \n");
#if defined(USE_ELEMENTAL)
    ElGlobalArraysPrint_d (eldga, g_A);
#else
    GA_Print (g_A);
#endif   
    // distribution information
    int my_lo[2] = {-1, -1}, my_hi[2] = {-1, -2};
    int my_ndim, my_dims[2];
#if defined(USE_ELEMENTAL)
    ElGlobalArraysInquire_d( eldga, g_A, &my_ndim, my_dims );
    ElGlobalArraysDistribution_d (eldga, g_A, rank, my_lo, my_hi);
#else
    int type;
    NGA_Inquire(g_A, &type, &my_ndim, my_dims);
    NGA_Distribution(g_A, rank, my_lo, my_hi);
#endif

    if (rank == 0)
    {
	printf ("ndim: %d\n", my_ndim);
	printf ("dims: %d %d\n", my_dims[0], my_dims[1]);
    }
    printf ("rank #[%d] range --> lo: (%d, %d) and hi: (%d, %d)\n", 
	    rank, my_lo[0], my_lo[1], my_hi[0], my_hi[1]);
    
#if defined(USE_ELEMENTAL)
    ElGlobalArraysDestroy_d (eldga, g_A);
#else
    GA_Destroy(g_A);    
#endif
    
    free (map);

#if defined(USE_ELEMENTAL)
    ElGlobalArraysTerminate_d( eldga );
    ElGlobalArraysDestruct_d( eldga );
    ElFinalize();
#else
    GA_Terminate();
    MPI_Finalize();
#endif
}
