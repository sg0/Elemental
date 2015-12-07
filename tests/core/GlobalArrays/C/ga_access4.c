#include<stdio.h>
#include<stdlib.h>
#include <assert.h>

#if defined(USE_ELEMENTAL)
#include "El.h"
ElGlobalArrays_d eldga;
#else
#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#endif

#define NITERS 	10
#define SIZE	10

int main(int argc, char **argv)
{
    int rank, nprocs;
    int g_A, g_B, ld, lo[2], hi[2];
    double val=2., *A;
    int *map, dims[2], block[2];
    /*
    int dims[2] = { 114, 114 };
    int block[2];
    block[0] = nprow;
    block[1] = npcol;
    int map[4] = { 0, 56, 0, 56 };
    int dims[2] = { 10, 10 };
    int block[2] = { 2, 2 };
    int map[4] = { 0, 5, 0, 5 };
    */
#if defined(USE_ELEMENTAL)
    ElInitialize( &argc, &argv );
    ElMPICommRank( MPI_COMM_WORLD, &rank );
    ElMPICommSize( MPI_COMM_WORLD, &nprocs );
    // instantiate el::global array
    ElGlobalArraysConstruct_d( &eldga );
    // initialize global arrays
    ElGlobalArraysInitialize_d( eldga );
    typedef ElInt ga_nbhdl_t;
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
    const int nprow = block[0];
    const int npcol = block[1];

    assert (nprocs == (nprow * npcol));

    double value = 1.;
    ga_nbhdl_t nbnb;

    // create irregular GA and initialize
#if defined(USE_ELEMENTAL)
    ElGlobalArraysCreateIrreg_d( eldga, 2, dims, "array_A", block, map, &g_A );
    ElGlobalArraysDuplicate_d( eldga, g_A, "array_B", &g_B );
    ElGlobalArraysFill_d (eldga, g_A, &val);
#else
    g_A = NGA_Create_irreg(C_DBL, 2, dims, "array_A", block, map);
    g_B = GA_Duplicate (g_A, "array_B");
    GA_Fill (g_A, &val);
#endif

    // start iterations
    for (int i = 0; i < NITERS; i++)
    {
	if (rank == 0)
	    printf ("[%d] GA A: \n", i);

#if defined(USE_ELEMENTAL)
	ElGlobalArraysPrint_d (eldga, g_A);
#else
	GA_Print (g_A);
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

	// acc data
#if defined(USE_ELEMENTAL)
	ElGlobalArraysNBAccumulate_d( eldga, g_B, lo, hi, A, &ld, &value, &nbnb );
#else
	NGA_NbAcc (g_B, lo, hi, A, &ld, &value, &nbnb);
#endif

#if defined(USE_ELEMENTAL)
	ElGlobalArraysRelease_d (eldga, g_A, lo, hi);
#else
	NGA_Release (g_A, lo, hi);
#endif

#if defined(USE_ELEMENTAL)
	ElGlobalArraysNBWait_d( eldga, &nbnb );
#else
	NGA_NbWait (&nbnb);
#endif   

	if (rank == 0)
	    printf ("[%d] GA B: \n", i);

#if defined(USE_ELEMENTAL)
	ElGlobalArraysSync_d( eldga );
#else
	GA_Sync ();
#endif 

#if defined(USE_ELEMENTAL)
	ElGlobalArraysPrint_d (eldga, g_B);
#else
	GA_Print (g_B);
#endif
    } // end of NITERS

#if defined(USE_ELEMENTAL)
    ElGlobalArraysDestroy_d (eldga, g_A);
    ElGlobalArraysDestroy_d (eldga, g_B);
#else
    GA_Destroy(g_A);    
    GA_Destroy(g_B);    
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
