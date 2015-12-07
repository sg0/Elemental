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

#define NITERS 10

int main(int argc, char **argv)
{
    int rank, nprocs;
    int g_A, ld, lo[2], hi[2];
    double val=2., *A;
    /*
    int dims[2] = { 114, 114 };
    int block[2];
    block[0] = nprow;
    block[1] = npcol;
    int map[4] = { 0, 56, 0, 56 };
    */
    int dims[2] = { 10, 10 };
    int block[2] = { 2, 2 };
    int map[4] = { 0, 5, 0, 5 };
    const int nprow = block[0];
    const int npcol = block[1];
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

    assert (nprocs == (nprow * npcol));

    // create irregular GA and initialize
#if defined(USE_ELEMENTAL)
    ElGlobalArraysCreateIrreg_d (eldga, 2, dims, "array_A", block, map, &g_A);
    ElGlobalArraysFill_d (eldga, g_A, &val);
#else
    g_A = NGA_Create_irreg(C_DBL, 2, dims, "array_A", block, map);
    GA_Fill(g_A, &val);
#endif

    // start iterations
    for (int i = 0; i < NITERS; i++)
    {
	// access local chunk and update
#if defined(USE_ELEMENTAL)
	ElGlobalArraysDistribution_d( eldga, g_A, rank, lo, hi );
	ElGlobalArraysAccess_d (eldga, g_A, lo, hi, &A, &ld);
#else    
	NGA_Distribution(g_A, rank, lo, hi);
	NGA_Access (g_A, lo, hi, &A, &ld);
#endif
	// distribution information
	if (i == 0)
	    printf ("rank #[%d] range --> lo: (%d, %d) and hi: (%d, %d)\n", 
		    rank, lo[0], lo[1], hi[0], hi[1]);

	if (rank == 0)
	    printf ("ldim = %d\n", ld);
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
	ElGlobalArraysReleaseUpdate_d (eldga, g_A, lo, hi);
#else
	NGA_Release_update (g_A, lo, hi);
#endif

#if defined(USE_ELEMENTAL)
	ElGlobalArraysSync_d (eldga);
#else
	GA_Sync ();
#endif   
	if (rank == 0)
	    printf ("Local portion of GA after calling release\n");

	for (int k = 0; k < nprocs; k++)
	{
	    if (k == rank)
	    {
		printf ("rank #%d\n", k);
		for (int j = 0; j < width; j++)
		{
		    for (int i = 0; i < height; i++)
		    {
			printf ("%f ", A[i + j * ld]);
		    }
		    printf ("\n");
		}
	    }
	    MPI_Barrier (MPI_COMM_WORLD);
	}

	if (rank == 0)
	    printf ("AFTER: \n");

#if defined(USE_ELEMENTAL)
	ElGlobalArraysPrint_d (eldga, g_A);
#else
	GA_Print (g_A);
#endif
    } // end of NITERS

#if defined(USE_ELEMENTAL)
    ElGlobalArraysDestroy_d (eldga, g_A);
#else
    GA_Destroy(g_A);    
#endif
    
#if defined(USE_ELEMENTAL)
    ElGlobalArraysTerminate_d( eldga );
    ElGlobalArraysDestruct_d( eldga );
    ElFinalize();
#else
    GA_Terminate();
    MPI_Finalize();
#endif
}
