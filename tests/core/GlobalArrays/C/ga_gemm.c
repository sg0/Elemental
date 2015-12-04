/*
 * Test Program for GA
 * This is to test GA_Add (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --used to duplicate and generate one more global array.., handle 'g_A' to 'g_B'
 * 
 *
 */
#include <stdio.h>
#include <stdlib.h>

#if defined(USE_ELEMENTAL)
#include "El.h"
ElGlobalArrays_d eldga;
#else
#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#endif

#define DIM 2
#define SIZE 5

#if defined(USE_ELEMENTAL)
void GA_Error (const char *str, int item)
{
    printf ("%s: %d\n", str, item);
    MPI_Abort (MPI_COMM_WORLD, item);
}
#endif

int main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, g_C, *local_C=NULL, dims[DIM]={SIZE,SIZE};
  double val1=5., val2=4., alpha=3., beta=2.;

#if defined(USE_ELEMENTAL)
  // initialize Elemental (which will initialize MPI)
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

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
#endif

  // create global arrays
#if defined(USE_ELEMENTAL)
  ElGlobalArraysCreate_d( eldga, DIM, dims, "array_A", NULL, &g_A );
#else
  g_A = NGA_Create(C_DBL, DIM, dims, "array_A", NULL);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysDuplicate_d( eldga, g_A, "array_B", &g_B );
  ElGlobalArraysDuplicate_d( eldga, g_A, "array_C", &g_C );
#else  
  g_B = GA_Duplicate(g_A, "array_B");
  g_C = GA_Duplicate(g_A, "array_C");
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysFill_d( eldga, g_A, &val1 );
  ElGlobalArraysFill_d( eldga, g_B, &val2 );
#else
  GA_Fill(g_A, &val1);
  GA_Fill(g_B, &val2);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysDgemm_d( eldga, 'N', 'T', SIZE, SIZE, SIZE, alpha, g_A, g_B, beta, g_C );
#else
  GA_Dgemm('N', 'T', SIZE, SIZE, SIZE, alpha, g_A, g_B, beta, g_C);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_d( eldga );
#else
  GA_Sync();
#endif

  if (rank == 0)
      printf ("alpha = %f and beta = %f\n", alpha, beta);
#if defined(USE_ELEMENTAL)
  ElGlobalArraysPrint_d( eldga, g_A );
  ElGlobalArraysPrint_d( eldga, g_B );
  ElGlobalArraysPrint_d( eldga, g_C );
#else 
  GA_Print(g_A);
  GA_Print(g_B);
  GA_Print(g_C);
#endif
 
  // Check
#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_d( eldga );
#else
  GA_Sync();
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysDestroy_d( eldga, g_A );
  ElGlobalArraysDestroy_d( eldga, g_B );
  ElGlobalArraysDestroy_d( eldga, g_C );
#else
  GA_Destroy(g_A);
  GA_Destroy(g_B);
  GA_Destroy(g_C);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysTerminate_d( eldga );
  // call el::global arrays destructor
  ElGlobalArraysDestruct_d( eldga );
  ElFinalize();
#else
  GA_Terminate();
  MPI_Finalize();
#endif

  return 0;
}
