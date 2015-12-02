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
ElGlobalArrays_i eliga;
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
  int g_A, g_B, dims[DIM]={SIZE,SIZE}, val1=5, val2=4;
  int lo[DIM]={SIZE-SIZE,SIZE-SIZE}, hi[DIM]={SIZE-1,SIZE-1}, ld=SIZE;

#if defined(USE_ELEMENTAL)
  // initialize Elemental (which will initialize MPI)
  ElInitialize( &argc, &argv );
  ElMPICommRank( MPI_COMM_WORLD, &rank );
  ElMPICommSize( MPI_COMM_WORLD, &nprocs );
  // instantiate el::global array
  ElGlobalArraysConstruct_i( &eliga );
  // initialize global arrays
  ElGlobalArraysInitialize_i( eliga );
#else
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
#endif

  // create global arrays
#if defined(USE_ELEMENTAL)
  ElGlobalArraysCreate_i( eliga, DIM, dims, "array_A", NULL, &g_A );
#else
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysDuplicate_i( eliga, g_A, "array_B", &g_B );
#else  
  g_B = GA_Duplicate(g_A, "array_B");
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysFill_i( eliga, g_A, &val1 );
  ElGlobalArraysFill_i( eliga, g_B, &val2 );
#else
  GA_Fill(g_A, &val1);
  GA_Fill(g_B, &val2);
#endif

  int dot_AB = -99;
#if defined(USE_ELEMENTAL)
  ElGlobalArraysDot_i( eliga, g_A, g_B, &dot_AB );
#else
  dot_AB = GA_Idot(g_A, g_B);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_i( eliga );
#else
  GA_Sync();
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysPrint_i( eliga, g_A );
  ElGlobalArraysPrint_i( eliga, g_B );
#else 
  GA_Print(g_A);
  GA_Print(g_B);
#endif
 
  // Check
#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_i( eliga );
#else
  GA_Sync();
#endif

  if(rank==0)
      printf ("Integer dot product of g_A and g_B: %d\n", dot_AB);

#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_i( eliga );
#else
  GA_Sync();
#endif

  if(rank==0)
    printf("Test Completed \n");

#if defined(USE_ELEMENTAL)
  ElGlobalArraysDestroy_i( eliga, g_A );
  ElGlobalArraysDestroy_i( eliga, g_B );
#else
  GA_Destroy(g_A);
  GA_Destroy(g_B);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysTerminate_i( eliga );
  // call el::global arrays destructor
  ElGlobalArraysDestruct_i( eliga );
  ElFinalize();
#else
  GA_Terminate();
  MPI_Finalize();
#endif

  return 0;
}
