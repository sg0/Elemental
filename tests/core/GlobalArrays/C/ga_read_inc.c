/*
 * Test Program for GA
 * This is to test NGA_Read_inc
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
  int rank, nprocs;
  int g_A, dims[DIM]={SIZE,SIZE}; 
  int lo[DIM]={SIZE-SIZE,SIZE-SIZE}, hi[DIM]={SIZE-1,SIZE-1}, ld=SIZE;
  int val=5;

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

  MA_init(C_DBL, 1000, 1000);

  GA_Initialize();
#endif

  // create global arrays
#if defined(USE_ELEMENTAL)
  ElGlobalArraysCreate_i( eliga, DIM, dims, "array_A", NULL, &g_A );
#else
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysFill_i( eliga, g_A, &val );
#else
  GA_Fill(g_A, &val);
#endif

  if (rank == 0) printf ("Initially: \n");
#if defined(USE_ELEMENTAL)
  ElGlobalArraysPrint_i( eliga, g_A );
#else 
  GA_Print(g_A);
#endif

  int subscript[2];
  subscript[0] = 2;
  subscript[1] = 2;
  long inc = 4;
  long prev = -99;
#if defined(USE_ELEMENTAL)
  ElGlobalArraysReadIncrement_i( eliga, g_A, subscript, inc, (ElInt *)&prev );
#else
  prev = NGA_Read_inc(g_A, subscript, inc);
#endif

  if (rank == 0) printf ("After read increment: \n");
#if defined(USE_ELEMENTAL)
  ElGlobalArraysPrint_i( eliga, g_A );
#else 
  GA_Print(g_A);
#endif

  printf ("[%d] Previous value at position {%d, %d}: %d\n", 
	  rank, subscript[0], subscript[1], prev);
 
  // Check
#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_i( eliga );
#else
  GA_Sync();
#endif

  if(rank==0)
    printf("Test Completed \n");

#if defined(USE_ELEMENTAL)
  ElGlobalArraysDestroy_i( eliga, g_A );
#else
  GA_Destroy(g_A);
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
