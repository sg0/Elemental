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
  int g_A, g_B, g_C, *local_C=NULL, dims[DIM]={SIZE,SIZE}, val1=5, val2=4, alpha=3, beta=2;
  int clo[DIM]={SIZE-SIZE,SIZE-SIZE}, chi[DIM]={SIZE-1,SIZE-1}, ld=SIZE;

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

  local_C=(int*)malloc(SIZE*SIZE*sizeof(int));

  // create global arrays
#if defined(USE_ELEMENTAL)
  ElGlobalArraysCreate_i( eliga, DIM, dims, "array_A", NULL, &g_A );
#else
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysDuplicate_i( eliga, g_A, "array_B", &g_B );
  ElGlobalArraysDuplicate_i( eliga, g_A, "array_C", &g_C );
#else  
  g_B = GA_Duplicate(g_A, "array_B");
  g_C = GA_Duplicate(g_A, "array_C");
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysFill_i( eliga, g_A, &val1 );
  ElGlobalArraysFill_i( eliga, g_B, &val2 );
#else
  GA_Fill(g_A, &val1);
  GA_Fill(g_B, &val2);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysAdd_i( eliga, &alpha, g_A, &beta, g_B, g_C );
#else
  GA_Add(&alpha, g_A, &beta, g_B, g_C);
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_i( eliga );
#else
  GA_Sync();
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysPrint_i( eliga, g_A );
  ElGlobalArraysPrint_i( eliga, g_B );
  ElGlobalArraysPrint_i( eliga, g_C );
#else 
  GA_Print(g_A);
  GA_Print(g_B);
  GA_Print(g_C);
#endif
 
  // Check
#if defined(USE_ELEMENTAL)
  ElGlobalArraysSync_i( eliga );
#else
  GA_Sync();
#endif

  if(rank==0)
  { 
#if defined(USE_ELEMENTAL)
      ElGlobalArraysGet_i( eliga, g_C, clo, chi, local_C, &ld );
#else
      NGA_Get(g_C, clo, chi, local_C, &ld);
#endif
      for(i=0; i<SIZE; i++)
	  for(j=0; j<SIZE; j++)
	      if(local_C[i+j*ld]!=(alpha*val1)+(beta*val2)) 
		  GA_Error ("ERROR: Values not matching", -99);
  }

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
  ElGlobalArraysDestroy_i( eliga, g_C );
#else
  GA_Destroy(g_A);
  GA_Destroy(g_B);
  GA_Destroy(g_C);
#endif

  free (local_C);

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
