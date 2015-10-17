/*
 * Test Program for GA
 * This is to test GA_Duplicate (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --helps to duplicate the content from handle 'g_A' 
 * Here used GA_Inquire to verify that g_A hanle returns the right values of created_array
 */

#include<stdio.h>

#if defined(USE_ELEMENTAL)
#include "El.h"
ElGlobalArrays_d eldga;
#define MAXDUP 100
#else
#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#endif

#define DIM 2

int main(int argc, char **argv)
{
  int rank, nprocs, i;
  int g_A, g_B;
#if defined(USE_ELEMENTAL)
  int g_C;
#endif
  int dims[DIM]={5,5}, dims2[DIM], ndim2, type2, dims3[DIM], ndim3, type3;
  double value=5, val2=4;

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

  MA_init(C_DBL, 1000, 1000);

  GA_Initialize();
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysCreate_d( eldga, DIM, dims, "array_A", &g_A );
  ElGlobalArraysFill_d( eldga, g_A, &value );
   
  ElGlobalArraysDuplicate_d( eldga, g_A, "array_B", &g_B );
  ElGlobalArraysPrint_d( eldga, g_A );
  ElGlobalArraysSync_d( eldga );
   
  ElGlobalArraysFill_d( eldga, g_B, &val2 );
  ElGlobalArraysPrint_d( eldga, g_B );
  
  // new array C
  int dims4[2];
  dims4[0] = 1;
  dims4[1] = dims[1];
  ElGlobalArraysCreate_d( eldga, DIM, dims4, "array_C", &g_C );
  ElGlobalArraysFill_d( eldga, g_C, &value );
  // create multiple duplicates
  int g[MAXDUP];
  for (int i = 0; i < MAXDUP; i++)
      ElGlobalArraysDuplicate_d( eldga, g_C,  "C Dup arrays", &g[i] );
#else
  g_A = NGA_Create(C_DBL, DIM, dims, "array_A", NULL);
  GA_Fill(g_A, &value);

  g_B = GA_Duplicate(g_A, "array_B");
  GA_Print(g_A);
  GA_Sync();
  GA_Fill(g_B, &val2);
  GA_Print(g_B);
#endif
  if(rank==0)
    {
#if defined(USE_ELEMENTAL)
      ElGlobalArraysInquire_d( eldga, g_A, &ndim2, dims2 );
      ElGlobalArraysInquire_d( eldga, g_A, &ndim3, dims3 );
#else
      NGA_Inquire(g_A, &type2, &ndim2, dims2);
      NGA_Inquire(g_A, &type3, &ndim3, dims3);
#endif
#if defined(USE_ELEMENTAL)
      printf(" %d \n", ndim2);
      if(ndim2!=ndim3 )
	printf("ERROR : Not equal \n");

      if(ndim2==ndim3 )
	printf("Equal \n");

      printf ("duplicate array handles: \n");
      for (int i = 0; i < MAXDUP; i++)
	  printf (" %d ", g[i]);
      printf ("\n\n");
#else
      printf(" %d -- %d,,\n", type2, ndim2);
      if(type2!=type3 || ndim2!=ndim3 )
	printf("ERROE : \n");

      if(type2==type3 || ndim2==ndim3 )
	printf("ERROE : Equal \n");
#endif
      for(i=0; i<DIM; i++)
	{
	  if(dims2[i]!=dims3[i])
	    printf("ERROE : \n");
	  printf("%d: %d[ %d] ...* \n", rank, i, dims2[i]);
	}
    }

  if(rank == 1)
    printf("Test Completed \n");

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
