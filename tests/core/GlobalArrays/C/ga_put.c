/*
 * Test Program for GA
 * This is to test GA_Put (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A' 
 */

#include<stdio.h>
#include<stdlib.h>

#if defined(USE_ELEMENTAL)
#include "El.h"
ElGlobalArrays_i eliga;
#else
#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"
#endif

#define DIM 2
#define SIZE 5

#if defined(USE_ELEMENTAL)
void GA_Error (const char *str, int item)
{
    printf ("%s: %d\n", str, item);
    exit(-1);
}
#endif

int main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A;
#if defined(USE_ELEMENTAL)
  // el::ga supports pointer to buffer
  int *local_A=NULL, *local_B=NULL;
#else
  int **local_A=NULL, **local_B=NULL;
#endif
  int dims[DIM]={SIZE,SIZE}, dims2[DIM], lo[DIM]={SIZE-SIZE,SIZE-SIZE}, hi[DIM]={SIZE-1,SIZE-1}, ld=5, value=5;

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

#if defined(USE_ELEMENTAL)
  local_A=(int*)malloc(SIZE*SIZE*sizeof(int));
  for(i=0; i<SIZE; i++)
      for(j=0; j<SIZE; j++) local_A[i*SIZE+j]=rand()%10;

  local_B=(int*)malloc(SIZE*SIZE*sizeof(int));
  for(i=0; i<SIZE; i++)
      for(j=0; j<SIZE; j++) local_B[i*SIZE+j]=0;
#else
  local_A=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE; i++)
    {
      local_A[i]=(int*)malloc(SIZE*sizeof(int));
      for(j=0; j<SIZE; j++) local_A[i][j]=rand()%10;
    }

  local_B=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE; i++)
    {
      local_B[i]=(int*)malloc(SIZE*sizeof(int));
      for(j=0; j<SIZE; j++) local_B[i][j]=rand()%10;
    }
#endif

#if defined(USE_ELEMENTAL)
  ElGlobalArraysCreate_i( eliga, DIM, dims, "array_A", &g_A );
  ElGlobalArraysPrint_i( eliga, g_A );
  // put data
  ElGlobalArraysPut_i( eliga, g_A, lo, hi, local_A, &ld );
  // sync is required because El::GA::Put enforces local
  // completion, not remote
  ElGlobalArraysSync_i( eliga );
  // get -- does not require a sync when get returns,
  // data transfer to local buffer is complete
  ElGlobalArraysGet_i( eliga, g_A, lo, hi, local_B, &ld );
  ElGlobalArraysPrint_i( eliga, g_A );
#else
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  //GA_Fill(g_A, &value);
  //  GA_Zero(g_A);
  GA_Print(g_A);

  NGA_Put(g_A, lo, hi, local_A, &ld);
  
  NGA_Get(g_A, lo, hi, local_B, &ld);

  GA_Sync();
  GA_Print(g_A);
#endif

  if(rank==0)
    {
#if defined(USE_ELEMENTAL)
// by default, elemental is column ordered storage
      printf(" local_B \n");
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    printf("%d ", local_B[j*SIZE+i]);
	  printf("\n");
	}

      printf("\n");
      printf(" local_A \n");
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    printf("%d ", local_A[j*SIZE+i]);
	  printf("\n");
	}
      
      printf("\n");
      for(j=0; j<SIZE; j++)
	{
	  for(i=0; i<SIZE; i++)
	    {
	      if(local_B[j*SIZE+i]!=local_A[j*SIZE+i])
		printf("ERROR : in passing values \n");
	      /* there is erroe in the above piece of code 
	       * have to find method to solve it
	       */
	    }
	}
#else
      printf(" local_B \n");
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    printf("%d ", local_B[i][j]);
	  printf("\n");
	}

      printf("\n");
      printf(" local_A \n");
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    printf("%d ", local_A[i][j]);
	  printf("\n");
	}
      
      printf("\n");
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    {
	      if(local_B[i][j]!=local_A[i][j])
		printf("ERROR : in passing values \n");
	      /* there is erroe in the above piece of code 
	       * have to find method to solve it
	       */
	    }
	}
#endif
    }

#if defined(USE_ELEMENTAL)
    free (local_A);
    free (local_B);
#else
  for(i=0; i<SIZE; i++)
    {
      free(local_A[i]);
      free(local_B[i]);
    }
    free (local_A);
    free (local_B);
#endif

    //  GA_Print(g_A);
  //..GA_Sync();
#if defined(USE_ELEMENTAL)
  ElGlobalArraysDestroy_i( eliga, g_A );
#else
  GA_Destroy(g_A);
#endif
  if(rank == 0)
#if defined(USE_ELEMENTAL)
    printf ("OK.\n");
#else
    GA_PRINT_MSG();
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
}
