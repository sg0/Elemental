/*
 * Test Program for GA
 * This is to test GA_NbAcc/NbPut
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#if defined(USE_ELEMENTAL)
#include "El.h"
ElGlobalArrays_i eliga;
#else
#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#endif

#define DIM	2
#define SIZE 	5
#define NITERS 	10

#if defined(USE_ELEMENTAL)
void GA_Error (const char *str, int item)
{
    printf ("%s: %d\n", str, item);
    MPI_Abort( MPI_COMM_WORLD, item );
}
#endif

int main(int argc, char **argv)
{
  int rank, nprocs;
  int g_A;
  int *local_A=NULL, *local_B=NULL, *output_A=NULL;
  int dims[DIM]={SIZE,SIZE}, dims2[DIM], lo[DIM]={SIZE-SIZE,SIZE-SIZE}, hi[DIM]={SIZE-1,SIZE-1}, ld=SIZE;
  int value=SIZE;

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

  local_A=(int*)malloc(SIZE*SIZE*sizeof(int));
  output_A=(int*)malloc(SIZE*SIZE*sizeof(int));
  memset (output_A, 0, SIZE*SIZE*sizeof(int));
  for(int j=0; j<SIZE; j++)
      for(int i=0; i<SIZE; i++) local_A[i+j*ld]=(i + j);
      //for(int i=0; i<SIZE; i++) local_A[i+j*ld]=(rand()%10);

  local_B=(int*)malloc(SIZE*SIZE*sizeof(int));
  memset (local_B, 0, SIZE*SIZE*sizeof(int));

  // nb handle
#if defined(USE_ELEMENTAL)
  typedef ElInt ga_nbhdl_t;
#endif
  ga_nbhdl_t nbnb;

#if defined(USE_ELEMENTAL)
  ElGlobalArraysCreate_i( eliga, DIM, dims, "array_A", NULL, &g_A );
  ElGlobalArraysFill_i( eliga, g_A, &value );
#else
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  GA_Fill(g_A, &value);
#endif

  if (rank == 0) printf ("Initial global array:\n");
#if defined(USE_ELEMENTAL)
  ElGlobalArraysPrint_i( eliga, g_A );
#else
  GA_Print(g_A);
#endif

  for (int i = 0; i < NITERS; i++)
  {
      // acc data
#if defined(USE_ELEMENTAL)
      ElGlobalArraysNBAccumulate_i( eliga, g_A, lo, hi, local_A, &ld, &value, &nbnb );
#else
      NGA_NbAcc(g_A, lo, hi, local_A, &ld, &value, &nbnb);
#endif
      
      // updated output
      MPI_Reduce (local_A, output_A, SIZE*SIZE, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

#if defined(USE_ELEMENTAL)
      ElGlobalArraysNBWait_i( eliga, &nbnb );
#else
      NGA_NbWait (&nbnb);
#endif 

      // get
      if (rank == 0) printf ("Get in iter #%d\n", i);
#if defined(USE_ELEMENTAL)
      ElGlobalArraysSync_i( eliga );
      ElGlobalArraysGet_i( eliga, g_A, lo, hi, local_B, &ld );
      ElGlobalArraysPrint_i( eliga, g_A );
#else
      GA_Sync();
      NGA_Get(g_A, lo, hi, local_B, &ld);
      GA_Print(g_A);
#endif
  } // end of iters

  if(rank==0)
    {
      printf(" Alpha (multiplier): %d\n", value);
      printf(" Original local buffer (before accumulation): \n");

      for(int i=0; i<SIZE; i++)
	{
	  for(int j=0; j<SIZE; j++)
	    printf("%d ", local_A[i*ld+j]);
	  printf("\n");
	}
      printf("\n");
      printf(" Get returns: \n");
      for(int i=0; i<SIZE; i++)
	{
	  for(int j=0; j<SIZE; j++)
	    printf("%d ", local_B[i*ld + j]);
	  printf("\n");
	}

      printf("\n");
      for(int i=0; i<SIZE; i++)
	{
	  for(int j=0; j<SIZE; j++)
	    {
	      if(local_B[i*ld+j]!=(value + (NITERS * value * (output_A[i*ld+j]))))
		  GA_Error("ERROR", -99);
	    }
	}
    }
#if defined(USE_ELEMENTAL)
  ElGlobalArraysDestroy_i( eliga, g_A );
#else
  GA_Destroy(g_A);
#endif
  if(rank == 0)
    printf ("OK. Test passed\n");

    free (local_A);
    free (local_B);
    free (output_A);

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
