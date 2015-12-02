#ifndef HAVE_STDIO_H
#define HAVE_STDIO_H
#endif
#ifndef HAVE_MATH_H
#define HAVE_MATH_H
#endif

#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif

#ifdef HAVE_STDIO_H
#   include <stdio.h>
#endif
#ifdef HAVE_MATH_H
#   include <math.h>
#endif

#if defined(USE_ELEMENTAL)
#include "El.h"
ElGlobalArrays_d eldga;
#else
#include "ga.h"
#include "macdecls.h"
#include "mp3.h"
#endif

#define N 100            /* dimension of matrices */
int me, nproc;

#if defined(USE_ELEMENTAL)
void GA_Error (const char *str, int item)
{
    printf ("%s: %d\n", str, item);
    exit(-1);
}
#endif

void do_work()
{
int ONE=1 ;   /* useful constants */
int g_a, g_b;
int n=N;

#if defined(USE_ELEMENTAL)
#else
int type=MT_F_DBL;
int me=GA_Nodeid(), nproc=GA_Nnodes();
#endif
int i, row;
int dims[2]={N,N};
int lo[2], hi[2];

/* Note: on all current platforms DoublePrecision == double */
double buf[N], err, alpha, beta;

#if defined(USE_ELEMENTAL)
     ElGlobalArraysCreate_d( eldga, 2, dims, "A", NULL, &g_a );
     if(me==0)printf("Creating matrix A\n");
#else
     if(me==0)printf("Creating matrix A\n");
     g_a = NGA_Create(type, 2, dims, "A", NULL);
     if(!g_a) GA_Error("create failed: A",n); 
#endif
     if(me==0)printf("OK\n");
     if(me==0)printf("Creating matrix B\n");
     /* create matrix B  so that it has dims and distribution of A*/

#if defined(USE_ELEMENTAL)
     ElGlobalArraysDuplicate_d( eldga, g_a, "B", &g_b );
     double zero = 0.0;
     ElGlobalArraysFill_d( eldga, g_a, &zero );
#else
     g_b = GA_Duplicate(g_a, "B");
     if(! g_b) GA_Error("duplicate failed",n); 
     if(me==0)printf("OK\n");

     GA_Zero(g_a);   /* zero the matrix */
#endif     

     if(me==0)printf("Initializing matrix A\n");
     /* fill in matrix A with random values in range 0.. 1 */ 
     lo[1]=0; hi[1]=n-1;
     for(row=me; row<n; row+= nproc){
         /* each process works on a different row in MIMD style */
         lo[0]=hi[0]=row;   
         for(i=0; i<n; i++) buf[i]=sin((double)i + 0.1*(row+1));
#if defined(USE_ELEMENTAL)
         ElGlobalArraysPut_d( eldga, g_a, lo, hi, buf, &ONE );
#else    
         NGA_Put(g_a, lo, hi, buf, &ONE);
#endif
     }

#if defined(USE_ELEMENTAL)
     ElGlobalArraysSync_d( eldga );
#endif

     if(me==0)printf("Symmetrizing matrix A\n");
#if defined(USE_ELEMENTAL)
     ElGlobalArraysSymmetrize_d( eldga, g_a );
#else
     GA_Symmetrize(g_a);   /* symmetrize the matrix A = 0.5*(A+A') */
#endif

     /* check if A is symmetric */ 
     if(me==0)printf("Checking if matrix A is symmetric\n");
#if defined(USE_ELEMENTAL)
     ElGlobalArraysTranspose_d( eldga, g_a, g_b );
#else
     GA_Transpose(g_a, g_b); /* B=A' */
#endif
     alpha=1.; beta=-1.;

#if defined(USE_ELEMENTAL)
     ElGlobalArraysAdd_d( eldga, &alpha, g_a, &beta, g_b, g_b );
     ElGlobalArraysDot_d( eldga, g_b, g_b, &err );
#else
     GA_Add(&alpha, g_a, &beta, g_b, g_b);  /* B= A - B */
     err= GA_Ddot(g_b, g_b);
#endif 
     if(me==0)printf("Error=%f\n",(double)err);
     
     if(me==0)printf("\nChecking atomic accumulate \n");

#if defined(USE_ELEMENTAL)
     ElGlobalArraysFill_d( eldga, g_a, &zero );
#else
     GA_Zero(g_a);   /* zero the matrix */
#endif 
     for(i=0; i<n; i++) buf[i]=(double)i;

     /* everybody accumulates to the same location/row */
     alpha = 1.0;
     row = n/2;
     lo[0]=hi[0]=row;
     lo[1]=0; hi[1]=n-1;
#if defined(USE_ELEMENTAL)
     ElGlobalArraysAccumulate_d( eldga, g_a, lo, hi, buf, &ONE, &alpha );
     ElGlobalArraysSync_d( eldga );
#else    
     NGA_Acc(g_a, lo, hi, buf, &ONE, &alpha );
     GA_Sync();
#endif
     if(me==0){ /* node 0 is checking the result */

#if defined(USE_ELEMENTAL)
        ElGlobalArraysGet_d( eldga, g_a, lo, hi, buf, &ONE );
#else    
        NGA_Get(g_a, lo, hi, buf,&ONE);
#endif
        for(i=0; i<n; i++) if(buf[i] != (double)nproc*i)
           GA_Error("failed: column=",i);
        printf("PASSED\n\n");
     }
     
#if defined(USE_ELEMENTAL)
     ElGlobalArraysDestroy_d( eldga, g_a );
     ElGlobalArraysDestroy_d( eldga, g_b );
#else
     GA_Destroy(g_a);
     GA_Destroy(g_b);
#endif
}
     
int main(argc, argv)
int argc;
char **argv;
{
#if defined(USE_ELEMENTAL)
  // initialize Elemental (which will initialize MPI)
  ElInitialize( &argc, &argv );
  ElMPICommRank( MPI_COMM_WORLD, &me );
  ElMPICommSize( MPI_COMM_WORLD, &nproc );
  // instantiate el::global array
  ElGlobalArraysConstruct_d( &eldga );
  // initialize global arrays
  ElGlobalArraysInitialize_d( eldga );
#else
    int heap=20000, stack=20000;
    MP_INIT(argc,argv);

    GA_INIT(argc,argv);                            /* initialize GA */
    me=GA_Nodeid(); 
    nproc=GA_Nnodes();
    if(me==0) {
       if(GA_Uses_fapi())GA_Error("Program runs with C array API only",1);
       printf("Using %ld processes\n",(long)nproc);
       fflush(stdout);
    }
#endif

#if defined(USE_ELEMENTAL)
#else
    heap /= nproc;
    stack /= nproc;
    if(! MA_init(MT_F_DBL, stack, heap)) 
       GA_Error("MA_init failed",stack+heap);  /* initialize memory allocator*/ 
#endif
    do_work();

#if defined(USE_ELEMENTAL)
    ElGlobalArraysSync_d( eldga );
#endif
    if(me==0)printf("Terminating ..\n");

#if defined(USE_ELEMENTAL)
    ElGlobalArraysTerminate_d( eldga );
    // call el::global arrays destructor
    ElGlobalArraysDestruct_d( eldga );
    ElFinalize();
#else
    GA_Terminate();
    MP_FINALIZE();
#endif
    return 0;
}
