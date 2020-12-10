#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#if defined(SINGLE) || defined(MIXED)
   typedef float real;
#  define f_floor floorf
#  ifdef USE_ERFC_HASTINGS
#    define f_erfc(x) erfcf_hastings(x)*exp2a
#  else
#    define f_erfc(x) erfcf(x)
#  endif
#  define f_sqrt  sqrtf
#  define f_inv(x) 1./x
#  define f_exp expf
#  define ti_eps 1e-5
#  define prec_eps 2.384185791015625e-7
#  define r_Format " %f"
#else
   typedef double real;
#  define f_floor floor
#  define f_erfc  erfc
#  define f_sqrt  sqrt
#  define f_inv(x) 1./x
#  define ti_eps 1e-12
#  define prec_eps 4.440892098500626e-16
#  define r_Format " %lf"
#  define r10_Format " %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf "
#  define i10_Format " %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i "
#  define f_exp exp
#endif

#ifdef MIXED
   typedef double realm;
#elif defined(SINGLE)
   typedef float realm;
#else
   typedef double realm;
#endif


#define EXTERN_C_BEG extern "C" {
#define EXTERN_C_END }

#define restrict __restrict__
#define WARP_SIZE  32
#define BLOCK_DIM 128

#define One 1
#define Onef 1.0
#define Two 2
#define Twof 2.0
#define Half 0.5

const int RED_BUFF_SIZE (1<<13);
const unsigned int ALL_LANES( 0xffffffff );

typedef struct Real3 { real x,y,z; }real3;
typedef struct Real3m { realm x,y,z; }real3m;

typedef struct Real6 { real x,y,z,xx,yy,zz; }real6;
typedef struct Rpole_elt { real c,dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz; }rpole_elt;


// This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
// the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
// error of 1.5e-7.
__device__
inline static float erfcf_hastings(float x)
{
   //float exp2a = expf(-x * x);
   float t = 1.0f / (1.0f + 0.3275911f * x);
   return (0.254829592f + (-0.284496736f + (1.421413741f + (-1.453152027f + 1.061405429f * t) * t) * t) * t) * t;//*exp2a;
}

/* ==
   Check Kernel launch result
   ==
*/
inline static void gpuAssert(cudaError_t code, const char* file, int line, int abort=1, const cudaError_t success=cudaSuccess){
   if (code!=success){
      fprintf(stderr,"GPUAssert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
inline static void gpuAssert(cusolverStatus_t code, const char* file, int line, int abort=1){
   if (code!=CUSOLVER_STATUS_SUCCESS){
      fprintf(stderr,"GPUAssert: cuSolver Error : %d %s %d\n", code, file, line);
      if (abort) exit(code);
   }
}
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__, cudaSuccess); }
#define gpuErrchkSolver(ans) { gpuAssert((ans), __FILE__, __LINE__); }


#define cudaKernelMaxGridSize(gS,kernel,bS,dynSMem)                                               \
   gpuErrchk(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&gS,kernel,bS,dynSMem))                \
   gS *= devProp.multiProcessorCount;



#endif
