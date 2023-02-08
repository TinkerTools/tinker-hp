#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>


// C++11
#ifdef __cplusplus
#   if __cplusplus < 201103L
#      error Must enable C++11.
#   endif
#endif

#if defined(__INTEL_COMPILER)
#   define TINKER_ICPC

#elif defined(__PGIC__)
#   define TINKER_PGI

#elif defined(__clang__)
#   define TINKER_CLANG
// xcode clang is different from llvm clang
#   ifdef __apple_build_version__
#      define TINKER_APPLE_CLANG
#   else
#      define TINKER_LLVM_CLANG
#   endif

#elif defined(__GNUC__)
#   define TINKER_GCC
#   if __GNUC__ <= 4 && __GNUC_MINOR__ <= 8
#      warning Your default GNU version is 4.8 where C++11 is incomplete.
#   endif

#endif

// Suppress Warnings
#ifdef TINKER_ICPC
// #161: unrecognized #pragma
#   pragma warning disable 161
#endif
#ifdef TINKER_CLANG
#   pragma clang diagnostic ignored "-Wextern-c-compat"
#endif

/// \ingroup cpp_syntax
/// `if constexpr` has been added to C++ since C++17.
/// `if CONSTEXPR` expands to `if constexpr` if this feature is supported.
/// Otherwise it expands to `if`.
#if __cplusplus >= 201703L && defined(__cpp_if_constexpr)
#   define CONSTEXPR constexpr
#else
#   define CONSTEXPR
#endif

/// \ingroup cpp_syntax
/// Reduces the "unused variable" warnings from the compiler.
#ifdef __has_cpp_attribute
#   if __has_cpp_attribute(maybe_unused)
#      define MAYBE_UNUSED [[maybe_unused]]
#   else
#      define MAYBE_UNUSED
#   endif
#elif defined(TINKER_ICPC) || defined(TINKER_CLANG) || defined(TINKER_GCC)
#   define MAYBE_UNUSED __attribute__((unused))
#else
#   define MAYBE_UNUSED
#endif

#define EXTERN_C_BEG extern "C" {
#define EXTERN_C_END }

#if defined(SINGLE) || defined(MIXED)
   typedef float real;
#  define f_floor floorf
#  ifdef USE_ERFC_HASTINGS
#    define f_erfc(x) erfcf_hastings(x)*exp2a
#  else
#    define f_erfc(x) erfcf(x)
#  endif
#  define f_sqrt  sqrtf
#  define f_pow   powf
#  define f_inv(x) (1/static_cast<float>(x))
#  define f_exp expf
#  define f_cos cosf
#  define f_sin sinf
#  define f_min fminf
#  define f_max fmaxf
#  define ti_eps 1e-5
#  define prec_eps 2.384185791015625e-7
#  define r_Format " %f"
#else
   typedef double real;
#  define f_floor floor
#  define f_erfc  erfc
#  define f_sqrt  sqrt
#  define f_pow   pow
#  define f_inv(x) (1/static_cast<double>(x))
#  define f_exp exp
#  define f_cos cosf
#  define f_sin sin
#  define f_min fmin
#  define f_max fmax
#  define ti_eps 1e-12
#  define prec_eps 4.440892098500626e-16
#  define r_Format " %lf"
#  define r10_Format " %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf "
#  define i10_Format " %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i "
#endif

#ifdef MIXED
   typedef double realm;
#elif defined(SINGLE)
   typedef float realm;
#else
   typedef double realm;
#endif


/// \ingroup cpp_syntax
/// Expands to \c __restrict__, which is a common C++ extension.
#ifdef __cplusplus
#define restrict __restrict__
#endif

#define WARP_SIZE  32
#define BLOCK_DIM  4*WARP_SIZE

const int RED_BUFF_SIZE (1<<13);
const unsigned int ALL_LANES( 0xffffffff );

typedef struct Real3 { real x,y,z; }real3;
typedef struct Real3m { realm x,y,z; }real3m;

typedef struct Real6 { real x,y,z,xx,yy,zz; }real6;
typedef struct Rpole_elt { real c,dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz; }rpole_elt;

namespace CodePrm {
const int  One    = 1;
const int  Two    = 2;
const real Half   = 0.5;
const real sqrtpi =1.772453850905516027;
const real zeror  = 0.0;
const real oner   = 1.0;
const real twor   = 2.0;
const realm zerom = 0.0;
const realm Halfm = 0.5;
const realm onem  = 1.0;
const realm twom  = 2.0;
const int SCAL    = 1<<1;
const int EWALD   = 1<<8;
const int THOLE   = 1<<9;
const int DIRDAMP = 1<<10;
const int CHGPEN  = 1<<11;
}


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

__device__
inline static real f_mmax( const real& val0, const real& val1 ){
   real mmax = f_min(val0,val1);
   return mmax==CodePrm::zeror ? f_max(val0,val1) : mmax;
}
__device__
inline static real f_mmax( const real& val0, real&& val1 ){
   real mmax = f_min(val0,val1);
   return mmax==CodePrm::zeror ? f_max(val0,val1) : mmax;
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
