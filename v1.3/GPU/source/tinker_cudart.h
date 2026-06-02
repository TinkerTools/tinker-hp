#ifndef TINKER_CUDART_H
#define TINKER_CUDART_H

#ifndef TINKER_CUF
#   define f_abs abs
#   define f_sign sign
#   define f_floor floor
#   define f_sqrt sqrt
#   define M_subroutine subroutine
#   define M_function  function
#   define f_inv(x) 1/x
#   define f_erfc(x) erfc(x)
#   define f_exp(x) exp(x)
#   define f_min(x) min(x)
#   define f_max(x) max(x)
#   define WRITE_C(x)
#else
#   define M_subroutine attributes(device) subroutine
#   define M_function attributes(device) function
#   define WRITE_C(x) x
#   if !defined(SINGLE) && !defined(MIXED)
#      define f_sqrt(x) sqrt(x)
#      define f_exp(x) exp(x)
#   else
#      if defined(USE_ERFC_HASTINGS)
#         define f_erfc(x) erfcf_hastings(x)*exp2a
#      endif
#   endif
#endif

#endif
