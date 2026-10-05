#ifndef TINKER_MACRO_H
#define TINKER_MACRO_H
#include "tinker_precision.h"
#include "tinker_types.h"

#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
#  define TINKERHP_LIBM_SIGN copysignf
#  define TINKERHP_LIBM_ABS  fabsf
#  define TINKERHP_LIBM_SQRT sqrtf
#  define TINKERHP_LIBM_FLOOR floorf
#  define TINKERHP_LIBM_ERFC erfcf
#  define TINKERHP_LIBM_EXP  expf
#  define TINKERHP_LIBM_MIN  fminf
#  define TINKERHP_LIBM_MAX  fmaxf
#  define TINKERHP_DEVICE_INV __frcp_rn
#else
#  define TINKERHP_LIBM_SIGN copysign
#  define TINKERHP_LIBM_ABS  fabs
#  define TINKERHP_LIBM_SQRT sqrt
#  define TINKERHP_LIBM_FLOOR floor
#  define TINKERHP_LIBM_ERFC erfc
#  define TINKERHP_LIBM_EXP  exp
#  define TINKERHP_LIBM_MIN  fmin
#  define TINKERHP_LIBM_MAX  fmax
#  define TINKERHP_DEVICE_INV __drcp_rn
#endif
#endif
