#ifndef TINKER_PRECISION_H
#define TINKER_PRECISION_H
#  ifdef TINKER_CUDART_H
#    warning( "tinker_cudart.h  should not be included before tinker_precision.h")
#  endif
#  define TINKER_SINGLE_PREC 0
#  define TINKER_MIXED_PREC  0
#  define TINKER_DOUBLE_PREC 0
#  ifdef SINGLE
#    undef  TINKER_SINGLE_PREC
#    define TINKER_SINGLE_PREC 1
#    define t_p 4
#    define r_p 4
#    define MPI_TPREC MPI_REAL4
#    define MPI_RPREC MPI_REAL4
#    define MPI_SHARED_TYPE 4_MPI_ADDRESS_KIND
#    define MPI_SHMRED_TYPE 4_MPI_ADDRESS_KIND
#    define M_PPTRS SPPTRS
#    define M_PPTRF SPPTRF
#    define M_gesv  sgesv
#    define M_gesvm sgesv
#  else
#    ifdef MIXED
#      undef  TINKER_MIXED_PREC
#      define TINKER_MIXED_PREC 1
#      define t_p 4
#      define r_p 8
#      define MPI_TPREC MPI_REAL4
#      define MPI_RPREC MPI_REAL8
#      define MPI_SHARED_TYPE 4_MPI_ADDRESS_KIND
#      define MPI_SHMRED_TYPE 8_MPI_ADDRESS_KIND
#      define M_PPTRS SPPTRS
#      define M_PPTRF SPPTRF
#      define M_gesv  sgesv
#      define M_gesvm dgesv
#    else
#      undef  TINKER_DOUBLE_PREC
#      define TINKER_DOUBLE_PREC 1
#      define t_p 8
#      define r_p 8
#      define MPI_TPREC MPI_REAL8
#      define MPI_RPREC MPI_REAL8
#      define MPI_SHARED_TYPE 8_MPI_ADDRESS_KIND
#      define MPI_SHMRED_TYPE 8_MPI_ADDRESS_KIND
#      define M_PPTRS DPPTRS
#      define M_PPTRF DPPTRF
#      define M_gesv  dgesv
#      define M_gesvm dgesv
#    endif
#  endif
#endif

#if ( TINKER_SINGLE_PREC + TINKER_MIXED_PREC +   \
      TINKER_DOUBLE_PREC ) != 1
#   error find error in PRECISION macro !!!
#endif

#  if TINKER_MIXED_PREC
#  else
#    ifdef USE_DETERMINISTIC_REDUCTION
#      warning 'Disabling Fixed point Arithmetic ! only to be \
#               Used with Mixed Precision Build -DMIXED'
#      undef USE_DETERMINISTIC_REDUCTION
#    endif
#  endif

#if defined(_CUDA) && defined(USE_NVSHMEM)
#  ifndef USE_NVSHMEM_CUDA
#    define USE_NVSHMEM_CUDA
#  endif
#endif
