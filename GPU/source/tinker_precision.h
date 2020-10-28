#ifndef TINKER_PRECISION_H
#define TINKER_PRECISION_H
#  ifdef TINKER_CUDART_H
#    warning( "tinker_cudart.h  should not be included before tinker_precision.h")
#  endif
#  ifdef SINGLE
#    define t_p 4
#    define r_p 4
#    define MPI_TPREC MPI_REAL4
#    define MPI_RPREC MPI_REAL4
#    define MPI_SHARED_TYPE 4_MPI_ADDRESS_KIND
#    define MPI_SHMRED_TYPE 4_MPI_ADDRESS_KIND
#    define M_PPTRS SPPTRS
#    define M_PPTRF SPPTRF
#    define M_dgesv sgesv
#  else
#    ifdef MIXED
#      define t_p 4
#      define r_p 8
#      define MPI_TPREC MPI_REAL4
#      define MPI_RPREC MPI_REAL8
#      define MPI_SHARED_TYPE 4_MPI_ADDRESS_KIND
#      define MPI_SHMRED_TYPE 8_MPI_ADDRESS_KIND
#      define M_PPTRS SPPTRS
#      define M_PPTRF SPPTRF
#      define M_dgesv sgesv
#    else
#      define t_p 8
#      define r_p 8
#      define MPI_TPREC MPI_REAL8
#      define MPI_RPREC MPI_REAL8
#      define MPI_SHARED_TYPE 8_MPI_ADDRESS_KIND
#      define MPI_SHMRED_TYPE 8_MPI_ADDRESS_KIND
#      define M_PPTRS DPPTRS
#      define M_PPTRF DPPTRF
#      define M_dgesv dgesv
#    endif
#  endif
#endif

#if defined(_CUDA) && defined(USE_NVSHMEM)
#  ifndef USE_NVSHMEM_CUDA
#    define USE_NVSHMEM_CUDA
#  endif
#endif
