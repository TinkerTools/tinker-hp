#ifndef TINKER_TYPES_H
#define TINKER_TYPES_H
#  ifndef TINKER_PRECISION_H
#    error include tinker_precision.h before tinker_types.h
#  endif

#  define md_p r_p
#  define en_p r_p
#  define Md_p(x) x ## _re_p
#  define En_p(x) x ## _re_p
#  define mdyn_rtyp real(md_p)
#  define ener_rtyp real(en_p)
#  define MPI_MDTYP MPI_RPREC

#  ifdef USE_DETERMINISTIC_REDUCTION
#    undef mdyn_rtyp
#    undef ener_rtyp
#    undef MPI_MDTYP
#    define mdyn_rtyp integer(8)
#    define ener_rtyp integer(8)
#    define MPI_MDTYP MPI_INTEGER8
#  endif

#define __concat(x,y) x ## y
#define __concat2(x,y,z) x ## y ## z
#define __char(x) #x

#define __Cat(x,y) __concat(x,y)
#define __Cat2(x,y,z) __concat2(x,y,z)
#define __cuCat2(x,y,z) __concat2(x,y,z)
#define __Char(x) __char(x)

#define __TINKER_FATAL__ call fatal1(__FILE__,__LINE__)

#define __use_grd__   1
#define __use_ene__   2
#define __use_vir__   4
#define __use_act__   8
#define __use_sca__  16

#define __use_mpi__          1
#define __radepsOpt__        2
#define __use_softcore__     4
#define __use_groups__       8
#define __use_gamd__        64
#define __use_polymer__    128
#define __use_shortRange__ 256
#define __use_longRange__  512
#define __use_ewald__     1024
#define __use_lambdadyn__ 2048

#endif
