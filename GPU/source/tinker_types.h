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
#endif
