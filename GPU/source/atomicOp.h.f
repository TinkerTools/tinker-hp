#ifndef ATOMICOP_H_F
#define ATOMICOP_H_F

c      interface atomic_add1
c      attributes(device) subroutine atomic_add_d( dat,val,rstat )
c      real(8)        :: dat
c      real(8)        :: val,rstat
c!DIR$ ignore_tkr (rd) dat
c!!pgi$ ignore_tkr (rd) dat
c      end subroutine
c      end interface

      interface atomic_add
        module procedure atomic_add_d
        module procedure atomic_add_s
        module procedure atomic_add_m
#ifdef USE_DETERMINISTIC_REDUCTION
        module procedure atomic_add_f
        module procedure atomic_add_f1
        module procedure atomic_add_f2
#endif
      end interface

      interface atomic_sub
        module procedure atomic_sub_d
        module procedure atomic_sub_s
        module procedure atomic_sub_m
#ifdef USE_DETERMINISTIC_REDUCTION
        module procedure atomic_sub_f
        module procedure atomic_sub_f1
        module procedure atomic_add_f2
#endif
      end interface

#endif
