#ifndef ATOMICOP_H_F
#warning("__FILE__ header undetected !! common interface unavalaible") 
#endif
#ifndef ATOMICOP_INC
#define ATOMICOP_INC
#include "convert.f.inc"

#ifdef TINKER_CUF
      !
      !  CUDA Section
      !
      attributes(device) subroutine atomic_add_i( dat,val )
      implicit none
      integer,device:: dat
      integer,intent(in)::val
      integer rstat
      rstat = atomicAdd(dat, val)
      end subroutine
      attributes(device) subroutine atomic_add_s( dat,val,rstat )
      implicit none
      real(4),device :: dat
      real(4),intent(in)::val
      real(4) :: rstat
      rstat = atomicAdd(dat, val)
      end subroutine
      attributes(device) subroutine atomic_add_d( dat,val,rstat )
      implicit none
      real(8),device :: dat
      real(8),intent(in)::val
      real(8) :: rstat
      rstat = atomicAdd(dat, val)
      end subroutine
      attributes(device) subroutine atomic_add_c( dat,val )
      implicit none
      real(t_p),device:: dat
      real(t_p),intent(in)::val
      real(t_p) stat
      stat = atomicAdd( dat,val )
      end subroutine
      attributes(device) subroutine atomic_add_m( dat,val )
      implicit none
      real(r_p),device:: dat
      real(t_p),intent(in)::val
      real(r_p) stat
      stat = atomicAdd(dat, real(val,r_p))
      end subroutine
      attributes(device) subroutine atomic_add_m1( dat,val )
      implicit none
      real(r_p),device:: dat
      real(r_p),intent(in)::val
      real(r_p) stat
      stat = atomicAdd(dat, val)
      end subroutine
      attributes(device) subroutine atomic_add_f( dat,val )
      implicit none
      mdyn_rtyp,device:: dat
      mdyn_rtyp,intent(in)::val
      mdyn_rtyp stat
      stat = atomicAdd(dat, val)
      end subroutine
      attributes(device) subroutine atomic_add_f1( dat,val )
      implicit none
      mdyn_rtyp,device:: dat
      real(t_p),intent(in)::val
      mdyn_rtyp stat
      stat = atomicAdd(dat, tp2mdr(val))
      end subroutine
      attributes(device) subroutine atomic_add_f2( dat,val )
      implicit none
      mdyn_rtyp,device:: dat
      real(r_p),intent(in)::val
      mdyn_rtyp stat
      stat = atomicAdd(dat, rp2mdr(val))
      end subroutine

      attributes(device) subroutine atomic_sub_i( dat,val )
      implicit none
      integer,device:: dat
      integer,intent(in)::val
      integer rstat
      rstat = atomicSub(dat, val)
      end subroutine
      attributes(device) subroutine atomic_sub_d( dat,val )
      implicit none
      real(4),device:: dat
      real(4),intent(in)::val
      real(4) rstat
      rstat = atomicSub(dat, val)
      end subroutine
      attributes(device) subroutine atomic_sub_s( dat,val )
      implicit none
      real(8),device:: dat
      real(8),intent(in)::val
      real(8) rstat
      rstat = atomicSub(dat, val)
      end subroutine
      attributes(device) subroutine atomic_sub_m( dat,val )
      implicit none
      real(r_p),device:: dat
      real(t_p),intent(in)::val
      real(r_p) rstat
      rstat = atomicSub(dat, real(val,r_p))
      end subroutine
      attributes(device) subroutine atomic_sub_m1( dat,val )
      implicit none
      real(r_p),device:: dat
      real(r_p),intent(in)::val
      real(r_p) stat
      stat = atomicSub(dat, val)
      end subroutine
      attributes(device) subroutine atomic_sub_f( dat,val )
      implicit none
      mdyn_rtyp,device:: dat
      mdyn_rtyp,intent(in)::val
      mdyn_rtyp stat
      stat = atomicSub(dat, val)
      end subroutine
      attributes(device) subroutine atomic_sub_f1( dat,val )
      implicit none
      mdyn_rtyp,device:: dat
      real(t_p),intent(in)::val
      mdyn_rtyp stat
      stat = atomicSub(dat, tp2mdr(val))
      end subroutine
      attributes(device) subroutine atomic_sub_f2( dat,val )
      implicit none
      mdyn_rtyp,device:: dat
      real(r_p),intent(in)::val
      mdyn_rtyp stat
      stat = atomicSub(dat, rp2mdr(val))
      end subroutine
#else
      !
      !   OpenACC Section
      !
      subroutine atomic_add_s( dat,val )
!$acc routine
      implicit none
      real(4) dat
      real(4),intent(in)::val
!$acc atomic
      dat = dat + val
      end subroutine
      subroutine atomic_add_d( dat,val )
!$acc routine
      implicit none
      real(8) dat
      real(8),intent(in)::val
!$acc atomic
      dat = dat + val
      end subroutine
      subroutine atomic_add_m( dat,val )
!$acc routine
      implicit none
      real(8) dat
      real(4),intent(in)::val
      real(8) val1
      val1 = val
!$acc atomic
      dat = dat + val1
      end subroutine
      subroutine atomic_add_f( dat,val )
!$acc routine
      implicit none
      mdyn_rtyp dat
      mdyn_rtyp,intent(in)::val
!$acc atomic
      dat = dat + val
      end subroutine
      subroutine atomic_add_f1( dat,val )
!$acc routine
      implicit none
      mdyn_rtyp dat
      real(t_p),intent(in)::val
      mdyn_rtyp val1
      val1 = tp2mdr(val)
!$acc atomic
      dat = dat + val1
      end subroutine
      subroutine atomic_add_f2( dat,val )
!$acc routine
      implicit none
      mdyn_rtyp dat
      real(r_p),intent(in)::val
      mdyn_rtyp val1
      val1 = rp2mdr(val)
!$acc atomic
      dat = dat + val1
      end subroutine


      subroutine atomic_sub_d( dat,val )
!$acc routine
      implicit none
      real(4) dat
      real(4),intent(in)::val
!$acc atomic
      dat = dat - val
      end subroutine
      subroutine atomic_sub_s( dat,val )
!$acc routine
      implicit none
      real(8) dat
      real(8),intent(in)::val
!$acc atomic
      dat = dat - val
      end subroutine
      subroutine atomic_sub_m( dat,val )
!$acc routine
      implicit none
      real(8) dat
      real(4),intent(in)::val
      real(8) val1
      val1 = val
!$acc atomic
      dat = dat - val1
      end subroutine
      subroutine atomic_sub_f( dat,val )
!$acc routine
      implicit none
      mdyn_rtyp dat
      mdyn_rtyp,intent(in)::val
!$acc atomic
      dat = dat - val
      end subroutine
      subroutine atomic_sub_f1( dat,val )
!$acc routine
      implicit none
      mdyn_rtyp dat
      real(t_p),intent(in)::val
      mdyn_rtyp val1
      val1 = tp2mdr(val)
!$acc atomic
      dat = dat - val1
      end subroutine
      subroutine atomic_sub_f2( dat,val )
!$acc routine
      implicit none
      mdyn_rtyp dat
      real(r_p),intent(in)::val
      mdyn_rtyp val1
      val1 = rp2mdr(val)
!$acc atomic
      dat = dat - val1
      end subroutine
#endif

#endif
