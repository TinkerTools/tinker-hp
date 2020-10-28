#include "tinker_cudart.h"

      M_subroutine 
     &            elj1_couple(rik2,xr,yr,zr,rv,eps,cut2
     &                  ,cut,off,e,ded)
!$acc routine
        use tinheader ,only: ti_p,re_p
        use tintypes  ,only: real3
#if defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
        use utilcu    ,only: f_sqrt
#endif
        implicit none
        real(t_p)  ,intent(in ):: rik2,xr,yr,zr,rv,eps,cut2
     &             ,cut,off
        real(t_p)  ,intent(out):: e
        type(real3),intent(out):: ded
        real(t_p) p6,rik,taper,dtaper,de
        real(t_p) r,r2,r3,rinv

        rik = f_sqrt(rik2)
        r   = rv*rv/rik2
        p6  = r*r*r
        e   = eps * (p6-2.0_ti_p)*p6
        de  = eps * (p6-1.0_ti_p)*p6 * (-12.0_ti_p/rik)
c
c     use energy switching if near the cutoff distance
c
        if (rik2 .gt. cut2) then
           rinv   = 1.0_ti_p/(cut - off)
           r      = (rik - off) * rinv
           r2     = r * r
           r3     = r2 * r
           taper  = r3 * (6*r2 - 15*r + 10)
           dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;

           de     = e*dtaper + de*taper
           e      = e * taper
        end if
c
c     find the chain rule terms for derivative components
c
        de    = de / rik
        ded%x = de * xr
        ded%y = de * yr
        ded%z = de * zr
      end

      M_subroutine 
     &            elj3_couple(rik2,xr,yr,zr,rv,eps,cut2
     &                  ,cut,off,e)
!$acc routine
        use tinheader ,only: ti_p,re_p
        use tintypes  ,only: real3
#if defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
        use utilcu    ,only: f_sqrt
#endif
        implicit none
        real(t_p),intent(in ):: rik2,xr,yr,zr,rv
     &           ,eps,cut2,cut,off
        real(t_p),intent(out):: e
        real(t_p) p6,rik,taper,dtaper
        real(t_p) r,r2,r3,rinv

        r   = rv*rv/rik2
        p6  = r*r*r
        e   = eps * (p6-2.0_ti_p)*p6
c
c     use energy switching if near the cutoff distance
c
        if (rik2 .gt. cut2) then
           rinv   = 1.0_ti_p/(cut - off)
           r      = (f_sqrt(rik2) - off) * rinv
           r2     = r * r
           r3     = r2 * r
           taper  = r3 * (6*r2 - 15*r + 10)
           dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;

           e      = e * taper
        end if
      end

      M_subroutine
     &            eljshortlong1_couple(rik2,xr,yr,zr,rv,eps,cut2
     &                  ,off,shortheal,c0,c1,c2,c3,c4,c5,e,ded,mode)
!$acc routine
        use interfaces,only: long_mode
        use tinheader ,only: ti_p,re_p
        use tintypes  ,only: real3
#if defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
        use utilcu    ,only: f_sqrt
#endif
        implicit none
        integer    ,intent(in ):: mode
        real(t_p)  ,intent(in ):: rik2,xr,yr,zr,rv,eps,cut2
     &             ,c0,c1,c2,c3,c4,c5,off,shortheal
        real(t_p)  ,intent(out):: e
        type(real3),intent(out):: ded
        real(t_p) p6,p12,rik,rik3,rik4,rik5,taper,dtaper,de,s,ds

        rik = f_sqrt(rik2)
        p6  = rv**6 / rik2**3
        p12 = p6 * p6
        e   = eps * (p12-2.0*p6)
        de  = eps * (p12-p6) * (-12.0/rik)
c
c     use energy switching if near the cutoff distance
c
        call switch_respa_inl(rik,off,shortheal,s,ds)
        if (mode.eq.long_mode) then
           e  = (1-s)*e
           de = -e*ds-s*de
c
c     use energy switching if near the cutoff distance
c
           if (rik2 .gt. cut2) then
              rik3  = rik2 * rik
              rik4  = rik2 * rik2
              rik5  = rik2 * rik3
              taper = c5*rik5 + c4*rik4 + c3*rik3
     &              + c2*rik2 + c1*rik  + c0
             dtaper = 5.0*c5*rik4 + 4.0*c4*rik3
     &              + 3.0*c3*rik2 + 2.0*c2*rik + c1
              de    = e*dtaper + de*taper
              e     = e * taper
           end if
        else
           e   = e * s
           de  = e*ds + de*s
        end if
c
c     find the chain rule terms for derivative components
c
        de    = de / rik
        ded%x = de * xr
        ded%y = de * yr
        ded%z = de * zr
      end
