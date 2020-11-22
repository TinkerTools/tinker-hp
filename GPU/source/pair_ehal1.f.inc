#include "tinker_cudart.h"

        M_subroutine
     &             ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                         ,cut2,cut,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,muti,mutk
     &                         ,e,dedx,dedy,dedz)
!$acc routine
#if  defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
           use utilcu    ,only: f_sqrt
#endif
           use tinheader ,only: ti_p,re_p
           implicit none
           real(ti_p),intent(in) ::xpos,ypos,zpos,rik2
           real(ti_p) rv2,eps2,vscale
           real(ti_p),intent(in) ::cut2,cut,off,ghal,dhal
           real(ti_p),intent(in) ::scexp,vlambda,scalpha
           logical   ,intent(in) ::muti,mutk
           real(ti_p),intent(out)::e
           real(ti_p),intent(out)::dedx,dedy,dedz

           real(ti_p) rik,rinv,r,r2,r3
           !real(ti_p) rik3,rik4,rik5,rik6,rik7,rv7,rv7orho
           real(ti_p) dtau,gtau,tau,tau7,de,taper,dtaper
           real(ti_p) rho,rho6,rho7
           real(ti_p) scal,s1,s2,t1,t2,dt1drho,dt2drho
c
c          compute the energy contribution for this interaction
c
#ifndef TINKER_NO_MUTATE
           if (muti .neqv. mutk) then
              rik  = f_sqrt(rik2)
              rho  = rik / rv2
              rho6 = rho**6
              rho7 = rho6 * rho
              eps2 = eps2 * vlambda**scexp *vscale
              scal = scalpha * (1.0_ti_p-vlambda)**2
              s1   = 1.0_ti_p / (scal+(rho+dhal)**7)
              s2   = 1.0_ti_p / (scal+rho7+ghal)
              t1   = (1.0_ti_p+dhal)**7 * s1
              t2   = (1.0_ti_p+ghal) * s2
              dt1drho = -7.0_ti_p*(rho+dhal)**6 * t1 * s1
              dt2drho = -7.0_ti_p*rho6 * t2 * s2
              e    = eps2 * t1 * (t2-2.0_ti_p)
              de   = eps2 * (dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv2
           else
#endif
              rik  = f_sqrt(rik2)
              rho  = rik / rv2
              rho6 = rho**6
              rho7 = rho6 * rho
              eps2 = eps2 * vscale
              s1   = 1.0_ti_p / (rho+dhal)**7
              s2   = 1.0_ti_p / (rho7+ghal)
              t1   = (1.0_ti_p+dhal)**7 * s1
              t2   = (1.0_ti_p+ghal) * s2
              dt1drho = -7.0_ti_p*(rho+dhal)**6 * t1 * s1
              dt2drho = -7.0_ti_p*rho6 * t2 * s2
              e    = eps2 * t1 * (t2-2.0_ti_p)
              de   = eps2 * (dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv2

              ! Inadequate for single precision compute (rik**6 overflow risk)
c             rik     =  f_sqrt(rik2)
c             rik6    =  rik2 * rik2 * rik2
c             rik7    =  rik2 * rik2 * rik2 * rik
c             dtau    =  (rik + dhal * rv2)**(-1) !tau / (dhal + 1.0)
c             rv7     =  rv2**7
c             rv7orho =  rv7 / (rik7 + ghal * rv7)
c             tau7    =  (dtau * (dhal + 1.0))**7
c             !tau7    =  tau ** 7.0
c             gtau    =  eps2 * tau7 * rik6
c    &                 * (ghal + 1.0) * rv7orho * rv7orho * vscale
c             e       =  eps2 * tau7 * rv7
c    &                 * ((ghal + 1.0) * rv7orho - 2.0) * vscale
c             de      = - 7.0 * (dtau*e + gtau)
#ifndef TINKER_NO_MUTATE
           end if
#endif
c
c          use energy switching if near the cutoff distance
c
           if(rik2 > cut2) then ! mask energy switch
              ! Single precision overflow due to rik6
c             taper  =  c5 * rik2*rik2*rik + c4 * rik2*rik2
c    &                + c3 * rik2*rik      + c2 * rik2
c    &                + c1 * rik           + c0
c             dtaper =  5.0 * c5 * rik2*rik2
c    &                + 4.0 * c4 * rik2*rik
c    &                + 3.0 * c3 * rik2
c    &                + 2.0 * c2 * rik
c    &                +       c1

              rinv   = 1.0_ti_p/(cut - off)
              r      = (rik - off) * rinv
              r2     = r * r
              r3     = r2 * r
              taper  = r3 * (6*r2 - 15*r + 10)
              dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;

              de     = ( e * dtaper + de * taper ) / rik
              e      =   e * taper
           else
              de = de / rik
           endif

           dedx  =  de * xpos
           dedy  =  de * ypos
           dedz  =  de * zpos

        end

        M_subroutine
     &          ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                         ,cut2,shortcut
     &                         ,scexp,vlambda,scalpha,muti,mutk
     &                         ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
!$acc routine
#if  defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
           use utilcu    ,only: f_sqrt
#endif
           use tinheader ,only: ti_p,re_p
           implicit none
           real(ti_p),intent(in) ::xpos,ypos,zpos,rik2
           real(ti_p) rv2,eps2,vscale
           real(ti_p),intent(in) ::cut2,ghal,dhal
     &                            ,shortcut,shortheal
           real(ti_p),intent(in) ::scexp,vlambda,scalpha
           logical   ,intent(in) ::muti,mutk
           real(ti_p),intent(out)::e
           real(ti_p),intent(out)::dedx,dedy,dedz

           real(ti_p) rik,rik3,rik4,rik5,rik6,rik7,rv7,rv7orho
           real(ti_p) dtau,gtau,tau,tau7,de,s,ds
           real(ti_p) rho,rho6,rho7
           real(ti_p) scal,s1,s2,t1,t2,dt1drho,dt2drho
c
c          compute the energy contribution for this interaction
c
#ifndef TINKER_NO_MUTATE
           if (muti .neqv. mutk) then
              rik  = f_sqrt(rik2)
              rho  = rik / rv2
              rho6 = rho**6
              rho7 = rho6 * rho
              eps2 = eps2 * vlambda**scexp *vscale
              scal = scalpha * (1.0_ti_p-vlambda)**2
              s1   = 1.0_ti_p / (scal+(rho+dhal)**7)
              s2   = 1.0_ti_p / (scal+rho7+ghal)
              t1   = (1.0_ti_p+dhal)**7 * s1
              t2   = (1.0_ti_p+ghal) * s2
              dt1drho = -7.0_ti_p*(rho+dhal)**6 * t1 * s1
              dt2drho = -7.0_ti_p*rho6 * t2 * s2
              e    = eps2 * t1 * (t2-2.0_ti_p)
              de   = eps2 * (dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv2
           else
#endif
              rik  = f_sqrt(rik2)
              rho  = rik / rv2
              rho6 = rho**6
              rho7 = rho6 * rho
              eps2 = eps2 * vscale
              s1   = 1.0_ti_p / (rho+dhal)**7
              s2   = 1.0_ti_p / (rho7+ghal)
              t1   = (1.0_ti_p+dhal)**7 * s1
              t2   = (1.0_ti_p+ghal) * s2
              dt1drho = -7.0_ti_p*(rho+dhal)**6 * t1 * s1
              dt2drho = -7.0_ti_p*rho6 * t2 * s2
              e    = eps2 * t1 * (t2-2.0_ti_p)
              de   = eps2 * (dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv2
#ifndef TINKER_NO_MUTATE
           end if
#endif
c
c          use energy switching if near the cutoff distance
c
           call switch_respa_inl(rik,shortcut,shortheal,s,ds)

           e     =   e* s
           de    = ( e*ds + de*s ) / rik

           dedx  =  de * xpos
           dedy  =  de * ypos
           dedz  =  de * zpos
        end

        M_subroutine
     &          ehal1_couple_long(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                         ,cut2,cut,off,shortcut
     &                         ,scexp,vlambda,scalpha,muti,mutk
     &                         ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
!$acc routine
#if  defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
           use utilcu    ,only: f_sqrt
#endif
           use tinheader ,only: ti_p,re_p
           implicit none
           real(ti_p),intent(in) ::xpos,ypos,zpos,rik2
           real(ti_p) rv2,eps2,vscale
           real(ti_p),intent(in) ::cut2,cut,off,ghal,dhal
     &                            ,shortcut,shortheal
           real(ti_p),intent(in) ::scexp,vlambda,scalpha
           logical   ,intent(in) ::muti,mutk
           real(ti_p),intent(out)::e
           real(ti_p),intent(out)::dedx,dedy,dedz

c          real(ti_p) rik,rik3,rik4,rik5,rik6,rik7,rv7,rv7orho
c          real(ti_p) dtau,gtau,tau,tau7
           real(ti_p) de,taper,dtaper,s,ds
           real(ti_p) rik,rinv,r,r2,r3
           real(ti_p) rho,rho6,rho7
           real(ti_p) scal,s1,s2,t1,t2,dt1drho,dt2drho
c
c          compute the energy contribution for this interaction
c
#ifndef TINKER_NO_MUTATE
           if (muti .neqv. mutk) then
              rik  = f_sqrt(rik2)
              rho  = rik / rv2
              rho6 = rho**6
              rho7 = rho6 * rho
              eps2 = eps2 * vlambda**scexp *vscale
              scal = scalpha * (1.0_ti_p-vlambda)**2
              s1   = 1.0_ti_p / (scal+(rho+dhal)**7)
              s2   = 1.0_ti_p / (scal+rho7+ghal)
              t1   = (1.0_ti_p+dhal)**7 * s1
              t2   = (1.0_ti_p+ghal) * s2
              dt1drho = -7.0_ti_p*(rho+dhal)**6 * t1 * s1
              dt2drho = -7.0_ti_p*rho6 * t2 * s2
              e    = eps2 * t1 * (t2-2.0_ti_p)
              de   = eps2 * (dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv2
           else
#endif
              rik  = f_sqrt(rik2)
              rho  = rik / rv2
              rho6 = rho**6
              rho7 = rho6 * rho
              eps2 = eps2 * vscale
              s1   = 1.0_ti_p / (rho+dhal)**7
              s2   = 1.0_ti_p / (rho7+ghal)
              t1   = (1.0_ti_p+dhal)**7 * s1
              t2   = (1.0_ti_p+ghal) * s2
              dt1drho = -7.0_ti_p*(rho+dhal)**6 * t1 * s1
              dt2drho = -7.0_ti_p*rho6 * t2 * s2
              e    = eps2 * t1 * (t2-2.0_ti_p)
              de   = eps2 * (dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv2
#ifndef TINKER_NO_MUTATE
           end if
#endif

           call switch_respa_inl(rik,shortcut,shortheal,s,ds)
           e       = (1-s)*e
           de      = - e*ds - s*de
c
c          use energy switching if near the cutoff distance
c
           if(rik2 > cut2) then ! mask energy switch
              rinv   = 1.0_ti_p/(cut - off)
              r      = (rik - off) * rinv
              r2     = r * r
              r3     = r2 * r
              taper  = r3 * (6*r2 - 15*r + 10)
              dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;

              de     = ( e * dtaper + de * taper ) / rik
              e      =   e * taper
           else
              de = de / rik
           endif

           dedx  =  de * xpos
           dedy  =  de * ypos
           dedz  =  de * zpos
        end
