#ifndef PAIR_EHAL_INC
#define PAIR_EHAL_INC
#include "tinker_cudart.h"
#include "convert.f.inc"
#include "switch_respa.f.inc"

#if 1
#  define GHAL 0.12_ti_p
#  define DHAL 0.07_ti_p
#  define SCEXP 5.0_ti_p
#  define SCAPLHA 0.7_ti_p
#else
#  define GHAL ghal
#  define DHAL dhal
#  define SCEXP scexp
#  define SCALPHA scalpha
#endif

        M_subroutine
     &             ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,vscale
     &                         ,cut2,rinv,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,mutik
     &                         ,e,dedx,dedy,dedz
     &                         ,ver,fea)
!$acc routine
#if  defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
           use utilcu    ,only: f_sqrt
#endif
           use tinheader ,only: ti_p,re_p
           implicit none
           integer   ,intent(in) ::ver,fea
           real(ti_p),intent(in) ::xpos,ypos,zpos,rik2
           real(ti_p) rv2,eps2,vscale
           real(ti_p),intent(in) ::cut2,rinv,off,ghal,dhal
           real(ti_p),intent(in) ::scexp,vlambda,scalpha
           integer(1),intent(in) ::mutik
           real(ti_p),intent(out)::e
           real(ti_p),intent(out)::dedx,dedy,dedz

           integer(1) one1
           real(ti_p) rik,r,r2,r3
           !real(ti_p) rik3,rik4,rik5,rik6,rik7,rv7,rv7orho
           real(ti_p) dtau,gtau,tau,tau7,de,taper,dtaper
           real(ti_p) rho,rho6,rho7
           real(ti_p) scal,s1,s2,t1,t2,dt1drho,dt2drho
           parameter(one1=1)
c
c          compute the energy contribution for this interaction
c
           rik  = f_sqrt(rik2)
           rho  = rik / rv2
           rho6 = rho**6
           rho7 = rho6 * rho
           if ( iand(fea,__use_softcore__).gt.0 ) then  !CONST TEST
              eps2 = eps2 * vscale
     &             * merge(vlambda**scexp,1.0_ti_p,mutik.eq.one1)
              scal = scalpha * (1.0_ti_p-vlambda)**2
           else
              eps2 = eps2 * vscale
              scal = 0.0
           end if
           s1   = (scal+(rho+dhal)**7)**(-1)
           s2   = (scal+rho7+ghal)**(-1)
           t1   = (1.0_ti_p+dhal)**7 * s1
           t2   = (1.0_ti_p+ghal) * s2
           if (iand(ver,__use_ene__).gt.0)  ! CONST TEST
     &        e    = eps2 * t1 * (t2-2.0_ti_p)
           if (iand(ver,__use_grd__).gt.0) then ! CONST TEST
              dt1drho = -7.0_ti_p*(rho+dhal)**6 * t1 * s1
              dt2drho = -7.0_ti_p*rho6 * t2 * s2
              de   = eps2 * (dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv2
           end if
           ! Inadequate for single precision compute (rik**6 overflow risk)
c          rik     =  f_sqrt(rik2)
c          rik6    =  rik2 * rik2 * rik2
c          rik7    =  rik2 * rik2 * rik2 * rik
c          dtau    =  (rik + dhal * rv2)**(-1) !tau / (dhal + 1.0)
c          rv7     =  rv2**7
c          rv7orho =  rv7 / (rik7 + ghal * rv7)
c          tau7    =  (dtau * (dhal + 1.0))**7
c          !tau7    =  tau ** 7.0
c          gtau    =  eps2 * tau7 * rik6
c    &              * (ghal + 1.0) * rv7orho * rv7orho * vscale
c          e       =  eps2 * tau7 * rv7
c    &              * ((ghal + 1.0) * rv7orho - 2.0) * vscale
c          de      = - 7.0 * (dtau*e + gtau)
c          if (iand(fea,__use_shortRange__+__use_longRange__).gt.0) then !CONST TEST
c             call switch_respa_inl(rik,shortcut,shortheal,s,ds)
c             if (iand(ver,__use_grd__).gt.0) de = ( e*ds + de*s )/rik
c             if (iand(ver,__use_ene__).gt.0) e  = e*s
c          end if
c
c          use energy switching if near the cutoff distance
c
           if (iand(fea,__use_shortRange__).eq.0) then ! CONST TEST (not shortrange)
           if (rik2 > cut2) then ! mask energy switch
              ! Single precision overflow due to rik6
c             taper  =  c5 * rik2*rik2*rik + c4 * rik2*rik2
c    &                + c3 * rik2*rik      + c2 * rik2
c    &                + c1 * rik           + c0
c             dtaper =  5.0 * c5 * rik2*rik2
c    &                + 4.0 * c4 * rik2*rik
c    &                + 3.0 * c3 * rik2
c    &                + 2.0 * c2 * rik
c    &                +       c1

              r      = (rik - off)*rinv
              r2     = r * r
              r3     = r2 * r
              taper  = r3 * (6*r2 - 15*r + 10)
              dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;

              if (iand(ver,__use_grd__).gt.0) de = ( e*dtaper + de*taper ) / rik !CONST TEST
              if (iand(ver,__use_ene__).gt.0)  e = e * taper  ! CONST TEST
           else
              if (iand(ver,__use_grd__).gt.0)de = de / rik ! CONST TEST
           end if
           end if

           if (iand(ver,__use_grd__).gt.0) then   ! CONST TEST
              dedx  =  de * xpos
              dedy  =  de * ypos
              dedz  =  de * zpos
           end if
        end

        M_subroutine
     &             duo_hal(xpos,ypos,zpos,rik2,rv2,eps2,vscale,fgrp
     &                    ,rinv,cut2,shortcut,off,shortheal,ghal,dhal
     &                    ,scexp,vlambda,scalpha,dt1lam,dt2lam,mutik
     &                    ,ulamdyn,ugrp,delambdav,e,dedx,dedy,dedz
     &                    ,ver,fea)
!$acc routine
#if  defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
           use utilcu    ,only: f_sqrt
#endif
           use tinheader ,only: ti_p,re_p,one1
           implicit none
           integer   ,intent(in):: ver,fea
           integer(1),intent(in):: mutik
           logical   ,intent(in):: ulamdyn,ugrp
           real(t_p) ,intent(in):: xpos,ypos,zpos,rik2,ghal,dhal
     &               ,rv2,vscale,cut2,shortcut,off,shortheal,fgrp
     &               ,rinv,scexp,vlambda,scalpha,dt1lam,dt2lam
           real(t_p)  eps2
           real(t_p) ,intent(out)::delambdav,e,dedx,dedy,dedz

           integer    ene,grd,shr,lgr,sfc,lbd,grp
           real(t_p)  rik,r,r2,r3
     &               ,dtau,gtau,tau,tau7,de,taper,dtaper
     &               ,rho,rho6,rho7,s,ds
     &               ,scal,s1,s2,t1,t2,dt1drho,dt2drho
           !real(t_p) rik3,rik4,rik5,rik6,rik7,rv7,rv7orho
           parameter (ene=__use_ene__
     &               ,grd=__use_grd__
     &               ,sfc=__use_softcore__
     &               ,lbd=__use_lambdadyn__
     &               ,grp=__use_groups__
     &               ,shr=__use_shortRange__
     &               ,lgr=__use_longRange__
     &               )
c          real(t_p) ds1dlambda,ds2dlambda,dt1dlambda,dt2dlambda
c
c          compute the energy contribution for this interaction
c
           rik  = f_sqrt(rik2)
           rho  = rik / rv2
           rho6 = rho**6
           rho7 = rho6 * rho
           IF ( iand(fea,sfc).ne.0.and.mutik.eq.one1 ) THEN  ! Apply a Soft-Core on interaction
              eps2 = eps2 * vscale
              scal = scalpha*(1.0-vlambda)**2
              r2   = vlambda**scexp
              s1   = (scal+(rho+DHAL)**7)**(-1)
              s2   = (scal+rho7+GHAL)**(-1)
              t1   = (1.0_ti_p+DHAL)**7 * s1
              t2   = (1.0_ti_p+GHAL) * s2
              e    = eps2*r2* t1* (t2-2.0_ti_p)
              if (IAND(fea,lbd).NE.0.and.ulamdyn) then ! Lambda dynamic
                delambdav = eps2*scexp*vlambda**(scexp-1.0)*t1*(t2-2.0)+
     &                      eps2*r2*
     &                   (dt1lam*s1*s1*(t2-2.0)+t1*dt2lam*s2*s2)
              end if
              IF (iand(ver,grd).ne.0) THEN
                 dt1drho = -7.0_ti_p*(rho+DHAL)**6 * t1 * s1
                 dt2drho = -7.0_ti_p*rho6 * t2 * s2
                 de   = eps2*r2* (dt1drho*(t2-2.0_ti_p)+t1*dt2drho)/ rv2
              END IF
           ELSE
              eps2 = eps2 * vscale
              s1   = ((rho+DHAL)**7)**(-1)
              s2   = (rho7+GHAL)**(-1)
              t1   = (1.0_ti_p+DHAL)**7 * s1
              t2   = (1.0_ti_p+GHAL) * s2
              e    = eps2* t1* (t2-2.0_ti_p)
              IF (iand(fea,lbd).ne.0) delambdav = 0
              IF (iand(ver,grd).ne.0) THEN
                 dt1drho = -7.0_ti_p*(rho+DHAL)**6 * t1 * s1
                 dt2drho = -7.0_ti_p*rho6 * t2 * s2
                 de   = eps2* (dt1drho*(t2-2.0_ti_p)+t1*dt2drho)/ rv2
              END IF
           END IF

           IF (iand(fea,shr+lgr).ne.0) THEN ! Apply range shift
              call switch_respa_inl(rik,shortcut,shortheal,s,ds)
              IF (iand(fea,shr).ne.0)  THEN ! Short Range
                 IF (iand(ver,grd).ne.0) de = ( e*ds + de*s )/rik
                                          e = e*s
              ELSE
                 IF (iand(ver,grd).ne.0) de = ( -e*ds + (1-s)*de )
                                          e = (1-s)*e
              END IF
           END IF
c
c     scale the interaction based on its group membership
c
           if (IAND(fea,grp).NE.0.and.ugrp) then
               e =  e * fgrp
              IF (IAND(ver,grd).NE.0) de = de * fgrp
           end if
c
c          use energy switching if near the cutoff distance
c
           IF (IAND(fea,shr).EQ.0) THEN ! CONST TEST (not short range)
c
           if (rik2 > cut2) then ! mask energy switch
              r      = (rik - off)*rinv
              r2     = r * r
              r3     = r2 * r
              taper  = r3 * (6*r2 - 15*r + 10)
              dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;

              IF (IAND(ver,grd).NE.0) de =(e*dtaper + de*taper )/rik
              IF (IAND(ver,ene).NE.0)  e = e * taper
           else
              IF (IAND(ver,grd).NE.0) de = de / rik
           end if
c
           END IF

           IF (iand(ver,grd).ne.0) THEN  ! Compute Gradient
              dedx  =  de * xpos
              dedy  =  de * ypos
              dedz  =  de * zpos
           END IF
        end
#endif
