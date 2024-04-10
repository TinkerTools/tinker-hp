#ifndef PAIR_ELJ_INC
#define PAIR_ELJ_INC
#include "tinker_cudart.h"
#include "convert.f.inc"
#include "switch_respa.f.inc"

c
c     ------------------------------------------------------
c     Lennard-Jones pairwise interaction computation routine
c     ------------------------------------------------------
c
      M_subroutine
     &           duo_lj_(rik2,rik,rik_i,rv,eps,cut2
     &                  ,rinv,off,sheal,ugrp,fgrp,mutik
     &                  ,sck,sct,scs,scalpha,galpha,ulamdyn
     &                  ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                  ,delambdav,e,de,ver,fea)
!$acc routine
        use tinheader ,only: one1
        use tintypes  ,only: real3
#if defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
        use utilcu    ,only: f_sqrt
#endif
        implicit none
        logical  ,intent(in ):: ulamdyn, ugrp
        integer  ,intent(in ):: ver,fea
        integer(1),intent(in ):: mutik
        real(t_p),intent(in ):: rik2,rik,rik_i,rv,eps,cut2
     &           ,rinv,off,sheal,fgrp,sct,sck,scs
     &           ,scalpha,lambda,vlambda,lambdavt,lambdavt1
     &           ,galpha,glamb
        real(t_p),intent(out):: e,de,delambdav

        integer    ene,grd,shr,sfc,sca,grp,lgr,lbd
        real(t_p)  p6,taper,dtaper
        real(t_p)  r,r3,ds,s,rho,rhok,rvk
     &            ,gsc,dgscrho,rvogsc6,dgsclambda,evdw,devdwgsc
        parameter(ene=__use_ene__
     &           ,grd=__use_grd__
     &           ,shr=__use_shortRange__
     &           ,lgr=__use_longRange__
     &           ,sfc=__use_softcore__
     &           ,grp=__use_groups__
     &           ,sca=__use_sca__
     &           ,lbd=__use_lambdadyn__)

        IF (iand(fea,sfc).ne.0) delambdav=0.0
        IF (IAND(fea,sfc).NE.0.and.mutik.eq.one1) then
           rho  = rik/rv
           rhok = rho**sck
           rvk  = rv**sck
           !Softcore expression 
           gsc  = rv *( galpha * glamb**scs + rhok )**(1.0/sck)
           rvogsc6 = (rv/gsc)**6

           !Softcore derivatives w.r.t rho, r and lambda_v 
           dgscrho = rvk * (rho**(sck-1.0))* gsc ** (1.0-sck)

           dgsclambda = -(scs/sck)*rvk*galpha*(glamb**(scs-1.0))*
     &       gsc**(1.0-sck)
c
           evdw = eps*rvogsc6*(rvogsc6 - 2.0)

c           if((lambda.gt.0.0).and.(lambda.lt.1.0)) then
              devdwgsc = eps*12.0*rvogsc6*(1.0-rvogsc6)/gsc
c           else
c              devdwgsc = 0.0
c           end if
c
           e  =  lambdavt*evdw
           IF (iand(ver,grd).ne.0)de = (lambdavt*dgscrho*devdwgsc)/rv
           delambdav = lambdavt1*evdw + lambdavt*dgsclambda*devdwgsc
        ELSE
           r   = rv*rv/rik2
           p6  = r*r*r
           e   = eps * (p6-2.0)*p6
           IF (iand(ver,grd).ne.0)de = eps * (p6-1.0)*p6*(-12.0*rik_i)
        END IF
c
c       use energy switching if near the cutoff distance
c
        IF (iand(fea,shr).eq.0) THEN
           if (rik2 .gt. cut2) then
              r      = (rik - off) * rinv
              !r3     = r*r*r
              taper  = r*r*r * (6*r*r - 15.0*r + 10.0)
              dtaper = 30* (r*(1.0-r))*(r*(1.0-r)) *rinv;

              IF (iand(ver,grd).gt.0) de = e*dtaper + de*taper
                                       e = e* taper
           end if
        END IF
c
c     scale the interaction based on its group membership
c
        if (IAND(fea,grp).NE.0.and.ugrp) then
            e =  e * fgrp
           IF (iand(ver,grd).ne.0)de = de * fgrp
        end if
c
c     use energy switching if close the cutoff distance (at short range)
c
        IF (iand(fea,shr+lgr).ne.0)
     &     call switch_respa_inl(rik,off,sheal,s,ds)

        IF (iand(fea,lgr).ne.0) THEN
           IF(iand(ver,grd).ne.0) de = -e*ds + (1.0-s)*de
           IF(iand(ver,ene).ne.0)  e = (1.0-s)*e
           IF (iand(ver,grd).ne.0.and.ulamdyn) THEN
             delambdav = delambdav*(1.0-s)
           END IF
        ELSE if(iand(fea,shr).ne.0) THEN
           IF(iand(ver,grd).ne.0) de = e*ds + de*s
           IF(iand(ver,ene).ne.0) e  = e * s
           IF (iand(ver,grd).ne.0.and.ulamdyn) THEN
             delambdav = delambdav*s
           END IF
        END IF
      end

      M_subroutine
     &            duo_lj(rik2,xr,yr,zr,rv,eps,cut2
     &                  ,rinv,off,sheal,ugrp,fgrp,mutik
     &                  ,sck,sct,scs,scalpha,galpha,ulamdyn
     &                  ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                  ,delambdav,e,ded,ver,fea)
!$acc routine
        use tinheader ,only: ti_p,re_p
        use tintypes  ,only: real3
#if defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
        use utilcu    ,only: f_sqrt
#endif
        implicit none
        logical  ,intent(in ):: ulamdyn, ugrp
        integer  ,intent(in ):: ver,fea
        integer(1),intent(in ):: mutik
        real(t_p),intent(in ):: rik2,xr,yr,zr,rv,eps,cut2
     &           ,rinv,off,sheal,fgrp,sct,sck,scs
     &           ,scalpha,lambda,vlambda,lambdavt,lambdavt1
     &           ,galpha,glamb
        real(t_p),intent(out):: e,delambdav
        type(real3),intent(out):: ded

        integer grd
        real(t_p) rik,rik_i,de
        parameter( grd=__use_grd__ )
        rik   = f_sqrt(rik2)
        rik_i = rik**(-1)

        call duo_lj_(rik2,rik,rik_i,rv,eps,cut2
     &              ,rinv,off,sheal,ugrp,fgrp,mutik
     &              ,sck,sct,scs,scalpha,galpha,ulamdyn
     &              ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &              ,delambdav,e,de,ver,fea)
c
c       find the chain rule terms for derivative components
c
        IF (iand(ver,grd).ne.0) THEN
           de    = de * rik_i
           ded%x = de * xr
           ded%y = de * yr
           ded%z = de * zr
        END IF
      end subroutine
#endif
