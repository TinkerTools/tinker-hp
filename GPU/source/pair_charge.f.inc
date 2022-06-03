#ifndef PAIR_CHARGE_INC
#define PAIR_CHARGE_INC
#include "tinker_cudart.h"
#include "convert.f.inc"
#include "switch_respa.f.inc"
      M_subroutine
     &          charge_couple(r2,xr,yr,zr,ebuffer,fik_,fik,aewald
     &                       ,cscale,mutik,use_lambdadyn
     &                       ,shortheal,chgshortcut
     &                       ,elambda,delambdae,e,ded,ver,fea)
!$acc routine
        use tinheader ,only: ti_p
        use math      ,only: sqrtpi
        use tintypes  ,only: real3,mdyn3_r
#ifdef TINKER_CUF
        use utilcu ,only: f_erfc
#  if defined(SINGLE)||defined(MIXED)
        use utilcu ,only: f_sqrt,f_exp
#  endif
#endif
        implicit none
        integer(1) ,intent(in ):: mutik
        integer    ,intent(in ):: fea,ver
        logical    ,intent(in ):: use_lambdadyn
        real(t_p)  ,intent(in ):: r2,xr,yr,zr,ebuffer
     &             ,fik,fik_,aewald,cscale,chgshortcut,shortheal
     &             ,elambda
        real(t_p)  ,intent(inout):: e,delambdae
        type(real3),intent(out):: ded

        integer(1)  onei1,twoi1
        integer     ene,grd,shr,lgr,sca,lbd,ace
        real(t_p)   r,rb,rb2,rew,de,erfterm,exp2a,s,ds
        parameter  (ene=__use_ene__, onei1=1, twoi1=2
     &             ,grd=__use_grd__
     &             ,shr=__use_shortRange__
     &             ,lgr=__use_longRange__
     &             ,sca=__use_sca__
     &             ,lbd=__use_lambdadyn__
     &             ,ace=shr+lgr+lbd
     &             )

        r   = f_sqrt(r2)
        rb  = r + ebuffer
        rb2 = rb * rb

        IF (iand(ver,sca).eq.0) THEN
           rew     = aewald * r
           exp2a   = f_exp(-rew**2)
           erfterm = f_erfc (rew)
           IF (iand(fea,ace).ne.0.or.iand(ver,ene).ne.0)
     &        e = ((fik/rb)*(erfterm))
           IF (iand(ver,grd).ne.0)de = -fik * ((erfterm)/rb2
     &                + (2.0*aewald/sqrtpi)*exp2a/rb)
        ELSE
           IF (iand(fea,ace).ne.0.or.iand(ver,ene).ne.0)
     &        e = ((fik/rb)*(cscale))
           IF (iand(ver,grd).ne.0)de = -fik     * (cscale)/rb2
        END IF

        IF (iand(fea,shr+lgr).ne.0) THEN
           call switch_respa_inl(r,chgshortcut,shortheal,s,ds)
           IF (iand(fea,shr).ne.0) THEN
              de = de*s + e*ds
              IF (iand(fea,ace).ne.0) e = e*s
           ELSE
              de = de *(1.0-s) - e*ds
              IF (iand(fea,ace).ne.0) e = e *(1.0-s)
           END IF
        END IF
c
c       lambdadyn
c
        IF (iand(fea,lbd).ne.0) then

        if (use_lambdadyn) then
           delambdae = 0.0
           if      (mutik.eq.twoi1) then
              IF (iand(ver,sca).ne.0) THEN
                 delambdae = 2.0*elambda*(fik_/rb)*erfterm
              ELSE
                 delambdae = 2.0*elambda*(fik_/rb)*cscale
              END IF
              IF     (iand(fea,shr).ne.0) THEN
                 delambdae = delambdae*s
              ELSE IF(iand(fea,lgr).ne.0) THEN
                 delambdae = delambdae*(1.0-s)
              END IF
           else if (mutik.eq.onei1) then
              IF (iand(ver,sca).eq.0) THEN
                 delambdae = (fik_/rb)*erfterm
              ELSE
                 delambdae = (fik_/rb)*cscale
              END IF
           end if
        end if

        END IF
c
c     form the chain rule terms for derivative expressions
c
        IF (iand(ver,grd).ne.0) THEN
           de     = de / r
           ded%x  = de * xr
           ded%y  = de * yr
           ded%z  = de * zr
        END IF
      end subroutine

      M_subroutine
     &            duo_chg1(r,ebuffer,fik,aewald,aewaldop,cscale
     &                    ,e,de,ver,fea)
!$acc routine
        use tinheader ,only: ti_p
        use math      ,only: sqrtpi
        use tintypes  ,only: real3,mdyn3_r
#ifdef TINKER_CUF
        use utilcu ,only: f_erfc,f_inv
#  if defined(SINGLE)||defined(MIXED)
        use utilcu ,only: f_sqrt,f_exp
#  endif
#endif
        implicit none
        real(t_p)  ,intent(in ):: ebuffer
     &             ,fik,aewald,aewaldop,cscale
        integer    ,intent(in ):: ver,fea
        real(t_p)  ,intent(inout):: e
        real(t_p)  ,intent(out):: de

        integer     ene,grd,sca
        real(t_p)   r,rb,rew,erfterm,exp2a
        parameter( ene=__use_ene__
     &           , grd=__use_grd__
     &           , sca=__use_sca__
     &           )

        rb  = (r + ebuffer)**(-1)
        if (iand(ver,sca).eq.0) then
           rew     = aewald * r
           exp2a   = f_exp(-rew**2)
           erfterm = f_erfc(rew)
           if(iand(ver,ene).ne.0)  e  = e + (fik*rb*erfterm)
           if(iand(ver,grd).ne.0) de  = -fik*rb*(erfterm*rb
     &                                + aewaldop*exp2a)
        else
           if(iand(ver,ene).ne.0)  e  = e + (fik*rb*cscale)
           if(iand(ver,grd).ne.0) de  = -fik*cscale*rb*rb
        end if
      end subroutine

#endif
