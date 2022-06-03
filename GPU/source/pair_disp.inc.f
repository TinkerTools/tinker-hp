#ifndef PAIR_DISP_INC
#define PAIR_DISP_INC
#include "tinker_precision.h"
#include "tinker_cudart.h"
#include "convert.f.inc"
#include "switch_respa.f.inc"
c
c
      M_subroutine
     &            duo_disp(r2,xr,yr,zr,ai,ci,ak,ck,rinv,dispshortcut,off
     &                    ,shortheal,mutik,vlambda,i_grp,fgrp,dspscale
     &                    ,e,ded,ver,fea)
      use tintypes  ,only: real3
#if defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
      use utilcu    ,only: f_sqrt,f_exp
#endif
      implicit none
      integer(1) ,intent(in):: mutik
      integer    ,intent(in):: ver,fea,i_grp
      real(t_p)  ,intent(in):: r2,xr,yr,zr,ai,ci,ak,ck,vlambda,rinv
     &           ,dispshortcut,off,shortheal,fgrp
      real(t_p)  ,intent(out):: e
      type(real3),intent(out):: ded

      integer(1)  onei1
      integer     ene,grd,grp,shr,lgr,sfc
      real(t_p)   r,r6,di,dk,di2,di3,expi,expk
     &           ,ai2,ai3,ak2,ak3,dk2,dk3,ti,tk,ti2,tk2
     &           ,damp3,damp5,ddamp,damp
     &           ,di4,di5,de,r1,r22,r3,taper,dtaper,s,ds

      parameter  (ene=__use_ene__, onei1=1
     &           ,grd=__use_grd__
     &           ,shr=__use_shortRange__
     &           ,lgr=__use_longRange__
     &           ,sfc=__use_softcore__
     &           ,grp=__use_groups__)

      r  = f_sqrt(r2)
      r6 = r2**3
      e  = -ci * ck / r6
      IF (iand(ver,grd).ne.0) de = -6.0 * e / r
c
c     find the damping factor for the dispersion interaction
c
      di   = ai * r
      dk   = ak * r
      di2  = di * di
      di3  = di * di2
      expi = f_exp(-di)
      expk = f_exp(-dk)
      if (ai .ne. ak) then
         ai2 = ai * ai
        !ai3 = ai * ai2
         ak2 = ak * ak
        !ak3 = ak * ak2
         dk2 = dk * dk
         dk3 = dk * dk2
         ti  = ak2 / (ak2-ai2)
         tk  = ai2 / (ai2-ak2)
         ti2 = ti * ti
         tk2 = tk * tk
         damp3 = 1.0 - ti2*(1.0+di+0.5*di2)*expi
     &               - tk2*(1.0+dk+0.5*dk2)*expk
     &               - 2.0*ti2*tk*(1.0+di)*expi
     &               - 2.0*tk2*ti*(1.0+dk)*expk
         damp5 = 1.0 - ti2*(1.0+di+0.5*di2+di3/6.0)*expi
     &               - tk2*(1.0+dk+0.5*dk2+dk3/6.0)*expk
     &               - 2.0*ti2*tk*(1.0+di+di2/3.0)*expi
     &               - 2.0*tk2*ti*(1.0+dk+dk2/3.0)*expk
         ddamp = 0.25 * di2 * ti2 * ai * expi
     &                * (r*ai+4.0*tk-1.0)
     &         + 0.25 * dk2 * tk2 * ak * expk
     &                * (r*ak+4.0*ti-1.0)
      else
         di4   = di2 * di2
         di5   = di2 * di3
         damp3 = 1.0 - (1.0+di +0.5*di2+7.0*di3/48.0+di4/48.0)*expi
         damp5 = 1.0 - (1.0+di +0.5*di2+    di3/6.0 +di4/24.0
     &                     +di5/144.0)*expi
         ddamp = ai * expi * (di5-3.0*di3-3.0*di2) / 96.0
      end if
      damp = 1.5*damp5 - 0.5*damp3
 
      !set use of lambda scaling for decoupling or annihilation
      IF (iand(fea,sfc).gt.0) THEN; if(mutik.eq.onei1) then
         IF (iand(ver,grd).ne.0) de = de *vlambda
                                  e =  e *vlambda
      end if; END IF
 
      !apply damping and scaling factors for this interaction
      if (iand(ver,grd).ne.0) de=(de*damp**2+2.0*e*damp*ddamp)*dspscale
                               e=( e*damp**2                 )*dspscale
 
      !use energy switching if near the cutoff distance
      IF (iand(fea,shr).eq.0) THEN
         if (r22 .gt. cut2) then
            r1     = (r - off)*rinv
            r22    = r1 * r1
            r3     = r22 * r1
            taper  = r3 * (6.0*r22 - 15*r1 + 10.0)
            dtaper = 30.0 * (r1*(1.0-r1))*(r1*(1.0-r1)) *rinv;

            IF (iand(ver,grd).ne.0) de = e*dtaper + de*taper
                                     e = e*taper
         end if
      END IF
 
      !scale the interaction based on its group membership
      IF (iand(fea,grp).gt.0) THEN; if(i_grp) then
         e  =  e *fgrp
         if (iand(ver,grd).ne.0) de = de *fgrp
      end if; END IF

      IF (iand(fea,shr+lgr).gt.0)
     &   call switch_respa(r,dispshortcut,shortheal,s,ds)

      IF (iand(fea,shr).gt.0) THEN
         IF (iand(ver,grd).ne.0) de = de *s + e *ds
         IF (iand(ver,ene).ne.0)  e =  e *s
      ELSE IF(iand(fea,lgr).gt.0) THEN
         IF (iand(ver,grd).ne.0) de = de *(1.0-s) - e *ds
         IF (iand(ver,ene).ne.0)  e =  e *(1.0-s)
      END IF

      IF (iand(ver,grd).ne.0) THEN
         ded%x = de * xr/r
         ded%y = de * yr/r
         ded%z = de * zr/r
      END IF
      end subroutine
c
c     "duo_disp_real" evaluates the real space portion of the Ewald
c     summation energy and gradient due to damped dispersion
c     interaction
c
      M_subroutine
     &     duo_disp_real(r2,xr,yr,zr,ai,ci,ak,ck
     &                  ,aewald,dspscale,mutik,vlambda,i_grp,fgrp
     &                  ,dispshortcut,shortheal,e,ded,ver,fea)
      use tintypes  ,only: real3
#if defined(TINKER_CUF) && (defined(SINGLE)||defined(MIXED))
      use utilcu    ,only: f_sqrt,f_exp
#endif
      implicit none
      integer(1) ,intent(in):: mutik
      integer    ,intent(in):: ver,fea,i_grp
      real(t_p)  ,intent(in):: r2,xr,yr,zr,ai,ci,ak,ck,aewald,vlambda
     &           ,fgrp,dspscale,dispshortcut,shortheal
      real(t_p)  ,intent(out):: e
      type(real3),intent(out):: ded

      integer(1)  onei1
      integer     ene,grd,sca,grp,shr,lgr,sfc
      real(t_p)   r,r6,ralpha2,term,expterm,r7
     &           ,di,di2,di3,dk,expi,expk,ai2,ak2,dk2,dk3,ti,ti2,tk,tk2
     &           ,damp,damp3,damp5,ddamp,di4,di5,scale,s,ds
      parameter  (ene=__use_ene__, onei1=1
     &           ,grd=__use_grd__, sca=__use_sca__
     &           ,shr=__use_shortRange__
     &           ,lgr=__use_longRange__
     &           ,sfc=__use_softcore__
     &           ,grp=__use_groups__
     &           )

      r    = f_sqrt(r2)
      r6   = r2**3
      IF (iand(ver,sca).eq.0) then
         ralpha2 = r2 * aewald**2
         term    = 1.0 + ralpha2 + 0.5*ralpha2**2
         expterm = f_exp(-ralpha2)
      END IF
c
c     find the damping factor for the dispersion interaction
c
      r7   = r6 * r
      di   = ai * r
      di2  = di * di
      di3  = di * di2
      dk   = ak * r
      expi = f_exp(-di)
      expk = f_exp(-dk)
      if (ai .ne. ak) then
         ai2   = ai * ai
         ak2   = ak * ak
         dk2   = dk * dk
         dk3   = dk * dk2
         ti    = ak2 / (ak2-ai2)
         tk    = ai2 / (ai2-ak2)
         ti2   = ti * ti
         tk2   = tk * tk
         damp3 = 1.0 - ti2*(1.0+di+0.5*di2)*expi
     &              - tk2*(1.0+dk+0.5*dk2)*expk
     &              - 2.0*ti2*tk*(1.0+di)*expi
     &              - 2.0*tk2*ti*(1.0+dk)*expk
         damp5 = 1.0 - ti2*(1.0+di+0.5*di2+di3/6.0)*expi
     &              - tk2*(1.0+dk+0.5*dk2+dk3/6.0)*expk
     &              - 2.0*ti2*tk*(1.0+di+di2/3.0)*expi
     &              - 2.0*tk2*ti*(1.0+dk+dk2/3.0)*expk
         ddamp = 0.25 * di2 * ti2 * ai * expi * (r*ai+4.0*tk-1.0)
     &         + 0.25 * dk2 * tk2 * ak * expk * (r*ak+4.0*ti-1.0)
      else
         di4   = di2 * di2
         di5   = di2 * di3
         damp3 = 1.0 - (1.0+di+0.5*di2+7.0*di3/48.0+di4/48.0)*expi
         damp5 = 1.0 - (1.0+di+0.5*di2+    di3/6.0 +di4/24.0
     &                 +di5/144.0)*expi
         ddamp = ai * expi * (di5-3.0*di3-3.0*di2) / 96.0
      end if
      damp  = 1.5*damp5 - 0.5*damp3
c
c     apply damping and scaling factors for this interaction
c
      scale = merge(damp**2,dspscale*damp**2,iand(ver,sca).eq.0)
      IF (iand(fea,grp).ne.0) THEN
         if (i_grp)  scale = scale * fgrp
      END IF
c
c     set use of lambda scaling for decoupling or annihilation
c
      IF (iand(fea,sfc).ne.0) THEN; if (mutik.eq.onei1) then
         scale = scale * vlambda
         ddamp = ddamp * vlambda
      end if; END IF

      IF (iand(ver,sca).ne.0) THEN
          e  = -ci * ck * scale / r6
         de  = -6.0*e/r
     &         -2.0*ci*ck*dspscale*damp*ddamp/r6
      ELSE
          e  = -ci * ck * (expterm*term+scale-1.0) / r6
         de  = -6.0*e/r - ci*ck*(-(ralpha2**3)* expterm/r)/r6
     &         -2.0*ci*ck*dspscale*damp*ddamp/r6
      END IF

      IF (iand(fea,shr+lgr).gt.0)
     &   call switch_respa(r,dispshortcut,shortheal,s,ds)

      IF     (iand(fea,shr).ne.0) THEN
         IF  (iand(ver,grd).ne.0) de = de *s + e *ds
         IF  (iand(ver,ene).ne.0)  e =  e *s
      ELSE IF(iand(fea,lgr).ne.0) THEN
         IF  (iand(ver,grd).ne.0) de = de *(1.0-s) - e *ds
         IF  (iand(ver,ene).ne.0)  e =  e *(1.0-s)
      END IF

      IF (iand(ver,grd).ne.0) THEN
         ded%x = de * xr/r
         ded%y = de * yr/r
         ded%z = de * zr/r
      END IF
      end subroutine
#endif
