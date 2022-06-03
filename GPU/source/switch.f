c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine switch  --  get switching function coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "switch" sets the coeffcients used by the fifth and seventh
c     order polynomial switching functions for spherical cutoffs
c
c
#include "tinker_precision.h"
      subroutine switch (mode)
      use cutoff
      use shunt
      use tinheader
      implicit none
      real(t_p) denom,term
      real(t_p) off3,off4,off5
      real(t_p) off6,off7
      real(t_p) cut3,cut4,cut5
      real(t_p) cut6,cut7
      character*11 mode
c
c
c     get the switching window for the current potential type
c
      if (mode(1:3) .eq. 'VDW') then
         off = vdwcut
         cut = vdwtaper
      elseif (mode(1:8) .eq. 'SHORTVDW') then
         off = vdwshortcut
      else if (mode(1:6) .eq. 'REPULS') then
         off = repcut
         cut = reptaper
      elseif (mode(1:11) .eq. 'SHORTREPULS') then
         off = repshortcut
      else if (mode(1:4) .eq. 'DISP') then
         off = dispcut
         cut = disptaper
      elseif (mode(1:10) .eq. 'SHORTDISP') then
         off = dispshortcut
      else if (mode(1:6) .eq. 'CHARGE') then
         off = chgcut
         cut = chgtaper
      else if (mode(1:6) .eq. 'CHGTRN') then
         off = ctrncut
         cut = ctrntaper
      elseif (mode(1:11) .eq. 'SHORTCHGTRN') then
         off = ctrnshortcut
      else if (mode(1:5) .eq. 'MPOLE') then
         off = mpolecut
         cut = mpoletaper
      else if (mode(1:5) .eq. 'EWALD') then
         off = ewaldcut
         cut = ewaldcut
      else if (mode(1:10) .eq. 'SHORTEWALD') then
         off = ewaldshortcut
         cut = ewaldshortcut
      else if (mode(1:5) .eq. 'DEWALD') then
         off = ewaldcut
         cut = ewaldcut
      else if (mode(1:10) .eq. 'SHORTDEWALD') then
         off = ewaldshortcut
         cut = ewaldshortcut
      else
         off = min(vdwcut,chgcut,mpolecut)
         cut = min(vdwtaper,chgtaper,mpoletaper)
      end if
cc
cc     test for replicate periodic boundaries at this cutoff
cc
c      call replica (off)
c
c     set switching coefficients to zero for truncation cutoffs
c
      c0 = 0.0_ti_p
      c1 = 0.0_ti_p
      c2 = 0.0_ti_p
      c3 = 0.0_ti_p
      c4 = 0.0_ti_p
      c5 = 0.0_ti_p
      f0 = 0.0_ti_p
      f1 = 0.0_ti_p
      f2 = 0.0_ti_p
      f3 = 0.0_ti_p
      f4 = 0.0_ti_p
      f5 = 0.0_ti_p
      f6 = 0.0_ti_p
      f7 = 0.0_ti_p
c
c     store the powers of the switching window cutoffs
c
      off2 = off * off
      off3 = off2 * off
      off4 = off2 * off2
      off5 = off2 * off3
      off6 = off3 * off3
      off7 = off3 * off4
      cut2 = cut * cut
      cut3 = cut2 * cut
      cut4 = cut2 * cut2
      cut5 = cut2 * cut3
      cut6 = cut3 * cut3
      cut6 = cut3 * cut3
      cut7 = cut3 * cut4
c
c     get 5th degree multiplicative switching function coefficients
c
      if (cut .lt. off) then
         denom = (off-cut)**5
         c0 = off*off2 * (off2-5.0_ti_p*off*cut+10.0_ti_p*cut2) / denom
         c1 = -30.0_ti_p * off2*cut2 / denom
         c2 = 30.0_ti_p * (off2*cut+off*cut2) / denom
         c3 = -10.0_ti_p * (off2+4.0_ti_p*off*cut+cut2) / denom
         c4 = 15.0_ti_p * (off+cut) / denom
         c5 = -6.0_ti_p / denom
      end if
c
c     get 7th degree additive switching function coefficients
c
      if (cut.lt.off .and. mode(1:6).eq.'CHARGE') then
         term = 9.3_ti_p * cut*off / (off-cut)
         denom = cut7 - 7.0_ti_p*cut6*off + 21.0_ti_p*cut5*off2
     &              - 35.0_ti_p*cut4*off3 + 35.0_ti_p*cut3*off4
     &              - 21.0_ti_p*cut2*off5 + 7.0_ti_p*cut*off6 - off7
         denom = term * denom
         f0 = cut3*off3 * (-39.0_ti_p*cut+64.0_ti_p*off) / denom
         f1 = cut2*off2
     &            * (117.0_ti_p*cut2-100.0_ti_p*cut*off-192.0_ti_p*off2)
     &            / denom
         f2 = cut*off * (-117.0_ti_p*cut3-84.0_ti_p*cut2*off
     &                   +534.0_ti_p*cut*off2+192.0_ti_p*off3) / denom
         f3 = (39.0_ti_p*cut4+212.0_ti_p*cut3*off-450.0_ti_p*cut2*off2
     &            -612.0_ti_p*cut*off3-64.0_ti_p*off4) / denom
         f4 = (-92.0_ti_p*cut3+66.0_ti_p*cut2*off
     &            +684.0_ti_p*cut*off2+217.0_ti_p*off3) / denom
         f5 = (42.0_ti_p*cut2-300.0_ti_p*cut*off-267.0_ti_p*off2)/denom
         f6 = (36.0_ti_p*cut+139.0_ti_p*off) / denom
         f7 = -25.0_ti_p / denom
      end if
      end
c
c
c
      subroutine switch_emtp(doder,r,rin,rout,s,ds)
      use tinheader
      implicit none
      real(t_p) r,rin,rout,s,ds
      logical doder
c
      if (r.le.rin) then
        s = 1.0_ti_p
      else if (r.gt.rin.and.r.le.rout) then
        s = (rout**2-r**2)**2*(rout**2+2*r**2-3*rin**2)/
     $          (rout**2-rin**2)**3
      else if (r.gt.rout) then
        s = 0.0_ti_p
      end if
      if (doder) then
        if (r.le.rin) then
          ds = 0.0_ti_p
        else if (r.gt.rin.and.r.le.rout) then
          ds = 4*r*(3*r**4-3*r**2*rout**2-3*rin**2*r**2+
     $          3*rin**2*rout**2)/
     $         (rout**2-rin**2)**3
        else if (r.gt.rout) then
          ds = 0.0_ti_p
        end if
      end if
      end
c
c
      subroutine switch_respa(r,rc,lambda,s,ds)
      use tinheader ,only: ti_p
      implicit none
      real(t_p),intent(in)   :: r,rc,lambda
      real(t_p),intent(inout):: s,ds
      real(t_p) u,du
c
      if (r.le.(rc-lambda)) then
        s = 1.0_ti_p
        ds = 0.0_ti_p
      else if (((rc-lambda).le.r).and.(r.le.rc)) then
        u = (1.0_ti_p/lambda)*(r-rc+lambda)
        du = 1.0_ti_p/lambda
        s = 1 + u**3*(15*u - 6*u**2 - 10)
        ds = 3*du*u**2*(15*u - 6*u**2 - 10) + 
     $   u**3*(15*du - 12*du*u)
      else if (r.ge.rc) then
        s = 0.0_ti_p
        ds = 0.0_ti_p
      end if
      end
c
