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
      subroutine switch (mode)
      use cutoff
      use shunt
      implicit none
      real*8 denom,term
      real*8 off3,off4,off5
      real*8 off6,off7
      real*8 cut3,cut4,cut5
      real*8 cut6,cut7
      character*6 mode
c
c
c     get the switching window for the current potential type
c
      if (mode(1:3) .eq. 'VDW') then
         off = vdwcut
         cut = vdwtaper
      else if (mode(1:6) .eq. 'CHARGE') then
         off = chgcut
         cut = chgtaper
      else if (mode(1:5) .eq. 'MPOLE') then
         off = mpolecut
         cut = mpoletaper
      else if (mode(1:5) .eq. 'EWALD') then
         off = ewaldcut
         cut = ewaldcut
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
      c0 = 0.0d0
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      c4 = 0.0d0
      c5 = 0.0d0
      f0 = 0.0d0
      f1 = 0.0d0
      f2 = 0.0d0
      f3 = 0.0d0
      f4 = 0.0d0
      f5 = 0.0d0
      f6 = 0.0d0
      f7 = 0.0d0
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
         c0 = off*off2 * (off2-5.0d0*off*cut+10.0d0*cut2) / denom
         c1 = -30.0d0 * off2*cut2 / denom
         c2 = 30.0d0 * (off2*cut+off*cut2) / denom
         c3 = -10.0d0 * (off2+4.0d0*off*cut+cut2) / denom
         c4 = 15.0d0 * (off+cut) / denom
         c5 = -6.0d0 / denom
      end if
c
c     get 7th degree additive switching function coefficients
c
      if (cut.lt.off .and. mode(1:6).eq.'CHARGE') then
         term = 9.3d0 * cut*off / (off-cut)
         denom = cut7 - 7.0d0*cut6*off + 21.0d0*cut5*off2
     &              - 35.0d0*cut4*off3 + 35.0d0*cut3*off4
     &              - 21.0d0*cut2*off5 + 7.0d0*cut*off6 - off7
         denom = term * denom
         f0 = cut3*off3 * (-39.0d0*cut+64.0d0*off) / denom
         f1 = cut2*off2
     &           * (117.0d0*cut2-100.0d0*cut*off-192.0d0*off2) / denom
         f2 = cut*off * (-117.0d0*cut3-84.0d0*cut2*off
     &                   +534.0d0*cut*off2+192.0d0*off3) / denom
         f3 = (39.0d0*cut4+212.0d0*cut3*off-450.0d0*cut2*off2
     &            -612.0d0*cut*off3-64.0d0*off4) / denom
         f4 = (-92.0d0*cut3+66.0d0*cut2*off
     &            +684.0d0*cut*off2+217.0d0*off3) / denom
         f5 = (42.0d0*cut2-300.0d0*cut*off-267.0d0*off2) / denom
         f6 = (36.0d0*cut+139.0d0*off) / denom
         f7 = -25.0d0 / denom
      end if
      return
      end
c
c
c
      subroutine switch_emtp(doder,r,rin,rout,s,ds)
      implicit none
      real*8 r,rin,rout,s,ds
      logical doder
c
      if (r.le.rin) then
        s = 1.0d0
      else if (r.gt.rin.and.r.le.rout) then
        s = (rout**2-r**2)**2*(rout**2+2*r**2-3*rin**2)/
     $          (rout**2-rin**2)**3
      else if (r.gt.rout) then
        s = 0.0d0
      end if
      if (doder) then
        if (r.le.rin) then
          ds = 0.0d0
        else if (r.gt.rin.and.r.le.rout) then
          ds = 4*r*(-3*r**2+rin**2)/
     $         (rout**2-rin**2)**3
        else if (r.gt.rout) then
          ds = 0.0d0
        end if
      end if
      return
      end
