c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine evcorr  --  long range vdw energy correction  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "evcorr" computes the long range van der Waals correction
c     to the energy via numerical integration
c
c     literature reference:
c
c     M. P. Allen and D. J. Tildesley, "Computer Simulation of
c     Liquids", Oxford University Press, 1987, Section 2.8
c
c
#include "tinker_precision.h"
      subroutine evcorr (elrc)
      use sizes
      use bound
      use boxes
      use math
      use shunt
      use tinheader,only:ti_p
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,it,kt
      integer nstep,ndelta
      real(t_p) zero,one,two
      real(r_p) elrc,etot
      real(t_p) range,rdelta
      real(t_p) termi,termik
      real(t_p) e,eps
      real(t_p) offset,taper
      real(t_p) rv,rv2,rv6,rv7
      real(t_p) r,r2,r3,r4
      real(t_p) r5,r6,r7
      real(t_p) p,p6,p12
      real(t_p) rho,tau,tau7
      character*10 mode
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p)
c
c
c     zero out the long range van der Waals correction
c
      elrc = 0
c
c     only applicable if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     set number of steps and range for numerical integration
c
      nstep = 2
      range = 100.0_ti_p
      ndelta = int(real(nstep,t_p)*(range-cut))
      rdelta = (range-cut) / real(ndelta,t_p)
      offset = cut - 0.5_ti_p*rdelta
c
c     find the van der Waals energy via double loop search
c
      elrc = zero
      do i = 1, nvt
         it = ivt(i)
         termi = two * pi * real(jvt(i),t_p)
         do k = 1, nvt
            kt = ivt(k)
            termik = termi * real(jvt(k),t_p)
            rv = radmin(kt,it)
            eps = epsilon(kt,it)
            rv2 = rv * rv
            rv6 = rv2 * rv2 * rv2
            rv7 = rv6 * rv
            etot = zero
            do j = 1, ndelta
               r = offset + real(j,t_p)*rdelta
               r2 = r * r
               r3 = r2 * r
               r6 = r3 * r3
               r7 = r6 * r
               e = zero
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  p6 = rv6 / r6
                  p12 = p6 * p6
                  e = eps * (p12 - two*p6)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
               rho = r7 + ghal*rv7
               tau = (dhal+one) / (r+dhal*rv)
               tau7 = tau**7
               e = eps * rv7 * tau7 * ((ghal+one)*rv7/rho-two)
c               else if (vdwtyp.eq.'BUCKINGHAM' .or.
c     &                  vdwtyp.eq.'MM3-HBOND') then
c                  p = sqrt(rv2/r2)
c                  p6 = rv6 / r6
c                  expterm = abuck * exp(-bbuck/p)
c                  e = eps * (expterm - cbuck*p6)
               end if
               if (r .lt. off) then
                  r4 = r2 * r2
                  r5 = r2 * r3
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  e = e * (one-taper)
               end if
               etot = etot + e*rdelta*r2
            end do
            elrc = elrc + termik*etot
         end do
      end do
      elrc = elrc / volbox
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine evcorr1  --  long range vdw energy & virial  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "evcorr1" computes the long range van der Waals correction
c     to the energy and virial via numerical integration
c
c     literature reference:
c
c     M. P. Allen and D. J. Tildesley, "Computer Simulation of
c     Liquids", Oxford University Press, 1987, Section 2.8
c
c
      subroutine evcorr1 (elrc,vlrc)
      use sizes
      use bound
      use boxes
      use math
      use shunt
      use tinheader ,only: ti_p
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,it,kt
      integer nstep,ndelta
      real(t_p) zero,one,two
      real(r_p) elrc,vlrc
      real(t_p) etot,vtot
      real(t_p) range,rdelta
      real(t_p) termi,termik
      real(t_p) e,de,eps
      real(t_p) offset
      real(t_p) taper,dtaper
      real(t_p) rv,rv2,rv6,rv7
      real(t_p) r,r2,r3,r4
      real(t_p) r5,r6,r7
      real(t_p) p,p6,p12
      real(t_p) rho,tau,tau7
      real(t_p) dtau,gtau
      real(t_p) rvterm,expterm
      character*10 mode
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p)
c
c
c     zero out the long range van der Waals corrections
c
      elrc = 0
      vlrc = 0
c
c     only applicable if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     set number of steps and range for numerical integration
c
      nstep = 2
      range = 100.0_ti_p
      ndelta = int(real(nstep,t_p)*(range-cut))
      rdelta = (range-cut) / real(ndelta,t_p)
      offset = cut - 0.5_ti_p*rdelta
c
c     find the van der Waals energy via double loop search
c
      do i = 1, nvt
         it = ivt(i)
         termi = two * pi * real(jvt(i),t_p)
         do k = 1, nvt
            kt = ivt(k)
            termik = termi * real(jvt(k),t_p)
            rv = radmin(kt,it)
            eps = epsilon(kt,it)
            rv2 = rv * rv
            rv6 = rv2 * rv2 * rv2
            rv7 = rv6 * rv
            etot = zero
            vtot = zero
            do j = 1, ndelta
               r = offset + real(j,t_p)*rdelta
               r2 = r * r
               r3 = r2 * r
               r6 = r3 * r3
               r7 = r6 * r
               e = zero
               de = zero
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  p6 = rv6 / r6
                  p12 = p6 * p6
                  e = eps * (p12 - two*p6)
                  de = eps * (p12-p6) * (-12.0_ti_p/r)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
              rho = r7 + ghal*rv7
              tau = (dhal+one) / (r+dhal*rv)
                  tau7 = tau**7
                  dtau = tau / (dhal+one)
                  gtau = eps*tau7*r6*(ghal+one)*(rv7/rho)**2
              e = eps * rv7 * tau7 * ((ghal+one)*rv7/rho-two)
               de = -7.0_ti_p * (dtau*e+gtau)
c               else if (vdwtyp.eq.'BUCKINGHAM' .or.
c     &                  vdwtyp.eq.'MM3-HBOND') then
c                  p = sqrt(rv2/r2)
c                  p6 = rv6 / r6
c                  rvterm = -bbuck / rv
c                  expterm = abuck * exp(-bbuck/p)
c                  e = eps * (expterm - cbuck*p6)
c                  de = eps * (rvterm*expterm+6.0_ti_p*cbuck*p6/r)
               end if
               if (r .lt. off) then
                  r4 = r2 * r2
                  r5 = r2 * r3
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0_ti_p*c5*r4 + 4.0_ti_p*c4*r3
     &                        + 3.0_ti_p*c3*r2 + two*c2*r + c1
                  de = de*(one-taper) - e*dtaper
                  e = e*(one-taper)
               end if
               etot = etot + e*rdelta*r2
               vtot = vtot + de*rdelta*r3
            end do
            elrc = elrc + termik*etot
            vlrc = vlrc + termik*vtot
         end do
      end do
      elrc = elrc / volbox
      vlrc = vlrc / (3.0_re_p*volbox)
      return
      end
c
c
c
      subroutine evcorr1gpu (elrc,vlrc)
      use sizes
      use bound
      use boxes
      use math
      use shunt
      use tinheader ,only: ti_p
      use utilgpu   ,only: def_queue
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,it,kt
      integer nstep,ndelta
      integer vdwtyp_i
      real(t_p) zero,one,two
      real(r_p) elrc,vlrc
      real(t_p) etot,vtot
      real(t_p) range,rdelta
      real(t_p) termi,termik
      real(t_p) e,de,eps
      real(t_p) offset
      real(t_p) taper,dtaper
      real(t_p) rv,rv2,rv6,rv7
      real(t_p) r,r2,r3,r4
      real(t_p) r5,r6,r7
      real(t_p) p,p6,p12
      real(t_p) rho,tau,tau7
      real(t_p) dtau,gtau
      real(t_p) rvterm,expterm
      character*10 mode
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p)

c
c     only applicable if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     zero out the long range van der Waals corrections
c
!$acc data present(elrc,vlrc) async(def_queue)
!$acc kernels async(def_queue)
      elrc = zero
      vlrc = zero
!$acc end kernels
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     set number of steps and range for numerical integration
c
      if (vdwtyp .eq. 'LENNARD-JONES') then
         vdwtyp_i=0
      else if (vdwtyp .eq. 'BUFFERED-14-7') then
         vdwtyp_i=1
      end if

      nstep = 2
      range = 100.0_ti_p
      ndelta = int(real(nstep,t_p)*(range-cut))
      rdelta = (range-cut) / real(ndelta,t_p)
      offset = cut - 0.5_ti_p*rdelta
c
c     find the van der Waals energy via double loop search
c
!$acc parallel present(elrc,vlrc,radmin,epsilon,ivt,jvt) 
!$acc&         async(def_queue)
!$acc loop collapse(3) reduction(+:elrc,vlrc)
      do i = 1, nvt
         do k = 1, nvt
            do j = 1, ndelta
               it      = ivt(i)
               termi   = two * pi * real(jvt(i),t_p)
               kt      = ivt(k)
               termik  = termi * real(jvt(k),t_p)
               rv      = radmin(kt,it)
               eps     = epsilon(kt,it)
               rv2     = rv * rv
               rv6     = rv2 * rv2 * rv2
               rv7     = rv6 * rv
               r       = offset + real(j,t_p)*rdelta
               r2      = r * r
               r3      = r2 * r
               r6      = r3 * r3
               r7      = r6 * r
               e       = zero
               de      = zero
               if (vdwtyp_i .eq. 0) then
                  p6   = rv6 / r6
                  p12  = p6 * p6
                  e    = eps * (p12 - two*p6)
                  de   = eps * (p12-p6) * (-12.0_ti_p/r)
               else if (vdwtyp_i .eq. 1) then
                  rho  = r7 + ghal*rv7
                  tau  = (dhal+one) / (r+dhal*rv)
                  tau7 = tau**7
                  dtau = tau / (dhal+one)
                  gtau = eps*tau7*r6*(ghal+one)*(rv7/rho)**2
                  e    = eps * rv7 * tau7 * ((ghal+one)*rv7/rho-two)
                  de   = -7.0_ti_p * (dtau*e+gtau)
               end if
               if (r .lt. off) then
                  r4   = r2 * r2
                  r5   = r2 * r3
                  taper  = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0_ti_p*c5*r4 + 4.0_ti_p*c4*r3
     &                   + 3.0_ti_p*c3*r2 + two*c2*r + c1
                  de   = de*(one-taper) - e*dtaper
                  e    = e*(one-taper)
               end if 
               elrc    = elrc + termik*(e*rdelta*r2)
               vlrc    = vlrc + termik*(de*rdelta*r3)
            end do
         end do
      end do
      elrc = elrc / volbox
      vlrc = vlrc / (3.0_re_p*volbox)
!$acc end parallel
!$acc end data
      return
      end
