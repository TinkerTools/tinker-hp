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
c     Liquids, 2nd Ed.", Oxford University Press, 2017, Section 2.8
c
c
#include "tinker_precision.h"
      subroutine evcorr (elrc)
      use atmtyp
      use atoms
      use bound
      use boxes
      use domdec   ,only:rank
      use math
      use mutant
      use potent
      use shunt
      use inform   ,only:deb_path
      use tinheader,only:ti_p
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,it,kt
      integer nstep,ndelta
      integer vdwtyp_i
      integer, allocatable :: mvt(:)
      real(t_p) zero,one,two
      real(r_p) elrc,etot
      real(t_p) range,rdelta
      real(t_p) fi,fk,fim,fkm,fik
      real(t_p) e,eps,vlam1
      real(t_p) offset,taper
      real(t_p) rv,rv2,rv6,rv7
      real(t_p) r,r2,r3,r4,rc
      real(t_p) r5,r6,r7
      real(t_p) cik,p,p6,p12
      real(t_p) t1,t2,ri
      real(t_p) rho,tau,tau7
      real(t_p) expterm
      character*10 mode
      enum,bind(C)
      enumerator LENNARD_JONES,BUFFERED_14_7
      end enum
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p)
c
c
c     zero out the long range van der Waals correction
c
!$acc serial async present(elrc)
      elrc = 0
!$acc end serial
c
c     only the master computes the vdw correction
c
      if (.not.use_bounds.or.rank.ne.0)  return
      if (deb_path) write(*,12) nvt
 12   format(3x,'evcorr',3x,'(nvt:',I8,')')
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
      if (vdwtyp.eq.'LENNARD-JONES') then
         vdwtyp_i=LENNARD_JONES
      else if (vdwtyp.eq.'BUFFERED-14-7') then
         vdwtyp_i=BUFFERED_14_7
      else
         return
      end if
c
c     set number of steps and range for numerical integration
c
      nstep = 2
      range = 100.0_ti_p
      ndelta = int(real(nstep,t_p)*(range-cut))
      rdelta = (range-cut) / real(ndelta,t_p)
      offset = cut - 0.5_ti_p*rdelta
      vlam1 = 1.0_ti_p - vlambda
      rc     = 1.0/(cut-off)
c
c     perform dynamic allocation of some local arrays
c
      allocate (mvt(n))
c
c     count the number of types and their frequencies
c
      nvt = 0
      do i = 1, n
         if (use_vdw)  it = jvdw(i)
         do k = 1, nvt
            if (ivt(k) .eq. it) then
               jvt(k) = jvt(k) + 1
               if (mut(i))  mvt(k) = mvt(k) + 1
               goto 10
            end if
         end do
         nvt = nvt + 1
         ivt(nvt) = it
         jvt(nvt) = 1
         mvt(nvt) = 0
         if (mut(i))  mvt(nvt) = 1
   10    continue
      end do

c
c     find the correction energy via double loop search
c
!$acc parallel loop collapse(2) async copyin(mvt)
!$acc&         present(ivt,jvt,radmin,epsilon,elrc)
      do i = 1, nvt
         do k = 1, nvt
            if (k.lt.i) cycle
            it  = ivt(i)
            kt  = ivt(k)
            fi  = 4.0_ti_p * pi * real(jvt(i),t_p)
            fim = 4.0_ti_p * pi * real(mvt(i),t_p)
            fk  = real(jvt(k),t_p)
            fkm = real(mvt(k),t_p)
c
c     set decoupling or annihilation for intraligand interactions
c
            if (vcouple .eq. 0) then
               fik = fi*fk - vlam1*(fim*(fk-fkm)+(fi-fim)*fkm)
            else
               fik = vlambda*fi*fk + vlam1*(fi-fim)*(fk-fkm)
            end if
            if (k .eq. i)  fik = 0.5_ti_p * fik
            rv  = radmin(kt,it)
            eps = epsilon(kt,it)
            etot = zero
c           rv2 = rv * rv
c           rv6 = rv2 * rv2 * rv2
c           rv7 = rv6 * rv
!$acc loop seq
            do j = 1, ndelta
               r = offset + real(j,t_p)*rdelta
c              r2 = r * r
c              r3 = r2 * r
c              r6 = r3 * r3
c              r7 = r6 * r
c              e = zero
               if (vdwtyp_i .eq. LENNARD_JONES) then
                  p6  = (rv / r)**6
                  e   = eps * (p6 - two)*p6
               else if (vdwtyp_i .eq. BUFFERED_14_7) then
                 rho  = r/rv
                 t1   = ((1.0_ti_p+dhal) / (rho+dhal))**7
                 t2   =  (1.0_ti_p+ghal) * ((rho**7+ghal)**(-1))
                 e    = eps * t1 * (t2-2.0_ti_p)
c               else if (vdwtyp.eq.'BUCKINGHAM' .or.
c     &                  vdwtyp.eq.'MM3-HBOND') then
c                  p = sqrt(rv2/r2)
c                  p6 = rv6 / r6
c                  expterm = abuck * exp(-bbuck/p)
c                  e = eps * (expterm - cbuck*p6)
               end if
               if (r .lt. off) then
                  ri    = (r-off)*rc
                  r2    = ri*ri
                  taper = r2*ri*(6*r2 - 15*ri + 10)
                  e     = e * (one - taper)
               end if
               etot = etot + e*rdelta*r*r
            end do
            elrc = elrc + fik*etot
         end do
      end do
!$acc serial async present(elrc)
      elrc = elrc / volbox
!$acc end serial
c
c     perform deallocation of some local arrays
c
      deallocate (mvt)
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
c     Liquids, 2nd Ed.", Oxford University Press, 2017, Section 2.8
c
c
      subroutine evcorr1 (elrc,vlrc)
      use atmtyp
      use atoms
      use bound
      use boxes
      use domdec   ,only:rank
      use inform   ,only:deb_path
      use math
      use mutant
      use potent
      use shunt
      use tinheader ,only:ti_p,re_p
      use utilgpu   ,only: def_queue
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,it,kt
      integer nstep,ndelta
      integer, allocatable :: mvt(:)
      real(t_p) zero,one,two
      real(r_p) elrc,vlrc
      real(t_p) etot,vtot
      real(t_p) range,rdelta
      real(t_p) fi,fk,fim,fkm,fik
      real(t_p) e,de,eps
      real(t_p) offset,vlam1
      real(t_p) taper,dtaper
      real(t_p) rv,rv2,rv6,rv7
      real(t_p) r,r2,r3,r4
      real(t_p) r5,r6,r7,ri
      real(t_p) cik,p,p6,p12
      real(t_p) rho,tau,tau7
      real(t_p) dtau,gtau
      real(t_p) rvterm,expterm
      character*10 mode
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p)
 12   format(2x,'evcorr1')
c
c
c     zero out the long range van der Waals corrections
c
      elrc = 0
      vlrc = 0
c
c     only the master computes the vdw correction
c
      if (.not. use_bounds.or.rank.ne.0)  return
      if (deb_path) write(*,12) 
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
      vlam1 = 1.0_ti_p - vlambda
c
c     perform dynamic allocation of some local arrays
c
      allocate (mvt(n))
c
c     count the number of vdw types and their frequencies
c
      nvt = 0
      do i = 1, n
         if (use_vdw)  it = jvdw(i)
         do k = 1, nvt
            if (ivt(k) .eq. it) then
               jvt(k) = jvt(k) + 1
               if (mut(i))  mvt(k) = mvt(k) + 1
               goto 10
            end if
         end do
         nvt = nvt + 1
         ivt(nvt) = it
         jvt(nvt) = 1
         mvt(nvt) = 0
         if (mut(i))  mvt(nvt) = 1
   10    continue
      end do
c
c     find the van der Waals energy via double loop search
c
      do i = 1, nvt
         do k = 1, nvt
            if (k.lt.i) cycle
         it = ivt(i)
         fi = 4.0_ti_p * pi * real(jvt(i),t_p)
         fim = 4.0_ti_p * pi * real(mvt(i),t_p)
            kt = ivt(k)
            fk = real(jvt(k),t_p)
            fkm = real(mvt(k),t_p)
c
c     set decoupling or annihilation for intraligand interactions
c
            if (vcouple .eq. 0) then
               fik = fi*fk - vlam1*(fim*(fk-fkm)+(fi-fim)*fkm)
            else
               fik = vlambda*fi*fk + vlam1*(fi-fim)*(fk-fkm)
            end if
            if (k .eq. i)  fik = 0.5_ti_p * fik
            rv = radmin(kt,it)
            eps = epsilon(kt,it)
            rv2 = rv * rv
            rv6 = rv2 * rv2 * rv2
            rv7 = rv6 * rv
            etot = 0.0_ti_p
            vtot = 0.0_ti_p
            do j = 1, ndelta
               r = offset + real(j,t_p)*rdelta
               r2 = r * r
               r3 = r2 * r
               r6 = r3 * r3
               r7 = r6 * r
               e = 0.0_ti_p
               de = 0.0_ti_p
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  p6 = rv6 / r6
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0_ti_p*p6)
                  de = eps * (p12-p6) * (-12.0_ti_p/r)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
                  rv7 = rv**7
                  rho = r7 + ghal*rv7
                  tau = (dhal+1.0_ti_p) / (r + dhal*rv)
                  tau7 = tau**7
                  dtau = tau / (dhal+1.0_ti_p)
                  gtau = eps*tau7*r6*(ghal+1.0_ti_p)*(rv7/rho)**2
                  e = eps*tau7*rv7*((ghal+1.0_ti_p)*rv7/rho-2.0_ti_p)
                  de = -7.0_ti_p * (dtau*e+gtau)
c               else if (vdwtyp.eq.'BUCKINGHAM' .or.
c     &                  vdwtyp.eq.'MM3-HBOND') then
c                  p = sqrt(rv2/r2)
c                  p6 = rv6 / r6
c                  rvterm = -bbuck / rv
c                  expterm = abuck * exp(-bbuck/p)
c                  e = eps * (expterm - cbuck*p6)
c                  de = eps * (rvterm*expterm+6.0d0*cbuck*p6/r)
               end if
               if (r .lt. off) then
                  r4 = r2 * r2
                  r5 = r2 * r3
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0_ti_p*c5*r4 + 4.0_ti_p*c4*r3
     &                        + 3.0_ti_p*c3*r2 + 2.0_ti_p*c2*r + c1
                  de = de*(1.0_ti_p-taper) - e*dtaper
                  e = e*(1.0_ti_p-taper)
               end if
               print*,e,taper,r2
               etot = etot + e*rdelta*r2
               vtot = vtot + de*rdelta*r3
            end do
            elrc = elrc + fik*etot
            vlrc = vlrc + fik*vtot
         end do
      end do
      elrc = elrc / volbox
      vlrc = vlrc / (3.0_re_p*volbox)
c
c     perform deallocation of some local arrays
c
      deallocate (mvt)
      return
      end
c
c     Device version of evcorr1
c
      subroutine evcorr1gpu (elrc,vlrc)
      use atmtyp
      use atoms
      use bound
      use boxes
      use domdec    ,only: rank
      use inform    ,only: deb_path
      use math
      use mutant
      use potent
      use shunt
      use tinMemory ,only: prmem_request
      use tinheader ,only: ti_p,re_p
      use utilgpu   ,only: def_queue
      use vdw
      use vdwpot
      implicit none
      real(r_p),intent(inout):: elrc,vlrc

      integer i,j,k,it,kt
      integer nstep,ndelta
      integer vdwtyp_i
      integer, save, allocatable :: mvt(:)
      logical, save :: f_in=.true.
      real(t_p) zero,one,two
      real(r_p) etot,vtot
      real(t_p) range,rdelta
      real(t_p) fi,fk,fim,fkm,fik
      real(t_p) e,de,eps
      real(t_p) offset,vlam1
      real(t_p) t1,t2,s1,s2
      real(t_p) dt1drho,dt2drho
      real(t_p) taper,dtaper
      real(t_p) rv,rv2,rv6,rv7
      real(t_p) r,r2,r3,r4,ri,ri2,ri3
      real(t_p) r5,r6,r7,rinv
      real(t_p) cik,p,p6,p12
      real(t_p) rho,rho6,tau,tau7
      real(t_p) dtau,gtau
      real(t_p) rvterm,expterm
      character*10 mode
      enum,bind(C)
      enumerator LENNARD_JONES,BUFFERED_14_7
      end enum
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p)
c
c
c     zero out the long range van der Waals corrections
c
!$acc serial async(def_queue) present(elrc,vlrc)
      elrc = 0.0_re_p
      vlrc = 0 0_re_p
!$acc end serial
c
c     only the master computes the vdw correction
c
      if (.not.use_bounds.or.rank.ne.0)  return
      if (deb_path) write(*,12) nvt
 12   format(3x,'evcorr1gpu',3x,'(nvt:',I8,')')
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
      if (vdwtyp.eq.'LENNARD-JONES') then
         vdwtyp_i=LENNARD_JONES
      else if (vdwtyp.eq.'BUFFERED-14-7') then
         vdwtyp_i=BUFFERED_14_7
      else
         return
      end if
c
c     set number of steps and range for numerical integration
c
      nstep = 2
      range = 100.0_ti_p
      ndelta = int(real(nstep,t_p)*(range-cut))
      rdelta = (range-cut) / real(ndelta,t_p)
      offset = cut - 0.5_ti_p*rdelta
      rinv   = 1.0_ti_p/(cut-off)
      vlam1  = 1.0_ti_p - vlambda

      if (f_in) then
      ! Request memory for mvt
      call prmem_request(mvt,n)
c
c     count the number of vdw types and their frequencies
c
      nvt = 0
      do i = 1, n
         if (use_vdw)  it = jvdw(i)
         do k = 1, nvt
            if (ivt(k) .eq. it) then
               jvt(k) = jvt(k) + 1
               if (mut(i))  mvt(k) = mvt(k) + 1
               goto 10
            end if
         end do
         nvt = nvt + 1
         ivt(nvt) = it
         jvt(nvt) = 1
         mvt(nvt) = 0
         if (mut(i))  mvt(nvt) = 1
   10    continue
      end do
!$acc update device(mvt)

      f_in=.false.
      end if
c
c     find the van der Waals energy via double loop search
c
!$acc parallel loop collapse(3) async(def_queue)
!$acc&         reduction(+:elrc,vlrc)
!$acc&         present(elrc,vlrc,radmin,epsilon,ivt,jvt,mvt)
      do i = 1, nvt
         do k = 1, nvt
            do j = 1, ndelta
               if (k.lt.i) cycle
               it  = ivt(i)
               kt  = ivt(k)
               fi  = 4.0_ti_p * pi * real(jvt(i),t_p)
               fim = 4.0_ti_p * pi * real(mvt(i),t_p)
               fk  = real(jvt(k),t_p)
               fkm = real(mvt(k),t_p)
c
c     set decoupling or annihilation for intraligand interactions
c
               if (vcouple .eq. 0) then
                  fik = fi*fk - vlam1*(fim*(fk-fkm)+(fi-fim)*fkm)
               else
                  fik = vlambda*fi*fk + vlam1*(fi-fim)*(fk-fkm)
               end if
               if (k .eq. i)  fik = 0.5_ti_p * fik
               r   = offset + real(j,t_p)*rdelta
               rv  = radmin(kt,it)
               eps = epsilon(kt,it)
               etot = 0.0_re_p
               vtot = 0.0_re_p
c              rv2 = rv * rv
c              rv6 = rv2 * rv2 * rv2
c              rv7 = rv6 * rv
c              r2 = r * r
c              r3 = r2 * r
c              r6 = r3 * r3
c              r7 = r6 * r
               if (vdwtyp_i .eq. LENNARD_JONES) then
                  p6  = (rv / r)**6
                  e    = eps * (p6 - two)*p6
                  de   = eps * (p6 - one)*p6 * (-12.0_ti_p/r)
               else if (vdwtyp_i .eq. BUFFERED_14_7) then
                  rho  = r / rv
                  rho6 = rho**6
                  s1   = ((rho + dhal)**7)**(-1)
                  s2   = (rho6*rho + ghal)**(-1)
                  t1   = (1.0_ti_p + dhal)**7 *s1
                  t2   = (1.0_ti_p + ghal)*s2
                  dt1drho = -7.0_ti_p*(rho+dhal)**6 *t1*s1
                  dt2drho = -7.0_ti_p*rho6*t2*s2
                  e    = eps*t1*(t2 - 2.0_ti_p)
                  de   = eps*(dt1drho*(t2-2.0_ti_p)+t1*dt2drho) / rv
c              else if (vdwtyp.eq.'BUCKINGHAM' .or.
c     &                 vdwtyp.eq.'MM3-HBOND') then
c                 p = sqrt(rv2/r2)
c                 p6 = rv6 / r6
c                 rvterm = -bbuck / rv
c                 expterm = abuck * exp(-bbuck/p)
c                 e = eps * (expterm - cbuck*p6)
c                 de = eps * (rvterm*expterm+6.0d0*cbuck*p6/r)
               end if
               if (r .lt. off) then
                  ri     = (r - off) * rinv
                  ri2    = ri * ri
                  ri3    = ri2 * ri
                  taper  = ri3 * (6*ri2 - 15*ri + 10)
                  dtaper = 30* (ri*(1.0-ri))*(ri*(1.0-ri)) *rinv;

                  de     = de*(one-taper) - e*dtaper
                  e      =  e*(one-taper)
               end if
               etot = (e*rdelta*r*r)
               vtot = (de*rdelta*r*r*r)
               elrc = elrc + fik*(etot)
               vlrc = vlrc + fik*(vtot)
            end do
         end do
      end do

!$acc serial async(def_queue) present(elrc,vlrc)
      elrc = elrc / volbox
      vlrc = vlrc / (3.0_re_p*volbox)
!$acc end serial

      end
