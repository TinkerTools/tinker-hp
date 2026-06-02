!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine ehal3  --  buffered 14-7 vdw energy & analysis  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "ehal3" calculates the buffered 14-7 van der Waals energy
!     and partitions the energy among the atoms
!
!
#include "tinker_macro.h"
subroutine ehal
   use energi
   use potent
   use vdwpot
   implicit none
   real(t_p) elrc
   character*11 mode
   logical      fullrange
!
   fullrange = .not.(use_vdwshort.or.use_vdwlong)
!
!     choose the method for summing over pairwise interactions
!
   if (use_vdwshort) then
      call ehalshort0c
   else if (use_vdwlong) then
      call ehallong0c
   else
      call ehal0c
   end if
!
!     apply long range van der Waals correction if desired
!
   if (use_vcorr.and.(fullrange.or.use_vdwlong)) then
      mode = 'VDW'
      call evcorr (mode,elrc)
      ev = ev + elrc
   end if
   return
end
!
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine ehal0c  --  buffered 14-7 analysis via list  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "ehal0c" calculates the buffered 14-7 van der Waals energy
!     and also partitions the energy among the atoms using a
!     pairwise neighbor list
!
!
subroutine ehal0c
   use atmlst
   use atoms
   use bound
   use couple
   use domdec
   use energi
   use group
   use inter
   use molcul
   use mutant
   use neigh
   use shunt
   use tinheader ,only:ti_p,re_p
   use usage
   use vdw
   use vdwpot
   implicit none
   integer i,j,iglob,kglob,kbis
   integer ii,iv,it,iivdw,inl
   integer kk,kv,kt
   integer, allocatable :: iv14(:)
   real(t_p) e,eps,rdn
   real(t_p) rv,rv7
   real(t_p) xi,yi,zi
   real(t_p) xr,yr,zr
   real(t_p) rho,tau,taper
   real(t_p) scal,t1,t2
   real(t_p) rik,rik2,rik3
   real(t_p) rik4,rik5,rik7
   real(t_p), allocatable :: xred(:)
   real(t_p), allocatable :: yred(:)
   real(t_p), allocatable :: zred(:)
   real(t_p), allocatable :: vscale(:)
   real(t_p) :: fgrp
   logical proceed,usei
   logical muti,mutk
   character*10 mode
!
!
!     zero out the van der Waals energy and partitioning terms
!
   ev = 0.0_re_p
!
!     perform dynamic allocation of some local arrays
!
   allocate (iv14(n))
   allocate (xred(nbloc))
   allocate (yred(nbloc))
   allocate (zred(nbloc))
   allocate (vscale(n))
!
!     set arrays needed to scale connected atom interactions
!
   vscale = 1.0_ti_p
   iv14 = 0
!
!     set the coefficients for the switching function
!
   mode = 'VDW'
   call switch (mode)
!
!     apply any reduction factor to the atomic coordinates
!
   do ii = 1, nvdwbloc
      iivdw = vdwglob(ii)
      iglob = ivdw(iivdw)
      i = loc(iglob)
      iv = ired(iglob)
      rdn = kred(iglob)
      xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
      yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
      zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
   end do
!
!
!     find the van der Waals energy via neighbor list search
!
   do ii = 1, nvdwlocnl
      iivdw = vdwglobnl(ii)
      iglob = ivdw(iivdw)
      i = loc(iglob)
      iv = ired(iglob)
      it = jvdw(iglob)
      xi = xred(i)
      yi = yred(i)
      zi = zred(i)
      usei = (use(iglob) .or. use(iv))
      muti = mut(iglob)
!
!     set interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = v2scale
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = v3scale
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = v4scale
         iv14(i14(j,iglob)) = iglob
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = v5scale
      end do
!
!     decide whether to compute the current interaction
!
      do kk = 1, nvlst(ii)
         kglob = vlst(kk,ii)
         kbis = loc(kglob)
         kv = ired(kglob)
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob) .or. use(kv))
         if(use_group) call groups(fgrp,iglob,kglob,0,0,0,0)
!
!     compute the energy contribution for this interaction
!
         if (proceed) then
            kt = jvdw(kglob)
            xr = xi - xred(kbis)
            yr = yi - yred(kbis)
            zr = zi - zred(kbis)
            if (use_bounds) call image (xr,yr,zr)
            rik2 = xr*xr + yr*yr + zr*zr
!
!     check for an interaction distance less than the cutoff
!
            if (rik2 .le. off2) then
               rik = sqrt(rik2)
               rv  = radmin(kt,it)
               eps = epsilon(kt,it)
               if (iv14(kglob) .eq. iglob) then
                  rv = radmin4(kt,it)
                  eps = epsilon4(kt,it)
               end if
               eps = eps * vscale(kglob)
!
!     get the interaction energy, via soft core if necessary
!
               if ((muti .and. .not.mutk) .or.&
                  &(mutk .and. .not.muti)) then
                  rho = rik / rv
                  eps = eps * vlambda**scexp
                  scal = scalpha * (1.0_ti_p-vlambda)**2
                  t1 = (1.0_ti_p+dhal)**7 / (scal+(rho+dhal)**7)
                  t2 = (1.0_ti_p+ghal) / (scal+rho**7+ghal)
                  e = eps * t1 * (t2-2.0_ti_p)
               else
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0_ti_p) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7&
                     &* ((ghal+1.0_ti_p)*rv7/rho-2.0_ti_p)
               end if
!
!     use energy switching if near the cutoff distance
!
               if (rik2 .gt. cut2) then
                  rik3 = rik2 * rik
                  rik4 = rik2 * rik2
                  rik5 = rik2 * rik3
                  taper = c5*rik5 + c4*rik4 + c3*rik3&
                     &+ c2*rik2 + c1*rik + c0
                  e = e * taper
               end if
!
!     scale the interaction by the group factor
!
               if (use_group) e = e * fgrp
!
!     increment the overall van der Waals energy components
!
               if (e .ne. 0.0_ti_p) then
                  ev = ev + e
               end if
!
!     increment the total intermolecular energy
!
               if (molcule(iglob) .ne. molcule(kglob)) then
                  einter = einter + e
               end if
!
            end if
         end if
      end do
!
!     reset interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = 1.0_ti_p
      end do
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (iv14)
   deallocate (xred)
   deallocate (yred)
   deallocate (zred)
   deallocate (vscale)
   return
end
!
!     ###############################################################################
!     ##                                                                           ##
!     ##  subroutine ehalshort0c  --  short range buffered 14-7 analysis via list  ##
!     ##                                                                           ##
!     ###############################################################################
!
!
!     "ehalshort0c" calculates the short range buffered 14-7 van der Waals energy
!     and also partitions the energy among the atoms using a
!     pairwise neighbor list
!
!
subroutine ehalshort0c
   use atmlst
   use atoms
   use bound
   use couple
   use cutoff
   use domdec
   use energi
   use group
   use inter
   use molcul
   use mutant
   use neigh
   use shunt
   use tinheader ,only:ti_p,re_p
   use usage
   use vdw
   use vdwpot
   implicit none
   integer i,j,iglob,kglob,kbis
   integer ii,iv,it,iivdw
   integer kk,kv,kt
   integer, allocatable :: iv14(:)
   real(t_p) e,eps,rdn
   real(t_p) rv,rv7
   real(t_p) xi,yi,zi
   real(t_p) xr,yr,zr
   real(t_p) rho,tau
   real(t_p) scal,t1,t2
   real(t_p) rik,rik2
   real(t_p) rik7
   real(t_p) s,ds
   real(t_p), allocatable :: xred(:)
   real(t_p), allocatable :: yred(:)
   real(t_p), allocatable :: zred(:)
   real(t_p), allocatable :: vscale(:)
   real(t_p) :: fgrp
   logical proceed,usei
   logical muti,mutk
   character*10 mode
!
!
!     zero out the van der Waals energy and partitioning terms
!
   ev = 0.0_re_p
!
!     perform dynamic allocation of some local arrays
!
   allocate (iv14(n))
   allocate (xred(nbloc))
   allocate (yred(nbloc))
   allocate (zred(nbloc))
   allocate (vscale(n))
!
!     set arrays needed to scale connected atom interactions
!
   vscale = 1.0_ti_p
   iv14 = 0
!
!     set the coefficients for the switching function
!
   mode = 'SHORTVDW'
   call switch (mode)
!
!     apply any reduction factor to the atomic coordinates
!
   do ii = 1, nvdwbloc
      iivdw = vdwglob(ii)
      iglob = ivdw(iivdw)
      i = loc(iglob)
      iv = ired(iglob)
      rdn = kred(iglob)
      xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
      yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
      zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
   end do
!
!
!     find the van der Waals energy via neighbor list search
!
   do ii = 1, nvdwlocnl
      iivdw = vdwglobnl(ii)
      iglob = ivdw(iivdw)
      i = loc(iglob)
      iv = ired(iglob)
      it = jvdw(iglob)
      xi = xred(i)
      yi = yred(i)
      zi = zred(i)
      usei = (use(iglob) .or. use(iv))
      muti = mut(iglob)
!
!     set interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = v2scale
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = v3scale
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = v4scale
         iv14(i14(j,iglob)) = iglob
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = v5scale
      end do
!
!     decide whether to compute the current interaction
!
      do kk = 1, nshortvlst(ii)
         kglob = shortvlst(kk,ii)
         kbis = loc(kglob)
         kv = ired(kglob)
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob) .or. use(kv))
         if(use_group) call groups(fgrp,iglob,kglob,0,0,0,0)
!
!     compute the energy contribution for this interaction
!
         if (proceed) then
            kt = jvdw(kglob)
            xr = xi - xred(kbis)
            yr = yi - yred(kbis)
            zr = zi - zred(kbis)
            if (use_bounds) call image (xr,yr,zr)
            rik2 = xr*xr + yr*yr + zr*zr
!
!     check for an interaction distance less than the cutoff
!
            if (rik2 .le. off2) then
               rik = sqrt(rik2)
               rv  = radmin(kt,it)
               eps = epsilon(kt,it)
               if (iv14(kglob) .eq. iglob) then
                  rv = radmin4(kt,it)
                  eps = epsilon4(kt,it)
               end if
               eps = eps * vscale(kglob)
!
!     get the interaction energy, via soft core if necessary
!
               if ((muti .and. .not.mutk) .or.&
                  &(mutk .and. .not.muti)) then
                  rho = rik / rv
                  eps = eps * vlambda**scexp
                  scal = scalpha * (1.0_ti_p-vlambda)**2
                  t1 = (1.0_ti_p+dhal)**7 / (scal+(rho+dhal)**7)
                  t2 = (1.0_ti_p+ghal) / (scal+rho**7+ghal)
                  e = eps * t1 * (t2-2.0_ti_p)
               else
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0_ti_p) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7&
                     &* ((ghal+1.0_ti_p)*rv7/rho-2.0_ti_p)
               end if
!
!     scale the interaction by the group factor
!
               if(use_group) e=e*fgrp
!
!     use energy switching if near the cutoff distance
!
               if (rik2 .gt. (off-shortheal)*(off-shortheal)) then
                  call switch_respa(rik,off,shortheal,s,ds)
                  e = e * s
               end if
!
!     increment the overall van der Waals energy components
!
               if (e .ne. 0.0_ti_p) then
                  ev = ev + e
               end if
!
!     increment the total intermolecular energy
!
               if (molcule(iglob) .ne. molcule(kglob)) then
                  einter = einter + e
               end if
!
            end if
         end if
      end do
!
!     reset interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = 1.0_ti_p
      end do
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (iv14)
   deallocate (xred)
   deallocate (yred)
   deallocate (zred)
   deallocate (vscale)
   return
end
!
!     ###############################################################################
!     ##                                                                           ##
!     ##  subroutine ehallong0c  --  long range buffered 14-7 analysis via list    ##
!     ##                                                                           ##
!     ###############################################################################
!
!
!     "ehallong0c" calculates the long range buffered 14-7 van der Waals energy
!     and also partitions the energy among the atoms using a
!     pairwise neighbor list
!
!
subroutine ehallong0c
   use atmlst
   use atoms
   use bound
   use couple
   use cutoff
   use domdec
   use energi
   use group
   use inter
   use molcul
   use mutant
   use neigh
   use shunt
   use tinheader ,only:ti_p,re_p
   use usage
   use vdw
   use vdwpot
   implicit none
   integer i,j,iglob,kglob,kbis
   integer ii,iv,it,iivdw
   integer kk,kv,kt
   integer, allocatable :: iv14(:)
   real(t_p) e,eps,rdn
   real(t_p) rv,rv7
   real(t_p) xi,yi,zi
   real(t_p) xr,yr,zr
   real(t_p) rho,tau
   real(t_p) scal,t1,t2
   real(t_p) rik,rik2
   real(t_p) rik7
   real(t_p) s,ds,vdwshortcut2
   real(t_p), allocatable :: xred(:)
   real(t_p), allocatable :: yred(:)
   real(t_p), allocatable :: zred(:)
   real(t_p), allocatable :: vscale(:)
   real(t_p) fgrp
   logical proceed,usei
   logical muti,mutk
   character*10 mode
!
!
!     zero out the van der Waals energy and partitioning terms
!
   ev = 0.0_re_p
!
!     perform dynamic allocation of some local arrays
!
   allocate (iv14(n))
   allocate (xred(nbloc))
   allocate (yred(nbloc))
   allocate (zred(nbloc))
   allocate (vscale(n))
!
!     set arrays needed to scale connected atom interactions
!
   vscale = 1.0_ti_p
   iv14 = 0
!
!     set the coefficients for the switching function
!
   mode = 'VDW'
   call switch (mode)
   vdwshortcut2 = (vdwshortcut-shortheal)**2
!
!     apply any reduction factor to the atomic coordinates
!
   do ii = 1, nvdwbloc
      iivdw = vdwglob(ii)
      iglob = ivdw(iivdw)
      i = loc(iglob)
      iv = ired(iglob)
      rdn = kred(iglob)
      xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
      yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
      zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
   end do
!
!
!     find the van der Waals energy via neighbor list search
!
   do ii = 1, nvdwlocnl
      iivdw = vdwglobnl(ii)
      iglob = ivdw(iivdw)
      i = loc(iglob)
      iv = ired(iglob)
      it = jvdw(iglob)
      xi = xred(i)
      yi = yred(i)
      zi = zred(i)
      usei = (use(iglob) .or. use(iv))
      muti = mut(iglob)
!
!     set interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = v2scale
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = v3scale
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = v4scale
         iv14(i14(j,iglob)) = iglob
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = v5scale
      end do
!
!     decide whether to compute the current interaction
!
      do kk = 1, nvlst(ii)
         kglob = vlst(kk,ii)
         kbis = loc(kglob)
         kv = ired(kglob)
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob) .or. use(kv))
         if(use_group) call groups(fgrp,iglob,kglob,0,0,0,0)
!
!     compute the energy contribution for this interaction
!
         if (proceed) then
            kt = jvdw(kglob)
            xr = xi - xred(kbis)
            yr = yi - yred(kbis)
            zr = zi - zred(kbis)
            if (use_bounds) call image (xr,yr,zr)
            rik2 = xr*xr + yr*yr + zr*zr
!
!     check for an interaction distance less than the cutoff
!
            if ((rik2 .le. off2).and.(rik2.ge.vdwshortcut2)) then
               rik = sqrt(rik2)
               rv  = radmin(kt,it)
               eps = epsilon(kt,it)
               if (iv14(kglob) .eq. iglob) then
                  rv = radmin4(kt,it)
                  eps = epsilon4(kt,it)
               end if
               eps = eps * vscale(kglob)
!
!     get the interaction energy, via soft core if necessary
!
               if ((muti .and. .not.mutk) .or.&
                  &(mutk .and. .not.muti)) then
                  rho = rik / rv
                  eps = eps * vlambda**scexp
                  scal = scalpha * (1.0_ti_p-vlambda)**2
                  t1 = (1.0_ti_p+dhal)**7 / (scal+(rho+dhal)**7)
                  t2 = (1.0_ti_p+ghal) / (scal+rho**7+ghal)
                  e = eps * t1 * (t2-2.0_ti_p)
               else
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0_ti_p) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7&
                     &* ((ghal+1.0_ti_p)*rv7/rho-2.0_ti_p)
               end if
!
!     scale the interaction based on its group membership
!
               if (use_group) then
                  e = e * fgrp
               end if
!
!     use energy switching if close the cutoff distance (at short range)
!
               call switch_respa(rik,vdwshortcut,shortheal,s,ds)
               e = (1-s)*e
!
!     use energy switching if near the cutoff distance
!
               if (rik2 .gt. (off-shortheal)*(off-shortheal)) then
                  call switch_respa(rik,off,shortheal,s,ds)
                  e = e * s
               end if
!
!     increment the overall van der Waals energy components
!
               if (e .ne. 0.0_ti_p) then
                  ev = ev + e
               end if
!
!     increment the total intermolecular energy
!
               if (molcule(iglob) .ne. molcule(kglob)) then
                  einter = einter + e
               end if
!
            end if
         end if
      end do
!
!     reset interaction scaling coefficients for connected atoms
!
      do j = 1, n12(iglob)
         vscale(i12(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n13(iglob)
         vscale(i13(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n14(iglob)
         vscale(i14(j,iglob)) = 1.0_ti_p
      end do
      do j = 1, n15(iglob)
         vscale(i15(j,iglob)) = 1.0_ti_p
      end do
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (iv14)
   deallocate (xred)
   deallocate (yred)
   deallocate (zred)
   deallocate (vscale)
   return
end
