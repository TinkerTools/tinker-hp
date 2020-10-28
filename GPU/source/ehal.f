c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3  --  buffered 14-7 vdw energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal3" calculates the buffered 14-7 van der Waals energy
c     and partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      subroutine ehal
      use energi
      use potent
      use vdwpot
      implicit none
      real(t_p) elrc
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_vdwshort) then
        call ehalshort0c
      else if (use_vdwlong) then
        call ehallong0c
      else
        call ehal0c
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr (elrc)
         ev = ev + elrc
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ehal0c  --  buffered 14-7 analysis via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ehal0c" calculates the buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
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
      logical proceed,usei
      logical muti,mutk
      character*10 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      ev = 0.0_re_p
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1.0_ti_p
      iv14 = 0
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
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
c
c
c     find the van der Waals energy via neighbor list search
c
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
c
c     set interaction scaling coefficients for connected atoms
c
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
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(ii)
            kglob = vlst(kk,ii)
            kbis = loc(kglob)
            kv = ired(kglob)
            mutk = mut(kglob)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)
               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv  = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
c
c     get the interaction energy, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
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
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0_ti_p)*rv7/rho-2.0_ti_p)
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0_ti_p) then
                     ev = ev + e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
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
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c     ###############################################################################
c     ##                                                                           ##
c     ##  subroutine ehalshort0c  --  short range buffered 14-7 analysis via list  ##
c     ##                                                                           ##
c     ###############################################################################
c
c
c     "ehalshort0c" calculates the short range buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine ehalshort0c
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff
      use domdec
      use energi
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
      logical proceed,usei
      logical muti,mutk
      character*10 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      ev = 0.0_re_p
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1.0_ti_p
      iv14 = 0
c
c     set the coefficients for the switching function
c
      mode = 'SHORTVDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
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
c
c
c     find the van der Waals energy via neighbor list search
c
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
c
c     set interaction scaling coefficients for connected atoms
c
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
c
c     decide whether to compute the current interaction
c
         do kk = 1, nshortvlst(ii)
            kglob = shortvlst(kk,ii)
            kbis = loc(kglob)
            kv = ired(kglob)
            mutk = mut(kglob)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)
               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv  = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
c
c     get the interaction energy, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
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
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0_ti_p)*rv7/rho-2.0_ti_p)
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. (off-shortheal)*(off-shortheal)) then
                     call switch_respa(rik,off,shortheal,s,ds)
                     e = e * s
                  end if
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0_ti_p) then
                     ev = ev + e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
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
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c     ###############################################################################
c     ##                                                                           ##
c     ##  subroutine ehallong0c  --  long range buffered 14-7 analysis via list    ##
c     ##                                                                           ##
c     ###############################################################################
c
c
c     "ehallong0c" calculates the long range buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine ehallong0c
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff
      use domdec
      use energi
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
      logical proceed,usei
      logical muti,mutk
      character*10 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      ev = 0.0_re_p
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1.0_ti_p
      iv14 = 0
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
      vdwshortcut2 = (vdwshortcut-shortheal)**2
c
c     apply any reduction factor to the atomic coordinates
c
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
c
c
c     find the van der Waals energy via neighbor list search
c
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
c
c     set interaction scaling coefficients for connected atoms
c
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
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(ii)
            kglob = vlst(kk,ii)
            kbis = loc(kglob)
            kv = ired(kglob)
            mutk = mut(kglob)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)
               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if ((rik2 .le. off2).and.(rik2.ge.vdwshortcut2)) then
                  rik = sqrt(rik2)
                  rv  = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
c
c     get the interaction energy, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
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
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0_ti_p)*rv7/rho-2.0_ti_p)
                  end if
c
c     use energy switching if close the cutoff distance (at short range)
c
                  call switch_respa(rik,vdwshortcut,shortheal,s,ds)
                  e = (1-s)*e
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. (off-shortheal)*(off-shortheal)) then
                     call switch_respa(rik,off,shortheal,s,ds)
                     e = e * s
                  end if
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0_ti_p) then
                     ev = ev + e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
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
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
