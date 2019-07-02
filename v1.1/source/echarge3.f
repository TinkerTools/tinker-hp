c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3  --  charge-charge energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms
c
c
      subroutine echarge3
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      call echarge3c
c
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3c  --  Ewald charge analysis via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3c" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation
c
c
      subroutine echarge3c
      use action
      use analyz
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use domdec
      use energi
      use ewald
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use usage
      use mpi
      implicit none
      integer i,j,k,iglob,iichg
      integer ii,kk,kkk,inl,kglob,kkchg
      integer in,kn
      real*8 e,efull
      real*8 f,fi,fik
      real*8 fs,fgrp
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      logical header,huge
      character*6 mode
      external erfc
c
c
c     zero out the Ewald summation energy and partitioning
c
      nec = 0
      ec = 0.0d0
      aec = 0.0d0
c
      if (nion .eq. 0)  return
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        call ecrecip
        if (use_pmecore) return
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      cscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
      do ii = 1, nionloc
         iichg = chgglob(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         e = fs * pchg(iichg)**2
         ec = ec + e
         nec = nec + 1
         aec(i) = aec(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nionloc
            iichg = chgglob(ii)
            iglob = iion(iichg)
            i = loc(iglob)
            xd = xd + pchg(iichg)*x(iglob)
            yd = yd + pchg(iichg)*y(iglob)
            zd = zd + pchg(iichg)*z(iglob)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
         nec = nec + 1
         do ii = 1, nionloc
            iichg = chgglob(ii)
            iglob = iion(iichg)
            i = loc(iglob)
            aec(i) = aec(i) + e/dble(nion)
         end do
      end if
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         inl = chglocnl(iichg)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         fi = f * pchg(iichg)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            cscale(i12(j,iglob)) = c2scale
         end do
         do j = 1, n13(iglob)
            cscale(i13(j,iglob)) = c3scale
         end do
         do j = 1, n14(iglob)
            cscale(i14(j,iglob)) = c4scale
         end do
         do j = 1, n15(iglob)
            cscale(i15(j,iglob)) = c5scale
         end do
         do kkk = 1, nelst(inl)
            kkchg = elst(kkk,inl)
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            k = loc(kglob)
c
c     compute the energy contribution for this interaction
c
            xr = xi - x(kglob)
            yr = yi - y(kglob)
            zr = zi - z(kglob)
c
c     find energy for interactions within real space cutoff
c
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rb = r + ebuffer
               fik = fi * pchg(kkchg)
               rew = aewald * r
               erfterm = erfc (rew)
               scale = cscale(kglob)
               scaleterm = scale - 1.0d0
               e = (fik/rb) * (erfterm+scaleterm)
               ec = ec + e
c
c     increment the overall charge-charge energy component
c
               efull = (fik/rb) * scale
               if (efull .ne. 0.0d0) then
                  nec = nec + 1
                  aec(i) = aec(i) + 0.5d0*efull
                  aec(k) = aec(k) + 0.5d0*efull
                  if (molcule(iglob) .ne. molcule(kglob))
     &               einter = einter + efull
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            cscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            cscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            cscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            cscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
