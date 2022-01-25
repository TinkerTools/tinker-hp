c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge  --  charge-charge energy & analysis   ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms
c
c
      subroutine echarge
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      call echarge0c
c
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge0c  --  Ewald charge analysis via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge0c" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation
c
c
      subroutine echarge0c
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
      integer i,iglob,iichg
      integer ii
      real*8 e
      real*8 f
      real*8 fs
      real*8 xd,yd,zd
      external erfc
c
c
c     zero out the Ewald summation energy and partitioning
c
      ec = 0.0d0
c
      if (nion .eq. 0)  return
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        if (use_crec) then
          call ecrecip
        end if
        if (use_pmecore) return
      end if

      if (use_cself) then
c
c     compute the Ewald self-energy term over all the atoms
c
        f = electric / dielec
        fs = -f * aewald / sqrtpi
        do ii = 1, nionloc
           iichg = chgglob(ii)
           iglob = iion(iichg)
           i = loc(iglob)
           e = fs * pchg(iichg)**2
           ec = ec + e
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
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_creal) then
          call ecreal0d
        end if
      end if
      return
      end
c
c     "ecreal0d" evaluates the real space portion of the Ewald sum
c     energy due to atomic charge interactions, 
c     using a pairwise neighbor list
c
c     if longrange, calculates just the long range part
c     if shortrange, calculates just the short range part
c
      subroutine ecreal0d
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use cutoff
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
      integer i,j,k,iglob,iichg,nnelst
      integer ii,kkk,kglob,kkchg
      real*8 e,efull
      real*8 f,fi,fik
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 fgrp
      real*8, allocatable :: cscale(:)
      real*8 s,ds,cshortcut2,facts
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName
      external erfc

c     compute the short, long, or full real space part of the Ewald summation
      shortrange = use_cshortreal
      longrange  = use_clong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'ecrealshort3d'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'ecreallong3d'
         mode        = 'EWALD'
      else
         RoutineName = 'ecreal3d'
         mode        = 'EWALD'
      endif

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
      call switch (mode)
      cshortcut2 = (chgshortcut - shortheal) ** 2

c
c     compute the real space portion of the Ewald summation
c
      MAINLOOP:
     &do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i = loc(iglob)
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
         nnelst = merge(nshortelst(ii),
     &                  nelst     (ii),
     &                  shortrange
     &                 )
         do kkk = 1, nnelst
            kkchg = merge(shortelst(kkk,ii),
     &                    elst     (kkk,ii),
     &                    shortrange
     &                   )
            if (kkchg.eq.0) cycle
            kglob = iion(kkchg)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
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
            testcut = merge(r2 .le. off2.and.r2.ge.cshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
c           if (r2 .le. off2) then
            if (testcut) then
               r = sqrt(r2)
               rb = r + ebuffer
               fik = fi * pchg(kkchg)
               rew = aewald * r
               erfterm = erfc (rew)
               scale = cscale(kglob)
               if (use_group)  scale = scale * fgrp
               scaleterm = scale - 1.0d0
               e = (fik/rb) * (erfterm+scaleterm)
c
c     use energy switching if near the cutoff distance
c     at short or long range
c
               if(shortrange .or. longrange)
     &            call switch_respa(r,chgshortcut,shortheal,s,ds)

               if(shortrange) then
                  facts =         s
               else if(longrange) then
                  facts = 1.0d0 - s
               endif
               ec = ec + e * facts
c
c     increment the overall charge-charge energy component
c
               efull = (fik/rb) * scale
               if (efull .ne. 0.0d0) then
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
      end do MAINLOOP
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
