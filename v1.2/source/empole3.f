c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions,
c     and partitions the energy among the atoms
c
      subroutine empole3
      use energi
      use potent
      use mpi
      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empole3c
c
      return
      end
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3c  --  Ewald multipole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3c" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among the atoms
c
c
      subroutine empole3c
      use sizes
      use action
      use analyz
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use mpi
      implicit none
      integer i,ii,iglob,iipole,ierr
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
      aem = 0d0
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $ then
        if (use_mrec) then
          call emrecip
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_mreal) then
          call emreal3d
        end if

        if (use_mself) then
c
c     compute the self-energy part of the Ewald summation
c
          term = 2.0d0 * aewald * aewald
          fterm = -f * aewald / sqrtpi
          do ii = 1, npoleloc
             iipole = poleglob(ii)
             iglob = ipole(iipole)
             i = loc(iglob)
             ci = rpole(1,iipole)
             dix = rpole(2,iipole)
             diy = rpole(3,iipole)
             diz = rpole(4,iipole)
             qixx = rpole(5,iipole)
             qixy = rpole(6,iipole)
             qixz = rpole(7,iipole)
             qiyy = rpole(9,iipole)
             qiyz = rpole(10,iipole)
             qizz = rpole(13,iipole)
             cii = ci*ci
             dii = dix*dix + diy*diy + diz*diz
             qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &                + qixx*qixx + qiyy*qiyy + qizz*qizz
             e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
             em = em + e
             nem = nem + 1
             aem(i) = aem(i) + e
          end do
c
c       compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                i = loc(iglob)
                dix = rpole(2,iipole)
                diy = rpole(3,iipole)
                diz = rpole(4,iipole)
                xd = xd + dix + rpole(1,iipole)*x(iglob)
                yd = yd + diy + rpole(1,iipole)*y(iglob)
                zd = zd + diz + rpole(1,iipole)*z(iglob)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               term = (2.0d0/3.0d0) * f * (pi/volbox)
               e = term * (xd*xd+yd*yd+zd*zd)
               em = em + e
               nem = nem + 1
             end if
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                i = loc(iglob)
                aem(i) = aem(i) + e/dble(npole)
             end do
          end if
        end if
      end if
      return
      end
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emreal3d  --  real space mpole analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emreal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions, and partitions
c     the energy among the atoms using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
c     if longrange, calculates just the long range part
c     if shortrange, calculates just the short range part

      subroutine emreal3d
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
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
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis,nnelst
      integer ii,kkk,iipole,kkpole
      real*8 e,efull,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 bn(0:4)
      real*8 fgrp,scale
      real*8 s,ds,mpoleshortcut2
      real*8 facts
      logical testcut,shortrange,longrange,fullrange
      real*8, allocatable :: mscale(:)
      logical header,huge
      character*10 mode
      character*80 :: RoutineName
      external erfc


c     compute the short, long, or full real space part of the summation
      shortrange = use_mpoleshortreal
      longrange  = use_mpolelong
      fullrange  = .not.(shortrange.or.longrange)
      if (shortrange) then 
         RoutineName = 'emrealshort3d'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'emreallong3d'
         mode        = 'EWALD'
      else
         RoutineName = 'emreal3d'
         mode        = 'EWALD'
      endif

      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      mscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
!     mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc  (iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         ci = rpole(1,iipole)
         dix = rpole(2,iipole)
         diy = rpole(3,iipole)
         diz = rpole(4,iipole)
         qixx = rpole(5,iipole)
         qixy = rpole(6,iipole)
         qixz = rpole(7,iipole)
         qiyy = rpole(9,iipole)
         qiyz = rpole(10,iipole)
         qizz = rpole(13,iipole)
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = m2scale
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = m3scale
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = m4scale
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         nnelst = merge(nshortelst(ii),
     &                  nelst     (ii),
     &                  shortrange
     &                 )
c        do kkk = 1, nelst(ii)
         do kkk = 1, nnelst
c           kkpole = elst(kkk,ii)
            kkpole = merge(shortelst(kkk,ii),
     &                     elst     (kkk,ii),
     &                     shortrange
     &                   )
            kglob = ipole(kkpole)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            kbis = loc(kglob)
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.mpoleshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
c           if (r2 .le. off2) then
            if (testcut) then
               r = sqrt(r2)
               ck = rpole(1,kkpole)
               dkx = rpole(2,kkpole)
               dky = rpole(3,kkpole)
               dkz = rpole(4,kkpole)
               qkxx = rpole(5,kkpole)
               qkxy = rpole(6,kkpole)
               qkxz = rpole(7,kkpole)
               qkyy = rpole(9,kkpole)
               qkyz = rpole(10,kkpole)
               qkzz = rpole(13,kkpole)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the full undamped energy for this interaction
c
               scale = mscale(kglob)
               if (use_group)  scale = scale * fgrp

               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = scale * efull
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  aem(i) = aem(i) + efull
                  if (molcule(iglob) .ne. molcule(kglob))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
c
               scalekk = 1.0d0 - scale
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               if(shortrange .or. longrange)
     &            call switch_respa(r,mpoleshortcut,shortheal,s,ds)
c
c     fix the s factor, depending on the range
c
               if(shortrange) then
                  facts  =         s
               else if(longrange) then
                  facts  = 1.0d0 - s
               else
                  facts  = 1.0d0
               endif

               em = em + facts * e
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  iglob,name(iglob),kglob,name(kglob),
     &             r,efull
   30             format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
