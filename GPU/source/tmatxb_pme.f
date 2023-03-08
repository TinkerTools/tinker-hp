c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_precision.h"
      subroutine tmatxb_pme(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use atmlst
      use atoms
      use bound
      use chgpen
      use couple
      use domdec
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      implicit none
      integer i,nrhs,iglob,kglob,iipole,kkpole,kkpoleloc,nnelst
      integer ipoleloc
      real(t_p)  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr2
      real(t_p) rr3,rr5
      real(t_p) rr3ik,rr5ik
      real(t_p) scalek
      real(t_p) dmp3,dmp5
      real(t_p) fid(3),fkd(3)
      real(t_p) fip(3),fkp(3)
      real(t_p) dmpik(7),dmpe(7)
      real(t_p) dlocal(6)
      real(t_p) duix,duiy,duiz
      real(t_p) puix,puiy,puiz
      real(t_p) dukx,duky,dukz
      real(t_p) pukx,puky,pukz
      real(t_p) alphai,alphak
      real(t_p) pol
      real(t_p), allocatable, target :: uscale(:)
      real(t_p), pointer     :: wscale(:)
      logical dodiag
      logical shortrange
      integer j, ii, kkk, irhs
      real(t_p) cutoff2
      character*11 mode
      character*80 :: RoutineName

      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='tmatxb_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='tmatxb_pme'
         mode = 'EWALD'
      endif
c
c     initialize the result vector
c
      efi = 0
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
      wscale => uscale
c
c     set arrays needed to scale connected atom interactions
c
      uscale(:) = 1.0
c
c     gather some parameters, then set up the damping factors.
c
      call switch (mode)

      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole   = poleglobnl(ii)
        iglob    = ipole     (iipole)
        i        = loc       (iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
        if (use_chgpen) then
           alphai = palpha(iipole)
        end if
c
c
c     set exclusion coefficients for connected atoms
c
        if (use_chgpen) then
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = w2scale
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = w3scale
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = w4scale
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = w5scale
        end do
        else
        do j = 1, np11(iglob)
           uscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           uscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           uscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           uscale(ip14(j,iglob)) = u4scale
        end do
        end if
        if (shortrange) then
          nnelst = nshortelst(ii)
        else
          nnelst = nelst(ii)
        end if
        do kkk = 1, nnelst
          if (shortrange) then
            kkpole = shortelst(kkk,ii)
          else
            kkpole = elst(kkk,ii)
          end if
          kglob = ipole(kkpole)
          kkpoleloc = poleloc(kkpole)
          if (kkpoleloc.eq.0) cycle
          xr = x(kglob) - x(iglob)
          yr = y(kglob) - y(iglob)
          zr = z(kglob) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2 .le. off2) then
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
            r = sqrt(r2)
            rr1 = 1.0 / r
            rr2 = rr1 * rr1
            rr3 = rr2 * rr1
            rr5 = 3.0 * rr2 * rr3
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)
            if (use_chgpen) then
               alphak = palpha(kkpole)
            end if
c
c     calculate real space Ewald error function damping
c
            call dampewald (7,r,r2,1.0,dmpe)
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               call dampthole2 (iipole,kkpole,5,r,dmpik)
               scalek = uscale(kglob)
               dmp3 = dmpe(3) - (1.0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0-scalek*dmpik(5))*rr5
               dlocal(1) = -dmp3 + dmp5*xr*xr
               dlocal(2) = dmp5*xr*yr
               dlocal(3) = dmp5*xr*zr
               dlocal(4) = -dmp3 + dmp5*yr*yr
               dlocal(5) = dmp5*yr*zr
               dlocal(6) = -dmp3 + dmp5*zr*zr
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               call dampmut (r,alphai,alphak,dmpik)
               scalek = wscale(kglob)
               rr3ik = dmpe(3) - (1.0-scalek*dmpik(3))*rr3
               rr5ik = dmpe(5) - (1.0-scalek*dmpik(5))*rr5
               dlocal(1) = -rr3ik + rr5ik*xr*xr
               dlocal(2) = rr5ik*xr*yr
               dlocal(3) = rr5ik*xr*zr
               dlocal(4) = -rr3ik + rr5ik*yr*yr
               dlocal(5) = rr5ik*yr*zr
               dlocal(6) = -rr3ik + rr5ik*zr*zr
            end if
            fid(1) = dlocal(1)*dukx+dlocal(2)*duky+dlocal(3)*dukz
            fid(2) = dlocal(2)*dukx+dlocal(4)*duky+dlocal(5)*dukz
            fid(3) = dlocal(3)*dukx+dlocal(5)*duky+dlocal(6)*dukz
            fkd(1) = dlocal(1)*duix+dlocal(2)*duiy+dlocal(3)*duiz
            fkd(2) = dlocal(2)*duix+dlocal(4)*duiy+dlocal(5)*duiz
            fkd(3) = dlocal(3)*duix+dlocal(5)*duiy+dlocal(6)*duiz

            fip(1) = dlocal(1)*pukx+dlocal(2)*puky+dlocal(3)*pukz
            fip(2) = dlocal(2)*pukx+dlocal(4)*puky+dlocal(5)*pukz
            fip(3) = dlocal(3)*pukx+dlocal(5)*puky+dlocal(6)*pukz
            fkp(1) = dlocal(1)*puix+dlocal(2)*puiy+dlocal(3)*puiz
            fkp(2) = dlocal(2)*puix+dlocal(4)*puiy+dlocal(5)*puiz
            fkp(3) = dlocal(3)*puix+dlocal(5)*puiy+dlocal(6)*puiz
            do j = 1, 3
               efi(j,1,ipoleloc) = efi(j,1,ipoleloc) - fid(j)
               efi(j,1,kkpoleloc) = efi(j,1,kkpoleloc) - fkd(j)
               efi(j,2,ipoleloc) = efi(j,2,ipoleloc) - fip(j)
               efi(j,2,kkpoleloc) = efi(j,2,kkpoleloc) - fkp(j)
            end do
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        if (use_chgpen) then
           do j = 1, n12(iglob)
              wscale(i12(j,iglob)) = 1.0
           end do
           do j = 1, n13(iglob)
              wscale(i13(j,iglob)) = 1.0
           end do
           do j = 1, n14(iglob)
              wscale(i14(j,iglob)) = 1.0
           end do
           do j = 1, n15(iglob)
              wscale(i15(j,iglob)) = 1.0
           end do
        else
           do j = 1, np11(iglob)
              uscale(ip11(j,iglob)) = 1.0
           end do
           do j = 1, np12(iglob)
              uscale(ip12(j,iglob)) = 1.0
           end do
           do j = 1, np13(iglob)
              uscale(ip13(j,iglob)) = 1.0
           end do
           do j = 1, np14(iglob)
              uscale(ip14(j,iglob)) = 1.0
           end do
        end if
      end do

      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npoleloc
          iipole = poleglob(i)
c
c     if no polarisability, take a negligeable value to allow convergence
c
          if (polarity(iipole).eq.0.0) then
             pol = tinypol ** (-1)
          else
             pol  = polarity(iipole) ** (-1)
          endif
          do irhs = 1, nrhs
            do j = 1, 3
              efi(j,irhs,i) = efi(j,irhs,i) +
     $           mu(j,irhs,i)*pol
            end do
          end do
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      return
      end
