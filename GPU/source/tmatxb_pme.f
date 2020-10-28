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
      use sizes
      use atoms
      use atmlst
      use domdec
      use ewald
      use group
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use tinheader ,only: ti_p,re_p
      implicit none
      integer i,nrhs,iglob,kglob,kbis,iipole,kkpole,kkpoleloc
      integer ipoleloc
      real(t_p)  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      logical dodiag
      integer j, ii, kkk, irhs, inl
      real(t_p)  dx, dy, dz, d, d2, damp, expdamp, pgamma,
     $  scale3, scale5, pdi, pti 
      real(t_p) ralpha, alsq2, alsq2n, bfac, exp2a
      real(t_p) rr3, rr5, dukx, duky, dukz, pukx, puky, pukz,
     $  puir, pukr, duir, dukr
      real(t_p) duix, duiy, duiz, puix, puiy, puiz
      real(t_p) bn(0:3), fid(3), fip(3), fimd(3), fimp(3)
      real(t_p) fkd(3), fkp(3), fkmd(3), fkmp(3)
      real(t_p)  zero, one, f50
      real(t_p), allocatable :: dscale(:)
      real(t_p)  cutoff2
      save    zero, one, f50
      data    zero/0._ti_p/, one/1._ti_p/, f50/50._ti_p/
      character*10 mode
c
c     initialize the result vector
c
      do i = 1, npolebloc
        do irhs = 1, nrhs
          do j = 1, 3
            efi(j,irhs,i) = zero
          end do
        end do
      end do
      allocate (dscale(n))
      do i = 1, n
         dscale(i) = 1.0_ti_p
      end do
c
c     gather some parameters, then set up the damping factors.
c
      mode = 'EWALD'
      call switch (mode)
      if(rank.eq.0.and.tinkerdebug)write(*,*) 'tmatxb_pme'
      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        pdi = pdamp(iipole)
        pti = thole(iipole)
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = u4scale
        end do
c
        do kkk = 1,nelst(ii)
          kkpole = elst(kkk,ii)
          kglob = ipole(kkpole)
          kkpoleloc = poleloc(kkpole)
          if (kkpoleloc.eq.0) cycle
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.cutoff2) then
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
            d  = sqrt(d2)
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)
c
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0_ti_p * aewald**2
            alsq2n = 0.0_ti_p
            if (aewald .gt. 0.0_ti_p)
     &        alsq2n = 1.0_ti_p / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 2
              bfac = real(j+j-1,t_p)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / d2
            end do
c
            scale3 = dscale(kglob)
            scale5 = dscale(kglob)
            damp = pdi*pdamp(kkpole)
            if (damp.ne.zero) then
              pgamma = min(pti,thole(kkpole))
              damp = -pgamma*(d/damp)**3
              if (damp .gt. -f50) then
                expdamp = exp(damp)
                scale3 = scale3 * (one - expdamp)
                scale5 = scale5 * (one - expdamp*(one - damp))
              end if
            end if
c
c     compute the field.
c
            rr3 = (1.0_ti_p-scale3) / (d*d2)
            rr5 = 3.0_ti_p * (1.0_ti_p-scale5) / (d*d2*d2)
            duir = dx*duix + dy*duiy + dz*duiz
            dukr = dx*dukx + dy*duky + dz*dukz
            puir = dx*puix + dy*puiy + dz*puiz
            pukr = dx*pukx + dy*puky + dz*pukz
            fimd(1) = -bn(1)*dukx + bn(2)*dukr*dx
            fimd(2) = -bn(1)*duky + bn(2)*dukr*dy
            fimd(3) = -bn(1)*dukz + bn(2)*dukr*dz
            fkmd(1) = -bn(1)*duix + bn(2)*duir*dx
            fkmd(2) = -bn(1)*duiy + bn(2)*duir*dy
            fkmd(3) = -bn(1)*duiz + bn(2)*duir*dz
            fimp(1) = -bn(1)*pukx + bn(2)*pukr*dx
            fimp(2) = -bn(1)*puky + bn(2)*pukr*dy
            fimp(3) = -bn(1)*pukz + bn(2)*pukr*dz
            fkmp(1) = -bn(1)*puix + bn(2)*puir*dx
            fkmp(2) = -bn(1)*puiy + bn(2)*puir*dy
            fkmp(3) = -bn(1)*puiz + bn(2)*puir*dz
            fid(1) = -rr3*dukx + rr5*dukr*dx
            fid(2) = -rr3*duky + rr5*dukr*dy
            fid(3) = -rr3*dukz + rr5*dukr*dz
            fkd(1) = -rr3*duix + rr5*duir*dx
            fkd(2) = -rr3*duiy + rr5*duir*dy
            fkd(3) = -rr3*duiz + rr5*duir*dz
            fip(1) = -rr3*pukx + rr5*pukr*dx
            fip(2) = -rr3*puky + rr5*pukr*dy
            fip(3) = -rr3*pukz + rr5*pukr*dz
            fkp(1) = -rr3*puix + rr5*puir*dx
            fkp(2) = -rr3*puiy + rr5*puir*dy
            fkp(3) = -rr3*puiz + rr5*puir*dz
            efi(1,1,ipoleloc) = efi(1,1,ipoleloc) - fimd(1) + fid(1)
            efi(2,1,ipoleloc) = efi(2,1,ipoleloc) - fimd(2) + fid(2)
            efi(3,1,ipoleloc) = efi(3,1,ipoleloc) - fimd(3) + fid(3)
            efi(1,2,ipoleloc) = efi(1,2,ipoleloc) - fimp(1) + fip(1)
            efi(2,2,ipoleloc) = efi(2,2,ipoleloc) - fimp(2) + fip(2)
            efi(3,2,ipoleloc) = efi(3,2,ipoleloc) - fimp(3) + fip(3)
c
            efi(1,1,kkpoleloc) = efi(1,1,kkpoleloc) - fkmd(1) + fkd(1)
            efi(2,1,kkpoleloc) = efi(2,1,kkpoleloc) - fkmd(2) + fkd(2)
            efi(3,1,kkpoleloc) = efi(3,1,kkpoleloc) - fkmd(3) + fkd(3)
            efi(1,2,kkpoleloc) = efi(1,2,kkpoleloc) - fkmp(1) + fkp(1)
            efi(2,2,kkpoleloc) = efi(2,2,kkpoleloc) - fkmp(2) + fkp(2)
            efi(3,2,kkpoleloc) = efi(3,2,kkpoleloc) - fkmp(3) + fkp(3)
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = 1.0_ti_p
        end do
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
          if (polarity(iipole) .eq. 0_ti_p) then
             do irhs = 1, nrhs
               do j = 1, 3
                 efi(j,irhs,i) = efi(j,irhs,i) +
     $              mu(j,irhs,i)/tinypol
               end do
             end do
          else 
             do irhs = 1, nrhs
               do j = 1, 3
                 efi(j,irhs,i) = efi(j,irhs,i) +
     $              mu(j,irhs,i)/polarity(iipole)
               end do
             end do
          end if
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
