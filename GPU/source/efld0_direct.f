c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_precision.h"
      subroutine efld0_direct(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use couple
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use tinheader
      implicit none
      integer i,iglob,kglob,nrhs,ipoleloc
      real(t_p)  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole,kbis
      real(t_p)  dx, dy, dz, d, d2, d3, d5, pdi, pti, ck,
     $  dkx, dky, dkz, dkr, qkxx, qkyy, qkzz, qkxy, qkxz, qkyz,
     $  qkx, qky, qkz, qkr
      real(t_p) dix, diy, diz, dir
      real(t_p) qixx, qiyy, qizz, qixy, qixz, qiyz,
     $  qix, qiy, qiz, qir, ci
      real(t_p) damp, pgamma, expdamp, scale3, scale5, scale7
      real(t_p) drr3,drr5,drr7
      real(t_p) prr3,prr5,prr7
      real(t_p) dsc3,dsc5,dsc7
      real(t_p) psc3,psc5,psc7
      real(t_p) zero, pt6, one, f50
      real(t_p) ralpha, alsq2, alsq2n, bfac, exp2a
      real(t_p) bn(0:3), fim(3), fid(3), fip(3)
      real(t_p) fkd(3), fkp(3), fkm(3)
      real(t_p) erfc, cutoff2
      real(t_p), allocatable :: dscale(:)
      real(t_p), allocatable :: pscale(:)
      save   zero, pt6, one, f50
      data   zero/0._ti_p/, pt6/0.6_ti_p/, one/1._ti_p/,f50/50._ti_p/
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
      mode = 'EWALD'
      call switch (mode)
      if(rank.eq.0.and.tinkerdebug)write(*,*) 'efld0_direct'
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      do i = 1, n
         dscale(i) = 1.0_ti_p
         pscale(i) = 1.0_ti_p
      end do
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        if ((i.eq.0).or.(i.gt.nbloc)) then
          write(iout,1000)
          cycle
        end if
        pdi = pdamp(iipole)
        pti = thole(iipole)
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
           pscale(i12(j,iglob)) = p2scale
        end do
        do j = 1, n13(iglob)
           pscale(i13(j,iglob)) = p3scale
        end do
        do j = 1, n14(iglob)
           pscale(i14(j,iglob)) = p4scale
           do k = 1, np11(iglob)
              if (i14(j,iglob) .eq. ip11(k,iglob))
     &           pscale(i14(j,iglob)) = p4scale * p41scale
           end do
        end do
        do j = 1, n15(iglob)
           pscale(i15(j,iglob)) = p5scale
        end do
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = d1scale
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = d2scale
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = d3scale
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = d4scale
        end do
c
        do kkk = 1, nelst(ii)
          kkpole = elst(kkk,ii)
          kbis = poleloc(kkpole)
          kglob = ipole(kkpole)
          if ((kbis.eq.0).or.(kbis.gt.npolebloc)) then
            write(iout,1000)
            cycle
          end if
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.cutoff2) then
            d  = sqrt(d2)
c
c     calculate the error function damping terms
c
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0_ti_p * aewald**2
            alsq2n = 0.0_ti_p
            if (aewald .gt. 0.0_ti_p)
     &        alsq2n = 1.0_ti_p / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 3
              bfac = real(j+j-1,t_p)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / d2
            end do
c
            d3   = d*d2
            d5   = d3*d2
            ck   = rpole(1,kkpole)
            dkx  = rpole(2,kkpole)
            dky  = rpole(3,kkpole)
            dkz  = rpole(4,kkpole)
            qkxx = rpole(5,kkpole)
            qkxy = rpole(6,kkpole)
            qkxz = rpole(7,kkpole)
            qkyy = rpole(9,kkpole)
            qkyz = rpole(10,kkpole)
            qkzz = rpole(13,kkpole)
            damp = pdi*pdamp(kkpole)
            scale3 = one
            scale5 = one
            scale7 = one
            if (damp.ne.zero) then
              pgamma = min(pti,thole(kkpole))
              damp = -pgamma*(d/damp)**3
              if (damp.gt.-F50) then
                expdamp = exp(damp)
                scale3 = one - expdamp
                scale5 = one - expdamp*(one-damp)
                scale7 = one - expdamp
     &                      *(one-damp + pt6*damp**2)
              end if
            end if
            dsc3 = scale3 * dscale(kglob)
            dsc5 = scale5 * dscale(kglob)
            dsc7 = scale7 * dscale(kglob)
            psc3 = scale3 * pscale(kglob)
            psc5 = scale5 * pscale(kglob)
            psc7 = scale7 * pscale(kglob)
            drr3 = (1.0_ti_p-dsc3) / (d*d2)
            drr5 = 3.0_ti_p * (1.0_ti_p-dsc5) / (d*d2*d2)
            drr7 = 15.0_ti_p * (1.0_ti_p-dsc7) / (d*d2*d2*d2)
            prr3 = (1.0_ti_p-psc3) / (d*d2)
            prr5 = 3.0_ti_p * (1.0_ti_p-psc5) / (d*d2*d2)
            prr7 = 15.0_ti_p * (1.0_ti_p-psc7) / (d*d2*d2*d2)
c
c     compute some intermediate quantities
c
            dir = dix*dx + diy*dy + diz*dz
            qix = qixx*dx + qixy*dy + qixz*dz
            qiy = qixy*dx + qiyy*dy + qiyz*dz
            qiz = qixz*dx + qiyz*dy + qizz*dz
            qir = qix*dx + qiy*dy + qiz*dz
            dkr = dkx*dx + dky*dy + dkz*dz
            qkx = qkxx*dx + qkxy*dy + qkxz*dz
            qky = qkxy*dx + qkyy*dy + qkyz*dz
            qkz = qkxz*dx + qkyz*dy + qkzz*dz
            qkr = qkx*dx + qky*dy + qkz*dz
c
            fim(1) = -dx*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dkx + 2.0_ti_p*bn(2)*qkx
            fim(2) = -dy*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dky + 2.0_ti_p*bn(2)*qky
            fim(3) = -dz*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dkz + 2.0_ti_p*bn(2)*qkz
            fkm(1) = dx*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*dix - 2.0_ti_p*bn(2)*qix
            fkm(2) = dy*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*diy - 2.0_ti_p*bn(2)*qiy
            fkm(3) = dz*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*diz - 2.0_ti_p*bn(2)*qiz
            fid(1) = -dx*(drr3*ck-drr5*dkr+drr7*qkr)
     &                   - drr3*dkx + 2.0_ti_p*drr5*qkx
            fid(2) = -dy*(drr3*ck-drr5*dkr+drr7*qkr)
     &                  - drr3*dky + 2.0_ti_p*drr5*qky
            fid(3) = -dz*(drr3*ck-drr5*dkr+drr7*qkr)
     &                   - drr3*dkz + 2.0_ti_p*drr5*qkz
            fkd(1) = dx*(drr3*ci+drr5*dir+drr7*qir)
     &                  - drr3*dix - 2.0_ti_p*drr5*qix
            fkd(2) = dy*(drr3*ci+drr5*dir+drr7*qir)
     &                  - drr3*diy - 2.0_ti_p*drr5*qiy
            fkd(3) = dz*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0_ti_p*drr5*qiz
            fip(1) = -dx*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dkx + 2.0_ti_p*prr5*qkx
            fip(2) = -dy*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dky + 2.0_ti_p*prr5*qky
            fip(3) = -dz*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dkz + 2.0_ti_p*prr5*qkz
            fkp(1) = dx*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*dix - 2.0_ti_p*prr5*qix
            fkp(2) = dy*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*diy - 2.0_ti_p*prr5*qiy
            fkp(3) = dz*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*diz - 2.0_ti_p*prr5*qiz
            ef(1,1,ipoleloc) = ef(1,1,ipoleloc) + fim(1) - fid(1)
            ef(2,1,ipoleloc) = ef(2,1,ipoleloc) + fim(2) - fid(2)
            ef(3,1,ipoleloc) = ef(3,1,ipoleloc) + fim(3) - fid(3)
            ef(1,2,ipoleloc) = ef(1,2,ipoleloc) + fim(1) - fip(1)
            ef(2,2,ipoleloc) = ef(2,2,ipoleloc) + fim(2) - fip(2)
            ef(3,2,ipoleloc) = ef(3,2,ipoleloc) + fim(3) - fip(3)
c
            ef(1,1,kbis) = ef(1,1,kbis) + fkm(1) - fkd(1)
            ef(2,1,kbis) = ef(2,1,kbis) + fkm(2) - fkd(2)
            ef(3,1,kbis) = ef(3,1,kbis) + fkm(3) - fkd(3)
            ef(1,2,kbis) = ef(1,2,kbis) + fkm(1) - fkp(1)
            ef(2,2,kbis) = ef(2,2,kbis) + fkm(2) - fkp(2)
            ef(3,2,kbis) = ef(3,2,kbis) + fkm(3) - fkp(3)
          end if
        end do
        do j = 1, n12(iglob)
           pscale(i12(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, n13(iglob)
           pscale(i13(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, n14(iglob)
           pscale(i14(j,iglob)) = 1.0_ti_p
        end do
        do j = 1, n15(iglob)
           pscale(i15(j,iglob)) = 1.0_ti_p
        end do
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
c
      deallocate (dscale)
      deallocate (pscale)
c
      return
      end
