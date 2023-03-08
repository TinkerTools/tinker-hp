c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_macro.h"
      subroutine efld0_direct(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use bound
      use chgpen
      use couple
      use cutoff
      use domdec
      use ewald
      use group
      use iounit
      use inform ,only: deb_Path
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
      integer i,iglob,kglob,nrhs,ipoleloc,nnelst
      real(t_p)  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole,kbis
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr2
      real(t_p) rr3,rr5,rr7
      real(t_p) rr3i,rr5i,rr7i
      real(t_p) rr3k,rr5k,rr7k
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qiyy,qizz
      real(t_p) qixy,qixz,qiyz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkyy,qkzz
      real(t_p) qkxy,qkxz,qkyz
      real(t_p) dir,dkr
      real(t_p) qix,qiy,qiz,qir
      real(t_p) qkx,qky,qkz,qkr
      real(t_p) fid(3), fip(3)
      real(t_p) fkd(3), fkp(3)
      real(t_p) cutoff2
      real(t_p) corei,corek
      real(t_p) vali,valk
      real(t_p) alphai,alphak
      real(t_p) dmp3,dmp5,dmp7
      real(t_p) dmpi(7),dmpk(7)
      real(t_p) dmpik(7),dmpe(7)
      real(t_p) scalek
      real(t_p), allocatable :: dscale(:)
      real(t_p), allocatable :: pscale(:)
      logical shortrange
      character*11 mode
      character*80 :: RoutineName
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='efld0_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='efld0_direct'
         mode = 'EWALD'
      endif
      if (deb_Path) write(*,*)
     &     ' efld0_direct', use_polarshortreal
c
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      dscale = 1.0d0
      pscale = 1.0d0

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
        ci  = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
        qixx = rpole(5,iipole)
        qixy = rpole(6,iipole)
        qixz = rpole(7,iipole)
        qiyy = rpole(9,iipole)
        qiyz = rpole(10,iipole)
        qizz = rpole(13,iipole)
        if (use_chgpen) then
           corei = pcore(iipole)
           vali = pval(iipole)
           alphai = palpha(iipole)
        end if
c
c     set exclusion coefficients for connected atoms
c
        if (dpequal) then
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
              dscale(i12(j,iglob)) = pscale(i12(j,iglob))
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
              dscale(i13(j,iglob)) = pscale(i13(j,iglob))
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
              dscale(i14(j,iglob)) = pscale(i14(j,iglob))
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
              dscale(i15(j,iglob)) = pscale(i15(j,iglob))
           end do
        else
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
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
        end if
c
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
          kbis = poleloc(kkpole)
          kglob = ipole(kkpole)
          if ((kbis.eq.0).or.(kbis.gt.npolebloc)) then
            write(iout,1000)
            cycle
          end if
          xr = x(kglob) - x(iglob)
          yr = y(kglob) - y(iglob)
          zr = z(kglob) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2.le.cutoff2) then
            r = sqrt(r2)
            rr1 = 1.0d0 / r
            rr2 = rr1 * rr1
            rr3 = rr2 * rr1
            rr5 = 3.0d0 * rr2 * rr3
            rr7 = 5.0d0 * rr2 * rr5
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
c
c     intermediates involving moments and separation distance
c
            dir = dix*xr + diy*yr + diz*zr
            qix = qixx*xr + qixy*yr + qixz*zr
            qiy = qixy*xr + qiyy*yr + qiyz*zr
            qiz = qixz*xr + qiyz*yr + qizz*zr
            qir = qix*xr + qiy*yr + qiz*zr
            dkr = dkx*xr + dky*yr + dkz*zr
            qkx = qkxx*xr + qkxy*yr + qkxz*zr
            qky = qkxy*xr + qkyy*yr + qkyz*zr
            qkz = qkxz*xr + qkyz*yr + qkzz*zr
            qkr = qkx*xr + qky*yr + qkz*zr
c
c     calculate real space Ewald error function damping
c
            call dampewald (7,r,r2,1.0d0,dmpe)
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               call dampthole (iipole,kkpole,7,r,dmpik)
               scalek = dscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fid(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fid(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fid(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkd(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkd(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkd(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz
               scalek = pscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fip(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fip(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fip(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkp(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkp(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkp(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kkpole)
               valk = pval(kkpole)
               alphak = palpha(kkpole)
               call dampdir (r,alphai,alphak,dmpi,dmpk)
               scalek = dscale(kglob)
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fid(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fid(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fid(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkd(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkd(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkd(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
               scalek = pscale(kglob)
               rr3 = rr2 * rr1
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fip(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fip(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fip(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkp(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkp(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkp(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
            end if
c
c     increment the field at each site due to this interaction
c
            do j = 1, 3
               ef(j,1,ipoleloc) = ef(j,1,ipoleloc) + fid(j)
               ef(j,1,kbis) = ef(j,1,kbis) + fkd(j)
               ef(j,2,ipoleloc) = ef(j,2,ipoleloc) + fip(j)
               ef(j,2,kbis) = ef(j,2,kbis) + fkp(j)
            end do
          end if
        end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = 1.0d0
               dscale(i12(j,iglob)) = 1.0d0
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = 1.0d0
               dscale(i13(j,iglob)) = 1.0d0
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = 1.0d0
               dscale(i14(j,iglob)) = 1.0d0
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = 1.0d0
               dscale(i15(j,iglob)) = 1.0d0
            end do
         else
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = 1.0d0
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = 1.0d0
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = 1.0d0
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = 1.0d0
            end do
            do j = 1, np11(iglob)
               dscale(ip11(j,iglob)) = 1.0d0
            end do
            do j = 1, np12(iglob)
               dscale(ip12(j,iglob)) = 1.0d0
            end do
            do j = 1, np13(iglob)
               dscale(ip13(j,iglob)) = 1.0d0
            end do
            do j = 1, np14(iglob)
               dscale(ip14(j,iglob)) = 1.0d0
            end do
         end if
      end do
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
