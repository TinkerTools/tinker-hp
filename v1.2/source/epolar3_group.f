c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar3_group  --  "group" polarization           ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar3_group" calculates the dipole polarization energy
c     using a subset of atoms in gaz phase
c
c
      subroutine epolar3_group
      use action
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use group
      use inform
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use virial
      use mpi
      implicit none
      integer i,iglob,iipole,iglobgroup
      real*8 a(3,3)
c
      if (deb_Path) write(iout,*), 'epolar3_group '
c
c
c     zero out the polarization energy and derivatives
c
      epgroup = 0.0d0
      if (allocated(uindgroup)) deallocate(uindgroup)
      allocate (uindgroup(3,npolegroup))
      if (allocated(uinpgroup)) deallocate(uinpgroup)
      allocate (uinpgroup(3,npolegroup))

c
      if (npolegroup .eq. 0)  return
c
c     rotates all the multipoles of the group
c
      do i = 1, npolegroup
         iglobgroup = ipolegroup(i)
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
         call rotmat (iipole,iglob,a)
         call rotsite (iipole,a)
      end do
      

      call newinduce_group
c
      call epreal3_group

      return
      end

c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal3_group  --  group polarization energy    ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal3_group" evaluates the 
c     polarization energy due to dipole polarization
c     of a subset of atoms
c

      subroutine epreal3_group
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use chgpot
      use chgpen
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
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      implicit none
      integer j,k,iglob,kglob
      integer ii,kk,kkk,iipole,kkpole
      integer iipolegroup,iglobgroup,kglobgroup
      real*8 e,f
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3
      real*8 sr3,sr5,sr7
      real*8 rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8 dmpik(7)
      real*8 scalek
      real*8, allocatable :: pscale(:)
c
      if (deb_Path) write(iout,*), 'epreal3_group '
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      pscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
c
c     compute the dipole polarization energy component
c
      do ii = 1, npolelocgroup
         iipolegroup = poleglobgroup(ii)
         iglobgroup = ipolegroup(iipolegroup)
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
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
         uix = uindgroup(1,iipolegroup)
         uiy = uindgroup(2,iipolegroup)
         uiz = uindgroup(3,iipolegroup)
         if (use_chgpen) then
            corei = pcore(iipole)
            vali = pval(iipole)
            alphai = palpha(iipole)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = p2scale
            do k = 1, np11(iglob)
               if (i12(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i12(j,iglob)) = p2iscale
            end do
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = p3scale
            do k = 1, np11(iglob)
               if (i13(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i13(j,iglob)) = p3iscale
            end do
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
               if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4iscale
            end do
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
            do k = 1, np11(iglob)
               if (i15(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i15(j,iglob)) = p5iscale
            end do
         end do
         do kkk = iglobgroup+1, npolegroup
            kglobgroup = ipolegroup(kkk)
            kglob = globglobgroup(kglobgroup)
            kkpole = pollist(kglob)
            kk = polelocgroup(kkk)
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            r2 = xr*xr + yr*yr + zr*zr
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
            ukx = uindgroup(1,kkk)
            uky = uindgroup(2,kkk)
            ukz = uindgroup(3,kkk)
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
            diu = dix*ukx + diy*uky + diz*ukz
            qiu = qix*ukx + qiy*uky + qiz*ukz
            uir = uix*xr + uiy*yr + uiz*zr
            dku = dkx*uix + dky*uiy + dkz*uiz
            qku = qkx*uix + qky*uiy + qkz*uiz
            ukr = ukx*xr + uky*yr + ukz*zr
c
c     find the energy value for Thole polarization damping
c
            if (use_thole) then
               call dampthole (iipole,kkpole,7,r,dmpik)
               scalek = pscale(kglob)
               rr3 = f / (r*r2)
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               sr3 = scalek * dmpik(3) * rr3
               sr5 = scalek * dmpik(5) * rr5
               sr7 = scalek * dmpik(7) * rr7
               term1 = ck*uir - ci*ukr + diu + dku
               term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
               term3 = uir*qkr - ukr*qir
               e = term1*sr3 + term2*sr5 + term3*sr7
c
c     find the energy value for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kkpole)
               valk = pval(kkpole)
               alphak = palpha(kkpole)
               call dampdir (r,alphai,alphak,dmpi,dmpk)
               scalek = pscale(kglob)
               rr3 = f * scalek / (r*r2)
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr3i = dmpi(3) * rr3
               rr5i = dmpi(5) * rr5
               rr7i = dmpi(7) * rr7
               rr3k = dmpk(3) * rr3
               rr5k = dmpk(5) * rr5
               rr7k = dmpk(7) * rr7
               rr3 = f / (r*r2)
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr3 = scalek*rr3
               e = uir*(corek*rr3+valk*rr3k)
     &                - ukr*(corei*rr3+vali*rr3i)
     &                + diu*rr3i + dku*rr3k
     &                + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                - dkr*uir*rr5k - dir*ukr*rr5i
     &                + qkr*uir*rr7k - qir*ukr*rr7i
            end if
c
c     compute the energy contribution for this interaction
c
            epgroup = epgroup + e
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
