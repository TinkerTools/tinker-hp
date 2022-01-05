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
c
      subroutine epreal3_group
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use group
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
      use virial
      use mpi
      implicit none
      integer j,k,iglob,kglob
      integer ii,kk,kkk,iipole,kkpole
      integer iglobgroup,kglobgroup
      real*8 e,f
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 usr3,usr5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 uixp,uiyp,uizp
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 ukxp,ukyp,ukzp
      real*8 dri,drk
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrri,qrrk
      real*8 duik,quik
      real*8 uri,urip
      real*8 urk,urkp

      real*8 term1,term2,term3
      real*8 scaled,scalep,scaleu
      real*8 rc3(3),rc5(3),rc7(3)

      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)

c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
c
c     set exclusion coefficients and arrays to store fields
c
      pscale = 1.0d0
      dscale = 1.0d0
      uscale = 1.0d0

c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
c
c
c     compute the dipole polarization gradient components
c
      do ii = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(ii))
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
         pdi = pdamp(iipole)
         pti = thole(iipole)
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
         uix = uindgroup(1,iglobgroup)
         uiy = uindgroup(2,iglobgroup)
         uiz = uindgroup(3,iglobgroup)
         uixp = uinpgroup(1,iglobgroup)
         uiyp = uinpgroup(2,iglobgroup)
         uizp = uinpgroup(3,iglobgroup)
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
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
         end do
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = d1scale
            uscale(ip11(j,iglob)) = u1scale
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = d2scale
            uscale(ip12(j,iglob)) = u2scale
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = d3scale
            uscale(ip13(j,iglob)) = u3scale
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = d4scale
            uscale(ip14(j,iglob)) = u4scale
         end do
c
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
            ukxp = uinpgroup(1,kkk)
            ukyp = uinpgroup(2,kkk)
            ukzp = uinpgroup(3,kkk)
c
c     get reciprocal distance terms for this interaction
c
            rr1 = f / r
            rr3 = rr1 / r2
            rr5 = 3.0d0 * rr3 / r2
            rr7 = 5.0d0 * rr5 / r2
            rr9 = 7.0d0 * rr7 / r2
c
c     apply Thole polarization damping to scale factors
c
            sc3 = 1.0d0
            sc5 = 1.0d0
            sc7 = 1.0d0
            do j = 1, 3
               rc3(j) = 0.0d0
               rc5(j) = 0.0d0
               rc7(j) = 0.0d0
            end do
            damp = pdi * pdamp(kkpole)
            if (damp .ne. 0.0d0) then
              pgamma = min(pti,thole(kk))
              if (pgamma .eq. 0.0d0) then
                 pgamma = max(pti,thole(kk))
              end if
              if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                 damp = pgamma * (r/damp)**3
                 if (damp .lt. 50.0d0) then
                    expdamp = exp(-damp)
                    sc3 = 1.0d0 - expdamp
                    sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                    sc7 = 1.0d0 - expdamp*(1.0d0+damp
     &                               +0.6d0*damp**2)
                    temp3 = damp * expdamp * rr5
                    temp5 = 3.0d0 * damp / r2
                    temp7 = (-1.0d0+3.0d0*damp) / r2
                    rc3(1) = xr * temp3
                    rc3(2) = yr * temp3
                    rc3(3) = zr * temp3
                    rc5(1) = rc3(1) * temp5
                    rc5(2) = rc3(2) * temp5
                    rc5(3) = rc3(3) * temp5
                    rc7(1) = rc5(1) * temp7
                    rc7(2) = rc5(2) * temp7
                    rc7(3) = rc5(3) * temp7
                 end if
              end if
            end if
cc
c     intermediates involving moments and distance separation
c
            dri = dix*xr + diy*yr + diz*zr
            drk = dkx*xr + dky*yr + dkz*zr
            qrix = qixx*xr + qixy*yr + qixz*zr
            qriy = qixy*xr + qiyy*yr + qiyz*zr
            qriz = qixz*xr + qiyz*yr + qizz*zr
            qrkx = qkxx*xr + qkxy*yr + qkxz*zr
            qrky = qkxy*xr + qkyy*yr + qkyz*zr
            qrkz = qkxz*xr + qkyz*yr + qkzz*zr
            qrri = qrix*xr + qriy*yr + qriz*zr
            qrrk = qrkx*xr + qrky*yr + qrkz*zr
            uri = uix*xr + uiy*yr + uiz*zr
            urk = ukx*xr + uky*yr + ukz*zr
            urip = uixp*xr + uiyp*yr + uizp*zr
            urkp = ukxp*xr + ukyp*yr + ukzp*zr
            duik = dix*ukx + diy*uky + diz*ukz
     &                + dkx*uix + dky*uiy + dkz*uiz
            quik = qrix*ukx + qriy*uky + qriz*ukz
     &                - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
            term1 = ck*uri - ci*urk + duik
            term2 = 2.0d0*quik - uri*drk - dri*urk
            term3 = uri*qrrk - urk*qrri
c
c     intermediates involving Thole damping and scale factors
c
            scaled = dscale(kglob)
            scalep = pscale(kglob)
            scaleu = uscale(kglob)
c
c     compute the full undamped energy for this interaction
c
            if (molcule(iglob) .ne. molcule(kglob))
     &         einter = einter + e
c
c     intermediates involving Thole damping and scale factors
c
            psr3 =  scalep* sc3*rr3
            psr5 =  scalep* sc5*rr5
            psr7 =  scalep* sc7*rr7
            dsr3 =  scaled* sc3*rr3
            dsr5 =  scaled* sc5*rr5
            dsr7 =  scaled* sc7*rr7
            usr3 =  scaleu* sc3*rr3
            usr5 =  scaleu* sc5*rr5
c
c     compute the energy contribution for this interaction
c
            e = term1*psr3 + term2*psr5 + term3*psr7

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
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = 1.0d0
            uscale(ip11(j,iglob)) = 1.0d0
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = 1.0d0
            uscale(ip12(j,iglob)) = 1.0d0
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = 1.0d0
            uscale(ip13(j,iglob)) = 1.0d0
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = 1.0d0
            uscale(ip14(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end
