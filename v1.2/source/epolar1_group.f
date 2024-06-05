c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1_group  --  "group" polarization           ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1_group" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using
c     a subset of atoms in gaz phase
c
c
      subroutine epolar1_group
      use action
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
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
      if (deb_Path) write(iout,*), 'epolar1_group '
c
c
c     zero out the polarization energy and derivatives
c
      epgroup = 0.0d0
      vir_group = 0d0
      if (allocated(depgroup)) deallocate(depgroup)
      allocate (depgroup(3,natgroup))
      depgroup = 0d0
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
      call epreal1_group
c
c     communicate the forces
c
      call commforcesgroup

      return
      end
      
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal1_group  --  group polarization derivs    ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal1_group" evaluates the 
c     polarization energy and gradient due to dipole polarization
c     of a subset of atoms
c
c
      subroutine epreal1_group
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
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
      use virial
      use mpi
      implicit none
      integer i,j,k,iglob,kglob
      integer ii,kk,kkk,iipole,kkpole
      integer iipolegroup,iglobgroup,kglobgroup
      integer ilocgroup,klocgroup
      integer ix,iy,iz
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
      real*8 term4,term5
      real*8 term6,term7,term8
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 tep(3)
      real*8 tix3,tiy3,tiz3
      real*8 tix5,tiy5,tiz5
      real*8 tkx3,tky3,tkz3
      real*8 tkx5,tky5,tkz5
      real*8 tuir,tukr
      real*8 tixx,tiyy,tizz
      real*8 tixy,tixz,tiyz
      real*8 tkxx,tkyy,tkzz
      real*8 tkxy,tkxz,tkyz
      real*8 scaled,scalep,scaleu

      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
c
      if (deb_Path) write(iout,*), 'eprea1_group '
c

c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (ufld(3,npolegroup))
      allocate (dufld(6,npolegroup))
c
c     set exclusion coefficients and arrays to store fields
c
      pscale = 1.0d0
      dscale = 1.0d0
      uscale = 1.0d0
      ufld = 0.0d0
      dufld = 0.0d0

c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
c
c
c     compute the dipole polarization gradient components
c
      do ii = 1, npolelocgroup
         iipolegroup = poleglobgroup(ii)
         iglobgroup = ipolegroup(iipolegroup)
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
         ilocgroup = locgroup(iglobgroup)
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
         uix = uindgroup(1,iipolegroup)
         uiy = uindgroup(2,iipolegroup)
         uiz = uindgroup(3,iipolegroup)
         uixp = uinpgroup(1,iipolegroup)
         uiyp = uinpgroup(2,iipolegroup)
         uizp = uinpgroup(3,iipolegroup)
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
         do kkk = iipolegroup+1, npolegroup
            kglobgroup = ipolegroup(kkk)
            kglob = globglobgroup(kglobgroup)
            kkpole = pollist(kglob)
            kk = polelocgroup(kkk)
            klocgroup = locgroup(kglobgroup)
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

            tix3 = psr3*ukx + dsr3*ukxp
            tiy3 = psr3*uky + dsr3*ukyp
            tiz3 = psr3*ukz + dsr3*ukzp
            tkx3 = psr3*uix + dsr3*uixp
            tky3 = psr3*uiy + dsr3*uiyp
            tkz3 = psr3*uiz + dsr3*uizp
            tuir = -psr5*urk - dsr5*urkp
            tukr = -psr5*uri - dsr5*urip
            ufld(1,ii) = ufld(1,ii) + tix3 + xr*tuir
            ufld(2,ii) = ufld(2,ii) + tiy3 + yr*tuir
            ufld(3,ii) = ufld(3,ii) + tiz3 + zr*tuir
            ufld(1,kk) = ufld(1,kk) + tkx3 + xr*tukr
            ufld(2,kk) = ufld(2,kk) + tky3 + yr*tukr
            ufld(3,kk) = ufld(3,kk) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
            tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
            tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
            tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
            tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
            tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
            tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
            tuir = -psr7*urk - dsr7*urkp
            tukr = -psr7*uri - dsr7*urip
            dufld(1,ii) = dufld(1,ii) + xr*tix5 + xr*xr*tuir
            dufld(2,ii) = dufld(2,ii) + xr*tiy5 + yr*tix5
     &                      + 2.0d0*xr*yr*tuir
            dufld(3,ii) = dufld(3,ii) + yr*tiy5 + yr*yr*tuir
            dufld(4,ii) = dufld(4,ii) + xr*tiz5 + zr*tix5
     &                      + 2.0d0*xr*zr*tuir
            dufld(5,ii) = dufld(5,ii) + yr*tiz5 + zr*tiy5
     &                      + 2.0d0*yr*zr*tuir
            dufld(6,ii) = dufld(6,ii) + zr*tiz5 + zr*zr*tuir
            dufld(1,kk) = dufld(1,kk) - xr*tkx5 - xr*xr*tukr
            dufld(2,kk) = dufld(2,kk) - xr*tky5 - yr*tkx5
     &                      - 2.0d0*xr*yr*tukr
            dufld(3,kk) = dufld(3,kk) - yr*tky5 - yr*yr*tukr
            dufld(4,kk) = dufld(4,kk) - xr*tkz5 - zr*tkx5
     &                      - 2.0d0*xr*zr*tukr
            dufld(5,kk) = dufld(5,kk) - yr*tkz5 - zr*tky5
     &                      - 2.0d0*yr*zr*tukr
            dufld(6,kk) = dufld(6,kk) - zr*tkz5 - zr*zr*tukr
c
c     get the field gradient for direct polarization force
c
            term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
            term2 = (sc3+sc5)*rr5*xr - rc3(1)
            term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
            term4 = 2.0d0 * sc5 * rr5
            term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
            term6 = xr * (sc7*rr9*xr-rc7(1))
            tixx = ci*term1 + dix*term2 - dri*term3
     &                - qixx*term4 + qrix*term5 - qrri*term6
     &                + (qriy*yr+qriz*zr)*sc7*rr7
            tkxx = ck*term1 - dkx*term2 + drk*term3
     &                - qkxx*term4 + qrkx*term5 - qrrk*term6
     &                + (qrky*yr+qrkz*zr)*sc7*rr7
            term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
            term2 = (sc3+sc5)*rr5*yr - rc3(2)
            term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
            term4 = 2.0d0 * sc5 * rr5
            term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
            term6 = yr * (sc7*rr9*yr-rc7(2))
            tiyy = ci*term1 + diy*term2 - dri*term3
     &                - qiyy*term4 + qriy*term5 - qrri*term6
     &                + (qrix*xr+qriz*zr)*sc7*rr7
            tkyy = ck*term1 - dky*term2 + drk*term3
     &                - qkyy*term4 + qrky*term5 - qrrk*term6
     &                + (qrkx*xr+qrkz*zr)*sc7*rr7
            term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
            term2 = (sc3+sc5)*rr5*zr - rc3(3)
            term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
            term4 = 2.0d0 * sc5 * rr5
            term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
            term6 = zr * (sc7*rr9*zr-rc7(3))
            tizz = ci*term1 + diz*term2 - dri*term3
     &                - qizz*term4 + qriz*term5 - qrri*term6
     &                + (qrix*xr+qriy*yr)*sc7*rr7
            tkzz = ck*term1 - dkz*term2 + drk*term3
     &                - qkzz*term4 + qrkz*term5 - qrrk*term6
     &                + (qrkx*xr+qrky*yr)*sc7*rr7
            term2 = sc3*rr5*xr - rc3(1)
            term1 = yr * term2
            term3 = sc5 * rr5 * yr
            term4 = yr * (sc5*rr7*xr-rc5(1))
            term5 = 2.0d0 * sc5 * rr5
            term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
            term7 = 2.0d0 * sc7 * rr7 * yr
            term8 = yr * (sc7*rr9*xr-rc7(1))
            tixy = -ci*term1 + diy*term2 + dix*term3
     &                - dri*term4 - qixy*term5 + qriy*term6
     &                + qrix*term7 - qrri*term8
            tkxy = -ck*term1 - dky*term2 - dkx*term3
     &                + drk*term4 - qkxy*term5 + qrky*term6
     &                + qrkx*term7 - qrrk*term8
            term2 = sc3*rr5*xr - rc3(1)
            term1 = zr * term2
            term3 = sc5 * rr5 * zr
            term4 = zr * (sc5*rr7*xr-rc5(1))
            term5 = 2.0d0 * sc5 * rr5
            term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
            term7 = 2.0d0 * sc7 * rr7 * zr
            term8 = zr * (sc7*rr9*xr-rc7(1))
            tixz = -ci*term1 + diz*term2 + dix*term3
     &                - dri*term4 - qixz*term5 + qriz*term6
     &                + qrix*term7 - qrri*term8
            tkxz = -ck*term1 - dkz*term2 - dkx*term3
     &                + drk*term4 - qkxz*term5 + qrkz*term6
     &                + qrkx*term7 - qrrk*term8
            term2 = sc3*rr5*yr - rc3(2)
            term1 = zr * term2
            term3 = sc5 * rr5 * zr
            term4 = zr * (sc5*rr7*yr-rc5(2))
            term5 = 2.0d0 * sc5 * rr5
            term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
            term7 = 2.0d0 * sc7 * rr7 * zr
            term8 = zr * (sc7*rr9*yr-rc7(2))
            tiyz = -ci*term1 + diz*term2 + diy*term3
     &                - dri*term4 - qiyz*term5 + qriz*term6
     &                + qriy*term7 - qrri*term8
            tkyz = -ck*term1 - dkz*term2 - dky*term3
     &                + drk*term4 - qkyz*term5 + qrkz*term6
     &                + qrky*term7 - qrrk*term8
c
c     get the dEd/dR terms for Thole direct polarization force
c
           depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &               - tkxx*uixp - tkxy*uiyp - tkxz*uizp
           depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &               - tkxy*uixp - tkyy*uiyp - tkyz*uizp
           depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &               - tkxz*uixp - tkyz*uiyp - tkzz*uizp
           frcx = scaled * depx
           frcy = scaled * depy
           frcz = scaled * depz
c
c     get the dEp/dR terms for Thole direct polarization force
c
           depx = tixx*ukx + tixy*uky + tixz*ukz
     &               - tkxx*uix - tkxy*uiy - tkxz*uiz
           depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &               - tkxy*uix - tkyy*uiy - tkyz*uiz
           depz = tixz*ukx + tiyz*uky + tizz*ukz
     &               - tkxz*uix - tkyz*uiy - tkzz*uiz
           frcx = frcx + scalep*depx
           frcy = frcy + scalep*depy
           frcz = frcz + scalep*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
           term1 = (sc3+sc5) * rr5
           term2 = term1*xr - rc3(1)
           term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
           tixx = uix*term2 + uri*term3
           tkxx = ukx*term2 + urk*term3
           term2 = term1*yr - rc3(2)
           term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
           tiyy = uiy*term2 + uri*term3
           tkyy = uky*term2 + urk*term3
           term2 = term1*zr - rc3(3)
           term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
           tizz = uiz*term2 + uri*term3
           tkzz = ukz*term2 + urk*term3
           term1 = sc5 * rr5 * yr
           term2 = sc3*rr5*xr - rc3(1)
           term3 = yr * (sc5*rr7*xr-rc5(1))
           tixy = uix*term1 + uiy*term2 - uri*term3
           tkxy = ukx*term1 + uky*term2 - urk*term3
           term1 = sc5 * rr5 * zr
           term3 = zr * (sc5*rr7*xr-rc5(1))
           tixz = uix*term1 + uiz*term2 - uri*term3
           tkxz = ukx*term1 + ukz*term2 - urk*term3
           term2 = sc3*rr5*yr - rc3(2)
           term3 = zr * (sc5*rr7*yr-rc5(2))
           tiyz = uiy*term1 + uiz*term2 - uri*term3
           tkyz = uky*term1 + ukz*term2 - urk*term3
           depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &               + tkxx*uixp + tkxy*uiyp + tkxz*uizp
           depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &               + tkxy*uixp + tkyy*uiyp + tkyz*uizp
           depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &               + tkxz*uixp + tkyz*uiyp + tkzz*uizp
           frcx = frcx + scaleu*depx
           frcy = frcy + scaleu*depy
           frcz = frcz + scaleu*depz
c
c     increment gradient and virial due to Cartesian forces
c
            depgroup(1,ilocgroup) = depgroup(1,ilocgroup) + frcx
            depgroup(2,ilocgroup) = depgroup(2,ilocgroup) + frcy
            depgroup(3,ilocgroup) = depgroup(3,ilocgroup) + frcz
            depgroup(1,klocgroup) = depgroup(1,klocgroup) - frcx
            depgroup(2,klocgroup) = depgroup(2,klocgroup) - frcy
            depgroup(3,klocgroup) = depgroup(3,klocgroup) - frcz
            vxx = -xr * frcx
            vxy = -0.5d0 * (yr*frcx+xr*frcy)
            vxz = -0.5d0 * (zr*frcx+xr*frcz)
            vyy = -yr * frcy
            vyz = -0.5d0 * (zr*frcy+yr*frcz)
            vzz = -zr * frcz
            vir_group(1,1) = vir_group(1,1) + vxx
            vir_group(2,1) = vir_group(2,1) + vxy
            vir_group(3,1) = vir_group(3,1) + vxz
            vir_group(1,2) = vir_group(1,2) + vxy
            vir_group(2,2) = vir_group(2,2) + vyy
            vir_group(3,2) = vir_group(3,2) + vyz
            vir_group(1,3) = vir_group(1,3) + vxz
            vir_group(2,3) = vir_group(2,3) + vyz
            vir_group(3,3) = vir_group(3,3) + vzz
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
c     torque is induced field and gradient cross permanent moments
c
      do i = 1, npolegroup
         iglob = globglobgroup(ipolegroup(i))
         iipole = pollist(iglob)
         ii = polelocgroup(i)
         dix = rpole(2,iipole)
         diy = rpole(3,iipole)
         diz = rpole(4,iipole)
         qixx = rpole(5,iipole)
         qixy = rpole(6,iipole)
         qixz = rpole(7,iipole)
         qiyy = rpole(9,iipole)
         qiyz = rpole(10,iipole)
         qizz = rpole(13,iipole)
         tep(1) = diz*ufld(2,ii) - diy*ufld(3,ii)
     &               + qixz*dufld(2,ii) - qixy*dufld(4,ii)
     &               + 2.0d0*qiyz*(dufld(3,ii)-dufld(6,ii))
     &               + (qizz-qiyy)*dufld(5,ii)
         tep(2) = dix*ufld(3,ii) - diz*ufld(1,ii)
     &               - qiyz*dufld(2,ii) + qixy*dufld(5,ii)
     &               + 2.0d0*qixz*(dufld(6,ii)-dufld(1,ii))
     &               + (qixx-qizz)*dufld(4,ii)
         tep(3) = diy*ufld(1,ii) - dix*ufld(2,ii)
     &               + qiyz*dufld(4,ii) - qixz*dufld(5,ii)
     &               + 2.0d0*qixy*(dufld(1,ii)-dufld(3,ii))
     &               + (qiyy-qixx)*dufld(2,ii)
         call torque_group(iipole,tep,fix,fiy,fiz,depgroup)
         iz = zaxis(iipole)
         ix = xaxis(iipole)
         iy = abs(yaxis(iipole))
         if (iz .eq. 0)  iz = iglob
         if (ix .eq. 0)  ix = iglob
         if (iy .eq. 0)  iy = iglob
         xiz = x(iz) - x(iglob)
         yiz = y(iz) - y(iglob)
         ziz = z(iz) - z(iglob)
         xix = x(ix) - x(iglob)
         yix = y(ix) - y(iglob)
         zix = z(ix) - z(iglob)
         xiy = x(iy) - x(iglob)
         yiy = y(iy) - y(iglob)
         ziy = z(iy) - z(iglob)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir_group(1,1) = vir_group(1,1) + vxx
         vir_group(2,1) = vir_group(2,1) + vxy
         vir_group(3,1) = vir_group(3,1) + vxz
         vir_group(1,2) = vir_group(1,2) + vxy
         vir_group(2,2) = vir_group(2,2) + vyy
         vir_group(3,2) = vir_group(3,2) + vyz
         vir_group(1,3) = vir_group(1,3) + vxz
         vir_group(2,3) = vir_group(2,3) + vyz
         vir_group(3,3) = vir_group(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (ufld)
      deallocate (dufld)
      return
      end

