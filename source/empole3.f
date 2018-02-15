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
c
      subroutine empole3
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'potent.i'
      include 'mpif.h'
      real*8 time0,time1
c
       time0 = mpi_wtime()
c
c     choose the method for summing over multipole interactions
c
      call empole3c
      time1 = mpi_wtime()
c
c     zero out energy terms and analysis which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
      end if
      if (.not. use_polar) then
         ep = 0.0d0
      end if
      return
      end
c
c     subroutine empole3c : computes the atomic multipole and dipole polarizability
c     interaction energy using a particle mesh Ewald summation
c
      subroutine empole3c
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'potent.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,ii,iglob,iipole
      integer ierr
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
      eintra = 0.0d0
c
c     zero out the multipole and polarization energies
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      if (use_pmecore) then
        call newinduce_pme
      else
        call newinduce_pme2
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  call emrecip
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        call ereal3c(eintra)
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
           uix = uind(1,iipole)
           uiy = uind(2,iipole)
           uiz = uind(3,iipole)
           cii = ci*ci
           dii = dix*dix + diy*diy + diz*diz
           qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &              + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
           uii = dix*uix + diy*uiy + diz*uiz
           e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
           ei = fterm * term * uii / 3.0d0
           nem = nem + 1
           nep = nep + 1
           em = em + e
           ep = ep + ei
           aem(i) = aem(i) + e
           aep(i) = aep(i) + ei
        end do
c
c       compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0d0
           yd = 0.0d0
           zd = 0.0d0
           xu = 0.0d0
           yu = 0.0d0
           zu = 0.0d0
           do ii = 1, npoleloc
              iipole = poleglob(ii)
              iglob = ipole(iipole)
              i = loc(iglob)
              dix = rpole(2,iipole)
              diy = rpole(3,iipole)
              diz = rpole(4,iipole)
              uix = uind(1,iipole)
              uiy = uind(2,iipole)
              uiz = uind(3,iipole)
              xd = xd + dix + rpole(1,iipole)*x(iglob)
              yd = yd + diy + rpole(1,iipole)*y(iglob)
              zd = zd + diz + rpole(1,iipole)*z(iglob)
              xu = xu + uix
              yu = yu + uiy
              zu = zu + uiz
           end do
           term = (2.0d0/3.0d0) * f * (pi/volbox)
           nem = nem + 1
           nep = nep + 1
           em = em + term*(xd*xd+yd*yd+zd*zd)
           ep = ep + term*(xd*xu+yd*yu+zd*zu)
        end if
c
c     ntermolecular energy is total minus intramolecular part
c
        einter = einter + em + ep - eintra
      end if
      return
      end
c
c
c
c      subroutine ereal3c   real space mpole analysis
c
c
c     "ereal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and dipole
c     polarizability and partitions the energy among the atoms
c
c
      subroutine ereal3c(eintra)
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'mpif.h'
      include 'openmp.i'
      integer i,j,k,m,inl
      integer ii,kkk,iipole,kkpole
      integer iglob, kglob
      real*8 e,ei,eintra
      real*8 f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8 eptemp
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      logical muse,puse
      character*6 mode
      external erfc
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      eptemp = 0.0d0
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      mscale = 1.0d0
      pscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         inl = polelocnl(iipole)
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
         uix = uind(1,iipole)
         uiy = uind(2,iipole)
         uiz = uind(3,iipole)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = m2scale
            pscale(i12(j,iglob)) = p2scale
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = m3scale
            pscale(i13(j,iglob)) = p3scale
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = m4scale
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
                if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = m5scale
            pscale(i15(j,iglob)) = p5scale
         end do
         do kkk = 1,nelst(inl)
           kglob = elst(kkk,inl)
           kkpole = pollist(kglob)
           xr = x(kglob) - x(iglob)
           yr = y(kglob) - y(iglob)
           zr = z(kglob) - z(iglob)
           call image (xr,yr,zr)
           r2 = xr*xr + yr* yr + zr*zr
           if (r2 .le. off2) then
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
              ukx = uind(1,kkpole)
              uky = uind(2,kkpole)
              ukz = uind(3,kkpole)
c
c     calculate the error function damping terms
c
              ralpha = aewald * r
              bn(0) = erfc(ralpha) / r
              alsq2 = 2.0d0 * aewald**2
              alsq2n = 0.0d0
              if (aewald .gt. 0.0d0)
     &           alsq2n = 1.0d0 / (sqrtpi*aewald)
              exp2a = exp(-ralpha**2)
              do m = 1, 4
                 bfac = dble(m+m-1)
                 alsq2n = alsq2 * alsq2n
                 bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
              end do
c
c     construct some intermediate quadrupole values
c
              qix = qixx*xr + qixy*yr + qixz*zr
              qiy = qixy*xr + qiyy*yr + qiyz*zr
              qiz = qixz*xr + qiyz*yr + qizz*zr
              qkx = qkxx*xr + qkxy*yr + qkxz*zr
              qky = qkxy*xr + qkyy*yr + qkyz*zr
              qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
              sc(2) = dix*dkx + diy*dky + diz*dkz
              sc(3) = dix*xr + diy*yr + diz*zr
              sc(4) = dkx*xr + dky*yr + dkz*zr
              sc(5) = qix*xr + qiy*yr + qiz*zr
              sc(6) = qkx*xr + qky*yr + qkz*zr
              sc(7) = qix*dkx + qiy*dky + qiz*dkz
              sc(8) = qkx*dix + qky*diy + qkz*diz
              sc(9) = qix*qkx + qiy*qky + qiz*qkz
              sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                    + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
              sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                    + diy*uky + uiz*dkz + diz*ukz
              sci(3) = uix*xr + uiy*yr + uiz*zr
              sci(4) = ukx*xr + uky*yr + ukz*zr
              sci(7) = qix*ukx + qiy*uky + qiz*ukz
              sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
              gl(0) = ci*ck
              gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
              gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                   + 2.0d0*(sc(7)-sc(8)+sc(10))
              gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
              gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
              gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
              gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                    - sc(3)*sci(4)
              gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
              e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &               + gl(3)*bn(3) + gl(4)*bn(4)
              ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
              rr1 = 1.0d0 / r
              rr3 = rr1 / r2
              rr5 = 3.0d0 * rr3 / r2
              rr7 = 5.0d0 * rr5 / r2
              rr9 = 7.0d0 * rr7 / r2
              scale3 = 1.0d0*pscale(kglob)
              scale5 = 1.0d0*pscale(kglob)
              scale7 = 1.0d0*pscale(kglob)
              damp = pdi * pdamp(kkpole)
              if (damp .ne. 0.0d0) then
                 pgamma = min(pti,thole(kkpole))
                 damp = -pgamma * (r/damp)**3
                 if (damp .gt. -50.0d0) then
                    expdamp = exp(damp)
                    scale3 = scale3 * (1.0d0-expdamp)
                    scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                    scale7 = scale7 * (1.0d0-(1.0d0-damp
     &                                   +0.6d0*damp**2)*expdamp)
                 end if
              end if
              efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                  + gl(3)*rr7 + gl(4)*rr9
              eifix = gli(1)*rr3*(1.0d0-scale3)
     &                   + gli(2)*rr5*(1.0d0-scale5)
     &                   + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
              e = e - efix*(1.0d0-mscale(kglob))
              ei = ei - eifix
              e = f * e
              ei = 0.5d0 * f * ei
c
c     increment the overall multipole and polarization energies
c
              muse = use_mpole
              puse = use_polar
              if (muse)  nem = nem + 1
              if (puse)  nep = nep + 1
              em = em + e
              eptemp = eptemp + ei
              aem(i) = aem(i) + e
              aep(i) = aep(i) + ei
c
c     increment the total intramolecular energy
c
              efix = f * efix * mscale(kglob)
              eifix = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                   + gli(3)*rr7*scale7
              eifix = 0.5d0 * f * eifix
              if (molcule(iglob) .eq. molcule(kglob)) then
                 eintra = eintra + efix + eifix
              end if
           end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0d0
            pscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0d0
            pscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0d0
            pscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0d0
            pscale(i15(j,iglob)) = 1.0d0
         end do
      end do
      ep = ep + eptemp
c
      deallocate (mscale)
      deallocate (pscale)
      return
      end
