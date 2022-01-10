c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      use domdec
      use energi
      use potent
      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empole1c
c
c     zero out energy and derivative terms which are not in use
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
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1c
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use cutoff
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use timestat
      use virial
      use mpi
      implicit none
      integer i,j,ii
      integer iipole,iglob
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
      real*8 time0,time1
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      dem = 0d0
      if (npole .eq. 0)  return
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
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
        time0 = mpi_wtime()
        call emrecip1
        time1 = mpi_wtime()
        timerec = timerec + time1 - time0
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        time0 = mpi_wtime()
        call emreal1c
        if ((use_emtp).or.(use_emtporig)) call dercorrecemtp(mpolecut)
        time1 = mpi_wtime()
        timereal = timereal + time1 - time0
c
c     compute the Ewald self-energy term over all the atoms
c
        term = 2.0d0 * aewald * aewald
        fterm = -f * aewald / sqrtpi
        do i = 1, npoleloc
           iipole = poleglob(i)
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
     &              + qixx*qixx + qiyy*qiyy + qizz*qizz
           e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
           em = em + e
        end do
c
c       compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0d0
           yd = 0.0d0
           zd = 0.0d0
           do i = 1, npoleloc
              iipole = poleglob(i)
              iglob = ipole(iipole)
              xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
              yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
              zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
           end do
           term = (2.0d0/3.0d0) * f * (pi/volbox)
           em = em + term*(xd*xd+yd*yd+zd*zd)
           do ii = 1, npoleloc
              iipole = poleglob(ii)
              iglob = ipole(iipole)
              i = loc(iglob)
              dem(1,i) = dem(1,i) + 2.0d0*term*rpole(1,iipole)*xd
              dem(2,i) = dem(2,i) + 2.0d0*term*rpole(1,iipole)*yd
              dem(3,i) = dem(3,i) + 2.0d0*term*rpole(1,iipole)*zd
           end do
           xdfield = -2.0d0 * term * xd
           ydfield = -2.0d0 * term * yd
           zdfield = -2.0d0 * term * zd
           do i = 1, npoleloc
              iipole = poleglob(i)
              trq(1) = rpole(3,iipole)*zdfield - rpole(4,iipole)*ydfield
              trq(2) = rpole(4,iipole)*xdfield - rpole(2,iipole)*zdfield
              trq(3) = rpole(2,iipole)*ydfield - rpole(3,iipole)*xdfield
              call torque(iipole,trq,frcx,frcy,frcz,dem)
           end do
c
c       boundary correction to virial due to overall cell dipole
c
           xd = 0.0d0
           yd = 0.0d0
           zd = 0.0d0
           xq = 0.0d0
           yq = 0.0d0
           zq = 0.0d0
           do i = 1, npoleloc
              iipole = poleglob(i)
              iglob = ipole(iipole)
              xd = xd + rpole(2,iipole)
              yd = yd + rpole(3,iipole)
              zd = zd + rpole(4,iipole)
              xq = xq + rpole(1,iipole)*x(iglob)
              yq = yq + rpole(1,iipole)*y(iglob)
              zq = zq + rpole(1,iipole)*z(iglob)
           end do
           xv = xd * xq
           yv = yd * yq
           zv = zd * zq
           vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                        + xq*xq + yq*yq + zq*zq)
           vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
           vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
           vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
           vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
           vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
           vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
           vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
           vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
           vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
        end if
      end if
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreal1c
      use atmlst
      use atoms
      use atmtyp
      use bound
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      use ewald
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      use virial
      use mpi
      implicit none
      integer i,j,k,iglob,kglob,kbis
      integer ii,kk,kkk,iipole,kkpole
      integer inl
      integer iax,iay,iaz
      real*8 e,efull,de,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 bn(0:5)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      character*6 mode
      external erfc
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate')
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,nbloc))
c
c     initialize connected atom scaling and torque arrays
c
      mscale = 1.0d0
      tem = 0.0d0
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
         if ((i.eq.0).or.(i.gt.nbloc)) then
           write(iout,1000)
           cycle
         end if
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
         do kkk = 1, nelst(ii)
            kkpole = elst(kkk,ii)
            kglob = ipole(kkpole)
            kbis = loc(kglob)
            if (kbis.eq.0) then
              write(iout,1000)
              cycle
            end if
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the full energy without any Ewald scaling
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = efull * mscale(kglob)
               if (molcule(iglob) .ne. molcule(kglob))
     &            einter = einter + efull
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0d0 - mscale(kglob)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
               rr11 = bn(5) - scalekk*rr11
c
c     compute the energy contributions for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
c
c     calculate intermediate terms for force and torque
c
               de = term1*rr3 + term2*rr5 + term3*rr7
     &                 + term4*rr9 + term5*rr11
               term1 = -ck*rr3 + drk*rr5 - qrrk*rr7
               term2 = ci*rr3 + dri*rr5 + qrri*rr7
               term3 = 2.0d0 * rr5
               term4 = 2.0d0 * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0d0 * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0d0 * rr7
c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c
c     increment force-based gradient and torque on first site
c
               dem(1,i) = dem(1,i) + frcx
               dem(2,i) = dem(2,i) + frcy
               dem(3,i) = dem(3,i) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kbis) = dem(1,kbis) - frcx
               dem(2,kbis) = dem(2,kbis) - frcy
               dem(3,kbis) = dem(3,kbis) - frcz
               tem(1,kbis) = tem(1,kbis) + ttmk(1)
               tem(2,kbis) = tem(2,kbis) + ttmk(2)
               tem(3,kbis) = tem(3,kbis) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
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
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         call torque (iipole,tem(1,i),fix,fiy,fiz,dem)
         iaz = zaxis(iipole)
         iax = xaxis(iipole)
         iay = yaxis(iipole)
         if (iaz .eq. 0)  iaz = iglob
         if (iax .eq. 0)  iax = iglob
         if (iay .eq. 0)  iay = iglob
         xiz = x(iaz) - x(iglob)
         yiz = y(iaz) - y(iglob)
         ziz = z(iaz) - z(iglob)
         xix = x(iax) - x(iglob)
         yix = y(iax) - y(iglob)
         zix = z(iax) - z(iglob)
         xiy = x(iay) - x(iglob)
         yiy = y(iay) - y(iglob)
         ziy = z(iay) - z(iglob)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
         vxz = zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine emrecip1  --  PME recip multipole energy & derivs  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to multipoles
c
      subroutine emrecip1
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use pme
      use math
      use mpole
      use potent
      use timestat
      use virial
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer i,j,k,ii,iipole,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm,f
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8, allocatable :: cmp(:,:),fmp(:,:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer nprocloc,commloc,rankloc,proc
      real*8 time0,time1
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_BEAD
        rankloc  = rank
      end if
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     indices into the electrostatic field array
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(cphirec)) deallocate (cphirec)
      allocate (cphirec(10,max(npolerecloc,1)))
      cphirec = 0d0
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0d0
c
c     dynamic allocation of local arrays
c
      allocate (cmp(10,npolerecloc))
      allocate (fmp(10,npolerecloc))
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vxy = 0.0d0
      vxz = 0.0d0
      vyy = 0.0d0
      vyz = 0.0d0
      vzz = 0.0d0
c
c     zero out pme grid
c
      qgridin_2d = 0d0
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         cmp(1,i) = rpole(1,iipole)
         cmp(2,i) = rpole(2,iipole)
         cmp(3,i) = rpole(3,iipole)
         cmp(4,i) = rpole(4,iipole)
         cmp(5,i) = rpole(5,iipole)
         cmp(6,i) = rpole(9,iipole)
         cmp(7,i) = rpole(13,iipole)
         cmp(8,i) = 2.0d0 * rpole(6,iipole)
         cmp(9,i) = 2.0d0 * rpole(7,iipole)
         cmp(10,i) = 2.0d0 * rpole(10,iipole)
         call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
         call bspline_fill_site(iglob,i)
      end do
c
c     assign permanent multipoles to PME grid and perform
c     the 3-D FFT forward transformation
c
      time0 = mpi_wtime()
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call grid_mpole_site(iglob,i,fmp(1,i))
      end do
      time1 = mpi_wtime()
      timegrid1 = timegrid1 + time1-time0
c
c     MPI : Begin reception
c
      time0 = mpi_wtime()
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,
     $   commloc,reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'

      do i = 1, nrec_recep
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1) + 
     $    qgridmpi(:,:,:,:,i) 
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1 - time0
c
c     Perform 3-D FFT forward transform
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
c     initialize variables required for the scalar summation
c
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
      time0 = mpi_wtime()
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0d0
      end if
c
c     make the scalar summation over reciprocal lattice
c
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
          do k1 = istart2(rankloc+1),iend2(rankloc+1)
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
            r1 = dble(m1)
            r2 = dble(m2)
            r3 = dble(m3)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0d0
            if (term .gt. -50.0d0) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
     &          + qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-
     $          jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
               eterm = 0.5d0 * electric * expterm * struc2
               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
               vxx = vxx + h1*h1*vterm - eterm
               vxy = vxy + h1*h2*vterm
               vxz = vxz + h1*h3*vterm
               vyy = vyy + h2*h2*vterm - eterm
               vyz = vyz + h2*h3*vterm
               vzz = vzz + h3*h3*vterm - eterm
             qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
           end if
 10        continue
          end do
        end do
      end do
c
c     save the virial for use in polarization computation
c
      vmxx = vxx
      vmxy = vxy
      vmxz = vxz
      vmyy = vyy
      vmyz = vyz
      vmzz = vzz
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5d0 * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e = f * expterm * struc2
           em = em + e
        end if
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, ksize2(rankloc+1)
         do j = 1, jsize2(rankloc+1)
           do i = 1, isize2(rankloc+1)
              term = qfac_2d(i,j,k)
              qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
              qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
            end do
         end do
      end do
      time1 = mpi_wtime()
      timescalar = timescalar + time1-time0
c
c     perform 3-D FFT backward transform and get potential
c
      time0 = mpi_wtime()
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1-time0
c
c     MPI : Begin reception
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,commloc,req2send(i),
     $   ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(req2rec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(req2send(i),status,ierr)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1-time0

      time0 = mpi_wtime()
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        call fphi_mpole_site(iglob,i)
      end do
      time1 = mpi_wtime()
      timegrid2 = timegrid2 + time1-time0
      do i = 1, npolerecloc
         do j = 1, 20
            fphirec(j,i) = electric * fphirec(j,i)
         end do
        call fphi_to_cphi_site (fphirec(1,i),cphirec(1,i))
      end do
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphirec(k,i)
            f1 = f1 + fmp(k,i)*fphirec(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphirec(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphirec(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         ii = locrec1(iglob)
         demrec(1,ii) = demrec(1,ii) + h1
         demrec(2,ii) = demrec(2,ii) + h2
         demrec(3,ii) = demrec(3,ii) + h3
      end do
      e = 0.5d0 * e
      em = em + e
c
c     distribute torques into the permanent multipole gradient
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         trq(1) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
         trq(2) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
         trq(3) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
         call torque_rec (iipole,trq,fix,fiy,fiz,demrec)
      end do
c
c     permanent multipole contribution to the internal virial
c
      do i = 1, npolerecloc
         vxx = vxx - cmp(2,i)*cphirec(2,i) - 2.0d0*cmp(5,i)*cphirec(5,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(9,i)*cphirec(9,i)
         vxy = vxy - 0.5d0*(cmp(3,i)*cphirec(2,i)+cmp(2,i)*cphirec(3,i))
     &            -(cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &            -0.5d0*cmp(8,i)*(cphirec(5,i)+cphirec(6,i))
     &            -0.5d0*(cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
         vxz = vxz - 0.5d0*(cmp(4,i)*cphirec(2,i)+cmp(2,i)*cphirec(4,i))
     &            -(cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &            -0.5d0*cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &            -0.5d0*(cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
         vyy = vyy - cmp(3,i)*cphirec(3,i) - 2.0d0*cmp(6,i)*cphirec(6,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(10,i)*cphirec(10,i)
         vyz = vyz - 0.5d0*(cmp(4,i)*cphirec(3,i)+cmp(3,i)*cphirec(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &            - 0.5d0*cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &            - 0.5d0*(cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
         vzz = vzz - cmp(4,i)*cphirec(4,i) - 2.0d0*cmp(7,i)*cphirec(7,i)
     &            - cmp(9,i)*cphirec(9,i) - cmp(10,i)*cphirec(10,i)
      end do
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      return
      end

c
c     "dercorrecemtp" computes the derivatives of the correction term that has to be added to the expression
c     of the electrostatic energy to include charge penetration effect
c
      subroutine dercorrecemtp(cutcorrec)
      use atmtyp
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use deriv
      use domdec
      use mpole
      use mplpot
      use polgrp
      use polpot
      use potent
      use shunt
      use usage
      use virial
      use neigh
      implicit none
      integer i,ii,iz,ix,iy,j,kkk,inl
      integer izloc,ixloc,iyloc
      integer iipole,iglob,kkpole,kglob,kbis
      integer k,kk,kz,kx,ky,l,m,n1,n2,l3
      integer kzloc,kxloc,kyloc
      real*8  cutcorrec,xsize,ysize,zsize,xsize2,ysize2,zsize2
      real*8 xmove,ymove,zmove
      real*8  f,fgrp,ci,dix,diy,diz
      real*8  ck,dkx,dky,dkz
      real*8  xr,yr,zr,r2,r,r3,r4,rr1,rr3,rr5,rr7
      real*8  alpha1,alpha2,beta1,beta2,rvdw1,rvdw2,etaemtp
      real*8  eta1,eta2
      real*8  sc,ecorr1,ecorr2,fm,e,ecorrtot
      real*8  sc1,sc2,sc11,sc22,scqir,scqkr
      real*8  ecorr11,ecorr12,ecorr13,ecorr14,ecorr15,ecorr21
      real*8  ecorr3,ecorr31,ecorr32,de31(3),de32(3),phi
      real*8  decorr21, d(3), di(3,3), dx(3), qit(3,3,3)
      real*8  qir(3),qkr(3),q(3,3),qr(3,3)
      real*8  rot(3,3),dri(3,3,3),driz(3,3,3),drix(3,3,3),driy(3,3,3)
      real*8  drk(3,3,3),drkz(3,3,3),drkx(3,3,3),drky(3,3,3)
      real*8  de(3),de11(3),de12(3)
      real*8  ftm(3),trqi(3),trqix(3),trqiy(3),trqiz(3)
      real*8  trqk(3),trqkx(3),trqky(3),trqkz(3)
      real*8  s,ds,rin,rout
      real*8  vxx,vyx,vzx,vyy,vzy,vzz
      real*8  xiz,yiz,ziz,xix,yix,zix,xiy,yiy,ziy
      real*8  xkz,ykz,zkz,xkx,ykx,zkx,xky,yky,zky
      real*8  qixx,qixy,qixz
      real*8  qiyy,qiyz,qizz
      real*8  qkxx,qkxy,qkxz
      real*8  qkyy,qkyz,qkzz
      real*8  qix,qiy,qiz
      real*8  qkx,qky,qkz
      real*8, allocatable :: mscale(:)
      logical usei,usek,proceed
      character*6 mode
c
      e = 0.0d0
      ecorrtot = 0.0d0
      mode = 'MPOLE'
      call switch (mode)
      rin = 1.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      if (npole .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         inl = polelocnl(iipole)
         iglob = glob(i)
         iz = zaxis(iipole)
         if (iz.ne.0) izloc = loc(iz)
         ix = xaxis(iipole)
         if (ix.ne.0) ixloc = loc(ix)
         iy = yaxis(iipole)
         if (iy.ne.0) iyloc = loc(iy)
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
c
c     set interaction scaling coefficients for connected atoms
c
         f = electric / dielec
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
         do kkk = 1, nelst(inl)
           kkpole = elst(kkk,inl)
           kglob = ipole(kkpole)
           kbis = loc(kglob)
           kz = zaxis(kkpole)
           if (kz.ne.0) kzloc = loc(kz)
           kx = xaxis(kkpole)
           if (kx.ne.0) kxloc = loc(kx)
           ky = yaxis(kkpole)
           if (ky.ne.0) kyloc = loc(ky)
           e  = 0.0d0
           ecorr1  = 0.0d0
           ecorr2  = 0.0d0
           vxx = 0.0d0
           vyx = 0.0d0
           vzx = 0.0d0
           vyy = 0.0d0
           vzy = 0.0d0
           vzz = 0.0d0
c
c     compute the energy correction for this interaction
c
           call aclear(3,de)
           xr = x(kglob) - x(iglob)
           yr = y(kglob) - y(iglob)
           zr = z(kglob) - z(iglob)
           call image (xr,yr,zr)
           dx(1) = xr
           dx(2) = yr
           dx(3) = zr
           r2 = xr*xr + yr*yr + zr*zr
           if (r2 .le. (cutcorrec)*(cutcorrec)) then
              r = sqrt(r2)
              rr1 = 1d0/r
              r3 = r2*r
              rr5 = 3.0d0 / (r2*r3)
              rr7 = 5.0d0 * rr5 / r2
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
c     construct some intermediate quadrupole values
c
              qix = qixx*xr + qixy*yr + qixz*zr
              qiy = qixy*xr + qiyy*yr + qiyz*zr
              qiz = qixz*xr + qiyz*yr + qizz*zr
              qkx = qkxx*xr + qkxy*yr + qkxz*zr
              qky = qkxy*xr + qkyy*yr + qkyz*zr
              qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
              qir(1) = qixx*xr + qixy*yr + qixz*zr
              qir(2) = qixy*xr + qiyy*yr + qiyz*zr
              qir(3) = qixz*xr + qiyz*yr + qizz*zr
              scqir  = qir(1)*xr + qir(2)*yr + qir(3)*zr
              qkr(1) = qkxx*xr + qkxy*yr + qkxz*zr
              qkr(2) = qkxy*xr + qkyy*yr + qkyz*zr
              qkr(3) = qkxz*xr + qkyz*yr + qkzz*zr
              scqkr = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
              if (use_emtp) then
                alpha1 = alphapen(iipole)
                beta1  = betapen(iipole)
                eta1   = gammapen(iipole)
c
                alpha2 = alphapen(kkpole)
                beta2  = betapen(kkpole)
                eta2   = gammapen(kkpole)
              else if (use_emtporig) then
                alpha1 = gamma1pen/vdwpen(iglob)
                beta1  = deltapen/vdwpen(iglob)
                eta1   = khipen*2/(vdwpen(iglob)+vdwpen(kglob))
c
                alpha2 = gamma1pen/vdwpen(kglob)
                beta2  = deltapen/vdwpen(kglob)
                eta2   = eta1
              end if
c
c    charge-charge correction term
c
              ecorr11 = valemtp(atomic(iglob))*(valemtp(atomic(kglob))
     &  -ck)*exp(-alpha2*r)
              ecorr12 = valemtp(atomic(kglob))*(valemtp(atomic(iglob))
     &  -ci)*exp(-alpha1*r)
              ecorr13 = (valemtp(atomic(iglob))-ci)*
     &  (valemtp(atomic(kglob))-ck)*exp(-(beta1+beta2)*r)
              ecorr14 = -(valemtp(atomic(iglob))-ci)*
     &  (valemtp(atomic(kglob))-ck)*exp(-beta1*r)
              ecorr15 = -(valemtp(atomic(iglob))-ci)*
     &  (valemtp(atomic(kglob))-ck)*exp(-beta2*r)
              ecorr1 = (ecorr11+ecorr12+ecorr13+ecorr14+ecorr15)/r
c
              de(1) = xr*(alpha2*ecorr11 +
     $         alpha1*ecorr12 + (beta1+beta2)*ecorr13 +
     $         beta1*ecorr14 + beta2*ecorr15)/r2 +
     $         xr*(ecorr11 + ecorr12 + ecorr13 + ecorr14 +
     $         ecorr15)/r3
              de(2) = yr*(alpha2*ecorr11 +
     $         alpha1*ecorr12 + (beta1+beta2)*ecorr13 +
     $         beta1*ecorr14 + beta2*ecorr15)/r2 +
     $         yr*(ecorr11 + ecorr12 + ecorr13 + ecorr14 +
     $         ecorr15)/r3
              de(3) = zr*(alpha2*ecorr11 +
     $         alpha1*ecorr12 + (beta1+beta2)*ecorr13 +
     $         beta1*ecorr14 + beta2*ecorr15)/r2 +
     $         zr*(ecorr11 + ecorr12 + ecorr13 + ecorr14 +
     $         ecorr15)/r3
c
c    charge-permanent dipole correction term
c
              sc11 = dkx*xr + dky*yr + dkz*zr
              sc1 = exp(-eta1*r)*(valemtp(atomic(iglob))-ci)
              sc22 = dix*xr + diy*yr + diz*zr
              sc2 = exp(-eta2*r)*(valemtp(atomic(kglob))-ck)
              ecorr2 = -(sc1*sc11-sc2*sc22)/r3
c
c    get switchung function terms
c
              call switch_emtp(.true.,r,cutcorrec-rin,cutcorrec,s,ds)
c
              fm = f*mscale(kglob)
              e = s*(ecorr1 + ecorr2)
c
              de11(1) = -(sc1/r3)*(sc11*eta1*xr/r +
     $        3*sc11*xr/r2 - dkx)
              de11(2) = -(sc1/r3)*(sc11*eta1*yr/r +
     $        3*sc11*yr/r2 - dky)
              de11(3) = -(sc1/r3)*(sc11*eta1*zr/r +
     $        3*sc11*zr/r2 - dkz)
              de12(1) = (sc2/r3)*(sc22*eta2*xr/r +
     $        3*sc22*xr/r2 - dix)
              de12(2) = (sc2/r3)*(sc22*eta2*yr/r +
     $        3*sc22*yr/r2 - diy)
              de12(3) = (sc2/r3)*(sc22*eta2*zr/r +
     $        3*sc22*zr/r2 - diz)
c
              de31 = 0d0
              de32 = 0d0
              ecorr3 = 0d0
              if (use_emtporig) then
c
c    charge-quadrupole correction term
c
              phi = omegapen*2/(vdwpen(iglob)+vdwpen(kglob))
              ecorr31 = (valemtp(atomic(iglob))-ci)*exp(-phi*r)
     &          *scqkr*rr5
              ecorr32 = (valemtp(atomic(kglob))-ck)*exp(-phi*r)
     &          *scqir*rr5
              ecorr3  = ecorr31 + ecorr32
              e = e + s*ecorr3
              de31(1) = ecorr31*(phi*xr*rr1)-2*qkr(1)*
     &        (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
     &      + (valemtp(atomic(iglob))-ci)*exp(-phi*r)*scqkr*rr7*xr
              de31(2) = ecorr31*(phi*yr*rr1)-2*qkr(2)*
     &        (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
     &      + (valemtp(atomic(iglob))-ci)*exp(-phi*r)*scqkr*rr7*yr
              de31(3) = ecorr31*(phi*zr*rr1)-2*qkr(3)*
     &        (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
     &      + (valemtp(atomic(iglob))-ci)*exp(-phi*r)*scqkr*rr7*zr

              de32(1) = ecorr32*(phi*xr*rr1)-2*qir(1)*
     &        (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
     &      + (valemtp(atomic(kglob))-ck)*exp(-phi*r)*scqir*rr7*xr
              de32(2) = ecorr32*(phi*yr*rr1)-2*qir(2)*
     &        (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
     &      + (valemtp(atomic(kglob))-ck)*exp(-phi*r)*scqir*rr7*yr
              de32(3) = ecorr32*(phi*zr*rr1)-2*qir(3)*
     &        (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
     &      + (valemtp(atomic(kglob))-ck)*exp(-phi*r)*scqir*rr7*zr
              end if
c
              ftm(1) = fm*s*
     $        (de(1)+de11(1)+de12(1)+de31(1)+de32(1)) - 
     $        fm*xr*ds*(ecorr1+ecorr2+ecorr3)/r
              ftm(2) = fm*s*
     $        (de(2)+de11(2)+de12(2)+de31(2)+de32(2)) - 
     $        fm*yr*ds*(ecorr1+ecorr2+ecorr3)/r
              ftm(3) = fm*s*
     $        (de(3)+de11(3)+de12(3)+de31(3)+de32(3)) - 
     $        fm*zr*ds*(ecorr1+ecorr2+ecorr3)/r
c
              dem(1,i) = dem(1,i) + ftm(1)
              dem(2,i) = dem(2,i) + ftm(2)
              dem(3,i) = dem(3,i) + ftm(3)
c
              vxx = vxx - xr*ftm(1)
              vyx = vyx - yr*ftm(1)
              vzx = vzx - zr*ftm(1)
              vyy = vyy - yr*ftm(2)
              vzy = vzy - zr*ftm(2)
              vzz = vzz - zr*ftm(3)
c
              dem(1,kbis) = dem(1,kbis) - ftm(1)
              dem(2,kbis) = dem(2,kbis) - ftm(2)
              dem(3,kbis) = dem(3,kbis) - ftm(3)
c
c    compute the contribution of the derivatives of the rotation matrixes
c
              trqi = 0d0
              trqix = 0d0
              trqiy = 0d0
              trqiz = 0d0
              trqk = 0d0
              trqkx = 0d0
              trqky = 0d0
              trqkz = 0d0
c
               call derrot(iipole,.true.,iglob,iz,ix,iy,rot,dri,driz,
     $          drix,driy)
c
               call aclear(3,d)
               d(1) = pole(2,iipole)
               d(2) = pole(3,iipole)
               d(3) = pole(4,iipole)
               q(1,1)  = pole(5,iipole)
               q(1,2)  = pole(6,iipole)
               q(1,3)  = pole(7,iipole)
               q(2,1)  = pole(8,iipole) 
               q(2,2)  = pole(9,iipole) 
               q(2,3)  = pole(10,iipole) 
               q(3,1)  = pole(11,iipole) 
               q(3,2)  = pole(12,iipole) 
               q(3,3)  = pole(13,iipole) 
               qr(1,1) = rpole(5,iipole)
               qr(1,2) = rpole(6,iipole)
               qr(1,3) = rpole(7,iipole)
               qr(2,1) = rpole(8,iipole) 
               qr(2,2) = rpole(9,iipole) 
               qr(2,3) = rpole(10,iipole) 
               qr(3,1) = rpole(11,iipole) 
               qr(3,2) = rpole(12,iipole) 
               qr(3,3) = rpole(13,iipole) 
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + dri(1,l,m)*d(m)
                   di(2,l) = di(2,l) + dri(2,l,m)*d(m)
                   di(3,l) = di(3,l) + dri(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(dri(l3,l,n1)
     $                         *rot(m,n2)
     $                         + dri(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               trqi(1) = fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $          exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
               trqi(2) = fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $          exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
               trqi(3) = fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $          exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
c
               if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqi(1) = trqi(1) + fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqi(2) = trqi(2) + fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqi(3) = trqi(3) + fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                      end do
                 end do
               end if
c
               dem(1,i) = dem(1,i) + trqi(1)
               dem(2,i) = dem(2,i) + trqi(2)
               dem(3,i) = dem(3,i) + trqi(3)
c
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + driz(1,l,m)*d(m)
                   di(2,l) = di(2,l) + driz(2,l,m)*d(m)
                   di(3,l) = di(3,l) + driz(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(driz(l3,l,n1)
     $                         *rot(m,n2)
     $                         + driz(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               if (iz.ne.0) then
                 trqiz(1) =
     $            fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
                 trqiz(2) =
     $            fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
                 trqiz(3) =
     $            fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
c
                  if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqiz(1) =trqiz(1) +fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqiz(2) =trqiz(2) +fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqiz(3) =trqiz(3) +fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                      end do
                    end do
                 end if
                 xiz = x(iz) - x(iglob)
                 yiz = y(iz) - y(iglob)
                 ziz = z(iz) - z(iglob)
                 vxx = vxx + xiz*trqiz(1)
                 vyx = vyx + yiz*trqiz(1)
                 vzx = vzx + ziz*trqiz(1)
                 vyy = vyy + yiz*trqiz(2)
                 vzy = vzy + ziz*trqiz(2)
                 vzz = vzz + ziz*trqiz(3)
                 dem(1,izloc) = dem(1,izloc) + trqiz(1)
                 dem(2,izloc) = dem(2,izloc) + trqiz(2)
                 dem(3,izloc) = dem(3,izloc) + trqiz(3)
               end if
c
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + drix(1,l,m)*d(m)
                   di(2,l) = di(2,l) + drix(2,l,m)*d(m)
                   di(3,l) = di(3,l) + drix(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(drix(l3,l,n1)
     $                         *rot(m,n2)
     $                         + drix(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               if (ix.ne.0) then
                 trqix(1) =
     $            fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
                 trqix(2) =
     $            fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
                 trqix(3) =
     $            fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
c
                  if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqix(1) = trqix(1)+fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqix(2) = trqix(2)+fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqix(3) = trqix(3)+fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                      end do
                    end do
                  end if
                 xix = x(ix) - x(iglob)
                 yix = y(ix) - y(iglob)
                 zix = z(ix) - z(iglob)
                 vxx = vxx + xix*trqix(1)
                 vyx = vyx + yix*trqix(1)
                 vzx = vzx + zix*trqix(1)
                 vyy = vyy + yix*trqix(2)
                 vzy = vzy + zix*trqix(2)
                 vzz = vzz + zix*trqix(3)
                 dem(1,ixloc) = dem(1,ixloc) + trqix(1)
                 dem(2,ixloc) = dem(2,ixloc) + trqix(2)
                 dem(3,ixloc) = dem(3,ixloc) + trqix(3)
               end if
c
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + driy(1,l,m)*d(m)
                   di(2,l) = di(2,l) + driy(2,l,m)*d(m)
                   di(3,l) = di(3,l) + driy(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(driy(l3,l,n1)
     $                         *rot(m,n2)
     $                         + driy(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               if (iy.ne.0) then
                 trqiy(1) =
     $            fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
                 trqiy(2) =
     $            fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
                 trqiy(3) =
     $            fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $            exp(-eta2*r)*(valemtp(atomic(kglob))-ck)/r3
c
                    if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqiy(1) = trqiy(1)+fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqiy(2) = trqiy(2)+fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                        trqiy(3) = trqiy(3)+fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(kglob))-ck)*exp(-phi*r)*rr5
                      end do
                    end do
                  end if
                 xiy = x(iy) - x(iglob)
                 yiy = y(iy) - y(iglob)
                 ziy = z(iy) - z(iglob)
                 vxx = vxx + xiy*trqiy(1)
                 vyx = vyx + yiy*trqiy(1)
                 vzx = vzx + ziy*trqiy(1)
                 vyy = vyy + yiy*trqiy(2)
                 vzy = vzy + ziy*trqiy(2)
                 vzz = vzz + ziy*trqiy(3)
                 dem(1,iyloc) = dem(1,iyloc) + trqiy(1)
                 dem(2,iyloc) = dem(2,iyloc) + trqiy(2)
                 dem(3,iyloc) = dem(3,iyloc) + trqiy(3)
               end if
c

               call derrot(kkpole,.true.,kglob,kz,kx,ky,rot,drk,drkz,
     $          drkx,drky)
c
               call aclear(3,d)
               d(1) = pole(2,kkpole)
               d(2) = pole(3,kkpole)
               d(3) = pole(4,kkpole)
               q(1,1)  = pole(5,kkpole)
               q(1,2)  = pole(6,kkpole)
               q(1,3)  = pole(7,kkpole)
               q(2,1)  = pole(8,kkpole) 
               q(2,2)  = pole(9,kkpole) 
               q(2,3)  = pole(10,kkpole) 
               q(3,1)  = pole(11,kkpole) 
               q(3,2)  = pole(12,kkpole) 
               q(3,3)  = pole(13,kkpole) 
               qr(1,1) = rpole(5,kkpole)
               qr(1,2) = rpole(6,kkpole)
               qr(1,3) = rpole(7,kkpole)
               qr(2,1) = rpole(8,kkpole) 
               qr(2,2) = rpole(9,kkpole) 
               qr(2,3) = rpole(10,kkpole) 
               qr(3,1) = rpole(11,kkpole) 
               qr(3,2) = rpole(12,kkpole) 
               qr(3,3) = rpole(13,kkpole) 
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + drk(1,l,m)*d(m)
                   di(2,l) = di(2,l) + drk(2,l,m)*d(m)
                   di(3,l) = di(3,l) + drk(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(drk(l3,l,n1)
     $                         *rot(m,n2)
     $                         + drk(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               trqk(1) =
     $          -fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $          exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
               trqk(2) =
     $          -fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $          exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
               trqk(3) =
     $          -fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $          exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
c
                 if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqk(1) = trqk(1)+fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqk(2) = trqk(2)+fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqk(3) = trqk(3)+fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                      end do
                    end do
                end if
c
               dem(1,kbis) = dem(1,kbis) + trqk(1)
               dem(2,kbis) = dem(2,kbis) + trqk(2)
               dem(3,kbis) = dem(3,kbis) + trqk(3)
c
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + drkz(1,l,m)*d(m)
                   di(2,l) = di(2,l) + drkz(2,l,m)*d(m)
                   di(3,l) = di(3,l) + drkz(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(drkz(l3,l,n1)
     $                         *rot(m,n2)
     $                         + drkz(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               if (kz.ne.0) then
                trqkz(1) =
     $             -fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $             exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
                trqkz(2) =
     $             -fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $             exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
                trqkz(3) =
     $            -fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $            exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
c
                if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqkz(1) = trqkz(1)+fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqkz(2) = trqkz(2)+fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqkz(3) = trqkz(3)+fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                      end do
                    end do
                 end if
                 xkz = x(kz) - x(kglob)
                 ykz = y(kz) - y(kglob)
                 zkz = z(kz) - z(kglob)
                 vxx = vxx + xkz*trqkz(1)
                 vyx = vyx + ykz*trqkz(1)
                 vzx = vzx + zkz*trqkz(1)
                 vyy = vyy + ykz*trqkz(2)
                 vzy = vzy + zkz*trqkz(2)
                 vzz = vzz + zkz*trqkz(3)
                 dem(1,kzloc) = dem(1,kzloc) + trqkz(1)
                 dem(2,kzloc) = dem(2,kzloc) + trqkz(2)
                 dem(3,kzloc) = dem(3,kzloc) + trqkz(3)
               end if
c
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + drkx(1,l,m)*d(m)
                   di(2,l) = di(2,l) + drkx(2,l,m)*d(m)
                   di(3,l) = di(3,l) + drkx(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(drkx(l3,l,n1)
     $                         *rot(m,n2)
     $                         + drkx(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               if (kx.ne.0) then
                 trqkx(1) =
     $            -fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $            exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
                 trqkx(2) =
     $            -fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $            exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
                 trqkx(3) =
     $            -fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $            exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
c
                 if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqkx(1) = trqkx(1)+fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqkx(2) = trqkx(2)+fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqkx(3) = trqkx(3)+fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                      end do
                    end do
                 end if
                 xkx = x(kx) - x(kglob)
                 ykx = y(kx) - y(kglob)
                 zkx = z(kx) - z(kglob)
                 vxx = vxx + xkx*trqkx(1)
                 vyx = vyx + ykx*trqkx(1)
                 vzx = vzx + zkx*trqkx(1)
                 vyy = vyy + ykx*trqkx(2)
                 vzy = vzy + zkx*trqkx(2)
                 vzz = vzz + zkx*trqkx(3)
                 dem(1,kxloc) = dem(1,kxloc) + trqkx(1)
                 dem(2,kxloc) = dem(2,kxloc) + trqkx(2)
                 dem(3,kxloc) = dem(3,kxloc) + trqkx(3)
               end if
c
               call aclear(9,di)
               call aclear(27,qit)
               do l = 1, 3
                 do m = 1, 3
                   di(1,l) = di(1,l) + drky(1,l,m)*d(m)
                   di(2,l) = di(2,l) + drky(2,l,m)*d(m)
                   di(3,l) = di(3,l) + drky(3,l,m)*d(m)
                   do n1 = 1,3
                      do n2 = 1,3
                         do l3 = 1,3
                            qit(l3,l,m) = qit(l3,l,m)
     $                         + q(n1,n2)*(drky(l3,l,n1)
     $                         *rot(m,n2)
     $                         + drky(l3,m,n2)*rot(l,n1))
                         end do
                      end do
                   end do
                 end do
               end do
               if (ky.ne.0) then
                 trqky(1) =
     $              -fm*s*(di(1,1)*xr + di(1,2)*yr + di(1,3)*zr)*
     $              exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
                 trqky(2) =
     $              -fm*s*(di(2,1)*xr + di(2,2)*yr + di(2,3)*zr)*
     $              exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
                 trqky(3) =
     $              -fm*s*(di(3,1)*xr + di(3,2)*yr + di(3,3)*zr)*
     $              exp(-eta1*r)*(valemtp(atomic(iglob))-ci)/r3
c
                  if (use_emtporig) then
c
c    contribution due to the quadrupoles
c
                    do l = 1, 3
                      do m = 1, 3
                        trqky(1) = trqky(1)+fm*s*qit(1,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqky(2) = trqky(2)+fm*s*qit(2,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                        trqky(3) = trqky(3)+fm*s*qit(3,l,m)*dx(l)*dx(m)*
     $                (valemtp(atomic(iglob))-ci)*exp(-phi*r)*rr5
                      end do
                    end do
                  end if
                 xky = x(ky) - x(kglob)
                 yky = y(ky) - y(kglob)
                 zky = z(ky) - z(kglob)
                 vxx = vxx + xky*trqky(1)
                 vyx = vyx + yky*trqky(1)
                 vzx = vzx + zky*trqky(1)
                 vyy = vyy + yky*trqky(2)
                 vzy = vzy + zky*trqky(2)
                 vzz = vzz + zky*trqky(3)
                 dem(1,ky) = dem(1,ky) + trqky(1)
                 dem(2,ky) = dem(2,ky) + trqky(2)
                 dem(3,ky) = dem(3,ky) + trqky(3)
               end if
c
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
c
c
c     apply the energy adjustments for scaled interactions
c
               e = fm * e
c
c     increment the total electrostatic energy
c
               em = em + e
               ecorrtot = ecorrtot + e
c
          end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
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
