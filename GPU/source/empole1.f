c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      subroutine empole1
      use domdec
      use energi
      use potent
      use tinheader ,only:ti_p,re_p
      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empole1c
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0_re_p
      end if
      if (.not. use_polar) then
         ep = 0.0_re_p
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
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use timestat
      use tinheader ,only:ti_p,re_p
      use virial
      use mpi
      implicit none
      integer i,j,ii
      integer iipole,iglob,ierr
      real(t_p) e,f
      real(t_p) term,fterm
      real(t_p) cii,dii,qii
      real(t_p) xd,yd,zd
      real(t_p) xq,yq,zq
      real(t_p) xv,yv,zv,vterm
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) xdfield,ydfield
      real(t_p) zdfield
      real(t_p) trq(3),frcx(3)
      real(t_p) frcy(3),frcz(3)
      real(t_p) time0,time1
c
c
c     zero out the atomic multipole energy and derivatives
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'empole1c'
      em = 0.0_re_p
      dem = 0_re_p
      if (npole .eq. 0)  return
!$acc update host(vir,dem)
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
     $   then
        if (use_mrec) then 
        call timer_enter( timer_rec )
        call emrecip1
        call timer_exit( timer_rec,quiet_timers )
        end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        call timer_enter( timer_real )
        if (use_mreal) then
          if (use_mpoleshortreal) then
c            write(*,*) 'gogo shortreal'
            call emrealshort1c
          else if (use_mpolelong) then
c            write(*,*) 'gogo longreal'
            call emreallong1c
          else
c          write(*,*) 'gogo real'
            call emreal1c
          end if
        end if
c
        if (use_mself) then
c          write(*,*) 'gogo emself'
c
c     compute the Ewald self-energy term over all the atoms
c
        term = 2.0_ti_p * aewald * aewald
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
           qii = 2.0_ti_p*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &              + qixx*qixx + qiyy*qiyy + qizz*qizz
           e = fterm * (cii + term*(dii/3.0_ti_p+
     &                             2.0_ti_p*term*qii/5.0_ti_p))
           em = em + e
        end do
c
c       compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0_ti_p
           yd = 0.0_ti_p
           zd = 0.0_ti_p
           do i = 1, npoleloc
              iipole = poleglob(i)
              iglob = ipole(iipole)
              xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
              yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
              zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
           end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
             if (rank.eq.0) then
               em = em + term*(xd*xd+yd*yd+zd*zd)
             end if
           do ii = 1, npoleloc
              iipole = poleglob(ii)
              iglob = ipole(iipole)
              i = loc(iglob)
              dem(1,i) = dem(1,i) + 2.0_ti_p*term*rpole(1,iipole)*xd
              dem(2,i) = dem(2,i) + 2.0_ti_p*term*rpole(1,iipole)*yd
              dem(3,i) = dem(3,i) + 2.0_ti_p*term*rpole(1,iipole)*zd
           end do
           xdfield = -2.0_ti_p * term * xd
           ydfield = -2.0_ti_p * term * yd
           zdfield = -2.0_ti_p * term * zd
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
           xd = 0.0_ti_p
           yd = 0.0_ti_p
           zd = 0.0_ti_p
           xq = 0.0_ti_p
           yq = 0.0_ti_p
           zq = 0.0_ti_p
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
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xq,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yq,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zq,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               xv = xd * xq
               yv = yd * yq
               zv = zd * zq
               vterm = term*(xd*xd+ yd*yd+ zd*zd+ 2.0_ti_p*(xv+yv+zv)
     &                            + xq*xq + yq*yq + zq*zq)
               vir(1,1) = vir(1,1) + 2.0_ti_p*term*(xq*xq+xv) + vterm
               vir(2,1) = vir(2,1) + 2.0_ti_p*term*(xq*yq+xv)
               vir(3,1) = vir(3,1) + 2.0_ti_p*term*(xq*zq+xv)
               vir(1,2) = vir(1,2) + 2.0_ti_p*term*(yq*xq+yv)
               vir(2,2) = vir(2,2) + 2.0_ti_p*term*(yq*yq+yv) + vterm
               vir(3,2) = vir(3,2) + 2.0_ti_p*term*(yq*zq+yv)
               vir(1,3) = vir(1,3) + 2.0_ti_p*term*(zq*xq+zv)
               vir(2,3) = vir(2,3) + 2.0_ti_p*term*(zq*yq+zv)
               vir(3,3) = vir(3,3) + 2.0_ti_p*term*(zq*zq+zv) + vterm
             end if
          end if
        end if
        call timer_exit( timer_real,quiet_timers )
      end if
!$acc update device(vir,dem)
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
      use shunt
      use tinheader ,only:ti_p,re_p
      use virial
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis
      integer ii,kkk,iipole,kkpole

      integer iax,iay,iaz
      real(t_p) e,efull,de,f
      real(t_p) bfac,erfc
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) scalekk
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) xix,yix,zix
      real(t_p) xiy,yiy,ziy
      real(t_p) xiz,yiz,ziz
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9,rr11
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) dikx,diky,dikz
      real(t_p) dirx,diry,dirz
      real(t_p) dkrx,dkry,dkrz
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrixr,qriyr,qrizr
      real(t_p) qrkxr,qrkyr,qrkzr
      real(t_p) qrrx,qrry,qrrz
      real(t_p) qikrx,qikry,qikrz
      real(t_p) qkirx,qkiry,qkirz
      real(t_p) qikrxr,qikryr,qikrzr
      real(t_p) qkirxr,qkiryr,qkirzr
      real(t_p) diqkx,diqky,diqkz
      real(t_p) dkqix,dkqiy,dkqiz
      real(t_p) diqkxr,diqkyr,diqkzr
      real(t_p) dkqixr,dkqiyr,dkqizr
      real(t_p) dqiqkx,dqiqky,dqiqkz
      real(t_p) dri,drk,qrri,qrrk
      real(t_p) diqrk,dkqri
      real(t_p) dik,qik,qrrik
      real(t_p) term1,term2,term3
      real(t_p) term4,term5,term6
      real(t_p) frcx,frcy,frcz
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) ttmi(3),ttmk(3)
      real(t_p) fix(3),fiy(3),fiz(3)
      real(t_p) bn(0:5)
      real(t_p), allocatable :: mscale(:)
      real(t_p), allocatable :: tem(:,:)
      character*10 mode

 1000 format(' Warning, system moved too much since last neighbor list',
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
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'emreal1c'
      mscale = 1.0_ti_p
      tem = 0.0_ti_p
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
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
               rr9 = 7.0_ti_p * rr7 / r2
               rr11 = 9.0_ti_p * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0_ti_p * aewald**2
               alsq2n = 0.0_ti_p
               if (aewald .gt. 0.0_ti_p)
     &            alsq2n = 1.0_ti_p / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = real(j+j-1,t_p)
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
               qik = 2.0_ti_p*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
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
     &                     - 2.0_ti_p*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0_ti_p*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0_ti_p*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0_ti_p*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0_ti_p*qrrik
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
               scalekk = 1.0_ti_p - mscale(kglob)
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
               term3 = 2.0_ti_p * rr5
               term4 = 2.0_ti_p * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0_ti_p * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0_ti_p * rr7
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
            mscale(i12(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0_ti_p
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
         if (iaz .le. 0)  iaz = iglob
         if (iax .le. 0)  iax = iglob
         if (iay .le. 0)  iay = iglob
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
c     OpenMP directives for the major loop structure
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
      use tinheader ,only:ti_p,re_p
      use virial
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer i,j,k,ii,iipole,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real(t_p) zero,one,two,half
      real(t_p) e,eterm,f
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) f1,f2,f3
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) vterm,struc2
      real(t_p) trq(3),fix(3)
      real(t_p) fiy(3),fiz(3)
      real(t_p), allocatable :: cmp(:,:),fmp(:,:)
      real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer nprocloc,commloc,rankloc,proc
      real(t_p) time0,time1

      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p, half=0.5_ti_p)

      parameter( ! indices into the electrostatic field array
     &  deriv1=(/ 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /),
     &  deriv2=(/ 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /),
     &  deriv3=(/ 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /))
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'emrecip1'

      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
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
      cphirec = 0_ti_p
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0_ti_p
c
c     dynamic allocation of local arrays
c
      allocate (cmp(10,npolerecloc))
      allocate (fmp(10,npolerecloc))
c
c     zero out the temporary virial accumulation variables
c
      vxx = zero
      vxy = zero
      vxz = zero
      vyy = zero
      vyz = zero
      vzz = zero
c
c     zero out pme grid
c
      qgridin_2d = 0_ti_p
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
         cmp(8,i) = two * rpole(6,iipole)
         cmp(9,i) = two * rpole(7,iipole)
         cmp(10,i) = two * rpole(10,iipole)
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
     $   n3mpimax,MPI_TPREC,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,proc,tag,
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
           qfac_2d(1,1,1) = zero
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
            r1 = real(m1,t_p)
            r2 = real(m2,t_p)
            r3 = real(m3,t_p)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = zero
            if (term .gt. -50.0_ti_p) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (one-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = zero
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
     &          + qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-
     $          jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2
               eterm = half * electric * expterm * struc2
               vterm = (two/hsq) * (one-term) * eterm
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
           expterm = half * pi / xbox
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
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_TPREC,prec_recep(i),tag,commloc,req2send(i),
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
      e = zero
      do i = 1, npolerecloc
         f1 = zero
         f2 = zero
         f3 = zero
         do k = 1, 10
            e = e + fmp(k,i)*fphirec(k,i)
            f1 = f1 + fmp(k,i)*fphirec(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphirec(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphirec(deriv3(k),i)
         end do
         f1 = real(nfft1,t_p) * f1
         f2 = real(nfft2,t_p) * f2
         f3 = real(nfft3,t_p) * f3
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
      e = half * e
      em = em + e
c
c     distribute torques into the permanent multipole gradient
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         trq(1) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
     &               + two*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
         trq(2) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
     &               + two*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
         trq(3) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
     &               + two*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
         call torque_rec (iipole,trq,fix,fiy,fiz,demrec)
      end do
c
c     permanent multipole contribution to the internal virial
c
      do i = 1, npolerecloc
         vxx = vxx - cmp(2,i)*cphirec(2,i) - two*cmp(5,i)*cphirec(5,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(9,i)*cphirec(9,i)
         vxy = vxy - half*(cmp(3,i)*cphirec(2,i)+cmp(2,i)*cphirec(3,i))
     &            -(cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &            -half*cmp(8,i)*(cphirec(5,i)+cphirec(6,i))
     &            -half*(cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
         vxz = vxz - half*(cmp(4,i)*cphirec(2,i)+cmp(2,i)*cphirec(4,i))
     &            -(cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &            -half*cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &            -half*(cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
         vyy = vyy - cmp(3,i)*cphirec(3,i) - two*cmp(6,i)*cphirec(6,i)
     &            - cmp(8,i)*cphirec(8,i) - cmp(10,i)*cphirec(10,i)
         vyz = vyz - half*(cmp(4,i)*cphirec(3,i)+cmp(3,i)*cphirec(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &            - half*cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &            - half*(cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
         vzz = vzz - cmp(4,i)*cphirec(4,i) - two*cmp(7,i)*cphirec(7,i)
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
c     ##################################################################################
c     ##                                                                              ##
c     ##  subroutine emrealshort1c  --  Ewald short range real space derivs via list  ##
c     ##                                                                              ##
c     ##################################################################################
c
c
c     "emrealshort1c" evaluates the short range real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emrealshort1c
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
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
      use shunt
      use tinheader ,only:ti_p,re_p
      use virial
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis
      integer ii,kkk,iipole,kkpole

      integer iax,iay,iaz
      real(t_p) e,efull,de,f
      real(t_p) bfac
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) scalekk
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) xix,yix,zix
      real(t_p) xiy,yiy,ziy
      real(t_p) xiz,yiz,ziz
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9,rr11
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) dikx,diky,dikz
      real(t_p) dirx,diry,dirz
      real(t_p) dkrx,dkry,dkrz
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrixr,qriyr,qrizr
      real(t_p) qrkxr,qrkyr,qrkzr
      real(t_p) qrrx,qrry,qrrz
      real(t_p) qikrx,qikry,qikrz
      real(t_p) qkirx,qkiry,qkirz
      real(t_p) qikrxr,qikryr,qikrzr
      real(t_p) qkirxr,qkiryr,qkirzr
      real(t_p) diqkx,diqky,diqkz
      real(t_p) dkqix,dkqiy,dkqiz
      real(t_p) diqkxr,diqkyr,diqkzr
      real(t_p) dkqixr,dkqiyr,dkqizr
      real(t_p) dqiqkx,dqiqky,dqiqkz
      real(t_p) dri,drk,qrri,qrrk
      real(t_p) diqrk,dkqri
      real(t_p) dik,qik,qrrik
      real(t_p) term1,term2,term3
      real(t_p) term4,term5,term6
      real(t_p) frcx,frcy,frcz
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) ttmi(3),ttmk(3)
      real(t_p) fix(3),fiy(3),fiz(3)
      real(t_p) bn(0:5)
      real(t_p) s,ds
      real(t_p), allocatable :: mscale(:)
      real(t_p), allocatable :: tem(:,:)
      character*10 mode

 1000 format(' Warning, system moved too much since last neighbor list',
     $   ' update, try lowering nlupdate')
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'emrealshort1c'

c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,nbloc))
c
c     initialize connected atom scaling and torque arrays
c
      mscale = 1.0_ti_p
      tem = 0.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'SHORTEWALD'
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
         do kkk = 1, nshortelst(ii)
            kkpole = shortelst(kkk,ii)
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
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
               rr9 = 7.0_ti_p * rr7 / r2
               rr11 = 9.0_ti_p * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0_ti_p * aewald**2
               alsq2n = 0.0_ti_p
               if (aewald .gt. 0.0_ti_p)
     &            alsq2n = 1.0_ti_p / (sqrtpi*aewald)
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
               qik = 2.0_ti_p*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
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
     &                     - 2.0_ti_p*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0_ti_p*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0_ti_p*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0_ti_p*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0_ti_p*qrrik
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
               scalekk = 1.0_ti_p - mscale(kglob)
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

               call switch_respa(r,off,shortheal,s,ds)
c               call switch_respa(r,off,4.0_ti_p,s,ds)
c               write(*,*) 's = ',s,'ds = ',ds

               em = em + e*s
c
c     calculate intermediate terms for force and torque
c
               de = term1*rr3 + term2*rr5 + term3*rr7
     &                 + term4*rr9 + term5*rr11

               term1 = -ck*rr3 + drk*rr5 - qrrk*rr7
               term2 = ci*rr3 + dri*rr5 + qrri*rr7
               term3 = 2.0_ti_p * rr5
               term4 = 2.0_ti_p * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0_ti_p * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0_ti_p * rr7

c
c     compute the force components for this interaction
c
               frcx = s*(de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx))
     &              - ds*xr*e/r
               frcy = s*(de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry))
     &              - ds*yr*e/r
               frcz = s*(de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz))
     &              - ds*zr*e/r
c
c     compute the torque components for this interaction
c

               ttmi(1) = s*(-rr3*dikx + term1*dirx 
     &                      + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx))
               ttmi(2) = s*(-rr3*diky + term1*diry 
     &                      + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry))
               ttmi(3) = s*(-rr3*dikz + term1*dirz 
     &                      + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz))
               ttmk(1) = s*(rr3*dikx + term2*dkrx 
     &                      - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx))
               ttmk(2) = s*(rr3*diky + term2*dkry 
     &                      - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry))
               ttmk(3) = s*(rr3*dikz + term2*dkrz 
     &                      - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz))
cc
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
            mscale(i12(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0_ti_p
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
         if (iaz .le. 0)  iaz = iglob
         if (iax .le. 0)  iax = iglob
         if (iay .le. 0)  iay = iglob
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
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c     ################################################################################
c     ##                                                                            ##
c     ##  subroutine emreallong1c  --  Ewald long range real space derivs via list  ##
c     ##                                                                            ##
c     ################################################################################
c
c
c     "emreallong1c" evaluates the long range part real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreallong1c
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
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
      use shunt
      use tinheader ,only:ti_p,re_p
      use virial
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis
      integer ii,kkk,iipole,kkpole

      integer iax,iay,iaz
      real(t_p) e,efull,de,f
      real(t_p) bfac,erfc
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) scalekk
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) xix,yix,zix
      real(t_p) xiy,yiy,ziy
      real(t_p) xiz,yiz,ziz
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9,rr11
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) dikx,diky,dikz
      real(t_p) dirx,diry,dirz
      real(t_p) dkrx,dkry,dkrz
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrixr,qriyr,qrizr
      real(t_p) qrkxr,qrkyr,qrkzr
      real(t_p) qrrx,qrry,qrrz
      real(t_p) qikrx,qikry,qikrz
      real(t_p) qkirx,qkiry,qkirz
      real(t_p) qikrxr,qikryr,qikrzr
      real(t_p) qkirxr,qkiryr,qkirzr
      real(t_p) diqkx,diqky,diqkz
      real(t_p) dkqix,dkqiy,dkqiz
      real(t_p) diqkxr,diqkyr,diqkzr
      real(t_p) dkqixr,dkqiyr,dkqizr
      real(t_p) dqiqkx,dqiqky,dqiqkz
      real(t_p) dri,drk,qrri,qrrk
      real(t_p) diqrk,dkqri
      real(t_p) dik,qik,qrrik
      real(t_p) term1,term2,term3
      real(t_p) term4,term5,term6
      real(t_p) frcx,frcy,frcz
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) ttmi(3),ttmk(3)
      real(t_p) fix(3),fiy(3),fiz(3)
      real(t_p) bn(0:5)
      real(t_p) s,ds,mpoleshortcut2
      real(t_p), allocatable :: mscale(:)
      real(t_p), allocatable :: tem(:,:)
      character*10 mode

 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate')
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'emreallong1c'

c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,nbloc))
c
c     initialize connected atom scaling and torque arrays
c
      mscale = 1.0_ti_p
      tem = 0.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
      mpoleshortcut2 = (mpoleshortcut-shortheal)**2
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
            if ((r2 .le. off2).and.(r2.ge.mpoleshortcut2)) then
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
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
               rr9 = 7.0_ti_p * rr7 / r2
               rr11 = 9.0_ti_p * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0_ti_p * aewald**2
               alsq2n = 0.0_ti_p
               if (aewald .gt. 0.0_ti_p)  
     &            alsq2n = 1.0_ti_p / (sqrtpi*aewald)
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
               qik = 2.0_ti_p*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
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
     &                     - 2.0_ti_p*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0_ti_p*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0_ti_p*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0_ti_p*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0_ti_p*qrrik
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
               scalekk = 1.0_ti_p - mscale(kglob)
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

               call switch_respa(r,mpoleshortcut,shortheal,s,ds)

               em = em + (1-s)*e
c
c     calculate intermediate terms for force and torque
c
               de = term1*rr3 + term2*rr5 + term3*rr7
     &                 + term4*rr9 + term5*rr11
               term1 = -ck*rr3 + drk*rr5 - qrrk*rr7
               term2 = ci*rr3 + dri*rr5 + qrri*rr7
               term3 = 2.0_ti_p * rr5
               term4 = 2.0_ti_p * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0_ti_p * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0_ti_p * rr7
c
c     compute the force components for this interaction
c
               frcx = (1-s)*(de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx))
     &                + ds*xr*e/r
               frcy = (1-s)*(de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry))
     &                + ds*yr*e/r
               frcz = (1-s)*(de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz))
     &                + ds*zr*e/r
c
c     compute the torque components for this interaction
c
               ttmi(1) = (1-s)*(-rr3*dikx + term1*dirx 
     &                      + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx))
               ttmi(2) = (1-s)*(-rr3*diky + term1*diry 
     &                      + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry))
               ttmi(3) = (1-s)*(-rr3*dikz + term1*dirz 
     &                      + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz))
               ttmk(1) = (1-s)*(rr3*dikx + term2*dkrx 
     &                      - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx))
               ttmk(2) = (1-s)*(rr3*diky + term2*dkry 
     &                      - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry))
               ttmk(3) = (1-s)*(rr3*dikz + term2*dkrz 
     &                      - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz))
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
            mscale(i12(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0_ti_p
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
         if (iaz .le. 0)  iaz = iglob
         if (iax .le. 0)  iax = iglob
         if (iay .le. 0)  iay = iglob
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
c     OpenMP directives for the major loop structure
c
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
