c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar3  --  induced dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar3" calculates the induced dipole polarization energy,
c     and partitions the energy among atoms
c
c
#include "tinker_precision.h"
      subroutine epolar3
      use polpot
      implicit none
c
c     choose the method for summing over polarization interactions
c
      if (polalg.eq.3) then
        !FIXME
        !call epolar3tcg
      else
        call epolar3c
      end if
      return
      end
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3c  --  Ewald polarization analysis; list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3c" calculates the polarization energy and analysis with
c     respect to Cartesian coordinates using particle mesh Ewald and
c     a neighbor list
c
c
      subroutine epolar3c
      use action
      use sizes
      use analyz
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use polar
      use polpot
      use potent
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,ii,ierr
      integer iipole,iglob
      real(t_p) e,f
      real(t_p) term,fterm
      real(t_p) dix,diy,diz
      real(t_p) uix,uiy,uiz,uii
      real(t_p) xd,yd,zd
      real(t_p) xu,yu,zu
c
c
c     zero out the dipole polarization energy and components
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'epolar3c'
      nep = 0
      ep = 0.0_re_p
      aep = 0_ti_p
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      if (use_polarshortreal) then
        call newinduce_shortreal
      else if (use_pmecore) then
        if (polalg.eq.5) then
          call dcinduce_pme
        else
          call newinduce_pme
        end if
      else
        if (polalg.eq.5) then
          call dcinduce_pme2
        else
          call newinduce_pme2
        end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        if (use_prec) then
          call eprecip
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_preal) then
          if (use_polarshortreal) then
            call eprealshort3d
          else
            call epreal3d
          end if
        end if

        if (use_pself) then
c
c     compute the Ewald self-energy term over all the atoms
c
          term = 2.0_ti_p * aewald * aewald
          fterm = -f * aewald / sqrtpi
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
             uii = dix*uix + diy*uiy + diz*uiz
             e = fterm * term * uii / 3.0_ti_p
             ep = ep + e
             nep = nep + 1
             aep(i) = aep(i) + e
          end do
c
c         compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0_ti_p
             yd = 0.0_ti_p
             zd = 0.0_ti_p
             xu = 0.0_ti_p
             yu = 0.0_ti_p
             zu = 0.0_ti_p
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                i = loc(iglob)
                xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
                yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
                zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
                xu = xu + uind(1,iipole)
                yu = yu + uind(2,iipole)
                zu = zu + uind(3,iipole)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xu,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yu,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zu,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
               ep = ep + term*(xd*xu+yd*yu+zd*zu)
               nep = nep + 1
             end if
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                aep(i) = aep(i) + e/real(npole,t_p)
             end do
            end if
        end if
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epreal3d  --  real space polar analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epreal3d" calculates the induced dipole polarization energy
c     and analysis using particle mesh Ewald and a neighbor list
c
c
      subroutine epreal3d
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use chgpot
      use couple
      use domdec
      use energi
      use ewald
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
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,j,k
      integer ii,kk,kkk,iipole,kkpole
      integer iglob,kglob
      real(t_p) e,efull,f
      real(t_p) damp,expdamp
      real(t_p) erfc,bfac
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) pdi,pti,pgamma
      real(t_p) sc3,sc5,sc7
      real(t_p) psc3,psc5,psc7
      real(t_p) psr3,psr5,psr7
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1
      real(t_p) rr3,rr5,rr7
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) uix,uiy,uiz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) ukx,uky,ukz
      real(t_p) dri,drk,uri,urk
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrri,qrrk
      real(t_p) duik,quik
      real(t_p) term1,term2,term3
      real(t_p) bn(0:3)
      real(t_p), allocatable :: pscale(:)
      logical header,huge
      character*10 mode
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'epreal3d'

c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      pscale = 1.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5_ti_p * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
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
         uix = uind(1,iipole)
         uiy = uind(2,iipole)
         uiz = uind(3,iipole)
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
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kkpole = elst(kkk,ii)
            kglob = ipole(kkpole)
            k = loc(kglob)
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
               ukx = uind(1,kkpole)
               uky = uind(2,kkpole)
               ukz = uind(3,kkpole)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
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
               do j = 1, 3
                  bfac = real(j+j-1,t_p)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 3
                  bn(j) = f * bn(j)
               end do
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0_ti_p
               sc5 = 1.0_ti_p
               sc7 = 1.0_ti_p
               damp = pdi * pdamp(kkpole)
               if (damp .ne. 0.0_ti_p) then
                  pgamma = min(pti,thole(kkpole))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0_ti_p) then
                     expdamp = exp(damp)
                     sc3 = 1.0_ti_p - expdamp
                     sc5 = 1.0_ti_p - (1.0_ti_p-damp)*expdamp
                     sc7 = 1.0_ti_p - (1.0_ti_p-damp+0.6_ti_p*damp**2)
     &                                    *expdamp
                  end if
               end if
c
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
               duik = dix*ukx + diy*uky + diz*ukz
     &                   + dkx*uix + dky*uiy + dkz*uiz
               quik = qrix*ukx + qriy*uky + qriz*ukz
     &                   - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0_ti_p*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     intermediates involving Thole damping and scale factors
c
               psr3 = rr3 * sc3 * pscale(kglob)
               psr5 = rr5 * sc5 * pscale(kglob)
               psr7 = rr7 * sc7 * pscale(kglob)
c
c     compute the full undamped energy for this interaction
c
               efull = term1*psr3 + term2*psr5 + term3*psr7
               if (efull .ne. 0.0_ti_p) then
                  nep = nep + 1
                  aep(i) = aep(i) + efull
                  if (molcule(iglob) .ne. molcule(kglob))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
c
               psc3 = 1.0_ti_p - sc3*pscale(kglob)
               psc5 = 1.0_ti_p - sc5*pscale(kglob)
               psc7 = 1.0_ti_p - sc7*pscale(kglob)
               psr3 = bn(1) - psc3*rr3
               psr5 = bn(2) - psc5*rr5
               psr7 = bn(3) - psc7*rr7
c
c     compute the energy contribution for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
               ep = ep + e
c
c     print message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0_ti_p)
               if ((debug.and.efull.ne.0.0_ti_p)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  iglob,name(iglob),kk,name(kglob),
     &              r,efull
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine eprealshort3d  --  real space polar analysis via list  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "epreal3d" calculates the short range real space induced dipole polarization energy
c     and analysis using particle mesh Ewald and a neighbor list
c
c
      subroutine eprealshort3d
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use ewald
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
      use sizes
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,j,k
      integer ii,kk,kkk,iipole,kkpole
      integer iglob,kglob
      real(t_p) e,efull,f
      real(t_p) damp,expdamp
      real(t_p) bfac
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) pdi,pti,pgamma
      real(t_p) sc3,sc5,sc7
      real(t_p) psc3,psc5,psc7
      real(t_p) psr3,psr5,psr7
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1
      real(t_p) rr3,rr5,rr7
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) uix,uiy,uiz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) ukx,uky,ukz
      real(t_p) dri,drk,uri,urk
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) qrri,qrrk
      real(t_p) duik,quik
      real(t_p) term1,term2,term3
      real(t_p) bn(0:3)
      real(t_p) s,ds
      real(t_p), allocatable :: pscale(:)
      logical header,huge
      character*10 mode
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'eprealshort3d'
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      pscale = 1.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5_ti_p * electric / dielec
      mode = 'MPOLE'
      mode = 'SHORTEWALD'
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
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
         uix = uind(1,iipole)
         uiy = uind(2,iipole)
         uiz = uind(3,iipole)
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
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kkpole = elst(kkk,ii)
            kglob = ipole(kkpole)
            k = loc(kglob)
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
               ukx = uind(1,kkpole)
               uky = uind(2,kkpole)
               ukz = uind(3,kkpole)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
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
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 3
                  bn(j) = f * bn(j)
               end do
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0_ti_p
               sc5 = 1.0_ti_p
               sc7 = 1.0_ti_p
               damp = pdi * pdamp(kkpole)
               if (damp .ne. 0.0_ti_p) then
                  pgamma = min(pti,thole(kkpole))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0_ti_p) then
                     expdamp = exp(damp)
                     sc3 = 1.0_ti_p - expdamp
                     sc5 = 1.0_ti_p - (1.0_ti_p-damp)*expdamp
                     sc7 = 1.0_ti_p - (1.0_ti_p-damp+0.6_ti_p*damp**2)
     &                                    *expdamp
                  end if
               end if
c
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
               duik = dix*ukx + diy*uky + diz*ukz
     &                   + dkx*uix + dky*uiy + dkz*uiz
               quik = qrix*ukx + qriy*uky + qriz*ukz
     &                   - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0_ti_p*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     intermediates involving Thole damping and scale factors
c
               psr3 = rr3 * sc3 * pscale(kglob)
               psr5 = rr5 * sc5 * pscale(kglob)
               psr7 = rr7 * sc7 * pscale(kglob)
c
c     compute the full undamped energy for this interaction
c
               efull = term1*psr3 + term2*psr5 + term3*psr7
               if (efull .ne. 0.0_ti_p) then
                  nep = nep + 1
                  aep(i) = aep(i) + efull
                  if (molcule(iglob) .ne. molcule(kglob))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
c
               psc3 = 1.0_ti_p - sc3*pscale(kglob)
               psc5 = 1.0_ti_p - sc5*pscale(kglob)
               psc7 = 1.0_ti_p - sc7*pscale(kglob)
               psr3 = bn(1) - psc3*rr3
               psr5 = bn(2) - psc5*rr5
               psr7 = bn(3) - psc7*rr7
c
c     compute the energy contribution for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
               call switch_respa(r,off,shortheal,s,ds)
               ep = ep + e*s
c
c     print message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0_ti_p)
               if ((debug.and.efull.ne.0.0_ti_p)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  iglob,name(iglob),kk,name(kglob),
     &              r,efull
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     subroutine totaldipole: get the total dipole moment of the system (permanent+induced)
c
      subroutine totaldipole
      use atmlst
      use atoms
      use boxes
      use charge
      use domdec
      use iounit
      use mpole
      use polar
      use potent
      use tinheader ,only:ti_p,re_p
      use units
      use mpi
      implicit none
      integer i,iipole,iglob,ierr
      integer iichg
      real(t_p) q,xr,yr,zr
      real(t_p) dipx,dipy,dipz
      real(t_p) mux,muy,muz,mudx,mudy,mudz,mupx,mupy,mupz
 1000 format(/'x dipolar moment : ',F14.5)
 1010 format(/'y dipolar moment : ',F14.5)
 1020 format(/'z dipolar moment : ',F14.5)
c
      dipx = 0.0_ti_p
      dipy = 0.0_ti_p
      dipz = 0.0_ti_p
      if (use_mpole) then
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          xr = x(iglob)
          yr = y(iglob)
          zr = z(iglob)
          q = rpole(1,iipole)
          mux = rpole(2,iipole)
          muy = rpole(3,iipole)
          muz = rpole(4,iipole)
          mudx = uind(1,iipole)
          mudy = uind(2,iipole)
          mudz = uind(3,iipole)
          mupx = uinp(1,iipole)
          mupy = uinp(2,iipole)
          mupz = uinp(3,iipole)
          dipx = dipx + q*xr + mux + 0.5*(mudx+mupx)
          dipy = dipy + q*yr + muy + 0.5*(mudy+mupy)
          dipz = dipz + q*zr + muz + 0.5*(mudz+mupz)
        end do
      else if (use_charge) then
        do i = 1, nionloc
          iichg = chgglob(i)
          iglob = iion(iichg)
          xr = x(iglob)
          yr = y(iglob)
          zr = z(iglob)
          q = pchg(iichg)
          dipx = dipx + q*xr 
          dipy = dipy + q*yr 
          dipz = dipz + q*zr 
        end do
      end if

      dipx = debye*dipx
      dipy = debye*dipy
      dipz = debye*dipz

      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,dipx,1,MPI_TPREC,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,dipy,1,MPI_TPREC,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,dipz,1,MPI_TPREC,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      else 
        call MPI_REDUCE(dipx,dipx,1,MPI_TPREC,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(dipy,dipy,1,MPI_TPREC,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(dipz,dipz,1,MPI_TPREC,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      end if
      if (rank.eq.0) write(iout,1000) dipx 
      if (rank.eq.0) write(iout,1010) dipy 
      if (rank.eq.0) write(iout,1020) dipz 
      return
      end
