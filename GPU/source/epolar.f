c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar  --  induced dipole polarization energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar" calculates the polarization energy due to induced
c     dipole interactions
c
c
#include "tinker_macro.h"
      subroutine epolar
      use polpot
      use group
      use tinheader,only: ti_p
      implicit none
      if(use_group .and. wgrp(1,2).eq.0._ti_p) then
        return
      endif
c
      if (polalg.eq.3) then
        !FIXME
        !call epolar3tcg
      else
        call epolar0c
      end if
      return
      end
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar0c  --  Ewald polarization derivs via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar0c" calculates the dipole polarization energy with respect
c     to Cartesian coordinates using particle mesh Ewald summation and
c     a neighbor list
c
c
      subroutine epolar0c
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use group
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
      integer ii,iglob,iipole,ierr
      real(t_p) e,f,term,fterm
      real(t_p) dix,diy,diz
      real(t_p) uix,uiy,uiz,uii
      real(t_p) xd,yd,zd
      real(t_p) xu,yu,zu
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0_re_p
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
            call epreal0c
        end if

        if (use_pself) then
c
c     compute the Ewald self-energy term over all the atoms
c
          term = 2.0_ti_p * aewald * aewald
          fterm = -f * aewald / sqrtpi
          do ii = 1, npoleloc
             iipole = poleglob(ii)
             dix = rpole(2,iipole)
             diy = rpole(3,iipole)
             diz = rpole(4,iipole)
             uix = uind(1,iipole)
             uiy = uind(2,iipole)
             uiz = uind(3,iipole)
             uii = dix*uix + diy*uiy + diz*uiz
             e = fterm * term * uii / 3.0_ti_p
             ep = ep + e
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
             end if
          end if
        end if
      end if

      !if (use_group .and. use_group_polar) call switch_group(.false.)

      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal0c  --  real space polar energy via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal0c" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine epreal0c
      use sizes
      use atmlst
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use ewald
      use group
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use tinheader ,only: ti_p
      use mpi
      implicit none
      integer i,j,k,iglob,kglob,kbis,nnelst
      integer ii,kkk,iipole,kkpole
      real(t_p) e,f
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) sr3,sr5,sr7
      real(t_p) r,r2,rr3,rr5,rr7
      real(t_p) rr3i,rr5i,rr7i
      real(t_p) rr3k,rr5k,rr7k
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) uix,uiy,uiz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) ukx,uky,ukz
      real(t_p) dir,diu,qiu,uir
      real(t_p) dkr,dku,qku,ukr
      real(t_p) qix,qiy,qiz,qir
      real(t_p) qkx,qky,qkz,qkr
      real(t_p) corei,corek
      real(t_p) vali,valk
      real(t_p) alphai,alphak
      real(t_p) term1,term2,term3
      real(t_p) dmpi(7),dmpk(7)
      real(t_p) dmpik(7),dmpe(7)
      real(t_p) scalek
      real(t_p), allocatable :: pscale(:)
      logical shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName

c     compute the short, or full real space part of the summation
      shortrange = use_polarshortreal
      longrange  = .false.
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'eprealshort0c'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'epreallong0c'
         mode        = 'EWALD'
      else
         RoutineName = 'epreal0c'
         mode        = 'MPOLE'
      endif
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      pscale = 1.0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5 * electric / dielec
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
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
c
c     evaluate all sites within the cutoff distance
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
            kglob = ipole(kkpole)
            kbis = loc(kglob)
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
c     calculate real space Ewald error function damping
c
               call dampewald (7,r,r2,f,dmpe)
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  call dampthole (iipole,kkpole,7,r,dmpik)
                  scalek = pscale(kglob)
                  rr3 = f/ (r*r2)
                  rr5 = 3.0 * rr3 / r2
                  rr7 = 5.0 * rr5 / r2
                  sr3 = scalek*dmpik(3) * rr3
                  sr5 = scalek*dmpik(5) * rr5
                  sr7 = scalek*dmpik(7) * rr7
                  sr3 = dmpe(3) - rr3 + sr3
                  sr5 = dmpe(5) - rr5 + sr5
                  sr7 = dmpe(7) - rr7 + sr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0*(qiu-qku) - uir*dkr - dir*ukr
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
                  rr5 = 3.0 * rr3 / r2
                  rr7 = 5.0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  rr3 = f / (r*r2)
                  rr5 = 3.0 * rr3 / r2
                  rr7 = 5.0 * rr5 / r2
                  rr3i = dmpe(3) - rr3 + rr3i
                  rr5i = dmpe(5) - rr5 + rr5i
                  rr7i = dmpe(7) - rr7 + rr7i
                  rr3k = dmpe(3) - rr3 + rr3k
                  rr5k = dmpe(5) - rr5 + rr5k
                  rr7k = dmpe(7) - rr7 + rr7k
                  rr3 = dmpe(3) - (1.0-scalek)*rr3
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
               !if(use_group) then
               !   call groups(fgrp,iglob,kglob,0,0,0,0)
               !   sc3 = sc3 * fgrp
               !   sc5 = sc5 * fgrp
               !   sc7 = sc7 * fgrp
               !endif
c
c     compute the energy contribution for this interaction
c
               ep = ep + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = 1.0
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = 1.0
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = 1.0
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = 1.0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine eprecip  --  PME recip space polarization energy  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "eprecip" evaluates the reciprocal space portion of particle
c     mesh Ewald summation energy due to dipole polarization
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine eprecip
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use fft
      use inform   ,only: deb_Path
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use tinheader ,only: ti_p,re_p
      use mpi
      implicit none
      integer iipole
      integer commloc
      integer nprocloc,rankloc
      integer i,j,k,iglob
      real(t_p) e
      real(t_p) f
      real(t_p) expterm
      real(t_p) struc2
      real(t_p) a(3,3),ftc(10,10)
      real(t_p) fuind(3)
c
      if (deb_Path) write(*,*) 'eprecip'
      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc =  comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = COMM_TINKER
      end if
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 0.000001)  return
      f = electric / dielec
cc
cc     get the fractional to Cartesian transformation matrix
cc
      call frac_to_cart (ftc)
cc
cc     initialize variables required for the scalar summation
cc
c      pterm = (pi/aewald)**2
c      volterm = pi * volbox
c      nff = nfft1 * nfft2
c      nf1 = (nfft1+1) / 2
c      nf2 = (nfft2+1) / 2
c      nf3 = (nfft3+1) / 2
cc
cc     remove scalar sum virial from prior multipole 3-D FFT
cc
c      if (.not. use_mpole) then
c         call bspline_fill
c         call table_fill
cc
cc     assign only the permanent multipoles to the PME grid
cc     and perform the 3-D FFT forward transformation
cc
c         do i = 1, npole
c            cmp(1,i) = rpole(1,i)
c            cmp(2,i) = rpole(2,i)
c            cmp(3,i) = rpole(3,i)
c            cmp(4,i) = rpole(4,i)
c            cmp(5,i) = rpole(5,i)
c            cmp(6,i) = rpole(9,i)
c            cmp(7,i) = rpole(13,i)
c            cmp(8,i) = 2.0_ti_p * rpole(6,i)
c            cmp(9,i) = 2.0_ti_p * rpole(7,i)
c            cmp(10,i) = 2.0_ti_p * rpole(10,i)
c         end do
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
cc
cc     make the scalar summation over reciprocal lattice
cc
c         do i = 1, ntot-1
c            k3 = i/nff + 1
c            j = i - (k3-1)*nff
c            k2 = j/nfft1 + 1
c            k1 = j - (k2-1)*nfft1 + 1
c            m1 = k1 - 1
c            m2 = k2 - 1
c            m3 = k3 - 1
c            if (k1 .gt. nf1)  m1 = m1 - nfft1
c            if (k2 .gt. nf2)  m2 = m2 - nfft2
c            if (k3 .gt. nf3)  m3 = m3 - nfft3
c            r1 = real(m1,t_p)
c            r2 = real(m2,t_p)
c            r3 = real(m3,t_p)
c            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c            hsq = h1*h1 + h2*h2 + h3*h3
c            term = -pterm * hsq
c            expterm = 0.0_ti_p
c            if (term .gt. -50.0_ti_p) then
c               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c               expterm = exp(term) / denom
c               if (.not. use_bounds) then
c                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
c               else if (octahedron) then
c                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
c               end if
c            end if
c            qfac(k1,k2,k3) = expterm
c         end do
cc
cc     account for zeroth grid point for nonperiodic system
cc
c         qfac(1,1,1) = 0.0_ti_p
c         if (.not. use_bounds) then
c            expterm = 0.5_ti_p * pi / xbox
c            qfac(1,1,1) = expterm
c         end if
cc
cc     complete the transformation of the PME grid
cc
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  term = qfac(i,j,k)
c                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
c                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
c               end do
c            end do
c         end do
cc
cc     perform 3-D FFT backward transform and get potential
cc
c         call fftback
c         call fphi_mpole (fphi)
c      end if
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = real(nfft1,t_p) * recip(i,1)
         a(2,i) = real(nfft2,t_p) * recip(i,2)
         a(3,i) = real(nfft3,t_p) * recip(i,3)
      end do
      e = 0_ti_p
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         fuind = 0_ti_p
         do j = 1, 3
            fuind(j) = a(j,1)*uind(1,iipole) + a(j,2)*uind(2,iipole)
     &                      + a(j,3)*uind(3,iipole)
         end do
         do k = 1, 3
            e = e + fuind(k)*fphirec(k+1,i)
         end do
      end do
      e = 0.5_ti_p * electric*  e
      ep = ep + e
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5_ti_p * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e = f * expterm * struc2
           ep = ep + e
        end if
      end if
c
      return
      end
