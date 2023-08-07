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
      subroutine epolar3
      use polpot
      implicit none
c
c     choose the method for summing over polarization interactions
c
      if (polalg.eq.3) then
        call epolar3tcg
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
      use analyz
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use group
      use math
      use mpole
      use polar
      use polpot
      use potent
      use mpi
      implicit none
      integer i,ii,ierr
      integer iipole,iglob
      real*8 e,f
      real*8 term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the dipole polarization energy and components
c
      nep = 0
      ep = 0.0d0
      aep = 0d0
      if (npole .eq. 0)  return
      aewald = apewald
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
          call epreal3d
        end if

        if (use_pself) then
c
c     compute the Ewald self-energy term over all the atoms
c
          term = 2.0d0 * aewald * aewald
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
             e = fterm * term * uii / 3.0d0
             ep = ep + e
             nep = nep + 1
             aep(i) = aep(i) + e
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
                xd = xd + rpole(2,iipole) + rpole(1,iipole)*x(iglob)
                yd = yd + rpole(3,iipole) + rpole(1,iipole)*y(iglob)
                zd = zd + rpole(4,iipole) + rpole(1,iipole)*z(iglob)
                xu = xu + uind(1,iipole)
                yu = yu + uind(2,iipole)
                zu = zu + uind(3,iipole)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,xu,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yu,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zu,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               term = (2.0d0/3.0d0) * f * (pi/volbox)
               e = term*(xd*xu+yd*yu+zd*zu)
               ep = ep + e
               nep = nep + 1
             end if
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                aep(i) = aep(i) + e/dble(npole)
             end do
          end if
        end if
      end if
c
c     get group polarization if necessary
c
      if (use_group) call switch_group_ene

c
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
c     if shortrange, calculates just the short range part
c
      subroutine epreal3d
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
      integer i,j,k,iglob,kglob,kbis,nnelst
      integer ii,kkk,iipole,kkpole
      real*8 e,efull,f
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
      real*8 dmpik(7),dmpe(7)
      real*8 scalek
      logical shortrange,longrange,fullrange
      real*8, allocatable :: pscale(:)
      logical header,huge
      character*11 mode
      character*80 :: RoutineName



c     compute the short, or full real space part of the summation
      shortrange = use_polarshortreal
      longrange  = .false.
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'eprealshort3d'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'epreallong3d'
         mode        = 'EWALD'
      else
         RoutineName = 'epreal3d'
         mode        = 'EWALD'
      endif
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
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  sr3 = scalek * dmpik(3) * rr3
                  sr5 = scalek * dmpik(5) * rr5
                  sr7 = scalek * dmpik(7) * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  efull = term1*sr3 + term2*sr5 + term3*sr7
                  sr3 = dmpe(3) - rr3 + sr3
                  sr5 = dmpe(5) - rr5 + sr5
                  sr7 = dmpe(7) - rr7 + sr7
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
                  efull = uir*(corek*rr3+valk*rr3k)
     &                       - ukr*(corei*rr3+vali*rr3i)
     &                       + diu*rr3i + dku*rr3k
     &                       + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                       - dkr*uir*rr5k - dir*ukr*rr5i
     &                       + qkr*uir*rr7k - qir*ukr*rr7i
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpe(3) - rr3 + rr3i
                  rr5i = dmpe(5) - rr5 + rr5i
                  rr7i = dmpe(7) - rr7 + rr7i
                  rr3k = dmpe(3) - rr3 + rr3k
                  rr5k = dmpe(5) - rr5 + rr5k
                  rr7k = dmpe(7) - rr7 + rr7k
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     compute the energy contribution for this interaction
c
               ep = ep + e
               if (efull .ne. 0.0d0) then
                  nep = nep + 1
                  aep(i) = aep(i) + 0.5d0*efull
                  aep(kbis) = aep(kbis) + 0.5d0*efull
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + efull
                  end if
               end if
c
c     print message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 10.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),kglob,name(kglob),r,efull
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
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
      use units
      use mpi
      use spectra, only: compute_dipole
      implicit none
      integer i,iipole,iglob,ierr
      integer iichg
      real*8 q,xr,yr,zr
      real*8 dip(3),dipind(3)
      real*8 dipx,dipy,dipz
      real*8 mux,muy,muz,mudx,mudy,mudz,mupx,mupy,mupz
 1000 format(/'x dipolar moment : ',F14.5)
 1010 format(/'y dipolar moment : ',F14.5)
 1020 format(/'z dipolar moment : ',F14.5)
c
      call compute_dipole(dip,dipind,.TRUE.)

      dipx = debye*(dip(1) + dipind(1))
      dipy = debye*(dip(2) + dipind(2))
      dipz = debye*(dip(3) + dipind(3))

      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,dipx,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,dipy,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,dipz,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      else 
        call MPI_REDUCE(dipx,dipx,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(dipy,dipy,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(dipz,dipz,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      end if
      if (rank.eq.0) write(iout,1000) dipx 
      if (rank.eq.0) write(iout,1010) dipy 
      if (rank.eq.0) write(iout,1020) dipz 
      return
      end
