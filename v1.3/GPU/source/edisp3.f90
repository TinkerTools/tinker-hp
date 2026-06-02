!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine edisp3  --  damped dispersion energy & analysis  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "edisp3" calculates the dispersion energy; also partitions
!     the energy among the atoms
!
!     literature reference:
!
!     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
!     Grounded Damped Dispersion Model with Particle Mesh Ewald
!     Summation", Journal of Chemical Physics, 149, 084115 (2018)
!
!
subroutine edisp3
   use analyz
   use atoms
   use domdec
   use dsppot
   use energi
   use ewald
   use inform
   use iounit
   use potent
   implicit none
   integer i
   real*8 elrc,aelrc
   character*11 mode
!
!
!     choose the method for summing over pairwise interactions
!
   if (use_dewald) then
      call edisp3d
   else
      call edisp3b
   end if
!
!     apply long range dispersion correction if desired
!
   if (use_dcorr .and. .not.use_dewald) then
      mode = 'DISP'
      call evcorr (mode,elrc)
      edsp = edsp + elrc
      aelrc = elrc / dble(n)
      do i = 1, nbloc
         aedsp(i) = aedsp(i) + aelrc
      end do
      if (verbose .and. elrc.ne.0.0d0) then
         if (digits .ge. 8) then
            write (iout,10)  elrc
10          format (/,' Long-Range Dispersion :',9x,f16.8)
         else if (digits .ge. 6) then
            write (iout,20)  elrc
20          format (/,' Long-Range Dispersion :',9x,f16.6)
         else
            write (iout,30)  elrc
30          format (/,' Long-Range Dispersion :',9x,f16.4)
         end if
      end if
   end if
   return
end
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine edisp3b  --  damp dispersion analysis via list  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "edisp3b" calculates the damped dispersion potential energy
!     and also partitions the energy among the atomsusing a pairwise
!     neighbor list
!
!
subroutine edisp3b
   use action
   use analyz
   use atmlst
   use atmtyp
   use atoms
   use bound
   use boxes
   use couple
   use cell
   use cutoff
   use disp
   use domdec
   use dsppot
   use energi
   use group
   use inform
   use inter
   use iounit
   use molcul
   use mutant
   use neigh
   use potent
   use shunt
   use usage
   implicit none
   integer i,j,k,iglob,kglob,kbis,nnvlst
   integer ii,iidisp,kk,kkk
   real*8 xi,yi,zi
   real*8 xr,yr,zr
   real*8 e,fgrp
   real*8 ci,ck
   real*8 r,r2,r3
   real*8 r4,r5,r6
   real*8 ai,ai2
   real*8 ak,ak2
   real*8 di,di2,di3
   real*8 di4,di5
   real*8 dk,dk2,dk3
   real*8 ti,ti2
   real*8 tk,tk2
   real*8 expi,expk
   real*8 damp3,damp5
   real*8 damp,taper
   real*8, allocatable :: dspscale(:)
   real*8 s,ds,dispshortcut2
   logical proceed,usei
   logical muti,mutk
   logical header,huge
   logical testcut,shortrange,longrange,fullrange
   character*11 mode
   character*80 :: RoutineName
!
!     compute the short, long, or full real space part of the Ewald summation
   shortrange = use_dispshort
   longrange  = use_displong
   fullrange  = .not.(shortrange.or.longrange)

   if (shortrange) then
      RoutineName = 'edispshort3b'
      mode        = 'SHORTDISP'
   else if (longrange) then
      RoutineName = 'edisplong3b'
      mode        = 'DISP'
   else
      RoutineName = 'edisp3b'
      mode        = 'DISP'
   endif
!
!     zero out the dispersion interaction
!
   edsp = 0.0d0
   if (ndisp .eq. 0) return
!
!     perform dynamic allocation of some local arrays
!
   allocate (dspscale(n))
!
!     initialize connected atom exclusion coefficients
!
   dspscale = 1.0d0
!
!     set cutoff and switching coefficients
!
   call switch (mode)
   dispshortcut2 = (dispshortcut - shortheal) ** 2
!
!     print header information if debug output was requested
!
   header = .true.
   if (debug .and. ndisp.ne.0) then
      header = .false.
      write (iout,10)
10    format (/,' Individual Dispersion Interactions :',&
         &//,' Type',14x,'Atom Names',15x,'Distance',&
         &8x,'Energy',/)
   end if
!
!     find the damped dispersion energy via neighbor list search
!
   do ii = 1, ndisplocnl
      iidisp = dispglobnl(ii)
      iglob = idisp(iidisp)
      i = loc(iglob)
      ci = csix(iidisp)
      ai = adisp(iidisp)
      usei = use(iglob)
      muti = mut(iglob)
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
!
!     set exclusion coefficients for connected atoms
!
      do j = 1, n12(iglob)
         dspscale(i12(j,iglob)) = dsp2scale
      end do
      do j = 1, n13(iglob)
         dspscale(i13(j,iglob)) = dsp3scale
      end do
      do j = 1, n14(iglob)
         dspscale(i14(j,iglob)) = dsp4scale
      end do
      do j = 1, n15(iglob)
         dspscale(i15(j,iglob)) = dsp5scale
      end do
!
!     decide whether to compute the current interaction
!
      nnvlst = merge(nshortvlst(ii),&
         &nvlst     (ii),&
         &shortrange&
         &)
      do kkk = 1, nnvlst
         kk = merge(shortvlst(kkk,ii),&
            &vlst     (kkk,ii),&
            &shortrange&
            &)
         kglob = idisp(kk)
         kbis = loc(kglob)
         if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
         ck = csix(kk)
         ak = adisp(kk)
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob))
         if (proceed) then
            xr = xi - x(kglob)
            yr = yi - y(kglob)
            zr = zi - z(kglob)
            if (use_bounds) call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.dispshortcut2,&
               &r2 .le. off2,&
               &longrange&
               &)
            if (testcut) then
               r = sqrt(r2)
               r6 = r2**3
               e = -ci * ck / r6
!
!     find the damping factor for the dispersion interaction
!
               di = ai * r
               di2 = di * di
               di3 = di * di2
               dk = ak * r
               expi = exp(-di)
               expk = exp(-dk)
               if (ai .ne. ak) then
                  ai2 = ai * ai
                  ak2 = ak * ak
                  dk2 = dk * dk
                  dk3 = dk * dk2
                  ti = ak2 / (ak2-ai2)
                  ti2 = ti * ti
                  tk = ai2 / (ai2-ak2)
                  tk2 = tk * tk
                  damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi&
                     &- tk2*(1.0d0+dk+0.5d0*dk2)*expk&
                     &- 2.0d0*ti2*tk*(1.0d0+di)*expi&
                     &- 2.0d0*tk2*ti*(1.0d0+dk)*expk
                  damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2&
                     &+di3/6.0d0)*expi&
                     &- tk2*(1.0d0+dk+0.5d0*dk2&
                     &+dk3/6.0d0)*expk&
                     &- 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi&
                     &- 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
               else
                  di4 = di2 * di2
                  di5 = di2 * di3
                  damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2&
                     &+7.0d0*di3/48.0d0+di4/48.0d0)*expi
                  damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2&
                     &+di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
               end if
               damp = 1.5d0*damp5 - 0.5d0*damp3
!
!     apply damping and scaling factors for this interaction
!
               e = e * damp**2 * dspscale(kglob)
!
!     scale the interaction based on its group membership
!
               if (use_group)  e = e * fgrp
!
!     set use of lambda scaling for decoupling or annihilation
!
               if (muti .or. mutk) then
                  if (vcouple .eq. 1) then
                     e = e * vlambda
                  else if (.not.muti .or. .not.mutk) then
                     e = e * vlambda
                  end if
               end if
!
!     use energy switching if near the cutoff distance
!
               if(longrange.or.fullrange) then
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3&
                        &+ c2*r2 + c1*r + c0
                     e = e * taper
                  end if
               end if

               if(shortrange .or. longrange)&
                  &call switch_respa(r,dispshortcut,shortheal,s,ds)

               if(shortrange) then
                  e  =   e * s
               else if(longrange) then
                  e  =   e * (1.0d0 - s)
               endif
!
!     increment the overall dispersion energy components
!
               if (e .ne. 0.0d0) then
                  edsp = edsp + e
                  nedsp = nedsp + 1
                  aedsp(i) = aedsp(i) + 0.5d0*e
                  aedsp(kbis) = aedsp(kbis) + 0.5d0*e
               end if
!
!     increment the total intermolecular energy
!
               if (molcule(iglob) .ne. molcule(kglob)) then
                  einter = einter + e
               end if
!
!     print a message if the energy of this interaction is large
!
               huge = (abs(e) .gt. 4.0d0)
               if ((debug.and.e.ne.0.0d0)&
                  &.or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
20                   format (/,' Individual Dispersion',&
                        &' Interactions :',&
                        &//,' Type',14x,'Atom Names',&
                        &15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,e
30                format (' Disper',4x,2(i7,'-',a3),&
                     &9x,f10.4,2x,f12.4)
               end if
            end if
         end if
      end do
!
!     reset exclusion coefficients for connected atoms
!
      do j = 1, n12(iglob)
         dspscale(i12(j,iglob)) = 1.0d0
      end do
      do j = 1, n13(iglob)
         dspscale(i13(j,iglob)) = 1.0d0
      end do
      do j = 1, n14(iglob)
         dspscale(i14(j,iglob)) = 1.0d0
      end do
      do j = 1, n15(iglob)
         dspscale(i15(j,iglob)) = 1.0d0
      end do
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (dspscale)
   return
end
!
!
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine edisp3d  --  Ewald dispersion analysis via list  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "edisp3d" calculates the damped dispersion energy and analysis
!     using particle mesh Ewald summation and a neighbor list
!
!
subroutine edisp3d
   use action
   use analyz
   use atmlst
   use atoms
   use disp
   use domdec
   use energi
   use ewald
   use pme
   use potent
   implicit none
   integer i,ii,iidisp
   real*8 e
!
!
!     zero out the damped dispersion energy and partitioning terms
!
   nedsp = 0
   edsp = 0.0d0
   aedsp = 0.0d0
   if (ndisp .eq. 0)  return
!
!     set Ewald coefficient
!
   aewald = adewald
!
!
!     compute the reciprocal space part of the Ewald summation
!
   if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))&
      &then
      if (use_disprec) then
         call edrecip
      end if
   end if
!
!     compute the real space portion of the Ewald summation
!
   if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))&
      &then
      if (use_dispreal) then
         call edreal3d
      end if
!
!     compute the self-energy portion of the Ewald summation
!
      if (use_dispself) then
         do ii = 1, ndisploc
            iidisp = dispglob(ii)
            i = idisp(iidisp)
            e = csix(iidisp)**2 * aewald**6 / 12.0d0
            edsp = edsp + e
            nedsp = nedsp + 1
            aedsp(ii) = aedsp(ii) + e
         end do
      end if
   end if
   return
end
!
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine edreal3d  --  real space disp analysis via list  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "edreal3d" evaluated the real space portion of the damped
!     dispersion energy and analysis using Ewald and a neighbor list
!
!
subroutine edreal3d
   use action
   use analyz
   use atmlst
   use atmtyp
   use atoms
   use bound
   use boxes
   use couple
   use cell
   use cutoff
   use disp
   use domdec
   use dsppot
   use energi
   use ewald
   use group
   use inform
   use inter
   use iounit
   use molcul
   use mutant
   use neigh
   use potent
   use shunt
   use usage
   implicit none
   integer i,j,iglob,kglob,kbis,nnvlst
   integer ii,iidisp,kk,kkk
   real*8 xi,yi,zi
   real*8 xr,yr,zr
   real*8 e,efull,fgrp
   real*8 ci,ck
   real*8 r,r2,r6
   real*8 ai,ai2
   real*8 ak,ak2
   real*8 di,di2,di3
   real*8 di4,di5
   real*8 dk,dk2,dk3
   real*8 ti,ti2
   real*8 tk,tk2
   real*8 expi,expk
   real*8 ralpha2
   real*8 expa,term
   real*8 damp3,damp5
   real*8 damp,scale
   real*8, allocatable :: dspscale(:)
   real*8 s,ds,dispshortcut2
   logical proceed,usei
   logical muti,mutk
   logical header,huge
   logical testcut,shortrange,longrange,fullrange
   character*11 mode
   character*80 :: RoutineName

!     compute the short, long, or full real space part of the Ewald summation
   shortrange = use_dispshortreal
   longrange  = use_displong
   fullrange  = .not.(shortrange.or.longrange)

   if (shortrange) then
      RoutineName = 'edrealshort3d'
      mode        = 'SHORTDEWALD'
   else if (longrange) then
      RoutineName = 'edreallong3d'
      mode        = 'DEWALD'
   else
      RoutineName = 'edreal3d'
      mode        = 'DEWALD'
   endif
!
!
!     perform dynamic allocation of some local arrays
!
   allocate (dspscale(n))
!
!     initialize connected atom exclusion coefficients
!
   dspscale = 1.0d0
!
!     set switching coefficients
!
   call switch (mode)
   dispshortcut2 = (dispshortcut - shortheal) ** 2
!
!     print header information if debug output was requested
!
   header = .true.
   if (debug .and. ndisp.ne.0) then
      header = .false.
      write (iout,10)
10    format (/,' Individual Dispersion Interactions :',&
         &//,' Type',14x,'Atom Names',15x,'Distance',&
         &8x,'Energy',/)
   end if
!
!
!     compute the real space portion of the Ewald summation
!
   do ii = 1, ndisplocnl
      iidisp = dispglobnl(ii)
      iglob = idisp(iidisp)
      i = loc(iglob)
      ci = csix(iidisp)
      ai = adisp(iidisp)
      usei = use(iglob)
      muti = mut(iglob)
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
!
!     set exclusion coefficients for connected atoms
!
      do j = 1, n12(iglob)
         dspscale(i12(j,iglob)) = dsp2scale
      end do
      do j = 1, n13(iglob)
         dspscale(i13(j,iglob)) = dsp3scale
      end do
      do j = 1, n14(iglob)
         dspscale(i14(j,iglob)) = dsp4scale
      end do
      do j = 1, n15(iglob)
         dspscale(i15(j,iglob)) = dsp5scale
      end do
      nnvlst = merge(nshortvlst(ii),&
         &nvlst     (ii),&
         &shortrange&
         &)
!
!     decide whether to compute the current interaction
!
      do kkk = 1, nnvlst
         kk = merge(shortvlst(kkk,ii),&
            &vlst     (kkk,ii),&
            &shortrange&
            &)
         kglob = idisp(kk)
         if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
         kbis = loc(kglob)
         ck = csix(kk)
         ak = adisp(kk)
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob))
!
!     compute the energy contribution for this interaction
!
         if (proceed) then
            xr = xi - x(kglob)
            yr = yi - y(kglob)
            zr = zi - z(kglob)
            if (use_bounds) call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.dispshortcut2,&
               &r2 .le. off2,&
               &longrange&
               &)
            if (testcut) then
               r = sqrt(r2)
               r6 = r2**3
               e = -ci * ck / r6
               ralpha2 = r2 * aewald**2
               term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
               expa = exp(-ralpha2) * term
!
!     find the damping factor for the dispersion interaction
!
               di = ai * r
               di2 = di * di
               di3 = di * di2
               dk = ak * r
               expi = exp(-di)
               expk = exp(-dk)
               if (ai .ne. ak) then
                  ai2 = ai * ai
                  ak2 = ak * ak
                  dk2 = dk * dk
                  dk3 = dk * dk2
                  ti = ak2 / (ak2-ai2)
                  ti2 = ti * ti
                  tk = ai2 / (ai2-ak2)
                  tk2 = tk * tk
                  damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi&
                     &- tk2*(1.0d0+dk+0.5d0*dk2)*expk&
                     &- 2.0d0*ti2*tk*(1.0d0+di)*expi&
                     &- 2.0d0*tk2*ti*(1.0d0+dk)*expk
                  damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2&
                     &+di3/6.0d0)*expi&
                     &- tk2*(1.0d0+dk+0.5d0*dk2&
                     &+dk3/6.0d0)*expk&
                     &- 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi&
                     &- 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
               else
                  di4 = di2 * di2
                  di5 = di2 * di3
                  damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2&
                     &+7.0d0*di3/48.0d0+di4/48.0d0)*expi
                  damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2&
                     &+di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
               end if
               damp = 1.5d0*damp5 - 0.5d0*damp3
!
!     apply damping and scaling factors for this interaction
!
               scale = dspscale(kglob) * damp**2
               if (use_group)  scale = scale * fgrp
!
!     set use of lambda scaling for decoupling or annihilation
!
               if (muti .or. mutk) then
                  if (vcouple .eq. 1) then
                     scale = scale * vlambda
                  else if (.not.muti .or. .not.mutk) then
                     scale = scale * vlambda
                  end if
               end if

               if(shortrange .or. longrange)&
                  &call switch_respa(r,dispshortcut,shortheal,s,ds)

               if(shortrange) then
                  e  =   e * s
               else if(longrange) then
                  e  =   e * (1.0d0 - s)
               endif

!
!     compute the full undamped energy for this interaction
!
               efull = e * scale

               if (efull .ne. 0.0d0) then
                  nedsp = nedsp + 1
                  aedsp(i) = aedsp(i) + 0.5d0*efull
                  aedsp(kbis) = aedsp(kbis) + 0.5d0*efull
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + efull
                  end if
               end if
!
!     increment the overall dispersion energy component
!
               scale = scale - 1.0d0
               e = e * (expa+scale)
               edsp = edsp + e
!
!     print a message if the energy of this interaction is large
!
               huge = (abs(efull) .gt. 4.0d0)
               if ((debug.and.efull.ne.0.0d0)&
                  &.or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
20                   format (/,' Individual Dispersion',&
                        &' Interactions :',&
                        &//,' Type',14x,'Atom Names',&
                        &15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  iglob,name(iglob),kglob,&
                     &name(kglob),r,efull
30                format (' Disper',4x,2(i7,'-',a3),9x,&
                     &f10.4,2x,f12.4)
               end if
            end if
         end if
      end do
!
!     reset exclusion coefficients for connected atoms
!
      do j = 1, n12(iglob)
         dspscale(i12(j,iglob)) = 1.0d0
      end do
      do j = 1, n13(iglob)
         dspscale(i13(j,iglob)) = 1.0d0
      end do
      do j = 1, n14(iglob)
         dspscale(i14(j,iglob)) = 1.0d0
      end do
      do j = 1, n15(iglob)
         dspscale(i15(j,iglob)) = 1.0d0
      end do
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (dspscale)
   return
end
