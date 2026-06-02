!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine echgtrn3  --  charge transfer energy & analysis  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "echgtrn3" calculates the charge transfer energy; also partitions
!     the energy among the atoms
!
!
!> @brief 
!> calculates the charge transfer energy; also partitions
!> the energy among the atoms
!> @param no params
subroutine echgtrn3
   use inform
   use iounit
   implicit none
!
   if (deb_Path) write(iout,*), 'echgtrn3 '
!
!
!     choose method for summing over charge transfer interactions
!
   call echgtrn3c
   return
end
!
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine echgtrn3c  --  neighbor list chgtrn analysis  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "echgtrn3c" calculates the charge transfer interaction energy
!     and also partitions the energy among the atoms using a pairwise
!     neighbor list
!
!
!> @brief 
!> calculates the charge transfer interaction energy
!> and also partitions the energy among the atoms using a pairwise
!> neighbor list
!> @param no params
subroutine echgtrn3c
   use action
   use analyz
   use atmtyp
   use atmlst
   use atoms
   use bound
   use chgpot
   use chgtrn
   use cell
   use couple
   use ctrpot
   use cutoff
   use domdec
   use energi
   use group
   use inform
   use inter
   use iounit
   use molcul
   use mplpot
   use mpole
   use mutant
   use neigh
   use potent
   use shunt
   use usage
   implicit none
   integer i,j,iipole,iglob,kglob,nnelst
   integer ii,kkk,kbis,kkpole
   real*8 e,f,fgrp
   real*8 r,r2,r3
   real*8 r4,r5
   real*8 xi,yi,zi
   real*8 xr,yr,zr
   real*8 chgi,chgk
   real*8 chgik
   real*8 alphai,alphak
   real*8 expi,expk
   real*8 expik
   real*8 taper
   real*8 s,ds,ctrnshortcut2,facts
   real*8, allocatable :: mscale(:)
   logical proceed,usei
   logical header,huge
   logical muti,mutk
   logical testcut,shortrange,longrange,fullrange
   character*11 mode
   character*80 :: RoutineName
!
   if (deb_Path) write(iout,*), 'echgtrn3c '
!
   shortrange = use_chgtrnshort
   longrange  = use_chgtrnlong
   fullrange  = .not.(shortrange.or.longrange)

   if (shortrange) then
      RoutineName = 'echgtrnshort3c'
      mode        = 'SHORTCHGTRN'
   else if (longrange) then
      RoutineName = 'echgtrnlong3c'
      mode        = 'CHGTRN'
   else
      RoutineName = 'echgtrn3c'
      mode        = 'CHGTRN'
   endif
!
!     zero out the charge transfer energy and partitioning terms
!
   nect = 0
   ect = 0.0d0
   aect = 0.0d0
   if (nct .eq. 0)  return
!
!     perform dynamic allocation of some local arrays
!
   allocate (mscale(n))
!
!     initialize connected atom exclusion coefficients
!
   mscale = 1.0d0
!
!     set the coefficients for the switching function
!
   call switch (mode)
   ctrnshortcut2 = (ctrnshortcut-shortheal)**2
!
!     set conversion factor
!
   f = electric / dielec
!
!     print header information if debug output was requested
!
   header = .true.
   if (debug .and. nct.ne.0) then
      header = .false.
      write (iout,10)
10    format (/,' Individual Charge Transfer Interactions :',&
      &//,' Type',14x,'Atom Names',15x,'Distance',&
      &8x,'Energy',/)
   end if
!
!     compute the charge transfer energy
!
   do ii = 1, npolelocnl
      iipole = poleglobnl(ii)
      iglob = ipole(iipole)
      i = loc(iglob)
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      chgi = chgct(iipole)
      alphai = dmpct(iipole)
      if (alphai .eq. 0.0d0)  alphai = 1000.0d0
      usei = use(iglob)
      muti = mut(iglob)
!
!     set exclusion coefficients for connected atoms
!
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
!
!     evaluate all sites within the cutoff distance
!
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
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob))
         if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
         if (proceed) then
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.ctrnshortcut2,&
            &r2 .le. off2,&
            &longrange&
            &)
            if (testcut) then
               r = sqrt(r2)
               chgk = chgct(kkpole)
               alphak = dmpct(kkpole)
               if (alphak .eq. 0.0d0)  alphak = 1000.0d0
               if (ctrntyp .eq. 'SEPARATE') then
                  expi = exp(-alphai*r)
                  expk = exp(-alphak*r)
                  e = -chgi*expk - chgk*expi
               else
                  chgik = sqrt(abs(chgi*chgk))
                  expik = exp(-0.5d0*(alphai+alphak)*r)
                  e = -chgik * expik
               end if
               e = f * e * mscale(kglob)
!
!     apply lambda scaling for interaction annihilation 
!
               if (muti .or. mutk) then
                  e = e * elambda
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
!
!     scale the interaction based on its group membership
!
               if (use_group)  e = e * fgrp

               if(shortrange .or. longrange)&
               &call switch_respa(r,ctrnshortcut,shortheal,s,ds)

               if(shortrange) then
                  facts =         s
               else if(longrange) then
                  facts = 1.0d0 - s
               else
                  facts = 1.0d0
               endif

               e  = e * facts
!
!     increment the overall charge transfer energy components
!
               if (e .ne. 0.0d0) then
                  nect = nect + 1
                  ect = ect + e
                  aect(i) = aect(i) + 0.5d0*e
                  aect(kbis) = aect(kbis) + 0.5d0*e
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
               huge = (e .gt. 10.0d0)
               if ((debug.and.e.ne.0.0d0)&
               &.or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
20                   format (/,' Individual Charge Transfer',&
                     &' Interactions :',&
                     &//,' Type',14x,'Atom Names',&
                     &15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  iglob,name(iglob),kglob,&
                  &name(kglob),r,e
30                format (' ChgTrn',4x,2(i7,'-',a3),&
                  &9x,f10.4,2x,f12.4)
               end if
            end if
         end if
      end do
!
!     reset exclusion coefficients for connected atoms
!
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
!
!     perform deallocation of some local arrays
!
   deallocate (mscale)
   return
end
