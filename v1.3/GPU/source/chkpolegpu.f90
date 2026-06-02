!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine chkpole  --  check multipoles at chiral sites  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "chkpole" inverts atomic multipole moments as necessary
!     at sites with chiral local reference frame definitions
!
!
#include "tinker_precision.h"
subroutine chkpolegpu(init)
   use sizes
   use atmlst
   use atoms
   use mpole
   use utilgpu,only:openacc_abort,rec_queue
   use tinheader,only: ti_p
   implicit none
   integer k,ii,iipole
   integer ia,ib,ic,id
   real(t_p) xad,yad,zad
   real(t_p) xbd,ybd,zbd
   real(t_p) xcd,ycd,zcd
   real(t_p) c1,c2,c3,vol
   logical check
   logical,intent(in):: init

   if (init) call openacc_abort(" init is not supposed to be true"//&
      &" in chkpolegpu routine")
!
!
!     loop over multipole sites testing for chirality inversion
!
!$acc parallel loop present(poleglobnl,ipole,xaxis,yaxis,zaxis,pole) &
!$acc         async(rec_queue)
   do ii = 1, npolelocnl
      iipole = poleglobnl(ii)
      check  = .true.
      if (ipolaxe(iipole).ne.Ax_Z_Then_X) check = .false.
      if (yaxis(iipole) .eq. 0)  check = .false.
      if (check) then
         k   = yaxis(iipole)
         ia  = ipole(iipole)
         ib  = zaxis(iipole)
         ic  = xaxis(iipole)
         id  = abs(k)
!
!     compute the signed parallelpiped volume at chiral site
!
         xad = x(ia)   - x(id)
         yad = y(ia)   - y(id)
         zad = z(ia)   - z(id)
         xbd = x(ib)   - x(id)
         ybd = y(ib)   - y(id)
         zbd = z(ib)   - z(id)
         xcd = x(ic)   - x(id)
         ycd = y(ic)   - y(id)
         zcd = z(ic)   - z(id)
         c1  = ybd*zcd - zbd*ycd
         c2  = ycd*zad - zcd*yad
         c3  = yad*zbd - zad*ybd
         vol = xad*c1  + xbd*c2 + xcd*c3
!
!     invert atomic multipole components involving the y-axis
!
         if (k.lt.0.and.vol.gt.0.0_ti_p .or.&
            &k.gt.0.and.vol.lt.0.0_ti_p) then
            yaxis(iipole) = -k
            pole( 3,iipole) = -pole( 3,iipole)
            pole( 6,iipole) = -pole( 6,iipole)
            pole( 8,iipole) = -pole( 8,iipole)
            pole(10,iipole) = -pole(10,iipole)
            pole(12,iipole) = -pole(12,iipole)
         end if
      end if
   end do
end

subroutine chkpolegpu_group()
   use sizes
   use atmlst
   use atoms
   use mpole
   use group
   use utilgpu  ,only:openacc_abort,rec_queue
   use tinheader,only: ti_p
   implicit none
   integer k,ii,iipole
   integer ia,ib,ic,id
   real(t_p) xad,yad,zad
   real(t_p) xbd,ybd,zbd
   real(t_p) xcd,ycd,zcd
   real(t_p) c1,c2,c3,vol
   logical check

!
!     loop over multipole sites testing for chirality inversion
!
!$acc parallel loop present(globpolegroup,ipole,xaxis,yaxis,zaxis,pole) &
!$acc         async(rec_queue)
   do ii = 1, npolegroup
      iipole = globpolegroup(ii)
      check  = .true.
      if (ipolaxe(iipole).ne.Ax_Z_Then_X) check = .false.
      if (yaxis(iipole) .eq. 0)  check = .false.
      if (check) then
         k   = yaxis(iipole)
         ia  = ipole(iipole)
         ib  = zaxis(iipole)
         ic  = xaxis(iipole)
         id  = abs(k)
!
!     compute the signed parallelpiped volume at chiral site
!
         xad = x(ia)   - x(id)
         yad = y(ia)   - y(id)
         zad = z(ia)   - z(id)
         xbd = x(ib)   - x(id)
         ybd = y(ib)   - y(id)
         zbd = z(ib)   - z(id)
         xcd = x(ic)   - x(id)
         ycd = y(ic)   - y(id)
         zcd = z(ic)   - z(id)
         c1  = ybd*zcd - zbd*ycd
         c2  = ycd*zad - zcd*yad
         c3  = yad*zbd - zad*ybd
         vol = xad*c1  + xbd*c2 + xcd*c3
!
!     invert atomic multipole components involving the y-axis
!
         if (k.lt.0.and.vol.gt.0.0_ti_p .or.&
            &k.gt.0.and.vol.lt.0.0_ti_p) then
            yaxis(iipole) = -k
            pole( 3,iipole) = -pole( 3,iipole)
            pole( 6,iipole) = -pole( 6,iipole)
            pole( 8,iipole) = -pole( 8,iipole)
            pole(10,iipole) = -pole(10,iipole)
            pole(12,iipole) = -pole(12,iipole)
         end if
      end if
   end do
end

