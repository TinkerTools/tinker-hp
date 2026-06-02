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
subroutine chkpole(init)
   use sizes
   use tinheader
   use atmlst
   use atoms
   use domdec ,only: nbloc
   use mpole
   implicit none
   integer k,ii,iipole,iglob,i,nloop_
   integer ia,ib,ic,id
   real(t_p) xad,yad,zad
   real(t_p) xbd,ybd,zbd
   real(t_p) xcd,ycd,zcd
   real(t_p) c1,c2,c3,vol
   logical check
   logical,intent(in):: init
!
   if (init) then
      nloop_ = npole
   else
      nloop_ = npolelocnl
   end if
!
!     loop over multipole sites testing for chirality inversion
!
   do ii = 1, nloop_
      if (init) then
         iipole = ii
      else
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc  (iglob)
         if ((i.le.0).or.(i.gt.nbloc)) cycle
      end if
      check = .true.
      if (polaxe(iipole) .ne. 'Z-then-X')  check = .false.
      if (yaxis(iipole) .eq. 0)  check = .false.
      if (check) then
         k  = yaxis(iipole)
         ia = ipole(iipole)
         ib = zaxis(iipole)
         ic = xaxis(iipole)
         id = abs(k)
!
!     compute the signed parallelpiped volume at chiral site
!
         xad = x(ia) - x(id)
         yad = y(ia) - y(id)
         zad = z(ia) - z(id)
         xbd = x(ib) - x(id)
         ybd = y(ib) - y(id)
         zbd = z(ib) - z(id)
         xcd = x(ic) - x(id)
         ycd = y(ic) - y(id)
         zcd = z(ic) - z(id)
         c1  = ybd*zcd - zbd*ycd
         c2  = ycd*zad - zcd*yad
         c3  = yad*zbd - zad*ybd
         vol = xad*c1 + xbd*c2 + xcd*c3
!
!     invert atomic multipole components involving the y-axis
!
         if (k.lt.0.and.vol.gt.0.0_ti_p .or.&
             k.gt.0.and.vol.lt.0.0_ti_p) then
            yaxis(iipole) = -k
            pole(3,iipole) = -pole(3,iipole)
            pole(6,iipole) = -pole(6,iipole)
            pole(8,iipole) = -pole(8,iipole)
            pole(10,iipole) = -pole(10,iipole)
            pole(12,iipole) = -pole(12,iipole)
         end if
      end if
   end do
end
