!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine eurey1  --  bond stretch energy & derivatives  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "eurey1" calculates the Urey-Bradley interaction energy and
!     its first derivatives with respect to Cartesian coordinates
!
!
#include "tinker_macro.h"
subroutine eurey1
   use atmlst
   use atoms
   use bound
   use deriv
   use domdec
   use energi
   use group
   use tinheader ,only:ti_p,re_p
   use urey
   use urypot
   use usage
   use virial
   implicit none
   integer i,ia,ic,iurey
   integer ialoc,icloc
   real(t_p) e,de,ideal,force
   real(t_p) dt,dt2,deddt
   real(t_p) dedx,dedy,dedz
   real(t_p) vxx,vyy,vzz
   real(t_p) vyx,vzx,vzy
   real(t_p) xac,yac,zac,rac
   real(t_p) fgrp
   integer iga,igc
   logical proceed
!
!$acc update host(deub,vir)
!
!     zero out the Urey-Bradley energy and first derivatives
!
   eub = 0.0_ti_p
!
!     calculate the Urey-Bradley 1-3 energy and first derivatives
!
   do iurey = 1, nureyloc
      i = ureyglob(iurey)
      ia = iury(1,i)
      ic = iury(3,i)
      ialoc = loc(ia)
      icloc = loc(ic)
      ideal = ul(i)
      force = uk(i)
!
!     decide whether to compute the current interaction
!
      proceed = .true.
      if (proceed)  proceed = (use(ia) .or. use(ic))
!
!     compute the value of the 1-3 distance deviation
!
      if (proceed) then
         xac = x(ia) - x(ic)
         yac = y(ia) - y(ic)
         zac = z(ia) - z(ic)
         if (use_polymer)  call image (xac,yac,zac)
         rac = sqrt(xac*xac + yac*yac + zac*zac)
         dt = rac - ideal
         dt2 = dt * dt
         e = ureyunit * force * dt2 * (1.0_ti_p+cury*dt+qury*dt2)
         deddt = 2.0_ti_p * ureyunit * force * dt&
            &* (1.0_ti_p+1.5_ti_p*cury*dt+2.0_ti_p*qury*dt2)

         if(use_group) then
            iga=grplist(ia)
            igc=grplist(ic)
            fgrp = wgrp(iga+1,igc+1)
            e = e*fgrp
            deddt = deddt*fgrp
         endif
!
!     compute chain rule terms needed for derivatives
!
         de = deddt / rac
         dedx = de * xac
         dedy = de * yac
         dedz = de * zac
!
!     increment the total Urey-Bradley energy and first derivatives
!
         eub = eub + e
         deub(1,ialoc) = deub(1,ialoc) + dedx
         deub(2,ialoc) = deub(2,ialoc) + dedy
         deub(3,ialoc) = deub(3,ialoc) + dedz
!
         deub(1,icloc) = deub(1,icloc) - dedx
         deub(2,icloc) = deub(2,icloc) - dedy
         deub(3,icloc) = deub(3,icloc) - dedz
!
!     increment the internal virial tensor components
!
         vxx = xac * dedx
         vyx = yac * dedx
         vzx = zac * dedx
         vyy = yac * dedy
         vzy = zac * dedy
         vzz = zac * dedz
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vyx
         vir(3,1) = vir(3,1) + vzx
         vir(1,2) = vir(1,2) + vyx
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vzy
         vir(1,3) = vir(1,3) + vzx
         vir(2,3) = vir(2,3) + vzy
         vir(3,3) = vir(3,3) + vzz
      end if
   end do
!$acc update device(deub,vir)
   return
end
