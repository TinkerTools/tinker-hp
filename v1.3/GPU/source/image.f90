!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine image  --  compute the minimum image distance  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "image" takes the components of pairwise distance between
!     two points in a periodic box and converts to the components
!     of the minimum image distance
!
!
#include "tinker_precision.h"
subroutine image (xr,yr,zr)
   use sizes
   use boxes
   use cell
   implicit none
   real(t_p) xr,yr,zr,cel
!
!
!     for orthogonal lattice, find the desired image directly
!
   if (orthogonal) then
      if (abs(xr) .gt. xcell2) then
         cel = sign(xcell,xr)
         xr  = xr - cel*floor((abs(xr)+xcell2)/xcell)
      end if
      if (abs(yr) .gt. ycell2) then
         cel = sign(ycell,yr)
         yr  = yr - cel*floor((abs(yr)+ycell2)/ycell)
      end if
      if (abs(zr) .gt. zcell2) then
         cel = sign(zcell,zr)
         zr  = zr - cel*floor((abs(zr)+zcell2)/zcell)
      end if
!
!     for monoclinic lattice, convert "xr" and "zr" to
!     fractional coordinates, find desired image and then
!     translate fractional coordinates back to Cartesian
!
   else if (monoclinic) then
      zr = zr / beta_sin
      xr = xr - zr*beta_cos
      do while (abs(xr) .gt. xcell2)
         xr = xr - sign(xcell,xr)
      end do
      do while (abs(yr) .gt. ycell2)
         yr = yr - sign(ycell,yr)
      end do
      do while (abs(zr) .gt. zcell2)
         zr = zr - sign(zcell,zr)
      end do
      xr = xr + zr*beta_cos
      zr = zr * beta_sin
!
!     for triclinic lattice, convert pairwise components to
!     fractional coordinates, find desired image and then
!     translate fractional coordinates back to Cartesian
!
   else if (triclinic) then
      zr = zr / gamma_term
      yr = (yr - zr*beta_term) / gamma_sin
      xr = xr - yr*gamma_cos - zr*beta_cos
      do while (abs(xr) .gt. xcell2)
         xr = xr - sign(xcell,xr)
      end do
      do while (abs(yr) .gt. ycell2)
         yr = yr - sign(ycell,yr)
      end do
      do while (abs(zr) .gt. zcell2)
         zr = zr - sign(zcell,zr)
      end do
      xr = xr + yr*gamma_cos + zr*beta_cos
      yr = yr*gamma_sin + zr*beta_term
      zr = zr * gamma_term
!
!     for truncated octahedron, use orthogonal box equations,
!     then perform extra tests to remove corner pieces
!
   else if (octahedron) then
      do while (abs(xr) .gt. xcell2)
         xr = xr - sign(xcell,xr)
      end do
      do while (abs(yr) .gt. ycell2)
         yr = yr - sign(ycell,yr)
      end do
      do while (abs(zr) .gt. zcell2)
         zr = zr - sign(zcell,zr)
      end do
      if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
         xr = xr - sign(xcell2,xr)
         yr = yr - sign(ycell2,yr)
         zr = zr - sign(zcell2,zr)
      end if
   end if
   return
end
!
subroutine image_acc(xr,yr,zr)
!$acc routine seq
   use boxes,only: orthogonal,octahedron&
      &,box34
   use cell
   use tinheader ,only: ti_p
   implicit none
   real(t_p) xr,yr,zr
   real(t_p) cel

   if (orthogonal) then
      if (abs(xr) .gt. xcell2) then
         cel = sign(xcell,xr)
         xr  = xr - cel*floor((abs(xr)+xcell2)/xcell)
      end if
      if (abs(yr) .gt. ycell2) then
         cel = sign(ycell,yr)
         yr  = yr - cel*floor((abs(yr)+ycell2)/ycell)
      end if
      if (abs(zr) .gt. zcell2) then
         cel = sign(zcell,zr)
         zr  = zr - cel*floor((abs(zr)+zcell2)/zcell)
      end if
   else if (octahedron) then
      if (abs(xr) .gt. xcell2)&
         &xr  = xr - sign(xcell,xr)&
         &*floor((abs(xr)-xcell2)*i_xcell + 1.0_ti_p)
      if (abs(yr) .gt. ycell2)&
         &yr  = yr - sign(ycell,yr)&
         &*floor((abs(yr)-ycell2)*i_ycell + 1.0_ti_p)
      if (abs(zr) .gt. zcell2)&
         &zr  = zr - sign(zcell,zr)&
         &*floor((abs(zr)-zcell2)*i_zcell + 1.0_ti_p)
      if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
         xr = xr - sign(xcell2,xr)
         yr = yr - sign(ycell2,yr)
         zr = zr - sign(zcell2,zr)
      end if
   end if
end
!
subroutine image1d(val,vcell,vcell2)
!$acc routine seq
   use tinheader,only:ti_p
   implicit none
   real(ti_p),intent(inout):: val
   real(ti_p),intent(in) :: vcell,vcell2

   if (abs(val) .gt. vcell2)&
      &val  = val - sign(vcell,val)*floor((abs(val)+vcell2)/vcell)
end subroutine
!
!    subroutine imagefrac: takes a vector in fractional coordinate and
!    get the associated minimum image convention vector in fractional coordinate
!

subroutine imagefrac(ar,br,cr)
   use boxes
   use cell
   implicit none
   real(t_p) ar,br,cr, cel
   if (abs(ar) .gt. xcell2) &
      ar = ar - sign(xcell,ar)*floor((abs(ar)+xcell2)/xcell)
   if (abs(br) .gt. ycell2) &
      br = br - sign(ycell,br)*floor((abs(br)+ycell2)/ycell)
   if (abs(cr) .gt. zcell2) &
      cr = cr - sign(zcell,cr)*floor((abs(cr)+zcell2)/zcell)
end
!
!    subroutine imagecell: takes a vector in "cell" coordinate and
!    get the associated minimum image convention vector in "cell" coordinate
!

subroutine imagecell(n_ar,n_br,n_cr,na,nb,nc)
   use boxes
   use cell
   implicit none
   integer n_ar,n_br,n_cr,na,nb,nc
   do while (abs(n_ar) .gt. int(na/2))
      n_ar = n_ar - sign(na,n_ar)
   end do
   do while (abs(n_br) .gt. int(nb/2))
      n_br = n_br - sign(nb,n_br)
   end do
   do while (abs(n_cr) .gt. int(nc/2))
      n_cr = n_cr - sign(nc,n_cr)
   end do
end
!
!     subroutine ctfvec: takes a vector in cartesian coordinates and put
!     it in fractional coordinates
!
subroutine ctfvec(xr,yr,zr,ar,br,cr)
use boxes
   implicit none
   real(t_p) xr,yr,zr
   real(t_p) ar,br,cr
   ar = ctfmat(1,1)*xr + ctfmat(1,2)*yr + ctfmat(1,3)*zr
   ar = xbox*ar
   br = ctfmat(2,1)*xr + ctfmat(2,2)*yr + ctfmat(2,3)*zr
   br = ybox*br
   cr = ctfmat(3,1)*xr + ctfmat(3,2)*yr + ctfmat(3,3)*zr
   cr = zbox*cr
end

!subroutine ftcvec!: takes a vector in fractional coordinates and put it in cartesian coordinates

subroutine ftcvec(ar,br,cr,xr,yr,zr)
   use boxes
   implicit none
   real(t_p) xr,yr,zr
   real(t_p) ar,br,cr
   real(t_p) artemp,brtemp,crtemp
   artemp = ar/xbox
   brtemp = br/ybox
   crtemp = cr/zbox
   xr = ftcmat(1,1)*artemp + ftcmat(2,1)*brtemp + ftcmat(3,1)*crtemp
   yr = ftcmat(1,2)*artemp + ftcmat(2,2)*brtemp + ftcmat(3,2)*crtemp
   zr = ftcmat(1,3)*artemp + ftcmat(2,3)*brtemp + ftcmat(3,3)*crtemp
end
