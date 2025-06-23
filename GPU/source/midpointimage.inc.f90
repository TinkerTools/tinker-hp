#include "tinker_cudart.h"

M_subroutine&
   &midpointimage_inl(xk,yk,zk,xr,yr,zr)
!$acc routine
use tinheader,only : ti_p,ti_eps,prec_eps
#ifdef TINKER_CUF
use utilcu,only:xcell,xcell2,ycell,ycell2,zcell,zcell2&
   &,i_xcell,i_ycell,i_zcell,eps_cell,octahedron&
   &,box34,f_abs,f_floor,f_sign
#else
use boxes  ,only: box34,octahedron
use cell
#endif
implicit none
real(ti_p) xk,yk,zk
real(ti_p) xr,yr,zr
!$acc routine(image_acc)
!
!     inline image compute by hand
!     Much efficient on GPU
!
if (f_abs(xr) .gt. xcell2)&
   &xr  = xr - f_sign(xcell,xr)&
   &*f_floor((f_abs(xr)-xcell2)*i_xcell+1.0_ti_p)
if (f_abs(yr) .gt. ycell2)&
   &yr  = yr - f_sign(ycell,yr)&
   &*f_floor((f_abs(yr)-ycell2)*i_ycell+1.0_ti_p)
if (f_abs(zr) .gt. zcell2)&
   &zr  = zr - f_sign(zcell,zr)&
   &*f_floor((f_abs(zr)-zcell2)*i_zcell+1.0_ti_p)
if (octahedron) then
   if (f_abs(xr)+f_abs(yr)+f_abs(zr).gt.box34) then
      xr = xr - f_sign(xcell2,xr)
      yr = yr - f_sign(ycell2,yr)
      zr = zr - f_sign(zcell2,zr)
   end if
end if
!
!     definition of the middle point between i and k atoms
!
xk = xk + 0.5*xr
yk = yk + 0.5*yr
zk = zk + 0.5*zr
!
!     Compute image of the middle point
!
if (f_abs(xk) .gt. xcell2)&
   &xk  = xk - f_sign(xcell,xk)&
   &*f_floor((f_abs(xk)-xcell2)*i_xcell+1.0_ti_p)
if (f_abs(yk) .gt. ycell2)&
   &yk  = yk - f_sign(ycell,yk)&
   &*f_floor((f_abs(yk)-ycell2)*i_ycell+1.0_ti_p)
if (f_abs(zk) .gt. zcell2)&
   &zk  = zk - f_sign(zcell,zk)&
   &*f_floor((f_abs(zk)-zcell2)*i_zcell+1.0_ti_p)
if (octahedron) then
   if (f_abs(xk)+f_abs(yk)+f_abs(zk).gt.box34) then
      xk = xk - f_sign(xcell2,xk)
      yk = yk - f_sign(ycell2,yk)
      zk = zk - f_sign(zcell2,zk)
   end if
end if
!
!     Adjust mid point position if necessary
!
if ((xcell2-f_abs(xk)).lt.eps_cell) xk= xk- f_sign(5*eps_cell,xk)
if ((ycell2-f_abs(yk)).lt.eps_cell) yk= yk- f_sign(5*eps_cell,yk)
if ((zcell2-f_abs(zk)).lt.eps_cell) zk= zk- f_sign(5*eps_cell,zk)

end subroutine

M_subroutine&
   &midpointimage1_inl(xk,yk,zk,xr,yr,zr)
!$acc routine
use tinheader,only : ti_p,ti_eps,prec_eps
#ifdef TINKER_CUF
use utilcu,only:xcell,xcell2,ycell,ycell2,zcell,zcell2&
   &,i_xcell,i_ycell,i_zcell,eps_cell,octahedron&
   &,box34,f_abs,f_floor,f_sign
#else
use boxes  ,only: box34,octahedron
use cell
#endif
implicit none
real(ti_p) xk,yk,zk
real(ti_p) xr,yr,zr
!$acc routine(image_acc)
!
!     definition of the middle point between i and k atoms
!
xk = xk + 0.5*xr
yk = yk + 0.5*yr
zk = zk + 0.5*zr
!
!     Adjust mid point position if necessary
!
if ((xcell2-f_abs(xk)).lt.eps_cell) xk= xk- f_sign(5*eps_cell,xk)
if ((ycell2-f_abs(yk)).lt.eps_cell) yk= yk- f_sign(5*eps_cell,yk)
if ((zcell2-f_abs(zk)).lt.eps_cell) zk= zk- f_sign(5*eps_cell,zk)
end subroutine

#ifndef TINKER_CUF
subroutine midpoint_inl(xk,yk,zk,xr,yr,zr,docompute)
!$acc routine seq
   use tinheader,only:ti_p,ti_eps
   use domdec,only:xbegproc,ybegproc,zbegproc,&
      &xendproc,yendproc,zendproc,rank
   use cell
   implicit none
   real(ti_p),value:: xk,yk,zk
   real(ti_p) xr,yr,zr
   logical docompute
!
   docompute = .false.
   call image_inl(xr,yr,zr)
!
!     definition of the middle point between i and k atoms
!
   xk = xk + xr/2
   yk = yk + yr/2
   zk = zk + zr/2
   call image_inl(xk,yk,zk)
!
!     Adjust mid point position if necessary
!
   if ((xcell2-f_abs(xk)).lt.eps_cell) xk= xk- f_sign(5*eps_cell,xk)
   if ((ycell2-f_abs(yk)).lt.eps_cell) yk= yk- f_sign(5*eps_cell,yk)
   if ((zcell2-f_abs(zk)).lt.eps_cell) zk= zk- f_sign(5*eps_cell,zk)

   if   ((zk.ge.zbegproc(rank+1)).and.(zk.lt.zendproc(rank+1))&
      &.and.(yk.ge.ybegproc(rank+1)).and.(yk.lt.yendproc(rank+1))&
      &.and.(xk.ge.xbegproc(rank+1)).and.(xk.lt.xendproc(rank+1)))&
      &then
      docompute = .true.
   end if
end
#endif
