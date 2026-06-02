#include "tinker_cudart.h"

#ifndef ORTHOGONAL_BOX_SHAPE_ONLY
#define INC_ALL_BOX_SHAPE
#endif

M_subroutine&
   midpointimage_inl(xk,yk,zk,xr,yr,zr)
!$acc routine
use tinheader,only : ti_p,ti_eps,prec_eps
#ifdef TINKER_CUF
use utilcu
#else
use boxes
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
#ifdef INC_ALL_BOX_SHAPE
if (monoclinic.or.triclinic) then
   call ctf_inl(xr,yr,zr)
   call ctf_inl(xk,yk,zk)
end if
#endif

!call image_orth_inl(xr,yr,zr)
if (f_abs(xr).gt. xcell2)&
   xr = xr - f_sign(xcell,xr) *f_floor((f_abs(xr)-xcell2)*i_xcell+1.0_ti_p)
if (f_abs(yr).gt. ycell2)&
   yr = yr - f_sign(ycell,yr) *f_floor((f_abs(yr)-ycell2)*i_ycell+1.0_ti_p)
if (f_abs(zr).gt. zcell2)&
   zr = zr - f_sign(zcell,zr) *f_floor((f_abs(zr)-zcell2)*i_zcell+1.0_ti_p)

#ifdef INC_ALL_BOX_SHAPE
if (octahedron) then
   if (f_abs(xr)+f_abs(yr)+f_abs(zr).gt.box34) then
      xr = xr - f_sign(xcell2,xr)
      yr = yr - f_sign(ycell2,yr)
      zr = zr - f_sign(zcell2,zr)
   end if
end if
#endif
!
!     definition of the middle point between i and k atoms
!
xk = xk + 0.5*xr
yk = yk + 0.5*yr
zk = zk + 0.5*zr
!
!     Compute image of the middle point
!
call image_orth_inl(xk,yk,zk)

#ifdef INC_ALL_BOX_SHAPE
if      (monoclinic.or.triclinic) then
   call ftc_inl(xr,yr,zr)
else if (octahedron) then
   if (f_abs(xk)+f_abs(yk)+f_abs(zk).gt.box34) then
      xk = xk - f_sign(xcell2,xk)
      yk = yk - f_sign(ycell2,yk)
      zk = zk - f_sign(zcell2,zk)
   end if
end if
#endif
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
use utilcu
#else
use boxes
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
   use boxes
   use cell
   implicit none
   real(ti_p),value:: xk,yk,zk
   real(ti_p) xr,yr,zr
   logical docompute
!
   docompute = .false.
#ifdef INC_ALL_BOX_SHAPE
if (monoclinic.or.triclinic) then
   call ctf_inl(xr,yr,zr)
   call ctf_inl(xk,yk,zk)
end if
#endif

   call image_orth_inl(xr,yr,zr)
!
!  definition of the middle point between i and k atoms
!
   xk = xk + xr/2
   yk = yk + yr/2
   zk = zk + zr/2
   call image_inl(xk,yk,zk)
#ifdef INC_ALL_BOX_SHAPE
if (monoclinic.or.triclinic) call ftc_inl(xr,yr,zr)
#endif
!
!  Adjust mid point position if necessary
!
   if ((xcell2-f_abs(xk)).lt.eps_cell) xk= xk- f_sign(5*eps_cell,xk)
   if ((ycell2-f_abs(yk)).lt.eps_cell) yk= yk- f_sign(5*eps_cell,yk)
   if ((zcell2-f_abs(zk)).lt.eps_cell) zk= zk- f_sign(5*eps_cell,zk)

   if   ((zk.ge.zbegproc(rank+1)).and.(zk.lt.zendproc(rank+1))&
   &.and.(yk.ge.ybegproc(rank+1)).and.(yk.lt.yendproc(rank+1))&
   &.and.(xk.ge.xbegproc(rank+1)).and.(xk.lt.xendproc(rank+1))) docompute = .true.
end subroutine
#endif
