#ifndef IMAGE_INC
#define IMAGE_INC

#include "tinker_cudart.h"

#ifndef ORTHOGONAL_BOX_SHAPE_ONLY
#define INC_OCTAHEDRON_BOX_SHAPE
#endif

      M_subroutine
     &            image_inl(xr,yr,zr)
!$acc routine seq
      use tinheader,only:ti_p,prec_eps
#ifdef TINKER_CUF
      use utilcu
#else
      use boxes ,only: octahedron,box34
      use cell
#endif
      implicit none
      real(ti_p),intent(inout):: xr,yr,zr

      if (f_abs(xr) .gt. xcell2)
     &   xr  = xr - f_sign(xcell,xr)
     &             *f_floor((f_abs(xr)-xcell2)*i_xcell + 1.0_ti_p)
      if (f_abs(yr) .gt. ycell2)
     &   yr  = yr - f_sign(ycell,yr)
     &             *f_floor((f_abs(yr)-ycell2)*i_ycell + 1.0_ti_p)
      if (f_abs(zr) .gt. zcell2)
     &   zr  = zr - f_sign(zcell,zr)
     &             *f_floor((f_abs(zr)-zcell2)*i_zcell + 1.0_ti_p)

#ifdef INC_OCTAHEDRON_BOX_SHAPE
      if (octahedron) then
         if (f_abs(xr)+f_abs(yr)+f_abs(zr) .gt. box34) then
            xr = xr - f_sign(xcell2,xr)
            yr = yr - f_sign(ycell2,yr)
            zr = zr - f_sign(zcell2,zr)
         end if
      end if
#endif
      end subroutine image_inl


      M_subroutine
     &            imagem_inl(xr,yr,zr,xbox2,ybox2,zbox2,box3)
!$acc routine seq
      use tinheader,only:ti_p,re_p,prec_eps
#ifdef TINKER_CUF
      use utilcu
#else
      use boxes ,only: octahedron,box34
      use cell
#endif
      implicit none
      real(re_p),intent(in)::xbox2,ybox2,zbox2,box3
      real(re_p),intent(inout):: xr,yr,zr

      if (abs(xr) .gt. xbox2) xr = xr - sign(2*xbox2,xr)
      if (abs(yr) .gt. ybox2) yr = yr - sign(2*ybox2,yr)
      if (abs(zr) .gt. zbox2) zr = zr - sign(2*zbox2,zr)

#ifdef INC_OCTAHEDRON_BOX_SHAPE
      if (octahedron) then
         if (abs(xr)+abs(yr)+abs(zr).gt.box3) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
#endif
      end subroutine imagem_inl


      M_subroutine
     &            image_orth_inl(xr,yr,zr)
!$acc routine seq
      use tinheader,only:ti_p
#ifdef TINKER_CUF
      use utilcu
#else
      use cell
#endif
      implicit none
      real(ti_p),intent(inout):: xr,yr,zr

      if (f_abs(xr) .gt. xcell2)  xr = xr - f_sign(xcell,xr)
      if (f_abs(yr) .gt. ycell2)  yr = yr - f_sign(ycell,yr)
      if (f_abs(zr) .gt. zcell2)  zr = zr - f_sign(zcell,zr)

      end subroutine image_orth_inl

      ! Complete Octahedron box image
      ! after inboxing in the cubic box
      !--------------------------------
      M_subroutine
     &            imageOctaComplete(xr,yr,zr)
!$acc routine seq
      use tinheader,only:ti_p,prec_eps
#ifdef TINKER_CUF
      use utilcu
#else
      use boxes ,only: octahedron,box34
      use cell
#endif
      real(t_p),intent(inout)::xr,yr,zr
      if (f_abs(xr)+f_abs(yr)+f_abs(zr) .gt. box34) then
         xr = xr - f_sign(xcell2,xr)
         yr = yr - f_sign(ycell2,yr)
         zr = zr - f_sign(zcell2,zr)
      end if
      end subroutine imageOctaComplete


      !TODO Extract erfcf_hastings from this file
      real(4) M_function
     &                  erfcf_hastings(x) result(res)
      implicit none
      real(4),intent(in) ::x
      !real(4),intent(out)::res
      real(4) t
 
      !float exp2a = f_expf(-x * x)
      t   = (1.0 + 0.3275911 * x)**(-1)
      res =  (0.254829592 + (-0.284496736 + 
     &       (1.421413741 + 
     &       (-1.453152027 + 1.061405429 * t) * t) * t) * t) * t !*exp2a
      end function
 

      subroutine image1d_inl(val,vcell,vcell2)
!$acc routine seq
      use tinheader,only:ti_p,prec_eps
#ifdef TINKER_CUF
      use utilcu   ,only:f_floor,f_sign
#endif
      implicit none
      real(ti_p),intent(inout):: val
      real(ti_p),intent(in) :: vcell,vcell2

      if (abs(val) .gt. vcell2)
     &   val  = val - f_sign(vcell,val)
     &               *f_floor((abs(val)-vcell2)/vcell + 1.0_ti_p)
      end subroutine image1d_inl

#endif
