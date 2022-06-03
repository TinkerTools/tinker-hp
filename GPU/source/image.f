c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  compute the minimum image distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c
#include "tinker_precision.h"
      subroutine image (xr,yr,zr)
      use sizes
      use boxes
      use cell
      implicit none
      real(t_p) xr,yr,zr,cel
c
c
c     for orthogonal lattice, find the desired image directly
c
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
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
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
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
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
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
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
c
      subroutine image_acc(xr,yr,zr) 
!$acc routine seq
      use boxes,only: orthogonal,octahedron
     &         ,box34
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
         if (abs(xr) .gt. xcell2)
     &      xr  = xr - sign(xcell,xr)
     &                *floor((abs(xr)-xcell2)*i_xcell + 1.0_ti_p)
         if (abs(yr) .gt. ycell2)
     &      yr  = yr - sign(ycell,yr)
     &                *floor((abs(yr)-ycell2)*i_ycell + 1.0_ti_p)
         if (abs(zr) .gt. zcell2)
     &      zr  = zr - sign(zcell,zr)
     &                *floor((abs(zr)-zcell2)*i_zcell + 1.0_ti_p)
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xcell2,xr)
            yr = yr - sign(ycell2,yr)
            zr = zr - sign(zcell2,zr)
         end if
      end if
      end
c
      subroutine image1d(val,vcell,vcell2)
!$acc routine seq
      use tinheader,only:ti_p
      implicit none
      real(ti_p),intent(inout):: val
      real(ti_p),intent(in) :: vcell,vcell2

      if (abs(val) .gt. vcell2)
     &   val  = val - sign(vcell,val)*floor((abs(val)+vcell2)/vcell)
      end subroutine
c
cold  subroutine imagevec (pos,n)
cold  use sizes
cold  use cell
cold  implicit none
cold  integer n,i
cold  real(t_p) pos(3,n)
cold  real(t_p) xr,yr,zr
cold
cold  do i = 1, n
cold    xr = pos(1,i)
cold    yr = pos(2,i)
cold    zr = pos(3,i)
cold
cold    for orthogonal lattice, find the desired image directly
cold
cold    if (orthogonal) then
cold    do while (abs(xr) .gt. xcell2)
cold       xr = xr - sign(xcell,xr)
cold    end do
cold    do while (abs(yr) .gt. ycell2)
cold       yr = yr - sign(ycell,yr)
cold    end do
cold    do while (abs(zr) .gt. zcell2)
cold       zr = zr - sign(zcell,zr)
cold    end do
cold    end if
cold    pos(1,i) = xr
cold    pos(2,i) = yr 
cold    pos(3,i) = zr
cold  end do
cold   do while (any(abs(pos(1,1:n)).gt.xcell2))
cold      where (    abs(pos(1,1:n)).gt.xcell2)
cold         pos(1,1:n) = pos(1,1:n) -sign(xcell,pos(1,1:n))
cold      end where
cold   enddo
cold   do while (any(abs(pos(2,1:n)).gt.ycell2))
cold      where (    abs(pos(2,1:n)).gt.ycell2)
cold         pos(2,1:n) = pos(2,1:n) -sign(ycell,pos(2,1:n))
cold      end where
cold   enddo
cold   do while (any(abs(pos(3,1:n)).gt.zcell2))
cold      where (    abs(pos(3,1:n)).gt.zcell2)
cold         pos(3,1:n) = pos(3,1:n) -sign(zcell,pos(3,1:n))
cold      end where
cold   enddo
cold  return
cold  end
