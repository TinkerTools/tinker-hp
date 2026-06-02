!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##                                                             ##
!     ##  subroutine image2  --  compute the minimum image distance  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "image" takes the components of pairwise distance between
!     two points in a periodic box and converts to the components
!     of the minimum image distance
!
!
!> @brief 
!> takes the components of pairwise distance between
!> two points in a periodic box and converts to the components
!> of the minimum image distance
!> @param[in] xr: first component of the distance vector
!> @param[in] yr: second component of the distance vector
!> @param[in] zr: third component of the distance vector
subroutine image (xr,yr,zr)
   use sizes
   use boxes
   use cell
   implicit none
   real*8 xr,yr,zr,cel
!
!
!     for orthogonal lattice, find the desired image directly
!
   if (orthogonal) then
      if (abs(xr) .gt. xcell2) then
         cel      = sign(xcell,xr)
         xr = xr - cel*floor(abs(xr)/xcell)
         if ((abs(xr)) .gt. xcell2)&
         &xr = xr-sign(xcell,xr)
      end if
      if (abs(yr) .gt. ycell2) then
         cel      = sign(ycell,yr)
         yr = yr - cel*floor(abs(yr)/ycell)
         if ((abs(yr)) .gt. ycell2)&
         &yr = yr-sign(ycell,yr)
      end if
      if (abs(zr) .gt. zcell2) then
         cel      = sign(zcell,zr)
         zr = zr - cel*floor(abs(zr)/zcell)
         if ((abs(zr)) .gt. zcell2)&
         &zr = zr-sign(zcell,zr)
      end if
!
   else if (monoclinic) then
      zr = zr / beta_sin
      xr = xr - zr*beta_cos
      if (abs(xr) .gt. xcell2) then
         cel      = sign(xcell,xr)
         xr = xr - cel*floor(abs(xr)/xcell)
         if ((abs(xr)) .gt. xcell2)&
         &xr = xr-sign(xcell,xr)
      end if
      if (abs(yr) .gt. ycell2) then
         cel      = sign(ycell,yr)
         yr = yr - cel*floor(abs(yr)/ycell)
         if ((abs(yr)) .gt. ycell2)&
         &yr = yr-sign(ycell,yr)
      end if
      if (abs(zr) .gt. zcell2) then
         cel      = sign(zcell,zr)
         zr = zr - cel*floor(abs(zr)/zcell)
         if ((abs(zr)) .gt. zcell2)&
         &zr = zr-sign(zcell,zr)
      end if
!      
!      do while (abs(xr) .gt. xcell2)
!         xr = xr - sign(xcell,xr)
!      end do
!      do while (abs(yr) .gt. ycell2)
!         yr = yr - sign(ycell,yr)
!      end do
!      do while (abs(zr) .gt. zcell2)
!         zr = zr - sign(zcell,zr)
!      end do
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
      if (abs(xr) .gt. xcell2) then
         cel      = sign(xcell,xr)
         xr = xr - cel*floor(abs(xr)/xcell)
         if ((abs(xr)) .gt. xcell2)&
         &xr = xr-sign(xcell,xr)
      end if
      if (abs(yr) .gt. ycell2) then
         cel      = sign(ycell,yr)
         yr = yr - cel*floor(abs(yr)/ycell)
         if ((abs(yr)) .gt. ycell2)&
         &yr = yr-sign(ycell,yr)
      end if
      if (abs(zr) .gt. zcell2) then
         cel      = sign(zcell,zr)
         zr = zr - cel*floor(abs(zr)/zcell)
         if ((abs(zr)) .gt. zcell2)&
         &zr = zr-sign(zcell,zr)
      end if
!      do while (abs(xr) .gt. xcell2)
!         xr = xr - sign(xcell,xr)
!      end do
!      do while (abs(yr) .gt. ycell2)
!         yr = yr - sign(ycell,yr)
!      end do
!      do while (abs(zr) .gt. zcell2)
!         zr = zr - sign(zcell,zr)
!      end do
      xr = xr + yr*gamma_cos + zr*beta_cos
      yr = yr*gamma_sin + zr*beta_term
      zr = zr * gamma_term
!
!     for truncated octahedron, use orthogonal box equations,
!     then perform extra tests to remove corner pieces
!
   else if (octahedron) then
      do while (abs(xr) .gt. xbox2)
         xr = xr - sign(xbox,xr)
      end do
      do while (abs(yr) .gt. ybox2)
         yr = yr - sign(ybox,yr)
      end do
      do while (abs(zr) .gt. zbox2)
         zr = zr - sign(zbox,zr)
      end do
      if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
         xr = xr - sign(xbox2,xr)
         yr = yr - sign(ybox2,yr)
         zr = zr - sign(zbox2,zr)
      end if
   end if
   return
end
!
!> @brief 
!> takes the components of pairwise distance between
!> two vector of points in a periodic box and converts to the components
!> of the minimum image distance
!> @param[in] pos: vector of pairwise distance
!> @param[in] n: size of the vector
subroutine imagevec (pos,n)
   use sizes
   use cell
   implicit none
   integer n,i
   real*8 pos(3,n)
   real*8 xr,yr,zr
!
   do i = 1, n
      xr = pos(1,i)
      yr = pos(2,i)
      zr = pos(3,i)
!
!       for orthogonal lattice, find the desired image directly
!
!       if (orthogonal) then
      do while (abs(xr) .gt. xcell2)
         xr = xr - sign(xcell,xr)
      end do
      do while (abs(yr) .gt. ycell2)
         yr = yr - sign(ycell,yr)
      end do
      do while (abs(zr) .gt. zcell2)
         zr = zr - sign(zcell,zr)
      end do
!       end if
      pos(1,i) = xr
      pos(2,i) = yr
      pos(3,i) = zr
   end do
!      do while (any(abs(pos(1,1:n)).gt.xcell2))
!         where (    abs(pos(1,1:n)).gt.xcell2)
!            pos(1,1:n) = pos(1,1:n) -sign(xcell,pos(1,1:n))
!         end where
!      enddo
!      do while (any(abs(pos(2,1:n)).gt.ycell2))
!         where (    abs(pos(2,1:n)).gt.ycell2)
!            pos(2,1:n) = pos(2,1:n) -sign(ycell,pos(2,1:n))
!         end where
!      enddo
!      do while (any(abs(pos(3,1:n)).gt.zcell2))
!         where (    abs(pos(3,1:n)).gt.zcell2)
!            pos(3,1:n) = pos(3,1:n) -sign(zcell,pos(3,1:n))
!         end where
!      enddo
   return
end
!
!    subroutine imagefrac: takes a vector in fractional coordinate and
!    get the associated minimum image convention vector in fractional coordinate
!

!> @brief 
!> takes a vector in fractional coordinate and
!> get the associated minimum image convention vector in fractional coordinate
!> @param[in] ar: first component of the distance vector
!> @param[in] br: second component of the distance vector
!> @param[in] cr: third component of the distance vector
subroutine imagefrac(ar,br,cr)
use boxes
use cell
implicit none
real*8 ar,br,cr,cel
if (abs(ar) .gt. xcell2) then
   cel      = sign(xcell,ar)
   ar = ar - cel*floor(abs(ar)/xcell)
   if ((abs(ar)) .gt. xcell2)&
   &ar = ar-sign(xcell,ar)
end if
if (abs(br) .gt. ycell2) then
   cel      = sign(ycell,br)
   br = br - cel*floor(abs(br)/ycell)
   if ((abs(br)) .gt. ycell2)&
   &br = br-sign(ycell,br)
end if
if (abs(cr) .gt. zcell2) then
   cel      = sign(zcell,cr)
   cr = cr - cel*floor(abs(cr)/zcell)
   if ((abs(cr)) .gt. zcell2)&
   &cr = cr-sign(zcell,cr)
end if
return
end
!
!     subroutine ctfvec: takes a vector in cartesian coordinates and put
!     it in fractional coordinates
!
!> @brief 
!> takes a vector in cartesian coordinates and put
!> it in fractional coordinates
!> @param[in] xr: first cartesian component of the distance vector
!> @param[in] yr: second cartesian component of the distance vector
!> @param[in] zr: third cartesian component of the distance vector
!> @param[in] ar: first fractional component of the distance vector
!> @param[in] br: second fractional component of the distance vector
!> @param[in] cr: third fractional component of the distance vector
subroutine ctfvec(xr,yr,zr,ar,br,cr)
use boxes
implicit none
real*8 xr,yr,zr
real*8 ar,br,cr
ar = ctfmat(1,1)*xr + ctfmat(1,2)*yr + ctfmat(1,3)*zr
ar = xbox*ar
br = ctfmat(2,1)*xr + ctfmat(2,2)*yr + ctfmat(2,3)*zr
br = ybox*br
cr = ctfmat(3,1)*xr + ctfmat(3,2)*yr + ctfmat(3,3)*zr
cr = zbox*cr
return
end

!subroutine ftcvec!: takes a vector in fractional coordinates and put it in cartesian coordinates

!> @brief 
!> takes a vector in fractional coordinates and put
!> it in cartesian coordinates
!> @param[in] xr: first cartesian component of the distance vector
!> @param[in] yr: second cartesian component of the distance vector
!> @param[in] zr: third cartesian component of the distance vector
!> @param[in] ar: first fractional component of the distance vector
!> @param[in] br: second fractional component of the distance vector
!> @param[in] cr: third fractional component of the distance vector
subroutine ftcvec(ar,br,cr,xr,yr,zr)
use boxes
implicit none
real*8 xr,yr,zr
real*8 ar,br,cr
real*8 artemp,brtemp,crtemp
artemp = ar/xbox
brtemp = br/ybox
crtemp = cr/zbox
xr = ftcmat(1,1)*artemp + ftcmat(2,1)*brtemp + ftcmat(3,1)*crtemp
yr = ftcmat(1,2)*artemp + ftcmat(2,2)*brtemp + ftcmat(3,2)*crtemp
zr = ftcmat(1,3)*artemp + ftcmat(2,3)*brtemp + ftcmat(3,3)*crtemp
return
end
