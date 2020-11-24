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
      subroutine image (xr,yr,zr)
      use sizes
      use boxes
      use cell
      implicit none
      real*8 xr,yr,zr
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
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

c     #################################################################
c     ##                                                             ##
c     ##  subroutine image2  --  compute the minimum image distance  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c
      subroutine image2 (xr,yr,zr)
      use sizes
      use boxes
      use cell
      implicit none
      real*8 xr,yr,zr,cel
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         if (abs(xr) .gt. xcell2) then
            cel      = sign(xcell,xr)
            xr = xr - cel*floor(abs(xr)/xcell)
            if ((abs(xr)) .gt. xcell2)
     &         xr = xr-sign(xcell,xr)
         end if
         if (abs(yr) .gt. xcell2) then
            cel      = sign(ycell,yr)
            yr = yr - cel*floor(abs(yr)/ycell)
            if ((abs(yr)) .gt. ycell2)
     &         yr = yr-sign(ycell,yr)
         end if
         if (abs(zr) .gt. zcell2) then
            cel      = sign(zcell,zr)
            zr = zr - cel*floor(abs(zr)/zcell)
            if ((abs(zr)) .gt. zcell2)
     &         zr = zr-sign(zcell,zr)
         end if
c         do while (abs(xr) .gt. xcell2)
c            xr = xr - sign(xcell,xr)
c         end do
c         do while (abs(yr) .gt. ycell2)
c            yr = yr - sign(ycell,yr)
c         end do
c         do while (abs(zr) .gt. zcell2)
c            zr = zr - sign(zcell,zr)
c         end do
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
c
      subroutine imagevec (pos,n)
      use sizes
      use cell
      implicit none
      integer n,i
      real*8 pos(3,n)
      real*8 xr,yr,zr
c
      do i = 1, n
        xr = pos(1,i)
        yr = pos(2,i)
        zr = pos(3,i)
c
c       for orthogonal lattice, find the desired image directly
c
c       if (orthogonal) then
        do while (abs(xr) .gt. xcell2)
           xr = xr - sign(xcell,xr)
        end do
        do while (abs(yr) .gt. ycell2)
           yr = yr - sign(ycell,yr)
        end do
        do while (abs(zr) .gt. zcell2)
           zr = zr - sign(zcell,zr)
        end do
c       end if
        pos(1,i) = xr
        pos(2,i) = yr 
        pos(3,i) = zr
      end do
c      do while (any(abs(pos(1,1:n)).gt.xcell2))
c         where (    abs(pos(1,1:n)).gt.xcell2)
c            pos(1,1:n) = pos(1,1:n) -sign(xcell,pos(1,1:n))
c         end where
c      enddo
c      do while (any(abs(pos(2,1:n)).gt.ycell2))
c         where (    abs(pos(2,1:n)).gt.ycell2)
c            pos(2,1:n) = pos(2,1:n) -sign(ycell,pos(2,1:n))
c         end where
c      enddo
c      do while (any(abs(pos(3,1:n)).gt.zcell2))
c         where (    abs(pos(3,1:n)).gt.zcell2)
c            pos(3,1:n) = pos(3,1:n) -sign(zcell,pos(3,1:n))
c         end where
c      enddo
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine imagevec2  --  compute the minimum image distance      ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "imagevec" takes the components of pairwise distances between
c     two points in a periodic box and converts to the components
c     of the minimum image distances 
c
      subroutine imagevec2_acc(pos2,n,size_pos)
      use cell
      
      implicit none
      integer n,size_pos
      real*8 pos2(size_pos,3)
      integer j
      real*8 cel

      do j=1,n
         if (abs(pos2(j,1)) .gt. xcell2) then
            cel      = sign(xcell,pos2(j,1))
            pos2(j,1) = pos2(j,1) -cel*floor(abs(pos2(j,1))/xcell)
            if (abs(pos2(j,1)) .gt. xcell2)
     &         pos2(j,1) = pos2(j,1)-sign(xcell,pos2(j,1))
         end if
      end do
      do j=1,n
         if (abs(pos2(j,2)) .gt. ycell2) then
            cel      = sign(ycell,pos2(j,2))
            pos2(j,2) = pos2(j,2) -cel*floor(abs(pos2(j,2))/ycell)
            if (abs(pos2(j,2)) .gt. ycell2)
     &         pos2(j,2) = pos2(j,2)-sign(ycell,pos2(j,2))
         end if
      end do
      do j=1,n
         if (abs(pos2(j,3)) .gt. zcell2) then
            cel      = sign(zcell,pos2(j,3))
            pos2(j,3) = pos2(j,3) -cel*floor(abs(pos2(j,3))/zcell)
            if (abs(pos2(j,3)) .gt. zcell2)
     &         pos2(j,3) = pos2(j,3)-sign(zcell,pos2(j,3))
         end if
      end do
      end
