c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine imagevec2  --  compute the minimum image distance  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "imagevec" takes the components of pairwise distances between
c     two points in a periodic box and converts to the components
c     of the minimum image distances
c
      subroutine imagevec2 (pos2,n)
      use sizes
      use boxes
      use cell

      implicit none
      integer n,i
      real*8 pos2(n,3)
      !DIR$ ASSUME_ALIGNED pos2:64
       do while (any(abs(pos2(1:n,1)).gt.xcell2))
          where (    abs(pos2(1:n,1)).gt.xcell2)
             pos2(:,1) = pos2(:,1) -sign(xcell,pos2(:,1))
          end where
       enddo
       do while (any(abs(pos2(1:n,2)).gt.ycell2))
          where (    abs(pos2(1:n,2)).gt.ycell2)
             pos2(:,2) = pos2(:,2) -sign(ycell,pos2(:,2))
          end where
       enddo
       do while (any(abs(pos2(1:n,3)).gt.zcell2))
          where (    abs(pos2(1:n,3)).gt.zcell2)
             pos2(:,3) = pos2(:,3) -sign(zcell,pos2(:,3))
          end where
       enddo
      return
      end
