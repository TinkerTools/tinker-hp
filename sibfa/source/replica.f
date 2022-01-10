c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine replica  --  periodicity via cell replication  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "replica" decides between images and replicates for generation
c     of periodic boundary conditions, and sets the cell replicate
c     list if the replicates method is to be used
c
c     only for dispersion for now
c
c
      subroutine replica (cutoff)
      use bound
      use boxes
      use cell
      use inform
      use iounit
      implicit none
      integer i,j,k
      integer nx,ny,nz
      real*8 cutoff,maximage
      real*8 xlimit,ylimit,zlimit
c
c     only necessary if periodic boundaries are in use
c
      ncell = 1
      if (.not. use_bounds)  return
c
c     find the maximum sphere radius inscribed in periodic box
c
c      if (orthogonal) then
         xlimit = xbox2
         ylimit = ybox2
         zlimit = zbox2
c      else if (monoclinic) then
c         xlimit = xbox2 * beta_sin
c         ylimit = ybox2
c         zlimit = zbox2 * beta_sin
c      else if (triclinic) then
c         xlimit = xbox2 * beta_sin * gamma_sin
c         ylimit = ybox2 * gamma_sin
c         zlimit = zbox2 * beta_sin
c      else if (octahedron) then
c         xlimit = (sqrt(3.0d0)/4.0d0) * xbox
c         ylimit = xlimit
c         zlimit = xlimit
c      end if
      maximage = min(xlimit,ylimit,zlimit)
cc
cc     use replicate method to handle cutoffs too large for images
cc
c      if (cutoff .le. maximage) then
c         use_replica = .false.
c      else
c         use_replica = .true.
c      end if
cc
cc     truncated octahedron cannot use the replicates method
cc
c      if (octahedron .and. use_replica) then
c         write (iout,10)
c   10    format (/,' REPLICA  --  Truncated Octahedron',
c     &              ' cannot be Replicated')
c         call fatal
c      end if
c
c     find the number of replicates needed based on cutoff
c
      nx = int(cutoff/xlimit)
      ny = int(cutoff/ylimit)
      nz = int(cutoff/zlimit)
      if (cutoff .gt. dble(nx)*xlimit)  nx = nx + 1
      if (cutoff .gt. dble(ny)*ylimit)  ny = ny + 1
      if (cutoff .gt. dble(nz)*zlimit)  nz = nz + 1
      if (nx .lt. 1)  nx = 1
      if (ny .lt. 1)  ny = 1
      if (nz .lt. 1)  nz = 1
cc
cc     set the replicated cell length and the half width
cc
c      xcell = dble(nx) * xbox
c      ycell = dble(ny) * ybox
c      zcell = dble(nz) * zbox
c      xcell2 = 0.5d0 * xcell
c      ycell2 = 0.5d0 * ycell
c      zcell2 = 0.5d0 * zcell
c
c     perform dynamic allocation of some global arrays
c
      ncell = nx*ny*nz
      if (allocated(icell)) then
         if (size(icell) .lt. 3*ncell) then
            deallocate (icell)
            allocate (icell(3,ncell))
         end if
      else
         allocate (icell(3,ncell))
      end if
c
c     assign indices to the required cell replicates
c
      ncell = 0
      do k = 0, nz-1
         do j = 0, ny-1
            do i = 0, nx-1
               ncell = ncell + 1
               icell(1,ncell) = i
               icell(2,ncell) = j
               icell(3,ncell) = k
            end do
         end do
      end do
c
c     print a message indicating the number of replicates used
c
      if (debug .and. ncell.gt.1) then
         if (max(nx,ny,nz) .lt. 10) then
            write (iout,20)  nx,ny,nz
   20       format (/,' REPLICA  --  Period Boundaries via',i2,' x',
     &                 i2,' x',i2,' Cell Replicate Set')
         else if (max(nx,ny,nz) .lt. 100) then
            write (iout,30)  nx,ny,nz
   30       format (/,' REPLICA  --  Period Boundaries via',i3,' x',
     &                 i3,' x',i3,' Cell Replicate Set')
         else
            write (iout,40)  nx,ny,nz
   40       format (/,' REPLICA  --  Period Boundaries via',i4,' x',
     &                 i4,' x',i4,' Cell Replicate Set')
         end if
      end if
      return
      end
