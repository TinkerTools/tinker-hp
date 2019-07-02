c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kinetic  --  compute kinetic energy components  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kinetic" computes the total kinetic energy and kinetic energy
c     contributions to the pressure tensor by summing over velocities
c
c
      subroutine kinetic (eksum,ekin,temp)
      use atmtyp
      use atoms
      use bath
      use domdec
      use group
      use mdstuf
      use moldyn
      use units
      use usage
      use mpi
      implicit none
      integer i,j,k,ierr,iglob
      integer start,stop
      real*8 eksum,weigh,temp
      real*8 term,value
      real*8 xr,yr,zr
      real*8 x2,y2,z2
      real*8 xcm,ycm,zcm
      real*8 ekin(3,3)
      real*8 inert(3,3)
c
c
c     zero out the total kinetic energy and its outer product
c
      eksum = 0.0d0
      ekin = 0.0d0
c
c     get the total kinetic energy and tensor for atomic sites
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            term = 0.5d0 * mass(iglob) / convert
            do j = 1, 3
               do k = 1, 3
                  value = term * v(j,iglob) * v(k,iglob)
                  ekin(k,j) = ekin(k,j) + value
               end do
            end do
         end if
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,ekin,9,MPI_REAL8,MPI_SUM,
     $  MPI_COMM_WORLD,ierr)
      eksum = ekin(1,1) + ekin(2,2) + ekin(3,3)
c
      if (isobaric .and. barostat.eq.'BUSSI') then
         term = dble(nfree) * gasconst * kelvin * taupres * taupres
         value = 0.5d0 * term * eta * eta
         do j = 1, 3
            ekin(j,j) = ekin(j,j) + value/3.0d0
         end do
         eksum = eksum + value
      end if
c
c     set the instantaneous temperature from total kinetic energy
c
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
      return
      end
