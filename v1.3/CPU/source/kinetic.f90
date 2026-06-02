!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine kinetic  --  compute kinetic energy components  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "kinetic" computes the total kinetic energy and kinetic energy
!     contributions to the pressure tensor by summing over velocities
!
!
!> @brief 
!> computes the total kinetic energy and kinetic energy
!> contributions to the pressure tensor by summing over velocities
!> @param no params
subroutine kinetic (eksum,ekin,temp)
   use atmtyp
   use atoms
   use bath
   use domdec
   use group
   use inform
   use iounit
   use mdstuf
   use moldyn
   use units
   use usage
   use mpi
   use qtb, only: corr_fact_qtb
   implicit none
   integer i,j,k,ierr,iglob

   real*8 eksum,temp
   real*8 term,value



   real*8 ekin(3,3)
!
   if (deb_Path) write(iout,*), 'kinetic '
!

!
!
!     zero out the total kinetic energy and its outer product
!
   eksum = 0.0d0
   ekin = 0.0d0
!
!     get the total kinetic energy and tensor for atomic sites
!
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

   do j=1,3; do k=1,3
         ekin(k,j)=ekin(k,j)/sqrt(corr_fact_qtb(k)*corr_fact_qtb(j))
      enddo; enddo

   call MPI_ALLREDUCE(MPI_IN_PLACE,ekin,9,MPI_REAL8,MPI_SUM,&
   &COMM_TINKER,ierr)
   eksum = ekin(1,1) + ekin(2,2) + ekin(3,3)
!
   if (isobaric .and. barostat.eq.'BUSSI') then
      term = dble(nfree) * gasconst * kelvin * taupres * taupres
      value = 0.5d0 * term * eta * eta
      do j = 1, 3
         ekin(j,j) = ekin(j,j) + value/3.0d0
      end do
      eksum = eksum + value
   end if
!
!     set the instantaneous temperature from total kinetic energy
!
   temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
   return
end
