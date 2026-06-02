!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module deriv  --  Cartesian coordinate derivative components  ##
!     ##                                                                ##
!     ####################################################################
!
!
module deriv
   implicit none
   real*8, allocatable :: desum(:,:) !<total energy Cartesian coordinate derivatives
   real*8, allocatable :: deb(:,:) !<bond stretch Cartesian coordinate derivatives
   real*8, allocatable :: dea(:,:) !<angle bend Cartesian coordinate derivatives
   real*8, allocatable :: deba(:,:) !<stretch-bend Cartesian coordinate derivatives
   real*8, allocatable :: deub(:,:) !<Urey-Bradley Cartesian coordinate derivatives
   real*8, allocatable :: deaa(:,:) !<angle-angle Cartesian coordinate derivatives
   real*8, allocatable :: deopb(:,:) !<out-of-plane bend Cartesian coordinate derivatives
   real*8, allocatable :: deopd(:,:) !<out-of-plane distance Cartesian coordinate derivatives
   real*8, allocatable :: deid(:,:) !<improper dihedral Cartesian coordinate derivatives
   real*8, allocatable :: det(:,:) !<torsional Cartesian coordinate derivatives
   real*8, allocatable :: dept(:,:) !<pi-orbital torsion Cartesian coordinate derivatives
   real*8, allocatable :: deit(:,:) !<improper torsion Cartesian coordinate derivatives
   real*8, allocatable :: deat(:,:) !<angle-torsion Cartesian coordinate derivatives
   real*8, allocatable :: debt(:,:) !<stretch-torsion Cartesian coordinate derivatives
   real*8, allocatable :: dett(:,:) !<torsion-torsion Cartesian coordinate derivatives
   real*8, allocatable :: dev(:,:) !<van der Waals Cartesian coordinate derivatives
   real*8, allocatable :: dec(:,:) !<charge-charge Cartesian coordinate derivatives
   real*8, allocatable :: der(:,:) !<repulsion Cartesian coordinate derivatives
   real*8, allocatable :: dedsp(:,:) !<dispersion Cartesian coordinate derivatives
   real*8, allocatable :: dect(:,:) !<charge transfer Cartesian coordinate derivatives
   real*8, allocatable :: dedsprec(:,:) !<reciprocal dispersion Cartesian coordinate derivatives
   real*8, allocatable :: dem(:,:) !<multipole Cartesian coordinate derivatives
   real*8, allocatable :: dep(:,:) !<polarization Cartesian coordinate derivatives
   real*8, allocatable :: deg(:,:) !<geometric restraint Cartesian coordinate derivatives
   real*8, allocatable :: dex(:,:) !<extra energy term Cartesian coordinate derivatives
   real*8, allocatable :: decrec(:,:) !<reciprocal charge-charge Cartesian coordinate derivatives
   real*8, allocatable :: demrec(:,:) !<reciprocal multipole Cartesian coordinate derivatives
   real*8, allocatable :: deprec(:,:) !<reciprocal polarization Cartesian coordinate derivatives
   real*8, allocatable :: debond(:,:) !bonded energy Cartesian coordinate derivatives
   real*8, allocatable :: desave(:,:) !<stored Cartesian coordinate derivatives
   real*8, allocatable :: desmd(:,:) !<extra smd energy term Cartesian coordinate derivatives
   real*8 :: delambda !<hamiltonian derivative with respect to lambda (to be sent to colvar)
   real*8 :: delambdae !<hamiltonian derivative with respect to elambda
   real*8 :: delambdav !<hamiltonian derivative with respect to vlambda
   real*8 :: delambdaesave !<stored hamiltonian derivative with respect to elambda
   real*8 :: delambdavsave !<stored hamiltonian derivative with respect to vlambda
   real*8 :: d2edlambda2 !<hamiltonian double derivative with respect to lambda (to be sent to colvar)
   real*8 :: d2edlambdae2 !<hamiltonian double derivative with respect to elambda (electrostatic interactions)
   real*8 :: d2edlambdav2 !<hamiltonian double derivative with respect to vlambda (vdw interactions)

   real*8 :: dlambdaelambda !<derivative of elambda with respect to lambda
   real*8 :: dlambdavlambda !<derivative of vlambda with respect to lambda
   real*8, allocatable  :: dxdelambda(:,:) !<hamiltonian double derivative with respect to x and lambda (to be sent to colvar)
   real*8, allocatable :: dxdelambdae(:,:) !<hamiltonian double derivative with respect to x and elambda (electrostatic interactions)
   real*8, allocatable :: dxdelambdav(:,:) !<hamiltonian double derivative with respect to x and vlambda (vdw interactions)
   logical :: dotstgrad !<flag when the main program is testgrad (communication of the forces one by one)

   logical :: abortall !< flag when all the processes receive abort instruction
   integer :: inte(2)
   integer :: cBond
   integer :: cNBond
   integer :: cSNBond
   integer :: cDef
   enum,bind(C)
      enumerator idBond,idSNBond,idNBond
   end enum
   parameter( cBond=1,cSNBond=2,cNBond=4,cDef=5 )
   save

contains

!> @brief 
!> reset all the arrays of reciprocal forces
!> @param no params
   subroutine resetForcesRec
      implicit none
      if(allocated(demrec)) demrec = 0.0d0
      if(allocated(decrec)) decrec = 0.0d0
      if(allocated(deprec)) deprec = 0.0d0
      if(allocated(dedsprec)) dedsprec = 0.0d0
   end subroutine resetForcesRec


!> @brief 
!> get the min, the max and the L1 norm of a vector
!> @param no params
   subroutine minmaxone1( mi,ma,on,vector,sz,name )
      use atoms
      use domdec
      use inform
      use mpi
      implicit none
      integer sz
      integer,parameter::lgli=10
      real(8) mi,ma,on
      real*8 vector(*)
      character(*),optional,intent(in)::name
      integer i,j,i1,i2,iglob,cap,cap1,cap2
      integer gli(lgli,nproc)
      real(8) val

      abortall = .false.
      if (present(name)) then
         cap = 1
         gli = 0
         if (name.eq.'devi') then
            do i = 1, sz/3
               cap2  = 0
               iglob = glob(i)
               do j  = 1,3
                  val   = vector((i-1)*3+j)
                  if (abs(val).gt.90.0) then
                     print*,j,iglob,rank,x(iglob),y(iglob)&
                     &,z(iglob),val
                     abort=.true.
                     cap2 = cap2 + 1
                  end if
               end do
               if (cap2.gt.0) then
                  cap1 = cap
                  cap  = cap + 1
                  if (cap1.le.lgli) gli(cap1,rank+1) = iglob
               end if
            end do
            do i = 1, nproc
               if (rank.eq.i-1) abortall=abort
               call MPI_BCAST(abortall,1,MPI_LOGICAL,i-1,COMM_TINKER,i1)
               if (abortall) then
                  abort = .true.
                  call MPI_AllGather(MPI_IN_PLACE,lgli,MPI_DATATYPE_NULL&
                  &,gli,lgli,MPI_INT,COMM_TINKER,i1)
                  exit
               end if
            end do

            if (abort) then

               cap  = 0
               cap1 = 0
               do i1 = 1, nproc; do i2 = 1,lgli
                     if ( gli(i2,i1).ne.0 ) then
                        if (cap.eq.0) then
                           cap = 1
                           inte(cap) = gli(i2,i1)
                        else if (cap.gt.0.and.abs(gli(i2,i1)-inte(1)).lt.5)&
                        &then
                           cap = cap + 1
                           inte(2) = gli(i2,i1)
                        else
                           cap = cap + 1
                        end if
                     end if
                  end do; end do

               if (cap.ne.2.and.rank.eq.0) then
                  print*,' more than one interactions found '&
                  &,rank,cap
                  do i = 1,nproc
                     do j = 1,lgli
                        if (gli(j,i).ne.0) write(*,'(I10,$)') gli(j,i)
                     end do
                     print*
                  end do
               end if

            end if
         end if

      end if

      do i = 1, sz
         val = vector(i)
         mi = min( mi,val )
         ma = max( ma,val )
         on = on + abs(val)
      end do
   end subroutine

end
