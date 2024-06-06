c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module deriv  --  Cartesian coordinate derivative components  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     desum   total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     deba    stretch-bend Cartesian coordinate derivatives
c     deub    Urey-Bradley Cartesian coordinate derivatives
c     deaa    angle-angle Cartesian coordinate derivatives
c     deopb   out-of-plane bend Cartesian coordinate derivatives
c     deopd   out-of-plane distance Cartesian coordinate derivatives
c     deid    improper dihedral Cartesian coordinate derivatives
c     deit    improper torsion Cartesian coordinate derivatives
c     det     torsional Cartesian coordinate derivatives
c     dept    pi-orbital torsion Cartesian coordinate derivatives
c     debt    stretch-torsion Cartesian coordinate derivatives
c     deat    angle-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     der     repulsion Cartesian coordinate derivatives
c     dedsp   dispersion Cartesian coordinate derivatives
c     dedsprec   reciprocal dispersion Cartesian coordinate derivatives
c     dect    charge transfer Cartesian coordinate derivatives
c     dec     charge-charge Cartesian coordinate derivatives
c     decrec  reciprocal charge-charge Cartesian coordinate derivatives
c     dem     multipole Cartesian coordinate derivatives
c     demrec  reciprocal multipole Cartesian coordinate derivatives
c     dep     polarization Cartesian coordinate derivatives
c     dep     reciprocal polarization Cartesian coordinate derivatives
c     deg     geometric restraint Cartesian coordinate derivatives
c     dex     extra energy term Cartesian coordinate derivatives
c     desave  stored Cartesian coordinate derivatives
c     desmd   extra smd energy term Cartesian coordinate derivatives
c
c     Lambda-dynamics derivatives
c
c     delambda           hamiltonian derivative with respect to lambda (to be sent to colvar)
c     delambdae          hamiltonian derivative with respect to elambda
c     delambdav          hamiltonian derivative with respect to vlambda
c     delambdaesave      stored hamiltonian derivative with respect to elambda
c     delambdavsave      stored hamiltonian derivative with respect to vlambda
c     dlambdaelambda     derivative of elambda with respect to lambda
c     dlambdavlambda     derivative of vlambda with respect to lambda     
c
c     Orthogonal Space Random Walk - note x stands for Cartesian coordinates
c     dxdelambda         hamiltonian double derivative with respect to x and lambda (to be sent to colvar)
c     dxdelambdae        hamiltonian double derivative with respect to x and elambda (electrostatic interactions)
c     dxdelambdav        hamiltonian double derivative with respect to x and vlambda (vdw interactions)
c     d2edlambda2         hamiltonian double derivative with respect to lambda (to be sent to colvar)
c     d2edlambdae2        hamiltonian double derivative with respect to elambda (electrostatic interactions)
c     d2edlambdav2        hamiltonian double derivative with respect to vlambda (vdw interactions)
c
c     dotstgrad : flag when the main program is testgrad (communication
c      of the forces one by one)
c
c
      module deriv
      implicit none
      real*8, allocatable :: desum(:,:),deb(:,:),dea(:,:),deba(:,:)
      real*8, allocatable :: deub(:,:),deaa(:,:),deopb(:,:),deopd(:,:)
      real*8, allocatable :: deid(:,:),det(:,:),dept(:,:),deit(:,:)
      real*8, allocatable :: deat(:,:)
      real*8, allocatable :: debt(:,:),dett(:,:),dev(:,:),dec(:,:)
      real*8, allocatable :: der(:,:),dedsp(:,:),dect(:,:)
      real*8, allocatable :: dedsprec(:,:)
      real*8, allocatable :: dem(:,:),dep(:,:)
      real*8, allocatable :: deg(:,:),dex(:,:)
      real*8, allocatable :: decrec(:,:),demrec(:,:),deprec(:,:)
      real*8, allocatable :: debond(:,:),desave(:,:)
      real*8, allocatable :: desmd(:,:)
      real*8 :: delambda,delambdae,delambdav
      real*8 :: delambdaesave,delambdavsave
      real*8 :: d2edlambda2,d2edlambdae2,d2edlambdav2
      real*8 :: dlambdaelambda, dlambdavlambda
      real*8, allocatable  :: dxdelambda(:,:)
      real*8, allocatable :: dxdelambdae(:,:), dxdelambdav(:,:)
      logical dotstgrad
      logical abortall
      integer inte(2)
      integer cBond,cNBond,cSNBond,cDef
      enum,bind(C)
      enumerator idBond,idSNBond,idNBond
      end enum
      parameter( cBond=1,cSNBond=2,cNBond=4,cDef=5 )
      save

      contains

      subroutine resetForcesRec
      implicit none
        if(allocated(demrec)) demrec = 0.0d0
        if(allocated(decrec)) decrec = 0.0d0
        if(allocated(deprec)) deprec = 0.0d0
        if(allocated(dedsprec)) dedsprec = 0.0d0
      end subroutine resetForcesRec


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
                  print*,j,iglob,rank,x(iglob),y(iglob)
     &                  ,z(iglob),val
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
                  call MPI_AllGather(MPI_IN_PLACE,lgli,MPI_DATATYPE_NULL
     $                              ,gli,lgli,MPI_INT,COMM_TINKER,i1)
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
                  else if (cap.gt.0.and.abs(gli(i2,i1)-inte(1)).lt.5)
     &                 then
                     cap = cap + 1
                     inte(2) = gli(i2,i1)
                  else
                     cap = cap + 1
                  end if
               end if
            end do; end do
            
            if (cap.ne.2.and.rank.eq.0) then
               print*,' more than one interactions found '
     &               ,rank,cap
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
