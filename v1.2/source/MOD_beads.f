c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module beads   --  pimd variables                              ##
c     ##                                                                 ##
c     #####################################################################
c
c
      module beads
      use dcdmod
      implicit none
      logical :: path_integral_md=.FALSE.
      integer :: nbeads,nbeadsloc   

      character(len=80) :: save_beads

      logical :: contract
      logical :: centroid_longrange
      logical :: centroid_recip
      logical :: polar_allbeads
      logical :: cay_correction
      integer :: nbeads_ctr, nbeadsloc_ctr  

      logical :: pi_comm_report
      integer :: pi_comm_unit
      logical :: PITIMER
      
      integer :: max_polymer_index=0
      integer :: ipolymer_loaded=-1
      integer :: ibead_loaded=-1
      logical :: add_buffer_pi


      integer, allocatable :: ilocpi_beg(:)
      integer, allocatable :: ilocpi_end(:)
      integer, allocatable :: nlocpis(:)
      integer :: nlocpi,nlocpiproc

      real*8  :: pi_start_isobaric
      real*8  :: lambda_trpmd
      logical :: default_lambda_trpmd=.TRUE.
      real*8 :: gyrpi
      real*8 :: temppi,temppi_cl,temppi_centroid
      real*8 :: epotpi_loc,etotpi_loc
      real*8 :: eksumpi_loc,ekinpi_loc(3,3)
      real*8 :: ekdynpi, epotpi, eintrapi, einterpi
      real*8 :: ekprim,ekvir,presvir
      real*8 :: ekvir_fast
      real*8 :: ekcentroid
      real*8 :: eintrapi_loc
      real*8 :: einterpi_loc
      real*8 :: virpi(3,3),ekinpi(3,3),ekinpi_fast(3,3)
      real*8 :: stresspi(3,3)
      real*8 :: dedvpi,dedvintrapi,dedvinterpi
      real*8 :: dippi(3), dipindpi(3)
      real*8 :: delambdavpi=0.d0, delambdaepi=0.d0 ! lambda dynamics
      
      logical :: pot_gathered=.TRUE.

      real(8), allocatable :: eigmat(:,:)
      real(8), allocatable ::  omkpi(:)
      
      real(8), allocatable :: contractor_mat(:,:)

      TYPE POLYMER_COMM_TYPE
        integer :: nbeads
        integer :: nbeadsproc
        integer,allocatable :: nbeadsloc(:)
        integer,allocatable :: ibead_beg(:)
        integer,allocatable :: ibead_end(:)
        integer,allocatable :: iproc(:)
        integer :: index

        logical :: allocated=.FALSE.

        real*8, allocatable :: pos(:,:,:),vel(:,:,:)
        real*8, allocatable :: eigpos(:,:,:),eigvel(:,:,:)
        real*8, allocatable :: forces(:,:,:)
        real*8, allocatable :: eigforces(:,:,:)
        real*8, allocatable :: forces_slow(:,:,:)
        
        real*8, allocatable :: epot(:)
        real*8, allocatable :: vir(:,:,:)

        real*8, allocatable :: dip(:,:)
        real*8, allocatable :: dipind(:,:)

        real*8 :: gyr

        type(dcdinfo_t),allocatable :: dcdinfo(:)
      END TYPE     

      save

      contains

      subroutine allocate_polymer(polymer,nu,allocate_vel
     &      ,allocate_forces,allocate_slow)
        use atoms, only: n
        use domdec, only: nproc_polymer
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
        integer, intent(in) :: nu
        LOGICAL, intent(in) :: allocate_vel, allocate_slow
        LOGICAL, intent(in) :: allocate_forces
        integer :: i

        call deallocate_polymer(polymer)

        max_polymer_index = max_polymer_index +1
        polymer%index = max_polymer_index

        polymer%nbeads=nu
        polymer%nbeadsproc = int(nu/nproc_polymer)
        allocate(polymer%ibead_beg(nproc_polymer))
        allocate(polymer%ibead_end(nproc_polymer))
        allocate(polymer%nbeadsloc(nproc_polymer))
        
        do i=1,nproc_polymer
          polymer%ibead_beg(i)=1+(i-1)*polymer%nbeadsproc
          polymer%ibead_end(i)=i*polymer%nbeadsproc
        enddo
        polymer%ibead_end(nproc_polymer)=nu
        do i=1,nproc_polymer
          polymer%nbeadsloc(i)=polymer%ibead_end(i)
     &                        -polymer%ibead_beg(i)+1
        enddo
        if(nproc_polymer>1) then
          allocate(polymer%iproc(nu))
          do i=1,nu
            polymer%iproc(i)=
     &         min(nproc_polymer,int((i-1)/polymer%nbeadsproc)+1)
          enddo
        endif

        allocate(polymer%epot(nu))
        allocate(polymer%vir(3,3,nu))
        allocate(polymer%dip(3,nu))
        allocate(polymer%dipind(3,nu))
        polymer%epot(:)=0.0d0
        polymer%vir(:,:,:)=0.0d0
        polymer%dip(:,:)=0.0d0
        polymer%dipind(:,:)=0.0d0

        allocate(polymer%pos(3,n,nu))
        if(allocate_vel) allocate(polymer%vel(3,n,nu))
        if(allocate_forces) allocate(polymer%forces(3,n,nu))
        if(allocate_slow) allocate(polymer%forces_slow(3,n,nu))

        polymer%allocated=.true.

      end subroutine allocate_polymer

      subroutine deallocate_polymer(polymer)
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer

        if(.not. polymer%allocated) return
        polymer%allocated=.false.

        if(allocated(polymer%pos)) then
          deallocate(polymer%pos)
        endif
        if(allocated(polymer%eigpos)) then
          deallocate(polymer%eigpos)
        endif
        if(allocated(polymer%vel)) then
          deallocate(polymer%vel)
        endif
        if(allocated(polymer%eigvel)) then
          deallocate(polymer%eigvel)
        endif
        if(allocated(polymer%forces)) then
          deallocate(polymer%forces)
        endif  
        if(allocated(polymer%eigforces)) then
          deallocate(polymer%eigforces)
        endif  

        if(allocated(polymer%forces_slow)) then
          deallocate(polymer%forces_slow)
        endif        

        if(allocated(polymer%epot)) then
          deallocate(polymer%epot)
        endif      
        
        if(allocated(polymer%vir)) then
          deallocate(polymer%vir)
        endif  

        if(allocated(polymer%dip)) then
          deallocate(polymer%dip)
        endif

        if(allocated(polymer%dipind)) then
          deallocate(polymer%dipind)
        endif 

        if(allocated(polymer%nbeadsloc)) then
          deallocate(polymer%nbeadsloc)
        endif 

        if(allocated(polymer%ibead_beg)) then
          deallocate(polymer%ibead_beg)
        endif 

        if(allocated(polymer%ibead_end)) then
          deallocate(polymer%ibead_end)
        endif 

        if(allocated(polymer%iproc)) then
          deallocate(polymer%iproc)
        endif 
        
      end subroutine deallocate_polymer

      subroutine initialize_normal_modes(polymer,verbose)
      use domdec
      use atoms
      use ascii, only: int_to_str
      use mpi
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: verbose
      integer :: k,ierr,nu,i,j
  
      nu=polymer%nbeads

      if(allocated(polymer%eigpos)) then
        deallocate(polymer%eigvel,polymer%eigpos)
      endif
      allocate(polymer%eigpos(3,n,nu))
      allocate(polymer%eigvel(3,n,nu))

      if(allocated(polymer%eigforces)) then
        deallocate(polymer%eigforces)
      endif
      allocate(polymer%eigforces(3,n,nu))

      do k=1,nu; do i=1,n; do j=1,3
        polymer%eigpos(j,i,k)=0.d0
        polymer%eigvel(j,i,k)=0.d0
        polymer%eigforces(j,i,k)=0.d0
      enddo ;enddo; enddo

      if(ranktot==0) write(0,*) "normal modes initialized"
      end subroutine initialize_normal_modes

      subroutine update_nlocpi(nloc)
        use domdec, only: nproc_polymer,rank_polymer
        implicit none
        integer, intent(in) :: nloc
        integer :: i

        if(.not.allocated(ilocpi_beg)) then
          allocate(ilocpi_beg(nproc_polymer))
          allocate(ilocpi_end(nproc_polymer))
          allocate(nlocpis(nproc_polymer))
        endif

        nlocpiproc = int(nloc/nproc_polymer)
        do i=1,nproc_polymer
          ilocpi_beg(i)=1+(i-1)*nlocpiproc
          ilocpi_end(i)=i*nlocpiproc
        enddo
        ilocpi_end(nproc_polymer)=nloc
        do i=1,nproc_polymer
          nlocpis(i)=ilocpi_end(i) - ilocpi_beg(i) + 1
        enddo
        nlocpi = nlocpis(rank_polymer+1)

      end subroutine update_nlocpi

      subroutine allocpi(polymer,polymer_ctr)
      use angle
      use atoms
      use bath
      use bitor
      use domdec
      use molcul
      use neigh
      use sizes
      use pitors
      use potent
      use tors
      use uprior
      use units
      use mdstuf
      use mpi
      implicit none
      integer ibead,i,ierr,j,k,contraction_level
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer,polymer_ctr
      real(8), allocatable :: WORK(:)
      real(8), allocatable :: WORK_ctr(:),eigMat_ctr(:,:),omkpi_ctr(:) 

      call update_nlocpi(nloc)
      
      if (allocated(eigmat)) then
          deallocate(eigmat)
      end if
      allocate(eigmat(nbeads,nbeads))

      if (allocated(omkpi)) then
          deallocate(omkpi)
      end if
      allocate(omkpi(nbeads))        

      if(ranktot==0) then
        allocate(WORK(3*nbeads))  
        eigmat=0
        DO i=1,nbeads-1
          eigmat(i,i)=2
          eigmat(i+1,i)=-1
          eigmat(i,i+1)=-1
        ENDDO
        eigmat(1,nbeads)=-1
        eigmat(nbeads,1)=-1
        eigmat(nbeads,nbeads)=2
        call DSYEV('V','U',nbeads,eigMat,nbeads, 
     $       omkpi,WORK,3*nbeads,ierr)
        omkpi(1)=0
        omkpi(:)=sqrt(omkpi)*(nbeads*boltzmann*kelvin/hbar_planck)
        do i=1,nbeads
          if(eigmat(i,1)<0) eigmat(i,:)=-eigmat(i,:)
        enddo
        deallocate(WORK)
      endif
      call MPI_BCAST(eigmat,nbeads**2,MPI_REAL8
     &    ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(omkpi,nbeads,MPI_REAL8
     &    ,0,MPI_COMM_WORLD,ierr)


      if (contract .and. nbeads_ctr>1) then
        if (allocated(contractor_mat)) then
          deallocate(contractor_mat)
        endif
        allocate(contractor_mat(nbeads_ctr,nbeads))
        if(ranktot==0) then
          allocate(eigmat_ctr(nbeads_ctr,nbeads_ctr))
          allocate(omkpi_ctr(nbeads_ctr))
          allocate(WORK_ctr(3*nbeads_ctr))
          eigmat_ctr=0
          do i=1,nbeads_ctr-1
            eigmat_ctr(i,i)=2
            eigmat_ctr(i+1,i)=-1
            eigmat_ctr(i,i+1)=-1
          enddo
          eigmat_ctr(1,nbeads_ctr)=-1
          eigmat_ctr(nbeads_ctr,1)=-1
          eigmat_ctr(nbeads_ctr,nbeads_ctr)=2
          call DSYEV('V','U',nbeads_ctr,eigMat_ctr,nbeads_ctr, 
     $        omkpi_ctr,WORK_ctr,3*nbeads_ctr,ierr)
          do i=1,nbeads_ctr
            if(eigMat_ctr(i,1)<0) eigMat_ctr(i,:)=-eigMat_ctr(i,:)
          enddo
          
          contractor_mat(:,:) = 0._8
          do i=1,nbeads_ctr ; do j=1,nbeads            
            do k=1,nbeads_ctr
              contractor_mat(i,j) = contractor_mat(i,j)
     &          + eigmat_ctr(i,k)*eigmat(j,k)
            enddo
          enddo; enddo
          contractor_mat=contractor_mat*sqrt(nbeads_ctr*1._8/nbeads)
          deallocate(WORK_ctr,eigMat_ctr,omkpi_ctr)
        endif
        call MPI_BCAST(contractor_mat,nbeads*nbeads_ctr,MPI_REAL8
     &    ,0,MPI_COMM_WORLD,ierr)
      endif
      
      call allocate_polymer(polymer,nbeads
     &                ,.TRUE., .TRUE., .TRUE. )
      nbeadsloc = polymer%nbeadsloc(rank_polymer+1)
      call initialize_normal_modes(polymer,.TRUE.)

      if(contract .and. nbeads_ctr>1) then        
        call allocate_polymer(polymer_ctr,nbeads_ctr
     &                    ,.FALSE.,.FALSE.,.TRUE.)
        nbeadsloc_ctr = polymer_ctr%nbeadsloc(rank_polymer+1)
      endif

      end subroutine allocpi

      subroutine update_normal_modes_pi(polymer)
      use atoms, only: n
      use domdec
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu
      real*8 :: sqrtnuinv

        nu=polymer%nbeads
        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)      
        sqrtnuinv = 1./sqrt(real(nu,8))
        DO k=1,nu; DO iloc=ilocbeg,ilocend 
          i=glob(iloc)
          polymer%eigpos(:,i,k) = 0.
          polymer%eigvel(:,i,k) = 0.
          do ibead=1,nu
            polymer%eigpos(:,i,k) = polymer%eigpos(:,i,k) 
     &          + eigmat(ibead,k)*polymer%pos(:,i,ibead)
            polymer%eigvel(:,i,k) = polymer%eigvel(:,i,k) 
     &          + eigmat(ibead,k)*polymer%vel(:,i,ibead)
          enddo
          polymer%eigpos(:,i,k) = polymer%eigpos(:,i,k)
     &                               *sqrtnuinv
          polymer%eigvel(:,i,k) = polymer%eigvel(:,i,k)
     &                               *sqrtnuinv
        ENDDO ; ENDDO 


      end subroutine update_normal_modes_pi

      subroutine set_eigforces_pi(polymer,forces)
      use atoms, only: n
      use domdec
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real*8, intent(in) :: forces(:,:,:)
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu
      real*8 :: sqrtnuinv

        nu=polymer%nbeads
        sqrtnuinv = 1./sqrt(real(nu,8))
        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)      
        DO k=1,nu; DO iloc=ilocbeg,ilocend 
          i=glob(iloc)
          polymer%eigforces(:,i,k) = 0.d0
          do ibead=1,nu
            polymer%eigforces(:,i,k) = polymer%eigforces(:,i,k) 
     &          + eigmat(ibead,k)*forces(:,i,ibead)*sqrtnuinv
          enddo                           
        ENDDO ; ENDDO 

      end subroutine set_eigforces_pi

      subroutine update_direct_space_pi(polymer)
      use atoms, only: n
      use domdec
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu
      real*8 :: sqrtnu

        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)
        nu=polymer%nbeads
        sqrtnu = sqrt(real(nu,8))
        DO ibead=1,nu; DO iloc=ilocbeg,ilocend 
          i=glob(iloc)
          polymer%pos(:,i,ibead) = 0.d0
          polymer%vel(:,i,ibead) = 0.d0
          do k=1,nu
            polymer%pos(:,i,ibead) = polymer%pos(:,i,ibead) 
     &          + eigmat(ibead,k)*polymer%eigpos(:,i,k)
            polymer%vel(:,i,ibead) = polymer%vel(:,i,ibead) 
     &          + eigmat(ibead,k)*polymer%eigvel(:,i,k)
          enddo
          polymer%pos(:,i,ibead) = polymer%pos(:,i,ibead)
     &                               *sqrtnu
          polymer%vel(:,i,ibead) = polymer%vel(:,i,ibead)
     &                               *sqrtnu
        ENDDO ; ENDDO

      end subroutine update_direct_space_pi

      subroutine contract_polymer(polymer,polymer_ctr)
        use domdec
        implicit none
        type(POLYMER_COMM_TYPE), intent(in) :: polymer
        type(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
        integer :: i,j,k,ibead,iloc,ilocbeg,ilocend
        integer :: nu,nuctr
        
        nuctr=polymer_ctr%nbeads
        nu=polymer%nbeads
        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)
        DO ibead=1,nuctr; DO iloc=ilocbeg,ilocend 
          i=glob(iloc)
          polymer_ctr%pos(:,i,ibead) = 0.d0 
          do k=1,nu
            polymer_ctr%pos(:,i,ibead) = polymer_ctr%pos(:,i,ibead) 
     &          + contractor_mat(ibead,k)*polymer%pos(:,i,k)
          enddo
        ENDDO ;ENDDO
      end subroutine contract_polymer

      subroutine project_forces_ctr(polymer,polymer_ctr,slow)
      use atoms, only: n
      use domdec
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer,polymer_ctr
      logical, intent(in), optional :: slow
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu,nuctr
      real*8 :: nuratio
      logical :: slow_

        slow_ = .false.
        if(present(slow)) slow_=slow

        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)
        nu=polymer%nbeads
        nuctr=polymer_ctr%nbeads
        nuratio=real(nu,8)/real(nuctr,8)
        if(slow_) then
          DO ibead=1,nu; do k=1,nuctr; DO iloc=ilocbeg,ilocend 
            i=glob(iloc)
            polymer%forces_slow(:,i,ibead)=
     &        polymer%forces_slow(:,i,ibead)
     &           + contractor_mat(k,ibead)*nuratio
     &                *polymer_ctr%forces_slow(:,i,k)
          ENDDO ;ENDDO; ENDDO
        else
          DO ibead=1,nu; do k=1,nuctr; DO iloc=ilocbeg,ilocend 
            i=glob(iloc)
            polymer%forces(:,i,ibead)=polymer%forces(:,i,ibead)
     &           + contractor_mat(k,ibead)*nuratio
     &                *polymer_ctr%forces_slow(:,i,k)
          ENDDO ;ENDDO; ENDDO
        endif

      end subroutine project_forces_ctr

      subroutine load_bead(polymer,ibead,force_load)
        use atoms
        use domdec
        implicit none
        type(POLYMER_COMM_TYPE), intent(in) :: polymer
        integer, intent(in) :: ibead
        logical, intent(in), optional :: force_load
        integer :: i,iloc
        logical :: already_loaded,force_load_

        force_load_=.false.
        if(present(force_load)) force_load_=force_load

        already_loaded = (ibead_loaded == ibead 
     &     .and. ipolymer_loaded == polymer%index)
        if(already_loaded .and. (.not. force_load_)) return

        ibead_loaded = ibead
        ipolymer_loaded = polymer%index
        if(ibead == 0) then
          DO iloc=1,nbloc
            i=glob(iloc)
            x(i)=polymer%eigpos(1,i,1)
            y(i)=polymer%eigpos(2,i,1)
            z(i)=polymer%eigpos(3,i,1)
          ENDDO
        else
          DO iloc=1,nbloc
            i=glob(iloc)
            x(i)=polymer%pos(1,i,ibead)
            y(i)=polymer%pos(2,i,ibead)
            z(i)=polymer%pos(3,i,ibead)
          ENDDO
        endif

      end subroutine load_bead

      subroutine stopwatchpi(timer,barrier,reset)
      use mpi
      implicit none
      real(8), intent(inout) :: timer
      logical, intent(in) :: barrier,reset
      real(8), save :: time0, time1
      integer :: ierr

      if(barrier) call MPI_Barrier(MPI_COMM_WORLD,ierr)
      if(reset) then
        time0 = mpi_wtime()
        timer = time0
      else
        time1 = mpi_wtime()
        timer = timer + time1 - time0
        time0 = time1
      endif

      end subroutine stopwatchpi

      subroutine kinetic_pi_fast(polymer)
      use units
      use atoms
      use mdstuf
      use domdec
      use atmtyp
      use mpi
      use bath, only: kelvin
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real*8 :: factor1,rnu,wnu2,ek0
      integer :: i,j,k,l,ibead,iloc,ilocbeg,nu

      nu = polymer%nbeads
      if (nu==1) return
      rnu = real(nu,8)
      factor1=0.5d0/convert
      wnu2=(rnu*boltzmann*kelvin/hbar_planck)**2
      ek0 = 3*n*boltzmann*kelvin*factor1

      ekvir_fast=0.d0
      ekinpi_fast=0.d0

      ilocbeg = ilocpi_beg(rank_polymer+1)
      DO k=2,nu; DO iloc=1,nlocpi
        i=glob(iloc+ilocbeg-1)
        DO j=1,3
          ekvir_fast = ekvir_fast + polymer%eigpos(j,i,k)
     &               *polymer%eigforces(j,i,k)
          DO l=1,3
            ekinpi_fast(l,j) = ekinpi_fast(l,j)
     &         + polymer%eigpos(j,i,k)*polymer%eigforces(l,i,k)
          ENDDO
        ENDDO
      ENDDO ; ENDDO

      end subroutine kinetic_pi_fast


      subroutine kinetic_pi(polymer)
      use units
      use atoms
      use mdstuf
      use domdec
      use atmtyp
      use mpi
      use bath, only: kelvin
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real*8 :: factor1,rnu,wnu2,ek0,hbar2beta,lth2
      integer :: i,j,k,l,ibead,iloc,ilocbeg,nu


      nu = polymer%nbeads
      rnu = real(nu,8)
      factor1=0.5d0/convert
      wnu2=(rnu*boltzmann*kelvin/hbar_planck)**2
      ek0 = 3*n*boltzmann*kelvin*factor1
      hbar2beta=hbar_planck**2/(boltzmann*kelvin)

      ekprim=0.d0
      ekvir=0.d0
      ekdynpi=0.d0
      Ekcentroid=0.d0
      ekinpi=0.d0
      gyrpi=0.d0

      ilocbeg = ilocpi_beg(rank_polymer+1)

      DO iloc=1,nlocpi
        i=glob(iloc+ilocbeg-1)
        DO j=1,3
          Ekcentroid=Ekcentroid+mass(i)*polymer%eigvel(j,i,1)**2
        ENDDO
      ENDDO

      DO k=2,nu; DO iloc=1,nlocpi
        i=glob(iloc+ilocbeg-1)
        lth2=hbar2beta/mass(i)
        DO j=1,3
          gyrpi=gyrpi+
     &     (polymer%pos(j,i,k)-polymer%eigpos(j,i,1))**2/lth2
          ekdynpi = ekdynpi + mass(i)*polymer%eigvel(j,i,k)**2
          ekvir = ekvir + polymer%eigpos(j,i,k)
     &               *polymer%eigforces(j,i,k)
          ekprim = ekprim + mass(i)*( omkpi(k)
     &              *polymer%eigpos(j,i,k) )**2
          DO l=1,3
            ekinpi(l,j) = ekinpi(l,j)
     &         + polymer%eigpos(j,i,k)*polymer%eigforces(l,i,k)
          ENDDO
        ENDDO
      ENDDO ; ENDDO

      if(nu==1) then
        gyrpi=0.d0
        Ekcentroid = Ekcentroid*factor1
        ekdynpi = ekcentroid
        ekprim = ekcentroid
        ekvir = Ekcentroid
        ekvir_fast = 0.0d0  
      else
        ekinpi(:,:) = - 0.5*(ekinpi(:,:)+ekinpi_fast(:,:))
        do j=1,3
          ekinpi(j,j) = ekinpi(j,j) + ek0
        enddo
        ! symmetrize ekinpi
        do j=1,3; do l=1,j-1
          ekinpi(j,l) = 0.5*(ekinpi(j,l)+ekinpi(l,j))
          ekinpi(l,j) = ekinpi(j,l)
        enddo; enddo

        gyrpi=sqrt(hbar2beta*gyrpi/real(3*n*nu,8))
        Ekcentroid = Ekcentroid*factor1
        ekdynpi = ekdynpi*factor1 + ekcentroid
        ekprim = ekdynpi - ekprim*factor1 !/rnu
        ekvir = Ekcentroid - 0.5d0*(ekvir+ekvir_fast) !/rnu
!        ekprim = rnu*ek0 - ekprim*factor1
!        ekvir =  ek0 - 0.5d0*(ekvir+ekvir_save) 
        ekvir_fast = 0.0d0 
        ekinpi_fast(:,:) = 0.0d0 

      endif

      end subroutine kinetic_pi

      subroutine reduce_observables_pi(polymer)
      use atoms
      use units
      use mdstuf
      use domdec
      use atmtyp
      use mpi
      use boxes
      use bath, only:anisotrop
      implicit none
      type(POLYMER_COMM_TYPE), intent(in) :: polymer
      integer :: i,j,k,ibead,iloc,ilocbeg,ierr
      real*8 :: factor2,factor,buffer(12)
      integer,save :: ncomm

      factor2=2.d0/(3*n*gasconst)
      factor=prescon/volbox

      if(nproctot==1) then
        temppi = factor2 * ekvir
        temppi_cl = factor2 *ekdynpi 
        temppi_centroid = factor2*ekcentroid
        dedvpi = (virpi(1,1)+virpi(2,2)+virpi(3,3))/(3.d0*volbox)
        presvir = prescon*( -dedvpi + 2.d0*ekvir
     &                              /(3.d0*volbox) )
        stresspi(:,:)=factor*(2*ekinpi(:,:)-virpi(:,:))
        return
      endif

      ncomm=6
      buffer(1)=ekdynpi
      buffer(2)=ekcentroid
      buffer(3)=ekprim
      buffer(4)=ekvir
      buffer(5)=epotpi
      buffer(6)=(virpi(1,1)+virpi(2,2)+virpi(3,3))/(3.d0*volbox)
      if (anisotrop) then
        ncomm=12
        buffer(7)=stresspi(1,1)
        buffer(8)=stresspi(2,2)
        buffer(9)=stresspi(3,3)
        buffer(10)=0.5*(stresspi(1,2)+stresspi(2,1))
        buffer(11)=0.5*(stresspi(1,3)+stresspi(3,1))
        buffer(12)=0.5*(stresspi(2,3)+stresspi(3,2))
      endif

      if(ranktot==0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer
     &       ,ncomm,MPI_REAL8
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      else
        call MPI_REDUCE(buffer,buffer
     &        ,ncomm,MPI_REAL8
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif

      if(ranktot==0) then
        ekdynpi=buffer(1)
        ekcentroid=buffer(2)
        ekprim=buffer(3)
        ekvir=buffer(4)
        epotpi=buffer(5)
        dedvpi=buffer(6)
        temppi = factor2 * ekvir
        temppi_cl = factor2 *ekdynpi
        temppi_centroid = factor2*ekcentroid  
        presvir = prescon*( -dedvpi + 2.d0*ekvir
     &                              /(3.d0*volbox) )
        if(anisotrop) then
          stresspi(1,1)=buffer(7)
          stresspi(2,2)=buffer(8)
          stresspi(3,3)=buffer(9)
          stresspi(1,2)=buffer(10)
          stresspi(1,3)=buffer(11)
          stresspi(2,3)=buffer(12)
          stresspi(2,1)=buffer(10)
          stresspi(3,1)=buffer(11)
          stresspi(3,2)=buffer(12)
        endif
      endif

      end subroutine reduce_observables_pi

      subroutine reset_observables(init)
        implicit none
        logical, intent(in) :: init
        epotpi_loc = 0
        eksumpi_loc = 0
        ekinpi_loc = 0
        eintrapi_loc = 0
        einterpi_loc = 0
        temppi = 0
        temppi_cl = 0
        temppi_centroid = 0
        eintrapi = 0
        einterpi = 0 
        epotpi = 0
        dedvintrapi = 0
        dedvinterpi = 0
        dedvpi = 0
        Ekcentroid = 0
        ekdynpi = 0
        ekvir = 0
        ekvir_fast = 0 
        ekprim = 0
        presvir = 0
        virpi(:,:) = 0
        ekinpi(:,:) = 0
        ekinpi_fast(:,:) = 0
        stresspi(:,:) = 0
        dippi(:) = 0
        dipindpi(:) = 0
        gyrpi = 0

      end subroutine reset_observables

      subroutine read_bead_config()
      use keys
      use domdec
      implicit none
      character*20  :: keyword
      character*240 :: record
      character*240 :: string
      character*240 :: xyzfile
      integer :: next,i,ios

c
c     check for keywords containing any altered parameters
c
      nbeads = 1 
      nbeads_ctr = 0
      pi_start_isobaric=-1.
      centroid_longrange=.FALSE.
      polar_allbeads=.TRUE.
      PITIMER=.FALSE.
      cay_correction=.FALSE.
      save_beads="ALL"
      add_buffer_pi=.FALSE.
      lambda_trpmd=1.0d0
      default_lambda_trpmd=.TRUE.

      do i = 1, nkey
        ios = 0
        next = 1
        record = keyline(i)
        call gettext (record,keyword,next)
        call upcase (keyword)
        string = record(next:240)
        
        select case (trim(keyword))
        case ('NBEADS')
          read (string,*,iostat=ios) nbeads
        case ('NBEADS_CTR')
          read (string,*,iostat=ios) nbeads_ctr
        case ('PI_COMM_REPORT')
          pi_comm_report=.TRUE.
        case ('SAVE_BEADS')
          call getword (record,save_beads,next)
          call upcase (save_beads)
        case ('PITIMER')
          PITIMER=.TRUE.
        case ('PI_START_ISOBARIC')
          read (string,*,iostat=ios) pi_start_isobaric
        case ('POLAR_CENTROID')
          polar_allbeads = .FALSE.
        case ('CENTROID_LONGRANGE')
          centroid_longrange = .TRUE.
        case ('CAY_CORR')
          cay_correction = .TRUE.
        case ('ADD_BUFFER_PI')
          add_buffer_pi=.TRUE.
        case ('LAMBDA_TRPMD')
          read (string,*,iostat=ios) lambda_trpmd
          default_lambda_trpmd=.FALSE.
        end select

        if (ios /= 0) then
          write(*,*) "Warning: keyword ",trim(keyword)
     &         ," not correctly read!"
        endif
      end do

      
      if(nbeads_ctr>nbeads) then
        write(0,*) "INPUT ERROR: nbeads_ctr must be lower than nbeads"
        call fatal
      endif

      contract = nbeads_ctr > 0 .AND. nbeads_ctr<nbeads
      centroid_recip=(contract.and.nbeads_ctr==1)
     &                     .OR. centroid_longrange

      if(.not. contract) nbeads_ctr = 0

      end subroutine read_bead_config


      end module
