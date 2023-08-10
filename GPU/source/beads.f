#include "tinker_precision.h"
#include "tinker_types.h"
      submodule(beads) sub_beads

      contains

      module subroutine allocate_polymer(polymer,nu,allocate_vel
     &      ,allocate_forces,allocate_slow)
        use atoms, only: n
        use domdec, only: nproc_polymer
        use tinheader, only: re_p
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
        integer, intent(in) :: nu
        LOGICAL, intent(in) :: allocate_vel, allocate_slow
        LOGICAL, intent(in) :: allocate_forces
        integer :: i

        call deallocate_polymer(polymer)
!$acc enter data create(polymer) async

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
!$acc enter data copyin(polymer%ibead_beg,polymer%ibead_end
!$acc&        ,polymer%nbeadsloc) async
        if(nproc_polymer>1) then
          allocate(polymer%iproc(nu))
          do i=1,nu
            polymer%iproc(i)=
     &         min(nproc_polymer,int((i-1)/polymer%nbeadsproc)+1)
          enddo
!$acc enter data copyin(polymer%iproc) async
        endif

        allocate(polymer%epot(nu))
        allocate(polymer%vir(3,3,nu))
        allocate(polymer%dip(3,nu))
        allocate(polymer%dipind(3,nu))
        polymer%epot(:)=0.0_re_p
        polymer%vir(:,:,:)=0.0_re_p
        polymer%dip(:,:)=0.0_re_p
        polymer%dipind(:,:)=0.0_re_p

        allocate(polymer%pos(3,n,nu))
!$acc enter data async create(polymer%pos) 
!$acc&    copyin(polymer%epot,polymer%vir
!$acc&          ,polymer%dip, polymer%dipind )

        
        if(allocate_vel) then
          allocate(polymer%vel(3,n,nu))
!$acc enter data create(polymer%vel) async
        endif
        
        if(allocate_forces) then
          allocate(polymer%forces(3,n,nu))
!$acc enter data create(polymer%forces) async
        endif
        
        if(allocate_slow) then
          allocate(polymer%forces_slow(3,n,nu))
!$acc enter data create(polymer%forces_slow) async
        endif

!$acc wait

        polymer%allocated=.true.

      end subroutine allocate_polymer

      module subroutine deallocate_polymer(polymer)
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer

        if(.not. polymer%allocated) return

        polymer%allocated=.false.

        if(allocated(polymer%pos)) then
!$acc exit data delete(polymer%pos) async
          deallocate(polymer%pos)
        endif
        if(allocated(polymer%eigpos)) then
!$acc exit data delete(polymer%eigpos) async
          deallocate(polymer%eigpos)
        endif
        if(allocated(polymer%vel)) then
!$acc exit data delete(polymer%vel) async
          deallocate(polymer%vel)
        endif
        if(allocated(polymer%eigvel)) then
!$acc exit data delete(polymer%eigvel) async
          deallocate(polymer%eigvel)
        endif
        if(allocated(polymer%forces)) then
!$acc exit data delete(polymer%forces) async
          deallocate(polymer%forces)
        endif  
        if(allocated(polymer%eigforces)) then
!$acc exit data delete(polymer%eigforces) async
          deallocate(polymer%eigforces)
        endif  

        if(allocated(polymer%forces_slow)) then
!$acc exit data delete(polymer%forces_slow) async
          deallocate(polymer%forces_slow)
        endif  

        if(allocated(polymer%epot)) then
!$acc exit data delete(polymer%epot) async
          deallocate(polymer%epot)
        endif      
        
        if(allocated(polymer%vir)) then
!$acc exit data delete(polymer%vir) async
          deallocate(polymer%vir)
        endif  

        if(allocated(polymer%dip)) then
!$acc exit data delete(polymer%dip) async
          deallocate(polymer%dip)
        endif

        if(allocated(polymer%dipind)) then
!$acc exit data delete(polymer%dipind) async
          deallocate(polymer%dipind)
        endif 

        if(allocated(polymer%nbeadsloc)) then
!$acc exit data delete(polymer%nbeadsloc) async
          deallocate(polymer%nbeadsloc)
        endif 

        if(allocated(polymer%ibead_beg)) then
!$acc exit data delete(polymer%ibead_beg) async
          deallocate(polymer%ibead_beg)
        endif 

        if(allocated(polymer%ibead_end)) then
!$acc exit data delete(polymer%ibead_end) async
          deallocate(polymer%ibead_end)
        endif 

        if(allocated(polymer%iproc)) then
!$acc exit data delete(polymer%iproc) async
          deallocate(polymer%iproc)
        endif 

!$acc exit data delete(polymer) async
!$acc wait
      end subroutine deallocate_polymer

      module subroutine initialize_normal_modes(polymer,verbose)
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
!$acc exit data delete(polymer%eigpos,polymer%eigvel) async
        deallocate(polymer%eigvel,polymer%eigpos)
      endif
      allocate(polymer%eigpos(3,n,nu))
      allocate(polymer%eigvel(3,n,nu))
!$acc enter data create(polymer%eigpos,polymer%eigvel) async


      if(allocated(polymer%eigforces)) then
!$acc exit data delete(polymer%eigforces) async
        deallocate(polymer%eigforces)
      endif
      allocate(polymer%eigforces(3,n,nu))
!$acc enter data create(polymer%eigforces) async

!$acc parallel loop collapse(3) async default(present)
      do k=1,nu; do i=1,n; do j=1,3
        polymer%eigpos(j,i,k)=0.d0
        polymer%eigvel(j,i,k)=0.d0
        polymer%eigforces(j,i,k)=0.d0
      enddo ;enddo; enddo

!$acc wait

      if(ranktot==0) write(0,*) "normal modes initialized"
      end subroutine initialize_normal_modes

      module subroutine update_nlocpi(nloc)
        use domdec, only: nproc_polymer,rank_polymer
        implicit none
        integer, intent(in) :: nloc
        integer :: i

        if(.not.allocated(ilocpi_beg)) then
          allocate(ilocpi_beg(nproc_polymer))
          allocate(ilocpi_end(nproc_polymer))
          allocate(nlocpis(nproc_polymer))
!$acc enter data create(ilocpi_beg,ilocpi_end,nlocpis) async
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
!$acc update device(ilocpi_beg,ilocpi_end,nlocpis) async

      end subroutine update_nlocpi

      module subroutine allocpi()
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
      real(8), allocatable :: WORK(:)
      real(8), allocatable :: WORK_ctr(:),eigMat_ctr(:,:),omkpi_ctr(:) 

      call update_nlocpi(nloc)
!$acc wait      
      if (allocated(eigmat)) then
!$acc exit data delete(eigmat)
          deallocate(eigmat)
      end if
      allocate(eigmat(nbeads,nbeads))
!$acc enter data create(eigmat)

      if (allocated(omkpi)) then
!$acc exit data delete(omkpi)
          deallocate(omkpi)
      end if
      allocate(omkpi(nbeads))        
!$acc enter data create(omkpi)

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
        !if(eigmat(1,1)<0) eigmat=-eigmat
        deallocate(WORK)
      endif
      call MPI_BCAST(eigmat,nbeads**2,MPI_REAL8
     &    ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(omkpi,nbeads,MPI_REAL8
     &    ,0,MPI_COMM_WORLD,ierr)
!$acc update device(eigmat,omkpi)


      if (contract .and. nbeads_ctr>1) then
        if (allocated(contractor_mat)) then
!$acc exit data delete(contractor_mat)
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
          !if(eigmat_ctr(1,1)<0) eigmat_ctr=-eigmat_ctr    
          
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
!$acc enter data copyin(contractor_mat)
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

      module subroutine update_normal_modes_pi(polymer)
      use atoms, only: n
      use domdec
      use inform, only: deb_Path
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu
      real(r_p) :: sqrtnuinv

        if (deb_Path) then
!$acc wait
          write(*,*) '   >> update_normal_modes_pi'
        endif
        nu=polymer%nbeads
        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)      
        sqrtnuinv = 1./sqrt(real(nu,r_p))
!$acc parallel loop collapse(3) async default(present)
        DO k=1,nu; DO iloc=ilocbeg,ilocend ; DO j=1,3
          i=glob(iloc)
          polymer%eigpos(j,i,k) = 0.
          polymer%eigvel(j,i,k) = 0.
!$acc loop seq
          do ibead=1,nu
            polymer%eigpos(j,i,k) = polymer%eigpos(j,i,k) 
     &          + eigmat(ibead,k)*polymer%pos(j,i,ibead)
            polymer%eigvel(j,i,k) = polymer%eigvel(j,i,k) 
     &          + eigmat(ibead,k)*polymer%vel(j,i,ibead)
          enddo
          polymer%eigpos(j,i,k) = polymer%eigpos(j,i,k)
     &                               *sqrtnuinv
          polymer%eigvel(j,i,k) = polymer%eigvel(j,i,k)
     &                               *sqrtnuinv
        ENDDO ; ENDDO ; ENDDO 

        if (deb_Path) then
!$acc wait
          write(*,*) '   << update_normal_modes_pi'
        endif


      end subroutine update_normal_modes_pi

      module subroutine set_eigforces_pi(polymer,forces)
      use atoms, only: n
      use domdec
      use inform, only: deb_Path
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real(r_p), intent(in) :: forces(:,:,:)
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu
      real(r_p) :: sqrtnuinv

        if (deb_Path) then
!$acc wait
          write(*,*) '   >> set_eigforces_pi'
        endif
        nu=polymer%nbeads
        sqrtnuinv = 1./sqrt(real(nu,r_p))
        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)      
!$acc parallel loop collapse(3) async default(present)
        DO k=1,nu; DO iloc=ilocbeg,ilocend ; DO j=1,3
          i=glob(iloc)
          polymer%eigforces(j,i,k) = 0.
!$acc loop seq
          do ibead=1,nu
            polymer%eigforces(j,i,k) = polymer%eigforces(j,i,k) 
     &          + eigmat(ibead,k)*forces(j,i,ibead)*sqrtnuinv
          enddo                           
        ENDDO ; ENDDO ; ENDDO 

        if (deb_Path) then
!$acc wait
          write(*,*) '   << set_eigforces_pi'
        endif

      end subroutine set_eigforces_pi

      module subroutine update_direct_space_pi(polymer)
      use atoms, only: n
      use domdec
      use inform, only: deb_Path
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu
      real(r_p) :: sqrtnu

        if (deb_Path) then
!$acc wait
          write(*,*) '   >> update_direct_space_pi'
        endif

        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)
        nu=polymer%nbeads
        sqrtnu=sqrt(real(nu,r_p))
!$acc parallel loop collapse(3) async default(present)
        DO ibead=1,nu; DO iloc=ilocbeg,ilocend ; DO j=1,3
          i=glob(iloc)
          polymer%pos(j,i,ibead) = 0.
          polymer%vel(j,i,ibead) = 0.
!$acc loop seq
          do k=1,nu
            polymer%pos(j,i,ibead) = polymer%pos(j,i,ibead) 
     &          + eigmat(ibead,k)*polymer%eigpos(j,i,k)
            polymer%vel(j,i,ibead) = polymer%vel(j,i,ibead) 
     &          + eigmat(ibead,k)*polymer%eigvel(j,i,k)
          enddo
          polymer%pos(j,i,ibead) = polymer%pos(j,i,ibead)
     &                               *sqrtnu
          polymer%vel(j,i,ibead) = polymer%vel(j,i,ibead)
     &                               *sqrtnu
        ENDDO ; ENDDO ; ENDDO

        if (deb_Path) then
!$acc wait
          write(*,*) '   << update_direct_space_pi'
        endif
      end subroutine update_direct_space_pi

      module subroutine contract_polymer(polymer,polymer_ctr)
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
!$acc parallel loop collapse(3) async default(present)
        DO ibead=1,nuctr; DO iloc=ilocbeg,ilocend ; DO j=1,3
          i=glob(iloc)
          polymer_ctr%pos(j,i,ibead) = 0.d0 
!$acc loop seq
          do k=1,nu
            polymer_ctr%pos(j,i,ibead) = polymer_ctr%pos(j,i,ibead) 
     &          + contractor_mat(ibead,k)*polymer%pos(j,i,k)
          enddo
        ENDDO ; ENDDO ;ENDDO
      end subroutine contract_polymer

      module subroutine project_forces_ctr(polymer,polymer_ctr,slow)
      use atoms, only: n
      use domdec
      use inform, only: deb_Path
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer,polymer_ctr
      logical, intent(in), optional :: slow
      integer :: i,j,k,ibead,iloc,ilocbeg,ilocend,nu,nuctr
      real(r_p) :: nuratio
      logical :: slow_

        slow_ = .false.
        if(present(slow)) slow_=slow

        if (deb_Path) write(*,*) '   >> project_forces_ctr'
        ilocbeg=ilocpi_beg(rank_polymer+1)
        ilocend=ilocpi_end(rank_polymer+1)
        nu=polymer%nbeads
        nuctr=polymer_ctr%nbeads
        nuratio=real(nu,r_p)/real(nuctr,r_p)
        if(slow_) then
!$acc parallel loop collapse(3) async default(present)
          DO ibead=1,nu; DO iloc=ilocbeg,ilocend ; DO j=1,3
            i=glob(iloc)
!$acc loop seq
            do k=1,nuctr
              polymer%forces_slow(j,i,ibead)=
     &           polymer%forces_slow(j,i,ibead)
     &           + contractor_mat(k,ibead)*nuratio
     &                *polymer_ctr%forces_slow(j,i,k)
            enddo
          ENDDO ; ENDDO ;ENDDO
        else
!$acc parallel loop collapse(3) async default(present)
          DO ibead=1,nu; DO iloc=ilocbeg,ilocend ; DO j=1,3
            i=glob(iloc)
!$acc loop seq
            do k=1,nuctr
              polymer%forces(j,i,ibead)=polymer%forces(j,i,ibead)
     &           + contractor_mat(k,ibead)*nuratio
     &                *polymer_ctr%forces_slow(j,i,k)
            enddo
          ENDDO ; ENDDO ;ENDDO
        endif

        if (deb_Path) write(*,*) '   << project_forces_ctr'
      end subroutine project_forces_ctr

      module subroutine load_bead(polymer,ibead,force_load)
        use atomsMirror
        use domdec
        use inform, only: deb_Path
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
        if(deb_path) then
!$acc wait
          write(*,*) '  load_bead',ibead
        endif

        ibead_loaded = ibead
        ipolymer_loaded = polymer%index
        if(ibead == 0) then
!$acc parallel loop async default(present)
          DO iloc=1,nbloc
            i=glob(iloc)
            x(i)=polymer%eigpos(1,i,1)
            y(i)=polymer%eigpos(2,i,1)
            z(i)=polymer%eigpos(3,i,1)
          ENDDO
        else
!$acc parallel loop async default(present)
          DO iloc=1,nbloc
            i=glob(iloc)
            x(i)=polymer%pos(1,i,ibead)
            y(i)=polymer%pos(2,i,ibead)
            z(i)=polymer%pos(3,i,ibead)
          ENDDO
        endif

        call reCast_position


      end subroutine load_bead

      module subroutine stopwatchpi(timer,barrier,reset)
      use mpi
      implicit none
      real(8), intent(inout) :: timer
      logical, intent(in) :: barrier,reset
      real(8), save :: time0, time1
      integer :: ierr

!$acc wait
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

      module subroutine kinetic_pi_fast(polymer)
      use units
      use atoms, only: n
      use mdstuf
      use domdec
      use atmtyp
      use mpi
      use bath, only: kelvin
      use inform, only: deb_Path
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real(r_p) :: factor1,rnu,wnu2,ek0,hbar2beta,lth2
      integer :: i,j,k,ibead,iloc,ilocbeg,nu
      real(r_p),save:: ekin11,ekin12,ekin13
      real(r_p),save:: ekin21,ekin22,ekin23
      real(r_p),save:: ekin31,ekin32,ekin33
      real(r_p),term
      logical,save::f_in=.true.

      
      nu = polymer%nbeads
      if(nu==1) return

      if (deb_Path) then
!$acc wait
        write(*,*) '   >> kinetic_pi_fast'
      endif

      if (f_in) then
         f_in=.false.
         ekin11=0;ekin22=0;ekin33=0;
         ekin21=0;ekin31=0;ekin12=0;
         ekin32=0;ekin13=0;ekin23=0;
!$acc enter data copyin(ekin11,ekin22,ekin33
!$acc&        ,ekin21,ekin31,ekin12,ekin32,ekin13,ekin23)
      end if

!$acc data present(ekvir_fast
!$acc&   ,ekin11,ekin22,ekin33
!$acc&   ,ekin21,ekin31,ekin12,ekin32,ekin13,ekin23
!$acc&   , ekinpi_fast) async
      ilocbeg = ilocpi_beg(rank_polymer+1)
!$acc parallel loop collapse(3) async default(present)
      DO k=2,nu; DO iloc=1,nlocpi ; DO j=1,3
        i=glob(iloc+ilocbeg-1)   
        ekvir_fast = ekvir_fast + polymer%eigpos(j,i,k)
     &              *polymer%eigforces(j,i,k)
        if (j.eq.1) then
          term=polymer%eigpos(1,i,k)
          ekin11=ekin11+term*polymer%eigforces(1,i,k)
          ekin21=ekin21+term*polymer%eigforces(2,i,k)
          ekin31=ekin31+term*polymer%eigforces(3,i,k)
        else if (j.eq.2) then
          term=polymer%eigpos(2,i,k)
          ekin12=ekin12+term*polymer%eigforces(1,i,k)
          ekin22=ekin22+term*polymer%eigforces(2,i,k)
          ekin32=ekin32+term*polymer%eigforces(3,i,k)
        else
          term=polymer%eigpos(3,i,k)
          ekin13=ekin13+term*polymer%eigforces(1,i,k)
          ekin23=ekin23+term*polymer%eigforces(2,i,k)
          ekin33=ekin33+term*polymer%eigforces(3,i,k)
        end if
      ENDDO ; ENDDO ; ENDDO  

!$acc serial async default(present)
        ekinpi_fast(1,1) = ekin11
        ekinpi_fast(2,1) = ekin21
        ekinpi_fast(3,1) = ekin31
        ekinpi_fast(1,2) = ekin12
        ekinpi_fast(2,2) = ekin22
        ekinpi_fast(3,2) = ekin32
        ekinpi_fast(1,3) = ekin13
        ekinpi_fast(2,3) = ekin23
        ekinpi_fast(3,3) = ekin33

        !---zero out the total kinetic energy and its outer product
        ekin11 = 0.0_re_p
        ekin22 = 0.0_re_p
        ekin33 = 0.0_re_p
        ekin21 = 0.0_re_p
        ekin31 = 0.0_re_p
        ekin12 = 0.0_re_p
        ekin32 = 0.0_re_p
        ekin13 = 0.0_re_p
        ekin23 = 0.0_re_p
!$acc end serial 
!$acc end data

      if (deb_Path) then
!$acc wait
        write(*,*) '   << kinetic_pi_fast'
      endif

      end subroutine kinetic_pi_fast


      module subroutine kinetic_pi(polymer)
      use units
      use atoms, only: n
      use mdstuf
      use domdec
      use atmtyp
      use mpi
      use bath, only: kelvin
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real(r_p) :: factor1,rnu,wnu2,ek0,hbar2beta,lth2
      integer :: i,j,k,ibead,iloc,ilocbeg,nu
      real(r_p),save:: ekin11,ekin12,ekin13
      real(r_p),save:: ekin21,ekin22,ekin23
      real(r_p),save:: ekin31,ekin32,ekin33
      real(r_p) :: term
      logical,save::f_in=.true.

      if (f_in) then
         f_in=.false.
         ekin11=0;ekin22=0;ekin33=0;
         ekin21=0;ekin31=0;ekin12=0;
         ekin32=0;ekin13=0;ekin23=0;
!$acc enter data copyin(ekin11,ekin22,ekin33
!$acc&        ,ekin21,ekin31,ekin12,ekin32,ekin13,ekin23)
      end if


!$acc data present(Ekcentroid,ekdynpi,ekprim,ekvir
!$acc&   ,temppi,temppi_cl,temppi_centroid,ekvir_fast,gyrpi
!$acc&   ,ekin11,ekin22,ekin33
!$acc&   ,ekin21,ekin31,ekin12,ekin32,ekin13,ekin23
!$acc&   ,ekinpi, ekinpi_fast) async
      nu = polymer%nbeads
      rnu = real(nu,r_p)
      factor1=0.5d0/convert
      wnu2=(rnu*boltzmann*kelvin/hbar_planck)**2
      ek0 = n*boltzmann*kelvin*factor1

      hbar2beta=hbar_planck**2/(boltzmann*kelvin)

!$acc serial async
      ekprim=0.d0
      ekvir=0.d0
      ekdynpi=0.d0
      Ekcentroid=0.d0
      gyrpi=0.d0
!$acc end serial  


      ilocbeg = ilocpi_beg(rank_polymer+1)
!$acc parallel loop collapse(3) async default(present)
      DO k=1,nu; DO iloc=1,nlocpi ; DO j=1,3
        i=glob(iloc+ilocbeg-1)   
        lth2=hbar2beta/mass(i)
        gyrpi=gyrpi+
     &     (polymer%pos(j,i,k)-polymer%eigpos(j,i,1))**2/lth2
        if(k==1) then
        Ekcentroid=Ekcentroid+mass(i)*polymer%eigvel(j,i,1)**2
        else
        ekdynpi = ekdynpi + mass(i)*polymer%eigvel(j,i,k)**2
        ekvir = ekvir + polymer%eigpos(j,i,k)
     &              *polymer%eigforces(j,i,k)
        ekprim = ekprim + mass(i)*( omkpi(k)
     &             *polymer%eigpos(j,i,k) )**2
        if (j.eq.1) then
          term=polymer%eigpos(1,i,k)
          ekin11=ekin11+term*polymer%eigforces(1,i,k)
          ekin21=ekin21+term*polymer%eigforces(2,i,k)
          ekin31=ekin31+term*polymer%eigforces(3,i,k)
        else if (j.eq.2) then
          term=polymer%eigpos(2,i,k)
          ekin12=ekin12+term*polymer%eigforces(1,i,k)
          ekin22=ekin22+term*polymer%eigforces(2,i,k)
          ekin32=ekin32+term*polymer%eigforces(3,i,k)
        else
          term=polymer%eigpos(3,i,k)
          ekin13=ekin13+term*polymer%eigforces(1,i,k)
          ekin23=ekin23+term*polymer%eigforces(2,i,k)
          ekin33=ekin33+term*polymer%eigforces(3,i,k)
        end if
        endif
      ENDDO ; ENDDO ; ENDDO  

      if(nu==1) then
!$acc serial async
        gyrpi=0.d0
        Ekcentroid = Ekcentroid*factor1
        ekdynpi = ekcentroid
        ekprim = ekcentroid
        ekvir = Ekcentroid
        ekvir_fast = 0.0_re_p  
!$acc end serial 
      else
!$acc serial async
        ekinpi(1,1) = ek0 - 0.5*(ekin11+ekinpi_fast(1,1))
        ekinpi(2,1) =     - 0.5*(ekin21+ekinpi_fast(2,1))
        ekinpi(3,1) =     - 0.5*(ekin31+ekinpi_fast(3,1))
        ekinpi(1,2) =     - 0.5*(ekin12+ekinpi_fast(1,2))
        ekinpi(2,2) = ek0 - 0.5*(ekin22+ekinpi_fast(2,2))
        ekinpi(3,2) =     - 0.5*(ekin32+ekinpi_fast(3,2))
        ekinpi(1,3) =     - 0.5*(ekin13+ekinpi_fast(1,3))
        ekinpi(2,3) =     - 0.5*(ekin23+ekinpi_fast(2,3))
        ekinpi(3,3) = ek0 - 0.5*(ekin33+ekinpi_fast(3,3))

        ekinpi(2,1)=0.5*(ekinpi(2,1)+ekinpi(1,2))
        ekinpi(1,2)=ekinpi(2,1)
        ekinpi(3,1)=0.5*(ekinpi(3,1)+ekinpi(1,3))
        ekinpi(1,3)=ekinpi(3,1)
        ekinpi(3,2)=0.5*(ekinpi(3,2)+ekinpi(2,3))
        ekinpi(2,3)=ekinpi(3,2)

        !---zero out the total kinetic energy and its outer product
        ekin11 = 0.0_re_p
        ekin22 = 0.0_re_p
        ekin33 = 0.0_re_p
        ekin21 = 0.0_re_p
        ekin31 = 0.0_re_p
        ekin12 = 0.0_re_p
        ekin32 = 0.0_re_p
        ekin13 = 0.0_re_p
        ekin23 = 0.0_re_p
        ekinpi_fast(1,1) = 0.0_re_p
        ekinpi_fast(2,1) = 0.0_re_p
        ekinpi_fast(3,1) = 0.0_re_p
        ekinpi_fast(1,2) = 0.0_re_p
        ekinpi_fast(2,2) = 0.0_re_p
        ekinpi_fast(3,2) = 0.0_re_p
        ekinpi_fast(1,3) = 0.0_re_p
        ekinpi_fast(2,3) = 0.0_re_p
        ekinpi_fast(3,3) = 0.0_re_p

        gyrpi=sqrt(hbar2beta*gyrpi/real(3*n*nu,r_p))
        Ekcentroid = Ekcentroid*factor1
        ekdynpi = ekdynpi*factor1 + ekcentroid
        ekprim = ekdynpi - ekprim*factor1 !/rnu
        ekvir = Ekcentroid - 0.5d0*(ekvir +ekvir_fast) !/rnu
!        ekprim = rnu*ek0 - ekprim*factor1
!        ekvir =  ek0 - 0.5d0*(ekvir+ekvir_fast) 
        ekvir_fast = 0.0_re_p  
!$acc end serial 

      endif

!$acc end data


      end subroutine kinetic_pi

      module subroutine reduce_observables_pi(polymer)
      use atoms, only :n
      use units
      use mdstuf
      use domdec
      use atmtyp
      use mpi
      use boxes
      use bath, only:anisotrop
      implicit none
      type(POLYMER_COMM_TYPE), intent(in) :: polymer
      logical, save:: f_in=.true.
      integer :: i,j,k,ibead,iloc,ilocbeg,ierr
      real(r_p) :: factor2,factor
      integer,save :: ncomm

      factor2=2.d0/(3*n*gasconst)
      factor=prescon/volbox

      if(nproctot==1) then
!$acc serial async present(temppi,temppi_cl,ekvir,ekdynpi
!$acc&      ,dedvpi,virpi,presvir,stresspi,ekinpi)
        temppi = factor2 * ekvir
        temppi_cl = factor2 *ekdynpi 
        temppi_centroid = factor2*ekcentroid
        dedvpi = (virpi(1,1)+virpi(2,2)+virpi(3,3))/(3.d0*volbox)
        presvir = prescon*( -dedvpi + 2.d0*ekvir
     &                              /(3.d0*volbox) )
        stresspi(:,:)=factor*(2*ekinpi(:,:)-virpi(:,:))
!$acc end serial
        return
      endif

      if(f_in) then
        if (anisotrop) then
          ncomm=12
        else
          ncomm=6
        endif
        allocate(bufferpi(ncomm))
!$acc enter data create(bufferpi) async
        f_in = .false.
      endif

!$acc serial async present(ekdynpi,ekprim,ekvir
!$acc& ,Ekcentroid,epotpi,virpi,volbox)
        bufferpi(1)=ekdynpi
        bufferpi(2)=ekcentroid
        bufferpi(3)=ekprim
        bufferpi(4)=ekvir
        bufferpi(5)=epotpi
        bufferpi(6)=(virpi(1,1)+virpi(2,2)+virpi(3,3))/(3.d0*volbox)
        if (anisotrop) then
          bufferpi(7)=stresspi(1,1)
          bufferpi(8)=stresspi(2,2)
          bufferpi(9)=stresspi(3,3)
          bufferpi(10)=0.5*(stresspi(1,2)+stresspi(2,1))
          bufferpi(11)=0.5*(stresspi(1,3)+stresspi(3,1))
          bufferpi(12)=0.5*(stresspi(2,3)+stresspi(3,2))
        endif
!$acc end serial
!$acc wait
!$acc host_data use_device(bufferpi)
        if(ranktot==0) then
          call MPI_REDUCE(MPI_IN_PLACE,bufferpi
     &       ,ncomm,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
          call MPI_REDUCE(bufferpi,bufferpi
     &        ,ncomm,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
!$acc end host_data

        if(ranktot==0) then
!$acc serial async present(ekdynpi,ekprim,ekvir
!$acc& ,Ekcentroid,epotpi,virpi,volbox
!$acc&  ,temppi,temppi_cl,temppi_centroid,dedvpi,presvir
!$acc&  ,stresspi)
          ekdynpi=bufferpi(1)
          ekcentroid=bufferpi(2)
          ekprim=bufferpi(3)
          ekvir=bufferpi(4)
          epotpi=bufferpi(5)
          dedvpi=bufferpi(6)
          temppi = factor2 * ekvir
          temppi_cl = factor2 *ekdynpi
          temppi_centroid = factor2*ekcentroid  
          presvir = prescon*( -dedvpi + 2.d0*ekvir
     &                              /(3.d0*volbox) )
          if(anisotrop) then
            stresspi(1,1)=bufferpi(7)
            stresspi(2,2)=bufferpi(8)
            stresspi(3,3)=bufferpi(9)
            stresspi(1,2)=bufferpi(10)
            stresspi(1,3)=bufferpi(11)
            stresspi(2,3)=bufferpi(12)
            stresspi(2,1)=bufferpi(10)
            stresspi(3,1)=bufferpi(11)
            stresspi(3,2)=bufferpi(12)
          endif
!$acc end serial

        endif

      end subroutine reduce_observables_pi

      module subroutine reset_observables(init)
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

        if(init) then
!$acc enter data copyin(epotpi_loc,eksumpi_loc,etotpi_loc,ekinpi_loc
!$acc&           ,eintrapi_loc,einterpi_loc,temppi,temppi_cl
!$acc&           ,eintrapi,einterpi,epotpi,temppi_centroid
!$acc&           ,dedvintrapi,dedvinterpi,dedvpi
!$acc&           ,Ekcentroid,ekdynpi,ekvir,ekprim
!$acc&            ,presvir,ekvir_fast,virpi,ekinpi,ekinpi_fast,stresspi
!$acc&            ,dippi, dipindpi,gyrpi)
        else
!$acc update device(epotpi_loc,eksumpi_loc,etotpi_loc,ekinpi_loc
!$acc&           ,eintrapi_loc,einterpi_loc,temppi,temppi_cl
!$acc&           ,eintrapi,einterpi,epotpi,temppi_centroid
!$acc&           ,dedvintrapi,dedvinterpi,dedvpi
!$acc&           ,Ekcentroid,ekdynpi,ekvir,ekprim
!$acc&            ,presvir,ekvir_fast,virpi,ekinpi,ekinpi_fast,stresspi
!$acc&            ,dippi, dipindpi, gyrpi)
        endif

      end subroutine reset_observables

      module subroutine read_bead_config()
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
      if (contract .and. nbeads_ctr==1) then
        centroid_longrange=.FALSE.
        centroid_recip=.TRUE.
      elseif(centroid_longrange) then
        centroid_recip=.TRUE.
      endif

      if(.not. contract) nbeads_ctr = 0

      end subroutine read_bead_config

      end submodule