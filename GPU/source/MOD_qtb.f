c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ########################################################################
c     ##                                                                    ##
c     ##  module QTB  --  parameters and arrays for Quantum Thermal Bath    ##
c     ##                                                                    ##
c     ########################################################################
c     omegacut : cutting frequency for the QTB random force
c     domega : frequency discretization in THz (ps-1)
c     omegasmear : length of the smearing window near omegacut
c     skipseg : number to seg to be skipped (because the system has not
c     vad : past velocities used for correlation computations
c
c     n_of_type: number of atoms of each type
c     adqtb_type: map the atoms to their adqtb type
c     ntype: number of different types
c
c     Literature reference: 
c     J-L Barrat & D. Rodney JStatPhys (2011)
c     and H Dammak PRL 103, 190601 (2009)
c
#include "tinker_macro.h"
      module qtb 
      use spectra, only: compteur_type
      implicit none
c      integer, allocatable :: repartnoise(:)
      integer :: skipseg
      real(r_p) :: omegasmear,gammar_min
      real(r_p), allocatable :: Htilde(:,:,:)
      real(r_p), allocatable :: corr_pot_ratio(:,:)
      logical :: noQTB
      logical :: corr_pot
      logical :: register_spectra
      logical :: qtb_thermostat=.FALSE.
      logical :: qtb_separate_fft=.FALSE.
      logical :: QTB_verbose=.FALSE.
      integer :: qtb_batch_size
      type(compteur_type) :: compteur, compteur_bis
      integer :: plan_qtb, nloc_plan=-1,nfft_plan=-1
      character*11 :: adqtb_optimizer
      character*11,allocatable :: adqtb_optimizer_type(:)

      integer :: adqtb_avg_adapt=1
      real(r_p) :: adqtb_tau_avg_default, adqtb_tau_adapt
      logical :: adqtb_correct_bias, adqtb_win_corr
      real(r_p),allocatable :: adqtb_tau_avg(:),A_gamma(:)

      logical :: isobaric_save,use_piston_save
      logical :: piqtb
c
c     adQTB specific variables
c
      logical adaptive_qtb
      real(r_p) :: A_gamma_default
      real(r_p) :: corr_fact_qtb(3)
      logical :: use_corr_fact_qtb
      logical :: corr_fact_converged
      real(r_p), allocatable :: gamma_type(:,:,:)
      real(r_p), allocatable :: gamma_type_init(:,:,:)
      real(r_p), allocatable :: mCvv_average_type(:,:,:)
      real(r_p), allocatable :: mCvv_average_sum(:,:,:)
      real(r_p), allocatable :: Cvf_average_type(:,:,:)
      real(r_p), allocatable :: mCvv_mean(:,:,:)
      real(r_p), allocatable :: Cvf_mean(:,:,:)
      real(r_p), allocatable :: dFDR_mean(:,:,:)
      real(r_p), allocatable :: dFDR_average_type(:,:,:)
      integer, allocatable :: adqtb_type(:),n_of_type(:)
      integer, allocatable :: adqtb_type_map(:)
      character*3, allocatable :: adqtb_name(:)
      integer :: ntype
      logical :: save_gamma_history=.FALSE.

      logical :: adqtb_start_classical=.FALSE.
      logical :: corr_pot_corr_fact=.TRUE.

      real(r_p) :: adqtb_smooth_default
      real(r_p), allocatable :: adqtb_smooth_type(:)

      abstract interface
          subroutine qtb_optimizer_type(itype,mCvv,Cvf,dFDR
     &   ,qtb_therm)
        implicit none
        integer, intent(in) :: itype
        real(r_p), intent(in) :: dFDR(:,:)
        real(r_p), intent(in) :: mCvv(:,:),Cvf(:,:)
        logical,intent(in) :: qtb_therm
        end subroutine qtb_optimizer_type
      end interface

      type qtbdat_t
        real(4), allocatable :: noise(:,:,:)
        real(4), allocatable :: fad(:,:,:), vad(:,:,:)
c        logical, allocatable :: proc_has_data(:)
c        logical :: currently_loc
      end type 

      type qtbopt_t
        procedure(qtb_optimizer_type), pointer, nopass :: p
      end type

      type(qtbdat_t), allocatable :: qtbdata(:)
      logical, allocatable :: qtbdata_allocated(:)
      type(qtbopt_t), allocatable :: adapt_gamma(:)

      save 
      private
      public :: qtbinit,qtbrandom
     &          ,adaptive_qtb
     &          , A_gamma_default, noQTB, corr_pot, skipseg
     &          , corr_fact_qtb, register_spectra
     &          ,save_gamma_history, qtb_thermostat
     &          ,reassignqtb
     &          ,qtb_batch_size, use_corr_fact_qtb
     &          ,adqtb_start_classical
     &          ,adqtb_optimizer,adqtb_correct_bias
     &          ,adqtb_tau_avg_default, adqtb_tau_adapt
     &          ,adqtb_avg_adapt
     &          ,adqtb_smooth_default
     &          ,qtbdata, QTB_update_state
     &          ,adqtb_win_corr,corr_pot_corr_fact
     &          ,get_qtb_compteur, piqtb,QTB_verbose

      contains

c
c      subroutine qtbinit: initialize the QTB
c      
      subroutine qtbinit(dt)
      use ascii, only: int_to_str
      use atmtyp
      use atoms
      use bath
      use beads
      use domdec
      use langevin
      use math
      use moldyn
      use mdstuf, only: integrate
      use sizes
      use usage
      use units
      use mpi
      use spectra
#ifdef _OPENACC
      use cufft
#endif
      use random_mod
      use tinMemory,only:prmem_request,prmem_requestm
      implicit none
      real(r_p), intent(in) ::  dt
      integer i,j,k,l,iglob,ierr,ibead
      integer compteur_lines
      real(r_p), allocatable :: r_memory(:)
      real(r_p) :: freq,tau_tmp
      logical :: adqtb_restart,corr_pot_found
      character*120 :: gamma_restart_file
      real(r_p) :: Tseg_file, gamma_file, omegacut_file,omega
      integer, allocatable :: tmp(:)
      character :: bufchar
      logical :: corr_fact_restart
      integer :: ios,uu
      character*11 :: opttype
      character*3 :: numberbeads
      character*240 :: exten

      call read_qtb_keys()

      if(register_spectra .and. ranktot==0) 
     &   write(*,*) 'Saving QTB spectra'
      if (adaptive_qtb) register_spectra=.TRUE.
      call initialize_spectra(dt)

      if(allocated(gamma_friction)) then
!$acc exit data delete(gamma_friction)
        deallocate(gamma_friction)
      endif
      allocate(gamma_friction(nbeads))
      gamma_friction(:)=gamma

      piqtb=nbeads>1
      call set_kubo(piqtb)
      if(noQTB) kubo=.FALSE.
      if(piqtb) then
        adqtb_start_classical=.FALSE.
        if(ranktot==1 .and. omegacut <= 1.5*omkpi(nbeads)) then
          write(*,*) "Warning: omegacut should be larger than",
     &       " 1.5*omkpi(nbeads)=",1.5*omkpi(nbeads)*cm1," cm-1"
        endif
        if(default_lambda_trpmd) lambda_trpmd=0.2d0
        do k=2,nbeads
          gamma_friction(k)=sqrt(gamma**2 + (lambda_trpmd*omkpi(k))**2)
        enddo
      endif
!$acc enter data copyin(gamma_friction(:)) async

      call compteur%initialize(nseg,startsavespec)
      call compteur_bis%initialize(adqtb_avg_adapt,1)
      
      if(ranktot==0) write(*,'(A,I6,A,f5.3,A)') 
     &     "Initializing QTB with nseg=",nseg
     &      ," , dt=",dt_spec*1000.d0," fs"
      omegasmear = pi/dt/100.
      if(adqtb_tau_avg_default<=0) then
        adqtb_tau_avg_default=10*Tseg
      endif
c
c     perform dynamic allocation of some pointer arrays
c
      !call prmem_request(noise,3,nloc,2*nseg,async=.TRUE.)
      !call prmem_request(fad,3,nloc,nseg,async=.TRUE.)

      allocate(qtbdata(n))
      allocate(qtbdata_allocated(n))
!$acc update host(glob) async
!$acc enter data create(qtbdata) async
      do i=1,n
        qtbdata_allocated(i)=.FALSE.
c        qtbdata(i)%currently_loc=.FALSE.
      enddo

!$acc wait
      do i=1,nloc
        iglob=glob(i)
        call allocate_qtbdata(iglob)
c        qtbdata(iglob)%currently_loc=.TRUE.
      enddo

c      if (allocated(repartnoise)) deallocate (repartnoise)
c      allocate(repartnoise(n))

c
c     Order the atom type for the adQTB
c
      allocate(n_of_type(maxtyp))
      allocate(adqtb_type(n))    
      n_of_type(:)=0 
      do i=1,n
        k=type(i)
        if((k.gt.0).and.(atomic(i)/=0)) then
          n_of_type(k)=n_of_type(k)+1 
        endif
      enddo
      allocate(tmp(maxtyp))  
      do k=1,maxtyp
        if((n_of_type(k)>0)) then
          ntype=ntype+1
          tmp(k)=ntype
        endif
      enddo
      allocate(adqtb_name(ntype))
      allocate(adqtb_type_map(ntype))
      deallocate(n_of_type)
      allocate(n_of_type(ntype))
      adqtb_type_map = -1
      adqtb_name(:) = '   '
      n_of_type(:) = 0
      do i=1,n
        k = type(i)
        if((k.gt.0).and.(atomic(i)/=0)) then
          adqtb_type(i)=tmp(k)
          adqtb_type_map(tmp(k))=type(i)
          adqtb_name(tmp(k))=name(i)
          n_of_type(tmp(k))=n_of_type(tmp(k))+1
        endif
      enddo 
!$acc enter data copyin(adqtb_type)
      if(ranktot==0) 
     &   write(*,*) "Number of adQTB types:",ntype

      if(ranktot==0 .AND. use_corr_fact_qtb) then
        write(*,*) 'WARNING - Using Kinetic correction'
        write(*,*) 'corr_fact_qtb=', corr_fact_qtb(1)
      endif

c
c         Potential correction
c
      allocate(corr_pot_ratio(nad,nbeads))        
      corr_pot_ratio(:,:)=0d0
      if(corr_pot) then
        if(ranktot==0) write(*,*)
     &     'WARNING -  Using potential correction'
        inquire (file='corr_pot.dat', exist=corr_pot_found)
        if (corr_pot_found) then
          open(315,file='corr_pot.dat')
          read(315,*) bufchar,Tseg_file, gamma_file,omegacut_file
          if((Tseg_file.ne.Tseg).or.(omegacut_file.ne.omegacut) 
     &            .or.(gamma_file.ne.gamma)) then
                  ! NOTE: dangerous comparison of real types
c            write(*,*) 'ERROR, NOT THE SAME PARAMETERS WHILE USING
c     & THE POTENTIAL CORRECTION'
c            call fatal
            corr_pot_found=.FALSE.
          else  
            do i=1,nad
              read(315,*) omega,corr_pot_ratio(i,:)
            enddo
          endif
          close(315)
        endif
        if(.NOT.corr_pot_found) then
          call generate_corr_pot()
        endif  
      endif
c
c     Initialize gammar
c
      allocate(gamma_type(nad,ntype,nbeads))
      allocate(gamma_type_init(nad,ntype,nbeads))

      inquire (file='gamma_restart.out', exist=adqtb_restart)
      if (adqtb_restart) then
        if(ranktot==0) write(*,*) 
     &     "Using gamma_restart.out to initialize QTB"
        ! read gamma restart file
        open(newunit=uu, file='gamma_restart.out')
        read(uu,*)
        do j=1,nad
          read(uu,*) freq,gamma_type_init(j,:,:)
        enddo
        close(uu)
      else
        ! initialize gamma_type to gamma
        do k=1,nbeads
          gamma_type_init(:,:,k)=gamma_friction(k)
        enddo
        if(adqtb_start_classical .and. (.not. noQTB)) then
          if(ranktot==0) 
     &     write(*,*) "adQTB start classical"
          do i=1,nad
            gamma_type_init(i,:,:)=gamma_type_init(i,:,:)*kubo_fact(i)
          enddo
        endif
      endif
      gamma_type = gamma_type_init

      if (piqtb) then
        !!! TEMPORARY
        use_corr_fact_qtb=.TRUE.
        corr_fact_qtb(:)=1.d0
        corr_fact_converged=.TRUE.
      endif
      if(.not.use_corr_fact_qtb) then
        inquire (file='corr_fact_qtb_restart.out'
     &     , exist=corr_fact_restart)
        if (corr_fact_restart) then
          open(newunit=uu, file='corr_fact_qtb_restart.out')
          read(uu,*) corr_fact_qtb, corr_fact_converged
          close(uu)
        endif
      endif


c
c     Kernel creation
c
      call prmem_requestm(Htilde,3*nseg,ntype,nbeads,async=.TRUE.)
      call adHnoise(.TRUE.)

c
c     Generate colored noise
c
!$acc wait
#ifdef _OPENACC
      if(.not. host_rand_platform) then
        do i=1,nloc
          iglob = glob(i)
          call normalgpuR4(qtbdata(iglob)%noise(1,1,1),3*2*nseg*nbeads)
        enddo
      endif
#endif
      if (host_rand_platform) then
        do i=1,nloc
          iglob = glob(i)
          do ibead=1,nbeads; do k=1,2*nseg; do j=1,3
            qtbdata(iglob)%noise(j,k,ibead)=normal()
          enddo;enddo; enddo          
!$acc update device(qtbdata(iglob)%noise)
        enddo        
      end if
      
      if(qtb_batch_size>0) qtb_separate_fft=.TRUE.
      if(qtb_separate_fft .and. rank==0) then
        !if(qtb_batch_size<=0) then 
        !  write(0,*) "Error: QTB_BATCH_SIZE is not properly defined!"
        !  call fatal
        !endif
        write(*,*) "Warning: QTB_SEPARATE_FFT saves memory",
     &       " but may be slower (batch size: "
     &        //int_to_str(qtb_batch_size)//")"
      endif

      !call update_fft_plan_qtb(3*nseg) 
      call convolseg()
c
c     Initialisation for adQTB-r
c
      if (adaptive_qtb) then 
        call initialize_adaptation_parameters()
        gammar_min = gamma/10.
        if (adqtb_restart) then
          skipseg=1
        endif

        if(rank==0) then
          allocate(adapt_gamma(ntype))
          do i=1,ntype
            opttype = adqtb_optimizer_type(adqtb_type_map(i))
            SELECT CASE(opttype)
            CASE("SIMPLE")
              !write(0,*) "SIMPLE adQTB optimizer"
              adapt_gamma(i)%p => adapt_gamma_simple
            CASE("RATIO")            
              adapt_gamma(i)%p => adapt_gamma_ratio
            CASE DEFAULT
              write(0,*) "Error: unknown adqtb optimizer "
     &           //opttype
              call fatal
            END SELECT
          enddo
        endif
        
      endif
        
      if (register_spectra) then
        allocate(mCvv_average_type(nad,ntype,nbeads))
        allocate(mCvv_average_sum(nad,3,nbeads))
        allocate(Cvf_average_type(nad,ntype,nbeads))
        allocate(dFDR_average_type(nad,ntype,nbeads))
        allocate(mCvv_mean(nad,ntype,nbeads))
        allocate(Cvf_mean(nad,ntype,nbeads))
        allocate(dFDR_mean(nad,ntype,nbeads))
        mCvv_mean=0.d0
        Cvf_mean=0.d0
        dFDR_mean=0.d0

        if(skipseg>compteur%startsave) compteur%startsave=skipseg
      endif

      isobaric_save=isobaric
      use_piston_save=use_piston
      if(register_spectra .AND. (.NOT.
     &   (noQTB .OR. corr_fact_restart .OR.
     &   use_corr_fact_qtb .OR. corr_fact_converged ))) then   
        ! WAIT FOR corr_fact_qtb TO BE COMPUTED TO ACTIVATE PISTON   
        isobaric=.FALSE.
        use_piston=.FALSE.
        if(ranktot==0 .and. isobaric_save) then
          write(*,*) "Warning: QTB piston deactivated until "
     &          //"corr_fact_qtb is computed for the first time (at "
     &          ,real(Tseg)," ps)"  
        endif
      endif

      end subroutine qtbinit

      subroutine read_qtb_keys()
      use keys
      use domdec
      use spectra
      implicit none
      character*20  :: keyword
      character*240 :: record
      character*240 :: string
      character*240 :: xyzfile
      integer :: next,i, ios

      qtb_verbose=.FALSE.

      A_gamma_default=0.5
      adqtb_start_classical=.false.
      skipseg=3
      register_spectra=.false.
      save_gamma_history=.FALSE.
      
      piqtb_classical_centroid=.FALSE.
      brieuc_piqtb=.FALSE.

      noQTB=.FALSE.
      qtb_batch_size=-1

      corr_pot=.true.
      corr_pot_corr_fact=.TRUE.
      corr_fact_qtb=1.0d0
      use_corr_fact_qtb=.FALSE.

      adqtb_optimizer="SIMPLE"
      adqtb_tau_avg_default=-1.
      adqtb_tau_adapt=-1.
      adqtb_correct_bias=.FALSE.
      adqtb_win_corr=.FALSE.
      adqtb_avg_adapt=1
      adqtb_smooth_default=-1

!      nstep_therm=0
      do i = 1, nkey
        ios = 0
        next = 1
        record = keyline(i)
        call gettext (record,keyword,next)
        call upcase (keyword)
        string = record(next:240)
        
        select case (trim(keyword))
        case ('A_GAMMA')
          read (string,*,iostat=ios) a_gamma_default
        case ('SKIPSEG')
          read (string,*,iostat=ios) skipseg
        case ('NOQTB')
          noQTB = .true.
        case ('GAMMA_HISTORY')
          save_gamma_history = .true.
        case ('QTB_VERBOSE')
          QTB_verbose=.TRUE.
        case ('CORR_FACT_QTB')
          read (string,*,iostat=ios) corr_fact_qtb(1)
          corr_fact_qtb(:)=corr_fact_qtb(1)
          use_corr_fact_qtb=.TRUE.
        case ('CORR_FACT_NO_POT')
          corr_pot_corr_fact=.FALSE.
        case ('NO_CORR_POT')
          corr_pot = .false.
        case ('REGISTER_SPECTRA')
          register_spectra = .true.
        case ('QTB_BATCH_SIZE')
          read (string,*,iostat=ios) qtb_batch_size
        case ('ADQTB_OPTIMIZER')
          call getword (record,adqtb_optimizer,next)
          call upcase (adqtb_optimizer)
        case ('ADQTB_START_CL')
          adqtb_start_classical = .true.
        case ('ADQTB_TAU_AVG')
          read (string,*,iostat=ios) adqtb_tau_avg_default
        case ('ADQTB_TAU_ADAPT')
          read (string,*,iostat=ios) adqtb_tau_adapt
        case ('ADQTB_BIAS_CORR')
          adqtb_correct_bias=.TRUE.
        case ('ADQTB_AVG_ADAPT')
          read (string,*,iostat=ios) adqtb_avg_adapt
        case ('ADQTB_WIN_CORR')
          adqtb_win_corr=.TRUE.
        case ('ADQTB_SMOOTH')
          read (string,*,iostat=ios) adqtb_smooth_default
        case ('PIQTB_CL_CENTROID')
          piqtb_classical_centroid=.TRUE.
        case ('PIQTB_BRIEUC')
          brieuc_piqtb=.TRUE.
        end select

        if (ios /= 0) then
          write(*,*) "Warning: keyword ",trim(keyword)
     &         ," not correctly read!"
        endif
      end do

      end subroutine read_qtb_keys

      subroutine initialize_adaptation_parameters()
      use keys
      use units
      implicit none
      integer :: i,k
      character*20 keyword
      character*240 dynfile
      character*240 record
      character*240 string
      character*11 opttype
      real(r_p) :: tmp,smooth
      integer next,ios

      allocate(A_gamma(maxtyp))
      allocate(adqtb_tau_avg(maxtyp))
      allocate(adqtb_smooth_type(maxtyp))
      allocate(adqtb_optimizer_type(maxtyp))
      adqtb_smooth_type(:)=adqtb_smooth_default
      A_gamma(:)=A_gamma_default
      adqtb_tau_avg(:)=adqtb_tau_avg_default
      adqtb_optimizer_type(:)=adqtb_optimizer

      do i = 1, nkey
        next = 1
        record = keyline(i)
        call gettext (record,keyword,next)
        call upcase (keyword)
        string = record(next:240)

        if (trim(keyword) /= 'ADQTB_MODIFIER') cycle

        smooth=-1
        tmp=-1
        read(string,*,iostat=ios) k,opttype,tmp,smooth

        if(ios/=0) then
          write(*,*) "Error: ADQTB_MODIFIER keyword is not",
     &         " properly defined"
          call fatal
        endif

        if (smooth<0.d0) smooth=adqtb_smooth_default

        call upcase(opttype)
        if (opttype=="RATIO") then
          if (tmp<0.d0) tmp=adqtb_tau_avg_default
          adqtb_tau_avg(k)=tmp
        elseif (opttype=="SIMPLE") then
          if (tmp<0.d0) tmp=A_gamma_default
          A_gamma(k)=tmp
        else
          write(*,*) "Error: unknown adqtb optimizer ",opttype
          call fatal
        endif
        adqtb_optimizer_type(k)=opttype

        write(*,'(A,I,3x,A,2F10.3)') "ADQTB_MODIFIER"
     &         ,k,opttype,tmp,smooth
      enddo

      adqtb_smooth_type(:)=adqtb_smooth_type(:)/cm1

      end subroutine initialize_adaptation_parameters

      subroutine allocate_qtbdata(iglob)
        use domdec
        use spectra
        use beads
        implicit none
        integer, intent(in) :: iglob
        integer :: i,j,k

        allocate(qtbdata(iglob)%noise(3,2*nseg,nbeads))
        allocate(qtbdata(iglob)%fad(3,nseg,nbeads))
        if(register_spectra) then
          allocate(qtbdata(iglob)%vad(3,nseg,nbeads))
!$acc enter data create(qtbdata(iglob)%noise,qtbdata(iglob)%fad
!$acc&                 ,qtbdata(iglob)%vad) async

!$acc parallel loop collapse(3) async
!$acc&  present(qtbdata(iglob)%vad)
          do k=1,nbeads; do i=1,nseg; do j=1,3
            qtbdata(iglob)%vad(j,i,k)=0.
          enddo; enddo; enddo
        else
!$acc enter data create(qtbdata(iglob)%noise,qtbdata(iglob)%fad) async
        endif
c        allocate(qtbdata(iglob)%proc_has_data(nproc))
c        qtbdata(iglob)%proc_has_data(:)=.FALSE.
c        qtbdata(iglob)%proc_has_data(rank+1)=.TRUE.
        qtbdata_allocated(iglob)=.TRUE.


      end subroutine allocate_qtbdata

      subroutine deallocate_qtbdata(iglob)
        implicit none
        integer, intent(in) :: iglob

        if(register_spectra) then
!$acc exit data delete(qtbdata(iglob)%noise,qtbdata(iglob)%fad
!$acc&                 ,qtbdata(iglob)%vad) async
          deallocate(qtbdata(iglob)%vad)
        else
!$acc exit data delete(qtbdata(iglob)%noise,qtbdata(iglob)%fad) async
        endif
        deallocate(qtbdata(iglob)%noise)
        deallocate(qtbdata(iglob)%fad)
c        deallocate(qtbdata(iglob)%proc_has_data)
        qtbdata_allocated(iglob)=.FALSE.

      end subroutine deallocate_qtbdata

c
c     subroutine adHnoise: compute the kernel for generating the colored noise (fill Htilde)
c
      subroutine adHnoise(write_kernel)
      use spectra
      use atmtyp
      use atomsMirror
      use bath
      use domdec
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use beads, only: nbeads
      use mpi
      implicit none
      integer :: i,j,k,l,iglob, ibead,uu
      logical, intent(in) :: write_kernel
      real(r_p) :: C_omega   !sinus cardinal
      real(r_p) :: f_omega   !smooth cutoff (Fermi function)
      real(r_p) :: theta_tilde
      real(r_p) :: h,u
      real(r_p) :: omega
      real(r_p) :: t
      real(r_p) :: g_ratio,gammak

      if(noQTB) then
        !classical kernel (equipartition of energy)
        Htilde(:,:,:)=sqrt(boltzmann*kelvin)
      else
        Htilde(:,:,:)=0.
        !quantum kernel (Energy propto Theta(omega,kT))
        do ibead=1,nbeads; do i=1,ntype
          !zero frequency term
          gammak=gamma_friction(ibead)
          Htilde(1,i,ibead)=sqrt(boltzmann*kelvin*
     &                    (1.0d0+corr_pot_ratio(1,ibead)))
     &                    *gamma_type(1,i,ibead)/gammak
          do k=2,(3*nseg)/2
            omega=(k-1)*domega   
            theta_tilde=QTB_kernel(omega,piqtb,ibead)*boltzmann*kelvin
            C_omega=(1-2*exp(-gammak*dt_spec)*cos(omega*dt_spec)
     &               +exp(-2*gammak*dt_spec)
     &            )/((gammak**2+omega**2)*dt_spec**2)
            f_omega=1.0d0/(1.0d0+exp((omega-omegacut)/omegasmear))
            g_ratio=1.d0
            if(k.le.nad) then
                g_ratio=gamma_type(k,i,ibead)/gammak
     &            *(1.0d0+corr_pot_ratio(k,ibead))
            endif
            Htilde(k,i,ibead)=sqrt(abs(
     &         theta_tilde*f_omega*C_omega*g_ratio))
            Htilde(3*nseg-k+2,i,ibead)=Htilde(k,i,ibead)  !symmetrize kernel
          enddo
        enddo; enddo
      endif

!$acc update device(Htilde) async

      if (write_kernel) then
        if(piqtb) then
          open(newunit=uu,file="PIQTB_kernel.out")
          do k=2,nad
            omega=(k-1)*domega   
            f_omega=1.0d0/(1.0d0+exp((omega-omegacut)/omegasmear))
            write(uu,'(f17.4)',advance="no") omega*cm1
            do ibead=1,nbeads
              theta_tilde=QTB_kernel(omega,piqtb,ibead)
              write(uu,'(f17.4)',advance="no") theta_tilde*f_omega
            enddo
            write(uu,*)
          enddo
          close(uu)
        endif
      endif

      end subroutine adHnoise

c
c     subroutine convolseg: generate the colored noise by performing a convolution
c                             between the white noise and the kernel Htilde
c
      subroutine convolseg()
      use atmtyp
      use atoms
      use bath
      use beads
      use domdec
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use mpi
      use inform
      use random_mod
#ifdef _OPENACC
      use openacc
      use cufft
      use utilgpu, only: rec_queue
#else
      use glassman
#endif      
      use spectra
      implicit none
      integer :: i,j,k,l,ierr,ibatch,ibeg,iend,nloc_,ibegbeg,iendend
      integer :: iglob,temp_size,nbatch,batch_size_last
      complex(4), allocatable :: s_fft(:,:)

c#ifdef _OPENACC
c      if(nproc>1) then
c!$acc update host(glob) async
c      endif
c#endif

      if (deb_Path) write(*,*) '   >> convolseg'

      if(qtb_separate_fft) then
        nloc_=nloc
        if (piqtb) nloc_=nlocpi
        nbatch=nloc_/(qtb_batch_size)
        batch_size_last=MOD(nloc_,qtb_batch_size)
        temp_size=qtb_batch_size
      else
        temp_size=nloc
        if (piqtb) temp_size=nlocpi
      endif
      if(ranktot==0 .and. QTB_verbose) 
     &     write(*,*) "--- Generating new colored noise",
     &       " segment ---"

      if(piqtb) then
        ibegbeg=ilocpi_beg(rank_polymer+1)
        iendend=ilocpi_end(rank_polymer+1)
      else
        ibegbeg=1
        iendend=nloc
      endif

      call update_fft_plan_qtb(3*nseg)  
      allocate(s_fft(3*nseg,3*temp_size*nbeads))
!$acc enter data create(s_fft) async
      if(qtb_separate_fft) then
        do ibatch=1,nbatch
          ibeg=(ibatch-1)*qtb_batch_size+ibegbeg
          iend=ibatch*qtb_batch_size
          call compute_convolution_batch(ibeg,iend,s_fft) 
        enddo

        if(batch_size_last>0) then
          ibeg=nbatch*qtb_batch_size+ibegbeg
          iend=iendend 
          call compute_convolution_batch(ibeg,iend,s_fft)
        endif
      else
        call compute_convolution_batch(ibegbeg,iendend,s_fft) 
      endif
!$acc exit data delete(s_fft) async
      deallocate(s_fft)

      if (piqtb) call broadcast_piqtb_noise()

      !UPDATE NOISE REPARTITION
c      if(nproc>1) then
c!$acc wait
c        repartnoise(:) = 0 
c        do i=1,nloc 
c          repartnoise(glob(i)) = rank
c        enddo        
c        !call MPI_BARRIER(COMM_TINKER,ierr)
cc
cc     communicate rank that generated qtb noise
cc 
c        call MPI_ALLREDUCE(MPI_IN_PLACE,repartnoise,n,MPI_INT,MPI_SUM,
c     $   COMM_TINKER,ierr) 
c
c      endif
      if (deb_Path) write(*,*) '   << convolseg'

      end subroutine convolseg

      subroutine compute_convolution_batch(ibeg,iend,s_fft)  
      use atmtyp
      use atoms
      use bath
      use beads
      use domdec
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use mpi
      use random_mod
      use inform
      use ascii, only: int_to_str
#ifdef _OPENACC
      use cufft
      use utilgpu  ,only: rec_queue
#else
      use glassman
#endif      
      use spectra
      implicit none
      integer :: i,j,k,l,ierr,ii
      integer :: iglob,ibead,nbatch_
      complex(4), intent(inout) :: s_fft(:,:)
      integer, intent(in) :: ibeg,iend
      real(4) :: rloc,rprev
      real(t_p), allocatable :: r(:)
      complex(t_p), allocatable :: s_fft1d(:),workfft(:)

      if (deb_Path) then
!$acc wait 
        write(*,*) '     >> compute_convolution_batch('
     &    //int_to_str(ibeg)//' => '//int_to_str(iend),
     &    ','//int_to_str(rank)//','//int_to_str(rank_polymer)//')'  
      endif


      nbatch_=(iend-ibeg+1)
      !GENERATE FRESH NOISE FOR THE NEXT SEGMENT
      allocate(r(3*nseg*nbatch_*nbeads))
!$acc enter data create(r) async
#ifdef _OPENACC
      if(.NOT. host_rand_platform) then
        call normalgpu(r,3*nseg*nbatch_*nbeads
     &   ,0.,1.,rec_queue)
      endif
#endif
      if(host_rand_platform) then
        call normalvec(r,3*nseg*nbatch_*nbeads)
!$acc update device(r) async
      endif

      !COPY PREVIOUS SEGMENTS IN WORKING ARRAY AND SHIFT
!$acc parallel loop collapse(4) default(present) async
      do ibead=1,nbeads; do i=ibeg,iend; do k=1,nseg;  do j=1,3
        iglob=glob(i)
        rloc=r(j+3*(k-1)+3*nseg*(i-ibeg)+3*nseg*nbatch_*(ibead-1))
        rprev=qtbdata(iglob)%noise(j,k+nseg,ibead)
        ii=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
        s_fft(k,ii)=cmplx(qtbdata(iglob)%noise(j,k,ibead),0.)
        s_fft(k+nseg,ii) = cmplx(rprev,0.)
        s_fft(k+2*nseg,ii) = cmplx(rloc,0.)
        qtbdata(iglob)%noise(j,k,ibead)=rprev
        qtbdata(iglob)%noise(j,k+nseg,ibead)=rloc
      enddo; enddo; enddo; enddo
!$acc exit data delete(r) async

      !FOURIER TRANSFORM => frequency space
#ifdef _OPENACC
!$acc host_data use_device(s_fft)
      ierr = cufftExecC2C(plan_qtb,s_fft,s_fft,CUFFT_FORWARD)
!$acc end host_data
      if(ierr/=0) then
        write(0,*) "first FFT failed in convolseg with error",ierr
        call fatal
      endif
#else
      allocate(s_fft1d(3*nseg))
      allocate(workfft(3*nseg))
      do ibead=1,nbeads; do i=1,nloc; do j=1,3
        ii=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
        s_fft1d(:)=s_fft(:,ii)
        call SPCFFT(s_fft1d,3*nseg,-1,workfft)
        s_fft(:,ii)=s_fft1d(:)
      enddo; enddo; enddo
#endif

      !MULTIPLY THE NOISE BY THE KERNEL IN FOURIER SPACE
!$acc parallel loop collapse(4) default(present) async
      do ibead=1,nbeads; do i=ibeg,iend ; do j=1,3; do k=1,3*nseg
        iglob=glob(i)
        ii=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
        s_fft(k,ii)=s_fft(k,ii)*Htilde(k,adqtb_type(iglob),ibead)
      enddo; enddo; enddo ; enddo

      !FOURIER TRANSFORM => back to time
#ifdef _OPENACC
!$acc host_data use_device(s_fft)
      ierr = cufftExecC2C(plan_qtb,s_fft,s_fft,CUFFT_INVERSE)
!$acc end host_data
      if(ierr/=0) then
        write(0,*) "second FFT failed in convolseg with error",ierr
        call fatal
      endif
#else
      do ibead=1,nbeads;do i=ibeg,iend ; do j=1,3
        ii=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
        s_fft1d(:)=s_fft(:,ii)
        call SPCFFT(s_fft1d,3*nseg,1,workfft)
        s_fft(:,ii)=s_fft1d(:)
      enddo; enddo; enddo
      deallocate(s_fft1d)
      deallocate(workfft)
#endif
      !USE THE MIDDLE SEGMENT FOR THE RANDOM FORCE
!$acc parallel loop collapse(4) default(present) async
      do ibead=1,nbeads; do i=ibeg,iend; do j=1,3; do k=1,nseg
        iglob=glob(i)
        ii=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
        qtbdata(iglob)%fad(j,k,ibead)=
     &     real(s_fft(nseg+k,ii),t_p)/(3*nseg)
     &        *sqrt(2.d0*mass(iglob)*gamma_friction(ibead)/dt_spec)
      enddo; enddo; enddo; enddo


      if (deb_Path) then
!$acc wait
        write(*,*) '     << compute_convolution_batch('
     &    //int_to_str(ibeg)//' => '//int_to_str(iend),
     &    ','//int_to_str(ranktot)//')'  
      endif

      end subroutine compute_convolution_batch


c
c     subroutine qtbrandom: propagate O block (from BAOAB) using the QTB random force
c                            
      subroutine qtbrandom()
      use atmtyp
      use atoms
      use bath
      use domdec
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use mpi
      use spectra
#ifdef _OPENACC
      use cufft
#endif
      implicit none
      integer i,j,k,l,ierr
      integer iglob,inoise
      real(r_p) a1

      if(piqtb .and. ranktot==0) then
        write(*,*) "Warning: qtbrandom not implemented for piqtb"
        call fatal
      endif

      inoise=compteur%i
      a1 = exp(-gamma*dt_spec)
c
c     save velocity and random force if you need spectra 
c     apply QTB random force
c      
      if (register_spectra) then

!$acc parallel loop  collapse(2) async default(present)
        do i=1,nloc ; do j=1,3
          iglob = glob(i)
          qtbdata(iglob)%vad(j,inoise,1)=real(sqrt(a1)*v(j,iglob)
     $        +0.5*dt_spec*qtbdata(iglob)%fad(j,inoise,1)/mass(iglob),4)
          v(j,iglob)=a1*v(j,iglob) 
     &     +dt_spec*qtbdata(iglob)%fad(j,inoise,1)/mass(iglob)
        enddo ;enddo

      else

!$acc parallel loop collapse(2) async default(present)
        do i=1,nloc; do j=1,3        
            iglob=glob(i)
            v(j,iglob)=a1*v(j,iglob) 
     &       +dt_spec*qtbdata(iglob)%fad(j,inoise,1)/mass(iglob)  
        enddo; enddo
        
      end if

      call QTB_update_state()

      end subroutine qtbrandom

      function get_qtb_compteur() result(inoise)
      implicit none
      integer :: inoise
      inoise=compteur%i
      end function get_qtb_compteur

!      subroutine next_qtb_kick(kick,inoise)
!      use atmtyp
!      use atoms
!      use bath
!      use beads
!      use domdec
!      use langevin
!      use math
!      use mdstuf
!      use moldyn
!      use usage
!      use units
!      use mpi
!      use spectra
!#ifdef _OPENACC
!      use cufft
!#endif
!      implicit none
!      integer, intent(out) :: inoise
!      real(t_p), intent(inout) :: kick(*)
!      
!      integer i,j,k,l,ierr,ii,ibeg,iend,nloc_
!      integer iglob, ibead
!      real(r_p) a1
!
!      inoise=compteur%i
!
!      if(piqtb) then
!        ibeg=ilocpi_beg(rank_polymer+1)
!        iend=ilocpi_end(rank_polymer+1)
!        nloc_=nlocpi
!      else
!        ibeg=1
!        iend=nloc
!        nloc_=nloc
!      endif
!c
!c     apply QTB random force
!c
!!$acc parallel loop collapse(3) async default(present)
!      do ibead=1,nbeads; do i=ibeg,iend; do j=1,3        
!          iglob=glob(i)
!          ii=j+3*(i-ibeg)+3*nloc_*(ibead-1)
!          kick(ii)=dt_spec*qtbdata(iglob)%fad(j,inoise,ibead)
!     &                         /mass(iglob)  
!      enddo; enddo ; enddo
!
!      if (.not.register_spectra) inoise=-1
!
!      end subroutine next_qtb_kick

      subroutine QTB_update_state()
      use atmtyp
      use atoms
      use bath
      use domdec
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use mpi
      use spectra
      implicit none
c
c     At the end of each segment compute the the new gamma and generate
c     new noise
c
      call compteur%update()
      if (.not.compteur%reset_step) return
      
      if(register_spectra .AND. (compteur%i_avg >= skipseg)) then
        if(isobaric_save .AND. (.NOT. isobaric)) then
          isobaric=.TRUE.
          use_piston=.TRUE.
          if(ranktot==0  .and. QTB_verbose) 
     &         write(*,*) "--- Activating QTB Piston ---"
        endif
        if(ranktot==0  .and. QTB_verbose) 
     &     write(*,*) "--- Writing QTB spectra ---"
        call QTB_spectra_adaptation_gamma()
      endif

      call cleanup_qtbdata()
      call convolseg()

      end subroutine QTB_update_state


c
c     subroutine QTB_spectra_adaptation_gamma: compute Cvv/Cvf spectra and DeltaFDT
c                                              and adapt gammar if adaptive
      subroutine QTB_spectra_adaptation_gamma()
      use ascii, only: int_to_str
      use atmtyp
      use atoms
      use bath
      use beads
      use domdec
      use deconvolution
      use energi
      use freeze
      use katoms, only: weight
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use mpi
      use spectra
      use glassman
      use inform, only: deb_Path
      implicit none
      integer :: i,j,k,l,ierr,ierr1,ierr2,idof,nloc_
      integer :: iglob, itype,nbatch,batch_size_last
      integer :: ibegbeg,iendend, ibead
      integer :: ibatch,iend,ibeg,temp_size,uu
      complex(4),  allocatable :: sv(:,:),sf(:,:)
      real(r_p),  allocatable :: mCvv_type(:,:,:,:), Cvf_type(:,:,:)
      real(r_p), allocatable :: dFDR_type(:,:,:)
      real(r_p), allocatable :: mCvv_sum(:,:,:)
      real(r_p) :: a1,gamma_disp,norm,gamma_ref
      real(r_p) :: b_avg, avg_correct,m,temp,mcvvsum_tmp
      logical :: qtb_therm      

      if (deb_Path) write(*,*) '   >> QTB_spectra_adaptation_gamma'

      call update_fft_plan_qtb(3*nseg)
c
c    reset averages if still in thermalization 
c
      qtb_therm=compteur%i_avg .le. compteur%startsave
      if (qtb_therm) then
        mCvv_average_type(:,:,:)=0d0
        mCvv_average_sum(:,:,:)=0d0
        Cvf_average_type(:,:,:)=0d0
        dFDR_average_type(:,:,:)=0d0
      endif

c     Allocation of the different used for the adqtb
      allocate(mCvv_type(3*nseg,ntype,3,nbeads))
      allocate(Cvf_type(3*nseg,ntype,nbeads))
      allocate(dFDR_type(nad,ntype,nbeads))
      allocate(mCvv_sum(nad,3,nbeads))
      
      if (piqtb) call reduce_piqtb_vad()
      !call gather_vad()

      if(qtb_separate_fft) then
        nloc_=nloc
        if (piqtb) nloc_=nlocpi
        nbatch=nloc_/qtb_batch_size
        batch_size_last=MOD(nloc_,qtb_batch_size)
        temp_size=qtb_batch_size
      else
        temp_size=nloc
        if (piqtb) temp_size=nlocpi
      endif

      if(piqtb) then
        ibegbeg=ilocpi_beg(rank_polymer+1)
        iendend=ilocpi_end(rank_polymer+1)
      else
        ibegbeg=1
        iendend=nloc
      endif

      allocate(sv(3*nseg,3*temp_size*nbeads))
      allocate(sf(3*nseg,3*temp_size*nbeads))
!$acc data create(sv,sf) copyout(mCvv_type,Cvf_type) async

!$acc parallel loop collapse(3) async
      do ibead=1,nbeads; do i=1,ntype; do k=1,3*nseg
!$acc loop seq
        do j=1,3
          mCvv_type(k,i,j,ibead)=0.0_re_p
        enddo
        Cvf_type(k,i,ibead)=0.0_re_p
      enddo; enddo; enddo

      if(qtb_separate_fft) then
        do ibatch=1,nbatch
          ibeg=(ibatch-1)*qtb_batch_size+ibegbeg
          iend=ibatch*qtb_batch_size
          call compute_QTB_spectra_batch(ibeg,iend,3*nseg
     &       ,sv,sf,mCvv_type,Cvf_type)
        enddo
        if(batch_size_last>0) then
          ibeg=nbatch*qtb_batch_size+ibegbeg
          iend=iendend 
          call compute_QTB_spectra_batch(ibeg,iend,3*nseg
     &       ,sv,sf,mCvv_type,Cvf_type)
        endif
      else
        call compute_QTB_spectra_batch(ibegbeg,iendend,3*nseg
     &     ,sv,sf,mCvv_type,Cvf_type)
      endif
!$acc end data
!$acc wait
      deallocate(sv,sf)
c
c     the master gets the sums of Cvv and Cvf contributions
c
c      if(piqtb .and. nproc_polymer>1) then
c        ! allreduce on the polymer communicator 
c        ! so that only spatial comm is needed now
c        call reduce_qtb_spectrum(mCvv_type,nbeads*ntype*3*nseg*3
c     &                            ,COMM_POLYMER,rank_polymer,.TRUE.)
c        call reduce_qtb_spectrum(Cvf_type,nbeads*ntype*3*nseg
c     &                            ,COMM_POLYMER,rank_polymer,.TRUE.)
c      endif

      if(nproctot>1) then 
        call reduce_qtb_spectrum(mCvv_type,nbeads*ntype*3*nseg*3
     &                              ,MPI_COMM_WORLD,ranktot,.FALSE.)
        call reduce_qtb_spectrum(Cvf_type,nbeads*ntype*3*nseg
     &                             ,MPI_COMM_WORLD,ranktot,.FALSE.)
      endif

      if (ranktot.eq.0) then
        do ibead=1,nbeads; do i=1,ntype
          m=weight(adqtb_type_map(i))

          mCvv_type(:,i,:,ibead)=mCvv_type(:,i,:,ibead)/n_of_type(i)
          Cvf_type(:,i,ibead)=Cvf_type(:,i,ibead)/n_of_type(i)

          if(adqtb_win_corr) then
            do j=1,3
              call correct_window_spectrum(mCvv_type(:,i,j,ibead))
            enddo
            call correct_window_spectrum(Cvf_type(:,i,ibead))
          endif

          do k=1,nad
          mcvvsum_tmp=SUM(mCvv_type(k,i,:,ibead))
          !compute Delta FDR
          dFDR_type(k,i,ibead)=mcvvsum_tmp*gamma_type(k,i,ibead)
     &       -Cvf_type(k,i,ibead) !/gamma_type(k,i)

          ! update intermediate averages for adaptation
          temp=mcvvsum_tmp - mCvv_mean(k,i,ibead)
          mCvv_mean(k,i,ibead)=mCvv_mean(k,i,ibead)+temp/compteur_bis%i
          temp=Cvf_type(k,i,ibead) - Cvf_mean(k,i,ibead)
          Cvf_mean(k,i,ibead)=Cvf_mean(k,i,ibead)+temp/compteur_bis%i
          temp=dFDR_type(k,i,ibead) - dFDR_mean(k,i,ibead)
          dFDR_mean(k,i,ibead)=dFDR_mean(k,i,ibead)+temp/compteur_bis%i
          enddo          
        enddo; enddo
c
c      adapt gamma
c
        if(adaptive_qtb) then         
          call compteur_bis%update()  
          if (compteur_bis%reset_step) then 
             if(ranktot==0 .and. QTB_verbose) then
              write(*,*) "--- Adapting gammar for adQTB ---"
            endif  
            do i=1,ntype
              call adapt_gamma(i)%p(i,mCvv_mean(:,i,:),Cvf_mean(:,i,:)
     &             ,dFDR_mean(:,i,:),qtb_therm)
            enddo
          endif
        endif

        if(.not.adqtb_win_corr) then
          ! apply window correction if not already done
          ! (just for printing, not for adaptation)
          do ibead=1,nbeads; do i=1,ntype
            do j=1,3
              call correct_window_spectrum(mCvv_type(:,i,j,ibead))
            enddo
            call correct_window_spectrum(Cvf_type(:,i,ibead))
            do k=1,nad
              !recompute Delta FDR after correction
              dFDR_type(k,i,ibead)=SUM(mCvv_type(k,i,:,ibead))
     &              *gamma_type(k,i,ibead) - Cvf_type(k,i,ibead)
            enddo
          enddo; enddo
        endif
        
        ! update averages
        mCvv_sum(:,:,:)=0.d0
        do ibead=1,nbeads; do i=1,ntype; do k=1,nad
          mCvv_sum(k,:,ibead)=mCvv_sum(k,:,ibead)
     &            +mCvv_type(k,i,:,ibead)*n_of_type(i)
          mCvv_average_type(k,i,ibead)=mCvv_average_type(k,i,ibead) 
     &                       +SUM(mCvv_type(k,i,:,ibead))
          Cvf_average_type(k,i,ibead)=Cvf_average_type(k,i,ibead) 
     &            +Cvf_type(k,i,ibead)/gamma_type(k,i,ibead) 
          dFDR_average_type(k,i,ibead)=dFdr_average_type(k,i,ibead) 
     &              +dFdr_type(k,i,ibead)
        enddo; enddo; enddo
        mCvv_sum=mCvv_sum/real(n,r_p)
        mCvv_average_sum=mCvv_average_sum+mCvv_sum
        mCvv_sum=mCvv_average_sum/compteur%nsample

c       WRITE SPECTRA AND GAMMAR
        do ibead=1,nbeads; do i=1,ntype
          m=weight(adqtb_type_map(i))
          if(piqtb) then
            open(newunit=uu, file="QTB_spectra_"//trim(adqtb_name(i))
     &           //"_"//int_to_str(adqtb_type_map(i))
     &           //'_mode'//int_to_str(ibead,3)
     &           //".out")
          else
            open(newunit=uu, file="QTB_spectra_"//trim(adqtb_name(i))
     &           //"_"//int_to_str(adqtb_type_map(i))
     &           //".out")
          endif
          write(uu,'(A)') 
     &        "#Frequency   mCvv   Cvf/gammar  dFDT  gammar"
          do l=1,nad
            write(uu,'(5E16.8)') domega*(l-1)*cm1
     &         ,mCvv_average_type(l,i,ibead)/compteur%nsample
     &         ,Cvf_average_type(l,i,ibead)/compteur%nsample
     &         ,dFDR_average_type(l,i,ibead)/compteur%nsample   
     &         ,gamma_type(l,i,ibead)
          enddo
          close(uu)
        enddo; enddo

        if(adaptive_qtb) then
          open(newunit=uu, file='gamma_restart.out')  
          call write_frequency_type_header(uu)
          do l=1,nad
            write(uu,'('//int_to_str(1+ntype*nbeads)//'E16.8)') 
     &              domega*(l-1)*cm1,gamma_type(l,:,:)
          enddo
          close(uu) 
          
          if(save_gamma_history) then
            open(newunit=uu, file='gamma_history.out',access='append')
            call write_frequency_type_header(uu)
            do l=1,nad
              write(uu,'('//int_to_str(1+ntype*nbeads)//'E16.8)') 
     &              domega*(l-1)*cm1,gamma_type(l,:,:)
            enddo
            close(uu)
          endif
              
        endif
        
      end if

      if(adaptive_qtb) then
c     the master sends the adapted gamma_r to everybody
        call MPI_BCAST(gamma_type,ntype*nad*nbeads,MPI_RPREC
     &        ,0,MPI_COMM_WORLD,ierr)
c     recompute the noise kernel with the new gamma
        call adHnoise(.FALSE.)
      endif

c
c       COMPUTE corr_fact_qtb from deconvolution
c
      if(use_corr_fact_qtb .OR. noQTB 
     & .OR. corr_fact_converged
     & .OR. (adqtb_start_classical .AND. (.not. adaptive_qtb)) 
     &  ) return
     
      ! TODO: corr_fact for piqtb
      call compute_corr_fact(mCvv_sum)

      
      if (deb_Path) write(*,*) '   << QTB_spectra_adaptation_gamma'

      
      end subroutine QTB_spectra_adaptation_gamma

      subroutine reduce_qtb_spectrum(spectrum,size,comm,rank,allreduce)
      use mpi
      implicit none
      integer, intent(in) :: size,comm,rank
      logical, intent(in) :: allreduce
      real(r_p), intent(inout) :: spectrum(*)
      integer :: ierr

      if (allreduce) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,spectrum,
     &       size,MPI_RPREC,MPI_SUM,
     &       comm,ierr)
      else
        if (rank.eq.0) then
          call MPI_REDUCE(MPI_IN_PLACE,spectrum,
     &       size,MPI_RPREC,MPI_SUM,0,
     &       comm,ierr)
        else
          call MPI_REDUCE(spectrum,spectrum,
     &       size,MPI_RPREC,MPI_SUM,0,
     &       comm,ierr)
        end if
      endif

      end subroutine reduce_qtb_spectrum

      subroutine compute_QTB_spectra_batch(ibeg,iend,nfft,sv,sf
     &   ,mCvv_type, Cvf_type)
      use atmtyp
      use atoms
      use bath
      use beads
      use domdec
      use energi
      use freeze
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use mpi
      use spectra
#ifdef _OPENACC
      use cufft
#else
      use glassman
#endif
      implicit none
      integer :: i,j,k,l,ierr,ierr1,ierr2,idof
      integer :: iglob, itype, ibead,nbatch_
      integer, intent(in) :: ibeg,iend,nfft
      complex(4),  intent(inout) :: sv(:,:),sf(:,:)
      real(r_p), intent(inout) :: mCvv_type(:,:,:,:), Cvf_type(:,:,:)
#ifndef _OPENACC
      complex(t_p), allocatable :: s_fft1d(:),workfft(:)
#endif

      nbatch_=(iend-ibeg+1)

      !PUT VELOCITIES AND FORCE IN TEMPORARY COMPLEX ARRAY
!$acc parallel loop collapse(4) async default(present)
      do ibead=1,nbeads; do i=ibeg,iend; do j=1,3; do k=1,nfft  
          iglob=glob(i)
          idof=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
          if(k<=nseg) then     
            sv(k,idof)=cmplx(qtbdata(iglob)%vad(j,k,ibead),0.)
            sf(k,idof)=cmplx(qtbdata(iglob)%fad(j,k,ibead),0.)
          else
            sv(k,idof)=cmplx(0.,0.)
            sf(k,idof)=cmplx(0.,0.)
          endif        
      enddo; enddo; enddo; enddo

      !FOURIER TRANSFORM
#ifdef _OPENACC
!$acc host_data use_device(sv,sf)
      ierr1 = cufftExecC2C(plan_qtb,sv,sv,CUFFT_FORWARD)
      ierr2 = cufftExecC2C(plan_qtb,sf,sf,CUFFT_FORWARD)
!$acc end host_data
      if(ierr1/=0 .or. ierr2/=0) then
        write(0,*) "FFT failed in QTB_spectra_adaptation_gamma"
        call fatal
      endif
#else
      allocate(s_fft1d(nfft))
      allocate(workfft(nfft))
      do ibead=1,nbeads; do i=ibeg,iend; do j=1,3
        idof=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
        s_fft1d(:)=sv(:,idof)
        call SPCFFT(s_fft1d,nfft,-1,workfft)
        sv(:,idof)=s_fft1d(:)
        s_fft1d(:)=sf(:,idof)
        call SPCFFT(s_fft1d,nfft,-1,workfft)
        sf(:,idof)=s_fft1d(:)
      enddo; enddo; enddo
      deallocate(s_fft1d)
      deallocate(workfft)
#endif

      !COMPUTE Cvv AND Cvf FROM FOURIER TRANSFORMED TRAJECTORIES
!$acc parallel loop collapse(4) default(present) async
      do ibead=1,nbeads; do i=ibeg,iend; do j=1,3; do k=1,nfft   
          idof=j+3*(i-ibeg)+3*nbatch_*(ibead-1)
          iglob=glob(i)
          itype=adqtb_type(iglob)
!$acc atomic update
          mCvv_type(k,itype,j,ibead)=mCvv_type(k,itype,j,ibead)
     &         +real(abs(sv(k,idof)),r_p)**2
     &            *mass(iglob)*dt_spec/real(nfft,r_p)
!$acc atomic update
          Cvf_type(k,itype,ibead)=Cvf_type(k,itype,ibead)
     &         +real(sv(k,idof)*conjg(sf(k,idof)),r_p)
     &            *dt_spec/real(nfft,r_p)     
         
      enddo; enddo; enddo; enddo   

      end subroutine compute_QTB_spectra_batch
      

      subroutine write_frequency_type_header(unit)
        use ascii, only: int_to_str
        use beads, only: nbeads
        implicit none
        integer, intent(in) :: unit
        integer :: i,ibead

        write(unit,'(A)',advance='no') '#Frequency' 
        do ibead=1,nbeads;do i=1,ntype
          write(unit,'(6x,A)',advance='no') 
     &       "Type_"//trim(adqtb_name(i))
     &        //"_"//int_to_str(adqtb_type_map(i))
        enddo; enddo
        write(unit,*)

      end subroutine write_frequency_type_header

      subroutine update_fft_plan_qtb(nfft) 
        use domdec, only:nloc,rank
#ifdef _OPENACC
        use cufft
        use utilgpu, only: rec_queue
        use openacc
#endif
        use spectra, only: nseg
        use beads
        implicit none
        integer, intent(in) :: nfft
        integer, pointer :: pnull
        integer :: ierr,nloc_plan_new,nloc_

        nloc_=nloc
        if (piqtb) nloc_=nlocpi 

        if(qtb_separate_fft) then
          if(qtb_batch_size>nloc_) then
            qtb_separate_fft=.FALSE.
            nloc_plan_new=nloc_
            if(rank==0) 
     &         write(*,*) "QTB_BATCH_SIZE>nloc: "
     &           ,"switching off QTB_SEPARATE_FFT"
          else
            nloc_plan_new=qtb_batch_size
          endif          
        else
          nloc_plan_new=nloc_
        endif

#ifdef _OPENACC
        if(nloc_plan_new==nloc_plan .AND. nfft_plan==nfft) return
        if(nloc_plan>=0) ierr = cufftDestroy(plan_qtb) 
        pnull => null()
        nloc_plan = nloc_plan_new
        nfft_plan=nfft
        ierr = cufftPlanMany(plan_qtb,1,(/nfft/)
     &          ,pnull,0,0
     &          ,pnull,0,0
     &          ,CUFFT_C2C, 3*nloc_plan*nbeads )
        ierr=ierr+cufftSetStream(plan_qtb
     &     ,acc_get_cuda_stream(rec_queue))  
#endif
      end subroutine update_fft_plan_qtb

c-------------------------------------------------------------
c       COMMUNICATION ROUTINES

      subroutine reassignqtb(nloc_save,max_data_recv
     &      ,nneig_recep, n_data_recv, iglob_recv, pneig_recep
     &      ,nneig_send , n_data_send, iglob_send, pneig_send 
     &      ,max_atoms_recv, max_atoms_send )
        use atomsMirror
        use mpi
        use domdec
        use spectra
        use inform, only: deb_Path
        use beads
        implicit none
        integer, intent(in) :: nloc_save, max_data_recv
        integer,intent(in) :: nneig_recep, n_data_recv(nneig_recep)
     &                     , iglob_recv(2*max_data_recv,nneig_recep)
     &                     , pneig_recep(:)
        integer,intent(in) :: nneig_send, n_data_send(0:nneig_send)
     &                     , iglob_send(1:2*nloc_save,0:nneig_send)
     &                     , pneig_send(:)
        integer, intent(in) :: max_atoms_recv, max_atoms_send
        real(4), allocatable :: data_send(:,:,:,:,:)
     &                        , data_recv(:,:,:,:,:)
        integer :: nseg_comm,i,j,ineighbor,ierr,iseg,iglob
        integer :: req_data_send(nneig_send),req_data_recv(nneig_recep)
        integer :: pn,ii,jj,kk,ibead
        !integer :: n_data_s(nneig_send), n_data_r(nneig_recep)
        !integer :: iglob_s(max_atoms_send,nneig_send)
        !integer :: iglob_r(max_atoms_recv,nneig_recep)
        !integer :: max_atoms_s, max_atoms_r
        integer :: nneig_s,nneig_r
        !logical :: has_data_send(max_atoms_send,nneig_send)
        !logical :: has_data_recv(max_atoms_recv,nneig_recep)
        !integer :: n_ask_send(nneig_send), n_ask_recv(nneig_recep)
        logical :: was_allocated
        
        if (deb_Path) then
!$acc wait
          write(*,*) '   >> reassignqtb'
        endif

        nseg_comm=3*nseg
        if (register_spectra) nseg_comm=4*nseg
        req_data_send (:) = MPI_REQUEST_NULL
        req_data_recv (:) = MPI_REQUEST_NULL

!$acc update host(iglob_send,iglob_recv) async
!$acc wait

        !ASK IF NEIGHBORS HAVE DATA
c        n_ask_send(:)=0
c        ii=0
c        do ineighbor=1,nneig_send
c          pn=pneig_send(ineighbor)
c          do i=1,n_data_send(ineighbor)
c            iglob=iglob_send(2*(i-1)+1,ineighbor)
c            ! IF I THINK THAT pn DOES NOT HAVE DATA, IT MEANS WE NEVER cCOMMUNICATED (THIS ATOMS) DURING SEGMENT
c            ! => I SHOULD ASK THEM IF THEY HAVE DATA OR NOT
c            if(.not.qtbdata(iglob)%proc_has_data(pn+1)) then
c              n_ask_send(ineighbor)=n_ask_send(ineighbor) + 1
c            endif
c          enddo
c          if(n_ask_send(ineighbor)>0) then
c            ii=ii+1
c            call MPI_IRECV(has_data_send(1,ineighbor) 
c     &      ,n_ask_send(ineighbor), MPI_LOGICAL, pn, 0
c     &      , COMM_TINKER, req_data_send(ii),ierr)
c          endif
c        enddo
c
c        n_ask_recv(:)=0 
c        max_atoms_r=0
c        n_data_r(:)=0
c        jj=0
        do ineighbor=1,nneig_recep
          pn=pneig_recep(ineighbor)         
          do i=1,n_data_recv(ineighbor)
            iglob=iglob_recv(2*(i-1)+1,ineighbor)  
            !was_allocated=qtbdata_allocated(iglob)       
            if(.not.qtbdata_allocated(iglob)) then
              !allocate qtbdata for receiving atom for which I don't have data yet
              !n_data_r(ineighbor)=n_data_r(ineighbor)+1
              !iglob_r(n_data_r(ineighbor),ineighbor)=iglob
              call allocate_qtbdata(iglob) 
            endif
            ! IF I THINK THAT pn DOES NOT HAVE DATA, IT MEANS WE NEVER COMMUNICATED (THIS ATOMS) DURING SEGMENT
            ! => I SHOULD TELL THEM IF I HAVE DATA OR NOT
c            if(.not.qtbdata(iglob)%proc_has_data(pn+1)) then
c              n_ask_recv(ineighbor)=n_ask_recv(ineighbor) + 1
c              has_data_recv(n_ask_recv(ineighbor),ineighbor)
c     &           =was_allocated
c              qtbdata(iglob)%proc_has_data(pn+1)=.TRUE.
c            endif            
c            qtbdata(iglob)%currently_loc=.TRUE.
          enddo
c          if(n_data_r(ineighbor)>max_atoms_r)
c     &        max_atoms_r=n_data_r(ineighbor)
c          if(n_ask_recv(ineighbor)>0) then
c            jj=jj+1
c            call MPI_ISEND(has_data_recv(1,ineighbor) 
c     &      ,n_ask_recv(ineighbor), MPI_LOGICAL, pn, 0
c     &      , COMM_TINKER, req_data_recv(jj),ierr)
c          endif
        enddo

c        call MPI_WAITALL(ii,req_data_send 
c     &                  ,MPI_STATUSES_IGNORE, ierr)
c        call MPI_WAITALL(jj,req_data_recv 
c     &                  ,MPI_STATUSES_IGNORE, ierr)


        !FILTER OUT THE PROCS WHO ALREADY HAVE DATA
c        max_atoms_s=0
c        n_data_s(:)=0
c        do ineighbor=1,nneig_send 
c          !ii=0
c          pn=pneig_send(ineighbor)
c          do i=1,n_data_send(ineighbor)
c            iglob=iglob_send(2*(i-1)+1,ineighbor)
c            !if(.not.qtbdata(iglob)%proc_has_data(pn+1)) then
c            !  ii=ii+1
c            !  if(.not. has_data_send(ii,ineighbor)) then
c            !    !send data if neigh does not have data yet
c            !    n_data_s(ineighbor)=n_data_s(ineighbor)+1
c            !    iglob_s(n_data_s(ineighbor),ineighbor)=iglob
c            !  endif              
c            !endif
cc            qtbdata(iglob)%proc_has_data(pn+1)=.TRUE.
c            !mark atom as not local
c            qtbdata(iglob)%currently_loc=.FALSE.
cc            call MPI_Isend(qtbdata(iglob)%proc_has_data
cc     &        ,nproc, MPI_LOGICAL
cc     &        ,pneig_send(ineighbor),1,COMM_TINKER
cc     &        ,req_data_send(ineighbor), ierr)
c
c          enddo
cc          if(n_data_s(ineighbor)>max_atoms_s)
cc     &        max_atoms_s=n_data_s(ineighbor)
c        enddo

        req_data_send (:) = MPI_REQUEST_NULL
        req_data_recv (:) = MPI_REQUEST_NULL

        allocate(data_send(3,nseg_comm,nbeads
     &               ,max_atoms_send,nneig_send))
        allocate(data_recv(3,nseg_comm,nbeads
     &               ,max_atoms_recv,nneig_recep))
!$acc data create(data_send,data_recv) 
!$acc&   copyin(n_data_send,n_data_recv) async

        !put data to send from qtbdata to buffer
!$acc parallel loop async default(present)
!$acc&      vector_length(512)
        do ineighbor=1,nneig_send 
!$acc loop vector collapse(4)
          do i=1,n_data_send(ineighbor); do ibead=1,nbeads
          do iseg=1,nseg; do j=1,3
            iglob=iglob_send(2*(i-1)+1,ineighbor)
            data_send(j,iseg,ibead,i,ineighbor)
     &         =qtbdata(iglob)%fad(j,iseg,ibead)
            data_send(j,iseg+nseg,ibead,i,ineighbor)
     &         =qtbdata(iglob)%noise(j,iseg,ibead)
            data_send(j,iseg+2*nseg,ibead,i,ineighbor)
     &         =qtbdata(iglob)%noise(j,iseg+nseg,ibead)
            if (register_spectra)
     &        data_send(j,iseg+3*nseg,ibead,i,ineighbor)
     &          =qtbdata(iglob)%vad(j,iseg,ibead)
c              data_send(j,iseg,i,ineighbor)
c     &           =qtbdata(iglob)%fad(j,iseg,1)
          enddo; enddo ; enddo; enddo
        enddo

!$acc wait
!$acc host_data use_device(data_recv,data_send)
        nneig_r=0
        do ineighbor = 1, nneig_recep
          if(n_data_recv(ineighbor)>0) then
            nneig_r=nneig_r+1
            call MPI_Irecv(data_recv(1,1,1,1,ineighbor)
     &        ,nseg_comm*3*n_data_recv(ineighbor)*nbeads
     &        ,MPI_REAL4, pneig_recep(ineighbor), 2, COMM_TINKER
     &        ,req_data_recv(nneig_r), ierr)
          endif
        end do

        nneig_s=0
        do ineighbor = 1, nneig_send
          if(n_data_send(ineighbor)>0) then
            nneig_s=nneig_s+1
            call MPI_Isend(data_send(1,1,1,1,ineighbor)
     &        ,nseg_comm*3*n_data_send(ineighbor)*nbeads
     &        ,MPI_REAL4, pneig_send(ineighbor), 2, COMM_TINKER
     &        ,req_data_send(nneig_s), ierr)
          endif
        end do
!$acc end host_data

        ! Wait all communications
        call MPI_Waitall(nneig_r,req_data_recv
     &                  ,MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(nneig_s ,req_data_send 
     &                  ,MPI_STATUSES_IGNORE, ierr)


        !put data received from buffer to qtbdata
!$acc parallel loop async default(present)
!$acc&         vector_length(512)
        do ineighbor=1,nneig_recep
!$acc loop vector collapse(4)
          do i=1,n_data_recv(ineighbor); do ibead=1,nbeads
          do iseg=1,nseg; do j=1,3
            iglob=iglob_recv(2*(i-1)+1,ineighbor)
            qtbdata(iglob)%fad(j,iseg,ibead)
     &         =data_recv(j,iseg,ibead,i,ineighbor)
            qtbdata(iglob)%noise(j,iseg,ibead)
     &         =data_recv(j,iseg+nseg,ibead,i,ineighbor)
            qtbdata(iglob)%noise(j,iseg+nseg,ibead)
     &         =data_recv(j,iseg+2*nseg,ibead,i,ineighbor)
            if (register_spectra)
     &         qtbdata(iglob)%vad(j,iseg,ibead)
     &           =data_recv(j,iseg+3*nseg,ibead,i,ineighbor)
c             qtbdata(iglob)%fad(j,iseg,1)
c     &         =data_recv(j,iseg,i,ineighbor)
          enddo; enddo; enddo ; enddo
        enddo

        do ineighbor=1,nneig_send 
          pn=pneig_send(ineighbor)
          do i=1,n_data_send(ineighbor)
            iglob=iglob_send(2*(i-1)+1,ineighbor)
            call deallocate_qtbdata(iglob)
          enddo
        enddo

!$acc end data

      deallocate(data_send,data_recv)

        if (deb_Path) then
!$acc wait
          write(*,*) '   << reassignqtb'
        endif
      end subroutine reassignqtb

      subroutine cleanup_qtbdata()
      use atoms
      use beads
      use spectra
      use domdec
      use inform, only: deb_Path
      implicit none
      integer :: i,j,k,iglob,ibead

      if(nproc==1) return

      if (deb_Path) write(*,*) '   >> cleanup_qtbdata'

      !call gather_noise()

      if(register_spectra) then
!$acc parallel loop async default(present) collapse(4)
        do ibead=1,nbeads; do i=1,nloc; do k=1,nseg; do j=1,3
          iglob=glob(i)
          qtbdata(iglob)%vad(j,k,ibead)=0.
        enddo; enddo; enddo; enddo
      endif

c      do i=1,n
c        if(qtbdata_allocated(i)) then 
c          if(.not. qtbdata(i)%currently_loc) then
c            qtbdata(i)%proc_has_data(:)=.FALSE.
c            qtbdata(i)%proc_has_data(rank+1)=.TRUE.
c          else
c            call deallocate_qtbdata(i)
c          endif
c        endif
c      enddo
c!$acc update host(glob)
c      do i=1,nloc
c        iglob=glob(i)
c        if(.not. qtbdata(iglob)%currently_loc) then
c          write(0,*) rank,"atom",iglob,"should be local"
c          call fatal
c        endif
c      enddo

      if (deb_Path) write(*,*) '   << cleanup_qtbdata'
      end subroutine cleanup_qtbdata

c      subroutine gather_vad()
c      use atoms
c      use domdec
c      use spectra
c      use mpi
c      use inform, only: deb_Path
c      implicit none
c      integer :: i,j,k,iproc,ii,jj,ierr,iglob,nsendtot
c      real(4), allocatable :: data_send(:,:,:)
c      real(4), allocatable :: data_recv(:,:,:)
c      integer, allocatable :: iglob_send(:)
c      integer, allocatable :: iproc_send(:),iproc_recv(:)
c      integer :: nproc_send, maxsend
c      integer :: nproc_recv, maxrecv,nsend_capture
c      integer :: req_s(nproc),req_r(nproc),req_size(2*nproc)
c      integer :: nsend(nproc),nrecv(nproc),procnum(nproc)
c     &         ,nsend_bis(nproc)
c      integer, parameter :: mptag_size=0,mptag_data=1
c
c      if(nproc==1) return
c      if (deb_Path) write(*,*) '   >> gather_vad'
c
c!$acc update host(repart) async
c      ! COUNT DATA TO SEND 
c      nsend(:)=0
c      allocate(iglob_send(n))
c
c      nsendtot=0
c!$acc wait
c      do i=1,n
c        iproc=repart(i)
c        if(qtbdata_allocated(i) 
c     &   .AND. iproc/=rank) then          
c          nsend(iproc+1)=nsend(iproc+1)+1
c          nsendtot=nsendtot+1
c          iglob_send(nsendtot)=i
c        endif
c      enddo
c      nproc_send=COUNT(nsend>0)
c      maxsend=maxval(nsend)
c
c      !SEND AND RECEIVE SIZES
c      req_size(:)=MPI_REQUEST_NULL
c      ii=0
c      procnum(:)=-1
c      nrecv(:)=0
c      allocate(iproc_send(nproc_send))
c      do iproc=1,nproc
c        if(iproc-1==rank) CYCLE
c        CALL MPI_ISEND(nsend(iproc),1,MPI_INTEGER
c     &   , iproc-1, mptag_size , COMM_TINKER, req_size(iproc),ierr)
c        CALL MPI_IRECV(nrecv(iproc),1,MPI_INTEGER
c     &   , iproc-1, mptag_size ,COMM_TINKER,req_size(nproc+iproc),ierr)
c        if(nsend(iproc)==0) CYCLE
c        !TWO WAY MAP OF PROCS TO SEND DATA
c        ii=ii+1
c        procnum(iproc)=ii
c        iproc_send(ii)=iproc
c      enddo
c
c       !WAIT FOR SIZES
c      call MPI_Waitall(2*nproc,req_size
c     &       ,MPI_STATUSES_IGNORE, ierr)
c
c      !ORGANIZE DATA TO SEND
c      nsend_bis(:)=0
c      maxrecv=maxval(nrecv)
c      nproc_recv=COUNT(nrecv>0)
c      allocate(data_send(3*nseg+1,maxsend,nproc_send))
c      allocate(data_recv(3*nseg+1,maxrecv,nproc_recv))
c
c!$acc data create(data_send,data_recv) 
c!$acc& copyin(iglob_send,procnum,nsend_bis,nrecv)
c
c!$acc parallel loop async default(present) 
c      do i=1,nsendtot
c        iglob=iglob_send(i)
c        iproc=repart(iglob)+1
c!$acc atomic capture
c        nsend_bis(iproc)=nsend_bis(iproc)+1
c        nsend_capture=nsend_bis(iproc)
c!$acc end atomic
c        data_send(1,nsend_capture,procnum(iproc))=iglob
c!$acc loop seq collapse(2)
c        do k=1,nseg; do j=1,3
c          data_send(3*(k-1)+j+1,nsend_capture,procnum(iproc))
c     &       =qtbdata(iglob)%vad(j,k)
c        enddo; enddo
c      enddo
c
c      !RECEIVE DATA
c      req_r(:)=MPI_REQUEST_NULL
c      req_s(:)=MPI_REQUEST_NULL
c
c      allocate(iproc_recv(nproc_recv))
c      jj=0
c
c!$acc wait
c!$acc host_data use_device(data_send, data_recv) 
c      ! SEND DATA    
c      do ii=1,nproc_send 
c        iproc=iproc_send(ii)
c        CALL MPI_ISEND(data_send(1,1,ii),(3*nseg+1)*nsend(iproc)
c     &   ,MPI_REAL4, iproc-1
c     &   ,mptag_data, COMM_TINKER,req_s(ii),ierr)
c      enddo
c      
c      do iproc=1,nproc
c        if(iproc-1==rank .OR. nrecv(iproc)==0) CYCLE
c        jj=jj+1
c        iproc_recv(jj)=iproc
c        CALL MPI_IRECV(data_recv(1,1,jj),(3*nseg+1)*nrecv(iproc)
c     &   ,MPI_REAL4, iproc-1
c     &   ,mptag_data,COMM_TINKER,req_r(jj),ierr)
c      enddo
c!$acc end host_data
c
c      !WAIT FOR DATA
c      call MPI_Waitall(ii,req_s
c     &       ,MPI_STATUSES_IGNORE, ierr)
c      call MPI_Waitall(jj,req_r
c     &       ,MPI_STATUSES_IGNORE, ierr)
c
c      
c      !REDUCE DATA
c!$acc parallel loop async default(present)
c!$acc&         vector_length(512) copyin(iproc_recv)
c      do ii=1,nproc_recv
c!$acc loop vector collapse(3)
c        do i=1,nrecv(iproc_recv(ii))
c          do k=1,nseg; do j=1,3
c            iglob=nint(data_recv(1,i,ii))
c!$acc atomic update
c            qtbdata(iglob)%vad(j,k)=qtbdata(iglob)%vad(j,k)
c     &         +data_recv(3*(k-1)+j+1,i,ii)
c          enddo; enddo
c        enddo
c      enddo
c!$acc end data
c
c      deallocate(data_recv, data_send)
c
c      if (deb_Path) write(*,*) '  << gather_vad'
c
c      end subroutine gather_vad
c
c      subroutine gather_noise()
c      use atoms
c      use domdec
c      use spectra
c      use mpi
c      use inform, only: deb_Path
c      implicit none
c      integer :: i,j,k,iproc_s,iproc_r,ierr,nseg_comm
c      integer :: maxsend,maxrecv,nproc_send,nproc_recv
c      integer :: nsendtot,ii,jj,iproc,iglob,nsend_capture
c      integer :: req_s(nproc),req_r(nproc)
c      integer :: nsend(nproc),nrecv(nproc),procnum(nproc)
c     &             ,nsend_bis(nproc)
c      integer, allocatable :: iglob_send(:)
c      integer, allocatable :: iproc_send(:),iproc_recv(:)
c      real(4), allocatable :: data_send(:,:,:),data_recv(:,:,:)
c
c      if(nproc==1) return
c      if (deb_Path) write(*,*) '   >> gather_noise'
c
c!$acc update host(repart) async
c      nseg_comm=2*nseg
c
c      !COUNT DATA TO SEND/RECEIVE
c      nrecv(:)=0
c      nsend(:)=0
c      nsendtot=0
c      allocate(iglob_send(n))
c!$acc wait
c      do i=1,n
c        iproc_s=repartnoise(i)+1
c        iproc_r=repart(i)+1
c        if(iproc_r-1==rank .AND. iproc_s-1/=rank) then
c          !RECEIVE IF I HAVE AN ATOM FOR WHICH
c          ! I DID NOT GENERATE THE NOISE
c          nrecv(iproc_s)=nrecv(iproc_s)+1
c        elseif(iproc_s-1==rank .AND. iproc_r-1/=rank) then
c          !SEND IF I GENERATED THE NOISE FOR AN ATOM
c          ! I DON'T HAVE
c          nsend(iproc_r)=nsend(iproc_r)+1
c          nsendtot=nsendtot+1
c          iglob_send(nsendtot)=i
c        endif
c      enddo
c      maxsend=maxval(nsend)
c      maxrecv=maxval(nrecv)
c      nproc_send=COUNT(nsend>0)
c      nproc_recv=COUNT(nrecv>0)
c
c      !MAP OF PROCS TO SEND/RECEIVE DATA
c      ii=0
c      jj=0
c      procnum(:)=-1
c      allocate(iproc_send(nproc_send))
c      allocate(iproc_recv(nproc_recv))
c      iproc_send(:)=-1
c      iproc_recv(:)=-1
c      do iproc=1,nproc
c        if(iproc-1==rank) CYCLE
c        if(nsend(iproc)>0) then
c          ii=ii+1
c          procnum(iproc)=ii
c          iproc_send(ii)=iproc
c        endif
c        if(nrecv(iproc)>0) then
c          jj=jj+1
c          iproc_recv(jj)=iproc
c        endif
c      enddo
c
c      !ORGANIZE DATA TO SEND
c      allocate(data_send(3*nseg_comm+1,maxsend,nproc_send))
c      allocate(data_recv(3*nseg_comm+1,maxrecv,nproc_recv))
c      nsend_bis(:)=0
c!$acc data create(data_send,data_recv) 
c!$acc& copyin(nrecv,iproc_recv,iglob_send,procnum,nsend_bis)
c
c!$acc parallel loop async default(present) 
c      do i=1,nsendtot
c        iglob=iglob_send(i)
c        iproc=repart(iglob)+1
c!$acc atomic capture
c        nsend_bis(iproc)=nsend_bis(iproc)+1
c        nsend_capture=nsend_bis(iproc)
c!$acc end atomic
c        data_send(1,nsend_capture,procnum(iproc))=iglob
c!$acc loop seq collapse(2)
c        do k=1,2*nseg; do j=1,3
c          data_send(3*(k-1)+j+1,nsend_capture,procnum(iproc))
c     &       =qtbdata(iglob)%noise(j,k)
c        enddo; enddo
c      enddo      
c      
c      ! SEND DATA    
c      req_s(:)=MPI_REQUEST_NULL  
c      req_r(:)=MPI_REQUEST_NULL
c!$acc wait
c!$acc host_data use_device(data_send,data_recv)
c      do ii=1,nproc_send 
c        iproc=iproc_send(ii)
c        CALL MPI_ISEND(data_send(1,1,ii),(3*nseg_comm+1)*nsend(iproc)
c     &   ,MPI_REAL4, iproc-1
c     &   ,0, COMM_TINKER,req_s(ii),ierr)
c      enddo
c
c      do jj=1,nproc_recv
c        iproc=iproc_recv(jj)
c        CALL MPI_IRECV(data_recv(1,1,jj),(3*nseg_comm+1)*nrecv(iproc)
c     &   ,MPI_REAL4, iproc-1
c     &   ,0,COMM_TINKER,req_r(jj),ierr)
c      enddo
c!$acc end host_data
c
c      !WAIT FOR DATA
c      call MPI_Waitall(nproc_send,req_s,MPI_STATUSES_IGNORE,ierr)          
c      call MPI_Waitall(nproc_recv,req_r,MPI_STATUSES_IGNORE,ierr)      
c
c      !ORGANIZE RECEIVED DATA
c!$acc parallel loop async default(present)
c!$acc&         vector_length(512)
c      do ii=1,nproc_recv
c!$acc loop vector collapse(3)
c        do i=1,nrecv(iproc_recv(ii))
c          do k=1,2*nseg; do j=1,3
c            iglob=nint(data_recv(1,i,ii))
c            qtbdata(iglob)%noise(j,k)=data_recv(3*(k-1)+j+1,i,ii)
c          enddo; enddo
c        enddo
c      enddo
c!$acc end data
c
c!$acc exit data delete(data_recv,data_send) async
c      deallocate(data_recv,data_send)
c      if (deb_Path) write(*,*) '   << gather_noise'
c
c      end subroutine gather_noise

c-------------------------------------------------------------
c     adQTB OPTIMIZERS
      subroutine adapt_gamma_simple(itype,mCvv,Cvf,dFDR
     &   ,qtb_therm)
        use ascii, only: int_to_str
        use spectra
        use bath
        use beads
        use units
        use langevin
        implicit none
        integer, intent(in) :: itype
        real(r_p), intent(in) :: dFDR(:,:)
        real(r_p), intent(in) :: mCvv(:,:),Cvf(:,:)
        logical,intent(in) :: qtb_therm
        real(r_p) :: a1,gamma_ref,norm,gamma_disp,normv(ntype)
        real(r_p), allocatable, save :: dfdr_m(:,:,:),dfdr_norm(:,:)
        real(r_p) :: b1,bias_correct1,smooth,dFDRs(nad,nbeads)
        integer, allocatable,save :: cnt_type(:)
        integer :: cnt
        integer :: i,j,l,ibead
        real(r_p) :: dfdr_norm_instant

        if(.not. allocated(dfdr_m)) then
          allocate(dfdr_m(nad,nbeads,ntype))
          dfdr_m(:,:,:)=0.d0
        endif
        if(.not. allocated(dfdr_norm)) then
          allocate(dfdr_norm(nbeads,ntype))
          dfdr_norm(:,:)=0
        endif
        if(.not. allocated(cnt_type)) then
          allocate(cnt_type(ntype))
          cnt_type(:)=0
        endif
        cnt_type(itype)=cnt_type(itype)+1
        cnt=cnt_type(itype)

        dFDRs(1:nad,:)=dFDR(1:nad,:)
        smooth=adqtb_smooth_type(adqtb_type_map(itype))
        if(smooth>0) then
          do ibead=1,nbeads
            call window_smooth_spectrum(dFDRs(:,ibead),smooth)
          enddo
        endif

        if(adqtb_tau_adapt>0) then
          b1 = exp(-Tseg*adqtb_avg_adapt/adqtb_tau_adapt)
        else
          b1=0.d0
        endif
        dfdr_m(:,:,itype)=b1*dfdr_m(:,:,itype) + (1.d0-b1)*dFDRs
        bias_correct1=1.d0 - b1**cnt

        do ibead=1,nbeads
          dfdr_norm_instant=NORM2(dfdr_m(:,ibead,itype)/bias_correct1)
          dfdr_norm(ibead,itype)=dfdr_norm(ibead,itype) 
     &      + (dfdr_norm_instant-dfdr_norm(ibead,itype))/cnt

          !a1=A_gamma*dt_spec*nseg/(boltzmann*kelvin)
          a1=dt_spec*nseg*gamma_friction(ibead)
     &           *A_gamma(adqtb_type_map(itype))
     &           /dfdr_norm(ibead,itype)  
          do j=1,nad  
            gamma_type(j,itype,ibead)=gamma_type(j,itype,ibead) !*dFDR_type(j,i) )
     &                      -a1*dfdr_m(j,ibead,itype)/bias_correct1
          enddo
        enddo

        gamma_type = max(gammar_min,gamma_type)

      end subroutine adapt_gamma_simple

      subroutine adapt_gamma_ratio(itype,mCvv,Cvf,dFDR
     &   ,qtb_therm)
        use spectra
        use atoms
        use langevin
        use ascii, only: int_to_str
        use katoms, only: weight
        use units
        use bath
        use beads
        implicit none
        integer, intent(in) :: itype
        real(r_p), intent(in) :: dFDR(:,:)
        real(r_p), intent(in) :: mCvv(:,:),Cvf(:,:)
        logical,intent(in) :: qtb_therm
        real(r_p) :: b1,tau,sqrtm,w,b_avg,avg_correct,a1
        real(r_p) :: avg_correct_prev
        real(r_p) :: adapt_correct,tmp
        real(r_p) :: gamma_new(nad),disp,smooth
        real(r_p), allocatable, save :: gamma_curr(:,:,:)
        real(r_p), allocatable, save :: Cvf_m(:,:,:)
        real(r_p), allocatable, save :: mcvv_m(:,:,:)
        integer, allocatable, save :: nadapt(:)
        integer :: i,j,l,ibead
        integer, allocatable, save :: cnt_type(:)
        real(r_p) :: mCvvs(nad,nbeads), Cvfs(nad,nbeads)
        integer :: cnt
        real(r_p) :: tau_avg

        if(.not. allocated(cnt_type)) then
          allocate(cnt_type(ntype))
          cnt_type(:)=0
          allocate(gamma_curr(nad,ntype,nbeads))
          allocate(mcvv_m(nad,nbeads,ntype))
          allocate(Cvf_m(nad,nbeads,ntype))
          allocate(nadapt(ntype))
          nadapt=0
          mcvv_m=0.d0
          Cvf_m=0.d0
          if(adqtb_correct_bias) then
            gamma_curr=0.d0
          else
            gamma_curr=gamma_type
          endif
        endif
        cnt_type(itype)=cnt_type(itype)+1
        cnt=cnt_type(itype)

        Cvfs(1:nad,:)=Cvf(1:nad,:)
        mCvvs(1:nad,:)=mCvv(1:nad,:)
        smooth=adqtb_smooth_type(adqtb_type_map(itype))
        if(smooth>0) then
          do ibead=1,nbeads
            call window_smooth_spectrum(Cvfs(:,ibead),smooth)
            call window_smooth_spectrum(mCvvs(:,ibead),smooth)
          enddo
        endif
        
        w=weight(adqtb_type_map(itype))
        sqrtm=merge(sqrt(w),1.d0,w>0)   
        tau=adqtb_tau_avg(adqtb_type_map(itype))
        if(tau>0) then
          b_avg = exp(-Tseg*adqtb_avg_adapt/tau)
          avg_correct = 1.d0 - b_avg**cnt
          if(cnt>1) then
            avg_correct_prev = 1.d0 
     &          - b_avg**(cnt-1)
          else
            avg_correct_prev=1.d0
          endif
        else
          tau=0.
          b_avg=0.d0
          avg_correct=1.d0
          avg_correct_prev=1.d0
        endif 
        do ibead=1,nbeads
          do j=1,nad
            mCvv_m(j,ibead,itype)=mCvv_m(j,ibead,itype)*b_avg
     &          + (1.d0-b_avg)*mCvvs(j,ibead)
            Cvf_m(j,ibead,itype)=b_avg*Cvf_m(j,ibead,itype) 
     &          + (1.d0-b_avg)*Cvfs(j,ibead)
          enddo
        
        !write(*,*) itype,Tseg*adqtb_avg_adapt*cnt,tau
c         if(Tseg*adqtb_avg_adapt*cnt>=tau) then
c           if(nadapt(itype)==0) then
c             write(*,*) "Started adapting Type_"
c     &         //trim(adqtb_name(itype))
c     &         //"_"//int_to_str(adqtb_type_map(itype))
c           endif
          
          gamma_new=Cvf_m(:,ibead,itype)/avg_correct
     &            /(mCvv_m(:,ibead,itype)/avg_correct+1.d-8)
          if(Tseg*adqtb_avg_adapt*cnt<tau) then
            b1 = Tseg*adqtb_avg_adapt*cnt/tau
            gamma_new = gamma_new*b1
     &             +gamma_type_init(:,itype,ibead)*(1.d0-b1)
          endif

          where(gamma_new<gammar_min)
            gamma_new=gammar_min
          endwhere    
          if(adqtb_tau_adapt>0 .and.
     &        Tseg*adqtb_avg_adapt*cnt>=tau)
     &    then  
            nadapt(itype)=nadapt(itype)+1          
            tau=adqtb_tau_adapt *sqrtm 
c     &           *real(n,8)/real(n_of_type(itype),8)
            b1 = exp(-Tseg*adqtb_avg_adapt/tau)            
            if(adqtb_correct_bias) then
              adapt_correct = (1.d0 - b1**nadapt(itype))
            else
              adapt_correct = 1.d0
            endif            
            !write(*,*) "adapt",itype," with tau=",tau,"ps"
          else            
            b1=0.d0
            adapt_correct = 1.
          endif            
          
          gamma_curr(:,itype,ibead) = b1*gamma_curr(:,itype,ibead) 
     &              + (1.d0 - b1)*gamma_new
          gamma_type(:,itype,ibead)=gamma_curr(:,itype,ibead)
     &                               /adapt_correct
c         endif  
        enddo

      end subroutine adapt_gamma_ratio

      subroutine window_smooth_spectrum(s,dw)
        use spectra
        implicit none 
        real(r_p), intent(inout) :: s(:)
        real(r_p), intent(in) :: dw
        integer :: i,j,jstart,jend
        real(r_p) :: w,w0,norm,b1
        real(r_p), allocatable :: s_in(:)

        allocate(s_in(size(s)))
        s_in(:)=s(:)
        do i=1,nad
          w0=(i-1)*domega
          jstart = max(1,int((w0-dw)/domega))
          jend = min(nad,1+int((w0+dw)/domega))
          norm=0
          s(i)=0.d0
          do j=jstart,jend  
          !do j=1,nad
            w=((j-1)*domega-w0)/dw
            b1=merge(exp(-1.d0/(1+w)), 0.d0, 1+w>0)
            b1=merge(b1*exp(-1.d0/(1-w)), 0.d0, 1-w>0)
            norm=norm+b1
            s(i) = s(i) + s_in(j)*b1
          enddo 
          s(i)=s(i)/norm
        enddo
        deallocate(s_in)

      end subroutine window_smooth_spectrum

      subroutine generate_corr_pot()
        use deconvolution
        use domdec
        use spectra
        use mpi
        use units
        use bath
        use beads
        use ascii, only: int_to_str
        use langevin
        implicit none
        real(r_p), allocatable :: omsym(:),s0sym(:)
        real(r_p), allocatable :: sout(:),srec(:)
        real(r_p) :: f_omega,w
        integer :: k,ierr,ibead

        if(rank==0) then
          if(ranktot==0) write(*,*) "--- Generating corr_pot.dat ---"
          allocate(omsym(2*nad-1),s0sym(2*nad-1))
          allocate(sout(2*nad-1),srec(2*nad-1))
          do ibead=1,nbeads
            omsym(nad)=0.d0
            s0sym(nad)=1.d0
            do k=1,nad-1
              w=k*domega
              omsym(nad+k)=w
              omsym(nad-k)=-w
              f_omega=QTB_kernel(w,piqtb,ibead)
     &           /(1.0d0+exp((w-omegacut)/omegasmear))
              s0sym(nad+k)=f_omega !/kubo_fact(k)
              s0sym(nad-k)=f_omega !/kubo_fact(k)
            enddo
            call deconvolute_spectrum(s0sym,omsym
     &         ,niter_deconvolution,.true.,real(1.0d-10,r_p)
     &         ,gamma_friction(ibead),boltzmann*kelvin
     &         ,kernel_lorentz_pot,sout,srec,verbose=.FALSE.)
            do k=1,nad
              corr_pot_ratio(k,ibead)=(sout(nad+k-1)-s0sym(nad+k-1))
     &                              /s0sym(nad+k-1)
            enddo
          enddo
          open(315,file='corr_pot.dat')
          write(315,'(A,3E17.8)') "# ",Tseg, gamma,omegacut
          do k=1,nad
            write(315,'('//int_to_str(1+nbeads)//'E17.8)')
     &         (k-1)*domega*cm1, corr_pot_ratio(k,:)
          enddo
          close(315)
        endif
        call MPI_BCAST(corr_pot_ratio,nbeads*nad,MPI_RPREC
     &     ,0,COMM_TINKER,ierr)
      end subroutine generate_corr_pot

      subroutine compute_corr_fact(mCvv_in)
        use deconvolution
        use bath
        use langevin
        use spectra
        use mpi
        use domdec
        use units
        implicit none
        real(r_p), intent(in) :: mCvv_in(:,:,:)
        real(r_p) :: mCvv(size(mCvv_in,1),size(mCvv_in,2))
        real(r_p), allocatable :: s0sym(:),omsym(:),sdec(:),srec(:)
        real(r_p), allocatable :: corr(:)
        real(r_p) :: rec_ratio,rec_ratio_min
        integer, save :: isame=0
        real(r_p), save :: corr_fact_prev(3)=1.0d0
        real(r_p), parameter :: thr=1.d-4
        integer :: k,ierr,jend,j

        if (piqtb) return

        if(rank==0) then

          if (anisotrop) then
            jend=3
            mCvv=mCvv_in(:,:,1)
          else
            mCvv(:,1)=sum(mCvv_in(:,1:3,1),dim=2)
            jend=1
          endif

          allocate(s0sym(2*nad-1))
          allocate(omsym(2*nad-1))
          allocate(sdec(2*nad-1))
          allocate(srec(2*nad-1))
          allocate(corr(nad))
          corr(:)=1.0d0
          if(corr_pot_corr_fact .and. corr_pot) then
            corr(1:nad) = 1.0d0+corr_pot_ratio(1:nad,1)
          endif
          rec_ratio_min=1.
          do j=1,jend 
            omsym(nad)=0.d0
            s0sym(nad)=mCvv(1,j)/corr(1)
            do k=1,nad-1
              omsym(nad+k)=k*domega
              omsym(nad-k)=-k*domega
              s0sym(nad+k)=mCvv(k+1,j)*kubo_fact(k+1)/corr(k+1)
              s0sym(nad-k)=mCvv(k+1,j)*kubo_fact(k+1)/corr(k+1)
            enddo
            call deconvolute_spectrum(s0sym,omsym
     &           ,niter_deconvolution,.false.,real(1.0d-10,r_p)
     &           ,gamma,boltzmann*kelvin,kernel_lorentz
     &          ,sdec,srec,verbose=.FALSE.)
  
            s0sym(nad)=s0sym(nad)*corr(1)
            srec(nad)=srec(nad)*corr(1)
            do k=1,nad-1
              s0sym(nad+k)=s0sym(nad+k)/kubo_fact(k+1)*corr(k+1)
              s0sym(nad-k)=s0sym(nad-k)/kubo_fact(k+1)*corr(k+1)
              srec(nad+k)=srec(nad+k)/kubo_fact(k+1)*corr(k+1)
              srec(nad-k)=srec(nad-k)/kubo_fact(k+1)*corr(k+1)
              sdec(nad+k)=sdec(nad+k)/kubo_fact(k+1)
              sdec(nad-k)=sdec(nad-k)/kubo_fact(k+1)
            enddo
            rec_ratio=sum(s0sym)/sum(srec)
            corr_fact_qtb(j)=sum(s0sym)/sum(sdec)
            if (rec_ratio<rec_ratio_min) then
              rec_ratio_min=rec_ratio
            endif
          enddo
          if(rec_ratio_min<0.95) then
            write(*,*) "Warning: reconvolution error is too high"
     &        //", corr_fact_qtb was not updated" 
            corr_fact_qtb=corr_fact_prev
          endif
          if (.not. anisotrop) corr_fact_qtb(1:3)=corr_fact_qtb(1)
          if (QTB_verbose) then
            write(*,'(A,3F7.3,A,1F5.3,A,1E14.7)') "corr_fact_qtb= "
     &      ,corr_fact_qtb(:),"  reconv_ratio= ",rec_ratio_min
     &      ,"  diff_to_prev",maxval(abs(corr_fact_qtb-corr_fact_prev))
          endif
        endif

        CALL MPI_BCAST(corr_fact_qtb,3,MPI_RPREC,0,COMM_TINKER,ierr)
        if(maxval(abs(corr_fact_qtb-corr_fact_prev))<thr) then
          isame=isame+1
        else
          isame=0
        endif
        if(isame>10) then
          if(rank==0 .and. QTB_verbose) write(*,*) 
     &       "corr_fact_qtb is converged",
     &      " (it did not change for 10 consecutive segments)."
          corr_fact_converged=.TRUE.
        endif
        corr_fact_prev=corr_fact_qtb

        if(rank==0) then
          open(315,file='corr_fact_qtb_restart.out')
          write(315,'(3f15.6,L)') corr_fact_qtb(:),corr_fact_converged
          close(315)
        endif
      end subroutine compute_corr_fact


      subroutine broadcast_piqtb_noise()
      use beads
      use mpi
      use beads
      use domdec
      use inform, only: deb_path
      use spectra
      implicit none
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk,l
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      integer :: nbeadslocmax, nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      integer :: cshift,nelt
      integer :: send_vel_
      real(4), allocatable :: bufferpi_s(:,:,:,:)
      real(4), allocatable :: bufferpi_r(:,:,:,:,:)

      if(nproc_polymer==1 .or. nproc==1) return
      reqs(:)=MPI_REQUEST_NULL

      if(deb_path) then
!$acc wait
        write(*,*) '>> broadcast_piqtb_noise'
      endif

      nelt=3*nseg

      maxloc=nlocpis(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(bufferpi_s(3,nelt,nbeads,maxloc))
      allocate(bufferpi_r(3,nelt,nbeads,maxloc
     &     ,nproc_polymer))
      
!$acc enter data create(bufferpi_r,bufferpi_s) async

        !PUT DATA INTO BUFFER ARRAYS

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)

!$acc parallel loop collapse(3) async
        DO k=1,nbeads; DO iloc=1,ilocend-ilocbeg+1; do l=1,3*nseg
          i=glob(ilocbeg+iloc-1)
          if (l<=2*nseg) then
!$acc loop seq
            DO j=1,3
            bufferpi_s(j,l,k,iloc)=qtbdata(i)%noise(j,l,k)
            ENDDO
          else
!$acc loop seq
            DO j=1,3
            bufferpi_s(j,l,k,iloc)=qtbdata(i)%fad(j,l-2*nseg,k)
            ENDDO
          endif
        ENDDO; ENDDO ; ENDDO

!$acc wait
!$acc host_data use_device(bufferpi_s,bufferpi_r)
      ii=0
      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        call MPI_ISEND(bufferpi_s
     &     ,3*nelt*(ilocend-ilocbeg+1)*nbeads
     &      ,MPI_REAL4
     &     ,iproc-1, mptag_data, COMM_POLYMER, reqs(ii),ierr)
      ENDDO

      DO iproc=1,nproc_polymer
        if(rank_polymer==iproc-1) CYCLE
        ii=ii+1
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        call MPI_IRECV(bufferpi_r(1,1,1,1,iproc)
     &   ,3*nelt*(ilocend-ilocbeg+1)*nbeads
     &   ,MPI_REAL4,iproc-1,mptag_data, COMM_POLYMER
     &   ,reqs(ii),ierr)
      ENDDO
!$acc end host_data

      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

!$acc parallel loop collapse(3) async
        DO k=1,nbeads; DO iloc=1,nloc; DO l=1,3*nseg
          iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
          if(rank_polymer==iproc-1) CYCLE
          ilocbeg=ilocpi_beg(iproc)
          ilocend=ilocpi_end(iproc)
          i=glob(iloc)
          kk=iloc-ilocbeg +1
          if (l<=2*nseg) then
!$acc loop seq
            DO j=1,3
            qtbdata(i)%noise(j,l,k)=bufferpi_r(j,l,k,kk,iproc)
            ENDDO
          else
!$acc loop seq
            DO j=1,3
            qtbdata(i)%fad(j,l-2*nseg,k)=bufferpi_r(j,l,k,kk,iproc)
            ENDDO
          endif
        ENDDO; ENDDO; ENDDO

c!$acc wait
c      call MPI_BARRIER(COMM_POLYMER,ierr)

!$acc exit data delete(bufferpi_s,bufferpi_r)
      deallocate(bufferpi_s,bufferpi_r)

      if(deb_path) then
!$acc wait
        write(*,*) '<< broadcast_piqtb_noise'
      endif

      end subroutine broadcast_piqtb_noise

      subroutine reduce_piqtb_vad()
      use mpi
      use beads
      use domdec
      use inform, only: deb_path
      use spectra
      implicit none
      integer :: ibead,k,i,j,iglob,iproc,maxsend,nsend,maxloc
      integer :: maxrecv,nproc_recv
      integer :: isend,ii,ierr,ns,nr,iloc,kk,l
      integer, parameter :: mptag_size=0,mptag_data=1
      integer :: reqs(2*nproc_polymer)
      LOGICAL :: flag
      integer :: nbeadslocmax, nlocproc
      integer :: ibeadbeg,ibeadend,nu
      integer :: ilocend,ilocbeg
      integer :: cshift,nelt,ineighbor
      integer :: send_vel_
      real(4), allocatable :: buffer(:,:,:,:)

      if(nproc_polymer==1 .or. nproc==1) return
      reqs(:)=MPI_REQUEST_NULL

      if(deb_path) then
!$acc wait
        write(*,*) '>> reduce_piqtb_vad'
      endif

      nelt=nseg

      maxloc=nlocpis(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(buffer(3,nelt,nloc,nbeads))
!$acc data create(buffer) async

        !PUT DATA INTO BUFFER ARRAYS

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)

!$acc parallel loop collapse(4) async
      do k=1,nbeads; do i=1,nloc; do l=1,nseg; do j=1,3
        iglob=glob(i)
        buffer(j,l,i,k)=qtbdata(iglob)%vad(j,l,k)
      enddo; enddo; enddo; enddo

!$acc wait
!$acc host_data use_device(buffer)
      call MPI_ALLREDUCE(MPI_IN_PLACE,buffer
     &       ,3*nelt*nbeads*nloc 
     &       ,MPI_REAL4,MPI_SUM
     &       ,COMM_POLYMER,ierr)
!$acc end host_data

!$acc parallel loop collapse(4) async
      do k=1,nbeads; do i=1,nloc; do l=1,nseg; do j=1,3
        iglob=glob(i)
        qtbdata(iglob)%vad(j,l,k)=buffer(j,l,i,k)
      enddo; enddo; enddo; enddo

!$acc end data
      deallocate(buffer)

      if(deb_path) then
!$acc wait
        write(*,*) '<< reduce_piqtb_vad'
      endif

      end subroutine reduce_piqtb_vad

      end module qtb
