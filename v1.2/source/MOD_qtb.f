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
      module qtb 
      use spectra, only: compteur_type
      implicit none
c      integer, allocatable :: repartnoise(:)
      integer :: skipseg
      real*8 :: omegasmear,gammar_min
      real*8, allocatable :: Htilde(:,:,:)
      real*8, allocatable :: corr_pot_ratio(:,:)
      logical :: noQTB
      logical :: corr_pot
      logical :: register_spectra
      logical :: qtb_thermostat=.FALSE.
      type(compteur_type) :: compteur,compteur_bis
      integer :: plan_qtb, nloc_plan=-1
      logical :: QTB_verbose=.FALSE.
      logical :: piqtb
c
c     adQTB specific variables
c
      character*11 :: adqtb_optimizer
      character*11,allocatable :: adqtb_optimizer_type(:)

      integer :: adqtb_avg_adapt=1
      real*8 :: adqtb_tau_avg_default, adqtb_tau_adapt
      logical :: adqtb_correct_bias, adqtb_win_corr
      real*8,allocatable :: adqtb_tau_avg(:),A_gamma(:)

      logical adaptive_qtb
      real*8 :: A_gamma_default
      real*8 :: corr_fact_qtb(3)=1.d0
      logical :: use_corr_fact_qtb
      logical :: corr_fact_converged=.FALSE.
      logical :: corr_pot_corr_fact=.TRUE.
      real*8, allocatable :: gamma_type(:,:,:)
      real*8, allocatable :: gamma_type_init(:,:,:)
      real*8, allocatable :: mCvv_average_type(:,:,:)
      real*8, allocatable :: mCvv_average_sum(:,:,:)
      real*8, allocatable :: Cvf_average_type(:,:,:)
      real*8, allocatable :: mCvv_mean(:,:,:)
      real*8, allocatable :: Cvf_mean(:,:,:)
      real*8, allocatable :: dFDR_mean(:,:,:)
      real*8, allocatable :: dFDR_average_type(:,:,:)
      integer, allocatable :: adqtb_type(:),n_of_type(:)
      character*3, allocatable :: adqtb_name(:)
      integer, allocatable :: adqtb_type_map(:)
      integer :: ntype
      logical :: save_gamma_history=.FALSE.

      logical :: adqtb_start_classical=.FALSE.

      real*8 :: adqtb_smooth_default
      real*8, allocatable :: adqtb_smooth_type(:)

      abstract interface
          subroutine qtb_optimizer_type(itype,mCvv,Cvf,dFDR
     &   ,qtb_therm)
        implicit none
        integer, intent(in) :: itype
        real*8, intent(in) :: dFDR(:,:)
        real*8, intent(in) :: mCvv(:,:),Cvf(:,:)
        logical,intent(in) :: qtb_therm
        end subroutine qtb_optimizer_type
      end interface


      type qtbdat_t
        real, allocatable :: noise(:,:,:)
        real, allocatable :: fad(:,:,:), vad(:,:,:)
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
     &          , use_corr_fact_qtb
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
      use spectra
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
      implicit none
      real*8, intent(in) ::  dt
      integer i,j,k,l,iglob,ierr,ibead,uu
      integer compteur_lines
      real*8 :: freq
      logical :: adqtb_restart,corr_pot_found
      character*120 :: gamma_restart_file
      real*8 :: Tseg_file, gamma_file, omegacut_file,omega
      integer, allocatable :: tmp(:)
      character :: bufchar
      logical :: corr_fact_restart
      character*11 :: opttype
      character*3 :: numberbeads
      character*240 :: exten

      interface
         function normal ()
         real*8 normal
         end function
      end interface

      call read_qtb_keys()
      
      if (adaptive_qtb) register_spectra=.TRUE.
      call initialize_spectra(dt)

      if(allocated(gamma_friction)) then
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
      allocate(qtbdata(n))
      allocate(qtbdata_allocated(n))
      do i=1,n
        qtbdata_allocated(i)=.FALSE.
c        qtbdata(i)%currently_loc=.FALSE.
      enddo

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
      adqtb_name(:) = '   '
      adqtb_type_map = -1
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
      if(ranktot==0)
     &   write(0,*) "Number of adQTB types:",ntype

      if(ranktot==0 .AND. use_corr_fact_qtb) then
        write(*,*) 'WARNING - Using Kinetic correction'
        write(*,*) 'corr_fact_qtb=', corr_fact_qtb(1)
      endif
c
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
          open(newunit=uu,file='corr_pot.dat')
          read(uu,*) bufchar,Tseg_file, gamma_file,omegacut_file
          if((Tseg_file.ne.Tseg).or.(omegacut_file.ne.omegacut) 
     &            .or.(gamma_file.ne.gamma)) then
                  ! NOTE: dangerous comparison of real types
c            write(*,*) 'ERROR, NOT THE SAME PARAMETERS WHILE USING
c     & THE POTENTIAL CORRECTION'
c            call fatal
            corr_pot_found=.FALSE.
          else  
            do i=1,nad
              read(uu,*) omega,corr_pot_ratio(i,:)
            enddo
          endif
          close(uu)
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
        if(ranktot==0)
     &   write(*,*) "Using gamma_restart.out to initialize QTB"
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
      allocate(Htilde(3*nseg,ntype,nbeads))
      call adHnoise(.TRUE.)
c
c     Generate colored noise
c
      do i=1,nloc
        iglob = glob(i)
        do ibead=1,nbeads;do k=1,2*nseg; do j=1,3
          qtbdata(iglob)%noise(j,k,ibead)=normal()
        enddo; enddo;enddo
      enddo        
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

      do i = 1, nkey
        ios = 0
        next = 1
        record = keyline(i)
        call gettext (record,keyword,next)
        call upcase (keyword)
        string = record(next:240)
        
        select case (trim(keyword))
        case ('A_GAMMA')
          read (string,*, iostat=ios) a_gamma_default
        case ('SKIPSEG')
          read (string,*, iostat=ios) skipseg
        case ('NOQTB')
          noQTB = .true.
        case ('GAMMA_HISTORY')
          save_gamma_history = .true.
        case ('QTB_VERBOSE')
          QTB_verbose=.TRUE.
        case ('CORR_FACT_QTB')
          read (string,*, iostat=ios) corr_fact_qtb(1)
          corr_fact_qtb(:)=corr_fact_qtb(1)
          use_corr_fact_qtb=.TRUE.
        case ('CORR_FACT_NO_POT')
          corr_pot_corr_fact=.FALSE.
        case ('NO_CORR_POT')
          corr_pot = .false.
        case ('REGISTER_SPECTRA')
          register_spectra = .true.
        case ('ADQTB_OPTIMIZER')
          call getword (record,adqtb_optimizer,next)
          call upcase (adqtb_optimizer)
        case ('ADQTB_START_CL')
          adqtb_start_classical = .true.
        case ('ADQTB_TAU_AVG')
          read (string,*, iostat=ios) adqtb_tau_avg_default
        case ('ADQTB_TAU_ADAPT')
          read (string,*, iostat=ios) adqtb_tau_adapt
        case ('ADQTB_BIAS_CORR')
          adqtb_correct_bias=.TRUE.
        case ('ADQTB_AVG_ADAPT')
          read (string,*, iostat=ios) adqtb_avg_adapt
        case ('ADQTB_WIN_CORR')
          adqtb_win_corr=.TRUE.
        case ('ADQTB_SMOOTH')
          read (string,*, iostat=ios) adqtb_smooth_default
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
      real*8 :: tmp,smooth
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

        if (keyword(1:15) .eq. 'ADQTB_MODIFIER ') then
          smooth=-1
          tmp=-1
          read(string,*,iostat=ios) k,opttype,tmp,smooth

          if(ios/=0) then
            write(*,*) "Error: ADQTB_MODIFIER keyword is not",
     &          " properly defined"
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
        endif
    
      enddo

      adqtb_smooth_type(:)=adqtb_smooth_type(:)/cm1

      end subroutine initialize_adaptation_parameters

      subroutine allocate_qtbdata(iglob)
        use domdec
        use spectra
        use beads
        implicit none
        integer, intent(in) :: iglob
        integer :: i,j

        allocate(qtbdata(iglob)%noise(3,2*nseg,nbeads))
        allocate(qtbdata(iglob)%fad(3,nseg,nbeads))
        if(register_spectra) then
          allocate(qtbdata(iglob)%vad(3,nseg,nbeads))
          qtbdata(iglob)%vad(:,:,:)=0.
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
          deallocate(qtbdata(iglob)%vad)
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
      implicit none
      logical, intent(in) :: write_kernel
      integer :: i,j,k,l,iglob,ibead,uu
      real*8 :: C_omega   !sinus cardinal
      real*8 :: f_omega   !smooth cutoff (Fermi function)
      real*8 :: theta_tilde
      real*8 :: h
      real*8 :: omega
      real*8 :: t,u
      real*8 :: g_ratio, gammak

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
      use spectra
      implicit none
      integer :: i,j,k,l,ierr,nloc_,ibegbeg,iendend, ibead
      integer :: iglob
      complex(FFT1D_s), allocatable :: s_fft(:),workfft(:)
      interface
         function normal ()
         real*8 normal
         end function
      end interface

      allocate(s_fft(3*nseg),workfft(3*nseg))

      if(ranktot==0 .and. QTB_verbose) 
     &     write(*,*) "--- Generating new colored noise segment ---"

      if(piqtb) then
        ibegbeg=ilocpi_beg(rank_polymer+1)
        iendend=ilocpi_end(rank_polymer+1)
      else
        ibegbeg=1
        iendend=nloc
      endif


      do ibead=1,nbeads
      do i=ibegbeg,iendend; iglob=glob(i); do j=1,3
        do k=1,nseg
          !COPY THE LAST SEGMENTS INTO TEMP ARRAY
          s_fft(k)=cmplx(qtbdata(iglob)%noise(j,k,ibead)
     &                  ,0.,kind=FFT1D_s)
          s_fft(k+nseg)=cmplx(qtbdata(iglob)%noise(j,k+nseg,ibead)
     &                  ,0.,kind=FFT1D_s)
          !SHIFT THE LAST SEGMENT
          qtbdata(iglob)%noise(j,k,ibead)
     &       =qtbdata(iglob)%noise(j,k+nseg,ibead)
          !GENERATE FRESH NOISE FOR THE NEXT SEGMENT
          qtbdata(iglob)%noise(j,k+nseg,ibead)=normal()
          !COPY THE NEW SEGMENT INTO TEMP ARRAY
          s_fft(k+2*nseg)=cmplx(qtbdata(iglob)%noise(j,k+nseg,ibead)
     &                     ,0.,kind=FFT1D_s)
        enddo
        !FOURIER TRANSFORM => to frequency
        call SPCFFT(s_fft,3*nseg,-1,workfft)
        !MULTIPLY THE NOISE BY THE KERNEL IN FOURIER SPACE
        s_fft(:)=s_fft(:)*Htilde(:,adqtb_type(iglob),ibead)
        !INV FOURIER TRANSFORM => back to time
        call SPCFFT(s_fft,3*nseg,1,workfft)
        !STORE THE RANDOM FORCE
        qtbdata(iglob)%fad(j,:,ibead)=
     &     real(s_fft(nseg+1:2*nseg))/(3*nseg)
     &      *sqrt(2.d0*mass(iglob)*gamma_friction(ibead)/dt_spec)
      enddo ; enddo   
      enddo    

      deallocate(s_fft)

      if (piqtb) call broadcast_piqtb_noise()

c      !UPDATE NOISE REPARTITION
c      if(nproc>1) then
c        repartnoise(:) = 0 
c        do i=1,nloc
c          repartnoise(glob(i)) = rank
c        enddo        
c        call MPI_BARRIER(COMM_TINKER,ierr)
c        call MPI_ALLREDUCE(MPI_IN_PLACE,repartnoise,n,MPI_INT,MPI_SUM,
c     $   COMM_TINKER,ierr) 
c      endif

      end subroutine convolseg


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
      implicit none
      integer i,j,k,l,ierr
      integer iglob,inoise
      real*8 a1

      if(piqtb .and. ranktot==0) then
        write(*,*) "Warning: qtbrandom not implemented for piqtb"
        call fatal
      endif

      inoise=compteur%i
      a1 = exp(-gamma*dt_spec)
c
c     save velocity and random force if you need spectra 
c      
      if (register_spectra) then
        do i=1,nloc ; iglob=glob(i);  do j=1,3
          qtbdata(iglob)%vad(j,inoise,1) = sqrt(a1)*v(j,iglob)
     $        +0.5d0*dt_spec*qtbdata(iglob)%fad(j,inoise,1)/mass(iglob)
        enddo ;enddo
      end if

c
c     apply QTB random force
c
      do i=1,nloc ; iglob=glob(i);  do j=1,3
          v(j,iglob)=a1*v(j,iglob) 
     $     +dt_spec*qtbdata(iglob)%fad(j,inoise,1)/mass(iglob)  
      enddo; enddo

      call QTB_update_state()

      end subroutine qtbrandom

      function get_qtb_compteur() result(inoise)
      implicit none
      integer :: inoise
      inoise=compteur%i
      end function get_qtb_compteur

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
          if(ranktot==0 .and. QTB_verbose) 
     &         write(*,*) "--- Activating QTB Piston ---"
        endif
        if(ranktot==0 .and. QTB_verbose)
     &      write(*,*) "--- Writing QTB spectra ---"
        call QTB_spectra_adaptation_gamma()
      endif

      call cleanup_qtbdata()
      call convolseg()

      end subroutine QTB_update_state

!
!     subroutine QTB_spectra_adaptation_gamma: compute Cvv/Cvf spectra and DeltaFDT
!                                              and adapt gammar if adaptive_qtb
      subroutine QTB_spectra_adaptation_gamma()
      use atmtyp
      use atoms
      use bath
      use beads
      use domdec
      use energi
      use freeze
      use katoms
      use langevin
      use math
      use mdstuf
      use moldyn
      use usage
      use units
      use mpi
      use spectra
      use ascii, only: int_to_str
      implicit none
      integer :: i,j,k,l,ierr,ierr1,ierr2,idof,ibead
      integer :: iglob, itype,uu,ibegbeg,iendend
      logical :: qtb_therm
      real*8,  allocatable :: mCvv_type(:,:,:,:), Cvf_type(:,:,:)
      real*8, allocatable :: dFDR_type(:,:,:),mCvv_sum(:,:,:)
      real*8 :: a1,fftfact,temp,mcvvsum_tmp,m
      complex(FFT1D_s), allocatable :: sv(:),sf(:),workfft(:)

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

      allocate(sv(3*nseg))
      allocate(sf(3*nseg))
      allocate(workfft(3*nseg))
      allocate(mCvv_type(3*nseg,ntype,3,nbeads))
      allocate(Cvf_type(3*nseg,ntype,nbeads))
      allocate(dFDR_type(nad,ntype,nbeads))
      allocate(mCvv_sum(nad,3,nbeads))
      mCvv_type=0.0d0
      Cvf_type=0.0d0

      if (piqtb) then
        call reduce_piqtb_vad()
        ibegbeg=ilocpi_beg(rank_polymer+1)
        iendend=ilocpi_end(rank_polymer+1)
      else
        ibegbeg=1
        iendend=nloc
      endif



      !call gather_vad()

      fftfact=dt_spec/real(3*nseg,8)
      do ibead=1,nbeads; do i=ibegbeg,iendend  
       iglob=glob(i)
       itype=adqtb_type(iglob)
       do j=1,3
         !PUT VELOCITIES AND FORCE IN TEMP ARRAY
         sv(1:nseg)=cmplx(qtbdata(iglob)%vad(j,:,ibead),0.,kind=FFT1D_s)
         sv(nseg+1:3*nseg)=cmplx(0.d0,0.d0,kind=FFT1D_s)
         sf(1:nseg)=cmplx(qtbdata(iglob)%fad(j,:,ibead),0.,kind=FFT1D_s)
         sf(nseg+1:3*nseg)=cmplx(0.d0,0.d0,kind=FFT1D_s)
         !FOURIER TRANSFORMS
         call SPCFFT(sv,3*nseg,-1,workfft)
         call SPCFFT(sf,3*nseg,-1,workfft)
         !COMPUTE Cvv AND Cvf FROM FOURIER TRANSFORMED TRAJECTORIES
         mCvv_type(:,itype,j,ibead)=mCvv_type(:, itype,j,ibead)
     &     +abs(sv(:))**2*mass(iglob)*fftfact
         Cvf_type(:,itype,ibead)=Cvf_type(:,itype,ibead)
     &     +real(sv(:)*conjg(sf(:)),8)*fftfact
       enddo
      enddo; enddo

      deallocate(sv,sf,workfft)

c
c     the master gets the sums of Cvv and Cvf contributions
c
      if(nproctot>1) then 
        call reduce_qtb_spectrum(mCvv_type,nbeads*ntype*3*nseg*3
     &                              ,MPI_COMM_WORLD,ranktot,.FALSE.)
        call reduce_qtb_spectrum(Cvf_type,nbeads*ntype*3*nseg
     &                             ,MPI_COMM_WORLD,ranktot,.FALSE.)
      endif

      if (ranktot.eq.0) then 
        ! Compute Delta FDR
        mCvv_sum(:,:,:)=0.d0
        do ibead=1,nbeads; do i=1,ntype; 
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
     &                     -Cvf_type(k,i,ibead)
          ! update intermediate averages for adaptation
            temp=mcvvsum_tmp - mCvv_mean(k,i,ibead)
            mCvv_mean(k,i,ibead) = mCvv_mean(k,i,ibead)
     &          + temp/compteur_bis%i
            temp=Cvf_type(k,i,ibead) - Cvf_mean(k,i,ibead)
            Cvf_mean(k,i,ibead) = Cvf_mean(k,i,ibead) 
     &          + temp/compteur_bis%i
            temp=dFDR_type(k,i,ibead) - dFDR_mean(k,i,ibead)
            dFDR_mean(k,i,ibead) = dFDR_mean(k,i,ibead) 
     &          + temp/compteur_bis%i
          enddo
        enddo; enddo

c
c      adapt gammar 
c
        if(adaptive_qtb) then
          call compteur_bis%update()  
          if (compteur_bis%reset_step) then 
            if(ranktot.eq.0 .and. QTB_verbose) then
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
     &          *gamma_type(k,i,ibead) - Cvf_type(k,i,ibead)
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
        mCvv_sum=mCvv_sum/real(n,8)
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
     &       "#Frequency   mCvv   Cvf/gammar  dFDT  gammar"
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
        call MPI_BCAST(gamma_type,ntype*nad*nbeads,MPI_REAL8
     &        ,0,MPI_COMM_WORLD,ierr)
c     recompute the noise kernel with the new gamma
        call adHnoise(.FALSE.)
      endif

c
c       COMPUTE corr_fact_qtb from deconvolution
c
      if(use_corr_fact_qtb .OR. noQTB 
     &   .OR. corr_fact_converged 
     &   .OR. (adqtb_start_classical .AND. (.not. adaptive_qtb))
     & ) return
      
      ! TODO: corr_fact for piqtb
      call compute_corr_fact(mCvv_sum)

      end subroutine QTB_spectra_adaptation_gamma

      subroutine write_frequency_type_header(unit)
        use ascii, only: int_to_str
        use beads, only: nbeads
        implicit none
        integer, intent(in) :: unit
        integer :: i, ibead

        write(unit,'(A)',advance='no') '#Frequency' 
        do ibead=1,nbeads; do i=1,ntype
          write(unit,'(6x,A)',advance='no') 
     &       "Type_"//trim(adqtb_name(i))
     &        //"_"//int_to_str(adqtb_type_map(i))
        enddo; enddo
        write(unit,*)

      end subroutine write_frequency_type_header

      subroutine reassignqtb(nloc_save,max_data_recv
     &      ,nneig_recep, n_data_recv, iglob_recv, pneig_recep
     &      ,nneig_send , n_data_send, iglob_send, pneig_send 
     &      ,max_atoms_recv, max_atoms_send )
        use atoms
        use mpi
        use domdec, only: rank,nproc,COMM_TINKER
        use spectra
        use beads
        implicit none
        integer, intent(in) :: nloc_save, max_data_recv
        integer,intent(in) :: nneig_recep, n_data_recv(nneig_recep)
     &                     , iglob_recv(max_data_recv,nneig_recep)
     &                     , pneig_recep(:)
        integer,intent(in) :: nneig_send, n_data_send(0:nneig_send)
     &                     , iglob_send(1:nloc_save,0:nneig_send)
     &                     , pneig_send(:)
        integer, intent(in) :: max_atoms_recv, max_atoms_send
        real, allocatable ::data_send(:,:,:,:,:),data_recv(:,:,:,:,:)
        integer :: nseg_comm,i,j,ineighbor,ierr,iseg,iglob
        integer :: req_data_send(nneig_send),req_data_recv(nneig_recep)
        integer :: pn,ii,jj,kk,ibead
c        integer :: n_data_s(nneig_send), n_data_r(nneig_recep)
c        integer :: iglob_s(max_atoms_send,nneig_send)
c        integer :: iglob_r(max_atoms_recv,nneig_recep)
c        integer :: max_atoms_s, max_atoms_r
        integer :: nneig_s,nneig_r
c        logical :: has_data_send(max_atoms_send,nneig_send)
c        logical :: has_data_recv(max_atoms_recv,nneig_recep)
c        integer :: n_ask_send(nneig_send), n_ask_recv(nneig_recep)
        logical :: was_allocated

        nseg_comm=3*nseg
        if (register_spectra) nseg_comm=4*nseg
        !nseg_comm=nseg
c        req_data_send (:) = MPI_REQUEST_NULL
c        req_data_recv (:) = MPI_REQUEST_NULL

c        !ASK IF NEIGHBORS HAVE DATA
c        n_ask_send(:)=0
c        ii=0
c        do ineighbor=1,nneig_send
c          pn=pneig_send(ineighbor)
c          do i=1,n_data_send(ineighbor)
c            iglob=iglob_send(i,ineighbor)
c            ! IF I THINK THAT pn DOES NOT HAVE DATA, IT MEANS WE NEVER COMMUNICATED (THIS ATOMS) DURING SEGMENT
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
            iglob=iglob_recv(i,ineighbor)  
            was_allocated=qtbdata_allocated(iglob)       
            if(.not.qtbdata_allocated(iglob)) then
              !allocate qtbdata for receiving atom for which I don't have data yet
c              n_data_r(ineighbor)=n_data_r(ineighbor)+1
c              iglob_r(n_data_r(ineighbor),ineighbor)=iglob
              call allocate_qtbdata(iglob) 
            endif
c            ! IF I THINK THAT pn DOES NOT HAVE DATA, IT MEANS WE NEVER COMMUNICATED (THIS ATOMS) DURING SEGMENT
c            ! => I SHOULD TELL THEM IF I HAVE DATA OR NOT
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
c
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
c            iglob=iglob_send(i,ineighbor)
c            !if(.not.qtbdata(iglob)%proc_has_data(pn+1)) then
c            !  ii=ii+1
c            !  if(.not. has_data_send(ii,ineighbor)) then
c            !    !send data if neigh does not have data yet
c            !    n_data_s(ineighbor)=n_data_s(ineighbor)+1
c            !    iglob_s(n_data_s(ineighbor),ineighbor)=iglob
c            !  endif              
c            !endif
c            qtbdata(iglob)%proc_has_data(pn+1)=.TRUE.
c            !mark atom as not local
c            qtbdata(iglob)%currently_loc=.FALSE.
c          enddo
c          if(n_data_s(ineighbor)>max_atoms_s)
c     &        max_atoms_s=n_data_s(ineighbor)
c        enddo

        req_data_send (:) = MPI_REQUEST_NULL
        req_data_recv (:) = MPI_REQUEST_NULL

        allocate(data_send(3,nseg_comm,nbeads
     &               ,max_atoms_send,nneig_send))
        allocate(data_recv(3,nseg_comm,nbeads
     &               ,max_atoms_recv,nneig_recep))
        do ineighbor=1,nneig_send 
          do i=1,n_data_send(ineighbor)
            iglob=iglob_send(i,ineighbor)
            data_send(:,1:nseg,:,i,ineighbor)
     &         =qtbdata(iglob)%fad(:,:,:)
            data_send(:,nseg+1:3*nseg,:,i,ineighbor)
     &         =qtbdata(iglob)%noise(:,:,:)
            if (register_spectra)
     &        data_send(:,3*nseg+1:4*nseg,:,i,ineighbor)
     &         =qtbdata(iglob)%vad(:,:,:)
          enddo
        enddo

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

        ! Wait all communications
        call MPI_Waitall(nneig_r,req_data_recv
     &                  ,MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(nneig_s ,req_data_send 
     &                  ,MPI_STATUSES_IGNORE, ierr)


        !put data received from buffer to qtbdata
        do ineighbor=1,nneig_recep
          do i=1,n_data_recv(ineighbor)
            iglob=iglob_recv(i,ineighbor)
            qtbdata(iglob)%fad(:,:,:)
     &        =data_recv(:,1:nseg,:,i,ineighbor)
            qtbdata(iglob)%noise(:,:,:)
     &        =data_recv(:,nseg+1:3*nseg,:,i,ineighbor)
            if (register_spectra)
     &       qtbdata(iglob)%vad(:,:,:)
     &        =data_recv(:,3*nseg+1:4*nseg,:,i,ineighbor)
          enddo
        enddo

        do ineighbor=1,nneig_send 
          do i=1,n_data_send(ineighbor)
            iglob=iglob_send(i,ineighbor)
            call deallocate_qtbdata(iglob)
          enddo
        enddo

      deallocate(data_send,data_recv)

      end subroutine reassignqtb

      subroutine reduce_qtb_spectrum(spectrum,size,comm,rank,allreduce)
      use mpi
      implicit none
      integer, intent(in) :: size,comm,rank
      logical, intent(in) :: allreduce
      real*8, intent(inout) :: spectrum(*)
      integer :: ierr

      if (allreduce) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,spectrum,
     &       size,MPI_REAL8,MPI_SUM,
     &       comm,ierr)
      else
        if (rank.eq.0) then
          call MPI_REDUCE(MPI_IN_PLACE,spectrum,
     &       size,MPI_REAL8,MPI_SUM,0,
     &       comm,ierr)
        else
          call MPI_REDUCE(spectrum,spectrum,
     &       size,MPI_REAL8,MPI_SUM,0,
     &       comm,ierr)
        end if
      endif

      end subroutine reduce_qtb_spectrum

      subroutine cleanup_qtbdata()
      use atoms
      use spectra
      use domdec
      implicit none
      integer :: i,j,k,iglob

      if(nproc==1) return

      if(register_spectra) then
        do i=1,nloc
          iglob=glob(i)
          qtbdata(iglob)%vad(:,:,:)=0.
        enddo
      endif

c      do i=1,n
c        if(qtbdata_allocated(i)) then 
c          if(qtbdata(i)%currently_loc) then
c            qtbdata(i)%proc_has_data(:)=.FALSE.
c            qtbdata(i)%proc_has_data(rank+1)=.TRUE.
c          else
c            call deallocate_qtbdata(i)
c          endif
c        endif
c      enddo
c      do i=1,nloc
c        iglob=glob(i)
c        if(.not. qtbdata(iglob)%currently_loc) then
c          write(0,*) rank,"atom",iglob,"should be local"
c          call fatal
c        endif
c      enddo

      end subroutine cleanup_qtbdata

       subroutine generate_corr_pot()
        use deconvolution
        use domdec
        use spectra
        use mpi
        use units
        use bath
        use langevin
        use beads
        use ascii, only: int_to_str
        implicit none
        real*8, allocatable :: omsym(:),s0sym(:)
        real*8, allocatable :: sout(:),srec(:)
        real*8 :: f_omega,w
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
     &         ,niter_deconvolution,.true.
     &         ,1.0d-10,gamma_friction(ibead)
     &         ,boltzmann*kelvin
     &         ,kernel_lorentz_pot,sout,srec,verbose=.FALSE.)
            do k=1,nad
              corr_pot_ratio(k,ibead)=(sout(nad+k-1)-s0sym(nad+k-1))
     &                          /s0sym(nad+k-1)
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
        call MPI_BCAST(corr_pot_ratio,nbeads*nad,MPI_REAL8
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
        real*8, intent(in) :: mCvv_in(:,:,:)
        real*8 :: mCvv(size(mCvv_in,1),size(mCvv_in,2))
        real*8, allocatable :: s0sym(:),omsym(:),sdec(:),srec(:)
        real*8, allocatable :: corr(:)
        real*8 :: rec_ratio,rec_ratio_min
        integer :: k,ierr,jend,j
        integer, save :: isame=0
        real*8, save :: corr_fact_prev(3)=1.0d0
        real*8, parameter :: thr=1.d-4

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
     &           ,niter_deconvolution,.false.,1.0d-10
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

        CALL MPI_BCAST(corr_fact_qtb,3,MPI_REAL8,0,COMM_TINKER,ierr)
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
        real*8, intent(in) :: dFDR(:,:)
        real*8, intent(in) :: mCvv(:,:),Cvf(:,:)
        logical,intent(in) :: qtb_therm
        real*8 :: a1,gamma_ref,norm,gamma_disp,normv(ntype)
        real*8, allocatable, save :: dfdr_m(:,:,:),dfdr_norm(:,:)
        real*8 :: b1,bias_correct1,smooth,dFDRs(nad,nbeads)
        integer, allocatable,save :: cnt_type(:)
        integer :: cnt
        integer :: i,j,l,ibead
        real*8 :: dfdr_norm_instant

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
        real*8, intent(in) :: dFDR(:,:)
        real*8, intent(in) :: mCvv(:,:),Cvf(:,:)
        logical,intent(in) :: qtb_therm
        real*8 :: b1,tau,sqrtm,w,b_avg,avg_correct,a1
        real*8 :: avg_correct_prev
        real*8 :: adapt_correct,tmp
        real*8 :: gamma_new(nad),disp,smooth
        real*8, allocatable, save :: gamma_curr(:,:,:)
        real*8, allocatable, save :: Cvf_m(:,:,:)
        real*8, allocatable, save :: mcvv_m(:,:,:)
        integer, allocatable, save :: nadapt(:)
        integer :: i,j,l,ibead
        integer, allocatable, save :: cnt_type(:)
        real*8 :: mCvvs(nad,nbeads), Cvfs(nad,nbeads)
        integer :: cnt
        real*8 :: tau_avg

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
        real*8, intent(inout) :: s(:)
        real*8, intent(in) :: dw
        integer :: i,j,jstart,jend
        real*8 :: w,w0,norm,b1
        real*8, allocatable :: s_in(:)

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

c-------------------------------------------------------------
      subroutine broadcast_piqtb_noise()
      use beads
      use mpi
      use beads
      use domdec
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

      DO iloc=1,ilocend-ilocbeg+1
        i=glob(ilocbeg+iloc-1)
        bufferpi_s(:,1:2*nseg,:,iloc)=qtbdata(i)%noise(:,:,:)
        bufferpi_s(:,2*nseg+1:3*nseg,:,iloc)=qtbdata(i)%fad(:,:,:)
      ENDDO

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

      call MPI_WAITALL(ii,reqs,MPI_STATUSES_IGNORE,ierr)

      DO iloc=1,nloc
        iproc=min(nproc_polymer,int((iloc-1)/nlocproc)+1)
        if(rank_polymer==iproc-1) CYCLE
        ilocbeg=ilocpi_beg(iproc)
        ilocend=ilocpi_end(iproc)
        i=glob(iloc)
        kk=iloc-ilocbeg +1
        qtbdata(i)%noise(:,:,:)=bufferpi_r(:,1:2*nseg,:,kk,iproc)
        qtbdata(i)%fad(:,:,:)=bufferpi_r(:,2*nseg+1:3*nseg,:,kk,iproc)
      ENDDO

      deallocate(bufferpi_s,bufferpi_r)

      end subroutine broadcast_piqtb_noise

      subroutine reduce_piqtb_vad()
      use mpi
      use beads
      use domdec
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

      nelt=nseg

      maxloc=nlocpis(nproc_polymer)
      nlocproc=int(nloc/nproc_polymer)

      allocate(buffer(3,nelt,nloc,nbeads))

        !PUT DATA INTO BUFFER ARRAYS

      ilocbeg=ilocpi_beg(rank_polymer+1)
      ilocend=ilocpi_end(rank_polymer+1)

      do k=1,nbeads; do i=1,nloc; do l=1,nseg; do j=1,3
        iglob=glob(i)
        buffer(j,l,i,k)=qtbdata(iglob)%vad(j,l,k)
      enddo; enddo; enddo; enddo

      call MPI_ALLREDUCE(MPI_IN_PLACE,buffer
     &       ,3*nelt*nbeads*nloc 
     &       ,MPI_REAL4,MPI_SUM
     &       ,COMM_POLYMER,ierr)

      do k=1,nbeads; do i=1,nloc; do l=1,nseg; do j=1,3
        iglob=glob(i)
        qtbdata(iglob)%vad(j,l,k)=buffer(j,l,i,k)
      enddo; enddo; enddo; enddo

      deallocate(buffer)

      end subroutine reduce_piqtb_vad

      end module qtb
