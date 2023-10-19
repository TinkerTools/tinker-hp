c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine baoarespabpi  --  BAOAB-respa Langevin PIMD step     ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "baoabrespapi" performs a single PIMD time step
c     via the BAOAB-respa recursion formula 
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c
      
      subroutine baoabrespapi(polymer,polymer_ctr,istep,dt)
      use bath
      use atoms
      use atmtyp
      use beads
      use domdec
      use inform
      use iounit
      use mdstuf
      use moldyn
      use timestat
      use mpi
      use spectra
      use deriv
      use units
      use bath
      use langevin
      use utilbaoabpi
      use cutoff
      use potent
      use virial
      use polar
      use moldyn
      use utilbaoabpi
      use commstuffpi
      use spectra
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
      integer, intent(in) :: istep
      real*8, intent(in) :: dt
      integer :: i,k
      real(8),save::timebegin=0d0
      real(8),save::timefull=0d0
      real*8 :: dt2,sqrtnu,dta
      integer :: iloc
      integer :: ibeadbeg
      logical :: only_long
      real*8 :: dip(3,nalt),dipind(3,nalt)
      interface
        subroutine baoabrespapi_fast(polymer
     &     ,istep,dt,nsteps,epot,vir_,dip_)
          import :: POLYMER_COMM_TYPE 
          TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
          integer, intent(in) :: istep,nsteps
          real*8, intent(in) :: dt
          real*8, intent(inout) :: epot,vir_(3,3)
          real*8 :: dip_(3,nsteps)
        end subroutine baoabrespapi_fast
      end interface

      dt2  = 0.5d0 * dt
      dta   = dt / dshort
      sqrtnu=sqrt(real(nbeads,8))        

      if(PITIMER) call stopwatchpi(timebegin,.true.,.true.)

      !! B SLOW !!
      call apply_B_PI(polymer,dt2)
      !! LOAD PREVIOUS FAST FORCES ON EIGFORCES !!
      call set_eigforces_pi(polymer,polymer%forces)

      if(isobaric) then
        call plangevinpi(polymer,dt2,-1,'A')
        call plangevinpi(polymer,dt,-1,'O')
      endif

      !! INTERNAL BAOAB STEPS FOR FAST EVOLVING FORCES !!
      call baoabrespapi_fast(polymer
     &   ,istep,dta,nalt,epotpi,virpi,dip)

      if(isobaric) call plangevinpi(polymer,dt2,-1,'A')
      if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

      call full_gradient_prepare_pi(polymer,polymer_ctr,istep,.true.)
      if(PITIMER) call stopwatchpi(timecom,.true.,.false.)

      if(centroid_recip) then
        !!! COMPUTE FORCES ON CENTROID!!!  
        if(allocated(derivs_centroid)) deallocate(derivs_centroid)
        allocate(derivs_centroid(3,nbloc))
        call load_bead(polymer,0)
        only_long = centroid_longrange .and.
     &        (.not. (contract.and.nbeads_ctr==1) )
        call compute_gradslow_centroid(eslow,derivs_centroid
     &              ,vir_ctr,only_long, polar_allbeads)
        if(rank_polymer == 0) then
          epotpi = epotpi + eslow
          virpi(:,:) = virpi(:,:) + vir_ctr(:,:)
        endif
        if(PITIMER) call stopwatchpi(timecontract,.false.,.false.)
        call commforcesrespa1(derivs_centroid,2)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
        DO k=1,nbeads; DO iloc=1,nloc
          i=glob(iloc)
          polymer%forces_slow(:,i,k) = 
     &         -derivs_centroid(:,iloc)
        ENDDO ; ENDDO

      else

        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
        DO k=1,nbeads; DO iloc=1,nloc
          i=glob(iloc)
          polymer%forces_slow(:,i,k) = 0.d0
        ENDDO ; ENDDO    

      endif

      if(.not. contract) then
          !! COMPUTE SLOW FORCES ON FULL POLYMER !! 
        if(allocated(derivs)) deallocate(derivs)
        allocate(derivs(3,nbloc,nbeadsloc))
        call compute_grad_beads(eslow,derivs,vir_ctr,polymer
     &         ,.not.centroid_longrange,.TRUE.,.FALSE. ! LONG, INT, SHORT
     &        ,polar_allbeads)
        epotpi = epotpi + eslow
        virpi(:,:) = virpi(:,:) + vir_ctr(:,:)
        if(PITIMER) call stopwatchpi(timegradpi,.false.,.false.)

        if(.not.centroid_recip)
     &       call commforcesrecpi(derivs,derivsRec,nbeadsloc)
        call commforcesblocpi(derivs,nbeadsloc,.FALSE.)

        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
        DO k=1,nbeadsloc; DO iloc=1,nloc
          i=glob(iloc)
          polymer%forces_slow(:,i,k+ibeadbeg-1) =
     &        polymer%forces_slow(:,i,k+ibeadbeg-1)       
     &        -derivs(:,iloc,k)
        ENDDO ; ENDDO      

        call comm_for_normal_modes(polymer,polymer%forces_slow )
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
      endif
        

      if(contract .and. nbeads_ctr>1) then
        !! COMPUTE SLOW FORCES ON CONTRACTED POLYMER !! 
        if(allocated(derivs_ctr)) deallocate(derivs_ctr)
        allocate(derivs_ctr(3,nbloc,nbeadsloc_ctr))
        call compute_grad_beads(eslow,derivs_ctr,vir_ctr,polymer_ctr
     &      ,.not.centroid_longrange,.TRUE.,.FALSE. ! LONG, INT, SHORT
     &     ,polar_allbeads)
        epotpi = epotpi + eslow
        virpi(:,:) = virpi(:,:) + vir_ctr(:,:)
        if(PITIMER) call stopwatchpi(timecontract,.false.,.false.)
        
        if(.not.centroid_recip) 
     &       call commforcesrecpi(derivs_ctr,derivsRec,nbeadsloc_ctr)
        call commforcesblocpi(derivs_ctr,nbeadsloc_ctr,.FALSE.)
        call commforces_ctr(polymer_ctr,derivs_ctr)
        call project_forces_ctr(polymer,polymer_ctr,.TRUE.)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
      endif

      call set_eigforces_pi(polymer,polymer%forces_slow)
      call apply_B_PI(polymer,dt2)
      if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

      call kinetic_pi(polymer)
      call reduce_observables_pi(polymer)
      if(ir) then
        dip(:,nalt)=dippi
        dipind(:,nalt)=dipindpi
        call save_dipole_respa(dip,dipind)
      endif
      if(isobaric) call plangevinpi(polymer,dt,-1,'B')

      call mdstatpi(istep,dt)
      call mdsavebeads (istep,dt,polymer)

      if(PITIMER) call stopwatchpi(timeobs,.false.,.false.)

      if(PITIMER) then
        timefull=timefull+mpi_wtime()-timebegin

        if (mod(istep,iprint).eq.0) then
          if(verbose .and. ranktot.eq.0) then
            write(*,*)
            write(*,*) "### PIMD TIMERS ###"
c                write(*,'(A,f10.7,A,f7.3,A)') '      Time prepare: '     
c     &      , timereass/iprint , ' ( ',100*timereass/timefull,' %)'
c                write(*,'(A,f10.7,A,f7.3,A)') ' Time normal modes: '
c     &           , timepush/iprint, ' ( ',100*timepush/timefull,' %)'
                write(*,'(A,f10.7,A,f7.3,A)') 'Time grad contract: '
     &  , timecontract/iprint , ' ( ',100*timecontract/timefull,' %)'
                write(*,'(A,f10.7,A,f7.3,A)') '    Time grad slow: '
     &     , timegradpi/iprint, ' ( ',100*timegradpi/timefull,' %)'
                write(*,'(A,f10.7,A,f7.3,A)') '    Time grad fast: '
     &            , timegradfastpi/iprint
     &                      ,' ( ',100*timegradfastpi/timefull,' %)'
                write(*,'(A,f10.7,A,f7.3,A)') '  Time observables: '
     &           , timeobs/iprint, ' ( ',100*timeobs/timefull,' %)'
                write(*,'(A,f10.7,A,f7.3,A)') '        Time comms: '
     &           ,timecom/iprint, ' ( ',100*timecom/timefull,' %)'
                write(*,'(A,f10.7,A,f7.3,A)') '        Time baoab: '
     &           , timeaoa/iprint, ' ( ',100*timeaoa/timefull,' %)'
                write(*,'(A,f10.7)')          '         Time step: '
     &           , timefull/iprint
                write(*,*) "###################"
                write(*,*)
          endif
          timefull=0.d0
          timepush=0d0
          timereass=0d0
          timeinit=0d0
          timecom=0.d0
          timeaoa=0.d0
          timegradpi=0.d0
          timegradfastpi=0.d0
          timeobs=0.d0
          timecontract=0.d0
          timesavenl=0.d0
        endif
      endif

      end subroutine baoabrespapi

      subroutine baoabrespapi_fast(polymer,
     &   istep,dta,nsteps,epot,vir_,dip)
      use bath
      use atoms
      use atmtyp
      use beads
      use domdec
      use inform
      use iounit
      use mdstuf
      use moldyn
      use timestat
      use mpi
      use spectra
      use deriv
      use units
      use bath
      use langevin
      use utilbaoabpi
      use cutoff
      use potent
      use virial
      use polar
      use moldyn
      use utilbaoabpi
      use commstuffpi
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer, intent(in) :: istep,nsteps
      real*8, intent(in) :: dta
      real*8, intent(inout) :: epot,vir_(3,3)
      real*8, intent(inout) :: dip(3,nsteps)
      integer :: i,k,ialt
      real*8 :: dta2
      integer :: iloc
      real*8 :: rsteps   
      integer :: ibeadbeg

      dta2 = 0.5d0*dta
      rsteps=real(nsteps,8)
      !! ALLOCATE NOISE VECTOR !!
      vir_(:,:) = 0.d0
      if(allocated(derivs)) deallocate(derivs)
      allocate(derivs(3,nbloc,nbeadsloc))

      !! START FAST EVOLVING LOOP !!
      fastloop: do ialt = 1, nsteps
        !! APPLY FIRST PART OF BAOAB INTEGRATION !! 
        call apply_B_PI(polymer,dta2)
        call apply_AOA_PI(polymer,dta)
        if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

        call fast_gradient_prepare_pi(polymer,istep)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)

        !! GRADFAST
        call compute_grad_beads(epot,derivs,vir_ctr,polymer
     &       ,.FALSE.,.FALSE.,.TRUE. ! LONG, INT, SHORT
     &       ,.FALSE.)
        vir_(:,:) = vir_(:,:) + vir_ctr(:,:)/rsteps
        if(ir) then
          dip(:,ialt)=dippi
          !dipind(:,ialt)=dipindpi
        endif
        if(PITIMER) call stopwatchpi(timegradfastpi,.false.,.false.)

        !!! COMMUNICATE BEADS FORCES (NEWTON LAW)!!!
        call commforcesblocpi(derivs,nbeadsloc,.true.)

        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
        DO k=1,nbeadsloc; DO iloc=1,nloc
          i=glob(iloc)
          polymer%forces(:,i,k+ibeadbeg-1) = -derivs(:,iloc,k)
        ENDDO ; ENDDO   

        call comm_for_normal_modes(polymer,polymer%forces)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)

        call set_eigforces_pi(polymer,polymer%forces)
        !! PROPAGATE VELOCITIES USING FAST GRADIENT !!
        call apply_B_PI(polymer,dta2)

      enddo fastloop

      call kinetic_pi_fast(polymer)

      end subroutine baoabrespapi_fast

