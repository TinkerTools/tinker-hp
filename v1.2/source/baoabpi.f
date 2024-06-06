c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine baoabpi  --  BAOAB Langevin PIMD step                ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "baoabpi" performs a single PIMD time step
c     via the BAOAB recursion formula 
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c
      subroutine baoabpi(polymer,polymer_ctr,istep,dt)
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
      use cutoff
      use potent
      use virial
      use polar
      use moldyn
      use utilbaoabpi
      use commstuffpi
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer_ctr
      integer, intent(in) :: istep
      real*8, intent(in) :: dt
      integer :: i,j,k
      real(8),save::timebegin=0d0
      real(8),save::timefull=0d0
      real*8 :: dt2,sqrtnu
      integer :: iloc
      integer :: ibeadbeg
      logical :: only_long
c
      if (deb_Path) write(iout,*), 'baoabpi '
c

      dt2=0.5d0*dt
      sqrtnu=sqrt(real(nbeads,8))
      

      if(PITIMER) call stopwatchpi(timebegin,.true.,.true.)

      !! APPLY FIRST PART OF BAOAB INTEGRATION !! 
      call apply_B_PI(polymer,dt2)
      if(isobaric) then
        call plangevinpi(polymer,dt2,-1,'A')
        call plangevinpi(polymer,dt,-1,'O')
      endif
      call apply_AOA_PI(polymer,dt)
      if(isobaric) call plangevinpi(polymer,dt2,-1,'A')
      if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)
      
      call full_gradient_prepare_pi(polymer,polymer_ctr,istep,.false.)
      if(PITIMER) call stopwatchpi(timecom,.false.,.false.)

      !!! COMPUTE FORCES ON CENTROID!!!  
      if(centroid_recip) then
        if(allocated(derivs_centroid)) deallocate(derivs_centroid)
        allocate(derivs_centroid(3,nbloc))
        call load_bead(polymer,0)
        only_long = centroid_longrange .and.
     &      (.not. (contract.and.nbeads_ctr==1) )
        call compute_gradslow_centroid(eslow,derivs_centroid
     &            ,vir_ctr,only_long, polar_allbeads)
        if(PITIMER) call stopwatchpi(timecontract,.false.,.false.)
        call commforcesrespa1(derivs_centroid,2)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
      endif

      !!! COMPUTE FORCES ON BEADS!!! 
      if(allocated(derivs)) deallocate(derivs)
      allocate(derivs(3,nbloc,nbeadsloc))
      call compute_grad_beads(epotpi,derivs,virpi,polymer
     &         ,.not. (centroid_longrange .or. contract)  !LONG 
     &         ,.not. contract ! INT
     &         , .TRUE.     !SHORT 
     &         ,polar_allbeads)

      if(PITIMER) call stopwatchpi(timegradpi,.false.,.false.)
      if(.not.(centroid_recip .or. contract))
     &    call commforcesrecpi(derivs,derivsRec,nbeadsloc)
      call commforcesblocpi(derivs,nbeadsloc,.FALSE.)  

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      DO k=1,nbeadsloc; DO iloc=1,nloc
        i=glob(iloc)
        polymer%forces(:,i,k+ibeadbeg-1) = -derivs(:,iloc,k)
      ENDDO ; ENDDO

      !! ADD CENTROID CONTRIBUTION TO ENERGY AND FORCES !!
      if(centroid_recip) then 
        if(rank_polymer==0) then
          epotpi = epotpi + eslow
          virpi(:,:) = virpi(:,:) + vir_ctr(:,:)
        endif  
        DO k=1,nbeadsloc; DO iloc=1,nloc
          i=glob(iloc)
          polymer%forces(:,i,k+ibeadbeg-1) = 
     &       polymer%forces(:,i,k+ibeadbeg-1)
     &                   -derivs_centroid(:,iloc) 
        ENDDO ; ENDDO
      endif


      call comm_for_normal_modes(polymer,polymer%forces)
      if(PITIMER) call stopwatchpi(timecom,.true.,.false.)

      !!! COMPUTE FORCES ON  CONTRACTED BEADS!!!  
      if(contract .and. nbeads_ctr>1) then
        if(allocated(derivs_ctr)) deallocate(derivs_ctr)
        allocate(derivs_ctr(3,nbloc,nbeadsloc_ctr))
        call compute_grad_beads(eslow,derivs_ctr,vir_ctr,polymer_ctr
     &       ,.not.centroid_longrange, .TRUE., .FALSE.  ![LONG, INT  , SHORT]
     &       ,polar_allbeads)
        epotpi = epotpi + eslow
        virpi(:,:) = virpi(:,:) + vir_ctr(:,:)
      if(PITIMER) call stopwatchpi(timecontract,.false.,.false.)
      
        if(.not.centroid_recip)
     &     call commforcesrecpi(derivs_ctr,derivsRec,nbeadsloc_ctr)
        call commforcesblocpi(derivs_ctr,nbeadsloc_ctr,.FALSE.)
        call commforces_ctr(polymer_ctr,derivs_ctr)
        call project_forces_ctr(polymer,polymer_ctr)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
      endif

      call set_eigforces_pi(polymer,polymer%forces)
      call centroid_force_corrections(polymer)
      call apply_B_PI(polymer,dt2)
      if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

      call kinetic_pi(polymer)
      call reduce_observables_pi(polymer)
      if(ir) call save_dipole_traj(dippi,dipindpi) 
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
c              write(*,'(A,f10.7,A,f7.3,A)') '      Time prepare: '     
c     &    , timereass/iprint , ' ( ',100*timereass/timefull,' %)'
c              write(*,'(A,f10.7,A,f7.3,A)') ' Time normal modes: '
c     &         , timepush/iprint, ' ( ',100*timepush/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') 'Time grad contract: '
     & ,timecontract/iprint , ' ( ',100*timecontract/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '   Time grad beads: '
     &         , timegradpi/iprint, ' ( ',100*timegradpi/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '  Time observables: '
     &         , timeobs/iprint, ' ( ',100*timeobs/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '        Time comms: '
     &         ,timecom/iprint, ' ( ',100*timecom/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '        Time baoab: '
     &         , timeaoa/iprint, ' ( ',100*timeaoa/timefull,' %)'
              write(*,'(A,f10.7)')          '         Time step: '
     &         , timefull/iprint
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
          timeobs=0.d0
          timecontract=0.d0
          timesavenl=0.d0
        endif
      endif

      end subroutine baoabpi

