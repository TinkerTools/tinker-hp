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
#include "tinker_precision.h"
      
      subroutine baoabrespapi(istep,dt)
      use bath
      use atomsMirror
      use atmtyp
      use beads
      use domdec
      use energi
      use inform
      use iounit
      use mdstuf
      use moldyn
      use timestat
      use mpi
      use spectra
      use deriv
      use utilgpu
      use random_mod
      use units
      use bath
      use langevin
      use utilbaoabpi
      use cutoff
      use utils
      use potent
      use virial
      use polar
      use moldyn
      use utilbaoabpi
      use commstuffpi
      use spectra
      implicit none
      integer, intent(in) :: istep
      real(r_p), intent(in) :: dt
      integer :: ibead,i,j,k,iglob,ierr,kk
      real(8),save::time0=0d0,time1=0d0,time00=0d0
      real(8),save::time01=0d0, timebegin=0d0
      real(8),save::timefull=0d0
      logical :: add_eslow
      real(r_p) :: dt2,sqrtnu,factor1,factor2,dta,dta2
      real(r_p) :: scale,eigx0,eigv0,a1,a2
      integer :: nnoise,iloc,iqtb
      real(r_p) :: gammak,nuratio     
      integer :: ilocbeg,ilocend,ibeadbeg
      logical :: only_long
      real(r_p), allocatable, save :: dip(:,:),dipind(:,:)
      interface
        subroutine baoabrespapi_fast(istep,dta,nsteps,epot,vir_
     &   ,dip_,dipind_)
        implicit none
        integer, intent(in) :: istep,nsteps
        real(r_p), intent(in) :: dta
        real(r_p), intent(inout) :: epot,vir_(3,3)
        real(r_p), intent(inout) :: dip_(:,:),dipind_(:,:)
        end subroutine baoabrespapi_fast
      end interface

      if(istep==1) then
        allocate(dip(3,nalt),dipind(3,nalt))
!$acc enter data create(eslow,vir_ctr,dip,dipind) async
        stepint=nalt
        stepfast=nalt2
      endif

      dt2  = 0.5_re_p * dt
      dta   = dt / dshort

      sqrtnu=sqrt(real(nbeads,r_p))        

      if(PITIMER) call stopwatchpi(timebegin,.true.,.true.)

      !! B SLOW !!
      call apply_B_PI(polymer,dt2)
      !! LOAD PREVIOUS FAST FORCES ON EIGFORCES !!
      call set_eigforces_pi(polymer,polymer%forces)

      if(isobaric) then
        call plangevinpi(polymer,dt2,istep,'A')
        call plangevinpi(polymer,dt,istep,'O')
      endif

      !! INTERNAL BAOAB STEPS FOR FAST EVOLVING FORCES !!
      call baoabrespapi_fast(istep,dta,nalt,epotpi,virpi
     &  ,dip,dipind )

      if(isobaric) call plangevinpi(polymer,dt2,istep,'A')
      if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

      call full_gradient_prepare_pi(polymer,polymer_ctr,istep,.true.)
      if(PITIMER) call stopwatchpi(timecom,.true.,.false.)

      if(centroid_recip) then
        !!! COMPUTE FORCES ON CENTROID!!!  
        call prmem_requestm(derivs_centroid,3,nbloc,async=.TRUE.)
        call load_bead(polymer,0)
        only_long = centroid_longrange .and.
     &        (.not. (contract.and.nbeads_ctr==1) )
        call compute_gradslow_centroid(eslow,derivs_centroid
     &              ,vir_ctr,only_long, polar_allbeads)
        if(rank_polymer == 0) then
!$acc serial async present(epotpi,eslow,virpi,vir_ctr)
          epotpi = epotpi + eslow
          do j=1,3; do i=1,3
            virpi(i,j) = virpi(i,j) + vir_ctr(i,j)
          enddo; enddo
!$acc end serial
        endif
        if(PITIMER) call stopwatchpi(timecontract,.false.,.false.)
        call commforcesrespa1(derivs_centroid,2)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
!$acc parallel loop collapse(3) async default(present)
        DO k=1,nbeads; DO iloc=1,nloc ; DO j=1,3
          i=glob(iloc)
          polymer%forces_slow(j,i,k) = 
     &         -derivs_centroid(j,iloc)
        ENDDO ; ENDDO ; ENDDO

        if (deb_Energy) call info_energy(rank)
        if (deb_Force)  then
          call info_forces(cNBond)
          if(ranktot==0) write(*,*) 'were Forces for centroid:'
        endif
        if (deb_Atom)   call info_minmax_pva

      else

        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
!$acc parallel loop collapse(3) async default(present)
        DO k=1,nbeads; DO iloc=1,nloc ; DO j=1,3
          i=glob(iloc)
          polymer%forces_slow(j,i,k) = 0._re_p
        ENDDO ; ENDDO ; ENDDO    

      endif

      if(.not. contract) then
          !! COMPUTE SLOW FORCES ON FULL POLYMER !! 
        call prmem_requestm(derivs,3,nbloc,nbeadsloc,async=.true.)
        call compute_grad_beads(eslow,derivs,vir_ctr,polymer
     &         ,.not.centroid_longrange,.TRUE.,.FALSE. ! LONG, INT, SHORT
     &        ,polar_allbeads)
!$acc serial async present(epotpi,eslow,virpi,vir_ctr)
        epotpi = epotpi + eslow
        do j=1,3; do i=1,3
          virpi(i,j) = virpi(i,j) + vir_ctr(i,j)
        enddo; enddo
!$acc end serial
        if(PITIMER) call stopwatchpi(timegrad,.false.,.false.)

        if(.not.centroid_recip)
     &       call commforcesrecpi(derivs,derivsRec,nbeadsloc)
        call commforcesblocpi(derivs,nbeadsloc,.FALSE.)

        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
!$acc parallel loop collapse(3) async default(present)
        DO k=1,nbeadsloc; DO iloc=1,nloc ; DO j=1,3
          i=glob(iloc)
          polymer%forces_slow(j,i,k+ibeadbeg-1) =
     &        polymer%forces_slow(j,i,k+ibeadbeg-1)       
     &        -derivs(j,iloc,k)
        ENDDO ; ENDDO ; ENDDO      

        call comm_for_normal_modes(polymer,polymer%forces_slow )
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
      endif
        

      if(contract .and. nbeads_ctr>1) then
        !! COMPUTE SLOW FORCES ON CONTRACTED POLYMER !! 
        call prmem_requestm(derivs_ctr,3,nbloc
     &    ,nbeadsloc_ctr,async=.true.)
        call compute_grad_beads(eslow,derivs_ctr,vir_ctr,polymer_ctr
     &      ,.not.centroid_longrange,.TRUE.,.FALSE. ! LONG, INT, SHORT
     &     ,polar_allbeads)
!$acc serial async present(epotpi,eslow,virpi,vir_ctr)
        epotpi = epotpi + eslow
        do j=1,3; do i=1,3
          virpi(i,j) = virpi(i,j) + vir_ctr(i,j)
        enddo; enddo
!$acc end serial
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

      if (calc_e.or.use_virial) then
        call kinetic_pi(polymer)
        call reduce_observables_pi(polymer)
      endif
      
      if(ir) then
!$acc serial async present(dip,dippi,dipind,dipindpi)
        dip(:,nalt)=dippi
        dipind(:,nalt)=dipindpi
!$acc end serial
!$acc update host(dip,dipind) async
!$acc wait
        call save_dipole_respa(dip,dipind)
      endif
      if(isobaric) call plangevinpi(polymer,dt,istep,'B')

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
     &           , timegrad/iprint, ' ( ',100*timegrad/timefull,' %)'
                write(*,'(A,f10.7,A,f7.3,A)') '    Time grad fast: '
     &      , timegradfast/iprint, ' ( ',100*timegradfast/timefull,' %)'
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
            timegrad=0.d0
            timegradfast=0.d0
            timeobs=0.d0
            timecontract=0.d0
            timesavenl=0.d0
          endif
        endif

      end subroutine baoabrespapi

      subroutine baoabrespapi_fast(istep,dta,nsteps,epot,vir_
     &   ,dip,dipind)
      use bath
      use atomsMirror
      use atmtyp
      use beads
      use domdec
      use energi
      use inform
      use iounit
      use mdstuf
      use moldyn
      use timestat
      use mpi
      use spectra
      use deriv
      use utilgpu
      use random_mod
      use units
      use bath
      use langevin
      use utilbaoabpi
      use cutoff
      use utils
      use potent
      use virial
      use polar
      use moldyn
      use utilbaoabpi
      use commstuffpi
      implicit none
      integer, intent(in) :: istep,nsteps
      real(r_p), intent(in) :: dta
      real(r_p), intent(inout) :: epot,vir_(3,3)
      real(r_p), intent(inout) :: dip(:,:),dipind(:,:)
      integer :: ibead,i,j,k,iglob,ierr,kk,ialt
      real(8),save::time0=0d0,time1=0d0,time00=0d0
      real(8),save::time01=0d0, timebegin=0d0
      real(8),save::timefull=0d0
      logical :: add_eslow
      real(r_p) :: dt2,sqrtnu,factor1,factor2,dta2
      real(r_p) :: scale,eigx0,eigv0,a1,a2
      integer :: nnoise,iloc,iqtb
      real(r_p) :: gammak,nuratio,rsteps   
      integer :: ilocbeg,ilocend,ibeadbeg
      logical :: only_long

      dta2 = 0.5d0*dta
      rsteps=real(nsteps,r_p)
      !! ALLOCATE NOISE VECTOR !!
      call prmem_request(noise,3*nlocpi*nbeads,async=.true.)
!$acc serial async present(vir_)
      do j=1,3; do i=1,3
        vir_(i,j) = 0.0_re_p
      enddo; enddo
!$acc end serial

      !! START FAST EVOLVING LOOP !!
      fastloop: do ialt = 1, nsteps
        if(deb_path) then
!$acc wait
          write(*,*) '  >> fastloop',ialt
        endif

        !! APPLY FIRST PART OF BAOAB INTEGRATION !! 
        call apply_B_PI(polymer,dta2)
        call apply_AOA_PI(polymer,dta)
        if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

        call fast_gradient_prepare_pi(polymer,istep)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)

        !! GRADFAST
        call prmem_requestm(derivs,3,nbloc,nbeadsloc,async=.true.)
        call compute_grad_beads(epot,derivs,vir_ctr,polymer
     &       ,.FALSE.,.FALSE.,.TRUE. ! LONG, INT, SHORT
     &       ,.FALSE.)
!$acc serial async present(vir_,vir_ctr)
        do j=1,3; do i=1,3
          vir_(i,j) = vir_(i,j) + vir_ctr(i,j)/rsteps
        enddo; enddo
!$acc end serial
        if(ir) then
!$acc serial async present(dip,dippi)
          dip(:,ialt)=dippi(:)
          !dipind(:,ialt)=dipindpi
!$acc end serial
        endif
        if(PITIMER) call stopwatchpi(timegradfast,.false.,.false.)

        !!! COMMUNICATE BEADS FORCES (NEWTON LAW)!!!
        call commforcesblocpi(derivs,nbeadsloc,.true.)

        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
!$acc parallel loop collapse(3) async default(present)
        DO k=1,nbeadsloc; DO iloc=1,nloc ; DO j=1,3
          i=glob(iloc)
          polymer%forces(j,i,k+ibeadbeg-1) = -derivs(j,iloc,k)
        ENDDO ; ENDDO ; ENDDO   

        call comm_for_normal_modes(polymer,polymer%forces)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)

        call set_eigforces_pi(polymer,polymer%forces)
        !! PROPAGATE VELOCITIES USING FAST GRADIENT !!
        call apply_B_PI(polymer,dta2)

      enddo fastloop

      if (calc_e .or. use_virial) then
        call kinetic_pi_fast(polymer)
      endif

      end subroutine baoabrespapi_fast

