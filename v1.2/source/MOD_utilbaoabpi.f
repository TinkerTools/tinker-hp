c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################################
c     ##                                                                                   ##
c     ##  module baoabopi  --  BAOAB Path Integral Langevin molecular dynamics step    ##
c     ##                                                                                   ##
c     #######################################################################################
c
c
c     "baoabpi" performs a single molecular dynamics time step
c     via the BAOAB recursion formula
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c
      module utilbaoabpi        
      implicit none

      real*8,allocatable, save::derivs(:,:,:)
      real*8,allocatable, save::derivs_ctr(:,:,:)
      real*8,allocatable, save::derivs_centroid(:,:)
      real*8,allocatable, save::derivsRec(:,:,:)
      real*8 :: eslow,vir_ctr(3,3)

      real(8), save :: timepush=0.d0, timeinit=0.d0
      real(8), save :: timereass=0.d0, timegradpi=0.d0
      real(8), save :: timesavenl=0.d0, timecom=0.d0
      real(8), save :: timecontract=0d0,timegradfastpi=0.d0
      real(8), save :: timeobs=0d0,timeaoa=0d0

      contains

      subroutine apply_B_PI(polymer,dt)
      use beads
      use units, only: convert
      use domdec
      use atmtyp, only: mass
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real*8, intent(in) :: dt
      integer :: nu,iloc, i,k,ilocbeg

      nu=polymer%nbeads
      ilocbeg = ilocpi_beg(rank_polymer+1)
      do k=1,nu; do iloc=1,nlocpi
        i=glob(iloc+ilocbeg-1)
        polymer%eigvel(:,i,k)=polymer%eigvel(:,i,k)
     &      +dt*convert*polymer%eigforces(:,i,k)/mass(i)
      enddo; enddo

      end subroutine apply_B_PI

      subroutine apply_AOA_PI(polymer,dt)
        use atoms
        use atmtyp
        use beads
        use domdec
        use inform
        use bath
        use qtb
        use langevin, only: gamma_friction
        use moldyn
        use units, only: boltzmann,convert
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
        real*8, intent(in) :: dt
        integer :: i,j,k,iloc,nu
        real*8 :: dt2, a1,a2,Rn
        integer :: ilocbeg, inoise
        interface
          function normal()
          real*8 normal
          end function
        end interface

        !! apply AOA for all beads local atoms !!
        nu=polymer%nbeads
        dt2=0.5d0*dt
        ilocbeg = ilocpi_beg(rank_polymer+1)

        !! CENTROID !!
        a1 = exp(-gamma_friction(1)*dt)
        a2 = sqrt((1.-a1**2)*boltzmann*kelvin)
        if (piqtb) then
          inoise = get_qtb_compteur()
          do iloc=1,nlocpi
            i=glob(iloc+ilocbeg-1)
            do j=1,3
              polymer%eigpos(j,i,1)=polymer%eigpos(j,i,1) 
     &           + dt2*polymer%eigvel(j,i,1)
              if(register_spectra) then
              qtbdata(i)%vad(j,inoise,1)=real(
     &           sqrt(a1)*polymer%eigvel(j,i,1)
     &          +0.5*dt*qtbdata(i)%fad(j,inoise,1)
     &                  /mass(i),4)
              endif
              polymer%eigvel(j,i,1)=a1*polymer%eigvel(j,i,1)
     &          + dt*qtbdata(i)%fad(j,inoise,1)
     &                /mass(i)
              polymer%eigpos(j,i,1)=polymer%eigpos(j,i,1) 
     &           + dt2*polymer%eigvel(j,i,1)
            enddo
          enddo
        else
          do iloc=1,nlocpi
            i=glob(iloc+ilocbeg-1)
            do j=1,3
              Rn = normal()
              polymer%eigpos(j,i,1)=polymer%eigpos(j,i,1) 
     &           + dt2*polymer%eigvel(j,i,1)
              polymer%eigvel(j,i,1)=polymer%eigvel(j,i,1)*a1 
     &                 + a2*Rn/sqrt(mass(i))
              polymer%eigpos(j,i,1)=polymer%eigpos(j,i,1) 
     &           + dt2*polymer%eigvel(j,i,1)
            enddo
          enddo
        endif

        if(nu==1) return
        !! FLUCTUATION MODES !!

        if(cay_correction) then
          call apply_A_cayhalf(polymer,dt)
        else
          call apply_A_fluct(polymer,dt2)
        endif

        if (piqtb) then
          do k=2,nu
            !! BLOCK O fluctuations!!
            a1     = exp(-gamma_friction(k)*dt)
            a2 = sqrt((1.-a1**2)*boltzmann*kelvin)
            do iloc=1,nlocpi
              i=glob(iloc+ilocbeg-1)
              do j=1,3
                if(register_spectra) then
                qtbdata(i)%vad(j,inoise,k)=real(
     &           s  qrt(a1)*polymer%eigvel(j,i,k)
     &            +0.5*dt*qtbdata(i)%fad(j,inoise,k)
     &                  /mass(i),4)
                endif
                polymer%eigvel(j,i,k)=a1*polymer%eigvel(j,i,k)
     &            + dt*qtbdata(i)%fad(j,inoise,k)
     &                  /mass(i)
              enddo
            enddo
          enddo
        else
          do k=2,nu
            !! BLOCK O fluctuations!!
            a1     = exp(-gamma_friction(k)*dt)
            a2 = sqrt((1.-a1**2)*boltzmann*kelvin)
            do iloc=1,nlocpi
              i=glob(iloc+ilocbeg-1)
              do j=1,3
                Rn = normal()
                polymer%eigvel(j,i,k)=polymer%eigvel(j,i,k)*a1 
     &               + a2*Rn/sqrt(mass(i))
              enddo
            enddo
          enddo
        endif

        if(cay_correction) then
          call apply_A_cayhalf(polymer,dt)
        else
          call apply_A_fluct(polymer,dt2)
        endif

        if(piqtb) call QTB_update_state()

      end subroutine apply_AOA_PI

      subroutine apply_A_cayhalf(polymer,dt)
        use atoms
        use atmtyp
        use beads
        use domdec
        use inform
        use bath
        use langevin
        use moldyn
        use units, only: boltzmann,convert
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
        real*8, intent(in) :: dt
        integer :: i,j,k,iloc,nu
        real*8 :: eigx0,eigv0,cayfact
        integer :: ilocbeg

        nu=polymer%nbeads
        ilocbeg = ilocpi_beg(rank_polymer+1)
        !! BLOCK A fluctuations CAY CORRECTION!!
        do k=2,nu 
          cayfact = sqrt(1.d0/(4.d0+(omkpi(k)*dt)**2))
          do iloc=1,nlocpi
            i=glob(iloc+ilocbeg-1)
            do j=1,3
              eigx0=polymer%eigpos(j,i,k)*cayfact
              eigv0=polymer%eigvel(j,i,k)*cayfact
              polymer%eigpos(j,i,k)=2*eigx0 + dt*eigv0
              polymer%eigvel(j,i,k)=-dt*omkpi(k)**2*eigx0+2*eigv0
            enddo
          enddo
        enddo

      end subroutine apply_A_cayhalf

      subroutine apply_A_fluct(polymer,dt)
        use atoms
        use atmtyp
        use beads
        use domdec
        use inform
        use bath
        use langevin
        use moldyn
        use units, only: boltzmann,convert
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
        real*8, intent(in) :: dt
        integer :: i,j,k,iloc,nu
        real*8 :: eigx0,eigv0
        integer :: ilocbeg

        nu=polymer%nbeads
        ilocbeg = ilocpi_beg(rank_polymer+1)
        !! BLOCK A fluctuations exact integration!!
        do k=2,nu
          do iloc=1,nlocpi
            i=glob(iloc+ilocbeg-1)
            do j=1,3
              eigx0=polymer%eigpos(j,i,k)
              eigv0=polymer%eigvel(j,i,k)
              polymer%eigpos(j,i,k)=eigx0*cos(omkpi(k)*dt)
     $               +eigv0*sin(omkpi(k)*dt)/omkpi(k)
              polymer%eigvel(j,i,k)=eigv0*cos(omkpi(k)*dt)
     $               -eigx0*sin(omkpi(k)*dt)*omkpi(k)
            enddo
          enddo
        enddo

      end subroutine apply_A_fluct

c----------------------------------------------------------------------
c       GRADIENT SUBROUTINES

      subroutine compute_gradslow_centroid(epot,derivs_,vir_
     &       ,only_long,remove_polarshort)
      use atoms
      use beads
      use cutoff
      use deriv
      use domdec
      use energi
      use moldyn
      use potent
      use virial
      use spectra
      use polpot
#ifdef COLVARS
      use colvars
#endif
#ifdef PLUMED
      use plumed
#endif
      implicit none
      real*8, intent(inout) :: epot,vir_(3,3)
      real*8, intent(inout) :: derivs_(3,nbloc)
      logical, intent(in) :: only_long ! 0 = slow, 1 = int+slow
      logical, intent(in) :: remove_polarshort
      integer i,j
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend 
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_angtor
      logical save_tortor,save_geom
      logical save_extra
      logical save_vdw,save_charge
      logical save_mpole
      logical save_mrec
      logical save_prec,save_crec
      logical save_disprec
      logical save_disp,save_chgtrn,save_repuls
      integer save_polalg,save_tcgorder
      real*8 save_tcgomega
      logical save_tcgprec,save_tcgguess,save_tcgpeek
      real*8 epshort,virshort(3,3)
      real*8 diploc(3)
#ifdef COLVARS
      logical save_colvars
#endif
#ifdef PLUMED
      logical save_plumed
#endif

      derivs_(:,:)=0.0d0
c
c     save the original state of fast-evolving potentials
c
      save_bond     = use_bond
      save_angle    = use_angle
      save_strbnd   = use_strbnd
      save_urey     = use_urey
      save_angang   = use_angang
      save_opbend   = use_opbend 
      save_opdist   = use_opdist
      save_improp   = use_improp
      save_imptor   = use_imptor
      save_tors     = use_tors
      save_pitors   = use_pitors
      save_strtor   = use_strtor
      save_angtor   = use_angtor
      save_tortor   = use_tortor
      save_geom     = use_geom
      save_extra    = use_extra
c
c     turn off fast-evolving valence potential energy terms
c
      use_bond      = .false.
      use_angle     = .false.
      use_strbnd    = .false.
      use_urey      = .false.
      use_angang    = .false.
      use_opbend    = .false.
      use_opdist    = .false.
      use_improp    = .false.
      use_imptor    = .false.
      use_tors      = .false.
      use_pitors    = .false.
      use_strtor    = .false.
      use_angtor    = .false.
      use_tortor    = .false.
      use_geom      = .false.
      use_extra     = .false.

      if(only_long .and. remove_polarshort) then
c     FIRST PASS TO SUBSTRACT POLAR SHORTRANGE

      save_vdw = use_vdw
      save_charge = use_charge
      save_mpole = use_mpole
      save_disp = use_disp
      save_chgtrn = use_chgtrn
      save_repuls = use_repuls
      
      save_crec    = use_crec
      save_mrec    = use_mrec
      save_prec    = use_prec
      save_disprec = use_disprec
      save_polalg = polalg
      save_tcgorder = tcgorder
      save_tcgprec = tcgprec
      save_tcgguess = tcgguess
      save_tcgpeek = tcgpeek
      save_tcgomega = tcgomega
#ifdef COLVARS
      save_colvars = use_colvars
#endif
#ifdef PLUMED
      save_plumed = lplumed
#endif

      use_vdw = .false.
      use_charge = .false.
      !use_mpole = .false.
      use_disp = .false.
      use_chgtrn = .false.
      use_repuls = .false.

      use_crec     = .false.
      use_mrec     = .false.
      use_prec     = .false.
      use_disprec  = .false.
      use_cself    = .false.
      use_mself    = .false.
      use_pself    = .false.
      use_dispself = .false.
      use_cshortreal     = .true.
      use_mpoleshortreal = .true.
      use_vdwshort       = .true.
      use_polarshortreal = .true.
      use_repulsshort = .true.
      use_dispshort = .true.
      use_dispshortreal = .true.
      use_chgtrnshort = .true.
      polalg = polalgshort
      tcgorder = tcgordershort
      tcgprec = tcgprecshort
      tcgguess = tcgguessshort
      tcgpeek = tcgpeekshort
      tcgomega = tcgomegashort
#ifdef COLVARS
      use_colvars = .false.
#endif
#ifdef PLUMED
      lplumed = .false.
#endif

      call gradient(epshort,derivs_)

      use_vdw = save_vdw
      use_charge = save_charge
      !use_mpole = save_mpole
      use_disp = save_disp
      use_chgtrn = save_chgtrn
      use_repuls = save_repuls

      use_crec   = save_crec
      use_mrec   = save_mrec
      use_prec   = save_prec
      use_disprec = save_disprec
      use_cself  = .true.
      use_mself  = .true.
      use_pself  = .true.
      use_dispself = .true.
      use_cshortreal     = .false.
      use_mpoleshortreal = .false.
      use_vdwshort       = .false.
      use_polarshortreal = .false.
      use_repulsshort = .false.
      use_dispshort = .false.
      use_dispshortreal = .false.
      use_chgtrnshort = .false.

      polalg = save_polalg
      tcgorder = save_tcgorder
      tcgprec = save_tcgprec
      tcgguess = save_tcgguess
      tcgpeek = save_tcgpeek
#ifdef COLVARS
      use_colvars = save_colvars
#endif
#ifdef PLUMED
      lplumed = save_plumed
#endif

      DO i=1,nbloc; do j=1,3
        derivs_(j,i) = -dep(j,i)
      ENDDO ; ENDDO
      epshort = ep
      virshort(:,:) = vir(:,:) 

      endif

c     SECOND PASS TO ADD LONGRANGE

      if(only_long) then
      use_vdwlong   = .true.
      use_clong     = .true.
      use_mpolelong = .true.
      use_repulslong = .true.
      use_displong = .true.
      use_chgtrnlong = .true.
      endif
c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (epot,derivs_)
      vir_(:,:) = vir(:,:) 
      if(ir .and. rank_polymer==0) then
        call compute_dipole(diploc,dipindpi,full_dipole)
      else
        dipindpi=0
      endif

      if(only_long .and. remove_polarshort) then
        epot = epot - epshort
        vir_(:,:) = vir_(:,:)-virshort(:,:)
      endif
c
c     restore the original state of fast-evolving potentials
c
      use_bond      = save_bond
      use_angle     = save_angle
      use_strbnd    = save_strbnd
      use_urey      = save_urey
      use_angang    = save_angang
      use_opbend    = save_opbend 
      use_opdist    = save_opdist
      use_improp    = save_improp
      use_imptor    = save_imptor
      use_tors      = save_tors
      use_pitors    = save_pitors
      use_strtor    = save_strtor
      use_angtor    = save_angtor
      use_tortor    = save_tortor
      use_geom      = save_geom
      use_extra     = save_extra

      if(only_long) then
      use_vdwlong   = .false.
      use_clong     = .false.
      use_mpolelong = .false.
      use_repulslong = .false.
      use_displong = .false.
      use_chgtrnlong = .false.
      endif


      end subroutine compute_gradslow_centroid

      subroutine compute_grad_beads(epot,derivs_,vir_,polymer
     &   ,do_long,do_int,do_short,compute_polar)
      use atoms
      use beads
      use cutoff
      use domdec
      use deriv
      use energi
      use polpot
      use potent
      use commstuffpi
      use virial
      use mpi
      use spectra
#ifdef COLVARS
      use colvars
#endif
#ifdef PLUMED
      use plumed
#endif
      implicit none
      real*8, intent(inout) ::  epot,vir_(3,3)
      real*8, intent(inout) :: derivs_(:,:,:)
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: do_long, do_int, do_short
      logical, intent(in) :: compute_polar
      real*8 save_tcgomega
      integer i,j,k,ibead
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend 
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_angtor
      logical save_tortor,save_geom
      logical save_extra
      logical save_mrec
      logical save_prec,save_crec,save_disprec
      logical save_smdvel, save_smdfor, save_list
      logical save_polar,save_vdw,save_charge,save_mpole
      logical save_repuls,save_disp,save_chgtrn,save_solv
      integer save_polalg,save_tcgorder
      logical save_tcgprec,save_tcgguess,save_tcgpeek
#ifdef COLVARS
      logical save_colvars
#endif
#ifdef PLUMED
      logical save_plumed
#endif
      real*8 nu
      integer nbeadsproc,ibeadbeg,ibeadend
      logical :: fullpotential,bonded_only


      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1

      do k=1,nbeadsproc; do i=1,nbloc; do j=1,3
        derivs_(j,i,k) = 0.0d0
      enddo; enddo; enddo
      if(.not. ANY([do_short,do_int,do_long])) return

      nu=real(polymer%nbeads,8)
      fullpotential = do_short .and. do_int .and. do_long
      bonded_only=.not.(do_long.or.do_int)

      if(do_long .and. (.not. do_int)) then
        if(ranktot==0) then
          write(*,*) "ERROR: compute_grad_beads can't do "
     &              ,"long range without int range"
          call fatal
        endif
      endif

      if(do_long) then
        call resetForcesRec
        if(allocated(derivsRec)) deallocate(derivsRec)
        allocate(derivsRec(3,nlocrec2,nbeadsproc))
        derivsRec(:,:,:) = 0.0d0
      endif
c
c     save the original state of fast-evolving potentials
c

      if(.not.do_short)then
        ! deactivate bonded terms
        save_bond   = use_bond
        save_angle  = use_angle
        save_strbnd = use_strbnd
        save_urey   = use_urey
        save_angang = use_angang
        save_opbend = use_opbend 
        save_opdist = use_opdist
        save_improp = use_improp
        save_imptor = use_imptor
        save_tors   = use_tors
        save_pitors = use_pitors
        save_strtor = use_strtor
        save_angtor = use_angtor
        save_tortor = use_tortor
        save_geom   = use_geom
        save_extra  = use_extra
        use_bond   = .false.
        use_angle  = .false.
        use_strbnd = .false.
        use_urey   = .false.
        use_angang = .false.
        use_opbend = .false.
        use_opdist = .false.
        use_improp = .false.
        use_imptor = .false.
        use_tors   = .false.
        use_pitors = .false.
        use_strtor = .false.
        use_angtor = .false.
        use_tortor = .false.
        use_geom   = .false.
        use_extra  = .false.
      endif

      if(bonded_only) then
        ! deactivate non-bonded terms
        save_vdw = use_vdw
        save_charge = use_charge
        save_mpole = use_mpole
        save_polar = use_polar
        save_repuls = use_repuls
        save_disp   = use_disp
        save_chgtrn = use_chgtrn
        save_solv   = use_solv
        save_list = use_list
        save_smdvel = use_smd_velconst
        save_smdfor = use_smd_forconst
#ifdef COLVARS
        save_colvars = use_colvars
#endif
#ifdef PLUMED
        save_plumed = lplumed
#endif
c
c     turn off slow-evolving nonbonded potential energy terms
c
        use_vdw = .false.
        use_charge = .false.
        use_mpole = .false.
        use_polar = .false.
        use_list = .false.
        use_solv   = .false.
        use_repuls = .false.
        use_disp   = .false.
        use_chgtrn = .false.
        use_smd_velconst = .false.
        use_smd_forconst = .false.
#ifdef COLVARS
        use_colvars = .false.
#endif
#ifdef PLUMED
        lplumed = .false.
#endif

      elseif(do_int .and.(.not. do_long)) then
        ! deactivate long-range terms (keep short-range non-bonded)
        save_crec   = use_crec
        save_mrec   = use_mrec
        save_prec   = use_prec
        save_disprec  = use_disprec
        save_polalg = polalg
        save_tcgorder = tcgorder
        save_tcgprec  = tcgprec
        save_tcgguess = tcgguess
        save_tcgpeek  = tcgpeek
        save_tcgomega = tcgomega
        save_smdvel = use_smd_velconst
        save_smdfor = use_smd_forconst
        save_polar = use_polar
#ifdef COLVARS
        save_colvars = use_colvars
#endif
#ifdef PLUMED
        save_plumed = lplumed
#endif
c  
c       turn off fast-evolving valence potential energy terms
c  
        use_crec   = .false.
        use_mrec   = .false.
        use_prec   = .false.
        use_disprec= .false.
        use_cself  = .false.
        use_mself  = .false.
        use_pself  = .false.
        use_dispself  = .false.

        use_polar  = compute_polar

        use_cshortreal     = .true.
        use_mpoleshortreal = .true.
        use_vdwshort       = .true.
        use_polarshortreal = compute_polar
        use_repulsshort    = .true.
        use_dispshort      = .true.
        use_dispshortreal  = .true.
        use_chgtrnshort    = .true.
        polalg     = polalgshort
        tcgorder   = tcgordershort
        tcgprec    = tcgprecshort
        tcgguess   = tcgguessshort
        tcgpeek    = tcgpeekshort
        tcgomega   = tcgomegashort
        use_smd_velconst   = .false.
        use_smd_forconst   = .false.
#ifdef COLVARS
        use_colvars = .false.
#endif
#ifdef PLUMED
        lplumed = .false.
#endif

      endif

      epot = 0.0d0
      vir_(:,:)=0.0d0
      polymer%epot(:)=0.0d0
      if (ir) then
        polymer%dip(:,:)=0.0d0
        polymer%dipind(:,:)=0.0d0
        if(do_short) dippi=0.0d0
        if(do_long) dipindpi=0.0d0
      endif

      do k=1,nbeadsproc
        ibead = k+ibeadbeg-1
        call load_bead(polymer, ibead)
        call gradient (polymer%epot(k),derivs_(:,:,k))
        if(do_long) call save_deriv_rec(derivsRec(:,:,k))

        if(ir) then
          if(bonded_only) then
            use_mpole = save_mpole
            use_charge = save_charge
            !if(save_mpole .and. full_dipole) call rotpole
          endif
          call compute_dipole(polymer%dip(:,k),polymer%dipind(:,k)
     &        ,full_dipole.and.(.not. bonded_only))
          if(do_short) dippi=dippi+polymer%dip(:,k)/nu
          if(do_long) dipindpi=dipindpi+polymer%dipind(:,k)/nu
          if(bonded_only) then
            use_mpole=.FALSE.
            use_charge=.FALSE.
          endif
        endif
        
        epot = epot+polymer%epot(k)/nu
        polymer%vir(:,:,k)=vir(:,:)
        vir_(:,:)=vir_(:,:)+vir(:,:)/nu
      enddo

c
c     restore the original state of fast-evolving potentials
c
      if(.not.do_short)then
        use_bond   = save_bond
        use_angle  = save_angle
        use_strbnd = save_strbnd
        use_urey   = save_urey
        use_angang = save_angang
        use_opbend = save_opbend 
        use_opdist = save_opdist
        use_improp = save_improp
        use_imptor = save_imptor
        use_tors   = save_tors
        use_pitors = save_pitors
        use_strtor = save_strtor
        use_angtor = save_angtor
        use_tortor = save_tortor
        use_geom   = save_geom
        use_extra  = save_extra
      endif

      if(bonded_only)then
        use_vdw = save_vdw
        use_charge = save_charge
        use_mpole = save_mpole
        use_polar = save_polar
        use_repuls = save_repuls
        use_disp   = save_disp
        use_chgtrn = save_chgtrn
        use_solv   = save_solv
        use_list = save_list
        use_smd_velconst = save_smdvel
        use_smd_forconst = save_smdfor
#ifdef COLVARS
        use_colvars = save_colvars
#endif
#ifdef PLUMED
        lplumed = save_plumed
#endif

      elseif(do_int .and.(.not. do_long)) then

        use_crec   = save_crec
        use_mrec   = save_mrec
        use_prec   = save_prec
        use_disprec= save_disprec
        use_cself  = .true.
        use_mself  = .true.
        use_pself  = .true.
        use_dispself = .true.
        use_cshortreal     = .false.
        use_mpoleshortreal = .false.
        use_vdwshort       = .false.
        use_polarshortreal = .false.
        use_repulsshort = .false.
        use_dispshort = .false.
        use_dispshortreal = .false.
        use_chgtrnshort = .false.
        polalg     = save_polalg
        tcgorder   = save_tcgorder
        tcgprec    = save_tcgprec
        tcgguess   = save_tcgguess
        tcgpeek    = save_tcgpeek
        use_smd_velconst = save_smdvel
        use_smd_forconst = save_smdfor
        use_polar = save_polar
#ifdef COLVARS
        use_colvars = save_colvars
#endif
#ifdef PLUMED
        lplumed = save_plumed
#endif
      endif

      end subroutine compute_grad_beads

c----------------------------------------------------------------------- 
c     BAROSTAT SUBROUTINES

      subroutine plangevinpi(polymer,dt,istep,bloc)
      use boxes
      use bath
      use beads
      use langevin
      use domdec
      use mpi
      use units
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      REAL*8, intent(in) :: dt
      integer, intent(in) :: istep
      character, intent(in) :: bloc
      integer :: iloc,ierr,i
      real*8 :: factor,a1p,a2p,scale(3),tmp
      real*8 :: normal

      SELECT CASE(bloc)
      CASE('B')
        if(ranktot>0) return
        if (anisotrop) then
          do i=1,3
            factor=convert/masspiston
            aextbox(i)=factor*(extvol*(stresspi(i,i)-atmsph)/prescon 
     &            +gasconst*kelvin)
            vextbox(i)=vextbox(i) + dt*aextbox(i)
          enddo
        else
          factor=3.d0*convert/masspiston
          aextvol=factor*(
     &          extvol*(presvir-atmsph)/prescon 
     &          +gasconst*kelvin )
          vextvol=vextvol + dt*aextvol
        endif
      CASE('O')
        if(ranktot>0) return
        a1p = exp(-gammapiston*dt)
        a2p = sqrt((1-a1p**2)*boltzmann*kelvin/masspiston)
        if (anisotrop) then
          do i=1,3
            vextbox(i)=a1p*vextbox(i)+a2p*normal()
          enddo
          temppiston = masspiston*sum(vextbox**2)/(3.*boltzmann) 
        else
          vextvol = vextvol*a1p + a2p*normal()
          temppiston = masspiston*vextvol**2/boltzmann
        endif
      CASE('A')
         if (anisotrop) then
          if(ranktot==0) scale(:) = exp(dt*vextbox(:))
          call MPI_BCAST(scale,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        else
          if(ranktot==0) tmp=exp(dt*vextvol)
          call MPI_BCAST(tmp,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          scale(:)=tmp
        endif
        extvol = extvol*scale(1)*scale(2)*scale(3)
        call rescale_box(istep,scale)
        DO iloc=1,nloc
          i=glob(iloc)
          polymer%eigpos(:,i,1)=polymer%eigpos(:,i,1)*scale(:)
          polymer%eigvel(:,i,1)=polymer%eigvel(:,i,1)/scale(:)
          !polymer%eigforces(:,i,1)=polymer%eigforces(:,i,1)/scale
        ENDDO
      CASE DEFAULT
        if(ranktot==0) then
          write(*,*) 'ERROR: unknown block in plangevinpi'
          write(*,*) 'bloc = ',bloc
          call fatal
        endif
      END SELECT

      
      end subroutine plangevinpi

      subroutine pscalepi (polymer,dt,istep)
      use atoms
      use bath
      use boxes
      use domdec
      use math
      use usage
      use beads
      use mpi
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      REAL*8, intent(in) :: dt
      integer, intent(in) :: istep
      integer i,iloc,ierr
      real*8 scale,third
      parameter(third = 1.0d0 / 3.0d0)
c
c     find the isotropic scale factor for constant pressure
c
      if(nproctot>1) then
        call MPI_BCAST(presvir,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      endif
      scale = (1.0d0 + (dt*compress/taupres)
     &                    *(presvir-atmsph))**third
c
c     modify the current periodic box dimension values
c
         xbox = xbox * scale
         ybox = ybox * scale
         zbox = zbox * scale
c
c     propagate the new box dimensions to other lattice values
c
         call lattice
c
c     couple to pressure bath via atom scaling in Cartesian space
c
        DO iloc=1,nloc
          i=glob(iloc)
          polymer%eigpos(:,i,1)=polymer%eigpos(:,i,1)*scale
          polymer%eigvel(:,i,1)=polymer%eigvel(:,i,1)/scale
        ENDDO 
c
c   also rescale xbegproc, xendproc...
c
         call ddpme3dnpt(scale,istep)
c
      end subroutine pscalepi

      end module utilbaoabpi

