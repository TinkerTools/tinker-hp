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
#include "tinker_macro.h"
      module utilbaoabpi        
      implicit none

      real(r_p),allocatable, save::derivs(:,:,:)
      real(r_p),allocatable, save::derivs_ctr(:,:,:)
      real(r_p),allocatable, save::derivs_centroid(:,:)
      real(t_p), allocatable, save:: noise(:)
      real(r_p),allocatable, save::derivsRec(:,:,:)
      real(r_p) :: eslow,vir_ctr(3,3)

      real(8), save :: timepush=0.d0, timeinit=0.d0
      real(8), save :: timereass=0.d0, timegrad=0.d0
      real(8), save :: timesavenl=0.d0, timecom=0.d0
      real(8), save :: timecontract=0d0,timegradfast=0.d0
      real(8), save :: timeobs=0d0,timeaoa=0d0

      contains

#include "convert.f.inc"
      subroutine apply_B_PI(polymer,dt)
      use beads
      use units, only: convert
      use domdec
      use atmtyp, only: mass
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      real(r_p), intent(in) :: dt
      integer :: nu,iloc, i,j,k,ilocbeg

      nu=polymer%nbeads
      ilocbeg = ilocpi_beg(rank_polymer+1)
!$acc parallel loop collapse(3) default(present) async
      do k=1,nu; do iloc=1,nlocpi; do j=1,3
        i=glob(iloc+ilocbeg-1)
        !! BLOCK B slow !!
        polymer%eigvel(j,i,k)=polymer%eigvel(j,i,k)
     &      +dt*convert*polymer%eigforces(j,i,k)/mass(i)
      enddo; enddo; enddo

      end subroutine apply_B_PI

      subroutine apply_AOA_PI(polymer,dt)
        use atomsMirror
        use atmtyp
        use beads
        use domdec
        use inform
        use random_mod
        use bath
        use qtb
        use langevin
        use moldyn
        use commstuffpi
        use units, only: boltzmann,convert
        implicit none
        TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
        real(r_p), intent(in) :: dt
        integer :: i,j,k,iloc,kk,nu
        real(r_p) :: dt2, eigx0,eigv0,gammak,a1,a2,cayfact
        integer :: ilocbeg,inoise

        nu=polymer%nbeads
        dt2=0.5d0*dt
        if (piqtb) then
          inoise = get_qtb_compteur()
        else
      !! Generate noise for all beads local atoms !!
          call normalarray(noise,3*nlocpi*nu)
        endif

      !! apply AOA for all beads local atoms !!
        ilocbeg = ilocpi_beg(rank_polymer+1)
!$acc parallel loop collapse(3) default(present) async
        do k=1,nu; do iloc=1,nlocpi; do j=1,3
          i=glob(iloc+ilocbeg-1)
          !! BLOCK B !!
c          polymer%eigvel(j,i,k)=polymer%eigvel(j,i,k)
c     &       +dt2*convert*polymer%eigforces(j,i,k)/mass(i)
          if(k==1) then
            !! BLOCK A centroid !!
            polymer%eigpos(j,i,1)=polymer%eigpos(j,i,1) 
     &       + dt2*polymer%eigvel(j,i,1)
          elseif(cay_correction) then
            !! BLOCK A fluctuations CAY CORRECTION!!
            cayfact = sqrt(1.d0/(4.d0+(omkpi(k)*dt)**2))
            eigx0=polymer%eigpos(j,i,k)*cayfact
            eigv0=polymer%eigvel(j,i,k)*cayfact
            polymer%eigpos(j,i,k)=2*eigx0 + dt*eigv0
            polymer%eigvel(j,i,k)=-dt*omkpi(k)**2*eigx0+2*eigv0
          else
            !! BLOCK A fluctuations !!
            eigx0=polymer%eigpos(j,i,k)
            eigv0=polymer%eigvel(j,i,k)
            polymer%eigpos(j,i,k)=eigx0*cos(omkpi(k)*dt2)
     $            +eigv0*sin(omkpi(k)*dt2)/omkpi(k)
            polymer%eigvel(j,i,k)=eigv0*cos(omkpi(k)*dt2)
     $           -eigx0*sin(omkpi(k)*dt2)*omkpi(k)
          endif

          !! BLOCK O !!
          a1=exp(-gamma_friction(k)*dt)
          if(piqtb) then
            if(register_spectra) then
              qtbdata(i)%vad(j,inoise,k)=real(
     &           sqrt(a1)*polymer%eigvel(j,i,k)
     &          +0.5*dt*qtbdata(i)%fad(j,inoise,k)
     &                  /mass(i),4)
            endif
            polymer%eigvel(j,i,k)=a1*polymer%eigvel(j,i,k)
     &        + dt*qtbdata(i)%fad(j,inoise,k)
     &              /mass(i)
          else
            kk=(k-1)*3*nlocpi
            a2 = sqrt((1.-a1**2)*boltzmann
     &                    *kelvin/mass(i))
            polymer%eigvel(j,i,k)=polymer%eigvel(j,i,k)*a1 
     &               + a2*noise(j+3*(iloc-1)+kk)
          endif


          if(k==1) then
            !! BLOCK A centroid !!
            polymer%eigpos(j,i,1)=polymer%eigpos(j,i,1) 
     &       + dt2*polymer%eigvel(j,i,1)
          elseif(cay_correction) then
            !! BLOCK A fluctuations CAY CORRECTION!!
            cayfact = sqrt(1.d0/(4.d0+(omkpi(k)*dt)**2))
            eigx0=polymer%eigpos(j,i,k)*cayfact
            eigv0=polymer%eigvel(j,i,k)*cayfact
            polymer%eigpos(j,i,k)=2*eigx0 + dt*eigv0
            polymer%eigvel(j,i,k)=-dt*omkpi(k)**2*eigx0+2*eigv0
          else
            !! BLOCK A fluctuations !!
            eigx0=polymer%eigpos(j,i,k)
            eigv0=polymer%eigvel(j,i,k)
            polymer%eigpos(j,i,k)=eigx0*cos(omkpi(k)*dt2)
     $            +eigv0*sin(omkpi(k)*dt2)/omkpi(k)
            polymer%eigvel(j,i,k)=eigv0*cos(omkpi(k)*dt2)
     $           -eigx0*sin(omkpi(k)*dt2)*omkpi(k)
          endif
        enddo; enddo; enddo

        if(piqtb) call QTB_update_state()


      end subroutine apply_AOA_PI

c----------------------------------------------------------------------
c       GRADIENT SUBROUTINES

      subroutine compute_gradslow_centroid(epot,derivs_,vir_
     &       ,only_long,remove_polarshort)
      use atomsMirror
      use beads
      use cutoff
      use deriv
      use domdec
      use energi
      use moldyn
      use potent
      use virial
      use tinheader
      use utilgpu
      use inform
      use utils
      use tinheader
      use spectra
      implicit none
      real(r_p), intent(inout) :: epot,vir_(3,3)
      real(r_p), intent(inout) :: derivs_(3,nbloc)
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
      logical save_mpole,save_polar
      logical save_list
      logical save_smdvel, save_smdfor
      logical save_mrec
      logical save_prec,save_crec,save_disprec
      integer save_polalg,save_tcgorder
      logical save_tcgprec,save_tcgguess,save_tcgpeek
      logical save_disp,save_repuls,save_chgtrn
      real(t_p) save_tcgomega
      real(r_p) epshort,virshort(3,3)
      real(r_p) diploc(3)

      call set_to_zero1m (derivs_,3*nbloc,rec_queue)

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

      bonded_l = .false.

      if(only_long) call detach_mpolar_pointer

      if(only_long .and. remove_polarshort) then
!$acc enter data create(epshort,virshort) async
c     FIRST PASS TO SUBSTRACT POLAR SHORTRANGE

      save_vdw = use_vdw
      save_charge = use_charge
      save_mpole = use_mpole
      save_disp = use_disp
      save_repuls = use_repuls
      save_chgtrn = use_chgtrn
      
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

      use_vdw = .false.
      use_charge = .false.
      !use_mpole = .false.

      use_crec      = .false.
      use_mrec      = .false.
      use_prec      = .false.
      use_disprec   = .false.
      use_cself     = .false.
      use_mself     = .false.
      use_pself     = .false.
      use_dispself  = .false.
      use_cshortreal     = .true.
      use_mpoleshortreal = .true.
      use_vdwshort       = .true.
      use_polarshortreal = .true.
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
      shortnonbonded_l   = .true.

      call gradient(epshort,derivs_)

      use_vdw = save_vdw
      use_charge = save_charge
      use_mpole = save_mpole

      use_crec       = save_crec
      use_mrec       = save_mrec
      use_prec       = save_prec
      use_disprec    = save_disprec
      use_cself      = .true.
      use_mself      = .true.
      use_pself      = .true.
      use_dispself   = .true.
      use_cshortreal     = .false.
      use_mpoleshortreal = .false.
      use_vdwshort       = .false.
      use_polarshortreal = .false.
      use_repulsshort    = .false.
      use_dispshort      = .false.
      use_dispshortreal  = .false.
      use_chgtrnshort    = .false.
      polalg     = save_polalg
      tcgorder   = save_tcgorder
      tcgprec    = save_tcgprec
      tcgguess   = save_tcgguess
      tcgpeek    = save_tcgpeek
      tcgomega   = save_tcgomega
      use_smd_velconst = save_smdvel
      use_smd_forconst = save_smdfor
      shortnonbonded_l   = .false.

c!$acc parallel loop collapse(2) async default(present)
c      DO i=1,nbloc; do j=1,3
c        derivs_(j,i) = -mdr2md(dep(j,i))
c      ENDDO ; ENDDO
      call set_to_zero1m (derivs_,3*nbloc,rec_queue)
      if (calc_e.or.use_virial) then
!$acc serial present(epshort,ep,virshort,virsave) async
        epshort = ep
        do j=1,3; do i=1,3
          virshort(i,j) = virsave(i,j)
          virsave(i,j) = 0.0
        enddo; enddo
!$acc end serial
      endif

      endif

c     SECOND PASS TO ADD LONGRANGE

      if(only_long) then
      use_vdwlong   = .true.
      use_clong     = .true.
      use_repulslong= .true.
      use_displong  = .true.
      use_mpolelong = .true.
      use_chgtrnlong= .true.
      endif
c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (epot,derivs_)
!$acc serial async present(vir_,vir) 
      do j=1,3; do i=1,3
        vir_(i,j) = vir(i,j)
      enddo; enddo
!$acc end serial
      if(ir .and. rank_polymer==0) then
!$acc data create(diploc) async
        call compute_dipole(diploc,dipindpi,full_dipole)
!$acc end data
      else
!$acc serial async
        dipindpi(:)=0
!$acc end serial
      endif

      if(only_long .and. remove_polarshort) then
        call remove_desave(derivs_)
        if (calc_e.or.use_virial) then
!$acc serial async present(epot,epshort,vir_,virshort) 
          epot = epot - epshort
          do j=1,3; do i=1,3
            vir_(i,j) = vir_(i,j)-virshort(i,j)
          enddo; enddo
!$acc end serial
        endif

!$acc exit data delete(epshort,virshort) async
      endif

      if (deb_Energy) call info_energy(rank)

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

      bonded_l = .true.

      if(only_long) then
      use_vdwlong   = .false.
      use_clong     = .false.
      use_mpolelong = .false.
      use_repulslong= .false.
      use_displong  = .false.
      use_chgtrnlong= .false.
      call attach_mpolar_pointer
      endif


      end subroutine compute_gradslow_centroid

      subroutine compute_grad_beads(epot,derivs_,vir_,polymer
     &   ,do_long,do_int,do_short,compute_polar)
      use atomsMirror
      use ascii, only: int_to_str
      use beads
      use cutoff
      use domdec
      use deriv
      use energi
      use inform
      use polpot
      use potent
      use utilgpu
      use utils
      use mpi
      use inform, only: deb_path
      use deriv, only: zero_forces_rec
      use commstuffpi
      use virial
      use spectra
      implicit none
      real(r_p), intent(inout) ::  epot,vir_(3,3)
      real(r_p), intent(inout) :: derivs_(:,:,:)
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: do_long, do_int, do_short
      logical, intent(in) :: compute_polar
      real(t_p) save_tcgomega
      integer i,j,k,iloc,ibead
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
      logical save_prec,save_crec
      integer save_polalg,save_tcgorder
      logical save_tcgprec,save_tcgguess,save_tcgpeek
      logical save_smdvel, save_smdfor, save_list
      logical save_polar,save_vdw,save_charge,save_mpole
      logical save_repuls,save_disp,save_chgtrn,save_solv,save_disprec
      real(r_p) eloc,nu
      integer nbeadsproc,ibeadbeg,ibeadend ,ierr
      real*8 time0,time1
      logical :: fullpotential,bonded_only
      real(r_p) diploc(3),dipindloc(3)

      if (deb_Path) write(*,*) '   >> compute_grad_beads'

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1

      call set_to_zero1m(derivs_,3*nbloc*nbeadsproc,rec_queue)
      if(.not. ANY([do_short,do_int,do_long])) return
      
      nu=real(polymer%nbeads,r_p)
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
        call zero_forces_rec
        call prmem_requestm(derivsRec,3,nlocrec2
     &     ,nbeadsproc,async=.true.)
        call set_to_zero1m(derivsRec,3*nlocrec2*nbeadsproc,rec_queue)
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

        bonded_l=.false.
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
c
c     turn off slow-evolving nonbonded potential energy terms
c
        use_vdw = .false.
        use_charge = .false.
        use_mpole = .false.
        use_polar = .false.
        use_solv   = .false.
        use_repuls = .false.
        use_disp   = .false.
        use_chgtrn = .false.
        use_list = .false.
        use_smd_velconst = .false.
        use_smd_forconst = .false.

        nonBonded_l = .false.

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
        use_dispself       = .false.

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

        if(.not.compute_polar) call detach_mpolar_pointer

      endif
c
c     get energy and gradient for slow-evolving potential terms
c

        do k=1,nbeadsproc
          ibead = k+ibeadbeg-1
          if (deb_Path) write(*,*) ranktot,' >> forcebead',ibead
          call load_bead(polymer, ibead)
          call gradient (polymer%epot(k),derivs_(:,:,k))
          !if(do_long) call save_deriv_rec(derivsRec(:,:,k))
          if(do_long) call add_forces_rec1(derivsRec(:,:,k))

          if(ir) then
            if(bonded_only) then
              use_mpole = save_mpole
              use_charge = save_charge
              !if(save_mpole .and. full_dipole) call rotpolegpu
            endif
            call compute_dipole(polymer%dip(:,k),polymer%dipind(:,k)
     &            ,full_dipole.and.(.not. bonded_only))
            if(bonded_only) then
              use_mpole=.FALSE.
              use_charge=.FALSE.
            endif
          endif
          
          if (use_virial) then
!$acc serial async default(present) present(vir)
            polymer%vir(:,:,k)=vir(:,:)
!$acc end serial
          end if


          if (deb_Energy) call info_energy(rank)
          if (deb_Force)  then
            if(.not.do_short) then
             call info_forces(cNBond)
            elseif(bonded_only) then
              call info_forces(cBond)
            else
              call info_forces(cDef)
            endif
            if(rank==0) write(*,*) 'were Forces for bead:'
     &         ,k+ibeadbeg-1
          endif
          if (deb_Atom)   call info_minmax_pva
          
          if (deb_Path) write(*,*) ranktot,' << forcebead',k+ibeadbeg-1
c          if(.not.do_short) then
c          call minmaxone(derivs_(:,:,k),3*nbloc
c     &       ,"d"//int_to_str(k+ibeadbeg-1))
c          call minmaxone(dev,3*nbloc
c     &       ,"dev"//int_to_str(k+ibeadbeg-1))
c          call minmaxone(dem,3*nbloc
c     &       ,"dem"//int_to_str(k+ibeadbeg-1))
c
c          endif
        enddo

        if(use_virial.or.calc_e.or.ir) then
!$acc serial async default(present)
!$acc& present(vir_,dippi,dipindpi,epot)
          vir_(:,:) = 0.0_re_p
          if (use_virial) then
            do k=1,nbeadsproc
              vir_(:,:)=vir_(:,:)+polymer%vir(:,:,k)/nu
              polymer%vir(:,:,k)=0.0_re_p
            enddo
          endif
          if(ir) then
            if(do_short) then
              dippi(:)=0.0_re_p
              do k=1,nbeadsproc
                dippi(:)=dippi(:)+polymer%dip(:,k)/nu
              enddo
            endif
            if(do_long) then
              dipindpi(:)=0.0_re_p
              do k=1,nbeadsproc
                dipindpi(:)=dipindpi(:)+polymer%dipind(:,k)/nu
              enddo
            endif
          endif
          epot=0.0_re_p
          if(calc_e) then
            do k=1,nbeadsproc
              epot=epot+polymer%epot(k)/nu
              polymer%epot(k)=0.0_re_p
            enddo
          endif
!$acc end serial
        endif

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

        bonded_l=.true.
      endif

      if(bonded_only)then
        use_vdw = save_vdw
        use_charge = save_charge
        use_mpole = save_mpole
        use_polar = save_polar
        use_disp   = save_disp
        use_chgtrn = save_chgtrn
        use_repuls = save_repuls
        use_solv   = save_solv
        use_list = save_list
        use_smd_velconst = save_smdvel
        use_smd_forconst = save_smdfor

        nonBonded_l = .true.

      elseif(do_int .and.(.not. do_long)) then
        if(.not.compute_polar) call attach_mpolar_pointer

        use_crec   = save_crec
        use_mrec   = save_mrec
        use_prec   = save_prec
        use_disprec= save_disprec
        use_cself  = .true.
        use_mself  = .true.
        use_pself  = .true.
        use_dispself       = .true.
        use_cshortreal     = .false.
        use_mpoleshortreal = .false.
        use_vdwshort       = .false.
        use_polarshortreal = .false.
        use_repulsshort    = .false.
        use_dispshort      = .false.
        use_dispshortreal  = .false.
        use_chgtrnshort    = .false.
        polalg     = save_polalg
        tcgorder   = save_tcgorder
        tcgprec    = save_tcgprec
        tcgguess   = save_tcgguess
        tcgpeek    = save_tcgpeek
        use_smd_velconst = save_smdvel
        use_smd_forconst = save_smdfor
        use_polar = save_polar
      endif

      if (deb_Path) then
!$acc wait
        write(*,*) '   << compute_grad_beads'
      endif

      end subroutine compute_grad_beads
c-----------------------------------------------------------------------
c    ENERGY SUBROUTINES
      subroutine compute_eslow_centroid(epot,only_long
     &   ,remove_polarshort)
      use atomsMirror
      use beads
      use cutoff
      use deriv
      use domdec
      use energi
      use moldyn
      use potent
      use virial
      use tinheader
      use utilgpu
      use inform
      use utils
      use tinheader
      use spectra
      implicit none
      real(r_p), intent(inout) :: epot
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
      logical save_mpole,save_polar
      logical save_list
      logical save_smdvel, save_smdfor
      logical save_mrec
      logical save_prec,save_crec,save_disprec
      integer save_polalg,save_tcgorder
      logical save_tcgprec,save_tcgguess,save_tcgpeek
      logical save_mlpot
      real(t_p) save_tcgomega
      real(r_p) epshort,virshort(3,3)
      real(r_p) diploc(3)
      interface
       function energy
        real(r_p) energy
       end function energy
      end interface

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
      save_mlpot    = use_mlpot
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
      use_mlpot     = .false.

      bonded_l = .false.

      if(only_long) call detach_mpolar_pointer

      if(only_long .and. remove_polarshort) then
c     FIRST PASS TO SUBSTRACT POLAR SHORTRANGE

      save_vdw = use_vdw
      save_charge = use_charge
      save_mpole = use_mpole
      
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

      use_vdw = .false.
      use_charge = .false.
      !use_mpole = .false.

      use_crec   = .false.
      use_mrec   = .false.
      use_prec   = .false.
      use_disprec= .false.
      use_cself  = .false.
      use_mself  = .false.
      use_pself  = .false.
      use_dispself       = .false.
      use_cshortreal     = .true.
      use_mpoleshortreal = .true.
      use_vdwshort       = .true.
      use_polarshortreal = .true.
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
      shortnonbonded_l   = .true.

      epshort = energy()

      use_vdw = save_vdw
      use_charge = save_charge
      use_mpole = save_mpole

      use_crec   = save_crec
      use_mrec   = save_mrec
      use_prec   = save_prec
      use_disprec= save_disprec
      use_cself  = .true.
      use_mself  = .true.
      use_pself  = .true.
      use_dispself       = .true.
      use_cshortreal     = .false.
      use_mpoleshortreal = .false.
      use_vdwshort       = .false.
      use_polarshortreal = .false.
      use_repulsshort    = .false.
      use_dispshort      = .false.
      use_dispshortreal  = .false.
      use_chgtrnshort    = .false.
      polalg     = save_polalg
      tcgorder   = save_tcgorder
      tcgprec    = save_tcgprec
      tcgguess   = save_tcgguess
      tcgpeek    = save_tcgpeek
      tcgomega   = save_tcgomega
      use_smd_velconst = save_smdvel
      use_smd_forconst = save_smdfor
      shortnonbonded_l   = .false.

!$acc update host(ep) async
!$acc wait
      epshort = ep

      endif

c     SECOND PASS TO ADD LONGRANGE

      if(only_long) then
      use_vdwlong   = .true.
      use_clong     = .true.
      use_repulslong= .true.
      use_displong  = .true.
      use_mpolelong = .true.
      use_chgtrnlong= .true.
      endif
c
c     get energy and gradient for slow-evolving potential terms
c
      epot = energy()

      if(deb_energy) call info_energy(0)

      if(only_long .and. remove_polarshort) then
        epot = epot - epshort
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
      use_mlpot     = save_mlpot

      bonded_l = .true.

      if(only_long) then
      use_vdwlong   = .false.
      use_clong     = .false.
      use_mpolelong = .false.
      use_repulslong= .false.
      use_displong  = .false.
      use_chgtrnlong= .false.
      call attach_mpolar_pointer
      endif


      end subroutine compute_eslow_centroid

      subroutine compute_energy_beads(epot,polymer
     &   ,do_long,do_int,do_short,compute_polar)
      use atomsMirror
      use ascii, only: int_to_str
      use beads
      use cutoff
      use domdec
      use deriv
      use energi
      use inform
      use polpot
      use potent
      use utilgpu
      use utils
      use mpi
      use inform, only: deb_path
      use deriv, only: zero_forces_rec
      use commstuffpi
      use virial
      use spectra
      implicit none
      real(r_p), intent(inout) ::  epot
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      logical, intent(in) :: do_long, do_int, do_short
      logical, intent(in) :: compute_polar
      real(t_p) save_tcgomega
      integer i,j,k,iloc,ibead
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
      logical save_prec,save_crec
      integer save_polalg,save_tcgorder
      logical save_tcgprec,save_tcgguess,save_tcgpeek
      logical save_smdvel, save_smdfor, save_list
      logical save_polar,save_vdw,save_charge,save_mpole
      logical save_repuls,save_disp,save_chgtrn,save_solv,save_disprec
      logical save_mlpot
      real(r_p) eloc,nu
      integer nbeadsproc,ibeadbeg,ibeadend ,ierr
      real*8 time0,time1
      logical :: fullpotential,bonded_only
      real(r_p) diploc(3),dipindloc(3)
      interface
       function energy
        real(r_p) energy
       end function energy
      end interface

      if (deb_Path) write(*,*) '   >> compute_grad_beads'

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1

      if(.not. ANY([do_short,do_int,do_long])) return
      nu=real(polymer%nbeads,r_p)
      fullpotential = do_short .and. do_int .and. do_long
      bonded_only=.not.(do_long.or.do_int)

      if(do_long .and. (.not. do_int)) then
        if(ranktot==0) then
          write(*,*) "ERROR: compute_grad_beads can't do "
     &              ,"long range without int range"
          call fatal
        endif
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
        save_mlpot  = use_mlpot
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
        use_mlpot  = .false.

        bonded_l=.false.
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
c
c     turn off slow-evolving nonbonded potential energy terms
c
        use_vdw = .false.
        use_charge = .false.
        use_mpole = .false.
        use_polar = .false.
        use_solv   = .false.
        use_repuls = .false.
        use_disp   = .false.
        use_chgtrn = .false.
        use_list = .false.
        use_smd_velconst = .false.
        use_smd_forconst = .false.

        nonBonded_l = .false.

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
        use_dispself       = .false.

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

        if(.not.compute_polar) call detach_mpolar_pointer

      endif
c
c     get energies
c
        epot=0.0_re_p
        do k=1,nbeadsproc
          ibead = k+ibeadbeg-1
          if (deb_Path) write(*,*) ranktot,' >> energybead',ibead
          call load_bead(polymer, ibead)
          polymer%epot(k) = energy()
          epot=epot+polymer%epot(k)/nu

          if (deb_energy) call info_energy(0)

          if (deb_Path) write(*,*) ranktot,' << energybead',ibead
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
        use_mlpot  = save_mlpot

        bonded_l=.true.
      endif

      if(bonded_only)then
        use_vdw = save_vdw
        use_charge = save_charge
        use_mpole = save_mpole
        use_polar = save_polar
        use_disp   = save_disp
        use_chgtrn = save_chgtrn
        use_repuls = save_repuls
        use_solv   = save_solv
        use_list = save_list
        use_smd_velconst = save_smdvel
        use_smd_forconst = save_smdfor

        nonBonded_l = .true.

      elseif(do_int .and.(.not. do_long)) then
        if(.not.compute_polar) call attach_mpolar_pointer

        use_crec   = save_crec
        use_mrec   = save_mrec
        use_prec   = save_prec
        use_disprec= save_disprec
        use_cself  = .true.
        use_mself  = .true.
        use_pself  = .true.
        use_dispself       = .true.
        use_cshortreal     = .false.
        use_mpoleshortreal = .false.
        use_vdwshort       = .false.
        use_polarshortreal = .false.
        use_repulsshort    = .false.
        use_dispshort      = .false.
        use_dispshortreal  = .false.
        use_chgtrnshort    = .false.
        polalg     = save_polalg
        tcgorder   = save_tcgorder
        tcgprec    = save_tcgprec
        tcgguess   = save_tcgguess
        tcgpeek    = save_tcgpeek
        use_smd_velconst = save_smdvel
        use_smd_forconst = save_smdfor
        use_polar = save_polar
      endif

      if (deb_Path) then
!$acc wait
        write(*,*) '   << compute_energy_beads'
      endif

      end subroutine compute_energy_beads

c----------------------------------------------------------------------- 
c     BAROSTAT SUBROUTINES

      subroutine pmontepi (polymer,epot,temp)
      use atmlst
      use atmtyp
      use atomsMirror
      use bath
      use boxes
      use domdec
      use energi
      use group
      use math
      use mdstuf
      use molcul
      use moldyn
      use inform     ,only: deb_Path,deb_Energy,deb_Force,verbose
      use random_mod
      use units
      use usage
      use sizes      ,only: tinkerdebug
      use mpi
      use beads
      implicit none
      type(POLYMER_COMM_TYPE),intent(inout) :: polymer
      real(r_p), intent(inout) :: epot,temp
      integer i,j,k,iglob,ierr
      integer start,stop
      real(r_p) term
      real(r_p) kt,expterm
      real(r_p) third,weigh
      real(r_p) step
      real(r_p) scale
      real(r_p) eold
      real(r_p) rnd6
      real(r_p) xcm,ycm,zcm
      real(r_p) vxcm,vycm,vzcm
      real(r_p) volold,cosine
      real(r_p) dpot,dpv,dkin
      real(r_p) xmove,ymove,zmove
      real(r_p) vxmove,vymove,vzmove
      real(r_p) xboxold,yboxold,zboxold
      real(t_p) alphaold,betaold,gammaold
      real(t_p) temp3(3,3)
      real(t_p) hbox(3,3)
      real(t_p) ascale(3,3)
      real(r_p), allocatable :: posold(:,:)
      real(r_p), allocatable :: vold(:,:)
      real(t_p) valrand
      integer nbeadsproc,ibeadbeg,ibeadend,ibead 
      logical dotrial
      logical isotropic
      parameter(third = 1.0_re_p / 3.0_re_p)
      interface 
        function energy ()
        real(r_p) energy
        end function
      end interface
c
c     decide whether to attempt a box size change at this step
c
      dotrial = .false.
      if (ranktot.eq.0) then
        valrand = random()
      end if
      call MPI_BCAST(valrand,1,MPI_TPREC,0,MPI_COMM_WORLD,ierr)
c
      if (valrand .lt. 1.0_ti_p/real(voltrial,t_p)) dotrial=.true.
c
c     set constants and decide on type of trial box size change
c
      if (.not.dotrial) return

      ibeadbeg=polymer%ibead_beg(rank_polymer+1)
      ibeadend=polymer%ibead_end(rank_polymer+1)
      nbeadsproc = ibeadend-ibeadbeg+1

      if (deb_Path) write(*,'(A,2x,2F12.6)') 
     &     ' __montecarlo barostat__'
     &      ,valrand,1.0_ti_p/real(voltrial,t_p)
        
      isotropic = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (posold(3,n))
      allocate (vold(3,n))
!$acc data create(posold,vold,eold,dpot)
!$acc&     present(x,y,z,v,glob,mass,molmass,kmol,use,imol
!$acc&            ,epot,temp) async

!$acc update host(temp) async

      eold = 0.d0
      do k=1,nbeadsproc
        call load_bead(polymer, k+ibeadbeg-1)
        eold = eold +  energy()
      enddo
      eold = eold/real(polymer%nbeads,r_p)

!$acc wait
      call MPI_ALLREDUCE(MPI_IN_PLACE,eold,1,MPI_RPREC,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)

      kt = gasconst * temp
      if (isothermal)  kt = gasconst * kelvin
c
c     save the system state prior to trial box size change
c
      xboxold  = xbox
      yboxold  = ybox
      zboxold  = zbox
      alphaold = alpha
      betaold  = beta
      gammaold = gamma
      volold   = volbox
      !eold     = epot
c
c     for the isotropic case, change the lattice lengths uniformly
c
      if (rank.eq.0) then
        valrand = random()
      end if
      call MPI_BCAST(valrand,1,MPI_TPREC,0,MPI_COMM_WORLD,ierr)
      step   = volmove * (2.0_re_p*valrand-1.0_re_p)
      volbox = volbox + step
      scale  = (volbox/volold)**third
      xbox   = xbox * scale
      ybox   = ybox * scale
      zbox   = zbox * scale
           print*," _scale_ ",scale
!$acc update device(volbox,xbox,ybox,zbox) async
      call lattice
!$acc parallel loop async collapse(2)
      do i = 1, nbloc; do j=1,3
        iglob = glob(i)
        posold(j,iglob) = polymer%eigpos(j,iglob,1)
        vold(j,iglob) = polymer%eigvel(j,iglob,1)
        polymer%eigpos(j,iglob,1)=polymer%eigpos(j,iglob,1)*scale
        polymer%eigvel(j,iglob,1)=polymer%eigvel(j,iglob,1)/scale
      end do; enddo;
!$acc parallel loop async collapse(3)
      do k=1, nbeads; do i = 1, nbloc; do j=1,3
        iglob = glob(i)
        ibead = ibeadbeg + k - 1
        polymer%pos(j,iglob,ibead) = polymer%pos(j,iglob,ibead) 
     &     + (polymer%eigpos(j,iglob,1) - posold(j,iglob))
      end do; enddo; enddo
c             call ddpme3dnpt(scale,0)
c
c     get the potential energy and PV work changes for trial move
c
      epot = 0.d0
      do k=1,nbeadsproc
        call load_bead(polymer, k+ibeadbeg-1)
        epot = epot +  energy()
      enddo
      epot = epot/real(polymer%nbeads,r_p)

      if(nproctot>1)
     $   call MPI_ALLREDUCE(MPI_IN_PLACE,epot,1,MPI_RPREC,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)

      dpot = epot - eold
      dpv = atmsph * (volbox-volold) / prescon
      write(*,'(A,3f12.2)') "dpot=",dpot,epot,eold
c
c     estimate the kinetic energy change as an ideal gas term
c
      dkin = real(nuse,r_p) * kt * log(volold/volbox)
c
c     compute the kinetic energy change from the actual velocities;
c     scale the kinetic energy change to match virial pressure
c
c        dkin = 0.0_re_p
c        do i = 1, n
c           if (use(i)) then
c              term = 0.5_re_p * mass(i) / convert
c              do j = 1, 3
c                 dkin = dkin + term*(v(j,i)**2-vold(j,i)**2)
c              end do
c           end if
c        end do
c        dkin = 0.907_re_p * dkin
c
c     acceptance ratio from Epot change, Ekin change and PV work
c
      term = -(dpot+dpv+dkin) / kt
      expterm = exp(term)
      write(*,'(A,4f10.4)') "expected acc ratio:"
     &   ,expterm, exp(-dpot/kt), exp(-dpv/kt), exp(-dkin/kt)
c
c     reject the step, and restore values prior to trial change
c
 15      format(A,F20.6,D24.12,4F20.6)
      if (ranktot.eq.0) then
        valrand = random()
      end if
      call MPI_BCAST(valrand,1,MPI_TPREC,0,MPI_COMM_WORLD,ierr)
      if (valrand .gt. expterm) then
          if (ranktot.eq.0.and.tinkerdebug.gt.0)
     &      write(*,15) ' Reject montecarlo',valrand
     &                 ,expterm,epot,eold,dpot,dkin
        epot = eold
        xbox = xboxold
        ybox = yboxold
        zbox = zboxold
!$acc update device(xbox,ybox,zbox,epot) async
        call lattice
!$acc parallel loop async collapse(3)
        do k=1, nbeads; do i = 1, nbloc; do j=1,3
          iglob = glob(i)
          ibead = ibeadbeg + k - 1
          polymer%pos(j,iglob,ibead) = polymer%pos(j,iglob,ibead) 
     &     - (polymer%eigpos(j,iglob,1) - posold(j,iglob))
        end do; enddo; enddo
!$acc parallel loop async collapse(2)
        do i = 1, nbloc; do j=1,3
          iglob = glob(i)
          polymer%eigpos(j,iglob,1)=posold(j,iglob)
          polymer%eigvel(j,iglob,1)=vold(j,iglob)
        end do; enddo;
      else
!$acc update device(epot) async
        if (ranktot.eq.0.and.tinkerdebug.gt.0)
     &      write(*,15) ' Accept montecarlo',valrand
     &                 ,expterm,epot,eold,dpot,dkin
        if (ranktot.eq.0.and.verbose.and.tinkerdebug.eq.0)
     &      write(*,*) 'Applied montecarlo barostat'
        !rescale domain decomposition related stuff
        call ddpme3dnpt(scale,0)
      end if
c
c     perform deallocation of some local arrays
c
!$acc end data 
      deallocate (posold)
      deallocate (vold)
      end subroutine pmontepi

      subroutine plangevinpi(polymer,dt,istep,bloc)
      use boxes
      use bath
      use beads
      use langevin
      use domdec
      use mpi
      use random_mod
      use units
      use inform, only: deb_Path
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      REAL(r_p), intent(in) :: dt
      integer, intent(in) :: istep
      character, intent(in) :: bloc
      integer :: iloc,ierr,i,j
      real(r_p) :: dt2,factor,a1p,a2p,tmp,scale(3)
      real(r_p) :: scalex,scaley,scalez

      if(deb_path) then
!$acc wait
        write(*,*) '  >> plangevinpi block '//bloc
      endif

      SELECT CASE(bloc)
      CASE('B')
        if(ranktot>0) return
!$acc update host(presvir,stresspi) async
!$acc wait
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
          call MPI_BCAST(scale,3,MPI_RPREC,0,MPI_COMM_WORLD,ierr)
        else
          if(ranktot==0) tmp=exp(dt*vextvol)
          call MPI_BCAST(tmp,1,MPI_RPREC,0,MPI_COMM_WORLD,ierr)
          scale(:)=tmp
        endif
        extvol = extvol*scale(1)*scale(2)*scale(3)
        call rescale_box(istep,scale)
        scalex=scale(1); scaley=scale(2); scalez=scale(3)
!$acc parallel loop async
!$acc& present(polymer%eigpos,polymer%eigvel)
        DO iloc=1,nloc
          i=glob(iloc)
          polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)*scalex
          polymer%eigvel(1,i,1)=polymer%eigvel(1,i,1)/scalex
          polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)*scaley
          polymer%eigvel(2,i,1)=polymer%eigvel(2,i,1)/scaley
          polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)*scalez
          polymer%eigvel(3,i,1)=polymer%eigvel(3,i,1)/scalez
          !polymer%eigforces(j,i,1)=polymer%eigforces(j,i,1)/scale
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
      use atomsMirror
      use bath
      use boxes
      use domdec
      use math
      use inform    ,only:deb_Path
      use tinheader ,only:ti_p,re_p
      use usage
      use beads
      use mpi
      implicit none
      TYPE(POLYMER_COMM_TYPE), intent(inout) :: polymer
      REAL(r_p), intent(in) :: dt
      integer, intent(in) :: istep
      integer i,j,k,iglob,iloc,ierr
      integer start,stop
      real(r_p) pres
      real(t_p) weigh,cosine
      real(r_p) scale,third
      real(r_p) xcm,xmove
      real(r_p) ycm,ymove
      real(r_p) zcm,zmove
      real(r_p) stress(3,3)
      real(t_p) temp(3,3)
      real(t_p) hbox(3,3)
      real(t_p) ascale(3,3)
      parameter(third = 1.0_re_p / 3.0_re_p)
c
c
c     find the isotropic scale factor for constant pressure
c
!$acc update host(presvir) async
!$acc wait
      if(nproctot>1) then
        call MPI_BCAST(presvir,1,MPI_RPREC,0,MPI_COMM_WORLD,ierr)
      endif
      scale = (1.0_re_p + (dt*compress/taupres)
     &                    *(presvir-atmsph))**third
c
c     modify the current periodic box dimension values
c
         xbox = xbox * scale
         ybox = ybox * scale
         zbox = zbox * scale
!$acc update device(xbox,ybox,zbox) async
c
c     propagate the new box dimensions to other lattice values
c
         call lattice

         if (deb_Path) then
 13      format(A,8F16.6)
         print 13,'pscale',scale,dt,compress,pres,atmsph,xbox,ybox,zbox
         end if
c
c     couple to pressure bath via atom scaling in Cartesian space
c
!$acc parallel loop collapse(2) async
!$acc& present(polymer%eigpos,polymer%eigvel)
        DO iloc=1,nloc  ; DO j=1,3
          i=glob(iloc)
          polymer%eigpos(j,i,1)=polymer%eigpos(j,i,1)*scale
          polymer%eigvel(j,i,1)=polymer%eigvel(j,i,1)/scale
        ENDDO ; ENDDO 
c
c   also rescale xbegproc, xendproc...
c
         call ddpme3dnpt(scale,istep)
c
      end subroutine pscalepi

      end module utilbaoabpi

