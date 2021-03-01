c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
#include "tinker_precision.h"
      subroutine mechanic
      use atoms,only:n
      use domdec
      use deriv ,only: ConfigPotential
      use inform
      use iounit
      use interfaces
      use neigh ,only: ineignl
      use potent
      use sizes ,only: tinkerdebug
      use tinMemory ,only: print_memory_usage
      use utils     ,only: Tinker_Wait
      use vdwpot
#ifdef _OPENACC
      use random_mod
      use utilgpu,only:rec_queue,thrust_cache_init
      use openacc
#endif
      use mpi
      implicit none
c
c     set the bonded connectivity lists and active atoms
c
      call kewald
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds(.true.)
      call angles(.true.)
      call torsions(.true.)
      call bitors(.true.)
      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c      call polymer
c
c     assign atom types, classes and other atomic information
c
      call katom
c
c     assign atoms to molecules and set the atom groups
c
      call molecule(.true.)
#ifdef _OPENACC
c
c     initiate Thrust cache memory
c
      call thrust_cache_init(n)
#endif
      if (nproc.eq.1) then
!$acc update host(ineignl) async
      end if
!$acc wait
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
c
c     assign bond, angle and cross term potential parameters
c
      call clean_shared_space(0)
c      call orbital
      if (use_bond .or. use_strbnd .or. use_strtor
     $    .or. (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))
     $                 call kbond
      if (use_angle .or. use_strbnd .or. use_angang) call kangle
      if (use_strbnd)  call kstrbnd(.true.)
      if (use_urey)    call kurey(.true.)
      if (use_angang)  call kangang(.true.)
c
c     assign out-of-plane deformation potential parameters
c
      if (use_angle .or. use_opbend)  call kopbend(.true.)
      if (use_angle .or. use_opdist)  call kopdist(.true.)
      if (use_improp)  call kimprop(.true.)
      if (use_imptor)  call kimptor(.true.)
c
c     assign torsion and torsion cross term potential parameters
c
      if (use_tors .or. use_strtor .or. use_tortor)  call ktors
      if (use_pitors)  call kpitors(.true.)
      if (use_strtor)  call kstrtor(.true.)
      if (use_tortor)  call ktortor(.true.)
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_charge)                call kcharge(.true.,0)
      if (use_vdw)                   call kvdw   (.true.,0)
      if (use_mpole .or. use_polar .or. use_solv) call kmpole(.true.,0)
      if (use_polar .or. use_mpole)  call kpolar(.true.,0)

      if (use_geom)    call kgeom(.true.) ! assign restraint parameters
c
!$acc wait
      call clean_shared_space(1)
      call initmpipmegpu
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     set holonomic constrains
c
      call shakeup(.true.)
c
c     SMD parametrization
c
      call ksmd(.true.)
c
c     aMD parametrization
c
      call kamd(.true.)
c
c     associate pointers routine
c
      call configure_routine
c
c     Set Potential
c
      call ConfigPotential
c
c     construct scaling factors
c
      call build_scaling_factor(0)
#ifdef _OPENACC
c
c     init cuda rand engine
c
      call init_curand_engine(n,rec_queue)
#endif
      if (tinkerdebug.gt.0) call print_field
c
c     quit if essential parameter information is missing
c
!$acc wait
      if (abort) then
         if (rank.eq.0) write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      end

      !TODO Delete this routine
      subroutine mHalte
      use domdec,only:rank
      use tinMemory
      use utils,only: Tinker_Wait
      implicit none
      call print_memory_usage;call sleep(1);call Tinker_Wait(9,rank)
      end subroutine

      !TODO Put this routine in an appropriate region
      subroutine print_field
      use atoms ,only:n
      use angle ,only:nangle
      use angang,only:nangang
      use bond  ,only:nbond
      use bitor ,only:nbitor
      use charge,only:nion
      use domdec
      use improp,only:niprop
      use imptor,only:nitors
      use kgeoms
      use mpole ,only:npole
      use mutant,only:nmut
      use strbnd,only:nstrbnd
      use opbend,only:nopbend
      use opdist,only:nopdist
      use potent
      use polar ,only:npolar
      use mpi
      use pitors,only:npitors
      use strtor,only:nstrtor
      use tors  ,only:ntors
      use tortor,only:ntortor
      use urey  ,only:nurey
      use vdw   ,only:nvdw
      implicit none
 12   format(1x,A15,2x,I8)
 13   format(5x,11("-"),2x,8("+"))
 14   format(1x,A15,2x,I8,3x,I8)
 15   format(1x,A15,2x,6I5)

      if (rank.eq.0) then
      if (use_bond  ) write(*,12) 'BOND'  ,nbond
      if (use_angle ) write(*,12) 'ANGLE' ,nangle
      if (use_strbnd) write(*,12) 'STRBND',nstrbnd
      if (use_urey  ) write(*,12) 'UREY'  ,nurey
      if (use_angang) write(*,12) 'ANGANG',nangang
      if (use_opbend) write(*,12) 'OPBEND',nopbend
      if (use_opdist) write(*,12) 'OPDIST',nopdist
      if (use_improp) write(*,12) 'IMPROP',niprop
      if (use_imptor) write(*,12) 'IMPTOR',nitors
      if (use_tors  ) write(*,12) 'TORS'  ,ntors
      if(nbitor.ne.0) write(*,12) 'BITORS',nbitor
      if (use_pitors) write(*,12) 'PITORS',npitors
      if (use_strtor) write(*,12) 'STRTOR',nstrtor
      if (use_tortor) write(*,12) 'TORTOR',ntortor
      if (use_charge) write(*,12) 'CHARGE',nion
      if (use_vdw   ) write(*,12) 'VDW'   ,nvdw
      if (use_mpole ) write(*,12) 'MPOLE' ,npole
      if (use_polar ) write(*,12) 'POLAR' ,npolar
      if (nmut.ne.0 ) write(*,12) 'MUTATION',nmut
      if (use_geom  ) then
         write(*,15) 'GEOM RESTRAINS',
     &               npfix,ndfix,nafix,ntfix,ngfix,nchir
         write(*,*)  '     USE_BASIN ',use_basin
         write(*,*)  '     USE_WALL  ',use_wall
      end if
      if (use_extra ) write(*,12) 'EXTRA'
      if (use_smd_forconst.or.use_smd_velconst)
     &   write(*,12) rank, 'SMD'
      if (use_gamd)  write(*,12) 'GAMD'

      print 13
      if (PotentialAll)         write(*,12) 'ALL POTENTIAL',n
      if (PotentialAmoeba)      write(*,12) 'AMOEBA',n
      if (PotentialAmoeba18)    write(*,12) 'AMOEBA18',n
      if (PotentialCharmm)      write(*,12) 'CHARMM',n
      if (PotentialWaterAmoeba) write(*,12) 'WATERAMOEBA',n
      if (PotentialWaterCharmm) write(*,12) 'WATERCHARMM',n
      end if

      end subroutine

      subroutine clean_shared_space(config)
      use domdec,only:rank
      use potent
      use utils,only: Tinker_Wait
      implicit none
      integer,intent(in):: config
      logical,save:: save_bond,save_angle
      logical,save:: save_strbnd,save_urey
      logical,save:: save_angang,save_opbend
      logical,save:: save_opdist,save_improp
      logical,save:: save_imptor,save_tors
      logical,save:: save_pitors,save_strtor
      logical,save:: save_tortor,save_geom
      logical,save:: save_metal,save_extra
      logical,save:: save_vdw,save_charge
      logical,save:: save_dipole
      logical,save:: save_mpole,save_polar
      logical,save:: save_rxnfld,save_solv
      logical,save:: save_list
      logical :: first_in=.true.
      logical passed,remove
      integer save_field,freeSpace
      parameter (save_field=0, freeSpace=1)

      if (config.eq.save_field) then

      save_vdw    = use_vdw
      save_charge = use_charge
      save_mpole  = use_mpole
      save_polar  = use_polar
      save_solv   = use_solv
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
      save_tortor = use_tortor
      save_geom   = use_geom
      save_extra  = use_extra

      else if (config.eq.freeSpace) then

      ! kangle
      passed =save_angle.or.use_strbnd.or.use_angang
      remove =.not.(use_angle .or. use_strbnd .or. use_angang)
      if (passed.and.remove) call delete_data_kangle
      ! kopbend
      passed = save_angle.or.save_opbend
      remove = .not.(use_opbend)
      if (passed.and.remove) call delete_data_kopbend
      ! ktors
      passed = save_tors.or.save_strtor.or.save_tortor
      remove = .not.(use_tors.or.use_strtor.or.use_tortor)
      if (passed.and.remove) call delete_data_ktors
      ! kmpole
      passed = save_solv.or.save_mpole.or.save_polar
      remove = .not.(use_mpole)
      if (passed.and.remove) call delete_data_mpole
      ! kpolar
      passed = save_polar.or.save_mpole
      remove = .not.(use_polar)
      if (passed.and.remove) call delete_data_polar

      else

 11   format("   WARNING !!!",/," clean_shared_space routine;",
     &   "Unknown configuration call;",/," Check source code" )
         print 11
         call Tinker_Wait(1,rank)

      end if

      end subroutine
c
c     subroutine mechanicstep : fill the array parameters between two time steps
c
      subroutine mechanicstep(istep)
      use domdec,only:nproc,rank,ndir,nrec
      use potent
      use timestat
      implicit none
      integer istep

      call timer_enter( timer_param )
      if (nproc.eq.1.or.(use_pmecore.and.ndir.eq.1.and.nrec.eq.1)) then
         call ksmd_recom
         call timer_exit( timer_param )
         return
      end if
c
!$acc wait
c     call molecule(.false.)
      call bonds   (.false.)
      call angles  (.false.)
      call torsions(.false.)
      call bitors  (.false.)
!$acc wait
      if (use_charge)  call kcharge(.false.,istep)
      if (use_mpole)   call kmpole (.false.,istep)
      if (use_polar)   call kpolar (.false.,istep)
      if (use_vdw)     call kvdw   (.false.,istep)
      if (use_strbnd)  call kstrbnd(.false.)
      if (use_urey)    call kurey  (.false.)
      if (use_angang)  call kangang(.false.)
      if (use_opbend)  call kopbend(.false.)
      if (use_opdist)  call kopdist(.false.)
      if (use_improp)  call kimprop(.false.)
      if (use_imptor)  call kimptor(.false.)
      if (use_pitors)  call kpitors(.false.)
      if (use_strtor)  call kstrtor(.false.)
      if (use_tortor)  call ktortor(.false.)
      if (use_geom)    call kgeom  (.false.)
      if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
!$acc wait
      call initmpipmegpu
      call build_scaling_factor(istep)
c
c     set holonomic constrains
c
      call shakeup(.false.)
      call timer_exit( timer_param )
      end subroutine
c
c     subroutine mechanicsteprespa: fill the array parameters between two time steps
c
      subroutine mechanicsteprespa(istep,fast)
      use domdec,only: nproc,rank,ndir,nrec
      use potent
      use timestat
      implicit none
      integer istep
      logical fast

      call timer_enter( timer_param )
      if (nproc.eq.1.or.(use_pmecore.and.ndir.eq.1.and.nrec.eq.1)) then
         if (fast) call ksmd_recom
         call timer_exit( timer_param )
         return
      end if
c
c     call molecule(.false.)
      if (fast) then
!$acc wait
         call bonds   (.false.)  ! Alloc(maxvalue*nbloc) glob(:nloc) -> bndglob(:nbondloc)
         call angles  (.false.)  ! Alloc(4*nbloc) glob(:nloc) -> angleglob(:nangleloc)
         call torsions(.false.)  ! Alloc(6*nbloc) bndglob(:nbondloc) -> torsglob(:ntorsloc)
         call bitors  (.false.)  ! Alloc(8*nbloc) angleglob(:nangleloc) -> bitorsglob(:nbitorloc)
!Wait for nbondloc, nangleloc, ntorsloc, nbitorloc
!$acc wait
c        print*,use_strbnd,use_urey,use_angang
c        print*,use_angle,use_opbend,use_opdist
c        print*,use_improp,use_imptor,use_pitors
c        print*,use_strtor,use_tortor,use_geom
         if (use_strbnd) call kstrbnd(.false.)  ! Alloc(nangleloc) angleglob(:nangleloc) -> strbndglob(:nstrbndloc)
         if (use_urey)   call kurey  (.false.)  ! Alloc(nangleloc) angleglob(:nangleloc) -> ureyglob(:nureyloc)
         if (use_angang) call kangang(.false.)  ! no ACC
         if (use_opbend) call kopbend(.false.)  ! Alloc(nangleloc) angleglob(:nangleloc) -> opbendglob(:nopbendloc)
         if (use_opdist) call kopdist(.false.)  ! no ACC
         if (use_improp) call kimprop(.false.)  ! no ACC
         if (use_imptor) call kimptor(.false.)  ! no ACC
         if (use_pitors) call kpitors(.false.)  ! Alloc(ntorsloc/!\) bndglob(:nbondloc) -> pitorsglob(npitorsloc)
         if (use_strtor) call kstrtor(.false.)  ! Alloc(ntorsloc) torsglob(:ntorsloc) -> strtorglob(:nstrtorloc) + nstrtor
         if (use_tortor) call ktortor(.false.)  ! Alloc(nbitorloc) bitorsglob(:nbitorloc) -> tortorglob(:ntortorloc)
         if (use_geom)   call kgeom  (.false.)  ! no ACC
         if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
!Wait for array lenghts nstrbndloc,nureyloc,nopbendloc,npitorsloc,nstrtorloc,ntortorloc
!$acc wait
      else
         if (use_charge) call kcharge(.false.,istep)
         if (use_mpole)  call kmpole (.false.,istep)
         if (use_polar)  call kpolar (.false.,istep)
         if (use_vdw)    call kvdw   (.false.,istep)
!$acc wait
         call initmpipmegpu
         call build_scaling_factor(istep)
c
c     set holonomic constrains
c
         call shakeup(.false.)
      end if
      call timer_exit( timer_param )
      end subroutine
c
c     subroutine mechanicsteprespa1: fill the array parameters between two time steps
c
      subroutine mechanicsteprespa1(istep,rule)
      use domdec ,only: rank,nproc,ndir,nrec
      use iounit
      use moldyn ,only: stepint,nalt
      use potent
      use timestat
      implicit none
      integer istep,rule
 1000 format(' illegal rule in mechanicsteprespa1.')
c
c     rule = 0: fast part of the forces
c     rule = 1: intermediate part of the forces
c     rule = 2: slow part of the forces
c
      call timer_enter( timer_param )
      ! Handle fusion between empole & polar routine
      if (rule.eq.1.and.use_polar) then
         if (stepint.eq.1.and.nalt.ne.1) then
            call attach_mpolar_pointer
         else if (stepint.eq.nalt) then
            call detach_mpolar_pointer
         end if
      end if

      if (nproc.eq.1.or.(use_pmecore.and.ndir.eq.1)) then
         if (rule.eq.0) call ksmd_recom
         call timer_exit( timer_param )
         return
      end if

c      call molecule(.false.)
      if (rule.eq.0) then
!$acc wait
        call bonds   (.false.)
        call angles  (.false.)
        call torsions(.false.)
        call bitors  (.false.)
!$acc wait
        if (use_strbnd) call kstrbnd(.false.)
        if (use_urey)   call kurey  (.false.)
        if (use_angang) call kangang(.false.)
        if (use_opbend) call kopbend(.false.)
        if (use_opdist) call kopdist(.false.)
        if (use_improp) call kimprop(.false.)
        if (use_imptor) call kimptor(.false.)
        if (use_pitors) call kpitors(.false.)
        if (use_strtor) call kstrtor(.false.)
        if (use_tortor) call ktortor(.false.)
        if (use_geom)   call kgeom  (.false.)
        if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
!$acc wait
      else if (rule.eq.1) then
        if (use_charge) call kcharge(.false.,istep)
        if (use_mpole)  call kmpole (.false.,istep)
        if (use_polar)  call kpolar (.false.,istep)
        if (use_vdw)    call kvdw   (.false.,istep)
!$acc wait
      else if (rule.eq.2) then
        if (use_charge) call kcharge(.false.,istep)
        if (use_mpole)  call kmpole (.false.,istep)
        if (use_polar)  call kpolar (.false.,istep)
        if (use_vdw)    call kvdw   (.false.,istep)
!$acc wait
        call build_scaling_factor(istep)
        call initmpipmegpu
c       call updategrid
c
c     set holonomic constrains
c
        call shakeup(.false.)
      else
         if (rank.eq.0) write(iout,1000) 
      end if
      call timer_exit( timer_param )
      end

      subroutine init_sub_config
      use atoms      ,only: n
      use interfaces
      use mdstuf
      use neigh      ,only: defaultlbuffer,defaultlbuffer1,lbuffer
      use utilgpu
      use sizes
      implicit none
      integer devtyp
      enum,bind(C)
        enumerator dev_quadro,dev_tesla
        enumerator dev_geforce
        enumerator dev_default
      end enum
      character(256) devname
      real(t_p) buff_t

      if (sub_config.eq.-1) then
         sub_config = itrf_legacy
#ifdef _OPENACC
         devname = devProp%name
         call upcase(devname)
         if (index(devname,'GEFORCE').gt.0) then
            devtyp = dev_geforce
         else if (index(devname,'QUADRO').gt.0) then
            devtyp = dev_quadro
         else if (index(devname,'TESLA').gt.0) then
            devtyp = dev_tesla
         else
            devtyp = dev_default
         end if
         if (devtyp.eq.dev_quadro.or.devtyp.eq.dev_tesla) then
            if ( n.gt.10000 ) then
               sub_config = itrf_adapted
#if TINKER_DOUBLE_PREC
               if (lbuffer.eq.defaultlbuffer1.or.
     &             lbuffer.eq.defaultlbuffer) then
                   if (integrate.eq.'RESPA1'.or.
     &                 integrate.eq.'BAOABRESPA1') then
                      buff_t = merge( 0.5,1.0,(n.gt.95000) )
                      call update_lbuffer(buff_t)
                   else
                      buff_t = merge( 0.7,0.4,(n.gt.95000) )
                      call update_lbuffer(buff_t)
                   end if
               end if
#endif
            else
#if TINKER_DOUBLE_PREC
               sub_config = itrf_legacy
#else
               sub_config = itrf_adapted
#endif
            end if
         else if (devtyp.eq.dev_geforce) then
#if TINKER_DOUBLE_PREC
            sub_config = itrf_legacy
#else
            sub_config = itrf_adapted
#endif
         else
            sub_config = itrf_adapted
         end if
#endif
         if (tinkerdebug.gt.0.and.rank.eq.0) then
 12         format( /,'  --- run mode :  LEGACY  ',/ )
 13         format( /,'  --- run mode :  ADAPTED ',/ )
            if      (sub_config.eq.itrf_legacy) then
               write(*,'(5x,a)') trim(devname)
               write(*,12)
            else if (sub_config.eq.itrf_adapted) then
               write(*,'(5x,a)') trim(devname)
               write(*,13)
            end if
         end if
      end if
      end subroutine

      subroutine configure_routine
      use interfaces
      use domdec ,only: nproc,nrec
      use precompute_pole,only: tmatxb_pme_compute,
     &    precompute_solvpole,precompute_mpole
      use potent, only: use_polar,use_mpole,use_vdw
     &          , use_charge,use_pmecore
      use polar , only: use_mpolar_ker
      use neigh , only: vlst_enable,vlst2_enable
     &          , mlst_enable,mlst2_enable
     &          , clst_enable,clst2_enable
      use sizes , only: tinkerdebug
      use vdwpot, only: vdwtyp
      implicit none
      integer configure,shft
      integer,parameter:: hexa_len=16,sh_len=4
#if _OPENACC
      procedure(tmatxb_pme_core1):: tmatxb_pme_core_v4
#endif

      call init_sub_config
      call init_routine_pointers

      vlst_enable  = .false.
      vlst2_enable = .false.
      mlst_enable  = .false.
      mlst2_enable = .false.
      clst_enable  = .false.
      clst2_enable = .false.

      shft = ishft(sub_config,-conf_vlist*sh_len)
      configure=iand(shft,hexa_len-1)
      if (btest(configure,list_verlet)) vlst_enable = .true.
#ifdef _OPENACC
      if (btest(configure,list_block_atoms)) vlst2_enable = .true.
#endif

      shft = ishft(sub_config,-conf_mlist*sh_len)
      configure=iand(shft,hexa_len-1)
      if (btest(configure,list_verlet)) mlst_enable = .true.
#ifdef _OPENACC
      if (btest(configure,list_block_atoms)) mlst2_enable = .true.
#endif

      if (use_vdw) then

         if (vdwtyp.eq.'LENNARD-JONES') then

         shft     = ishft(sub_config,-conf_elj*sh_len)
         configure=  iand(shft,hexa_len-1)
         eljsl1c_p => eljshortlong1cgpu
         elj1c_p    => elj1cgpu
         elj3c_p    => elj3cgpu
         if      (configure.eq.conf_elj) then
            eljsl1c_p  => eljshortlong1cgpu
            elj1c_p    => elj1c
            elj3c_p    => elj3c
            vlst_enable =.true.
         else if (configure.eq.conf_elj1gpu) then
            eljsl1c_p  => eljshortlong1cgpu
            elj1c_p    => elj1cgpu
            elj3c_p    => elj3cgpu
            vlst_enable =.true.
#ifdef _OPENACC
         else if (configure.eq.conf_elj1cu) then
            eljsl1c_p   => eljshortlong1cgpu
            elj1c_p     => elj1c_cu
            elj3c_p     => elj3c_cu
            vlst2_enable =.true.
#endif
         else
            vlst_enable =.true.
         end if

         end if

         if (vdwtyp.eq.'BUFFERED-14-7') then

         shft     = ishft(sub_config,-conf_ehal*sh_len)
         configure=iand(shft,hexa_len-1)
         ehalshort1c_p => ehalshort1cgpu
         ehallong1c_p  => ehallong1cgpu
         ehal3c_p  => ehal3cgpu
         ehalshortlong3c_p => ehalshortlong3cgpu
         if      (configure.eq.conf_ehal1) then
            ehal1c_p => ehal1cgpu1
            vlst_enable =.true.
         else if (configure.eq.conf_ehal2) then
            ehal1c_p => ehal1cgpu2
            vlst_enable =.true.
#ifdef _OPENACC
         else if (configure.eq.conf_ehalcu) then
            ehal1c_p => ehal1c_cu
            ehalshort1c_p => ehalshortlong1c_cu
            ehallong1c_p  => ehalshortlong1c_cu
            ehal3c_p  => ehal3c_cu
            ehalshortlong3c_p => ehalshortlong3c_cu
            vlst2_enable =.true.
#endif
         else
            ehal1c_p => ehal1cgpu2
            vlst_enable =.true.
         end if

         end if

      end if

      if (use_charge) then
         shft     = ishft(sub_config,-conf_echg*sh_len)
         configure=  iand(shft,hexa_len-1)
         if      (configure.eq.conf_ecreal1d)    then
            ecreal1d_p          => ecreal1d
            ecrealshortlong1d_p => ecrealshortlong1dgpu
            ecreal3d_p          => ecreal3dgpu
            clst_enable = .true.
         else if (configure.eq.conf_ecreal1dgpu) then
            ecreal1d_p          => ecreal1dgpu
            ecrealshortlong1d_p => ecrealshortlong1dgpu
            ecreal3d_p          => ecreal3dgpu
            clst_enable = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_ecreal1d_cu) then
            ecreal1d_p          => ecreal1d_cu
            ecrealshortlong1d_p => ecrealshortlong1dgpu
            ecreal3d_p          => ecreal3d_cu
            clst2_enable = .true.
#endif
         else
            ecreal1d_p          => ecreal1dgpu
            ecrealshortlong1d_p => ecrealshortlong1dgpu
            ecreal3d_p          => ecreal3dgpu
            clst_enable = .true.
         end if
      end if

      if (use_mpole) then
         shft = ishft(sub_config,-conf_mpole*sh_len)
         configure=iand(shft,hexa_len-1)
         emreal1c_p     => emreal1cgpu
         emrealshort1c_p=> emrealshort1cgpu
         emreallong1c_p => emreallong1cgpu
         emrealshortlong1c_core_p => emrealshortlong1c_core
         emrealshortlong3d_p      => emrealshortlong3d
         if      (configure.eq.conf_emreal1c_1) then
            emreal1c_core_p => emreal1c_core1
            emreal3d_p => emreal3dgpu
            mlst_enable = .true.
         else if (configure.eq.conf_emreal1c_2) then
            emreal1c_core_p => emreal1c_core2
            emreal3d_p => emreal3dgpu
            mlst_enable = .true.
         else if (configure.eq.conf_emreal1c_pre.or.precompute_mpole)
     &        then
            emreal1c_core_p => emreal1c_core3
            emreal3d_p => emreal3dgpu
            mlst_enable = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_emreal1c_4) then
            emreal1c_core_p => emreal1c_core4
            emrealshortlong1c_core_p => emrealshortlong1c_core2
            emreal3d_p => emreal3d_cu
            emrealshortlong3d_p => emrealshortlong3d_cu
            mlst2_enable = .true.
         else if (configure.eq.conf_emreal1c_5) then
            emreal1c_core_p => emreal1c_core5
            emrealshortlong1c_core_p => emrealshortlong1c_core2
            emreal3d_p => emreal3d_cu
            emrealshortlong3d_p => emrealshortlong3d_cu
            mlst2_enable = .true.
#endif
         else
            emreal1c_core_p => emreal1c_core2
            emreal3d_p => emreal3dgpu
            mlst_enable = .true.
         end if
      end if

      if (use_polar) then
         shft = ishft(sub_config,-conf_efld0*sh_len)
         configure=iand(shft,hexa_len-1)
         otf_dc_efld0_directgpu_p => otf_dc_efld0_directgpu2
         if      (configure.eq.conf_efld0_directgpu1) then
            efld0_directgpu_p => efld0_directgpu
            mlst_enable = .true.
         else if (configure.eq.conf_efld0_directgpu2) then
            efld0_directgpu_p => efld0_directgpu2
            mlst_enable = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_efld0_directgpu3) then
            efld0_directgpu_p => efld0_directgpu3
        otf_dc_efld0_directgpu_p => otf_dc_efld0_directgpu3
            mlst2_enable = .true.
#endif
         else
            efld0_directgpu_p => efld0_directgpu2
         otf_dc_efld0_directgpu_p => otf_dc_efld0_directgpu2
            mlst_enable = .true.
         end if
         efld0_directgpu_p1  => efld0_directgpu_p

         shft = ishft(sub_config,-conf_tmat*sh_len)
         configure=iand(shft,hexa_len-1)
         tmatxb_p              => tmatxb_pmegpu
         otf_dc_tmatxb_pme_core_p => otf_dc_tmatxb_pme_core2
         if      (configure.eq.conf_tmatxb_ref) then
            nullify(tmatxb_p)
            tmatxb_p          => tmatxb_pmegpu1
            mlst_enable = .true.
         else if (configure.eq.conf_tmatxb_pre.or.
     &            precompute_solvpole) then
            tmatxb_pme_core_p => tmatxb_pme_compute
            precompute_solvpole = .true.
            mlst_enable = .true.
         else if (configure.eq.conf_tmatxb_c1) then
            tmatxb_pme_core_p => tmatxb_pme_core1
            mlst_enable = .true.
         else if (configure.eq.conf_tmatxb_c2) then
            tmatxb_pme_core_p => tmatxb_pme_core2
            mlst_enable = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_tmatxb_c3) then
              tmatxb_pme_core_p => tmatxb_pme_core3
       otf_dc_tmatxb_pme_core_p => otf_dc_tmatxb_pme_core3
            mlst2_enable = .true.
         else if (configure.eq.conf_tmatxb_v4) then
              tmatxb_pme_core_p => tmatxb_pme_core_v4
       otf_dc_tmatxb_pme_core_p => otf_dc_tmatxb_pme_core3
            mlst2_enable = .true.
#endif
         else
            tmatxb_pme_core_p => tmatxb_pme_core2
            mlst_enable = .true.
         end if
         tmatxb_p1  => tmatxb_p

         shft = ishft(sub_config,-conf_polar*sh_len)
         configure=iand(shft,hexa_len-1)
         epreal1c_p => epreal1cgpu
         if      (configure.eq.conf_epreal1c_1) then
            epreal1c_core_p => epreal1c_core1
            epreal3d_p      => epreal3dgpu
            mlst_enable     = .true.
         else if (configure.eq.conf_epreal1c_2) then
            epreal1c_core_p => epreal1c_core2
            epreal3d_p      => epreal3dgpu
            mlst_enable     = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_epreal1c_3) then
            epreal1c_core_p => epreal1c_core3
            epreal3d_p      => epreal3d_cu
            mlst2_enable    = .true.
            use_mpolar_ker  = .true.
            call attach_mpolar_pointer
#endif
         else
            epreal1c_core_p => epreal1c_core2
            epreal3d_p      => epreal3dgpu
            mlst_enable     = .true.
         end if
      end if

      shft = ishft(sub_config,-conf_fphi*sh_len)
      configure=iand(shft,hexa_len-1)
      if      (configure.eq.conf_fphi_acc) then
         fphi_uind_site1_p => fphi_uind_sitegpu1
         fphi_uind_site2_p => fphi_uind_sitegpu2
         fphi_mpole_site_p => fphi_mpole_sitegpu
         grid_pchg_force_p => grid_pchg_force
#ifdef _OPENACC
      else if (configure.eq.conf_fphi_cuda) then
         fphi_uind_site1_p => fphi_uind_sitecu1
         fphi_uind_site2_p => fphi_uind_sitecu2
         fphi_mpole_site_p => fphi_mpole_sitecu
         grid_pchg_force_p => grid_pchg_forcecu
#endif
      else
         fphi_uind_site1_p => fphi_uind_sitegpu1
         fphi_uind_site2_p => fphi_uind_sitegpu2
         fphi_mpole_site_p => fphi_mpole_sitegpu
         grid_pchg_force_p => grid_pchg_force
      end if

      shft = ishft(sub_config,-conf_grid*sh_len)
      configure=iand(shft,hexa_len-1)
      if      (configure.eq.conf_grid_acc) then
         grid_uind_site_p  => grid_uind_sitegpu
         grid_mpole_site_p => grid_mpole_sitegpu
         grid_pchg_site_p  => grid_pchg_sitegpu
#ifdef _OPENACC
      else if (configure.eq.conf_grid_cuda) then
         if (nproc.eq.1) then
            grid_uind_site_p  => grid_uind_sitecu
            grid_mpole_site_p => grid_mpole_sitecu
            grid_pchg_site_p  => grid_pchg_sitecu
         else if (use_pmecore.and.nrec.eq.1) then
            grid_uind_site_p  => grid_uind_sitecu
            grid_mpole_site_p => grid_mpole_sitecu
            grid_pchg_site_p  => grid_pchg_sitecu
         else
            grid_uind_site_p  => grid_uind_sitegpu
            grid_mpole_site_p => grid_mpole_sitegpu
            grid_pchg_site_p  => grid_pchg_sitegpu
         end if
#endif
      else
         grid_uind_site_p  => grid_uind_sitegpu
         grid_mpole_site_p => grid_mpole_sitegpu
         grid_pchg_site_p  => grid_pchg_sitegpu
      end if

c     print '(/,A21,Z16,I30)', 'routine config number ', 
c    &      sub_config,sub_config

      call config_cover_routines
      call disable_commreal_study
      end

c
c     Associate routine pointers related to 
c           MPI Async Overlapping computation (Real V Rec)
c
      subroutine config_cover_routines
      use interfaces
      use utilgpu
      use inform
      use polpot ,only: polalg
      use potent ,only: use_pmecore
      implicit none
      logical,parameter:: OK=.true.

      ! Association
      if (dir_queue.ne.rec_queue) then

         if (use_pmecore) then
            if (verbose.and.rank.eq.0) then
 13            format(5x,'--- Disabling Async Comput ---',/,
     &         5x,'Not Suited with pme-proc option')
               print 13
            end if
#ifdef _OPENACC
            call finalize_async_recover
#endif
            return
         end if

         if (nproc.eq.1.and.OK) then
 14         format(5x,'--- Disabling Async Comput ---',/,
     &             2x,'Affect performances during sequential run')
            print 14
#ifdef _OPENACC
            call finalize_async_recover
#endif
            return
         end if

         emreal1c_cp     => emreal1c_p
         emreallong1c_cp => emreallong1c_p
         emreal1c_p      => tinker_void_sub
         emreallong1c_p  => tinker_void_sub

         efld0_direct_cp   => efld0_directgpu_p
         efld0_directgpu_p => efld0_direct_void

         if (polalg.eq.1) then
            tmatxb_cp    => tmatxb_p
            tmatxb_p     => tmatxb_void
         end if

         epreal1c_cp     => epreal1c_p
         epreal1c_p      => tinker_void_sub

      end if

      end subroutine

      subroutine attach_mpolar_pointer
      use interfaces
      use polar ,only: use_mpolar_ker
      implicit none

#ifdef _OPENACC
      if (use_mpolar_ker) then
      epreal1c_core_p => mpreal1c_core
      emreal1c_p      => tinker_void_sub
      emrealshort1c_p => tinker_void_sub
      end if
#endif

      end subroutine

      subroutine detach_mpolar_pointer
      use interfaces
      use polar ,only: use_mpolar_ker
      implicit none
      integer shft,configure
      integer,parameter:: hexa_len=16,sh_len=4

#ifdef _OPENACC
      if (use_mpolar_ker.and.
     &   associated(epreal1c_core_p,mpreal1c_core)) 
     &   then
         shft = ishft(sub_config,-conf_polar*sh_len)
         configure=iand(shft,hexa_len-1)
         if      (configure.eq.conf_epreal1c_1) then
            epreal1c_core_p => epreal1c_core1
         else if (configure.eq.conf_epreal1c_2) then
            epreal1c_core_p => epreal1c_core2
         else if (configure.eq.conf_epreal1c_3) then
            epreal1c_core_p => epreal1c_core3
         else
            epreal1c_core_p => epreal1c_core2
         end if
         emreal1c_p      => emreal1cgpu
         emrealshort1c_p => emrealshort1cgpu
      end if
#endif
      end subroutine

      subroutine disable_commreal_study
      use domdec  ,only: nproc,rank
      use interfaces
      use inform  ,only: verbose
      use potent  ,only: use_mpole,use_polar
      use utilcomm,only: no_commdir
#ifdef _OPENACC
      use utilcu  ,only: cu_update_balanced_comput
#endif
      implicit none
      logical,parameter::accept=.false.

      if (nproc.gt.6.and.accept.and.
     &   (associated(efld0_directgpu_p,efld0_directgpu3).or.
     &    associated(efld0_direct_cp,efld0_directgpu3)).and.
     &    associated(tmatxb_pme_core_p,tmatxb_pme_core3))
     &   then
         no_commdir=.true.
 12      format(/,
     &   ' *****  Bypassing most of real space communications *****')
         if (rank.eq.0.and.verbose) write(*,12)

         ! Reset poleglobnl for each process
         if (use_mpole)  call kmpole (.false.,0)
         if (use_polar)  call kpolar (.false.,0)
      end if
#ifdef _OPENACC
      call cu_update_balanced_comput(.not.no_commdir)
#endif
      end subroutine
