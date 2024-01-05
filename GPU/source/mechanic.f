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
#include "tinker_macro.h"
      subroutine mechanic
      use atoms,only:n
      use cutoff 
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
      use polpot
      implicit none
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds   (.true.)
      call angles  (.true.)
      call torsions(.true.)
      call bitors  (.true.)
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
c      call orbital
c
c     assign electrostatic and dispersion Ewald sum parameters
c
      if (.not.use_ani_only) call kewald
c
      call clean_shared_space(0)
c
c     assign bond, angle and cross term potential parameters
c
      call kbond
      call kmlpot(.true.,0)
      if (use_angle .or. use_strbnd .or. use_angang)
     &                 call kangle
      if (use_strbnd)  call kstrbnd(.true.)
      if (use_urey  )  call kurey  (.true.)
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
      if (use_tors .or. use_strtor .or. use_tortor)
     &                 call ktors
      if (use_pitors)  call kpitors(.true.)
      if (use_strtor)  call kstrtor(.true.)
      if (use_angtor)  call kangtor(.true.)
      if (use_tortor)  call ktortor(.true.)
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_charge)  call kcharge(.true.,0)
      if (use_vdw   )  call kvdw   (.true.,0)
      if (use_mpole .or. use_polar .or. use_chgtrn .or. use_solv)
     &                 call kmpole (.true.,0)
      if (use_polar .or. use_chgtrn)
     &                 call kpolar (.true.,0)

      if (use_chgtrn)  call kchgtrn(.true.,0)
      if (use_mpole.or.use_polar.or.use_chgtrn) then
         if (nproc.eq.1) then
                       call kpolar_reset1(0)
         else
                       call kpolar_reset (0)
         end if
      end if
      if (use_chgflx)  call kchgflx
c
c     assign restraint parameters
c
      if (use_geom  )  call kgeom  (.true.)
c
c     assign repulsion and dispersion parameters
c
      if (use_repuls)  call krepel 
      if (use_disp  )  call kdisp  (.true.,0) 
c
!$acc wait
      call clean_shared_space(1)
      call dev_sort_global_buffers(0)
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
#ifdef _OPENACC
c
c     init cuda rand engine
c
      call init_curand_engine(n,rec_queue)
#endif
c
c     associate pointers routine
c
      call configure_routine
c
c     construct scaling factors
c
      call build_scaling_factor(0)
c
c     Set Potential
c
      call ConfigPotential
c
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
      use angtor,only:nangtor
      use bond  ,only:nbond
      use bitor ,only:nbitor
      use cflux ,only:naflx,nbflx
      use chgtrn,only:nct
      use chgpen,only:ncp
      use disp  ,only:ndisp
      use charge,only:nion
      use domdec,only:ranktot
      use group ,only:ngrp,use_group
      use improp,only:niprop
      use imptor,only:nitors
      use kgeoms,only:npfix,ndfix,nafix,ntfix,ngfix,nchir
     &          ,use_basin,use_wall
      use mpole ,only:npole
      use mutant,only:nmut
      use strbnd,only:nstrbnd
      use opbend,only:nopbend
      use opdist,only:nopdist
      use polar ,only:npolar
      use polpot,only:use_thole,use_dirdamp
      use potent
      use pitors,only:npitors
      use repel ,only:nrep
      use strtor,only:nstrtor
      use tors  ,only:ntors
      use tortor,only:ntortor
      use urey  ,only:nurey
      use vdw   ,only:nvdw
      implicit none
      character*10 prop
      logical u_mpole,u_mpolec,u_polar,u_polarc
 12   format(1x,A15,2x,I8)
 13   format(5x,11("-"),2x,8("+"))
 14   format(1x,A15,2x,I8,3x,I8)
 17   format(1x,A15,2x,I8,2x,A)

      if (ranktot.eq.0) then
         prop = merge(' DIRDAMP','        ',use_dirdamp)
         u_mpole  = use_mpole.and..not.use_chgpen
         u_mpolec = use_mpole.and.use_chgpen
         u_polar  = use_polar.and.use_thole
         u_polarc = use_polar.and.use_chgpen.and..not.use_dirdamp

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
      if (use_angtor) write(*,12) 'ANGTOR',nangtor
      if (use_tortor) write(*,12) 'TORTOR',ntortor
      if (use_charge) write(*,12) 'CHARGE',nion
      if (use_chgtrn) write(*,12) 'CHGTRN',nct
      if (use_chgflx) write(*,14) 'CHGFLX',naflx,nbflx
      if (use_vdw   ) write(*,12) 'VDW'   ,nvdw
      if (use_disp  ) write(*,12) 'DISP'  ,ndisp
      if (use_repuls) write(*,12) 'REPULS',nrep
      if (u_mpole   ) write(*,12) 'MPOLE' ,npole
      if (u_mpolec  ) write(*,12) 'MPOLE CHGPEN',npole
      if (u_polar   ) write(*,17) 'POLAR THOLE ',npolar,prop
      if (u_polarc  ) write(*,12) 'POLAR CHGPEN',ncp
      if (nmut.ne.0 ) write(*,12) 'MUTATION',nmut
      if (use_group ) write(*,12) 'GROUP' ,ngrp
      if (use_geom  ) then
 15      format(1x,A15,2x,I7,5I5)
 16      format(1x,15x,2x,A7,5A5)
         write(*,16) 'pos','dist','angl','tors','grpd','chir'
         write(*,15) 'GEOM RESTRAINS',
     &               npfix,ndfix,nafix,ntfix,ngfix,nchir
         write(*,*)  '     USE_BASIN ',use_basin
         write(*,*)  '     USE_WALL  ',use_wall
      end if
      if (use_extra ) write(*,12) 'EXTRA'
      if (use_smd_forconst.or.use_smd_velconst)
     &   write(*,12) 'SMD'
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

      subroutine zero_field
      use atoms ,only:n
      use angle ,only:nangle,nangleloc
      use angang,only:nangang,nangangloc
      use angtor,only:nangtor,nangtorloc
      use bond  ,only:nbond,nbondloc
      use bitor ,only:nbitor,nbitorloc
      use charge,only:nion,nionloc
      use domdec
      use improp,only:niprop,niproploc
      use imptor,only:nitors,nitorsloc
      use kgeoms
      use mpole ,only:npole,npoleloc
      use mutant,only:nmut
      use strbnd,only:nstrbnd,nstrbndloc
      use opbend,only:nopbend,nopbendloc
      use opdist,only:nopdist,nopdistloc
      use potent
      use polar ,only:npolar
      use mpi
      use pitors,only:npitors,npitorsloc
      use strtor,only:nstrtor,nstrtorloc
      use tors  ,only:ntors,ntorsloc
      use tortor,only:ntortor,ntortorloc
      use urey  ,only:nurey,nureyloc
      use vdw   ,only:nvdw,nvdwloc
      implicit none

      nbond  =0;  nbondloc  =0;
      nangle =0;  nangleloc =0;
      nstrbnd=0;  nstrbndloc=0;
      nurey  =0;  nureyloc  =0;
      nangang=0;  nangangloc=0;
      nopbend=0;  nopbendloc=0;
      nopdist=0;  nopdistloc=0;
      niprop =0;  niproploc =0;
      nitors =0;  nitorsloc =0;
      ntors  =0;  ntorsloc  =0;
      nbitor =0;  nbitorloc =0;
      npitors=0;  npitorsloc=0;
      nstrtor=0;  nstrtorloc=0;
      nangtor=0;  nangtorloc=0;
      ntortor=0;  ntortorloc=0;
      nion   =0;  nionloc   =0;
      nvdw   =0;  nvdwloc   =0;
      npole  =0;  npoleloc  =0;
      npolar =0;
      nmut   =0;

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
      logical,save:: save_angtor,save_tortor,save_geom
      logical,save:: save_metal,save_extra
      logical,save:: save_vdw,save_charge
      logical,save:: save_chgtrn,save_disp,save_repuls
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
      save_disp   = use_disp
      save_repuls = use_repuls
      save_charge = use_charge
      save_mpole  = use_mpole
      save_polar  = use_polar
      save_chgtrn = use_chgtrn
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
      save_angtor = use_angtor
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
      remove = .not.(use_mpole.or.use_polar.or.use_chgtrn)
      if (passed.and.remove) call dealloc_shared_mpole
      ! kpolar
      passed = save_polar.or.save_mpole
      remove = .not.(use_polar)
      if (passed.and.remove) call dealloc_shared_polar

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
      if (use_mpole.or.use_polar.or.use_chgtrn)
     &                 call kpolar_reset(istep)
      if (use_vdw   )  call kvdw   (.false.,istep)
      if (use_disp  )  call kdisp  (.false.,istep)
      if (use_strbnd)  call kstrbnd(.false.)
      if (use_urey  )  call kurey  (.false.)
      if (use_angang)  call kangang(.false.)
      if (use_opbend)  call kopbend(.false.)
      if (use_opdist)  call kopdist(.false.)
      if (use_improp)  call kimprop(.false.)
      if (use_imptor)  call kimptor(.false.)
      if (use_pitors)  call kpitors(.false.)
      if (use_strtor)  call kstrtor(.false.)
      if (use_angtor)  call kangtor(.false.)
      if (use_tortor)  call ktortor(.false.)
      if (use_geom  )  call kgeom  (.false.)
      if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
      if (use_mlpot)   call kmlpot_reset(istep)
!$acc wait
      call dev_sort_global_buffers(istep)
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
         if (use_urey  ) call kurey  (.false.)  ! Alloc(nangleloc) angleglob(:nangleloc) -> ureyglob(:nureyloc)
         if (use_angang) call kangang(.false.)  ! no ACC
         if (use_opbend) call kopbend(.false.)  ! Alloc(nangleloc) angleglob(:nangleloc) -> opbendglob(:nopbendloc)
         if (use_opdist) call kopdist(.false.)  ! no ACC
         if (use_improp) call kimprop(.false.)
         if (use_imptor) call kimptor(.false.)
         if (use_pitors) call kpitors(.false.)  ! Alloc(ntorsloc/!\) bndglob(:nbondloc) -> pitorsglob(npitorsloc)
         if (use_strtor) call kstrtor(.false.)  ! Alloc(ntorsloc) torsglob(:ntorsloc) -> strtorglob(:nstrtorloc) + nstrtor
         if (use_angtor) call kangtor(.false.)
         if (use_tortor) call ktortor(.false.)  ! Alloc(nbitorloc) bitorsglob(:nbitorloc) -> tortorglob(:ntortorloc)
         if (use_geom  ) call kgeom  (.false.)  ! no ACC
         if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
!Wait for array lenghts nstrbndloc,nureyloc,nopbendloc,npitorsloc,nstrtorloc,ntortorloc
!$acc wait
      else
         if (use_charge) call kcharge(.false.,istep)
         if (use_mpole.or.use_polar.or.use_chgtrn)
     &                   call kpolar_reset(istep)
         if (use_vdw   ) call kvdw   (.false.,istep)
         if (use_disp  ) call kdisp  (.false.,istep)
         if (use_mlpot)  call kmlpot_reset(istep)
!$acc wait
         call dev_sort_global_buffers(istep)
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
        if (use_urey  ) call kurey  (.false.)
        if (use_angang) call kangang(.false.)
        if (use_opbend) call kopbend(.false.)
        if (use_opdist) call kopdist(.false.)
        if (use_improp) call kimprop(.false.)
        if (use_imptor) call kimptor(.false.)
        if (use_pitors) call kpitors(.false.)
        if (use_strtor) call kstrtor(.false.)
        if (use_angtor) call kangtor(.false.)
        if (use_tortor) call ktortor(.false.)
        if (use_geom  ) call kgeom  (.false.)
        if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
!$acc wait
      else if (rule.eq.1) then
        if (use_charge) call kcharge(.false.,istep)
        if (use_mpole.or.use_polar.or.use_chgtrn)
     &                  call kpolar_reset(istep)
        if (use_vdw   ) call kvdw   (.false.,istep)
!$acc wait
      else if (rule.eq.2) then
        if (use_charge) call kcharge(.false.,istep)
        if (use_mpole.or.use_polar.or.use_chgtrn)
     &                  call kpolar_reset(istep)
        if (use_vdw   ) call kvdw   (.false.,istep)
        if (use_disp  ) call kdisp  (.false.,istep)
        if (use_mlpot)  call kmlpot_reset(istep)
!$acc wait
        call dev_sort_global_buffers(istep)
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
      use potent     ,only: use_chgpen
      use polpot     ,only: use_thole, use_dirdamp
      use utilgpu
      use sizes
      use beads, only: centroid_recip
      implicit none
      integer devtyp,devtyp1
      enum,bind(C)
        enumerator dev_quadro,dev_tesla,dev_geforce
        enumerator dev_other
      end enum
      enum,bind(C)
        enumerator dev_fp32,dev_fp64
      end enum
      character(256) devname
      real(t_p) buff_t

      if (sub_config.eq.-1) then
         sub_config = itrf_legacy
         devname(:) = ""
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
            devtyp = dev_other
         end if
         if (index(devname,'100').gt.0 .or. 
     &       index(devname,'TITAN V').gt.0 ) then
            devtyp1 = dev_fp64
         else
            devtyp1 = dev_fp32
         end if

         if (devtyp1.eq.dev_fp64) then
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
                      buff_t = merge( 0.4,0.7,(n.gt.95000) )
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
         else if (devtyp1.eq.dev_fp32) then
#if TINKER_DOUBLE_PREC
            sub_config = itrf_legacy
#else
            sub_config = itrf_adapted
#endif
         else
            sub_config = itrf_adapted
         end if
         if (use_chgpen.or.use_dirdamp)
     &      sub_config = itrf_adapted
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
      use atoms  ,only: n
      use couple ,only: maxn12_
      use interfaces
      use charge ,only: nion
      use domdec ,only: nproc,nrec
      use group  ,only: use_group
      use inform ,only: tinEssai,deb_Path,app_id,dynamic_a
      use precompute_pole,only: tmatxb_pme_compute,polar_precomp
      use potent ,only: use_polar,use_mpole,use_vdw
     &           , use_charge, use_pmecore, use_lambdadyn, fuse_chglj
     &           , use_chgpen, use_chgtrn, use_disp, use_repuls
     &           , use_chgflx
      use polar  ,only: use_mpolar_ker
      use polpot ,only: use_dirdamp, use_thole
     &           , u1scale,u2scale,u3scale,u4scale
      use neigh  ,only: vlst_enable,vlst2_enable
     &           , mlst_enable,mlst2_enable
     &           , clst_enable,clst2_enable
      use sizes  ,only: tinkerdebug
#ifdef _OPENACC
      use utilgpu,only: allow_precomp
#endif
      use vdw    ,only: nvdw
      use vdwpot ,only: vdwtyp
      implicit none
      integer configure,shft
      integer,parameter:: hexa_len=16,sh_len=4
      real(t_p) puscale
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
         elj1c_p    => elj1c_ac
         elj3c_p    => elj3c_ac
         if      (configure.eq.conf_elj) then
            vlst_enable =.true.
         else if (configure.eq.conf_elj1gpu) then
            vlst_enable =.true.
#ifdef _OPENACC
         else if (configure.eq.conf_elj1cu.and.maxn12_.lt.5) then
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
         ehal1c_p => ehal1c_ac
         ehal3c_p => ehal3c_ac
         if      (configure.eq.conf_ehal1) then
            vlst_enable =.true.
         else if (configure.eq.conf_ehal2) then
            vlst_enable =.true.
#ifdef _OPENACC
         else if (configure.eq.conf_ehalcu) then
            ehal1c_p => ehal1c_cu
            ehal3c_p  => ehal3c_cu
            vlst2_enable =.true.
#endif
         else
            vlst_enable =.true.
         end if

         end if

      end if

      if (use_disp) then

         vlst_enable = .true.
         vlst2_enable = .true.

      end if

      if (use_charge) then
         shft     = ishft(sub_config,-conf_echg*sh_len)
         configure=  iand(shft,hexa_len-1)
         if      (configure.eq.conf_ecreal1d)    then
            ecreal1d_p          => ecreal1d
            ecreal3d_p          => ecreal3dgpu
            clst_enable = .true.
         else if (configure.eq.conf_ecreal1dgpu) then
            ecreal1d_p          => ecreal1dgpu
            ecreal3d_p          => ecreal3dgpu
            clst_enable = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_ecreal1d_cu) then
            ecreal1d_p          => ecreal1d_cu
            ecreal3d_p          => ecreal3d_cu
            clst2_enable = .true.

            if (vdwtyp.eq.'LENNARD-JONES'.and.maxn12_.lt.5
     &         .and.associated(elj1c_p,elj1c_cu).and..not.use_group
     &         .and.n.eq.nvdw.and..not.use_lambdadyn
     &         .and.app_id.eq.dynamic_a.and.btest(tinEssai,0)) then
               fuse_chglj   = .true.
               vlst2_enable = .false.
               elj1c_p      => tinker_void_sub
               ecreal1d_p   => ecreal_lj1_cu
               call kcharge(.true.,0)
            end if
#endif
         else
            ecreal1d_p          => ecreal1dgpu
            ecreal3d_p          => ecreal3dgpu
            clst_enable = .true.
         end if
      end if

      if (use_mpole) then
         shft = ishft(sub_config,-conf_mpole*sh_len)
         configure=iand(shft,hexa_len-1)
         emreal1c_p     => emreal1c_
         emreal1ca_p    => emreal1ca_ac
         emreal3d_p     => emreal3d_ac
         if      (configure.eq.conf_emreal1c_0) then
            mlst_enable = .true.
         else if (configure.eq.conf_emreal1c_1) then
            mlst_enable = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_emreal1c_2) then
            emreal1ca_p => emreal1ca_cu
            emreal3d_p => emreal3d_cu
            mlst2_enable = .true.
         else if (configure.eq.conf_emreal1c_3) then
            emreal1ca_p => emreal1c_core4
            emreal3d_p => emreal3d_cu
            mlst2_enable = .true.
#endif
         else
            mlst_enable = .true.
         end if
      end if

      if (use_polar) then
         shft = ishft(sub_config,-conf_efld0*sh_len)
         configure=iand(shft,hexa_len-1)
         otf_dc_efld0_directgpu_p => otf_dc_efld0_directgpu2
         if      (configure.eq.conf_efld0_directgpu1) then
            efld0_directgpu_p => efld0_directgpu2
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
         shft = merge(3*n/(2*nproc),n,nproc.gt.1)
         puscale = u1scale*u2scale*u3scale*u4scale
         if      (configure.eq.conf_tmatxb_ref) then
            nullify(tmatxb_p)
            tmatxb_p          => tmatxb_pmegpu1
            mlst_enable = .true.
         else if (polar_precomp.and.use_thole.and..not.use_lambdadyn
     &           .and.puscale.eq.1.0.and.
     &           (configure.eq.conf_tmatxb_pre.or.
     &             (sub_config.eq.itrf_legacy
#ifdef _OPENACC
     &              .and.allow_precomp(shft)
#endif
     &           ))) then
            tmatxb_pme_core_p => tmatxb_pme_compute
            mlst_enable = .true.
   33       format(/,' --- Tinker-HP: Enabling precomputation in'
     &            ,'polarisation solver')
            if (tinkerdebug.gt.0) print 33
         else if (configure.eq.conf_tmatxb_c1) then
            tmatxb_pme_core_p => tmatxb_pme_core1
            mlst_enable = .true.
         else if (configure.eq.conf_tmatxb_c3) then
            tmatxb_pme_core_p => tmatxb_pme_core2
            mlst_enable = .true.
#ifdef _OPENACC
         else if (configure.eq.conf_tmatxb_c2) then
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
         if (.not.associated(tmatxb_pme_core_p,tmatxb_pme_compute))
     &      polar_precomp = .false.

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
            use_mpolar_ker  =
     &      merge(.true.,.false.,use_mpole.and.tinkerdebug.le.0
     &           .and..not.use_lambdadyn.and..not.use_chgpen
     &           .and..not.use_dirdamp.and..not.use_group
     &           .and..not.use_chgflx)
            !TODO Add group
            call attach_mpolar_pointer
#endif
         else
            epreal1c_core_p => epreal1c_core2
            epreal3d_p      => epreal3dgpu
            mlst_enable     = .true.
         end if
      end if

      if (use_chgtrn) then
#ifdef _OPENACC
         mlst2_enable = .true.
#else
         mlst_enable  = .true.
#endif
      end if

      if (use_disp) then
#ifdef _OPENACC
         vlst2_enable = .true.
#else
         vlst_enable = .false.
#endif
      end if

      shft = ishft(sub_config,-conf_fphi*sh_len)
      configure=iand(shft,hexa_len-1)
      if      (configure.eq.conf_fphi_acc) then
         fphi_uind_site1_p => fphi_uind_sitegpu1
         fphi_uind_site2_p => fphi_uind_sitegpu2
         fphi_mpole_site_p => fphi_mpole_sitegpu
         grid_pchg_force_p => grid_pchg_force
         grid_disp_force_p => grid_disp_force
#ifdef _OPENACC
      else if (configure.eq.conf_fphi_cuda) then
         fphi_uind_site1_p => fphi_uind_sitecu1
         fphi_uind_site2_p => fphi_uind_sitecu2
         fphi_mpole_site_p => fphi_mpole_sitecu
         grid_pchg_force_p => grid_pchg_forcecu
         grid_disp_force_p => grid_disp_forcecu
#endif
      else
         fphi_uind_site1_p => fphi_uind_sitegpu1
         fphi_uind_site2_p => fphi_uind_sitegpu2
         fphi_mpole_site_p => fphi_mpole_sitegpu
         grid_pchg_force_p => grid_pchg_force
         grid_disp_force_p => grid_disp_force
      end if

      shft = ishft(sub_config,-conf_grid*sh_len)
      configure=iand(shft,hexa_len-1)
      if      (configure.eq.conf_grid_acc) then
         grid_uind_site_p  => grid_uind_sitegpu
         grid_mpole_site_p => grid_mpole_sitegpu
         grid_pchg_site_p  => grid_pchg_sitegpu
         grid_disp_site_p  => grid_disp_sitegpu
#ifdef _OPENACC
      else if (configure.eq.conf_grid_cuda) then
         if (nproc.eq.1) then
            grid_uind_site_p  => grid_uind_sitecu
            grid_mpole_site_p => grid_mpole_sitecu
            grid_pchg_site_p  => grid_pchg_sitecu
            grid_disp_site_p  => grid_disp_sitecu
         else if (use_pmecore.and.nrec.eq.1) then
            grid_uind_site_p  => grid_uind_sitecu
            grid_mpole_site_p => grid_mpole_sitecu
            grid_pchg_site_p  => grid_pchg_sitecu
            grid_disp_site_p  => grid_disp_sitecu
         else
            grid_uind_site_p  => grid_uind_sitegpu
            grid_mpole_site_p => grid_mpole_sitegpu
            grid_pchg_site_p  => grid_pchg_sitegpu
            grid_disp_site_p  => grid_disp_sitegpu
         end if
#endif
      else
         grid_uind_site_p  => grid_uind_sitegpu
         grid_mpole_site_p => grid_mpole_sitegpu
         grid_pchg_site_p  => grid_pchg_sitegpu
         grid_disp_site_p  => grid_disp_sitegpu
      end if

      shft     = ishft(sub_config,-conf_conv*sh_len)
      configure= iand (shft,hexa_len-1)
      if      (configure.eq.conf_conv_acc) then
         pme_conv_p => pme_conv_gpu
#ifdef _OPENACC
      else if (configure.eq.conf_conv_cuda) then
         pme_conv_p => pme_conv_cu
#endif
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
      use polpot
      use potent ,only: use_pmecore
      implicit none
      logical,parameter:: OK=.false.

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

         ecreal1d_cp     => ecreal1d_p
         ecreal1d_p      => tinker_void_sub

         emreal1c_cp     => emreal1c_p
         emreal1c_p      => tinker_void_sub

         efld0_direct_cp   => efld0_directgpu_p
         efld0_directgpu_p => efld0_direct_void

         if (polalg.eq.pcg_SId.or.polalg.eq.step_pcg_SId
     &   .or.polalg.eq.step_pcg_short_SId) then
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
         emreal1c_p      => emreal1c_
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
      use tinheader ,only: isrel_build
      implicit none
      logical,parameter::accept=.false.

      if (nproc.gt.7.and.accept.and.use_polar.and.
     &   .not.isrel_build.and.
     &   (associated(efld0_directgpu_p,efld0_directgpu3).or.
     &    associated(efld0_direct_cp,efld0_directgpu3)).and.
     &    associated(tmatxb_pme_core_p,tmatxb_pme_core3))
     &   then

         ! Enable ExtraComputation
         no_commdir=.true.
 12      format(/,
     &   ' *****  Bypassing most of real space communications *****')
         if (rank.eq.0.and.verbose) write(*,12)

         ! Take into account neighbor extension
         call ddpme3d
         call reassignpme(.true.)

         ! Reset poleglobnl for each process
         if (use_mpole)  call kmpole (.false.,0)
         if (use_polar)  call kpolar (.false.,0)
      end if
#ifdef _OPENACC
      call cu_update_balanced_comput(.not.no_commdir)
#endif
      end subroutine
