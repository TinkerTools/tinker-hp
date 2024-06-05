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
      subroutine mechanic
      use domdec
      use inform
      use iounit
      use potent
      use mpi
      implicit none
c
      if (deb_Path) write(iout,*), 'mechanic '
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c     assign atom types, classes and other atomic information
c
      call katom
c
c     assign atoms to molecules and set the atom groups
c
      call molecule(.true.)
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
c      call orbital
c
c     assign electrostatic and dispersion Ewald sum parameters
c
      call kewald
c
c     assign bond, angle and cross term potential parameters
c
      call kbond
      call kangle
      call kstrbnd
      call kurey
      call kangang
c
c     assign out-of-plane deformation potential parameters
c
      call kopbend
      call kopdist
      call kimprop
      call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      call ktors
      call kpitors
      call kstrtor
      call kangtor
      call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      call kcharge
      call kvdw
      call kmpole
      call kpolar

      call kchgtrn
      call kchgflx
c
c      if (use_polar) call initmpipme
c
c     assign repulsion and dispersion parameters
c
      call krepel 
      call kdisp 
c
c     assign restraint parameters
c
      call kgeom
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     set holonomic constrains
c
      call shakeup
c
c     SMD parametrization
c
      call ksmd(.true.)
c
c     quit if essential parameter information is missing
c
      if (abort) then
         if (rank.eq.0) write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
c
c     subroutine mechanic_init_para: initialize parallel parameters after domain decomposition
c
c
      subroutine mechanic_init_para
      use inform
      use iounit
      use potent
      implicit none
c
      if (deb_Path) write(iout,*), 'mechanic_init_para '
c
c
      call bonds_update
      call angles_update
      call torsions_update
      call bitors_update

      call kewald_dd_init
      
      if (use_strbnd) call kstrbnd_update
      if (use_urey)   call kurey_update
      if (use_angang)  call kangang_update

      if (use_opbend)  call kopbend_update
      if (use_opdist) call kopdist_update
      if (use_improp) call kimprop_update
      if (use_imptor) call kimptor_update

      if (use_pitors) call kpitors_update
      if (use_strtor) call kstrtor_update
      if (use_angtor) call kangtor_update
      if (use_tortor) call ktortor_update

      if (use_charge) call kcharge_update(0)
      if (use_vdw)    call kvdw_update(0)
      if (use_mpole)  call kmpole_update(0)
      if (use_polar)  call kpolar_update(0)

      if (use_chgtrn) call kchgtrn_update(0)

      if (use_polar)  call initmpipme

      if (use_disp)   call kdisp_update(0)

      if (use_geom)   call kgeom_update

      call shakeup_update

      if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
c
      return
      end
c
c     subroutine mechanic_up_para: update parallel parameters
c
c
      subroutine mechanic_up_para(istep)
      use inform
      use iounit
      use potent
      implicit none
      integer istep
c
      if (deb_Path) write(iout,*), 'mechanic_up_para '
c
c
      call bonds_update
      call angles_update
      call torsions_update
      call bitors_update

c      call molecule_update
      
      if (use_strbnd) call kstrbnd_update
      if (use_urey)   call kurey_update
      if (use_angang)  call kangang_update

      if (use_opbend)  call kopbend_update
      if (use_opdist) call kopdist_update
      if (use_improp) call kimprop_update
      if (use_imptor) call kimptor_update

      if (use_pitors) call kpitors_update
      if (use_strtor) call kstrtor_update
      if (use_angtor) call kangtor_update
      if (use_tortor) call ktortor_update

      if (use_charge) call kcharge_update(istep)
      if (use_vdw)    call kvdw_update(istep)
      if (use_mpole)  call kmpole_update(istep)
      if (use_polar)  call kpolar_update(istep)

      if (use_chgtrn) call kchgtrn_update(istep)

      if (use_polar)  call initmpipme

      if (use_disp)   call kdisp_update(istep)

      if (use_geom)   call kgeom_update

      call shakeup_update

      if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
c
      return
      end

c
c     subroutine mechanic_up_para_respa: update parallel parameters between two respa steps
c
c
      subroutine mechanic_up_para_respa(istep,fast)
      use inform
      use iounit
      use potent
      implicit none
      logical fast
      integer istep
c
      if (deb_Path) write(iout,*), 'mechanic_up_para_respa '
c
c
      if (fast) then
        call bonds_update
        call angles_update
        call torsions_update
        call bitors_update
        if (use_strbnd) call kstrbnd_update
        if (use_urey) call kurey_update
        if (use_angang) call kangang_update
        if (use_opbend)  call kopbend_update
        if (use_opdist)  call kopdist_update
        if (use_improp)  call kimprop_update
        if (use_imptor)  call kimptor_update
        if (use_pitors)  call kpitors_update
        if (use_strtor)  call kstrtor_update
        if (use_angtor)  call kangtor_update
        if (use_tortor)  call ktortor_update
        if (use_geom)  call kgeom_update
        if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
      else
        if (use_charge) call kcharge_update(istep)
        if (use_mpole) call kmpole_update(istep)
        if (use_polar) call kpolar_update(istep)
        if (use_vdw) call kvdw_update(istep)
        if (use_disp) call kdisp_update(istep)
        if (use_chgtrn) call kchgtrn_update(istep)
        if (use_polar) call initmpipme
c
c     set holonomic constrains
c
        call shakeup_update
      end if
      return
      end
c
c     subroutine mechanic_up_para_respa1: update parameters between two time steps
c
      subroutine mechanic_up_para_respa1(istep,rule)
      use domdec
      use inform
      use iounit
      use potent
      implicit none
      integer istep,rule
 1000 format(' illegal rule in mechanicsteprespa1.')
c
      if (deb_Path) write(iout,*), 'mechanic_up_para_respa1 '
c
c
c     rule = 0: fast part of the forces
c     rule = 1: intermediate part of the forces
c     rule = 2: slow part of the forces
c
c      call molecule(.false.)
      if (rule.eq.0) then
        call bonds_update
        call angles_update
        call torsions_update
        call bitors_update
        if (use_strbnd) call kstrbnd_update
        if (use_urey) call kurey_update
        if (use_angang) call kangang_update
        if (use_opbend)  call kopbend_update
        if (use_opdist)  call kopdist_update
        if (use_improp)  call kimprop_update
        if (use_imptor)  call kimptor_update
        if (use_pitors)  call kpitors_update
        if (use_strtor)  call kstrtor_update
        if (use_angtor)  call kangtor_update
        if (use_tortor)  call ktortor_update
        if (use_geom)  call kgeom_update
        if (use_smd_velconst .or. use_smd_forconst) call ksmd(.false.)
      else if (rule.eq.1) then
        if (use_charge) call kcharge_update(istep)
        if (use_mpole) call kmpole_update(istep)
        if (use_polar) call kpolar_update(istep)
        if (use_vdw) call kvdw_update(istep)
        if (use_chgtrn) call kchgtrn_update(istep)
      else if (rule.eq.2) then
        if (use_charge) call kcharge_update(istep)
        if (use_mpole) call kmpole_update(istep)
        if (use_polar) call kpolar_update(istep)
        if (use_vdw) call kvdw_update(istep)
        if (use_disp) call kdisp_update(istep)
        if (use_chgtrn) call kchgtrn_update(istep)
        if ((istep.ne.-1).and.use_polar) call initmpipme
c
c     set holonomic constrains
c
        call shakeup_update
      else
         if (rank.eq.0) write(iout,1000) 
      end if
      return
      end

c
