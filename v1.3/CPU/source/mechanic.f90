!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine mechanic  --  initialize molecular mechanics  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "mechanic" sets up needed parameters for the potential energy
!     calculation and reads in many of the user selectable options
!
!
!> @brief 
!> sets up needed parameters for the potential energy
!> calculation and reads in many of the user selectable options
!> @param no params
subroutine mechanic
   use domdec
   use inform
   use iounit
   use potent
   use mpi
   implicit none
!
   if (deb_Path) write(iout,*), 'mechanic '
!
!
!     set the bonded connectivity lists and active atoms
!
   call attach
   call active
!
!     find bonds, angles, torsions, bitorsions and small rings
!
   call bonds
   call angles
   call torsions
   call bitors
   call rings
!
!     get the base force field from parameter file and keyfile
!
   call field
!
!     assign atom types, classes and other atomic information
!
   call katom
!
!     assign atoms to molecules and set the atom groups
!
   call molecule(.true.)
   call cluster
!
!     find any pisystem atoms, bonds and torsional angles
!
!      call orbital
!
!     assign electrostatic and dispersion Ewald sum parameters
!
   call kewald
!
!     assign bond, angle and cross term potential parameters
!
   call kbond
   call kangle
   call kstrbnd
   call kurey
   call kangang
!
!     assign out-of-plane deformation potential parameters
!
   call kopbend
   call kopdist
   call kimprop
   call kimptor
!
!     assign torsion and torsion cross term potential parameters
!
   call ktors
   call kpitors
   call kstrtor
   call kangtor
   call ktortor
!
!     assign van der Waals and electrostatic potential parameters
!
   call kcharge
   call kvdw
   call kmpole
   call kpolar

   call kchgtrn
   call kchgflx
!
!      if (use_polar) call initmpipme
!
!     assign repulsion and dispersion parameters
!
   call krepel
   call kdisp
!
!     assign restraint parameters
!
   call kgeom
!
!     set hybrid parameter values for free energy perturbation
!
   call mutate
!
!     set holonomic constrains
!
   call shakeup
!
!     SMD parametrization
!
   call ksmd(.true.)
!
!     quit if essential parameter information is missing
!
   if (abort) then
      if (rank.eq.0) write (iout,10)
10    format (/,' MECHANIC  --  Some Required Potential Energy',&
      &' Parameters are Undefined')
      call fatal
   end if
   return
end
!
!     subroutine mechanic_init_para: initialize parallel parameters after domain decomposition
!
!
!> @brief 
!> initialize parallel parameters after domain decomposition
!> @param no params
subroutine mechanic_init_para
   use inform
   use iounit
   use potent
   implicit none
!
   if (deb_Path) write(iout,*), 'mechanic_init_para '
!
!
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
!
   return
end
!
!     subroutine mechanic_up_para: update parallel parameters
!
!
!> @brief 
!> update parallel parameters
!> @param[in]  istep: index of timestep
subroutine mechanic_up_para(istep)
   use inform
   use iounit
   use potent
   implicit none
   integer istep
!
   if (deb_Path) write(iout,*), 'mechanic_up_para '
!
!
   call bonds_update
   call angles_update
   call torsions_update
   call bitors_update

!      call molecule_update

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
!
   return
end

!
!     subroutine mechanic_up_para_respa: update parallel parameters between two respa steps
!
!
!> @brief 
!> update parallel parameters between two respa steps
!> @param no params
subroutine mechanic_up_para_respa(istep,fast)
   use inform
   use iounit
   use potent
   implicit none
   logical fast
   integer istep
!
   if (deb_Path) write(iout,*), 'mechanic_up_para_respa '
!
!
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
!
!     set holonomic constrains
!
      call shakeup_update
   end if
   return
end
!
!     subroutine mechanic_up_para_respa1: update parameters between two time steps
!
!> @brief 
!> update parameters between two time steps of respa1 integrator
!> @param[in]  istep: index of timestep
!> @param[in]  rule: rule depending on short, intermediate or large timestep
subroutine mechanic_up_para_respa1(istep,rule)
   use domdec
   use inform
   use iounit
   use potent
   implicit none
   integer istep,rule
1000 format(' illegal rule in mechanicsteprespa1.')
!
   if (deb_Path) write(iout,*), 'mechanic_up_para_respa1 '
!
!
!     rule = 0: fast part of the forces
!     rule = 1: intermediate part of the forces
!     rule = 2: slow part of the forces
!
!      call molecule(.false.)
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
!
!     set holonomic constrains
!
      call shakeup_update
   else
      if (rank.eq.0) write(iout,1000)
   end if
   return
end

!
