!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine field  --  get the potential energy functions  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "field" sets the force field potential energy functions from
!     a parameter file and modifications specified in a keyfile
!
!
!> @brief 
!> sets the force field potential energy functions from
!> a parameter file and modifications specified in a keyfile
!> @param no params
subroutine field
   use inform
   use iounit
   use keys
   use potent
   implicit none
   integer i
   character*240 record
!
!
   if (deb_Path) write(iout,*), 'field '
!
!
!     set the default values for the active potentials
!
   use_bond = .true.
   use_angle = .true.
   use_strbnd = .true.
   use_urey = .true.
   use_angang = .true.
   use_opbend = .true.
   use_opdist = .true.
   use_improp = .true.
   use_imptor = .true.
   use_tors = .true.
   use_pitors = .true.
   use_angtor = .true.
   use_strtor = .true.
   use_tortor = .true.
   use_vdw = .true.
   use_vdwshort = .false.
   use_vdwlong = .false.
   use_charge = .true.
   use_cshortreal = .false.
   use_clong = .false.
   use_creal = .true.
   use_crec = .true.
   use_cself = .true.
   use_mpole = .true.
   use_mpoleshortreal = .false.
   use_mpolelong = .false.
   use_mreal = .true.
   use_mrec = .true.
   use_mself = .true.
   use_polar = .true.
   use_preal = .true.
   use_prec = .true.
   use_pself = .true.
   use_solv = .false.
   use_geom = .true.
   use_extra = .true.
   use_repulsshort = .false.
   use_repulslong = .false.
   use_repuls = .true.
   use_disp = .true.
   use_dispshortreal = .false.
   use_displong = .false.
   use_dispreal = .true.
   use_disprec = .true.
   use_dispself = .true.
   use_chgtrn = .true.
   use_chgflx = .true.
   use_lambdadyn = .false.
!
!     read the potential energy force field parameter file
!
   call getprm
!
!     check keywords for potential function control parameters
!
   do i = 1, nkey
      record = keyline(i)
      call prmkey (record)
   end do
   return
end
