!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module mdstate --  control of molecular dynamics trajectory  ##
!     ##                                                               ##
!     ###################################################################
!
!
#include "tinker_macro.h"
!
!     track_mds : tracking MD state switch
!     ms_back_p : MD State backup period
!     md_state  : new type to hold an MD State
!     ms        : md_state type for other
!
!
module mdstate
   implicit none
   logical track_mds,fw_mds
   integer ms_back_p
   character(5) :: energy_format = 'F15.4'
   logical :: has_energy_fmt = .false.
   type md_state
      logical isothermal,isobaric
      integer istep,randseed,n
      real(r_p) kelvin,atmsph
      real(r_p) dt
      real(r_p) xbox,ybox,zbox
      character*64 dumpdyn,dumpdat
      logical  lddyn,lddat
      real(r_p),allocatable::state(:)
   end type
   type(md_state), target:: ms(2)
   real(r_p)     ,pointer:: ms_x(:),ms_y(:),ms_z(:)&
      &,ms_v(:),ms_a(:),ms_alt(:)

   interface
      module subroutine set_energy_format(total_energy)
         real(r_p), intent(in) :: total_energy
      end subroutine
      module subroutine write_energy_hdr(iout, include_temp_cl, include_pressure)
         integer, intent(in) :: iout
         logical, intent(in) :: include_temp_cl, include_pressure
      end subroutine
      module subroutine write_energy_row(iout, istep, total_energy, potential_energy, &
         &kinetic_energy, temperature, temp_cent, pressure, include_temp_cl, &
         &include_pressure)
         integer, intent(in) :: iout, istep
         real(r_p), intent(in) :: total_energy, potential_energy, kinetic_energy
         real(r_p), intent(in) :: temperature, temp_cent, pressure
         logical, intent(in) :: include_temp_cl, include_pressure
      end subroutine
      module subroutine mds_init
      end subroutine
      module subroutine mds_save
      end subroutine
      module subroutine mds_prt
      end subroutine
   end interface

end module
