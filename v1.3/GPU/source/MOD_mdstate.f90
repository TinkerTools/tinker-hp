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
      module subroutine mds_init
      end subroutine
      module subroutine mds_save
      end subroutine
      module subroutine mds_prt
      end subroutine
   end interface

end module
