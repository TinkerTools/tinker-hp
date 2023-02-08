c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module mdstate --  control of molecular dynamics trajectory  ##
c     ##                                                               ##
c     ###################################################################
c
c
#include "tinker_macro.h"
c
c     track_mds : tracking MD state switch
c     ms_back_p : MD State backup period
c     md_state  : new type to hold an MD State
c     ms        : md_state type for other
c
c
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
      real(r_p)     ,pointer:: ms_x(:),ms_y(:),ms_z(:)
     &              ,ms_v(:),ms_a(:),ms_alt(:)

      interface
       module subroutine mds_init
       end subroutine
       module subroutine mds_save
       end subroutine
       module subroutine mds_prt
       end subroutine
      end interface

      end module
