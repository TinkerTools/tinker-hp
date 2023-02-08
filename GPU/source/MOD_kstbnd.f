c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module kstbnd  --  forcefield parameters for stretch-bend  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnsb   maximum number of stretch-bend parameter entries
c
c     stbn     force constant parameters for stretch-bend terms
c     ksb      integer signature of atom classes for stretch-bend terms
c     ksb_sys  integer signature of atom classes for stretch-bend terms of the system simulated
c
c
#include "tinker_macro.h"
      module kstbnd
      implicit none
      integer maxnsb
      parameter (maxnsb=2000)
      integer(8) ksb(maxnsb)
      integer(8) ksb_sys(0:maxnsb)
      real(t_p) stbn(2,maxnsb)
      save
!$acc declare create(ksb,ksb_sys)
      end
