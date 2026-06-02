!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #######################################################################
!     ##                                                                   ##
!     ##  module kurybr  --  forcefield parameters for Urey-Bradley terms  ##
!     ##                                                                   ##
!     #######################################################################
!
!
!     maxnu   maximum number of Urey-Bradley parameter entries
!
!     ucon    force constant parameters for Urey-Bradley terms
!     dst13   ideal 1-3 distance parameters for Urey-Bradley terms
!     ku      integer signature of atom classes for Urey-Bradley terms
!     ku_sys  integer signature of atom classes for Urey-Bradley terms of the simulated system
!
!
#include "tinker_macro.h"
module kurybr
   implicit none
   integer maxnu,maxnups,maxnuq
   parameter (maxnu=2000)
   parameter (maxnups=500)
   parameter (maxnuq=500)
   real(t_p) uconps(maxnups), dst13ps(maxnups)
   real(t_p) uconq(maxnuq), dst13q(maxnuq)
   real(t_p) ucon(maxnu),dst13(maxnu)
   integer(8) ku(maxnu),kups(maxnups),kuq(maxnuq)
   integer(8) ku_sys(0:maxnu+maxnups+maxnuq)
   save
!$acc declare create(ku,kups,ku_sys,kuq)
end
