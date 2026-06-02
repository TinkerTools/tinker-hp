!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module kcflux -- charge flux term forcefield parameters  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     maxncfb   maximum number of bond stretch charge flux entries
!     maxncfa   maximum number of angle bend charge flux entries
!
!     cflb      charge flux over stretching of a bond length
!     cfla      charge flux over bending of a bond angle
!     cflab     charge flux over asymmetric bond within an angle
!     kcfb      string of atom classes for bond stretch charge flux
!     kcfa      string of atom classes for angle bend charge flux
!
!
#include "tinker_macro.h"
module kcflux
   implicit none
   integer maxncfb
   integer maxncfa
   parameter (maxncfb=2000)
   parameter (maxncfa=2000)
   real(t_p) cflb(maxncfb)
   real(t_p) cfla(2,maxncfa)
   real(t_p) cflab(2,maxncfa)
   integer(8) kcfb(maxncfb), kcfa(maxncfa)
end
