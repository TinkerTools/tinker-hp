!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ## module  khbond  --  forcefield parameters for H-bonding terms  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     maxnhb   maximum number of hydrogen bonding pair entries
!
!     radhb    radius parameter for hydrogen bonding pairs
!     epshb    well depth parameter for hydrogen bonding pairs
!     khb      string of atom types for hydrogen bonding pairs
!
!
#include "tinker_macro.h"
module khbond
   implicit none
   integer maxnhb
   parameter (maxnhb=500)
   real(t_p) radhb(maxnhb),epshb(maxnhb)
   character*8 khb(maxnhb)
   save
end
