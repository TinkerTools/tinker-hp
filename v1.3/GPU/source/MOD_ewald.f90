!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module ewald  --  parameters and options for Ewald summation  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     aewald     Ewald convergence coefficient value (Ang-1)
!     ewald_tiny tiniest Ewald convergence coefficient value (Ang-1)
!     aeewald    Ewald convergence coefficient value (Ang-1) for electrostatics
!     apewald    Ewald convergence coefficient value (Ang-1) for polarization
!     adwald     Ewald convergence coefficient value (Ang-1) for dispersion
!     boundary   Ewald boundary condition; none, tinfoil or vacuum
!
!
#include "tinker_macro.h"
module ewald
   implicit none
   real(t_p) aewald,adewald,aeewald,apewald,ewald_tiny
   character*7 boundary
   parameter( ewald_tiny=1d-6 )
end
