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
!
module ewald
   implicit none
   real*8 :: aewald !<Ewald convergence coefficient value (Ang-1)
   real*8 :: adewald !<Ewald convergence coefficient value (Ang-1) for dispersion
   real*8 :: aeewald !<Ewald convergence coefficient value (Ang-1) for electrostatics
   real*8 :: apewald !<Ewald convergence coefficient value (Ang-1) for polarization
   character*7 :: boundary !<Ewald boundary condition; none, tinfoil or vacuum
   save
end
