!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ########################################################################
!     ##                                                                    ##
!     ##  module langevin  --  parameters and arrays for langevin dynamics  ##
!     ##                                                                    ##
!     ########################################################################
!
!
module langevin
   implicit none
   real*8 gamma !<friction parameter in ps-1
   real*8, allocatable :: gamma_friction(:) !<friction parameter in ps-1 (PIMD)
   real*8, allocatable :: Rn(:,:) !<white noise
   save
end
