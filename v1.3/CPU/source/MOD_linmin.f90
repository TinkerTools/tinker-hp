!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module linmin  --  parameters for line search minimization  ##
!     ##                                                              ##
!     ##################################################################
!
!
module  linmin
   implicit none
   integer :: intmax !<maximum number of interpolations during line search
   real*8 :: stpmin !<minimum step length in current line search direction
   real*8 :: stpmax !<maximum step length in current line search direction
   real*8 :: cappa !<stringency of line search (0=tight < cappa < 1=loose)
   real*8 :: slpmax !<projected gradient above which stepsize is reduced
   real*8 :: angmax !<maximum angle between search direction and -gradient
   save
end
