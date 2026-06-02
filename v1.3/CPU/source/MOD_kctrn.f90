!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module kctrn  --  charge transfer forcefield parameters  ##
!     ##                                                           ##
!     ###############################################################
!
!
module kctrn
   use sizes
   implicit none
   real*8 :: ctchg(maxtyp) !<charge transfer magnitude for each atom class
   real*8 :: ctdmp(maxtyp) !<alpha charge transfer parameter for each atom class
   save
end
