!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kopbnd  --  forcefield parameters for out-of-plane bend  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module kopbnd
   implicit none
   integer :: maxnopb !<maximum number of out-of-plane bending entries
   parameter (maxnopb=500)
   logical, allocatable :: jopb(:)
   real*8 :: opbn(maxnopb) !<force constant parameters for out-of-plane bending
   character*16 :: kopb(maxnopb) !<string of atom classes for out-of-plane bending
end
