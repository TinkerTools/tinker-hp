!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kvdwpr  --  forcefield parameters for special vdw terms  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module kvdwpr
   implicit none
   integer :: maxnvp !<maximum number of special van der Waals pair entries
   parameter (maxnvp=500)
   real*8 :: radpr(maxnvp) !<radius parameter for special van der Waals pairs
   real*8 :: epspr(maxnvp) !<well depth parameter for special van der Waals pairs
   character*8 :: kvpr(maxnvp) !<string of atom classes for special van der Waals pairs
   save
end
