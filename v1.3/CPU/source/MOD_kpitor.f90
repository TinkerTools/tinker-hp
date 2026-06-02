!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kpitor  --  forcefield parameters for pi-orbit torsions  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module kpitor
   implicit none
   integer maxnpt !<maximum number of pi-orbital torsion parameter entries
   parameter (maxnpt=500)
   real*8 ptcon(maxnpt) !<force constant parameters for pi-orbital torsions
   character*8 kpt(maxnpt) !<string of atom classes for pi-orbital torsion terms
   save
end
