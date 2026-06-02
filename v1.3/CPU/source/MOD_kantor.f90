!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  module kantor  --  angle-torsion forcefield parameters  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     maxnat   !<maximum number of angle-torsion parameter entries
!
!     atcon    !<torsional amplitude parameters for angle-torsion
!     kat      !<string of atom classes for angle-torsion terms
!
!
module kantor
   implicit none
   integer :: maxnat !<maximum number of angle-torsion parameter entries
   parameter (maxnat=500)
   real*8 :: atcon(6,maxnat) !<torsional amplitude parameters for angle-torsion
   character*16 :: kat(maxnat) !<string of atom classes for angle-torsion terms
   save
end
