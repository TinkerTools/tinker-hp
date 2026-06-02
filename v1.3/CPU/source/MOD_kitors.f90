!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kitors  --  forcefield parameters for improper torsions  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module kitors
   implicit none
   integer ::  maxnti !<maximum number of improper torsion parameter entries
   parameter (maxnti=500)
   real*8 :: ti1(2,maxnti) !<torsional parameters for improper 1-fold rotation
   real*8 :: ti2(2,maxnti) !<torsional parameters for improper 2-fold rotation
   real*8 :: ti3(2,maxnti) !<torsional parameters for improper 3-fold rotation
   character*16 :: kti(maxnti) !<string of atom classes for improper torsional parameters
   save
end
