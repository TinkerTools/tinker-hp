!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module ksttor  --  forcefield parameters for stretch-torsions  ##
!     ##                                                                 ##
!     #####################################################################
!
!
module ksttor
   implicit none
   integer :: maxnbt !<maximum number of stretch-torsion parameter entries
   parameter (maxnbt=500)
   real*8 :: btcon(9,maxnbt) !<force constant parameters for stretch-torsion
   character*16 :: kbt(maxnbt) !<string of atom classes for stretch-torsion terms
   save
end
