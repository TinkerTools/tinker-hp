!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kiprop  --  forcefield parameters for improper dihedral  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!
module kiprop
   implicit none
   integer :: maxndi !<maximum number of improper dihedral parameter entries
   parameter (maxndi=500)
   real*8 :: dcon(maxndi) !<force constant parameters for improper dihedrals
   real*8 :: tdi(maxndi) !<ideal dihedral angle values for improper dihedrals
   character*16 :: kdi(maxndi) !<string of atom classes for improper dihedral angles
   save
end
