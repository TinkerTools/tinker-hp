!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module katoms    --  forcefield parameters for the atom types  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     weight     !<average atomic mass of each atom type
!     atmcls     !<atom class number for each of the atom types
!     atmnum     !<atomic number for each of the atom types
!     ligand     !<number of atoms to be attached to each atom type
!     symbol     !<modified atomic symbol for each atom type
!     describe   !<string identifying each of the atom types
!
!
module katoms
   use sizes
   implicit none
   integer :: atmcls(maxtyp) !<atom class number for each of the atom types
   integer :: atmnum(maxtyp) !<atomic number for each of the atom types
   integer :: ligand(maxtyp) !<number of atoms to be attached to each atom type
   real*8 :: weight(maxtyp) !<average atomic mass of each atom type
   character*3 :: symbol(maxtyp) !<modified atomic symbol for each atom type
   character*24 :: describe(maxtyp) !<string identifying each of the atom types
   save
end
