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
!     weight     average atomic mass of each atom type
!     atmcls     atom class number for each of the atom types
!     atmnum     atomic number for each of the atom types
!     ligand     number of atoms to be attached to each atom type
!     symbol     modified atomic symbol for each atom type
!     describe   string identifying each of the atom types
!
!
#include "tinker_macro.h"
module katoms
   use sizes ,only: maxtyp
   implicit none
   integer atmcls(maxtyp),atmnum(maxtyp)
   integer ligand(maxtyp)
   real(r_p) weight(maxtyp)
   character*3 symbol(maxtyp)
   character*24 describe(maxtyp)
   save
end
