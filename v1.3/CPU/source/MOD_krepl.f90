!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module krepl  --  Pauli repulsion forcefield parameters  ##
!     ##                                                           ##
!     ###############################################################
!
!
!
module krepl
   use sizes
   implicit none
   real*8 :: prsiz(maxtyp) !<Pauli repulsion size value for each atom class
   real*8 :: prdmp(maxtyp) !<alpha Pauli repulsion parameter for each atom class
   real*8 :: prele(maxtyp) !<number of valence electrons for each atom class
   save
end
