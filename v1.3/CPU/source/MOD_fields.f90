!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module fields  --  molecular mechanics force field description  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     biotyp       !<force field atom type of each biopolymer type
!     forcefield   !<string used to describe the current forcefield
!
!
module fields
   use sizes
   implicit none
   integer :: biotyp(maxbio) !<force field atom type of each biopolymer type
   character*20 :: forcefield !<string used to describe the current forcefield
   save
end
