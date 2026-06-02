!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #######################################################################
!     ##                                                                   ##
!     ##  module pitors  --  pi-orbital torsions in the current structure  ##
!     ##                                                                   ##
!     #######################################################################
!
!
module pitors
   implicit none
   integer :: npitors !<total number of pi-orbital torsional interactions
   integer :: npitorsloc !<local number of pi-orbital torsional interactions
   integer :: winipit !<window object corresponding to ipit
   integer :: winnbpitors !<window object corresponding to nbpitors
   integer :: winkpit !<window object corresponding to kpit
   integer, pointer :: ipit(:,:) !<numbers of the atoms in each pi-orbital torsion
   integer, pointer :: nbpitors(:) !<number of pi-orbital torsional interactions before each atom
   real*8, pointer :: kpit(:) !<2-fold pi-orbital torsional force constants
   save
end
