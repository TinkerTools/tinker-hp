!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module urey  --  Urey-Bradley interactions in the structure  ##
!     ##                                                               ##
!     ###################################################################
!
!
!
module urey
   implicit none
   integer :: nurey !<total number of Urey-Bradley terms in the system
   integer :: nureyloc !<local number of Urey-Bradley terms in the system
   integer, pointer :: iury(:,:) !<numbers of the atoms in each Urey-Bradley interaction
   integer, pointer :: nburey(:) !<numbers of Urey-Bradley interactions before each atom
   real*8, pointer ::  uk(:) !<Urey-Bradley force constants (kcal/mole/Ang**2)
   real*8, pointer :: ul(:) !<ideal 1-3 distance values in Angstroms
   integer :: winiury !<window object corresponding to iury
   integer :: winnburey !<window object corresponding to nburey
   integer :: winuk !<window object corresponding to uk
   integer :: winul !<window object corresponding to ul
   integer :: winureytyp !<window object corresponding to ureytyp
   character*8, pointer :: ureytyp(:) !<type of urey-bradley interaction
   save
end
