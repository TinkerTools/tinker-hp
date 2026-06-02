!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module atmlst    --  local geometry terms involving each atom  ##
!     ##                                                                 ##
!     #####################################################################
!
!
module atmlst
   implicit none
   integer :: winbndlist !<window object corresponding to bndlist
   integer :: winanglist !<window object corresponding to anglist
   integer :: winbalist !<window object corresponding to balist
   integer, pointer :: bndlist(:,:) !<list of the bond numbers involving each atom
   integer, pointer :: anglist(:,:) !<list of the angle numbers centered on each atom
   integer, pointer :: balist(:,:) !<numbers of the bonds comprising each angle
   integer, allocatable :: bndglob(:) !<local - global bond  correspondance
   integer, allocatable :: angleglob(:) !<local - global angle correspondance
   integer, allocatable :: torsglob(:) !<local - global torsion correspondance
   integer, allocatable :: bitorsglob(:) !<local - global bitorsion correspondance      
   integer, allocatable :: strbndglob(:) !<local - global strech bending correspondance 
   integer, allocatable :: ureyglob(:) !<local - global urey bradley correspondance
   integer, allocatable :: angangglob(:)  !<local - global angle correspondance
   integer, allocatable :: opbendglob(:) !<local - global out of plane bending correspondance 
   integer, allocatable :: opdistglob(:)  !<local - global out of plane distance correspondance
   integer, allocatable :: impropglob(:)  !<local - global improper dihedral correspondance
   integer, allocatable :: imptorglob(:)  !<local - global improper torsion correspondance
   integer, allocatable :: pitorsglob(:)  !<local - global pi torsion correspondance
   integer, allocatable :: strtorglob(:)  !<local - global strech torsion correspondance
   integer, allocatable :: angtorglob(:)  !<local - global angle-torsion correspondance
   integer, allocatable :: tortorglob(:)  !<local - global torsion torsion correspondance
   integer, allocatable :: vdwglob(:) !<local - global vdw correspondance
   integer, allocatable :: poleglob(:)  !<local - global direct multipole correspondance
   integer, allocatable :: polerecglob(:)  !<local - global reciprocal multipole correspondance
   integer, allocatable :: dispglob(:) !<local - global dispersion correspondance
   integer, allocatable :: disprecglob(:)  !<local - global reciprocal dispersion correspondance
   integer, allocatable :: chgglob(:)  !<local - global direct charge correspondance
   integer, allocatable :: chgrecglob(:)  !<local - global reciprocal charge correspondance
   integer, allocatable :: molculeglob(:)  !<local - global molecule correspondance
   integer, allocatable :: npfixglob(:)  !<local - global position restrains correspondance
   integer, allocatable :: ndfixglob(:)  !<local - global distance restrains correspondance
   integer, allocatable :: nafixglob(:)  !<local - global angle restrains correspondance
   integer, allocatable :: ntfixglob(:)  !<local - global torsion restrains correspondance
   integer, allocatable :: ngfixglob(:)  !<local - global group restrains correspondance
   integer, allocatable :: nchirglob(:)  !<local - global chiral restrains correspondance
   integer, allocatable :: ratglob(:)  !<local - <global constrains correspondance
   integer, allocatable :: chgglobnl(:)  !<localnl - global direct charge correspondance
   integer, allocatable :: vdwglobnl(:)  !<localnl - global vdw correspondance
   integer, allocatable :: poleglobnl(:)  !<localnl - global direct multipole correspondance
   integer, allocatable :: dispglobnl(:)  !<localnl - global dispersion correspondance
   save
end
