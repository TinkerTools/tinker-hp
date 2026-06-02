!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module kgeoms  --  parameters for the geometrical restraints  ##
!     ##                                                                ##
!     ####################################################################
!
!
!
module kgeoms
   implicit none
   integer :: npfix !<number of position restraints to be applied
   integer :: ndfix !<number of distance restraints to be applied
   integer :: nafix !<number of angle restraints to be applied
   integer :: ntfix !<number of torsion restraints to be applied
   integer :: ngfix !<number of group restraints to be applied
   integer :: nchir !<number of chirality restraints to be applied
   integer :: npfixloc !<number of local position restraints to be applied
   integer :: ndfixloc !<number of local distance restraints to be applied
   integer :: nafixloc !<number of local angle restraints to be applied
   integer :: ntfixloc !<number of local torsion restraints to be applied
   integer :: ngfixloc !<number of local group restraints to be applied
   integer :: nchirloc !<number of local chirality restraints to be applied
   integer :: winipfix !< window object corresponding to ipfix
   integer :: winkpfix !< window object corresponding to kpfix
   integer :: winidfix !< window object corresponding to idfix
   integer :: winiafix !< window object corresponding to iapfix
   integer :: winitfix !< window object corresponding to itfix
   integer :: winigfix !< window object corresponding to igfix
   integer :: winichir !< window object corresponding to ichir
   integer :: winxpfix !< window object corresponding to xpfix
   integer :: winypfix !< window object corresponding to ypfix
   integer :: winzpfix !< window object corresponding to zpfix
   integer :: winpfix !< window object corresponding to pfix
   integer :: windfix !< window object corresponding to dfix
   integer :: winafix !< window object corresponding to afix
   integer :: wintfix !< window object corresponding to tfix
   integer :: wingfix !< window object corresponding to gfix
   integer :: winchir !< window object corresponding to chir
   integer, pointer :: ipfix(:) !<atom number involved in each position restraint
   integer, pointer :: kpfix(:,:) !<flags to use x-, y-, z-coordinate position restraints
   integer, pointer :: idfix(:,:) !<atom numbers defining each distance restraint
   integer, pointer :: iafix(:,:) !<atom numbers defining each angle restraint
   integer, pointer :: itfix(:,:) !<atom numbers defining each torsional restraint
   integer, pointer :: igfix(:,:) !<group numbers defining each group distance restraint
   integer, pointer :: ichir(:,:) !<atom numbers defining each chirality restraint
   logical :: use_basin !<logical flag governing use of Gaussian basin
   logical :: use_wall !<logical flag governing use of droplet boundary
   real*8, pointer :: xpfix(:) !<x-coordinate target for each restrained position
   real*8, pointer :: ypfix(:) !<y-coordinate target for each restrained position
   real*8, pointer :: zpfix(:) !<z-coordinate target for each restrained position
   real*8, pointer :: pfix(:,:) !<force constant and flat-well range for each position
   real*8, pointer :: dfix(:,:) !<force constant and target range for each distance
   real*8, pointer :: afix(:,:) !<force constant and target range for each angle
   real*8, pointer :: tfix(:,:) !<force constant and target range for each torsion
   real*8, pointer :: gfix(:,:) !<force constant and target range for each group distance
   real*8, pointer :: chir(:,:) !<force constant and target range for chiral centers
   real*8 :: depth !<depth of shallow Gaussian basin restraint
   real*8 :: width !<exponential width coefficient of Gaussian basin
   real*8 :: rwall !<radius of spherical droplet boundary restraint
   save
end
