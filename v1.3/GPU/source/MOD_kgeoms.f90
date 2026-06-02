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
!     xpfix      x-coordinate target for each restrained position
!     winxpfix    window object corresponding to xpfix
!     ypfix      y-coordinate target for each restrained position
!     winypfix    window object corresponding to ypfix
!     zpfix      z-coordinate target for each restrained position
!     winzpfix    window object corresponding to zpfix
!     pfix       force constant and flat-well range for each position
!     winpfix    window object corresponding to pfix
!     dfix       force constant and target range for each distance
!     windfix    window object corresponding to dfix
!     afix       force constant and target range for each angle
!     winafix    window object corresponding to afix
!     tfix       force constant and target range for each torsion
!     wintfix    window object corresponding to tfix
!     gfix       force constant and target range for each group distance
!     wingfix    window object corresponding to gfix
!     chir       force constant and target range for chiral centers
!     winchir    window object corresponding to chir
!     depth      depth of shallow Gaussian basin restraint
!     width      exponential width coefficient of Gaussian basin
!     rwall      radius of spherical droplet boundary restraint
!     npfix      number of position restraints to be applied
!     ipfix      atom number involved in each position restraint
!     winipfix    window object corresponding to ipfix
!     kpfix      flags to use x-, y-, z-coordinate position restraints
!     winkpfix    window object corresponding to kpfix
!     ndfix      number of distance restraints to be applied
!     idfix      atom numbers defining each distance restraint
!     winidfix    window object corresponding to idfix
!     nafix      number of angle restraints to be applied
!     iafix      atom numbers defining each angle restraint
!     winiafix    window object corresponding to iafix
!     ntfix      number of torsional restraints to be applied
!     itfix      atom numbers defining each torsional restraint
!     winitfix    window object corresponding to itfix
!     ngfix      number of group distance restraints to be applied
!     igfix      group numbers defining each group distance restraint
!     winigfix    window object corresponding to igfix
!     nchir      number of chirality restraints to be applied
!     ichir      atom numbers defining each chirality restraint
!     winichir    window object corresponding to ichir
!     use_basin  logical flag governing use of Gaussian basin
!     use_wall   logical flag governing use of droplet boundary
!
!
#include "tinker_macro.h"
module kgeoms
   implicit none
   integer npfix
   integer ndfix
   integer nafix
   integer ntfix
   integer ngfix
   integer nchir
   integer npfixloc,ndfixloc,nafixloc,ntfixloc,ngfixloc,nchirloc
   integer, pointer :: ipfix(:)
   integer, pointer :: kpfix(:,:), idfix(:,:), iafix(:,:)
   integer, pointer :: itfix(:,:), igfix(:,:), ichir(:,:)
   integer :: winipfix,winkpfix,winidfix,winiafix
   integer :: winitfix,winigfix,winichir
   integer :: winxpfix,winypfix,winzpfix
   integer :: wintfix,wingfix,winchir
   integer :: winpfix,windfix,winafix
   real(t_p), pointer ::  xpfix(:),ypfix(:),zpfix(:)
   real(t_p), pointer ::  pfix(:,:),dfix(:,:),afix(:,:)
   real(t_p), pointer ::  tfix(:,:),gfix(:,:),chir(:,:)
   real(t_p) depth,width,rwall
   logical use_basin,use_wall
end
