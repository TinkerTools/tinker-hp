c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module kgeoms  --  parameters for the geometrical restraints  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     xpfix      x-coordinate target for each restrained position
c     winxpfix    window object corresponding to xpfix
c     ypfix      y-coordinate target for each restrained position
c     winypfix    window object corresponding to ypfix
c     zpfix      z-coordinate target for each restrained position
c     winzpfix    window object corresponding to zpfix
c     pfix       force constant and flat-well range for each position
c     winpfix    window object corresponding to pfix
c     dfix       force constant and target range for each distance
c     windfix    window object corresponding to dfix
c     afix       force constant and target range for each angle
c     winafix    window object corresponding to afix
c     tfix       force constant and target range for each torsion
c     wintfix    window object corresponding to tfix
c     gfix       force constant and target range for each group distance
c     wingfix    window object corresponding to gfix
c     chir       force constant and target range for chiral centers
c     winchir    window object corresponding to chir
c     depth      depth of shallow Gaussian basin restraint
c     width      exponential width coefficient of Gaussian basin
c     rwall      radius of spherical droplet boundary restraint
c     npfix      number of position restraints to be applied
c     ipfix      atom number involved in each position restraint
c     winipfix    window object corresponding to ipfix
c     kpfix      flags to use x-, y-, z-coordinate position restraints
c     winkpfix    window object corresponding to kpfix
c     ndfix      number of distance restraints to be applied
c     idfix      atom numbers defining each distance restraint
c     winidfix    window object corresponding to idfix
c     nafix      number of angle restraints to be applied
c     iafix      atom numbers defining each angle restraint
c     winiafix    window object corresponding to iafix
c     ntfix      number of torsional restraints to be applied
c     itfix      atom numbers defining each torsional restraint
c     winitfix    window object corresponding to itfix
c     ngfix      number of group distance restraints to be applied
c     igfix      group numbers defining each group distance restraint
c     winigfix    window object corresponding to igfix
c     nchir      number of chirality restraints to be applied
c     ichir      atom numbers defining each chirality restraint
c     winichir    window object corresponding to ichir
c     use_basin  logical flag governing use of Gaussian basin
c     use_wall   logical flag governing use of droplet boundary
c
c
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
