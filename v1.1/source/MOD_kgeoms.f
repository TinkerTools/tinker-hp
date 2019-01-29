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
c     ypfix      y-coordinate target for each restrained position
c     zpfix      z-coordinate target for each restrained position
c     pfix       force constant and flat-well range for each position
c     dfix       force constant and target range for each distance
c     afix       force constant and target range for each angle
c     tfix       force constant and target range for each torsion
c     gfix       force constant and target range for each group distance
c     chir       force constant and target range for chiral centers
c     depth      depth of shallow Gaussian basin restraint
c     width      exponential width coefficient of Gaussian basin
c     rwall      radius of spherical droplet boundary restraint
c     npfix      number of position restraints to be applied
c     ipfix      atom number involved in each position restraint
c     kpfix      flags to use x-, y-, z-coordinate position restraints
c     ndfix      number of distance restraints to be applied
c     idfix      atom numbers defining each distance restraint
c     nafix      number of angle restraints to be applied
c     iafix      atom numbers defining each angle restraint
c     ntfix      number of torsional restraints to be applied
c     itfix      atom numbers defining each torsional restraint
c     ngfix      number of group distance restraints to be applied
c     igfix      group numbers defining each group distance restraint
c     nchir      number of chirality restraints to be applied
c     ichir      atom numbers defining each chirality restraint
c     use_basin  logical flag governing use of Gaussian basin
c     use_wall   logical flag governing use of droplet boundary
c
c
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
      real*8, pointer ::  xpfix(:),ypfix(:),zpfix(:)
      real*8, pointer ::  pfix(:,:),dfix(:,:),afix(:,:)
      real*8, pointer ::  tfix(:,:),gfix(:,:),chir(:,:)
      real*8 depth,width,rwall
      logical use_basin,use_wall
      save
      end
