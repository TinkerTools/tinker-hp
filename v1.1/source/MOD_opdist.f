c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module opdist  --  out-of-plane distances in current structure  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     nopdist   total number of out-of-plane distances in the system
c     nopdistloc   local number of out-of-plane distances in the system
c     iopd      numbers of the atoms in each out-of-plane distance
c     nbopdist number of angle used in out-of-plane distance before each atom
c     opdk      force constant values for out-of-plane distance
c
c
      module opdist
      implicit none
      integer nopdist,nopdistloc
      integer, pointer :: iopd(:,:),nbopdist(:)
      real*8, pointer ::  opdk(:)
      save
      end
