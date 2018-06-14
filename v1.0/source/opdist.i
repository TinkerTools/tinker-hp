c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  opdist.i  --  out-of-plane distances in current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     opdk      force constant values for out-of-plane distance
c     nopdist   total number of out-of-plane distances in the system
c     iopb      numbers of the atoms in each out-of-plane distance
c
c
      integer nopdist,nopdistloc
      integer, pointer :: iopd(:,:),nbopdist(:)
      real*8, pointer ::  opdk(:)
      common /opdist/ opdk,nopdist,iopd,nopdistloc,nbopdist
