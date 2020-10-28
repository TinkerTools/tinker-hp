c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module strtor  --  stretch-torsions in the current structure  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     kst       1-, 2- and 3-fold stretch-torsion force constants
c     winkst    window object corresponding to kst
c     nstrtor   total number of stretch-torsion interactions
c     nstrtorloc   local number of stretch-torsion interactions
c     nbstrtor   number of stretch-torsion interactions before each atom
c     winnbstrtor    window object corresponding to nbstrtor
c     ist       torsion and bond numbers used in stretch-torsion
c     winist    window object corresponding to ist
c
c
#include "tinker_precision.h"
      module strtor
      implicit none
      integer nstrtor,nstrtorloc
      integer, pointer :: ist(:,:),nbstrtor(:)
      real(t_p), pointer :: kst(:,:)
      integer :: winist,winnbstrtor,winkst
      end
