c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module repel  --  Pauli repulsion for current structure  ##        
c     ##                                                           ##
c     ###############################################################
c
c
c     nrep      total number of repulsion sites in the system
c     nreploc    local number of repulsion sites in the system
c     sizpr     Pauli repulsion size parameter value at each site
c     dmppr     Pauli repulsion alpha damping value at each site
c     elepr     Pauli repulsion valence electrons at each site
c     winsizepr window object corresponding to sizepr
c     windmppr window object corresponding to damppr
c     winelepr window object corresponding to elepr
c
c
#include "tinker_macro.h"
      module repel
      implicit none
      integer nrep,nreploc
      real(t_p), pointer :: sizpr(:)
      real(t_p), pointer :: dmppr(:)
      real(t_p), pointer :: elepr(:)
      integer :: winsizpr,windmppr,winelepr
      end
