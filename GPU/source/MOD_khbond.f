c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ## module  khbond  --  forcefield parameters for H-bonding terms  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     maxnhb   maximum number of hydrogen bonding pair entries
c
c     radhb    radius parameter for hydrogen bonding pairs
c     epshb    well depth parameter for hydrogen bonding pairs
c     khb      string of atom types for hydrogen bonding pairs
c
c
#include "tinker_precision.h"
      module khbond
      implicit none
      integer maxnhb
      parameter (maxnhb=500)
      real(t_p) radhb(maxnhb),epshb(maxnhb)
      character*8 khb(maxnhb)
      save
      end
