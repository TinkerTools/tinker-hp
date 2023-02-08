c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module urypot  --  specifics of Urey-Bradley functional form  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     cury       cubic coefficient in Urey-Bradley potential
c     qury       quartic coefficient in Urey-Bradley potential
c     ureyunit   convert Urey-Bradley energy to kcal/mole
c
c
#include "tinker_macro.h"
      module urypot
      implicit none
      real(t_p) cury,qury
      real(t_p) ureyunit
      save
      end
