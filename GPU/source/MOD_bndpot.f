c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module bndpot  --  specifics of bond stretch functional forms  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     cbnd      cubic coefficient in bond stretch potential
c     qbnd      quartic coefficient in bond stretch potential
c     bndunit   convert bond stretch energy to kcal/mole
c     bndtyp    type of bond stretch potential energy function
c
c
#include "tinker_precision.h"
      module bndpot
      implicit none
      real(t_p) cbnd,qbnd
      real(t_p) bndunit
      character*8 bndtyp
      integer bndtyp_i
      enum,bind(C)
      enumerator BND_HARMONIC
      enumerator BND_MORSE
      enumerator BND_NO_TYPE
      end enum
      save
      end
