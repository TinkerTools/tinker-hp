c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module katoms    --  forcefield parameters for the atom types  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     weight     average atomic mass of each atom type
c     atmcls     atom class number for each of the atom types
c     atmnum     atomic number for each of the atom types
c     ligand     number of atoms to be attached to each atom type
c     symbol     modified atomic symbol for each atom type
c     describe   string identifying each of the atom types
c
c
#include "tinker_precision.h"
      module katoms
      use sizes ,only: maxtyp
      implicit none
      integer atmcls(maxtyp),atmnum(maxtyp)
      integer ligand(maxtyp)
      real(r_p) weight(maxtyp)
      character*3 symbol(maxtyp)
      character*24 describe(maxtyp)
      save
      end
