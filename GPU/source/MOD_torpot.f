c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module torpot  --  specifics of torsional functional forms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     idihunit  convert improper dihedral energy to kcal/mole
c     itorunit  convert improper torsion amplitudes to kcal/mole
c     torsunit  convert torsional parameter amplitudes to kcal/mole
c     ptorunit  convert pi-orbital torsion energy to kcal/mole
c     storunit  convert stretch-torsion energy to kcal/mole
c     ttorunit  convert stretch-torsion energy to kcal/mole
c
c
#include "tinker_precision.h"
      module torpot
      implicit none
      real(t_p) idihunit,itorunit,torsunit
      real(t_p) ptorunit,storunit,ttorunit
      save
      end
