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
c     atorunit  convert angle-torsion energy to kcal/mole
c     ttorunit  convert stretch-torsion energy to kcal/mole
c
c
      module torpot
      implicit none
      real*8 idihunit,itorunit,torsunit
      real*8 ptorunit,storunit,ttorunit
      real*8 atorunit
      save
      end
