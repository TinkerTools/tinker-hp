c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module krepl  --  Pauli repulsion forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     prsiz     Pauli repulsion size value for each atom class
c     prdmp     alpha Pauli repulsion parameter for each atom class
c     prele     number of valence electrons for each atom class
c
c
      module krepl
      use sizes
      implicit none
      real*8 prsiz(maxtyp),prdmp(maxtyp),prele(maxtyp)
      save
      end
