c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  bond.i  --  covalent bonds in the current structure  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     bk      bond stretch force constants (kcal/mole/Ang**2)
c     bl      ideal bond length values in Angstroms
c     nbond   total number of bond stretches in the system
c     ibnd    numbers of the atoms in each bond stretch
c
c
      integer nbond,nbondloc
      integer, pointer :: ibnd(:,:)
      real*8, pointer ::  bk(:),bl(:)
      common /bond/  bk,bl,nbond,ibnd,nbondloc


