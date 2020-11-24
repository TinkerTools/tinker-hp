c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module bond  --  covalent bonds in the current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     bk      bond stretch force constants (kcal/mole/Ang**2)
c     winbk    window object corresponding to bk
c     bl      ideal bond length values in Angstroms
c     winbl    window object corresponding to bk
c     nbond   total number of bond stretches in the system
c     nbondloc   local number of bond stretches in the system
c     ibnd    numbers of the atoms in each bond stretch
c     winibnd    window object corresponding to ibnd
c
c
      module bond
      implicit none
      integer nbond,nbondloc
      integer, pointer :: ibnd(:,:)
      integer :: winibnd
      real*8, pointer ::  bk(:),bl(:)
      integer :: winbk,winbl
      save
      end
