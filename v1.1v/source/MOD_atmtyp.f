c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module atmtyp  --  atomic properties for each current atom  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     mass      atomic weight for each atom in the system
c     tag       integer atom labels from input coordinates file
c     class     atom class number for each atom in the system
c     atomic    atomic number for each atom in the system
c     valence   valence number for each atom in the system
c     name      atom name for each atom in the system
c     story     descriptive type for each atom in system
c
c
      module atmtyp
      implicit none
      integer, pointer :: tag(:),class(:)
      integer, pointer :: atomic(:),valence(:)
      real*8, pointer :: mass(:)
      character*3, pointer :: name(:)
      character*24, pointer :: story(:)
      save
      end
