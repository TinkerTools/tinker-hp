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
c     winmass    window object corresponding to mass
c     tag       integer atom labels from input coordinates file
c     wintag    window object corresponding to tag
c     class     atom class number for each atom in the system
c     winclass    window object corresponding to class
c     atomic    atomic number for each atom in the system
c     winatomic    window object corresponding to atomic
c     valence   valence number for each atom in the system
c     winvalence    window object corresponding to valence
c     name      atom name for each atom in the system
c     winname    window object corresponding to name
c     story     descriptive type for each atom in system
c     winstory    window object corresponding to story
c
c
      module atmtyp
      implicit none
      integer, allocatable :: tag(:)
      integer, pointer :: class(:)
      integer, pointer :: atomic(:),valence(:)
      integer :: winclass,winatomic,winvalence
      real*8, pointer :: mass(:)
      character*3, allocatable :: name(:)
      character*24, pointer :: story(:)
      integer :: winmass,winstory
      save
      end
