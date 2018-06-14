c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  atmtyp.i  --  atomic properties for each current atom  ##
c     ##                                                         ##
c     #############################################################
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
      integer, pointer :: tag(:),class(:)
      integer, pointer :: atomic(:),valence(:)
      real*8, pointer :: mass(:)
      character*3, pointer :: name(:)
      character*24, pointer :: story(:)
      common /atmtyp/ mass,tag,class,
     &                atomic,valence,name,story
