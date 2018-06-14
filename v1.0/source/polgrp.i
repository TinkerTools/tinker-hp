c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  polgrp.i  --  polarizable site group connectivity lists   ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxp11   maximum number of atoms in a polarization group
c     maxp12   maximum number of atoms in groups 1-2 to an atom
c     maxp13   maximum number of atoms in groups 1-3 to an atom
c     maxp14   maximum number of atoms in groups 1-4 to an atom
c
c     np11     number of atoms in polarization group of each atom
c     ip11     atom numbers of atoms in same group as each atom
c     np12     number of atoms in groups 1-2 to each atom
c     ip12     atom numbers of atoms in groups 1-2 to each atom
c     np13     number of atoms in groups 1-3 to each atom
c     ip13     atom numbers of atoms in groups 1-3 to each atom
c     np14     number of atoms in groups 1-4 to each atom
c     ip14     atom numbers of atoms in groups 1-4 to each atom
c
c
      integer maxp11,maxp12
      integer maxp13,maxp14
      parameter (maxp11=40)
      parameter (maxp12=40)
      parameter (maxp13=40)
      parameter (maxp14=40)
      integer, pointer :: np11(:),ip11(:,:),np12(:),ip12(:,:)
      integer, pointer :: np13(:),ip13(:,:),np14(:),ip14(:,:)
      common /polgrp/ np11,ip11,
     &                np12,ip12,
     &                np13,ip13,
     &                np14,ip14
