c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module polgrp  --  polarizable site group connectivity lists   ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     maxp11   maximum number of atoms in a polarization group
c     maxp12   maximum number of atoms in groups 1-2 to an atom
c     maxp13   maximum number of atoms in groups 1-3 to an atom
c     maxp14   maximum number of atoms in groups 1-4 to an atom
c     maxscalp maximum number of atoms in scaled polar interaction per atom
c
c     np11     number of atoms in polarization group of each atom
c     winnp11    window object corresponding to np11
c     ip11     atom numbers of atoms in same group as each atom
c     winip11    window object corresponding to ip11
c     np12     number of atoms in groups 1-2 to each atom
c     winnp12    window object corresponding to np12
c     ip12     atom numbers of atoms in groups 1-2 to each atom
c     winip12    window object corresponding to ip12
c     np13     number of atoms in groups 1-3 to each atom
c     winnp13    window object corresponding to np13
c     ip13     atom numbers of atoms in groups 1-3 to each atom
c     winip13    window object corresponding to ip13
c     np14     number of atoms in groups 1-4 to each atom
c     winnp14    window object corresponding to np14
c     ip14     atom numbers of atoms in groups 1-4 to each atom
c     winip14    window object corresponding to ip14
c
c
      module polgrp
      use sizes
      implicit none
      integer maxp11,maxp12
      integer maxp13,maxp14
      integer maxscalp
      parameter (maxp11=200)
      parameter (maxp12=200)
      parameter (maxp13=200)
      parameter (maxp14=200)
      parameter (maxscalp=480)
      integer, pointer :: np11(:),ip11(:,:),np12(:),ip12(:,:)
      integer winnp11,winip11,winnp12,winip12
      integer, pointer :: np13(:),ip13(:,:),np14(:),ip14(:,:)
      integer winnp13,winip13,winnp14,winip14
      save
      end
