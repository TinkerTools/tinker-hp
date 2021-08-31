c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module couple  --  near-neighbor atom connectivity lists   ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxn13   maximum number of atoms 1-3 connected to an atom
c     maxn14   maximum number of atoms 1-4 connected to an atom
c     maxn15   maximum number of atoms 1-5 connected to an atom
c
c     n12      number of atoms directly bonded to each atom
c     winn12    window object corresponding to n12
c     i12      atom numbers of atoms 1-2 connected to each atom
c     wini12    window object corresponding to i12
c     n13      number of atoms in a 1-3 relation to each atom
c     winn13    window object corresponding to n13
c     i13      atom numbers of atoms 1-3 connected to each atom
c     wini13    window object corresponding to i13
c     n14      number of atoms in a 1-4 relation to each atom
c     winn14    window object corresponding to n14
c     i14      atom numbers of atoms 1-4 connected to each atom
c     wini14    window object corresponding to i14
c     n15      number of atoms in a 1-5 relation to each atom
c     winn15    window object corresponding to n15
c     i15      atom numbers of atoms 1-5 connected to each atom
c     wini15    window object corresponding to i15
c
c
      module couple
      use sizes
      implicit none
      integer maxn13,maxn14,maxn15
      parameter (maxn13=3*maxvalue)
      parameter (maxn14=3*maxvalue)
      parameter (maxn15=3*maxvalue)
      integer, allocatable :: n12(:),i12(:,:)
      integer, pointer :: n13(:),i13(:,:)
      integer, pointer ::  n14(:),i14(:,:),n15(:),i15(:,:)
      integer :: winn13,wini13
      integer :: winn14,wini14,winn15,wini15
      save
      end
