c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  couple.i  --  near-neighbor atom connectivity lists   ##
c     ##                                                        ##
c     ############################################################
c
c
c     maxn13   maximum number of atoms 1-3 connected to an atom
c     maxn14   maximum number of atoms 1-4 connected to an atom
c     maxn15   maximum number of atoms 1-5 connected to an atom
c
c     n12      number of atoms directly bonded to each atom
c     i12      atom numbers of atoms 1-2 connected to each atom
c     n13      number of atoms in a 1-3 relation to each atom
c     i13      atom numbers of atoms 1-3 connected to each atom
c     n14      number of atoms in a 1-4 relation to each atom
c     i14      atom numbers of atoms 1-4 connected to each atom
c     n15      number of atoms in a 1-5 relation to each atom
c     i15      atom numbers of atoms 1-5 connected to each atom
c
c
      integer maxn13,maxn14,maxn15
      parameter (maxn13=3*maxval)
      parameter (maxn14=3*maxval)
      parameter (maxn15=3*maxval)
      integer, pointer :: n12(:),i12(:,:),n13(:),i13(:,:)
      integer, pointer ::  n14(:),i14(:,:),n15(:),i15(:,:)
      common /couple/ n12,i12,n13,
     &                i13,n14,i14,
     &                n15,i15
