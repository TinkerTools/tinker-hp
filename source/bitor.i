c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  bitor.i  --  bitorsions within the current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     nbitor  total number of bitorsions in the system
c     ibitor  numbers of the atoms in each bitorsion
c
c
      integer nbitor,nbitorloc
      integer, pointer :: ibitor(:,:),nbbitors(:)
      common /bitor/ nbitor,ibitor,nbitorloc,nbbitors
