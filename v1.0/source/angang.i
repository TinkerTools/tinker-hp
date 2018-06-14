c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  angang.i  --  angle-angle terms in current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     kaa       force constant for angle-angle cross terms
c     nangang   total number of angle-angle interactions
c     iaa       angle numbers used in each angle-angle term
c
c
      integer nangang,nangangloc
      integer, pointer :: iaa(:,:),nbangang(:)
      real*8, pointer ::  kaa(:)
      common /angang/ kaa,nangang,iaa,nbangang
