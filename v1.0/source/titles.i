c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  titles.i  --  title for the current molecular system  ##
c     ##                                                        ##
c     ############################################################
c
c
c     ltitle   length in characters of the nonblank title string
c     title    title used to describe the current structure
c
c
      integer ltitle
      character*120 title
      common /titles/ ltitle,title
