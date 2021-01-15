c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module titles  --  title for the current molecular system  ##
c     ##                                                             ##
c     #################################################################
c
c
c     ltitle   length in characters of the nonblank title string
c     title    title used to describe the current structure
c
c
      module titles
      implicit none
      integer ltitle
      character*240 title
      save
      end
