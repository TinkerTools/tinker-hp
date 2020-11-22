c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module cell  --  periodic boundaries using replicated cells  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     xcell    length of the a-axis of the complete replicated cell
c     ycell    length of the b-axis of the complete replicated cell
c     zcell    length of the c-axis of the complete replicated cell
c     xcell2   half the length of the a-axis of the replicated cell
c     ycell2   half the length of the b-axis of the replicated cell
c     zcell2   half the length of the c-axis of the replicated cell
c
c
      module cell
      implicit none
      real*8 xcell,ycell,zcell
      real*8 xcell2,ycell2,zcell2
      save
      end
