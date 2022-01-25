c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module mplpot  --  specifics of atomic multipole functions  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     m2scale   scale factor for 1-2 multipole energy interactions
c     m3scale   scale factor for 1-3 multipole energy interactions
c     m4scale   scale factor for 1-4 multipole energy interactions
c     m5scale   scale factor for 1-5 multipole energy interactions
c     pentyp       type of penetration damping (NONE, GORDON1, GORDON2)
c
c
      module mplpot
      implicit none
      real*8 m2scale,m3scale
      real*8 m4scale,m5scale
      character*7 pentyp
      save
      end
