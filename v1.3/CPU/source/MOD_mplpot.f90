!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module mplpot  --  specifics of atomic multipole functions  ##
!     ##                                                              ##
!     ##################################################################
!
!
!
module mplpot
   implicit none
   real*8 :: m2scale !<Scale factor for 1-2 multipole energy interactions
   real*8 :: m3scale !<Scale factor for 1-3 multipole energy interactions
   real*8 :: m4scale !<Scale factor for 1-4 multipole energy interactions
   real*8 :: m5scale !<Scale factor for 1-5 multipole energy interactions
   character*7 :: pentyp !<type of penetration damping (NONE, GORDON1, GORDON2)
   save
end
