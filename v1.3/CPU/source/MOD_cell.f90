!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module cell  --  periodic boundaries using replicated cells  ##
!     ##                                                               ##
!     ###################################################################
!
!
module cell
   implicit none
   real*8 :: xcell !<length of the a-axis of the complete replicated cell
   real*8 :: ycell !<length of the b-axis of the complete replicated cell
   real*8 :: zcell !<length of the c-axis of the complete replicated cell
   real*8 :: xcell2 !<half the length of the a-axis of the replicated cell
   real*8 :: ycell2 !<half the length of the b-axis of the replicated cell
   real*8 :: zcell2 !<half the length of the c-axis of the replicated cell
   save
end
