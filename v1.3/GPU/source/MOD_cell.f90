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
!     xcell    length of the a-axis of the complete replicated cell
!     i_[xyz]cell  inverse length of the a-axis of the complete replicated cell
!     ycell    length of the b-axis of the complete replicated cell
!     zcell    length of the c-axis of the complete replicated cell
!     xcell2   half the length of the a-axis of the replicated cell
!     ycell2   half the length of the b-axis of the replicated cell
!     zcell2   half the length of the c-axis of the replicated cell
!     eps_cell contains an epsilon value such that max(xcell,ycell,zcell)+eps_cell > max(xcell,ycell,zcell)
!
!
#include "tinker_macro.h"
module cell
   implicit none
   real(t_p) xcell,ycell,zcell
   real(t_p) xcell2,ycell2,zcell2
   real(t_p) i_xcell,i_ycell,i_zcell
   real(t_p) eps_cell

!$acc declare create(xcell,ycell,zcell)
!$acc declare create(xcell2,ycell2,zcell2)
!$acc declare create(i_xcell,i_ycell,i_zcell)
!$acc declare create(eps_cell)
end
