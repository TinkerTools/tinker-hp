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
c     i_[xyz]cell  inverse length of the a-axis of the complete replicated cell
c     ycell    length of the b-axis of the complete replicated cell
c     zcell    length of the c-axis of the complete replicated cell
c     xcell2   half the length of the a-axis of the replicated cell
c     ycell2   half the length of the b-axis of the replicated cell
c     zcell2   half the length of the c-axis of the replicated cell
c     eps_cell contains an epsilon value such that max(xcell,ycell,zcell)+eps_cell > max(xcell,ycell,zcell) 
c
c
#include "tinker_precision.h"
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
