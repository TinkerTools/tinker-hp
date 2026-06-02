!
!
!     ##########################################################
!     ##  COPYRIGHT (C) 2020 by Chengwen Liu & Jay W. Ponder  ##
!     ##                 All Rights Reserved                  ##
!     ##########################################################
!
!     #############################################################
!     ##                                                         ##
!     ##  module cflux  --  charge flux terms in current system  ##
!     ##                                                         ##
!     #############################################################
!
!
module cflux
   implicit none
   integer :: nbflx !<total number of bond charge flux interactions
   integer :: naflx !<total number of angle charge flux interactions
   integer :: winbflx !<window object corresponding to bflx
   integer :: winaflx !<window object corresponding to aflx
   integer :: winabflx !<window object corresponding to abflx
   real*8, pointer :: bflx(:) !<bond stretching charge flux constant (electrons/Ang)
   real*8, pointer :: aflx(:,:) !<angle bending charge flux constant (electrons/radian)
   real*8, pointer :: abflx(:,:) !<asymmetric stretch charge flux constant (electrons/Ang)
   save
end
