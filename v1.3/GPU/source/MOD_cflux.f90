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
!     nbflx   total number of bond charge flux interactions
!     naflx   total number of angle charge flux interactions
!     bflx    bond stretching charge flux constant (electrons/Ang)
!     winbflx window object corresponding to bflx
!     aflx    angle bending charge flux constant (electrons/radian)
!     winaflx window object corresponding to aflx
!     abflx   asymmetric stretch charge flux constant (electrons/Ang)
!     winabflx window object corresonding to abflx
!
!
#include "tinker_macro.h"
module cflux
   use mdevice, only: channel_t
   implicit none
   integer nbflx
   integer naflx
   real(t_p), pointer :: bflx(:)
   real(t_p), pointer :: aflx(:,:)
   real(t_p), pointer :: abflx(:,:)
   real(t_p), allocatable, target :: cflux_pot(:),cflux_potrec(:)
   integer winbflx,winaflx,winabflx

   interface
      module subroutine adflux1(dcf,de,channel)
         real(r_p) dcf(*),de(*)
         type(channel_t), intent(in) :: channel
      end subroutine
      module subroutine adflux2(dcf,de,channel)
         mdyn_rtyp dcf(*),de(*)
         type(channel_t), intent(in) :: channel
      end subroutine
   end interface

   interface
      module subroutine dcflux1 (pot,dcf,channel)
         real(t_p) pot(*)
         real(r_p) dcf(*)
         type(channel_t), intent(in) :: channel
      end subroutine
      module subroutine dcflux2 (pot,dcf,channel)
         real(t_p),intent(in   ):: pot(*)
         mdyn_rtyp,intent(inout):: dcf(*)
         type(channel_t), intent(in) :: channel
      end subroutine
   end interface

   interface
      module subroutine release_cflux_workspaces(queue)
         integer,intent(in):: queue
      end subroutine
   end interface
end
