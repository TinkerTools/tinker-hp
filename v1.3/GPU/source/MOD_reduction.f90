!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module reduction  --  GPU reduction storage               ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
module reduction
   use tinheader, only: ti_p,re_p
   use lifecycle_observer, only: LifecycleReduction, record_lifecycle_event
   implicit none

   integer,parameter :: RED_BUFF_SIZE=ishft(1,15)
   integer   nred_buff(RED_BUFF_SIZE)
   real(t_p) vred_buff(6*RED_BUFF_SIZE)
   real(t_p) vred_buf1(6*RED_BUFF_SIZE)
   real(r_p) ered_buff(RED_BUFF_SIZE)
   real(r_p),target :: lam_buff(RED_BUFF_SIZE)
   ener_rtyp ered_buf1(RED_BUFF_SIZE)
   logical :: reduction_data_on_device=.false.

   interface reduce_buffer
      module procedure reduce_buffer_m
   end interface
   interface reduce_energy_virial
      module procedure reduce_energy_virial0
      module procedure reduce_energy_virial1
      module procedure reduce_energy_virial2
   end interface
   interface reduce_energy_action
      module procedure reduce_energy_action0
      module procedure reduce_energy_action1
   end interface

contains

   subroutine create_reduction_buffers(queue)
      implicit none
      integer,intent(in):: queue

      if (.not.reduction_data_on_device) then
!$acc enter data create(ered_buff,ered_buf1,lam_buff,vred_buff,vred_buf1,nred_buff)
         reduction_data_on_device=.true.
      end if
      call zero_mdred_buffers(queue)
   end subroutine

   subroutine delete_reduction_buffers(queue)
      implicit none
      integer,intent(in):: queue

      if (.not.reduction_data_on_device) return
      call record_lifecycle_event(LifecycleReduction)
!$acc wait(queue)
!$acc exit data delete(ered_buff,ered_buf1,lam_buff,vred_buff,vred_buf1,nred_buff)
      reduction_data_on_device=.false.
   end subroutine

   subroutine reduce_buffer_m(r_buff,sz_buff,redVal,queue)
      implicit none
      integer,intent(in):: sz_buff,queue
      real(r_p) r_buff(*)
      real(r_p),intent(out):: redVal
      integer i

!$acc parallel loop present(redVal,r_buff) async(queue)
      do i = 1,sz_buff
         redVal = redVal + r_buff(i)
         r_buff(i) = 0
      end do
   end subroutine

   subroutine reduce_energy_virial0(energy,vxx,vxy,vxz,vyy,vyz,vzz,queue)
      implicit none
      integer,intent(in):: queue
      real(r_p) energy,vxx,vxy,vxz,vyy,vyz,vzz
      integer i

!$acc parallel loop async(queue) present(energy,vxx,vxy,vxz,vyy,vyz,vzz,ered_buff,vred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_buff(i)
         vxx = vxx + real(vred_buff(i),r_p)
         vxy = vxy + real(vred_buff(RED_BUFF_SIZE+i),r_p)
         vxz = vxz + real(vred_buff(2*RED_BUFF_SIZE+i),r_p)
         vyy = vyy + real(vred_buff(3*RED_BUFF_SIZE+i),r_p)
         vyz = vyz + real(vred_buff(4*RED_BUFF_SIZE+i),r_p)
         vzz = vzz + real(vred_buff(5*RED_BUFF_SIZE+i),r_p)
         ered_buff(i) = 0
         vred_buff(i) = 0
         vred_buff(RED_BUFF_SIZE+i) = 0
         vred_buff(2*RED_BUFF_SIZE+i) = 0
         vred_buff(3*RED_BUFF_SIZE+i) = 0
         vred_buff(4*RED_BUFF_SIZE+i) = 0
         vred_buff(5*RED_BUFF_SIZE+i) = 0
      end do
   end subroutine

   subroutine reduce_energy_virial1(energy,vxx,vxy,vxz,vyy,vyz,vzz,ered_b,queue)
      implicit none
      integer,intent(in):: queue
      ener_rtyp,intent(inout):: ered_b(RED_BUFF_SIZE)
      ener_rtyp energy
      real(r_p) vxx,vxy,vxz,vyy,vyz,vzz
      integer i

!$acc parallel loop async(queue) present(energy,vxx,vxy,vxz,vyy,vyz,vzz,ered_b,vred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_b(i)
         vxx = vxx + real(vred_buff(i),r_p)
         vxy = vxy + real(vred_buff(RED_BUFF_SIZE+i),r_p)
         vxz = vxz + real(vred_buff(2*RED_BUFF_SIZE+i),r_p)
         vyy = vyy + real(vred_buff(3*RED_BUFF_SIZE+i),r_p)
         vyz = vyz + real(vred_buff(4*RED_BUFF_SIZE+i),r_p)
         vzz = vzz + real(vred_buff(5*RED_BUFF_SIZE+i),r_p)
         ered_b(i) = 0
         vred_buff(i) = 0
         vred_buff(RED_BUFF_SIZE+i) = 0
         vred_buff(2*RED_BUFF_SIZE+i) = 0
         vred_buff(3*RED_BUFF_SIZE+i) = 0
         vred_buff(4*RED_BUFF_SIZE+i) = 0
         vred_buff(5*RED_BUFF_SIZE+i) = 0
      end do
   end subroutine

   subroutine reduce_energy_virial2(energy,vxx,vxy,vxz,vyy,vyz,vzz,ered_b,vred_b,queue)
      implicit none
      integer,intent(in):: queue
      real(r_p),intent(inout):: ered_b(RED_BUFF_SIZE)
      real(t_p),intent(inout):: vred_b(6*RED_BUFF_SIZE)
      real(r_p) energy,vxx,vxy,vxz,vyy,vyz,vzz
      integer i

!$acc parallel loop async(queue) present(energy,vxx,vxy,vxz,vyy,vyz,vzz,ered_b,vred_b)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_b(i)
         vxx = vxx + real(vred_b(i),r_p)
         vxy = vxy + real(vred_b(RED_BUFF_SIZE+i),r_p)
         vxz = vxz + real(vred_b(2*RED_BUFF_SIZE+i),r_p)
         vyy = vyy + real(vred_b(3*RED_BUFF_SIZE+i),r_p)
         vyz = vyz + real(vred_b(4*RED_BUFF_SIZE+i),r_p)
         vzz = vzz + real(vred_b(5*RED_BUFF_SIZE+i),r_p)
         ered_b(i) = 0
         vred_b(i) = 0
         vred_b(RED_BUFF_SIZE+i) = 0
         vred_b(2*RED_BUFF_SIZE+i) = 0
         vred_b(3*RED_BUFF_SIZE+i) = 0
         vred_b(4*RED_BUFF_SIZE+i) = 0
         vred_b(5*RED_BUFF_SIZE+i) = 0
      end do
   end subroutine

   subroutine reduce_energy_action0(energy,n,queue)
      implicit none
      integer,intent(in):: queue
      real(r_p) energy
      integer n,i

!$acc parallel loop async(queue) present(energy,n,ered_buff,nred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_buff(i)
         n = n + nred_buff(i)
         ered_buff(i) = 0
         nred_buff(i) = 0
      end do
   end subroutine

   subroutine reduce_energy_action1(energy,n,ered_b,queue)
      implicit none
      integer,intent(in):: queue
      ener_rtyp ered_b(RED_BUFF_SIZE)
      ener_rtyp energy
      integer n,i

!$acc parallel loop async(queue) present(energy,n,ered_b,nred_buff)
      do i = 1,RED_BUFF_SIZE
         energy = energy + ered_b(i)
         n = n + nred_buff(i)
         ered_b(i) = 0
         nred_buff(i) = 0
      end do
   end subroutine

   subroutine zero_mdred_buffers(queue)
      implicit none
      integer,intent(in):: queue
      integer i

!$acc parallel loop async(queue) default(present)
      do i = 1,6*RED_BUFF_SIZE
         if (i < RED_BUFF_SIZE+1) then
            ered_buff(i) = 0
            ered_buf1(i) = 0
            nred_buff(i) = 0
         end if
         vred_buff(i) = 0
         vred_buf1(i) = 0
      end do
   end subroutine

   subroutine zero_evir_red_buffer(queue)
      implicit none
      integer,intent(in):: queue
      integer i

!$acc parallel loop async(queue) present(ered_buff,ered_buf1,vred_buff)
      do i = 1,6*RED_BUFF_SIZE
         if (i < RED_BUFF_SIZE+1) then
            ered_buff(i) = 0
            ered_buf1(i) = 0
         end if
         vred_buff(i) = 0
      end do
   end subroutine

   subroutine zero_en_red_buffer(queue)
      implicit none
      integer,intent(in):: queue
      integer i

!$acc parallel loop async(queue) present(ered_buff,ered_buf1,nred_buff)
      do i = 1,RED_BUFF_SIZE
         ered_buff(i) = 0
         ered_buf1(i) = 0
         nred_buff(i) = 0
      end do
   end subroutine

end module reduction
