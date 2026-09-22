!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  subroutine prtdynbeads  --  output of MD restart information  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     "prtdynbeads" writes out the information needed to restart a
!     molecular dynamics trajectory to an external disk file while
!     using PIMD
!
!
#include "tinker_precision.h"
subroutine prtdynbeads(polymer,ibead_beg,ibead_end)
   use atomsMirror, only: n
   use atmtyp, only: mass
   use units, only: convert
   use beads
   use boxes
   use files
   use restart, only: write_restart
   use status, only: failed, status_t
   use titles
   use timestat ,only:timer_io,timer_enter,timer_exit,quiet_timers
   implicit none
   type(POLYMER_COMM_TYPE), intent(in) :: polymer
   integer, intent(in) :: ibead_beg,ibead_end
   integer i,ibead
   real(r_p) :: wt
   character*3 numberbeads
   character*240 dynfile
   real(r_p), allocatable :: a_bead(:,:),aalt_bead(:,:)
   type(status_t) :: status_value
!
   allocate(a_bead(3,n),aalt_bead(3,n))
   call timer_enter(timer_io)
!
!     update an existing restart file or open a new one
!
   do ibead = ibead_beg,ibead_end

      write(numberbeads, '(i3.3)') ibead
      dynfile = filename(1:leng)//'_beads'//numberbeads//'.dyn'
      do i = 1,n
         wt = convert/mass(i)
         a_bead(:,i) = polymer%forces(:,i,ibead)*wt
         if (allocated(polymer%forces_slow)) then
            aalt_bead(:,i) = polymer%forces_slow(:,i,ibead)*wt
         else
            aalt_bead(:,i) = 0.0_re_p
         endif
      enddo
      call write_restart(trim(dynfile),title,xbox,ybox,zbox,alpha,beta,gamma,&
         &polymer%pos(1,:,ibead),polymer%pos(2,:,ibead),polymer%pos(3,:,ibead),&
         &polymer%vel(:,:,ibead),a_bead,aalt_bead,status_value)
      if (failed(status_value)) then
         write (6,'(/,a)') trim(status_value%message)
         call fatal
      endif

   enddo
   deallocate(a_bead,aalt_bead)
   call timer_exit(timer_io)
end
