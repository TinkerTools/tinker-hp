!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module mdtraj  --  trajectory persistence interface      ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
module mdtraj
   use dcdmod ,only: dcdinfo_t
   use status ,only: status_t
   implicit none

   integer,parameter :: arc_backend = 1
   integer,parameter :: dcd_backend = 2

   type trajectory_t
      ! Copy policy: noncopyable; use the lifecycle procedures only.
      integer :: backend = 0
      integer :: frame_count = 0
      logical :: is_active = .false.
      logical :: append_mode = .false.
      logical :: owns_dcd_file = .false.
      character(len=240) :: filename = ''
      type(dcdinfo_t) :: dcdinfo
   end type trajectory_t

   type frame_view_t
      ! Copy policy: value-copy of non-owning call-time aliases.
      integer :: atoms = 0
      logical :: host_global = .false.
      real(r_p),pointer :: x(:) => null(), y(:) => null(), z(:) => null()
   end type frame_view_t

   type(trajectory_t) :: dynamic_trajectory

   interface
      module subroutine open_trajectory(writer,backend,filename,append,&
         &status_value)
         type(trajectory_t),intent(inout) :: writer
         integer,intent(in) :: backend
         character(*),intent(in) :: filename
         logical,intent(in) :: append
         type(status_t),intent(out) :: status_value
      end subroutine open_trajectory

      module subroutine write_frame(writer,frame,step,dt,status_value)
         type(trajectory_t),intent(inout) :: writer
         type(frame_view_t),intent(in) :: frame
         integer,intent(in) :: step
         real(r_p),intent(in) :: dt
         type(status_t),intent(out) :: status_value
      end subroutine write_frame

      module subroutine write_snapshot(filename,frame,status_value,append)
         character(*),intent(in) :: filename
         type(frame_view_t),intent(in) :: frame
         type(status_t),intent(out) :: status_value
         logical,intent(in),optional :: append
      end subroutine write_snapshot

      module subroutine close_trajectory(writer,status_value)
         type(trajectory_t),intent(inout) :: writer
         type(status_t),intent(out) :: status_value
      end subroutine close_trajectory
   end interface
end module mdtraj
