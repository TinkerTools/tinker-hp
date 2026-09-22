!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
submodule (mdtraj) mdtraj_operations
   use dcdmod ,only: dcdio_write,dcdfile_close
   use status ,only: clear_status,set_status,status_invalid_argument,status_ok
   implicit none
   integer,external :: freeunit

   interface
      subroutine prtxyz(ixyz,x_frame,y_frame,z_frame,atom_count_frame)
         integer,intent(in) :: ixyz
         real(r_p),optional,intent(in) :: x_frame(:),y_frame(:),z_frame(:)
         integer,optional,intent(in) :: atom_count_frame
      end subroutine prtxyz
   end interface
contains

   module subroutine open_trajectory(writer,backend,filename,append,&
      &status_value)
      type(trajectory_t),intent(inout) :: writer
      integer,intent(in) :: backend
      character(*),intent(in) :: filename
      logical,intent(in) :: append
      type(status_t),intent(out) :: status_value

      call clear_status(status_value)
      if (backend /= arc_backend .and. backend /= dcd_backend) then
         call set_status(status_value,status_invalid_argument,&
            &'Trajectory backend must be ARC or DCD')
         return
      end if
      if (len_trim(filename) == 0 .or. writer%is_active) then
         call set_status(status_value,status_invalid_argument,&
            &'Trajectory filename is empty or a session is already active')
         return
      end if

      writer%backend = backend
      writer%filename = trim(filename)
      writer%append_mode = append
      writer%is_active = .true.
      writer%frame_count = 0
      writer%owns_dcd_file = .false.
      writer%dcdinfo%idcd = 0
      writer%dcdinfo%nframes = 0
      writer%dcdinfo%natoms = 0
      call set_status(status_value,status_ok,'Trajectory session opened')
   end subroutine open_trajectory

   module subroutine write_frame(writer,frame,step,dt,status_value)
      type(trajectory_t),intent(inout) :: writer
      type(frame_view_t),intent(in) :: frame
      integer,intent(in) :: step
      real(r_p),intent(in) :: dt
      type(status_t),intent(out) :: status_value
      integer :: unit,ios
      logical :: exist

      call clear_status(status_value)
      if (.not.writer%is_active .or. frame%atoms <= 0 .or. &
         &.not.frame%host_global .or. .not.associated(frame%x) .or. &
         &.not.associated(frame%y) .or. .not.associated(frame%z) .or. &
         &size(frame%x) /= frame%atoms .or. size(frame%y) /= frame%atoms .or. &
         &size(frame%z) /= frame%atoms .or. step < 0 .or. dt <= 0.0) then
         call set_status(status_value,status_invalid_argument,&
            &'Trajectory frame is not a complete host frame')
         return
      end if

      if (writer%backend == dcd_backend) then
         call dcdio_write(writer%dcdinfo,writer%filename,step,dt,&
            &writer%append_mode,frame%x,frame%y,frame%z)
         writer%owns_dcd_file = .true.
      else
         inquire(file=trim(writer%filename),exist=exist)
         unit = freeunit ()
         if (exist .and. writer%append_mode) then
            open(unit=unit,file=trim(writer%filename),status='old',&
               &position='append',action='write',iostat=ios)
         else
            open(unit=unit,file=trim(writer%filename),status='new',&
               &action='write',iostat=ios)
         end if
         if (ios /= 0) then
            call set_status(status_value,status_invalid_argument,&
               &'Unable to open ARC trajectory')
            return
         end if
         call prtxyz(unit,frame%x,frame%y,frame%z,frame%atoms)
         close(unit=unit,iostat=ios)
         if (ios /= 0) then
            call set_status(status_value,status_invalid_argument,&
               &'Unable to close ARC trajectory')
            return
         end if
      end if
      writer%frame_count = writer%frame_count + 1
      call set_status(status_value,status_ok,'Trajectory frame written')
   end subroutine write_frame

   module subroutine write_snapshot(filename,frame,status_value)
      character(*),intent(in) :: filename
      type(frame_view_t),intent(in) :: frame
      type(status_t),intent(out) :: status_value
      integer :: unit,ios

      call clear_status(status_value)
      if (len_trim(filename) == 0 .or. frame%atoms <= 0 .or. &
         &.not.frame%host_global .or. .not.associated(frame%x) .or. &
         &.not.associated(frame%y) .or. .not.associated(frame%z) .or. &
         &size(frame%x) /= frame%atoms .or. size(frame%y) /= frame%atoms .or. &
         &size(frame%z) /= frame%atoms) then
         call set_status(status_value,status_invalid_argument,&
            &'Snapshot needs a name and a complete host frame')
         return
      end if

      unit = freeunit ()
      open(unit=unit,file=trim(filename),status='new',action='write',&
         &iostat=ios)
      if (ios /= 0) then
         call set_status(status_value,status_invalid_argument,&
            &'Unable to open XYZ snapshot')
         return
      end if
      call prtxyz(unit,frame%x,frame%y,frame%z,frame%atoms)
      close(unit=unit,iostat=ios)
      if (ios /= 0) then
         call set_status(status_value,status_invalid_argument,&
            &'Unable to close XYZ snapshot')
         return
      end if
      call set_status(status_value,status_ok,'XYZ snapshot written')
   end subroutine write_snapshot

   module subroutine close_trajectory(writer,status_value)
      type(trajectory_t),intent(inout) :: writer
      type(status_t),intent(out) :: status_value

      call clear_status(status_value)
      if (writer%owns_dcd_file) call dcdfile_close(writer%dcdinfo)
      if (allocated(writer%dcdinfo%filename)) &
         &deallocate(writer%dcdinfo%filename)
      writer%backend = 0
      writer%frame_count = 0
      writer%is_active = .false.
      writer%append_mode = .false.
      writer%owns_dcd_file = .false.
      writer%filename = ''
      writer%dcdinfo%idcd = 0
      call set_status(status_value,status_ok,'Trajectory already closed')
   end subroutine close_trajectory
end submodule mdtraj_operations
