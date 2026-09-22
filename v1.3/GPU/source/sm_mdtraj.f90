!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
submodule (mdtraj) mdtraj_operations
   use atmtyp
   use atoms     ,only: type,pbcunwrap,pbcWrapIdx
   use bound
   use boxes
   use cell      ,only: xcell,ycell,zcell
   use couple
   use dcdmod ,only: dcdio_close,dcdio_write
   use inform ,only: digits
   use status ,only: clear_status,failed,set_status,status_invalid_argument,status_ok
   use tinheader
   use titles
   use timestat, only: timer_enter,timer_exit,timer_io
   implicit none
   integer,external :: freeunit
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
      if (.not.writer%is_active .or. .not.valid_frame(frame) .or. &
         &step < 0 .or. dt <= 0.0) then
         call set_status(status_value,status_invalid_argument,&
            &'Trajectory frame is not a complete host frame')
         return
      end if

      if (writer%backend == dcd_backend) then
         call dcdio_write(writer%dcdinfo,writer%filename,step,dt,&
            &writer%append_mode,frame%x,frame%y,frame%z,status_value)
         writer%owns_dcd_file = writer%dcdinfo%idcd /= 0
         if (failed(status_value)) return
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
         call write_xyz(unit,frame)
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

   module subroutine write_snapshot(filename,frame,status_value,append)
      character(*),intent(in) :: filename
      type(frame_view_t),intent(in) :: frame
      type(status_t),intent(out) :: status_value
      logical,intent(in),optional :: append
      integer :: unit,ios
      logical :: append_mode,exist

      call clear_status(status_value)
      if (len_trim(filename) == 0 .or. .not.valid_frame(frame)) then
         call set_status(status_value,status_invalid_argument,&
            &'Snapshot needs a name and a complete host frame')
         return
      end if

      unit = freeunit ()
      append_mode = .false.
      if (present(append)) append_mode = append
      inquire(file=trim(filename),exist=exist)
      if (append_mode .and. exist) then
         open(unit=unit,file=trim(filename),status='old',position='append',&
            &action='write',iostat=ios)
      else
         open(unit=unit,file=trim(filename),status='replace',action='write',&
            &iostat=ios)
      end if
      if (ios /= 0) then
         call set_status(status_value,status_invalid_argument,&
            &'Unable to open XYZ snapshot')
         return
      end if
      call write_xyz(unit,frame)
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
      if (writer%owns_dcd_file) call dcdio_close(writer%dcdinfo,status_value)
      if (failed(status_value)) return
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

   subroutine write_xyz(unit,frame)
      type(frame_view_t),intent(in) :: frame
      integer,intent(in) :: unit
      integer :: i,k,atom_count,size,crdsiz
      real(r_p) :: crdmin,crdmax,xi,yi,zi
      character*2 :: atmc,crdc,digc
      character*25 :: fstr

      atom_count = frame%atoms
      call timer_enter(timer_io)
      atmc = 'i6'
      if (atom_count .ge. 100000) atmc = 'i7'
      if (atom_count .ge. 1000000) atmc = 'i8'
      crdmin = 0.0_ti_p
      crdmax = 0.0_ti_p
      do i = 1, atom_count
         crdmin = min(crdmin,frame%x(i),frame%y(i),frame%z(i))
         crdmax = max(crdmax,frame%x(i),frame%y(i),frame%z(i))
      end do
      crdsiz = 6
      if (crdmin .le. -1000.0_re_p) crdsiz = 7
      if (crdmax .ge. 10000.0_re_p) crdsiz = 7
      if (crdmin .le. -10000.0_re_p) crdsiz = 8
      if (crdmax .ge. 100000.0_re_p) crdsiz = 8
      crdsiz = crdsiz + max(6,digits)
      size = 0
      call numeral(crdsiz,crdc,size)
      if (digits .le. 6) then
         digc = '6 '
      else if (digits .le. 8) then
         digc = '8'
      else
         digc = '10'
      end if

      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (unit,fstr(1:4)) atom_count
      else
         fstr = '('//atmc//',2x,a)'
         write (unit,fstr(1:9)) atom_count,title(1:ltitle)
      end if
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (unit,fstr) xbox,ybox,zbox,alpha,beta,gamma
      end if

      fstr = '('//atmc//',2x,a3,3f'//crdc//&
         &'.'//digc//',i6,8'//atmc//')'
      if (pbcunwrap) then
         do i = 1, atom_count
            xi = frame%x(i)
            yi = frame%y(i)
            zi = frame%z(i)
            xi = xi + int(pbcWrapIdx(4*(i-1)+1))*xcell
            yi = yi + int(pbcWrapIdx(4*(i-1)+2))*ycell
            zi = zi + int(pbcWrapIdx(4*(i-1)+3))*zcell
            write (unit,fstr) i,name(i),xi,yi,zi,type(i),&
               &(i12(k,i),k=1,n12(i))
         end do
      else
         do i = 1, atom_count
            write (unit,fstr) i,name(i),frame%x(i),frame%y(i),frame%z(i),&
               &type(i),(i12(k,i),k=1,n12(i))
         end do
      end if
      call timer_exit(timer_io)
   end subroutine write_xyz

   logical function valid_frame(frame)
      type(frame_view_t),intent(in) :: frame

      valid_frame = frame%atoms > 0 .and. frame%host_global .and. &
         &associated(frame%x) .and. associated(frame%y) .and. &
         &associated(frame%z) .and. size(frame%x) == frame%atoms .and. &
         &size(frame%y) == frame%atoms .and. size(frame%z) == frame%atoms
   end function valid_frame
end submodule mdtraj_operations
