!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas in Austin
!
#include "tinker_macro.h"
submodule (restart) restart_codec
   use status,only:clear_status,set_status,status_invalid_argument,status_ok
   implicit none
contains

   module subroutine write_restart(filename,title,xbox,ybox,zbox,alpha,beta,&
      &gamma,x,y,z,v,a,aalt,status_value)
      character(*) ,intent(in) :: filename,title
      real(r_p)   ,intent(in) :: xbox,ybox,zbox
      real(t_p)   ,intent(in) :: alpha,beta,gamma
      real(r_p)   ,intent(in) :: x(:),y(:),z(:),v(:,:),a(:,:),aalt(:,:)
      type(status_t),intent(out) :: status_value
      character(len=240) :: record
      integer :: i,ios,unit,atom_count
      logical :: opened,ok

      call clear_status(status_value)
      opened = .false.
      ok = .false.
      atom_count = size(x)
      if (len_trim(filename) == 0 .or. atom_count <= 0 .or. &
         &size(y) /= atom_count .or. size(z) /= atom_count .or. &
         &size(v,1) /= 3 .or. size(v,2) /= atom_count .or. &
         &size(a,1) /= 3 .or. size(a,2) /= atom_count .or. &
         &size(aalt,1) /= 3 .or. size(aalt,2) /= atom_count) then
         call codec_error(status_value,'Invalid DYN state for writing')
         return
      end if

      open (newunit=unit,file=trim(filename),status='replace',action='write',&
         &form='formatted',iostat=ios,iomsg=record)
      if (ios /= 0) then
         call codec_error(status_value,'Unable to open DYN file for writing: '//trim(record))
         return
      end if
      opened = .true.

      write (unit,'(a)',iostat=ios,iomsg=record) ' Number of Atoms and Title :'
      if (ios /= 0) goto 900
      write (unit,'(i8,2x,a)',iostat=ios,iomsg=record) atom_count,trim(title)
      if (ios /= 0) goto 900
      write (unit,'(a)',iostat=ios,iomsg=record) ' Periodic Box Dimensions :'
      if (ios /= 0) goto 900
      write (unit,'(3d26.16)',iostat=ios,iomsg=record) xbox,ybox,zbox
      if (ios /= 0) goto 900
      write (unit,'(3d26.16)',iostat=ios,iomsg=record) alpha,beta,gamma
      if (ios /= 0) goto 900
      write (unit,'(a)',iostat=ios,iomsg=record) ' Current Atomic Positions :'
      if (ios /= 0) goto 900
      do i = 1,atom_count
         write (unit,'(3d26.16)',iostat=ios,iomsg=record) x(i),y(i),z(i)
         if (ios /= 0) goto 900
      end do
      write (unit,'(a)',iostat=ios,iomsg=record) ' Current Atomic Velocities :'
      if (ios /= 0) goto 900
      do i = 1,atom_count
         write (unit,'(3d26.16)',iostat=ios,iomsg=record) v(:,i)
         if (ios /= 0) goto 900
      end do
      write (unit,'(a)',iostat=ios,iomsg=record) ' Current Atomic Accelerations :'
      if (ios /= 0) goto 900
      do i = 1,atom_count
         write (unit,'(3d26.16)',iostat=ios,iomsg=record) a(:,i)
         if (ios /= 0) goto 900
      end do
      write (unit,'(a)',iostat=ios,iomsg=record) ' Alternate Atomic Accelerations :'
      if (ios /= 0) goto 900
      do i = 1,atom_count
         write (unit,'(3d26.16)',iostat=ios,iomsg=record) aalt(:,i)
         if (ios /= 0) goto 900
      end do
      ok = .true.

900   if (opened) then
         close (unit=unit,iostat=ios,iomsg=record)
         if (ok .and. ios /= 0) then
            call codec_error(status_value,'Unable to close DYN file after writing: '//trim(record))
            ok = .false.
         end if
      end if
      if (ok) then
         call set_status(status_value,status_ok,'DYN restart written without aalt2')
      else if (status_value%code == status_ok) then
         call codec_error(status_value,'Unable to write DYN restart: '//trim(record))
      end if
   end subroutine write_restart

   module subroutine read_restart(filename,expected_atoms,xbox,ybox,zbox,alpha,&
      &beta,gamma,x,y,z,v,a,aalt,aalt2,aalt2_present,status_value)
      character(*) ,intent(in) :: filename
      integer      ,intent(in) :: expected_atoms
      real(r_p)    ,intent(out) :: xbox,ybox,zbox
      real(t_p)    ,intent(out) :: alpha,beta,gamma
      real(r_p)    ,intent(inout) :: x(:),y(:),z(:),v(:,:),a(:,:),aalt(:,:),&
         &aalt2(:,:)
      logical      ,intent(out) :: aalt2_present
      type(status_t),intent(out) :: status_value
      real(r_p),allocatable :: x_read(:),y_read(:),z_read(:),v_read(:,:),&
         &a_read(:,:),aalt_read(:,:),aalt2_read(:,:)
      real(r_p) :: xbox_read,ybox_read,zbox_read
      real(r_p) :: angle_values(3)
      real(t_p) :: alpha_read,beta_read,gamma_read
      character(len=240) :: record
      character(len=240) :: message
      character(len=240) :: close_message
      integer :: i,ios,unit,atom_count,close_ios
      logical :: opened,ok

      call clear_status(status_value)
      aalt2_present = .false.
      message = ''
      opened = .false.
      ok = .false.
      if (len_trim(filename) == 0 .or. expected_atoms <= 0 .or. &
         &size(x) /= expected_atoms .or. size(y) /= expected_atoms .or. &
         &size(z) /= expected_atoms .or. size(v,1) /= 3 .or. &
         &size(v,2) /= expected_atoms .or. size(a,1) /= 3 .or. &
         &size(a,2) /= expected_atoms .or. size(aalt,1) /= 3 .or. &
         &size(aalt,2) /= expected_atoms .or. size(aalt2,1) /= 3 .or. &
         &size(aalt2,2) /= expected_atoms) then
         call codec_error(status_value,'Invalid destination for DYN restart')
         return
      end if

      open (newunit=unit,file=trim(filename),status='old',action='read',&
         &form='formatted',iostat=ios,iomsg=record)
      if (ios /= 0) then
         call codec_error(status_value,'Unable to open DYN file for reading: '//trim(record))
         return
      end if
      opened = .true.

      read (unit,'(a240)',iostat=ios,iomsg=message)
      if (ios /= 0) goto 900
      read (unit,'(i8)',iostat=ios,iomsg=message) atom_count
      if (ios /= 0) goto 900
      if (atom_count /= expected_atoms) then
         call codec_error(status_value,'DYN restart has an incorrect atom count')
         goto 900
      end if

      read (unit,'(a240)',iostat=ios,iomsg=message)
      if (ios /= 0) goto 900
      read (unit,'(3d26.16)',iostat=ios,iomsg=message) xbox_read,ybox_read,zbox_read
      if (ios /= 0) goto 900
      read (unit,'(3d26.16)',iostat=ios,iomsg=message) angle_values(1),&
         &angle_values(2),angle_values(3)
      if (ios /= 0) goto 900
      alpha_read = real(angle_values(1),t_p)
      beta_read = real(angle_values(2),t_p)
      gamma_read = real(angle_values(3),t_p)

      allocate(x_read(expected_atoms),y_read(expected_atoms),z_read(expected_atoms))
      allocate(v_read(3,expected_atoms),a_read(3,expected_atoms),aalt_read(3,expected_atoms))
      read (unit,'(a240)',iostat=ios,iomsg=message)
      if (ios /= 0) goto 900
      do i = 1,expected_atoms
         read (unit,'(3d26.16)',iostat=ios,iomsg=message) x_read(i),y_read(i),z_read(i)
         if (ios /= 0) goto 900
      end do
      read (unit,'(a240)',iostat=ios,iomsg=message)
      if (ios /= 0) goto 900
      do i = 1,expected_atoms
         read (unit,'(3d26.16)',iostat=ios,iomsg=message) v_read(:,i)
         if (ios /= 0) goto 900
      end do
      read (unit,'(a240)',iostat=ios,iomsg=message)
      if (ios /= 0) goto 900
      do i = 1,expected_atoms
         read (unit,'(3d26.16)',iostat=ios,iomsg=message) a_read(:,i)
         if (ios /= 0) goto 900
      end do
      read (unit,'(a240)',iostat=ios,iomsg=message)
      if (ios /= 0) goto 900
      do i = 1,expected_atoms
         read (unit,'(3d26.16)',iostat=ios,iomsg=message) aalt_read(:,i)
         if (ios /= 0) goto 900
      end do

      read (unit,'(a240)',iostat=ios,iomsg=message)
      if (ios > 0) goto 900
      if (ios == 0) then
         allocate(aalt2_read(3,expected_atoms))
         do i = 1,expected_atoms
            read (unit,'(3d26.16)',iostat=ios,iomsg=message) aalt2_read(:,i)
            if (ios /= 0) goto 900
         end do
         aalt2_present = .true.
      end if
      ok = .true.

900   if (opened) then
         close (unit=unit,iostat=close_ios,iomsg=close_message)
         if (ok .and. close_ios /= 0) then
            call codec_error(status_value,'Unable to close DYN file after reading: '//trim(close_message))
            ok = .false.
         end if
      end if
      if (ok) then
         x = x_read
         y = y_read
         z = z_read
         v = v_read
         a = a_read
         aalt = aalt_read
         if (aalt2_present) aalt2 = aalt2_read
         xbox = xbox_read
         ybox = ybox_read
         zbox = zbox_read
         alpha = alpha_read
         beta = beta_read
         gamma = gamma_read
         call set_status(status_value,status_ok,'DYN restart read')
      else if (status_value%code == status_ok) then
         call codec_error(status_value,'Malformed DYN restart: '//trim(message))
      end if
      if (allocated(x_read)) deallocate(x_read,y_read,z_read,v_read,a_read,aalt_read)
      if (allocated(aalt2_read)) deallocate(aalt2_read)
   end subroutine read_restart

   subroutine codec_error(status_value,message)
      type(status_t),intent(out) :: status_value
      character(*),intent(in) :: message

      call set_status(status_value,status_invalid_argument,trim(message))
   end subroutine codec_error

end submodule restart_codec
