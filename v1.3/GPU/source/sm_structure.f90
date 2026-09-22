!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
submodule (structure) structure_loader
   use atmtyp, only: tag, name
   use atoms, only: n, nloop, n_pe, pbcWrapIdx, pbcWrap, x, y, z, xold, yold, zold, type
   use atomsMirror, only: xm => x, ym => y, zm => z, atomsmirror_init, &
      &reCast_position, download_position
   use bound, only: use_bounds
   use boxes
   use couple, only: n12, i12
   use iso_c_binding, only: c_f_pointer, c_loc
   use inform, only: app_id, dynamic_a, pimd_a, abort
   use nvshmem, only: npes
   use output, only: coordtype
   use sizes, only: maxvalue, maxbnd, maxang, maxtors, maxbitor
   use status, only: clear_status, failed, set_status, status_invalid_argument, &
      &status_missing_argument
   use tinMemory, only: mhostonly, prmem_request
   use tinheader, only: ti_p, re_p
   use titles, only: title, ltitle
   implicit none

contains

   module subroutine load_initial(configuration, structure_value, status_value)
      implicit none
      type(config_t), intent(in) :: configuration
      type(structure_t), intent(out) :: structure_value
      type(status_t), intent(out) :: status_value
      character(len=240) :: xyzfile
      integer :: unit, ios
      logical :: exists

      call clear_status(status_value)
      call clear_structure(structure_value)
      if (len_trim(configuration%xyz_path) == 0) then
         call set_status(status_value, status_missing_argument, &
            'Configuration supplied no coordinate file')
         return
      end if

      xyzfile = configuration%xyz_path
      call suffix(xyzfile, 'xyz', 'old')
      inquire(file=trim(xyzfile), exist=exists)
      if (.not.exists) then
         call set_status(status_value, status_invalid_argument, &
            'Unable to Find the Cartesian Coordinates File')
         return
      end if
      open(newunit=unit,file=trim(xyzfile),status='old',action='read',iostat=ios)
      if (ios /= 0) then
         call set_status(status_value, status_invalid_argument, &
            'Unable to Open the Cartesian Coordinates File')
         return
      end if

      structure_value%source_path = xyzfile
      call read_xyz(unit, structure_value, status_value)
      close(unit)
      if (failed(status_value)) then
         call clear_structure(structure_value)
         return
      end if

      call publish_legacy_state(structure_value)
   end subroutine load_initial

   module subroutine read_frame(unit, structure_value, status_value)
      implicit none
      integer, intent(in) :: unit
      type(structure_t), intent(inout) :: structure_value
      type(status_t), intent(out) :: status_value
      type(structure_t) :: next_structure
      logical :: at_end

      call read_xyz(unit, next_structure, status_value, at_end)
      if (at_end) then
         abort = .true.
         return
      end if
      if (failed(status_value)) return
      if (next_structure%atom_count /= structure_value%atom_count) then
         call set_status(status_value, status_invalid_argument, &
            'Trajectory Frame Atom Count Does Not Match the Initial Structure')
         return
      end if

      next_structure%source_path = structure_value%source_path
      call publish_legacy_state(next_structure)
      structure_value = next_structure
   end subroutine read_frame

   subroutine read_xyz(unit, structure_value, status_value, at_end)
      implicit none
      integer, intent(in) :: unit
      type(structure_t), intent(inout) :: structure_value
      type(status_t), intent(out) :: status_value
      logical, intent(out), optional :: at_end
      character(len=240) :: record
      character(len=3) :: atom_name
      integer :: i, j, k, ios, atom_count, atom_tag, atom_type
      integer :: neighbor(maxvalue)
      integer, allocatable :: remap(:)
      logical :: has_box
      real(re_p) :: xpos, ypos, zpos
      real(re_p) :: box_length(3)
      real(ti_p) :: box_angle(3)

      call clear_status(status_value)
      if (present(at_end)) at_end = .false.
      call next_record(unit, record, ios)
      if (ios /= 0) then
         if (present(at_end)) then
            at_end = .true.
            return
         end if
         call set_status(status_value, status_invalid_argument, &
            'Coordinate File Does Not Contain an Atom Count')
         return
      end if
      read(record,*,iostat=ios) atom_count
      if (ios /= 0 .or. atom_count <= 0) then
         call set_status(status_value, status_invalid_argument, &
            'Coordinate File Contains an Invalid Atom Count')
         return
      end if

      structure_value%atom_count = atom_count
      call extract_title(record, structure_value%title)
      allocate(structure_value%x(atom_count), structure_value%y(atom_count), &
         &structure_value%z(atom_count), structure_value%tag(atom_count), &
         &structure_value%atom_type(atom_count), structure_value%atom_name(atom_count), &
         &structure_value%neighbor_count(atom_count), &
         &structure_value%connectivity(maxvalue,atom_count))
      structure_value%x = 0.0_re_p
      structure_value%y = 0.0_re_p
      structure_value%z = 0.0_re_p
      structure_value%tag = 0
      structure_value%atom_type = 0
      structure_value%atom_name = ' '
      structure_value%neighbor_count = 0
      structure_value%connectivity = 0
      structure_value%box_length = 0.0_re_p
      structure_value%box_angle = 0.0_ti_p

      has_box = .false.
      do i = 1, atom_count
         do
            call next_record(unit, record, ios)
            if (ios /= 0) then
               call set_status(status_value, status_invalid_argument, &
                  'Coordinate File Ends Before All Atoms Are Read')
               return
            end if
            neighbor = 0
            read(record,*,iostat=ios) atom_tag, atom_name, xpos, ypos, zpos, &
               &atom_type, (neighbor(j),j=1,maxvalue)
            if (ios <= 0 .and. verify(trim(atom_name), &
                  &'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789') == 0) exit
            if (i /= 1 .or. has_box) then
               call set_status(status_value, status_invalid_argument, &
                  'Malformed Atom Record in Coordinate File')
               return
            end if
            read(record,*,iostat=ios) box_length, box_angle
            if (ios /= 0) then
               call set_status(status_value, status_invalid_argument, &
                  'Malformed Periodic Box Record in Coordinate File')
               return
            end if
            structure_value%box_length = box_length
            structure_value%box_angle = box_angle
            has_box = .true.
         end do
         structure_value%tag(i) = atom_tag
         structure_value%atom_name(i) = atom_name
         structure_value%x(i) = xpos
         structure_value%y(i) = ypos
         structure_value%z(i) = zpos
         structure_value%atom_type(i) = atom_type
         do j = 1, maxvalue
            if (neighbor(j) == 0) exit
            if (neighbor(j) < 0) then
               call set_status(status_value, status_invalid_argument, &
                  'Negative Atom Connectivity in Coordinate File')
               return
            end if
            structure_value%neighbor_count(i) = j
            structure_value%connectivity(j,i) = neighbor(j)
         end do
      end do
      structure_value%has_bounds = has_box

      do i = 1, atom_count
         if (structure_value%tag(i) < 1) then
            call set_status(status_value, status_invalid_argument, &
               'Atom Labels Must Be Positive and Unique')
            return
         end if
      end do
      allocate(remap(maxval(structure_value%tag)))
      remap = 0
      do i = 1, atom_count
         if (remap(structure_value%tag(i)) /= 0) then
            call set_status(status_value, status_invalid_argument, &
               'Atom Labels Must Be Positive and Unique')
            return
         end if
         remap(structure_value%tag(i)) = i
      end do
      do i = 1, atom_count
         do j = 1, structure_value%neighbor_count(i)
            k = structure_value%connectivity(j,i)
            if (k < 1 .or. k > size(remap)) then
               call set_status(status_value, status_invalid_argument, &
                  'Atom Connectivity References an Unknown Atom')
               return
            end if
            if (remap(k) == 0) then
               call set_status(status_value, status_invalid_argument, &
                  'Atom Connectivity References an Unknown Atom')
               return
            end if
            structure_value%connectivity(j,i) = remap(k)
         end do
         call sort_neighbors(structure_value%neighbor_count(i), &
            &structure_value%connectivity(:,i))
      end do
      do i = 1, atom_count
         do j = 1, structure_value%neighbor_count(i)
            k = structure_value%connectivity(j,i)
            if (.not.any(structure_value%connectivity(1:structure_value%neighbor_count(k),k) == i)) then
               call set_status(status_value, status_invalid_argument, &
                  'Coordinate Connectivity Must Be Bidirectional')
               return
            end if
         end do
      end do
      structure_value%max_neighbor_count = maxval(structure_value%neighbor_count)
   end subroutine read_xyz

   subroutine publish_legacy_state(structure_value)
      implicit none
      type(structure_t), intent(in) :: structure_value
      integer :: i
      logical :: use_wrap

      n = structure_value%atom_count
      call prmem_request(x,n)
      call prmem_request(y,n)
      call prmem_request(z,n)
      call prmem_request(xold,n)
      call prmem_request(yold,n)
      call prmem_request(zold,n)
      call prmem_request(type,n)
      call prmem_request(tag,n,config=mhostonly)
      call atomsmirror_init

      if (associated(name)) deallocate(name)
      allocate(name(n))
      do i = 1, n
         name(i) = structure_value%atom_name(i)
      end do
      type = structure_value%atom_type
      do i = 1, n
         tag(i) = i
      end do
      xold = 0.0_ti_p
      yold = 0.0_ti_p
      zold = 0.0_ti_p
      do i = 1, n
         xm(i) = structure_value%x(i)
         ym(i) = structure_value%y(i)
         zm(i) = structure_value%z(i)
      end do
!$acc update device(xm,ym,zm,xold,yold,zold,type)
      call reCast_position
!$acc wait
      call download_position

!$acc wait
      if (allocated(n12)) then
!$acc exit data delete(n12)
         deallocate(n12)
      end if
      if (allocated(i12)) then
!$acc exit data delete(i12)
         deallocate(i12)
      end if
      allocate(n12(n), i12(maxvalue,n))
      n12 = structure_value%neighbor_count
      i12 = structure_value%connectivity
!$acc enter data create(n12,i12)
!$acc update device(n12(:),i12(:,:))

      nloop = merge(n,(int(n/16)+1)*16,(mod(n,16) == 0))
      if (npes > 0) n_pe = ceiling(1.0*n/npes)
      maxbnd = 4*n
      maxang = 6*n
      maxtors = 18*n
      maxbitor = 54*n

      use_wrap = app_id == dynamic_a .or. app_id == pimd_a
      if (use_wrap) then
         call prmem_request(pbcWrap,n)
         call c_f_pointer(c_loc(pbcWrap),pbcWrapIdx,[4*n])
         pbcWrap = 0
!$acc update device(pbcWrap)
      end if

      use_bounds = structure_value%has_bounds
      if (structure_value%has_bounds) then
         xbox = structure_value%box_length(1)
         ybox = structure_value%box_length(2)
         zbox = structure_value%box_length(3)
         alpha = structure_value%box_angle(1)
         beta = structure_value%box_angle(2)
         gamma = structure_value%box_angle(3)
!$acc update device(xbox,ybox,zbox,alpha,beta,gamma)
      end if
      coordtype = 'CARTESIAN'
      title = structure_value%title
      ltitle = len_trim(title)
      abort = .false.
   end subroutine publish_legacy_state

   subroutine next_record(unit, record, ios)
      implicit none
      integer, intent(in) :: unit
      character(*), intent(out) :: record
      integer, intent(out) :: ios

      do
         read(unit,'(a)',iostat=ios) record
         if (ios /= 0 .or. len_trim(record) /= 0) return
      end do
   end subroutine next_record

   subroutine extract_title(record, value)
      implicit none
      character(*), intent(in) :: record
      character(*), intent(out) :: value
      integer :: first, last

      first = 1
      do while (first <= len(record) .and. (record(first:first) == ' ' .or. &
         &record(first:first) == achar(9)))
         first = first + 1
      end do
      do while (first <= len(record) .and. record(first:first) /= ' ' .and. &
         &record(first:first) /= achar(9))
         first = first + 1
      end do
      if (first > len(record)) then
         value = ''
         return
      end if
      last = len_trim(record)
      value = adjustl(record(first:last))
   end subroutine extract_title

   subroutine sort_neighbors(count, values)
      implicit none
      integer, intent(in) :: count
      integer, intent(inout) :: values(:)
      integer :: i, j, value

      do i = 1, count - 1
         do j = i + 1, count
            if (values(j) < values(i)) then
               value = values(i)
               values(i) = values(j)
               values(j) = value
            end if
         end do
      end do
   end subroutine sort_neighbors

   subroutine clear_structure(structure_value)
      implicit none
      type(structure_t), intent(out) :: structure_value

      structure_value%source_path = ''
      structure_value%title = ''
      structure_value%atom_count = 0
      structure_value%max_neighbor_count = 0
      structure_value%has_bounds = .false.
      structure_value%box_length = 0.0_re_p
      structure_value%box_angle = 0.0_ti_p
   end subroutine clear_structure

end submodule structure_loader
