!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###########################################################
!     ##                                                       ##
!     ##  subroutine active  --  set the list of active atoms  ##
!     ##                                                       ##
!     ###########################################################
!
!
!     "active" sets the list of atoms that are used during
!     each potential energy function calculation
!
!
!> @brief 
!> sets the list of atoms that are used during each potential energy function calculation
!> @param no params
subroutine active
   use atoms
   use domdec
   use keys
   use inform
   use iounit
   use usage
   use mpi
   implicit none
   integer i,j,next,ierr
   integer nmobile,nfixed
   integer center,nsphere
   integer, allocatable :: mobile(:)
   integer, allocatable :: fixed(:)
   real*8 xcenter,ycenter,zcenter
   real*8 radius,radius2,dist2
   character*20 keyword
   character*240 record
   character*240 string
!
   if (deb_Path) write(iout,*), 'active '
!
!     deallocate global pointers if necessary
!
   call dealloc_shared_active
!
!     allocate global pointers
!
   call alloc_shared_active
!
   if (hostrank.ne.0) goto 90
!
!     perform dynamic allocation of some local arrays
!
   allocate (mobile(n))
   allocate (fixed(n))
!
!     set defaults for the numbers and lists of active atoms
!
   nuse = n
   use(0) = .false.
   do i = 1, n
      use(i) = .true.
   end do
   nmobile = 0
   nfixed = 0
   do i = 1, n
      mobile(i) = 0
      fixed(i) = 0
   end do
   nsphere = 0
!
!     get any keywords containing active atom parameters
!
   do j = 1, nkey
      next = 1
      record = keyline(j)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:240)
!
!     get any lists of atoms whose coordinates are active
!
      if (keyword(1:7) .eq. 'ACTIVE ') then
         read (string,*,err=10,end=10)  (mobile(i),i=nmobile+1,n)
10       continue
         do while (mobile(nmobile+1) .ne. 0)
            nmobile = nmobile + 1
            mobile(nmobile) = max(-n,min(n,mobile(nmobile)))
         end do
!
!     get any lists of atoms whose coordinates are inactive
!
      else if (keyword(1:9) .eq. 'INACTIVE ') then
         read (string,*,err=20,end=20)  (fixed(i),i=nfixed+1,n)
20       continue
         do while (fixed(nfixed+1) .ne. 0)
            nfixed = nfixed + 1
            fixed(nfixed) = max(-n,min(n,fixed(nfixed)))
         end do
!
!     get the center and radius of the sphere of active atoms
!
      else if (keyword(1:7) .eq. 'SPHERE ') then
         center = 0
         xcenter = 0.0d0
         ycenter = 0.0d0
         zcenter = 0.0d0
         radius = 0.0d0
         read (string,*,err=30,end=30)  xcenter,ycenter,&
         &zcenter,radius
30       continue
         if (radius .eq. 0.0d0) then
            read (string,*,err=60,end=60)  center,radius
            xcenter = x(center)
            ycenter = y(center)
            zcenter = z(center)
         end if
         nsphere = nsphere + 1
         if (nsphere .eq. 1) then
            nuse = 0
            do i = 1, n
               use(i) = .false.
            end do
            if (verbose) then
               write (iout,40)
40             format (/,' Active Site Spheres used to',&
               &' Select Active Atoms :',&
               &//,3x,'Atom Center',11x,'Coordinates',&
               &12x,'Radius',6x,'# Active Atoms')
            end if
         end if
         radius2 = radius * radius
         do i = 1, n
            if (.not. use(i)) then
               dist2 = (x(i)-xcenter)**2 + (y(i)-ycenter)**2&
               &+ (z(i)-zcenter)**2
               if (dist2 .le. radius2) then
                  nuse = nuse + 1
                  use(i) = .true.
               end if
            end if
         end do
         if (verbose.and.(rank.eq.0)) then
            write (iout,50)  center,xcenter,ycenter,&
            &zcenter,radius,nuse
50          format (2x,i8,6x,3f9.2,2x,f9.2,7x,i8)
         end if
60       continue
      end if
   end do
!
!     set active atoms to those marked as not inactive
!
   i = 1
   do while (fixed(i) .ne. 0)
      if (fixed(i) .gt. 0) then
         j = fixed(i)
         if (use(j)) then
            use(fixed(i)) = .false.
            nuse = nuse - 1
         end if
         i = i + 1
      else
         do j = abs(fixed(i)), abs(fixed(i+1))
            if (use(j)) then
               use(j) = .false.
               nuse = nuse - 1
            end if
         end do
         i = i + 2
      end if
   end do
!
!     set active atoms to only those marked as active
!
   i = 1
   do while (mobile(i) .ne. 0)
      if (i .eq. 1) then
         nuse = 0
         do j = 1, n
            use(j) = .false.
         end do
      end if
      if (mobile(i) .gt. 0) then
         j = mobile(i)
         if (.not. use(j)) then
            use(j) = .true.
            nuse = nuse + 1
         end if
         i = i + 1
      else
         do j = abs(mobile(i)), abs(mobile(i+1))
            if (.not. use(j)) then
               use(j) = .true.
               nuse = nuse + 1
            end if
         end do
         i = i + 2
      end if
   end do
!
!     use logical array to set the list of active atoms
!
   j = 0
   do i = 1, n
      if (use(i)) then
         j = j + 1
         iuse(j) = i
      end if
   end do
!
!     output the final list of the active atoms
!
   if ((debug .and. nuse.gt.0 .and. nuse.lt.n).and.(rank.eq.0)) then
      write (iout,70)
70    format (/,' List of Active Atoms for Energy',&
      &' Calculations :',/)
      write (iout,80)  (iuse(i),i=1,nuse)
80    format (3x,10i7)
   end if
!
!     perform deallocation of some local arrays
!
   deallocate (mobile)
   deallocate (fixed)
90 continue
   call MPI_BCAST(nuse,1,MPI_INT,0,hostcomm,ierr)
   return
end
!
!     subroutine dealloc_shared_active : deallocate shared memory pointers for active
!     parameter arrays
!
!> @brief 
!> deallocate shared memory pointers for active parameter arrays
!> @param no params
subroutine dealloc_shared_active
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use usage
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   if (associated(use)) then
      CALL MPI_Win_shared_query(winuse, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winuse,ierr)
   end if
   if (associated(iuse)) then
      CALL MPI_Win_shared_query(winiuse, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winiuse,ierr)
   end if
   return
end
!
!
!     subroutine alloc_shared_active : allocate shared memory pointers for active
!     parameter arrays
!
!> @brief 
!> allocate shared memory pointers for active parameter arrays
!> @param no params
subroutine alloc_shared_active
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use atoms
   use domdec
   use usage
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   integer :: arrayshape(1)
!
!     use
!
   arrayshape=(/n+1/)
   if (hostrank == 0) then
      windowsize = int(n+1,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winuse, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winuse, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,use,arrayshape)
   use(0:) => use
!
!     iuse
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winiuse, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winiuse, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,iuse,arrayshape)
   return
end
