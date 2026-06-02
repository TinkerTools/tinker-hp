!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine kangang  --  angle-angle parameter assignment  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "kangang" assigns the parameters for angle-angle cross term
!     interactions and processes new or changed parameter values
!
!
!> @brief 
!> assigns the parameters for angle-angle cross term
!> interactions and processes new or changed parameter values
!> @param no params
subroutine kangang
   use angang
   use angle
   use atmlst
   use atmtyp
   use atoms
   use couple
   use domdec
   use inform
   use iounit
   use kanang
   use keys
   use potent
   use tors
   implicit none
   integer i,j,k,m,next
   integer it,ia,ic
   integer nang,jang,kang
   real*8 fa,faa,aak(3)
   logical header
   character*20 keyword
   character*240 record
   character*240 string
!
   if (deb_Path) write(iout,*), 'kangang '
!
!
!     process keywords containing angle-angle parameters
!
   header = .true.
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      if (keyword(1:7) .eq. 'ANGANG ') then
         it = 0
         do j = 1, 3
            aak(j) = 0.0d0
         end do
         string = record(next:240)
         read (string,*,err=10,end=10)  it,(aak(j),j=1,3)
10       continue
         if (.not. silent) then
            if (header) then
               header = .false.
               if (rank.eq.0) write (iout,20)
20             format (/,' Additional Angle-Angle Parameters :',&
               &//,5x,'Atom Class',8x,'K(AA) 1',4x,&
               &'K(AA) 2',4x,'K(AA) 3',/)
            end if
            if (rank.eq.0) write (iout,30)  it,(aak(j),j=1,3)
30          format (9x,i3,7x,3f11.3)
         end if
         do j = 1, 3
            anan(j,it) = aak(j)
         end do
      end if
   end do
!
!     assign the angle-angle parameters for each angle pair
!
   nangang = 0
   do i = 1, n
      nang = n12(i) * (n12(i)-1) / 2
      it = class(i)
      nbangang(i) = nangang
      do j = 1, nang-1
         jang = anglist(j,i)
         ia = iang(1,jang)
         ic = iang(3,jang)
         m = 1
         if (atomic(ia) .le. 1)  m = m + 1
         if (atomic(ic) .le. 1)  m = m + 1
         fa = anan(m,it)
         do k = j+1, nang
            kang = anglist(k,i)
            ia = iang(1,kang)
            ic = iang(3,kang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            faa = fa * anan(m,it)
            if (faa .ne. 0.0d0) then
               nangang = nangang + 1
            end if
         end do
      end do
   end do
!
!     deallocate global pointers if necessary
!
   call dealloc_shared_angang
!
!     allocate global pointers
!
   call alloc_shared_angang
!
   nangangloc = 0
   do i = 1, n
      nang = n12(i) * (n12(i)-1) / 2
      it = class(i)
      do j = 1, nang-1
         jang = anglist(j,i)
         ia = iang(1,jang)
         ic = iang(3,jang)
         m = 1
         if (atomic(ia) .le. 1)  m = m + 1
         if (atomic(ic) .le. 1)  m = m + 1
         fa = anan(m,it)
         do k = j+1, nang
            kang = anglist(k,i)
            ia = iang(1,kang)
            ic = iang(3,kang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            faa = fa * anan(m,it)
            if (faa .ne. 0.0d0) then
               nangangloc = nangangloc + 1
               iaa(1,nangangloc) = jang
               iaa(2,nangangloc) = kang
               kaa(nangangloc) = faa
            end if
         end do
      end do
   end do
!
!     turn off the angle-angle potential if it is not used
!
   if (nangang .eq. 0)  use_angang = .false.
   return
end
!
!     subroutine kangang_update: update local angang
!
!> @brief 
!> update local angang parameters
!> @param no params
subroutine kangang_update
   use angang
   use angle
   use atmlst
   use atmtyp
   use atoms
   use couple
   use domdec
   use inform
   use iounit
   use kanang
   use keys
   use potent
   use tors
   implicit none
   integer i,j,k,m
   integer it,ia,ic
   integer iglob
   integer nang,jang,kang,nangangloc1,angangcount
   real*8 fa,faa
!
   if (deb_Path) write(iout,*), 'kangang_update '
!
!
   if (allocated(angangglob)) deallocate(angangglob)
   allocate (angangglob(15*nloc))
   nangangloc = 0
   do i = 1, nloc
      iglob = glob(i)
      nang = n12(iglob) * (n12(iglob)-1) / 2
      it = class(iglob)
      angangcount = nbangang(iglob)
      nangangloc1 = 0
      do j = 1, nang-1
         jang = anglist(j,iglob)
         ia = iang(1,jang)
         ic = iang(3,jang)
         m = 1
         if (atomic(ia) .le. 1)  m = m + 1
         if (atomic(ic) .le. 1)  m = m + 1
         fa = anan(m,it)
         do k = j+1, nang
            kang = anglist(k,iglob)
            ia = iang(1,kang)
            ic = iang(3,kang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            faa = fa * anan(m,it)
            if (faa .ne. 0.0d0) then
               nangangloc = nangangloc + 1
               nangangloc1 = nangangloc1 + 1
               angangglob(nangangloc) = angangcount + nangangloc1
            end if
         end do
      end do
   end do
!
   return
end
!
!     subroutine dealloc_shared_angang : deallocate shared memory pointers for angang
!     parameter arrays
!
!> @brief 
!> deallocate shared memory pointers for angang
!> parameter arrays
!> @param no params
subroutine dealloc_shared_angang
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use angang
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
!
   if (associated(iaa)) then
      CALL MPI_Win_shared_query(winiaa, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winiaa,ierr)
   end if
   if (associated(kaa)) then
      CALL MPI_Win_shared_query(winkaa, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winkaa,ierr)
   end if
   return
end
!
!     subroutine alloc_shared_angang : allocate shared memory pointers for angang
!     parameter arrays
!
!> @brief 
!> allocate shared memory pointers for angang
!> parameter arrays
!> @param no params
subroutine alloc_shared_angang
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use sizes
   use angang
   use domdec
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   integer :: arrayshape(1),arrayshape2(2)
!
!     iaa
!
   arrayshape2=(/2,nangang/)
   if (hostrank == 0) then
      windowsize = int(2*nangang,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winiaa, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winiaa, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,iaa,arrayshape2)
!
!     kaa
!
   arrayshape=(/nangang/)
   if (hostrank == 0) then
      windowsize = int(nangang,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winkaa, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winkaa, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,kaa,arrayshape)
   return
end
