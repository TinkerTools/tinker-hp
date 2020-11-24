c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine attach  --  setup of connectivity arrays  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "attach" generates lists of 1-3, 1-4 and 1-5 connectivities
c     starting from the previously determined list of attached
c     atoms (ie, 1-2 connectivity)
c
c
      subroutine attach
      use sizes
      use atoms
      use couple
      use domdec
      use iounit
      implicit none
      integer i,j,k,m
      integer jj,kk
c
c     allocate global arrays
c
      call alloc_shared_attach
c
c       only master of the node fill the arrays
c
      if (hostrank.ne.0) return
      n13 = 0
      i13 = 0
      n14 = 0
      i14 = 0
      n15 = 0
      i15 = 0
c
c     loop over all atoms finding all the 1-3 relationships;
c     note "n12" and "i12" have already been setup elsewhere
c
      do i = 1, n
         n13(i) = 0
         do j = 1, n12(i)
            jj = i12(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (kk .eq. i)  goto 10
               do m = 1, n12(i)
                  if (kk .eq. i12(m,i))  goto 10
               end do
               n13(i) = n13(i) + 1
               i13(n13(i),i) = kk
   10          continue
            end do
         end do
         if ((n13(i) .gt. maxn13).and.(rank.eq.0)) then
            write (iout,20)  i
   20       format (/,' ATTACH  --  Too many 1-3 Connected Atoms',
     &                 ' Attached to Atom',i6)
            call fatal
         end if
         call sort (n13(i),i13(:,i))
      end do
c
c
c     loop over all atoms finding all the 1-4 relationships
c
      do i = 1, n
         n14(i) = 0
         do j = 1, n13(i)
            jj = i13(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (kk .eq. i)  goto 30
               do m = 1, n12(i)
                  if (kk .eq. i12(m,i))  goto 30
               end do
               do m = 1, n13(i)
                  if (kk .eq. i13(m,i))  goto 30
               end do
               n14(i) = n14(i) + 1
               i14(n14(i),i) = kk
   30          continue
            end do
         end do
         if ((n14(i) .gt. maxn14).and.(rank.eq.0)) then
            write (iout,40)  i
   40       format (/,' ATTACH  --  Too many 1-4 Connected Atoms',
     &                 ' Attached to Atom',i6)
            call fatal
         end if
         call sort (n14(i),i14(:,i))
      end do
c
c     loop over all atoms finding all the 1-5 relationships
c
      do i = 1, n
         n15(i) = 0
         do j = 1, n14(i)
            jj = i14(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (kk .eq. i)  goto 50
               do m = 1, n12(i)
                  if (kk .eq. i12(m,i))  goto 50
               end do
               do m = 1, n13(i)
                  if (kk .eq. i13(m,i))  goto 50
               end do
               do m = 1, n14(i)
                  if (kk .eq. i14(m,i))  goto 50
               end do
               n15(i) = n15(i) + 1
               i15(n15(i),i) = kk
   50          continue
            end do
         end do
         if ((n15(i) .gt. maxn15).and.(rank.eq.0)) then
            write (iout,60)  i
   60       format (/,' ATTACH  --  Too many 1-5 Connected Atoms',
     &                 ' Attached to Atom',i6)
            call fatal
         end if
         call sort (n15(i),i15(:,i))
      end do
      return
      end
c
c     subroutine alloc_shared_attach : allocate shared memory pointers for attach
c     parameter arrays
c
      subroutine alloc_shared_attach
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use couple
      use domdec
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
c      if (associated(i13)) deallocate(i13)
c      if (associated(n13)) deallocate(n13)
c      if (associated(i14)) deallocate(i14)
c      if (associated(n14)) deallocate(n14)
c      if (associated(i15)) deallocate(i15)
c      if (associated(n15)) deallocate(n15)
      if (associated(i13)) then
        CALL MPI_Win_shared_query(wini13, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wini13,ierr)
      end if
      if (associated(n13)) then
        CALL MPI_Win_shared_query(winn13, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winn13,ierr)
      end if
      if (associated(i14)) then
        CALL MPI_Win_shared_query(wini14, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wini14,ierr)
      end if
      if (associated(n14)) then
        CALL MPI_Win_shared_query(winn14, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winn14,ierr)
      end if
      if (associated(i15)) then
        CALL MPI_Win_shared_query(wini15, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wini15,ierr)
      end if
      if (associated(n15)) then
        CALL MPI_Win_shared_query(winn15, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winn15,ierr)
      end if
c
c    i13
c
      arrayshape2=(/maxn13,n/)
      if (hostrank == 0) then
        windowsize = int(maxn13*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wini13, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wini13, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,i13,arrayshape2)
c
c    n13
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winn13, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winn13, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,n13,arrayshape)
c
c    i14
c
      arrayshape2=(/maxn14,n/)
      if (hostrank == 0) then
        windowsize = int(maxn14*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wini14, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wini14, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,i14,arrayshape2)
c
c    n14
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winn14, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winn14, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,n14,arrayshape)
c
c    i15
c
      arrayshape2=(/maxn15,n/)
      if (hostrank == 0) then
        windowsize = int(maxn15*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wini15, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wini15, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,i15,arrayshape2)
c
c    n14
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winn15, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winn15, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,n15,arrayshape)
      return
      end
