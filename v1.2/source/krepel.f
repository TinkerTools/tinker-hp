c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine krepel  --  Pauli repulsion term assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "krepel" assigns the size values, exponential parameter and
c     number of valence electrons for Pauli repulsion interactions
c     and processes any new or changed values for these parameters
c
c
      subroutine krepel
      use atoms
      use atmtyp
      use inform
      use iounit
      use krepl
      use keys
      use potent
      use repel
      use sizes
      implicit none
      integer i,k
      integer ia,ic,next
      real*8 spr,apr,epr
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
      if (deb_Path) write(iout,*), 'krepel '
c
c
c     process keywords containing Pauli repulsion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'REPULSION ') then
            k = 0
            spr = 0.0d0
            apr = 0.0d0
            epr = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  spr,apr,epr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Pauli Repulsion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',15x,'Size',11x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (k .le. maxclass) then
                  prsiz(k) = spr
                  prdmp(k) = apr
                  prele(k) = -abs(epr)
                  if (.not. silent) then
                     write (iout,30)  k,spr,apr,epr
   30                format (6x,i6,7x,2f15.4,f15.3)
                  end if
               else
                  write (iout,40)
   40             format (/,' KREPEL  --  Too many Pauli Repulsion',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     deallocate global pointers if necessary
c
      call dealloc_shared_rep
c
c     allocate global pointers
c
      call alloc_shared_rep
c
c     assign the repulsion size, alpha and valence parameters 
c     
      do i = 1, n
         ic = class(i)
         sizpr(i) = prsiz(ic)
         dmppr(i) = prdmp(ic)
         elepr(i) = prele(ic)
      end do
c
c     process keywords containing atom specific Pauli repulsion
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'REPULSION ') then
            ia = 0
            spr = 0.0d0
            apr = 0.0d0
            epr = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,spr,apr,epr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Pauli Repulsion Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',17x,'Size',12x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,spr,apr,epr
   60             format (6x,i6,7x,2f15.4,f15.3)
               end if
               sizpr(ia) = spr
               dmppr(ia) = apr
               elepr(ia) = -abs(epr)
            end if
   70       continue
         end if
      end do
c
c     remove zero and undefined repulsion sites from the list
c
      nrep = 0
      do i = 1, n
         if (sizpr(i) .ne. 0.0d0) then
            nrep = nrep + 1
         end if
      end do
c
c     turn off the Pauli repulsion potential if not used
c
      if (nrep .eq. 0)  then
        use_repuls = .false.
      end if
      return
      end
c
c     subroutine dealloc_shared_rep : deallocate shared memory pointers for repulsion
c     parameter arrays
c
      subroutine dealloc_shared_rep
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use repel
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
c
      if (associated(sizpr)) then
        CALL MPI_Win_shared_query(winsizpr, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winsizpr,ierr)
      end if
      if (associated(dmppr)) then
        CALL MPI_Win_shared_query(windmppr, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(windmppr,ierr)
      end if
      if (associated(elepr)) then
        CALL MPI_Win_shared_query(winelepr, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winelepr,ierr)
      end if
      return
      end
c
c     subroutine alloc_shared_rep : allocate shared memory pointers for repulsion
c     parameter arrays
c
      subroutine alloc_shared_rep
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use domdec
      use repel
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1)
c
c     sizpr
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winsizpr, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winsizpr, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,sizpr,arrayshape)
c
c     dmppr
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, windmppr, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(windmppr, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,dmppr,arrayshape)
c
c     elepr
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winelepr, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winelepr, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,elepr,arrayshape)
      return
      end
