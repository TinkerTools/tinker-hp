c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kangang  --  angle-angle parameter assignment  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kangang" assigns the parameters for angle-angle cross term
c     interactions and processes new or changed parameter values
c
c
      subroutine kangang(init)
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
      integer iglob
      integer nang,jang,kang,nangangloc1,angangcount
      real*8 fa,faa,aak(3)
      logical header
      character*20 keyword
      character*120 record
      character*120 string
      logical init
c
      if (init) then
c
c     process keywords containing angle-angle parameters
c
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
              string = record(next:120)
              read (string,*,err=10,end=10)  it,(aak(j),j=1,3)
   10         continue
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Angle-Angle Parameters :',
     &                      //,5x,'Atom Class',8x,'K(AA) 1',4x,
     &                         'K(AA) 2',4x,'K(AA) 3',/)
                 end if
                 if (rank.eq.0) write (iout,30)  it,(aak(j),j=1,3)
   30            format (9x,i3,7x,3f11.3)
              end if
              do j = 1, 3
                 anan(j,it) = aak(j)
              end do
           end if
        end do
c
c       assign the angle-angle parameters for each angle pair
c
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
c
c       allocate global arrays
c
        call alloc_shared_angang
c        if (associated(iaa)) deallocate (iaa)
c        allocate (iaa(2,nangang))
c        if (associated(kaa)) deallocate (kaa)
c        allocate (kaa(nangang))
c
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
c
c       turn off the angle-angle potential if it is not used
c
        if (nangang .eq. 0)  use_angang = .false.
      end if
      if (allocated(angangglob)) deallocate(angangglob)
      allocate (angangglob(ntors))
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
c
      return
      end
c
c     subroutine alloc_shared_angang : allocate shared memory pointers for angang
c     parameter arrays
c
      subroutine alloc_shared_angang
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use angang
      use domdec
      use mpi
      implicit none
      integer :: win,win2
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
      if (associated(iaa)) deallocate(iaa)
      if (associated(kaa)) deallocate(kaa)
c
c     iaa
c
      arrayshape2=(/2,nangang/)
      if (hostrank == 0) then
        windowsize = int(2*nangang,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iaa,arrayshape2)
c
c     kaa
c
      arrayshape=(/nangang/)
      if (hostrank == 0) then
        windowsize = int(nangang,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,kaa,arrayshape)
      return
      end
