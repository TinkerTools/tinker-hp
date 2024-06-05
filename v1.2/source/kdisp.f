c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kdisp  --  dispersion parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kdisp" assigns C6 coefficients and damping parameters for
c     dispersion interactions and processes any new or changed
c     values for these parameters
c
c
      subroutine kdisp
      use atoms
      use atmlst
      use atmtyp
      use cutoff
      use disp
      use domdec
      use dsppot
      use inform
      use iounit
      use kdsp
      use keys
      use neigh
      use pme
      use potent
      use sizes
      implicit none
      integer i,k
      integer ia,ic,next
      real*8 cs,adsp
      real*8 csixi
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
      if (deb_Path) write(iout,*), 'kdisp '
c
c    
c
c     process keywords containing damped dispersion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'DISPERSION ') then
            k = 0
            cs = 0.0d0
            adsp = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  cs,adsp
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Damped Dispersion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',16x,'C6',12x,'Damp',/)
               end if
               if (k .le. maxclass) then
                  dspsix(k) = cs
                  dspdmp(k) = adsp
                  if (.not. silent) then
                     write (iout,30)  k,cs,adsp
   30                format (6x,i6,7x,f15.4,f15.4)
                  end if
               else
                  write (iout,40)
   40             format (/,' KDISP  --  Too many Damped',
     &                       ' Dispersion Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     deallocate global pointers if necessary
c
      call dealloc_shared_disp
c
c     allocate global pointers
c
      call alloc_shared_disp
c
c     assign the dispersion C6 values and alpha parameters 
c
      do i = 1, n
         ic = class(i)
         csix(i) = dspsix(ic)
         adisp(i) = dspdmp(ic)
      end do
c
c     process keywords containing atom specific dispersion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'DISPERSION ') then
            ia = 0
            cs = 0.0d0
            adsp = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,cs,adsp
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Dispersion Values for',
     &                       ' Specific Atoms :',
     &                    //,8x,'Atom',19x,'C6',12x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,cs,adsp
   60             format (6x,i6,7x,f15.4,f15.4)
               end if
               csix(ia) = cs
               adisp(ia) = adsp
            end if
   70       continue
         end if
      end do
c 
c     remove zero and undefined dispersion sites from the list
c     
      ndisp = 0
      do i = 1, n
         if (csix(i) .ne. 0.0d0) then 
            nbdisp(i) = ndisp
            ndisp = ndisp + 1
            idisp(ndisp) = i
            displist(i) = ndisp
            csix(ndisp) = csix(i)
            adisp(ndisp) = adisp(i)
         end if
      end do
c
c     compute pairwise sum of C6 coefficients needed for PME
c
      csixpr = 0.0d0
      if (use_dewald) then
         do i = 1, ndisp
            csixi = csix(i)
            do k = 1, ndisp
               csixpr = csixpr + csixi*csix(k)
            end do
         end do
      end if
c
c     turn off the dispersion potential if not used
c
      if (ndisp .eq. 0)  then
        use_disp = .false.
        use_dlist = .false.
      end if

      if (allocated(displocnl)) deallocate(displocnl)
      allocate (displocnl(ndisp))
      if (use_dewald) then
        if (allocated(disprecglob)) deallocate(disprecglob)
        allocate (disprecglob(n))
        if (allocated(disprecloc)) deallocate(disprecloc)
        allocate (disprecloc(n))
      end if
      if (allocated(dispglob)) deallocate(dispglob)
      allocate (dispglob(n))

      return
      end
c
c     subroutine kdisp_update
c
      subroutine kdisp_update(istep)
      use atoms
      use atmlst
      use atmtyp
      use cutoff
      use disp
      use domdec
      use dsppot
      use inform
      use iounit
      use kdsp
      use keys
      use neigh
      use pme
      use potent
      use sizes
      implicit none
      integer i,istep,iglob,iproc,idisploc
      integer dispcount
      integer modnl
      real*8 d
c
      if (deb_Path) write(iout,*), 'kdisp_update '
c
c
c
c     remove zero and undefined dispersion sites from the list
c       
      ndisploc = 0
      do i = 1, nloc
        iglob = glob(i)
        dispcount = nbdisp(iglob)
         if (csix(iglob) .ne. 0.0d0) then 
            ndisploc = ndisploc + 1
            dispglob(ndisploc) = dispcount + 1
         end if
      end do

      ndispbloc = ndisploc
      do iproc = 1, n_recep2
        do i = 1, domlen(p_recep2(iproc)+1)
          iglob = glob(bufbeg(p_recep2(iproc)+1)+i-1)
          dispcount = nbdisp(iglob)
          if (csix(iglob) .ne. 0.0d0) then 
            ndispbloc = ndispbloc + 1
            dispglob(ndispbloc) = dispcount + 1
          end if
        end do
      end do

      if (use_dewald) then
        ndisprecloc = 0
        do i = 1, nlocrec
           iglob = globrec(i)
           idisploc = displist(iglob)
           if (idisploc.eq.0) cycle
           ndisprecloc = ndisprecloc + 1
           disprecglob(ndisprecloc) = idisploc
           disprecloc(idisploc) = ndisprecloc
        end do
      end if
c
      modnl = mod(istep,ineigup)
      if (istep.eq.-1) return
      if (modnl.ne.0) return
      if (allocated(dispglobnl)) deallocate(dispglobnl)
      allocate (dispglobnl(nlocnl))

      ndisplocnl = 0
      do i = 1, nlocnl
        iglob = ineignl(i)
        dispcount = nbdisp(iglob)
        if (csix(iglob) .ne. 0.0d0) then
          call distprocpart(iglob,rank,d,.true.)
          if (repart(iglob).eq.rank) d = 0.0d0
            if (d*d.le.(dbuf2/4)) then
              ndisplocnl = ndisplocnl + 1
              dispglobnl(ndisplocnl) = dispcount + 1
              displocnl(dispcount+1) = ndisplocnl
            end if
        end if
      end do

      return
      end
c
c     subroutine dealloc_shared_disp : deallocate shared memory pointers for dispersion
c     parameter arrays
c
      subroutine dealloc_shared_disp
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use disp
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
c
      if (associated(idisp)) then
        CALL MPI_Win_shared_query(winidisp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winidisp,ierr)
      end if
      if (associated(csix)) then
        CALL MPI_Win_shared_query(wincsix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wincsix,ierr)
      end if
      if (associated(adisp)) then
        CALL MPI_Win_shared_query(winadisp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winadisp,ierr)
      end if
      if (associated(displist)) then
        CALL MPI_Win_shared_query(windisplist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(windisplist,ierr)
      end if
      if (associated(nbdisp)) then
        CALL MPI_Win_shared_query(winnbdisp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbdisp,ierr)
      end if
      return
      end
c
c     subroutine alloc_shared_disp : allocate shared memory pointers for dispersion
c     parameter arrays
c
      subroutine alloc_shared_disp
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use domdec
      use disp
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1)
c
c     idisp
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
     $  hostcomm, baseptr, winidisp, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winidisp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,idisp,arrayshape)
c
c     csix
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
     $  hostcomm, baseptr, wincsix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wincsix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,csix,arrayshape)
c
c     adisp
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
     $  hostcomm, baseptr, winadisp, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winadisp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,adisp,arrayshape)
c
c     displist
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
     $  hostcomm, baseptr, windisplist, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(windisplist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,displist,arrayshape)
c
c     nbdisp
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
     $  hostcomm, baseptr, winnbdisp, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbdisp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbdisp,arrayshape)
      return
      end
