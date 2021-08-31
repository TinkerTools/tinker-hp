c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine bonds  --  locate and store covalent bonds  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bonds" finds the total number of covalent bonds and
c     stores the atom numbers of the atoms defining each bond
c
c
      subroutine bonds(init)
      use atmlst
      use atoms
      use bond
      use couple
      use domdec
      use iounit
      implicit none
      integer i,j,k,m
      integer iglob
      real*8 xi,yi,zi,xk,yk,zk
      logical init,docompute
   10              format (/,' BONDS  --  Too many Bonds; Increase',
     &                        ' MAXBND')
c
      if (init) then
c
c
c     loop over all atoms, storing the atoms in each bond
c
        nbond = 0
        do i = 1, n
           do j = 1, n12(i)
              k = i12(j,i)
              if (i.lt.k) then
                nbond = nbond + 1
                if (nbond .gt. 4*n) then
                   if (rank.eq.0) write (iout,10)
                   call fatal
                end if
              end if
           end do
        end do
c
c       deallocate global pointers if necessary
c
        call dealloc_shared_bond
c
c       allocate global pointers
c
        call alloc_shared_bond
c
        nbondloc = 0
        do i = 1, n
           do j = 1, n12(i)
              k = i12(j,i)
              if (i.lt.k) then
                nbondloc = nbondloc + 1
                ibnd(1,nbondloc) = i
                ibnd(2,nbondloc) = k
                bndlist(j,i) = nbondloc
                do m = 1, n12(k)
                   if (i .eq. i12(m,k)) then
                      bndlist(m,k) = nbondloc
                      goto 20
                   end if
                end do
   20           continue
              end if
           end do
        end do
      end if
      if (allocated(bndglob)) deallocate(bndglob)
      allocate (bndglob(maxvalue*nbloc))
      nbondloc = 0
      do i = 1, nloc
         iglob = glob(i)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         do j = 1, n12(iglob)
            k = i12(j,iglob)
            xk = x(k)
            yk = y(k)
            zk = z(k)
            call halfcell(xi,yi,zi,xk,yk,zk,docompute)
            if (docompute) then
              nbondloc = nbondloc + 1
              bndglob(nbondloc) = bndlist(j,iglob)
            end if
         end do
      end do
      return
      end
c
c     subroutine dealloc_shared_bond : deallocate shared memory pointers for bond
c     parameter arrays
c
      subroutine dealloc_shared_bond
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use atmlst
      use bond
      use pitors
      use tors
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      if (associated(bndlist)) then
        CALL MPI_Win_shared_query(winbndlist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winbndlist,ierr)
      end if
      if (associated(bk)) then
        CALL MPI_Win_shared_query(winbk, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winbk,ierr)
      end if
      if (associated(bl)) then
        CALL MPI_Win_shared_query(winbl, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winbl,ierr)
      end if
      if (associated(ibnd)) then
        CALL MPI_Win_shared_query(winibnd, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winibnd,ierr)
      end if
      if (associated(nbtors)) then
        CALL MPI_Win_shared_query(winnbtors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbtors,ierr)
      end if
      if (associated(nbpitors)) then
        CALL MPI_Win_shared_query(winnbpitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbpitors,ierr)
      end if
      return
      end
c
c
c     subroutine alloc_shared_bond : allocate shared memory pointers for bond
c     parameter arrays
c
      subroutine alloc_shared_bond
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atmlst
      use atoms
      use bond
      use domdec
      use pitors
      use tors
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
c     bk
c
      arrayshape=(/nbond/)
      if (hostrank == 0) then
        windowsize = int(nbond,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winbk, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winbk, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,bk,arrayshape)
c
c     bl
c
      arrayshape=(/nbond/)
      if (hostrank == 0) then
        windowsize = int(nbond,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winbl, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winbl, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,bl,arrayshape)
c
c     bndlist
c
      arrayshape2=(/8,n/)
      if (hostrank == 0) then
        windowsize = int(8*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winbndlist, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winbndlist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,bndlist,arrayshape2)
c
c    ibnd
c
      arrayshape2=(/2,nbond/)
      if (hostrank == 0) then
        windowsize = int(2*nbond,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winibnd, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winibnd, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ibnd,arrayshape2)
c
c    nbtors
c
      arrayshape=(/nbond/)
      if (hostrank == 0) then
        windowsize = int(nbond,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbtors, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbtors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbtors,arrayshape)
c
c    nbpitors
c
      arrayshape=(/nbond/)
      if (hostrank == 0) then
        windowsize = int(nbond,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbpitors, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbpitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbpitors,arrayshape)
      return
      end
