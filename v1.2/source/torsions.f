c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine torsions  --  locate and store torsions  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "torsions" finds the total number of dihedral angles and
c     the numbers of the four atoms defining each dihedral angle
c
c
      subroutine torsions(init)
      use atmlst
      use atoms
      use bond
      use couple
      use domdec
      use iounit
      use tors
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
      integer ibond,torscount,ntorsloc1
      logical init
c
      if (init) then
c
c     loop over all bonds, storing the atoms in each torsion
c
        ntors = 0
        do i = 1, nbond
           ib = ibnd(1,i)
           ic = ibnd(2,i)
           do j = 1, n12(ib)
              ia = i12(j,ib)
              if (ia .ne. ic) then
                 do k = 1, n12(ic)
                    id = i12(k,ic)
                    if (id.ne.ib .and. id.ne.ia) then
                       ntors = ntors + 1
                       if (ntors .gt. 8*n) then
                          if (rank.eq.0) write (iout,10)
   10                     format (/,' TORSIONS  --  Too many Torsional',
     &                               ' Angles; Increase MAXTORS')
                          call fatal
                       end if
                    end if
                 end do
              end if
           end do
        end do
c
c       deallocate global pointers if necessary
c
        call dealloc_shared_torsions
c
c       allocate global pointers
c
        call alloc_shared_torsions
c
        ntorsloc = 0
        do i = 1, nbond
           ib = ibnd(1,i)
           ic = ibnd(2,i)
           nbtors(i) = ntorsloc
           do j = 1, n12(ib)
              ia = i12(j,ib)
              if (ia .ne. ic) then
                 do k = 1, n12(ic)
                    id = i12(k,ic)
                    if (id.ne.ib .and. id.ne.ia) then
                       ntorsloc = ntorsloc + 1
                       itors(1,ntorsloc) = ia
                       itors(2,ntorsloc) = ib
                       itors(3,ntorsloc) = ic
                       itors(4,ntorsloc) = id
                    end if
                 end do
              end if
           end do
        end do
      end if
      if (allocated(torsglob)) deallocate (torsglob)
      allocate (torsglob(6*nbloc))
      ntorsloc = 0
      do i = 1, nbondloc
        ibond = bndglob(i)
        torscount = nbtors(ibond)
        ib = ibnd(1,ibond)
        ic = ibnd(2,ibond)
        ntorsloc1  = 0
        do j = 1, n12(ib)
           ia = i12(j,ib)
           if (ia .ne. ic) then
              do k = 1, n12(ic)
                 id = i12(k,ic)
                 if (id.ne.ib .and. id.ne.ia) then
                    ntorsloc = ntorsloc + 1
                    ntorsloc1 = ntorsloc1 + 1
                    torsglob(ntorsloc) = torscount + ntorsloc1
                 end if
              end do
           end if
        end do
      end do
      return
      end
c
c     subroutine dealloc_shared_torsions : deallocate shared memory pointers for torsions
c     parameter arrays
c
      subroutine dealloc_shared_torsions
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use angtor
      use improp
      use imptor
      use pitors
      use strtor
      use tors
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
c
      if (associated(itors)) then
        CALL MPI_Win_shared_query(winitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winitors,ierr)
      end if
      if (associated(tors1)) then
        CALL MPI_Win_shared_query(wintors1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintors1,ierr)
      end if
      if (associated(tors2)) then
        CALL MPI_Win_shared_query(wintors2, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintors2,ierr)
      end if
      if (associated(tors3)) then
        CALL MPI_Win_shared_query(wintors3, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintors3,ierr)
      end if
      if (associated(tors4)) then
        CALL MPI_Win_shared_query(wintors4, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintors4,ierr)
      end if
      if (associated(tors5)) then
        CALL MPI_Win_shared_query(wintors5, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintors5,ierr)
      end if
      if (associated(tors6)) then
        CALL MPI_Win_shared_query(wintors6, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintors6,ierr)
      end if
      if (associated(kprop)) then
        CALL MPI_Win_shared_query(winkprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winkprop,ierr)
      end if
      if (associated(vprop)) then
        CALL MPI_Win_shared_query(winvprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winvprop,ierr)
      end if
      if (associated(iiprop)) then
        CALL MPI_Win_shared_query(winiiprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winiiprop,ierr)
      end if
      if (associated(itors1)) then
        CALL MPI_Win_shared_query(winitors1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winitors1,ierr)
      end if
      if (associated(itors2)) then
        CALL MPI_Win_shared_query(winitors2, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winitors2,ierr)
      end if
      if (associated(itors3)) then
        CALL MPI_Win_shared_query(winitors3, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winitors3,ierr)
      end if
      if (associated(iitors)) then
        CALL MPI_Win_shared_query(winiitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winiitors,ierr)
      end if
      if (associated(kpit)) then
        CALL MPI_Win_shared_query(winkpit, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winkpit,ierr)
      end if
      if (associated(ipit)) then
        CALL MPI_Win_shared_query(winipit, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winipit,ierr)
      end if
      if (associated(ist)) then
        CALL MPI_Win_shared_query(winist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winist,ierr)
      end if
      if (associated(kst)) then
        CALL MPI_Win_shared_query(winkst, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winkst,ierr)
      end if
      if (associated(nbstrtor)) then
        CALL MPI_Win_shared_query(winnbstrtor, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbstrtor,ierr)
      end if
      if (associated(iat)) then
        CALL MPI_Win_shared_query(winiat, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winiat,ierr)
      end if
      if (associated(kant)) then
        CALL MPI_Win_shared_query(winkant, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winkant,ierr)
      end if
      if (associated(nbangtor)) then
        CALL MPI_Win_shared_query(winnbangtor, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbangtor,ierr)
      end if
      return
      end
c
c
c     subroutine alloc_shared_torsions : allocate shared memory pointers for torsions
c     parameter arrays
c
      subroutine alloc_shared_torsions
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use angtor
      use domdec
      use improp
      use imptor
      use pitors
      use strtor
      use tors
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
c     itors
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winitors, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,itors,arrayshape2)
c
c     tors1
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintors1, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintors1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,tors1,arrayshape2)
c
c     tors2
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintors2, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintors2, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,tors2,arrayshape2)
c
c     tors3
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintors3, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintors3, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,tors3,arrayshape2)
c
c     tors4
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintors4, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintors4, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,tors4,arrayshape2)
c
c     tors5
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintors5, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintors5, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,tors5,arrayshape2)
c
c     tors6
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintors6, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintors6, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,tors6,arrayshape2)
c
c     kprop
c
      arrayshape=(/ntors/)
      if (hostrank == 0) then
        windowsize = int(ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winkprop, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winkprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,kprop,arrayshape)
c
c     vprop
c
      arrayshape=(/ntors/)
      if (hostrank == 0) then
        windowsize = int(ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winvprop, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winvprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,vprop,arrayshape)
c
c     iiprop
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winiiprop, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winiiprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iiprop,arrayshape2)
c
c     itors1
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winitors1, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winitors1, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,itors1,arrayshape2)
c
c     itors2
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winitors2, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winitors2, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,itors2,arrayshape2)
c
c     itors3
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winitors3, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winitors3, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,itors3,arrayshape2)
c
c     iitors
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winiitors, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winiitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iitors,arrayshape2)
c
c     kpit
c
      arrayshape=(/ntors/)
      if (hostrank == 0) then
        windowsize = int(ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winkpit, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winkpit, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,kpit,arrayshape)
c
c     ipit
c
      arrayshape2=(/6,ntors/)
      if (hostrank == 0) then
        windowsize = int(6*ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winipit, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winipit, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ipit,arrayshape2)
c
c     ist
c
      arrayshape2=(/4,ntors/)
      if (hostrank == 0) then
        windowsize = int(4*ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winist, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ist,arrayshape2)
c
c     kst
c
      arrayshape2=(/9,ntors/)
      if (hostrank == 0) then
        windowsize = int(9*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winkst, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winkst, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,kst,arrayshape2)
c
c     nbstrtor
c
      arrayshape=(/ntors/)
      if (hostrank == 0) then
        windowsize = int(ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbstrtor, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbstrtor, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbstrtor,arrayshape)
c
c     iat
c
      arrayshape2=(/3,ntors/)
      if (hostrank == 0) then
        windowsize = int(3*ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winiat, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winiat, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iat,arrayshape2)
c
c     kant
c
      arrayshape2=(/6,ntors/)
      if (hostrank == 0) then
        windowsize = int(6*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winkant, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winkant, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,kant,arrayshape2)
c
c     nbangtor
c
      arrayshape=(/ntors/)
      if (hostrank == 0) then
        windowsize = int(ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbangtor, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbangtor, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbangtor,arrayshape)
      return
      end
