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
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'bond.i'
      include 'couple.i'
      include 'iounit.i'
      include 'tors.i'
      include 'openmp.i'
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
c       allocate arrays
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
      if (associated(torsglob)) deallocate (torsglob)
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
c
c     subroutine alloc_shared_bond : allocate shared memory pointers for torsions
c     parameter arrays
c
      subroutine alloc_shared_torsions
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use mpi
      implicit none
      include 'sizes.i'
      include 'tors.i'
      include 'improp.i'
      include 'imptor.i'
      include 'pitors.i'
      include 'strtor.i'
      include 'openmp.i'
      integer :: win,win2
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
      if (associated(itors)) deallocate (itors)
      if (associated(tors1)) deallocate (tors1)
      if (associated(tors2)) deallocate (tors2)
      if (associated(tors3)) deallocate (tors3)
      if (associated(tors4)) deallocate (tors4)
      if (associated(tors5)) deallocate (tors5)
      if (associated(tors6)) deallocate (tors6)
      if (associated(kprop)) deallocate (kprop)
      if (associated(vprop)) deallocate (vprop)
      if (associated(iiprop)) deallocate (iiprop)
      if (associated(itors1)) deallocate (itors1)
      if (associated(itors2)) deallocate (itors2)
      if (associated(itors3)) deallocate (itors3)
      if (associated(iitors)) deallocate (iitors)
      if (associated(iiprop)) deallocate (iiprop)
      if (associated(kpit)) deallocate (kpit)
      if (associated(ipit)) deallocate (ipit)
      if (associated(ist)) deallocate (ist)
      if (associated(kst)) deallocate (kst)
      if (associated(nbstrtor)) deallocate (nbstrtor)
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iitors,arrayshape2)
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iiprop,arrayshape2)
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ipit,arrayshape2)
c
c     ist
c
      arrayshape2=(/2,ntors/)
      if (hostrank == 0) then
        windowsize = int(2*ntors,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
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
      CALL C_F_POINTER(baseptr,ist,arrayshape2)
c
c     kst
c
      arrayshape2=(/3,ntors/)
      if (hostrank == 0) then
        windowsize = int(3*ntors,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
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
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbstrtor,arrayshape)
      return
      end
