!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##########################################################
!     ##                                                      ##
!     ##  subroutine bitors  --  locate and store bitorsions  ##
!     ##                                                      ##
!     ##########################################################
!
!
!     "bitors" finds the total number of bitorsions, pairs of
!     overlapping dihedral angles, and the numbers of the five
!     atoms defining each bitorsion
!
!
!> @brief 
!>"bitors" finds the total number of bitorsions, pairs of
!>overlapping dihedral angles, and the numbers of the five
!>atoms defining each bitorsion
!> @param no params
subroutine bitors
   use angle
   use atmlst
   use bitor
   use couple
   use domdec
   use inform
   use iounit
   implicit none
   integer i,j,k
   integer ia,ib,ic,id,ie
!
   if (deb_Path) write(iout,*), 'bitors '
!
!
!     loop over all angles, storing the atoms in each bitorsion
!
   nbitor = 0
   do i = 1, nangle
      ib = iang(1,i)
      ic = iang(2,i)
      id = iang(3,i)
      do j = 1, n12(ib)
         ia = i12(j,ib)
         if (ia.ne.ic .and. ia.ne.id) then
            do k = 1, n12(id)
               ie = i12(k,id)
               if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
                  nbitor = nbitor + 1
                  if (nbitor .gt. maxbitor) then
                     if (rank.eq.0) write (iout,10)
10                   format (/,' BITORS  --  Too many Adjacent',&
                     &' Torsions; Increase MAXBITOR')
                     call fatal
                  end if
               end if
            end do
         end if
      end do
   end do
!
!     deallocate global pointers if necessary
!
   call dealloc_shared_bitors
!
!     allocate global pointers
!
   call alloc_shared_bitors
!
   nbitorloc = 0
   do i = 1, nangle
      ib = iang(1,i)
      ic = iang(2,i)
      id = iang(3,i)
      nbbitors(i) = nbitorloc
      do j = 1, n12(ib)
         ia = i12(j,ib)
         if (ia.ne.ic .and. ia.ne.id) then
            do k = 1, n12(id)
               ie = i12(k,id)
               if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
                  nbitorloc = nbitorloc + 1
                  ibitor(1,nbitorloc) = ia
                  ibitor(2,nbitorloc) = ib
                  ibitor(3,nbitorloc) = ic
                  ibitor(4,nbitorloc) = id
                  ibitor(5,nbitorloc) = ie
               end if
            end do
         end if
      end do
   end do
   return
end
!
!     subroutine bitors_update : update local bitors
!
!> @brief 
!> update local bitors
!> @param no params
subroutine bitors_update
   use angle
   use atmlst
   use bitor
   use couple
   use domdec
   use inform
   use iounit
   implicit none
   integer i,j,k,iangle,nbitorloc1,bitorscount
   integer ia,ib,ic,id,ie
!
   if (deb_Path) write(iout,*), 'bitors_update '
!
!
   if (allocated(bitorsglob)) deallocate(bitorsglob)
   allocate (bitorsglob(8*nbloc))
   nbitorloc = 0
   do i = 1, nangleloc
      iangle = angleglob(i)
      bitorscount = nbbitors(iangle)
      ib = iang(1,iangle)
      ic = iang(2,iangle)
      id = iang(3,iangle)
      nbitorloc1 = 0
      do j = 1, n12(ib)
         ia = i12(j,ib)
         if (ia.ne.ic .and. ia.ne.id) then
            do k = 1, n12(id)
               ie = i12(k,id)
               if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
                  nbitorloc = nbitorloc + 1
                  nbitorloc1 = nbitorloc1 + 1
                  bitorsglob(nbitorloc) = bitorscount + nbitorloc1
               end if
            end do
         end if
      end do
   end do
   return
end
!
!     subroutine dealloc_shared_bitors : deallocate shared memory pointers for bitors
!     parameter arrays
!
!> @brief 
!>deallocate shared memory pointers for bitors
!>parameter arrays
!> @param no params
subroutine dealloc_shared_bitors
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use bitor
   use tortor
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
!
   if (associated(ibitor)) then
      CALL MPI_Win_shared_query(winibitor, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winibitor,ierr)
   end if
   if (associated(itt)) then
      CALL MPI_Win_shared_query(winitt, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winitt,ierr)
   end if
   if (associated(nbtortor)) then
      CALL MPI_Win_shared_query(winnbtortor, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winnbtortor,ierr)
   end if
   return
end
!
!
!     subroutine alloc_shared_bitors : allocate shared memory pointers for bitors
!     parameter arrays
!
!> @brief 
!>allocate shared memory pointers for bitors parameter arrays
!> @param no params
subroutine alloc_shared_bitors
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use sizes
   use bitor
   use domdec
   use tortor
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   integer :: arrayshape(1),arrayshape2(2)
!
!     ibitor
!
   arrayshape2=(/5,nbitor/)
   if (hostrank == 0) then
      windowsize = int(5*nbitor,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winibitor, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winibitor, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,ibitor,arrayshape2)
!
!     itt
!
   arrayshape2=(/3,nbitor/)
   if (hostrank == 0) then
      windowsize = int(3*nbitor,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winitt, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winitt, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,itt,arrayshape2)
!
!     nbtortor
!
   arrayshape=(/nbitor/)
   if (hostrank == 0) then
      windowsize = int(3*nbitor,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winnbtortor, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winnbtortor, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,nbtortor,arrayshape)
   return
end
