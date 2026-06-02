!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##########################################################
!     ##                                                      ##
!     ##  subroutine molecule  --  assign atoms to molecules  ##
!     ##                                                      ##
!     ##########################################################
!
!
!     "molecule" counts the molecules, assigns each atom to
!     its molecule and computes the mass of each molecule
!
!
!> @brief 
!> counts the molecules, assigns each atom to
!> its molecule and computes the mass of each molecule
!> @param[in] init: flag true if first call
subroutine molecule(init)
   use sizes
   use atmlst
   use atmtyp
   use atoms
   use couple
   use domdec
   use inform
   use iounit
   use molcul
   use mpi
   implicit none
   integer i,j,k,ii
   integer mi,mj,mk
   integer iglob,jglob,ierr
   integer, allocatable :: list(:)
   logical init
!
   if (deb_Path) write(iout,*), 'molecule '
!
!
   if (init) then
!
!     deallocate global pointers if necessary
!
      call dealloc_shared_mol
!
!     allocate global pointers
!
      call alloc_shared_mol
!
!       only master of the node fill the arrays
!
      if (hostrank.ne.0) goto 30
!
!
!       zero number of molecules and molecule membership list
!
      nmol = 0
      do i = 1, n
         molcule(i) = 0
      end do
!
!       assign each atom to its respective molecule
!
      do i = 1, n
         if (molcule(i) .eq. 0) then
            nmol = nmol + 1
            molcule(i) = nmol
         end if
         mi = molcule(i)
         do ii = 1, n12(i)
            j = i12(ii,i)
            mj = molcule(j)
            if (mj .eq. 0) then
               molcule(j) = mi
            else if (mi .lt. mj) then
               nmol = nmol - 1
               do k = 1, n
                  mk = molcule(k)
                  if (mk .eq. mj) then
                     molcule(k) = mi
                  else if (mk .gt. mj) then
                     molcule(k) = mk - 1
                  end if
               end do
            else if (mi .gt. mj) then
               nmol = nmol - 1
               do k = 1, n
                  mk = molcule(k)
                  if (mk .eq. mi) then
                     molcule(k) = mj
                  else if (mk .gt. mi) then
                     molcule(k) = mk - 1
                  end if
               end do
               mi = mj
            end if
         end do
      end do
!
!       perform dynamic allocation of some local arrays
!
      allocate (list(n))
!
!       pack atoms of each molecule into a contiguous indexed list
!
      do i = 1, n
         list(i) = molcule(i)
      end do
      call sort3 (n,list,kmol)
!
!       find the first and last atom in each molecule
!
      k = 1
      imol(1,1) = 1
      do i = 2, n
         j = list(i)
         if (j .ne. k) then
            imol(2,k) = i - 1
            k = j
            imol(1,k) = i
         end if
      end do
      imol(2,nmol) = n
!
!       perform deallocation of some local arrays
!
      deallocate (list)
!
!       sort the list of atoms in each molecule by atom number
!
      do i = 1, nmol
         k = imol(2,i) - imol(1,i) + 1
         call sort (k,kmol(imol(1,i)))
      end do
!
!       if all atomic masses are zero, set them all to unity
!
      do i = 1, n
         if (mass(i) .ne. 0.0d0)  goto 10
      end do
      do i = 1, n
         mass(i) = 1.0d0
      end do
10    continue
!
!       compute the mass of each molecule and the total mass
!
      totmass = 0.0d0
      do i = 1, nmol
         molmass(i) = 0.0d0
         do k = imol(1,i), imol(2,i)
            molmass(i) = molmass(i) + mass(kmol(k))
         end do
         totmass = totmass + molmass(i)
      end do
!
!
   end if
   if (allocated(molculeglob)) deallocate(molculeglob)
   allocate (molculeglob(nbloc))
   nmoleloc = 0
   do i = 1, nbloc
      iglob = glob(i)
      do j = 1, i-1
         jglob = glob(j)
         if (molcule(iglob).eq.molcule(jglob)) then
            goto 20
         end if
      end do
      nmoleloc = nmoleloc + 1
      molculeglob(nmoleloc) = molcule(iglob)
20    continue
   end do
30 call MPI_BCAST(nmol,1,MPI_INT,0,hostcomm,ierr)
!
   return
end
!
!     subroutine dealloc_shared_mol : deallocate shared memory pointers for molecule
!     parameter arrays
!
!> @brief 
!> deallocate shared memory pointers for molecule
!> parameter arrays
!> @param no params
subroutine dealloc_shared_mol
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use molcul
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
!
   if (associated(molcule)) then
      CALL MPI_Win_shared_query(winmolcule, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winmolcule,ierr)
   end if
   if (associated(kmol)) then
      CALL MPI_Win_shared_query(winkmol, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winkmol,ierr)
   end if
   if (associated(imol)) then
      CALL MPI_Win_shared_query(winimol, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winimol,ierr)
   end if
   if (associated(molmass)) then
      CALL MPI_Win_shared_query(winmolmass, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winmolmass,ierr)
   end if
   return
end
!
!     subroutine alloc_shared_mol : allocate shared memory pointers for molecule
!     parameter arrays
!
!> @brief 
!> allocate shared memory pointers for molecule
!> parameter arrays
!> @param no params
subroutine alloc_shared_mol
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use sizes
   use atoms
   use domdec
   use molcul
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   integer :: arrayshape(1),arrayshape2(2)
!
!     molcule
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
   &hostcomm, baseptr, winmolcule, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winmolcule, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,molcule,arrayshape)
!
!    kmol
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
   &hostcomm, baseptr, winkmol, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winkmol, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,kmol,arrayshape)
!
!    imol
!
   arrayshape2=(/2,n/)
   if (hostrank == 0) then
      windowsize = int(2*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winimol, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winimol, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,imol,arrayshape2)
!
!     molmass
!
   arrayshape=(/n/)
   if (hostrank == 0) then
      windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winmolmass, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winmolmass, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,molmass,arrayshape)
!
   return
end
