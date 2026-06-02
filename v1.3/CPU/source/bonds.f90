!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine bonds  --  locate and store covalent bonds  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "bonds" finds the total number of covalent bonds and
!     stores the atom numbers of the atoms defining each bond
!
!
!> @brief 
!> finds the total number of covalent bonds and
!> stores the atom numbers of the atoms defining each bond displays the total potential energy
!> @param no params
subroutine bonds
   use atmlst
   use atoms
   use bond
   use couple
   use domdec
   use inform
   use iounit
   implicit none
   integer i,j,k,m
10 format (/,' BONDS  --  Too many Bonds; Increase',&
   &' MAXBND')
!
!
   if (deb_Path) write(iout,*), 'bonds '
!
!
!     loop over all atoms, storing the atoms in each bond
!
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
!
!     deallocate global pointers if necessary
!
   call dealloc_shared_bond
!
!     allocate global pointers
!
   call alloc_shared_bond
!
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
20          continue
         end if
      end do
   end do
   return
end
!
!     subroutine bonds_update : update local bonds
!
!> @brief 
!> update local bonds
!> @param no params
subroutine bonds_update
   use atmlst
   use atoms
   use bond
   use couple
   use domdec
   use inform
   use iounit
   implicit none
   integer i,j,k,iglob
   logical docompute
   real*8 xk,yk,zk,xi,yi,zi
!
   if (deb_Path) write(iout,*), 'bonds_update '
!
!
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
!
!     subroutine dealloc_shared_bond : deallocate shared memory pointers for bond
!     parameter arrays
!
!> @brief 
!> deallocate shared memory pointers for bond
!> arameter arrays
!> @param no params
subroutine dealloc_shared_bond
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use atmlst
   use bond
   use bndpot
   use pitors
   use tors
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   if (associated(bndlist)) then
      CALL MPI_Win_shared_query(winbndlist, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winbndlist,ierr)
   end if
   if (associated(bk)) then
      CALL MPI_Win_shared_query(winbk, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winbk,ierr)
   end if
   if (associated(bl)) then
      CALL MPI_Win_shared_query(winbl, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winbl,ierr)
   end if
   if (associated(ba)) then
      CALL MPI_Win_shared_query(winba, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winba,ierr)
   end if
   if (associated(bndtyp)) then
      CALL MPI_Win_shared_query(winbndtyp, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winbndtyp,ierr)
   end if
   if (associated(ibnd)) then
      CALL MPI_Win_shared_query(winibnd, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winibnd,ierr)
   end if
   if (associated(nbtors)) then
      CALL MPI_Win_shared_query(winnbtors, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winnbtors,ierr)
   end if
   if (associated(nbpitors)) then
      CALL MPI_Win_shared_query(winnbpitors, 0, windowsize, disp_unit,&
      &baseptr, ierr)
      CALL MPI_Win_free(winnbpitors,ierr)
   end if
   return
end
!
!
!     subroutine alloc_shared_bond : allocate shared memory pointers for bond
!     parameter arrays
!
!> @brief 
!> allocate shared memory pointers for bond
!> parameter arrays
!> @param no params
subroutine alloc_shared_bond
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use sizes
   use atmlst
   use atoms
   use bond
   use bndpot
   use domdec
   use pitors
   use tors
   use mpi
   implicit none
   INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
   INTEGER :: disp_unit,ierr
   TYPE(C_PTR) :: baseptr
   integer :: arrayshape(1),arrayshape2(2)
!
!     bk
!
   arrayshape=(/nbond/)
   if (hostrank == 0) then
      windowsize = int(nbond,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winbk, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winbk, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,bk,arrayshape)
!
!     bl
!
   arrayshape=(/nbond/)
   if (hostrank == 0) then
      windowsize = int(nbond,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winbl, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winbl, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,bl,arrayshape)
!
!     ba
!
   arrayshape=(/nbond/)
   if (hostrank == 0) then
      windowsize = int(nbond,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winba, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winba, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,ba,arrayshape)
!
!     bndtyp
!
   arrayshape=(/nbond/)
   if (hostrank == 0) then
      windowsize = int(nbond,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winbndtyp, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winbndtyp, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,bndtyp,arrayshape)
!
!     bndlist
!
   arrayshape2=(/8,n/)
   if (hostrank == 0) then
      windowsize = int(8*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winbndlist, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winbndlist, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,bndlist,arrayshape2)
!
!    ibnd
!
   arrayshape2=(/2,nbond/)
   if (hostrank == 0) then
      windowsize = int(2*nbond,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winibnd, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winibnd, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,ibnd,arrayshape2)
!
!    nbtors
!
   arrayshape=(/nbond/)
   if (hostrank == 0) then
      windowsize = int(nbond,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winnbtors, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winnbtors, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,nbtors,arrayshape)
!
!    nbpitors
!
   arrayshape=(/nbond/)
   if (hostrank == 0) then
      windowsize = int(nbond,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
   else
      windowsize = 0_MPI_ADDRESS_KIND
   end if
   disp_unit = 1
!
!    allocation
!
   CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,&
   &hostcomm, baseptr, winnbpitors, ierr)
   if (hostrank /= 0) then
      CALL MPI_Win_shared_query(winnbpitors, 0, windowsize, disp_unit,&
      &baseptr, ierr)
   end if
!
!    association with fortran pointer
!
   CALL C_F_POINTER(baseptr,nbpitors,arrayshape)
   return
end
