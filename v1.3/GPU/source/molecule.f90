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
#include "tinker_macro.h"
subroutine molecule(init)
   use sizes
   use atmlst
   use atmtyp
   use atoms
   use couple
   use domdec
   use inform
   use molcul
   use mpi
#ifdef _OPENACC
   !use thrust
#endif
   use utilgpu
   implicit none
   integer i,j,k,ii
   integer mi,mj,mk
   integer iglob,jglob,ierr
   integer iimol,jmol
   integer captur
   integer, allocatable :: list(:)
   logical init
   logical not_find
!
   if (init) then
!
!     allocate global arrays
!
      if (deb_Path) print*,'molecule init'
      call alloc_shared_mol
!
!       only master of the node fill the arrays
!
      if (hostrank.ne.0) goto 20
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
         if (mass(i) .ne. 0.0_ti_p)  goto 10
      end do
      do i = 1, n
         mass(i) = 1.0_ti_p
      end do
10    continue
!
!       compute the mass of each molecule and the total mass
!
      totmass = 0.0_ti_p
      do i = 1, nmol
         molmass(i) = 0.0_ti_p
         do k = imol(1,i), imol(2,i)
            molmass(i) = molmass(i) + mass(kmol(k))
         end do
         totmass = totmass + molmass(i)
      end do
20    call MPI_barrier(hostcomm,ierr)
      call MPI_BCAST(nmol,1,MPI_INT,0,hostcomm,ierr)
      call update_device_mol
      if (nproc.gt.1) return
!
   end if

   call prmem_request(molculeglob,nmol,async=.true.)

   if (nbloc.eq.n) then

      ! Copy molecule in molculeglob
      nmoleloc = nmol
!$acc parallel loop default(present) async
      do i = 1,nmol
         molculeglob(i) = i
      end do

   else

!!$acc data present(nmoleloc,molecule,glob)
!!$acc serial async
!         nmoleloc = 0
!!$acc end serial
!!$acc parallel loop present(moleculeglob)
!         do i = 1,nmol
!            molculeglob(i)=0
!         end do
!!$acc parallel loop present(molculeglob) async
!         do i = 1,nbloc
!            iglob = glob(i)
!            iimol = molcule(iglob)
!            molculeglob(iimol) = iimol
!         end do
!!$acc host_data use_device(moleculeglob)
!         call thrust_remove_zero_async_int(moleculeglob,nmol
!     &             ,rec_stream)
!!$acc end host_data
!!$acc parallel loop present(moleculeglob) async
!         do i = 1, nmol
!            if (moleculeglob(i).ne.0) then
!               nmoleloc = nmoleloc+1
!            end if
!         end do
!!$acc end data


!$acc serial present(nmoleloc) async
      nmoleloc = 0
!$acc end serial
!$acc parallel loop present(molcule,glob,molculeglob) &
!$acc         vector_length(32) &
!$acc         present(nmoleloc) &
!$acc         private(not_find) async
      do i = 1, nbloc
         iglob = glob(i)
         iimol = molcule(iglob)
         not_find = .true.
!$acc loop vector
         do j = 1, i-1
            jmol = molcule(glob(j))
            if (iimol.eq.jmol) then
!$acc atomic write
               not_find=.false.
            end if
         end do
!$acc loop vector
         do j = 1,1
            if (not_find) then
!$acc atomic capture
               nmoleloc = nmoleloc + 1
               captur   = nmoleloc
!$acc end atomic
               molculeglob(captur) = molcule(iglob)
            end if
         end do
      end do
!$acc update host(nmoleloc) async

   end if
!
end

subroutine update_device_mol
   use atmtyp
   use domdec ,only: rank
   use molcul
   use sizes  ,only: tinkerdebug
   use tinMemory
   implicit none

#ifdef _OPENACC
12 format(2x,'upload_device_mol')
   if(rank.eq.0.and.tinkerdebug) print 12
#endif
!$acc enter data create(nmoleloc)
!$acc update device(mass)
!$acc update device(molcule,imol,kmol,molmass)
end subroutine
!
!     subroutine alloc_shared_mol : allocate shared memory pointers for molecule
!     parameter arrays
!
subroutine alloc_shared_mol
   use sizes
   use atoms
   use domdec
   use molcul
   use mpi
   use tinMemory
   implicit none
!
   if (associated(kmol).and.n.eq.size(kmol)) return !Exit condition

   call shmem_request(molcule,winmolcule,[n],config=mhostacc)
   call shmem_request(kmol,winkmol,[n],config=mhostacc)
   call shmem_request(imol,winimol,[2,n],config=mhostacc)
   call shmem_request(molmass,winmolmass,[n],config=mhostacc)

end
