c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine bitors  --  locate and store bitorsions  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "bitors" finds the total number of bitorsions, pairs of
c     overlapping dihedral angles, and the numbers of the five
c     atoms defining each bitorsion
c
c
#include "tinker_precision.h"
      subroutine bitors(init)
      use angle
      use atmlst
      use bitor
      use couple
      use domdec
      use inform ,only: deb_Path
      use iounit
      use utilgpu
      use sizes
      implicit none
      integer i,j,k,iangle,nbitorloc1,bitorscount
      integer ia,ib,ic,id,ie
      logical init
      integer nbitorloc_capture
#ifdef USE_NVSHMEM_CUDA
      integer ind,ipe
#endif
c
      if (init) then
c
c     loop over all angles, storing the atoms in each bitorsion
c
        if (deb_Path) print*,'bitors init'
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
   10                     format (/,' BITORS  --  Too many Adjacent',
     &                               ' Torsions; Increase MAXBITOR')
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
        call alloc_shared_bitors
c
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

        call upload_device_bitors
      end if

      call prmem_request(bitorsglob,8*nbloc,async=.true.)
!$acc data present(nbitorloc,nangleloc) async
!$acc serial async
      nbitorloc = 0
!$acc end serial
!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc& present(angleglob,nbbitors,n12,i12,bitorsglob)
#else
!$acc& present(angleglob,nbbitors,iang,n12,i12,bitorsglob)
#endif
      do i = 1, nangleloc
         iangle = angleglob(i)
         bitorscount = nbbitors(iangle)
#ifdef USE_NVSHMEM_CUDA
         ipe =     (iangle-1)/nangle_pe
         ind = mod((iangle-1),nangle_pe) +1
         ib = d_iang(ipe)%pel(1,ind)
         ic = d_iang(ipe)%pel(2,ind)
         id = d_iang(ipe)%pel(3,ind)
#else
         ib = iang(1,iangle)
         ic = iang(2,iangle)
         id = iang(3,iangle)
#endif
         nbitorloc1 = 0
!$acc loop seq
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id) then
!$acc loop seq
               do k = 1, n12(id)
                  ie = i12(k,id)
                  if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
!$acc atomic capture
                     nbitorloc = nbitorloc + 1
                     nbitorloc_capture = nbitorloc
!$acc end atomic
                     nbitorloc1 = nbitorloc1 + 1
                     bitorsglob(nbitorloc_capture) = 
     &                      bitorscount + nbitorloc1
                  end if
               end do
            end if
         end do
      end do
!$acc update host(nbitorloc) async
!$acc end data

      end

      subroutine upload_device_bitors
      use bitor
      use tortor
      use domdec ,only: rank,hostcomm
      use inform ,only: deb_Path
      use mpi    ,only: MPI_BARRIER
      use nvshmem
      use tinMemory
      implicit none
      integer c_size
      integer ierr
#ifdef _OPENACC

 12   format(2x,'upload_device_bitors')
      if(deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
!$acc update device(nbbitors)
!$acc update device(ibitor)
!$acc enter data copyin(nbitorloc)
#endif
      end subroutine
c
c     subroutine alloc_shared_bitors : allocate shared memory pointers for bitors
c     parameter arrays
c
      subroutine alloc_shared_bitors
      use sizes
      use bitor
      use tortor
      use tinMemory
      implicit none

c     if(associated(itt).and.size(itt).eq.3*nbitor) return

      call shmem_request(ibitor,  winibitor,[5,nbitor],config=mhostacc)
      call shmem_request(itt,     winitt,   [3,nbitor],config=mhostacc)
      call shmem_request(nbtortor,winnbtortor,[nbitor],config=mhostacc)
      end
