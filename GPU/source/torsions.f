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
#include "tinker_precision.h"
      subroutine torsions(init)
      use atmlst
      use atoms
      use bond
      use couple
      use domdec
      use inform
      use iounit
      use nvshmem
      use tors
      use utilgpu
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
      integer ind,ipe
      integer ibond,torscount,ntorsloc1
      integer ntorsloc_capture
      logical init
c
      if (init) then
c
c     loop over all bonds, storing the atoms in each torsion
c
        if(deb_Path) print*,'torsions init'
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

#ifdef USE_NVSHMEM_CUDA
        ntors_pe = value_pe(ntors)
#endif

        call update_device_torsions
      end if

      call prmem_request(torsglob,6*nbloc,async=.true.)
!     Wait for nbondloc:
!     CPU scalar might be used in loop bounds
!$acc wait

!$acc data present(ntorsloc,nbondloc)
!$acc serial async
      ntorsloc = 0
!$acc end serial
!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc& default(present) deviceptr(d_nbtors,d_ibnd)
#else
!$acc& present(bndglob, nbtors, ibnd, n12, i12) 
#endif
      do i = 1, nbondloc
        ibond = bndglob(i)
#ifdef USE_NVSHMEM_CUDA
        ipe =     (ibond-1)/nbond_pe
        ind = mod((ibond-1),nbond_pe) +1
        torscount = d_nbtors(ipe)%pel(ind)
        ib  = d_ibnd(ipe)%pel(1,ind)
        ic  = d_ibnd(ipe)%pel(2,ind)
#else
        torscount = nbtors(ibond)
        ib  = ibnd(1,ibond)
        ic  = ibnd(2,ibond)
#endif
        ntorsloc1  = 0
!$acc loop seq
        do j = 1, n12(ib)
           ia = i12(j,ib)
           if (ia .ne. ic) then
!$acc loop seq
              do k = 1, n12(ic)
                 id = i12(k,ic)
                 if (id.ne.ib .and. id.ne.ia) then
!$acc atomic capture
                    ntorsloc = ntorsloc + 1
                    ntorsloc_capture = ntorsloc
!$acc end atomic
                    ntorsloc1 = ntorsloc1 + 1
                    torsglob(ntorsloc_capture) = torscount + ntorsloc1
                 end if
              end do
           end if
        end do
      end do
!$acc update host(ntorsloc) async
!$acc end data

      end subroutine

      subroutine update_device_torsions
      use bond
      use domdec,only: rank,hostcomm
      use mpi   ,only: MPI_BARRIER
      use nvshmem
      use sizes ,only:tinkerdebug
      use tors
      use tinMemory
      implicit none
      integer ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_torsions')
      if(rank.eq.0.and.tinkerdebug) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif

#ifdef USE_NVSHMEM_CUDA
      if (nbond.ne.0) then
      call shmem_update_device(nbtors,size(nbtors),
     &     dst=c_nbtors(mype)%pel,nd=nbond_pe,config=mnvshonly)
c     call nvshmem_check_data( nbtors,c_nbtors(mype)%pel,nbond,
c    &     nbond_pe,1)
      end if

      if (ntors.ne.0) then
      call shmem_update_device(itors,size(itors),
     &     dst=c_itors(mype)%pel,nd=4*ntors_pe,config=mnvshonly)
      end if
#else
!$acc update device(nbtors)
!$acc update device(itors)
#endif
!$acc enter data create(ntorsloc)
      end subroutine
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
      use mpi
      use pitors
      use strtor
      use tors
      use tinMemory
      implicit none
c
      !Exit condition
      !if (associated(kprop).and.ntors.eq.size(kprop)) return

#ifdef USE_NVSHMEM_CUDA
      call shmem_request(itors, winitors, [4,ntors],
     &     c_itors, d_itors, config=mhostnvsh)
      call shmem_request(tors1, wintors1, [4,ntors],
     &     c_tors1, d_tors1, config=mhostnvsh)
      call shmem_request(tors2, wintors2, [4,ntors],
     &     c_tors2, d_tors2, config=mhostnvsh)
      call shmem_request(tors3, wintors3, [4,ntors],
     &     c_tors3, d_tors3, config=mhostnvsh)
      call shmem_request(tors4, wintors4, [4,ntors],
     &     c_tors4, d_tors4, config=mhostnvsh)
      call shmem_request(tors5, wintors5, [4,ntors],
     &     c_tors5, d_tors5, config=mhostnvsh)
      call shmem_request(tors6, wintors6, [4,ntors],
     &     c_tors6, d_tors6, config=mhostnvsh)

      ! self association (OpenAcc visibility)
      d_itors => d_itors
      d_tors1 => d_tors1
      d_tors2 => d_tors2
      d_tors3 => d_tors3
      d_tors4 => d_tors4
      d_tors5 => d_tors5
      d_tors6 => d_tors6
#else
      ! tors parameters data
      call shmem_request(itors,winitors,[4,ntors],config=mhostacc)
      call shmem_request(tors1,wintors1,[4,ntors],config=mhostacc)
      call shmem_request(tors2,wintors2,[4,ntors],config=mhostacc)
      call shmem_request(tors3,wintors3,[4,ntors],config=mhostacc)
      call shmem_request(tors4,wintors4,[4,ntors],config=mhostacc)
      call shmem_request(tors5,wintors5,[4,ntors],config=mhostacc)
      call shmem_request(tors6,wintors6,[4,ntors],config=mhostacc)
#endif
      ! improp parameters data
      call shmem_request( kprop, winkprop,  [ntors],config=mhostacc)
      call shmem_request( vprop, winvprop,  [ntors],config=mhostacc)
      call shmem_request(iiprop,winiiprop,[4,ntors],config=mhostacc)

      call shmem_request(itors1,winitors1,[4,ntors],config=mhostacc)
      call shmem_request(itors2,winitors2,[4,ntors],config=mhostacc)
      call shmem_request(itors3,winitors3,[4,ntors],config=mhostacc)
      call shmem_request(iitors,winiitors,[4,ntors],config=mhostacc)
      !TODO 1.2 Warn louis about iiprop's double allocation in Tinker-HP 1.2
      call shmem_request(kpit,winkpit,  [ntors],config=mhostacc)
      call shmem_request(ipit,winipit,[6,ntors],config=mhostacc)

      ! kstrtor parameters data
      call shmem_request(ist,     winist,    [4,ntors],config=mhostacc)
      call shmem_request(kst,     winkst,    [9,ntors],config=mhostacc)
      call shmem_request(nbstrtor,winnbstrtor, [ntors],config=mhostacc)

      ! angtor parameters data
      call shmem_request( iat,       winiat, [4,ntors],config=mhostacc)
      call shmem_request(kant,      winkant, [6,ntors],config=mhostacc)
      call shmem_request(nbangtor,winnbangtor, [ntors],config=mhostacc)
      end
