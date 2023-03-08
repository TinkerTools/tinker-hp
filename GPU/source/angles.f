c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine angles  --  locate and store bond angles  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "angles" finds the total number of bond angles and stores
c     the atom numbers of the atoms defining each angle; for
c     each angle to a trivalent central atom, the third bonded
c     atom is stored for use in out-of-plane bending
c
c
#include "tinker_precision.h"
      subroutine angles(init)
      use angle
      use atmlst
      use atoms
      use couple
      use domdec
      use iounit
      use inform  ,only: deb_Path
      use nvshmem
#ifdef _OPENACC
      use thrust
#endif
      use utilgpu
      implicit none
      integer i,j,k,m,iglob
      integer ia,ib,ic
      integer ipe,ind,ianglst
      logical init
      integer nangleloc_capture

      if (init) then
c
c     loop over all atoms, storing the atoms in each bond angle
c
        if(deb_Path) print*,'angles init'
        nangle = 0
        do i = 1, n
           do j = 1, n12(i)-1
              do k = j+1, n12(i)
                nangle = nangle + 1
                if (nangle .gt. 4*n) then
                   if (rank.eq.0) write (iout,10)
   10              format (/,' ANGLES  --  Too many Bond Angles;',
     &                        ' Increase MAXANG')
                   call fatal
                end if
              end do
           end do
        end do
c
c       allocate arrays
c
        call alloc_shared_angles
c
        nangleloc = 0
        do i = 1, n
           m = 0
           do j = 1, n12(i)-1
              do k = j+1, n12(i)
                nangleloc = nangleloc + 1
                m = m + 1
                anglist(m,i) = nangleloc
                iang(1,nangleloc) = i12(j,i)
                iang(2,nangleloc) = i
                iang(3,nangleloc) = i12(k,i)
                iang(4,nangleloc) = 0
              end do
           end do
c
c       set the out-of-plane atom for angles at trivalent centers
c
           if (n12(i) .eq. 3) then
              iang(4,nangleloc) = i12(1,i)
              iang(4,nangleloc-1) = i12(2,i)
              iang(4,nangleloc-2) = i12(3,i)
           end if
        end do
c
c       store the numbers of the bonds comprising each bond angle
c
        do i = 1, nangle
           ia = iang(1,i)
           ib = iang(2,i)
           ic = iang(3,i)
           do k = 1, n12(ib)
              if (i12(k,ib) .eq. ia)  balist(1,i) = bndlist(k,ib)
              if (i12(k,ib) .eq. ic)  balist(2,i) = bndlist(k,ib)
           end do
        end do

        call MPI_BARRIER(hostcomm,i)
        nangle_pe=value_pe(nangle)
        call upload_device_angles
        call prmem_request(angleloc,nangle)
      end if

      call prmem_request(angleglob,4*nbloc,async=.true.)

!$acc data present(nangleloc)
!$acc serial async
      nangleloc = 0
!$acc end serial
!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc& default(present) deviceptr(d_anglist)
#else
!$acc& present(anglist,angleloc,angleglob, glob, n12)
#endif
      do i = 1, nloc
        m = 0
        iglob = glob(i)
!$acc loop seq
        do j = 1, n12(iglob)-1
!$acc loop seq
          do k = j+1, n12(iglob)
!$acc atomic capture
            nangleloc = nangleloc + 1
            nangleloc_capture = nangleloc
!$acc end atomic
            m = m + 1
#ifdef USE_NVSHMEM_CUDA
            ipe     =     (iglob-1)/n_pe
            ind     = mod((iglob-1),n_pe) +1
            ianglst = d_anglist(ipe)%pel(m,ind)
            angleglob(nangleloc_capture) = ianglst
            angleloc (ianglst)           = nangleloc_capture
#else
            angleglob(nangleloc_capture) = anglist(m,iglob)
            angleloc (anglist(m,iglob))  = nangleloc_capture
#endif
          end do
        end do
      end do
!$acc update host(nangleloc) async
!$acc end data

   20 continue
      end

      subroutine upload_device_angles
      use angle
      use atmlst
      use atoms
      use domdec
      use nvshmem
      use sizes
      use tinMemory
      implicit none
      integer c_size

#ifdef _OPENACC
 12   format(2x,'upload_device_angles')
      if(rank.eq.0.and.tinkerdebug) print 12
#endif
#ifdef USE_NVSHMEM_CUDA
      if (n.ne.0) then
      c_size = size(c_anglist(mype)%pel)
      call shmem_update_device(anglist,size(anglist),
     &     dst=c_anglist(mype)%pel,nd=c_size,config=mnvshonly)
      end if

      if (nangle.ne.0) then
      c_size = size(c_iang(mype)%pel)
      call shmem_update_device(iang,size(iang),
     &     dst=c_iang(mype)%pel,nd=c_size,config=mnvshonly)
      end if
#else
!$acc update device(anglist,balist,iang)
#endif
!$acc enter data copyin(nangleloc)
      end subroutine
c
c
c     subroutine alloc_shared_angles : allocate shared memory pointers for angles
c     parameter arrays
c
      subroutine alloc_shared_angles
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use angang
      use angle
      use angpot
      use atmlst
      use atoms
      use bitor
      use domdec
      use improp
      use imptor
      use opbend
      use opdist
      use strbnd
      use tinMemory
      use urey
      use mpi
      implicit none
c
c     if (associated(ak).and.size(ak).eq.nangle
c    &   .and.16*n.eq.size(anglist)) return ! Exit condition

#ifdef USE_NVSHMEM_CUDA
      call shmem_request(anglist, winanglist, [16,n],
     &     c_anglist, d_anglist, config=mhostnvsh)
      call shmem_request(iang,    winiang,    [4,nangle],
     &     c_iang, d_iang, config=mhostnvsh)

      ! self association (OpenAcc visibility)
      d_anglist => d_anglist
      d_iang    => d_iang
#else
      call shmem_request(anglist, winanglist,[16,n],   config=mhostacc)
      call shmem_request(iang,    winiang, [4,nangle], config=mhostacc)
#endif
      call shmem_request(balist,  winbalist, [2,6*n],  config=mhostacc)
      call shmem_request(ak,      winak,     [nangle], config=mhostacc)
      call shmem_request(anat,    winanat,   [nangle], config=mhostacc)
      call shmem_request(afld,    winafld,   [nangle], config=mhostacc)
      call shmem_request(angtyp,  winangtyp, [nangle], config=mhostonly)
      call shmem_request(angtypI, winangtypI,[nangle], config=mhostacc)

      call shmem_request(isb,     winisb,  [3,nangle], config=mhostacc)
      call shmem_request(sbk,     winsbk,  [2,nangle], config=mhostacc)
      call shmem_request(nbstrbnd,winnbstrbnd,[nangle],config=mhostacc)

      call shmem_request(uk,      winuk,     [nangle], config=mhostacc)
      call shmem_request(ul,      winul,     [nangle], config=mhostacc)
      call shmem_request(iury,    winiury, [3,nangle], config=mhostacc)
      call shmem_request(nburey,  winnburey,  [nangle],config=mhostacc)

      call shmem_request(nbbitors,winnbbitors,[nangle],config=mhostacc)
      call shmem_request(nbangang,winnbangang,     [n],config=mhostonly)
      call shmem_request(nbopbend,winnbopbend,[nangle],config=mhostacc)
      call shmem_request(nbopdist,winnbopdist,     [n],config=mhostonly)
      call shmem_request(nbimprop,winnbimprop,     [n],config=mhostacc)
      call shmem_request(nbimptor,winnbimptor,     [n],config=mhostacc)
      end
c
c
c     subroutine dealloc_shared_angles : deallocate shared memory pointers for angles
c     parameter arrays
c
      subroutine dealloc_shared_angles
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use angang
      use angle
      use angpot
      use atmlst
      use atoms
      use bitor
      use domdec
      use improp
      use imptor
      use opbend
      use opdist
      use strbnd
      use tinMemory
      use urey
      use mpi
      implicit none
c
c     if (associated(ak).and.size(ak).eq.nangle
c    &   .and.16*n.eq.size(anglist)) return ! Exit condition

#ifdef USE_NVSHMEM_CUDA
      call shmem_request(anglist, winanglist, [0,0],
     &     c_anglist, d_anglist, config=mhostnvsh)
      call shmem_request(iang,    winiang,    [0,0],
     &     c_iang, d_iang, config=mhostnvsh)

      ! self association (OpenAcc visibility)
      d_anglist => d_anglist
      d_iang    => d_iang
#else
      call shmem_request(anglist, winanglist,[0,0],   config=mhostacc)
      call shmem_request(iang,    winiang, [0,0], config=mhostacc)
#endif
      call shmem_request(balist,  winbalist, [0,0*n], config=mhostacc)
      call shmem_request(ak,      winak,     [0], config=mhostacc)
      call shmem_request(anat,    winanat,   [0], config=mhostacc)
      call shmem_request(afld,    winafld,   [0], config=mhostacc)
      call shmem_request(angtyp,  winangtyp, [0], config=mhostonly)
      call shmem_request(angtypI, winangtypI,[0], config=mhostacc)

      call shmem_request(isb,     winisb,  [0,0], config=mhostacc)
      call shmem_request(sbk,     winsbk,  [0,0], config=mhostacc)
      call shmem_request(nbstrbnd,winnbstrbnd,[0],config=mhostacc)

      call shmem_request(uk,      winuk,     [0], config=mhostacc)
      call shmem_request(ul,      winul,     [0], config=mhostacc)
      call shmem_request(iury,    winiury, [0,0], config=mhostacc)
      call shmem_request(nburey,  winnburey, [0], config=mhostacc)

      call shmem_request(nbbitors,winnbbitors,[0],config=mhostacc)
      call shmem_request(nbangang,winnbangang,[0],config=mhostonly)
      call shmem_request(nbopbend,winnbopbend,[0],config=mhostacc)
      call shmem_request(nbopdist,winnbopdist,[0],config=mhostonly)
      call shmem_request(nbimprop,winnbimprop,[0],config=mhostacc)
      call shmem_request(nbimptor,winnbimptor,[0],config=mhostacc)
      end
