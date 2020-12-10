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
#include "tinker_precision.h"
      module bound_internal
      contains
      !routine that says whether an interaction between two particules has 
      !to be computed within the current domain or not 
      !(dd half cell method, Newton's3rd law)
#include "image.f.inc"
      subroutine halfcell(xi,yi,zi,xj,yj,zj,docompute,use_bounds)
          implicit none
          real(t_p) xr,yr,zr
          real(t_p) xi,yi,zi
          real(t_p) xj,yj,zj
          logical,intent(in):: use_bounds
          logical docompute
!$acc routine seq
          xr = xi - xj
          yr = yi - yj
          zr = zi - zj
          if (use_bounds) call image_inl(xr,yr,zr)
          docompute = (xr > 0)
     &               .or. (xr == 0 .and. yr > 0)
     &               .or. (xr == 0 .and. yr == 0 .and. zr > 0)
      end
      end
      
      subroutine bonds(init)
      use atmlst
      use atoms
      use bond
      use bound
      use couple
      use domdec
      use iounit
      use nvshmem
      use utilgpu
      use tinMemory
      use bound_internal
      implicit none
      integer i,j,k,m
      integer iglob,kglob
      integer ipe,ind,d_size,iibndl
      real(t_p) xi,yi,zi,xk,yk,zk
      logical init,docompute
      integer nbondloc_capture
   10              format (/,' BONDS  --  Too many Bonds; Increase',
     &                        ' MAXBND')
c
      if (init) then
c
c
c     loop over all atoms, storing the atoms in each bond
c
        if(rank.eq.0.and.tinkerdebug) print*,'bonds init'
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
c       allocate arrays
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
        call MPI_BARRIER(hostcomm,i)

        call upload_device_bond
      end if

      call prmem_request(bndglob,maxvalue*nbloc,async=.true.)

!$acc data present(nbondloc)
!$acc serial async
      nbondloc = 0
!$acc end serial
!$acc parallel loop default(none) async
#ifdef USE_NVSHMEM_CUDA
!$acc& present(bndglob,glob,x,y,z,n12,i12) deviceptr(d_bndlist)
#else
!$acc& present(bndglob, glob, x, y, z, n12, i12, bndlist)
#endif
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
            call halfcell(xi,yi,zi,xk,yk,zk,docompute,use_bounds)
            if (docompute) then
!$acc atomic capture
              nbondloc = nbondloc + 1
              nbondloc_capture = nbondloc
!$acc end atomic
#ifdef USE_NVSHMEM_CUDA
              ipe    =     (iglob-1)/n_pe
              ind    = mod((iglob-1),n_pe) +1
              iibndl = d_bndlist(ipe)%pel(j,ind)
#else
              iibndl = bndlist(j,iglob)
#endif
              bndglob(nbondloc_capture) = iibndl
            end if
         end do
      end do
!$acc update host(nbondloc) async
!$acc end data
      end

      subroutine upload_device_bond
      use atmlst
      use atoms
      use bond
      use domdec,only:rank
      use nvshmem
      use sizes,only:tinkerdebug
      use tinMemory
      implicit none
      integer d_size

#ifdef _OPENACC
 12   format(2x,'upload_device_bond')
      if(rank.eq.0.and.tinkerdebug) print 12
#endif

#ifdef USE_NVSHMEM_CUDA
      nbond_pe = value_pe(nbond)
      if (nbond.ne.0) then
      d_size   = size(c_ibnd(mype)%pel)
      !print*,'bonds.f i',rank,d_size,size(ibnd)
      call shmem_update_device(ibnd,size(ibnd),
     &     dst=c_ibnd(mype)%pel,nd=d_size,config=mhostnvsh)
      end if

      if (n.ne.0) then
      d_size = size(c_bndlist(mype)%pel)
      call shmem_update_device(bndlist,size(bndlist),
     &     dst=c_bndlist(mype)%pel,nd=d_size,config=mhostnvsh)
      end if
#else
!$acc update device(ibnd(:,:),bndlist(:,:))
#endif
!$acc enter data copyin(nbondloc)

      end subroutine

      subroutine delete_data_bonds
      use atmlst
      use bond
      use domdec
      use sizes
      implicit none
      end subroutine
c
c     subroutine alloc_shared_bond : allocate shared memory pointers for bond
c     parameter arrays
c
      subroutine alloc_shared_bond
      use atmlst
      use atoms ,only:n
      use domdec
      use bond
      use pitors
      use tors
      use tinMemory
      use utils
      implicit none
 
c     if (associated(bl)     .and.nbond.eq.size(bl) .and.
c    &    associated(bndlist).and.  8*n.eq.size(bndlist)) return

#ifdef USE_NVSHMEM_CUDA
      call shmem_request(ibnd    ,winibnd    ,[2,nbond], c_ibnd, d_ibnd,
     &     mhostnvsh)
      call shmem_request(bndlist ,winbndlist ,[8,n]    , c_bndlist,
     &     d_bndlist, mhostnvsh)
      call shmem_request(bk      ,winbk      ,[nbond]  , c_bk, d_bk,
     &     mhostnvsh)
      call shmem_request(bl      ,winbl      ,[nbond]  , c_bl, d_bl,
     &     mhostnvsh)
      call shmem_request(nbtors, winnbtors, [nbond], c_nbtors,
     &     d_nbtors, config=mhostnvsh)
      call shmem_request(nbpitors,winnbpitors,[nbond] , c_nbpitors,
     &     d_nbpitors,  config=mhostnvsh)
#else
      call shmem_request(bk      ,winbk      ,[nbond]  ,config=mhostacc)
      call shmem_request(bl      ,winbl      ,[nbond]  ,config=mhostacc)
      call shmem_request(bndlist ,winbndlist ,[8,n]    ,config=mhostacc)
      call shmem_request(ibnd    ,winibnd    ,[2,nbond],config=mhostacc)
      call shmem_request(nbtors  ,winnbtors  ,[nbond]  ,config=mhostacc)
      call shmem_request(nbpitors,winnbpitors,[nbond]  ,config=mhostacc)
#endif

#ifdef USE_NVSHMEM_CUDA
      ! Explicit Reassociate to make device pointer data accessible by OpenACC
      i2dDPC_expl => d_ibnd
      d_ibnd      => i2dDPC_expl
      i2dDPC_expl => d_bndlist
      d_bndlist   => i2dDPC_expl
      rDPC_expl   => d_bl
      d_bl        => rDPC_expl
      rDPC_expl   => d_bk
      d_bk        => rDPC_expl
      ! self association for openAcc visibility
      d_nbtors    => d_nbtors
      d_nbpitors  => d_nbpitors
#endif
      end
