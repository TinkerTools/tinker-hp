c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine kpitors  --  find pi-orbital torsion parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "kpitors" assigns pi-orbital torsion parameters to torsions
c     needing them, and processes any new or changed values
c
c
#include "tinker_precision.h"
      subroutine kpitors(init)
      use atmlst
      use atmtyp
      use bond
      use couple
      use domdec
      use inform
      use iounit
      use keys
      use kpitor
      use nvshmem
      use pitors
      use potent
      use tinheader
      use tors
      use utils
      use utilgpu
      implicit none
      integer i,j,npt
      integer ia,ib
      integer ibond,pitorscount,npitorsloc1
      integer ita,itb,temp
      integer ipe,ind
      integer size,next
      integer :: isys=0
      integer*8 pt
      real(t_p) tp
      logical header
      character*4 pa,pb
      character*8 blank
      character*20 keyword
      character*120 record
      character*120 string
      logical init
      integer npitorsloc_capture
c
c     blank = '        '
      if (init) then
c
c     process keywords containing pi-orbital torsion parameters
c
        if(rank.eq.0.and.tinkerdebug) print*,'kpitors init'
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'PITORS ') then
              ia = 0
              ib = 0
              tp = 0.0_ti_p
              string = record(next:120)
              read (string,*,err=10,end=10)  ia,ib,tp
   10         continue
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Pi-Orbital Torsion',
     &                         ' Parameters :',
     &                      //,5x,'Atom Classes',7x,'2-Fold',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,tp
   30            format (6x,2i4,4x,f12.3)
              end if
c             size = 4
c             call numeral (ia,pa,size)
c             call numeral (ib,pb,size)
              if (ia .le. ib) then
                 call front_convert_base(0,ia,ib,pt)
c                pt = pa//pb
              else
                 call front_convert_base(0,ib,ia,pt)
c                pt = pb//pa
              end if
              do j = 1, maxnpt
                 if (kpt(j).eq.-1 .or. kpt(j).eq.pt) then
                    kpt(j) = pt
                    ptcon(j) = tp
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KPITORS  --  Too many Pi-Orbital Torsion',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        npt = maxnpt
        do i = maxnpt, 1, -1
           if (kpt(i) .eq. -1)  npt = i - 1
        end do
c
c       assign pi-orbital torsion parameters as required
c
        npitors = 0
        if (npt .ne. 0) then
           do i = 1, nbond
              ia = ibnd(1,i)
              ib = ibnd(2,i)
              nbpitors(i) = npitors
              if (n12(ia).eq.3 .and. n12(ib).eq.3) then
                 ita = class(ia)
                 itb = class(ib)
c                size = 4
c                call numeral (ita,pa,size)
c                call numeral (itb,pb,size)
                 if (ita .le. itb) then
                    call front_convert_base(0,ita,itb,pt)
c                   pt = pa//pb
                 else
                    call front_convert_base(0,itb,ita,pt)
c                   pt = pb//pa
                 end if
                 do j = 1, npt
                    if (kpt(j) .eq. pt) then
                       npitors = npitors + 1
                       kpit(npitors) = ptcon(j)
                       ipit(1,npitors) = i12(1,ia)
                       ipit(2,npitors) = i12(2,ia)
                       ipit(3,npitors) = ia
                       ipit(4,npitors) = ib
                       ipit(5,npitors) = i12(1,ib)
                       ipit(6,npitors) = i12(2,ib)
                       if (i12(1,ia) .eq. ib)
     &                    ipit(1,npitors) = i12(3,ia)
                       if (i12(2,ia) .eq. ib)
     &                    ipit(2,npitors) = i12(3,ia)
                       if (i12(1,ib) .eq. ia)
     &                    ipit(5,npitors) = i12(3,ib)
                       if (i12(2,ib) .eq. ia)
     &                    ipit(6,npitors) = i12(3,ib)
                       if (.not.is_find8(kpt_sys(1),isys,pt)) then
                          isys  = isys + 1
                          kpt_sys(isys) = pt
                       end if
                       goto 60
                    end if
                 end do
              end if
   60         continue
           end do
           kpt_sys(0) = isys
           !print*,'npt  ', npt,isys
        end if
c
c       turn off the pi-orbital torsion potential if it is not used
c
        if (npitors .eq. 0) use_pitors = .false.

        if (use_pitors) then
           ! Send data to device
           call upload_device_kpitors
        else
           call delete_data_kpitors
           return
        end if

      end if
c
      npt = size_i8_to_i(kpt_sys(0))
c     npt = maxnpt
c     do i = maxnpt, 1, -1
c        if (kpt(i) .eq. -1)  npt = i - 1
c     end do
c
!Wait for ntorsloc
!GLITCHY : size is `ntorsloc`, but loop in kernel is `1, nbondloc` ?
!$acc wait
      call prmem_request(pitorsglob,ntorsloc,async=.true.)
      
!$acc data present(npitorsloc, nbondloc) async
!$acc serial async
      npitorsloc = 0
!$acc end serial
      if (npt .ne. 0) then
!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc& default(present) deviceptr(d_nbpitors)
#else
!$acc& present(bndglob,ibnd,nbpitors,n12,class,pitorsglob,
!$acc&     kpt_sys)
#endif
         do i = 1, nbondloc
            ibond = bndglob(i)
#ifdef USE_NVSHMEM_CUDA
            ipe  = (ibond-1)/nbond_pe
            ind  = mod((ibond-1),nbond_pe) +1
            ia   = d_ibnd(ipe)%pel(1,ind)
            ib   = d_ibnd(ipe)%pel(2,ind)
            pitorscount = d_nbpitors(ipe)%pel(ind)
#else
            ia   = ibnd(1,ibond)
            ib   = ibnd(2,ibond)
            pitorscount = nbpitors(ibond)
#endif
            npitorsloc1 = 0
            if (n12(ia).eq.3 .and. n12(ib).eq.3) then
               ita = class(ia)
               itb = class(ib)
               call front_convert_base(0,min(ita,itb),max(ita,itb),pt)
               do j = 1, npt
                  if (kpt_sys(j) .eq. pt) then
!$acc atomic capture
                     npitorsloc = npitorsloc + 1
                     npitorsloc_capture = npitorsloc
!$acc end atomic
                     npitorsloc1 = npitorsloc1 + 1
                     pitorsglob(npitorsloc_capture) = pitorscount 
     &                                              + npitorsloc1
                     exit
                  end if
               end do
            end if
         end do
      end if
!$acc update host(npitorsloc) async
!$acc end data
      end
c
c     Upload Pitors shared data parameters on Device
c
      subroutine upload_device_kpitors
      use bond
      use domdec,only: rank,hostcomm
      use kpitor
      use mpi   ,only: MPI_BARRIER
      use nvshmem
      use pitors
      use sizes
      use tors
      use tinMemory
      implicit none
      integer ierr
#ifdef _OPENACC
 12   format(2x,'upload_device_kpitors')
      if(rank.eq.0.and.tinkerdebug) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif

#ifdef USE_NVSHMEM_CUDA
      call shmem_update_device(nbpitors,size(nbpitors),
     &     dst=c_nbpitors(mype)%pel,nd=nbond_pe,config=mnvshonly)
#else
!$acc update device(nbpitors)
#endif
!$acc enter data copyin(npitorsloc)
!$acc update device(kpit,ipit)
!$acc update device(kpt,kpt_sys)
      end subroutine

      subroutine delete_data_kpitors
      use domdec,only: rank
      use pitors
      use sizes ,only: tinkerdebug
      use tinMemory
      use tors
      implicit none

 12   format(2x,'delete_data_kpitors')
      if(rank.eq.0.and.tinkerdebug) print 12

!$acc exit data delete(npitorsloc)
      call shmem_request(kpit,winkpit,[0]  ,config=mhostacc)
      call shmem_request(ipit,winipit,[0,0],config=mhostacc)
      end subroutine
