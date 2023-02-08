c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kstrtor  --  find stretch-torsion parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kstrtor" assigns stretch-torsion parameters to torsions
c     needing them, and processes any new or changed values
c
c
#include "tinker_macro.h"
      module kstrtor_inl
        interface fron_conv_base_inl
          module procedure front_convert_base3
          module procedure front_convert_base5
        end interface
        contains
        include "fb_conv_base.f.inc"
      end module

      subroutine kstrtor(init)
      use atmlst
      use atmtyp
      use couple
      use domdec
      use inform
      use iounit
      use keys
      use ksttor
      use kstrtor_inl
      use potent
      use tinheader
      use tinMemory,only:prmem_request
#ifdef _OPENACC
      use thrust
      use utilgpu  ,only:rec_stream
#endif
      use utils    ,only:is_find8
      use sizes    ,only:tinkerdebug
      use strtor
      use tors
      implicit none
      integer i,j,k,nbt,nbt_sys
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer iitors,strtorcount,nstrtorloc1
      integer nstrtorloc_cap
      integer next,zero
      !integer size
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      integer :: isys=0
      integer*8 pt,pt0
      real(t_p) bt1,bt2,bt3
      real(t_p) bt4,bt5,bt6
      real(t_p) bt7,bt8,bt9
      logical header,swap
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*16 blank
      character*20 keyword
      character*240 record
      character*240 string
      parameter(zero=0)
      logical,intent(in):: init
c
      blank = '                '
      zeros = '0000'
      if (init) then
c
c     process keywords containing stretch-torsion parameters
c
        if (deb_Path) print*,'kstrtor init'
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:8) .eq. 'STRTORS ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              bt1 = 0.0_ti_p
              bt2 = 0.0_ti_p
              bt3 = 0.0_ti_p
              bt4 = 0.0_ti_p
              bt5 = 0.0_ti_p
              bt6 = 0.0_ti_p
              bt7 = 0.0_ti_p
              bt8 = 0.0_ti_p
              bt9 = 0.0_ti_p
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,id,bt1,bt2,bt3,
     &                                     bt4,bt5,bt6,bt7,bt8,bt9
   10         continue
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20              format (/,' Additional Stretch-Torsion Parameters :',
     &                    //,5x,'Atom Classes',10x,'Stretch',
     &                       9x,'1-Fold',6x,'2-Fold',6x,'3-Fold',/)
                 end if
                 if (rank.eq.0) write (iout,30) ia,ib,ic,id,bt1,bt2,bt3,
     &                          bt4,bt5,bt6,bt7,bt8,bt9
   30            format (2x,4i4,8x,'1st Bond',3x,3f12.3,
     &                 /,26x,'2nd Bond',3x,3f12.3,
     &                 /,26x,'3rd Bond',3x,3f12.3)
              end if
c             size = 4
c             call numeral (ita,pa,size)
c             call numeral (itb,pb,size)
c             call numeral (itc,pc,size)
c             call numeral (itd,pd,size)
              if (itb .lt. itc) then
c                pt = pa//pb//pc//pd
                 call fron_conv_base_inl(zero,ita,itb,itc,itd,pt)
              else if (itc .lt. itb) then
c                pt = pd//pc//pb//pa
                 call fron_conv_base_inl(zero,itd,itc,itb,ita,pt)
              else if (ita .le. itd) then
c                pt = pa//pb//pc//pd
                 call fron_conv_base_inl(zero,ita,itb,itc,itd,pt)
              else if (itd .lt. ita) then
c                pt = pd//pc//pb//pa
                 call fron_conv_base_inl(zero,itd,itc,itb,ita,pt)
              end if
              do j = 1, maxnbt
                 if (kbt(j).eq.-1 .or. kbt(j).eq.pt) then
                    kbt(j) = pt
                    btcon(4,j) = bt4
                    btcon(5,j) = bt5
                    btcon(6,j) = bt6
                    if (swap) then
                       btcon(1,j) = bt7
                       btcon(2,j) = bt8
                       btcon(3,j) = bt9
                       btcon(7,j) = bt1
                       btcon(8,j) = bt2
                       btcon(9,j) = bt3
                    else
                       btcon(1,j) = bt1
                       btcon(2,j) = bt2
                       btcon(3,j) = bt3
                       btcon(7,j) = bt7
                       btcon(8,j) = bt8
                       btcon(9,j) = bt9
                    end if
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KSTRTOR  --  Too many Stretch-Torsion',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        nbt = maxnbt
        do i = maxnbt, 1, -1
           if (kbt(i) .eq. -1)  nbt = i - 1
        end do
c
c       assign the stretch-torsion parameters for each torsion
c
        nstrtor = 0
        if (nbt .ne. 0) then
           do i = 1, ntors
              ia = itors(1,i)
              ib = itors(2,i)
              ic = itors(3,i)
              id = itors(4,i)
              nbstrtor(i) = nstrtor
              ita = class(ia)
              itb = class(ib)
              itc = class(ic)
              itd = class(id)
c             size = 4
c             call numeral (ita,pa,size)
c             call numeral (itb,pb,size)
c             call numeral (itc,pc,size)
c             call numeral (itd,pd,size)
              if (itb .lt. itc) then
c                pt = pa//pb//pc//pd
                 call fron_conv_base_inl(zero,ita,itb,itc,itd,pt)
                 call fron_conv_base_inl(zero,zero,itb,itc,zero,pt0)
                 swap = .false.
              else if (itc .lt. itb) then
c                pt = pd//pc//pb//pa
                 call fron_conv_base_inl(zero,itd,itc,itb,ita,pt)
                 call fron_conv_base_inl(zero,zero,itc,itb,zero,pt0)
                 swap = .true.
              else if (ita .le. itd) then
c                pt = pa//pb//pc//pd
                 call fron_conv_base_inl(zero,ita,itb,itc,itd,pt)
                 call fron_conv_base_inl(zero,zero,itb,itc,zero,pt0)
                 swap = .false.
              else if (itd .lt. ita) then
c                pt = pd//pc//pb//pa
                 call fron_conv_base_inl(zero,itd,itc,itb,ita,pt)
                 call fron_conv_base_inl(zero,zero,itc,itb,zero,pt0)
                 swap = .true.
              end if
              do j = 1, nbt
                 if (kbt(j) .eq. pt) then
                    nstrtor = nstrtor + 1
                    kst(4,nstrtor) = btcon(4,j)
                    kst(5,nstrtor) = btcon(5,j)
                    kst(6,nstrtor) = btcon(6,j)
                    if (swap) then
                       kst(1,nstrtor) = btcon(7,j)
                       kst(2,nstrtor) = btcon(8,j)
                       kst(3,nstrtor) = btcon(9,j)
                       kst(7,nstrtor) = btcon(1,j)
                       kst(8,nstrtor) = btcon(2,j)
                       kst(9,nstrtor) = btcon(3,j)
                    else
                       kst(1,nstrtor) = btcon(1,j)
                       kst(2,nstrtor) = btcon(2,j)
                       kst(3,nstrtor) = btcon(3,j)
                       kst(7,nstrtor) = btcon(7,j)
                       kst(8,nstrtor) = btcon(8,j)
                       kst(9,nstrtor) = btcon(9,j)
                    end if
                    ! construct the system class parameter
                    if (.not.is_find8(kbt_sys(1),isys,kbt(j))) then
                       isys = isys +1
                       kbt_sys(isys) = kbt(j)
                    end if
                    ist(1,nstrtor) = i
                    do k = 1, n12(ia)
                       if (i12(k,ia) .eq. ib) then
                          ist(2,nstrtor) = bndlist(k,ia)
                          goto 60
                       endif
                    end do
   60               continue
                    do k = 1, n12(ib)
                       if (i12(k,ib) .eq. ic) then
                          ist(3,nstrtor) = bndlist(k,ib)
                          goto 70
                       end if
                    end do
   70               continue
                    do k = 1, n12(ic)
                       if (i12(k,ic) .eq. id) then
                          ist(4,nstrtor) = bndlist(k,ic)
                          goto 100
                       end if
                    end do
                 end if
              end do
              do j = 1, nbt
                 if (kbt(j) .eq. pt0) then
                    nstrtor = nstrtor + 1
                    kst(4,nstrtor) = btcon(4,j)
                    kst(5,nstrtor) = btcon(5,j)
                    kst(6,nstrtor) = btcon(6,j)
                    if (swap) then
                       kst(1,nstrtor) = btcon(7,j)
                       kst(2,nstrtor) = btcon(8,j)
                       kst(3,nstrtor) = btcon(9,j)
                       kst(7,nstrtor) = btcon(1,j)
                       kst(8,nstrtor) = btcon(2,j)
                       kst(9,nstrtor) = btcon(3,j)
                    else
                       kst(1,nstrtor) = btcon(1,j)
                       kst(2,nstrtor) = btcon(2,j)
                       kst(3,nstrtor) = btcon(3,j)
                       kst(7,nstrtor) = btcon(7,j)
                       kst(8,nstrtor) = btcon(8,j)
                       kst(9,nstrtor) = btcon(9,j)
                    end if
                    ! construct the system class parameter
                    if (.not.is_find8(kbt_sys(1),isys,kbt(j))) then
                       isys = isys +1
                       kbt_sys(isys) = kbt(j)
                    end if
                    ist(1,nstrtor) = i
                    do k = 1, n12(ia)
                       if (i12(k,ia) .eq. ib) then
                          ist(2,nstrtor) = bndlist(k,ia)
                          goto 80
                       endif
                    end do
   80               continue
                    do k = 1, n12(ib)
                       if (i12(k,ib) .eq. ic) then
                          ist(3,nstrtor) = bndlist(k,ib)
                          goto 90
                       end if
                    end do
   90               continue
                    do k = 1, n12(ic)
                       if (i12(k,ic) .eq. id) then
                          ist(4,nstrtor) = bndlist(k,ic)
                          goto 100
                       end if
                    end do
                 end if
              end do
  100         continue
           end do
           kbt_sys(0) = isys
        end if
c
c       turn off the stretch-torsion potential if it is not used
c
        if (nstrtor .eq. 0)  use_strtor = .false.
c
c     Upload data Stretch Bend data on Device
c
        if (deb_Path) print*, 'nbt',nbt,int(kbt_sys(0),4)
        if (use_strtor) then
           !print*,nstrtor,isys
           call upload_device_kstrtor
        else
           call delete_data_kstrtor
           return
        end if
      end if
c     nbt = maxnbt
c     do i = maxnbt, 1, -1
c        if (kbt(i) .eq. -1)  nbt = i - 1
c     end do
      nbt = size_i8_to_i(kbt_sys(0))

      call prmem_request(strtorglob,ntorsloc,async=.true.)

!$acc serial present(nstrtorloc) async
      nstrtorloc = 0
!$acc end serial

      if (nbt .ne. 0) then
!$acc parallel loop 
#ifdef USE_NVSHMEM_CUDA
!$acc&         present(nstrtorloc,nstrtor,strtorglob,
#else
!$acc&         present(nstrtorloc,nstrtor,strtorglob,itors,
#endif
!$acc&   torsglob,nbstrtor) async
         do i = 1, ntorsloc
            iitors = torsglob(i)
#ifdef USE_NVSHMEM_CUDA
            ipe = (iitors-1)/ntors_pe
            ind = mod((iitors-1),ntors_pe) +1
            ia = d_itors(ipe)%pel(1,ind)
            ib = d_itors(ipe)%pel(2,ind)
            ic = d_itors(ipe)%pel(3,ind)
            id = d_itors(ipe)%pel(4,ind)
#else
            ia = itors(1,iitors)
            ib = itors(2,iitors)
            ic = itors(3,iitors)
            id = itors(4,iitors)
#endif
            strtorcount = nbstrtor(iitors)
            nstrtorloc1 = 0
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
c           size = 4
c           call numeral (ita,pa,size)
c           call numeral (itb,pb,size)
c           call numeral (itc,pc,size)
c           call numeral (itd,pd,size)
            if (itb .lt. itc) then
c              pt = pa//pb//pc//pd
               call fron_conv_base_inl(zero,ita,itb,itc,itd,pt)
               call fron_conv_base_inl(zero,zero,itb,itc,zero,pt0)
            else if (itc .lt. itb) then
c              pt = pd//pc//pb//pa
               call fron_conv_base_inl(zero,itd,itc,itb,ita,pt)
               call fron_conv_base_inl(zero,zero,itc,itb,zero,pt0)
            else if (ita .le. itd) then
c              pt = pa//pb//pc//pd
               call fron_conv_base_inl(zero,ita,itb,itc,itd,pt)
               call fron_conv_base_inl(zero,zero,itb,itc,zero,pt0)
            else if (itd .lt. ita) then
c              pt = pd//pc//pb//pa
               call fron_conv_base_inl(zero,itd,itc,itb,ita,pt)
               call fron_conv_base_inl(zero,zero,itc,itb,zero,pt0)
            end if
            do j = 1, nbt
               if (kbt_sys(j) .eq. pt .or. kbt_sys(j).eq. pt0) then
!$acc atomic capture
                  nstrtorloc = nstrtorloc + 1
                  nstrtorloc_cap = nstrtorloc
!$acc end atomic
                  nstrtorloc1 = nstrtorloc1 + 1
                  strtorglob(nstrtorloc_cap) = strtorcount + nstrtorloc1
                  exit
               end if
            end do
         end do
!$acc update host(nstrtorloc) async
      end if

      end

      subroutine upload_device_kstrtor
      use domdec,only: rank,hostcomm
      use ksttor
      use mpi   ,only: MPI_BARRIER
      use sizes ,only: tinkerdebug
      use strtor
      implicit none
      integer ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_kstrtor')
      if(rank.eq.0.and.tinkerdebug) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif

!$acc enter data copyin(nstrtorloc,nstrtor)
!$acc update device(ist,kst,nbstrtor)
!$acc update device(kbt,kbt_sys)
      end subroutine

      subroutine delete_data_kstrtor
      use domdec,only: rank
      use ksttor
      use sizes ,only: tinkerdebug
      use strtor
      use tinMemory
      implicit none

 12   format(2x,'delete_data_kstrtor')
      if(rank.eq.0.and.tinkerdebug) print 12

!$acc exit data delete(nstrtorloc,nstrtor)
      call shmem_request(ist,     winist,    [0,0],config=mhostacc)
      call shmem_request(kst,     winkst,    [0,0],config=mhostacc)
      call shmem_request(nbstrtor,winnbstrtor, [0],config=mhostacc)
      end subroutine
