!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine kchgflx  --  charge flux parameter assignment  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "kchgflx" assigns bond stretch and angle bend charge flux
!     correction values and processes any new or changed values
!     for these parameters
!
!
#include "tinker_macro.h"
module kchgflx_inl
  contains
  subroutine host_distribute(nelem,off,n)
  use domdec
  use mpi
  implicit none
  integer nelem,off,n
  integer nhost,i,left
  call MPI_COMM_SIZE(hostcomm,nhost,i)
  if (nhost.eq.1) then
     n    = nelem
     off  = 0
  else
     n    = nelem / nhost
     off  = hostrank*n
     left = modulo(nelem,nhost)
     if (hostrank.eq.nhost-1) n = n+left
  end if
  !write(0,'(A,3I4,3I10)') '>>>', rank, hostrank, nhost, off, n, left
  end subroutine
#include "fb_conv_base.inc.f90"
end module

subroutine kchgflx
   use sizes
   use angle
   use atmlst
   use atoms
   use atmtyp
   use bond
   use cflux
   use couple
   use domdec
   use inform
   use iounit
   use kangs
   use kbonds
   use kcflux
   use kchgflx_inl
   use keys
   use mpi
   use potent
   use tinheader
   use usage
   implicit none
   integer i,j,k
   integer ia,ib,ic
   integer ita,itb,itc
   integer na,nb
   integer size,next
   integer ierr,ndis,off
   integer(8),parameter:: initi=-1
   real*8 fc,bd,cfb
   real*8 cfa1,cfa2
   real*8 cfb1,cfb2
   logical headerb
   logical headera
   character*4 pa,pb,pc
   character*8 blank8
   integer(8) pt,pt2,pt3
   character*8 ptab,ptbc
   character*12 blank12
   character*20 keyword
   character*240 record
   character*240 string

   if (deb_Path) print*, 'kchgflx'
!
!     allocate global pointers
!
   call alloc_shared_chgflx
!
!     process keywords containing charge flux parameters
!
   blank8 = '        '
   blank12 = '            '
   size = 4
   headerb = .true.
   headera = .true.
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      if (keyword(1:9) .eq. 'BNDCFLUX ') then
         ia = 0
         ib = 0
         cfb = 0.0d0
         string = record(next:240)
         read (string,*,err=10,end=10) ia,ib,cfb
10       continue
         if (headerb .and. .not.silent) then
            headerb = .false.
            write (iout,20)
20          format (/,' Additional Bond Charge Flux Parameters :',&
               &//,5x,'Atom Classes',19x,'K(CFB)',/)
         end if
         if (.not. silent) then
            write (iout,30)  ia,ib,cfb
30          format (6x,2i4,13x,f15.6)
         end if
         !call numeral (ia,pa,size)
         !call numeral (ib,pb,size)
         if (ia .le. ib) then
            call front_convert_base3(0,ia,ib,pt2)
         else
            call front_convert_base3(0,ib,ia,pt2)
         end if
         do j = 1, maxncfb
            if (kcfb(j).eq.initi .or. kcfb(j).eq.pt2) then
               kcfb(j) = pt2
               if (ia .lt. ib) then
                  cflb(j) = cfb
               else if (ib .lt. ia) then
                  cflb(j) = -cfb
               else
                  cflb(j) = 0.0d0
                  write (iout,40)
40                format (/,' KCHGFLX  --  Bond Charge Flux for',&
                     &' Identical Classes Set to Zero')
               end if
               goto 50
            end if
         end do
50       continue
      else if (keyword(1:9) .eq. 'ANGCFLUX ') then
         ia = 0
         ib = 0
         ic = 0
         cfa1 = 0.0d0
         cfa2 = 0.0d0
         cfb1 = 0.0d0
         cfb2 = 0.0d0
         string = record(next:240)
         read (string,*,err=60,end=60) ia,ib,ic,cfa1,cfa2,cfb1,cfb2
60       continue
         if (headera .and. .not.silent) then
            headera = .false.
            write (iout,70)
70          format (/,' Additional Angle Charge Flux Parameters :',&
               &//,5x,'Atom Classes',10x,'K(CFA1)',&
               &7x,'K(CFA2)',7x,'K(CFB1)',7x,'K(CFB2)',/)
         end if
         if (.not. silent) then
            write (iout,80)  ia,ib,ic,cfa1,cfa2,cfb1,cfb2
80          format (4x,3i4,4x,4f14.6)
         end if
         !call numeral (ia,pa,size)
         !call numeral (ib,pb,size)
         !call numeral (ic,pc,size)
         if (ia .le. ic) then
            call front_convert_base3(ia,ib,ic,pt3)
            !pt3 = pa//pb//pc
         else
            call front_convert_base3(ic,ib,ia,pt3)
            !pt3 = pc//pb//pa
         end if
         do j = 1, maxncfa
            if (kcfa(j).eq.initi .or. kcfa(j).eq.pt3) then
               kcfa(j) = pt3
               cfla(1,j) = cfa1
               cfla(2,j) = cfa2
               cflab(1,j) = cfb1
               cflab(2,j) = cfb2
               goto 90
            end if
         end do
90       continue
      end if
   end do
!
!     determine the total number of forcefield parameters
!
   nb = maxncfb
   do i = maxncfb, 1, -1
      if (kcfb(i) .eq. initi)  nb = i - 1
   end do
   na = maxncfa
   do i = maxncfa, 1, -1
      if (kcfa(i) .eq. initi)  na = i - 1
   end do
!
!     assign bond charge flux parameters for each bond
!
   if (hostrank.ne.0) goto 100
   nbflx = 0
   do i = 1, nbond
      ia = ibnd(1,i)
      ib = ibnd(2,i)
      ita = class(ia)
      itb = class(ib)
      size = 4
      !call numeral (ita,pa,size)
      !call numeral (itb,pb,size)
      if (ita .le. itb) then
         call front_convert_base3(0,ia,ib,pt2)
         !pt2 = pa//pb
      else
         call front_convert_base3(0,ib,ia,pt2)
         !pt2 = pb//pa
      end if
      bflx(i) = 0.0d0
      do j = 1, nb
         if (kcfb(j) .eq. pt2) then
            nbflx = nbflx + 1
            if (ita .le. itb) then
               bflx(i) = cflb(j)
            else
               bflx(i) = -cflb(j)
            end if
         end if
      end do
   end do
!
!    assign angle charge flux parameters for each angle
!
   naflx = 0
   do i = 1, nangle
      ia = iang(1,i)
      ib = iang(2,i)
      ic = iang(3,i)
      ita = class(ia)
      itb = class(ib)
      itc = class(ic)
      !call numeral (ita,pa,size)
      !call numeral (itb,pb,size)
      !call numeral (itc,pc,size)
      if (ita .le. itc) then
         call front_convert_base3(ia,ib,ic,pt3)
         !pt3 = pa//pb//pc
      else
         call front_convert_base3(ic,ib,ia,pt3)
         !pt3 = pc//pb//pa
      end if
      !if (ita .le. itb) then
      !   call front_convert_base3(0,ia,ib,ptab)
      !   !ptab = pa//pb
      !else
      !   call front_convert_base3(0,ib,ia,ptab)
      !   !ptab = pb//pa
      !end if
      !if (itb .le. itc) then
      !   call front_convert_base3(0,ib,ic,ptbc)
      !   !ptbc = pb//pc
      !else
      !   call front_convert_base3(0,ic,ib,ptbc)
      !   !ptbc = pc//pb
      !end if
      aflx(1,i)  = 0.0d0
      aflx(2,i)  = 0.0d0
      abflx(1,i) = 0.0d0
      abflx(2,i) = 0.0d0
      do j = 1, na
         if (kcfa(j) .eq. pt3) then
            naflx = naflx + 1
            aflx(1,i) = cfla(1,j)
            aflx(2,i) = cfla(2,j)
            abflx(1,i) = cflab(1,j)
            abflx(2,i) = cflab(2,j)
         end if
      end do
   end do
100 call MPI_BARRIER(hostcomm,ierr)
   call MPI_BCAST(nbflx,1,MPI_INT,0,hostcomm,ierr)
   call MPI_BCAST(naflx,1,MPI_INT,0,hostcomm,ierr)
!
!     turn off bond and angle charge flux term if not used
!
   if (nbflx.eq.0 .and. naflx.eq.0)  use_chgflx = .false.
   if (.not.use_charge .and. .not.use_mpole&
      &.and. .not.use_polar)  use_chgflx = .false.
   if (use_chgflx.and.use_pmecore) then
      if (ranktot.eq.0) then
         write(iout,*) "Charge Flux not compatible with PMECORES"
         __TINKER_FATAL__
      end if
   end if

   if (use_chgflx) then
      call upload_dev_shd_chgflx
   else
      call dealloc_shared_chgflx
   end if
end

subroutine upload_dev_shd_chgflx
   use cflux
   implicit none
!$acc update device(aflx,bflx,abflx) async
end subroutine
!
!     subroutine dealloc_shared_chgflx :
!     deallocate shared memory pointers for chgflx parameter arrays
!
subroutine dealloc_shared_chgflx
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use cflux
   use inform    ,only: deb_Path
   use tinMemory
!
12 format(2x,'dealloc_shared_chgflx')
   if (deb_Path) print 12

   if (.not.associated(aflx)) return
   call shmem_request(aflx, winaflx, [2,0],config=mhostacc)
   call shmem_request(abflx,winabflx,[2,0],config=mhostacc)
   call shmem_request(bflx, winbflx,   [0],config=mhostacc)
end
!
!     subroutine alloc_shared_chgflx :
!     allocate shared memory pointers for chgflx parameter arrays
!
subroutine alloc_shared_chgflx
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
   use angle
   use atoms
   use bond
   use cflux
   use domdec
   use mpi
   use sizes
   use tinMemory
   implicit none
!
   if (associated(aflx).and.size(aflx).eq.2*nangle) return

   call shmem_request(aflx, winaflx, [2,nangle],config=mhostacc)
   call shmem_request(abflx,winabflx,[2,nangle],config=mhostacc)
   call shmem_request(bflx, winbflx,    [nbond],config=mhostacc)
end
