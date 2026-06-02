!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine krepel  --  Pauli repulsion term assignment  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "krepel" assigns the size values, exponential parameter and
!     number of valence electrons for Pauli repulsion interactions
!     and processes any new or changed values for these parameters
!
!
#include "tinker_macro.h"
subroutine krepel
   use atoms
   use atmtyp
   use domdec ,only: hostrank,hostcomm
   use inform
   use iounit
   use krepl
   use keys
   use mpi
   use potent
   use repel
   use sizes
   use tinheader
   implicit none
   integer i,k,ierr
   integer ia,ic,next
   real*8 spr,apr,epr
   logical header
   character*20 keyword
   character*240 record
   character*240 string
!
   if (deb_Path) print*,'krepel'
!
!     process keywords containing Pauli repulsion parameters
!
   header = .true.
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      if (keyword(1:10) .eq. 'REPULSION ') then
         k = 0
         spr = 0.0_ti_p
         apr = 0.0_ti_p
         epr = 0.0_ti_p
         call getnumb (record,k,next)
         string = record(next:240)
         read (string,*,err=10,end=10)  spr,apr,epr
10       continue
         if (k .gt. 0) then
            if (header .and. .not.silent) then
               header = .false.
               write (iout,20)
20             format (/,' Additional Pauli Repulsion',&
                  &' Parameters :',&
                  &//,5x,'Atom Class',15x,'Size',11x,'Damp',&
                  &8x,'Valence'/)
            end if
            if (k .le. maxclass) then
               prsiz(k) = spr
               prdmp(k) = apr
               prele(k) = -abs(epr)
               if (.not. silent) then
                  write (iout,30)  k,spr,apr,epr
30                format (6x,i6,7x,2f15.4,f15.3)
               end if
            else
               write (iout,40)
40             format (/,' KREPEL  --  Too many Pauli Repulsion',&
                  &' Parameters')
               abort = .true.
            end if
         end if
      end if
   end do
!
!     allocate global pointers
!
   call alloc_shared_rep
   if (hostrank.ne.0) goto 100
!
!     assign the repulsion size, alpha and valence parameters
!
   do i = 1, n
      ic = class(i)
      sizpr(i) = prsiz(ic)
      dmppr(i) = prdmp(ic)
      elepr(i) = prele(ic)
   end do
!
!     process keywords containing atom specific Pauli repulsion
!
   header = .true.
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      if (keyword(1:10) .eq. 'REPULSION ') then
         ia = 0
         spr = 0.0_ti_p
         apr = 0.0_ti_p
         epr = 0.0_ti_p
         string = record(next:240)
         read (string,*,err=70,end=70)  ia,spr,apr,epr
         if (ia.lt.0 .and. ia.ge.-n) then
            ia = -ia
            if (header .and. .not.silent) then
               header = .false.
               write (iout,50)
50             format (/,' Additional Pauli Repulsion Values',&
                  &' for Specific Atoms :',&
                  &//,8x,'Atom',17x,'Size',12x,'Damp',&
                  &8x,'Valence'/)
            end if
            if (.not. silent) then
               write (iout,60)  ia,spr,apr,epr
60             format (6x,i6,7x,2f15.4,f15.3)
            end if
            sizpr(ia) = spr
            dmppr(ia) = apr
            elepr(ia) = -abs(epr)
         end if
70       continue
      end if
   end do

100 call MPI_BARRIER(hostcomm,ierr)
!
!     remove zero and undefined repulsion sites from the list
!
   nrep = 0
   do i = 1, n
      if (sizpr(i) .ne. 0.0_ti_p) then
         nrep = nrep + 1
      end if
   end do
!
!     turn off the Pauli repulsion potential if not used
!
   if (nrep .eq. 0)  then
      use_repuls = .false.
      call dealloc_shared_rep
   else
      call upload_dev_shr_repel
   end if
end

subroutine upload_dev_shr_repel
   use repel
   implicit none
!$acc update device(sizpr,dmppr,elepr) async
end subroutine
!
!     subroutine dealloc_shared_rep : deallocate shared memory pointers for repulsion
!     parameter arrays
!
subroutine dealloc_shared_rep
   use inform
   use repel
   use tinMemory
   implicit none
!
   if(.not.associated(sizpr)) return
   if(deb_Path) print*, '  dealloc_shared_rep'
   call shmem_request(sizpr, winsizpr, [0] ,config=mhostacc)
   call shmem_request(dmppr, windmppr, [0] ,config=mhostacc)
   call shmem_request(elepr, winelepr, [0] ,config=mhostacc)
end
!
!     subroutine alloc_shared_rep : allocate shared memory pointers for repulsion
!     parameter arrays
!
subroutine alloc_shared_rep
   use sizes
   use atoms
   use repel
   use tinMemory
   implicit none
!
   if (associated(sizpr).and.size(sizpr).eq.n) return
   call shmem_request(sizpr, winsizpr, [n] ,config=mhostacc)
   call shmem_request(dmppr, windmppr, [n] ,config=mhostacc)
   call shmem_request(elepr, winelepr, [n] ,config=mhostacc)
end
