!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine kchgtrn  --  charge transfer term assignment  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "kchgtrn" assigns charge magnitude and damping parameters for
!     charge transfer interactions and processes any new or changed
!     values for these parameters
!
!
#include "tinker_macro.h"
subroutine kchgtrn(init,istep)
   use atmlst
   use atoms
   use atmtyp
   use chgpen
   use chgtrn
   use cutoff
   use domdec
   use inform
   use iounit
   use kctrn
   use keys
   use mplpot
   use mpole
   use neigh
   use pme
   use polar
   use polpot
   use potent
   use sizes
   use mpi
   implicit none
   integer istep,modnl
   integer i,k,ierr
   integer iproc,iglob,polecount,iipole
   integer ia,ic,next
   real*8 d
   real*8 chtrn,actrn
   logical header
   character*20 keyword
   character*240 record
   character*240 string
   logical init
!
   if (init) then

      if (deb_Path) print*, 'kchgtrn init'
!
!     allocate global pointers
!
      call alloc_shared_chgct

      if (hostrank.ne.0) goto 100
!
!     process keywords containing charge transfer parameters
!
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHGTRN ') then
            k = 0
            chtrn = 0.0d0
            actrn = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  chtrn,actrn
10          continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
20                format (/,' Additional Charge Transfer',&
                     &' Parameters :',&
                     &//,5x,'Atom Class',13x,'Charge',11x,'Damp',/)
               end if
               if (k .le. maxclass) then
                  ctchg(k) = chtrn
                  ctdmp(k) = actrn
                  if (.not. silent) then
                     write (iout,30)  k,chtrn,actrn
30                   format (6x,i6,7x,f15.4,f15.4)
                  end if
               else
                  write (iout,40)
40                format (/,' KCHGTRN  --  Too many Charge',&
                     &' Transfer Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
!
!
!     assign the charge transfer charge and alpha parameters
!
      do i = 1, n
         ic = class(i)
         chgct(i) = ctchg(ic)
         dmpct(i) = ctdmp(ic)
      end do
!
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHGTRN ') then
            ia = 0
            chtrn = 0.0d0
            actrn = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,chtrn,actrn
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
50                format (/,' Additional Charge Transfer Values',&
                     &' for Specific Atoms :',&
                     &//,8x,'Atom',16x,'Charge',11x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,chtrn,actrn
60                format (6x,i6,7x,f15.4,f15.4)
               end if
               chgct(ia) = chtrn
               dmpct(ia) = actrn
            end if
70          continue
         end if
      end do
!
!     remove zero or undefined electrostatic sites from the list
!
      if (use_chgtrn) then
         npole = 0
         ncp = 0
         npolar = 0
         nct = 0
         do i = 1, n
            if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0 .or.&
               &chgct(i).ne. 0.0d0 .or. dmpct(i).ne.0.0d0) then
               nbpole(i) = npole
               npole = npole + 1
               ipole(npole) = i
               pollist(i) = npole
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               yaxis(npole) = yaxis(i)
               polaxe(npole) = polaxe(i)
               do k = 1, maxpole
                  pole(k,npole) = pole(k,i)
               end do
               mono0(npole) = pole(1,i)
               if (palpha(i) .ne. 0.0d0)  ncp = ncp + 1
               pcore(npole) = pcore(i)
               pval(npole) = pval(i)
               pval0(npole) = pval(i)
               palpha(npole) = palpha(i)
               if (polarity(i) .ne. 0.0d0) then
                  npolar = npolar + 1
!                   douind(i) = .true.
               end if
               if (tholed(i) .ne. 0.0d0)  use_tholed = .true.
               polarity(npole) = polarity(i)
               thole(npole) = thole(i)
               tholed(npole) = tholed(i)
               pdamp(npole) = pdamp(i)
               if (chgct(i).ne.0.0d0 .or. dmpct(i).ne.0.0d0) then
                  nct = nct + 1
               end if
               chgct(npole) = chgct(i)
               dmpct(npole) = dmpct(i)
            end if
         end do
         call FromPolaxe2Ipolaxe ! Brodcast 'nZ_onlyglob' after this call
      end if
100   call MPI_BARRIER(hostcomm,ierr)
      call MPI_BCAST(npole ,1,MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(npolar,1,MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(nct   ,1,MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(ncp   ,1,MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(nZ_Onlyglob,1,    MPI_INT,0,hostcomm,ierr)
      call MPI_BCAST(use_tholed,1,MPI_LOGICAL,0,hostcomm,ierr)
!
!     test multipoles at chiral sites and invert if necessary
!
      if (use_chgtrn) call chkpole(.true.)
!
!     turn off individual electrostatic potentials if not used
!
      if (npole  .eq.0) use_mpole  = .false.
      if (ncp    .ne.0) use_chgpen = .true.
      if (ncp    .ne.0) use_thole  = .false.
      if (npolar .eq.0) use_polar  = .false.
      if (nct    .eq.0) use_chgtrn = .false.
      if (use_tholed ) use_thole  = .true.

      call upload_dev_shr_chgct
      if (use_chgtrn) then
         use_mlist = .true.
      else
         call dealloc_shared_chgct
         return
      end if

   end if
end

subroutine upload_dev_shr_chgct
   use chgpen
   use chgtrn
   use mpole
   use polar
   implicit none

!$acc update device(chgct,dmpct)
!$acc update device(nbpole,ipole,pollist,xaxis,yaxis,zaxis,ipolaxe &
!$acc      ,pole,pcore,pval,pval0,palpha,polarity,thole,tholed,pdamp)

end subroutine
!
!     subroutine dealloc_shared_chgct : deallocate shared memory pointers for charge transfer
!     parameter arrays
!
subroutine dealloc_shared_chgct
   use chgpen
   use chgtrn
   use inform
   use tinMemory
   implicit none

   if (deb_Path) print*, '  dealloc_shared_chgct'
   call shmem_request(chgct, winchgct,  [0], config=mhostacc)
   call shmem_request(dmpct, windmpct,  [0], config=mhostacc)
end
!
!     subroutine alloc_shared_chgct : allocate shared memory pointers for charge transfer
!     parameter arrays
!
subroutine alloc_shared_chgct
   use atoms
   use chgpen
   use chgtrn
   use tinMemory
   implicit none

   if (associated(chgct).and.size(chgct).eq.n) return
   call shmem_request(chgct, winchgct,  [n], config=mhostacc)
   call shmem_request(dmpct, windmpct,  [n], config=mhostacc)
end
