!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine kangang  --  angle-angle parameter assignment  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "kangang" assigns the parameters for angle-angle cross term
!     interactions and processes new or changed parameter values
!
!
#include "tinker_macro.h"
subroutine kangang(init)
   use angang
   use angle
   use atmlst
   use atmtyp
   use atoms
   use couple
   use domdec
   use inform
   use iounit
   use kanang
   use keys
   use potent
   use tinheader
   use tinMemory
   use tors
   implicit none
   integer i,j,k,m,next
   integer it,ia,ic
   integer iglob
   integer nang,jang,kang,nangangloc1,angangcount
   real(t_p) fa,faa,aak(3)
   logical header
   character*20 keyword
   character*240 record
   character*240 string
   logical init
!
   if (init) then
!
!     process keywords containing angle-angle parameters
!
      if(deb_Path) print*,'kangang'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'ANGANG ') then
            it = 0
            do j = 1, 3
               aak(j) = 0.0_ti_p
            end do
            string = record(next:240)
            read (string,*,err=10,end=10)  it,(aak(j),j=1,3)
10          continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  if (rank.eq.0) write (iout,20)
20                format (/,' Additional Angle-Angle Parameters :',&
                     &//,5x,'Atom Class',8x,'K(AA) 1',4x,&
                     &'K(AA) 2',4x,'K(AA) 3',/)
               end if
               if (rank.eq.0) write (iout,30)  it,(aak(j),j=1,3)
30             format (9x,i3,7x,3f11.3)
            end if
            do j = 1, 3
               anan(j,it) = aak(j)
            end do
         end if
      end do
!
!       assign the angle-angle parameters for each angle pair
!
      nangang = 0
      do i = 1, n
         nang = n12(i) * (n12(i)-1) / 2
         it = class(i)
         nbangang(i) = nangang
         do j = 1, nang-1
            jang = anglist(j,i)
            ia = iang(1,jang)
            ic = iang(3,jang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            fa = anan(m,it)
            do k = j+1, nang
               kang = anglist(k,i)
               ia = iang(1,kang)
               ic = iang(3,kang)
               m = 1
               if (atomic(ia) .le. 1)  m = m + 1
               if (atomic(ic) .le. 1)  m = m + 1
               faa = fa * anan(m,it)
               if (faa .ne. 0.0_ti_p) then
                  nangang = nangang + 1
               end if
            end do
         end do
      end do
!
!       allocate global arrays
!
      call alloc_shared_angang
!
      nangangloc = 0
      do i = 1, n
         nang = n12(i) * (n12(i)-1) / 2
         it = class(i)
         do j = 1, nang-1
            jang = anglist(j,i)
            ia = iang(1,jang)
            ic = iang(3,jang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            fa = anan(m,it)
            do k = j+1, nang
               kang = anglist(k,i)
               ia = iang(1,kang)
               ic = iang(3,kang)
               m = 1
               if (atomic(ia) .le. 1)  m = m + 1
               if (atomic(ic) .le. 1)  m = m + 1
               faa = fa * anan(m,it)
               if (faa .ne. 0.0_ti_p) then
                  nangangloc = nangangloc + 1
                  iaa(1,nangangloc) = jang
                  iaa(2,nangangloc) = kang
                  kaa(nangangloc) = faa
               end if
            end do
         end do
      end do
!
!       turn off the angle-angle potential if it is not used
!
      if (nangang .eq. 0)  use_angang = .false.

      if (.not.use_angang) then
         call dealloc_shared_angang  ! Clean angang memory space
         return
      end if

   end if
   call prmem_request(angangglob,ntorsloc,config=mhostonly)
   nangangloc = 0
   do i = 1, nloc
      iglob = glob(i)
      nang = n12(iglob) * (n12(iglob)-1) / 2
      it = class(iglob)
      angangcount = nbangang(iglob)
      nangangloc1 = 0
      do j = 1, nang-1
         jang = anglist(j,iglob)
         ia = iang(1,jang)
         ic = iang(3,jang)
         m = 1
         if (atomic(ia) .le. 1)  m = m + 1
         if (atomic(ic) .le. 1)  m = m + 1
         fa = anan(m,it)
         do k = j+1, nang
            kang = anglist(k,iglob)
            ia = iang(1,kang)
            ic = iang(3,kang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            faa = fa * anan(m,it)
            if (faa .ne. 0.0_ti_p) then
               nangangloc = nangangloc + 1
               nangangloc1 = nangangloc1 + 1
               angangglob(nangangloc) = angangcount + nangangloc1
            end if
         end do
      end do
   end do
!
   return
end
!
!     subroutine dealloc_shared_angang : deallocate shared memory pointers for angang
!     parameter arrays
!
subroutine dealloc_shared_angang
   use sizes
   use angang
   use domdec
   use tinMemory
   implicit none
!
   call shmem_request(iaa,winiaa,[2,0])
   call shmem_request(kaa,winkaa,  [0])
end
!
!     subroutine alloc_shared_angang : allocate shared memory pointers for angang
!     parameter arrays
!
subroutine alloc_shared_angang
   use sizes
   use angang
   use domdec
   use tinMemory
   implicit none
!
   if (associated(kaa).and.nangang.eq.size(kaa)) return

   call shmem_request(iaa,winiaa,[2,nangang])
   call shmem_request(kaa,winkaa,  [nangang])
end
