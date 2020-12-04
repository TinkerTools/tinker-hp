c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine katom  --  atom type parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "katom" assigns an atom type definitions to each atom in
c     the structure and processes any new or changed values
c
c     literature reference:
c
c     K. A. Feenstra, B. Hess and H. J. C. Berendsen, "Improving
c     Efficiency of Large Time-Scale Molecular Dynamics Simulations
c     of Hydrogen-Rich Systems", Journal of Computational Chemistry,
c     8, 786-798 (1999)  [heavy hydrogen reweighting]
c
#include "tinker_precision.h"
      subroutine katom
      use atmtyp
      use atoms
      use couple
      use domdec
      use keys
      use inform
      use iounit
      use katoms
      use mpi,only:mpi_barrier
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,j,k,next,ierr
      integer nh
      integer cls,atn,lig
      real(t_p) wght
      real(r_p) hmax,hmass
      real(r_p) sum,dmin,dmass
      logical header,heavy
      character*3 symb
      character*20 keyword
      character*24 notice
      character*120 record
      character*120 string
c
c     allocate global arrays
c
      if (deb_Path) print*,'katoms'
      call alloc_shared_katom
c
c     process keywords containing atom type parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            cls = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0_ti_p
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            if (cls .eq. 0)  cls = k
            atmcls(k) = cls
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:120)
            read (string,*,err=40,end=40)  atn,wght,lig
            if (k.ge.1 .and. k.le.maxtyp) then
               if (header .and. .not.silent) then
                  header = .false.
                  if (rank.eq.0) write (iout,10)
   10             format (/,' Additional Atom Type Parameters :',
     &                    //,5x,'Type  Class  Symbol  Description',
     &                       15x,'Atomic',4x,'Mass',3x,'Valence',/)
               end if
               symbol(k) = symb
               describe(k) = notice
               atmnum(k) = atn
               weight(k) = wght
               ligand(k) = lig
               if (.not. silent) then
                  if (rank.eq.0) then
                    write (iout,20)  k,cls,symb,notice,atn,wght,lig
                  end if
   20             format (2x,i6,1x,i6,5x,a3,3x,a24,i6,f11.3,i6)
               end if
            else if (k .ge. maxtyp) then
               if (rank.eq.0) write (iout,30)
   30          format (/,' KATOM   --  Too many Atom Types;',
     &                    ' Increase MAXTYP')
               abort = .true.
            end if
   40       continue
         end if
      end do
c
c     transfer atom type values to individual atoms
c
      do i = 1, n
         k = type(i)
         if (k .eq. 0) then
            class(i) = 0
            atomic(i) = 0
            mass(i) = 0.0_ti_p
            valence(i) = 0
            story(i) = 'Undefined Atom Type     '
         else
            if (symbol(k) .ne. '   ')  name(i) = symbol(k)
            class(i) = atmcls(k)
            atomic(i) = atmnum(k)
            mass(i) = weight(k)
            valence(i) = ligand(k)
            story(i) = describe(k)
         end if
      end do
c
c     repartition hydrogen masses to use "heavy" hydrogens
c
      heavy = .false.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:15) .eq. 'HEAVY-HYDROGEN ') then
            heavy = .true.
         end if
      end do
      if (heavy) then
         if (hostrank.ne.0) goto 11
         hmax = 4.0d0
         do i = 1, n
            nh = 0
            sum = mass(i)
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  nh = nh + 1
                  sum = sum + mass(k)
               end if
            end do
            hmass = min(hmax,sum/dble(nh+1))
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  dmass = hmass - mass(k)
                  mass(k) = hmass
                  mass(i) = mass(i) - dmass
               end if
            end do
         end do
         do i = 1, n
            if (mass(i) .lt. hmax) then
               dmass = hmax - mass(i)
               dmin = hmax + dmass
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (mass(k) .gt. dmin) then
                     mass(k) = mass(k) - dmass
                     mass(i) = hmax
                     goto 50
                  end if
               end do
               do j = 1, n13(i)
                  k = i13(j,i)
                  if (mass(k) .gt. dmin) then
                     mass(k) = mass(k) - dmass
                     mass(i) = hmax
                     goto 50
                  end if
               end do
   50          continue
            end if
         end do
 11      call MPI_BARRIER(hostcomm,ierr)
      end if
c
c     process keywords containing atom types for specific atoms
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0_ti_p
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:120)
            read (string,*,err=80,end=80)  atn,wght,lig
            if (k.lt.0 .and. k.ge.-n) then
               if (header .and. .not.silent) then
                  header = .false.
                  if (rank.eq.0) write (iout,60)
   60             format (/,' Additional Atom Types for',
     &                       ' Specific Atoms :',
     &                    //,5x,'Atom  Class  Symbol  Description',
     &                       15x,'Atomic',4x,'Mass',3x,'Valence',/)
               end if
               k = -k
               if (cls .eq. 0)  cls = k
               class(k) = cls
               name(k) = symb
               story(k) = notice
               atomic(k) = atn
               mass(k) = wght
               valence(k) = lig
               if ((.not. silent).and.(rank.eq.0)) then
                  write (iout,70)  k,cls,symb,notice,atn,wght,lig
   70             format (2x,i6,1x,i6,5x,a3,3x,a24,i6,f11.3,i6)
               end if
            end if
   80       continue
         end if
      end do
c
c     check for presence of undefined atom types or classes
c
      header = .true.
      do i = 1, n
         k = type(i)
         cls = class(i)
         if (k.lt.1 .or. k.gt.maxtyp
     &          .or. cls.lt.1 .or. cls.gt.maxclass) then
            abort = .true.
            if (header) then
               header = .false.
               if (rank.eq.0) write (iout,90)
   90          format (/,' Undefined Atom Types or Classes :',
     &                 //,' Type',10x,'Atom Number',5x,'Atom Type',
     &                    5x,'Atom Class',/)
            end if
            if (rank.eq.0) write (iout,100)  i,k,cls
  100       format (' Atom',12x,i5,10x,i5,10x,i5)
         end if
      end do
      call mpi_barrier(hostcomm,ierr)
      call update_device_katom
c
c     check the number of atoms attached to each atom
c
      header = .true.
      do i = 1, n
         if (n12(i) .ne. valence(i)) then
            if (header) then
               header = .false.
               if (rank.eq.0) write (iout,110)
  110          format (/,' Atoms with an Unusual Number of Attached',
     &                    ' Atoms :',
     &                 //,' Type',11x,'Atom Name',6x,'Atom Type',7x,
     &                    'Expected',4x,'Found',/)
            end if
            if (rank.eq.0) then
               write (iout,120)  i,name(i),type(i),valence(i),
     $          n12(i)
            end if
  120       format (' Valence',7x,i5,'-',a3,8x,i5,10x,i5,5x,i5)
         end if
      end do
      return
      end

      subroutine update_device_katom
      use atmtyp
      use atmlst
      use inform,only: deb_Path
      implicit none
#ifdef _OPENACC
      if (deb_Path) write(*,*) ' update_device_katom'
!$acc update device(mass,class,atomic)
#endif
      end subroutine
c
c     subroutine alloc_shared_katom : allocate shared memory pointers for katom
c     parameter arrays
c
      subroutine alloc_shared_katom
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atmtyp
      use atoms
      use domdec
      use mpi
      use tinMemory
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
      if (associated(atomic).and.n.eq.size(atomic)) return !Exit condidion

      call shmem_int_req(class,winclass,[n],config=mhostacc)
      call shmem_request(atomic,winatomic,[n],config=mhostacc)
      call shmem_requestm(mass,winmass,[n],config=mhostacc)
      call shmem_request(valence,winvalence,[n])

      if (associated(story)) then
         s_shmem = s_shmem - sizeof(story)
         call MPI_Win_shared_query(winstory, 0, windowsize,
     $        disp_unit, baseptr, ierr)
         call MPI_Win_free(winstory,ierr)
      end if
c
c     story
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*24_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c     allocation
c
      call MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $     hostcomm, baseptr, winstory, ierr)
      if (hostrank /= 0) then
         call MPI_Win_shared_query(winstory, 0, windowsize, disp_unit,
     $        baseptr, ierr)
      end if
c
c     association with fortran pointer
c
      CALL C_F_POINTER(baseptr,story,arrayshape)
      s_shmem = s_shmem + sizeof(story)
      end
