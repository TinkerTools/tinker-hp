c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kdisp  --  dispersion parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kdisp" assigns C6 coefficients and damping parameters for
c     dispersion interactions and processes any new or changed
c     values for these parameters
c
c
#include "tinker_precision.h"
      subroutine kdisp(init,istep)
      use atoms
      use atmlst
      use atmtyp
      use cutoff
      use disp
      use domdec
      use dsppot
      use inform
      use iounit
      use kdsp
      use keys
      use neigh
      use pme
      use potent
      use sizes
      use tinMemory
      use utilgpu
      implicit none
      integer i,k,istep,iglob,iproc,idisploc
      integer ia,ic,next
      integer dispcount,capt,idml,ibbeg
      integer modnl
      real(t_p) cs,adsp
      real(t_p) csixi
      real(t_p) d
      logical header
      character*20 keyword
      character*240 record
      character*240 string
      logical init
!$acc routine(distprocpart1)
c
      if (init) then
c
c     process keywords containing damped dispersion parameters
c
        if (deb_Path) print*, 'kdisp'
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:11) .eq. 'DISPERSION ') then
              k = 0
              cs = 0.0d0
              adsp = 0.0d0
              call getnumb (record,k,next)
              string = record(next:240)
              read (string,*,err=10,end=10)  cs,adsp
   10         continue
              if (k .gt. 0) then
                 if (header .and. .not.silent) then
                    header = .false.
                    write (iout,20)
   20               format (/,' Additional Damped Dispersion',
     &                         ' Parameters :',
     &                      //,5x,'Atom Class',16x,'C6',12x,'Damp',/)
                 end if
                 if (k .le. maxclass) then
                    dspsix(k) = cs
                    dspdmp(k) = adsp
                    if (.not. silent) then
                       write (iout,30)  k,cs,adsp
   30                  format (6x,i6,7x,f15.4,f15.4)
                    end if
                 else
                    write (iout,40)
   40               format (/,' KDISP  --  Too many Damped',
     &                         ' Dispersion Parameters')
                    abort = .true.
                 end if
              end if
           end if
        end do
c
c       allocate global pointers
c
        call alloc_shared_disp
c
c       assign the dispersion C6 values and alpha parameters 
c
        do i = 1, n
           ic = class(i)
           csix(i) = dspsix(ic)
           adisp(i) = dspdmp(ic)
        end do
c
c       process keywords containing atom specific dispersion parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:11) .eq. 'DISPERSION ') then
              ia = 0
              cs = 0.0d0
              adsp = 0.0d0
              string = record(next:240)
              read (string,*,err=70,end=70)  ia,cs,adsp
              if (ia.lt.0 .and. ia.ge.-n) then
                 ia = -ia
                 if (header .and. .not.silent) then
                    header = .false.
                    write (iout,50)
   50               format (/,' Additional Dispersion Values for',
     &                         ' Specific Atoms :',
     &                      //,8x,'Atom',19x,'C6',12x,'Damp',/)
                 end if
                 if (.not. silent) then
                    write (iout,60)  ia,cs,adsp
   60               format (6x,i6,7x,f15.4,f15.4)
                 end if
                 csix(ia) = cs
                 adisp(ia) = adsp
              end if
   70         continue
           end if
        end do
c 
c       remove zero and undefined dispersion sites from the list
c       
        ndisp = 0
        do i = 1, n
           if (csix(i) .ne. 0.0d0) then 
              nbdisp(i)    = ndisp
              ndisp        = ndisp + 1
              idisp(ndisp) = i
              displist(i)  = ndisp
              !csix(ndisp) = csix(i)
              !adisp(ndisp) = adisp(i)
           end if
        end do
c
c       compute pairwise sum of C6 coefficients needed for PME
c
        csixpr = 0.0d0
        if (use_dewald) then
           do i = 1, n
              csixi = csix(i)
              if (csixi.ne.0.0) then
                 do k = 1, n
                    csixpr = csixpr + csixi*csix(k)
                 end do
              end if
           end do
        end if
c
c       turn off the dispersion potential if not used
c
        if (ndisp .eq. 0)  then
          use_disp = .false.
          use_dlist = .false.
          call dealloc_shared_disp
          return
        else
          call upload_dev_shr_disp
        end if

        call prmem_request(displocnl,ndisp,async=.false.)
        if (use_dewald) then
           call prmem_request(disprecglob,n,async=.false.)
           call prmem_request(disprecloc ,n,async=.false.)
        end if
      end if
c
      call prmem_request(dispglob,nbloc,async=.true.)
!$acc data copy(ndisploc,ndispbloc,ndisprecloc,ndisplocnl) async
c 
c     remove zero and undefined dispersion sites from the list
c       
      modnl = mod(istep,ineigup)
!$acc serial async
      ndisploc  = 0
      ndispbloc = ndisploc
      if (use_dewald) ndisprecloc = 0
      if (modnl.eq.0.and.istep.lt.0) ndisplocnl=0
!$acc end serial
!$acc parallel loop async default(present)
      do i = 1, nloc
        iglob = glob(i)
        dispcount = nbdisp(iglob)
         if (csix(iglob) .ne. 0.0_ti_p) then 
!$acc atomic capture
            ndisploc = ndisploc + 1
            capt     = ndisploc
!$acc end atomic
            dispglob(capt) = dispcount + 1
         end if
      end do

      do iproc = 1, n_recep2
         idml  = domlen(p_recep2(iproc)+1)
         ibbeg = bufbeg(p_recep2(iproc)+1)
!$acc parallel loop async default(present)
        do i = 1, idml
          iglob = glob(ibbeg+i-1)
          dispcount = nbdisp(iglob)
          if (csix(iglob) .ne. 0.0_ti_p) then 
!$acc atomic capture
            ndispbloc = ndispbloc + 1
            capt      = ndispbloc
!$acc end atomic
            dispglob(capt) = dispcount + 1
          end if
        end do
      end do
!$acc serial async
      ndispbloc = ndispbloc + ndisploc
!$acc end serial

      if (use_dewald) then
!$acc parallel loop async default(present)
        do i = 1, nlocrec
           iglob = globrec(i)
           idisploc = displist(iglob)
           if (idisploc.eq.0) cycle
!$acc atomic capture
           ndisprecloc = ndisprecloc + 1
           capt        = ndisprecloc
!$acc end atomic
           disprecglob(capt)    = idisploc
           disprecloc(idisploc) = capt
        end do
      end if
c
      if (modnl.eq.0.and.istep.ne.-1) then

      call prmem_request(dispglobnl,nlocnl,async=.true.)
!$acc parallel loop async default(present)
      do i = 1, nlocnl
        iglob     = ineignl(i)
        dispcount = nbdisp(iglob)
        if (csix(iglob) .ne. 0.0_ti_p) then
          call distprocpart1(iglob,rank,d,.true.,x,y,z)
          if (repart(iglob).eq.rank) d = 0.0_ti_p
            if (d*d.le.(dbuf2/4)) then
!$acc atomic capture
              ndisplocnl = ndisplocnl + 1
              capt       = ndisplocnl
!$acc end atomic
              dispglobnl(capt)       = dispcount + 1
              displocnl(dispcount+1) = capt
            end if
        end if
      end do

      end if
!$acc end data
      end

      subroutine upload_dev_shr_disp
      use disp
      implicit none
!$acc update device(idisp,csix,adisp,displist,nbdisp) async
      end subroutine
c
c     subroutine dealloc_shared_disp : deallocate shared memory pointers for dispersion
c     parameter arrays
c
      subroutine dealloc_shared_disp
      use disp
      use inform    ,only: deb_Path
      use tinMemory
      implicit none

      if (deb_Path) print*, "  dealloc_shared_disp"
      call shmem_request(   idisp,   winidisp,[0], config=mhostacc)
      call shmem_request(    csix,    wincsix,[0], config=mhostacc)
      call shmem_request(   adisp,   winadisp,[0], config=mhostacc)
      call shmem_request(displist,windisplist,[0], config=mhostacc)
      call shmem_request(  nbdisp,  winnbdisp,[0], config=mhostacc)

      end
c
c     subroutine alloc_shared_disp : allocate shared memory pointers for dispersion
c     parameter arrays
c
      subroutine alloc_shared_disp
      use atoms
      use disp
      use tinMemory
      implicit none

      call shmem_request( idisp,  winidisp,[n], config=mhostacc)
      call shmem_request(  csix,   wincsix,[n], config=mhostacc)
      call shmem_request( adisp,  winadisp,[n], config=mhostacc)
      call shmem_request(displist,  windisplist,[n], config=mhostacc)
      call shmem_request(nbdisp, winnbdisp, [n], config=mhostacc)

      end
