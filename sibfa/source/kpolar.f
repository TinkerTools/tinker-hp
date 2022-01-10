c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c
      subroutine kpolar(init,istep)
      use atmlst
      use atoms
      use domdec
      use inform
      use iounit
      use keys
      use kpolr
      use mpole
      use neigh
      use polar
      use potent
      use mpi
      use uprior
      implicit none
      integer istep,modnl,ierr
      integer i,j,k
      integer iproc,iglob,polecount,iipole
      integer npg,next
      integer pg(maxvalue)
      real*8 pol,thl
      real*8 sixth
      real*8 d
      logical header
      character*20 keyword
      character*120 record
      character*120 string
      logical init
c
      if (init) then
c
c     allocate global arrays
c
        call alloc_shared_polar
c
        if (.not.(use_polar)) return
        if (hostrank.ne.0) goto 90
c
c       process keywords containing polarizability parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:9) .eq. 'POLARIZE ') then
              k = 0
              pol = 0.0d0
              thl = -1.0d0
              do j = 1, maxvalue
                 pg(j) = 0
              end do
              call getnumb (record,k,next)
              string = record(next:120)
              read (string,*,err=10,end=10) pol,thl,(pg(j),j=1,maxvalue)
   10         continue
              if (k .gt. 0) then
                 if (header .and. .not.silent) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Atomic Dipole',
     &                         ' Polarizability Parameters :',
     &                      //,5x,'Atom Type',11x,'Alpha',8x,
     &                         'Damp',5x,'Group Atom Types'/)
                 end if
                 if (k .le. maxtyp) then
                    polr(k) = pol
                    athl(k) = thl
                    do j = 1, maxvalue
                       pgrp(j,k) = pg(j)
                       if (pg(j) .eq. 0) then
                          npg = j - 1
                          goto 30
                       end if
                    end do
   30               continue
                    if (.not. silent) then
                       if (rank.eq.0) write (iout,40)  k,pol,thl,
     &                    (pg(j),j=1,npg)
   40                  format (4x,i6,10x,f10.3,2x,f10.3,7x,20i5)
                    end if
                 else
                    if (rank.eq.0) write (iout,50)
   50               format (/,' KPOLAR  --  Too many Dipole',
     &                         ' Polarizability Parameters')
                    abort = .true.
                 end if
              end if
           end if
        end do
c
c       find and store the atomic dipole polarizability parameters
c
        do i = 1, n
           polarity(i) = polr(type(i))
           thole(i) = athl(type(i))
        end do
c
c       process keywords containing atom specific polarizabilities
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:9) .eq. 'POLARIZE ') then
              k = 0
              pol = 0.0d0
              thl = 0.0d0
              call getnumb (record,k,next)
              if (k.lt.0 .and. k.ge.-n) then
                 k = -k
                 string = record(next:120)
                 read (string,*,err=60,end=60)  pol,thl
   60            continue
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,70)
   70               format (/,' Additional Dipole Polarizabilities',
     &                         ' for Specific Atoms :',
     &                      //,6x,'Atom',15x,'Alpha',8x,'Damp',/)
                 end if
                 if (.not. silent) then
                    if (rank.eq.0) write (iout,80)  k,pol,thl
   80               format (4x,i6,10x,f10.3,2x,f10.3)
                 end if
                 polarity(k) = pol
                 thole(k) = thl
              end if
           end if
        end do
c
c         assign polarization group connectivity of each atom
c
        call polargrp
c
c       remove zero and undefined polarizable sites from the list
c
        npole = 0
        npolar = 0
        do i = 1, n
           if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
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
              if (polarity(i) .ne. 0.0d0)  npolar = npolar + 1
              polarity(npole) = polarity(i)
              thole(npole) = thole(i)
              if (use_emtp) then
                alphapen(npole) = alphapen(i)
                betapen(npole) = betapen(i)
                gammapen(npole) = gammapen(i)
              end if
           end if
        end do
c
c       set the values used in the scaling of the polarizability
c
        sixth = 1.0d0 / 6.0d0
        do i = 1, npole
           if (thole(i) .eq. 0.0d0) then
              pdamp(i) = 0.0d0
           else
              pdamp(i) = polarity(i)**sixth
           end if
        end do
 90     call MPI_BARRIER(hostcomm,ierr)
        call MPI_BCAST(npole,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(npolar,1,MPI_INT,0,hostcomm,ierr)
c
c       test multipoles at chiral sites and invert if necessary
c
        call chkpole
c
c       turn off polarizable multipole potential if it is not used
c
        if (npole .eq. 0)  use_mpole = .false.
        if (npolar .eq. 0)  use_polar = .false.
c
c  allocate predictor arrays
c
        if (use_polar) then
          if (allocated(udalt))  deallocate (udalt)
          if (allocated(upalt))  deallocate (upalt)
          allocate (udalt(maxualt,3,n))
          allocate (upalt(maxualt,3,n))
c
c         set the Gear predictor binomial coefficients
c
          gear(1) = 6.0d0
          gear(2) = -15.0d0
          gear(3) = 20.0d0
          gear(4) = -15.0d0
          gear(5) = 6.0d0
          gear(6) = -1.0d0
          gear(7) = 0.0d0
c
c         set always stable predictor-corrector (ASPC) coefficients
c
          aspc(1) = 22.0d0 / 7.0d0
          aspc(2) = -55.0d0 / 14.0d0
          aspc(3) = 55.0d0 / 21.0d0
          aspc(4) = -22.0d0 / 21.0d0
          aspc(5) = 5.0d0 / 21.0d0
          aspc(6) = -1.0d0 / 42.0d0
          aspc(7) = 0.0d0
c
c         initialize prior values of induced dipole moments
c
          nualt = 0
          udalt = 0.0d0
        end if
c
      end if
      npoleloc = 0
      do i = 1, nloc
         iglob = glob(i)
         iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
         if (iipole.eq.0) cycle
         polecount = nbpole(iglob)
         if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0.0d0) then
            npoleloc = npoleloc + 1
            poleglob(npoleloc) = polecount + 1
            poleloc(polecount+1) = npoleloc
         end if
      end do
      bufbegpole(rank+1) = 1
      domlenpole(rank+1) = npoleloc
      npolebloc = npoleloc
      do iproc = 1, n_recep1
        if (domlen(p_recep1(iproc)+1).ne.0) then
          bufbegpole(p_recep1(iproc)+1) = npolebloc + 1
        else
          bufbegpole(p_recep1(iproc)+1) = 1
        end if
        do i = 1, domlen(p_recep1(iproc)+1)
          iglob = glob(bufbeg(p_recep1(iproc)+1)+i-1)
          iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
          if (iipole.eq.0) cycle
          polecount = nbpole(iglob)
          if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0.0d0) then
            npolebloc = npolebloc + 1
            poleglob(npolebloc) = polecount + 1
            poleloc(polecount+1) = npolebloc
          end if
        end do
        if (domlen(p_recep1(iproc)+1).ne.0) then
          domlenpole(p_recep1(iproc)+1) = 
     $      npolebloc-bufbegpole(p_recep1(iproc)+1)+1
        else
          domlenpole(p_recep1(iproc)+1) = 0
        end if
      end do
c
c
      npolerecloc = 0
      do i = 1, nlocrec
         iglob = globrec(i)
         iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
         if (iipole.eq.0) cycle
         polecount = nbpole(iglob)
         if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0.0d0) then
            npolerecloc = npolerecloc + 1
            polerecglob(npolerecloc) = polecount + 1
            polerecloc(polecount+1) = npolerecloc
         end if
      end do
      domlenpolerec(rank+1) = npolerecloc
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return
      if (allocated(poleglobnl)) deallocate(poleglobnl)
      allocate (poleglobnl(n))
c
      npolelocnl = 0
      do i = 1, nlocnl
        iglob = ineignl(i)
        iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
        if (iipole.eq.0) cycle
        polecount = nbpole(iglob)
        if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0.0d0) then
          call distprocpart(iglob,rank,d,.true.)
          if (repart(iglob).eq.rank) d = 0.0d0
          if (d*d.le.(mbuf2/4)) then
            npolelocnl = npolelocnl + 1
            poleglobnl(npolelocnl) = polecount + 1
            polelocnl(polecount+1) = npolelocnl
          end if
        end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polargrp  --  polarization group connectivity  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polargrp" generates members of the polarization group of
c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
c     connectivities
c
c
      subroutine polargrp
      use sizes
      use atoms
      use couple
      use domdec
      use inform
      use iounit
      use kpolr
      use polgrp
      implicit none
      integer maxlist,maxkeep
      parameter (maxkeep=100)
      parameter (maxlist=1000)
      integer i,j,k
      integer it,jt
      integer jj,kk
      integer start,stop
      integer nlist,nkeep
      integer keep(maxkeep)
      integer list(maxlist)
      integer, allocatable :: mask(:)
      logical done
c
c
c     find the directly connected group members for each atom
c
      do i = 1, n
         np11(i) = 1
         ip11(1,i) = i
         it = type(i)
         do j = 1, n12(i)
            jj = i12(j,i)
            jt = type(jj)
            do k = 1, maxvalue
               kk = pgrp(k,it)
               if (kk .eq. 0)  goto 20
               if (pgrp(k,it) .eq. jt) then
                  np11(i) = np11(i) + 1
                  if (np11(i) .le. maxp11) then
                     ip11(np11(i),i) = jj
                  else
                     if (rank.eq.0) write (iout,10)
   10                format (/,' POLARGRP  --  Too many Atoms',
     &                          ' in Polarization Group')
                     abort = .true.
                     goto 30
                  end if
               end if
            end do
   20       continue
         end do
      end do
   30 continue
c
c     perform dynamic allocation of some local arrays
c
      allocate (mask(n))
c
c     find any other group members for each atom in turn
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         done = .false.
         start = 1
         stop = np11(i)
         do j = start, stop
            jj = ip11(j,i)
            if (jj .lt. i) then
               done = .true.
               np11(i) = np11(jj)
               do k = 1, np11(i)
                  ip11(k,i) = ip11(k,jj)
               end do
            else
               mask(jj) = i
            end if
         end do
         do while (.not. done)
            done = .true.
            do j = start, stop
               jj = ip11(j,i)
               do k = 1, np11(jj)
                  kk = ip11(k,jj)
                  if (mask(kk) .ne. i) then
                     np11(i) = np11(i) + 1
                     if (np11(i) .le. maxp11) then
                        ip11(np11(i),i) = kk
                     else
                        if (rank.eq.0) write (iout,40)
   40                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                        goto 50
                     end if
                     mask(kk) = i
                  end if
               end do
            end do
            if (np11(i) .ne. stop) then
               done = .false.
               start = stop + 1
               stop = np11(i)
            end if
         end do
         call sort (np11(i),ip11(1,i))
      end do
   50 continue
c
c     loop over atoms finding all the 1-2 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         nkeep = 0
         do j = 1, np11(i)
            jj = ip11(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (mask(kk) .ne. i) then
                  nkeep = nkeep + 1
                  keep(nkeep) = kk
               end if
            end do
         end do
         nlist = 0
         do j = 1, nkeep
            jj = keep(j)
            do k = 1, np11(jj)
               kk = ip11(k,jj)
               nlist = nlist + 1
               list(nlist) = kk
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp12) then
            np12(i) = nlist
            do j = 1, nlist
               ip12(j,i) = list(j)
            end do
         else
            if (rank.eq.0) write (iout,60)
   60       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-2 Polarization Group')
            abort = .true.
            goto 70
         end if
      end do
   70 continue
c
c     loop over atoms finding all the 1-3 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np12(i)
            jj = ip12(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp13) then
            np13(i) = nlist
            do j = 1, nlist
               ip13(j,i) = list(j)
            end do
         else
            if (rank.eq.0) write (iout,80)
   80       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-3 Polarization Group')
            abort = .true.
            goto 90
         end if
      end do
   90 continue
c
c     loop over atoms finding all the 1-4 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         do j = 1, np13(i)
            jj = ip13(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np13(i)
            jj = ip13(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp14) then
            np14(i) = nlist
            do j = 1, nlist
               ip14(j,i) = list(j)
            end do
         else
            if (rank.eq.0) write (iout,100)
  100       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-4 Polarization Group')
            abort = .true.
            goto 110
         end if
      end do
  110 continue
c
c     perform deallocation of some local arrays
c
      deallocate (mask)
      return
      end
c
c     subroutine alloc_shared_polar : allocate shared memory pointers for polar
c     parameter arrays
c
      subroutine alloc_shared_polar
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use domdec
      use polar
      use mpi
      implicit none
      integer :: win,win2
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
      if (associated(polarity)) deallocate (polarity)
      if (associated(thole)) deallocate (thole)
      if (associated(pdamp)) deallocate (pdamp)
c
c     polarity
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,polarity,arrayshape)
c
c     thole
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,thole,arrayshape)
c
c     pdamp
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pdamp,arrayshape)
      return
      end
