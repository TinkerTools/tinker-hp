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
      use chgpen
      use cutoff
      use domdec
      use inform
      use iounit
      use keys
      use kpolr
      use mpole
      use neigh
      use pme
      use polar
      use polpot
      use potent
      use mpi
      use uprior
      implicit none
      integer istep,modnl,ierr
      integer i,j,k,it
      integer iproc,iglob,polecount,iipole
      integer npg,next
      integer pg(maxvalue)
      real*8 pol,thl
      real*8 sixth
      real*8 d
      logical header
      character*20 keyword
      character*240 record
      character*240 string
      logical init
c
      if (init) then
c
c     deallocate global pointers if necessary
c
        call dealloc_shared_polar
c
c     allocate global pointers
c
        call alloc_shared_polar
c
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
              string = record(next:240)
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
                 string = record(next:240)
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
c       find and store the atomic dipole polarizability parameters
c
        sixth = 1.0d0 / 6.0d0
        npolar = n
        do i = 1, n
           polarity(i) = 0.0d0
           thole(i) = 0.0d0
           dirdamp(i) = 0.0d0
           pdamp(i) = 0.0d0
           it = type(i)
           if (it .ne. 0) then
              polarity(i) = polr(it)
              thole(i) = athl(it)
              dirdamp(i) = ddir(it)
              if (thole(i) .eq. 0.0d0) then
                pdamp(i) = 0.0d0
              else
                pdamp(i) = polarity(i)**sixth
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
        if ((use_polar .or. use_repuls) .and. .not.use_chgtrn) then
          npole = 0
          npolar = 0
          ncp = 0
          ipole = 0
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
                mono0(npole) = pole(1,i)
                if (polarity(i) .ne. 0.0d0)  npolar = npolar + 1
                if (dirdamp(i) .ne. 0.0d0)  use_dirdamp = .true.
                polarity(npole) = polarity(i)
                thole(npole) = thole(i)
                dirdamp(npole) = dirdamp(i)
                pdamp(npole) = pdamp(i)
                if (palpha(i) .ne. 0.0d0)  ncp = ncp + 1
                pcore(npole) = pcore(i)
                pval(npole) = pval(i)
                pval0(npole) = pval(i)
                palpha(npole) = palpha(i)
             end if
          end do
        end if
 90     call MPI_BARRIER(hostcomm,ierr)
        call MPI_BCAST(npole,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(npolar,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(ncp,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(xaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(yaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(zaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(use_dirdamp,1,MPI_LOGICAL,0,hostcomm,ierr)
c
c       test multipoles at chiral sites and invert if necessary
c
        if (use_polar .and. .not.use_chgtrn)  call chkpole(.true.)
c
c       initialization for TCG and omega fit
c
        if ((polalg.eq.3).and.tcgpeek) then 
           poleps = 0.00000001
        end if
c
c       turn off polarizable multipole potential if it is not used
c
        if (npole .eq. 0)  then
          use_mpole = .false.
        end if
        if (npolar .eq. 0)  use_polar = .false.
        if (ncp .ne. 0)  use_chgpen = .true.
        if (ncp .ne. 0)  use_thole = .false.
        if (use_dirdamp)  use_thole = .true.
c
c  allocate predictor arrays
c
        if (use_polar) then
          if (allocated(udalt))  deallocate (udalt)
          if (allocated(upalt))  deallocate (upalt)
          allocate (udalt(maxualt,3,n))
          allocate (upalt(maxualt,3,n))
          if (allocated(udshortalt))  deallocate (udshortalt)
          if (allocated(upshortalt))  deallocate (upshortalt)
          allocate (udshortalt(maxualt,3,n))
          allocate (upshortalt(maxualt,3,n))
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
c       copy original polarizability values that won't change during mutation
c
        polarity_orig = polarity
c
      end if

      if ((use_polar .or. use_repuls) .and. .not.use_chgtrn) then

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
            if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0.0d0)
     $        then
              npolebloc = npolebloc + 1
              poleglob(npolebloc) = polecount + 1
              poleloc(polecount+1) = npolebloc
            end if
          end do
          if (domlen(p_recep1(iproc)+1).ne.0) then
            domlenpole(p_recep1(iproc)+1) = 
     $        npolebloc-bufbegpole(p_recep1(iproc)+1)+1
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

        do i = nlocrec+1, nlocrec2
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
        npolerecloc = domlenpolerec(rank+1)
c
        if ((polalg.eq.3).and.(tcgpeek).and.(tcgomegafit)) then
          if (mod((istep-1), omegafitfreq).eq.0) then
            omegafitstep = .true.
          end if
        end if
c
c  deallocate/reallocate B-spline arrays
c
        if (allocated(thetai1)) deallocate (thetai1)
        if (allocated(thetai2)) deallocate (thetai2)
        if (allocated(thetai3)) deallocate (thetai3)
        allocate (thetai1(4,bsorder,nlocrec))
        allocate (thetai2(4,bsorder,nlocrec))
        allocate (thetai3(4,bsorder,nlocrec))
c
        modnl = mod(istep,ineigup)
        if (istep.eq.-1) return
        if (modnl.ne.0) return
        if (allocated(poleglobnl)) deallocate(poleglobnl)
        allocate (poleglobnl(nlocnl))
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
      end if
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
         call sort (np11(i),ip11(:,i))
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
c     subroutine dealloc_shared_polar : deallocate shared memory pointers for polar
c     parameter arrays
c
      subroutine dealloc_shared_polar
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use polar
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
c
      if (associated(polarity)) then
        CALL MPI_Win_shared_query(winpolarity, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpolarity,ierr)
      end if
      if (associated(polarity_orig)) then
        CALL MPI_Win_shared_query(winpolarity_orig, 0, windowsize,
     $  disp_unit,baseptr, ierr)
        CALL MPI_Win_free(winpolarity_orig,ierr)
      end if
      if (associated(thole)) then
        CALL MPI_Win_shared_query(winthole, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winthole,ierr)
      end if
      if (associated(pdamp)) then
        CALL MPI_Win_shared_query(winpdamp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpdamp,ierr)
      end if
      if (associated(dirdamp)) then
        CALL MPI_Win_shared_query(windirdamp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(windirdamp,ierr)
      end if
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
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1)
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
     $  hostcomm, baseptr, winpolarity, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpolarity, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,polarity,arrayshape)
c
c     polarity_orig
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
     $  hostcomm, baseptr, winpolarity_orig, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpolarity_orig, 0, windowsize,
     $  disp_unit,baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,polarity_orig,arrayshape)
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
     $  hostcomm, baseptr, winthole, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winthole, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winpdamp, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpdamp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pdamp,arrayshape)
c
c     dirdamp
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
     $  hostcomm, baseptr, windirdamp, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(windirdamp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,dirdamp,arrayshape)
      return
      end
!===================================================
!     sub omegafit 
!===================================================
! Objective : find the omega_fit
! Every 'N' timesteps, where N has to be thought,
!  - compute polarization using CG + 10-8 convergence,
!  - use epolar3(mu_CG) as an energy reference 
!  - compute polarization using TCG,
!  - dichotomy-search the right omega minimizing
!      f(omega) = ep(mu_CG) - ep(mu_TCG(omega))
      subroutine omegafit
      use atmlst
      use atoms
      use chgpot
      use domdec
      use energi
      use mpole
      use polpot
      use sizes
      use polar
      use potent
      use mpi
      implicit none

      real*8, allocatable, dimension(:,:,:) :: 
     $       musave
      real*8 :: omega_min, omega_max,
     $           o1, o2, cvg_crit,f, ep1, ep2, sprod,
     $           cnp, cr, efinal, epnp, epsave
      real*8 :: rms1,rms2
      logical :: cvged, dofit
      integer :: n_max_iter, i, iipole, j, ierr

      if (allocated(residue)) deallocate(residue)
      if (allocated(munp))    deallocate(munp)
      if (allocated(efres))   deallocate(efres)
      if (allocated(musave))   deallocate(musave)
      allocate(munp(3,npolebloc))
      allocate(residue(3,npolebloc))
      allocate(efres(3,npolebloc))
      allocate(musave(3,2,npoleloc))


c      omegafitstep = .true.
!      write(*,*)  'tcgomega',tcgomega
      omega_min = -5d0
      omega_max = 5d0
      n_max_iter = 20
      f = electric/dielec
c
      !-2. Recover muTCG NOTPEEKED, residue, field from TCG (via
      !      epolar1tcg) in 'munp', 'residue', 'efres'

      call epolar
c
      do i = 1, npoleloc
        iipole = poleglob(i)
        musave(:,1,i) = uind(:,iipole)
        musave(:,2,i) = uinp(:,iipole)
      end do
      !-1. Recover the energy from tightly converged CG polarization
      !     (via epolar0c)
      epsave = ep
      polalg = 1
      poleps = 0.0000000001
      call epolar
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,ep,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      else
        call MPI_REDUCE(ep,ep,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      end if
      epcg = ep
      rms1 = 0d0
      rms2 = 0d0
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
           rms1 = rms1 + (musave(j,1,i)-uind(j,iipole))**2
           rms2 = rms2 + (musave(j,2,i)-uinp(j,iipole))**2
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,rms1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rms2,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      rms1 = sqrt(rms1/(3*npole))
      rms2 = sqrt(rms2/(3*npole))
      if (rank.eq.0) write(*,*) 'rms init = ',max(rms1,rms2)
c

c
      polalg = 3
      ep = epsave
      cvg_crit = abs(epcg*.0001d0)
      if (cvg_crit .le. 1d-6) cvg_crit = 1d-6

      !0. diagvec r to get alpha.r
      if (.not.(tcgprec)) call diagvec(1, residue,residue)

      !4. Dichotomy !

      ! <E,mupeek> = <munp,E> + omega*<alpha.r,E>
      cnp = sprod(3*npoleloc,munp, efres)
      cr  = sprod(3*npoleloc,residue, efres)
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,cnp,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(MPI_IN_PLACE,cr,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      else
        call MPI_REDUCE(cnp,cnp,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
        call MPI_REDUCE(cr,cr,1,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      end if

      o1 = omega_min
      o2 = omega_max
      ep1 = -.5d0*f*(cnp + o1*cr)
      ep2 = -.5d0*f*(cnp + o2*cr)
      epnp = -.5d0*f*cnp

!      if (rank.eq.0) then
!         dofit=.true.
!         cvged=.false.
!        ! If <residue,E> is too small, don't omegafit
!        if (abs(.5d0*f*cr) .le. cvg_crit) then
!           write(*,*) "Peek step too small, residue energy :",-.5d0*f*cr
!           write(*,*)  
!           dofit = .false.
!        end if
!        
!        if (dofit) then
!           do i = 1, n_max_iter
!              if (abs(epcg-ep1).le. cvg_crit) then 
!                 tcgomega = o1
!                 cvged=.true.
!                 efinal = ep1
!                 exit
!              else if (abs(epcg-ep2).le. cvg_crit) then 
!                 tcgomega = o2
!                 cvged=.true.
!                 efinal = ep2
!                 exit
!              end if 
!              
!              if ((abs(epcg-ep1) .ge. abs(epcg-ep2))) then
!                 o1 = o1 + (o2-o1)*.5d0
!                 ep1 = -.5d0*f*(cnp + o1*cr)
!              else if (abs(epcg-ep1).le.abs(epcg-ep2)) then
!                 o2 = o2 - (o2-o1)*.5d0
!                 ep2 = -.5d0*f*(cnp + o2*cr)
!              end if
!           end do
!        end if
!
!        if (cvged) then
!         write(*,160) 'Omega_fit :', tcgomega, ' (conv. at ', 
!     $                       cvg_crit, ' kcal/mole)'
!         write(*,161)  'Initial (nopeek) energy :', epnp, ' kcal/mole'
!         write(*,161)  'Reference (CG) pol. energy :', epcg,' kcal/mole'
!         write(*,161)  'Peeked TCG pol. energy :', efinal, ' kcal/mole'
!        else 
!           write(*,*)  'Omega fitting could not converge.'
!           write(*,162) 'Omega_fit (not converged):', tcgomega
!           !STOP
!        end if
!      end if


      if (rank.eq.0) then
         dofit = .true.
         ! If <residue,E> is too small, don't omegafit
         if (abs(.5d0*f*cr) .le. cvg_crit) then
            dofit = .false.
            cvged = .false.
         end if

         if (dofit) then
            cvged = .false.
            do i = 1, n_max_iter
               if (abs(epcg-ep1).le. cvg_crit) then 
                  tcgomega = o1
                  cvged=.true.
                  efinal = ep1
                  exit
               else if (abs(epcg-ep2).le. cvg_crit) then 
                  tcgomega = o2
                  cvged=.true.
                  efinal = ep2
                  exit
               end if 
               
               if ((abs(epcg-ep1) .ge. abs(epcg-ep2))) then
                  o1 = o1 + (o2-o1)*.5d0
                  ep1 = -.5d0*f*(cnp + o1*cr)
               else if (abs(epcg-ep1).le.abs(epcg-ep2)) then
                  o2 = o2 - (o2-o1)*.5d0
                  ep2 = -.5d0*f*(cnp + o2*cr)
               end if
            end do
         end if

         if (.not. dofit) then
            write(*,*) "< mu_peek, Efi > to small to fit omega", 
     $                 "(residue energy : ", -.5d0*f*cr, ")"
         end if
         if (cvged) then 
            write(*,160) 'Omega_fit :', tcgomega, ' (conv. at ', 
     $                          cvg_crit, ' kcal/mole)'
            write(*,161) 'Initial (nopeek) energy :', epnp, ' kcal/mole'
            write(*,161) 'Reference (CG) pol. energy :', epcg,
     $                        ' kcal/mole'
            write(*,161) 'Peeked TCG pol. energy :', efinal,' kcal/mole'
         else 
            if (.not. dofit) then
               write(*,*)  'Omega fitting could not converge (cr too',
     $                     'small)'
            else if (dofit) then
               write(*,*)  'Omega fitting could not converge (unable',
     $                     ' to fit)'
            end if
            write(*,162) 'Omega_fit (not converged):', tcgomega
         end if
      end if







c
c     send the converged value of omega to everybody
c
      call MPI_BCAST(tcgomega,1,MPI_REAL8,0,COMM_TINKER,ierr)
c
c     store converged dipole values for rms
c
      do i = 1, npoleloc
        iipole = poleglob(i)
        musave(:,1,i) = uind(:,iipole)
        musave(:,2,i) = uinp(:,iipole)
      end do

      call epolar

      rms1 = 0d0
      rms2 = 0d0
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
           rms1 = rms1 + (musave(j,1,i)-uind(j,iipole))**2
           rms2 = rms2 + (musave(j,2,i)-uinp(j,iipole))**2
        end do
      end do
      rms1 = rms1
      call MPI_ALLREDUCE(MPI_IN_PLACE,rms1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rms2,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      rms1 = sqrt(rms1/(3*npole))
      rms2 = sqrt(rms2/(3*npole))
      if (rank.eq.0) then
         write(*,*) 'rms fin = ',max(rms1,rms2)
         if (max(rms1, rms2) .ge. 0.01) then
            write(*,*) 'WARNING: dipoles rmse greater than 0.01'
         end if
      end if
       

160   format('   ', a11, f6.3, a12, g9.2, a11)
161   format('   ', a28, g10.3, a11)
162   format('   ', a28, g10.3, a11)

      omegafitstep = .false.

      deallocate(munp,residue,efres,musave)

      return
      end
