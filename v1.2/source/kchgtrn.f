c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kchgtrn  --  charge transfer term assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kchgtrn" assigns charge magnitude and damping parameters for
c     charge transfer interactions and processes any new or changed
c     values for these parameters
c
c
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
c
      if (init) then
c
c     deallocate global pointers if necessary
c
        call dealloc_shared_chgct
c
c     allocate global pointers
c
        call alloc_shared_chgct

        if (hostrank.ne.0) goto 100
c
c     process keywords containing charge transfer parameters
c
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
   10         continue
              if (k .gt. 0) then
                 if (header .and. .not.silent) then
                    header = .false.
                    write (iout,20)
   20               format (/,' Additional Charge Transfer',
     &                         ' Parameters :',
     &                     //,5x,'Atom Class',13x,'Charge',11x,'Damp',/)
                 end if
                 if (k .le. maxclass) then
                    ctchg(k) = chtrn
                    ctdmp(k) = actrn
                    if (.not. silent) then
                       write (iout,30)  k,chtrn,actrn
   30                  format (6x,i6,7x,f15.4,f15.4)
                    end if
                 else
                    write (iout,40)
   40               format (/,' KCHGTRN  --  Too many Charge',
     &                         ' Transfer Parameters')
                    abort = .true.
                 end if
              end if
           end if
        end do
c
c
c     assign the charge transfer charge and alpha parameters 
c     
        do i = 1, n
           ic = class(i)
           chgct(i) = ctchg(ic)
           dmpct(i) = ctdmp(ic)
        end do
c
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
   50               format (/,' Additional Charge Transfer Values',
     &                         ' for Specific Atoms :',
     &                      //,8x,'Atom',16x,'Charge',11x,'Damp',/)
                 end if
                 if (.not. silent) then
                    write (iout,60)  ia,chtrn,actrn
   60               format (6x,i6,7x,f15.4,f15.4)
                 end if
                 chgct(ia) = chtrn
                 dmpct(ia) = actrn
              end if
   70         continue
           end if
        end do
c
c     remove zero or undefined electrostatic sites from the list
c
        if (use_chgtrn) then
          npole = 0
          ncp = 0
          npolar = 0
          nct = 0
          do i = 1, n
             if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0 .or.
     &              chgct(i).ne. 0.0d0 .or. dmpct(i).ne.0.0d0) then
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
c                   douind(i) = .true.
                end if
                if (dirdamp(i) .ne. 0.0d0)  use_dirdamp = .true.
                polarity(npole) = polarity(i)
                thole(npole) = thole(i)
                dirdamp(npole) = dirdamp(i)
                pdamp(npole) = pdamp(i)
                if (chgct(i).ne.0.0d0 .or. dmpct(i).ne.0.0d0) then
                   nct = nct + 1
                end if
                chgct(npole) = chgct(i)
                dmpct(npole) = dmpct(i)
             end if
          end do
        end if
 100    call MPI_BARRIER(hostcomm,ierr)
        call MPI_BCAST(npole,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(npolar,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(nct,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(ncp,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(xaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(yaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(zaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(use_dirdamp,1,MPI_LOGICAL,0,hostcomm,ierr)
c
c     test multipoles at chiral sites and invert if necessary
c
        if (use_chgtrn) call chkpole(.true.)
c
c     turn off individual electrostatic potentials if not used
c
        if (npole .eq. 0)  use_mpole = .false.
        if (ncp .ne. 0)  use_chgpen = .true.
        if (ncp .ne. 0)  use_thole = .false.
        if (use_dirdamp)  use_thole = .true.
        if (npolar .eq. 0)  use_polar = .false.
        if (nct .eq. 0)  use_chgtrn = .false.
      end if

      if (.not. use_chgtrn) return
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
c
c     also deallocate/reallocate B-spline arrays
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
      return
      end
c
c     subroutine dealloc_shared_chgct : deallocate shared memory pointers for charge transfer
c     parameter arrays
c
      subroutine dealloc_shared_chgct
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use chgtrn
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
c
      if (associated(chgct)) then
        CALL MPI_Win_shared_query(winchgct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winchgct,ierr)
      end if
      if (associated(dmpct)) then
        CALL MPI_Win_shared_query(windmpct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(windmpct,ierr)
      end if
      return
      end
c
c     subroutine alloc_shared_chgct : allocate shared memory pointers for charge transfer
c     parameter arrays
c
      subroutine alloc_shared_chgct
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use domdec
      use chgtrn
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1)
c
c     chgct
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
     $  hostcomm, baseptr, winchgct, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winchgct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,chgct,arrayshape)
c
c     dmpct
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
     $  hostcomm, baseptr, windmpct, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(windmpct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,dmpct,arrayshape)
c
      return
      end
