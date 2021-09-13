c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kcharge  --  assign partial charge parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kcharge" assigns partial charges to the atoms within
c     the structure and processes any new or changed values
c
c
      subroutine kcharge(init,istep)
      use atmlst
      use atmtyp
      use atoms
      use charge
      use chgpot
      use couple
      use cutoff
      use domdec
      use fields
      use keys
      use kchrge
      use inform
      use iounit
      use neigh
      use potent
      use pme
      use mpi
      implicit none
      integer i,j,k,m,iglob,ionloc
      integer ia,next,ierr,istep,iproc
      integer modnl,ioncount
      integer, allocatable :: list(:)
      integer, allocatable :: nc12(:)
      real*8 cg,d
      logical header
      character*20 keyword
      character*240 record
      character*240 string
      logical init
c
      if (init) then
c
c       deallocate global pointers if necessary
c
        call dealloc_shared_chg
c
c       allocate global pointers
c
        call alloc_shared_chg
        if (hostrank.ne.0) goto 1000
c
c       process keywords containing partial charge parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'CHARGE ') then
              ia = 0
              cg = 0.0d0
              string = record(next:240)
              read (string,*,err=40,end=40)  ia,cg
              if (ia .gt. 0) then
                 if (header .and. .not.silent) then
                    header = .false.
                    write (iout,10)
   10               format (/,' Additional Atomic Partial Charge',
     &                         ' Parameters :',
     &                      //,5x,'Atom Type',10x,'Charge',/)
                 end if
                 if (ia .le. maxtyp) then
                    chg(ia) = cg
                    if (.not. silent) then
                       write (iout,20)  ia,cg
   20                  format (4x,i6,8x,f12.4)
                    end if
                 else
                    write (iout,30)
   30               format (/,' KCHARGE  --  Too many Partial Charge',
     &                         ' Parameters')
                    abort = .true.
                 end if
              end if
   40         continue
           end if
        end do
c
c       find and store all the atomic partial charges
c
        do i = 1, n
           pchg(i) = chg(type(i))
        end do
c
c       use special charge parameter assignment method for MMFF
c
        if (forcefield .eq. 'MMFF94')  call kchargem
c
c       process keywords containing atom specific partial charges
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'CHARGE ') then
              ia = 0
              cg = 0.0d0
              string = record(next:240)
              read (string,*,err=70,end=70)  ia,cg
              if (ia.lt.0 .and. ia.ge.-n) then
                 ia = -ia
                 if (header .and. .not.silent) then
                    header = .false.
                    write (iout,50)
   50               format (/,' Additional Partial Charges for',
     &                         ' Specific Atoms :',
     &                      //,6x,'Atom',14x,'Charge',/)
                 end if
                 if (.not. silent) then
                    write (iout,60)  ia,cg
   60               format (4x,i6,8x,f12.4)
                 end if
                 pchg(ia) = cg
              end if
   70         continue
           end if
        end do
c
c       perform dynamic allocation of some local arrays
c
        allocate (list(n))
        allocate (nc12(n))
c
c       remove zero partial charges from the list of charges
c
        nion = 0
        do i = 1, n
           list(i) = 0
           if (pchg(i) .ne. 0.0d0) then
              nbchg(i) = nion
              nion = nion + 1
              iion(nion) = i
              jion(nion) = i
              kion(nion) = i
              pchg(nion) = pchg(i)
              list(i) = nion
           end if
        end do
c
c       perform deallocation of some local arrays
c
        chglist = list
        deallocate (list)
        deallocate (nc12)
 1000   call MPI_BARRIER(hostcomm,ierr)
        call MPI_BCAST(nion,1,MPI_INT,0,hostcomm,ierr)
c
c       turn off charge-charge and charge-dipole terms if not used
c
        if (nion .eq. 0) then
           use_charge = .false.
           use_clist = .false.
        end if
        if (.not.(use_charge)) return
        if (allocated(chglocnl)) deallocate(chglocnl)
        allocate (chglocnl(n))
        if (allocated(chgloc)) deallocate(chgloc)
        allocate (chgloc(n))
        if (allocated(chgrecloc)) deallocate(chgrecloc)
        allocate (chgrecloc(n))
      end if
c
      if (allocated(chgglob)) deallocate(chgglob)
      allocate (chgglob(nbloc))
      if (allocated(chgrecglob)) deallocate(chgrecglob)
      allocate (chgrecglob(nlocrec2))
c
c       remove zero partial charges from the list of local charges
c
      nionloc = 0
      do i = 1, nloc
         iglob = glob(i)
         ionloc = chglist(iglob)
         if (ionloc.eq.0) cycle
          ioncount = nbchg(iglob)
          nionloc = nionloc + 1
          chgglob(nionloc) = ioncount + 1
          chgloc(ioncount+1) = nionloc
      end do
      domlenpole(rank+1) = nionloc
      nionbloc = nionloc
      do iproc = 1, n_recep1
        if (domlen(p_recep1(iproc)+1).ne.0) then
          bufbegpole(p_recep1(iproc)+1) = nionbloc + 1
        else
          bufbegpole(p_recep1(iproc)+1) = 1
        end if
        do i = 1, domlen(p_recep1(iproc)+1)
          iglob = glob(bufbeg(p_recep1(iproc)+1)+i-1)
          ionloc = chglist(iglob)
          if (ionloc.eq.0) cycle
          ioncount = nbchg(iglob)
          nionbloc = nionbloc + 1
          chgglob(nionbloc) = ioncount + 1
          chgloc(ioncount+1) = nionbloc
        end do
        if (domlen(p_recep1(iproc)+1).ne.0) then
          domlenpole(p_recep1(iproc)+1) = 
     $      nionbloc-bufbegpole(p_recep1(iproc)+1)+1
        else
          domlenpole(p_recep1(iproc)+1) = 0
        end if
      end do
c
      nionrecloc = 0
      do i = 1, nlocrec
         iglob = globrec(i)
         ionloc = chglist(iglob)
         if (ionloc.eq.0) cycle
         ioncount = nbchg(iglob)
         nionrecloc = nionrecloc + 1
         chgrecglob(nionrecloc) = ioncount + 1
         chgrecloc(ioncount+1) = nionrecloc
      end do
      domlenpolerec(rank+1) = nionrecloc
      do i = nlocrec+1, nlocrec2
         iglob = globrec(i)
         ionloc = chglist(iglob)
         if (ionloc.eq.0) cycle
         ioncount = nbchg(iglob)
         nionrecloc = nionrecloc + 1
         chgrecglob(nionrecloc) = ioncount + 1
         chgrecloc(ioncount+1) = nionrecloc
      end do
      nionrecloc = domlenpolerec(rank+1)
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
      if (allocated(chgglobnl)) deallocate(chgglobnl)
      allocate (chgglobnl(nlocnl))
c
      nionlocnl = 0
      do i = 1, nlocnl
        iglob = ineignl(i)
        ionloc = chglist(iglob)
        if (ionloc.eq.0) cycle
        ioncount = nbchg(iglob)
        call distprocpart(iglob,rank,d,.true.)
        if (repart(iglob).eq.rank) d = 0.0d0
        if (d*d.le.(cbuf2/4)) then
          nionlocnl = nionlocnl + 1
          chgglobnl(nionlocnl) = ioncount + 1
          chglocnl(ioncount + 1) = nionlocnl
        end if
      end do
c
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kchargem  --  assign MMFF charge parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kchargem" assigns partial charges to the atoms according to
c     the Merck Molecular Force Field (MMFF)
c
c
      subroutine kchargem
      use sizes
      use atmtyp
      use atoms
      use charge
      use couple
      use merck
      implicit none
      integer i,j,k,m
      integer it,kt,bt
      integer ic,kc
      real*8, allocatable :: pchg0(:)
      logical emprule
c
c
c     set and store MMFF base atomic partial charge values
c
      do i = 1, n
         it = type(i)
         pchg(i) = 0.0d0
         if (it .eq. 107)  pchg(i) = -0.5d0
         if (it .eq. 113) then
            pchg(i) = 0.0d0
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt .eq. 185)  pchg(i) = -0.5d0
            end do
         end if
         if (it .eq. 114)  pchg(i) = -1.0d0 / 3.0d0
         if (it .eq. 115)  pchg(i) = -3.0d0
         if (it .eq. 116)  pchg(i) = -0.5d0
         if (it .eq. 118)  pchg(i) = -0.5d0
         if (it .eq. 119)  pchg(i) = -2.0d0 / 3.0d0
         if (it .eq. 121)  pchg(i) = -0.25d0
         if (it .eq. 123)  pchg(i) = 1.0d0
         if (it .eq. 124)  pchg(i) = -1.0d0
         if (it .eq. 125)  pchg(i) = -1.0d0
         if (it .eq. 154)  pchg(i) = 1.0d0
         if (it .eq. 156)  pchg(i) = 1.0d0
         if (it .eq. 159)  pchg(i) = 1.0d0
         if (it .eq. 160)  pchg(i) = 1.0d0
         if (it .eq. 161)  pchg(i) = 0.5d0
         if (it .eq. 162)  pchg(i) = 1.0d0 / 3.0d0
         if (it .eq. 165)  pchg(i) = 1.0d0
         if (it .eq. 168) then
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt.eq.168 .or. kt.eq.142)  pchg(i) = 1.0d0
            end do
         end if
         if (it .eq. 169)  pchg(i) = -1.0d0
         if (it .eq. 182)  pchg(i) = -0.5d0
         if (it .eq. 183) then
            pchg(i) = -1.0d0
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt .eq. 87)  pchg(i) = -0.5d0
            end do
         end if
         if (it .eq. 195)  pchg(i) = 1.0d0
         if (it .eq. 196)  pchg(i) = 1.0d0
         if (it .eq. 197)  pchg(i) = 1.0d0
         if (it .eq. 201)  pchg(i) = 2.0d0
         if (it .eq. 202)  pchg(i) = 3.0d0
         if (it .eq. 203)  pchg(i) = -1.0d0
         if (it .eq. 204)  pchg(i) = -1.0d0
         if (it .eq. 205)  pchg(i) = -1.0d0
         if (it .eq. 206)  pchg(i) = 1.0d0
         if (it .eq. 207)  pchg(i) = 1.0d0
         if (it .eq. 208)  pchg(i) = 1.0d0
         if (it .eq. 209)  pchg(i) = 2.0d0
         if (it .eq. 210)  pchg(i) = 2.0d0
         if (it .eq. 211)  pchg(i) = 2.0d0
         if (it .eq. 212)  pchg(i) = 1.0d0
         if (it .eq. 213)  pchg(i) = 2.0d0
         if (it .eq. 214)  pchg(i) = 2.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (pchg0(n))
c
c     modify MMFF base charges using a bond increment scheme
c
      do i = 1, n
         pchg0(i) = pchg(i)
      end do
      do i = 1, n
         it = type(i)
         ic = class(i)
         if (pchg0(i).lt.0.0d0 .or. it.eq.162) then
            pchg(i) = (1.0d0-crd(ic)*fcadj(ic)) * pchg0(i)
         end if
         do j = 1, n12(i)
            k = i12(j,i)
            kt = type(k)
            kc = class(k)
            if (pchg0(k).lt.0.0d0 .or. kt.eq.162) then
               pchg(i) = pchg(i) + fcadj(kc)*pchg0(k)
            end if
            bt = 0
            do m = 1, nlignes
               if ((i.eq.bt_1(m,1) .and. i12(j,i).eq.bt_1(m,2)).or.
     &                (i12(j,i).eq.bt_1(m,1) .and. i.eq.bt_1(m,2))) then
                  bt = 1
               end if
            end do
            emprule = .false.
            if (bt .eq. 1) then
               pchg(i) = pchg(i) + bci_1(kc,ic)
               if (bci_1(kc,ic) .eq. 1000.0d0) then
                  emprule = .true.
                  goto 10
               end if
            else if (bt .eq. 0) then
               pchg(i) = pchg(i) + bci(kc,ic)
               if (bci(kc,ic) .eq. 1000.0d0) then
                  emprule = .true.
                  goto 10
               end if
            end if
         end do
   10    continue
         if (emprule) then
            pchg(i) = (1.0d0-crd(ic)*fcadj(ic)) * pchg0(i)
            do j = 1, n12(i)
               k = i12(j,i)
               kc = class(k)
               pchg(i) = pchg(i) + fcadj(kc)*pchg0(i12(j,i))
            end do
            do j = 1, n12(i)
               k = i12(j,i)
               kc = class(k)
               bt = 0
               do k = 1, nlignes
                  if ((i.eq.bt_1(k,1) .and.
     &                      i12(j,i).eq.bt_1(k,2)) .or.
     &                   (i12(j,i).eq.bt_1(k,1) .and.
     &                      i.eq.bt_1(k,2))) then
                     bt = 1
                  end if
               end do
               if (bt .eq. 1) then
                  if (bci_1(kc,ic) .eq. 1000.0d0) then
                     pchg(i) = pchg(i) + pbci(ic) - pbci(kc)
                  else
                     pchg(i) = pchg(i) + bci_1(kc,ic)
                  end if
               else if (bt .eq. 0) then
                  if (bci(kc,ic) .eq. 1000.0d0) then
                     pchg(i) = pchg(i) + pbci(ic) - pbci(kc)
                  else
                     pchg(i) = pchg(i) + bci(kc,ic)
                  end if
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pchg0)
      return
      end
c
c     subroutine dealloc_shared_chg : deallocate shared memory pointers for chg 
c     parameter arrays
c
      subroutine dealloc_shared_chg
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use charge
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
c
      if (associated(iion)) then
        CALL MPI_Win_shared_query(winiion, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winiion,ierr)
      end if
      if (associated(jion)) then
        CALL MPI_Win_shared_query(winjion, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winjion,ierr)
      end if
      if (associated(kion)) then
        CALL MPI_Win_shared_query(winkion, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winkion,ierr)
      end if
      if (associated(pchg)) then
        CALL MPI_Win_shared_query(winpchg, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpchg,ierr)
      end if
      if (associated(nbchg)) then
        CALL MPI_Win_shared_query(winnbchg, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbchg,ierr)
      end if
      if (associated(chglist)) then
        CALL MPI_Win_shared_query(winchglist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winchglist,ierr)
      end if
      return
      end
c
c     subroutine alloc_shared_chg : allocate shared memory pointers for chg 
c     parameter arrays
c
      subroutine alloc_shared_chg
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use charge
      use domdec
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
c     iion 
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winiion, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winiion, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iion,arrayshape)
c
c     jion
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winjion, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winjion, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,jion,arrayshape)
c
c     kion
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winkion, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winkion, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,kion,arrayshape)
c
c     pchg
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
     $  hostcomm, baseptr, winpchg, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpchg, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pchg,arrayshape)
c
c     nbchg
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbchg, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbchg, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbchg,arrayshape)
c
c     chglist
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winchglist, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winchglist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,chglist,arrayshape)
      return
      end
