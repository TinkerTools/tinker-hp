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
#include "tinker_macro.h"
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
      use mpi
      use neigh
      use potent
      use tinMemory
#ifdef _OPENACC
      use thrust
#endif
      use utilgpu
      implicit none
      integer i,j,k,m,iglob,ionloc
      integer ia,next,ierr,istep,iproc
      integer idomlen,ibufbeg,icap
      integer modnl,nion_s,n0ion
      integer, allocatable :: list(:)
      integer, allocatable :: nc12(:)
      real(t_p) cg,d,dcut2
      logical header
      character*20 keyword
      character*240 record
      character*240 string
      logical init
!$acc routine(distprocpart1)
c
      if (init) then
c
c       deallocate global pointers if necessary
c
        if (deb_Path) print*,'kcharge'
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
              cg = 0.0_ti_p
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
           pchg(i) = merge(chg(type(i)),0.0_ti_p,type(i).ne.0)
           pchg0(i) = 0.0_ti_p
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
              cg = 0.0_ti_p
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
        n0ion= 0
        do i = 1, n
           list(i) = 0
           if (pchg(i).ne.0.0_ti_p.or.fuse_chglj) then
              nbchg(i) = nion
              nion = nion + 1
              iion(nion) = i
              jion(nion) = i
              kion(nion) = i
              pchg(nion) = pchg(i)
             pchg0(nion) = pchg(i)
              list(i) = nion
           else
              iion(n-n0ion) = i
              jion(n-n0ion) = i
              kion(n-n0ion) = i
              !list(i)       = n-n0ion
              list(i)       = 0
              n0ion = n0ion + 1
           end if
        end do
        do i = nion+1,n
           pchg(i) = 0.0_ti_p
          pchg0(i) = 0.0_ti_p
        end do
c
c       copy original charge values that won't change during mutation
c
        if (use_lambdadyn) then
           pchg_orig(:) = pchg(:)
        end if
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
           use_clist  = .false.
        end if

        ! Update of delete shared data
        if (use_charge) then
           call upload_device_kcharge
           !allocate global arrays
           call prmem_request(chglocnl ,n)
           call prmem_request(chgloc   ,n)
           call prmem_request(chgrecloc,n)
        else
           call dealloc_shared_chg
           return
        end if
      end if
c
      call prmem_request(chgglob,nbloc,async=.true.)
      call prmem_request(chgrecglob,nbloc,async=.true.)
c
c     remove zero partial charges from the list of local charges
c
      nionloc = 0

!$acc data present(glob,globrec,chglist,chgglob,chgloc,chgrecglob,
!$acc&     chgrecloc)
!$acc parallel loop copy(nionloc)
      do i = 1, nloc
         iglob  = glob(i)
         ionloc = chglist(iglob)
         if (ionloc.eq.0) cycle
!$acc atomic capture
          nionloc = nionloc + 1
          icap    = nionloc
!$acc end atomic
          chgglob(icap) = ionloc
c         chgloc(ionloc)= icap
      end do
      domlenpole(rank+1) = nionloc
      nionbloc = nionloc
c     nionlocloop is nionloc if nionloc is a multiple of 16, or the first one greater
      nionlocloop = merge( nionloc,
     &                     (int(nionloc/16)+1)*16,
     &                     (mod(nionloc,16).eq.0))

      do iproc = 1, n_recep1
        if (domlen(p_recep1(iproc)+1).ne.0) then
          nion_s = nionbloc + 1
          bufbegpole(p_recep1(iproc)+1) = nion_s
        else
          nion_s = 1
          bufbegpole(p_recep1(iproc)+1) = nion_s
        end if
        idomlen = domlen(p_recep1(iproc)+1)
        ibufbeg = bufbeg(p_recep1(iproc)+1)
!$acc parallel loop copy(nionbloc)
        do i = 1, idomlen
          iglob  = glob(ibufbeg+i-1)
          ionloc = chglist(iglob)
          if (ionloc.eq.0) cycle
!$acc atomic capture
          nionbloc = nionbloc + 1
          icap     = nionbloc
!$acc end atomic
          chgglob(icap)  = ionloc
c         chgloc(ionloc) = icap
        end do
        idomlen = nionbloc - nion_s + 1
        domlenpole(p_recep1(iproc)+1) = idomlen 
      end do
c
      nionrecloc = 0
!$acc parallel loop copy(nionrecloc)
      do i = 1, nlocrec
         iglob  = globrec(i)
         ionloc = chglist(iglob)
         if (ionloc.eq.0) cycle
!$acc atomic capture
         nionrecloc = nionrecloc + 1
         icap       = nionrecloc
!$acc end atomic
         chgrecglob(icap)  = ionloc
c        chgrecloc(ionloc) = icap
      end do
      domlenpolerec(rank+1) = nionrecloc
c
!$acc end data
      modnl = mod(istep,ineigup)
      if (modnl.ne.0.or.istep.lt.0) return

      call prmem_request(chgglobnl,nlocnl,async=.true.)
c
      dcut2 = merge(max(vbuf2,cbuf2),cbuf2,fuse_chglj)
      nionlocnl = 0
!$acc parallel loop copy(nionlocnl)
!$acc&         present(ineignl,chglist,repart,chgglobnl,chglocnl)
      do i = 1, nlocnl
        iglob  = ineignl(i)
        ionloc = chglist(iglob)

        if (ionloc.eq.0) cycle
        call distprocpart1(iglob,rank,d,.true.,x,y,z)
        if (repart(iglob).eq.rank) d = 0.0_ti_p
        if (d*d.le.(dcut2/4)) then
!$acc atomic capture
          nionlocnl = nionlocnl + 1
          icap      = nionlocnl
!$acc end atomic
          chgglobnl(icap)  = ionloc
          chglocnl(ionloc) = icap
        end if

      end do

c     nionlocnlloop is nionlocnl if nionlocnl is a multiple of 16, or the first one greater
      nionlocnlloop = merge( nionlocnl,
     &                     (int(nionlocnl/16)+1)*16,
     &                     (mod(nionlocnl,16).eq.0))

      end

c
c     Upload Charge data on device
c
      subroutine upload_device_kcharge
      use atoms
      use charge
      use domdec,only: rank,hostcomm
      use mpi   ,only: MPI_BARRIER
      use sizes
      use potent,only: use_lambdadyn
      use tinMemory
      implicit none
      integer ierr
#ifdef _OPENACC
 12   format(2x,'upload_device_kcharge')
      if(rank.eq.0.and.tinkerdebug) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif
!$acc update device(iion,pchg,chglist,nbchg)
      if (use_lambdadyn) then
!$acc update device(pchg_orig)
      end if
      end subroutine

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
      use tinheader
      implicit none
      integer i,j,k,m
      integer it,kt,bt
      integer ic,kc
      real(t_p), allocatable :: pbase(:)
      logical emprule
c
c
c     set and store MMFF base atomic partial charge values
c
      do i = 1, n
         it = type(i)
         pchg(i) = 0.0_ti_p
         if (it .eq. 107)  pchg(i) = -0.5_ti_p
         if (it .eq. 113) then
            pchg(i) = 0.0_ti_p
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt .eq. 185)  pchg(i) = -0.5_ti_p
            end do
         end if
         if (it .eq. 114)  pchg(i) = -1.0_ti_p / 3.0_ti_p
         if (it .eq. 115)  pchg(i) = -3.0_ti_p
         if (it .eq. 116)  pchg(i) = -0.5_ti_p
         if (it .eq. 118)  pchg(i) = -0.5_ti_p
         if (it .eq. 119)  pchg(i) = -2.0_ti_p / 3.0_ti_p
         if (it .eq. 121)  pchg(i) = -0.25_ti_p
         if (it .eq. 123)  pchg(i) = 1.0_ti_p
         if (it .eq. 124)  pchg(i) = -1.0_ti_p
         if (it .eq. 125)  pchg(i) = -1.0_ti_p
         if (it .eq. 154)  pchg(i) = 1.0_ti_p
         if (it .eq. 156)  pchg(i) = 1.0_ti_p
         if (it .eq. 159)  pchg(i) = 1.0_ti_p
         if (it .eq. 160)  pchg(i) = 1.0_ti_p
         if (it .eq. 161)  pchg(i) = 0.5_ti_p
         if (it .eq. 162)  pchg(i) = 1.0_ti_p / 3.0_ti_p
         if (it .eq. 165)  pchg(i) = 1.0_ti_p
         if (it .eq. 168) then
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt.eq.168 .or. kt.eq.142)  pchg(i) = 1.0_ti_p
            end do
         end if
         if (it .eq. 169)  pchg(i) = -1.0_ti_p
         if (it .eq. 182)  pchg(i) = -0.5_ti_p
         if (it .eq. 183) then
            pchg(i) = -1.0_ti_p
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt .eq. 87)  pchg(i) = -0.5_ti_p
            end do
         end if
         if (it .eq. 195)  pchg(i) = 1.0_ti_p
         if (it .eq. 196)  pchg(i) = 1.0_ti_p
         if (it .eq. 197)  pchg(i) = 1.0_ti_p
         if (it .eq. 201)  pchg(i) = 2.0_ti_p
         if (it .eq. 202)  pchg(i) = 3.0_ti_p
         if (it .eq. 203)  pchg(i) = -1.0_ti_p
         if (it .eq. 204)  pchg(i) = -1.0_ti_p
         if (it .eq. 205)  pchg(i) = -1.0_ti_p
         if (it .eq. 206)  pchg(i) = 1.0_ti_p
         if (it .eq. 207)  pchg(i) = 1.0_ti_p
         if (it .eq. 208)  pchg(i) = 1.0_ti_p
         if (it .eq. 209)  pchg(i) = 2.0_ti_p
         if (it .eq. 210)  pchg(i) = 2.0_ti_p
         if (it .eq. 211)  pchg(i) = 2.0_ti_p
         if (it .eq. 212)  pchg(i) = 1.0_ti_p
         if (it .eq. 213)  pchg(i) = 2.0_ti_p
         if (it .eq. 214)  pchg(i) = 2.0_ti_p
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (pbase(n))
c
c     modify MMFF base charges using a bond increment scheme
c
      do i = 1, n
         pbase(i) = pchg(i)
      end do
      do i = 1, n
         it = type(i)
         ic = class(i)
         if (pbase(i).lt.0.0_ti_p .or. it.eq.162) then
            pchg(i) = (1.0_ti_p-crd(ic)*fcadj(ic)) * pbase(i)
         end if
         do j = 1, n12(i)
            k = i12(j,i)
            kt = type(k)
            kc = class(k)
            if (pbase(k).lt.0.0_ti_p .or. kt.eq.162) then
               pchg(i) = pchg(i) + fcadj(kc)*pbase(k)
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
               if (bci_1(kc,ic) .eq. 1000.0_ti_p) then
                  emprule = .true.
                  goto 10
               end if
            else if (bt .eq. 0) then
               pchg(i) = pchg(i) + bci(kc,ic)
               if (bci(kc,ic) .eq. 1000.0_ti_p) then
                  emprule = .true.
                  goto 10
               end if
            end if
         end do
   10    continue
         if (emprule) then
            pchg(i) = (1.0_ti_p-crd(ic)*fcadj(ic)) * pbase(i)
            do j = 1, n12(i)
               k = i12(j,i)
               kc = class(k)
               pchg(i) = pchg(i) + fcadj(kc)*pbase(i12(j,i))
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
                  if (bci_1(kc,ic) .eq. 1000.0_ti_p) then
                     pchg(i) = pchg(i) + pbci(ic) - pbci(kc)
                  else
                     pchg(i) = pchg(i) + bci_1(kc,ic)
                  end if
               else if (bt .eq. 0) then
                  if (bci(kc,ic) .eq. 1000.0_ti_p) then
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
      deallocate (pbase)
      return
      end
c
c     subroutine dealloc_shared_chg : deallocate shared memory pointers for chg 
c     parameter arrays
c
      subroutine dealloc_shared_chg
      use atoms
      use charge
      use domdec
      use inform,only: deb_Path
      use potent,only: use_lambdadyn
      use sizes
      use tinMemory
      implicit none
 
 12   format(2x,'dealloc_shared_chg')
      if(deb_Path) print 12
      call shmem_request(iion,   winiion,   [0], config=mhostacc)
      call shmem_request(jion,   winjion,   [0])
      call shmem_request(kion,   winkion,   [0])
      call shmem_request(pchg,   winpchg,   [0], config=mhostacc)
      if (use_lambdadyn) then
         call shmem_request(pchg_orig,winpchg_orig,[0],config=mhostacc)
      end if
      call shmem_request(pchg0, winpchg0,   [0], config=mhostacc)
      call shmem_request(nbchg, winnbchg,   [0], config=mhostacc)
      call shmem_request(chglist,winchglist,[0], config=mhostacc)

      end subroutine
c
c     subroutine alloc_shared_chg : allocate shared memory pointers for chg 
c     parameter arrays
c
      subroutine alloc_shared_chg
      use sizes
      use atoms
      use charge
      use domdec
      use potent,only: use_lambdadyn
      use tinMemory
      implicit none
c
      if (associated(iion).and.n.eq.size(iion)) return

      call shmem_request(iion,   winiion,   [n], config=mhostacc)
      call shmem_request(jion,   winjion,   [n])
      call shmem_request(kion,   winkion,   [n])
      call shmem_request(pchg,   winpchg,   [n], config=mhostacc)
      if (use_lambdadyn) then
         call shmem_request(pchg_orig,winpchg_orig,[n],config=mhostacc)
      end if
      call shmem_request(pchg0, winpchg0,   [n], config=mhostacc)
      call shmem_request(nbchg,  winnbchg,  [n], config=mhostacc)
      call shmem_request(chglist,winchglist,[n], config=mhostacc)

      end
