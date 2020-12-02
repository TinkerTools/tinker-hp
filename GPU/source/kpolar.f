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
#include "tinker_precision.h"
      subroutine kpolar(init,istep)
      use atmlst
      use atoms
      use cutoff
      use domdec
      use inform
      use iounit
      use,intrinsic::iso_c_binding
      use keys
      use kpolr
      use mpole
      use mdstuf ,only: integrate
      use neigh
      use nvshmem
      use polar
      use polpot
      use potent
      use mpi
      use uprior
      use tinheader ,only:ti_p,re_p
      use tinMemory
      implicit none
      integer istep,modnl,ierr
      integer i,j,k,ii,maxdr
      integer iproc,iglob,polecount,iipole
      integer ipr,npole_cap,ibufbeg
      integer,save:: nlocnl_save=0
      integer npg,next
      integer pg(maxvalue)
      integer winalt
      real(t_p),parameter:: zero=0.0_ti_p
      real(t_p) pol,thl
      real(t_p) sixth
      real(t_p) d
      real(t_p),pointer:: buffer3d(:,:,:)
      real(8) m1,m2,m3
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
        if (deb_Path) print*,'kpolar'
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
              pol = zero
              thl = -1.0_ti_p
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
           else if (keyword(1:11) .eq. 'INTEGRATOR ') then
              call getword (record,integrate,next)
              call upcase (integrate)
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
              pol = zero
              thl = zero
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
  90    call polargrp
        if (hostrank.ne.0) goto 100
c
c       remove zero and undefined polarizable sites from the list
c
        npole  = 0
        npolar = 0
        do i = 1, n
           if (polsiz(i).ne.0 .or. polarity(i).ne.zero) then
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
              if (polarity(i) .ne. zero)  npolar = npolar + 1
              polarity(npole) = polarity(i)
              thole(npole) = thole(i)
              if (use_emtp) then
                alphapen(npole) = alphapen(i)
                betapen(npole) = betapen(i)
                gammapen(npole) = gammapen(i)
              end if
           end if
        end do
        call FromPolaxe2Ipolaxe ! Brodcast 'nZ_onlyglob' after this call
c
c       Display if some poles does not have polarity
c
        npolar_ne_npole = (npolar.ne.npole)
        if (npolar.ne.npole) write(*,99) npole, npolar
 99        format("There are some poles without polarity :",/,
     &            3x,I12, "poles for ", I12, "polar atoms.")
c
c       set the values used in the scaling of the polarizability
c
        sixth = 1.0_ti_p / 6.0_ti_p
        do i = 1, npole
           if (thole(i) .eq. zero) then
              pdamp(i) = zero
           else
              pdamp(i) = polarity(i)**sixth
           end if
        end do
 100    call MPI_BARRIER(hostcomm,ierr)
        call MPI_BCAST(npole,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(npolar,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(integrate,11,MPI_CHAR,0,hostcomm,ierr)
        call MPI_BCAST(nZ_Onlyglob,1,MPI_INT,0,hostcomm,ierr)
c
c       test multipoles at chiral sites and invert if necessary
c
        call chkpole(.true.)
c
c       turn off polarizable multipole potential if it is not used
c
        if (npole .eq. 0)  then
           use_mpole = .false.
           use_mlist = .false.
!$acc update device(use_mpole)
           endif
        if (npolar .eq. 0) then
        use_polar = .false.
!$acc update device(use_polar)
        end if
c
c       initialization for TCG and omega fit
c
        if ((polalg.eq.3).and.tcgpeek) then 
           !TODO 1.2  Ask louis about this value for single prec
           poleps = 0.00000001
        end if

        if (use_polar) then
           call upload_device_polar(1)
        else
           return
        end if
c
c  allocate predictor arrays
c
        if (use_polar) then
           if (deb_Path) then
              call mem_get(m1,m2)
           end if
#ifdef USE_NVSHMEM_CUDA
           call shmem_request(buffer3d,winalt,[maxualt,3,n],
     &                c_udalt,d_udalt,config=mnvshonly)
           call shmem_request(buffer3d,winalt,[maxualt,3,n],
     &                c_upalt,d_upalt,config=mnvshonly)
           call shmem_request(buffer3d,winalt,[3,2,n],
     &                c_altbuf,d_altbuf,config=mnvshonly)
           d_udalt  => d_udalt
           d_upalt  => d_upalt
           d_altbuf => d_altbuf
#else
           call prmem_request(udalt,maxualt,3,n)
           call prmem_request(upalt,maxualt,3,n)
#endif
           if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1')
     &        .and.(.not.use_pmecore.or.use_pmecore.and.rank.lt.ndir))
     &        then
#ifdef USE_NVSHMEM_CUDA
              call shmem_request(buffer3d,winalt,[maxualt,3,n],
     &              c_udshortalt,d_udshortalt,config=mnvshonly)
              call shmem_request(buffer3d,winalt,[maxualt,3,n],
     &              c_upshortalt,d_upshortalt,config=mnvshonly)
              d_udshortalt => d_udshortalt
              d_upshortalt => d_upshortalt
#else
              call prmem_request(udshortalt,maxualt,3,n)
              call prmem_request(upshortalt,maxualt,3,n)
#endif
           end if
           if (deb_Path) then
              call mem_get(m1,m3)
 14   format(" Polar Predictor workSpace size",F9.3," Mio per process")
              print 14, m3-m2
           end if
c
c          set the Gear predictor binomial coefficients
c
           gear(1) = 6.0_ti_p
           gear(2) = -15.0_ti_p
           gear(3) = 20.0_ti_p
           gear(4) = -15.0_ti_p
           gear(5) = 6.0_ti_p
           gear(6) = -1.0_ti_p
           gear(7) = 0.0_ti_p
c
c          set always stable predictor-corrector (ASPC) coefficients
c
           aspc(1) = 22.0_ti_p / 7.0_ti_p
           aspc(2) = -55.0_ti_p / 14.0_ti_p
           aspc(3) = 55.0_ti_p / 21.0_ti_p
           aspc(4) = -22.0_ti_p / 21.0_ti_p
           aspc(5) = 5.0_ti_p / 21.0_ti_p
           aspc(6) = -1.0_ti_p / 42.0_ti_p
           aspc(7) = 0.0_ti_p
!$acc update device(gear,aspc)
c
c          initialize prior values of induced dipole moments
c
#ifdef USE_NVSHMEM_CUDA
!$acc parallel loop collapse(3) deviceptr(d_upalt,d_udalt)
           do i = 1, n_pe
              do j = 1, 3
                 do k = 1, maxualt
                    d_udalt(mype)%pel(k,j,i) = 0.0_ti_p
                    d_upalt(mype)%pel(k,j,i) = 0.0_ti_p
                 end do
              end do
           end do
#else
!$acc parallel loop collapse(3) default(present)
           do i = 1, npole
              do j = 1, 3
                 do k = 1, maxualt
                    udalt(k,j,i) = 0.0_ti_p
                    upalt(k,j,i) = 0.0_ti_p
                 end do
              end do
           end do
#endif

           if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1')
     &        .and.(.not.use_pmecore.or.use_pmecore.and.rank.lt.ndir))
     &        then
#ifdef USE_NVSHMEM_CUDA
!$acc parallel loop collapse(3) deviceptr(d_upshortalt)
              do i = 1, n_pe
                 do j = 1, 3
                    do k = 1, maxualt
                       d_udshortalt(mype)%pel(k,j,i) = 0.0_ti_p
                       d_upshortalt(mype)%pel(k,j,i) = 0.0_ti_p
                    end do
                 end do
              end do
#else
!$acc parallel loop collapse(3) default(present)
              do i = 1, npole
                 do j = 1, 3
                    do k = 1, maxualt
                       upshortalt(k,j,i) = 0.0_ti_p
                       udshortalt(k,j,i) = 0.0_ti_p
                    end do
                 end do
              end do
#endif
           end if
        end if
c
!$acc enter data copyin(nZ_onlyloc)
      end if

      if (nproc.eq.1) then
         call kpolar_reset1(istep)
      else
         call kpolar_reset(istep)
      end if

      end subroutine

      subroutine kpolar_reset(istep)
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
      use utilcomm,only:do_not_commpole
#ifdef _OPENACC
      use thrust
#endif
      use utilgpu
      implicit none
      integer,intent(in)::istep
      integer i,j,iipole,iglob,iproc
      integer polecount,npole_cap,modnl
      integer ipr,ibufbeg
      real(t_p) d
!$acc routine(distprocpart1)

!$acc data present(glob,pollist,poleglob,poleloc,nbpole,bufbegpole,
!$acc&  domlenpole,domlenpolerec,polerecglob,poleloc,
!$acc&  polerecloc)
!$acc&     present(npoleloc,npolebloc,npolerecloc,npolerecloc_old)
!$acc&     async(rec_queue)

      !TODO clean this routine
!$acc serial async(rec_queue)
      npoleloc = 0
!$acc end serial
!$acc parallel loop async(rec_queue)
      do i = 1, nloc
         iglob = glob(i)
         iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
         if (iipole.eq.0) cycle
         polecount = nbpole(iglob)
         if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
!$acc atomic capture
            npoleloc = npoleloc + 1
            npole_cap = npoleloc
!$acc end atomic
            poleglob(npole_cap) = polecount + 1
            poleloc(polecount+1) = npole_cap
         end if
      end do
!$acc serial async(rec_queue)
      bufbegpole(rank+1) = 1
      domlenpole(rank+1) = npoleloc
      npolebloc = npoleloc
!$acc end serial
#ifndef _OPENACC
      do iproc = 1, n_recep1
        ipr = p_recep1(iproc)+1
        ibufbeg = bufbeg(ipr)
        if (domlen(ipr).ne.0) then
!$acc serial async
          bufbegpole(ipr) = npolebloc + 1
!$acc end serial
        else
!$acc serial async
          bufbegpole(ipr) = 1
!$acc end serial
        end if
!$acc serial loop async
        do i = 1, domlen(ipr)
          iglob  = glob(ibufbeg+i-1)
          iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
          if (iipole.eq.0) cycle
          polecount = nbpole(iglob)
          if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
!$acc atomic capture
            npolebloc = npolebloc + 1
            npole_cap = npolebloc
!$acc end atomic
            poleglob(npole_cap) = polecount + 1
            poleloc(polecount+1) = npole_cap
          end if
        end do
        if (domlen(ipr).ne.0) then
!$acc serial async
          domlenpole(ipr) = npolebloc-bufbegpole(ipr)+1
!$acc end serial
        else
!$acc serial async
          domlenpole(ipr) = 0
!$acc end serial
        end if
      end do
#else
      do_not_commpole = .false.
#endif
c!$acc update host(poleglob(:),poleloc(:),npoleloc,npolebloc,
c!$acc&  domlenpole,bufbegpole) async
c
!$acc serial async(rec_queue)
      npolerecloc_old = npolerecloc
      npolerecloc  = 0
!$acc end serial
!$acc parallel loop async(rec_queue)
      do i = 1, nlocrec
         iglob  = globrec(i)
         iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
         if (iipole.eq.0) cycle
         polecount = nbpole(iglob)
         if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
!$acc atomic capture
            npolerecloc = npolerecloc + 1
            npole_cap   = npolerecloc
!$acc end atomic
            polerecglob(npole_cap) = polecount + 1
            polerecloc(polecount+1) = npole_cap
         end if
      end do
c!$acc update device(polerecglob(:),polerecloc(:)) async
!$acc serial async(rec_queue)
      domlenpolerec(rank+1) = npolerecloc
!$acc end serial
c
!$acc update host(domlenpole,domlenpolerec,bufbegpole,
!$acc& npoleloc,npolebloc,npolerecloc,
!$acc& polerecglob) async
!$acc end data
c
      call update_gang_rec(npolerecloc)
c
      npolelocloop = merge( npoleloc,
     &                    (int(npoleloc/16)+1)*16,
     &                    (mod(npoleloc,16).eq.0))
c
#ifdef _OPENACC
!$acc wait
!$acc host_data use_device(poleglob,polerecglob)
c     call thrust_sort(poleglob,npolebloc)
      call thrust_sort(polerecglob,npolerecloc,rec_stream)
!$acc end host_data
c!$acc parallel loop async(rec_queue)
c      do i = 1,npolebloc
c         poleloc(poleglob(i)) = i
c      end do
!$acc parallel loop present(polerecloc,polerecglob)
!$acc&         async(rec_queue)
      do i = 1,npolerecloc
         polerecloc(polerecglob(i)) = i
      end do
#endif
c
c     comput nZ_Onlyloc if necessary
c
      if (nZ_Onlyglob.ne.0) then
         call compute_nZ_Onlyaxe
      end if
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0.or.istep.lt.0) return
c
      call prmem_request(poleglobnl,nlocnl,async=.true.)
c
!$acc serial async(rec_queue) present(npolelocnl)
      npolelocnl = 0
!$acc end serial
!$acc parallel loop async(rec_queue)
!$acc&         present(npolelocnl,nbpole,ineignl,pollist,
!$acc&  polsiz,repart,poleglobnl,polelocnl,x,y,z)
      do i = 1, nlocnl
        iglob  = ineignl(i)
        iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
        if (iipole.eq.0) cycle
        polecount = nbpole(iglob)
        if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
          call distprocpart1(iglob,rank,d,.true.,x,y,z)
          if (repart(iglob).eq.rank) d = 0
          if (d*d.le.(mbuf2/4)) then
!$acc atomic capture
            npolelocnl = npolelocnl + 1
            npole_cap  = npolelocnl
!$acc end atomic
            poleglobnl(npole_cap) = polecount + 1
            polelocnl(polecount+1) = npole_cap
          end if
        end if
      end do
!$acc update host(npolelocnl) async(rec_queue)
#ifdef _OPENACC
!$acc wait(rec_queue)
!$acc host_data use_device(poleglobnl)
      call thrust_sort(poleglobnl,npolelocnl,rec_stream)
!$acc end host_data
!$acc parallel loop present(polelocnl,poleglobnl) async(rec_queue)
      do i = 1,npolelocnl
         polelocnl(poleglobnl(i)) = i
      end do
#endif

c
c     npolelocnlloop is npolelocnl if npolelocnl is a multiple of 16
c                    or the first one greater
c     npolelocnlloop = merge( npolelocnl,
c    &                     (int(npolelocnl/16)+1)*16,
c    &                     (mod(npolelocnl,16).eq.0))

c     npoleblocloop = merge( npolebloc,
c    &                    (int(npolebloc/16)+1)*16,
c    &                    (mod(npolebloc,16).eq.0))
      end

      subroutine kpolar_reset1(istep)
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
      use utilcomm,only:do_not_commpole
      use utilgpu
      implicit none
      integer,intent(in)::istep
      integer i,j,iipole,iglob,iproc
      integer polecount,npole_cap,modnl
      integer ipr,ibufbeg
      real(t_p) d

      do_not_commpole = .true.
      npoleloc = 0
!$acc data present(poleglob,poleloc,polerecglob,polerecloc)

      do i = 1, nloc
         iglob = glob(i)
         iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
         if (iipole.eq.0) cycle
         polecount = nbpole(iglob)
         if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
            npoleloc = npoleloc + 1
            poleglob(npoleloc) = polecount + 1
            poleloc(polecount+1) = npoleloc
         end if
      end do
      bufbegpole(rank+1) = 1
      domlenpole(rank+1) = npoleloc
      npolebloc = npoleloc

      do iproc = 1, n_recep1
        ipr = p_recep1(iproc)+1
        ibufbeg = bufbeg(ipr)
        if (domlen(ipr).ne.0) then
          bufbegpole(ipr) = npolebloc + 1
        else
          bufbegpole(ipr) = 1
        end if
        do i = 1, domlen(ipr)
          iglob  = glob(ibufbeg+i-1)
          iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
          if (iipole.eq.0) cycle
          polecount = nbpole(iglob)
          if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
            npolebloc = npolebloc + 1
            poleglob(npolebloc) = polecount + 1
            poleloc(polecount+1) = npolebloc
          end if
        end do
        if (domlen(ipr).ne.0) then
          domlenpole(ipr) = npolebloc-bufbegpole(ipr)+1
        else
          domlenpole(ipr) = 0
        end if
      end do
!$acc update device(poleglob,poleloc) async
c
      npolerecloc_old = npolerecloc
      npolerecloc  = 0
      do i = 1, nlocrec
         iglob  = globrec(i)
         iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
         if (iipole.eq.0) cycle
         polecount = nbpole(iglob)
         if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
            npolerecloc = npolerecloc + 1
            polerecglob(npolerecloc) = polecount + 1
            polerecloc(polecount+1) = npolerecloc
         end if
      end do
      domlenpolerec(rank+1) = npolerecloc
!$acc update device(polerecglob,polerecloc) async
c
      call update_gang_rec(npolerecloc)
c
      npolelocloop = merge( npoleloc,
     &                    (int(npoleloc/16)+1)*16,
     &                    (mod(npoleloc,16).eq.0))
c
!$acc end data
c
c     comput nZ_Onlyloc if necessary
c
      if (nZ_Onlyglob.ne.0 .and.nproc.ne.1) then
         call compute_nZ_Onlyaxe
      end if
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0.and.istep.ge.0) return
c
      call prmem_request(poleglobnl,nlocnl,async=.true.)
c     if (allocated(poleglobnl)) deallocate(poleglobnl)
c     allocate (poleglobnl(nlocnl))

!$acc data present(poleglobnl,polelocnl)
c
      npolelocnl = 0
      do i = 1, nlocnl
        iglob  = ineignl(i)
        iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
        if (iipole.eq.0) cycle
        polecount = nbpole(iglob)
        if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
          call distprocpart(iglob,rank,d,.true.)
          if (repart(iglob).eq.rank) d = 0
          if (d*d.le.(mbuf2/4)) then
            npolelocnl = npolelocnl + 1
            poleglobnl(npolelocnl) = polecount + 1
            polelocnl(polecount+1) = npolelocnl
          end if
        end if
      end do
!$acc update device(poleglobnl,polelocnl) async
!$acc end data
c
c     npolelocnlloop is npolelocnl if npolelocnl is a multiple of 16
c                    or the first one greater
c     npolelocnlloop = merge( npolelocnl,
c    &                     (int(npolelocnl/16)+1)*16,
c    &                     (mod(npolelocnl,16).eq.0))

c     npoleblocloop = merge( npolebloc,
c    &                    (int(npolebloc/16)+1)*16,
c    &                    (mod(npolebloc,16).eq.0))

      end

      subroutine compute_nZ_Onlyaxe
      use atmlst,only : poleglob,polerecglob
      use mpole, only : nZ_Onlyloc,npoleloc,npolerecloc,ipolaxe
      implicit none
      integer i,k
      integer:: dire=0,reci=0

!$acc data present(nZ_onlyloc,poleglob,polerecglob)
      dire = 0
      reci = 0
!$acc enter data copyin(dire,reci) async
!$acc parallel loop async present(dire)
      do i = 1,npoleloc
         k = poleglob(i)
         if (btest(ipolaxe(k),3)) dire = dire +1
      end do
!$acc parallel loop async present(reci)
      do i = 1,npolerecloc
         k = polerecglob(i)
         if (btest(ipolaxe(k),3)) reci = reci +1
      end do
!$acc serial async
      nZ_Onlyloc = max(dire,reci)
!$acc end serial
!$acc update host(nZ_Onlyloc) async
!$acc exit data delete(dire,reci) async
!$acc end data

      end subroutine

      ! Communicate to neigbhor process poleglob configuration
      ! This will simplify direct fiels communication when 
      ! poleglob is being construct parallel on device
      subroutine orderPole
      use atoms   ,only : n
      use atmlst  ,only : poleglob
      use domdec  ,only : nproc,ndir,rank,pbig_send,nbig_send,pbig_recep
     &            , nbig_recep,bufbeg,bufbegpole,domlen,domlenpole
     &            , glob,nbloc,comm_dir,COMM_TINKER
      use inform  ,only : deb_Path
      use mpole   ,only : poleloc,pollist,polsiz,nbpole,npoleloc
     &                  , npolebloc,nz_Onlyglob,npole
      use potent  ,only : use_pmecore
      use timestat,only : timer_enter,timer_exit,timer_eneig
     &            , quiet_timers
      use utilcomm,only : reqs_poleglob,reqr_poleglob
     &            , do_not_commpole
      use utilgpu ,only : rec_queue
      use mpi
      implicit none
      integer ierr,i,j,iproc,ipr,ibufbeg,idomlen
      integer iipole
      integer commloc,tag0,status(MPI_STATUS_SIZE)
      integer npolelocmpi(nproc)
      integer mpi_status1(MPI_STATUS_SIZE)
      parameter(tag0=0)

      if (ndir.eq.1.or.do_not_commpole) return
      if (use_pmecore.and.rank.ge.ndir) return

      if (use_pmecore) then
         commloc = comm_dir
      else
         commloc = COMM_TINKER
      end if
      call timer_enter( timer_eneig )
      if (deb_Path) write(*,'(3X,A)') '>>  orderPole'

      ! TODO Optimize this routine for multi-time step
      !  --- Short Range

      ! update switch
      do_not_commpole = .true.
c
c     Initiate poleglob sending
c
!$acc host_data use_device(poleglob)
      do i = 1, nbig_send
         call MPI_ISEND(poleglob(1),npoleloc,MPI_INT,pbig_send(i)
     &                 ,tag0,commloc,reqs_poleglob(i),ierr)
      end do

      if (n.eq.npole) then

        domlenpole(:) = domlen(:)
        bufbegpole(:) = bufbeg(:)
        npolebloc     = nbloc

      else

        ! Intercept communication size
        do i = 1, nbig_recep
           ipr = pbig_recep(i) + 1
           call MPI_Probe(pbig_recep(i), tag0, commloc,
     &                    mpi_status1, ierr)
           call MPI_Get_count( mpi_status1, MPI_INT,
     &                         domlenpole(ipr), ierr )
           ! Compute starting process index of
           ! CSR Storage
           if (domlen(ipr) .ne. 0) then
              if (i.eq.1) then
                 bufbegpole(ipr) = 1 + npoleloc
                 npolebloc = npolebloc + domlenpole(ipr)
              else
                 bufbegpole(ipr) = bufbegpole(pbig_recep(i-1)+1) +
     &                             domlenpole(pbig_recep(i-1)+1)
                 npolebloc = npolebloc + domlenpole(ipr)
              end if
           else
              bufbegpole(ipr) = 1
              domlenpole(ipr) = 0
           end if
        end do

      end if
c
c     Begin Message reception
c
      do i = 1, nbig_recep
         ipr     = pbig_recep(i) + 1
         idomlen = domlenpole(ipr)
         ibufbeg = bufbegpole(ipr)
         call MPI_IRECV(poleglob(ibufbeg),idomlen,MPI_INT,pbig_recep(i)
     &                 ,tag0,commloc,reqr_poleglob(i),ierr)
      end do
!$acc end host_data
c
c     Wait for poleglob communication request to end
c     And Fill the missing part of it
c
!$acc data present(poleglob,poleloc)
      do i = 1, nbig_recep
         ipr     = pbig_recep(i) + 1
         idomlen = domlenpole(ipr)
         ibufbeg = bufbegpole(ipr) - 1
         call MPI_WAIT(reqr_poleglob(i),status,ierr)
         ! Finish npoleloc construction
!$acc parallel loop async(rec_queue)
         do j = 1, idomlen
            iipole = poleglob(ibufbeg+j)
            poleloc(iipole) = ibufbeg+j
         end do
         call MPI_WAIT(reqs_poleglob(i),status,ierr)
      end do
!$acc end data

      !print*,rank,domlenpole,bufbegpole
      if (deb_Path) write(*,'(3X,A)') '<<  orderPole'
      call timer_exit ( timer_eneig,quiet_timers )
      end subroutine
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
      integer jj,kk,ierr
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
      if (hostrank.ne.0) goto 120
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
  120 call MPI_BARRIER(hostcomm,ierr)
c     call upload_device_polar(0)
      end

c
c     Update polar's data on device
c
      subroutine upload_device_polar(config)
      use domdec
      use inform ,only: deb_Path
      use mpole
      use mpi    ,only: MPI_BARRIER
      use polar
      use polgrp
      use sizes
      use tinMemory
      implicit none
      integer,intent(in)::config
      integer up_polargrp , up_other
      parameter( up_polargrp=0, up_other=1 )
      integer ierr
#ifdef _OPENACC
 12   format(2x,'upload_device_polar')
      if(deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif
      if (config.eq.up_polargrp) then
c!$acc update device(ip11(:,:),np11(:),ip12(:,:),np12(:),
c!$acc&              ip13(:,:),np13(:),ip14(:,:),np14(:))
      else
!$acc update device(thole,pdamp,polarity)

!$acc update device(xaxis,yaxis,zaxis)
!$acc update device(ipolaxe,pole,nbpole,ipole)
!$acc update device(pollist)
      end if

      end subroutine

      subroutine delete_data_polar
      use domdec
      use inform ,only: deb_Path
      use mpole
      use polar
      use polgrp
      use sizes
      use tinMemory
      implicit none

 12   format(2x,'delete_data_polar')
      if (deb_Path) print 12

      call shmem_request(polarity,winpolarity,[0],config=mhostacc)
      call shmem_request(thole,   winthole,   [0],config=mhostacc)
      call shmem_request(pdamp,   winpdamp,   [0],config=mhostacc)

      end subroutine
c
c     subroutine alloc_shared_polar : allocate shared memory pointers for polar
c     parameter arrays
c
      subroutine alloc_shared_polar
      use sizes
      use atoms
      use domdec
      use polar
      use mpi
      use tinMemory
      implicit none

      if (associated(thole).and.n.eq.size(thole)) return
c
      call shmem_request(polarity,winpolarity,[n],config=mhostacc)
      call shmem_request(thole,   winthole,   [n],config=mhostacc)
      call shmem_request(pdamp,   winpdamp,   [n],config=mhostacc)
      end
