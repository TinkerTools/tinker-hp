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
      use chgpen
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
      integer i,j,k,ii,it,maxdr
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
      character*240 record
      character*240 string
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
           else if (keyword(1:11) .eq. 'INTEGRATOR ') then
              call getword (record,integrate,next)
              call upcase (integrate)
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
              pol = zero
              thl = zero
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
        do i = 1, n
           polarity(i) = 0.0_ti_p
           thole(i)    = 0.0_ti_p
           dirdamp(i)  = 0.0_ti_p
           pdamp(i)    = 0.0_ti_p
           it = type(i)
           if (it .ne. 0) then
              polarity(i) = polr(it)
              thole(i) = athl(it)
              dirdamp(i) = ddir(it)
              if (thole(i) .eq. 0.0_ti_p) then
                pdamp(i) = 0.0_ti_p
              else
                pdamp(i) = polarity(i)**sixth
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
        if ((use_polar .or. use_repuls)) then
          npole  = 0
          npolar = 0
          ncp    = 0
          ipole  = 0
          do i = 1, n
             if (polsiz(i).ne.0 .or. polarity(i).ne.0.0_ti_p) then
                nbpole(i) = npole
                npole     = npole + 1
                ipole(npole) = i
                pollist(i)   = npole
                zaxis (npole)= zaxis(i)
                xaxis (npole)= xaxis(i)
                yaxis (npole)= yaxis(i)
                polaxe(npole)= polaxe(i)
                do k = 1, maxpole
                   pole(k,npole) = pole(k,i)
                end do
                mono0(npole) = pole(1,i)
                if (polarity(i) .ne. 0.0_ti_p)  npolar = npolar + 1
                if (dirdamp(i) .ne. 0.0_ti_p)  use_dirdamp = .true.
                polarity(npole) = polarity(i)
                thole   (npole) = thole(i)
                dirdamp (npole) = dirdamp(i)
                pdamp   (npole) = pdamp(i)
                if (palpha(i) .ne. 0.0_ti_p)  ncp = ncp + 1
                pcore (npole) = pcore(i)
                pval  (npole) = pval(i)
                pval0 (npole) = pval(i)
                palpha(npole) = palpha(i)
             end if
          end do
        end if
        call FromPolaxe2Ipolaxe ! Brodcast 'nZ_onlyglob' after this call
c
c       Display if some poles does not have polarity
c
        npolar_ne_npole = (npolar.ne.npole)
        if (npolar.ne.npole) write(*,99) npole, npolar
 99        format("There are some poles without polarity :",/,
     &            3x,I12, "poles for ", I12, "polar atoms.")

 100    call MPI_BARRIER(hostcomm,ierr)
        call MPI_BCAST(npole, 1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(npolar,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(ncp,   1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(integrate,  11,   MPI_CHAR,0,hostcomm,ierr)
        call MPI_BCAST(nZ_Onlyglob, 1,    MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(use_dirdamp, 1,MPI_LOGICAL,0,hostcomm,ierr)
c
c       test multipoles at chiral sites and invert if necessary
c
        if (use_polar.and..not.use_chgtrn) call chkpole(.true.)
c
c       copy original polarizability values that won't change during mutation
c
        if (use_lambdadyn) then
           polarity_orig(:) = polarity(:)
        end if
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
           use_mlist = .false.
!$acc update device(use_mpole)
        endif
        if (npolar .eq. 0) then
           use_polar = .false.
!$acc update device(use_polar)
        end if
        if (ncp .ne. 0)  use_chgpen = .true.
        if (ncp .ne. 0)  use_thole  = .false.
        if (use_dirdamp) use_thole  = .true.
c
c       initialization for TCG and omega fit
c
        if ((polalg.eq.3).and.tcgpeek) then 
           !TODO 1.2  Ask louis about this value for single prec
           poleps = 0.00000001
        end if
!$acc enter data copyin(nZ_onlyloc)

        if (use_polar) then
           call upload_device_polar(1)
        else
           return
        end if
c
c  allocate predictor arrays
c
        if (use_polar.and.use_pred) then
           if (deb_Path) then
              call mem_get(m1,m2)
           end if
           lalt   = maxualt+1
           lshalt = maxualt+1
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
           call prmem_request(udalt,3,n,maxualt)
           call prmem_request(upalt,3,n,maxualt)
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
              call prmem_request(udshortalt,3,n,maxualt)
              call prmem_request(upshortalt,3,n,maxualt)
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
           do k = 1, maxualt
              do i = 1, npole
                 do j = 1, 3
                    udalt(j,i,k) = 0.0_ti_p
                    upalt(j,i,k) = 0.0_ti_p
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
              do k = 1, maxualt
                 do j = 1, 3
                    do i = 1, npole
                       upshortalt(j,i,k) = 0.0_ti_p
                       udshortalt(j,i,k) = 0.0_ti_p
                    end do
                 end do
              end do
#endif
           end if
        end if
c
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
      use utilcomm,only:do_not_commpole,no_commdir
#ifdef _OPENACC
      use thrust
#endif
      use utilgpu
      implicit none
      integer,intent(in)::istep
      integer i,j,iipole,iglob,iproc
      integer polecount,npole_cap,modnl
      integer ipr,ibufbeg
      integer npolerecloc2
      real(t_p) d,distcut2
!$acc routine(distprocpart1)

!$acc data present(glob,pollist,poleglob,poleloc,nbpole,bufbegpole,
!$acc&  domlenpole,domlenpolerec,polerecglob,poleloc,polerecloc)
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
          bufbegpole(ipr) = npolebloc + 1
        else
          bufbegpole(ipr) = 1
        end if
        do i = 1, domlen(ipr)
          iglob  = glob(ibufbeg+i-1)
          iipole = pollist(iglob)
          !skip atom if it is not in the multipole list
          if (iipole.eq.0) cycle
          polecount = nbpole(iglob)
          if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0) then
            npolebloc = npolebloc + 1
            npole_cap = npolebloc
            poleglob(npole_cap) = polecount + 1
            poleloc(polecount+1) = npole_cap
          end if
        end do
        if (domlen(ipr).ne.0) then
          domlenpole(ipr) = npolebloc-bufbegpole(ipr)+1
        else
          domlenpole(ipr) = 0
        end if
      end do
#else
      do_not_commpole = .false.
#endif
c
!$acc serial async(rec_queue)
      npolerecloc_old = npolerecloc
      npolerecloc  = 0
!$acc end serial
!$acc parallel loop async(rec_queue)
      do i = 1, nlocrec
         iglob  = globrec(i)
         iipole = pollist(iglob)
 
         !skip atom if it is not in the multipole list
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
!$acc parallel loop async(rec_queue)
      do i = nlocrec+1, nlocrec2
         iglob = globrec(i)
         iipole = pollist(iglob)
 
         !skip atom if it is not in the multipole list
         if (iipole.eq.0) cycle
         polecount = nbpole(iglob)
         if (polsiz(iglob) .ne. 0 .or. polarity(iipole).ne.0.0d0) then
!$acc atomic capture
            npolerecloc = npolerecloc + 1
            npole_cap   = npolerecloc
!$acc end atomic
            polerecglob(npole_cap) = polecount + 1
            polerecloc(polecount+1) = npole_cap
         end if
      end do
c
!$acc serial async(rec_queue)
      !npolerecloc2= npolerecloc
      npolerecloc = domlenpolerec(rank+1)
!$acc end serial
c
!$acc update host(domlenpole,domlenpolerec,bufbegpole
!$acc&    ,npoleloc,npolebloc,npolerecloc
!$acc&    ,polerecglob) async(rec_queue)
!$acc end data
c
c     npolelocloop = merge( npoleloc,
c    &                    (int(npoleloc/16)+1)*16,
c    &                    (mod(npoleloc,16).eq.0))
c
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0.or.istep.lt.0) return
c
      call prmem_request(poleglobnl,nlocnl,async=.true.)
c
!$acc serial async(rec_queue) present(npolelocnl)
      npolelocnl = 0
!$acc end serial
      if (no_commdir) then
         distcut2 = mbuf2
      else
         distcut2 = mbuf2/4
      end if

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
          if (d*d.le.distcut2) then
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
c
c     npolelocnlloop is npolelocnl if npolelocnl is a multiple of 16
c                    or the first one greater
c     npolelocnlloop = merge( npolelocnl,
c    &                     (int(npolelocnl/16)+1)*16,
c    &                     (mod(npolelocnl,16).eq.0))
c     npoleblocloop  = merge( npolebloc,
c    &                     (int(npolebloc/16)+1)*16,
c    &                     (mod(npolebloc,16).eq.0))
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
!$acc update device(polerecglob,polerecloc) async
c
      call update_gang_rec(npolerecloc)
c
      npolelocloop = merge( npoleloc,
     &                    (int(npoleloc/16)+1)*16,
     &                    (mod(npoleloc,16).eq.0))
c
c     comput nZ_Onlyloc if necessary
c
      if (nZ_Onlyglob.ne.0) call compute_nZ_Onlyaxe
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0.and.istep.ge.0) return
c
      call prmem_request(poleglobnl,nlocnl,async=.true.)
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
      use domdec,only : rank,Bdecomp1d
      use mpole ,only : nZ_Onlyloc,npolebloc,npolerecloc,ipolaxe
     &          ,ipole,Ax_Z_Only
      use pme   ,only : GridDecomp1d
      implicit none
      integer i,k
      integer looplen,reci

!$acc data create(reci) present(nZ_onlyloc) async

!$acc serial async
      nZ_Onlyloc = 0
      reci       = 0
!$acc end serial
!$acc parallel loop async default(present)
      do i = 1,npolebloc
         k = poleglob(i)
         if (ipolaxe(k).eq.Ax_Z_Only) nZ_Onlyloc = nZ_Onlyloc + 1
      end do
      if (.not.(Bdecomp1d.and.GridDecomp1d)) then
!$acc parallel loop async default(present)
         do i = 1,npolerecloc
            k = polerecglob(i)
            if (ipolaxe(k).eq.Ax_Z_Only) reci = reci+1
         end do
      end if
!$acc serial async
      nZ_Onlyloc = max(reci,nZ_Onlyloc)
!$acc end serial
!$acc update host(nZ_Onlyloc) async
!$acc wait
! -- wait for download
! -- to be used in ~rotpolegpu~ which is 2 calls away

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
      use utilgpu ,only : rec_queue,dir_queue
      use sizes   ,only : tinkerdebug
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

c
c     comput nZ_Onlyloc if necessary
c
      if (nZ_Onlyglob.ne.0) then
         call compute_nZ_Onlyaxe
      end if

      !print*,rank,domlenpole,bufbegpole
      if (btest(tinkerdebug,0)) call MPI_BARRIER(commloc,ierr)
      if (deb_Path)  write(*,'(3X,A)') '<<  orderPole'

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
      parameter (maxkeep=maxp11*maxvalue)
      parameter (maxlist=max(maxp12*maxp12,maxp12*maxp13))
      integer i,j,k
      integer it,jt
      integer jj,kk,ierr
      integer start,stop
      integer nlist,nkeep
      integer keep(maxkeep)
      integer list(maxlist)
      integer, allocatable :: mask(:)
      logical done
      integer max11,max12,max13,max14
      integer min11,min12,min13,min14
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

      if (tinkerdebug.gt.0) then
         max11=0;max12=0;max13=0;max14=0
         min11=0;min12=0;min13=0;min14=0
         do i = 1,n
            max11 = max(max11,np11(i))
            min11 = min(min11,np11(i))
            max12 = max(max12,np12(i))
            min12 = min(min12,np12(i))
            max13 = max(max13,np13(i))
            min13 = min(min13,np13(i))
            max14 = max(max14,np14(i))
            min14 = min(min14,np14(i))
         end do
         if (rank.eq.0) then
 71   format(10x,'min',8x,'max')
 72   format(A9,I4,7X,I4)
            print 71
            print 72, 'np11',min11,max11
            print 72, 'np12',min12,max12
            print 72, 'np13',min13,max13
            print 72, 'np14',min14,max14
         end if
      end if
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
      use chgpen
      use domdec
      use inform ,only: deb_Path
      use mpole
      use mpi    ,only: MPI_BARRIER
      use polar
      use polgrp
      use potent ,only: use_lambdadyn
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
!$acc update device(thole,pdamp,dirdamp,polarity)
      if (use_lambdadyn) then
!$acc update device(polarity_orig)
      end if

!$acc update device(xaxis,yaxis,zaxis
!$acc&       ,ipolaxe,ipole,pole,pollist,nbpole
!$acc&       ,mono0,pcore,pval,pval0,palpha)
      end if

      end subroutine

      subroutine dealloc_shared_polar
      use domdec
      use inform ,only: deb_Path
      use mpole
      use polar
      use polgrp
      use potent ,only: use_lambdadyn
      use sizes
      use tinMemory
      implicit none

 12   format(2x,'dealloc_shared_polar')
      if (deb_Path) print 12

      call shmem_request(polarity,winpolarity,[0],config=mhostacc)
      if (use_lambdadyn) then
         call shmem_request(polarity_orig,winpolarity_orig,[0]
     &                     ,config=mhostacc)
      end if
      call shmem_request(thole,   winthole,   [0],config=mhostacc)
      call shmem_request(pdamp,   winpdamp,   [0],config=mhostacc)
      call shmem_request(dirdamp, windirdamp, [0],config=mhostacc)

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
      use potent ,only: use_lambdadyn
      use mpi
      use tinMemory
      implicit none

      if (associated(thole).and.n.eq.size(thole)) return
c
      call shmem_request(polarity,winpolarity,[n],config=mhostacc)
      if (use_lambdadyn) then
         call shmem_request(polarity_orig,winpolarity_orig,[n]
     &                     ,config=mhostacc)
      end if
      call shmem_request(thole,   winthole,   [n],config=mhostacc)
      call shmem_request(pdamp,   winpdamp,   [n],config=mhostacc)
      call shmem_request(dirdamp, windirdamp, [n],config=mhostacc)
      end
