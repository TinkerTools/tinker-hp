#include "tinker_macro.h"
      submodule(mdstate) submdstate
      use atoms
      use argue
      use bath
      use boxes
      use cell
      use domdec     ,only: nproc,rank,nloc
      use energi     ,only: etot_ave,etot_std
      use files
      use keys
      use inform     ,only: verbose
      use mdstuf
      use moldyn
      use random_mod ,only: get_randseed,random
     &               ,samplevec,get_pickcount
#ifdef _OPENACC
     &               ,randomgpu,normalgpu
#endif
      use titles
      use utils
      use utilgpu
      implicit none
      integer i,ndat,n_save
      character*120 mddat(20)
      integer(8) npik,npikD,npikDr
      logical respa_l,respa1_l
      real(r_p) ,pointer:: ms_alt2(:)
      contains

      subroutine mds_checkdyn
      implicit none
      integer i,itmp,nc,s,next
      integer(8) ii,pack_,nn,nnD
      real(r_p) rtmp,valrandom
      logical ltmp
      character*64 keyword,string

 10   format('mds_check Issue detected with ',a
     &      ,'parameter (given/loaded)',2I14)
 11   format('mds_check Issue detected with ',a
     &      ,'parameter (given/loaded)',2F14.7)
 12   format('mds_check Issue detected with ',a
     &      ,'parameter (given/loaded)',2L4)
 13   format(32('-'))

      nc    = 0
      npik  = 0
      npikD = 0
      npikDr= 0

      do i = 1,ndat
         next=1
         call gettext( mddat(i),keyword,next )
         next=next+1
         if      ( keyword(1:3).eq.'N ') then
            read(mddat(i)(next:),*,err=66,end=66) itmp
            if (itmp.ne.n) then;
               print 10,'N ',n,itmp; nc=nc+1;
            end if

         else if ( keyword(1:7).eq.'RANDOM ') then
            read(mddat(i)(next:),'(I14,3I19)',err=66,end=66) itmp
     &          ,npik,npikD,npikDr
            if(itmp.ne.ms(1)%randseed) then
               print 10,'RANDOM ',ms(1)%randseed,itmp; nc=nc+1;
               !print*, npik,npikD
            end if

         else if ( keyword(1:12).eq.'ISOTHERMAL ') then
            read(mddat(i)(next:),'(L)',err=66,end=66) ltmp
            if (ltmp.neqv.isothermal) then
               print 12,'ISOTHERMAL ',isothermal,ltmp; nc=nc+1
            end if

         else if ( keyword(1:9).eq.'ISOBARIC ') then
            read(mddat(i)(next:),'(L)',err=66,end=66) ltmp
            if (ltmp.neqv.isobaric) then
               print 12,'ISOBARIC ',isobaric,ltmp; nc=nc+1;
            end if

         else if ( keyword(1:3).eq.'DT ') then
            read(mddat(i)(next:),*,err=66,end=66) rtmp
            if (abs(rtmp-ms(1)%dt).gt.1e-4) then
               print 11,'DT ',ms(1)%dt,rtmp; nc=nc+1;
            end if

         else if ( keyword(1:14).eq.'ETOT_MEAN_STD ') then
            read(mddat(i)(next:),*,err=66,end=66)
     &           etot_ave,etot_std

         else if ( keyword(1:7).eq.'KELVIN ') then
            read(mddat(i)(next:),*,err=66,end=66) rtmp
            if (abs(rtmp-kelvin).gt.1e-4) then
               print 11,'KELVIN ',kelvin,rtmp; nc=nc+1;
            end if

         else if ( keyword(1:9).eq.'PRESSURE ') then
            read(mddat(i)(next:),*,err=66,end=66) rtmp
            if (abs(rtmp-atmsph).gt.1e-4) then
               print 11,'PRESSURE ',kelvin,rtmp; nc=nc+1;
            end if
         end if
      end do

      if (nc.ne.0) then
         write(0,13)
         __TINKER_FATAL__
      end if

      if (npik.ne.0.or.npikD.ne.0) then
 16      format(" -- Aligning Random Generator to Save State --")
         if (rank.eq.0.and.verbose) print 16
#ifdef _OPENACC
      ! Align the device random generator to the saved state

      ! Random Device alignment
      s     = 3*n
      pack_ = npikDr/(s*20)
      do i = 1,20
         do ii = 1,pack_; call randomgpu(samplevec,s); end do
      end do
      pack_ = mod(npikDr,int(s*20,8))
      do i = 1,pack_/s; call randomgpu(samplevec,s); end do
      s     = mod(pack_,int(s,8))
      if (s.gt.0) call randomgpu(samplevec,s)

      ! Normal Device alignment
      s     = 3*n
      pack_ = npikD/(s*20)
      if (rank.eq.0.and.verbose) write(*,'(a)',advance='no') ' <'
      do i = 1,20
         do ii = 1,pack_; call normalgpu(samplevec,s); end do
         if (rank.eq.0.and.verbose) write(*,'(a)',advance='no') '|'
      end do
      pack_ = mod(npikD,int(s*20,8))
      do i = 1,pack_/s; call normalgpu(samplevec,s); end do
      s     = mod(pack_,int(s,8))
      if (s.gt.0) call normalgpu(samplevec,s)
      if (rank.eq.0.and.verbose) write(*,'(a2)',advance='no') '> '
#endif
      ! CPU Alignment
      pack_ = (npik)/20
      if (rank.eq.0.and.verbose) write(*,'(a)',advance='no') '<'
      ! Align the host random generator to the saved state
      do i = 1,20
         do ii = 1,pack_; itmp = random(); end do
         if (rank.eq.0.and.verbose) write(*,'(a)',advance='no') '|'
      end do
      pack_ = mod(npik,20)
      do ii = 1,pack_; itmp = random(); end do
      end if
!$acc wait
      if (rank.eq.0.and.verbose) write(*,'(a)') '>'

      !call get_pickcount(npik,npikD)
      !print*, npik,npikD

 66   continue

      end subroutine

      module subroutine mds_init
      implicit none
      integer,parameter:: id=1
      integer ios,unt,i,freeunit,s

 25   format(a120)
 10   format(" Code Erreur:",I3," Detected while reading ",A)

      call fetchkey('TRACK-MDSTATE',ms_back_p,-1)
      track_mds = (ms_back_p.gt.0)
      fw_mds    = existkey('FWRITE-MDSTATE')
      if (ms_back_p.lt.0.or.nproc.gt.1) return

      respa_l  = (index(integrate,'RESPA').gt.0)
      respa1_l = (index(integrate,'RESPA1').gt.0)
      s        = merge( 12,9, respa_l)
      s        = merge(s+3,s,respa1_l)
      n_save           = 0
      ms(1)%n          = n
      ms(1)%istep      = 0
      ms(1)%randseed   = get_randseed()
      ms(1)%isothermal = isothermal
      ms(1)%isobaric   = isobaric
      ms(1)%xbox       = xbox
      ms(1)%ybox       = ybox
      ms(1)%zbox       = zbox
      ms(1)%kelvin     = kelvin
      ms(1)%atmsph     = atmsph
      read(arg(3),*,iostat=ios) ms(1)%dt  ! Fetch timestep among arguments
      ms(1)%dt         = 1d-3*ms(1)%dt
      ms(1)%dumpdyn    = filename(1:leng)//'.mdx.dyn'
      ms(1)%dumpdat    = filename(1:leng)//'.mdx.dat'
      call version(ms(1)%dumpdyn,'old')
      call version(ms(1)%dumpdat,'old')
      inquire(file=ms(1)%dumpdyn,exist=ms(1)%lddyn)
      inquire(file=ms(1)%dumpdat,exist=ms(1)%lddat)
c
c     Compare MD actual State against saved State
c
      if (ms(1)%lddat) then
         unt = freeunit()
         open(unit=unt,file=ms(1)%dumpdat,status='old')
         rewind(unit=unt)
         ios = 0
         ndat= 0
         do while( ios.eq.0 )
            read(unt,25,iostat=ios) mddat(ndat+1)
            ndat = ndat+1
         end do
         if (ios.gt.0) then
            write(0,10) ios, trim(ms(1)%dumpdat)
            __TINKER_FATAL__
         end if
         close(unit=unt)
 36      format(/,' ---Tinker-HP: Resuming MD from a Save State')
         if (rank.eq.0.and.verbose) print 36
         call mds_checkdyn
      end if

      if (verbose.and.nproc.gt.1) then
 20      format(" MD State Tracking is not supported in parallel run")
         if(rank.eq.0) write(0,20)
         __TINKER_FATAL__
      end if

      if (.not.track_mds) return
      call version(ms(1)%dumpdyn,'new')
      call prmem_requestm(ms(1)%state,s*n)

      ms_x(1:s*n) => ms(id)%state(    1:)
      ms_y(1:n)   => ms(id)%state(  n+1:)
      ms_z(1:n)   => ms(id)%state(2*n+1:)
      ms_v(1:3*n) => ms(id)%state(3*n+1:)
      ms_a(1:3*n) => ms(id)%state(6*n+1:)
      if (respa_l)
     &   ms_alt(1:3*n)  => ms(id)%state(9*n+1:)
      if (respa1_l)
     &   ms_alt2(1:3*n) => ms(id)%state(12*n+1:)

      if (rank.eq.0.and.verbose) then
 46      format(' ---Tinker-HP: Backup MD State every ',I0,' steps')
         print 46, ms_back_p
      end if

      end subroutine

      module subroutine mds_save
      implicit none
      integer i,j

      n_save      = n_save + 1
      ms(1)%istep = step_c
      ms(1)%xbox  = xbox
      ms(1)%ybox  = ybox
      ms(1)%zbox  = zbox
      call get_pickcount(npik,npikD,npikDr)

      if (rank.eq.0.and.verbose.and.n_save.lt.5) then
 16      format(" -- Tinker-HP : Saving Actual MD State")
         print 16
      end if

!$acc parallel loop async default(present)
      do i = 1,n
         ms_x(i) = x(i)
         ms_y(i) = y(i)
         ms_z(i) = z(i)
         do j = 1,3
            ms_v(3*(i-1)+j) = v(j,i) 
            ms_a(3*(i-1)+j) = a(j,i) 
         end do
      end do
#if TINKER_SINGLE_PREC
      __TINKER_FATAL__
#else
      if( respa_l)call mem_move_r8(ms_alt ,aalt ,int(3*n,8),rec_stream)
      if(respa1_l)call mem_move_r8(ms_alt2,aalt2,int(3*n,8),rec_stream)
#endif
      end subroutine

      module subroutine mds_prt
      implicit none
      integer unt,idat,i,k,freeunit
      character*40 fstr
      character*2 atmc

      if (n_save.eq.0) return
!$acc update host(ms_x) async
      unt = freeunit ()

      if (rank.eq.0.and.verbose) then
 11      format(" ---Tinker-HP !! Dumping Last MD State saved in ",A)
         write(0,11) ms(1)%dumpdat
      end if
c
c     Dump State file
c
      if (ms(1)%lddat) then
         open   (unit=unt,file=ms(1)%dumpdat,status='old')
         rewind (unit=unt)
      else
         open   (unit=unt,file=ms(1)%dumpdat,status='new')
      end if
 200  format('-- THIS FILE IS NOT SUPPOSED TO BE MANUALLY EDITED --')
 20   format('N ',I0)
 21   format('RANDOM ',I14,3I19)
 22   format('STEP ',I0)
 23   format('ISOTHERMAL ',L1)
 24   format('ISOBARIC ',L1)
 25   format('DT ',F26.16)
 26   format('ETOT_MEAN_STD ',2F26.16)
 29   format('KELVIN ',F26.16)
 30   format('PRESSURE ',F26.8)
      write(unt,200)
      write(unt,20) ms(1)%n
      if (ms(1)%isothermal.or.ms(1)%isobaric) then
         write(unt,21) ms(1)%randseed, npik, npikD, npikDr
      end if
      write(unt,22) ms(1)%istep
      write(unt,23) ms(1)%isothermal
      write(unt,24) ms(1)%isobaric
      write(unt,25) ms(1)%dt
      write(unt,26) etot_ave,etot_std
      write(unt,29) ms(1)%kelvin
      write(unt,30) ms(1)%atmsph
      close(unit=unt)

      if (ms(1)%lddyn) then
         open (unit=unt,file=ms(1)%dumpdyn,status='old')
         rewind (unit=unt)
      else
         open (unit=unt,file=ms(1)%dumpdyn,status='new')
      end if
c
c     Dump dynfile
c
      fstr = '('' Number of Atoms and Title :'')'
      write (unt,fstr(1:32))
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (unt,fstr(1:4))  n
      else
         fstr = '('//atmc//',2x,a)'
         write (unt,fstr(1:9))  n,title(1:ltitle)
      end if
c
c     save the periodic box edge lengths and angles
c
      fstr = '('' Periodic Box Dimensions :'')'
      write (unt,fstr(1:30))
      fstr = '(3d26.16)'
      write (unt,fstr(1:9))  ms(1)%xbox,ms(1)%ybox,ms(1)%zbox
      write (unt,fstr(1:9))  alpha,beta,gamma
c
c     save the atomic positions, velocities and accelerations
c
!$acc wait
      fstr = '('' Current Atomic Positions :'')'
      write (unt,fstr(1:31))
      fstr = '(3d26.16)'
      do i = 1, n
         write (unt,fstr(1:9))  ms_x(i),ms_y(i),ms_z(i)
      end do
      fstr = '('' Current Atomic Velocities :'')'
      write (unt,fstr(1:32))
      fstr = '(3d26.16)'
      do i = 1, n
         k = 3*(i-1)
         write (unt,fstr(1:9))  ms_v(1+k),ms_v(2+k),ms_v(3+k)
      end do
      fstr =  '('' Current Atomic Accelerations :'')'
      write (unt,fstr(1:36))
      fstr = '(3d26.16)'
      do i = 1, n
         k = 3*(i-1)
         write (unt,fstr(1:9))  ms_a(1+k),ms_a(2+k),ms_a(3+k)
      end do
      if (respa_l) then
      fstr =  '('' Alternate Atomic Accelerations :'')'
      write (unt,fstr(1:38))
      fstr = '(3d26.16)'
      do i = 1, n
         k = 3*(i-1)
         write (unt,fstr(1:9))  ms_alt(1+k),ms_alt(2+k),ms_alt(3+k)
      end do
      end if
      if (respa1_l) then
      fstr =  '('' Alternate 2 Atomic Accelerations :'')'
      write (unt,fstr(1:40))
      fstr = '(3d26.16)'
      do i = 1, n
         k = 3*(i-1)
         write (unt,fstr(1:9))  ms_alt2(1+k),ms_alt2(2+k),ms_alt2(3+k)
      end do
      end if
      close (unit=unt)
      end subroutine

      end submodule
