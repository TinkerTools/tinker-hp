c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kurey  --  Urey-Bradley parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kurey" assigns the force constants and ideal distances
c     for the Urey-Bradley 1-3 interactions; also processes any
c     new or changed parameter values
c
c
#include "tinker_macro.h"
      subroutine kurey(init)
      use angle
      use atmlst
      use atmtyp
      use domdec
      use inform
      use iounit
      use keys
      use kurybr
      use potent
      use utils
      use urey
      use utilgpu  
      use urypot
      implicit none
      integer i,j,nu,nups,nutot,nuq
      integer ia,ib,ic
      integer ita,itb,itc
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      integer size,next,temp
      integer iangle,ureycount,nureyloc1,nureyloc_capture
      integer::isys=0
      integer*8 pt
      real(t_p) bb,tt
      logical header
      character*4 pa,pb,pc
      character*20 keyword
      character*240 record
      character*240 string
      logical init,max_reach
c
      if (init) then
c
c     process keywords containing Urey-Bradley parameters
c
        if(deb_Path) print*,'kurey'
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:9) .eq. 'UREYBRAD ') then
              ia = 0
              ib = 0
              ic = 0
              bb = 0.0_ti_p
              tt = 0.0_ti_p
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,bb,tt
   10         continue
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Urey-Bradley Parameters :',
     &                      //,5x,'Atom Classes',8x,'K(UB)',5x,
     &                         'Distance',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,bb,tt
   30            format (4x,3i4,2x,f12.3,f12.4)
              end if
c             size = 4
c             call numeral (ia,pa,size)
c             call numeral (ib,pb,size)
c             call numeral (ic,pc,size)
              if (ia .le. ic) then
                 call front_convert_base(ia,ib,ic,pt)
c                pt = pa//pb//pc
              else
                 call front_convert_base(ic,ib,ia,pt)
c                pt = pc//pb//pa
              end if
              max_reach =.true.
              do j = 1, maxnu
                 if (ku(j).eq. -1 .or. ku(j).eq.pt) then
                    ku(j)    = pt
                    ucon(j)  = bb
                    dst13(j) = tt
                    max_reach=.false.
                    exit
                 end if
              end do
              if (max_reach) then
                if (rank.eq.0) write (iout,*)
     &         ' KUREY  --  Too many Urey-Bradley',
     &                   ' Interaction Parameters'
                abort = .true.
              endif
           end if
        end do
c
c     process keywords containing Angle repulsion parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'ANGREP ') then
              ia = 0
              ib = 0
              ic = 0
              bb = 0.0d0
              tt = 0.0d0
              string = record(next:240)
              read (string,*)  ia,ib,ic,bb,tt
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,'(A,5x,A,8x,A,5x,A)')
     &                  ' Additional Angle repulsion Parameters :',
     &                     'Atom Classes','D','1/b'
                 end if
                 if (rank.eq.0) write (iout,'(4x,3i4,2x,f12.3,f12.4)')
     &               ia,ib,ic,bb,tt
              end if
              size = 4
c              call numeral (ia,pa,size)
c              call numeral (ib,pb,size)
c              call numeral (ic,pc,size)
              if (ia .le. ic) then
                 call front_convert_base(ia,ib,ic,pt)
c                pt = pa//pb//pc
              else
                 call front_convert_base(ic,ib,ia,pt)
c                pt = pc//pb//pa
              end if
              max_reach =.true.
              do j = 1, maxnups
                 if (kups(j).eq.-1 .or. kups(j).eq.pt) then
                    kups(j) = pt
                    uconps(j) = bb
                    dst13ps(j) = tt
                    max_reach=.false.
                    exit
                 end if
              end do
              if (max_reach) then
                if (rank.eq.0) write (iout,*)
     &              ' KUREY  --  Too many Angle repulsion',
     &                   ' Interaction Parameters'
                abort = .true.
              end if
           end if
        end do
c
c     process keywords containing Quartic Urey parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:12) .eq. 'UREYQUARTIC ') then
              ia = 0
              ib = 0
              ic = 0
              bb = 0.0d0
              tt = 0.0d0
              string = record(next:240)
              read (string,*)  ia,ib,ic,bb,tt
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,'(A,5x,A,8x,A,5x,A)')
     &                  ' Additional Urey Quartic Parameters :',
     &                     'Atom Classes','D','1/b'
                 end if
                 if (rank.eq.0) write (iout,'(4x,3i4,2x,f12.3,f12.4)')
     &               ia,ib,ic,bb,tt
              end if
              size = 4
c              call numeral (ia,pa,size)
c              call numeral (ib,pb,size)
c              call numeral (ic,pc,size)
              if (ia .le. ic) then
                 call front_convert_base(ia,ib,ic,pt)
c                pt = pa//pb//pc
              else
                 call front_convert_base(ic,ib,ia,pt)
c                pt = pc//pb//pa
              end if
              max_reach =.true.
              do j = 1, maxnuq
                 if (kuq(j).eq.-1 .or. kuq(j).eq.pt) then
                    kuq(j) = pt
                    uconq(j) = bb
                    dst13q(j) = tt
                    max_reach=.false.
                    exit
                 end if
              end do
              if (max_reach) then
                if (rank.eq.0) write (iout,*)
     &              ' KUREY  --  Too many Quartic Urey',
     &                   ' Interaction Parameters'
                abort = .true.
              end if
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        nu = maxnu
        nups = maxnups
        nuq = maxnuq
        do i = maxnu, 1, -1
           if (ku(i) .eq. -1)  nu = i - 1
        end do
        do i = maxnups, 1, -1
           if (kups(i) .eq. -1)  nups = i - 1
        end do
        do i = maxnups, 1, -1
           if (kuq(i) .eq. -1)  nuq = i - 1
        end do
        nutot = nu + nups + nuq

        if (nups > 0 .or. nuq > 0) then
          disable_fuse_bonded = .true.
        endif
c
c       assign the Urey-Bradley parameters for each angle
c
        nurey = 0
        if (nutot .ne. 0) then
           do i = 1, nangle
              ia = iang(1,i)
              ib = iang(2,i)
              ic = iang(3,i)
              nburey(i) = nurey
              ita = class(ia)
              itb = class(ib)
              itc = class(ic)
c             size = 4
c             call numeral (ita,pa,size)
c             call numeral (itb,pb,size)
c             call numeral (itc,pc,size)
              if (ita .le. itc) then
                 call front_convert_base(ita,itb,itc,pt)
c                pt = pa//pb//pc
              else
                 call front_convert_base(itc,itb,ita,pt)
c                pt = pc//pb//pa
              end if
              do j = 1, nu
                 if (ku(j) .eq. pt) then
                    nurey = nurey + 1
                    iury(1,nurey) = ia
                    iury(2,nurey) = ib
                    iury(3,nurey) = ic
                    uk(nurey) = ucon(j)
                    ul(nurey) = dst13(j)
                    ureytypI(nurey) = UREY_BRAD
                    ! ku_sys construction
                    if (.not.is_find8(ku_sys(1),isys,pt)) then
                       isys = isys + 1
                       ku_sys(isys) = pt
                    end if
                    goto 60
                 end if
              end do
              do j = 1, nups
                 if (kups(j) .eq. pt) then
                    nurey = nurey + 1
                    iury(1,nurey) = ia
                    iury(2,nurey) = ib
                    iury(3,nurey) = ic
                    ureytypI(nurey) = UREY_ANGREP

                    uk(nurey) = uconps(j)
                    ul(nurey) = dst13ps(j)
                    if (.not.is_find8(ku_sys(1),isys,pt)) then
                       isys = isys + 1
                       ku_sys(isys) = pt
                    end if
                    goto 60
                 end if
              end do
              do j = 1, nuq
                 if (kuq(j) .eq. pt) then
                    nurey = nurey + 1
                    iury(1,nurey) = ia
                    iury(2,nurey) = ib
                    iury(3,nurey) = ic
                    ureytypI(nurey) = UREY_QUARTIC

                    uk(nurey) = uconq(j)
                    ul(nurey) = dst13q(j)
                    if (.not.is_find8(ku_sys(1),isys,pt)) then
                       isys = isys + 1
                       ku_sys(isys) = pt
                    end if
                    goto 60
                 end if
              end do
   60         continue
           end do
           ku_sys(0) = isys
           !print*,'nu   ',nu,isys
        end if
c
c       turn off the Urey-Bradley potential if it is not used
c
        if (nurey .eq. 0) use_urey = .false.
c
c     Upload or delete Urey data to device
c
        if (use_urey) then
           call upload_device_kurey
        else
           call delete_data_kurey
           return
        end if
      end if

      nu = size_i8_to_i(ku_sys(0))
!Wait for nangleloc
!$acc wait
      call prmem_request(ureyglob, nangleloc, async=.false.)
      
!$acc data present(nureyloc,nangleloc)
!$acc serial async
      nureyloc = 0
!$acc end serial
!$acc parallel loop async 
#ifdef USE_NVSHMEM_CUDA
!$acc& present(angleglob, nburey, class, ureyglob)
#else
!$acc& present(angleglob, nburey, iang, class, ureyglob)
#endif
      do i = 1, nangleloc
        iangle    = angleglob(i)
        ureycount = nburey(iangle)
#ifdef USE_NVSHMEM_CUDA
        ipe       =     (iangle-1)/nangle_pe
        ind       = mod((iangle-1),nangle_pe) +1
        ia        = d_iang(ipe)%pel(1,ind)
        ib        = d_iang(ipe)%pel(2,ind)
        ic        = d_iang(ipe)%pel(3,ind)
#else
        ia        = iang(1,iangle)
        ib        = iang(2,iangle)
        ic        = iang(3,iangle)
#endif
        ita       = class(ia)
        itb       = class(ib)
        itc       = class(ic)
        call front_convert_base3(min(ita,itc),itb,max(ita,itc),pt)
        nureyloc1 = 0
        do j = 1, nu
           if (ku_sys(j) .eq. pt) then
!$acc atomic capture
              nureyloc  = nureyloc + 1
              nureyloc_capture = nureyloc
!$acc end atomic
              nureyloc1 = nureyloc1 + 1
              ureyglob(nureyloc_capture) = ureycount + nureyloc1
              exit
           end if
        end do
      end do
!$acc update host(nureyloc) async
!$acc end data

 80   continue
      end

      subroutine upload_device_kurey
      use domdec,only: rank,hostcomm
      use inform,only: deb_Path
      use kurybr
      use mpi   ,only: MPI_BARRIER
      use tinMemory
      use urey
      use urypot
      implicit none
      integer ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_kurey')
      if(deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif
!$acc update device(ul,uk,iury,nburey,ureytypI)
!$acc update device(ku,ku_sys)
!$acc enter data copyin(nureyloc)
      end subroutine

      subroutine delete_data_kurey
      use domdec,only:rank
      use inform,only: deb_Path
      use kurybr
      use tinMemory
      use urey
      use urypot
      implicit none

 12   format(2x,'delete_data_kurey')
      if(deb_Path) print 12

      call shmem_request(uk,    winuk,     [0],config=mhostacc)
      call shmem_request(ul,    winul,     [0],config=mhostacc)
      call shmem_request(iury,  winiury, [0,0],config=mhostacc)
      call shmem_request(nburey,winnburey, [0],config=mhostacc)
      call shmem_request(ureytypI,winureytypI, [0],config=mhostacc)
      end subroutine
