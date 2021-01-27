c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kopbend  --  out-of-plane bending parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kopbend" assigns the force constants for out-of-plane bends
c     at trigonal centers via Wilson-Decius-Cross or Allinger angles;
c     also processes any new or changed parameter values
c
c
#include "tinker_precision.h"
      subroutine kopbend(init)
      use angle
      use angpot
      use atmlst
      use atmtyp
      use couple
      use domdec
      use fields
      use inform
      use iounit
      use keys
      use kopbnd
      use opbend
      use potent
      use utils
      use usage
      use utilgpu
      implicit none
      integer i,j,it,iangle
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      integer nopbendloc1,opbendcount,nopbendloc_capture
      integer nopb,size,temp
      integer next,number
      integer:: isys=0
      integer*8 pt,pt0,pt1
      real(t_p) fopb
      logical header,done
c      logical, allocatable :: jopb(:)
      character*4 pa,pb,pc,pd
      character*4 zero4
      character*8 zero8
      character*16 blank
      character*20 keyword
      character*240 record
      character*240 string
      logical init
c
c     blank = '                '
c     zero4 = '0000'
c     zero8 = '00000000'
      if (init) then
c
c
c     process keywords containing out-of-plane bend parameters
c
        if(deb_Path) print*,'kopbend init'
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'OPBEND ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              fopb = 0.0_ti_p
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,id,fopb
   10         continue
c             size = 4
c             call numeral (ia,pa,size)
c             call numeral (ib,pb,size)
c             call numeral (ic,pc,size)
c             call numeral (id,pd,size)
              if (ic .le. id) then
                 call front_convert_base(0,ia,ib,ic,id,pt)
c                pt = pa//pb//pc//pd
              else
                 call front_convert_base(0,ia,ib,id,ic,pt)
c                pt = pa//pb//pd//pc
              end if
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Out-of-Plane Bend',
     &                         ' Parameters :',
     &                      //,5x,'Atom Classes',19x,'K(OPB)',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,id,fopb
   30            format (4x,4i4,10x,f12.3)
              end if
              size = 4
              do j = 1, maxnopb
                 if (kopb(j).eq.-1 .or. kopb(j).eq.pt) then
                    kopb(j) = pt
                    opbn(j) = fopb
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KOPBEND --  Too many Out-of-Plane',
     &                   ' Angle Bending Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       use special out-of-plane bend parameter assignment for MMFF
c
        if (forcefield .eq. 'MMFF94') then
           call kopbendm
           call upload_device_kopbend(0)
           return
        end if
c
c     determine the total number of forcefield parameters
c
        nopb = maxnopb
        do i = maxnopb, 1, -1
           if (kopb(i) .eq. -1)  nopb = i - 1
        end do
c
c       perform dynamic allocation of some local arrays
c
        call prmem_request(jopb,maxclass)
c
c       make list of atom classes using out-of-plane bending
c
        do i = 1, maxclass
           jopb(i) = .false.
        end do
        do i = 1, maxnopb
           if (kopb(i) .eq. -1)  goto 60
           call back_convert_base(ita,itb,it,itc,itd,kopb(i))
c          it = number(kopb(i)(5:8))
           jopb(it) = .true.
        end do
   60   continue
c
c       allocate arrays
c
        call alloc_shared_opbend
c
c       assign out-of-plane bending parameters for each angle
c
        nopbend = 0
        if (nopb .ne. 0) then
           header = .true.
           do i = 1, nangle
              ib = iang(2,i)
              nbopbend(i) = nopbend
              itb = class(ib)
              if (jopb(itb) .and. n12(ib).eq.3) then
                 ia = iang(1,i)
                 ita = class(ia)
                 ic = iang(3,i)
                 itc = class(ic)
                 id = iang(4,i)
                 itd = class(id)
c                size = 4
c                call numeral (ita,pa,size)
c                call numeral (itb,pb,size)
c                call numeral (itc,pc,size)
c                call numeral (itd,pd,size)
                 if (ita .le. itc) then
                    call front_convert_base(0,itd,itb,ita,itc,pt)
c                   pt = pd//pb//pa//pc
                 else
                    call front_convert_base(0,itd,itb,itc,ita,pt)
c                   pt = pd//pb//pc//pa
                 end if
c                pt1 = pd//pb//zero8
c                pt0 = zero4//pb//zero8
                 call front_convert_base(0,itd,itb,0,0,pt1)
                 call front_convert_base(0,  0,itb,0,0,pt0)
                 done = .false.
                 do j = 1, nopb
                    if (kopb(j).eq.pt .or.kopb(j).eq.pt1.or.
     &                  kopb(j).eq.pt0) then
                       nopbend = nopbend + 1
                       iopb(nopbend) = i
                       opbk(nopbend) = opbn(j)
                       done = .true.
                       if (.not.is_find8(kopb_sys(1),isys,kopb(j))) then
                          isys = isys + 1
                          kopb_sys(isys) = kopb(j)
                       end if
                       goto 70
                    end if
                 end do
   70            continue
                 if (use_opbend .and. .not.done) then
                    if (use(ia) .or. use(ib) .or. use(ic) .or. use(id))
     &                 abort = .true.
                    if (header) then
                       header = .false.
                       if (rank.eq.0) write (iout,80)
   80                  format (/,' Undefined Out-of-Plane Bend',
     &                            ' Parameters :',
     &                         //,' Type',24x,'Atom Names',24x,
     &                            'Atom Classes',/)
                    end if
                   if (rank.eq.0) write (iout,90)  id,name(id),ib,
     &                              name(ib),ia,name(ia),
     &                              ic,name(ic),itd,itb,ita,itc
   90               format (' Angle-OP',3x,4(i6,'-',a3),5x,4i5)
                 end if
              else
                 iang(4,i) = ib
              end if
           end do
           kopb_sys(0) = isys
           if (deb_Path) print*,'nopb ',nopb,isys
           call upload_device_kopbend(0)
        end if
cc
cc       perform deallocation of some local arrays
cc
c        deallocate (jopb)
c
c
c       turn off the out-of-plane bending term if it is not used
c
        if (nopbend .eq. 0)  use_opbend = .false.
        call upload_device_kopbend(1)
      end if

      nopb = size_i8_to_i(kopb_sys(0))
!Wait for nangleloc
!$acc wait
      call prmem_request(opbendglob,nangleloc,async=.false.)

!$acc data present(nopbendloc,nangleloc)
!$acc serial async
      nopbendloc = 0
!$acc end serial
!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc& present(angleglob, nbopbend, class, jopb, n12)
#else
!$acc& present(angleglob, nbopbend, iang, class, jopb, n12)
#endif
      do i = 1, nangleloc
         iangle      = angleglob(i)
         opbendcount = nbopbend(iangle)
#ifdef USE_NVSHMEM_CUDA
         ipe  =     (iangle-1)/nangle_pe
         ind  = mod((iangle-1),nangle_pe) +1
         ib   = d_iang(ipe)%pel(2,ind)
#else
         ib   = iang(2,iangle)
#endif
         itb  = class(ib)
         if (jopb(itb) .and. n12(ib).eq.3) then
#ifdef USE_NVSHMEM_CUDA
            ia   = d_iang(ipe)%pel(1,ind)
            ic   = d_iang(ipe)%pel(3,ind)
            id   = d_iang(ipe)%pel(4,ind)
#else
            ia   = iang(1,iangle)
            ic   = iang(3,iangle)
            id   = iang(4,iangle)
#endif
            ita  = class(ia)
            itc  = class(ic)
            itd  = class(id)
            if (ita .gt. itc) then
               temp = ita; ita = itc; itc = temp;
            end if
            call front_convert_base(0,itd,itb,ita,itc,pt )
            call front_convert_base(0,itd,itb,  0,  0,pt1)
            call front_convert_base(0,  0,itb,  0,  0,pt0)
            nopbendloc1 = 0
            do j = 1, kopb_sys(0)
               if (kopb_sys(j).eq.pt  .or.
     &             kopb_sys(j).eq.pt1 .or.
     &             kopb_sys(j).eq.pt0) then
!$acc atomic capture
                  nopbendloc  = nopbendloc  + 1
                  nopbendloc_capture  = nopbendloc
!$acc end atomic
                  nopbendloc1 = nopbendloc1 + 1
                  opbendglob(nopbendloc_capture) = opbendcount 
     &                  + nopbendloc1
                  exit
               end if
            end do
         end if
      end do
!$acc update host(nopbendloc) async
!$acc end data
      end

      subroutine upload_device_kopbend(config)
      use angle
      use angpot
      use domdec,only:rank,hostcomm
      use fields
      use kopbnd
      use mpi   ,only: MPI_BARRIER
      use nvshmem
      use opbend
      use sizes,only:tinkerdebug
      use tinMemory
      implicit none
      integer,intent(in)::config
      integer c_size,ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_kopbend')
      if(rank.eq.0.and.tinkerdebug) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif

      if (config.eq.0) then
!$acc update device(iopb,opbk,nbopbend)
         if (forcefield .eq. 'MMFF94') return
#ifdef USE_NVSHMEM_CUDA
      c_size = size(c_iang(mype)%pel)
      call shmem_update_device(iang,size(iang),
     &     dst=c_iang(mype)%pel,nd=c_size,config=mnvshonly)
#else
!$acc update device(iang)
#endif
!$acc update device(kopb_sys)
      else if (config.eq.1) then
!$acc update device(angtyp)
!$acc update device(jopb)
!$acc enter data copyin(nopbendloc)
      else
 14      format('WARNING ! upload_device_kopbend ! Unknown',
     &          ' configuration',I4)
         print 14, config
      end if
      end subroutine

      subroutine delete_data_kopbend
      use angle
      use angpot
      use domdec,only:rank
      use fields
      use kopbnd
      use opbend
      use sizes,only:tinkerdebug
      use tinMemory
      implicit none

 12   format(2x,'delete_data_kopbend')
      if(rank.eq.0.and.tinkerdebug) print 12
      call shmem_request(opbk,winopbk,[0],config=mhostacc)
      call shmem_request(iopb,winiopb,[0],config=mhostacc)
      if (forcefield .eq. 'MMFF94') return
!$acc exit data delete(nopbendloc)
      end subroutine

c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine kopbendm  --  MMFF out-of-plane bend parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "kopbendm" assigns the force constants for out-of-plane bends
c     according to the Merck Molecular Force Field (MMFF)
c
c
      subroutine kopbendm
      use angle
      use atmtyp
      use atoms
      use kopbnd
      use merck
      use opbend
      use utils
      implicit none
      integer i,j,m
      integer nopb,size
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer itta,ittb
      integer ittc,ittd
      integer:: isys=0
      integer*8 pt
      character*4 pa,pb,pc,pd
      character*16 blank
c
c
c     determine the total number of forcefield parameters
c
      blank = '                '
      nopb = maxnopb
      do i = maxnopb, 1, -1
         if (kopb(i) .eq. -1)  nopb = i - 1
      end do
c
c     assign MMFF out-of-plane bending parameter values
c
      nopbend = 0
      if (nopb .ne. 0) then
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            id = iang(4,i)
            itta = type(ia)
            ittb = type(ib)
            ittc = type(ic)
            ittd = type(id)
            m = 0
   10       continue
            m = m + 1
            if (m .eq. 1) then
               ita = eqclass(itta,1)
               itb = eqclass(ittb,1)
               itc = eqclass(ittc,1)
               itd = eqclass(ittd,1)
            else if (m .eq. 2) then
               ita = eqclass(itta,2)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,2)
               itd = eqclass(ittd,2)
            else if (m .eq. 3) then
               ita = eqclass(itta,3)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,3)
               itd = eqclass(ittd,3)
            else if (m .eq. 4) then
               ita = eqclass(itta,4)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,4)
               itd = eqclass(ittd,4)
            else if (m .eq. 5) then
               ita = eqclass(itta,5)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,5)
               itd = eqclass(ittd,5)
            end if
            if (ia.ne.0 .and. ib.ne.0 .and. ic.ne.0 .and. id.ne.0) then
               if (m .gt. 5) then
                  nopbend = nopbend + 1
                  iopb(nopbend) = i
                  opbk(nopbend) = 0.0
               else
c                 size = 4
c                 call numeral (ita,pa,size)
c                 call numeral (itb,pb,size)
c                 call numeral (itc,pc,size)
c                 call numeral (itd,pd,size)
                  if (itd.le.ita .and. itd.le.itc) then
                     if (ita .le. itc) then
                        call front_convert_base(0,itd,itb,ita,itc,pt)
c                       pt = pd//pb//pa//pc
                     else
                        call front_convert_base(0,itd,itb,itc,ita,pt)
c                       pt = pd//pb//pc//pa
                     end if
                  else if (ita.le.itc .and. ita.le.itd) then
                     if (itd .le. itc) then
                        call front_convert_base(0,ita,itb,itd,itc,pt)
c                       pt = pa//pb//pd//pc
                     else
                        call front_convert_base(0,ita,itb,itc,itd,pt)
c                       pt = pa//pb//pc//pd
                     end if
                  else if (itc.le.ita .and. itc.le.itd) then
                     if (ita .le. itd) then
                        call front_convert_base(0,itc,itb,ita,itd,pt)
c                       pt = pc//pb//pa//pd
                     else
                        call front_convert_base(0,itc,itb,itd,ita,pt)
c                       pt = pc//pb//pd//pa
                     end if
                  end if
                  do j = 1, nopb
                     if (kopb(j) .eq. pt) then
                        nopbend = nopbend + 1
                        iopb(nopbend) = i
                        opbk(nopbend) = opbn(j)
                        if (.not.is_find8(kopb_sys(1),isys,pt)) then
                           isys = isys + 1
                           kopb_sys(isys) = pt
                        end if
                        goto 20
                     end if
                  end do
                  if (class(ib).eq.8 .or. class(ib).eq.17 .or.
     &                class(ib).eq.26 .or. class(ib).eq.43 .or.
     &                class(ib).eq.49 .or. class(ib).eq.73 .or.
     &                class(ib).eq.82) then
                     nopbend = nopbend + 1
                     iopb(nopbend) = i
                     opbk(nopbend) = 0.0
                     goto 20
                  end if
                  goto 10
   20             continue
               end if
            end if
         end do
         kopb_sys(0) = isys
         !print*,'nopb ',nopb,isys
!$acc update device (kopb_sys) async
      end if

c
c     turn off the out-of-plane bending term if it is not used
c
c      if (nopbend .eq. 0)  use_opbend = .false.
      return
      end
c
c     subroutine alloc_shared_opbend : allocate shared memory pointers for opbend
c     parameter arrays
c
      subroutine alloc_shared_opbend
      use angle  ,only: nangle
      use opbend
      use tinMemory
      implicit none

      if (associated(opbk).and.nangle.eq.size(opbk)) return

      call shmem_request(opbk,winopbk,[nangle],config=mhostacc)
      call shmem_request(iopb,winiopb,[nangle],config=mhostacc)
      end
