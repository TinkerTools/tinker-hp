c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine kimprop  --  improper dihedral parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "kimprop" assigns potential parameters to each improper
c     dihedral in the structure and processes any changed values
c
c
#include "tinker_macro.h"
      subroutine kimprop(init)
      use atmlst
      use atmtyp
      use atoms
      use couple
      use domdec
      use improp
      use inform
      use iounit
      use keys
      use kiprop
      use potent
#ifdef _OPENACC
      use thrust
      use utilgpu  ,only:rec_stream
#endif
      use utils
      use tinheader ,only:ti_p,re_p
      use tinMemory
      use tors
      implicit none
      integer i,j,k
      integer ia,ib,ic,id,ie
      integer ita,itb,itc,itd,ite
      integer iglob,impropcount,niproploc1
      integer isize,next,capt
      integer,save:: ndis,ndi
      real(t_p) tk,tv,symm
      logical header,done
      integer pa,pb,pc,pd
      character*12 zeros
      character*16 blank
      integer(8) pt0,pt1,pti
      integer(8) pt2,pt3
      integer(8) pt(6)
      character*20 keyword
      character*240 record
      character*240 string
      logical init
c
      blank = '                '
      zeros = '000000000000'
      if (init) then
         ndis = 0
c
c     process keywords containing improper dihedral parameters
c
        if(deb_Path) write(*,*) 'kimprop init'

        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:9) .eq. 'IMPROPER ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              tk = 0.0_ti_p
              tv = 0.0_ti_p
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,id,tk,tv
   10         continue
              !isize = 4
              !call numeral (ia,pa,isize)
              !call numeral (ib,pb,isize)
              !call numeral (ic,pc,isize)
              !call numeral (id,pd,isize)
              !pti = pa//pb//pc//pd
              call front_convert_base(0,ia,ib,ic,id,pti)
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Improper Dihedral',
     &                         ' Parameters :',
     &                      //,5x,'Atom Classes',20x,'K(ID)',
     &                         7x,'Angle',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,id,tk,tv
   30            format (4x,4i4,10x,2f12.3)
              end if
              do j = 1, maxndi
                 if (kdi(j).eq.-1 .or. kdi(j).eq.pti) then
                    kdi(j) = pti
                    dcon(j) = tk
                    tdi(j) = tv
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KIMPROP  --  Too many Improper Dihedral',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        ndi = maxndi
        do i = maxndi, 1, -1
           if (kdi(i) .eq. -1)  ndi = i - 1
        end do
c
c       assign improper dihedral parameters for each improper angle;
c       multiple symmetrical parameters are given partial weights
c
        niprop = 0
        if (ndi .ne. 0) then
           do i = 1, n
              if (n12(i) .eq. 3) then
                 ia = i
                 ib = i12(1,i)
                 ic = i12(2,i)
                 id = i12(3,i)
                 nbimprop(i) = niprop
                 ita = class(ia)
                 itb = class(ib)
                 itc = class(ic)
                 itd = class(id)
                 !isize = 4
                 !call numeral (ita,pa,isize)
                 !call numeral (itb,pb,isize)
                 !call numeral (itc,pc,isize)
                 !call numeral (itd,pd,isize)
                 !pt(1) = pa//pb//pc//pd
                 !pt(2) = pa//pb//pd//pc
                 !pt(3) = pa//pc//pb//pd
                 !pt(4) = pa//pc//pd//pb
                 !pt(5) = pa//pd//pb//pc
                 !pt(6) = pa//pd//pc//pb
                 call front_convert_base(0,ita,itb,itc,itd,pt(1))
                 call front_convert_base(0,ita,itb,itd,itc,pt(2))
                 call front_convert_base(0,ita,itc,itb,itd,pt(3))
                 call front_convert_base(0,ita,itc,itd,itb,pt(4))
                 call front_convert_base(0,ita,itd,itb,itc,pt(5))
                 call front_convert_base(0,ita,itd,itc,itd,pt(6))
                 !pt3 = pa//pb//zeros//zeros
                 !pt2 = pa//pc//zeros//zeros
                 !pt1 = pa//pd//zeros//zeros
                 !pt0 = pa//zeros
                 call front_convert_base(0,ita,itb,0,0,pt3)
                 call front_convert_base(0,ita,itc,0,0,pt2)
                 call front_convert_base(0,ita,itd,0,0,pt1)
                 call front_convert_base(0,ita,  0,0,0,pt0)
                 symm = 1.0_ti_p
                 if (itb.eq.itc .or. itb.eq.itd .or. itc.eq.itd)
     &              symm = 2.0_ti_p
                 if (itb.eq.itc .and. itb.eq.itd .and. itc.eq.itd)
     &              symm = 6.0_ti_p
                 done = .false.
                 do j = 1, ndi
                    call back_convert_base(ite,pa,itb,itc,itd,kdi(j))
                    if (ita .eq. pa) then
                       do k = 1, 6
                          if (kdi(j) .eq. pt(k)) then
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ic
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = id
                             else if (k .eq. 4) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 5) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             else if (k .eq. 6) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = ib
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                             done = .true.
                             ! Build system kdi
                             if (.not.is_find8(kdi_sys(1),ndis,kdi(j)))
     &                          then
                                ndis = ndis + 1
                                kdi_sys(ndis) = kdi(j)
                             end if
                          end if
                       end do
                    end if
                 end do
                 if (.not. done) then
                    do j = 1, ndi
                       if (kdi(j) .eq. pt1) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          done = .true.
                          ! Build system kdi
                          if(.not.is_find8(kdi_sys(1),ndis,pt1))then
                             ndis = ndis + 1
                             kdi_sys(ndis) = kdi(j)
                          end if
                       else if (kdi(j) .eq. pt2) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          done = .true.
                          ! Build system kdi
                          if(.not.is_find8(kdi_sys(1),ndis,pt2))then
                             ndis = ndis + 1
                             kdi_sys(ndis) = kdi(j)
                          end if
                       else if (kdi(j) .eq. pt3) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          done = .true.
                          ! Build system kdi
                          if(.not.is_find8(kdi_sys(1),ndis,pt3))then
                             ndis = ndis + 1
                             kdi_sys(ndis) = kdi(j)
                          end if
                       end if
                    end do
                 end if
                 if (.not. done) then
                    do j = 1, ndi
                       if (kdi(j) .eq. pt0) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          ! Build system kdi
                          if(.not.is_find8(kdi_sys(1),ndis,pt0))then
                             ndis = ndis + 1
                             kdi_sys(ndis) = kdi(j)
                          end if
                       end if
                    end do
                 end if
              end if
           end do
        end if
c
c       turn off the improper dihedral potential if it is not used
c
        if (niprop .eq. 0)  use_improp = .false.

        if (.not.use_improp) then
           call delete_data_kimprop
           !!! Exit ---------------
           return
        end if

        call update_device_kimprop
c
c       determine the total number of forcefield parameters
c
        ndi = maxndi
        do i = maxndi, 1, -1
           if (kdi(i).eq.-1)  ndi = i - 1
        end do

        ! store system forcefield parameter's number
        kdi_sys(0) = ndis
        if (deb_Path) write(*,'(A,2I6)') '  ndi',ndi,ndis

      end if

      call prmem_request(impropglob,6*nbloc,async=.true.)
!$acc serial async present(niproploc)
      niproploc = 0
!$acc end serial
      if (ndi .ne. 0) then

!$acc parallel loop gang vector async private(pt)
!$acc&         present(glob,nbimprop,n12,i12,class,impropglob,kdi_sys)
!$acc&         present(niproploc)
         do i = 1, nloc
            iglob       = glob(i)
            impropcount = nbimprop(iglob)
            if (n12(iglob) .eq. 3) then
               ia  = iglob
               ib  = i12(1,iglob)
               ic  = i12(2,iglob)
               id  = i12(3,iglob)
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
               call front_convert_base(0,ita,itb,itc,itd,pt(1))
               call front_convert_base(0,ita,itb,itd,itc,pt(2))
               call front_convert_base(0,ita,itc,itb,itd,pt(3))
               call front_convert_base(0,ita,itc,itd,itb,pt(4))
               call front_convert_base(0,ita,itd,itb,itc,pt(5))
               call front_convert_base(0,ita,itd,itc,itd,pt(6))
               call front_convert_base(0,ita,itb,0,0,pt3)
               call front_convert_base(0,ita,itc,0,0,pt2)
               call front_convert_base(0,ita,itd,0,0,pt1)
               call front_convert_base(0,ita,  0,0,0,pt0)
               done = .false.
               niproploc1 = 0
               do j = 1, ndis
                  call back_convert_base(ite,pa,itb,itc,itd,kdi_sys(j))
                  if (ita .eq. pa) then
                     do k = 1, 6
                        if (kdi_sys(j) .eq. pt(k)) then
!$acc atomic capture
                           niproploc  = niproploc + 1
                           capt       = niproploc
!$acc end atomic
                           niproploc1 = niproploc1 + 1
                           impropglob(capt)=impropcount +niproploc1
                           done = .true.
                        end if
                     end do
                  end if
               end do
               if (.not.done) then
                  do j = 1, ndis
                     if (kdi_sys(j) .eq. pt1) then
                        do k = 1, 3
!$acc atomic capture
                           niproploc  = niproploc + 1
                           capt       = niproploc
!$acc end atomic
                           niproploc1 = niproploc1 + 1
                           impropglob(capt)=impropcount +niproploc1
                        end do
                        done = .true.
                     else if (kdi_sys(j) .eq. pt2) then
                        do k = 1, 3
!$acc atomic capture
                           niproploc  = niproploc + 1
                           capt       = niproploc
!$acc end atomic
                           niproploc1 = niproploc1 + 1
                           impropglob(capt)=impropcount +niproploc1
                        end do
                        done = .true.
                     else if (kdi_sys(j) .eq. pt3) then
                        do k = 1, 3
!$acc atomic capture
                           niproploc  = niproploc + 1
                           capt       = niproploc
!$acc end atomic
                           niproploc1 = niproploc1 + 1
                           impropglob(capt)=impropcount +niproploc1
                        end do
                        done = .true.
                     end if
                  end do
               end if
               if (.not. done) then
                  do j = 1, ndi
                     if (kdi_sys(j) .eq. pt0) then
                        do k = 1, 3
!$acc atomic capture
                           niproploc  = niproploc + 1
                           capt       = niproploc
!$acc end atomic
                           niproploc1 = niproploc1 + 1
                           impropglob(capt)=impropcount +niproploc1
                        end do
                     end if
                  end do
               end if
            end if
         end do
!$acc update host(niproploc) async
      end if
c
      end

      subroutine update_device_kimprop
      use domdec,only: rank,hostcomm
      use improp
      use inform
      use kiprop
      use mpi   ,only: MPI_BARRIER
      use tinMemory
      implicit none
      integer ierr

#ifdef _OPENACC
12    format(2x,'update_device_kimprop')
      if (deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif
!$acc update device(kdi_sys)
!$acc update device(iiprop,nbimprop,kprop,vprop)
!$acc enter data create(niproploc)
      end subroutine

      subroutine delete_data_kimprop
      use improp
      use inform
      use kiprop
      use tinMemory
      implicit none

12    format(2x,'delete_data_kimprop')
      if (deb_Path) print 12
      call shmem_request(nbimprop,winnbimprop,  [0],config=mhostacc)
      call shmem_request(   kprop,   winkprop,  [0],config=mhostacc)
      call shmem_request(   vprop,   winvprop,  [0],config=mhostacc)
      call shmem_request(  iiprop,  winiiprop,[0,0],config=mhostacc)
      end subroutine
