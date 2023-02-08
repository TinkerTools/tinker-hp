c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine kimptor  --  improper torsion parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "kimptor" assigns torsional parameters to each improper
c     torsion in the structure and processes any changed values
c
c
#include "tinker_macro.h"
      subroutine kimptor(init)
      use atmlst
      use atmtyp
      use atoms
      use couple
      use domdec
      use keys
      use imptor
      use inform
      use iounit
      use kitors
      use math
      use potent
      use tinheader ,only:ti_p,re_p
      use tinMemory
      use tors
      use utils
      implicit none
      integer i,j,k,nti,ntis
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer bta,btb,btc,btd,bte
      integer iglob,imptorcount,nitorsloc1
      integer isize,next,capt
      integer ft(3)
      real(t_p) angle,symm
      real(t_p) vt(3),st(3)
      logical header,done
c     character*4 pa,pb,pc,pd
      integer,parameter:: zero=0
      integer(8) pt0,pt1,pti
      integer(8) pt2,pt3
      integer(8) pt(6)
      character*20 keyword
      character*240 record
      character*240 string
      logical init
c
      if (init) then
c
c     process keywords containing improper torsion parameters
c
        if(rank.eq.0.and.tinkerdebug) print*,'kimptor init'
        header = .true.
        ntis = 0

        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:8) .eq. 'IMPTORS ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              do j = 1, 3
                 vt(j) = 0.0_ti_p
                 st(j) = 0.0_ti_p
                 ft(j) = 0
              end do
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,3)
   10         continue
c             isize = 4
c             call numeral (ia,pa,isize)
c             call numeral (ib,pb,isize)
c             call numeral (ic,pc,isize)
c             call numeral (id,pd,isize)
c             pti = pa//pb//pc//pd
              call front_convert_base(0,ia,ib,ic,id,pti)
              call torphase (ft,vt,st)
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20             format (/,' Additional Improper Torsion Parameters :',
     &                    //,5x,'Atom Classes',15x,'1-Fold',12x,
     &                       '2-Fold',12x,'3-Fold',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,id,
     &               (vt(j),st(j),j=1,3)
   30            format (4x,4i4,2x,3(f11.3,f7.1))
              end if
              do j = 1, maxnti
                 if (kti(j).eq.-1 .or. kti(j).eq.pti) then
                    kti(j) = pti
                    ti1(1,j) = vt(1)
                    ti1(2,j) = st(1)
                    ti2(1,j) = vt(2)
                    ti2(2,j) = st(2)
                    ti3(1,j) = vt(3)
                    ti3(2,j) = st(3)
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KIMPTOR  --  Too many Improper Torsion',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        nti = maxnti
        do i = maxnti, 1, -1
           if (kti(i) .eq. -1)  nti = i - 1
        end do
c
c       assign improper torsional parameters for each improper torsion;
c       multiple symmetrical parameters are given partial weights
c
        nitors = 0
        if (nti .ne. 0) then
           do i = 1, n
              if (n12(i) .eq. 3) then
                 ia = i12(1,i)
                 ib = i12(2,i)
                 ic = i
                 id = i12(3,i)
                 nbimptor(i) = nitors
                 ita = class(ia)
                 itb = class(ib)
                 itc = class(ic)
                 itd = class(id)
c                isize = 4
c                call numeral (ita,pa,isize)
c                call numeral (itb,pb,isize)
c                call numeral (itc,pc,isize)
c                call numeral (itd,pd,isize)
c                pt(1) = pa//pb//pc//pd
c                pt(2) = pb//pa//pc//pd
c                pt(3) = pa//pd//pc//pb
c                pt(4) = pd//pa//pc//pb
c                pt(5) = pb//pd//pc//pa
c                pt(6) = pd//pb//pc//pa
                 call front_convert_base(0,ita,itb,itc,itd,pt(1))
                 call front_convert_base(0,itb,ita,itc,itd,pt(2))
                 call front_convert_base(0,ita,itd,itc,itb,pt(3))
                 call front_convert_base(0,itd,ita,itc,itb,pt(4))
                 call front_convert_base(0,itb,itd,itc,ita,pt(5))
                 call front_convert_base(0,itd,itb,itc,ita,pt(6))
                 call front_convert_base(0,zero,zero,itc,itd,pt3)
                 call front_convert_base(0,zero,zero,itc,itb,pt2)
                 call front_convert_base(0,zero,zero,itc,ita,pt1)
                 call front_convert_base(0,zero,zero,itc,zero,pt0)
                 symm = 1.0_ti_p
                 if (ita.eq.itb .or. ita.eq.itd .or. itb.eq.itd)
     &              symm = 2.0_ti_p
                 if (ita.eq.itb .and. ita.eq.itd .and. itb.eq.itd)
     &              symm = 6.0_ti_p
                 done = .false.
                 do j = 1, nti
                    call back_convert_base(bte,bta,btb,btc,btd,kti(j))
                    if ( btc .eq. itc) then
                       do k = 1, 6
                          if (kti(j) .eq. pt(k)) then
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = id
                             else if (k .eq. 3) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ib
                             else if (k .eq. 4) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             else if (k .eq. 5) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 6) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = ia
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                             done = .true.
                             ! Build system kti
                             if (.not.is_find8(kti_sys(1),ntis,kti(j)))
     &                          then
                                ntis = ntis + 1
                                kti_sys(ntis) = kti(j)
                             end if
                          end if
                       end do
                    end if
                 end do
                 if (.not. done) then
                    do j = 1, nti
                       if (kti(j) .eq. pt1) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          done = .true.
                          ! Build system kti
                          if (.not.is_find8(kti_sys(1),ntis,kti(j)))
     &                       then
                             ntis = ntis + 1
                             kti_sys(ntis) = kti(j)
                          end if
                       else if (kti(j) .eq. pt2) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          done = .true.
                          ! Build system kti
                          if (.not.is_find8(kti_sys(1),ntis,kti(j)))
     &                       then
                             ntis = ntis + 1
                             kti_sys(ntis) = kti(j)
                          end if
                       else if (kti(j) .eq. pt3) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          done = .true.
                          ! Build system kti
                          if (.not.is_find8(kti_sys(1),ntis,kti(j)))
     &                       then
                             ntis = ntis + 1
                             kti_sys(ntis) = kti(j)
                          end if
                       end if
                    end do
                 end if
                 if (.not. done) then
                    do j = 1, nti
                       if (kti(j) .eq. pt0) then
                          symm = 3.0_ti_p
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          ! Build system kti
                          if (.not.is_find8(kti_sys(1),ntis,kti(j)))
     &                       then
                             ntis = ntis + 1
                             kti_sys(ntis) = kti(j)
                          end if
                       end if
                    end do
                 end if
              end if
           end do
        end if
c
c       find the cosine and sine of the phase angle for each torsion
c
        do i = 1, nitors
           angle = itors1(2,i) / radian
           itors1(3,i) = cos(angle)
           itors1(4,i) = sin(angle)
           angle = itors2(2,i) / radian
           itors2(3,i) = cos(angle)
           itors2(4,i) = sin(angle)
           angle = itors3(2,i) / radian
           itors3(3,i) = cos(angle)
           itors3(4,i) = sin(angle)
        end do
c
c       turn off the improper torsional potential if it is not used
c
        if (nitors .eq. 0)  use_imptor = .false.

        if (use_imptor) then
           call update_device_kimptor
        else
           call delete_data_kimptor
           return
        end if

        ! store system forcefield parameter's number
        kti_sys(0) = ntis
        if (deb_Path) write(*,'(A,2I6)') '  nti',nti,ntis

      end if
      ntis = kti_sys(0)

      call prmem_request(imptorglob,6*nbloc,config=mhostacc)
!$acc serial async present(nitorsloc)
      nitorsloc = 0
!$acc end serial
      if (nti .ne. 0) then

!$acc parallel loop gang vector private(pt)
!$acc&         default(present) present(nitorsloc)
         do i = 1, nloc
            iglob = glob(i)
            imptorcount = nbimptor(iglob)
            nitorsloc1 = 0
            if (n12(iglob) .eq. 3) then
               ia = i12(1,iglob)
               ib = i12(2,iglob)
               ic = iglob
               id = i12(3,iglob)
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
c              isize = 4
c              call numeral (ita,pa,isize)
c              call numeral (itb,pb,isize)
c              call numeral (itc,pc,isize)
c              call numeral (itd,pd,isize)
               call front_convert_base(0,ita,itb,itc,itd,pt(1))
               call front_convert_base(0,itb,ita,itc,itd,pt(2))
               call front_convert_base(0,ita,itd,itc,itb,pt(3))
               call front_convert_base(0,itd,ita,itc,itb,pt(4))
               call front_convert_base(0,itb,itd,itc,ita,pt(5))
               call front_convert_base(0,itd,itb,itc,ita,pt(6))
               call front_convert_base(0,zero,zero,itc,itd,pt3)
               call front_convert_base(0,zero,zero,itc,itb,pt2)
               call front_convert_base(0,zero,zero,itc,ita,pt1)
               call front_convert_base(0,zero,zero,itc,0  ,pt0)
               symm = 1.0_ti_p
               done = .false.
               do j = 1, ntis
                  call back_convert_base(bte,bta,btb,btc,btd,kti_sys(j))
                  if (btc .eq. itc) then
                     do k = 1, 6
                        if (kti_sys(j) .eq. pt(k)) then
!$acc atomic capture
                           nitorsloc = nitorsloc + 1
                           capt      = nitorsloc
!$acc end atomic
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(capt)=imptorcount +nitorsloc1
                           done = .true.
                        end if
                     end do
                  end if
               end do
               if (.not. done) then
                  do j = 1, ntis
                     if (kti_sys(j) .eq. pt1) then
                        do k = 1, 3
!$acc atomic capture
                           nitorsloc = nitorsloc + 1
                           capt      = nitorsloc
!$acc end atomic
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(capt)=imptorcount +nitorsloc1
                           done = .true.
                        end do
                        done = .true.
                     else if (kti_sys(j) .eq. pt2) then
                        do k = 1, 3
!$acc atomic capture
                           nitorsloc = nitorsloc + 1
                           capt      = nitorsloc
!$acc end atomic
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(capt)=imptorcount +nitorsloc1
                        end do
                        done = .true.
                     else if (kti_sys(j) .eq. pt3) then
                        do k = 1, 3
!$acc atomic capture
                           nitorsloc = nitorsloc + 1
                           capt      = nitorsloc
!$acc end atomic
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(capt)=imptorcount +nitorsloc1
                        end do
                        done = .true.
                     end if
                  end do
               end if
               if (.not. done) then
                  do j = 1, ntis
                     if (kti_sys(j) .eq. pt0) then
                        do k = 1, 3
!$acc atomic capture
                           nitorsloc = nitorsloc + 1
                           capt      = nitorsloc
!$acc end atomic
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(capt)=imptorcount +nitorsloc1
                        end do
                     end if
                  end do
               end if
            end if
         end do
      end if
!$acc update host(nitorsloc) async
      end

      subroutine update_device_kimptor
      use domdec,only: rank,hostcomm
      use imptor
      use inform
      use kitors
      use mpi   ,only: MPI_BARRIER
      use tinMemory
      implicit none
      integer ierr

#ifdef _OPENACC
12    format(2x,'update_device_kimptor')
      if (deb_Path) write(*,12)
      call MPI_BARRIER(hostcomm,ierr)
#endif
!$acc update device(kti_sys)
!$acc update device(iitors,itors1,itors2,itors3,nbimptor)
!$acc enter data create(nitorsloc)
      end subroutine

      subroutine delete_data_kimptor
      use imptor
      use inform
      use kitors
      use tinMemory
      implicit none

12    format(2x,'delete_data_kimptor')
      if (deb_Path) print 12
      call shmem_request(nbimptor,winnbimptor,[0],config=mhostacc)
      call shmem_request(itors1  ,winitors1,[0,0],config=mhostacc)
      call shmem_request(itors2  ,winitors2,[0,0],config=mhostacc)
      call shmem_request(itors3  ,winitors3,[0,0],config=mhostacc)
      call shmem_request(iitors  ,winiitors,[0,0],config=mhostacc)
      end subroutine
