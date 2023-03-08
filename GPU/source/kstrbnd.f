c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kstrbnd  --  assign stretch-bend parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kstrbnd" assigns parameters for stretch-bend interactions
c     and processes new or changed parameter values
c
c
#include "tinker_macro.h"
      subroutine kstrbnd(init)
      use angle
      use atmlst
      use atmtyp
      use couple
      use domdec
      use fields
      use inform
      use iounit
      use keys
      use kstbnd
      use potent
      use strbnd
      use tinheader
      use utils
      use utilgpu
      implicit none
      integer i,j,k,nsb
      integer ia,ib,ic
      integer ita,itb,itc
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      integer nba,nbc
      integer next
      integer iangle,strbndcount,nstrbndloc1
      integer*8 ipt
      integer::isys=0
      real(t_p) sb1,sb2,temp
      logical header
      character*4 pa,pb,pc
      character*12 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
      logical init
      integer nstrbndloc_capture
c
      blank = '            '
      if (init) then
c
c     process keywords containing stretch-bend parameters
c
        if(deb_Path) print*,'kstrbnd init'
        isys   = 0
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'STRBND ') then
              ia = 0
              ib = 0
              ic = 0
              sb1 = 0.0_ti_p
              sb2 = 0.0_ti_p
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,sb1,sb2
   10         continue
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Stretch-Bend Parameters :',
     &                      //,5x,'Atom Classes',6x,'K(SB)-1',5x,
     &                         'K(SB)-2',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,sb1,sb2
   30            format (4x,3i4,2x,2f12.3)
              end if
c             size = 4
c             call numeral (ia,pa,size)
c             call numeral (ib,pb,size)
c             call numeral (ic,pc,size)
              if (ia .le. ic) then
                 call front_convert_base(ia,ib,ic,ipt)
c                pt = pa//pb//pc
              else
                 call front_convert_base(ic,ib,ia,ipt)
c                pt = pc//pb//pa
                 temp = sb1; sb1 = sb2; sb2 = temp
              end if
              do j = 1, maxnsb
                 if (ksb(j).eq.-1 .or. ksb(j).eq.ipt) then
                    ksb(j) = ipt
                    stbn(1,j) = sb1
                    stbn(2,j) = sb2
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KSTRBND  --  Too many Stretch-Bend',
     &                   ' Interaction Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        nsb = maxnsb
        do i = maxnsb, 1, -1
           if (ksb(i) .eq. -1)  nsb = i - 1
        end do
c
c       use special stretch-bend parameter assignment method for MMFF
c
        if (forcefield .eq. 'MMFF94') then
           call kstrbndm
           goto 80  !Exit routine
        end if
c
c     assign the stretch-bend parameters for each angle
c
        nstrbnd = 0
        if (nsb .ne. 0) then
           do i = 1, nangle
              ia = iang(1,i)
              ib = iang(2,i)
              ic = iang(3,i)
              nbstrbnd(i) = nstrbnd
              ita = class(ia)
              itb = class(ib)
              itc = class(ic)
c             size = 4
c             call numeral (ita,pa,size)
c             call numeral (itb,pb,size)
c             call numeral (itc,pc,size)
              if (ita .le. itc) then
                 call front_convert_base(ita,itb,itc,ipt)
c                pt = pa//pb//pc
              else
                 call front_convert_base(itc,itb,ita,ipt)
c                pt = pc//pb//pa
              end if
              do j = 1, nsb
                 if (ksb(j) .eq. ipt) then
                    nstrbnd = nstrbnd + 1
                    do k = 1, n12(ib)
                       if (i12(k,ib) .eq. ia)  nba = bndlist(k,ib)
                       if (i12(k,ib) .eq. ic)  nbc = bndlist(k,ib)
                    end do
                    isb(1,nstrbnd) = i
                    isb(2,nstrbnd) = nba
                    isb(3,nstrbnd) = nbc
                    if (ita .le. itc) then
                       sbk(1,nstrbnd) = stbn(1,j)
                       sbk(2,nstrbnd) = stbn(2,j)
                    else
                       sbk(1,nstrbnd) = stbn(2,j)
                       sbk(2,nstrbnd) = stbn(1,j)
                    end if
                    if (.not.is_find8(ksb_sys(1),isys,ipt)) then
                       isys          = isys + 1
                       ksb_sys(isys) = ipt
                    end if
                    goto 60
                 end if
              end do
   60         continue
           end do
           ksb_sys(0) = isys
        end if
        !print*,'nsb  ',nsb,isys
c
c       turn off the stretch-bend potential if it is not used
c
        if (nstrbnd .eq. 0) use_strbnd = .false.
c
c       Update strbnd shared memory data on device
c
        if (use_strbnd) then
           call upload_device_kstrbnd
        else
           call delete_data_kstrbnd
           return
        end if
      end if

      nsb = size_i8_to_i(ksb_sys(0))
! Wait for nangleloc
!$acc wait
      call prmem_request(strbndglob, nangleloc, async=.false.)
!$acc data present(nstrbndloc,nangleloc)
!$acc serial async
      nstrbndloc = 0
!$acc end serial
!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc& present(angleglob, nbstrbnd, class, ksb_sys)
#else
!$acc& present(angleglob, nbstrbnd, iang, class, ksb_sys)
#endif
      do i = 1, nangleloc
        iangle      = angleglob(i)
        strbndcount = nbstrbnd(iangle)
#ifdef USE_NVSHMEM_CUDA
        ipe         =     (iangle-1)/nangle_pe
        ind         = mod((iangle-1),nangle_pe) +1
        ia          = d_iang(ipe)%pel(1,ind)
        ib          = d_iang(ipe)%pel(2,ind)
        ic          = d_iang(ipe)%pel(3,ind)
#else               
        ia          = iang(1,iangle)
        ib          = iang(2,iangle)
        ic          = iang(3,iangle)
#endif              
        ita         = class(ia)
        itb         = class(ib)
        itc         = class(ic)
        call front_convert_base3(min(ita,itc),itb,max(ita,itc),ipt)
        nstrbndloc1 = 0
!$acc loop seq    
        do j = 1, nsb
           if (ksb_sys(j) .eq. ipt) then
!$acc atomic capture
              nstrbndloc = nstrbndloc + 1
              nstrbndloc_capture = nstrbndloc
!$acc end atomic
              nstrbndloc1 = nstrbndloc1 + 1
              strbndglob(nstrbndloc_capture) = strbndcount + nstrbndloc1
              exit
           end if
        end do
      end do
!$acc update host(nstrbndloc) async
!$acc end data
!$acc wait
80    continue      
      end subroutine

      subroutine upload_device_kstrbnd
      use domdec,only: rank,hostcomm
      use fields,only:forcefield
      use inform,only:deb_Path
      use kstbnd
      use mpi   ,only: MPI_BARRIER
      use sizes ,only:tinkerdebug
      use strbnd
      implicit none
      integer ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_kstrbnd')
      if(deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif

!$acc update device(isb,sbk)
      if (forcefield .eq. 'MMFF94') return
!$acc update device(nbstrbnd)
!$acc update device(ksb_sys)
!$acc enter data copyin(nstrbndloc)
      end subroutine

      subroutine delete_data_kstrbnd
      use domdec,only:rank
      use fields,only:forcefield
      use inform,only:deb_Path
      use kstbnd
      use sizes ,only:tinkerdebug
      use strbnd
      use tinMemory
      implicit none

 12   format(2x,'delete_data_kstrbnd')
      if(deb_Path) print 12

      call shmem_request(isb,  winisb,  [0,0], config=mhostacc)
      call shmem_request(sbk,  winsbk,  [0,0], config=mhostacc)
      if (forcefield .eq. 'MMFF94') return
      call shmem_request(nbstrbnd,winnbstrbnd,[0], config=mhostacc)
      end subroutine
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kstrbndm  --  assign MMFF str-bnd parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kstrbndm" assigns parameters for stretch-bend interactions
c     according to the Merck Molecular Force Field (MMFF)
c
c     Note: "stbnt" is the MMFF Stretch-Bend Type for angle "a-b-c",
c     where atom "a" has a smaller class number than atom "c"
c
c     if the BT of a-b = 1, then stbnt = 1
c     if the BT of b-c = 1, then stbnt = 2
c     if both = 1, then stbnt = 3
c     if 4-membered ring, then stbnt = 4
c     if 3-membered ring, then stbnt = 5
c     if 3-membered ring with BT of a-b = 1, then stbnt = 6
c     if 3-membered ring with BT of b-c = 1, then stbnt = 7
c     if 3-membered ring with BT of both = 1, then stbnt = 8
c     if 4-membered ring with BT of a-b = 1, then stbnt = 9
c     if 4-membered ring with BT of b-c = 1, then stbnt = 10
c     if 4-membered ring with BT of both = 1, then stbnt = 11
c     else, if all BT = 0 and no small ring, then stbnt = 0
c
c     literature references:
c
c     T. A. Halgren, "Merck Molecular Force Field. I. Basis, Form,
c     Scope, Parametrization, and Performance of MMFF94", Journal of
c     Computational Chemistry, 17, 490-519 (1995)
c
c     T. A. Halgren, "Merck Molecular Force Field. V. Extension of
c     MMFF94 Using Experimental Data, Additional Computational Data,
c     and Empirical Rules", Journal of Computational Chemistry, 17,
c     616-641 (1995)
c
c
      subroutine kstrbndm
      use sizes
      use angle
      use atmlst
      use atmtyp
      use couple
      use merck
      use potent
      use ring
      use strbnd
      use tinheader
      implicit none
      integer i,j,k,l,m
      integer ia,ib,ic
      integer ita,itb,itc
      integer ina,inb,inc
      integer ira,irb,irc
      integer nb1,nb2
      integer stbnt,ab,bc
      logical ring3,ring4
c
c
c     assign stretch-bend parameters for each angle
c
      nstrbnd = 0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
c
c     stretch-bend interactions are omitted for linear angles
c
         if (lin(class(ib)) .eq. 0) then
            ita = class (ia)
            itb = class (ib)
            itc = class (ic)
            ina = atomic(ia)
            inb = atomic(ib)
            inc = atomic(ic)
            sbk(1,nstrbnd+1) = 0.0_ti_p
            sbk(2,nstrbnd+1) = 0.0_ti_p
            do k = 1, n12(ib)
               if (i12(k,ib) .eq. ia)  nb1 = bndlist(k,ib)
               if (i12(k,ib) .eq. ic)  nb2 = bndlist(k,ib)
            end do
            stbnt = 0
            ab = 0
            bc = 0
c
c     check if the atoms belong to a single 3- or 4-membered ring
c
            ring3 = .false.
            ring4 = .false.
            do j = 1, nring3
               do k = 1, 3
                  if (ia .eq. iring3(k,j)) then
                     do l = 1, 3
                        if (ib .eq. iring3(l,j)) then
                           do m = 1, 3
                              if (ic .eq. iring3(m,j))
     &                           ring3 = .true.
                           end do
                        end if
                     end do
                  end if
               end do
            end do
            if (.not. ring3) then
               do j = 1, nring4
                  do k = 1, 4
                     if (ia .eq. iring4(k,j)) then
                        do l = 1, 4
                           if (ib .eq. iring4(l,j)) then
                              do m = 1, 4
                                 if (ic .eq. iring4(m,j))
     &                              ring4 = .true.
                              end do
                           end if
                        end do
                     end if
                  end do
               end do
            end if
c
c     determine the MMFF stretch-bend type for the current angle
c
            if (ita .lt. itc) then
               do j = 1, nlignes
                  if (((ia.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ia.eq.bt_1(j,2)))) then
                     ab = 1
                  end if
                  if (((ic.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ic.eq.bt_1(j,2)))) then
                     bc = 1
                  end if
               end do
               if (ab.eq.1 .and. bc.eq.0)  stbnt = 1
               if (ab.eq.0 .and. bc.eq.1)  stbnt = 2
               if (ab.eq.1 .and. bc.eq.1)  stbnt = 3
               if (stbnt.eq.0 .AND. ring3) then
                  stbnt = 5
               else if (stbnt.eq.1 .and. ring3) then
                  stbnt = 6
               else if (stbnt.eq.2 .and. ring3) then
                  stbnt = 7
               else if (stbnt.eq.3 .and. ring3) then
                  stbnt = 8
               else if (stbnt.eq.0 .and. ring4) then
                  stbnt = 4
               else if (stbnt.eq.1 .and. ring4) then
                  stbnt = 9
               else if (stbnt.eq.2 .and. ring4) then
                  stbnt = 10
               else if (stbnt.eq.3 .and. ring4) then
                  stbnt = 11
               end if
            else if (ita .gt. itc) then
               do j = 1, nlignes
                  if (((ia.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ia.eq.bt_1(j,2)))) then
                     ab = 1
                  end if
                  if (((ic.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ic.eq.bt_1(j,2)))) then
                     bc = 1
                  end if
               end do
               if (ab.eq.1 .and. bc.eq.0)  stbnt = 2
               if (ab.eq.0 .and. bc.eq.1)  stbnt = 1
               if (ab.eq.1 .and. bc.eq.1)  stbnt = 3
               if (stbnt.eq.0 .and. ring3) then
                  stbnt = 5
               else if (stbnt.eq.1 .and. ring3) then
                  stbnt = 6
               else if (stbnt.eq.2 .and. ring3) then
                  stbnt = 7
               else if (stbnt.eq.3 .and. ring3) then
                  stbnt = 8
               else if (stbnt.eq.0 .and. ring4) then
                  stbnt = 4
               else if (stbnt.eq.1 .and. ring4) then
                  stbnt = 9
               else if (stbnt.eq.2 .and. ring4) then
                  stbnt = 10
               else if (stbnt.eq.3 .and. ring4) then
                  stbnt = 11
               end if
            else if (ita .eq. itc) then
               do j = 1, nlignes
                  if (((ic.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ic.eq.bt_1(j,2)))) then
                     bc = 1
                  end if
                  if (((ia.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ia.eq.bt_1(j,2)))) then
                     ab = 1
                  end if
               end do
               if (ab.eq.1 .and. bc.eq.0)  stbnt = 1
               if (ab.eq.0 .and. bc.eq.1)  stbnt = 2
               if (ab.eq.1 .and. bc.eq.1)  stbnt = 3
               if (stbnt.eq.0 .and. ring3) then
                  stbnt = 5
               else if (stbnt.eq.1 .and. ring3) then
                  stbnt = 6
               else if (stbnt.eq.2 .and. ring3) then
                  stbnt = 7
               else if (stbnt.eq.3 .and. ring3) then
                  stbnt = 8
               else if (stbnt.eq.0 .and. ring4) then
                  stbnt = 4
               else if (stbnt.eq.1 .and. ring4) then
                  stbnt = 9
               else if (stbnt.eq.2 .and. ring4) then
                  stbnt = 10
               else if (stbnt.eq.3 .and. ring4) then
                  stbnt = 11
               end if
            end if
c
c     find the periodic table row for the atoms in the angle
c
            if (ina .eq. 1)  ira = 0
            if (ina.ge.3 .and. ina.le.10)  ira = 1
            if (ina.ge.11 .and. ina.le.18)  ira = 2
            if (ina.ge.19 .and. ina.le.36)  ira = 3
            if (ina.ge.37 .and. ina.le.54)  ira = 4
            if (inb .eq. 1)  irb = 0
            if (inb.ge.3 .and. inb.le.10)  irb = 1
            if (inb.ge.11 .and. inb.le.18)  irb = 2
            if (inb.ge.19 .and. inb.le.36)  irb = 3
            if (inb.ge.37 .and. inb.le.54)  irb = 4
            if (inc .eq. 1)  irc = 0
            if (inc.ge.3 .and. inc.le.10)  irc = 1
            if (inc.ge.11 .and. inc.le.18)  irc = 2
            if (inc.ge.19 .and. inc.le.36)  irc = 3
            if (inc.ge.37 .and. inc.le.54)  irc = 4
c
c     assign parameters via explicit values or empirical rules
c
            if (stbnt .eq. 11) then
               if ((stbn_abc11(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba11(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc11(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba11(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 10) then
               if ((stbn_abc10(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba10(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc10(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba10(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 9) then
               if ((stbn_abc9(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba9(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc9(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba9(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 8) then
               if ((stbn_abc8(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba3(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc8(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba8(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 7) then
               if ((stbn_abc7(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba7(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc7(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba7(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 6) then
               if ((stbn_abc6(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba3(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc6(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba6(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 5) then
               if (((stbn_abc5(ita,itb,itc).ne.1000.0_ti_p) .and.
     &              (stbn_cba3(ita,itb,itc).ne.1000.0_ti_p))
     &            .or. (ita.eq.22.and.itb.eq.22.and.itc.eq.22)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc5(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba5(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 4) then
               if ((stbn_abc4(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba4(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc4(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba4(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 3) then
               if ((stbn_abc3(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba3(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc3(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba3(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 2) then
              if ((stbn_abc2(ita,itb,itc).ne.1000.0_ti_p) .and.
     &            (stbn_cba2(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc2(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba2(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 1) then
               if ((stbn_abc1(ita,itb,itc).ne.1000.0_ti_p) .and.
     &             (stbn_cba1(ita,itb,itc).ne.1000.0_ti_p)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc1(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba1(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 0) then
               if (((stbn_abc(ita,itb,itc) .ne. 1000.0_ti_p) .and.
     &              (stbn_cba(ita,itb,itc) .ne. 1000.0_ti_p))
     &            .or. (ita.eq.12.AND.itb.eq.20.AND.itc.eq.20)
     &            .or. (ita.eq.20.AND.itb.eq.20.AND.itc.eq.12)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            end if
         end if
      end do
c
c     turn off the stretch-bend potential if it is not used
c
      if (nstrbnd .eq. 0) then
         use_strbnd = .false.
         return
      end if
      if (use_strbnd) then
         call upload_device_kstrbnd
      else
         call delete_data_kstrbnd
      end if
      end
