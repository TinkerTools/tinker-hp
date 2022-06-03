c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kangle  --  angle bend parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kangle" assigns the force constants and ideal angles for
c     the bond angles; also processes new or changed parameters
c
c
#include "tinker_precision.h"
      subroutine kangle
      use angle
      use angpot
      use atmtyp
      use couple
      use domdec
      use fields
      use kangs
      use keys
      use inform
      use iounit
      use potent
      use tinheader
      use usage
      implicit none
      integer i,j
      integer ia,ib,ic
      integer ita,itb,itc
      integer na,na5,na4
      integer na3,naf,nap,ierr
      integer jen,ih,nh
      integer next,size
      integer minat,iring
      real(t_p) fc,an,pr
      real(t_p) an1,an2,an3
      logical header,done
      logical use_ring
      character*4 pa,pb,pc
      character*6 label
      character*12 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing angle bending parameters
c
      if(deb_Path) print*,'kangle'
      blank = '         '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:6) .eq. 'ANGLE ')  iring = 0
         if (keyword(1:7) .eq. 'ANGLE5 ')  iring = 5
         if (keyword(1:7) .eq. 'ANGLE4 ')  iring = 4
         if (keyword(1:7) .eq. 'ANGLE3 ')  iring = 3
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an1 = 0.0_ti_p
            an2 = 0.0_ti_p
            an3 = 0.0_ti_p
            jen = 0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,fc,an1,an2,an3
   10       continue
            if (an2.ne.0.0_ti_p .or. an3.ne.0.0_ti_p)  jen = 1
            if (.not. silent) then
               if (header) then
                  header = .false.
                  if (rank.eq.0) write (iout,20)
   20             format (/,' Additional Angle Bending Parameters :',
     &                    //,5x,'Atom Classes',9x,'K(B)',7x,'Angle',/)
               end if
               if (iring .eq. 0) then
                  if (jen .eq. 0) then
                     if (rank.eq.0) write (iout,30)  ia,ib,ic,fc,an1
   30                format (4x,3i4,2x,2f12.3)
                  else if (an1 .ne. 0.0_ti_p) then
                     if (rank.eq.0) write (iout,40)  ia,ib,ic,fc,an1
   40                format (4x,3i4,2x,2f12.3,3x,'0-H''s')
                  end if
                  if (an2 .ne. 0.0_ti_p) then
                     if (rank.eq.0) write (iout,50)  ia,ib,ic,fc,an2
   50                format (4x,3i4,2x,2f12.3,3x,'1-H''s')
                  end if
                  if (an3 .ne. 0.0_ti_p) then
                     if (rank.eq.0) write (iout,60)  ia,ib,ic,fc,an3
   60                format (4x,3i4,2x,2f12.3,3x,'2-H''s')
                  end if
               else
                  if (iring .eq. 5)  label = '5-Ring'
                  if (iring .eq. 4)  label = '4-Ring'
                  if (iring .eq. 3)  label = '3-Ring'
                  if (jen .eq. 0) then
                     if (rank.eq.0) write (iout,70)  ia,ib,ic,fc,an1,
     $                label
   70                format (4x,3i4,2x,2f12.3,3x,a6)
                  else if (an1 .ne. 0.0_ti_p) then
                     if (rank.eq.0) write (iout,80)  ia,ib,ic,fc,an1,
     $                label
   80                format (4x,3i4,2x,2f12.3,3x,a6,3x,'0-H''s')
                  end if
                  if (an2 .ne. 0.0_ti_p) then
                     if (rank.eq.0) write (iout,90)  ia,ib,ic,fc,an2,
     $               label
   90                format (4x,3i4,2x,2f12.3,3x,a6,3x,'1-H''s')
                  end if
                  if (an3 .ne. 0.0_ti_p) then
                     if (rank.eq.0) write (iout,100)  ia,ib,ic,fc,an3,
     $               label
  100                format (4x,3i4,2x,2f12.3,3x,a6,3x,'2-H''s')
                  end if
               end if
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            if (iring .eq. 0) then
               do j = 1, maxna
                  if (ka(j).eq.blank .or. ka(j).eq.pt) then
                     ka(j) = pt
                     acon(j) = fc
                     ang(1,j) = an1
                     ang(2,j) = an2
                     ang(3,j) = an3
                     goto 150
                  end if
               end do
               if (rank.eq.0) write (iout,110)
  110          format (/,' KANGLE  --  Too many Bond Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 5) then
               do j = 1, maxna5
                  if (ka5(j).eq.blank .or. ka5(j).eq.pt) then
                     ka5(j) = pt
                     acon5(j) = fc
                     ang5(1,j) = an1
                     ang5(2,j) = an2
                     ang5(3,j) = an3
                     goto 150
                  end if
               end do
               if (rank.eq.0) write (iout,120)
  120          format (/,' KANGLE  --  Too many 5-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 4) then
               do j = 1, maxna4
                  if (ka4(j).eq.blank .or. ka4(j).eq.pt) then
                     ka4(j) = pt
                     acon4(j) = fc
                     ang4(1,j) = an1
                     ang4(2,j) = an2
                     ang4(3,j) = an3
                     goto 150
                  end if
               end do
               if (rank.eq.0) write (iout,130)
  130          format (/,' KANGLE  --  Too many 4-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 3) then
               do j = 1, maxna3
                  if (ka3(j).eq.blank .or. ka3(j).eq.pt) then
                     ka3(j) = pt
                     acon3(j) = fc
                     ang3(1,j) = an1
                     ang3(2,j) = an2
                     ang3(3,j) = an3
                     goto 150
                  end if
               end do
               if (rank.eq.0) write (iout,140)
  140          format (/,' KANGLE  --  Too many 3-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            end if
  150       continue
         end if
      end do
c
c     process keywords containing in-plane angle bending parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'ANGLEP ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an1 = 0.0_ti_p
            an2 = 0.0_ti_p
            jen = 0
            string = record(next:240)
            read (string,*,err=160,end=160)  ia,ib,ic,fc,an1,an2
  160       continue
            if (an2 .ne. 0.0_ti_p)  jen = 1
            if (.not. silent) then
               if (header) then
                  header = .false.
                  if (rank.eq.0) write (iout,170)
  170             format (/,' Additional In-Plane Angle Bending',
     &                       ' Parameters :',
     &                    //,5x,'Atom Classes',13x,'K(B)',10x,'Angle',/)
               end if
               if (jen .eq. 0) then
                  if (rank.eq.0) write (iout,180)  ia,ib,ic,fc,an1
  180             format (4x,3i4,3x,2f15.3)
               else if (an1 .ne. 0.0_ti_p) then
                  if (rank.eq.0) write (iout,190)  ia,ib,ic,fc,an1
  190             format (4x,3i4,3x,2f15.3,3x,'0-H''s')
               end if
               if (an2 .ne. 0.0_ti_p) then
                  if (rank.eq.0) write (iout,200)  ia,ib,ic,fc,an2
  200             format (4x,3i4,3x,2f15.3,3x,'1-H''s')
               end if
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, maxnap
               if (kap(j).eq.blank .or. kap(j).eq.pt) then
                  kap(j) = pt
                  aconp(j) = fc
                  angp(1,j) = an1
                  angp(2,j) = an2
                  goto 230
               end if
            end do
            write (iout,220)
  220       format (/,' KANGLE  --  Too many In-Plane Angle',
     &                    ' Bending Parameters')
            abort = .true.
  230       continue
         end if
      end do
c
c     process keywords containing Fourier angle bending parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:7) .eq. 'ANGLEF ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an = 0.0_ti_p
            pr = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=240,end=240)  ia,ib,ic,fc,an,pr
  240       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  if (rank.eq.0) write (iout,250)
  250             format (/,' Additional Fourier Angle Bending',
     &                       ' Parameters :',
     &                    //,5x,'Atom Classes',9x,'K(B)',7x,'Shift',
     &                       6x,'Period',/)
               end if
               if (rank.eq.0) write (iout,260)  ia,ib,ic,fc,an,pr
  260          format (4x,3i4,2x,3f12.3)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, maxnaf
               if (kaf(j).eq.blank .or. kaf(j).eq.pt) then
                  kaf(j) = pt
                  aconf(j) = fc
                  angf(1,j) = an
                  angf(2,j) = pr
                  goto 280
               end if
            end do
            if (rank.eq.0) write (iout,270)
  270       format (/,' KANGLE  --  Too many Fourier Angle',
     &                    ' Bending Parameters')
            abort = .true.
  280       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      na = maxna
      na5 = maxna5
      na4 = maxna4
      na3 = maxna3
      nap = maxnap
      naf = maxnaf
      do i = maxna, 1, -1
         if (ka(i) .eq. blank)  na = i - 1
      end do
      do i = maxna5, 1, -1
         if (ka5(i) .eq. blank)  na5 = i - 1
      end do
      do i = maxna4, 1, -1
         if (ka4(i) .eq. blank)  na4 = i - 1
      end do
      do i = maxna3, 1, -1
         if (ka3(i) .eq. blank)  na3 = i - 1
      end do
      do i = maxnap, 1, -1
         if (kap(i) .eq. blank)  nap = i - 1
      end do
      do i = maxnaf, 1, -1
         if (kaf(i) .eq. blank)  naf = i - 1
      end do
      use_ring = .false.
      if (min(na5,na4,na3) .ne. 0)  use_ring = .true.
c
c     set generic parameters for use with any number of hydrogens
c
      do i = 1, na
         if (ang(2,i).eq.0.0_ti_p .and. ang(3,i).eq.0.0_ti_p) then
            ang(2,i) = ang(1,i)
            ang(3,i) = ang(1,i)
         end if
      end do
      do i = 1, na5
         if (ang5(2,i).eq.0.0_ti_p .and. ang5(3,i).eq.0.0_ti_p) then
            ang5(2,i) = ang5(1,i)
            ang5(3,i) = ang5(1,i)
         end if
      end do
      do i = 1, na4
         if (ang4(2,i).eq.0.0_ti_p .and. ang4(3,i).eq.0.0_ti_p) then
            ang4(2,i) = ang4(1,i)
            ang4(3,i) = ang4(1,i)
         end if
      end do
      do i = 1, na3
         if (ang3(2,i).eq.0.0_ti_p .and. ang3(3,i).eq.0.0_ti_p) then
            ang3(2,i) = ang3(1,i)
            ang3(3,i) = ang3(1,i)
         end if
      end do
      do i = 1, nap
         if (angp(2,i) .eq. 0.0d0) then
            angp(2,i) = angp(1,i)
         end if
      end do
c
c     use special angle parameter assignment method for MMFF
c
      if (forcefield .eq. 'MMFF94') then
         call kanglem
         goto 400  !to the end
      end if
c
c     assign ideal bond angle and force constant for each angle
c
      header = .true.
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pt = pa//pb//pc
         else
            pt = pc//pb//pa
         end if
         ak(i) = 0.0_ti_p
         anat(i) = 0.0_ti_p
         afld(i) = 0.0_ti_p
         angtyp(i) = 'HARMONIC'
        angtypI(i) = ANG_HARMONIC
         done = .false.
c
c     count number of non-angle hydrogens on the central atom
c
         nh = 1
         do j = 1, n12(ib)
            ih = i12(j,ib)
            if (ih.ne.ia .and. ih.ne.ic .and. atomic(ih).eq.1)
     &         nh = nh + 1
         end do
c
c     make a check for bond angles contained inside small rings
c
         iring = 0
         if (use_ring) then
            call chkring (iring,ia,ib,ic,0)
            if (iring .eq. 6)  iring = 0
            if (iring.eq.5 .and. na5.eq.0)  iring = 0
            if (iring.eq.4 .and. na4.eq.0)  iring = 0
            if (iring.eq.3 .and. na3.eq.0)  iring = 0
         end if
c
c     assign angle bending parameters for bond angles
c
         if (iring .eq. 0) then
            do j = 1, na
               if (ka(j).eq.pt .and. ang(nh,j).ne.0.0_ti_p) then
                  ak(i) = acon(j)
                  anat(i) = ang(nh,j)
                  done = .true.
                  goto 290
               end if
            end do
c
c     assign bending parameters for 5-membered ring angles
c
         else if (iring .eq. 5) then
            do j = 1, na5
               if (ka5(j).eq.pt .and. ang5(nh,j).ne.0.0_ti_p) then
                  ak(i) = acon5(j)
                  anat(i) = ang5(nh,j)
                  done = .true.
                  goto 290
               end if
            end do
c
c     assign bending parameters for 4-membered ring angles
c
         else if (iring .eq. 4) then
            do j = 1, na4
               if (ka4(j).eq.pt .and. ang4(nh,j).ne.0.0_ti_p) then
                  ak(i) = acon4(j)
                  anat(i) = ang4(nh,j)
                  done = .true.
                  goto 290
               end if
            end do
c
c     assign bending parameters for 3-membered ring angles
c
         else if (iring .eq. 3) then
            do j = 1, na3
               if (ka3(j).eq.pt .and. ang3(nh,j).ne.0.0_ti_p) then
                  ak(i) = acon3(j)
                  anat(i) = ang3(nh,j)
                  done = .true.
                  goto 290
               end if
            end do
         end if
c
c     assign in-plane angle bending parameters for bond angles
c
         if (.not.done .and. n12(ib).eq.3) then
            do j = 1, nap
               if (kap(j).eq.pt .and. angp(nh,j).ne.0.0_ti_p) then
                  ak(i) = aconp(j)
                  anat(i) = angp(nh,j)
                  angtyp(i) = 'IN-PLANE'
                 angtypI(i) = ANG_IN_PLANE
                  done = .true.
                  goto 290
               end if
            end do
         end if
c
c     assign Fourier angle bending parameters for bond angles
c
         if (.not. done) then
            do j = 1, naf
               if (kaf(j) .eq. pt) then
                  ak(i) = aconf(j)
                  anat(i) = angf(1,j)
                  afld(i) = angf(2,j)
                  angtyp(i) = 'FOURIER'
                 angtypI(i) = ANG_FOURIER
                  done = .true.
                  goto 290
               end if
            end do
         end if
c
c     warning if suitable angle bending parameter not found
c
  290    continue
         minat = min(atomic(ia),atomic(ib),atomic(ic))
         if (minat .eq. 0)  done = .true.
         if (use_angle .and. .not.done) then
            if (use(ia) .or. use(ib) .or. use(ic))  abort = .true.
            if (header) then
               header = .false.
               if (rank.eq.0) write (iout,300)
  300          format (/,' Undefined Angle Bending Parameters :',
     &                 //,' Type',18x,'Atom Names',19x,
     &                    'Atom Classes',/)
            end if
            label = 'Angle '
            if (iring .eq. 5)  label = '5-Ring'
            if (iring .eq. 4)  label = '4-Ring'
            if (iring .eq. 3)  label = '3-Ring'
            write (iout,310)  label,ia,name(ia),ib,name(ib),
     &                        ic,name(ic),ita,itb,itc
  310       format (1x,a6,5x,3(i6,'-',a3),7x,3i5)
         end if
      end do
c
c     turn off the angle bending potential if it is not used
c
      if (nangle .eq. 0)  use_angle = .false.
  400 continue
      call MPI_BARRIER(hostcomm,ierr)
      call upload_device_kangle
      end

      subroutine upload_device_kangle
      use angle
      use angpot
      use domdec,only:rank
      use sizes,only:tinkerdebug
      use tinMemory
      implicit none
#ifdef _OPENACC
 12   format(2x,'upload_device_kangle')
      if(rank.eq.0.and.tinkerdebug) print 12
#endif
!$acc update device(afld,anat,ak,angtypI)
      end subroutine

      subroutine delete_data_kangle
      use angle
      use angpot
      use domdec,only:rank
      use potent
      use sizes,only:tinkerdebug
      use tinMemory
      implicit none
 12   format(2x,'delete_data_kangle')
      if(rank.eq.0.and.tinkerdebug) print 12

      call shmem_request(ak,     winak,     [0], config=mhostacc)
      call shmem_request(anat,   winanat,   [0], config=mhostacc)
      call shmem_request(afld,   winafld,   [0], config=mhostacc)
      call shmem_request(angtyp, winangtyp, [0], config=mhostonly)
      call shmem_request(angtypI,winangtypI,[0], config=mhostacc)
      end subroutine
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kanglem  --  MMFF angle parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kanglem" assigns the force constants and ideal angles for
c     bond angles according to the Merck Molecular Force Field (MMFF)
c
c     literature reference:
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
      subroutine kanglem
      use angle
      use angpot
      use atmtyp
      use atoms
      use bond
      use merck
      use potent
      use ring
      use tinheader
      implicit none
      integer i,j,k,l,m
      integer ia,ib,ic
      integer ita,itb,itc
      integer ina,inb,inc
      integer itta,ittb,ittc
      integer bnd_ab,bnd_bc
      integer at,minat
      integer mclass
      real(t_p) d,beta
      real(t_p) z2(100),c(100)
      logical done
      logical ring3,ring4
c
c
c     set empirical rule parameters for some common elements
c
      do i = 1, 100
         z2(i) = 1000.0_ti_p
         c(i) = 1000.0_ti_p
      end do
      z2(1) = 1.395_ti_p
      z2(5) = 0.0_ti_p
      z2(6) = 2.494_ti_p
      z2(7) = 2.711_ti_p
      z2(8) = 3.045_ti_p
      z2(9) = 2.847_ti_p
      z2(14) = 2.350_ti_p
      z2(15) = 2.350_ti_p
      z2(16) = 2.980_ti_p
      z2(17) = 2.909_ti_p
      z2(35) = 3.017_ti_p
      z2(33) = 0.0_ti_p
      z2(53) = 3.086_ti_p
      c(1) = 0.0_ti_p
      c(5) = 0.704_ti_p
      c(6) = 1.016_ti_p
      c(7) = 1.113_ti_p
      c(8) = 1.337_ti_p
      c(9) = 0.0_ti_p
      c(14) = 0.811_ti_p
      c(15) = 1.068_ti_p
      c(16) = 1.249_ti_p
      c(17) = 1.078_ti_p
      c(35) = 0.0_ti_p
      c(33) = 0.825_ti_p
      c(53) = 0.0_ti_p
c
c     assign MMFF bond angle and force constant for each angle
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itta = type(ia)
         ittb = type(ib)
         ittc = type(ic)
         ina = atomic(ia)
         inb = atomic(ib)
         inc = atomic(ic)
c
c     set angle index value, accounting for MMFF bond type = 1
c
         at = 0
         do j = 1, nlignes
            if ((ia.eq.bt_1(j,1) .and. ib.eq.bt_1(j,2)) .or.
     &          (ib.eq.bt_1(j,1) .and. ia.eq.bt_1(j,2))) then
               at = at + 1
            end if
            if ((ic.eq.bt_1(j,1) .and. ib.eq.bt_1(j,2)) .or.
     &          (ib.eq.bt_1(j,1) .and. ic.eq.bt_1(j,2))) then
               at = at + 1
            end if
         end do
c
c     determine if the atoms belong to a 3- or 4-membered ring
c
         ring3 = .false.
         ring4 = .false.
         do j = 1, nring3
            do k = 1, 3
               if (ia .eq. iring3(k,j)) then
                  do l = 1, 3
                     if (ib .eq. iring3(l,j)) then
                        do m = 1, 3
                           if (ic .eq. iring3(m,j))  ring3 = .true.
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
                              if (ic .eq. iring4(m,j))  ring4 = .true.
                           end do
                        end if
                     end do
                  end if
               end do
            end do
         end if
c
c     set special index value when 3- or 4-rings are present
c
         if (at.eq.0 .and. ring4) then
            at = 4
         else if (at.eq.1 .and. ring4) then
            at = 7
         else if (at.eq.2 .and. ring4) then
            at = 8
         else if (at.eq.0 .and. ring3) then
            at = 3
         else if (at.eq.1 .and. ring3) then
            at = 5
         else if (at.eq.2 .and. ring3) then
            at = 6
         end if
c
c     setup the atom class equivalencies assignment
c
         mclass = 0
   10    continue
         mclass = mclass + 1
         if (mclass .eq. 1) then
            ita = eqclass(itta,1)
            itb = eqclass(ittb,1)
            itc = eqclass(ittc,1)
         else if (mclass .eq. 2) then
            ita = eqclass(itta,2)
            itb = eqclass(ittb,2)
            itc = eqclass(ittc,2)
         else if (mclass .eq. 3) then
            ita = eqclass(itta,3)
            itb = eqclass(ittb,2)
            itc = eqclass(ittc,3)
         else if (mclass .eq. 4) then
            ita = eqclass(itta,4)
            itb = eqclass(ittb,2)
            itc = eqclass(ittc,4)
         else if (mclass .eq. 5) then
            ita = eqclass(itta,5)
            itb = eqclass(ittb,2)
            itc = eqclass(ittc,5)
         end if
         if (mclass .gt. 5) then
            goto 20
         else
            if (at .eq. 0) then
               ak(i) = mmff_ka(ita,itb,itc)
               anat(i) = mmff_ang0(ita,itb,itc)
            else if (at .eq. 1) then
               ak(i) = mmff_ka1(ita,itb,itc)
               anat(i) = mmff_ang1(ita,itb,itc)
            else if (at .eq. 2) then
               ak(i) = mmff_ka2(ita,itb,itc)
               anat(i) = mmff_ang2(ita,itb,itc)
            else if (at .eq. 3) then
               ak(i) = mmff_ka3(ita,itb,itc)
               anat(i) = mmff_ang3(ita,itb,itc)
            else if (at .eq. 4) then
               ak(i) = mmff_ka4(ita,itb,itc)
               anat(i) = mmff_ang4(ita,itb,itc)
            else if (at .eq. 5) then
               ak(i) = mmff_ka5(ita,itb,itc)
               anat(i) = mmff_ang5(ita,itb,itc)
            else if (at .eq. 6) then
               ak(i) = mmff_ka6(ita,itb,itc)
               anat(i) = mmff_ang6(ita,itb,itc)
            else if (at .eq. 7) then
               ak(i) = mmff_ka7(ita,itb,itc)
               anat(i) = mmff_ang7(ita,itb,itc)
            else if (at .eq. 8) then
               ak(i) = mmff_ka8(ita,itb,itc)
               anat(i) = mmff_ang8(ita,itb,itc)
            end if
c
c     use empirical rule to calculate the force constant
c
            if (mclass .eq. 5) then
               if (z2(ina) .eq. 1000.0_ti_p)  goto 20
               if (z2(inb) .eq. 1000.0_ti_p)  goto 20
               if (z2(inc) .eq. 1000.0_ti_p)  goto 20
               if (c(ina) .eq. 1000.0_ti_p)  goto 20
               if (c(inb) .eq. 1000.0_ti_p)  goto 20
               if (c(inc) .eq. 1000.0_ti_p)  goto 20
               do k = 1, nbond
                  if ((min(ia,ib).eq.ibnd(1,k)) .and.
     &                (max(ia,ib).eq.ibnd(2,k))) then
                     bnd_ab = k
                  end if
                  if ((min(ic,ib).eq.ibnd(1,k)) .and.
     &                (max(ic,ib).eq.ibnd(2,k))) then
                     bnd_bc = k
                  end if
               end do
               d = (bl(bnd_ab)-bl(bnd_bc))**2
     &                / (bl(bnd_ab)+bl(bnd_bc))**2
               beta = 1.0_ti_p
               if (ring4)  beta = 0.85_ti_p
               if (ring3)  beta = 0.05_ti_p
               ak(i) = beta*1.75_ti_p*z2(ina)*z2(inc)*c(inb)
     &                 / ((0.01745329252_ti_p*anat(i))**2
     &                      *(bl(bnd_ab)+bl(bnd_bc))*exp(2.0_ti_p*d))
            end if
            done = .true.
            if (ak(i) .eq. 1000.0_ti_p)  done = .false.
            if (anat(i) .eq. 1000.0_ti_p)  done = .false.
            if (.not. done)  goto 10
            goto 20
         end if
c
c     use empirical rule for ideal angle and force constant
c
   20    continue
         minat = min(ina,inb,inc)
         if (minat .eq. 0)  done = .true.
         if (.not. done) then
            if (use_angle) then
               anat(i) = 120.0_ti_p
               if (crd(itb) .eq. 4)  anat(i) = 109.45_ti_p
               if (crd(itb) .eq. 2) then
                  if (inb .eq. 8) then
                     anat(i) = 105.0_ti_p
                  else if (inb .gt. 10) then
                     anat(i) = 95.0_ti_p
                  else if (lin(itb) .eq. 1) then
                     anat(i) = 180.0_ti_p
                  end if
               end if
               if (crd(itb).eq.3 .and. val(itb).eq.3
     &                .and. mltb(itb).eq.0) then
                  if (inb .eq. 7) then
                     anat(i) = 107.0_ti_p
                  else
                     anat(i) = 92.0_ti_p
                  end if
               end if
               if (ring3)  anat(i) = 60.0_ti_p
               if (ring4)  anat(i) = 90.0_ti_p
               do k = 1, nbond
                  if ((min(ia,ib).eq.ibnd(1,k)) .and.
     &                (max(ia,ib).eq.ibnd(2,k))) then
                     bnd_ab = k
                  end if
                  if ((min(ic,ib).eq.ibnd(1,k)) .and.
     &                (max(ic,ib).eq.ibnd(2,k))) then
                     bnd_bc = k
                  end if
               end do
               d = (bl(bnd_ab)-bl(bnd_bc))**2
     &                / (bl(bnd_ab)+bl(bnd_bc))**2
               beta = 1.0_ti_p
               if (ring4)  beta = 0.85_ti_p
               if (ring3)  beta = 0.05_ti_p
               ak(i) = beta*1.75_ti_p*z2(ina)*z2(inc)*c(inb)
     &                 / ((0.01745329252_ti_p*anat(i))**2
     &                      *(bl(bnd_ab)+bl(bnd_bc))*exp(2.0_ti_p*d))
            end if
         end if
         angtyp(i) = 'HARMONIC'
        angtypI(i) = ANG_HARMONIC
         if (anat(i) .eq. 180.0_ti_p) then
            angtyp(i) = 'LINEAR'
           angtypI(i) = ANG_LINEAR
         end if
      end do
c
c     turn off the angle bending potential if it is not used
c
      if (nangle .eq. 0)  use_angle = .false.
!$acc update device(beta)
      end
