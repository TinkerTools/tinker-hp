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
      use usage
      implicit none
      integer i,j
      integer ia,ib,ic
      integer ita,itb,itc
      integer na,na5,na4
      integer na3,naf,nap,naps
      integer jen,ih,nh
      integer next,size
      integer minat,iring
      real*8 fc,an,pr
      real*8 an1,an2,an3
      logical header,done
      logical use_ring
      character*4 pa,pb,pc
      character*6 label
      character*12 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
      logical :: max_reach
c
c
c     process keywords containing angle bending parameters
c
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
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            jen = 0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,fc,an1,an2,an3
   10       continue
            if (an2.ne.0.0d0 .or. an3.ne.0.0d0)  jen = 1
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
                  else if (an1 .ne. 0.0d0) then
                     if (rank.eq.0) write (iout,40)  ia,ib,ic,fc,an1
   40                format (4x,3i4,2x,2f12.3,3x,'0-H''s')
                  end if
                  if (an2 .ne. 0.0d0) then
                     if (rank.eq.0) write (iout,50)  ia,ib,ic,fc,an2
   50                format (4x,3i4,2x,2f12.3,3x,'1-H''s')
                  end if
                  if (an3 .ne. 0.0d0) then
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
                  else if (an1 .ne. 0.0d0) then
                     if (rank.eq.0) write (iout,80)  ia,ib,ic,fc,an1,
     $                label
   80                format (4x,3i4,2x,2f12.3,3x,a6,3x,'0-H''s')
                  end if
                  if (an2 .ne. 0.0d0) then
                     if (rank.eq.0) write (iout,90)  ia,ib,ic,fc,an2,
     $               label
   90                format (4x,3i4,2x,2f12.3,3x,a6,3x,'1-H''s')
                  end if
                  if (an3 .ne. 0.0d0) then
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
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            jen = 0
            string = record(next:240)
            read (string,*,err=160,end=160)  ia,ib,ic,fc,an1,an2
  160       continue
            if (an2 .ne. 0.0d0)  jen = 1
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
               else if (an1 .ne. 0.0d0) then
                  if (rank.eq.0) write (iout,190)  ia,ib,ic,fc,an1
  190             format (4x,3i4,3x,2f15.3,3x,'0-H''s')
               end if
               if (an2 .ne. 0.0d0) then
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
            fc = 0.0d0
            an = 0.0d0
            pr = 0.0d0
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
c     process keywords containing Partridge-Schwenke angle bending parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:8) .eq. 'ANGLEPS ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an = 0.0d0
            pr = 0.0d0
            string = record(next:240)
            read (string,*)  ia,ib,ic,fc,an,pr
            if (.not. silent) then
               if (header) then
                  header = .false.
                  if (rank.eq.0) write (iout,'(A,5x,A,9x,A,7x,A,6x,A)')
     &                  ' Additional Partridge-Scwhenke Angle ',
     &                       ' Bending Parameters :',
     &                    'Atom Classes','theta_e','r_e','beta'
               end if
               if (rank.eq.0) write (iout,'(4x,3i4,2x,3f12.3)')  
     &           ia,ib,ic,fc,an,pr
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
            max_reach = .true.
            do j = 1, maxnaps
               if (kaps(j).eq.blank .or. kaps(j).eq.pt) then
                  kaps(j) = pt
                  angps(1,j) = fc
                  angps(2,j) = an
                  angps(3,j) = pr
                  max_reach=.false.
                  exit
               end if
            end do
            if(max_reach) then
              if (rank.eq.0) write (iout,*)
     &          ' KANGLE  --  Too many P-S Angle',
     &                    ' Bending Parameters'
              abort = .true.
            endif
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
      naps = maxnaps
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
      do i = maxnaps, 1, -1
         if (kaps(i) .eq. blank)  naps = i - 1
      end do
      use_ring = .false.
      if (min(na5,na4,na3) .ne. 0)  use_ring = .true.

      if (naps > 0) call initialize_angleps()
c
c     set generic parameters for use with any number of hydrogens
c
      do i = 1, na
         if (ang(2,i).eq.0.0d0 .and. ang(3,i).eq.0.0d0) then
            ang(2,i) = ang(1,i)
            ang(3,i) = ang(1,i)
         end if
      end do
      do i = 1, na5
         if (ang5(2,i).eq.0.0d0 .and. ang5(3,i).eq.0.0d0) then
            ang5(2,i) = ang5(1,i)
            ang5(3,i) = ang5(1,i)
         end if
      end do
      do i = 1, na4
         if (ang4(2,i).eq.0.0d0 .and. ang4(3,i).eq.0.0d0) then
            ang4(2,i) = ang4(1,i)
            ang4(3,i) = ang4(1,i)
         end if
      end do
      do i = 1, na3
         if (ang3(2,i).eq.0.0d0 .and. ang3(3,i).eq.0.0d0) then
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
         return
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
         ak(i) = 0.0d0
         anat(i) = 0.0d0
         afld(i) = 0.0d0
         angtyp(i) = 'HARMONIC'
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
               if (ka(j).eq.pt .and. ang(nh,j).ne.0.0d0) then
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
               if (ka5(j).eq.pt .and. ang5(nh,j).ne.0.0d0) then
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
               if (ka4(j).eq.pt .and. ang4(nh,j).ne.0.0d0) then
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
               if (ka3(j).eq.pt .and. ang3(nh,j).ne.0.0d0) then
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
               if (kap(j).eq.pt .and. angp(nh,j).ne.0.0d0) then
                  ak(i) = aconp(j)
                  anat(i) = angp(nh,j)
                  angtyp(i) = 'IN-PLANE'
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
                  done = .true.
                  goto 290
               end if
            end do
         end if
c
c     assign Partridge-Schwenke angle bending parameters
c
         if (.not. done) then
            do j = 1, naps
               if (kaps(j) .eq. pt) then
                  anat(i) = angps(1,j) !theta_e
                  afld(i) = angps(2,j) !r_e
                    ak(i) = angps(3,j) !beta
                  angtyp(i) = 'ANGLEPS'
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
      return
      end
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
      use domdec
      use merck
      use potent
      use ring
      use mpi
      implicit none
      integer i,j,k,l,m
      integer ia,ib,ic
      integer ita,itb,itc
      integer ina,inb,inc
      integer itta,ittb,ittc
      integer bnd_ab,bnd_bc
      integer at,minat
      integer mclass
      integer ierr
      real*8 d,beta
      real*8 z2(100),c(100)
      logical done
      logical ring3,ring4
c
c
c     set empirical rule parameters for some common elements
c
      if (hostrank.ne.0) goto 30
      do i = 1, 100
         z2(i) = 1000.0d0
         c(i) = 1000.0d0
      end do
      z2(1) = 1.395d0
      z2(5) = 0.0d0
      z2(6) = 2.494d0
      z2(7) = 2.711d0
      z2(8) = 3.045d0
      z2(9) = 2.847d0
      z2(14) = 2.350d0
      z2(15) = 2.350d0
      z2(16) = 2.980d0
      z2(17) = 2.909d0
      z2(35) = 3.017d0
      z2(33) = 0.0d0
      z2(53) = 3.086d0
      c(1) = 0.0d0
      c(5) = 0.704d0
      c(6) = 1.016d0
      c(7) = 1.113d0
      c(8) = 1.337d0
      c(9) = 0.0d0
      c(14) = 0.811d0
      c(15) = 1.068d0
      c(16) = 1.249d0
      c(17) = 1.078d0
      c(35) = 0.0d0
      c(33) = 0.825d0
      c(53) = 0.0d0
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
               if (z2(ina) .eq. 1000.0d0)  goto 20
               if (z2(inb) .eq. 1000.0d0)  goto 20
               if (z2(inc) .eq. 1000.0d0)  goto 20
               if (c(ina) .eq. 1000.0d0)  goto 20
               if (c(inb) .eq. 1000.0d0)  goto 20
               if (c(inc) .eq. 1000.0d0)  goto 20
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
               beta = 1.0d0
               if (ring4)  beta = 0.85d0
               if (ring3)  beta = 0.05d0
               ak(i) = beta*1.75d0*z2(ina)*z2(inc)*c(inb)
     &                 / ((0.01745329252d0*anat(i))**2
     &                      *(bl(bnd_ab)+bl(bnd_bc))*exp(2.0d0*d))
            end if
            done = .true.
            if (ak(i) .eq. 1000.0d0)  done = .false.
            if (anat(i) .eq. 1000.0d0)  done = .false.
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
               anat(i) = 120.0d0
               if (crd(itb) .eq. 4)  anat(i) = 109.45d0
               if (crd(itb) .eq. 2) then
                  if (inb .eq. 8) then
                     anat(i) = 105.0d0
                  else if (inb .gt. 10) then
                     anat(i) = 95.0d0
                  else if (lin(itb) .eq. 1) then
                     anat(i) = 180.0d0
                  end if
               end if
               if (crd(itb).eq.3 .and. val(itb).eq.3
     &                .and. mltb(itb).eq.0) then
                  if (inb .eq. 7) then
                     anat(i) = 107.0d0
                  else
                     anat(i) = 92.0d0
                  end if
               end if
               if (ring3)  anat(i) = 60.0d0
               if (ring4)  anat(i) = 90.0d0
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
               beta = 1.0d0
               if (ring4)  beta = 0.85d0
               if (ring3)  beta = 0.05d0
               ak(i) = beta*1.75d0*z2(ina)*z2(inc)*c(inb)
     &                 / ((0.01745329252d0*anat(i))**2
     &                      *(bl(bnd_ab)+bl(bnd_bc))*exp(2.0d0*d))
            end if
         end if
         angtyp(i) = 'HARMONIC'
         if (anat(i) .eq. 180.0d0)  angtyp(i) = 'LINEAR'
      end do
 30   call MPI_BARRIER(hostcomm,ierr)
c
c     turn off the angle bending potential if it is not used
c
      if (nangle .eq. 0)  use_angle = .false.
      return
      end

c
c     Initialize parameters for Partridge-Scwhenke angle terms
c
      subroutine initialize_angleps()
      use angpot
      use units
      implicit none
      integer :: i
      real*8 :: cbasis(245),ccore(245),crest(245)
      real*8, parameter :: f5z = 0.99967788500000d0
     &                    ,fbasis = 0.15860145369897d0
     &                    ,fcore = -1.6351695982132d0
     &                    ,frest = 1d0

      idx_ps(:,1) = (/ 
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 
     &   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
     &   2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &   3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 
     &   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
     &   4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 
     &   4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 
     &   6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
     &   6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 
     &   5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 
     &   7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 
     &   6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 
     &   9, 9, 9, 9, 9 /)
      idx_ps(:,2) = (/ 
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
     &   2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 
     &   2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 
     &   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
     &   2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 
     &   3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 
     &   1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &   2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 
     &   4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 
     &   2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 
     &   4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 
     &   1, 1, 1, 1, 1 /)
      idx_ps(:,3) = (/ 
     &   1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5, 
     &   6, 7, 8, 9,10,11,12,13,14, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 
     &   12,13, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 2, 3, 4, 5,
     &   6, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 1, 
     &   2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 
     &   11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8,
     &   9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8, 
     &   9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
     &   1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 
     &   3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 
     &   7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 
     &   4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 
     &   3, 4, 5, 6, 7/)
      c5z_ps(:)= (/ 
     &  4.2278462684916D+04, 4.5859382909906D-02, 9.4804986183058D+03,
     &  7.5485566680955D+02, 1.9865052511496D+03, 4.3768071560862D+02,
     &  1.4466054104131D+03, 1.3591924557890D+02,-1.4299027252645D+03,
     &  6.6966329416373D+02, 3.8065088734195D+03,-5.0582552618154D+02,
     &  -3.2067534385604D+03, 6.9673382568135D+02, 1.6789085874578D+03,
     &  -3.5387509130093D+03,-1.2902326455736D+04,-6.4271125232353D+03,
     &  -6.9346876863641D+03,-4.9765266152649D+02,-3.4380943579627D+03,
     &  3.9925274973255D+03,-1.2703668547457D+04,-1.5831591056092D+04,
     &  2.9431777405339D+04, 2.5071411925779D+04,-4.8518811956397D+04,
     &  -1.4430705306580D+04, 2.5844109323395D+04,-2.3371683301770D+03,
     &  1.2333872678202D+04, 6.6525207018832D+03,-2.0884209672231D+03,
     &  -6.3008463062877D+03, 4.2548148298119D+04, 2.1561445953347D+04,
     &  -1.5517277060400D+05, 2.9277086555691D+04, 2.6154026873478D+05,
     &  -1.3093666159230D+05,-1.6260425387088D+05, 1.2311652217133D+05,
     &  -5.1764697159603D+04, 2.5287599662992D+03, 3.0114701659513D+04,
     &  -2.0580084492150D+03, 3.3617940269402D+04, 1.3503379582016D+04,
     &  -1.0401149481887D+05,-6.3248258344140D+04, 2.4576697811922D+05,
     &  8.9685253338525D+04,-2.3910076031416D+05,-6.5265145723160D+04,
     &  8.9184290973880D+04,-8.0850272976101D+03,-3.1054961140464D+04,
     &  -1.3684354599285D+04, 9.3754012976495D+03,-7.4676475789329D+04,
     &  -1.8122270942076D+05, 2.6987309391410D+05, 4.0582251904706D+05,
     &  -4.7103517814752D+05,-3.6115503974010D+05, 3.2284775325099D+05,
     &  1.3264691929787D+04, 1.8025253924335D+05,-1.2235925565102D+04,
     &  -9.1363898120735D+03,-4.1294242946858D+04,-3.4995730900098D+04,
     &  3.1769893347165D+05, 2.8395605362570D+05,-1.0784536354219D+06,
     &  -5.9451106980882D+05, 1.5215430060937D+06, 4.5943167339298D+05,
     &  -7.9957883936866D+05,-9.2432840622294D+04, 5.5825423140341D+03,
     &  3.0673594098716D+03, 8.7439532014842D+04, 1.9113438435651D+05,
     &  -3.4306742659939D+05,-3.0711488132651D+05, 6.2118702580693D+05,
     &  -1.5805976377422D+04,-4.2038045404190D+05, 3.4847108834282D+05,
     &  -1.3486811106770D+04, 3.1256632170871D+04, 5.3344700235019D+03,
     &  2.6384242145376D+04, 1.2917121516510D+05,-1.3160848301195D+05,
     &  -4.5853998051192D+05, 3.5760105069089D+05, 6.4570143281747D+05,
     &  -3.6980075904167D+05,-3.2941029518332D+05,-3.5042507366553D+05,
     &  2.1513919629391D+03, 6.3403845616538D+04, 6.2152822008047D+04,
     &  -4.8805335375295D+05,-6.3261951398766D+05, 1.8433340786742D+06,
     &  1.4650263449690D+06,-2.9204939728308D+06,-1.1011338105757D+06,
     &  1.7270664922758D+06, 3.4925947462024D+05,-1.9526251371308D+04,
     &  -3.2271030511683D+04,-3.7601575719875D+05, 1.8295007005531D+05,
     &  1.5005699079799D+06,-1.2350076538617D+06,-1.8221938812193D+06,
     &  1.5438780841786D+06,-3.2729150692367D+03, 1.0546285883943D+04,
     &  -4.7118461673723D+04,-1.1458551385925D+05, 2.7704588008958D+05,
     &  7.4145816862032D+05,-6.6864945408289D+05,-1.6992324545166D+06,
     &  6.7487333473248D+05, 1.4361670430046D+06,-2.0837555267331D+05,
     &  4.7678355561019D+05,-1.5194821786066D+04,-1.1987249931134D+05,
     &  1.3007675671713D+05, 9.6641544907323D+05,-5.3379849922258D+05,
     &  -2.4303858824867D+06, 1.5261649025605D+06, 2.0186755858342D+06,
     &  -1.6429544469130D+06,-1.7921520714752D+04, 1.4125624734639D+04,
     &  -2.5345006031695D+04, 1.7853375909076D+05,-5.4318156343922D+04,
     &  -3.6889685715963D+05, 4.2449670705837D+05, 3.5020329799394D+05,
     &  9.3825886484788D+03,-8.0012127425648D+05, 9.8554789856472D+04,
     &  4.9210554266522D+05,-6.4038493953446D+05,-2.8398085766046D+06,
     &  2.1390360019254D+06, 6.3452935017176D+06,-2.3677386290925D+06,
     &  -3.9697874352050D+06,-1.9490691547041D+04, 4.4213579019433D+04,
     &  1.6113884156437D+05,-7.1247665213713D+05,-1.1808376404616D+06,
     &  3.0815171952564D+06, 1.3519809705593D+06,-3.4457898745450D+06,
     &  2.0705775494050D+05,-4.3778169926622D+05, 8.7041260169714D+03,
     &  1.8982512628535D+05,-2.9708215504578D+05,-8.8213012222074D+05,
     &  8.6031109049755D+05, 1.0968800857081D+06,-1.0114716732602D+06,
     &  1.9367263614108D+05, 2.8678295007137D+05,-9.4347729862989D+04,
     &  4.4154039394108D+04, 5.3686756196439D+05, 1.7254041770855D+05,
     &  -2.5310674462399D+06,-2.0381171865455D+06, 3.3780796258176D+06,
     &  7.8836220768478D+05,-1.5307728782887D+05,-3.7573362053757D+05,
     &  1.0124501604626D+06, 2.0929686545723D+06,-5.7305706586465D+06,
     &  -2.6200352535413D+06, 7.1543745536691D+06,-1.9733601879064D+04,
     &  8.5273008477607D+04, 6.1062454495045D+04,-2.2642508675984D+05,
     &  2.4581653864150D+05,-9.0376851105383D+05,-4.4367930945690D+05,
     &  1.5740351463593D+06, 2.4563041445249D+05,-3.4697646046367D+03,
     &  -2.1391370322552D+05, 4.2358948404842D+05, 5.6270081955003D+05,
     &  -8.5007851251980D+05,-6.1182429537130D+05, 5.6690751824341D+05,
     &  -3.5617502919487D+05,-8.1875263381402D+02,-2.4506258140060D+05,
     &  2.5830513731509D+05, 6.0646114465433D+05,-6.9676584616955D+05,
     &  5.1937406389690D+05, 1.7261913546007D+05,-1.7405787307472D+04,
     &  -3.8301842660567D+05, 5.4227693205154D+05, 2.5442083515211D+06,
     &  -1.1837755702370D+06,-1.9381959088092D+06,-4.0642141553575D+05,
     &  1.1840693827934D+04,-1.5334500255967D+05, 4.9098619510989D+05,
     &  6.1688992640977D+05, 2.2351144690009D+05,-1.8550462739570D+06,
     &  9.6815110649918D+03,-8.1526584681055D+04,-8.0810433155289D+04,
     &  3.4520506615177D+05, 2.5509863381419D+05,-1.3331224992157D+05,
     &  -4.3119301071653D+05,-5.9818343115856D+04, 1.7863692414573D+03,
     &  8.9440694919836D+04,-2.5558967650731D+05,-2.2130423988459D+04,
     &  4.4973674518316D+05,-2.2094939343618D+05/)
      cbasis(:)= (/
     &  6.9770019624764D-04,-2.4209870001642D+01, 1.8113927151562D+01,
     &  3.5107416275981D+01,-5.4600021126735D+00,-4.8731149608386D+01,
     &  3.6007189184766D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  -7.7178474355102D+01,-3.8460795013977D+01,-4.6622480912340D+01,
     &  5.5684951167513D+01, 1.2274939911242D+02,-1.4325154752086D+02,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00,-6.0800589055949D+00,
     &  8.6171499453475D+01,-8.4066835441327D+01,-5.8228085624620D+01,
     &  2.0237393793875D+02, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  3.3525582670313D+02, 7.0056962392208D+01,-4.5312502936708D+01,
     &  -3.0441141194247D+02, 2.8111438108965D+02, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-1.2983583774779D+02, 3.9781671212935D+01,
     &  -6.6793945229609D+01,-1.9259805675433D+02, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-8.2855757669957D+02,-5.7003072730941D+01,
     &  -3.5604806670066D+01, 9.6277766002709D+01, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 8.8645622149112D+02,-7.6908409772041D+01,
     &  6.8111763314154D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  2.5090493428062D+02,-2.3622141780572D+02, 5.8155647658455D+02,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 2.8919570295095D+03,
     &  -1.7871014635921D+02,-1.3515667622500D+02, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-3.6965613754734D+03, 2.1148158286617D+02,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00,-1.4795670139431D+03,
     &  3.6210798138768D+02, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  -5.3552886800881D+03, 3.1006384016202D+02, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 1.6241824368764D+03, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 4.3764909606382D+03, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 1.0940849243716D+03, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 3.0743267832931D+03, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00/)
      ccore(:)= (/
     &  2.4332191647159D-02,-2.9749090113656D+01, 1.8638980892831D+01,
     &  -6.1272361746520D+00, 2.1567487597605D+00,-1.5552044084945D+01,
     &  8.9752150543954D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  -3.5693557878741D+02,-3.0398393196894D+00,-6.5936553294576D+00,
     &  1.6056619388911D+01, 7.8061422868204D+01,-8.6270891686359D+01,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00,-3.1688002530217D+01,
     &  3.7586725583944D+01,-3.2725765966657D+01,-5.6458213299259D+00,
     &  2.1502613314595D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  5.2789943583277D+02,-4.2461079404962D+00,-2.4937638543122D+01,
     &  -1.1963809321312D+02, 2.0240663228078D+02, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-6.2574211352272D+02,-6.9617539465382D+00,
     &  -5.9440243471241D+01, 1.4944220180218D+01, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-1.2851139918332D+03,-6.5043516710835D+00,
     &  4.0410829440249D+01,-6.7162452402027D+01, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 1.0031942127832D+03, 7.6137226541944D+01,
     &  -2.7279242226902D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  -3.3059000871075D+01, 2.4384498749480D+01,-1.4597931874215D+02,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 1.6559579606045D+03,
     &  1.5038996611400D+02,-7.3865347730818D+01, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-1.9738401290808D+03,-1.4149993809415D+02,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00,-1.2756627454888D+02,
     &  4.1487702227579D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  -1.7406770966429D+03,-9.3812204399266D+01, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-1.1890301282216D+03, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 2.3723447727360D+03, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-1.0279968223292D+03, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 5.7153838472603D+02, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00/)
      crest(:)=(/
     &  0.0000000000000D+00,-4.7430930170000D+00,-1.4422132560000D+01,
     &  -1.8061146510000D+01, 7.5186735000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  -2.7962099800000D+02, 1.7616414260000D+01,-9.9741392630000D+01,
     &  7.1402447000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00,-7.8571336480000D+01,
     &  5.2434353250000D+01, 7.7696745000000D+01, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  1.7799123760000D+02, 1.4564532380000D+02, 2.2347226000000D+02,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-4.3823284100000D+02,-7.2846553000000D+02,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00,-2.6752313750000D+02, 3.6170310000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     &  0.0000000000000D+00, 0.0000000000000D+00/)

      
      do i=1,245
        c5z_ps(i)=f5z*c5z_ps(i)+fbasis*cbasis(i)+fcore*ccore(i)
     &       +frest*crest(i)
      enddo
      c5z_ps(1) = c5z_ps(1)*2d0
      c5z_ps(:) = c5z_ps(:)*hartree/efreq

      end
