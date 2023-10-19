c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kbond  --  bond stretch parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kbond" assigns a force constant and ideal bond length
c     to each bond in the structure and processes any new or
c     changed parameter values
c
c
      subroutine kbond
      use atmtyp
      use bond
      use domdec
      use fields
      use inform
      use iounit
      use kbonds
      use keys
      use potent
      use usage
      use bndpot
      implicit none
      integer i,j
      integer ia,ib,ita,itb
      integer nb,nb5,nb4,nb3
      integer nbm,nbm4
      integer size,next
      integer minat,iring
      real*8 fc,bd,balpha
      logical header,done
      logical use_ring
      character*4 pa,pb
      character*6 label
      character*8 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
      logical :: max_reach
c
c     process keywords containing bond stretch parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:5) .eq. 'BOND ')  iring = 0
         if (keyword(1:6) .eq. 'BOND5 ')  iring = 5
         if (keyword(1:6) .eq. 'BOND4 ')  iring = 4
         if (keyword(1:6) .eq. 'BOND3 ')  iring = 3
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,fc,bd
   10       continue
            if (.not. silent) then
               if (header) then
                header = .false.
                if (rank.eq.0) write (iout,20)
   20           format (/,' Additional Bond Stretching Parameters :',
     &                  //,5x,'Atom Classes',9x,'K(S)',6x,'Length',/)
               end if
               if (iring .eq. 0) then
                  if (rank.eq.0) write (iout,30)  ia,ib,fc,bd
   30             format (6x,2i4,4x,f12.3,f12.4)
               else
                  if (iring .eq. 5)  label = '5-Ring'
                  if (iring .eq. 4)  label = '4-Ring'
                  if (iring .eq. 3)  label = '3-Ring'
                  if (rank.eq.0) write (iout,40)  ia,ib,fc,bd,label
   40             format (6x,2i4,4x,f12.3,f12.4,3x,a6)
               end if
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            if (iring .eq. 0) then
               max_reach = .true.
               do j = 1, maxnb
                  if (kb(j).eq.blank .or. kb(j).eq.pt) then
                     kb(j) = pt
                     bcon(j) = fc
                     blen(j) = bd
                     max_reach = .false.
                     exit
                  end if
               end do
               if (max_reach) then
                  if (rank.eq.0) write (iout,*)
     &              ' KBOND  --  Too many Bond Stretching',
     &                       ' Parameters'
                  abort = .true.
               end if
            else if (iring .eq. 5) then
               max_reach = .true.
               do j = 1, maxnb5
                  if (kb5(j).eq.blank .or. kb5(j).eq.pt) then
                     kb5(j) = pt
                     bcon5(j) = fc
                     blen5(j) = bd
                     max_reach = .false.
                     exit
                  end if
               end do
               if (max_reach) then
                  if (rank.eq.0) write (iout,*)
     &               ' KBOND  --  Too many 5-Ring Stretching',
     &                       ' Parameters'
                  abort = .true.
               end if
            else if (iring .eq. 4) then
               max_reach = .true.
               do j = 1, maxnb4
                  if (kb4(j).eq.blank .or. kb4(j).eq.pt) then
                     kb4(j) = pt
                     bcon4(j) = fc
                     blen4(j) = bd
                     max_reach = .false.
                     exit
                  end if
               end do
               if(max_reach) then
                  if (rank.eq.0) write (iout,*)
     &             ' KBOND  --  Too many 4-Ring Stretching',
     &                       ' Parameters'
                  abort = .true.
               endif
            else if (iring .eq. 3) then
               max_reach = .true.
               do j = 1, maxnb3
                  if (kb3(j).eq.blank .or. kb3(j).eq.pt) then
                     kb3(j) = pt
                     bcon3(j) = fc
                     blen3(j) = bd
                     max_reach = .false.
                     exit
                  end if
               end do
               if(max_reach) then
                  if (rank.eq.0) write (iout,*)
     &               ' KBOND  --  Too many 3-Ring Stretching',
     &                       ' Parameters'
                  abort = .true.
               endif
            end if
         end if
      end do
c
c     process keywords containing MORSE stretch parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:6) .eq. 'MORSE ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            balpha = 2.d0
            string = record(next:240)
            read (string,*)  ia,ib,fc,bd,balpha
            if (.not. silent) then
               if (header) then
                header = .false.
                if (rank.eq.0) write (iout,'(A,5x,A,9x,A,6x,A)') 
     &             ' Additional Morse Stretching Parameters :',
     &             'Atom Classes','K(S)','Length'
               end if
               if (rank.eq.0) 
     &               write (iout,'(6x,2i4,4x,f12.3,f12.4,f12.4)')
     &                       ia,ib,fc,bd,balpha
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, maxnbm
              max_reach=.true.
              if (kbm(j).eq.blank .or. kbm(j).eq.pt) then
                  kbm(j) = pt
                  bmor(1,j) = fc
                  bmor(2,j) = bd
                  bmor(3,j) = balpha
                  max_reach = .false.
                  exit
              end if
            end do
            if (max_reach) then
               if (rank.eq.0) write (iout,*)  
     &           'KBOND  --  Too many Morse Stretching Parameters'
               abort = .true.
            endif
         end if
      end do
c
c     process keywords containing MORSE4 stretch parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:7) .eq. 'MORSE4 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            balpha = 2.d0
            string = record(next:240)
            read (string,*)  ia,ib,fc,bd,balpha
            if (.not. silent) then
               if (header) then
                header = .false.
                if (rank.eq.0) write (iout,'(A,5x,A,9x,A,6x,A)') 
     &             ' Additional Morse4 Stretching Parameters :',
     &             'Atom Classes','K(S)','Length'
               end if
               if (rank.eq.0) 
     &               write (iout,'(6x,2i4,4x,f12.3,f12.4,f12.4)')
     &                       ia,ib,fc,bd,balpha
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, maxnbm4
              max_reach=.true.
              if (kbm4(j).eq.blank .or. kbm4(j).eq.pt) then
                  kbm4(j) = pt
                  bmor4(1,j) = fc
                  bmor4(2,j) = bd
                  bmor4(3,j) = balpha
                  max_reach = .false.
                  exit
              end if
            end do
            if (max_reach) then
               if (rank.eq.0) write (iout,*)  
     &           'KBOND  --  Too many Morse4 Stretching Parameters'
               abort = .true.
            endif
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nb = maxnb
      nbm = maxnbm
      nbm4 = maxnbm4
      nb5 = maxnb5
      nb4 = maxnb4
      nb3 = maxnb3
      do i = maxnb, 1, -1
         if (kb(i) .eq. blank)  nb = i - 1
      end do
      do i = maxnb5, 1, -1
         if (kb5(i) .eq. blank)  nb5 = i - 1
      end do
      do i = maxnb4, 1, -1
         if (kb4(i) .eq. blank)  nb4 = i - 1
      end do
      do i = maxnb3, 1, -1
         if (kb3(i) .eq. blank)  nb3 = i - 1
      end do
      do i = maxnbm, 1, -1
         if (kbm(i) .eq. blank)  nbm = i - 1
      end do
      do i = maxnbm4, 1, -1
         if (kbm4(i) .eq. blank)  nbm4 = i - 1
      end do
      use_ring = .false.
      if (min(nb5,nb4,nb3) .ne. 0)  use_ring = .true.
c
c     use special bond parameter assignment method for MMFF
c
      if (forcefield .eq. 'MMFF94') then
         call kbondm
         return
      end if
c
c     assign ideal bond length and force constant for each bond
c
      header = .true.
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt = pa//pb
         else
            pt = pb//pa
         end if
         bk(i) = 0.0d0
         bl(i) = 0.0d0
         ba(i) = 2.0d0
         bndtyp(i) = bndtyp_default
         done = .false.
c
c     make a check for bonds contained inside small rings
c
         iring = 0
         if (use_ring) then
            call chkring (iring,ia,ib,0,0)
            if (iring .eq. 6)  iring = 0
            if (iring.eq.5 .and. nb5.eq.0)  iring = 0
            if (iring.eq.4 .and. nb4.eq.0)  iring = 0
            if (iring.eq.3 .and. nb3.eq.0)  iring = 0
         end if
c
c     assign bond stretching parameters for each bond
c
         if (iring .eq. 0) then
            do j = 1, nb
               if (kb(j) .eq. pt) then
                  bk(i) = bcon(j)
                  bl(i) = blen(j)
                  done = .true.
                  goto 130
               end if
            end do
c
c     assign stretching parameters for 5-membered ring bonds
c
         else if (iring .eq. 5) then
            do j = 1, nb5
               if (kb5(j) .eq. pt) then
                  bk(i) = bcon5(j)
                  bl(i) = blen5(j)
                  done = .true.
                  goto 130
               end if
            end do
c
c     assign stretching parameters for 4-membered ring bonds
c
         else if (iring .eq. 4) then
            do j = 1, nb4
               if (kb4(j) .eq. pt) then
                  bk(i) = bcon4(j)
                  bl(i) = blen4(j)
                  done = .true.
                  goto 130
               end if
            end do
c
c     assign stretching parameters for 3-membered ring bonds
c
         else if (iring .eq. 3) then
            do j = 1, nb3
               if (kb3(j) .eq. pt) then
                  bk(i) = bcon3(j)
                  bl(i) = blen3(j)
                  done = .true.
                  goto 130
               end if
            end do
         end if
c
c     assign stretching parameters for Morse bonds
c
         if (.not. done) then
            do j = 1, nbm
               if (kbm(j) .eq. pt) then
                  bk(i) = bmor(1,j)
                  bl(i) = bmor(2,j)
                  ba(i) = bmor(3,j)
                  bndtyp(i) = 'MORSE'
                  done = .true.
                  goto 130
               end if
            end do
         end if
c
c     assign stretching parameters for Morse4 bonds
c
         if (.not. done) then
            do j = 1, nbm4
               if (kbm4(j) .eq. pt) then
                  bk(i) = bmor4(1,j)
                  bl(i) = bmor4(2,j)
                  ba(i) = bmor4(3,j)
                  bndtyp(i) = 'MORSE4'
                  done = .true.
                  goto 130
               end if
            end do
         end if
c
c     warning if suitable bond stretching parameter not found
c
  130    continue
         minat = min(atomic(ia),atomic(ib))
         if (minat .eq. 0)  done = .true.
         if (use_bond .and. .not.done) then
            if (use(ia) .or. use(ib))  abort = .true.
            if (header) then
               header = .false.
               if (rank.eq.0) write (iout,140)
  140          format (/,' Undefined Bond Stretching Parameters :',
     &                 //,' Type',13x,'Atom Names',11x,
     &                    'Atom Classes',/)
            end if
            label = 'Bond  '
            if (iring .eq. 5)  label = '5-Ring'
            if (iring .eq. 4)  label = '4-Ring'
            if (iring .eq. 3)  label = '3-Ring'
            if (rank.eq.0) write (iout,150)
     $        label,ia,name(ia),ib,name(ib),ita
     $        ,itb
  150       format (1x,a6,5x,i6,'-',a3,i6,'-',a3,7x,2i5)
         end if
      end do
c
c     check for electronegativity bond length corrections
c
      if (hostrank.eq.0) call keneg
c
c     turn off the bond stretch potential if it is not used
c
      if (nbond .eq. 0)  use_bond = .false.
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine keneg  --  assign electronegativity parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "keneg" applies primary and secondary electronegativity bond
c     length corrections to applicable bond parameters
c
c     note this version does not scale multiple corrections to the
c     same bond by increasing powers of 0.62 as in MM3
c
c
      subroutine keneg
      use angle
      use atmlst
      use atmtyp
      use bond
      use couple
      use kbonds
      use keys
      use inform
      use iounit
      use tors
      implicit none
      integer i,j,k,m,nel
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer size,next
      real*8 dl,factor
      logical header
      character*4 pa,pb,pc,pd
      character*12 blank
      character*12 pt,pt1,pt2
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing electronegativity parameters
c
      blank = '            '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'ELECTNEG ') then
            ia = 0
            ib = 0
            ic = 0
            dl = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,dl
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Electronegativity',
     &                       ' Parameters :',
     &                    //,5x,'Atom Classes',18x,'dLength',/)
               end if
               write (iout,30)  ia,ib,ic,dl
   30          format (4x,3i4,14x,f12.4)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            pt = pa//pb//pc
            do j = 1, maxnel
               if (kel(j).eq.blank .or. kel(j).eq.pt) then
                  kel(j) = pt
                  dlen(j) = dl
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KENEG  --  Too many Electronegativity',
     &                    ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nel = maxnel
      do i = 1, maxnel
         if (kel(i) .eq. blank) then
            nel = i - 1
            goto 60
         end if
      end do
   60 continue
c
c     check angles for primary electronegativity corrections
c
      if (nel .ne. 0) then
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
            pt1 = pa//pb//pc
            pt2 = pc//pb//pa
c
c     search the parameter set for a match to either bond
c
            do j = 1, nel
               if (kel(j) .eq. pt1) then
                  do k = 1, n12(ia)
                     if (i12(k,ia) .eq. ib) then
                        m = bndlist(k,ia)
                        bl(m) = bl(m) + dlen(j)
                     end if
                  end do
                  goto 70
               else if (kel(j) .eq. pt2) then
                  do k = 1, n12(ic)
                     if (i12(k,ic) .eq. ib) then
                        m = bndlist(k,ic)
                        bl(m) = bl(m) + dlen(j)
                     end if
                  end do
                  goto 70
               end if
            end do
   70       continue
         end do
c
c     check torsions for secondary electronegativity corrections
c
         factor = 0.4d0
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            pt1 = pa//pb//pd
            pt2 = pd//pc//pa
c
c     turn off electronegativity effect for attached hydrogens
c
            if (atomic(id) .le. 1)  pt1 = blank
            if (atomic(ia) .le. 1)  pt2 = blank
c
c     search the parameter set for a match to either bond
c
            do j = 1, nel
               if (kel(j) .eq. pt1) then
                  do k = 1, n12(ia)
                     if (i12(k,ia) .eq. ib) then
                        m = bndlist(k,ia)
                        bl(m) = bl(m) + factor*dlen(j)
                     end if
                  end do
                  goto 80
               else if (kel(j) .eq. pt2) then
                  do k = 1, n12(id)
                     if (i12(k,id) .eq. ic) then
                        m = bndlist(k,id)
                        bl(m) = bl(m) + factor*dlen(j)
                     end if
                  end do
                  goto 80
               end if
            end do
   80       continue
         end do
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine kbondm  --  assign MMFF bond stretch parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "kbondm" assigns a force constant and ideal bond length to
c     each bond according to the Merck Molecular Force Field (MMFF)
c
c     literature reference:
c
c     R. Blom and A. Haaland, "A Modification of the Schomaker-Stevenson
c     Rule for Prediction of Single Bond Distances", Journal of
c     Molecular Structure, 128, 21-27 (1985)
c
c
      subroutine kbondm
      use sizes
      use atmtyp
      use bond
      use keys
      use merck
      implicit none
      integer i,j
      integer ia,ib,ita,itb
      integer next,minat
      integer list(20)
      real*8 khia,khib,cst
      real*8 rad0a,rad0b
      logical header,done
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     get single bonds that could be double (MMFF bond type=1)
c
      nlignes = 0
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:12) .eq. 'MMFF-PIBOND ') then
            do j = 1, 20
               list(j) = 0
            end do
            string = record(next:240)
            read (string,*,err=10,end=10)  (list(j),j=1,20)
   10       continue
            do j = 1, 20, 2
               if (list(j).ne.0 .and. list(j+1).ne.0) then
                  nlignes = nlignes + 1
                  bt_1(nlignes,1) = list(j)
                  bt_1(nlignes,2) = list(j+1)
               else
                  goto 20
               end if
            end do
   20       continue
         end if
      end do
c
c     assign MMFF bond length and force constant values
c
      header = .true.
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         if (ia .le. ib) then
            do j = 1, nlignes
               if (ia.eq.bt_1(j,1) .and. ib.eq.bt_1(j,2)) then
                  bk(i) = mmff_kb1(ita,itb)
                  bl(i) = mmff_b1(ita,itb)
                  done = .true.
                  if (bk(i) .eq. 1000.0d0)  done = .false.
                  if (bl(i) .eq. 1000.0d0)  done = .false.
                  goto 30
               end if
            end do
            bk(i) = mmff_kb(ita,itb)
            bl(i) = mmff_b0(ita,itb)
            done = .true.
            if (bk(i) .eq. 1000.0d0)  done = .false.
            if (bl(i) .eq. 1000.0d0)  done = .false.
            goto 30
         else if (ib .le. ia) then
            do j = 1, nlignes
               if (ib.eq.bt_1(j,1) .and. ia.eq.bt_1(j,2)) then
                  bk(i) = mmff_kb1(itb,ita)
                  bl(i) = mmff_b1(itb,ita)
                  done = .true.
                  if (bk(i) .eq. 1000.0d0)  done = .false.
                  if (bl(i) .eq. 1000.0d0)  done = .false.
                  goto 30
               end if
            end do
            bk(i) = mmff_kb(itb,ita)
            bl(i) = mmff_b0(itb,ita)
            done = .true.
            if (bk(i) .eq. 1000.0d0)  done = .false.
            if (bl(i) .eq. 1000.0d0)  done = .false.
            goto 30
         end if
c
c     estimate missing bond parameters via an empirical rule
c
   30    continue
         minat = min(atomic(ia),atomic(ib))
         if (minat .eq. 0)  done = .true.
         if (.not. done) then
            khia = paulel(atomic(ia))
            khib = paulel(atomic(ib))
            rad0a = rad0(atomic(ia))
            rad0b = rad0(atomic(ib))
            cst = 0.085d0
            if (atomic(ia).eq.1 .or. atomic(ib).eq.1)  cst = 0.05d0
            bl(i) = rad0a + rad0b - cst*abs(khia-khib)**1.4d0
            bk(i) = kbref(atomic(ia),atomic(ib))
     &                 * (r0ref(atomic(ia),atomic(ib))/bl(i))**6
         end if
      end do
c
c     turn off the bond stretch potential if it is not used
c
c      if (nbond .eq. 0)  use_bond = .false.
      return
      end
