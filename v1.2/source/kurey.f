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
      subroutine kurey
      use angle
      use atmlst
      use atmtyp
      use domdec
      use inform
      use iounit
      use keys
      use kurybr
      use potent
      use urey   
      implicit none
      integer i,j,nu,nups,nuq,nutot
      integer ia,ib,ic
      integer ita,itb,itc
      integer size,next
      real*8 bb,tt
      logical header
      character*4 pa,pb,pc
      character*12 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
      logical max_reach
c
      blank = '            '
c
c     process keywords containing Urey-Bradley parameters
c
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
            bb = 0.0d0
            tt = 0.0d0
            string = record(next:240)
            read (string,*)  ia,ib,ic,bb,tt
            if (.not. silent) then
               if (header) then
                  header = .false.
                  if (rank.eq.0) write (iout,20)
   20             format (/,' Additional Urey-Bradley Parameters :',
     &                    //,5x,'Atom Classes',8x,'K(UB)',5x,
     &                       'Distance',/)
               end if
               if (rank.eq.0) write (iout,30)  ia,ib,ic,bb,tt
   30          format (4x,3i4,2x,f12.3,f12.4)
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
            do j = 1, maxnu
               if (ku(j).eq.blank .or. ku(j).eq.pt) then
                  ku(j) = pt
                  ucon(j) = bb
                  dst13(j) = tt
                  max_reach = .false.
                  exit
               end if
            end do
            if (max_reach) then
              if (rank.eq.0) write (iout,*)
     &          ' KUREY  --  Too many Urey-Bradley',
     &                 ' Interaction Parameters'
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
     &                ' Additional Angle repulsion Parameters :',
     &                   'Atom Classes','D','1/b'
               end if
               if (rank.eq.0) write (iout,'(4x,3i4,2x,f12.3,f12.4)')
     &             ia,ib,ic,bb,tt
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
            max_reach =.true.
            do j = 1, maxnups
               if (kups(j).eq.blank .or. kups(j).eq.pt) then
                  kups(j) = pt
                  uconps(j) = bb
                  dst13ps(j) = tt
                  max_reach=.false.
                  exit
               end if
            end do
            if (max_reach) then
              if (rank.eq.0) write (iout,*)
     &            ' KUREY  --  Too many Angle repulsion',
     &                 ' Interaction Parameters'
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
     &                ' Additional Angle repulsion Parameters :',
     &                   'Atom Classes','D','1/b'
               end if
               if (rank.eq.0) write (iout,'(4x,3i4,2x,f12.3,f12.4)')
     &             ia,ib,ic,bb,tt
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
            max_reach =.true.
            do j = 1, maxnuq
               if (kuq(j).eq.blank .or. kuq(j).eq.pt) then
                  kuq(j) = pt
                  uconq(j) = bb
                  dst13q(j) = tt
                  max_reach=.false.
                  exit
               end if
            end do
            if (max_reach) then
              if (rank.eq.0) write (iout,*)
     &            ' KUREY  --  Too many Quartic Urey',
     &                 ' Interaction Parameters'
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
         if (ku(i) .eq. blank)  nu = i - 1
      end do
      do i = maxnups, 1, -1
         if (kups(i) .eq. blank)  nups = i - 1
      end do
      do i = maxnuq, 1, -1
         if (kuq(i) .eq. blank)  nuq = i - 1
      end do
      nutot = nu + nups + nuq
c
c     assign the Urey-Bradley parameters for each angle
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
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, nu
               if (ku(j) .eq. pt) then
                  nurey = nurey + 1
                  iury(1,nurey) = ia
                  iury(2,nurey) = ib
                  iury(3,nurey) = ic
                  ureytyp(nurey) = 'UREYBRAD'
                  uk(nurey) = ucon(j)
                  ul(nurey) = dst13(j)
                  goto 60
               end if
            end do
            do j = 1, nups
               if (kups(j) .eq. pt) then
                  nurey = nurey + 1
                  iury(1,nurey) = ia
                  iury(2,nurey) = ib
                  iury(3,nurey) = ic
                  ureytyp(nurey) = 'ANGREP'
                  uk(nurey) = uconps(j)
                  ul(nurey) = dst13ps(j)
                  goto 60
               end if
            end do
            do j = 1, nuq
               if (kuq(j) .eq. pt) then
                  nurey = nurey + 1
                  iury(1,nurey) = ia
                  iury(2,nurey) = ib
                  iury(3,nurey) = ic
                  ureytyp(nurey) = 'UREYQUAR'
                  uk(nurey) = uconq(j)
                  ul(nurey) = dst13q(j)
                  goto 60
               end if
            end do
   60       continue
         end do
      end if
c
c     turn off the Urey-Bradley potential if it is not used
c
      if (nurey .eq. 0)  use_urey = .false.
      return
      end
c
c     subroutine kurey_update: update local urey bradley
c
      subroutine kurey_update
      use angle
      use atmlst
      use atmtyp
      use domdec
      use inform
      use iounit
      use keys
      use kurybr
      use potent
      use urey   
      implicit none
      integer i,j,nu,nups,nuq,nutot
      integer ia,ib,ic
      integer ita,itb,itc
      integer size
      integer iangle,ureycount,nureyloc1
      character*4 pa,pb,pc
      character*12 blank,pt
c
      blank = '            '

      nu = maxnu
      nups = maxnups
      nuq = maxnuq
      do i = maxnu, 1, -1
         if (ku(i) .eq. blank)  nu = i - 1
      end do
      do i = maxnups, 1, -1
          if (kups(i) .eq. blank)  nups = i - 1
      end do
       do i = maxnuq, 1, -1
          if (kuq(i) .eq. blank)  nuq = i - 1
      end do
      nutot = nu + nups + nuq
      
      if (allocated(ureyglob)) deallocate(ureyglob)
      allocate (ureyglob(nangleloc))
      nureyloc = 0
      do i = 1, nangleloc
        iangle = angleglob(i)
        ureycount = nburey(iangle)
        ia = iang(1,iangle)
        ib = iang(2,iangle)
        ic = iang(3,iangle)
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
        nureyloc1 = 0
        do j = 1, nu
           if (ku(j) .eq. pt) then
              nureyloc = nureyloc + 1
              nureyloc1 = nureyloc1 + 1
              ureyglob(nureyloc) = ureycount + nureyloc1
              goto 70
           end if
        end do
        do j = 1, nups
           if (kups(j) .eq. pt) then
              nureyloc = nureyloc + 1
              nureyloc1 = nureyloc1 + 1
              ureyglob(nureyloc) = ureycount + nureyloc1
              goto 70
           end if
        end do
        do j = 1, nuq
           if (kuq(j) .eq. pt) then
              nureyloc = nureyloc + 1
              nureyloc1 = nureyloc1 + 1
              ureyglob(nureyloc) = ureycount + nureyloc1
              goto 70
           end if
        end do
 70     continue
      end do
      return
      end
