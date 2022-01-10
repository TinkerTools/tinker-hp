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
      use urey   
      implicit none
      integer i,j,nu
      integer ia,ib,ic
      integer ita,itb,itc
      integer size,next
      integer iangle,ureycount,nureyloc1
      real*8 bb,tt
      logical header
      character*4 pa,pb,pc
      character*12 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
      logical init
c
      blank = '            '
      if (init) then
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
              string = record(next:120)
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
              size = 4
              call numeral (ia,pa,size)
              call numeral (ib,pb,size)
              call numeral (ic,pc,size)
              if (ia .le. ic) then
                 pt = pa//pb//pc
              else
                 pt = pc//pb//pa
              end if
              do j = 1, maxnu
                 if (ku(j).eq.blank .or. ku(j).eq.pt) then
                    ku(j) = pt
                    ucon(j) = bb
                    dst13(j) = tt
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KUREY  --  Too many Urey-Bradley',
     &                   ' Interaction Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        nu = maxnu
        do i = maxnu, 1, -1
           if (ku(i) .eq. blank)  nu = i - 1
        end do
c
c       assign the Urey-Bradley parameters for each angle
c
        nurey = 0
        if (nu .ne. 0) then
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
                    uk(nurey) = ucon(j)
                    ul(nurey) = dst13(j)
                    goto 60
                 end if
              end do
   60         continue
           end do
        end if
c
c       turn off the Urey-Bradley potential if it is not used
c
        if (nurey .eq. 0)  use_urey = .false.
      end if
      nu = maxnu
      do i = maxnu, 1, -1
         if (ku(i) .eq. blank)  nu = i - 1
      end do
      if (allocated(ureyglob)) deallocate(ureyglob)
      allocate (ureyglob(nangle))
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
 70     continue
      end do
      return
      end
