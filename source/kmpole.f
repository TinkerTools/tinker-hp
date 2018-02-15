c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kmpole  --  multipole parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kmpole" assigns atomic multipole moments to the atoms of
c     the structure and processes any new or changed values
c
c
      subroutine kmpole(init,istep)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kmulti.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'potent.i'
      include 'units.i'
      include 'neigh.i'
      include 'openmp.i'
      include 'mpif.h'
      integer istep,modnl,ierr
      integer i,j,k,l,m,ia
      integer iproc,iglob,polecount
      integer ji,ki,li
      integer it,jt,kt,lt
      integer imp,nmp
      integer size,next
      integer number
      integer kz,kx,ky
      integer ztyp,xtyp,ytyp
      integer, allocatable :: mpt(:)
      integer, allocatable :: mpz(:)
      integer, allocatable :: mpx(:)
      integer, allocatable :: mpy(:)
      real*8 mpl(13)
      real*8 emtp1,emtp2,emtp3
      real*8 d
      logical header,path
      character*4 pa,pb,pc,pd
      character*8 axt
      character*16 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
      logical init
c
      if (init) then
c
c     allocate global arrays
c
        if (associated(rpole)) deallocate (rpole)
        allocate (rpole(13,n))
        if (associated(polelocnl)) deallocate(polelocnl)
        allocate (polelocnl(n))
        if (associated(poleloc)) deallocate(poleloc)
        allocate (poleloc(n))
        if (associated(polerecloc)) deallocate(polerecloc)
        allocate (polerecloc(n))
        if (.not.associated(uind)) allocate (uind(3,n))
        if (.not.associated(uinp)) allocate (uinp(3,n))
        uind = 0d0
        uinp = 0d0
c
c     allocate global arrays
c
        call alloc_shared_mpole
c
c       count the number of existing multipole parameters
c
        if (hostrank.ne.0) goto 230
        blank = '                '
        nmp = maxnmp
        do i = maxnmp, 1, -1
           if (kmp(i) .eq. blank)  nmp = i - 1
        end do
c
c       find and count new multipole parameters in the keyfile
c
        imp = 0
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:10) .eq. 'MULTIPOLE ') then
              k = 0
              string = record(next:120)
              read (string,*,err=10,end=10)  k,kz,kx,ky,mpl(1)
              goto 40
   10         continue
              read (string,*,err=20,end=20)  k,kz,kx,mpl(1)
              goto 40
   20         continue
              read (string,*,err=30,end=30)  k,kz,mpl(1)
              goto 40
   30         continue
              read (string,*,err=50,end=50)  k,mpl(1)
   40         continue
              if (k .gt. 0) then
                 record = keyline(i+1)
                 read (record,*,err=50,end=50)  mpl(2),mpl(3),mpl(4)
                 record = keyline(i+2)
                 read (record,*,err=50,end=50)  mpl(5)
                 record = keyline(i+3)
                 read (record,*,err=50,end=50)  mpl(8),mpl(9)
                 record = keyline(i+4)
                 read (record,*,err=50,end=50)  mpl(11),mpl(12),mpl(13)
                 imp = imp + 1
              end if
   50         continue
           end if
         end do
c
c       check for too many combined parameter values
c
        nmp = nmp + imp
        if (nmp .gt. maxnmp) then
           if (rank.eq.0) write (iout,60)
   60      format (/,' KMPOLE  --  Too many Atomic Multipole',
     &                ' Parameters')
           abort = .true.
        end if
c
c       move existing parameters to make room for new values
c
        if (imp .ne. 0) then
           do j = nmp, imp+1, -1
              k = j - imp
              kmp(j) = kmp(k)
              mpaxis(j) = mpaxis(k)
              do m = 1, 13
                 multip(m,j) = multip(m,k)
              end do
           end do
        end if
c
c       process keywords containing atomic multipole parameters
c
        imp = 0
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:10) .eq. 'MULTIPOLE ') then
              k = 0
              kz = 0
              kx = 0
              ky = 0
              axt = 'Z-then-X'
              do j = 1, 13
                 mpl(j) = 0.0d0
              end do
              string = record(next:120)
              read (string,*,err=70,end=70)  k,kz,kx,ky,mpl(1)
              goto 100
   70         continue
              ky = 0
              read (string,*,err=80,end=80)  k,kz,kx,mpl(1)
              goto 100
   80         continue
              kx = 0
              read (string,*,err=90,end=90)  k,kz,mpl(1)
              goto 100
   90         continue
              kz = 0
              read (string,*,err=130,end=130)  k,mpl(1)
  100         continue
              if (k .gt. 0) then
               if (kz .eq. 0)  axt = 'None'
               if (kz.ne.0 .and. kx.eq.0)  axt = 'Z-Only'
               if (kz.lt.0 .or. kx.lt.0)  axt = 'Bisector'
               if (kx.lt.0 .and. ky.lt.0)  axt = 'Z-Bisect'
               if (max(kz,kx,ky) .lt. 0)  axt = '3-Fold'
               kz = abs(kz)
               kx = abs(kx)
               ky = abs(ky)
               record = keyline(i+1)
               read (record,*,err=130,end=130)  mpl(2),mpl(3),mpl(4)
               record = keyline(i+2)
               read (record,*,err=130,end=130)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=130,end=130)  mpl(8),mpl(9)
               record = keyline(i+4)
               read (record,*,err=130,end=130)  mpl(11),mpl(12),mpl(13)
               mpl(6) = mpl(8)
               mpl(7) = mpl(11)
               mpl(10) = mpl(12)
               if (header .and. .not.silent) then
                  header = .false.
                  if (rank.eq.0) write (iout,110)
  110             format (/,' Additional Atomic Multipole Parameters :',
     &                    //,5x,'Atom Type',5x,'Coordinate Frame',
     &                       ' Definition',8x,'Multipole Moments')
               end if
               if (.not. silent) then
                  if (rank.eq.0) write (iout,120)  k,kz,kx,ky,axt,
     &                             (mpl(j),j=1,5),
     &                             mpl(8),mpl(9),(mpl(j),j=11,13)
  120             format (/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,
     &                       /,48x,3f9.5,/,48x,f9.5,
     &                       /,48x,2f9.5,/,48x,3f9.5)
               end if
               size = 4
               call numeral (k,pa,size)
               call numeral (kz,pb,size)
               call numeral (kx,pc,size)
               call numeral (ky,pd,size)
               pt = pa//pb//pc//pd
               imp = imp + 1
               kmp(imp) = pt
               mpaxis(imp) = axt
               do j = 1, 13
                  multip(j,imp) = mpl(j)
               end do
              end if
  130         continue
           end if
        end do
c
c    find and count new sibfacp parameters in the keyfile
c
        do i = 1, nkey
          next = 1
          record = keyline(i)
          call gettext (record,keyword,next)
          call upcase (keyword)
          if (keyword(1:8) .eq. 'SIBFACP ') then
            ia = 0
            emtp1 = 0.0d0
            emtp2 = 0.0d0
            emtp3 = 0.0d0
            string = record(next:120)
            read (string,*,err=220,end=220) ia,emtp1,emtp2,emtp3
            if (ia.ne.0) then
              sibfacp(1,ia) = emtp1
              sibfacp(2,ia) = emtp2
              sibfacp(3,ia) = emtp3
            end if
  220       continue
          end if
        end do
c
c       zero out local axes, multipoles and polarization attachments
c
        npole = n
        do i = 1, n
           pollist(i) = 0
           zaxis(i) = 0
           xaxis(i) = 0
           yaxis(i) = 0
           polaxe(i) = '        '
           do j = 1, 13
              pole(j,i) = 0.0d0
           end do
           np11(i) = 0
           np12(i) = 0
           np13(i) = 0
           np14(i) = 0
           if (use_emtp) then
             alphapen(i) = 0.0d0
             betapen(i)  = 0.0d0
             gammapen(i) = 0.0d0
           end if
        end do
c
c     assign sibfacp parameters
c
        if (use_emtp) then
          do i = 1, n
            it = type(i)
            do k = 1, maxtyp
              if (it.eq.k) then
                alphapen(i) = sibfacp(1,k)
                betapen(i)  = sibfacp(2,k)
                gammapen(i) = sibfacp(3,k)
              end if
            end do
          end do
        end if
c
c       perform dynamic allocation of some local arrays
c
        allocate (mpt(maxnmp))
        allocate (mpz(maxnmp))
        allocate (mpx(maxnmp))
        allocate (mpy(maxnmp))
c
c       store the atom types associated with each parameter
c
        do i = 1, nmp
           mpt(i) = number(kmp(i)(1:4))
           mpz(i) = number(kmp(i)(5:8))
           mpx(i) = number(kmp(i)(9:12))
           mpy(i) = number(kmp(i)(13:16))
        end do
c
c       assign multipole parameters via only 1-2 connected atoms
c
        do i = 1, n
           it = type(i)
           do imp = 1, nmp
              if (it .eq. mpt(imp)) then
                 ztyp = mpz(imp)
                 xtyp = mpx(imp)
                 ytyp = mpy(imp)
                 do j = 1, n12(i)
                    ji = i12(j,i)
                    jt = type(ji)
                    if (jt .eq. ztyp) then
                       do k = 1, n12(i)
                          ki = i12(k,i)
                          kt = type(ki)
                          if (kt.eq.xtyp .and. ki.ne.ji) then
                             if (ytyp .eq. 0) then
                                zaxis(i) = ji
                                xaxis(i) = ki
                                polaxe(i) = mpaxis(imp)
                                do m = 1, 13
                                   pole(m,i) = multip(m,imp)
                                end do
                                goto 140
                             end if
                             do l = 1, n12(i)
                                li = i12(l,i)
                                lt = type(li)
                                if (lt.eq.ytyp .and. li.ne.ji
     &                                 .and. li.ne.ki) then
                                   zaxis(i) = ji
                                   xaxis(i) = ki
                                   yaxis(i) = li
                                   polaxe(i) = mpaxis(imp)
                                   do m = 1, 13
                                      pole(m,i) = multip(m,imp)
                                   end do
                                   goto 140
                                end if
                             end do
                          end if
                       end do
                    end if
                 end do
              end if
           end do
c
c       assign multipole parameters via 1-2 and 1-3 connected atoms
c
           do imp = 1, nmp
              if (it .eq. mpt(imp)) then
                 ztyp = mpz(imp)
                 xtyp = mpx(imp)
                 ytyp = mpy(imp)
                 do j = 1, n12(i)
                    ji = i12(j,i)
                    jt = type(ji)
                    if (jt .eq. ztyp) then
                       do k = 1, n13(i)
                          ki = i13(k,i)
                          kt = type(ki)
                          path = .false.
                          do m = 1, n12(ki)
                             if (i12(m,ki) .eq. ji)  path = .true.
                          end do
                          if (kt.eq.xtyp .and. path) then
                             if (ytyp .eq. 0) then
                                zaxis(i) = ji
                                xaxis(i) = ki
                                polaxe(i) = mpaxis(imp)
                                do m = 1, 13
                                   pole(m,i) = multip(m,imp)
                                end do
                                goto 140
                             end if
                             do l = 1, n13(i)
                                li = i13(l,i)
                                lt = type(li)
                                path = .false.
                                do m = 1, n12(li)
                                   if (i12(m,li) .eq. ji)  path = .true.
                                end do
                                if (lt.eq.ytyp .and. li.ne.ki
     &                                 .and. path) then
                                   zaxis(i) = ji
                                   xaxis(i) = ki
                                   yaxis(i) = li
                                   polaxe(i) = mpaxis(imp)
                                   do m = 1, 13
                                      pole(m,i) = multip(m,imp)
                                   end do
                                   goto 140
                                end if
                             end do
                          end if
                       end do
                    end if
                 end do
              end if
           end do
c
c       assign multipole parameters via only a z-defining atom
c
           do imp = 1, nmp
              if (it .eq. mpt(imp)) then
                 ztyp = mpz(imp)
                 xtyp = mpx(imp)
                 ytyp = mpy(imp)
                 do j = 1, n12(i)
                    ji = i12(j,i)
                    jt = type(ji)
                    if (jt .eq. ztyp) then
                       if (xtyp .eq. 0) then
                          zaxis(i) = ji
                          polaxe(i) = mpaxis(imp)
                          do m = 1, 13
                             pole(m,i) = multip(m,imp)
                          end do
                          goto 140
                       end if
                    end if
                 end do
              end if
           end do
c
c       assign multipole parameters via no connected atoms
c
           do imp = 1, nmp
              if (it .eq. mpt(imp)) then
                 ztyp = mpz(imp)
                 xtyp = mpx(imp)
                 ytyp = mpy(imp)
                 if (ztyp .eq. 0) then
                    polaxe(i) = mpaxis(imp)
                    do m = 1, 13
                       pole(m,i) = multip(m,imp)
                    end do
                    goto 140
                 end if
              end if
           end do
  140      continue
        end do
c
c       perform deallocation of some local arrays
c
        deallocate (mpt)
        deallocate (mpz)
        deallocate (mpx)
        deallocate (mpy)
c
c       process keywords with multipole parameters for specific atoms
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:10) .eq. 'MULTIPOLE ') then
              k = 0
              kz = 0
              kx = 0
              ky = 0
              axt = 'Z-then-X'
              do j = 1, 13
                 mpl(j) = 0.0d0
              end do
              string = record(next:120)
              read (string,*,err=150,end=150)  k,kz,kx,ky,mpl(1)
              goto 180
  150         continue
              ky = 0
              read (string,*,err=160,end=160)  k,kz,kx,mpl(1)
              goto 180
  160         continue
              kx = 0
              read (string,*,err=170,end=170)  k,kz,mpl(1)
              goto 180
  170         continue
              kz = 0
              read (string,*,err=210,end=210)  k,mpl(1)
  180         continue
              if (k.lt.0 .and. k.ge.-n) then
                k = -k
                if (kz .eq. 0)  axt = 'None'
                if (kz.ne.0 .and. kx.eq.0)  axt = 'Z-Only'
                if (kz.lt.0 .or. kx.lt.0)  axt = 'Bisector'
                if (kx.lt.0 .and. ky.lt.0)  axt = 'Z-Bisect'
                if (max(kz,kx,ky) .lt. 0)  axt = '3-Fold'
                kz = abs(kz)
                kx = abs(kx)
                ky = abs(ky)
                record = keyline(i+1)
                read (record,*,err=210,end=210)  mpl(2),mpl(3),mpl(4)
                record = keyline(i+2)
                read (record,*,err=210,end=210)  mpl(5)
                record = keyline(i+3)
                read (record,*,err=210,end=210)  mpl(8),mpl(9)
                record = keyline(i+4)
                read (record,*,err=210,end=210)  mpl(11),mpl(12),mpl(13)
                mpl(6) = mpl(8)
                mpl(7) = mpl(11)
                mpl(10) = mpl(12)
                if (header .and. .not.silent) then
                   header = .false.
                   if (rank.eq.0) write (iout,190)
  190              format (/,' Additional Atomic Multipoles',
     &                        ' for Specific Atoms :',
     &                     //,6x,'Atom',9x,'Coordinate Frame',
     &                        ' Definition',8x,'Multipole Moments')
                end if
                if (.not. silent) then
                   if (rank.eq.0) write (iout,200)  k,kz,kx,ky,axt,
     &                               (mpl(j),j=1,5),
     &                               mpl(8),mpl(9),(mpl(j),j=11,13)
  200              format (/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,
     &                        /,48x,3f9.5,/,48x,f9.5,
     &                        /,48x,2f9.5,/,48x,3f9.5)
                end if
                zaxis(k) = kz
                xaxis(k) = kx
                yaxis(k) = ky
                polaxe(k) = axt
                do j = 1, 13
                   pole(j,k) = mpl(j)
                end do
              end if
  210         continue
           end if
        end do
c
c       convert the dipole and quadrupole moments to Angstroms,
c       quadrupole divided by 3 for use as traceless values
c
        do i = 1, n
           do k = 2, 4
              pole(k,i) = pole(k,i) * bohr
           end do
           do k = 5, 13
              pole(k,i) = pole(k,i) * bohr**2 / 3.0d0
           end do
        end do
 230    call MPI_BARRIER(hostcomm,ierr)
c
c       get the order of the multipole expansion at each site
c
        do i = 1, n
           size = 0
           do k = 1, maxpole
              if (pole(k,i) .ne. 0.0d0)  size = max(k,size)
           end do
           if (size .gt. 4) then
              size = 13
           else if (size .gt. 1) then
              size = 4
           end if
           polsiz(i) = size
        end do
c
c       if polarization not used, zero out induced dipoles
c
        if (.not. use_polar) then
           do i = 1, n
              do j = 1, 3
                 uind(j,i) = 0.0d0
                 uinp(j,i) = 0.0d0
c                 uinds(j,i) = 0.0d0
c                 uinps(j,i) = 0.0d0
              end do
           end do
        end if
c
c       remove any zero or undefined atomic multipoles
c
        if (.not.use_polar .and. .not.use_solv) then
           npole = 0
           do i = 1, n
              if (polsiz(i) .ne. 0) then
                 nbpole(i) = npole
                 npole = npole + 1
                 if (hostrank.eq.0) then
                   ipole(npole) = i
                   pollist(i) = npole
                   zaxis(npole) = zaxis(i)
                   xaxis(npole) = xaxis(i)
                   yaxis(npole) = yaxis(i)
                   polaxe(npole) = polaxe(i)
                   do j = 1, maxpole
                      pole(j,npole) = pole(j,i)
                   end do
                   if (use_emtp) then
                     alphapen(npole) = alphapen(i)
                     betapen(npole) = betapen(i)
                     gammapen(npole) = gammapen(i)
                   end if
                 end if
              end if
           end do
c
c       test multipoles at chiral sites and invert if necessary
c
           call chkpole
c
c       turn off the atomic multipole potential if it is not used
c
           if (npole .eq. 0)  use_mpole = .false.
        end if
      end if
c
      if (associated(poleglob)) deallocate(poleglob)
      allocate (poleglob(nbloc))
      if (associated(polerecglob)) deallocate(polerecglob)
      allocate (polerecglob(nlocrec))
c
c       remove any zero or undefined atomic multipoles
c
      if (.not.use_polar .and. .not.use_solv) then
         npoleloc = 0
         do i = 1, nloc
            iglob = glob(i)
            polecount = nbpole(iglob)
            if (polsiz(iglob) .ne. 0) then
               npoleloc = npoleloc + 1
               poleglob(npoleloc) = polecount + 1
               poleloc(polecount+1) = npoleloc
            end if
         end do
        npolebloc = npoleloc
        do iproc = 1, n_recep1
          if (domlen(p_recep1(iproc)+1).ne.0) then
            bufbegpole(p_recep1(iproc)+1) = npolebloc + 1
          else
            bufbegpole(p_recep1(iproc)+1) = 1
          end if
          do i = 1, domlen(p_recep1(iproc)+1)
            iglob = glob(bufbeg(p_recep1(iproc)+1)+i-1)
            polecount = nbpole(iglob)
            if (polsiz(iglob) .ne. 0.0d0) then
              npolebloc = npolebloc + 1
              poleglob(npolebloc) = polecount + 1
              poleloc(polecount+1) = npolebloc
            end if
          end do
        end do
c
c  deal with case where ewald summation is used
c
        npolerecloc = 0
        do i = 1, nlocrec
           iglob = globrec(i)
           polecount = nbpole(iglob)
           if (polsiz(iglob) .ne. 0) then
              npolerecloc = npolerecloc + 1
              polerecglob(npolerecloc) = polecount + 1
              polerecloc(polecount+1) = npolerecloc
           end if
        end do
c
        modnl = mod(istep,ineigup)
        if (modnl.ne.0) return
        if (associated(poleglobnl)) deallocate(poleglobnl)
        allocate (poleglobnl(nlocnl))
c
        npolelocnl = 0
        do i = 1, nlocnl
          iglob = ineignl(i)
          polecount = nbpole(iglob)
          if (polsiz(iglob) .ne. 0 ) then
            call distprocpart(iglob,rank,d,.true.)
            if (repart(iglob).eq.rank) d = 0.0d0
            if (d*d.le.(mbuf2/4)) then
              npolelocnl = npolelocnl + 1
              poleglobnl(npolelocnl) = polecount + 1
              polelocnl(polecount+1) = npolelocnl
            end if
          end if
        end do
      end if
c
      return
      end
c
c     subroutine alloc_shared_mpole : allocate shared memory pointers for mpole
c     parameter arrays
c
      subroutine alloc_shared_mpole
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use mpi
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polgrp.i'
      include 'potent.i'
      include 'openmp.i'
      integer :: win,win2
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
      if (associated(polsiz)) deallocate (polsiz)
      if (associated(zaxis)) deallocate (zaxis)
      if (associated(yaxis)) deallocate (yaxis)
      if (associated(xaxis)) deallocate (xaxis)
      if (associated(pole)) deallocate (pole)
      if (associated(polaxe)) deallocate (polaxe)
      if (associated(ipole)) deallocate (ipole)
      if (associated(pollist)) deallocate (pollist)
      if (associated(nbpole)) deallocate (nbpole)
      if (use_emtp) then
        if (associated(alphapen)) deallocate (alphapen)
        if (associated(betapen)) deallocate (betapen)
        if (associated(gammapen)) deallocate (gammapen)
      end if
c
c     polsiz
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,polsiz,arrayshape)
c
c     zaxis
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,zaxis,arrayshape)
c
c     yaxis
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,yaxis,arrayshape)
c
c     xaxis
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,xaxis,arrayshape)
c
c     pole
c
      arrayshape2=(/13,n/)
      if (hostrank == 0) then
        windowsize = int(13*n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pole,arrayshape2)
c
c     polaxe
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,polaxe,arrayshape)
c
c     ipole
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ipole,arrayshape)
c
c     pollist
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pollist,arrayshape)
c
c     nbpole
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbpole,arrayshape)
      if (use_emtp) then
c
c     alphapen
c
        arrayshape=(/n/)
        if (hostrank == 0) then
          windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
        else
          windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
c
c    allocation
c
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $    MPI_INFO_NULL,hostcomm, baseptr, win, ierr)
        if (hostrank /= 0) then
          CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $    baseptr, ierr)
        end if
c
c    association with fortran pointer
c
        CALL C_F_POINTER(baseptr,alphapen,arrayshape)
c
c     betapen
c
        arrayshape=(/n/)
        if (hostrank == 0) then
          windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
        else
          windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
c
c    allocation
c
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $    MPI_INFO_NULL,hostcomm, baseptr, win, ierr)
        if (hostrank /= 0) then
          CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $    baseptr, ierr)
        end if
c
c    association with fortran pointer
c
        CALL C_F_POINTER(baseptr,betapen,arrayshape)
c
c     gammapen
c
        arrayshape=(/n/)
        if (hostrank == 0) then
          windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
        else
          windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
c
c    allocation
c
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $    MPI_INFO_NULL,hostcomm, baseptr, win, ierr)
        if (hostrank /= 0) then
          CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $    baseptr, ierr)
        end if
c
c    association with fortran pointer
c
        CALL C_F_POINTER(baseptr,gammapen,arrayshape)
      end if
c
      if (associated(ip11)) deallocate(ip11)
      if (associated(np11)) deallocate (np11)
      if (associated(ip12)) deallocate (ip12)
      if (associated(np12)) deallocate (np12)
      if (associated(ip13)) deallocate (ip13)
      if (associated(np13)) deallocate (np13)
      if (associated(ip14)) deallocate (ip14)
      if (associated(np14)) deallocate (np14)
c
c     ip11
c
      arrayshape2=(/maxp11,n/)
      if (hostrank == 0) then
        windowsize = int(maxp11*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ip11,arrayshape2)
c
c     np11
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,np11,arrayshape)
c
c     ip12
c
      arrayshape2=(/maxp12,n/)
      if (hostrank == 0) then
        windowsize = int(maxp12*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ip12,arrayshape2)
c
c     np12
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,np12,arrayshape)
c
c     ip13
c
      arrayshape2=(/maxp13,n/)
      if (hostrank == 0) then
        windowsize = int(maxp13*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ip13,arrayshape2)
c
c     np13
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,np13,arrayshape)
c
c     ip14
c
      arrayshape2=(/maxp14,n/)
      if (hostrank == 0) then
        windowsize = int(maxp14*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ip14,arrayshape2)
c
c     np14
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,np14,arrayshape)
      return
      end
