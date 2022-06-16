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
      use sizes
      use atmlst
      use atmtyp
      use atoms
      use chgpen
      use couple
      use cutoff
      use domdec
      use inform
      use iounit
      use kcpen
      use keys
      use kmulti
      use mpole
      use neigh
      use pme
      use polar
      use polgrp
      use potent
      use units
      use mpi
      implicit none
      integer istep,modnl,ierr
      integer i,j,k,l,m,ia,iipole
      integer iproc,iglob,polecount
      integer ji,ki,li
      integer it,jt,kt,lt
      integer ic,imp,nmp
      integer size,next
      integer number
      integer kz,kx,ky
      integer ztyp,xtyp,ytyp
      integer, allocatable :: mpt(:)
      integer, allocatable :: mpz(:)
      integer, allocatable :: mpx(:)
      integer, allocatable :: mpy(:)
      real*8 pel,pal
      real*8 mpl(13)
      real*8 d
      logical header,path
      character*4 pa,pb,pc,pd
      character*8 axt
      character*16 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
      logical init
c
      if (init) then
c
c     allocate global arrays
c
        if (allocated(rpole)) deallocate (rpole)
        allocate (rpole(13,n))
        if (allocated(xaxis)) deallocate(xaxis)
        allocate (xaxis(n))
        if (allocated(yaxis)) deallocate(yaxis)
        allocate (yaxis(n))
        if (allocated(zaxis)) deallocate(zaxis)
        allocate (zaxis(n))
        if (allocated(polelocnl)) deallocate(polelocnl)
        allocate (polelocnl(n))
        if (allocated(poleloc)) deallocate(poleloc)
        allocate (poleloc(n))
        if (allocated(polerecloc)) deallocate(polerecloc)
        allocate (polerecloc(n))
        if (allocated(poleglob)) deallocate(poleglob)
        allocate (poleglob(n))
        if (allocated(polerecglob)) deallocate(polerecglob)
        allocate (polerecglob(n))
        if (allocated(uind)) deallocate (uind)
        if (allocated(uinp)) deallocate (uinp)
        allocate (uind(3,n))
        allocate (uinp(3,n))
        uind = 0d0
        uinp = 0d0

        vmxx = 0d0
        vmxy = 0d0
        vmxz = 0d0
        vmyy = 0d0
        vmyz = 0d0
        vmzz = 0d0
c
c       deallocate global pointers if necessary
c
        call dealloc_shared_mpole
c
c       allocate global pointers
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
              string = record(next:240)
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
              string = record(next:240)
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
c       zero out local axes, multipoles and polarization attachments
c
        npole = 0
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
        end do
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
              string = record(next:240)
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
c 230    call MPI_BARRIER(hostcomm,ierr)
c
c       get the order of the multipole expansion at each site
c
        npole = n
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
c     find new charge penetration parameters in the keyfile
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'CHGPEN ') then
              k = 0
              pel = 0.0d0
              pal = 0.0d0
              string = record(next:240)
              read (string,*,err=340,end=340)  k,pel,pal
              cpele(k) = abs(pel)
              cpalp(k) = pal
              if (header .and. .not.silent) then
                 header = .false.
                 write (iout,320)
  320            format (/,' Additional Charge Penetration Parameters :',
     &                   //,5x,'Atom Class',11x,'Core Chg',11x,'Damp',/)
              end if
              if (.not. silent) then
                 write (iout,330)  k,pel,pal
  330            format (6x,i6,7x,f15.3,f15.4)
              end if
  340         continue
           end if
        end do
c
c     assign the charge penetration charge and alpha parameters 
c     
        ncp = 0
        do i = 1, n
           pcore(i) = 0.0d0
           pval(i) = pole(1,i)
           pval0(i) = pval(i)
           palpha(i) = 0.0d0
           ic = class(i)
           if (ic .ne. 0) then
              pcore(i) = cpele(ic)
              pval(i) = pole(1,i) - cpele(ic)
              pval0(i) = pval(i)
              palpha(i) = cpalp(ic)
           end if
        end do
c
c     process keywords with charge penetration for specific atoms
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:7) .eq. 'CHGPEN ') then
              k = 0
              pel = 0.0d0
              pal = 0.0d0
              string = record(next:240)
              read (string,*,err=370,end=370)  k,pel,pal
              if (k.lt.0 .and. k.ge.-n) then
                 k = -k
                 pcore(k) = abs(pel)
                 pval(k) = pole(1,k) - abs(pel)
                 palpha(k) = pal
                 if (header .and. .not.silent) then
                    header = .false.
                    write (iout,350)
  350               format (/,' Additional Charge Penetration',
     &                         ' for Specific Atoms :',
     &                      //,5x,'Atom',17x,'Core Chg',11x,'Damp',/)
                 end if
                 if (.not. silent) then
                    write (iout,360)  k,pel,pal
  360               format (6x,i6,7x,f15.3,f15.4)
                 end if
              end if
  370         continue
           end if
        end do
c
c       remove any zero or undefined atomic multipoles
c
        if ((.not.use_polar).and.(.not.use_chgtrn)) then
           npole = 0
           ncp = 0
           do i = 1, n
              if (polsiz(i) .ne. 0) then
                 nbpole(i) = npole
                 npole = npole + 1
                 ipole(npole) = i
                 pollist(i) = npole
                 zaxis(npole) = zaxis(i)
                 xaxis(npole) = xaxis(i)
                 yaxis(npole) = yaxis(i)
                 polaxe(npole) = polaxe(i)
                 do j = 1, maxpole
                    pole(j,npole) = pole(j,i)
                 end do
                 mono0(npole) = pole(1,i)
                 if (palpha(i) .ne. 0.0d0)  ncp = ncp + 1
                 pcore(npole) = pcore(i)
                 pval(npole) = pval(i)
                 pval0(npole) = pval(i)
                 palpha(npole) = palpha(i)
              end if
           end do
         end if
c
c       test multipoles at chiral sites and invert if necessary
c
        if (use_mpole .and. .not.use_polar .and. .not.use_chgtrn)
     &   call chkpole(.true.)
c
 230    call MPI_BARRIER(hostcomm,ierr)
        call MPI_BCAST(npole,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(ncp,1,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(xaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(yaxis,n,MPI_INT,0,hostcomm,ierr)
        call MPI_BCAST(zaxis,n,MPI_INT,0,hostcomm,ierr)
c
        if (npole .eq.0)  then
          use_mpole = .false.
        end if
        if (ncp .ne. 0)  use_chgpen = .true.
c
c       if polarization not used, zero out induced dipoles
c
        if (.not. use_polar) then
           do i = 1, n
              do j = 1, 3
                 uind(j,i) = 0.0d0
                 uinp(j,i) = 0.0d0
              end do
           end do
        end if
c
c       copy original multipole values that won't change during mutation
c
        pole_orig = pole
      end if
c
c
c       remove any zero or undefined atomic multipoles
c
      if (.not.use_polar) then
         npoleloc = 0
         do i = 1, nloc
            iglob = glob(i)
            iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
            if (iipole.eq.0) cycle
            polecount = nbpole(iglob)
            if (polsiz(iglob) .ne. 0) then
               npoleloc = npoleloc + 1
               poleglob(npoleloc) = polecount + 1
               poleloc(polecount+1) = npoleloc
            end if
         end do
        domlenpole(rank+1) = npoleloc
        npolebloc = npoleloc
        do iproc = 1, n_recep1
          if (domlen(p_recep1(iproc)+1).ne.0) then
            bufbegpole(p_recep1(iproc)+1) = npolebloc + 1
          else
            bufbegpole(p_recep1(iproc)+1) = 1
          end if
          do i = 1, domlen(p_recep1(iproc)+1)
            iglob = glob(bufbeg(p_recep1(iproc)+1)+i-1)
            iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
            if (iipole.eq.0) cycle
            polecount = nbpole(iglob)
            if (polsiz(iglob) .ne. 0.0d0) then
              npolebloc = npolebloc + 1
              poleglob(npolebloc) = polecount + 1
              poleloc(polecount+1) = npolebloc
            end if
          end do
          if (domlen(p_recep1(iproc)+1).ne.0) then
            domlenpole(p_recep1(iproc)+1) = 
     $        npolebloc-bufbegpole(p_recep1(iproc)+1)+1
          else
            domlenpole(p_recep1(iproc)+1) = 0
          end if
        end do
c
        npolerecloc = 0
        do i = 1, nlocrec
           iglob = globrec(i)
           iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
           if (iipole.eq.0) cycle
           polecount = nbpole(iglob)
           if (polsiz(iglob) .ne. 0) then
              npolerecloc = npolerecloc + 1
              polerecglob(npolerecloc) = polecount + 1
              polerecloc(polecount+1) = npolerecloc
           end if
        end do
        domlenpolerec(rank+1) = npolerecloc
c
        do i = nlocrec+1, nlocrec2
           iglob = globrec(i)
           iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
           if (iipole.eq.0) cycle
           polecount = nbpole(iglob)
           if (polsiz(iglob) .ne. 0 ) then
              npolerecloc = npolerecloc + 1
              polerecglob(npolerecloc) = polecount + 1
              polerecloc(polecount+1) = npolerecloc
           end if
        end do
        npolerecloc = domlenpolerec(rank+1)
c
        modnl = mod(istep,ineigup)
        if (istep.eq.-1) return
        if (modnl.ne.0) return
        if (allocated(poleglobnl)) deallocate(poleglobnl)
        allocate (poleglobnl(nlocnl))
c
        npolelocnl = 0
        do i = 1, nlocnl
          iglob = ineignl(i)
          iipole = pollist(iglob)
c
c   skip atom if it is not in the multipole list
c
          if (iipole.eq.0) cycle
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
c     subroutine dealloc_shared_mpole : deallocate shared memory pointers for mpole
c     parameter arrays
c
      subroutine dealloc_shared_mpole
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use chgpen
      use mpole
      use polgrp
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
c
      if (associated(polsiz)) then
        CALL MPI_Win_shared_query(winpolsiz, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpolsiz,ierr)
      end if
      if (associated(pole)) then
        CALL MPI_Win_shared_query(winpole, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpole,ierr)
      end if
      if (associated(pole_orig)) then
        CALL MPI_Win_shared_query(winpole_orig,0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpole_orig,ierr)
      end if
      if (associated(mono0)) then
        CALL MPI_Win_shared_query(winmono0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winmono0,ierr)
      end if
      if (associated(polaxe)) then
        CALL MPI_Win_shared_query(winpolaxe, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpolaxe,ierr)
      end if
      if (associated(ipole)) then
        CALL MPI_Win_shared_query(winipole, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winipole,ierr)
      end if
      if (associated(pollist)) then
        CALL MPI_Win_shared_query(winpollist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpollist,ierr)
      end if
      if (associated(nbpole)) then
        CALL MPI_Win_shared_query(winnbpole, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbpole,ierr)
      end if
      if (associated(pcore)) then
        CALL MPI_Win_shared_query(winpcore, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpcore,ierr)
      end if
      if (associated(pval)) then
        CALL MPI_Win_shared_query(winpval, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpval,ierr)
      end if
      if (associated(pval0)) then
        CALL MPI_Win_shared_query(winpval0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpval0,ierr)
      end if
      if (associated(palpha)) then
        CALL MPI_Win_shared_query(winpalpha, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpalpha,ierr)
      end if
      if (associated(ip11)) then
        CALL MPI_Win_shared_query(winip11, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winip11,ierr)
      end if
      if (associated(np11)) then
        CALL MPI_Win_shared_query(winnp11, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnp11,ierr)
      end if
      if (associated(ip12)) then
        CALL MPI_Win_shared_query(winip12, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winip12,ierr)
      end if
      if (associated(np12)) then
        CALL MPI_Win_shared_query(winnp12, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnp12,ierr)
      end if
      if (associated(ip13)) then
        CALL MPI_Win_shared_query(winip13, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winip13,ierr)
      end if
      if (associated(np13)) then
        CALL MPI_Win_shared_query(winnp13, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnp13,ierr)
      end if
      if (associated(ip14)) then
        CALL MPI_Win_shared_query(winip14, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winip14,ierr)
      end if
      if (associated(np14)) then
        CALL MPI_Win_shared_query(winnp14, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnp14,ierr)
      end if
      return
      end
c
c     subroutine alloc_shared_mpole : allocate shared memory pointers for mpole
c     parameter arrays
c
      subroutine alloc_shared_mpole
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use chgpen
      use domdec
      use mpole
      use polgrp
      use potent
      use mpi
      implicit none
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
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
     $  hostcomm, baseptr, winpolsiz, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpolsiz, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,polsiz,arrayshape)
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
     $  hostcomm, baseptr, winpole, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpole, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pole,arrayshape2)
c
c     pole_orig
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
     $  hostcomm, baseptr, winpole_orig, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpole_orig,0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pole_orig,arrayshape2)
c
c     mono0
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
     $  hostcomm, baseptr, winmono0, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winmono0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,mono0,arrayshape)
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
     $  hostcomm, baseptr, winpolaxe, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpolaxe, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winipole, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winipole, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winpollist, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpollist, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winnbpole, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbpole, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbpole,arrayshape)
c
c     pcore
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
     $  hostcomm, baseptr, winpcore, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpcore, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pcore,arrayshape)
c
c     pval
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
     $  hostcomm, baseptr, winpval, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpval, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pval,arrayshape)
c
c     pval0
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
     $  hostcomm, baseptr, winpval0, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpval0, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pval0,arrayshape)
c
c     palpha
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
     $  hostcomm, baseptr, winpalpha, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpalpha, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,palpha,arrayshape)
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
     $  hostcomm, baseptr, winip11, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winip11, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winnp11, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnp11, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winip12, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winip12, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winnp12, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnp12, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winip13, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winip13, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winnp13, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnp13, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winip14, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winip14, 0, windowsize, disp_unit,
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
     $  hostcomm, baseptr, winnp14, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnp14, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,np14,arrayshape)
      return
      end
