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
#include "tinker_precision.h"
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
      use utilgpu
      use tinMemory
      use mpi
      implicit none
      integer istep,modnl,ierr
      integer i,j,k,l,m,ia,iipole
      integer iproc,iglob,polecount
      integer npole_cap,ipr1,ibufbeg
      integer sizepole
      integer ji,ki,li
      integer it,jt,kt,lt
      integer ic,imp,nmp
      integer size_t,next
      integer number
      integer kz,kx,ky
      integer ztyp,xtyp,ytyp
      !integer:: n3_f=0,nbis=0,nzbi=0,nzon=0,nzth=0,none=0
      integer, allocatable :: mpt(:)
      integer, allocatable :: mpz(:)
      integer, allocatable :: mpx(:)
      integer, allocatable :: mpy(:)
      real(t_p) pel,pal
      real(t_p) mpl(13)
      real(t_p) d
      logical header,path
      character*4 pa,pb,pc,pd
      character*8 axt
      character*16 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
      logical init
!$acc routine(distprocpart1)
c
      if (init) then
c
c     allocate global arrays
c
        if (deb_Path) print*,'kmpole init'

        call prmem_request(rpole,13,n)
        call prmem_request(uind ,03,n)
        call prmem_request(uinp ,03,n)
        call prmem_request(polelocnl ,n)
        call prmem_request(poleloc   ,n)
        call prmem_request(polerecloc,n)

        uinp_p => uinp
        uind_p => uind

        vmxx = 0d0
        vmxy = 0d0
        vmxz = 0d0
        vmyy = 0d0
        vmyz = 0d0
        vmzz = 0d0
!$acc update host(vmxx,vmxy,vmxz,vmyy,vmyz,vmzz)

!$acc parallel loop collapse(2)
!$acc&         present(uind,uinp)
        do i = 1,n
           do j = 1,3
              uind(j,i) = 0
              uinp(j,i) = 0
              if(i.eq.1.and.j.eq.1) then
                 uind(1,1) = 1
                 uinp(1,1) = 1
              end if
           end do
        end do
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
                 mpl(j) = 0.0_ti_p
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
               size_t = 4
               call numeral (k,pa,size_t)
               call numeral (kz,pb,size_t)
               call numeral (kx,pc,size_t)
               call numeral (ky,pd,size_t)
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
              pole(j,i) = 0.0_ti_p
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
                 mpl(j) = 0.0_ti_p
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
              pole(k,i) = pole(k,i) * bohr**2 / 3.0_ti_p
           end do
        end do
c 230    call MPI_BARRIER(hostcomm,ierr)
c
c       get the order of the multipole expansion at each site
c
        npole = n
        do i = 1, n
           size_t = 0
           do k = 1, maxpole
              if (pole(k,i) .ne. 0.0_ti_p)  size_t = max(k,size_t)
           end do
           if (size_t .gt. 4) then
              size_t = 13
           else if (size_t .gt. 1) then
              size_t = 4
           end if
           polsiz(i) = size_t
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
              pel = 0.0_ti_p
              pal = 0.0_ti_p
              string = record(next:240)
              read (string,*,err=340,end=340)  k,pel,pal
              cpele(k) = abs(pel)
              cpalp(k) = pal
              if (header .and. .not.silent) then
                 header = .false.
                 write (iout,320)
  320           format (/,' Additional Charge Penetration Parameters :',
     &                  //,5x,'Atom Class',11x,'Core Chg',11x,'Damp',/)
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
           pcore(i)  = 0.0_ti_p
           pval(i)   = pole(1,i)
           pval0(i)  = pval(i)
           palpha(i) = 0.0_ti_p
           ic = class(i)
           if (ic .ne. 0) then
              pcore(i)  = cpele(ic)
              pval(i)   = pole(1,i) - cpele(ic)
              pval0(i)  = pval(i)
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
              pel = 0.0_ti_p
              pal = 0.0_ti_p
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
                 nbpole(i)     = npole
                 npole         = npole + 1
                 ipole(npole)  = i
                 pollist(i)    = npole
                 zaxis(npole)  = zaxis(i)
                 xaxis(npole)  = xaxis(i)
                 yaxis(npole)  = yaxis(i)
                 polaxe(npole) = polaxe(i)
                 do j = 1, maxpole
                    pole(j,npole) = pole(j,i)
                 end do
                 mono0(npole)  = pole(1,i)
                 if (palpha(i) .ne. 0.0)  ncp = ncp + 1
                 pcore(npole)  = pcore(i)
                 pval(npole)   = pval(i)
                 pval0(npole)  = pval(i)
                 palpha(npole) = palpha(i)
              end if
           end do
           call FromPolaxe2Ipolaxe ! Brodcast 'nZ_onlyglob' after this call
        end if
c
c       copy original multipole values that won't change during mutation
c
        if (use_lambdadyn) then
           pole_orig(1:13,1:n) = pole(1:13,1:n)
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

        if (.not. use_polar)
     &     call MPI_BCAST(nZ_Onlyglob,1,MPI_INT,0,hostcomm,ierr)
        !TODO 1.2 Why xaxis need to be allocatable
c
c       Update use_mpole switch
c
        if (npole.eq.0)  then
           use_mpole = .false.
           use_mlist = .false.
!$acc update device(use_mpole)
        end if

        call upload_device_mpole
c
c       if polarization not used, zero out induced dipoles
c
        if (.not. use_polar) then
!$acc parallel loop collapse(2) 
!$acc&         present(uind,uinp) async
           do i = 1, n
              do j = 1, 3
                 uind(j,i) = 0.0_ti_p
                 uinp(j,i) = 0.0_ti_p
              if(i.eq.1.and.j.eq.1) then
                 uind(1,1) = 1
                 uinp(1,1) = 1
              end if
              end do
           end do
        end if
      end if
c
      call prmem_request(poleglob   ,   nbloc,async=.true.)
      call prmem_request(polerecglob,nlocrec2,async=.true.)
c
c       remove any zero or undefined atomic multipoles
c
      !TODO Clean this part
      if (.not.use_mpole.or.use_polar) return

         if (nproc.eq.1) then
            call kpolar_reset1(istep)
         else
            call kpolar_reset(istep)
         end if

      end

      subroutine FromPolaxe2Ipolaxe
      use atoms  ,only: n
      use domdec ,only: rank
      use inform ,only: deb_Path
      use mpole  ,only: polaxe,ipolaxe,nZ_Onlyglob
     &           ,Ax_3_Fold,Ax_Bisector,Ax_Z_Bisect,Ax_Z_Only
     &           ,Ax_Z_Then_X,Ax_None
      use sizes  ,only: tinkerdebug
      implicit none
      integer i
      integer n3_f,nbis,nzbi,nzon,nzth,none

      !                  !!!   WARNING    !!!
      ! This routine is called normally called by first process of
      ! each comput node able to share memory
      ! If this is the case, do not forget to broadcast 'nZ_Onlyglob'
      ! after a process synchronisation call
c
c       Translate polaxe to integer to acc GPU comput
c
        nZ_Onlyglob = 0
        do i=1,n
           select case (polaxe(i))
              case ('3-Fold'  ) ; ipolaxe(i) = Ax_3_Fold
              case ('Bisector') ; ipolaxe(i) = Ax_Bisector
              case ('Z-Bisect') ; ipolaxe(i) = Ax_Z_Bisect
              case ('Z-Only'  ) ; 
                      ipolaxe(i) = Ax_Z_Only
                      nZ_Onlyglob = nZ_Onlyglob+1
              case ('Z-then-X') ; ipolaxe(i) = Ax_Z_Then_X
              case ('None'    ) ; ipolaxe(i) = Ax_None
              case default      ; ipolaxe(i) = Ax_None
           endselect
        end do
c       if (nZ_Onlyglob.ne.0.and.rank.eq.0)
c    &     print*,'nZ_Only axe : ',nZ_Onlyglob

        if (tinkerdebug.gt.0) then
           n3_f = 0
           nbis = 0
           nzbi = 0
           nzon = 0
           nzth = 0
           none = 0
           do i=1,n
              select case (polaxe(i))
                 case ('3-Fold'  );n3_f =n3_f +1
                 case ('Bisector');nbis =nbis +1  
                 case ('Z-Bisect');nzbi =nzbi +1  
                 case ('Z-Only'  );nzon =nzon +1  
                 case ('Z-then-X');nzth =nzth +1   
                 case ('None'    );none =none +1 
                 case default     ;none =none +1 
              end select
           end do
 13     format('   3-Fold   Bisect Z-Bisect   Z-Only Z-then-X     None')
 14     format(6I9)
           write(*,*) " >>>>>   type   polar axe      Num   >>>>>"
           write(*,13)
           write(*,14) n3_f,nbis,nzbi,nzon,nzth,none
           write(*,*) " <<<<< "
        end if
      end subroutine

      subroutine upload_device_mpole
      use chgpen
      use domdec ,only: rank,hostcomm
     &           ,bufbegpole,domlenpole,domlenpolerec
      use inform ,only: deb_Path
      use mpi    ,only: MPI_BARRIER
      use mpole
      use polar
      use potent ,only: use_lambdadyn
      use polgrp
      use sizes
      implicit none
      integer ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_mpole')
      if(deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif
!$acc update device(nbpole,polsiz
!$acc&       ,ipole,xaxis,yaxis,zaxis
!$acc&       ,ipolaxe,pollist,pole
!$acc&       ,mono0,pcore,pval,pval0,palpha)
      if (use_lambdadyn) then
!$acc  update device(pole_orig)
      end if

!$acc enter data copyin(npoleloc,npolebloc,npolerecloc
!$acc&      ,npolerecloc_old,npolelocnl,bufbegpole
!$acc&      ,domlenpole,domlenpolerec)

      end subroutine

      subroutine dealloc_shared_mpole
      use chgpen
      use domdec
      use inform ,only: deb_Path
      use mpole
      use polar
      use polgrp
      use potent
      use sizes
      use tinMemory
      implicit none

 12   format(2x,'dealloc_shared_mpole')
      if(deb_Path) print 12

      call shmem_request(polsiz, winpolsiz, [0],config=mhostacc)
      call shmem_request(zaxis,  winzaxis,  [0],config=mhostacc)
      call shmem_request(yaxis,  winyaxis,  [0],config=mhostacc)
      call shmem_request(xaxis,  winxaxis,  [0],config=mhostacc)
      call shmem_request(pole,   winpole,[13,0],config=mhostacc)
      if (use_lambdadyn) then
      call shmem_request(pole_orig,winpole_orig,[13,0],config=mhostacc)
      end if
      call shmem_request(mono0,   winmono0, [0],config=mhostacc)
      call shmem_request(polaxe, winpolaxe, [0])
      call shmem_request(ipolaxe,winipolaxe,[0],config=mhostacc)
      call shmem_request(ipole,  winipole,  [0],config=mhostacc)
      call shmem_request(pollist,winpollist,[0],config=mhostacc)
      call shmem_request(nbpole, winnbpole, [0],config=mhostacc)
      call shmem_request(pcore,   winpcore, [0],config=mhostacc)
      call shmem_request(pval,     winpval, [0],config=mhostacc)
      call shmem_request(pval0,   winpval0, [0],config=mhostacc)
      call shmem_request(palpha, winpalpha, [0],config=mhostacc)

      call shmem_request(ip11,winip11,[maxp11,0])
      call shmem_request(np11,winnp11,[0])
      call shmem_request(ip12,winip12,[maxp12,0])
      call shmem_request(np12,winnp12,[0])
      call shmem_request(ip13,winip13,[maxp13,0])
      call shmem_request(np13,winnp13,[0])
      call shmem_request(ip14,winip14,[maxp14,0])
      call shmem_request(np14,winnp14,[0])

      end subroutine
c
c     subroutine alloc_shared_mpole : allocate shared memory pointers for mpole
c     parameter arrays
c
      subroutine alloc_shared_mpole
      use sizes
      use atoms
      use chgpen
      use domdec
      use mpole
      use polgrp
      use potent
      use mpi
      use tinMemory
      implicit none
c
      if (associated(ipole).and.n.eq.size(pole)) return

      call shmem_request(polsiz, winpolsiz, [n],config=mhostacc)
      call shmem_request(zaxis,   winzaxis, [n],config=mhostacc)
      call shmem_request(yaxis,   winyaxis, [n],config=mhostacc)
      call shmem_request(xaxis,   winxaxis, [n],config=mhostacc)
      call shmem_request(pole,   winpole,[13,n],config=mhostacc)
      if (use_lambdadyn) then
      call shmem_request(pole_orig,winpole_orig,[13,n],config=mhostacc)
      end if
      call shmem_request(mono0,   winmono0, [n],config=mhostacc)
      call shmem_request(polaxe, winpolaxe, [n])
      call shmem_request(ipolaxe,winipolaxe,[n],config=mhostacc)
      call shmem_request(ipole,   winipole, [n],config=mhostacc)
      call shmem_request(pollist,winpollist,[n],config=mhostacc)
      call shmem_request(nbpole, winnbpole, [n],config=mhostacc)
      call shmem_request(pcore,   winpcore, [n],config=mhostacc)
      call shmem_request(pval,     winpval, [n],config=mhostacc)
      call shmem_request(pval0,   winpval0, [n],config=mhostacc)
      call shmem_request(palpha, winpalpha, [n],config=mhostacc)
c
      call shmem_request(ip11,winip11,[maxp11,n])
      call shmem_request(np11,winnp11,[n])
      call shmem_request(ip12,winip12,[maxp12,n])
      call shmem_request(np12,winnp12,[n])
      call shmem_request(ip13,winip13,[maxp13,n])
      call shmem_request(np13,winnp13,[n])
      call shmem_request(ip14,winip14,[maxp14,n])
      call shmem_request(np14,winnp14,[n])
      end
