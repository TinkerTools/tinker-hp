c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kvdw  --  van der Waals parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kvdw" assigns the parameters to be used in computing the
c     van der Waals interactions and processes any new or changed
c     values for these parameters
c
c
#include "tinker_precision.h"
      subroutine kvdw(init,istep)
      use atmlst
      use atmtyp
      use atoms
      use couple
      use cutoff
      use domdec
      use fields
      use keys
      use inform
      use iounit
      use khbond
      use kvdwpr
      use kvdws
      use math
      use merck
      use neigh
      use potent
      use tinheader ,only:ti_p,re_p
      use tinMemory ,only:prmem_request
      use vdw
      use vdwpot
      implicit none
      integer istep,modnl,ierr
      integer i,k,it
      integer ia,ib,next
      integer size,number
      integer iglob,vdwcount,iproc
      integer nvdw_cap,ipr2,ibufbeg
      real(t_p) zero,one,two,four
      real(t_p) rd,ep,rdn,gik
      real(t_p), allocatable :: srad(:)
      real(t_p), allocatable :: srad4(:)
      real(t_p), allocatable :: seps(:)
      real(t_p), allocatable :: seps4(:)
      logical header
      character*4 pa,pb
      character*8 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
      logical init
      parameter(zero=0.0_ti_p,  one=1.0_ti_p,
     &           two=2.0_ti_p, four=4.0_ti_p)
c
      blank = '        '
      if (init) then
cc
cc     allocate global arrays
cc
c        call alloc_shared_vdw
c        if (hostrank.ne.0) goto 1000
c
c       process keywords containing van der Waals parameters
c
        if (deb_Path) print*,'kvdw init'
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:4) .eq. 'VDW ') then
              call getnumb (record,k,next)
              if (k.ge.1 .and. k.le.maxclass) then
                 rd = rad(k)
                 ep = eps(k)
                 rdn = reduct(k)
                 string = record(next:120)
                 read (string,*,err=10,end=10)  rd,ep,rdn
   10            continue
                 if (header .and. .not.silent) then
                    header = .false.
                    if (vdwindex .eq. 'CLASS') then
                     if (rank.eq.0) write (iout,20)
   20                format (/,' Additional van der Waals Parameters :',
     &                      //,5x,'Atom Class',10x,'Size',6x,
     &                            'Epsilon',5x,'Reduction',/)
                    else
                     if (rank.eq.0) write (iout,30)
   30                format (/,' Additional van der Waals Parameters :',
     &                         //,5x,'Atom Type',11x,'Size',6x,
     &                            'Epsilon',5x,'Reduction',/)
                    end if
                 end if
                 rad(k) = rd
                 eps(k) = ep
                 reduct(k) = rdn
                 if (.not. silent) then
                    if (rank.eq.0) write (iout,40)  k,rd,ep,rdn
   40               format (4x,i6,8x,2f12.4,f12.3)
                 end if
              else if (k .gt. maxclass) then
                 if (rank.eq.0) write (iout,50)  maxclass
   50            format (/,' KVDW  --  Only Atom Classes through',i4,
     &                      ' are Allowed')
                 abort = .true.
              end if
           end if
        end do
c
c       process keywords containing 1-4 van der Waals parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:6) .eq. 'VDW14 ') then
              call getnumb (record,k,next)
              if (k.ge.1 .and. k.le.maxclass) then
                 rd = rad4(k)
                 ep = eps4(k)
                 string = record(next:120)
                 read (string,*,err=60,end=60)  rd,ep
   60            continue
                 if (header .and. .not.silent) then
                    header = .false.
                    if (vdwindex .eq. 'CLASS') then
                       if (rank.eq.0) write (iout,70)
   70                  format (/,' Additional 1-4 van der Waals',
     &                            ' Parameters :',
     &                         //,5x,'Atom Class',10x,'Size',6x,
     &                            'Epsilon',/)
                    else
                       if (rank.eq.0) write (iout,80)
   80                  format (/,' Additional 1-4 van der Waals',
     &                            ' Parameters :',
     &                         //,5x,'Atom Type',11x,'Size',6x,
     &                            'Epsilon',/)
                    end if
                 end if
                 rad4(k) = rd
                 eps4(k) = ep
                 if (.not. silent) then
                    if (rank.eq.0) write (iout,90)  k,rd,ep
   90               format (4x,i6,8x,2f12.4)
                 end if
              else if (k .gt. maxclass) then
                 if (rank.eq.0) write (iout,100)  maxclass
  100            format (/,' KVDW  --  Only Atom Classes through',i4,
     &                      ' are Allowed')
                 abort = .true.
              end if
           end if
        end do
c
c       process keywords containing specific pair vdw parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:6) .eq. 'VDWPR ') then
              ia = 0
              ib = 0
              rd = zero
              ep = zero
              string = record(next:120)
              read (string,*,err=150,end=150)  ia,ib,rd,ep
              if (header .and. .not.silent) then
                 header = .false.
                 if (vdwindex .eq. 'CLASS') then
                    if (rank.eq.0) write (iout,110)
  110               format (/,' Additional van der Waals Parameters',
     &                         ' for Specific Pairs :',
     &                      //,5x,'Atom Classes',6x,'Size Sum',
     &                         4x,'Epsilon',/)
                 else
                    if (rank.eq.0) write (iout,120)
  120               format (/,' Additional van der Waals Parameters',
     &                         ' for Specific Pairs :',
     &                      //,5x,'Atom Types',8x,'Size Sum',
     &                         4x,'Epsilon',/)
                 end if
              end if
              if (.not. silent) then
                 if (rank.eq.0) write (iout,130)  ia,ib,rd,ep
  130            format (6x,2i4,4x,2f12.4)
              end if
              size = 4
              call numeral (ia,pa,size)
              call numeral (ib,pb,size)
              if (ia .le. ib) then
                 pt = pa//pb
              else
                 pt = pb//pa
              end if
              do k = 1, maxnvp
                 if (kvpr(k).eq.blank .or. kvpr(k).eq.pt) then
                    kvpr(k) = pt
                    radpr(k) = rd
                    epspr(k) = ep
                    goto 150
                 end if
              end do
              if (rank.eq.0) write (iout,140)
  140         format (/,' KVDW  --  Too many Special VDW Pair',
     &                   ' Parameters')
              abort = .true.
  150         continue
           end if
        end do
c
c       process keywords containing hydrogen bonding vdw parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:6) .eq. 'HBOND ') then
              ia = 0
              ib = 0
              rd = zero
              ep = zero
              string = record(next:120)
              read (string,*,err=200,end=200)  ia,ib,rd,ep
              if (header .and. .not.silent) then
                 header = .false.
                 if (vdwindex .eq. 'CLASS') then
                    if (rank.eq.0) write (iout,160)
  160               format (/,' Additional van der Waals Hydrogen',
     &                         ' Bonding Parameters :',
     &                      //,5x,'Atom Classes',6x,'Size Sum',
     &                         4x,'Epsilon',/)
                 else
                    if (rank.eq.0) write (iout,170)
  170               format (/,' Additional van der Waals Hydrogen',
     &                         ' Bonding Parameters :',
     &                      //,5x,'Atom Types',8x,'Size Sum',
     &                         4x,'Epsilon',/)
                 end if
              end if
              if (.not. silent) then
                 if (rank.eq.0) write (iout,180)  ia,ib,rd,ep
  180            format (6x,2i4,4x,2f12.4)
              end if
              size = 4
              call numeral (ia,pa,size)
              call numeral (ib,pb,size)
              if (ia .le. ib) then
                 pt = pa//pb
              else
                 pt = pb//pa
              end if
              do k = 1, maxnvp
                 if (khb(k).eq.blank .or. khb(k).eq.pt) then
                    khb(k) = pt
                    radhb(k) = rd
                    epshb(k) = ep
                    goto 200
                 end if
              end do
              if (rank.eq.0) write (iout,190)
  190         format (/,' KVDW  --  Too many Hydrogen Bonding Pair',
     &                   ' Parameters')
              abort = .true.
  200         continue
           end if
        end do
c
c     allocate global arrays
c
        call alloc_shared_vdw
        if (hostrank.ne.0) goto 1000
c
c       use atom class or type as index into vdw parameters
c
        k = 0
        do i = 1, n
           jvdw(i) = class(i)
           if (vdwindex .eq. 'TYPE')  jvdw(i) = type(i)
           k = max(k,jvdw(i))
        end do
        if (k .gt. maxclass) then
           if (rank.eq.0) write (iout,210)
  210      format (/,' KVDW  --  Unable to Index VDW Parameters;',
     &                ' Increase MAXCLASS')
           abort = .true.
        end if
c
c       count the number of vdw types and their frequencies
c
        nvt = 0
        do i = 1, n
           it = jvdw(i)
           do k = 1, nvt
              if (ivt(k) .eq. it) then
                 jvt(k) = jvt(k) + 1
                 goto 220
              end if
           end do
           nvt = nvt + 1
           ivt(nvt) = it
           jvt(nvt) = 1
  220      continue
        end do
 1000   call MPI_BARRIER(hostcomm,ierr)
c
c       perform dynamic allocation of some local arrays
c
        allocate (srad(maxtyp))
        allocate (srad4(maxtyp))
        allocate (seps(maxtyp))
        allocate (seps4(maxtyp))
c
c       get the vdw radii and well depths for each atom type
c
        do i = 1, maxtyp
           if (rad4(i) .eq. zero)  rad4(i) = rad(i)
           if (eps4(i) .eq. zero)  eps4(i) = eps(i)
           if (radtyp .eq. 'SIGMA') then
              rad(i) = twosix * rad(i)
              rad4(i) = twosix * rad4(i)
           end if
           if (radsiz .eq. 'DIAMETER') then
              rad(i) = 0.5_ti_p * rad(i)
              rad4(i) = 0.5_ti_p * rad4(i)
           end if
           srad(i) = sqrt(rad(i))
           eps(i) = abs(eps(i))
           seps(i) = sqrt(eps(i))
           srad4(i) = sqrt(rad4(i))
           eps4(i) = abs(eps4(i))
           seps4(i) = sqrt(eps4(i))
        end do
c
c       use combination rules to set pairwise vdw radii sums
c
        do i = 1, maxclass
           do k = i, maxclass
              if (radrule(1:6) .eq. 'MMFF94') then
                 if (i .ne. k) then
                    if (DA(i).eq.'D' .or. DA(k).eq.'D') then
                       rd = 0.5_ti_p * (rad(i)+rad(k))
                    else
                       gik = (rad(i)-rad(k))/(rad(i)+rad(k))
                       rd = 0.5_ti_p * (rad(i)+rad(k))
     &                     *(one+0.2_ti_p*(one-exp(-12.0_ti_p*gik*gik)))
                    end if
                 else
                    rd = rad(i)
                 end if
              else if (rad(i).eq.zero .and. rad(k).eq.zero) then
                 rd = zero
              else if (radrule(1:10) .eq. 'ARITHMETIC') then
                 rd = rad(i) + rad(k)
              else if (radrule(1:9) .eq. 'GEOMETRIC') then
                 rd = two * (srad(i) * srad(k))
              else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
                 rd = two*(rad(i)**3+rad(k)**3)/(rad(i)**2+rad(k)**2)
              else
                 rd = rad(i) + rad(k)
              end if
              radmin(i,k) = rd
              radmin(k,i) = rd
           end do
        end do
c
c       use combination rules to set pairwise well depths
c
        do i = 1, maxclass
           do k = i, maxclass
              if (epsrule(1:6) .eq. 'MMFF94') then
                 ep = 181.16_ti_p*G(i)*G(k)*alph(i)*alph(k)
     &                   / ((sqrt(alph(i)/Nn(i))+sqrt(alph(k)/Nn(k)))
     &                                *radmin(i,k)**6)
                 if (i .eq. k)  eps(i) = ep
              else if (eps(i).eq.zero .and. eps(k).eq.zero) then
                 ep = zero
              else if (epsrule(1:10) .eq. 'ARITHMETIC') then
                 ep = 0.5_ti_p * (eps(i) + eps(k))
              else if (epsrule(1:9) .eq. 'GEOMETRIC') then
                 ep = seps(i) * seps(k)
              else if (epsrule(1:8) .eq. 'HARMONIC') then
                 ep = two * (eps(i)*eps(k)) / (eps(i)+eps(k))
              else if (epsrule(1:3) .eq. 'HHG') then
                 ep = four * (eps(i)*eps(k)) / (seps(i)+seps(k))**2
              else
                 ep = seps(i) * seps(k)
              end if
              epsilon(i,k) = ep
              epsilon(k,i) = ep
           end do
        end do
c
c       use combination rules to set pairwise 1-4 vdw radii sums
c
        do i = 1, maxclass
           do k = i, maxclass
              if (radrule(1:6) .eq. 'MMFF94') then
                 if (i .ne. k) then
                    if (DA(i).eq.'D' .or. DA(k).eq.'D') then
                       rd = 0.5_ti_p * (rad(i)+rad(k))
                    else
                       gik = (rad(i)-rad(k))/(rad(i)+rad(k))
                       rd = 0.5_ti_p * (rad(i)+rad(k))
     &                     *(one+0.2_ti_p*(one-exp(-12.0_ti_p*gik*gik)))
                    end if
                 else
                    rd = rad(i)
                 end if
              else if (rad4(i).eq.zero .and. rad4(k).eq.zero) then
                 rd = zero
              else if (radrule(1:10) .eq. 'ARITHMETIC') then
                 rd = rad4(i) + rad4(k)
              else if (radrule(1:9) .eq. 'GEOMETRIC') then
                 rd = two * (srad4(i) * srad4(k))
              else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
                 rd = two * (rad4(i)**3+rad4(k)**3)
     &                           / (rad4(i)**2+rad4(k)**2)
              else
                 rd = rad4(i) + rad4(k)
              end if
              radmin4(i,k) = rd
              radmin4(k,i) = rd
           end do
        end do
c
c       use combination rules to set pairwise 1-4 well depths
c
        do i = 1, maxclass
           do k = i, maxclass
              if (epsrule(1:6) .eq. 'MMFF94') then
                 ep = 181.16_ti_p*G(i)*G(k)*alph(i)*alph(k)
     &                   / ((sqrt(alph(i)/Nn(i))+sqrt(alph(k)/Nn(k)))
     &                                *radmin(i,k)**6)
                 if (i .eq. k)  eps4(i) = ep
              else if (eps4(i).eq.zero .and. eps4(k).eq.zero) then
                 ep = zero
              else if (epsrule(1:10) .eq. 'ARITHMETIC') then
                 ep = 0.5_ti_p * (eps4(i) + eps4(k))
              else if (epsrule(1:9) .eq. 'GEOMETRIC') then
                 ep = seps4(i) * seps4(k)
              else if (epsrule(1:8) .eq. 'HARMONIC') then
                 ep = two * (eps4(i)*eps4(k)) / (eps4(i)+eps4(k))
              else if (epsrule(1:3) .eq. 'HHG') then
                 ep = four * (eps4(i)*eps4(k)) / (seps4(i)+seps4(k))**2
              else
                 ep = seps4(i) * seps4(k)
              end if
              epsilon4(i,k) = ep
              epsilon4(k,i) = ep
           end do
        end do
c
c       perform deallocation of some local arrays
c
        deallocate (srad)
        deallocate (srad4)
        deallocate (seps)
        deallocate (seps4)
c
c       use reduced values for MMFF donor-acceptor pairs
c
        if (forcefield .eq. 'MMFF94') then
           do i = 1, maxclass
              do k = i, maxclass
                 if ((da(i).eq.'D' .and. da(k).eq.'A') .or.
     &               (da(i).eq.'A' .and. da(k).eq.'D')) then
                    epsilon(i,k) = epsilon(i,k) * 0.5_ti_p
                    epsilon(k,i) = epsilon(k,i) * 0.5_ti_p
                    radmin(i,k) = radmin(i,k) * 0.8_ti_p
                    radmin(k,i) = radmin(k,i) * 0.8_ti_p
                    epsilon4(i,k) = epsilon4(i,k) * 0.5_ti_p
                    epsilon4(k,i) = epsilon4(k,i) * 0.5_ti_p
                    radmin4(i,k) = radmin4(i,k) * 0.8_ti_p
                    radmin4(k,i) = radmin4(k,i) * 0.8_ti_p
                 end if
              end do
           end do
        end if
c
c       vdw reduction factor information for each individual atom
c
        do i = 1, n
           kred(i) = reduct(jvdw(i))
           if (n12(i).ne.1 .or. kred(i).eq.zero) then
              ired(i) = i
           else
              ired(i) = i12(1,i)
           end if
        end do
c
c       radii and well depths for special atom class pairs
c
        do i = 1, maxnvp
           if (kvpr(i) .eq. blank)  goto 230
           ia = number(kvpr(i)(1:4))
           ib = number(kvpr(i)(5:8))
           if (rad(ia) .eq. zero)  rad(ia) = 0.001_ti_p
           if (rad(ib) .eq. zero)  rad(ib) = 0.001_ti_p
           if (radtyp .eq. 'SIGMA')  radpr(i) = twosix * radpr(i)
           radmin(ia,ib) = radpr(i)
           radmin(ib,ia) = radpr(i)
           epsilon(ia,ib) = abs(epspr(i))
           epsilon(ib,ia) = abs(epspr(i))
           radmin4(ia,ib) = radpr(i)
           radmin4(ib,ia) = radpr(i)
           epsilon4(ia,ib) = abs(epspr(i))
           epsilon4(ib,ia) = abs(epspr(i))
        end do
  230   continue
c
c       radii and well depths for hydrogen bonding pairs
c
        if (vdwtyp .eq. 'MM3-HBOND') then
           do i = 1, maxclass
              do k = 1, maxclass
                 radhbnd(k,i) = zero
                 epshbnd(k,i) = zero
              end do
           end do
           do i = 1, maxnhb
              if (khb(i) .eq. blank)  goto 240
              ia = number(khb(i)(1:4))
              ib = number(khb(i)(5:8))
              if (rad(ia) .eq. zero)  rad(ia) = 0.001_ti_p
              if (rad(ib) .eq. zero)  rad(ib) = 0.001_ti_p
              if (radtyp .eq. 'SIGMA')  radhb(i) = twosix * radhb(i)
              radhbnd(ia,ib) = radhb(i)
              radhbnd(ib,ia) = radhb(i)
              epshbnd(ia,ib) = abs(epshb(i))
              epshbnd(ib,ia) = abs(epshb(i))
           end do
  240      continue
        end if
c
c       set coefficients for Gaussian fit to eps=1 and radmin=1
c
        if (vdwtyp .eq. 'GAUSSIAN') then
           if (gausstyp .eq. 'LJ-4') then
              ngauss = 4
              igauss(1,1) = 846706.7_ti_p
              igauss(2,1) = 15.464405_ti_p * twosix**2
              igauss(1,2) = 2713.651_ti_p
              igauss(2,2) = 7.346875_ti_p * twosix**2
              igauss(1,3) = -9.699172_ti_p
              igauss(2,3) = 1.8503725_ti_p * twosix**2
              igauss(1,4) = -0.7154420_ti_p
              igauss(2,4) = 0.639621_ti_p * twosix**2
           else if (gausstyp .eq. 'LJ-2') then
              ngauss = 2
              igauss(1,1) = 14487.1_ti_p
              igauss(2,1) = 9.05148_ti_p * twosix**2
              igauss(1,2) = -5.55338_ti_p
              igauss(2,2) = 1.22536_ti_p * twosix**2
           else if (gausstyp .eq. 'MM3-2') then
              ngauss = 2
              igauss(1,1) = 2438.886_ti_p
              igauss(2,1) = 9.342616_ti_p
              igauss(1,2) = -6.197368_ti_p
              igauss(2,2) = 1.564486_ti_p
           else if (gausstyp .eq. 'MM2-2') then
              ngauss = 2
              igauss(1,1) = 3423.562_ti_p
              igauss(2,1) = 9.692821_ti_p
              igauss(1,2) = -6.503760_ti_p
              igauss(2,2) = 1.585344_ti_p
           else if (gausstyp .eq. 'IN-PLACE') then
              ngauss = 2
              igauss(1,1) = 500.0_ti_p
              igauss(2,1) = 6.143_ti_p
              igauss(1,2) = -18.831_ti_p
              igauss(2,2) = 2.209_ti_p
           end if
        end if
c
c     remove zero-sized atoms from the list of local vdw sites
c
        nvdw = 0
        do i = 1, n
           if (rad(jvdw(i)) .ne. zero) then
              nbvdw(i) = nvdw
              nvdw = nvdw + 1
              ivdw(nvdw) = i
           end if
        end do
c
c       turn off the van der Waals potential if it is not used
c
        if (nvdw .eq. 0)  then
          use_vdw = .false.
          use_vlist = .false.
        end if

        if (use_vdw) then
           call upload_device_kvdw
           call prmem_request(vdwlocnl,nvdw)
           call compress_vdwData
        else
           call delete_data_vdw
           return
        end if

      end if
c
      if (nproc.eq.1) then
         call kvdw_reset1(istep)
      else
         call kvdw_reset(istep)
      end if

      end

      subroutine compress_vdwData
      use atoms  ,only: n
      use domdec ,only: hostrank,hostcomm
      use mpi
      use tinMemory ,only:shmem_request,mhostacc
      use sizes     ,only: maxclass
      use vdw
      implicit none
      integer i,j,key,class_i,class_j,posi,posj
      integer ierr
      integer sort_class(n),jvdw_key(n)
      integer scan_class(maxclass+1)
      real(t_p) rad_val

      call shmem_request(jvdw_c,winjvdw_c,[n],config=mhostacc)
      if (hostrank.ne.0) goto 30

      scan_class(:) = 0
      do i = 1,n
         sort_class(i) = jvdw(i)
      end do
      call sort3 (n,sort_class,jvdw_key)  !kort by key
c
c     Count number of system class
c     and renumber them from one
c
      nvdwclass = 1
      jvdw_c(jvdw_key(1)) = nvdwclass
      scan_class(nvdwclass+1) = 1

      do i = 2,n
         key = jvdw_key(i)
         if (sort_class(i).ne.sort_class(i-1)) nvdwclass=nvdwclass+1
         jvdw_c(key) = nvdwclass
         scan_class(nvdwclass+1) = scan_class(nvdwclass+1) + 1
      end do

      ! Scan class to pin point class changing in sorted array
      scan_class(1) = 1
      do i = 2,nvdwclass+1
         scan_class(i) = scan_class(i-1) + scan_class(i)
      end do

 30   continue
      call MPI_BCAST(nvdwclass,1,MPI_INT,0,hostcomm,ierr)
      ! Allocate memory according to number of classes
      call shmem_request(radmin_c  ,winradmin_c  ,[nvdwclass**2],
     &     config=mhostacc)
      call shmem_request(radmin4_c ,winradmin4_c ,[nvdwclass**2],
     &     config=mhostacc)
      call shmem_request(epsilon_c ,winepsilon_c ,[nvdwclass**2],
     &     config=mhostacc)
      call shmem_request(epsilon4_c,winepsilon4_c,[nvdwclass**2],
     &     config=mhostacc)

c
c     Compress radmin,etc,...
c     With renumbered classes
c
      if (hostrank.eq.0) then
         do i = 1,nvdwclass
            class_i = sort_class(scan_class(i))
            do j = 1,nvdwclass
               class_j = sort_class(scan_class(j))
               rad_val = radmin(class_j,class_i) 
               radmin_c((j-1)*nvdwclass+i) = rad_val
               radmin_c((i-1)*nvdwclass+j) = rad_val
               rad_val = radmin4(class_j,class_i) 
               radmin4_c((j-1)*nvdwclass+i) = rad_val
               radmin4_c((i-1)*nvdwclass+j) = rad_val
               rad_val = epsilon(class_j,class_i)
               epsilon_c((j-1)*nvdwclass+i) = rad_val
               epsilon_c((i-1)*nvdwclass+j) = rad_val
               rad_val = epsilon4(class_j,class_i)
               epsilon4_c((j-1)*nvdwclass+i) = rad_val
               epsilon4_c((i-1)*nvdwclass+j) = rad_val
            end do
         end do
      end if

      call mpi_barrier(hostcomm,ierr)
!$acc update device(jvdw_c,radmin_c,radmin4_c,epsilon_c,epsilon4_c)
      end subroutine

      subroutine kvdw_reset(istep)
      use atoms  ,only : x,y,z
      use atmlst ,only : vdwglob,vdwglobnl
      use domdec ,only : glob,p_recep2,domlen,nlocnl,repart,rank,
     &                   nbloc,nloc,bufbeg,n_recep2
      use kvdws  ,only : rad
      use neigh  ,only : ineigup,vbuf2,ineignl
      use utilgpu,only : prmem_request,rec_queue
#ifdef _OPENACC
     &                  ,rec_stream
      use thrust
#endif
      use vdw    ,only : nbvdw,jvdw,nvdwbloc,nvdwloc,nvdwlocnl,
     &                   vdwlocnl
      implicit none
      integer,intent(in) :: istep
      integer i,iglob,iproc,modnl
      integer vdwcount,nvdw_cap,ipr2,ibufbeg
      real(t_p) d
!$acc routine(distprocpart1)

      call prmem_request(vdwglob,nbloc,async=.true.)
!$acc data present(vdwglob,glob,nbvdw,rad,jvdw)
!$acc&     present(nvdwbloc,nvdwloc)
!$acc&     async(rec_queue)
c
c     remove zero-sized atoms from the list of vdw sites
c
!$acc serial async(rec_queue)
      nvdwloc = 0
!$acc end serial
c
!$acc parallel loop async(rec_queue)
      do i = 1, nloc
         iglob = glob(i)
         vdwcount = nbvdw(iglob)
         if (rad(jvdw(iglob)) .ne. 0) then
!$acc atomic capture
            nvdwloc = nvdwloc + 1
            nvdw_cap = nvdwloc
!$acc end atomic
            vdwglob(nvdw_cap) = vdwcount + 1
         end if
      end do
c
!$acc serial async(rec_queue)
      nvdwbloc = nvdwloc
!$acc end serial
c
      do iproc = 1, n_recep2
        ipr2    = p_recep2(iproc)+1 
        ibufbeg = bufbeg(ipr2)
!$acc parallel loop async(rec_queue)
        do i = 1, domlen(ipr2)
          iglob = glob(ibufbeg+i-1)
          vdwcount = nbvdw(iglob)
          if (rad(jvdw(iglob)) .ne. 0) then
!$acc atomic capture
            nvdwbloc = nvdwbloc + 1
            nvdw_cap = nvdwbloc
!$acc end atomic
            vdwglob(nvdw_cap) = vdwcount + 1
          end if
        end do
      end do
!$acc update host(nvdwloc,nvdwbloc) async(rec_queue)

c#ifdef _OPENACC
c!$acc wait(rec_queue)
c!$acc host_data use_device(vdwglob)
c      call thrust_sort(vdwglob,nvdwbloc)
c!$acc end host_data
c#endif
c
!$acc end data
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return
 
      call prmem_request(vdwglobnl,nlocnl,async=.true.)

!$acc serial async(rec_queue) present(nvdwlocnl)
      nvdwlocnl = 0
!$acc end serial
c
!$acc parallel loop async(rec_queue)
!$acc&    present(ineignl,nbvdw,rad,jvdw,repart,vdwglobnl,
!$acc&  vdwlocnl,nvdwlocnl,x,y,z)
      do i = 1, nlocnl
        iglob    = ineignl(i)
        vdwcount = nbvdw(iglob)
        if (rad(jvdw(iglob)) .ne. 0) then
          call distprocpart1(iglob,rank,d,.true.,x,y,z)
          if (repart(iglob).eq.rank) d = 0
          if (d*d.le.(vbuf2/4)) then
!$acc atomic capture
             nvdwlocnl = nvdwlocnl + 1
             nvdw_cap  = nvdwlocnl
!$acc end atomic
             vdwglobnl(nvdw_cap)  = vdwcount + 1
             vdwlocnl(vdwcount+1) = nvdw_cap
          end if
        end if
      end do
!$acc update host(nvdwlocnl) async(rec_queue)

#ifdef _OPENACC
!$acc wait
!$acc host_data use_device(vdwglobnl)
      call thrust_sort(vdwglobnl,nvdwlocnl,rec_stream)
!$acc end host_data
!$acc parallel loop present(vdwglobnl,vdwlocnl) async(rec_queue)
      do i =1, nvdwlocnl
         nvdw_cap = vdwglobnl(i)
         vdwlocnl(nvdw_cap) = i
      end do
#endif
      end
c
c     Reconfigure vdw glob on host and send them on device
c
      subroutine kvdw_reset1(istep)
      use atmlst,only : vdwglob,vdwglobnl
      use domdec,only : glob,p_recep2,domlen,nlocnl,repart,rank,
     &                  nbloc,nloc,bufbeg,n_recep2
      use kvdws ,only : rad
      use neigh ,only : ineigup,vbuf2,ineignl
      use tinheader,only: ti_p
      use tinMemory,only: prmem_request
      use vdw   ,only : nbvdw,jvdw,nvdwbloc,nvdwloc,nvdwlocnl,
     &                  vdwlocnl,nvdwblocloop
      implicit none
      integer,intent(in) :: istep
      integer i,iglob,iproc,modnl
      integer vdwcount,nvdw_cap,ipr2,ibufbeg
      real(t_p) d

!$acc wait
      call prmem_request(vdwglob,nbloc,async=.false.)
c
c     remove zero-sized atoms from the list of vdw sites
c
      nvdwloc = 0
c
      do i = 1, nloc
         iglob = glob(i)
         vdwcount = nbvdw(iglob)
         if (rad(jvdw(iglob)) .ne. 0) then
            nvdwloc = nvdwloc + 1
            vdwglob(nvdwloc) = vdwcount + 1
         end if
      end do
c
      nvdwbloc = nvdwloc
c
      do iproc = 1, n_recep2
        ipr2    = p_recep2(iproc)+1 
        ibufbeg = bufbeg(ipr2)
        do i = 1, domlen(ipr2)
          iglob = glob(ibufbeg+i-1)
          vdwcount = nbvdw(iglob)
          if (rad(jvdw(iglob)) .ne. 0) then
            nvdwbloc = nvdwbloc + 1
            vdwglob(nvdwbloc) = vdwcount + 1
          end if
        end do
      end do

c     nvdwblocloop is nvdwbloc if nvdwbloc is a multiple of 16, or the first one greater
      nvdwblocloop = merge( nvdwbloc,
     &                     (int(nvdwbloc/16)+1)*16,
     &                     (mod(nvdwbloc,16).eq.0))
!$acc data present(vdwglob)
!$acc update device(vdwglob(:))
!$acc end data

c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return

      call prmem_request(vdwglobnl,nlocnl,async=.true.)
c     if (allocated(vdwglobnl)) deallocate(vdwglobnl)
c     allocate (vdwglobnl(nlocnl))

c
      nvdwlocnl = 0
      do i = 1, nlocnl
        iglob    = ineignl(i)
        vdwcount = nbvdw(iglob)
        if (rad(jvdw(iglob)) .ne. 0.0_ti_p) then
          call distprocpart(iglob,rank,d,.true.)
          if (repart(iglob).eq.rank) d = 0.0_ti_p
          if (d*d.le.(vbuf2/4)) then
             nvdwlocnl = nvdwlocnl + 1
             vdwglobnl(nvdwlocnl)  = vdwcount + 1
             vdwlocnl(vdwcount+1) = nvdwlocnl
          end if
        end if
      end do
!$acc data present(vdwglobnl(:),vdwlocnl(:))
!$acc update device(vdwglobnl(:),vdwlocnl(:))
!$acc end data

      end


      subroutine upload_device_kvdw
      use domdec,only: rank,hostcomm
      use inform,only: deb_Path
      use kvdws
      use mpi   ,only: MPI_BARRIER
      use tinMemory
      use vdw
      implicit none
      integer ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_vdw')
      if (deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif
!$acc update device(jvdw)
!$acc update device(radmin,radmin4,epsilon,epsilon4)
!$acc update device(ired,kred)
!$acc update device(ivt,jvt)
!$acc update device(ivdw,nbvdw)
!$acc enter data copyin(nvdwloc,nvdwbloc,nvdwlocnl)
!$acc enter data copyin(rad)
      end subroutine

      subroutine delete_data_vdw
      use domdec,only: rank
      use inform,only: deb_Path
      use kvdws
      use tinMemory
      use vdw
      implicit none
 
 12   format(2x,'delete_data_vdw')
      if(deb_Path) print 12

!$acc exit data delete(nvdwloc,nvdwbloc,nvdwlocnl)
!$acc exit data delete(rad)

      call shmem_request(jvdw, winjvdw, [0], config=mhostacc)
      call shmem_request(ivdw, winivdw, [0], config=mhostacc)
      call shmem_request(ired, winired, [0], config=mhostacc)
      call shmem_request(kred, winkred, [0], config=mhostacc)
      call shmem_request(ivt,  winivt,  [0], config=mhostacc)
      call shmem_request(jvt,  winjvt,  [0], config=mhostacc)
      call shmem_request(nbvdw,winnbvdw,[0], config=mhostacc)

      call shmem_request(radmin,  winradmin,  [0,0], config=mhostacc)
      call shmem_request(epsilon, winepsilon, [0,0], config=mhostacc)
      call shmem_request(radmin4, winradmin4, [0,0], config=mhostacc)
      call shmem_request(epsilon4,winepsilon4,[0,0], config=mhostacc)
      call shmem_request(radhbnd, winradhbnd, [0,0])
      call shmem_request(epshbnd, winepshbnd, [0,0])

      end subroutine
c
c     subroutine alloc_shared_vdw : allocate shared memory pointers for vdw
c     parameter arrays
c
      subroutine alloc_shared_vdw
      use atoms
      use domdec
      use mpi
      use sizes
      use tinMemory
      use vdw
      implicit none
      integer shape2(2)

      if (associated(jvdw).and.n.eq.size(jvdw)) return

      shape2=(/maxclass,maxclass/)

      call shmem_request(jvdw, winjvdw, [n], config=mhostacc)
      call shmem_request(ivdw, winivdw, [n], config=mhostacc)
      call shmem_request(ired, winired, [n], config=mhostacc)
      call shmem_request(kred, winkred, [n], config=mhostacc)
      call shmem_request(ivt,  winivt,  [n], config=mhostacc)
      call shmem_request(jvt,  winjvt,  [n], config=mhostacc)
      call shmem_request(nbvdw,winnbvdw,[n], config=mhostacc)

      call shmem_request(radmin,  winradmin,  shape2, config=mhostacc)
      call shmem_request(epsilon, winepsilon, shape2, config=mhostacc)
      call shmem_request(radmin4, winradmin4, shape2, config=mhostacc)
      call shmem_request(epsilon4,winepsilon4,shape2, config=mhostacc)
      call shmem_request(radhbnd, winradhbnd, shape2)
      call shmem_request(epshbnd, winepshbnd, shape2)

      end
