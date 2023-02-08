c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ktortor  --  tors-tors parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ktortor" assigns torsion-torsion parameters to adjacent
c     torsion pairs and processes any new or changed values
c
c
#include "tinker_macro.h"
      subroutine ktortor(init)
      use atmlst
      use atmtyp
      use bitor
      use domdec
      use inform
      use iounit
      use keys
      use ktrtor
      use potent
      use tortor
#ifdef _OPENACC
      use thrust
#endif
      use utils
      use utilgpu
      implicit none
      integer i,j,k,m
      integer ia,ib,ic,id,ie
      integer ita,itb,itc,itd,ite
      integer size,next,ntt
      integer nx,ny,nxy
      integer isys
      integer iibitor,tortorcount,ntortorloc1
      integer*8 pt,pt1,pt2
      real(t_p) eps
      real(t_p) tx(maxtgrd2)
      real(t_p) ty(maxtgrd2)
      real(t_p) tf(maxtgrd2)
      real(t_p) bs(0:maxtgrd)
      real(t_p) cs(0:maxtgrd)
      real(t_p) ds(0:maxtgrd)
      real(t_p) tmp1(0:maxtgrd)
      real(t_p) tmp2(0:maxtgrd)
      real(t_p) tmp3(0:maxtgrd)
      real(t_p) tmp4(0:maxtgrd)
      real(t_p) tmp5(0:maxtgrd)
      real(t_p) tmp6(0:maxtgrd)
      real(t_p) tmp7(0:maxtgrd)
      logical header,cyclic
      character*4 pa,pb,pc,pd,pe
      character*20 blank
      character*20 keyword
      character*240 record
      character*240 string
      logical init
      integer ntortorloc_capture
c
      if (init) then
c
c     process keywords containing torsion-torsion parameters
c
        if (deb_Path) print*,'ktortor init'
        blank = '                    '
        header = .true.
        isys   = 0
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:8) .eq. 'TORTORS ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              ie = 0
              nx = 0
              ny = 0
              nxy = 0
              do j = 1, maxtgrd2
                 tx(j) = 0.0_ti_p
                 ty(j) = 0.0_ti_p
                 tf(j) = 0.0_ti_p
              end do
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,id,ie,nx,ny
              nxy = nx * ny
              do j = 1, nxy
                 record = keyline(i+j)
                 read (record,*,err=10,end=10)  tx(j),ty(j),tf(j)
              end do
   10         continue
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,'Additional Torsion-Torsion Parameters :',
     &                      //,5x,'Atom Classes',12x,'GridSize1',
     &                         5x,'GridSize2',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,id,ie,nx,ny
   30            format (1x,5i4,6x,i8,6x,i8)
              end if
c             size = 4
c             call numeral (ia,pa,size)
c             call numeral (ib,pb,size)
c             call numeral (ic,pc,size)
c             call numeral (id,pd,size)
c             call numeral (ie,pe,size)
c             pt = pa//pb//pc//pd//pe
              call front_convert_base(ia,ib,ic,id,ie,pt)
              do j = 1, maxntt
                 if (ktt(j).eq.-1 .or. ktt(j).eq.pt) then
                    ktt(j) = pt
                    nx = nxy
                    call sort9 (nx,tx)
                    ny = nxy
                    call sort9 (ny,ty)
                    tnx(j) = nx
                    tny(j) = ny
                    do k = 1, nx
                       ttx(k,j) = tx(k)
                    end do
                    do k = 1, ny
                       tty(k,j) = ty(k)
                    end do
                    do k = 1, nxy
                       tbf(k,j) = tf(k)
                    end do
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KTORTOR  --  Too many Torsion-Torsion',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        ntt = maxntt
        do i = maxntt, 1, -1
           if (ktt(i) .eq. -1)  ntt = i - 1
        end do
c
c       check whether each torsion-torsion parameter is periodic;
c       assumes the "tbf" array is sorted with both indices in
c       increasing order and the first index changing most rapidly
c
        do i = 1, ntt
           cyclic = .true.
           eps = 0.000001_ti_p
           nx = tnx(i) - 1
           ny = tny(i) - 1
           if (abs(abs(ttx(1,i)-ttx(tnx(i),i))-360.0_ti_p) .gt. eps)
     &        cyclic = .false.
           if (abs(abs(tty(1,i)-tty(tny(i),i))-360.0_ti_p) .gt. eps)
     &        cyclic = .false.
           if (cyclic) then
              do j = 1, tny(i)
                 k = (j-1)*tnx(i) + 1
                 if (abs(tbf(k,i)-tbf(k+nx,i)) .gt. eps) then
                    if (rank.eq.0) write (iout,60)  tbf(k,i),tbf(k+nx,i)
   60               format (/,' KTORTOR  --  Warning, Unequal Tor-Tor',
     &                          ' Values',3x,2f12.5)
                 end if
              end do
              k = ny * tnx(i)
              do j = 1, tnx(i)
                 if (abs(tbf(j,i)-tbf(j+k,i)) .gt. eps) then
                    if (rank.eq.0) write (iout,70)  tbf(j,i),tbf(j+k,i)
   70               format (/,' KTORTOR  --  Warning, Unequal Tor-Tor',
     &                          ' Values',3x,2f12.5)
                 end if
              end do
           end if
c
c       spline fit the derivatives about the first torsion
c
           do j = 1, tnx(i)
              tmp1(j-1) = ttx(j,i)
           end do
           m = 0
           do j = 1, tny(i)
              do k = 1, tnx(i)
                 tmp2(k-1) = tbf(m+k,i)
              end do
              if (cyclic) then
                 call cspline (nx,tmp1,tmp2,bs,cs,ds,tmp3,
     &                           tmp4,tmp5,tmp6,tmp7)
              else
                 call nspline (nx,tmp1,tmp2,bs,cs,tmp3,
     &                           tmp4,tmp5,tmp6,tmp7)
              end if
              do k = 1, tnx(i)
                 tbx(m+k,i) = bs(k-1)
              end do
              m = m + tnx(i)
           end do
c
c       spline fit the derivatives about the second torsion
c
           do j = 1, tny(i)
              tmp1(j-1) = tty(j,i)
           end do
           m = 1
           do j = 1, tnx(i)
              do k = 1, tny(i)
                 tmp2(k-1) = tbf(m+(k-1)*tnx(i),i)
              end do
              if (cyclic) then
                 call cspline (ny,tmp1,tmp2,bs,cs,ds,tmp3,
     &                           tmp4,tmp5,tmp6,tmp7)
              else
                 call nspline (ny,tmp1,tmp2,bs,cs,tmp3,
     &                           tmp4,tmp5,tmp6,tmp7)
              end if
              do k = 1, tny(i)
                 tby(m+(k-1)*tnx(i),i) = bs(k-1)
              end do
              m = m + 1
           end do
c
c       spline fit the cross derivatives about both torsions
c
           m = 1
           do j = 1, tnx(i)
              do k = 1, tny(i)
                 tmp2(k-1) = tbx(m+(k-1)*tnx(i),i)
              end do
              if (cyclic) then
                 call cspline (ny,tmp1,tmp2,bs,cs,ds,tmp3,
     &                            tmp4,tmp5,tmp6,tmp7)
              else
                 call nspline (ny,tmp1,tmp2,bs,cs,tmp3,
     &                           tmp4,tmp5,tmp6,tmp7)
              end if
              do k = 1, tny(i)
                 tbxy(m+(k-1)*tnx(i),i) = bs(k-1)
              end do
              m = m + 1
           end do
        end do
c
c       assign torsion-torsion parameters for each bitorsion
c
        ntortor = 0
        do i = 1, nbitor
           ia = ibitor(1,i)
           ib = ibitor(2,i)
           ic = ibitor(3,i)
           id = ibitor(4,i)
           ie = ibitor(5,i)
           nbtortor(i) = ntortor
           ita = class(ia)
           itb = class(ib)
           itc = class(ic)
           itd = class(id)
           ite = class(ie)
c          size = 4
c          call numeral (ita,pa,size)
c          call numeral (itb,pb,size)
c          call numeral (itc,pc,size)
c          call numeral (itd,pd,size)
c          call numeral (ite,pe,size)
c          pt1 = pa//pb//pc//pd//pe
c          pt2 = pe//pd//pc//pb//pa
           call front_convert_base(ita,itb,itc,itd,ite,pt1)
           call front_convert_base(ite,itd,itc,itb,ita,pt2)
c
c       find parameters for this torsion-torsion interaction
c
           do j = 1, ntt
              if (ktt(j) .eq. pt1) then
                 ntortor = ntortor + 1
                 itt(1,ntortor) = i
                 itt(2,ntortor) = j
                 itt(3,ntortor) = 1
                 if (.not.is_find8(ktt_sys(1),isys,pt1)) then
                    isys = isys +1
                    ktt_sys(isys) = pt1
                 end if
                 goto 80
              else if (ktt(j) .eq. pt2) then
                 ntortor = ntortor + 1
                 itt(1,ntortor) = i
                 itt(2,ntortor) = j
                 itt(3,ntortor) = -1
                 if (.not.is_find8(ktt_sys(1),isys,pt2)) then
                    isys = isys +1
                    ktt_sys(isys) = pt2
                 end if
                 goto 80
              end if
           end do
   80      continue
        end do
        ktt_sys(0) = isys

c
c       turn off the torsion-torsion potential if it is not used
c
        if (ntortor .eq. 0) use_tortor = .false.
c
c     Upload Tortor shared data parameters on Device
c
        if (use_tortor) then
           call upload_device_ktortor
        else
           call delete_data_ktortor
           return
        end if

      end if

      ntt = size_i8_to_i(ktt_sys(0))
!Wait for nbitorloc
!$acc wait
      call prmem_request(tortorglob,nbitorloc,async=.false.)

!$acc data present(ntortorloc,nbitorloc)
!$acc serial async
      ntortorloc = 0
!$acc end serial
!$acc parallel loop async
!$acc& present(bitorsglob, ibitor, class, nbtortor, ktt_sys, tortorglob)
      do i = 1, nbitorloc
         iibitor = bitorsglob(i)
         ia = ibitor(1,iibitor)
         ib = ibitor(2,iibitor)
         ic = ibitor(3,iibitor)
         id = ibitor(4,iibitor)
         ie = ibitor(5,iibitor)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         ite = class(ie)
         call front_convert_base(ita,itb,itc,itd,ite,pt1)
         call front_convert_base(ite,itd,itc,itb,ita,pt2)
         tortorcount = nbtortor(iibitor)
         ntortorloc1 = 0
c     find parameters for this torsion-torsion interaction
         do j = 1, ntt
            if ( (ktt_sys(j) .eq. pt1) 
     &       .or.(ktt_sys(j) .eq. pt2)) then
!$acc atomic capture
               ntortorloc = ntortorloc + 1
               ntortorloc_capture = ntortorloc
!$acc end atomic
               ntortorloc1 = ntortorloc1 + 1
               tortorglob(ntortorloc_capture) = tortorcount +ntortorloc1
               exit
            end if
         end do
      end do
!$acc update host(ntortorloc) async
!$acc end data

      end

      subroutine upload_device_ktortor
      use bitor
      use domdec,only: rank,hostcomm
      use ktrtor
      use mpi   ,only: MPI_BARRIER
      use sizes ,only: tinkerdebug
      use tortor
      implicit none
      integer ierr

#ifdef _OPENACC
 12   format(2x,'upload_device_kpitors')
      if(rank.eq.0.and.tinkerdebug) print 12
      call MPI_BARRIER(hostcomm,ierr)
      call malloc_device_data
#endif
!$acc update device(tbxy,tbf,tbx,tby,ttx,tty,tnx,tny)
!$acc update device(itt)
!$acc update device(nbtortor)
!$acc enter data copyin(ktt_sys)
!$acc enter data copyin(ntortorloc)
      end subroutine

      subroutine delete_data_ktortor
      use bitor
      use domdec,only: rank
      use ktrtor
      use sizes ,only: tinkerdebug
      use tinMemory
      use tortor
      implicit none

 12   format(2x,'delete_data_kpitors')
      if(rank.eq.0.and.tinkerdebug) print 12

      call free_device_data
      call shmem_request(itt,     winitt,   [0,0],config=mhostacc)
      call shmem_request(nbtortor,winnbtortor,[0],config=mhostacc)
!$acc exit data delete(ktt_sys)
!$acc exit data delete(ntortorloc)
      end subroutine
