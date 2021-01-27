
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kgeom  --  restraint term parameter assignment  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kgeom" asisgns parameters for geometric restraint terms
c     to be included in the potential energy calculation
c
c
#include "tinker_precision.h"
      subroutine kgeom(init)
      use atmlst
      use argue  ,only: arg
      use atmtyp
      use atoms
      use bound
      use couple
      use domdec
      use files  ,only: filename,leng
      use group
      use iounit
      use keys
      use kgeoms
      use molcul
      use USampling ,only: init_USampling,US_enable
      use potent
      use utilgpu,only: prmem_request
      implicit none
      integer i,j,k,ncap
      integer ip,next
      integer ia,ib,ic,id
      integer, allocatable :: rpos(:)
      integer nrpos,l,sizegroup,n_capture
      real(t_p) p1,p2,p3,p4,p5
      real(t_p) d1,d2,d3
      real(t_p) a1,a2,a3
      real(t_p) t1,t2,t3
      real(t_p) g1,g2,g3
      real(t_p) xr,yr,zr
      real(t_p) xcm,ycm,zcm
      real(t_p) geometry,weigh
      real(t_p) xad,yad,zad
      real(t_p) xbd,ybd,zbd
      real(t_p) xcd,ycd,zcd
      real(t_p) c1,c2,c3
      real(t_p) vol,ratio
      logical exist,keep
      logical intermol
      character*1 letter
      character*20 keyword
      character*240 record
      character*240 string
      logical init,docompute
      real(t_p) xa,ya,za,xb,yb,zb
      real(t_p) pos(3,4)
      real(t_p), allocatable :: posgroup(:,:)
c
!$acc update host(x,y,z) async
!$acc wait
      if (init) then
c
c     allocate global arrays
c
        if (rank.eq.0.and.tinkerdebug) print*,'kgeom'
        call alloc_shared_geom
c
c       set the default values for the restraint variables
c
        npfix = 0
        ndfix = 0
        nafix = 0
        ntfix = 0
        ngfix = 0
        nchir = 0
        depth = 0.0_ti_p
        width = 0.0_ti_p
        use_basin = .false.
        use_wall = .false.
c
c       allocate local arrays
c
        allocate (rpos(n))
        nrpos = 0
        do i = 1, n
           rpos(i) = 0
        end do
c
c       search the keywords for restraint parameters
c
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           string = record(next:240)
c
c       get atom restrained to a specified position range
c
           if (keyword(1:18) .eq. 'RESTRAIN-POSITION ') then
c              read (string,*,err=10,end=10)  (rpos(l),l=nrpos+1,n)
c   10         continue
c              do while (rpos(nrpos+1) .ne. 0)
c                 nrpos = nrpos + 1
c                 rpos(nrpos) = max(-n,min(n,rpos(nrpos)))
c              end do
              p1 = 0.0_ti_p
              p2 = 0.0_ti_p
              p3 = 0.0_ti_p
              p4 = 0.0_ti_p
              p5 = 0.0_ti_p
              next = 1
              call getword (string,letter,next)
              if (letter .eq. ' ') then
                 read (string,*,err=10,end=10)  ip,p1,p2,p3,p4,p5
   10            continue
                 if (p4 .eq. 0.0_ti_p)  p4 = 100.0_ti_p
                 npfix = npfix + 1
                 ipfix(npfix) = ip
                 kpfix(1,npfix) = 1
                 kpfix(2,npfix) = 1
                 kpfix(3,npfix) = 1
                 xpfix(npfix) = p1
                 ypfix(npfix) = p2
                 zpfix(npfix) = p3
                 pfix(1,npfix) = p4
                 pfix(2,npfix) = p5
              else
                 call upcase (letter)
                 read (string,*,err=20,end=20)  ip
                 string = string(next:240)
                 read (string,*,err=20,end=20)  p1,p2,p3
   20            continue
                 if (p2 .eq. 0.0_ti_p)  p2 = 100.0_ti_p
                 npfix = npfix + 1
                 ipfix(npfix) = ip
                 kpfix(1,npfix) = 0
                 kpfix(2,npfix) = 0
                 kpfix(3,npfix) = 0
                 if (letter .eq. 'X') then
                    kpfix(1,npfix) = 1
                    xpfix(npfix) = p1
                 else if (letter .eq. 'Y') then
                    kpfix(2,npfix) = 1
                    ypfix(npfix) = p1
                 else if (letter .eq. 'Z') then
                    kpfix(3,npfix) = 1
                    zpfix(npfix) = p1
                 end if
                 pfix(1,npfix) = p2
                 pfix(2,npfix) = p3
              end if
c
c       get list of atoms restrained at their initial position (equilibration phase)
c
           else if (keyword(1:14) .eq. 'RESTRAIN-LIST ') then
              read (string,*,err=11,end=11)  (rpos(l),l=nrpos+1,n)
   11         continue
              do while (rpos(nrpos+1) .ne. 0)
                 nrpos = nrpos + 1
                 rpos(nrpos) = max(-n,min(n,rpos(nrpos)))
              end do
c
c       restrain backbone atoms at their initial position (equilibration phase)
c
           else if (keyword(1:18) .eq. 'RESTRAIN-BACKBONE ') then
             do j = 1, n
               if (name(j).eq.'CA') then
                 npfix = npfix + 1
                 ipfix(npfix) = j
                 kpfix(1,npfix) = 1
                 kpfix(2,npfix) = 1
                 kpfix(3,npfix) = 1
                 xpfix(npfix) = x(j) 
                 ypfix(npfix) = y(j) 
                 zpfix(npfix) = z(j) 
                 pfix(1,npfix) = 100.0_ti_p
                 pfix(2,npfix) = 0.0_ti_p
               end if
             end do
c
c       get atoms restrained to a specified distance range
c
           else if (keyword(1:18) .eq. 'RESTRAIN-DISTANCE ') then
              d1 = 100.0_ti_p
              d2 = 0.0_ti_p
              d3 = 0.0_ti_p
              exist = .false.
              read (string,*,err=40,end=40)  ia,ib,d1,d2
              exist = .true.
   40         continue
              read (string,*,err=50,end=50)  ia,ib,d1,d2,d3
   50         continue
              if (.not. exist) then
                 xr = x(ia) - x(ib)
                 yr = y(ia) - y(ib)
                 zr = z(ia) - z(ib)
                 intermol = (molcule(ia) .ne. molcule(ib))
                 if (use_bounds .and. intermol)  call image (xr,yr,zr)
                 d2 = sqrt(xr*xr + yr*yr + zr*zr)
              end if
              if (d3 .eq. 0.0_ti_p)  d3 = d2
              ndfix = ndfix + 1
              idfix(1,ndfix) = ia
              idfix(2,ndfix) = ib
              dfix(1,ndfix) = d1
              dfix(2,ndfix) = d2
              dfix(3,ndfix) = d3
c
c       get atoms restrained to a specified angle range
c
           else if (keyword(1:15) .eq. 'RESTRAIN-ANGLE ') then
              a1 = 10.0_ti_p
              a2 = 0.0_ti_p
              a3 = 0.0_ti_p
              exist = .false.
              read (string,*,err=70,end=70)  ia,ib,ic,a1,a2
              exist = .true.
   70         continue
              read (string,*,err=80,end=80)  ia,ib,ic,a1,a2,a3
   80         continue
              if (.not. exist)  a2 = geometry (ia,ib,ic,0)
              if (a3 .eq. 0.0_ti_p)  a3 = a2
              nafix = nafix + 1
              iafix(1,nafix) = ia
              iafix(2,nafix) = ib
              iafix(3,nafix) = ic
              afix(1,nafix) = a1
              afix(2,nafix) = a2
              afix(3,nafix) = a3
c
c       get atoms restrained to a specified torsion range
c
           else if (keyword(1:17).eq.'RESTRAIN-TORSION ') then
              t1 = 1.0_ti_p
              t2 = 0.0_ti_p
              t3 = 0.0_ti_p
              exist = .false.
              read (string,*,err=100,end=100)  ia,ib,ic,id,t1,t2
              exist = .true.
  100         continue
              read (string,*,err=110,end=110)  ia,ib,ic,id,t1,t2,t3
              exist = .true.
  110         continue
              if (.not. exist)  t2 = geometry (ia,ib,ic,id)
              if (t3 .eq. 0.0_ti_p)  t3 = t2
              do while (t2 .gt. 180.0_ti_p)
                 t2 = t2 - 360.0_ti_p
              end do
              do while (t2 .lt. -180.0_ti_p)
                 t2 = t2 + 360.0_ti_p
              end do
              do while (t3 .gt. 180.0_ti_p)
                 t3 = t3 - 360.0_ti_p
              end do
              do while (t3 .lt. -180.0_ti_p)
                 t3 = t3 + 360.0_ti_p
              end do
              ntfix = ntfix + 1
              itfix(1,ntfix) = ia
              itfix(2,ntfix) = ib
              itfix(3,ntfix) = ic
              itfix(4,ntfix) = id
              tfix(1,ntfix) = t1
              tfix(2,ntfix) = t2
              tfix(3,ntfix) = t3
c
c       get groups restrained to a specified distance range
c
           else if (keyword(1:16) .eq. 'RESTRAIN-GROUPS ') then
              g1 = 100.0_ti_p
              g2 = 0.0_ti_p
              g3 = 0.0_ti_p
              exist = .false.
              read (string,*,err=130,end=130)  ia,ib,g1,g2
              exist = .true.
  130         continue
              read (string,*,err=140,end=140)  ia,ib,g1,g2,g3
  140         continue
              if (.not. exist) then
                 xcm = 0.0_ti_p
                 ycm = 0.0_ti_p
                 zcm = 0.0_ti_p
                 do j = igrp(1,ia), igrp(2,ia)
                    k = kgrp(j)
                    weigh = mass(k)
                    xcm = xcm + x(k)*weigh
                    ycm = ycm + y(k)*weigh
                    zcm = zcm + z(k)*weigh
                 end do
                 weigh = max(1.0_ti_p,grpmass(ia))
                 xr = xcm / weigh
                 yr = ycm / weigh
                 zr = zcm / weigh
                 xcm = 0.0_ti_p
                 ycm = 0.0_ti_p
                 zcm = 0.0_ti_p
                 do j = igrp(1,ib), igrp(2,ib)
                    k = kgrp(j)
                    weigh = mass(k)
                    xcm = xcm + x(k)*weigh
                    ycm = ycm + y(k)*weigh
                    zcm = zcm + z(k)*weigh
                 end do
                 weigh = max(1.0_ti_p,grpmass(ib))
                 xr = xr - xcm/weigh
                 yr = yr - ycm/weigh
                 zr = zr - zcm/weigh
                 intermol = (molcule(kgrp(igrp(1,ia))) .ne.
     &                       molcule(kgrp(igrp(1,ib))))
                 if (use_bounds .and. intermol)  call image (xr,yr,zr)
                 g2 = sqrt(xr*xr + yr*yr + zr*zr)
              end if
              if (g3 .eq. 0.0_ti_p)  g3 = g2
              ngfix = ngfix + 1
              igfix(1,ngfix) = ia
              igfix(2,ngfix) = ib
              gfix(1,ngfix) = g1
              gfix(2,ngfix) = g2
              gfix(3,ngfix) = g3
c
c       maintain chirality as found in the original input structure
c
           else if (keyword(1:18) .eq. 'ENFORCE-CHIRALITY ') then
              do j = 1, n
                 if (n12(j) .eq. 4) then
                    ia = i12(1,j)
                    ib = i12(2,j)
                    ic = i12(3,j)
                    id = i12(4,j)
                    keep = .true.
                    if (n12(ia) .eq. 1) then
                       if (type(ia) .eq. type(ib))  keep = .false.
                       if (type(ia) .eq. type(ic))  keep = .false.
                       if (type(ia) .eq. type(id))  keep = .false.
                    else if (n12(ib) .eq. 1) then
                       if (type(ib) .eq. type(ic))  keep = .false.
                       if (type(ib) .eq. type(id))  keep = .false.
                    else if (n12(ic) .eq. 1) then
                       if (type(ic) .eq. type(id))  keep = .false.
                    end if
                    if (keep) then
                       nchir = nchir + 1
                       ichir(1,nchir) = ia
                       ichir(2,nchir) = ib
                       ichir(3,nchir) = ic
                       ichir(4,nchir) = id
                       xad = x(ia) - x(id)
                       yad = y(ia) - y(id)
                       zad = z(ia) - z(id)
                       xbd = x(ib) - x(id)
                       ybd = y(ib) - y(id)
                       zbd = z(ib) - z(id)
                       xcd = x(ic) - x(id)
                       ycd = y(ic) - y(id)
                       zcd = z(ic) - z(id)
                       c1 = ybd*zcd - zbd*ycd
                       c2 = ycd*zad - zcd*yad
                       c3 = yad*zbd - zad*ybd
                       vol = xad*c1 + xbd*c2 + xcd*c3
                       ratio = abs(vol/(xad*xbd*xcd))
                       chir(1,nchir) = 10.0_ti_p
                       if (ratio .gt. 0.1_ti_p) then
                          chir(2,nchir) = 0.5_ti_p * vol
                          chir(3,nchir) = 2.0_ti_p * vol
                       else
                          chir(2,nchir) = -2.0_ti_p * abs(vol)
                          chir(3,nchir) = 2.0_ti_p * abs(vol)
                       end if
                    end if
                 end if
              end do
c
c       setup any shallow Gaussian basin restraint between atoms
c
           else if (keyword(1:6) .eq. 'BASIN ') then
              depth = 0.0_ti_p
              width = 0.0_ti_p
              read (string,*,err=160,end=160)  depth,width
  160         continue
              use_basin = .true.
              if (depth .eq. 0.0_ti_p)  use_basin = .false.
              if (width .eq. 0.0_ti_p)  use_basin = .false.
              if (depth .gt. 0.0_ti_p)  depth = -depth
c
c       setup any spherical droplet restraint between atoms
c
           else if (keyword(1:5) .eq. 'WALL ') then
              rwall = 0_ti_p
              read (string,*,err=170,end=170)  rwall
  170         continue
              if (rwall .gt. 0.0_ti_p)  use_wall = .true.
           end if
        end do
c
c       turn on the geometric restraint potential if it is used
c
c
c     set restrained list of atoms at their starting position
c
        i = 1
        do while (rpos(i) .ne. 0)
           if (i .eq. 1) then
              npfix = 0
           end if
           if (rpos(i) .gt. 0) then
              j = rpos(i)
              npfix = npfix + 1
              ipfix(npfix) = j
              kpfix(1,npfix) = 1
              kpfix(2,npfix) = 1
              kpfix(3,npfix) = 1
              xpfix(npfix) = x(j) 
              ypfix(npfix) = y(j) 
              zpfix(npfix) = z(j) 
              pfix(1,npfix) = 100.0_ti_p
              pfix(2,npfix) = 0.0_ti_p
              i = i + 1
           else
              do j = abs(rpos(i)), abs(rpos(i+1))
                  npfix = npfix + 1
                  ipfix(npfix) = j
                  kpfix(1,npfix) = 1
                  kpfix(2,npfix) = 1
                  kpfix(3,npfix) = 1
                  xpfix(npfix) = x(j) 
                  ypfix(npfix) = y(j) 
                  zpfix(npfix) = z(j) 
                  pfix(1,npfix) = 100.0_ti_p
                  pfix(2,npfix) = 0.0_ti_p
              end do
              i = i + 2
           end if
        end do
        use_geom = .false.
        if (npfix .ne. 0)  use_geom = .true.
        if (ndfix .ne. 0)  use_geom = .true.
        if (nafix .ne. 0)  use_geom = .true.
        if (ntfix .ne. 0)  use_geom = .true.
        if (ngfix .ne. 0)  use_geom = .true.
        if (nchir .ne. 0)  use_geom = .true.
        if (use_basin)  use_geom = .true.
        if (use_wall)  use_geom = .true.
c
c       deallocate (rpos)
        call MPI_barrier(hostcomm,i)

        ! Init US for Fred if required
        if (use_geom.and.US_enable) call init_USampling(nproc,arg)
c
c     Upload Geometric constrain data
c

        if (use_geom) then
           call upload_device_kgeom
        else
           call delete_data_kgeom
           return
        end if
      end if
c

      !TODO 1.2 Move this part to device. Allow mpi comput
      call prmem_request(npfixglob,nloc,async=.true.)
      call prmem_request(ndfixglob,nloc,async=.true.)
      call prmem_request(nafixglob,nloc,async=.true.)
      call prmem_request(ntfixglob,nloc,async=.true.)
      call prmem_request(ngfixglob,nloc,async=.true.)
      call prmem_request(nchirglob,nloc,async=.true.)
c

      npfixloc = 0
      if (npfix.ne.0) then
!$acc parallel loop async
!$acc&         copy(npfixloc) present(ipfix,repart,npfixglob)
      do i = 1, npfix
        ip = ipfix(i)
        if (repart(ip).eq.rank) then
!$acc atomic capture
          npfixloc = npfixloc + 1
          ncap     = npfixloc
!$acc end atomic
          npfixglob(ncap) = i
        end if
      end do
      end if

c
      if (ndfix.ne.0) then
      ndfixloc = 0
      do i = 1, ndfix
          ia = idfix(1,i)
          ib = idfix(2,i)
          xa = x(ia)
          ya = y(ia)
          za = z(ia)
          xb = x(ib)
          yb = y(ib)
          zb = z(ib)
          call midpoint(xa,ya,za,xb,yb,zb,docompute)
          if (docompute) then
            ndfixloc = ndfixloc + 1
            ndfixglob(ndfixloc) = i
          end if
      end do
!$acc update device(ndfixglob) async
      end if

c
      nafixloc = 0
      if (nafix.ne.0) then
      do i = 1, nafix
        ia = iafix(1,i)
        pos(1,1) = x(ia)
        pos(2,1) = y(ia)
        pos(3,1) = z(ia)
        ib = iafix(2,i)
        pos(1,2) = x(ib)
        pos(2,2) = y(ib)
        pos(3,2) = z(ib)
        ic = iafix(3,i)
        pos(1,3) = x(ic)
        pos(2,3) = y(ic)
        pos(3,3) = z(ic)
        call midpointgroup(pos,3,docompute)
        if (docompute) then
          nafixloc = nafixloc + 1
          nafixglob(nafixloc) = i
        end if
      end do
!$acc update device(nafixglob) async
      end if

c
      ntfixloc = 0
      if (ntfix.ne.0) then
      do i = 1, ntfix
        ia = itfix(1,i)
        pos(1,1) = x(ia)
        pos(2,1) = y(ia)
        pos(3,1) = z(ia)
        ib = itfix(2,i)
        pos(1,2) = x(ib)
        pos(2,2) = y(ib)
        pos(3,2) = z(ib)
        ic = itfix(3,i)
        pos(1,3) = x(ic)
        pos(2,3) = y(ic)
        pos(3,3) = z(ic)
        id = itfix(4,i)
        pos(1,4) = x(id)
        pos(2,4) = y(id)
        pos(3,4) = z(id)
        call midpointgroup(pos,4,docompute)
        if (docompute) then
          ntfixloc = ntfixloc + 1
          ntfixglob(ntfixloc) = i
        end if
      end do
!$acc update device(ntfixglob) async
      end if

c
      ngfixloc = 0
      if (ngfix.ne.0) then
      do i = 1, ngfix
        ia = igfix(1,i)
        sizegroup = igrp(2,ia)-igrp(1,ia)+1
        allocate (posgroup(3,sizegroup))
        l = 0
        do j = igrp(1,ia),igrp(2,ia)
          k = kgrp(j)
          l = l+1
          posgroup(1,l) = x(k)
          posgroup(2,l) = y(k)
          posgroup(3,l) = z(k)
        end do
        call midpointgroup(posgroup,sizegroup,docompute)
        if (docompute) then
          ngfixloc = ngfixloc + 1
          ngfixglob(ngfixloc) = i
        end if
        deallocate (posgroup)
      end do
!$acc update device(ngfixglob) async
      end if

c
      nchirloc = 0
      if (nchir.ne.0) then
      do i = 1, nchir
        ia = ichir(1,i)
        pos(1,1) = x(ia)
        pos(2,1) = y(ia)
        pos(3,1) = z(ia)
        ib = ichir(2,i)
        pos(1,2) = x(ib)
        pos(2,2) = y(ib)
        pos(3,2) = z(ib)
        ic = ichir(3,i)
        pos(1,3) = x(ic)
        pos(2,3) = y(ic)
        pos(3,3) = z(ic)
        id = ichir(4,i)
        pos(1,4) = x(id)
        pos(2,4) = y(id)
        pos(3,4) = z(id)
        call midpointgroup(pos,4,docompute)
        if (docompute) then
          nchirloc = nchirloc + 1
          nchirglob(nchirloc) = i
        end if
      end do
!$acc update device(nchirglob) async
      end if

      end

      subroutine upload_device_kgeom
      use domdec,only: rank,hostcomm
      use mpi   ,only: MPI_BARRIER
      use kgeoms
      use vdw
      use sizes
      use tinMemory
      implicit none
      integer ierr
#ifdef _OPENACC
 12   format(2x,'upload_device_kgeom')
      if(rank.eq.0.and.tinkerdebug) print 12
      call MPI_BARRIER(hostcomm,ierr)
#endif
!$acc update device(ipfix,xpfix,ypfix,zpfix) async
!$acc update device(kpfix,pfix,dfix,tfix,gfix) async
!$acc update device(afix,iafix,idfix,itfix,igfix) async
!$acc update device(chir,ichir) async
      end subroutine

      subroutine delete_data_kgeom
      use domdec
      use kgeoms
      use sizes
      use tinMemory
      implicit none

 12   format(2x,'delete_data_kgeom')
      if(rank.eq.0.and.tinkerdebug) print 12

      call shmem_request(xpfix,winxpfix,  [0], config=mhostacc)
      call shmem_request(ypfix,winypfix,  [0], config=mhostacc)
      call shmem_request(zpfix,winzpfix,  [0], config=mhostacc)
      call shmem_request(pfix, winpfix, [2,0], config=mhostacc)
      call shmem_request(dfix, windfix, [3,0], config=mhostacc)
      call shmem_request(afix, winafix, [3,0], config=mhostacc)
      call shmem_request(tfix, wintfix, [3,0], config=mhostacc)
      call shmem_request(gfix, wingfix, [3,0], config=mhostacc)
      call shmem_request(chir, winchir, [3,0], config=mhostacc)
      call shmem_request(ipfix,winipfix,  [0], config=mhostacc)
      call shmem_request(kpfix,winkpfix,[3,0], config=mhostacc)
      call shmem_request(idfix,winidfix,[2,0], config=mhostacc)
      call shmem_request(iafix,winiafix,[3,0], config=mhostacc)
      call shmem_request(itfix,winitfix,[4,0], config=mhostacc)
      call shmem_request(igfix,winigfix,[2,0], config=mhostacc)
      call shmem_request(ichir,winichir,[4,0], config=mhostacc)

      end subroutine
c
c     subroutine alloc_shared_geom : allocate shared memory pointers for geometric restraint
c     parameter arrays
c
      subroutine alloc_shared_geom
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use atoms
      use domdec
      use kgeoms
      use vdw
      use mpi
      use tinMemory
      implicit none

      if (associated(xpfix).and.n.eq.size(ipfix)) return
c
      call shmem_request(xpfix,winxpfix,  [n], config=mhostacc)
      call shmem_request(ypfix,winypfix,  [n], config=mhostacc)
      call shmem_request(zpfix,winzpfix,  [n], config=mhostacc)
      call shmem_request(pfix, winpfix, [2,n], config=mhostacc)
      call shmem_request(dfix, windfix, [3,n], config=mhostacc)
      call shmem_request(afix, winafix, [3,n], config=mhostacc)
      call shmem_request(tfix, wintfix, [3,n], config=mhostacc)
      call shmem_request(gfix, wingfix, [3,n], config=mhostacc)
      call shmem_request(chir, winchir, [3,n], config=mhostacc)
      call shmem_request(ipfix,winipfix,  [n], config=mhostacc)
      call shmem_request(kpfix,winkpfix,[3,n], config=mhostacc)
      call shmem_request(idfix,winidfix,[2,n], config=mhostacc)
      call shmem_request(iafix,winiafix,[3,n], config=mhostacc)
      call shmem_request(itfix,winitfix,[4,n], config=mhostacc)
      call shmem_request(igfix,winigfix,[2,n], config=mhostacc)
      call shmem_request(ichir,winichir,[4,n], config=mhostacc)

      end
