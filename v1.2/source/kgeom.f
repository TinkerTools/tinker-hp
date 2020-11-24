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
      subroutine kgeom(init)
      use atmlst
      use atmtyp
      use atoms
      use bound
      use couple
      use domdec
      use group
      use iounit
      use keys
      use kgeoms
      use molcul
      use potent
      implicit none
      integer i,j,k
      integer ip,next
      integer ia,ib,ic,id

      integer l,sizegroup
      real*8 p1,p2,p3,p4,p5
      real*8 d1,d2,d3
      real*8 a1,a2,a3
      real*8 t1,t2,t3
      real*8 g1,g2,g3
      real*8 xr,yr,zr
      real*8 xcm,ycm,zcm
      real*8 geometry,weigh
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 c1,c2,c3
      real*8 vol,ratio
      logical exist,keep
      logical intermol
      character*1 letter
      character*20 keyword
      character*120 record
      character*120 string
      logical init,docompute
      real*8 xa,ya,za,xb,yb,zb
      real*8 pos(3,4)
      real*8, allocatable :: posgroup(:,:)
c
      if (init) then
c
c     allocate global arrays
c
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
        depth = 0.0d0
        width = 0.0d0
        use_basin = .false.
        use_wall = .false.
cc
cc       allocate local arrays
cc
c        allocate (rpos(n))
c        nrpos = 0
c        do i = 1, n
c           rpos(i) = 0
c        end do
c
c       search the keywords for restraint parameters
c
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           string = record(next:120)
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
              p1 = 0.0d0
              p2 = 0.0d0
              p3 = 0.0d0
              p4 = 0.0d0
              p5 = 0.0d0
              next = 1
              call getword (string,letter,next)
              if (letter .eq. ' ') then
                 read (string,*,err=10,end=10)  ip,p1,p2,p3,p4,p5
   10            continue
                 if (p4 .eq. 0.0d0)  p4 = 100.0d0
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
                 string = string(next:120)
                 read (string,*,err=20,end=20)  p1,p2,p3
   20            continue
                 if (p2 .eq. 0.0d0)  p2 = 100.0d0
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
c       get atoms restrained to a specified distance range
c
           else if (keyword(1:18) .eq. 'RESTRAIN-DISTANCE ') then
              d1 = 100.0d0
              d2 = 0.0d0
              d3 = 0.0d0
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
              if (d3 .eq. 0.0d0)  d3 = d2
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
              a1 = 10.0d0
              a2 = 0.0d0
              a3 = 0.0d0
              exist = .false.
              read (string,*,err=70,end=70)  ia,ib,ic,a1,a2
              exist = .true.
   70         continue
              read (string,*,err=80,end=80)  ia,ib,ic,a1,a2,a3
   80         continue
              if (.not. exist)  a2 = geometry (ia,ib,ic,0)
              if (a3 .eq. 0.0d0)  a3 = a2
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
              t1 = 1.0d0
              t2 = 0.0d0
              t3 = 0.0d0
              exist = .false.
              read (string,*,err=100,end=100)  ia,ib,ic,id,t1,t2
              exist = .true.
  100         continue
              read (string,*,err=110,end=110)  ia,ib,ic,id,t1,t2,t3
              exist = .true.
  110         continue
              if (.not. exist)  t2 = geometry (ia,ib,ic,id)
              if (t3 .eq. 0.0d0)  t3 = t2
              do while (t2 .gt. 180.0d0)
                 t2 = t2 - 360.0d0
              end do
              do while (t2 .lt. -180.0d0)
                 t2 = t2 + 360.0d0
              end do
              do while (t3 .gt. 180.0d0)
                 t3 = t3 - 360.0d0
              end do
              do while (t3 .lt. -180.0d0)
                 t3 = t3 + 360.0d0
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
              g1 = 100.0d0
              g2 = 0.0d0
              g3 = 0.0d0
              exist = .false.
              read (string,*,err=130,end=130)  ia,ib,g1,g2
              exist = .true.
  130         continue
              read (string,*,err=140,end=140)  ia,ib,g1,g2,g3
  140         continue
              if (.not. exist) then
                 xcm = 0.0d0
                 ycm = 0.0d0
                 zcm = 0.0d0
                 do j = igrp(1,ia), igrp(2,ia)
                    k = kgrp(j)
                    weigh = mass(k)
                    xcm = xcm + x(k)*weigh
                    ycm = ycm + y(k)*weigh
                    zcm = zcm + z(k)*weigh
                 end do
                 weigh = max(1.0d0,grpmass(ia))
                 xr = xcm / weigh
                 yr = ycm / weigh
                 zr = zcm / weigh
                 xcm = 0.0d0
                 ycm = 0.0d0
                 zcm = 0.0d0
                 do j = igrp(1,ib), igrp(2,ib)
                    k = kgrp(j)
                    weigh = mass(k)
                    xcm = xcm + x(k)*weigh
                    ycm = ycm + y(k)*weigh
                    zcm = zcm + z(k)*weigh
                 end do
                 weigh = max(1.0d0,grpmass(ib))
                 xr = xr - xcm/weigh
                 yr = yr - ycm/weigh
                 zr = zr - zcm/weigh
                 intermol = (molcule(kgrp(igrp(1,ia))) .ne.
     &                       molcule(kgrp(igrp(1,ib))))
                 if (use_bounds .and. intermol)  call image (xr,yr,zr)
                 g2 = sqrt(xr*xr + yr*yr + zr*zr)
              end if
              if (g3 .eq. 0.0d0)  g3 = g2
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
                       chir(1,nchir) = 10.0d0
                       if (ratio .gt. 0.1d0) then
                          chir(2,nchir) = 0.5d0 * vol
                          chir(3,nchir) = 2.0d0 * vol
                       else
                          chir(2,nchir) = -2.0d0 * abs(vol)
                          chir(3,nchir) = 2.0d0 * abs(vol)
                       end if
                    end if
                 end if
              end do
c
c       setup any shallow Gaussian basin restraint between atoms
c
           else if (keyword(1:6) .eq. 'BASIN ') then
              depth = 0.0d0
              width = 0.0d0
              read (string,*,err=160,end=160)  depth,width
  160         continue
              use_basin = .true.
              if (depth .eq. 0.0d0)  use_basin = .false.
              if (width .eq. 0.0d0)  use_basin = .false.
              if (depth .gt. 0.0d0)  depth = -depth
c
c       setup any spherical droplet restraint between atoms
c
           else if (keyword(1:5) .eq. 'WALL ') then
              rwall = 0.0d0
              read (string,*,err=170,end=170)  rwall
  170         continue
              if (rwall .gt. 0.0d0)  use_wall = .true.
           end if
        end do
c
c       turn on the geometric restraint potential if it is used
c
c
c     set restrained atoms 
c
c        i = 1
c        do while (rpos(i) .ne. 0)
c           if (i .eq. 1) then
c              npfix = 0
c           end if
c           if (rpos(i) .gt. 0) then
c              j = rpos(i)
c              npfix = npfix + 1
c              ipfix(npfix) = j
c              kpfix(1,npfix) = 1
c              kpfix(2,npfix) = 1
c              kpfix(3,npfix) = 1
c              xpfix(npfix) = x(j) 
c              ypfix(npfix) = y(j) 
c              zpfix(npfix) = z(j) 
c              pfix(1,npfix) = 100.0d0
c              pfix(2,npfix) = 0.0d0
c              i = i + 1
c           else
c              do j = abs(rpos(i)), abs(rpos(i+1))
c                  npfix = npfix + 1
c                  ipfix(npfix) = j
c                  kpfix(1,npfix) = 1
c                  kpfix(2,npfix) = 1
c                  kpfix(3,npfix) = 1
c                  xpfix(npfix) = x(j) 
c                  ypfix(npfix) = y(j) 
c                  zpfix(npfix) = z(j) 
c                  pfix(1,npfix) = 100.0d0
c                  pfix(2,npfix) = 0.0d0
c              end do
c              i = i + 2
c           end if
c        end do
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
c        deallocate (rpos)
      end if
c
      if (allocated(npfixglob)) deallocate(npfixglob)
      allocate (npfixglob(nloc))
      if (allocated(ndfixglob)) deallocate(ndfixglob)
      allocate (ndfixglob(nloc))
      if (allocated(nafixglob)) deallocate(nafixglob)
      allocate (nafixglob(nloc))
      if (allocated(ntfixglob)) deallocate(ntfixglob)
      allocate (ntfixglob(nloc))
      if (allocated(ngfixglob)) deallocate(ngfixglob)
      allocate (ngfixglob(nloc))
      if (allocated(nchirglob)) deallocate(nchirglob)
      allocate (nchirglob(nloc))
c
      npfixloc = 0
      do i = 1, npfix
        ip = ipfix(i)
        if (repart(ip).eq.rank) then
          npfixloc = npfixloc + 1
          npfixglob(npfixloc) = i
        end if
      end do
c
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
c
      nafixloc = 0
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
c
      ntfixloc = 0
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
c
      ngfixloc = 0
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
c
      nchirloc = 0
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
      return
      end
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
      implicit none

      integer(KIND=MPI_ADDRESS_KIND) :: windowsize
      integer :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
c      if (associated(xpfix)) deallocate (xpfix)
c      if (associated(ypfix)) deallocate (ypfix)
c      if (associated(zpfix)) deallocate (zpfix)
c      if (associated(pfix)) deallocate (pfix)
c      if (associated(dfix)) deallocate (dfix)
c      if (associated(afix)) deallocate (afix)
c      if (associated(tfix)) deallocate (tfix)
c      if (associated(gfix)) deallocate (gfix)
c      if (associated(chir)) deallocate (chir)
c      if (associated(ipfix)) deallocate (ipfix)
c      if (associated(kpfix)) deallocate (kpfix)
c      if (associated(idfix)) deallocate (idfix)
c      if (associated(iafix)) deallocate (iafix)
c      if (associated(itfix)) deallocate (itfix)
c      if (associated(igfix)) deallocate (igfix)
c      if (associated(ichir)) deallocate (ichir)
      if (associated(xpfix)) then
        CALL MPI_Win_shared_query(winxpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winxpfix,ierr)
      end if
      if (associated(ypfix)) then
        CALL MPI_Win_shared_query(winypfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winypfix,ierr)
      end if
      if (associated(zpfix)) then
        CALL MPI_Win_shared_query(winzpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winzpfix,ierr)
      end if
      if (associated(pfix)) then
        CALL MPI_Win_shared_query(winpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winpfix,ierr)
      end if
      if (associated(dfix)) then
        CALL MPI_Win_shared_query(windfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(windfix,ierr)
      end if
      if (associated(afix)) then
        CALL MPI_Win_shared_query(winafix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winafix,ierr)
      end if
      if (associated(tfix)) then
        CALL MPI_Win_shared_query(wintfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wintfix,ierr)
      end if
      if (associated(gfix)) then
        CALL MPI_Win_shared_query(wingfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(wingfix,ierr)
      end if
      if (associated(chir)) then
        CALL MPI_Win_shared_query(winchir, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winchir,ierr)
      end if
      if (associated(ipfix)) then
        CALL MPI_Win_shared_query(winipfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winipfix,ierr)
      end if
      if (associated(kpfix)) then
        CALL MPI_Win_shared_query(winkpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winkpfix,ierr)
      end if
      if (associated(idfix)) then
        CALL MPI_Win_shared_query(winidfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winidfix,ierr)
      end if
      if (associated(iafix)) then
        CALL MPI_Win_shared_query(winiafix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winiafix,ierr)
      end if
      if (associated(itfix)) then
        CALL MPI_Win_shared_query(winitfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winitfix,ierr)
      end if
      if (associated(igfix)) then
        CALL MPI_Win_shared_query(winigfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winigfix,ierr)
      end if
      if (associated(ichir)) then
        CALL MPI_Win_shared_query(winichir, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winichir,ierr)
      end if
      
c
c     xpfix
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
     $  hostcomm, baseptr, winxpfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winxpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,xpfix,arrayshape)
c
c     ypfix
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
     $  hostcomm, baseptr, winypfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winypfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ypfix,arrayshape)
c
c     zpfix
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
     $  hostcomm, baseptr, winzpfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winzpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,zpfix,arrayshape)
c
c     pfix
c
      arrayshape2=(/2,n/)
      if (hostrank == 0) then
        windowsize = int(2*n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winpfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,pfix,arrayshape2)
c
c     dfix
c
      arrayshape2=(/3,n/)
      if (hostrank == 0) then
        windowsize = int(3*n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, windfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(windfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,dfix,arrayshape2)
c
c     afix
c
      arrayshape2=(/3,n/)
      if (hostrank == 0) then
        windowsize = int(3*n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winafix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winafix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,afix,arrayshape2)
c
c     tfix
c
      arrayshape2=(/3,n/)
      if (hostrank == 0) then
        windowsize = int(3*n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wintfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wintfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,tfix,arrayshape2)
c
c     gfix
c
      arrayshape2=(/3,n/)
      if (hostrank == 0) then
        windowsize = int(3*n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, wingfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(wingfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,gfix,arrayshape2)
c
c     chir 
c
      arrayshape2=(/3,n/)
      if (hostrank == 0) then
        windowsize = int(3*n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winchir, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winchir, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,chir,arrayshape2)
c
c     ipfix
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
     $  hostcomm, baseptr, winipfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winipfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ipfix,arrayshape)
c
c     kpfix
c
      arrayshape2=(/3,n/)
      if (hostrank == 0) then
        windowsize = int(3*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winkpfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winkpfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,kpfix,arrayshape2)
c
c     idfix
c
      arrayshape2=(/2,n/)
      if (hostrank == 0) then
        windowsize = int(2*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winidfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winidfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,idfix,arrayshape2)
c
c     iafix
c
      arrayshape2=(/3,n/)
      if (hostrank == 0) then
        windowsize = int(3*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winiafix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winiafix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iafix,arrayshape2)
c
c     itfix
c
      arrayshape2=(/4,n/)
      if (hostrank == 0) then
        windowsize = int(4*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winitfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winitfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,itfix,arrayshape2)
c
c     igfix
c
      arrayshape2=(/2,n/)
      if (hostrank == 0) then
        windowsize = int(2*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winigfix, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winigfix, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,igfix,arrayshape2)
c
c     ichir
c
      arrayshape2=(/4,n/)
      if (hostrank == 0) then
        windowsize = int(4*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winichir, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winichir, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ichir,arrayshape2)
c
      return
      end
