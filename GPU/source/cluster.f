
#include "tinker_precision.h"
      module clusterInl
      contains
#include "image.f.inc"
      end module

      subroutine cluster1
      use atmlst
      use atoms
      use bound
      use cell
      use couple
      use cutoff
      use domdec
      use divcon
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use tinheader ,only:ti_p
      use tinMemory
      implicit none
      integer i,ii,k,iipole,iglob,tk,step
      integer curdim
      integer, allocatable :: tmp(:)

      if (.not. allocated(means))then
        if (use_pmecore) then
          km = n/(natprblk*ndir)
        else
          km = n/(natprblk*nproc)
        end if
        step = npoleloc/(km+1)
        allocate (means(3,km))
        allocate (oldmeans(3,km))
        allocate(npergrp(km))
        tk=1
c
c       initialize the means of our blocks
c
        do ii = 1, km
          iipole = poleglob(tk)
          iglob = ipole(iipole)
          means(1,ii) = x(iglob)
          means(2,ii) = y(iglob)
          means(3,ii) = z(iglob)
          npergrp(ii) = 1
          tk = tk + step
        end do
        allocate(grplst(n))
      end if
c
c     build group list for all atoms
c
      do i = 1, n
        grplst(i) = -1
      end do

c      if(rank .eq. 0)print*,km," blocks per proc.",step

      call kmeans

      maxdim = 0
      allocate(tmp(km))

      do k = 1, km
        if(npergrp(k) .gt. maxdim) maxdim = npergrp(k)
        tmp(k) = 0
      end do
      if (.not. allocated(atmofst)) allocate(atmofst(n))
      call iclear(n,atmofst)
C
C     intra kblock list
C
      do ii = 1, npoleloc
        iipole = poleglob(ii)
        iglob = ipole(iipole)
        k = grplst(iglob)
        tmp(k) = tmp(k) + 1
        atmofst(iglob) = tmp(k)
      end do

C
C     Kofst is where each Zmat begins in the array 
C
      maxdim = npergrp(1)
      if (.not. allocated(kofst)) allocate(kofst(km))
      call iclear(km,kofst)
      kofst(1)=0
      do i = 2, km
        if(npergrp(i) .gt. maxdim) maxdim = npergrp(i)
        kofst(i) = kofst(i-1)+npergrp(i-1)*npergrp(i-1)*9 
      end do

      maxdim = maxdim*3
      curdim = kofst(km)+npergrp(km)*npergrp(km)*9
      if (.not. allocated(zmat) .or. zdim .lt. curdim) then
        if(allocated(zmat)) deallocate(zmat)
        allocate(zmat(curdim))
        zdim = curdim
      end if
      call aclear(curdim,zmat)
c
c     allocate some memory and clear the arrays:
c
      deallocate(tmp)
      return
      end


      subroutine kmeans
      use atmlst
      use atoms
      use bound
      use cell
      use couple
      use cutoff
      use divcon
      use domdec
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use tinheader
      implicit none
      integer i,ii,k,iipole,iglob,chk,tk
      integer iter,lngth,lngth_inc,loca
      integer, allocatable :: tknblist(:)
      real(t_p)  change,rad,reinit,xr,yr,zr,dist
      real(t_p)  mindist,inc,xi,yi,zi
      real(t_p)  xmin,ymin,zmin,xmax,ymax,zmax

      xmin = xbegproc(rank+1)
      ymin = ybegproc(rank+1)
      zmin = zbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymax = yendproc(rank+1)
      zmax = zendproc(rank+1)

      do k = 1, km
        oldmeans(1,k) = means(1,k)
        oldmeans(2,k) = means(2,k)
        oldmeans(3,k) = means(3,k)
      end do

      lngth = 20*npoleloc
      lngth_inc = 20*npoleloc
      allocate(knblist(lngth))
      allocate(point(npoleloc+1))
      rad = cut2
      reinit = sqrt(rad)/2.0_ti_p


      change = 9.0d9
      iter = 0
      do i = 1, 1000
c
c   begin neighbors list
c
        if(iter .eq. 0 .or. inc/km .gt. reinit) then
C          if(iter .ne. 0)print*,"reinit"
          inc = 0.0_ti_p
          loca = 1
          do ii = 1, npoleloc
            iipole = poleglob(ii)
            iglob = ipole(iipole)
            point(ii) = loca
            chk=0
            do k=1,km
              xr = x(iglob) - means(1,k)
              yr = y(iglob) - means(2,k)
              zr = z(iglob) - means(3,k)
              call image(xr,yr,zr)
              dist = xr*xr + yr*yr + zr*zr   
              if(dist .lt. rad)then
                knblist(loca)=k
                loca = loca + 1
                chk = chk + 1
                if(loca .EQ. lngth) then
                  allocate(tknblist(lngth+lngth_inc))
                  tknblist(1:lngth)=knblist
                  lngth = lngth+lngth_inc
                  deallocate(knblist)
                  allocate(knblist(lngth))
                  knblist=tknblist
                  deallocate(tknblist)
                end if
              end if
            end do
            if( chk .eq. 0) then
C              print*,"unable to find closest"
              mindist = 9.0d10
              do k=1,km
                xr = x(iglob) - means(1,k)
                yr = y(iglob) - means(2,k)
                zr = z(iglob) - means(3,k)
                dist = xr*xr + yr*yr + zr*zr
                call image(xr,yr,zr)
                if(dist .LT. mindist) then
                  mindist = dist
                  tk = k
                end if
              end do
              knblist(loca) = tk
              loca = loca + 1
              if(loca .eq. lngth) then
                allocate(tknblist(lngth+lngth_inc))
                tknblist(1:lngth)=knblist
                lngth = lngth+lngth_inc
                deallocate(knblist)
                allocate(knblist(lngth))
                knblist=tknblist
                deallocate(tknblist)
              end if
            end if
          end do
          point(npoleloc+1) = loca - 1
        end if

c
c   end neighbors list begin K-means
c
        do ii = 1, npoleloc
          iipole = poleglob(ii)
          iglob = ipole(iipole)
          mindist = 9.0d20
          do k = point(ii), point(ii+1)
            xr = x(iglob) - means(1,knblist(k))
            yr = y(iglob) - means(2,knblist(k))
            zr = z(iglob) - means(3,knblist(k))
            call image(xr,yr,zr)
            dist = xr*xr + yr*yr + zr*zr
            if(dist .lt. mindist) then
              mindist = dist
              grplst(iglob) = knblist(k)
            end if
          end do
        end do

        do k = 1, km
          npergrp(k) = 0
          xr = 0.0_ti_p
          yr = 0.0_ti_p
          zr = 0.0_ti_p
          do ii = 1, npoleloc
            iipole = poleglob(ii)
            iglob = ipole(iipole)
            if(grplst(iglob) .eq. k) then
              xi = x(iglob) - means(1,k)
              yi = y(iglob) - means(2,k)
              zi = z(iglob) - means(3,k)
              call image(xi,yi,zi)
              xr = xr + xi
              yr = yr + yi
              zr = zr + zi
              npergrp(k) = npergrp(k) + 1 
            end if
          end do
          means(1,k) = xr/REAL(npergrp(k)) + means(1,k)
          means(2,k) = yr/REAL(npergrp(k)) + means(2,k)
          means(3,k) = zr/REAL(npergrp(k)) + means(3,k)
          if(means(1,k) .lt. xmin)
     & means(1,k) = means(1,k) + xcell
          if(means(2,k) .lt. ymin) 
     & means(2,k) = means(2,k) + ycell
          if(means(3,k) .lt. zmin)
     & means(3,k) = means(3,k) + zcell
          if(means(1,k) .gt. xmax)
     & means(1,k) = means(1,k) - xcell
          if(means(2,k) .gt. ymax)
     & means(2,k) = means(2,k) - ycell
          if(means(3,k) .gt. zmax)
     & means(3,k) = means(3,k) - zcell
        end do

        change = 0.0_ti_p
        do k = 1, km
          change = change + abs(means(1,k)-oldmeans(1,k))
          change = change + abs(means(2,k)-oldmeans(2,k))
          change = change + abs(means(3,k)-oldmeans(3,k))
          oldmeans(1,k) = means(1,k)
          oldmeans(2,k) = means(2,k)
          oldmeans(3,k) = means(3,k)
        end do 
        inc = inc + change
        iter = iter + 1
        if(change .lt. 1.0d-10) exit
      end do


      deallocate(knblist)
      deallocate(point)
      return
      end


      subroutine cluster2
      use atmlst
      use atoms
      use bound
      use cell
      use clusterInl
      use couple
      use cutoff
      use divcon
      use domdec
      use iounit
      use inform ,only: deb_Path
#ifdef _OPENACC
      use interfaces ,only: initcuSolver
#endif
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use sizes ,only: tinkerdebug
#ifdef _OPENACC
      use thrust
#endif
      use tinheader
      use tinMemory
      use utils
      use utilgpu ,only: rec_queue
#ifdef _OPENACC
     &            , rec_stream
#endif
      implicit none
      integer km1,km2
      integer i,ii,xp,yp,zp,k
      integer iipole,iglob
      integer curdim,mdim,capture,igrp
      real(t_p)  xmin,ymin,zmin,xmax,ymax,zmax
      real(t_p)  xl,yl,zl,xs,ys,zs,pos
      real(t_p)  xi,yi,zi
      integer,allocatable:: Trigrp(:)
      logical,save:: first_in=.true.

      if (first_in) then
         if (allocated(grplst).or.allocated(npergrp).or.
     &       allocated(atmofst)) then
 11      format( " Error !!! Source Code alteration !!! ",/,3x,
     &           " Variables should not be allocated at this point !")
             write(*,11) 
             call fatal
         end if
         if (use_pmecore) then
           km = n/(natprblk*ndir)
         else
           km = n/(natprblk*nproc)
         end if
         call pcluster
#ifdef _OPENACC
         call initcuSolver(rec_stream)
#endif
12       format(1x,"Cluster2: Iproc",I5,3x,"K grid",3I6,3x,"km",I6)
         if (tinkerdebug) print 12 ,rank,dcx,dcy,dcz,km
         call prmem_request(grplst,n)
         call prmem_request(npergrp,km)
         call prmem_request(atmofst,n)
         first_in=.false.
      end if

      if (deb_Path) print*,'cluster2'
      xmin = xbegproc(rank+1)
      ymin = ybegproc(rank+1)
      zmin = zbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymax = yendproc(rank+1)
      zmax = zendproc(rank+1)

      xl = abs(xmax - xmin)
      yl = abs(ymax - ymin)
      zl = abs(zmax - zmin)

      xs = xl/real(dcx,t_p)
      ys = yl/real(dcy,t_p)
      zs = zl/real(dcz,t_p)
      maxdim = 0

      !npergrp = 0
      call set_to_zero1_int(npergrp,km,rec_queue)
!$acc parallel loop default(present) async
      do i = 1, n
         grplst(i) = -1
      end do

!$acc parallel loop default(present) async
      do i = 1, npoleloc
        iipole = poleglob(i)
        iglob  = ipole(iipole)
        xi     = x(iglob)
        yi     = y(iglob)
        zi     = z(iglob) 
        call image_inl(xi,yi,zi)

c       pos = xmax
c       grplst(iglob) = 1
c       do xp = 1, dcx
c         if(xi .gt. pos) then
c           grplst(iglob) = grplst(iglob) + (xp-1)
c           goto 10
c         end if
c         pos = pos - xs
c       end do
c10   continue

c       pos = ymax
c       do yp = 1, dcy
c         if(yi .gt. pos) then
c           grplst(iglob) = grplst(iglob) + (yp-1)*dcx
c           goto 20
c         end if
c         pos = pos - ys
c       end do
c20   continue

c       pos = zmax
c       do zp = 1, dcz
c         if(zi .gt. pos) then
c           grplst(iglob) = grplst(iglob) + (zp-1)*dcx*dcy
c           goto 30
c         end if
c         pos = pos - zs
c       end do
c30   continue

        xp   = int((xmax-xi)/xs)
        yp   = int((ymax-yi)/ys)
        zp   = int((zmax-zi)/zs)
        if (xp==dcx) xp = dcx-1
        if (yp==dcy) yp = dcy-1
        if (zp==dcz) zp = dcz-1
        igrp = 1 + (xp) + (yp)*dcx + (zp)*dcx*dcy
        grplst(iglob) = igrp
!$acc atomic capture
        npergrp(igrp) = npergrp(igrp) + 1
        capture       = npergrp(igrp)
!$acc end atomic
        maxdim        = max(maxdim,capture)
      end do
      
c     do k=1,km
c       if(npergrp(k) .gt. maxdim) maxdim = npergrp(k)
c     end do

      !call iclear(n,atmofst)
      call set_to_zero1_int(atmofst,n,rec_queue)
      call prmem_request(klst,maxdim,km,async=.true.)

      call set_to_zero1_int(npergrp,km,rec_queue)
c
c     intra kblock list
c
!$acc parallel loop default(present) async
      do ii = 1, npoleloc
        iipole     = poleglob(ii)
        iglob      = ipole(iipole)
        k          = grplst(iglob)
!$acc atomic capture
        npergrp(k) = npergrp(k) + 1     
        capture    = npergrp(k)
!$acc end atomic
        klst(capture,k) = ii
        atmofst(iglob)  = capture
      end do

      if (.not. allocated(kofst)) then
         call prmem_request(kofst,km,async=.true.)
      end if
      !call iclear(km,kofst)
      call set_to_zero1_int(kofst,km,rec_queue)

#ifdef _OPENACC
      allocate(Trigrp(km))
!$acc enter data create(Trigrp) async
!$acc parallel loop default(present) async
      do i = 1,km
         mdim = npergrp(i)*3
         Trigrp(i) = mdim*(mdim+1)/2
      end do
!$acc host_data use_device(Trigrp,kofst)
      call thrust_exclusive_scan(Trigrp,km,kofst,rec_stream)
!$acc end host_data
!$acc update host(kofst,npergrp) async
!$acc exit data delete(Trigrp) async
      deallocate(Trigrp)
#else
      kofst(1)=0
      do i = 2, km
        mdim     = npergrp(i-1)*3
        kofst(i) = kofst(i-1) + mdim*(mdim+1)/2 
      end do
#endif

!$acc wait
      maxdim = maxdim*3
      mdim   = npergrp(km)*3
      curdim = kofst(km) + mdim*(mdim+1)/2 

13    format(1x,"Cluster2: Zmat length",I10," Maxdim",2I10)
      if (deb_Path)
     &   print 13,curdim,maxdim,maxdim**2
      if (.not. allocated(zmat) .or. zdim .lt. curdim) then
        call prmem_request(zmat,curdim,async=.true.)
        zdim = curdim
      end if
      call set_to_zero1(zmat,curdim,rec_queue)

      end


      subroutine pcluster
      use atmlst
      use atoms
      use bound
      use boxes
      use cell
      use couple
      use cutoff
      use divcon
      use domdec
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use tinheader
      implicit none
      real(t_p)  xmin,ymin,zmin,xmax,ymax,zmax
      integer num,istep,bkm,sft,nstp
      integer key(3),list2(3),list3(3),diff,nwdiff
      real(t_p) list1(3)
      integer n1,n2,n3,i,j,res,kadd

      xmin = xbegproc(rank+1)
      ymin = ybegproc(rank+1)
      zmin = zbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymax = yendproc(rank+1)
      zmax = zendproc(rank+1)

c
c    Sort the axis by size
c
      list1(1) = abs(xmax - xmin)
      list1(2) = abs(ymax - ymin)
      list1(3) = abs(zmax - zmin)
      call sort2(3,list1,key)
c
c    normailze the axis and get difference between the longest and the two shorter
c
      list1(1) = list1(1)/list1(3)
      list1(2) = list1(2)/list1(3)
      list1(3) = 1.0_ti_p

c
c     get incriments 
c
      list3(1) = 1
      list3(2) = max(floor(list1(2)/list1(1)),1)
      list3(3) = max(floor(list1(3)/list1(1)),1)
      

c
c     get a starting point
c

      list2(3) = max(floor(km**(1.0/3.0)),1)
      list2(2) = max(floor(list1(2)*km**(1.0/3.0)),1)
      list2(1) = max(floor(list1(1)*km**(1.0/3.0)),1)

      diff = abs(km - list2(3)*list2(2)*list2(1))

c
c    build up and try to get as close as possible to
c    the desired km while keeping the number of
c    segments of each axis proportional to its length
c
      
      do while(.true.)
        do i = 1, list3(3)
          list2(3) = list2(3) + 1
          nwdiff = abs(km - list2(3)*list2(2)*list2(1))
          if(nwdiff .gt. diff)then
            list2(3) = list2(3) - 1
            goto 10
          end if
          diff = nwdiff
        end do  

        do i = 1, list3(2)
          list2(2) = list2(2) + 1
          nwdiff = abs(km - list2(3)*list2(2)*list2(1))
          if(nwdiff .gt. diff)then
            list2(2) = list2(2) - 1
            goto 10
          end if
          diff = nwdiff
        end do 

        do i = 1, list3(1)
          list2(1) = list2(1) + 1
          nwdiff = abs(km - list2(3)*list2(2)*list2(1))
          if(nwdiff .gt. diff)then
            list2(1) = list2(1) - 1
            goto 10
          end if
          diff = nwdiff
        end do     
      end do

 10     continue
        list2(3) = MAX(list2(3),1)
        list2(2) = MAX(list2(2),1)
        list2(1) = MAX(list2(1),1)
C      print*,"Requested block size", npoleloc/km
        km = list2(3)*list2(2)*list2(1)
C      print*,"Used block size", npoleloc/km


      if (key(3).eq.1) then
        if (key(2).eq.2) then
          dcx = list2(3)
          dcy = list2(2)
          dcz = list2(1)
        else
          dcx = list2(3)
          dcy = list2(1)
          dcz = list2(2)
        end if
      else if (key(3).eq.2) then
        if (key(2).eq.1) then
          dcx = list2(2)
          dcy = list2(3)
          dcz = list2(1)
        else
          dcx = list2(1)
          dcy = list2(3)
          dcz = list2(2)
        end if
      else
        if (key(2).eq.1) then
          dcx = list2(2)
          dcy = list2(1)
          dcz = list2(3)
        else
          dcx = list2(1)
          dcy = list2(2)
          dcz = list2(3)
        end if
      end if

      return
      end

c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cluster  --  set user-defined groups of atoms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cluster" gets the partitioning of the system into groups
c     and stores a list of the group to which each atom belongs
c
c
      subroutine cluster
      use atoms
      use atmtyp
      use bound
      use domdec,only:rank
      use group
      use inform
      use iounit
      use keys
      use molcul
      use tinheader
      use tinMemory
      implicit none
      integer i,j,k
      integer next,size
      integer gnum,ga,gb
      integer, allocatable :: list(:)
      real(t_p) wg
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (deb_Path) print*,'cluster'

      if (.not. allocated(igrp)) 
     &   call prmem_request(igrp,2,maxgrp,ncst=0)
      if (.not. allocated(grpmass)) 
     &   call prmem_request(grpmass,maxgrp,nst=0)
      if (.not. allocated(wgrp))
     &   allocate (wgrp(0:maxgrp,0:maxgrp))
      call prmem_request (   kgrp,n)
      call prmem_request (grplist,n)
c
c     set defaults for the group atom list and weight options
c
      use_group = .false.
      use_intra = .false.
      use_inter = .false.
      ngrp = 0
      kgrp = 0
      grplist = 0
      do i = 1, maxgrp
         igrp(1,i) = 1
         igrp(2,i) = 0
      end do
      wgrp = 1.0_ti_p
c
c     perform dynamic allocation of some local arrays
c
      size = max(20,n)
      allocate (list(size))
c
c     get any keywords containing atom group definitions
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'GROUP ') then
            use_group = .true.
            gnum = 0
            do i = 1, 20
               list(i) = 0
            end do
            call getnumb (record,gnum,next)
            if (gnum .gt. maxgrp) then
               write (iout,10)
   10          format (/,' CLUSTER  --  Too many Atom Groups;',
     &                    ' Increase MAXGRP')
               call fatal
            end if
            string = record(next:240)
            read (string,*,err=20,end=20)  (list(i),i=1,20)
   20       continue
            i = 1
            do while (list(i) .ne. 0)
               if (list(i) .gt. 0) then
                  grplist(list(i)) = gnum
                  i = i + 1
               else
                  do k = abs(list(i)), abs(list(i+1))
                     grplist(k) = gnum
                  end do
                  i = i + 2
               end if
            end do
c
c     get any keywords with weights for group interactions
c
         else if (keyword(1:15) .eq. 'GROUP-MOLECULE ') then
            use_group = .true.
            use_inter = .true.
            use_intra = .false.
            if (nmol .gt. maxgrp) then
               write (iout,30)
   30          format (/,' CLUSTER  --  Too many Atom Groups;',
     &                    ' Increase MAXGRP')
               call fatal
            end if
            do i = 1, nmol
               do k = imol(1,i), imol(2,i)
                  grplist(kmol(k)) = i
               end do
            end do
c
c     get any keywords with weights for group interactions
c
         else if (keyword(1:13) .eq. 'GROUP-SELECT ') then
            ga = 0
            gb = 0
            wg = -1.0_ti_p
            string = record(next:240)
            read (string,*,err=40,end=40)  ga,gb,wg
   40       continue
            if (wg .lt. 0.0_ti_p)  wg = 1.0_ti_p
            wgrp(ga,gb) = wg
            wgrp(gb,ga) = wg
            use_inter = .false.
c
c     get keywords to select common sets of group interactions
c
         else if (keyword(1:12) .eq. 'GROUP-INTRA ') then
            use_intra = .true.
            use_inter = .false.
         else if (keyword(1:12) .eq. 'GROUP-INTER ') then
            use_inter = .true.
            use_intra = .false.
         end if
      end do
c
c     pack atoms of each group into a contiguous indexed list
c
      if (use_group) then
         do i = 1, n
            list(i) = grplist(i)
         end do
         call sort3 (n,list,kgrp)
c
c     find the first and last atom in each of the groups
c
         k = list(1)
         igrp(1,k) = 1
         do i = 1, n
            j = list(i)
            if (j .ne. k) then
               igrp(2,k) = i - 1
               igrp(1,j) = i
               k = j
            end if
            ngrp = max(j,ngrp)
         end do
         igrp(2,j) = n

c
c     sort the list of atoms in each group by atom number
c
         do i = 0, ngrp
            size = igrp(2,i) - igrp(1,i) + 1
            if (igrp(1,i) .ne. 0)
     &         call sort (size,kgrp(igrp(1,i)))
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     use only intragroup or intergroup interactions if selected
c
      if (use_intra) then
         do i = 0, ngrp
            do j = 0, ngrp
               wgrp(j,i) = 0.0_ti_p
            end do
            wgrp(i,i) = 1.0_ti_p
         end do
      end if
      if (use_inter) then
         do i = 0, ngrp
            do j = 0, ngrp
               wgrp(j,i) = 1.0_ti_p
            end do
            wgrp(i,i) = 0.0_ti_p
         end do
      end if
c
c     disable consideration of interactions with any empty groups
c
      do i = 0, ngrp
         size = igrp(2,i) - igrp(1,i) + 1
         if (size .eq. 0) then
            do j = 0, ngrp
               wgrp(j,i) = 0.0_ti_p
               wgrp(i,j) = 0.0_ti_p
            end do
         end if
      end do
cc
cc     turn off bounds and replicas for intragroup calculations
cc
c      if (use_intra) then
c         use_bounds = .false.
c         use_replica = .false.
c         call cutoffs
c      end if
c
c     compute the total mass of all atoms in each group
c
      do i = 1, ngrp
         grpmass(i) = 0.0_ti_p
         do j = igrp(1,i), igrp(2,i)
            grpmass(i) = grpmass(i) + mass(kgrp(j))
         end do
      end do
!$acc update device(grplist)
!$acc update device(grpmass,kgrp,igrp)

c     output the final list of atoms in each group
c
      if (debug .and. use_group) then
         do i = 1, ngrp
            size = igrp(2,i) - igrp(1,i) + 1
            if (size .ne. 0) then
               write (iout,50)  i
   50          format (/,' List of Atoms in Group',i3,' :',/)
               write (iout,60)  (kgrp(j),j=igrp(1,i),igrp(2,i))
   60          format (3x,10i7)
            end if
         end do
      end if
c
c     output the weights for intragroup and intergroup interactions
c
      if (debug .and. use_group) then
         header = .true.
         do i = 0, ngrp
            do j = i, ngrp
               if (wgrp(j,i) .ne. 0.0_ti_p) then
                  if (header) then
                     header = .false.
                     write (iout,70)
   70                format (/,' Active Sets of Intra- and InterGroup',
     &                          ' Interactions :',
     &                       //,11x,'Groups',15x,'Type',14x,'Weight',/)
                  end if
                  if (i .eq. j) then
                     write (iout,80)  i,j,wgrp(j,i)
   80                format (5x,2i6,12x,'IntraGroup',5x,f12.4)
                  else
                     write (iout,90)  i,j,wgrp(j,i)
   90                format (5x,2i6,12x,'InterGroup',5x,f12.4)
                  end if
               end if
            end do
         end do
      end if
      return
      end
