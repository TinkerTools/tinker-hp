
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
      implicit none
      integer kbis,kkpole,kglob,kkk,inl,ipoleloc,extradd,km1,km2
      integer i,ii,k,iipole,iglob,chk,tk,step,maxd,extra,maxc,minc
      integer iter,lngth,lngth_inc,loca,maxiter,curdim,n3,mindim
      integer, allocatable :: tmp(:)
      real*8  change,rad,reinit,xr,yr,zr,dist
      real*8  mindist,inc,xi,yi,zi
      logical alfit

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
      implicit none
      integer kbis,kkpole,kglob,kkk,inl,ipoleloc
      integer i,ii,k,iipole,iglob,chk,tk,step,maxd,maxc
      integer iter,lngth,lngth_inc,loca,maxiter,curdim
      integer, allocatable :: tknblist(:)
      integer, allocatable :: tmp(:)
      real*8  change,rad,reinit,xr,yr,zr,dist
      real*8  mindist,inc,xi,yi,zi
      real*8  xmin,ymin,zmin,xmax,ymax,zmax

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
      reinit = sqrt(rad)/2.0d0


      change = 9.0d9
      iter = 0
      do i = 1, 1000
c
c   begin neighbors list
c
        if(iter .eq. 0 .or. inc/km .gt. reinit) then
C          if(iter .ne. 0)print*,"reinit"
          inc = 0.0d0
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
          xr = 0.0d0
          yr = 0.0d0
          zr = 0.0d0
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

        change = 0.0d0
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
      use potent
      use shunt
      implicit none
      integer km1,km2
      integer i,ii,xp,yp,zp,k
      integer iipole,iglob
      integer curdim,proc,mdim
      real*8  xmin,ymin,zmin,xmax,ymax,zmax
      real*8  xl,yl,zl,xs,ys,zs,pos
      real*8  xi,yi,zi

      xmin = xbegproc(rank+1)
      ymin = ybegproc(rank+1)
      zmin = zbegproc(rank+1)
      xmax = xendproc(rank+1)
      ymax = yendproc(rank+1)
      zmax = zendproc(rank+1)

      xl = abs(xmax - xmin)
      yl = abs(ymax - ymin)
      zl = abs(zmax - zmin)


      if (.not. allocated(grplst))then
        if (use_pmecore) then
          km = n/(natprblk*ndir)
        else
          km = n/(natprblk*nproc)
        end if
        call pcluster
c        print*,"proc",rank,"K grid",dcx,dcy,dcz
        allocate(grplst(n))
        allocate(npergrp(km))
      end if

      xs = xl/REAL(dcx)
      ys = yl/REAL(dcy)
      zs = zl/REAL(dcz)

      npergrp = 0
      grplst = -1

      do i = 1, npoleloc
        iipole = poleglob(i)
        iglob = ipole(iipole)
        xi = x(iglob)
        yi = y(iglob)
        zi = z(iglob) 
        call image(xi,yi,zi)

        pos = xmax
        grplst(iglob) = 1
        do xp = 1, dcx
          if(xi .gt. pos) then
            grplst(iglob) = grplst(iglob) + (xp-1)
            goto 10
          end if
          pos = pos - xs
        end do
 10   continue

        pos = ymax
        do yp = 1, dcy
          if(yi .gt. pos) then
            grplst(iglob) = grplst(iglob) + (yp-1)*dcx
            goto 20
          end if
          pos = pos - ys
        end do
 20   continue

        pos = zmax
        do zp = 1, dcz
          if(zi .gt. pos) then
            grplst(iglob) = grplst(iglob) + (zp-1)*dcx*dcy
            goto 30
          end if
          pos = pos - zs
        end do
 30   continue

        npergrp(grplst(iglob)) = npergrp(grplst(iglob)) + 1
      end do
      
      maxdim=0
      do k=1,km
        if(npergrp(k) .gt. maxdim) maxdim = npergrp(k)
      end do

      if (.not. allocated(atmofst)) allocate(atmofst(n))
      call iclear(n,atmofst)
      if(allocated(klst)) deallocate(klst)
      allocate (klst(maxdim,km))
      npergrp = 0
      
C
C     intra kblock list
C
      do ii = 1, npoleloc
        iipole = poleglob(ii)
        iglob = ipole(iipole)
        k = grplst(iglob)
        npergrp(k) = npergrp(k) + 1     
        klst(npergrp(k),k) = ii
        atmofst(iglob) = npergrp(k)
      end do


      if (.not. allocated(kofst)) allocate(kofst(km))
      call iclear(km,kofst)
      kofst(1)=0
      do i = 2, km
        mdim = npergrp(i-1)*3
        kofst(i) = kofst(i-1) + mdim*(mdim+1)/2 
      end do

      maxdim = maxdim*3
      mdim = npergrp(km)*3
      curdim = kofst(km) + mdim*(mdim+1)/2 
C      print*,"Zmat length",curdim
      if (.not. allocated(zmat) .or. zdim .lt. curdim) then
        if(allocated(zmat)) deallocate(zmat)
        allocate(zmat(curdim))
        zdim = curdim
      end if
      call aclear(curdim,zmat)
     
      return
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
      implicit none
      real*8  xmin,ymin,zmin,xmax,ymax,zmax
      integer num,istep,bkm,sft,nstp
      integer key(3),list2(3),list3(3),diff,nwdiff
      real*8 list1(3)
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
      list1(3) = 1.0d0

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
      use group
      use inform
      use iounit
      use keys
      use molcul
      implicit none
      integer i,j,k
      integer next,size
      integer gnum,ga,gb
      integer, allocatable :: list(:)
      real*8 wg
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(igrp))  allocate (igrp(2,0:maxgrp))
      if (.not. allocated(grpmass))  allocate (grpmass(0:maxgrp))
      if (.not. allocated(wgrp))  allocate (wgrp(0:maxgrp,0:maxgrp))
      if (allocated(kgrp))  deallocate (kgrp)
      if (allocated(grplist))  deallocate (grplist)
      allocate (kgrp(n))
      allocate (grplist(n))
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
      wgrp = 1.0d0
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
            wg = -1.0d0
            string = record(next:240)
            read (string,*,err=40,end=40)  ga,gb,wg
   40       continue
            if (wg .lt. 0.0d0)  wg = 1.0d0
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
               wgrp(j,i) = 0.0d0
            end do
            wgrp(i,i) = 1.0d0
         end do
      end if
      if (use_inter) then
         do i = 0, ngrp
            do j = 0, ngrp
               wgrp(j,i) = 1.0d0
            end do
            wgrp(i,i) = 0.0d0
         end do
      end if
c
c     disable consideration of interactions with any empty groups
c
      do i = 0, ngrp
         size = igrp(2,i) - igrp(1,i) + 1
         if (size .eq. 0) then
            do j = 0, ngrp
               wgrp(j,i) = 0.0d0
               wgrp(i,j) = 0.0d0
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
         grpmass(i) = 0.0d0
         do j = igrp(1,i), igrp(2,i)
            grpmass(i) = grpmass(i) + mass(kgrp(j))
         end do
      end do
c
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
               if (wgrp(j,i) .ne. 0.0d0) then
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
