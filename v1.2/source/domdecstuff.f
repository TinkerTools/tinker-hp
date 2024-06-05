c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine drivermpi  --  driver for MPI related quantities ##
c     ##                            (3d spatial decomposition)        ##
c     ##                                                              ##
c     ##################################################################
c
cc
c
      subroutine drivermpi
      use atoms
      use domdec
      use inform
      use iounit
      use potent
      use mpi
      implicit none
      integer iproc, ierr
      integer total_group, direct_group, rec_group
      integer, allocatable :: direct_rank(:)
c
      if (deb_Path) write(iout,*), 'drivermpi '
c
c
      ndir = nproc - nrec
c
c     MPI : get the atoms repartition over the processes
c
c     deallocate global arrays
c
      if (allocated(glob)) deallocate(glob)
      if (allocated(loc)) deallocate(loc)
      if (allocated(repart)) deallocate(repart)
      if (allocated(domlen)) deallocate(domlen)
      if (allocated(domlenpole)) deallocate(domlenpole)
      if (allocated(domlenpolerec)) deallocate(domlenpolerec) 
      if (allocated(zbegproc)) deallocate(zbegproc)
      if (allocated(zendproc)) deallocate(zendproc)
      if (allocated(ybegproc)) deallocate(ybegproc)
      if (allocated(yendproc)) deallocate(yendproc)
      if (allocated(xbegproc)) deallocate(xbegproc)
      if (allocated(xendproc)) deallocate(xendproc)
      if (allocated(p_recep1)) deallocate(p_recep1)
      if (allocated(p_recepshort1)) deallocate(p_recepshort1)
      if (allocated(p_send1)) deallocate(p_send1)
      if (allocated(p_sendshort1)) deallocate(p_sendshort1)
      if (allocated(p_recep2)) deallocate(p_recep2)
      if (allocated(p_recepshort2)) deallocate(p_recepshort2)
      if (allocated(p_send2)) deallocate(p_send2)
      if (allocated(p_sendshort2)) deallocate(p_sendshort2)
      if (allocated(pneig_recep)) deallocate(pneig_recep)
      if (allocated(pneig_send)) deallocate(pneig_send)
      if (allocated(pbig_recep)) deallocate(pbig_recep)
      if (allocated(pbigshort_recep)) deallocate(pbigshort_recep) 
      if (allocated(pbig_send)) deallocate(pbig_send)
      if (allocated(pbigshort_send)) deallocate(pbigshort_send)
      if (allocated(precdir_recep)) deallocate(precdir_recep) 
      if (allocated(precdir_send)) deallocate(precdir_send)
      if (allocated(precdir_recep1)) deallocate(precdir_recep1)
      if (allocated(precdir_recep2)) deallocate(precdir_recep2)
      if (allocated(precdir_send1)) deallocate(precdir_send1)
      if (allocated(precdir_send2)) deallocate(precdir_send2)
      if (allocated(globrec)) deallocate(globrec)
      if (allocated(locrec)) deallocate(locrec)
      if (allocated(bufbegrec)) deallocate(bufbegrec)
      if (allocated(bufbegpole)) deallocate(bufbegpole)
      if (allocated(bufbeg)) deallocate(bufbeg)
      if (use_pmecore) then
        if (nrec.eq.0) then
          if (rank.eq.0) write(iout,*)
     $     'no cores assigned to compute reciprocal space contribution'
           call fatal
        end if
        if (nproc-nrec.lt.1) then
          if (rank.eq.0) write(iout,*)
     $     'not enough cores to compute reciprocal space contribution'
           call fatal
        end if
c
c    allocate global arrays
c
        allocate (glob(n))
        allocate (loc(n))
        allocate (globrec(n))
        allocate (locrec(n))
        allocate (repart(n))
        allocate (domlen(nproc))
        allocate (domlenpole(nproc))
        allocate (domlenpolerec(nproc))
        allocate (zbegproc(nproc))
        allocate (zendproc(nproc))
        allocate (ybegproc(nproc))
        allocate (yendproc(nproc))
        allocate (xbegproc(nproc))
        allocate (xendproc(nproc))
        allocate(p_recep1(nproc))
        allocate(p_recepshort1(nproc))
        allocate(p_send1(nproc))
        allocate(p_sendshort1(nproc))
        allocate(p_recep2(nproc))
        allocate(p_recepshort2(nproc))
        allocate(p_send2(nproc))
        allocate(p_sendshort2(nproc))
        allocate(pneig_recep(nproc))
        allocate(pneig_send(nproc))
        allocate(precdir_recep(nproc))
        allocate(precdir_send(nproc))
        allocate(precdir_recep1(nproc))
        allocate(precdir_send1(nproc))
        allocate(precdir_recep2(nproc))
        allocate(precdir_send2(nproc))
        allocate(bufbegpole(nproc))
        allocate(bufbeg(nproc))
        allocate(bufbegrec(nproc))
        allocate(pbig_recep(nproc))
        allocate(pbigshort_recep(nproc))
        allocate(pbig_send(nproc))
        allocate(pbigshort_send(nproc))
c
c
c     build the two mpi groups for the computation of the multipole interactions
c     (direct and reciprocal space) and associated communicators
c
        allocate (direct_rank(0:nproc-1))
        call MPI_Comm_group(COMM_TINKER, total_group, ierr)
        do iproc = 0, ndir - 1
          direct_rank(iproc) = iproc
        end do
        call MPI_Group_incl(total_group,ndir,direct_rank,
     $   direct_group,ierr)
        call MPI_Group_excl(total_group,ndir,direct_rank,
     $   rec_group,ierr)
        call MPI_Comm_create(COMM_TINKER,direct_group,comm_dir,ierr)
        call MPI_Comm_create(COMM_TINKER,rec_group,comm_rec,ierr)
        if (rank.le.ndir-1) then
          call MPI_COMM_RANK(comm_dir,rank_bis,ierr)
          CALL MPI_Comm_split_type(comm_dir, MPI_COMM_TYPE_SHARED, 0,
     $      MPI_INFO_NULL, hostcomm,ierr)
          CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
        else
          call MPI_COMM_RANK(comm_rec,rank_bis,ierr)
          CALL MPI_Comm_split_type(comm_rec, MPI_COMM_TYPE_SHARED, 0,
     $      MPI_INFO_NULL, hostcomm,ierr)
          CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
        end if
        deallocate (direct_rank)
c
c       call the dd load balancing routine
c
        call ddpme3d
      else
        allocate (glob(n))
        allocate (loc(n))
        allocate (repart(n))
        allocate (domlen(nproc))
        allocate (domlenpole(nproc))
        allocate (domlenpolerec(nproc))
        allocate (zbegproc(nproc))
        allocate (zendproc(nproc))
        allocate (ybegproc(nproc))
        allocate (yendproc(nproc))
        allocate (xbegproc(nproc))
        allocate (xendproc(nproc))
        allocate(p_recep1(nproc))
        allocate(p_recepshort1(nproc))
        allocate(p_send1(nproc))
        allocate(p_sendshort1(nproc))
        allocate(p_recep2(nproc))
        allocate(p_recepshort2(nproc))
        allocate(p_send2(nproc))
        allocate(p_sendshort2(nproc))
        allocate(pneig_recep(nproc))
        allocate(pneig_send(nproc))
        allocate(pbig_recep(nproc))
        allocate(pbigshort_recep(nproc))
        allocate(pbig_send(nproc))
        allocate(pbigshort_send(nproc))
        allocate(precdir_recep(nproc))
        allocate(precdir_send(nproc))
        allocate(precdir_recep1(nproc))
        allocate(precdir_send1(nproc))
        allocate(precdir_recep2(nproc))
        allocate(precdir_send2(nproc))
        allocate (globrec(n))
        allocate (locrec(n))
        allocate (bufbegrec(nproc))
        allocate (bufbegpole(nproc))
        allocate (bufbeg(nproc))
        call ddpme3d
      end if
      return
      end
c
c     subroutine allocstep: deallocate arrays and reallocate them with proper size
c     (memory distribution)
c
      subroutine allocstep
      use deriv
      use domdec
      use inform
      use iounit
      implicit none
c
      if (deb_Path) write(iout,*), 'allocstep '
c
c
      if (allocated(desum)) deallocate (desum)
      allocate (desum(3,nbloc))
      desum = 0d0
      if (allocated(decrec)) deallocate (decrec)
      allocate (decrec(3,nlocrec2))
      decrec = 0d0
      if (allocated(demrec)) deallocate (demrec)
      allocate (demrec(3,nlocrec2))
      demrec = 0d0
      if (allocated(deprec)) deallocate (deprec)
      allocate (deprec(3,nlocrec2))
      deprec = 0d0
      if (allocated(dedsprec)) deallocate (dedsprec)
      allocate (dedsprec(3,nlocrec2))
      dedsprec = 0d0
c
      if (allocated(deb)) deallocate (deb)
      allocate (deb(3,nbloc))
      if (allocated(dea)) deallocate (dea)
      allocate (dea(3,nbloc))
      if (allocated(deba)) deallocate (deba)
      allocate (deba(3,nbloc))
      if (allocated(deub)) deallocate (deub)
      allocate (deub(3,nbloc))
      if (allocated(deaa)) deallocate (deaa)
      allocate (deaa(3,nbloc))
      if (allocated(deopb)) deallocate (deopb)
      allocate (deopb(3,nbloc))
      if (allocated(deopd)) deallocate (deopd)
      allocate (deopd(3,nbloc))
      if (allocated(deid)) deallocate (deid)
      allocate (deid(3,nbloc))
      if (allocated(deit)) deallocate (deit)
      allocate (deit(3,nbloc))
      if (allocated(det)) deallocate (det)
      allocate (det(3,nbloc))
      if (allocated(dept)) deallocate (dept)
      allocate (dept(3,nbloc))
      if (allocated(deat)) deallocate (deat)
      allocate (deat(3,nbloc))
      if (allocated(debt)) deallocate (debt)
      allocate (debt(3,nbloc))
      if (allocated(dett)) deallocate (dett)
      allocate (dett(3,nbloc))
      if (allocated(dev)) deallocate (dev)
      allocate (dev(3,nbloc))
      if (allocated(der)) deallocate (der)
      allocate (der(3,nbloc))
      if (allocated(dedsp)) deallocate (dedsp)
      allocate (dedsp(3,nbloc))
      if (allocated(dect)) deallocate (dect)
      allocate (dect(3,nbloc))
      if (allocated(dec)) deallocate (dec)
      if (allocated(dec)) deallocate (dec)
      allocate (dec(3,nbloc))
      if (allocated(dem)) deallocate (dem)
      allocate (dem(3,nbloc))
      if (allocated(dep)) deallocate (dep)
      allocate (dep(3,nbloc))
      if (allocated(deg)) deallocate (deg)
      allocate (deg(3,nbloc))
      if (allocated(dex)) deallocate (dex)
      allocate (dex(3,nbloc))
      if (allocated(desmd)) deallocate (desmd)
      allocate (desmd(3,nbloc))
      if (allocated(debond)) deallocate (debond)
      allocate (debond(3,nbloc))
c
      return
      end
c
c     subroutine allocstepsrespa: deallocate arrays and reallocate them with proper size
c     (memory distribution)
c
      subroutine allocsteprespa(fast)
      use deriv
      use domdec
      use inform
      use iounit
      implicit none
      logical fast
c
      if (deb_Path) write(iout,*), 'allocsteprespa '
c
c
      if (allocated(desum)) deallocate (desum)
      allocate (desum(3,nbloc))
      if (allocated(deb)) deallocate (deb)
      allocate (deb(3,nbloc))
      if (allocated(dea)) deallocate (dea)
      allocate (dea(3,nbloc))
      if (allocated(deba)) deallocate (deba)
      allocate (deba(3,nbloc))
      if (allocated(deub)) deallocate (deub)
      allocate (deub(3,nbloc))
      if (allocated(deaa)) deallocate (deaa)
      allocate (deaa(3,nbloc))
      if (allocated(deopb)) deallocate (deopb)
      allocate (deopb(3,nbloc))
      if (allocated(deopd)) deallocate (deopd)
      allocate (deopd(3,nbloc))
      if (allocated(deid)) deallocate (deid)
      allocate (deid(3,nbloc))
      if (allocated(deit)) deallocate (deit)
      allocate (deit(3,nbloc))
      if (allocated(det)) deallocate (det)
      allocate (det(3,nbloc))
      if (allocated(dept)) deallocate (dept)
      allocate (dept(3,nbloc))
      if (allocated(deat)) deallocate (deat)
      allocate (deat(3,nbloc))
      if (allocated(debt)) deallocate (debt)
      allocate (debt(3,nbloc))
      if (allocated(dett)) deallocate (dett)
      allocate (dett(3,nbloc))
      if (allocated(decrec)) deallocate (decrec)
      allocate (decrec(3,nlocrec2))
      if (allocated(demrec)) deallocate (demrec)
      allocate (demrec(3,nlocrec2))
      if (allocated(deprec)) deallocate (deprec)
      allocate (deprec(3,nlocrec2))
      if (allocated(dedsprec)) deallocate (dedsprec)
      allocate (dedsprec(3,nlocrec2))
      if (allocated(debond)) deallocate (debond)
      allocate (debond(3,nbloc))
c
      if (.not.(fast)) then
        decrec = 0d0
        demrec = 0d0
        deprec = 0d0
        dedsprec = 0d0
      end if
c
      if (allocated(dev)) deallocate (dev)
      allocate (dev(3,nbloc))
      if (allocated(der)) deallocate (der)
      allocate (der(3,nbloc))
      if (allocated(dedsp)) deallocate (dedsp)
      allocate (dedsp(3,nbloc))
      if (allocated(dect)) deallocate (dect)
      allocate (dect(3,nbloc))
      if (allocated(dec)) deallocate (dec)
      allocate (dec(3,nbloc))
      if (allocated(dem)) deallocate (dem)
      allocate (dem(3,nbloc))
      if (allocated(dep)) deallocate (dep)
      allocate (dep(3,nbloc))
      if (allocated(deg)) deallocate (deg)
      allocate (deg(3,nbloc))
      if (allocated(dex)) deallocate (dex)
      allocate (dex(3,nbloc))
      if (allocated(desmd)) deallocate (desmd)
      allocate (desmd(3,nbloc))
c
      return
      end
c
c     subroutine ddnumber : get the number of subdivision along each axis for
c     3d spatial decomposition
c
      subroutine ddnumber(num,istep)
      use boxes
      use domdec
      use inform
      use iounit
      use keys
      implicit none
      integer num,istep
      integer key(3),list2(3)
      integer, allocatable :: d(:)
      real*8 list1(3)
      integer n1,n2,n3,i,res,next
      character*20 keyword
      character*240 record
      character*240 string
 10   format('Nx = ',I5,2x,'Ny = ',I5,2x,'Nz = ',I5,2x)
 11   format('User defined 3D decompostion ','Nx = ',I5,2x,
     $    'Ny = ',I5,2x,'Nz = ',I5,2x)
 12   format('User defined 3D decomposition not compatible with number',
     $    ' of cores ','Nx*Ny*Nz = ',I5,2x,'number of procs = ',I5,2x)
 13   format('The program will impose the 3D decomposition')
c
c
c     check for keywords containing domain decomposition parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:9) .eq. 'DECOMP3D ') then
           read (string,*,err=20,end=20)  nxdd,nydd,nzdd
 20        continue
           if (nxdd*nydd*nzdd.eq.nproc) then
             if (istep.eq.0.and.verbose) then
               if (rank.eq.0) write(iout,11) nxdd,nydd,nzdd
             end if
             return
           else
             if (rank.eq.0) then
               write(iout,12) nxdd*nydd*nzdd,nproc
               write(iout,13) 
             end if
           end if
         end if
      end do
c
      allocate (d(num))
      d = 0
c
c    Sort the axis by size
c
      list1(1) = xbox
      list1(2) = ybox
      list1(3) = zbox
      call sort2(3,list1,key)
c
c     get prime decomposition of number
c
      call prime(num,d,i)
      i = i-1
      if (i.eq.1) then
        if (key(3).eq.1) then
          nxdd = num
          nydd = 1
          nzdd = 1
        else if (key(3).eq.2) then
          nxdd = 1
          nydd = num
          nzdd = 1
        else
          nxdd = 1
          nydd = 1
          nzdd = num
        end if
      else if (i.eq.2) then
        if (key(3).eq.1) then
          if (key(2).eq.2) then
            nxdd = d(2)
            nydd = d(1)
            nzdd = 1
          else
            nxdd = d(2)
            nydd = 1
            nzdd = d(1)
          end if
        else if (key(3).eq.2) then
          if (key(2).eq.1) then
            nxdd = d(1)
            nydd = d(2)
            nzdd = 1
          else
            nxdd = 1
            nydd = d(2)
            nzdd = d(1)
          end if
        else
          if (key(3).eq.1) then
            nxdd = d(1)
            nydd = 1
            nzdd = d(2)
          else
            nxdd = 1
            nydd = d(1)
            nzdd = d(2)
          end if
        end if
      else
c
        n1 = floor(num**(1.0/3.0))
        do i = 0, n1-2
          res = mod(num,n1-i)
          if (res.eq.0) goto 30
        end do
 30     continue
        n1 = n1-i
        n2 = floor((num/n1)**(1.0/2.0))
        do i = 0, n2-2
          res = mod(num/n1,n2-i)
          if (res.eq.0) goto 40
        end do
 40     continue
        n2 = n2 - i
        n3 = num/(n1*n2)
        list2(1) = n1
        list2(2) = n2
        list2(3) = n3
        call sort(3,list2)
c
        if (list2(1).eq.1) then
c
c      try dividing first by a smaller number
c
          n1 = floor(num**(1.0/3.0)) - 1
          if (n1.eq.0) goto 70
          do i = 0, n1-2
            res = mod(num,n1-i)
            if (res.eq.0) goto 50
          end do
 50       continue
          n1 = n1-i
          n2 = floor((num/n1)**(1.0/2.0))
          do i = 0, n2-2
            res = mod(num/n1,n2-i)
            if (res.eq.0) goto 60
          end do
 60       continue
          n2 = n2 - i
          n3 = num/(n1*n2)
          list2(1) = n1
          list2(2) = n2
          list2(3) = n3
          call sort(3,list2)
        end if
c
 70     if (key(3).eq.1) then
          if (key(2).eq.2) then
            nxdd = list2(3)
            nydd = list2(2)
            nzdd = list2(1)
          else
            nxdd = list2(3)
            nydd = list2(1)
            nzdd = list2(2)
          end if
        else if (key(3).eq.2) then
          if (key(2).eq.1) then
            nxdd = list2(2)
            nydd = list2(3)
            nzdd = list2(1)
          else
            nxdd = list2(1)
            nydd = list2(3)
            nzdd = list2(2)
          end if
        else
          if (key(2).eq.1) then
            nxdd = list2(2)
            nydd = list2(1)
            nzdd = list2(3)
          else
            nxdd = list2(1)
            nydd = list2(2)
            nzdd = list2(3)
          end if
        end if
      end if
      if (istep.eq.0.and.verbose) then
        if (rank.eq.0) write(iout,*) '3D Domain Decomposition'
        if (rank.eq.0) write(iout,10) nxdd,nydd,nzdd
      end if
      deallocate (d)
      return
      end
c
c
c     subroutine distproc : get the minimum distance between two 3d domains
c
c     nb : to be modified for non cubic unit cells
      subroutine distproc(iproc1,iproc2,dist,do3d)
      use domdec
      implicit none
      integer iproc1,iproc2
      integer i,j
      real*8 dist,dist1,dist2,dist3,dist4
      real*8 disttemp,disttemp2
      real*8 x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      real*8 xtemp(8),ytemp(8),ztemp(8)
      real*8 xtempbis(8),ytempbis(8),ztempbis(8)
      real*8 xr,yr,zr
      logical do3d
      dist = 10000.0d0
c
      if (do3d) then
        x1 = xbegproc(iproc1+1)
        x2 = xendproc(iproc1+1)
        x3 = xbegproc(iproc2+1)
        x4 = xendproc(iproc2+1)
        y1 = ybegproc(iproc1+1)
        y2 = yendproc(iproc1+1)
        y3 = ybegproc(iproc2+1)
        y4 = yendproc(iproc2+1)
        z1 = zbegproc(iproc1+1)
        z2 = zendproc(iproc1+1)
        z3 = zbegproc(iproc2+1)
        z4 = zendproc(iproc2+1)
c
c         first case, "same" x,y
c
        if ((((x1.le.x3).and.(x2.ge.x3)).or.((x1.ge.x3).and.(x1.le.x4)))
     $   .and. (((y1.le.y3).and.(y2.ge.y3)).or.((y1.ge.y3)
     $   .and.(y1.le.y4))))
     $     then
          dist1 = z1-z4
          call image(0.0d0,0.0d0,dist1)
          dist2 = z3-z2
          call image(0.0d0,0.0d0,dist2)
          dist3 = z1-z3
          call image(0.0d0,0.0d0,dist3)
          dist4 = z2-z4
          call image(0.0d0,0.0d0,dist4)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist3 = abs(dist3)
          dist4 = abs(dist4)
          dist = min(dist1,dist2,dist3,dist4)
c
c         second case, "same" x,z
c
        else if ((((x1.le.x3).and.(x2.ge.x3)).or.((x1.ge.x3)
     $   .and.(x1.le.x4)))
     $   .and. (((z1.le.z3).and.(z2.ge.z3)).or.
     $   ((z1.ge.z3).and.(z1.le.z4))))
     $     then
          dist1 = y1-y4
          call image(0.0d0,dist1,0.0d0)
          dist2 = y3-y2
          call image(0.0d0,dist2,0.0d0)
          dist3 = y1-y3
          call image(0.0d0,dist3,0.0d0)
          dist4 = y2-y4
          call image(0.0d0,dist4,0.0d0)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist3 = abs(dist3)
          dist4 = abs(dist4)
          dist = min(dist1,dist2,dist3,dist4)
c
c         third case, "same" y,z
c
        else if ((((y1.le.y3).and.(y2.ge.y3)).or.((y1.ge.y3)
     $   .and.(y1.le.y4)))
     $   .and. (((z1.le.z3).and.(z2.ge.z3)).or.
     $   ((z1.ge.z3).and.(z1.le.z4))))
     $     then
          dist1 = x1-x4
          call image(dist1,0.0d0,0.0d0)
          dist2 = x3-x2
          call image(dist2,0.0d0,0.0d0)
          dist3 = x1-x3
          call image(dist3,0.0d0,0.0d0)
          dist4 = x2-x4
          call image(dist4,0.0d0,0.0d0)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist3 = abs(dist3)
          dist4 = abs(dist4)
          dist = min(dist1,dist2,dist3,dist4)
c
c    along one "edge"
c
        else if ((x1.le.x3).and.(x2.ge.x3)) then
          dist = 1000.0d0
          xtemp(1) = x3
          ytemp(1) = y1
          ztemp(1) = z1
          xtemp(2) = x3
          ytemp(2) = y2
          ztemp(2) = z1
          xtemp(3) = x3
          ytemp(3) = y1
          ztemp(3) = z2
          xtemp(4) = x3
          ytemp(4) = y2
          ztemp(4) = z2
          xtempbis(1) = x3
          ytempbis(1) = y3
          ztempbis(1) = z3
          xtempbis(2) = x3
          ytempbis(2) = y3
          ztempbis(2) = z4
          xtempbis(3) = x3
          ytempbis(3) = y4
          ztempbis(3) = z3
          xtempbis(4) = x3
          ytempbis(4) = y4
          ztempbis(4) = z4
          do i = 1, 4
            do j = 1, 4
              xr = xtemp(i) - xtempbis(j)
              yr = ytemp(i) - ytempbis(j)
              zr = ztemp(i) - ztempbis(j)
              call image(xr,yr,zr)
              disttemp2 = xr*xr + yr*yr + zr*zr
              disttemp = sqrt(disttemp2)
              if (disttemp.le.dist) dist = disttemp
            end do
          end do
        else if ((x1.le.x4).and.(x2.ge.x4)) then
          dist = 1000.0d0
          xtemp(1) = x4
          ytemp(1) = y1
          ztemp(1) = z1
          xtemp(2) = x4
          ytemp(2) = y2
          ztemp(2) = z1
          xtemp(3) = x4
          ytemp(3) = y1
          ztemp(3) = z2
          xtemp(4) = x4
          ytemp(4) = y2
          ztemp(4) = z2
          xtempbis(1) = x4
          ytempbis(1) = y3
          ztempbis(1) = z3
          xtempbis(2) = x4
          ytempbis(2) = y3
          ztempbis(2) = z4
          xtempbis(3) = x4
          ytempbis(3) = y4
          ztempbis(3) = z3
          xtempbis(4) = x4
          ytempbis(4) = y4
          ztempbis(4) = z4
          do i = 1, 4
            do j = 1, 4
              xr = xtemp(i) - xtempbis(j)
              yr = ytemp(i) - ytempbis(j)
              zr = ztemp(i) - ztempbis(j)
              call image(xr,yr,zr)
              disttemp2 = xr*xr + yr*yr + zr*zr
              disttemp = sqrt(disttemp2)
              if (disttemp.le.dist) dist = disttemp
            end do
          end do
        else if ((y1.le.y3).and.(y2.ge.y3)) then
          dist = 1000.0d0
          xtemp(1) = x1
          ytemp(1) = y3
          ztemp(1) = z1
          xtemp(2) = x2
          ytemp(2) = y3
          ztemp(2) = z1
          xtemp(3) = x1
          ytemp(3) = y3
          ztemp(3) = z2
          xtemp(4) = x2
          ytemp(4) = y3
          ztemp(4) = z2
          xtempbis(1) = x3
          ytempbis(1) = y3
          ztempbis(1) = z3
          xtempbis(2) = x3
          ytempbis(2) = y3
          ztempbis(2) = z4
          xtempbis(3) = x4
          ytempbis(3) = y3
          ztempbis(3) = z3
          xtempbis(4) = x4
          ytempbis(4) = y3
          ztempbis(4) = z4
          do i = 1, 4
            do j = 1, 4
              xr = xtemp(i) - xtempbis(j)
              yr = ytemp(i) - ytempbis(j)
              zr = ztemp(i) - ztempbis(j)
              call image(xr,yr,zr)
              disttemp2 = xr*xr + yr*yr + zr*zr
              disttemp = sqrt(disttemp2)
              if (disttemp.le.dist) dist = disttemp
            end do
          end do
        else if ((y1.le.y4).and.(y2.ge.y4)) then
          dist = 1000.0d0
          xtemp(1) = x1
          ytemp(1) = y4
          ztemp(1) = z1
          xtemp(2) = x2
          ytemp(2) = y4
          ztemp(2) = z1
          xtemp(3) = x1
          ytemp(3) = y4
          ztemp(3) = z2
          xtemp(4) = x2
          ytemp(4) = y4
          ztemp(4) = z2
          xtempbis(1) = x3
          ytempbis(1) = y4
          ztempbis(1) = z3
          xtempbis(2) = x3
          ytempbis(2) = y4
          ztempbis(2) = z4
          xtempbis(3) = x4
          ytempbis(3) = y4
          ztempbis(3) = z3
          xtempbis(4) = x4
          ytempbis(4) = y4
          ztempbis(4) = z4
          do i = 1, 4
            do j = 1, 4
              xr = xtemp(i) - xtempbis(j)
              yr = ytemp(i) - ytempbis(j)
              zr = ztemp(i) - ztempbis(j)
              call image(xr,yr,zr)
              disttemp2 = xr*xr + yr*yr + zr*zr
              disttemp = sqrt(disttemp2)
              if (disttemp.le.dist) dist = disttemp
            end do
          end do
        else if ((z1.le.z3).and.(z2.ge.z3)) then
          dist = 1000.0d0
          xtemp(1) = x1
          ytemp(1) = y1
          ztemp(1) = z3
          xtemp(2) = x2
          ytemp(2) = y1
          ztemp(2) = z3
          xtemp(3) = x1
          ytemp(3) = y2
          ztemp(3) = z3
          xtemp(4) = x2
          ytemp(4) = y2
          ztemp(4) = z3
          xtempbis(1) = x3
          ytempbis(1) = y3
          ztempbis(1) = z3
          xtempbis(2) = x3
          ytempbis(2) = y4
          ztempbis(2) = z3
          xtempbis(3) = x4
          ytempbis(3) = y3
          ztempbis(3) = z3
          xtempbis(4) = x4
          ytempbis(4) = y4
          ztempbis(4) = z3
          do i = 1, 4
            do j = 1, 4
              xr = xtemp(i) - xtempbis(j)
              yr = ytemp(i) - ytempbis(j)
              zr = ztemp(i) - ztempbis(j)
              call image(xr,yr,zr)
              disttemp2 = xr*xr + yr*yr + zr*zr
              disttemp = sqrt(disttemp2)
              if (disttemp.le.dist) dist = disttemp
            end do
          end do
        else if ((z1.le.z4).and.(z2.ge.z4)) then
          dist = 1000.0d0
          xtemp(1) = x1
          ytemp(1) = y1
          ztemp(1) = z4
          xtemp(2) = x2
          ytemp(2) = y1
          ztemp(2) = z4
          xtemp(3) = x1
          ytemp(3) = y2
          ztemp(3) = z4
          xtemp(4) = x2
          ytemp(4) = y2
          ztemp(4) = z4
          xtempbis(1) = x3
          ytempbis(1) = y3
          ztempbis(1) = z4
          xtempbis(2) = x3
          ytempbis(2) = y4
          ztempbis(2) = z4
          xtempbis(3) = x4
          ytempbis(3) = y3
          ztempbis(3) = z4
          xtempbis(4) = x4
          ytempbis(4) = y4
          ztempbis(4) = z4
          do i = 1, 4
            do j = 1, 4
              xr = xtemp(i) - xtempbis(j)
              yr = ytemp(i) - ytempbis(j)
              zr = ztemp(i) - ztempbis(j)
              call image(xr,yr,zr)
              disttemp2 = xr*xr + yr*yr + zr*zr
              disttemp = sqrt(disttemp2)
              if (disttemp.le.dist) dist = disttemp
            end do
          end do

        else
c
c       on a "corner"
c
          xtemp(1) = xbegproc(iproc1+1)
          ytemp(1) = ybegproc(iproc1+1)
          ztemp(1) = zbegproc(iproc1+1)
c
          xtemp(2) = xbegproc(iproc1+1)
          ytemp(2) = ybegproc(iproc1+1)
          ztemp(2) = zendproc(iproc1+1)
c
          xtemp(3) = xbegproc(iproc1+1)
          ytemp(3) = yendproc(iproc1+1)
          ztemp(3) = zbegproc(iproc1+1)
c
          xtemp(4) = xbegproc(iproc1+1)
          ytemp(4) = yendproc(iproc1+1)
          ztemp(4) = zendproc(iproc1+1)
c
          xtemp(5) = xendproc(iproc1+1)
          ytemp(5) = ybegproc(iproc1+1)
          ztemp(5) = zbegproc(iproc1+1)
c
          xtemp(6) = xendproc(iproc1+1)
          ytemp(6) = ybegproc(iproc1+1)
          ztemp(6) = zendproc(iproc1+1)
c
          xtemp(7) = xendproc(iproc1+1)
          ytemp(7) = yendproc(iproc1+1)
          ztemp(7) = zbegproc(iproc1+1)
c
          xtemp(8) = xendproc(iproc1+1)
          ytemp(8) = yendproc(iproc1+1)
          ztemp(8) = zendproc(iproc1+1)
c
          xtempbis(1) = xbegproc(iproc2+1)
          ytempbis(1) = ybegproc(iproc2+1)
          ztempbis(1) = zbegproc(iproc2+1)
c
          xtempbis(2) = xbegproc(iproc2+1)
          ytempbis(2) = ybegproc(iproc2+1)
          ztempbis(2) = zendproc(iproc2+1)
c
          xtempbis(3) = xbegproc(iproc2+1)
          ytempbis(3) = yendproc(iproc2+1)
          ztempbis(3) = zbegproc(iproc2+1)
c
          xtempbis(4) = xbegproc(iproc2+1)
          ytempbis(4) = yendproc(iproc2+1)
          ztempbis(4) = zendproc(iproc2+1)
c
          xtempbis(5) = xendproc(iproc2+1)
          ytempbis(5) = ybegproc(iproc2+1)
          ztempbis(5) = zbegproc(iproc2+1)
c
          xtempbis(6) = xendproc(iproc2+1)
          ytempbis(6) = ybegproc(iproc2+1)
          ztempbis(6) = zendproc(iproc2+1)
c
          xtempbis(7) = xendproc(iproc2+1)
          ytempbis(7) = yendproc(iproc2+1)
          ztempbis(7) = zbegproc(iproc2+1)
c
          xtempbis(8) = xendproc(iproc2+1)
          ytempbis(8) = yendproc(iproc2+1)
          ztempbis(8) = zendproc(iproc2+1)
          dist = 1000.0d0
          do i = 1, 8
            do j = 1, 8
              xr = xtemp(i) - xtempbis(j)
              yr = ytemp(i) - ytempbis(j)
              zr = ztemp(i) - ztempbis(j)
              call image(xr,yr,zr)
              disttemp2 = xr*xr + yr*yr + zr*zr
              disttemp = sqrt(disttemp2)
              if (disttemp.le.dist) dist = disttemp
            end do
          end do
        end if
      else
        dist1 = zbegproc(iproc1+1)-zendproc(iproc2+1)
        xr = 0.0d0
        yr = 0.0d0
        call image(xr,yr,dist1)
        dist2 = zbegproc(iproc2+1)-zendproc(iproc1+1)
        xr = 0.0d0
        yr = 0.0d0
        call image(xr,yr,dist2)
        dist1 = abs(dist1)
        dist2 = abs(dist2)
        dist = min(dist1,dist2)
      end if
      return
      end
c
c     subroutine distprocpart : get the minimum distance between a 3d domain and an atom
c
c     to be modified for non cubic unit cells
      subroutine distprocpart(i,iproc,dist,do3d)
      use atoms
      use domdec
      implicit none
      integer iproc
      integer i,j
      real*8 dist,dist1,dist2
      real*8 disttemp,disttemp2
      real*8 x1,x2,x3,y1,y2,y3,z1,z2,z3
      real*8 xtemp(8),ytemp(8),ztemp(8)
      real*8 xr,yr,zr
      logical do3d
c
      x3 = x(i)
      y3 = y(i)
      z3 = z(i)
      call image(x3,y3,z3)
c
      if (do3d) then
        x1 = xbegproc(iproc+1)
        x2 = xendproc(iproc+1)
        y1 = ybegproc(iproc+1)
        y2 = yendproc(iproc+1)
        z1 = zbegproc(iproc+1)
        z2 = zendproc(iproc+1)
c
c       deal with atom exactly on the boundary of the proc's domain
c
c       on the "x" boundary
c
        if (((y1.le.y3).and.(y2.ge.y3)).and.((z1.le.z3).and.(z2.ge.z3))
     $   .and.((x1.eq.x3).or.(x2.eq.x3))) then
          dist = 0d0
          return
        end if
c
c       on the "y" boundary
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((z1.le.z3).and.(z2.ge.z3))
     $   .and.((y1.eq.y3).or.(y2.eq.y3))) then
          dist = 0d0
          return
        end if
c
c       on the "z" boundary
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3))
     $   .and.((z1.eq.z3).or.(z2.eq.z3))) then
          dist = 0d0
          return
        end if
c
c         first case, "same" x,y
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3)))
     $     then
          dist1 = z3-z2
          call image(0.0d0,0.0d0,dist1)
          dist2 = z1-z3
          call image(0.0d0,0.0d0,dist2)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist = min(dist1,dist2)
c
c         second case, "same" x,z
c
        else if (((x1.le.x3).and.(x2.ge.x3)).and.
     $    ((z1.le.z3).and.(z2.ge.z3)))
     $     then
          dist1 = y3-y2
          call image(0.0d0,dist1,0.0d0)
          dist2 = y1-y3
          call image(0.0d0,dist2,0.0d0)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist = min(dist1,dist2)
c
c         third case, "same" y,z
c
        else if (((y1.le.y3).and.(y2.ge.y3)).and.
     $    ((z1.le.z3).and.(z2.ge.z3)))
     $     then
          dist1 = x3-x2
          call image(dist1,0.0d0,0.0d0)
          dist2 = x1-x3
          call image(dist2,0.0d0,0.0d0)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist = min(dist1,dist2)
c
c     along one "edge"
c
        else if ((x1.le.x3).and.(x2.ge.x3)) then
          dist = 1000.0d0
          xtemp(1) = x3
          ytemp(1) = y1
          ztemp(1) = z1
          xtemp(2) = x3
          ytemp(2) = y2
          ztemp(2) = z1
          xtemp(3) = x3
          ytemp(3) = y1
          ztemp(3) = z2
          xtemp(4) = x3
          ytemp(4) = y2
          ztemp(4) = z2
          do j = 1, 4
            xr = x3 - xtemp(j)
            yr = y3 - ytemp(j)
            zr = z3 - ztemp(j)
            call image(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
          end do
        else if ((y1.le.y3).and.(y2.ge.y3)) then
          dist = 1000.0d0
          xtemp(1) = x1
          ytemp(1) = y3
          ztemp(1) = z1
          xtemp(2) = x2
          ytemp(2) = y3
          ztemp(2) = z1
          xtemp(3) = x1
          ytemp(3) = y3
          ztemp(3) = z2
          xtemp(4) = x2
          ytemp(4) = y3
          ztemp(4) = z2
          do j = 1, 4
            xr = x3 - xtemp(j)
            yr = y3 - ytemp(j)
            zr = z3 - ztemp(j)
            call image(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
          end do
        else if ((z1.le.z3).and.(z2.ge.z3)) then
          dist = 1000.0d0
          xtemp(1) = x1
          ytemp(1) = y1
          ztemp(1) = z3
          xtemp(2) = x2
          ytemp(2) = y1
          ztemp(2) = z3
          xtemp(3) = x1
          ytemp(3) = y2
          ztemp(3) = z3
          xtemp(4) = x2
          ytemp(4) = y2
          ztemp(4) = z3
          do j = 1, 4
            xr = x3 - xtemp(j)
            yr = y3 - ytemp(j)
            zr = z3 - ztemp(j)
            call image(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
          end do
c
       else
c
c       on a "corner"
c
          dist = 1000.0d0
          xtemp(1) = xbegproc(iproc+1)
          ytemp(1) = ybegproc(iproc+1)
          ztemp(1) = zbegproc(iproc+1)
c
          xtemp(2) = xbegproc(iproc+1)
          ytemp(2) = ybegproc(iproc+1)
          ztemp(2) = zendproc(iproc+1)
c
          xtemp(3) = xbegproc(iproc+1)
          ytemp(3) = yendproc(iproc+1)
          ztemp(3) = zbegproc(iproc+1)
c
          xtemp(4) = xbegproc(iproc+1)
          ytemp(4) = yendproc(iproc+1)
          ztemp(4) = zendproc(iproc+1)
c
          xtemp(5) = xendproc(iproc+1)
          ytemp(5) = ybegproc(iproc+1)
          ztemp(5) = zbegproc(iproc+1)
c
          xtemp(6) = xendproc(iproc+1)
          ytemp(6) = ybegproc(iproc+1)
          ztemp(6) = zendproc(iproc+1)
c
          xtemp(7) = xendproc(iproc+1)
          ytemp(7) = yendproc(iproc+1)
          ztemp(7) = zbegproc(iproc+1)
c
          xtemp(8) = xendproc(iproc+1)
          ytemp(8) = yendproc(iproc+1)
          ztemp(8) = zendproc(iproc+1)
c
          do j = 1, 8
            xr = x3 - xtemp(j)
            yr = y3 - ytemp(j)
            zr = z3 - ztemp(j)
            call image(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
          end do
        end if
      else
        dist1 = zbegproc(iproc+1)-z3
        xr = 0.0d0
        yr = 0.0d0
        call image(xr,yr,dist1)
        dist2 = z3-zendproc(iproc+1)
        xr = 0.0d0
        yr = 0.0d0
        call image(xr,yr,dist2)
        dist1 = abs(dist1)
        dist2 = abs(dist2)
        dist = min(dist1,dist2)
      end if
      return
      end
c     subroutine ddpme: domain decomposition load balancing
c     assign atom sites to MPI processes based on a domain decomposition
c
c
      subroutine ddpme3d
      use atoms
      use boxes
      use cell
      use cutoff
      use domdec
      use inform
      use iounit
      use keys
      use neigh
      use potent
      use mpi
      implicit none
      integer i,j,k
      integer nprocloc,rankloc,iproc,iproc1
      real*8 xr,yr,zr
      real*8 dist
      real*8 mbuf,vbuf,neigbuf,bigbuf
      real*8 mshortbuf,vshortbuf,bigshortbuf
      real*8 eps1,eps2
      real*8, allocatable :: xbegproctemp(:),ybegproctemp(:)
      real*8, allocatable :: zbegproctemp(:)
      real*8, allocatable :: xendproctemp(:),yendproctemp(:)
      real*8, allocatable :: zendproctemp(:)
      integer p,q,r,numneig
      integer temp_x,temp_y,temp_z,tempproc
      integer, allocatable :: neigproc(:,:),numneigproc(:),filledproc(:)
c
      integer, allocatable :: nrecep1(:),precep1(:,:)
      integer, allocatable :: nrecep2(:),precep2(:,:)
      integer, allocatable :: nbigrecep(:),pbigrecep(:,:)
c
      integer, allocatable :: nrecepshort1(:),precepshort1(:,:)
      integer, allocatable :: nrecepshort2(:),precepshort2(:,:)
      integer, allocatable :: nbigshortrecep(:),pbigshortrecep(:,:)
c
      integer, allocatable :: nneigrecep(:),pneigrecep(:,:)
c
 1000 format(' Warning, less than 10 atoms on process number',I6,x,
     $   ' number of cores may be too high compared to the number of '
     $    'atoms')
c
      if (deb_Path) write(iout,*), 'ddpme3d '
c
c
      eps1 = 10d-10
      eps2 = 10d-8
      nbloc = 0
      nloc  = 0
      nlocrec = 0
      if (use_pmecore) then
        nprocloc = ndir
        rankloc  = rank_bis
      else
        nprocloc = nproc
        rankloc  = rank
      end if
c
      allocate (nrecep1(nproc))
      allocate (precep1(nproc,nproc))
      allocate (nrecep2(nproc))
      allocate (precep2(nproc,nproc))
      allocate (nbigrecep(nproc))
      allocate (pbigrecep(nproc,nproc))
c
c     also deal with short range non bonded forces (respa-1 like integrators)
c
      allocate (nrecepshort1(nproc))
      allocate (precepshort1(nproc,nproc))
      allocate (nrecepshort2(nproc))
      allocate (precepshort2(nproc,nproc))
      allocate (nbigshortrecep(nproc))
      allocate (pbigshortrecep(nproc,nproc))
c
      allocate (nneigrecep(nproc))
      allocate (pneigrecep(nproc,nproc))
      nrecep1 = 0
      nrecep2 = 0
      nbigrecep = 0
      nneigrecep = 0
      precep1 = 0
      precep2 = 0
      pbigrecep = 0
c
      nrecepshort1 = 0
      nrecepshort2 = 0
      nbigshortrecep = 0
      precepshort1 = 0
      precepshort2 = 0
      pbigshortrecep = 0
c
      pneigrecep = 0
c
      allocate (xbegproctemp(nproc))
      allocate (ybegproctemp(nproc))
      allocate (xendproctemp(nproc))
      allocate (yendproctemp(nproc))
      allocate (zbegproctemp(nproc))
      allocate (zendproctemp(nproc))
      allocate (neigproc(26,nproc))
      allocate (numneigproc(nproc))
      allocate (filledproc(nproc))
      neigproc = 0
      numneigproc = 0
      filledproc = 0
      xbegproc = 0d0
      ybegproc = 0d0
      zbegproc = 0d0
      xbegproctemp = 0d0
      ybegproctemp = 0d0
      zbegproctemp = 0d0
      xendproc = 0d0
      yendproc = 0d0
      zendproc = 0d0
      xendproctemp = 0d0
      yendproctemp = 0d0
      zendproctemp = 0d0
      repart = -1
      domlen = 0
      glob = 0
      loc = 0
c
c     get the number of subdivision along each axis for dd
c
      call ddnumber(nprocloc,0)
c
      nx_box = xbox/nxdd
      ny_box = ybox/nydd
      nz_box = zbox/nzdd
      do i = 0, nxdd-1
        xbegproctemp(i+1) = -xbox2 + i*nx_box
        xendproctemp(i+1) = -xbox2 + (i+1)*nx_box
      end do
      do i = 0, nydd-1
        ybegproctemp(i+1) = -ybox2 + i*ny_box
        yendproctemp(i+1) = -ybox2 + (i+1)*ny_box
      end do
      do i = 0, nzdd-1
        zbegproctemp(i+1) = -zbox2 + i*nz_box
        zendproctemp(i+1) = -zbox2 + (i+1)*nz_box
      end do
c
c     assign processes
c
      do k = 1, nzdd
        do j = 1, nydd
          do i = 1, nxdd
              iproc = (k-1)*nydd*nxdd+(j-1)*nxdd+i
              xbegproc(iproc) = xbegproctemp(i)
              xendproc(iproc) = xendproctemp(i)
              ybegproc(iproc) = ybegproctemp(j)
              yendproc(iproc) = yendproctemp(j)
              zbegproc(iproc) = zbegproctemp(k)
              zendproc(iproc) = zendproctemp(k)
          end do
        end do
      end do
c
c     count number of particules per domain
c
      domlen = 0
      do i = 1, n
        xr = x(i)
        yr = y(i)
        zr = z(i)
        call image(xr,yr,zr)
        if (abs(xr-xbox2).lt.eps1) xr = xr-eps2
        if (abs(yr-ybox2).lt.eps1) yr = yr-eps2
        if (abs(zr-zbox2).lt.eps1) zr = zr-eps2
        do iproc = 0, nprocloc-1
          if ((zr.ge.zbegproc(iproc+1)).and.
     $      (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
     $      .and.(yr.lt.yendproc(iproc+1)).and.(xr.ge.xbegproc(iproc+1))
     $      .and.(xr.lt.xendproc(iproc+1))) then
             repart(i) = iproc
            domlen(repart(i)+1) = domlen(repart(i)+1) + 1
          end if
        end do
      end do
c
c     iterative procedure to change the size of domains for load balancing
c
c   get neighbouring processes
c
      do k = 1, nzdd
        do j = 1, nydd
          do i = 1, nxdd
            iproc = (k-1)*nydd*nxdd+(j-1)*nxdd+i
            numneig = 0
            filledproc = 0
            filledproc(iproc) = 1
            do p = -1,1
              do q = -1,1
                do r = -1,1
                  if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
c
                  temp_x = p+i
                  temp_y = q+j-1
                  temp_z = r+k-1
                  if ((i.eq.1).and.(p.eq.-1)) temp_x = nxdd
                  if ((i.eq.nxdd).and.(p.eq.1)) temp_x = 1
                  if ((j.eq.1).and.(q.eq.-1)) temp_y = nydd-1
                  if ((j.eq.nydd).and.(q.eq.1)) temp_y = 0
                  if ((k.eq.1).and.(r.eq.-1)) temp_z = nzdd-1
                  if ((k.eq.nzdd).and.(r.eq.1)) temp_z = 0
                  tempproc = temp_z*nydd*nxdd+temp_y*nxdd+
     $              temp_x
                  if (filledproc(tempproc).eq.1) goto 10
                  filledproc(tempproc) = 1
                  numneig = numneig+1
                  neigproc(numneig,iproc) = tempproc
 10             continue
                end do
              end do
            end do
            numneigproc(iproc) = numneig
          end do
        end do
      end do
c
      if ((use_pmecore).and.(rank.le.ndir-1)) then
        nloc = domlen(rankloc+1)
c
c       check for low atom number in each atom for 'extreme' parallelism
c
        if ((nloc.lt.10).and.(nproc.gt.32)) then
          write(iout,1000) rank   
        end if
      else if ((use_pmecore).and.(rank.gt.ndir-1)) then
        nloc = 0
      else
        nloc = domlen(rankloc+1)
c
c       check for low atom number in each atom for 'extreme' parallelism
c
        if ((nloc.lt.10).and.(nproc.gt.32)) then
          write(iout,1000) rank   
        end if
      end if
c
c     get the processes to receive data from 
c
      p_send1 = 0
      p_recep1 = 0
      p_send2 = 0
      p_recep2 = 0
      pbig_send = 0
      pbig_recep = 0
      n_recep1 = 0
      n_recep2 = 0
      nneig_recep = 0
      nbig_recep = 0
c
      p_sendshort1 = 0
      p_recepshort1 = 0
      p_sendshort2 = 0
      p_recepshort2 = 0
      pbigshort_send = 0
      pbigshort_recep = 0
      n_recepshort1 = 0
      n_recepshort2 = 0
      nbigshort_recep = 0
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) then
        goto 80
      end if
c
c     choose cutoff depending on electrostatic interaction
c
      if ((use_mpole).or.(use_polar)) then
        mbuf = sqrt(mbuf2)
        mshortbuf = sqrt(mshortbuf2)
      else if (use_charge) then
        mbuf = sqrt(cbuf2)
        mshortbuf = sqrt(cshortbuf2)
      else
        mbuf = 0.0d0
        mshortbuf = 0.0d0
      end if
c
      vbuf = sqrt(vbuf2)+2.0d0
      vshortbuf = sqrt(vshortbuf2)+2.0d0
      neigbuf = lbuffer
c
c     get maximum cutoff value
c
      bigbuf = max(mbuf,vbuf,ddcut)
      bigshortbuf = max(mshortbuf,vshortbuf,ddcut)
c
      do iproc = 0, nprocloc-1
        do iproc1 = 0, nprocloc-1
          if (iproc.eq.iproc1) goto 70
          call distproc(iproc1,iproc,dist,.true.)
          if (dist.le.(mbuf/2)) then
            nrecep1(iproc1+1) = nrecep1(iproc1+1)+1
            precep1(nrecep1(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(vbuf/2)) then
            nrecep2(iproc1+1) = nrecep2(iproc1+1)+1
            precep2(nrecep2(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(bigbuf/2)) then
            nbigrecep(iproc1+1) = nbigrecep(iproc1+1)+1
            pbigrecep(nbigrecep(iproc1+1),iproc1+1) = iproc
          end if
c
c    also deal with short range non bonded interactions (respa-1 like integrators)
c
          if (dist.le.(mshortbuf/2)) then
            nrecepshort1(iproc1+1) = nrecepshort1(iproc1+1)+1
            precepshort1(nrecepshort1(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(vshortbuf/2)) then
            nrecepshort2(iproc1+1) = nrecepshort2(iproc1+1)+1
            precepshort2(nrecepshort2(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(bigshortbuf/2)) then
            nbigshortrecep(iproc1+1) = nbigshortrecep(iproc1+1)+1
            pbigshortrecep(nbigshortrecep(iproc1+1),iproc1+1) = iproc
          end if
c
          if (dist.le.neigbuf) then
            nneigrecep(iproc1+1) = nneigrecep(iproc1+1)+1
            pneigrecep(nneigrecep(iproc1+1),iproc1+1) = iproc
          end if
 70     continue
        end do
      end do
c
      n_recep1 = nrecep1(rankloc+1)
      p_recep1 = precep1(:,rankloc+1)
      n_recep2 = nrecep2(rankloc+1)
      p_recep2 = precep2(:,rankloc+1)
      nbig_recep = nbigrecep(rankloc+1)
      pbig_recep = pbigrecep(:,rankloc+1)
c
      n_recepshort1 = nrecepshort1(rankloc+1)
      p_recepshort1 = precepshort1(:,rankloc+1)
      n_recepshort2 = nrecepshort2(rankloc+1)
      p_recepshort2 = precepshort2(:,rankloc+1)
      nbigshort_recep = nbigshortrecep(rankloc+1)
      pbigshort_recep = pbigshortrecep(:,rankloc+1)
c      write(*,*) 'nbigshort_recep = ',nbigshort_recep
c      write(*,*) 'nbig_recep = ',nbig_recep
c      write(*,*) 'n_recepshort1 = ',n_recepshort1
c      write(*,*) 'n_recep1 = ',n_recep1
c      write(*,*) 'n_recepshort2 = ',n_recepshort2
c      write(*,*) 'n_recep2 = ',n_recep2
c
      nneig_recep = nneigrecep(rankloc+1)
      pneig_recep = pneigrecep(:,rankloc+1)
c
      n_send1 = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, nrecep1(iproc+1)
            if (precep1(i,iproc+1).eq.rankloc) then
              n_send1 = n_send1 + 1
              p_send1(n_send1) = iproc
            end if
          end do
        end if
      end do
c
c     get the processes to send data to
c
      n_send2 = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, nrecep2(iproc+1)
            if (precep2(i,iproc+1).eq.rankloc) then
              n_send2 = n_send2 + 1
              p_send2(n_send2) = iproc
            end if
          end do
        end if
      end do
c
c     get the processes to send data to
c
      nbig_send = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, nbigrecep(iproc+1)
            if (pbigrecep(i,iproc+1).eq.rankloc) then
              nbig_send = nbig_send + 1
              pbig_send(nbig_send) = iproc
            end if
          end do
        end if
      end do
c
      n_sendshort1 = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, nrecepshort1(iproc+1)
            if (precepshort1(i,iproc+1).eq.rankloc) then
              n_sendshort1 = n_sendshort1 + 1
              p_sendshort1(n_sendshort1) = iproc
            end if
          end do
        end if
      end do
c
c     get the processes to send data to
c
      n_sendshort2 = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, nrecepshort2(iproc+1)
            if (precepshort2(i,iproc+1).eq.rankloc) then
              n_sendshort2 = n_sendshort2 + 1
              p_sendshort2(n_send2) = iproc
            end if
          end do
        end if
      end do
c
c     get the processes to send data to
c
      nbigshort_send = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, nbigshortrecep(iproc+1)
            if (pbigshortrecep(i,iproc+1).eq.rankloc) then
              nbigshort_send = nbigshort_send + 1
              pbigshort_send(nbigshort_send) = iproc
            end if
          end do
        end if
      end do
c
      nneig_send = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, nneigrecep(iproc+1)
            if (pneigrecep(i,iproc+1).eq.rankloc) then
              nneig_send = nneig_send + 1
              pneig_send(nneig_send) = iproc
            end if
          end do
        end if
      end do
c
 80   call orderbuffer(.true.)
      deallocate (filledproc)
      deallocate (numneigproc)
      deallocate (neigproc)
      deallocate (xbegproctemp)
      deallocate (xendproctemp)
      deallocate (ybegproctemp)
      deallocate (yendproctemp)
      deallocate (zbegproctemp)
      deallocate (zendproctemp)
      deallocate (nrecep1)
      deallocate (precep1)
      deallocate (nrecep2)
      deallocate (precep2)
      deallocate (nbigrecep)
      deallocate (pbigrecep)
      deallocate (nrecepshort1)
      deallocate (precepshort1)
      deallocate (nrecepshort2)
      deallocate (precepshort2)
      deallocate (nbigshortrecep)
      deallocate (pbigshortrecep)
      deallocate (nneigrecep)
      deallocate (pneigrecep)
      return
      end
c
c
c     subroutine halfcell : routine that says whether an interaction between two particules
c     has to be computed within the current domain or not (dd half cell method, Newton's
c     3rd law)
c
c
      subroutine halfcell(xi,yi,zi,xj,yj,zj,docompute)
      use bound
      implicit none
      real*8 xr,yr,zr
      real*8 xi,yi,zi
      real*8 xj,yj,zj
      logical docompute
c
      docompute = .false.
c
      xr = xi - xj
      yr = yi - yj
      zr = zi - zj
      if (use_bounds) call image(xr,yr,zr)
      if (xr.lt.0.0d0) then
        docompute = .true.
      else if ((xr.eq.0.0d0).and.(yr.lt.0.0d0)) then
        docompute = .true.
      else if ((xr.eq.0.0d0).and.(yr.eq.0.0d0).and.
     $    (zr.lt.0.0d0)) then
        docompute = .true.
      end if
      return
      end
c
c     subroutine midpoint : routine that says whether an interaction between two particules
c     has to be computed within the current domain or not (dd midpoint method)
c
c
      subroutine midpoint(xi,yi,zi,xk,yk,zk,docompute)
      use cell
      use domdec
      implicit none
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 xrmid,yrmid,zrmid
      real*8 eps1,eps2
      logical docompute
c
      eps1 = 10d-10
      eps2 = 10d-8
      docompute = .false.
c
c      call image(xi,yi,zi)
c      call image(xk,yk,zk)
      xr = xi - xk
      yr = yi - yk
      zr = zi - zk
      call image(xr,yr,zr)
c     
c     definition of the middle point between i and k atoms
c
      xrmid = xk + xr/2
      yrmid = yk + yr/2
      zrmid = zk + zr/2
      call image(xrmid,yrmid,zrmid)
      if (abs(xrmid-xcell2).lt.eps1) xrmid = xrmid-eps2
      if (abs(yrmid-ycell2).lt.eps1) yrmid = yrmid-eps2
      if (abs(zrmid-zcell2).lt.eps1) zrmid = zrmid-eps2
      if ((zrmid.ge.zbegproc(rank+1)).and.
     $  (zrmid.lt.zendproc(rank+1)).and.(yrmid.ge.ybegproc(rank+1))
     $  .and.(yrmid.lt.yendproc(rank+1))
     $  .and.(xrmid.ge.xbegproc(rank+1))
     $  .and.(xrmid.lt.xendproc(rank+1))) then
        docompute = .true.
      end if
      return
      end
c
c     subroutine midpointimage : routine that says whether an interaction between two particules
c     has to be computed within the current domain or not (dd midpoint method), also returns
c     minimum image of the distance vector
c
c
      subroutine midpointimage(xi,yi,zi,xk,yk,zk,xr,yr,zr,docompute)
      use cell
      use domdec
      implicit none
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 xrmid,yrmid,zrmid
      real*8 eps1,eps2
      logical docompute
      eps1 = 10d-10
      eps2 = 10d-8
c
      docompute = .false.
c
      xr = xi - xk
      yr = yi - yk
      zr = zi - zk
      call image(xr,yr,zr)
c     
c     definition of the middle point between i and k atoms
c
      xrmid = xk + xr/2
      yrmid = yk + yr/2
      zrmid = zk + zr/2
      call image(xrmid,yrmid,zrmid)
      if (abs(xrmid-xcell2).lt.eps1) xrmid = xrmid-eps2
      if (abs(yrmid-ycell2).lt.eps1) yrmid = yrmid-eps2
      if (abs(zrmid-zcell2).lt.eps1) zrmid = zrmid-eps2
      if ((zrmid.ge.zbegproc(rank+1)).and.
     $  (zrmid.lt.zendproc(rank+1)).and.(yrmid.ge.ybegproc(rank+1))
     $  .and.(yrmid.lt.yendproc(rank+1))
     $  .and.(xrmid.ge.xbegproc(rank+1))
     $  .and.(xrmid.lt.xendproc(rank+1))) then
        docompute = .true.
      end if
      return
      end
c
c
c
c     subroutine midpointgroup : routine that says whether an interaction between a
c     group of particules
c     has to be computed within the current domain or not (dd midpoint method, Newton's
c     3rd law)
c
c
      subroutine midpointgroup(pos,nb,docompute)
      use domdec
      implicit none
      integer nb,i
      real*8 pos(3,nb),posdir(3,nb)
      real*8 xrmid,yrmid,zrmid
      logical docompute
c
      docompute = .false.

      do i = 1, nb
        call image(pos(1,i),pos(2,i),pos(3,i))
      end do
      do i = 2, nb
        posdir(1,i) = pos(1,i) - pos(1,1)
        posdir(2,i) = pos(2,i) - pos(2,1)
        posdir(3,i) = pos(3,i) - pos(3,1)
        call image(posdir(1,i),posdir(2,i),posdir(3,i))
      end do
c     
c     definition of the middle point between all the atoms
c
      xrmid = pos(1,1) 
      yrmid = pos(2,1) 
      zrmid = pos(3,1) 
      do i = 2, nb
        xrmid = xrmid + posdir(1,i)/nb 
        yrmid = yrmid + posdir(2,i)/nb 
        zrmid = zrmid + posdir(3,i)/nb 
      end do
c
      call image(xrmid,yrmid,zrmid)
      if ((zrmid.ge.zbegproc(rank+1)).and.
     $  (zrmid.lt.zendproc(rank+1)).and.(yrmid.ge.ybegproc(rank+1))
     $  .and.(yrmid.lt.yendproc(rank+1))
     $  .and.(xrmid.ge.xbegproc(rank+1))
     $  .and.(xrmid.lt.xendproc(rank+1))) then
        docompute = .true.
      end if
      return
      end
c
c     subroutine ddpme3dnpt : rescale geomtry of the domains and recompute the related quantities
c     for communications
c
      subroutine ddpme3dnpt(scaleiso,istep)
      real*8 scaleiso
      integer istep
      real*8 scale(3)

      scale(:)=scaleiso
      call ddpme3dnptaniso(scale,istep)

      end subroutine ddpme3dnpt

      subroutine ddpme3dnptaniso(scale,istep)
      use cutoff
      use domdec
      use neigh
      use potent
      use mpi
      implicit none
      integer nprocloc,rankloc,commloc,iproc
      integer i,status(MPI_STATUS_SIZE),tag,ierr
      integer istep,modnl
      real*8 mbuf,vbuf,neigbuf,bigbuf
      real*8 mshortbuf,vshortbuf,bigshortbuf
      real*8 dist, scale(3)
      integer, allocatable :: bufcount(:),buffer(:,:)
      integer, allocatable :: reqsend(:),reqrec(:)
cc
      modnl = mod(istep,ineigup)
cc
      if (use_pmecore) then
        nprocloc = ndir
        commloc  = comm_dir
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
      do iproc = 1, nprocloc
        xbegproc(iproc) = scale(1)*xbegproc(iproc)
        xendproc(iproc) = scale(1)*xendproc(iproc)
        ybegproc(iproc) = scale(2)*ybegproc(iproc)
        yendproc(iproc) = scale(2)*yendproc(iproc)
        zbegproc(iproc) = scale(3)*zbegproc(iproc)
        zendproc(iproc) = scale(3)*zendproc(iproc)
      end do
      if (modnl.ne.0) return
c
      allocate (reqsend(nprocloc))
      allocate (reqrec(nprocloc))
      allocate (bufcount(nprocloc))
      allocate (buffer(nprocloc,nprocloc))
c
c     get the processes to receive data from
c
      p_send1 = 0
      p_recep1 = 0
      p_send2 = 0
      p_recep2 = 0
      pbig_send = 0
      pbig_recep = 0
      n_recep1 = 0
      n_recep2 = 0
      nneig_recep = 0
      nbig_recep = 0
c
      p_sendshort1 = 0
      p_recepshort1 = 0
      p_sendshort2 = 0
      p_recepshort2 = 0
      pbigshort_send = 0
      pbigshort_recep = 0
      n_recepshort1 = 0
      n_recepshort2 = 0
      nbigshort_recep = 0

      if ((use_pmecore).and.(rank.gt.ndir-1)) then
        goto 80
      end if
c
c     choose cutoff depending on electrostatic interaction
c
      if ((use_mpole).or.(use_polar)) then
        mbuf = sqrt(mbuf2)
        mshortbuf = sqrt(mshortbuf2)
      else if (use_charge) then
        mbuf = sqrt(cbuf2)
        mshortbuf = sqrt(cshortbuf2)
      else
        mbuf = 0.0d0
        mshortbuf = 0d0
      end if
c
      vbuf = sqrt(vbuf2)+2.0d0
      vshortbuf = sqrt(vshortbuf2)+2.0d0
      neigbuf = lbuffer
c
c     get maximum cutoff value
c
      bigbuf = max(mbuf+2d0,vbuf,ddcut)
      bigshortbuf = max(mshortbuf+2d0,vshortbuf,ddcut)
c
      do iproc = 0, nprocloc-1
          if (iproc.eq.rank) goto 70
          call distproc(rank,iproc,dist,.true.)
          if (dist.le.(mbuf/2)) then
            n_recep1 = n_recep1+1
            p_recep1(n_recep1) = iproc
          end if
          if (dist.le.(vbuf/2)) then
            n_recep2 = n_recep2+1
            p_recep2(n_recep2) = iproc
          end if
          if (dist.le.(bigbuf/2)) then
            nbig_recep = nbig_recep+1
            pbig_recep(nbig_recep) = iproc
          end if
c
c    also deal with short range non bonded interactions (respa-1 like integrators)
c
          if (dist.le.(mshortbuf/2)) then
            n_recepshort1 = n_recepshort1+1
            p_recepshort1(n_recepshort1) = iproc
          end if
          if (dist.le.(vshortbuf/2)) then
            n_recepshort2 = n_recepshort2+1
            p_recepshort2(n_recepshort2) = iproc
          end if
          if (dist.le.(bigshortbuf/2)) then
            nbigshort_recep = nbigshort_recep+1
            pbigshort_recep(nbigshort_recep) = iproc
          end if
c
          if (dist.le.neigbuf) then
            nneig_recep = nneig_recep+1
            pneig_recep(nneig_recep) = iproc
          end if
 70     continue
      end do
c
      n_send1 = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(n_recep1,1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(p_recep1,n_recep1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              n_send1 = n_send1 + 1
              p_send1(n_send1) = iproc
            end if
          end do
        end if
      end do
c
c     get the processes to send data to
c
      n_send2 = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(n_recep2,1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(p_recep2,n_recep2,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              n_send2 = n_send2 + 1
              p_send2(n_send2) = iproc
            end if
          end do
        end if
      end do
c
      nneig_send = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(nneig_recep,1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(pneig_recep,nneig_recep,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              nneig_send = nneig_send + 1
              pneig_send(nneig_send) = iproc
            end if
          end do
        end if
      end do
c
      nbig_send = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(nbig_recep,1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(pbig_recep,nbig_recep,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              nbig_send = nbig_send + 1
              pbig_send(nbig_send) = iproc
            end if
          end do
        end if
      end do
c
c
c
      n_sendshort1 = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(n_recepshort1,1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(p_recepshort1,n_recepshort1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              n_sendshort1 = n_sendshort1 + 1
              p_sendshort1(n_sendshort1) = iproc
            end if
          end do
        end if
      end do
c
c     get the processes to send data to
c
      n_sendshort2 = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(n_recepshort2,1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(p_recepshort2,n_recepshort2,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              n_sendshort2 = n_sendshort2 + 1
              p_sendshort2(n_sendshort2) = iproc
            end if
          end do
        end if
      end do
c
      nbigshort_send = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(nbigshort_recep,1,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(pbigshort_recep,nbigshort_recep,MPI_INT,iproc,
     $     tag,commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              nbigshort_send = nbigshort_send + 1
              pbigshort_send(nbigshort_send) = iproc
            end if
          end do
        end if
      end do
cc
c
 80   call orderbuffer(.false.)
      call orderbufferrec

c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (bufcount)
      return
      end
