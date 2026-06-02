!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine drivermpi  --  driver for MPI related quantities ##
!     ##                            (3d spatial decomposition)        ##
!     ##                                                              ##
!     ##################################################################
!
!c
!
!> @brief 
!> driver to initialize 3d domain decomposition
!> @param no params
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
!
   if (deb_Path) write(iout,*), 'drivermpi '
!
!
   ndir = nproc - nrec
!
!     MPI : get the atoms repartition over the processes
!
!     deallocate global arrays
!
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
   if (allocated(cbegproc)) deallocate(cbegproc)
   if (allocated(cendproc)) deallocate(cendproc)
   if (allocated(bbegproc)) deallocate(bbegproc)
   if (allocated(bendproc)) deallocate(bendproc)
   if (allocated(abegproc)) deallocate(abegproc)
   if (allocated(aendproc)) deallocate(aendproc)
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
         if (rank.eq.0) write(iout,*)&
         &'no cores assigned to compute reciprocal space contribution'
         call fatal
      end if
      if (nproc-nrec.lt.1) then
         if (rank.eq.0) write(iout,*)&
         &'not enough cores to compute reciprocal space contribution'
         call fatal
      end if
!
!    allocate global arrays
!
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
      allocate (cbegproc(nproc))
      allocate (cendproc(nproc))
      allocate (bbegproc(nproc))
      allocate (bendproc(nproc))
      allocate (abegproc(nproc))
      allocate (aendproc(nproc))
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
!
!
!     build the two mpi groups for the computation of the multipole interactions
!     (direct and reciprocal space) and associated communicators
!
      allocate (direct_rank(0:nproc-1))
      call MPI_Comm_group(COMM_TINKER, total_group, ierr)
      do iproc = 0, ndir - 1
         direct_rank(iproc) = iproc
      end do
      call MPI_Group_incl(total_group,ndir,direct_rank,&
      &direct_group,ierr)
      call MPI_Group_excl(total_group,ndir,direct_rank,&
      &rec_group,ierr)
      call MPI_Comm_create(COMM_TINKER,direct_group,comm_dir,ierr)
      call MPI_Comm_create(COMM_TINKER,rec_group,comm_rec,ierr)
      if (rank.le.ndir-1) then
         call MPI_COMM_RANK(comm_dir,rank_bis,ierr)
         CALL MPI_Comm_split_type(comm_dir, MPI_COMM_TYPE_SHARED, 0,&
         &MPI_INFO_NULL, hostcomm,ierr)
         CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
      else
         call MPI_COMM_RANK(comm_rec,rank_bis,ierr)
         CALL MPI_Comm_split_type(comm_rec, MPI_COMM_TYPE_SHARED, 0,&
         &MPI_INFO_NULL, hostcomm,ierr)
         CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
      end if
      deallocate (direct_rank)
!
!       call the dd load balancing routine
!
      call ddpme3d
   else
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
      allocate (cbegproc(nproc))
      allocate (cendproc(nproc))
      allocate (bbegproc(nproc))
      allocate (bendproc(nproc))
      allocate (abegproc(nproc))
      allocate (aendproc(nproc))
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
      allocate (bufbegrec(nproc))
      allocate (bufbegpole(nproc))
      allocate (bufbeg(nproc))
      call ddpme3d
   end if
   return
end
!
!     subroutine allocstep: deallocate arrays and reallocate them with proper size
!     (memory distribution)
!
!> @brief 
!> subroutine allocstep: deallocate arrays and reallocate them with proper size
!> (memory distribution)
!> @param no params
subroutine allocstep
   use deriv
   use domdec
   use inform
   use iounit
   implicit none
!
   if (deb_Path) write(iout,*), 'allocstep '
!
!
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
!
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
!
   return
end
!
!     subroutine allocstepsrespa: deallocate arrays and reallocate them with proper size
!     (memory distribution)
!
!> @brief 
!> subroutine allocstepsrespa: deallocate arrays and reallocate them with proper size
!> (memory distribution)
!> @param no params
subroutine allocsteprespa(fast)
   use deriv
   use domdec
   use inform
   use iounit
   implicit none
   logical fast
!
   if (deb_Path) write(iout,*), 'allocsteprespa '
!
!
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
!
   if (.not.(fast)) then
      decrec = 0d0
      demrec = 0d0
      deprec = 0d0
      dedsprec = 0d0
   end if
!
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
!
   return
end
!
!     subroutine ddnumber : get the number of subdivision along each axis for
!     3d spatial decomposition
!
!> @brief 
!> get the number of subdivision along each axis for
!> 3d spatial decomposition
!> @param no params
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
10 format('Nx = ',I5,2x,'Ny = ',I5,2x,'Nz = ',I5,2x)
11 format('User defined 3D decompostion ','Nx = ',I5,2x,&
   &'Ny = ',I5,2x,'Nz = ',I5,2x)
12 format('User defined 3D decomposition not compatible with number',&
   &' of cores ','Nx*Ny*Nz = ',I5,2x,'number of procs = ',I5,2x)
13 format('The program will impose the 3D decomposition')
!
!
!     check for keywords containing domain decomposition parameters
!
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:240)
      if (keyword(1:9) .eq. 'DECOMP3D ') then
         read (string,*,err=20,end=20)  nxdd,nydd,nzdd
20       continue
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
!
   allocate (d(num))
   d = 0
!
!    Sort the axis by size
!
   list1(1) = xbox
   list1(2) = ybox
   list1(3) = zbox
   call sort2(3,list1,key)
!
!     get prime decomposition of number
!
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
!
      n1 = floor(num**(1.0/3.0))
      do i = 0, n1-2
         res = mod(num,n1-i)
         if (res.eq.0) goto 30
      end do
30    continue
      n1 = n1-i
      n2 = floor((num/n1)**(1.0/2.0))
      do i = 0, n2-2
         res = mod(num/n1,n2-i)
         if (res.eq.0) goto 40
      end do
40    continue
      n2 = n2 - i
      n3 = num/(n1*n2)
      list2(1) = n1
      list2(2) = n2
      list2(3) = n3
      call sort(3,list2)
!
      if (list2(1).eq.1) then
!
!      try dividing first by a smaller number
!
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
!
70    if (key(3).eq.1) then
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
!
!
!> @brief 
!> 3d domain decomposition
!> @param no params
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
   integer nprocloc,rankloc,iproc
   real*8 xr,yr,zr
   real*8 mbuf,vbuf,neigbuf,bigbuf
   real*8 mshortbuf,vshortbuf,bigshortbuf
   real*8 eps1,eps2
   real*8, allocatable :: abegproctemp(:),bbegproctemp(:)
   real*8, allocatable :: cbegproctemp(:)
   real*8, allocatable :: aendproctemp(:),bendproctemp(:)
   real*8, allocatable :: cendproctemp(:)
   real*8 ar,br,cr
   real*8 amin,amax
   real*8 bmin,bmax
   real*8 cmin,cmax
   real*8 nnorm
   real*8 a(3),b(3),c(3),nn(3),o(3),pp(3)
   integer na_cell,nb_cell,nc_cell
   integer n_ar,n_br,n_cr
   integer n_abeg,n_bbeg,n_cbeg
   integer na_buf,nb_buf,nc_buf
   integer p,q,r
   integer proc,tempa,tempb,tempc
!
1000 format(' Warning, less than 10 atoms on process number',I6,x,&
   &' number of cores may be too high compared to the number of '&
   &'atoms')
!
   if (deb_Path) write(iout,*), 'ddpme3d '
!
!
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
!
   nneig_recep = 0
   pneig_recep = 0
   nbig_recep = 0
   pbig_recep = 0
   nbigshort_recep = 0
   pbigshort_recep = 0
!
   allocate (abegproctemp(nproc))
   allocate (aendproctemp(nproc))
   allocate (bbegproctemp(nproc))
   allocate (bendproctemp(nproc))
   allocate (cbegproctemp(nproc))
   allocate (cendproctemp(nproc))
   abegproc = 0d0
   bbegproc = 0d0
   cbegproc = 0d0
   abegproctemp = 0d0
   bbegproctemp = 0d0
   cbegproctemp = 0d0
   aendproc = 0d0
   bendproc = 0d0
   cendproc = 0d0
   aendproctemp = 0d0
   bendproctemp = 0d0
   cendproctemp = 0d0
   repart = -1
   domlen = 0
   glob = 0
   loc = 0
!
!     get the number of subdivision along each axis for dd
!
   call ddnumber(nprocloc,0)
!
   na_box = xbox/nxdd
   nb_box = ybox/nydd
   nc_box = zbox/nzdd
   amin = -xbox2
   amax = xbox2
   bmin = -ybox2
   bmax = ybox2
   cmin = -zbox2
   cmax = zbox2
   do i = 0, nxdd-1
      abegproctemp(i+1) = -xbox2 + i*na_box
      aendproctemp(i+1) = -xbox2 + (i+1)*na_box
   end do
   do i = 0, nydd-1
      bbegproctemp(i+1) = -ybox2 + i*nb_box
      bendproctemp(i+1) = -ybox2 + (i+1)*nb_box
   end do
   do i = 0, nzdd-1
      cbegproctemp(i+1) = -zbox2 + i*nc_box
      cendproctemp(i+1) = -zbox2 + (i+1)*nc_box
   end do
!
!     assign processes
!
   do k = 1, nzdd
      do j = 1, nydd
         do i = 1, nxdd
            iproc = (k-1)*nydd*nxdd+(j-1)*nxdd+i
            abegproc(iproc) = abegproctemp(i)
            aendproc(iproc) = aendproctemp(i)
            bbegproc(iproc) = bbegproctemp(j)
            bendproc(iproc) = bendproctemp(j)
            cbegproc(iproc) = cbegproctemp(k)
            cendproc(iproc) = cendproctemp(k)
         end do
      end do
   end do
!
!     count number of particules per domain
!
   domlen = 0
   do i = 1, n
      xr = x(i)
      yr = y(i)
      zr = z(i)
!
!       get fractional coordinates
!
      call ctfvec(xr,yr,zr,ar,br,cr)
!
!       put particle in unit cell
!
      call imagefrac(ar,br,cr)
      if (abs(ar-xbox2).lt.eps1) ar = ar-eps2
      if (abs(br-ybox2).lt.eps1) br = br-eps2
      if (abs(cr-zbox2).lt.eps1) cr = cr-eps2
      n_ar = int((ar+xbox2)/na_box)
      n_br = int((br+ybox2)/nb_box)
      n_cr = int((cr+zbox2)/nc_box)
      iproc = n_ar + n_br*nxdd + n_cr*nydd*nxdd
      repart(i) = iproc
      domlen(repart(i)+1) = domlen(repart(i)+1) + 1
   end do
!
! get distance between consecutive planes (parallel to unit cell) corresponding
! to cells
!
a(1) = lvec(1,1)
a(2) = lvec(1,2)
a(3) = lvec(1,3)
b(1) = lvec(2,1)
b(2) = lvec(2,2)
b(3) = lvec(2,3)
c(1) = lvec(3,1)
c(2) = lvec(3,2)
c(3) = lvec(3,3)
!lena = abs(amax-amin)
nn(1) = b(2)*c(3)-b(3)*c(2)
nn(2) = b(3)*c(1)-b(1)*c(3)
nn(3) = b(1)*c(2)-b(2)*c(1)
nnorm = sqrt(nn(1)**2+nn(2)**2+nn(3)**2)
nn(1) = nn(1)/nnorm
nn(2) = nn(2)/nnorm
nn(3) = nn(3)/nnorm
o(1) = 1d0
o(2) = 0d0
o(3) = 0d0
call ftcvec(o(1),o(2),o(3),pp(1),pp(2),pp(3))
scala = abs(nn(1)*pp(1)+nn(2)*pp(2)+nn(3)*pp(3))
na_cell = na_box
nn(1) = a(2)*c(3)-a(3)*c(2)
nn(2) = a(3)*c(1)-a(1)*c(3)
nn(3) = a(1)*c(2)-a(2)*c(1)
nnorm = sqrt(nn(1)**2+nn(2)**2+nn(3)**2)
nn(1) = nn(1)/nnorm
nn(2) = nn(2)/nnorm
nn(3) = nn(3)/nnorm
o(1) = 0d0
o(2) = 1d0
o(3) = 0d0
call ftcvec(o(1),o(2),o(3),pp(1),pp(2),pp(3))
scalb = abs(nn(1)*pp(1)+nn(2)*pp(2)+nn(3)*pp(3))
nb_cell = nb_box
nn(1) = a(2)*b(3)-a(3)*b(2)
nn(2) = a(3)*b(1)-a(1)*b(3)
nn(3) = a(1)*b(2)-a(2)*b(1)
nnorm = sqrt(nn(1)**2+nn(2)**2+nn(3)**2)
nn(1) = nn(1)/nnorm
nn(2) = nn(2)/nnorm
nn(3) = nn(3)/nnorm
o(1) = 0d0
o(2) = 0d0
o(3) = 1d0
call ftcvec(o(1),o(2),o(3),pp(1),pp(2),pp(3))
scalc = abs(nn(1)*pp(1)+nn(2)*pp(2)+nn(3)*pp(3))
nc_cell = nc_box
!
!     choose cutoff depending on electrostatic interaction
!
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
!
   vbuf = sqrt(vbuf2)+2.0d0
!
!  take some margin because of torques
!
   mbuf = sqrt(mbuf2)+2.0d0
   vshortbuf = sqrt(vshortbuf2)+2.0d0
   mshortbuf = sqrt(mshortbuf2)+2.0d0
   neigbuf = 2.0d0!lbuffer
!
!     get maximum cutoff value
!
   bigbuf = max(mbuf,vbuf,ddcut)
   bigshortbuf = max(mshortbuf,vshortbuf,ddcut)
!
! how many cells do we need as neighbours in each direction and for each cutoff
!
! coordinates of "cells" in fractional space
!
!   get neighbouring processes
!
   nbig_recep = 0
   pbig_recep = 0
   nbig_send = 0
   pbig_send = 0

   n_recep2 = 0
   p_recep2 = 0
   n_send2 = 0
   p_send2 = 0

   n_recep1 = 0
   p_recep1 = 0
   n_send1 = 0
   p_send1 = 0

   nneig_recep = 0
   pneig_recep = 0
   nneig_send = 0
   pneig_send = 0

   n_recepshort1 = 0
   p_recepshort1 = 0
   n_sendshort1 = 0
   p_sendshort1 = 0

   n_recepshort2 = 0
   p_recepshort2 = 0
   n_sendshort2 = 0
   p_sendshort2 = 0

   nbigshort_recep = 0
   pbigshort_recep = 0
   nbigshort_send = 0
   pbigshort_send = 0
  
   if ((use_pmecore).and.(rank.le.ndir-1)) then
      nloc = domlen(rankloc+1)
!
!       check for low atom number in each atom for 'extreme' parallelism
!
      if ((nloc.lt.10).and.(nproc.gt.32)) then
         write(iout,1000) rank
      end if
   else if ((use_pmecore).and.(rank.gt.ndir-1)) then
      nloc = 0
   else
      nloc = domlen(rankloc+1)
!
!       check for low atom number in each atom for 'extreme' parallelism
!
      if ((nloc.lt.10).and.(nproc.gt.32)) then
         write(iout,1000) rank
      end if
   end if

   if ((use_pmecore).and.(rank.gt.ndir-1)) then
      goto 180
   end if

   n_abeg = int((0.5*(abegproc(rank+1)+aendproc(rank+1))+xbox2)/na_box)
   n_bbeg = int((0.5*(bbegproc(rank+1)+bendproc(rank+1))+ybox2)/nb_box)
   n_cbeg = int((0.5*(cbegproc(rank+1)+cendproc(rank+1))+zbox2)/nc_box)

!
!   neighboring processes: bigbuf
!
na_buf = ceiling(bigbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(bigbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(bigbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 10
         do i = 1, nbig_recep
           if (proc.eq.pbig_recep(i)) goto 10
         end do
         nbig_recep = nbig_recep+1
         pbig_recep(nbig_recep) = proc
 10      continue 
       end do
     end do
   end do
!
!   neighboring processes: mbuf
!
na_buf = ceiling(mbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(mbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(mbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 20
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 20
         do i = 1, n_recep1
           if (proc.eq.p_recep1(i)) goto 20
         end do
         n_recep1 = n_recep1+1
         p_recep1(n_recep1) = proc
 20      continue 
       end do
     end do
   end do
!
!   neighboring processes: vbuf
!
na_buf = ceiling(vbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(vbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(vbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 30
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 30
         do i = 1, n_recep2
           if (proc.eq.p_recep2(i)) goto 30
         end do
         n_recep2 = n_recep2+1
         p_recep2(n_recep2) = proc
 30      continue 
       end do
     end do
   end do
!
!   neighboring processes: mshortbuf
!
na_buf = ceiling(mshortbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(mshortbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(mshortbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 40
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 40
         do i = 1, n_recepshort1
           if (proc.eq.p_recepshort1(i)) goto 40
         end do
         n_recepshort1 = n_recepshort1+1
         p_recepshort1(n_recepshort1) = proc
 40      continue 
       end do
     end do
   end do
!
!   neighboring processes: vshortbuf
!
na_buf = ceiling(vshortbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(vshortbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(vshortbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 50
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 50
         do i = 1, n_recepshort2
           if (proc.eq.p_recepshort2(i)) goto 50
         end do
         n_recepshort2 = n_recepshort2+1
         p_recepshort2(n_recepshort2) = proc
 50      continue 
       end do
     end do
   end do
!
!   neighboring processes: bigshortbuf
!
na_buf = ceiling(bigshortbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(bigshortbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(bigshortbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 60
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 60
         do i = 1, nbigshort_recep
           if (proc.eq.pbigshort_recep(i)) goto 60
         end do
         nbigshort_recep = nbigshort_recep+1
         pbigshort_recep(nbigshort_recep) = proc
 60      continue 
       end do
     end do
   end do
!
!   neighboring processes: neigbuf
!
na_buf = ceiling(neigbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(neigbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(neigbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 70
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 70
         do i = 1, nneig_recep
           if (proc.eq.pneig_recep(i)) goto 70
         end do
         nneig_recep = nneig_recep+1
         pneig_recep(nneig_recep) = proc
 70      continue 
       end do
     end do
   end do

   n_send1 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recep1
            if (p_recep1(i).eq.iproc) then
               n_send1 = n_send1 + 1
               p_send1(n_send1) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   n_send2 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recep2
            if (p_recep2(i).eq.iproc) then
               n_send2 = n_send2 + 1
               p_send2(n_send2) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   nbig_send = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, nbig_recep
            if (pbig_recep(i).eq.iproc) then
               nbig_send = nbig_send + 1
               pbig_send(nbig_send) = iproc
            end if
         end do
      end if
   end do
!
   n_sendshort1 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recepshort1
            if (p_recepshort1(i).eq.iproc) then
               n_sendshort1 = n_sendshort1 + 1
               p_sendshort1(n_sendshort1) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   n_sendshort2 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recepshort2
            if (p_recepshort2(i).eq.iproc) then
               n_sendshort2 = n_sendshort2 + 1
               p_sendshort2(n_send2) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   nbigshort_send = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, nbigshort_recep
            if (pbigshort_recep(i).eq.iproc) then
               nbigshort_send = nbigshort_send + 1
               pbigshort_send(nbigshort_send) = iproc
            end if
         end do
      end if
   end do
!
   nneig_send = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, nneig_recep
            if (pneig_recep(i).eq.iproc) then
               nneig_send = nneig_send + 1
               pneig_send(nneig_send) = iproc
            end if
         end do
      end if
   end do
!if (rank.eq.0) then
!      write(*,*) 'nbigshort_send = ',nbigshort_send
!      write(*,*) 'nbig_send = ',nbig_send
!      write(*,*) 'n_sendshort1 = ',n_sendshort1
!      write(*,*) 'n_send1 = ',n_send1
!      write(*,*) 'n_sendshort2 = ',n_sendshort2
!      write(*,*) 'n_send2 = ',n_send2
!end if
!
180 call orderbuffer(.true.)
    
   deallocate (abegproctemp)
   deallocate (aendproctemp)
   deallocate (bbegproctemp)
   deallocate (bendproctemp)
   deallocate (cbegproctemp)
   deallocate (cendproctemp)
   return
end
!
!
!     subroutine halfcell : routine that says whether an interaction between two particules
!     has to be computed within the current domain or not (dd half cell method, Newton's
!     3rd law)
!
!
!> @brief 
!> routine that says whether an interaction between two particules
!> has to be computed within the current domain or not (dd half cell method, Newton's
!> 3rd law)
!> @param no params
subroutine halfcell(xi,yi,zi,xj,yj,zj,docompute)
   use bound
   implicit none
   real*8 xr,yr,zr
   real*8 xi,yi,zi
   real*8 xj,yj,zj
   logical docompute
!
   docompute = .false.
!
   xr = xi - xj
   yr = yi - yj
   zr = zi - zj
   if (use_bounds) call image(xr,yr,zr)
   if (xr.lt.0.0d0) then
      docompute = .true.
   else if ((xr.eq.0.0d0).and.(yr.lt.0.0d0)) then
      docompute = .true.
   else if ((xr.eq.0.0d0).and.(yr.eq.0.0d0).and.&
   &(zr.lt.0.0d0)) then
      docompute = .true.
   end if
   return
end
!
!     subroutine midpoint : routine that says whether an interaction between two particules
!     has to be computed within the current domain or not (dd midpoint method)
!
!
!> @brief 
!> routine that says whether an interaction between two particules
!> has to be computed within the current domain or not (dd midpoint method)
!> @param no params
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
!
   eps1 = 10d-10
   eps2 = 10d-8
   docompute = .false.
!
!      call image(xi,yi,zi)
!      call image(xk,yk,zk)
   xr = xi - xk
   yr = yi - yk
   zr = zi - zk
   call image(xr,yr,zr)
!
!     definition of the middle point between i and k atoms
!
   xrmid = xk + xr/2
   yrmid = yk + yr/2
   zrmid = zk + zr/2
   call image(xrmid,yrmid,zrmid)
   if (abs(xrmid-xcell2).lt.eps1) xrmid = xrmid-eps2
   if (abs(yrmid-ycell2).lt.eps1) yrmid = yrmid-eps2
   if (abs(zrmid-zcell2).lt.eps1) zrmid = zrmid-eps2
   if ((zrmid.ge.zbegproc(rank+1)).and.&
   &(zrmid.lt.zendproc(rank+1)).and.(yrmid.ge.ybegproc(rank+1))&
   &.and.(yrmid.lt.yendproc(rank+1))&
   &.and.(xrmid.ge.xbegproc(rank+1))&
   &.and.(xrmid.lt.xendproc(rank+1))) then
      docompute = .true.
   end if
   return
end
!
!     subroutine ddbuild: build lists of cells for a given cutoff (cell-list
!     method) and also nlocnl
!
!> @brief 
!> build lists of cells for a given cutoff (cell-list
!> method) and also nlocnl
!> build lists of neighbors for a given cutoff
!> @param[in] cut: cutoff 
!> @param[in] nrecep: number of neighboring procs 
!> @param[in] precep: list of neighboring procs 
!> @param[in] init: logical true if first call 
subroutine ddbuild(cut,nrecep,precep,init)
use atoms
use boxes
use cutoff
use domdec
use iounit
use neigh
use potent
use mpi
implicit none
integer i,j,k,l,iglob,iloc,count
integer p,q,r
integer :: nrecep
integer, dimension(nproc) :: precep
integer ncell,icell,jcell,ncelltemp
integer, allocatable :: lcell_neig(:),indcelltemp(:)
integer na,nb,nc
integer n_xr,n_yr,n_zr
integer n_amin,n_bmin,n_cmin
integer n_amax,n_bmax,n_cmax
integer tempx,tempy,tempz
integer imin,imax,jmin,jmax,kmin,kmax
integer itemp,jtemp,ktemp
real*8 lena_cell,lenb_cell,lenc_cell
real*8 eps1,eps2,eps3
real*8 cut
real*8 amin,amax,bmin,bmax,cmin,cmax
real*8 xr,yr,zr,ar,br,cr
real*8 lena
real*8 lenb
real*8 lenc
logical init,testdist


eps1 = 10d-10
eps2 = 10d-8
eps3 = 10d-5

!
amin = -xbox2
amax =  xbox2
bmin = -ybox2
bmax =  ybox2
cmin = -zbox2
cmax =  zbox2

lena = abs(amax-amin)
lenb = abs(bmax-bmin)
lenc = abs(cmax-cmin)

na = max(1,int(3d0*lena*scala/(cut)))
nb = max(1,int(3d0*lenb*scalb/(cut)))
nc = max(1,int(3d0*lenc*scalc/(cut)))
!write(*,*) 'rank = ',rank,'na = ',na,'nb = ',nb,'nc = ',nc
!write(*,*) 'rank = ',rank,'lena = ',lena,amin,amax
!write(*,*) 'rank = ',rank,'nbig_recep = ',nbig_recep

lena_cell = lena/na
lenb_cell = lenb/nb
lenc_cell = lenc/nc


ncell = na*nb*nc

!write(*,*) 'ncell = ',ncell
!write(*,*) 'na_cut = ',na_cut,'nb_cut = ',nb_cut,'nc_cut = ',nc_cut
!write(*,*) 'scala = ',scala,'scalb = ',scalb,'scalc = ',scalc
!end if
if (allocated(cell_len)) deallocate (cell_len)
allocate (cell_len(ncell))
cell_len = 0
if (allocated(repartcell)) deallocate (repartcell)
allocate (repartcell(n))
if (allocated(indcell)) deallocate (indcell)
allocate (indcell(n))
if (allocated(bufbegcell)) deallocate (bufbegcell)
allocate (bufbegcell(ncell))
allocate (indcelltemp(n))
indcelltemp = 0

if (init) then
if (allocated(lcell_neig)) deallocate(lcell_neig)
allocate(lcell_neig(ncell))
lcell_neig = 0
if (allocated(lcell_nl)) deallocate(lcell_nl)
allocate(lcell_nl(ncell))
lcell_nl = 0
if (allocated(neigcell)) deallocate (neigcell)
allocate (neigcell(400,ncell))
if (allocated(numneigcell)) deallocate (numneigcell)
allocate (numneigcell(ncell))



ar = abegproc(rank+1)
br = bbegproc(rank+1)
cr = cbegproc(rank+1)
n_amin = int((ar+xbox2)/lena_cell)
n_bmin = int((br+ybox2)/lenb_cell)
n_cmin = int((cr+zbox2)/lenc_cell)

ar = aendproc(rank+1)
br = bendproc(rank+1)
cr = cendproc(rank+1)
n_amax = int((ar+xbox2)/lena_cell)
n_bmax = int((br+ybox2)/lenb_cell)
n_cmax = int((cr+zbox2)/lenc_cell)

!
! ok if distance between cells below list-buffer
testdist = (cut/3d0.ge.lbuffer)
if (testdist.eq..false.) then
  if (rank.eq.0) write(iout,*) 'warning: margin taken for domain decomposition too small, cutoff divided by 3 has to be greater than list-buffer'
end if
!
imin = n_amin-2
imax = n_amax+2
jmin = n_bmin-2
jmax = n_bmax+2
kmin = n_cmin-2
kmax = n_cmax+2

ncell_nl = 0
lcell_nl = 0

do k = kmin,kmax
  do j = jmin,jmax
    do i = imin,imax
      itemp = modulo(i,na)
      jtemp = modulo(j,nb)
      ktemp = modulo(k,nc)
      icell = ktemp*nb*na+jtemp*na+itemp+1
      do l = 1, ncell_nl
        if (icell.eq.lcell_nl(l)) goto 20
      end do
      ncell_nl = ncell_nl+1
      lcell_nl(ncell_nl) = icell
 20     continue
    end do
  end do
end do
!write(*,*) 'ncell_nl = ',ncell_nl

!
!


do k = 0, nc-1
  do j = 0, nb-1
    do i = 0, na-1
  ncelltemp = 0
  lcell_neig = 0
  icell = k*nb*na+j*na+i+1
  do p = -3,3
    do q = -3,3
       do r = -3,3
       if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
!       if (abs(p)+abs(q)+abs(r).gt.6) goto 10
       tempx = modulo(i+p,na) 
       tempy = modulo(j+q,nb) 
       tempz = modulo(k+r,nc) 
       jcell = tempx + tempy*na + tempz*nb*na + 1
       do l = 1, ncell
         if (jcell.eq.lcell_neig(l)) goto 10
         if (jcell.eq.icell) goto 10
       end do
       ncelltemp = ncelltemp + 1
       lcell_neig(ncelltemp) = jcell
       neigcell(ncelltemp,icell) = jcell
10 continue
      end do
    end do
  end do
  numneigcell(icell) = ncelltemp
  end do
 end do
end do

end if


repartcell = -1
do i = 1, nbloc
  iglob = glob(i)
  xr = x(iglob)
  yr = y(iglob)
  zr = z(iglob)
!
!       get fractional coordinates
!
  call ctfvec(xr,yr,zr,ar,br,cr)
!
!       put particle in unit cell
!
  call imagefrac(ar,br,cr)

  if (abs(ar-xbox).lt.eps1) ar = ar-eps2
  if (abs(br-ybox).lt.eps1) br = br-eps2
  if (abs(cr-zbox).lt.eps1) cr = cr-eps2

  n_xr = int((ar+xbox2)/lena_cell)
  n_yr = int((br+ybox2)/lenb_cell)
  n_zr = int((cr+zbox2)/lenc_cell)
  icell = n_xr + n_yr*na + n_zr*nb*na+1
  repartcell(iglob) = icell
  cell_len(icell) = cell_len(icell) + 1
  indcelltemp(iglob) = cell_len(icell)
end do
!
bufbegcell = 0
bufbegcell(1) = 1
count = cell_len(1)
do icell = 2, ncell
   if (cell_len(icell).ne.0) then
      bufbegcell(icell) = count + 1
   else
      bufbegcell(icell) = 1
   end if
   count = count + cell_len(icell)
end do

do i = 1, nbloc
   iglob = glob(i)
   icell = repartcell(iglob)
   iloc  = bufbegcell(icell) + indcelltemp(iglob) - 1
   indcell(iloc) = iglob
end do
!
!  filter out "nlocnl" atoms
!
if (.not.allocated(ineignl)) allocate (ineignl(n))
ineignl = 0

nlocnl = 0
do icell = 1, ncell_nl 
  jcell = lcell_nl(icell)
  do i = 1, cell_len(jcell)
    iloc = bufbegcell(jcell) + i - 1
    iglob = indcell(iloc)
    nlocnl = nlocnl + 1
    ineignl(nlocnl) = iglob
  end do
end do
!write(*,*) 'rank = ',rank,'nlocnl = ',nlocnl
!
deallocate (indcelltemp)

return
end
!
!     subroutine midpointimage : routine that says whether an interaction between two particules
!     has to be computed within the current domain or not (dd midpoint method), also returns
!     minimum image of the distance vector
!
!
!> @brief 
!> routine that says whether an interaction between two particules
!> has to be computed within the current domain or not (dd midpoint method), also returns
!> minimum image of the distance vector
!> @param no params
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
!
   docompute = .false.
!
   xr = xi - xk
   yr = yi - yk
   zr = zi - zk
   call image(xr,yr,zr)
!
!     definition of the middle point between i and k atoms
!
   xrmid = xk + xr/2
   yrmid = yk + yr/2
   zrmid = zk + zr/2
   call image(xrmid,yrmid,zrmid)
   if (abs(xrmid-xcell2).lt.eps1) xrmid = xrmid-eps2
   if (abs(yrmid-ycell2).lt.eps1) yrmid = yrmid-eps2
   if (abs(zrmid-zcell2).lt.eps1) zrmid = zrmid-eps2
   if ((zrmid.ge.zbegproc(rank+1)).and.&
   &(zrmid.lt.zendproc(rank+1)).and.(yrmid.ge.ybegproc(rank+1))&
   &.and.(yrmid.lt.yendproc(rank+1))&
   &.and.(xrmid.ge.xbegproc(rank+1))&
   &.and.(xrmid.lt.xendproc(rank+1))) then
      docompute = .true.
   end if
   return
end
!
!> @brief 
!> routine that says whether an interaction between two particules
!> has to be computed within the current domain or not (dd midpoint method), also returns
!> minimum image of the distance vector
!> @param no params
   subroutine midpointimagetriclinic(xi,yi,zi,xk,yk,zk,xr,yr,zr,docompute)
   use boxes
   use cell
   use domdec
   implicit none
   real*8 xi,yi,zi
   real*8 xk,yk,zk
   real*8 xr,yr,zr
   real*8 ai,bi,ci
   real*8 ak,bk,ck
   real*8 ar,br,cr
   real*8 armid,brmid,crmid
   real*8 eps1,eps2
   logical docompute
   integer n_ar,n_br,n_cr,iproc
   eps1 = 10d-10
   eps2 = 10d-8
!
   docompute = .false.
   call ctfvec(xk,yk,zk,ak,bk,ck)
   call ctfvec(xi,yi,zi,ai,bi,ci)
!
   ar = ai - ak 
   br = bi - bk 
   cr = ci - ck 
   call imagefrac(ar,br,cr)
   call ftcvec(ar,br,cr,xr,yr,zr)
!   call image(xr,yr,zr)
!  
!  definition of the middle point between i and k atoms
!
   armid = ak + ar/2
   brmid = bk + br/2
   crmid = ck + cr/2
   call imagefrac(armid,brmid,crmid)
   if (abs(armid-xcell2).lt.eps1) armid = armid-eps2
   if (abs(brmid-ycell2).lt.eps1) brmid = brmid-eps2
   if (abs(crmid-zcell2).lt.eps1) crmid = crmid-eps2
   n_ar = int((armid+xbox2)/na_box)
   n_br = int((brmid+ybox2)/nb_box)
   n_cr = int((crmid+zbox2)/nc_box)
   iproc = n_ar + n_br*nxdd + n_cr*nydd*nxdd
   if (rank.eq.iproc) docompute = .true.
!
   return
   end
!
!
!
!     subroutine midpointgroup : routine that says whether an interaction between a
!     group of particules
!     has to be computed within the current domain or not (dd midpoint method, Newton's
!     3rd law)
!
!
!> @brief 
!> routine that says whether an interaction between a
!> group of particules
!> has to be computed within the current domain or not (dd midpoint method, Newton's
!> 3rd law)
!> @param no params
subroutine midpointgroup(pos,nb,docompute)
   use domdec
   implicit none
   integer nb,i
   real*8 pos(3,nb),posdir(3,nb)
   real*8 xrmid,yrmid,zrmid
   logical docompute
!
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
!
!     definition of the middle point between all the atoms
!
   xrmid = pos(1,1)
   yrmid = pos(2,1)
   zrmid = pos(3,1)
   do i = 2, nb
      xrmid = xrmid + posdir(1,i)/nb
      yrmid = yrmid + posdir(2,i)/nb
      zrmid = zrmid + posdir(3,i)/nb
   end do
!
   call image(xrmid,yrmid,zrmid)
   if ((zrmid.ge.zbegproc(rank+1)).and.&
   &(zrmid.lt.zendproc(rank+1)).and.(yrmid.ge.ybegproc(rank+1))&
   &.and.(yrmid.lt.yendproc(rank+1))&
   &.and.(xrmid.ge.xbegproc(rank+1))&
   &.and.(xrmid.lt.xendproc(rank+1))) then
      docompute = .true.
   end if
   return
end

!> @brief 
!> driver for update of domain decomposition after rescaling of volume due to
!> barostat
!> @param[in] scaleiso: scaling factor of every direction
!> @param[in] istep: timestep number
subroutine ddpme3dnpt(scaleiso,istep)
   real*8 scaleiso
   integer istep
   real*8 scale(3)

   scale(:)=scaleiso
   call ddpme3dnptaniso(scale,istep)

end subroutine ddpme3dnpt

!> @brief 
!> update of domain decomposition after rescaling of volume due to
!> barostat
!> @param[in] scale: scaling factor of every direction
!> @param[in] istep: timestep number
subroutine ddpme3dnptaniso(scale,istep)
   use boxes
   use cutoff
   use domdec
   use neigh
   use potent
   use mpi
   implicit none
   integer nprocloc,rankloc,commloc,iproc
   integer i
   integer istep,modnl
   real*8 mbuf,vbuf,neigbuf,bigbuf
   real*8 mshortbuf,vshortbuf,bigshortbuf
   real*8 scale(3)
   integer n_abeg,n_bbeg,n_cbeg
   integer na_buf,nb_buf,nc_buf
   integer p,q,r
   integer proc,tempa,tempb,tempc
   integer, allocatable :: bufcount(:),buffer(:,:)
   integer, allocatable :: reqsend(:),reqrec(:)
!c
   modnl = mod(istep,ineigup)
!c
   if (use_pmecore) then
      nprocloc = ndir
      commloc  = comm_dir
      rankloc  = rank_bis
   else
      nprocloc = nproc
      commloc  = COMM_TINKER
      rankloc  = rank
   end if
!
   do iproc = 1, nprocloc
      abegproc(iproc) = scale(1)*abegproc(iproc)
      aendproc(iproc) = scale(1)*aendproc(iproc)
      bbegproc(iproc) = scale(2)*bbegproc(iproc)
      bendproc(iproc) = scale(2)*bendproc(iproc)
      cbegproc(iproc) = scale(3)*cbegproc(iproc)
      cendproc(iproc) = scale(3)*cendproc(iproc)
   end do
   if (modnl.ne.0) return
!
   allocate (reqsend(nprocloc))
   allocate (reqrec(nprocloc))
   allocate (bufcount(nprocloc))
   allocate (buffer(nprocloc,nprocloc))
!
!     get the processes to receive data from
!
   nbig_recep = 0
   pbig_recep = 0
   nbig_send = 0
   pbig_send = 0

   n_recep2 = 0
   p_recep2 = 0
   n_send2 = 0
   p_send2 = 0

   n_recep1 = 0
   p_recep1 = 0
   n_send1 = 0
   p_send1 = 0

   nneig_recep = 0
   pneig_recep = 0
   nneig_send = 0
   pneig_send = 0

   n_recepshort1 = 0
   p_recepshort1 = 0
   n_sendshort1 = 0
   p_sendshort1 = 0

   n_recepshort2 = 0
   p_recepshort2 = 0
   n_sendshort2 = 0
   p_sendshort2 = 0

   nbigshort_recep = 0
   pbigshort_recep = 0
   nbigshort_send = 0
   pbigshort_send = 0

   if ((use_pmecore).and.(rank.gt.ndir-1)) then
      goto 80
   end if
!
!     choose cutoff depending on electrostatic interaction
!
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
!
   vbuf = sqrt(vbuf2)+2.0d0
!
!  take some margin because of torques
!
   mbuf = sqrt(mbuf2)+2.0d0
   vshortbuf = sqrt(vshortbuf2)+2.0d0
   mshortbuf = sqrt(mshortbuf2)+2.0d0
   neigbuf = 2.0d0!lbuffer
!
!     get maximum cutoff value
!
   bigbuf = max(mbuf,vbuf,ddcut)
   bigshortbuf = max(mshortbuf,vshortbuf,ddcut)
!!
   n_abeg = int((0.5*(abegproc(rank+1)+aendproc(rank+1))+xbox2)/na_box)
   n_bbeg = int((0.5*(bbegproc(rank+1)+bendproc(rank+1))+ybox2)/nb_box)
   n_cbeg = int((0.5*(cbegproc(rank+1)+cendproc(rank+1))+zbox2)/nc_box)

!
!   neighboring processes: bigbuf
!
na_buf = ceiling(bigbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(bigbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(bigbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 10
         do i = 1, nbig_recep
           if (proc.eq.pbig_recep(i)) goto 10
         end do
         nbig_recep = nbig_recep+1
         pbig_recep(nbig_recep) = proc
 10      continue 
       end do
     end do
   end do
!
!   neighboring processes: mbuf
!
na_buf = ceiling(mbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(mbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(mbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 20
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 20
         do i = 1, n_recep1
           if (proc.eq.p_recep1(i)) goto 20
         end do
         n_recep1 = n_recep1+1
         p_recep1(n_recep1) = proc
 20      continue 
       end do
     end do
   end do
!
!   neighboring processes: vbuf
!
na_buf = ceiling(vbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(vbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(vbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 30
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 30
         do i = 1, n_recep2
           if (proc.eq.p_recep2(i)) goto 30
         end do
         n_recep2 = n_recep2+1
         p_recep2(n_recep2) = proc
 30      continue 
       end do
     end do
   end do
!
!   neighboring processes: mshortbuf
!
na_buf = ceiling(mshortbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(mshortbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(mshortbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 40
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 40
         do i = 1, n_recepshort1
           if (proc.eq.p_recepshort1(i)) goto 40
         end do
         n_recepshort1 = n_recepshort1+1
         p_recepshort1(n_recepshort1) = proc
 40      continue 
       end do
     end do
   end do
!
!   neighboring processes: vshortbuf
!
na_buf = ceiling(vshortbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(vshortbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(vshortbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 50
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 50
         do i = 1, n_recepshort2
           if (proc.eq.p_recepshort2(i)) goto 50
         end do
         n_recepshort2 = n_recepshort2+1
         p_recepshort2(n_recepshort2) = proc
 50      continue 
       end do
     end do
   end do
!
!   neighboring processes: bigshortbuf
!
na_buf = ceiling(bigshortbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(bigshortbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(bigshortbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 60
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 60
         do i = 1, nbigshort_recep
           if (proc.eq.pbigshort_recep(i)) goto 60
         end do
         nbigshort_recep = nbigshort_recep+1
         pbigshort_recep(nbigshort_recep) = proc
 60      continue 
       end do
     end do
   end do
!
!   neighboring processes: neigbuf
!
na_buf = ceiling(neigbuf*nxdd/(2d0*scala*xbox))
nb_buf = ceiling(neigbuf*nydd/(2d0*scalb*ybox))
nc_buf = ceiling(neigbuf*nzdd/(2d0*scalc*zbox))

   do p = -na_buf,na_buf
     do q = -nb_buf,nb_buf
        do r = -nc_buf,nc_buf
         if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 70
         tempa = modulo(n_abeg+p,nxdd)
         tempb = modulo(n_bbeg+q,nydd)
         tempc = modulo(n_cbeg+r,nzdd)
         proc = tempa + tempb*nxdd + tempc*nydd*nxdd
         if (proc.eq.rank) goto 70
         do i = 1, nneig_recep
           if (proc.eq.pneig_recep(i)) goto 70
         end do
         nneig_recep = nneig_recep+1
         pneig_recep(nneig_recep) = proc
 70      continue 
       end do
     end do
   end do
   n_send1 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recep1
            if (p_recep1(i).eq.iproc) then
               n_send1 = n_send1 + 1
               p_send1(n_send1) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   n_send2 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recep2
            if (p_recep2(i).eq.iproc) then
               n_send2 = n_send2 + 1
               p_send2(n_send2) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   nbig_send = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, nbig_recep
            if (pbig_recep(i).eq.iproc) then
               nbig_send = nbig_send + 1
               pbig_send(nbig_send) = iproc
            end if
         end do
      end if
   end do
!
   n_sendshort1 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recepshort1
            if (p_recepshort1(i).eq.iproc) then
               n_sendshort1 = n_sendshort1 + 1
               p_sendshort1(n_sendshort1) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   n_sendshort2 = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, n_recepshort2
            if (p_recepshort2(i).eq.iproc) then
               n_sendshort2 = n_sendshort2 + 1
               p_sendshort2(n_send2) = iproc
            end if
         end do
      end if
   end do
!
!     get the processes to send data to
!
   nbigshort_send = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, nbigshort_recep
            if (pbigshort_recep(i).eq.iproc) then
               nbigshort_send = nbigshort_send + 1
               pbigshort_send(nbigshort_send) = iproc
            end if
         end do
      end if
   end do
!
   nneig_send = 0
!
   do iproc = 0, nprocloc-1
      if (iproc.ne.rankloc) then
         do i = 1, nneig_recep
            if (pneig_recep(i).eq.iproc) then
               nneig_send = nneig_send + 1
               pneig_send(nneig_send) = iproc
            end if
         end do
      end if
   end do
80 call orderbuffer(.false.)
   call orderbufferrec

   return
end
