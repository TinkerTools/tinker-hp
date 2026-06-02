!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"
module domdecstuff_inl
contains
#include "image.inc.f90"
end module

!     ##################################################################
!     ##                                                              ##
!     ##  subroutine drivermpi  --  driver for MPI related quantities ##
!     ##                            (3d spatial decomposition)        ##
!     ##                                                              ##
!     ##################################################################
!
!
subroutine drivermpi
   use atoms
   use domdec
   use iounit
   use potent
   use mpi
#ifdef _OPENACC
   use utilcu,only: copy_data_to_cuda_env
#endif
   use utilcomm
   use tinMemory
   implicit none
   integer iproc, ierr
   integer total_group, direct_group, rec_group
   integer, allocatable :: direct_rank(:)
!
   ndir = nproc - nrec
!
!     MPI : get the atoms repartition over the processes
!
!     allocate global arrays
!
   call prmem_request(glob,n)
   call prmem_request(loc,n)
   call prmem_request(globrec,n)
   call prmem_request(locrec,n)
   call prmem_request(repart,n)
   call prmem_request(domlen,nproc,config=mhostonly)
   call prmem_request(domlenpole,nproc,config=mhostonly)
   call prmem_request(domlenpolerec,nproc,config=mhostonly)
   call prmem_request(p_recepshort1,nproc,config=mhostonly)
   call prmem_request(p_recepshort2,nproc,config=mhostonly)
   call prmem_request(p_send1,nproc,config=mhostonly)
   call prmem_request(p_send2,nproc,config=mhostonly)
   call prmem_request(p_sendshort1,nproc,config=mhostonly)
   call prmem_request(p_sendshort2,nproc,config=mhostonly)
   call prmem_request(pneig_recep,nproc,config=mhostonly)
   call prmem_request(pneig_send,nproc ,config=mhostonly)
   call prmem_request(precdir_recep ,nproc,config=mhostonly)
   call prmem_request(precdir_send  ,nproc,config=mhostonly)
   call prmem_request(precdir_recep1,nproc,config=mhostonly)
   call prmem_request(precdir_send1 ,nproc,config=mhostonly)
   call prmem_request(precdir_recep2,nproc,config=mhostonly)
   call prmem_request(precdir_send2 ,nproc,config=mhostonly)
   call prmem_request(bufbegpole,nproc,config=mhostonly)
   call prmem_request(bufbeg    ,nproc,config=mhostonly)
   call prmem_request(bufbegrec ,nproc,config=mhostonly)
   call prmem_request(ptorque_recep     ,nproc,config=mhostonly)
   call prmem_request(ptorqueshort_recep,nproc,config=mhostonly)
   call prmem_request(ptorque_send      ,nproc,config=mhostonly)
   call prmem_request(ptorqueshort_send ,nproc,config=mhostonly)
   call prmem_request(pbig_recep     ,nproc,config=mhostonly)
   call prmem_request(pbigshort_recep,nproc,config=mhostonly)
   call prmem_request(pbig_send      ,nproc,config=mhostonly)
   call prmem_request(pbigshort_send ,nproc,config=mhostonly)
   call prmem_request(reqs_dird_a    ,nproc,config=mhostonly)
   call prmem_request(reqr_dird_a    ,nproc,config=mhostonly)
   call prmem_request(reqs_dirdir    ,nproc,config=mhostonly)
   call prmem_request(reqr_dirdir    ,nproc,config=mhostonly)
   call prmem_request(reqs_recdir    ,nproc,config=mhostonly)
   call prmem_request(reqr_recdir    ,nproc,config=mhostonly)
   call prmem_request(reqs_dirrec    ,nproc,config=mhostonly)
   call prmem_request(reqr_dirrec    ,nproc,config=mhostonly)
   call prmem_request(reqs_recdirsolv,nproc,config=mhostonly)
   call prmem_request(reqr_recdirsolv,nproc,config=mhostonly)
   call prmem_request(reqs_poleglob  ,nproc,config=mhostonly)
   call prmem_request(reqr_poleglob  ,nproc,config=mhostonly)
   call prmem_request(cbegproc, nproc)
   call prmem_request(cendproc, nproc)
   call prmem_request(bbegproc, nproc)
   call prmem_request(bendproc, nproc)
   call prmem_request(abegproc, nproc)
   call prmem_request(aendproc, nproc)
   ! TODO Remove declare directive
   if (allocated(zbegproc)) deallocate(zbegproc)
   if (allocated(zendproc)) deallocate(zendproc)
   if (allocated(ybegproc)) deallocate(ybegproc)
   if (allocated(yendproc)) deallocate(yendproc)
   if (allocated(xbegproc)) deallocate(xbegproc)
   if (allocated(xendproc)) deallocate(xendproc)
   if (allocated(p_recep1)) deallocate(p_recep1)
   if (allocated(p_recep2)) deallocate(p_recep2)
   allocate (zbegproc(nproc))
   allocate (zendproc(nproc))
   allocate (ybegproc(nproc))
   allocate (yendproc(nproc))
   allocate (xbegproc(nproc))
   allocate (xendproc(nproc))
   allocate(p_recep1(nproc))
   allocate(p_recep2(nproc))
!
14 format('no cores assigned to compute reciprocal space ',&
      &'contribution')
   if (use_pmecore) then
!$acc update device(use_pmecore)
      if (nrec.eq.0) then
         if (rank.eq.0) write(iout,14)
         call fatal
      end if
      if (nproc-nrec.lt.1) then
         if (rank.eq.0) write(iout,14)
         call fatal
      end if
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
      else
         call MPI_COMM_RANK(comm_rec,rank_bis,ierr)
      end if
!$acc update device(rank_bis)
      deallocate (direct_rank)
   end if
#ifdef _OPENACC
   call copy_data_to_cuda_env(ndir,0)
#endif
!
!     call the dd load balancing routine
!
   call ddpme3d
end
!
!     subroutine allocstep: deallocate arrays and reallocate them with proper size
!     (memory distribution)
!
subroutine allocstep
   use deriv    ,only: mem_alloc_deriv
   use inform   ,only: deb_Path
   use timestat
   use tinMemory,only: mem_get
   implicit none
   real(8) m1,m2,m3
!
   call timer_enter( timer_clear )
   if (deb_Path) call mem_get(m1,m2)

   call mem_alloc_deriv
!
   if (deb_Path) then
      call mem_get(m1,m3)
      if ( m3-m2.ne.0.0 ) then
12       format(" Rank ",I3,"; Forces memory diff",F9.3," Mio")
         print 12, rank, m3-m2
      end if
   end if
!
   call timer_exit( timer_clear,quiet_timers )
end
!
!     subroutine allocstepsrespa: deallocate arrays and reallocate them with proper size
!     (memory distribution)
!
subroutine allocsteprespa(fast)
   use deriv
   use domdec
   use inform ,only: deb_Path
   use potent ,pa=>PotentialAll
   use timestat
   use tinheader
   use tinMemory
   implicit none
   logical fast
   integer i,j
   integer,save:: nb=0,nbr=0
   real(8) m1,m2,m3
!
!     Optimise reallocation of direct force array
!
   if (deb_Path) call mem_get(m1,m2)
   call timer_enter( timer_clear )

   if (fast) then
      call mem_alloc_deriv(cBond)
   else
      call mem_alloc_deriv(cNBond)
   end if

   if (deb_Path) then
      call mem_get(m1,m3)
      if ( m3-m2.ne.0.0 ) then
12       format(" Rank ",I3,"; nbloc ", I10&
            &,"; Forces memory diff",F9.3," Mio")
         print 12, rank,nbloc, m3-m2
      end if
   end if
!
   call timer_exit( timer_clear,quiet_timers )
end
!
!     subroutine ddnumber : get the number of subdivision along each axis for
!     3d spatial decomposition
!
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
   real(r_p) list1(3)
   integer n1,n2,n3,i,res,next
   character*20 keyword
   character*240 record
   character*240 string
10 format('Nx = ',I5,2x,'Ny = ',I5,2x,'Nz = ',I5,2x)
11 format('User defined 3D decompostion ','Nx = ',I5,2x,&
      &'Ny = ',I5,2x,'Nz = ',I5,2x)
12 format('User defined 3D decomposition not compatible with number',&
      &' of cores ','Nx*Ny*Nz = ',I5,2x,'number of procs = ',I5,2x)
14 format('User privileged 1D decompostion ','Nx = ',I5,2x,&
      &'Ny = ',I5,2x,'Nz = ',I5,2x)
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
            if (istep.eq.0.and.verbose.and.dd_verbose) then
               if (ranktot.eq.0) write(iout,11) nxdd,nydd,nzdd
               dd_verbose = .false.
            end if
            return
         else
            if (ranktot.eq.0) then
               write(iout,12) nxdd*nydd*nzdd,nproc
               write(iout,13)
            end if
         end if
      else if (keyword(1:9).eq.'DECOMP1D ' .and. ndir.ne.1) then
         nxdd = 1
         nydd = 1
         nzdd = ndir
         if (ranktot.eq.0.and.istep.eq.0.and.verbose.and.dd_verbose)&
            &write(iout,14)  nxdd,nydd,nzdd
         dd_verbose = .false.
         return
      end if
   end do
!
   ! Force 1Decomp for 2 process
   if (ndir.eq.2) then
      nxdd = 1
      nydd = 1
      nzdd = 2
      if (ranktot.eq.0.and.verbose.and.dd_verbose) then
         write(iout,*) '3D Domain Decomposition'
         write(iout,10) nxdd,nydd,nzdd
         dd_verbose = .false.
      end if
      return
   end if
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
   if (istep.eq.0.and.verbose.and.ndir.gt.1.and.dd_verbose) then
      if (ranktot.eq.0) write(iout,*) '3D Domain Decomposition'
      if (ranktot.eq.0) write(iout,10) nxdd,nydd,nzdd
      dd_verbose = .false.
   end if
   deallocate (d)
end
!
!     subroutine distproc : get the minimum distance between two 3d domains
!
!     nb : to be modified for non cubic unit cells
subroutine distproc(iproc1,iproc2,dist,do3d)
   use domdec
   use tinheader
   implicit none
   integer iproc1,iproc2
   integer i,j
   real(t_p) dist,dist1,dist2,dist3,dist4
   real(t_p) disttemp,disttemp2
   real(t_p) x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
   real(t_p) xtemp(8),ytemp(8),ztemp(8)
   real(t_p) xtempbis(8),ytempbis(8),ztempbis(8)
   real(t_p) xr,yr,zr
   logical do3d
   dist = 10000.0
!
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
!
!         first case, "same" x,y
!
      if ((((x1.le.x3).and.(x2.ge.x3)).or.((x1.ge.x3).and.(x1.le.x4)))&
         &.and. (((y1.le.y3).and.(y2.ge.y3)).or.((y1.ge.y3)&
         &.and.(y1.le.y4))))&
         &then
         dist1 = z1-z4
         call image(0.0_ti_p,0.0_ti_p,dist1)
         dist2 = z3-z2
         call image(0.0_ti_p,0.0_ti_p,dist2)
         dist3 = z1-z3
         call image(0.0_ti_p,0.0_ti_p,dist3)
         dist4 = z2-z4
         call image(0.0_ti_p,0.0_ti_p,dist4)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist3 = abs(dist3)
         dist4 = abs(dist4)
         dist = min(dist1,dist2,dist3,dist4)
!
!         second case, "same" x,z
!
      else if ((((x1.le.x3).and.(x2.ge.x3)).or.((x1.ge.x3)&
         &.and.(x1.le.x4)))&
         &.and. (((z1.le.z3).and.(z2.ge.z3)).or.&
         &((z1.ge.z3).and.(z1.le.z4))))&
         &then
         dist1 = y1-y4
         call image(0.0_ti_p,dist1,0.0_ti_p)
         dist2 = y3-y2
         call image(0.0_ti_p,dist2,0.0_ti_p)
         dist3 = y1-y3
         call image(0.0_ti_p,dist3,0.0_ti_p)
         dist4 = y2-y4
         call image(0.0_ti_p,dist4,0.0_ti_p)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist3 = abs(dist3)
         dist4 = abs(dist4)
         dist = min(dist1,dist2,dist3,dist4)
!
!         third case, "same" y,z
!
      else if ((((y1.le.y3).and.(y2.ge.y3)).or.((y1.ge.y3)&
         &.and.(y1.le.y4)))&
         &.and. (((z1.le.z3).and.(z2.ge.z3)).or.&
         &((z1.ge.z3).and.(z1.le.z4))))&
         &then
         dist1 = x1-x4
         call image(dist1,0.0_ti_p,0.0_ti_p)
         dist2 = x3-x2
         call image(dist2,0.0_ti_p,0.0_ti_p)
         dist3 = x1-x3
         call image(dist3,0.0_ti_p,0.0_ti_p)
         dist4 = x2-x4
         call image(dist4,0.0_ti_p,0.0_ti_p)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist3 = abs(dist3)
         dist4 = abs(dist4)
         dist = min(dist1,dist2,dist3,dist4)
!
!    along one "edge"
!
      else if ((x1.le.x3).and.(x2.ge.x3)) then
         dist = 1000.0_ti_p
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
         dist = 1000.0_ti_p
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
         dist = 1000.0_ti_p
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
         dist = 1000.0_ti_p
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
         dist = 1000.0_ti_p
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
         dist = 1000.0_ti_p
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
!
!       on a "corner"
!
         xtemp(1) = xbegproc(iproc1+1)
         ytemp(1) = ybegproc(iproc1+1)
         ztemp(1) = zbegproc(iproc1+1)
!
         xtemp(2) = xbegproc(iproc1+1)
         ytemp(2) = ybegproc(iproc1+1)
         ztemp(2) = zendproc(iproc1+1)
!
         xtemp(3) = xbegproc(iproc1+1)
         ytemp(3) = yendproc(iproc1+1)
         ztemp(3) = zbegproc(iproc1+1)
!
         xtemp(4) = xbegproc(iproc1+1)
         ytemp(4) = yendproc(iproc1+1)
         ztemp(4) = zendproc(iproc1+1)
!
         xtemp(5) = xendproc(iproc1+1)
         ytemp(5) = ybegproc(iproc1+1)
         ztemp(5) = zbegproc(iproc1+1)
!
         xtemp(6) = xendproc(iproc1+1)
         ytemp(6) = ybegproc(iproc1+1)
         ztemp(6) = zendproc(iproc1+1)
!
         xtemp(7) = xendproc(iproc1+1)
         ytemp(7) = yendproc(iproc1+1)
         ztemp(7) = zbegproc(iproc1+1)
!
         xtemp(8) = xendproc(iproc1+1)
         ytemp(8) = yendproc(iproc1+1)
         ztemp(8) = zendproc(iproc1+1)
!
         xtempbis(1) = xbegproc(iproc2+1)
         ytempbis(1) = ybegproc(iproc2+1)
         ztempbis(1) = zbegproc(iproc2+1)
!
         xtempbis(2) = xbegproc(iproc2+1)
         ytempbis(2) = ybegproc(iproc2+1)
         ztempbis(2) = zendproc(iproc2+1)
!
         xtempbis(3) = xbegproc(iproc2+1)
         ytempbis(3) = yendproc(iproc2+1)
         ztempbis(3) = zbegproc(iproc2+1)
!
         xtempbis(4) = xbegproc(iproc2+1)
         ytempbis(4) = yendproc(iproc2+1)
         ztempbis(4) = zendproc(iproc2+1)
!
         xtempbis(5) = xendproc(iproc2+1)
         ytempbis(5) = ybegproc(iproc2+1)
         ztempbis(5) = zbegproc(iproc2+1)
!
         xtempbis(6) = xendproc(iproc2+1)
         ytempbis(6) = ybegproc(iproc2+1)
         ztempbis(6) = zendproc(iproc2+1)
!
         xtempbis(7) = xendproc(iproc2+1)
         ytempbis(7) = yendproc(iproc2+1)
         ztempbis(7) = zbegproc(iproc2+1)
!
         xtempbis(8) = xendproc(iproc2+1)
         ytempbis(8) = yendproc(iproc2+1)
         ztempbis(8) = zendproc(iproc2+1)
         dist = 1000.0_ti_p
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
      xr = 0.0_ti_p
      yr = 0.0_ti_p
      call image(xr,yr,dist1)
      dist2 = zbegproc(iproc2+1)-zendproc(iproc1+1)
      xr = 0.0_ti_p
      yr = 0.0_ti_p
      call image(xr,yr,dist2)
      dist1 = abs(dist1)
      dist2 = abs(dist2)
      dist = min(dist1,dist2)
   end if
   return
end
!
!     subroutine distprocpart : get the minimum distance between a 3d domain and an atom
!
!     to be modified for non cubic unit cells
subroutine distprocpart(i,iproc,dist,do3d)
   use tinheader
   use atoms
   use cell
   use domdecstuff_inl
   use domdec,only:xbegproc,xendproc,ybegproc,yendproc,&
      &zbegproc,zendproc
   implicit none
   integer iproc
   integer i,j
   real(t_p) dist,dist1,dist2
   real(t_p) disttemp,disttemp2
   real(t_p) x1,x2,x3,y1,y2,y3,z1,z2,z3
   real(t_p) xtemp(8),ytemp(8),ztemp(8)
   real(t_p) xr,yr,zr
!     real(t_p):: zero=0
   logical do3d
!
   x3 = x(i)
   y3 = y(i)
   z3 = z(i)
   call image_inl(x3,y3,z3)
!
   if (do3d) then
      x1 = xbegproc(iproc+1)
      x2 = xendproc(iproc+1)
      y1 = ybegproc(iproc+1)
      y2 = yendproc(iproc+1)
      z1 = zbegproc(iproc+1)
      z2 = zendproc(iproc+1)
!
!       deal with atom exactly on the boundary of the proc's domain
!
!       on the "x" boundary
!
      if (((y1.le.y3).and.(y2.ge.y3)).and.((z1.le.z3).and.(z2.ge.z3))&
         &.and.((x1.eq.x3).or.(x2.eq.x3))) then
         dist = 0.0_ti_p
         return
      end if
!
!       on the "y" boundary
!
      if (((x1.le.x3).and.(x2.ge.x3)).and.((z1.le.z3).and.(z2.ge.z3))&
         &.and.((y1.eq.y3).or.(y2.eq.y3))) then
         dist = 0.0_ti_p
         return
      end if
!
!       on the "z" boundary
!
      if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3))&
         &.and.((z1.eq.z3).or.(z2.eq.z3))) then
         dist = 0.0_ti_p
         return
      end if
!
!         first case, "same" x,y
!
      if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3)))&
         &then
         dist1 = z3-z2
         call image1d_inl(dist1,zcell,zcell2)
         dist2 = z1-z3
         call image1d_inl(dist2,zcell,zcell2)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist = min(dist1,dist2)
!
!         second case, "same" x,z
!
      else if (((x1.le.x3).and.(x2.ge.x3)).and.&
         &((z1.le.z3).and.(z2.ge.z3)))&
         &then
         dist1 = y3-y2
         call image1d_inl(dist1,ycell,ycell2)
         dist2 = y1-y3
         call image1d_inl(dist2,ycell,ycell2)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist = min(dist1,dist2)
!
!         third case, "same" y,z
!
      else if (((y1.le.y3).and.(y2.ge.y3)).and.&
         &((z1.le.z3).and.(z2.ge.z3)))&
         &then
         dist1 = x3-x2
         call image1d_inl(dist1,xcell,xcell2)
         dist2 = x1-x3
         call image1d_inl(dist2,xcell,xcell2)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist = min(dist1,dist2)
!
!     along one "edge"
!
      else if ((x1.le.x3).and.(x2.ge.x3)) then
         dist = 1000.0_ti_p
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
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
      else if ((y1.le.y3).and.(y2.ge.y3)) then
         dist = 1000.0_ti_p
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
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
      else if ((z1.le.z3).and.(z2.ge.z3)) then
         dist = 1000.0_ti_p
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
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
!
      else
!
!       on a "corner"
!
         dist = 1000.0_ti_p
         xtemp(1) = x1
         ytemp(1) = y1
         ztemp(1) = z1
!
         xtemp(2) = x1
         ytemp(2) = y1
         ztemp(2) = z2
!
         xtemp(3) = x1
         ytemp(3) = y2
         ztemp(3) = z1
!
         xtemp(4) = x1
         ytemp(4) = y2
         ztemp(4) = z2
!
         xtemp(5) = x2
         ytemp(5) = y1
         ztemp(5) = z1
!
         xtemp(6) = x2
         ytemp(6) = y1
         ztemp(6) = z2
!
         xtemp(7) = x2
         ytemp(7) = y2
         ztemp(7) = z1
!
         xtemp(8) = x2
         ytemp(8) = y2
         ztemp(8) = z2
!
         do j = 1, 8
            xr = x3 - xtemp(j)
            yr = y3 - ytemp(j)
            zr = z3 - ztemp(j)
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
      end if
   else
      dist1 = zbegproc(iproc+1)-z3
      !xr = 0.0_ti_p
      !yr = 0.0_ti_p
      call image1d_inl(dist1,zcell,zcell2)
      dist2 = z3-zendproc(iproc+1)
      !xr = 0.0_ti_p
      !yr = 0.0_ti_p
      call image1d_inl(dist2,zcell,zcell2)
      dist1 = abs(dist1)
      dist2 = abs(dist2)
      dist = min(dist1,dist2)
   end if
end
!
!     Same a dispprocpart except for the interface
!
subroutine distprocpart1(i,iproc,dist,do3d,x,y,z)
!$acc routine
   use tinheader
   use atoms ,only:n
   use cell
   use domdecstuff_inl
   use domdec,only:xbegproc,xendproc,ybegproc,yendproc,&
      &zbegproc,zendproc
   implicit none
   integer iproc
   integer i,j
   real(t_p) x(n),y(n),z(n)
   real(t_p) dist,dist1,dist2
   real(t_p) disttemp,disttemp2
   real(t_p) x1,x2,x3,y1,y2,y3,z1,z2,z3
   real(t_p) xtemp(8),ytemp(8),ztemp(8)
   real(t_p) xr,yr,zr
!     real(t_p):: zero=0
   logical do3d
!
   x3 = x(i)
   y3 = y(i)
   z3 = z(i)
   call image_inl(x3,y3,z3)
!
   if (do3d) then
      dist = 1000.0_ti_p
      x1 = xbegproc(iproc+1)
      x2 = xendproc(iproc+1)
      y1 = ybegproc(iproc+1)
      y2 = yendproc(iproc+1)
      z1 = zbegproc(iproc+1)
      z2 = zendproc(iproc+1)
!       deal with atom exactly on the boundary of the proc's domain
!
!       on the "x" boundary
!
      if (((y1.le.y3).and.(y2.ge.y3)).and.((z1.le.z3).and.(z2.ge.z3))&
         &.and.((x1.eq.x3).or.(x2.eq.x3))) then
         dist = 0.0_ti_p
         return
      end if
!
!       on the "y" boundary
!
      if (((x1.le.x3).and.(x2.ge.x3)).and.((z1.le.z3).and.(z2.ge.z3))&
         &.and.((y1.eq.y3).or.(y2.eq.y3))) then
         dist = 0.0_ti_p
         return
      end if
!
!       on the "z" boundary
!
      if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3))&
         &.and.((z1.eq.z3).or.(z2.eq.z3))) then
         dist = 0.0_ti_p
         return
      end if
!
!         first case, "same" x,y
!
      if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3)))&
         &then
         dist1 = z3-z2
         call image1d_inl(dist1,zcell,zcell2)
         dist2 = z1-z3
         call image1d_inl(dist2,zcell,zcell2)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist = min(dist1,dist2)
!
!         second case, "same" x,z
!
      else if (((x1.le.x3).and.(x2.ge.x3)).and.&
         &((z1.le.z3).and.(z2.ge.z3)))&
         &then
         dist1 = y3-y2
         call image1d_inl(dist1,ycell,ycell2)
         dist2 = y1-y3
         call image1d_inl(dist2,ycell,ycell2)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist = min(dist1,dist2)
!
!         third case, "same" y,z
!
      else if (((y1.le.y3).and.(y2.ge.y3)).and.&
         &((z1.le.z3).and.(z2.ge.z3)))&
         &then
         dist1 = x3-x2
         call image1d_inl(dist1,xcell,xcell2)
         dist2 = x1-x3
         call image1d_inl(dist2,xcell,xcell2)
         dist1 = abs(dist1)
         dist2 = abs(dist2)
         dist = min(dist1,dist2)
!
!     along one "edge"
!
      else if ((x1.le.x3).and.(x2.ge.x3)) then
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
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
      else if ((y1.le.y3).and.(y2.ge.y3)) then
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
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
      else if ((z1.le.z3).and.(z2.ge.z3)) then
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
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
!
      else
!
!       on a "corner"
!
         xtemp(1) = x1
         ytemp(1) = y1
         ztemp(1) = z1
!
         xtemp(2) = x1
         ytemp(2) = y1
         ztemp(2) = z2
!
         xtemp(3) = x1
         ytemp(3) = y2
         ztemp(3) = z1
!
         xtemp(4) = x1
         ytemp(4) = y2
         ztemp(4) = z2
!
         xtemp(5) = x2
         ytemp(5) = y1
         ztemp(5) = z1
!
         xtemp(6) = x2
         ytemp(6) = y1
         ztemp(6) = z2
!
         xtemp(7) = x2
         ytemp(7) = y2
         ztemp(7) = z1
!
         xtemp(8) = x2
         ytemp(8) = y2
         ztemp(8) = z2
!
         do j = 1, 8
            xr = x3 - xtemp(j)
            yr = y3 - ytemp(j)
            zr = z3 - ztemp(j)
            call image_inl(xr,yr,zr)
            disttemp2 = xr*xr + yr*yr + zr*zr
            disttemp = sqrt(disttemp2)
            if (disttemp.le.dist) dist = disttemp
         end do
      end if
   else
      dist1 = zbegproc(iproc+1)-z3
      !xr = 0.0_ti_p
      !yr = 0.0_ti_p
      call image1d_inl(dist1,zcell,zcell2)
      dist2 = z3-zendproc(iproc+1)
      !xr = 0.0_ti_p
      !yr = 0.0_ti_p
      call image1d_inl(dist2,zcell,zcell2)
      dist1 = abs(dist1)
      dist2 = abs(dist2)
      dist = min(dist1,dist2)
   end if
end

subroutine build_domain_delimiters
   use boxes
   use cell
   use domdec
   implicit none
   integer   i,j,k,iproc
   real(t_p), dimension(nproc):: xbegproctemp,ybegproctemp&
      &,zbegproctemp,xendproctemp,yendproctemp,zendproctemp

   nx_box = xbox/nxdd
   ny_box = ybox/nydd
   nz_box = zbox/nzdd
   na_box = xbox/nxdd
   nb_box = ybox/nydd
   nc_box = zbox/nzdd
!
!     assign processes
!
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

            abegproc(iproc) = xbegproctemp(i)
            aendproc(iproc) = xendproctemp(i)
            bbegproc(iproc) = ybegproctemp(j)
            bendproc(iproc) = yendproctemp(j)
            cbegproc(iproc) = zbegproctemp(k)
            cendproc(iproc) = zendproctemp(k)
         end do
      end do
   end do
!$acc update device(xbegproc,xendproc, &
!$acc ybegproc,yendproc,zbegproc,zendproc) async
end subroutine
!
!     subroutine ddpme: domain decomposition load balancing
!     assign atom sites to MPI processes based on a domain decomposition
!
!
subroutine ddpme3d
   use atoms
   use boxes
   use cell
   use cutoff
   use domdec
   use inform
   use iounit
   use keys
   use mpi
   use mlpot
   use neigh
   use potent
   use mpi
   use tinheader
   use utilcomm ,only: no_commdir
   implicit none
   integer i,j,k,ierr
   integer nprocloc,rankloc,iproc
   real(t_p) xr,yr,zr
   real(t_p) mbuf,vbuf,neigbuf,bigbuf,anibuf
   real(t_p) mshortbuf,vshortbuf,bigshortbuf
   real(t_p) eps1,eps2
   real(t_p), allocatable :: abegproctemp(:),bbegproctemp(:)
   real(t_p), allocatable :: cbegproctemp(:)
   real(t_p), allocatable :: aendproctemp(:),bendproctemp(:)
   real(t_p), allocatable :: cendproctemp(:)
   real(t_p) ar,br,cr
   real(t_p) amin,amax
   real(t_p) bmin,bmax
   real(t_p) cmin,cmax
   real(t_p) nnorm
   real(t_p) a(3),b(3),c(3),nn(3),o(3),pp(3)
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
   nbloc   = 0
   nloc    = 0
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
   ! Get decomposition dimension
   if (nxdd.eq.1.and.nydd.eq.1) then
      Bdecomp1d=.true.
      Bdecomp2d=.false.
      Bdecomp3d=.false.
   else if (nxdd.eq.1) then
      Bdecomp1d=.false.
      Bdecomp2d=.true.
      Bdecomp3d=.false.
   else
      Bdecomp1d=.false.
      Bdecomp2d=.false.
      Bdecomp3d=.true.
   end if
!
   call build_domain_delimiters
!
   eps1   =  5*xcell2*prec_eps
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
      if (abs(ar-xbox2).lt.eps1) ar = ar-0.05*sign(real(nx_box,t_p),ar)
      if (abs(br-ybox2).lt.eps1) br = br-0.05*sign(real(ny_box,t_p),br)
      if (abs(cr-zbox2).lt.eps1) cr = cr-0.05*sign(real(nz_box,t_p),cr)
      n_ar = int((ar+xbox2)/na_box)
      n_br = int((br+ybox2)/nb_box)
      n_cr = int((cr+zbox2)/nc_box)
      iproc = n_ar + n_br*nxdd + n_cr*nydd*nxdd
      repart(i) = iproc
      domlen(repart(i)+1) = domlen(repart(i)+1) + 1
   end do
!$acc update device(repart(:)) async
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

   if (orthogonal.or.Octahedron) then
      scala = 1.0
      scalb = 1.0
      scalc = 1.0
   end if

   if (nproc.gt.1) then
      call MPI_Bcast(scala,1,MPI_RPREC,0,COMM_TINKER,ierr)
      call MPI_Bcast(scalb,1,MPI_RPREC,0,COMM_TINKER,ierr)
      call MPI_Bcast(scalc,1,MPI_RPREC,0,COMM_TINKER,ierr)
   end if

   if (ranktot.eq.0.and.tinkerdebug.gt.0) then
      33 format(' ddpme3d:: Box Form Factor  -x',F7.3,' -y',F7.3,' -z',F7.3)
      print 33, scala, scalb, scalc
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
      mshortbuf = 0.0d0
   end if
!
   vbuf = sqrt(vbuf2)+lbuffer
!
!  take some margin because of torques
!
   mbuf = sqrt(mbuf2)+lbuffer
   vshortbuf = sqrt(vshortbuf2)+lbuffer
   mshortbuf = sqrt(mshortbuf2)+lbuffer
   anibuf    = merge(MLpot_rfield + lbuffer, 0.0_ti_p, use_mlpot)
   neigbuf = lbuffer
!
!     get maximum cutoff value
!
   bigbuf = max(mbuf,vbuf,ddcut,anibuf)
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
   if ( rank.eq.0.and.tinkerdebug.gt.0.and.nproc.gt.1 ) then
46    format(a21,I3,';',$)
47    format(I4,$)
      write(*,46) 'nbigshort_recep    = ',nbigshort_recep
      write(*,47) ( pbigshort_recep(i),i=1,nbigshort_recep)
      write(*,*)
      write(*,46) 'nbig_recep         = ',nbig_recep
      write(*,47) (pbig_recep(i),i=1,nbig_recep)
      write(*,*)
      write(*,46) 'nneig_recep        = ',nneig_recep
      write(*,47) (pneig_recep(i),i=1,nneig_recep)
      write(*,*)
      write(*,46) 'n_recepshort1      = ',n_recepshort1
      write(*,47) (p_recepshort1(i),i=1,n_recepshort1)
      write(*,*)
      write(*,46) 'n_recep1           = ',n_recep1
      write(*,47) (p_recep1(i),i=1,n_recep1)
      write(*,*)
      write(*,46) 'n_recepshort2      = ',n_recepshort2
      write(*,47) (p_recepshort2(i),i=1,n_recepshort2)
      write(*,*)
      write(*,46) 'n_recep2           = ',n_recep2
      write(*,47) (p_recep2(i),i=1,n_recep2)
      write(*,*)
      write(*,46) 'ntorqueshort_recep = ',ntorqueshort_recep
      write(*,47) (ptorqueshort_recep(i),i=1,ntorqueshort_recep)
      write(*,*)
      write(*,46) 'ntorque_recep      = ',ntorque_recep
      write(*,47) (ptorque_recep(i),  i=1,ntorque_recep)
      write(*,*)
   end if
!
180 call orderbuffer(.true.)
    
   deallocate (abegproctemp)
   deallocate (aendproctemp)
   deallocate (bbegproctemp)
   deallocate (bendproctemp)
   deallocate (cbegproctemp)
   deallocate (cendproctemp)
end
!
!     subroutine ddbuild: build lists of neighbors for a given cutoff
!
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
   integer na1,nb1,nc
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

   na1 = max(1,int(3d0*lena*scala/(cut)))
   nb1 = max(1,int(3d0*lenb*scalb/(cut)))
   nc = max(1,int(3d0*lenc*scalc/(cut)))
   !write(*,*) 'rank = ',rank,'na1 = ',na1,'nb1 = ',nb1,'nc = ',nc
   !write(*,*) 'rank = ',rank,'lena = ',lena,amin,amax
   !write(*,*) 'rank = ',rank,'nbig_recep = ',nbig_recep

   lena_cell = lena/na1
   lenb_cell = lenb/nb1
   lenc_cell = lenc/nc


   ncell = na1*nb1*nc

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
         itemp = modulo(i,na1)
         jtemp = modulo(j,nb1)
         ktemp = modulo(k,nc)
         icell = ktemp*nb1*na1+jtemp*na1+itemp+1
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
     do j = 0, nb1-1
       do i = 0, na1-1
     ncelltemp = 0
     lcell_neig = 0
     icell = k*nb1*na1+j*na1+i+1
     do p = -3,3
       do q = -3,3
          do r = -3,3
          if ((p.eq.0).and.(q.eq.0).and.(r.eq.0)) goto 10
   !       if (abs(p)+abs(q)+abs(r).gt.6) goto 10
          tempx = modulo(i+p,na1) 
          tempy = modulo(j+q,nb1) 
          tempz = modulo(k+r,nc) 
          jcell = tempx + tempy*na1 + tempz*nb1*na1 + 1
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
     icell = n_xr + n_yr*na1 + n_zr*nb1*na1+1
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

subroutine AtomDebRepart(ierr)
   use atoms
   use domdec
   use mpi
   use potent
   implicit none
   integer nd,nr,ierr
   integer comm_d,comm_r

   if (use_pmecore) then
      comm_d = comm_dir
      comm_r = comm_rec
   else
      comm_d = COMM_TINKER
      comm_r = COMM_TINKEr
   end if

   call MPI_AllReduce(nloc,nd,1,MPI_INT,MPI_SUM,comm_d,ierr)
   call MPI_AllReduce(nlocrec,nr,1,MPI_INT,MPI_SUM,comm_r,ierr)

11 format("An issue has been detected during reassign process")
12 format(A,' > nloc ',I10,' ntot ',3I10)

   if (nd.ne.n) then
      if (rank.eq.0) print 11
      write(*,12) " direct space",nloc,n,nd; ierr=1;
   end if

   if (rank.eq.0.And.nr.ne.n) write(*,12) " rec    space",nlocrec,n
end subroutine

subroutine AtomDebLocation
   use atomsMirror
   use cell
   use domdec
   use inform
   use mpi
   implicit none
   integer i,j
   real(r_p) boundx,boundy,boundz

   if (deb_Path) print*,'AtomDebLocation'
   boundx = 1.5*xcell2
   boundy = 1.5*ycell2
   boundz = 1.5*zcell2
!$acc wait
!$acc parallel loop present(x,y,z)
   do i = 1, n
      if (abs(x(i)).ge.boundx) print*, "atomi",i,"out of bound",boundx
      if (abs(y(i)).ge.boundy) print*, "atomi",i,"out of bound",boundy
      if (abs(z(i)).ge.boundz) print*, "atomi",i,"out of bound",boundz
   end do
end subroutine
!
!     subroutine midpoint : routine that says whether an interaction between two particules
!     has to be computed within the current domain or not (dd midpoint method)
!
!
subroutine midpoint(xi,yi,zi,xk,yk,zk,docompute)
   use cell
   use domdec
   implicit none
   real(t_p) xi,yi,zi
   real(t_p) xk,yk,zk
   real(t_p) xr,yr,zr
   real(t_p) xrmid,yrmid,zrmid
   logical docompute
!
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
   if (xcell2-abs(xrmid).lt.eps_cell)&
      &xrmid = xrmid-sign(4*eps_cell,xrmid)
   if (ycell2-abs(yrmid).lt.eps_cell)&
      &yrmid = yrmid-sign(4*eps_cell,yrmid)
   if (zcell2-abs(zrmid).lt.eps_cell)&
      &zrmid = zrmid-sign(4*eps_cell,zrmid)

   if   ((zrmid.ge.zbegproc(rank+1)).and.(zrmid.lt.zendproc(rank+1))&
      &.and.(yrmid.ge.ybegproc(rank+1)).and.(yrmid.lt.yendproc(rank+1))&
      &.and.(xrmid.ge.xbegproc(rank+1)).and.(xrmid.lt.xendproc(rank+1)))&
      &then
      docompute = .true.
   end if
   return
end

!
!     subroutine midpointimage : routine that says whether an interaction between two particules
!     has to be computed within the current domain or not (dd midpoint method), also returns
!     minimum image of the distance vector
!
!
subroutine midpointimage(xi,yi,zi,xk,yk,zk,xr,yr,zr,docompute)
   use cell
   use domdec
   implicit none
   real(t_p) xi,yi,zi
   real(t_p) xk,yk,zk
   real(t_p) xr,yr,zr
   real(t_p) xrmid,yrmid,zrmid
   logical docompute
!
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
   if ((xcell2-abs(xrmid)).lt.eps_cell)&
      &xrmid = xrmid-sign(4*eps_cell,xrmid)
   if ((ycell2-abs(yrmid)).lt.eps_cell)&
      &yrmid = yrmid-sign(4*eps_cell,xrmid)
   if ((zcell2-abs(zrmid)).lt.eps_cell)&
      &zrmid = zrmid-sign(4*eps_cell,xrmid)
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
!
subroutine midpointimagetriclinic(xi,yi,zi,xk,yk,zk,xr,yr,zr,docompute)
   use boxes
   use cell
   use domdec
   implicit none
   real(t_p) xi,yi,zi
   real(t_p) xk,yk,zk
   real(t_p) xr,yr,zr
   real(t_p) ai,bi,ci
   real(t_p) ak,bk,ck
   real(t_p) ar,br,cr
   real(t_p) armid,brmid,crmid
   real(t_p) eps1,eps2
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
end
!
!
!     subroutine midpointgroup : routine that says whether an interaction between a
!     group of particules
!     has to be computed within the current domain or not (dd midpoint method, Newton's
!     3rd law)
!
!
subroutine midpointgroup(pos,nb,docompute)
   use domdec
   implicit none
   integer nb,i
   real(t_p) pos(3,nb),posdir(3,nb)
   real(t_p) xrmid,yrmid,zrmid
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
   !TODO 1.2 I think there's a check missing
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
!
!     subroutine ddpme3dnpt : rescale geomtry of the domains and recompute the related quantities
!     for communications
!
subroutine ddpme3dnpt(scaleiso,istep)
   real(r_p) scaleiso
   integer istep
   real(r_p) scale(3)

   scale(:)=scaleiso
   call ddpme3dnptaniso(scale,istep)

end subroutine ddpme3dnpt

subroutine ddpme3dnptaniso(scale,istep)
   use boxes
   use cutoff
   use domdec
   use neigh
   use potent
   use mpi
   use tinheader
   implicit none
   integer nprocloc,rankloc,commloc,iproc
   integer i
   integer istep,modnl
   real(t_p) mbuf,vbuf,neigbuf,bigbuf
   real(t_p) mshortbuf,vshortbuf,bigshortbuf
   real(r_p) scale(3)
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
   call build_domain_delimiters
!
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

   if (nproc.eq.1) return

   if (use_pmecore) then
      if (ndir.eq.1.and.nrec.eq.1) return
      if (ndir.eq.1.and.nrec.gt.1) goto 80
      if (rank.gt.ndir-1) goto 80
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
      mbuf = 0.0_ti_p
      mshortbuf = 0.0_ti_p
   end if
!
   vbuf = sqrt(vbuf2)+lbuffer
!
!  take some margin because of torques
!
   mbuf = sqrt(mbuf2)+lbuffer
   vshortbuf = sqrt(vshortbuf2)+lbuffer
   mshortbuf = sqrt(mshortbuf2)+lbuffer
   neigbuf = lbuffer
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
!
80 call orderbuffer_gpu(.false.)
   call orderbufferrec_gpu

end
