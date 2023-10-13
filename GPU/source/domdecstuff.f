c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_macro.h"
      module domdecstuff_inl
        contains
#include "image.f.inc"
      end module

c     ##################################################################
c     ##                                                              ##
c     ##  subroutine drivermpi  --  driver for MPI related quantities ##
c     ##                            (3d spatial decomposition)        ##
c     ##                                                              ##
c     ##################################################################
c
c
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
c
      ndir = nproc - nrec
c
c     MPI : get the atoms repartition over the processes
c
c     allocate global arrays
c
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
c
 14   format('no cores assigned to compute reciprocal space ',
     &       'contribution')
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
     $        direct_group,ierr)
         call MPI_Group_excl(total_group,ndir,direct_rank,
     $        rec_group,ierr)
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
c
c     call the dd load balancing routine
c
      call ddpme3d
      end
c
c     subroutine allocstep: deallocate arrays and reallocate them with proper size
c     (memory distribution)
c
      subroutine allocstep
      use deriv    ,only: mem_alloc_deriv
      use inform   ,only: deb_Path
      use timestat
      use tinMemory,only: mem_get
      implicit none
      real(8) m1,m2,m3
c
      call timer_enter( timer_clear )
      if (deb_Path) call mem_get(m1,m2)

      call mem_alloc_deriv
c
      if (deb_Path) then
         call mem_get(m1,m3)
         if ( m3-m2.ne.0.0 ) then
 12   format(" Rank ",I3,"; Forces memory diff",F9.3," Mio")
            print 12, rank, m3-m2
         end if
      end if
c
      call timer_exit( timer_clear,quiet_timers )
      end
c
c     subroutine allocstepsrespa: deallocate arrays and reallocate them with proper size
c     (memory distribution)
c
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
c
c     Optimise reallocation of direct force array
c
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
 12   format(" Rank ",I3,"; nbloc ", I10
     &      ,"; Forces memory diff",F9.3," Mio")
            print 12, rank,nbloc, m3-m2
         end if
      end if
c
      call timer_exit( timer_clear,quiet_timers )
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
      real(r_p) list1(3)
      integer n1,n2,n3,i,res,next
      character*20 keyword
      character*240 record
      character*240 string
 10   format('Nx = ',I5,2x,'Ny = ',I5,2x,'Nz = ',I5,2x)
 11   format('User defined 3D decompostion ','Nx = ',I5,2x,
     $    'Ny = ',I5,2x,'Nz = ',I5,2x)
 12   format('User defined 3D decomposition not compatible with number',
     $    ' of cores ','Nx*Ny*Nz = ',I5,2x,'number of procs = ',I5,2x)
 14   format('User privileged 1D decompostion ','Nx = ',I5,2x,
     $    'Ny = ',I5,2x,'Nz = ',I5,2x)
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
           if (ranktot.eq.0.and.istep.eq.0.and.verbose.and.dd_verbose)
     &        write(iout,14)  nxdd,nydd,nzdd
              dd_verbose = .false.
           return
         end if
      end do
c
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
      if (istep.eq.0.and.verbose.and.ndir.gt.1.and.dd_verbose) then
        if (ranktot.eq.0) write(iout,*) '3D Domain Decomposition'
        if (ranktot.eq.0) write(iout,10) nxdd,nydd,nzdd
        dd_verbose = .false.
      end if
      deallocate (d)
      end
c
c     subroutine distproc : get the minimum distance between two 3d domains
c
c     nb : to be modified for non cubic unit cells
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
c
c         second case, "same" x,z
c
        else if ((((x1.le.x3).and.(x2.ge.x3)).or.((x1.ge.x3)
     $   .and.(x1.le.x4)))
     $   .and. (((z1.le.z3).and.(z2.ge.z3)).or.
     $   ((z1.ge.z3).and.(z1.le.z4))))
     $     then
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
c
c         third case, "same" y,z
c
        else if ((((y1.le.y3).and.(y2.ge.y3)).or.((y1.ge.y3)
     $   .and.(y1.le.y4)))
     $   .and. (((z1.le.z3).and.(z2.ge.z3)).or.
     $   ((z1.ge.z3).and.(z1.le.z4))))
     $     then
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
c
c    along one "edge"
c
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
c
c     subroutine distprocpart : get the minimum distance between a 3d domain and an atom
c
c     to be modified for non cubic unit cells
      subroutine distprocpart(i,iproc,dist,do3d)
      use tinheader
      use atoms
      use cell
      use domdecstuff_inl
      use domdec,only:xbegproc,xendproc,ybegproc,yendproc,
     &                zbegproc,zendproc
      implicit none
      integer iproc
      integer i,j
      real(t_p) dist,dist1,dist2
      real(t_p) disttemp,disttemp2
      real(t_p) x1,x2,x3,y1,y2,y3,z1,z2,z3
      real(t_p) xtemp(8),ytemp(8),ztemp(8)
      real(t_p) xr,yr,zr
c     real(t_p):: zero=0
      logical do3d
c
      x3 = x(i)
      y3 = y(i)
      z3 = z(i)
      call image_inl(x3,y3,z3)
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
     $     .and.((x1.eq.x3).or.(x2.eq.x3))) then
           dist = 0.0_ti_p
           return
        end if
c
c       on the "y" boundary
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((z1.le.z3).and.(z2.ge.z3))
     $     .and.((y1.eq.y3).or.(y2.eq.y3))) then
           dist = 0.0_ti_p
           return
        end if
c
c       on the "z" boundary
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3))
     $     .and.((z1.eq.z3).or.(z2.eq.z3))) then
           dist = 0.0_ti_p
           return
        end if
c
c         first case, "same" x,y
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3)))
     $     then
          dist1 = z3-z2
          call image1d_inl(dist1,zcell,zcell2)
          dist2 = z1-z3
          call image1d_inl(dist2,zcell,zcell2)
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
          call image1d_inl(dist1,ycell,ycell2)
          dist2 = y1-y3
          call image1d_inl(dist2,ycell,ycell2)
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
          call image1d_inl(dist1,xcell,xcell2)
          dist2 = x1-x3
          call image1d_inl(dist2,xcell,xcell2)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist = min(dist1,dist2)
c
c     along one "edge"
c
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
c
        else
c
c       on a "corner"
c
          dist = 1000.0_ti_p
          xtemp(1) = x1
          ytemp(1) = y1
          ztemp(1) = z1
c
          xtemp(2) = x1
          ytemp(2) = y1
          ztemp(2) = z2
c
          xtemp(3) = x1
          ytemp(3) = y2
          ztemp(3) = z1
c
          xtemp(4) = x1
          ytemp(4) = y2
          ztemp(4) = z2
c
          xtemp(5) = x2
          ytemp(5) = y1
          ztemp(5) = z1
c
          xtemp(6) = x2
          ytemp(6) = y1
          ztemp(6) = z2
c
          xtemp(7) = x2
          ytemp(7) = y2
          ztemp(7) = z1
c
          xtemp(8) = x2
          ytemp(8) = y2
          ztemp(8) = z2
c
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
c
c     Same a dispprocpart except for the interface
c
      subroutine distprocpart1(i,iproc,dist,do3d,x,y,z)
!$acc routine
      use tinheader
      use atoms ,only:n
      use cell
      use domdecstuff_inl
      use domdec,only:xbegproc,xendproc,ybegproc,yendproc,
     &                zbegproc,zendproc
      implicit none
      integer iproc
      integer i,j
      real(t_p) x(n),y(n),z(n)
      real(t_p) dist,dist1,dist2
      real(t_p) disttemp,disttemp2
      real(t_p) x1,x2,x3,y1,y2,y3,z1,z2,z3
      real(t_p) xtemp(8),ytemp(8),ztemp(8)
      real(t_p) xr,yr,zr
c     real(t_p):: zero=0
      logical do3d
c
      x3 = x(i)
      y3 = y(i)
      z3 = z(i)
      call image_inl(x3,y3,z3)
c
      if (do3d) then
        dist = 1000.0_ti_p
        x1 = xbegproc(iproc+1)
        x2 = xendproc(iproc+1)
        y1 = ybegproc(iproc+1)
        y2 = yendproc(iproc+1)
        z1 = zbegproc(iproc+1)
        z2 = zendproc(iproc+1)
c       deal with atom exactly on the boundary of the proc's domain
c
c       on the "x" boundary
c
        if (((y1.le.y3).and.(y2.ge.y3)).and.((z1.le.z3).and.(z2.ge.z3))
     $   .and.((x1.eq.x3).or.(x2.eq.x3))) then
          dist = 0.0_ti_p
          return
        end if
c
c       on the "y" boundary
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((z1.le.z3).and.(z2.ge.z3))
     $   .and.((y1.eq.y3).or.(y2.eq.y3))) then
          dist = 0.0_ti_p
          return
        end if
c
c       on the "z" boundary
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3))
     $   .and.((z1.eq.z3).or.(z2.eq.z3))) then
          dist = 0.0_ti_p
          return
        end if
c
c         first case, "same" x,y
c
        if (((x1.le.x3).and.(x2.ge.x3)).and.((y1.le.y3).and.(y2.ge.y3)))
     $     then
          dist1 = z3-z2
          call image1d_inl(dist1,zcell,zcell2)
          dist2 = z1-z3
          call image1d_inl(dist2,zcell,zcell2)
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
          call image1d_inl(dist1,ycell,ycell2)
          dist2 = y1-y3
          call image1d_inl(dist2,ycell,ycell2)
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
          call image1d_inl(dist1,xcell,xcell2)
          dist2 = x1-x3
          call image1d_inl(dist2,xcell,xcell2)
          dist1 = abs(dist1)
          dist2 = abs(dist2)
          dist = min(dist1,dist2)
c
c     along one "edge"
c
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
c
        else
c
c       on a "corner"
c
          xtemp(1) = x1
          ytemp(1) = y1
          ztemp(1) = z1
c
          xtemp(2) = x1
          ytemp(2) = y1
          ztemp(2) = z2
c
          xtemp(3) = x1
          ytemp(3) = y2
          ztemp(3) = z1
c
          xtemp(4) = x1
          ytemp(4) = y2
          ztemp(4) = z2
c
          xtemp(5) = x2
          ytemp(5) = y1
          ztemp(5) = z1
c
          xtemp(6) = x2
          ytemp(6) = y1
          ztemp(6) = z2
c
          xtemp(7) = x2
          ytemp(7) = y2
          ztemp(7) = z1
c
          xtemp(8) = x2
          ytemp(8) = y2
          ztemp(8) = z2
c
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
      real(t_p), dimension(nproc):: xbegproctemp,ybegproctemp
     &         ,zbegproctemp,xendproctemp,yendproctemp,zendproctemp

      nx_box = xbox/nxdd
      ny_box = ybox/nydd
      nz_box = zbox/nzdd
c
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

      !xendproctemp(nxdd) = xbox2
      !yendproctemp(nydd) = ybox2
      !yendproctemp(nzdd) = zbox2
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
!$acc update device(xbegproc,xendproc,
!$acc& ybegproc,yendproc,zbegproc,zendproc) async
      end subroutine
c
c     subroutine ddpme: domain decomposition load balancing
c     assign atom sites to MPI processes based on a domain decomposition
c
c
      subroutine ddpme3d
      use atoms
      use ani
      use boxes
      use cell
      use cutoff
      use domdec
      use iounit
      use sizes
      use neigh
      use potent
      use mpi
      use tinheader
      use utilcomm ,only: no_commdir
      implicit none
      integer i,j,k
      integer nprocloc,rankloc,iproc,iproc1
      real(t_p) xr,yr,zr
      real(t_p) dist,distRatio
      real(t_p) mbuf,vbuf,torquebuf,neigbuf,bigbuf
      real(t_p) mshortbuf,vshortbuf,torqueshortbuf,bigshortbuf
      real(t_p) anibuf,distRatio1
      real(t_p) eps1
      integer p,q,r,numneig
      integer temp_x,temp_y,temp_z,tempproc
      integer, allocatable :: neigproc(:,:),numneigproc(:),filledproc(:)
c
      integer, allocatable :: nrecep1(:),precep1(:,:)
      integer, allocatable :: nrecep2(:),precep2(:,:)
      integer, allocatable :: ntorquerecep(:),ptorquerecep(:,:)
      integer, allocatable :: nbigrecep(:),pbigrecep(:,:)
c
      integer, allocatable :: nrecepshort1(:),precepshort1(:,:)
      integer, allocatable :: nrecepshort2(:),precepshort2(:,:)
      integer, allocatable :: ntorqueshortrecep(:),
     $                        ptorqueshortrecep(:,:)
      integer, allocatable :: nbigshortrecep(:),pbigshortrecep(:,:)
      integer, allocatable :: nneigrecep(:),pneigrecep(:,:)
c
 1000 format(' Warning, less than 10 atoms on process number',I6,1x,
     $   ' number of cores may be too high compared to the number of ',
     $   'atoms')
c
      nbloc = 0
      nloc  = 0
      !nlocrec = 0
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
      allocate (ntorquerecep(nproc))
      allocate (ptorquerecep(nproc,nproc))
      allocate (nbigrecep(nproc))
      allocate (pbigrecep(nproc,nproc))
c
c     also deal with short range non bonded forces (respa-1 like integrators)
c
      allocate (nrecepshort1(nproc))
      allocate (precepshort1(nproc,nproc))
      allocate (nrecepshort2(nproc))
      allocate (precepshort2(nproc,nproc))
      allocate (ntorqueshortrecep(nproc))
      allocate (ptorqueshortrecep(nproc,nproc))
      allocate (nbigshortrecep(nproc))
      allocate (pbigshortrecep(nproc,nproc))
c
      allocate (nneigrecep(nproc))
      allocate (pneigrecep(nproc,nproc))
      nrecep1      = 0
      nrecep2      = 0
      ntorquerecep = 0
      nbigrecep    = 0
      nneigrecep   = 0
      precep1      = 0
      precep2      = 0
      ptorquerecep = 0
      pbigrecep    = 0
c
      nrecepshort1      = 0
      nrecepshort2      = 0
      ntorqueshortrecep = 0
      nbigshortrecep    = 0
      precepshort1      = 0
      precepshort2      = 0
      ptorqueshortrecep = 0
      pbigshortrecep    = 0
c
      pneigrecep = 0
c
      allocate (neigproc(26,nproc))
      allocate (numneigproc(nproc))
      allocate (filledproc(nproc))
      neigproc     = 0
      numneigproc  = 0
      filledproc   = 0
      xbegproc     = 0_ti_p
      ybegproc     = 0_ti_p
      zbegproc     = 0_ti_p
      xendproc     = 0_ti_p
      yendproc     = 0_ti_p
      zendproc     = 0_ti_p
      repart       = 0
      domlen       = 0
      glob         = 0
      loc          = 0
!$acc update device(glob,loc) async
c
c     get the number of subdivision along each axis for dd
c
      call ddnumber(nprocloc,0)

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
c
      call build_domain_delimiters
      eps1   =  5*xcell2*prec_eps
c
c     count number of particules per domain
c
      domlen = 0
      do i = 1, n
        xr = x(i)
        yr = y(i)
        zr = z(i)
        call image(xr,yr,zr)
        ! Avoid box limits
        if ((xcell2-abs(xr)).lt.eps1) xr = xr-0.05*sign(nx_box,xr)
        if ((ycell2-abs(yr)).lt.eps1) yr = yr-0.05*sign(ny_box,yr)
        if ((zcell2-abs(zr)).lt.eps1) zr = zr-0.05*sign(nz_box,zr)
        do iproc = 0, nprocloc-1
          if ((zr.ge.zbegproc(iproc+1)).and.(zr.lt.zendproc(iproc+1))
     &   .and.(yr.ge.ybegproc(iproc+1)).and.(yr.lt.yendproc(iproc+1))
     &   .and.(xr.ge.xbegproc(iproc+1)).and.(xr.lt.xendproc(iproc+1)))
     &    then
             repart(i) = iproc
             domlen(repart(i)+1) = domlen(repart(i)+1) + 1
          end if
        end do
      end do
!$acc update device(repart(:)) async
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

      if (rank.eq.0.and.nproc.gt.1.and.tinkerdebug.gt.0) then
         print*,'---  ddpme3d'
 14      format(A9)
 15      format(F9.3)
 16      format(I9)
         write(*,14,advance='no') 'xbegproc '
         do i = 1,nproc
            write(*,15,advance='no') xbegproc(i)
         end do
         write(*,*)
         write(*,14,advance='no') 'xendproc '
         do i = 1,nproc
            write(*,15,advance='no') xendproc(i)
         end do
         write(*,*)
         write(*,14,advance='no') 'ybegproc '
         do i = 1,nproc
            write(*,15,advance='no') ybegproc(i)
         end do
         write(*,*)
         write(*,14,advance='no') 'yendproc '
         do i = 1,nproc
            write(*,15,advance='no') yendproc(i)
         end do
         write(*,*)
         write(*,14,advance='no') 'zbegproc '
         do i = 1,nproc
            write(*,15,advance='no') zbegproc(i)
         end do
         write(*,*)
         write(*,14,advance='no') 'zendproc '
         do i = 1,nproc
            write(*,15,advance='no') zendproc(i)
         end do
         write(*,*)
         write(*,14,advance='no') 'domlen '
         do i = 1,nproc
            write(*,16,advance='no') domlen(i)
         end do
         write(*,*)
      end if

cc
cc     warning bug in iterative load balancing
cc
cc
cc     change size with neighbors
cc
c      exchgsizex = nx_box/20
c      nitmax = 0!3
c      ndiff  = n/(20*nprocloc)
cc
c      do it = 1, nitmax
c        do iproc = 1, nprocloc
c          x2 = xendproc(iproc)
c          y1 = ybegproc(iproc)
c          z1 = zbegproc(iproc)
cc
cc    check neighbors
cc
c          do i = 1, numneigproc(iproc)
c            iproc1 = neigproc(i,iproc)
c            x3 = xbegproc(iproc1)
c            y3 = ybegproc(iproc1)
c            z3 = zbegproc(iproc1)
cc
c            if ((x2.eq.x3).and.(y1.eq.y3).and.(z1.eq.z3)) then
c              if ((domlen(iproc1).ge.domlen(iproc)+ndiff)) then
c                xendproc(iproc)   = xendproc(iproc) +
c     $            exchgsizex
c                xbegproc(iproc1)   = xbegproc(iproc1) + exchgsizex
c              end if
c            end if
c          end do
c        end do
cc
cc     count number of particules per domain
cc
c        domlen = 0
c        do i = 1, n
c          xr = x(i)
c          yr = y(i)
c          zr = z(i)
c          call image(xr,yr,zr)
c          if ((xcell2-abs(xr)).lt.eps1) xr = xr-sign(eps2,xr)
c          if ((ycell2-abs(yr)).lt.eps1) yr = yr-sign(eps2,yr)
c          if ((zcell2-abs(zr)).lt.eps1) zr = zr-sign(eps2,zr)
c          do iproc = 0, nprocloc-1
c            if ((zr.ge.zbegproc(iproc+1)).and.
c     $        (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
c     $        .and.(yr.lt.yendproc(iproc+1))
c     $        .and.(xr.ge.xbegproc(iproc+1))
c     $        .and.(xr.lt.xendproc(iproc+1))) then
c               repart(i) = iproc
c              domlen(repart(i)+1) = domlen(repart(i)+1) + 1
c            end if
c          end do
c        end do
c      end do
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
      p_send1      = 0
      p_recep1     = 0
      p_send2      = 0
      p_recep2     = 0
      ptorque_send = 0
      ptorque_recep = 0
      pbig_send     = 0
      pbig_recep    = 0
      n_recep1      = 0
      n_recep2      = 0
      ntorque_recep = 0
      nneig_recep   = 0
      nbig_recep    = 0
c
      p_sendshort1  = 0
      p_recepshort1 = 0
      p_sendshort2  = 0
      p_recepshort2 = 0
      ptorqueshort_send  = 0
      ptorqueshort_recep = 0
      pbigshort_send  = 0
      pbigshort_recep = 0
      n_recepshort1   = 0
      n_recepshort2   = 0
      ntorqueshort_recep = 0
      nbigshort_recep    = 0
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) then
        goto 80
      end if
c
c     choose cutoff depending on electrostatic interaction
c
      if (use_mpole) then
        mbuf = sqrt(mbuf2)
        mshortbuf = sqrt(mshortbuf2)
      else if (use_charge) then
        mbuf = sqrt(cbuf2)
        mshortbuf = sqrt(cshortbuf2)
      else
        mbuf = 0.0_ti_p
        mshortbuf = 0.0_ti_p
      end if
c
      vbuf = sqrt(vbuf2)+2.0_ti_p
      vshortbuf = sqrt(vshortbuf2)+2.0_ti_p
      torquebuf = mbuf + lbuffer
      torqueshortbuf = mshortbuf + lbuffer
      torquebuf2 = torquebuf*torquebuf 
      torqueshortbuf2 = torqueshortbuf*torqueshortbuf 
      anibuf  = 0
      if (use_mlpot) then
         anibuf  = MLpotcut + lbuffer
      end if
      neigbuf = anibuf
c
c     get maximum cutoff value
c
      bigbuf = max(torquebuf,vbuf,ddcut,anibuf)
      bigshortbuf = max(torqueshortbuf,vshortbuf,ddcut)

      ! Set the distannce ratio ragarding direct space neighbor process
      if (no_commdir) then
         distRatio  = 1
         distRatio1 = 1
      else
         distRatio  = 0.5
         distRatio1 = 0.5
      end if
      if (use_mlpot) then
         distRatio1 = 1
      end if
      distRatio1 = max(distRatio1,distRatio)
c
      do iproc = 0, nprocloc-1
        do iproc1 = 0, nprocloc-1
          if (iproc.eq.iproc1) goto 70
          call distproc(iproc1,iproc,dist,.true.)
          if (dist.le.(mbuf*distRatio)) then
            nrecep1(iproc1+1) = nrecep1(iproc1+1)+1
            precep1(nrecep1(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(vbuf/2)) then
            nrecep2(iproc1+1) = nrecep2(iproc1+1)+1
            precep2(nrecep2(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(torquebuf*distRatio)) then
            ntorquerecep(iproc1+1) = ntorquerecep(iproc1+1)+1
            ptorquerecep(ntorquerecep(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(bigbuf*distRatio1)) then
            nbigrecep(iproc1+1) = nbigrecep(iproc1+1)+1
            pbigrecep(nbigrecep(iproc1+1),iproc1+1) = iproc
          end if
c
c    also deal with short range non bonded interactions (respa-1 like integrators)
c
          if (dist.le.(mshortbuf*distRatio)) then
            nrecepshort1(iproc1+1) = nrecepshort1(iproc1+1)+1
            precepshort1(nrecepshort1(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(vshortbuf/2)) then
            nrecepshort2(iproc1+1) = nrecepshort2(iproc1+1)+1
            precepshort2(nrecepshort2(iproc1+1),iproc1+1) = iproc
          end if
          if (dist.le.(torquebuf*distRatio)) then
            ntorqueshortrecep(iproc1+1) = ntorqueshortrecep(iproc1+1)+1
            ptorqueshortrecep(ntorqueshortrecep(iproc1+1),iproc1+1) = 
     $         iproc
          end if
          if (dist.le.(bigshortbuf*distRatio)) then
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
      ntorque_recep = ntorquerecep(rankloc+1)
      ptorque_recep = ptorquerecep(:,rankloc+1)
      nbig_recep = nbigrecep(rankloc+1)
      pbig_recep = pbigrecep(:,rankloc+1)
c
      n_recepshort1 = nrecepshort1(rankloc+1)
      p_recepshort1 = precepshort1(:,rankloc+1)
      n_recepshort2 = nrecepshort2(rankloc+1)
      p_recepshort2 = precepshort2(:,rankloc+1)
      ntorqueshort_recep = ntorqueshortrecep(rankloc+1)
      ptorqueshort_recep = ptorqueshortrecep(:,rankloc+1)
      nbigshort_recep = nbigshortrecep(rankloc+1)
      pbigshort_recep = pbigshortrecep(:,rankloc+1)
c
      nneig_recep = nneigrecep(rankloc+1)
      pneig_recep = pneigrecep(:,rankloc+1)
!$acc update device(p_recep1,p_recep2)
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
      ntorque_send = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, ntorquerecep(iproc+1)
            if (ptorquerecep(i,iproc+1).eq.rankloc) then
              ntorque_send = ntorque_send + 1
              ptorque_send(ntorque_send) = iproc
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
      ntorqueshort_send = 0
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          do i = 1, ntorqueshortrecep(iproc+1)
            if (ptorqueshortrecep(i,iproc+1).eq.rankloc) then
              ntorqueshort_send = ntorqueshort_send + 1
              ptorqueshort_send(ntorqueshort_send) = iproc
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
      if ( rank.eq.0.and.tinkerdebug.gt.0.and.nproc.gt.1 ) then
 46       format(a21,I3,';',$)
 47       format(I3,$)
          write(*,46) 'nbigshort_recep    = ',nbigshort_recep
          write(*,47) ( pbigshort_recep(i),i=1,nbigshort_recep)
          write(*,*)
          write(*,46) 'nbig_recep         = ',nbig_recep
          write(*,47) (pbig_recep(i),i=1,nbig_recep)
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
c
 80   call orderbuffer(.true.)
      deallocate (filledproc)
      deallocate (numneigproc)
      deallocate (neigproc)
      deallocate (nrecep1)
      deallocate (precep1)
      deallocate (nrecep2)
      deallocate (precep2)
      deallocate (nbigrecep)
      deallocate (pbigrecep)
      deallocate (ntorquerecep)
      deallocate (ptorquerecep)
      deallocate (nrecepshort1)
      deallocate (precepshort1)
      deallocate (nrecepshort2)
      deallocate (precepshort2)
      deallocate (nbigshortrecep)
      deallocate (pbigshortrecep)
      deallocate (ntorqueshortrecep)
      deallocate (ptorqueshortrecep)
      deallocate (nneigrecep)
      deallocate (pneigrecep)
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

 11   format("An issue has been detected during reassign process")
 12   format(A,' > nloc ',I10,' ntot ',3I10)

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
c
c     subroutine midpoint : routine that says whether an interaction between two particules
c     has to be computed within the current domain or not (dd midpoint method)
c
c
      subroutine midpoint(xi,yi,zi,xk,yk,zk,docompute)
      use cell
      use domdec
      implicit none
      real(t_p) xi,yi,zi
      real(t_p) xk,yk,zk
      real(t_p) xr,yr,zr
      real(t_p) xrmid,yrmid,zrmid
      logical docompute
c
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
      if (xcell2-abs(xrmid).lt.eps_cell)
     &   xrmid = xrmid-sign(4*eps_cell,xrmid)
      if (ycell2-abs(yrmid).lt.eps_cell)
     &   yrmid = yrmid-sign(4*eps_cell,yrmid)
      if (zcell2-abs(zrmid).lt.eps_cell)
     &   zrmid = zrmid-sign(4*eps_cell,zrmid)

      if   ((zrmid.ge.zbegproc(rank+1)).and.(zrmid.lt.zendproc(rank+1))
     $ .and.(yrmid.ge.ybegproc(rank+1)).and.(yrmid.lt.yendproc(rank+1))
     $ .and.(xrmid.ge.xbegproc(rank+1)).and.(xrmid.lt.xendproc(rank+1)))
     $   then
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
      real(t_p) xi,yi,zi
      real(t_p) xk,yk,zk
      real(t_p) xr,yr,zr
      real(t_p) xrmid,yrmid,zrmid
      logical docompute
c
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
      if ((xcell2-abs(xrmid)).lt.eps_cell)
     &   xrmid = xrmid-sign(4*eps_cell,xrmid)
      if ((ycell2-abs(yrmid)).lt.eps_cell)
     &   yrmid = yrmid-sign(4*eps_cell,xrmid)
      if ((zcell2-abs(zrmid)).lt.eps_cell)
     &   zrmid = zrmid-sign(4*eps_cell,xrmid)
      if ((zrmid.ge.zbegproc(rank+1)).and.
     &    (zrmid.lt.zendproc(rank+1)).and.(yrmid.ge.ybegproc(rank+1))
     &  .and.(yrmid.lt.yendproc(rank+1))
     &  .and.(xrmid.ge.xbegproc(rank+1))
     &  .and.(xrmid.lt.xendproc(rank+1))) then
        docompute = .true.
      end if
      return
      end
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
      real(t_p) pos(3,nb),posdir(3,nb)
      real(t_p) xrmid,yrmid,zrmid
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
      !TODO 1.2 I think there's a check missing
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
      real(r_p) scaleiso
      integer istep
      real(r_p) scale(3)

      scale(:)=scaleiso
      call ddpme3dnptaniso(scale,istep)

      end subroutine ddpme3dnpt

      subroutine ddpme3dnptaniso(scale,istep)
      use cutoff
      use domdec
      use neigh
      use potent
      use mpi
      use tinheader
      implicit none
      integer nprocloc,rankloc,commloc,iproc
      integer i,status(MPI_STATUS_SIZE),tag,ierr
      integer istep,modnl
      real(t_p) mbuf,vbuf,torquebuf,neigbuf,bigbuf
      real(t_p) mshortbuf,vshortbuf,torqueshortbuf,bigshortbuf
      real(t_p) dist
      real(r_p) scale(3)
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
      call build_domain_delimiters

      if (modnl.ne.0) return
c
      allocate (reqsend(nprocloc))
      allocate (reqrec(nprocloc))
      allocate (bufcount(nprocloc))
      allocate (buffer(nprocloc,nprocloc))
c
c     get the processes to receive data from : for commstep and for commtorque
c
      p_send1    = 0
      p_recep1   = 0
      p_send2    = 0
      p_recep2   = 0
      ptorque_send  = 0
      ptorque_recep = 0
      pbig_send  = 0
      pbig_recep = 0
      n_recep1   = 0
      n_recep2   = 0
      ntorque_recep = 0
      nneig_recep   = 0
      nbig_recep    = 0
c
      p_sendshort1   = 0
      p_recepshort1  = 0
      p_sendshort2   = 0
      p_recepshort2  = 0
      ptorqueshort_send  = 0
      ptorqueshort_recep = 0
      pbigshort_send  = 0
      pbigshort_recep = 0
      n_recepshort1   = 0
      n_recepshort2   = 0
      ntorqueshort_recep = 0
      nbigshort_recep = 0

      if (nproc.eq.1) return

      if (use_pmecore) then
        if (ndir.eq.1.and.nrec.eq.1) return
        if (ndir.eq.1.and.nrec.gt.1) goto 80
        if (rank.gt.ndir-1) goto 80
      end if
c
c     choose cutoff depending on electrostatic interaction
c
      if (use_mpole) then
        mbuf = sqrt(mbuf2)
        mshortbuf = sqrt(mshortbuf2)
      else if (use_charge) then
        mbuf = sqrt(cbuf2)
        mshortbuf = sqrt(cshortbuf2)
      else
        mbuf = 0.0_ti_p
        mshortbuf = 0.0_ti_p
      end if
c
      vbuf = sqrt(vbuf2)+2.0_ti_p
      vshortbuf = sqrt(vshortbuf2)+2.0_ti_p
      torquebuf = mbuf + lbuffer
      torqueshortbuf = mshortbuf + lbuffer
      torquebuf2 = torquebuf*torquebuf
      torqueshortbuf2 = torqueshortbuf*torqueshortbuf
      neigbuf = lbuffer
c
c     get maximum cutoff value
c
      bigbuf = max(torquebuf,vbuf,ddcut)
      bigshortbuf = max(torqueshortbuf,vshortbuf,ddcut)
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
          if (dist.le.(torquebuf/2)) then
            ntorque_recep = ntorque_recep+1
            ptorque_recep(ntorque_recep) = iproc
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
          if (dist.le.(torquebuf/2)) then
            ntorqueshort_recep = ntorqueshort_recep+1
            ptorqueshort_recep(ntorqueshort_recep) = iproc
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
!$acc update device(p_recep1,p_recep2) async
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
      ntorque_send = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(ntorque_recep,1,MPI_INT,iproc,tag,
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
          call MPI_ISEND(ptorque_recep,ntorque_recep,MPI_INT,iproc,tag,
     $     commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              ntorque_send = ntorque_send + 1
              ptorque_send(ntorque_send) = iproc
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
      ntorqueshort_send = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(bufcount(iproc+1),1,MPI_INT,
     $     iproc,tag,commloc,reqrec(iproc+1),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(ntorqueshort_recep,1,MPI_INT,iproc,tag,
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
          call MPI_ISEND(ptorqueshort_recep,ntorqueshort_recep,MPI_INT,
     $     iproc,tag,commloc,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              ntorqueshort_send = ntorqueshort_send + 1
              ptorqueshort_send(ntorqueshort_send) = iproc
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
c      if (bigbuf.eq.torquebuf) then
c        nbig_recep = ntorque_recep
c        pbig_recep = ptorque_recep
c        nbig_send = ntorque_send
c        pbig_send = ptorque_send
c      else
c        nbig_recep = n_recep2
c        pbig_recep = p_recep2
c        nbig_send = n_send2
c        pbig_send = p_send2
c      end if
cc
c
 80   call orderbuffer_gpu(.false.)
      call orderbufferrec_gpu

c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (bufcount)
      return
      end
