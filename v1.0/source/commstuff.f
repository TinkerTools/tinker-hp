c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin

c
c     subroutine commforces : if necessary communicate some forces between processes
c
      subroutine commforces(derivs)
      implicit none
      include 'timestat.i'
      include 'mpif.h'
      real*8 derivs(3,*)
      real*8 time0,time1,time2
c
      time0 = mpi_wtime()
      call commforcesrec(derivs)
      time1 = mpi_wtime()
      timedirreccomm = timedirreccomm + time1-time0
      call commforcesbloc(derivs)
      time2 = mpi_wtime()
      timebondedvdw = timebondedvdw + time2-time1
      timeforcescomm = timeforcescomm + time2-time0
      return
      end
c
c     subroutine commforcesrespa : communicate some forces between processes
c     when respa integrator is used
c
      subroutine commforcesrespa(derivs,fast)
      implicit none
      include 'timestat.i'
      include 'mpif.h'
      real*8 derivs(3,*)
      real*8 time0,time1,time2
      logical fast
c
      if (fast) then
        time0 = mpi_wtime()
        call commforcesbonded(derivs)
        time1 = mpi_wtime()
        timebondedvdw = timebondedvdw + time1-time0
        timeforcescomm = timeforcescomm + time1-time0
      else
        time0 = mpi_wtime()
        call commforcesrec(derivs)
        time1 = mpi_wtime()
        timedirreccomm = timedirreccomm + time1-time0
        call commforcesbloc(derivs)
        time2 = mpi_wtime()
        timebondedvdw = timebondedvdw + time2-time1
        timeforcescomm = timeforcescomm + time2-time0
      end if
      return
      end
c
c     subroutine reduceen : make the necessary reductions over the processes to get the
c     total values of the energies
c
      subroutine reduceen(epot)
      implicit none
      include 'energi.i'
      include 'inter.i'
      include 'virial.i'
      include 'potent.i'
      include 'openmp.i'
      include 'mpif.h'
      integer ierr,commloc
      real*8 epot
      real*8 buffer1(5),buffer2(15)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = MPI_COMM_WORLD
      end if
c
      buffer1(1) = em
      buffer1(2) = ep
      buffer1(3) = epot
      buffer1(4) = esum
      buffer1(5) = einter 
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer1,5,MPI_REAL8,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      else
        call MPI_REDUCE(buffer1,buffer1,5,MPI_REAL8,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      end if
      em = buffer1(1)
      ep = buffer1(2)
      epot = buffer1(3)
      esum = buffer1(4)
      einter = buffer1(5)
c
c   MPI: get virial
c
      call MPI_ALLREDUCE(MPI_IN_PLACE,vir,9,MPI_REAL8,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)
c
      if ((use_pmecore).and.(rank.gt.(ndir-1))) return
c
      buffer2(1) = eba
      buffer2(2) = ea
      buffer2(3) = eb
      buffer2(4) = ev
      buffer2(5) = eub
      buffer2(6) = eaa
      buffer2(7) = eopb
      buffer2(8) = eopd
      buffer2(9) = eid
      buffer2(10) = eit
      buffer2(11) = et
      buffer2(12) = ept
      buffer2(13) = ebt
      buffer2(14) = ett
      buffer2(15) = eg
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer2,15,MPI_REAL8,MPI_SUM,0,
     $     commloc,ierr)
      else
        call MPI_REDUCE(buffer1,buffer2,15,MPI_REAL8,MPI_SUM,0,
     $     commloc,ierr)
      end if
c
      eba  = buffer2(1) 
      ea   = buffer2(2) 
      eb   = buffer2(3) 
      ev   = buffer2(4) 
      eub  = buffer2(5) 
      eaa  = buffer2(6) 
      eopb = buffer2(7) 
      eopd = buffer2(8) 
      eid  = buffer2(9) 
      eit  = buffer2(10)
      et   = buffer2(11)
      ept  = buffer2(12)
      ebt  = buffer2(13)
      ett  = buffer2(14)
      eg   = buffer2(15)
c
      return
      end
c
c     subroutine commforcesbonded : deal with communications of some
c     bonded forces modifier avec nlocnl ?
c
      subroutine commforcesbonded(derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'potent.i'
      include 'deriv.i'
      include 'mpif.h'
      integer i,j,tag,ierr,iproc
      integer iloc
      integer, allocatable :: req(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     allocate some arrays
c
      allocate (req(nproc*nproc))
      allocate (buffer(3,max(1,nloc),nneig_send))
      allocate (buffers(3,nbloc))
c      buffer = 0d0
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nneig_send
        tag = nproc*rank + pneig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pneig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      do iproc = 1, nneig_recep
        do i = 1, domlen(pneig_recep(iproc)+1)
          iloc = bufbeg(pneig_recep(iproc)+1)+i-1
          buffers(1,iloc) = debond(1,iloc)
          buffers(2,iloc) = debond(2,iloc)
          buffers(3,iloc) = debond(3,iloc)
        end do
      end do
c
c     MPI : begin sending
c
      do i = 1, nneig_recep
        tag = nproc*pneig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pneig_recep(i)+1)),
     $   3*domlen(pneig_recep(i)+1),
     $   MPI_REAL8,pneig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nneig_send
        tag = nproc*rank + pneig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nneig_recep
        tag = nproc*pneig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nneig_send
        do j = 1, nloc
          derivs(1,j) = derivs(1,j) + buffer(1,j,i)
          derivs(2,j) = derivs(2,j) + buffer(2,j,i)
          derivs(3,j) = derivs(3,j) + buffer(3,j,i)
        end do
      end do
c
      deallocate (req)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c
c     subroutine commforcesbloc : deal with communications of forces by bloc
c     between two time steps
c
      subroutine commforcesbloc(derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'openmp.i'
      include 'potent.i'
      include 'mpif.h'
      integer i,tag,ierr
      integer, allocatable :: req(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     allocate some arrays
c
      allocate (req(nproc*nproc))
      allocate (buffer(3,max(1,nloc),nbig_send))
c      buffer = 0d0
      allocate (buffers(3,nbloc))
c      buffers = 0d0
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = desum(:,(nloc+1):nbloc)
c
c     MPI : begin sending
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        derivs(:,1:nloc) = derivs(:,1:nloc) + buffer(:,1:nloc,i)
      end do
c
      deallocate (req)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commforcesrec : deal with communications of reciprocal forces
c
      subroutine commforcesrec(derivs)
      implicit none
      include 'sizes.i'
      include 'openmp.i'
      include 'deriv.i'
      include 'mpif.h'
      integer i,j,tag,ierr,iproc,iglob,ilocrec
      integer, allocatable :: req(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 derivs(3,*)
      integer status(MPI_STATUS_SIZE)
c
c     allocate some arrays
c
      allocate (req(nproc*nproc))
      allocate (buffer(3,max(1,nloc),nrecdir_send))
c      buffer = 0d0
      allocate (buffers(3,max(1,nblocrec)))
c      buffers = 0d0
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
        tag = nproc*rank + precdir_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,precdir_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
c
c     MPI : move in buffer
c
      buffers(:,(1):nblocrec)=torquerec(:,(1):nblocrec)+
     $           torquerecp(:,(1):nblocrec) +
     $  demrec(:,(1):nblocrec) + deprec(:,(1):nblocrec)
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
        tag = nproc*precdir_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $   3*domlen(precdir_recep(i)+1),
     $   MPI_REAL8,precdir_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
        end if
      end do
c
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
        tag = nproc*rank + precdir_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
        tag = nproc*precdir_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
          do j = 1, nloc
            derivs(1,j) = derivs(1,j) + buffer(1,j,i)
            derivs(2,j) = derivs(2,j) + buffer(2,j,i)
            derivs(3,j) = derivs(3,j) + buffer(3,j,i)
          end do
        else
          do j = 1, nloc
            iglob = glob(j)
            ilocrec = locrec1(iglob)
            derivs(1,j) = derivs(1,j) + torquerec(1,ilocrec) +
     $        torquerecp(1,ilocrec) + demrec(1,ilocrec)
     $      + deprec(1,ilocrec)
            derivs(2,j) = derivs(2,j) + torquerec(2,ilocrec) +
     $        torquerecp(2,ilocrec) + demrec(2,ilocrec)
     $      + deprec(2,ilocrec)
            derivs(3,j) = derivs(3,j) + torquerec(3,ilocrec) +
     $        torquerecp(3,ilocrec) + demrec(3,ilocrec)
     $      + deprec(3,ilocrec)
          end do
        end if
      end do
c
c    Communicate separately the direct and reciprocal torques if testgrad is
c    being used
c
      if (dotstgrad) then
c
c     MPI : begin reception in buffer
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $       MPI_REAL8,precdir_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep(iproc)+1)
              ilocrec = bufbegrec(precdir_recep(iproc)+1)+i-1
              buffers(1,ilocrec) =demrec(1,ilocrec)+torquerec(1,ilocrec)
              buffers(2,ilocrec) =demrec(2,ilocrec)+torquerec(2,ilocrec)
              buffers(3,ilocrec) =demrec(3,ilocrec)+torquerec(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $       3*domlen(precdir_recep(i)+1),
     $       MPI_REAL8,precdir_recep(i),tag,MPI_COMM_WORLD,
     $       req(tag),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
c
c       MPI : move in global arrays
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            do j = 1, nloc
              iglob = glob(j)
              dem(1,j) = dem(1,j) + buffer(1,j,i)
              dem(2,j) = dem(2,j) + buffer(2,j,i)
              dem(3,j) = dem(3,j) + buffer(3,j,i)
            end do
          else
            do j = 1, nloc
              iglob = glob(j)
              ilocrec = locrec1(iglob)
              dem(1,j) = dem(1,j)+demrec(1,ilocrec)+torquerec(1,ilocrec)
              dem(2,j) = dem(2,j)+demrec(2,ilocrec)+torquerec(2,ilocrec)
              dem(3,j) = dem(3,j)+demrec(3,ilocrec)+torquerec(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $       MPI_REAL8,precdir_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep(iproc)+1)
              ilocrec = bufbegrec(precdir_recep(iproc)+1)+i-1
              buffers(1,ilocrec)=deprec(1,ilocrec)+torquerecp(1,ilocrec)
              buffers(2,ilocrec)=deprec(2,ilocrec)+torquerecp(2,ilocrec)
              buffers(3,ilocrec)=deprec(3,ilocrec)+torquerecp(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $       3*domlen(precdir_recep(i)+1),
     $       MPI_REAL8,precdir_recep(i),tag,MPI_COMM_WORLD,
     $       req(tag),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
c
c       MPI : move in global arrays
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            do j = 1, nloc
              iglob = glob(j)
              dep(1,j) = dep(1,j) + buffer(1,j,i)
              dep(2,j) = dep(2,j) + buffer(2,j,i)
              dep(3,j) = dep(3,j) + buffer(3,j,i)
            end do
          else
            do j = 1, nloc
              iglob = glob(j)
              ilocrec = locrec1(iglob)
              dep(1,j) =dep(1,j)+deprec(1,ilocrec)+torquerecp(1,ilocrec)
              dep(2,j) =dep(2,j)+deprec(2,ilocrec)+torquerecp(2,ilocrec)
              dep(3,j) =dep(3,j)+deprec(3,ilocrec)+torquerecp(3,ilocrec)
            end do
          end if
        end do
      end  if
      deallocate (req)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine ddpme: domain decomposition load balancing
c     assign atom sites to MPI processes based on a domain decomposition
c
c
      subroutine ddpme3d
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'cell.i'
      include 'cutoff.i'
      include 'openmp.i'
      include 'mpif.h'
      include 'neigh.i'
      include 'potent.i'
      integer nx, ny, nz, i,j,k
      integer nprocloc,rankloc,iproc,iproc1
      integer it,nitmax
      real*8  ndiff
      real*8  exchgsizex
      real*8 xr,yr,zr
      real*8 x2,x3,y1,y3,z1,z3
      real*8 dist
      real*8 mbuf,vbuf,torquebuf,neigbuf,bigbuf
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
      integer, allocatable :: ntorquerecep(:),ptorquerecep(:,:)
      integer, allocatable :: nneigrecep(:),pneigrecep(:,:)
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
      allocate (ntorquerecep(nproc))
      allocate (ptorquerecep(nproc,nproc))
      allocate (nneigrecep(nproc))
      allocate (pneigrecep(nproc,nproc))
      nrecep1 = 0
      nrecep2 = 0
      ntorquerecep = 0
      nneigrecep = 0
      precep1 = 0
      precep2 = 0
      ptorquerecep = 0
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
      repart = 0
      domlen = 0
      glob = 0
      loc = 0
c
c     get the number of subdivision along each axis for dd
c
      call ddnumber(nprocloc,nx,ny,nz,0)
c
      nx_box = xbox/nx
      ny_box = ybox/ny
      nz_box = zbox/nz
      do i = 0, nx-1
        xbegproctemp(i+1) = -xcell2 + i*nx_box
        xendproctemp(i+1) = -xcell2 + (i+1)*nx_box
      end do
      do i = 0, ny-1
        ybegproctemp(i+1) = -ycell2 + i*ny_box
        yendproctemp(i+1) = -ycell2 + (i+1)*ny_box
      end do
      do i = 0, nz-1
        zbegproctemp(i+1) = -zcell2 + i*nz_box
        zendproctemp(i+1) = -zcell2 + (i+1)*nz_box
      end do
c
c     assign processes
c
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
              iproc = (k-1)*ny*nx+(j-1)*nx+i
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
        if (abs(xr-xcell2).lt.eps1) xr = xr-eps2
        if (abs(yr-ycell2).lt.eps1) yr = yr-eps2
        if (abs(zr-zcell2).lt.eps1) zr = zr-eps2
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
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            iproc = (k-1)*ny*nx+(j-1)*nx+i
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
                  if ((i.eq.1).and.(p.eq.-1)) temp_x = nx
                  if ((i.eq.nx).and.(p.eq.1)) temp_x = 1
                  if ((j.eq.1).and.(q.eq.-1)) temp_y = ny-1
                  if ((j.eq.ny).and.(q.eq.1)) temp_y = 0
                  if ((k.eq.1).and.(r.eq.-1)) temp_z = nz-1
                  if ((k.eq.nz).and.(r.eq.1)) temp_z = 0
                  tempproc = temp_z*ny*nx+temp_y*nx+
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
c     change size with neighbors
c
      exchgsizex = nx_box/20
      nitmax = 3 
      ndiff  = n/(20*nprocloc)
c
      do it = 1, nitmax
        do iproc = 1, nprocloc
          x2 = xendproc(iproc)
          y1 = ybegproc(iproc)
          z1 = zbegproc(iproc)
c
c    check neighbors
c
          do i = 1, numneigproc(iproc)
            iproc1 = neigproc(i,iproc)
            x3 = xbegproc(iproc1)
            y3 = ybegproc(iproc1)
            z3 = zbegproc(iproc1)
c
            if ((x2.eq.x3).and.(y1.eq.y3).and.(z1.eq.z3)) then
              if ((domlen(iproc1).ge.domlen(iproc)+ndiff)) then
                xendproc(iproc)   = xendproc(iproc) +
     $            exchgsizex
                xbegproc(iproc1)   = xbegproc(iproc1) + exchgsizex
              end if
            end if
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
          if (abs(xr-xcell2).lt.eps1) xr = xr-eps2
          if (abs(yr-ycell2).lt.eps1) yr = yr-eps2
          if (abs(zr-zcell2).lt.eps1) zr = zr-eps2
          do iproc = 0, nprocloc-1
            if ((zr.ge.zbegproc(iproc+1)).and.
     $        (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
     $        .and.(yr.lt.yendproc(iproc+1))
     $        .and.(xr.ge.xbegproc(iproc+1))
     $        .and.(xr.lt.xendproc(iproc+1))) then
               repart(i) = iproc
              domlen(repart(i)+1) = domlen(repart(i)+1) + 1
            end if
          end do
        end do
      end do
c
      if ((use_pmecore).and.(rank.le.ndir-1)) then
        nloc = domlen(rankloc+1)
      else if ((use_pmecore).and.(rank.gt.ndir-1)) then
        nloc = 0
      else
        nloc = domlen(rankloc+1)
      end if
c
c     get the processes to receive data from : for commstep and for commtorque
c
      p_send1 = 0
      p_recep1 = 0
      p_send2 = 0
      p_recep2 = 0
      ptorque_send = 0
      ptorque_recep = 0
      pbig_send = 0
      pbig_recep = 0
      n_recep1 = 0
      n_recep2 = 0
      ntorque_recep = 0
      nneig_recep = 0
      nbig_recep = 0
      if ((use_pmecore).and.(rank.gt.ndir-1)) then
        goto 80
      end if
c
      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2)+2.0d0
      torquebuf = mbuf + lbuffer
      torquebuf2 = torquebuf*torquebuf 
      neigbuf = lbuffer
c
c     get maximum cutoff value
c
      if (torquebuf.gt.(vbuf)) then
        bigbuf = torquebuf
      else
        bigbuf = vbuf
      end if
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
          if (dist.le.(torquebuf/2)) then
            ntorquerecep(iproc1+1) = ntorquerecep(iproc1+1)+1
            ptorquerecep(ntorquerecep(iproc1+1),iproc1+1) = iproc
          end if
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
      if (bigbuf.eq.torquebuf) then
        nbig_recep = ntorque_recep
        pbig_recep = ptorque_recep
        nbig_send = ntorque_send
        pbig_send = ptorque_send
      else
        nbig_recep = n_recep2
        pbig_recep = p_recep2
        nbig_send = n_send2
        pbig_send = p_send2
      end if
c
 80   call orderbuffer(.true.,domlen)
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
      deallocate (ntorquerecep)
      deallocate (ptorquerecep)
      deallocate (nneigrecep)
      deallocate (pneigrecep)
      return
      end
c
c     subroutine commforceststgrad : deal with communications of forces one by one
c     for testgrad
c
      subroutine commforceststgrad(detot)
      implicit none
      include 'sizes.i'
      include 'openmp.i'
      include 'deriv.i'
      include 'mpif.h'
      integer i,j,k,tag,ierr
      integer, allocatable :: req(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 detot(3,*)
      integer status(MPI_STATUS_SIZE)
c
c     allocate some arrays
c
      allocate (req(nproc*nproc))
c
      allocate (buffer(3,max(1,nloc),nbig_send))
      buffer = 0d0
      allocate (buffers(3,nbloc))
      buffers = 0d0
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deb(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deb(1,j) = deb(1,j) + buffer(1,j,i)
          deb(2,j) = deb(2,j) + buffer(2,j,i)
          deb(3,j) = deb(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dea(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dea(1,j) = dea(1,j) + buffer(1,j,i)
          dea(2,j) = dea(2,j) + buffer(2,j,i)
          dea(3,j) = dea(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deba(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deba(1,j) = deba(1,j) + buffer(1,j,i)
          deba(2,j) = deba(2,j) + buffer(2,j,i)
          deba(3,j) = deba(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deub(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deub(1,j) = deub(1,j) + buffer(1,j,i)
          deub(2,j) = deub(2,j) + buffer(2,j,i)
          deub(3,j) = deub(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deaa(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deaa(1,j) = deaa(1,j) + buffer(1,j,i)
          deaa(2,j) = deaa(2,j) + buffer(2,j,i)
          deaa(3,j) = deaa(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deopb(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deopb(1,j) = deopb(1,j) + buffer(1,j,i)
          deopb(2,j) = deopb(2,j) + buffer(2,j,i)
          deopb(3,j) = deopb(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deopd(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deopd(1,j) = deopd(1,j) + buffer(1,j,i)
          deopd(2,j) = deopd(2,j) + buffer(2,j,i)
          deopd(3,j) = deopd(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deid(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deid(1,j) = deid(1,j) + buffer(1,j,i)
          deid(2,j) = deid(2,j) + buffer(2,j,i)
          deid(3,j) = deid(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deit(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deit(1,j) = deit(1,j) + buffer(1,j,i)
          deit(2,j) = deit(2,j) + buffer(2,j,i)
          deit(3,j) = deit(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = det(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          det(1,j) = det(1,j) + buffer(1,j,i)
          det(2,j) = det(2,j) + buffer(2,j,i)
          det(3,j) = det(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dept(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dept(1,j) = dept(1,j) + buffer(1,j,i)
          dept(2,j) = dept(2,j) + buffer(2,j,i)
          dept(3,j) = dept(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = debt(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          debt(1,j) = debt(1,j) + buffer(1,j,i)
          debt(2,j) = debt(2,j) + buffer(2,j,i)
          debt(3,j) = debt(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dett(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dett(1,j) = dett(1,j) + buffer(1,j,i)
          dett(2,j) = dett(2,j) + buffer(2,j,i)
          dett(3,j) = dett(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dev(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dev(1,j) = dev(1,j) + buffer(1,j,i)
          dev(2,j) = dev(2,j) + buffer(2,j,i)
          dev(3,j) = dev(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deg(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deg(1,j) = deg(1,j) + buffer(1,j,i)
          deg(2,j) = deg(2,j) + buffer(2,j,i)
          deg(3,j) = deg(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = torquedir(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dem(1,j) = dem(1,j) + buffer(1,j,i)
          dem(2,j) = dem(2,j) + buffer(2,j,i)
          dem(3,j) = dem(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = torquedirp(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dep(1,j) = dep(1,j) + buffer(1,j,i)
          dep(2,j) = dep(2,j) + buffer(2,j,i)
          dep(3,j) = dep(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dem(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dem(1,j) = dem(1,j) + buffer(1,j,i)
          dem(2,j) = dem(2,j) + buffer(2,j,i)
          dem(3,j) = dem(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dep(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dep(1,j) = dep(1,j) + buffer(1,j,i)
          dep(2,j) = dep(2,j) + buffer(2,j,i)
          dep(3,j) = dep(3,j) + buffer(3,j,i)
          detot(1,j) = detot(1,j) + buffer(1,j,i)
          detot(2,j) = detot(2,j) + buffer(2,j,i)
          detot(3,j) = detot(3,j) + buffer(3,j,i)
        end do
      end do
      deallocate (req)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine allreduceen : make the necessary allreductions over the processes to get the
c     total values of the energies
c
      subroutine allreduceen(epot)
      implicit none
      include 'energi.i'
      include 'inter.i'
      include 'mpif.h'
      integer ierr
      real*8 epot
c
       call MPI_ALLREDUCE(MPI_IN_PLACE,eba,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ea,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eb,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,em,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ep,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ev,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eub,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eaa,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eopb,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eopd,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eid,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eit,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,et,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ept,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ebt,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ett,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,esum,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,epot,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eg,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,einter,1,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,ierr)
      return
      end
c
c     subroutine distproc : get the minimum distance between two 3d domains
c
c     nb : to be modified for non cubic unit cells
      subroutine distproc(iproc1,iproc2,dist,do3d)
      implicit none
      include 'openmp.i'
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
c     subroutine distproc : get the minimum distance between a 3d domain and an atom
c
c     to be modified for non cubic unit cells
      subroutine distprocpart(i,iproc,dist,do3d)
      implicit none
      include 'atoms.i'
      include 'openmp.i'
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
c
c     subroutine orderbuffer : get the indexes of the particules in the neighboring processes
c
      subroutine orderbuffer(init,domlenold)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'cutoff.i'
      include 'openmp.i'
      include 'mpif.h'
      integer, allocatable :: counttemp(:),ind1temp(:)
      integer, allocatable :: globtemp(:)
      integer domlenold(*)
      integer count
      integer i,iproc,iproc1,tag,req(nproc*nproc),ierr
      integer iloc,iglob
      integer status(MPI_STATUS_SIZE)
      logical init
c
      allocate (counttemp(nproc))
      allocate (ind1temp(n))
      allocate (globtemp(n))
      globtemp = glob
c
c      ind1temp = 0
      counttemp = 0
c
c     get the domlen values of the coupled real space processes
c
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
        tag = nproc*rank + precdir_recep(iproc) + 1
        call MPI_IRECV(domlen(precdir_recep(iproc)+1),1,MPI_INT,
     $   precdir_recep(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
        tag = nproc*precdir_send(iproc) + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,precdir_send(iproc),
     $   tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
        tag = nproc*rank + precdir_recep(iproc) + 1
        call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
        tag = nproc*precdir_send(iproc) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     get the domlen values of the coupled real space processes
c
      do iproc = 1, nbig_recep
        tag = nproc*rank + pbig_recep(iproc) + 1
        call MPI_IRECV(domlen(pbig_recep(iproc)+1),1,MPI_INT,
     $   pbig_recep(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, nbig_send
        tag = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,pbig_send(iproc),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, nbig_recep
        tag = nproc*rank + pbig_recep(iproc) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do iproc = 1, nbig_send
        tag = nproc*pbig_send(iproc) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
      if (init) then
        do i = 1, n
c
c       check if atom is in the domain or in one of the neighboring domains
c
          counttemp(repart(i)+1) = counttemp(repart(i)+1) + 1
          ind1temp(i) = counttemp(repart(i)+1)
          if (repart(i).eq.rank) then
            glob(ind1temp(i)) = i
            loc(i) = ind1temp(i)
          end if
        end do
      else
        do i = 1, nloc
          iglob = globtemp(i)
          counttemp(repart(iglob)+1) = counttemp(repart(iglob)+1) + 1
          ind1temp(iglob) = counttemp(repart(iglob)+1)
          if (repart(iglob).eq.rank) then
            glob(ind1temp(iglob)) = iglob
            loc(iglob) = ind1temp(iglob)
          end if
        end do
        do iproc = 1, nneig_recep
          do i = 1, domlenold(pneig_recep(iproc)+1)
            iloc = bufbeg(pneig_recep(iproc)+1)+i-1
            iglob = globtemp(iloc)
            counttemp(repart(iglob)+1) = counttemp(repart(iglob)+1) + 1
            ind1temp(iglob) = counttemp(repart(iglob)+1)
            if (repart(iglob).eq.rank) then
              glob(ind1temp(iglob)) = iglob
              loc(iglob) = ind1temp(iglob)
            end if
          end do
        end do
      end if
c
      bufbeg(rank+1) = 1
      count = domlen(rank+1)
      do iproc = 1, nbig_recep
        if (domlen(pbig_recep(iproc)+1).ne.0) then
          bufbeg(pbig_recep(iproc)+1) = count + 1
        else
          bufbeg(pbig_recep(iproc)+1) = 1
        end if
        count = count + domlen(pbig_recep(iproc)+1)
      end do
      nbloc = count
      nblocrecdir = count
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).eq.rank) goto 10
        do iproc1 = 1, nbig_recep
          if (pbig_recep(iproc1).eq.precdir_recep(iproc)) goto 10
        end do
        if (domlen(precdir_recep(iproc)+1).ne.0) then
          bufbeg(precdir_recep(iproc)+1) = count + 1
        else
          bufbeg(precdir_recep(iproc)+1) = 1
        end if
        count = count + domlen(precdir_recep(iproc)+1)
  10    continue
      end do
      nblocrecdir = count
c
c     send and receive the indexes of the neighboring processes
c
      do iproc = 1, nbig_recep
        tag = nproc*rank + pbig_recep(iproc) + 1
        call MPI_IRECV(glob(bufbeg(pbig_recep(iproc)+1)),
     $   domlen(pbig_recep(iproc)+1),MPI_INT,pbig_recep(iproc),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, nbig_send
        tag = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,pbig_send(iproc),
     $   tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, nbig_recep
        tag = nproc*rank + pbig_recep(iproc) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do iproc = 1, nbig_send
        tag = nproc*pbig_send(iproc) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
      do iproc = 1, nbig_recep
        do i = 1, domlen(pbig_recep(iproc)+1)
          loc(glob(bufbeg(pbig_recep(iproc)+1)+i-1)) =
     $      bufbeg(pbig_recep(iproc)+1)+i-1
        end do
      end do
c
c     also get the indexes of recdir neighboring processes
c
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
          tag = nproc*rank + precdir_recep(iproc) + 1
          call MPI_IRECV(glob(bufbeg(precdir_recep(iproc)+1)),
     $     domlen(precdir_recep(iproc)+1),MPI_INT,
     $     precdir_recep(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
          tag = nproc*precdir_send(iproc) + rank + 1
          call MPI_ISEND(glob,domlen(rank+1),MPI_INT,
     $     precdir_send(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
          tag = nproc*rank + precdir_recep(iproc) + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
          tag = nproc*precdir_send(iproc) + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
          do i = 1, domlen(precdir_recep(iproc)+1)
            loc(glob(bufbeg(precdir_recep(iproc)+1)+i-1)) =
     $        bufbeg(precdir_recep(iproc)+1)+i-1
          end do
        end if
      end do
c
      deallocate (globtemp)
      deallocate (ind1temp)
      deallocate (counttemp)
      return
      end
c
c
c     subroutine orderbufferrec : get the indexes of the particules in the reciprocal processes
c
      subroutine orderbufferrec
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'potent.i'
      integer, allocatable :: counttemp(:),ind1temp(:)
      integer count,rankloc
      integer i,iproc
      integer iloc,iglob,iproc1
c
      allocate (counttemp(nproc))
      allocate (ind1temp(n))
c
      if (use_pmecore) then
        rankloc = rank_bis
      else
        rankloc = rank
      end if
c
c      ind1temp = 0 
      counttemp = 0
      domlenrec = 0
c
      if ((.not.(use_pmecore)).or.((use_pmecore).and.(rank.gt.ndir-1)))
     $ then
        do iproc = 1, nrecdir_recep
          do i = 1, domlen(precdir_recep(iproc)+1)
            iloc = bufbeg(precdir_recep(iproc)+1)+i-1
            iglob = glob(iloc)
c
c       check if atom is in the domain or in one of the neighboring domains
c
            counttemp(repartrec(iglob)+1) =
     $        counttemp(repartrec(iglob)+1) + 1
            ind1temp(iglob) = counttemp(repartrec(iglob)+1)
            domlenrec(repartrec(iglob)+1)=
     $       domlenrec(repartrec(iglob)+1) + 1
            if (repartrec(iglob).eq.rankloc) then
              globrec(ind1temp(iglob)) = iglob
              locrec(iglob) = ind1temp(iglob)
            end if
          end do
        end do
      end if
c
      if (use_pmecore) then
        nlocrec = domlenrec(rank_bis+1)
      else
        nlocrec = domlenrec(rank+1)
      end if
c
      count = 0
      bufbegrec(precdir_recep(1)+1) = 1
      if (nrecdir_recep.gt.0) then
        count = domlen(precdir_recep(1)+1)
      end if
      do iproc = 2, nrecdir_recep
        if (domlen(precdir_recep(iproc)+1).ne.0) then
          bufbegrec(precdir_recep(iproc)+1) = count + 1
        else
          bufbegrec(precdir_recep(iproc)+1) = 1
        end if
        count = count + domlen(precdir_recep(iproc)+1)
      end do
      nblocrec = count
c
      do iproc = 1, nrecdir_recep
        do i = 1, domlen(precdir_recep(iproc)+1)
          globrec1(bufbegrec(precdir_recep(iproc)+1)+i-1) =
     $       glob(bufbeg(precdir_recep(iproc)+1)+i-1)
          locrec1(globrec1(bufbegrec(precdir_recep(iproc)+1)+i-1)) =
     $      bufbegrec(precdir_recep(iproc)+1)+i-1
        end do
      end do
c
      deallocate (ind1temp)
      deallocate (counttemp)
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
      implicit none
      include 'bound.i'
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
c     has to be computed within the current domain or not (dd midpoint method, Newton's
c     3rd law)
c
c
      subroutine midpoint(xi,yi,zi,xk,yk,zk,docompute)
      implicit none
      include 'openmp.i'
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 xrmid,yrmid,zrmid
      logical docompute
c
      docompute = .false.
c
      call image(xi,yi,zi)
      call image(xk,yk,zk)
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
c     subroutine midpointgroup : routine that says whether an interaction between a
c     group of particules
c     has to be computed within the current domain or not (dd midpoint method, Newton's
c     3rd law)
c
c
      subroutine midpointgroup(pos,nb,docompute)
      implicit none
      include 'openmp.i'
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
      subroutine ddpme3dnpt(scale,istep)
      implicit none
      include 'openmp.i'
      include 'mpif.h'
      include 'neigh.i'
      include 'potent.i'
      integer nprocloc,rankloc,commloc,iproc
      integer i,status(MPI_STATUS_SIZE),tag,ierr
      integer istep,modnl
      real*8 mbuf,vbuf,torquebuf,neigbuf,bigbuf
      real*8 dist, scale
      integer, allocatable :: bufcount(:),buffer(:,:),req(:)
cc
      modnl = mod(istep,ineigup)
cc
      if (use_pmecore) then
        nprocloc = ndir
        commloc  = comm_dir
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = MPI_COMM_WORLD
        rankloc  = rank
      end if
c
      do iproc = 1, nprocloc
        xbegproc(iproc) = scale*xbegproc(iproc)
        xendproc(iproc) = scale*xendproc(iproc)
        ybegproc(iproc) = scale*ybegproc(iproc)
        yendproc(iproc) = scale*yendproc(iproc)
        zbegproc(iproc) = scale*zbegproc(iproc)
        zendproc(iproc) = scale*zendproc(iproc)
      end do
      if (modnl.ne.0) return
c
      allocate (req(nprocloc*nprocloc))
      allocate (bufcount(nprocloc))
      allocate (buffer(nprocloc,nprocloc))
c
c     get the processes to receive data from : for commstep and for commtorque
c
      p_send1 = 0
      p_recep1 = 0
      p_send2 = 0
      p_recep2 = 0
      ptorque_send = 0
      ptorque_recep = 0
      pbig_send = 0
      pbig_recep = 0
      n_recep1 = 0
      n_recep2 = 0
      ntorque_recep = 0
      nneig_recep = 0
      nbig_recep = 0
      if ((use_pmecore).and.(rank.gt.ndir-1)) then
        goto 80
      end if
c
      mbuf = sqrt(mbuf2)
      vbuf = sqrt(vbuf2)+2.0d0
      torquebuf = mbuf + lbuffer
      torquebuf2 = torquebuf*torquebuf
      neigbuf = lbuffer
c
c     get maximum cutoff value
c
      if (torquebuf.gt.vbuf) then
        bigbuf = torquebuf
      else
        bigbuf = vbuf
      end if
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
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(n_recep1,1,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(p_recep1,n_recep1,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
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
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(ntorque_recep,1,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(ptorque_recep,ntorque_recep,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
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
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(n_recep2,1,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(p_recep2,n_recep2,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
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
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(nneig_recep,1,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     receive the buffers and get the information to send
c
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_IRECV(buffer(1,iproc+1),bufcount(iproc+1),MPI_INT,
     $     iproc,tag,commloc,req(tag),ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_ISEND(pneig_recep,nneig_recep,MPI_INT,iproc,tag,
     $     commloc,req(tag),ierr)
        end if
      end do
      do iproc = 0, nprocloc-1
        if (iproc.ne.rankloc) then
          tag = nprocloc*rankloc + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nprocloc*iproc + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
          do i = 1, bufcount(iproc+1)
            if (buffer(i,iproc+1).eq.rankloc) then
              nneig_send = nneig_send + 1
              pneig_send(nneig_send) = iproc
            end if
          end do
        end if
      end do
      if (bigbuf.eq.torquebuf) then
        nbig_recep = ntorque_recep
        pbig_recep = ptorque_recep
        nbig_send = ntorque_send
        pbig_send = ptorque_send
      else
        nbig_recep = n_recep2
        pbig_recep = p_recep2
        nbig_send = n_send2
        pbig_send = p_send2
      end if
cc
c
 80   call orderbuffer(.false.,domlen)
      call orderbufferrec

c
      deallocate (req)
      deallocate (buffer)
      deallocate (bufcount)
      return
      end
c
c
c     subroutine orderbufferrespa : get the indexes of the particules in the neighboring processes with respa integrator
c
      subroutine orderbufferrespa(init,fast,domlenold)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'mpif.h'
      integer, allocatable :: counttemp(:),ind1temp(:)
      integer, allocatable :: globtemp(:)
      integer domlenold(*)
      integer count
      integer i,iproc,iproc1,tag,req(nproc*nproc),ierr
      integer iloc,iglob
      integer status(MPI_STATUS_SIZE)
      integer, allocatable :: ptemp_recep(:),ptemp_send(:)
      integer ntemp_recep,ntemp_send
      logical init,fast
c
      allocate (ptemp_recep(nproc))
      allocate (ptemp_send(nproc))
      if (fast) then
        ptemp_recep = pneig_recep
        ptemp_send = pneig_send
        ntemp_recep = nneig_recep
        ntemp_send = nneig_send
      else
        ptemp_recep = pbig_recep
        ptemp_send = pbig_send
        ntemp_recep = nbig_recep
        ntemp_send = nbig_send
      end if
c
      allocate (counttemp(nproc))
      allocate (ind1temp(n))
      allocate (globtemp(n))
      globtemp = glob
c
c      ind1temp = 0
      counttemp = 0
c
c     get the domlen values of the coupled real space processes
c
      if ((.not.(fast))) then
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
          tag = nproc*rank + precdir_recep(iproc) + 1
          call MPI_IRECV(domlen(precdir_recep(iproc)+1),1,MPI_INT,
     $     precdir_recep(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
          end if
        end do
        do iproc = 1, nrecdir_send
          if (precdir_send(iproc).ne.rank) then
          tag = nproc*precdir_send(iproc) + rank + 1
          call MPI_ISEND(domlen(rank+1),1,MPI_INT,precdir_send(iproc),
     $     tag,MPI_COMM_WORLD,req(tag),ierr)
          end if
        end do
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
          tag = nproc*rank + precdir_recep(iproc) + 1
          call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
        do iproc = 1, nrecdir_send
          if (precdir_send(iproc).ne.rank) then
          tag = nproc*precdir_send(iproc) + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
      end if
c
c     get the domlen values of the coupled real space processes
c
      do iproc = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(iproc) + 1
        call MPI_IRECV(domlen(ptemp_recep(iproc)+1),1,MPI_INT,
     $   ptemp_recep(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, ntemp_send
        tag = nproc*ptemp_send(iproc) + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,ptemp_send(iproc),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(iproc) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do iproc = 1, ntemp_send
        tag = nproc*ptemp_send(iproc) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
      if (init) then
        do i = 1, n
c
c       check if atom is in the domain or in one of the neighboring domains
c
          counttemp(repart(i)+1) = counttemp(repart(i)+1) + 1
          ind1temp(i) = counttemp(repart(i)+1)
          if (repart(i).eq.rank) then
            glob(ind1temp(i)) = i
            loc(i) = ind1temp(i)
          end if
        end do
      else
        do i = 1, nloc
          iglob = globtemp(i)
          counttemp(repart(iglob)+1) = counttemp(repart(iglob)+1) + 1
          ind1temp(iglob) = counttemp(repart(iglob)+1)
          if (repart(iglob).eq.rank) then
            glob(ind1temp(iglob)) = iglob
            loc(iglob) = ind1temp(iglob)
          end if
        end do
        do iproc = 1, nneig_recep
          do i = 1, domlenold(pneig_recep(iproc)+1)
            iloc = bufbeg(pneig_recep(iproc)+1)+i-1
            iglob = globtemp(iloc)
            counttemp(repart(iglob)+1) = counttemp(repart(iglob)+1) + 1
            ind1temp(iglob) = counttemp(repart(iglob)+1)
            if (repart(iglob).eq.rank) then
              glob(ind1temp(iglob)) = iglob
              loc(iglob) = ind1temp(iglob)
            end if
          end do
        end do
      end if
c
      bufbeg(rank+1) = 1
      count = domlen(rank+1)
      do iproc = 1, ntemp_recep
        if (domlen(ptemp_recep(iproc)+1).ne.0) then
          bufbeg(ptemp_recep(iproc)+1) = count + 1
        else
          bufbeg(ptemp_recep(iproc)+1) = 1
        end if
        bufbeg(ptemp_recep(iproc)+1) = count + 1
        count = count + domlen(ptemp_recep(iproc)+1)
      end do
      nbloc = count
      if ((.not.(fast))) then
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).eq.rank) goto 10
          do iproc1 = 1, nbig_recep
            if (pbig_recep(iproc1).eq.precdir_recep(iproc)) goto 10
          end do
          if (domlen(precdir_recep(iproc)+1).ne.0) then
            bufbeg(precdir_recep(iproc)+1) = count + 1
          else
            bufbeg(precdir_recep(iproc)+1) = 1
          end if
          count = count + domlen(precdir_recep(iproc)+1)
  10      continue
        end do
        nblocrecdir = count
      end if
c
c     send and receive the indexes of the neighboring processes
c
      do iproc = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(iproc) + 1
        call MPI_IRECV(glob(bufbeg(ptemp_recep(iproc)+1)),
     $   domlen(ptemp_recep(iproc)+1),MPI_INT,ptemp_recep(iproc),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, ntemp_send
        tag = nproc*ptemp_send(iproc) + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,ptemp_send(iproc),
     $   tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
      do iproc = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(iproc) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do iproc = 1, ntemp_send
        tag = nproc*ptemp_send(iproc) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do

      do iproc = 1, ntemp_recep
        do i = 1, domlen(ptemp_recep(iproc)+1)
          loc(glob(bufbeg(ptemp_recep(iproc)+1)+i-1)) =
     $      bufbeg(ptemp_recep(iproc)+1)+i-1
        end do
      end do
c
c     also get the indexes of recdir neighboring processes
c
      if ((.not.(fast))) then
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            tag = nproc*rank + precdir_recep(iproc) + 1
            call MPI_IRECV(glob(bufbeg(precdir_recep(iproc)+1)),
     $       domlen(precdir_recep(iproc)+1),MPI_INT,
     $       precdir_recep(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
          end if
        end do
        do iproc = 1, nrecdir_send
          if (precdir_send(iproc).ne.rank) then
            tag = nproc*precdir_send(iproc) + rank + 1
            call MPI_ISEND(glob,domlen(rank+1),MPI_INT,
     $       precdir_send(iproc),tag,MPI_COMM_WORLD,req(tag),ierr)
          end if
        end do
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            tag = nproc*rank + precdir_recep(iproc) + 1
            call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
        do iproc = 1, nrecdir_send
          if (precdir_send(iproc).ne.rank) then
            tag = nproc*precdir_send(iproc) + rank + 1
            call MPI_WAIT(req(tag),status,ierr)
          end if
        end do

        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep(iproc)+1)
              loc(glob(bufbeg(precdir_recep(iproc)+1)+i-1)) =
     $          bufbeg(precdir_recep(iproc)+1)+i-1
            end do
          end if
        end do
      end if
c
      deallocate (ptemp_recep)
      deallocate (ptemp_send)
      deallocate (globtemp)
      deallocate (ind1temp)
      deallocate (counttemp)
      return
      end
