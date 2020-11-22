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
c
      subroutine drivermpi
      use atoms
      use domdec
      use iounit
      use potent
      use mpi
      implicit none
      integer iproc, ierr
      integer total_group, direct_group, rec_group
      integer, allocatable :: direct_rank(:)
c
      ndir = nproc - nrec
c
c     MPI : get the atoms repartition over the processes
c
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
        if (.not. allocated(glob)) allocate (glob(n))
        if (.not. allocated(loc)) allocate (loc(n))
        if (.not. allocated(globrec)) allocate (globrec(n))
        if (.not. allocated(locrec)) allocate (locrec(n))
        if (.not. allocated(globrec1)) allocate (globrec1(n))
        if (.not. allocated(locrec1)) allocate (locrec1(n))
        if (.not.allocated(repart)) allocate (repart(n))
        if (.not.allocated(domlen)) allocate (domlen(nproc))
        if (.not.allocated(domlenpole)) allocate (domlenpole(nproc))
        if (.not.allocated(domlenpolerec)) 
     $     allocate (domlenpolerec(nproc))
        if (.not.allocated(zbegproc)) allocate (zbegproc(nproc))
        if (.not.allocated(zendproc)) allocate (zendproc(nproc))
        if (.not.allocated(ybegproc)) allocate (ybegproc(nproc))
        if (.not.allocated(yendproc)) allocate (yendproc(nproc))
        if (.not.allocated(xbegproc)) allocate (xbegproc(nproc))
        if (.not.allocated(xendproc)) allocate (xendproc(nproc))
        if (.not.allocated(p_recep1))   allocate(p_recep1(nproc))
        if (.not.allocated(p_send1))  allocate(p_send1(nproc))
        if (.not.allocated(p_recep2))   allocate(p_recep2(nproc))
        if (.not.allocated(p_send2))  allocate(p_send2(nproc))
        if (.not.allocated(pneig_recep))   allocate(pneig_recep(nproc))
        if (.not.allocated(pneig_send))  allocate(pneig_send(nproc))
        if (.not.allocated(precdir_recep))
     $     allocate(precdir_recep(nproc))
        if (.not.allocated(precdir_send)) allocate(precdir_send(nproc))
        if (.not.allocated(precdir_recep1))
     $     allocate(precdir_recep1(nproc))
        if (.not.allocated(precdir_send1))
     $     allocate(precdir_send1(nproc))
        if (.not. allocated(bufbegpole)) allocate (bufbegpole(nproc))
        if (.not. allocated(bufbeg)) allocate (bufbeg(nproc))
        if (.not. allocated(bufbegrec)) allocate (bufbegrec(nproc))
        if (.not.allocated(ptorque_recep))
     $      allocate(ptorque_recep(nproc))
        if (.not.allocated(ptorque_send)) allocate(ptorque_send(nproc))
        if (.not.allocated(pbig_recep)) allocate(pbig_recep(nproc))
        if (.not.allocated(pbig_send)) allocate(pbig_send(nproc))
c
c
c     build the two mpi groups for the computation of the multipole interactions
c     (direct and reciprocal space) and associated communicators
c
        allocate (direct_rank(0:nproc-1))
        call MPI_Comm_group(MPI_COMM_WORLD, total_group, ierr)
        do iproc = 0, ndir - 1
          direct_rank(iproc) = iproc
        end do
        call MPI_Group_incl(total_group,ndir,direct_rank,
     $   direct_group,ierr)
        call MPI_Group_excl(total_group,ndir,direct_rank,
     $   rec_group,ierr)
        call MPI_Comm_create(MPI_COMM_WORLD,direct_group,comm_dir,ierr)
        call MPI_Comm_create(MPI_COMM_WORLD,rec_group,comm_rec,ierr)
        if (rank.le.ndir-1) then
          call MPI_COMM_RANK(comm_dir,rank_bis,ierr)
        else
          call MPI_COMM_RANK(comm_rec,rank_bis,ierr)
        end if
c
c       call the dd load balancing routine
c
        call ddpme3d
      else
        if (.not.allocated(glob)) allocate (glob(n))
        if (.not.allocated(loc)) allocate (loc(n))
        if (.not.allocated(repart)) allocate (repart(n))
        if (.not.allocated(domlen)) allocate (domlen(nproc))
        if (.not.allocated(domlenpole)) allocate (domlenpole(nproc))
        if (.not.allocated(domlenpolerec)) 
     $     allocate (domlenpolerec(nproc))
        if (.not.allocated(zbegproc)) allocate (zbegproc(nproc))
        if (.not.allocated(zendproc)) allocate (zendproc(nproc))
        if (.not.allocated(ybegproc)) allocate (ybegproc(nproc))
        if (.not.allocated(yendproc)) allocate (yendproc(nproc))
        if (.not.allocated(xbegproc)) allocate (xbegproc(nproc))
        if (.not.allocated(xendproc)) allocate (xendproc(nproc))
        if (.not.allocated(p_recep1))   allocate(p_recep1(nproc))
        if (.not.allocated(p_send1))  allocate(p_send1(nproc))
        if (.not.allocated(p_recep2))   allocate(p_recep2(nproc))
        if (.not.allocated(p_send2))  allocate(p_send2(nproc))
        if (.not.allocated(pneig_recep))  allocate(pneig_recep(nproc))
        if (.not.allocated(pneig_send))  allocate(pneig_send(nproc))
        if (.not.allocated(pbig_recep))  allocate(pbig_recep(nproc))
        if (.not.allocated(pbig_send))  allocate(pbig_send(nproc))
        if (.not.allocated(precdir_recep))
     $        allocate(precdir_recep(nproc))
        if (.not.allocated(precdir_send)) allocate(precdir_send(nproc))
        if (.not.allocated(precdir_recep1))
     $     allocate(precdir_recep1(nproc))
        if (.not.allocated(precdir_send1))
     $     allocate(precdir_send1(nproc))
        if (.not.allocated(ptorque_recep))
     $      allocate(ptorque_recep(nproc))
        if (.not.allocated(ptorque_send)) allocate(ptorque_send(nproc))
        if (.not. allocated(globrec))
     $       allocate (globrec(n))
        if (.not. allocated(locrec)) allocate (locrec(n))
        if (.not. allocated(globrec1)) allocate (globrec1(n))
        if (.not. allocated(locrec1)) allocate (locrec1(n))
        if (.not. allocated(bufbegrec)) allocate (bufbegrec(nproc))
        if (.not. allocated(bufbegpole)) allocate (bufbegpole(nproc))
        if (.not. allocated(bufbeg)) allocate (bufbeg(nproc))
        call ddpme3d
      end if
      return
      end
c
c     subroutine commposrec: communicate positions for reciprocal space
c
      subroutine commposrec
      use atoms
      use domdec
      use mpi
      implicit none
      integer i,iproc
      integer iglob,iloc
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:)
      real*8, allocatable :: buffers(:,:)
      integer status(MPI_STATUS_SIZE)
      integer tag,ierr
c
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(3,nblocrecdir))
      allocate (buffers(3,nloc))
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
         tag = nproc*rank + precdir_recep(i) + 1
         call MPI_IRECV(buffer(1,bufbeg(precdir_recep(i)+1)),
     $    3*domlen(precdir_recep(i)+1),
     $    MPI_REAL8,precdir_recep(i),tag,MPI_COMM_WORLD,reqrec(i),ierr)
        end if
      end do
c
c     MPI : move in buffer
c
      do i = 1, nloc
        iglob = glob(i)
        buffers(1,i) = x(iglob)
        buffers(2,i) = y(iglob)
        buffers(3,i) = z(iglob)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
         tag = nproc*precdir_send(i) + rank + 1
         call MPI_ISEND(buffers,3*nloc,
     $    MPI_REAL8,precdir_send(i),tag,MPI_COMM_WORLD,
     $    reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
         call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
         call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
         do i = 1, domlen(precdir_recep(iproc)+1)
           iloc = bufbeg(precdir_recep(iproc)+1)+i-1
           iglob = glob(iloc)
           x(iglob) = buffer(1,iloc)
           y(iglob) = buffer(2,iloc)
           z(iglob) = buffer(3,iloc)
         end do
        end if
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commpos : communicate positions to neighboring processes  
c
      subroutine commpos
      use atoms
      use freeze
      use domdec
      use potent
      use mpi
      implicit none
      integer i,tag,ierr,iproc
      integer iloc,iglob
      integer, allocatable :: reqrec(:),reqsend(:)
      real*8, allocatable :: buffer(:,:)
      real*8, allocatable :: buffers(:,:)
      integer status(MPI_STATUS_SIZE)
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(3,max(1,nbloc)))
      allocate (buffers(3,max(1,nloc)))
c
c     communicate positions
c
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_recep
        tag = nproc*rank + pbig_recep(i) + 1
        call MPI_IRECV(buffer(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      do i = 1, nloc
        iglob = glob(i)
        buffers(1,i) = x(iglob)
        buffers(2,i) = y(iglob)
        buffers(3,i) = z(iglob)
      end do
c
c     MPI : begin sending
c
      do i = 1, nbig_send
        tag = nproc*pbig_send(i) + rank + 1
        call MPI_ISEND(buffers,3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,MPI_COMM_WORLD,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, nbig_recep
        do i = 1, domlen(pbig_recep(iproc)+1)
          iloc = bufbeg(pbig_recep(iproc)+1)+i-1
          iglob = glob(iloc)
          x(iglob) = buffer(1,iloc)
          y(iglob) = buffer(2,iloc)
          z(iglob) = buffer(3,iloc)
        end do
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine reassign: deal with particles changing of domain between two time steps
c
      subroutine reassign
      use atoms
      use bound
      use cell
      use domdec
      use freeze
      use moldyn
      use neigh
      use potent
      use mpi
      implicit none
      integer i
      integer iproc,ierr,iglob
      integer iloc
      real*8 xr,yr,zr,eps1,eps2
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: count2(:)
      integer, allocatable :: buf3(:), buflen3(:), bufbeg3(:)
      integer, allocatable :: buflen4(:),bufbeg4(:)
      real*8, allocatable :: buffers(:,:),buffer(:,:)
      real*8, allocatable :: buffers1(:,:),buffer1(:,:)
      integer, allocatable :: buf3bis(:,:)
      integer, allocatable :: glob1(:)
      integer rankloc,commloc,nprocloc,proc,j,jglob
      integer nloc1
      eps1 = 1.0d-10
      eps2 = 1.0d-8
c
      if (use_pmecore) then
        nprocloc = ndir
        commloc  = comm_dir
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = MPI_COMM_WORLD
        rankloc  = rank
      end if
      allocate (glob1(n))
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
c
      allocate (count2(nproc))
      allocate (buflen3(nproc))
      allocate (bufbeg3(nproc))
      allocate (buf3bis(nloc,nproc))
      count2 = 0
      buflen3 = 0
      bufbeg3 = 0
c
c     get particules that changed of domain
c
      do i = 1, nloc
        iglob = glob(i)
        xr = x(iglob)
        yr = y(iglob)
        zr = z(iglob)
        if (use_bounds) call image(xr,yr,zr)
        if (abs(xr-xcell2).lt.eps1) xr = xr-eps2
        if (abs(yr-ycell2).lt.eps1) yr = yr-eps2
        if (abs(zr-zcell2).lt.eps1) zr = zr-eps2
        do iproc = 0, nprocloc-1
          if (iproc.eq.rank) cycle
          if ((zr.ge.zbegproc(iproc+1)).and.
     $     (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
     $    .and.(yr.lt.yendproc(iproc+1)).and.(xr.ge.xbegproc(iproc+1))
     $    .and.(xr.lt.xendproc(iproc+1))) then
            buflen3(iproc+1) = buflen3(iproc+1)+1
            buf3bis(buflen3(iproc+1),iproc+1) = iglob
          end if
        end do
      end do
      bufbeg3(1) = 1
      do iproc = 1, nproc-1
        bufbeg3(iproc+1) = bufbeg3(iproc)+buflen3(iproc)
      end do
c
      allocate (buf3(bufbeg3(nproc)+buflen3(nproc)))
      buf3 = 0
c
      do iproc = 1, nproc
        buf3(bufbeg3(iproc):(bufbeg3(iproc)+buflen3(iproc)-1)) = 
     $    buf3bis(1:buflen3(iproc),iproc)
      end do
c
      if (nprocloc.eq.1) return
c
      if ((use_pmecore).and.(rank.ge.ndir)) goto 20
c
c     get size of buffers
c
      allocate (buflen4(nprocloc))
      buflen4 = 0
      allocate (bufbeg4(nprocloc))
      bufbeg4 = 0
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buflen4(proc+1),1,MPI_INT,
     $   proc,tag,commloc,reqrec(iproc),ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_send(iproc)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buflen3(proc+1),1,MPI_INT,
     $   proc,tag,commloc,reqsend(iproc),ierr)
      end do
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_recep(iproc)
        call MPI_WAIT(reqsend(iproc),status,ierr)
      end do
c
      bufbeg4(pneig_recep(1)+1) = 1
      do iproc = 2, nneig_recep
        bufbeg4(pneig_recep(iproc)+1) = bufbeg4(pneig_recep(iproc-1)+1)
     $    +buflen4(pneig_recep(iproc-1)+1)
      end do
c
c     communicate the corresponding indexes, positions, speed and accelerations
c
      proc = pneig_send(nneig_send)
      allocate (buffers(10,bufbeg3(proc+1)+buflen3(proc+1)))
      proc = pneig_recep(nneig_recep)
      allocate (buffer(10,bufbeg4(proc+1)+buflen4(proc+1)))
c
c     Begin reception 
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buffer(1,bufbeg4(proc+1)),10*buflen4(proc+1),
     $    MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrec(i),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        do j = 0, buflen3(proc+1)-1
          jglob = buf3(bufbeg3(proc+1)+j)
          buffers(1,bufbeg3(proc+1)+j) = x(jglob)
          buffers(2,bufbeg3(proc+1)+j) = y(jglob)
          buffers(3,bufbeg3(proc+1)+j) = z(jglob)
          buffers(4,bufbeg3(proc+1)+j) = v(1,jglob)
          buffers(5,bufbeg3(proc+1)+j) = v(2,jglob)
          buffers(6,bufbeg3(proc+1)+j) = v(3,jglob)
          buffers(7,bufbeg3(proc+1)+j) = a(1,jglob)
          buffers(8,bufbeg3(proc+1)+j) = a(2,jglob)
          buffers(9,bufbeg3(proc+1)+j) = a(3,jglob)
          buffers(10,bufbeg3(proc+1)+j) = jglob
        end do
      end do
c
c     send the positions, the velocities and the accelerations
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buffers(1,bufbeg3(proc+1)),10*buflen3(proc+1),
     $   MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqsend(i),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        call MPI_WAIT(reqrec(i),status,ierr)
        do j = 0, buflen4(proc+1)-1
          iglob = int(buffer(10,bufbeg4(proc+1)+j))
          x(iglob) = buffer(1,bufbeg4(proc+1)+j)
          y(iglob) = buffer(2,bufbeg4(proc+1)+j)
          z(iglob) = buffer(3,bufbeg4(proc+1)+j)
          v(1,iglob) = buffer(4,bufbeg4(proc+1)+j)
          v(2,iglob) = buffer(5,bufbeg4(proc+1)+j)
          v(3,iglob) = buffer(6,bufbeg4(proc+1)+j)
          a(1,iglob) = buffer(7,bufbeg4(proc+1)+j)
          a(2,iglob) = buffer(8,bufbeg4(proc+1)+j)
          a(3,iglob) = buffer(9,bufbeg4(proc+1)+j)
        end do
      end do
c
c     if rattle is used, also send the old coordinates
c
      if (use_rattle) then

        proc = pneig_send(nneig_send)
        allocate (buffers1(4,bufbeg3(proc+1)+buflen3(proc+1)))
        proc = pneig_recep(nneig_recep)
        allocate (buffer1(4,bufbeg4(proc+1)+buflen4(proc+1)))
c
c     Begin reception 
c
        do i = 1, nneig_recep
          proc = pneig_recep(i)
          tag = nprocloc*rankloc + proc + 1
          call MPI_IRECV(buffer1(1,bufbeg4(proc+1)),4*buflen4(proc+1),
     $      MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrec(i),ierr)
        end do
        do i = 1, nneig_send
          proc = pneig_send(i)
          do j = 0, buflen3(proc+1)-1
            jglob = buf3(bufbeg3(proc+1)+j)
            buffers1(1,bufbeg3(proc+1)+j) = jglob
            buffers1(2,bufbeg3(proc+1)+j) = xold(jglob)
            buffers1(3,bufbeg3(proc+1)+j) = yold(jglob)
            buffers1(4,bufbeg3(proc+1)+j) = zold(jglob)
          end do
        end do
c
c       send the old positions
c
        do i = 1, nneig_send
          proc = pneig_send(i)
          tag = nprocloc*proc + rankloc + 1
          call MPI_ISEND(buffers1(1,bufbeg3(proc+1)),4*buflen3(proc+1),
     $     MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqsend(i),ierr)
        end do
        do i = 1, nneig_send
          proc = pneig_send(i)
          call MPI_WAIT(reqsend(i),status,ierr)
        end do
        do i = 1, nneig_recep
          proc = pneig_recep(i)
          call MPI_WAIT(reqrec(i),status,ierr)
          do j = 0, buflen4(proc+1)-1
            iglob = int(buffer1(1,bufbeg4(proc+1)+j))
            xold(iglob) = buffer1(2,bufbeg4(proc+1)+j)
            yold(iglob) = buffer1(3,bufbeg4(proc+1)+j)
            zold(iglob) = buffer1(4,bufbeg4(proc+1)+j)
          end do
        end do
        deallocate (buffer1)
        deallocate (buffers1)
      end if

c
c     reorder indexes accordingly and build local repart array
c
      repart = 0
c
c     remove atoms that left the domain 
c
      nloc1 = 0
      glob1 = 0
      do i = 1, nloc
        iglob = glob(i)
        do j = 1, bufbeg3(nproc)+buflen3(nproc)-1
          jglob = buf3(j)
          if (iglob.eq.jglob) goto 10
        end do
        nloc1 = nloc1 + 1
        glob1(nloc1) = iglob
        loc(iglob) = nloc1
        repart(iglob) = rank
 10     continue
      end do
c
c     add atoms that entered the domain
c
      proc = pneig_recep(nneig_recep)+1
      do j = 1, bufbeg4(proc)+buflen4(proc)-1
        jglob = int(buffer(10,j))
        nloc1 = nloc1 + 1
        glob1(nloc1) = jglob
        loc(jglob) = nloc1
        repart(jglob) = rank
      end do
      deallocate (buffer)
      deallocate (buffers)
c
      nloc = nloc1
      domlen(rank+1) = nloc
      glob = glob1
c
 20   call orderbuffer(.false.)
c
      deallocate (buf3bis)
      deallocate (glob1)
      deallocate (count2)
      deallocate (buf3)
      deallocate (buflen3)
      deallocate (bufbeg3)
      deallocate (reqsend)
      deallocate (reqrec)
      return
      end
c
c     subroutine allocstep: deallocate arrays and reallocate them with proper size
c     (memory distribution)
c
      subroutine allocstep
      use deriv
      use domdec
      implicit none
c
      if (allocated(desum)) deallocate (desum)
      allocate (desum(3,nbloc))
      desum = 0d0
      if (allocated(decrec)) deallocate (decrec)
      allocate (decrec(3,nblocrec))
      decrec = 0d0
      if (allocated(demrec)) deallocate (demrec)
      allocate (demrec(3,nblocrec))
      demrec = 0d0
      if (allocated(deprec)) deallocate (deprec)
      allocate (deprec(3,nblocrec))
      deprec = 0d0
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
      if (allocated(debt)) deallocate (debt)
      allocate (debt(3,nbloc))
      if (allocated(dett)) deallocate (dett)
      allocate (dett(3,nbloc))
      if (allocated(dev)) deallocate (dev)
      allocate (dev(3,nbloc))
      if (allocated(dec)) deallocate (dec)
      allocate (dec(3,nbloc))
      if (allocated(dem)) deallocate (dem)
      allocate (dem(3,nbloc))
      if (allocated(dep)) deallocate (dep)
      allocate (dep(3,nblocloop))
!!     if (allocated(dep1)) deallocate (dep1)
!!     allocate (dep1(nblocloop))
!!     if (allocated(dep2)) deallocate (dep2)
!!     allocate (dep2(nblocloop))
!!     if (allocated(dep3)) deallocate (dep3)
!!     allocate (dep3(nblocloop))
      if (allocated(deg)) deallocate (deg)
      allocate (deg(3,nbloc))
      if (allocated(dex)) deallocate (dex)
      allocate (dex(3,nbloc))
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
      implicit none
      logical fast
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
      if (allocated(debt)) deallocate (debt)
      allocate (debt(3,nbloc))
      if (allocated(dett)) deallocate (dett)
      allocate (dett(3,nbloc))
      if (allocated(decrec)) deallocate (decrec)
      allocate (decrec(3,nblocrec))
      if (allocated(demrec)) deallocate (demrec)
      allocate (demrec(3,nblocrec))
      if (allocated(deprec)) deallocate (deprec)
      allocate (deprec(3,nblocrec))
      if (allocated(debond)) deallocate (debond)
      allocate (debond(3,nbloc))
c
      if (.not.(fast)) then
        decrec = 0d0
        demrec = 0d0
        deprec = 0d0
      end if
c
      if (allocated(dev)) deallocate (dev)
      allocate (dev(3,nbloc))
      if (allocated(dec)) deallocate (dec)
      allocate (dec(3,nbloc))
      if (allocated(dem)) deallocate (dem)
      allocate (dem(3,nbloc))
      if (allocated(dep)) deallocate (dep)
      allocate (dep(3,nblocloop))
!!     if (allocated(dep1)) deallocate (dep1)
!!     allocate (dep1(nblocloop))
!!     if (allocated(dep2)) deallocate (dep2)
!!     allocate (dep2(nblocloop))
!!     if (allocated(dep3)) deallocate (dep3)
!!     allocate (dep3(nblocloop))

      if (allocated(deg)) deallocate (deg)
      allocate (deg(3,nbloc))
      if (allocated(dex)) deallocate (dex)
      allocate (dex(3,nbloc))
c
      return
      end
c
c     subroutine sendall : everybody gets all the positions, speed and accel of the system
c
      subroutine sendall
      use atoms
      use domdec
      use moldyn
      use mpi
      implicit none
      integer i,iglob,iloc,iproc,tag,ierr
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: postemp(:,:)
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (postemp(9,n))
      postemp = 0d0
c
c     put the arrays to be sent in the right order
c
      do i = 1, nloc
        iglob = glob(i)
        postemp(1,i) = x(iglob)
        postemp(2,i) = y(iglob)
        postemp(3,i) = z(iglob)
        postemp(4,i) = v(1,iglob)
        postemp(5,i) = v(2,iglob)
        postemp(6,i) = v(3,iglob)
        postemp(7,i) = a(1,iglob)
        postemp(8,i) = a(2,iglob)
        postemp(9,i) = a(3,iglob)
      end do
c
c     get the sizes of the domains
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(domlen(iproc+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        call MPI_WAIT(reqrec(iproc+1),status,ierr)
        call MPI_WAIT(reqsend(iproc+1),status,ierr)
        end if
      end do
c
      bufbeg(rank+1) = 1
      count = domlen(rank+1)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          bufbeg(iproc+1) = count + 1
          count = count + domlen(iproc+1)
        end if
      end do
c
c     get the indexes
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(glob(bufbeg(iproc+1)),domlen(iproc+1),MPI_INT,
     $   iproc,tag,MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     get the positions
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(postemp(1,bufbeg(iproc+1)),9*domlen(iproc+1),
     $   MPI_REAL8,iproc,tag,MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,9*nloc,MPI_REAL8,iproc,tag,
     $    MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     put it back in the global arrays
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          do i = 1, domlen(iproc+1)
            iloc = bufbeg(iproc+1)+i-1
            iglob = glob(iloc)
            x(iglob) = postemp(1,iloc)
            y(iglob) = postemp(2,iloc)
            z(iglob) = postemp(3,iloc)
            v(1,iglob) = postemp(4,iloc)
            v(2,iglob) = postemp(5,iloc)
            v(3,iglob) = postemp(6,iloc)
            a(1,iglob) = postemp(7,iloc)
            a(2,iglob) = postemp(8,iloc)
            a(3,iglob) = postemp(9,iloc)
          end do
        end if
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (postemp)
      return
      end
c
c     subroutine ddnumber : get the number of subdivision along each axis for
c     3d spatial decomposition
c
      subroutine ddnumber(num,nx,ny,nz,istep)
      use boxes
      use domdec
      use inform
      use iounit
      implicit none
      integer num,istep
      integer key(3),list2(3)
      integer, allocatable :: d(:)
      real*8 list1(3)
      integer nx,ny,nz,n1,n2,n3,i,res
 10   format('Nx = ',I5,2x,'Ny = ',I5,2x,'Nz = ',I5,2x)
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
          nx = num
          ny = 1
          nz = 1
        else if (key(3).eq.2) then
          nx = 1
          ny = num
          nz = 1
        else
          nx = 1
          ny = 1
          nz = num
        end if
      else if (i.eq.2) then
        if (key(3).eq.1) then
          if (key(2).eq.2) then
            nx = d(2)
            ny = d(1)
            nz = 1
          else
            nx = d(2)
            ny = 1
            nz = d(1)
          end if
        else if (key(3).eq.2) then
          if (key(2).eq.1) then
            nx = d(1)
            ny = d(2)
            nz = 1
          else
            nx = 1
            ny = d(2)
            nz = d(1)
          end if
        else
          if (key(3).eq.1) then
            nx = d(1)
            ny = 1
            nz = d(2)
          else
            nx = 1
            ny = d(1)
            nz = d(2)
          end if
        end if
      else
c
        n1 = floor(num**(1.0/3.0))
        do i = 0, n1-2
          res = mod(num,n1-i)
          if (res.eq.0) goto 20
        end do
 20     continue
        n1 = n1-i
        n2 = floor((num/n1)**(1.0/2.0))
        do i = 0, n2-2
          res = mod(num/n1,n2-i)
          if (res.eq.0) goto 30
        end do
 30     continue
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
          if (n1.eq.0) goto 60
          do i = 0, n1-2
            res = mod(num,n1-i)
            if (res.eq.0) goto 40
          end do
 40       continue
          n1 = n1-i
          n2 = floor((num/n1)**(1.0/2.0))
          do i = 0, n2-2
            res = mod(num/n1,n2-i)
            if (res.eq.0) goto 50
          end do
 50       continue
          n2 = n2 - i
          n3 = num/(n1*n2)
          list2(1) = n1
          list2(2) = n2
          list2(3) = n3
          call sort(3,list2)
        end if
c
 60     if (key(3).eq.1) then
          if (key(2).eq.2) then
            nx = list2(3)
            ny = list2(2)
            nz = list2(1)
          else
            nx = list2(3)
            ny = list2(1)
            nz = list2(2)
          end if
        else if (key(3).eq.2) then
          if (key(2).eq.1) then
            nx = list2(2)
            ny = list2(3)
            nz = list2(1)
          else
            nx = list2(1)
            ny = list2(3)
            nz = list2(2)
          end if
        else
          if (key(2).eq.1) then
            nx = list2(2)
            ny = list2(1)
            nz = list2(3)
          else
            nx = list2(1)
            ny = list2(2)
            nz = list2(3)
          end if
        end if
      end if
      if (istep.eq.0.and.verbose) then
        if (rank.eq.0) write(iout,*) '3D Domain Decomposition'
        if (rank.eq.0) write(iout,10) nx,ny,nz
      end if
      deallocate (d)
      return
      end
c
c     subroutine reassignrespa: deal with particles changing of domain between two time steps
c
      subroutine reassignrespa(fast,ialt,nalt)
      use atoms
      use bound
      use cell
      use domdec
      use moldyn
      use neigh
      use potent
      use mpi
      implicit none
      integer i,iproc,ierr,iglob,iloc
      real*8 xr,yr,zr,eps1,eps2
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: count2(:)
      integer, allocatable :: buf3(:), buflen3(:), bufbeg3(:)
      integer, allocatable :: buflen4(:),bufbeg4(:)
      real*8, allocatable :: buffers(:,:),buffer(:,:)
      integer, allocatable :: buf3bis(:,:)
      integer, allocatable :: glob1(:)
      integer rankloc,commloc,nprocloc,proc,j,jglob
      integer nloc1
      integer ialt,nalt
      logical fast
      eps1 = 1.0d-10
      eps2 = 1.0d-8
c
      if (ialt.ne.nalt) return
c
      if (use_pmecore) then
        nprocloc = ndir
        commloc  = comm_dir
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = MPI_COMM_WORLD
        rankloc  = rank
      end if
      allocate (glob1(n))
c
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
      allocate (count2(nproc))
      allocate (buflen3(nproc))
      allocate (bufbeg3(nproc))
      allocate (buf3bis(nloc,nproc))
      count2 = 0
      buflen3 = 0
      bufbeg3 = 0
c
c     get particules that changed of domain
c
      do i = 1, nloc
        iglob = glob(i)
        xr = x(iglob)
        yr = y(iglob)
        zr = z(iglob)
        if (use_bounds) call image(xr,yr,zr)
        if (abs(xr-xcell2).lt.eps1) xr = xr-eps2
        if (abs(yr-ycell2).lt.eps1) yr = yr-eps2
        if (abs(zr-zcell2).lt.eps1) zr = zr-eps2
        do iproc = 0, nprocloc-1
          if (iproc.eq.rank) cycle
          if ((zr.ge.zbegproc(iproc+1)).and.
     $     (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
     $    .and.(yr.lt.yendproc(iproc+1)).and.(xr.ge.xbegproc(iproc+1))
     $    .and.(xr.lt.xendproc(iproc+1))) then
            buflen3(iproc+1) = buflen3(iproc+1)+1
            buf3bis(buflen3(iproc+1),iproc+1) = iglob
          end if
        end do
      end do
      bufbeg3(1) = 1
      do iproc = 1, nproc-1
        bufbeg3(iproc+1) = bufbeg3(iproc)+buflen3(iproc)
      end do
c
      allocate (buf3(bufbeg3(nproc)+buflen3(nproc)))
      buf3 = 0
c
      do iproc = 1, nproc
        buf3(bufbeg3(iproc):(bufbeg3(iproc)+buflen3(iproc)-1)) = 
     $    buf3bis(1:buflen3(iproc),iproc)
      end do
c
      if (nprocloc.eq.1) return
c
      if ((use_pmecore).and.(rank.ge.ndir)) goto 20
c
c     get size of buffers
c
      allocate (buflen4(nprocloc))
      buflen4 = 0
      allocate (bufbeg4(nprocloc))
      bufbeg4 = 0
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buflen4(proc+1),1,MPI_INT,
     $   proc,tag,commloc,reqrec(iproc),ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_send(iproc)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buflen3(proc+1),1,MPI_INT,
     $   proc,tag,commloc,reqsend(iproc),ierr)
      end do
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nneig_send
        tag = nprocloc*proc + rankloc + 1
        call MPI_WAIT(reqsend(iproc),status,ierr)
      end do
c
      bufbeg4(pneig_recep(1)+1) = 1
      do iproc = 2, nneig_recep
        bufbeg4(pneig_recep(iproc)+1) = bufbeg4(pneig_recep(iproc-1)+1)
     $    +buflen4(pneig_recep(iproc-1)+1)
      end do
c
c     communicate the corresponding indexes, positions, speed and accelerations
c
      proc = pneig_send(nneig_send)
      allocate (buffers(10,bufbeg3(proc+1)+buflen3(proc+1)))
      proc = pneig_recep(nneig_recep)
      allocate (buffer(10,bufbeg4(proc+1)+buflen4(proc+1)))
c
c     Begin reception 
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buffer(1,bufbeg4(proc+1)),10*buflen4(proc+1),
     $    MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrec(i),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        do j = 0, buflen3(proc+1)-1
          jglob = buf3(bufbeg3(proc+1)+j)
          buffers(1,bufbeg3(proc+1)+j) = x(jglob)
          buffers(2,bufbeg3(proc+1)+j) = y(jglob)
          buffers(3,bufbeg3(proc+1)+j) = z(jglob)
          buffers(4,bufbeg3(proc+1)+j) = v(1,jglob)
          buffers(5,bufbeg3(proc+1)+j) = v(2,jglob)
          buffers(6,bufbeg3(proc+1)+j) = v(3,jglob)
          buffers(7,bufbeg3(proc+1)+j) = a(1,jglob)
          buffers(8,bufbeg3(proc+1)+j) = a(2,jglob)
          buffers(9,bufbeg3(proc+1)+j) = a(3,jglob)
          buffers(10,bufbeg3(proc+1)+j) = jglob
        end do
      end do
c
c     send the positions, the velocities and the accelerations
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buffers(1,bufbeg3(proc+1)),10*buflen3(proc+1),
     $   MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqsend(i),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        tag = nprocloc*proc + rankloc + 1
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_WAIT(reqrec(i),status,ierr)
        do j = 0, buflen4(proc+1)-1
          iglob = int(buffer(10,bufbeg4(proc+1)+j))
          x(iglob) = buffer(1,bufbeg4(proc+1)+j)
          y(iglob) = buffer(2,bufbeg4(proc+1)+j)
          z(iglob) = buffer(3,bufbeg4(proc+1)+j)
          v(1,iglob) = buffer(4,bufbeg4(proc+1)+j)
          v(2,iglob) = buffer(5,bufbeg4(proc+1)+j)
          v(3,iglob) = buffer(6,bufbeg4(proc+1)+j)
          a(1,iglob) = buffer(7,bufbeg4(proc+1)+j)
          a(2,iglob) = buffer(8,bufbeg4(proc+1)+j)
          a(3,iglob) = buffer(9,bufbeg4(proc+1)+j)
        end do
      end do
c
c     reorder indexes accordingly and build local repart array
c
      repart = 0
c
c     remove atoms that left the domain 
c
      nloc1 = 0
      glob1 = 0
      do i = 1, nloc
        iglob = glob(i)
        do j = 1, bufbeg3(nproc)+buflen3(nproc)-1
          jglob = buf3(j)
          if (iglob.eq.jglob) goto 10
        end do
        nloc1 = nloc1 + 1
        glob1(nloc1) = iglob
        loc(iglob) = nloc1
        repart(iglob) = rank
 10     continue
      end do
c
c     add atoms that entered the domain
c
      proc = pneig_recep(nneig_recep)+1
      do j = 1, bufbeg4(proc)+buflen4(proc)-1
        jglob = int(buffer(10,j))
        nloc1 = nloc1 + 1
        glob1(nloc1) = jglob
        loc(jglob) = nloc1
        repart(jglob) = rank
      end do
c
      nloc = nloc1
      domlen(rank+1) = nloc
      glob = glob1
c
 20   call orderbuffer(.false.)
c
      deallocate (buf3bis)
      deallocate (glob1)
      deallocate (count2)
      deallocate (buf3)
      deallocate (buflen3)
      deallocate (bufbeg3)
      deallocate (reqsend)
      deallocate (reqrec)
      return
      end
c
c     subroutine commposrespa : deal with communications of positions between two time steps
c     with respa integrator
c
      subroutine commposrespa(fast)
      use atoms
      use domdec
      use potent
      use mpi
      implicit none
      integer i,tag,ierr,iproc
      integer iloc,iglob
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:)
      real*8, allocatable :: buffers(:,:)
      integer status(MPI_STATUS_SIZE)
      integer, allocatable :: ptemp_recep(:),ptemp_send(:)
      integer ntemp_recep,ntemp_send
      logical fast
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
      allocate (ptemp_recep(nproc))
      allocate (ptemp_send(nproc))
c

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
c
c     allocate some arrays
c
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (buffer(3,nbloc))
      allocate (buffers(3,nloc))
c
c     communicate positions
c
c
c     MPI : begin reception in buffer
c
      do i = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(i) + 1
        call MPI_IRECV(buffer(1,bufbeg(ptemp_recep(i)+1)),
     $   3*domlen(ptemp_recep(i)+1),
     $   MPI_REAL8,ptemp_recep(i),tag,MPI_COMM_WORLD,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      do i = 1, nloc
        iglob = glob(i)
        buffers(1,i) = x(iglob)
        buffers(2,i) = y(iglob)
        buffers(3,i) = z(iglob)
      end do
c
c     MPI : begin sending
c
      do i = 1, ntemp_send
        tag = nproc*ptemp_send(i) + rank + 1
        call MPI_ISEND(buffers,3*nloc,
     $   MPI_REAL8,ptemp_send(i),tag,MPI_COMM_WORLD,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, ntemp_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, ntemp_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, ntemp_recep
        do i = 1, domlen(ptemp_recep(iproc)+1)
          iloc = bufbeg(ptemp_recep(iproc)+1)+i-1
          iglob = glob(iloc)
          x(iglob) = buffer(1,iloc)
          y(iglob) = buffer(2,iloc)
          z(iglob) = buffer(3,iloc)
        end do
      end do
c
      deallocate (ptemp_recep)
      deallocate (ptemp_send)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine sendallpos : everybody gets all the positions
c
      subroutine sendallpos
      use atoms
      use domdec
      use mpi
      implicit none
      integer i,iglob,iloc,iproc,tag,ierr
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: postemp(:,:)
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (postemp(3,n))
      postemp = 0d0
c
c     put the arrays to be sent in the right order
c
      do i = 1, nloc
        iglob = glob(i)
        postemp(1,i) = x(iglob)
        postemp(2,i) = y(iglob)
        postemp(3,i) = z(iglob)
      end do
c
c     get the sizes of the domains
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(domlen(iproc+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        call MPI_WAIT(reqsend(iproc+1),status,ierr)
        call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
      bufbeg(rank+1) = 1
      count = domlen(rank+1)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          if (domlen(iproc+1).eq.0) then
            bufbeg(iproc+1) = 1
          else
            bufbeg(iproc+1) = count + 1
          end if
          count = count + domlen(iproc+1)
        end if
      end do
c
c     get the indexes
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(glob(bufbeg(iproc+1)),domlen(iproc+1),MPI_INT,
     $   iproc,tag,MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     get the positions
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(postemp(1,bufbeg(iproc+1)),3*domlen(iproc+1),
     $   MPI_REAL8,iproc,tag,MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,3*nloc,MPI_REAL8,iproc,tag,
     $    MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     put it back in the global arrays
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          do i = 1, domlen(iproc+1)
            iloc = bufbeg(iproc+1)+i-1
            iglob = glob(iloc)
            x(iglob) = postemp(1,iloc)
            y(iglob) = postemp(2,iloc)
            z(iglob) = postemp(3,iloc)
          end do
        end if
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (postemp)
      return
      end
c
c     subroutine sendvecmin : deal with communications during minimization routine
c
      subroutine sendvecmin(xx)
      use atoms
      use domdec
      use mpi
      implicit none
      integer i,iglob,iloc,iproc,tag,ierr,j
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: temp(:,:)
      real*8 xx(3*n)
c
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (temp(3,n))
      temp = 0d0
c
c     put the arrays to be sent in the right order
c
      do i = 1, nloc
        iglob = glob(i)
        do j = 1, 3
          temp(j,i) = xx(3*(iglob-1)+j)
        end do
      end do
c
c     get the sizes of the domains
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(domlen(iproc+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        call MPI_WAIT(reqsend(iproc+1),status,ierr)
        call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
      bufbeg(rank+1) = 1
      count = domlen(rank+1)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          if (domlen(iproc+1).eq.0) then
            bufbeg(iproc+1) = 1
          else
            bufbeg(iproc+1) = count + 1
          end if
          count = count + domlen(iproc+1)
        end if
      end do
c
c     get the indexes
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(glob(bufbeg(iproc+1)),domlen(iproc+1),MPI_INT,
     $   iproc,tag,MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     get the values of xx
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(temp(1,bufbeg(iproc+1)),3*domlen(iproc+1),
     $   MPI_REAL8,iproc,tag,MPI_COMM_WORLD,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(temp,3*nloc,MPI_REAL8,iproc,tag,
     $    MPI_COMM_WORLD,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     put it back in the global arrays
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          do i = 1, domlen(iproc+1)
            iloc = bufbeg(iproc+1)+i-1
            iglob = glob(iloc)
            xx(3*(iglob-1)+1) = temp(1,iloc)
            xx(3*(iglob-1)+2) = temp(2,iloc)
            xx(3*(iglob-1)+3) = temp(3,iloc)
          end do
        end if
      end do
c
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (temp)
      return
      end
