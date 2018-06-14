
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
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'openmp.i'
      include 'potent.i'
      include 'mpif.h'
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
        if (.not. associated(glob)) allocate (glob(n))
        if (.not. associated(loc)) allocate (loc(n))
        if (.not. associated(globrec)) allocate (globrec(n))
        if (.not. associated(locrec)) allocate (locrec(n))
        if (.not. associated(globrec1)) allocate (globrec1(n))
        if (.not. associated(locrec1)) allocate (locrec1(n))
        if (.not.associated(repart)) allocate (repart(n))
        if (.not.associated(domlen)) allocate (domlen(nproc))
        if (.not.associated(zbegproc)) allocate (zbegproc(nproc))
        if (.not.associated(zendproc)) allocate (zendproc(nproc))
        if (.not.associated(ybegproc)) allocate (ybegproc(nproc))
        if (.not.associated(yendproc)) allocate (yendproc(nproc))
        if (.not.associated(xbegproc)) allocate (xbegproc(nproc))
        if (.not.associated(xendproc)) allocate (xendproc(nproc))
        if (.not.associated(p_recep1))   allocate(p_recep1(nproc))
        if (.not.associated(p_send1))  allocate(p_send1(nproc))
        if (.not.associated(p_recep2))   allocate(p_recep2(nproc))
        if (.not.associated(p_send2))  allocate(p_send2(nproc))
        if (.not.associated(pneig_recep))   allocate(pneig_recep(nproc))
        if (.not.associated(pneig_send))  allocate(pneig_send(nproc))
        if (.not.associated(precdir_recep))
     $     allocate(precdir_recep(nproc))
        if (.not.associated(precdir_send)) allocate(precdir_send(nproc))
        if (.not.associated(precdir_recep1))
     $     allocate(precdir_recep1(nproc))
        if (.not.associated(precdir_send1))
     $     allocate(precdir_send1(nproc))
        if (.not. associated(bufbegpole)) allocate (bufbegpole(nproc))
        if (.not. associated(bufbeg)) allocate (bufbeg(nproc))
        if (.not. associated(bufbegrec)) allocate (bufbegrec(nproc))
        if (.not.associated(ptorque_recep))
     $      allocate(ptorque_recep(nproc))
        if (.not.associated(ptorque_send)) allocate(ptorque_send(nproc))
        if (.not.associated(pbig_recep)) allocate(pbig_recep(nproc))
        if (.not.associated(pbig_send)) allocate(pbig_send(nproc))
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
        if (.not.associated(glob)) allocate (glob(n))
        if (.not.associated(loc)) allocate (loc(n))
        if (.not.associated(repart)) allocate (repart(n))
        if (.not.associated(domlen)) allocate (domlen(nproc))
        if (.not.associated(zbegproc)) allocate (zbegproc(nproc))
        if (.not.associated(zendproc)) allocate (zendproc(nproc))
        if (.not.associated(ybegproc)) allocate (ybegproc(nproc))
        if (.not.associated(yendproc)) allocate (yendproc(nproc))
        if (.not.associated(xbegproc)) allocate (xbegproc(nproc))
        if (.not.associated(xendproc)) allocate (xendproc(nproc))
        if (.not.associated(p_recep1))   allocate(p_recep1(nproc))
        if (.not.associated(p_send1))  allocate(p_send1(nproc))
        if (.not.associated(p_recep2))   allocate(p_recep2(nproc))
        if (.not.associated(p_send2))  allocate(p_send2(nproc))
        if (.not.associated(pneig_recep))  allocate(pneig_recep(nproc))
        if (.not.associated(pneig_send))  allocate(pneig_send(nproc))
        if (.not.associated(pbig_recep))  allocate(pbig_recep(nproc))
        if (.not.associated(pbig_send))  allocate(pbig_send(nproc))
        if (.not.associated(precdir_recep))
     $        allocate(precdir_recep(nproc))
        if (.not.associated(precdir_send)) allocate(precdir_send(nproc))
        if (.not.associated(precdir_recep1))
     $     allocate(precdir_recep1(nproc))
        if (.not.associated(precdir_send1))
     $     allocate(precdir_send1(nproc))
        if (.not.associated(ptorque_recep))
     $      allocate(ptorque_recep(nproc))
        if (.not.associated(ptorque_send)) allocate(ptorque_send(nproc))
        if (.not. associated(globrec))
     $       allocate (globrec(n))
        if (.not. associated(locrec)) allocate (locrec(n))
        if (.not. associated(globrec1)) allocate (globrec1(n))
        if (.not. associated(locrec1)) allocate (locrec1(n))
        if (.not. associated(bufbegrec)) allocate (bufbegrec(nproc))
        if (.not. associated(bufbegpole)) allocate (bufbegpole(nproc))
        if (.not. associated(bufbeg)) allocate (bufbeg(nproc))
        call ddpme3d
      end if
      return
      end
c
c     subroutine commposrec: communicate positions for reciprocal space
c
      subroutine commposrec
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,iproc
      integer iglob,iloc
      integer, allocatable :: req(:)
      real*8, allocatable :: buffer(:,:)
      real*8, allocatable :: buffers(:,:)
      integer status(MPI_STATUS_SIZE)
      integer tag,ierr
c
c
c     allocate some arrays
c
      allocate (req(nproc*nproc))
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
     $    MPI_REAL8,precdir_recep(i),tag,MPI_COMM_WORLD,req(tag),ierr)
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
     $    req(tag),ierr)
        end if
      end do
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
         tag = nproc*rank + precdir_recep(i) + 1
         call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
         tag = nproc*precdir_send(i) + rank + 1
         call MPI_WAIT(req(tag),status,ierr)
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
      deallocate (req)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commpos : communicate positions to neighboring processes  
c
      subroutine commpos
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'potent.i'
      include 'mpif.h'
      integer i,tag,ierr,iproc
      integer iloc,iglob
      integer, allocatable :: req(:)
      real*8, allocatable :: buffer(:,:)
      real*8, allocatable :: buffers(:,:)
      integer status(MPI_STATUS_SIZE)
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     allocate some arrays
c
      allocate (req(nproc*nproc))
      allocate (buffer(3,nbloc))
      allocate (buffers(3,nloc))
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
     $   MPI_REAL8,pbig_recep(i),tag,MPI_COMM_WORLD,req(tag),ierr)
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
     $   req(tag),ierr)
      end do
c
      do i = 1, nbig_recep
        tag = nproc*rank + pbig_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nbig_send
        tag = nproc*pbig_send(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
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
      deallocate (req)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine reassign: deal with particles changing of domain between two time steps
c
      subroutine reassign
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'openmp.i'
      include 'potent.i'
      include 'mpif.h'
      include 'neigh.i'
      integer i
      integer iproc,ierr,iglob
      integer iloc
      real*8 xr,yr,zr,eps1,eps2
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reparttemp(:)
      integer, allocatable :: req(:)
      integer, allocatable :: count2(:)
      integer, allocatable :: buf3(:), buflen3(:), bufbeg3(:)
      integer, allocatable :: domlentemp(:)
      integer nprocloc
      eps1 = 1.0d-10
      eps2 = 1.0d-8
c
      if (use_pmecore) then
        nprocloc = ndir
      else
        nprocloc = nproc
      end if
      allocate (domlentemp(nproc))
c
      allocate (req(nproc*nproc))
      allocate (reparttemp(nblocrecdir))
      reparttemp = 0
      domlentemp = domlen
c
      allocate (count2(nproc))
      allocate (buflen3(nproc))
      allocate (bufbeg3(nproc))
      count2 = 0
      buflen3 = 0
      bufbeg3 = 0
c
c     recompute number of atoms per domain
c
      nlocold = nloc
      nblocold = nbloc
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
          if ((zr.ge.zbegproc(iproc+1)).and.
     $     (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
     $    .and.(yr.lt.yendproc(iproc+1)).and.(xr.ge.xbegproc(iproc+1))
     $    .and.(xr.lt.xendproc(iproc+1))) then
             reparttemp(i) = iproc
          end if
        end do
      end do
c
c     every process sends his portion of reparttemp
c
      do i = 1, nneig_recep
        tag = nproc*rank + pneig_recep(i) + 1
        call MPI_IRECV(reparttemp(bufbeg(pneig_recep(i)+1)),
     $   domlen(pneig_recep(i)+1),MPI_INT,pneig_recep(i),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     begin sending
c
      do i = 1, nneig_send
        tag = nproc*pneig_send(i) + rank + 1
        call MPI_ISEND(reparttemp,
     $  domlen(rank+1),MPI_INT,pneig_send(i),tag,MPI_COMM_WORLD,
     $  req(tag),ierr)
      end do
c
      do i = 1, nneig_recep
        tag = nproc*rank + pneig_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nneig_send
        tag = nproc*pneig_send(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
        tag = nproc*rank + precdir_recep(i) + 1
        call MPI_IRECV(reparttemp(bufbeg(precdir_recep(i)+1)),
     $   domlen(precdir_recep(i)+1),MPI_INT,precdir_recep(i),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
c
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
        tag = nproc*precdir_send(i) + rank + 1
        call MPI_ISEND(reparttemp,
     $  domlen(rank+1),MPI_INT,precdir_send(i),tag,MPI_COMM_WORLD,
     $  req(tag),ierr)
        end if
      end do
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
        tag = nproc*rank + precdir_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
        tag = nproc*precdir_send(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
      domlen = 0
c
      do i = 1, nloc
        iglob = glob(i)
        repart(iglob) = reparttemp(i)
        domlen(repart(iglob)+1) = domlen(repart(iglob)+1)+1
      end do
c
      do iproc = 1, nneig_recep
        do i = 1, domlentemp(pneig_recep(iproc)+1)
        iloc = bufbeg(pneig_recep(iproc)+1)+i-1
        iglob = glob(iloc)
          repart(iglob) = reparttemp(iloc)
          domlen(repart(iglob)+1) = domlen(repart(iglob)+1)+1
        end do
      end do
      do iproc = 1, nrecdir_recep
        do i = 1, domlentemp(precdir_recep(iproc)+1)
          iloc = bufbeg(precdir_recep(iproc)+1)+i-1
          iglob = glob(iloc)
          repart(iglob) = reparttemp(iloc)
        end do
      end do
c
c    get the positions to send (change of domains)
c
      do i = 1, nloc
        iglob = glob(i)
        if (rank.ne.repart(iglob)) then
          buflen3(repart(iglob)+1) = buflen3(repart(iglob)+1)+1
        end if
      end do
      bufbeg3(1) = 1
      do iproc = 1, nproc-1
        bufbeg3(iproc+1) = bufbeg3(iproc)+buflen3(iproc)
      end do
c
      allocate (buf3(bufbeg3(nproc)+buflen3(nproc)))
      buf3 = 0
c
      do i = 1, nloc
        iglob = glob(i)
        if (rank.ne.repart(iglob)) then
          buf3(bufbeg3(repart(iglob)+1)+count2(repart(iglob)+1))=
     $      iglob
          count2(repart(iglob)+1) = count2(repart(iglob)+1) + 1
        end if
      end do
c
      call orderbuffer(.false.,domlentemp)
      nloc = domlen(rank+1)
c
c
c     send the positions, the speed and the accelerations
c
      if ((.not.(use_pmecore)).or.(rank.le.ndir-1)) then
        call commreassign(buflen3,bufbeg3,buf3)
      end if
c
      deallocate (domlentemp)
      deallocate (count2)
      deallocate (buf3)
      deallocate (buflen3)
      deallocate (bufbeg3)
      deallocate (reparttemp)
      deallocate (req)
      return
      end
c
c     subroutine commreassign: send and receive positions of atoms that changed of domain
c
      subroutine commreassign(buflen3,bufbeg3,buf3)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'moldyn.i'
      include 'potent.i'
      include 'mpif.h'
      integer i,j,iproc,tag,ierr,proc
      integer status(MPI_STATUS_SIZE)
      integer buflen3(*),bufbeg3(*),buf3(*)
      integer, allocatable :: buf4(:),buflen4(:),bufbeg4(:)
      integer, allocatable :: req(:),req2(:)
      real*8, allocatable :: buffers(:,:),buffer(:,:)
      integer nprocloc,commloc,rankloc
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
c
      if (nprocloc.eq.1) return
c
      allocate (req(nprocloc*nprocloc))
      allocate (req2(nprocloc*nprocloc))
c
c     deal with Direct-Recip communications
c
      allocate (buflen4(nprocloc))
      buflen4 = 0
      allocate (bufbeg4(nprocloc))
      bufbeg4 = 0
c
c     send and receive sizes of the buffers
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buflen4(proc+1),1,MPI_INT,
     $   proc,tag,commloc,req(tag),ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_send(iproc)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buflen3(proc+1),1,MPI_INT,
     $   proc,tag,commloc,req(tag),ierr)
      end do
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        tag = nprocloc*rankloc + proc + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_recep(iproc)
        tag = nprocloc*proc + rankloc + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
      bufbeg4(pneig_recep(1)+1) = 1
      do iproc = 2, nneig_recep
        bufbeg4(pneig_recep(iproc)+1) = bufbeg4(pneig_recep(iproc-1)+1)
     $    +buflen4(pneig_recep(iproc-1)+1)
      end do
c
      proc = pneig_recep(nneig_recep)
      allocate (buf4(bufbeg4(proc+1)+buflen4(proc+1)))
      buf4 = 0
c
c     send and receive list of corresponding indexes
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buf4(bufbeg4(proc+1)),buflen4(proc+1),
     $   MPI_INT,proc,tag,MPI_COMM_WORLD,req2(tag),ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_send(iproc)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buf3(bufbeg3(proc+1)),buflen3(proc+1),
     $   MPI_INT,proc,tag,MPI_COMM_WORLD,req2(tag),ierr)
      end do
c
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        tag = nprocloc*rankloc + proc + 1
        call MPI_WAIT(req2(tag),status,ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_send(iproc)
        tag = nprocloc*proc + rankloc + 1
        call MPI_WAIT(req2(tag),status,ierr)
      end do
c
c     now do the actual communications
c
      proc = pneig_send(nneig_send)
      allocate (buffers(9,bufbeg3(proc+1)+buflen3(proc+1)))
      proc = pneig_recep(nneig_recep)
      allocate (buffer(9,bufbeg4(proc+1)+buflen4(proc+1)))
c
c     Begin reception of positions
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buffer(1,bufbeg4(proc+1)),9*buflen4(proc+1),
     $    MPI_REAL8,proc,tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        do j = 0, buflen3(proc+1)-1
          buffers(1,bufbeg3(proc+1)+j) = x(buf3(bufbeg3(proc+1)+j))
          buffers(2,bufbeg3(proc+1)+j) = y(buf3(bufbeg3(proc+1)+j))
          buffers(3,bufbeg3(proc+1)+j) = z(buf3(bufbeg3(proc+1)+j))
          buffers(4,bufbeg3(proc+1)+j) = v(1,buf3(bufbeg3(proc+1)+j))
          buffers(5,bufbeg3(proc+1)+j) = v(2,buf3(bufbeg3(proc+1)+j))
          buffers(6,bufbeg3(proc+1)+j) = v(3,buf3(bufbeg3(proc+1)+j))
          buffers(7,bufbeg3(proc+1)+j) = a(1,buf3(bufbeg3(proc+1)+j))
          buffers(8,bufbeg3(proc+1)+j) = a(2,buf3(bufbeg3(proc+1)+j))
          buffers(9,bufbeg3(proc+1)+j) = a(3,buf3(bufbeg3(proc+1)+j))
        end do
      end do
c
c     send the positions, the velocities and the accelerations
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buffers(1,bufbeg3(proc+1)),9*buflen3(proc+1),
     $   MPI_REAL8,proc,tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        tag = nprocloc*proc + rankloc + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_WAIT(req(tag),status,ierr)
        do j = 0, buflen4(proc+1)-1
          x(buf4(bufbeg4(proc+1)+j)) = buffer(1,bufbeg4(proc+1)+j)
          y(buf4(bufbeg4(proc+1)+j)) = buffer(2,bufbeg4(proc+1)+j)
          z(buf4(bufbeg4(proc+1)+j)) = buffer(3,bufbeg4(proc+1)+j)
          v(1,buf4(bufbeg4(proc+1)+j)) = buffer(4,bufbeg4(proc+1)+j)
          v(2,buf4(bufbeg4(proc+1)+j)) = buffer(5,bufbeg4(proc+1)+j)
          v(3,buf4(bufbeg4(proc+1)+j)) = buffer(6,bufbeg4(proc+1)+j)
          a(1,buf4(bufbeg4(proc+1)+j)) = buffer(7,bufbeg4(proc+1)+j)
          a(2,buf4(bufbeg4(proc+1)+j)) = buffer(8,bufbeg4(proc+1)+j)
          a(3,buf4(bufbeg4(proc+1)+j)) = buffer(9,bufbeg4(proc+1)+j)
        end do
      end do
c
c
c      deallocate (buffers)
c      deallocate (buffer)
      deallocate (req)
      deallocate (req2)
c      deallocate (buf4)
      deallocate (buflen4)
      deallocate (bufbeg4)
      return
      end
c
c     subroutine allocstep: deallocate arrays and reallocate them with proper size
c     (memory distribution)
c
      subroutine allocstep
      implicit none
      include 'deriv.i'
      include 'openmp.i'
c
      if (associated(torquedir)) deallocate (torquedir)
      allocate (torquedir(3,nbloc))
      torquedir = 0d0
      if (associated(torquedirp)) deallocate (torquedirp)
      allocate (torquedirp(3,nbloc))
      torquedirp = 0d0
      if (associated(desum)) deallocate (desum)
      allocate (desum(3,nbloc))
      desum = 0d0
      if (associated(torquerec)) deallocate (torquerec)
      allocate (torquerec(3,nblocrec))
      torquerec = 0d0
      if (associated(torquerecp)) deallocate (torquerecp)
      allocate (torquerecp(3,nblocrec))
      torquerecp = 0d0
      if (associated(demrec)) deallocate (demrec)
      allocate (demrec(3,nblocrec))
      demrec = 0d0
      if (associated(deprec)) deallocate (deprec)
      allocate (deprec(3,nblocrec))
      deprec = 0d0
c
      if (associated(deb)) deallocate (deb)
      allocate (deb(3,nbloc))
      if (associated(dea)) deallocate (dea)
      allocate (dea(3,nbloc))
      if (associated(deba)) deallocate (deba)
      allocate (deba(3,nbloc))
      if (associated(deub)) deallocate (deub)
      allocate (deub(3,nbloc))
      if (associated(deaa)) deallocate (deaa)
      allocate (deaa(3,nbloc))
      if (associated(deopb)) deallocate (deopb)
      allocate (deopb(3,nbloc))
      if (associated(deopd)) deallocate (deopd)
      allocate (deopd(3,nbloc))
      if (associated(deid)) deallocate (deid)
      allocate (deid(3,nbloc))
      if (associated(deit)) deallocate (deit)
      allocate (deit(3,nbloc))
      if (associated(det)) deallocate (det)
      allocate (det(3,nbloc))
      if (associated(dept)) deallocate (dept)
      allocate (dept(3,nbloc))
      if (associated(debt)) deallocate (debt)
      allocate (debt(3,nbloc))
      if (associated(dett)) deallocate (dett)
      allocate (dett(3,nbloc))
      if (associated(dev)) deallocate (dev)
      allocate (dev(3,nbloc))
      if (associated(dem)) deallocate (dem)
      allocate (dem(3,nbloc))
      if (associated(dep)) deallocate (dep)
      allocate (dep(3,nbloc))
      if (associated(deg)) deallocate (deg)
      allocate (deg(3,nbloc))
      if (associated(debond)) deallocate (debond)
      allocate (debond(3,nbloc))
c
      return
      end
c
c     subroutine allocstepsrespa: deallocate arrays and reallocate them with proper size
c     (memory distribution)
c
      subroutine allocsteprespa(fast)
      implicit none
      include 'deriv.i'
      include 'openmp.i'
      logical fast
c
      if (associated(desum)) deallocate (desum)
      allocate (desum(3,nbloc))
      if (associated(deb)) deallocate (deb)
      allocate (deb(3,nbloc))
      if (associated(dea)) deallocate (dea)
      allocate (dea(3,nbloc))
      if (associated(deba)) deallocate (deba)
      allocate (deba(3,nbloc))
      if (associated(deub)) deallocate (deub)
      allocate (deub(3,nbloc))
      if (associated(deaa)) deallocate (deaa)
      allocate (deaa(3,nbloc))
      if (associated(deopb)) deallocate (deopb)
      allocate (deopb(3,nbloc))
      if (associated(deopd)) deallocate (deopd)
      allocate (deopd(3,nbloc))
      if (associated(deid)) deallocate (deid)
      allocate (deid(3,nbloc))
      if (associated(deit)) deallocate (deit)
      allocate (deit(3,nbloc))
      if (associated(det)) deallocate (det)
      allocate (det(3,nbloc))
      if (associated(dept)) deallocate (dept)
      allocate (dept(3,nbloc))
      if (associated(debt)) deallocate (debt)
      allocate (debt(3,nbloc))
      if (associated(dett)) deallocate (dett)
      allocate (dett(3,nbloc))
      if (associated(torquedir)) deallocate (torquedir)
      allocate (torquedir(3,nbloc))
      if (associated(torquedirp)) deallocate (torquedirp)
      allocate (torquedirp(3,nbloc))
      if (associated(torquerec)) deallocate (torquerec)
      allocate (torquerec(3,nblocrec))
      if (associated(torquerecp)) deallocate (torquerecp)
      allocate (torquerecp(3,nblocrec))
      if (associated(demrec)) deallocate (demrec)
      allocate (demrec(3,nblocrec))
      if (associated(deprec)) deallocate (deprec)
      allocate (deprec(3,nblocrec))
      if (associated(debond)) deallocate (debond)
      allocate (debond(3,nbloc))
c
      if (.not.(fast)) then
        torquedir = 0d0
        torquedirp = 0d0
        torquerec = 0d0
        torquerecp = 0d0
        demrec = 0d0
        deprec = 0d0
      end if
c
      if (associated(dev)) deallocate (dev)
      allocate (dev(3,nbloc))
      if (associated(dem)) deallocate (dem)
      allocate (dem(3,nbloc))
      if (associated(dep)) deallocate (dep)
      allocate (dep(3,nbloc))
      if (associated(deg)) deallocate (deg)
      allocate (deg(3,nbloc))
c
      return
      end
c
c     subroutine sendall : everybody gets all the positions, speed and accel of the system
c
      subroutine sendall
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'moldyn.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,iglob,iloc,iproc,tag,ierr
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: req(:)
      real*8, allocatable :: postemp(:,:)
c
      allocate (req(nproc*nproc))
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
     $   MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_WAIT(req(tag),status,ierr)
        tag = nproc*iproc + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
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
     $   iproc,tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nproc*iproc + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     get the positions
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(postemp(1,bufbeg(iproc+1)),9*domlen(iproc+1),
     $   MPI_REAL8,iproc,tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,9*nloc,MPI_REAL8,iproc,tag,
     $    MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nproc*iproc + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
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
      deallocate (req)
      deallocate (postemp)
      return
      end
c
c     subroutine ddnumber : get the number of subdivision along each axis for
c     3d spatial decomposition
c
      subroutine ddnumber(num,nx,ny,nz,istep)
      implicit none
      include 'boxes.i'
      include 'inform.i'
      include 'iounit.i'
      include 'openmp.i'
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
c     subroutine reassignrespa: deal with change of domain with respa integrator
c
      subroutine reassignrespa(fast,ialt,nalt)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'openmp.i'
      include 'potent.i'
      include 'mpif.h'
      include 'neigh.i'
      integer i
      integer iproc,ierr,iglob
      integer iloc
      real*8 xr,yr,zr
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reparttemp(:)
      integer, allocatable :: req(:)
      integer, allocatable :: count2(:)
      integer, allocatable :: buf3(:), buflen3(:), bufbeg3(:)
      integer, allocatable :: domlentemp(:)
      integer nprocloc
      integer ialt,nalt
      logical fast
c
      if (use_pmecore) then
        nprocloc = ndir
      else
        nprocloc = nproc
      end if
      allocate (domlentemp(nproc))
c
      allocate (req(nproc*nproc))
c      allocate (reparttemp(nblocrecdir))
      allocate (reparttemp(n))
      reparttemp = 0
      domlentemp = domlen
c
      allocate (count2(nproc))
      allocate (buflen3(nproc))
      allocate (bufbeg3(nproc))
      count2 = 0
      buflen3 = 0
      bufbeg3 = 0
c
c     recompute number of atoms per domain
c
      nlocold = nloc
      nblocold = nbloc
      do i = 1, nloc
        iglob = glob(i)
        xr = x(iglob)
        yr = y(iglob)
        zr = z(iglob)
        if (use_bounds) call image(xr,yr,zr)
        do iproc = 0, nprocloc-1
          if ((zr.ge.zbegproc(iproc+1)).and.
     $     (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
     $    .and.(yr.lt.yendproc(iproc+1)).and.(xr.ge.xbegproc(iproc+1))
     $    .and.(xr.lt.xendproc(iproc+1))) then
             reparttemp(i) = iproc
          end if
        end do
      end do
c
c     every process sends his portion of reparttemp
c
      do i = 1, nneig_recep
        tag = nproc*rank + pneig_recep(i) + 1
        call MPI_IRECV(reparttemp(bufbeg(pneig_recep(i)+1)),
     $   domlen(pneig_recep(i)+1),MPI_INT,pneig_recep(i),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     begin sending
c
      do i = 1, nneig_send
        tag = nproc*pneig_send(i) + rank + 1
        call MPI_ISEND(reparttemp,
     $  domlen(rank+1),MPI_INT,pneig_send(i),tag,MPI_COMM_WORLD,
     $  req(tag),ierr)
      end do
c
      do i = 1, nneig_recep
        tag = nproc*rank + pneig_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nneig_send
        tag = nproc*pneig_send(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
      if ((.not.(fast))) then
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
          tag = nproc*rank + precdir_recep(i) + 1
          call MPI_IRECV(reparttemp(bufbeg(precdir_recep(i)+1)),
     $     domlen(precdir_recep(i)+1),MPI_INT,precdir_recep(i),tag,
     $     MPI_COMM_WORLD,req(tag),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
          tag = nproc*precdir_send(i) + rank + 1
          call MPI_ISEND(reparttemp,
     $    domlen(rank+1),MPI_INT,precdir_send(i),tag,MPI_COMM_WORLD,
     $    req(tag),ierr)
          end if
        end do
c
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
          tag = nproc*rank + precdir_recep(i) + 1
          call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
          tag = nproc*precdir_send(i) + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
          end if
        end do
      end if
c
      domlen = 0
c
      do i = 1, nloc
        iglob = glob(i)
        repart(iglob) = reparttemp(i)
        domlen(repart(iglob)+1) = domlen(repart(iglob)+1)+1
      end do
c
      do iproc = 1, nneig_recep
        do i = 1, domlentemp(pneig_recep(iproc)+1)
          iloc = bufbeg(pneig_recep(iproc)+1)+i-1
          iglob = glob(iloc)
          repart(iglob) = reparttemp(iloc)
          domlen(repart(iglob)+1) = domlen(repart(iglob)+1)+1
        end do
      end do
      if ((.not.(fast))) then
        do iproc = 1, nrecdir_recep
          do i = 1, domlentemp(precdir_recep(iproc)+1)
            iloc = bufbeg(precdir_recep(iproc)+1)+i-1
            iglob = glob(iloc)
            repart(iglob) = reparttemp(iloc)
          end do
        end do
      end if
c
c    get the positions to send (change of domains)
c
      do i = 1, nloc
        iglob = glob(i)
        if (rank.ne.repart(iglob)) then
          buflen3(repart(iglob)+1) = buflen3(repart(iglob)+1)+1
        end if
      end do
      bufbeg3(1) = 1
      do iproc = 1, nproc-1
        bufbeg3(iproc+1) = bufbeg3(iproc)+buflen3(iproc)
      end do
c
      allocate (buf3(bufbeg3(nproc)+buflen3(nproc)))
      buf3 = 0
c
      do i = 1, nloc
        iglob = glob(i)
        if (rank.ne.repart(iglob)) then
          buf3(bufbeg3(repart(iglob)+1)+count2(repart(iglob)+1))=
     $      iglob
          count2(repart(iglob)+1) = count2(repart(iglob)+1) + 1
        end if
      end do
c
      call orderbufferrespa(.false.,(ialt.ne.nalt),domlentemp)
      nloc = domlen(rank+1)
c
c
c     send the positions, the speed and the accelerations
c
      if ((.not.(use_pmecore)).or.(rank.le.ndir-1)) then
        call commreassign(buflen3,bufbeg3,buf3)
      end if
c
      deallocate (domlentemp)
      deallocate (count2)
      deallocate (buf3)
      deallocate (buflen3)
      deallocate (bufbeg3)
      deallocate (reparttemp)
      deallocate (req)
      return
      end
c
c     subroutine commposrespa : deal with communications of positions between two time steps
c     with respa integrator
c
      subroutine commposrespa(fast)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'potent.i'
      include 'mpif.h'
      integer i,tag,ierr,iproc
      integer iloc,iglob
      integer, allocatable :: req(:)
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
      allocate (req(nproc*nproc))
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
     $   MPI_REAL8,ptemp_recep(i),tag,MPI_COMM_WORLD,req(tag),ierr)
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
     $   req(tag),ierr)
      end do
c
      do i = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, ntemp_send
        tag = nproc*ptemp_send(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
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
      deallocate (req)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine sendallpos : everybody gets all the positions
c
      subroutine sendallpos
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,iglob,iloc,iproc,tag,ierr
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: req(:)
      real*8, allocatable :: postemp(:,:)
c
      allocate (req(nproc*nproc))
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
     $   MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_WAIT(req(tag),status,ierr)
        tag = nproc*iproc + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
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
     $   iproc,tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nproc*iproc + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     get the positions
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(postemp(1,bufbeg(iproc+1)),3*domlen(iproc+1),
     $   MPI_REAL8,iproc,tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,3*nloc,MPI_REAL8,iproc,tag,
     $    MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nproc*iproc + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
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
      deallocate (req)
      deallocate (postemp)
      return
      end
c
c     subroutine sendvecmin : deal with communications during minimization routine
c
      subroutine sendvecmin(xx)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,iglob,iloc,iproc,tag,ierr,j
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: req(:)
      real*8, allocatable :: temp(:,:)
      real*8 xx(3*n)
c
      allocate (req(nproc*nproc))
      req = MPI_REQUEST_NULL
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
     $   MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_WAIT(req(tag),status,ierr)
        tag = nproc*iproc + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
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
     $   iproc,tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nproc*iproc + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
        end if
      end do
c
c     get the values of xx
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(temp(1,bufbeg(iproc+1)),3*domlen(iproc+1),
     $   MPI_REAL8,iproc,tag,MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(temp,3*nloc,MPI_REAL8,iproc,tag,
     $    MPI_COMM_WORLD,req(tag),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_WAIT(req(tag),status,ierr)
          tag = nproc*iproc + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
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
      deallocate (req)
      deallocate (temp)
      return
      end
