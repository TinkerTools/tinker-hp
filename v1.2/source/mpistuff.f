c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
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
     $    MPI_REAL8,precdir_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $    MPI_REAL8,precdir_send(i),tag,COMM_TINKER,
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,
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
c     subroutine commposshort : communicate positions to short range neighboring processes  
c
      subroutine commposshort
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
      do i = 1, nbigshort_recep
        tag = nproc*rank + pbigshort_recep(i) + 1
        call MPI_IRECV(buffer(1,bufbeg(pbigshort_recep(i)+1)),
     $   3*domlen(pbigshort_recep(i)+1),
     $   MPI_REAL8,pbigshort_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
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
      do i = 1, nbigshort_send
        tag = nproc*pbigshort_send(i) + rank + 1
        call MPI_ISEND(buffers,3*nloc,
     $   MPI_REAL8,pbigshort_send(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbigshort_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbigshort_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, nbigshort_recep
        do i = 1, domlen(pbigshort_recep(iproc)+1)
          iloc = bufbeg(pbigshort_recep(iproc)+1)+i-1
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
        commloc  = COMM_TINKER
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
     $    MPI_REAL8,proc,tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,proc,tag,COMM_TINKER,reqsend(i),ierr)
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
     $      MPI_REAL8,proc,tag,COMM_TINKER,reqrec(i),ierr)
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
     $     MPI_REAL8,proc,tag,COMM_TINKER,reqsend(i),ierr)
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
      repart = -1
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
     $   COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqsend(iproc+1),ierr)
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
     $   iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
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
     $   MPI_REAL8,iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,9*nloc,MPI_REAL8,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
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
c     subroutine reassignrespa: deal with particles changing of domain between two time steps
c
c     subroutine reassignrespa(fast,ialt,naltloc)
      subroutine reassignrespa(ialt,naltloc)
      use atoms
      use bound
      use cell
      use domdec
      use moldyn
      use neigh
      use potent
      use mpi
      implicit none
      integer i,iproc,ierr,iglob
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
      integer ialt,naltloc

      eps1 = 1.0d-10
      eps2 = 1.0d-8
c
      if (ialt.ne.naltloc) return
c
      if (use_pmecore) then
        nprocloc = ndir
        commloc  = comm_dir
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
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
     $    MPI_REAL8,proc,tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,proc,tag,COMM_TINKER,reqsend(i),ierr)
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
      repart = -1
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
     $   MPI_REAL8,ptemp_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,ptemp_send(i),tag,COMM_TINKER,
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
     $   COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqsend(iproc+1),ierr)
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
     $   iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
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
     $   MPI_REAL8,iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,3*nloc,MPI_REAL8,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
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
     $   COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqsend(iproc+1),ierr)
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
     $   iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
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
     $   MPI_REAL8,iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(temp,3*nloc,MPI_REAL8,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
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
c
c     subroutine commforces : communicate forces between processes
c
      subroutine commforces(derivs)
      use timestat
      use mpi
      implicit none
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
      use timestat
      use mpi
      implicit none

      real*8 derivs(3,*)
      real*8 time0,time1,time2
      logical fast
c
      if (fast) then
        time0 = mpi_wtime()
        call commforcesbonded(derivs)
        time1 = mpi_wtime()
        timeforcescomm = timeforcescomm + time1-time0
      else
        time0 = mpi_wtime()
        call commforcesrec(derivs)
        time1 = mpi_wtime()
        timedirreccomm = timedirreccomm + time1-time0
        call commforcesbloc(derivs)
        time2 = mpi_wtime()
        timeforcescomm = timeforcescomm + time2-time0
      end if
      return
      end
c
c     subroutine commforcesrespa1 : communicate some forces between processes
c     when respa1 integrator is used
c
      subroutine commforcesrespa1(derivs,rule)
      use domdec
      use iounit
      use timestat
      use mpi
      implicit none
      integer rule
      real*8 derivs(3,*)
      real*8 time0,time1,time2
 1000 format(' illegal rule in commforcesrespa1.')
c
c     rule = 0: fast part of the forces
c     rule = 1: intermediate part of the forces
c     rule = 2: slow part of the forces
c
      if (rule.eq.0) then
        time0 = mpi_wtime()
        call commforcesbonded(derivs)
        time1 = mpi_wtime()
        timeforcescomm = timeforcescomm + time1-time0
      else if (rule.eq.1) then
        time0 = mpi_wtime()
        call commforcesshortbloc(derivs)
        time1 = mpi_wtime()
        timeforcescomm = timeforcescomm + time1-time0
      else if (rule.eq.2) then
        time0 = mpi_wtime()
        call commforcesrec(derivs)
        time1 = mpi_wtime()
        timedirreccomm = timedirreccomm + time1-time0
        call commforcesbloc(derivs)
        time2 = mpi_wtime()
        timeforcescomm = timeforcescomm + time2-time0
      else 
         if (rank.eq.0) write(iout,1000) 
      end if
      return
      end
c
c     subroutine reduceen : make the necessary reductions over the processes to get the
c     total values of the energies
c
      subroutine reduceen(epot)
      use domdec
      use energi
      use inter
      use potent
      use virial
      use mpi
      implicit none
      integer ierr,commloc
      real*8 epot
      real*8 buffer1(6),buffer2(17)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      buffer1(1) = ec
      buffer1(2) = em
      buffer1(3) = ep
      buffer1(4) = epot
      buffer1(5) = esum
      buffer1(6) = einter 
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer1,6,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      else
        call MPI_REDUCE(buffer1,buffer1,6,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      end if
      ec = buffer1(1)
      em = buffer1(2)
      ep = buffer1(3)
      epot = buffer1(4)
      esum = buffer1(5)
      einter = buffer1(6)
c
c   MPI: get virial
c
      call MPI_ALLREDUCE(MPI_IN_PLACE,vir,9,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
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
      buffer2(16) = ex
      buffer2(17) = eat
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer2,17,MPI_REAL8,MPI_SUM,0,
     $     commloc,ierr)
      else
        call MPI_REDUCE(buffer2,buffer2,17,MPI_REAL8,MPI_SUM,0,
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
      ex   = buffer2(16)
      eat  = buffer2(17)
c
      return
      end
c
c     subroutine commforcesbonded : deal with communications of some
c     bonded forces modifier avec nlocnl ?
c
      subroutine commforcesbonded(derivs)
      use sizes
      use atoms
      use deriv
      use domdec
      use potent
      use mpi
      implicit none
      integer i,j,tag,ierr,iproc
      integer iloc
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
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
     $   MPI_REAL8,pneig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pneig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nneig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nneig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
      deallocate (reqsend)
      deallocate (reqrec)
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
      use sizes
      use atoms
      use deriv
      use domdec
      use potent
      use mpi
      implicit none
      integer i,tag,ierr
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        derivs(:,1:nloc) = derivs(:,1:nloc) + buffer(:,1:nloc,i)
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commforcesshortbloc : deal with communications of short range forces by bloc
c     between two time steps
c
      subroutine commforcesshortbloc(derivs)
      use sizes
      use atoms
      use deriv
      use domdec
      use potent
      use mpi
      implicit none
      integer i,tag,ierr
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if ((use_pmecore).and.(rank.gt.ndir-1)) return
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(3,max(1,nloc),nbigshort_send))
c      buffer = 0d0
      allocate (buffers(3,nbloc))
c      buffers = 0d0
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nbigshort_send
        tag = nproc*rank + pbigshort_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbigshort_send(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = desum(:,(nloc+1):nbloc)
c
c     MPI : begin sending
c
      do i = 1, nbigshort_recep
        tag = nproc*pbigshort_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbigshort_recep(i)+1)),
     $   3*domlen(pbigshort_recep(i)+1),
     $   MPI_REAL8,pbigshort_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbigshort_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbigshort_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbigshort_send
        derivs(:,1:nloc) = derivs(:,1:nloc) + buffer(:,1:nloc,i)
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commforcesrec : deal with communications of reciprocal forces
c
      subroutine commforcesrec(derivs)
      use sizes
      use deriv
      use domdec
      use mpi
      implicit none
      integer i,j,tag,ierr,iproc,iglob,ilocrec
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 derivs(3,*)
      integer status(MPI_STATUS_SIZE)
c
c     allocate some arrays
c
c      allocate (req(nproc*nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
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
     $   MPI_REAL8,precdir_send(i),tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
c
c     MPI : move in buffer
c
      buffers(:,(1):nblocrec)= demrec(:,(1):nblocrec) + 
     $  deprec(:,(1):nblocrec)+decrec(:,(1):nblocrec)
c
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
        tag = nproc*precdir_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $   3*domlen(precdir_recep(i)+1),
     $   MPI_REAL8,precdir_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
        call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
        call MPI_WAIT(reqsend(i),status,ierr)
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
            derivs(1,j) =  derivs(1,j)       + demrec(1,ilocrec)
     $                   + deprec(1,ilocrec) + decrec(1,ilocrec)
            derivs(2,j) =  derivs(2,j)       + demrec(2,ilocrec)
     $                   + deprec(2,ilocrec) + decrec(2,ilocrec)
            derivs(3,j) =  derivs(3,j)       + demrec(3,ilocrec)
     $                   + deprec(3,ilocrec) + decrec(3,ilocrec)
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
     $      MPI_REAL8,precdir_send(i),tag,COMM_TINKER,reqrec(i),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep(iproc)+1)
              ilocrec = bufbegrec(precdir_recep(iproc)+1)+i-1
              buffers(1,ilocrec) =demrec(1,ilocrec)
              buffers(2,ilocrec) =demrec(2,ilocrec)
              buffers(3,ilocrec) =demrec(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $       3*domlen(precdir_recep(i)+1),
     $       MPI_REAL8,precdir_recep(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
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
              dem(1,j) = dem(1,j)+demrec(1,ilocrec)
              dem(2,j) = dem(2,j)+demrec(2,ilocrec)
              dem(3,j) = dem(3,j)+demrec(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_REAL8,precdir_send(i),tag,COMM_TINKER,reqrec(i),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep(iproc)+1)
              ilocrec = bufbegrec(precdir_recep(iproc)+1)+i-1
              buffers(1,ilocrec)=deprec(1,ilocrec)
              buffers(2,ilocrec)=deprec(2,ilocrec)
              buffers(3,ilocrec)=deprec(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $       3*domlen(precdir_recep(i)+1),
     $       MPI_REAL8,precdir_recep(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
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
              dep(1,j) =dep(1,j)+deprec(1,ilocrec)
              dep(2,j) =dep(2,j)+deprec(2,ilocrec)
              dep(3,j) =dep(3,j)+deprec(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_REAL8,precdir_send(i),tag,COMM_TINKER,reqrec(i),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep
          if (precdir_recep(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep(iproc)+1)
              ilocrec = bufbegrec(precdir_recep(iproc)+1)+i-1
              buffers(1,ilocrec) = decrec(1,ilocrec)
              buffers(2,ilocrec) = decrec(2,ilocrec)
              buffers(3,ilocrec) = decrec(3,ilocrec)
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $       3*domlen(precdir_recep(i)+1),
     $       MPI_REAL8,precdir_recep(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep
          if (precdir_recep(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
          end if
        end do
c
c       MPI : move in global arrays
c
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            do j = 1, nloc
              iglob = glob(j)
              dec(1,j) = dec(1,j) + buffer(1,j,i)
              dec(2,j) = dec(2,j) + buffer(2,j,i)
              dec(3,j) = dec(3,j) + buffer(3,j,i)
            end do
          else
            do j = 1, nloc
              iglob = glob(j)
              ilocrec = locrec1(iglob)
              dec(1,j) = dec(1,j) + decrec(1,ilocrec)
              dec(2,j) = dec(2,j) + decrec(2,ilocrec)
              dec(3,j) = dec(3,j) + decrec(3,ilocrec)
            end do
          end if
        end do
      end  if
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commforceststgrad : deal with communications of forces one by one
c     for testgrad
c
      subroutine commforceststgrad(detot)
      use sizes
      use deriv
      use domdec
      use mpi
      implicit none
      integer i,j,tag,ierr
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:,:)
      real*8, allocatable :: buffers(:,:)
      real*8 detot(3,*)
      integer status(MPI_STATUS_SIZE)
c
c     allocate some arrays
c
c      allocate (req(nproc*nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqrec(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = deat(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          deat(1,j) = deat(1,j) + buffer(1,j,i)
          deat(2,j) = deat(2,j) + buffer(2,j,i)
          deat(3,j) = deat(3,j) + buffer(3,j,i)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqrec(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dex(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dex(1,j) = dex(1,j) + buffer(1,j,i)
          dex(2,j) = dex(2,j) + buffer(2,j,i)
          dex(3,j) = dex(3,j) + buffer(3,j,i)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers(:,(nloc+1):nbloc) = dec(:,(nloc+1):nbloc)
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $   3*domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nbig_send
        do j = 1, nloc
          dec(1,j) = dec(1,j) + buffer(1,j,i)
          dec(2,j) = dec(2,j) + buffer(2,j,i)
          dec(3,j) = dec(3,j) + buffer(3,j,i)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
     $   MPI_REAL8,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nbig_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nbig_recep
        call MPI_WAIT(reqsend(i),status,ierr)
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
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine allreduceen : make the necessary allreductions over the processes to get the
c     total values of the energies
c
      subroutine allreduceen(epot)
      use domdec
      use energi
      use inter
      use mpi
      implicit none
      integer ierr
      real*8 epot
c
       call MPI_ALLREDUCE(MPI_IN_PLACE,eba,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ea,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eb,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ec,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,em,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ep,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ev,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eub,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eaa,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eopb,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eopd,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eid,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eit,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,et,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ept,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ebt,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eat,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ett,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,esum,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,epot,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eg,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ex,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,einter,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
      return
      end
c
c     subroutine orderbuffer : get the indexes of the particules in the neighboring processes
c
      subroutine orderbuffer(init)
      use sizes
      use atoms
      use cutoff
      use domdec
      use mpi
      implicit none
      integer, allocatable :: counttemp(:),ind1temp(:)
      integer count
      integer i,iproc,iproc1,tag,ierr

      integer, allocatable :: reqsend(:),reqrec(:)
      integer status(MPI_STATUS_SIZE)
      logical init
c
      allocate (counttemp(nproc))
      allocate (ind1temp(n))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
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
     $   precdir_recep(iproc),tag,COMM_TINKER,reqrec(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
        tag = nproc*precdir_send(iproc) + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,precdir_send(iproc),
     $   tag,COMM_TINKER,reqsend(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
        call MPI_WAIT(reqrec(iproc),status,ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
        call MPI_WAIT(reqsend(iproc),status,ierr)
        end if
      end do
c
c     get the domlen values of the coupled real space processes
c
      do iproc = 1, nbig_recep
        tag = nproc*rank + pbig_recep(iproc) + 1
        call MPI_IRECV(domlen(pbig_recep(iproc)+1),1,MPI_INT,
     $   pbig_recep(iproc),tag,COMM_TINKER,reqrec(iproc),ierr)
      end do
      do iproc = 1, nbig_send
        tag = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,pbig_send(iproc),tag,
     $   COMM_TINKER,reqsend(iproc),ierr)
      end do
      do iproc = 1, nbig_recep
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nbig_send
        call MPI_WAIT(reqsend(iproc),status,ierr)
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
     $   COMM_TINKER,reqrec(iproc),ierr)
      end do
      do iproc = 1, nbig_send
        tag = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,pbig_send(iproc),
     $   tag,COMM_TINKER,reqsend(iproc),ierr)
      end do
      do iproc = 1, nbig_recep
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nbig_send
        call MPI_WAIT(reqsend(iproc),status,ierr)
      end do
c
      do iproc = 1, nbig_recep
        do i = 1, domlen(pbig_recep(iproc)+1)
          loc(glob(bufbeg(pbig_recep(iproc)+1)+i-1)) =
     $      bufbeg(pbig_recep(iproc)+1)+i-1
          repart(glob(bufbeg(pbig_recep(iproc)+1)+i-1)) = 
     $     pbig_recep(iproc)
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
     $     precdir_recep(iproc),tag,COMM_TINKER,reqrec(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
          tag = nproc*precdir_send(iproc) + rank + 1
          call MPI_ISEND(glob,domlen(rank+1),MPI_INT,
     $     precdir_send(iproc),tag,COMM_TINKER,reqsend(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
          call MPI_WAIT(reqrec(iproc),status,ierr)
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
          call MPI_WAIT(reqsend(iproc),status,ierr)
        end if
      end do
c
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
          do i = 1, domlen(precdir_recep(iproc)+1)
            loc(glob(bufbeg(precdir_recep(iproc)+1)+i-1)) =
     $        bufbeg(precdir_recep(iproc)+1)+i-1
            repart(glob(bufbeg(precdir_recep(iproc)+1)+i-1)) = 
     $       precdir_recep(iproc)
          end do
        end if
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (ind1temp)
      deallocate (counttemp)
      return
      end
c
c
c     subroutine orderbufferrec : get the indexes of the particules in the reciprocal processes
c
      subroutine orderbufferrec
      use atoms
      use domdec
      use potent
      implicit none
      integer, allocatable :: counttemp(:),ind1temp(:)
      integer count,rankloc
      integer i,iproc
      integer iloc,iglob
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
c    subroutine commdirdir : deal with communications of direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles 
c        - 2: wait for the communications to be done
c
      subroutine commdirdir(nrhs,rule,mu,reqrec,reqsend)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
      integer :: reqrec(nproc),reqsend(nproc)
      real*8 mu(3,nrhs,npolebloc)
 1000 format(' illegal rule in commdirdir.')
c
      if (rule.eq.0) then
c
c     MPI : begin reception
c
        do i = 1, n_recep1
          tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   COMM_TINKER,reqrec(i),ierr)
        end do
      else if (rule.eq.1) then
c
c     MPI : begin sending
c
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
          call MPI_ISEND(mu,3*nrhs*npoleloc,
     $     MPI_REAL8,p_send1(i),tag,COMM_TINKER,
     $     reqsend(i),ierr)
        end do
      else if (rule.eq.2) then
        do i = 1, n_recep1
           call MPI_WAIT(reqrec(i),status,ierr)
         end do
         do i = 1, n_send1
           call MPI_WAIT(reqsend(i),status,ierr)
         end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c    subroutine commdirdirshort : deal with communications of short range direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles 
c        - 2: wait for the communications to be done
c
      subroutine commdirdirshort(nrhs,rule,mu,reqrec,reqsend)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
      integer :: reqrec(nproc),reqsend(nproc)
      real*8 mu(3,nrhs,npolebloc)
 1000 format(' illegal rule in commdirdir.')
c
      if (rule.eq.0) then
c
c     MPI : begin reception
c
        do i = 1, n_recepshort1
          tag = nproc*rank + p_recepshort1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recepshort1(i)+1)),
     $   3*nrhs*domlen(p_recepshort1(i)+1),MPI_REAL8,p_recepshort1(i),
     $   tag,COMM_TINKER,reqrec(i),ierr)
        end do
      else if (rule.eq.1) then
c
c     MPI : begin sending
c
        do i = 1, n_sendshort1
          tag = nproc*p_sendshort1(i) + rank + 1
          call MPI_ISEND(mu,3*nrhs*npoleloc,
     $     MPI_REAL8,p_sendshort1(i),tag,COMM_TINKER,
     $     reqsend(i),ierr)
        end do
      else if (rule.eq.2) then
        do i = 1, n_recepshort1
           call MPI_WAIT(reqrec(i),status,ierr)
         end do
         do i = 1, n_sendshort1
           call MPI_WAIT(reqsend(i),status,ierr)
         end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c
c    subroutine commrecdir : deal with communications of  reciprocal fields
c
c    rule determines what to do:
c        - 0: start reception of reciprocal fields
c        - 1: start sending of reciprocal fields
c        - 2: wait for the communications to be done
c
      subroutine commrecdirfields(rule,efrec,ef,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: reqrecdirrec(nproc),reqrecdirsend(nproc)
      real*8 efrec(10,max(1,npolerecloc))
      real*8 ef(10,max(1,npoleloc))
      real*8 buffermpi1(10,max(npoleloc,1))
      real*8 buffermpi2(10,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdirfields.')
c
      if (rule.eq.0) then
c
c     Begin the reception of the reciprocal fields
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpi1(1,bufbeg1(proc+1)),
     $        10*buflen1(proc+1),
     $        MPI_REAL8,proc,tag,COMM_TINKER,reqrecdirrec(i),ierr)
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            do j = 0, buflen2(proc+1)-1
              call amove(10,efrec(1,polerecloc(
     $       buf2(bufbeg2(proc+1)+j))),buffermpi2(1,bufbeg2(proc+1)+j))
            end do
          end if
        end do
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpi2(1,bufbeg2(proc+1)),
     $       10*buflen2(proc+1),MPI_REAL8,proc,tag,COMM_TINKER,
     $       reqrecdirsend(i),ierr)
          end if
        end do
      else if (rule.eq.2) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirsend(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirrec(i),status,ierr)
            do j = 0, buflen1(proc+1)-1
              call amove(10,buffermpi1(1,bufbeg1(proc+1)+j),
     $          ef(1,poleloc(buf1(bufbeg1(proc+1)+j))))
            end do
          end if
        end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c    subroutine commrecdirsolv : deal with communications of reciprocal fields in the solver
c
c    rule determines what to do:
c        - 0: start reception of reciprocal fields
c        - 1: start sending of reciprocal fields
c        - 2: wait for the communications to be done
c
      subroutine commrecdirsolv(nrhs,rule,dipfieldbis,dipfield,
     $ buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: reqrecdirrec(nproc),reqrecdirsend(nproc)
      real*8 dipfield(3,nrhs,max(1,npoleloc))
      real*8 dipfieldbis(3,nrhs,max(1,npolerecloc))
      real*8 buffermpi1(3,nrhs,max(npoleloc,1))
      real*8 buffermpi2(3,nrhs,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdirdip.')
c
      if (rule.eq.0) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
           tag = nproc*rank + proc + 1
           call MPI_IRECV(buffermpi1(1,1,bufbeg1(proc+1)),
     $       3*nrhs*buflen1(proc+1),
     $       MPI_REAL8,proc,tag,COMM_TINKER,reqrecdirrec(i),ierr)
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,dipfieldbis(1,1,polerecloc(
     $          buf2(bufbeg2(proc+1)+j))),
     $          buffermpi2(1,1,bufbeg2(proc+1)+j))
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpi2(1,1,bufbeg2(proc+1)),
     $       3*nrhs*buflen2(proc+1),
     $       MPI_REAL8,proc,tag,COMM_TINKER,reqrecdirsend(i),ierr)
          end if
        end do
      else if (rule.eq.2) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirsend(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirrec(i),status,ierr)
            do j = 0, buflen1(proc+1)-1
              call amove(3*nrhs,buffermpi1(1,1,bufbeg1(proc+1)+j),
     $          dipfield(1,1,poleloc(buf1(bufbeg1(proc+1)+j))))
            end do
          end if
        end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c
c    subroutine commrecdirdipsolv : deal with communications of dipoles for PME 
c
c    rule determines what to do:
c        - 0: start reception of dipoles
c        - 1: start sending of dipoles 
c        - 2: wait for the communications to be done
c
      subroutine commrecdirdip(nrhs,rule,diprec,dip,buffermpimu1,
     $ buffermpimu2,req2rec,req2send)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: req2rec(nproc),req2send(nproc)
      real*8 dip(3,nrhs,max(1,npoleloc))
      real*8 diprec(3,nrhs,max(1,npolerecloc))
c      real*8 buffermpimu1(3,nrhs,max(npoleloc,1))
c      real*8 buffermpimu2(3,nrhs,max(1,npolerecloc))
      real*8 buffermpimu1(3,nrhs,*)
      real*8 buffermpimu2(3,nrhs,*)
 1000 format(' illegal rule in commrecdirdip.')
c
      if (rule.eq.0) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),3*nrhs*
     $       buflen2(proc+1),MPI_REAL8,proc,tag,COMM_TINKER,
     $       req2rec(i),ierr)
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            do j = 0, buflen1(proc+1)-1
              call amove(3*nrhs,dip(1,1,poleloc(
     $          buf1(bufbeg1(proc+1)+j))),
     $          buffermpimu1(1,1,bufbeg1(proc+1)+j))
            end do
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),3*nrhs*
     $       buflen1(proc+1),MPI_REAL8,proc,tag,COMM_TINKER,
     $       req2send(i),ierr)
          end if
        end do
      else if (rule.eq.2) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(req2send(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(req2rec(i),status,ierr)
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $          diprec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
            end do
          end if
        end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c
c     subroutine commfield : communicate some direct fields (Newton's third law)
c
      subroutine commfield(nrhs,ef)
      use atoms
      use domdec
      use mpole
      use potent
      use mpi
      implicit none
      integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer, allocatable :: reqrec(:), reqsend(:)
      real*8 ef(3,nrhs,*)
      real*8, allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_REAL8,
     $   p_send1(i),tag,commloc,reqrec(i),ierr)
      end do
c
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlenpole(p_recep1(i)+1),MPI_REAL8,p_recep1(i),
     $  tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
      do i = 1, n_send1
        do j = 1, npoleloc
          do k = 1, nrhs
            ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
            ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
            ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
          end do
        end do
      end do
c
      deallocate (buffer)
      deallocate (reqrec)
      deallocate (reqsend)
      return
      end
c
c     subroutine commfieldshort : communicate some short range direct fields (Newton's third law)
c
      subroutine commfieldshort(nrhs,ef)
      use atoms
      use domdec
      use mpole
      use potent
      use mpi
      implicit none
      integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer, allocatable :: reqrec(:), reqsend(:)
      real*8 ef(3,nrhs,*)
      real*8, allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      allocate (buffer(3,nrhs,max(npoleloc,1),n_sendshort1))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
      do i = 1, n_sendshort1
        tag = nproc*rank + p_sendshort1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_REAL8,
     $   p_sendshort1(i),tag,commloc,reqrec(i),ierr)
      end do
c
      do i = 1, n_recepshort1
        tag = nproc*p_recepshort1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recepshort1(i)+1)),
     $   3*nrhs*domlenpole(p_recepshort1(i)+1),MPI_REAL8,
     $   p_recepshort1(i),tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, n_sendshort1
        tag = nproc*rank + p_sendshort1(i) + 1
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_recepshort1
        tag = nproc*p_recepshort1(i) + rank + 1
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
      do i = 1, n_sendshort1
        do j = 1, npoleloc
          do k = 1, nrhs
            ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
            ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
            ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
          end do
        end do
      end do
c
      deallocate (buffer)
      deallocate (reqrec)
      deallocate (reqsend)
      return
      end
