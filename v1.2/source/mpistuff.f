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
      use potent
      use mpi
      implicit none
      integer i,iproc
      integer iglob,iloc
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:),buffer2(:,:)
      real*8, allocatable :: buffers(:,:),buffers2(:,:)
      integer status(MPI_STATUS_SIZE)
      integer tag,ierr
      integer rankloc,commloc,nprocloc
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(3,nblocrecdir))
      buffer = 0d0
      allocate (buffer2(3,nlocrec2))
      buffer2 = 0d0
      allocate (buffers(3,nloc))
      buffers = 0d0
      allocate (buffers2(3,nlocrec))
      buffers2 = 0d0
c
c     first do dir rec communications
c
c     MPI : begin reception in buffer
c

      do i = 1, nrecdir_recep2
        if (precdir_recep2(i).ne.rank) then
         tag = nproc*rank + precdir_recep2(i) + 1
         call MPI_IRECV(buffer(1,bufbeg(precdir_recep2(i)+1)),
     $    3*domlen(precdir_recep2(i)+1),
     $    MPI_REAL8,precdir_recep2(i),tag,COMM_TINKER,reqrec(i),ierr)
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
      do i = 1, nrecdir_send2
        if (precdir_send2(i).ne.rank) then
         tag = nproc*precdir_send2(i) + rank + 1
         call MPI_ISEND(buffers,3*nloc,
     $    MPI_REAL8,precdir_send2(i),tag,COMM_TINKER,
     $    reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nrecdir_recep2
        if (precdir_recep2(i).ne.rank) then
         call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nrecdir_send2
        if (precdir_send2(i).ne.rank) then
         call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc).ne.rank) then
         do i = 1, domlen(precdir_recep2(iproc)+1)
           iloc = bufbeg(precdir_recep2(iproc)+1)+i-1
           iglob = glob(iloc)
           x(iglob) = buffer(1,iloc)
           y(iglob) = buffer(2,iloc)
           z(iglob) = buffer(3,iloc)
         end do
        end if
      end do
c
c     then do rec rec communications
c
c
c     MPI : begin reception in buffer
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(buffer2(1,bufbegrec(prec_recep(i)+1)),
     $   3*domlenrec(prec_recep(i)+1),
     $   MPI_REAL8,prec_recep(i),tag,commloc,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      do i = 1, nlocrec
        iglob = globrec(i)
        buffers2(1,i) = x(iglob)
        buffers2(2,i) = y(iglob)
        buffers2(3,i) = z(iglob)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(buffers2,3*nlocrec,
     $   MPI_REAL8,prec_send(i),tag,commloc,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, nrec_recep
        do i = 1, domlenrec(prec_recep(iproc)+1)
          iloc = bufbegrec(prec_recep(iproc)+1)+i-1
          iglob = globrec(iloc)
          x(iglob) = buffer2(1,iloc)
          y(iglob) = buffer2(2,iloc)
          z(iglob) = buffer2(3,iloc)
        end do
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      deallocate (buffer2)
      deallocate (buffers2)
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
#ifdef COLVARS
      use colvars
#endif
      use domdec
      use freeze
      use mdstuf
      use moldyn
      use neigh
      use potent
      use mpi
      use qtb, only: qtb_thermostat, reassignqtb
      implicit none
      integer i,iproc,ierr,iglob
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
      integer max_atoms_recv,max_atoms_send
      integer, allocatable :: iglob_send(:,:), iglob_recv(:,:)
      integer, allocatable :: n_data_send(:), n_data_recv(:)
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
      allocate (buflen3(nprocloc))
      allocate (bufbeg3(nprocloc))
      allocate (buf3bis(nloc,nprocloc))
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
c        if (use_bounds) call image(xr,yr,zr)
        call image(xr,yr,zr)
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
      do iproc = 1, nprocloc-1
        bufbeg3(iproc+1) = bufbeg3(iproc)+buflen3(iproc)
      end do
c
      allocate (buf3(bufbeg3(nprocloc)+buflen3(nprocloc)))
      buf3 = 0
c
      do iproc = 1, nprocloc
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
      allocate(n_data_send(0:nneig_send))
      n_data_send(:)=0
      allocate(n_data_recv(nneig_recep))
      n_data_recv(:)=0
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
        n_data_recv(iproc)=buflen4(proc+1)
      end do
      do iproc = 1, nneig_send
        proc = pneig_recep(iproc)
        call MPI_WAIT(reqsend(iproc),status,ierr)
        n_data_send(iproc)=buflen3(proc+1)
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

      max_atoms_send=maxval(n_data_send)
      allocate(iglob_send(nloc,0:nneig_send))
      iglob_send(:,:)=0
      do i = 1, nneig_send
        proc = pneig_send(i)
        do j = 0, buflen3(proc+1)-1
          jglob = buf3(bufbeg3(proc+1)+j)
          iglob_send(j+1,i)=jglob
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

      max_atoms_recv=maxval(n_data_recv)
      allocate(iglob_recv(max_atoms_recv,nneig_recep))
      iglob_recv(:,:)=0
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        call MPI_WAIT(reqrec(i),status,ierr)
        do j = 0, buflen4(proc+1)-1
          iglob = int(buffer(10,bufbeg4(proc+1)+j))
          iglob_recv(j+1,i)=iglob
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
c     also send the pbcwrap indexes
c     (unnecessary since it is synchronized globally in mdsave)
c      call commpbcwrapindex(bufbeg3,buf3,buflen3,bufbeg4,buflen4)
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

#ifdef COLVARS
c
c     if mts+colvars is used, also send aalt and aalt2 values
c
      if (mts.and.use_colvars) then

        proc = pneig_send(nneig_send)
        allocate (buffers1(7,bufbeg3(proc+1)+buflen3(proc+1)))
        proc = pneig_recep(nneig_recep)
        allocate (buffer1(7,bufbeg4(proc+1)+buflen4(proc+1)))
c
c     Begin reception 
c
        do i = 1, nneig_recep
          proc = pneig_recep(i)
          tag = nprocloc*rankloc + proc + 1
          call MPI_IRECV(buffer1(1,bufbeg4(proc+1)),7*buflen4(proc+1),
     $      MPI_REAL8,proc,tag,COMM_TINKER,reqrec(i),ierr)
        end do
        do i = 1, nneig_send
          proc = pneig_send(i)
          do j = 0, buflen3(proc+1)-1
            jglob = buf3(bufbeg3(proc+1)+j)
            buffers1(1,bufbeg3(proc+1)+j) = jglob
            buffers1(2,bufbeg3(proc+1)+j) = aalt(1,jglob)
            buffers1(3,bufbeg3(proc+1)+j) = aalt(2,jglob)
            buffers1(4,bufbeg3(proc+1)+j) = aalt(3,jglob)
            buffers1(5,bufbeg3(proc+1)+j) = aalt2(1,jglob)
            buffers1(6,bufbeg3(proc+1)+j) = aalt2(2,jglob)
            buffers1(7,bufbeg3(proc+1)+j) = aalt2(3,jglob)
          end do
        end do
c
c       send the old positions
c
        do i = 1, nneig_send
          proc = pneig_send(i)
          tag = nprocloc*proc + rankloc + 1
          call MPI_ISEND(buffers1(1,bufbeg3(proc+1)),7*buflen3(proc+1),
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
            aalt(1,iglob) = buffer1(2,bufbeg4(proc+1)+j)
            aalt(2,iglob) = buffer1(3,bufbeg4(proc+1)+j)
            aalt(3,iglob) = buffer1(4,bufbeg4(proc+1)+j)
            aalt2(1,iglob) = buffer1(5,bufbeg4(proc+1)+j)
            aalt2(2,iglob) = buffer1(6,bufbeg4(proc+1)+j)
            aalt2(3,iglob) = buffer1(7,bufbeg4(proc+1)+j)
          end do
        end do
        deallocate (buffer1)
        deallocate (buffers1)
      end if
#endif

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
        do j = 1, bufbeg3(nprocloc)+buflen3(nprocloc)-1
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

      if(qtb_thermostat) call reassignqtb(nloc,max_atoms_recv
     &      ,nneig_recep, n_data_recv, iglob_recv, pneig_recep 
     &      ,nneig_send , n_data_send, iglob_send, pneig_send 
     &      ,max_atoms_recv, max_atoms_send)
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
      integer, intent(in) :: ialt,naltloc
      if (ialt.ne.naltloc) return
      call reassign
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
      real*8 time0,time1
c
      time0 = mpi_wtime()
      call commforcesrec(derivs)
      time1 = mpi_wtime()
      timecommforcesrec = timecommforcesrec + time1-time0
      time0 = mpi_wtime()
      call commforcesbloc(derivs)
      time1 = mpi_wtime()
      timecommforcesreal = timecommforcesreal + time1-time0
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
      real*8 time0,time1
      logical fast
c
      if (fast) then
        time0 = mpi_wtime()
        call commforcesbonded(derivs)
        time1 = mpi_wtime()
        timecommforcesreal = timecommforcesreal + time1-time0
      else
        time0 = mpi_wtime()
        call commforcesrec(derivs)
        time1 = mpi_wtime()
        timecommforcesrec = timecommforcesrec + time1-time0
        time0 = mpi_wtime()
        call commforcesbloc(derivs)
        time1 = mpi_wtime()
        timecommforcesreal = timecommforcesreal + time1-time0
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
      real*8 time0,time1
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
        timecommforcesreal = timecommforcesreal + time1-time0
      else if (rule.eq.1) then
        time0 = mpi_wtime()
        call commforcesshortbloc(derivs)
        time1 = mpi_wtime()
        timecommforcesreal = timecommforcesreal + time1-time0
      else if (rule.eq.2) then
        time0 = mpi_wtime()
        call commforcesrec(derivs)
        time1 = mpi_wtime()
        timecommforcesrec = timecommforcesrec + time1-time0
        time0 = mpi_wtime()
        call commforcesbloc(derivs)
        time1 = mpi_wtime()
        timecommforcesreal = timecommforcesreal + time1-time0
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
      real*8 buffer1(7),buffer2(20)
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
      if (use_dewald) then
        buffer1(7) = edsp
      else
        buffer1(7) = 0d0
      end if
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer1,7,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      else
        call MPI_REDUCE(buffer1,buffer1,7,MPI_REAL8,MPI_SUM,0,
     $     COMM_TINKER,ierr)
      end if
      ec = buffer1(1)
      em = buffer1(2)
      ep = buffer1(3)
      epot = buffer1(4)
      esum = buffer1(5)
      einter = buffer1(6)
      if (use_dewald) edsp = buffer1(7)
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
      buffer2(17) = er
      buffer2(18) = ect
      buffer2(19) = eat
      if (.not.use_dewald) then
        buffer2(20) = edsp
      else
        buffer2(20) = 0d0
      end if
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer2,20,MPI_REAL8,MPI_SUM,0,
     $     commloc,ierr)
      else
        call MPI_REDUCE(buffer2,buffer2,20,MPI_REAL8,MPI_SUM,0,
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
      er   = buffer2(17)
      ect  = buffer2(18)
      eat  = buffer2(19)
      if (.not.use_dewald) edsp = buffer2(20)
c      call info_energy(0)
c
      return
      end

      subroutine info_energy(rank)
      use energi
      implicit none
      integer,intent(in):: rank
      real(8) ebonded

      if (rank.eq.0) then
!$acc wait
 20   format ( 40('-'))
         print 20
 30   format (1x,A,F18.6)
!$acc update host(esum,eb,ea,eba,eub,eaa,eopb,eopd,
!$acc&                  eid,eit,et,ept,ebt,ett,eg,ex,eat,
!$acc&                  ev,ec,ecrec,em,emrec,ep,eprec,esave,ensmd,
!$acc&                  eDaMD,ePaMD,eW1aMD,eW2aMD,
!$acc&                  ev_r,ec_r,em_r,ep_r,eb_r)

         ebonded = 
     &   eb+ea+eba+eub+eaa+eid+eit+et+ept+ebt+ett+eat+eopb+eopd+eg+ex
         if (eb   /=real(0,8)) print 30, 'eb     = ',eb
         !if (eb_r /=tp2enr(0.0)) print 30, 'eb_r = ',enr2en(eb_r)
         if (ea   /=real(0,8)) print 30, 'ea     = ',ea 
         if (eba  /=real(0,8)) print 30, 'eba    = ',eba
         if (eub  /=real(0,8)) print 30, 'eub    = ',eub
         if (eaa  /=real(0,8)) print 30, 'eaa    = ',eaa
         if (eid  /=real(0,8)) print 30, 'eid    = ',eid
         if (eit  /=real(0,8)) print 30, 'eit    = ',eit
         if (et   /=real(0,8)) print 30, 'et     = ',et
         if (ept  /=real(0,8)) print 30, 'ept    = ',ept
         if (ebt  /=real(0,8)) print 30, 'ebt    = ',ebt
         if (ett  /=real(0,8)) print 30, 'ett    = ',ett
         if (eat  /=real(0,8)) print 30, 'eat    = ',ett
         if (eopb /=real(0,8)) print 30, 'eopb   = ',eopb
         if (eopd /=real(0,8)) print 30, 'eopd   = ',eopd
         if (eg   /=real(0,8)) print 30, 'eg     = ',eg
         if (ex   /=real(0,8)) print 30, 'ex     = ',ex
         if (ebonded/=real(0,8)) print 30, 'ebonded =',ebonded
         if (ec   /=real(0,8)) print 30, 'ec     = ',ec
         if (ev   /=real(0,8)) print 30, 'ev     = ',ev
         if (em   /=real(0,8)) print 30, 'em     = ',em
         if (ep   /=real(0,8)) print 30, 'ep     = ',ep
         !if (esum/=real(0,8)) print 30, 'esum    = ',esum
      end if
!$acc wait
      end subroutine
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
      use potent
      use mpi
      implicit none
      integer i,j,tag,ierr,iproc,iglob,iloc
      integer jloc,jglob,jlocrec
      integer rankloc,commloc,nprocloc
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:,:),buffer2(:,:,:)
      real*8, allocatable :: buffers(:,:),buffers2(:,:)
      real*8, allocatable :: temprec(:,:)
      real*8 derivs(3,*)
      integer status(MPI_STATUS_SIZE)

      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(3,max(1,nloc),nrecdir_send1))
      allocate (buffers(3,max(1,nblocrecdir)))
      allocate (buffer2(3,nlocrec2,nrec_send))
      allocate (buffers2(3,max(1,nlocrec2)))
      allocate (temprec(3,nlocrec2))
c
c     first do rec rec communications
c
c
c     MPI : begin reception in buffer
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(buffer2(1,1,i),3*nlocrec,
     $   MPI_REAL8,prec_send(i),tag,commloc,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      buffers2(:,(nlocrec+1):nlocrec2)= demrec(:,(nlocrec+1):nlocrec2)+
     $  deprec(:,(nlocrec+1):nlocrec2)+ decrec(:,(nlocrec+1):nlocrec2)+
     $  dedsprec(:,(nlocrec+1):nlocrec2)
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(buffers2(1,bufbegrec(prec_recep(i)+1)),
     $   3*domlenrec(prec_recep(i)+1),
     $   MPI_REAL8,prec_recep(i),tag,commloc,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, nrec_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      temprec = demrec + deprec + decrec + dedsprec 
      do i = 1, nrec_send
        do j = 1, nlocrec
          temprec(1,j) = temprec(1,j) + buffer2(1,j,i)
          temprec(2,j) = temprec(2,j) + buffer2(2,j,i)
          temprec(3,j) = temprec(3,j) + buffer2(3,j,i)
        end do
      end do

c
c     then do rec dir communications
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
        tag = nproc*rank + precdir_send1(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_REAL8,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
c
c     Move in array
c
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
          do j = 1, domlen(precdir_recep1(i)+1)
            jloc = bufbeg(precdir_recep1(i)+1)+j-1
            jglob = glob(jloc)
            if (repartrec(jglob).eq.rankloc) then
              jlocrec = locrec(jglob)
              buffers(1,jloc) = temprec(1,jlocrec)
              buffers(2,jloc) = temprec(2,jlocrec)
              buffers(3,jloc) = temprec(3,jlocrec)
            else
              buffers(1,jloc) = 0d0
              buffers(2,jloc) = 0d0
              buffers(3,jloc) = 0d0
            end if
          end do
        end if
      end do
c
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
        tag = nproc*precdir_recep1(i) + rank + 1
        call MPI_ISEND(buffers(1,bufbeg(precdir_recep1(i)+1)),
     $   3*domlen(precdir_recep1(i)+1),
     $   MPI_REAL8,precdir_recep1(i),tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
        call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        if (precdir_recep1(i).ne.rank) then
        call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nrecdir_send1
        if (precdir_send1(i).ne.rank) then
          do j = 1, nloc
            derivs(1,j) = derivs(1,j) + buffer(1,j,i)
            derivs(2,j) = derivs(2,j) + buffer(2,j,i)
            derivs(3,j) = derivs(3,j) + buffer(3,j,i)
          end do
        else
          do j = 1, nlocrec
            iglob = globrec(j)
            iloc = loc(iglob)
            if (repart(iglob).eq.rank) then
              derivs(1,iloc) = derivs(1,iloc) + temprec(1,j)
              derivs(2,iloc) = derivs(2,iloc) + temprec(2,j)
              derivs(3,iloc) = derivs(3,iloc) + temprec(3,j)
            end if
          end do
        end if
      end do
c
c    Communicate separately the direct and reciprocal torques if testgrad is
c    being used
c
      if (dotstgrad) then
c
c     first do rec rec communications
c
c     MPI : begin reception in buffer
c
        do i = 1, nrec_send
          tag = nprocloc*rankloc + prec_send(i) + 1
          call MPI_IRECV(buffer2(1,1,i),3*nlocrec,
     $     MPI_REAL8,prec_send(i),tag,commloc,reqrec(i),ierr)
        end do
c
c     MPI : move in buffer
c
        buffers2(:,(nlocrec+1):nlocrec2)=demrec(:,(nlocrec+1):nlocrec2)
c
        do i = 1, nrec_recep
          tag = nprocloc*prec_recep(i) + rankloc + 1
          call MPI_ISEND(buffers2(1,bufbegrec(prec_recep(i)+1)),
     $     3*domlenrec(prec_recep(i)+1),
     $     MPI_REAL8,prec_recep(i),tag,commloc,
     $     reqsend(i),ierr)
        end do
c
        do i = 1, nrec_send
          call MPI_WAIT(reqrec(i),status,ierr)
        end do
        do i = 1, nrec_recep
          call MPI_WAIT(reqsend(i),status,ierr)
        end do
c
c     MPI : move in global arrays
c
        temprec = demrec

        do i = 1, nrec_send
          do j = 1, nlocrec
            temprec(1,j) = temprec(1,j) + buffer2(1,j,i)
            temprec(2,j) = temprec(2,j) + buffer2(2,j,i)
            temprec(3,j) = temprec(3,j) + buffer2(3,j,i)
          end do
        end do
c
c     MPI : begin reception in buffer
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            tag = nproc*rank + precdir_send1(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_REAL8,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep1
          if (precdir_recep1(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep1(iproc)+1)
              jloc = bufbeg(precdir_recep1(iproc)+1)+i-1
              jglob = glob(jloc)
              if (repartrec(jglob).eq.rankloc) then
                jlocrec = locrec(jglob)
                buffers(1,jloc) = temprec(1,jlocrec)
                buffers(2,jloc) = temprec(2,jlocrec)
                buffers(3,jloc) = temprec(3,jlocrec)
              else
                buffers(1,jloc) = 0d0
                buffers(2,jloc) = 0d0
                buffers(3,jloc) = 0d0
              end if
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            tag = nproc*precdir_recep1(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbeg(precdir_recep1(i)+1)),
     $       3*domlen(precdir_recep1(i)+1),
     $       MPI_REAL8,precdir_recep1(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
          end if
        end do
c
c       MPI : move in global arrays
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            do j = 1, nloc
              iglob = glob(j)
              dem(1,j) = dem(1,j) + buffer(1,j,i)
              dem(2,j) = dem(2,j) + buffer(2,j,i)
              dem(3,j) = dem(3,j) + buffer(3,j,i)
            end do
          else
            do j = 1, nlocrec
              iglob = globrec(j)
              iloc = loc(iglob)
              if (repart(iglob).eq.rank) then
                dem(1,iloc) = dem(1,iloc)+temprec(1,j)
                dem(2,iloc) = dem(2,iloc)+temprec(2,j)
                dem(3,iloc) = dem(3,iloc)+temprec(3,j)
              end if
            end do
          end if
        end do
c
c
c     first do rec rec communications
c
c     MPI : begin reception in buffer
c
        do i = 1, nrec_send
          tag = nprocloc*rankloc + prec_send(i) + 1
          call MPI_IRECV(buffer2(1,1,i),3*nlocrec,
     $     MPI_REAL8,prec_send(i),tag,commloc,reqrec(i),ierr)
        end do
c
c     MPI : move in buffer
c
        buffers2(:,(nlocrec+1):nlocrec2)=deprec(:,(nlocrec+1):nlocrec2)
c
        do i = 1, nrec_recep
          tag = nprocloc*prec_recep(i) + rankloc + 1
          call MPI_ISEND(buffers2(1,bufbegrec(prec_recep(i)+1)),
     $     3*domlenrec(prec_recep(i)+1),
     $     MPI_REAL8,prec_recep(i),tag,commloc,
     $     reqsend(i),ierr)
        end do
c
        do i = 1, nrec_send
          call MPI_WAIT(reqrec(i),status,ierr)
        end do
        do i = 1, nrec_recep
          call MPI_WAIT(reqsend(i),status,ierr)
        end do
c
c     MPI : move in global arrays
c
        temprec = deprec

        do i = 1, nrec_send
          do j = 1, nlocrec
            temprec(1,j) = temprec(1,j) + buffer2(1,j,i)
            temprec(2,j) = temprec(2,j) + buffer2(2,j,i)
            temprec(3,j) = temprec(3,j) + buffer2(3,j,i)
          end do
        end do
c
c     MPI : begin reception in buffer
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            tag = nproc*rank + precdir_send1(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_REAL8,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep1
          if (precdir_recep1(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep1(iproc)+1)
              jloc = bufbeg(precdir_recep1(iproc)+1)+i-1
              jglob = glob(jloc)
              if (repartrec(jglob).eq.rankloc) then
                jlocrec = locrec(jglob)
                buffers(1,jloc) = temprec(1,jlocrec)
                buffers(2,jloc) = temprec(2,jlocrec)
                buffers(3,jloc) = temprec(3,jlocrec)
              else
                buffers(1,jloc) = 0d0
                buffers(2,jloc) = 0d0
                buffers(3,jloc) = 0d0
              end if
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            tag = nproc*precdir_recep1(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbeg(precdir_recep1(i)+1)),
     $       3*domlen(precdir_recep1(i)+1),
     $       MPI_REAL8,precdir_recep1(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
          end if
        end do
c
c       MPI : move in global arrays
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            do j = 1, nloc
              iglob = glob(j)
              dep(1,j) = dep(1,j) + buffer(1,j,i)
              dep(2,j) = dep(2,j) + buffer(2,j,i)
              dep(3,j) = dep(3,j) + buffer(3,j,i)
            end do
          else
            do j = 1, nlocrec
              iglob = globrec(j)
              iloc = loc(iglob)
              if (repart(iglob).eq.rank) then
                dep(1,iloc) = dep(1,iloc)+temprec(1,j)
                dep(2,iloc) = dep(2,iloc)+temprec(2,j)
                dep(3,iloc) = dep(3,iloc)+temprec(3,j)
              end if
            end do
          end if
        end do
c
c     first do rec rec communications
c
c     MPI : begin reception in buffer
c
        do i = 1, nrec_send
          tag = nprocloc*rankloc + prec_send(i) + 1
          call MPI_IRECV(buffer2(1,1,i),3*nlocrec,
     $     MPI_REAL8,prec_send(i),tag,commloc,reqrec(i),ierr)
        end do
c
c     MPI : move in buffer
c
        buffers2(:,(nlocrec+1):nlocrec2)=decrec(:,(nlocrec+1):nlocrec2)
c
        do i = 1, nrec_recep
          tag = nprocloc*prec_recep(i) + rankloc + 1
          call MPI_ISEND(buffers2(1,bufbegrec(prec_recep(i)+1)),
     $     3*domlenrec(prec_recep(i)+1),
     $     MPI_REAL8,prec_recep(i),tag,commloc,
     $     reqsend(i),ierr)
        end do
c
        do i = 1, nrec_send
          call MPI_WAIT(reqrec(i),status,ierr)
        end do
        do i = 1, nrec_recep
          call MPI_WAIT(reqsend(i),status,ierr)
        end do
c
c     MPI : move in global arrays
c
        do i = 1, nrec_send
          do j = 1, nlocrec
            temprec(1,j) = temprec(1,j) + buffer2(1,j,i)
            temprec(2,j) = temprec(2,j) + buffer2(2,j,i)
            temprec(3,j) = temprec(3,j) + buffer2(3,j,i)
          end do
        end do
        temprec = decrec
c
c     MPI : begin reception in buffer
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            tag = nproc*rank + precdir_send1(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_REAL8,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep1
          if (precdir_recep1(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep1(iproc)+1)
              jloc = bufbeg(precdir_recep1(iproc)+1)+i-1
              jglob = glob(jloc)
              if (repartrec(jglob).eq.rankloc) then
                jlocrec = locrec(jglob)
                buffers(1,jloc) = temprec(1,jlocrec)
                buffers(2,jloc) = temprec(2,jlocrec)
                buffers(3,jloc) = temprec(3,jlocrec)
              else
                buffers(1,jloc) = 0d0
                buffers(2,jloc) = 0d0
                buffers(3,jloc) = 0d0
              end if
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            tag = nproc*precdir_recep1(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbeg(precdir_recep1(i)+1)),
     $       3*domlen(precdir_recep1(i)+1),
     $       MPI_REAL8,precdir_recep1(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
          end if
        end do
c
c       MPI : move in global arrays
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            do j = 1, nloc
              iglob = glob(j)
              dec(1,j) = dec(1,j) + buffer(1,j,i)
              dec(2,j) = dec(2,j) + buffer(2,j,i)
              dec(3,j) = dec(3,j) + buffer(3,j,i)
            end do
          else
            do j = 1, nlocrec
              iglob = globrec(j)
              iloc = loc(iglob)
              if (repart(iglob).eq.rank) then
                dec(1,iloc) = dec(1,iloc)+temprec(1,j)
                dec(2,iloc) = dec(2,iloc)+temprec(2,j)
                dec(3,iloc) = dec(3,iloc)+temprec(3,j)
              end if
            end do
          end if
        end do
c
c     first do rec rec communications
c
c     MPI : begin reception in buffer
c
        do i = 1, nrec_send
          tag = nprocloc*rankloc + prec_send(i) + 1
          call MPI_IRECV(buffer2(1,1,i),3*nlocrec,
     $     MPI_REAL8,prec_send(i),tag,commloc,reqrec(i),ierr)
        end do
c
c     MPI : move in buffer
c
        buffers2(:,(nlocrec+1):nlocrec2)=dedsprec(:,(nlocrec+1):
     $     nlocrec2)
c
        do i = 1, nrec_recep
          tag = nprocloc*prec_recep(i) + rankloc + 1
          call MPI_ISEND(buffers2(1,bufbegrec(prec_recep(i)+1)),
     $     3*domlenrec(prec_recep(i)+1),
     $     MPI_REAL8,prec_recep(i),tag,commloc,
     $     reqsend(i),ierr)
        end do
c
        do i = 1, nrec_send
          call MPI_WAIT(reqrec(i),status,ierr)
        end do
        do i = 1, nrec_recep
          call MPI_WAIT(reqsend(i),status,ierr)
        end do
c
c     MPI : move in global arrays
c
        temprec = dedsprec

        do i = 1, nrec_send
          do j = 1, nlocrec
            temprec(1,j) = temprec(1,j) + buffer2(1,j,i)
            temprec(2,j) = temprec(2,j) + buffer2(2,j,i)
            temprec(3,j) = temprec(3,j) + buffer2(3,j,i)
          end do
        end do
c
c     MPI : begin reception in buffer
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            tag = nproc*rank + precdir_send1(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_REAL8,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
          end if
        end do
c
c       MPI : move in buffer
c
        do iproc = 1, nrecdir_recep1
          if (precdir_recep1(iproc).ne.rank) then
            do i = 1, domlen(precdir_recep1(iproc)+1)
              jloc = bufbeg(precdir_recep1(iproc)+1)+i-1
              jglob = glob(jloc)
              if (repartrec(jglob).eq.rankloc) then
                jlocrec = locrec(jglob)
                buffers(1,jloc) = temprec(1,jlocrec)
                buffers(2,jloc) = temprec(2,jlocrec)
                buffers(3,jloc) = temprec(3,jlocrec)
              else
                buffers(1,jloc) = 0d0
                buffers(2,jloc) = 0d0
                buffers(3,jloc) = 0d0
              end if
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            tag = nproc*precdir_recep1(i) + rank + 1
            call MPI_ISEND(buffers(1,bufbeg(precdir_recep1(i)+1)),
     $       3*domlen(precdir_recep1(i)+1),
     $       MPI_REAL8,precdir_recep1(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
          end if
        end do
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          if (precdir_recep1(i).ne.rank) then
            call MPI_WAIT(reqsend(i),status,ierr)
          end if
        end do
c
c       MPI : move in global arrays
c
        do i = 1, nrecdir_send1
          if (precdir_send1(i).ne.rank) then
            do j = 1, nloc
              iglob = glob(j)
              dedsp(1,j) = dedsp(1,j) + buffer(1,j,i)
              dedsp(2,j) = dedsp(2,j) + buffer(2,j,i)
              dedsp(3,j) = dedsp(3,j) + buffer(3,j,i)
            end do
          else
            do j = 1, nlocrec
              iglob = globrec(j)
              iloc = loc(iglob)
              if (repart(iglob).eq.rank) then
                dedsp(1,j) = dedsp(1,j)+temprec(1,j)
                dedsp(2,j) = dedsp(2,j)+temprec(2,j)
                dedsp(3,j) = dedsp(3,j)+temprec(3,j)
              end if
            end do
          end if
        end do
      end  if

      deallocate (temprec)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      deallocate (buffer2)
      deallocate (buffers2)
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
      buffers(:,(nloc+1):nbloc) = der(:,(nloc+1):nbloc)
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
          der(1,j) = der(1,j) + buffer(1,j,i)
          der(2,j) = der(2,j) + buffer(2,j,i)
          der(3,j) = der(3,j) + buffer(3,j,i)
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
      buffers(:,(nloc+1):nbloc) = dedsp(:,(nloc+1):nbloc)
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
          dedsp(1,j) = dedsp(1,j) + buffer(1,j,i)
          dedsp(2,j) = dedsp(2,j) + buffer(2,j,i)
          dedsp(3,j) = dedsp(3,j) + buffer(3,j,i)
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
      buffers(:,(nloc+1):nbloc) = dect(:,(nloc+1):nbloc)
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
          dect(1,j) = dect(1,j) + buffer(1,j,i)
          dect(2,j) = dect(2,j) + buffer(2,j,i)
          dect(3,j) = dect(3,j) + buffer(3,j,i)
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
       call MPI_ALLREDUCE(MPI_IN_PLACE,er,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,edsp,1,MPI_REAL8,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ect,1,MPI_REAL8,MPI_SUM,
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
c      call info_energy(0)
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
c     get the domlen values of the coupled real space processes (dir-recip neighbors but not present in the previous list)
c
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc).ne.rank) then
        tag = nproc*rank + precdir_recep2(iproc) + 1
        call MPI_IRECV(domlen(precdir_recep2(iproc)+1),1,MPI_INT,
     $   precdir_recep2(iproc),tag,COMM_TINKER,reqrec(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_send2
        if (precdir_send2(iproc).ne.rank) then
        tag = nproc*precdir_send2(iproc) + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,precdir_send2(iproc),
     $   tag,COMM_TINKER,reqsend(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc).ne.rank) then
        call MPI_WAIT(reqrec(iproc),status,ierr)
        end if
      end do
      do iproc = 1, nrecdir_send2
        if (precdir_send2(iproc).ne.rank) then
        call MPI_WAIT(reqsend(iproc),status,ierr)
        end if
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
      do iproc = 1, nrecdir_recep1
        if (precdir_recep1(iproc).eq.rank) goto 10
        do iproc1 = 1, nbig_recep
          if (pbig_recep(iproc1).eq.precdir_recep1(iproc)) goto 10
        end do
        if (domlen(precdir_recep1(iproc)+1).ne.0) then
          bufbeg(precdir_recep1(iproc)+1) = count + 1
        else
          bufbeg(precdir_recep1(iproc)+1) = 1
        end if
        count = count + domlen(precdir_recep1(iproc)+1)
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
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc).ne.rank) then
          tag = nproc*rank + precdir_recep2(iproc) + 1
          call MPI_IRECV(glob(bufbeg(precdir_recep2(iproc)+1)),
     $     domlen(precdir_recep2(iproc)+1),MPI_INT,
     $     precdir_recep2(iproc),tag,COMM_TINKER,reqrec(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_send2
        if (precdir_send2(iproc).ne.rank) then
          tag = nproc*precdir_send2(iproc) + rank + 1
          call MPI_ISEND(glob,domlen(rank+1),MPI_INT,
     $     precdir_send2(iproc),tag,COMM_TINKER,reqsend(iproc),ierr)
        end if
      end do
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc).ne.rank) then
          call MPI_WAIT(reqrec(iproc),status,ierr)
        end if
      end do
      do iproc = 1, nrecdir_send2
        if (precdir_send2(iproc).ne.rank) then
          call MPI_WAIT(reqsend(iproc),status,ierr)
        end if
      end do
c
      do iproc = 1, nrecdir_recep2
        if (precdir_recep2(iproc).ne.rank) then
          do i = 1, domlen(precdir_recep2(iproc)+1)
            loc(glob(bufbeg(precdir_recep2(iproc)+1)+i-1)) =
     $        bufbeg(precdir_recep2(iproc)+1)+i-1
            repart(glob(bufbeg(precdir_recep2(iproc)+1)+i-1)) = 
     $       precdir_recep2(iproc)
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
      use mpi
      use pme
      implicit none
      integer, allocatable :: counttemp(:),ind1temp(:)
      integer count,rankloc,nprocloc,commloc
      integer i,iproc
      integer iloc,iglob,ierr,tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqsend(:),reqrec(:)
c
      allocate (counttemp(nproc))
      allocate (ind1temp(n))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
c
      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc = comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = COMM_TINKER
      end if
c
c      ind1temp = 0 
      counttemp = 0
      domlenrec = 0
c
      if ((.not.(use_pmecore)).or.((use_pmecore).and.(rank.gt.ndir-1)))
     $ then
        do iproc = 1, nrecdir_recep1
          do i = 1, domlen(precdir_recep1(iproc)+1)
            iloc = bufbeg(precdir_recep1(iproc)+1)+i-1
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
c     get the domlen values of the coupled reciprocal space processes
c
      do iproc = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(iproc) + 1
        call MPI_IRECV(domlenrec(prec_recep(iproc)+1),1,MPI_INT,
     $   prec_recep(iproc),tag,commloc,reqrec(iproc),ierr)
      end do
      do iproc = 1, nrec_send
        tag = nprocloc*prec_send(iproc) + rankloc + 1
        call MPI_ISEND(nlocrec,1,MPI_INT,prec_send(iproc),tag,
     $   commloc,reqsend(iproc),ierr)
      end do
      do iproc = 1, nrec_recep
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nrec_send
        call MPI_WAIT(reqsend(iproc),status,ierr)
      end do
c
      bufbegrec(rankloc+1) = 1
      count = nlocrec

      nlocrec2 = count
      count = nlocrec2

      do iproc = 1, nrec_recep
        if (domlenrec(prec_recep(iproc)+1).ne.0) then
          bufbegrec(prec_recep(iproc)+1) = count + 1
        else
          bufbegrec(prec_recep(iproc)+1) = 1
        end if
        count = count + domlenrec(prec_recep(iproc)+1)
 10     continue
      end do
      nlocrec2 = count
c
c     send and receive the indexes of the neighboring reciprocal processes
c

      do iproc = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(iproc) + 1
        call MPI_IRECV(globrec(bufbegrec(prec_recep(iproc)+1)),
     $   domlenrec(prec_recep(iproc)+1),MPI_INT,prec_recep(iproc),tag,
     $   commloc,reqrec(iproc),ierr)
      end do
      do iproc = 1, nrec_send
        tag = nprocloc*prec_send(iproc) + rankloc + 1
        call MPI_ISEND(globrec,domlenrec(rankloc+1),MPI_INT,
     $   prec_send(iproc),
     $   tag,commloc,reqsend(iproc),ierr)
      end do
      do iproc = 1, nrec_recep
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nrec_send
        call MPI_WAIT(reqsend(iproc),status,ierr)
      end do
c
      do iproc = 1, nrec_recep
        do i = 1, domlenrec(prec_recep(iproc)+1)
          locrec(globrec(bufbegrec(prec_recep(iproc)+1)+i-1)) =
     $      bufbegrec(prec_recep(iproc)+1)+i-1
          repartrec(globrec(bufbegrec(prec_recep(iproc)+1)+i-1)) = 
     $     prec_recep(iproc)
        end do
      end do

      deallocate (ind1temp)
      deallocate (counttemp)
      deallocate (reqsend)
      deallocate (reqrec)
c
c     also allocate spline arrays
c
      if (allocated(thetai1))  deallocate (thetai1)
      if (allocated(thetai2))  deallocate (thetai2)
      if (allocated(thetai3))  deallocate (thetai3)
      allocate (thetai1(4,bsorder,nlocrec))
      allocate (thetai2(4,bsorder,nlocrec))
      allocate (thetai3(4,bsorder,nlocrec))

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
c
c     subroutine commchgflx : communicate increment to modify charges to neighboring processes  
c
      subroutine commchgflx(pdelta)
      use atmlst
      use atoms
      use charge
      use domdec
      use mpole
      use mpi
      use potent
      implicit none
      integer i,tag,ierr,iproc
      integer iloc,iglob,j
      integer, allocatable :: reqrec(:),reqsend(:)
      real*8, allocatable :: buffer(:,:),buffer2(:)
      real*8, allocatable :: buffers(:)
      real*8 pdelta(*)
      integer status(MPI_STATUS_SIZE)
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(max(1,nloc),nbig_send))
      allocate (buffers(max(1,nbloc)))
c
      do i = 1, nbig_recep
        do j = 1, domlen(pbig_recep(i)+1)
          iloc = bufbeg(pbig_recep(i)+1)+j-1
          iglob = glob(iloc)
          buffers(iloc) = pdelta(iglob)
        end do
      end do
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
        call MPI_IRECV(buffer(1,i),nloc,MPI_REAL8,pbig_send(i),
     $   tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
        call MPI_ISEND(buffers(bufbeg(pbig_recep(i)+1)),
     $   domlen(pbig_recep(i)+1),
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
      do iproc = 1, nbig_send
        do i = 1, nloc
          iglob = glob(i)
          pdelta(iglob) = pdelta(iglob) + buffer(i,iproc)
        end do
      end do
c
c     then send accumulated values to neighbors
c
      deallocate (buffer)
      allocate (buffer2(max(1,nbloc)))
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_recep
        tag = nproc*rank + pbig_recep(i) + 1
        call MPI_IRECV(buffer2(bufbeg(pbig_recep(i)+1)),
     $   domlen(pbig_recep(i)+1),
     $   MPI_REAL8,pbig_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      do i = 1, nloc
        iglob = glob(i)
        buffers(i) = pdelta(iglob)
      end do
c
c     MPI : begin sending
c
      do i = 1, nbig_send
        tag = nproc*pbig_send(i) + rank + 1
        call MPI_ISEND(buffers,nloc,
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
          pdelta(iglob) = buffer2(iloc)
        end do
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer2)
      deallocate (buffers)
      return
      end
c
c     subroutine commpot : communicate potential to neighboring processes (chgflux)
c
c     rule = 0: communicate to pbig_* list of processes
c     rule = 1: communicate to pneig_* list of processes
c
c
      subroutine commpot(pot,rule)
      use atmlst
      use atoms
      use charge
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer i,tag,ierr,iproc
      integer iloc
      integer rule
      integer, allocatable :: reqrec(:),reqsend(:)
      integer :: n_recep, n_send
      integer, allocatable :: p_send(:),p_recep(:)
      real*8, allocatable :: buffer(:)
      real*8, allocatable :: buffers(:)
      real*8 pot(*)
      integer status(MPI_STATUS_SIZE)
 1000 format(' illegal rule in commpot.')
c
      allocate (p_send(nproc))
      allocate (p_recep(nproc))
c
      if (rule.eq.0) then
        n_recep = nbig_recep
        n_send  = nbig_send
        p_send  = pbig_send
        p_recep = pbig_recep
      else if (rule.eq.1) then
        n_recep = nneig_recep
        n_send  = nneig_send
        p_send  = pneig_send
        p_recep = pneig_recep
      else
        if (rank.eq.0) write(iout,1000)
      end if

c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(max(1,nbloc)))
      allocate (buffers(max(1,nloc)))
c
c     MPI : begin reception in buffer
c
      do i = 1, n_recep
        tag = nproc*rank + p_recep(i) + 1
        call MPI_IRECV(buffer(bufbeg(p_recep(i)+1)),
     $   domlen(p_recep(i)+1),
     $   MPI_REAL8,p_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
      do i = 1, nloc
        buffers(i) = pot(i)
      end do
c
c     MPI : begin sending
c
      do i = 1, n_send
        tag = nproc*p_send(i) + rank + 1
        call MPI_ISEND(buffers,nloc,
     $   MPI_REAL8,p_send(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, n_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do iproc = 1, n_recep
        do i = 1, domlen(p_recep(iproc)+1)
          iloc = bufbeg(p_recep(iproc)+1)+i-1
          pot(iloc) = buffer(iloc)
        end do
      end do
c
      deallocate (p_send)
      deallocate (p_recep)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commpotsum : communicate and sum potential to neighboring processes (chgflux)
c
c     rule = 0: communicate to pbig_* list of processes
c     rule = 1: communicate to pneig_* list of processes
c
c
      subroutine commpotsum(pot,rule)
      use atmlst
      use atoms
      use charge
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer i,tag,ierr
      integer rule
      integer, allocatable :: reqrec(:),reqsend(:)
      integer :: n_recep, n_send
      integer, allocatable :: p_send(:),p_recep(:)
      real*8, allocatable :: buffer(:,:)
      real*8, allocatable :: buffers(:)
      real*8 pot(*)
      integer status(MPI_STATUS_SIZE)
 1000 format(' illegal rule in commpot.')
c
      allocate (p_send(nproc))
      allocate (p_recep(nproc))
c
      if (rule.eq.0) then
        n_recep = nbig_recep
        n_send  = nbig_send
        p_send  = pbig_send
        p_recep = pbig_recep
      else if (rule.eq.1) then
        n_recep = nneig_recep
        n_send  = nneig_send
        p_send  = pneig_send
        p_recep = pneig_recep
      else
        if (rank.eq.0) write(iout,1000)
      end if

c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(max(1,nloc),n_send))
      allocate (buffers(max(1,nbloc)))
c
c     MPI : move in buffer
c
      buffers((nloc+1):nbloc) = pot((nloc+1):nbloc)
c
      do i = 1, n_send
        tag = nproc*rank + p_send(i) + 1
        call MPI_IRECV(buffer(1,i),nloc,MPI_REAL8,p_send(i),
     $   tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, n_recep
        tag = nproc*p_recep(i) + rank + 1
        call MPI_ISEND(buffers(bufbeg(p_recep(i)+1)),
     $   domlen(p_recep(i)+1),
     $   MPI_REAL8,p_recep(i),tag,COMM_TINKER,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, n_send
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_recep
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, n_send
        pot(1:nloc) = pot(1:nloc) + buffer(1:nloc,i)
      end do
c
c
      deallocate (p_send)
      deallocate (p_recep)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (buffer)
      deallocate (buffers)
      return
      end
c
c     subroutine commpotrec : communicate reciprocal pot to local
c
      subroutine commpotrec(pot,potrec)
      use sizes
      use deriv
      use domdec
      use potent
      use mpi
      implicit none
      integer i,j,tag,ierr
      integer jloc,jglob,jlocrec
      integer, allocatable :: reqsend(:),reqrec(:)
      real*8, allocatable :: buffer(:,:)
      real*8, allocatable :: buffers(:)
      real*8 pot(*),potrec(*)
      integer status(MPI_STATUS_SIZE)
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (buffer(max(1,nloc),nrecdir_send1))
      buffer = 0d0
      allocate (buffers(max(1,nblocrecdir)))
      buffers = 0d0
c
c     do rec dir communications
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_send1
        tag = nproc*rank + precdir_send1(i) + 1
        call MPI_IRECV(buffer(1,i),nloc,
     $   MPI_REAL8,precdir_send1(i),tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     Move in array
c
      do i = 1, nrecdir_recep1
        do j = 1, domlen(precdir_recep1(i)+1)
          jloc = bufbeg(precdir_recep1(i)+1)+j-1
          jglob = glob(jloc)
          if (repartrec(jglob).eq.rank) then
            jlocrec = locrec(jglob)
            buffers(jloc) = potrec(jlocrec)
          end if
        end do
      end do
c
      do i = 1, nrecdir_recep1
        tag = nproc*precdir_recep1(i) + rank + 1
        call MPI_ISEND(buffers(bufbeg(precdir_recep1(i)+1)),
     $   domlen(precdir_recep1(i)+1),
     $   MPI_REAL8,precdir_recep1(i),tag,COMM_TINKER,reqsend(i),ierr)
      end do
c
      do i = 1, nrecdir_send1
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrecdir_recep1
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     MPI : move in global arrays
c
      do i = 1, nrecdir_send1
        do j = 1, nloc
          pot(j) = pot(j) + buffer(j,i)
        end do
      end do

      return
      end
c
c     subroutine commpbcwrapindex: communicate local pbcwrap index during reassign
c
      subroutine commpbcwrapindex(bufbeg3,buf3,buflen3,bufbeg4,
     $  buflen4)
      use atoms
      use domdec
      use mpi
      use potent
      implicit none
      integer tag,status(MPI_STATUS_SIZE),i,j,proc
      integer nprocloc,rankloc,commloc
      integer iglob,jglob,ierr
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, dimension (*) :: buflen3,bufbeg3,buflen4,bufbeg4
      integer, dimension (*) :: buf3
      integer, allocatable :: buffers(:,:),buffer(:,:)
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
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
c
c
c     communicate the corresponding indexes
c
      proc = pneig_send(nneig_send)
      allocate (buffers(4,bufbeg3(proc+1)+buflen3(proc+1)))
      proc = pneig_recep(nneig_recep)
      allocate (buffer(4,bufbeg4(proc+1)+buflen4(proc+1)))
c
c     Begin reception 
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buffer(1,bufbeg4(proc+1)),4*buflen4(proc+1),
     $    MPI_INT,proc,tag,COMM_TINKER,reqrec(i),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        do j = 0, buflen3(proc+1)-1
          jglob = buf3(bufbeg3(proc+1)+j)
          buffers(1,bufbeg3(proc+1)+j) = pbcwrapindex(1,jglob)
          buffers(2,bufbeg3(proc+1)+j) = pbcwrapindex(2,jglob)
          buffers(3,bufbeg3(proc+1)+j) = pbcwrapindex(3,jglob)
          buffers(4,bufbeg3(proc+1)+j) = jglob
        end do
      end do
c
c     send the pbc indexes
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buffers(1,bufbeg3(proc+1)),4*buflen3(proc+1),
     $   MPI_INT,proc,tag,COMM_TINKER,reqsend(i),ierr)
      end do
      do i = 1, nneig_send
        proc = pneig_send(i)
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        call MPI_WAIT(reqrec(i),status,ierr)
        do j = 0, buflen4(proc+1)-1
          iglob = int(buffer(4,bufbeg4(proc+1)+j))
          pbcwrapindex(1,iglob) = buffer(1,bufbeg4(proc+1)+j)
          pbcwrapindex(2,iglob) = buffer(2,bufbeg4(proc+1)+j)
          pbcwrapindex(3,iglob) = buffer(3,bufbeg4(proc+1)+j)
        end do
      end do
      return
      end
