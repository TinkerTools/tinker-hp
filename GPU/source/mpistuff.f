c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_precision.h"
c
c     subroutine commposrec: communicate positions for reciprocal space
c
      subroutine commposrec
      use atomsMirror
      use domdec
      use inform   ,only: deb_Path
      use mpi
      use timestat ,only: timer_enter,timer_exit,timer_commstep
     &             ,quiet_timers
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2
      implicit none
      integer i,iproc
      integer iglob,iloc
      integer  ,allocatable:: reqsend(:),reqrec(:)
      real(r_p),pointer    :: buffer(:,:),buffers(:,:)
      integer status(MPI_STATUS_SIZE)
      integer tag,ierr,ibufbeg
c
      if (nproc.eq.1) return
      call timer_enter( timer_commstep )
      if (deb_Path) write(*,*) '   >> commposrec'
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      call prmem_requestm(buffMpi_p1,3*nblocrecdir,async=.false.)
      call prmem_requestm(buffMpi_p2,3*nloc,async=.false.)
      buffer (1:3,1:nblocrecdir)=>buffMpi_p1(1:3*nblocrecdir)
      buffers(1:3,1:nloc)       =>buffMpi_p2(1:3*nloc)
c
c     MPI : begin reception in buffer
c
!$acc data present(buffer,buffers)
!$acc&     present(x,y,z,glob)

!$acc host_data use_device(buffer,buffers)
      do i = 1, nrecdir_recep
        if (precdir_recep(i).ne.rank) then
          tag = nproc*rank + precdir_recep(i) + 1
          call MPI_IRECV(buffer(1,bufbeg(precdir_recep(i)+1)),
     $         3*domlen(precdir_recep(i)+1),MPI_RPREC,
     $         precdir_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
c
c     MPI : move in buffer
c
!$acc parallel loop async
      do i = 1, nloc
         iglob        = glob(i)
         buffers(1,i) = x(iglob)
         buffers(2,i) = y(iglob)
         buffers(3,i) = z(iglob)
      end do
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, nrecdir_send
        if (precdir_send(i).ne.rank) then
          tag = nproc*precdir_send(i) + rank + 1
          call MPI_ISEND(buffers,3*nloc,MPI_RPREC,precdir_send(i),
     &         tag,COMM_TINKER,reqsend(i),ierr)
        end if
      end do
!$acc end host_data
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
        ibufbeg = bufbeg(precdir_recep(iproc)+1)
        if (precdir_recep(iproc).ne.rank) then
!$acc parallel loop async
         do i = 1, domlen(precdir_recep(iproc)+1)
           iloc     = ibufbeg+i-1
           iglob    = glob(iloc)
           x(iglob) = buffer(1,iloc)
           y(iglob) = buffer(2,iloc)
           z(iglob) = buffer(3,iloc)
         end do
        end if
      end do
c
!$acc end data
!$acc wait
      nullify(buffer)
      nullify(buffers)
      deallocate (reqsend)
      deallocate (reqrec)

      if (deb_Path) write(*,*) '   << commposrec'
      call timer_exit( timer_commstep,quiet_timers )
      end
c
c     subroutine commpos : communicate positions to neighboring processes
c
      subroutine commpos
      use atomsMirror
      use domdec
      use inform   ,only: deb_Path
      use freeze
      use mpi
      use potent
      use timestat ,only: timer_enter,timer_exit,timer_commstep
     &             ,quiet_timers
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2
      implicit none
      integer i,tag,ierr,iproc
      integer iloc,iglob
      real(r_p),pointer:: buffer(:,:),buffers(:,:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer status(MPI_STATUS_SIZE)
c
      if (nproc.eq.1.or.(use_pmecore).and.(rank.gt.ndir-1)) return
      call timer_enter( timer_commstep )
      if (deb_Path) write(*,*) '   >> commpos'
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      call prmem_requestm(buffMpi_p1,3*nbloc,async=.false.)
      call prmem_requestm(buffMpi_p2,3*nloc,async=.false.)
      buffer (1:3,1:nbloc)=>buffMpi_p1(1:3*nbloc)
      buffers(1:3,1:nloc) =>buffMpi_p2(1:3*nloc)
!$acc enter data attach(buffer,buffers)
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
     $   MPI_RPREC,pbig_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_RPREC,pbig_send(i),tag,COMM_TINKER,
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
!$acc update device(x(:),y(:),z(:)) async
c
!$acc exit data detach(buffer,buffers) async
      nullify(buffer)
      nullify(buffers)
      deallocate (reqsend)
      deallocate (reqrec)
!$acc wait

      if (deb_Path) write(*,*) '   << commpos'
      call timer_exit( timer_commstep,quiet_timers )
      end
c
c     subroutine commposshort : communicate positions to short range neighboring processes
c
      subroutine commposshort
      use atomsMirror
      use domdec
      use freeze
      use inform   ,only: deb_Path
      use mpi
      use potent
      use timestat ,only: timer_enter,timer_exit,timer_commstep
     &             ,quiet_timers
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2
      implicit none
      integer i,tag,ierr,iproc
      integer iloc,iglob,ibufbeg,idomlen
      real(r_p),pointer:: buffer(:,:),buffers(:,:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer status(MPI_STATUS_SIZE)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return

      call timer_enter( timer_commstep )
      if (deb_Path) write(*,*) '   >> commposshort'
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      call prmem_requestm(buffMpi_p1,3*nbloc,async=.false.)
      call prmem_requestm(buffMpi_p2,3*nloc,async=.false.)
      buffer (1:3,1:nbloc)=>buffMpi_p1(1:3*nbloc)
      buffers(1:3,1:nloc) =>buffMpi_p2(1:3*nloc)
!$acc enter data attach(buffer,buffers)
      !FIXME Offload this part on device
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
     $   MPI_RPREC,pbigshort_recep(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_RPREC,pbigshort_send(i),tag,COMM_TINKER,
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
!$acc exit data detach(buffer,buffers)
      nullify(buffer)
      nullify(buffers)
      deallocate (reqsend)
      deallocate (reqrec)

      if (deb_Path) write(*,*) '   << commposshort'
      call timer_exit( timer_commstep,quiet_timers )
      end
c
c     subroutine reassign: deal with particles changing of domain between two time steps
c
      subroutine reassign
      use atomsMirror
      use bound
      use cell
      use domdec
      use freeze
      use moldyn
      use neigh
      use potent
      use timestat ,only: timer_enter,timer_exit,timer_eneig
     &             ,quiet_timers
      use tinheader,only: ti_eps
      use mpi
      implicit none
      integer i
      integer iproc,ierr,iglob

      real(t_p) xr,yr,zr,eps1,eps2
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: count2(:)
      integer, allocatable :: buf3(:), buflen3(:), bufbeg3(:)
      integer, allocatable :: buflen4(:),bufbeg4(:)
      real(r_p), allocatable :: buffers(:,:),buffer(:,:)
      real(r_p), allocatable :: buffers1(:,:),buffer1(:,:)
      integer, allocatable :: buf3bis(:,:)
      integer, allocatable :: glob1(:)
      integer rankloc,commloc,nprocloc,proc,j,jglob
      integer nloc1
#if (defined(SINGLE)||defined(MIXED))
      parameter(eps1 = 1d1*ti_eps, eps2 = 1d3*ti_eps)
#else
      parameter(eps1 = 1d2*ti_eps, eps2 = 1d4*ti_eps)
#endif
c
      if (nproc.eq.1) return

      call timer_enter( timer_eneig )
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

      count2  = 0
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
          if ((zr.ge.zbegproc(iproc+1)).and.(zr.lt.zendproc(iproc+1))
     $   .and.(yr.ge.ybegproc(iproc+1)).and.(yr.lt.yendproc(iproc+1))
     $   .and.(xr.ge.xbegproc(iproc+1)).and.(xr.lt.xendproc(iproc+1)))
     $    then
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
      if (nprocloc.eq.1) then
         call timer_exit( timer_eneig,quiet_timers )
         return
      end if
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
        tag  = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buflen4(proc+1),1,MPI_INT,
     $   proc,tag,commloc,reqrec(iproc),ierr)
      end do
      do iproc = 1, nneig_send
        proc = pneig_send(iproc)
        tag  = nprocloc*proc + rankloc + 1
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
     $       MPI_RPREC,proc,tag,COMM_TINKER,reqrec(i),ierr)
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
     $   MPI_RPREC,proc,tag,COMM_TINKER,reqsend(i),ierr)
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
!$acc data present(x,y,z,v,a)
!$acc update device(x,y,z,v,a) async
!$acc end data
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
     $         MPI_RPREC,proc,tag,COMM_TINKER,reqrec(i),ierr)
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
     $         MPI_RPREC,proc,tag,COMM_TINKER,reqsend(i),ierr)
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
!$acc data present(xold,yold,zold)
!$acc update device(xold,yold,zold) async
!$acc end data
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
        glob1(nloc1)  = iglob
        loc(iglob)    = nloc1
        repart(iglob) = rank
 10     continue
      end do
c
c     add atoms that entered the domain
c
      proc = pneig_recep(nneig_recep)+1
      do j = 1, bufbeg4(proc)+buflen4(proc)-1
         jglob         = int(buffer(10,j))
         nloc1         = nloc1 + 1
         glob1(nloc1)  = jglob
         loc(jglob)    = nloc1
         repart(jglob) = rank
      end do
!$acc update device(loc(:),repart(:)) async
      deallocate (buffer)
      deallocate (buffers)
c
      nloc = nloc1
      domlen(rank+1) = nloc
      glob = glob1
!$acc update device(glob) async
c
 20   call orderbuffer(.false.)
c
      call timer_exit( timer_eneig,quiet_timers )
      deallocate (buf3bis)
      deallocate (glob1)
      deallocate (count2)
      deallocate (buf3)
      deallocate (buflen3)
      deallocate (bufbeg3)
      deallocate (reqsend)
      deallocate (reqrec)
!$acc wait
      end

c
c     Assign all atoms to their real space domain
c
      subroutine AllDirAssign
      use atoms
      use cell
      use domdec
      use inform
      use potent
      use timestat ,only: timer_eneig,quiet_timers,timer_enter
     &             ,timer_exit
      implicit none
      integer i,iproc
      integer nprocloc,rankloc
      real(t_p) xr,yr,zr
!$acc routine(image_acc)

      if (nproc.eq.1.or.use_pmecore.and.ndir.eq.1) return
      if (deb_Path) print*, 'AllDirAssign'
      call timer_enter(timer_eneig)

      if (use_pmecore) then
        nprocloc = ndir
        rankloc  = rank_bis
      else
        nprocloc = nproc
        rankloc  = rank
      end if

!$acc parallel loop async default(present)
      do i = 1,n
         glob(i)   = 0
         repart(i) = 0
         loc(i)    = 0
      end do

      domlen(:) = 0

!$acc parallel loop async default(present) copy(domlen)
      do i = 1, n
        xr = x(i)
        yr = y(i)
        zr = z(i)
        call image_acc(xr,yr,zr)
        if (abs(xr-xcell2).lt.eps_cell) xr = xr-5*eps_cell
        if (abs(yr-ycell2).lt.eps_cell) yr = yr-5*eps_cell
        if (abs(zr-zcell2).lt.eps_cell) zr = zr-5*eps_cell
        do iproc = 0, nprocloc-1
          if ((zr.ge.zbegproc(iproc+1)).and.(zr.lt.zendproc(iproc+1))
     $   .and.(yr.ge.ybegproc(iproc+1)).and.(yr.lt.yendproc(iproc+1))
     $   .and.(xr.ge.xbegproc(iproc+1)).and.(xr.lt.xendproc(iproc+1)))
     &       then
             repart(i) = iproc
!$acc atomic
             domlen(iproc+1) = domlen(iproc+1) + 1
          end if
        end do
      end do

!$acc wait
      if ((use_pmecore).and.(rank.le.ndir-1)) then
        nloc = domlen(rankloc+1)
      else if ((use_pmecore).and.(rank.gt.ndir-1)) then
        nloc = 0
      else
        nloc = domlen(rankloc+1)
      end if

      call orderbuffer_gpu(.true.)
      call timer_exit(timer_eneig,quiet_timers)
      end subroutine

c
c     Reciproqual space domain assignment of all atoms
c
      subroutine AllRecAssign
        use atoms  ,only: x,y,z,n
        use boxes  ,only: recip
        use chunks ,only: grdoff
        use domdec ,only: rank,rank_bis,ndir,nrec,nproc
     &             ,repartrec
        use fft    ,only: kstart1,kend1,jstart1,jend1
        use inform ,only: deb_Path
        use pme      ,only: pme_eps,nfft1,nfft2,nfft3,bsorder
        use potent   ,only: use_pmecore
        use tinheader,only: ti_p
        use timestat ,only: timer_eneig,quiet_timers,timer_enter
     &               ,timer_exit
        implicit none
        integer i,j,k,ifr,iproc,nprocloc,rankloc
        real(t_p) xi,yi,zi,w,fr

        if (use_pmecore) then
          nprocloc = ndir
          rankloc  = rank_bis
          if (nrec.eq.1) return
        else
          nprocloc = nproc
          rankloc  = rank
          if (nproc.eq.1) return
        end if
        if (deb_Path) print*, 'AllRecAssign'
        call timer_enter(timer_eneig)

!$acc parallel loop async default(present)
        do i = 1, n
          xi  = x(i)
          yi  = y(i)
          zi  = z(i)
          w   = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
          fr  = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
          ifr = int(fr-pme_eps)
          !if (ifr-fr.eq.0.0.and.ifr>0) ifr=ifr-1
          k   = ifr- bsorder + grdoff
          if (k .lt. 1)  k = k + nfft3
          w   = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
          fr  = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
          ifr = int(fr-pme_eps)
          !if (ifr-fr.eq.0.0.and.ifr>0) ifr=ifr-1
          j   = ifr- bsorder + grdoff
          if (j.lt.1) j = j + nfft2
!$acc loop seq
          do iproc = 0, nprocloc-1
            if (((k.ge.kstart1(iproc+1)).and.(k.le.kend1(iproc+1))).and.
     $       ((j.ge.jstart1(iproc+1)).and.(j.le.jend1(iproc+1)))) then
              repartrec(i) = iproc
            end if
          end do
        end do
        call orderbufferrec_gpu
        call timer_exit(timer_eneig,quiet_timers)
      end subroutine

      subroutine reassignrespa(ialt,nalt)
        use atomsMirror ,only : x,y,z,n
        use bound ,only : use_bounds
        use cell  ,only : xcell2, ycell2, zcell2, eps_cell
        use domdec,only : nloc,nproc,glob,comm_dir,rank,rank_bis,ndir,
     &                    xbegproc,ybegproc,zbegproc,
     &                    xendproc,yendproc,zendproc,
     &                    repart, loc, domlen, COMM_TINKER,
     &                    nneig_send, nneig_recep,
     &                    pneig_recep, pneig_send
        use inform,only : deb_Path,tindPath
        use moldyn,only : a,v
        use potent,only : use_pmecore
        use timestat ,only: timer_enter,timer_exit,timer_eneig
     &               ,quiet_timers
        use tinheader,only: ti_eps
        use mpi
        use sizes    ,only: tinkerdebug
        implicit none
        integer, intent(in) :: ialt, nalt
        integer :: nprocloc, commloc, rankloc ! local (pme) mpi info
        real(t_p) :: xr, yr, zr ! adjusted position
        integer :: local_neighbor(nloc) ! Local domain for iloc. 0 is local rank
        type t_elt
          real(r_p) :: x,y,z,vx,vy,vz,ax,ay,az
        end type
        type(t_elt),pointer :: d
        integer :: iglob_send(nloc,0:nneig_send) ! Global index to send to each proc
        integer, allocatable :: iglob_recv(:,:) ! Global index to recieve from procs
        type(t_elt),allocatable,target :: data_send(:,:), data_recv(:,:) ! Data (x,y,...) to exchange
        integer  n_data_send(0:nneig_send), n_data_recv(nneig_recep)
        integer req_iglob_send(nneig_send),req_iglob_recv(nneig_recep)
        integer  req_data_send(nneig_send), req_data_recv(nneig_recep)
        integer  mpi_status1(MPI_STATUS_SIZE)
        integer i, ierr, iglob, iloc, iproc, ineighbor
        integer nloc_save, n_data_send_capture,max_data_recv
        integer nloc_capture
        integer :: max_data_recv_save=0
#if  (defined(SINGLE)||defined(MIXED))
        real(t_p), parameter:: eps1= 10*ti_eps, eps2= 1d2*ti_eps
#else
        real(t_p), parameter:: eps1= 1d2*ti_eps, eps2= 1d4*ti_eps
#endif
!$acc routine(image_acc) seq

        if (ialt/=nalt) return
        if (use_pmecore) then
           nprocloc = ndir
           commloc  = comm_dir
           rankloc  = rank_bis
           if(rank >= ndir) then
              call timer_enter( timer_eneig )
              goto 20
           end if
        else
           nprocloc = nproc
           commloc  = COMM_TINKER
           rankloc  = rank
        end if
        if (nproc.eq.1) return
        if (use_pmecore.and.ndir.eq.1) then
           call timer_enter( timer_eneig )
           goto 20
        end if

        if (deb_Path) write(*,'(A,2I4)') '  -- reassignrespa ',ialt,nalt
        call timer_enter( timer_eneig )

!$acc wait
!$acc data copyin(pneig_recep, pneig_send) async
!$acc& present(repart,loc,glob,x,y,z,a,v)

!$acc enter data create(local_neighbor) async

        ! Compute domain owning each particle (local_neighbor(:))
!$acc kernels default(present) async
        local_neighbor(:) = 0
!$acc end kernels
!$acc parallel loop default(present) async
        do i = 1, nloc
          iglob = glob(i)
          xr = x(iglob)
          yr = y(iglob)
          zr = z(iglob)
          if (use_bounds) call image_acc(xr,yr,zr)
          if (abs(xr-xcell2).lt.eps_cell) xr = xr-eps2
          if (abs(yr-ycell2).lt.eps_cell) yr = yr-eps2
          if (abs(zr-zcell2).lt.eps_cell) zr = zr-eps2
!$acc loop seq
          do ineighbor = 1, nneig_send
            iproc = pneig_send(ineighbor)+1
            if( (iproc /= rank+1)
     &    .and. (xr>=xbegproc(iproc)) .and. (xr<xendproc(iproc))
     &    .and. (yr>=ybegproc(iproc)) .and. (yr<yendproc(iproc))
     &    .and. (zr>=zbegproc(iproc)) .and. (zr<zendproc(iproc)) ) then
              local_neighbor(i) = ineighbor
              exit
            endif
          end do
        end do
        ! Filter global indexes to send in iglob_send
!$acc enter data create(n_data_send, iglob_send) async
!$acc kernels default(present) async
        n_data_send(:) = 0
!$acc end kernels
!$acc parallel loop default(present) async
        do i = 1, nloc
          ineighbor = local_neighbor(i)
!$acc atomic capture
          n_data_send(ineighbor) = n_data_send(ineighbor) + 1
          n_data_send_capture = n_data_send(ineighbor)
!$acc end atomic
          iglob_send(n_data_send_capture,ineighbor) = glob(i)
        end do

!$acc exit data delete(local_neighbor) async

!$acc update host(n_data_send) async
!Wait for n_data_send
!$acc wait
        ! Allocate and fill data to send
        allocate( data_send( maxval(n_data_send(1:nneig_send)),
     &                       nneig_send))
!$acc enter data create(data_send) async
!$acc parallel loop default(present) async
!$acc& present( iglob_send )
        do ineighbor = 1, nneig_send
!$acc loop vector
          do i = 1, n_data_send(ineighbor)
            iglob = iglob_send(i,ineighbor)
            d => data_send( i, ineighbor )
            d%x = x(iglob)
            d%y = y(iglob)
            d%z = z(iglob)
            d%vx = v(1,iglob)
            d%vy = v(2,iglob)
            d%vz = v(3,iglob)
            d%ax = a(1,iglob)
            d%ay = a(2,iglob)
            d%az = a(3,iglob)
          end do
        end do

        ! Initiate send data
        req_iglob_send(:) = MPI_REQUEST_NULL
        req_data_send(:) = MPI_REQUEST_NULL
!Wait for iglob_send and data_send
!$acc wait
!$acc host_data use_device(iglob_send, data_send)
        do ineighbor = 1, nneig_send
          call MPI_Isend(iglob_send(1,ineighbor),n_data_send(ineighbor)
     &            ,MPI_INT,   pneig_send(ineighbor), 0, COMM_TINKER,
     &            req_iglob_send(ineighbor), ierr)
        end do
!$acc end host_data

        ! Probe for sizes
        n_data_recv(:) = 0
        do ineighbor = 1, nneig_recep
          call MPI_Probe(pneig_recep(ineighbor), 0, COMM_TINKER,
     &                   mpi_status1, ierr)
          call MPI_Get_count( mpi_status1, MPI_INT,
     &                        n_data_recv(ineighbor), ierr )
        end do

        ! Allocate and recv data
        max_data_recv = maxval(n_data_recv)
        allocate( iglob_recv( max_data_recv, nneig_recep) )
        allocate( data_recv( max_data_recv, nneig_recep) )
!$acc enter data create(iglob_recv,data_recv) async
        req_iglob_recv(:) = MPI_REQUEST_NULL
        req_data_recv(:) = MPI_REQUEST_NULL
!Wait for iglob_recv and data_recv
!$acc wait
!$acc host_data use_device(iglob_recv,data_recv,data_send)
        do ineighbor = 1, nneig_recep
          call MPI_Irecv(iglob_recv(1,ineighbor),n_data_recv(ineighbor)
     &            ,MPI_INT, pneig_recep(ineighbor), 0, COMM_TINKER,
     &            req_iglob_recv(ineighbor), ierr)
          call MPI_Irecv(data_recv(1,ineighbor),9*n_data_recv(ineighbor)
     &            ,MPI_RPREC, pneig_recep(ineighbor), 1, COMM_TINKER,
     &            req_data_recv(ineighbor), ierr)
        end do

        ! Wait for iglob exchange before sending datas
        ! Overlaping those comms like before disrupts mpi with
        ! multi-node running
        call MPI_Waitall(nneig_recep,req_iglob_recv,MPI_STATUSES_IGNORE,
     &                   ierr)
        call MPI_Waitall(nneig_send ,req_iglob_send,MPI_STATUSES_IGNORE,
     &                   ierr)

        do ineighbor = 1, nneig_send
          call MPI_Isend(data_send(1,ineighbor),9*n_data_send(ineighbor)
     &            ,MPI_RPREC, pneig_send(ineighbor), 1, COMM_TINKER,
     &            req_data_send(ineighbor), ierr)
        end do
!$acc end host_data

        ! Wait all communications
        call MPI_Waitall(nneig_recep,req_data_recv ,MPI_STATUSES_IGNORE,
     &                   ierr)
        call MPI_Waitall(nneig_send ,req_data_send ,MPI_STATUSES_IGNORE,
     &                   ierr)
!$acc exit data delete(data_send) async

        ! Rebuild glob(:) with old local atoms

        nloc = n_data_send(0)
!$acc enter data copyin(nloc) async

!$acc parallel loop default(present) async
!$acc&         present(iglob_send)
        do i = 1, n
           if (i.lt.nloc+1)
     &        glob(i) = iglob_send(i, 0)
           loc(i)    = 0
           repart(i) = 0
        end do

!$acc exit data delete(iglob_send,n_data_send) async


        ! Add new atoms to glob(:) and add data
!$acc parallel loop vector_length(512) async
!$acc& default(present) present(nloc) copyin(n_data_recv)
        do ineighbor = 1, nneig_recep
!$acc loop vector
          do i = 1, n_data_recv(ineighbor)
            iglob = iglob_recv(i,ineighbor)
!$acc atomic capture
            nloc = nloc+1
            nloc_capture = nloc
!$acc end atomic
            glob(nloc_capture) = iglob
            d => data_recv(i,ineighbor)
            x(iglob)   = d%x
            y(iglob)   = d%y
            z(iglob)   = d%z
            v(1,iglob) = d%vx
            v(2,iglob) = d%vy
            v(3,iglob) = d%vz
            a(1,iglob) = d%ax
            a(2,iglob) = d%ay
            a(3,iglob) = d%az
          end do
        end do
!$acc exit data copyout(nloc) async
!Wait for nloc (is needed with pgi 19.4 because nloc is the loop bound in next kernel)
!$acc wait

        nloc_save = nloc
        if (btest(tinkerdebug,tindPath).and.
     &     max_data_recv.gt.max_data_recv_save) then
           max_data_recv_save = max_data_recv
 16     format('iRank',I4,'; ',I6,' atoms(max) reassigned; nloc=',I10)
           print 16,rank,max_data_recv_save,nloc
        end if
 17     format('reassign nloc ',I9,' max_data_recv',I9)
        if (deb_Path) write(*,17) nloc,max_data_recv
        if (tinkerdebug.gt.0.and.nproc.gt.1) then
           call AtomDebRepart
        end if

!$acc parallel loop default(present)
        do i = 1, nloc
          iglob         = glob(i)
          loc(iglob)    = i
          repart(iglob) = rank
        end do

!$acc exit data delete(iglob_recv,data_recv) async
        domlen(rank+1) = nloc

!$acc end data
!copyin(pneig_recep, pneig_send)

20      call orderbuffer_gpu(.false.)
!!$acc update host(loc,glob) async
        call timer_exit( timer_eneig,quiet_timers )
      end subroutine

c
c     subroutine reassignrespa: deal with particles changing of domain between two time steps
c
      subroutine reassignrespa_old(ialt,nalt)
      use atomsMirror
      use bound
      use cell
      use domdec
      use moldyn, nalt2=>nalt
      use neigh
      use potent
      use mpi
      use timestat ,only: timer_enter,timer_exit,timer_eneig
     &             ,quiet_timers
      use tinheader
      implicit none
      integer i,iproc,ierr,iglob,iloc
      real(t_p) xr,yr,zr,eps1,eps2
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: count2(:)
      integer, allocatable :: buf3(:), buflen3(:), bufbeg3(:)
      integer, allocatable :: buflen4(:),bufbeg4(:)
      real(r_p), allocatable :: buffers(:,:),buffer(:,:)
      integer, allocatable :: buf3bis(:,:)
      integer, allocatable :: glob1(:)
      integer rankloc,commloc,nprocloc,proc,j,jglob,jglob_max
      integer nloc1
      integer ialt,nalt
      integer buflen3_capture,nloc1_capture
#if (defined(SINGLE)||defined(MIXED))
      parameter (eps1= 10*ti_eps, eps2= 1d2*ti_eps)
#else
      parameter (eps1= 1d2*ti_eps, eps2= 1d4*ti_eps)
#endif
!$acc routine(image_acc) seq
c
      if (ialt.ne.nalt.or.nproc.eq.1) return
c
      call timer_enter( timer_eneig )
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
      count2  = 0

!$acc wait

!$acc data async
!$acc& create(buflen3,bufbeg3,buf3bis)

!$acc kernels async
      buflen3(1:nproc) = 0
!$acc end kernels
c
c     get particules that changed of domain
c
!$acc parallel loop default(present) async
      do i = 1, nloc
        iglob = glob(i)
        xr = x(iglob)
        yr = y(iglob)
        zr = z(iglob)
        if (use_bounds) call image_acc(xr,yr,zr)
        if (abs(xr-xcell2).lt.eps_cell) xr = xr-eps2
        if (abs(yr-ycell2).lt.eps_cell) yr = yr-eps2
        if (abs(zr-zcell2).lt.eps_cell) zr = zr-eps2
!$acc loop seq
        do iproc = 0, nprocloc-1
          if (iproc.eq.rank) cycle
          if ((zr.ge.zbegproc(iproc+1)).and.
     $     (zr.lt.zendproc(iproc+1)).and.(yr.ge.ybegproc(iproc+1))
     $    .and.(yr.lt.yendproc(iproc+1)).and.(xr.ge.xbegproc(iproc+1))
     $    .and.(xr.lt.xendproc(iproc+1))) then
!$acc atomic capture
            buflen3(iproc+1) = buflen3(iproc+1)+1
            buflen3_capture = buflen3(iproc+1)
!$acc end atomic
            buf3bis(buflen3_capture,iproc+1) = iglob
          end if
        end do
      end do

!$acc kernels default(present) async
      bufbeg3(1) = 1
      bufbeg3(2:) = 0
!$acc end kernels

!$acc parallel loop default(present) async
      do iproc = 1, nproc-1
        bufbeg3(iproc+1) = bufbeg3(iproc)+buflen3(iproc)
      end do

!$acc update host(bufbeg3,buflen3) async
!$acc wait
      allocate (buf3(bufbeg3(nproc)+buflen3(nproc)))
!$acc data create(buf3) async

!$acc kernels async
      buf3(:) = 0
!$acc end kernels

!$acc parallel loop default(present) async
      do iproc = 1, nproc
        buf3(bufbeg3(iproc):(bufbeg3(iproc)+buflen3(iproc)-1)) =
     $    buf3bis(1:buflen3(iproc),iproc)
      end do

      ! GLITCHY : gotos
      if (nprocloc.eq.1) goto 25
      if ((use_pmecore).and.(rank.ge.ndir)) goto 20
c
c     get size of buffers
c
      allocate (buflen4(nprocloc))
      allocate (bufbeg4(nprocloc))

!$acc data create(buflen4,bufbeg4) async

!$acc kernels async
      buflen4(:) = 0
      bufbeg4(:) = 0
!$acc end kernels

!$acc wait
!$acc host_data use_device(buflen4)
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buflen4(proc+1),1,MPI_INT,
     $   proc,tag,commloc,reqrec(iproc),ierr)
      end do
!$acc end host_data
!$acc host_data use_device(buflen3)
      do iproc = 1, nneig_send
        proc = pneig_send(iproc)
        tag = nprocloc*proc + rankloc + 1
        call MPI_ISEND(buflen3(proc+1),1,MPI_INT,
     $   proc,tag,commloc,reqsend(iproc),ierr)
      end do
!$acc end host_data
      do iproc = 1, nneig_recep
        proc = pneig_recep(iproc)
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nneig_send
        tag = nprocloc*proc + rankloc + 1
        call MPI_WAIT(reqsend(iproc),status,ierr)
      end do
!$acc serial copyin(pneig_recep) async
      bufbeg4(pneig_recep(1)+1) = 1
!$acc end serial
!$acc parallel loop default(present) async
!$acc& copyin(pneig_recep)
      do iproc = 2, nneig_recep
        bufbeg4(pneig_recep(iproc)+1) = bufbeg4(pneig_recep(iproc-1)+1)
     $    +buflen4(pneig_recep(iproc-1)+1)
      end do
c
c     communicate the corresponding indexes, positions, speed and accelerations
c
!$acc update host(bufbeg4,buflen4) async
!$acc wait
      proc = pneig_send(nneig_send)
      allocate (buffers(10,bufbeg3(proc+1)+buflen3(proc+1)))
      proc = pneig_recep(nneig_recep)
      allocate (buffer(10,bufbeg4(proc+1)+buflen4(proc+1)))
!$acc data create(buffer, buffers) async
c
c     Begin reception
c
!$acc host_data use_device(buffer)
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_IRECV(buffer(1,bufbeg4(proc+1)),10*buflen4(proc+1),
     $    MPI_RPREC,proc,tag,COMM_TINKER,reqrec(i),ierr)
      end do
!$acc end host_data
      do i = 1, nneig_send
        proc = pneig_send(i)
!$acc parallel loop default(present) async
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
!$acc wait
!$acc host_data use_device(buffers)
      do i = 1, nneig_send
         proc = pneig_send(i)
         tag = nprocloc*proc + rankloc + 1
         call MPI_ISEND(buffers(1,bufbeg3(proc+1)),
     $        10*buflen3(proc+1),MPI_RPREC,proc,tag,
     $        COMM_TINKER,reqsend(i),ierr)
      end do

!$acc end host_data
      do i = 1, nneig_send
        proc = pneig_send(i)
        tag = nprocloc*proc + rankloc + 1
        call MPI_WAIT(reqsend(i),status,ierr)
      end do

      do i = 1, nneig_recep
        proc = pneig_recep(i)
        tag = nprocloc*rankloc + proc + 1
        call MPI_WAIT(reqrec(i),status,ierr)
!$acc parallel loop default(present) async
        do j = 0, buflen4(proc+1)-1
          iglob = int(buffer(10,bufbeg4(proc+1)+j))
          x(  iglob) = buffer(1,bufbeg4(proc+1)+j)
          y(  iglob) = buffer(2,bufbeg4(proc+1)+j)
          z(  iglob) = buffer(3,bufbeg4(proc+1)+j)
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

c
c     remove atoms that left the domain
c

!$acc data create(nloc1,glob1) async
!present(loc, repart)

!$acc kernels async
      nloc1 = 0
      glob1(:) = 0
      repart(:) = 0
!$acc end kernels

      jglob_max = bufbeg3(nproc)+buflen3(nproc)-1
!$acc update device(loc) async

!$acc parallel loop default(present) async
      do i = 1, nloc
        iglob = glob(i)
!$acc loop seq
        do j = 1, jglob_max
          jglob = buf3(j)
          if (iglob.eq.jglob) goto 10
        end do
!$acc atomic capture
        nloc1 = nloc1 + 1
        nloc1_capture = nloc1
!$acc end atomic
        glob1(nloc1_capture) = iglob
        loc(iglob) = nloc1_capture
        repart(iglob) = rank
 10     continue
      end do
c
c     add atoms that entered the domain
c
      proc = pneig_recep(nneig_recep)+1
      jglob_max = bufbeg4(proc)+buflen4(proc)-1
!$acc parallel loop default(present) async
      do j = 1, jglob_max
        jglob = int(buffer(10,j))
!$acc atomic capture
        nloc1 = nloc1 + 1
        nloc1_capture = nloc1
!$acc end atomic
        glob1(nloc1_capture) = jglob
        loc(jglob) = nloc1_capture
        repart(jglob) = rank
      end do

!$acc serial present(nloc1) copyout(nloc) async
      nloc = nloc1
!$acc end serial
!$acc wait
      domlen(rank+1) = nloc
!$acc kernels async
      glob(:) = glob1(:)
!$acc end kernels

!$acc end data
!create(nloc1,glob1)

!$acc end data
!create(buffer, buffers)

!$acc end data
!create(buflen4,bufbeg4)

!$acc update host(loc,repart,glob) async
!$acc wait
c
 20   call orderbuffer(.false.)
 25   continue
!$acc end data
!create(buf3)

!$acc end data
!create(buflen3,bufbeg3,buf3bis)
      deallocate (buf3bis)
      deallocate (glob1)
      deallocate (count2)
      deallocate (buf3)
      deallocate (buflen3)
      deallocate (bufbeg3)
      deallocate (reqsend)
      deallocate (reqrec)
!$acc update host(x,y,z,a,v) async
!$acc wait
      call timer_exit( timer_eneig,quiet_timers )
      end
c
c     subroutine commposrespa : deal with communications of positions between two time steps
c     with respa integrator
c
      subroutine commposrespa(fast)
      use atomsMirror
      use domdec
      use mpi
      use inform   ,only: deb_Path
      use potent
      use timestat ,only: timer_enter,timer_exit,timer_commstep
     &             ,quiet_timers
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2
      implicit none
      integer i,tag,ierr,iproc
      integer iloc,iglob
      real(r_p),pointer:: buffer(:,:),buffers(:,:)
      integer,dimension(nproc):: reqsend,reqrec
      integer,dimension(nproc):: ptemp_recep,ptemp_send
      integer status(MPI_STATUS_SIZE)
      integer ntemp_recep,ntemp_send,ibufbeg
      logical fast
c
      if (ndir.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return

      if (deb_Path) write(*,*) '   >> commposrespa',fast
      call timer_enter( timer_commstep )
c
      if (fast) then
        ptemp_recep = pneig_recep
        ptemp_send  = pneig_send
        ntemp_recep = nneig_recep
        ntemp_send  = nneig_send
      else
        ptemp_recep = pbig_recep
        ptemp_send  = pbig_send
        ntemp_recep = nbig_recep
        ntemp_send  = nbig_send
      end if
c
c     Reallocate some buffers if necessary
c
      call prmem_requestm(buffMpi_p1,3*nbloc,async=.false.)
      call prmem_requestm(buffMpi_p2,3*nloc,async=.false.)
      buffer (1:3,1:nbloc)=>buffMpi_p1(1:3*nbloc)
      buffers(1:3,1:nloc) =>buffMpi_p2(1:3*nloc)

!$acc data present(buffer,buffers)
!$acc&     present(glob,x,y,z)
c
c     communicate positions
c
c
c     MPI : begin reception in buffer
c
!$acc host_data use_device(buffer,buffers)
      do i = 1, ntemp_recep
        tag = nproc*rank + ptemp_recep(i) + 1
        call MPI_IRECV(buffer(1,bufbeg(ptemp_recep(i)+1)),
     $       3*domlen(ptemp_recep(i)+1),MPI_RPREC,ptemp_recep(i),
     $       tag,COMM_TINKER,reqrec(i),ierr)
      end do
c
c     MPI : move in buffer
c
!$acc parallel loop async
      do i = 1, nloc
        iglob        = glob(i)
        buffers(1,i) = x(iglob)
        buffers(2,i) = y(iglob)
        buffers(3,i) = z(iglob)
      end do
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, ntemp_send
        tag = nproc*ptemp_send(i) + rank + 1
        call MPI_ISEND(buffers,3*nloc,
     $       MPI_RPREC,ptemp_send(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
      end do
!$acc end host_data
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
        ibufbeg = bufbeg(ptemp_recep(iproc)+1)
!$acc parallel loop async
        do i = 1, domlen(ptemp_recep(iproc)+1)
          iloc     = ibufbeg+i-1
          iglob    = glob(iloc)
          x(iglob) = buffer(1,iloc)
          y(iglob) = buffer(2,iloc)
          z(iglob) = buffer(3,iloc)
        end do
      end do
c
!$acc end data
!$acc wait
      nullify(buffer)
      nullify(buffers)

      if (deb_Path) write(*,*) '   << commposrespa'
      call timer_exit( timer_commstep,quiet_timers )
      end
c
c     subroutine sendall : everybody gets all the positions, speed and accel of the system
c
      subroutine sendall
      use atomsMirror
      use domdec
      use moldyn
      use mpi
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,iglob,iloc,iproc,tag,ierr
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: reqsend(:),reqrec(:)
      real(r_p), allocatable :: postemp(:,:)
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (postemp(9,n))
      postemp = 0_ti_p
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
!$acc update device(glob) async
c
c     get the positions
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(postemp(1,bufbeg(iproc+1)),9*domlen(iproc+1),
     $   MPI_RPREC,iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,9*nloc,MPI_RPREC,iproc,tag,
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
!$acc update device(x,y,z,v,a) async
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (postemp)
!$acc wait
      return
      end
c
c     subroutine sendallpos : everybody gets all the positions
c
      subroutine sendallpos
      use atomsMirror
      use domdec
      use mpi
      implicit none
      integer i,iglob,iloc,iproc,tag,ierr
      integer idomlen,ibufbeg
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: reqsend(:),reqrec(:)
      real(r_p), allocatable :: postemp(:,:)
c
      if (nproc.eq.1) return

      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (postemp(3,n))
!$acc data create(postemp) present(x,y,z,glob) async
c
c     put the arrays to be sent in the right order
c
!$acc parallel loop collapse(2) async
      do i = 1, n
         do iloc = 1, 3
            postemp(iloc,i) = 0
         end do
      end do

!$acc parallel loop async
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
     $       COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,iproc,tag,
     $       COMM_TINKER,reqsend(iproc+1),ierr)
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
!$acc host_data use_device(glob)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(glob(bufbeg(iproc+1)),domlen(iproc+1),MPI_INT,
     $       iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
!$acc wait
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $       COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
!$acc end host_data

      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     get the positions
c
!$acc host_data use_device(postemp)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(postemp(1,bufbeg(iproc+1)),3*domlen(iproc+1),
     $       MPI_RPREC,iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(postemp,3*nloc,MPI_RPREC,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
!$acc end host_data

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
          ibufbeg = bufbeg(iproc+1)
!$acc parallel loop async
          do i = 1, domlen(iproc+1)
            iloc  = ibufbeg+i-1
            iglob = glob(iloc)
            x(iglob) = postemp(1,iloc)
            y(iglob) = postemp(2,iloc)
            z(iglob) = postemp(3,iloc)
          end do
        end if
      end do
c
!$acc end data
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (postemp)
      end
c
c     subroutine sendvecmin : deal with communications during minimization routine
c
      subroutine sendvecmin(xx)
      use atomsMirror
      use domdec
      use mpi
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,iglob,iloc,iproc,tag,ierr,j
      integer ibufbeg
      integer status(MPI_STATUS_SIZE),count
      integer, allocatable :: reqsend(:),reqrec(:)
      real(r_p), allocatable :: temp(:,:)
      real(r_p) xx(3*n)
c
      if (nproc.eq.1) return

      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (temp(3,n))
!$acc data create(temp) present(xx) async

!$acc parallel loop collapse(2) async
      do i = 1, n
         do j = 1, 3
            temp(j,i) = 0.0_re_p
         end do
      end do
c
c     put the arrays to be sent in the right order
c
!$acc parallel loop collapse(2) async
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
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
!$acc host_data use_device(glob)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(glob(bufbeg(iproc+1)),domlen(iproc+1),MPI_INT,
     $   iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
!$acc wait
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,iproc,tag,
     $    COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
!$acc end host_data
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
c     get the values of xx
c
!$acc host_data use_device(temp)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(temp(1,bufbeg(iproc+1)),3*domlen(iproc+1),
     $       MPI_RPREC,iproc,tag,COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(temp,3*nloc,MPI_RPREC,iproc,tag,
     $       COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
!$acc end host_data
c
c     put it back in the global arrays
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          ibufbeg = bufbeg(iproc+1)
!$acc parallel loop present(glob,xx,temp) async
          do i = 1, domlen(iproc+1)
            iloc  = ibufbeg+i-1
            iglob = glob(iloc)
            xx(3*(iglob-1)+1) = temp(1,iloc)
            xx(3*(iglob-1)+2) = temp(2,iloc)
            xx(3*(iglob-1)+3) = temp(3,iloc)
          end do
        end if
      end do
c
!$acc end data
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (temp)
      end
c
c     subroutine commforces : communicate forces between processes
c
      subroutine commforces(derivs)
      use domdec,only:nbloc,nproc
      use timestat
      use mpi
      implicit none
      real(r_p) derivs(3,nbloc)
      real*8 time0,time1,time2
c
      call timer_enter( timer_fcomm )
      call commforcesrec(derivs)
      call commforcesbloc(derivs)
      call timer_exit( timer_fcomm,quiet_timers )
      end
c
c     subroutine commforcesrespa : communicate some forces between processes
c     when respa integrator is used
c
      subroutine commforcesrespa(derivs,fast)
      use domdec,only:nbloc,nproc
      use inform,only:deb_Path
      use timestat
      use mpi
      implicit none
      integer ialt,nalt
      real(r_p) derivs(3,nbloc)
      real*8 time0,time1,time2
      logical fast
c
      call timer_enter( timer_fcomm )
      if (deb_Path) write(*,*) ' -- commforcesrespa'
      if (fast) then
        call commforcesbonded(derivs)
      else
        call commforcesrec(derivs)
        call commforcesbloc(derivs)
      end if
      call timer_exit( timer_fcomm,quiet_timers )
      end
c
c     subroutine commforcesrespa1 : communicate some forces between processes
c     when respa1 integrator is used
c
      subroutine commforcesrespa1(derivs,rule)
      use domdec
      use iounit
      use inform,only:deb_Path
      use timestat
      use mpi
      implicit none
      integer rule
      real(r_p) derivs(3,*)
      real*8 time0,time1,time2
 1000 format(' illegal rule in commforcesrespa1.')
c
c     rule = 0: fast part of the forces
c     rule = 1: intermediate part of the forces
c     rule = 2: slow part of the forces
c
      call timer_enter( timer_fcomm )
      if (deb_Path) write(*,*) ' -- commforcesrespa1'
      if (rule.eq.0) then
        call commforcesbonded(derivs)
      else if (rule.eq.1) then
        call commforcesshortbloc(derivs)
      else if (rule.eq.2) then
        call commforcesrec(derivs)
        call commforcesbloc(derivs)
      else
         if (rank.eq.0) write(iout,1000)
      end if
      call timer_exit( timer_fcomm,quiet_timers )
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
      use timestat
      implicit none
      integer ierr,commloc
      real(r_p) epot
      real(r_p) buffer1(6),buffer2(16)
      real(r_p) vir_copy(3,3)
c
      if (ndir.eq.1) return

      call timer_enter(timer_reduceen)
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
!$acc data create(buffer1,buffer2,vir_copy)
!$acc&     present_or_copy(epot)
!$acc&     present(eb,eba,eub,eopb,et,ept,ett,ebt,ea
!$acc&           ,eaa,eopd,eid,eit,ec,ev,em,ep,eg,ex
!$acc&           ,esum)
c
!$acc serial async
      buffer1(1) = ec
      buffer1(2) = em
      buffer1(3) = ep
      buffer1(4) = epot
      buffer1(5) = esum
      buffer1(6) = einter
!$acc end serial
c
!$acc host_data use_device(buffer1)
!$acc wait
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer1,6,MPI_RPREC,MPI_SUM,0,
     $       COMM_TINKER,ierr)
      else
        call MPI_REDUCE(buffer1,buffer1,6,MPI_RPREC,MPI_SUM,0,
     $       COMM_TINKER,ierr)
      end if
!$acc end host_data
c
!$acc serial async
      ec     = buffer1(1)
      em     = buffer1(2)
      ep     = buffer1(3)
      epot   = buffer1(4)
      esum   = buffer1(5)
      einter = buffer1(6)
!$acc end serial
c
c   MPI: get virial
c
!$acc kernels async
      vir_copy(:,:) = vir(:,:)
!$acc end kernels
c
!$acc host_data use_device(vir_copy)
!$acc wait
      call MPI_ALLREDUCE(MPI_IN_PLACE,vir_copy,9,MPI_RPREC,MPI_SUM,
     $                   COMM_TINKER,ierr)
!$acc end host_data
c
!$acc kernels async
      vir(:,:) = vir_copy(:,:)
!$acc end kernels
c
      if ((use_pmecore).and.(rank.gt.(ndir-1))) goto 10
c
!$acc serial async
      buffer2(1)  = eba
      buffer2(2)  = ea
      buffer2(3)  = eb
      buffer2(4)  = ev
      buffer2(5)  = eub
      buffer2(6)  = eaa
      buffer2(7)  = eopb
      buffer2(8)  = eopd
      buffer2(9)  = eid
      buffer2(10) = eit
      buffer2(11) = et
      buffer2(12) = ept
      buffer2(13) = ebt
      buffer2(14) = ett
      buffer2(15) = eg
      buffer2(16) = ex
!$acc end serial
c
!$acc host_data use_device(buffer2)
!$acc wait
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,buffer2,16,MPI_RPREC,MPI_SUM,0,
     $       commloc,ierr)
      else
        call MPI_REDUCE(buffer2,buffer2,16,MPI_RPREC,MPI_SUM,0,
     $       commloc,ierr)
      end if
!$acc end host_data
c
!$acc serial async
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
!$acc end serial

  10  continue
!$acc end data
      call timer_exit(timer_reduceen,quiet_timers)
      end
c
c     subroutine commforcesbonded : deal with communications of some
c     bonded forces modifier avec nlocnl ?
c
      subroutine commforcesbonded(derivs)
      use sizes
      use deriv     ,only: debond=>desum
      use domdec
      use inform    ,only: deb_Path
      use potent
      use mpi
      use timestat
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1
      implicit none
      integer i,j,k,tag,ierr,iproc
      integer iloc,iloc_beg,idomlen
      integer sz
      integer reqsend(nproc),reqrec(nproc)
      real(r_p), pointer :: buffer(:,:,:)
      real(r_p) derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if ((nproc.eq.1).or.(use_pmecore.and.rank.gt.ndir-1)) return
      if (deb_Path) write(*,*) '   >> commforcesbonded'
      call timer_enter ( timer_dirbondfcomm )
c
c     allocate some arrays
c
      sz = 3*max(1,nloc)*nneig_send
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nneig_send) => buffMpi_p1(1:sz)
!$acc enter data attach(buffer)
c      buffer = 0_ti_p
c
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nneig_send
         tag = nproc*rank + pneig_send(i) + 1
!$acc host_data use_device(buffer)
         call MPI_IRECV(buffer(1,1,i),3*nloc,
     $        MPI_RPREC,pneig_send(i),tag,
     $        COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, nneig_recep
        tag = nproc*pneig_recep(i) + rank + 1
!$acc host_data use_device(debond)
        call MPI_ISEND(debond(1,bufbeg(pneig_recep(i)+1)),
     $       3*domlen(pneig_recep(i)+1),MPI_RPREC,pneig_recep(i),
     $       tag,COMM_TINKER,reqsend(i),ierr)
!$acc end host_data
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
!$acc parallel loop gang vector collapse(2) async
!$acc&         present(derivs,buffer)
      do j = 1, nloc
         do k = 1, 3
            do i = 1, nneig_send
               derivs(k,j) = derivs(k,j) + buffer(k,j,i)
            end do
         end do
      end do
!$acc exit data detach(buffer) async

      if (deb_Path) write(*,*) '   << commforcesbonded'
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end

      ! In place version of commforcesbonded
      subroutine commforcesbonded1(derivs)
      use sizes
      use domdec
      use inform    ,only: deb_Path
      use potent
      use mpi
      use timestat
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1
      implicit none
      integer i,j,k,tag,ierr,iproc
      integer iloc,iloc_beg,idomlen
      integer sz
      integer reqsend(nproc),reqrec(nproc)
      real(r_p), pointer :: buffer(:,:,:)
      real(r_p) derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if ((nproc.eq.1).or.(use_pmecore.and.rank.gt.ndir-1)) return
      !if (deb_Path) write(*,*) '   >> commforcesbonded1'
      call timer_enter ( timer_dirbondfcomm )
c
c     allocate some arrays
c
      sz = 3*max(1,nloc)*nneig_send
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nneig_send) => buffMpi_p1(1:sz)
!$acc enter data attach(buffer)
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nneig_send
         tag = nproc*rank + pneig_send(i) + 1
!$acc host_data use_device(buffer)
         call MPI_IRECV(buffer(1,1,i),3*nloc,
     $        MPI_RPREC,pneig_send(i),tag,
     $        COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, nneig_recep
        tag = nproc*pneig_recep(i) + rank + 1
!$acc host_data use_device(derivs)
        call MPI_ISEND(derivs(1,bufbeg(pneig_recep(i)+1)),
     $       3*domlen(pneig_recep(i)+1),MPI_RPREC,pneig_recep(i),
     $       tag,COMM_TINKER,reqsend(i),ierr)
!$acc end host_data
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
!$acc parallel loop gang vector collapse(2) async
!$acc&         present(derivs,buffer)
      do j = 1, nloc
         do k = 1, 3
            do i = 1, nneig_send
               derivs(k,j) = derivs(k,j) + buffer(k,j,i)
            end do
         end do
      end do

!$acc exit data detach(buffer) async
      !if (deb_Path) write(*,*) '   << commforcesbonded1'
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end
c
c
c     subroutine commforcesbloc : deal with communications of forces by bloc
c     between two time steps
c
      subroutine commforcesbloc(derivs)
      use deriv
      use domdec
      use inform,only: deb_Path
      use mpi
      use potent
      use sizes
      use timestat
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1
      implicit none
      integer i,j,k,tag,ierr,sz
      integer reqsend(nproc),reqrec(nproc)
      real(r_p), pointer :: buffer(:,:,:)
      real(r_p) derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
      if (deb_Path) write(*,*) '   >> commforcesbloc'
      call timer_enter( timer_dirbondfcomm )
c
c     allocate some arrays
c
      sz = 3*max(1,nloc)*nbig_send
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nbig_send) => buffMpi_p1(1:sz)
!$acc enter data attach(buffer)
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
         tag = nproc*rank + pbig_send(i) + 1
!$acc host_data use_device(buffer)
         call MPI_IRECV(buffer(1,1,i),3*nloc,
     $        MPI_RPREC,pbig_send(i),tag,
     $        COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, nbig_recep
         tag = nproc*pbig_recep(i) + rank + 1
!$acc host_data use_device(desum)
         call MPI_ISEND(desum(1,bufbeg(pbig_recep(i)+1)),
     $        3*domlen(pbig_recep(i)+1),MPI_RPREC,pbig_recep(i),tag,
     $        COMM_TINKER,reqsend(i),ierr)
!$acc end host_data
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
!$acc parallel loop gang vector collapse(2) async
!$acc&         present(derivs,buffer)
      do j = 1, nloc
         do k = 1, 3
            do i = 1, nbig_send
              derivs(k,j) = derivs(k,j) + buffer(k,j,i)
            end do
         end do
      end do
c
!$acc exit data detach(buffer) async
      nullify(buffer)

      if (deb_Path) write(*,*) '   << commforcesbloc'
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end
      ! In place version of commforcesbloc
      subroutine commforcesbloc1(derivs)
      use deriv
      use domdec
      use inform,only: deb_Path
      use mpi
      use potent
      use sizes
      use timestat
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1
      implicit none
      integer i,j,k,tag,ierr,sz
      integer reqsend(nproc),reqrec(nproc)
      real(r_p), pointer :: buffer(:,:,:)
      real(r_p) derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
      !if (deb_Path) write(*,*) '   >> commforcesbloc1'
      call timer_enter( timer_dirbondfcomm )
c
c     allocate some arrays
c
      sz = 3*max(1,nloc)*nbig_send
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nbig_send) => buffMpi_p1(1:sz)
!$acc enter data attach(buffer)
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
         tag = nproc*rank + pbig_send(i) + 1
!$acc host_data use_device(buffer)
         call MPI_IRECV(buffer(1,1,i),3*nloc,
     $        MPI_RPREC,pbig_send(i),tag,
     $        COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : begin sending
c
!$acc wait
      do i = 1, nbig_recep
         tag = nproc*pbig_recep(i) + rank + 1
!$acc host_data use_device(derivs)
         call MPI_ISEND(derivs(1,bufbeg(pbig_recep(i)+1)),
     $        3*domlen(pbig_recep(i)+1),MPI_RPREC,pbig_recep(i),tag,
     $        COMM_TINKER,reqsend(i),ierr)
!$acc end host_data
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
!$acc parallel loop gang vector collapse(2) async
!$acc&         present(derivs,buffer)
      do j = 1, nloc
         do k = 1, 3
            do i = 1, nbig_send
              derivs(k,j) = derivs(k,j) + buffer(k,j,i)
            end do
         end do
      end do
c
!$acc exit data detach(buffer) async
      nullify(buffer)

      !if (deb_Path) write(*,*) '   << commforcesbloc1'
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end
c
c     subroutine commforcesshortbloc : deal with communications of short range forces by bloc
c     between two time steps
c
      subroutine commforcesshortbloc(derivs)
      use deriv
      use domdec
      use inform,only: deb_Path
      use mpi
      use potent
      use sizes
      use timestat
      use utilgpu,only: rec_queue
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1
      implicit none
      integer i,j,k,tag,ierr,sz
      integer, allocatable :: reqsend(:),reqrec(:)
      real(r_p), pointer :: buffer(:,:,:)
c     real(r_p), allocatable :: buffers(:,:)
      real(r_p) derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
      if (deb_Path) write(*,*) '   >> commforcesshortbloc'
      call timer_enter( timer_dirbondfcomm )
c
c     allocate some arrays
c
      sz = 3*max(1,nloc)*nbigshort_send
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nbigshort_send) => buffMpi_p1(1:sz)
!$acc enter data attach(buffer)

c
c     communicate forces
c
c     MPI : begin reception in buffer
c

!$acc host_data use_device(buffer)
      do i = 1, nbigshort_send
        tag = nproc*rank + pbigshort_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $       MPI_RPREC,pbigshort_send(i),tag,
     $       COMM_TINKER,reqrec(i),ierr)
      end do
!$acc end host_data
c
c     MPI : begin sending
c
!$acc wait(rec_queue)
!$acc host_data use_device(desum)
      do i = 1, nbigshort_recep
        tag = nproc*pbigshort_recep(i) + rank + 1
        call MPI_ISEND(desum(1,bufbeg(pbigshort_recep(i)+1)),
     $       3*domlen(pbigshort_recep(i)+1),
     $       MPI_RPREC,pbigshort_recep(i),tag,
     $       COMM_TINKER,reqsend(i),ierr)
      end do
!$acc end host_data
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

!$acc parallel loop gang vector collapse(2) async(rec_queue)
!$acc&         present(buffer,derivs)
      do j = 1,nloc
         do k = 1,3
            do i = 1, nbigshort_send
              derivs(k,j) = derivs(k,j) + buffer(k,j,i)
            end do
         end do
      end do
c
!$acc exit data detach(buffer) async(rec_queue)
      deallocate (reqsend)
      deallocate (reqrec)

      if (deb_Path) write(*,*) '   << commforcesshortbloc'
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end
      ! In-place version of commforcesshortbloc
      subroutine commforcesshortbloc1(derivs)
      use deriv
      use domdec
      use inform,only: deb_Path
      use mpi
      use potent
      use sizes
      use timestat
      use utilgpu,only: rec_queue
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1
      implicit none
      integer i,j,k,tag,ierr,sz
      integer, allocatable :: reqsend(:),reqrec(:)
      real(r_p), pointer :: buffer(:,:,:)
c     real(r_p), allocatable :: buffers(:,:)
      real(r_p) derivs(3,nbloc)
      integer status(MPI_STATUS_SIZE)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
      !if (deb_Path) write(*,*) '   >> commforcesshortbloc1'
      call timer_enter( timer_dirbondfcomm )
c
c     allocate some arrays
c
      sz = 3*max(1,nloc)*nbigshort_send
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nbigshort_send) => buffMpi_p1(1:sz)
!$acc enter data attach(buffer)

c
c     communicate forces
c
c     MPI : begin reception in buffer
c

!$acc host_data use_device(buffer)
      do i = 1, nbigshort_send
        tag = nproc*rank + pbigshort_send(i) + 1
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $       MPI_RPREC,pbigshort_send(i),tag,
     $       COMM_TINKER,reqrec(i),ierr)
      end do
!$acc end host_data
c
c     MPI : begin sending
c
!$acc wait(rec_queue)
!$acc host_data use_device(derivs)
      do i = 1, nbigshort_recep
        tag = nproc*pbigshort_recep(i) + rank + 1
        call MPI_ISEND(derivs(1,bufbeg(pbigshort_recep(i)+1)),
     $       3*domlen(pbigshort_recep(i)+1),
     $       MPI_RPREC,pbigshort_recep(i),tag,
     $       COMM_TINKER,reqsend(i),ierr)
      end do
!$acc end host_data
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

!$acc parallel loop gang vector collapse(2) async(rec_queue)
!$acc&         present(buffer,derivs)
      do j = 1,nloc
         do k = 1,3
            do i = 1, nbigshort_send
              derivs(k,j) = derivs(k,j) + buffer(k,j,i)
            end do
         end do
      end do
c
!$acc exit data detach(buffer) async(rec_queue)
      deallocate (reqsend)
      deallocate (reqrec)

      !if (deb_Path) write(*,*) '   << commforcesshortbloc1'
      call timer_exit( timer_dirbondfcomm,quiet_timers )
      end
c
c
c     subroutine commforce_one : deal with communications of forces by bloc
c     between two time steps
c
      subroutine commforce_one(force,rule)
      use sizes
      use deriv
      use domdec
      use potent
      use mpi
      implicit none
      integer i,j
c
      integer  ,intent(in)   ::rule
      real(r_p),intent(inout):: force(3,nbloc)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
c
c     allocate some arrays
c
      if      (rule.eq.0) then
         call commforcesbonded1(force)
      else if (rule.eq.1) then
         call commforcesshortbloc1(force)
      else if (rule.eq.2) then
         call commforcesbloc1(force)
      else
 12     format('   WARNING  !!!'
     &       ,/,'--- Unknown config for commforce_one routine ---')
         print 12
      end if
      end
c
c     subroutine commforcesrec : deal with communications of reciprocal forces
c
      subroutine commforcerecdir_one(dir,rec,fsum)
      use sizes
      use domdec
      use mpi
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_p1
      implicit none
      integer i,j,k,tag,ierr,iproc,iglob,ilocrec,ibufbeg
      integer sz
      integer reqsend(nproc),reqrec(nproc)
      real(r_p), pointer :: buffer(:,:,:)
      real(r_p) dir(3,nbloc),rec(3,nblocrec)
      real(r_p) fsum(3,nloc)
      integer status(MPI_STATUS_SIZE)

      !case sequential
      if (nproc.eq.1) then
!$acc parallel loop async present(fsum,dir,rec)
         do i = 1,nloc
            do j = 1,3
               fsum(j,i) = dir(j,i)+rec(j,i)
            end do
         end do
         return
      end if
c
c     allocate some arrays
c
      sz = 3*max(1,nloc)*nrecdir_send
      call prmem_requestm(buffMpi_p1,sz,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nrecdir_send) => buffMpi_p1(1:sz)
!$acc enter data attach(buffer)

c     allocate (buffers(3,max(1,nblocrec)))
c
c
c     buffer  = 0_re_p
c     buffers = 0_re_p
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_send
         if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
!$acc host_data use_device(buffer)
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $           MPI_RPREC,precdir_send(i),tag,
     $           COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
         end if
      end do
c
c     MPI : move in buffer
c
!$acc wait
      do i = 1, nrecdir_recep
         if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
!$acc host_data use_device(rec)
            call MPI_ISEND(rec(1,bufbegrec(precdir_recep(i)+1)),
     $           3*domlen(precdir_recep(i)+1),MPI_RPREC,
     $           precdir_recep(i),tag,COMM_TINKER,reqsend(i),ierr)
!$acc end host_data
         end if
      end do
c
      ! save direct forces to fsum
!$acc parallel loop collapse(2) present(fsum,dir) async
      do i = 1,nloc
         do j = 1, 3
            fsum(j,i) = dir(j,i)
         end do
      end do

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
!$acc data present(buffer,rec,fsum,locrec1,glob)
      do i = 1, nrecdir_send
         if (precdir_send(i).ne.rank) then
!$acc parallel loop collapse(2) async
            do j = 1, nloc
               do k = 1, 3
                  fsum(k,j) = fsum(k,j) + buffer(k,j,i)
               end do
            end do
         else
!$acc parallel loop collapse(2) async
            do j = 1, nloc
               do k = 1, 3
                  ilocrec   = locrec1(glob(j))
                  fsum(k,j) = fsum(k,j) + rec(k,ilocrec)
               end do
            end do
         end if
      end do
!$acc end data
c
!$acc exit data detach(buffer) async
      end subroutine

      subroutine commforcesrec(derivs)
      use deriv
      use domdec
      use mpi
      use potent ,pa=>PotentialAll
      use inform ,only: deb_Path
      use sizes
      use timestat ,only: timer_enter,timer_exit,quiet_timers
     &             ,timer_dirreccomm,timer_fmanage
      use tinMemory,only: prmem_requestm
      use utilcomm ,only: buffMpi_p1,buffMpi_p2
      implicit none
      integer i,j,k,tag,ierr,iproc,iglob,ilocrec
      integer sz1,sz2
      integer reqsend(nproc),reqrec(nproc)
      real(r_p) derivs(3,nbloc)
      real(r_p),pointer:: buffer(:,:,:)
      real(r_p),pointer:: buffers(:,:)
      integer status(MPI_STATUS_SIZE)

      !case sequential
      if (nproc.eq.1) then
         call timer_enter( timer_fmanage )
         if (pa) then
!$acc parallel loop collapse(2)
!$acc&         default(present) async
         do i = 1, nloc
            do j = 1, 3
               derivs(j,i) = derivs(j,i) + demrec(j,i) +
     $                       deprec(j,i) + decrec(j,i)
            end do
         end do
         else 
            if (use_polar) then
!$acc parallel loop collapse(2)
!$acc&         default(present) async
               do i = 1, nloc
                  do j = 1, 3
                     derivs(j,i) = derivs(j,i) + demrec(j,i) +
     $                             deprec(j,i)
                  end do
               end do
            else if (use_mpole) then
!$acc parallel loop collapse(2)
!$acc&         default(present) async
               do i = 1, nloc
                  do j = 1, 3
                     derivs(j,i) = derivs(j,i) + demrec(j,i)
                  end do
               end do
            end if
            if (use_charge) then
!$acc parallel loop collapse(2)
!$acc&         default(present) async
               do i = 1, nloc
                  do j = 1, 3
                     derivs(j,i) = derivs(j,i) + decrec(j,i)
                  end do
               end do
            end if
         end if
         call timer_exit ( timer_fmanage,quiet_timers )

         if (.not.dotstgrad) return
         goto 100
      end if

      if (deb_Path) write(*,*) '   >> commforcesrec'
      call timer_enter( timer_dirreccomm )
c
c     allocate some arrays
c
      sz1 = 3*max(1,nloc)*nrecdir_send
      sz2 = 3*max(1,nblocrec)
      call prmem_requestm(buffMpi_p1,sz1,async=.false.)
      call prmem_requestm(buffMpi_p2,sz2,async=.false.)
      buffer(1:3,1:max(1,nloc),1:nrecdir_send) => buffMpi_p1(1:sz1)
      buffers(1:3,1:max(1,nblocrec))           => buffMpi_p2(1:sz2)
!$acc enter data attach(buffer,buffers)

c
c     buffer  = 0_re_p
c     buffers = 0_re_p
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nrecdir_send
         if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
!$acc host_data use_device(buffer)
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $           MPI_RPREC,precdir_send(i),tag,
     $           COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
         end if
      end do
c
c     MPI : move in buffer
c
!$acc data present(buffer,buffers)
c
      if (pa) then
!$acc parallel loop collapse(2) async
!$acc&       present(demrec,deprec,decrec)
         do i = 1, nblocrec
            do j = 1, 3
               buffers(j,i)= demrec(j,i) + deprec(j,i) + decrec(j,i)
            end do
         end do
      else
         if (use_polar) then
!$acc parallel loop collapse(2) async
!$acc&       present(demrec,deprec)
            do i = 1, nblocrec
               do j = 1, 3
                  buffers(j,i)= demrec(j,i) + deprec(j,i)
               end do
            end do
         else if (use_mpole) then
!$acc parallel loop collapse(2) async
!$acc&       present(demrec)
            do i = 1, nblocrec
               do j = 1, 3
                  buffers(j,i)= demrec(j,i)
               end do
            end do
         end if
         if (use_charge) then
!$acc parallel loop collapse(2) async
!$acc&       present(decrec)
            do i = 1, nblocrec
               do j = 1, 3
                  buffers(j,i)=  decrec(j,i)
               end do
            end do
         end if
      end if
c
!$acc wait
      do i = 1, nrecdir_recep
         if (precdir_recep(i).ne.rank) then
            tag = nproc*precdir_recep(i) + rank + 1
            sz1 = 3*domlen(precdir_recep(i)+1)
!$acc host_data use_device(buffers)
            call MPI_ISEND(buffers(1,bufbegrec(precdir_recep(i)+1)),
     $           sz1,MPI_RPREC,
     $           precdir_recep(i),tag,COMM_TINKER,reqsend(i),ierr)
!$acc end host_data
            !write(0,*) 'iRank',rank,'send',sz1,'to ',precdir_recep(i)
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
!$acc parallel loop collapse(2) async
!$acc&      present(derivs)
            do j = 1, nloc
               do k = 1, 3
                  derivs(k,j) = derivs(k,j) + buffer(k,j,i)
               end do
            end do
         else
            if (use_polar.and..not.use_charge) then
!$acc parallel loop collapse(2) async
!$acc&      present(glob,locrec1,derivs,demrec,deprec)
            do j = 1, nloc
               do k = 1, 3
                  iglob       = glob(j)
                  ilocrec     = locrec1(iglob)
                  derivs(k,j) =  derivs(k,j) + demrec(k,ilocrec) +
     $                     deprec(k,ilocrec)
               end do
            end do
            else if (use_mpole.and..not.use_charge) then
!$acc parallel loop collapse(2) async
!$acc&      present(glob,locrec1,derivs,demrec)
            do j = 1, nloc
               do k = 1, 3
                  iglob       = glob(j)
                  ilocrec     = locrec1(iglob)
                  derivs(k,j) =  derivs(k,j) + demrec(k,ilocrec)
               end do
            end do
            else if (.not.use_polar.and..not.use_mpole.and.use_charge)
     &         then
!$acc parallel loop collapse(2) async
!$acc&      present(glob,locrec1,derivs,decrec)
            do j = 1, nloc
               do k = 1, 3
                  iglob       = glob(j)
                  ilocrec     = locrec1(iglob)
                  derivs(k,j) =  derivs(k,j) + decrec(k,ilocrec)
               end do
            end do
            else
!$acc parallel loop collapse(2) async
!$acc&      present(glob,locrec1,derivs,demrec,deprec,decrec)
            do j = 1, nloc
               do k = 1, 3
                  iglob       = glob(j)
                  ilocrec     = locrec1(iglob)
                  derivs(k,j) =  derivs(k,j) + demrec(k,ilocrec) +
     $                     deprec(k,ilocrec) + decrec(k,ilocrec)
               end do
            end do
            end if
         end if
      end do
!$acc end data

      call timer_exit( timer_dirreccomm,quiet_timers )
c
c
c    Communicate separately the direct and reciprocal torques if testgrad is
c    being used
c
 100  continue
      if (dotstgrad) then
        ! case sequential
        if (nproc.eq.1) return
        call timer_enter( timer_dirreccomm )
c
c     MPI : begin reception in buffer
c
        if (use_mpole) then
        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_RPREC,precdir_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $       MPI_RPREC,precdir_recep(i),tag,COMM_TINKER,
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
        
        end if
c
        if (use_polar) then

        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_RPREC,precdir_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $       MPI_RPREC,precdir_recep(i),tag,COMM_TINKER,
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

        end if
c
        if (use_charge) then

        do i = 1, nrecdir_send
          if (precdir_send(i).ne.rank) then
            tag = nproc*rank + precdir_send(i) + 1
            call MPI_IRECV(buffer(1,1,i),3*nloc,
     $      MPI_RPREC,precdir_send(i),tag,COMM_TINKER,reqrec(i),ierr)
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
     $       MPI_RPREC,precdir_recep(i),tag,COMM_TINKER,
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

        end if
        call timer_exit( timer_dirreccomm,quiet_timers )
      end if
!$acc exit data detach(buffer,buffers)

      if (deb_Path) write(*,*) '   << commforcesrec'
      end
c
c     subroutine Icommforceststgrad : deal with communications of
c                                     one force (de) for testgrad
c
      subroutine Icommforceststgrad(detot,de,buffer,buffers,
     &                                       reqrec,reqsend)
      use sizes
      use domdec
      use mpi
      implicit none
      real(r_p),intent(inout)::detot(3,*),de(3,*)
      real(r_p),intent(  out)::buffer(3,max(1,nloc),nbig_send)
      real(r_p),intent(  out)::buffers(3,nbloc)
      integer,  intent(  out)::reqsend(nproc),reqrec(nproc)
      integer i,j,k,tag,ierr
      integer status(MPI_STATUS_SIZE)
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
!$acc data present(detot,de,buffer,buffer)
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
!$acc host_data use_device(buffer)
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $   MPI_RPREC,pbig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : move in buffer
c
!$acc parallel loop collapse(2)
      do i = nloc+1,nbloc
         do j = 1, 3
            buffers(j,i) = de(j,i)
         end do
      end do
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
!$acc host_data use_device(buffers)
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $       3*domlen(pbig_recep(i)+1),
     $       MPI_RPREC,pbig_recep(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
!$acc end host_data
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
!$acc parallel loop collapse(2)
        do j = 1, nloc
           do k = 1, 3
             de   (k,j) = de   (k,j) + buffer(k,j,i)
             detot(k,j) = detot(k,j) + buffer(k,j,i)
          end do
        end do
      end do
!$acc end data
      end subroutine

c
c     subroutine commforceststgrad : deal with communications of forces one by one
c     for testgrad
c
      subroutine commforceststgrad(detot)
      use sizes
      use deriv
      use domdec
      use mpi
      use tinheader
      implicit none
      integer i,j,k,tag,ierr
      integer, allocatable :: reqsend(:),reqrec(:)
      real(r_p), allocatable :: buffer(:,:,:)
      real(r_p), allocatable :: buffers(:,:)
      real(r_p) detot(3,*)
      integer status(MPI_STATUS_SIZE)
c     FIXME portage en cours
c
c     allocate some arrays
c
c      allocate (req(nproc*nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
c
      allocate (buffer(3,max(1,nloc),nbig_send))
      allocate (buffers(3,nbloc))
!$acc data present(detot)
!$acc&     present(deb)
!$acc&     create(buffer,buffers)

!$acc kernels
      buffer(1:3,1:nloc,1:nbig_send) = 0.0_ti_p
      buffers(1:3,1:nbloc)  = 0.0_ti_p
!$acc end kernels
c
c     communicate forces
c
c     MPI : begin reception in buffer
c
      do i = 1, nbig_send
        tag = nproc*rank + pbig_send(i) + 1
!$acc host_data use_device(buffer)
        call MPI_IRECV(buffer(1,1,i),3*nloc,
     $       MPI_RPREC,pbig_send(i),tag,
     $       COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
      end do
c
c     MPI : move in buffer
c
!$acc parallel loop collapse(2)
      do i = nloc+1,nbloc
         do j = 1, 3
            buffers(j,i) = deb(j,i)
         end do
      end do
c
      do i = 1, nbig_recep
        tag = nproc*pbig_recep(i) + rank + 1
!$acc host_data use_device(buffers)
        call MPI_ISEND(buffers(1,bufbeg(pbig_recep(i)+1)),
     $       3*domlen(pbig_recep(i)+1),
     $       MPI_RPREC,pbig_recep(i),tag,COMM_TINKER,
     $       reqsend(i),ierr)
!$acc end host_data
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
!$acc parallel loop collapse(2)
        do j = 1, nloc
           do k = 1, 3
             deb  (k,j) = deb  (k,j) + buffer(k,j,i)
             detot(k,j) = detot(k,j) + buffer(k,j,i)
          end do
        end do
      end do

      call Icommforceststgrad(detot,dea  ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deba ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deub ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deaa ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deopb,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deopd,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deid ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deit ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,det  ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,dept ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,debt ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,dett ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,dev  ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,dec  ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,dem  ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,dep  ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,deg  ,buffer,buffers,
     &                                    reqrec,reqsend)
      call Icommforceststgrad(detot,dex  ,buffer,buffers,
     &                                    reqrec,reqsend)

!$acc update host(buffer,buffers)
!$acc update host(detot(1:3,1:nbloc))
!$acc end data
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
      use domdec ,only: COMM_TINKER,nproc,ndir
      use energi
      use inter
      use mpi
      use potent ,only: use_pmecore
      implicit none
      integer ierr
      real(t_p) epot
c
      if (nproc.eq.1.or.use_pmecore.and.ndir.eq.1) return

!$acc host_data use_device(eba,ea,eb,ec,em,ep,ev,eub,eaa,eopb,eopd,
!$acc&          eid,eit,et,ept,ebt,ett,esum,epot,eg,ex,einter)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eba,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ea,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eb,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ec,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,em,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ep,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ev,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eub,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eaa,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eopb,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eopd,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eid,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eit,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,et,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ept,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ebt,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ett,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,esum,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,epot,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,eg,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,ex,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE,einter,1,MPI_RPREC,MPI_SUM,
     $    COMM_TINKER,ierr)
!$acc end host_data
      end
c
c     subroutine orderbuffer : get the indexes of the particules in the neighboring processes
c
      subroutine orderbuffer(init)
      use sizes
      use atoms  ,only: n
      use inform ,only: deb_Path
      use cutoff
      use domdec
      use mpi
      implicit none
      integer, allocatable :: counttemp(:),ind1temp(:)
      integer count
      integer i,iproc,iproc1,tag,ierr
      integer iloc,iglob
      integer, allocatable :: reqsend(:),reqrec(:)
      integer status(MPI_STATUS_SIZE)
      logical init
c
      allocate (counttemp(nproc))
      allocate (ind1temp(n))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      if (deb_Path) write(*,*) '   >> orderbuffer',init
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
c     nvdwblocloop is nvdwbloc if nvdwbloc is a multiple of 16, or the first one greater
      nblocloop = merge( nbloc,
     &                     (int(nbloc/16)+1)*16,
     &                     (mod(nbloc,16).eq.0))

      nblocrecdir = count
!!$acc update device(nbloc,nblocloop) async
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
!$acc update device(repart,loc,glob) async
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (ind1temp)
      deallocate (counttemp)

      if (deb_Path) write(*,*) '   << orderbuffer'
      end

c
c     Orderbuffer with openacc directives
c
      subroutine orderbuffer_gpu(init)
      use sizes
      use atoms  , only : n
      use cutoff
      use domdec
      use inform , only : deb_Path
      use mpi
      use utilgpu, only : openacc_abort
      implicit none
      integer :: counttemp(nproc)
      integer count
      integer i,iproc,iproc1,tag,ierr
      integer iloc,iglob
      integer idomlen,iprec,ibeg
      integer :: reqsend(nproc),reqrec(nproc)
      integer status(MPI_STATUS_SIZE)
      integer cap,cpt
      logical init

      ! This approach MPI_Allgather is fine when total number of prcocess
      ! doesn't exceed a certain amount
      ! Consider using classical method with more than 64 MPI process

      if (deb_Path) write(*,*) '   >> orderbuffer_gpu',init
      ! Get domlen values of all processes
      call MPI_Allgather( MPI_IN_PLACE, 1, MPI_DATATYPE_NULL,
     &                    domlen      , 1, MPI_INT,
     &                    COMM_TINKER, ierr )

      if (init) then
        cpt = 0
!$acc parallel loop present(repart,glob,loc) copyin(cpt)
        do i = 1, n
          !check if atom is in the domain or in one of the neighboring domains
          if (repart(i).eq.rank) then
!$acc atomic capture
             cpt       = cpt + 1
             cap       = cpt
!$acc end atomic
             glob(cap) = i
             loc(i)    = cap
          end if
        end do
      end if

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
c     nvdwblocloop is nvdwbloc if nvdwbloc is a multiple of 16, or the first one greater
      nblocloop = merge( nbloc,
     &                     (int(nbloc/16)+1)*16,
     &                     (mod(nbloc,16).eq.0))

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

      !send and receive the indexes of the neighboring processes
      do iproc = 1, nbig_recep
        tag = nproc*rank + pbig_recep(iproc) + 1
!$acc host_data use_device(glob)
        call MPI_IRECV(glob(bufbeg(pbig_recep(iproc)+1)),
     $   domlen(pbig_recep(iproc)+1),MPI_INT,pbig_recep(iproc),tag,
     $   COMM_TINKER,reqrec(iproc),ierr)
!$acc end host_data
      end do
!$acc wait
      do iproc = 1, nbig_send
        tag = nproc*pbig_send(iproc) + rank + 1
!$acc host_data use_device(glob)
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,pbig_send(iproc),
     $   tag,COMM_TINKER,reqsend(iproc),ierr)
!$acc end host_data
      end do
      do iproc = 1, nbig_recep
        call MPI_WAIT(reqrec(iproc),status,ierr)
      end do
      do iproc = 1, nbig_send
        call MPI_WAIT(reqsend(iproc),status,ierr)
      end do

c
c     also get the indexes of recdir neighboring processes
c
      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
          tag = nproc*rank + precdir_recep(iproc) + 1
!$acc host_data use_device(glob)
          call MPI_IRECV(glob(bufbeg(precdir_recep(iproc)+1)),
     $     domlen(precdir_recep(iproc)+1),MPI_INT,
     $     precdir_recep(iproc),tag,COMM_TINKER,reqrec(iproc),ierr)
!$acc end host_data
        end if
      end do
      do iproc = 1, nrecdir_send
        if (precdir_send(iproc).ne.rank) then
          tag = nproc*precdir_send(iproc) + rank + 1
!$acc host_data use_device(glob)
          call MPI_ISEND(glob,domlen(rank+1),MPI_INT,
     $     precdir_send(iproc),tag,COMM_TINKER,reqsend(iproc),ierr)
!$acc end host_data
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

      !Update loc & repart
      idomlen = domlen(rank+1)+1

!$acc data present(repart,glob,loc)

!$acc parallel loop async
      do i = idomlen, nblocrecdir
        loc(glob(i)) = i
      end do

      do iproc = 1, nbig_recep
        iprec   = pbig_recep(iproc)
        idomlen = domlen(iprec+1)
        ibeg    = bufbeg(iprec+1)
!$acc parallel loop async
        do i = 1, idomlen
           repart(glob(ibeg+i-1)) = iprec
        end do
      end do

      do iproc = 1, nrecdir_recep
        if (precdir_recep(iproc).ne.rank) then
          iprec   = precdir_recep(iproc)
          idomlen = domlen(iprec+1)
          ibeg    = bufbeg(iprec+1)
!$acc parallel loop async
          do i = 1, idomlen
            repart(glob(ibeg+i-1)) = iprec
          end do
        end if
      end do
!$acc end data

      if (deb_Path) write(*,*) '   << orderbuffer_gpu'
      end
c
c
c     subroutine orderbufferrec : get the indexes of the particules in the reciprocal processes
c
      subroutine orderbufferrec
      use atoms  , only : n
      use domdec
      use inform , only : deb_Path
      use potent
      implicit none
      integer counttemp(nproc)
      integer,allocatable::ind1temp(:)
      integer count,rankloc
      integer i,iproc
      integer iloc,iglob
      integer iprec,ind1
c
      if (deb_Path) write(*,*) '   >> orderbufferrec'
c     allocate (ind1temp(n))
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
!$acc wait
c
      if ((.not.(use_pmecore)).or.((use_pmecore).and.(rank.gt.ndir-1)))
     $  then
        do iproc = 1, nrecdir_recep
          do i = 1, domlen(precdir_recep(iproc)+1)
            iloc  = bufbeg(precdir_recep(iproc)+1)+i-1
            iglob = glob(iloc)
            domlenrec(repartrec(iglob)+1) =
     $          domlenrec(repartrec(iglob)+1) + 1
            ind1  = domlenrec(repartrec(iglob)+1)
c
c       check if atom is in the domain or in one of the neighboring domains
c
            if (repartrec(iglob).eq.rankloc) then
              globrec(ind1) = iglob
              locrec(iglob) = ind1
            end if
          end do
        end do
!$acc update device(globrec) async
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
      do iproc = 1, nrecdir_recep
        do i = 1, domlen(precdir_recep(iproc)+1)
          globrec1(bufbegrec(precdir_recep(iproc)+1)+i-1) =
     $       glob(bufbeg(precdir_recep(iproc)+1)+i-1)
          locrec1(globrec1(bufbegrec(precdir_recep(iproc)+1)+i-1)) =
     $      bufbegrec(precdir_recep(iproc)+1)+i-1
        end do
      end do
!$acc update device(locrec1)
c
c     deallocate (ind1temp)

      if (deb_Path) write(*,*) '   << orderbufferrec'
      end

      subroutine orderbufferrec_gpu
      use atoms  , only : n
      use domdec
      use inform , only : deb_Path
      use potent
      implicit none
c     integer, allocatable :: counttemp(:),ind1temp(:)
      integer count,rankloc
      integer i,iproc
      integer idomlen,ibeg,ibeg1,iprec,irep,ind1
      integer iloc,iglob
c
      if (deb_Path) write(*,*) '   >> orderbufferrec_gpu'
c
      if (use_pmecore) then
        rankloc = rank_bis
      else
        rankloc = rank
      end if
c
!$acc data present(domlenrec,globrec,globrec1
!$acc&         ,repartrec,locrec,locrec1,glob)

!$acc parallel loop async
      do i = 1,nproc
         domlenrec(i)=0
      end do
c
      if ((.not.(use_pmecore)).or.((use_pmecore).and.(rank.gt.ndir-1)))
     $ then
        do iproc = 1, nrecdir_recep
          iprec   = precdir_recep(iproc)
          idomlen = domlen(iprec+1)
          ibeg    = bufbeg(iprec+1)
!$acc parallel loop async
          do i = 1, idomlen
            iloc  = ibeg+i-1
            iglob = glob(iloc)
            irep  = repartrec(iglob)+1
!$acc atomic capture
            domlenrec(irep) = domlenrec(irep) + 1
            ind1            = domlenrec(irep)
!$acc end atomic
c
c       check if atom is in the domain or in one of the neighboring domains
c
            if (irep.eq.rankloc+1) then
              globrec(ind1) = iglob
              locrec(iglob) = ind1
            end if
          end do
        end do
      end if
!$acc update host(domlenrec) async
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
      do iproc  = 1, nrecdir_recep
        iprec   = precdir_recep(iproc)
        idomlen = domlen(iprec+1)
        ibeg    = bufbegrec(iprec+1)
        ibeg1   = bufbeg(iprec+1)
!$acc parallel loop async
        do i = 1, idomlen
          globrec1(ibeg+i-1) = glob(ibeg1+i-1)
          locrec1(globrec1(ibeg+i-1)) = ibeg+i-1
        end do
      end do
!$acc end data
c
!$acc wait
      ! Update nlocrec
      if (use_pmecore) then
        nlocrec = domlenrec(rank_bis+1)
      else
        nlocrec = domlenrec(rank+1)
      end if

      if (deb_Path) write(*,*) '   << orderbufferrec_gpu'
      end
c
c    subroutine commdirdir : deal with communications of direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles
c        - 2: wait for the communications to be done
c
      subroutine commdirdirgpu(nrhs,rule,mu,reqrec1,reqsend1)
      use domdec
      use iounit
      use inform  ,only: deb_Path
      use mpole
      use mpi
      use utilgpu ,only:rec_queue
      use sizes   ,only:tinkerdebug
      use utilcomm,only:reqsend=>reqs_dirdir,reqrec=>reqr_dirdir
     &            ,no_commdir
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,tag0,i
      logical lexit
      integer idomlen,ibufbeg,ipr
      integer :: reqrec1(nproc),reqsend1(nproc)
c     real(t_p) ef(3,nrhs,npolebloc)
      real(t_p) mu(3,nrhs,npolebloc)
      parameter(tag0=0,lexit=.false.)
 1000 format(' illegal rule in commdirdirgpu.')
 41   format(7x,'>> ',A20,   3x,'recv')
 42   format(7x,   3x,A20,' >>','send')
 43   format(7x,'<< ',A20,   3x,'wait')

      if ((n_recep1.eq.0.and.n_send1.eq.0).or.lexit) return

      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
         if (deb_Path) write(0,41) 'commdirdirgpu       '
c        if (btest(tinkerdebug,0).and.nproc.eq.4) then
c22         format('___',A6,I5,3x,9I5)
c           write(0,22)
c    &      'domlenpole',rank,npoleloc,p_send1(1:4),p_recep1(1:4)
c        end if
c
c     MPI : begin reception
c
         do i = 1, n_recep1
            tag     = nproc*rank + p_recep1(i) + 1
            ibufbeg = bufbegpole(p_recep1(i)+1)
            idomlen = domlenpole(p_recep1(i)+1)
!$acc host_data use_device(mu)
            call MPI_IRECV(mu(1,1,ibufbeg),3*nrhs*idomlen
     &          ,MPI_TPREC,p_recep1(i),tag0
     &          ,COMM_TINKER,reqrec(i),ierr)
!$acc end host_data
         end do
      else if (rule.eq.1) then
c
c     MPI : begin sending
c
!$acc wait(rec_queue)
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
!$acc host_data use_device(mu)
          call MPI_ISEND(mu,3*nrhs*npoleloc
     &        ,MPI_TPREC,p_send1(i),tag0
     &        ,COMM_TINKER,reqsend(i),ierr)
!$acc end host_data
        end do
         if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
         if (deb_Path) write(*,42) 'commdirdirgpu       '
      else if (rule.eq.2) then
         do i = 1, n_send1
            call MPI_WAIT(reqsend(i),status,ierr)
         end do
         do i = 1, n_recep1
            call MPI_WAIT(reqrec(i),status,ierr)
         end do
         if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
         if (deb_Path) write(*,43) 'commdirdirgpu       '
      else
         if (rank.eq.0) write(iout,1000)
         call fatal
      end if
      call timer_exit( timer_polsolvcomm,quiet_timers )
      end

      subroutine commdirdir(nrhs,rule,mu,reqrec,reqsend)
      use domdec
      use iounit
      use inform  ,only: deb_Path
      use mpole
      use mpi
      use sizes   ,only: tinkerdebug
      use utilcomm,only: no_commdir
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
      integer :: reqrec(nproc),reqsend(nproc)
c     real(t_p) ef(3,nrhs,npolebloc)
      real(t_p) mu(3,nrhs,npolebloc)
 1000 format(' illegal rule in commdirdir.')
c
      if (no_commdir) return

      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
c
c     MPI : begin reception
c
         do i = 1, n_recep1
            tag = nproc*rank + p_recep1(i) + 1
            call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     &           3*nrhs*domlen(p_recep1(i)+1),MPI_TPREC,p_recep1(i),
     &           tag,COMM_TINKER,reqrec(i),ierr)
         end do
      else if (rule.eq.1) then
c
c     MPI : begin sending
c
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
          call MPI_ISEND(mu,3*nrhs*npoleloc,
     &         MPI_TPREC,p_send1(i),tag,COMM_TINKER,
     &         reqsend(i),ierr)
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
      call timer_exit( timer_polsolvcomm,quiet_timers )
      end
c
c    subroutine commdirdirshort : deal with communications of short range direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles
c        - 2: wait for the communications to be done
c
      subroutine commdirdirshort(nrhs,rule,mu,reqrec1,reqsend1)
      use domdec
      use iounit
      use inform  ,only: deb_Path
      use mpole
      use mpi
      use utilgpu ,only:def_queue
      use utilcomm,only:reqsend=>reqs_dirdir,reqrec=>reqr_dirdir
     &            ,no_commdir
      use sizes
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
      integer ibufbeg
      integer :: reqrec1(nproc),reqsend1(nproc)
      real(t_p)  mu(3,nrhs,npolebloc)
 1000 format(' illegal rule in commdirdir.')
 41   format(7x,'>> ',A20,   3x,'recv')
 42   format(7x,   3x,A20,' >>','send')
 43   format(7x,'<< ',A20,   3x,'wait')
c
      if ((n_recepshort1.eq.0.and.n_sendshort1.eq.0)
     &   .or.no_commdir) return

      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
         if (deb_Path) write(0,41) 'commdirdirshort     '
c        if (btest(tinkerdebug,0).and.nproc.eq.4) then
c22         format('___',A6,I5,3x,10I7)
c           write(0,22)
c    &      'domlen',rank,npoleloc,domlen(1:4),npolebloc,bufbegpole(1:4)
c        end if
c
c     MPI : begin reception
c
!$acc host_data use_device(mu)
        do i = 1, n_recepshort1
          tag = nproc*rank + p_recepshort1(i) + 1
          ibufbeg = bufbegpole(p_recepshort1(i)+1)
          call MPI_IRECV(mu(1,1,ibufbeg),
     $         3*nrhs*domlenpole(p_recepshort1(i)+1),MPI_TPREC,
     $         p_recepshort1(i),tag,COMM_TINKER,reqrec(i),ierr)
        end do
!$acc end host_data
      else if (rule.eq.1) then
c
c     MPI : begin sending
c
!$acc wait(def_queue)
!$acc host_data use_device(mu)
        do i = 1, n_sendshort1
          tag = nproc*p_sendshort1(i) + rank + 1
          call MPI_ISEND(mu,3*nrhs*npoleloc,
     $         MPI_TPREC,p_sendshort1(i),tag,COMM_TINKER,
     $         reqsend(i),ierr)
        end do
!$acc end host_data
         if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
         if (deb_Path) write(*,42) 'commdirdirshort     '
      else if (rule.eq.2) then
         do i = 1, n_recepshort1
            call mpi_wait(reqrec(i),status,ierr)
         end do
         do i = 1, n_sendshort1
            call mpi_wait(reqsend(i),status,ierr)
         end do
         if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
         if (deb_Path) write(*,43) 'commdirdirshort     '
      else
         if (rank.eq.0) write(iout,1000)
         call fatal
      end if
      call timer_exit( timer_polsolvcomm,quiet_timers )
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
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      implicit none
      integer rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: reqrecdirrec(nproc),reqrecdirsend(nproc)
      real(t_p) efrec(10,max(1,npolerecloc))
      real(t_p) ef(10,max(1,npoleloc))
      real(t_p) buffermpi1(10,max(npoleloc,1))
      real(t_p) buffermpi2(10,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdirfields.')
c
      call timer_enter( timer_polsolvcomm )
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
     $        MPI_TPREC,proc,tag,COMM_TINKER,reqrecdirrec(i),ierr)
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
     $       10*buflen2(proc+1),MPI_TPREC,proc,tag,COMM_TINKER,
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
      call timer_exit( timer_polsolvcomm,quiet_timers )
      end
      subroutine commrecdirfieldsgpu(rule,efrec,ef,buffermpi1,buffermpi2
     $ ,reqrecdirrec,reqrecdirsend)
      use domdec
      use iounit
      use mpole
      use mpi
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      implicit none
      integer rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,k,proc
      integer ibufbeg
      integer lenbuf
      integer :: reqrecdirrec(nproc),reqrecdirsend(nproc)
      real(t_p) efrec(10,max(1,npolerecloc))
      real(t_p) ef(10,max(1,npoleloc))
      real(t_p) buffermpi1(10,max(npoleloc,1))
      real(t_p) buffermpi2(10,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdirfieldsgpu.')
 41   format(7x,'>> ',A20,   3x,'recv')
 42   format(7x,   3x,A20,' >>','send')
 43   format(7x,'<< ',A20,   3x,'wait')

      if (nrecdir_send1.eq.nrecdir_recep1.and.nrecdir_recep1.le.1)
     &   return
c
      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
c
c     Begin the reception of the reciprocal fields
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
!$acc host_data use_device(buffermpi1)
            call MPI_IRECV(buffermpi1(1,bufbeg1(proc+1)),
     $             10*buflen1(proc+1),MPI_TPREC,proc,tag,
     &                  COMM_TINKER,reqrecdirrec(i),ierr)
!$acc end host_data
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
             lenbuf = buflen2(proc+1)-1
             ibufbeg= bufbeg2(proc+1)
!$acc parallel loop gang collapse(2) async
!$acc&         default(present)
            do j = 0, lenbuf
               do k = 1, 10
               buffermpi2(k,ibufbeg+j)
     &               = efrec(k,polerecloc(buf2(ibufbeg+j)))
               end do
            end do
          end if
        end do
!$acc wait
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
!$acc host_data use_device(buffermpi2)
            call MPI_ISEND(buffermpi2(1,bufbeg2(proc+1)),
     $             10*buflen2(proc+1),MPI_TPREC,proc,tag,
     $                 COMM_TINKER,reqrecdirsend(i),ierr)
!$acc end host_data
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
            lenbuf = buflen1(proc+1)-1
            ibufbeg= bufbeg1(proc+1)
            call MPI_WAIT(reqrecdirrec(i),status,ierr)
!$acc parallel loop gang collapse(2) async
!$acc&         default(present)
            do j = 0, lenbuf
              do k = 1, 10
                ef(k,poleloc(buf1(ibufbeg+j))) = buffermpi1(k,ibufbeg+j)
              end do
            end do
          end if
        end do
!$acc wait
      else
        if (rank.eq.0) write(iout,1000)
        call fatal
      end if
      call timer_exit( timer_polsolvcomm,quiet_timers )
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
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: reqrecdirrec(nproc),reqrecdirsend(nproc)
      real(t_p) dipfield(3,nrhs,max(1,npoleloc))
      real(t_p) dipfieldbis(3,nrhs,max(1,npolerecloc))
      real(t_p) buffermpi1(3,nrhs,max(npoleloc,1))
      real(t_p) buffermpi2(3,nrhs,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdisolv.')
c
      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
           tag = nproc*rank + proc + 1
           call MPI_IRECV(buffermpi1(1,1,bufbeg1(proc+1)),
     $       3*nrhs*buflen1(proc+1),
     $       MPI_TPREC,proc,tag,COMM_TINKER,reqrecdirrec(i),ierr)
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
     $       MPI_TPREC,proc,tag,COMM_TINKER,reqrecdirsend(i),ierr)
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
      call timer_exit( timer_polsolvcomm,quiet_timers )
      end
      subroutine commrecdirsolvgpu(nrhs,rule,dipfieldbis,dipfield,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
      use domdec
      use iounit
      use inform  ,only:deb_Path
      use mpole
      use mpi
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      use sizes   ,only:tinkerdebug
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,k,l,proc
      integer lenbuf,ibufbeg
      integer reqrecdirrec(nproc),reqrecdirsend(nproc)
      real(t_p) dipfield(3,nrhs,max(1,npoleloc))
      real(t_p) dipfieldbis(3,nrhs,max(1,npolerecloc))
      real(t_p) buffermpi1(3,nrhs,max(npoleloc,1))
      real(t_p) buffermpi2(3,nrhs,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdirsolvgpu.')
 41   format(7x,'>> ',A20,   3x,'recv')
 42   format(7x,   3x,A20,' >>','send')
 43   format(7x,'<< ',A20,   3x,'wait')
      !nrhs==2 is a parameter

      if (nrecdir_send1.eq.nrecdir_recep1.and.nrecdir_recep1.le.1)
     &   return
c
      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
           tag = nproc*rank + proc + 1
!$acc host_data use_device(buffermpi1)
           call MPI_IRECV(buffermpi1(1,1,bufbeg1(proc+1)),
     $          3*nrhs*buflen1(proc+1),MPI_TPREC,proc,tag,
     $                   COMM_TINKER,reqrecdirrec(i),ierr)
!$acc end host_data
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
             lenbuf = buflen2(proc+1)-1
            ibufbeg = bufbeg2(proc+1)
!$acc parallel loop gang collapse(3) async
!$acc&         default(present)
             do j = 0, lenbuf
                do k = 1, 2
                   do l = 1, 3
                     buffermpi2(l,k,ibufbeg+j)= 
     &               dipfieldbis(l,k,polerecloc(buf2(ibufbeg+j)))
                   end do
                end do
             end do
          end if
        end do
c
!$acc wait
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
!$acc host_data use_device(buffermpi2)
            call MPI_ISEND(buffermpi2(1,1,bufbeg2(proc+1)),
     $           3*nrhs*buflen2(proc+1),MPI_TPREC,proc,tag,
     $                   COMM_TINKER,reqrecdirsend(i),ierr)
!$acc end host_data
          end if
        end do
        if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
        if (deb_Path) write(*,42) 'commrecdirsolvgpu   '
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
            lenbuf = buflen1(proc+1)-1
           ibufbeg = bufbeg1(proc+1)
!$acc parallel loop gang collapse(3) async
!$acc&         default(present)
            do j = 0, lenbuf
              do k = 1, 2
                 do l = 1, 3
                    dipfield(l,k,poleloc(buf1(ibufbeg+j)))
     $                       = buffermpi1(l,k,ibufbeg+j)
                 end do
              end do
            end do
          end if
        end do
!$acc wait
        if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
        if (deb_Path) write(*,43) 'commrecdirsolvgpu   '
      else
        if (rank.eq.0) write(iout,1000)
        call fatal
      end if
      call timer_exit( timer_polsolvcomm,quiet_timers )
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
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: req2rec(nproc),req2send(nproc)
      real(t_p) dip(3,nrhs,max(1,npoleloc))
      real(t_p) diprec(3,nrhs,max(1,npolerecloc))
c      real(t_p) buffermpimu1(3,nrhs,max(npoleloc,1))
c      real(t_p) buffermpimu2(3,nrhs,max(1,npolerecloc))
      real(t_p) buffermpimu1(3,nrhs,*)
      real(t_p) buffermpimu2(3,nrhs,*)
 1000 format(' illegal rule in commrecdirdip.')
c
      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),3*nrhs*
     $       buflen2(proc+1),MPI_TPREC,proc,tag,COMM_TINKER,
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
     $       buflen1(proc+1),MPI_TPREC,proc,tag,COMM_TINKER,
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
      call timer_exit( timer_polsolvcomm,quiet_timers )
      return
      end
      subroutine commrecdirdipgpu(nrhs,rule,diprec,dip,buffermpimu1,
     $           buffermpimu2,req2rec,req2send)
      use domdec
      use iounit
      use inform  ,only:deb_path
      use mpole
      use mpi
      use timestat,only:timer_enter,timer_exit,timer_polsolvcomm
     &            ,quiet_timers
      use sizes   ,only:tinkerdebug
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,k,l,proc
      integer lenbuf,ibufbeg
      integer :: req2rec(nproc),req2send(nproc)
      real(t_p) dip(3,nrhs,max(1,npoleloc))
      real(t_p) diprec(3,nrhs,max(1,npolerecloc))
c      real(t_p) buffermpimu1(3,nrhs,max(npoleloc,1))
c      real(t_p) buffermpimu2(3,nrhs,max(1,npolerecloc))
      real(t_p) buffermpimu1(3,nrhs,*)
      real(t_p) buffermpimu2(3,nrhs,*)
 1000 format(' illegal rule in commrecdirdipgpu.')
 41   format(7x,'>> ',A20,   3x,'recv')
 42   format(7x,   3x,A20,' >>','send')
 43   format(7x,'<< ',A20,   3x,'wait')
      
      if (nrecdir_send1.eq.nrecdir_recep1.and.nrecdir_recep1.le.1)
     &   return
c
      call timer_enter( timer_polsolvcomm )
      if (rule.eq.0) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
!$acc host_data use_device(buffermpimu2)
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),
     $             3*nrhs*buflen2(proc+1),MPI_TPREC,proc,tag,
     $                           COMM_TINKER,req2rec(i),ierr)
!$acc end host_data
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            lenbuf = buflen1(proc+1)-1
           ibufbeg = bufbeg1(proc+1)
!$acc parallel loop gang collapse(3) async
!$acc&         default(present)
            do j = 0,lenbuf
               do k = 1, 2
                 do l = 1, 3
                    buffermpimu1(l,k,ibufbeg+j)
     &                  = dip(l,k,poleloc(buf1(ibufbeg+j)))
                 end do
               end do
            end do
          end if
        end do
!$acc wait
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
!$acc host_data use_device(buffermpimu1)
            call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),
     $             3*nrhs*buflen1(proc+1),MPI_TPREC,proc,tag,
     $                          COMM_TINKER,req2send(i),ierr)
!$acc end host_data
          end if
        end do
        if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
        if (deb_Path) write(*,41) 'commrecdirdipgpu    '
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
            lenbuf = buflen2(proc+1)-1
           ibufbeg = bufbeg2(proc+1)
!$acc parallel loop gang collapse(3) async
!$acc&         default(present)
            do j = 0, lenbuf
               do k = 1, 2
                 do l = 1, 3
                   diprec(l,k,polerecloc(buf2(ibufbeg+j)))
     &                  = buffermpimu2(l,k,ibufbeg+j)
                 end do
              end do
            end do
          end if
        end do
!$acc wait
        if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
        if (deb_Path) write(*,42) 'commrecdirdipgpu    '
      else
        if (rank.eq.0) write(iout,1000)
        call fatal
      end if
      call timer_exit( timer_polsolvcomm,quiet_timers )
      end
c
c
c     subroutine commfield : communicate some direct fields (Newton's third law)
c
      subroutine commfield(nrhs,ef)
      use domdec
      use mpi
      use mpole
      use potent
      use timestat,only:timer_enter,timer_exit,timer_polfieldcomm
     &            ,quiet_timers
      implicit none
      integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer, allocatable :: reqrec(:), reqsend(:)
      real(t_p) ef(3,nrhs,*)
      real(t_p), allocatable :: buffer(:,:,:,:)
c
      call timer_enter( timer_polfieldcomm )
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = MPI_COMM_WORLD
      end if
c
      allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     $   p_send1(i),tag,commloc,reqrec(i),ierr)
      end do
c
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlenpole(p_recep1(i)+1),MPI_TPREC,p_recep1(i),
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
      call timer_exit( timer_polfieldcomm,quiet_timers )
      end
c
c     subroutine commfield : communicate some direct fields (Newton's third law)
c
      subroutine commfieldgpu(nrhs,ef)
      use domdec
      use mpi
      use mpole
      use potent
      use inform  ,only: deb_Path
      use sizes   ,only: tinkerdebug
      use timestat,only:timer_enter,timer_exit,timer_polfieldcomm
     &            ,quiet_timers
      use utilcomm,only:buff_field,no_commdir
      use utilgpu ,only:dir_queue,prmem_request
      implicit none
      integer  ,intent(in)    :: nrhs
      real(t_p),intent(inout) :: ef(3,nrhs,max(1,npolebloc))
      integer i,i1,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer ibufbeg
      integer :: reqrec(nproc), reqsend(nproc)
      !real(t_p)  :: buffer(3,nrhs,max(npoleloc,1),n_send1)
c
 41   format(7x,'>> ',A20,   3x)
 42   format(7x,   3x,A20,' >>')
 43   format(7x,'<< ',A20,   3x)

      if ( no_commdir ) return
      if (n_send1.gt.0.or.n_recep1.gt.0) then
      call timer_enter( timer_polfieldcomm )
      if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
      if (deb_Path) write(0,41) 'commfieldgpu        '
c        if (btest(tinkerdebug,0).and.nproc.eq.4) then
c22         format('_domlenpole ',I5,I7,3x,5I7,' _')
c           write(0,22) rank,npoleloc,domlenpole(1:4)
c        end if

         if (use_pmecore) then
            commloc = comm_dir
         else
            commloc = MPI_COMM_WORLD
         end if

         call prmem_request(buff_field,3,nrhs,
     &        max(npoleloc,1),n_send1,async=.true.)
c
        do i = 1, n_send1
           tag = nproc*rank + p_send1(i) + 1
!$acc host_data use_device(buff_field)
           call MPI_IRECV(buff_field(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     &          p_send1(i),tag,commloc,reqrec(i),ierr)
!$acc end host_data
        end do
c
!$acc wait(dir_queue)
        do i = 1, n_recep1
           tag = nproc*p_recep1(i) + rank + 1
           ibufbeg = bufbegpole(p_recep1(i)+1)
!$acc host_data use_device(ef)
           call MPI_ISEND(ef(1,1,ibufbeg),
     &          3*nrhs*domlenpole(p_recep1(i)+1),MPI_TPREC,
     &          p_recep1(i),tag,commloc,reqsend(i),ierr)
!$acc end host_data
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
!$acc parallel loop gang vector collapse(3)
!$acc&         present(ef,buff_field) async(dir_queue)
        do j = 1, npoleloc
           do k = 1, nrhs
              do i1=1,3
                 do i = 1, n_send1
                ef(i1,k,j) = ef(i1,k,j) + buff_field(i1,k,j,i)
                 end do
              end do
           end do
        end do

      if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
      if (deb_Path) write(*,43) 'commfieldgpu        '
      call timer_exit( timer_polfieldcomm,quiet_timers )
      end if
      end
c
c     subroutine commfieldshort : communicate some short range direct fields (newton's third law)
c
      subroutine commfieldshort(nrhs,ef)
      use domdec
      use mpi
      use mpole
      use potent
      use inform  ,only: deb_Path
      use sizes   ,only: tinkerdebug
      use timestat,only:timer_enter,timer_exit,timer_polfieldcomm
     &            ,quiet_timers
      use utilcomm ,only: buffer=>buff_field,no_commdir
      use utilgpu  ,only: prmem_request,def_queue
      implicit none
      integer nrhs,i,j,k,l,tag,status(mpi_status_size)
      integer ierr,ibufbeg,commloc
      integer, allocatable :: reqrec(:), reqsend(:)
      real(t_p) ef(3,nrhs,*)
c
      if (no_commdir) return
      if (n_sendshort1.eq.0.and.n_recepshort1.eq.0) return
 41   format(7x,'>> ',A20,   3x)
 42   format(7x,   3x,A20,' >>')
 43   format(7x,'<< ',A20,   3x)

      call timer_enter( timer_polfieldcomm )
      if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
      if (deb_Path) write(0,41) 'commfieldshort      '
c        if (btest(tinkerdebug,0).and.nproc.eq.4) then
c22         format('_',A6,I5,3x,5I7,'_ ')
c           write(0,22) 'npoleloc',rank,npoleloc
c        end if

      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      call prmem_request(buffer,3,nrhs,max(npoleloc,1),n_sendshort1,
     &                    async=.true.)
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
!$acc host_data use_device(buffer)
      do i = 1, n_sendshort1
        tag = nproc*rank + p_sendshort1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     $       p_sendshort1(i),tag,commloc,reqrec(i),ierr)
      end do
!$acc end host_data
c
!$acc wait(def_queue)
!$acc host_data use_device(ef)
      do i = 1, n_recepshort1
        tag     = nproc*p_recepshort1(i) + rank + 1
        ibufbeg = bufbegpole(p_recepshort1(i)+1)
        call MPI_ISEND(ef(1,1,ibufbeg),
     $       3*nrhs*domlenpole(p_recepshort1(i)+1),MPI_TPREC,
     $       p_recepshort1(i),tag,commloc,reqsend(i),ierr)
      end do
!$acc end host_data
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
!$acc parallel loop gang vector collapse(3)
!$acc&         default(present) async(def_queue)
      do j = 1, npoleloc
        do k = 1, nrhs
          do l = 1,3
            do i = 1, n_sendshort1
              ef(l,k,j) = ef(l,k,j) + buffer(l,k,j,i)
            end do
          end do
        end do
      end do
c
      deallocate (reqrec)
      deallocate (reqsend)

      if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)
      if (deb_Path) write(*,43) 'commfieldshort      '
      call timer_exit( timer_polfieldcomm,quiet_timers )
      end
