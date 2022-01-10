c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rattle  --  RATTLE distance constraint method  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rattle" implements the first portion of the RATTLE algorithm
c     by correcting atomic positions and half-step velocities to
c     maintain interatomic distance and absolute spatial constraints
c
c     literature reference:
c
c     H. C. Andersen, "RATTLE: A Velocity Version of the SHAKE
c     Algorithm for Molecular Dynamics Calculations", Journal of
c     Computational Physics, 52, 24-34 (1983)
c
c
      subroutine rattle (dt)
      use atmlst
      use atmtyp
      use atoms
      use domdec
      use freeze
      use inform
      use iounit
      use moldyn
      use usage
      use mpi
      implicit none
      integer i,j,k,iloc,iglob,ierr
      integer ia,ib,ialoc,ibloc,mode
      integer niter,maxiter
      integer start,stop
      real*8 dt,eps,sor
      real*8 xr,yr,zr
      real*8 xo,yo,zo
      real*8 xv,yv,zv
      real*8 dot,rma,rmb
      real*8 weigh,dist2
      real*8 delta,term
      real*8 xterm,yterm,zterm
      real*8, allocatable :: displace(:,:)
      logical done
      logical, allocatable :: moved(:)
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (moved(nbloc))
      allocate (update(nbloc))
      allocate (displace(3,nbloc))
      displace = 0d0
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, nbloc
         iglob = glob(i)
         if (use(iglob)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     initialize arrays for MPI 
c
      call initmpirattle
      call commrattleinit
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      sor = 1.25d0
      eps = rateps
c
c     apply RATTLE to distances and half-step velocity values
c
      niter = 0
      done = .false.
      do while (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         displace = 0d0
         do iloc = 1, nratloc
            i = ratglob(iloc) 
            ia = irat(1,i)
            ib = irat(2,i)
            ialoc = loc(ia)
            ibloc = loc(ib)
            if (moved(ialoc) .or. moved(ibloc)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (ratimage(i))  call image (xr,yr,zr)
               dist2 = xr**2 + yr**2 + zr**2
               delta = krat(i)**2 - dist2
               if (abs(delta) .gt. eps) then
                  done = .false.
                  update(ialoc) = .true.
                  update(ibloc) = .true.
                  xo = xold(ib) - xold(ia)
                  yo = yold(ib) - yold(ia)
                  zo = zold(ib) - zold(ia)
                  if (ratimage(i))  call image (xo,yo,zo)
                  dot = xr*xo + yr*yo + zr*zo
                  rma = 1.0d0 / mass(ia)
                  rmb = 1.0d0 / mass(ib)
                  term = sor * delta / (2.0d0 * (rma+rmb) * dot)
                  xterm = xo * term
                  yterm = yo * term
                  zterm = zo * term
                  x(ia) =  x(ia) - xterm*rma
                  y(ia) =  y(ia) - yterm*rma
                  z(ia) =  z(ia) - zterm*rma
                  x(ib) =  x(ib) + xterm*rmb
                  y(ib) =  y(ib) + yterm*rmb
                  z(ib) =  z(ib) + zterm*rmb
                  displace(1,ialoc) =  displace(1,ialoc) - xterm*rma
                  displace(2,ialoc) =  displace(2,ialoc) - yterm*rma
                  displace(3,ialoc) =  displace(3,ialoc) - zterm*rma
                  displace(1,ibloc) =  displace(1,ibloc) + xterm*rmb
                  displace(2,ibloc) =  displace(2,ibloc) + yterm*rmb
                  displace(3,ibloc) =  displace(3,ibloc) + zterm*rmb
                  rma = rma / dt
                  rmb = rmb / dt
                  v(1,ia) = v(1,ia) - xterm*rma
                  v(2,ia) = v(2,ia) - yterm*rma
                  v(3,ia) = v(3,ia) - zterm*rma
                  v(1,ib) = v(1,ib) + xterm*rmb
                  v(2,ib) = v(2,ib) + yterm*rmb
                  v(3,ib) = v(3,ib) + zterm*rmb
               end if
            end if
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,done,1,MPI_LOGICAL,MPI_LAND,
     $    COMM_BEAD,ierr)
c
         do i = 1, nbloc
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
      call commrattleend
c
c    faire routine qui envoie chaque partie des contraintes convergees a tous les voisins
c
c     perform deallocation of some local arrays
c
      deallocate (moved)
      deallocate (update)
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' RATTLE  --  Warning, Distance Constraints',
     &              ' not Satisfied')
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' RATTLE   --  Distance Constraints met at',i6,
     &              ' Iterations')
      end if
c
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rattle2  --  RATTLE atom velocity corrections  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rattle2" implements the second portion of the RATTLE algorithm
c     by correcting the full-step velocities in order to maintain
c     interatomic distance constraints
c
c
      subroutine rattle2 (dt)
      use atmlst
      use atmtyp
      use atoms
      use domdec
      use group
      use freeze
      use inform
      use iounit
      use moldyn
      use units
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,iloc,iglob,ierr
      integer ia,ib,ialoc,ibloc,mode
      integer niter,maxiter
      integer start,stop
      real*8 dt,eps,sor
      real*8 xr,yr,zr
      real*8 xv,yv,zv
      real*8 dot,rma,rmb
      real*8 weigh,vterm,term
      real*8 xterm,yterm,zterm
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: displace(:,:)
      real*8 virtemp(3,3)
      logical done
      logical, allocatable :: moved(:)
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (moved(nbloc))
      allocate (update(nbloc))
      allocate (displace(3,nbloc))
      displace = 0d0
c
      virtemp = 0d0
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, nbloc
         iglob = glob(i)
         if (use(iglob)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     initialize arrays for MPI 
c
      call initmpirattle
      call commrattleinit
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      niter = 0
      done = .false.
      sor = 1.25d0
      eps = rateps / dt
      vterm = 2.0d0 / (dt*convert)
c
c     apply the RATTLE algorithm to correct the velocities
c
      do while (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         displace = 0d0
         do iloc = 1, nratloc
            i = ratglob(iloc)
            ia = irat(1,i)
            ib = irat(2,i)
            ialoc = loc(ia)
            ibloc = loc(ib)
            if (moved(ialoc) .or. moved(ibloc)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (ratimage(i))  call image (xr,yr,zr)
               xv = v(1,ib) - v(1,ia)
               yv = v(2,ib) - v(2,ia)
               zv = v(3,ib) - v(3,ia)
               dot = xr*xv + yr*yv + zr*zv
               rma = 1.0d0 / mass(ia)
               rmb = 1.0d0 / mass(ib)
               term = -dot / ((rma+rmb) * krat(i)**2)
               if (abs(term) .gt. eps) then
                  done = .false.
                  update(ialoc) = .true.
                  update(ibloc) = .true.
                  term = sor * term
                  xterm = xr * term
                  yterm = yr * term
                  zterm = zr * term
                  v(1,ia) = v(1,ia) - xterm*rma
                  v(2,ia) = v(2,ia) - yterm*rma
                  v(3,ia) = v(3,ia) - zterm*rma
                  v(1,ib) = v(1,ib) + xterm*rmb
                  v(2,ib) = v(2,ib) + yterm*rmb
                  v(3,ib) = v(3,ib) + zterm*rmb
                  displace(1,ialoc) =  displace(1,ialoc) - xterm*rma
                  displace(2,ialoc) =  displace(2,ialoc) - yterm*rma
                  displace(3,ialoc) =  displace(3,ialoc) - zterm*rma
                  displace(1,ibloc) =  displace(1,ibloc) + xterm*rmb
                  displace(2,ibloc) =  displace(2,ibloc) + yterm*rmb
                  displace(3,ibloc) =  displace(3,ibloc) + zterm*rmb
c
c     increment the internal virial tensor components
c
                  xterm = xterm * vterm
                  yterm = yterm * vterm
                  zterm = zterm * vterm
                  vxx = xr * xterm
                  vyx = yr * xterm
                  vzx = zr * xterm
                  vyy = yr * yterm
                  vzy = zr * yterm
                  vzz = zr * zterm
                  virtemp(1,1) = virtemp(1,1) - vxx
                  virtemp(2,1) = virtemp(2,1) - vyx
                  virtemp(3,1) = virtemp(3,1) - vzx
                  virtemp(1,2) = virtemp(1,2) - vyx
                  virtemp(2,2) = virtemp(2,2) - vyy
                  virtemp(3,2) = virtemp(3,2) - vzy
                  virtemp(1,3) = virtemp(1,3) - vzx
                  virtemp(2,3) = virtemp(2,3) - vzy
                  virtemp(3,3) = virtemp(3,3) - vzz
               end if
            end if
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,done,1,MPI_LOGICAL,MPI_LAND,
     $    COMM_BEAD,ierr)
c
         do i = 1, nbloc
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
      call commrattleend
      call MPI_ALLREDUCE(MPI_IN_PLACE,virtemp,9,MPI_REAL8,MPI_SUM,
     $ COMM_BEAD,ierr)
      vir = vir + virtemp
c
c
c     perform deallocation of some local arrays
c
      deallocate (moved)
      deallocate (update)
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' RATTLE2  --  Warning, Velocity Constraints',
     &              ' not Satisfied')
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' RATTLE2  --  Velocity Constraints met at',i6,
     &              ' Iterations')
      end if
c
      return
      end
c
c    subroutine initmpirattle : build the arrays to communicate positions of atoms
c    during the calculation of the constrains using the rattle algorithm
c
c
      subroutine initmpirattle
      use atmlst
      use atoms
      use domdec
      use freeze
      use mpole
      use mpi
      implicit none
      integer ierr,iloc,iglob,j,proc
      integer i,iproc,tag,ia,ib
      integer count1
      integer status(MPI_STATUS_SIZE)
      integer, allocatable :: count(:)
      integer, allocatable :: reqrec(:),reqsend(:)
      allocate (count(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
c
      count = 0 
c
c     deal with constrains communications
c
      if (allocated(bufrat1)) deallocate (bufrat1)
      allocate (bufrat1(nbloc))
      if (allocated(bufrat2)) deallocate (bufrat2)
      allocate (bufrat2(nbloc))
      if (allocated(buflenrat1)) deallocate (buflenrat1)
      allocate (buflenrat1(nproc))
      buflenrat1 = 0
      if (allocated(buflenrat2)) deallocate (buflenrat2)
      allocate (buflenrat2(nproc))
      buflenrat2 = 0
      if (allocated(bufbegrat1)) deallocate (bufbegrat1)
      allocate (bufbegrat1(nproc))
      bufbegrat1 = 0
      if (allocated(bufbegrat2)) deallocate (bufbegrat2)
      allocate (bufbegrat2(nproc))
      bufbegrat2 = 0
c
      do iloc = 1, nratloc
        i = ratglob(iloc)
        ia = irat(1,i)
        ib = irat(2,i)
        if (repart(ia).ne.rank) then
          buflenrat2(repart(ia)+1) = buflenrat2(repart(ia)+1)+1
        end if
        if (repart(ib).ne.rank) then
          buflenrat2(repart(ib)+1) = buflenrat2(repart(ib)+1)+1
        end if
      end do
      count1 = 0
      do iproc = 1, nneig_recep
        if (pneig_recep(iproc).ne.rank) then
          if (buflenrat2(pneig_recep(iproc)+1).ne.0) then
            bufbegrat2(pneig_recep(iproc)+1) = count1 + 1
          else
            bufbegrat2(pneig_recep(iproc)+1) = 1
          end if
          count1 = count1 + buflenrat2(pneig_recep(iproc)+1)
        end if
      end do
c
      do iloc = 1, nratloc
        i = ratglob(iloc)
        ia = irat(1,i)
        ib = irat(2,i)
        if (repart(ia).ne.rank) then
          bufrat2(bufbegrat2(repart(ia)+1)+count(repart(ia)+1))=
     $      ia
          count(repart(ia)+1) = count(repart(ia)+1) + 1
        end if
        if (repart(ib).ne.rank) then
          bufrat2(bufbegrat2(repart(ib)+1)+count(repart(ib)+1))=
     $      ib
          count(repart(ib)+1) = count(repart(ib)+1) + 1
        end if
      end do
c
c     send and receive sizes of the buffers
c
       do i = 1, nneig_send
         if (pneig_send(i).ne.rank) then
          tag = nproc*rank + pneig_send(i) + 1
          call MPI_IRECV(buflenrat1(pneig_send(i)+1),1,MPI_INT,
     $   pneig_send(i),tag,COMM_BEAD,reqrec(i),ierr)
        end if
      end do
      do i = 1, nneig_recep
        if (pneig_recep(i).ne.rank) then
          tag = nproc*pneig_recep(i) + rank + 1
          call MPI_ISEND(buflenrat2(pneig_recep(i)+1),1,MPI_INT,
     $     pneig_recep(i),tag,COMM_BEAD,reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nneig_send
        if (pneig_send(i).ne.rank) then
          tag = nproc*rank + pneig_send(i) + 1
          call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nneig_recep
        if (pneig_recep(i).ne.rank) then
          tag = nproc*pneig_recep(i) + rank + 1
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
      count1 = 0
      do iproc = 1, nneig_send
        if (pneig_send(iproc).ne.rank) then
          if (buflenrat1(pneig_send(iproc)+1).ne.0) then
            bufbegrat1(pneig_send(iproc)+1) = count1 + 1
          else
            bufbegrat1(pneig_send(iproc)+1) = 1
          end if
          count1 = count1 + buflenrat1(pneig_send(iproc)+1)
        end if
      end do
c
c     send and receive list of corresponding indexes
c
      do i = 1, nneig_send
        if (pneig_send(i).ne.rank) then
          tag = nproc*rank + pneig_send(i) + 1
          call MPI_IRECV(bufrat1(bufbegrat1(pneig_send(i)+1)),
     $     buflenrat1(pneig_send(i)+1),
     $     MPI_INT,pneig_send(i),tag,COMM_BEAD,reqrec(i),ierr)
        end if
      end do
      do i = 1, nneig_recep
        if (pneig_recep(i).ne.rank) then
          tag = nproc*pneig_recep(i) + rank + 1
          call MPI_ISEND(bufrat2(bufbegrat2(pneig_recep(i)+1)),
     $     buflenrat2(pneig_recep(i)+1),MPI_INT,pneig_recep(i),tag,
     $     COMM_BEAD,reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nneig_send
        if (pneig_send(i).ne.rank) then
          tag = nproc*rank + pneig_send(i) + 1
          call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nneig_recep
        if (pneig_recep(i).ne.rank) then
          tag = nproc*pneig_recep(i) + rank + 1
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (count)
      return
      end
c
c     subroutine commrattleinit : deal with necessary comms before rattle iterations
c
      subroutine commrattleinit
      use atoms
      use domdec
      use freeze
      use moldyn
      use mpi
      implicit none
      integer i,j,proc,ierr,iglob,iloc
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqrec(:),reqsend(:)
      real*8, allocatable :: buffermpi1(:,:),buffermpi2(:,:)
c
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (buffermpi1(9,max(nbloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpi2(9,max(nbloc,1)))
      buffermpi2 = 0d0
c
c     Begin the reception of the constrained atoms
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_IRECV(buffermpi1(1,bufbegrat2(proc+1)),
     $      9*buflenrat2(proc+1),
     $      MPI_REAL8,proc,tag,COMM_BEAD,reqrec(i),ierr)
        end if
      end do
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        if (proc.ne.rank) then
          do j = 0, buflenrat1(proc+1)-1
            iglob = bufrat1(bufbegrat1(proc+1)+j)
            buffermpi2(1,bufbegrat1(proc+1)+j) = 
     $          x(iglob)
            buffermpi2(2,bufbegrat1(proc+1)+j) = 
     $          y(iglob)
            buffermpi2(3,bufbegrat1(proc+1)+j) = 
     $          z(iglob)
            buffermpi2(4,bufbegrat1(proc+1)+j) = 
     $          v(1,iglob)
            buffermpi2(5,bufbegrat1(proc+1)+j) = 
     $          v(2,iglob)
            buffermpi2(6,bufbegrat1(proc+1)+j) = 
     $          v(3,iglob)
            buffermpi2(7,bufbegrat1(proc+1)+j) = 
     $          xold(iglob)
            buffermpi2(8,bufbegrat1(proc+1)+j) = 
     $          yold(iglob)
            buffermpi2(9,bufbegrat1(proc+1)+j) = 
     $          zold(iglob)
          end do
        end if
      end do
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_ISEND(buffermpi2(1,bufbegrat1(proc+1)),
     $     9*buflenrat1(proc+1),MPI_REAL8,proc,tag,COMM_BEAD,
     $     reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_WAIT(reqrec(i),status,ierr)
          do j = 0, buflenrat2(proc+1)-1
            iglob = bufrat2(bufbegrat2(proc+1)+j)
            x(iglob) = buffermpi1(1,bufbegrat2(proc+1)+j)
            y(iglob) = buffermpi1(2,bufbegrat2(proc+1)+j)
            z(iglob) = buffermpi1(3,bufbegrat2(proc+1)+j)
            v(1,iglob) = buffermpi1(4,bufbegrat2(proc+1)+j)
            v(2,iglob) = buffermpi1(5,bufbegrat2(proc+1)+j)
            v(3,iglob) = buffermpi1(6,bufbegrat2(proc+1)+j)
            xold(iglob) = buffermpi1(7,bufbegrat2(proc+1)+j)
            yold(iglob) = buffermpi1(8,bufbegrat2(proc+1)+j)
            zold(iglob) = buffermpi1(9,bufbegrat2(proc+1)+j)
          end do
        end if
      end do
c
      deallocate (buffermpi2)
      deallocate (buffermpi1)
      deallocate (reqsend)
      deallocate (reqrec)
      return
      end
c
c     subroutine commrattleend : deal with necessary comms after rattle iterations
c
      subroutine commrattleend
      use atoms
      use domdec
      use freeze
      use moldyn
      use mpi
      implicit none
      integer i,j,proc,ierr,iglob,iloc
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqrec(:),reqsend(:)
      real*8, allocatable :: buffermpi1(:,:),buffermpi2(:,:)
c
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (buffermpi1(9,max(nbloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpi2(9,max(nbloc,1)))
      buffermpi2 = 0d0
c
c     Begin the reception of the constrained atoms
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_IRECV(buffermpi2(1,bufbegrat1(proc+1)),
     $      9*buflenrat1(proc+1),
     $      MPI_REAL8,proc,tag,COMM_BEAD,reqrec(i),ierr)
        end if
      end do
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        if (proc.ne.rank) then
          do j = 0, buflenrat2(proc+1)-1
            iglob = bufrat2(bufbegrat2(proc+1)+j)
            buffermpi1(1,bufbegrat2(proc+1)+j) = 
     $          x(iglob)
            buffermpi1(2,bufbegrat2(proc+1)+j) = 
     $          y(iglob)
            buffermpi1(3,bufbegrat2(proc+1)+j) = 
     $          z(iglob)
            buffermpi1(4,bufbegrat2(proc+1)+j) = 
     $          v(1,iglob)
            buffermpi1(5,bufbegrat2(proc+1)+j) = 
     $          v(2,iglob)
            buffermpi1(6,bufbegrat2(proc+1)+j) = 
     $          v(3,iglob)
            buffermpi1(7,bufbegrat2(proc+1)+j) = 
     $          xold(iglob)
            buffermpi1(8,bufbegrat2(proc+1)+j) = 
     $          yold(iglob)
            buffermpi1(9,bufbegrat2(proc+1)+j) = 
     $          zold(iglob)
          end do
        end if
      end do
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_ISEND(buffermpi1(1,bufbegrat2(proc+1)),
     $     9*buflenrat2(proc+1),MPI_REAL8,proc,tag,COMM_BEAD,
     $     reqsend(i),ierr)
        end if
      end do
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
      do i = 1, nneig_send
        proc = pneig_send(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_WAIT(reqrec(i),status,ierr)
          do j = 0, buflenrat1(proc+1)-1
            iglob = bufrat1(bufbegrat1(proc+1)+j)
            x(iglob) = buffermpi2(1,bufbegrat1(proc+1)+j)
            y(iglob) = buffermpi2(2,bufbegrat1(proc+1)+j)
            z(iglob) = buffermpi2(3,bufbegrat1(proc+1)+j)
            v(1,iglob) = buffermpi2(4,bufbegrat1(proc+1)+j)
            v(2,iglob) = buffermpi2(5,bufbegrat1(proc+1)+j)
            v(3,iglob) = buffermpi2(6,bufbegrat1(proc+1)+j)
            xold(iglob) = buffermpi2(7,bufbegrat1(proc+1)+j)
            yold(iglob) = buffermpi2(8,bufbegrat1(proc+1)+j)
            zold(iglob) = buffermpi2(9,bufbegrat1(proc+1)+j)
          end do
        end if
      end do
c
      deallocate (buffermpi2)
      deallocate (buffermpi1)
      deallocate (reqsend)
      deallocate (reqrec)
      return
      end
