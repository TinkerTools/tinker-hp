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
#include "tinker_precision.h"
      subroutine rattle (dt)
      use atmlst
      use atmtyp
      use atomsMirror
      use domdec
      use freeze
      use molcul
      use inform
      use iounit
      use moldyn
      use tinheader ,only: ti_p,re_p
      use usage
      use mpi
      implicit none
      integer i,j,k,iloc,iglob,ierr,imol_
      integer ia,ib,ialoc,ibloc,mode
      integer niter,maxiter
      real(r_p) dt,eps
      real(t_p) sor
      real(t_p) xr,yr,zr
      real(t_p) xo,yo,zo
      real(t_p) dot,rma,rmb
      real(t_p) weigh,dist2
      real(t_p) delta,term
      real(t_p) xterm,yterm,zterm
      real(t_p), allocatable :: displace(:,:)
      integer :: n_not_done,ii,nratmol_
      logical done
      logical, allocatable :: moved(:)
      logical, allocatable :: update(:)
c
c#ifdef _OPENACC
c 15   format(' Rattle Feature is for now unavailable !!!'
c     &    ,/,'   Please Target host build(CPU) to benefit from it')
c      write(0,15)
c      call fatal
c#endif

      if (nproc > 1) then
         write (0,*) 'Error: RATTLE is not implemented for for nproc>1'
         call fatal
      end if

c
c     initialize arrays for MPI 
c
c      call initmpirattle
c      call commrattleinit
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      sor = 1.25_ti_p
      eps = rateps
c
c     apply RATTLE to distances and half-step velocity values
c
      niter = 0
      done = .false.
      !do while (.not.done .and. niter.lt.maxiter)
      !   niter = niter + 1
      !   n_not_done = 0
      !   ! do iloc = 1, nratloc
      !   !    i = ratglob(iloc) 

!$acc parallel loop default(present) async
      do imol_ = 1, nmol
         if (nratmol(imol_) == 0) cycle
!$acc loop seq
         do niter = 1,maxiter
           n_not_done=0
!$acc loop seq
           do ii = 1,nratmol(imol_)
            i = iratmol(ii,imol_)
            ia = irat(1,i)
            ib = irat(2,i)
            if(.not.(use(ia) .or. use(ib))) cycle
            
            xr = x(ib) - x(ia)
            yr = y(ib) - y(ia)
            zr = z(ib) - z(ia)
            !if (ratimage(i))  call image (xr,yr,zr)
            dist2 = xr**2 + yr**2 + zr**2
            delta = krat(i)**2 - dist2
            if(abs(delta) <= eps) cycle

            n_not_done = n_not_done + 1
            xo = xold(ib) - xold(ia)
            yo = yold(ib) - yold(ia)
            zo = zold(ib) - zold(ia)
            !if (ratimage(i))  call image (xo,yo,zo)
            dot = xr*xo + yr*yo + zr*zo
            rma = 1.0_ti_p / mass(ia)
            rmb = 1.0_ti_p / mass(ib)
            term = sor * delta / (2.0_ti_p * (rma+rmb) * dot)
            xterm = xo * term
            yterm = yo * term
            zterm = zo * term
            x(ia) =  x(ia) - xterm*rma
            y(ia) =  y(ia) - yterm*rma
            z(ia) =  z(ia) - zterm*rma
            x(ib) =  x(ib) + xterm*rmb
            y(ib) =  y(ib) + yterm*rmb
            z(ib) =  z(ib) + zterm*rmb
            rma = rma / dt
            rmb = rmb / dt
            v(1,ia) = v(1,ia) - xterm*rma
            v(2,ia) = v(2,ia) - yterm*rma
            v(3,ia) = v(3,ia) - zterm*rma
            v(1,ib) = v(1,ib) + xterm*rmb
            v(2,ib) = v(2,ib) + yterm*rmb
            v(3,ib) = v(3,ib) + zterm*rmb
           end do
           if (n_not_done == 0) exit
         end do
      end do
c!$acc wait
c         !write(*,*) 'n_not_done = ',n_not_done
c         done = (n_not_done == 0)
cc         call MPI_ALLREDUCE(MPI_IN_PLACE,done,1,MPI_LOGICAL,MPI_LAND,
cc     $    COMM_TINKER,ierr)
cc
c      end do


      call reCast_position

c      call commrattleend
c
c     write information on the number of iterations needed
c
c      if (niter .eq. maxiter) then
c         write (iout,10)
c   10    format (/,' RATTLE  --  Warning, Distance Constraints',
c     &              ' not Satisfied')
c         call fatal
c      if (debug) then
c         write (iout,20)  niter
c   20    format (' RATTLE   --  Distance Constraints met at',i6,
c     &              ' Iterations')
c      end if
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
      use atomsMirror
      use domdec
      use group
      use freeze
      use inform
      use iounit
      use moldyn
      use tinheader ,only: ti_p,re_p
      use units
      use usage
      use virial
      use molcul
      use mpi
      implicit none
      integer i,j,k,iloc,iglob,ierr
      integer ia,ib,ialoc,ibloc,mode
      integer niter,maxiter
      integer start,stop
      real(r_p) dt,eps
      real(t_p) sor
      real(t_p) xr,yr,zr
      real(t_p) xv,yv,zv
      real(t_p) dot,rma,rmb
      real(t_p) weigh,vterm,term
      real(t_p) xterm,yterm,zterm
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      real(r_p), save :: gxx,gyy,gzz
      real(r_p), save :: gyx,gzx,gzy
      real(r_p) gxx_loc,gyy_loc,gzz_loc
      real(r_p) gyx_loc,gzx_loc,gzy_loc
      logical done
      integer :: n_not_done,imol_,ii
      logical, save :: f_in=.TRUE.
c
c#ifdef _OPENACC
c 15   format(' Rattle Feature is for now unavailable !!!'
c     &    ,/,' > Please Target host(CPU) build to benefit from it')
c      write(0,15)
c      call fatal
c#endif
c
c     initialize arrays for MPI 
c
c      call initmpirattle
c      call commrattleinit
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      niter = 0
      done = .false.
      sor = 1.25_ti_p
      eps = rateps / dt
      vterm = 2.0_ti_p / (dt*convert)
c
c     apply the RATTLE algorithm to correct the velocities
c
c      do while (.not.done .and. niter.lt.maxiter)
c         niter = niter + 1
c         n_not_done = 0

      if (f_in) then
         f_in = .false.
        gxx = 0.0_re_p
        gyy = 0.0_re_p
        gzz = 0.0_re_p
        gyx = 0.0_re_p
        gzx = 0.0_re_p
        gzy = 0.0_re_p
!$acc enter data async copyin(gxx,gyy,gzz,gyx,gzx,gzy)
      endif

!$acc parallel loop default(present)  async
!$acc& present(gxx,gyy,gzz,gyx,gzx,gzy)
!$acc& reduction(+:gxx,gyy,gzz,gyx,gzx,gzy)
      do imol_ = 1, nmol
         if (nratmol(imol_) == 0) cycle
         gxx_loc = 0.0_re_p
         gyy_loc = 0.0_re_p
         gzz_loc = 0.0_re_p
         gyx_loc = 0.0_re_p
         gzx_loc = 0.0_re_p
         gzy_loc = 0.0_re_p
!$acc loop seq
         do niter = 1, maxiter
           n_not_done=0
!$acc loop seq
           do ii = 1,nratmol(imol_)
            i = iratmol(ii,imol_)
            ia = irat(1,i)
            ib = irat(2,i)
            if(.not.(use(ia) .or. use(ib))) cycle

            xr = x(ib) - x(ia)
            yr = y(ib) - y(ia)
            zr = z(ib) - z(ia)
            !if (ratimage(i))  call image (xr,yr,zr)
            xv = v(1,ib) - v(1,ia)
            yv = v(2,ib) - v(2,ia)
            zv = v(3,ib) - v(3,ia)
            dot = xr*xv + yr*yv + zr*zv
            rma = 1.0_ti_p / mass(ia)
            rmb = 1.0_ti_p / mass(ib)
            term = -dot / ((rma+rmb) * krat(i)**2)
            if (abs(term) <= eps) cycle

            n_not_done = n_not_done + 1
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
            gxx_loc = gxx_loc + vxx
            gyy_loc = gyy_loc + vyy
            gzz_loc = gzz_loc + vzz
            gyx_loc = gyx_loc + vyx
            gzx_loc = gzx_loc + vzx
            gzy_loc = gzy_loc + vzy
           end do
           if (n_not_done == 0) exit
         end do

         gxx = gxx + gxx_loc
         gyy = gyy + gyy_loc
         gzz = gzz + gzz_loc
         gyx = gyx + gyx_loc
         gzx = gzx + gzx_loc
         gzy = gzy + gzy_loc
      end do
c!$acc wait
         !write(*,*) 'n_not_done 2 = ',n_not_done
         ! done = (n_not_done == 0)
c         call MPI_ALLREDUCE(MPI_IN_PLACE,done,1,MPI_LOGICAL,MPI_LAND,
c     $    COMM_TINKER,ierr)
c
      ! end do
c      call commrattleend
c      call MPI_ALLREDUCE(MPI_IN_PLACE,virtemp,9,MPI_TPREC,MPI_SUM,
c     $ COMM_TINKER,ierr)

      if (use_virial) then
!$acc serial async present(vir,gxx,gyy,gzz,gyx,gzx,gzy)
         vir(1,1) = vir(1,1) - gxx
         vir(2,1) = vir(2,1) - gyx
         vir(3,1) = vir(3,1) - gzx
         vir(1,2) = vir(1,2) - gyx
         vir(2,2) = vir(2,2) - gyy
         vir(3,2) = vir(3,2) - gzy
         vir(1,3) = vir(1,3) - gzx
         vir(2,3) = vir(2,3) - gzy
         vir(3,3) = vir(3,3) - gzz
         gxx = 0.0_re_p
         gyy = 0.0_re_p
         gzz = 0.0_re_p
         gyx = 0.0_re_p
         gzx = 0.0_re_p
         gzy = 0.0_re_p
!$acc end serial
c!$acc update host(vir) async
      endif




c
c     write information on the number of iterations needed
c
c      if (niter .eq. maxiter) then
c         write (iout,10)
c   10    format (/,' RATTLE2  --  Warning, Velocity Constraints',
c     &              ' not Satisfied')
c         call fatal
c      else if (debug) then
c         write (iout,20)  niter
c   20    format (' RATTLE2  --  Velocity Constraints met at',i6,
c     &              ' Iterations')
c      end if
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
      use atomsMirror
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
     $   pneig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
      do i = 1, nneig_recep
        if (pneig_recep(i).ne.rank) then
          tag = nproc*pneig_recep(i) + rank + 1
          call MPI_ISEND(buflenrat2(pneig_recep(i)+1),1,MPI_INT,
     $     pneig_recep(i),tag,COMM_TINKER,reqsend(i),ierr)
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
     $     MPI_INT,pneig_send(i),tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
      do i = 1, nneig_recep
        if (pneig_recep(i).ne.rank) then
          tag = nproc*pneig_recep(i) + rank + 1
          call MPI_ISEND(bufrat2(bufbegrat2(pneig_recep(i)+1)),
     $     buflenrat2(pneig_recep(i)+1),MPI_INT,pneig_recep(i),tag,
     $     COMM_TINKER,reqsend(i),ierr)
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
      use tinheader ,only: ti_p,re_p
      implicit none
      integer i,j,proc,ierr,iglob,iloc
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqrec(:),reqsend(:)
      real(t_p), allocatable :: buffermpi1(:,:),buffermpi2(:,:)
c
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (buffermpi1(9,max(nbloc,1)))
      buffermpi1 = 0_ti_p
      allocate (buffermpi2(9,max(nbloc,1)))
      buffermpi2 = 0_ti_p
c
c     Begin the reception of the constrained atoms
c
      do i = 1, nneig_recep
        proc = pneig_recep(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_IRECV(buffermpi1(1,bufbegrat2(proc+1)),
     $      9*buflenrat2(proc+1),
     $      MPI_TPREC,proc,tag,COMM_TINKER,reqrec(i),ierr)
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
     $     9*buflenrat1(proc+1),MPI_TPREC,proc,tag,COMM_TINKER,
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
!$acc update device(x(:),y(:),z(:),v,xold,yold,zold)
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
      integer i,j,proc,ierr,iglob
      integer tag,status(MPI_STATUS_SIZE)
      integer, allocatable :: reqrec(:),reqsend(:)
      real(t_p), allocatable :: buffermpi1(:,:),buffermpi2(:,:)
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
     $      MPI_TPREC,proc,tag,COMM_TINKER,reqrec(i),ierr)
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
     $     9*buflenrat2(proc+1),MPI_TPREC,proc,tag,COMM_TINKER,
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
