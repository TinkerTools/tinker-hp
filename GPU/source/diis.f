c      
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine diis                                      ##
c     ##                                                       ##
c     ###########################################################
c
#include "tinker_precision.h"
      subroutine diis(ndsmax,n,xdiis,ediis,bmat,nmat,reqdiis,comm)
      implicit none
c
c     driver to perform direct inversion in iterative subsspace
c     extrapolation.
c
c
      integer ndsmax, n, comm
      integer nmat
      real(t_p)  xdiis(n,ndsmax), ediis(n,ndsmax)
      real(t_p)  bmat(ndsmax+1,ndsmax+1)
      real(t_p)  bmat_buf(2:nmat-10,2:nmat-10)
      integer reqdiis(*)
      integer j, k, l, jst, jnd, quo, res
c
c1000 format(' DIIS will restart the extrapolation.')
 1010 format(' Restarting from ',I4)
c
c     if needed, resize the matrix and restart the extrapolation.
c
      if(nmat.ge.ndsmax) then
        if (nmat>10) then 
!$acc data create(bmat_buf) present(bmat,xdiis,ediis) async

        quo = (nmat-10)/10

!$acc parallel loop collapse(2) async
        do j = 2, nmat - 10
          do k = 2, nmat - 10
            bmat_buf(j,k) = bmat(j+10,k+10)
          end do
        end do
!$acc parallel loop collapse(2) async
        do j = 2, nmat - 10
          do k = 2, nmat - 10
            bmat(j,k) = bmat_buf(j,k)
          end do
        end do

        do l = 0, quo
          jst = l*10+1
          jnd = (l+1)*10
          if (l.eq.quo) jnd = nmat-10
!$acc parallel loop collapse(2) async
          do j = jst, jnd
            do k = 1, n
               xdiis(k,j) = xdiis(k,j+10)
               ediis(k,j) = ediis(k,j+10)
            end do
            !call amove(n,xdiis(1,j+10),xdiis(1,j))
            !call amove(n,ediis(1,j+10),ediis(1,j))
          end do
        end do

!$acc end data
        end if
        write(6,1010) nmat - 10
        nmat = nmat - 10
      end if
c
c     build the diis matrix
c
      call makeb(ndsmax,n,nmat,ediis,bmat,reqdiis,comm)
      nmat = nmat + 1
      return
c
      end
c
      subroutine extrap(n,nmat,x,cex,xex)
c
c     perform the DIIS extrapolation.
c
      implicit none
      integer n, nmat
      real(t_p)  x(n,*), cex(*), xex(*)
      integer j, k

!$acc parallel loop collapse(2) default(present) async
      do j = 1, nmat
        do k = 1, n
!$acc atomic
          xex(k) = xex(k) + cex(j+1)*x(k,j)
        end do
      end do
      end
c
      subroutine makeb(ndsmax,n,nvec,e,b,reqdiis,comm)
      use mpi
      implicit none 
c
c   build or update pulay's diis matrix. 
c
      integer,parameter::ti_p=t_p
      integer ndsmax, n, nvec, j, ierr, reqdiis(*), comm
      real(t_p) e(n,*), b(ndsmax+1,*), zero, one 
      real(t_p) sprod  ! Function for scalar product written in linalg.f
      parameter ( zero=0.0_ti_p, one=1.0_ti_p )
c
      if(nvec.eq.1) then 
c
c   initialize:
c
!$acc serial present(b) async
        b(1,1) = zero
        b(1,2) = one
        b(2,1) = one
!$acc end serial
        if (n.gt.0) call sprodgpu(n,e(1,1),e(1,1),b(2,2))
!$acc wait
!$acc host_data use_device(b)
        call MPI_ALLREDUCE(MPI_IN_PLACE,b(2,2),1,MPI_TPREC,MPI_SUM,
     $        comm,ierr)
!$acc end host_data
      else
c
c   update the lagrangian line:
c
!$acc serial async
        b(nvec+1,1) = one
        b(1,nvec+1) = one
!$acc end serial
c
c   calculate the new matrix elements:
c
        do j = 1, nvec - 1
          if (n.gt.0) call sprodgpu(n,e(1,j),e(1,nvec),b(nvec+1,j+1))
          if (n.gt.0) call sprodgpu(n,e(1,j),e(1,nvec),b(j+1,nvec+1))
!$acc wait
!$acc host_data use_device(b)
          call MPI_ALLREDUCE(MPI_IN_PLACE,b(nvec+1,j+1),1,MPI_TPREC,
     $      MPI_SUM,comm,ierr)
          call MPI_ALLREDUCE(MPI_IN_PLACE,b(j+1,nvec+1),1,MPI_TPREC,
     $      MPI_SUM,comm,ierr)
!$acc end host_data
        end do
        if(n.gt.0) call sprodgpu(n,e(1,nvec),e(1,nvec),b(nvec+1,nvec+1))
!$acc wait
!$acc host_data use_device(b)
        call MPI_ALLREDUCE(MPI_IN_PLACE,b(nvec+1,nvec+1),1,
     $   MPI_TPREC,MPI_SUM,comm,ierr)
!$acc end host_data

      endif
      end
c
      subroutine no_com_diis(ndsmax,n,xdiis,ediis
     $ ,bmat,nmat)
      implicit none
c
c     driver to perform direct inversion in iterative subsspace
c     extrapolation.
c
c
      integer ndsmax, n, nmat
      real(t_p)  xdiis(n,ndsmax), ediis(n,ndsmax)
      real(t_p)  bmat(ndsmax+1,ndsmax+1)
      real(t_p)  bmat_buf(2:nmat-10,2:nmat-10)
      integer j, k, l, jst, jnd, quo, res
c
c1000 format(' DIIS will restart the extrapolation.')
 1010 format(' Restarting from ',I4)
c
c     if needed, resize the matrix and restart the extrapolation.
c
      if(nmat.ge.ndsmax) then
        if (nmat>10) then 
!$acc data create(bmat_buf) present(bmat,xdiis,ediis) async

        quo = (nmat-10)/10

!$acc parallel loop collapse(2) async
        do j = 2, nmat - 10
          do k = 2, nmat - 10
            bmat_buf(j,k) = bmat(j+10,k+10)
          end do
        end do
!$acc parallel loop collapse(2) async
        do j = 2, nmat - 10
          do k = 2, nmat - 10
            bmat(j,k) = bmat_buf(j,k)
          end do
        end do

        do l = 0, quo
          jst = l*10+1
          jnd = (l+1)*10
          if (l.eq.quo) jnd = nmat-10
!$acc parallel loop collapse(2) async default(present)
          do j = jst, jnd
            do k = 1, n
               xdiis(k,j) = xdiis(k,j+10)
               ediis(k,j) = ediis(k,j+10)
            end do
            !call amove(n,xdiis(1,j+10),xdiis(1,j))
            !call amove(n,ediis(1,j+10),ediis(1,j))
          end do
        end do

!$acc end data
        end if
        write(6,1010) nmat - 10
        nmat = nmat - 10
      end if

c
c     build the diis matrix
c
      call no_com_makeb(ndsmax,n,nmat,ediis,bmat)
      nmat = nmat + 1
      return
c
      end
c
      subroutine no_com_makeb(ndsmax,n,nvec,e,b)
      use mpi
      implicit none 
c
c   build or update pulay's diis matrix. 
c
      integer,parameter::ti_p=t_p
      integer ndsmax, n, nvec, j, ierr
      real(t_p) e(n,*), b(ndsmax+1,*), zero, one
      real(t_p) sprod  ! Function for scalar product written in linalg.f
      parameter ( zero=0.0_ti_p, one=1.0_ti_p )
c
      if(nvec.eq.1) then 
c
c   initialize:
c
!$acc serial async present(b)
        b(1,1) = zero
        b(1,2) = one
        b(2,1) = one
!$acc end serial
        if (n.gt.0) call sprodgpu(n,e(1,1),e(1,1),b(2,2))
      else
c
c   update the lagrangian line:
c
!$acc serial async present(b)
        b(nvec+1,1) = one
        b(1,nvec+1) = one
!$acc end serial
c
c   calculate the new matrix elements:
c
        do j = 1, nvec - 1
          if (n.gt.0) call sprodgpu(n,e(1,j),e(1,nvec),b(nvec+1,j+1))
          if (n.gt.0) call sprodgpu(n,e(1,j),e(1,nvec),b(j+1,nvec+1))
        end do
        if (n.gt.0)call sprodgpu(n,e(1,nvec),e(1,nvec),b(nvec+1,nvec+1))
      endif
      return
      end
