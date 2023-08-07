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
      subroutine diis(ndsmax,n,xdiis,ediis,bmat,nmat,reqdiis,comm)
      implicit none
c
c     driver to perform direct inversion in iterative subsspace
c     extrapolation.
c
c
      integer ndsmax, n, nmat, comm
      real*8  xdiis(n,ndsmax), ediis(n,ndsmax), bmat(ndsmax+1,ndsmax+1)
      integer j, k 
      integer reqdiis(*)
c
 1000 format(' DIIS will restart the extrapolation.')
 1010 format(' Restarting from ',I4)
c
c     if needed, resize the matrix and restart the extrapolation.
c
      if(nmat.ge.ndsmax) then
        do j = 2, nmat - 10
          do k = 2, nmat - 10
            bmat(j,k) = bmat(j+10,k+10)
          end do
        end do
        do j = 1, nmat - 10
          call amove(n,xdiis(1,j+10),xdiis(1,j))
          call amove(n,ediis(1,j+10),ediis(1,j))
        end do
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
      implicit none
c
c     perform the DIIS extrapolation.
c
      integer n, nmat
      real*8  x(n,*), cex(*), xex(*)
      integer j, k
      do j = 1, nmat
        do k = 1, n
          xex(k) = xex(k) + cex(j+1)*x(k,j)
        end do
      end do
      return
      end
c
      subroutine makeb(ndsmax,n,nvec,e,b,reqdiis,comm)
      use mpi
      implicit none 
c
c   build or update pulay's diis matrix. 
c
      integer ndsmax, n, nvec, j, i, ierr, reqdiis(*), comm
      real*8 e(n,*), b(ndsmax+1,*), zero, one, sprod
      data zero/0.d0/, one/1.d0/
c
c
c
      if(nvec.eq.1) then 
c
c   initialize:
c
        b(1,1) = zero
        b(1,2) = one
        b(2,1) = one
        if (n.gt.0) b(2,2) = sprod(n,e(1,1),e(1,1))
        call MPI_IALLREDUCE(MPI_IN_PLACE,b(2,2),1,MPI_REAL8,MPI_SUM,
     $        comm,reqdiis(1),ierr)
      else
c
c   update the lagrangian line:
c
        b(nvec+1,1) = one
        b(1,nvec+1) = one
c
c   calculate the new matrix elements:
c
        do 100 j = 1, nvec - 1
          if (n.gt.0) b(nvec+1,j+1) = sprod(n,e(1,j),e(1,nvec))
          if (n.gt.0) b(j+1,nvec+1) = sprod(n,e(1,j),e(1,nvec))
          call MPI_IALLREDUCE(MPI_IN_PLACE,b(nvec+1,j+1),1,MPI_REAL8,
     $      MPI_SUM,comm,reqdiis(2*j-1),ierr)
          call MPI_IALLREDUCE(MPI_IN_PLACE,b(j+1,nvec+1),1,MPI_REAL8,
     $      MPI_SUM,comm,reqdiis(2*j),ierr)
 100      continue
        if (n.gt.0) b(nvec+1,nvec+1) = sprod(n,e(1,nvec),e(1,nvec))
        call MPI_IALLREDUCE(MPI_IN_PLACE,b(nvec+1,nvec+1),1,
     $   MPI_REAL8,MPI_SUM,comm,reqdiis(2*nvec-1),ierr)
      endif
c
      return
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
      real*8  xdiis(n,ndsmax), ediis(n,ndsmax), bmat(ndsmax+1,ndsmax+1)
      integer j, k 
c
 1000 format(' DIIS will restart the extrapolation.')
 1010 format(' Restarting from ',I4)
c
c     if needed, resize the matrix and restart the extrapolation.
c
      if(nmat.ge.ndsmax) then
        do j = 2, nmat - 10
          do k = 2, nmat - 10
            bmat(j,k) = bmat(j+10,k+10)
          end do
        end do
        do j = 1, nmat - 10
          call amove(n,xdiis(1,j+10),xdiis(1,j))
          call amove(n,ediis(1,j+10),ediis(1,j))
        end do
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
      integer ndsmax, n, nvec, j, ierr
      real*8 e(n,*), b(ndsmax+1,*), zero, one, sprod
      data zero/0.d0/, one/1.d0/
c
      if(nvec.eq.1) then 
c
c   initialize:
c
        b(1,1) = zero
        b(1,2) = one
        b(2,1) = one
        if (n.gt.0) b(2,2) = sprod(n,e(1,1),e(1,1))
      else
c
c   update the lagrangian line:
c
        b(nvec+1,1) = one
        b(1,nvec+1) = one
c
c   calculate the new matrix elements:
c
        do 100 j = 1, nvec - 1
          if (n.gt.0) b(nvec+1,j+1) = sprod(n,e(1,j),e(1,nvec))
          if (n.gt.0) b(j+1,nvec+1) = sprod(n,e(1,j),e(1,nvec))
 100      continue
        if (n.gt.0) b(nvec+1,nvec+1) = sprod(n,e(1,nvec),e(1,nvec))
      endif
      return
      end
