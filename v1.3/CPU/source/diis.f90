!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###########################################################
!     ##                                                       ##
!     ##  subroutine diis                                      ##
!     ##                                                       ##
!     ###########################################################
!
!> @brief 
!> compute quantities to run diis extrapolation
!> @param[in] ndsmax: max number of iterations to extrapolate on
!> @param[in] n: size of matrix
!> @param[in] xdiis: matrix of previous solutions
!> @param[in] ediis: matrix of previous errors
!> @param[in] bmat: B diis matrix
!> @param[in] nmat: current number of iterations to extrapolate on
!> @param[in] reqdiis: array of request (MPI)
!> @param[in] comm: MPI communicator
subroutine diis(ndsmax,n,xdiis,ediis,bmat,nmat,reqdiis,comm)
   use inform
   use iounit
   implicit none
!
!     driver to perform direct inversion in iterative subsspace
!     extrapolation.
!
!
   integer ndsmax, n, nmat, comm
   real*8  xdiis(n,ndsmax), ediis(n,ndsmax), bmat(ndsmax+1,ndsmax+1)
   integer j, k
   integer reqdiis(*)
!
1000 format(' DIIS will restart the extrapolation.')
1010 format(' Restarting from ',I4)
!
   if (deb_Path) write(iout,*), 'diis '
!
!
!     if needed, resize the matrix and restart the extrapolation.
!
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
!
!     build the diis matrix
!
   call makeb(ndsmax,n,nmat,ediis,bmat,reqdiis,comm)
   nmat = nmat + 1
   return
!
end
!
!> @brief 
!> compute run diis extrapolation
!> @param[in] n: size of matrix
!> @param[in] nmat: number of iterations to extrapolate on
!> @param[in] x: vector of previous solutions
!> @param[in] cex: extrapolation coefficients
!> @param[in] xex: extrapolated solution
subroutine extrap(n,nmat,x,cex,xex)
   implicit none
!
!     perform the DIIS extrapolation.
!
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
!
!> @brief 
!> compute b matrix for diis extrapolation
!> @param[in] ndsmax: max number of iterations to extrapolate on
!> @param[in] n: size of matrix
!> @param[in] nvec: current number of iterations to extrapolate on
!> @param[in] e: matrix of previous errors
!> @param[in] b: B diis matrix
!> @param[in] reqdiis: array of request (MPI)
!> @param[in] comm: MPI communicator
subroutine makeb(ndsmax,n,nvec,e,b,reqdiis,comm)
   use mpi
   implicit none
!
!   build or update pulay's diis matrix.
!
   integer ndsmax, n, nvec, j, ierr, reqdiis(*), comm
   real*8 e(n,*), b(ndsmax+1,*), zero, one, sprod
   data zero/0.d0/, one/1.d0/
!
!
!
   if(nvec.eq.1) then
!
!   initialize:
!
      b(1,1) = zero
      b(1,2) = one
      b(2,1) = one
      if (n.gt.0) b(2,2) = sprod(n,e(1,1),e(1,1))
      call MPI_IALLREDUCE(MPI_IN_PLACE,b(2,2),1,MPI_REAL8,MPI_SUM,&
      &comm,reqdiis(1),ierr)
   else
!
!   update the lagrangian line:
!
      b(nvec+1,1) = one
      b(1,nvec+1) = one
!
!   calculate the new matrix elements:
!
      do 100 j = 1, nvec - 1
         if (n.gt.0) b(nvec+1,j+1) = sprod(n,e(1,j),e(1,nvec))
         if (n.gt.0) b(j+1,nvec+1) = sprod(n,e(1,j),e(1,nvec))
         call MPI_IALLREDUCE(MPI_IN_PLACE,b(nvec+1,j+1),1,MPI_REAL8,&
         &MPI_SUM,comm,reqdiis(2*j-1),ierr)
         call MPI_IALLREDUCE(MPI_IN_PLACE,b(j+1,nvec+1),1,MPI_REAL8,&
         &MPI_SUM,comm,reqdiis(2*j),ierr)
100   continue
      if (n.gt.0) b(nvec+1,nvec+1) = sprod(n,e(1,nvec),e(1,nvec))
      call MPI_IALLREDUCE(MPI_IN_PLACE,b(nvec+1,nvec+1),1,&
      &MPI_REAL8,MPI_SUM,comm,reqdiis(2*nvec-1),ierr)
   endif
!
   return
end
!
subroutine no_com_diis(ndsmax,n,xdiis,ediis&
&,bmat,nmat)
   implicit none
!
!     driver to perform direct inversion in iterative subsspace
!     extrapolation.
!
!
   integer ndsmax, n, nmat
   real*8  xdiis(n,ndsmax), ediis(n,ndsmax), bmat(ndsmax+1,ndsmax+1)
   integer j, k
!
1000 format(' DIIS will restart the extrapolation.')
1010 format(' Restarting from ',I4)
!
!     if needed, resize the matrix and restart the extrapolation.
!
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

!
!     build the diis matrix
!
   call no_com_makeb(ndsmax,n,nmat,ediis,bmat)
   nmat = nmat + 1
   return
!
end
!
subroutine no_com_makeb(ndsmax,n,nvec,e,b)
   use mpi
   implicit none
!
!   build or update pulay's diis matrix.
!
   integer ndsmax, n, nvec, j
   real*8 e(n,*), b(ndsmax+1,*), zero, one, sprod
   data zero/0.d0/, one/1.d0/
!
   if(nvec.eq.1) then
!
!   initialize:
!
      b(1,1) = zero
      b(1,2) = one
      b(2,1) = one
      if (n.gt.0) b(2,2) = sprod(n,e(1,1),e(1,1))
   else
!
!   update the lagrangian line:
!
      b(nvec+1,1) = one
      b(1,nvec+1) = one
!
!   calculate the new matrix elements:
!
      do 100 j = 1, nvec - 1
         if (n.gt.0) b(nvec+1,j+1) = sprod(n,e(1,j),e(1,nvec))
         if (n.gt.0) b(j+1,nvec+1) = sprod(n,e(1,j),e(1,nvec))
100   continue
      if (n.gt.0) b(nvec+1,nvec+1) = sprod(n,e(1,nvec),e(1,nvec))
   endif
   return
end
