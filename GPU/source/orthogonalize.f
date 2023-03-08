c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_precision.h"
#if TINKER_DOUBLE_PREC
#define M_sbev dsbev
#define M_cublasGemm cublasDgemm
#elif TINKER_MIXED_PREC
#define M_sbev dsbev
#define M_cublasGemm cublasSgemm
#else
#define M_sbev ssbev
#define M_cublasGemm cublasSgemm
#endif
      
      ! Compute the value of the characteristical polynome of a
      ! tri-diagonal sub-matrix stored in Diag and SubDiag
      real(r_p) function pol_c(Diag, SubDiag, i, lambda) result(pol_val)
      implicit none
      real(r_p),intent(in) :: Diag(*),SubDiag(*)
      integer  ,intent(in) :: i
      real(r_p),intent(in) :: lambda
      integer k
      real(r_p) pol_s,new_val

      if (i.eq.0) then
         pol_val = 1.0
      else if (i.eq.1) then
         pol_val = Diag(1) - lambda
      else
         pol_s   = 1.0
         pol_val = Diag(1) - lambda
         do k = 2, i
            new_val = (Diag(k)-lambda)*pol_val - SubDiag(k-1)**2*pol_s
            pol_s   = pol_val
            pol_val = new_val
         end do
      end if

      end function

      real(r_p) function pol_cn(Diag,SubDiag,lambda,p0,p1)
     &          result(pval)
      implicit none
      real(r_p),intent(in):: Diag,SubDiag
      real(r_p),intent(in):: lambda
      real(r_p),intent(in):: p0,p1

      pval = (Diag-lambda)*p1 - (SubDiag**2)*p0
      end function

      integer function Nracines(Diag, SubDiag, i, mu)
      implicit none
      real(r_p),intent(in) :: Diag(*),SubDiag(*)
      integer  ,intent(in) :: i
      real(r_p),intent(in) :: mu
      integer,parameter::pos=31
      integer k
      real(8) p0, p1, p2
      real(8),parameter::zero=0.0
      logical tinker_isnan_m
      logical new_l,old_l

      Nracines = 0
      p0       = 1
      p1       = Diag(1) - mu

      if (p1.lt.0) Nracines=1
      old_l    = merge(.false.,.true.,p1.lt.zero)

      do k = 2, i
         p2 = (Diag(k)-mu)*p1 - SubDiag(k-1)**2*p0
         p0 = p1
         p1 = p2
         new_l = merge(.false.,.true.,p2.lt.zero)
         if (p2.eq.zero) new_l = old_l
         if (new_l.neqv.old_l) Nracines = Nracines + 1
c        if (i.eq.1450.and.mu.eq.3.125d0) then
c41         format(I6,d20.6,2F11.6,I6)
c           print 41, k,p2,Diag(k),SubDiag(k-1)**2,Nracines
c        end if
         old_l = new_l

         ! Verify sequence stability and compute an estimation
         if (abs(p2).gt.1.0d300.or.abs(p2).lt.1.0d-300) then
            Nracines = i*Nracines/k
            exit
         end if
      end do

      if (tinker_isnan_m(p2).and.Nracines.ne.0) Nracines=i

      end function

      ! Compute Approximate EigenValues of Tri-Diagonal Matrix following
      ! the Givens bissection method
      ! (Use a dichotomic approach)
      subroutine bissectionGivens1(Diag,SubDiag,n,i,ri,si,tol, Eig)
      use inform
      implicit none
      real(r_p),intent(in) :: Diag(n),SubDiag(n-1)
      integer  ,intent(in) :: i,n
      real(r_p),intent(in) :: ri,si,tol
      real(r_p),intent(out):: Eig(*)
      integer   k,ii,it,stride
      real(r_p) r,s,d,norm
      integer Nracines
      integer Nrac

      if (deb_path) write(*,*) 'bissectionGivens1'
      stride = n/(i-1)

      do k = 0, i-1
         r    = min(ri,si)
         s    = max(ri,si)
         d    = (ri+si)/2.0
         norm = s-r
         ii   = k*stride
         it   = 0
         !if (k.ne.0) ii = n-i+k
         if (k.eq.i-1) ii = n-1
         do while( norm.gt.tol )
            it   = it+1
            Nrac = Nracines(Diag,SubDiag,n,d)
            if (Nrac.ge.ii+1) then
               s = d
            else
               r = d
            end if
            d    = (r+s)/2.0
            norm = s-r

c12         format(A,I3,2I6,3F16.7)
c           write(*,12) 'bissG',it,ii,Nrac,r,s,norm
         end do
         Eig(k+1) = (r+s)/2.0
      end do

      end subroutine

      subroutine MatxVec(Rin,Rout)
      use atmlst
      use domdec
      use ewald
      use iounit
      use inform    ,only: deb_Path,abort,minmaxone
      use interfaces,only: tmatxb_pmegpu,tmatxb_p
      use math
      use mpole
      use orthogonalize,only: MdiagITmat,lz_AQ,lz_AQ1,l_idx
     &                 ,EigenVec_l,saveResAmat_l
      use polar
      use polar_temp,only: murec,dipfield,dipfieldbis,diag
      use polpot
      use potent    ,only: use_polarshortreal
      use sizes
      use timestat
      use units
      use utils
      use utilcomm  ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
      use utilgpu
      implicit none
      integer i,j,k,nrhs
      real(t_p),contiguous:: Rin(:,:,:),Rout(:,:,:)
      integer  ,dimension(nproc):: reqsend,reqrec,req2send,req2rec
      real(t_p) term
      parameter(nrhs=2)

      ! Extension communications
      if (use_polarshortreal) then
         if (.not.skpPcomm) then
         call commdirdirshort(nrhs,1,Rin,reqrec,reqsend)
         call commdirdirshort(nrhs,0,Rin,reqrec,reqsend)
         call commdirdirshort(nrhs,2,Rin,reqrec,reqsend)
         end if
      else
         call commdirdirgpu(nrhs,0,Rin,reqrec,reqsend)
         call commdirdirgpu(nrhs,1,Rin,reqrec,reqsend)
         call commdirdirgpu(nrhs,2,Rin,reqrec,reqsend)
         call commrecdirdipgpu(nrhs,0,murec,Rin,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
         call commrecdirdipgpu(nrhs,1,murec,Rin,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
         call commrecdirdipgpu(nrhs,2,murec,Rin,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
      end if

      ! Compute matrix Vector product
c
      call timer_enter( timer_realdip )
      call tmatxb_p(2,.true.,Rin,Rout)
      call timer_exit ( timer_realdip,quiet_timers )
c
      if (.not.use_polarshortreal) then
      call timer_enter( timer_recdip )
      call tmatxbrecipgpu(Rin,murec,nrhs,dipfield,dipfieldbis)
      call timer_exit ( timer_recdip,quiet_timers )
      end if

      ! Completion communications
      if (use_polarshortreal) then
         call commfieldshort(nrhs,Rout)
      else
         call commfieldgpu(nrhs,Rout)
         call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,req2send,req2rec)
         call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,req2send,req2rec)
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,req2send,req2rec)
      end if

      ! Store result in tall matrix
      if (.not.use_polarshortreal) then
         term = (4.0_ti_p/3.0_ti_p)*aewald**3 / sqrtpi
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,nrhs; do i=1,3
            Rout(i,j,k) = Rout(i,j,k) + dipfield(i,j,k)
     &                  - term*Rin(i,j,k)
         end do; end do; end do
      end if

      if (EigenVec_l.and.saveResAmat_l) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,nrhs; do i=1,3
            if (btest(j,0)) then; lz_AQ(i,k,l_idx) = Rout(i,j,k)
            else;                lz_AQ1(i,k,l_idx) = Rout(i,j,k)
            end if
         end do; end do; end do
      end if

      ! Apply diag preconditionner
      if (MdiagITmat.ne.0) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,nrhs; do i=1,3
            Rout(i,j,k) = diag(k)*Rout(i,j,k)
         end do; end do; end do
      end if

      end subroutine

      subroutine InnerProd( veca, vecb, sza, Ip )
      use domdec
      use mpi
      implicit none
      integer  ,intent(in) :: sza
      real(t_p),intent(in) :: veca(*),vecb(*)
      real(r_p),intent(out):: Ip
      integer i,j,k

      Ip = 0
!$acc parallel loop async default(present)
      do i = 1,sza
         Ip = Ip + veca(i)*vecb(i)
      end do

      if (nproc.gt.1) then
!$acc wait
         call MPI_ALLREDUCE(MPI_IN_PLACE,Ip,1,MPI_RPREC,MPI_SUM
     &                     ,COMM_TINKER,i)
      end if
      end subroutine

      subroutine InnerProdInl( veca, vecb, sza, Ip )
      use domdec
      use mpi
      implicit none
      integer  ,intent(in) :: sza
      real(t_p),intent(in) :: veca(*),vecb(*)
      real(r_p),intent(out):: Ip
      integer i
!$acc routine vector

      Ip = 0
!$acc loop vector
      do i = 1,sza
         Ip = Ip + veca(i)*vecb(i)
      end do

      end subroutine

      subroutine InnerMProdInl( veca, vecb, diag, sza, Ip )
      use domdec
      use mpi
      implicit none
      integer  ,intent(in) :: sza
      real(t_p),intent(in) :: veca(*),vecb(*),diag(*)
      real(r_p),intent(out):: Ip
      integer i
!$acc routine vector

      Ip = 0
!$acc loop vector
      do i = 1,sza
         Ip = Ip + veca(i)*vecb(i)/diag((i-1)/3+1)
      end do

      end subroutine

      subroutine InnerProdEf( veca, vecb, sza, Ip1, Ip2 )
      use domdec
      use mpi
      use orthogonalize,only: debg
      use utilgpu ,only: rec_queue
      implicit none
      integer  dim3,nrhs
      parameter(dim3=3,nrhs=2)
      integer  ,intent(in) :: sza
      real(t_p),intent(in) :: veca(dim3,nrhs,*),vecb(dim3,nrhs,*)
      real(r_p),intent(out):: Ip1,Ip2
      real(r_p) Ip(nrhs)
      integer i,j,k

      Ip1= 0; Ip2= 0
!$acc parallel loop async(rec_queue) collapse(3) default(present)
      do k=1,sza; do j=1,nrhs; do i=1,dim3
         if (btest(j,0)) then
            Ip1 = Ip1 + veca(i,j,k)*vecb(i,j,k)
         else
            Ip2 = Ip2 + veca(i,j,k)*vecb(i,j,k)
         end if
      end do; end do; end do

      if (nproc.gt.1) then
!$acc wait
         Ip(1)=Ip1; Ip(2)=Ip2
         call MPI_ALLREDUCE(MPI_IN_PLACE,Ip,2,MPI_RPREC,MPI_SUM
     &                     ,COMM_TINKER,i)
         Ip1=Ip(1); Ip2=Ip(2)
      end if

      end subroutine

      subroutine InnerMProdEf( veca, vecb, sza, Ip1, Ip2 )
      use domdec
      use mpi
      use mpole      ,only: npoleloc
      use orthogonalize,only: debg
      use polar_temp ,only: diag
      use polar      ,only: tinypol
      use sizes      ,only: tinkerdebug
      implicit none
      integer  dim3,nrhs
      parameter(dim3=3,nrhs=2)
      integer  ,intent(in) :: sza
      real(t_p),intent(in) :: veca(dim3,nrhs,*),vecb(dim3,nrhs,*)
      real(r_p),intent(out):: Ip1,Ip2
      real(t_p) tmp,hugepole
      real(r_p) Ip(nrhs)
      integer i,j,k
      parameter(hugepole=1d5)

      if (tinkerdebug.gt.0.and.sza.gt.npoleloc) then
         print*, 'InnerMProdEf Issue  !! ',
     &           'size > npoleloc ', npoleloc,sza
      end if

      Ip1= 0; Ip2= 0
!$acc parallel loop async collapse(3) default(present)
      do k=1,sza; do j=1,nrhs; do i=1,dim3
         if (btest(j,0)) then
            Ip1 = Ip1 + veca(i,j,k)*vecb(i,j,k)/diag(k)
         else
            Ip2 = Ip2 + veca(i,j,k)*vecb(i,j,k)/diag(k)
         end if
      end do; end do; end do

      if (nproc.gt.1) then
!$acc wait
         Ip(1)=Ip1; Ip(2)=Ip2
         call MPI_ALLREDUCE(MPI_IN_PLACE,Ip,2,MPI_RPREC,MPI_SUM
     &                     ,COMM_TINKER,i)
         Ip1=Ip(1); Ip2=Ip(2)
      end if

      end subroutine

      subroutine lanzcos_init(k0,compEv,mu_)
      use domdec
      use orthogonalize
      use mpole
      use mpi
#ifdef _OPENACC
      use interfaces ,only: initcuSolver
#endif
      use inform     ,only: minmaxone
      use polar_temp
      use tinMemory
      use utilgpu
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
      implicit none
      integer  ,parameter :: nrhs=2
      integer  ,intent(in):: k0,compEv
      real(t_p),intent(in):: mu_(3,nrhs,*)
      real(r_p) l2_norm,l2_norm1
      real(t_p) temp
      integer i,j,k,ierr,sz
      integer(8) sz_
      logical EigV

      maxKD     = k0
      krylovDim = k0
      nEigenV   = compEv
      EigenVec_l= merge(.true.,.false.,nEigenV.gt.0)

      if (nEigenV.gt.krylovDim) then
 61      format(' ERROR ! --- Krylov Space dimension is lower than the'
     &   ,' Eigen vector number ---',/,2I4 ) 
         write(0,61) krylovDim, nEigenV
      end if

      i = 3*k0 + nEigenV**2
      call prmem_requestm(lz_SS_,nrhs*i)
      call prmem_request (lz_WorkSpace,12*npolebloc+2*maxKD)
      call prmem_request (res,3,nrhs,npoleloc)
      call prmem_request (dipfield,3,nrhs,max(npoleloc,1))
      call prmem_request (dipfieldbis,3,nrhs,max(npolerecloc,1))

      if (nproc.gt.1) then
         call prmem_request(buffermpi1  ,3,nrhs,max(npoleloc,1))
         call prmem_request(buffermpimu1,3,nrhs,max(npoleloc,1))
         call prmem_request(buffermpi2  ,3,nrhs,max(npolerecloc,1))
         call prmem_request(buffermpimu2,3,nrhs,max(npolerecloc,1))
      end if

      call Tinker_shellEnv("ORTH",OrthBase,0)
      call Tinker_shellEnv("MINV_TMAT",MdiagITmat,0)

      if(debg.and.rank.eq.0) then
         if (OrthBase.ne.0) print*,' < ----  Enable Base'
     &       ,' Orthoginalisation ---> '
         if (MdiagITmat.ne.0)
     &       print*, ' < ---- Enable M^{-1} * A ----  > '
      end if

      ! lz_SS   [ Vp  | Vp1 |  lz_Q  |  lz_Q1  | lz_AQ | lz_AQ1 | lz_Z | lz_Z1 | MatE_ | MatE1_ ]
      !                     | TVp | TVp1 |

      ! lz_SS_  [ lz_T | lz_T1 | Eigval | EigVal1 |  MatE | MatE1 ]

      if (nEigenV.gt.0) then
         i = 3*npoleloc*maxKD
         j = 3*npoleloc*compEv
         sz = nrhs*(2*    i   +2*j+maxKD**2+nEigenV**2)
        sz_ = nrhs*(2*int(i,8)+2*j+maxKD**2+nEigenV**2)
         if (sz_.gt.int(ishft(1,30),8)) then
 32         format(__FILE__,I4,A,I0,/,A)
            write(0,*) __LINE__,' maximum memory size reached!! '
     &      ,sz_,'  Update your allocator' 
            call fatal
         end if

         call prmem_request(lz_SS,sz)
         call prmem_request(Evec_,nEigenV,nrhs)
         call prmem_request(MIdk,nEigenV,nEigenV)
         call prmem_request(Ipivot,maxKD)
         call prmem_requestm(Evec,nEigenV,nrhs)

           Vp(1:3,1:npoleloc,1:compEv)=>lz_SS(  1:  j); k=j
          Vp1(1:3,1:npoleloc,1:compEv)=>lz_SS(k+1:k+j); k=2*j
          TVp(1:3,1:npoleloc,1:compEv)=>lz_SS(k+1:k+j); k=3*j
         TVp1(1:3,1:npoleloc,1:compEv)=>lz_SS(k+1:k+j); k=4*j
         lz_Q(1:3,1:npoleloc,1:maxKD) =>lz_SS(k+1:k+i); k=k+i
        lz_Q1(1:3,1:npoleloc,1:maxKD) =>lz_SS(k+1:k+i); k=k+i
        lz_AQ(1:3,1:npoleloc,1:maxKD) =>lz_SS(k+1:k+i); k=k+i
       lz_AQ1(1:3,1:npoleloc,1:maxKD) =>lz_SS(k+1:k+i); k=k+i
      lz_Z(1:maxKD,1:maxKD)=>lz_SS(k+1:k+maxKD**2);  k=k+maxKD**2
      lz_Z1(1:maxKD,1:maxKD)=>lz_SS(k+1:k+maxKD**2); k=k+maxKD**2
        MatE_(1:nEigenV,1:nEigenV)=>lz_SS(k+1:k+nEigenV**2)
         k=k+nEigenV**2
       MatE1_(1:nEigenV,1:nEigenV)=>lz_SS(k+1:k+nEigenV**2)


         k = nrhs*3*k0
         MatE (1:nEigenV,1:nEigenV)=> lz_SS_(k+1:k+nEigenV**2)
         k = k + nEigenV**2
         MatE1(1:nEigenV,1:nEigenV)=> lz_SS_(k+1:k+nEigenV**2)
#ifdef _OPENACC
         call init_cublas_handle
         call initcuSolver(rec_stream)
#endif
      end if

      lz_V(1:3,1:2,1:npolebloc,1:2) => lz_WorkSpace(1:12*npolebloc)
       lz_T(1:2*k0) => lz_SS_(1:2*k0);     k =     2*k0
      lz_T1(1:2*k0) => lz_SS_(k+1:k+2*k0); k = k + 2*k0
       EigVal(1:k0) => lz_SS_(k+1:k+k0);   k = k +   k0
      EigVal1(1:k0) => lz_SS_(k+1:k+k0)

      do i = 1,2*k0
         lz_T (i) = 0
         lz_T1(i) = 0
      end do

      if (polgsf.ne.0) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k = 1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
               lz_V(i,j,k,1) = mu_(i,j,k)*diag(k)
            else
               lz_V(i,j,k,1) = mu_(i,j,k)*diag(k)
            end if
         end do; end do; end do
      else
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k = 1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
               lz_V(i,j,k,1) = mu_(i,j,k)
            else
               lz_V(i,j,k,1) = mu_(i,j,k)
            end if
         end do; end do; end do
      end if

      if (MdiagITmat.ne.0) then
         call InnerMProdEf(lz_V,lz_V,npoleloc,l2_norm,l2_norm1)
      else
         call InnerProdEf(lz_V,lz_V,npoleloc,l2_norm,l2_norm1)
      end if
!$acc wait
      l2_norm  = real(1,r_p)/sqrt(l2_norm);
      l2_norm1 = real(1,r_p)/sqrt(l2_norm1);

      ! Normalisation
      if (EigenVec_l) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k = 1,npoleloc; do j=1,2; do i=1,3
         if (btest(j,0)) then
         temp         = lz_V(i,j,k,1)*real(l2_norm,t_p)
         lz_Q (i,k,1) = temp
         else
         temp         = lz_V(i,j,k,1)*real(l2_norm1,t_p)
         lz_Q1(i,k,1) = temp
         end if
         lz_V(i,j,k,1) = temp
      end do; end do; end do
      else
!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k = 1,npoleloc; do j=1,2; do i=1,3
         if (btest(j,0)) then
            lz_V(i,j,k,1) = lz_V(i,j,k,1)*real(l2_norm,t_p)
         else
            lz_V(i,j,k,1) = lz_V(i,j,k,1)*real(l2_norm1,t_p)
         end if
      end do; end do; end do
      end if

      end subroutine

      subroutine lanzcos_iteration
      use domdec
      use inform  ,only: deb_path,minmaxone
      use mpi
      use mpole
      use orthogonalize
      use polar_temp   ,only: res,diag
      use utilgpu ,only: rec_queue,Tinker_shellEnv
      use sizes   ,only: tinkerdebug
      implicit none
      integer  i,j,k,kk,it,ierr
      integer lossOrthIt,notho
      real(t_p) temp
     &         ,nrmbi,nrmbi1,alphaa,alphaa1,betaa,betaa1
      real(r_p) alpha,alpha1,beta,beta1,beta_s,beta1_s
     &         ,nrm_b,nrm_b1,db_scal,db_scal1
      real(r_p) alpha_v(2),beta_v(2)
      real(r_p),parameter:: tol=1d-12

#if TINKER_DOUBLE_PREC
      real(t_p),parameter:: eps=1d-6
#else
      real(t_p),parameter:: eps=1e-4
#endif
      real(t_p),pointer:: scalV(:)
      character*8 cfg
      parameter(lossOrthIt=100)

      interface
      subroutine MatxVec(Rin,Rout)
      real(t_p),contiguous:: Rin(:,:,:),Rout(:,:,:)
      end subroutine
      end interface

      if (deb_path) write(*,*) 'lanzcos_iteration'
      scalV(1:krylovDim)=>
     &         lz_WorkSpace(12*npolebloc+1:12*npolebloc+krylovDim)
      notho = 0
      beta  = 0; beta1  = 0
      l_idx = 1
      saveResAmat_l = .true.

      call MatxVec(lz_V(:,:,:,1),res)

      if (MdiagITmat.ne.0) then
         call InnerMProdEf(lz_V(1:,1:,1:,1),res,npoleloc,alpha,alpha1)
      else
         call InnerProdEf (lz_V(1:,1:,1:,1),res,npoleloc,alpha,alpha1)
      end if
!$acc wait
      alphaa = alpha; alphaa1=alpha1;

      if (MdiagITmat.ne.0) then
!$acc    parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j = 1,2; do i = 1,3
            if (btest(j,0)) then
               temp  = res(i,j,k) - alphaa*lz_V(i,j,k,1)
               res(i,j,k) = temp
               beta  = beta + temp**2/diag(k)
            else
               temp  = res(i,j,k) - alphaa1*lz_V(i,j,k,1)
               res(i,j,k) = temp
               beta1 = beta1 + temp**2/diag(k)
            end if
         end do; end do; end do
      else
!$acc    parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j = 1,2; do i = 1,3
            if (btest(j,0)) then
               temp  = res(i,j,k) - alphaa*lz_V(i,j,k,1)
               res(i,j,k) = temp
               beta  = beta + temp**2
            else
               temp  = res(i,j,k) - alphaa1*lz_V(i,j,k,1)
               res(i,j,k) = temp
               beta1 = beta1 + temp**2
            end if
         end do; end do; end do
      end if
!$acc wait
      if (nproc.gt.1) then
         beta_v(1)=beta; beta_v(2)=beta1
         call MPI_ALLREDUCE(MPI_IN_PLACE,beta_v,2,MPI_RPREC,MPI_SUM,
     &        COMM_TINKER,ierr)
         beta=beta_v(1); beta1=beta_v(2)
      end if

      beta  = sqrt(beta); beta1=sqrt(beta1)
      nrmbi = real(1,r_p)/beta
      nrmbi1= real(1,r_p)/beta1

      if (tinkerdebug.gt.0) then
         if(MdiagITmat.ne.0) then
            call InnerMProdEf(lz_V(1:,1:,1:,1),res,npoleloc,db_scal
     &                       ,db_scal1)
         else
            call InnerProdEf(lz_V(1:,1:,1:,1),res,npoleloc,db_scal
     &                      ,db_scal1)
         end if
!$acc wait
         if (rank.eq.0) then
 14      format(A,I5,2F12.7,2D16.7)
         print 14, 'it lanzcos',1,alpha,alpha1,db_scal,beta
         end if
      end if

      lz_T (1) = alpha
      lz_T1(1) = alpha1

      do kk = 2, krylovDim

         if (beta.lt.tol) exit
         lz_T(maxKD+kk-1) = beta
        lz_T1(maxKD+kk-1) = beta1
         l_idx = kk

         if (EigenVec_l) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
               lz_V(i,j,k,2) = lz_V(i,j,k,1)
               temp          = res(i,j,k)*nrmbi
               lz_V(i,j,k,1) = temp
                lz_Q(i,k,kk) = temp
            else
               lz_V(i,j,k,2) = lz_V(i,j,k,1)
               temp          = res(i,j,k)*nrmbi1
               lz_V(i,j,k,1) = temp
               lz_Q1(i,k,kk) = temp
            end if
         end do; end do; end do
         else
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
               lz_V(i,j,k,2) = lz_V(i,j,k,1)
               lz_V(i,j,k,1) = res(i,j,k)*nrmbi
            else
               lz_V(i,j,k,2) = lz_V(i,j,k,1)
               lz_V(i,j,k,1) = res(i,j,k)*nrmbi1
            end if
         end do; end do; end do
         end if

         if (tinkerdebug.gt.0) then  ! (DEBUG INFO)

         if(MdiagITmat.ne.0) then
         call InnerMProdEf(lz_V(1:,1:,1:,2),lz_V,npoleloc,db_scal
     &                    ,db_scal1)
         else
         call InnerProdEf(lz_V(1:,1:,1:,2),lz_V,npoleloc,db_scal
     &                   ,db_scal1)
         end if

c        if (EigenVec_l.and.kk.gt.lossOrthIt) then
c           !cfg=merge('verb    ','lz_Q    ',kk.gt.286.and.kk.lt.291)
c           cfg='lz_Q'
c           call chkOrthogonality(lz_Q(1,1,1),3*npoleloc,kk
c    &                           ,3*npoleloc,cfg,'lz_Q    ')
c        end if
         end if  ! (DEBUG INFO)

         call MatxVec(lz_V(:,:,:,1),res)

         if (MdiagITmat.ne.0) then
            call InnerMProdEf(lz_V(1:,1:,1:,1),res,npoleloc,alpha
     &                       ,alpha1)
         else
            call InnerProdEf (lz_V(1:,1:,1:,1),res,npoleloc,alpha
     &                       ,alpha1)
         end if
!$acc wait

         lz_T(kk) = alpha
         lz_T1(kk) = alpha1
         !beta_s=beta; beta1_s=beta1
         betaa =beta; betaa1 =beta1
         beta  = 0;   beta1  = 0

         alphaa= alpha; alphaa1= alpha1

         if (MdiagITmat.ne.0) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
            do k=1,npoleloc; do j=1,2; do i=1,3
               if (btest(j,0)) then
               temp = res(i,j,k) - alphaa*lz_V(i,j,k,1)
     &                           - betaa *lz_V(i,j,k,2)
               res(i,j,k) = temp
               beta  = beta + temp**2/diag(k)
               else
               temp  = res(i,j,k) - alphaa1*lz_V(i,j,k,1)
     &                            - betaa1 *lz_V(i,j,k,2)
               res(i,j,k) = temp
               beta1 = beta1 + temp**2/diag(k)
               end if
            end do; end do; end do
         else
!$acc parallel loop collapse(3) default(present) async(rec_queue)
            do k=1,npoleloc; do j=1,2; do i=1,3
               if (btest(j,0)) then
               temp = res(i,j,k) - alphaa*lz_V(i,j,k,1)
     &                           - betaa *lz_V(i,j,k,2)
               res(i,j,k) = temp
               beta  = beta + temp**2
               else
               temp  = res(i,j,k) - alphaa1*lz_V(i,j,k,1)
     &                            - betaa1 *lz_V(i,j,k,2)
               res(i,j,k) = temp
               beta1 = beta1 + temp**2
               end if
            end do; end do; end do
         end if

!$acc wait
         if (nproc.gt.1) then
            beta_v(1)=beta; beta_v(2)=beta1
            call MPI_ALLREDUCE(MPI_IN_PLACE,beta_v,2,MPI_RPREC
     &          ,MPI_SUM,COMM_TINKER,ierr)
            beta=beta_v(1); beta1=beta_v(2)
         end if
         beta  = sqrt(beta); beta1=sqrt(beta1)
         nrmbi = real(real(1,r_p)/beta,t_p)
         nrmbi1= real(real(1,r_p)/beta1,t_p)

         if (OrthBase.ne.0.and.EigenVec_l.and.kk.gt.lossOrthIt) then  ! FULL ORTHOGONALISATION
!$acc parallel loop gang vector_length(256) default(present) async
            do it = 1,kk-02
               alpha1 = 0
!$acc loop vector collapse(2)
               do k= 1,npoleloc; do i = 1,3
                  alpha1 = alpha1+ res(i,1,k)*lz_Q(i,k,it)
               end do; end do
!$acc loop vector
               do i =1, 256; if(i.eq.1) then
                  scalV(it) = real(alpha1,t_p)
               end if; end do
            end do
            if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(scalV)
               call MPI_ALLREDUCE(MPI_IN_PLACE,scalV,kk-02,MPI_TPREC
     &                           ,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
            end if
!$acc parallel loop gang vector_length(512) default(present) async
            do it = 1,kk-02
               alphaa = scalV(it)
               if (abs(alphaa).le.eps) cycle
!$acc loop vector collapse(2)
               do k= 1,npoleloc; do i = 1,3
                  temp = alphaa*lz_Q(i,k,it)
!$acc atomic
                  res(i,1,k) = res(i,1,k) - temp 
               end do; end do
            end do
            nrm_b = 0
!$acc parallel loop async collapse(2) default(present)
            do k = 1,npoleloc; do i =1,3
               nrm_b = nrm_b + res(i,1,k)**2
            end do; end do
!$acc wait
            if (nproc.gt.1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE,nrm_b,1,MPI_RPREC
     &                        ,MPI_SUM,COMM_TINKER,ierr)
            end if
            nrm_b = sqrt(nrm_b)
            nrmbi = nrm_b
            beta  = nrm_b

            if (tinkerdebug.gt.0) then
               notho = 0
!$acc wait
!$acc parallel loop async default(present)
               do i =1,kk-02
                  if (abs(scalV(i)).gt.eps) then
                     notho=notho+1
                     !if(kk.gt.285.and.kk.lt.290) print*, i,scalV(i)
                  end if
               end do
!$acc wait
c!$acc update host(scalV)
c               call printVec("scalV",scalV,kk-40,precd=12,prec=15)
            end if
         end if   ! (ORTHOGONALISATION PROCESS)

         if (tinkerdebug.gt.0.and.rank.eq.0) then
 13      format(A,I5,4F12.7,2D16.7,2x,'northo(',I0,')')
         print 13, 'it lanzcos',kk,alpha,alpha1,betaa,betaa1
     &           ,db_scal,db_scal1,notho
         end if

      end do

c     if (krylovDim.gt.120.and.OrthBase.and.EigenVec_l) then
c        call OrthogonalizeBase1(lz_Q,3*npoleloc,3*npoleloc
c    &       ,krylovDim,120)
c     end if
      if (tinkerdebug.gt.0) then
 15   format(A,I4,A,I4,A)
      if (kk.ne.krylovDim+1) write(*,15) 'early termination of',
     &   'lanzcos ',kk,'k0(',krylovDim,')'
         if (EigenVec_l) then
            if( MdiagITmat.ne.0) then
            call chkMOrthogonality_(lz_Q,3*npoleloc,krylovDim
     &          ,3*npoleloc,'verb    ')
            call chkMOrthogonality_(lz_Q1,3*npoleloc,krylovDim
     &          ,3*npoleloc,'verb    ')
            else
            call chkOrthogonality(lz_Q(1:,1:,1),3*npoleloc,krylovDim
     &          ,3*npoleloc,'lz_Q    ')
            call chkOrthogonality(lz_Q1(1:,1:,1),3*npoleloc,krylovDim
     &          ,3*npoleloc,'lz_Q1   ')
            end if
         end if
      end if

      end subroutine

      subroutine chkOrthogonality_r4(Mat,nl,nc,ldM,name)
      use domdec ,only: nproc
      implicit none
      integer  ,intent(in):: nc,nl,ldM
      real(4)  ,intent(in):: Mat(ldM,*)
      character(8),intent(in) :: name
      integer i,j,k,ii,vl,num
      real(r_p) scal,scal_max,scal_min,scal_zer,eps,sol
      logical disp
      parameter(vl=512,eps=1d-5)
      !TODO Add parallel version
      if (nproc.gt.1) return

      num=0
      disp=merge(.true.,.false.,name.eq.'verb    ')
      scal_max=-huge(scal_max)
      scal_min= huge(scal_min)
      scal_zer= huge(scal_zer)
!$acc data copyin(num,scal_max,scal_min,scal_zer) async
!$acc parallel loop gang vector_length(vl) collapse(2)
!$acc&         default(present) async
      do i = 1,nc; do j = 1,nc
         if (i.gt.j) cycle
         scal = 0
!$acc loop vector
         do k = 1,nl
            scal = scal + Mat(k,i)*Mat(k,j)
         end do
!$acc loop vector
         do k = 1,vl; if (k.eq.1) then
            sol = merge(1.0,0.0,i.eq.j)
            if (abs(scal-sol).gt.eps) then
               num= num + 1
               scal_min= min(scal_min,scal)
               scal_max= max(scal_max,scal)
               scal_zer= min(scal_zer,abs(scal))
               if (disp) print*,'orth',i,j,scal
            end if
         end if; end do
      end do; end do;
!$acc serial async copyin(name)
      if (num.gt.0) then
         print*,"IssueOrth ",name,num,real(scal_min,4),real(scal_max,4)
     &         ,real(scal_zer,4)
      end if
!$acc end serial
!$acc end data
      end subroutine

      subroutine chkMOrthogonality_(Mat,nl,nc,ldM,name)
      use domdec ,only: nproc
      use polar_temp,only: diag
      implicit none
      integer  ,intent(in):: nc,nl,ldM
      real(t_p),intent(in):: Mat(ldM,*)
      character(8),intent(in) :: name
      integer i,j,k,ii,vl,num
      real(r_p) scal,scal_max,scal_min,scal_zer,eps,sol
      logical disp
#if TINKER_DOUBLE_PREC
      parameter(vl=512,eps=1d-8)
#else
      parameter(vl=512,eps=1d-7)
#endif
      !TODO Add parallel version
      if (nproc.gt.1) return

      num =0
      disp=merge(.true.,.false.,name.eq.'verb    ')
      scal_max=-huge(scal_max)
      scal_min= huge(scal_min)
      scal_zer= huge(scal_zer)
!$acc data copyin(num,scal_max,scal_min,scal_zer) async
!$acc parallel loop gang vector_length(vl) collapse(2)
!$acc&         default(present) async
      do i = 1,nc; do j = 1,nc
         if (i.gt.j) cycle
         scal = 0
!$acc loop vector
         do k = 1,nl
            scal = scal + Mat(k,i)*Mat(k,j)/diag((k-1)/3+1)
         end do
!$acc loop vector
         do k = 1,vl; if (k.eq.1) then
            sol = merge(1.0,0.0,i.eq.j)
            if (abs(scal-sol).gt.eps) then
               num= num + 1
               scal_min= min(scal_min,scal)
               scal_max= max(scal_max,scal)
               scal_zer= min(scal_zer,abs(scal))
               if (disp) print*,'orthM',i,j,scal
            end if
         end if; end do
      end do; end do;
!$acc serial async copyin(name)
      if (num.gt.0) then
         print*,"IssueOrthm ",name,num,real(scal_min,4),real(scal_max,4)
     &         ,real(scal_zer,4)
      end if
!$acc end serial
!$acc end data
      end subroutine

      subroutine chkOrthogonality_r8(Mat,nl,nc,ldM,name)
      use domdec ,only: nproc
      implicit none
      integer,intent(in):: nc,nl,ldM
      real(8),intent(in):: Mat(ldM,*)
      character(8),intent(in) :: name
      integer i,j,k,ii,vl,num
      real(r_p) scal,scal_max,scal_zer,scal_min,eps,sol
      logical disp
      parameter(vl=512,eps=1d-8)

      !TODO Add parallel version
      if (nproc.gt.1) return

      num =0
      disp=merge(.true.,.false.,name.eq.'verb    ')
      scal_max=-huge(scal_max)
      scal_min= huge(scal_min)
      scal_zer= huge(scal_zer)
!$acc data copyin(num,scal_max,scal_min,scal_zer) async
!$acc parallel loop gang vector_length(vl) collapse(2)
!$acc&         default(present) async
      do i = 1,nc; do j = 1,nc
         if (i.gt.j) cycle
         scal = 0
!$acc loop vector
         do k = 1,nl
            scal = scal + Mat(k,i)*Mat(k,j)
         end do
!$acc loop vector
         do k = 1,vl; if (k.eq.1) then
            sol = merge(1.0,0.0,i.eq.j)
            if (abs(scal-sol).gt.eps) then
               num = num + 1
               scal_min= min(scal_min,scal)
               scal_max= max(scal_max,scal)
               scal_zer= min(scal_zer,abs(scal))
               if (disp) print*,'orth ',i,j,scal
            end if
         end if; end do
      end do; end do;
!$acc serial async copyin(name)
      if (num.gt.0) then
         print*,"IssueOrth  ",name,num,real(scal_min,4),real(scal_max,4)
     &         ,real(scal_zer,4)
      end if
!$acc end serial
!$acc end data

      end subroutine

      subroutine OrthogonalizeBase1(Base,ldb,nl,nv,start0)
      use domdec
      use orthogonalize
      use mpi
      use sizes,only: tinkerdebug
      implicit none
      integer  ,intent(in)   :: nl,nv,ldb,start0
      real(t_p),intent(inout):: Base(ldb,*)
      integer i,j,k,i1,j1,vl,ierr,str
      real(r_p) a_scal, a_nrm
      real(t_p) aScal, aNrm
      real(t_p),pointer:: v_scal(:)
      real(t_p),parameter::eps=real(1d-8,t_p)
      parameter(vl=256)

      v_scal(1:nv) => lz_WorkSpace(1:nv)

!$acc host_data use_device(v_scal,Base)
      do i = 1,start0
         if (nproc.eq.1) then
!$acc parallel loop gang vector_length(vl) async
!$acc&         deviceptr(v_scal,Base)
         do j = start0+1,nv  ! Inner product <[1..start0] | i>
            a_scal=0
!$acc loop vector
            do k = 1,nl; a_scal= a_scal+ Base(k,j)*Base(k,i); end do
            aScal=a_scal
!$acc loop vector
            do k = 1,nl; Base(k,j)= Base(k,j)- aScal*Base(k,i); end do
         end do
         else
!$acc parallel loop gang vector_length(vl) async
!$acc&         deviceptr(v_scal,Base)
         do j = start0+1,nv  ! Inner product <[1..start0] | i>
            a_scal=0
!$acc loop vector
            do k = 1,nl; a_scal= a_scal+ Base(k,j)*Base(k,i); end do
!$acc loop vector
            do k = 1,vl; if (k.eq.1) then;
               v_scal(j)=real(a_scal,t_p)
            end if; end do;
         end do
         if (nproc.gt.1) then
!$acc wait
            call MPI_ALLREDUCE(MPI_IN_PLACE,v_scal(start0+1),nv-start0 ! Parallel Inner Product
     &                        ,MPI_TPREC,MPI_SUM,COMM_TINKER,ierr)
         end if
!$acc parallel loop gang vector_length(vl) async
!$acc&         deviceptr(v_scal,Base)
         do j = start0+1,nv  ! Orthogonalize latest vectors with vec i
            aScal=v_scal(j)
            if (abs(aScal).le.eps) cycle
!$acc loop vector
            do k = 1,nl; Base(k,j)= Base(k,j)- aScal*Base(k,i); end do
         end do
         end if
      end do
!$acc end host_data
      call OrthogonalizeBase(Base(1,start0),ldb,nl,nv-start0+1)

      end subroutine

      subroutine OrthogonalizeBase(Base,ldb,nl,nv)
      use domdec
      use orthogonalize
      use mpi
      use sizes,only: tinkerdebug
      implicit none
      integer  ,intent(in)   :: nl,nv,ldb
      real(t_p),intent(inout):: Base(ldb,*)
      integer i,j,k,i1,j1,vl,ierr
      real(r_p) a_scal, a_nrm
      real(t_p) aScal, aNrm
      real(t_p),parameter::eps=real(1d-8,t_p)
      real(t_p),pointer:: v_scal(:)
      parameter(vl=256)

      v_scal(1:nv) => lz_WorkSpace(1:nv)

!$acc host_data use_device(v_scal,Base)

      a_nrm=0
!$acc parallel loop async
!$acc&         deviceptr(Base)
      do j = 1,nl; a_nrm= a_nrm+ Base(j,1)**2; end do ! Compute norm
!$acc wait
      if (nproc.gt.1) call MPI_ALLREDUCE(MPI_IN_PLACE,a_nrm,1
     &                    ,MPI_RPREC,MPI_SUM,COMM_TINKER,ierr)
      a_nrm = sqrt(a_nrm) ! Comp norm 2
      aNrm  = a_nrm       ! Cast
!$acc parallel loop async
!$acc&         deviceptr(Base)
      do j = 1,nl; Base(j,1)= Base(j,1)/aNrm; end do ! Normalize

      do i = 1,nv-1
!$acc parallel loop gang vector_length(vl) async
!$acc&         deviceptr(v_scal,Base)
         do j = i+1,nv  ! Inner product <[i+1..nv] | i>
            a_scal=0
!$acc loop vector
            do k = 1,nl; a_scal= a_scal+ Base(k,j)*Base(k,i); end do
!$acc loop vector
            do k = 1,vl; if (k.eq.1) then;
               v_scal(j)=real(a_scal,t_p)
            end if; end do;
         end do
         if (nproc.gt.1) then
!$acc wait
            call MPI_ALLREDUCE(MPI_IN_PLACE,v_scal(i+1),nv-i ! Parallel Inner Product
     &                        ,MPI_TPREC,MPI_SUM,COMM_TINKER,ierr)
         end if
!$acc parallel loop gang vector_length(vl) async
!$acc&         deviceptr(v_scal,Base)
         do j = i+1,nv  ! Orthogonalize latest vectors with vec i
            aScal=v_scal(j)
            if (abs(aScal).le.eps) cycle
!$acc loop vector
            do k = 1,nl; Base(k,j)= Base(k,j)- aScal*Base(k,i); end do
         end do

         i1 = i+1
         a_nrm = 0
!$acc parallel loop async
!$acc&         deviceptr(Base)
         do j = 1,nl; a_nrm= a_nrm+ Base(j,i1)**2; end do ! Compute norm
!$acc wait
         if (nproc.gt.1) call MPI_ALLREDUCE(MPI_IN_PLACE,a_nrm,1
     &                       ,MPI_RPREC,MPI_SUM,COMM_TINKER,ierr)
         a_nrm = sqrt(a_nrm) ! Comp norm 2
         aNrm  = a_nrm       ! Cast
         if(tinkerdebug.gt.0) print*, 'orthogonalize',i,a_nrm
!$acc parallel loop async
!$acc&         deviceptr(Base)
         do j = 1,nl; Base(j,i1)= Base(j,i1)/aNrm; end do ! Normalize
      end do

!$acc end host_data

      end subroutine

      subroutine OrthogonalizeLastVec1(Base,ldb,nl,nv)
      use domdec
      use orthogonalize
      use mpi
      use sizes     ,only: tinkerdebug
      use tinheader ,only: re_p
      implicit none
      integer  ,intent(in)   :: nl,nv,ldb
      real(t_p),intent(inout):: Base(ldb,*)
      integer i,j,k,i1,j1,vl,ierr,str
      real(r_p) a_scal, a_nrm
      real(t_p) aScal, aNrm
      real(t_p),pointer:: v_scal(:)
      real(t_p),parameter::eps=real(1d-8,t_p)
      parameter(vl=256)

      v_scal(1:nv) => lz_WorkSpace(1:nv)

!$acc host_data use_device(v_scal,Base)
!$acc parallel loop gang vector_length(vl) async
!$acc&         deviceptr(v_scal,Base)
      do j = 1,nv  ! Inner product <[1..start0] | i>
         a_scal=0
!$acc loop vector
         do k = 1,nl; a_scal= a_scal+ Base(k,j)*Base(k,nv); end do
!$acc loop vector
         do k = 1,vl; if (k.eq.1) then;
            v_scal(j)=real(a_scal,t_p)
         end if; end do;
      end do
      if (nproc.gt.1) then
!$acc wait
         call MPI_ALLREDUCE(MPI_IN_PLACE,v_scal,nv ! Parallel Inner Product
     &                     ,MPI_TPREC,MPI_SUM,COMM_TINKER,ierr)
      end if
!$acc parallel loop gang vector_length(vl) async
!$acc&         deviceptr(v_scal,Base)
      do j = 1,nv-1  ! Inner product <[1..start0] | i>
         aScal=v_scal(j)/v_scal(nv)
!$acc loop vector
         do k = 1,nl;
!$acc atomic
            Base(k,nv)= Base(k,nv)- aScal*Base(k,j);
         end do
      end do
!$acc end host_data

      call InnerProd(Base(1,nv),Base(1,nv),nl,a_nrm)
!$acc wait
      aNrm = real(1.0_re_p/sqrt(a_nrm),t_p)
      call ascale(nl,aNrm,Base(1,nv),Base(1,nv))

      end subroutine

      subroutine lejaOrdering(array,n)
      use orthogonalize
      implicit none
      integer  ,intent(in)::n
      real(r_p),intent(inout)::array(n)
      integer i,j,k,ind
      real(r_p) elt, maxv
      real(r_p) work(n)

      ! Get First value
      call getMaxr(array,n,maxv,ind)
      if (ind.ne.1) then
         array(ind)= array(1)
         array(1)  = maxv
      end if
      work(:) = 1

      do i = 2,n
         elt = array(i-1)
         do j = i,n ! Compute sequence
            work(j) = work(j)* abs(array(j)-elt)
         end do
         call getMaxr(work(i),n-i+1,maxv,ind)
         ind = ind + i-1
         if (ind.ne.i) then ! Switch array and workSpace
            elt        = array(ind)
            array(ind) = array(i)
            array(i)   = elt
            elt       = work(ind)
            work(ind) = work(i)
            work(i)   = elt
         end if
      end do

      end subroutine

#ifdef _OPENACC
      subroutine init_cublas_handle
      use orthogonalize
      use utilgpu
      integer ierr

      if (init_cublas_handle_l) then
         return
      else
         init_cublas_handle_l=.true.
      end if

      ierr = cublasCreate(cBhandle)
      ierr = ierr + cublasSetStream(cBhandle, rec_stream)
      if (ierr.ne.0) then
         print*, 'cuBlas init fails ',ierr
      end if
      end subroutine
#endif

      subroutine getEigenVec
      use domdec ,only: rank,nproc,COMM_TINKER
      use mpole
      use mpi
      use orthogonalize
      use inform ,only: deb_path,minmaxone
      use tinmemory, only: prmem_request
      use random_mod
      use sizes  ,only: tinkerdebug
      implicit none
      integer kd,ldAB,i,info,j,lda,i1,i2
      real(t_p) alpha,beta
      parameter( kd=1, ldAB=kd+1, alpha=1, beta=0 )
      real(r_p)  work(3*krylovDim)
      real(r_p) AB(ldAB,krylovDim),Z(krylovDim,krylovDim)
      real(r_p) nrm,nrm1,scal,scal1,scal2

      if (deb_path) print*, ' getEigenVec'

#if 0
      call prmem_requestm (lz_WS_,2*maxKD**2)
      i = maxKD**2
       lz_MT(1:maxKD,1:maxKD) =>lz_WS_(1:i)
      lz_MT1(1:maxKD,1:maxKD) =>lz_WS_(i+1:2*i)
#endif

      do i = 1,krylovDim
         AB(1,i) = lz_T(i)
         AB(2,i) = lz_T(maxKD+i)
      end do

#if !TINKERHP_REL_BUILD
      call M_sbev('V','L',krylovDim,kd,AB,ldAB
     &           ,EigVal,Z,maxKD,work,info)
      lz_Z(1:krylovDim,1:krylovDim) = Z(:,:)
#endif
      if (info.ne.0) then
         print*, ' M_sbev fails with error',info 
      end if


      do i = 1,krylovDim
         AB(1,i) = lz_T1(i)
         AB(2,i) = lz_T1(maxKD+i)
      end do

#if !TINKERHP_REL_BUILD
      call M_sbev('V','L',krylovDim,kd,AB,ldAB
     &           ,EigVal1,Z,maxKD,work,info)
#endif
      lz_Z1(1:krylovDim,1:krylovDim) = Z(:,:)

      if (.not.EigenVec_l) return

 12   format('eig val',50F10.6,/)
 13   format('EigV',I4,' nrm2(',2F9.5,') scal(',2F11.7,')',2I3,F11.7)

!$acc update device(EigVal,lz_Z,lz_Z1) async

#ifdef _OPENACC
      lda = 3*npoleloc
!$acc host_data use_device(lz_Q,lz_Q1,lz_Z,lz_Z1,Vp,Vp1)
      info= M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                  ,lda,nEigenV,krylovDim
     &                  ,alpha,lz_Q,lda,lz_Z,maxKD
     &                  ,beta ,  Vp,lda)
      info= info+ M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                  ,lda,nEigenV,krylovDim
     &                  ,alpha,lz_Q1,lda,lz_Z1,maxKD
     &                  ,beta ,  Vp1,lda)
!$acc end host_data
      if (info.ne.CUBLAS_STATUS_SUCCESS) then
         print '(A,A,I4,A,I0)'
     &        ,' M_cublasGemm ',__FILE__,__LINE__
     &        ,' fails with error ',info
      end if

      if (tinkerdebug.gt.0) then
      call chkOrthogonality(lz_Z,krylovDim,krylovDim
     &     ,maxKD,"lz_Z    ")
      call chkOrthogonality(lz_Z1,krylovDim,krylovDim
     &     ,maxKD,"lz_Z1   ")
         if (MdiagITmat) then
      call chkMOrthogonality_(Vp(1,1,1),3*npoleloc,nEigenV
     &     ,3*npoleloc,"Vp_prev ")
      call chkMOrthogonality_(Vp1(1,1,1),3*npoleloc,nEigenV
     &     ,3*npoleloc,"Vp1_prev ")
         else
      call chkOrthogonality(Vp(1,1,1),3*npoleloc,nEigenV
     &     ,3*npoleloc,"Vp_prev ")
      call chkOrthogonality(Vp1(1,1,1),3*npoleloc,nEigenV
     &     ,3*npoleloc,"Vp1_prev ")
         end if
      end if

      if (krylovDim.gt.100.or.MdiagITmat.and.OrthBase) then
         call OrthogonalizeBase(Vp(1,1,1),3*npoleloc,3*npoleloc,nEigenV)
      end if

      call calcEigenVectorImage
      call BuildEMat

      if (tinkerdebug.gt.0) then
         if (nproc.eq.1) then
         if (MdiagITmat) then
         call chkMOrthogonality_(Vp(1,1,1),3*npoleloc,nEigenV
     &        ,3*npoleloc,"Vp  orth")
         call chkMOrthogonality_(Vp1(1,1,1),3*npoleloc,nEigenV
     &        ,3*npoleloc,"Vp1 orth")
         else
         call chkOrthogonality(Vp(1,1,1),3*npoleloc,nEigenV
     &        ,3*npoleloc,"Vp  orth ")
         call chkOrthogonality(Vp1(1,1,1),3*npoleloc,nEigenV
     &        ,3*npoleloc,"Vp1 orth ")
         end if
!$acc wait
         else
!!$acc update host(Vp,lz_Q) async
!$acc wait
         do i = 1,nEigenV
            nrm= 0; scal= 0; nrm1= 0; scal1= 2; scal2=0;
            i1   = int(nEigenV*random()) +1
            i2   = int(nEigenV*random()) +1
            if (i1.eq.i) i1 = mod(i1+1,nEigenV) +1
            if (i2.eq.i) i2 = mod(i2+2,nEigenV) +1
!$acc parallel loop default(present)
            do j = 1,lda
               nrm  = nrm   + Vp  (j,1,i)**2
               nrm1 = nrm1  + lz_Q(j,1,i)**2
               scal = scal  + Vp  (j,1,i)*Vp(j,1,i1)
              scal2 = scal2 + Vp  (j,1,i)*Vp(j,1,i2)
              scal1 = min(scal1 , real(abs(lz_Q(j,1,i)),r_p))
            end do
            !call minmaxone(Vp(1,1,i),3*npoleloc,'Vp')
            call MPI_ALLREDUCE(MPI_IN_PLACE,nrm,1,MPI_RPREC,MPI_SUM,
     &                         COMM_TINKER,info)
            call MPI_ALLREDUCE(MPI_IN_PLACE,scal,1,MPI_RPREC,MPI_SUM,
     &                         COMM_TINKER,info)
            if (rank.eq.0) print 13, i,nrm,nrm1,scal,scal2,i1,i2,scal1
         end do
         end if
      end if
#endif

      end subroutine

      subroutine calcEigenVectorImage
      use domdec
      use inform     ,only: deb_path
#if TINKER_MIXED_PREC && defined(_OPENACC)
      use interfaces ,only: cuPOTRI,cuPOTRF,cuGESV
     &                ,cuPOTRIm
#endif
      use mpi
      use mpole      ,only: npoleloc,npolebloc
      use orthogonalize
#ifdef _OPENACC
      use utilgpu    ,only: rec_stream
#endif
      use sizes      ,only: tinkerdebug
      implicit none
      integer i,j,k,it,it1
      integer lda
      real(t_p),pointer:: work(:,:,:,:)
      !real(t_p) eps
      integer MdiagITmat_s
      interface
      subroutine MatxVec(Rin,Rout)
      real(t_p),contiguous:: Rin(:,:,:),Rout(:,:,:)
      end subroutine
      end interface

      if (deb_path) print*,'  calcEigenVectorImage'
      lz_V(1:3,1:2,1:npolebloc,1:2) => lz_WorkSpace(1:12*npolebloc)

      MdiagITmat_s  = MdiagITmat
      MdiagITmat    = 0
      saveResAmat_l = .false.

c      do it = 1,nEigenV
c         !eps = -huge(eps)
c!$acc parallel loop async collapse(3) default(present)
c         do k = 1,npoleloc; do j = 1,2; do i = 1,3
c            if (btest(j,0)) then; lz_V(i,j,k,1) = Vp (i,k,it)
c            else;                 lz_V(i,j,k,1) = Vp1(i,k,it); end if
c         end do; end do; end do
c         call Matxvec(lz_V(:,:,:,1),lz_V(:,:,:,2))
c!$acc parallel loop async collapse(3) default(present)
c         do k = 1,npoleloc; do j=1,2; do i = 1,3
c            if (btest(j,0)) then; TVp (i,k,it) = lz_V(i,j,k,2)
c            else;                 TVp1(i,k,it) = lz_V(i,j,k,2); end if
c            !eps  = max(eps,abs(lz_V(i,1,k,1)-lz_V(i,2,k,2)))
c         end do; end do; end do
cc!$acc wait
cc         print*, 'diff Tmat vec',it,eps
c      end do

      MdiagITmat=MdiagITmat_s

#if _OPENACC
!$acc host_data use_device(lz_AQ,lz_AQ1,lz_Z,lz_Z1,TVp,TVp1)
      lda = 3*npoleloc
      it = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                 ,lda,nEigenV,krylovDim
     &                 ,1.0_ti_p,lz_AQ,lda,lz_Z,maxKD
     &                 ,0.0_ti_p,TVp,lda)
      it = it+ M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                 ,lda,nEigenV,krylovDim
     &                 ,1.0_ti_p,lz_AQ1,lda,lz_Z1,maxKD
     &                 ,0.0_ti_p,  TVp1,lda)
!$acc end host_data
      if (it.ne.CUBLAS_STATUS_SUCCESS) then
         print '(A,A,I4,A,I0)'
     &         ,' M_cublasGemm ',__FILE__,__LINE__
     &         ,' fails with error ',it
      end if
#endif

      end subroutine

      subroutine BuildEMat
      use domdec
      use inform     ,only: deb_path
#if defined(_OPENACC)
      use interfaces ,only: cuPOTRI,cuPOTRF,cuGESV
#if TINKER_MIXED_PREC
     &               ,cuPOTRIm
#endif
#endif
      use mpi
      use mpole      ,only: npoleloc,npolebloc
      use orthogonalize
#ifdef _OPENACC
      use utilgpu    ,only: rec_stream
#endif
      use sizes      ,only: tinkerdebug
      implicit none
      integer i,j,k,it,it1
      real(t_p),pointer:: work(:,:,:,:)
      real(t_p) eps
      integer sd,sdec,ndec,ndec1,stdec,ledec,dim1
      logical MdiagITmat_s
      real(r_p) tot1,tot2,tmp
      real(t_p) tot
      parameter(sdec=256)

      call setMIdk

      dim1  = 3*npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)
c
c     MatE =  Vp^{T} * Tmat * Vp
c
!$acc parallel loop gang vector_length(sdec) collapse(2)
!$acc&         default(present) private(tot1) async
      do it = 1,nEigenV; do it1 = 1,nEigenV;

         if (it1.lt.it) cycle
         tot1  = 0
         tot2  = 0
         !stdec = (sd-1)*sdec+1
         !ledec = min(sd*sdec,dim1)
!$acc loop vector
         do k = 1,dim1
            tot1 = tot1 + Vp (k,1,it1)*TVp (k,1,it)
            tot2 = tot2 + Vp1(k,1,it1)*TVp1(k,1,it)
         end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!!$acc atomic
            MatE(it1,it) = MatE(it1,it) + tot1
            MatE1(it1,it) = MatE1(it1,it) + tot2
            if(it.ne.it1) MatE (it,it1) = MatE (it,it1) + tot1
            if(it.ne.it1) MatE1(it,it1) = MatE1(it,it1) + tot2
         end if; end do;

      end do; end do

      if (nproc.gt.1) then
!$acc wait
      call MPI_ALLREDUCE(MPI_IN_PLACE,MatE,2*nEigenV**2,MPI_RPREC
     &                  ,MPI_SUM,COMM_TINKER,it)
      end if

      ! Casting operation in both MatE and MatE1
!$acc parallel loop async default(present)
      do i=1,2*nEigenV**2
         MatE_(i,1) = MatE(i,1)
         !MatE1_(i,1) = MatE(i,1)
      end do

      ! MatE Cholesky Inversion
#ifdef _OPENACC
      ! Cholesky Factorisation & Inversion
!$acc host_data use_device(MatE_,MatE1_)
      call cuPOTRF(nEigenV,MatE_,nEigenV,rec_stream)
      call cuPOTRI(nEigenV,MatE_,nEigenV,rec_stream)
      call cuPOTRF(nEigenV,MatE1_,nEigenV,rec_stream)
      call cuPOTRI(nEigenV,MatE1_,nEigenV,rec_stream)
!$acc end host_data

       ! LU Decomposition & Inversion
c!$acc host_data use_device(MatE_,MatE1_,Ipivot,MIdk)
c      call cuGESV(nEigenV,nEigenV,MatE_,nEigenV,Ipivot,MIdk,nEigenV
c     &           ,rec_stream)
c!$acc parallel loop async collapse(2) deviceptr(MatE_,MIdk)
c      do i = 1,nEigenV**2
c         do j= 1,nEigenV
c            MatE_(j,i) = MIdk(j,i)
c            Midk (j,i) = merge(1,0,(i.eq.j))
c         end do
c      end do
c      call cuGESV(nEigenV,nEigenV,MatE1_,nEigenV,Ipivot,MIdk,nEigenV
c     &           ,rec_stream)
c!$acc parallel loop async collapse(2) deviceptr(MatE_,MIdk)
c      do i = 1,nEigenV**2
c         do j= 1,nEigenV
c            MatE1_(j,i) = MIdk(j,i)
c            Midk (j,i) = merge(1,0,(i.eq.j))
c         end do
c      end do
c!$acc end host_data

#if TINKER_DOUBLE_PREC
      if (tinkerdebug.gt.0) then
      tot1=1.0; tot2=0.0
!$acc host_data use_device(MatE,MatE_,MatE1,MatE1_,MIdk)
      it = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                 ,nEigenV,nEigenV,nEigenV
     &                 ,tot1,MatE,nEigenV,MatE_,nEigenV
     &                 ,tot2,MIdk,nEigenV)
!$acc end host_data
      if (it.ne.CUBLAS_STATUS_SUCCESS) then
         print '(A,A,I4,A,I0)'
     &         ,' M_cublasGemm ',__FILE__,__LINE__
     &         ,' fails with error ',it
      end if
      tot1 = 1000
!$acc wait
!$acc update host(MatE_,MatE,MatE1_,MatE1,MIdk)
      do i = 1,nEigenV; tot1=min(tot1,MatE(i,i)); end do;
      if (.true.) then
         call printVec(" MatE   ", MatE, nEigenV**2, prec=12, precd=8
     &        ,period=min(nEigenV,20))
         call printVec(" MatE1  ", MatE1, nEigenV**2, prec=12, precd=8
     &        ,period=min(nEigenV,20))
      end if
c     call printVec(" MatE_ ",MatE_,nEigenV**2, prec=12, precd=8,
c    &      period=min(nEigenV,20))
      call chkOrthogonality(MIdk,nEigenV,nEigenV,nEigenV,'E inv   ')
      print*, 'min E diag', tot1
      end if
#endif
#endif

      end subroutine
c
c     Compute projection of Vin in the deflated subspace
c
      subroutine projectorOpe(Vin,Vout,invert)
      use domdec
      use orthogonalize
      use mpole
      use mpi
      use inform
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)
      integer  ,intent(in)   :: invert

      integer ii,i,j,k,dim1
      integer sd,sdec,ndec,ndec1,stdec,ledec
      real(r_p) tot1,tot2,tmp
      real(t_p) tot,tot_
      parameter(sdec=64)

      dim1  = npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)

!$acc kernels async present(Evec)
      Evec(1:2*nEigenV,1) = 0
!$acc end kernels

!$acc parallel loop gang vector_length(sdec) collapse(2)
!$acc&         default(present) private(tot1,tot2) async
      do ii = 1,nEigenV; do sd = 1,ndec

            tot1  = 0
            tot2  = 0
         stdec = (sd-1)*sdec+1
         ledec = min(sd*sdec,dim1)
!$acc loop vector collapse(2)
         do k = stdec,ledec; do i = 1,3
            tot1 = tot1 + Vp (i,k,ii)*Vin(i,1,k)
            tot2 = tot2 + Vp1(i,k,ii)*Vin(i,2,k)
         end do; end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!$acc atomic
            Evec(ii,1) = Evec(ii,1) + tot1
!$acc atomic
            Evec(ii,2) = Evec(ii,2) + tot2
         end if; end do;

      end do; end do

      if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(Evec)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Evec,nEigenV*2,MPI_RPREC 
     &       ,MPI_SUM,COMM_TINKER,i)
!$acc end host_data
      end if

!$acc kernels async default(present)
      Evec_(1:2*nEigenV,1) = Evec(1:2*nEigenV,1)
!$acc end kernels

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0;
         if(btest(j,0)) then;
!$acc loop seq
            do ii = 1,nEigenV
              tot = tot + Vp (i,k,ii)*Evec_(ii,j)
            end do
         else
!$acc loop seq
            do ii = 1,nEigenV
              tot = tot + Vp1(i,k,ii)*Evec_(ii,j)
            end do
         end if
         Vout(i,j,k) = Vin(i,j,k) - tot
      end do; end do; end do;

      call minmaxone(Vout(1:,1,1),6*npoleloc,'Vout')

      end subroutine

      subroutine ProjectorPGeneric(Vin,Vout)
      use domdec  ,only: COMM_TINKER,rank,nproc
      use orthogonalize
      use mpi
      use mpole
      use inform
      use tinheader
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)

      integer ii,i,j,k,ierr
      integer sd,sdec,ndec,ndec1,stdec,ledec,dim1
      real(r_p) tot1,tot2,tmp
      real(t_p) tot,tot_
      parameter(sdec=64)

      dim1  = npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)

      print*, 'ProjectorPGeneric not updated yet !!!'

!$acc kernels async present(Evec)
      Evec(1:2*nEigenV,1) = 0
!$acc end kernels

!$acc parallel loop gang vector_length(sdec) collapse(2)
!$acc&         default(present) private(tot1,tot2) async
      do ii = 1,nEigenV; do sd = 1,ndec

            tot1  = 0
            tot2  = 0
         stdec = (sd-1)*sdec+1
         ledec = min(sd*sdec,dim1)
!$acc loop vector collapse(2)
         do k = stdec,ledec; do i = 1,3
            tot1 = tot1 + Vp(i,k,ii)*Vin(i,1,k)
            tot2 = tot2 + Vp(i,k,ii)*Vin(i,2,k)
         end do; end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!$acc atomic
            Evec(ii,1) = Evec(ii,1) + tot1
!$acc atomic
            Evec(ii,2) = Evec(ii,2) + tot2
         end if; end do;

      end do; end do

      if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(Evec)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Evec,nEigenV*2,MPI_RPREC 
     &       ,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
      end if

!$acc kernels async default(present)
      Evec_(1:2*nEigenV,1) = Evec(1:2*nEigenV,1)
!$acc end kernels

#ifdef _OPENACC
!$acc host_data use_device(MatE_,Evec_,MIdk)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                   ,nEigenV,2,nEigenV
     &                   ,1.0_ti_p,MatE_,nEigenV,Evec_,nEigenV
     &                   ,0.0_ti_p,MIdk,nEigenV)
!$acc end host_data
#endif

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0
!$acc loop seq
         do ii = 1,nEigenV
            tot = tot + TVp(i,k,ii)*MIdk(ii,j)
         end do
         Vout(i,j,k) = Vin(i,j,k) - tot
      end do; end do; end do;

      end subroutine

      subroutine ProjectorPTransGeneric(Vin,Vout)
      use domdec  ,only: COMM_TINKER,rank,nproc
      use orthogonalize
      use mpi
      use mpole
      use inform
      use tinheader
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)

      integer ii,i,j,k,ierr
      integer sd,sdec,ndec,ndec1,stdec,ledec,dim1
      real(r_p) tot1,tot2,tmp
      real(t_p) tot,tot_
      parameter(sdec=64)

      dim1  = npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)

!$acc kernels async present(Evec)
      Evec(1:2*nEigenV,1) = 0
!$acc end kernels

!$acc parallel loop gang vector_length(sdec) collapse(2)
!$acc&         default(present) private(tot1,tot2) async
      do ii = 1,nEigenV; do sd = 1,ndec

            tot1 = 0
            tot2 = 0
         stdec = (sd-1)*sdec+1
         ledec = min(sd*sdec,dim1)
!$acc loop vector collapse(2)
         do k = stdec,ledec; do i = 1,3
            tot1 = tot1 + TVp (i,k,ii)*Vin(i,1,k)
            tot2 = tot2 + TVp1(i,k,ii)*Vin(i,2,k)
         end do; end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!$acc atomic
            Evec(ii,1) = Evec(ii,1) + tot1
!$acc atomic
            Evec(ii,2) = Evec(ii,2) + tot2
         end if; end do;

      end do; end do

      if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(Evec)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Evec,nEigenV*2,MPI_RPREC 
     &       ,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
      end if

      if (.false.) then
!$acc parallel loop collapse(2) present(Evec_,Evec,EigVal) async
      do j = 1,2; do i = 1,nEigenV
         if (btest(j,0)) then; MIdk(i,j) = Evec(i,j)/EigVal (i)
         else;                 MIdk(i,j) = Evec(i,j)/EigVal1(i); end if
      end do; end do
      else
!$acc kernels async default(present)
      Evec_(1:2*nEigenV,1) = Evec(1:2*nEigenV,1)
!$acc end kernels

#ifdef _OPENACC
!$acc host_data use_device(MatE_,MatE1_,Evec_,MIdk)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_T,CUBLAS_OP_N
     &                   ,nEigenV,2,nEigenV
     &                   ,1.0_ti_p,MatE_,nEigenV,Evec_,nEigenV
     &                   ,0.0_ti_p,MIdk,nEigenV)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_T,CUBLAS_OP_N
     &                   ,nEigenV,2,nEigenV
     &                   ,1.0_ti_p,MatE1_,nEigenV,Evec_(1,2),nEigenV
     &                   ,0.0_ti_p,MIdk(1,2),nEigenV)
!$acc end host_data
#endif
      end if

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0;
         if (btest(j,0)) then
!$acc loop seq
            do ii = 1,nEigenV
               tot = tot + Vp(i,k,ii)*MIdk(ii,j)
            end do
         else
!$acc loop seq
            do ii = 1,nEigenV
               tot = tot + Vp1(i,k,ii)*MIdk(ii,j)
            end do
         end if
         Vout(i,j,k) = Vin(i,j,k) - tot
      end do; end do; end do;

      end subroutine
c
c     Compute Vout = alpha*Vout + Vp*EigVal^{-1}*Vp^{T}
c
      subroutine ApplyQxV (Vin,Vout,alpha)
      use domdec  ,only: COMM_TINKER,rank,nproc
      use inform  ,only: minmaxone
      use orthogonalize
      use mpi
      use mpole
      use tinheader
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)
      real(t_p),intent(in)   :: alpha
      integer ii,i,j,k,ierr
      integer sd,sdec,ndec,ndec1,stdec,ledec,dim1
      real(r_p) tot1,tot2,tmp
      real(t_p) tot,tot_
      parameter(sdec=64)

      dim1  = npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)

!$acc kernels present(Evec) async
      Evec(1:2*nEigenV,1) = 0
!$acc end kernels

!$acc parallel loop gang vector_length(sdec) collapse(2)
!$acc&         default(present) private(tot1,tot2) async
      do ii = 1,nEigenV; do sd = 1,ndec

         tot1  = 0
         tot2  = 0
         stdec = (sd-1)*sdec+1
         ledec = min(sd*sdec,dim1)
!$acc loop vector collapse(2)
         do k = stdec,ledec; do i = 1,3
            tot1 = tot1 + Vp (i,k,ii)*Vin(i,1,k)
            tot2 = tot2 + Vp1(i,k,ii)*Vin(i,2,k)
         end do; end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!$acc atomic
            Evec(ii,1) = Evec(ii,1) + tot1
!$acc atomic
            Evec(ii,2) = Evec(ii,2) + tot2
         end if; end do

      end do; end do

      if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(Evec)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Evec,nEigenV*2,MPI_RPREC 
     &       ,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
      end if

      if (.false.) then
!$acc parallel loop collapse(2) present(Evec_,Evec,EigVal) async
      do j = 1,2; do i = 1,nEigenV
         if (btest(j,0)) then; Evec_(i,j) = Evec(i,j)/EigVal (i)
         else;                 Evec_(i,j) = Evec(i,j)/EigVal1(i); end if
      end do; end do
      else

!$acc kernels async default(present)
      Evec_(1:2*nEigenV,1) = Evec(1:2*nEigenV,1)
!$acc end kernels
#ifdef _OPENACC
!$acc host_data use_device(MatE_,MatE1_,Evec_,MIdk)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                   ,nEigenV,1,nEigenV
     &                   ,1.0_ti_p,MatE_,nEigenV,Evec_,nEigenV
     &                   ,0.0_ti_p,MIdk,nEigenV)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                   ,nEigenV,1,nEigenV
     &                   ,1.0_ti_p,MatE1_,nEigenV,Evec_(1,2),nEigenV
     &                   ,0.0_ti_p,MIdk(1,2),nEigenV)
!$acc end host_data
#endif

!$acc kernels async default(present)
      Evec_(1:2*nEigenV,1) = MIdk(1:2*nEigenV,1)
!$acc end kernels

      end if

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0;
         if(btest(j,0)) then
!$acc loop seq
            do ii = 1,nEigenV
               tot =tot  +Vp (i,k,ii)*Evec_(ii,j)
            end do
         else
!$acc loop seq
            do ii = 1,nEigenV
               tot =tot  +Vp1(i,k,ii)*Evec_(ii,j)
            end do
         end if
         Vout(i,j,k) = alpha*Vout(i,j,k) + tot
      end do; end do; end do;

      end subroutine

      subroutine diagPrecnd(Vin,Vout)
      use mpole
      use polar_temp ,only: diag
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,npoleloc)
      real(t_p),intent(inout):: Vout(3,2,npoleloc)
      integer i,j,k

!$acc parallel loop collapse(3) async default(present)
      do k=1,npoleloc; do j=1,2; do i=1,3
         Vout(i,j,k) = diag(k)*Vin(i,j,k)
      end do; end do; end do

      end subroutine

      subroutine adaptedDeflat2 (Vin,Vout,Minv)
      use inform  ,only: minmaxone
      use orthogonalize
      use mpole
      use tinheader
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,npoleloc)
      real(t_p),intent(inout):: Vout(3,2,npoleloc)
      real(t_p),intent(in)   :: Minv(npoleloc)
      integer i
      real(t_p),pointer:: work(:,:,:)
      work(1:3,1:2,1:npoleloc)=> lz_WorkSpace(6*npoleloc+1:12*npoleloc)

      if(debg) print*, 'adaptedDeflat2'
      call diagPrecnd(Vin,work)
      if(debg) call minmaxone(work(1:,1,1),6*npoleloc,'work')
      if (.false.) then
         call projectorOpe(work,work,0)
      else
         call ProjectorPTransGeneric(work,work)
      end if
      if(debg) call minmaxone(work(1:,1,1),6*npoleloc,'work')
      call ApplyQxV(Vin,work,1.0_ti_p)
      if(debg) call minmaxone(work(1:,1,1),6*npoleloc,'work')
      if(debg) print*, '-----------------------------------------------'

!$acc parallel loop async default(present)
      do i=1,6*npoleloc
         Vout(i,1,1) = work(i,1,1)
      end do

      end subroutine

      subroutine restartDeflat
      use orthogonalize
      implicit none
      !call calcEigenVectorImage
      wKD = 1
      end subroutine

      subroutine loadKrylovBase(mu,Amu,beta_,beta1_)
      use orthogonalize
      use mpole
      use utilgpu  ,only: rec_queue
      implicit none
      real(t_p) mu(3,2,*)
      real(t_p) Amu(3,2,*)
      real(r_p) beta_,beta1_
      integer i,j,k
      real(t_p) beta,beta1

      ! Check out if maximum capacity has been reached
      if (wKD.gt.maxKD) return
      if (.true.) print*, 'loadKrylovBase',wKD

      if (wKD.eq.1) then
!$acc parallel loop collapse(3) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
            lz_Q (i,k,wKD) =  mu(i,j,k)
            lz_AQ(i,k,wKD) = Amu(i,j,k)
            else
            lz_Q1 (i,k,wKD) =  mu(i,j,k)
            lz_AQ1(i,k,wKD) = Amu(i,j,k)
            end if
         end do; end do; end do
      else
         beta = real(beta_ ,t_p)
         beta1= real(beta1_,t_p)
!$acc parallel loop collapse(3) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
            lz_Q  (i,k,wKD) =  mu(i,j,k)
            lz_AQ (i,k,wKD) = Amu(i,j,k) - beta*lz_AQ(i,k,wKD-1)
            else
            lz_Q1 (i,k,wKD) =  mu(i,j,k)
            lz_AQ1(i,k,wKD) = Amu(i,j,k) - beta1*lz_AQ1(i,k,wKD-1)
            end if
         end do; end do; end do
      end if

      wKD = wKD +1
      end subroutine

      subroutine FinalizeKrylovBase
      use orthogonalize
      use domdec
      use inform        ,only: deb_path
      use mpole
      use mpi
      use sizes         ,only: tinkerdebug
      implicit none
      integer ii,i,j,k,ierr
      integer sd,sdec,ndec,ndec1,stdec,ledec,dim1
      real(r_p) tot1,tot2,tmp
      real(t_p) tot,tot_
      parameter(sdec=64)

      if (wKD.lt.4) return

      if (deb_path) print*,'   FinalizeKrylovBase'
      krylovDim = wKD - 1

      if (tinkerdebug.gt.0.and.EigenVec_l) then
         call chkOrthogonality(lz_Q(1:,1:,1),3*npoleloc,krylovDim
     &       ,3*npoleloc,'lz_Q    ')
         call chkMOrthogonality_(lz_Q(1:,1:,1),3*npoleloc,krylovDim
     &       ,3*npoleloc,'lz_Q    ')
      end if

      dim1  = npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)

!$acc kernels present(MIdk) async
      MIdk(1:nEigenV**2,1) = 0
!$acc end kernels

!$acc parallel loop gang vector_length(sdec) collapse(3)
!$acc&         default(present) private(tot1,tot2) async
      do ii = 1,krylovDim; do j=1,krylovDim; do sd = 1,ndec

         tot1  = 0
         tot2  = 0
         stdec = (sd-1)*sdec+1
         ledec = min(sd*sdec,dim1)
!$acc loop vector collapse(2)
         do k = stdec,ledec; do i = 1,3
            tot1 = tot1 + lz_Q(i,k,ii)*lz_Q(i,k,j)
            !tot2 = tot2 + lz_Q(i,k,ii)*lz_Q(i,k,j)
         end do; end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!$acc atomic
            MIdk(ii,j) = MIdk(ii,j) + real(tot1,t_p)
c!$acc atomic
c            Evec(ii,2) = Evec(ii,2) + tot2
         end if; end do

      end do; end do; end do

      if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(MIdk)
         call MPI_ALLREDUCE(MPI_IN_PLACE,MIdk,nEigenV**2,MPI_TPREC 
     &       ,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
      end if

!$acc wait
!$acc update host(MIdk)
      call printVec('MIdk ch0',MIdk,nEigenV**2
     &             ,prec=13,precd=8,period=nEigenV)

      end subroutine

      subroutine polarEingenTest
      use atoms
      use domdec
      use mpole
      use orthogonalize
      use polar_temp ,only: mu
      use utilgpu    ,only: Tinker_shellEnv
      use tinheader  ,only: re_p
      implicit none
      integer sz_init,sz_max,i,jump,nlop,ire
      integer,parameter::ncount=20
      real(r_p) Eig(ncount)
      real(r_p) ri,si,sat
      real(r_p),parameter:: tol=1d-6
      integer Nracines

      ri = 0.0

      call Tinker_shellEnv('KRYLOV_INIT',sz_init,50)
      call Tinker_shellEnv('KRYLOV_JUMP',jump,3*n/nlop)
      call Tinker_shellEnv('EIG_MAX',si,15.0_re_p)
      call Tinker_shellEnv('EIG_NLOOP',nlop,1)
      sz_max = jump*(nlop-1)+sz_init

      call lanzcos_init(sz_max,0,mu)
      call lanzcos_iteration

      do i = 1,nlop
      call bissectionGivens1(lz_T,lz_T(sz_max+1:),sz_init
     &                      ,ncount,ri,si,tol,Eig)

 12   format(A,I7,A,20F10.6)
      write(*,12) 'EigenVal k(',sz_init,')',Eig(1:ncount)
      sz_init = min(sz_init + jump, 3*n)
      end do

      end subroutine

      subroutine polarEingenVal
      use atoms
      use domdec
      use mpole
      use orthogonalize
      use polar_temp ,only: mu
      use sizes      ,only: tinkerdebug
      use utilgpu    ,only: Tinker_shellEnv
      use random_mod
      implicit none
      integer i,sz_init
      integer ncount
      real(r_p) ri,si,sat
      real(r_p),allocatable:: Eig(:)
      real(r_p),parameter:: tol=1d-6
      integer Nracines
      character*20 ::fmt

      ri = 0.0
      si = 25

      call Tinker_shellEnv("KRYLOV_DIM",sz_init, min(40,3*n/100))
      call Tinker_shellEnv("NEIGEN_VEC",nEigenV, min(5,sz_init))

      call lanzcos_init(sz_init,nEigenV,mu)
      call lanzcos_iteration
      allocate(Eig(sz_init))

      if (EigenVec_l) then
         call getEigenVec
      else
         ncount = sz_init
         call getEigenVec
c        call bissectionGivens1(lz_T,lz_T(sz_init+1),sz_init
c    &                         ,ncount,ri,si,tol,EigVal)
      end if
      Eig = EigVal

      if (tinkerdebug.gt.0.and.rank.eq.0)
#if TINKER_DOUBLE_PREC
     &    call printVec('EigenVal__',Eigval,min(150,sz_init)
     &         ,prec=13,precd=9)
#else
     &    call printVec('EigenVal__',Eigval,min(150,sz_init)
     &         ,prec=11,precd=7)
#endif
      do i = 1,sz_init ! SomeHow stmv Eigenvalues disrupts BAOABRESPA1 long range
         if (Eig(i).gt.5) Eig(i) = 3.41 + i*0.005
      end do

      if (EigenVec_l) then; call lejaOrdering(Eig,sz_init)
      else; call lejaOrdering(Eig,ncount); end if

      call load_shift_attr(Eig,lz_T(sz_init+1),sz_init)

      end subroutine

