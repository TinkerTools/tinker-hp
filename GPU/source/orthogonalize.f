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
      use inform    ,only: deb_Path,abort
      use interfaces,only: tmatxb_pmegpu,tmatxb_p
      use math
      use mpole
      use polar
      use polar_temp,only: murec
     &              , dipfield,dipfieldbis
      use polpot
      use potent    ,only: use_polarshortreal
      use sizes
      use timestat
      use units
      use utils
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
      use utilgpu
      implicit none
      real(t_p),contiguous:: Rin(:,:,:),Rout(:,:,:)
      integer   i,j,k,nrhs
      integer,dimension(nproc):: reqsend,reqrec,req2send,req2rec
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
!!$acc parallel loop collapse(3) default(present) async(rec_queue)
!      do k=1,npoleloc; do j=1,nrhs; do i=1,3
!         Rout(i,j,k) = polarity(poleglob(k))*Rout(i,j,k)
!      end do; end do; end do
      end subroutine

      subroutine lanzcos_init(k0,compEv,mu)
      use domdec
      use orthogonalize
      use mpole
      use mpi
      use interfaces ,only: initcuSolver
      use polar_temp
      use tinMemory
      use utilgpu
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
      implicit none
      integer  ,intent(in):: k0,compEv
      real(t_p),intent(in):: mu(3,2,*)
      real(r_p) l2_norm,l2_norm1
      real(t_p) temp
      integer i,j,k,ierr
      logical EigV

      krylovdim = k0
      nEigenV   = compEv
      EigenVec_l= merge(.true.,.false.,nEigenV.gt.0)

      i = 3*(k0+betw) + nEigenV**2

      call prmem_requestm(lz_StoreSpace_rp,i)
      call prmem_request (lz_WorkSpace,12*npolebloc+2*krylovdim)
      call prmem_request (res,3,2,npoleloc)
      call prmem_request (dipfield,3,2,max(npoleloc,1))
      call prmem_request (dipfieldbis,3,2,max(npolerecloc,1))

      if (nproc.gt.1) then
         call prmem_request(buffermpi1  ,3,2,max(npoleloc,1))
         call prmem_request(buffermpimu1,3,2,max(npoleloc,1))
         call prmem_request(buffermpi2  ,3,2,max(npolerecloc,1))
         call prmem_request(buffermpimu2,3,2,max(npolerecloc,1))
      end if

      ! lz_StoreSpace     [ Vp  | lz_Q | MatE_ ]
      !                         | TVp
      ! lz_StoreSpace_rp  [ lz_T  | Eigval | MatE ]

      if (nEigenV.gt.0) then
         i = 3*npoleloc*(krylovdim+compEv)
         j = 3*npoleloc*compEv
         call prmem_request(lz_StoreSpace,i+nEigenV**2)
         call prmem_request(Evec_,nEigenV,2)
         call prmem_request(MIdk,nEigenV,nEigenV)
         call prmem_requestm(Evec,nEigenV,2)
           Vp(1:3,1:npoleloc,1:compEv)=>lz_StoreSpace(1:j)
          TVp(1:3,1:npoleloc,1:compEv)=>lz_StoreSpace(1+j:2*j)
         lz_Q(1:3,1:npoleloc,1:krylovdim)=>lz_StoreSpace(j+1:i)
         k = 3*k0+betw
         MatE (1:nEigenV,1:nEigenV)=> lz_StoreSpace_rp(k+1:k+nEigenV**2)
         MatE_(1:nEigenV,1:nEigenV)=> lz_StoreSpace   (i+1:i+nEigenV**2)
#ifdef _OPENACC
         call init_cublas_handle
         call initcuSolver(rec_stream)
#endif
      end if

      lz_V(1:3,1:2,1:npolebloc,1:2) => lz_WorkSpace(1:12*npolebloc)
      lz_T(1:2*k0) => lz_StoreSpace_rp(1:2*k0)
      EigVal(1:k0) => lz_StoreSpace_rp(2*k0+1+betw:3*k0+betw)

      do i = 1,2*k0
         lz_T(i) = 0
      end do

      l2_norm = 0; l2_norm1=0
!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k=1,npoleloc; do j=1,2; do i=1,3;
         if (btest(j,0)) then
            l2_norm = l2_norm + mu(i,j,k)**2
         else
            l2_norm1 = l2_norm1 + mu(i,j,k)**2
         end if
      end do; end do; end do;
!$acc wait
      if (nproc.gt.1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE,l2_norm,1,MPI_RPREC,MPI_SUM,
     &                   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,l2_norm1,1,MPI_RPREC,MPI_SUM,
     &                   COMM_TINKER,ierr)
      end if
      l2_norm=sqrt(l2_norm); l2_norm1=sqrt(l2_norm1);

      ! Normalisation
      if (EigenVec_l) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k = 1,npoleloc; do j=1,2; do i=1,3
         temp          = mu(i,1,k)/l2_norm 
         lz_V(i,j,k,1) = temp
         if (j.eq.1) lz_Q(i,k,1) = temp
      end do; end do; end do
      else
!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k = 1,npoleloc; do j=1,2; do i=1,3
         if (btest(j,0)) then
            lz_V(i,j,k,1) = mu(i,1,k)/l2_norm
         else
            lz_V(i,j,k,1) = mu(i,1,k)/l2_norm
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
      use polar_temp   ,only: res
      use utilgpu ,only: rec_queue,Tinker_shellEnv
      use sizes   ,only: tinkerdebug
      implicit none
      integer  i,j,k,kk,it,ierr
      integer lossOrthIt,makeOrth,notho
      real(t_p) temp,nrmb,alphaa
      real(r_p) alpha,alpha1,beta,beta1,beta_s,beta1_s,nrm_b,l2_norm
      real(r_p) alpha_v(2),beta_v(2)
      real(r_p),parameter:: tol=1d-12
#if TINKER_DOUBLE_PREC
      real(t_p),parameter:: eps=1d-1
#else
      real(t_p),parameter:: eps=1e-1
#endif
      real(t_p),pointer:: scalV(:)
      character*8 cfg
      parameter(lossOrthIt=200)

      interface
      subroutine MatxVec(Rin,Rout)
      real(t_p),contiguous:: Rin(:,:,:),Rout(:,:,:)
      end subroutine
      end interface

      if (deb_path) write(*,*) 'lanzcos_iteration'
      call Tinker_shellEnv("ORTH",makeOrth,0)
      scalV(1:krylovdim)=>
     &         lz_WorkSpace(12*npolebloc+1:12*npolebloc+krylovdim)
      notho = 0

      call MatxVec(lz_V(:,:,:,1),res)
      alpha = 0; alpha1 = 0
      beta  = 0; beta1  = 0
!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k=1,npoleloc; do j=1,2; do i=1,3
         if (btest(j,0)) then
         alpha = alpha + res(i,j,k)*lz_V(i,j,k,1)
c        else
c        alpha1 = alpha1 + res(i,j,k)*lz_V(i,j,k,1)
         end if
      end do; end do; end do
!$acc wait

      if (nproc.gt.1) then
         alpha_v(1) = alpha; alpha_v(2) = alpha1;
         call MPI_ALLREDUCE(MPI_IN_PLACE,alpha_v,2,MPI_RPREC,MPI_SUM,
     &        COMM_TINKER,ierr)
         alpha = alpha_v(1); alpha1 = alpha_v(2)
      end if

!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k=1,npoleloc; do j = 1,2; do i = 1,3
         if (btest(j,0)) then
            temp = res(i,j,k) - alpha*lz_V(i,j,k,1)
            res(i,j,k) = temp
            beta  = beta + temp**2
c        else
c           temp = res(i,j,k) - alpha1*lz_V(i,j,k,1)
c           res(i,j,k) = temp
c           beta1  = beta1 + temp**2
         end if
      end do; end do; end do
!$acc wait

      if (nproc.gt.1) then
         beta_v(1)=beta; beta_v(2)=beta1
         call MPI_ALLREDUCE(MPI_IN_PLACE,beta_v,2,MPI_RPREC,MPI_SUM,
     &        COMM_TINKER,ierr)
         beta=beta_v(1); beta1=beta_v(2)
      end if
      beta = sqrt(beta); beta1=sqrt(beta1)
      nrmb = beta
      lz_T(1) = alpha

      do kk = 2, krylovDim

         if (beta.lt.tol) exit
         lz_T(krylovDim+kk-1) = beta

         if (EigenVec_l) then
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
               lz_V(i,j,k,2) = lz_V(i,j,k,1)
               temp          = res(i,j,k)/nrmb
               lz_V(i,j,k,1) = temp
                lz_Q(i,k,kk) = temp
c           else
c              lz_V(i,j,k,2) = lz_V(i,j,k,1)
c              temp          = res(i,j,k)/nrmb
c              lz_V(i,j,k,1) = temp
            end if
         end do; end do; end do
         else
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
               lz_V(i,j,k,2) = lz_V(i,j,k,1)
               lz_V(i,j,k,1) = res(i,j,k)/nrmb
c           else
c              lz_V(i,j,k,2) = lz_V(i,j,k,1)
c              lz_V(i,j,k,1) = res(i,j,k)/nrmb
            end if
         end do; end do; end do
         end if

         if (tinkerdebug.gt.0) then  ! (DEBUG INFO)
         l2_norm=0
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,1; do i=1,3
            l2_norm = l2_norm + lz_V(i,1,k,2)*lz_V(i,1,k,1)
         end do; end do; end do
!$acc wait
         if (nproc.gt.1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE,l2_norm,1,MPI_RPREC,MPI_SUM,
     &           COMM_TINKER,ierr)
         end if

c        if (tinkerdebug.gt.0.and.EigenVec_l.and.kk.gt.100) then
c           !cfg=merge('verb    ','lz_Q    ',kk.gt.286.and.kk.lt.291)
c           cfg='lz_Q'
c           call chkOrthogonality(lz_Q(1,1,1),3*npoleloc,kk
c    &                           ,3*npoleloc,cfg)
c        end if
         end if  ! (DEBUG INFO)

         call MatxVec(lz_V(:,:,:,1),res)
         alpha = 0; alpha1 = 0
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
               alpha = alpha + res(i,j,k)*lz_V(i,j,k,1)
c           else
c              alpha1 = alpha1 + res(i,j,k)*lz_V(i,j,k,1)
            end if
         end do; end do; end do
!$acc wait
         if (nproc.gt.1) then
            alpha_v(1) = alpha; alpha_v(2) = alpha1;
            call MPI_ALLREDUCE(MPI_IN_PLACE,alpha_v,2,MPI_RPREC,MPI_SUM,
     &           COMM_TINKER,ierr)
            alpha = alpha_v(1); alpha1 = alpha_v(2)
         end if
         lz_T(kk) = alpha
         alpha1   = alpha
         beta_s=beta; beta1_s=beta1
         beta  = 0; beta1  = 0

         alphaa = alpha
!$acc parallel loop collapse(3) default(present) async(rec_queue)
         do k=1,npoleloc; do j=1,2; do i=1,3
            if (btest(j,0)) then
            temp = res(i,j,k) - alphaa*lz_V(i,j,k,1)
     &                        - beta_s*lz_V(i,j,k,2)
            res(i,j,k) = temp
            beta  = beta + temp**2
c           else
c           temp  = res(i,j,k) - alpha1*lz_V(i,j,k,1)
c    &                        - beta1_s*lz_V(i,j,k,2)
c           res(i,j,k) = temp
c           beta1 = beta1 + temp**2
            end if
         end do; end do; end do
!$acc wait
         if (nproc.gt.1) then
            beta_v(1)=beta; beta_v(2)=beta1
            call MPI_ALLREDUCE(MPI_IN_PLACE,beta_v,2,MPI_RPREC,MPI_SUM,
     &           COMM_TINKER,ierr)
            beta=beta_v(1); beta1=beta_v(2)
         end if
         beta = sqrt(beta); beta1=sqrt(beta1)
         nrmb = beta

         if (makeOrth.and.EigenVec_l.and.kk.gt.lossOrthIt) then  ! FULL ORTHOGONALISATION
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
            nrmb  = nrm_b
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

!$acc parallel loop async collapse(2) default(present)
         do k = 1,npoleloc; do i=1,3
            res(i,2,k) = res(i,1,k)
         end do; end do

         if (tinkerdebug.gt.0.and.rank.eq.0) then
 13      format(A,I8,2F14.7,D15.7,2x,'northo(',I0,')')
         print 13, 'it lanzcos',kk,alpha,beta,l2_norm,notho
         end if

      end do

c     if (krylovDim.gt.120.and.makeOrth.and.EigenVec_l) then
c        call OrthogonalizeBase1(lz_Q,3*npoleloc,3*npoleloc
c    &       ,krylovDim,120)
c     end if
      if (tinkerdebug.gt.0) then
 14   format(A,2I4)
      if (kk.ne.krylovDim+1) write(*,14) 'early termination of',
     &   'lanzcos ',kk,'k0(',krylovDim,')'
         if (EigenVec_l) then
            call chkOrthogonality(lz_Q,3*npoleloc,krylovDim
     &          ,3*npoleloc,'lz_Q')
         end if
      end if

      end subroutine

      subroutine chkOrthogonality_(Mat,nl,nc,ldM,name)
      use domdec ,only: nproc
      implicit none
      integer,intent(in):: nc,nl,ldM
      real(t_p),intent(in):: Mat(ldM,*)
      character(8) name
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
!$acc serial async
      if (num.gt.0) then
         print*,"IssueOrth ",num,real(scal_min,4),real(scal_max,4)
     &         ,real(scal_zer,4)
      end if
!$acc end serial
!$acc end data
      end subroutine

      subroutine chkOrthogonality_r(Mat,nl,nc,ldM)
      use domdec ,only: nproc
      implicit none
      integer,intent(in):: nc,nl,ldM
      real(r_p),intent(in):: Mat(ldM,*)
      integer i,j,k,ii,vl,num
      real(r_p) scal,scal_max,scal_zer,scal_min,eps,sol
      parameter(vl=512,eps=1d-8)

      !TODO Add parallel version
      if (nproc.gt.1) return

      num=0
      scal_max=-huge(scal_max)
      scal_min= huge(scal_min)
      scal_zer= 1d6
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
            end if
         end if; end do
      end do; end do;
!$acc serial async
      if (num.gt.0) then
         print*,"Issue in chkOrth",num,scal_min,scal_zer,scal_max
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
      real(r_p) W(krylovDim), work(3*krylovDim)
      real(r_p) AB(ldAB,krylovDim), Z(krylovDim,krylovDim)
      real(t_p),pointer:: Z_(:,:)
      real(r_p) nrm,nrm1,scal,scal1,scal2

      if (deb_path) print*, ' getEigenVec'

      call prmem_request (lz_WorkSpace,krylovDim**2)
      Z_(1:krylovDim,1:krylovDim) => lz_WorkSpace(1:krylovDim**2)
      do i = 1,krylovDim
         AB(1,i) = lz_T(i)
         AB(2,i) = lz_T(krylovDim+i)
      end do

      call M_sbev('V','L',krylovDim,kd,AB,ldAB
     &           ,EigVal,Z,krylovDim,work,info)
      Z_(:,:) = Z(:,:)

      if (info.ne.0) then
         print*, ' M_sbev fails with error',info 
      end if

      if (.not.EigenVec_l) return

 12   format('eig val',50F10.6,/)
 13   format('EigV',I4,' nrm2(',2F9.5,') scal(',2F11.7,')',2I3,F11.7)

      !if(rank.eq.0) write(*,12) W(1:50)
!$acc update device(EigVal,Z_) async

#ifdef _OPENACC
      lda = 3*npoleloc
!$acc host_data use_device(lz_Q,Z_,Vp)
      info= M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                  ,lda,nEigenV,krylovDim
     &                  ,alpha,lz_Q,lda,Z_,krylovDim
     &                  ,beta ,  Vp,lda)
!$acc end host_data
      if (info.ne.CUBLAS_STATUS_SUCCESS) then
         print '(A,A,I4,A,I0)'
     &        ,' M_cublasGemm ',__FILE__,__LINE__
     &        ,' fails with error ',info
      end if

      if (nproc.eq.1) then
      call chkOrthogonality(Vp(1,1,1),3*npoleloc,nEigenV
     &     ,3*npoleloc,'Vp')
      end if

      if (krylovdim.gt.100) then
         call OrthogonalizeBase(Vp(1,1,1),3*npoleloc,3*npoleloc,nEigenV)
      end if

      call calcEigenVectorImage

      if (tinkerdebug.gt.0) then
         if (nproc.eq.1) then
         call chkOrthogonality(Vp(1,1,1),3*npoleloc,nEigenV
     &        ,3*npoleloc,'Vp')
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
      use interfaces ,only: cuPOTRI,cuPOTRF,cuGESV
#if TINKER_MIXED_PREC
     &                ,cuPOTRIm
#endif
      use mpi
      use mpole      ,only: npoleloc,npolebloc
      use orthogonalize
      use utilgpu    ,only: rec_stream
      use sizes      ,only: tinkerdebug
      implicit none
      integer i,j,k,it,it1
      real(t_p),pointer:: work(:,:,:,:)
      real(t_p) eps
      integer sd,sdec,ndec,ndec1,stdec,ledec,dim1
      real(r_p) tot1,tot2,tmp
      real(t_p) tot
      parameter(sdec=256)
      interface
      subroutine MatxVec(Rin,Rout)
      real(t_p),contiguous:: Rin(:,:,:),Rout(:,:,:)
      end subroutine
      end interface

      if (deb_path) print*,'  calcEigenVectorImage'
      lz_V(1:3,1:2,1:npolebloc,1:2) => lz_WorkSpace(1:12*npolebloc)

      do it = 1,nEigenV
         eps = -huge(eps)
!$acc parallel loop async default(present)
         do k = 1,npoleloc; do j = 1,2; do i = 1,3
            lz_V(i,j,k,1) = Vp(i,k,it)
         end do; end do; end do
         call Matxvec(lz_V(:,:,:,1),lz_V(:,:,:,2))
!$acc parallel loop async default(present)
         do k = 1,npoleloc; do i = 1,3
            TVp(i,k,it) = lz_V(i,1,k,2)
            !eps  = max(eps,abs(lz_V(i,1,k,1)-lz_V(i,2,k,2)))
         end do; end do
c!$acc wait
c         print*, 'diff Tmat vec',it,eps
      end do

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
         !stdec = (sd-1)*sdec+1
         !ledec = min(sd*sdec,dim1)
!$acc loop vector
         do k = 1,dim1
            tot1 = tot1 + Vp(k,1,it1)*TVp(k,1,it)
         end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!!$acc atomic
            MatE(it1,it) = MatE(it1,it) + tot1
            MatE(it,it1) = MatE(it,it1) + tot1
         end if; end do;

      end do; end do

      if (nproc.gt.1) then
!$acc wait
      call MPI_ALLREDUCE(MPI_IN_PLACE,MatE,nEigenV**2,MPI_RPREC
     &                  ,MPI_SUM,COMM_TINKER,it)
      end if

      ! Casting operation
!$acc parallel loop async default(present)
      do i=1,nEigenV**2
         MatE_(i,1) = MatE(i,1)
      end do

      ! MatE Cholesky Inversion
#ifdef _OPENACC
      ! Cholesky Factorisation & Inversion
!$acc host_data use_device(MatE_)
      call cuPOTRF(nEigenV,MatE_,nEigenV,rec_stream)
      call cuPOTRI(nEigenV,MatE_,nEigenV,rec_stream)
!$acc end host_data

       ! LU Decomposition & Inversion
c!$acc host_data use_device(MatE_,Ipivot,MIdk)
c      call cuGESV(nEigenV,nEigenV,MatE_,nEigenV,Ipivot,MIdk,nEigenV
c     &           ,rec_stream)
c!$acc parallel loop async deviceptr(MatE_,MIdk)
c      do i = 1,nEigenV**2
c         MatE_(i,1) = MIdk(i,1)
c      end do
c!$acc end host_data
#endif

#if TINKER_DOUBLE_PREC
      if (tinkerdebug.gt.0) then
      tot1=1.0; tot2=0.0
!$acc host_data use_device(MatE,MatE_,MIdk)
      it = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                  ,nEigenV,nEigenV,nEigenV
     &                  ,tot1,MatE,nEigenV,MatE_,nEigenV
     &                  ,tot2,MIdk,nEigenV)
!$acc end host_data
      if (it.ne.CUBLAS_STATUS_SUCCESS) then
         print '(A,A,I4,A,I0)'
     &         ,' M_cublasGemm ',__FILE__,__LINE__
     &         ,' fails with error ',it
      end if
!$acc wait
!$acc update host(MatE_,MatE,MIdk)
c     call printVec(" MatE  ",MatE ,nEigenV**2, prec=12, precd=8,
c    &      period=min(nEigenV,20))
c     call printVec(" MatE_ ",MatE_,nEigenV**2, prec=12, precd=8,
c    &      period=min(nEigenV,20))
      call chkOrthogonality(MIdk,nEigenV,nEigenV,nEigenV,'nEigenV')
      end if
#endif

      end subroutine
c
c     Compute projection of Vin in the deflated subspace
c
      subroutine projectorOpe(Vin,Vout,invert)
      use orthogonalize
      use mpole
      use inform
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)
      integer  ,intent(in)   :: invert

      integer ii,i,j,k,dim1
      integer sd,sdec,ndec,ndec1,stdec,ledec
      real(r_p) tot1,tot2,tmp
      real(t_p) tot
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

      if (invert.eq.0) then
!$acc kernels async default(present)
      Evec_(1:2*nEigenV,1) = Evec(1:2*nEigenV,1)
!$acc end kernels
      else
!$acc parallel loop collapse(2) present(Evec_,Evec,EigVal) async
      do j = 1,2; do i = 1,nEigenV
         Evec_(i,j) = Evec(i,j)/EigVal(i)
      end do; end do
      end if

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0
!$acc loop seq
         do ii = 1,nEigenV
            tot = tot + Vp(i,k,ii)*Evec_(ii,j)
         end do
         Vout(i,j,k) = Vin(i,j,k) - tot
      end do; end do; end do;

      !call minmaxone(Vin,6*npoleloc,'Vout')

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
      real(t_p) tot
      real(t_p),pointer::work(:,:,:)
      parameter(sdec=64)

      dim1  = npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)
      work(1:3,1:2,1:npoleloc)=>lz_WorkSpace(1:6*npoleloc)

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

!$acc host_data use_device(MatE_,Evec_,MIdk)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                   ,nEigenV,2,nEigenV
     &                   ,1.0_ti_p,MatE_,nEigenV,Evec_,nEigenV
     &                   ,0.0_ti_p,MIdk,nEigenV)
!$acc end host_data

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
      real(t_p) tot
      real(t_p),pointer::work(:,:,:)
      parameter(sdec=64)

      dim1  = npoleloc
      ndec1 = dim1/sdec
      ndec  = merge(ndec1+1,ndec1,mod(dim1,sdec).gt.0)
      work(1:3,1:2,1:npoleloc)=>lz_WorkSpace(1:6*npoleloc)

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
            tot1 = tot1 + TVp(i,k,ii)*Vin(i,1,k)
            tot2 = tot2 + TVp(i,k,ii)*Vin(i,2,k)
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

!$acc host_data use_device(MatE_,Evec_,MIdk)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_T,CUBLAS_OP_N
     &                   ,nEigenV,2,nEigenV
     &                   ,1.0_ti_p,MatE_,nEigenV,Evec_,nEigenV
     &                   ,0.0_ti_p,MIdk,nEigenV)
!$acc end host_data

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0
!$acc loop seq
         do ii = 1,nEigenV
            tot = tot + Vp(i,k,ii)*MIdk(ii,j)
         end do
         Vout(i,j,k) = Vin(i,j,k) - tot
      end do; end do; end do;
      !call minmaxone(work,6*npoleloc,'work')

      end subroutine
c
c     Compute Vout = alpha*Vout + Vp*EigVal^{-1}*Vp^{T}
c
      subroutine invertDefaultQmat (Vin,Vout,alpha)
      use inform  ,only: minmaxone
      use orthogonalize
      use mpole
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)
      real(t_p),intent(in)   :: alpha
      integer ii,i,j,k,dim1
      integer sd,sdec,ndec,ndec1,stdec,ledec
      real(r_p) tot1,tot2,tmp
      real(t_p) tot
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
            tot1 = tot1 + Vp(i,k,ii)*Vin(i,1,k)
            tot2 = tot2 + Vp(i,k,ii)*Vin(i,2,k)
         end do; end do

!$acc loop vector
         do k = 1,sdec; if (k.eq.1) then
!$acc atomic
            Evec(ii,1) = Evec(ii,1) + tot1
!$acc atomic
            Evec(ii,2) = Evec(ii,2) + tot2
         end if; end do

      end do; end do

!$acc parallel loop collapse(2) present(Evec_,Evec,EigVal) async
      do j = 1,2; do i = 1,nEigenV
         Evec_(i,j) = Evec(i,j)/EigVal(i)
      end do; end do

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0
!$acc loop seq
         do ii = 1,nEigenV
            tot = tot + Vp(i,k,ii)*Evec_(ii,j) 
         end do
         Vout(i,j,k) = alpha*Vout(i,j,k) + tot
      end do; end do; end do;

      end subroutine

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
      real(t_p) tot
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
            tot1 = tot1 + Vp(i,k,ii)*Vin(i,1,k)
            tot2 = tot2 + Vp(i,k,ii)*Vin(i,2,k)
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

!$acc kernels async default(present)
      Evec_(1:2*nEigenV,1) = Evec(1:2*nEigenV,1)
!$acc end kernels

!$acc host_data use_device(MatE_,Evec_,MIdk)
      ierr = M_cublasGemm(cBhandle,CUBLAS_OP_N,CUBLAS_OP_N
     &                   ,nEigenV,2,nEigenV
     &                   ,1.0_ti_p,MatE_,nEigenV,Evec_,nEigenV
     &                   ,0.0_ti_p,MIdk,nEigenV)
!$acc end host_data

!$acc parallel loop collapse(3) default(present) async
      do k = 1,npoleloc; do j = 1,2; do i = 1,3;
         tot = 0
!$acc loop seq
         do ii = 1,nEigenV
            tot = tot + Vp(i,k,ii)*MIdk(ii,j) 
         end do
         Vout(i,j,k) = alpha*Vout(i,j,k) + tot
      end do; end do; end do;

      end subroutine

      subroutine diagPrecnd(Vin,Vout,Minv)
      use mpole
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)
      real(t_p),intent(in)   :: Minv(npoleloc)
      integer i,j,k

!$acc parallel loop collapse(3) async default(present)
      do k=1,npoleloc; do j=1,2; do i=1,3
         Vout(i,j,k) = Minv(k)*Vin(i,j,k)
      end do; end do; end do

      end subroutine

      subroutine adaptedDeflat2 (Vin,Vout,Minv)
      use inform  ,only: minmaxone
      use orthogonalize
      use mpole
      use tinheader
      implicit none
      real(t_p),intent(in)   :: Vin (3,2,*)
      real(t_p),intent(inout):: Vout(3,2,*)
      real(t_p),intent(in)   :: Minv(npoleloc)
      integer i
      real(t_p),pointer:: work(:,:,:)
      work(1:3,1:2,1:npoleloc)=> lz_WorkSpace(6*npoleloc+1:12*npoleloc)

      print*, 'adaptedDeflat2'
      call diagPrecnd(Vin,work,Minv)
      call minmaxone(work,6*npoleloc,'work')
      call ProjectorPTransGeneric(work,work)
      call minmaxone(work,6*npoleloc,'work')
      call ApplyQxV(Vin,work,1.0_ti_p)
      call minmaxone(work,6*npoleloc,'work')
      print*, '------------------------------------------------'

!$acc parallel loop async default(present)
      do i=1,6*npoleloc
         Vout(i,1,1) = work(i,1,1)
      end do

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
      call bissectionGivens1(lz_T,lz_T(sz_max+1),sz_init
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

