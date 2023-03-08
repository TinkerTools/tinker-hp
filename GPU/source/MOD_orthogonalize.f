c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################################
c     ##                                                                                ##
c     ##  module orthogonalize  --  gather routines and data for                        ##
c     ##                            matrix orthogonalization and projection operators   ##
c     ##                                                                                ##
c     ####################################################################################
c
c     krylovDim   holds the krylov subspace dimension
c     nEigenV     holds the number of EigenVector to be computed
c     EigenVec_l  computes the EigenVectors to Tmat ? yes : no
c     lz_WorkSpace   holds the workspace for lanzcos iterations
c     lz_SS[*][_]    holds the storage space for module subroutine -- (*)instances  (_)inc precision
c     lz_T[*]        holds tridiagonal matrix computed from lanzcos iterations (*)dipoles
c
#include "tinker_macro.h"
c
      module orthogonalize
#ifdef _OPENACC
      use cublas_v2
#endif
      integer,parameter:: orthMaxInstances=2
      integer krylovDim,maxKD,wKD,nEigenV
      integer betw,l_idx
      integer MdiagITmat,OrthBase
      logical EigenVec_l,saveResAmat_l
      logical debg
      integer  ,allocatable,target:: Ipivot(:)
      real(t_p),allocatable,target:: lz_WorkSpace(:)
      real(t_p),allocatable,target:: lz_SS(:),lz_SS1(:)
      real(r_p),allocatable,target:: lz_SS_(:)
     &         , lz_SS1_(:), lz_WS_(:)
      real(t_p),allocatable,target:: Evec_(:,:), MIdk(:,:)
      real(r_p),allocatable,target:: Evec (:,:)

      real(r_p),pointer :: lz_T(:),lz_T1(:),lz_MT(:,:),lz_MT1(:,:)
     &         , EigVal(:), Eigval1(:)
      real(t_p),pointer :: lz_V(:,:,:,:)
      real(t_p),pointer :: lz_Q(:,:,:), Vp(:,:,:), TVp(:,:,:)
     &         , lz_AQ(:,:,:),lz_Q1(:,:,:),Vp1(:,:,:),TVp1(:,:,:)
     &         , lz_AQ1(:,:,:)
     &         , lz_Z(:,:),lz_Z1(:,:)
      real(r_p),pointer :: MatE (:,:),MatE1 (:,:)
      real(t_p),pointer :: MatE_(:,:),MatE1_(:,:)

      parameter(betw=32,debg=.false.)

#ifdef _OPENACC
      type(cublashandle):: cBhandle
      logical :: init_cublas_handle_l=.false.
#endif

      interface printVec
         module procedure printVec1_
         module procedure printVec1_2d
#if TINKER_MIXED_PREC
         module procedure printVec1_r
         module procedure printVec1_r_2d
#endif
      end interface

      interface chkOrthogonality
        subroutine chkOrthogonality_r4(Mat,nl,nc,ldM,name)
        implicit none
        integer,intent(in):: nc,nl,ldM
        real(4),intent(in):: Mat(ldM,*)
        character(8),intent(in) :: name
!DIR$ ignore_tkr (r) Mat
        end subroutine
        subroutine chkOrthogonality_r8(Mat,nl,nc,ldM,name)
        implicit none
        integer,intent(in):: nc,nl,ldM
        real(8),intent(in):: Mat(ldM,*)
        character(8),intent(in) :: name
!DIR$ ignore_tkr (r) Mat
        end subroutine
      end interface

      interface
      real(r_p) function pol_c(Diag, SubDiag, i, lambda) result(pol_val)
      real(r_p),intent(in) :: Diag(*),SubDiag(*)
      integer  ,intent(in) :: i
      real(r_p),intent(in) :: lambda
      end function
      end interface

      contains

      subroutine getMaxr(array,n,val,ind)
      implicit none
      real(r_p) array(*)
      integer,intent(in)::n
      real(r_p),intent(out):: val
      integer,intent(out):: ind
      integer i

      val = array(1)
      ind = 1
      do i = 2,n
         if (array(i).gt.val) then
            val = array(i)
            ind = i
         end if
      end do

      end subroutine

      subroutine printVec1_(name,Vec,siz,prec,precD,period)
      implicit none
      character(*)         :: name
      real(t_p),intent(in) :: Vec(*)
      integer,intent(in)   :: siz
      integer,optional     :: prec,precD,period
      integer i,begi,last,prec_,precd_,period_, nloop
      character*15::fmt,fmt1

      prec_=10; precd_=6; period_=10;
      if (present(prec)  ) prec_  = prec
      if (present(precD) ) precd_ = precD
      if (present(period)) period_= period
      nloop = ceiling(siz/(period_*1.0d0))

 11   format("(6X,"I0,"I",I0,")")
 12   format("(I6,",I0,"F",I0,".",I0,")")
      write(fmt1,11) period_,prec_
      write(fmt,12) period_,prec_,precd_

      write(*,*) name
      write(*,fmt1) [ (i,i=1,period_) ]
      do i = 1,nloop
         begi = (i-1)*period_+1
         last = min(i*period_,siz)
         write(*,fmt) i,Vec(begi:last)
      end do

      end subroutine

      subroutine printVec1_2d(name,Vec,siz,prec,precD,period)
      implicit none
      character(*)         :: name
      real(t_p),intent(in) :: Vec(1,*)
      integer,intent(in)   :: siz
      integer,optional     :: prec,precD,period
      integer i,begi,last,prec_,precd_,period_, nloop
      character*15::fmt,fmt1

      prec_=10; precd_=6; period_=10;
      if (present(prec)  ) prec_  = prec
      if (present(precD) ) precd_ = precD
      if (present(period)) period_= period
      nloop = ceiling(siz/(period_*1.0d0))

 11   format("(6X,"I0,"I",I0,")")
 12   format("(I6,",I0,"F",I0,".",I0,")")
      write(fmt1,11) period_,prec_
      write(fmt,12) period_,prec_,precd_

      write(*,*) name
      write(*,fmt1) [ (i,i=1,period_) ]
      do i = 1,nloop
         begi = (i-1)*period_+1
         last = min(i*period_,siz)
         write(*,fmt) i,Vec(1,begi:last)
      end do

      end subroutine

      subroutine printVec1_r(name,Vec,siz,prec,precD,period)
      implicit none
      character(*)         :: name
      real(r_p),intent(in) :: Vec(*)
      integer,intent(in)   :: siz
      integer,optional     :: prec,precD,period
      integer i,begi,last,prec_,precd_,period_, nloop
      character*15::fmt,fmt1

      prec_=10; precd_=6; period_=10;
      if (present(prec)  ) prec_  = prec
      if (present(precD) ) precd_ = precD
      if (present(period)) period_= period
      nloop = ceiling(siz/(period_*1.0d0))

 11   format("(6X,"I0,"I",I0,")")
 12   format("(I6,",I0,"F",I0,".",I0,")")
      write(fmt1,11) period_,prec_
      write(fmt,12) period_,prec_,precd_

      write(*,*) name
      write(*,fmt1) [ (i,i=1,period_) ]
      do i = 1,nloop
         begi = (i-1)*period_+1
         last = min(i*period_,siz)
         write(*,fmt) i,Vec(begi:last)
      end do

      end subroutine

      subroutine printVec1_r_2d(name,Vec,siz,prec,precD,period)
      implicit none
      character(*)         :: name
      real(r_p),intent(in) :: Vec(1,*)
      integer,intent(in)   :: siz
      integer,optional     :: prec,precD,period
      integer i,begi,last,prec_,precd_,period_, nloop
      character*15::fmt,fmt1

      prec_=10; precd_=6; period_=10;
      if (present(prec)  ) prec_  = prec
      if (present(precD) ) precd_ = precD
      if (present(period)) period_= period
      nloop = ceiling(siz/(period_*1.0d0))

 11   format("(6X,"I0,"I",I0,")")
 12   format("(I6,",I0,"F",I0,".",I0,")")
      write(fmt1,11) period_,prec_
      write(fmt,12) period_,prec_,precd_

      write(*,*) name
      write(*,fmt1) [ (i,i=1,period_) ]
      do i = 1,nloop
         begi = (i-1)*period_+1
         last = min(i*period_,siz)
         write(*,fmt) i,Vec(1,begi:last)
      end do

      end subroutine

      subroutine setMIdk
      implicit none
      integer i,j

!$acc parallel loop collapse(2) async default(present)
      do i = 1,nEigenV; do j = 1,nEigenV
         if (i.eq.j) then
            MIdk(i,i) = 1
         else
            MIdk(j,i) = 0
         end if
         MatE (j,i) = 0
         MatE_(j,i) = 0
         MatE1 (j,i) = 0
         MatE1_(j,i) = 0
      end do; end do

      end subroutine

      end module
