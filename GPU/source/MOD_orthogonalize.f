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
c     lz_StoreSpace*  holds the sotorage space for module subroutine
c     lz_T        holds tridiagonal matrix computed from lanzcos iterations
c
#include "tinker_precision.h"
c
      module orthogonalize
#ifdef _OPENACC
      use cublas_v2
#endif
      integer,parameter:: orthMaxInstances=2
      integer krylovDim,nEigenV
      integer betw
      logical EigenVec_l
      integer  ,target:: Ipivot(150)
      real(t_p),allocatable,target:: lz_WorkSpace(:)
      real(t_p),allocatable,target:: lz_StoreSpace(:),lz_StoreSpace1(:)
      real(r_p),allocatable,target:: lz_StoreSpace_rp(:)
     &         , lz_StoreSpace1_rp(:)
      real(t_p),allocatable,target:: Evec_(:,:), MIdk(:,:)
      real(r_p),allocatable,target:: Evec (:,:)

      real(r_p),pointer :: lz_T(:), EigVal(:)
      real(t_p),pointer :: lz_V(:,:,:,:)
      real(t_p),pointer :: lz_Q(:,:,:),Vp(:,:,:),TVp(:,:,:)
      real(r_p),pointer :: MatE(:,:)
      real(t_p),pointer :: MatE_(:,:)

      parameter(betw=32)

#ifdef _OPENACC
      type(cublashandle):: cBhandle
      logical :: init_cublas_handle_l=.false.
#endif
!$acc declare create(Ipivot)

      interface printVec
         module procedure printVec1_
#if TINKER_MIXED_PREC
         module procedure printVec1_r
#endif
      end interface

      interface chkOrthogonality
        subroutine chkOrthogonality_(Mat,nl,nc,ldM,name)
        implicit none
        integer,intent(in):: nc,nl,ldM
        real(t_p),intent(in):: Mat(ldM,*)
        character(*) name
!DIR$ ignore_tkr (r) Mat
        end subroutine
#if TINKER_MIXED_PREC
        subroutine chkOrthogonality_r(Mat,nl,nc,ldM)
        implicit none
        integer,intent(in):: nc,nl,ldM
        real(r_p),intent(in):: Mat(ldM,*)
!DIR$ ignore_tkr (r) Mat
        end subroutine
#endif
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
      end do; end do

      end subroutine

      end module
