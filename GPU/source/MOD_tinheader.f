c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module tinker header  --  General Stuff for Global Handle    ##
c     ##                                                               ##
c     ###################################################################
c
c
c     tp  define by the preprocessor as tinker-HP parameter for floating presicion
c     rp  define by the preprocessor as tinker-HP parameter for mixed precision
c     MPI_TYPE  define by the preprocessor as tinker-HP parameter for MPI_exchange
c     ti_eps    contains an espilon value for Tinker
c     prec_eps  smallest value depending on the precision such that (2+prec_esp>2)
c     prec1_eps  smallest value depending on the precision such that (1+prec_esp>1)

#include "tinker_precision.h"
      module tinheader
      !use mpi,only:MPI_REAL4,MPI_REAL8
      implicit none
      integer ti_p,re_p
      !integer MPI_TYPE,MPI_RTYPE
      real(t_p) ti_eps,prec_eps,prec1_eps
      real(r_p) precm_eps

      parameter(ti_p=t_p)
      parameter(re_p=r_p)
      parameter(prec_eps =2*epsilon(ti_eps))
      parameter(prec1_eps=epsilon(ti_eps))
      parameter(precm_eps=2*epsilon(precm_eps))
      !parameter(MPI_TYPE=MPI_TPREC)
      !parameter(MPI_RTYPE=MPI_RPREC)

#if (defined(SINGLE) || defined(MIXED))
      parameter(ti_eps=1d-5)
#else
      parameter(ti_eps=1d-12)
#endif

      end module

      module tinTypes
c
c     Scalar derive type definition
c
      type real3
         real(t_p) x,y,z
      end type real3
      type real3_red
         real(r_p) x,y,z
      end type
      type real6
         real(t_p) x,y,z
         real(t_p) xx,yy,zz
      end type real6
      type real6_red
         real(t_p) x,y,z
         real(t_p) xx,yy,zz
      end type
      type cross
         real(t_p) xx,xy,xz,yy,yz,zz
      end type
      type cross_red
         real(r_p) xx,xy,xz,yy,yz,zz
      end type
      type iframe
         integer ist,ien
         integer jst,jen
         integer kst,ken
      end type iframe
      type rframe
         real(t_p) xst,xen
         real(t_p) yst,yen
         real(t_p) zst,zen
      end type rframe
      type int3
         integer i,j,k
      end type int3
      type rpole_elt 
        real(t_p) c,dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz
      end type

      type IntPointerContainer
        integer,pointer:: pel(:)
      end type
      type Int2dPointerContainer
        integer,pointer:: pel(:,:)
      end type
      type RealPointerContainer
        real(t_p),pointer:: pel(:)
      end type
      type Real2dPointerContainer
        real(t_p),pointer:: pel(:,:)
      end type
#ifdef _CUDA
      type IntDevPointerContainer
        integer,device,pointer:: pel(:)
      end type
      type Int2dDevPointerContainer
        integer,device,pointer:: pel(:,:)
      end type
      type LogicalDevPointerContainer
        logical,device,pointer:: pel(:)
      end type
      type Logical1DevPointerContainer
        logical(1),device,pointer:: pel(:)
      end type
      type Char8DevPointerContainer
        character(8),device,pointer:: pel(:)
      end type
      type RealDevPointerContainer
        real(t_p),device,pointer:: pel(:)
      end type
      type Real2dDevPointerContainer
        real(t_p),device,pointer:: pel(:,:)
      end type
      type Real3dDevPointerContainer
        real(t_p),device,pointer:: pel(:,:,:)
      end type
#endif

      interface assignment (=)
         module procedure assign_real_to_real3
         module procedure assign_real_to_real3_red
      end interface

      contains 

      subroutine assign_real_to_real3( left, right )
!$acc routine
      type(real3) ,intent(out):: left
      real(t_p)   ,intent(in) :: right

      left%x = right
      left%y = right
      left%z = right
      end subroutine
      subroutine assign_real_to_real3_red( left, right )
!$acc routine
      type(real3_red),intent(out):: left
      real(r_p)      ,intent(in) :: right

      left%x = right
      left%y = right
      left%z = right
      end subroutine

      end module
