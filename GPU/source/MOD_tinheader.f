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
#include "tinker_types.h"
      module tinheader
      !use mpi,only:MPI_REAL4,MPI_REAL8
      implicit none
      integer ti_p,re_p
      integer i_init
      integer(1) zero1,one1,two1,three1
      integer(4) zeroi4,onei4,twoi4
      integer(8) zero8,one8
      !integer MPI_TYPE,MPI_RTYPE
      real(t_p) ti_eps,prec_eps,prec1_eps
      real(r_p) precm_eps
      real(t_p) a_init
      real(t_p)  zeror,oner,twor,halfr
      real(r_p)  zerom,onem,twom,halfm
      mdyn_rtyp  zeromd,onemd,twomd

      parameter (
     &           zero1=0,one1=1,two1=2,three1=3
     &          ,zeroi4=0,onei4=1,twoi4=2
     &          ,zero8=0,one8=1
     &          ,zeror=0,oner=1.0,twor=2.0,halfr=0.5
     &          ,zerom=0,onem=1.0,twom=2.0,halfm=0.5
     &          ,zeromd=0,onemd=1.0,twomd=2.0
     &          ,ti_p=t_p
     &          ,re_p=r_p
     &          ,i_init=-1
     &          ,a_init=-1.0
     &          ,prec_eps =2*epsilon(ti_eps)
     &          ,prec1_eps=  epsilon(ti_eps)
     &          ,precm_eps=2*epsilon(precm_eps)
     &          )
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
      type mdyn3_r
         mdyn_rtyp x,y,z
      end type
      type real6
         real(t_p) x,y,z
         real(t_p) xx,yy,zz
      end type real6
      type real7
       real(t_p) x,y,z
       real(t_p) xx,yy,zz,pa
      end type
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
