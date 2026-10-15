!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module tinker header  --  General Stuff for Global Handle    ##
!     ##                                                               ##
!     ###################################################################
!
!
!     tp  define by the preprocessor as tinker-HP parameter for floating presicion
!     rp  define by the preprocessor as tinker-HP parameter for mixed precision
!     MPI_TYPE  define by the preprocessor as tinker-HP parameter for MPI_exchange
!     ti_eps    contains an espilon value for Tinker
!     prec_eps  smallest value depending on the precision such that (2+prec_esp>2)
!     prec1_eps  smallest value depending on the precision such that (1+prec_esp>1)

#include "tinker_macro.h"
module tinheader
   !use mpi,only:MPI_REAL4,MPI_REAL8
   implicit none
   integer    ti_p,re_p
   integer    i_init
   integer(1) zero1,one1,two1,three1
   integer(4) zeroi4,onei4,twoi4
   integer(8) zero8,one8
   logical    isrel_build
   !integer MPI_TYPE,MPI_RTYPE
   real(t_p)  ti_eps,prec_eps,prec1_eps
   real(r_p)  precm_eps
   real(t_p)  a_init
   real(t_p)  zeror,oner,twor,halfr
   real(r_p)  zerom,onem,twom,halfm
   mdyn_rtyp  zeromd,onemd,twomd
   real(t_p)  Infinity
   real(r_p)  Infinity_r

   parameter (&
      &zero1=0,one1=1,two1=2,three1=3&
      &,zeroi4=0,onei4=1,twoi4=2&
      &,zero8=0,one8=1&
      &,zeror=0,oner=1.0,twor=2.0,halfr=0.5&
      &,zerom=0,onem=1.0,twom=2.0,halfm=0.5&
      &,zeromd=0,onemd=1.0,twomd=2.0&
#if (defined(SINGLE) || defined(MIXED))
      &,ti_eps=1d-5&
#else
      &,ti_eps=1d-12&
#endif
      &,ti_p  =t_p&
      &,re_p  =r_p&
      &,Infinity=4*huge(0.0_ti_p)&
      &,Infinity_r=4*huge(0.0_re_p)&
      &,i_init=-1&
      &,a_init=-1.0&
      &,prec_eps =2*epsilon(zeror)&
      &,prec1_eps=  epsilon(zeror)&
      &,precm_eps=2*epsilon(zerom)&
#if TINKERHP_REL_BUILD
      &,isrel_build=.true.&
#else
      &,isrel_build=.false.&
#endif
      &)
   !parameter(MPI_TYPE=MPI_TPREC)
   !parameter(MPI_RTYPE=MPI_RPREC)

end module

module tinTypes
!
!     Scalar derive type definition
!
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

   type dcdinfo_t
      integer(kind=4) :: idcd,nframes,istart,iend,nevery
      integer(kind=8) :: nframes_pos,iend_pos,curr_pos
      integer(kind=8) :: filesize, framesize
      integer(kind=4) :: natoms
      real(kind=4) :: timestep
      character(len=80), allocatable :: titlesdcd(:)
      character(:), allocatable :: filename
   end type dcdinfo_t

   type trajectory_t
      ! Copy policy: noncopyable; use the lifecycle procedures only.
      integer :: backend = 0
      integer :: frame_count = 0
      logical :: is_active = .false.
      logical :: append_mode = .false.
      logical :: owns_dcd_file = .false.
      character(len=240) :: filename = ''
      type(dcdinfo_t) :: dcdinfo
   end type trajectory_t


   ! Copy policy: noncopyable; owns bead workspaces and trajectory lifecycles.
   TYPE POLYMER_COMM_TYPE
      integer :: nbeads
      integer :: nbeadsproc
      integer,allocatable :: nbeadsloc(:)
      integer,allocatable :: ibead_beg(:)
      integer,allocatable :: ibead_end(:)
      integer,allocatable :: iproc(:)
      integer :: index

      logical :: allocated=.FALSE.

      real(r_p), allocatable :: pos(:,:,:),vel(:,:,:)
      real(r_p), allocatable :: eigpos(:,:,:),eigvel(:,:,:)
      real(r_p), allocatable :: forces(:,:,:)
      real(r_p), allocatable :: eigforces(:,:,:)
      real(r_p), allocatable :: forces_slow(:,:,:)

      real(r_p), allocatable :: epot(:)
      real(r_p), allocatable :: vir(:,:,:)

      real(r_p), allocatable :: dip(:,:)
      real(r_p), allocatable :: dipind(:,:)

      real(r_p) :: gyr

      type(dcdinfo_t),allocatable :: dcdinfo(:)
      type(trajectory_t),allocatable :: dcd_trajectory(:)
   END TYPE

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

end module
