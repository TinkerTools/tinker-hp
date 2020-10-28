c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################################
c     ##                                                                             ##
c     ##  module interfaces  -- Tinker-HP interfaces functions and subroutines       ##
c     ##                                                                             ##
c     #################################################################################
c
#include "tinker_precision.h"
      module interfaces
      use iso_c_binding ,only: c_int,c_float,c_double
#ifdef _OPENACC
      use cudafor ,only:cuda_stream_kind
#endif
      enum, bind(C) !:: typelist
         enumerator list_verlet,list_block_block,list_block_atoms
      end enum
      enum ,bind(C)
         enumerator conf_ehal1    !0
         enumerator conf_ehal2    !1
         enumerator conf_ehalcu   !2
      end enum
      enum, bind(C)
         enumerator conf_elj1
         enumerator conf_elj1gpu
         enumerator conf_elj1cu
      end enum
      enum, bind(C)
         enumerator conf_ecreal1d
         enumerator conf_ecreal1dgpu
         enumerator conf_ecreal1d_cu
      end enum
      enum, bind(C) !:: proc_emreal1c
         enumerator conf_emreal1c_1     !0
         enumerator conf_emreal1c_2     !1
         enumerator conf_emreal1c_pre   !2  ( Precompute Routine )
         enumerator conf_emreal1c_4     !3  ( CUDA C )
         enumerator conf_emreal1c_5     !4  ( CUDA Fortran )
      end enum
      enum, bind(C) !:: proc_elfd
         enumerator conf_efld0_directgpu1
         enumerator conf_efld0_directgpu2
         enumerator conf_efld0_directgpu3
      end enum
      enum, bind(C)
         enumerator conf_tmatxb_ref
         enumerator conf_tmatxb_c1
         enumerator conf_tmatxb_pre
         enumerator conf_tmatxb_c2
         enumerator conf_tmatxb_c3
      end enum
      enum, bind(C) !:: proc_epreal1c
         enumerator conf_epreal1c_1
         enumerator conf_epreal1c_2
         enumerator conf_epreal1c_3
      end enum
      enum, bind(C)
         enumerator conf_fphi_acc
         enumerator conf_fphi_cuda
      end enum
      enum, bind(C)
         enumerator conf_grid_acc
         enumerator conf_grid_cuda
      end enum
      enum, bind(C) !:: config
         enumerator conf_any         !0
         enumerator conf_vlist       !1
         enumerator conf_mlist       !2
         enumerator conf_clist       !3
         enumerator ::conf_ehal=4    !4
         enumerator conf_elj         !5
         enumerator conf_echg        !6
         enumerator ::conf_mpole=8   !8
         enumerator conf_efld0       !9
         enumerator conf_tmat        !a
         enumerator conf_polar       !b
         enumerator conf_fphi        !c
         enumerator conf_grid        !d
      end enum

      enum, bind(C)
         enumerator PrDynamic
         enumerator PrAnalyze
         enumerator PrMinimize
         enumerator PrTestgrad
      end enum

#if (defined(SINGLE) | defined(MIXED))
                                !fedcba9876543210!
      integer(8) :: sub_config=Z'0011242402220000'
#else
                                !fedcba9876543210!
      integer(8) :: sub_config=Z'0000131101110000'
#endif
      ! parameter for long range interactions comput
      integer,parameter :: short_mode=0
      ! parameter for short range interactions comput
      integer,parameter :: long_mode=1

      integer ProgramID


!  #############################################################################
!  SECTION
!                Neighbor list routines Interfaces
!
!  #############################################################################

      interface
        subroutine build_cell_list2(ineignl,cell_order,
     &              nlocnl,buf,bx_cell,by_cell,bz_cell)
          integer,intent(in):: nlocnl
          integer,intent(in):: ineignl(:)
          integer,intent(inout)::cell_order(:)
          integer,intent(inout)::bx_cell,by_cell,bz_cell
          real(t_p),intent(in)::buf
        end subroutine
      end interface

      interface
         subroutine reorder_nblist(nblist,nneig,nneig_cut,nloc,cutoff2,
     &                             ktype,glob)
            integer,intent(in) :: nloc
            integer,intent(inout) :: nblist(:,:)
            integer,intent(out) :: nneig_cut(:)
            integer,intent(in) :: nneig(:)
            integer,intent(in) :: glob(:),ktype(:)
            real(t_p),intent(in) ::  cutoff2
         end subroutine
      end interface

      interface
        subroutine set_ElecData_CellOrder(rebuild_nl)
        logical,intent(in):: rebuild_nl
        end subroutine
      end interface

      interface
        subroutine pre_process_adjacency_matrix
     &             (nblock,nlocnlb_pair)
        integer,intent(in):: nblock
        integer,intent(out):: nlocnlb_pair
        end subroutine
      end interface

      procedure(tinker_void_sub) ,pointer:: vlist_block_p
     &           ,mlist_block_p,clist_block_p
     &           ,vlistcell_p,mlistcell_p,clistcell_p



!  #############################################################################
!  SECTION
!                Ehal routines Interfaces
!
!  #############################################################################
      interface
        subroutine ehal1c_correct_scaling(xred,yred,zred,
     &             vxx,vxy,vxz,vyy,vyz,vzz)
          real(t_p),intent(in):: xred(:)
          real(t_p),intent(in):: yred(:)
          real(t_p),intent(in):: zred(:)
          real(r_p)  vxx,vxy,vxz,vyy,vyz,vzz
        end subroutine
        subroutine ehal1c_correct_scaling_shortlong(xred,yred,zred,
     &             vxx,vxy,vxz,vyy,vyz,vzz,mode)
          real(t_p),intent(in):: xred(:)
          real(t_p),intent(in):: yred(:)
          real(t_p),intent(in):: zred(:)
          real(r_p)  vxx,vxy,vxz,vyy,vyz,vzz
          character*10,intent(in):: mode
        end subroutine
        subroutine ehal3c_correct_scaling(xred,yred,zred)
          real(t_p),intent(in):: xred(:)
          real(t_p),intent(in):: yred(:)
          real(t_p),intent(in):: zred(:)
        end subroutine
        subroutine ehalshortlong3c_correct_scaling(xred,yred,zred,mode)
          real(t_p),intent(in):: xred(:)
          real(t_p),intent(in):: yred(:)
          real(t_p),intent(in):: zred(:)
          character*10,intent(in):: mode
        end subroutine
      end interface
      interface
        subroutine ehal3c_cu
        end subroutine
        subroutine ehalshortlong3c_cu
        end subroutine
      end interface

      ! Lennard-Jones subroutines
      interface 
        subroutine elj1c
        end subroutine
        subroutine elj1cgpu
        end subroutine
        subroutine elj1c_cu
        end subroutine
        subroutine eljshortlong1cgpu
        end subroutine
        subroutine elj3c
        end subroutine
        subroutine elj3cgpu
        end subroutine
        subroutine elj3c_cu
        end subroutine
      end interface
      interface
        subroutine elj1_scaling(xred,yred,zred,
     &             vxx,vxy,vxz,vyy,vyz,vzz)
        real(t_p),intent(in):: xred(:)
        real(t_p),intent(in):: yred(:)
        real(t_p),intent(in):: zred(:)
        real(r_p)  vxx,vxy,vxz
        real(r_p)  vyy,vyz,vzz
        end subroutine
        subroutine elj1shortlong_scaling(xred,yred,zred,
     &             vxx,vxy,vxz,vyy,vyz,vzz)
        real(t_p),intent(in):: xred(:)
        real(t_p),intent(in):: yred(:)
        real(t_p),intent(in):: zred(:)
        real(r_p)  vxx,vxy,vxz
        real(r_p)  vyy,vyz,vzz
        end subroutine
        subroutine elj3_scaling(xred,yred,zred)
        real(t_p),intent(in):: xred(:)
        real(t_p),intent(in):: yred(:)
        real(t_p),intent(in):: zred(:)
        end subroutine
      end interface

      procedure(tinker_void_sub):: ehal1cgpu1,ehal1cgpu2,ehal1c_cu
     &                           , ehalshort1c,ehal1c,ehallong1c
     &                           , ehalshort1cgpu,ehalshortlong1c_cu
     &                           , ehallong1cgpu
     &                           , ehal3c,ehal3cgpu,ehalshortlong3cgpu
      procedure(tinker_void_sub) ,pointer:: ehal1c_p,ehalshort1c_p
     &                           ,ehallong1c_p,ehal3c_p
     &                           ,ehalshortlong3c_p
      procedure(elj1cgpu)        ,pointer::elj1c_p
      procedure(elj1cgpu)        ,pointer::eljsl1c_p
      procedure(elj3cgpu)        ,pointer::elj3c_p



!  #############################################################################
!  SECTION
!                Tmatxb routines Interfaces
!
!  #############################################################################
      interface
        subroutine tmatxb_pmevec(nrhs,dodiag,mu,ef)
           integer  ,intent(in) :: nrhs
           logical  ,intent(in) :: dodiag
           real(t_p)            :: mu(3,2,*)
           real(t_p)            :: ef(3,2,*)
        end subroutine
        subroutine tmatxb_pmegpu(nrhs,dodiag,mu,ef)
           integer  ,intent(in) :: nrhs
           logical  ,intent(in) :: dodiag
           real(t_p),intent(in) :: mu(:,:,:)
           real(t_p),intent(out):: ef(:,:,:)
        end subroutine
        subroutine tmatxb_pmegpu1(nrhs,dodiag,mu,ef)
           integer  ,intent(in) :: nrhs
           logical  ,intent(in) :: dodiag
           real(t_p),intent(in) :: mu(:,:,:)
           real(t_p),intent(out):: ef(:,:,:)
        end subroutine
        subroutine tmatxb_pme_core1(mu,efi)
           real(t_p),intent(in)   ::  mu(:,:,:)
           real(t_p),intent(inout):: efi(:,:,:)
        end subroutine
        subroutine tmatxb_pme_core2(mu,efi)
           real(t_p),intent(in)   ::  mu(:,:,:)
           real(t_p),intent(inout):: efi(:,:,:)
        end subroutine
        subroutine otf_dc_tmatxb_pme_core2(mu,efi)
           real(t_p) ,intent(in)   :: mu (:,:,:)
           real(t_p) ,intent(inout):: efi(:,:,:)
        end subroutine
        subroutine otf_dc_tmatxb_pme_core3(mu,efi)
           real(t_p) ,intent(in)   :: mu (:,:,:)
           real(t_p) ,intent(inout):: efi(:,:,:)
        end subroutine
        subroutine tmatxb_pme_core3(mu,efi)
           real(t_p),intent(in)   ::  mu(:,:,:)
           real(t_p),intent(inout):: efi(:,:,:)
        end subroutine
        subroutine tmatxb_correct_interactions(mu,efi)
           real(t_p),intent(in)   ::  mu(:,:,:)
           real(t_p),intent(inout):: efi(:,:,:)
        end subroutine
        subroutine otf_dc_tmatxb_correct_interactions(mu,efi)
           real(t_p) ,intent(in) :: mu (:,:,:)
           real(t_p) ,intent(out):: efi(:,:,:)
        end subroutine
      end interface
      procedure(tmatxb_pmegpu)   ,pointer:: tmatxb_p
      procedure(tmatxb_pme_core1),pointer:: tmatxb_pme_core_p
      procedure(otf_dc_tmatxb_pme_core2),pointer:: 
     &          otf_dc_tmatxb_pme_core_p



!  #############################################################################
!  SECTION
!                FFt Grid routines Interfaces
!
!  #############################################################################
      interface fphi_uind_sitegpu
         subroutine fphi_uind_sitegpu1(fphi1,fphi2,fphi_sum)
            real(t_p),intent(out),dimension(:,:)::fphi1,fphi2
            real(t_p),intent(out),dimension(:,:)::fphi_sum
         end subroutine
         subroutine fphi_uind_sitegpu2(fphi1,fphi2)
            real(t_p),intent(out),dimension(:,:)::fphi1,fphi2
         end subroutine
      end interface
      interface fphi_uind_sitecu
         subroutine fphi_uind_sitecu2(fphi1,fphi2)
            real(t_p),intent(out),dimension(:,:)::fphi1,fphi2
         end subroutine
         subroutine fphi_uind_sitecu1(fphi1,fphi2,fphi_sum)
            real(t_p),intent(out),dimension(:,:)::fphi1,fphi2
            real(t_p),intent(out),dimension(:,:)::fphi_sum
         end subroutine
      end interface
      interface
         subroutine grid_uind_sitegpu(fuindvec,fuinpvec,qgrid2loc)
           real(t_p),intent(in),dimension(3,*) :: fuindvec,fuinpvec
           real(t_p) qgrid2loc(:,:,:,:,:)
         end subroutine
         subroutine grid_uind_sitecu(fuindvec,fuinpvec,qgrid2loc)
           real(t_p),intent(in),dimension(3,*) :: fuindvec,fuinpvec
           real(t_p) qgrid2loc(:,:,:,:,:)
         end subroutine
      end interface
      interface
         subroutine grid_mpole_sitegpu(fmpvec)
         real(t_p) fmpvec(10,*)
         end subroutine
         subroutine grid_mpole_sitecu(fmpvec)
         real(t_p) fmpvec(10,*)
         end subroutine
      end interface
      interface
         subroutine grid_pchg_sitegpu
         end subroutine
         subroutine grid_pchg_sitecu
         end subroutine
      end interface
      interface
         subroutine grid_pchg_force
         end subroutine
         subroutine grid_pchg_forcecu
         end subroutine
      end interface
      interface
         subroutine bspline_fill_sitegpu(config)
         integer,intent(in),optional::config
         end subroutine
      end interface
      procedure(tinker_void_sub) :: fphi_mpole_sitegpu,fphi_mpole_sitecu
      procedure(fphi_uind_sitegpu2),pointer:: fphi_uind_site2_p
      procedure(fphi_uind_sitegpu1),pointer:: fphi_uind_site1_p
      procedure(grid_uind_sitegpu) ,pointer:: grid_uind_site_p
      procedure(grid_mpole_sitegpu),pointer:: grid_mpole_site_p
      procedure(grid_pchg_sitegpu) ,pointer:: grid_pchg_site_p
      procedure(grid_pchg_force)   ,pointer:: grid_pchg_force_p
      procedure(tinker_void_sub)   ,pointer:: fphi_mpole_site_p


!  #############################################################################
!  SECTION
!                Charge routines Interfaces
!
!  #############################################################################
      interface
        subroutine ecreal1d
        end subroutine
        subroutine ecreal1dgpu
        end subroutine
#ifdef _CUDA
        subroutine ecreal1d_cu
        end subroutine
#endif
        subroutine ecrealshortlong1dgpu
        end subroutine
        subroutine ecreal_scaling
        end subroutine
      end interface

      interface
        subroutine ecreal3c
        end subroutine
        subroutine ecreal3dgpu
        end subroutine
#ifdef _CUDA
        subroutine ecreal3d_cu
        end subroutine
#endif
        subroutine ecreal3_scaling
        end subroutine
      end interface

      procedure(ecreal1dgpu),pointer:: ecreal1d_p
     &                      ,ecrealshortlong1d_p
      procedure(ecreal3dgpu),pointer:: ecreal3d_p


!  #############################################################################
!  SECTION
!                Empole routines Interfaces
!
!  #############################################################################
      interface
        subroutine commpoleglob(polebuff)
          integer,intent(inout)::polebuff(:)
        end subroutine
      end interface

      interface
        subroutine
     &  emreal_correct_interactions(tem,vxx,vxy,vxz,vyy,vyz,vzz)
           real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz
           real(t_p),intent(inout):: tem(:,:)
        end subroutine
        subroutine emreal_correct_interactions_shortlong
     &                          (tem,vxx,vxy,vxz,vyy,vyz,vzz)
           real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz
           real(t_p),intent(inout):: tem(:,:)
        end subroutine
        subroutine emreal3_correct_interactions
        end subroutine
        subroutine emreal3_correct_interactions_shortlong
        end subroutine
      end interface
      interface
        subroutine
     &  emreal1c_core1(tem,vxx,vxy,vxz,vyy,vyz,vzz)
           real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz
           real(t_p),intent(inout):: tem(3,*)
        end subroutine
      end interface
      interface
        subroutine emreal3d
        end subroutine
        subroutine emrealshort3d
        end subroutine
        subroutine emreallong3d
        end subroutine
        subroutine emreal3dgpu
        end subroutine
        subroutine emreal3d_cu
        end subroutine
        subroutine emrealshortlong3d
        end subroutine
        subroutine emrealshortlong3d_cu
        end subroutine
      end interface
      procedure(emreal1c_core1):: emreal1c_core2,emreal1c_core3
     &                           ,emreal1c_core4,emreal1c_core5
     &                           ,emrealshortlong1c_core
     &                           ,emrealshortlong1c_core2
      procedure(emreal1c_core1)  ,pointer:: emreal1c_core_p
     &                           ,emrealshortlong1c_core_p
      procedure(emreal3dgpu),pointer:: emreal3d_p
     &                      ,emrealshortlong3d_p



!  #############################################################################
!  SECTION
!               Induced dipoles and polarisation solver routines Interfaces
!
!  #############################################################################
      interface
        subroutine efld0_directgpu2(nrhs,ef)
           integer  ,intent(in)   :: nrhs
           real(t_p),intent(inout):: ef(:,:,:)
        end subroutine
        subroutine otf_dc_efld0_directgpu2(nrhs,ef)
           integer  ,intent(in)   :: nrhs
           real(t_p),intent(inout):: ef(:,:,:)
        end subroutine
        subroutine otf_dc_efld0_directgpu3(nrhs,ef)
           integer  ,intent(in)   :: nrhs
           real(t_p),intent(inout):: ef(:,:,:)
        end subroutine
        subroutine efld0_direct_correct_scaling(ef)
           real(t_p),intent(inout):: ef(:,:,:)
        end subroutine
      end interface

      interface
        subroutine inducepcg_pme2gpu(matvec,nrhs,precnd,ef,mu,murec)
           import tmatxb_pmegpu
           integer  ,intent(in)   :: nrhs
           logical  ,intent(in)   :: precnd
           real(t_p),intent(in)   :: ef   (:,:,:)
           real(t_p),intent(inout):: mu   (:,:,:)
           real(t_p),intent(inout):: murec(:,:,:)
           procedure(tmatxb_pmegpu)::matvec
        end subroutine
        subroutine inducepcg_pmegpu(matvec,nrhs,precnd,ef,mu,murec)
           import tmatxb_pmegpu
           integer  ,intent(in)   :: nrhs
           logical  ,intent(in)   :: precnd
           real(t_p),intent(in)   :: ef   (:,:,:)
           real(t_p),intent(inout):: mu   (:,:,:)
           real(t_p),intent(inout):: murec(:,:,:)
           procedure(tmatxb_pmegpu)::matvec
        end subroutine
        subroutine inducejac_pme2gpu(matvec,nrhs,dodiis,ef,mu,murec)
           import tmatxb_pmegpu
           integer  ,intent(in)   :: nrhs
           logical  ,intent(in)   :: dodiis
           real(t_p),intent(in)   :: ef   (:,:,:)
           real(t_p),intent(inout):: mu   (:,:,:)
           real(t_p),intent(inout):: murec(:,:,:)
           procedure(tmatxb_pmegpu)::matvec
        end subroutine
        subroutine inducejac_pmegpu(matvec,nrhs,dodiis,ef,mu,murec)
           import tmatxb_pmegpu
           integer  ,intent(in)   :: nrhs
           logical  ,intent(in)   :: dodiis
           real(t_p),intent(in)   :: ef   (:,:,:)
           real(t_p),intent(inout):: mu   (:,:,:)
           real(t_p),intent(inout):: murec(:,:,:)
           procedure(tmatxb_pmegpu)::matvec
        end subroutine
        subroutine inducepcg_shortrealgpu(matvec,nrhs,precnd,ef,mu)
           import tmatxb_pmegpu
           integer  ,intent(in)   :: nrhs
           logical  ,intent(in)   :: precnd
           real(t_p),intent(in)   :: ef   (:,:,:)
           real(t_p),intent(inout):: mu   (:,:,:)
           procedure(tmatxb_pmegpu)::matvec
        end subroutine
      end interface

!  #############################################################################
!  SECTION
!               Epolar routines Interfaces
!
!  #############################################################################
      interface
        subroutine epreal1c_core1(trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
           real(t_p) trqvec(3,*)
           real(r_p) vxx,vxy,vxz,vyy,vyz,vzz
        end subroutine
      end interface

      interface
        subroutine epreal1c_correct_scale
     &           (trqvec,vxx,vxy,vxz,vyy,vyz,vzz)
           real(t_p),intent(inout):: trqvec(:,:)
           real(r_p),intent(inout):: vxx,vxy,vxz,vyy,vyz,vzz
        end subroutine
        subroutine epreal3c_correct_scale
        end subroutine
      end interface

      interface
        subroutine epreal3d
        end subroutine
        subroutine epreal3dgpu
        end subroutine
        subroutine epreal3d_cu
        end subroutine
      end interface
      procedure(epreal1c_core1)::epreal1c_core2,epreal1c_core3
      procedure(epreal1c_core1),pointer:: epreal1c_core_p
      procedure(epreal3dgpu)   ,pointer:: epreal3d_p


!  #############################################################################
!  SECTION
!              Tinker CUDA(C/C++) Wrapper routines Interfaces
!
!  #############################################################################
#ifdef _CUDA
      interface
        subroutine cu_filter_lst_sparse
     &     (cell_glob,cell_scan,xred,yred,zred,matb_lst,
     &      nvdwlocnlb,nblock,matsize,nvdwlocnlb_pair,vbuf2,
     &      vblst,cell_len,xcell,ycell,zcell,xcell2,ycell2,zcell2,
     &      rec_stream)
     &       bind(C,name="filter_lst_sparse")
           import cuda_stream_kind
           integer,device:: cell_glob(*),cell_scan(*),matb_lst(*)
           integer,device:: cell_len(*),vblst(*)
           real(t_p),device:: xred(*),yred(*),zred(*)
           integer,value::nvdwlocnlb,nblock,matsize,nvdwlocnlb_pair
           integer(cuda_stream_kind),value:: rec_stream
           real(t_p),value::vbuf2
           real(t_p),value::xcell,ycell,zcell,xcell2,ycell2,zcell2
        end subroutine
      end interface

      interface 
        subroutine cu_tmatxb_pme 
     &            (ipole,pglob,ploc,ieblst,eblst,x,y,z
     &            ,pdamp,thole,polarity,mu,efi
     &            ,npolelocnlb,npolelocnlb_pair,npolebloc,n,nproc
     &            ,cut2,alsq2,alsq2n,aewald
     &            , xcell, ycell, zcell,xcell2,ycell2,zcell2
     &            ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &            ,stream)
     &            bind(C,name="cu_tmatxb_pme")
        import cuda_stream_kind
        integer,value,intent(in)::npolelocnlb,npolebloc,n,nproc
     &         ,npolelocnlb_pair
        integer(cuda_stream_kind),value::stream
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,xcell,xcell2,ycell,ycell2,zcell,zcell2
     &           ,cut2,alsq2,alsq2n,aewald
        integer,device::ipole(*),pglob(*),ploc(*),ieblst(*),eblst(*)
        real(t_p),device::pdamp(*),thole(*),polarity(*),x(*),y(*),z(*)
        real(t_p),device:: mu(6,*)
        real(t_p),device:: efi(6,*)
        end subroutine
      end interface

      interface 
        subroutine cu_efld0_direct
     &            (ipole,pglob,ploc,ieblst,eblst,x,y,z
     &            ,pdamp,thole,polarity,rpole,efi
     &            ,npolelocnlb,npolelocnlb_pair,npolebloc,n,nproc
     &            ,cut2,alsq2,alsq2n,aewald
     &            , xcell, ycell, zcell,xcell2,ycell2,zcell2
     &            ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &            ,stream)
     &            bind(C,name="cu_efld0_direct")
        import cuda_stream_kind
        integer,value,intent(in)::npolelocnlb,npolebloc,n,nproc
     &         ,npolelocnlb_pair
        integer(cuda_stream_kind),value::stream
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,xcell,xcell2,ycell,ycell2,zcell,zcell2
     &           ,cut2,alsq2,alsq2n,aewald
        integer,device::ipole(*),pglob(*),ploc(*),ieblst(*),eblst(*)
        real(t_p),device::pdamp(*),thole(*),polarity(*),x(*),y(*),z(*)
     &           ,rpole(13,*)
        real(t_p),device:: efi(6,*)
        end subroutine
        subroutine cu_otfdc_efld0_direct
     &            (ipole,pglob,ploc,grplst,atmofst,npergrp,kofst
     &            ,ieblst,eblst,x,y,z,pdamp,thole,polarity,rpole
     &            ,efi,zmat
     &            ,npolelocnlb,npolelocnlb_pair,npolebloc,n,nproc
     &            ,cut2,alsq2,alsq2n,aewald
     &            , xcell, ycell, zcell,xcell2,ycell2,zcell2
     &            ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &            ,stream)
     &            bind(C,name="cu_otfdc_efld0_direct")
        import cuda_stream_kind
        integer,value,intent(in)::npolelocnlb,npolebloc,n,nproc
     &         ,npolelocnlb_pair
        integer(cuda_stream_kind),value::stream
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,xcell,xcell2,ycell,ycell2,zcell,zcell2
     &           ,cut2,alsq2,alsq2n,aewald
        integer,device::ipole(*),pglob(*),ploc(*),ieblst(*),eblst(*)
     &           ,grplst(*),atmofst(*),npergrp(*),kofst(*)
        real(t_p),device::pdamp(*),thole(*),polarity(*),x(*),y(*),z(*)
     &           ,rpole(13,*)
        real(t_p),device:: efi(6,*),zmat(*)
        end subroutine
      end interface

      interface
        subroutine cu_emreal1c
     &        ( ipole, pglob, loc, ieblst, eblst, x, y, z, rpole
     &        , dem, tem, em_buffer, vir_buffer
     &        , npolelocnlb, npolelocnlb_pair, npolebloc, n
     &        , off2, f, alsq2, alsq2n, aewald
     &        , xcell, ycell, zcell,xcell2,ycell2,zcell2
     &        ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &        , stream)
     &        bind(C,name="cu_emreal1c")
        import cuda_stream_kind
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair
        integer(cuda_stream_kind),value::stream
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,xcell,xcell2,ycell,ycell2,zcell,zcell2
     &           ,off2,alsq2,alsq2n,aewald,f
        integer,device::ipole(*),pglob(*),loc(*),ieblst(*),eblst(*)
        real(t_p),device:: x(*),y(*),z(*),rpole(13,*)
        real(t_p),device:: tem(3,*),vir_buffer(*)
        real(r_p),device:: dem(3,*),em_buffer(*)
        end subroutine
      end interface 

      interface
        subroutine C_init_env(devicenum,nproc,rank,debug)
     &             bind(C,name="C_init_env")
        integer,intent(in),value::devicenum,nproc,rank,debug
        end subroutine
        subroutine C_get_cell(xcell,ycell,zcell,eps_cell,octahedron
     &                       ,box34) bind(C,name="C_get_cell")
        real(t_p),intent(in),value::xcell,ycell,zcell,eps_cell
     &                      ,box34
        integer ,intent(in),value::octahedron
        end subroutine

        ! Interfaces to cuSolver wrapper
        subroutine initcuSolver( stream )
     &             bind(C,name="initcuSolverHandle")
        import cuda_stream_kind
        integer(cuda_stream_kind),value:: stream
        end subroutine
        subroutine cuPOTRF(n, A, lda, stream)
     &             bind(C,name="cuPOTRF_Wrapper")
        import cuda_stream_kind,c_int
        integer(c_int),value :: n,lda
        real(t_p)     ,device:: A(*)
        integer(cuda_stream_kind),value:: stream
        end subroutine
        subroutine cuPOTRS(n, A, lda, B, ldb, stream)
     &             bind(C,name="cuPOTRS_Wrapper")
        import cuda_stream_kind,c_int
        integer(c_int),value :: n,lda,ldb
        real(t_p)     ,device:: A(*),B(*)
        integer(cuda_stream_kind),value:: stream
        end subroutine
        subroutine cuGESV(n, nrhs, A, lda, Ipiv, B, ldb, stream)
     &             bind(C,name="cuGESV_Wrapper")
        import cuda_stream_kind,c_int
        integer(c_int),value :: n,lda,ldb,nrhs
        integer(c_int),device:: Ipiv(*)
        real(t_p)     ,device:: A(*),B(*)
        integer(cuda_stream_kind),value:: stream
        end subroutine
        subroutine destroycuSolver()
     &             bind(C,name="destroycuSolverHandle")
        end subroutine
      end interface
#endif



!  #############################################################################
!  SECTION
!                Scaling factor routines Interfaces
!
!  #############################################################################

      interface 
         subroutine vdw_scaling
         end subroutine
      end interface
      interface 
         subroutine vdw_scalingnl
         end subroutine
      end interface
      interface 
         subroutine mpole_scaling
         end subroutine
      end interface
      interface 
         subroutine polar_scaling
         end subroutine
      end interface
      interface 
         subroutine polar_scaling1
         end subroutine
      end interface

      procedure(polar_scaling),pointer:: polar_scaling_p

!  #############################################################################
!  SECTION
!                Torque routines Interfaces
!
!  #############################################################################
      interface torquegpu
         subroutine torquedgpu(trqvec,frx,fry,frz,de,extract)
            real(t_p),intent(in   ),dimension(:,:)::trqvec
            real(r_p),intent(inout),dimension(:,:)::de
            real(t_p),intent(out  ),dimension(:,:)::frx,fry,frz
            logical*1,intent(in   )               ::extract
         end subroutine
         subroutine torquergpu(n,pole,poleloc,trqvec,de,queue)
            integer  ,intent(in   )::n,queue
            integer  ,intent(in   )::pole(:),poleloc(:)
            real(t_p),intent(in   )::trqvec(:,:)
            real(r_p),intent(inout)::de(:,:)
         end subroutine
      end interface


      procedure(tinker_void_sub) :: vlist_block,mlist_block
     &           ,vlistcell,mlistcell,clistcell


      procedure(efld0_directgpu2):: efld0_directgpu,efld0_directgpu3
      procedure(otf_dc_efld0_directgpu2)::oft_dc_efld0_directgpu3

      procedure(efld0_directgpu2),pointer:: efld0_directgpu_p
      procedure(otf_dc_efld0_directgpu2),pointer:: 
     &          otf_dc_efld0_directgpu_p

!  #############################################################################
!  SECTION
!               Minimize routines Interfaces
!
!  #############################################################################
      interface
        subroutine optsave (ncycle,xx)
          integer ncycle
          real(r_p) xx(*)
        end subroutine
        subroutine lbfgs (n,x0,minimum,grdmin,fgvalue,optsave)
          integer n
          real(r_p) x0(*)
          real(r_p) minimum
          real(r_p) grdmin
          real(r_p) fgvalue
          external  fgvalue,optsave
        end subroutine
        subroutine search (n,f,g,x,p,f_move,angle,ncalls,
     &                            fgvalue,status)
          integer n,ncalls
          real(r_p) fgvalue
          real(r_p) f
          real(r_p) f_move
          real(r_p) angle
          real(r_p) x(*)
          real(r_p) g(*)
          real(r_p) p(*)
          character*9 status
          external fgvalue
        end subroutine
      end interface


      contains

      subroutine tinker_void_sub
      ! Do nothing
      end subroutine

      end module
