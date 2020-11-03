c
c     sorbonne university
c     washington university in saint louis
c     university of texas at austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1cu  --  buffered 14-7 energy & derivatives##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der waals energy and
c     its first derivatives with respect to cartesian coordinates
c
c     *_t texture variable destined to be attached to their target
#ifdef _CUDA
#define TINKER_CUF
#include  "tinker_precision.h"
      module ehal1cu
        use cudafor
        use tinheader,only: ti_p
        use sizes    ,only: maxclass
        use utilcu   ,only: nproc,all_lanes,VDW_BLOCK_DIM
        use utilgpu  ,only: BLOCK_SIZE,RED_BUFF_SIZE

        type(dim3) :: gridDim,blockDim
        integer  ,pointer,device::ired_t(:),cellv_glob_t(:)
     &           ,cellv_loc_t(:),loc_ired_t(:),vblst_t(:),ivblst_t(:)
     &           ,jvdw_t(:)
        real(t_p),pointer,texture::radmin_t(:,:),epsilon_t(:,:)
        real(t_p),pointer,device::
     &            kred_t(:),xred_t(:),yred_t(:),zred_t(:)

        contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "switch_respa.f.inc"
#include "pair_ehal1.f.inc"

        attributes(global) 
     &  subroutine ehal1_cu(xred,yred,zred,cellv_glob
     &           ,cellv_loc,loc_ired,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,dev,ev_buff,vir_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
     &           ,scexp,vlambda,scalpha,mut
#ifdef TINKER_DEBUG
     &           ,inter
#endif
     &                 )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
     &           ,cut2,cut,ghal,dhal,off2,off
     &           ,scexp,vlambda,scalpha
        logical  ,device :: mut(n)
        real(r_p),device :: ev_buff(RED_BUFF_SIZE)
        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb)
     &           ,vblst(nvdwlocnlb_pair*2),cellv_loc(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
        real(r_p),device,intent(inout)::dev(3,nbloc)
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
#endif

        integer ithread,iwarp,nwarp,ilane,istat,srclane
        integer dstlane,klane
        integer idx,kdx,kdx_,ii,j,i,iglob,it
        integer kglob,kbis,kt,kt_,kvloc
#ifdef TINKER_DEBUG
        integer kglob_
#endif
c       integer,shared:: ninte
        real(t_p) e
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz
        real(r_p) gxi,gyi,gzi
        real(r_p) gxk,gyk,gzk
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block
        logical muti
        logical,shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
c       ninte   = 0

c       if (ithread.eq.1) print*,'ehal1c_cu_in'
c    &   ,blockDim%x,gridDim%x,nwarp,nvdwlocnlb_pair,n
c    &   ,c0,c1,c2,c3,c3,cut2,ghal,dhal
c    &   ,ev,vxx,vxy,vxz,vyy

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Set Data to compute to zero
           ev_    = 0
           gxi    = 0
           gyi    = 0
           gzi    = 0
           gxk    = 0
           gyk    = 0
           gzk    = 0
           vxx_   = 0
           vxy_   = 0
           vxz_   = 0
           vyy_   = 0
           vyz_   = 0
           vzz_   = 0

           ! Load atom block i parameters
           idx    = (vblst(2*ii+1)-1)*warpsize + ilane 
           iglob  = cellv_glob(idx)
           i      = cellv_loc (idx)
c          ivloc  = loc_ired(idx)
           it     = jvdw(idx)
           xi     = xred(idx)
           yi     = yred(idx)
           zi     = zred(idx)
c          redi   = merge (1.0_ti_p,kred(iglob),(i.eq.ivloc))
           muti   = mut(iglob)

           ! Load atom block k parameters
           kdx    = (vblst(2*ii+2)-1)*warpsize + ilane
           kglob  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
c          kvloc  = loc_ired(kdx)
           kt     = jvdw(kdx)
           xk     = xred(kdx)
           yk     = yred(kdx)
           zk     = zred(kdx)
           same_block = (idx.ne.kdx)
c          redk   = merge (1.0_ti_p,kred(kglob),(kbis.eq.kvloc))
           mutk(threadIdx%x) = mut(kglob)
           !call syncthreads()

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
              kdx_    = __shfl(kdx  ,srclane)
              kt_     = __shfl(kt   ,srclane)
#ifdef TINKER_DEBUG
              kglob_  = __shfl(kglob,srclane)
#endif

              xpos    = xi - __shfl(xk,srclane)
              ypos    = yi - __shfl(yk,srclane)
              zpos    = zi - __shfl(zk,srclane)

              call image_inl(xpos,ypos,zpos)

              dedx = 0.0; dedy = 0.0; dedz = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                       ,same_block)
              if (do_pair.and.kdx_<=nvdwlocnl.and.rik2<=off2) then

                 rv2   =  radmin (kt_,it)
                 eps2  = epsilon (kt_,it)

                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                         ,cut2,cut,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,muti,mutk(klane)
     &                         ,e,dedx,dedy,dedz)

                 ev_   = ev_ + e

                 vxx_  = vxx_  + xpos * dedx
                 vxy_  = vxy_  + ypos * dedx
                 vxz_  = vxz_  + zpos * dedx
                 vyy_  = vyy_  + ypos * dedy
                 vyz_  = vyz_  + zpos * dedy
                 vzz_  = vzz_  + zpos * dedz

#ifdef TINKER_DEBUG
                 if (iglob<kglob_) then
                    istat = Atomicadd(inter(iglob) ,1)
                 else
                    istat = Atomicadd(inter(kglob_),1)
                 end if
#endif

              end if

              dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

              ! Increment gradient
              gxi = gxi + real(dedx,r_p)
              gyi = gyi + real(dedy,r_p)
              gzi = gzi + real(dedz,r_p)

              gxk = gxk + real(__shfl(dedx,dstlane),r_p)
              gyk = gyk + real(__shfl(dedy,dstlane),r_p)
              gzk = gzk + real(__shfl(dedz,dstlane),r_p)

           end do

           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )
 
           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi)
              istat   = atomicAdd (dev(2,i),gyi)
              istat   = atomicAdd (dev(3,i),gzi)
           end if

           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk)
              istat   = atomicSub (dev(2,kbis),gyk)
              istat   = atomicSub (dev(3,kbis),gzk)
           end if

           ! Increment virial term of van der Waals
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )

        end do
c       call syncthreads
c       if (ithread.eq.1) then
c          print*,ev
c          print*,vxx,vzz,vxz
c          print*,dev(1,1),dev(2,1),dev(3,1)
c          print*,dev(1,2),dev(2,2),dev(3,2)
c       end if

        end subroutine

        attributes(global)
     &  subroutine ehal3_cu2
     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,ev_buff,nev_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
     &           ,scexp,vlambda,scalpha,mut
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &           )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
     &           ,cut2,cut,ghal,dhal,off2,off
     &           ,scexp,vlambda,scalpha
        logical  ,device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        real(r_p),device::  ev_buff(RED_BUFF_SIZE)
        integer  ,device:: nev_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb),ivblst(nvdwlocnlb_pair)
     &           ,vblst(nvdwlocnlb_pair*(BLOCK_SIZE))
     &           ,cellv_loc(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
#ifdef TINKER_DEBUG
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kglob,kbis,kt,kt_,ist,nev_
        real(t_p) e,istat
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk,yk,zk,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz,devx,devy,devz
        logical do_pair,same_block,accept_mid
        logical muti
        logical,shared:: mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Set Data to compute to zero
           ev_    = 0
           nev_   = 0

           ! Load atom block i parameters
           iblock = ivblst(ii+1)
           if (iblock.eq.0) cycle
           idx    = (iblock-1)*warpsize + ilane 
           iglob  = cellv_glob(idx)
           i      = cellv_loc (idx)
           it     = jvdw  (idx)
           xi     = xred  (idx)
           yi     = yred  (idx)
           zi     = zred  (idx)
           muti   = mut(iglob)

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt     = jvdw  (kdx)
           xk     = xred  (kdx)
           yk     = yred  (kdx)
           zk     = zred  (kdx)
           same_block = (idx.ne.kdx)
           mutk(threadIdx%x) = mut(kglob)

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
              kglob_  = __shfl(kglob,srclane)
#endif
              kt_     = __shfl(kt   ,srclane)

              if (nproc.gt.1) then
                 xk_   = __shfl(xk,srclane)
                 yk_   = __shfl(yk,srclane)
                 zk_   = __shfl(zk,srclane)
                 xpos  = xi - xk_
                 ypos  = yi - yk_
                 zpos  = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
                 call image_inl(xpos,ypos,zpos)
              end if

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

                 rv2   =  radmin_t (kt_,it)
                 eps2  = epsilon_t (kt_,it)

                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                         ,cut2,cut,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,muti,mutk(klane)
     &                         ,e,dedx,dedy,dedz)

                 ev_   = ev_ + e
                 nev_  = nev_ + 1

#ifdef TINKER_DEBUG
#endif
              end if

              dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

           end do

           kt_ = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(kt_) , ev_ )
           istat = atomicAdd(nev_buff(kt_) ,nev_ )
 
        end do
        end subroutine


        attributes(global)
     &  subroutine ehal1_cu2
     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,dev,ev_buff,vir_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
     &           ,scexp,vlambda,scalpha,mut
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &           ,inter,rank
#endif
     &           )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
     &           ,cut2,cut,ghal,dhal,off2,off
     &           ,scexp,vlambda,scalpha
        logical  ,device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        real(r_p),device :: ev_buff (RED_BUFF_SIZE)
        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb),ivblst(nvdwlocnlb_pair)
     &           ,vblst(nvdwlocnlb_pair*(BLOCK_SIZE))
     &           ,cellv_loc(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
        real(r_p),device,intent(inout)::dev(3,nbloc)
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
        integer,value :: rank
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kglob,kbis,ist
        integer ,shared,dimension(VDW_BLOCK_DIM):: kt
        real(t_p) e,istat
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
        real(t_p) dedx,dedy,dedz
        real(r_p),shared,dimension(VDW_BLOCK_DIM)::
     &            gxi,gyi,gzi,gxk,gyk,gzk
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid
        logical muti
        logical,shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.
c       ninte   = 0

c       if (ithread.eq.1) print*,'ehal1c_cu2_in'
c    &   ,blockDim%x,gridDim%x,nwarp,nvdwlocnlb_pair,n
c    &   ,c0,c1,c2,c3,c3,cut2,ghal,dhal
c    &   ,ev,vxx,vxy,vxz,vyy,nproc
c    &   ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Set Data to compute to zero
           ev_    = 0
           gxi(threadIdx%x)    = 0
           gyi(threadIdx%x)    = 0
           gzi(threadIdx%x)    = 0
           gxk(threadIdx%x)    = 0
           gyk(threadIdx%x)    = 0
           gzk(threadIdx%x)    = 0
           vxx_   = 0
           vxy_   = 0
           vxz_   = 0
           vyy_   = 0
           vyz_   = 0
           vzz_   = 0

           ! Load atom block i parameters
           iblock = ivblst(ii+1)
           if (iblock.eq.0) cycle
           idx    = (iblock-1)*warpsize + ilane 
           iglob  = cellv_glob(idx)
           i      = cellv_loc (idx)
           it     = jvdw  (idx)
           xi     = xred  (idx)
           yi     = yred  (idx)
           zi     = zred  (idx)
#ifndef TINKER_NO_MUTATE
           muti   = mut(iglob)
#endif

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt(threadIdx%x) = jvdw(kdx)
           xk(threadIdx%x) = xred(kdx)
           yk(threadIdx%x) = yred(kdx)
           zk(threadIdx%x) = zred(kdx)
           same_block = (idx.ne.kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob)
#endif

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
              kglob_  = __shfl(kglob,srclane)
#endif

              if (nproc.gt.1) then
                 xk_   = xk(klane)
                 yk_   = yk(klane)
                 zk_   = zk(klane)
                 xpos  = xi - xk_
                 ypos  = yi - yk_
                 zpos  = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 xpos    = xi - xk(klane)
                 ypos    = yi - yk(klane)
                 zpos    = zi - zk(klane)
                 call image_inl(xpos,ypos,zpos)
              end if

              dedx = 0.0; dedy = 0.0; dedz = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                       ,same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

                 rv2   =  radmin_t (kt(klane),it)
                 eps2  = epsilon_t (kt(klane),it)

                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                         ,cut2,cut,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,muti,mutk(klane)
     &                         ,e,dedx,dedy,dedz)

                 ev_   = ev_ + e

#ifdef TINKER_DEBUG
                 if (iglob<kglob_) then
                    ist = Atomicadd(inter(iglob),1)
                 else
                    ist = Atomicadd(inter(kglob_),1)
                 end if
c                if (rank.eq.0) then
c                if (iglob.eq.12) then
c                   if (iglob<kglob_) then
c                      print*,kglob_, kdx_,ii+1,iblock,j,ilane,'i'
c                   else
c                      print*, kglob_,real(-xpos,4),real(-ypos,4),-zpos
c                   end if
c                end if
c                if (kglob_.eq.12) then
c                   if (iglob<kglob_) then
c                      print*, iglob, kdx_,ii+1,iblock,j,ilane,'i'
c                   else
c                      print*, iglob,real(-xpos,4),real(-ypos,4),-zpos
c                   endif
c                end if
c                end if
#endif

                 vxx_  = vxx_  + xpos * dedx
                 vxy_  = vxy_  + ypos * dedx
                 vxz_  = vxz_  + zpos * dedx
                 vyy_  = vyy_  + ypos * dedy
                 vyz_  = vyz_  + zpos * dedy
                 vzz_  = vzz_  + zpos * dedz


                 !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

                 ! Accumulate interaction gradient
                 gxi(threadIdx%x) = gxi(threadIdx%x) + real(dedx,r_p)
                 gyi(threadIdx%x) = gyi(threadIdx%x) + real(dedy,r_p)
                 gzi(threadIdx%x) = gzi(threadIdx%x) + real(dedz,r_p)

                 gxk(klane) = gxk(klane) + real(dedx,r_p)
                 gyk(klane) = gyk(klane) + real(dedy,r_p)
                 gzk(klane) = gzk(klane) + real(dedz,r_p)

              end if

           end do

           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )
 
           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi(threadIdx%x))
              istat   = atomicAdd (dev(2,i),gyi(threadIdx%x))
              istat   = atomicAdd (dev(3,i),gzi(threadIdx%x))
           end if

           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
           end if

           ! Increment virial term of van der Waals
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )

        end do
        end subroutine

c        attributes(global)
c     &  subroutine ehal1_cu3
c     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
c     &           ,ivblst,vblst,jvdw,epsilon,radmin
c     &           ,ired,kred,dev,ev_buff,vir_buff
c     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
c     &           ,nvdwclass
c     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
c     &           ,scexp,vlambda,scalpha,mut
c     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
c#ifdef TINKER_DEBUG
c     &           ,inter,rank
c#endif
c     &           )
c
c        implicit none
c        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
c     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
c        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
c     &           ,cut2,cut,ghal,dhal,off2,off
c     &           ,scexp,vlambda,scalpha
c        logical  ,device :: mut(n)
c        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
c     &           ,p_zbeg,p_zend
c        real(r_p),device :: ev_buff (RED_BUFF_SIZE)
c        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
c        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
c     &           ,loc_ired(nvdwlocnlb),ivblst(nvdwlocnlb_pair)
c     &           ,vblst(nvdwlocnlb_pair*(BLOCK_SIZE))
c     &           ,cellv_loc(nvdwlocnlb)
c        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
c        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
c     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
c        real(t_p),device,intent(in):: xred(nvdwlocnlb)
c     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
c        real(r_p),device,intent(inout)::dev(3,nbloc)
c#ifdef TINKER_DEBUG
c        integer,device:: inter(*)
c        integer,value :: rank
c#endif
c
c        integer ithread,iwarp,nwarp,ilane,srclane
c        integer dstlane,klane,iilane,iblock
c        integer idx,kdx,kdx_,kglob_,ii,j,i,i1
c        integer ,shared,dimension(VDW_BLOCK_DIM):: iglob,it
c        integer kbis,ist
c        integer ,shared,dimension(VDW_BLOCK_DIM):: kt,kglob,ikstat
c        real(t_p) e,istat
c        real(r_p) ev_
c        real(t_p) xk_,yk_,zk_,xpos,ypos,zpos
c        real(t_p) rik2,rv2,eps2
c        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
c     &           ,xi,yi,zi
c        real(t_p) dedx,dedy,dedz
c        real(r_p),shared,dimension(VDW_BLOCK_DIM)::
c     &            gxi,gyi,gzi,gxk,gyk,gzk
c        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
c        logical do_pair,same_block,accept_mid
c        logical,shared::mutk(VDW_BLOCK_DIM)
c     &         ,muti(VDW_BLOCK_DIM)
c
c        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
c        iwarp   = (ithread-1) / warpsize
c        nwarp   = blockDim%x*gridDim%x / warpsize
c        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
c        accept_mid = .true.
cc       ninte   = 0
c
cc       if (ithread.eq.1) print*,'ehal1c_cu2_in'
cc    &   ,blockDim%x,gridDim%x,nwarp,nvdwlocnlb_pair,n
cc    &   ,c0,c1,c2,c3,c3,cut2,ghal,dhal
cc    &   ,ev,vxx,vxy,vxz,vyy,nproc
cc    &   ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
c
c        do ii = iwarp, nvdwlocnlb_pair-1, nwarp
c
c           ! Set Data to compute to zero
c           ev_    = 0
c           gxi(threadIdx%x)    = 0
c           gyi(threadIdx%x)    = 0
c           gzi(threadIdx%x)    = 0
c           gxk(threadIdx%x)    = 0
c           gyk(threadIdx%x)    = 0
c           gzk(threadIdx%x)    = 0
c           vxx_   = 0
c           vxy_   = 0
c           vxz_   = 0
c           vyy_   = 0
c           vyz_   = 0
c           vzz_   = 0
c
c           ! Set interaction mask
c           ikstat(threadIdx%x) = 0
c
c           ! Load atom block i parameters
c           iblock = ivblst(ii+1)
c           if (iblock.eq.0) cycle
c           idx    = (iblock-1)*warpsize + ilane 
c           iglob(threadIdx%x)  = cellv_glob(idx)
c           i      = cellv_loc (idx)
c           it(threadIdx%x) = jvdw  (idx)
c           xi(threadIdx%x) = xred  (idx)
c           yi(threadIdx%x) = yred  (idx)
c           zi(threadIdx%x) = zred  (idx)
c#ifndef TINKER_NO_MUTATE
c           muti(threadIdx%x)   = mut(iglob(threadIdx%x))
c#endif
c
c           ! Load atom block k neighbor parameter
c           kdx    = vblst( ii*warpsize + ilane )
c           kglob(threadIdx%x)  = cellv_glob(kdx)
c           kbis   = cellv_loc(kdx)
c           kt(threadIdx%x) = jvdw(kdx)
c           xk(threadIdx%x) = xred(kdx)
c           yk(threadIdx%x) = yred(kdx)
c           zk(threadIdx%x) = zred(kdx)
c           same_block = (idx.ne.kdx)
c#ifndef TINKER_NO_MUTATE
c           mutk(threadIdx%x) = mut(kglob(threadIdx%x))
c#endif
c
c           ! Interact block i with block k
c           do i1 = 0, 1
c              srclane = iand( ilane+i1-1,warpsize-1 ) + 1
c              iilane  = threadIdx%x-ilane + srclane
c           do j = 0,warpsize-1
c              ist     = AtomicOr(ikstat(iilane),ishft(1,j))
c              if (btest(ist,j)) cycle
c              srclane = iand( ilane+j-1,warpsize-1 ) + 1
c              klane   = threadIdx%x-ilane + srclane
c              if (nproc.gt.1) then
c                 xk_   = xk(klane)
c                 yk_   = yk(klane)
c                 zk_   = zk(klane)
c                 xpos  = xi(iilane) - xk_
c                 ypos  = yi(iilane) - yk_
c                 zpos  = zi(iilane) - zk_
c                 call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
c                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
c     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
c     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
c                    accept_mid = .false.
c                 else
c                    accept_mid = .true.
c                 end if
c              else
c                 xpos    = xi(iilane) - xk(klane)
c                 ypos    = yi(iilane) - yk(klane)
c                 zpos    = zi(iilane) - zk(klane)
c                 call image_inl(xpos,ypos,zpos)
c              end if
c
c              dedx = 0.0; dedy = 0.0; dedz = 0.0;
c
c              rik2  = xpos**2 + ypos**2 + zpos**2
c
c              do_pair = merge(.true.,iglob(iilane).lt.kglob(klane)
c     &                       ,same_block)
c              if (do_pair.and.rik2<=off2
c     &           .and.accept_mid) then
c
c                 rv2   =  radmin_t (kt(klane),it(iilane))
c                 eps2  = epsilon_t (kt(klane),it(iilane))
c
c                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
c     &                         ,cut2,cut,off,ghal,dhal
c     &                         ,scexp,vlambda,scalpha,muti(iilane)
c     &                         ,mutk(klane),e,dedx,dedy,dedz)
c
c                 ev_   = ev_   + e
c
c                 vxx_  = vxx_  + xpos * dedx
c                 vxy_  = vxy_  + ypos * dedx
c                 vxz_  = vxz_  + zpos * dedx
c                 vyy_  = vyy_  + ypos * dedy
c                 vyz_  = vyz_  + zpos * dedy
c                 vzz_  = vzz_  + zpos * dedz
c
c
c                 !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1
c
c                 ! Accumulate interaction gradient
c                 gxi(iilane) = gxi(iilane) + real(dedx,r_p)
c                 gyi(iilane) = gyi(iilane) + real(dedy,r_p)
c                 gzi(iilane) = gzi(iilane) + real(dedz,r_p)
c
c                 gxk(klane)  = gxk(klane) + real(dedx,r_p)
c                 gyk(klane)  = gyk(klane) + real(dedy,r_p)
c                 gzk(klane)  = gzk(klane) + real(dedz,r_p)
c
c              end if
c
c           end do
c           end do
c
c           call syncwarp(ALL_LANES)
c           ist  = iand(ithread-1,RED_BUFF_SIZE-1) + 1
c           !increment the van der Waals energy
c           istat = atomicAdd( ev_buff(ist) ,ev_ )
c 
c           !increment the van der Waals derivatives
c           if (idx.le.nvdwlocnl) then
c              istat   = atomicAdd (dev(1,i),gxi(threadIdx%x))
c              istat   = atomicAdd (dev(2,i),gyi(threadIdx%x))
c              istat   = atomicAdd (dev(3,i),gzi(threadIdx%x))
c           end if
c
c           if (kdx.le.nvdwlocnl) then
c              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
c              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
c              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
c           end if
c
c           ! Increment virial term of van der Waals
c           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+ist),vxx_ )
c           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+ist),vxy_ )
c           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+ist),vxz_ )
c           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+ist),vyy_ )
c           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+ist),vyz_ )
c           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+ist),vzz_ )
c
c        end do
c        end subroutine

        attributes(global)
     &  subroutine ehalshortlong3_cu
     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,ev_buff,nev_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,shortcut2,off
     &           ,scexp,vlambda,scalpha,mut
     &           ,vdwshortcut,shortheal,ghal,dhal,use_short
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &           )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5,cut2,cut
     &           ,ghal,dhal,off2,off,shortcut2,vdwshortcut
     &           ,shortheal,scexp,vlambda,scalpha
        logical  ,device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        real(r_p),device :: ev_buff(RED_BUFF_SIZE)
        integer  ,device ::nev_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb),ivblst(nvdwlocnlb_pair)
     &           ,vblst(nvdwlocnlb_pair*(BLOCK_SIZE))
     &           ,cellv_loc(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
        logical  ,value,intent(in)::use_short
#ifdef TINKER_DEBUG
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kglob,kbis,kt,kt_,ist,nev_
        real(t_p) e,istat
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk,yk,zk,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz,devx,devy,devz
        logical do_pair,same_block,accept_mid
        logical muti
        logical,shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Set Data to compute to zero
           ev_    = 0
           nev_   = 0

           ! Load atom block i parameters
           iblock = ivblst(ii+1)
           if (iblock.eq.0) cycle
           idx    = (iblock-1)*warpsize + ilane 
           iglob  = cellv_glob(idx)
           i      = cellv_loc (idx)
           it     = jvdw  (idx)
           xi     = xred  (idx)
           yi     = yred  (idx)
           zi     = zred  (idx)
#ifndef TINKER_NO_MUTATE
           muti   = mut(iglob)
#endif

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt     = jvdw  (kdx)
           xk     = xred  (kdx)
           yk     = yred  (kdx)
           zk     = zred  (kdx)
           same_block = (idx.ne.kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob)
#endif

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
              kglob_  = __shfl(kglob,srclane)
#endif
              kt_     = __shfl(kt   ,srclane)

              if (nproc.gt.1) then
                 xk_   = __shfl(xk,srclane)
                 yk_   = __shfl(yk,srclane)
                 zk_   = __shfl(zk,srclane)
                 xpos  = xi - xk_
                 ypos  = yi - yk_
                 zpos  = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
                 call image_inl(xpos,ypos,zpos)
              end if

              dedx = 0.0; dedy = 0.0; dedz = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                       ,same_block)
              if (do_pair.and.rik2>=shortcut2.and.rik2<=off2
     &           .and.accept_mid) then

                 rv2   =  radmin_t (kt_,it)
                 eps2  = epsilon_t (kt_,it)

                 if (use_short) then
                 call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,off
     &                     ,scexp,vlambda,scalpha,muti,mutk(klane)
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
                 else
                 call ehal1_couple_long(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,cut,off,vdwshortcut
     &                     ,scexp,vlambda,scalpha,muti,mutk(klane)
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
                 end if

                 ev_   = ev_ + e
                 nev_  = nev_+ 1

#ifdef TINKER_DEBUG
#endif
              end if

              !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1
           end do

           kt_ = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(kt_) , ev_ )
           istat = atomicAdd(nev_buff(kt_) ,nev_ )
 
        end do
        end subroutine


        attributes(global)
     &  subroutine ehal1short_cu
     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,dev,ev_buff,vir_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,off2,off
     &           ,scexp,vlambda,scalpha,mut
     &           ,shortheal,ghal,dhal,use_short
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &           ,inter,rank
#endif
     &           )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
     &           ,cut2,ghal,dhal,off2,off,shortheal
     &           ,scexp,vlambda,scalpha
        logical  ,device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        real(r_p),device :: ev_buff (RED_BUFF_SIZE)
        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb),ivblst(nvdwlocnlb_pair)
     &           ,vblst(nvdwlocnlb_pair*(BLOCK_SIZE))
     &           ,cellv_loc(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
        real(r_p),device,intent(inout)::dev(3,nbloc)
        logical  ,value,intent(in)::use_short
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
        integer,value :: rank
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kbis,ist
        integer,shared,dimension(VDW_BLOCK_DIM):: kt,kglob
        real(t_p) e,istat
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz,devx,devy,devz
        real(r_p),shared,dimension(VDW_BLOCK_DIM):: 
     &            gxi,gyi,gzi,gxk,gyk,gzk
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid
        logical muti
        logical,shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.
c       ninte   = 0

c       if (ithread.eq.1) print*,'ehal1c_cu2_in'
c    &   ,blockDim%x,gridDim%x,nwarp,nvdwlocnlb_pair,n
c    &   ,c0,c1,c2,c3,c3,cut2,ghal,dhal
c    &   ,ev,vxx,vxy,vxz,vyy,nproc
c    &   ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Set Data to compute to zero
           ev_    = 0
           gxi(threadIdx%x)  = 0
           gyi(threadIdx%x)  = 0
           gzi(threadIdx%x)  = 0
           gxk(threadIdx%x)  = 0
           gyk(threadIdx%x)  = 0
           gzk(threadIdx%x)  = 0
           vxx_   = 0
           vxy_   = 0
           vxz_   = 0
           vyy_   = 0
           vyz_   = 0
           vzz_   = 0

           ! Load atom block i parameters
           iblock = ivblst(ii+1)
           if (iblock.eq.0) cycle
           idx    = (iblock-1)*warpsize + ilane 
           iglob  = cellv_glob(idx)
           i      = cellv_loc (idx)
           it     = jvdw  (idx)
           xi     = xred  (idx)
           yi     = yred  (idx)
           zi     = zred  (idx)
#ifndef TINKER_NO_MUTATE
           muti   = mut(iglob)
#endif

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob(threadIdx%x)  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt(threadIdx%x) = jvdw  (kdx)
           xk(threadIdx%x) = xred  (kdx)
           yk(threadIdx%x) = yred  (kdx)
           zk(threadIdx%x) = zred  (kdx)
           same_block = (idx.ne.kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob(threadIdx%x))
#endif

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
#endif

              if (nproc.gt.1) then
                 xk_   = xk(klane)
                 yk_   = yk(klane)
                 zk_   = zk(klane)
                 xpos  = xi - xk_
                 ypos  = yi - yk_
                 zpos  = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 xpos    = xi - xk(klane)
                 ypos    = yi - yk(klane)
                 zpos    = zi - zk(klane)
                 call image_inl(xpos,ypos,zpos)
              end if

              dedx = 0.0; dedy = 0.0; dedz = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                       ,same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

                 rv2   =  radmin_t (kt(klane),it)
                 eps2  = epsilon_t (kt(klane),it)

                 call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,off
     &                     ,scexp,vlambda,scalpha,muti,mutk(klane)
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

                 ev_   = ev_ + e

#ifdef TINKER_DEBUG
                 if (iglob<kglob(klane)) then
                    ist = Atomicadd(inter(iglob),1)
                 else
                    ist = Atomicadd(inter(kglob(klane)),1)
                 end if
c                if (rank.eq.0) then
c                if (iglob.eq.12) then
c                   if (iglob<kglob(klane)) then
c                      print*,kglob(klane), kdx_,ii+1,iblock,j,ilane,'i'
c                   else
c                      print*, kglob(klane),real(-xpos,4),real(-ypos,4),-zpos
c                   end if
c                end if
c                if (kglob(klane).eq.12) then
c                   if (iglob<kglob(klane)) then
c                      print*, iglob, kdx_,ii+1,iblock,j,ilane,'i'
c                   else
c                      print*, iglob,real(-xpos,4),real(-ypos,4),-zpos
c                   endif
c                end if
c                end if
#endif

                 vxx_  = vxx_  + xpos * dedx
                 vxy_  = vxy_  + ypos * dedx
                 vxz_  = vxz_  + zpos * dedx
                 vyy_  = vyy_  + ypos * dedy
                 vyz_  = vyz_  + zpos * dedy
                 vzz_  = vzz_  + zpos * dedz

                 !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

                 ! Resolve gradient
                 gxi(threadIdx%x) = gxi(threadIdx%x) + real(dedx,r_p)
                 gyi(threadIdx%x) = gyi(threadIdx%x) + real(dedy,r_p)
                 gzi(threadIdx%x) = gzi(threadIdx%x) + real(dedz,r_p)

                 gxk(klane) = gxk(klane) + real(dedx,r_p)
                 gyk(klane) = gyk(klane) + real(dedy,r_p)
                 gzk(klane) = gzk(klane) + real(dedz,r_p)

              end if

           end do

           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )
 
           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi(threadIdx%x))
              istat   = atomicAdd (dev(2,i),gyi(threadIdx%x))
              istat   = atomicAdd (dev(3,i),gzi(threadIdx%x))
           end if

           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
           end if

           ! Increment virial term of van der Waals
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )

        end do
        end subroutine

        attributes(global)
     &  subroutine ehal1long_cu
     &           (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,dev,ev_buff,vir_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,shortcut2
     &           ,scexp,vlambda,scalpha,mut
     &           ,vdwshortcut,shortheal,ghal,dhal,use_short
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &           ,inter,rank
#endif
     &           )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
     &           ,cut2,cut,ghal,dhal,off2,off,vdwshortcut,shortheal
     &           ,shortcut2
     &           ,scexp,vlambda,scalpha
        logical  ,device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        real(r_p),device :: ev_buff (RED_BUFF_SIZE)
        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb),ivblst(nvdwlocnlb_pair)
     &           ,vblst(nvdwlocnlb_pair*(BLOCK_SIZE))
     &           ,cellv_loc(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
        real(r_p),device,intent(inout)::dev(3,nbloc)
        logical  ,value,intent(in)::use_short
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
        integer,value :: rank
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kglob,kbis,ist
        integer,shared,dimension(VDW_BLOCK_DIM):: kt
        real(t_p) e,istat
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p),shared,dimension(VDW_BLOCK_DIM):: xk,yk,zk
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz
        real(r_p),shared,dimension(VDW_BLOCK_DIM)::
     &            gxi,gyi,gzi,gxk,gyk,gzk
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid
        logical muti
        logical,shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.
c       ninte   = 0

c       if (ithread.eq.1) print*,'ehal1c_cu2_in'
c    &   ,blockDim%x,gridDim%x,nwarp,nvdwlocnlb_pair,n
c    &   ,c0,c1,c2,c3,c3,cut2,ghal,dhal
c    &   ,ev,vxx,vxy,vxz,vyy,nproc
c    &   ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Set Data to compute to zero
           ev_    = 0
           gxi(threadIdx%x)  = 0
           gyi(threadIdx%x)  = 0
           gzi(threadIdx%x)  = 0
           gxk(threadIdx%x)  = 0
           gyk(threadIdx%x)  = 0
           gzk(threadIdx%x)  = 0
           vxx_   = 0
           vxy_   = 0
           vxz_   = 0
           vyy_   = 0
           vyz_   = 0
           vzz_   = 0

           ! Load atom block i parameters
           iblock = ivblst(ii+1)
           if (iblock.eq.0) cycle
           idx    = (iblock-1)*warpsize + ilane 
           iglob  = cellv_glob(idx)
           i      = cellv_loc (idx)
           it     = jvdw  (idx)
           xi     = xred  (idx)
           yi     = yred  (idx)
           zi     = zred  (idx)
#ifndef TINKER_NO_MUTATE
           muti   = mut(iglob)
#endif

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt(threadIdx%x) = jvdw  (kdx)
           xk(threadIdx%x) = xred  (kdx)
           yk(threadIdx%x) = yred  (kdx)
           zk(threadIdx%x) = zred  (kdx)
           same_block = (idx.ne.kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob)
#endif

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
              kglob_  = __shfl(kglob,srclane)
#endif

              if (nproc.gt.1) then
                 xk_   = xk(klane)
                 yk_   = yk(klane)
                 zk_   = zk(klane)
                 xpos  = xi - xk_
                 ypos  = yi - yk_
                 zpos  = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 xpos    = xi - xk(klane)
                 ypos    = yi - yk(klane)
                 zpos    = zi - zk(klane)
                 call image_inl(xpos,ypos,zpos)
              end if

              dedx = 0.0; dedy = 0.0; dedz = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                       ,same_block)
              if (do_pair.and.rik2>=shortcut2.and.rik2<=off2
     &           .and.accept_mid) then

                 rv2   =  radmin_t (kt(klane),it)
                 eps2  = epsilon_t (kt(klane),it)

                 call ehal1_couple_long(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,cut,off,vdwshortcut
     &                     ,scexp,vlambda,scalpha,muti,mutk(klane)
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

                 ev_   = ev_ + e

#ifdef TINKER_DEBUG
                 if (iglob<kglob_) then
                    ist = Atomicadd(inter(iglob),1)
                 else
                    ist = Atomicadd(inter(kglob_),1)
                 end if
c                if (rank.eq.0) then
c                if (iglob.eq.12) then
c                   if (iglob<kglob_) then
c                      print*,kglob_, kdx_,ii+1,iblock,j,ilane,'i'
c                   else
c                      print*, kglob_,real(-xpos,4),real(-ypos,4),-zpos
c                   end if
c                end if
c                if (kglob_.eq.12) then
c                   if (iglob<kglob_) then
c                      print*, iglob, kdx_,ii+1,iblock,j,ilane,'i'
c                   else
c                      print*, iglob,real(-xpos,4),real(-ypos,4),-zpos
c                   endif
c                end if
c                end if
#endif

                 vxx_  = vxx_  + xpos * dedx
                 vxy_  = vxy_  + ypos * dedx
                 vxz_  = vxz_  + zpos * dedx
                 vyy_  = vyy_  + ypos * dedy
                 vyz_  = vyz_  + zpos * dedy
                 vzz_  = vzz_  + zpos * dedz

                 !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

                 ! Accumulate gradient
                 gxi(threadIdx%x)= gxi(threadIdx%x) + real(dedx,r_p)
                 gyi(threadIdx%x)= gyi(threadIdx%x) + real(dedy,r_p)
                 gzi(threadIdx%x)= gzi(threadIdx%x) + real(dedz,r_p)

                 gxk(klane) = gxk(klane) + real(dedx,r_p)
                 gyk(klane) = gyk(klane) + real(dedy,r_p)
                 gzk(klane) = gzk(klane) + real(dedz,r_p)

              end if

           end do

           it = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )
 
           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi(threadIdx%x))
              istat   = atomicAdd (dev(2,i),gyi(threadIdx%x))
              istat   = atomicAdd (dev(3,i),gzi(threadIdx%x))
           end if

           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
           end if

           ! Increment virial term of van der Waals
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )

        end do
        end subroutine
        !TODO Instruct reduction kernels for energy and virial

        subroutine set_vdw_texture(kred,radmin,epsilon,xred,yred,zred
     &                   ,vblst,ivblst,loc_ired,jvdw,cellv_glob
     &                   ,cellv_loc,n,nvdwlocnl,nvdwlocnlb,nvdwclass
     &                   ,nvdwlocnlb_pair)
        implicit none
        integer  ,value ,intent(in):: n,nvdwlocnl,nvdwlocnlb
     &           ,nvdwlocnlb_pair,nvdwclass
        integer  ,target,device:: vblst(nvdwlocnlb_pair*(BLOCK_SIZE+2))
     &           ,ivblst(nvdwlocnlb_pair),jvdw(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb),cellv_glob(nvdwlocnlb)
     &           ,cellv_loc(nvdwlocnlb)
        real(t_p),target,device:: kred(n),xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb)
     &           ,zred(nvdwlocnlb),radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass)

        type(c_devptr)::base_devptr
        integer(cuda_count_kind) ::dpitch,spitch
        integer   istat
        real(t_p) rstat
        logical,save:: first_entry=.true.

c       base_devptr = c_devloc(loc_ired)
c       call c_f_pointer(base_devptr,loc_ired_t,(/nvdwlocnlb/))
c       base_devptr = c_devloc(jvdw)
c       call c_f_pointer(base_devptr,jvdw_t,(/nvdwlocnlb/))
c       base_devptr = c_devloc(xred)
c       call c_f_pointer(base_devptr,xred_p,(/nvdwlocnlb/))
c       base_devptr = c_devloc(yred)
c       call c_f_pointer(base_devptr,yred_p,(/nvdwlocnlb/))
c       base_devptr = c_devloc(zred)
c       call c_f_pointer(base_devptr,zred_p,(/nvdwlocnlb/))
c       base_devptr = c_devloc(vblst)
c       call c_f_pointer(base_devptr,vblst_t,
c    &               (/nvdwlocnlb_pair*(BLOCK_SIZE+2)/))
c       base_devptr = c_devloc(ivblst)
c       call c_f_pointer(base_devptr,ivblst_t,(/nvdwlocnlb_pair/))

        !xred_t   => xred
        !yred_t   => yred
        !zred_t   => zred
        !vblst_t  => vblst
        !ivblst_t => ivblst
        !jvdw_t   => jvdw
        !loc_ired_t => loc_ired

        if (.not.first_entry) return
        kred_t    => kred
        radmin_t  => radmin
        epsilon_t => epsilon
c       base_devptr = c_devloc(kred)
c       call c_f_pointer(base_devptr,kred_t,(/n/))
c       base_devptr = c_devloc(radmin)
c       call c_f_pointer(base_devptr,radmin_t,(/nvdwclass,nvdwclass/))
c       base_devptr = c_devloc(epsilon)
c       call c_f_pointer(base_devptr,epsilon_t,(/nvdwclass,nvdwclass/))

c       allocate (ired_p(n))
c       allocate (kred_p(n))
c       allocate (radmin_p(maxclass,maxclass))
c       allocate (epsilon_p(maxclass,maxclass))
c       istat = cudaMallocPitch(radmin_p,dpitch,maxclass,maxclass)
c       istat = cudaMallocPitch(epsilon_p,spitch,maxclass,maxclass)

c       print*,size(radmin_p),shape(epsilon_p),spitch,dpitch
c       dpitch = maxclass
c       spitch = dpitch
c       istat  = cudaMemcpy(ired_p,ired,n)
c       istat  = istat + cudaMemcpy(kred_p,kred,n)
c       istat  = istat + cudaMemcpy2d(radmin_p,dpitch,radmin,
c    &           maxclass,maxclass,maxclass)
c       istat  = istat + cudaMemcpy2d(epsilon_p,spitch,epsilon,
c    &           maxclass,maxclass,maxclass)

c       if (istat.ne.0) then
c          print*, 'Error allocating Cuda Fortran Array',istat
c       end if
        first_entry=.false.

        end subroutine

      end module

      subroutine ehal1c_kernel
      use ehal1cu
      implicit none
      end subroutine
#else
      ! For Portability issue
      subroutine void_ehalcu
      end
#endif
