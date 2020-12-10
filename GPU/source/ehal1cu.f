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
#include  "tinker_types.h"
      module ehal1cu
        use cudafor
        use tinheader,only: ti_p
        use sizes    ,only: maxclass
        use utilcu   ,only: nproc,all_lanes,VDW_BLOCK_DIM
     &               ,skipvdw12,cu_update_skipvdw12,use_virial
     &               ,vcouple
        use utilgpu  ,only: BLOCK_SIZE,RED_BUFF_SIZE
        use sizes    ,only: maxvalue
        use vdw      ,only: vdw_skipvdw12=>skipvdw12

        integer(1),private:: one1,two1
        integer  ,pointer,device::ired_t(:),cellv_glob_t(:)
     &           ,cellv_loc_t(:),loc_ired_t(:),vblst_t(:),ivblst_t(:)
     &           ,jvdw_t(:),i12_p(:,:)
        real(t_p),pointer,texture::radmin_t(:,:),epsilon_t(:,:)
        real(t_p),pointer,device::
     &            kred_t(:),xred_t(:),yred_t(:),zred_t(:)
        parameter(one1=1, two1=2)

        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
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
        integer(1),device :: mut(n)
        ener_rtyp,device :: ev_buff(RED_BUFF_SIZE)
        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::cellv_glob(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb)
     &           ,vblst(nvdwlocnlb_pair*2),cellv_loc(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
        mdyn_rtyp,device,intent(inout)::dev(3,nbloc)
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
#endif

        integer ithread,iwarp,nwarp,ilane,istat,srclane
        integer dstlane,klane
        integer idx,kdx,kdx_,ii,j,i,iglob,it
        integer kglob,kbis,kt,kt_,kvloc,k
        integer kglob_
        integer ai12(maxvalue)
        real(t_p) e
        ener_rtyp ev_
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz
        mdyn_rtyp gxi,gyi,gzi
        mdyn_rtyp gxk,gyk,gzk
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,ik12
        integer(1) muti,mutik
        integer(1),shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

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
           if (skipvdw12) then
             do k = 1,maxvalue
                ai12(k) = i12_p(k,iglob)
             end do
           end if

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
              kglob_  = __shfl(kglob,srclane)

              xpos    = xi - __shfl(xk,srclane)
              ypos    = yi - __shfl(yk,srclane)
              zpos    = zi - __shfl(zk,srclane)
              mutik   = muti + mutk(klane)

              if (skipvdw12) then
                 ik12    = .false.
                 do k = 1,maxvalue
                    if (ai12(k).eq.kglob_) ik12=.true.
                 end do
              end if

              call image_inl(xpos,ypos,zpos)

              dedx = 0.0; dedy = 0.0; dedz = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              if (ik12.and.skipvdw12) goto 10
              do_pair = merge(.true.,iglob.lt.kglob_,same_block)

              if (do_pair.and.kdx_<=nvdwlocnl.and.rik2<=off2) then

                 ! Annihilate
                 if (vcouple.and.mutik.eq.two1) mutik=one1
                 rv2   =  radmin (kt_,it)
                 eps2  = epsilon (kt_,it)

                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                         ,cut2,cut,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,mutik
     &                         ,e,dedx,dedy,dedz)

                 ev_   = ev_ + tp2enr(e)

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

 10           continue
              dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

              ! Increment gradient
              gxi = gxi + tp2mdr(dedx)
              gyi = gyi + tp2mdr(dedy)
              gzi = gzi + tp2mdr(dedz)

              gxk = gxk + tp2mdr(__shfl(dedx,dstlane))
              gyk = gyk + tp2mdr(__shfl(dedy,dstlane))
              gzk = gzk + tp2mdr(__shfl(dedz,dstlane))

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
        integer(1),device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        ener_rtyp,device::  ev_buff(RED_BUFF_SIZE)
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
        integer dstlane,klane,iblock,k
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kbis,kt_,ist,nev_
        integer ai12(maxvalue)
        real(t_p) e,istat
        ener_rtyp ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz,devx,devy,devz
        logical do_pair,same_block,accept_mid,ik12
        integer(1) muti,mutik
#ifndef TINKER_NO_MUTATE
        integer(1),shared:: mutk(VDW_BLOCK_DIM)
#endif
        integer,shared,dimension(VDW_BLOCK_DIM):: kglob,kt
        real(t_p) ,shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk

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
           do j = 1, maxvalue
              ai12(j) = i12_p(j,iglob)
           end do

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
           call syncwarp(ALL_LANES)
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
#endif
              kglob_  = kglob(klane)
              kt_     = kt(klane)
#ifndef TINKER_NO_MUTATE
              mutik   = muti + mutk(klane)
#endif

              if (skipvdw12) then
                 ! Look for 1-2 bond interaction
                 ! and skip it
                 ik12 = .false.
                 do k = 1,maxvalue
                    if (ai12(k).eq.kglob_) ik12=.true.
                 end do
                 if (ik12) goto 10
              end if

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

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.kglob_
     &                       ,same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

#ifndef TINKER_NO_MUTATE
                 ! Annihilate
                 if (vcouple.and.mutik.eq.two1) mutik=one1
#endif
                 rv2   =  radmin_t (kt_,it)
                 eps2  = epsilon_t (kt_,it)

                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                         ,cut2,cut,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,mutik
     &                         ,e,dedx,dedy,dedz)

                 ev_   = ev_ + tp2enr(e)
                 nev_  = nev_ + 1

#ifdef TINKER_DEBUG
#endif
              end if

 10           continue

              !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

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
        integer(1),device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
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
        mdyn_rtyp,device,intent(inout)::dev(3,nbloc)
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
        integer,value :: rank
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kbis,ist,k
        integer ,shared,dimension(VDW_BLOCK_DIM):: kt,kglob
        integer ai12(maxvalue)
        real(t_p) e,istat
        ener_rtyp ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
        real(t_p) dedx,dedy,dedz
        mdyn_rtyp,shared,dimension(VDW_BLOCK_DIM)::
     &            gxk,gyk,gzk
        mdyn_rtyp gxi,gyi,gzi
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid,ik12
        integer(1) muti,mutik
#ifndef TINKER_NO_MUTATE
        integer(1),shared::mutk(VDW_BLOCK_DIM)
#endif

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.
c       ninte   = 0

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob(threadIdx%x) = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt(threadIdx%x) = jvdw(kdx)
           xk(threadIdx%x) = xred(kdx)
           yk(threadIdx%x) = yred(kdx)
           zk(threadIdx%x) = zred(kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob(threadIdx%x))
#endif
           call syncwarp(ALL_LANES)

           ! Set Data to compute to zero
           ev_    = 0
           gxi                 = 0
           gyi                 = 0
           gzi                 = 0
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
           if (skipvdw12) then
             do k = 1,maxvalue
                ai12(k) = i12_p(k,iglob)
             end do
           end if

           same_block = (idx.ne.kdx)

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
#endif
              kglob_  = kglob(klane)
#ifndef TINKER_NO_MUTATE
              mutik   = muti + mutk(klane)
#endif

              if (skipvdw12) then
                 ! Look for 1-2 bond interaction
                 ! and skip it
                 ik12 =.false.
                 do k = 1, maxvalue
                    if (ai12(k).eq.kglob_) ik12=.true.
                 end do
                 if (ik12) cycle
              end if
                 
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

              do_pair = merge(.true.,iglob.lt.kglob(klane),same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

#ifndef TINKER_NO_MUTATE
                 ! Annihilation
                 if (vcouple.and.mutik.eq.two1) mutik=one1
#endif
                 rv2   =  radmin_t (kt(klane),it)
                 eps2  = epsilon_t (kt(klane),it)

                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                         ,cut2,cut,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,mutik
     &                         ,e,dedx,dedy,dedz)

                 ev_   = ev_ + tp2enr(e)

#ifdef TINKER_DEBUG
                 if (iglob<kglob_) then
                    ist = Atomicadd(inter(iglob),1)
                 else
                    ist = Atomicadd(inter(kglob_),1)
                 end if
#endif

                 vxx_  = vxx_  + xpos * dedx
                 vxy_  = vxy_  + ypos * dedx
                 vxz_  = vxz_  + zpos * dedx
                 vyy_  = vyy_  + ypos * dedy
                 vyz_  = vyz_  + zpos * dedy
                 vzz_  = vzz_  + zpos * dedz

                 !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

                 ! Accumulate interaction gradient
                 gxi = gxi + tp2mdr(dedx)
                 gyi = gyi + tp2mdr(dedy)
                 gzi = gzi + tp2mdr(dedz)

                 gxk(klane) = gxk(klane) + tp2mdr(dedx)
                 gyk(klane) = gyk(klane) + tp2mdr(dedy)
                 gzk(klane) = gzk(klane) + tp2mdr(dedz)

              end if

           end do

           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )

           ! Increment virial term of van der Waals
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )

           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi)
              istat   = atomicAdd (dev(2,i),gyi)
              istat   = atomicAdd (dev(3,i),gzi)
           end if

           call syncwarp(ALL_LANES)
           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
           end if

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
c        integer(1),device :: mut(n)
c        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
c     &           ,p_zbeg,p_zend
c        ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
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
c        mdyn_rtyp,device,intent(inout)::dev(3,nbloc)
c#ifdef TINKER_DEBUG
c        integer,device:: inter(*)
c        integer,value :: rank
c#endif
c
c        integer ithread,iwarp,nwarp,ilane,srclane
c        integer dstlane,klane,iilane,iblock
c        integer idx,kdx,kdx_,kglob_,ii,j,i,i1
c        integer ,shared,dimension(VDW_BLOCK_DIM):: iglob,it
c        integer kbis,ist,k
c        integer ,shared,dimension(VDW_BLOCK_DIM):: kt,kglob,ikstat
c        integer,shared:: ai12(maxvalue,VDW_BLOCK_DIM)
c        real(t_p) e,istat
c        ener_rtyp ev_
c        real(t_p) xk_,yk_,zk_,xpos,ypos,zpos
c        real(t_p) rik2,rv2,eps2
c        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
c     &           ,xi,yi,zi
c        real(t_p) dedx,dedy,dedz
c        mdyn_rtyp,shared,dimension(VDW_BLOCK_DIM)::
c     &            gxi,gyi,gzi,gxk,gyk,gzk
c        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
c        logical do_pair,same_block,accept_mid,ik12
c        integer(1),shared::mutk(VDW_BLOCK_DIM)
c     &            ,muti(VDW_BLOCK_DIM)
c        integer(1) mutik
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
cc    &   ,vxx,vxy,vxz,vyy,nproc
cc    &   ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
c
c        do ii = iwarp, nvdwlocnlb_pair-1, nwarp
c
c           ! Load atom block k neighbor parameter
c           kdx    = vblst( ii*warpsize + ilane )
c           kglob(threadIdx%x)  = cellv_glob(kdx)
c           kbis   = cellv_loc(kdx)
c           kt(threadIdx%x) = jvdw(kdx)
c           xk(threadIdx%x) = xred(kdx)
c           yk(threadIdx%x) = yred(kdx)
c           zk(threadIdx%x) = zred(kdx)
c#ifndef TINKER_NO_MUTATE
c           mutk(threadIdx%x) = mut(kglob(threadIdx%x))
c#endif
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
c           if (skipvdw12) then
c             do k = 1,maxvalue
c                ai12(kn,threadIdx%x) = i12_p(k,iglob)
c             end do
c           end if
c
c           same_block = (idx.ne.kdx)
c           call syncwarp(ALL_LANES)
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
c              kglob_  = kglob(klane)
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
c              if (skipvdw12) then
c                 ! Look for 1-2 bond interaction
c                 ! and skip it
c                 ik12 =.false.
c                 do k = 1, maxvalue
c                    if (ai12(k,iilane).eq.kglob_) ik12=.true.
c                 end do
c                 if (ik12) cycle
c              end if
c
c              mutik = mut(ilane) + mutk(klane)
c              dedx = 0.0; dedy = 0.0; dedz = 0.0;
c
c              rik2  = xpos**2 + ypos**2 + zpos**2
c
c              do_pair = merge(.true.,iglob(iilane).lt.kglob_,same_block)
c              if (do_pair.and.rik2<=off2
c     &           .and.accept_mid) then
c
c#ifndef TINKER_NO_MUTATE
c                 ! Annihilation
c                 if (vcouple.and.mutik.eq.two1) mutik=one1
c#endif
c                 rv2   =  radmin_t (kt(klane),it(iilane))
c                 eps2  = epsilon_t (kt(klane),it(iilane))
c
c                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
c     &                         ,cut2,cut,off,ghal,dhal
c     &                         ,scexp,vlambda,scalpha,mutik
c     &                         ,e,dedx,dedy,dedz)
c
c                 ev_   = ev_   + tp2enr(e)
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
c                 gxi(iilane) = gxi(iilane) + tp2mdr(dedx)
c                 gyi(iilane) = gyi(iilane) + tp2mdr(dedy)
c                 gzi(iilane) = gzi(iilane) + tp2mdr(dedz)
c
c                 gxk(klane)  = gxk(klane) + tp2mdr(dedx)
c                 gyk(klane)  = gyk(klane) + tp2mdr(dedy)
c                 gzk(klane)  = gzk(klane) + tp2mdr(dedz)
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
        integer(1),device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        ener_rtyp,device :: ev_buff(RED_BUFF_SIZE)
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
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it,k
        integer kbis,kt,kt_,ist,nev_
        integer,shared,dimension(VDW_BLOCK_DIM)::kglob
        integer ai12(maxvalue)
        real(t_p) e,istat
        ener_rtyp ev_
        real(t_p) xi,yi,zi,xk,yk,zk,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz,devx,devy,devz
        logical do_pair,same_block,accept_mid,ik12
        integer(1) muti,mutik
        integer(1),shared::mutk(VDW_BLOCK_DIM)

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
           if (skipvdw12) then
             do k = 1,maxvalue
                ai12(k) = i12_p(k,iglob)
             end do
           end if

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob(threadIdx%x)  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt     = jvdw  (kdx)
           xk     = xred  (kdx)
           yk     = yred  (kdx)
           zk     = zred  (kdx)
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
              kglob_  = kglob(klane)
              kt_     = __shfl(kt   ,srclane)
#ifndef TINKER_NO_MUTATE
              mutik   = muti + mutk(klane)
#endif

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

              do_pair = merge(.true.,iglob.lt.kglob_
     &                       ,same_block)

              if (skipvdw12) then
                 ! Look for 1-2 bond interaction
                 ! and skip it
                 ik12 =.false.
                 do k = 1, maxvalue
                    if (ai12(k).eq.kglob_) ik12=.true.
                 end do
                 if (ik12) cycle
              end if
                 
              if (do_pair.and.rik2>=shortcut2.and.rik2<=off2
     &           .and.accept_mid) then

#ifndef TINKER_NO_MUTATE
                 ! Annihilate
                 if (vcouple.and.mutik.eq.two1) mutik=one1
#endif
                 rv2   =  radmin_t (kt_,it)
                 eps2  = epsilon_t (kt_,it)

                 if (use_short) then
                 call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,off
     &                     ,scexp,vlambda,scalpha,mutik
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
                 else
                 call ehal1_couple_long(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,cut,off,vdwshortcut
     &                     ,scexp,vlambda,scalpha,mutik
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)
                 end if

                 ev_   = ev_ + tp2enr(e)
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
        integer(1),device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
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
        mdyn_rtyp,device,intent(inout)::dev(3,nbloc)
        logical  ,value,intent(in)::use_short
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
        integer,value :: rank
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kbis,ist,k
        integer,shared,dimension(VDW_BLOCK_DIM):: kt,kglob
        integer ai12(maxvalue)
        real(t_p) e,istat
        ener_rtyp ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz,devx,devy,devz
        mdyn_rtyp,shared,dimension(VDW_BLOCK_DIM)::gxk,gyk,gzk
        mdyn_rtyp gxi,gyi,gzi
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid,ik12
        integer(1) muti,mutik
        integer(1),shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob(threadIdx%x)  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt(threadIdx%x) = jvdw  (kdx)
           xk(threadIdx%x) = xred  (kdx)
           yk(threadIdx%x) = yred  (kdx)
           zk(threadIdx%x) = zred  (kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob(threadIdx%x))
#endif
           call syncwarp(ALL_LANES)

           ! Set Data to compute to zero
           ev_    = 0
           gxi               = 0
           gyi               = 0
           gzi               = 0
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
           if (skipvdw12) then
             do k = 1,maxvalue
                ai12(k) = i12_p(k,iglob)
             end do
           end if

           same_block = (idx.ne.kdx)

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
              kglob_  = kglob(klane)
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
#endif
#ifndef TINKER_NO_MUTATE
              mutik   = muti + mutk(klane)
#endif

              if (skipvdw12) then
                 ! Look for 1-2 bond interaction
                 ! and skip it
                 ik12 =.false.
                 do k = 1, maxvalue
                    if (ai12(k).eq.kglob_) ik12=.true.
                 end do
                 if (ik12) cycle
              end if

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

              do_pair = merge(.true.,iglob.lt.kglob_,same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

#ifndef TINKER_NO_MUTATE
                 ! Annihilate
                 if (vcouple.and.mutik.eq.two1) mutik=one1
#endif

                 rv2   =  radmin_t (kt(klane),it)
                 eps2  = epsilon_t (kt(klane),it)

                 call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,off
     &                     ,scexp,vlambda,scalpha,mutik
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

                 ev_   = ev_ + tp2enr(e)

#ifdef TINKER_DEBUG
                 if (iglob<kglob(klane)) then
                    ist = Atomicadd(inter(iglob),1)
                 else
                    ist = Atomicadd(inter(kglob(klane)),1)
                 end if
#endif

                 if (use_virial) then
                 vxx_  = vxx_  + xpos * dedx
                 vxy_  = vxy_  + ypos * dedx
                 vxz_  = vxz_  + zpos * dedx
                 vyy_  = vyy_  + ypos * dedy
                 vyz_  = vyz_  + zpos * dedy
                 vzz_  = vzz_  + zpos * dedz
                 end if

                 !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

                 ! Resolve gradient
                 gxi = gxi + tp2mdr(dedx)
                 gyi = gyi + tp2mdr(dedy)
                 gzi = gzi + tp2mdr(dedz)

                 gxk(klane) = gxk(klane) + tp2mdr(dedx)
                 gyk(klane) = gyk(klane) + tp2mdr(dedy)
                 gzk(klane) = gzk(klane) + tp2mdr(dedz)

              end if

           end do

           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )
 
           if (use_virial) then
           ! Increment virial term of van der Waals
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )
           end if

           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi)
              istat   = atomicAdd (dev(2,i),gyi)
              istat   = atomicAdd (dev(3,i),gzi)
           end if

           call syncwarp(ALL_LANES)
           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
           end if

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
        integer(1),device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
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
        mdyn_rtyp,device,intent(inout)::dev(3,nbloc)
        logical  ,value,intent(in)::use_short
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
        integer,value :: rank
#endif

        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kbis,ist,k
        integer,shared,dimension(VDW_BLOCK_DIM):: kt,kglob
        integer ai12(maxvalue)
        real(t_p) e,istat
        ener_rtyp ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p),shared,dimension(VDW_BLOCK_DIM):: xk,yk,zk
        real(t_p) rik2,rv2,eps2
        real(t_p) dedx,dedy,dedz
        mdyn_rtyp,shared,dimension(VDW_BLOCK_DIM)::gxk,gyk,gzk
        mdyn_rtyp gxi,gyi,gzi
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid,ik12
        integer(1) muti,mutik
        integer(1),shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob(threadIdx%x)  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt(threadIdx%x) = jvdw  (kdx)
           xk(threadIdx%x) = xred  (kdx)
           yk(threadIdx%x) = yred  (kdx)
           zk(threadIdx%x) = zred  (kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob(threadIdx%x))
#endif
           call syncwarp(ALL_LANES)

           ! Set Data to compute to zero
           ev_    = 0
           gxi               = 0
           gyi               = 0
           gzi               = 0
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
           if (skipvdw12) then
             do k = 1,maxvalue
                ai12(k) = i12_p(k,iglob)
             end do
           end if

           same_block = (idx.ne.kdx)

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
              kglob_  = kglob(klane)
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
#endif
#ifndef TINKER_NO_MUTATE
              mutik   = muti + mutk(klane)
#endif

              if (skipvdw12) then
                 ! Look for 1-2 bond interaction
                 ! and skip it
                 ik12 =.false.
                 do k = 1, maxvalue
                    if (ai12(k).eq.kglob_) ik12=.true.
                 end do
                 if (ik12) cycle
              end if
                 
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

              do_pair = merge(.true.,iglob.lt.kglob_,same_block)
              if (do_pair.and.rik2>=shortcut2.and.rik2<=off2
     &           .and.accept_mid) then

#ifndef TINKER_NO_MUTATE
                 ! Annihilate
                 if (vcouple.and.mutik.eq.two1) mutik=one1
#endif

                 rv2   =  radmin_t (kt(klane),it)
                 eps2  = epsilon_t (kt(klane),it)

                 call ehal1_couple_long(xpos,ypos,zpos,rik2,rv2,eps2
     &                     ,1.0_ti_p,cut2,cut,off,vdwshortcut
     &                     ,scexp,vlambda,scalpha,mutik
     &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

                 ev_   = ev_ + tp2enr(e)

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
                 gxi = gxi  + tp2mdr(dedx)
                 gyi = gyi  + tp2mdr(dedy)
                 gzi = gzi  + tp2mdr(dedz)

                 gxk(klane) = gxk(klane) + tp2mdr(dedx)
                 gyk(klane) = gyk(klane) + tp2mdr(dedy)
                 gzk(klane) = gzk(klane) + tp2mdr(dedz)

              end if

           end do

           it = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )
 
           ! Increment virial term of van der Waals
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )

           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi)
              istat   = atomicAdd (dev(2,i),gyi)
              istat   = atomicAdd (dev(3,i),gzi)
           end if

           call syncwarp(ALL_LANES)
           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
           end if

        end do
        end subroutine

        subroutine attach_vdwcu_data(i12,kred,radmin,epsilon
     &                   ,n,nvdwclass)
        implicit none
        integer  ,intent(in):: n,nvdwclass
        integer  ,intent(in),device,target:: i12(maxvalue,n)
        real(t_p),intent(in),device,target:: kred(n)
        real(t_p),device,target::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass)
        kred_t    => kred
        i12_p     => i12
        radmin_t  => radmin
        epsilon_t => epsilon
        end subroutine

c        subroutine set_vdw_texture(kred,radmin,epsilon,xred,yred,zred
c     &                   ,vblst,ivblst,loc_ired,jvdw,cellv_glob
c     &                   ,cellv_loc,n,nvdwlocnl,nvdwlocnlb,nvdwclass
c     &                   ,nvdwlocnlb_pair)
c        implicit none
c        integer  ,value ,intent(in):: n,nvdwlocnl,nvdwlocnlb
c     &           ,nvdwlocnlb_pair,nvdwclass
c        integer  ,target,device:: vblst(nvdwlocnlb_pair*(BLOCK_SIZE+2))
c     &           ,ivblst(nvdwlocnlb_pair),jvdw(nvdwlocnlb)
c     &           ,loc_ired(nvdwlocnlb),cellv_glob(nvdwlocnlb)
c     &           ,cellv_loc(nvdwlocnlb)
c        real(t_p),target,device:: kred(n),xred(nvdwlocnlb)
c     &           ,yred(nvdwlocnlb)
c     &           ,zred(nvdwlocnlb),radmin(nvdwclass,nvdwclass)
c     &           ,epsilon(nvdwclass,nvdwclass)
c
c        type(c_devptr)::base_devptr
c        integer(cuda_count_kind) ::dpitch,spitch
c        integer   istat
c        real(t_p) rstat
c        logical,save:: first_entry=.true.
c        ! TODO Clean interface
c
cc       base_devptr = c_devloc(loc_ired)
cc       call c_f_pointer(base_devptr,loc_ired_t,(/nvdwlocnlb/))
cc       base_devptr = c_devloc(jvdw)
cc       call c_f_pointer(base_devptr,jvdw_t,(/nvdwlocnlb/))
cc       base_devptr = c_devloc(xred)
cc       call c_f_pointer(base_devptr,xred_p,(/nvdwlocnlb/))
cc       base_devptr = c_devloc(yred)
cc       call c_f_pointer(base_devptr,yred_p,(/nvdwlocnlb/))
cc       base_devptr = c_devloc(zred)
cc       call c_f_pointer(base_devptr,zred_p,(/nvdwlocnlb/))
cc       base_devptr = c_devloc(vblst)
cc       call c_f_pointer(base_devptr,vblst_t,
cc    &               (/nvdwlocnlb_pair*(BLOCK_SIZE+2)/))
cc       base_devptr = c_devloc(ivblst)
cc       call c_f_pointer(base_devptr,ivblst_t,(/nvdwlocnlb_pair/))
c
c        !xred_t   => xred
c        !yred_t   => yred
c        !zred_t   => zred
c        !vblst_t  => vblst
c        !ivblst_t => ivblst
c        !jvdw_t   => jvdw
c        !loc_ired_t => loc_ired
c
c        if (.not.first_entry) return
c        kred_t    => kred
c        radmin_t  => radmin
c        epsilon_t => epsilon
cc       base_devptr = c_devloc(kred)
cc       call c_f_pointer(base_devptr,kred_t,(/n/))
cc       base_devptr = c_devloc(radmin)
cc       call c_f_pointer(base_devptr,radmin_t,(/nvdwclass,nvdwclass/))
cc       base_devptr = c_devloc(epsilon)
cc       call c_f_pointer(base_devptr,epsilon_t,(/nvdwclass,nvdwclass/))
c
cc       allocate (ired_p(n))
cc       allocate (kred_p(n))
cc       allocate (radmin_p(maxclass,maxclass))
cc       allocate (epsilon_p(maxclass,maxclass))
cc       istat = cudaMallocPitch(radmin_p,dpitch,maxclass,maxclass)
cc       istat = cudaMallocPitch(epsilon_p,spitch,maxclass,maxclass)
c
cc       print*,size(radmin_p),shape(epsilon_p),spitch,dpitch
cc       dpitch = maxclass
cc       spitch = dpitch
cc       istat  = cudaMemcpy(ired_p,ired,n)
cc       istat  = istat + cudaMemcpy(kred_p,kred,n)
cc       istat  = istat + cudaMemcpy2d(radmin_p,dpitch,radmin,
cc    &           maxclass,maxclass,maxclass)
cc       istat  = istat + cudaMemcpy2d(epsilon_p,spitch,epsilon,
cc    &           maxclass,maxclass,maxclass)
c
cc       if (istat.ne.0) then
cc          print*, 'Error allocating Cuda Fortran Array',istat
cc       end if
c        first_entry=.false.
c
c        end subroutine

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
