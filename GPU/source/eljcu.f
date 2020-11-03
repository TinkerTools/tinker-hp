c
c     sorbonne university
c     washington university in saint louis
c     university of texas at austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elj1cu  --  buffered 14-7 energy & derivatives ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der waals energy and
c     its first derivatives with respect to cartesian coordinates
c
#define TINKER_CUF
#include  "tinker_precision.h"
      module eljcu
        use cudafor
        use tinheader,only: ti_p
        use tintypes ,only: real3
        use sizes    ,only: maxclass
        use utilcu   ,only: ndir,VDW_BLOCK_DIM,ALL_LANES,use_virial
        use utilgpu  ,only: BLOCK_SIZE,RED_BUFF_SIZE

        contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "switch_respa.f.inc"
#include "pair_elj.f.inc"

        attributes(global)
     &  subroutine elj1_cu
     &            (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,dev,ev_buff,vir_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
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

        integer ithread,iwarp,nwarp,ilane,srclane,klane
        integer dstlane,iblock
        integer idx,kdx,kdx_,ii,j,i,iglob,it,ivloc
        integer kbis,kt_,kvloc,ist
        integer,shared,dimension(VDW_BLOCK_DIM)::kglob,kt
        real(t_p) e,istat
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
        real(t_p) rik2,rv2,eps2,redi,redk,redk_
        type(real3) ded
        real(t_p) devx,devy,devz
        real(r_p) gxi,gyi,gzi
        real(r_p),shared,dimension(VDW_BLOCK_DIM):: gxk,gyk,gzk
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.
c       ninte   = 0

c       if (ithread.eq.1) print*,'ehal1c_cu2_in'
c    &   ,blockDim%x,gridDim%x,nwarp,nvdwlocnlb_pair,n
c    &   ,c0,c1,c2,c3,c3,cut2,ghal,dhal
c    &   ,ev,vxx,vxy,vxz,vyy,ndir
c    &   ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Set Data to compute to zero
           ev_    = 0
           gxi    = 0
           gyi    = 0
           gzi    = 0
           gxk(threadIdx%x) = 0
           gyk(threadIdx%x) = 0
           gzk(threadIdx%x) = 0
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

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob(threadIdx%x) = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kt(threadIdx%x) = jvdw  (kdx)
           xk(threadIdx%x) = xred  (kdx)
           yk(threadIdx%x) = yred  (kdx)
           zk(threadIdx%x) = zred  (kdx)
           same_block = (idx.ne.kdx)

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
#endif
              kt_     = kt(klane)

              if (ndir.gt.1) then
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

              ded%x = 0.0; ded%y = 0.0; ded%z = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                       ,same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

                 rv2   =  radmin (kt_,it)
                 eps2  = epsilon (kt_,it)

                 call elj1_couple(rik2,xpos,ypos,zpos,rv2,eps2,cut2
     &                           ,cut,off,e,ded)

                 ev_   = ev_ + e

#ifdef TINKER_DEBUG
                 if (iglob<kglob(klane)) then
                    ist = Atomicadd(inter(iglob),1)
                 else
                    ist = Atomicadd(inter(kglob(klane)),1)
                 end if
#endif

                 if (use_virial) then
                 vxx_  = vxx_  + xpos * ded%x
                 vxy_  = vxy_  + ypos * ded%x
                 vxz_  = vxz_  + zpos * ded%x
                 vyy_  = vyy_  + ypos * ded%y
                 vyz_  = vyz_  + zpos * ded%y
                 vzz_  = vzz_  + zpos * ded%z
                 end if

              end if

              dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

              ! Accumulate interaction gradient
              gxi = gxi + real(ded%x,r_p)
              gyi = gyi + real(ded%y,r_p)
              gzi = gzi + real(ded%z,r_p)

              gxk(klane) = gxk(klane) + real(ded%x,r_p)
              gyk(klane) = gyk(klane) + real(ded%y,r_p)
              gzk(klane) = gzk(klane) + real(ded%z,r_p)

           end do

           call syncwarp(ALL_LANES)
           kt_   = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(kt_) ,ev_ )
 
           !increment the van der Waals derivatives
           if (idx.le.nvdwlocnl) then
              istat   = atomicAdd (dev(1,i),gxi)
              istat   = atomicAdd (dev(2,i),gyi)
              istat   = atomicAdd (dev(3,i),gzi)
           end if

           if (kdx.le.nvdwlocnl) then
              istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
              istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
              istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
           end if

           ! Increment virial term of van der Waals
           if (use_virial) then
           istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+kt_),vxx_ )
           istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+kt_),vxy_ )
           istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+kt_),vxz_ )
           istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+kt_),vyy_ )
           istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+kt_),vyz_ )
           istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+kt_),vzz_ )
           end if

        end do
        end subroutine

        attributes(global)
     &  subroutine elj3_cu
     &            (xred,yred,zred,cellv_glob,cellv_loc,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,ev_buff,nev_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass
     &           ,c0,c1,c2,c3,c4,c5,cut2,cut,off2,off,ghal,dhal
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &           )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
     &           ,cut2,cut,off2,off,ghal,dhal
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        real(r_p),device :: ev_buff (RED_BUFF_SIZE)
        integer  ,device :: nev_buff(RED_BUFF_SIZE)
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
        integer dstlane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it,ivloc
        integer kglob,kbis,kt,kt_,kvloc,ist
        integer nev_, istat
        real(t_p) e,rstat
        real(r_p) ev_
        real(t_p) xi,yi,zi,xk,yk,zk,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p) rik2,rv2,eps2,redi,redk,redk_
        logical do_pair,same_block,accept_mid

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.
c       ninte   = 0

c       if (ithread.eq.1) print*,'ehal1c_cu2_in'
c    &   ,blockDim%x,gridDim%x,nwarp,nvdwlocnlb_pair,n
c    &   ,c0,c1,c2,c3,c3,cut2,ghal,dhal
c    &   ,ev,vxx,vxy,vxz,vyy,ndir
c    &   ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

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
           ivloc  = loc_ired  (idx)
           it     = jvdw  (idx)
           xi     = xred  (idx)
           yi     = yred  (idx)
           zi     = zred  (idx)
           redi   = merge (1.0_ti_p,kred(iglob),(i.eq.ivloc))

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob  = cellv_glob(kdx)
           kbis   = cellv_loc(kdx)
           kvloc  = loc_ired  (kdx)
           kt     = jvdw  (kdx)
           xk     = xred  (kdx)
           yk     = yred  (kdx)
           zk     = zred  (kdx)
           same_block = (idx.ne.kdx)
           redk   = merge (1.0_ti_p,kred(kglob),(kbis.eq.kvloc))

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
              kglob_  = __shfl(kglob,srclane)
#endif
              kt_     = __shfl(kt   ,srclane)

              if (ndir.gt.1) then
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
                 xpos  = xi - __shfl(xk,srclane)
                 ypos  = yi - __shfl(yk,srclane)
                 zpos  = zi - __shfl(zk,srclane)
                 call image_inl(xpos,ypos,zpos)
              end if

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                       ,same_block)

              if (do_pair.and.rik2<=off2.and.accept_mid) then

                 rv2   =  radmin (kt_,it)
                 eps2  = epsilon (kt_,it)

                 call elj3_couple(rik2,xpos,ypos,zpos,rv2,eps2,cut2
     &                           ,cut,off,e)

                 ev_   = ev_ + e
                 nev_  = nev_+ 1

#ifdef TINKER_DEBUG
c                if (iglob<kglob_) then
c                   ist = Atomicadd(inter(iglob),1)
c                else
c                   ist = Atomicadd(inter(kglob_),1)
c                end if
#endif

              end if

              dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

           end do

           kt_   = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           rstat = atomicAdd( ev_buff(kt_) ,ev_ )
           istat = atomicAdd(nev_buff(kt_) ,nev_ )
 
        end do
        end subroutine

      end module
