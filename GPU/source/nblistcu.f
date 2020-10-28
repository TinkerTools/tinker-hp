#ifdef _CUDA
#include "tinker_precision.h"
#define TINKER_CUF
      module nblistcu

c     use utilcu    ,only: all_lanes
      use cudadevice
      use utilgpu   ,only: BLOCK_SIZE
c     use tinheader ,only: ti_p,ti_eps
      integer,parameter::grp=2

      contains
#include "image.f.inc"
#include "midpointimage.f.inc"

        attributes(global)
     &  subroutine filter_lst_mat(cell_glob,block_scan,x,y,z,
     &             matb_lst,nlocnlb,nblock,matsize,nblocknlb_pair,
     &             cutbuff2,lst,nlst)
        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb,
     &                     nblocknlb_pair
        real(t_p),value,intent(in)  :: cutbuff2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: lst(nblocknlb_pair*BLOCK_SIZE)
        integer,device,intent(inout):: nlst(nblock)
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp_b,iwarp_b,ilane
        integer bit_si,bit_sh,bit2,pos2,j,srclane
        integer idx,kdx,iscan,kglob,kneig,it1,it2
        integer flag,jth,nbits,oldnb,oldnb_,accept,bset
        integer jflag
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos
        integer ilane_set
        parameter( bit_si=32, bit_sh=5 )

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp_b = (threadIdx%x-1) / warpsize
        nwarp_b = blockDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do it1 = blockIdx%x-1,nblock-1,gridDim%x
           idx   = it1*warpsize + ilane
           iscan = block_scan(it1+1)*warpsize
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)

           do it2 = it1+iwarp_b, nblock-1, nwarp_b
              pos2  = ishft(it2,-bit_sh)+1
              bit2  = iand (it2,bit_si-1)
              if (.not.btest(matb_lst(it1*matsize+pos2),bit2)) cycle
              kdx   = it2*warpsize + ilane
              kglob = cell_glob(kdx)
              xk    = x(kdx)
              yk    = y(kdx)
              zk    = z(kdx)
              jflag = 0

              if (it1.eq.it2) then
                 lst (iscan+ilane) = kglob
              else

                 do j = 1, warpsize
                    srclane = j
                    xpos    = xi - __shfl(xk,srclane)
                    ypos    = yi - __shfl(yk,srclane)
                    zpos    = zi - __shfl(zk,srclane)
                    call image_inl(xpos,ypos,zpos)
                    if ( xpos**2+ypos**2+zpos**2.lt.cutbuff2 ) then
                       ilane_set = 1
                    else
                       ilane_set = 0
                    end if
                    accept = ballot(ilane_set)
                    if (accept.ne.0) jflag = ibset (jflag,j-1)
                 end do
                 nbits = __popc(jflag)
                 
                 jth  = 0
                 bset = 0
                 do while (bset.ne.ilane.and.jth.lt.warpsize)
                    if (btest(jflag,jth)) bset = bset + 1
                    jth = jth + 1
                 end do
                 kneig = __shfl(kglob,jth)
                 
                 if (ilane.eq.1) oldnb_= atomicAdd(nlst(it1+1),nbits)
                 oldnb = __shfl(oldnb_,1)
                 if (ilane.lt.nbits+1) lst(iscan+oldnb+ilane) = kneig
              end if
           end do
        end do
        end subroutine


        attributes(global)
     &  subroutine filter_lst_sparse
     &             (cell_glob,block_scan,x,y,z,
     &             matb_lst,nlocnlb,nblock,matsize,nblocknlb_pair,
     &             cutbuff2,lst,nlst)

        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nblocknlb_pair
        real(t_p),value,intent(in)  :: cutbuff2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nlst(nblock)
     &                ,lst(nblocknlb_pair*(BLOCK_SIZE+2))
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp,iwarp,ilane
        integer iblock,j,srclane
        integer idx,kdx,kdx_,iscan,kglob,kneig,ii
        integer flag,jth,nbits,oldnb,oldnb_,accept,bset
        integer jflag,ilane_set
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do ii = iwarp,nblocknlb_pair-1,nwarp
           iblock= lst(2*ii+1)
           idx   = (iblock-1)*warpsize + ilane
           iscan = block_scan(iblock)*warpsize + 2*nblocknlb_pair
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)

           kdx   = (lst(2*ii+2)-1)*warpsize + ilane
c          kglob = cell_glob(kdx)
           xk    = x(kdx)
           yk    = y(kdx)
           zk    = z(kdx)

           jth   = warpsize
           nbits = 0

           if (idx.eq.kdx) then
              lst(iscan+ilane) = kdx
           else
              do j = 1, warpsize
                 srclane = j
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
                 call image_inl(xpos,ypos,zpos)
                 ilane_set = merge(1,0,
     &                       (xpos**2+ypos**2+zpos**2.lt.cutbuff2))
                 accept  = ballot(ilane_set)
                 if (accept.ne.0) then
                    ! count the number of validation
                    nbits = nbits + 1
                    ! find corresponding position for ilane
                    if (nbits.eq.ilane) jth = j
                 end if
              end do

              if (ilane.eq.1) oldnb= atomicAdd(nlst(iblock),nbits)

              kdx_  = __shfl(kdx,jth)
              oldnb = __shfl(oldnb,1)

              if (ilane.lt.nbits+1) lst(iscan+oldnb+ilane) = kdx_

           end if
        end do
        end subroutine

        attributes(global)
     &  subroutine filter_lst_sparse1
     &             (cell_glob,block_scan,x,y,z,
     &             matb_lst,nlocnlb,nblock,matsize,nblocknlb_pair,
     &             cutbuff2,lst,nlst_l,nlst_r)

        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nblocknlb_pair
        real(t_p),value,intent(in)  :: cutbuff2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nlst_l(nblock),nlst_r(nblock)
     &                ,lst(nblocknlb_pair*(BLOCK_SIZE+2))
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp,iwarp,ilane
        integer iblock,j,srclane
        integer idx,kdx,kdx_l,kdx_r,iscan,iscan1
        integer nlbits,nrbits,kglob,kneig,ii
        integer flag,jlth,jrth,oldnbl,oldnbr,accept,bset
        integer jflag,ilane_set
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do ii = iwarp,nblocknlb_pair-1,nwarp
           iblock= lst(2*ii+1)
           idx   = (iblock-1)*warpsize + ilane
           iscan = block_scan(iblock)*warpsize + 2*nblocknlb_pair
           if (iblock<nblock) then
              iscan1= block_scan(iblock+1)*warpsize + 2*nblocknlb_pair
           end if
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)

           kdx   = (lst(2*ii+2)-1)*warpsize + ilane
c          kglob = cell_glob(kdx)
           xk    = x(kdx)
           yk    = y(kdx)
           zk    = z(kdx)

           jlth   = warpsize
           jrth   = 1
           nlbits = 0
           nrbits = 0

           if (idx.eq.kdx) then
              lst(iscan+ilane) = kdx
           else
              do j = 1, warpsize
                 srclane = j
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
                 call image_inl(xpos,ypos,zpos)
                 ilane_set = merge(1,0,
     &                       (xpos**2+ypos**2+zpos**2.lt.cutbuff2))
                 accept  = ballot(ilane_set)
                 if (accept.ne.0) then
                    if (__popc(accept)>grp) then
                       ! count the number of validation
                       nlbits = nlbits + 1
                       ! find corresponding position for ilane
                       if (nlbits.eq.ilane) jlth = j
                    else
                       if (warpsize-nrbits.eq.ilane) jrth = j
                       nrbits = nrbits + 1
                    end if
                 end if
              end do

              ! Post processing
              if ((nlbits+nrbits)>0) then
                 if (ilane.eq.1)
     &              oldnbl= atomicAdd(nlst_l(iblock),nlbits)
                 if (ilane.eq.warpsize)
     &              oldnbr = atomicAdd(nlst_r(iblock),nrbits)

                 kdx_l  = __shfl(kdx,jlth)
                 kdx_r  = __shfl(kdx,jrth)
                 oldnbl = __shfl(oldnbl,1)
                 oldnbr = __shfl(oldnbr,warpsize)

                 if      (ilane.lt.nlbits+1) then
                    lst(iscan+oldnbl+ilane) = kdx_l
                    !if (iblock.eq.1) print*,kdx_l,'l'
                 else if (ilane.gt.warpsize-nrbits) then
                    lst(iscan1-oldnbr+(ilane-warpsize)) = kdx_r
                 end if
              end if
           end if
        end do
        end subroutine

c
c       "filter_shortlst_sparse" 
c       Filter both nblist and shortnblist with the same adjacency matrix
c
        attributes(global)
     &  subroutine filter_shortlst_sparse1
     &             (cell_glob,block_scan,x,y,z,
     &              matb_lst,nlocnlb,nblock,matsize,nblocknlb_pair,
     &              shortcutbuff2,cutbufbeg2,cutbuff2,
     &              lst,nlst_l,nlst_r,lst1,nlst_l1,nlst_r1)

        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nblocknlb_pair
        real(t_p),value,intent(in)  :: cutbuff2,shortcutbuff2,cutbufbeg2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nlst_l(nblock),nlst_r(nblock)
     &                ,lst(nblocknlb_pair*(BLOCK_SIZE+2))
        integer,device,intent(inout):: nlst_l1(nblock),nlst_r1(nblock)
     &                ,lst1(nblocknlb_pair*(BLOCK_SIZE+2))
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp,iwarp,ilane
        integer iblock,j,srclane
        integer idx,kdx,kdx_l,kdx_r,iscan,iscan1
        integer nlbits,nrbits,nlbits1,nrbits1,kglob,kneig,ii
        integer flag,jlth,jrth,jlth1,jrth1,oldnbl,oldnbr
     &         ,oldnbl1,oldnbr1,accept,accept1
        integer jflag,ilane_set,ilane_set1
        real(t_p) r2,xi,yi,zi,xk,yk,zk,xpos,ypos,zpos

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do ii = iwarp,nblocknlb_pair-1,nwarp
           iblock= lst(2*ii+1)
           idx   = (iblock-1)*warpsize + ilane
           iscan = block_scan(iblock)*warpsize + 2*nblocknlb_pair
           if (iblock<nblock) then
              iscan1= block_scan(iblock+1)*warpsize + 2*nblocknlb_pair
           end if
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)

           kdx   = (lst(2*ii+2)-1)*warpsize + ilane
c          kglob = cell_glob(kdx)
           xk    = x(kdx)
           yk    = y(kdx)
           zk    = z(kdx)

           jlth    = warpsize
           jrth    = 1
           nlbits  = 0
           nrbits  = 0
           nlbits1 = 0
           nrbits1 = 0

           if (idx.eq.kdx) then
               lst(iscan+ilane) = kdx
              lst1(iscan+ilane) = kdx
           else
              do j = 1, warpsize
                 srclane = j
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
                 call image_inl(xpos,ypos,zpos)
                 r2      = xpos**2+ypos**2+zpos**2
                 ilane_set  = merge(1,0,
     &                       (r2.ge.cutbufbeg2.and.r2.le.cutbuff2))
                 ilane_set1 = merge(1,0,(r2.le.shortcutbuff2))
                 accept   = ballot(ilane_set)
                 accept1  = ballot(ilane_set1)

                 if (accept.ne.0) then
                    if (__popc(accept)>grp) then
                       ! count the number of validation
                       nlbits = nlbits + 1
                       ! find corresponding position for ilane
                       if (nlbits.eq.ilane) jlth = j
                    else
                       if (warpsize-nrbits.eq.ilane) jrth = j
                       nrbits = nrbits + 1
                    end if
                 end if
                 if (accept1.ne.0) then !Same operation for short list
                    if (__popc(accept1)>grp) then
                       nlbits1 = nlbits1 + 1
                       if (nlbits1.eq.ilane) jlth1 = j
                    else
                       if (warpsize-nrbits1.eq.ilane) jrth1 = j
                       nrbits1 = nrbits1 + 1
                    end if
                 end if
              end do

              ! Post processing
              if ((nlbits+nrbits)>0) then
                 if (ilane.eq.1)
     &              oldnbl= atomicAdd(nlst_l(iblock),nlbits)
                 if (ilane.eq.warpsize)
     &              oldnbr = atomicAdd(nlst_r(iblock),nrbits)

                 kdx_l  = __shfl(kdx,jlth)
                 kdx_r  = __shfl(kdx,jrth)
                 oldnbl = __shfl(oldnbl,1)
                 oldnbr = __shfl(oldnbr,warpsize)

                 if      (ilane.lt.nlbits+1) then
                    lst(iscan+oldnbl+ilane) = kdx_l
                 else if (ilane.gt.warpsize-nrbits) then
                    lst(iscan1-oldnbr+(ilane-warpsize)) = kdx_r
                 end if
              end if
              ! Post-processing short list
              if ((nlbits1+nrbits1)>0) then
                 if (ilane.eq.1)
     &              oldnbl= atomicAdd(nlst_l1(iblock),nlbits1)
                 if (ilane.eq.warpsize)
     &              oldnbr = atomicAdd(nlst_r1(iblock),nrbits1)

                 kdx_l  = __shfl(kdx,jlth1)
                 kdx_r  = __shfl(kdx,jrth1)
                 oldnbl = __shfl(oldnbl,1)
                 oldnbr = __shfl(oldnbr,warpsize)

                 if      (ilane.lt.nlbits1+1) then
                    lst1(iscan+oldnbl+ilane) = kdx_l
                    !if (iblock.eq.1) print*,kdx_l,'l'
                 else if (ilane.gt.warpsize-nrbits1) then
                    lst1(iscan1-oldnbr+(ilane-warpsize)) = kdx_r
                 end if
              end if
           end if
        end do
        end subroutine

        attributes(global)
     &  subroutine filter_lst_sparse_mpi
     &             (cell_glob,block_scan,x,y,z,
     &             matb_lst,nlocnlb,nblock,matsize,nblocknlb_pair,
     &             cutbuff2,xbeg_p,xend_p,ybeg_p,yend_p,zbeg_p,zend_p,
     &             lst,nlst,idx_glob)

        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nblocknlb_pair
        real(t_p),value,intent(in)  :: cutbuff2,xbeg_p,xend_p
     &                 ,ybeg_p,yend_p,zbeg_p,zend_p
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nlst(nblock),idx_glob(nlocnlb)
     &                ,lst(nblocknlb_pair*(BLOCK_SIZE+2))
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp,iwarp,ilane,ilane_mid
        integer iblock,j,srclane
        integer idx,kdx,kdx_,iscan,kglob,kneig,ii
        integer flag,jth,nbits,oldnb,accept,accept_mid,bset
        integer jflag,ilane_set
        real(t_p) xi,yi,zi,xk,yk,zk,xk_,yk_,zk_,xpos,ypos,zpos

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do ii = iwarp,nblocknlb_pair-1,nwarp
           iblock= lst(2*ii+1)
           idx   = (iblock-1)*warpsize + ilane
           iscan = block_scan(iblock)*warpsize + 2*nblocknlb_pair
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)

           kdx   = (lst(2*ii+2)-1)*warpsize + ilane
c          kglob = cell_glob(kdx)
           xk    = x(kdx)
           yk    = y(kdx)
           zk    = z(kdx)

           jth   = warpsize
c          bset  = 0
           nbits = 0

           do j = 1, warpsize
              srclane = j
              xk_     = __shfl(xk,srclane)
              yk_     = __shfl(yk,srclane)
              zk_     = __shfl(zk,srclane)
              xpos    = xi - xk_
              ypos    = yi - yk_
              zpos    = zi - zk_
              call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
              if  ((zk.lt.zbeg_p).or.(zk.ge.zend_p)
     &         .or.(yk.lt.ybeg_p).or.(yk.ge.yend_p)
     &         .or.(xk.lt.xbeg_p).or.(xk.ge.xend_p)) then
                 ilane_mid = 1
                 !Tag idx
                 idx  = ibset(idx,warpsize-1)
              else
                 ilane_mid = 0
              end if

              ilane_set = merge(1,0,
     &                    (xpos**2+ypos**2+zpos**2.lt.cutbuff2))
              accept  = ballot(ilane_set)
              accept_mid = ballot (ilane_mid)
              if (accept.ne.0) then
                 ! count the number of validation
                 nbits = nbits + 1
                 ! find corresponding position for ilane
                 if (nbits.eq.ilane) jth = j
              end if
              if (accept_mid.ne.0.and.ilane.eq.srclane) 
     &           kdx = ibset(kdx,warpsize-1) !Tag kdx
           end do

           if (abs(idx).eq.abs(kdx)) then
              lst (iscan+ilane) = kdx
           else
              if (ilane.eq.1) oldnb= atomicAdd(nlst(iblock),nbits)

              kdx_  = __shfl(kdx,jth)
              oldnb = __shfl(oldnb,1)

              if (ilane.lt.nbits+1) lst(iscan+oldnb+ilane) = kdx_
           end if

           if (idx.lt.0)
     &        oldnb = atomicOr (idx_glob(abs(idx)),ishft(1,warpsize-1))

        end do
        end subroutine
      end module
#else
      subroutine void_nblist_cu
      implicit none
      !Do Nothing
      end subroutine
#endif
