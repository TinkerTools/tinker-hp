#ifdef _CUDA
#include "tinker_precision.h"
#define TINKER_CUF
      module nblistcu

      use cudadevice
      use neigh     ,only: n_bl=>normal_block,i_bl=>inside_block
     &              ,d_bl=>disc_block
      use utilcu    ,only: all_lanes,xcell,ycell,zcell
     &              ,xcell2,ycell2,zcell2,mincell2,f_sqrt,f_abs
      use utilgpu   ,only: BLOCK_SIZE,WARP_SIZE,WARP_SHFT,inf
c     use tinheader ,only: ti_p,ti_eps
      integer,parameter::REORDER_BSIZE=64
      integer,parameter::grp=2
      integer,parameter::InMID_BS=128
      integer,parameter::sztb=128 !leading dimension of trackb array field from nblst_t

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
        integer idx,kdx,kdx_,iscan,kglob,it1,it2
        integer flag,jth,nbits,oldnb,oldnb_,accept,bset
        integer(8) it
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos
        integer ilane_set
        parameter( bit_si=32, bit_sh=5 )

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        !iwarp_b = (threadIdx%x-1) / warpsize
        !nwarp_b = blockDim%x / warpsize
        iwarp_b = (ithread-1) / warpsize
        nwarp_b = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do it = iwarp_b,(int(nblock,8)*(nblock+1)/2)-1,nwarp_b
           it2 = (sqrt(1.0d0+8.0d0*it)-1)/2
           it1 = it - int(it2,8)*(it2+1)/2

              pos2  = ishft(it2,-bit_sh)+1
              bit2  = iand (it2,bit_si-1)
              if (.not.btest(matb_lst(it1*matsize+pos2),bit2)) cycle
           !if (it1.ge.nblock-1.and.ilane.eq.1) then
           !   print*, int(it),it1,it2,nblock*(nblock+1)/2
           !end if
        !do it1 = blockIdx%x-1,nblock-1,gridDim%x
           idx   = it1*warpsize + ilane
           iscan = block_scan(it1+1)*warpsize + 2*nblocknlb_pair
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)

           !do it2 = it1+iwarp_b, nblock-1, nwarp_b
              kdx   = it2*warpsize + ilane
              !kglob = cell_glob(kdx)
              xk    = x(kdx)
              yk    = y(kdx)
              zk    = z(kdx)

              if (it1.eq.it2) then
                 lst (iscan+ilane) = kdx
              else

                 jth   = warpsize
                 nbits = 0
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
                    if (accept.ne.0) then
                       ! count the number of validation
                       nbits = nbits + 1
                       ! find corresponding position for ilane
                       if (nbits.eq.ilane) jth = j
                    end if
                 end do

                 if (ilane.eq.1) oldnb= atomicAdd(nlst(it1+1),nbits)

                 kdx_  = __shfl(kdx,jth)
                 oldnb = __shfl(oldnb,1)

                 if (ilane.lt.nbits+1) lst(iscan+oldnb+ilane) = kdx_
              end if
           !end do
        end do
        end subroutine


        ! Extract from block-block list the block-atoms list
        attributes(global)
     &  subroutine filter_lst_sparse
     &             (cell_glob,block_scan,x,y,z,b_stat,b_rmid
     &             ,matb_lst,nlocnlb,nblock,matsize,nb_p,cutbuff2
     &             ,lst,blst,alst,nlst,lst1,nlst1,nablst,trackb)

        !nb_p              number of block pairs
        !nablst            number of neighbor atoms per block in block-atom list
        !trackb            track block id reservation (in block-atom list) per block
        !blst alst & nlst  block-atom list & size   
        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nb_p
        real(t_p),value,intent(in)  :: cutbuff2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock),b_stat(*)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nablst(nblock)
     &                ,trackb(sztb*nblock),lst(nb_p*2)
     &                ,alst(nb_p*(BLOCK_SIZE)),blst(nb_p),nlst(4)
     &                ,lst1(nb_p),nlst1
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z
        real(t_p),device,intent(in) :: b_rmid(*)

        integer ithread,nwarp,iwarp,ilane
        integer iblock,kblock,j,srclane
        integer idx,kdx,kdx_,iscan,kglob,kneig,ii
        integer flag,jth,nbits,oldnb,oldnb_,accept,bset
        integer ilane_set,total,stat_i,inc
        integer r_onb,oldbl,newbl,idx_l,nthr_ces,left_ces,rght_ces
        integer(1) info
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / WARP_SIZE
        nwarp   = (blockDim%x*gridDim%x) / WARP_SIZE
        ilane   = iand( threadIdx%x-1,WARP_SIZE-1 ) + 1

        do ii = iwarp,nb_p-1,nwarp
           iblock= lst(2*ii+1)
           kblock= lst(2*ii+2)

           stat_i = b_stat(iblock)
c          if (iblock.gt.nblock.or.kblock.gt.nblock.and.ilane.eq.1)
c    &        print*,'filter_lst_sparse: out of block'
c    &              ,iblock,kblock,ii

           ! Filter block-block list in lst1
           block; real(t_p) d2,l2b
           if (iblock.eq.kblock) then
              if (ilane.eq.1) then
                 oldnb= atomicAdd(nlst1,1)
                 lst1(2*oldnb+1) = iblock
                 lst1(2*oldnb+2) = kblock
              end if
              goto 30
           else if (stat_i.eq.i_bl) then
              xpos = b_rmid(4*(iblock-1)+1) - b_rmid(4*(kblock-1)+1)
              ypos = b_rmid(4*(iblock-1)+2) - b_rmid(4*(kblock-1)+2)
              zpos = b_rmid(4*(iblock-1)+3) - b_rmid(4*(kblock-1)+3)
              l2b  = b_rmid(4*(iblock-1)+4) + b_rmid(4*(kblock-1)+4)
              d2   = xpos**2+ypos**2+zpos**2
              if (d2.lt.cutbuff2) then
                 if (ilane.eq.1) then
                    oldnb= atomicAdd(nlst1,1)
                    lst1(2*oldnb+1) = iblock
                    lst1(2*oldnb+2) = kblock
                 end if
                 goto 30
              else if (f_sqrt(d2).gt. l2b+f_sqrt(cutbuff2)) then
                 goto 30
              end if
           else if (stat_i.eq.n_bl) then
              if (b_stat(kblock).eq.d_bl) then
                 info=1
              else
                 xpos = b_rmid(4*(iblock-1)+1) - b_rmid(4*(kblock-1)+1)
                 ypos = b_rmid(4*(iblock-1)+2) - b_rmid(4*(kblock-1)+2)
                 zpos = b_rmid(4*(iblock-1)+3) - b_rmid(4*(kblock-1)+3)
                 l2b  = b_rmid(4*(iblock-1)+4) + b_rmid(4*(kblock-1)+4)
                 call image1_inl(xpos,ypos,zpos,info)
                 d2   = xpos**2+ypos**2+zpos**2
                 zpos = f_sqrt(cutbuff2)
                 if (d2.lt.cutbuff2) then
                    if (ilane.eq.1) then
                       oldnb= atomicAdd(nlst1,1)
                       lst1(2*oldnb+1) = iblock
                       lst1(2*oldnb+2) = kblock
                    end if
                    goto 30
                 else if (f_sqrt(d2).gt.l2b+zpos) then
                    goto 30
                 else if (2*(mincell2-l2b).lt.zpos) then
                    info=1
                 end if
              end if
           end if
           end block

           ! Filter block-atom list
           idx   = (iblock-1)*WARP_SIZE + ilane
           kdx   = (kblock-1)*WARP_SIZE + ilane
           !iscan = block_scan(iblock)*WARP_SIZE + 2*nb_p
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)
           xk    = x(kdx)
           yk    = y(kdx)
           zk    = z(kdx)

           jth   = WARP_SIZE
           nbits = 0

           if (stat_i.eq.i_bl.or.
     &         stat_i.eq.n_bl.and..not.info) then
              do j = 1, WARP_SIZE
                 srclane = j
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
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
           else
              do j = 1, WARP_SIZE
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
           end if

           if (nbits.eq.0) then
              !if (ilane.eq.1) then
              !   oldbl = atomicadd(nlst(5),1)
              !end if
              goto 30
           end if

           ! Append block-atom list

           if (ilane.eq.1) then
              block; integer d_onb,ltb,rtb
              !Update number of atoms neighbor to block i
              oldnb = atomicAdd(nablst(iblock),nbits)
              r_onb = iand(oldnb,WARP_SIZE-1)
              if (r_onb.eq.0) r_onb = WARP_SIZE
              d_onb = ishft(oldnb,-WARP_SHFT)
              if (d_onb.ge.sztb) print*,'filter_lst_sparse1: trackb to',
     &            'small',iblock,oldnb,d_onb

              !control reservation on nlst
              inc   = merge(1,0,(r_onb+nbits.gt.WARP_SIZE))
              if (inc) then
                 newbl = atomicadd(nlst(1),inc)                  ! Get next block
                 if (r_onb.eq.WARP_SIZE) then
                    ltb   = sztb*(iblock-1)+d_onb +1
                    !rtb   = ltb
                    oldbl = atomicexch(trackb(ltb),newbl+1)   ! Mark reserved block
                    if (oldbl.ne.0) print*, 'Tag Issue 0',oldbl
                    oldbl = newbl
                 else
                    ltb   = sztb*(iblock-1)+d_onb +1
                    !rtb   = ltb+1
                    oldbl = atomicexch(trackb(ltb+1),newbl+1)   ! Mark reserved block
                    if (oldbl.ne.0) print*, 'Tag Issue 1',oldbl
                    oldbl = atomicand(trackb(ltb),ALL_LANES)
                    do while( oldbl.eq.0 )
                       oldbl = atomicand(trackb(ltb),ALL_LANES)
                    end do
                    oldbl = oldbl-1
                 end if
c                if (rtb.ne.ltb+1.and.r_onb.ne.WARP_SIZE)
c    &              print*,'trouble',oldnb,newbl,ltb,rtb
c                if (rtb.ne.ltb.and.r_onb.eq.WARP_SIZE)
c    &              print*,'ptrouble',oldnb,nbits,newbl,ltb,rtb
                 if (oldbl.eq.-1) print*,'severe Tag Issue'
     &                                  ,iblock,r_onb,ltb
              else
                 ! Get previous block
                 ltb   = sztb*(iblock-1)+d_onb +1
                 oldbl = atomicand(trackb(ltb),ALL_LANES)
                 do while( oldbl.eq.0 )
                    oldbl = atomicand(trackb(ltb),ALL_LANES)
                 end do
                 oldbl = oldbl-1
                 newbl = oldbl
              end if
              if (inc) blst(newbl+1) = iblock   ! Append block list 
              !if (iblock.EQ.1) then
              !   print*,oldnb,nbits,inc.eq.1,d_onb,oldbl,newbl
              !end if
              end block
           end if

           kdx_  = __shfl(kdx,jth)
           r_onb = __shfl(r_onb,1)
           newbl = __shfl(newbl,1)
           oldbl = __shfl(oldbl,1)

           nthr_ces = WARP_SIZE - r_onb
           left_ces = (oldbl+1)*WARP_SIZE - nthr_ces + ilane
           rght_ces = (newbl  )*WARP_SIZE - nthr_ces + ilane
           idx_l    = merge(left_ces,rght_ces,ilane.le.nthr_ces)
c          if (ilane.eq.1.and.inc.and.iblock.eq.1) print*,
c    &           newbl,oldbl,oldnb,nbits,r_onb,left_ces
           if (ilane.lt.nbits+1) alst(idx_l) = kdx_

 30        continue
        end do
        end subroutine

        attributes(global)
     &  subroutine filter_lsts_sparse
     &             (cell_glob,block_scan,x,y,z,b_stat,b_rmid
     &             ,matb_lst,nlocnlb,nblock,matsize,nb_p
     &             ,scutbuff2,lcbuff2,cutbuff2
     &             ,lst,blst,alst,nlst,lst1,nlst1
     &                 ,alst_,nablst,trackb)

        !nb_p              number of block pairs
        !nablst            number of neighbor atoms per block in block-atom list
        !trackb            track block id reservation (in block-atom list) per block
        !blst alst & nlst  block-atom list & size   
        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nb_p
        real(t_p),value,intent(in)  :: scutbuff2,lcbuff2,cutbuff2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock),b_stat(*)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nablst(nblock)
     &                ,trackb(2*sztb*nblock),lst(nb_p*2)
     &                ,alst(nb_p*(BLOCK_SIZE)),blst(2*nb_p),nlst(8)
     &                ,lst1(2*nb_p),nlst1(2)
     &                ,alst_(nb_p*(BLOCK_SIZE))
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z
        real(t_p),device,intent(in) :: b_rmid(*)

        integer ithread,nwarp,iwarp,ilane
        integer iblock,kblock,j,srclane
        integer idx,kdx,kdx_,iscan,kglob,kneig,ii
        integer flag,jth,jth1,nbits,nbits1,nbits_,accept,accept1
        integer ilane_set,ilane_set1,stat_i,inc,oldnb,oldnb_
        integer r_onb,oldbl,newbl,idx_l,nthr_ces,left_ces,rght_ces
        integer ofl,oft,ofn
        integer(1) info
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos,d2
        logical add 

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / WARP_SIZE
        nwarp   = (blockDim%x*gridDim%x) / WARP_SIZE
        ilane   = iand( threadIdx%x-1,WARP_SIZE-1 ) + 1
        ofl     = merge(0,nb_p,       ilane.eq.1)
        ofn     = merge(0,nblock,     ilane.eq.1)
        oft     = merge(0,sztb*nblock,ilane.eq.1)

        do ii = iwarp,nb_p-1,nwarp
           iblock = lst(2*ii+1)
           kblock = lst(2*ii+2)

           stat_i = b_stat(iblock)
           add    = .false.
           ! Filter block-block list in lst1
           block; real(t_p) l2b;
           if      (iblock.eq.kblock) then
              if (ilane.lt.3) then
                 oldnb= atomicAdd(nlst1(ilane),1)
                 lst1(ofl+2*oldnb+1) = iblock
                 lst1(ofl+2*oldnb+2) = kblock
              end if
              goto 30
           else if (stat_i.eq.i_bl) then
              xpos = b_rmid(4*(iblock-1)+1) - b_rmid(4*(kblock-1)+1)
              ypos = b_rmid(4*(iblock-1)+2) - b_rmid(4*(kblock-1)+2)
              zpos = b_rmid(4*(iblock-1)+3) - b_rmid(4*(kblock-1)+3)
              l2b  = b_rmid(4*(iblock-1)+4) + b_rmid(4*(kblock-1)+4)
              d2   = xpos**2+ypos**2+zpos**2
              if (d2.lt.cutbuff2) then
                 add = merge(d2.gt.lcbuff2,d2.lt.scutbuff2,ilane.eq.1)
                 if (ilane.lt.3.and.add) then
                    oldnb= atomicAdd(nlst1(ilane),1)
                    lst1(ofl+2*oldnb+1) = iblock
                    lst1(ofl+2*oldnb+2) = kblock
                 end if
              else if (f_sqrt(d2).gt. l2b+f_sqrt(cutbuff2)) then
                 goto 30
              end if
           else if (stat_i.eq.n_bl) then
              if (b_stat(kblock).eq.d_bl) then
                 info=1
              else
                 xpos = b_rmid(4*(iblock-1)+1) - b_rmid(4*(kblock-1)+1)
                 ypos = b_rmid(4*(iblock-1)+2) - b_rmid(4*(kblock-1)+2)
                 zpos = b_rmid(4*(iblock-1)+3) - b_rmid(4*(kblock-1)+3)
                 l2b  = b_rmid(4*(iblock-1)+4) + b_rmid(4*(kblock-1)+4)
                 call image1_inl(xpos,ypos,zpos,info)
                 d2   = xpos**2+ypos**2+zpos**2
                 zpos = f_sqrt(cutbuff2)
                 if (d2.lt.cutbuff2) then
                   add = merge(d2.gt.lcbuff2,d2.lt.scutbuff2,ilane.eq.1)
                    if (ilane.lt.3.and.add) then
                       oldnb= atomicAdd(nlst1(ilane),1)
                       lst1(ofl+2*oldnb+1) = iblock
                       lst1(ofl+2*oldnb+2) = kblock
                    end if
                 else if (f_sqrt(d2).gt. l2b+zpos) then
                    goto 30
                 else if (2*(mincell2-l2b).lt.zpos) then
                    info=1
                 end if
              end if
           end if
           end block
           ! Jump to next iteration if i & k blocks are accepted in both lists
           if (iand(ballot(add),3).eq.3) then
              goto 30
           end if

           ! Filter block-atom list
           idx   = (iblock-1)*WARP_SIZE + ilane
           kdx   = (kblock-1)*WARP_SIZE + ilane
           !iscan = block_scan(iblock)*WARP_SIZE + 2*nb_p
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)
           xk    = x(kdx)
           yk    = y(kdx)
           zk    = z(kdx)

           jth   = WARP_SIZE; jth1  = WARP_SIZE
           nbits = 0; nbits1= 0

           if (stat_i.eq.i_bl.or.
     &         stat_i.eq.n_bl.and..not.info) then
              do j = 1, WARP_SIZE
                 srclane = j
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
                 d2      = xpos**2+ypos**2+zpos**2
                 ilane_set  = merge(1,0
     &                     ,(d2.gt.lcbuff2.and.d2.lt.cutbuff2))
                 ilane_set1 = merge(1,0,(d2.lt.scutbuff2))
                 accept  = ballot(ilane_set)
                 accept1 = ballot(ilane_set1)
                 if (accept.ne.0) then
                    ! count the number of validation
                    nbits = nbits + 1
                    ! find corresponding position for ilane
                    if (nbits.eq.ilane) jth = j
                 end if
                 if (accept1.ne.0) then
                    ! count the number of validation
                    nbits1 = nbits1 + 1
                    ! find corresponding position for ilane
                    if (nbits1.eq.ilane) jth1 = j
                 end if
              end do
           else
              do j = 1, WARP_SIZE
                 srclane = j
                 xpos    = xi - __shfl(xk,srclane)
                 ypos    = yi - __shfl(yk,srclane)
                 zpos    = zi - __shfl(zk,srclane)
                 call image_inl(xpos,ypos,zpos)
                 d2      = xpos**2+ypos**2+zpos**2
                 ilane_set  = merge(1,0
     &                     ,(d2.gt.lcbuff2.and.d2.lt.cutbuff2))
                 ilane_set1 = merge(1,0,(d2.lt.scutbuff2))
                 accept  = ballot(ilane_set )
                 accept1 = ballot(ilane_set1)
                 if (accept.ne.0) then
                    ! count the number of validation
                    nbits = nbits + 1
                    ! find corresponding position for ilane
                    if (nbits.eq.ilane) jth = j
                 end if
                 if (accept1.ne.0) then
                    ! count the number of validation
                    nbits1 = nbits1 + 1
                    ! find corresponding position for ilane
                    if (nbits1.eq.ilane) jth1 = j
                 end if
              end do
           end if

           if (nbits.eq.0.and.nbits1.eq.0) then
              !if (ilane.eq.1) then
              !   oldbl = atomicadd(nlst(5),1)
              !end if
              goto 30
           end if

           nbits_ = merge(nbits,nbits1,ilane.eq.1)
           add    = (.not.add).and.(nbits_.gt.0).and.ilane.lt.3

           ! Append block-atom list

           if (add) then
              block; integer d_onb,ltb,rtb
              !Update number of atoms neighbor to block i
              oldnb = atomicAdd(nablst(ofn+iblock),nbits_)
              r_onb = iand(oldnb,WARP_SIZE-1)
              if (r_onb.eq.0) r_onb = WARP_SIZE
              d_onb = ishft(oldnb,-WARP_SHFT)
              if (d_onb.ge.sztb) print*,': trackb to',
     &            'small',iblock,oldnb,d_onb,ilane
           !if (iblock.eq.81.and.ilane.lt.3) then
           !   print*,oldnb,ilane,nbits_,kblock
           !end if

              !control reservation on nlst
              inc   = merge(1,0,(r_onb+nbits_.gt.WARP_SIZE))
              if (inc) then
                 newbl = atomicadd(nlst(ilane),inc)                  ! Get next block
                 if (r_onb.eq.WARP_SIZE) then
                    ltb   = sztb*(iblock-1)+d_onb +1
                    !rtb   = ltb
                    oldbl = atomicexch(trackb(oft+ltb),newbl+1)   ! Mark reserved block
                    if (oldbl.ne.0) print*,'Tag Issue 0',ilane,oft,oldbl
                    oldbl = newbl
                 else
                    ltb   = sztb*(iblock-1)+d_onb +1
                    !rtb   = ltb+1
                    oldbl = atomicexch(trackb(oft+ltb+1),newbl+1)   ! Mark reserved block
                    if (oldbl.ne.0) print*,'Tag Issue 1',ilane,oft,oldbl
                    oldbl = atomicand(trackb(oft+ltb),ALL_LANES)
                    do while( oldbl.eq.0 )
                       oldbl = atomicand(trackb(oft+ltb),ALL_LANES)
                    end do
                    oldbl = oldbl-1
                 end if
c                if (rtb.ne.ltb+1.and.r_onb.ne.WARP_SIZE)
c    &              print*,'trouble',oldnb,newbl,ltb,rtb
c                if (rtb.ne.ltb.and.r_onb.eq.WARP_SIZE)
c    &              print*,'ptrouble',oldnb,nbits_,newbl,ltb,rtb
                 if (oldbl.eq.-1) print*,'severe Tag Issue'
     &                                  ,iblock,r_onb,ltb
              else
                 ! Get previous block
                 ltb   = sztb*(iblock-1)+d_onb +1
                 oldbl = atomicand(trackb(oft+ltb),ALL_LANES)
                 do while( oldbl.eq.0 )
                    oldbl = atomicand(trackb(oft+ltb),ALL_LANES)
                 end do
                 oldbl = oldbl-1
                 newbl = oldbl
              end if
              if (inc) blst(ofl+newbl+1) = iblock   ! Append block list 
              !if (iblock.EQ.1) then
              !   print*,oldnb,nbits_,inc.eq.1,d_onb,oldbl,newbl
              !end if
              end block
           end if

           ! Add Long List
           block
           integer r_onb_,newbl_,oldbl_
           if (__shfl(merge(1,0,add),1)) then
           kdx_  = __shfl(kdx,jth)
           r_onb_= __shfl(r_onb,1)
           newbl_= __shfl(newbl,1)
           oldbl_= __shfl(oldbl,1)

           nthr_ces = WARP_SIZE - r_onb_
           left_ces = (oldbl_+1)*WARP_SIZE - nthr_ces + ilane
           rght_ces = (newbl_  )*WARP_SIZE - nthr_ces + ilane
           idx_l    = merge(left_ces,rght_ces,ilane.le.nthr_ces)
c          if (ilane.eq.1.and.inc.and.iblock.eq.1) print*,
c    &           newbl_,oldbl_,oldnb,nbits,r_onb_,left_ces
           if (ilane.lt.nbits+1) alst(idx_l) = kdx_
           end if

           ! Add Short List
           if (__shfl(merge(1,0,add),2)) then
           kdx_  = __shfl(kdx,jth1)
           r_onb_= __shfl(r_onb,2)
           newbl_= __shfl(newbl,2)
           oldbl_= __shfl(oldbl,2)

           nthr_ces = WARP_SIZE - r_onb_
           left_ces = (oldbl_+1)*WARP_SIZE - nthr_ces + ilane
           rght_ces = (newbl_  )*WARP_SIZE - nthr_ces + ilane
           idx_l    = merge(left_ces,rght_ces,ilane.le.nthr_ces)
c          if (ilane.eq.1.and.inc.and.iblock.eq.1) print*,
c    &           newbl_,oldbl_,oldnb,nbits1,r_onb_,left_ces
           if (ilane.lt.nbits1+1) alst_(idx_l) = kdx_
           end if
           end block

 30        continue
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
        integer ilane_set,total,total1
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
              total =0
              total1=0
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
                 if (accept.EQ.all_lanes) then
                    total  = total +1
                 end if
                 total1 = total1+__popc(accept)
              end do
              !if (lst(2*ii+2).eq.iblock+1.and.ilane.eq.1) then
              !   print*,iblock,total,total1,'ii'
              !end if

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
        integer ilane_set,ilane_set1
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
     &  subroutine filter_lst_sparse2
     &             (cell_glob,block_scan,x,y,z,
     &             matb_lst,nlocnlb,nblock,matsize,nblocknlb_pair,
     &             cutbuff2,lst,abpl,nlst_l,nlst_r)

        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nblocknlb_pair
        real(t_p),value,intent(in)  :: cutbuff2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nlst_l(nblock),nlst_r(nblock)
     &         ,lst(nblocknlb_pair*2),abpl(nblocknlb_pair*BLOCK_SIZE)
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp,iwarp,ilane
        integer iblock,j,srclane
        integer idx,kdx,kdx_l,kdx_r,iscan,iscan1
        integer nlbits,nrbits,kglob,kneig,ii
        integer flag,jlth,jrth,oldnbl,oldnbr,accept,bset
        integer ilane_set,total,total1
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do ii = iwarp,nblocknlb_pair-1,nwarp
           iblock= lst(2*ii+1)
           idx   = (iblock-1)*warpsize + ilane
           iscan = block_scan(iblock)*warpsize
           if (iblock<nblock) then
              iscan1= block_scan(iblock+1)*warpsize
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
              total =0
              total1=0
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
                 if (accept.EQ.all_lanes) then
                    total  = total +1
                 end if
                 total1 = total1+__popc(accept)
              end do
              !if (lst(2*ii+2).eq.iblock+1.and.ilane.eq.1) then
              !   print*,iblock,total,total1,'ii'
              !end if

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
                    abpl(iscan+oldnbl+ilane) = kdx_l
                    !if (iblock.eq.1) print*,kdx_l,'l'
                 else if (ilane.gt.warpsize-nrbits) then
                    abpl(iscan1-oldnbr+(ilane-warpsize)) = kdx_r
                 end if
              end if
           end if
        end do
        end subroutine

c
c       "filter_shortlst_sparse2"
c       Filter both nblist and shortnblist with the same adjacency matrix
c
        attributes(global)
     &  subroutine filter_shortlst_sparse2
     &             (cell_glob,block_scan,x,y,z,
     &              matb_lst,nlocnlb,nblock,matsize,nblocknlb_pair,
     &              shortcutbuff2,cutbufbeg2,cutbuff2,
     &              lst,abpl,nlst_l,nlst_r,abpl_1,nlst_l1,nlst_r1)

        implicit none
        integer,value ,intent(in)   :: nblock,matsize,nlocnlb
     &                ,nblocknlb_pair
        real(t_p),value,intent(in)  :: cutbuff2,shortcutbuff2,cutbufbeg2
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: nlst_l(nblock),nlst_r(nblock)
     &         ,lst(nblocknlb_pair*2),abpl(nblocknlb_pair*BLOCK_SIZE)
        integer,device,intent(inout):: nlst_l1(nblock),nlst_r1(nblock)
     &                ,abpl_1(nblocknlb_pair*BLOCK_SIZE)
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp,iwarp,ilane
        integer iblock,j,srclane
        integer idx,kdx,kdx_l,kdx_r,iscan,iscan1
        integer nlbits,nrbits,nlbits1,nrbits1,kglob,kneig,ii
        integer flag,jlth,jrth,jlth1,jrth1,oldnbl,oldnbr
     &         ,oldnbl1,oldnbr1,accept,accept1
        integer ilane_set,ilane_set1
        real(t_p) r2,xi,yi,zi,xk,yk,zk,xpos,ypos,zpos

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do ii = iwarp,nblocknlb_pair-1,nwarp
           iblock= lst(2*ii+1)
           idx   = (iblock-1)*warpsize + ilane
           iscan = block_scan(iblock)*warpsize
           if (iblock<nblock) then
              iscan1= block_scan(iblock+1)*warpsize
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
              abpl  (iscan+ilane) = kdx
              abpl_1(iscan+ilane) = kdx
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
                    abpl(iscan+oldnbl+ilane) = kdx_l
                 else if (ilane.gt.warpsize-nrbits) then
                    abpl(iscan1-oldnbr+(ilane-warpsize)) = kdx_r
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
                    abpl_1(iscan+oldnbl+ilane) = kdx_l
                    !if (iblock.eq.1) print*,kdx_l,'l'
                 else if (ilane.gt.warpsize-nrbits1) then
                    abpl_1(iscan1-oldnbr+(ilane-warpsize)) = kdx_r
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
        integer ilane_set
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

!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       Filter C1 list and extract pairwise atom nblist in a (npairs,2)
!       pattern
        attributes(global)
     &  subroutine filter_pairwise_atom_lst
     &             (cell_glob,block_scan,x,y,z
     &             ,matb_lst,n,nlocnlb,nblock,matsize,nb_pair
     &             ,cutbuff2,lst,blst,xbeg,xend,ybeg,yend,zbeg,zend)

        implicit none
        integer  ,value,intent(in)  :: nblock,matsize,n,nlocnlb
     &                 ,nb_pair
        real(t_p),value,intent(in)  :: cutbuff2,xbeg,xend,ybeg,yend
     &                 ,zbeg,zend
        integer,device,intent(in)   :: cell_glob(nlocnlb)
        integer,device,intent(in)   :: block_scan(nblock)
        integer,device,intent(in)   :: matb_lst(matsize*nblock)
        integer,device,intent(inout):: blst(nb_pair*2)
     &                ,lst(2*nb_pair*(BLOCK_SIZE**2))
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer ithread,nwarp,iwarp,ilane
        integer iblock,j,srclane
        integer idx,kdx,kdx_,iscan,kglob,kneig,ii
        integer flag,jth,nbits,oldnb,oldnb_,accept,bset
        integer pair_a,ilane_set
        real(t_p) xi,yi,zi,xk,yk,zk,xk_,yk_,zk_,xpos,ypos,zpos
        logical accept_mid

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        pair_a  = nb_pair*(BLOCK_SIZE**2)

        do ii = iwarp,nb_pair-1,nwarp
           iblock= blst(2*ii+1)
           idx   = (iblock-1)*warpsize + ilane
           !iscan = block_scan(iblock)*(warpsize**2)
           iscan = ii*BLOCK_SIZE**2
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)

           kdx   = (blst(2*ii+2)-1)*warpsize + ilane
c          kglob = cell_glob(kdx)
           xk    = x(kdx)
           yk    = y(kdx)
           zk    = z(kdx)

           jth   = warpsize
           nbits = 0

           !if (idx.eq.kdx) then
           !   do j = 1,warpsize
           !      kdx_ = kdx - ilane + j
           !      if (kdx_.gt.idx.and.kdx_.lt.n+1) then
           !         lst(       iscan+(j-1)*warpsize+ilane) = idx-1
           !         lst(pair_a+iscan+(j-1)*warpsize+ilane) = kdx_-1
           !      end if
           !   end do
           !else
              do j = 1, warpsize
                 srclane = j
                 kdx_    = kdx - ilane + j
                 xk_     = __shfl(xk,srclane)
                 yk_     = __shfl(yk,srclane)
                 zk_     = __shfl(zk,srclane)
                 xpos    = xi - xk_
                 ypos    = yi - yk_
                 zpos    = zi - zk_
                 call midpointimage_inl( xk_,yk_,zk_, xpos,ypos,zpos )
                 if ((zk_.lt.zbeg).or.(zk_.ge.zend)
     &           .or.(yk_.lt.ybeg).or.(yk_.ge.yend)
     &           .or.(xk_.lt.xbeg).or.(xk_.ge.xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
                 if (idx.eq.kdx) then ! Remove half
                    accept_mid = accept_mid.and.(kdx_.gt.idx)
                 end if
                 if ((xpos**2+ypos**2+zpos**2.lt.cutbuff2)
     &               .and.accept_mid.and.kdx_.lt.n+1) then
                    lst(       iscan+(j-1)*warpsize+ilane) = idx-1
                    lst(pair_a+iscan+(j-1)*warpsize+ilane) = kdx_-1
                 end if
c                ilane_set = merge(1,0,
c    &                       (xpos**2+ypos**2+zpos**2.lt.cutbuff2))
c                accept  = ballot(ilane_set)
c                if (accept.ne.0) then
c                   ! count the number of validation
c                   nbits = nbits + 1
c                   ! find corresponding position for ilane
c                   if (nbits.eq.ilane) jth = j
c                end if
              end do

c             if (ilane.eq.1) oldnb= atomicAdd(nlst(iblock),nbits)

c             kdx_  = __shfl(kdx,jth)
c             oldnb = __shfl(oldnb,1)

c             if (ilane.lt.nbits+1) lst(iscan+oldnb+ilane) = kdx_

           !end if
        end do
        end subroutine

        attributes(global)
     &  subroutine blocksat_kcu(a_kind,a_id,s_key,x,y,z,nb,na,nab,n
     &            ,cutbuf,len_xcell,len_ycell,len_zcell,eps_cell
     &            ,s_kind,sgl_id,so_x,so_y,so_z,b_stat,b_rmid)
        implicit none
        integer  ,value :: nb,na,nab,n
        real(t_p),value :: cutbuf,len_xcell,len_ycell,len_zcell
     &           ,eps_cell
        !integer  ,device,intent(in)::
        real(t_p),device,intent(in ):: x(*),y(*),z(*)
        real(t_p),device,intent(out):: so_x(*),so_y(*),so_z(*)
        integer  ,device:: b_stat(  nb),a_kind(*),a_id(*)
        integer  ,device:: s_key(*),s_kind(*),sgl_id(*)
        real(t_p),device:: b_rmid(4*nb)

        integer ithread,iwarp,owarp,nwarp,ilane,ib,idx
        integer ikind,iglob,attr
        real(t_p) xr,yr,zr,xs,ys,zs,xmid,ymid,zmid
        logical isIn,travers
        !real(t_p),shared:: minxr(InMID_BS),maxxr(InMID_BS)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1)            / warpsize
        owarp   = (threadIdx%x-1)        / warpsize
        nwarp   = (blockDim%x*gridDim%x) / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        ! Compute middle point of each block
        do ib = iwarp, nb-1, nwarp
           idx   = ib*BLOCK_SIZE + ilane
           idx   = merge(idx,na,idx.le.na)

           ikind = a_kind(s_key(idx)) 
           iglob = a_id (ikind)
           xr    = x    (iglob)
           yr    = y    (iglob)
           zr    = z    (iglob)
           xs=xr; ys=yr; zs=zr;  !Save position

           call image_inl(xr,yr,zr)
           if ((xcell2-xr).lt.eps_cell) xr = xr-0.05*len_xcell
           if ((ycell2-yr).lt.eps_cell) yr = yr-0.05*len_ycell
           if ((zcell2-zr).lt.eps_cell) zr = zr-0.05*len_zcell

           block;
           real(t_p) xlow,xhig,ylow,yhig,zlow,zhig,r2,r
           integer it,it1
           xlow = xr; xhig = xr;
           ylow = yr; yhig = yr;
           zlow = zr; zhig = zr;

           ! Min Max Warp Reduction
           it=1
           do while( it.lt.WARP_SIZE )
              xlow = min( xlow,__shfl_xor(xlow,it+1) )
              ylow = min( ylow,__shfl_xor(ylow,it+1) )
              zlow = min( zlow,__shfl_xor(zlow,it+1) )
              xhig = max( xhig,__shfl_xor(xhig,it+1) )
              yhig = max( yhig,__shfl_xor(yhig,it+1) )
              zhig = max( zhig,__shfl_xor(zhig,it+1) )
              it = it*2
           end do
           ! Extract block Mid point
           xmid = ( xlow+xhig )*0.5_ti_p
           ymid = ( ylow+yhig )*0.5_ti_p
           zmid = ( zlow+zhig )*0.5_ti_p

           ! Compute block Radius
           r2 = (xr-xmid)**2 + (yr-ymid)**2 + (zr-zmid)**2
           r  = f_sqrt(r2)
           it = 1
           !do while( it.lt.WARP_SIZE )
           !   r  = max( r,__shfl_xor(r,it+1) )
           !   it = it*2
           !end do
           r  = max( r,__shfl_xor(r, 1+1) )
           r  = max( r,__shfl_xor(r, 2+1) )
           r  = max( r,__shfl_xor(r, 4+1) )
           r  = max( r,__shfl_xor(r, 8+1) )
           r  = max( r,__shfl_xor(r,16+1) )
           !r  = max( r,__shfl_xor(r,32+1) )

           ! Check if block is deep inside the box
           isIn =  (xlow+xcell2.gt.cutbuf.and.xcell2-xhig.gt.cutbuf
     &        .and. ylow+ycell2.gt.cutbuf.and.ycell2-yhig.gt.cutbuf
     &        .and. zlow+zcell2.gt.cutbuf.and.zcell2-zhig.gt.cutbuf)

           ! Check if block volume is discontinuous
           it  = ballot( xr+xcell2.lt.len_xcell )
           it1 = ballot( xcell2-xr.lt.len_xcell )
           travers = merge( .true.,.false.,it.ne.0.and.it1.ne.0 )

           ! Check if block is irregular (use d_bl tag for it)
           xr=xs; yr=ys; zr=zs;
              ! Compute dist & warp comparison
           xr = f_abs( xr - __shfl(xs,WARP_SIZE/2) )
           yr = f_abs( yr - __shfl(ys,WARP_SIZE/2) )
           zr = f_abs( zr - __shfl(zs,WARP_SIZE/2) )
           it = ballot(xr.gt.xcell2)
           it = ior( it,ballot(yr.gt.ycell2) )
           it = ior( it,ballot(zr.gt.zcell2) )
           travers = travers.or.(it.ne.0)

           attr = merge(d_bl,merge(i_bl,n_bl,isIn),travers)

           if (attr.ne.d_bl) then
           ! Extract block physical center
           xlow = xs; xhig = xs;
           ylow = ys; yhig = ys;
           zlow = zs; zhig = zs;

           ! Min Max Warp Reduction
           it=1
           do while( it.lt.WARP_SIZE )
              xlow = min( xlow,__shfl_xor(xlow,it+1) )
              ylow = min( ylow,__shfl_xor(ylow,it+1) )
              zlow = min( zlow,__shfl_xor(zlow,it+1) )
              xhig = max( xhig,__shfl_xor(xhig,it+1) )
              yhig = max( yhig,__shfl_xor(yhig,it+1) )
              zhig = max( zhig,__shfl_xor(zhig,it+1) )
              it = it*2
           end do
           ! Extract block Mid point
           xmid = ( xlow+xhig )*0.5_ti_p
           ymid = ( ylow+yhig )*0.5_ti_p
           zmid = ( zlow+zhig )*0.5_ti_p
           end if

           ! Load block data
           if (ilane.eq.1) then
              b_stat   (ib+1)    = attr
              b_rmid(4*(ib  )+1) = xmid
              b_rmid(4*(ib  )+2) = ymid
              b_rmid(4*(ib  )+3) = zmid
              b_rmid(4*(ib  )+4) = r
           end if
           end block

           if (idx.gt.na) then
              ikind = nab; iglob = n
              xs = inf; ys = inf; zs = inf;
           end if
           ! Spatial ordering of positions and global index buffers
           s_kind(idx) = ikind
           sgl_id(idx) = iglob
           so_x  (idx) = xs
           so_y  (idx) = ys
           so_z  (idx) = zs

        end do
        end subroutine

!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       Reorder list
        attributes(global)
     &  subroutine reorder_list_block(ilst,lst,x,y,z,
     &             nb_pair,nblock,nlocnlb,cut2,cutbuff2)

        implicit none
        integer  ,value,intent(in)  :: nblock,nb_pair,nlocnlb
        real(t_p),value,intent(in)  :: cut2,cutbuff2
        integer,device,intent(inout):: ilst(nb_pair)
     &                ,lst(nb_pair*(BLOCK_SIZE))
        real(t_p),device,intent(in),dimension(nlocnlb)::x,y,z

        integer,shared:: atset(REORDER_BSIZE),accept(REORDER_BSIZE)
        integer,shared:: atsave(REORDER_BSIZE)
        integer,shared:: accloc(REORDER_BSIZE/warpsize)
        integer,shared:: rejloc(REORDER_BSIZE/warpsize)
        integer,shared:: accscan(REORDER_BSIZE/warpsize)
        integer,shared:: rejscan(REORDER_BSIZE/warpsize)
        integer,shared:: acctot,nal_blocks

        integer ilane,klane,ithread,iwarp,nwarp,nwarp_b,bwarp,srclane
        integer idx,kdx,ibs,nal_block1
        integer i,j,k,ib,iblock,acceptIdxst,ial_blocks
        integer accpr,iwrite,nrecover,offset
        real(t_p) xi,yi,zi,xk,yk,zk,xpos,ypos,zpos

        parameter(nwarp_b=REORDER_BSIZE/warpsize)
        ithread = threadIdx%x
        !iwarp   = (ithread-1) / warpsize
        !nwarp   = (blockDim%x*gridDim%x) / warpsize
        bwarp   = (threadIdx%x-1)/warpsize + 1
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        do ib = blockIdx%x,nblock,gridDim%x
           ! Deal with global
           if (ithread.eq.1) nal_blocks  = 0
           if (REORDER_BSIZE>warpsize) call syncthreads
           nal_block1 = 0

           ibs   = ib
           iblock= ilst(ibs)
           !Search index to the first block such that (ib==iblock)
           do while(iblock.lt.ib)
              ibs    = ibs + 1
              iblock = ilst(ibs)
           end do
           !Gather data on processing i block
           idx   = (ib-1)*warpsize + ilane
           xi    = x(idx)
           yi    = y(idx)
           zi    = z(idx)
           !k     = 0

           ! Expand neighbor block to every warp
           acceptIdxst = (ibs-1)*warpsize
           ibs    = ibs + bwarp -1

           ! Loop over ib neighbor blocks
           do while(nal_block1.eq.0)

              iblock = ilst(ibs)
              if (iblock.ne.ib.and.ilane.eq.1) i=atomicAdd(nal_blocks,1)
              !if (REORDER_BSIZE>warpsize) call syncthreads

              ! Reset global values
              if (ilane.eq.1) then
                 acctot     = 0
                 accscan(1) = 1
                 rejscan(1) = 1
                 accloc(bwarp) = 0
                 rejloc(bwarp) = 0
              end if
              accept(threadIdx%x) = 0
              !atset(threadIdx%x)  = -200
              !atsave(threadIdx%x) = -250

              if (iblock.eq.ib) then
              !Gather data on k block
              kdx   = lst( (ibs-1)*warpsize + ilane )
              xk    = x(kdx)
              yk    = y(kdx)
              zk    = z(kdx)

              if (idx.eq.kdx) then
                 accept(threadIdx%x) = 1
                 accpr=warpsize
              else
                 do j = 0, warpsize-1
                    srclane = iand( ilane-1+j,warpsize-1 ) + 1
                    klane   = threadIdx%x - ilane + srclane
                    xpos    = xi - __shfl(xk,srclane)
                    ypos    = yi - __shfl(yk,srclane)
                    zpos    = zi - __shfl(zk,srclane)
                    call image_inl(xpos,ypos,zpos)
                    if (xpos**2+ypos**2+zpos**2.lt.cut2) then
                       !accept(klane) = ior(accept(klane),ishft(1,ilane-1))
                       accept(klane) = 1
                    end if
                 end do
                 accpr=0
                 do j = 0, warpsize-1
                    srclane = iand( ilane-1+j,warpsize-1 ) + 1
                    klane   = threadIdx%x - ilane + srclane
                    if (accept(klane).ne.0) accpr=accpr+1
                 end do
              end if
              if (ilane.eq.1) then
                 accloc(bwarp)= accpr
                 rejloc(bwarp)= warpsize-accpr     
              end if
              end if

              if (REORDER_BSIZE>warpsize) call syncthreads
              if (ilane.eq.2.and.iblock.eq.ib) i=atomicAdd(acctot,accpr)
              if (ithread.eq.1) then
                 do j = 2, nwarp_b
                    accscan(j) = accscan(j-1) + accloc(j-1)
                    rejscan(j) = rejscan(j-1) + rejloc(j-1)
                 end do
              end if

              if (REORDER_BSIZE>warpsize) call syncthreads
              if (iblock.eq.ib) then
              if (ilane.eq.2) accloc(bwarp) = acctot
              call syncwarp(all_lanes)
              !Where to write reordered kdx
              if (accloc(bwarp).eq.blockDim%x) then
                 iwrite = ithread
              else if (accept(threadIdx%x).ne.0) then
                 iwrite = atomicAdd(accscan(bwarp),1)
              else
                 iwrite = atomicAdd(rejscan(bwarp),1) + accloc(bwarp)
              end if
              atset(iwrite) = kdx
              end if
              nal_block1 = nal_blocks ! save to register

              if (REORDER_BSIZE>warpsize) call syncthreads
              if (iblock.eq.ib) then
              ! Global reordering in lst
              iwrite = (ibs-bwarp)*warpsize
              if (iwrite.eq.acceptIdxst) then
                 lst(acceptIdxst+threadIdx%x) = atset(threadidx%x)
              else
                 if (ithread.le.accloc(bwarp)) then
                    atsave(ithread) = lst(acceptIdxst+ithread)
                    lst(acceptIdxst+ithread) = atset(ithread)
                    ! atoms to recover from list
                    nrecover=min(accloc(bwarp),iwrite-acceptIdxst)
                    if (ithread.le.nrecover) then
                       offset = max(iwrite,acceptIdxst+accloc(bwarp))
                       lst(offset+ithread)=atsave(ithread)
                    end if
                 else
                    ! Atoms to push at the end (Those rejected)
                    lst((ibs-1)*warpsize+ilane) = atset(threadIdx%x)
                 end if
              end if
              acceptIdxst = acceptIdxst + accloc(bwarp)
              end if

c             if (ib.eq.1.and.ithread.eq.1) then
c                print*,ibs,accloc(bwarp),acceptIdxst
c    &                 ,acceptIdxst-accloc(bwarp),iwrite,nal_block1
c             end if
              !if (REORDER_BSIZE>warpsize) call syncthreads

              !Increment index
              ibs = ibs + nwarp_b
              !k   = k +1
           end do
        end do
        end subroutine

      end module
#else
      subroutine void_nblist_cu
      implicit none
      !Do Nothing
      end subroutine
#endif
