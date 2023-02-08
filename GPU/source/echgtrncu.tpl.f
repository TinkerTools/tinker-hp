c
c     Chgtrn CUDA Kernel Template 
c
c     "echgtrncu" calculates the charge transfer energy and first
c     derivatives with respect to Cartesian coordinates on device
c     version __tver__ features __tfea__
c
      attributes(global) subroutine __Cat(echgtrn,__sufx__)
     &        ( ipole, pglob, loc, bapl, abpl, x, y, z, chgct, dmpct
     &        , grplist, wgrp
     &        , dect, e_buff, v_buff, n_buff
     &        , na, nb, nab, nbap, ctrntyp, ugrp
     &        , f, ctrnscut, loff2, off2, off, cut2, sheal, rinv
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &        , ctscal_ik, ctscal_val, ipole_, loc_, x_, y_, z_
     &        , n_ctscale
     &        )
      implicit none
      integer  ,intent(in),value:: na,nb,nab,nbap,ctrntyp,n_ctscale
      logical  ,intent(in),value:: ugrp
      real(t_p),intent(in),value:: f,ctrnscut,loff2,off2,off,cut2,sheal
     &         ,rinv,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      integer  ,intent(in),device:: ipole(*),pglob(*),loc(*),bapl(*)
     &         ,abpl(*),grplist(*),ctscal_ik(2,*),ipole_(*)
     &         ,loc_(*)
      real(t_p),intent(in),device:: x(*),y(*),z(*),chgct(*),dmpct(*)
     &         ,wgrp(maxgrp+1,*),x_(*),y_(*),z_(*),ctscal_val(*)
      integer  ,device:: n_buff(*)
      ener_rtyp,device:: e_buff(RED_BUFF_SIZE)
      mdyn_rtyp,device:: dect(3,*)
      real(t_p),device:: v_buff(*)
c
      integer   ithread,iwarp,nwarp,ilane,klane,srclane,lot
      integer   ii,j,i,kbis,iblock,idx,kdx,ver,fea
      integer   iipole,iglob,iploc,kpole,kglob,kploc,nect_
      real(t_p) xk_,yk_,zk_,d2
      real(t_p) vir_(6),e,fgrp
      ener_rtyp ect_
      real(t_p) chgi,chgk,alphai,alphak
      type(real3) posi,pos,frc
      type(real3),shared:: posk(BLOCK_DIM)
      type(mdyn3_r) frc_i,frc_k
      logical   do_pair,same_block,accept_mid
      parameter(ver=__tver__, fea=__tfea__)
c
      ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
      iwarp   = ((ithread-1) / warpsize) + 1
      nwarp   = blockDim%x*gridDim%x / warpsize
      ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
      accept_mid = .true.
c
c     Correct Scaling Interactions
c
      do ii = ithread, n_ctscale, blockDim%x*gridDim%x
      block;
      real(t_p) mscale
      integer kkpole
         iipole = ctscal_ik(1,ii)
         kkpole = ctscal_ik(2,ii)
         mscale =-ctscal_val(ii)
         iglob  = ipole_(iipole)
         kglob  = ipole_(kkpole)
         chgi   = chgct (iipole)
         alphai = dmpct (iipole)
         chgk   = chgct (kkpole)
         alphak = dmpct (kkpole)

         pos%x  = x_(kglob) - x_(iglob)
         pos%y  = y_(kglob) - y_(iglob)
         pos%z  = z_(kglob) - z_(iglob)
         call image_inl(pos%x,pos%y,pos%z)
         d2      = pos%x**2 + pos%y**2 + pos%z**2
#if __tfea__ & __use_longRange__
         if (loff2.gt.d2.or.d2.gt.off2) cycle  !apply cutoff
#else
         if (d2.gt.off2) cycle  !apply cutoff
#endif
         if (alphai.eq.0.0) alphai=1000.0
         if (alphak.eq.0.0) alphak=1000.0
#if __tver__ & __use_ene__
         ect_ = 0
#endif
#if __tfea__ & __use_groups__
         if (ugrp) call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
#endif
         ! compute mpole one interaction
         call duo_chgtrn(d2,pos%x,pos%y,pos%z,chgi,chgk,alphai,alphak
     &           ,f,mscale,fgrp,rinv,off,ctrnscut,cut2,sheal
     &           ,ctrntyp,ugrp,e,frc,ver+__use_sca__,fea)

         lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
#if __tver__ & __use_ene__
         ! update energy
         call atomic_add_f( e_buff(lot),tp2enr(e) )
#endif
#if __tver__ & __use_act__
         if (mscale.eq.-oner) call atomic_sub_i(n_buff(lot), 1)
#endif
#if __tver__ & __use_vir__
         ! increment the virial due to pairwise Cartesian forces c
         call atomic_add_c(v_buff(0*RED_BUFF_SIZE+lot),pos%x*frc%x)
         call atomic_add_c(v_buff(1*RED_BUFF_SIZE+lot),pos%y*frc%x)
         call atomic_add_c(v_buff(2*RED_BUFF_SIZE+lot),pos%z*frc%x)
         call atomic_add_c(v_buff(3*RED_BUFF_SIZE+lot),pos%y*frc%y)
         call atomic_add_c(v_buff(4*RED_BUFF_SIZE+lot),pos%z*frc%y)
         call atomic_add_c(v_buff(5*RED_BUFF_SIZE+lot),pos%z*frc%z)
#endif
#if __tver__ & __use_grd__
         i      = loc_(iglob)
         kbis   = loc_(kglob)
         ! increment force-based gradient and torque on first site
         call atomic_add_f1( dect(1,i)    ,-frc%x )
         call atomic_add_f1( dect(2,i)    ,-frc%y )
         call atomic_add_f1( dect(3,i)    ,-frc%z )
         call atomic_add_f1( dect(1,kbis) , frc%x )
         call atomic_add_f1( dect(2,kbis) , frc%y )
         call atomic_add_f1( dect(3,kbis) , frc%z )
#endif
      end block
      end do
c
c       Compute pairwise interaction
c
      do ii = iwarp, nbap, nwarp
      block;
         iblock = bapl(ii)
         if (iblock==0) cycle
         idx     = (iblock-1)*WARP_SIZE + ilane;
         iipole  = pglob(idx)
         iglob   = ipole(idx)
         posi%x  =     x(idx)
         posi%y  =     y(idx)
         posi%z  =     z(idx)
         chgi    = chgct(iipole)
         alphai  = dmpct(iipole)
         !  Load atom block k parameters
         kdx     = abpl( (ii-1)*WARP_SIZE+ ilane )
         kpole   = pglob(kdx)
         kglob   = ipole(kdx)
         posk(threadIdx%x)%x  = x(kdx)
         posk(threadIdx%x)%y  = y(kdx)
         posk(threadIdx%x)%z  = z(kdx)
         chgk    = chgct(kpole)
         alphak  = dmpct(kpole)
           !* set compute Data to 0
#if __tver__ & __use_ene__
         frc_i%x = 0;
         frc_i%y = 0;
         frc_i%z = 0;
         frc_k%x = 0;
         frc_k%y = 0;
         frc_k%z = 0;
#endif
#if __tver__ & __use_ene__
          ect_ = 0
#endif
#if __tver__ & __use_act__
         nect_ = 0
#endif
#if __tver__ & __use_ene__
         vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
         vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;
#endif
         if (alphai.eq.0.0) alphai=1000.0
         if (alphak.eq.0.0) alphak=1000.0

         same_block = (idx.ne.kdx)

         do j = 0,warpsize-1
            srclane  = iand( ilane+j-1,warpsize-1 ) + 1
            klane    = threadIdx%x-ilane + srclane
#if __tfea__ & __use_mpi__
            if (nproc.gt.1) then
               xk_   = posk(klane)%x
               yk_   = posk(klane)%y
               zk_   = posk(klane)%z
               pos%x = posi%x - xk_
               pos%y = posi%y - yk_
               pos%z = posi%z - zk_
               call midpointimage_inl(xk_,yk_,zk_,pos%x,pos%y,pos%z)
               if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &         .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &         .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                  accept_mid = .false.
               else
                  accept_mid = .true.
                  pos%x=-pos%x; pos%y=-pos%y; pos%z=-pos%z;
               end if
            else
#endif
               pos%x = posk(klane)%x - posi%x
               pos%y = posk(klane)%y - posi%y
               pos%z = posk(klane)%z - posi%z
               call image_inl(pos%x,pos%y,pos%z)
#if __tfea__ & __use_mpi__
            end if
#endif
            d2      = pos%x**2 + pos%y**2 + pos%z**2
            do_pair = merge(.true.,iglob.lt.kglob,same_block)
#if __tfea__ & __use_longRange__
            if (do_pair.and.loff2<=d2.and.d2<=off2.and.accept_mid)
#else
            if (do_pair.and.d2<=off2.and.accept_mid)
#endif
     &           then
#if __tfea__ & __use_groups__
                 if (ugrp)
     &              call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
#endif
#if __tfea__ & __use_lambdadyn__
                 if (ulamdyn) mutik = muti + mutk(klane)
#endif
                 ! compute one interaction
                 call duo_chgtrn(d2,pos%x,pos%y,pos%z,chgi,chgk,alphai
     &                   ,alphak,f,1.0,fgrp,rinv,off,ctrnscut,cut2,sheal
     &                   ,ctrntyp,ugrp,e,frc,ver,fea)
#if __tver__ & __use_act__
                 if (e.ne.0.0) nect_ = nect_ + 1
#endif
#if __tver__ & __use_ene__
                 ect_  = ect_ + tp2enr(e)
#endif
#if __tver__ & __use_grd__
                 frc_i%x = frc_i%x - tp2mdr(frc%x)
                 frc_i%y = frc_i%y - tp2mdr(frc%y)
                 frc_i%z = frc_i%z - tp2mdr(frc%z)
                 frc_k%x = frc_k%x + tp2mdr(frc%x)
                 frc_k%y = frc_k%y + tp2mdr(frc%y)
                 frc_k%z = frc_k%z + tp2mdr(frc%z)
#endif
#if __tver__ & __use_vir__
                 vir_(1) = vir_(1) + pos%x* frc%x
                 vir_(2) = vir_(2) + pos%y* frc%x
                 vir_(3) = vir_(3) + pos%z* frc%x
                 vir_(4) = vir_(4) + pos%y* frc%y
                 vir_(5) = vir_(5) + pos%z* frc%y
                 vir_(6) = vir_(6) + pos%z* frc%z
#endif
              end if

                lot   = iand(ilane,warpsize-1)+1
              kglob   = __shfl(kglob  , lot)
              chgk    = __shfl(chgk   , lot)
              alphak  = __shfl(alphak , lot)
#if __tver__ & __use_grd__
              frc_k%x = __shfl(frc_k%x, lot)
              frc_k%y = __shfl(frc_k%y, lot)
              frc_k%z = __shfl(frc_k%z, lot)
#endif
           end do

           lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
#if __tver__ & __use_ene__
           ! Update energy buffer
           call atomic_add_f(e_buff(lot), ect_)
#endif
#if __tver__ & __use_act__
           ! Update counter buffer
           call atomic_add_i(n_buff(lot), nect_)
#endif
#if __tver__ & __use_grd__
           i       = loc  (idx)
           kbis    = loc  (kdx)
           ! Update forces
           call atomic_add_f( dect(1,i   ), frc_i%x )
           call atomic_add_f( dect(2,i   ), frc_i%y )
           call atomic_add_f( dect(3,i   ), frc_i%z )
           call atomic_add_f( dect(1,kbis), frc_k%x )
           call atomic_add_f( dect(2,kbis), frc_k%y )
           call atomic_add_f( dect(3,kbis), frc_k%z )
#endif
#if __tver__ & __use_vir__
           ! Update virial buffer
           call atomic_add_c(v_buff(0*RED_BUFF_SIZE+lot),vir_(1))
           call atomic_add_c(v_buff(1*RED_BUFF_SIZE+lot),vir_(2))
           call atomic_add_c(v_buff(2*RED_BUFF_SIZE+lot),vir_(3))
           call atomic_add_c(v_buff(3*RED_BUFF_SIZE+lot),vir_(4))
           call atomic_add_c(v_buff(4*RED_BUFF_SIZE+lot),vir_(5))
           call atomic_add_c(v_buff(5*RED_BUFF_SIZE+lot),vir_(6))
#endif
        end block
        end do

      end subroutine
