c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1cu" : CUDA Template for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates on device
c     version __tver__ features __tfea__
c
c
      attributes(global) subroutine __Cat(emreal,__sufx__)
     &      ( ipole, pglob, loc, ieblst, eblst
     &      , npolelocnlb, npolelocnlb_pair, npolebloc, n
     &      , x, y, z, rpole
     &      , shortheal, scut, loff2, off2, f
     &      , alsq2, alsq2n, aewald
     &      , ugrp, grplist, wgrp
     &      , ulamdyn, elambda, mut, u_cflx, pot
     &      , dem, tem, em_buff, vir_buff, lam_buff, nem_buff
     &      , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &      , mcorrect_ik,mcorrect_scale,ipole_,loc_,x_,y_,z_
     &      , n_mscale
     &      )
      implicit none
      integer  ,value :: npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair,n_mscale
      logical  ,value :: ugrp,ulamdyn,u_cflx
      real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &         ,p_zbeg,p_zend,loff2,shortheal,elambda
     &         ,scut,off2,alsq2,alsq2n,aewald,f
      integer  ,device,intent(in)::ipole(*),pglob(*),loc(*)
     &         ,ieblst(*),eblst(*),grplist(*),ipole_(*),loc_(*)
     &         ,mcorrect_ik(2,*)
      integer(1),device,intent(in):: mut(n)
      real(t_p),device,intent(in):: x(*),y(*),z(*),x_(*),y_(*),z_(*)
     &         , rpole(13,*),wgrp(ngrp+1,*),mcorrect_scale(*)
      real(t_p),device:: tem(3,*),pot(*),vir_buff(*)
      real(r_p),device:: lam_buff(*)
      integer  ,device:: nem_buff(*)
      mdyn_rtyp,device:: dem(3,*)
      ener_rtyp,device:: em_buff(*)

      integer ithread,iwarp,nwarp,ilane,klane,srclane,lot
      integer ii,j,i,kbis,iblock,idx,kdx,ver,fea
      integer iipole,iglob,iploc,kpole,kglob,kploc,nem_
      integer(1) muti,mutik
      integer(1),shared :: mutk(BLOCK_DIM)
      real(t_p) xk_,yk_,zk_,d2
      ener_rtyp em_
      real(r_p) delambdae
      real(t_p) rstat,fgrp
      type(real3) posi,pos
      type(real3) frc
      type(mdyn3_r) frc_i
      type(mdyn3_r),shared::frc_k(BLOCK_DIM)
      type(real3)  ,shared::posk(BLOCK_DIM)
      type(real3)   ttmi
      type(real3)  ,shared:: ttmk(BLOCK_DIM)
      real(t_p)    ,shared:: potk(BLOCK_DIM)
      real(t_p) vir_(6),poti
      type(rpole_elt) ip
      type(rpole_elt),shared:: kp(BLOCK_DIM)
      logical do_pair,same_block,accept_mid
      parameter(ver = __tver__, fea=__tfea__)
c
      ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
      iwarp   = ((ithread-1) / warpsize) + 1
      nwarp   = blockDim%x*gridDim%x / warpsize
      ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
      accept_mid = .true.
c
c     Correct Scaling Interactions
c
      do ii = ithread, n_mscale, blockDim%x*gridDim%x
      block;
      real(t_p) mscale,delambdae_
      integer kkpole
         iipole = mcorrect_ik(1,ii)
         kkpole = mcorrect_ik(2,ii)
         mscale = mcorrect_scale(ii)
         iglob  = ipole_(iipole)
         kglob  = ipole_(kkpole)
         i      = loc_(iglob)
         kbis   = loc_(kglob)

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
         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)
         kp(threadIdx%x)%c   = rpole(01,kkpole)
         kp(threadIdx%x)%dx  = rpole(02,kkpole)
         kp(threadIdx%x)%dy  = rpole(03,kkpole)
         kp(threadIdx%x)%dz  = rpole(04,kkpole)
         kp(threadIdx%x)%qxx = rpole(05,kkpole)
         kp(threadIdx%x)%qxy = rpole(06,kkpole)
         kp(threadIdx%x)%qxz = rpole(07,kkpole)
         kp(threadIdx%x)%qyy = rpole(09,kkpole)
         kp(threadIdx%x)%qyz = rpole(10,kkpole)
         kp(threadIdx%x)%qzz = rpole(13,kkpole)
#if __tfea__ & __use_lambdadyn__
         mutik  = mut(iglob)+mut(kglob)
#endif
#if __tfea__ & __use_chgflx__
         if (u_cflx) then
            poti = 0
            potk(threadIdx%x) = 0
         end if
#endif
c20        continue
#if __tver__ & __use_ene__
         em_ = 0
#endif
#if __tver__ & __use_grd__
         call  real3_zr(ttmi)
         call  real3_zr(ttmk(threadIdx%x))
#endif
#if __tfea__ & __use_groups__
         if (ugrp) call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
#endif
         ! compute mpole one interaction
         call duo_mpole(d2,pos%x,pos%y,pos%z,ip,kp(threadIdx%x),mscale
     &           ,scut,shortheal,aewald,f,alsq2n,alsq2,ugrp,fgrp
     &           ,ulamdyn,mutik,elambda,u_cflx,poti,potk(threadIdx%x)
     &           ,delambdae_,em_,frc,ttmi,ttmk(threadIdx%x)
     &           ,ver+__use_sca__,fea)

         ! update energy
         lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
#if __tver__ & __use_ene__
         call atomic_add_f( em_buff(lot), em_ )
#endif
#if __tver__ & __use_act__
         if (mscale.eq.oner) call atomic_sub_i(nem_buff(lot), 1)
#endif
#if  __tfea__ & __use_lambdadyn__
         if (ulamdyn) call atomic_add_m(lam_buff(lot),delambdae_)
#endif
#if __tfea__ & __use_chgflx__
         if (u_cflx) then
            call atomic_add_c(pot(i),poti)
            call atomic_add_c(pot(kbis),potk(threadIdx%x))
         end if
#endif
#if __tver__ & __use_vir__
         ! increment the virial due to pairwise Cartesian forces c
         call atomic_add_c(vir_buff(0*RED_BUFF_SIZE+lot),-pos%x*frc%x)
         call atomic_add_c(vir_buff(1*RED_BUFF_SIZE+lot),-pos%y*frc%x)
         call atomic_add_c(vir_buff(2*RED_BUFF_SIZE+lot),-pos%z*frc%x)
         call atomic_add_c(vir_buff(3*RED_BUFF_SIZE+lot),-pos%y*frc%y)
         call atomic_add_c(vir_buff(4*RED_BUFF_SIZE+lot),-pos%z*frc%y)
         call atomic_add_c(vir_buff(5*RED_BUFF_SIZE+lot),-pos%z*frc%z)
#endif
#if __tver__ & __use_grd__
         ! increment force-based gradient and torque on first site
         call atomic_add_f1( dem(1,i)    , frc%x )
         call atomic_add_f1( dem(2,i)    , frc%y )
         call atomic_add_f1( dem(3,i)    , frc%z )
         call atomic_add_f1( dem(1,kbis) ,-frc%x )
         call atomic_add_f1( dem(2,kbis) ,-frc%y )
         call atomic_add_f1( dem(3,kbis) ,-frc%z )
         ! increment force-based gradient and torque on second site
         call atomic_add_c( tem(1,i)    , ttmi%x )
         call atomic_add_c( tem(2,i)    , ttmi%y )
         call atomic_add_c( tem(3,i)    , ttmi%z )
         call atomic_add_c( tem(1,kbis) , ttmk(threadIdx%x)%x )
         call atomic_add_c( tem(2,kbis) , ttmk(threadIdx%x)%y )
         call atomic_add_c( tem(3,kbis) , ttmk(threadIdx%x)%z )
#endif
        end block
        end do
c
c       Compute pairwise interaction
c
        do ii = iwarp, npolelocnlb_pair, nwarp
        block; real(t_p) delambdae_
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ip%c    = rpole( 1, iipole)
           ip%dx   = rpole( 2, iipole)
           ip%dy   = rpole( 3, iipole)
           ip%dz   = rpole( 4, iipole)
           ip%qxx  = rpole( 5, iipole)
           ip%qxy  = rpole( 6, iipole)
           ip%qxz  = rpole( 7, iipole)
           ip%qyy  = rpole( 9, iipole)
           ip%qyz  = rpole(10, iipole)
           ip%qzz  = rpole(13, iipole)
           muti    = mut(iglob)
           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kp(threadIdx%x)%c    = rpole( 1, kpole)
           kp(threadIdx%x)%dx   = rpole( 2, kpole)
           kp(threadIdx%x)%dy   = rpole( 3, kpole)
           kp(threadIdx%x)%dz   = rpole( 4, kpole)
           kp(threadIdx%x)%qxx  = rpole( 5, kpole)
           kp(threadIdx%x)%qxy  = rpole( 6, kpole)
           kp(threadIdx%x)%qxz  = rpole( 7, kpole)
           kp(threadIdx%x)%qyy  = rpole( 9, kpole)
           kp(threadIdx%x)%qyz  = rpole(10, kpole)
           kp(threadIdx%x)%qzz  = rpole(13, kpole)
           mutk(threadIdx%x) = mut(kglob)
           !* set compute Data to 0
#if __tver__ & __use_ene__
           frc_i%x = 0;
           frc_i%y = 0;
           frc_i%z = 0;
           frc_k(threadIdx%x)%x = 0;
           frc_k(threadIdx%x)%y = 0;
           frc_k(threadIdx%x)%z = 0;
           ttmi%x  = 0.0;
           ttmi%y  = 0.0;
           ttmi%z  = 0.0;
           ttmk(threadIdx%x)%x = 0.0;
           ttmk(threadIdx%x)%y = 0.0;
           ttmk(threadIdx%x)%z = 0.0;
#endif
#if __tfea__ & __use_chgflx__
           if (u_cflx) then
              poti = 0
              potk(threadIdx%x) = 0
           end if
#endif
#if __tver__ & __use_ene__
           em_ = 0
#endif
#if __tver__ & __use_act__
           nem_ = 0
#endif
#if __tfea__ & __use_lambdadyn__
           delambdae = 0
#endif
#if __tver__ & __use_ene__
           vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
           vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;
#endif
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
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
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
                 mutik = muti + mutk(klane)
                 if (.not.(mutik.eq.1.and.elambda.eq.0)) then

                 ! compute one interaction
                 call duo_mpole(d2,pos%x,pos%y,pos%z,ip,kp(klane),zeror
     &                   ,scut,shortheal,aewald,f,alsq2n,alsq2,ugrp,fgrp
     &                   ,ulamdyn,mutik,elambda,u_cflx,poti,potk(klane)
     &                   ,delambdae_,em_,frc,ttmi,ttmk(klane)
     &                   ,ver,fea)
#if __tfea__ & __use_lambdadyn__
                 if (ulamdyn) delambdae = delambdae + delambdae_
#endif
#if __tver__ & __use_act__
                 nem_ = nem_ + 1
#endif
#if __tver__ & __use_grd__
                 frc_i%x = frc_i%x + tp2mdr(frc%x)
                 frc_i%y = frc_i%y + tp2mdr(frc%y)
                 frc_i%z = frc_i%z + tp2mdr(frc%z)
                 !store in large container for mixed precision
                 frc_k(klane)%x = frc_k(klane)%x - tp2mdr(frc%x)
                 frc_k(klane)%y = frc_k(klane)%y - tp2mdr(frc%y)
                 frc_k(klane)%z = frc_k(klane)%z - tp2mdr(frc%z)
#endif
#if __tver__ & __use_vir__
                 vir_(1) = vir_(1) - pos%x * frc%x
                 vir_(2) = vir_(2) - pos%y * frc%x
                 vir_(3) = vir_(3) - pos%z * frc%x
                 vir_(4) = vir_(4) - pos%y * frc%y
                 vir_(5) = vir_(5) - pos%z * frc%y
                 vir_(6) = vir_(6) - pos%z * frc%z
#endif
                 end if
                end if

                lot = iand(ilane,warpsize-1)+1
              kglob = __shfl(kglob, lot)
           end do

           lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
#if __tver__ & __use_ene__
           ! Update energy buffer
           call atomic_add_f(em_buff(lot), em_)
#endif
#if __tfea__ & __use_lambdadyn__
           call atomic_add_m1(lam_buff(lot), delambdae)
#endif
#if __tver__ & __use_act__
           ! Update counter buffer
           call atomic_add_i(nem_buff(lot), nem_)
#endif
#if __tfea__ & __use_chgflx__
           if (u_cflx) then
              call atomic_add_c(pot(i),poti)
              call atomic_add_c(pot(kbis),potk(threadIdx%x))
           end if
#endif
#if __tver__ & __use_grd__
           ! Update forces
           call atomic_add_f( dem(1,i   ), frc_i%x )
           call atomic_add_f( dem(2,i   ), frc_i%y )
           call atomic_add_f( dem(3,i   ), frc_i%z )
           call atomic_add_f( dem(1,kbis), frc_k(threadIdx%x)%x )
           call atomic_add_f( dem(2,kbis), frc_k(threadIdx%x)%y )
           call atomic_add_f( dem(3,kbis), frc_k(threadIdx%x)%z )

           ! Update torque
           call atomic_add_c( tem(1,i)   ,ttmi%x )
           call atomic_add_c( tem(2,i)   ,ttmi%y )
           call atomic_add_c( tem(3,i)   ,ttmi%z )
           call atomic_add_c( tem(1,kbis),ttmk(threadIdx%x)%x )
           call atomic_add_c( tem(2,kbis),ttmk(threadIdx%x)%y )
           call atomic_add_c( tem(3,kbis),ttmk(threadIdx%x)%z )
#endif
#if __tver__ & __use_vir__
           ! Update virial buffer
           call atomic_add_c(vir_buff(0*RED_BUFF_SIZE+lot),vir_(1))
           call atomic_add_c(vir_buff(1*RED_BUFF_SIZE+lot),vir_(2))
           call atomic_add_c(vir_buff(2*RED_BUFF_SIZE+lot),vir_(3))
           call atomic_add_c(vir_buff(3*RED_BUFF_SIZE+lot),vir_(4))
           call atomic_add_c(vir_buff(4*RED_BUFF_SIZE+lot),vir_(5))
           call atomic_add_c(vir_buff(5*RED_BUFF_SIZE+lot),vir_(6))
#endif
        end block
        end do

        end subroutine
