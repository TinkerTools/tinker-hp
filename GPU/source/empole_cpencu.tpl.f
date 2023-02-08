c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1cu" : CUDA Template for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates on device
c     version __tver__ features __tfea__
c
      attributes(global) subroutine __Cat(emreal_cpen,__sufx__)
     &      ( ipole, pglob, loc, bapl, abpl, nab, nabp
     &      , nbloc, n, pentyp
     &      , x, y, z, pcore, pval, palpha, rpole
     &      , shortheal, scut, loff2, off2, f, aewald
     &      , ugrp, grplist, wgrp, ulamdyn, elambda, mut, ucflx
     &      , dem, tem, pot, em_buff, vir_buff, lam_buff, nem_buff
     &      , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &      , mcorrect_ik,mcorrect_scale,ipole_,loc_,x_,y_,z_
     &      , n_mscale
     &      )
      implicit none
      integer  ,value :: nab,nbloc,n
     &         ,nabp,n_mscale,pentyp
      logical  ,value :: ugrp,ulamdyn,ucflx
      real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &         ,p_zbeg,p_zend,loff2,shortheal,elambda
     &         ,scut,off2,aewald,f
      integer  ,device,intent(in)::ipole(*),pglob(*),loc(*)
     &         ,bapl(*),abpl(*),grplist(*),ipole_(*),loc_(*)
     &         ,mcorrect_ik(2,*)
      integer(1),device,intent(in):: mut(n)
      real(t_p),device,intent(in):: x(*),y(*),z(*),x_(*),y_(*),z_(*)
     &         ,pcore(*),pval(*),palpha(*),rpole(13,*)
     &         ,wgrp(maxgrp+1,*),mcorrect_scale(*)
      real(t_p),device:: tem(3,*),vir_buff(*),pot(*)
      real(r_p),device:: lam_buff(*)
      integer  ,device:: nem_buff(*)
      mdyn_rtyp,device:: dem(3,*)
      ener_rtyp,device:: em_buff(*)

      integer     ithread,iwarp,nwarp,ilane,klane,srclane,lot
      integer     ii,j,iblock,idx,kdx,ver,fea
      integer     iipole,iglob,kpole,kglob,nem_
      integer(1)  muti,mutik
      integer(1) ,shared :: mutk(BLOCK_DIM)
      logical     do_pair,same_block,accept_mid,deb
      real(t_p)   xk_,yk_,zk_,d2,e
      real(t_p)   rstat,corei,corek,vali,valk,alphai,alphak,poti,potk
     &           ,fgrp,vir_(6)
#if __tfea__ & __use_lambdadyn__
      real(r_p)   delambdae
#endif
      ener_rtyp   em_
      type(real3)     pi,d,frc,ttmi
      type(real3)    ,shared:: pk(BLOCK_DIM),ttmk(BLOCK_DIM)
      type(mdyn3_r)   frc_i
      type(mdyn3_r)  ,shared:: frc_k(BLOCK_DIM)
      type(rpole_elt) ip
      type(rpole_elt),shared:: kp(BLOCK_DIM)

      parameter(ver=__tver__, fea=__tfea__
     &         ,deb=.false.)
c
      ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
      iwarp   = ((ithread-1) / warpsize) + 1
      nwarp   = blockDim%x*gridDim%x / warpsize
      ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
      accept_mid = .true.
      !deb = .false.
c
c     Correct Scaling Interactions
c
      do ii = ithread, n_mscale, blockDim%x*gridDim%x
      block;
      real(t_p) mscale,delambdae_
      integer kkpole
         iipole = mcorrect_ik(1,ii)
         kkpole = mcorrect_ik(2,ii)
         mscale =-mcorrect_scale(ii)
         iglob  = ipole_(iipole)
         kglob  = ipole_(kkpole)

         d%x  = x_(kglob) - x_(iglob)
         d%y  = y_(kglob) - y_(iglob)
         d%z  = z_(kglob) - z_(iglob)
         call image_inl(d%x,d%y,d%z)
         d2      = d%x**2 + d%y**2 + d%z**2
#if __tfea__ & __use_longRange__
         if (loff2.gt.d2.or.d2.gt.off2) cycle  !apply cutoff
#else
         if (d2.gt.off2) cycle  !apply cutoff
#endif
         corei  = pcore (iipole)
         vali   = pval  (iipole)
         alphai = palpha(iipole)
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
         corek  = pcore (kkpole)
         valk   = pval  (kkpole)
         alphak = palpha(kkpole)
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
         if (ucflx) then
            poti = 0
            potk = 0
         end if
#endif
#if __tver__ & __use_grd__
         call  real3_zr(ttmi)
         call  real3_zr(ttmk(threadIdx%x))
#endif
#if __tfea__ & __use_groups__
         if (ugrp) call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
#endif
         ! compute mpole one interaction
         call duo_mpole_cpen(d2,d%x,d%y,d%z,ip,kp(threadIdx%x)
     &           ,mscale,scut,shortheal,aewald,f
     &           ,pentyp,corei,corek,vali,valk,alphai,alphak
     &           ,ugrp,fgrp,ucflx
     &           ,e,frc,ttmi,ttmk(threadIdx%x)
     &           ,poti,potk,ver+__use_sca__,fea,deb)

         ! update energy
         lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
#if __tver__ & __use_ene__
         call atomic_add_f( em_buff(lot), tp2enr(e) )
#endif
#if __tver__ & __use_act__
         if (mscale.eq.-oner) call atomic_sub_i(nem_buff(lot), 1)
#endif
#if  __tfea__ & __use_lambdadyn__
         if (ulamdyn) call atomic_add_m(lam_buff(lot),delambdae_)
#endif
#if __tver__ & __use_grd__
         associate(iloc=>iglob,kloc=>kglob)
         iloc   = loc_(iglob)
         kloc   = loc_(kglob)
         ! increment force-based gradient and torque on first site
         call atomic_add_f1( dem(1,iloc), frc%x )
         call atomic_add_f1( dem(2,iloc), frc%y )
         call atomic_add_f1( dem(3,iloc), frc%z )
         call atomic_add_f1( dem(1,kloc),-frc%x )
         call atomic_add_f1( dem(2,kloc),-frc%y )
         call atomic_add_f1( dem(3,kloc),-frc%z )
         ! increment force-based gradient and torque on second site
         call atomic_add_c( tem(1,iloc), ttmi%x )
         call atomic_add_c( tem(2,iloc), ttmi%y )
         call atomic_add_c( tem(3,iloc), ttmi%z )
         call atomic_add_c( tem(1,kloc), ttmk(threadIdx%x)%x )
         call atomic_add_c( tem(2,kloc), ttmk(threadIdx%x)%y )
         call atomic_add_c( tem(3,kloc), ttmk(threadIdx%x)%z )
#  if __tfea__ & __use_chgflx__
         if (ucflx) then
            call atomic_add_c( pot(iloc),poti )
            call atomic_add_c( pot(kloc),potk )
         end if
#  endif
#  if __tver__ & __use_vir__
         ! increment the virial due to pairwise Cartesian forces c
         call atomic_add_c(vir_buff(0*RED_BUFF_SIZE+lot),-d%x*frc%x)
         call atomic_add_c(vir_buff(1*RED_BUFF_SIZE+lot)
     &              ,-0.5*(d%y*frc%x+d%x*frc%y))
         call atomic_add_c(vir_buff(2*RED_BUFF_SIZE+lot)
     &              ,-0.5*(d%z*frc%x+d%x*frc%z))
         call atomic_add_c(vir_buff(3*RED_BUFF_SIZE+lot),-d%y*frc%y)
         call atomic_add_c(vir_buff(4*RED_BUFF_SIZE+lot)
     &              ,-0.5*(d%z*frc%y+d%y*frc%z))
         call atomic_add_c(vir_buff(5*RED_BUFF_SIZE+lot),-d%z*frc%z)
#  endif
         end associate
#endif
        end block
        end do
c
c       Compute pairwise interaction
c
        do ii = iwarp, nabp, nwarp
        block; real(t_p) delambdae_
           iblock = bapl(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           pi%x  =     x(idx)
           pi%y  =     y(idx)
           pi%z  =     z(idx)
           corei  = pcore (iipole)
           vali   = pval  (iipole)
           alphai = palpha(iipole)
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
#if __tfea__ & __use_lambdadyn__
           muti    = mut(iglob)
#endif
           !  Load atom block k parameters
           kdx     = abpl( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           pk(threadIdx%x)%x  = x(kdx)
           pk(threadIdx%x)%y  = y(kdx)
           pk(threadIdx%x)%z  = z(kdx)
           corek   = pcore (kpole)
           valk    = pval  (kpole)
           alphak  = palpha(kpole)
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
#if __tfea__ & __use_lambdadyn__
           mutk(threadIdx%x) = mut(kglob)
#endif
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
#if __tver__ & __use_ene__
           em_ = 0
#endif
#if __tver__ & __use_act__
           nem_ = 0
#endif
#if __tfea__ & __use_lambdadyn__
           delambdae = 0
#endif
#if __tfea__ & __use_chgflx__
           if (ucflx) then
              poti=0; potk=0;
           end if
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
                 xk_   = pk(klane)%x
                 yk_   = pk(klane)%y
                 zk_   = pk(klane)%z
                 d%x = pi%x - xk_
                 d%y = pi%y - yk_
                 d%z = pi%z - zk_
                 call midpointimage_inl(xk_,yk_,zk_,d%x,d%y,d%z)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                    d%x=-d%x; d%y=-d%y; d%z=-d%z;
                 end if
              else
#endif
                 d%x = pk(klane)%x - pi%x
                 d%y = pk(klane)%y - pi%y
                 d%z = pk(klane)%z - pi%z
                 call image_inl(d%x,d%y,d%z)
#if __tfea__ & __use_mpi__
              end if
#endif
              d2      = d%x**2 + d%y**2 + d%z**2
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
              !if (iglob.Eq.1.and.kglob.EQ.2) then
              !   deb=.true.
              !else
              !   deb=.false.
              !end if
#if __tfea__ & __use_lambdadyn__
                 if (ulamdyn) mutik = muti + mutk(klane)
#endif
                 ! compute one interaction
                 call duo_mpole_cpen
     &                   (d2,d%x,d%y,d%z,ip,kp(klane),oner
     &                   ,scut,shortheal,aewald,f
     &                   ,pentyp,corei,corek,vali,valk,alphai,alphak
     &                   ,ugrp,fgrp,ucflx
     &                   ,e,frc,ttmi,ttmk(klane),poti,potk,ver,fea,deb)
#if __tfea__ & __use_lambdadyn__
                 if (ulamdyn) delambdae = delambdae + delambdae_
#endif
#if __tver__ & __use_ene__
                 em_  = em_ + tp2enr(e)
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
#  if __tver__ & __use_vir__
                 vir_(1) = vir_(1) - d%x*frc%x
                 vir_(2) = vir_(2) - 0.5*(d%y*frc%x + d%x*frc%y)
                 vir_(3) = vir_(3) - 0.5*(d%z*frc%x + d%x*frc%z)
                 vir_(4) = vir_(4) - d%y*frc%y
                 vir_(5) = vir_(5) - 0.5*(d%z*frc%y + d%y*frc%z)
                 vir_(6) = vir_(6) - d%z*frc%z
#  endif
#endif
              end if

                lot = iand(ilane,warpsize-1)+1
              kglob = __shfl( kglob,lot )
              corek = __shfl(  corek,lot )
               valk = __shfl(   valk,lot )
             alphak = __shfl( alphak,lot )
#if __tfea__ & __use_chgflx__
              potk  = __shfl( potk,lot ) 
#endif
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
           associate(iloc=>iglob,kloc=>kglob)
#if __tver__ & __use_grd__
           iloc    = loc  (idx)
           kloc    = loc  (kdx)
           ! Update forces
           call atomic_add_f( dem(1,iloc), frc_i%x )
           call atomic_add_f( dem(2,iloc), frc_i%y )
           call atomic_add_f( dem(3,iloc), frc_i%z )
           call atomic_add_f( dem(1,kloc), frc_k(threadIdx%x)%x )
           call atomic_add_f( dem(2,kloc), frc_k(threadIdx%x)%y )
           call atomic_add_f( dem(3,kloc), frc_k(threadIdx%x)%z )

           ! Update torque
           call atomic_add_c( tem(1,iloc),ttmi%x )
           call atomic_add_c( tem(2,iloc),ttmi%y )
           call atomic_add_c( tem(3,iloc),ttmi%z )
           call atomic_add_c( tem(1,kloc),ttmk(threadIdx%x)%x )
           call atomic_add_c( tem(2,kloc),ttmk(threadIdx%x)%y )
           call atomic_add_c( tem(3,kloc),ttmk(threadIdx%x)%z )
#endif
#if __tfea__ & __use_chgflx__
           if (ucflx) then
              call atomic_add_c( pot(iloc),poti )
              call atomic_add_c( pot(kloc),potk )
           end if
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
           end associate
        end block
        end do

        end subroutine
