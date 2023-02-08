c
c      version __tver__  features __tfea__
c
      attributes(global)
     &subroutine __Cat(lj,__sufx__)
     &          (x,y,z,sgl_id,slc_id
     &          ,b2pl,bapl,abpl,jvdw,i12,mut,epsilon,radmin,radv,epsv
     &          ,d_x,d_y,d_z,ev_buff,vir_buff,nred_buff
     &          ,nb2p,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nab
     &          ,nvdwclass,radrule,epsrule
     &          ,cut2,cut,off2,off,loff2,sheal,rinv
     &          ,ugrp,grplist,wgrp
     &          ! lambdaDyn
     &          ,sck,sct,scs,scalpha,galpha,ulamdyn,lambda,vlambda,glamb
     &          ,lambdavt,lambdavt1,lam_buff
     &          ! Scaling factor
     &          ,x_,y_,z_,correct_ik,correct_scale,szc,n_vscale
     &          ,loc,jvdw_,radmin1,epsilon1,radmin4,epsilon4,dev
     &          ! Box data
     &          ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &          )

      implicit none
      ! Input Scalar
      integer  ,value,intent(in) :: nb2p,nvdwlocnlb_pair,n,nbloc
     &         ,nvdwlocnl,nab,nvdwclass,radrule,epsrule
     &         ,szc,n_vscale
      logical  ,value:: ugrp,ulamdyn
      real(t_p),value,intent(in) :: cut2,cut,rinv,off2,off,loff2,sheal
     &         ,sck,sct,scs,scalpha,galpha,lambda,vlambda,glamb
     &         ,lambdavt,lambdavt1
     &         ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      ! Input Vector
      integer  ,device,intent(in):: sgl_id(nab)
     &         ,bapl(nvdwlocnlb_pair),b2pl(nb2p)
     &         ,abpl(nvdwlocnlb_pair*(BLOCK_SIZE))
     &         ,slc_id(nab),i12(maxvalue,n),jvdw(nab)
     &         ,correct_ik(szc,2),loc(n),jvdw_(*),grplist(*)
      integer(1),device:: mut(n)
      real(t_p),device,intent(in):: radmin(nvdwclass,nvdwclass)
     &         ,epsilon(nvdwclass,nvdwclass),radv(*),epsv(*)
     &         ,correct_scale(*), wgrp(ngrp+1,*)
     &         ,radmin1(maxclass,*),epsilon1(maxclass,*)
     &         ,radmin4(maxclass,*),epsilon4(maxclass,*)
     &         ,x(nab),y(nab),z(nab),x_(n),y_(n),z_(n)
      ! Output Array
      integer  ,device :: nred_buff(  RED_BUFF_SIZE)
      ener_rtyp,device ::   ev_buff(  RED_BUFF_SIZE)
      real(t_p),device ::  vir_buff(6*RED_BUFF_SIZE)
      real(r_p),device ::  lam_buff(  RED_BUFF_SIZE)
      mdyn_rtyp,device :: d_x(*),d_y(*),d_z(*),dev(3,*)

      ! Local workspace
#if __tfea__ & (__use_softcore__+__use_lambdadyn__)
      integer(1) muti,mutk,mutik
#endif
      integer ithread,iwarp,nwarp,ilane
      integer ii,j
      integer idx,kdx,iglob,kglob
      integer ai12(4)
      logical do_pair
      real(t_p) xi,yi,zi,xpos,ypos,zpos,xk,yk,zk,rik2,one,zero
      real(t_p) fgrp
#if __tver__ & __use_ene__
      ener_rtyp ev_
#endif
#if __tver__ & __use_grd__
      mdyn_rtyp gxi,gyi,gzi
      mdyn_rtyp gxk,gyk,gzk
#endif
#if __tfea__ & __use_lambdadyn__
      real(r_p) delambdav0
#endif
#if __tver__ & __use_vir__
      real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
#endif
      parameter( zero=0.0
     &         , one =1.0
#if !(__tfea__ & __use_groups__)
     &         , fgrp=1.0
#endif
     &         )

      ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
      iwarp   = (ithread-1) / warpsize
      !nwarp   = blockDim%x*gridDim%x / warpsize
      ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
      idx     = iand(ithread-1,RED_BUFF_SIZE-1) + 1
c
c     Scaling factor correction
c
      associate(do_scale4=>do_pair)
      do ii = ithread,n_vscale,blockDim%x*gridDim%x
      block
      integer  it,kt
      real(t_p) e,vscale,rvi,epsi,sc4,delambdav_
      type(real3) ded
         iglob  = correct_ik(ii,1)
         kglob  = correct_ik(ii,2)
         vscale = correct_scale(ii)
         it     = jvdw_(iglob)
         kt     = jvdw_(kglob)
         sc4    = 0
         do_scale4= .false.
         if (vscale.lt.zero) then ! 1-4 Interactions
            sc4   = -vscale
            vscale= one
         end if

         xpos   = x_(iglob) - x_(kglob)
         ypos   = y_(iglob) - y_(kglob)
         zpos   = z_(iglob) - z_(kglob)
         call image_inl(xpos,ypos,zpos)
c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
         rik2   = xpos**2 + ypos**2 + zpos**2
#if __tfea__ & __use_longRange__
         if (rik2>off2.or.rik2<loff2) cycle
#else
         if (rik2>off2) cycle
#endif
20       continue
        !replace 1-4 interactions
         if (do_scale4) then
            rvi  = radmin4 (kt ,it)
            epsi = epsilon4(kt ,it)
         else
           rvi  =  radmin1(kt ,it)
           epsi = epsilon1(kt ,it)
         end if
#if __tfea__ & __use_softcore__
         mutik = mut(iglob) + mut(kglob)
         if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilation
#endif
#if __tfea__ & __use_groups__
         if (ugrp)
     &      call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
#endif
c
c     compute the energy contribution for this interaction
c
         call duo_lj(rik2,xpos,ypos,zpos,rvi,epsi*vscale,cut2
     &              ,rinv,off,sheal,ugrp,fgrp,mutik
     &              ,sck,sct,scs,scalpha,galpha,ulamdyn
     &              ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &              ,delambdav_,e,ded,__tver__+__use_sca__,__tfea__)

#if __tver__ & __use_ene__
         e = merge(e,-e,do_scale4)
#endif
#if __tver__ & __use_grd__
         if (.not.do_scale4) then
            ded%x=-ded%x; ded%y=-ded%y; ded%z=-ded%z
         end if
#endif
#if __tver__ & __use_ene__
         block; ener_rtyp fstat;
         fstat= atomicadd( ev_buff(idx),tp2enr(e) )
         end block
#endif
#if __tfea__ & __use_lambdadyn__
         if (ulamdyn) then
            delambdav_ = merge( delambdav_,-delambdav_,do_scale4 )
            call atomic_add_m( lam_buff(idx),delambdav_ )
         end if
#endif
#if __tver__ & __use_act__
         if (vscale.eq.one) call atomic_sub_i( nred_buff(idx),1 )
#endif
         associate( iloc=>iglob,kloc=>kglob )
#if __tfea__ & __use_mpi__
         iloc   = loc(iglob)
         kloc   = loc(kglob)
#endif
#if __tver__ & __use_grd__
         call atomic_sub_f1( dev(1,kloc), ded%x )
         call atomic_sub_f1( dev(2,kloc), ded%y )
         call atomic_sub_f1( dev(3,kloc), ded%z )

         call atomic_add_f1( dev(1,iloc), ded%x )
         call atomic_add_f1( dev(2,iloc), ded%y )
         call atomic_add_f1( dev(3,iloc), ded%z )
#endif
         end associate
#if __tver__ & __use_vir__
         if (use_virial) then
        !increment the total van der Waals energy 
         call atomic_add_c( vir_buff(0*RED_BUFF_SIZE+idx),xpos*ded%x )
         call atomic_add_c( vir_buff(1*RED_BUFF_SIZE+idx),ypos*ded%x )
         call atomic_add_c( vir_buff(2*RED_BUFF_SIZE+idx),zpos*ded%x )
         call atomic_add_c( vir_buff(3*RED_BUFF_SIZE+idx),ypos*ded%y )
         call atomic_add_c( vir_buff(4*RED_BUFF_SIZE+idx),zpos*ded%y )
         call atomic_add_c( vir_buff(5*RED_BUFF_SIZE+idx),zpos*ded%z )
         end if
#endif
         ! deal with 1-4 Interactions
         if (sc4.gt.0) then
            vscale  =  sc4
            do_scale4 = .true.
            sc4     = 0
            goto 20
         end if
      end block
      end do
      end associate
c
c       Compute pairwise
c
      do ii = iwarp, nb2p+nvdwlocnlb_pair-1
     &      , blockDim%x*gridDim%x / warpsize
      block
#if __tver__ & __use_act__
      integer nev_
#endif
#if (__tfea__ & __radepsOpt__ )
      real(t_p) rvi,rvk,epsi,epsk
#else
      integer it,kt
#endif
         ! Load atom block k neighbor parameter
         if (ii.lt.nb2p) then
            kdx = b2pl(ii*2+2)
            kdx = (kdx-1)*warpsize + ilane
         else
            kdx = abpl( (ii-nb2p)*warpsize+ilane )
            if(kdx.eq.0) kdx = nab
         end if
         kglob  = sgl_id(kdx)
         !kbis   = slc_id(kdx)
         xk     = x     (kdx)
         yk     = y     (kdx)
         zk     = z     (kdx)
#if (__tfea__ & __radepsOpt__ )
         rvk    = radv  (kdx)
         epsk   = epsv  (kdx)
#else
         kt     = jvdw  (kdx)
#endif
#if __tfea__ & __use_softcore__
         mutk   = mut(kglob)
#endif
         ! Load atom block i parameters
         if (ii.lt.nb2p) then
            idx = b2pl(ii*2+1)
         else
            idx = bapl(ii-nb2p+1)
         end if
         if (idx.eq.0) cycle
         idx    = (idx-1)*warpsize + ilane 
         iglob  = sgl_id(idx)
         !i      = slc_id (idx)
         xi     = x     (idx)
         yi     = y     (idx)
         zi     = z     (idx)
#if (__tfea__ & __radepsOpt__ )
         rvi    = radv  (idx)
         epsi   = epsv  (idx)
#else
         it     = jvdw  (idx)
#endif
#if __tfea__ & __use_softcore__
         muti   = mut(iglob)
#endif
        ! Set Data to compute to zero
#if __tver__ & __use_ene__
         ev_    = 0
#endif
#if __tver__ & __use_act__
         nev_   = 0
#endif
#if __tver__ & __use_grd__
         gxi  = 0; gxk = 0
         gyi  = 0; gyk = 0
         gzi  = 0; gzk = 0
#endif
#if __tver__ & __use_vir__
         vxx_ = 0; vyy_ = 0
         vxy_ = 0; vyz_ = 0
         vxz_ = 0; vzz_ = 0
#endif
#if __tfea__ & __use_lambdadyn__
         if (ulamdyn) delambdav0 = 0
#endif
         if (skipvdw12) then; block; integer k;
            do k = 1,4; ai12(k) = i12(k,iglob); end do;
         end block; end if

         ! Interact block i with block k
         do j = 0,warpsize-1
            !block; integer srclane,klane
            !srclane = iand( ilane+j-1,warpsize-1 ) + 1
            !klane   = threadIdx%x-ilane + srclane
            !end block

            do_pair = merge(.true.,iglob.lt.kglob,idx.ne.kdx)

            if (skipvdw12) then; block; logical ik11,ik12,ik13
               ik11 = merge(.false.,.true.,ai12(1).eq.kglob)
               ik12 = merge(.false.,.true.,ai12(2).eq.kglob)
               ik13 = merge(.false.,.true.,ai12(3).eq.kglob)
               ik12 = merge(.false.,  ik12,ai12(4).eq.kglob)
               do_pair = do_pair.and.ik11.and.ik12.and.ik13
            end block; end if
#if __tfea__ & __use_mpi__
            if (ndir.gt.1) then
               block
               real(t_p) xk_,yk_,zk_
               xk_   = xk
               yk_   = yk
               zk_   = zk
               xpos  = xi - xk_
               ypos  = yi - yk_
               zpos  = zi - zk_
               call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
               if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &         .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &         .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                  do_pair = do_pair.and. .false.
               else
                  do_pair = do_pair.and. .true.
               end if
               end block
            else
#endif
               xpos    = xi - xk
               ypos    = yi - yk
               zpos    = zi - zk
               call image_inl(xpos,ypos,zpos)
#if __tfea__ & __use_mpi__
            end if
#endif
            rik2  = xpos**2 + ypos**2 + zpos**2
#if __tfea__ & __use_longRange__
            do_pair = do_pair.and.rik2<=off2.and.rik2>=loff2
#else
            do_pair = do_pair.and.rik2<=off2
#endif

            if (do_pair) then; block
            real(t_p) e,rv2,eps2,delambdav_
            type(real3) ded
#if (__tfea__ & __radepsOpt__ )
               call get_rad(rvi,rvk,rv2,radrule)
               call get_eps(epsi,epsk,eps2,epsrule)
#else
               rv2   =  radmin (kt,it)
               eps2  = epsilon (kt,it)
#endif
#if __tfea__ & __use_softcore__
               mutik = muti + mutk
               if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilation
#endif
#if __tfea__ & __use_groups__
               if (ugrp)
     &            call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
#endif
               call duo_lj(rik2,xpos,ypos,zpos,rv2,eps2,cut2
     &                    ,rinv,off,sheal,ugrp,fgrp,mutik
     &                    ,sck,sct,scs,scalpha,galpha,ulamdyn
     &                    ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                    ,delambdav_,e,ded,__tver__,__tfea__)
#if __tver__ & __use_ene__
               ev_   = ev_ + tp2enr(e)
#endif
#if __tver__ & __use_act__
               nev_  = nev_+1
#endif
#if __tfea__ & __use_lambdadyn__
               delambdav0 = delambdav0 + delambdav_
#endif
#if __tver__ & __use_grd__
               ! Accumulate interaction gradient
               gxi = gxi + tp2mdr(ded%x)
               gyi = gyi + tp2mdr(ded%y)
               gzi = gzi + tp2mdr(ded%z)

               gxk = gxk + tp2mdr(ded%x)
               gyk = gyk + tp2mdr(ded%y)
               gzk = gzk + tp2mdr(ded%z)
#endif
#if __tver__ & __use_vir__
               if (use_virial) then
               vxx_  = vxx_  + xpos * ded%x
               vxy_  = vxy_  + ypos * ded%x
               vxz_  = vxz_  + zpos * ded%x
               vyy_  = vyy_  + ypos * ded%y
               vyz_  = vyz_  + zpos * ded%y
               vzz_  = vzz_  + zpos * ded%z
               end if
#endif
            end block
            end if

c 10        continue
            block; integer il
            il   = iand(ilane,warpsize-1)+1
            kglob= __shfl(kglob,il )
#if (__tfea__ & __radepsOpt__ )
            rvk  = __shfl(  rvk,il )
            epsk = __shfl( epsk,il )
#else
            kt   = __shfl(   kt,il )
#endif
#if __tfea__ & __use_softcore__
            mutk = __shfl( int(mutk),il )
#endif
            xk   = __shfl(   xk,il )
            yk   = __shfl(   yk,il )
            zk   = __shfl(   zk,il )
#if __tver__ & __use_grd__
            gxk  = __shfl(  gxk,il )
            gyk  = __shfl(  gyk,il )
            gzk  = __shfl(  gzk,il )
#endif
            end block
         end do

         associate(itr=>iglob)
         itr   = iand(ithread-1,RED_BUFF_SIZE-1) + 1
#if __tver__ & __use_ene__
         !increment the van der Waals energy
         call atomic_add_f( ev_buff(itr) ,ev_ )
#endif
#if __tfea__ & __use_lambdadyn__
         if (ulamdyn)
     &      call atomic_add_m1(lam_buff(itr) ,delambdav0 )
#endif
#if __tver__ & __use_act__
         !increment interaction counter
         call atomic_add_i(nred_buff(itr),nev_)
#endif
#if __tver__ & __use_vir__
         ! Increment virial term of van der Waals
         if (use_virial) then
         call atomic_Add_c( vir_buff(0*RED_BUFF_SIZE+itr),vxx_ )
         call atomic_Add_c( vir_buff(1*RED_BUFF_SIZE+itr),vxy_ )
         call atomic_Add_c( vir_buff(2*RED_BUFF_SIZE+itr),vxz_ )
         call atomic_Add_c( vir_buff(3*RED_BUFF_SIZE+itr),vyy_ )
         call atomic_Add_c( vir_buff(4*RED_BUFF_SIZE+itr),vyz_ )
         call atomic_Add_c( vir_buff(5*RED_BUFF_SIZE+itr),vzz_ )
         end if
#endif
         end associate
#if __tver__ & __use_grd__
         !increment the van der Waals derivatives
         if (idx.le.nvdwlocnl) then
            call atomic_add_f (d_x(idx),gxi)
            call atomic_add_f (d_y(idx),gyi)
            call atomic_add_f (d_z(idx),gzi)
         end if
         if (kdx.le.nvdwlocnl) then
            call atomic_sub_f (d_x(kdx),gxk)
            call atomic_sub_f (d_y(kdx),gyk)
            call atomic_sub_f (d_z(kdx),gzk)
         end if
#endif
      end block
      end do
      end subroutine
