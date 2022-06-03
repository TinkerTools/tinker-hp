c
c     Hal CUDA Kernel Template 
c
c     "ehal1" calculates the buffered 14-7 van der waals energy and
c     its first derivatives with respect to cartesian coordinates
c     version __tver__ features __tfea__
c
      attributes(global)
     &subroutine __Cat( hal,__sufx__ )
     &         (xred,yred,zred,sgl_id,slc_id,loc_ired
     &         ,vblst1,ivblst,vblst,b_stat,b_rmid,grplist
     &         ,jvdw,epsilon,radmin,radv,epsv,wgrp,ired,kred
     &         ,devx,devy,devz,ev_buff,vir_buff,nred_buff,lam_buff
     &         ,nb2p,nbap,n,nbloc,nvdwlocnl
     &         ,nab,nvdwclass,radrule,epsrule
     &         ,c0,c1,c2,c3,c4,c5,rinv,shortheal,ghal,dhal
     &         ,cut2,cut,lowOff2,lowOff,off2,off
     &         ,scexp,vlambda,scalpha,dt1lam,dt2lam,mut
     &         ,ulamdyn,ugrp
     &         ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,rank
     &         ! Scaling Factor Params
     &         ,n_vscale,sizvc,vcorrect_ik,vcorrect_scale,loc
     &         ,jvdw_,xred_,yred_,zred_
     &         ,radmin4,epsilon4,radmin_,epsilon_,dev
     &         )

      implicit none
      integer  ,value:: nb2p,nbap,rank,n,nbloc,nvdwlocnl,nab
     &         ,nvdwclass,radrule,epsrule,n_vscale,sizvc
      logical  ,value:: ulamdyn,ugrp
      real(t_p),value:: c0,c1,c2,c3,c4,c5,rinv,shortheal
     &         ,ghal,dhal,cut2,cut,lowOff2,lowOff,off2,off
     &         ,scexp,vlambda,scalpha,dt1lam,dt2lam
     &         ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      integer(1),device :: mut(n)
      integer  ,device :: nred_buff (RED_BUFF_SIZE)
      ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
      real(r_p),device :: lam_buff(RED_BUFF_SIZE)
      real(t_p),device :: vir_buff(RED_BUFF_SIZE)
      integer  ,device,intent(in)::sgl_id(nab)
     &         ,loc_ired(nab),ivblst(nbap)
     &         ,vblst(nbap*(BLOCK_SIZE))
     &         ,vblst1(nbap),b_stat(*)
     &         ,slc_id(nab),grplist(n),ired(n),jvdw(nab)
      real(t_p),device,intent(in):: radmin(nvdwclass,nvdwclass),radv(*)
     &         ,epsilon(nvdwclass,nvdwclass),epsv(*),kred(n),b_rmid(4,*)
     &         ,xred(nab),yred(nab),zred(nab),wgrp(maxgrp+1,*)
      mdyn_rtyp,device,intent(inout)::devx(nbloc),devy(nbloc)
     &         ,devz(nbloc)

      integer  ,device,intent(in):: vcorrect_ik(sizvc,2)
     &         ,jvdw_(n),loc(n)
      real(t_p),device,intent(in):: vcorrect_scale(sizvc)
     &         ,xred_(nbloc),yred_(nbloc),zred_(nbloc)
     &         ,radmin_(maxclass,*),radmin4(maxclass,*)
     &         ,epsilon_(maxclass,*),epsilon4(maxclass,*)
      mdyn_rtyp,device::dev(3,nbloc)

      integer   ithread,iwarp,nwarp,ilane,iblock,kb
     &         ,idx,kdx,ii,j,iglob,kglob,ai12(maxvalue)
      integer(1) mutik,info
      logical   do_pair,accept_mid
      real(t_p) e,one
#if __tver__ & __use_ene__
      ener_rtyp ev_
#endif
      real(t_p) xi,yi,zi,xpos,ypos,zpos,xk,yk,zk
     &         ,rik2,rv2,eps2
     &         ,fgrp
#if __tver__ & __use_vir__
     &         ,vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
#endif
#if __tver__ & __use_grd__
      mdyn_rtyp gxk,gyk,gzk
      mdyn_rtyp gxi,gyi,gzi
#endif
#if __tfea__ & __use_softcore__
      integer(1) muti,mutk
#else
      parameter(mutik=0)
#endif
      parameter(one=1.0)
#if !(__tfea__ & __use_mpi__)
      parameter(accept_mid=.true.)
#endif

      ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
      iwarp   = (ithread-1) / warpsize
      nwarp   = blockDim%x*gridDim%x / warpsize
      ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
#if (__tfea__ & __use_mpi__)
      accept_mid = .true.
#endif
c                ==========================================
c     ================================================================
c                    Compute scaling factor Correction               |
c     ================================================================
c                ==========================================
      block
      integer it,kt
      logical do_scale4
      real(t_p) vscale,vscale4
      real(t_p) dedx,dedy,dedz,delambdav_
      do ii = ithread,n_vscale,blockDim%x*gridDim%x
         iglob  = vcorrect_ik(ii,1)
         kglob  = vcorrect_ik(ii,2)
         vscale = vcorrect_scale(ii)
         it     = jvdw_(iglob)
         kt     = jvdw_(kglob)
#if __tfea__ & __use_softcore__
         mutik  = mut(iglob) + mut(kglob)
#endif
         associate(i=>iglob,kbis=>kglob)
         i      = loc(iglob)
         kbis   = loc(kglob)
         do_scale4 = .false.
         vscale4   = 0

         if (vscale.lt.0) then 
            vscale4 = -vscale
            vscale  = one
         end if
c
c        compute the energy contribution for this interaction
c
         xpos   = xred_(i) - xred_(kbis)
         ypos   = yred_(i) - yred_(kbis)
         zpos   = zred_(i) - zred_(kbis)
         call image_inl(xpos,ypos,zpos)
c
c        decide whether to compute the current interaction
c        and check for an interaction distance less than the cutoff
c
         rik2   = xpos**2 + ypos**2 + zpos**2
#if __tfea__ & __use_longRange__
         if (rik2<lowOff2.or.rik2>off2) cycle
#else
         if (rik2>off2) cycle
#endif
#if __tfea__ & __use_softcore__
         if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilate
#endif
#if __tfea__ & __use_groups__
         if (ugrp) call groups2_inl(fgrp,iglob,kglob,grplist,wgrp)
#endif
 20      continue
         if (do_scale4) then
           !replace 1-4 interactions
            rv2  = radmin4 (kt,it)
            eps2 = epsilon4(kt,it)
         else
            rv2  =  radmin_(kt,it)
            eps2 = epsilon_(kt,it)
         end if

         call duo_hal
     &       (xpos,ypos,zpos,rik2,rv2,eps2,vscale,fgrp
     &       ,rinv,cut2,lowOff,off,shortheal,ghal,dhal,scexp
     &       ,vlambda,scalpha,dt1lam,dt2lam,mutik,ulamdyn,ugrp
     &       ,delambdav_,e,dedx,dedy,dedz
     &       ,__tver__+__use_sca__,__tfea__)

         if (.not.do_scale4) then
#if __tver__ & __use_ene__
            e   = -e
#endif
#if __tver__ & __use_act__
            if (vscale.eq.one) call atomic_sub_i( nred_buff(j),1 )
#endif
#if __tver__ & __use_grd__
            dedx = -dedx; dedy = -dedy; dedz = -dedz;
#endif
         end if

         j = iand(ithread-1,RED_BUFF_SIZE-1)+1
#if __tver__ & __use_ene__
         !increment the total van der Waals energy 
         call atomic_add_f( ev_buff(j),tp2enr(e) )
#endif
#if __tfea__ & __use_lambdadyn__
         if (ulamdyn) then
            delambdav_ = merge( delambdav_,-delambdav_,do_scale4 )
            call atomic_add_m( lam_buff(j),delambdav_ )
         end if
#endif
#if __tver__ & __use_grd__
         call atomic_sub_f1( dev(1,kbis),dedx )
         call atomic_sub_f1( dev(2,kbis),dedy )
         call atomic_sub_f1( dev(3,kbis),dedz )

         call atomic_add_f1( dev(1,i),dedx )
         call atomic_add_f1( dev(2,i),dedy )
         call atomic_add_f1( dev(3,i),dedz )
#endif
#if __tver__ & __use_vir__
         ! Increment virial term of van der Waals
         if (use_virial) then
            call atomic_add_c( vir_buff(0*RED_BUFF_SIZE+j),xpos*dedx )
            call atomic_add_c( vir_buff(1*RED_BUFF_SIZE+j),ypos*dedx )
            call atomic_add_c( vir_buff(2*RED_BUFF_SIZE+j),zpos*dedx )
            call atomic_add_c( vir_buff(3*RED_BUFF_SIZE+j),ypos*dedy )
            call atomic_add_c( vir_buff(4*RED_BUFF_SIZE+j),zpos*dedy )
            call atomic_add_c( vir_buff(5*RED_BUFF_SIZE+j),zpos*dedz )
         end if
#endif
         end associate
c
         ! deal with 1-4 Interactions
         if (vscale4.gt.0) then
            vscale    =  vscale4
            do_scale4 = .true.
            vscale4   = 0
            goto 20
         end if
      end do
      end block
c                  ==========================================
c       ================================================================
c                        Compute Pairwise interactions                 |
c       ================================================================
c                  ==========================================
      do ii = iwarp, nb2p+nbap-1, nwarp
      block
#if __tver__ & __use_act__
      integer nev_
#endif
#if __tfea__ & __radepsOpt__
      real(t_p) rvi,rvk,epsi,epsk
#else
      integer it,kt
#endif
      real(t_p) dedx,dedy,dedz,delambdav_
#if __tfea__ & __use_lambdadyn__
      real(r_p) delambdav
#endif
         if (ii.lt.nb2p) then
            kb  = vblst1(ii*2+2)
         kdx    = (kb-1)*warpsize + ilane
         else
         kdx    = vblst( (ii-nb2p)*warpsize+ilane )
         if(kdx.eq.0) kdx = nab
         end if
         ! Load atom block k neighbor parameter
         kglob  = sgl_id(kdx)
#if __tfea__ & __radepsOpt__
         rvk    = radv(kdx)
         epsk   = epsv(kdx)
#else
         kt     = jvdw(kdx)
#endif
         xk     = xred(kdx)
         yk     = yred(kdx)
         zk     = zred(kdx)
#if __tfea__ & __use_softcore__
         mutk   = mut(kglob)
#endif
        ! Set Data to compute to zero
#if __tver__ & __use_ene__
        ev_    = 0
#endif
#if __tver__ & __use_act__
        nev_    = 0
#endif
#if __tver__ & __use_grd__
        gxi    = 0; gyi    = 0; gzi    = 0;
        gxk    = 0; gyk    = 0; gzk    = 0;
#endif
#if __tver__ & __use_vir__
        vxx_   = 0; vxy_   = 0; vxz_   = 0;
        vyy_   = 0; vyz_   = 0; vzz_   = 0;
#endif
#if __tfea__ & __use_lambdadyn__
        if (ulamdyn) delambdav = 0
#endif

        ! Load atom block i parameters
        if (ii.lt.nb2p) then
           iblock = vblst1(ii*2+1)
        else
           iblock = ivblst(ii-nb2p+1)
        end if
        if (iblock.eq.0) cycle
        idx    = (iblock-1)*warpsize + ilane 
        iglob  = sgl_id(idx)
#if __tfea__ & __radepsOpt__
        rvi    = radv(idx)
        epsi   = epsv(idx)
#else
        it     = jvdw(idx)
#endif
        xi     = xred(idx)
        yi     = yred(idx)
        zi     = zred(idx)
#if __tfea__ & __use_softcore__
        muti   = mut(iglob)
#endif
        ! Study blocks proximity
        j      = b_stat(iblock)
        info   = merge(0,1,j.eq.i_bl)
        if (info.and.ii.lt.nb2p.and.j.eq.n_bl) then
           j = b_stat(kb)
           info = merge(0,1,j.eq.i_bl)
           if (j.eq.n_bl) then
              xpos = b_rmid(1,iblock) - b_rmid(1,kb)
              ypos = b_rmid(2,iblock) - b_rmid(2,kb)
              zpos = b_rmid(3,iblock) - b_rmid(3,kb)
              info = 0
              !if(xpos**2+ypos**2+zpos**2.gt.off2) info=1
              if( f_abs(xpos).gt.xcell2 ) info=1
              if( f_abs(ypos).gt.ycell2 ) info=1
              if( f_abs(zpos).gt.zcell2 ) info=1
           end if
        end if

        if (skipvdw12) then; block; integer k
          do k = 1,maxvalue
             ai12(k) = i12_p(k,iglob)
          end do
        end block; end if

        ! Interact block i with block k
        do j = 0,warpsize-1
           !srclane = iand( ilane+j-1,warpsize-1 ) + 1
           !klane   = threadIdx%x-ilane + srclane
           if (skipvdw12) then; block;
              ! Look for 1-2 bond interaction
              ! and skip it
              logical ik12
              integer k
              ik12 =.false.
              do k = 1, maxvalue
                 if (ai12(k).eq.kglob) ik12=.true.
              end do
              if (ik12) goto 10
           end block; end if
           if (info) then
#if __tfea__ & __use_mpi__
              if (nproc.gt.1) then
                 block; real(t_p) xk_,yk_,zk_
                 xk_   = xk
                 yk_   = yk
                 zk_   = zk
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
                 end block
              else
#endif
                 xpos    = xi - xk
                 ypos    = yi - yk
                 zpos    = zi - zk
                 if(info) call image_inl(xpos,ypos,zpos)
#if __tfea__ & __use_mpi__
              end if
#endif
           else
#if __tfea__ & __use_mpi__
              if (nproc.gt.1) then
                 block; real(t_p) xk_,yk_,zk_
                 xk_   = xk
                 yk_   = yk
                 zk_   = zk
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
                 end block
              else
#endif
                 xpos    = xi - xk
                 ypos    = yi - yk
                 zpos    = zi - zk
#if __tfea__ & __use_mpi__
              end if
#endif
           end if
#if __tver__ & __use_grd__
           dedx = 0.0; dedy = 0.0; dedz = 0.0;
#endif
           rik2 = xpos**2 + ypos**2 + zpos**2

           do_pair = merge(iglob.lt.kglob,.true.,idx.eq.kdx)
     &                .and.accept_mid
#if __tfea__ & __use_longRange__
           do_pair = do_pair.and.(lowOff2<rik2.and.rik2<=off2)
#else
           do_pair = do_pair.and.rik2<=off2
#endif
           if (do_pair) then
#if __tfea__ & __use_softcore__
              mutik   = muti + mutk
              if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilation
#endif
#if __tfea__ & __use_groups__
              if (ugrp) call groups2_inl(fgrp,iglob,kglob,grplist,wgrp)
#endif
#if __tfea__ & __radepsOpt__ 
              call get_rad(rvi,rvk,rv2,radrule)
              call get_eps(epsi,epsk,eps2,epsrule)
#else
              rv2   =  radmin_t (kt,it)
              eps2  = epsilon_t (kt,it)
#endif
              call duo_hal
     &            (xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p,fgrp
     &            ,rinv,cut2,lowOff,off,shortheal,ghal,dhal,scexp
     &            ,vlambda,scalpha,dt1lam,dt2lam,mutik,ulamdyn,ugrp
     &            ,delambdav_,e,dedx,dedy,dedz,__tver__,__tfea__)
#if __tver__ & __use_ene__
              ev_   = ev_ + tp2enr(e)
#endif
#if __tver__ & __use_act__
              nev_  = nev_ + 1
#endif
#if __tfea__ & __use_lambdadyn__
              if (ulamdyn) delambdav = delambdav + delambdav_
#endif
#if __tver__ & __use_vir__
              if (use_virial) then
              vxx_  = vxx_  + xpos * dedx
              vxy_  = vxy_  + ypos * dedx
              vxz_  = vxz_  + zpos * dedx
              vyy_  = vyy_  + ypos * dedy
              vyz_  = vyz_  + zpos * dedy
              vzz_  = vzz_  + zpos * dedz
              end if
#endif
              !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1
#if __tver__ & __use_grd__
              ! Accumulate interaction gradient
              gxi = gxi + tp2mdr(dedx)
              gyi = gyi + tp2mdr(dedy)
              gzi = gzi + tp2mdr(dedz)

              gxk = gxk + tp2mdr(dedx)
              gyk = gyk + tp2mdr(dedy)
              gzk = gzk + tp2mdr(dedz)
#endif
           end if

 10        continue
           block; integer il
           il    = iand(ilane,warpsize-1)+1
           kglob = __shfl( kglob,il )
#if __tfea__ & __use_softcore__
           mutk  = __shfl( int(mutk) ,il )
#endif
#if __tfea__ & __radepsOpt__ 
           rvk   = __shfl(  rvk,il )
           epsk  = __shfl( epsk,il )
#else
           kt    = __shfl( kt   ,il )
#endif
           xk    = __shfl( xk   ,il )
           yk    = __shfl( yk   ,il )
           zk    = __shfl( zk   ,il )
#if __tver__ & __use_grd__
           gxk   = __shfl( gxk  ,il )
           gyk   = __shfl( gyk  ,il )
           gzk   = __shfl( gzk  ,il )
#endif
           end block

        end do

        j = iand(ithread-1,RED_BUFF_SIZE-1) + 1
#if __tver__ & __use_ene__
        !increment the van der Waals energy
        call atomic_add_f( ev_buff(j),ev_ )
#endif
#if __tfea__ & __use_lambdadyn__
        if (ulamdyn)
     &     call atomic_add_m1(lam_buff(j) ,delambdav )
#endif
#if __tver__ & __use_act__
        !increment the van der Waals energy
        call atomic_add_i( nred_buff(j),nev_ )
#endif
#if __tver__ & __use_grd__
        !associate(kbis=>kglob,i=>iglob)
        !kbis   = slc_id(kdx)
        !i      = slc_id (idx)
        !increment the van der Waals derivatives
        if (idx.le.nvdwlocnl) then
           call atomic_add_f( devx(idx),gxi )
           call atomic_add_f( devy(idx),gyi )
           call atomic_add_f( devz(idx),gzi )
        end if
        if (kdx.le.nvdwlocnl) then
           call atomic_sub_f( devx(kdx),gxk )
           call atomic_sub_f( devy(kdx),gyk )
           call atomic_sub_f( devz(kdx),gzk )
        end if
        !end associate
#endif
#if __tver__ & __use_vir__
        ! Increment virial term of van der Waals
        if (use_virial) then
        call atomic_add_c( vir_buff(0*RED_BUFF_SIZE+j),vxx_ )
        call atomic_add_c( vir_buff(1*RED_BUFF_SIZE+j),vxy_ )
        call atomic_add_c( vir_buff(2*RED_BUFF_SIZE+j),vxz_ )
        call atomic_add_c( vir_buff(3*RED_BUFF_SIZE+j),vyy_ )
        call atomic_add_c( vir_buff(4*RED_BUFF_SIZE+j),vyz_ )
        call atomic_add_c( vir_buff(5*RED_BUFF_SIZE+j),vzz_ )
        end if
#endif
      end block
      end do
      end subroutine

