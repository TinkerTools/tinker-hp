c
c     Disp CUDA Kernel Template 
c
c     "edispcu" calculates the damped dispersion energy and first
c     derivatives with respect to Cartesian coordinates
c
c     version __tver__ features __tfea__
c
      attributes(global)
     &subroutine __Cat( disp,__sufx__ )
     &         (x,y,z,sgl_id,slc_id,i12
     &         ,b2pl,bapl,abpl,b_stat,b_rmid
     &         ,dedspx,dedspy,dedspz,dedsp,ev_buff,vir_buff,nred_buff
     &         ,nb2p,nbap,n,na,nab
     &         ,rinv,aewald,shortheal,cut2,cut,sCut2,sCut,off2,off
     &         ,mut,vlambda,grplist,wgrp,i_grp
     &         ! Box data 
     &         ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,rank
     &         ! Scaling Factor Params
     &         ,n_scale,correct_ik,correct_scale,loc,x_,y_,z_
     &         )

      implicit none
      integer  ,value:: nb2p,nbap,n,na,nab,i_grp,rank
      real(t_p),value:: c0,c1,c2,c3,c4,c5,rinv,aewald,shortheal
     &         ,cut2,cut,sCut2,sCut,off2,off,vlambda
      real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend
     &         ,p_zbeg,p_zend
      integer(1),device,intent(in):: mut(n)
      integer  ,device,intent(in)::sgl_id(nab),slc_id(nab)
     &         ,b2pl(nbap),bapl(nbap),abpl(nbap*(BLOCK_SIZE))
     &         ,b_stat(*),i12(maxvalue,*),grplist(n)
      real(t_p),device,intent(in):: b_rmid(4,*),x(nab)
     &         ,y(nab),z(nab),wgrp(maxgrp+1,*)
      integer  ,device :: nred_buff (RED_BUFF_SIZE)
      ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
      real(t_p),device :: vir_buff(RED_BUFF_SIZE)
      mdyn_rtyp,device :: dedspx(na),dedspy(na),dedspz(na)
     &                    dedsp(3,nbloc)

      integer  ,value:: n_scale
      integer  ,device,intent(in):: correct_ik(2,*),loc(n)
      real(t_p),device,intent(in):: correct_scale(*)
     &         ,x_(n),y_(n),z_(n)

      integer ithread,iwarp,nwarp,ilane,FEA,VER
      integer iblock,kb
      integer idx,kdx,ii,j,iglob
      integer kglob
      integer ai12(4)
      real(t_p) e,one
#if __tver__ & __use_ene__
      ener_rtyp ed_
#endif
      real(t_p) xi,yi,zi,xk,yk,zk,xr,yr,zr,rik2
      real(t_p) ai,ci,ak,ck,fgrp
      type(real3) ded
#if __tver__ & __use_grd__
      mdyn_rtyp gxk,gyk,gzk
      mdyn_rtyp gxi,gyi,gzi
#endif
#if __tver__ & __use_vir__
      real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
#endif
      logical do_pair,accept_mid
      integer(1) mutik,info
#if __tfea__ & __use_softcore__
      integer(1) muti,mutk
#else
      parameter(mutik=0)
#endif
      parameter(one=1.0,FEA=__tfea__,VER=__tver__)
#if !(__tfea__ & __use_mpi__)
      parameter(accept_mid=.true.)
#endif
#if __tfea__ & __use_groups__
      parameter(fgrp=1.0)
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
      integer ver0
      logical do_scale4
      real(t_p) scale,vscale4
      parameter(ver0=ver+__use_sca__)
      do ii = ithread,n_scale,blockDim%x*gridDim%x
         iglob  = correct_ik(ii,1)
         kglob  = correct_ik(ii,2)
         scale  = correct_scale(ii)
#if __tfea__ & __use_softcore__
         mutik  = mut(iglob) + mut(kglob)
#endif
         associate(i=>iglob,kbis=>kglob)
         i      = loc(iglob)
         kbis   = loc(kglob)
         ci     = csix(iglob)
         ck     = csix(kglob)
         ai     = adisp(iglob)
         ak     = adisp(kglob)
c
c        compute the energy contribution for this interaction
c
         xr     = x_(i) - x_(kbis)
         yr     = y_(i) - y_(kbis)
         zr     = z_(i) - z_(kbis)
         call image_inl(xr,yr,zr)
c
c        decide whether to compute the current interaction
c        and check for an interaction distance less than the cutoff
c
         rik2   = xr**2 + yr**2 + zr**2
#if __tfea__ & __use_longRange__
         if (rik2<sCut2.or.rik2>off2) cycle
#else
         if (rik2>off2) cycle
#endif
#if __tfea__ & __use_softcore__
         if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilate
#endif
#if __tfea__ & __use_groups__
         if (i_grp) call groups2_inl(fgrp,iglob,kglob,grplist,wgrp)
#endif

#if __tfea__ & __use_ewald__
         call duo_disp_real(rik2,xr,yr,zr,ai,ci,ak,ck
     &                  ,aewald,scale,mutik,vlambda,i_grp,fgrp
     &                  ,sCut,shortheal,e,ded,ver0,fea)
#else
         call duo_disp(rik2,xr,yr,zr,ai,ci,ak,ck,rinv,sCut,off
     &                ,shortheal,mutik,vlambda,i_grp,fgrp,scale
     &                ,e,ded,ver0,fea)
#endif

#if __tver__ & __use_act__
         if (scale.eq.one) call atomic_sub_i( nred_buff(j),1 )
#endif
         j = iand(ithread-1,RED_BUFF_SIZE-1)+1
#if __tver__ & __use_ene__
         !increment the total energy 
         call atomic_add_f( ev_buff(j),tp2enr(e) )
#endif
#if __tver__ & __use_grd__
         call atomic_sub_f1( dedsp(1,kbis),ded%x )
         call atomic_sub_f1( dedsp(2,kbis),ded%y )
         call atomic_sub_f1( dedsp(3,kbis),ded%z )

         call atomic_add_f1( dedsp(1,i),ded%x )
         call atomic_add_f1( dedsp(2,i),ded%y )
         call atomic_add_f1( dedsp(3,i),ded%z )
#endif
#if __tver__ & __use_vir__
         ! Increment virial term
         if (use_virial) then
            call atomic_add_c( vir_buff(0*RED_BUFF_SIZE+j),xr*ded%x )
            call atomic_add_c( vir_buff(1*RED_BUFF_SIZE+j),yr*ded%x )
            call atomic_add_c( vir_buff(2*RED_BUFF_SIZE+j),zr*ded%x )
            call atomic_add_c( vir_buff(3*RED_BUFF_SIZE+j),yr*ded%y )
            call atomic_add_c( vir_buff(4*RED_BUFF_SIZE+j),zr*ded%y )
            call atomic_add_c( vir_buff(5*RED_BUFF_SIZE+j),zr*ded%z )
         end if
#endif
         end associate
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
      integer ned_
#endif
        ! Load atom block i parameters
        if (ii.lt.nb2p) then
           iblock = b2pl(ii*2+1)
        else
           iblock = bapl(ii-nb2p+1)
        end if
        if (iblock.eq.0) cycle
        idx    = (iblock-1)*warpsize + ilane 
        iglob  = sgl_id(idx)
        xi     = x(idx)
        yi     = y(idx)
        zi     = z(idx)
        ai     = adisp(iglob)
        ci     = csix (iglob)
#if __tfea__ & __use_softcore__
        muti   = mut(iglob)
#endif
        ! Load atom block k neighbor parameter
        if (ii.lt.nb2p) then
           kb  = b2pl(ii*2+2)
        kdx    = (kb-1)*warpsize + ilane
        else
        kdx    = abpl( (ii-nb2p)*warpsize+ilane )
        if(kdx.eq.0) kdx = nab
        end if
        kglob  = sgl_id(kdx)
        xk     = x(kdx)
        yk     = y(kdx)
        zk     = z(kdx)
        ak     = adisp(kglob)
        ck     = csix (kglob)
#if __tfea__ & __use_softcore__
        mutk   = mut(kglob)
#endif
        ! Set Data to compute to zero
#if __tver__ & __use_ene__
        ed_    = 0
#endif
#if __tver__ & __use_act__
        ned_    = 0
#endif
#if __tver__ & __use_grd__
        gxi    = 0; gyi    = 0; gzi    = 0;
        gxk    = 0; gyk    = 0; gzk    = 0;
#endif
#if __tver__ & __use_vir__
        vxx_   = 0; vxy_   = 0; vxz_   = 0;
        vyy_   = 0; vyz_   = 0; vzz_   = 0;
#endif
        ! Study blocks proximity
        j      = b_stat(iblock)
        info   = merge(0,1,j.eq.i_bl)
        if (info.and.ii.lt.nb2p.and.j.eq.n_bl) then
           j = b_stat(kb)
           info = merge(0,1,j.eq.i_bl)
           if (j.eq.n_bl) then
              xr   = b_rmid(1,iblock) - b_rmid(1,kb)
              yr   = b_rmid(2,iblock) - b_rmid(2,kb)
              zr   = b_rmid(3,iblock) - b_rmid(3,kb)
              info = 0
              if( f_abs(xr).gt.xcell2 ) info=1
              if( f_abs(yr).gt.ycell2 ) info=1
              if( f_abs(zr).gt.zcell2 ) info=1
           end if
        end if

        if (skipvdw12) then; block; integer k
           do k = 1,4
              ai12(k) = i12(k,iglob)
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
              do k = 1, 4
                 if (ai12(k).eq.kglob) ik12=.true.
              end do
              if (ik12) goto 10
           end block; end if
           if (info) then
#if __tfea__ & __use_mpi__
              if (nproc.gt.1) then
                 block; real(t_p) xk_,yk_,zk_
                 xk_  = xk
                 yk_  = yk
                 zk_  = zk
                 xr   = xi - xk_
                 yr   = yi - yk_
                 zr   = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xr,yr,zr)
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
                 xr   = xi - xk
                 yr   = yi - yk
                 zr   = zi - zk
                 if(info) call image_inl(xr,yr,zr)
#if __tfea__ & __use_mpi__
              end if
#endif
           else
#if __tfea__ & __use_mpi__
              if (nproc.gt.1) then
                 block; real(t_p) xk_,yk_,zk_
                 xk_  = xk
                 yk_  = yk
                 zk_  = zk
                 xr   = xi - xk_
                 yr   = yi - yk_
                 zr   = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xr,yr,zr)
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
                 xr   = xi - xk
                 yr   = yi - yk
                 zr   = zi - zk
#if __tfea__ & __use_mpi__
              end if
#endif
           end if
#if __tver__ & __use_grd__
           dedx = 0.0; dedy = 0.0; dedz = 0.0;
#endif
           rik2 = xr**2 + yr**2 + zr**2

           do_pair = merge(iglob.lt.kglob,.true.,idx.eq.kdx)
     &                .and.accept_mid
#if __tfea__ & __use_longRange__
           do_pair = do_pair.and.(sCut2<rik2.and.rik2<=off2)
#else
           do_pair = do_pair.and.rik2<=off2
#endif
           if (do_pair) then
#if __tfea__ & __use_groups__
              if (i_grp) call groups2_inl(fgrp,iglob,kglob,grplist,wgrp)
#endif
#if __tfea__ & __use_softcore__
              mutik   = muti + mutk
              if (vcouple.and.mutik.eq.two1) mutik=one1 ! Annihilation
#endif

              ! Compute paiwise dispersion interaction
#if __tfea__ & __use_ewald__
              call duo_disp_real(rik2,xr,yr,zr,ai,ci,ak,ck
     &                          ,aewald,scale,mutik,vlambda,i_grp,fgrp
     &                          ,sCut,shortheal,e,ded,ver0,fea)
#else
              call duo_disp(rik2,xr,yr,zr,ai,ci,ak,ck,rinv,sCut,off
     &                     ,shortheal,mutik,vlambda,i_grp,fgrp,scale
     &                     ,e,ded,ver,fea)
#endif

#if __tver__ & __use_ene__
              ed_   = ed_ + tp2enr(e)
#endif
#if __tver__ & __use_act__
              ned_  = ned_ + 1
#endif
#if __tver__ & __use_vir__
              if (use_virial) then
              vxx_  = vxx_  + xr * dedx
              vxy_  = vxy_  + yr * dedx
              vxz_  = vxz_  + zr * dedx
              vyy_  = vyy_  + yr * dedy
              vyz_  = vyz_  + zr * dedy
              vzz_  = vzz_  + zr * dedz
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
           ak    = __shfl( ak   ,il )
           ck    = __shfl( ck   ,il )
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
        call atomic_add_f( ev_buff(j),ed_ )
#endif
#if __tver__ & __use_grd__
        !associate(kbis=>kglob,i=>iglob)
        !kbis   = slc_id(kdx)
        !i      = slc_id (idx)
        !increment the van der Waals derivatives
        if (idx.le.na) then
           call atomic_add_f( dedspx(idx),gxi )
           call atomic_add_f( dedspy(idx),gyi )
           call atomic_add_f( dedspz(idx),gzi )
        end if
        if (kdx.le.na) then
           call atomic_sub_f( dedspx(kdx),gxk )
           call atomic_sub_f( dedspy(kdx),gyk )
           call atomic_sub_f( dedspz(kdx),gzk )
        end if
        !end associate
#endif
#if __tver__ & __use_act__
        !increment the van der Waals energy
        call atomic_add_i( nred_buff(j),nev_ )
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
