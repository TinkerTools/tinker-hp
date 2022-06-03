c
c      Template version of ~echg_lj1_kcu~ routine
c      Compute Fused Lennard-Jones and Point Charge interactions in the
c      same kernel
c      version __tver__ features __tfea__
c
        attributes(global)
     &  subroutine __Cat( echg_lj1_kcu, __sufx__ )
     &           (x_sp,y_sp,z_sp,glob_sp,loc_sp,ichg_
     &           ,iblst,blst,jvdw_sp,i12,epsilon,radmin,radv,epsv
     &           ,d_x,d_y,d_z,ev_buff,vir_buff
     &           ,nlocnlb_pair,n,nbloc,nlocnl,nlocnlb
     &           ,nvdwclass,radrule,epsrule
     &           ,vcut2,vcut,voff2,voff,rinv,sheal
     &           ,coff2,f,aewald,aewaldop,ebuffer,pchg
     &           ! Scaling factor
     &           ,x,y,z,correct_ik,vcorrect_scale,ccorrect_scale,szc
     &           ,n_scale,loc,jvdw
     &           ,radmin1,epsilon1,radmin4,epsilon4,de
     &           ! Box data
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &           )

        implicit none
        integer  ,value,intent(in):: nlocnlb_pair,n,nbloc
     &           ,nlocnl,nlocnlb,nvdwclass,radrule,epsrule
     &           ,szc,n_scale
        real(t_p),value,intent(in):: coff2,f,aewald,aewaldop,ebuffer
     &           ,vcut2,vcut,rinv,sheal,voff2,voff
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
        ener_rtyp,device :: ev_buff (  RED_BUFF_SIZE)
        real(t_p),device :: vir_buff(6*RED_BUFF_SIZE)
        integer  ,device,intent(in)::glob_sp(nlocnlb)
     &           ,iblst(nlocnlb_pair), ichg_(nlocnlb)
     &           , blst(nlocnlb_pair*(BLOCK_SIZE))
     &           ,loc_sp(nlocnlb),i12(maxvalue,n)
     &           ,correct_ik(szc,2),loc(n),jvdw(n)
        integer  ,device,intent(in)::jvdw_sp(nlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),radv(*),epsv(*)
     &           ,vcorrect_scale(*),ccorrect_scale(*)
     &           ,radmin1(maxclass,*),epsilon1(maxclass,*)
     &           ,radmin4(maxclass,*),epsilon4(maxclass,*)
        real(t_p),device,intent(in):: x_sp(nlocnlb)
     &           ,y_sp(nlocnlb),z_sp(nlocnlb),x(*),y(*),z(*),pchg(*)
        mdyn_rtyp,device,intent(inout)::d_x(*),d_y(*),d_z(*),de(3,*)

        integer ithread,iwarp,nwarp,ilane,ver,ver1,fea
        integer ii,j, no_scale, yes_scale
        integer idx,kdx,iglob,kglob
        logical do_pair
        real(t_p) xi,yi,zi,xpos,ypos,zpos,xk,yk,zk,rik2,rik
        real(t_p) fi,fk,zer,one
        real(t_p) e,dev_,dec_
        ener_rtyp ev_
        mdyn_rtyp gxi,gyi,gzi
        mdyn_rtyp gxk,gyk,gzk
#if __tver__ & __use_vir__
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
#endif
        parameter(no_scale=0,yes_scale=1,zer=0,one=1
     &           ,ver=__tver__
     &           ,ver1=ver+__use_sca__
     &           ,fea=__tfea__)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        !nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1

        idx     = iand(ithread-1,RED_BUFF_SIZE-1) + 1
c
c       Scaling factor correction
c       -------------------------
        associate(do_scale4=>do_pair)
        do ii = ithread,n_scale,blockDim%x*gridDim%x
        block
        real(t_p) rvi,epsi,rik_i
        type(real3) ded
           iglob  = correct_ik(ii,1)
           kglob  = correct_ik(ii,2)

           xpos   = x(iglob) - x(kglob)
           ypos   = y(iglob) - y(kglob)
           zpos   = z(iglob) - z(kglob)
           call image_inl(xpos,ypos,zpos)
c
c       decide whether to compute the current interaction
c       and check for an interaction distance less than the cutoff
c
           rik2   = xpos**2 + ypos**2 + zpos**2
#if __tver__ & __use_ene__
           e=0;
#endif
#if __tver__ & __use_grd__
           dev_=0; dec_=0
#endif
           rik    = f_sqrt(rik2)
           rik_i  = rik**(-1)

           if (rik2>voff2) goto 50
           block
           integer   it,kt
           real(t_p) vscale
           integer(1) mutik
           real(t_p)  fgrp,sck,sct,scs,scalpha,galpha
     &               ,lambda,vlambda,glamb,lambdavt,lambdavt1,delambdav_
           logical    ugrp,ulamdyn

              vscale = vcorrect_scale(ii)
              it     = jvdw(iglob)
              kt     = jvdw(kglob)

              if (vscale.lt.0) vscale=1 ! 1-4 Interactions
c20           continue
             !replace 1-4 interactions
              !if (do_scale4) then
              !   rvi  = radmin4 (kt ,it)
              !   epsi = epsilon4(kt ,it)
              !else
                 rvi  =  radmin1(kt ,it)
                 epsi = epsilon1(kt ,it)
              !end if
c
c       compute the energy contribution for this interaction
c
              call duo_lj_(rik2,rik,rik_i,rvi,epsi*vscale,vcut2
     &                    ,rinv,voff,sheal,ugrp,fgrp,mutik
     &                    ,sck,sct,scs,scalpha,galpha,ulamdyn
     &                    ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                    ,delambdav_,e,dev_,ver1,fea)

#if __tver__ & __use_ene__
              e    = -e
#endif
#if __tver__ & __use_grd__
              dev_ = -dev_
#endif
           end block

 50        continue

           if (rik2.le.coff2) then; block
           real(t_p) scale_f,fik
              scale_f =   ccorrect_scale(2*ii+1)
              fik     = f*ccorrect_scale(2*ii+2)
              call duo_chg1(rik,ebuffer,fik,aewald,aewaldop,scale_f
     &                     ,e,dec_,ver1,fea)
           end block
           end if

#if __tver__ & (__use_grd__+__use_vir__)
           block
           real(t_p) de_
              de_   = (dev_+dec_) * rik_i
              ded%x = de_ * xpos
              ded%y = de_ * ypos
              ded%z = de_ * zpos
           end block
#endif
#if __tver__ & __use_ene__
           block; ener_rtyp fstat;
           fstat= atomicadd( ev_buff(idx),tp2enr(e) )
           end block
#endif
           associate( iloc=>iglob,kloc=>kglob )
#if __tfea__ & __use_mpi__
           iloc   = loc(iglob)
           kloc   = loc(kglob)
#endif
#if __tver__ & __use_grd__
           call atomic_sub_f1( de(1,kloc), ded%x )
           call atomic_sub_f1( de(2,kloc), ded%y )
           call atomic_sub_f1( de(3,kloc), ded%z )

           call atomic_add_f1( de(1,iloc), ded%x )
           call atomic_add_f1( de(2,iloc), ded%y )
           call atomic_add_f1( de(3,iloc), ded%z )
#endif
           end associate
#if __tver__ & __use_vir__
           if (use_virial) then
          !increment the total van der Waals energy 
           vxx_= atomicadd( vir_buff(0*RED_BUFF_SIZE+idx),xpos*ded%x )
           vxx_= atomicadd( vir_buff(1*RED_BUFF_SIZE+idx),ypos*ded%x )
           vxx_= atomicadd( vir_buff(2*RED_BUFF_SIZE+idx),zpos*ded%x )
           vxx_= atomicadd( vir_buff(3*RED_BUFF_SIZE+idx),ypos*ded%y )
           vxx_= atomicadd( vir_buff(4*RED_BUFF_SIZE+idx),zpos*ded%y )
           vxx_= atomicadd( vir_buff(5*RED_BUFF_SIZE+idx),zpos*ded%z )
           end if
#endif
        end block
        end do
        end associate
c
c       Vdw Special 1-4 pairwise interactions
c       -------------------------------------
        do ii = ithread,n_scale,blockDim%x*gridDim%x
        block
        integer   it,kt
        real(t_p) vscale,rvi,epsi
        type(real3) ded
        integer(1) mutik
        real(t_p)  fgrp,sck,sct,scs,scalpha,galpha
     &            ,lambda,vlambda,glamb,lambdavt,lambdavt1,delambdav_
        logical    ugrp,ulamdyn
           iglob  = correct_ik(ii,1)
           kglob  = correct_ik(ii,2)
           vscale = vcorrect_scale(ii)
           if (vscale.gt.0) cycle
           it     = jvdw(iglob)
           kt     = jvdw(kglob)
           vscale = -vscale ! 1-4 Interactions

           xpos   = x(iglob) - x(kglob)
           ypos   = y(iglob) - y(kglob)
           zpos   = z(iglob) - z(kglob)
           call image_inl(xpos,ypos,zpos)
c
c       decide whether to compute the current interaction
c       and check for an interaction distance less than the cutoff
c
           rik2   = xpos**2 + ypos**2 + zpos**2
           if (rik2>voff2) cycle

           rvi  = radmin4 (kt ,it)
           epsi = epsilon4(kt ,it)
c
c       compute the energy contribution for this interaction
c
           call duo_lj(rik2,xpos,ypos,zpos,rvi,epsi*vscale,vcut2
     &                ,rinv,voff,sheal,ugrp,fgrp,mutik
     &                ,sck,sct,scs,scalpha,galpha,ulamdyn
     &                ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                ,delambdav_,e,ded,ver+__use_sca__,fea)

#if __tver__ & __use_ene__
           block; ener_rtyp fstat;
           fstat= atomicadd( ev_buff(idx),tp2enr(e) )
           end block
#endif
           associate( iloc=>iglob,kloc=>kglob )
#if __tfea__ & __use_mpi__
           iloc   = loc(iglob)
           kloc   = loc(kglob)
#endif
#if __tver__ & __use_grd__
           call atomic_sub_f1( de(1,kloc), ded%x )
           call atomic_sub_f1( de(2,kloc), ded%y )
           call atomic_sub_f1( de(3,kloc), ded%z )

           call atomic_add_f1( de(1,iloc), ded%x )
           call atomic_add_f1( de(2,iloc), ded%y )
           call atomic_add_f1( de(3,iloc), ded%z )
#endif
           end associate
#if __tver__ & __use_vir__
           if (use_virial) then
          !increment the total van der Waals energy 
           vxx_= atomicadd( vir_buff(0*RED_BUFF_SIZE+idx),xpos*ded%x )
           vxx_= atomicadd( vir_buff(1*RED_BUFF_SIZE+idx),ypos*ded%x )
           vxx_= atomicadd( vir_buff(2*RED_BUFF_SIZE+idx),zpos*ded%x )
           vxx_= atomicadd( vir_buff(3*RED_BUFF_SIZE+idx),ypos*ded%y )
           vxx_= atomicadd( vir_buff(4*RED_BUFF_SIZE+idx),zpos*ded%y )
           vxx_= atomicadd( vir_buff(5*RED_BUFF_SIZE+idx),zpos*ded%z )
           end if
#endif
        end block
        end do
c
c       Compute pairwise
c       ----------------
        do ii = iwarp, nlocnlb_pair-1
     &        , blockDim%x*gridDim%x / warpsize
        block
        integer it,kt,ai12(4)
        real(t_p) rvi,rvk,epsi,epsk,rik_i

           ! Load atom block k neighbor parameter
           kdx    =  blst( ii*warpsize + ilane )
           kglob  = glob_sp(kdx)
           xk     = x_sp  (kdx)
           yk     = y_sp  (kdx)
           zk     = z_sp  (kdx)
           fk     = pchg(ichg_(kdx))
#if (__tfea__ & __radepsOpt__ )
           rvk    = radv  (kdx)
           epsk   = epsv  (kdx)
#else
           kt     = jvdw_sp  (kdx)
#endif

           ! Set Data to compute to zero
#if __tver__ & __use_ene__
           ev_    = 0
#endif
#if __tver__ & __use_grd__
           gxi    = 0
           gyi    = 0
           gzi    = 0
           gxk    = 0
           gyk    = 0
           gzk    = 0
#endif
#if __tver__ & __use_vir__
           vxx_   = 0
           vxy_   = 0
           vxz_   = 0
           vyy_   = 0
           vyz_   = 0
           vzz_   = 0
#endif
           ! Load atom block i parameters
           idx = iblst(ii+1)
           if (idx.eq.0) cycle
           idx    = (idx-1)*warpsize + ilane 
           iglob  = glob_sp(idx)
           xi     = x_sp  (idx)
           yi     = y_sp  (idx)
           zi     = z_sp  (idx)
           fi     = f*pchg(ichg_(idx))
#if (__tfea__ & __radepsOpt__ )
           rvi    = radv  (idx)
           epsi   = epsv  (idx)
#else
           it     = jvdw_sp(idx)
#endif
           if (skipvdw12) then; block; integer k;
              do k = 1,4; ai12(k) = i12(k,iglob); end do
           end block; end if

           ! Interact block i with block k
           do j = 0,warpsize-1

              do_pair = merge(.true.,iglob.lt.kglob,idx.ne.kdx)
              if (skipvdw12) then; block; logical ik11,ik12
                 ik11 = merge(.false.,.true.,ai12(1).eq.kglob)
                 ik12 = merge(.false.,.true.,ai12(2).eq.kglob)
                 ik11 = merge(.false.,  ik11,ai12(3).eq.kglob)
                 ik11 = merge(.false.,  ik12,ai12(4).eq.kglob)
                 do_pair = do_pair.and.ik12.and.ik11
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
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    do_pair = do_pair.and..false.
                 else
                    do_pair = do_pair.and..true.
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
#if __tver__ & __use_ene__
              e   =0;
#endif
#if __tver__ & __use_grd__
              dev_=0; dec_=0;
#endif
              if ((.not.(rik2.lt.Inf)).or.rik2.gt.max(voff2,coff2)
     &           .or..not.do_pair) goto 20

              rik   = f_sqrt(rik2)
              rik_i = rik**(-1)

              computeInteraction:block
              integer(1) mutik
              real(t_p)  fgrp,sck,sct,scs,scalpha,galpha
     &                  ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                  ,delambdav_
              logical    ugrp,ulamdyn
              if (rik2<=voff2) then
              block
              real(t_p) rv2,eps2

#if (__tfea__ & __radepsOpt__ )
                 call get_rad(rvi,rvk,rv2,radrule)
                 call get_eps(epsi,epsk,eps2,epsrule)
#else
                 rv2   =  radmin (kt,it)
                 eps2  = epsilon (kt,it)
#endif
                 call duo_lj_(rik2,rik,rik_i,rv2,eps2,vcut2
     &                       ,rinv,voff,sheal,ugrp,fgrp,mutik
     &                       ,sck,sct,scs,scalpha,galpha,ulamdyn
     &                       ,lambda,vlambda,glamb,lambdavt,lambdavt1
     &                       ,delambdav_,e,dev_,ver,fea)
              end block
              end if

              if (rik2<=coff2) then
                 call duo_chg1(rik,ebuffer,fi*fk,aewald,aewaldop,oner
     &                        ,e,dec_,ver,fea)
              end if

#if __tver__ & __use_ene__
              ! Accumulate interaction energy
              ev_   = ev_ + tp2enr(e)
#endif
#if __tver__ & __use_grd__
              block
              real(t_p) de_
              type(real3) ded
                 de_   = (dev_+dec_) * rik_i
                 ded%x = de_ * xpos
                 ded%y = de_ * ypos
                 ded%z = de_ * zpos
#if    __tver__ & __use_vir__
                 if (use_virial) then
                    vxx_  = vxx_  + xpos*ded%x
                    vxy_  = vxy_  + ypos*ded%x
                    vxz_  = vxz_  + zpos*ded%x
                    vyy_  = vyy_  + ypos*ded%y
                    vyz_  = vyz_  + zpos*ded%y
                    vzz_  = vzz_  + zpos*ded%z
                 end if
#endif
                 ! Accumulate interaction gradient
                 gxi = gxi + tp2mdr(ded%x)
                 gyi = gyi + tp2mdr(ded%y)
                 gzi = gzi + tp2mdr(ded%z)

                 gxk = gxk + tp2mdr(ded%x)
                 gyk = gyk + tp2mdr(ded%y)
                 gzk = gzk + tp2mdr(ded%z)
              end block
#endif
              end block computeInteraction

 20           continue
              block; integer il
              il   = iand(ilane,warpsize-1)+1
              kglob= __shfl(kglob,il )
#if (__tfea__ & __radepsOpt__ )
              rvk  = __shfl(  rvk,il )
              epsk = __shfl( epsk,il )
#else
              kt   = __shfl(   kt,il )
#endif
              xk   = __shfl(   xk,il )
              yk   = __shfl(   yk,il )
              zk   = __shfl(   zk,il )
              fk   = __shfl(   fk,il )
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
           !i      = loc_sp (idx)
           !kbis   = loc_sp(kdx)
           if (idx.le.nlocnl) then
              call atomic_add_f (d_x(idx),gxi)
              call atomic_add_f (d_y(idx),gyi)
              call atomic_add_f (d_z(idx),gzi)
              !call atomic_add_f (dev(1,i),gxi)
              !call atomic_add_f (dev(2,i),gyi)
              !call atomic_add_f (dev(3,i),gzi)
           end if

           !call syncwarp(ALL_LANES)
           if (kdx.le.nlocnl) then
              call atomic_sub_f (d_x(kdx),gxk)
              call atomic_sub_f (d_y(kdx),gyk)
              call atomic_sub_f (d_z(kdx),gzk)
              !call atomic_sub_f (dev(1,kbis),gxk)
              !call atomic_sub_f (dev(2,kbis),gyk)
              !call atomic_sub_f (dev(3,kbis),gzk)
           end if
#endif
        end block
        end do
        end subroutine

