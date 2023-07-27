c
c     Hal CUDA Kernel Template 
c
c     "ecreal" evaluates the real space portion of the Ewald sum
c     energy and forces due to atomic charge interactions, using a neighbor list
c     version __tver__ features __tfea__
c
        attributes(global) subroutine __Cat(ecreal,__sufx__)
     &        ( iion, cglob, loc, ieblst, eblst
     &        , nionlocnlb, nionlocnlb_pair, nionbloc, n
     &        , x, y, z, pchg, pchg_, mut
     &        , off2,loff2,scut,shortheal,f,aewald,ebuffer
     &        , elambda, use_lambdadyn
     &        , dec, ec_buff, vir_buff, lam_buff, act_buff
     &        , use_group, grplist, wgrp
     &        , correct_ik,correct_scale,chglist,x_,y_,z_,loc_
     &        , n_scale,siz,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &        )
        implicit none
        integer   ,value::nionlocnlb,nionbloc,n
     &            ,nionlocnlb_pair,siz,n_scale
        real(t_p) ,value:: p_xbeg,p_xend,p_ybeg,p_yend
     &            ,p_zbeg,p_zend
     &            ,off2,loff2,scut,shortheal,aewald,f,ebuffer,elambda
        logical   ,value:: use_lambdadyn, use_group
        integer(1),device,intent(in):: mut(n)
        integer   ,device,intent(in)::iion(*),cglob(*),loc(*),loc_(*)
     &            ,ieblst(*),eblst(*),correct_ik(siz,2)
     &            ,grplist(*),chglist(*)
        real(t_p) ,device,intent(in):: x(*),y(*),z(*),pchg(*),pchg_(*)
     &            ,x_(*),y_(*),z_(*),correct_scale(*),wgrp(ngrp+1,*)
        integer   ,device:: act_buff(*)
        real(t_p) ,device:: vir_buff(*)
        real(r_p) ,device:: lam_buff(*)
        ener_rtyp ,device:: ec_buff(*)
        mdyn_rtyp ,device:: dec(3,*)

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,k,nadd
        integer iblock,idx,kdx
        integer iichg,iglob,icloc,kchg,kcloc
        integer location,ver,ver1,fea,RBS
        integer kglob
        real(t_p) xk_,yk_,zk_,d2,fi,fi_
        ener_rtyp ec_
        real(t_p) rstat,one_f
        type(real3) posi,pos
        type(mdyn3_r) frc_i
#if __tfea__ & (__use_lambdadyn__+__use_softcore__)
        integer(1)     ,shared:: mutk(BLOCK_DIM)
        real(t_p)      ,shared::  fk_(BLOCK_DIM)
#endif
        real(t_p)      ,shared::   fk(BLOCK_DIM)
        type(mdyn3_r)  ,shared::frc_k(BLOCK_DIM)
        type(real3)    ,shared:: posk(BLOCK_DIM)
        real(t_p) vir_(6)
        logical  do_pair,same_block,accept_mid
        parameter( ver  = __tver__, one_f=1.0
     &           , ver1 = __tver__+__use_sca__
     &           , fea  = __tfea__
     &           , RBS  = RED_BUFF_SIZE
     &           )

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.
c
c     compute scaling interactions
c
        do ii = ithread, n_scale, blockDim%x*gridDim%x
        block; 
        integer(1) mutik
        real(t_p)  scale_f,fik_,fik,delambdae_,e
        type(real3) ded
           iglob   = correct_ik(ii,1)
           kglob   = correct_ik(ii,2)
#if __tfea__ & __use_lambdadyn__
           if (use_lambdadyn) then
           scale_f = correct_scale(2*ii+1)
           fik     = f*pchg(chglist(iglob))*pchg(chglist(kglob))
           fik_    = f*pchg_(chglist(iglob))*pchg_(chglist(kglob))
           mutik   = mut(iglob)+mut(kglob)
           else
#endif
           scale_f =   correct_scale(2*ii+1)
           fik     = f*correct_scale(2*ii+2)
#if __tfea__ & __use_lambdadyn__
           end if
#endif
           posi%x  = x_(iglob)
           posi%y  = y_(iglob)
           posi%z  = z_(iglob)
           pos%x   = posi%x - x_(kglob)
           pos%y   = posi%y - y_(kglob)
           pos%z   = posi%z - z_(kglob)
           call image_inl (pos%x,pos%y,pos%z)
           d2 = pos%x**2 + pos%y**2 + pos%z**2
#if __tfea__ & __use_longRange__
           if (d2.gt.loff2 .and. d2.le.off2) then
#else
           if (d2.le.off2) then
#endif
#if __tfea__ & __use_groups__
              if (use_group) then; block
              real(t_p) fgrp
                 call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
                 scale_f = 1.0 - scale_f*(1.0-fgrp)
              end block; end if
#endif
             !correct pair
             !compute the pairwise contribution for this interaction
              call charge_couple(d2,pos%x,pos%y,pos%z,ebuffer
     &                          ,fik_,fik,aewald,scale_f,mutik
     &                          ,use_lambdadyn,shortheal,scut,elambda
     &                          ,delambdae_,e,ded,ver1,fea)

              block; integer lot
              lot   = iand( ithread-1,RBS-1 ) + 1
#if __tver__ & __use_ene__
              ! increment the overall energy
              call atomic_add_f( ec_buff(lot),tp2enr(e) )
#endif
#if __tver__ & __use_act__
              if(scale_f.eq.-1.0) call atomic_sub_i( act_buff(lot),1 )
#endif
#if __tfea__ & __use_lambdadyn__
              if (use_lambdadyn) 
     &           call atomic_add_m(lam_buff(lot),delambdae_)
#endif
#if __tver__ & __use_vir__
              ! increment the internal virial tensor components
              call atomic_add_c(vir_buff(0*RBS+lot),pos%x*ded%x)
              call atomic_add_c(vir_buff(1*RBS+lot),pos%y*ded%x)
              call atomic_add_c(vir_buff(2*RBS+lot),pos%z*ded%x)
              call atomic_add_c(vir_buff(3*RBS+lot),pos%y*ded%y)
              call atomic_add_c(vir_buff(4*RBS+lot),pos%z*ded%y)
              call atomic_add_c(vir_buff(5*RBS+lot),pos%z*ded%z)
#endif
              end block
#if __tver__ & __use_grd__
#if __tfea__ & __use_mpi__
              i  = loc_(iglob)
              k  = loc_(kglob)
#else
              i  = iglob
              k  = kglob
#endif
              ! increment derivative expressions
              call atomic_add_f1( dec(1,i),ded%x )
              call atomic_add_f1( dec(2,i),ded%y )
              call atomic_add_f1( dec(3,i),ded%z )
              call atomic_sub_f1( dec(1,k),ded%x )
              call atomic_sub_f1( dec(2,k),ded%y )
              call atomic_sub_f1( dec(3,k),ded%z )
#endif
           end if
        end block
        end do
c
c     compute the real space Ewald energy and first derivatives
c
        do ii = iwarp, nionlocnlb_pair, nwarp
        block;
        integer(1) muti
        real(r_p)  delambdae0
           iblock = ieblst(ii)
           if (iblock==0) cycle

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kchg    = cglob(kdx)
           kglob   = iion (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           fk  (threadIdx%x)    = pchg(kchg)
#if __tfea__ & (__use_lambdadyn__+__use_softcore__)
           if (use_lambdadyn) then
              mutk(threadIdx%x) = mut(kglob)
               fk_(threadIdx%x) = pchg_(kchg)
           end if
#endif
           !  Load atom block i parameters
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iichg   = cglob(idx)
           iglob   = iion (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           fi      = f*pchg(iichg)
#if __tfea__ & (__use_lambdadyn__+__use_softcore__)
           if (use_lambdadyn) then
              muti    = mut(iglob)
              fi_     = f*pchg_(iichg)
           end if
#endif
#if __tver__ & __use_grd__
           ! zero data to compute
           frc_i%x = 0.0;
           frc_i%y = 0.0;
           frc_i%z = 0.0;
           frc_k(threadIdx%x)%x = 0.0;
           frc_k(threadIdx%x)%y = 0.0;
           frc_k(threadIdx%x)%z = 0.0;
#endif
           !* set compute Data to 0
#if __tver__ & __use_ene__
           ec_ = 0
#endif
#if __tver__ & __use_act__
           nadd = 0
#endif
#if __tfea__ & __use_lambdadyn__
           delambdae0 = 0
#endif
#if __tver__ & __use_vir__
           vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
           vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;
#endif

           same_block = (idx.ne.kdx)

           do j = 0,warpsize-1
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#if __tfea__ & __use_mpi__
              if (ndir.gt.1) then
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
                 end if
              else
#endif
                 pos%x = posi%x - posk(klane)%x 
                 pos%y = posi%y - posk(klane)%y 
                 pos%z = posi%z - posk(klane)%z 
                 call image_inl(pos%x,pos%y,pos%z)
#if __tfea__ & __use_mpi__
              end if
#endif
              d2       = pos%x**2 + pos%y**2 + pos%z**2
              do_pair  = merge(.true.,iglob.lt.kglob
     &                       ,same_block)
#if __tfea__ & __use_longRange__
              if (do_pair.and.d2.ge.loff2.and.d2.le.off2.and.accept_mid)
     &           then
#else
              if (do_pair.and.d2.le.off2.and.accept_mid) then
#endif
              ! compute one interaction
              block
              integer(1) mutik
              real(t_p) e, delambdae_
              type(real3) frc
              real(t_p) fgrp
#if !(__tfea__ & __use_groups__)
              parameter( fgrp=1.0 )
#endif
#if __tfea__ & __use_groups__
                 if (use_group) then
                    call groups2_inl(fgrp,iglob,kglob,ngrp,grplist,wgrp)
                    fgrp = 1.0-fgrp
                 else
                    fgrp = 0.0
                 end if
#endif 
#if __tfea__ & (__use_lambdadyn__+__use_softcore__)
                 mutik = muti + mutk(klane)
#endif
                 call charge_couple(d2,pos%x,pos%y,pos%z,ebuffer
#if __tfea__ & (__use_lambdadyn__+__use_softcore__)
     &                      ,fi_*fk_(klane),fi*fk(klane),aewald
#else
     &                      ,fi_,fi*fk(klane),aewald
#endif
     &                      ,fgrp,mutik,use_lambdadyn,shortheal,scut
     &                      ,elambda,delambdae_,e,frc,ver,fea)

#if __tver__ & __use_ene__
                 ec_     =  ec_ + tp2enr(e)
#endif
#if __tver__ & __use_act__
                 nadd    = nadd +1
#endif
#if __tfea__ & __use_lambdadyn__
                 delambdae0 = delambdae0 + delambdae_
#endif
#if __tver__ & __use_vir__
                 vir_(1) = vir_(1) + pos%x * frc%x
                 vir_(2) = vir_(2) + pos%y * frc%x
                 vir_(3) = vir_(3) + pos%z * frc%x
                 vir_(4) = vir_(4) + pos%y * frc%y
                 vir_(5) = vir_(5) + pos%z * frc%y
                 vir_(6) = vir_(6) + pos%z * frc%z
#endif
#if __tver__ & __use_grd__
                 frc_i%x = frc_i%x + tp2mdr( frc%x )
                 frc_i%y = frc_i%y + tp2mdr( frc%y )
                 frc_i%z = frc_i%z + tp2mdr( frc%z )
                 frc_k(klane)%x = frc_k(klane)%x - tp2mdr( frc%x )
                 frc_k(klane)%y = frc_k(klane)%y - tp2mdr( frc%y )
                 frc_k(klane)%z = frc_k(klane)%z - tp2mdr( frc%z )
#endif
              end block
              end if
              block; integer il
              il   = iand(ilane,warpsize-1)+1
              kglob= __shfl(kglob,il )
              end block
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
#if __tfea__ & __use_lambdadyn__
           if (use_lambdadyn)
     &        call atomic_add_m1(lam_buff(location),delambdae0)
#endif
#if __tver__ & __use_act__
           call atomic_add_i( act_buff(location),nadd )
#endif
#if __tver__ & __use_ene__
           ! Update energy buffer
           call atomic_Add_f(ec_buff(location), ec_)
#endif
#if __tver__ & __use_vir__
           ! Update virial buffer
           call atomic_Add_c(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           call atomic_Add_c(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           call atomic_Add_c(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           call atomic_Add_c(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           call atomic_Add_c(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           call atomic_Add_c(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))
#endif
#if __tver__ & __use_grd__
           ! Update forces
           i       = loc  (idx)
           k       = loc  (kdx)
           call atomic_add_f( dec(1,i), frc_i%x )
           call atomic_add_f( dec(2,i), frc_i%y )
           call atomic_add_f( dec(3,i), frc_i%z )
           call atomic_add_f( dec(1,k), frc_k(threadIdx%x)%x )
           call atomic_add_f( dec(2,k), frc_k(threadIdx%x)%y )
           call atomic_add_f( dec(3,k), frc_k(threadIdx%x)%z )
#endif
        end block
        end do

        end subroutine
