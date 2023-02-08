c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "epolar_cpencu" : CUDA Template for calculation of the real space portion of the Ewald
c     summation energy and gradient due to dipole polarization
c     via a block-atom list with charge penetration
c     version __tver__ features __tfea__
c
c
#if !(__tfea__&__use_chgpen__)
#  error  __FILE__ is specific to charge penetration feature
#endif
        attributes(global) subroutine __Cat(epreal_cpen,__sufx__)
     &        ( ipole,pglob,loc,ploc,bapl,abpl
     &        , x,y,z,rpole,uind,uinp
     &        , pcore,pval,palpha
     &        , dep,trq,pot,ep_buff,vir_buff,nep_buff
     &        , na,nab,nabp,npolebloc,n,idec,pentyp
     &        , off2,f,ewald,ucflx
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &        , scal_ik,scal_val,n_scal,ipole_,polelocnl,loc_,x_,y_,z_
     &        )

        implicit  none
        integer  ,value,intent(in):: na,nab,npolebloc,n,nabp,n_scal
     &           ,pentyp,idec
        logical  ,value,intent(in):: ucflx
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &           ,off2,alsq2,alsq2n,ewald,f
        integer  ,device,intent(in)::ipole(*),pglob(*),loc(*),ploc(*)
     &           ,bapl(*),abpl(*),ipole_(*),loc_(*),polelocnl(*)
     &           ,scal_ik(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
     &           ,pcore(*),pval(*),palpha(*),uind(3,*),uinp(3,*),pot(*)
     &           ,scal_val(*),x_(*),y_(*),z_(*)
        integer  ,device:: nep_buff(*)
        real(t_p),device:: trq(3,*),vir_buff(*)
        mdyn_rtyp,device:: dep(3,*)
        ener_rtyp,device:: ep_buff(*)

        integer   ithread,iwarp,nwarp,ilane,klane,istat,srclane
        integer   beg,ii,j,i,k,ver,ver0,fea,iblock,idx,kdx,kdx_,lot
        integer   iipole,iglob,iploc,kpole,kploc,location,nep_
        logical   do_pair,same_block,accept_mid
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) ep_,vir_(6)
        real(t_p) corei,corek,vali,valk,alphai,alphak,poti,potk
        type(real3    ) pi,d,ui,uip,trqi,frc
        type(mdyn3_r  ) frc_i
        type(rpole_elt) ip

        integer        ,shared:: kglob(BLOCK_DIM)
        type(real3)    ,shared,dimension(BLOCK_DIM):: pk,uk,ukp,trqk
        type(mdyn3_r)  ,shared:: frc_k(BLOCK_DIM)
        type(rpole_elt),shared:: kp(BLOCK_DIM)

        parameter(ver=__tver__,ver0=ver+__use_sca__,fea=__tfea__)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        if (idec.eq.1) then
c
c     Loop over pairwise scaled interactions 
c
        do ii = ithread, n_scal, blockDim%x*gridDim%x; block
        real(t_p) dscal,pscal,wscal,e
           iipole  = scal_ik(2*(ii-1)+1)
           kpole   = scal_ik(2*(ii-1)+2)

           iglob   = ipole_(iipole)
           kglob(threadIdx%x)   = ipole_(kpole)

           d%x     = x_(kglob(threadIdx%x)) - x_(iglob)
           d%y     = y_(kglob(threadIdx%x)) - y_(iglob)
           d%z     = z_(kglob(threadIdx%x)) - z_(iglob)

           call image_inl(d%x,d%y,d%z)
           ! cutoff
           d2      = d%x**2 + d%y**2 + d%z**2
           if (d2>off2) cycle

           dscal   = scal_val(3*(ii-1)+1)
           pscal   = scal_val(3*(ii-1)+2)
           wscal   = scal_val(3*(ii-1)+3)

            corei  = pcore(   iipole)
             vali  = pval (   iipole)
           alphai  = palpha(  iipole)

           ip%c    = rpole( 1,iipole)
           ip%dx   = rpole( 2,iipole)
           ip%dy   = rpole( 3,iipole)
           ip%dz   = rpole( 4,iipole)
           ip%qxx  = rpole( 5,iipole)
           ip%qxy  = rpole( 6,iipole)
           ip%qxz  = rpole( 7,iipole)
           ip%qyy  = rpole( 9,iipole)
           ip%qyz  = rpole(10,iipole)
           ip%qzz  = rpole(13,iipole)

           ui%x    = uind ( 1,iipole)
           ui%y    = uind ( 2,iipole)
           ui%z    = uind ( 3,iipole)
           uip%x   = uinp ( 1,iipole)
           uip%y   = uinp ( 2,iipole)
           uip%z   = uinp ( 3,iipole)

            corek  = pcore(    kpole)
             valk  = pval (    kpole)
           alphak  = palpha(   kpole)

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

           uk(threadIdx%x)%x    = uind ( 1, kpole)
           uk(threadIdx%x)%y    = uind ( 2, kpole)
           uk(threadIdx%x)%z    = uind ( 3, kpole)
           ukp(threadIdx%x)%x   = uinp ( 1, kpole)
           ukp(threadIdx%x)%y   = uinp ( 2, kpole)
           ukp(threadIdx%x)%z   = uinp ( 3, kpole)

           ! zero outputs
#if (__tver__ & __use_ene__)
           e = 0;
#endif
#if ( __tver__ & __use_grd__ )
            frc%x=0;  frc%y=0;  frc%z=0;
           trqk%x=0; trqk%y=0; trqk%z=0;
           trqi%x=0; trqi%y=0; trqi%z=0;
#endif

           ! Compute polar interaction
           call duo_polar_cpen(ui,uip,ip
     &             ,uk(threadIdx%x),ukp(threadIdx%x),kp(threadIdx%x)
     &             ,d,d2,f,corei,corek,vali,valk,alphai,alphak
     &             ,dscal,pscal,wscal,ewald,pentyp,ucflx
     &             ,poti,potk,e,frc,trqi,trqk(threadIdx%x),ver0,fea)

           lot     = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
#if (__tver__ & __use_ene__)
           ! Increment energy
           call atomic_add_f( ep_buff(lot), tp2enr(e) )
#endif
#if ( __tver__ & __use_grd__ )
           i       = loc_ (iglob)
           k       = loc_ (kglob(threadIdx%x))
           iploc   = polelocnl(iipole)
           kploc   = polelocnl(kpole)
           ! Increment gradient due to Cartesian forces
           call atomic_add_f1( dep(1,i),-frc%x )
           call atomic_add_f1( dep(2,i),-frc%y )
           call atomic_add_f1( dep(3,i),-frc%z )
           call atomic_add_f1( dep(1,k), frc%x )
           call atomic_add_f1( dep(2,k), frc%y )
           call atomic_add_f1( dep(3,k), frc%z )
           ! Increment torque
           call atomic_add_c( trq(1,iploc),trqi%x )
           call atomic_add_c( trq(2,iploc),trqi%y )
           call atomic_add_c( trq(3,iploc),trqi%z )
           call atomic_add_c( trq(1,kploc),trqk(threadIdx%x)%x )
           call atomic_add_c( trq(2,kploc),trqk(threadIdx%x)%y )
           call atomic_add_c( trq(3,kploc),trqk(threadIdx%x)%z )
#endif
#if ( __tfea__ & __use_chgflx__ )
           call atomic_add_c( pot(i),poti )
           call atomic_add_c( pot(k),potk )
#endif
#if ( __tver__ & __use_vir__ )
           ! Increment virial
           call atomic_add_c(vir_buff(0*RED_BUFF_SIZE+lot),d%x*frc%x)
           call atomic_add_c(vir_buff(1*RED_BUFF_SIZE+lot)
     &                      ,d%y*frc%x+d%x*frc%y)
           call atomic_add_c(vir_buff(2*RED_BUFF_SIZE+lot)
     &                      ,d%z*frc%x+d%x*frc%z)
           call atomic_add_c(vir_buff(3*RED_BUFF_SIZE+lot),d%y*frc%y)
           call atomic_add_c(vir_buff(4*RED_BUFF_SIZE+lot)
     &                      ,d%z*frc%y+d%y*frc%z)
           call atomic_add_c(vir_buff(5*RED_BUFF_SIZE+lot),d%z*frc%z)
#endif
        end block; end do

        end if
c
c     Compute pairwise without scaling correction using nblist (C2)
c
        do ii = iwarp, nabp, nwarp; block
        real(t_p) e
           iblock = bapl(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           pi%x    =     x(idx)
           pi%y    =     y(idx)
           pi%z    =     z(idx)
            corei  = pcore(   iipole)
             vali  = pval (   iipole)
           alphai  = palpha(  iipole)
           ip%c    = rpole( 1,iipole)
           ip%dx   = rpole( 2,iipole)
           ip%dy   = rpole( 3,iipole)
           ip%dz   = rpole( 4,iipole)
           ip%qxx  = rpole( 5,iipole)
           ip%qxy  = rpole( 6,iipole)
           ip%qxz  = rpole( 7,iipole)
           ip%qyy  = rpole( 9,iipole)
           ip%qyz  = rpole(10,iipole)
           ip%qzz  = rpole(13,iipole)
             ui%x  = uind ( 1,iipole)
             ui%y  = uind ( 2,iipole)
             ui%z  = uind ( 3,iipole)
            uip%x  = uinp ( 1,iipole)
            uip%y  = uinp ( 2,iipole)
            uip%z  = uinp ( 3,iipole)

           !  Load atom block k parameters
           kdx     = abpl( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob(threadIdx%x)  = ipole(kdx)
            corek  = pcore(    kpole)
             valk  = pval (    kpole)
           alphak  = palpha(   kpole)
           pk(threadIdx%x)%x   = x(kdx)
           pk(threadIdx%x)%y   = y(kdx)
           pk(threadIdx%x)%z   = z(kdx)
           kp(threadIdx%x)%c   = rpole( 1, kpole)
           kp(threadIdx%x)%dx  = rpole( 2, kpole)
           kp(threadIdx%x)%dy  = rpole( 3, kpole)
           kp(threadIdx%x)%dz  = rpole( 4, kpole)
           kp(threadIdx%x)%qxx = rpole( 5, kpole)
           kp(threadIdx%x)%qxy = rpole( 6, kpole)
           kp(threadIdx%x)%qxz = rpole( 7, kpole)
           kp(threadIdx%x)%qyy = rpole( 9, kpole)
           kp(threadIdx%x)%qyz = rpole(10, kpole)
           kp(threadIdx%x)%qzz = rpole(13, kpole)

             uk(threadIdx%x)%x = uind ( 1, kpole)
             uk(threadIdx%x)%y = uind ( 2, kpole)
             uk(threadIdx%x)%z = uind ( 3, kpole)
            ukp(threadIdx%x)%x = uinp ( 1, kpole)
            ukp(threadIdx%x)%y = uinp ( 2, kpole)
            ukp(threadIdx%x)%z = uinp ( 3, kpole)

           !* set compute Data to 0
#if (__tver__ & __use_grd__)
           frc_i%x = 0;
           frc_i%y = 0;
           frc_i%z = 0;
           frc_k(threadIdx%x)%x = 0;
           frc_k(threadIdx%x)%y = 0;
           frc_k(threadIdx%x)%z = 0;
           trqi%x  = 0.0_ti_p;
           trqi%y  = 0.0_ti_p;
           trqi%z  = 0.0_ti_p;
           trqk(threadIdx%x)%x = 0.0_ti_p;
           trqk(threadIdx%x)%y = 0.0_ti_p;
           trqk(threadIdx%x)%z = 0.0_ti_p;
#endif
#if (__tver__ & __use_ene__)
           ep_ = 0
#endif
#if (__tver__ & __use_act__)
           nep_ = 0
#endif
#if ( __tfea__ & __use_chgflx__ )
           poti = 0
           potk = 0
#endif
#if (__tver__ & __use_vir__)
           vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
           vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;
#endif
           same_block = (idx.ne.kdx)

           do j = 0,warpsize-1
              srclane  = iand( ilane-1+j,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#if (__tfea__ & __use_mpi__)
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
#if (__tfea__ & __use_mpi__)
              end if
#endif
              d2      = d%x**2 + d%y**2 + d%z**2
              do_pair = merge(.true.,iglob.lt.kglob(klane),same_block)
              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call duo_polar_cpen(ui,uip,ip,uk(klane),ukp(klane)
     &                   ,kp(klane),d,d2,f,corei,corek,vali,valk
     &                   ,alphai,alphak,1.0,1.0,1.0,ewald,pentyp,ucflx
     &                   ,poti,potk,e,frc,trqi,trqk(klane)
     &                   ,ver,fea)

#if (__tver__ & __use_ene__)
                 ep_ = ep_ + tp2enr(e)
#endif
#if (__tver__ & __use_act__)
                 nep_ = nep_ + 1
#endif
#if (__tver__ & __use_vir__)
                 vir_(1) = vir_(1) + d%x*frc%x
                 vir_(2) = vir_(2) + 0.5*(d%y*frc%x + d%x*frc%y)
                 vir_(3) = vir_(3) + 0.5*(d%z*frc%x + d%x*frc%z)
                 vir_(4) = vir_(4) + d%y*frc%y
                 vir_(5) = vir_(5) + 0.5*(d%z*frc%y + d%y*frc%z)
                 vir_(6) = vir_(6) + d%z*frc%z
#endif
#if (__tver__ & __use_grd__)
                 frc_i%x = frc_i%x - tp2mdr(frc%x)
                 frc_i%y = frc_i%y - tp2mdr(frc%y)
                 frc_i%z = frc_i%z - tp2mdr(frc%z)
             frc_k(threadIdx%x)%x = frc_k(threadIdx%x)%x + tp2mdr(frc%x)
             frc_k(threadIdx%x)%y = frc_k(threadIdx%x)%y + tp2mdr(frc%y)
             frc_k(threadIdx%x)%z = frc_k(threadIdx%x)%z + tp2mdr(frc%z)
#endif
              end if

                lot = iand(ilane,warpsize-1)+1
               corek= __shfl( corek,lot )
                valk= __shfl(  valk,lot )
              alphak= __shfl(alphak,lot )
#if ( __tfea__ & __use_chgflx__ )
                potk= __shfl(  potk,lot )
#endif
           end do

           lot = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

#if (__tver__ & __use_ene__)
           call atomic_add_f(ep_buff(lot), tp2enr(ep_))
#endif
#if (__tver__ & __use_act__)
           call atomic_add_i(nep_buff(lot),nep_ )
#endif
#if (__tver__ & __use_grd__)
           iploc   = ploc (idx)
           kploc   = ploc (kdx)
           i       = loc  (idx)
           k       = loc  (kdx)
           if (idx.le.na) then
              call atomic_add_f( dep(1,i   ), frc_i%x )
              call atomic_add_f( dep(2,i   ), frc_i%y )
              call atomic_add_f( dep(3,i   ), frc_i%z )
              call atomic_add_c( trq(1,iploc),trqi%x )
              call atomic_add_c( trq(2,iploc),trqi%y )
              call atomic_add_c( trq(3,iploc),trqi%z )
           end if
           if (kdx.le.na) then
              call atomic_add_f( dep(1,k    ),frc_k(threadIdx%x)%x )
              call atomic_add_f( dep(2,k    ),frc_k(threadIdx%x)%y )
              call atomic_add_f( dep(3,k    ),frc_k(threadIdx%x)%z )
              call atomic_add_c( trq(1,kploc),trqk(threadIdx%x)%x )
              call atomic_add_c( trq(2,kploc),trqk(threadIdx%x)%y )
              call atomic_add_c( trq(3,kploc),trqk(threadIdx%x)%z )
           end if
#endif
#if ( __tfea__ & __use_chgflx__ )
           if (idx.le.na) call atomic_add_c( pot(i),poti )
           if (kdx.le.na) call atomic_add_c( pot(k),potk )
#endif
#if (__tver__ & __use_vir__)
           call atomic_add_c(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           call atomic_add_c(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           call atomic_add_c(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           call atomic_add_c(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           call atomic_add_c(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           call atomic_add_c(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))
#endif
        end block; end do
        end subroutine

