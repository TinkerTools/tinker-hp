#include "image.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"

c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine chkttor  --  check torsion-torsion chirality  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "chkttor" tests the attached atoms at a torsion-torsion central
c     site and inverts the angle values if the site is chiral
c
c     note that the sign convention used in this version is correct
c     for phi-psi torsion-torsion interactions as defined in the
c     AMOEBA protein force field; the code may need to be altered
c     for other uses of the torsion-torsion potential, and will not
c     correctly handle enantiomeric sugar rings in nucleic acids
c
c
      ! Differs from ~chkttor~ only through calling interface
      M_subroutine
     &   chkttor_(xia,xib,xic,xid,yia,yib,yic,yid,zia,zib,zic,zid
     &            ,sign,value1,value2)
!$acc routine seq
      implicit none
      real(t_p),intent(in)::xia,xib,xic,xid,yia,yib,yic,yid
     &         ,zia,zib,zic,zid
      real(t_p),intent(inout):: sign
      real(t_p),intent(inout):: value1
      real(t_p),intent(inout):: value2
      real(t_p) xac,yac,zac
      real(t_p) xbc,ybc,zbc
      real(t_p) xdc,ydc,zdc
      real(t_p) c1,c2,c3,vol

      yac = yia - yic
      zac = zia - zic
      ybc = yib - yic
      zbc = zib - zic
      ydc = yid - yic
      zdc = zid - zic
      c1  = ybc*zdc - zbc*ydc
      c2  = ydc*zac - zdc*yac
      c3  = yac*zbc - zac*ybc
      vol = (xia-xic)*c1 + (xib-xic)*c2 + (xid-xic)*c3
 
      !invert the angle values if chirality has an inverted sign
      if (vol .lt. 0.0) then
         sign   = -1.0
         value1 = -value1
         value2 = -value2
      end if
      end

      ! see original implementation in bicubic.f
      M_subroutine
     &             bcucof (y,y1,y2,y12,d1,d2,c)
!$acc routine seq
      implicit none
      real(t_p) d1,d2,d1d2
      real(t_p) y(4),y1(4),y2(4),y12(4),c(4,4)
      d1d2 = d1*d2

      y1(1) = d1 * y1(1)
      y1(2) = d1 * y1(2)
      y1(3) = d1 * y1(3)
      y1(4) = d1 * y1(4)
      y2(1) = d2 * y2(1)
      y2(2) = d2 * y2(2)
      y2(3) = d2 * y2(3)
      y2(4) = d2 * y2(4)
      y12(1) = d1d2 * y12(1)
      y12(2) = d1d2 * y12(2)
      y12(3) = d1d2 * y12(3)
      y12(4) = d1d2 * y12(4)

      c(1,1) =  y(1)
      c(1,2) = y2(1)
      c(1,3) = 3 * (-y(1) + y(4)) - (2*y2(1) + y2(4))
      c(1,4) = 2 * ( y(1) - y(4)) +    y2(1) + y2(4)

      c(2,1) =  y1(1)
      c(2,2) = y12(1)
      c(2,3) = 3 * (y1(4) - y1(1)) - (2*y12(1) + y12(4))
      c(2,4) = 2 * (y1(1) - y1(4)) +    y12(1) + y12(4)

      c(3,1) = 3 * ( y(2) -  y(1)) - (2* y1(1) +  y1(2))
      c(3,2) = 3 * (y2(2) - y2(1)) - (2*y12(1) + y12(2))
      c(3,3) = 9 *( y (1) -     y  (2) +     y  (3) -     y  (4)) 
     &       + 6 * y1 (1) + 3 * y1 (2) - 3 * y1 (3) - 6 * y1 (4) 
     &       + 6 * y2 (1) - 6 * y2 (2) - 3 * y2 (3) + 3 * y2 (4)
     &       + 4 * y12(1) + 2 * y12(2) +     y12(3) + 2 * y12(4)
      c(3,4) = 6 * (-y(1) +     y(2) -      y(3) +      y(4)) 
     &        -4 * y1(1) - 2 * y1(2) + 2 * y1(3) + 4 * y1(4) 
     &        -3 * y2(1) + 3 * y2(2) + 3 * y2(3) - 3 * y2(4)
     &        -2 *y12(1) -    y12(2) -    y12(3) - 2 *y12(4)

      c(4,1) = 2 * ( y(1) -  y(2)) +  y1(1) +  y1(2)
      c(4,2) = 2 * (y2(1) - y2(2)) + y12(1) + y12(2)
      c(4,3) = 6 * (-y(1) +  y(2)  -   y(3) +   y(4))
     &       + 3 * (y1(3) + y1(4)  -  y1(1) -  y1(2))
     &       + 2 * (2 * (y2(2) - y2(1)) + y2(3) - y2(4)) 
     &       - 2 * (y12(1) + y12(2)) - y12(3) - y12(4)
      c(4,4) = 4 * ( y(1) -  y (2) + y  (3) - y  (4))
     &       + 2 * (y1(1) + y1 (2) - y1 (3) - y1 (4))
     &       + 2 * (y2(1) - y2 (2) - y2 (3) + y2 (4))
     &       +     y12(1) + y12(2) + y12(3) + y12(4)
      end
c
      ! See original implementation in bicubic.f
      M_subroutine 
     &             bcuint1 (y,y1,y2,y12,x1l,x1u,x2l,x2u,
     &                       x1,x2,ansy,ansy1,ansy2,c)
!$acc routine seq
!pgi$r unroll
      implicit none
      integer,parameter::ti_p=t_p
      integer i
      real(t_p) x1,x1l,x1u
      real(t_p) x2,x2l,x2u
      real(t_p) t,u,ansy
      real(t_p) ansy1,ansy2
      real(t_p) y(4),y12(4)
      real(t_p) y1(4),y2(4)
      real(t_p) c(4,4)
 
      !get coefficients, then perform bicubic interpolation
      call bcucof (y,y1,y2,y12,x1u-x1l,x2u-x2l,c)

      t = (x1-x1l) / (x1u-x1l)
      u = (x2-x2l) / (x2u-x2l)
      ansy = 0.0_ti_p
      ansy1 = 0.0_ti_p
      ansy2 = 0.0_ti_p
      do i = 4, 1, -1
       ansy = t*ansy + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
       ansy1 = u*ansy1 +(3.0*c(4,i)*t+2.0*c(3,i))*t + c(2,i)
       ansy2 = t*ansy2 +(3.0*c(i,4)*u+2.0*c(i,3))*u + c(i,2)
      end do
      ansy1 = ansy1 / (x1u-x1l)
      ansy2 = ansy2 / (x2u-x2l)
      end
c
c       torsion-torsion interaction kernel energy term
c
      M_subroutine 
     &   ker_tortor
     &  (k,ia,ib,ic,id,ie,iitortor,ntortorloc,n,loc     ! Input
     &  ,i12,n12,atomic,typeAtm
     &  ,x,y,z,tnx,tny,ttx,tty,tbf,tbx,tby,tbxy,ttorunit,fgrp
     &  ,use_polymer,use_group,use_amd_dih,use_virial
#ifdef  _CUDA_ONLY
     &  ,ett,eDaMD,dettx,detty,dettz            ! Output
#else   
     &  ,ett,dett            ! Output
#endif  
     &  ,vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
     &  ,ver,fea)
      use couple ,only: maxvalue
      use ktrtor ,only: maxntt,maxtgrd,maxtgrd2
      use math   ,only: sqrtpi,radian
!$acc routine
      implicit none
      integer  ,intent(in)::
     &          k,ia,ib,ic,id,ie,iitortor,ntortorloc,n,ver,fea
     &         ,loc(n),i12(maxvalue,n),n12(n),atomic(*)
     &         ,typeAtm(*),tnx(maxntt),tny(maxntt)
      logical  ,intent(in):: use_amd_dih,use_group,use_virial
     &         ,use_polymer
      real(t_p),intent(in)::
     &          ttorunit,fgrp,x(*),y(*),z(*)
     &         ,tbf(maxtgrd2,maxntt),tbxy(maxtgrd2,maxntt)
     &         ,ttx(maxtgrd,maxntt),tty(maxtgrd,maxntt)
     &         ,tbx(maxtgrd2,maxntt),tby(maxtgrd2,maxntt)
#ifdef _CUDA_ONLY
      ener_rtyp,intent(inout):: ett
      real(r_p),intent(inout):: eDaMD
      mdyn_rtyp,intent(inout):: dettx(*),detty(*),dettz(*)
      real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
#else
      real(r_p),intent(inout):: ett
      real(r_p),intent(inout):: dett(*)
      real(r_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
#endif

      integer   pos1,pos2,nlo,xlo,ylo,nhi,nt
      integer   j,m,i1,k1,ia1
      integer   grd,ene,vrl,plm,grp,mpr
      real(t_p) e,sign
      real(t_p) angle1,angle2,value1,value2,cosine1,cosine2
      real(t_p) xia,yia,zia,xib,yib,zib,xic,yic,zic
      real(t_p) xid,yid,zid,xie,yie,zie
      real(t_p) xt,yt,zt,rt2,xu,yu,zu,ru2
      real(t_p) xv,yv,zv,rv2,rtru,rurv
      real(t_p) x1l,x1u,y1l,y1u
      real(t_p) xba,yba,zba,xdc,ydc,zdc,xcb,ycb,zcb,xed,yed,zed
      real(t_p) xia1,yia1,zia1,rcb,rdc
      real(t_p) dedang1,dedang2
      real(t_p) cw(4,4)
      real(t_p) ftt(4),ft12(4),ft1(4),ft2(4)
      real(t_p) vxx2,vyy2,vzz2,vyx2,vzx2,vzy2

      real(t_p) xca,yca,zca,xdb,ydb,zdb,xec,yec,zec
     &         ,dedxia,dedyia,dedzia,dedxib,dedyib,dedzib
     &         ,dedxic,dedyic,dedzic,dedxid,dedyid,dedzid
     &         ,dedxib2,dedyib2,dedzib2,dedxic2,dedyic2,dedzic2
     &         ,dedxid2,dedyid2,dedzid2,dedxie2,dedyie2,dedzie2
      real(t_p) dedxt,dedyt,dedzt,dedxu,dedyu,dedzu,dedxv,dedyv,dedzv
      parameter(
     &    grd=__use_grd__, ene=__use_ene__, vrl=__use_vir__
     &   ,plm=__use_polymer__, grp=__use_groups__, mpr=__use_mpi__
     &         )

      !compute the values of the torsional angles
      xia = x(ia); yia = y(ia); zia = z(ia)
      xib = x(ib); yib = y(ib); zib = z(ib)
      xic = x(ic); yic = y(ic); zic = z(ic)
      xid = x(id); yid = y(id); zid = z(id)
      xie = x(ie); yie = y(ie); zie = z(ie)
      xba = xib - xia; yba = yib - yia; zba = zib - zia
      xcb = xic - xib; ycb = yic - yib; zcb = zic - zib
      xdc = xid - xic; ydc = yid - yic; zdc = zid - zic
      xed = xie - xid; yed = yie - yid; zed = zie - zid

      IF (IAND(fea,plm).ne.0.and.use_polymer) THEN
         call image_inl (xba,yba,zba)
         call image_inl (xcb,ycb,zcb)
         call image_inl (xdc,ydc,zdc)
         call image_inl (xed,yed,zed)
      END IF
      xt   = yba*zcb - ycb*zba
      yt   = zba*xcb - zcb*xba
      zt   = xba*ycb - xcb*yba
      xu   = ycb*zdc - ydc*zcb
      yu   = zcb*xdc - zdc*xcb
      zu   = xcb*ydc - xdc*ycb
      rt2  = xt*xt + yt*yt + zt*zt
      ru2  = xu*xu + yu*yu + zu*zu
      rtru = sqrt(rt2 * ru2)
      xv   = ydc*zed - yed*zdc
      yv   = zdc*xed - zed*xdc
      zv   = xdc*yed - xed*ydc
      rv2  = xv*xv + yv*yv + zv*zv
      rurv = sqrt(ru2 * rv2)
      if (rtru.ne.0.0 .and. rurv.ne.0.0) then
         rcb     = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
         cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
         cosine1 = min(1.0,max(-1.0,cosine1))
         angle1  = radian * acos(cosine1)
         sign    = xba*xu + yba*yu + zba*zu
         if (sign .lt. 0.0)  angle1 = -angle1
         value1  = angle1
         rdc     = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
         cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
         cosine2 = min(1.0,max(-1.0,cosine2))
         angle2  = radian * acos(cosine2)
         sign    = xcb*xv + ycb*yv + zcb*zv
         if (sign .lt. 0.0)  angle2 = -angle2
         value2  = angle2

        !check for inverted chirality at the central atom
c           call chkttor (ib,ic,id,sign,value1,value2,x,y,z,
c    &                    typeAtm,atomic)
c
        !test for chirality at the central torsion-torsion site
         sign = 1.0
         if (n12(ic) .eq. 4) then
            j = 0
!$acc          loop seq
            do i1 = 1, 4
               m = i12(i1,ic)
               if (m.ne.ib.and. m.ne.id) then
                  if (j .eq. 0) then
                     j = m
                  else
                     k1 = m
                  end if
               end if
            end do
            ia1 = 0
            if (typeAtm(j ) .gt. typeAtm(k1))  ia1 = j
            if (typeAtm(k1) .gt. typeAtm(j))  ia1 = k1
            if (atomic(j)  .gt. atomic(k1))  ia1 = j
            if (atomic(k1) .gt. atomic(j))  ia1 = k1
            if (ia1 .ne. 0) then
               xia1 = x(ia1); yia1=y(ia1); zia1=z(ia1);
               call chkttor_(xia1,xib,xic,xid,yia1,yib,yic,yid
     &                      ,zia1,zib,zic,zid,sign,value1,value2)
            end if
         end if
c
c     use bicubic interpolation to compute spline values
c
         nlo = 1
         nhi = tnx(k)
         do while (nhi-nlo .gt. 1)
            nt = (nhi+nlo) / 2
            if (ttx(nt,k) .gt. value1) then
               nhi = nt
            else
               nlo = nt
            end if
         end do
         xlo = nlo
         nlo = 1
         nhi = tny(k)
         do while (nhi-nlo .gt. 1)
            nt = (nhi + nlo)/2
            if (tty(nt,k) .gt. value2) then
               nhi = nt
            else
               nlo = nt
            end if
         end do
         ylo  = nlo
         x1l  = ttx(xlo,k)
         x1u  = ttx(xlo+1,k)
         y1l  = tty(ylo,k)
         y1u  = tty(ylo+1,k)
         pos2 = ylo*tnx(k) + xlo
         pos1 = pos2 - tnx(k)
         ftt(1)  = tbf(pos1,k)
         ftt(2)  = tbf(pos1+1,k)
         ftt(3)  = tbf(pos2+1,k)
         ftt(4)  = tbf(pos2,k)
         ft1(1)  = tbx(pos1,k)
         ft1(2)  = tbx(pos1+1,k)
         ft1(3)  = tbx(pos2+1,k)
         ft1(4)  = tbx(pos2,k)
         ft2(1)  = tby(pos1,k)
         ft2(2)  = tby(pos1+1,k)
         ft2(3)  = tby(pos2+1,k)
         ft2(4)  = tby(pos2,k)
         ft12(1) = tbxy(pos1,k)
         ft12(2) = tbxy(pos1+1,k)
         ft12(3) = tbxy(pos2+1,k)
         ft12(4) = tbxy(pos2,k)
         call bcuint1 (ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u
     &                ,value1,value2,e,dedang1,dedang2
     &                ,cw)
         e = ttorunit * e
         dedang1 = sign * ttorunit * radian * dedang1
         dedang2 = sign * ttorunit * radian * dedang2
c
c     chain rule terms for first angle derivative components
c
         IF (IAND(fea,grp).NE.0.and.use_group) THEN
            e = e * fgrp
            dedang1 = dedang1 * fgrp
            dedang2 = dedang2 * fgrp
         END IF

         IF (IAND(ver,grd).NE.0) THEN
         xca = xic - xia
         yca = yic - yia
         zca = zic - zia
         xdb = xid - xib
         ydb = yid - yib
         zdb = zid - zib
         IF (IAND(fea,plm).ne.0.and.use_polymer) THEN
            call image_inl (xca,yca,zca)
            call image_inl (xdb,ydb,zdb)
         END IF
         dedxt =  dedang1 * (yt*zcb - ycb*zt) / (rt2*rcb)
         dedyt =  dedang1 * (zt*xcb - zcb*xt) / (rt2*rcb)
         dedzt =  dedang1 * (xt*ycb - xcb*yt) / (rt2*rcb)
         dedxu = -dedang1 * (yu*zcb - ycb*zu) / (ru2*rcb)
         dedyu = -dedang1 * (zu*xcb - zcb*xu) / (ru2*rcb)
         dedzu = -dedang1 * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c       compute first derivative components for first angle
c
         dedxia = zcb*dedyt - ycb*dedzt
         dedyia = xcb*dedzt - zcb*dedxt
         dedzia = ycb*dedxt - xcb*dedyt
         dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
         dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
         dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
         dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
         dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
         dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
         dedxid = zcb*dedyu - ycb*dedzu
         dedyid = xcb*dedzu - zcb*dedxu
         dedzid = ycb*dedxu - xcb*dedyu
c
c       chain rule terms for second angle derivative components
c
         xec = xie - xic
         yec = yie - yic
         zec = zie - zic
         IF (IAND(fea,plm).ne.0.and.use_polymer) THEN
            call image_inl (xdb,ydb,zdb)
            call image_inl (xec,yec,zec)
         END IF
         dedxu =  dedang2 * (yu*zdc - ydc*zu) / (ru2*rdc)
         dedyu =  dedang2 * (zu*xdc - zdc*xu) / (ru2*rdc)
         dedzu =  dedang2 * (xu*ydc - xdc*yu) / (ru2*rdc)
         dedxv = -dedang2 * (yv*zdc - ydc*zv) / (rv2*rdc)
         dedyv = -dedang2 * (zv*xdc - zdc*xv) / (rv2*rdc)
         dedzv = -dedang2 * (xv*ydc - xdc*yv) / (rv2*rdc)
c
c       compute first derivative components for second angle
c
         dedxib2 = zdc*dedyu - ydc*dedzu
         dedyib2 = xdc*dedzu - zdc*dedxu
         dedzib2 = ydc*dedxu - xdc*dedyu
         dedxic2 = ydb*dedzu - zdb*dedyu + zed*dedyv - yed*dedzv
         dedyic2 = zdb*dedxu - xdb*dedzu + xed*dedzv - zed*dedxv
         dedzic2 = xdb*dedyu - ydb*dedxu + yed*dedxv - xed*dedyv
         dedxid2 = zcb*dedyu - ycb*dedzu + yec*dedzv - zec*dedyv
         dedyid2 = xcb*dedzu - zcb*dedxu + zec*dedxv - xec*dedzv
         dedzid2 = ycb*dedxu - xcb*dedyu + xec*dedyv - yec*dedxv
         dedxie2 = zdc*dedyv - ydc*dedzv
         dedyie2 = xdc*dedzv - zdc*dedxv
         dedzie2 = ydc*dedxv - xdc*dedyv
         END IF
 
         !increment the torsion-torsion energy and gradient
         IF (IAND(ver,ene).ne.0) THEN
#ifdef _CUDA_ONLY
            ett = ett + tp2enr(e)
#else
            ett = ett + real(e,r_p)
#endif
         END IF

         IF (IAND(ver,grd).ne.0) THEN
         associate(ialoc=>ia,ibloc=>ib,icloc=>ic,idloc=>id,ieloc=>ie)
#ifdef _CUDA_ONLY
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
         ieloc = loc(ie)
         call atomic_add_f1( dettx(ibloc), dedxib+dedxib2 )
         call atomic_add_f1( detty(ibloc), dedyib+dedyib2 )
         call atomic_add_f1( dettz(ibloc), dedzib+dedzib2 )
c
         call atomic_add_f1( dettx(ialoc), dedxia )
         call atomic_add_f1( detty(ialoc), dedyia )
         call atomic_add_f1( dettz(ialoc), dedzia )
c
         call atomic_add_f1( dettx(icloc), dedxic+dedxic2 )
         call atomic_add_f1( detty(icloc), dedyic+dedyic2 )
         call atomic_add_f1( dettz(icloc), dedzic+dedzic2 )
c
         call atomic_add_f1( dettx(idloc), dedxid+dedxid2 )
         call atomic_add_f1( detty(idloc), dedyid+dedyid2 )
         call atomic_add_f1( dettz(idloc), dedzid+dedzid2 )
c
         call atomic_add_f1( dettx(ieloc), dedxie2 )
         call atomic_add_f1( detty(ieloc), dedyie2 )
         call atomic_add_f1( dettz(ieloc), dedzie2 )
#else
         ialoc = 3*(loc(ia)-1)
         ibloc = 3*(loc(ib)-1)
         icloc = 3*(loc(ic)-1)
         idloc = 3*(loc(id)-1)
         ieloc = 3*(loc(ie)-1)
         call atomic_add( dett(1+ibloc), dedxib+dedxib2 )
         call atomic_add( dett(2+ibloc), dedyib+dedyib2 )
         call atomic_add( dett(3+ibloc), dedzib+dedzib2 )
c
         call atomic_add( dett(1+ialoc), dedxia )
         call atomic_add( dett(2+ialoc), dedyia )
         call atomic_add( dett(3+ialoc), dedzia )
c
         call atomic_add( dett(1+icloc), dedxic+dedxic2 )
         call atomic_add( dett(2+icloc), dedyic+dedyic2 )
         call atomic_add( dett(3+icloc), dedzic+dedzic2 )
c
         call atomic_add( dett(1+idloc), dedxid+dedxid2 )
         call atomic_add( dett(2+idloc), dedyid+dedyid2 )
         call atomic_add( dett(3+idloc), dedzid+dedzid2 )
c
         call atomic_add( dett(1+ieloc), dedxie2 )
         call atomic_add( dett(2+ieloc), dedyie2 )
         call atomic_add( dett(3+ieloc), dedzie2 )
#endif
         end associate
         END IF

         !increment the internal virial tensor components
         IF (IAND(ver,vrl).ne.0.and.use_virial) THEN;
           vxx2 = xdc*(dedxid2+dedxie2) -xcb*dedxib2 +xed*dedxie2
           vyx2 = ydc*(dedxid2+dedxie2) -ycb*dedxib2 +yed*dedxie2
           vzx2 = zdc*(dedxid2+dedxie2) -zcb*dedxib2 +zed*dedxie2
           vyy2 = ydc*(dedyid2+dedyie2) -ycb*dedyib2 +yed*dedyie2
           vzy2 = zdc*(dedyid2+dedyie2) -zcb*dedyib2 +zed*dedyie2
           vzz2 = zdc*(dedzid2+dedzie2) -zcb*dedzib2 +zed*dedzie2

           ! Compiler Issue find here - test with next release
           vxx_ =vxx_ +vxx2 +xcb*(dedxic+dedxid) -xba*dedxia +xdc*dedxid
           vxy_ =vxy_ +vyx2 +ycb*(dedxic+dedxid) -yba*dedxia +ydc*dedxid
           vxz_ =vxz_ +vzx2 +zcb*(dedxic+dedxid) -zba*dedxia +zdc*dedxid
           vyy_ =vyy_ +vyy2 +ycb*(dedyic+dedyid) -yba*dedyia +ydc*dedyid
           vyz_ =vyz_ +vzy2 +zcb*(dedyic+dedyid) -zba*dedyia +zdc*dedyid
           vzz_ =vzz_ +vzz2 +zcb*(dedzic+dedzid) -zba*dedzia +zdc*dedzid
         END IF
#ifdef _CUDA_ONLY
         if (use_amd_dih) call atomic_add_m( eDaMD,e )
#endif
      end if
      end subroutine

