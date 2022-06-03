c
c     sorbonne university
c     washington university in saint louis
c     university of texas at austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module     eliaison1                                      ##
c     ##                                                            ##
c     ################################################################
c
c
#define TINKER_CUF
#include  "tinker_precision.h"
#include  "tinker_types.h"
#include  "tinker_cudart.h"
      module eliaisoncu
        use angpot   ,only: OPB_W_D_C,OPB_ALLINGER
     &               ,ANG_IN_PLANE,ANG_HARMONIC,ANG_LINEAR,ANG_FOURIER
        use bndpot   ,only: BND_HARMONIC,BND_MORSE,BND_NO_TYPE
        use couple   ,only: maxvalue
        use cudafor  ,only: atomicAdd
        use ktrtor   ,only: maxntt,maxtgrd,maxtgrd2
        use math     ,only: sqrtpi,radian
        use tinheader,only: ti_p
        use utilcu   ,only: use_virial
#if ~ TINKER_DOUBLE_PREC
     &               ,f_sqrt,f_exp
#endif
        use utilgpu  ,only: RED_BUFF_SIZE
        integer lia_Gs,lia_Bd
        parameter(lia_Bd=128)
        logical,parameter::calc_e=.false.

#include "atomicOp.h.f"

        contains
#include "image.f.inc"
#include "atomicOp.inc.f"

#define __tver__ (__use_grd__+__use_ene__)
#define __tfea__ (__use_mpi__+__use_gamd__)
#define __sufx__

c
        M_subroutine ker_pitors
     &              (ipitors,npitorsloc,ipit,kpit,ptorunit,iuse
     &              ,x,y,z,loc,useAll,use_amd_dih,use_vir
     &              ,ept,eDaMD,deptx,depty,deptz
     &              ,vxx_,vxy_,vxz_,vyy_,vyz_,vzz_)
        implicit none
        integer  ,intent(in)::
     &            ipitors,npitorsloc,ipit(6,*),iuse(*),loc(*)
        logical  ,intent(in):: useAll,use_amd_dih,use_vir
        real(t_p),intent(in):: ptorunit
     &           ,kpit(*),x(*),y(*),z(*)
        ener_rtyp :: ept
        real(r_p) :: eDaMD
        mdyn_rtyp :: deptx(*),depty(*),deptz(*)
        real(t_p) :: vxx_,vxy_,vxz_,vyy_,vyz_,vzz_

        integer ia,ib,ic,id,ie,ig
        real(t_p) e,dedphi
        real(t_p) xt,yt,zt,rt2
        real(t_p) xu,yu,zu,ru2
        real(t_p) xtu,ytu,ztu
        real(t_p) rdc,rtru
        real(t_p) v2,c2,s2
        real(t_p) phi2,dphi2
        real(t_p) sine,cosine
        real(t_p) sine2,cosine2
        real(t_p) xic,yic,zic
        real(t_p) xid,yid,zid
        real(t_p) xip,yip,zip
        real(t_p) xiq,yiq,ziq
        real(t_p) xad,yad,zad
        real(t_p) xbd,ybd,zbd
        real(t_p) xec,yec,zec
        real(t_p) xgc,ygc,zgc
        real(t_p) xcp,ycp,zcp
        real(t_p) xdc,ydc,zdc
        real(t_p) xqd,yqd,zqd
        real(t_p) xdp,ydp,zdp
        real(t_p) xqc,yqc,zqc
        real(t_p) dedxt,dedyt,dedzt
        real(t_p) dedxu,dedyu,dedzu
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        real(t_p) dedxid,dedyid,dedzid
        real(t_p) dedxie,dedyie,dedzie
        real(t_p) dedxig,dedyig,dedzig
        real(t_p) dedxip,dedyip,dedzip
        real(t_p) dedxiq,dedyiq,dedziq
        integer proceed

        ia = ipit(1,ipitors)
        ib = ipit(2,ipitors)
        ic = ipit(3,ipitors)
        id = ipit(4,ipitors)
        ie = ipit(5,ipitors)
        ig = ipit(6,ipitors)

        !decide whether to compute the current interaction
        proceed = merge(1,(iuse(ia).or.iuse(ib).or.iuse(ic)
     &                 .or.iuse(id).or.iuse(ie).or.iuse(ig)), useAll)

        !compute the value of the pi-orbital torsion angle
        if (proceed) then
           block
           real(t_p) xia,yia,zia,xib,yib,zib
           real(t_p) xie,yie,zie,xig,yig,zig
           xia = x(ia)
           yia = y(ia)
           zia = z(ia)
           xib = x(ib)
           yib = y(ib)
           zib = z(ib)
           xic = x(ic)
           yic = y(ic)
           zic = z(ic)
           xid = x(id)
           yid = y(id)
           zid = z(id)
           xie = x(ie)
           yie = y(ie)
           zie = z(ie)
           xig = x(ig)
           yig = y(ig)
           zig = z(ig)
           xad = xia - xid
           yad = yia - yid
           zad = zia - zid
           xbd = xib - xid
           ybd = yib - yid
           zbd = zib - zid
           xec = xie - xic
           yec = yie - yic
           zec = zie - zic
           xgc = xig - xic
           ygc = yig - yic
           zgc = zig - zic
           end block
#if __tfea__ & __use_polymer__
           if (use_polymer) then
              call image_inl (xad,yad,zad)
              call image_inl (xbd,ybd,zbd)
              call image_inl (xec,yec,zec)
              call image_inl (xgc,ygc,zgc)
           end if
#endif
           xip = yad*zbd - ybd*zad + xic
           yip = zad*xbd - zbd*xad + yic
           zip = xad*ybd - xbd*yad + zic
           xiq = yec*zgc - ygc*zec + xid
           yiq = zec*xgc - zgc*xec + yid
           ziq = xec*ygc - xgc*yec + zid
           xcp = xic - xip
           ycp = yic - yip
           zcp = zic - zip
           xdc = xid - xic
           ydc = yid - yic
           zdc = zid - zic
           xqd = xiq - xid
           yqd = yiq - yid
           zqd = ziq - zid
#if __tfea__ & __use_polymer__
           if (use_polymer) then
              call image_inl (xcp,ycp,zcp)
              call image_inl (xdc,ydc,zdc)
              call image_inl (xqd,yqd,zqd)
           end if
#endif
           xt = ycp*zdc - ydc*zcp
           yt = zcp*xdc - zdc*xcp
           zt = xcp*ydc - xdc*ycp
           xu = ydc*zqd - yqd*zdc
           yu = zdc*xqd - zqd*xdc
           zu = xdc*yqd - xqd*ydc
           xtu = yt*zu - yu*zt
           ytu = zt*xu - zu*xt
           ztu = xt*yu - xu*yt
           rt2 = xt*xt + yt*yt + zt*zt
           ru2 = xu*xu + yu*yu + zu*zu
           rtru = sqrt(rt2 * ru2)
           if (rtru .ne. 0.0_ti_p) then
              rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
              cosine = (xt*xu + yt*yu + zt*zu) / rtru
              sine = (xdc*xtu + ydc*ytu + zdc*ztu) / (rdc*rtru)
c
c         set the pi-orbital torsion parameters for this angle
c
              v2 = kpit(ipitors)
              c2 = -1.0_ti_p
              s2 = 0.0_ti_p
c
c         compute the multiple angle trigonometry and the phase terms
c
              cosine2 = cosine*cosine - sine*sine
              sine2 = 2.0_ti_p * cosine * sine
              phi2 = 1.0_ti_p + (cosine2*c2 + sine2*s2)
              dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
c
c         calculate pi-orbital torsion energy and master chain rule term
c
              e = ptorunit * v2 * phi2
              dedphi = ptorunit * v2 * dphi2
c
c         chain rule terms for first derivative components
c
              xdp = xid - xip
              ydp = yid - yip
              zdp = zid - zip
              xqc = xiq - xic
              yqc = yiq - yic
              zqc = ziq - zic
              dedxt = dedphi * (yt*zdc - ydc*zt) / (rt2*rdc)
              dedyt = dedphi * (zt*xdc - zdc*xt) / (rt2*rdc)
              dedzt = dedphi * (xt*ydc - xdc*yt) / (rt2*rdc)
              dedxu = -dedphi * (yu*zdc - ydc*zu) / (ru2*rdc)
              dedyu = -dedphi * (zu*xdc - zdc*xu) / (ru2*rdc)
              dedzu = -dedphi * (xu*ydc - xdc*yu) / (ru2*rdc)
c
c         compute first derivative components for pi-orbital angle
c
              dedxip = zdc*dedyt - ydc*dedzt
              dedyip = xdc*dedzt - zdc*dedxt
              dedzip = ydc*dedxt - xdc*dedyt
              dedxic = ydp*dedzt - zdp*dedyt + zqd*dedyu - yqd*dedzu
              dedyic = zdp*dedxt - xdp*dedzt + xqd*dedzu - zqd*dedxu
              dedzic = xdp*dedyt - ydp*dedxt + yqd*dedxu - xqd*dedyu
              dedxid = zcp*dedyt - ycp*dedzt + yqc*dedzu - zqc*dedyu
              dedyid = xcp*dedzt - zcp*dedxt + zqc*dedxu - xqc*dedzu
              dedzid = ycp*dedxt - xcp*dedyt + xqc*dedyu - yqc*dedxu
              dedxiq = zdc*dedyu - ydc*dedzu
              dedyiq = xdc*dedzu - zdc*dedxu
              dedziq = ydc*dedxu - xdc*dedyu
c
c         compute first derivative components for individual atoms
c
              dedxia = ybd*dedzip - zbd*dedyip
              dedyia = zbd*dedxip - xbd*dedzip
              dedzia = xbd*dedyip - ybd*dedxip
              dedxib = zad*dedyip - yad*dedzip
              dedyib = xad*dedzip - zad*dedxip
              dedzib = yad*dedxip - xad*dedyip
              dedxie = ygc*dedziq - zgc*dedyiq
              dedyie = zgc*dedxiq - xgc*dedziq
              dedzie = xgc*dedyiq - ygc*dedxiq
              dedxig = zec*dedyiq - yec*dedziq
              dedyig = xec*dedziq - zec*dedxiq
              dedzig = yec*dedxiq - xec*dedyiq
              dedxic = dedxic + dedxip - dedxie - dedxig
              dedyic = dedyic + dedyip - dedyie - dedyig
              dedzic = dedzic + dedzip - dedzie - dedzig
              dedxid = dedxid + dedxiq - dedxia - dedxib
              dedyid = dedyid + dedyiq - dedyia - dedyib
              dedzid = dedzid + dedziq - dedzia - dedzib
c
c       increment the total pi-orbital torsion energy and gradient
c
#if __tver__ & __use_ene__
              !call atomic_add_m( ept,e )
              ept = ept + tp2enr(e)
#endif
#if __tver__ & __use_grd__
              associate(ialoc=>ia,ibloc=>ib,icloc=>ic,idloc=>id
     &                 ,ieloc=>ie,igloc=>ig)
              ialoc = loc(ia)
              ibloc = loc(ib)
              icloc = loc(ic)
              idloc = loc(id)
              ieloc = loc(ie)
              igloc = loc(ig)

              call atomic_add_f1( deptx(icloc), dedxic )
              call atomic_add_f1( depty(icloc), dedyic )
              call atomic_add_f1( deptz(icloc), dedzic )
c
              call atomic_add_f1( deptx(ialoc), dedxia )
              call atomic_add_f1( depty(ialoc), dedyia )
              call atomic_add_f1( deptz(ialoc), dedzia )
c
              call atomic_add_f1( deptx(ibloc), dedxib )
              call atomic_add_f1( depty(ibloc), dedyib )
              call atomic_add_f1( deptz(ibloc), dedzib )
c
              call atomic_add_f1( deptx(idloc), dedxid )
              call atomic_add_f1( depty(idloc), dedyid )
              call atomic_add_f1( deptz(idloc), dedzid )
c
              call atomic_add_f1( deptx(ieloc), dedxie )
              call atomic_add_f1( depty(ieloc), dedyie )
              call atomic_add_f1( deptz(ieloc), dedzie )
c
              call atomic_add_f1( deptx(igloc), dedxig )
              call atomic_add_f1( depty(igloc), dedyig )
              call atomic_add_f1( deptz(igloc), dedzig )
              end associate
#endif
#if __tver__ & __use_vir__
              if (use_vir) then; block
                 real(t_p) vxterm,vyterm,vzterm
                 !increment the internal virial tensor components
                 vxterm = dedxid + dedxia + dedxib
                 vyterm = dedyid + dedyia + dedyib
                 vzterm = dedzid + dedzia + dedzib

                 vxx_ = vxx_ + xdc*vxterm + xcp*dedxip - xqd*dedxiq
                 vxy_ = vxy_ + ydc*vxterm + ycp*dedxip - yqd*dedxiq
                 vxz_ = vxz_ + zdc*vxterm + zcp*dedxip - zqd*dedxiq
                 vyy_ = vyy_ + ydc*vyterm + ycp*dedyip - yqd*dedyiq
                 vyz_ = vyz_ + zdc*vyterm + zcp*dedyip - zqd*dedyiq
                 vzz_ = vzz_ + zdc*vzterm + zcp*dedzip - zqd*dedziq
              end block; end if
#endif
#if __tfea__ & __use_gamd__
              if (use_amd_dih) call atomic_add_m( eDaMD,e )
#endif
           end if
        end if
        end subroutine

        ! see alco bicubic.f
        M_subroutine
     &               bcucof (y,y1,y2,y12,d1,d2,c)
!$acc   routine seq
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
     &         + 6 * y1 (1) + 3 * y1 (2) - 3 * y1 (3) - 6 * y1 (4) 
     &         + 6 * y2 (1) - 6 * y2 (2) - 3 * y2 (3) + 3 * y2 (4)
     &         + 4 * y12(1) + 2 * y12(2) +     y12(3) + 2 * y12(4)
        c(3,4) = 6 * (-y(1) +     y(2) -      y(3) +      y(4)) 
     &          -4 * y1(1) - 2 * y1(2) + 2 * y1(3) + 4 * y1(4) 
     &          -3 * y2(1) + 3 * y2(2) + 3 * y2(3) - 3 * y2(4)
     &          -2 *y12(1) -    y12(2) -    y12(3) - 2 *y12(4)

        c(4,1) = 2 * ( y(1) -  y(2)) +  y1(1) +  y1(2)
        c(4,2) = 2 * (y2(1) - y2(2)) + y12(1) + y12(2)
        c(4,3) = 6 * (-y(1) +  y(2)  -   y(3) +   y(4))
     &         + 3 * (y1(3) + y1(4)  -  y1(1) -  y1(2))
     &         + 2 * (2 * (y2(2) - y2(1)) + y2(3) - y2(4)) 
     &         - 2 * (y12(1) + y12(2)) - y12(3) - y12(4)
        c(4,4) = 4 * ( y(1) -  y (2) + y  (3) - y  (4))
     &         + 2 * (y1(1) + y1 (2) - y1 (3) - y1 (4))
     &         + 2 * (y2(1) - y2 (2) - y2 (3) + y2 (4))
     &         +     y12(1) + y12(2) + y12(3) + y12(4)
        end
c
        ! See also bicubic.f
        M_subroutine 
     &               bcuint1 (y,y1,y2,y12,x1l,x1u,x2l,x2u,
     &                         x1,x2,ansy,ansy1,ansy2,c)
!$acc   routine seq
!pgi$r  unroll
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
!$acc   routine(bcucof) seq
 
        !get coefficients, then perform bicubic interpolation
        call bcucof (y,y1,y2,y12,x1u-x1l,x2u-x2l,c)

        t = (x1-x1l) / (x1u-x1l)
        u = (x2-x2l) / (x2u-x2l)
        ansy = 0.0_ti_p
        ansy1 = 0.0_ti_p
        ansy2 = 0.0_ti_p
        do i = 4, 1, -1
         ansy = t*ansy + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
         ansy1 = u*ansy1 +(3.0_ti_p*c(4,i)*t+2.0_ti_p*c(3,i))*t + c(2,i)
         ansy2 = t*ansy2 +(3.0_ti_p*c(i,4)*u+2.0_ti_p*c(i,3))*u + c(i,2)
        end do
        ansy1 = ansy1 / (x1u-x1l)
        ansy2 = ansy2 / (x2u-x2l)
        end

        M_subroutine
     &     chkttor1 (xia,xib,xic,xid,yia,yib,yic,yid,zia,zib,zic,zid
     &              ,sign,value1,value2)
!$acc   routine seq
        implicit none
        real(t_p),intent(in)::xia,xib,xic,xid,yia,yib,yic,yid
     &           ,zia,zib,zic,zid
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
        if (vol .lt. 0.0_ti_p) then
           sign   = -1.0_ti_p
           value1 = -value1
           value2 = -value2
        end if
        end
c
c       torsion-torsion interaction kernel energy term
c
        M_subroutine ker_tortor
     &  (iitortor,ntortorloc,n,itt,ibitor,iuse,loc     ! Input
     &  ,i12,n12,atomic,typeAtm
     &  ,x,y,z,tnx,tny,ttx,tty,tbf,tbx,tby,tbxy,ttorunit
     &  ,use_polymer,use_amd_dih,useAll
     &  ,ett,eDaMD,dettx,detty,dettz            ! Output
     &  ,vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
     &  )
        use couple ,only: maxvalue
!$acc   routine
!$acc   routine(bcuint1)
        implicit none
        integer  ,intent(in)::
     &            iitortor,ntortorloc,n,itt(3,*),ibitor(5,*),iuse(*)
     &           ,loc(n),i12(maxvalue,n),n12(n),atomic(*)
     &           ,typeAtm(*),tnx(maxntt),tny(maxntt)
        logical  ,intent(in):: use_amd_dih,use_polymer,useAll
        real(t_p),intent(in)::
     &            ttorunit,x(*),y(*),z(*)
     &           ,tbf(maxtgrd2,maxntt),tbxy(maxtgrd2,maxntt)
     &           ,ttx(maxtgrd,maxntt),tty(maxtgrd,maxntt)
     &           ,tbx(maxtgrd2,maxntt),tby(maxtgrd2,maxntt)
        ener_rtyp,intent(inout):: ett
        real(r_p),intent(inout):: eDaMD
        mdyn_rtyp,intent(inout):: dettx(*),detty(*),dettz(*)
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_

        integer i,k
        integer pos1,pos2
        integer ia,ib,ic,id,ie
        integer nlo,nhi,nt
        integer xlo,ylo
        integer j,m,i1,k1,ia1
        integer proceed
        real(t_p) e,sign
        real(t_p) angle1,angle2
        real(t_p) value1,value2
        real(t_p) cosine1,cosine2
        real(t_p) xia,yia,zia,xib,yib,zib,xic,yic,zic
        real(t_p) xid,yid,zid,xie,yie,zie
        real(t_p) xt,yt,zt,rt2
        real(t_p) xu,yu,zu,ru2
        real(t_p) xv,yv,zv,rv2
        real(t_p) rtru
        real(t_p) rurv
        real(t_p) x1l,x1u
        real(t_p) y1l,y1u
        real(t_p) xba,yba,zba
        real(t_p) xdc,ydc,zdc
        real(t_p) xcb,ycb,zcb
        real(t_p) xed,yed,zed
        real(t_p) rcb,rdc
        real(t_p) dedang1,dedang2
        real(t_p) cw(4,4)
        real(t_p) ftt(4),ft12(4),ft1(4),ft2(4)

        i        = itt(1,iitortor)
        k        = itt(2,iitortor)
        if (itt(3,iitortor) .eq. 1) then
           ia = ibitor(1,i)
           ib = ibitor(2,i)
           ic = ibitor(3,i)
           id = ibitor(4,i)
           ie = ibitor(5,i)
        else
           ia = ibitor(5,i)
           ib = ibitor(4,i)
           ic = ibitor(3,i)
           id = ibitor(2,i)
           ie = ibitor(1,i)
        end if
 
        !decide whether to compute the current interaction
        proceed = merge(1,(iuse(ia).or.iuse(ib).or.iuse(ic)
     &                 .or.iuse(id).or.iuse(ie))  , useAll)
 
        !compute the values of the torsional angles
        if (proceed) then
           xia = x(ia)
           yia = y(ia)
           zia = z(ia)
           xib = x(ib)
           yib = y(ib)
           zib = z(ib)
           xic = x(ic)
           yic = y(ic)
           zic = z(ic)
           xid = x(id)
           yid = y(id)
           zid = z(id)
           xie = x(ie)
           yie = y(ie)
           zie = z(ie)
           xba = xib - xia
           yba = yib - yia
           zba = zib - zia
           xcb = xic - xib
           ycb = yic - yib
           zcb = zic - zib
           xdc = xid - xic
           ydc = yid - yic
           zdc = zid - zic
           xed = xie - xid
           yed = yie - yid
           zed = zie - zid
#if __tfea__ & __use_polymer__
           if (use_polymer) then
              call image_inl (xba,yba,zba)
              call image_inl (xcb,ycb,zcb)
              call image_inl (xdc,ydc,zdc)
              call image_inl (xed,yed,zed)
           end if
#endif
           xt = yba*zcb - ycb*zba
           yt = zba*xcb - zcb*xba
           zt = xba*ycb - xcb*yba
           xu = ycb*zdc - ydc*zcb
           yu = zcb*xdc - zdc*xcb
           zu = xcb*ydc - xdc*ycb
           rt2 = xt*xt + yt*yt + zt*zt
           ru2 = xu*xu + yu*yu + zu*zu
           rtru = sqrt(rt2 * ru2)
           xv = ydc*zed - yed*zdc
           yv = zdc*xed - zed*xdc
           zv = xdc*yed - xed*ydc
           rv2 = xv*xv + yv*yv + zv*zv
           rurv = sqrt(ru2 * rv2)
           if (rtru.ne.0.0_ti_p .and. rurv.ne.0.0_ti_p) then
              rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
              cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
              cosine1 = min(1.0_ti_p,max(-1.0_ti_p,cosine1))
              angle1 = radian * acos(cosine1)
              sign = xba*xu + yba*yu + zba*zu
              if (sign .lt. 0.0_ti_p)  angle1 = -angle1
              value1 = angle1
              rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
              cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
              cosine2 = min(1.0_ti_p,max(-1.0_ti_p,cosine2))
              angle2 = radian * acos(cosine2)
              sign = xcb*xv + ycb*yv + zcb*zv
              if (sign .lt. 0.0_ti_p)  angle2 = -angle2
              value2 = angle2
 
             !check for inverted chirality at the central atom
c             call chkttor (ib,ic,id,sign,value1,value2,x,y,z,
c    &                      typeAtm,atomic)
c
             !test for chirality at the central torsion-torsion site
              sign = 1.0_ti_p
              if (n12(ic) .eq. 4) then
                 j = 0
!$acc            loop seq
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
                    block
                    real(t_p) xia1,yia1,zia1
                    xia1 = x(ia1); yia1=y(ia1); zia1=z(ia1);
                    call chkttor1(xia1,xib,xic,xid,yia1,yib,yic,yid
     &                           ,zia1,zib,zic,zid,sign,value1,value2)
                    end block
                 end if
              end if
c
c       use bicubic interpolation to compute spline values
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
              ftt(1) = tbf(pos1,k)
              ftt(2) = tbf(pos1+1,k)
              ftt(3) = tbf(pos2+1,k)
              ftt(4) = tbf(pos2,k)
              ft1(1) = tbx(pos1,k)
              ft1(2) = tbx(pos1+1,k)
              ft1(3) = tbx(pos2+1,k)
              ft1(4) = tbx(pos2,k)
              ft2(1) = tby(pos1,k)
              ft2(2) = tby(pos1+1,k)
              ft2(3) = tby(pos2+1,k)
              ft2(4) = tby(pos2,k)
              ft12(1) = tbxy(pos1,k)
              ft12(2) = tbxy(pos1+1,k)
              ft12(3) = tbxy(pos2+1,k)
              ft12(4) = tbxy(pos2,k)
              call bcuint1 (ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u
     &                     ,value1,value2,e,dedang1,dedang2
     &                     ,cw)
              e = ttorunit * e
              dedang1 = sign * ttorunit * radian * dedang1
              dedang2 = sign * ttorunit * radian * dedang2
c
c       chain rule terms for first angle derivative components
c
              block
              real(t_p) xca,yca,zca,xdb,ydb,zdb,xec,yec,zec
     &                 ,dedxia,dedyia,dedzia,dedxib,dedyib,dedzib
     &                 ,dedxic,dedyic,dedzic,dedxid,dedyid,dedzid
     &                 ,dedxib2,dedyib2,dedzib2,dedxic2,dedyic2,dedzic2
     &                 ,dedxid2,dedyid2,dedzid2,dedxie2,dedyie2,dedzie2

              xca = xic - xia
              yca = yic - yia
              zca = zic - zia
              xdb = xid - xib
              ydb = yid - yib
              zdb = zid - zib
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xca,yca,zca)
                 call image_inl (xdb,ydb,zdb)
              end if
#endif
              block
              real(t_p) dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
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
              end block
c
c       chain rule terms for second angle derivative components
c
              xec = xie - xic
              yec = yie - yic
              zec = zie - zic
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xdb,ydb,zdb)
                 call image_inl (xec,yec,zec)
              end if
#endif
              block
              real(t_p) dedxu,dedyu,dedzu,dedxv,dedyv,dedzv
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
              dedxic2 = ydb*dedzu - zdb*dedyu
     &                + zed*dedyv - yed*dedzv
              dedyic2 = zdb*dedxu - xdb*dedzu
     &                + xed*dedzv - zed*dedxv
              dedzic2 = xdb*dedyu - ydb*dedxu
     &                + yed*dedxv - xed*dedyv
              dedxid2 = zcb*dedyu - ycb*dedzu
     &                + yec*dedzv - zec*dedyv
              dedyid2 = xcb*dedzu - zcb*dedxu
     &                + zec*dedxv - xec*dedzv
              dedzid2 = ycb*dedxu - xcb*dedyu
     &                + xec*dedyv - yec*dedxv
              dedxie2 = zdc*dedyv - ydc*dedzv
              dedyie2 = xdc*dedzv - zdc*dedxv
              dedzie2 = ydc*dedxv - xdc*dedyv
              end block
 
              !increment the torsion-torsion energy and gradient
 
#if __tver__ & __use_ene__
              !call atomic_add( ett,e )
              ett = ett + tp2enr(e)
#endif
#if __tver__ & __use_grd__
              associate(ialoc=>ia,ibloc=>ib,icloc=>ic
     &                 ,idloc=>id,ieloc=>ie)
              ialoc = loc(ia)
              ibloc = loc(ib)
              icloc = loc(ic)
              idloc = loc(id)
              ieloc = loc(ie)
c
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
              end associate
#endif
#if __tver__ & __use_vir__
              !increment the internal virial tensor components
              if (use_virial) then; block
              real(t_p) vxx2,vyy2,vzz2,vyx2,vzx2,vzy2
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
              end block; end if
#endif
#if __tfea__ & __use_gamd__
              if (use_amd_dih) call atomic_add( eDaMD,e )
#endif
              end block
           end if
        end if
        end subroutine
c
c       eliaison1_cu gathers all bonded routines in a single kernel
c
        attributes(global) subroutine __Cat(eliaison1_kcu,__sufx__)
     &  ( !bond
     &   bndglob,ibnd,bl,bk,typeAtm,nbondloc,bndunit,cbnd,qbnd,bndtyp_i
     &  ,use_polymer,deb
     &  ,off_bond,num_bond
     &    !urey
     &  ,ureyglob,nureyloc,iury,ul,uk,ureyunit,cury,qury,deub
     &  ,off_urey,num_urey
     &    !angle
     &  ,angleglob,nangleloc,angtypI,iang,anat,ak,afld,angunit
     &  ,cang,qang,pang,sang,dea
     &  ,off_angl,num_angl
     &    !opbend
     &  ,opbendglob,nopbendloc,iopb,opbk,opbunit,opbtypInt
     &  ,copb,qopb,popb,sopb,deopb
     &  ,off_opbd,num_opbd
     &    !strbnd
     &  ,strbndglob,nstrbndloc,isb,sbk,stbnunit,deba
     &  ,off_stbd,num_stbd
     &    !tors
     &  ,torsglob,ntorsloc,itors,tors1,tors2,tors3,tors4,tors5,tors6
     &  ,torsunit,et,det
     &  ,off_tors,num_tors
     &    !pitors
     &  ,pitorsglob,npitorsloc,ipit,kpit,ptorunit,ept,dept
     &  ,off_pitr,num_pitr
     &    !tortor
     &  ,tortorglob,ntortorloc,itt,ibitor,tnx,tny,ttx,tty,tbf,tbx,tby
     &  ,tbxy,i12,n12,atomic,ttorunit,ett,dett
     &  ,off_totr,num_totr
     &    !angtor
     &  ,angtorglob,nangtorloc,iat,kant,atorunit,deat
     &  ,off_antr,num_antr
     &    !improp
     &  ,impropglob,niproploc,iiprop,vprop,kprop,idihunit,eid,deid
     &  ,off_impr,num_impr
     &    !imptor
     &  ,imptorglob,nitorsloc,iitors,itors1,itors2,itors3
     &  ,itorunit,eit,deit
     &  ,off_impt,num_impt
     &    !aMD/GaMD
     &  ,aMDwattype,use_amd_wat1,use_amd_dih,eDaMD,eW1aMD,deW1aMD
     &
     &  ,x_,y_,z_,x,y,z,loc,iuse
     &  ,n,nbloc,grsz,useAll
     &  ,grx,gry,grz,gr_x,gr_y,gr_z,ev_buff,vir_buff
     &  )
        implicit none
        integer  ,value:: n,nbloc,nbondloc,nureyloc,nangleloc,nopbendloc
     &           ,nstrbndloc,ntorsloc,npitorsloc,ntortorloc,nangtorloc
     &           ,niproploc,nitorsloc
     &           ,bndtyp_i,opbtypInt,grsz
     &           ,off_bond,num_bond,off_urey,num_urey,off_angl,num_angl
     &           ,off_tors,num_tors,off_pitr,num_pitr,off_totr,num_totr
     &           ,off_antr,num_antr,off_impr,num_impr,off_impt,num_impt
     &           ,off_opbd,num_opbd,off_stbd,num_stbd
        logical  ,value:: useAll,use_polymer,use_amd_wat1,use_amd_dih
        real(t_p),value:: bndunit,cbnd,qbnd,ureyunit,angunit,opbunit
     &           ,stbnunit,torsunit,ptorunit,ttorunit,atorunit,idihunit
     &           ,itorunit,cang,qang,pang,sang
     &           ,cury,qury,copb,qopb,popb,sopb

        integer  ,device,intent(in),dimension(*)::
     &            bndglob,ureyglob,angleglob,opbendglob,strbndglob
     &           ,torsglob,pitorsglob,tortorglob,angtorglob,impropglob
     &           ,imptorglob,iopb,angtypI,typeAtm,atomic,loc,iuse
        integer  ,device,intent(in)             ::
     &            ibnd(2,*),iury(3,*),iang(4,*),isb(3,*),itors(4,*)
     &           ,ipit(6,*),itt(3,*),ibitor(5,*),iat(4,*),iiprop(4,*)
     &           ,iitors(4,*),tnx(maxntt),tny(maxntt)
     &           ,i12(maxvalue,n),n12(n)
     &           ,aMDwattype(2)
c       logical  ,device,intent(in):: iuse(0:n)
        real(t_p),device,intent(in),dimension(*)::
     &            x,y,z,bl,bk,ul,uk,opbk,anat,ak,afld,kpit,vprop,kprop
        real(t_p),device,intent(in)             ::
     &            sbk(2,*),tors1(4,*),tors2(4,*),tors3(4,*),tors4(4,*)
     &           ,tors5(4,*),tors6(4,*),itors1(4,*),itors2(4,*)
     &           ,itors3(4,*),kant(6,*),tbf(maxtgrd2,maxntt)
     &           ,ttx(maxtgrd,maxntt),tty(maxtgrd,maxntt)
     &           ,tbx(maxtgrd2,maxntt),tby(maxtgrd2,maxntt)
     &           ,tbxy(maxtgrd2,maxntt)
        real(r_p),device,intent(in),dimension(*)::
     &            x_,y_,z_

        real(r_p),device,intent(inout),dimension(*)::
     &            deb,deub,dea,deopb,deba,det,dept,dett,deat,deid,deit
     &           ,deW1aMD
        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
        real(r_p),device :: et,ept,ett,eid,eit,eDaMD,eW1aMD
        ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
        mdyn_rtyp,device,intent(inout),dimension(*)::grx,gry,grz
     &           ,gr_x,gr_y,gr_z

        integer ithread,stride,istat,it
        real(t_p) rstat
        real(r_p) rrstat
        ener_rtyp enrgy
        real(t_p) g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
        !real(t_p),shared:: cw(4,4,lia_Bd)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        !stride  = blockDim%x*gridDim%x
        enrgy   = 0
        g_vxx =0; g_vxy =0; g_vxz =0; g_vyy =0; g_vyz =0; g_vzz =0;

        bond:block
        integer   i,ia,ib,iaglob,ibglob,ibond
     &           ,ia1,ib1,ind,ipe
        real(t_p) e,de,ideal,force,expterm,bde ,dt,dt2,deddt
     &           ,rab,xab,yab,zab
        real(t_p) dedx,dedy,dedz

        !calculate the bond stretch energy and first derivatives
        do ibond = ithread, nbondloc, num_bond
           i     = bndglob(ibond)
#ifdef   USE_NVSHMEM_CUDA
           ipe   = (i-1)/nbond_pe
           ind   = mod((i-1),nbond_pe) +1
           ia    = d_ibnd(ipe)%pel(1,ind)
           ib    = d_ibnd(ipe)%pel(2,ind)
           ideal = d_bl  (ipe)%pel(ind)
           force = d_bk  (ipe)%pel(ind)
#else
           ia    = ibnd(1,i)
           ib    = ibnd(2,i)
           ideal = bl(i)
           force = bk(i)
#endif  
 
           !compute the value of the bond length deviation
           if ( useAll.or.iuse(ia).or.iuse(ib) ) then
              xab = x_(ia) - x_(ib)
              yab = y_(ia) - y_(ib)
              zab = z_(ia) - z_(ib)
#if __tfea__ & __use_polymer__
              if (use_polymer)  call image_inl (xab,yab,zab)
#endif
              rab = f_sqrt(xab*xab + yab*yab + zab*zab)
              dt  = rab - ideal
              if (dt.eq.0.0) cycle
c
c       harmonic potential uses Taylor expansion of Morse potential
c       through the fourth power of the bond length deviation
c
              if (bndtyp_i .eq. BND_HARMONIC) then
                 dt2 = dt * dt
                 e = bndunit * force * dt2 * (1.0_ti_p+cbnd*dt+qbnd*dt2)
                 deddt = 2.0_ti_p * bndunit * force * dt
     &                   * (1.0_ti_p+1.5_ti_p*cbnd*dt+2.0_ti_p*qbnd*dt2)
c
c       Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c       with the approximations alpha = sqrt(ForceConst/BDE) = -2
c       and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
              else if (bndtyp_i .eq. BND_MORSE) then
                 expterm = f_exp(-2.0_ti_p*dt)
                 bde   = 0.25_ti_p * bndunit * force
                 e     = bde * (1.0_ti_p-expterm)**2
                 deddt = 4.0_ti_p * bde * (1.0_ti_p-expterm) * expterm
              end if
 
              !compute chain rule terms needed for derivatives
              if (rab .eq. 0.0_ti_p) then
                 de = 0.0_ti_p
              else
                 de = deddt / rab
              end if
              dedx = de * xab
              dedy = de * yab
              dedz = de * zab
 
#if __tver__ & __use_ene__
             !increment the total Bond energy and first derivatives and virial tensor
              enrgy = enrgy + tp2enr(e) 
#endif
              associate(ialoc=>ia,ibloc=>ib)
#if __tver__ & __use_grd__
              ialoc = loc(ia)
              ibloc = loc(ib)
              call atomic_add_f1( grx(ialoc),dedx )
              call atomic_add_f1( gry(ialoc),dedy )
              call atomic_add_f1( grz(ialoc),dedz )
c
              call atomic_sub_f1( grx(ibloc),dedx )
              call atomic_sub_f1( gry(ibloc),dedy )
              call atomic_sub_f1( grz(ibloc),dedz )
#endif
#if __tver__ & __use_vir__
              if (use_virial) then
                 g_vxx = g_vxx + xab*dedx
                 g_vxy = g_vxy + yab*dedx
                 g_vxz = g_vxz + zab*dedx
                 g_vyy = g_vyy + yab*dedy
                 g_vyz = g_vyz + zab*dedy
                 g_vzz = g_vzz + zab*dedz
              end if
#endif
#if __tfea__ & __use_gamd__
              !aMD storage if waters are considered
              if (use_amd_wat1) then
              if ( typeAtm(ia)== aMDwattype(1)
     $        .or. typeAtm(ib)== aMDwattype(1)) then
                 rrstat = atomicAdd( eW1aMD,real(e,r_p) )
                 it = (ialoc-1)*3
                 call atomic_add_m( deW1aMD(it+1),dedx )
                 call atomic_add_m( deW1aMD(it+2),dedy )
                 call atomic_add_m( deW1aMD(it+3),dedz )
                 it = (ibloc-1)*3
                 call atomic_sub_m( deW1aMD(it+1),dedx )
                 call atomic_sub_m( deW1aMD(it+2),dedy )
                 call atomic_sub_m( deW1aMD(it+3),dedz )
              end if
              end if
#endif
              end associate
           end if
        end do
        end block bond


        if (ithread.gt.off_urey) then; urey:block
        integer i,ia,ic,iurey
        real(t_p) e,de,ideal,force
        real(t_p) dt,dt2,deddt
        real(t_p) dedx,dedy,dedz
        real(t_p) xac,yac,zac,rac
        !logical proceed
 
        !calculate the Urey-Bradley 1-3 energy and first derivatives
        do iurey = ithread, nureyloc+off_urey, num_urey
           i     = ureyglob(iurey-off_urey)
           ia    = iury(1,i)
           ic    = iury(3,i)
           ideal = ul(i)
           force = uk(i)
c
c       compute the value of the 1-3 distance deviation
c
           if ( useAll.or.iuse(ia).or.iuse(ic) ) then
              xac = x(ia) - x(ic)
              yac = y(ia) - y(ic)
              zac = z(ia) - z(ic)
#if __tfea__ & __use_polymer__
              if (use_polymer)  call image_inl (xac,yac,zac)
#endif
              rac = f_sqrt(xac*xac + yac*yac + zac*zac)
              dt  = rac - ideal
              dt2 = dt * dt
              e   = ureyunit * force * dt2 * (1.0_ti_p+cury*dt+qury*dt2)
              deddt = 2.0_ti_p * ureyunit * force * dt
     &             * (1.0_ti_p+1.5_ti_p*cury*dt+2.0_ti_p*qury*dt2)
 
             !compute chain rule terms needed for derivatives
              de    = deddt / rac
              dedx  = de * xac
              dedy  = de * yac
              dedz  = de * zac
 
#if __tver__ & __use_ene__
             !increment the total Urey-Bradley energy and first derivatives and virial tensor
              enrgy = enrgy + tp2enr(e)
#endif
#if __tver__ & __use_grd__
              associate(ialoc=>ia,icloc=>ic)
              ialoc  = loc(ia)
              icloc  = loc(ic)
              call atomic_add_f1( grx(ialoc),dedx )
              call atomic_add_f1( gry(ialoc),dedy )
              call atomic_add_f1( grz(ialoc),dedz )
c
              call atomic_sub_f1( grx(icloc),dedx )
              call atomic_sub_f1( gry(icloc),dedy )
              call atomic_sub_f1( grz(icloc),dedz )
              end associate
#endif
#if __tver__ & __use_vir__
              if (use_virial) then
                 g_vxx = g_vxx + xac*dedx
                 g_vxy = g_vxy + yac*dedx
                 g_vxz = g_vxz + zac*dedx
                 g_vyy = g_vyy + yac*dedy
                 g_vyz = g_vyz + zac*dedy
                 g_vzz = g_vzz + zac*dedz
              end if
#endif
           end if
        end do
        end block urey; end if


        if(ithread.gt.off_opbd) then; opbend:block
        integer i,iopbend,iopbendloc
        integer ia,ib,ic,id
#ifdef   USE_NVSHMEM_CUDA
        integer ipe,ind
#endif  
        real(t_p) e,angle1,force
        real(t_p) dot,cosine
        real(t_p) cc,ee,bkk2,term
        real(t_p) deddt,dedcos
        real(t_p) dt,dt2,dt3,dt4
        real(t_p) xab,yab,zab
        real(t_p) xcb,ycb,zcb
        real(t_p) xdb,ydb,zdb
        real(t_p) xad,yad,zad
        real(t_p) xcd,ycd,zcd
        real(t_p) rdb2,rad2,rcd2
        real(t_p) rab2,rcb2
        real(t_p) dccdxia,dccdyia,dccdzia
        real(t_p) dccdxic,dccdyic,dccdzic
        real(t_p) dccdxid,dccdyid,dccdzid
        real(t_p) deedxia,deedyia,deedzia
        real(t_p) deedxic,deedyic,deedzic
        real(t_p) deedxid,deedyid,deedzid
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        real(t_p) dedxid,dedyid,dedzid
        integer proceed
c
c       calculate the out-of-plane bending energy and derivatives
c
        do iopbendloc = ithread, nopbendloc+off_opbd, num_opbd
           iopbend = opbendglob(iopbendloc-off_opbd)
           i  = iopb(iopbend)
#ifdef   USE_NVSHMEM_CUDA
           ipe=     (i-1)/nangle_pe
           ind= mod((i-1),nangle_pe) +1
           ia = d_iang(ipe)%pel(1,ind)
           ib = d_iang(ipe)%pel(2,ind)
           ic = d_iang(ipe)%pel(3,ind)
           id = d_iang(ipe)%pel(4,ind)
#else
           ia = iang(1,i)
           ib = iang(2,i)
           ic = iang(3,i)
           id = iang(4,i)
#endif  
           force = opbk(iopbend)
 
          !decide whether to compute the current interaction
           proceed = merge(1,iuse(ia).or.iuse(ib).or.iuse(ic)
     &                   .or.iuse(id), useAll)
 
          !get the coordinates of the atoms at trigonal center
           if (proceed) then
              block
              real(t_p) xia,yia,zia,xib,yib,zib
              real(t_p) xic,yic,zic,xid,yid,zid
              xia = x(ia)
              yia = y(ia)
              zia = z(ia)
              xib = x(ib)
              yib = y(ib)
              zib = z(ib)
              xic = x(ic)
              yic = y(ic)
              zic = z(ic)
              xid = x(id)
              yid = y(id)
              zid = z(id)
c
c       compute the out-of-plane bending angle
c
              xab = xia - xib
              yab = yia - yib
              zab = zia - zib
              xcb = xic - xib
              ycb = yic - yib
              zcb = zic - zib
              xdb = xid - xib
              ydb = yid - yib
              zdb = zid - zib
              xad = xia - xid
              yad = yia - yid
              zad = zia - zid
              xcd = xic - xid
              ycd = yic - yid
              zcd = zic - zid
              end block
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xab,yab,zab)
                 call image_inl (xcb,ycb,zcb)
                 call image_inl (xdb,ydb,zdb)
                 call image_inl (xad,yad,zad)
                 call image_inl (xcd,ycd,zcd)
              end if
#endif
c
c       W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
              if (opbtypInt .eq. OPB_W_D_C ) then
                 rab2 = xab*xab + yab*yab + zab*zab
                 rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                 dot = xab*xcb+yab*ycb+zab*zcb
                 cc = rab2*rcb2 - dot*dot
c
c       Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
              else if (opbtypInt .eq. OPB_ALLINGER) then
                 rad2 = xad*xad + yad*yad + zad*zad
                 rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
                 dot = xad*xcd + yad*ycd + zad*zcd
                 cc = rad2*rcd2 - dot*dot
              end if
c
c       find the out-of-plane angle bending energy
c
              ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)
     &                + zdb*(xab*ycb-yab*xcb)
              rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
              if (rdb2.ne.0.0_ti_p .and. cc.ne.0.0_ti_p) then
                 bkk2 = rdb2 - ee*ee/cc
                 cosine = f_sqrt(bkk2/rdb2)
                 cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
                 angle1 = radian * acos(cosine)
                 dt = angle1
                 dt2 = dt * dt
                 dt3 = dt2 * dt
                 dt4 = dt2 * dt2
                 e = opbunit * force * dt2
     &                * (1.0_ti_p+copb*dt+qopb*dt2+popb*dt3+sopb*dt4)
                 deddt = opbunit * force * dt * radian
     &                * (2.0_ti_p + 3.0_ti_p*copb*dt + 4.0_ti_p*qopb*dt2
     &                 + 5.0_ti_p*popb*dt3 + 6.0_ti_p*sopb*dt4)
                 dedcos = -deddt * sign(1.0_ti_p,ee) / f_sqrt(cc*bkk2)
c
c       chain rule terms for first derivative components
c
                 if (opbtypInt .eq. OPB_W_D_C) then
                    term = ee / cc
                    dccdxia = (xab*rcb2-xcb*dot) * term
                    dccdyia = (yab*rcb2-ycb*dot) * term
                    dccdzia = (zab*rcb2-zcb*dot) * term
                    dccdxic = (xcb*rab2-xab*dot) * term
                    dccdyic = (ycb*rab2-yab*dot) * term
                    dccdzic = (zcb*rab2-zab*dot) * term
                    dccdxid = 0.0_ti_p
                    dccdyid = 0.0_ti_p
                    dccdzid = 0.0_ti_p
                 else if (opbtypInt .eq. OPB_ALLINGER) then
                    term = ee / cc
                    dccdxia = (xad*rcd2-xcd*dot) * term
                    dccdyia = (yad*rcd2-ycd*dot) * term
                    dccdzia = (zad*rcd2-zcd*dot) * term
                    dccdxic = (xcd*rad2-xad*dot) * term
                    dccdyic = (ycd*rad2-yad*dot) * term
                    dccdzic = (zcd*rad2-zad*dot) * term
                    dccdxid = -dccdxia - dccdxic
                    dccdyid = -dccdyia - dccdyic
                    dccdzid = -dccdzia - dccdzic
                 end if
                 term = ee / rdb2
                 deedxia = ydb*zcb - zdb*ycb
                 deedyia = zdb*xcb - xdb*zcb
                 deedzia = xdb*ycb - ydb*xcb
                 deedxic = yab*zdb - zab*ydb
                 deedyic = zab*xdb - xab*zdb
                 deedzic = xab*ydb - yab*xdb
                 deedxid = ycb*zab - zcb*yab + xdb*term
                 deedyid = zcb*xab - xcb*zab + ydb*term
                 deedzid = xcb*yab - ycb*xab + zdb*term
c
c       compute first derivative components for this angle
c
                 dedxia = dedcos * (dccdxia+deedxia)
                 dedyia = dedcos * (dccdyia+deedyia)
                 dedzia = dedcos * (dccdzia+deedzia)
                 dedxic = dedcos * (dccdxic+deedxic)
                 dedyic = dedcos * (dccdyic+deedyic)
                 dedzic = dedcos * (dccdzic+deedzic)
                 dedxid = dedcos * (dccdxid+deedxid)
                 dedyid = dedcos * (dccdyid+deedyid)
                 dedzid = dedcos * (dccdzid+deedzid)
                 dedxib = -dedxia - dedxic - dedxid
                 dedyib = -dedyia - dedyic - dedyid
                 dedzib = -dedzia - dedzic - dedzid
c
c       increment the out-of-plane bending energy and gradient and virial tensor
c
#if __tver__ & __use_ene__
                 !istat = atomicAdd( eopb, e )
                 enrgy = enrgy + tp2enr(e)
#endif
#if __tver__ & __use_grd__
                 associate(ialoc=>ia,ibloc=>ib,icloc=>ic,idloc=>id)
                 ialoc = loc(ia)
                 ibloc = loc(ib)
                 icloc = loc(ic)
                 idloc = loc(id)
                 call atomic_add_f1( grx(ibloc),dedxib )
                 call atomic_add_f1( gry(ibloc),dedyib )
                 call atomic_add_f1( grz(ibloc),dedzib )
c
                 call atomic_add_f1( grx(ialoc),dedxia )
                 call atomic_add_f1( gry(ialoc),dedyia )
                 call atomic_add_f1( grz(ialoc),dedzia )
c
                 call atomic_add_f1( grx(icloc),dedxic )
                 call atomic_add_f1( gry(icloc),dedyic )
                 call atomic_add_f1( grz(icloc),dedzic )
c
                 call atomic_add_f1( grx(idloc),dedxid )
                 call atomic_add_f1( gry(idloc),dedyid )
                 call atomic_add_f1( grz(idloc),dedzid )
                 end associate
#endif
#if __tver__ & __use_vir__
                 if (use_virial) then
                 g_vxx = g_vxx + xab*dedxia + xcb*dedxic + xdb*dedxid
                 g_vxy = g_vxy + yab*dedxia + ycb*dedxic + ydb*dedxid
                 g_vxz = g_vxz + zab*dedxia + zcb*dedxic + zdb*dedxid
                 g_vyy = g_vyy + yab*dedyia + ycb*dedyic + ydb*dedyid
                 g_vyz = g_vyz + zab*dedyia + zcb*dedyic + zdb*dedyid
                 g_vzz = g_vzz + zab*dedzia + zcb*dedzic + zdb*dedzid
                 end if
#endif
              end if
           end if
        end do
        end block opbend; end if


        if(ithread.gt.off_stbd) then; strbend:block
        integer i,j,k,istrbnd,istrbndloc
        integer ia,ib,ic
        integer ialoc,ibloc,icloc
#ifdef USE_NVSHMEM_CUDA
        integer ipe,ind
#endif
        real(t_p) e,dr1,dr2,dt
        real(t_p) angle1
        real(t_p) force1,force2
        real(t_p) dot,cosine
        real(t_p) xab,yab,zab
        real(t_p) xcb,ycb,zcb
        real(t_p) rab,rab2
        real(t_p) rcb,rcb2
        real(t_p) xp,yp,zp,rp
        real(t_p) term1,term2
        real(t_p) termr,term1t,term2t
        real(t_p) ddtdxia,ddtdyia,ddtdzia
        real(t_p) ddtdxic,ddtdyic,ddtdzic
        real(t_p) ddrdxia,ddrdyia,ddrdzia
        real(t_p) ddrdxic,ddrdyic,ddrdzic
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        integer proceed
c
c       calculate the stretch-bend energy and first derivatives
c
        do istrbndloc = ithread, nstrbndloc+off_stbd, num_stbd
           istrbnd = strbndglob(istrbndloc-off_stbd)
           i       = isb(1,istrbnd)
#ifdef   USE_NVSHMEM_CUDA
           ipe     =     (i-1)/nangle_pe
           ind     = mod((i-1),nangle_pe) +1
           ia      = d_iang(ipe)%pel(1,ind)
           ib      = d_iang(ipe)%pel(2,ind)
           ic      = d_iang(ipe)%pel(3,ind)
#else
           ia     = iang(1,i)
           ib     = iang(2,i)
           ic     = iang(3,i)
#endif  
           force1 = sbk(1,istrbnd)
           force2 = sbk(2,istrbnd)
c
c       get the coordinates of the atoms in the angle
c
           if (useAll.or.iuse(ia).or.iuse(ib).or.iuse(ic)) then
              block
              real(t_p) xia,yia,zia,xib,yib,zib,xic,yic,zic
              xia = x(ia)
              yia = y(ia)
              zia = z(ia)
              xib = x(ib)
              yib = y(ib)
              zib = z(ib)
              xic = x(ic)
              yic = y(ic)
              zic = z(ic)
c
c       compute the value of the bond angle
c
              xab = xia - xib
              yab = yia - yib
              zab = zia - zib
              xcb = xic - xib
              ycb = yic - yib
              zcb = zic - zib
              end block
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xab,yab,zab)
                 call image_inl (xcb,ycb,zcb)
              end if
#endif
              rab2 = xab*xab + yab*yab + zab*zab
              rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
              if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
                 rab = sqrt(rab2)
                 rcb = sqrt(rcb2)
                 xp = ycb*zab - zcb*yab
                 yp = zcb*xab - xcb*zab
                 zp = xcb*yab - ycb*xab
                 rp = sqrt(xp*xp + yp*yp + zp*zp)
                 rp = max(rp,0.001_ti_p)
                 dot = xab*xcb + yab*ycb + zab*zcb
                 cosine = dot / (rab*rcb)
                 cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
                 angle1 = radian * acos(cosine)
c
c       find chain rule terms for the bond angle deviation
c
                 dt    = angle1 - anat(i)
                 term1 = -radian / (rab2*rp)
                 term2 = radian / (rcb2*rp)
                 ddtdxia = term1 * (yab*zp-zab*yp)
                 ddtdyia = term1 * (zab*xp-xab*zp)
                 ddtdzia = term1 * (xab*yp-yab*xp)
                 ddtdxic = term2 * (ycb*zp-zcb*yp)
                 ddtdyic = term2 * (zcb*xp-xcb*zp)
                 ddtdzic = term2 * (xcb*yp-ycb*xp)
c
c       find chain rule terms for the bond length deviations
c
                 j = isb(2,istrbnd)
                 k = isb(3,istrbnd)
#ifdef   USE_NVSHMEM_CUDA
                 ipe = (j-1)/nbond_pe
                 ind = mod((j-1),nbond_pe) +1
                 dr1 = rab - d_bl(ipe)%pel(ind)
                 ipe = (k-1)/nbond_pe
                 ind = mod((k-1),nbond_pe) +1
                 dr2 = rcb - d_bl(ipe)%pel(ind)
#else
                 dr1 = rab - bl(j)
                 dr2 = rcb - bl(k)
#endif  
                 term1 = 1.0_ti_p / rab
                 term2 = 1.0_ti_p / rcb
                 ddrdxia = term1 * xab
                 ddrdyia = term1 * yab
                 ddrdzia = term1 * zab
                 ddrdxic = term2 * xcb
                 ddrdyic = term2 * ycb
                 ddrdzic = term2 * zcb
c
c       abbreviations used in defining chain rule terms
c
                 term1  = stbnunit * force1
                 term2  = stbnunit * force2
                 termr  = term1*dr1 + term2*dr2
                 term1t = term1 * dt
                 term2t = term2 * dt
c
c       get the energy and master chain rule terms for derivatives
c
                 e = termr * dt
                 dedxia = term1t*ddrdxia + termr*ddtdxia
                 dedyia = term1t*ddrdyia + termr*ddtdyia
                 dedzia = term1t*ddrdzia + termr*ddtdzia
                 dedxic = term2t*ddrdxic + termr*ddtdxic
                 dedyic = term2t*ddrdyic + termr*ddtdyic
                 dedzic = term2t*ddrdzic + termr*ddtdzic
                 dedxib = -dedxia - dedxic
                 dedyib = -dedyia - dedyic
                 dedzib = -dedzia - dedzic
 
#if __tver__ & __use_ene__
                !increment the total stretch-bend energy and derivatives
                 !stat = atomicAdd( eba, e )
                 enrgy = enrgy + tp2enr(e)
#endif
                 associate(ialoc=>ia,ibloc=>ib,icloc=>ic)
#if __tver__ & __use_grd__
                 ialoc  = loc(ia)
                 ibloc  = loc(ib)
                 icloc  = loc(ic)
                 call atomic_add_f1( grx(ialoc),dedxia )
                 call atomic_add_f1( gry(ialoc),dedyia )
                 call atomic_add_f1( grz(ialoc),dedzia )

                 call atomic_add_f1( grx(ibloc),dedxib )
                 call atomic_add_f1( gry(ibloc),dedyib )
                 call atomic_add_f1( grz(ibloc),dedzib )

                 call atomic_add_f1( grx(icloc),dedxic )
                 call atomic_add_f1( gry(icloc),dedyic )
                 call atomic_add_f1( grz(icloc),dedzic )
#endif
#if __tfea__ & __use_gamd__
                !aMD storage if waters are considered
                 if (use_amd_wat1) then
                 if (typeAtm(ia)==aMDwattype(1) 
     &          .or. typeAtm(ib)==aMDwattype(1)
     &          .or. typeAtm(ic)==aMDwattype(1)) then
                    rrstat = atomicAdd( eW1aMD, real(e,r_p) )
                    it = (ialoc-1)*3
                    call atomic_add_m( deW1aMD(it+1),dedxia )
                    call atomic_add_m( deW1aMD(it+2),dedyia )
                    call atomic_add_m( deW1aMD(it+3),dedzia )
                    it = (ibloc-1)*3
                    call atomic_add_m( deW1aMD(it+1),dedxib )
                    call atomic_add_m( deW1aMD(it+2),dedyib )
                    call atomic_add_m( deW1aMD(it+3),dedzib )
                    it = (icloc-1)*3
                    call atomic_add_m( deW1aMD(it+1),dedxic )
                    call atomic_add_m( deW1aMD(it+2),dedyic )
                    call atomic_add_m( deW1aMD(it+3),dedzic )
                 end if
                 end if
#endif
                 end associate
#if __tver__ & __use_vir__
                 if (use_virial) then
                   !increment the internal virial tensor components
                    g_vxx = g_vxx + xab*dedxia + xcb*dedxic
                    g_vxy = g_vxy + yab*dedxia + ycb*dedxic
                    g_vxz = g_vxz + zab*dedxia + zcb*dedxic
                    g_vyy = g_vyy + yab*dedyia + ycb*dedyic
                    g_vyz = g_vyz + zab*dedyia + zcb*dedyic
                    g_vzz = g_vzz + zab*dedzia + zcb*dedzic
                 end if
#endif
              end if
           end if
        end do
        end block strbend; end if


        if (ithread.gt.off_angl) then; angle: block
        integer i,ia,ib,ic,id
#ifdef USE_NVSHMEM_CUDA
        integer ipe,ind
#endif  
        integer iangle
        real(t_p) e,ideal,force
        real(t_p) fold,factor,dot
        real(t_p) cosine,sine
        real(t_p) angle1
        real(t_p) dt,dt2,dt3,dt4
        real(t_p) deddt,term
        real(t_p) terma,termc
        real(t_p) xia,yia,zia
        real(t_p) xib,yib,zib
        real(t_p) xic,yic,zic
        real(t_p) xid,yid,zid
        real(t_p) xab,yab,zab
        real(t_p) xcb,ycb,zcb
        real(t_p) xp,yp,zp,rp
        real(t_p) xad,yad,zad
        real(t_p) xbd,ybd,zbd
        real(t_p) xcd,ycd,zcd
        real(t_p) xip,yip,zip
        real(t_p) xap,yap,zap
        real(t_p) xcp,ycp,zcp
        real(t_p) rab2,rcb2
        real(t_p) rap2,rcp2
        real(t_p) xt,yt,zt
        real(t_p) rt2,ptrt2
        real(t_p) xm,ym,zm,rm
        real(t_p) delta,delta2
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        real(t_p) dedxid,dedyid,dedzid
        real(t_p) dedxip,dedyip,dedzip
        real(t_p) dpdxia,dpdyia,dpdzia
        real(t_p) dpdxic,dpdyic,dpdzic
        integer proceed
c
c       calculate the bond angle bending energy term
c
        do iangle = ithread, nangleloc+off_angl, num_angl
           i     = angleglob(iangle-off_angl)
#ifdef   USE_NVSHMEM_CUDA
           ipe   =     (i-1)/nangle_pe
           ind   = mod((i-1),nangle_pe) +1
           ia    = d_iang(ipe)%pel(1,ind)
           ib    = d_iang(ipe)%pel(2,ind)
           ic    = d_iang(ipe)%pel(3,ind)
           id    = d_iang(ipe)%pel(4,ind)
#else
           ia = iang(1,i)
           ib = iang(2,i)
           ic = iang(3,i)
           id = iang(4,i)
#endif  
           ideal = anat(i)
           force = ak(i)
 
          !decide whether to compute the current interaction
           if (angtypI(i) .eq. ANG_IN_PLANE) then
              proceed = merge
     &        (1,iuse(ia).or.iuse(ib).or.iuse(ic).or.iuse(id),useAll)
           else
              proceed = merge(1,iuse(ia).or.iuse(ib).or.iuse(ic),useAll)
           end if
 
          !get the coordinates of the atoms in the angle
           if (proceed) then  !Proceed
              xia = x(ia)
              yia = y(ia)
              zia = z(ia)
              xib = x(ib)
              yib = y(ib)
              zib = z(ib)
              xic = x(ic)
              yic = y(ic)
              zic = z(ic)
 
             !compute the bond angle bending energy and gradient
              if (angtypI(i) .ne. ANG_IN_PLANE) then
                 xab = xia - xib
                 yab = yia - yib
                 zab = zia - zib
                 xcb = xic - xib
                 ycb = yic - yib
                 zcb = zic - zib
#if __tfea__ & __use_polymer__
                 if (use_polymer) then
                    call image_inl (xab,yab,zab)
                    call image_inl (xcb,ycb,zcb)
                 end if
#endif
                 rab2 = xab*xab + yab*yab + zab*zab
                 rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                 if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
                    xp  = ycb*zab - zcb*yab
                    yp  = zcb*xab - xcb*zab
                    zp  = xcb*yab - ycb*xab
                    rp  = sqrt(xp*xp + yp*yp + zp*zp)
                    rp  = max(rp,0.000001_ti_p)
                    dot = xab*xcb + yab*ycb + zab*zcb
                    cosine = dot / sqrt(rab2*rcb2)
                    cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
                    angle1 = radian * acos(cosine)
c
c       get the energy and master chain rule term for derivatives
c
                    if (angtypI(i) .eq. ANG_HARMONIC) then
                       dt  = angle1 - ideal
                       dt2 = dt * dt
                       dt3 = dt2 * dt
                       dt4 = dt2 * dt2
                       e   = angunit * force * dt2
     &                   * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                     deddt = angunit * force * dt * radian
     &              * (2.0_ti_p + 3.0_ti_p*cang*dt + 4.0_ti_p*qang*dt2
     &                          + 5.0_ti_p*pang*dt3 + 6.0_ti_p*sang*dt4)
                    else if (angtypI(i) .eq. ANG_LINEAR) then
                       factor = 2.0_ti_p * angunit * radian**2
                       sine   = sqrt(1.0_ti_p-cosine*cosine)
                       e      = factor * force * (1.0_ti_p+cosine)
                       deddt  = -factor * force * sine
                    else if (angtypI(i) .eq. ANG_FOURIER) then
                       fold   = afld(i)
                       factor = 2.0_ti_p * angunit * (radian/fold)**2
                       cosine = cos((fold*angle1-ideal)/radian)
                       sine   = sin((fold*angle1-ideal)/radian)
                       e      =  factor * force * (1.0_ti_p+cosine)
                       deddt  = -factor * force * fold * sine
                    end if
c
c       compute derivative components for this interaction
c
                    terma = -deddt / (rab2*rp)
                    termc =  deddt / (rcb2*rp)
                    dedxia = terma * (yab*zp-zab*yp)
                    dedyia = terma * (zab*xp-xab*zp)
                    dedzia = terma * (xab*yp-yab*xp)
                    dedxic = termc * (ycb*zp-zcb*yp)
                    dedyic = termc * (zcb*xp-xcb*zp)
                    dedzic = termc * (xcb*yp-ycb*xp)
                    dedxib = -dedxia - dedxic
                    dedyib = -dedyia - dedyic
                    dedzib = -dedzia - dedzic
c
c       increment the total bond angle energy and derivatives
c
#if __tver__ & __use_ene__
                    !stat= atomicAdd( ea,e )
                    enrgy = enrgy + tp2enr(e)
#endif
#if __tver__ & __use_grd__
                    associate(ialoc=>ia,ibloc=>ib,icloc=>ic)
                    ialoc = loc(ia)
                    ibloc = loc(ib)
                    icloc = loc(ic)
                    call atomic_add_f1( grx(ialoc), dedxia )
                    call atomic_add_f1( gry(ialoc), dedyia )
                    call atomic_add_f1( grz(ialoc), dedzia )
c
                    call atomic_add_f1( grx(ibloc), dedxib )
                    call atomic_add_f1( gry(ibloc), dedyib )
                    call atomic_add_f1( grz(ibloc), dedzib )
c
                    call atomic_add_f1( grx(icloc), dedxic )
                    call atomic_add_f1( gry(icloc), dedyic )
                    call atomic_add_f1( grz(icloc), dedzic )
                    end associate
#endif
#if __tver__ & __use_vir__
                   !increment the internal virial tensor components
                    if (use_virial) then
                       g_vxx = g_vxx + xab*dedxia + xcb*dedxic
                       g_vxy = g_vxy + yab*dedxia + ycb*dedxic
                       g_vxz = g_vxz + zab*dedxia + zcb*dedxic
                       g_vyy = g_vyy + yab*dedyia + ycb*dedyic
                       g_vyz = g_vyz + zab*dedyia + zcb*dedyic
                       g_vzz = g_vzz + zab*dedzia + zcb*dedzic
                    end if
#endif
                 end if
c
c       compute the projected in-plane angle energy and gradient
c
              else
                 xid = x(id)
                 yid = y(id)
                 zid = z(id)
                 xad = xia - xid
                 yad = yia - yid
                 zad = zia - zid
                 xbd = xib - xid
                 ybd = yib - yid
                 zbd = zib - zid
                 xcd = xic - xid
                 ycd = yic - yid
                 zcd = zic - zid
#if __tfea__ & __use_polymer__
                 if (use_polymer) then
                    call image_inl (xad,yad,zad)
                    call image_inl (xbd,ybd,zbd)
                    call image_inl (xcd,ycd,zcd)
                 end if
#endif
                 xt = yad*zcd - zad*ycd
                 yt = zad*xcd - xad*zcd
                 zt = xad*ycd - yad*xcd
                 rt2 = xt*xt + yt*yt + zt*zt
                 delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
                 xip = xib + xt*delta
                 yip = yib + yt*delta
                 zip = zib + zt*delta
                 xap = xia - xip
                 yap = yia - yip
                 zap = zia - zip
                 xcp = xic - xip
                 ycp = yic - yip
                 zcp = zic - zip
#if __tfea__ & __use_polymer__
                 if (use_polymer) then
                    call image_inl (xap,yap,zap)
                    call image_inl (xcp,ycp,zcp)
                 end if
#endif
                 rap2 = xap*xap + yap*yap + zap*zap
                 rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
                 if (rap2.ne.0.0_ti_p .and. rcp2.ne.0.0_ti_p) then
                    xm = ycp*zap - zcp*yap
                    ym = zcp*xap - xcp*zap
                    zm = xcp*yap - ycp*xap
                    rm = sqrt(xm*xm + ym*ym + zm*zm)
                    rm = max(rm,0.000001_ti_p)
                    dot = xap*xcp + yap*ycp + zap*zcp
                    cosine = dot / sqrt(rap2*rcp2)
                    cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
                    angle1 = radian * acos(cosine)
c
c       get the energy and master chain rule term for derivatives
c
                    dt  = angle1 - ideal
                    dt2 = dt * dt
                    dt3 = dt2 * dt
                    dt4 = dt2 * dt2
                    e   = angunit * force * dt2
     &                * (1.0_ti_p+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                 deddt  = angunit * force * dt * radian
     &               * (2.0_ti_p + 3.0_ti_p*cang*dt + 4.0_ti_p*qang*dt2
     &                          + 5.0_ti_p*pang*dt3 + 6.0_ti_p*sang*dt4)
c
c       chain rule terms for first derivative components
c
                    terma = -deddt / (rap2*rm)
                    termc = deddt / (rcp2*rm)
                    dedxia = terma * (yap*zm-zap*ym)
                    dedyia = terma * (zap*xm-xap*zm)
                    dedzia = terma * (xap*ym-yap*xm)
                    dedxic = termc * (ycp*zm-zcp*ym)
                    dedyic = termc * (zcp*xm-xcp*zm)
                    dedzic = termc * (xcp*ym-ycp*xm)
                    dedxip = -dedxia - dedxic
                    dedyip = -dedyia - dedyic
                    dedzip = -dedzia - dedzic
c
c       chain rule components for the projection of the central atom
c
                    delta2 = 2.0_ti_p * delta
                    ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
                    term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
                    dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
                    term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
                    dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
                    term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
                    dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
                    term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
                    dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
                    term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
                    dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
                    term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
                    dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2
c
c       compute derivative components for this interaction
c
                    dedxia =  dedxia + dpdxia
                    dedyia =  dedyia + dpdyia
                    dedzia =  dedzia + dpdzia
                    dedxib =  dedxip
                    dedyib =  dedyip
                    dedzib =  dedzip
                    dedxic =  dedxic + dpdxic
                    dedyic =  dedyic + dpdyic
                    dedzic =  dedzic + dpdzic
                    dedxid = -dedxia - dedxib - dedxic
                    dedyid = -dedyia - dedyib - dedyic
                    dedzid = -dedzia - dedzib - dedzic
c
c       increment the total bond angle energy and derivatives
c
#if __tver__ & __use_ene__
                    !istat= atomicAdd( ea, e )
                    enrgy = enrgy + tp2enr(e)
#endif

                    associate(ialoc=>ia,ibloc=>ib,icloc=>ic,idloc=>id)
#if __tver__ & __use_grd__
                    ialoc = loc(ia)
                    ibloc = loc(ib)
                    icloc = loc(ic)
                    idloc = loc(id)
                    call atomic_add_f1( grx(ialoc), dedxia )
                    call atomic_add_f1( gry(ialoc), dedyia )
                    call atomic_add_f1( grz(ialoc), dedzia )
c
                    call atomic_add_f1( grx(ibloc), dedxib )
                    call atomic_add_f1( gry(ibloc), dedyib )
                    call atomic_add_f1( grz(ibloc), dedzib )
c
                    call atomic_add_f1( grx(icloc), dedxic )
                    call atomic_add_f1( gry(icloc), dedyic )
                    call atomic_add_f1( grz(icloc), dedzic )
c
                    call atomic_add_f1( grx(idloc), dedxid )
                    call atomic_add_f1( gry(idloc), dedyid )
                    call atomic_add_f1( grz(idloc), dedzid )
#endif
#if __tfea__ & __use_gamd__
                    !aMD storage if waters are considered
                    if (use_amd_wat1) then  ! Check amd
                    if (typeAtm(ia) == aMDwattype(1)
     &             .or. typeAtm(ib) == aMDwattype(1)
     &             .or. typeAtm(ic) == aMDwattype(1)) then   !Check type
                       rrstat= atomicAdd(eW1aMD, real(e,r_p) )
                       it= (ialoc-1)*3
                       call atomic_add_m( deW1aMD(it+1), dedxia )
                       call atomic_add_m( deW1aMD(it+2), dedyia )
                       call atomic_add_m( deW1aMD(it+3), dedzia )
                       it= (ibloc-1)*3
                       call atomic_add_m( deW1aMD(it+1), dedxib )
                       call atomic_add_m( deW1aMD(it+2), dedyib )
                       call atomic_add_m( deW1aMD(it+3), dedzib )
                       it= (icloc-1)*3
                       call atomic_add_m( deW1aMD(it+1), dedxic )
                       call atomic_add_m( deW1aMD(it+2), dedyic )
                       call atomic_add_m( deW1aMD(it+3), dedzic )
                    end if  !Check type
                    end if  !Check amd
#endif
#if __tver__ & __use_vir__
                    if (use_virial) then
                   !increment the internal virial tensor components
                    g_vxx = g_vxx + xad*dedxia + xbd*dedxib + xcd*dedxic
                    g_vxy = g_vxy + yad*dedxia + ybd*dedxib + ycd*dedxic
                    g_vxz = g_vxz + zad*dedxia + zbd*dedxib + zcd*dedxic
                    g_vyy = g_vyy + yad*dedyia + ybd*dedyib + ycd*dedyic
                    g_vyz = g_vyz + zad*dedyia + zbd*dedyib + zcd*dedyic
                    g_vzz = g_vzz + zad*dedzia + zbd*dedzib + zcd*dedzic
                    end if
#endif
                    end associate
                 end if
              end if   ! In_plane angle
           end if      ! Proceed
        end do
        end block angle; end if


        if (ithread.gt.off_tors) then; torsional: block
        integer i,ia,ib,ic,id
        integer itor
#ifdef   USE_NVSHMEM_CUDA
        integer ipe,ind
#endif  
        real(t_p) e,dedphi,rcb
        real(t_p) v1,v2,v3,v4,v5,v6
        real(t_p) c1,c2,c3,c4,c5,c6
        real(t_p) s1,s2,s3,s4,s5,s6
        real(t_p) cosine,cosine2
        real(t_p) cosine3,cosine4
        real(t_p) cosine5,cosine6
        real(t_p) sine,sine2,sine3
        real(t_p) sine4,sine5,sine6
        real(t_p) phi1,phi2,phi3
        real(t_p) phi4,phi5,phi6
        real(t_p) xia,yia,zia
        real(t_p) xib,yib,zib
        real(t_p) xic,yic,zic
        real(t_p) xid,yid,zid
        real(t_p) xba,yba,zba
        real(t_p) xdc,ydc,zdc
        real(t_p) xcb,ycb,zcb
        real(t_p) xca,yca,zca
        real(t_p) xdb,ydb,zdb
        real(t_p) xt,yt,zt,rt2
        real(t_p) xu,yu,zu,ru2
        real(t_p) xtu,ytu,ztu,rtru
        real(t_p) dphi1,dphi2,dphi3
        real(t_p) dphi4,dphi5,dphi6
        real(t_p) dedxt,dedyt,dedzt
        real(t_p) dedxu,dedyu,dedzu
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        real(t_p) dedxid,dedyid,dedzid
        integer proceed
c
c       calculate the torsional angle energy and first derivatives
c
        do itor = ithread, ntorsloc+off_tors, num_tors
           i   = torsglob(itor-off_tors)
#ifdef   USE_NVSHMEM_CUDA
           ipe = (i-1)/ntors_pe
           ind = mod((i-1),ntors_pe) +1
           ia  = d_itors(ipe)%pel(1,ind)
           ib  = d_itors(ipe)%pel(2,ind)
           ic  = d_itors(ipe)%pel(3,ind)
           id  = d_itors(ipe)%pel(4,ind)
#else
           ia = itors(1,i)
           ib = itors(2,i)
           ic = itors(3,i)
           id = itors(4,i)
#endif  
          !decide whether to compute the current interaction
           proceed = merge(1,(iuse(ia).or.iuse(ib).or.iuse(ic)
     &                    .or.iuse(id)), useAll)

          !compute the value of the torsional angle
           if (proceed) then
              xia = x(ia)
              yia = y(ia)
              zia = z(ia)
              xib = x(ib)
              yib = y(ib)
              zib = z(ib)
              xic = x(ic)
              yic = y(ic)
              zic = z(ic)
              xid = x(id)
              yid = y(id)
              zid = z(id)
              xba = xib - xia
              yba = yib - yia
              zba = zib - zia
              xcb = xic - xib
              ycb = yic - yib
              zcb = zic - zib
              xdc = xid - xic
              ydc = yid - yic
              zdc = zid - zic
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xba,yba,zba)
                 call image_inl (xcb,ycb,zcb)
                 call image_inl (xdc,ydc,zdc)
              end if
#endif
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
              rt2 = xt*xt + yt*yt + zt*zt
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
              if (rtru .ne. 0.0_ti_p) then  !Zero
                 rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
                 cosine = (xt*xu + yt*yu + zt*zu) / rtru
                 sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c       set the torsional parameters for this angle
c
#ifdef   USE_NVSHMEM_CUDA
                 v1 = d_tors1(ipe)%pel(1,ind)
                 c1 = d_tors1(ipe)%pel(3,ind)
                 s1 = d_tors1(ipe)%pel(4,ind)
                 v2 = d_tors2(ipe)%pel(1,ind)
                 c2 = d_tors2(ipe)%pel(3,ind)
                 s2 = d_tors2(ipe)%pel(4,ind)
                 v3 = d_tors3(ipe)%pel(1,ind)
                 c3 = d_tors3(ipe)%pel(3,ind)
                 s3 = d_tors3(ipe)%pel(4,ind)
                 v4 = d_tors4(ipe)%pel(1,ind)
                 c4 = d_tors4(ipe)%pel(3,ind)
                 s4 = d_tors4(ipe)%pel(4,ind)
                 v5 = d_tors5(ipe)%pel(1,ind)
                 c5 = d_tors5(ipe)%pel(3,ind)
                 s5 = d_tors5(ipe)%pel(4,ind)
                 v6 = d_tors6(ipe)%pel(1,ind)
                 c6 = d_tors6(ipe)%pel(3,ind)
                 s6 = d_tors6(ipe)%pel(4,ind)
#else
                 v1 = tors1(1,i)
                 c1 = tors1(3,i)
                 s1 = tors1(4,i)
                 v2 = tors2(1,i)
                 c2 = tors2(3,i)
                 s2 = tors2(4,i)
                 v3 = tors3(1,i)
                 c3 = tors3(3,i)
                 s3 = tors3(4,i)
                 v4 = tors4(1,i)
                 c4 = tors4(3,i)
                 s4 = tors4(4,i)
                 v5 = tors5(1,i)
                 c5 = tors5(3,i)
                 s5 = tors5(4,i)
                 v6 = tors6(1,i)
                 c6 = tors6(3,i)
                 s6 = tors6(4,i)
#endif  
c
c       compute the multiple angle trigonometry and the phase terms
c
                 cosine2 = cosine*cosine  - sine*sine
                 sine2   = 2.0_ti_p * cosine * sine
                 cosine3 = cosine*cosine2 - sine*sine2
                 sine3   = cosine*sine2   + sine*cosine2
                 cosine4 = cosine*cosine3 - sine*sine3
                 sine4   = cosine*sine3   + sine*cosine3
                 cosine5 = cosine*cosine4 - sine*sine4
                 sine5   = cosine*sine4   + sine*cosine4
                 cosine6 = cosine*cosine5 - sine*sine5
                 sine6   = cosine*sine5   + sine*cosine5
                 phi1  = 1.0_ti_p + (cosine *c1 + sine *s1)
                 phi2  = 1.0_ti_p + (cosine2*c2 + sine2*s2)
                 phi3  = 1.0_ti_p + (cosine3*c3 + sine3*s3)
                 phi4  = 1.0_ti_p + (cosine4*c4 + sine4*s4)
                 phi5  = 1.0_ti_p + (cosine5*c5 + sine5*s5)
                 phi6  = 1.0_ti_p + (cosine6*c6 + sine6*s6)
                 dphi1 =            (cosine *s1 - sine *c1)
                 dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
                 dphi3 = 3.0_ti_p * (cosine3*s3 - sine3*c3)
                 dphi4 = 4.0_ti_p * (cosine4*s4 - sine4*c4)
                 dphi5 = 5.0_ti_p * (cosine5*s5 - sine5*c5)
                 dphi6 = 6.0_ti_p * (cosine6*s6 - sine6*c6)
c
c       calculate torsional energy and master chain rule term
c
                 e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                         + v4*phi4 + v5*phi5 + v6*phi6)
                 dedphi = torsunit * (v1*dphi1 + v2*dphi2 + v3*dphi3
     &                              + v4*dphi4 + v5*dphi5 + v6*dphi6)
c
c       chain rule terms for first derivative components
c
                 xca = xic - xia
                 yca = yic - yia
                 zca = zic - zia
                 xdb = xid - xib
                 ydb = yid - yib
                 zdb = zid - zib
#if __tfea__ & __use_polymer__
                 if (use_polymer) then
                    call image_inl (xca,yca,zca)
                    call image_inl (xdb,ydb,zdb)
                 end if
#endif
                 dedxt =  dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
                 dedyt =  dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
                 dedzt =  dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
                 dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
                 dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
                 dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c       compute first derivative components for this angle
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
c       increment the total torsional angle energy and gradient
c
#if __tver__ & __use_ene__
                 !call atomic_add_m(et, e)
                 enrgy = enrgy + tp2enr(e)
#endif
#if __tver__ & __use_grd__
                 associate(ialoc=>ia,ibloc=>ib,icloc=>ic,idloc=>id)
                 ialoc = loc(ia)
                 ibloc = loc(ib)
                 icloc = loc(ic)
                 idloc = loc(id)
 
                 call atomic_add_f1( gr_x(ialoc), dedxia )
                 call atomic_add_f1( gr_y(ialoc), dedyia )
                 call atomic_add_f1( gr_z(ialoc), dedzia )
c
                 call atomic_add_f1( gr_x(ibloc), dedxib )
                 call atomic_add_f1( gr_y(ibloc), dedyib )
                 call atomic_add_f1( gr_z(ibloc), dedzib )
c
                 call atomic_add_f1( gr_x(icloc), dedxic )
                 call atomic_add_f1( gr_y(icloc), dedyic )
                 call atomic_add_f1( gr_z(icloc), dedzic )
c
                 call atomic_add_f1( gr_x(idloc), dedxid )
                 call atomic_add_f1( gr_y(idloc), dedyid )
                 call atomic_add_f1( gr_z(idloc), dedzid )
                 end associate
#endif
#if __tver__ & __use_vir__
          !increment the internal virial tensor components
           if (use_virial) then
           g_vxx =g_vxx +xcb*(dedxic+dedxid) -xba*dedxia +xdc*dedxid
           g_vxy =g_vxy +ycb*(dedxic+dedxid) -yba*dedxia +ydc*dedxid
           g_vxz =g_vxz +zcb*(dedxic+dedxid) -zba*dedxia +zdc*dedxid
           g_vyy =g_vyy +ycb*(dedyic+dedyid) -yba*dedyia +ydc*dedyid
           g_vyz =g_vyz +zcb*(dedyic+dedyid) -zba*dedyia +zdc*dedyid
           g_vzz =g_vzz +zcb*(dedzic+dedzid) -zba*dedzia +zdc*dedzid
           end if
#endif
#if __tfea__ & __use_gamd__
                 if (use_amd_dih) call atomic_add_m( eDaMD,e )
#endif
              end if  ! Zero
           end if     ! Proceed
        end do
        end block torsional; end if


        if(ithread.gt.off_pitr) then; pitors: block
        integer i,ipitors
c
c       calculate the pi-orbital torsion angle energy term
c
        do ipitors = ithread, npitorsloc+off_pitr, num_pitr
           i  = pitorsglob(ipitors-off_pitr)
           call ker_pitors
     &         (i,npitorsloc,ipit,kpit,ptorunit,iuse
     &         ,x,y,z,loc,useAll,use_amd_dih,use_virial
     &         ,enrgy,eDaMD,gr_x,gr_y,gr_z
     &         ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
        end do
        end block pitors; end if

        if(ithread.gt.off_totr) then; tortor:block
        integer itortor,iitortor
        !calculate the torsion-torsion interaction energy term
        do itortor = ithread, ntortorloc+off_totr, num_totr
           iitortor = tortorglob(itortor-off_totr)
           call ker_tortor
     &     (iitortor,ntortorloc,n,itt,ibitor,iuse,loc     ! Input
     &     ,i12,n12,atomic,typeAtm
     &     ,x,y,z,tnx,tny,ttx,tty,tbf,tbx,tby,tbxy,ttorunit
     &     ,use_polymer,use_amd_dih,useAll
     &     ,enrgy,eDaMD,gr_x,gr_y,gr_z            ! Output
     &     ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &     )
        end do
        end block tortor; end if


        if(ithread.gt.off_impr) then; improp:block
        integer i,ia,ib,ic,id
        integer ialoc,ibloc,icloc,idloc
        integer iimprop
        real(t_p) e,dedphi
        real(t_p) dt
        real(t_p) ideal,force
        real(t_p) cosine,sine
        real(t_p) rcb,angle
        real(t_p) xt,yt,zt
        real(t_p) xu,yu,zu
        real(t_p) xtu,ytu,ztu
        real(t_p) rt2,ru2,rtru
        real(t_p) xia,yia,zia,xib,yib,zib
        real(t_p) xic,yic,zic,xid,yid,zid
        real(t_p) xba,yba,zba,xcb,ycb,zcb
        real(t_p) xdc,ydc,zdc,xca,yca,zca
        real(t_p) xdb,ydb,zdb
        real(t_p) dedxt,dedyt,dedzt
        real(t_p) dedxu,dedyu,dedzu
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        real(t_p) dedxid,dedyid,dedzid
        integer proceed
c
c       calculate the improper dihedral angle energy term
c
        do iimprop = ithread, niproploc+off_impr, num_impr
           i     = impropglob(iimprop-off_impr)
           ia    = iiprop(1,i)
           ib    = iiprop(2,i)
           ic    = iiprop(3,i)
           id    = iiprop(4,i)
           ialoc = loc(ia)
           ibloc = loc(ib)
           icloc = loc(ic)
           idloc = loc(id)
 
          !decide whether to compute the current interaction
           proceed = merge(1,iuse(ia).or.iuse(ib).or.iuse(ic)
     &                   .or.iuse(id), useAll)
c
c       compute the value of the improper dihedral angle
c
           if (proceed) then
              xia = x(ia)
              yia = y(ia)
              zia = z(ia)
              xib = x(ib)
              yib = y(ib)
              zib = z(ib)
              xic = x(ic)
              yic = y(ic)
              zic = z(ic)
              xid = x(id)
              yid = y(id)
              zid = z(id)
              xba = xib - xia
              yba = yib - yia
              zba = zib - zia
              xcb = xic - xib
              ycb = yic - yib
              zcb = zic - zib
              xdc = xid - xic
              ydc = yid - yic
              zdc = zid - zic
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xba,yba,zba)
                 call image_inl (xcb,ycb,zcb)
                 call image_inl (xdc,ydc,zdc)
              end if
#endif
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
              rt2 = xt*xt + yt*yt + zt*zt
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
              if (rtru .ne. 0.0_ti_p) then
                 rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
                 cosine = (xt*xu + yt*yu + zt*zu) / rtru
                 sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
                 cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
                 angle = radian * acos(cosine)
                 if (sine .lt. 0.0_ti_p)  angle = -angle
c
c       set the improper dihedral parameters for this angle
c
                 ideal = vprop(i)
                 force = kprop(i)
                 if (abs(angle+ideal) .lt. abs(angle-ideal))
     &              ideal = -ideal
                 dt = angle - ideal
                 do while (dt .gt. 180.0_ti_p)
                    dt = dt - 360.0_ti_p
                 end do
                 do while (dt .lt. -180.0_ti_p)
                    dt = dt + 360.0_ti_p
                 end do
c
c       calculate improper energy and master chain rule term
c
                 e      = idihunit * force * dt**2
                 dedphi = 2.0_ti_p * radian * idihunit * force * dt
c
c       chain rule terms for first derivative components
c
                 xca = xic - xia
                 yca = yic - yia
                 zca = zic - zia
                 xdb = xid - xib
                 ydb = yid - yib
                 zdb = zid - zib
#if __tfea__ & __use_polymer__
                 if (use_polymer) then
                    call image_inl (xca,yca,zca)
                    call image_inl (xdb,ydb,zdb)
                 end if
#endif
                 dedxt =  dedphi * (yt*zcb - ycb*zt) /(rt2*rcb)
                 dedyt =  dedphi * (zt*xcb - zcb*xt) /(rt2*rcb)
                 dedzt =  dedphi * (xt*ycb - xcb*yt) /(rt2*rcb)
                 dedxu = -dedphi * (yu*zcb - ycb*zu) /(ru2*rcb)
                 dedyu = -dedphi * (zu*xcb - zcb*xu) /(ru2*rcb)
                 dedzu = -dedphi * (xu*ycb - xcb*yu) /(ru2*rcb)
c
c       compute first derivative components for this angle
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
c       calculate improper dihedral energy and derivatives
c
#if __tver__ & __use_ene__
                 !call atomicÃaddÃ§m( eid,e )
                 enrgy = enrgy + tp2enr(e)
#endif

                 call atomic_add_f1( grx(ialoc), dedxia )
                 call atomic_add_f1( gry(ialoc), dedyia )
                 call atomic_add_f1( grz(ialoc), dedzia )
c
                 call atomic_add_f1( grx(ibloc), dedxib )
                 call atomic_add_f1( gry(ibloc), dedyib )
                 call atomic_add_f1( grz(ibloc), dedzib )
c
                 call atomic_add_f1( grx(icloc), dedxic )
                 call atomic_add_f1( gry(icloc), dedyic )
                 call atomic_add_f1( grz(icloc), dedzic )
c
                 call atomic_add_f1( grx(idloc), dedxid )
                 call atomic_add_f1( gry(idloc), dedyid )
                 call atomic_add_f1( grz(idloc), dedzid )
#if __tver__ & __use_vir__
           !increment the internal virial tensor components
           if (use_virial) then
           g_vxx =g_vxx +xcb*(dedxic+dedxid) -xba*dedxia +xdc*dedxid
           g_vxy =g_vxy +ycb*(dedxic+dedxid) -yba*dedxia +ydc*dedxid
           g_vxz =g_vxz +zcb*(dedxic+dedxid) -zba*dedxia +zdc*dedxid
           g_vyy =g_vyy +ycb*(dedyic+dedyid) -yba*dedyia +ydc*dedyid
           g_vyz =g_vyz +zcb*(dedyic+dedyid) -zba*dedyia +zdc*dedyid
           g_vzz =g_vzz +zcb*(dedzic+dedzid) -zba*dedzia +zdc*dedzid
           end if
#endif
              end if
           end if
        end do
        end block improp; end if


        if(ithread.gt.off_impt) then; imptor:block
        integer i,ia,ib,ic,id
        integer ialoc,ibloc,icloc,idloc
        integer iimptor
        real(t_p) e,dedphi
        real(t_p) rcb
        real(t_p) xt,yt,zt
        real(t_p) xu,yu,zu
        real(t_p) xtu,ytu,ztu
        real(t_p) rt2,ru2,rtru
        real(t_p) v1,v2,v3
        real(t_p) c1,c2,c3
        real(t_p) s1,s2,s3
        real(t_p) sine,cosine
        real(t_p) sine2,cosine2
        real(t_p) sine3,cosine3
        real(t_p) phi1,phi2,phi3
        real(t_p) xia,yia,zia
        real(t_p) xib,yib,zib
        real(t_p) xic,yic,zic
        real(t_p) xid,yid,zid
        real(t_p) xba,yba,zba
        real(t_p) xcb,ycb,zcb
        real(t_p) xdc,ydc,zdc
        real(t_p) xca,yca,zca
        real(t_p) xdb,ydb,zdb
        real(t_p) dphi1,dphi2,dphi3
        real(t_p) dedxt,dedyt,dedzt
        real(t_p) dedxu,dedyu,dedzu
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        real(t_p) dedxid,dedyid,dedzid
        real(r_p) dedx,dedy,dedz
        real(t_p) vxx,vyy,vzz
        real(t_p) vyx,vzx,vzy
        integer proceed
        proceed = 1
c
c       calculate the improper torsional angle energy term
c
        do iimptor = ithread, nitorsloc+off_impt, num_impt
           i  = imptorglob(iimptor-off_impt)
           ia = iitors(1,i)
           ib = iitors(2,i)
           ic = iitors(3,i)
           id = iitors(4,i)
           ialoc = loc(ia)
           ibloc = loc(ib)
           icloc = loc(ic)
           idloc = loc(id)
 
          !decide whether to compute the current interaction
           proceed = merge(1,iuse(ia).or.iuse(ib).or.iuse(ic)
     &                   .or.iuse(id), useAll)
c
c       compute the value of the torsional angle
c
           if (proceed) then
              xia = x(ia)
              yia = y(ia)
              zia = z(ia)
              xib = x(ib)
              yib = y(ib)
              zib = z(ib)
              xic = x(ic)
              yic = y(ic)
              zic = z(ic)
              xid = x(id)
              yid = y(id)
              zid = z(id)
              xba = xib - xia
              yba = yib - yia
              zba = zib - zia
              xcb = xic - xib
              ycb = yic - yib
              zcb = zic - zib
              xdc = xid - xic
              ydc = yid - yic
              zdc = zid - zic
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xba,yba,zba)
                 call image_inl (xcb,ycb,zcb)
                 call image_inl (xdc,ydc,zdc)
              end if
#endif
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
              rt2 = xt*xt + yt*yt + zt*zt
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
              if (rtru .ne. 0.0_ti_p) then
                 rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
                 cosine = (xt*xu + yt*yu + zt*zu) / rtru
                 sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c       set the improper torsional parameters for this angle
c
                 v1 = itors1(1,i)
                 c1 = itors1(3,i)
                 s1 = itors1(4,i)
                 v2 = itors2(1,i)
                 c2 = itors2(3,i)
                 s2 = itors2(4,i)
                 v3 = itors3(1,i)
                 c3 = itors3(3,i)
                 s3 = itors3(4,i)
c
c       compute the multiple angle trigonometry and the phase terms
c
                 cosine2 = cosine*cosine - sine*sine
                 sine2 = 2.0_ti_p * cosine * sine
                 cosine3 = cosine*cosine2 - sine*sine2
                 sine3 = cosine*sine2 + sine*cosine2
                 phi1 = 1.0_ti_p + (cosine*c1 + sine*s1)
                 phi2 = 1.0_ti_p + (cosine2*c2 + sine2*s2)
                 phi3 = 1.0_ti_p + (cosine3*c3 + sine3*s3)
                 dphi1 = (cosine*s1 - sine*c1)
                 dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
                 dphi3 = 3.0_ti_p * (cosine3*s3 - sine3*c3)
c
c       calculate improper torsion energy and master chain rule term
c
                 e = itorunit * (v1*phi1+v2*phi2+v3*phi3)
                 dedphi = itorunit * (v1*dphi1+v2*dphi2+v3*dphi3)
c
c       chain rule terms for first derivative components
c
                 xca = xic - xia
                 yca = yic - yia
                 zca = zic - zia
                 xdb = xid - xib
                 ydb = yid - yib
                 zdb = zid - zib
#if __tfea__ & __use_polymer__
                 if (use_polymer) then
                    call image_inl (xca,yca,zca)
                    call image_inl (xdb,ydb,zdb)
                 end if
#endif
                 dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
                 dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
                 dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
                 dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
                 dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
                 dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c       compute first derivative components for this angle
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
c       increment the improper torsion energy and gradient
c
#if __tver__ & __use_ene__
                 !eit = eit + e
                 enrgy = enrgy + tp2enr(e)
#endif

                 call atomic_add_f1( grx(icloc), dedxic )
                 call atomic_add_f1( gry(icloc), dedyic )
                 call atomic_add_f1( grz(icloc), dedzic )
c
                 call atomic_add_f1( grx(ialoc), dedxia )
                 call atomic_add_f1( gry(ialoc), dedyia )
                 call atomic_add_f1( grz(ialoc), dedzia )
c
                 call atomic_add_f1( grx(ibloc), dedxib )
                 call atomic_add_f1( gry(ibloc), dedyib )
                 call atomic_add_f1( grz(ibloc), dedzib )
c
                 call atomic_add_f1( grx(idloc), dedxid )
                 call atomic_add_f1( gry(idloc), dedyid )
                 call atomic_add_f1( grz(idloc), dedzid )

#if __tver__ & __use_vir__
          !increment the internal virial tensor components
           if (use_virial) then
           g_vxx =g_vxx +xcb*(dedxic+dedxid) -xba*dedxia +xdc*dedxid
           g_vxy =g_vxy +ycb*(dedxic+dedxid) -yba*dedxia +ydc*dedxid
           g_vxz =g_vxz +zcb*(dedxic+dedxid) -zba*dedxia +zdc*dedxid
           g_vyy =g_vyy +ycb*(dedyic+dedyid) -yba*dedyia +ydc*dedyid
           g_vyz =g_vyz +zcb*(dedyic+dedyid) -zba*dedyia +zdc*dedyid
           g_vzz =g_vzz +zcb*(dedzic+dedzid) -zba*dedzia +zdc*dedzid
           end if
#endif
              end if
           end if
        end do
        end block imptor; end if

        if(ithread.gt.off_antr) then; angtor:block
        integer i,k,iangtor
        integer ia,ib,ic,id
        integer iiangtor
        real(t_p) e,e1,e2
        real(t_p) rcb,fgrp
        real(t_p) ddt,dedphi
        real(t_p) rt2,ru2,rtru
        real(t_p) rba2,rcb2,rdc2
        real(t_p) dot,dt
        real(t_p) xt,yt,zt
        real(t_p) xu,yu,zu
        real(t_p) xtu,ytu,ztu
        real(t_p) v1,v2,v3
        real(t_p) c1,c2,c3
        real(t_p) s1,s2,s3
        real(t_p) sine,cosine
        real(t_p) sine2,cosine2
        real(t_p) sine3,cosine3
        real(t_p) phi1,phi2,phi3
        real(t_p) dphi1,dphi2,dphi3
        real(t_p) angle1,cosang
        real(t_p) xba,yba,zba
        real(t_p) xcb,ycb,zcb
        real(t_p) xdc,ydc,zdc
        real(t_p) xca,yca,zca
        real(t_p) xdb,ydb,zdb
        real(t_p) terma,termb
        real(t_p) termc,termd
        real(t_p) dedxt,dedyt,dedzt
        real(t_p) dedxu,dedyu,dedzu
        real(t_p) dedxia,dedyia,dedzia
        real(t_p) dedxib,dedyib,dedzib
        real(t_p) dedxic,dedyic,dedzic
        real(t_p) dedxid,dedyid,dedzid
        integer proceed
c
c       calculate the angle-torsion energy and first derviatives
c
        do iangtor = ithread, nangtorloc+off_antr, num_antr
           iiangtor = angtorglob(iangtor-off_antr)
           i  = iat(1,iiangtor)
           ia = itors(1,i)
           ib = itors(2,i)
           ic = itors(3,i)
           id = itors(4,i)
c
c       decide whether to compute the current interaction
c
c          if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
           proceed = merge(1,(iuse(ia).or.iuse(ib).or.iuse(ic)
     &                    .or.iuse(id)) , useAll)
c
c       compute the value of the torsional angle
c
           if (proceed) then
              block
              real(t_p) xia,xib,xic,xid,yia,yib,yic,yid,zia,zib,zic,zid
              xia = x(ia)
              yia = y(ia)
              zia = z(ia)
              xib = x(ib)
              yib = y(ib)
              zib = z(ib)
              xic = x(ic)
              yic = y(ic)
              zic = z(ic)
              xid = x(id)
              yid = y(id)
              zid = z(id)
              xba = xib - xia
              yba = yib - yia
              zba = zib - zia
              xcb = xic - xib
              ycb = yic - yib
              zcb = zic - zib
              xdc = xid - xic
              ydc = yid - yic
              zdc = zid - zic
              xca = xic - xia
              yca = yic - yia
              zca = zic - zia
              xdb = xid - xib
              ydb = yid - yib
              zdb = zid - zib
              end block
#if __tfea__ & __use_polymer__
              if (use_polymer) then
                 call image_inl (xba,yba,zba)
                 call image_inl (xcb,ycb,zcb)
                 call image_inl (xdc,ydc,zdc)
              end if
#endif
              rba2 = xba*xba + yba*yba + zba*zba
              rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
              rdc2 = xdc*xdc + ydc*ydc + zdc*zdc
              if (min(rba2,rcb2,rdc2) .ne. 0.0_ti_p) then
                 xt = yba*zcb - ycb*zba
                 yt = zba*xcb - zcb*xba
                 zt = xba*ycb - xcb*yba
                 xu = ycb*zdc - ydc*zcb
                 yu = zcb*xdc - zdc*xcb
                 zu = xcb*ydc - xdc*ycb
                 xtu = yt*zu - yu*zt
                 ytu = zt*xu - zu*xt
                 ztu = xt*yu - xu*yt
                 rt2 = xt*xt + yt*yt + zt*zt
                 rt2 = max(rt2,0.000001_ti_p)
                 ru2 = xu*xu + yu*yu + zu*zu
                 ru2 = max(ru2,0.000001_ti_p)
                 rtru = sqrt(rt2*ru2)
#if __tfea__ & __use_polymer__
                 if (use_polymer) then
                    call image_inl (xca,yca,zca)
                    call image_inl (xdb,ydb,zdb)
                 end if
#endif
                 rcb = sqrt(rcb2)
                 cosine = (xt*xu + yt*yu + zt*zu) / rtru
                 sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c       compute multiple angle trigonometry and phase terms
c
                 c1 = tors1(3,i)
                 s1 = tors1(4,i)
                 c2 = tors2(3,i)
                 s2 = tors2(4,i)
                 c3 = tors3(3,i)
                 s3 = tors3(4,i)
                 cosine2 = cosine*cosine - sine*sine
                 sine2   = 2.0_ti_p * cosine * sine
                 cosine3 = cosine*cosine2 - sine*sine2
                 sine3   = cosine*sine2 + sine*cosine2
                 phi1  = 1.0_ti_p + (cosine*c1 + sine*s1)
                 phi2  = 1.0_ti_p + (cosine2*c2 + sine2*s2)
                 phi3  = 1.0_ti_p + (cosine3*c3 + sine3*s3)
                 dphi1 = (cosine*s1 - sine*c1)
                 dphi2 = 2.0_ti_p * (cosine2*s2 - sine2*c2)
                 dphi3 = 3.0_ti_p * (cosine3*s3 - sine3*c3)
c
c       set the angle-torsion parameters for the first angle
c
                 v1 = kant(1,iiangtor)
                 v2 = kant(2,iiangtor)
                 v3 = kant(3,iiangtor)
                 k  = iat(2,iiangtor)
                 dot = xba*xcb + yba*ycb + zba*zcb
                 cosang = -dot / sqrt(rba2*rcb2)
                 angle1 = radian * acos(cosang)
                 dt = angle1 - anat(k)
                 e1 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
                 dedphi = atorunit*dt*(v1*dphi1 + v2*dphi2 + v3*dphi3)
                 ddt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
cc
cc       scale the interaction based on its group membership
cc
c                 if (use_group) then
c                    e1 = e1 * fgrp
c                    dedphi = dedphi * fgrp
c                    ddt = ddt * fgrp
c                 end if
c
c       compute derivative components for this interaction
c
                 dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb)
                 dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb)
                 dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb)
                 dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb)
                 dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb)
                 dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb)
c
c       increment chain rule components for the first angle
c
                 terma = -ddt / (rba2*sqrt(rt2))
                 termc = ddt / (rcb2*sqrt(rt2))
                 dedxia = terma*(zba*yt-yba*zt) + zcb*dedyt - ycb*dedzt
                 dedyia = terma*(xba*zt-zba*xt) + xcb*dedzt - zcb*dedxt
                 dedzia = terma*(yba*xt-xba*yt) + ycb*dedxt - xcb*dedyt
                 dedxib = terma*(yba*zt-zba*yt) + termc*(zcb*yt-ycb*zt)
     &                       + yca*dedzt - zca*dedyt
     &                       + zdc*dedyu - ydc*dedzu
                 dedyib = terma*(zba*xt-xba*zt) + termc*(xcb*zt-zcb*xt)
     &                       + zca*dedxt - xca*dedzt
     &                       + xdc*dedzu - zdc*dedxu
                 dedzib = terma*(xba*yt-yba*xt) + termc*(ycb*xt-xcb*yt)
     &                       + xca*dedyt - yca*dedxt
     &                       + ydc*dedxu - xdc*dedyu
                 dedxic = termc*(ycb*zt-zcb*yt) + zba*dedyt
     &                       - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyic = termc*(zcb*xt-xcb*zt) + xba*dedzt
     &                       - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzic = termc*(xcb*yt-ycb*xt) + yba*dedxt
     &                       - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxid = zcb*dedyu - ycb*dedzu
                 dedyid = xcb*dedzu - zcb*dedxu
                 dedzid = ycb*dedxu - xcb*dedyu
c
c       get the angle-torsion values for the second angle
c
                 v1 = kant(4,iiangtor)
                 v2 = kant(5,iiangtor)
                 v3 = kant(6,iiangtor)
                 k = iat(3,iiangtor)
                 dot = xcb*xdc + ycb*ydc + zcb*zdc
                 cosang = acos(-dot / sqrt(rcb2*rdc2))
                 angle1 = radian * cosang
                 dt = angle1 - anat(k)
                 e2 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
                 dedphi = atorunit*dt*(v1*dphi1 + v2*dphi2 + v3*dphi3)
                 ddt = atorunit* radian* (v1*phi1 + v2*phi2 + v3*phi3)
cc
cc       scale the interaction based on its group membership
cc
c                 if (use_group) then
c                    e2 = e2 * fgrp
c                    dedphi = dedphi * fgrp
c                    ddt = ddt * fgrp
c                 end if
c
c       compute derivative components for this interaction
c
                 dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb)
                 dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb)
                 dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb)
                 dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb)
                 dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb)
                 dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb)
c
c       increment chain rule components for the second angle
c
                 termb  = -ddt / (rcb2*sqrt(ru2))
                 termd  = ddt / (rdc2*sqrt(ru2))
                 dedxia = dedxia + zcb*dedyt - ycb*dedzt
                 dedyia = dedyia + xcb*dedzt - zcb*dedxt
                 dedzia = dedzia + ycb*dedxt - xcb*dedyt
                 dedxib = dedxib + termb*(zcb*yu-ycb*zu) + yca*dedzt
     &                       - zca*dedyt + zdc*dedyu - ydc*dedzu
                 dedyib = dedyib + termb*(xcb*zu-zcb*xu) + zca*dedxt
     &                       - xca*dedzt + xdc*dedzu - zdc*dedxu
                 dedzib = dedzib + termb*(ycb*xu-xcb*yu) + xca*dedyt
     &                       - yca*dedxt + ydc*dedxu - xdc*dedyu
                 dedxic = dedxic + termb*(ycb*zu-zcb*yu)
     &                       + termd*(zdc*yu-ydc*zu) + zba*dedyt
     &                       - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyic = dedyic + termb*(zcb*xu-xcb*zu)
     &                       + termd*(xdc*zu-zdc*xu) + xba*dedzt
     &                       - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzic = dedzic + termb*(xcb*yu-ycb*xu)
     &                       + termd*(ydc*xu-xdc*yu) + yba*dedxt
     &                       - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxid = dedxid + termd*(ydc*zu-zdc*yu)
     &                       + zcb*dedyu - ycb*dedzu
                 dedyid = dedyid + termd*(zdc*xu-xdc*zu)
     &                       + xcb*dedzu - zcb*dedxu
                 dedzid = dedzid + termd*(xdc*yu-ydc*xu)
     &                       + ycb*dedxu - xcb*dedyu

#if __tver__ & __use_ene__
                !increment the angle-torsion energy and gradient
                 enrgy   = enrgy + tp2enr( e1+e2 )
                 !call atomic_add_m( eat,e )
#endif

#if __tver__ & __use_grd__
                 associate(ialoc=>ia,ibloc=>ib,icloc=>ic,idloc=>id)
                 ialoc = loc(ia)
                 ibloc = loc(ib)
                 icloc = loc(ic)
                 idloc = loc(id)
                 call atomic_add_f1( grx(ialoc), dedxia )
                 call atomic_add_f1( gry(ialoc), dedyia )
                 call atomic_add_f1( grz(ialoc), dedzia )
                 call atomic_add_f1( grx(ibloc), dedxib )
                 call atomic_add_f1( gry(ibloc), dedyib )
                 call atomic_add_f1( grz(ibloc), dedzib )
                 call atomic_add_f1( grx(icloc), dedxic )
                 call atomic_add_f1( gry(icloc), dedyic )
                 call atomic_add_f1( grz(icloc), dedzic )
                 call atomic_add_f1( grx(idloc), dedxid )
                 call atomic_add_f1( gry(idloc), dedyid )
                 call atomic_add_f1( grz(idloc), dedzid )
                 end associate
#endif
#if __tver__ & __use_vir__
           if (use_virial) then
          !increment the internal viriaens components
           g_vxx =g_vxx +xcb*(dedxic+dedxid) -xba*dedxia +xdc*dedxid
           g_vxy =g_vxy +ycb*(dedxic+dedxid) -yba*dedxia +ydc*dedxid
           g_vxz =g_vxz +zcb*(dedxic+dedxid) -zba*dedxia +zdc*dedxid
           g_vyy =g_vyy +ycb*(dedyic+dedyid) -yba*dedyia +ydc*dedyid
           g_vyz =g_vyz +zcb*(dedyic+dedyid) -zba*dedyia +zdc*dedyid
           g_vzz =g_vzz +zcb*(dedzic+dedzid) -zba*dedzia +zdc*dedzid
           end if
#endif
              end if
           end if
        end do
        end block angtor; end if

        it    = iand(ithread-1,RED_BUFF_SIZE-1) + 1
#if __tver__ & __use_ene__
        ! Increment energy term of bonded terms
        rrstat = atomicAdd( ev_buff(it) ,enrgy )
#endif
#if __tver__ & __use_vir__
        if (use_virial) then
        ! Increment virial term of bonded terms
        rstat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),g_vxx )
        rstat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),g_vxy )
        rstat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),g_vxz )
        rstat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),g_vyy )
        rstat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),g_vyz )
        rstat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),g_vzz )
        end if
#endif
        end subroutine eliaison1_kcu
      end module

      module eliaisoncu_dat
      use deriv, only: tdes_l
      use tinMemory,only: mipk
      implicit none
      mdyn_rtyp,device,pointer,dimension(:)::grx,gry,grz
     &         ,grtx,grty,grtz

      contains

      subroutine attach_p(buff1,stride,off_ws,opt)
      use atoms ,only:n
      use angle ,only:nangle,nangleloc
      use angang,only:nangang,nangangloc
      use angtor,only:nangtor,nangtorloc
      use bond  ,only:nbond,nbondloc
      use bitor ,only:nbitor,nbitorloc
      use charge,only:nion,nionloc
      !use domdec
      use improp,only:niprop,niproploc
      use imptor,only:nitors,nitorsloc
      use kgeoms
      use mpole ,only:npole,npoleloc
      use mutant,only:nmut
      use strbnd,only:nstrbnd,nstrbndloc
      use opbend,only:nopbend,nopbendloc
      use opdist,only:nopdist,nopdistloc
      !use potent
      use polar ,only:npolar
      !use mpi
      use pitors,only:npitors,npitorsloc
      use strtor,only:nstrtor,nstrtorloc
      use tors  ,only:ntors,ntorsloc
      use tortor,only:ntortor,ntortorloc
      use urey  ,only:nurey,nureyloc
      use vdw   ,only:nvdw,nvdwloc
      implicit none
      mdyn_rtyp,device,intent(in),target:: buff1(*)
      integer      ,intent(in):: stride,opt
      integer(mipk),intent(in):: off_ws
      integer(mipk) o1

      if (opt.eq.2) then
         if (tdes_l) then
            grx (1:stride) => buff1(0*stride+1:1*stride)
            gry (1:stride) => buff1(1*stride+1:2*stride)
            grz (1:stride) => buff1(2*stride+1:3*stride)
            grtx(1:stride) => buff1(0*stride+1:1*stride)
            grty(1:stride) => buff1(1*stride+1:2*stride)
            grtz(1:stride) => buff1(2*stride+1:3*stride)
         else
            o1  = off_ws
            grx (1:stride) => buff1(o1+0*stride+1:o1+1*stride)
            gry (1:stride) => buff1(o1+1*stride+1:o1+2*stride)
            grz (1:stride) => buff1(o1+2*stride+1:o1+3*stride)
            grtx(1:stride) => buff1(o1+3*stride+1:o1+4*stride)
            grty(1:stride) => buff1(o1+4*stride+1:o1+5*stride)
            grtz(1:stride) => buff1(o1+5*stride+1:o1+6*stride)
         end if
      end if

      end subroutine

      end module
