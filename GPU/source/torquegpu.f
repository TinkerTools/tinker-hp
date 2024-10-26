c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine torquegpu  --  convert single site torque to force  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "torque2gpu" takes the torque values on multiple sites defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame on the device
c
c     npolelocnl and npolelocnlloop came from mpole module
c
c     Matching between polaxe and ipolaxe
c     O  ->  'None'
c     1  ->  '3-Fold'
c     2  ->  'Bisector'
c     4  ->  'Z-bisect'
c     8  ->  'Z-Only'
c     16 ->  'Z-Then-X'
c
#include "tinker_macro.h"
      module torquegpu_inl
      use atmlst
      use atoms
      use deriv
      use domdec
      use inform ,only: deb_Path
      use mpole  ,only: npolelocnl,xaxis,yaxis,zaxis,ipole,ipolaxe
     &           ,Ax_3_Fold,Ax_Bisector,Ax_Z_Bisect,Ax_Z_Only
     &           ,Ax_Z_Then_X,Ax_None
      use sizes
      use mpi
      use utilgpu
      use timestat

      implicit none
#include "atomicOp.h.f"

      contains
#include "convert.f.inc"
#include "atomicOp.inc.f"

      end module

      subroutine torquedgpu (trqvec,frcx,frcy,frcz,de,extract)
      use torquegpu_inl
      implicit none
      integer i,iipole,iglob,iloc,j,ca
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer axetyp
      real(t_p) du,dv,dw,dt
      real(t_p) dot
      real(t_p) usiz,vsiz,wsiz
      real(t_p) psiz,rsiz,ssiz
      real(t_p) t1siz,t2siz
      real(t_p) uvsiz,uwsiz,vwsiz
      real(t_p) ursiz,ussiz
      real(t_p) vssiz,wssiz
      real(t_p) delsiz,dphiddel
      real(t_p) uvcos,uwcos,vwcos
      real(t_p) urcos,uscos
      real(t_p) vscos,wscos
      real(t_p) rwcos,wpcos
      real(t_p) upcos,vpcos
      real(t_p) rucos,rvcos
      real(t_p) ut1cos,ut2cos
      real(t_p) uvsin,uwsin,vwsin
      real(t_p) ursin,ussin
      real(t_p) vssin,wssin
      real(t_p) rwsin,wpsin
      real(t_p) upsin,vpsin
      real(t_p) rusin,rvsin
      real(t_p) ut1sin,ut2sin
      real(t_p) dphidu,dphidv,dphidw
      real(t_p) dphidr,dphids
      real(t_p) ux,uy,uz,vx,vy,vz,wx,wy,wz
      real(t_p) px,py,pz,rx,ry,rz,sx,sy,sz
      real(t_p) t1x,t1y,t1z,t2x,t2y,t2z
      real(t_p) uvx,uvy,uvz,uwx,uwy,uwz,vwx,vwy,vwz
      real(t_p) urx,ury,urz,usx,usy,usz
      real(t_p) vsx,vsy,vsz,wsx,wsy,wsz
      real(t_p) delx,dely,delz,epsx,epsy,epsz
      real(t_p) trq(3)
      real(t_p) del(3),eps(3)
      real(t_p),intent(in):: trqvec(:,:)
      real(t_p),intent(out):: frcx(:,:),frcy(:,:),frcz(:,:)
      real(md_p),intent(inout):: de(:,:)
      logical*1,intent(in):: extract !for trqvec depending on his size
!$acc routine(distprocpart1)
c
      if (deb_Path) write(*,'(3x,a)') 'torquedgpu'
      call timer_enter( timer_torque )

      if (tinkerdebug.gt.0) then
         j=0

!$acc parallel loop gang vector_length(128) async(dir_queue)
!$acc&         present(ipolaxe,loc,ipole,
!$acc&  xaxis,yaxis,zaxis,poleglobnl) copyin(j)
      do i =1, npolelocnl
        iipole = poleglobnl(i)
        axetyp = ipolaxe(iipole)

        if (axetyp .eq. Ax_None) cycle

        ia    = zaxis(iipole)
        ib    = ipole(iipole)  !iglob
        ic    = xaxis(iipole)
        id    = yaxis(iipole)
        ibloc = loc(ib) !iloc
        ialoc = -1; idloc=-1; icloc=-1;
        if (ia.gt.0) then; ialoc = loc(ia); else; ialoc=ibloc; endif
        if (ic.gt.0) then; icloc = loc(ic); else; icloc=ibloc; endif
        if (id.gt.0) then; idloc = loc(id); else; idloc=ibloc; endif

        if (ibloc.eq.0.or.ibloc.gt.nbloc) then
!$acc atomic capture
           j = j+1
           ca = j
!$acc end atomic
           if(ca.le.10) print*,'WARNING!!! torquedgpu Outside loc'
           print*,'b ',iipole,ib,ibloc,rank
        end if
        if (ialoc.eq.0.or.ialoc.gt.nbloc) then
!$acc atomic capture
           j = j+1
           ca = j
!$acc end atomic
           if(ca.le.10) print*,'WARNING!!! torquedgpu Outside loc'
           call distprocpart1(ib,rank,uz,.true.,x,y,z)
           print*,'a ',iipole,ia,ialoc,ib,ibloc,nbloc,uz,rank
        end if
        if (icloc.eq.0.or.icloc.gt.nbloc) then
           call distprocpart1(ib,rank,ux,.true.,x,y,z)
           print*,'c ',iipole,ic,icloc,ib,ibloc,nbloc,ux,rank
        end if
        if (idloc.eq.0.or.idloc.gt.nbloc) then
           call distprocpart1(id,rank,uy,.true.,x,y,z)
           print*,'d ',iipole,id,idloc,ib,ibloc,nbloc,uy,rank
        end if
      end do

      end if
c
c     zero out force components on local frame-defining atoms
c
!$acc parallel loop collapse(2) async(dir_queue)
!$acc&         present(frcx,frcy,frcz)
      do i = 1, npolelocnl
         do j = 1, 3
            frcz(j,i) = 0.0_ti_p
            frcx(j,i) = 0.0_ti_p
            frcy(j,i) = 0.0_ti_p
         end do
      end do
c
c     get the local frame type and the frame-defining atoms
c
!$acc parallel loop gang vector_length(64) async(dir_queue)
!$acc&         private(trq,del,eps)
!$acc&         present(trqvec,frcx,frcy,frcz,de)
!$acc&         present(x,y,z,ipolaxe,loc,ipole,
!$acc&  xaxis,yaxis,zaxis,poleglobnl)
      do i =1, npolelocnl
        iipole = poleglobnl(i)
        axetyp = ipolaxe(iipole)

        if (axetyp .eq. Ax_None) cycle

        ia    = zaxis(iipole)
        ib    = ipole(iipole)  !iglob
        ic    = xaxis(iipole)
        ibloc = loc(ib) !iloc
        ialoc = merge(loc(ia),nbloc,ia.gt.0)
        icloc = merge(loc(ic),nbloc,ic.gt.0)
        ! trqvec can be either npolelocnl or nbloc size
        iloc  = merge(ibloc,i,extract)
        ! Check axis local Id
        if (ialoc.eq.0.or.ialoc.gt.nbloc) cycle
        if (icloc.eq.0.or.icloc.gt.nbloc) cycle

        trq(1) = trqvec(1,iloc)
        trq(2) = trqvec(2,iloc)
        trq(3) = trqvec(3,iloc)
c
c     construct the three rotation axes for the local frame
c
        ux   = x(ia) - x(ib)
        uy   = y(ia) - y(ib)
        uz   = z(ia) - z(ib)
        usiz = sqrt(ux*ux + uy*uy + uz*uz)

        if (axetyp .ne. Ax_Z_Only) then
           vx   = x(ic) - x(ib)
           vy   = y(ic) - y(ib)
           vz   = z(ic) - z(ib)
           vsiz = sqrt(vx*vx + vy*vy + vz*vz)
        else
           vx   = 1.0_ti_p
           vy   = 0.0_ti_p
           vz   = 0.0_ti_p
           vsiz = 1.0_ti_p
           dot  = ux / usiz
           if (abs(dot) .gt. 0.866_ti_p) then
              vx = 0.0_ti_p
              vy = 1.0_ti_p
           end if
        end if

        if (axetyp.eq.Ax_Z_Bisect .or. axetyp.eq.Ax_3_Fold) then
           id    = yaxis(iipole)
           idloc = merge(loc(id),nbloc,id.gt.0)
           ! Check axis local Id
           if (idloc.eq.0.or.idloc.gt.nbloc) cycle

           wx = x(id) - x(ib)
           wy = y(id) - y(ib)
           wz = z(id) - z(ib)
        else
           wx = uy*vz - uz*vy
           wy = uz*vx - ux*vz
           wz = ux*vy - uy*vx
        end if

        wsiz  = sqrt(wx*wx + wy*wy + wz*wz)
        ux    = ux / usiz
        uy    = uy / usiz
        uz    = uz / usiz
        vx    = vx / vsiz
        vy    = vy / vsiz
        vz    = vz / vsiz
        wx    = wx / wsiz
        wy    = wy / wsiz
        wz    = wz / wsiz
c
c     find the perpendicular and angle for each pair of axes
c
        uvx   = vy*uz - vz*uy
        uvy   = vz*ux - vx*uz
        uvz   = vx*uy - vy*ux
        uvsiz = sqrt(uvx*uvx + uvy*uvy + uvz*uvz)
        uwx   = wy*uz - wz*uy
        uwy   = wz*ux - wx*uz
        uwz   = wx*uy - wy*ux
        uwsiz = sqrt(uwx*uwx + uwy*uwy + uwz*uwz)
        vwx   = wy*vz - wz*vy
        vwy   = wz*vx - wx*vz
        vwz   = wx*vy - wy*vx
        vwsiz = sqrt(vwx*vwx + vwy*vwy + vwz*vwz)

        uvx   = uvx / uvsiz
        uvy   = uvy / uvsiz
        uvz   = uvz / uvsiz
        uwx   = uwx / uwsiz
        uwy   = uwy / uwsiz
        uwz   = uwz / uwsiz
        vwx   = vwx / vwsiz
        vwy   = vwy / vwsiz
        vwz   = vwz / vwsiz
c
c     get sine and cosine of angles between the rotation axes
c
        uvcos = ux*vx + uy*vy + uz*vz
        uvsin = sqrt(abs(1.0_ti_p - uvcos*uvcos))
        uwcos = ux*wx + uy*wy + uz*wz
        uwsin = sqrt(abs(1.0_ti_p - uwcos*uwcos))
        vwcos = vx*wx + vy*wy + vz*wz
        vwsin = sqrt(abs(1.0_ti_p - vwcos*vwcos))
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
        dphidu = -trq(1)*ux - trq(2)*uy - trq(3)*uz
        dphidv = -trq(1)*vx - trq(2)*vy - trq(3)*vz
        dphidw = -trq(1)*wx - trq(2)*wy - trq(3)*wz

        if (axetyp .eq. Ax_Z_Bisect) then
c
c     build some additional axes needed for the Z-Bisect method
c
           rx    = vx + wx
           ry    = vy + wy
           rz    = vz + wz
           rsiz  = sqrt(rx*rx + ry*ry + rz*rz)
           sx    = uy*rz - uz*ry
           sy    = uz*rx - ux*rz
           sz    = ux*ry - uy*rx
           ssiz  = sqrt(sx*sx + sy*sy + sz*sz)
                 
           rsiz  = 1.0_ti_p / rsiz
           ssiz  = 1.0_ti_p / ssiz
           rx    = rx * rsiz
           ry    = ry * rsiz
           rz    = rz * rsiz
           sx    = sx * ssiz
           sy    = sy * ssiz
           sz    = sz * ssiz

           urx   = ry*uz - rz*uy
           ury   = rz*ux - rx*uz
           urz   = rx*uy - ry*ux
           ursiz = sqrt(urx*urx + ury*ury + urz*urz)
           usx   = sy*uz - sz*uy
           usy   = sz*ux - sx*uz
           usz   = sx*uy - sy*ux
           ussiz = sqrt(usx*usx + usy*usy + usz*usz)
           vsx   = sy*vz - sz*vy
           vsy   = sz*vx - sx*vz
           vsz   = sx*vy - sy*vx
           vssiz = sqrt(vsx*vsx + vsy*vsy + vsz*vsz)
           wsx   = sy*wz - sz*wy
           wsy   = sz*wx - sx*wz
           wsz   = sx*wy - sy*wx
           wssiz = sqrt(wsx*wsx + wsy*wsy + wsz*wsz)

           urx   = urx / ursiz
           ury   = ury / ursiz
           urz   = urz / ursiz
           usx   = usx / ussiz
           usy   = usy / ussiz
           usz   = usz / ussiz
           vsx   = vsx / vssiz
           vsy   = vsy / vssiz
           vsz   = vsz / vssiz
           wsx   = wsx / wssiz
           wsy   = wsy / wssiz
           wsz   = wsz / wssiz
c
c     get sine and cosine of angles between the rotation axes
c
           urcos = ux*rx + uy*ry + uz*rz
           ursin = sqrt(abs(1.0_ti_p - urcos*urcos))
           uscos = ux*sx + uy*sy + uz*sz
           ussin = sqrt(abs(1.0_ti_p - uscos*uscos))
           vscos = vx*sx + vy*sy + vz*sz
           vssin = sqrt(abs(1.0_ti_p - vscos*vscos))
           wscos = wx*sx + wy*sy + wz*sz
           wssin = sqrt(abs(1.0_ti_p - wscos*wscos))
c
c     compute the projection of v and w onto the ru-plane
c
           t1x   = vx - sx*vscos
           t1y   = vy - sy*vscos
           t1z   = vz - sz*vscos
           t2x   = wx - sx*wscos
           t2y   = wy - sy*wscos
           t2z   = wz - sz*wscos

           t1siz = sqrt(t1x*t1x+t1y*t1y+t1z*t1z)
           t2siz = sqrt(t2x*t2x+t2y*t2y+t2z*t2z)

           t1x   = t1x / t1siz
           t1y   = t1y / t1siz
           t1z   = t1z / t1siz
           t2x   = t2x / t2siz
           t2y   = t2y / t2siz
           t2z   = t2z / t2siz

           ut1cos = ux*t1x + uy*t1y + uz*t1z
           ut1sin = sqrt(abs(1.0_ti_p - ut1cos*ut1cos))
           ut2cos = ux*t2x + uy*t2y + uz*t2z
           ut2sin = sqrt(abs(1.0_ti_p - ut2cos*ut2cos))
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           dphids = -trq(1)*sx - trq(2)*sy - trq(3)*sz
        end if
c
c     force distribution for the Z-Only local coordinate method
c
        if (axetyp .eq. Ax_Z_Only) then

           du = uvx*dphidv/(usiz*uvsin)
     &        + uwx*dphidw/usiz
           call atomic_add(de(1,ialoc),real( du,md_p))
           call atomic_add(de(1,ibloc),real(-du,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           du = uvy*dphidv/(usiz*uvsin)
     &        + uwy*dphidw/usiz
           call atomic_add(de(2,ialoc),real( du,md_p))
           call atomic_add(de(2,ibloc),real(-du,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           du = uvz*dphidv/(usiz*uvsin)
     &        + uwz*dphidw/usiz
           call atomic_add(de(3,ialoc),real( du,md_p))
           call atomic_add(de(3,ibloc),real(-du,md_p))
           frcz(3,i)   = frcz(3,i)   + du
c
c     force distribution for the Z-then-X local coordinate method
c
        else if (axetyp .eq. Ax_Z_Then_X) then

           du =  uvx*dphidv/(usiz*uvsin) + uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           du =  uvy*dphidv/(usiz*uvsin) + uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           du =  uvz*dphidv/(usiz*uvsin) + uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
           frcz(3,i)   = frcz(3,i) + du
           frcx(3,i)   = frcx(3,i) + dv

c
c     force distribution for the Bisector local coordinate method
c
        else if (axetyp .eq. Ax_Bisector) then

           du =  uvx*dphidv/(usiz*uvsin) + 0.5_ti_p*uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwx*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           du =  uvy*dphidv/(usiz*uvsin) + 0.5_ti_p*uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwy*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           du =  uvz*dphidv/(usiz*uvsin) + 0.5_ti_p*uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwz*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
           frcz(3,i)   = frcz(3,i)   + du
           frcx(3,i)   = frcx(3,i)   + dv

c
c     force distribution for the Z-Bisect local coordinate method
c
        else if (axetyp .eq. Ax_Z_Bisect) then

           du = urx*dphidr/(usiz*ursin) 
     &        + usx*dphids/usiz
           dv = (vssin*sx-vscos*t1x)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sx-wscos*t2x)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,idloc),real(dw,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           frcy(1,i)   = frcy(1,i)   + dw
           du = ury*dphidr/(usiz*ursin) + usy*dphids/usiz
           dv = (vssin*sy-vscos*t1y)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sy-wscos*t2y)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,idloc),real(dw,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           frcy(2,i)   = frcy(2,i)   + dw
           du = urz*dphidr/(usiz*ursin) + usz*dphids/usiz
           dv = (vssin*sz-vscos*t1z)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sz-wscos*t2z)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,idloc),real(dw,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
           frcz(3,i)   = frcz(3,i)   + du
           frcx(3,i)   = frcx(3,i)   + dv
           frcy(3,i)   = frcy(3,i)   + dw

c
c     force distribution for the 3-Fold local coordinate method
c
        else if (axetyp .eq. Ax_3_Fold) then
           px     = ux + vx + wx
           py     = uy + vy + wy
           pz     = uz + vz + wz
           psiz   = sqrt(px*px + py*py + pz*pz) 
                  
           px     = px / psiz
           py     = py / psiz
           pz     = pz / psiz
                  
           wpcos  = wx*px + wy*py + wz*pz
           upcos  = ux*px + uy*py + uz*pz
           vpcos  = vx*px + vy*py + vz*pz
           wpsin  = sqrt(abs(1.0_ti_p - wpcos*wpcos))
           upsin  = sqrt(abs(1.0_ti_p - upcos*upcos))
           vpsin  = sqrt(abs(1.0_ti_p - vpcos*vpcos))
                  
           rx     = ux + vx
           ry     = uy + vy
           rz     = uz + vz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rwcos  = rx*wx + ry*wy + rz*wz
           rwsin  = sqrt(abs(1.0_ti_p - rwcos*rwcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*wz - rz*wy
           del(2) = rz*wx - rx*wz
           del(3) = rx*wy - ry*wx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*wz - del(3)*wy
           eps(2) = del(3)*wx - del(1)*wz
           eps(3) = del(1)*wy - del(2)*wx

!$acc loop seq
           do j = 1, 3
              dw = del(j)*dphidr/(wsiz*rwsin)
     &           + eps(j)*dphiddel*wpcos/(wsiz*psiz) 
              call atomic_add(de(j,idloc),real( dw,md_p))
              call atomic_add(de(j,ibloc),real(-dw,md_p))
              frcy(j,i)   = frcy(j,i)   + dw
           end do         
           rx     = vx + wx
           ry     = vy + wy
           rz     = vz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)

           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rucos  = rx*ux + ry*uy + rz*uz
           rusin  = sqrt(abs(1.0_ti_p - rucos*rucos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*uz - rz*uy
           del(2) = rz*ux - rx*uz
           del(3) = rx*uy - ry*ux
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*uz - del(3)*uy
           eps(2) = del(3)*ux - del(1)*uz
           eps(3) = del(1)*uy - del(2)*ux
!$acc loop seq
           do j = 1, 3
              du = del(j)*dphidr/(usiz*rusin)
     &           + eps(j)*dphiddel*upcos/(usiz*psiz)
              call atomic_add(de(j,ialoc),real( du,md_p))
              call atomic_add(de(j,ibloc),real(-du,md_p))
              frcz(j,i)   = frcz(j,i)   + du
           end do
           rx     = ux + wx
           ry     = uy + wy
           rz     = uz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rvcos  = rx*vx + ry*vy + rz*vz 
           rvsin  = sqrt(abs(1.0_ti_p - rvcos*rvcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*vz - rz*vy
           del(2) = rz*vx - rx*vz
           del(3) = rx*vy - ry*vx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*vz - del(3)*vy
           eps(2) = del(3)*vx - del(1)*vz
           eps(3) = del(1)*vy - del(2)*vx
!$acc loop seq
           do j = 1, 3
              dv = del(j)*dphidr/(vsiz*rvsin)
     &           + eps(j)*dphiddel*vpcos/(vsiz*psiz)
              call atomic_add(de(j,icloc),real( dv,md_p))
              call atomic_add(de(j,ibloc),real(-dv,md_p))
              frcx(j,i)   = frcx(j,i)   + dv
           end do
        end if
      end do

      call timer_exit( timer_torque )
      end

#ifdef USE_DETERMINISTIC_REDUCTION
      ! Included Fixed precision managment
      subroutine torquefgpu (trqvec,frcx,frcy,frcz,de,extract)
      use torquegpu_inl
      implicit none
      integer i,iipole,iglob,iloc,j
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer axetyp
      real(t_p) du,dv,dw,dt
      real(t_p) dot
      real(t_p) usiz,vsiz,wsiz
      real(t_p) psiz,rsiz,ssiz
      real(t_p) t1siz,t2siz
      real(t_p) uvsiz,uwsiz,vwsiz
      real(t_p) ursiz,ussiz
      real(t_p) vssiz,wssiz
      real(t_p) delsiz,dphiddel
      real(t_p) uvcos,uwcos,vwcos
      real(t_p) urcos,uscos
      real(t_p) vscos,wscos
      real(t_p) rwcos,wpcos
      real(t_p) upcos,vpcos
      real(t_p) rucos,rvcos
      real(t_p) ut1cos,ut2cos
      real(t_p) uvsin,uwsin,vwsin
      real(t_p) ursin,ussin
      real(t_p) vssin,wssin
      real(t_p) rwsin,wpsin
      real(t_p) upsin,vpsin
      real(t_p) rusin,rvsin
      real(t_p) ut1sin,ut2sin
      real(t_p) dphidu,dphidv,dphidw
      real(t_p) dphidr,dphids
      real(t_p) ux,uy,uz,vx,vy,vz,wx,wy,wz
      real(t_p) px,py,pz,rx,ry,rz,sx,sy,sz
      real(t_p) t1x,t1y,t1z,t2x,t2y,t2z
      real(t_p) uvx,uvy,uvz,uwx,uwy,uwz,vwx,vwy,vwz
      real(t_p) urx,ury,urz,usx,usy,usz
      real(t_p) vsx,vsy,vsz,wsx,wsy,wsz
      real(t_p) delx,dely,delz,epsx,epsy,epsz
      real(t_p) trq(3)
      real(t_p) del(3),eps(3)
      real(t_p),intent(in):: trqvec(:,:)
      real(t_p),intent(out):: frcx(:,:),frcy(:,:),frcz(:,:)
      mdyn_rtyp,intent(inout):: de(:,:)
      logical*1,intent(in):: extract !for trqvec depending on his size
c
      if (deb_Path) write(*,'(3x,a)') 'torquefgpu'
      call timer_enter( timer_torque )

c
c     zero out force components on local frame-defining atoms
c


!$acc parallel loop collapse(2) async(dir_queue)
!$acc&         present(frcx,frcy,frcz)
      do i = 1, npolelocnl
         do j = 1, 3
            frcz(j,i) = 0.0_ti_p
            frcx(j,i) = 0.0_ti_p
            frcy(j,i) = 0.0_ti_p
         end do
      end do
c
c     get the local frame type and the frame-defining atoms
c
!$acc parallel loop gang vector_length(64) async(dir_queue)
!$acc&         private(trq,del,eps)
!$acc&         present(trqvec,frcx,frcy,frcz,de)
!$acc&         present(x,y,z,ipolaxe,loc,ipole,
!$acc&  xaxis,yaxis,zaxis,poleglobnl)
      do i =1, npolelocnl
        iipole = poleglobnl(i)
        axetyp = ipolaxe(iipole)

        ia    = zaxis(iipole)
        ib    = ipole(iipole)  !iglob
        ic    = xaxis(iipole)
        ibloc = loc(ib) !iloc
        ialoc = merge(loc(ia),nbloc,ia.gt.0)
        icloc = merge(loc(ic),nbloc,ic.gt.0)
        ! trqvec can be either npolelocnl or nbloc size
        iloc  = merge(ibloc,i,extract)
        ! Check atom local Id
        if (ialoc.eq.0.or.ialoc.gt.nbloc) cycle
        if (icloc.eq.0.or.icloc.gt.nbloc) cycle

        trq(1) = trqvec(1,iloc)
        trq(2) = trqvec(2,iloc)
        trq(3) = trqvec(3,iloc)
c
c     construct the three rotation axes for the local frame
c
        ux   = x(ia) - x(ib)
        uy   = y(ia) - y(ib)
        uz   = z(ia) - z(ib)
        usiz = sqrt(ux*ux + uy*uy + uz*uz)

        if (axetyp .ne. Ax_Z_Only) then
           vx   = x(ic) - x(ib)
           vy   = y(ic) - y(ib)
           vz   = z(ic) - z(ib)
           vsiz = sqrt(vx*vx + vy*vy + vz*vz)
        else
           vx   = 1.0_ti_p
           vy   = 0.0_ti_p
           vz   = 0.0_ti_p
           vsiz = 1.0_ti_p
           dot  = ux / usiz
           if (abs(dot) .gt. 0.866_ti_p) then
              vx = 0.0_ti_p
              vy = 1.0_ti_p
           end if
        end if
        if (axetyp.eq.Ax_Z_Bisect .or. axetyp.eq.Ax_3_Fold) then
           id    = yaxis(iipole)
           idloc = merge(loc(id),nbloc,id.gt.0)
           if (idloc.eq.0.or.idloc.gt.nbloc) cycle

           wx = x(id) - x(ib)
           wy = y(id) - y(ib)
           wz = z(id) - z(ib)
        else
           wx = uy*vz - uz*vy
           wy = uz*vx - ux*vz
           wz = ux*vy - uy*vx
        end if

        wsiz  = sqrt(wx*wx + wy*wy + wz*wz)
        ux    = ux / usiz
        uy    = uy / usiz
        uz    = uz / usiz
        vx    = vx / vsiz
        vy    = vy / vsiz
        vz    = vz / vsiz
        wx    = wx / wsiz
        wy    = wy / wsiz
        wz    = wz / wsiz
c
c     find the perpendicular and angle for each pair of axes
c
        uvx   = vy*uz - vz*uy
        uvy   = vz*ux - vx*uz
        uvz   = vx*uy - vy*ux
        uvsiz = sqrt(uvx*uvx + uvy*uvy + uvz*uvz)
        uwx   = wy*uz - wz*uy
        uwy   = wz*ux - wx*uz
        uwz   = wx*uy - wy*ux
        uwsiz = sqrt(uwx*uwx + uwy*uwy + uwz*uwz)
        vwx   = wy*vz - wz*vy
        vwy   = wz*vx - wx*vz
        vwz   = wx*vy - wy*vx
        vwsiz = sqrt(vwx*vwx + vwy*vwy + vwz*vwz)

        uvx   = uvx / uvsiz
        uvy   = uvy / uvsiz
        uvz   = uvz / uvsiz
        uwx   = uwx / uwsiz
        uwy   = uwy / uwsiz
        uwz   = uwz / uwsiz
        vwx   = vwx / vwsiz
        vwy   = vwy / vwsiz
        vwz   = vwz / vwsiz
c
c     get sine and cosine of angles between the rotation axes
c
        uvcos = ux*vx + uy*vy + uz*vz
        uvsin = sqrt(abs(1.0_ti_p - uvcos*uvcos))
        uwcos = ux*wx + uy*wy + uz*wz
        uwsin = sqrt(abs(1.0_ti_p - uwcos*uwcos))
        vwcos = vx*wx + vy*wy + vz*wz
        vwsin = sqrt(abs(1.0_ti_p - vwcos*vwcos))
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
        dphidu = -trq(1)*ux - trq(2)*uy - trq(3)*uz
        dphidv = -trq(1)*vx - trq(2)*vy - trq(3)*vz
        dphidw = -trq(1)*wx - trq(2)*wy - trq(3)*wz

        if (axetyp .eq. Ax_Z_Bisect) then
c
c     build some additional axes needed for the Z-Bisect method
c
           rx    = vx + wx
           ry    = vy + wy
           rz    = vz + wz
           rsiz  = sqrt(rx*rx + ry*ry + rz*rz)
           sx    = uy*rz - uz*ry
           sy    = uz*rx - ux*rz
           sz    = ux*ry - uy*rx
           ssiz  = sqrt(sx*sx + sy*sy + sz*sz)
                 
           rsiz  = 1.0_ti_p / rsiz
           ssiz  = 1.0_ti_p / ssiz
           rx    = rx * rsiz
           ry    = ry * rsiz
           rz    = rz * rsiz
           sx    = sx * ssiz
           sy    = sy * ssiz
           sz    = sz * ssiz

           urx   = ry*uz - rz*uy
           ury   = rz*ux - rx*uz
           urz   = rx*uy - ry*ux
           ursiz = sqrt(urx*urx + ury*ury + urz*urz)
           usx   = sy*uz - sz*uy
           usy   = sz*ux - sx*uz
           usz   = sx*uy - sy*ux
           ussiz = sqrt(usx*usx + usy*usy + usz*usz)
           vsx   = sy*vz - sz*vy
           vsy   = sz*vx - sx*vz
           vsz   = sx*vy - sy*vx
           vssiz = sqrt(vsx*vsx + vsy*vsy + vsz*vsz)
           wsx   = sy*wz - sz*wy
           wsy   = sz*wx - sx*wz
           wsz   = sx*wy - sy*wx
           wssiz = sqrt(wsx*wsx + wsy*wsy + wsz*wsz)

           urx   = urx / ursiz
           ury   = ury / ursiz
           urz   = urz / ursiz
           usx   = usx / ussiz
           usy   = usy / ussiz
           usz   = usz / ussiz
           vsx   = vsx / vssiz
           vsy   = vsy / vssiz
           vsz   = vsz / vssiz
           wsx   = wsx / wssiz
           wsy   = wsy / wssiz
           wsz   = wsz / wssiz
c
c     get sine and cosine of angles between the rotation axes
c
           urcos = ux*rx + uy*ry + uz*rz
           ursin = sqrt(abs(1.0_ti_p - urcos*urcos))
           uscos = ux*sx + uy*sy + uz*sz
           ussin = sqrt(abs(1.0_ti_p - uscos*uscos))
           vscos = vx*sx + vy*sy + vz*sz
           vssin = sqrt(abs(1.0_ti_p - vscos*vscos))
           wscos = wx*sx + wy*sy + wz*sz
           wssin = sqrt(abs(1.0_ti_p - wscos*wscos))
c
c     compute the projection of v and w onto the ru-plane
c
           t1x   = vx - sx*vscos
           t1y   = vy - sy*vscos
           t1z   = vz - sz*vscos
           t2x   = wx - sx*wscos
           t2y   = wy - sy*wscos
           t2z   = wz - sz*wscos

           t1siz = sqrt(t1x*t1x+t1y*t1y+t1z*t1z)
           t2siz = sqrt(t2x*t2x+t2y*t2y+t2z*t2z)

           t1x   = t1x / t1siz
           t1y   = t1y / t1siz
           t1z   = t1z / t1siz
           t2x   = t2x / t2siz
           t2y   = t2y / t2siz
           t2z   = t2z / t2siz

           ut1cos = ux*t1x + uy*t1y + uz*t1z
           ut1sin = sqrt(abs(1.0_ti_p - ut1cos*ut1cos))
           ut2cos = ux*t2x + uy*t2y + uz*t2z
           ut2sin = sqrt(abs(1.0_ti_p - ut2cos*ut2cos))
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           dphids = -trq(1)*sx - trq(2)*sy - trq(3)*sz
        end if
c
c     force distribution for the Z-Only local coordinate method
c
        if (axetyp .eq. Ax_Z_Only) then

           du = uvx*dphidv/(usiz*uvsin)
     &        + uwx*dphidw/usiz
           call atomic_add(de(1,ialoc),( du))
           call atomic_add(de(1,ibloc),(-du))
           frcz(1,i)   = frcz(1,i)   + du
           du = uvy*dphidv/(usiz*uvsin)
     &        + uwy*dphidw/usiz
           call atomic_add(de(2,ialoc),( du))
           call atomic_add(de(2,ibloc),(-du))
           frcz(2,i)   = frcz(2,i)   + du
           du = uvz*dphidv/(usiz*uvsin)
     &        + uwz*dphidw/usiz
           call atomic_add(de(3,ialoc),( du))
           call atomic_add(de(3,ibloc),(-du))
           frcz(3,i)   = frcz(3,i)   + du
c
c     force distribution for the Z-then-X local coordinate method
c
        else if (axetyp .eq. Ax_Z_Then_X) then

           du =  uvx*dphidv/(usiz*uvsin) + uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(1,ialoc),(du))
           call atomic_add(de(1,icloc),(dv))
           call atomic_add(de(1,ibloc),(dt))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           du =  uvy*dphidv/(usiz*uvsin) + uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(2,ialoc),(du))
           call atomic_add(de(2,icloc),(dv))
           call atomic_add(de(2,ibloc),(dt))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           du =  uvz*dphidv/(usiz*uvsin) + uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(3,ialoc),(du))
           call atomic_add(de(3,icloc),(dv))
           call atomic_add(de(3,ibloc),(dt))
           frcz(3,i)   = frcz(3,i) + du
           frcx(3,i)   = frcx(3,i) + dv

c
c     force distribution for the Bisector local coordinate method
c
        else if (axetyp .eq. Ax_Bisector) then

           du =  uvx*dphidv/(usiz*uvsin) + 0.5_ti_p*uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwx*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(1,ialoc),(du))
           call atomic_add(de(1,icloc),(dv))
           call atomic_add(de(1,ibloc),(dt))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           du =  uvy*dphidv/(usiz*uvsin) + 0.5_ti_p*uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwy*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(2,ialoc),(du))
           call atomic_add(de(2,icloc),(dv))
           call atomic_add(de(2,ibloc),(dt))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           du =  uvz*dphidv/(usiz*uvsin) + 0.5_ti_p*uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwz*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(3,ialoc),(du))
           call atomic_add(de(3,icloc),(dv))
           call atomic_add(de(3,ibloc),(dt))
           frcz(3,i)   = frcz(3,i)   + du
           frcx(3,i)   = frcx(3,i)   + dv

c
c     force distribution for the Z-Bisect local coordinate method
c
        else if (axetyp .eq. Ax_Z_Bisect) then

           du = urx*dphidr/(usiz*ursin) 
     &        + usx*dphids/usiz
           dv = (vssin*sx-vscos*t1x)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sx-wscos*t2x)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(1,ialoc),(du))
           call atomic_add(de(1,icloc),(dv))
           call atomic_add(de(1,idloc),(dw))
           call atomic_add(de(1,ibloc),(dt))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           frcy(1,i)   = frcy(1,i)   + dw
           du = ury*dphidr/(usiz*ursin) + usy*dphids/usiz
           dv = (vssin*sy-vscos*t1y)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sy-wscos*t2y)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(2,ialoc),(du))
           call atomic_add(de(2,icloc),(dv))
           call atomic_add(de(2,idloc),(dw))
           call atomic_add(de(2,ibloc),(dt))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           frcy(2,i)   = frcy(2,i)   + dw
           du = urz*dphidr/(usiz*ursin) + usz*dphids/usiz
           dv = (vssin*sz-vscos*t1z)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sz-wscos*t2z)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(3,ialoc),(du))
           call atomic_add(de(3,icloc),(dv))
           call atomic_add(de(3,idloc),(dw))
           call atomic_add(de(3,ibloc),(dt))
           frcz(3,i)   = frcz(3,i)   + du
           frcx(3,i)   = frcx(3,i)   + dv
           frcy(3,i)   = frcy(3,i)   + dw

c
c     force distribution for the 3-Fold local coordinate method
c
        else if (axetyp .eq. Ax_3_Fold) then
           px     = ux + vx + wx
           py     = uy + vy + wy
           pz     = uz + vz + wz
           psiz   = sqrt(px*px + py*py + pz*pz) 
                  
           px     = px / psiz
           py     = py / psiz
           pz     = pz / psiz
                  
           wpcos  = wx*px + wy*py + wz*pz
           upcos  = ux*px + uy*py + uz*pz
           vpcos  = vx*px + vy*py + vz*pz
           wpsin  = sqrt(abs(1.0_ti_p - wpcos*wpcos))
           upsin  = sqrt(abs(1.0_ti_p - upcos*upcos))
           vpsin  = sqrt(abs(1.0_ti_p - vpcos*vpcos))
                  
           rx     = ux + vx
           ry     = uy + vy
           rz     = uz + vz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rwcos  = rx*wx + ry*wy + rz*wz
           rwsin  = sqrt(abs(1.0_ti_p - rwcos*rwcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*wz - rz*wy
           del(2) = rz*wx - rx*wz
           del(3) = rx*wy - ry*wx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*wz - del(3)*wy
           eps(2) = del(3)*wx - del(1)*wz
           eps(3) = del(1)*wy - del(2)*wx

!$acc loop seq
           do j = 1, 3
              dw = del(j)*dphidr/(wsiz*rwsin)
     &           + eps(j)*dphiddel*wpcos/(wsiz*psiz) 
              call atomic_add(de(j,idloc),( dw))
              call atomic_add(de(j,ibloc),(-dw))
              frcy(j,i)   = frcy(j,i)   + dw
           end do         
           rx     = vx + wx
           ry     = vy + wy
           rz     = vz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)

           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rucos  = rx*ux + ry*uy + rz*uz
           rusin  = sqrt(abs(1.0_ti_p - rucos*rucos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*uz - rz*uy
           del(2) = rz*ux - rx*uz
           del(3) = rx*uy - ry*ux
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*uz - del(3)*uy
           eps(2) = del(3)*ux - del(1)*uz
           eps(3) = del(1)*uy - del(2)*ux
!$acc loop seq
           do j = 1, 3
              du = del(j)*dphidr/(usiz*rusin)
     &           + eps(j)*dphiddel*upcos/(usiz*psiz)
              call atomic_add(de(j,ialoc),( du))
              call atomic_add(de(j,ibloc),(-du))
              frcz(j,i)   = frcz(j,i)   + du
           end do
           rx     = ux + wx
           ry     = uy + wy
           rz     = uz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rvcos  = rx*vx + ry*vy + rz*vz 
           rvsin  = sqrt(abs(1.0_ti_p - rvcos*rvcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*vz - rz*vy
           del(2) = rz*vx - rx*vz
           del(3) = rx*vy - ry*vx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*vz - del(3)*vy
           eps(2) = del(3)*vx - del(1)*vz
           eps(3) = del(1)*vy - del(2)*vx
!$acc loop seq
           do j = 1, 3
              dv = del(j)*dphidr/(vsiz*rvsin)
     &           + eps(j)*dphiddel*vpcos/(vsiz*psiz)
              call atomic_add(de(j,icloc),( dv))
              call atomic_add(de(j,ibloc),(-dv))
              frcx(j,i)   = frcx(j,i)   + dv
           end do
        end if
      end do

      call timer_exit( timer_torque )
      end

      subroutine torquerfgpu (natom,natombloc,polevec,ilocvec,
     &                       trqvec,de,queue)
      use torquegpu_inl
      implicit none
      integer,intent(in):: natom,natombloc
      integer,intent(in):: polevec(:),ilocvec(:)
      integer,intent(in):: queue
      real(t_p),intent(in)   ::trqvec(:,:)
      mdyn_rtyp,intent(inout):: de(:,:)

      integer i,iipole,j
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer axetyp,shade2
      real(t_p) du,dv,dw,dt
      real(t_p) usiz,vsiz,wsiz
      real(t_p) psiz,rsiz,ssiz
      real(t_p) t1siz,t2siz,dot
      real(t_p) uvsiz,uwsiz,vwsiz
      real(t_p) ursiz,ussiz
      real(t_p) vssiz,wssiz
      real(t_p) delsiz,dphiddel
      real(t_p) uvcos,uwcos,vwcos
      real(t_p) urcos,uscos
      real(t_p) vscos,wscos
      real(t_p) rwcos,wpcos
      real(t_p) upcos,vpcos
      real(t_p) rucos,rvcos
      real(t_p) ut1cos,ut2cos
      real(t_p) uvsin,uwsin,vwsin
      real(t_p) ursin,ussin
      real(t_p) vssin,wssin
      real(t_p) rwsin,wpsin
      real(t_p) upsin,vpsin
      real(t_p) rusin,rvsin
      real(t_p) ut1sin,ut2sin
      real(t_p) dphidu,dphidv,dphidw
      real(t_p) dphidr,dphids
      real(t_p) trq(3)
      real(t_p) ux,uy,uz,vx,vy,vz,wx,wy,wz
      real(t_p) px,py,pz,rx,ry,rz,sx,sy,sz
      real(t_p) t1x,t1y,t1z,t2x,t2y,t2z
      real(t_p) uvx,uvy,uvz,uwx,uwy,uwz,vwx,vwy,vwz
      real(t_p) urx,ury,urz,usx,usy,usz
      real(t_p) vsx,vsy,vsz,wsx,wsy,wsz
      real(t_p) delx,dely,delz,epsx,epsy,epsz
      real(t_p) del(3),eps(3)
      character*3 mode
c
      if (deb_Path) write(*,'(3x,a)') 'torquerfgpu'
      call timer_enter( timer_torque )
      shade2 = natombloc
c
c     get the local frame type and the frame-defining atoms
c
!$acc parallel loop gang vector_length(64) async(queue)
!$acc&         private(trq,del,eps)
!$acc&         present(polevec,ilocvec,trqvec,de)
!$acc&         present(x,y,z,ipolaxe,ipole,
!$acc&   xaxis,yaxis,zaxis)
      do i =1, natom
        iipole = polevec(i)
        axetyp = ipolaxe(iipole)

        if (axetyp .eq. Ax_None) cycle

        ia    = zaxis(iipole)
        ib    = ipole(iipole)  !iglob
        ic    = xaxis(iipole)
        id    = yaxis(iipole)
        ibloc = ilocvec(ib)  !iloc
        ialoc = merge(ilocvec(ia),shade2,ia.gt.0)
        icloc = merge(ilocvec(ic),shade2,ic.gt.0)

        ! Chech axis local Id
        if (ialoc.eq.0.or.ialoc.gt.shade2) cycle
        if (icloc.eq.0.or.icloc.gt.shade2) cycle

        trq(1) = trqvec(1,i)
        trq(2) = trqvec(2,i)
        trq(3) = trqvec(3,i)
c
c       construct the three rotation axes for the local frame
c
        ux   = x(ia) - x(ib)
        uy   = y(ia) - y(ib)
        uz   = z(ia) - z(ib)
        usiz = sqrt(ux*ux + uy*uy + uz*uz)

        if (axetyp .ne. Ax_Z_Only) then
           vx   = x(ic) - x(ib)
           vy   = y(ic) - y(ib)
           vz   = z(ic) - z(ib)
           vsiz = sqrt(vx*vx + vy*vy + vz*vz)
        else
           vx   = 1.0_ti_p
           vy   = 0.0_ti_p
           vz   = 0.0_ti_p
           vsiz = 1.0_ti_p
           dot  = ux / usiz
           if (abs(dot) .gt. 0.866_ti_p) then
              vx = 0.0_ti_p
              vy = 1.0_ti_p
           end if
        end if
        if (axetyp.eq.Ax_Z_Bisect .or. axetyp.eq.Ax_3_Fold) then
           idloc = merge(ilocvec(id),shade2,id.gt.0)
           if (idloc.eq.0.or.idloc.gt.shade2) cycle

           wx = x(id) - x(ib)
           wy = y(id) - y(ib)
           wz = z(id) - z(ib)
        else
           wx = uy*vz - uz*vy
           wy = uz*vx - ux*vz
           wz = ux*vy - uy*vx
        end if

        wsiz  = sqrt(wx*wx + wy*wy + wz*wz)
        ux    = ux / usiz
        uy    = uy / usiz
        uz    = uz / usiz
        vx    = vx / vsiz
        vy    = vy / vsiz
        vz    = vz / vsiz
        wx    = wx / wsiz
        wy    = wy / wsiz
        wz    = wz / wsiz
c
c       find the perpendicular and angle for each pair of axes
c
        uvx   = vy*uz - vz*uy
        uvy   = vz*ux - vx*uz
        uvz   = vx*uy - vy*ux
        uvsiz = sqrt(uvx*uvx + uvy*uvy + uvz*uvz)
        uwx   = wy*uz - wz*uy
        uwy   = wz*ux - wx*uz
        uwz   = wx*uy - wy*ux
        uwsiz = sqrt(uwx*uwx + uwy*uwy + uwz*uwz)
        vwx   = wy*vz - wz*vy
        vwy   = wz*vx - wx*vz
        vwz   = wx*vy - wy*vx
        vwsiz = sqrt(vwx*vwx + vwy*vwy + vwz*vwz)

        uvx   = uvx / uvsiz
        uvy   = uvy / uvsiz
        uvz   = uvz / uvsiz
        uwx   = uwx / uwsiz
        uwy   = uwy / uwsiz
        uwz   = uwz / uwsiz
        vwx   = vwx / vwsiz
        vwy   = vwy / vwsiz
        vwz   = vwz / vwsiz
c
c       get sine and cosine of angles between the rotation axes
c
        uvcos = ux*vx + uy*vy + uz*vz
        uvsin = sqrt(abs(1.0_ti_p - uvcos*uvcos))
        uwcos = ux*wx + uy*wy + uz*wz
        uwsin = sqrt(abs(1.0_ti_p - uwcos*uwcos))
        vwcos = vx*wx + vy*wy + vz*wz
        vwsin = sqrt(abs(1.0_ti_p - vwcos*vwcos))
c
c       negative of dot product of torque with unit vectors gives
c       result of infinitesimal rotation along these vectors
c
        dphidu = -trq(1)*ux - trq(2)*uy - trq(3)*uz
        dphidv = -trq(1)*vx - trq(2)*vy - trq(3)*vz
        dphidw = -trq(1)*wx - trq(2)*wy - trq(3)*wz

        if (axetyp .eq. Ax_Z_Bisect) then
c
c       build some additional axes needed for the Z-Bisect method
c
           rx    = vx + wx
           ry    = vy + wy
           rz    = vz + wz
           rsiz  = sqrt(rx*rx + ry*ry + rz*rz)
           sx    = uy*rz - uz*ry
           sy    = uz*rx - ux*rz
           sz    = ux*ry - uy*rx
           ssiz  = sqrt(sx*sx + sy*sy + sz*sz)
                 
           rsiz  = 1.0_ti_p / rsiz
           ssiz  = 1.0_ti_p / ssiz
           rx    = rx * rsiz
           ry    = ry * rsiz
           rz    = rz * rsiz
           sx    = sx * ssiz
           sy    = sy * ssiz
           sz    = sz * ssiz

           urx   = ry*uz - rz*uy
           ury   = rz*ux - rx*uz
           urz   = rx*uy - ry*ux
           ursiz = sqrt(urx*urx + ury*ury + urz*urz)
           usx   = sy*uz - sz*uy
           usy   = sz*ux - sx*uz
           usz   = sx*uy - sy*ux
           ussiz = sqrt(usx*usx + usy*usy + usz*usz)
           vsx   = sy*vz - sz*vy
           vsy   = sz*vx - sx*vz
           vsz   = sx*vy - sy*vx
           vssiz = sqrt(vsx*vsx + vsy*vsy + vsz*vsz)
           wsx   = sy*wz - sz*wy
           wsy   = sz*wx - sx*wz
           wsz   = sx*wy - sy*wx
           wssiz = sqrt(wsx*wsx + wsy*wsy + wsz*wsz)

           urx   = urx / ursiz
           ury   = ury / ursiz
           urz   = urz / ursiz
           usx   = usx / ussiz
           usy   = usy / ussiz
           usz   = usz / ussiz
           vsx   = vsx / vssiz
           vsy   = vsy / vssiz
           vsz   = vsz / vssiz
           wsx   = wsx / wssiz
           wsy   = wsy / wssiz
           wsz   = wsz / wssiz
c
c       get sine and cosine of angles between the rotation axes
c
           urcos = ux*rx + uy*ry + uz*rz
           ursin = sqrt(abs(1.0_ti_p - urcos*urcos))
           uscos = ux*sx + uy*sy + uz*sz
           ussin = sqrt(abs(1.0_ti_p - uscos*uscos))
           vscos = vx*sx + vy*sy + vz*sz
           vssin = sqrt(abs(1.0_ti_p - vscos*vscos))
           wscos = wx*sx + wy*sy + wz*sz
           wssin = sqrt(abs(1.0_ti_p - wscos*wscos))
c
c       compute the projection of v and w onto the ru-plane
c
           t1x   = vx - sx*vscos
           t1y   = vy - sy*vscos
           t1z   = vz - sz*vscos
           t2x   = wx - sx*wscos
           t2y   = wy - sy*wscos
           t2z   = wz - sz*wscos

           t1siz = sqrt(t1x*t1x+t1y*t1y+t1z*t1z)
           t2siz = sqrt(t2x*t2x+t2y*t2y+t2z*t2z)

           t1x   = t1x / t1siz
           t1y   = t1y / t1siz
           t1z   = t1z / t1siz
           t2x   = t2x / t2siz
           t2y   = t2y / t2siz
           t2z   = t2z / t2siz

           ut1cos = ux*t1x + uy*t1y + uz*t1z
           ut1sin = sqrt(abs(1.0_ti_p - ut1cos*ut1cos))
           ut2cos = ux*t2x + uy*t2y + uz*t2z
           ut2sin = sqrt(abs(1.0_ti_p - ut2cos*ut2cos))
c
c       negative of dot product of torque with unit vectors gives
c       result of infinitesimal rotation along these vectors
c
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           dphids = -trq(1)*sx - trq(2)*sy - trq(3)*sz
        end if
c
c       force distribution for the Z-Only local coordinate method
c
        if (axetyp .eq. Ax_Z_Only) then

           du = uvx*dphidv/(usiz*uvsin)
     &        + uwx*dphidw/usiz
           call atomic_add(de(1,ialoc),( du))
           call atomic_add(de(1,ibloc),(-du))
           du = uvy*dphidv/(usiz*uvsin)
     &        + uwy*dphidw/usiz
           call atomic_add(de(2,ialoc),( du))
           call atomic_add(de(2,ibloc),(-du))
           du = uvz*dphidv/(usiz*uvsin)
     &        + uwz*dphidw/usiz
           call atomic_add(de(3,ialoc),( du))
           call atomic_add(de(3,ibloc),(-du))
c
c       force distribution for the Z-then-X local coordinate method
c
        else if (axetyp .eq. Ax_Z_Then_X) then

           du =  uvx*dphidv/(usiz*uvsin) + uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(1,ialoc),(du))
           call atomic_add(de(1,icloc),(dv))
           call atomic_add(de(1,ibloc),(dt))
           du =  uvy*dphidv/(usiz*uvsin) + uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(2,ialoc),(du))
           call atomic_add(de(2,icloc),(dv))
           call atomic_add(de(2,ibloc),(dt))
           du =  uvz*dphidv/(usiz*uvsin) + uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(3,ialoc),(du))
           call atomic_add(de(3,icloc),(dv))
           call atomic_add(de(3,ibloc),(dt))
c
c       force distribution for the Bisector local coordinate method
c
        else if (axetyp .eq. Ax_Bisector) then

           du =  uvx*dphidv/(usiz*uvsin) + 0.5_ti_p*uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwx*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(1,ialoc),(du))
           call atomic_add(de(1,icloc),(dv))
           call atomic_add(de(1,ibloc),(dt))
           du =  uvy*dphidv/(usiz*uvsin) + 0.5_ti_p*uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwy*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(2,ialoc),(du))
           call atomic_add(de(2,icloc),(dv))
           call atomic_add(de(2,ibloc),(dt))
           du =  uvz*dphidv/(usiz*uvsin) + 0.5_ti_p*uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwz*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(3,ialoc),(du))
           call atomic_add(de(3,icloc),(dv))
           call atomic_add(de(3,ibloc),(dt))
c
c       force distribution for the Z-Bisect local coordinate method
c
        else if (axetyp .eq. Ax_Z_Bisect) then

           du = urx*dphidr/(usiz*ursin) 
     &        + usx*dphids/usiz
           dv = (vssin*sx-vscos*t1x)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sx-wscos*t2x)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(1,ialoc),(du))
           call atomic_add(de(1,icloc),(dv))
           call atomic_add(de(1,idloc),(dw))
           call atomic_add(de(1,ibloc),(dt))
           du = ury*dphidr/(usiz*ursin) + usy*dphids/usiz
           dv = (vssin*sy-vscos*t1y)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sy-wscos*t2y)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(2,ialoc),(du))
           call atomic_add(de(2,icloc),(dv))
           call atomic_add(de(2,idloc),(dw))
           call atomic_add(de(2,ibloc),(dt))
           du = urz*dphidr/(usiz*ursin) + usz*dphids/usiz
           dv = (vssin*sz-vscos*t1z)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sz-wscos*t2z)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(3,ialoc),(du))
           call atomic_add(de(3,icloc),(dv))
           call atomic_add(de(3,idloc),(dw))
           call atomic_add(de(3,ibloc),(dt))
c
c       force distribution for the 3-Fold local coordinate method
c
        else if (axetyp .eq. Ax_3_Fold) then
           px     = ux + vx + wx
           py     = uy + vy + wy
           pz     = uz + vz + wz
           psiz   = sqrt(px*px + py*py + pz*pz) 
                  
           px     = px / psiz
           py     = py / psiz
           pz     = pz / psiz
                  
           wpcos  = wx*px + wy*py + wz*pz
           upcos  = ux*px + uy*py + uz*pz
           vpcos  = vx*px + vy*py + vz*pz
           wpsin  = sqrt(abs(1.0_ti_p - wpcos*wpcos))
           upsin  = sqrt(abs(1.0_ti_p - upcos*upcos))
           vpsin  = sqrt(abs(1.0_ti_p - vpcos*vpcos))
                  
           rx     = ux + vx
           ry     = uy + vy
           rz     = uz + vz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rwcos  = rx*wx + ry*wy + rz*wz
           rwsin  = sqrt(abs(1.0_ti_p - rwcos*rwcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*wz - rz*wy
           del(2) = rz*wx - rx*wz
           del(3) = rx*wy - ry*wx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*wz - del(3)*wy
           eps(2) = del(3)*wx - del(1)*wz
           eps(3) = del(1)*wy - del(2)*wx
!$acc loop seq
           do j = 1, 3
              dw = del(j)*dphidr/(wsiz*rwsin)
     &           + eps(j)*dphiddel*wpcos/(wsiz*psiz) 
              call atomic_add(de(j,idloc),( dw))
              call atomic_add(de(j,ibloc),(-dw))
           end do         
           rx     = vx + wx
           ry     = vy + wy
           rz     = vz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)

           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rucos  = rx*ux + ry*uy + rz*uz
           rusin  = sqrt(abs(1.0_ti_p - rucos*rucos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*uz - rz*uy
           del(2) = rz*ux - rx*uz
           del(3) = rx*uy - ry*ux
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*uz - del(3)*uy
           eps(2) = del(3)*ux - del(1)*uz
           eps(3) = del(1)*uy - del(2)*ux
!$acc loop seq
           do j = 1, 3
              du = del(j)*dphidr/(usiz*rusin)
     &           + eps(j)*dphiddel*upcos/(usiz*psiz)
              call atomic_add(de(j,ialoc),( du))
              call atomic_add(de(j,ibloc),(-du))
           end do
           rx     = ux + wx
           ry     = uy + wy
           rz     = uz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rvcos  = rx*vx + ry*vy + rz*vz 
           rvsin  = sqrt(abs(1.0_ti_p - rvcos*rvcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*vz - rz*vy
           del(2) = rz*vx - rx*vz
           del(3) = rx*vy - ry*vx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*vz - del(3)*vy
           eps(2) = del(3)*vx - del(1)*vz
           eps(3) = del(1)*vy - del(2)*vx
!$acc loop seq
           do j = 1, 3
              dv = del(j)*dphidr/(vsiz*rvsin)
     &           + eps(j)*dphiddel*vpcos/(vsiz*psiz)
              call atomic_add(de(j,icloc),( dv))
              call atomic_add(de(j,ibloc),(-dv))
           end do
        end if
      end do

      call timer_exit( timer_torque )
      end
#endif
c
c     This routine is different form the first one only by his output
c
      subroutine torquergpu (natom,natombloc,polevec,ilocvec,
     &                       trqvec,de,queue)
      use torquegpu_inl
      implicit none
      integer i,iipole,j
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer shade2
      integer,intent(in):: natom,natombloc
      integer,intent(in):: polevec(:),ilocvec(:)
      integer,intent(in):: queue
      real(t_p) du,dv,dw,dt
      real(t_p) usiz,vsiz,wsiz
      real(t_p) psiz,rsiz,ssiz
      real(t_p) t1siz,t2siz,dot
      real(t_p) uvsiz,uwsiz,vwsiz
      real(t_p) ursiz,ussiz
      real(t_p) vssiz,wssiz
      real(t_p) delsiz,dphiddel
      real(t_p) uvcos,uwcos,vwcos
      real(t_p) urcos,uscos
      real(t_p) vscos,wscos
      real(t_p) rwcos,wpcos
      real(t_p) upcos,vpcos
      real(t_p) rucos,rvcos
      real(t_p) ut1cos,ut2cos
      real(t_p) uvsin,uwsin,vwsin
      real(t_p) ursin,ussin
      real(t_p) vssin,wssin
      real(t_p) rwsin,wpsin
      real(t_p) upsin,vpsin
      real(t_p) rusin,rvsin
      real(t_p) ut1sin,ut2sin
      real(t_p) dphidu,dphidv,dphidw
      real(t_p) dphidr,dphids
      real(t_p) trq(3)
      real(t_p) ux,uy,uz,vx,vy,vz,wx,wy,wz
      real(t_p) px,py,pz,rx,ry,rz,sx,sy,sz
      real(t_p) t1x,t1y,t1z,t2x,t2y,t2z
      real(t_p) uvx,uvy,uvz,uwx,uwy,uwz,vwx,vwy,vwz
      real(t_p) urx,ury,urz,usx,usy,usz
      real(t_p) vsx,vsy,vsz,wsx,wsy,wsz
      real(t_p) delx,dely,delz,epsx,epsy,epsz
      real(t_p) del(3),eps(3)
      real(t_p),intent(in)   ::trqvec(:,:)
      real(md_p),intent(inout):: de(:,:)
      character*3 mode
      integer axetyp
c
      if (deb_Path) write(*,'(3x,a)') 'torquergpu'
      call timer_enter( timer_torque )
      shade2 = natombloc

c
c     get the local frame type and the frame-defining atoms
c
!$acc parallel loop gang vector_length(64) async(queue)
!$acc&         private(trq,del,eps)
!$acc&         present(polevec,ilocvec,trqvec,de)
!$acc&         present(x,y,z,ipolaxe,ipole,
!$acc&   xaxis,yaxis,zaxis)
      do i =1, natom
        iipole = polevec(i)
        axetyp = ipolaxe(iipole)

        if (axetyp .eq. Ax_None) cycle
        ia    = zaxis(iipole)
        ib    = ipole(iipole)  !iglob
        ic    = xaxis(iipole)
        id    = yaxis(iipole)
        ibloc = ilocvec(ib)    !iloc
        ialoc = merge(ilocvec(ia),shade2,ia.gt.0)
        icloc = merge(ilocvec(ic),shade2,ic.gt.0)        
        ! Chech axis local Id
        if (ialoc.eq.0.or.ialoc.gt.shade2) cycle;
        if (icloc.eq.0.or.icloc.gt.shade2) cycle;

        trq(1) = trqvec(1,i)
        trq(2) = trqvec(2,i)
        trq(3) = trqvec(3,i)
c
c       construct the three rotation axes for the local frame
c
        ux   = x(ia) - x(ib)
        uy   = y(ia) - y(ib)
        uz   = z(ia) - z(ib)
        usiz = sqrt(ux*ux + uy*uy + uz*uz)

        if (axetyp .ne. Ax_Z_Only) then
           vx   = x(ic) - x(ib)
           vy   = y(ic) - y(ib)
           vz   = z(ic) - z(ib)
           vsiz = sqrt(vx*vx + vy*vy + vz*vz)
        else
           vx   = 1.0_ti_p
           vy   = 0.0_ti_p
           vz   = 0.0_ti_p
           vsiz = 1.0_ti_p
           dot  = ux / usiz
           if (abs(dot) .gt. 0.866_ti_p) then
              vx = 0.0_ti_p
              vy = 1.0_ti_p
           end if
        end if
        if (axetyp.eq.Ax_Z_Bisect .or. axetyp.eq.Ax_3_Fold) then
           idloc = merge(ilocvec(id),shade2,id.gt.0)
           if (idloc.eq.0.or.idloc.gt.shade2) cycle

           wx = x(id) - x(ib)
           wy = y(id) - y(ib)
           wz = z(id) - z(ib)
        else
           wx = uy*vz - uz*vy
           wy = uz*vx - ux*vz
           wz = ux*vy - uy*vx
        end if

        wsiz  = sqrt(wx*wx + wy*wy + wz*wz)
        ux    = ux / usiz
        uy    = uy / usiz
        uz    = uz / usiz
        vx    = vx / vsiz
        vy    = vy / vsiz
        vz    = vz / vsiz
        wx    = wx / wsiz
        wy    = wy / wsiz
        wz    = wz / wsiz
c
c       find the perpendicular and angle for each pair of axes
c
        uvx   = vy*uz - vz*uy
        uvy   = vz*ux - vx*uz
        uvz   = vx*uy - vy*ux
        uvsiz = sqrt(uvx*uvx + uvy*uvy + uvz*uvz)
        uwx   = wy*uz - wz*uy
        uwy   = wz*ux - wx*uz
        uwz   = wx*uy - wy*ux
        uwsiz = sqrt(uwx*uwx + uwy*uwy + uwz*uwz)
        vwx   = wy*vz - wz*vy
        vwy   = wz*vx - wx*vz
        vwz   = wx*vy - wy*vx
        vwsiz = sqrt(vwx*vwx + vwy*vwy + vwz*vwz)

        uvx   = uvx / uvsiz
        uvy   = uvy / uvsiz
        uvz   = uvz / uvsiz
        uwx   = uwx / uwsiz
        uwy   = uwy / uwsiz
        uwz   = uwz / uwsiz
        vwx   = vwx / vwsiz
        vwy   = vwy / vwsiz
        vwz   = vwz / vwsiz
c
c       get sine and cosine of angles between the rotation axes
c
        uvcos = ux*vx + uy*vy + uz*vz
        uvsin = sqrt(abs(1.0_ti_p - uvcos*uvcos))
        uwcos = ux*wx + uy*wy + uz*wz
        uwsin = sqrt(abs(1.0_ti_p - uwcos*uwcos))
        vwcos = vx*wx + vy*wy + vz*wz
        vwsin = sqrt(abs(1.0_ti_p - vwcos*vwcos))
c
c       negative of dot product of torque with unit vectors gives
c       result of infinitesimal rotation along these vectors
c
        dphidu = -trq(1)*ux - trq(2)*uy - trq(3)*uz
        dphidv = -trq(1)*vx - trq(2)*vy - trq(3)*vz
        dphidw = -trq(1)*wx - trq(2)*wy - trq(3)*wz

        if (axetyp .eq. Ax_Z_Bisect) then
c
c       build some additional axes needed for the Z-Bisect method
c
           rx    = vx + wx
           ry    = vy + wy
           rz    = vz + wz
           rsiz  = sqrt(rx*rx + ry*ry + rz*rz)
           sx    = uy*rz - uz*ry
           sy    = uz*rx - ux*rz
           sz    = ux*ry - uy*rx
           ssiz  = sqrt(sx*sx + sy*sy + sz*sz)
                 
           rsiz  = 1.0_ti_p / rsiz
           ssiz  = 1.0_ti_p / ssiz
           rx    = rx * rsiz
           ry    = ry * rsiz
           rz    = rz * rsiz
           sx    = sx * ssiz
           sy    = sy * ssiz
           sz    = sz * ssiz

           urx   = ry*uz - rz*uy
           ury   = rz*ux - rx*uz
           urz   = rx*uy - ry*ux
           ursiz = sqrt(urx*urx + ury*ury + urz*urz)
           usx   = sy*uz - sz*uy
           usy   = sz*ux - sx*uz
           usz   = sx*uy - sy*ux
           ussiz = sqrt(usx*usx + usy*usy + usz*usz)
           vsx   = sy*vz - sz*vy
           vsy   = sz*vx - sx*vz
           vsz   = sx*vy - sy*vx
           vssiz = sqrt(vsx*vsx + vsy*vsy + vsz*vsz)
           wsx   = sy*wz - sz*wy
           wsy   = sz*wx - sx*wz
           wsz   = sx*wy - sy*wx
           wssiz = sqrt(wsx*wsx + wsy*wsy + wsz*wsz)

           urx   = urx / ursiz
           ury   = ury / ursiz
           urz   = urz / ursiz
           usx   = usx / ussiz
           usy   = usy / ussiz
           usz   = usz / ussiz
           vsx   = vsx / vssiz
           vsy   = vsy / vssiz
           vsz   = vsz / vssiz
           wsx   = wsx / wssiz
           wsy   = wsy / wssiz
           wsz   = wsz / wssiz
c
c       get sine and cosine of angles between the rotation axes
c
           urcos = ux*rx + uy*ry + uz*rz
           ursin = sqrt(abs(1.0_ti_p - urcos*urcos))
           uscos = ux*sx + uy*sy + uz*sz
           ussin = sqrt(abs(1.0_ti_p - uscos*uscos))
           vscos = vx*sx + vy*sy + vz*sz
           vssin = sqrt(abs(1.0_ti_p - vscos*vscos))
           wscos = wx*sx + wy*sy + wz*sz
           wssin = sqrt(abs(1.0_ti_p - wscos*wscos))
c
c       compute the projection of v and w onto the ru-plane
c
           t1x   = vx - sx*vscos
           t1y   = vy - sy*vscos
           t1z   = vz - sz*vscos
           t2x   = wx - sx*wscos
           t2y   = wy - sy*wscos
           t2z   = wz - sz*wscos

           t1siz = sqrt(t1x*t1x+t1y*t1y+t1z*t1z)
           t2siz = sqrt(t2x*t2x+t2y*t2y+t2z*t2z)

           t1x   = t1x / t1siz
           t1y   = t1y / t1siz
           t1z   = t1z / t1siz
           t2x   = t2x / t2siz
           t2y   = t2y / t2siz
           t2z   = t2z / t2siz

           ut1cos = ux*t1x + uy*t1y + uz*t1z
           ut1sin = sqrt(abs(1.0_ti_p - ut1cos*ut1cos))
           ut2cos = ux*t2x + uy*t2y + uz*t2z
           ut2sin = sqrt(abs(1.0_ti_p - ut2cos*ut2cos))
c
c       negative of dot product of torque with unit vectors gives
c       result of infinitesimal rotation along these vectors
c
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           dphids = -trq(1)*sx - trq(2)*sy - trq(3)*sz
        end if
c
c       force distribution for the Z-Only local coordinate method
c
        if (axetyp .eq. Ax_Z_Only) then

           du = uvx*dphidv/(usiz*uvsin)
     &        + uwx*dphidw/usiz
           call atomic_add(de(1,ialoc),real( du,md_p))
           call atomic_add(de(1,ibloc),real(-du,md_p))
           du = uvy*dphidv/(usiz*uvsin)
     &        + uwy*dphidw/usiz
           call atomic_add(de(2,ialoc),real( du,md_p))
           call atomic_add(de(2,ibloc),real(-du,md_p))
           du = uvz*dphidv/(usiz*uvsin)
     &        + uwz*dphidw/usiz
           call atomic_add(de(3,ialoc),real( du,md_p))
           call atomic_add(de(3,ibloc),real(-du,md_p))
c
c       force distribution for the Z-then-X local coordinate method
c
        else if (axetyp .eq. Ax_Z_Then_X) then

           du =  uvx*dphidv/(usiz*uvsin) + uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           du =  uvy*dphidv/(usiz*uvsin) + uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           du =  uvz*dphidv/(usiz*uvsin) + uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
c
c       force distribution for the Bisector local coordinate method
c
        else if (axetyp .eq. Ax_Bisector) then

           du =  uvx*dphidv/(usiz*uvsin) + 0.5_ti_p*uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwx*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           du =  uvy*dphidv/(usiz*uvsin) + 0.5_ti_p*uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwy*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           du =  uvz*dphidv/(usiz*uvsin) + 0.5_ti_p*uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwz*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
c
c       force distribution for the Z-Bisect local coordinate method
c
        else if (axetyp .eq. Ax_Z_Bisect) then

           du = urx*dphidr/(usiz*ursin) 
     &        + usx*dphids/usiz
           dv = (vssin*sx-vscos*t1x)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sx-wscos*t2x)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,idloc),real(dw,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           du = ury*dphidr/(usiz*ursin) + usy*dphids/usiz
           dv = (vssin*sy-vscos*t1y)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sy-wscos*t2y)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,idloc),real(dw,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           du = urz*dphidr/(usiz*ursin) + usz*dphids/usiz
           dv = (vssin*sz-vscos*t1z)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sz-wscos*t2z)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,idloc),real(dw,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
c
c       force distribution for the 3-Fold local coordinate method
c
        else if (axetyp .eq. Ax_3_Fold) then
           px     = ux + vx + wx
           py     = uy + vy + wy
           pz     = uz + vz + wz
           psiz   = sqrt(px*px + py*py + pz*pz) 
                  
           px     = px / psiz
           py     = py / psiz
           pz     = pz / psiz
                  
           wpcos  = wx*px + wy*py + wz*pz
           upcos  = ux*px + uy*py + uz*pz
           vpcos  = vx*px + vy*py + vz*pz
           wpsin  = sqrt(abs(1.0_ti_p - wpcos*wpcos))
           upsin  = sqrt(abs(1.0_ti_p - upcos*upcos))
           vpsin  = sqrt(abs(1.0_ti_p - vpcos*vpcos))
                  
           rx     = ux + vx
           ry     = uy + vy
           rz     = uz + vz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rwcos  = rx*wx + ry*wy + rz*wz
           rwsin  = sqrt(abs(1.0_ti_p - rwcos*rwcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*wz - rz*wy
           del(2) = rz*wx - rx*wz
           del(3) = rx*wy - ry*wx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*wz - del(3)*wy
           eps(2) = del(3)*wx - del(1)*wz
           eps(3) = del(1)*wy - del(2)*wx
!$acc loop seq
           do j = 1, 3
              dw = del(j)*dphidr/(wsiz*rwsin)
     &           + eps(j)*dphiddel*wpcos/(wsiz*psiz) 
              call atomic_add(de(j,idloc),real( dw,md_p))
              call atomic_add(de(j,ibloc),real(-dw,md_p))
           end do         
           rx     = vx + wx
           ry     = vy + wy
           rz     = vz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)

           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rucos  = rx*ux + ry*uy + rz*uz
           rusin  = sqrt(abs(1.0_ti_p - rucos*rucos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*uz - rz*uy
           del(2) = rz*ux - rx*uz
           del(3) = rx*uy - ry*ux
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*uz - del(3)*uy
           eps(2) = del(3)*ux - del(1)*uz
           eps(3) = del(1)*uy - del(2)*ux
!$acc loop seq
           do j = 1, 3
              du = del(j)*dphidr/(usiz*rusin)
     &           + eps(j)*dphiddel*upcos/(usiz*psiz)
              call atomic_add(de(j,ialoc),real( du,md_p))
              call atomic_add(de(j,ibloc),real(-du,md_p))
           end do
           rx     = ux + wx
           ry     = uy + wy
           rz     = uz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rvcos  = rx*vx + ry*vy + rz*vz 
           rvsin  = sqrt(abs(1.0_ti_p - rvcos*rvcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*vz - rz*vy
           del(2) = rz*vx - rx*vz
           del(3) = rx*vy - ry*vx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*vz - del(3)*vy
           eps(2) = del(3)*vx - del(1)*vz
           eps(3) = del(1)*vy - del(2)*vx
!$acc loop seq
           do j = 1, 3
              dv = del(j)*dphidr/(vsiz*rvsin)
     &           + eps(j)*dphiddel*vpcos/(vsiz*psiz)
              call atomic_add(de(j,icloc),real( dv,md_p))
              call atomic_add(de(j,ibloc),real(-dv,md_p))
           end do
        end if
      end do

      call timer_exit( timer_torque )
      end


      subroutine torquegpu_group (trqvec,frcx,frcy,frcz,de)
      use torquegpu_inl
      use group
      use mpole, only: pollist
      implicit none
      real(t_p),intent(in):: trqvec(:,:)
      real(t_p),intent(out):: frcx(:,:),frcy(:,:),frcz(:,:)
      mdyn_rtyp,intent(inout):: de(:,:)

      integer i,iipole,iglob,iloc,j,ca
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer :: iglobgroup,iploc
      integer axetyp,iipolegroup
      real(t_p) du,dv,dw,dt
      real(t_p) dot
      real(t_p) usiz,vsiz,wsiz
      real(t_p) psiz,rsiz,ssiz
      real(t_p) t1siz,t2siz
      real(t_p) uvsiz,uwsiz,vwsiz
      real(t_p) ursiz,ussiz
      real(t_p) vssiz,wssiz
      real(t_p) delsiz,dphiddel
      real(t_p) uvcos,uwcos,vwcos
      real(t_p) urcos,uscos
      real(t_p) vscos,wscos
      real(t_p) rwcos,wpcos
      real(t_p) upcos,vpcos
      real(t_p) rucos,rvcos
      real(t_p) ut1cos,ut2cos
      real(t_p) uvsin,uwsin,vwsin
      real(t_p) ursin,ussin
      real(t_p) vssin,wssin
      real(t_p) rwsin,wpsin
      real(t_p) upsin,vpsin
      real(t_p) rusin,rvsin
      real(t_p) ut1sin,ut2sin
      real(t_p) dphidu,dphidv,dphidw
      real(t_p) dphidr,dphids
      real(t_p) ux,uy,uz,vx,vy,vz,wx,wy,wz
      real(t_p) px,py,pz,rx,ry,rz,sx,sy,sz
      real(t_p) t1x,t1y,t1z,t2x,t2y,t2z
      real(t_p) uvx,uvy,uvz,uwx,uwy,uwz,vwx,vwy,vwz
      real(t_p) urx,ury,urz,usx,usy,usz
      real(t_p) vsx,vsy,vsz,wsx,wsy,wsz
      real(t_p) delx,dely,delz,epsx,epsy,epsz
      real(t_p) trq(3)
      real(t_p) del(3),eps(3)
!$acc routine(distprocpart1)
c
      if (deb_Path) write(*,'(3x,a)') 'torquegpu_group'
      call timer_enter( timer_torque )

      if (tinkerdebug.gt.0) then
         j=0

!$acc parallel loop gang vector_length(128) async(dir_queue)
!$acc&         present(ipolaxe,loc,ipole,
!$acc&  xaxis,yaxis,zaxis,globglobgroup,ipolegroup,pollist,
!$acc&  locgroup,loclocgroup) copyin(j)
      do i =1, npolegroup
        iglob = globglobgroup(ipolegroup(i))
        iipole = pollist(iglob)
        axetyp = ipolaxe(iipole)

        if (axetyp .eq. Ax_None) cycle

        ia    = zaxis(iipole)
        ib    = ipole(iipole)  !iglob
        ic    = xaxis(iipole)
        id    = yaxis(iipole)
        ibloc = locgroup(loclocgroup(ib)) !iloc
        ialoc = -1; idloc=-1; icloc=-1;
        if (ia.gt.0) then; ialoc = locgroup(loclocgroup(ia))
          else; ialoc=ibloc; endif
        if (ic.gt.0) then; icloc = locgroup(loclocgroup(ic))
          else; icloc=ibloc; endif
        if (id.gt.0) then; idloc = locgroup(loclocgroup(id))
          else; idloc=ibloc; endif

        if (ibloc.eq.0.or.ibloc.gt.natgroup) then
!$acc atomic capture
           j = j+1
           ca = j
!$acc end atomic
           if(ca.le.10) print*,'WARNING!!! torquegpu_group Outside loc'
           print*,'b ',iipole,ib,ibloc,rank
        end if
        if (ialoc.eq.0.or.ialoc.gt.natgroup) then
!$acc atomic capture
           j = j+1
           ca = j
!$acc end atomic
           if(ca.le.10) print*,'WARNING!!! torquegpu_group Outside loc'
           call distprocpart1(ib,rank,uz,.true.,x,y,z)
           print*,'a ',iipole,ia,ialoc,ib,ibloc,nbloc,uz,rank
        end if
        if (icloc.eq.0.or.icloc.gt.natgroup) then
           call distprocpart1(ib,rank,ux,.true.,x,y,z)
           print*,'c ',iipole,ic,icloc,ib,ibloc,nbloc,ux,rank
        end if
        if (idloc.eq.0.or.idloc.gt.natgroup) then
           call distprocpart1(id,rank,uy,.true.,x,y,z)
           print*,'d ',iipole,id,idloc,ib,ibloc,nbloc,uy,rank
        end if
      end do

      end if
c
c     zero out force components on local frame-defining atoms
c
!$acc parallel loop collapse(2) async(dir_queue)
!$acc&         present(frcx,frcy,frcz)
      do i = 1, npolegroup
         do j = 1, 3
            frcz(j,i) = 0.0_ti_p
            frcx(j,i) = 0.0_ti_p
            frcy(j,i) = 0.0_ti_p
         end do
      end do
c
c     get the local frame type and the frame-defining atoms
c
!$acc parallel loop gang vector_length(32) async(dir_queue)
!$acc&         private(trq,del,eps)
!$acc&         present(trqvec,frcx,frcy,frcz,de)
!$acc&         present(x,y,z,ipolaxe,loc,globglobgroup,
!$acc&  xaxis,yaxis,zaxis,ipolegroup,pollist,ipole,
!$acc&  locgroup,loclocgroup)
      do i =1, npolegroup
        iglobgroup=ipolegroup(i)
        iglob = globglobgroup(iglobgroup)
        iipole = pollist(iglob)
        axetyp = ipolaxe(iipole)

        if (axetyp .eq. Ax_None) cycle

        ia    = zaxis(iipole)
        ib    = ipole(iipole)  !iglob
        ic    = xaxis(iipole)
        ibloc = locgroup(loclocgroup(ib)) !iloc
        ialoc = merge(locgroup(loclocgroup(ia)),natgroup,ia.gt.0)
        icloc = merge(locgroup(loclocgroup(ic)),natgroup,ic.gt.0)
        iloc=ibloc
        iploc=polelocgroup(i)

        ! Chech axis local Id
        if (ialoc.eq.0.or.ialoc.gt.natgroup) cycle
        if (icloc.eq.0.or.icloc.gt.natgroup) cycle

        trq(1) = trqvec(1,iploc)
        trq(2) = trqvec(2,iploc)
        trq(3) = trqvec(3,iploc)
c
c     construct the three rotation axes for the local frame
c
        ux   = x(ia) - x(ib)
        uy   = y(ia) - y(ib)
        uz   = z(ia) - z(ib)
        usiz = sqrt(ux*ux + uy*uy + uz*uz)

        if (axetyp .ne. Ax_Z_Only) then
           vx   = x(ic) - x(ib)
           vy   = y(ic) - y(ib)
           vz   = z(ic) - z(ib)
           vsiz = sqrt(vx*vx + vy*vy + vz*vz)
        else
           vx   = 1.0_ti_p
           vy   = 0.0_ti_p
           vz   = 0.0_ti_p
           vsiz = 1.0_ti_p
           dot  = ux / usiz
           if (abs(dot) .gt. 0.866_ti_p) then
              vx = 0.0_ti_p
              vy = 1.0_ti_p
           end if
        end if
        if (axetyp.eq.Ax_Z_Bisect .or. axetyp.eq.Ax_3_Fold) then
           id    = yaxis(iipole)
           idloc = merge(locgroup(loclocgroup(id)),natgroup,id.gt.0)
           if (idloc.eq.0.or.idloc.gt.natgroup) cycle
           wx = x(id) - x(ib)
           wy = y(id) - y(ib)
           wz = z(id) - z(ib)
        else
           wx = uy*vz - uz*vy
           wy = uz*vx - ux*vz
           wz = ux*vy - uy*vx
        end if

        wsiz  = sqrt(wx*wx + wy*wy + wz*wz)
        ux    = ux / usiz
        uy    = uy / usiz
        uz    = uz / usiz
        vx    = vx / vsiz
        vy    = vy / vsiz
        vz    = vz / vsiz
        wx    = wx / wsiz
        wy    = wy / wsiz
        wz    = wz / wsiz
c
c     find the perpendicular and angle for each pair of axes
c
        uvx   = vy*uz - vz*uy
        uvy   = vz*ux - vx*uz
        uvz   = vx*uy - vy*ux
        uvsiz = sqrt(uvx*uvx + uvy*uvy + uvz*uvz)
        uwx   = wy*uz - wz*uy
        uwy   = wz*ux - wx*uz
        uwz   = wx*uy - wy*ux
        uwsiz = sqrt(uwx*uwx + uwy*uwy + uwz*uwz)
        vwx   = wy*vz - wz*vy
        vwy   = wz*vx - wx*vz
        vwz   = wx*vy - wy*vx
        vwsiz = sqrt(vwx*vwx + vwy*vwy + vwz*vwz)

        uvx   = uvx / uvsiz
        uvy   = uvy / uvsiz
        uvz   = uvz / uvsiz
        uwx   = uwx / uwsiz
        uwy   = uwy / uwsiz
        uwz   = uwz / uwsiz
        vwx   = vwx / vwsiz
        vwy   = vwy / vwsiz
        vwz   = vwz / vwsiz
c
c     get sine and cosine of angles between the rotation axes
c
        uvcos = ux*vx + uy*vy + uz*vz
        uvsin = sqrt(abs(1.0_ti_p - uvcos*uvcos))
        uwcos = ux*wx + uy*wy + uz*wz
        uwsin = sqrt(abs(1.0_ti_p - uwcos*uwcos))
        vwcos = vx*wx + vy*wy + vz*wz
        vwsin = sqrt(abs(1.0_ti_p - vwcos*vwcos))
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
        dphidu = -trq(1)*ux - trq(2)*uy - trq(3)*uz
        dphidv = -trq(1)*vx - trq(2)*vy - trq(3)*vz
        dphidw = -trq(1)*wx - trq(2)*wy - trq(3)*wz

        if (axetyp .eq. Ax_Z_Bisect) then
c
c     build some additional axes needed for the Z-Bisect method
c
           rx    = vx + wx
           ry    = vy + wy
           rz    = vz + wz
           rsiz  = sqrt(rx*rx + ry*ry + rz*rz)
           sx    = uy*rz - uz*ry
           sy    = uz*rx - ux*rz
           sz    = ux*ry - uy*rx
           ssiz  = sqrt(sx*sx + sy*sy + sz*sz)
                 
           rsiz  = 1.0_ti_p / rsiz
           ssiz  = 1.0_ti_p / ssiz
           rx    = rx * rsiz
           ry    = ry * rsiz
           rz    = rz * rsiz
           sx    = sx * ssiz
           sy    = sy * ssiz
           sz    = sz * ssiz

           urx   = ry*uz - rz*uy
           ury   = rz*ux - rx*uz
           urz   = rx*uy - ry*ux
           ursiz = sqrt(urx*urx + ury*ury + urz*urz)
           usx   = sy*uz - sz*uy
           usy   = sz*ux - sx*uz
           usz   = sx*uy - sy*ux
           ussiz = sqrt(usx*usx + usy*usy + usz*usz)
           vsx   = sy*vz - sz*vy
           vsy   = sz*vx - sx*vz
           vsz   = sx*vy - sy*vx
           vssiz = sqrt(vsx*vsx + vsy*vsy + vsz*vsz)
           wsx   = sy*wz - sz*wy
           wsy   = sz*wx - sx*wz
           wsz   = sx*wy - sy*wx
           wssiz = sqrt(wsx*wsx + wsy*wsy + wsz*wsz)

           urx   = urx / ursiz
           ury   = ury / ursiz
           urz   = urz / ursiz
           usx   = usx / ussiz
           usy   = usy / ussiz
           usz   = usz / ussiz
           vsx   = vsx / vssiz
           vsy   = vsy / vssiz
           vsz   = vsz / vssiz
           wsx   = wsx / wssiz
           wsy   = wsy / wssiz
           wsz   = wsz / wssiz
c
c     get sine and cosine of angles between the rotation axes
c
           urcos = ux*rx + uy*ry + uz*rz
           ursin = sqrt(abs(1.0_ti_p - urcos*urcos))
           uscos = ux*sx + uy*sy + uz*sz
           ussin = sqrt(abs(1.0_ti_p - uscos*uscos))
           vscos = vx*sx + vy*sy + vz*sz
           vssin = sqrt(abs(1.0_ti_p - vscos*vscos))
           wscos = wx*sx + wy*sy + wz*sz
           wssin = sqrt(abs(1.0_ti_p - wscos*wscos))
c
c     compute the projection of v and w onto the ru-plane
c
           t1x   = vx - sx*vscos
           t1y   = vy - sy*vscos
           t1z   = vz - sz*vscos
           t2x   = wx - sx*wscos
           t2y   = wy - sy*wscos
           t2z   = wz - sz*wscos

           t1siz = sqrt(t1x*t1x+t1y*t1y+t1z*t1z)
           t2siz = sqrt(t2x*t2x+t2y*t2y+t2z*t2z)

           t1x   = t1x / t1siz
           t1y   = t1y / t1siz
           t1z   = t1z / t1siz
           t2x   = t2x / t2siz
           t2y   = t2y / t2siz
           t2z   = t2z / t2siz

           ut1cos = ux*t1x + uy*t1y + uz*t1z
           ut1sin = sqrt(abs(1.0_ti_p - ut1cos*ut1cos))
           ut2cos = ux*t2x + uy*t2y + uz*t2z
           ut2sin = sqrt(abs(1.0_ti_p - ut2cos*ut2cos))
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           dphids = -trq(1)*sx - trq(2)*sy - trq(3)*sz
        end if
c
c     force distribution for the Z-Only local coordinate method
c
        if (axetyp .eq. Ax_Z_Only) then

           du = uvx*dphidv/(usiz*uvsin)
     &        + uwx*dphidw/usiz
           call atomic_add(de(1,ialoc),real( du,md_p))
           call atomic_add(de(1,ibloc),real(-du,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           du = uvy*dphidv/(usiz*uvsin)
     &        + uwy*dphidw/usiz
           call atomic_add(de(2,ialoc),real( du,md_p))
           call atomic_add(de(2,ibloc),real(-du,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           du = uvz*dphidv/(usiz*uvsin)
     &        + uwz*dphidw/usiz
           call atomic_add(de(3,ialoc),real( du,md_p))
           call atomic_add(de(3,ibloc),real(-du,md_p))
           frcz(3,i)   = frcz(3,i)   + du
c
c     force distribution for the Z-then-X local coordinate method
c
        else if (axetyp .eq. Ax_Z_Then_X) then

           du =  uvx*dphidv/(usiz*uvsin) + uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           du =  uvy*dphidv/(usiz*uvsin) + uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           du =  uvz*dphidv/(usiz*uvsin) + uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin)
           dt = - du - dv
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
           frcz(3,i)   = frcz(3,i) + du
           frcx(3,i)   = frcx(3,i) + dv

c
c     force distribution for the Bisector local coordinate method
c
        else if (axetyp .eq. Ax_Bisector) then

           du =  uvx*dphidv/(usiz*uvsin) + 0.5_ti_p*uwx*dphidw/usiz
           dv = -uvx*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwx*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           du =  uvy*dphidv/(usiz*uvsin) + 0.5_ti_p*uwy*dphidw/usiz
           dv = -uvy*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwy*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           du =  uvz*dphidv/(usiz*uvsin) + 0.5_ti_p*uwz*dphidw/usiz
           dv = -uvz*dphidu/(vsiz*uvsin) + 0.5_ti_p*vwz*dphidw/vsiz
           dt = - du - dv
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
           frcz(3,i)   = frcz(3,i)   + du
           frcx(3,i)   = frcx(3,i)   + dv

c
c     force distribution for the Z-Bisect local coordinate method
c
        else if (axetyp .eq. Ax_Z_Bisect) then

           du = urx*dphidr/(usiz*ursin) 
     &        + usx*dphids/usiz
           dv = (vssin*sx-vscos*t1x)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sx-wscos*t2x)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(1,ialoc),real(du,md_p))
           call atomic_add(de(1,icloc),real(dv,md_p))
           call atomic_add(de(1,idloc),real(dw,md_p))
           call atomic_add(de(1,ibloc),real(dt,md_p))
           frcz(1,i)   = frcz(1,i)   + du
           frcx(1,i)   = frcx(1,i)   + dv
           frcy(1,i)   = frcy(1,i)   + dw
           du = ury*dphidr/(usiz*ursin) + usy*dphids/usiz
           dv = (vssin*sy-vscos*t1y)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sy-wscos*t2y)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(2,ialoc),real(du,md_p))
           call atomic_add(de(2,icloc),real(dv,md_p))
           call atomic_add(de(2,idloc),real(dw,md_p))
           call atomic_add(de(2,ibloc),real(dt,md_p))
           frcz(2,i)   = frcz(2,i)   + du
           frcx(2,i)   = frcx(2,i)   + dv
           frcy(2,i)   = frcy(2,i)   + dw
           du = urz*dphidr/(usiz*ursin) + usz*dphids/usiz
           dv = (vssin*sz-vscos*t1z)*dphidu
     &        / (vsiz*(ut1sin+ut2sin))
           dw = (wssin*sz-wscos*t2z)*dphidu
     &        / (wsiz*(ut1sin+ut2sin))
           dt = - du - dv - dw
           call atomic_add(de(3,ialoc),real(du,md_p))
           call atomic_add(de(3,icloc),real(dv,md_p))
           call atomic_add(de(3,idloc),real(dw,md_p))
           call atomic_add(de(3,ibloc),real(dt,md_p))
           frcz(3,i)   = frcz(3,i)   + du
           frcx(3,i)   = frcx(3,i)   + dv
           frcy(3,i)   = frcy(3,i)   + dw

c
c     force distribution for the 3-Fold local coordinate method
c
        else if (axetyp .eq. Ax_3_Fold) then
           px     = ux + vx + wx
           py     = uy + vy + wy
           pz     = uz + vz + wz
           psiz   = sqrt(px*px + py*py + pz*pz) 
                  
           px     = px / psiz
           py     = py / psiz
           pz     = pz / psiz
                  
           wpcos  = wx*px + wy*py + wz*pz
           upcos  = ux*px + uy*py + uz*pz
           vpcos  = vx*px + vy*py + vz*pz
           wpsin  = sqrt(abs(1.0_ti_p - wpcos*wpcos))
           upsin  = sqrt(abs(1.0_ti_p - upcos*upcos))
           vpsin  = sqrt(abs(1.0_ti_p - vpcos*vpcos))
                  
           rx     = ux + vx
           ry     = uy + vy
           rz     = uz + vz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rwcos  = rx*wx + ry*wy + rz*wz
           rwsin  = sqrt(abs(1.0_ti_p - rwcos*rwcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*wz - rz*wy
           del(2) = rz*wx - rx*wz
           del(3) = rx*wy - ry*wx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*wz - del(3)*wy
           eps(2) = del(3)*wx - del(1)*wz
           eps(3) = del(1)*wy - del(2)*wx

!$acc loop seq
           do j = 1, 3
              dw = del(j)*dphidr/(wsiz*rwsin)
     &           + eps(j)*dphiddel*wpcos/(wsiz*psiz) 
              call atomic_add(de(j,idloc),real( dw,md_p))
              call atomic_add(de(j,ibloc),real(-dw,md_p))
              frcy(j,i)   = frcy(j,i)   + dw
           end do         
           rx     = vx + wx
           ry     = vy + wy
           rz     = vz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)

           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rucos  = rx*ux + ry*uy + rz*uz
           rusin  = sqrt(abs(1.0_ti_p - rucos*rucos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*uz - rz*uy
           del(2) = rz*ux - rx*uz
           del(3) = rx*uy - ry*ux
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*uz - del(3)*uy
           eps(2) = del(3)*ux - del(1)*uz
           eps(3) = del(1)*uy - del(2)*ux
!$acc loop seq
           do j = 1, 3
              du = del(j)*dphidr/(usiz*rusin)
     &           + eps(j)*dphiddel*upcos/(usiz*psiz)
              call atomic_add(de(j,ialoc),real( du,md_p))
              call atomic_add(de(j,ibloc),real(-du,md_p))
              frcz(j,i)   = frcz(j,i)   + du
           end do
           rx     = ux + wx
           ry     = uy + wy
           rz     = uz + wz
           rsiz   = sqrt(rx*rx + ry*ry + rz*rz)
                  
           rx     = rx / rsiz
           ry     = ry / rsiz
           rz     = rz / rsiz

           rvcos  = rx*vx + ry*vy + rz*vz 
           rvsin  = sqrt(abs(1.0_ti_p - rvcos*rvcos))
           dphidr = -trq(1)*rx - trq(2)*ry - trq(3)*rz
           del(1) = ry*vz - rz*vy
           del(2) = rz*vx - rx*vz
           del(3) = rx*vy - ry*vx
           delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))

           del(1) = del(1) / delsiz
           del(2) = del(2) / delsiz
           del(3) = del(3) / delsiz

           dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
           eps(1) = del(2)*vz - del(3)*vy
           eps(2) = del(3)*vx - del(1)*vz
           eps(3) = del(1)*vy - del(2)*vx
!$acc loop seq
           do j = 1, 3
              dv = del(j)*dphidr/(vsiz*rvsin)
     &           + eps(j)*dphiddel*vpcos/(vsiz*psiz)
              call atomic_add(de(j,icloc),real( dv,md_p))
              call atomic_add(de(j,ibloc),real(-dv,md_p))
              frcx(j,i)   = frcx(j,i)   + dv
           end do
        end if
      end do

      call timer_exit( timer_torque )
      end
