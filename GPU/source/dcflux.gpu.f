c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dcflux  --  charge flux gradient chain rule  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dcflux" takes as input the electrostatic potential at each
c     atomic site and calculates gradient chain rule terms due to
c     charge flux coupling with bond stretching and angle bending
c
c     literature reference:
c
c     C. Liu, J.-P. Piquemal and P. Ren, "Implementation of Geometry-
c     Dependent Charge Flux into the Polarizable AMOEBA+ Potential",
c     Journal of Physical Chemistry Letters, 11, 419-426 (2020)
c
c
#include "tinker_macro.h"
      submodule(cflux) subdcflux

#include "atomicOp.h.f"

      contains
#include "image.f.inc"
#include "atomicOp.inc.f"
#include "convert.f.inc"

      module subroutine dcflux1 (pot,dcf)
      use atmlst
      use angle
      use atoms
      use bond
      use bound
      use cflux
      use domdec
      use inform
      use sizes
      use utilgpu
      use virial
      implicit none
      integer   i,ia,ib,ic,ialoc,ibloc,icloc,iangle,ibond
      real(t_p) xa,ya,za,xb,yb,zb,xc,yc,zc,xab,yab,zab,xcb,ycb,zcb
      real(t_p) rab,rab2,rab3,rcb,rcb2,rcb3
      real(t_p) dpot,dpota,dpotc
      real(t_p) ddqdx,ddqdy,ddqdz
      real(t_p) fx,fy,fz
      real(t_p) fxa1,fya1,fza1
      real(t_p) fxb1,fyb1,fzb1
      real(t_p) fxc1,fyc1,fzc1
      real(t_p) fxa2,fya2,fza2
      real(t_p) fxb2,fyb2,fzb2
      real(t_p) fxc2,fyc2,fzc2
      real(t_p) pb,pb1,pb2
      real(t_p) pa1,pa2
      real(t_p) eps,dot
      real(t_p) term,fterm
      real(t_p) termxa,termxc,termya,termyc,termza,termzc
      real(t_p) pot(*)
      real(r_p) dcf(*)
c
c     set tolerance for minimum distance and angle values
c
      parameter(eps = 0.0001)
      if (deb_Path) print*, "  dcflux1"
c
c     calculate the charge flux forces due to bond stretches
c
!$acc parallel loop async(def_queue) default(present)
      do ibond = 1, nbondloc
         i = bndglob(ibond)
         ia = ibnd(1,i)
         ialoc = loc(ia)
         ib = ibnd(2,i)
         ibloc = loc(ib)
         pb = bflx(i)
         xa = x(ia)
         ya = y(ia)
         za = z(ia)
         xb = x(ib)
         yb = y(ib)
         zb = z(ib)
         xab = xa - xb
         yab = ya - yb
         zab = za - zb
         if (use_polymer)  call image_inl (xab,yab,zab)
         rab2 = max(xab*xab+yab*yab+zab*zab,eps)
         pb = pb / sqrt(rab2)
         dpot = pot(ibloc) - pot(ialoc)
         ddqdx = pb * (xa-xb)
         ddqdy = pb * (ya-yb)
         ddqdz = pb * (za-zb)
         fx = dpot * ddqdx
         fy = dpot * ddqdy
         fz = dpot * ddqdz
         ialoc = (ialoc-1)*3
         ibloc = (ibloc-1)*3
         call atomic_add( dcf(ialoc+1),fx )
         call atomic_add( dcf(ialoc+2),fy )
         call atomic_add( dcf(ialoc+3),fz )
         call atomic_sub( dcf(ibloc+1),fx )
         call atomic_sub( dcf(ibloc+2),fy )
         call atomic_sub( dcf(ibloc+3),fz )
      end do
c
c     calculate the charge flux forces due to angle bends
c
!$acc parallel loop async(def_queue) default(present)
      do iangle = 1, nangleloc
         i     = angleglob(iangle)
         ia    = iang(1,i)
         ialoc = loc(ia)
         ib    = iang(2,i)
         ibloc = loc(ib)
         ic    = iang(3,i)
         icloc = loc(ic)
         pa1 = aflx(1,i)
         pa2 = aflx(2,i)
         pb1 = abflx(1,i)
         pb2 = abflx(2,i)
         xa  = x(ia)
         ya  = y(ia)
         za  = z(ia)
         xb  = x(ib)
         yb  = y(ib)
         zb  = z(ib)
         xc  = x(ic)
         yc  = y(ic)
         zc  = z(ic)
         xab = xa - xb
         yab = ya - yb
         zab = za - zb
         xcb = xc - xb
         ycb = yc - yb
         zcb = zc - zb
         if (use_polymer) then
            call image_inl (xab,yab,zab)
            call image_inl (xcb,ycb,zcb)
         end if
         rab2  = max(xab*xab+yab*yab+zab*zab,eps)
         rcb2  = max(xcb*xcb+ycb*ycb+zcb*zcb,eps)
         rab   = sqrt(rab2)
         rab3  = rab2 * rab
         rcb   = sqrt(rcb2)
         rcb3  = rcb2 * rcb
c
c     get terms corresponding to asymmetric bond stretches
c
         dpota = pot(ialoc) - pot(ibloc)
         dpotc = pot(icloc) - pot(ibloc)
         pb1   = dpota * pb1
         pb2   = dpotc * pb2
         fxa1  = pb2 * xab/rab
         fya1  = pb2 * yab/rab
         fza1  = pb2 * zab/rab
         fxc1  = pb1 * xcb/rcb
         fyc1  = pb1 * ycb/rcb
         fzc1  = pb1 * zcb/rcb
         fxb1  = -fxa1 - fxc1
         fyb1  = -fya1 - fyc1
         fzb1  = -fza1 - fzc1
c
c     get terms corresponding to bond angle bending
c
         dot    = xab*xcb + yab*ycb + zab*zcb
         term   = -rab*rcb / max(sqrt(rab2*rcb2-dot*dot),eps)
         fterm  = term*(dpota*pa1+dpotc*pa2)
         termxa = xcb/(rab*rcb) - xab*dot/(rab3*rcb)
         termya = ycb/(rab*rcb) - yab*dot/(rab3*rcb)
         termza = zcb/(rab*rcb) - zab*dot/(rab3*rcb)
         termxc = xab/(rab*rcb) - xcb*dot/(rab*rcb3)
         termyc = yab/(rab*rcb) - ycb*dot/(rab*rcb3)
         termzc = zab/(rab*rcb) - zcb*dot/(rab*rcb3)
         fxa2   = fterm * termxa
         fya2   = fterm * termya
         fza2   = fterm * termza
         fxc2   = fterm * termxc
         fyc2   = fterm * termyc
         fzc2   = fterm * termzc
         fxb2   = -fxa2 - fxc2
         fyb2   = -fya2 - fyc2
         fzb2   = -fza2 - fzc2
         ialoc  = (ialoc-1)*3
         ibloc  = (ibloc-1)*3
         icloc  = (icloc-1)*3
         call atomic_add( dcf(ialoc+1),fxa1+fxa2 )
         call atomic_add( dcf(ialoc+2),fya1+fya2 )
         call atomic_add( dcf(ialoc+3),fza1+fza2 )
         call atomic_add( dcf(ibloc+1),fxb1+fxb2 )
         call atomic_add( dcf(ibloc+2),fyb1+fyb2 )
         call atomic_add( dcf(ibloc+3),fzb1+fzb2 )
         call atomic_add( dcf(icloc+1),fxc1+fxc2 )
         call atomic_add( dcf(icloc+2),fyc1+fyc2 )
         call atomic_add( dcf(icloc+3),fzc1+fzc2 )
      end do
      end

      ! Same as dcflux2 except for _dcf_ type which is mdyn_rtyp
      module subroutine dcflux2 (pot,dcf)
      use atmlst
      use angle
      use atoms
      use bond
      use bound
      use cflux
      use domdec
      use inform
      use sizes
      use utilgpu
      use tinheader
      use virial
      implicit none
      real(t_p),intent(in   ):: pot(*)
      mdyn_rtyp,intent(inout):: dcf(*)
      integer   i,ia,ib,ic,ialoc,ibloc,icloc,iangle,ibond
      real(t_p) xa,ya,za,xb,yb,zb,xc,yc,zc,xab,yab,zab,xcb,ycb,zcb
      real(t_p) rab,rab2,rab3,rcb,rcb2,rcb3
      real(t_p) dpot,dpota,dpotc,ddqdx,ddqdy,ddqdz
      real(t_p) fx,fy,fz
      real(t_p) fxa1,fya1,fza1,fxb1,fyb1,fzb1
      real(t_p) fxc1,fyc1,fzc1,fxa2,fya2,fza2
      real(t_p) fxb2,fyb2,fzb2,fxc2,fyc2,fzc2
      real(t_p) pb,pb1,pb2,pa1,pa2
      real(t_p) eps,dot
      real(t_p) term,fterm
      real(t_p) termxa,termxc,termya,termyc,termza,termzc
c
c     set tolerance for minimum distance and angle values
c
      parameter(eps = 0.0001)
      if (deb_Path) print*, "  dcflux2"
c
c     calculate the charge flux forces due to bond stretches
c
!$acc parallel loop async(def_queue) default(present)
      do ibond = 1, nbondloc
         i     = bndglob(ibond)
         ia    = ibnd(1,i)
         ib    = ibnd(2,i)
         pb    = bflx(i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         xa    = x(ia)
         ya    = y(ia)
         za    = z(ia)
         xb    = x(ib)
         yb    = y(ib)
         zb    = z(ib)
         xab   = xa - xb
         yab   = ya - yb
         zab   = za - zb
         if (use_polymer)  call image_inl (xab,yab,zab)
         rab2  = max(xab*xab+yab*yab+zab*zab,eps)
         pb    = pb / sqrt(rab2)
         dpot  = pot(ibloc) - pot(ialoc)
         ddqdx = pb * (xa-xb)
         ddqdy = pb * (ya-yb)
         ddqdz = pb * (za-zb)
         fx    = dpot * ddqdx
         fy    = dpot * ddqdy
         fz    = dpot * ddqdz
         ialoc = (ialoc-1)*3
         ibloc = (ibloc-1)*3
         call atomic_add( dcf(ialoc+1),fx )
         call atomic_add( dcf(ialoc+2),fy )
         call atomic_add( dcf(ialoc+3),fz )
         call atomic_sub( dcf(ibloc+1),fx )
         call atomic_sub( dcf(ibloc+2),fy )
         call atomic_sub( dcf(ibloc+3),fz )
      end do
c
c     calculate the charge flux forces due to angle bends
c
!$acc parallel loop async(def_queue) default(present)
      do iangle = 1, nangleloc
         i     = angleglob(iangle)
         ia    = iang(1,i)
         ib    = iang(2,i)
         ic    = iang(3,i)
         pa1   = aflx(1,i)
         pa2   = aflx(2,i)
         pb1   = abflx(1,i)
         pb2   = abflx(2,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         xa  = x(ia)
         ya  = y(ia)
         za  = z(ia)
         xb  = x(ib)
         yb  = y(ib)
         zb  = z(ib)
         xc  = x(ic)
         yc  = y(ic)
         zc  = z(ic)
         xab = xa - xb
         yab = ya - yb
         zab = za - zb
         xcb = xc - xb
         ycb = yc - yb
         zcb = zc - zb
         if (use_polymer) then
            call image_inl (xab,yab,zab)
            call image_inl (xcb,ycb,zcb)
         end if
         rab2  = max(xab*xab+yab*yab+zab*zab,eps)
         rcb2  = max(xcb*xcb+ycb*ycb+zcb*zcb,eps)
         rab   = sqrt(rab2)
         rab3  = rab2 * rab
         rcb   = sqrt(rcb2)
         rcb3  = rcb2 * rcb
c
c     get terms corresponding to asymmetric bond stretches
c
         dpota = pot(ialoc) - pot(ibloc)
         dpotc = pot(icloc) - pot(ibloc)
         pb1   = dpota * pb1
         pb2   = dpotc * pb2
         fxa1  = pb2 * xab/rab
         fya1  = pb2 * yab/rab
         fza1  = pb2 * zab/rab
         fxc1  = pb1 * xcb/rcb
         fyc1  = pb1 * ycb/rcb
         fzc1  = pb1 * zcb/rcb
         fxb1  = -fxa1 - fxc1
         fyb1  = -fya1 - fyc1
         fzb1  = -fza1 - fzc1
c
c     get terms corresponding to bond angle bending
c
         dot    = xab*xcb + yab*ycb + zab*zcb
         term   = -rab*rcb / max(sqrt(rab2*rcb2-dot*dot),eps)
         fterm  = term*(dpota*pa1+dpotc*pa2)
         termxa = xcb/(rab*rcb) - xab*dot/(rab3*rcb)
         termya = ycb/(rab*rcb) - yab*dot/(rab3*rcb)
         termza = zcb/(rab*rcb) - zab*dot/(rab3*rcb)
         termxc = xab/(rab*rcb) - xcb*dot/(rab*rcb3)
         termyc = yab/(rab*rcb) - ycb*dot/(rab*rcb3)
         termzc = zab/(rab*rcb) - zcb*dot/(rab*rcb3)
         fxa2   = fterm * termxa
         fya2   = fterm * termya
         fza2   = fterm * termza
         fxc2   = fterm * termxc
         fyc2   = fterm * termyc
         fzc2   = fterm * termzc
         fxb2   = -fxa2 - fxc2
         fyb2   = -fya2 - fyc2
         fzb2   = -fza2 - fzc2
         ialoc  = (ialoc-1)*3
         ibloc  = (ibloc-1)*3
         icloc  = (icloc-1)*3
         call atomic_add( dcf(ialoc+1),fxa1+fxa2 )
         call atomic_add( dcf(ialoc+2),fya1+fya2 )
         call atomic_add( dcf(ialoc+3),fza1+fza2 )
         call atomic_add( dcf(ibloc+1),fxb1+fxb2 )
         call atomic_add( dcf(ibloc+2),fyb1+fyb2 )
         call atomic_add( dcf(ibloc+3),fzb1+fzb2 )
         call atomic_add( dcf(icloc+1),fxc1+fxc2 )
         call atomic_add( dcf(icloc+2),fyc1+fyc2 )
         call atomic_add( dcf(icloc+3),fzc1+fzc2 )
      end do
      end
c
c     "adflux" increments global gradient buffer _de_
c
      module subroutine adflux1(dcf,de)
      use atoms
      use atmlst
      use domdec
      use mpole
      use utilgpu
      use virial
      implicit none
      real(r_p) dcf(*),de(*)
      integer   ii,iipole,iglob,i
      real(t_p) xi,yi,zi
      real(r_p) frcx,frcy,frcz

!$acc parallel loop async(def_queue) default(present)
!$acc&         present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     reduction(+:g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do ii = 1, npolebloc
         iglob   = ipole(poleglob(ii))
         i       = (loc(iglob)-1)*3
         frcx    = dcf(i+1)
         frcy    = dcf(i+2)
         frcz    = dcf(i+3)
         de(i+1) = de (i+1) + frcx
         de(i+2) = de (i+2) + frcy
         de(i+3) = de (i+3) + frcz

         if (use_virial) then
            xi   = x(iglob)
            yi   = y(iglob)
            zi   = z(iglob)
            g_vxx= g_vxx + xi*(frcx)
            g_vxy= g_vxy + yi*(frcx)
            g_vxz= g_vxz + zi*(frcx)
            g_vyy= g_vyy + yi*(frcy)
            g_vyz= g_vyz + zi*(frcy)
            g_vzz= g_vzz + zi*(frcz)
         end if

         dcf(i+1) = 0
         dcf(i+2) = 0
         dcf(i+3) = 0
      end do
      end subroutine
      ! Same as adflux1 except for _de_ type which is mdyn_rtyp
      module subroutine adflux2(dcf,de)
      use atoms
      use atmlst
      use domdec
      use inform
      use mpole
      use utilgpu
      use virial
      implicit none
      mdyn_rtyp dcf(*)
      mdyn_rtyp  de(*)
      integer   ii,iipole,iglob,i
      real(t_p) xi,yi,zi
      mdyn_rtyp frcx,frcy,frcz

      if (use_virial) then

!$acc parallel loop async(def_queue) default(present)
!$acc&         present(g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     reduction(+:g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do ii = 1, npolebloc
         iglob   = ipole(poleglob(ii))
         i       = (loc(iglob)-1)*3
         frcx    = dcf(i+1)
         frcy    = dcf(i+2)
         frcz    = dcf(i+3)
         de(i+1) = de (i+1) + frcx
         de(i+2) = de (i+2) + frcy
         de(i+3) = de (i+3) + frcz

         xi      = x(iglob)
         yi      = y(iglob)
         zi      = z(iglob)
         g_vxx   = g_vxx + xi*mdr2md(frcx)
         g_vxy   = g_vxy + yi*mdr2md(frcx)
         g_vxz   = g_vxz + zi*mdr2md(frcx)
         g_vyy   = g_vyy + yi*mdr2md(frcy)
         g_vyz   = g_vyz + zi*mdr2md(frcy)
         g_vzz   = g_vzz + zi*mdr2md(frcz)

         dcf(i+1)= 0
         dcf(i+2)= 0
         dcf(i+3)= 0
      end do

      else

!$acc parallel loop async(def_queue) default(present)
      do ii = 1, npolebloc
         iglob   = ipole(poleglob(ii))
         i       = (loc(iglob)-1)*3
         frcx    = dcf(i+1)
         frcy    = dcf(i+2)
         frcz    = dcf(i+3)
         de(i+1) = de (i+1) + frcx
         de(i+2) = de (i+2) + frcy
         de(i+3) = de (i+3) + frcz

         dcf(i+1)= 0
         dcf(i+2)= 0
         dcf(i+3)= 0
      end do

      end if
      end subroutine

      end submodule
