c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1  --  buffered 14-7 energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      subroutine ehal1
      use deriv,only:dev
      use energi
      use potent
      use virial
      use vdwpot
      use mpi
      implicit none
      real(t_p) elrc,vlrc
c
c
c     choose the method for summing over pairwise interactions
c
      call ehal1c
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr1 (elrc,vlrc)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1c  --  buffered 14-7 vdw derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1c" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise neighbor list
c
c
      subroutine ehal1c
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use group
      use inter
      use iounit
      use molcul
      use mutant
      use neigh
      use potent
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis,iivdw,nnvlst
      integer ii,iv,it,ivloc
      integer kk,kv,kt,kvloc
      integer, allocatable :: iv14(:)
      real(t_p) e,de,eps,rdn,fgrp
      real(t_p) rv,rv7
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) redi,rediv
      real(t_p) redk,redkv
      real(t_p) dedx,dedy,dedz
      real(t_p) rho,rho6,rho7
      real(t_p) tau,tau7,scal
      real(t_p) s1,s2,t1,t2
      real(t_p) dt1drho,dt2drho
      real(t_p) dt1dlambda, dt2dlambda
      real(t_p) ds1dlambda, ds2dlambda
      real(t_p) dtau,gtau
      real(t_p) taper,dtaper
      real(t_p) rik,rik2,rik3
      real(t_p) rik4,rik5
      real(t_p) rik6,rik7
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      real(t_p) s,ds,vdwshortcut2,facts,factds
      real(t_p), allocatable :: xred(:)
      real(t_p), allocatable :: yred(:)
      real(t_p), allocatable :: zred(:)
      real(t_p), allocatable :: vscale(:)
      logical proceed,usei
      logical testcut,shortrange,longrange,fullrange
      logical muti,mutk,mutik
      character*11 mode
      character*80 :: RoutineName
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate VDW')

c     choose the method for summing over pairwise interactions
      shortrange = use_vdwshort
      longrange  = use_vdwlong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'ehalshort1c'
         mode        = 'SHORTVDW'
      else if (longrange) then
         RoutineName = 'ehallong1c'
         mode        = 'VDW'
      else
         RoutineName = 'ehal1c'
         mode        = 'VDW'
      endif

c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0
      dev = 0
      delambdav = 0
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1.0
      iv14 = 0
c
c     set the coefficients for the switching function
c
      call switch (mode)
      vdwshortcut2 = (vdwshortcut-shortheal)**2
c
c     apply any reduction factor to the atomic coordinates
c
      do ii = 1, nvdwbloc
         iivdw = vdwglob(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         rdn = kred(iglob)
         xr = x(iglob) - x(iv)
         yr = y(iglob) - y(iv)
         zr = z(iglob) - z(iv)
         if (use_polymer) call image(xr,yr,zr)
         xred(i) = rdn*xr + x(iv)
         yred(i) = rdn*yr + y(iv)
         zred(i) = rdn*zr + z(iv)
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         if ((i.eq.0).or.(i.gt.nbloc)) then
           write(iout,1000)
           cycle
         end if
         iv = ired(iglob)
         ivloc = loc(iv)
         redi = kred(iglob)
         rediv = 1.0 - redi
         it = jvdw(iglob)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(iglob) .or. use(iv))
         muti = mut(iglob)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = v2scale
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = v3scale
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = v4scale
            iv14(i14(j,iglob)) = iglob
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         if (shortrange) then
           nnvlst = nshortvlst(ii)
         else
           nnvlst = nvlst(ii)
         end if
         do kk = 1, nnvlst
            if (shortrange) then
              kglob = shortvlst(kk,ii)
            else
              kglob = vlst(kk,ii)
            end if
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            kbis = loc(kglob)
            kv = ired(kglob)
            kvloc = loc(kv)
            if ((kvloc.eq.0).or.(kbis.gt.nbloc).or.(kvloc.gt.nbloc))
     $         then
               write(iout,1000)
               cycle
            end if
            redk = kred(kglob)
            redkv = 1.0 - redk
            mutk = mut(kglob)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)

               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               testcut = merge(rik2 .le. off2.and.rik2.ge.vdwshortcut2,
     &                         rik2 .le. off2,
     &                         longrange
     &                        )
               if (testcut) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)

c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = .false.
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        mutik = .true.
                     else if (.not.muti .or. .not.mutk) then
                        mutik = .true.
                     end if
                  end if
c
c     get the energy and gradient, via soft core if necessary
c
                  if (mutik) then
                     rho = rik / rv
                     rho6 = rho**6
                     rho7 = rho6 * rho
c                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0-vlambda)**2
                     s1 = 1.0 / (scal+(rho+dhal)**7)
                     s2 = 1.0 / (scal+rho7+ghal)
                     t1 = (1.0+dhal)**7 * s1
                     t2 = (1.0+ghal) * s2
                     dt1drho = -7.0*(rho+dhal)**6 * t1 * s1
                     dt2drho = -7.0*rho6 * t2 * s2
                     e = eps * (vlambda**scexp)* t1 * (t2-2.0)
                     if (use_lambdadyn) then
                       ds1dlambda=2*scalpha*(1.0-vlambda)*s1*s1
                       ds2dlambda=2*scalpha*(1.0-vlambda)*s2*s2
                       dt1dlambda=(1.0+dhal)**7 *ds1dlambda
                       dt2dlambda=(1.0+ghal)*ds2dlambda
                       delambdav = delambdav +
     $                    eps*scexp*vlambda**(scexp-1)*t1*(t2-2.0)+
     $                    eps*(vlambda**scexp)*(dt1dlambda*(t2-2.0)+
     $                    t1*dt2dlambda)
                     end if
                     de = eps * (vlambda**scexp) *
     $                    (dt1drho*(t2-2.0)+t1*dt2drho) / rv
                  else
                     rv7 = rv**7
                     rik6 = rik2**3
                     rik7 = rik6 * rik
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0) / (rik + dhal*rv)
                     tau7 = tau**7
                     dtau = tau / (dhal+1.0)
                     gtau = eps*tau7*rik6*(ghal+1.0)*(rv7/rho)**2
                     e = eps*tau7*rv7*((ghal+1.0)*rv7/rho-2.0)
                     de = -7.0 * (dtau*e+gtau)
                  end if
c
c     use energy switching if near the cutoff distance
c     at short or long range
c
                  if(longrange.or.fullrange) then
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                        + c2*rik2 + c1*rik  + c0
                        dtaper =   5.0*c5*rik4 + 4.0*c4*rik3
     &                           + 3.0*c3*rik2 + 2.0*c2*rik + c1
                        de = e * dtaper + de * taper
                        e  = e * taper
                     end if
                  endif
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     use energy switching if close the cutoff distance (at short range)
c
                  if(shortrange .or. longrange)
     &               call switch_respa(rik,vdwshortcut,shortheal,s,ds)

                  if(shortrange) then
                     facts =          s
                     factds =        ds
                  else if(longrange) then
                     facts  = 1.0 - s
                     factds =      - ds
                  else
                     facts  = 1.0
                     factds = 0.0
                  endif

                  de = de * facts + e * factds
                  e  = e  * facts

                  ev = ev + e
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  if (iglob .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
c
                     dev(1,ivloc) = dev(1,ivloc) + dedx*rediv
                     dev(2,ivloc) = dev(2,ivloc) + dedy*rediv
                     dev(3,ivloc) = dev(3,ivloc) + dedz*rediv
                  end if
                  if (kglob .eq. kv) then
                     dev(1,kbis) = dev(1,kbis) - dedx
                     dev(2,kbis) = dev(2,kbis) - dedy
                     dev(3,kbis) = dev(3,kbis) - dedz
                  else
                     dev(1,kbis) = dev(1,kbis) - dedx*redk
                     dev(2,kbis) = dev(2,kbis) - dedy*redk
                     dev(3,kbis) = dev(3,kbis) - dedz*redk
c
                     dev(1,kvloc) = dev(1,kvloc) - dedx*redkv
                     dev(2,kvloc) = dev(2,kvloc) - dedy*redkv
                     dev(3,kvloc) = dev(3,kvloc) - dedz*redkv
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz

                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = 1.0
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = 1.0
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = 1.0
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = 1.0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)

      return
      end
