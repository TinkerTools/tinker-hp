c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey1" calculates the Urey-Bradley interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
      subroutine eurey1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'bound.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'urey.i'
      include 'urypot.i'
      include 'usage.i'
      include 'virial.i'
      include 'openmp.i'
      integer i,ia,ic,iurey
      integer ialoc,icloc
      real*8 e,de,ideal,force
      real*8 dt,dt2,deddt,fgrp
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 xac,yac,zac,rac
      logical proceed
c
c
c     zero out the Urey-Bradley energy and first derivatives
c
      eub = 0.0d0
c
c     calculate the Urey-Bradley 1-3 energy and first derivatives
c
      do iurey = 1, nureyloc
         i = ureyglob(iurey)
         ia = iury(1,i)
         ic = iury(3,i)
         ialoc = loc(ia)
         icloc = loc(ic)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ic,0,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ic))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xac = x(ia) - x(ic)
            yac = y(ia) - y(ic)
            zac = z(ia) - z(ic)
            if (use_polymer)  call image (xac,yac,zac)
            rac = sqrt(xac*xac + yac*yac + zac*zac)
            dt = rac - ideal
            dt2 = dt * dt
            e = ureyunit * force * dt2 * (1.0d0+cury*dt+qury*dt2)
            deddt = 2.0d0 * ureyunit * force * dt
     &                 * (1.0d0+1.5d0*cury*dt+2.0d0*qury*dt2)
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               deddt = deddt * fgrp
            end if
c
c     compute chain rule terms needed for derivatives
c
            de = deddt / rac
            dedx = de * xac
            dedy = de * yac
            dedz = de * zac
c
c     increment the total Urey-Bradley energy and first derivatives
c
            eub = eub + e
            deub(1,ialoc) = deub(1,ialoc) + dedx
            deub(2,ialoc) = deub(2,ialoc) + dedy
            deub(3,ialoc) = deub(3,ialoc) + dedz
c
            deub(1,icloc) = deub(1,icloc) - dedx
            deub(2,icloc) = deub(2,icloc) - dedy
            deub(3,icloc) = deub(3,icloc) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xac * dedx
            vyx = yac * dedx
            vzx = zac * dedx
            vyy = yac * dedy
            vzy = zac * dedy
            vzz = zac * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
      return
      end
