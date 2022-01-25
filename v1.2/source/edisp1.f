c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edisp1  --  damped dispersion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp1" calculates the damped dispersion energy and first
c     derivatives with respect to Cartesian coordinates
c
c     literature reference:
c
c     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
c     Grounded Damped Dispersion Model with Particle Mesh Ewald
c     Summation", Journal of Chemical Physics, 149, 084115 (2018)
c
c
      subroutine edisp1
      use dsppot
      use energi
      use ewald
      use potent
      use virial
      implicit none
      real*8 elrc,vlrc
      character*11 mode
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_dewald) then
        call edisp1d
      else
        call edisp1b
      end if
c
c     apply long range dispersion correction if desired
c
      if (use_dcorr .and. .not.use_dewald) then
         mode = 'DISP'
         call evcorr1 (mode,elrc,vlrc)
         edsp = edsp + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edisp1b  --  neighbor list dispersion derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edisp1b" calculates the damped dispersion energy and
c     derivatives with respect to Cartesian coordinates using
c     a pairwise neighbor list
c
c
      subroutine edisp1b
      use atoms
      use atmlst
      use bound
      use cell
      use couple
      use cutoff
      use deriv
      use disp
      use domdec
      use dsppot
      use energi
      use group
      use mutant
      use neigh
      use potent
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,iglob,kglob,kbis,nnvlst
      integer ii,iidisp,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ai2,ai3
      real*8 ak,ak2,ak3
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,ddamp
      real*8 taper,dtaper
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: dspscale(:)
      real*8 s,ds,dispshortcut2,facts,factds
      logical proceed,usei
      logical muti,mutk
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName
c
c     compute the short, long, or full real space part of the Ewald summation
      shortrange = use_dispshort
      longrange  = use_displong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'edispshort1b'
         mode        = 'SHORTDISP'
      else if (longrange) then
         RoutineName = 'edisplong1b'
         mode        = 'DISP'
      else
         RoutineName = 'edisp1b'
         mode        = 'DISP'
      endif
c
c     zero out the damped dispersion energy and derivatives
c
      edsp = 0.0d0
      dedsp = 0.0d0
      if (ndisp .eq. 0) return
c
c     perform dynamic allocation of some local arrays
c
      allocate (dspscale(n))
c
c     initialize connected atom exclusion coefficients
c
      dspscale = 1.0d0
c
c     set cutoff and switching coefficients
c
      call switch (mode)
      dispshortcut2 = (dispshortcut-shortheal)**2
c
c     find dispersion energy and derivatives via neighbor list
c
      do ii = 1, ndisplocnl
         iidisp = dispglobnl(ii)
         iglob = idisp(iidisp)
         i = loc(iglob)
         ci = csix(iidisp)
         ai = adisp(iidisp)
         usei = use(iglob)
         muti = mut(iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = dsp2scale
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = dsp3scale
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = dsp4scale
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = dsp5scale
         end do
c
c     decide whether to compute the current interaction
c
         nnvlst = merge(nshortvlst(ii),
     &                  nvlst     (ii),
     &                  shortrange
     &                 )
         do kkk = 1, nnvlst
            kk = merge(shortvlst(kkk,ii),
     &                    vlst     (kkk,ii),
     &                    shortrange
     &                   )
            kglob = idisp(kk)
            kbis = loc(kglob)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            ck = csix(kk)
            ak = adisp(kk)
            mutk = mut(kglob)
            proceed = (usei .or. use(kglob))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(kglob)
               yr = yi - y(kglob)
               zr = zi - z(kglob)
               if (use_bounds) call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               testcut = merge(r2 .le. off2.and.r2.ge.dispshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
               if (testcut) then
                  r = sqrt(r2)
                  r6 = r2**3
                  e = -ci * ck / r6
                  de = -6.0d0 * e / r
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ai3 = ai * ai2
                     ak2 = ak * ak
                     ak3 = ak * ak2
                     dk2 = dk * dk
                     dk3 = dk * dk2
                     ti = ak2 / (ak2-ai2)
                     tk = ai2 / (ai2-ak2)
                     ti2 = ti * ti
                     tk2 = tk * tk
                     damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2)*expk
     &                          - 2.0d0*ti2*tk*(1.0d0+di)*expi
     &                          - 2.0d0*tk2*ti*(1.0d0+dk)*expk
                     damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2
     &                                       +di3/6.0d0)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2
     &                                    +dk3/6.0d0)*expk
     &                          - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                          - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                     ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                          * (r*ai+4.0d0*tk-1.0d0)
     &                       + 0.25d0 * dk2 * tk2 * ak * expk
     &                            * (r*ak+4.0d0*ti-1.0d0)
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                          / 96.0d0
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     set use of lambda scaling for decoupling or annihilation
c
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        de = de * vlambda
                        e = e * vlambda
                     else if (.not.muti .or. .not.mutk) then
                        de = de * vlambda
                        e = e * vlambda
                     end if
                  end if
c
c     apply damping and scaling factors for this interaction
c
                  de = de*damp**2 + 2.0d0*e*damp*ddamp
                  e = e * damp**2
                  e = e * dspscale(kglob)
                  de = de * dspscale(kglob)
c
c     use energy switching if near the cutoff distance
c
                  if(longrange.or.fullrange) then
                    if (r2 .gt. cut2) then
                       r3 = r2 * r
                       r4 = r2 * r2
                       r5 = r2 * r3
                       taper = c5*r5 + c4*r4 + c3*r3
     &                            + c2*r2 + c1*r + c0
                       dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                             + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                       de = e*dtaper + de*taper
                       e = e * taper
                    end if
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if

                  if(shortrange .or. longrange)
     &               call switch_respa(r,chgshortcut,shortheal,s,ds)

                  if(shortrange) then
                     facts  =         s
                     factds =      + ds
                  else if(longrange) then
                     facts  = 1.0d0 - s
                     factds =      - ds
                  else
                     facts  = 1.0d0
                     factds = 0.0d0
                  endif

                  de = de * facts + e * factds
                  e  =  e * facts
c
c     increment the overall damped dispersion energy component
c
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  dedx = de * xr/r
                  dedy = de * yr/r
                  dedz = de * zr/r
                  dedsp(1,i) = dedsp(1,i) + dedx
                  dedsp(2,i) = dedsp(2,i) + dedy
                  dedsp(3,i) = dedsp(3,i) + dedz
                  dedsp(1,kbis) = dedsp(1,kbis) - dedx
                  dedsp(2,kbis) = dedsp(2,kbis) - dedy
                  dedsp(3,kbis) = dedsp(3,kbis) - dedz
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
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c                                                                        
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edisp1d  --  Ewald dispersion derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp1d" calculates the damped dispersion energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine edisp1d
      use atoms
      use atmlst
      use deriv
      use disp
      use domdec
      use energi
      use ewald
      use pme
      use potent
      implicit none
      integer i,iidisp
      real*8 term
c
c
c     zero out the damped dispersion energy and derivatives
c
      edsp = 0.0d0
      dedsp = 0d0
      if (ndisp .eq. 0)  return
      aewald = adewald
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        if (use_disprec) then
          call edrecip1
        end if
      end if
c
c     compute the real space portion of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_dispreal) then
          call edreal1d
        end if
c
c     compute the self-energy portion of the Ewald summation
c
        if (use_dispself) then
          do i = 1, ndisploc
             iidisp = dispglob(i)
             term = aewald**6 / 12.0d0
             edsp = edsp + term*csix(iidisp)*csix(iidisp)
          end do
        end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edreal1d  --  Ewald real disp derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to damped dispersion
c     interactions via a neighbor list
c
c
      subroutine edreal1d
      use atoms
      use atmlst
      use bound
      use boxes
      use couple
      use cell
      use cutoff
      use deriv
      use disp
      use domdec
      use dsppot
      use energi
      use ewald
      use group
      use mutant
      use neigh
      use potent
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,iglob,kglob,kbis,nnvlst
      integer ii,iidisp,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r6,r7
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,ddamp
      real*8 ralpha2,scale
      real*8 expterm,term
      real*8 expa,rterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 s,ds,dispshortcut2,facts,factds
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName

c     choose the method for summing over pairwise interactions
      shortrange = use_dispshortreal
      longrange  = use_displong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'edrealshort1c'
         mode        = 'SHORTDEWALD'
      else if (longrange) then
         RoutineName = 'edreallong1c'
         mode        = 'DEWALD'
      else
         RoutineName = 'ereal1c'
         mode        = 'DEWALD'
      endif
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dspscale(n))
c
c     initialize connected atom exclusion coefficients
c
      dspscale = 1.0d0
c
c     set cutoff and switching coefficients
c
      call switch (mode)
      dispshortcut2 = (dispshortcut-shortheal)**2
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisplocnl
         iidisp = dispglobnl(ii)
         iglob = idisp(iidisp)
         i = loc(iglob)
         ci = csix(iidisp)
         ai = adisp(iidisp)
         usei = use(iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         muti = mut(iglob)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = dsp2scale
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = dsp3scale
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = dsp4scale
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = dsp5scale
         end do
c
c     decide whether to compute the current interaction
c
         nnvlst = merge(nshortvlst(ii),
     &                  nvlst     (ii),
     &                  shortrange
     &                 )
         do kkk = 1, nnvlst
            kk = merge(shortvlst(kkk,ii),
     &                    vlst     (kkk,ii),
     &                    shortrange
     &                   )
            kglob = idisp(kk)
            kbis = loc(kglob)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            ck = csix(kk)
            ak = adisp(kk)
            mutk = mut(kglob)
            proceed = (usei .or. use(kglob))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(kglob)
               yr = yi - y(kglob)
               zr = zi - z(kglob)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               testcut = merge(r2 .le. off2.and.r2.ge.dispshortcut2,
     &                         r2 .le. off2,
     &                         longrange
     &                        )
               if (testcut) then
                  r6 = r2**3
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expterm = exp(-ralpha2)
                  expa = expterm * term
c
c     find the damping factor for the dispersion interaction
c
                  r = sqrt(r2)
                  r7 = r6 * r
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ak2 = ak * ak
                     dk2 = dk * dk
                     dk3 = dk * dk2
                     ti = ak2 / (ak2-ai2)
                     ti2 = ti * ti
                     tk = ai2 / (ai2-ak2)
                     tk2 = tk * tk
                     damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2)*expk
     &                          - 2.0d0*ti2*tk*(1.0d0+di)*expi
     &                          - 2.0d0*tk2*ti*(1.0d0+dk)*expk
                     damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2
     &                                       +di3/6.0d0)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2
     &                                    +dk3/6.0d0)*expk
     &                          - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                          - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                     ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                          * (r*ai+4.0d0*tk-1.0d0)
     &                       + 0.25d0 * dk2 * tk2 * ak * expk
     &                            * (r*ak+4.0d0*ti-1.0d0)
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                          / 96.0d0
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(kglob) * damp**2
                  if (use_group)  scale = scale * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        scale = scale * vlambda
                        ddamp = ddamp * vlambda
                     else if (.not.muti .or. .not.mutk) then
                        scale = scale * vlambda
                        ddamp = ddamp * vlambda
                     end if
                  end if

                  scale = scale - 1.0d0
                  e = -ci * ck * (expa+scale) / r6
                  rterm = -(ralpha2**3) * expterm / r
                  de = -6.0d0*e/r - ci*ck*rterm/r6
     &                    - 2.0d0*ci*ck*dspscale(kglob)*damp*ddamp/r6

                  if(shortrange .or. longrange)
     &               call switch_respa(r,chgshortcut,shortheal,s,ds)

                  if(shortrange) then
                     facts  =         s
                     factds =      + ds
                  else if(longrange) then
                     facts  = 1.0d0 - s
                     factds =      - ds
                  else
                     facts  = 1.0d0
                     factds = 0.0d0
                  endif

                  de = de * facts + e * factds
                  e  =  e * facts
c
c     increment the overall damped dispersion energy component
c
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  dedx = de * xr/r
                  dedy = de * yr/r
                  dedz = de * zr/r
                  dedsp(1,i) = dedsp(1,i) + dedx
                  dedsp(2,i) = dedsp(2,i) + dedy
                  dedsp(3,i) = dedsp(3,i) + dedz
                  dedsp(1,kbis) = dedsp(1,kbis) - dedx
                  dedsp(2,kbis) = dedsp(2,kbis) - dedy
                  dedsp(3,kbis) = dedsp(3,kbis) - dedz
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
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c                                                                        
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edrecip1  --  PME recip disp energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edrecip1" evaluates the reciprocal space portion of particle
c     mesh Ewald energy and gradient due to damped dispersion
c
c
      subroutine edrecip1
      use atmlst
      use boxes
      use bound
      use disp
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use math
      use mpi
      use pme
      use potent
      use virial
      implicit none
      integer i,j,k,ii,iidisp,iglob
      integer iproc
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff
      integer i0,iatm,igrd0
      integer it1,it2,it3
      integer j0,jgrd0
      integer k0,kgrd0
      integer istart,iend,jstart,jend,kstart,kend
      real*8 e,fi,denom
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 term,vterm
      real*8 expterm
      real*8 erfcterm
      real*8 hsq,struc2
      real*8 h,hhh,b,bfac
      real*8 term1,denom0
      real*8 fac1,fac2,fac3
      real*8 de1,de2,de3
      real*8 dn1,dn2,dn3
      real*8 dt1,dt2,dt3
      real*8 t1,t2,t3
      integer, allocatable :: req(:),reqbcast(:)
      integer status(MPI_STATUS_SIZE),tag,ierr,proc
      integer nprocloc,rankloc,commloc
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
c
      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc = comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = COMM_TINKER
      end if
c
c     set Ewald coefficient
c
      aewald = adewald
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     dynamic allocation of local arrays
c
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (req(nproc*nproc))
      allocate (reqbcast(nproc*nproc))
c
c
      do i = 1, ndisprecloc
        iidisp = disprecglob(i)
        iglob = idisp(iidisp)
        call bspline_fill_site(iglob,i)
      end do
c
      qgridin_2d = 0d0
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,commloc,req(tag),
     $   ierr)
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      do i = 1, ndisprecloc
        iidisp = disprecglob(i)
        iglob = idisp(iidisp)
        call grid_disp_site(iglob,i)
      enddo
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,
     $   commloc,req(tag),ierr)
      end do
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1) + 
     $    qgridmpi(:,:,:,:,i) 
      end do
c
c     perform the 3-D FFT forward transformation
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     use scalar sum to get the reciprocal space energy
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
         qgridout_2d(1,1,1,1) = 0.0d0
         qgridout_2d(2,1,1,1) = 0.0d0
      end if
      bfac = pi / aewald
      fac1 = 2.0d0*pi**(3.5d0)
      fac2 = aewald**3
      fac3 = -2.0d0*aewald*pi**2
      denom0 = (6.0d0*volbox)/(pi**1.5d0)
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      nff = nfft1 * nfft2
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
          do k1 = istart2(rankloc+1),iend2(rankloc+1)
           m1 = k1 - 1
           m2 = k2 - 1
           m3 = k3 - 1
           if (k1 .gt. nf1)  m1 = m1 - nfft1
           if (k2 .gt. nf2)  m2 = m2 - nfft2
           if (k3 .gt. nf3)  m3 = m3 - nfft3
           if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
           r1 = dble(m1)
           r2 = dble(m2)
           r3 = dble(m3)
           h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
           h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
           h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
           hsq = h1*h1 + h2*h2 + h3*h3
           h = sqrt(hsq)
           b = h*bfac
           hhh = h*hsq
           term = -b*b
           expterm = 0.0d0
           erfcterm = erfc(b)
           denom = denom0*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
           if (term .gt. -50.0d0) then
              expterm = exp(term)
              erfcterm = erfc(b)
              if (.not. use_bounds) then
                 expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
                 erfcterm = erfcterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
              else if (octahedron) then
                 if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
                 if (mod(m1+m2+m3,2) .ne. 0)  erfcterm = 0.0d0
              end if
              term1 = fac1*erfcterm*hhh + expterm*(fac2 + fac3*hsq)
              struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $  k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2 + 
     $  qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,
     $  k3-kstart2(rankloc+1)+1)**2
              e = -(term1 / denom) * struc2
              edsp = edsp + e
              vterm = 3.0d0*(fac1*erfcterm*h + 
     $          fac3*expterm)*struc2/denom
              vir(1,1) = vir(1,1) + h1*h1*vterm - e
              vir(2,1) = vir(2,1) + h1*h2*vterm
              vir(3,1) = vir(3,1) + h1*h3*vterm
              vir(1,2) = vir(1,2) + h2*h1*vterm
              vir(2,2) = vir(2,2) + h2*h2*vterm - e
              vir(3,2) = vir(3,2) + h2*h3*vterm
              vir(1,3) = vir(1,3) + h3*h1*vterm
              vir(2,3) = vir(2,3) + h3*h2*vterm
              vir(3,3) = vir(3,3) + h3*h3*vterm - e
           end if
           qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $      k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1) = 
     $      -(term1/denom) * qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $      k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)
           qgridout_2d(2,k1-istart2(rankloc+1)+1,
     $      k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1) = 
     $      -(term1/denom) * qgridout_2d(2,k1-istart2(rankloc+1)+1,
     $      k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)
 10        continue
          end do
        end do
      end do
c
c     perform the 3-D FFT backward transformation
c
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,commloc,reqbcast(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,commloc,reqbcast(tag),
     $   ierr)
      end do
c
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
c
c     get first derivatives of the reciprocal space energy 
c
      dn1 = dble(nfft1)
      dn2 = dble(nfft2)
      dn3 = dble(nfft3)
      do ii = 1, ndisprecloc
         iidisp = disprecglob(ii)
         iglob = idisp(iidisp)
         iatm = iglob
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         fi = csix(iidisp)
         de1 = 0.0d0
         de2 = 0.0d0
         de3 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-sign(nfft3,k0))/2
            t3 = thetai3(1,it3,ii)
            dt3 = dn3 * thetai3(2,it3,ii)
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-sign(nfft2,j0))/2
               t2 = thetai2(1,it2,ii)
               dt2 = dn2 * thetai2(2,it2,ii)
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-sign(nfft1,i0))/2
                  t1 = thetai1(1,it1,ii)
                  dt1 = dn1 * thetai1(2,it1,ii)
c
                  kstart = kstart1(rankloc+1)
                  kend = kend1(rankloc+1)
                  jstart = jstart1(rankloc+1)
                  jend = jend1(rankloc+1)
                  istart = istart1(rankloc+1)
                  iend = iend1(rankloc+1)

                  if (((k.ge.kstart).and.(k.le.kend)).and.
     $              ((j.ge.jstart).and.(j.le.jend)).and.
     $              ((i.ge.istart).and.(i.le.iend))) then
                    term = qgridin_2d(1,i-istart+1,j-jstart+1,
     $                   k-kstart+1,1)
                    goto 100
                  end if
c
                  do iproc = 1, nrec_send
                    proc = prec_send(iproc)
                    kstart = kstart1(proc+1)
                    kend = kend1(proc+1)
                    jstart = jstart1(proc+1)
                    jend = jend1(proc+1)
                    istart = istart1(proc+1)
                    iend = iend1(proc+1)
                    if (((k.ge.kstart).and.(k.le.kend)).and.
     $                ((j.ge.jstart).and.(j.le.jend)).and.
     $                ((i.ge.istart).and.(i.le.iend))) then
                     term=qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,
     $                 iproc+1)
                      goto 100
                    end if
                   end do
 100               continue

                  de1 = de1 + 2.0d0*term*dt1*t2*t3
                  de2 = de2 + 2.0d0*term*dt2*t1*t3
                  de3 = de3 + 2.0d0*term*dt3*t1*t2
               end do
            end do
         end do
         dedsprec(1,ii) = dedsprec(1,ii) + fi*(recip(1,1)*de1
     &                            +recip(1,2)*de2+recip(1,3)*de3)
         dedsprec(2,ii) = dedsprec(2,ii) + fi*(recip(2,1)*de1
     &                            +recip(2,2)*de2+recip(2,3)*de3)
         dedsprec(3,ii) = dedsprec(3,ii) + fi*(recip(3,1)*de1
     &                            +recip(3,2)*de2+recip(3,3)*de3)
      end do
c
c     account for the energy and virial correction terms
c
      term = csixpr * aewald**3 / denom0
      if (rankloc.eq.0) then
        edsp = edsp - term
        vir(1,1) = vir(1,1) + term
        vir(2,2) = vir(2,2) + term
        vir(3,3) = vir(3,3) + term
      end if
      return
      end
