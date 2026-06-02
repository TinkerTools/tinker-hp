!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!
!     "bspline_fill_site" finds B-spline coefficients and derivatives
!     for PME i-th atomic sites along the fractional coordinate axes
!
!
!> @brief 
!> finds B-spline coefficients and derivatives
!> for PME i-th atomic sites along the fractional coordinate axes
!> @param[in] i: atomic index
!> @param[in] impi: atomic index in the "parallel index"
subroutine bspline_fill_site(i,impi)
   use atoms
   use boxes
   use pme
   implicit none
   integer i,ifr,impi
   real*8 xi,yi,zi
   real*8 w,fr,eps
!
!     offset used to shift sites off exact lattice bounds
!
   eps = 1.0d-8
!
!     get the B-spline coefficients for the i-th atomic site
!
   xi = x(i)
   yi = y(i)
   zi = z(i)
   w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
   fr = dble(nfft1) * (w-anint(w)+0.5d0)
   ifr = int(fr-eps)
   w = fr - dble(ifr)
   igrid(1,i) = ifr - bsorder
   call bsplgen (w,thetai1(1,1,impi))
   w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
   fr = dble(nfft2) * (w-anint(w)+0.5d0)
   ifr = int(fr-eps)
   w = fr - dble(ifr)
   igrid(2,i) = ifr - bsorder
   call bsplgen (w,thetai2(1,1,impi))
   w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
   fr = dble(nfft3) * (w-anint(w)+0.5d0)
   ifr = int(fr-eps)
   w = fr - dble(ifr)
   igrid(3,i) = ifr - bsorder
   call bsplgen (w,thetai3(1,1,impi))
   return
end
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine bsplgen  --  B-spline coefficients for an atom  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "bsplgen" gets B-spline coefficients and derivatives for
!     a single PME atomic site along a particular direction
!
!
!> @brief 
!> "bsplgen" gets B-spline coefficients and derivatives for
!> a single PME atomic site along a particular direction
!> for PME i-th atomic sites along the fractional coordinate axes
!> @param[in] w: offset
!> @param[in] thetai: B-spline value
subroutine bsplgen (w,thetai)
   use pme
   use potent
   implicit none
   integer i,j,k
   integer level
   real*8 w,denom
   real*8 thetai(4,*)
   real*8 temp(maxorder,maxorder)
!
!
!     set B-spline depth for partial charges or multipoles
!
   level = 2
   if (use_mpole .or. use_polar)  level = 4
!
!     initialization to get to 2nd order recursion
!
   temp(2,2) = w
   temp(2,1) = 1.0d0 - w
!
!     perform one pass to get to 3rd order recursion
!
   temp(3,3) = 0.5d0 * w * temp(2,2)
   temp(3,2) = 0.5d0 * ((1.0d0+w)*temp(2,1)+(2.0d0-w)*temp(2,2))
   temp(3,1) = 0.5d0 * (1.0d0-w) * temp(2,1)
!
!     compute standard B-spline recursion to desired order
!
   do i = 4, bsorder
      k = i - 1
      denom = 1.0d0 / dble(k)
      temp(i,i) = denom * w * temp(k,k)
      do j = 1, i-2
         temp(i,i-j) = denom * ((w+dble(j))*temp(k,i-j-1)&
         &+(dble(i-j)-w)*temp(k,i-j))
      end do
      temp(i,1) = denom * (1.0d0-w) * temp(k,1)
   end do
!
!     get coefficients for the B-spline first derivative
!
   k = bsorder - 1
   temp(k,bsorder) = temp(k,bsorder-1)
   do i = bsorder-1, 2, -1
      temp(k,i) = temp(k,i-1) - temp(k,i)
   end do
   temp(k,1) = -temp(k,1)
!
!     get coefficients for the B-spline second derivative
!
   if (level .eq. 4) then
      k = bsorder - 2
      temp(k,bsorder-1) = temp(k,bsorder-2)
      do i = bsorder-2, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
      temp(k,bsorder) = temp(k,bsorder-1)
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
!
!     get coefficients for the B-spline third derivative
!
      k = bsorder - 3
      temp(k,bsorder-2) = temp(k,bsorder-3)
      do i = bsorder-3, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
      temp(k,bsorder-1) = temp(k,bsorder-2)
      do i = bsorder-2, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
      temp(k,bsorder) = temp(k,bsorder-1)
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
   end if
!
!     copy coefficients from temporary to permanent storage
!
   do i = 1, bsorder
      do j = 1, level
         thetai(j,i) = temp(bsorder-j+1,i)
      end do
   end do
   return
end
!
!
!       "grid_mpole_site" places the i-th fractional atomic multipole onto
!       the particle mesh Ewald grid
!
!
!> @brief 
!> places the i-th fractional atomic multipole onto
!> the particle mesh Ewald grid
!> @param[in] isite: atomic index
!> @param[in] impi: atomic index in the "parallel index"
!> @param[in] fmp: fractional atomic multipole
subroutine grid_mpole_site(isite,impi,fmp)
   use chunks
   use domdec
   use fft
   use pme
   use potent
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer i,j,k,m,impi,rankloc
   integer iproc,proc
   integer ii,jj,kk
   integer isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer nearpt(3)
   integer abound(6)
   real*8 v0,u0,t0
   real*8 v1,u1,t1
   real*8 v2,u2,t2
   real*8 term0,term1,term2
   real*8 fmp(10)
!
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
!
!       put the permanent multipole moments onto the grid
!
   iatm = isite
   nearpt(1) = igrid(1,iatm) + grdoff
   nearpt(2) = igrid(2,iatm) + grdoff
   nearpt(3) = igrid(3,iatm) + grdoff
   abound(1) = nearpt(1) - nlpts
   abound(2) = nearpt(1) + nrpts
   abound(3) = nearpt(2) - nlpts
   abound(4) = nearpt(2) + nrpts
   abound(5) = nearpt(3) - nlpts
   abound(6) = nearpt(3) + nrpts
   call adjust (offsetx,nfft1,1,abound(1),&
   &abound(2),1,nfft1)
   call adjust (offsety,nfft2,1,abound(3),&
   &abound(4),1,nfft2)
   call adjust (offsetz,nfft3,1,abound(5),&
   &abound(6),1,nfft3)
   do kk = abound(5), abound(6)
      k = kk
      m = k + offsetz
      if (k .lt. 1)  k = k + nfft3
      v0 = thetai3(1,m,impi)
      v1 = thetai3(2,m,impi)
      v2 = thetai3(3,m,impi)
      do jj = abound(3), abound(4)
         j = jj
         m = j + offsety
         if (j .lt. 1)  j = j + nfft2
         u0 = thetai2(1,m,impi)
         u1 = thetai2(2,m,impi)
         u2 = thetai2(3,m,impi)
         term0 = fmp(1)*u0*v0 + fmp(3)*u1*v0&
         &+ fmp(4)*u0*v1 + fmp(6)*u2*v0&
         &+ fmp(7)*u0*v2 + fmp(10)*u1*v1
         term1 = fmp(2)*u0*v0 + fmp(8)*u1*v0&
         &+ fmp(9)*u0*v1
         term2 = fmp(5) * u0 * v0
         do ii = abound(1), abound(2)
            i = ii
            m = i + offsetx
            if (i .lt. 1)  i = i + nfft1
            t0 = thetai1(1,m,impi)
            t1 = thetai1(2,m,impi)
            t2 = thetai1(3,m,impi)
!
            kstart = kstart1(rankloc+1)
            kend = kend1(rankloc+1)
            jstart = jstart1(rankloc+1)
            jend = jend1(rankloc+1)
            istart = istart1(rankloc+1)
            iend = iend1(rankloc+1)
            if (((k.ge.kstart).and.(k.le.kend)).and.&
            &((j.ge.jstart).and.(j.le.jend)).and.&
            &((i.ge.istart).and.(i.le.iend))) then
               qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1) =&
               &qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1)&
               &+ term0*t0 + term1*t1 + term2*t2
               goto 10
            end if
            do iproc = 1, nrec_send
               proc = prec_send(iproc)
               kstart = kstart1(proc+1)
               kend = kend1(proc+1)
               jstart = jstart1(proc+1)
               jend = jend1(proc+1)
               istart = istart1(proc+1)
               iend = iend1(proc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
                  qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                  &iproc+1) = qgridin_2d(1,i-istart+1,j-jstart+1,&
                  &k-kstart+1,iproc+1) + term0*t0 + term1*t1 +&
                  &term2*t2
                  goto 10
               end if
            end do
10          continue
         end do
      end do
   end do
   return
end
!
!       "grid_uind_site" places the fractional induced dipoles onto the
!       particle mesh Ewald grid
!
!
!> @brief 
!> places the fractional induced dipoles onto the
!> particle mesh Ewald grid
!> the particle mesh Ewald grid
!> @param[in] isite: atomic index
!> @param[in] impi: atomic index in the "parallel index"
!> @param[in] fuind: fractional d induced dipoles
!> @param[in] fuind: fractional p induced dipoles
!> @param[in] qgrid2loc: PME grid
subroutine grid_uind_site(isite,impi,fuind,fuinp,qgrid2loc)
   use chunks
   use domdec
   use fft
   use pme
   use potent
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer i,j,k,m,impi
   integer iproc,proc
   integer ii,jj,kk
   integer isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer nearpt(3)
   integer abound(6)
   real*8 v0,u0,t0
   real*8 v1,u1,t1
   real*8 term01,term11
   real*8 term02,term12
   real*8 fuind(3)
   real*8 fuinp(3)
   real*8 qgrid2loc(2,n1mpimax,n2mpimax,n3mpimax,nrec_send+1)
!
!       put the induced dipole moments onto the grid
!
   iatm = isite
   nearpt(1) = igrid(1,iatm) + grdoff
   nearpt(2) = igrid(2,iatm) + grdoff
   nearpt(3) = igrid(3,iatm) + grdoff
   abound(1) = nearpt(1) - nlpts
   abound(2) = nearpt(1) + nrpts
   abound(3) = nearpt(2) - nlpts
   abound(4) = nearpt(2) + nrpts
   abound(5) = nearpt(3) - nlpts
   abound(6) = nearpt(3) + nrpts
   call adjust (offsetx,nfft1,1,abound(1),&
   &abound(2),1,nfft1)
   call adjust (offsety,nfft2,1,abound(3),&
   &abound(4),1,nfft2)
   call adjust (offsetz,nfft3,1,abound(5),&
   &abound(6),1,nfft3)
   do kk = abound(5), abound(6)
      k = kk
      m = k + offsetz
      if (k .lt. 1)  k = k + nfft3
      v0 = thetai3(1,m,impi)
      v1 = thetai3(2,m,impi)
      do jj = abound(3), abound(4)
         j = jj
         m = j + offsety
         if (j .lt. 1)  j = j + nfft2
         u0 = thetai2(1,m,impi)
         u1 = thetai2(2,m,impi)
         term01 = fuind(2)*u1*v0&
         &+ fuind(3)*u0*v1
         term11 = fuind(1)*u0*v0
         term02 = fuinp(2)*u1*v0&
         &+ fuinp(3)*u0*v1
         term12 = fuinp(1)*u0*v0
         do ii = abound(1), abound(2)
            i = ii
            m = i + offsetx
            if (i .lt. 1)  i = i + nfft1
            t0 = thetai1(1,m,impi)
            t1 = thetai1(2,m,impi)
!
            if (use_pmecore) then
               kstart = kstart1(rank_bis+1)
               kend = kend1(rank_bis+1)
               jstart = jstart1(rank_bis+1)
               jend = jend1(rank_bis+1)
               istart = istart1(rank_bis+1)
               iend = iend1(rank_bis+1)
            else
               kstart = kstart1(rank+1)
               kend = kend1(rank+1)
               jstart = jstart1(rank+1)
               jend = jend1(rank+1)
               istart = istart1(rank+1)
               iend = iend1(rank+1)
            end if
            if (((k.ge.kstart).and.(k.le.kend)).and.&
            &((j.ge.jstart).and.(j.le.jend)).and.&
            &((i.ge.istart).and.(i.le.iend))) then
               qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,1) =&
               &qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,1)&
               &+ term01*t0 + term11*t1
               qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,1) =&
               &qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,1)&
               &+ term02*t0 + term12*t1
               goto 10
            end if
            do iproc = 1, nrec_send
               proc = prec_send(iproc)
               kstart = kstart1(proc+1)
               kend = kend1(proc+1)
               jstart = jstart1(proc+1)
               jend = jend1(proc+1)
               istart = istart1(proc+1)
               iend = iend1(proc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
                  qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &= qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,&
                  &iproc+1) + term01*t0 + term11*t1
                  qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,&
                  &iproc+1) = qgrid2loc(2,i-istart+1,j-jstart+1,&
                  &k-kstart+1,iproc+1) + term02*t0&
                  &+ term12*t1
                  goto 10
               end if
            end do
10          continue
!
         end do
      end do
   end do
   return
end
!
!       #################################################################
!       ##                                                             ##
!       ##  subroutine adjust  --  alter site bounds for the PME grid  ##
!       ##                                                             ##
!       #################################################################
!
!
!       "adjust" modifies site bounds on the PME grid and returns
!       an offset into the B-spline coefficient arrays
!
!
!> @brief 
!> modifies site bounds on the PME grid and returns
!> an offset into the B-spline coefficient arrays
!> @param no params
subroutine adjust (offset,nfft,nchk,amin,amax,cmin,cmax)
   implicit none
   integer offset
   integer nfft,nchk
   integer amin,amax
   integer cmin,cmax
!
!
!       modify grid offset and bounds for site at edge of chunk
!
   offset = 0
   if (nchk .ne. 1) then
      if (amin.lt.cmin .or. amax.gt.cmax) then
         if (amin.lt.1 .or. amax.gt.nfft) then
            if (cmin .eq. 1) then
               offset = 1 - amin
               amin = 1
            else if (cmax .eq. nfft) then
               amax = nfft
               amin = amin + nfft
            end if
         else
            if (cmin .gt. amin) then
               offset = cmin - amin
               amin = cmin
            else
               amax = cmax
            end if
         end if
      end if
   end if
   offset = offset + 1 - amin
   return
end
!
!
!       "fphi_mpole_site" extracts the permanent multipole potential on the i-th site from
!       the particle mesh Ewald grid
!
!
!> @brief 
!> "fphi_mpole_site" extracts the permanent multipole potential on the i-th site from
!> the particle mesh Ewald grid
!> @param[in] i: atomic index
!> @param[in] impi: atomic index in the "parallel index"
subroutine fphi_mpole_site(isite,impi)
   use domdec
   use fft
   use pme
   use potent
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer i,j,k,impi,rankloc
   integer iproc,proc
   integer isite,iatm
   integer i0,j0,k0
   integer it1,it2,it3
   integer igrd0,jgrd0,kgrd0
   real*8 v0,v1,v2,v3
   real*8 u0,u1,u2,u3
   real*8 t0,t1,t2,t3,tq
   real*8 tu00,tu10,tu01,tu20,tu11
   real*8 tu02,tu21,tu12,tu30,tu03
   real*8 tuv000,tuv100,tuv010,tuv001
   real*8 tuv200,tuv020,tuv002,tuv110
   real*8 tuv101,tuv011,tuv300,tuv030
   real*8 tuv003,tuv210,tuv201,tuv120
   real*8 tuv021,tuv102,tuv012,tuv111
!
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
!
   iatm = isite
   igrd0 = igrid(1,iatm)
   jgrd0 = igrid(2,iatm)
   kgrd0 = igrid(3,iatm)
   tuv000 = 0.0d0
   tuv001 = 0.0d0
   tuv010 = 0.0d0
   tuv100 = 0.0d0
   tuv200 = 0.0d0
   tuv020 = 0.0d0
   tuv002 = 0.0d0
   tuv110 = 0.0d0
   tuv101 = 0.0d0
   tuv011 = 0.0d0
   tuv300 = 0.0d0
   tuv030 = 0.0d0
   tuv003 = 0.0d0
   tuv210 = 0.0d0
   tuv201 = 0.0d0
   tuv120 = 0.0d0
   tuv021 = 0.0d0
   tuv102 = 0.0d0
   tuv012 = 0.0d0
   tuv111 = 0.0d0
   k0 = kgrd0
   do it3 = 1, bsorder
      k0 = k0 + 1
      k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
      v0 = thetai3(1,it3,impi)
      v1 = thetai3(2,it3,impi)
      v2 = thetai3(3,it3,impi)
      v3 = thetai3(4,it3,impi)
      tu00 = 0.0d0
      tu10 = 0.0d0
      tu01 = 0.0d0
      tu20 = 0.0d0
      tu11 = 0.0d0
      tu02 = 0.0d0
      tu30 = 0.0d0
      tu21 = 0.0d0
      tu12 = 0.0d0
      tu03 = 0.0d0
      j0 = jgrd0
      do it2 = 1, bsorder
         j0 = j0 + 1
         j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
         u0 = thetai2(1,it2,impi)
         u1 = thetai2(2,it2,impi)
         u2 = thetai2(3,it2,impi)
         u3 = thetai2(4,it2,impi)
         t0 = 0.0d0
         t1 = 0.0d0
         t2 = 0.0d0
         t3 = 0.0d0
         i0 = igrd0
         do it1 = 1, bsorder
            i0 = i0 + 1
            i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
!
            tq = 0.0d0
            kstart = kstart1(rankloc+1)
            kend = kend1(rankloc+1)
            jstart = jstart1(rankloc+1)
            jend = jend1(rankloc+1)
            istart = istart1(rankloc+1)
            iend = iend1(rankloc+1)
            if (((k.ge.kstart).and.(k.le.kend)).and.&
            &((j.ge.jstart).and.(j.le.jend)).and.&
            &((i.ge.istart).and.(i.le.iend))) then
               tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1)
               goto 10
            end if
            do iproc = 1, nrec_send
               proc = prec_send(iproc)
               kstart = kstart1(proc+1)
               kend = kend1(proc+1)
               jstart = jstart1(proc+1)
               jend = jend1(proc+1)
               istart = istart1(proc+1)
               iend = iend1(proc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
                  tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                  &iproc+1)
                  goto 10
               end if
            end do
10          continue
!
            t0 = t0 + tq*thetai1(1,it1,impi)
            t1 = t1 + tq*thetai1(2,it1,impi)
            t2 = t2 + tq*thetai1(3,it1,impi)
            t3 = t3 + tq*thetai1(4,it1,impi)
         end do
         tu00 = tu00 + t0*u0
         tu10 = tu10 + t1*u0
         tu01 = tu01 + t0*u1
         tu20 = tu20 + t2*u0
         tu11 = tu11 + t1*u1
         tu02 = tu02 + t0*u2
         tu30 = tu30 + t3*u0
         tu21 = tu21 + t2*u1
         tu12 = tu12 + t1*u2
         tu03 = tu03 + t0*u3
      end do
      tuv000 = tuv000 + tu00*v0
      tuv100 = tuv100 + tu10*v0
      tuv010 = tuv010 + tu01*v0
      tuv001 = tuv001 + tu00*v1
      tuv200 = tuv200 + tu20*v0
      tuv020 = tuv020 + tu02*v0
      tuv002 = tuv002 + tu00*v2
      tuv110 = tuv110 + tu11*v0
      tuv101 = tuv101 + tu10*v1
      tuv011 = tuv011 + tu01*v1
      tuv300 = tuv300 + tu30*v0
      tuv030 = tuv030 + tu03*v0
      tuv003 = tuv003 + tu00*v3
      tuv210 = tuv210 + tu21*v0
      tuv201 = tuv201 + tu20*v1
      tuv120 = tuv120 + tu12*v0
      tuv021 = tuv021 + tu02*v1
      tuv102 = tuv102 + tu10*v2
      tuv012 = tuv012 + tu01*v2
      tuv111 = tuv111 + tu11*v1
   end do
   fphirec(1,impi) = tuv000
   fphirec(2,impi) = tuv100
   fphirec(3,impi) = tuv010
   fphirec(4,impi) = tuv001
   fphirec(5,impi) = tuv200
   fphirec(6,impi) = tuv020
   fphirec(7,impi) = tuv002
   fphirec(8,impi) = tuv110
   fphirec(9,impi) = tuv101
   fphirec(10,impi) = tuv011
   fphirec(11,impi) = tuv300
   fphirec(12,impi) = tuv030
   fphirec(13,impi) = tuv003
   fphirec(14,impi) = tuv210
   fphirec(15,impi) = tuv201
   fphirec(16,impi) = tuv120
   fphirec(17,impi) = tuv021
   fphirec(18,impi) = tuv102
   fphirec(19,impi) = tuv012
   fphirec(20,impi) = tuv111
   return
end
!
!       "fphi_chg_site" extracts the permanent fixed charge potential on the i-th site from
!       the particle mesh Ewald grid
!
!
!> @brief 
!> "fphi_chg_site" extracts the permanent fixed charge potential on the i-th site from
!> the particle mesh Ewald grid
!> @param[in] isite: atomic index
!> @param[in] impi: atomic index in the "parallel index"
subroutine fphi_chg_site(isite,impi)
   use domdec
   use fft
   use pme
   use potent
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer i,j,k,impi,rankloc
   integer iproc,proc
   integer isite,iatm
   integer i0,j0,k0
   integer it1,it2,it3
   integer igrd0,jgrd0,kgrd0
   real*8 v0,v1
   real*8 u0,u1
   real*8 t0,t1,tq
   real*8 tu00,tu10,tu01
   real*8 tuv000,tuv100,tuv010,tuv001
!
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
!
   iatm = isite
   igrd0 = igrid(1,iatm)
   jgrd0 = igrid(2,iatm)
   kgrd0 = igrid(3,iatm)
   tuv000 = 0.0d0
   tuv001 = 0.0d0
   tuv010 = 0.0d0
   tuv100 = 0.0d0
   k0 = kgrd0
   do it3 = 1, bsorder
      k0 = k0 + 1
      k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
      v0 = thetai3(1,it3,impi)
      v1 = thetai3(2,it3,impi)
      tu00 = 0.0d0
      tu10 = 0.0d0
      tu01 = 0.0d0
      j0 = jgrd0
      do it2 = 1, bsorder
         j0 = j0 + 1
         j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
         u0 = thetai2(1,it2,impi)
         u1 = thetai2(2,it2,impi)
         t0 = 0.0d0
         t1 = 0.0d0
         i0 = igrd0
         do it1 = 1, bsorder
            i0 = i0 + 1
            i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
!
            tq = 0.0d0
            kstart = kstart1(rankloc+1)
            kend = kend1(rankloc+1)
            jstart = jstart1(rankloc+1)
            jend = jend1(rankloc+1)
            istart = istart1(rankloc+1)
            iend = iend1(rankloc+1)
            if (((k.ge.kstart).and.(k.le.kend)).and.&
            &((j.ge.jstart).and.(j.le.jend)).and.&
            &((i.ge.istart).and.(i.le.iend))) then
               tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1)
               goto 10
            end if
            do iproc = 1, nrec_send
               proc = prec_send(iproc)
               kstart = kstart1(proc+1)
               kend = kend1(proc+1)
               jstart = jstart1(proc+1)
               jend = jend1(proc+1)
               istart = istart1(proc+1)
               iend = iend1(proc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
                  tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                  &iproc+1)
                  goto 10
               end if
            end do
10          continue
!
            t0 = t0 + tq*thetai1(1,it1,impi)
            t1 = t1 + tq*thetai1(2,it1,impi)
         end do
         tu00 = tu00 + t0*u0
         tu10 = tu10 + t1*u0
         tu01 = tu01 + t0*u1
      end do
      tuv000 = tuv000 + tu00*v0
      tuv100 = tuv100 + tu10*v0
      tuv010 = tuv010 + tu01*v0
      tuv001 = tuv001 + tu00*v1
   end do
   fphirec(1,impi) = tuv000
   fphirec(2,impi) = tuv100
   fphirec(3,impi) = tuv010
   fphirec(4,impi) = tuv001
   return
end
!
!     "fphi_uind_site" extracts the induced dipole potential at the i-th site from
!     the particle mesh Ewald grid
!
!
!> @brief 
!> extracts the induced dipole potential at the i-th site from
!> the particle mesh Ewald grid
!> finds B-spline coefficients and derivatives
!> for PME i-th atomic sites along the fractional coordinate axes
!> @param[in] isite: atomic index
!> @param[in] impi: atomic index in the "parallel index"
!> @param[in] fdip_phi1: fractional coordinates of d induced fields
!> @param[in] fdip_phi2: fractional coordinates of p induced fields
!> @param[in] fdip_phi2: fractional coordinates of total induced fields
subroutine fphi_uind_site(isite,impi,fdip_phi1,fdip_phi2,&
&fdip_sum_phi)
   use domdec
   use fft
   use pme
   use potent
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer i,j,k,impi
   integer iproc,proc
   integer isite,iatm
   integer i0,j0,k0
   integer it1,it2,it3
   integer igrd0,jgrd0,kgrd0
   real*8 v0,v1,v2,v3
   real*8 u0,u1,u2,u3
   real*8 t0,t1,t2,t3
   real*8 t0_1,t0_2,t1_1,t1_2
   real*8 t2_1,t2_2,tq_1,tq_2
   real*8 tu00,tu10,tu01,tu20,tu11
   real*8 tu02,tu30,tu21,tu12,tu03
   real*8 tu00_1,tu01_1,tu10_1
   real*8 tu00_2,tu01_2,tu10_2
   real*8 tu20_1,tu11_1,tu02_1
   real*8 tu20_2,tu11_2,tu02_2
   real*8 tuv100_1,tuv010_1,tuv001_1
   real*8 tuv100_2,tuv010_2,tuv001_2
   real*8 tuv200_1,tuv020_1,tuv002_1
   real*8 tuv110_1,tuv101_1,tuv011_1
   real*8 tuv200_2,tuv020_2,tuv002_2
   real*8 tuv110_2,tuv101_2,tuv011_2
   real*8 tuv000,tuv100,tuv010,tuv001
   real*8 tuv200,tuv020,tuv002,tuv110
   real*8 tuv101,tuv011,tuv300,tuv030
   real*8 tuv003,tuv210,tuv201,tuv120
   real*8 tuv021,tuv102,tuv012,tuv111
   real*8 fdip_phi1(10)
   real*8 fdip_phi2(10)
   real*8 fdip_sum_phi(20)
!
   iatm = isite
   igrd0 = igrid(1,iatm)
   jgrd0 = igrid(2,iatm)
   kgrd0 = igrid(3,iatm)
   tuv100_1 = 0.0d0
   tuv010_1 = 0.0d0
   tuv001_1 = 0.0d0
   tuv200_1 = 0.0d0
   tuv020_1 = 0.0d0
   tuv002_1 = 0.0d0
   tuv110_1 = 0.0d0
   tuv101_1 = 0.0d0
   tuv011_1 = 0.0d0
   tuv100_2 = 0.0d0
   tuv010_2 = 0.0d0
   tuv001_2 = 0.0d0
   tuv200_2 = 0.0d0
   tuv020_2 = 0.0d0
   tuv002_2 = 0.0d0
   tuv110_2 = 0.0d0
   tuv101_2 = 0.0d0
   tuv011_2 = 0.0d0
   tuv000 = 0.0d0
   tuv001 = 0.0d0
   tuv010 = 0.0d0
   tuv100 = 0.0d0
   tuv200 = 0.0d0
   tuv020 = 0.0d0
   tuv002 = 0.0d0
   tuv110 = 0.0d0
   tuv101 = 0.0d0
   tuv011 = 0.0d0
   tuv300 = 0.0d0
   tuv030 = 0.0d0
   tuv003 = 0.0d0
   tuv210 = 0.0d0
   tuv201 = 0.0d0
   tuv120 = 0.0d0
   tuv021 = 0.0d0
   tuv102 = 0.0d0
   tuv012 = 0.0d0
   tuv111 = 0.0d0
   k0 = kgrd0
   do it3 = 1, bsorder
      k0 = k0 + 1
      k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
      v0 = thetai3(1,it3,impi)
      v1 = thetai3(2,it3,impi)
      v2 = thetai3(3,it3,impi)
      v3 = thetai3(4,it3,impi)
      tu00_1 = 0.0d0
      tu01_1 = 0.0d0
      tu10_1 = 0.0d0
      tu20_1 = 0.0d0
      tu11_1 = 0.0d0
      tu02_1 = 0.0d0
      tu00_2 = 0.0d0
      tu01_2 = 0.0d0
      tu10_2 = 0.0d0
      tu20_2 = 0.0d0
      tu11_2 = 0.0d0
      tu02_2 = 0.0d0
      tu00 = 0.0d0
      tu10 = 0.0d0
      tu01 = 0.0d0
      tu20 =0.0d0
      tu11 =0.0d0
      tu02 = 0.0d0
      tu30 = 0.0d0
      tu21 = 0.0d0
      tu12 = 0.0d0
      tu03 = 0.0d0
      j0 = jgrd0
      do it2 = 1, bsorder
         j0 = j0 + 1
         j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
         u0 = thetai2(1,it2,impi)
         u1 = thetai2(2,it2,impi)
         u2 = thetai2(3,it2,impi)
         u3 = thetai2(4,it2,impi)
         t0_1 = 0.0d0
         t1_1 = 0.0d0
         t2_1 = 0.0d0
         t0_2 = 0.0d0
         t1_2 = 0.0d0
         t2_2 = 0.0d0
         t3 = 0.0d0
         i0 = igrd0
         do it1 = 1, bsorder
            i0 = i0 + 1
            i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
!
            tq_1 = 0d0
            tq_2 = 0d0
            if (use_pmecore) then
               kstart = kstart1(rank_bis+1)
               kend = kend1(rank_bis+1)
               jstart = jstart1(rank_bis+1)
               jend = jend1(rank_bis+1)
               istart = istart1(rank_bis+1)
               iend = iend1(rank_bis+1)
            else
               kstart = kstart1(rank+1)
               kend = kend1(rank+1)
               jstart = jstart1(rank+1)
               jend = jend1(rank+1)
               istart = istart1(rank+1)
               iend = iend1(rank+1)
            end if
            if (((k.ge.kstart).and.(k.le.kend)).and.&
            &((j.ge.jstart).and.(j.le.jend)).and.&
            &((i.ge.istart).and.(i.le.iend))) then
               tq_1  = qgrid2in_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
               &1)
               tq_2  = qgrid2in_2d(2,i-istart+1,j-jstart+1,k-kstart+1,&
               &1)
               goto 10
            end if
            do iproc = 1, nrec_send
               proc = prec_send(iproc)
               kstart = kstart1(proc+1)
               kend = kend1(proc+1)
               jstart = jstart1(proc+1)
               jend = jend1(proc+1)
               istart = istart1(proc+1)
               iend = iend1(proc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
                  tq_1  = qgrid2in_2d(1,i-istart+1,j-jstart+1,&
                  &k-kstart+1,iproc+1)
                  tq_2  = qgrid2in_2d(2,i-istart+1,j-jstart+1,&
                  &k-kstart+1,iproc+1)
                  goto 10
               end if
            end do
10          continue
!
            t0_1 = t0_1 + tq_1*thetai1(1,it1,impi)
            t1_1 = t1_1 + tq_1*thetai1(2,it1,impi)
            t2_1 = t2_1 + tq_1*thetai1(3,it1,impi)
            t0_2 = t0_2 + tq_2*thetai1(1,it1,impi)
            t1_2 = t1_2 + tq_2*thetai1(2,it1,impi)
            t2_2 = t2_2 + tq_2*thetai1(3,it1,impi)
            t3 = t3 + (tq_1+tq_2)*thetai1(4,it1,impi)
         end do
         tu00_1 = tu00_1 + t0_1*u0
         tu10_1 = tu10_1 + t1_1*u0
         tu01_1 = tu01_1 + t0_1*u1
         tu20_1 = tu20_1 + t2_1*u0
         tu11_1 = tu11_1 + t1_1*u1
         tu02_1 = tu02_1 + t0_1*u2
         tu00_2 = tu00_2 + t0_2*u0
         tu10_2 = tu10_2 + t1_2*u0
         tu01_2 = tu01_2 + t0_2*u1
         tu20_2 = tu20_2 + t2_2*u0
         tu11_2 = tu11_2 + t1_2*u1
         tu02_2 = tu02_2 + t0_2*u2
         t0 = t0_1 + t0_2
         t1 = t1_1 + t1_2
         t2 = t2_1 + t2_2
         tu00 = tu00 + t0*u0
         tu10 = tu10 + t1*u0
         tu01 = tu01 + t0*u1
         tu20 = tu20 + t2*u0
         tu11 = tu11 + t1*u1
         tu02 = tu02 + t0*u2
         tu30 = tu30 + t3*u0
         tu21 = tu21 + t2*u1
         tu12 = tu12 + t1*u2
         tu03 = tu03 + t0*u3
      end do
      tuv100_1 = tuv100_1 + tu10_1*v0
      tuv010_1 = tuv010_1 + tu01_1*v0
      tuv001_1 = tuv001_1 + tu00_1*v1
      tuv200_1 = tuv200_1 + tu20_1*v0
      tuv020_1 = tuv020_1 + tu02_1*v0
      tuv002_1 = tuv002_1 + tu00_1*v2
      tuv110_1 = tuv110_1 + tu11_1*v0
      tuv101_1 = tuv101_1 + tu10_1*v1
      tuv011_1 = tuv011_1 + tu01_1*v1
      tuv100_2 = tuv100_2 + tu10_2*v0
      tuv010_2 = tuv010_2 + tu01_2*v0
      tuv001_2 = tuv001_2 + tu00_2*v1
      tuv200_2 = tuv200_2 + tu20_2*v0
      tuv020_2 = tuv020_2 + tu02_2*v0
      tuv002_2 = tuv002_2 + tu00_2*v2
      tuv110_2 = tuv110_2 + tu11_2*v0
      tuv101_2 = tuv101_2 + tu10_2*v1
      tuv011_2 = tuv011_2 + tu01_2*v1
      tuv000 = tuv000 + tu00*v0
      tuv100 = tuv100 + tu10*v0
      tuv010 = tuv010 + tu01*v0
      tuv001 = tuv001 + tu00*v1
      tuv200 = tuv200 + tu20*v0
      tuv020 = tuv020 + tu02*v0
      tuv002 = tuv002 + tu00*v2
      tuv110 = tuv110 + tu11*v0
      tuv101 = tuv101 + tu10*v1
      tuv011 = tuv011 + tu01*v1
      tuv300 = tuv300 + tu30*v0
      tuv030 = tuv030 + tu03*v0
      tuv003 = tuv003 + tu00*v3
      tuv210 = tuv210 + tu21*v0
      tuv201 = tuv201 + tu20*v1
      tuv120 = tuv120 + tu12*v0
      tuv021 = tuv021 + tu02*v1
      tuv102 = tuv102 + tu10*v2
      tuv012 = tuv012 + tu01*v2
      tuv111 = tuv111 + tu11*v1
   end do
   fdip_phi1(2) = tuv100_1
   fdip_phi1(3) = tuv010_1
   fdip_phi1(4) = tuv001_1
   fdip_phi1(5) = tuv200_1
   fdip_phi1(6) = tuv020_1
   fdip_phi1(7) = tuv002_1
   fdip_phi1(8) = tuv110_1
   fdip_phi1(9) = tuv101_1
   fdip_phi1(10) = tuv011_1
   fdip_phi2(2) = tuv100_2
   fdip_phi2(3) = tuv010_2
   fdip_phi2(4) = tuv001_2
   fdip_phi2(5) = tuv200_2
   fdip_phi2(6) = tuv020_2
   fdip_phi2(7) = tuv002_2
   fdip_phi2(8) = tuv110_2
   fdip_phi2(9) = tuv101_2
   fdip_phi2(10) = tuv011_2
   fdip_sum_phi(1) = tuv000
   fdip_sum_phi(2) = tuv100
   fdip_sum_phi(3) = tuv010
   fdip_sum_phi(4) = tuv001
   fdip_sum_phi(5) = tuv200
   fdip_sum_phi(6) = tuv020
   fdip_sum_phi(7) = tuv002
   fdip_sum_phi(8) = tuv110
   fdip_sum_phi(9) = tuv101
   fdip_sum_phi(10) = tuv011
   fdip_sum_phi(11) = tuv300
   fdip_sum_phi(12) = tuv030
   fdip_sum_phi(13) = tuv003
   fdip_sum_phi(14) = tuv210
   fdip_sum_phi(15) = tuv201
   fdip_sum_phi(16) = tuv120
   fdip_sum_phi(17) = tuv021
   fdip_sum_phi(18) = tuv102
   fdip_sum_phi(19) = tuv012
   fdip_sum_phi(20) = tuv111
!
   return
end
!
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine cmp_to_fmp  --  transformation of multipoles  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "cmp_to_fmp" transforms the atomic multipoles from Cartesian
!     to fractional coordinates
!
!
!> @brief 
!> transforms the atomic multipoles from Cartesian
!> to fractional coordinates
!> @param[in] cmp: atomic multipoles in cartesian coordinates
!> @param[in] fmp: atomic multipoles in fractional coordinates
subroutine cmp_to_fmp(cmp,fmp)
   use mpole
   implicit none
   integer i,j,k
   real*8 ctf(10,10)
   real*8 cmp(10,*)
   real*8 fmp(10,*)
!
!
!     find the matrix to convert Cartesian to fractional
!
   call cart_to_frac (ctf)
!
!     apply the transformation to get the fractional multipoles
!
   do i = 1, npole
      fmp(1,i) = ctf(1,1) * cmp(1,i)
      do j = 2, 4
         fmp(j,i) = 0.0d0
         do k = 2, 4
            fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
         end do
      end do
      do j = 5, 10
         fmp(j,i) = 0.0d0
         do k = 5, 10
            fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
         end do
      end do
   end do
   return
end
!
!     "cmp_to_fmp_site" transforms the ith atomic multipole from Cartesian
!     to fractional coordinates
!
!
!> @brief 
!> transforms the atomic multipoles from Cartesian
!> to fractional coordinates
!> @param[in] cmp: atomic multipoles in cartesian coordinates
!> @param[in] fmp: atomic multipoles in fractional coordinates
subroutine cmp_to_fmp_site(cmp,fmp)
   use sizes
   implicit none
   integer j,k
   real*8 ctf(10,10)
   real*8 cmp(*)
   real*8 fmp(*)
!
!
!     find the matrix to convert Cartesian to fractional
!
   call cart_to_frac (ctf)
!
!     apply the transformation to get the fractional multipoles
!
   fmp(1) = ctf(1,1) * cmp(1)
   do j = 2, 4
      fmp(j) = 0.0d0
      do k = 2, 4
         fmp(j) = fmp(j) + ctf(j,k)*cmp(k)
      end do
   end do
   do j = 5, 10
      fmp(j) = 0.0d0
      do k = 5, 10
         fmp(j) = fmp(j) + ctf(j,k)*cmp(k)
      end do
   end do
   return
end
!
!
!     ############################################################
!     ##                                                        ##
!     ##  subroutine cart_to_frac  --  Cartesian to fractional  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "cart_to_frac" computes a transformation matrix to convert
!     a multipole object in Cartesian coordinates to fractional
!
!     note the multipole components are stored in the condensed
!     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
!
!
!> @brief 
!> "cart_to_frac" computes a transformation matrix to convert
!> a multipole object in Cartesian coordinates to fractional
!> note the multipole components are stored in the condensed
!> order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
!> @param[out] ctf: transformation matrix
subroutine cart_to_frac (ctf)
   use boxes
   use pme
   implicit none
   integer i,j,k,m
   integer i1,i2
   integer qi1(6)
   integer qi2(6)
   real*8 a(3,3)
   real*8 ctf(10,10)
   data qi1  / 1, 2, 3, 1, 1, 2 /
   data qi2  / 1, 2, 3, 2, 3, 3 /
!
!
!     set the reciprocal vector transformation matrix
!
   do i = 1, 3
      a(1,i) = dble(nfft1) * recip(i,1)
      a(2,i) = dble(nfft2) * recip(i,2)
      a(3,i) = dble(nfft3) * recip(i,3)
   end do
!
!     get the Cartesian to fractional conversion matrix
!
   do i = 1, 10
      do j = 1, 10
         ctf(j,i) = 0.0d0
      end do
   end do
   ctf(1,1) = 1.0d0
   do i = 2, 4
      do j = 2, 4
         ctf(i,j) = a(i-1,j-1)
      end do
   end do
   do i1 = 1, 3
      k = qi1(i1)
      do i2 = 1, 6
         i = qi1(i2)
         j = qi2(i2)
         ctf(i1+4,i2+4) = a(k,i) * a(k,j)
      end do
   end do
   do i1 = 4, 6
      k = qi1(i1)
      m = qi2(i1)
      do i2 = 1, 6
         i = qi1(i2)
         j = qi2(i2)
         ctf(i1+4,i2+4) = a(k,i)*a(m,j) + a(k,j)*a(m,i)
      end do
   end do
   return
end
!
!
!> @brief 
!> transform electrostatic potential (and derivatives) from cartesian to
!fractional coordinatses
!> @param[in] fphi: potential in fractional coordinates
!> @param[in] cphi: potential in cartesian coordinates
subroutine fphi_to_cphi_site(fphi,cphi)
   use boxes
   use sizes
   implicit none
   integer j,k
   real*8 cphi(10)
   real*8 fphi(20)
!
!
!     find the matrix to convert fractional to Cartesian
!
   call frac_to_cart
!
!     apply the transformation to get the Cartesian potential
!
   cphi(1) = ftc(1,1) * fphi(1)
   do j = 2, 4
      cphi(j) = 0.0d0
      do k = 2, 4
         cphi(j) = cphi(j) + ftc(j,k)*fphi(k)
      end do
   end do
   do j = 5, 10
      cphi(j) = 0.0d0
      do k = 5, 10
         cphi(j) = cphi(j) + ftc(j,k)*fphi(k)
      end do
   end do
   return
end
!
!
!
!     ############################################################
!     ##                                                        ##
!     ##  subroutine frac_to_cart  --  fractional to Cartesian  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "frac_to_cart" computes a transformation matrix to convert
!     a multipole object in fraction coordinates to Cartesian
!
!     note the multipole components are stored in the condensed
!     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
!
!
!> @brief 
!> computes a transformation matrix to convert
!> a multipole object in fractional coordinates to cartesian
!> note the multipole components are stored in the condensed
!> order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
!> @param[out] ctf: transformation matrix
subroutine frac_to_cart
!      use sizes
   use boxes
   use pme
   implicit none
   integer i,j,k,m
   integer i1,i2
   integer qi1(6)
   integer qi2(6)
   real*8 a(3,3)
   data qi1  / 1, 2, 3, 1, 1, 2 /
   data qi2  / 1, 2, 3, 2, 3, 3 /
!
!
!     set the reciprocal vector transformation matrix
!
   do i = 1, 3
      a(i,1) = dble(nfft1) * recip(i,1)
      a(i,2) = dble(nfft2) * recip(i,2)
      a(i,3) = dble(nfft3) * recip(i,3)
   end do
!
!     get the fractional to Cartesian conversion matrix
!
   do i = 1, 10
      do j = 1, 10
         ftc(j,i) = 0.0d0
      end do
   end do
   ftc(1,1) = 1.0d0
   do i = 2, 4
      do j = 2, 4
         ftc(i,j) = a(i-1,j-1)
      end do
   end do
   do i1 = 1, 3
      k = qi1(i1)
      do i2 = 1, 3
         i = qi1(i2)
         ftc(i1+4,i2+4) = a(k,i) * a(k,i)
      end do
      do i2 = 4, 6
         i = qi1(i2)
         j = qi2(i2)
         ftc(i1+4,i2+4) = 2.0d0 * a(k,i) * a(k,j)
      end do
   end do
   do i1 = 4, 6
      k = qi1(i1)
      m = qi2(i1)
      do i2 = 1, 3
         i = qi1(i2)
         ftc(i1+4,i2+4) = a(k,i) * a(m,i)
      end do
      do i2 = 4, 6
         i = qi1(i2)
         j = qi2(i2)
         ftc(i1+4,i2+4) = a(k,i)*a(m,j) + a(m,i)*a(k,j)
      end do
   end do
   return
end
!
!     subroutine reassignpme: initialize arrays to compute reciprocal part of
!     electrostatic and polarization interactions, also reassign atoms on 'reciprocal cores'
!     between time steps
!
!> @brief 
!> initialize arrays to compute reciprocal part of
!> electrostatic and polarization interactions, also reassign atoms on 'reciprocal cores'
!> between time steps
!> @param[in] init: first call to the routine
subroutine reassignpme(init)
   use atoms
   use boxes
   use chunks
   use domdec
   use inform
   use iounit
   use fft
   use math
   use pme
   use potent
   use mpi
   implicit none
   real*8 disttemp,jtemp(4),ktemp(4),jtempbis(4),ktempbis(4)
   integer ierr,iglob,j
   integer i,k
   integer iproc,iproc1,temp
   integer proc1,proc2
   integer nprocloc,rankloc,nprocloc1
   integer, allocatable :: req(:)
   integer, allocatable :: counttemp(:),repartrectemp(:)
   integer, allocatable :: domlentemp(:)
   integer, allocatable :: nrecdirrecep(:),precdirrecep(:,:)
   integer, allocatable :: nrecdirreceptemp(:)
   integer, allocatable :: nrecrecep(:),precrecep(:,:)
   integer ifr,kbeg,kend,kmin,kmax,sizek,sizej
   integer jbeg,jend,jmin,jmax
   real*8 wbeg,wend
   integer tag, status(MPI_STATUS_SIZE)
   integer nlocrecold
   real*8 fr,xi,yi,zi,eps,w
   real*8 dist,xr,yr
   logical init
   logical interj,interk
!
   if (deb_Path) write(iout,*), 'reassignpme '
!
!
   allocate (req(nproc*nproc))
   eps = 1.0d-8
   allocate (counttemp(nproc))
   allocate (repartrectemp(n))
   allocate (domlentemp(nproc))
   counttemp = 0
!
!     get geometrical correspondance between processes
!
   if (init) then
      if (use_pmecore) then
         nprocloc = ndir
         rankloc  = rank_bis
      else
         nprocloc = nproc
         rankloc  = rank
      end if
      allocate (nrecdirrecep(nproc))
      allocate (precdirrecep(nproc,nproc))
      allocate (nrecdirreceptemp(nproc))
      allocate (nrecrecep(nproc))
      allocate (precrecep(nproc,nproc))
      nrecdirrecep = 0
      nrecdirreceptemp = 0
      precdirrecep = 0
      nrecrecep = 0
      precrecep = 0
!
      if (allocated(prec_send)) deallocate(prec_send)
      allocate (prec_send(nproc))
      prec_send = 0
      if (allocated(prec_recep)) deallocate(prec_recep)
      allocate (prec_recep(nproc))
      prec_recep = 0
      if (allocated(repartrec)) deallocate(repartrec)
      allocate (repartrec(n))
      if (allocated(domlenrec)) deallocate(domlenrec)
      allocate (domlenrec(nproc))
!
      precdir_recep = 0
      precdir_send = 0
      precdir_recep1 = 0
      precdir_send1 = 0
!
      if (use_pmecore) then
         nprocloc1 = nrec
      else
         nprocloc1 = nproc
      end if
      do iproc = 0, nprocloc1 - 1
         kmin = kstart1(iproc+1)
         kmax = kend1(iproc+1)
         jmin = jstart1(iproc+1)
         jmax = jend1(iproc+1)
         if (kmin.lt.1) kmin = kmin + nfft3
         if (kmax.gt.nfft3) kmax = kmax - nfft3
         if (jmin.lt.1) jmin = jmin + nfft2
         if (jmax.gt.nfft2) jmax = jmax - nfft2

!
         do iproc1 = 0, nprocloc - 1
!
!            wbeg = (zbegproc(iproc1+1))*recip(3,3)
!            wbeg = (cbegproc(iproc1+1))*recip(3,3)
            wbeg = (cbegproc(iproc1+1))/zbox
            fr = dble(nfft3) * (wbeg-anint(wbeg)+0.5d0)
            ifr = int(fr-eps)
            kbeg = ifr - bsorder + grdoff
!            wend = (zendproc(iproc1+1))*recip(3,3)
!            wend = (cendproc(iproc1+1))*recip(3,3)
            wend = (cendproc(iproc1+1))/zbox
            fr = dble(nfft3) * (wend-anint(wend)+0.5d0)
            ifr = int(fr-eps)
            kend = ifr - bsorder + grdoff

            sizek = abs(kend-kbeg)

            if (kbeg.lt.1) kbeg = kbeg + nfft3
            if (kbeg+sizek-1.le.nfft3) then
               kend = kbeg + sizek + 1
            else
               temp = nfft3 - kbeg + 1
               kend = sizek - temp + 1
            end if
!
!            wbeg = (ybegproc(iproc1+1))*recip(2,2)
!            wbeg = (bbegproc(iproc1+1))*recip(2,2)
            wbeg = (bbegproc(iproc1+1))/ybox
            fr = dble(nfft2) * (wbeg-anint(wbeg)+0.5d0)
            ifr = int(fr-eps)
            jbeg = ifr - bsorder + grdoff
!            wend = (yendproc(iproc1+1))*recip(2,2)
!            wend = (bendproc(iproc1+1))*recip(2,2)
            wend = (bendproc(iproc1+1))/ybox
            fr = dble(nfft2) * (wend-anint(wend)+0.5d0)
            ifr = int(fr-eps)
            jend = ifr - bsorder + grdoff

            sizej = abs(jend-jbeg)
            

            if (jbeg.lt.1) jbeg = jbeg + nfft2
            if (jbeg+sizej-1.le.nfft2) then
               jend = jbeg + sizej + 1
            else
               temp = nfft2 - jbeg + 1
               jend = sizej - temp + 1
            end if
!
!    first the "k" segments must intersect
!
            interk = .false.


            if ((kbeg.le.kend).and.(((kmin.ge.kbeg).and.(kmin.le.kend)).or.&
           &  ((kmax.ge.kbeg).and.(kmax.le.kend)).or.&
           &((kbeg.ge.kmin).and.(kbeg.le.kmax)).or.&
           &((kend.ge.kmin).and.(kend.le.kmax))&
           &)) then
             interk = .true.
            end if
            if ((kbeg.gt.kend).and.((kmax.ge.kbeg).or.&
           &  (kmax.le.kend).or.(kmin.ge.kbeg).or.&
           &(kmin.le.kend))) then
             interk = .true.
            end if
!
!    then the "j" segments must intersect
!
            interj = .false. 

            if ((jbeg.le.jend).and.(((jmin.ge.jbeg).and.(jmin.le.jend)).or.&
           &  ((jmax.ge.jbeg).and.(jmax.le.jend)).or.&
           &((jbeg.ge.jmin).and.(jbeg.le.jmax)).or.&
           &((jend.ge.jmin).and.(jend.le.jmax))&
           &)) then
             interj = .true.
            end if
            if ((jbeg.gt.jend).and.((jmax.ge.jbeg).or.&
           &  (kmax.le.jend).or.(jmin.ge.jbeg).or.&
           &(jmin.le.jend))) then
             interj = .true.
            end if


            if (interj.and.interk) then
              nrecdirrecep(iproc+1) = nrecdirrecep(iproc+1)+1
              precdirrecep(nrecdirrecep(iproc+1),iproc+1) = iproc1
            end if


         end do
      end do
      if (use_pmecore) then
         if (rank.gt.ndir-1) then
            nrecdir_recep = nrecdirrecep(rankloc+1)
            precdir_recep = precdirrecep(:,rankloc+1)
         end if
      else
         nrecdir_recep = nrecdirrecep(rankloc+1)
         precdir_recep = precdirrecep(:,rankloc+1)
      end if
!
!   Receive also the neighbors within a cutoff to rotate reciprocal multipoles
!
      if (use_pmecore) then
         nprocloc = ndir
         rankloc  = rank
      else
         nprocloc = nproc
         rankloc  = rank
      end if
!
      nrecdirreceptemp = nrecdirrecep
!
!      store this list (without the processes to rotate the reciprocal multipoles)
!
      nrecdir_recep1 = nrecdir_recep
      precdir_recep1 = precdir_recep
!
!     get the processes to send data to
!
      nrecdir_send1 = 0
!
      if (use_pmecore) then
         if (rank.le.ndir-1) then
            do iproc = 0, nrec-1
               do i = 1, nrecdirrecep(iproc+1)
                  if (precdirrecep(i,iproc+1).eq.rankloc) then
                     nrecdir_send1 = nrecdir_send1 + 1
                     precdir_send1(nrecdir_send1) = iproc+ndir
                  end if
               end do
            end do
         end if
      else
         do iproc = 0, nprocloc-1
            do i = 1, nrecdirrecep(iproc+1)
               if (precdirrecep(i,iproc+1).eq.rankloc) then
                  nrecdir_send1 = nrecdir_send1 + 1
                  precdir_send1(nrecdir_send1) = iproc
               end if
            end do
         end do
      end if
      if (use_pmecore) then
         if (rank.gt.ndir-1) then
            nrecdir_recep = nrecdirreceptemp(rankloc+1-ndir)
            nrecdirrecep = nrecdirreceptemp
            precdir_recep = precdirrecep(:,rankloc+1-ndir)
         end if
      else
         nrecdir_recep = nrecdirreceptemp(rankloc+1)
         precdir_recep = precdirrecep(:,rankloc+1)
      end if
!
!     get the processes to send data to
!
      nrecdir_send = 0
!
      if (use_pmecore) then
         if (rank.le.ndir-1) then
            do iproc = 0, nrec-1
               do i = 1, nrecdirreceptemp(iproc+1)
                  if (precdirrecep(i,iproc+1).eq.rankloc) then
                     nrecdir_send = nrecdir_send + 1
                     precdir_send(nrecdir_send) = iproc+ndir
                  end if
               end do
            end do
         end if
      else
         do iproc = 0, nprocloc-1
            do i = 1, nrecdirreceptemp(iproc+1)
               if (precdirrecep(i,iproc+1).eq.rankloc) then
                  nrecdir_send = nrecdir_send + 1
                  precdir_send(nrecdir_send) = iproc
               end if
            end do
         end do
      end if
!
!      Remove procs already in real space neighbor list and put them in another list
!
      nrecdir_recep2 = 0
      precdir_recep2 = 0
      do i = 1, nrecdir_recep
         proc1 = precdir_recep(i)
         do j = 1, nbig_recep
            proc2 = pbig_recep(j)
            if (proc1.eq.proc2) goto 20
         end do
         nrecdir_recep2 = nrecdir_recep2 + 1
         precdir_recep2(nrecdir_recep2) = proc1
20       continue
      end do
!c
      nrecdir_send2 = 0
      precdir_send2 = 0
      do i = 1, nrecdir_send
         proc1 = precdir_send(i)
         do j = 1, nbig_send
            proc2 = pbig_send(j)
            if (proc1.eq.proc2) goto 30
         end do
         nrecdir_send2 = nrecdir_send2 + 1
         precdir_send2(nrecdir_send2) = proc1
30       continue
      end do
!
!      Remove procs already in recdir space neighbor list and put them in another list
!
      if (use_pmecore) then
         nprocloc = nrec
         rankloc  = rank_bis
      else
         nprocloc = nproc
         rankloc  = rank
      end if
!
      if ((.not.(use_pmecore)).or.(rank.gt.ndir-1)) then
!
!     get the reciprocal processes to send and receive data from
!
!
         nrec_recep = 0
         do iproc = 0, nprocloc-1
            do iproc1 = 0, nprocloc-1
               dist = 1000.0d0
               if (iproc.eq.iproc1) goto 100
               jtemp(1) = jstart1(iproc+1)
               ktemp(1) = kstart1(iproc+1)
               jtemp(2) = jend1(iproc+1)
               ktemp(2) = kstart1(iproc+1)
               jtemp(3) = jend1(iproc+1)
               ktemp(3) = kend1(iproc+1)
               jtemp(4) = jstart1(iproc+1)
               ktemp(4) = kend1(iproc+1)
               jtempbis(1) = jstart1(iproc1+1)
               ktempbis(1) = kstart1(iproc1+1)
               jtempbis(2) = jend1(iproc1+1)
               ktempbis(2) = kstart1(iproc1+1)
               jtempbis(3) = jend1(iproc1+1)
               ktempbis(3) = kend1(iproc1+1)
               jtempbis(4) = jstart1(iproc1+1)
               ktempbis(4) = kend1(iproc1+1)
               do i = 1, 4
                  do j = 1, 4
                     xr = abs(jtemp(i) - jtempbis(j))
                     yr = abs(ktemp(i) - ktempbis(j))
!                     ar = abs(jtemp(i) - jtempbis(j))
!                     br = abs(ktemp(i) - ktempbis(j))
                     if (xr.ge.nfft2-2) xr = nfft2-xr
                     if (yr.ge.nfft3-2) yr = nfft3-yr
!                     if (ar.ge.nfft2-2) ar = nfft2-ar
!                     if (br.ge.nfft3-2) br = nfft3-br
!                     call ftcvec(0d0,ar,br,xr,yr,zr)
                     disttemp = sqrt(xr*xr+yr*yr)
!                     disttemp = sqrt(xr*xr+yr*yr+zr*zr)
                     if (disttemp.le.dist) dist = disttemp
                  end do
               end do
               if ((iproc1.ne.iproc).and.(dist.le.1.5*max(nlpts,nrpts)))&
               &then
                  nrecrecep(iproc+1) = nrecrecep(iproc+1)+1
                  precrecep(nrecrecep(iproc+1),iproc+1) = iproc1
               end if
100            continue
            end do
         end do
         nrec_recep = nrecrecep(rankloc+1)
         prec_recep = precrecep(:,rankloc+1)
!
!     get the processes to send data to
!
         nrec_send = 0
!
         do iproc = 0, nprocloc-1
            do i = 1, nrecrecep(iproc+1)
               if (precrecep(i,iproc+1).eq.rankloc) then
                  nrec_send = nrec_send + 1
                  prec_send(nrec_send) = iproc
               end if
            end do
         end do
!
         if (allocated(qgridin_2d)) deallocate (qgridin_2d)
         allocate (qgridin_2d(2,n1mpimax,n2mpimax,n3mpimax,&
         &nrec_send+1))
         if (allocated(qgrid2in_2d)) deallocate (qgrid2in_2d)
         allocate (qgrid2in_2d(2,n1mpimax,n2mpimax,n3mpimax,&
         &nrec_send+1))
         if (allocated(qgridout_2d)) deallocate (qgridout_2d)
         allocate (qgridout_2d(2,isize2(rankloc+1),jsize2(rankloc+1),&
         &ksize2(rankloc+1)))
         if (allocated(qgrid2out_2d)) deallocate (qgrid2out_2d)
         allocate (qgrid2out_2d(2,isize2(rankloc+1),jsize2(rankloc+1),&
         &ksize2(rankloc+1)))
      end if
!
      domlenrec = 0
      do i = 1, n
         xi = x(i)
         yi = y(i)
         zi = z(i)
         w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr = dble(nfft3) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         k = ifr- bsorder + grdoff
         if (k .lt. 1)  k = k + nfft3
         w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr = dble(nfft2) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         j = ifr- bsorder + grdoff
         if (j .lt. 1)  j = j + nfft2
         do iproc = 0, nprocloc-1
            if (((k.ge.kstart1(iproc+1)).and.(k.le.kend1(iproc+1))).and.&
            &((j.ge.jstart1(iproc+1)).and.(j.le.jend1(iproc+1)))) then
               repartrec(i) = iproc
            end if
         end do
      end do
      domlentemp = domlen
      call orderbuffer(.false.)
      call orderbufferrec
!
      deallocate (precdirrecep)
      deallocate (nrecdirrecep)
      deallocate (nrecdirreceptemp)
      deallocate (precrecep)
      deallocate (nrecrecep)
   end if
!
!
   if (use_pmecore) then
      nprocloc = nrec
      rankloc  = rank_bis
   else
      nprocloc = nproc
      rankloc  = rank
   end if
!
   nlocrecold = nlocrec
   repartrec = -1
   do i = 1, nloc
      iglob = glob(i)
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
      fr = dble(nfft3) * (w-anint(w)+0.5d0)
      ifr = int(fr-eps)
      k = ifr- bsorder + grdoff
      if (k .lt. 1)  k = k + nfft3
      w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
      fr = dble(nfft2) * (w-anint(w)+0.5d0)
      ifr = int(fr-eps)
      j = ifr- bsorder + grdoff
      if (j .lt. 1)  j = j + nfft2
      do iproc = 0, nprocloc-1
         if (((k.ge.kstart1(iproc+1)).and.(k.le.kend1(iproc+1))).and.&
         &((j.ge.jstart1(iproc+1)).and.(j.le.jend1(iproc+1)))) then
            repartrec(i) = iproc
         end if
      end do
   end do
!
!     send and receive repartrec
!
   do i = 1, nrecdir_recep1
      if (precdir_recep1(i).ne.rank) then
         tag = nproc*rank + precdir_recep1(i) + 1
         call MPI_IRECV(repartrec(bufbeg(precdir_recep1(i)+1)),&
         &domlen(precdir_recep1(i)+1),MPI_INT,precdir_recep1(i),tag,&
         &COMM_TINKER,req(tag),ierr)
      end if
   end do
!
   do i = 1, nrecdir_send1
      if (precdir_send1(i).ne.rank) then
         tag = nproc*precdir_send1(i) + rank + 1
         call MPI_ISEND(repartrec,&
         &nloc,MPI_INT,precdir_send1(i),tag,COMM_TINKER,&
         &req(tag),ierr)
      end if
   end do
!
   do i = 1, nrecdir_recep1
      if (precdir_recep1(i).ne.rank) then
         tag = nproc*rank + precdir_recep1(i) + 1
         call MPI_WAIT(req(tag),status,ierr)
      end if
   end do
   do i = 1, nrecdir_send1
      if (precdir_send1(i).ne.rank) then
         tag = nproc*precdir_send1(i) + rank + 1
         call MPI_WAIT(req(tag),status,ierr)
      end if
   end do
!
!     reorder repartrec
!
   call imove(nblocrecdir,repartrec,repartrectemp)
   do i = 1, nloc
      iglob = glob(i)
      repartrec(iglob) =  repartrectemp(i)
   end do
   do iproc = 1, nrecdir_recep1
      if (precdir_recep1(iproc).ne.rank) then
         do i = 1, domlen(precdir_recep1(iproc)+1)
            iglob = glob(bufbeg(precdir_recep1(iproc)+1)+i-1)
            repartrec(iglob) =&
            &repartrectemp(bufbeg(precdir_recep1(iproc)+1)+i-1)
         end do
      end if
   end do
   call orderbufferrec
!
   deallocate (domlentemp)
   deallocate (counttemp)
   deallocate (repartrectemp)
   deallocate (req)
   return
end
!
!     "grid_pchg_site" places the i-th fractional atomic charge onto
!     the particle mesh Ewald grid
!
!
!> @brief 
!> places the i-th fractional atomic charge onto
!> the particle mesh Ewald grid
!> @param[in] isite: atomic index
!> @param[in] impi: atomic index in the "parallel index"
!> @param[in] q: charge
subroutine grid_pchg_site(isite,impi,q)
   use chunks
   use domdec
   use fft
   use pme
   use potent
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer i,j,k,m,impi,rankloc
   integer iproc,proc
   integer ii,jj,kk
   integer isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer nearpt(3)
   integer abound(6)
   real*8 q
   real*8 v0,u0,t0
   real*8 term
!
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
!
!     put the permanent multipole moments onto the grid
!
   iatm = isite
   nearpt(1) = igrid(1,iatm) + grdoff
   nearpt(2) = igrid(2,iatm) + grdoff
   nearpt(3) = igrid(3,iatm) + grdoff
   abound(1) = nearpt(1) - nlpts
   abound(2) = nearpt(1) + nrpts
   abound(3) = nearpt(2) - nlpts
   abound(4) = nearpt(2) + nrpts
   abound(5) = nearpt(3) - nlpts
   abound(6) = nearpt(3) + nrpts
   call adjust (offsetx,nfft1,1,abound(1),&
   &abound(2),1,nfft1)
   call adjust (offsety,nfft2,1,abound(3),&
   &abound(4),1,nfft2)
   call adjust (offsetz,nfft3,1,abound(5),&
   &abound(6),1,nfft3)
   do kk = abound(5), abound(6)
      k = kk
      m = k + offsetz
      if (k .lt. 1)  k = k + nfft3
      v0 = thetai3(1,m,impi)
      do jj = abound(3), abound(4)
         j = jj
         m = j + offsety
         if (j .lt. 1)  j = j + nfft2
         u0 = thetai2(1,m,impi)
         term = q*u0*v0
         do ii = abound(1), abound(2)
            i = ii
            m = i + offsetx
            if (i .lt. 1)  i = i + nfft1
            t0 = thetai1(1,m,impi)
!
            kstart = kstart1(rankloc+1)
            kend = kend1(rankloc+1)
            jstart = jstart1(rankloc+1)
            jend = jend1(rankloc+1)
            istart = istart1(rankloc+1)
            iend = iend1(rankloc+1)
            if (((k.ge.kstart).and.(k.le.kend)).and.&
            &((j.ge.jstart).and.(j.le.jend)).and.&
            &((i.ge.istart).and.(i.le.iend))) then
               qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1) =&
               &qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1)&
               &+ term*t0
               goto 10
            end if
            do iproc = 1, nrec_send
               proc = prec_send(iproc)
               kstart = kstart1(proc+1)
               kend = kend1(proc+1)
               jstart = jstart1(proc+1)
               jend = jend1(proc+1)
               istart = istart1(proc+1)
               iend = iend1(proc+1)
               if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
                  qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                  &iproc+1) = qgridin_2d(1,i-istart+1,j-jstart+1,&
                  &k-kstart+1,iproc+1) + term*t0
                  goto 10
               end if
            end do
10          continue
         end do
      end do
   end do
   return
end
!
!     #######################################################################
!     ##                                                                   ##
!     ##  subroutine grid_disp_site  --  put dispersion sites on PME grid  ##
!     ##                                                                   ##
!     #######################################################################
!
!
!     "grid_disp_site" places the i-th damped dispersion coefficients onto
!     the particle mesh Ewald grid
!
!
!> @brief 
!> places the i-th damped dispersion coefficients onto
!> the particle mesh Ewald grid
!> @param[in] i: atomic index
!> @param[in] impi: atomic index in the "parallel index"
subroutine grid_disp_site(isite,impi)
   use atoms
   use disp
   use domdec
   use chunks
   use fft
   use pme
   use potent
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer istart3,iend3,jstart3,jend3,kstart3,kend3
   integer i,j,k,m,impi,rankloc
   integer ii,jj,kk
   integer iproc,proc
   integer isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer nearpt(3)
   integer abound(6)
   integer cbound(6)
   real*8 v0,u0,t0
   real*8 term

   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
   kstart = kstart1(rankloc+1)
   kend = kend1(rankloc+1)
   jstart = jstart1(rankloc+1)
   jend = jend1(rankloc+1)
   istart = istart1(rankloc+1)
   iend = iend1(rankloc+1)
!
   iatm = isite
   nearpt(1) = igrid(1,iatm) + grdoff
   nearpt(2) = igrid(2,iatm) + grdoff
   nearpt(3) = igrid(3,iatm) + grdoff
   abound(1) = nearpt(1) - nlpts
   abound(2) = nearpt(1) + nrpts
   abound(3) = nearpt(2) - nlpts
   abound(4) = nearpt(2) + nrpts
   abound(5) = nearpt(3) - nlpts
   abound(6) = nearpt(3) + nrpts
   call adjust (offsetx,nfft1,1,abound(1),&
   &abound(2),cbound(1),cbound(2))
   call adjust (offsety,nfft2,1,abound(3),&
   &abound(4),cbound(3),cbound(4))
   call adjust (offsetz,nfft3,1,abound(5),&
   &abound(6),cbound(5),cbound(6))
   do kk = abound(5), abound(6)
      k = kk
      m = k + offsetz
      if (k .lt. 1)  k = k + nfft3
      v0 = thetai3(1,m,impi) * csix(isite)
      do jj = abound(3), abound(4)
         j = jj
         m = j + offsety
         if (j .lt. 1)  j = j + nfft2
         u0 = thetai2(1,m,impi)
         term = v0 * u0
         do ii = abound(1), abound(2)
            i = ii
            m = i + offsetx
            if (i .lt. 1)  i = i + nfft1
            t0 = thetai1(1,m,impi)
!
            if (((k.ge.kstart).and.(k.le.kend)).and.&
            &((j.ge.jstart).and.(j.le.jend)).and.&
            &((i.ge.istart).and.(i.le.iend))) then
               qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1) =&
               &qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,1)&
               &+ term*t0
               goto 10
            end if

            do iproc = 1, nrec_send
               proc = prec_send(iproc)
               kstart3 = kstart1(proc+1)
               kend3 = kend1(proc+1)
               jstart3 = jstart1(proc+1)
               jend3 = jend1(proc+1)
               istart3 = istart1(proc+1)
               iend3 = iend1(proc+1)
               if (((k.ge.kstart3).and.(k.le.kend3)).and.&
               &((j.ge.jstart3).and.(j.le.jend3)).and.&
               &((i.ge.istart3).and.(i.le.iend3))) then
                  qgridin_2d(1,i-istart3+1,j-jstart3+1,k-kstart3+1,&
                  &iproc+1) = qgridin_2d(1,i-istart3+1,j-jstart3+1,&
                  &k-kstart3+1,iproc+1) + term*t0
                  goto 10
               end if
            end do
10          continue

         end do
      end do
   end do
!
   return
end
