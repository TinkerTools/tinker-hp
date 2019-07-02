c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
      function sprod(n,v,w)
      implicit real*8 (a-h,o-z)
      dimension v(*), w(*)
      sprod=0.d0
      do j=1,n
        sprod=sprod+v(j)*w(j)
      enddo
      return
      end
c
      subroutine aadd(n,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*)
      do j = 1, n
        c(j) = a(j) + b(j)
      enddo
      return
      end
c
      subroutine asub(n,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*)
      do j = 1, n
        c(j) = a(j) - b(j)
      enddo
      return
      end
c
      subroutine acasb(n,a,b,c,s)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*)
      do j = 1, n
        c(j) = a(j) + s*b(j)
      enddo
      return
      end
c
      subroutine acasb2(n,a,b1,b2,c,s1,s2)
      implicit real*8 (a-h,o-z)
      dimension a(*),b1(*),b2(*),c(*)
      do j = 1, n
        c(j) = a(j) + s1*b1(j) + s2*b2(j)
      enddo
      return
      end
c
      subroutine amove(n,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*)
      do j = 1, n
        b(j) = a(j)
      enddo
      return
      end
c
      subroutine imove(n,a,b)
      implicit integer (a-z)
      dimension a(*),b(*)
      do j = 1, n
        b(j) = a(j)
      enddo
      return
      end
c
      subroutine aneg(n,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*)
      do j = 1, n
        b(j) = -a(j)
      enddo
      return
      end
c
      subroutine aset2(n,s,a)
      implicit real*8 (a-h,o-z)
      dimension a(*)
      do j = 1, n
        a(j) = s
      enddo
      return
      end
c
      subroutine asasbc(n,a,b,c,s1,s2)
      implicit real*8 (a-h,o-z)
      dimension a(*), b(*), c(*)
      do j = 1, n
        c(j) = s1*a(j) + s2*b(j)
      enddo
      return
      end
c
      subroutine ascale(n,s,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(*), b(*)
      do j = 1, n
        b(j) = s*a(j)
      enddo
      return
      end
c
      subroutine iscale(n,s,a,b)
      implicit integer (a-z)
      dimension a(*), b(*)
      do j = 1, n
        b(j) = s*a(j)
      enddo
      return
      end
c
      integer function iarsum(istart,n,ia)
      implicit integer(a-z)
      dimension ia(*)
c
      sum = 0
      do 100 i = 1, n
  100   sum = sum + ia(istart+i-1)
      iarsum = sum
      return
      end
c
      function arrsum(n,a)
      implicit real*8(a-h,o-z)
      dimension a(*)
      arrsum = 0.d0
      do 100 j = 1, n
 100    arrsum = arrsum + a(j)
      return
      end
c
      subroutine aclear(n,a)
      implicit real*8 (a-h,o-z)
      dimension a(*)
      save zero
      data zero/0.0d0/
c
      do 10 j = 1, n
 10     a(j) = zero
      return
      end
c
      subroutine iclear(n,a)
      implicit integer (a-z)
      dimension a(*)
      save zero
      data zero/0/
c
      do 10 j = 1, n
 10     a(j) = zero
      return
      end
c
      subroutine ione(n,a)
      implicit integer (a-z)
      dimension a(*)
      save one
      data one/1/
c
      do 10 j = 1, n
 10     a(j) = one
      return
      end
c
      subroutine rmsvec(n,v,vrms,vmax)
      implicit real*8(a-h,o-z)
      dimension v(*)
      save zero
      data zero/0.0d0/
c
      vmax=zero
      vrms=zero
      do 100 j = 1, n
        if(abs(v(j)).gt.vmax) vmax=abs(v(j))
 100    vrms = vrms + (v(j)**2)
      vrms=sqrt(vrms/dble(n))
      return
      end
c
      subroutine ylmbas(lmax,x,vfac,basis,vplm,vcos,vsin)
      use math
      implicit none
c
c   calculates a basis of real spherical harmonics up to order lmax.
c   the routine computes as a first the generalized legendre polynomials
c   and then weighs them with sine and cosine coefficients.
c   this lasts are evaluated starting from cos(phi) = y/sqrt(1-z**2)
c   using chebyshev polynomials.
c
c   input:  lmax  ... maximum angular momentum of the spherical harmonics
c                     basis
c           x     ... unit vector (x \in s_2) where to evaluate the basis
c           vfac  ... scratch array dimensioned 2*nbasis. Contains the
c                     normalization factors for the spherical harmonics.
c
c   output: basis ... spherical harmonics at x. basis is dimensioned (lmax+1)^2
c                     and the harmonics are ordered for increasing l and m
c                     (i.e.: 0,0; 1,-1; 1,0; 1,1; 2,-2; ...)
c
c   scratch: vplm ... scratch array dimensioned (lmax + 1)^2. it hosts the
c                     generalized legendre polynomials.
c            vcos ... scratch array dimensioned lmax + 1. it hosts the cos(m*\phi)
c                     values
c            vsin ... scratch array dimensioned lmax + 1. it hosts the sin(m*\phi)
c                     values.
c
      integer lmax
      real*8 basis(*), x(*), vfac(2,*), vcos(*), vsin(*), vplm(*)
c
      integer l, m, ind
      real*8  zero, one, cthe, sthe, cphi, sphi, plm
      save zero, one
      data zero/0.0d0/, one/1.0d0/
c
c   get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
c   coordinates of x.
c
      cthe = x(3)
      sthe = sqrt(one - cthe*cthe)
      if(sthe.ne.zero) then
        cphi = x(1)/sthe
        sphi = x(2)/sthe
      else
        cphi = one
        sphi = zero
      endif
c
c   evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is
c   pointless if z = 1, as the only non vanishing terms will be the
c   ones with m=0.
c
      if(sthe.ne.zero) then
        call trgev(lmax,cphi,sphi,vcos,vsin)
      else
        do m = 1, lmax + 1
          vcos(m) = one
          vsin(m) = zero
        end do
      end if
c
c   evaluate the generalized legendre polynomials.
c
      call polleg(lmax,cthe,sthe,vplm)
c
c   now build the spherical harmonics. we will distinguish m=0,
c   m>0 and m<0:
c
      do l = 0, lmax
        ind = l**2 + l + 1
c       m = 0
        basis(ind) = vfac(1,ind)*vplm(ind)
        do m = 1, l
          plm = vplm(ind+m)
c         m > 0
          basis(ind+m) = vfac(1,ind+m)*plm*vcos(m+1)
c         m < 0
          basis(ind-m) = vfac(2,ind-m)*plm*vsin(m+1)
        end do
      end do
      return
      end
c
      subroutine dbasis(lmax,x,vfac,basis,dbas,vplm,vcos,vsin)
      use math
      implicit none
c
c   calculates a basis of real spherical harmonics up to order lmax.
c   the routine computes as a first the generalized legendre polynomials
c   and then weighs them with sine and cosine coefficients.
c   this lasts are evaluated starting from cos(phi) = y/sqrt(1-z**2)
c   using chebyshev polynomials.
c
c   input:  lmax  ... maximum angular momentum of the spherical harmonics
c                     basis
c           x     ... unit vector (x \in s_2) where to evaluate the basis
c           vfac  ... scratch array dimensioned lmax + 1. it hosts the first 2*lmax
c                     + 1 factorials. this is to be filled somewhere at the
c                     beginning of the program, as it is pointless to recompute
c                     it all the times. notice that 0! = vfac(1) according to
c                     fortran conventions.
c
c   output: basis ... spherical harmonics at x. basis is dimensioned (lmax+1)^2
c                     and the harmonics are ordered for increasing l and m
c                     (i.e.: 0,0; 1,-1; 1,0; 1,1; 2,-2; ...)
c           dbas  ... spherical harmonics cartesian gradient.
c
c   scratch: vplm ... scratch array dimensioned (lmax + 1)^2. it hosts the
c                     generalized legendre polynomials.
c            vcos ... scratch array dimensioned lmax + 1. it hosts the cos(m*\phi)
c                     values
c            vsin ... scratch array dimensioned lmax + 1. it hosts the sin(m*\phi)
c                     values.
c
      integer lmax
      real*8 basis(*), x(*), dbas(3,*), vfac(2,*), vcos(*), vsin(*),
     $  vplm(*)
      integer l, m, ind
      real*8  zero, pt5, one, cthe, sthe, cphi, sphi, et(3), ep(3), fln,
     $  dbthe, fnorm, plm, pp1, pm1, pp, dbphi
      save zero, pt5, one
      data zero/0.0d0/, pt5/0.5d0/, one/1.0d0/
c
c     get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
c     coordinates of x.
c
      cthe = x(3)
      sthe = sqrt(one - cthe*cthe)
      if(sthe.ne.zero) then
        cphi = x(1)/sthe
        sphi = x(2)/sthe
      else
        cphi = one
        sphi = zero
      endif
c
c     evaluate the dirivatives of theta and phi:
c
      et(1) = cthe*cphi
      et(2) = cthe*sphi
      et(3) = -sthe
      if(sthe.ne.zero) then
        ep(1) = -sphi/sthe
        ep(2) = cphi/sthe
        ep(3) = zero
      else
        ep(1) = zero
        ep(2) = zero
        ep(3) = zero
      endif
c
c     evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is
c     pointless if z = 1, as the only non vanishing terms will be the
c     ones with m=0.
c
      if(sthe.ne.zero) then
        call trgev(lmax,cphi,sphi,vcos,vsin)
      else
        do m = 1, lmax + 1
          vcos(m) = one
          vsin(m) = zero
        end do
      end if
c
c     evaluate the generalized legendre polynomials.
c
      call polleg(lmax,cthe,sthe,vplm)
c
c     now build the spherical harmonics. we will distinguish m=0,
c     m>0 and m<0:
c
      do l = 0, lmax
        ind = l**2 + l + 1
c
c     functions for m = 0
c
        fln = vfac(1,ind)
        basis(ind) = fln*vplm(ind)
c
c     gradient for m = 0
c
        if(l.gt.0) then
          dbthe = fln*vplm(ind+1)
          dbas(1,ind) = dbthe*et(1)
          dbas(2,ind) = dbthe*et(2)
          dbas(3,ind) = dbthe*et(3)
        else
          dbas(1,ind) = zero
          dbas(2,ind) = zero
          dbas(3,ind) = zero
        endif
c
        do m = 1, l
          fnorm = vfac(1,ind+m)
c
c     get p_l^m for the function, p_l^{m+1} and p_l^{m-1} for the
c     gradient
c
          plm = fnorm*vplm(ind+m)
          pp1 = zero
          if(m.lt.l) pp1 = -pt5*vplm(ind+m+1)
          pm1 = pt5*(dble(l+m)*dble(l-m+1)*vplm(ind+m-1))
          pp  = pp1 + pm1
c
c     m > 0
c
          dbthe = -fnorm*pp*vcos(m+1)
          dbphi = -dble(m)*fnorm*vplm(ind+m)*vsin(m+1)
          basis(ind+m)  = plm*vcos(m+1)
          dbas(1,ind+m) = dbthe*et(1) + dbphi*ep(1)
          dbas(2,ind+m) = dbthe*et(2) + dbphi*ep(2)
          dbas(3,ind+m) = dbthe*et(3)
c
c     m < 0
c
          dbthe = -fnorm*pp*vsin(m+1)
          dbphi = dble(m)*fnorm*vplm(ind+m)*vcos(m+1)
          basis(ind-m) = plm*vsin(m+1)
          dbas(1,ind-m) = dbthe*et(1) + dbphi*ep(1)
          dbas(2,ind-m) = dbthe*et(2) + dbphi*ep(2)
          dbas(3,ind-m) = dbthe*et(3)
        end do
      end do
      return
      end
c
      subroutine trgev(n,x,y,cx,sx)
      implicit none
      integer n
      real*8  x, y, cx(*), sx(*)
      integer i
      real*8  zero, one, two
      save zero, one, two
      data zero/0.0d0/, one/1.0d0/, two/2.0d0/
c
      cx(1) = one
      cx(2) = x
      sx(1) = zero
      sx(2) = y
      do i = 3, n+1
        cx(i) = two*x*cx(i-1) - cx(i-2)
        sx(i) = two*x*sx(i-1) - sx(i-2)
      end do
      return
      end
c
      subroutine polleg(lmax,x,y,plm)
      implicit none
c
c   computes the l,m associated legendre polynomial for -1 <= x <= 1
c   using the recurrence formula
c
c   (l-m)p(l,m) = x(2l-1)p(l-1,m) - (l+m-1)p(l-2,m)
c
      integer lmax
      real*8  x, y, plm(*)
      integer m, ind, ind2, l
      real*8 one, two, fact, pmm, somx2, pmm1, pmmo, pll
      save one, two
      data one/1.0d0/, two/2.0d0/
c
      fact  = one
      pmm   = one
      somx2 = y
c     somx2 = sqrt((one-x)*(one+x))
      do m = 0, lmax
c
c  do pmm for each m:
c
        ind  = (m+1)*(m+1)
        plm(ind) = pmm
        if(m.eq.lmax) return
        pmm1 = x*(two*dble(m)+one)*pmm
        ind2 = ind + 2*m + 2
        plm(ind2) = pmm1
        pmmo = pmm
        do l = m+2, lmax
          pll   = (x*(two*dble(l)-one)*pmm1 -
     $           (dble(l)+dble(m)-1)*pmm)/(dble(l)-dble(m))
          ind = l*l + l + 1
          plm(ind+m) = pll
          pmm  = pmm1
          pmm1 = pll
        end do
        pmm  = -pmmo*fact*somx2
        fact = fact + two
      end do
      return
      end
c
      subroutine derrot(ii,doder,i,iz,ix,iy,r,dri,driz,drix,driy)
      use atoms
      use mpole
      implicit none
c
c     given the index ii of a multipole, gives as output the indexes of
c     the atoms used to define the local frame and the derivatives of
c     the rotation matrix with respect to the position of such atoms.
c
      integer ii, i, iz, ix, iy
      real*8  r(3,3), dri(3,3,3), driz(3,3,3), drix(3,3,3), driy(3,3,3)
      logical doder
      integer j, k, l
      real*8  xxi, xxi2, eet, eet2, zze, zze2, uu, uu2, uu3, 
     $  vez, eex, eex2, zero, one
      real*8  xi(3), eta(3), zeta(3), u(3), v(3), ez(3), ex(3), ey(3)
      real*8  ezu(3,3), exv(3,3), exez(3,3), eyez(3,3), eyex(3,3),
     $  uri(3,3), uriz(3,3), urix(3,3), uriy(3,3), vri(3,3), vrix(3,3),
     $  vriy(3,3), exi(3,3), exiz(3,3), exix(3,3), exiy(3,3), eyi(3,3),
     $  eyiz(3,3), eyix(3,3), eyiy(3,3), ezi(3,3), eziz(3,3), ezix(3,3),
     $  eziy(3,3)
c
      save zero, one
      data zero/0.0d0/, one/1.0d0/
c
 1000 format('polaxe not recognized or not implemented.')
c
c     Get the reference atoms coordinates.
c
      i  = ipole(ii)
c
      if (polaxe(i).eq.'None') then
        r    = zero
        dri  = zero
        driz = zero
        drix = zero
        driy = zero
        return
      end if
c
      iz = zaxis(i)
      ix = xaxis(i)
      iy = yaxis(i)
      if (iz.ne.0) then
        xi(1)   = x(iz) - x(i)
        xi(2)   = y(iz) - y(i)
        xi(3)   = z(iz) - z(i)
        xxi2    = xi(1)*xi(1) + xi(2)*xi(2) + xi(3)*xi(3)
        xxi     = sqrt(xxi2)
      end if
      if (ix.ne.0) then
        eta(1)  = x(ix) - x(i)
        eta(2)  = y(ix) - y(i)
        eta(3)  = z(ix) - z(i)
        eet2    = eta(1)*eta(1) + eta(2)*eta(2) + eta(3)*eta(3)
        eet     = sqrt(eet2)
      end if
      if (iy.ne.0) then
        zeta(1) = x(iy) - x(i)
        zeta(2) = y(iy) - y(i)
        zeta(3) = z(iy) - z(i)
        zze2    = zeta(1)*zeta(1) + zeta(2)*zeta(2) + zeta(3)*zeta(3)
        zze     = sqrt(zze2)
      end if
c
c     We will write ez = u/|u|, ex = (v - (v*ez)ez)/(v - (v*ez)ez) and
c     ey = ez x ex for all the possible definitions.
c     Here, we define for each method the u and v vectors:
c
      if (polaxe(i).eq.'Z-then-X') then
        u(1) = xi(1)
        u(2) = xi(2)
        u(3) = xi(3)
        v(1) = eta(1)
        v(2) = eta(2)
        v(3) = eta(3)
      else if (polaxe(i).eq.'Bisector') then
        u(1) = eet*xi(1) + xxi*eta(1)
        u(2) = eet*xi(2) + xxi*eta(2)
        u(3) = eet*xi(3) + xxi*eta(3)
        v(1) = eta(1)
        v(2) = eta(2)
        v(3) = eta(3)
      else if (polaxe(i).eq.'Z-Bisect') then
        u(1) = xi(1)
        u(2) = xi(2)
        u(3) = xi(3)
        v(1) = zze*eta(1) + eet*zeta(1)
        v(2) = zze*eta(2) + eet*zeta(2)
        v(3) = zze*eta(3) + eet*zeta(3)
      else if (polaxe(i).eq.'3-Fold') then
        u(1) = eet*zze*xi(1) + xxi*zze*eta(1) + eet*zze*zeta(1)
        u(2) = eet*zze*xi(2) + xxi*zze*eta(2) + eet*zze*zeta(2)
        u(3) = eet*zze*xi(3) + xxi*zze*eta(3) + eet*zze*zeta(3)
        v(1) = eta(1)
        v(2) = eta(2)
        v(3) = eta(3)
      else
        write(6,1000)
        call fatal
      end if
      uu2 = u(1)*u(1) + u(2)*u(2) + u(3)*u(3)
      uu  = sqrt(uu2)
c
c     Assemble the three versors:
c
      ez(1) = u(1)/uu
      ez(2) = u(2)/uu
      ez(3) = u(3)/uu
      vez   = ez(1)*v(1) + ez(2)*v(2) + ez(3)*v(3)
      ex(1) = v(1) - vez*ez(1)
      ex(2) = v(2) - vez*ez(2)
      ex(3) = v(3) - vez*ez(3)
      eex2  = ex(1)*ex(1) + ex(2)*ex(2) + ex(3)*ex(3)
      eex   = sqrt(eex2)
      ex(1) = ex(1)/eex
      ex(2) = ex(2)/eex
      ex(3) = ex(3)/eex
      ey(1) = ez(2)*ex(3) - ez(3)*ex(2)
      ey(2) = ez(3)*ex(1) - ez(1)*ex(3)
      ey(3) = ez(1)*ex(2) - ez(2)*ex(1)
      do j = 1, 3
        r(j,1) = ex(j)
        r(j,2) = ey(j)
        r(j,3) = ez(j)
      end do
      if (.not. doder) return
c
c     clear everything.
c
      do j = 1, 3
        do k = 1, 3
          ezu(j,k)  = zero
          exv(j,k)  = zero
          exez(j,k) = zero
          eyez(j,k) = zero
          eyex(j,k) = zero
          uri(j,k)  = zero
          uriz(j,k) = zero
          urix(j,k) = zero
          uriy(j,k) = zero
          vri(j,k)  = zero
          vrix(j,k) = zero
          vriy(j,k) = zero
          exi(j,k)  = zero
          exiz(j,k) = zero
          exix(j,k) = zero
          exiy(j,k) = zero
          eyi(j,k)  = zero
          eyiz(j,k) = zero
          eyix(j,k) = zero
          eyiy(j,k) = zero
          ezi(j,k)  = zero
          eziz(j,k) = zero
          ezix(j,k) = zero
          eziy(j,k) = zero
          do l = 1, 3
            dri(j,k,l)  = zero
            driz(j,k,l) = zero
            drix(j,k,l) = zero
            driy(j,k,l) = zero
          end do
        end do
      end do
c
c     We will assemble the derivatives of the rotation matrices as
c     the product of two contributions, according to the chain rule.
c     The first part is the derivative of the versors wrt the u, v
c     vectors and is independent of the model; the second part is the
c     derivatives of u and v wrt the positions of the involved atoms
c     and is specific. We will assemble here the first part. For later
c     convenience we will also compute here dex/dez, dey/dez and
c     dey/dex.
c
      uu3  = uu*uu2
      do j = 1, 3
        ezu(j,j)  = one/uu
        exv(j,j)  = one/eex
        exez(j,j) = -vez/eex
        do k = 1, 3
          ezu(j,k)  = ezu(j,k) - u(j)*u(k)/uu3
          exv(j,k)  = exv(j,k) - ez(j)*ez(k)/eex - ex(j)*ex(k)/eex
          exez(j,k) = exez(j,k) + ex(j)*vez*v(k)/eex2 - ez(j)*v(k)/eex
        end do
      end do
      eyez(1,1) = zero
      eyez(1,2) =  ex(3)
      eyez(1,3) = -ex(2)
      eyez(2,1) = -ex(3)
      eyez(2,2) = zero
      eyez(2,3) =  ex(1)
      eyez(3,1) =  ex(2)
      eyez(3,2) = -ex(1)
      eyez(3,3) = zero
      eyex(1,1) = zero
      eyex(1,2) = -ez(3)
      eyex(1,3) =  ez(2)
      eyex(2,1) =  ez(3)
      eyex(2,2) = zero
      eyex(2,3) = -ez(1)
      eyex(3,1) = -ez(2)
      eyex(3,2) =  ez(1)
      eyex(3,3) = zero
c
c     We compute now all the chain rule contributions, in particular:
c       du/dr(i), du/dr(iz), du/dr(ix), du/dr(iy)
c       dv/dr(i), dv/dr(ix), dv/dr(iy)
c     which we will use to assemble the derivatives of the rotation
c     matrix.
c
      if (polaxe(i).eq.'Z-then-X') then
        do j = 1, 3
          uri(j,j)  = -one
          uriz(j,j) = one
          vri(j,j)  = -one
          vrix(j,j) = one
        end do
      else if (polaxe(i).eq.'Bisector') then
        do j = 1, 3
          uri(j,j)  = -xxi - eet
          uriz(j,j) = eet
          urix(j,j) = xxi
          vri(j,j)  = -one
          vrix(j,j) = one
          do k = 1, 3
            uri(j,k)  = uri(j,k) - eta(j)*xi(k)/xxi - xi(j)*eta(k)/eet
            uriz(j,k) = uriz(j,k) + eta(j)*xi(k)/xxi
            urix(j,k) = urix(j,k) + xi(j)*eta(k)/eet
          end do
        end do
      else if (polaxe(i).eq.'Z-Bisect') then
        do j = 1, 3
          uri(j,j)  = -one
          uriz(j,j) = one
          vri(j,j)  = -eet - zze
          vrix(j,j) = zze
          vriy(j,j) = eet
          do k = 1, 3
            vri(j,k)  = vri(j,k)  - eta(j)*zeta(k)/zze
     $                            - zeta(j)*eta(k)/eet
            vrix(j,k) = vrix(j,k) + zeta(j)*eta(k)/eet
            vriy(j,k) = vriy(j,k) + eta(j)*zeta(k)/zze
          end do
        end do
      else if (polaxe(i).eq.'3-Fold') then
        do j = 1, 3
          uri(j,j)  = -xxi*eet - xxi*zze - eet*zze
          uriz(j,j) = +eet*zze
          urix(j,j) = +xxi*zze
          uriy(j,j) = +xxi*eet
          vri(j,j)  = -one
          vrix(j,j) = one
          do k = 1, 3
            uri(j,k)  = uri(j,k)  - (zze*eta(j)+eet*zeta(j))*xi(k)/xxi
     $                            - (zze*xi(j)+xxi*zeta(j))*eta(k)/eet
     $                            - (eet*xi(j)+xxi*eta(j))*zeta(k)/zze
            uriz(j,k) = uriz(j,k) + (zze*eta(j)+eet*zeta(j))*xi(k)/xxi
            urix(j,k) = urix(j,k) + (zze*xi(j)+xxi*zeta(j))*eta(k)/eet
            uriy(j,k) = uriy(j,k) + (eet*xi(j)+xxi*eta(j))*zeta(k)/zze
          end do
        end do
      else
        write(6,1000)
        call fatal
      end if
c
c     We are, finally, ready to assemble the derivatives of each versor.
c     The code does not make distinctions anymore between the various
c     methods to define the local frame. This is not top efficient, but
c     the price to pay for a more efficient code would be a very messy
c     one: as the overall cost of the computation is small, order is
c     preferred.
c
c     ez derivatives:
c
      do l = 1, 3
        do j = 1, 3
          do k = 1, 3
            ezi(l,j)  = ezi(l,j)  + ezu(l,k)*uri(k,j)
            eziz(l,j) = eziz(l,j) + ezu(l,k)*uriz(k,j)
            ezix(l,j) = ezix(l,j) + ezu(l,k)*urix(k,j)
            eziy(l,j) = eziy(l,j) + ezu(l,k)*uriy(k,j)
          end do
        end do
      end do
c
c     ex derivatives:
c
      do l = 1, 3
        do j = 1, 3
          do k = 1, 3
            exi(l,j) =exi(l,j) +exv(l,k)*vri(k,j) +exez(l,k)*ezi(k,j)
            exiz(l,j)=exiz(l,j)                   +exez(l,k)*eziz(k,j)
            exix(l,j)=exix(l,j)+exv(l,k)*vrix(k,j)+exez(l,k)*ezix(k,j)
            exiy(l,j)=exiy(l,j)+exv(l,k)*vriy(k,j)+exez(l,k)*eziy(k,j)
          end do
        end do
      end do
c
c     ey derivatives:
c
      do l = 1, 3
        do j = 1, 3
          do k = 1, 3
            eyi(l,j) =eyi(l,j) +eyex(l,k)*exi(k,j) +eyez(l,k)*ezi(k,j)
            eyiz(l,j)=eyiz(l,j)+eyex(l,k)*exiz(k,j)+eyez(l,k)*eziz(k,j)
            eyix(l,j)=eyix(l,j)+eyex(l,k)*exix(k,j)+eyez(l,k)*ezix(k,j)
            eyiy(l,j)=eyiy(l,j)+eyex(l,k)*exiy(k,j)+eyez(l,k)*eziy(k,j)
          end do
        end do
      end do
c
c     Finally, assemble the derivatives of the rotation matrix.
c
      do j = 1, 3
        do k = 1, 3
          dri(k,j,1)  = exi(j,k)
          dri(k,j,2)  = eyi(j,k)
          dri(k,j,3)  = ezi(j,k)
          driz(k,j,1) = exiz(j,k)
          driz(k,j,2) = eyiz(j,k)
          driz(k,j,3) = eziz(j,k)
          drix(k,j,1) = exix(j,k)
          drix(k,j,2) = eyix(j,k)
          drix(k,j,3) = ezix(j,k)
          driy(k,j,1) = exiy(j,k)
          driy(k,j,2) = eyiy(j,k)
          driy(k,j,3) = eziy(j,k)
        end do
      end do
      return
      end
c
      logical function isnan(x)
      real*8 x
      isnan = .false.
      return
      end
