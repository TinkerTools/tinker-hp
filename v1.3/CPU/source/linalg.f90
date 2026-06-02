!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!> @brief 
!> computes scalar product between two vectors
!> @param[in] n: dimension of the vectors
!> @param[in] v: first vector
!> @param[in] w: second vector
real*8 function sprod(n,v,w)
   integer j,n
   real*8, dimension(*) :: v, w
   sprod=0.d0
   do j=1,n
      sprod=sprod+v(j)*w(j)
   enddo
   return
end
!
!> @brief 
!> computes sum of two vectors
!> @param[in] n: dimension of the vectors
!> @param[in] a: first vector
!> @param[in] b: second vector
!> @param[in] c: result of the sum
subroutine aadd(n,a,b,c)
   integer j,n
   real*8, dimension(*) :: a,b,c
   do j = 1, n
      c(j) = a(j) + b(j)
   enddo
   return
end
!
!
!> @brief 
!> copy a vector of reals in an other one
!> @param[in] n: dimension of the vectors
!> @param[in] a: first vector
!> @param[in] b: second vector
subroutine amove(n,a,b)
   integer j,n
   real*8, dimension(*) :: a,b
   do j = 1, n
      b(j) = a(j)
   enddo
   return
end
!
!> @brief 
!> copy a vector of integer in an other one
!> @param[in] n: dimension of the vectors
!> @param[in] a: first vector
!> @param[in] b: second vector
subroutine imove(n,a,b)
   integer j,n
   integer, dimension(*) :: a,b
   do j = 1, n
      b(j) = a(j)
   enddo
   return
end
!
!> @brief 
!> set values of a vector of reals to zero
!> @param[in] n: dimension of the vector
!> @param[in] a: vector
subroutine aclear(n,a)
   integer j,n
   real*8 zero
   real*8, dimension(*) :: a
   save zero
   data zero/0.0d0/
!
   do j = 1, n
      a(j) = zero
   enddo
   return
end
!
!> @brief 
!> set values of a vector of integers to zero
!> @param[in] n: dimension of the vector
!> @param[in] a: vector
subroutine iclear(n,a)
   integer j,n
   integer zero
   integer, dimension(*) :: a
   save zero
   data zero/0/
!
   do j = 1, n
      a(j) = zero
   enddo
   return
end
!
!> @brief 
!> given the index ii of a multipole, gives as output the indexes of
!> the atoms used to define the local frame and the derivatives of
!> the rotation matrix with respect to the position of such atoms.
!> @param[in] ii: index of a multipole
!> @param[in] doder: flag governing computation of derivatives
!> @param[in] i: index of the atomic site
!> @param[in] iz: index of the "iz" atom for local frame definition
!> @param[in] ix: index of the "ix" atom for local frame definition
!> @param[in] iy: index of the "iy" atom for local frame definition
!> @param[in] r: rotation matrix
!> @param[in] dri: derivatives of rotation matrix wrt i coordinates
!> @param[in] driz: derivatives of rotation matrix wrt iz coordinates
!> @param[in] drix: derivatives of rotation matrix wrt ix coordinates
!> @param[in] driy: derivatives of rotation matrix wrt iy coordinates
subroutine derrot(ii,doder,i,iz,ix,iy,r,dri,driz,drix,driy)
   use atoms
   use mpole
   implicit none
!
!     given the index ii of a multipole, gives as output the indexes of
!     the atoms used to define the local frame and the derivatives of
!     the rotation matrix with respect to the position of such atoms.
!
   integer ii, i, iz, ix, iy
   real*8  r(3,3), dri(3,3,3), driz(3,3,3), drix(3,3,3), driy(3,3,3)
   logical doder
   integer j, k, l
   real*8  xxi, xxi2, eet, eet2, zze, zze2, uu, uu2, uu3,&
   &vez, eex, eex2, zero, one
   real*8  xi(3), eta(3), zeta(3), u(3), v(3), ez(3), ex(3), ey(3)
   real*8  ezu(3,3), exv(3,3), exez(3,3), eyez(3,3), eyex(3,3),&
   &uri(3,3), uriz(3,3), urix(3,3), uriy(3,3), vri(3,3), vrix(3,3),&
   &vriy(3,3), exi(3,3), exiz(3,3), exix(3,3), exiy(3,3), eyi(3,3),&
   &eyiz(3,3), eyix(3,3), eyiy(3,3), ezi(3,3), eziz(3,3), ezix(3,3),&
   &eziy(3,3)
!
   save zero, one
   data zero/0.0d0/, one/1.0d0/
!
1000 format('polaxe not recognized or not implemented.')
!
!     Get the reference atoms coordinates.
!
   i  = ipole(ii)
!
   if (polaxe(ii).eq.'None') then
      r    = zero
      dri  = zero
      driz = zero
      drix = zero
      driy = zero
      return
   end if
!
   iz = zaxis(ii)
   if (iz.le.0) iz = i
   ix = xaxis(ii)
   if (ix.le.0) ix = i
   iy = yaxis(ii)
   if (iy.le.0) iy = i

   if (iz.ne.0) then
      xi(1)   = x(iz) - x(ii)
      xi(2)   = y(iz) - y(ii)
      xi(3)   = z(iz) - z(ii)
      xxi2    = xi(1)*xi(1) + xi(2)*xi(2) + xi(3)*xi(3)
      xxi     = sqrt(xxi2)
   end if
   if (ix.ne.0) then
      eta(1)  = x(ix) - x(ii)
      eta(2)  = y(ix) - y(ii)
      eta(3)  = z(ix) - z(ii)
      eet2    = eta(1)*eta(1) + eta(2)*eta(2) + eta(3)*eta(3)
      eet     = sqrt(eet2)
   end if
   if (iy.ne.0) then
      zeta(1) = x(iy) - x(ii)
      zeta(2) = y(iy) - y(ii)
      zeta(3) = z(iy) - z(ii)
      zze2    = zeta(1)*zeta(1) + zeta(2)*zeta(2) + zeta(3)*zeta(3)
      zze     = sqrt(zze2)
   end if
!
!     We will write ez = u/|u|, ex = (v - (v*ez)ez)/(v - (v*ez)ez) and
!     ey = ez x ex for all the possible definitions.
!     Here, we define for each method the u and v vectors:
!
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
!
!     Assemble the three versors:
!
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
!
!     clear everything.
!
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
!
!     We will assemble the derivatives of the rotation matrices as
!     the product of two contributions, according to the chain rule.
!     The first part is the derivative of the versors wrt the u, v
!     vectors and is independent of the model; the second part is the
!     derivatives of u and v wrt the positions of the involved atoms
!     and is specific. We will assemble here the first part. For later
!     convenience we will also compute here dex/dez, dey/dez and
!     dey/dex.
!
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
!
!     We compute now all the chain rule contributions, in particular:
!       du/dr(i), du/dr(iz), du/dr(ix), du/dr(iy)
!       dv/dr(i), dv/dr(ix), dv/dr(iy)
!     which we will use to assemble the derivatives of the rotation
!     matrix.
!
   if (polaxe(ii).eq.'Z-then-X') then
      do j = 1, 3
         uri(j,j)  = -one
         uriz(j,j) = one
         vri(j,j)  = -one
         vrix(j,j) = one
      end do
   else if (polaxe(ii).eq.'Bisector') then
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
   else if (polaxe(ii).eq.'Z-Bisect') then
      do j = 1, 3
         uri(j,j)  = -one
         uriz(j,j) = one
         vri(j,j)  = -eet - zze
         vrix(j,j) = zze
         vriy(j,j) = eet
         do k = 1, 3
            vri(j,k)  = vri(j,k)  - eta(j)*zeta(k)/zze&
            &- zeta(j)*eta(k)/eet
            vrix(j,k) = vrix(j,k) + zeta(j)*eta(k)/eet
            vriy(j,k) = vriy(j,k) + eta(j)*zeta(k)/zze
         end do
      end do
   else if (polaxe(ii).eq.'3-Fold') then
      do j = 1, 3
         uri(j,j)  = -xxi*eet - xxi*zze - eet*zze
         uriz(j,j) = +eet*zze
         urix(j,j) = +xxi*zze
         uriy(j,j) = +xxi*eet
         vri(j,j)  = -one
         vrix(j,j) = one
         do k = 1, 3
            uri(j,k)  = uri(j,k)  - (zze*eta(j)+eet*zeta(j))*xi(k)/xxi&
            &- (zze*xi(j)+xxi*zeta(j))*eta(k)/eet&
            &- (eet*xi(j)+xxi*eta(j))*zeta(k)/zze
            uriz(j,k) = uriz(j,k) + (zze*eta(j)+eet*zeta(j))*xi(k)/xxi
            urix(j,k) = urix(j,k) + (zze*xi(j)+xxi*zeta(j))*eta(k)/eet
            uriy(j,k) = uriy(j,k) + (eet*xi(j)+xxi*eta(j))*zeta(k)/zze
         end do
      end do
   else
      write(6,1000)
      call fatal
   end if
!
!     We are, finally, ready to assemble the derivatives of each versor.
!     The code does not make distinctions anymore between the various
!     methods to define the local frame. This is not top efficient, but
!     the price to pay for a more efficient code would be a very messy
!     one: as the overall cost of the computation is small, order is
!     preferred.
!
!     ez derivatives:
!
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
!
!     ex derivatives:
!
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
!
!     ey derivatives:
!
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
!
!     Finally, assemble the derivatives of the rotation matrix.
!
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
!
logical function isnan(x)
   real*8 x
   isnan = .false.
   return
end
