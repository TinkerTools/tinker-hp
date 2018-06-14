c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chkpole  --  check multipoles at chiral sites  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chkpole" inverts atomic multipole moments as necessary
c     at sites with chiral local reference frame definitions
c
c
      subroutine chkpole
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'mpole.i'
      integer k,ii,iipole
      integer ia,ib,ic,id
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 c1,c2,c3,vol
      logical check
c
c
c     loop over multipole sites testing for chirality inversion
c
      do ii = 1, npoleloc
         iipole = poleglob(ii)
         check = .true.
         if (polaxe(iipole) .ne. 'Z-then-X')  check = .false.
         if (yaxis(iipole) .eq. 0)  check = .false.
         if (check) then
            k = yaxis(iipole)
            ia = ipole(iipole)
            ib = zaxis(iipole)
            ic = xaxis(iipole)
            id = abs(k)
c
c     compute the signed parallelpiped volume at chiral site
c
            xad = x(ia) - x(id)
            yad = y(ia) - y(id)
            zad = z(ia) - z(id)
            xbd = x(ib) - x(id)
            ybd = y(ib) - y(id)
            zbd = z(ib) - z(id)
            xcd = x(ic) - x(id)
            ycd = y(ic) - y(id)
            zcd = z(ic) - z(id)
            c1 = ybd*zcd - zbd*ycd
            c2 = ycd*zad - zcd*yad
            c3 = yad*zbd - zad*ybd
            vol = xad*c1 + xbd*c2 + xcd*c3
c
c     invert atomic multipole components involving the y-axis
c
            if (k.lt.0.and.vol.gt.0.0d0 .or.
     &          k.gt.0.and.vol.lt.0.0d0) then
               yaxis(iipole) = -k
               pole(3,iipole) = -pole(3,iipole)
               pole(6,iipole) = -pole(6,iipole)
               pole(8,iipole) = -pole(8,iipole)
               pole(10,iipole) = -pole(10,iipole)
               pole(12,iipole) = -pole(12,iipole)
            end if
         end if
      end do
      return
      end
