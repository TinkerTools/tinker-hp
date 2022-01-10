c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine dispersion  --  dispersion energy                   ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "dispersion" calculates the repulsion energy
c
c
      subroutine edispersion
      use bound
      implicit none
      if (use_replica) then
        call edispersion0a
      else
        call edispersion0b
      end if
      return
      end
c
c     "dispersion0b" calculates the dispersion energy
c     
      subroutine edispersion0b
      use analyz
      use atmtyp
      use atmlst
      use atoms
      use bond
      use bound
      use cell
      use charge
      use chargetransfer
      use cutoff
      use dispersion
      use domdec
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mpole
      use neigh
      use shunt
      use usage
      use potent
      use katoms
      use repulsion
      implicit none
c
      integer i,j,k,iglob,jglob,ilp,jlp
      integer ilpnl,jj,jnl
      real*8 xi,yi,zi,xj,yj,zj,rij2,rij
      real*8 rij3,rij4,rij5
      real*8 xij,yij,zij
      real*8 rdispi,rdispj,mij,mij6,mij8,mij10
      real*8 bij,rdpi,rdpj,dij
      real*8 rij6,rij8,rij10
      real*8 dmp6,dmp8,dmp10,disp6,disp8,disp10
      real*8 ci,cj,gdi,gdj,tij,edisp
      real*8 xlp,ylp,zlp,xk,yk,zk
      real*8 xklp,yklp,zklp,rklp2,rklp
      real*8 xdisp,rsij
      real*8 expo6
      real*8 edisp1,xdisp1,exdisp1
      real*8 edisp2,xdisp2,exdisp2
      real*8 edisp3,xdisp3,exdisp3
      real*8 xkj,ykj,zkj,rkj2,rkj 
      real*8 rkj3,rkj4,rkj5
      real*8 bkj,bkjx,rdpk,dkj,rskj
      real*8 mkj,mkj6,mkj8,mkj10
      real*8 clp,gdk,tkj,clp1,clp2
      integer l
      real*8 xl,yl,zl
      real*8 xkl,ykl,zkl,rkl2,rkl
      real*8 bkl,rdpl,dkl,rskl
      real*8 mkl,mkl6,mkl8,mkl10
      real*8 ck,cl,gdl,tkl
      real*8 taper,rlp
      real*8 e
      logical docompute
      character*6 mode
c
      call rotlp
c
      disp6 = 0.0d0
      disp8 = 0.0d0
      disp10 = 0.0d0
      exdisp = 0.0d0
      xdisp = 0.0d0
      edisp = 0.0d0
      edisp1 = 0.0d0
      edisp2 = 0.0d0
      edisp3 = 0.0d0
      xdisp1 = 0.0d0
      xdisp2 = 0.0d0
      xdisp3 = 0.0d0
      exdisp1 = 0.0d0
      exdisp2 = 0.0d0
      exdisp3 = 0.0d0
c
c     set the coefficients for the switching function
c
      mode = 'DISP'
      call switch (mode)

c
c     atomic dispersion computation
c
      do i = 1, nlocnl
         iglob = ineignl(i)
         xi = x(iglob) 
         yi = y(iglob) 
         zi = z(iglob) 
         do j = 1, natatlst(i)
            jglob = atatlst(j,i)
            if (molcule(iglob) .ne. molcule(jglob)) then
               xj = x(jglob) 
               yj = y(jglob) 
               zj = z(jglob) 
               xij = xj - xi 
               yij = yj - yi 
               zij = zj - zi 
               if (use_bounds)  call image (xij,yij,zij)
               rij2 = xij*xij + yij*yij + zij*zij
c
               if (rij2.gt.(dispcut*dispcut)) cycle
c
               rij = sqrt(rij2) 
               bij =  gorbrep(iglob)*gorbrep(jglob)
               rdpi = vdwdisp1(iglob)
               rdpj = vdwdisp1(jglob)
               dij = ((rdpi+rdpj)*bdmp/rij)-1
               rsij = (rdpi + rdpj)*bdmp
               expo6 = admp6*dij
c              si dump pas assez efficace , > 1
               if (rij.ge.rsij) then
                 dmp6 = bij
                 dmp8 = bij
                 dmp10 = bij
               else
                 dmp6 = bij*exp(-expo6)
                 dmp8 = bij*exp(-admp8*dij)
                 dmp10 = bij*exp(-admp10*dij)
               end if
               mij = rij/(2*sqrt(rdpi*rdpj))
               mij6 = mij*mij*mij*mij*mij*mij 
               mij8 = mij6*mij*mij
               mij10 = mij8*mij*mij
               disp6 = c6disp*dmp6/mij6
               disp8 = c8disp*dmp8/mij8
               disp10 = c10disp*dmp10/mij10
c              
c                 exchange disp              
c
               ci = rpole(1,iglob) 
               cj = rpole(1,jglob)
               gdi = 1 - ci/valemtp(atomic(iglob)) 
               gdj = 1 - cj/valemtp(atomic(jglob)) 
               tij = bij*gdi*gdj
               xdisp1 = tij*cxd*exp(-axd*mij)
               e = (disp6 + disp8 + disp10)*scdp*
     $         facdispij*discof
               e = e + xdisp1*facdispij
c
c     use energy switching if near the cutoff distance
c
               if (rij2 .gt. cut2) then
                  rij3 = rij2 * rij
                  rij4 = rij2 * rij2
                  rij5 = rij2 * rij3
                  taper = c5*rij5 + c4*rij4 + c3*rij3
     &                       + c2*rij2 + c1*rij + c0
                  e = e * taper
                end if
               exdisp = exdisp + e
            end if
         end do
      end do  
c
c     lp-atom dispersion
c
      do ilpnl = 1, nlplocnl
         ilp = lpglobnl(ilpnl)
         k = lpatom(ilp)
         xlp = rlonepair(1,ilp)
         ylp = rlonepair(2,ilp)
         zlp = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do jj = 1, nlpatlst(ilpnl)
            j = lpatlst(jj,ilpnl)
            if (molcule(k) .ne. molcule(j)) then
               disp6 =0.0d0
               disp8 = 0.0d0
               disp10 = 0.0d0
               dmp6 =0.0d0
               dmp8 = 0.0d0
               dmp10 = 0.0d0
c               
               xj = x(j) 
               yj = y(j) 
               zj = z(j) 
               xkj = xj - xlp
               ykj = yj - ylp
               zkj = zj - zlp
               if (use_bounds)  call image (xkj,ykj,zkj)
               rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
c
               if (rkj2.gt.(dispcut*dispcut)) cycle
c
               rkj = sqrt(rkj2)
               bkj = gorbrep(j)
               bkjx = gorbrep(k)*gorbrep(j)
c
c   vdw radius is the vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp2(k)! + dincr_lprep(ilp)
               rdpj = vdwdisp1(j)
               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
               rskj = (rdpk + rdpj)*bdmp
               expo6 = -admp6*dkj
               if (rkj > rskj) expo6=0.0d0
               dmp6 = bkj*exp(expo6)
               dmp8 = bkj*exp(-admp8*dkj)
               dmp10 = bkj*exp(-admp10*dkj)
               mkj = rkj/(2*sqrt(rdpk*rdpj))
               mkj6 = mkj*mkj*mkj*mkj*mkj*mkj 
               mkj8 = mkj6*mkj*mkj
               mkj10 = mkj8*mkj*mkj
               clp = lpcharge(ilp)!2.0d0
               disp6 = -0.5*clp*colpa*c6disp*dmp6/mkj6
               disp8 = -0.5*clp*colpa*c8disp*dmp8/mkj8
               disp10 = -0.5*clp*colpa*c10disp*dmp10/mkj10
c
c              exchange disp
c
               cj = rpole(1,j)
               ck = rpole(1,k)
               gdk = 1 - ck/valemtp(atomic(k)) 
               gdj = 1 - cj/valemtp(atomic(j)) 
               tkj = bkjx*gdk*gdj
               e = (disp6 + disp8 + disp10)*facdispij*
     $            discof
               xdisp2 = 0.5*clp*tkj*cxdla*exp(-axdla*mkj)
               e = e + xdisp2*facdispij
c
c     use energy switching if near the cutoff distance
c
               if (rkj2 .gt. cut2) then
                  rkj3 = rkj2 * rkj
                  rkj4 = rkj2 * rkj2
                  rkj5 = rkj2 * rkj3
                  taper = c5*rkj5 + c4*rkj4 + c3*rkj3
     &                       + c2*rkj2 + c1*rkj + c0
                  e = e * taper
                end if
               exdisp = exdisp + e
            end if   
         end do 
      end do  
c
c      lone pair / lone pair dispersion term 
c
      do ilpnl = 1, nlplocnl
         ilp = lpglobnl(ilpnl)
         k = lpatom(ilp)    
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, nlplplst(ilpnl)
            jlp = lplplst(j,ilpnl)
            l = lpatom(jlp)
            if (molcule(k) .ne. molcule(l)) then
               disp6 =0.0d0
               disp8 = 0.0d0
               disp10 = 0.0d0
               xj = rlonepair(1,jlp)
               yj = rlonepair(2,jlp)
               zj = rlonepair(3,jlp)
               xl = x(l) 
               yl = y(l) 
               zl = z(l) 
               xij = xj - xi 
               yij = yj - yi 
               zij = zj - zi 
               if (use_bounds)  call image (xij,yij,zij)
               rij2 = xij*xij + yij*yij + zij*zij
c
               if (rij2.gt.(dispcut*dispcut)) cycle
c
               rij = sqrt(rij2)
               bkl = 1  
c
c    vdw radius is vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp3(k)! + dincr_lprep(ilp) 
               rdpl = vdwdisp3(l)! + dincr_lprep(jlp) 
c
               dkl = ((rdpk+rdpl)*bdmp/rij)-1
               rskl = (rdpk + rdpl)*bdmp
               expo6 = -admp6*dkl
               if (rij > rskl) expo6 = 0.0d0 
               dmp6 = bkl*exp(expo6)
               dmp8 = bkl*exp(-admp8*dkl)
               dmp10 = bkl*exp(-admp10*dkl)
               mkl = rij/(2*sqrt(rdpk*rdpl))
               mkl6 = mkl*mkl*mkl*mkl*mkl*mkl 
               mkl8 = mkl6*mkl*mkl
               mkl10 = mkl8*mkl*mkl
c
c          Dispersion term computation
c
               clp1 = lpcharge(ilp)
               clp2 = lpcharge(jlp)
               disp6 = -0.5*clp1*0.5*clp2*c6disp*dmp6/mkl6
               disp8 = -0.5*clp1*0.5*clp2*c8disp*dmp8/mkl8
               disp10 = -0.5*clp1*0.5*clp2*c10disp*dmp10/mkl10
c
c          Exchange dispersion                  
c
               bkl = gorbrep(k)*gorbrep(l)
               ck = rpole(1,k) 
               cl = rpole(1,l) 
               gdk = 1 - ck/valemtp(atomic(k))
               gdl = 1 - cl/valemtp(atomic(l))
               tkl = bkl*gdk*gdl
               xdisp3 = 0.5*clp1*0.5*clp2*tkl*cxdlp*exp(-axdlp*mkl) 
               e = (disp6 + disp8 + disp10)*colp*
     $           facdispij*discof
               e = e + xdisp3*facdispij
c
c     use energy switching if near the cutoff distance
c
               if (rij2 .gt. cut2) then
                  rij3 = rij2 * rij
                  rij4 = rij2 * rij2
                  rij5 = rij2 * rij3
                  taper = c5*rij5 + c4*rij4 + c3*rij3
     &                       + c2*rij2 + c1*rij + c0
                  e = e * taper
                end if
               exdisp = exdisp + e
            end if  
         end do
      end do
      einter = einter + exdisp
      return
      end
c
c     "dispersion3a" calculates the dispersion energy with a double loop (and replicas for large cutoff)
c     
      subroutine edispersion0a
      use analyz
      use atmtyp
      use atmlst
      use atoms
      use bond
      use bound
      use boxes
      use cell
      use charge
      use chargetransfer
      use cutoff
      use dispersion
      use domdec
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mpole
      use neigh
      use shunt
      use usage
      use potent
      use katoms
      use repulsion
      implicit none
c
      integer i,j,k,iglob,jglob,ilp,jlp
      integer m
      integer ilpnl,jj,jnl
      real*8 xi,yi,zi,xj,yj,zj,rij2,rij
      real*8 rij3,rij4,rij5
      real*8 xij,yij,zij
      real*8 rdispi,rdispj,mij,mij6,mij8,mij10
      real*8 bij,rdpi,rdpj,dij
      real*8 rij6,rij8,rij10
      real*8 dmp6,dmp8,dmp10,disp6,disp8,disp10
      real*8 ci,cj,gdi,gdj,tij,edisp
      real*8 xlp,ylp,zlp,xk,yk,zk
      real*8 xklp,yklp,zklp,rklp2,rklp
      real*8 xdisp,rsij
      real*8 expo6
      real*8 edisp1,xdisp1,exdisp1
      real*8 edisp2,xdisp2,exdisp2
      real*8 edisp3,xdisp3,exdisp3
      real*8 xkj,ykj,zkj,rkj2,rkj 
      real*8 rkj3,rkj4,rkj5
      real*8 bkj,bkjx,rdpk,dkj,rskj
      real*8 mkj,mkj6,mkj8,mkj10
      real*8 clp,gdk,tkj,clp1,clp2
      real*8 xl,yl,zl
      real*8 xkl,ykl,zkl,rkl2,rkl
      real*8 bkl,rdpl,dkl,rskl
      real*8 mkl,mkl6,mkl8,mkl10
      real*8 ck,cl,gdl,tkl
      real*8 rlp
      real*8 taper,e
      real*8 xmove,ymove,zmove
      integer l
      logical docompute
      character*6 mode
c
      call rotlp
c
      disp6 = 0.0d0
      disp8 = 0.0d0
      disp10 = 0.0d0
      exdisp = 0.0d0
      xdisp = 0.0d0
      edisp = 0.0d0
      edisp1 = 0.0d0
      edisp2 = 0.0d0
      edisp3 = 0.0d0
      xdisp1 = 0.0d0
      xdisp2 = 0.0d0
      xdisp3 = 0.0d0
      exdisp1 = 0.0d0
      exdisp2 = 0.0d0
      exdisp3 = 0.0d0

c
c     set the coefficients for the switching function
c
      mode = 'DISP'
      call switch (mode)
c
c     atomic dispersion computation
c
      do i = 1, nloc
         iglob = glob(i)
         xi = x(iglob) 
         yi = y(iglob) 
         zi = z(iglob) 
         do j = 1, n
            jglob = j
            if (jglob.le.iglob) cycle
            if (molcule(iglob) .ne. molcule(jglob)) then
               xj = x(jglob) 
               yj = y(jglob) 
               zj = z(jglob) 
               xij = xj - xi 
               yij = yj - yi 
               zij = zj - zi 
c               if (use_bounds)  call image (xij,yij,zij)
               rij2 = xij*xij + yij*yij + zij*zij
c
               if (rij2.gt.(dispcut*dispcut)) cycle
c
               rij = sqrt(rij2) 
               bij =  gorbrep(iglob)*gorbrep(jglob)
               rdpi = vdwdisp1(iglob)
               rdpj = vdwdisp1(jglob)
               dij = ((rdpi+rdpj)*bdmp/rij)-1
               rsij = (rdpi + rdpj)*bdmp
               expo6 = admp6*dij
c              si dump pas assez efficace , > 1
               if (rij.ge.rsij) then
                 dmp6 = bij
                 dmp8 = bij
                 dmp10 = bij
               else
                 dmp6 = bij*exp(-expo6)
                 dmp8 = bij*exp(-admp8*dij)
                 dmp10 = bij*exp(-admp10*dij)
               end if
               mij = rij/(2*sqrt(rdpi*rdpj))
               mij6 = mij*mij*mij*mij*mij*mij 
               mij8 = mij6*mij*mij
               mij10 = mij8*mij*mij
               disp6 = c6disp*dmp6/mij6
               disp8 = c8disp*dmp8/mij8
               disp10 = c10disp*dmp10/mij10
c              
c                 exchange disp              
c
               ci = rpole(1,iglob) 
               cj = rpole(1,jglob)
               gdi = 1 - ci/valemtp(atomic(iglob)) 
               gdj = 1 - cj/valemtp(atomic(jglob)) 
               tij = bij*gdi*gdj
               xdisp1 = tij*cxd*exp(-axd*mij)
               e = (disp6 + disp8 + disp10)*scdp*
     $         facdispij*discof
               e = e + xdisp1*facdispij
c
c     use energy switching if near the cutoff distance
c
               if (rij2 .gt. cut2) then
                  rij3 = rij2 * rij
                  rij4 = rij2 * rij2
                  rij5 = rij2 * rij3
                  taper = c5*rij5 + c4*rij4 + c3*rij3
     &                       + c2*rij2 + c1*rij + c0
                  e = e * taper
               end if
               exdisp = exdisp + e
            end if
         end do
      end do  
c
c     lp-atom dispersion
c
      do i = 1, nlploc
         ilp = lpglob(i)
         k = lpatom(ilp)
         xlp = rlonepair(1,ilp)
         ylp = rlonepair(2,ilp)
         zlp = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, n
            if (molcule(k) .ne. molcule(j)) then
               disp6 =0.0d0
               disp8 = 0.0d0
               disp10 = 0.0d0
               dmp6 =0.0d0
               dmp8 = 0.0d0
               dmp10 = 0.0d0
c               
               xj = x(j) 
               yj = y(j) 
               zj = z(j) 
               xkj = xj - xlp
               ykj = yj - ylp
               zkj = zj - zlp
c               if (use_bounds)  call image (xkj,ykj,zkj)
               rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
c
               if (rkj2.gt.(dispcut*dispcut)) cycle
c
               rkj = sqrt(rkj2)
               bkj = gorbrep(j)
               bkjx = gorbrep(k)*gorbrep(j)
c
c   vdw radius is the vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp2(k)! + dincr_lprep(ilp)
               rdpj = vdwdisp1(j)
               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
               rskj = (rdpk + rdpj)*bdmp
               expo6 = -admp6*dkj
               if (rkj > rskj) expo6=0.0d0
               dmp6 = bkj*exp(expo6)
               dmp8 = bkj*exp(-admp8*dkj)
               dmp10 = bkj*exp(-admp10*dkj)
               mkj = rkj/(2*sqrt(rdpk*rdpj))
               mkj6 = mkj*mkj*mkj*mkj*mkj*mkj 
               mkj8 = mkj6*mkj*mkj
               mkj10 = mkj8*mkj*mkj
               clp = lpcharge(ilp)!2.0d0
               disp6 = -0.5*clp*colpa*c6disp*dmp6/mkj6
               disp8 = -0.5*clp*colpa*c8disp*dmp8/mkj8
               disp10 = -0.5*clp*colpa*c10disp*dmp10/mkj10
c
c              exchange disp
c
               cj = rpole(1,j)
               ck = rpole(1,k)
               gdk = 1 - ck/valemtp(atomic(k)) 
               gdj = 1 - cj/valemtp(atomic(j)) 
               tkj = bkjx*gdk*gdj
               xdisp2 = 0.5*clp*tkj*cxdla*exp(-axdla*mkj)
               e = (disp6 + disp8 + disp10)*facdispij*
     $            discof
               e = e + xdisp2*facdispij
c
c     use energy switching if near the cutoff distance
c
               if (rkj2 .gt. cut2) then
                  rkj3 = rkj2 * rkj
                  rkj4 = rkj2 * rkj2
                  rkj5 = rkj2 * rkj3
                  taper = c5*rkj5 + c4*rkj4 + c3*rkj3
     &                       + c2*rkj2 + c1*rkj + c0
                  e = e * taper
               end if
               exdisp = exdisp + e
            end if   
         end do 
      end do  
c
c      lone pair / lone pair dispersion term 
c
      do i = 1, nlploc
         ilp = lpglob(i)
         k = lpatom(ilp)    
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, nlp
            jlp = j
            if (jlp.le.ilp) cycle
            l = lpatom(jlp)
            if (molcule(k) .ne. molcule(l)) then
               disp6 =0.0d0
               disp8 = 0.0d0
               disp10 = 0.0d0
               xj = rlonepair(1,jlp)
               yj = rlonepair(2,jlp)
               zj = rlonepair(3,jlp)
               xl = x(l) 
               yl = y(l) 
               zl = z(l) 
               xij = xj - xi 
               yij = yj - yi 
               zij = zj - zi 
c               if (use_bounds)  call image (xij,yij,zij)
               rij2 = xij*xij + yij*yij + zij*zij
c
               if (rij2.gt.(dispcut*dispcut)) cycle
c
               rij = sqrt(rij2)
               bkl = 1  
c
c    vdw radius is vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp3(k)! + dincr_lprep(ilp) 
               rdpl = vdwdisp3(l)! + dincr_lprep(jlp) 
c
               dkl = ((rdpk+rdpl)*bdmp/rij)-1
               rskl = (rdpk + rdpl)*bdmp
               expo6 = -admp6*dkl
               if (rij > rskl) expo6 = 0.0d0 
               dmp6 = bkl*exp(expo6)
               dmp8 = bkl*exp(-admp8*dkl)
               dmp10 = bkl*exp(-admp10*dkl)
               mkl = rij/(2*sqrt(rdpk*rdpl))
               mkl6 = mkl*mkl*mkl*mkl*mkl*mkl 
               mkl8 = mkl6*mkl*mkl
               mkl10 = mkl8*mkl*mkl
c
c          Dispersion term computation
c
               clp1 = lpcharge(ilp)
               clp2 = lpcharge(jlp)
               disp6 = -0.5*clp1*0.5*clp2*c6disp*dmp6/mkl6
               disp8 = -0.5*clp1*0.5*clp2*c8disp*dmp8/mkl8
               disp10 = -0.5*clp1*0.5*clp2*c10disp*dmp10/mkl10
c
c          Exchange dispersion                  
c
               bkl = gorbrep(k)*gorbrep(l)
               ck = rpole(1,k) 
               cl = rpole(1,l) 
               gdk = 1 - ck/valemtp(atomic(k))
               gdl = 1 - cl/valemtp(atomic(l))
               tkl = bkl*gdk*gdl
               xdisp3 = 0.5*clp1*0.5*clp2*tkl*cxdlp*exp(-axdlp*mkl) 
               e = (disp6 + disp8 + disp10)*colp*
     $           facdispij*discof
               e = e + xdisp3*facdispij
c
c     use energy switching if near the cutoff distance
c
               if (rij2 .gt. cut2) then
                  rij3 = rij2 * rij
                  rij4 = rij2 * rij2
                  rij5 = rij2 * rij3
                  taper = c5*rij5 + c4*rij4 + c3*rij3
     &                       + c2*rij2 + c1*rij + c0
                  e = e * taper
               end if
               exdisp = exdisp + e
            end if  
         end do
      end do
      einter = einter + exdisp
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
c     atomic dispersion computation
c
      do i = 1, nloc
         iglob = ineignl(i)
         xi = x(iglob) 
         yi = y(iglob) 
         zi = z(iglob) 
c
         do j = 1, n
            jglob = j
c            if (jglob.le.iglob) cycle
c
               do k = 2, ncell
c            if (molcule(iglob) .ne. molcule(jglob)) then
c
c     set the distance to translate along each cell axis
c
                 xmove = icell(1,k) * xbox
                 ymove = icell(2,k) * ybox
                 zmove = icell(3,k) * zbox
                 xj = x(jglob) 
                 yj = y(jglob) 
                 zj = z(jglob) 
                 xij = xj - xi 
                 yij = yj - yi 
                 zij = zj - zi 
                 xij = xij + xmove
                 yij = yij + ymove
                 zij = zij + zmove
c                 if (use_bounds)  call image (xij,yij,zij)
                 rij2 = xij*xij + yij*yij + zij*zij
c
                 if (rij2.gt.(dispcut*dispcut)) cycle
c
                 rij = sqrt(rij2) 
                 bij =  gorbrep(iglob)*gorbrep(jglob)
                 rdpi = vdwdisp1(iglob)
                 rdpj = vdwdisp1(jglob)
                 dij = ((rdpi+rdpj)*bdmp/rij)-1
                 rsij = (rdpi + rdpj)*bdmp
                 expo6 = admp6*dij
c                si dump pas assez efficace , > 1
                 if (rij.ge.rsij) then
                   dmp6 = bij
                   dmp8 = bij
                   dmp10 = bij
                 else
                   dmp6 = bij*exp(-expo6)
                   dmp8 = bij*exp(-admp8*dij)
                   dmp10 = bij*exp(-admp10*dij)
                 end if
                 mij = rij/(2*sqrt(rdpi*rdpj))
                 mij6 = mij*mij*mij*mij*mij*mij 
                 mij8 = mij6*mij*mij
                 mij10 = mij8*mij*mij
                 disp6 = c6disp*dmp6/mij6
                 disp8 = c8disp*dmp8/mij8
                 disp10 = c10disp*dmp10/mij10
c                
c                   exchange disp              
c
                 ci = rpole(1,iglob) 
                 cj = rpole(1,jglob)
                 gdi = 1 - ci/valemtp(atomic(iglob)) 
                 gdj = 1 - cj/valemtp(atomic(jglob)) 
                 tij = bij*gdi*gdj
                 xdisp1 = tij*cxd*exp(-axd*mij)
                 e = (disp6 + disp8 + disp10)*scdp*
     $           facdispij*discof
                 e = e + xdisp1*facdispij
c
c     use energy switching if near the cutoff distance
c
                 if (rij2 .gt. cut2) then
                    rij3 = rij2 * rij
                    rij4 = rij2 * rij2
                    rij5 = rij2 * rij3
                    taper = c5*rij5 + c4*rij4 + c3*rij3
     &                         + c2*rij2 + c1*rij + c0
                    e = e * taper
                 end if
                 exdisp = exdisp + e
               end do
c            end if
         end do
      end do  
c
c     lp-atom dispersion
c
      do i = 1, nlploc
         ilp = lpglob(i)
         k = lpatom(ilp)
         xlp = rlonepair(1,ilp)
         ylp = rlonepair(2,ilp)
         zlp = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, n
               do l = 2, ncell
c            if (molcule(iglob) .ne. molcule(jglob)) then
c
c     set the distance to translate along each cell axis
c
                 xmove = icell(1,l) * xbox
                 ymove = icell(2,l) * ybox
                 zmove = icell(3,l) * zbox
                 disp6 =0.0d0
                 disp8 = 0.0d0
                 disp10 = 0.0d0
                 dmp6 =0.0d0
                 dmp8 = 0.0d0
                 dmp10 = 0.0d0
c                 
                 xj = x(j) 
                 yj = y(j) 
                 zj = z(j) 
                 xkj = xj - xlp
                 ykj = yj - ylp
                 zkj = zj - zlp
                 xkj = xkj + xmove
                 ykj = ykj + ymove
                 zkj = zkj + zmove
c                 if (use_bounds)  call image (xkj,ykj,zkj)
                 rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
c
                 if (rkj2.gt.(dispcut*dispcut)) cycle
c
                 rkj = sqrt(rkj2)
                 bkj = gorbrep(j)
                 bkjx = gorbrep(k)*gorbrep(j)
c
c   vdw radius is the vdw of the carrier! plus increment of the lp
c
                 rdpk = vdwdisp2(k)! + dincr_lprep(ilp)
                 rdpj = vdwdisp1(j)
                 dkj = ((rdpk+rdpj)*bdmp/rkj)-1
                 rskj = (rdpk + rdpj)*bdmp
                 expo6 = -admp6*dkj
                 if (rkj > rskj) expo6=0.0d0
                 dmp6 = bkj*exp(expo6)
                 dmp8 = bkj*exp(-admp8*dkj)
                 dmp10 = bkj*exp(-admp10*dkj)
                 mkj = rkj/(2*sqrt(rdpk*rdpj))
                 mkj6 = mkj*mkj*mkj*mkj*mkj*mkj 
                 mkj8 = mkj6*mkj*mkj
                 mkj10 = mkj8*mkj*mkj
                 clp = lpcharge(ilp)!2.0d0
                 disp6 = -0.5*clp*colpa*c6disp*dmp6/mkj6
                 disp8 = -0.5*clp*colpa*c8disp*dmp8/mkj8
                 disp10 = -0.5*clp*colpa*c10disp*dmp10/mkj10
c
c                exchange disp
c
                 cj = rpole(1,j)
                 ck = rpole(1,k)
                 gdk = 1 - ck/valemtp(atomic(k)) 
                 gdj = 1 - cj/valemtp(atomic(j)) 
                 tkj = bkjx*gdk*gdj
                 xdisp2 = 0.5*clp*tkj*cxdla*exp(-axdla*mkj)
                 e = (disp6 + disp8 + disp10)*facdispij*
     $              discof
                 e = e + xdisp2*facdispij
c
c     use energy switching if near the cutoff distance
c
                 if (rkj2 .gt. cut2) then
                    rkj3 = rkj2 * rkj
                    rkj4 = rkj2 * rkj2
                    rkj5 = rkj2 * rkj3
                    taper = c5*rkj5 + c4*rkj4 + c3*rkj3
     &                         + c2*rkj2 + c1*rkj + c0
                    e = e * taper
                 end if
                 exdisp = exdisp + e
               end do
c            end if   
         end do 
      end do  
c
c      lone pair / lone pair dispersion term 
c
      do i = 1, nlploc
         ilp = lpglob(i)
         k = lpatom(ilp)    
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, nlp
            jlp = j
c            if (jlp.le.ilp) cycle
            l = lpatom(jlp)
c            if (molcule(k) .ne. molcule(l)) then
               do m = 2, ncell
c
c     set the distance to translate along each cell axis
c
                 xmove = icell(1,m) * xbox
                 ymove = icell(2,m) * ybox
                 zmove = icell(3,m) * zbox

                 disp6 =0.0d0
                 disp8 = 0.0d0
                 disp10 = 0.0d0
                 xj = rlonepair(1,jlp)
                 yj = rlonepair(2,jlp)
                 zj = rlonepair(3,jlp)
                 xl = x(l) 
                 yl = y(l) 
                 zl = z(l) 
                 xij = xj - xi 
                 yij = yj - yi 
                 zij = zj - zi 
                 xij = xij + xmove
                 yij = yij + ymove
                 zij = zij + zmove
c                 if (use_bounds)  call image (xij,yij,zij)
                 rij2 = xij*xij + yij*yij + zij*zij
c
                 if (rij2.gt.(dispcut*dispcut)) cycle
c
                 rij = sqrt(rij2)
                 bkl = 1  
c
c    vdw radius is vdw of the carrier! plus increment of the lp
c
                 rdpk = vdwdisp3(k)! + dincr_lprep(ilp) 
                 rdpl = vdwdisp3(l)! + dincr_lprep(jlp) 
c
                 dkl = ((rdpk+rdpl)*bdmp/rij)-1
                 rskl = (rdpk + rdpl)*bdmp
                 expo6 = -admp6*dkl
                 if (rij > rskl) expo6 = 0.0d0 
                 dmp6 = bkl*exp(expo6)
                 dmp8 = bkl*exp(-admp8*dkl)
                 dmp10 = bkl*exp(-admp10*dkl)
                 mkl = rij/(2*sqrt(rdpk*rdpl))
                 mkl6 = mkl*mkl*mkl*mkl*mkl*mkl 
                 mkl8 = mkl6*mkl*mkl
                 mkl10 = mkl8*mkl*mkl
c
c          Dispersion term computation
c
                 clp1 = lpcharge(ilp)
                 clp2 = lpcharge(jlp)
                 disp6 = -0.5*clp1*0.5*clp2*c6disp*dmp6/mkl6
                 disp8 = -0.5*clp1*0.5*clp2*c8disp*dmp8/mkl8
                 disp10 = -0.5*clp1*0.5*clp2*c10disp*dmp10/mkl10
c
c          Exchange dispersion                  
c
                 bkl = gorbrep(k)*gorbrep(l)
                 ck = rpole(1,k) 
                 cl = rpole(1,l) 
                 gdk = 1 - ck/valemtp(atomic(k))
                 gdl = 1 - cl/valemtp(atomic(l))
                 tkl = bkl*gdk*gdl
                 xdisp3 = 0.5*clp1*0.5*clp2*tkl*cxdlp*exp(-axdlp*mkl) 
                 e = (disp6 + disp8 + disp10)*colp*
     $             facdispij*discof
                 e = e + xdisp3*facdispij
c
c     use energy switching if near the cutoff distance
c
                 if (rij2 .gt. cut2) then
                    rij3 = rij2 * rij
                    rij4 = rij2 * rij2
                    rij5 = rij2 * rij3
                    taper = c5*rij5 + c4*rij4 + c3*rij3
     &                         + c2*rij2 + c1*rij + c0
                    e = e * taper
                 end if
                 exdisp = exdisp + e
               end do
c            end if  
         end do
      end do
      return
      end
