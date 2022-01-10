c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine erepulsion  --  repulsion energy                    ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "erepulsion" calculates the repulsion energy
c
c
      subroutine erepulsion
      implicit none
      call erepulsion0b
      return
      end
c
c     "erepulsion0b" calculates the repulsion energy
c     
      subroutine erepulsion0b
      use analyz
      use atmtyp
      use atmlst
      use atoms
      use bound
      use cell
      use charge
      use cutoff
      use domdec
      use energi
      use inform
      use inter
      use iounit
      use molcul
      use mpole
      use neigh
      use usage
      use bond
      use repulsion !! to complete : kac=tabulé 
      use chargetransfer
      use potent
      use katoms
      use kct
      implicit none
      integer i,j,ia,ib,jc,jd,inl,jj,ilpnl
      integer k,p
      real*8 xa,xb,xc,xd
      real*8 ya,yb,yc,yd
      real*8 za,zb,zc,zd
      real*8 xab,xac,xbd,xbc 
      real*8 yab,yac,ybd,ybc 
      real*8 zab,zac,zbd,zbc 
      real*8 rab2,rac2,rbd2,rbc2
      real*8 rab,rac,rbd,rbc 
      real*8 xcd,ycd,zcd
      real*8 xad,yad,zad
      real*8 rda2,rcd2,rda,rcd
      real*8 kac,kad,kbc,kbd
      real*8 bmac2,bmad2,bmbc2,bmbd2
      real*8 bmac,bmad,bmbc,bmbd
      real*8 pac,pda,pbc,pbd
      real*8 sasc,sasd,sbsc,sbsd
      real*8 sasc2,sasd2,sbsc2,sbsd2
      real*8 cos11,cos12,cos13,cos14
      real*8 cos15,cos16,cos17,cos18
      real*8 csa,csb,csc,csd
      real*8 cpa,cpb,cpc,cpd
      real*8 olap1ac,olap1ad,olap1bc,olap1bd
      real*8 olap1ac2,olap1ad2,olap1bc2,olap1bd2
      real*8 onab,oncd
      real*8 bxab,bxcd,byab,bycd,bzab,bzcd
      real*8 xdist1,ydist1,zdist1,dist12
      real*8 dist1,e11,e12
      real*8 ss11,ss12,ss13,ss14
      real*8 sp11,sp12,sp13,sp14
      real*8 ps11,ps12,ps13,ps14
      real*8 pp11,pp12,pp13,pp14
      real*8 olap1ss,olap1ps,olap1sp,olap1pp
      real*8 olap1,crep1,erep1,olap12
      real*8 spac,spad,spbc,spbd
      real*8 psac,psad,psbc,psbd
      real*8 erep2,erep3,erep22,erep32
      real*8 xlp,ylp,zlp,xk,yk,zk
      real*8 xalp,yalp,zalp,ralp2,ralp
      real*8 xklp,yklp,zklp,rklp2,rklp
      real*8 kalp,kblp,kak,kbk
      real*8 bmak2,bmak,pak,sask,sask2
      real*8 bmbk2,bmbk,pbk,sbsk,sbsk2
      real*8 xblp,yblp,zblp,rblp2,rblp
      real*8 xbk,ybk,zbk,rbk2,rbk
      real*8 cos21,cos22,sp21,sp22
      real*8 cos23,cos24,pp21,pp22
      real*8 onlp,dist22,dist2
      real*8 olap2ak,olap2ak2,olap2bk,olap2bk2,olap2,olap22
      real*8 e21,e22
      real*8 cpk,spak,spbk
      real*8 xak,yak,zak,rak2,rak
      real*8 xdist2,ydist2,zdist2
      real*8 cia,cib,cjc,cjd
      real*8 xi,yi,zi,xp,yp,zp
      real*8 xip,yip,zip,rip2,rip
      real*8 cilp,cjlp,cjk,cip
      real*8 xj,yj,zj
      real*8 xjk,yjk,zjk,rjk2,rjk
      real*8 xij,yij,zij,rij2,rij
      real*8 kij,bmij2,bmij,pij,sisj,sisj2
      real*8 cpi,cpj,cos31,cos32
      real*8 pp31,oni,onj
      real*8 bxip,byip,bzip,bxjk,byjk,bzjk
      real*8 xdist3,ydist3,zdist3,dist32,dist3
      real*8 e31,e32,olap3,olap32
      real*8 psak,psbk,ps21,ps22,csk
      real*8 ss21,ss22
      real*8 erepc1
      real*8 csi,csj,spij,psij,ss31,ps31,sp31
      real*8 xpk,ypk,zpk,rpk2,rpk
      real*8 vdwreplp,vdwreplpi,vdwreplpj
      real*8 erep12
      logical docompute
      integer ibond,jlp,ilp,jlpnl,ii,jnl
c
c
      erep1 = 0.0d0
      erep2 = 0.0d0
      erep3 = 0.0d0
      erep12 = 0.0d0
      erep22 = 0.0d0
      erep32 = 0.0d0
c
c     zero out the repulsion  energy and partitioning terms
c
      erep = 0.0d0
      call rotlp
c
c     calculate bond-bond repulsion term
c
      do inl = 1, nbondlocnl
         i = bondglobnl(inl)
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         cia = rpole(1,ia)
         cib = rpole(1,ib)
         xa = x(ia)  
         ya = y(ia)
         za = z(ia)
         xb = x(ib)
         yb = y(ib)
         zb = z(ib)
         xab = xb - xa
         yab = yb - ya
         zab = zb - za
         rab2 = xab*xab + yab*yab + zab*zab
         rab = sqrt(rab2)
         do jj = 1, nbondlst(inl) 
            j = bondlst(jj,inl)
            jc = ibnd(1,j)                     
            jd = ibnd(2,j)
            if (molcule(ia) .ne. molcule(jc)) then
            xc = x(jc)  
            yc = y(jc)
            zc = z(jc)
            xd = x(jd)
            yd = y(jd)
            zd = z(jd)
            xcd = xd - xc
            ycd = yd - yc
            zcd = zd - zc
            rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
            rcd = sqrt(rcd2)
            xad = xd - xa
            yad = yd - ya
            zad = zd - za
            if (use_bounds)  call image (xad,yad,zad)
            rda2 = xad*xad + yad*yad + zad*zad
c
            if (rda2.gt.(repcut*repcut)) cycle
c
            rda = sqrt(rda2)
            xbd = xd - xb
            ybd = yd - yb
            zbd = zd - zb
            if (use_bounds)  call image (xbd,ybd,zbd)
            rbd2 = xbd*xbd + ybd*ybd + zbd*zbd         
            rbd = sqrt(rbd2)
            xac = xc - xa
            yac = yc - ya
            zac = zc - za
            if (use_bounds)  call image (xac,yac,zac)
            rac2 = xac*xac + yac*yac + zac*zac         
            rac = sqrt(rac2)
            xbc = xc - xb
            ybc = yc - yb
            zbc = zc - zb
            if (use_bounds)  call image (xbc,ybc,zbc)
            rbc2 = xbc*xbc + ybc*ybc + zbc*zbc         
            rbc = sqrt(rbc2)
            cjc = rpole(1,jc)
            cjd = rpole(1,jd)
c                                                                                   
c     Get the proportionnality factor between two atoms for overlap AO s 
c
            kac = gorbrep(ia)*gorbrep(jc)
            kad = gorbrep(ia)*gorbrep(jd)
            kbc = gorbrep(ib)*gorbrep(jc)
            kbd = gorbrep(ib)*gorbrep(jd)
c
c      orbital s contribution computation 
c                 AO atom A on C 
c
            bmac2 = kac*((1-cia/valemtp(atomic(ia)))*         
     &                    (1-cjc/valemtp(atomic(jc))))
            bmac = sqrt(bmac2)
            pac = rac/(4*sqrt(vdwrep(ia)*vdwrep(jc)))
            sasc = bmac*exp(-alpha*pac)
            sasc2 = bmac*exp(-alpha2*pac)
c                     
c           AO atom A on D 
c
            bmad2 = kad*((1-cia/valemtp(atomic(ia)))*
     &              (1-cjd/valemtp(atomic(jd))))
            bmad = sqrt(bmad2)
            pda = rda/(4*sqrt(vdwrep(ia)*vdwrep(jd)))
            sasd = bmad*exp(-alpha*pda)
            sasd2 = bmad*exp(-alpha2*pda)
            
c           AO atom B on C                                        
c                                                             
            bmbc2 = kbc*((1-cib/valemtp(atomic(ib)))*       
     &              (1-cjc/valemtp(atomic(jc))))
            bmbc = sqrt(bmbc2)
            pbc = rbc/(4*sqrt(vdwrep(ib)*vdwrep(jc)))
            sbsc = bmbc*exp(-alpha*pbc)
            sbsc2 = bmbc*exp(-alpha2*pbc)
c
c          AO atom B on D                 
c
            bmbd2 = kbd*((1-cib/valemtp(atomic(ib)))*
     &               (1-cjd/valemtp(atomic(jd))))
            bmbd = sqrt(bmbd2)
            pbd = rbd/(4*sqrt(vdwrep(ib)*vdwrep(jd)))
            sbsd = bmbd*exp(-alpha*pbd)
            sbsd2 = bmbd*exp(-alpha2*pbd)
c
            csa = chyb(1,ia)
            csb = chyb(1,ib)
            csc = chyb(1,jc)
            csd = chyb(1,jd)
            cpa = chyb(2,ia)
            cpb = chyb(2,ib)
            cpc = chyb(2,jc)
            cpd = chyb(2,jd)
c
c           exponential prefactor ss -> sp
c
            psac = forb(atomic(ia),atomic(jc))
            psad = forb(atomic(ia),atomic(jd))
            psbc = forb(atomic(ib),atomic(jc))
            psbd = forb(atomic(ib),atomic(jd))
            spac = forb(atomic(jc),atomic(ia))
            spbc = forb(atomic(jc),atomic(ib))
            spad = forb(atomic(jd),atomic(ia)) 
            spbd = forb(atomic(jd),atomic(ib))
c
c           overlap bond-bond ss
c
            ss11 = csa*csc
            ss12 = csa*csd
            ss13 = csb*csc
            ss14 = csb*csd
c
c           orbital overlap p(AB) s(CD)
c
            cos11 = (xab*xac + yab*yac + zab*zac)/(rab*rac)
            cos12 = (xab*xad + yab*yad + zab*zad)/(rab*rda)
            cos13 = (-xab*xbc - yab*ybc - zab*zbc)/(rab*rbc)
            cos14 = (-xab*xbd - yab*ybd - zab*zbd)/(rab*rbd)
            cos15 = (-xcd*xac - ycd*yac - zcd*zac)/(rcd*rac)
            cos16 = (-xcd*xbc - ycd*ybc - zcd*zbc)/(rcd*rbc)
            cos17 = (xcd*xad + ycd*yad + zcd*zad)/(rcd*rda) 
            cos18 = (xcd*xbd + ycd*ybd + zcd*zbd)/(rcd*rbd)
            ps11 = cpa*csc*psac*cos11
            ps12 = cpa*csd*psad*cos12
            ps13 = cpb*csc*psbc*cos13
            ps14 = cpb*csd*psbd*cos14
c
c           orbital overlap s(AB) p(CD)
c
            sp11 = csa*cpc*spac*cos15
            sp12 = csa*cpd*spad*cos17
            sp13 = csb*cpc*spbc*cos16
            sp14 = csb*cpd*spbd*cos18
c
c           orbital overlap pp
c
            pp11 = cpa*cpc*2*cos11*cos15
            pp12 = cpa*cpd*2*cos12*cos17
            pp13 = cpb*cpc*2*cos13*cos16
            pp14 = cpb*cpd*2*cos14*cos18
c
c           final summation of the overlap terms
c
            olap1ac = (ss11 + ps11 + sp11 + pp11)*sasc
            olap1ad = (ss12 + ps12 + sp12 + pp12)*sasd
            olap1bc = (ss13 + ps13 + sp13 + pp13)*sbsc
            olap1bd = (ss14 + ps14 + sp14 + pp14)*sbsd
            olap1 = olap1ac + olap1ad + olap1bc + olap1bd
            olap1ac2 = (ss11 + ps11 + sp11 + pp11)*sasc2
            olap1ad2 = (ss12 + ps12 + sp12 + pp12)*sasd2
            olap1bc2 = (ss13 + ps13 + sp13 + pp13)*sbsc2
            olap1bd2 = (ss14 + ps14 + sp14 + pp14)*sbsd2
            olap12 = olap1ac2 + olap1ad2 + olap1bc2 + olap1bd2
c
c           bond occupation number 
c
            onab = 2.0d0
            oncd = 2.0d0
            bxab = (xa + xb)/2                           
            byab = (ya + yb)/2
            bzab = (za + zb)/2
            bxcd = (xc + xd)/2
            bycd = (yc + yd)/2
            bzcd = (zc + zd)/2
            xdist1 = bxab - bxcd
            ydist1 = byab - bycd 
            zdist1 = bzab - bzcd
            if (use_bounds)  call image (xdist1,ydist1,zdist1)
            dist12 = xdist1*xdist1 + ydist1*ydist1
     &               + zdist1*zdist1                    
            dist1 = sqrt(dist12) 
            e11 = onab*olap1*olap1*oncd/dist1
            e12 = onab*olap12*olap12*oncd/dist12
            erep = erep + e11*cvrep11 + e12*cvrep12
c
            end if
         end do
      end do   
c           
c     Compute Bond - Lone Pair repulsion
c
      do inl = 1, nbondlocnl
         i = bondglobnl(inl)
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         cia = rpole(1,ia)
         cib = rpole(1,ib)
         xa = x(ia)  
         ya = y(ia)
         za = z(ia)
         xb = x(ib)
         yb = y(ib)
         zb = z(ib)
         xab = xb - xa
         yab = yb - ya
         zab = zb - za
         rab2 = xab*xab + yab*yab + zab*zab
         rab = sqrt(rab2)
         do ii = 1, nbondlplst(inl)
            jlp = bondlplst(ii,inl)
            k = lpatom(jlp)       
            if (molcule(ia) .ne. molcule(k)) then  
             xlp = rlonepair(1,jlp)
             ylp = rlonepair(2,jlp)
             zlp = rlonepair(3,jlp)
             xk = x(k) 
             yk = y(k) 
             zk = z(k) 
             xalp = xlp - xa 
             yalp = ylp - ya 
             zalp = zlp - za 
             if (use_bounds)  call image (xalp,yalp,zalp)
             ralp2 = (xalp*xalp + yalp*yalp + zalp*zalp) 
c
             if (ralp2.gt.(repcut*repcut)) cycle
c
             ralp = sqrt(ralp2)
             xblp = xlp - xb 
             yblp = ylp - yb 
             zblp = zlp - zb 
             if (use_bounds)  call image (xblp,yblp,zblp)
             rblp2 = (xblp*xblp + yblp*yblp + zblp*zblp) 
             rblp = sqrt(rblp2)
             xklp = xlp - xk 
             yklp = ylp - yk 
             zklp = zlp - zk 
             rklp2 = xklp*xklp + yklp*yklp + zklp*zklp
             rklp = sqrt(rklp2)
             xak = xk - xa
             yak = yk - ya
             zak = zk - za
             if (use_bounds)  call image (xak,yak,zak)
             rak2 = xak*xak + yak*yak + zak*zak
             rak = sqrt(rak2)
             xbk = xk - xb
             ybk = yk - yb
             zbk = zk - zb
             if (use_bounds)  call image (xbk,ybk,zbk)
             rbk2 = xbk*xbk + ybk*ybk + zbk*zbk
             rbk = sqrt(rbk2)
             cjk = rpole(1,k)
             cjlp = lpcharge(jlp)
             kak = gorbrep(ia)*gorbrep(k)
             kbk = gorbrep(ib)*gorbrep(k)
            
             bmak2 = kak*((1-cia/valemtp(atomic(ia)))*
     &                (1-cjk/valemtp(atomic(k))))
             bmak = sqrt(bmak2)
c
c    vdw radius is vdw of the carrier plus increment of the lp
c
             vdwreplp = vdwrep(k) + dincr_lprep(jlp)
             pak = rak/(4*(sqrt(vdwrep(ia)*vdwreplp)))
             sask = bmak*exp(-alpha*pak)            
             sask2 = bmak*exp(-alpha2*pak)    

             bmbk2 = kbk*((1-cib/valemtp(atomic(ib)))*
     &                (1-cjk/valemtp(atomic(k))))
             bmbk = sqrt(bmbk2)
             pbk = rbk/(4*(sqrt(vdwrep(ib)*vdwreplp)))
             sbsk = bmbk*exp(-alpha*pbk) 
             sbsk2 = bmbk*exp(-alpha2*pbk) 
            
             csa = chyb(1,ia)
             csb = chyb(1,ib)
             csk = chyb(1,k)
             cpa = chyb(2,ia)
             cpb = chyb(2,ib)
             cpk = chyb(2,k)
              
             spak = forb(atomic(ia),atomic(k))
             spbk = forb(atomic(ib),atomic(k))
             psak = forb(atomic(k),atomic(ia))
             psbk = forb(atomic(k),atomic(ib))
            
             ss21 = csa*csk
             ss22 = csb*csk
            
             cos21 = - (xak*xklp + yak*yklp + zak*zklp)/(rak*rklp)
             cos22 = - (xbk*xklp + ybk*yklp + zbk*zklp)/(rbk*rklp)
            
             sp21 = csa*cpk*cos21*psak
             sp22 = csb*cpk*cos22*psbk
             
             cos23 = (xab*xak + yab*yak + zab*zak)/(rab*rak)  
             cos24 = - (xab*xbk + yab*ybk + zab*zbk)/(rab*rbk)  
            
             ps21 = cpa*csk*cos23*spak
             ps22 = cpb*csk*cos24*spbk
            
             pp21 = cpa*cpk*cos21*cos23*2 
             pp22 = cpb*cpk*cos22*cos24*2
c           
             onab = 2.0
             onlp = lpcharge(jlp)!2.0
c           
             bxab = (xa + xb)/2
             byab = (ya + yb)/2
             bzab = (za + zb)/2
c
             xdist2 = bxab - xlp
             ydist2 = byab - ylp
             zdist2 = bzab - zlp
             if (use_bounds)  call image (xdist2,ydist2,zdist2)
c       
             dist22 = xdist2*xdist2 + ydist2*ydist2 + zdist2*zdist2
             dist2 = sqrt(dist22)
c       
             olap2ak = (ss21 + sp21 + ps21 + pp21)*sask
             olap2bk = (ss22 + sp22 + ps22 + pp22)*sbsk
             olap2 = olap2ak + olap2bk
             olap2ak2 = (ss21 + sp21 + ps21 + pp21)*sask2
             olap2bk2 = (ss22 + sp22 + ps22 + pp22)*sbsk2
             olap22 = olap2ak2 + olap2bk2
c             
             e21 = onab*onlp*olap2*olap2/dist2
             e22 = onab*onlp*olap22*olap22/dist22
c
             erep = erep + cvrep21*e21 + cvrep22*e22 
            end if
         end do
      end do  
c
c     lone pair / lone pair repulsion      
c
      do ilpnl = 1, nlplocnl
         ilp = lpglobnl(ilpnl)
         p = lpatom(ilp)   
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xp = x(p)
         yp = y(p)
         zp = z(p)
         xip = xp - xi
         yip = yp - yi
         zip = zp - zi
         rip2 = xip*xip + yip*yip + zip*zip
         rip = sqrt(rip2)
         cilp = lpcharge(ilp)
         cip = rpole(1,p)
         do j = 1, nlplplst(ilpnl)
            jlp = lplplst(j,ilpnl)
            k = lpatom(jlp)
            if (molcule(p) .ne. molcule(k)) then
               xj = rlonepair(1,jlp)
               yj = rlonepair(2,jlp)
               zj = rlonepair(3,jlp)
               xk = x(k)
               yk = y(k)
               zk = z(k)
               xjk = xk - xj
               yjk = yk - yj
               zjk = zk - zj
               rjk2 = xjk*xjk + yjk*yjk + zjk*zjk 
               rjk = sqrt(rjk2) 
               xij = xj - xi
               yij = yj - yi
               zij = zj - zi
               if (use_bounds)  call image (xij,yij,zij)
               rij2 = xij*xij + yij*yij + zij*zij
c
               if (rij2.gt.(repcut*repcut)) cycle
c
               rij = sqrt(rij2)
               cjk = rpole(1,k)
               xpk = xk - xp
               ypk = yk - yp
               zpk = zk - zp
               if (use_bounds)  call image (xpk,ypk,zpk)
               rpk2 = xpk*xpk + ypk*ypk + zpk*zpk
               rpk = sqrt(rpk2)

               kij = gorbrep(p)*gorbrep(k)
               cjlp = lpcharge(jlp)

               bmij2 = kij*((1-cip/valemtp(atomic(p)))*
     &                 (1-cjk/valemtp(atomic(k))))
               bmij = sqrt(bmij2)
c
c    vdw radius is vdw of the carrier plus increment of the lp
c
               vdwreplpi = vdwrep(p) + dincr_lprep(ilp)
               vdwreplpj = vdwrep(k) + dincr_lprep(jlp)
               pij = rpk/(4*(sqrt(vdwreplpi*vdwreplpj)))
               sisj = bmij*exp(-alpha*pij)
               sisj2 = bmij*exp(-alpha2*pij) 

               csi = chyb(1,p)
               csj = chyb(1,k)

               cpi = chyb(2,p)
               cpj = chyb(2,k)
               spij = forb(atomic(p),atomic(k))
               psij = forb(atomic(k),atomic(p))

               cos31 = - (xip*xpk + yip*ypk + zip*zpk)/(rpk*rip)
               cos32 = (xjk*xpk + yjk*ypk + zjk*zpk)/(rjk*rpk) 

               ss31 = csi*csj
               sp31 = csi*cpj*spij*cos32
               ps31 = csi*cpj*psij*cos31
               pp31 = cpi*cpj*2*cos31*cos32

               oni = cilp!2.0
               onj = cjlp!2.0
               olap3 = (ss31 + sp31 + ps31 + pp31)*sisj
               olap32 = (ss31 + sp31 + ps31 + pp31)*sisj2

               e31 = oni*onj*olap3*olap3/rij
               e32 = oni*onj*olap32*olap32/rij2
               erep = erep + e31*cvrep31 + e32*cvrep32 

            end if 
         end do 
      end do     
      einter = einter + erep
      return
      end
