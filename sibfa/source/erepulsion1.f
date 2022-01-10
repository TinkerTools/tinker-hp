c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine erepulsion1  --  repulsion derivatives              ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "erepulsion3" calculates the repulsion energy
c
c
      subroutine erepulsion1
      implicit none

      call erepulsion1b
      return
      end
c
c     "erepulsion1b" calculates the repulsion energy and forces
c     
      subroutine erepulsion1b
      use action
      use analyz
      use atmtyp
      use atmlst
      use atoms
      use bound
      use charge
      use chargetransfer
      use cutoff
      use deriv
      use domdec
      use energi
      use kct
      use molcul
      use mpole
      use neigh
      use usage
      use bond
      use repulsion  
      use potent
      implicit none
      real*8 d
      integer i,j,ia,ib,jc,jd
      integer ii,jj,inl,ilpnl,jlpnl
      integer ialoc,ibloc,jcloc,jdloc
      integer k,kloc,ploc,ixlploc,izlploc
      real*8 xa,xb,xc,xd
      real*8 ya,yb,yc,yd
      real*8 za,zb,zc,zd
      real*8 xab,xac,xbd,xbc,xcb,xba,xca,xdc,xda,xdb
      real*8 yab,yac,ybd,ybc,ycb,yba,yca,ydc,yda,ydb
      real*8 zab,zac,zbd,zbc,zcb,zba,zca,zdc,zda,zdb
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
      real*8 xbk,ybk,zbk,rbk2,rbk,xkb,ykb,zkb
      real*8 cos21,cos22,sp21,sp22
      real*8 cos23,cos24,pp21,pp22
      real*8 onlp,dist22,dist2
      real*8 olap2ak,olap2ak2,olap2bk,olap2bk2,olap2,olap22
      real*8 e21,e22
      real*8 cpk,spak,spbk
      real*8 xak,yak,zak,rak2,rak,xka,yka,zka
      real*8 xdist2,ydist2,zdist2
      real*8 cia,cib,cjc,cjd
      real*8 testolap
      integer p
      real*8 xi,yi,zi,xp,yp,zp
      real*8 xip,yip,zip,rip2,rip,xpi,ypi,zpi
      real*8 cilp,cjlp,cjk,cip
      real*8 xj,yj,zj
      real*8 xjk,yjk,zjk,rjk2,rjk,xkj,ykj,zkj
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
      real*8 xpk,ypk,zpk,rpk2,rpk,xkp,ykp,zkp
      real*8 erep12
      real*8 scabac,scabad,scbabc,scbabd,sccdca,sccdcb,sccdad,scdcdb
      real*8 sckaklp,sckbklp,scabak,scbabk,scpipk,sckjkp
      real*8 di(3,3),dix(3,3),diz(3,3)
      real*8 di2(3,3),dix2(3,3),diz2(3,3)
      real*8 dcoslp1,dcoslp2,dcoslp3
      real*8 dlp1,dlp2,dlp3
      real*8 vdwreplp,vdwreplpi,vdwreplpj
      real*8 xr,yr,zr
      real*8, allocatable :: dcos1(:,:),dcos2(:,:),dcos3(:,:),dcos4(:,:)
      real*8, allocatable :: dcos5(:,:),dcos6(:,:),dcos7(:,:),dcos8(:,:)
      real*8, allocatable :: dcos21(:,:),dcos22(:,:),dcos23(:,:),
     $ dcos24(:,:)
      real*8, allocatable :: dcos31(:,:),dcos32(:,:)
      real*8, allocatable :: dolap1ac(:,:),dolap1ad(:,:),dolap1bc(:,:)
      real*8, allocatable :: dolap1bd(:,:),dolaptot(:,:)
      real*8, allocatable :: dolap1ac2(:,:),dolap1ad2(:,:),
     $ dolap1bc2(:,:)
      real*8, allocatable :: dolap1bd2(:,:),dolaptot2(:,:)
      real*8, allocatable :: dolap2ak(:,:),dolap2bk(:,:)
      real*8, allocatable :: dolap2ak2(:,:),dolap2bk2(:,:)
      real*8, allocatable :: ddist2(:,:)
      logical docompute
      integer ibond,jlp,ilp
      real*8 forcetemp(3)
c
      allocate (dcos1(3,4))
      allocate (dcos2(3,4))
      allocate (dcos3(3,4))
      allocate (dcos4(3,4))
      allocate (dcos5(3,4))
      allocate (dcos6(3,4))
      allocate (dcos7(3,4))
      allocate (dcos8(3,4))
c
      allocate (dolap1ac(3,4))
      allocate (dolap1ad(3,4))
      allocate (dolap1bc(3,4))
      allocate (dolap1bd(3,4))
      allocate (dolaptot(3,6))
c
      allocate (dolap1ac2(3,4))
      allocate (dolap1ad2(3,4))
      allocate (dolap1bc2(3,4))
      allocate (dolap1bd2(3,4))
      allocate (dolaptot2(3,6))

      allocate (dolap2ak(3,5))
      allocate (dolap2bk(3,5))
      allocate (dolap2ak2(3,5))
      allocate (dolap2bk2(3,5))
c
      allocate (dcos21(3,5))
      allocate (dcos22(3,5))
      allocate (dcos23(3,5))
      allocate (dcos24(3,5))
      dcos21 = 0d0
      dcos22 = 0d0
      dcos23 = 0d0
      dcos24 = 0d0
      allocate (ddist2(3,6))
      ddist2 = 0d0
c
      allocate (dcos31(3,6))
      allocate (dcos32(3,6))
      dcos31 = 0d0
      dcos32 = 0d0

c
      di = 0d0
      dix = 0d0
      diz = 0d0
      di2 = 0d0
      dix2 = 0d0
      diz2 = 0d0
c
      erep1 = 0.0d0
      erep2 = 0.0d0
      erep3 = 0.0d0
      erep12 = 0.0d0
      erep22 = 0.0d0
      erep32 = 0.0d0
      testolap = 0
c
c     zero out the repulsion  energy and partitioning terms
c
      nerep = 0
      erep = 0.0d0
      call rotlp
c
c     calculate bond-bond repulsion term
c
      do inl = 1, nbondlocnl
         i = bondglobnl(inl)
         if (inl.eq.0) cycle
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
         xba = -xab
         yba = -yab
         zba = -zab
         rab2 = xab*xab + yab*yab + zab*zab
         rab = sqrt(rab2)
         do jj = 1, nbondlst(inl) 
            j = bondlst(jj,inl)
c
c  correspondance between indexes and number (1,2,3,4)
c  1 -> ia, 2-> ib, 3 -> jc, 4 -> jd
c
            dcos1 = 0d0
            dcos2 = 0d0
            dcos3 = 0d0
            dcos4 = 0d0
            dcos5 = 0d0
            dcos6 = 0d0
            dcos7 = 0d0
            dcos8 = 0d0
            dolap1ac = 0d0
            dolap1ad = 0d0
            dolap1bc = 0d0
            dolap1bd = 0d0
            dolaptot = 0d0
            dolap1ac2 = 0d0
            dolap1ad2 = 0d0
            dolap1bc2 = 0d0
            dolap1bd2 = 0d0
            dolaptot2 = 0d0
            jc = ibnd(1,j)                     
            jd = ibnd(2,j)
            if (molcule(ia) .ne. molcule(jc)) then
            xc = x(jc)  
            yc = y(jc)
            zc = z(jc)
            xd = x(jd)
            yd = y(jd)
            zd = z(jd)
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
            if (dist12.gt.(repcut*repcut)) cycle
c
            xcd = xd - xc
            ycd = yd - yc
            zcd = zd - zc
            xdc = -xcd
            ydc = -ycd
            zdc = -zcd
            rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
            rcd = sqrt(rcd2)
            xad = xd - xa
            yad = yd - ya
            zad = zd - za
            if (use_bounds)  call image (xad,yad,zad)
            xda = -xad
            yda = -yad
            zda = -zad
            rda2 = xad*xad + yad*yad + zad*zad
c
            rda = sqrt(rda2)
            xbd = xd - xb
            ybd = yd - yb
            zbd = zd - zb
            if (use_bounds)  call image (xbd,ybd,zbd)
            xdb = -xbd
            ydb = -ybd
            zdb = -zbd
            rbd2 = xbd*xbd + ybd*ybd + zbd*zbd         
            rbd = sqrt(rbd2)
            xac = xc - xa
            yac = yc - ya
            zac = zc - za
            if (use_bounds)  call image (xac,yac,zac)
            xca = -xac
            yca = -yac
            zca = -zac
            rac2 = xac*xac + yac*yac + zac*zac         
            rac = sqrt(rac2)
            xbc = xc - xb
            ybc = yc - yb
            zbc = zc - zb
            if (use_bounds)  call image (xbc,ybc,zbc)
            xcb = -xbc
            ycb = -ybc
            zcb = -zbc
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
            scabac = xab*xac + yab*yac + zab*zac
c
            dcos1(1,1) =  - (xac+xab)/(rab*rac) +
     $       scabac*xab/(rac*rab**3) + scabac*xac/(rab*rac**3)
            dcos1(2,1) =  - (yac+yab)/(rab*rac) +
     $       scabac*yab/(rac*rab**3) + scabac*yac/(rab*rac**3)
            dcos1(3,1) =  - (zac+zab)/(rab*rac) +
     $       scabac*zab/(rac*rab**3) + scabac*zac/(rab*rac**3)
     
            dcos1(1,2) =  xac/(rab*rac) -(scabac*xab)/(rac*rab**3)
            dcos1(2,2) =  yac/(rab*rac) -(scabac*yab)/(rac*rab**3)
            dcos1(3,2) =  zac/(rab*rac) -(scabac*zab)/(rac*rab**3)
     
            dcos1(1,3) = xab/(rab*rac) -(scabac*xac)/(rab*rac**3)
            dcos1(2,3) = yab/(rab*rac) -(scabac*yac)/(rab*rac**3)
            dcos1(3,3) = zab/(rab*rac) -(scabac*zac)/(rab*rac**3)
c
            cos12 = (xab*xad + yab*yad + zab*zad)/(rab*rda)
            scabad = xab*xad + yab*yad + zab*zad
            dcos2(1,1) =  - (xad+xab)/(rab*rda) +
     $       scabad*xab/(rda*rab**3) + scabad*xad/(rab*rda**3)
            dcos2(2,1) =  - (yad+yab)/(rab*rda) +
     $       scabad*yab/(rda*rab**3) + scabad*yad/(rab*rda**3)
            dcos2(3,1) =  - (zad+zab)/(rab*rda) +
     $       scabad*zab/(rda*rab**3) + scabad*zad/(rab*rda**3)
c
            dcos2(1,2) =  xad/(rab*rda) -(scabad*xab)/(rda*rab**3)
            dcos2(2,2) =  yad/(rab*rda) -(scabad*yab)/(rda*rab**3)
            dcos2(3,2) =  zad/(rab*rda) -(scabad*zab)/(rda*rab**3)

            dcos2(1,4) =  xab/(rab*rda) -(scabad*xad)/(rab*rda**3)
            dcos2(2,4) =  yab/(rab*rda) -(scabad*yad)/(rab*rda**3)
            dcos2(3,4) =  zab/(rab*rda) -(scabad*zad)/(rab*rda**3)
c
            cos13 = (xba*xbc + yba*ybc + zba*zbc)/(rab*rbc)
            scbabc = xba*xbc + yba*ybc + zba*zbc
            dcos3(1,2) =  - (xba+xbc)/(rab*rbc) +
     $       scbabc*xba/(rbc*rab**3) + scbabc*xbc/(rab*rbc**3)
            dcos3(2,2) =  - (yba+ybc)/(rab*rbc) +
     $       scbabc*yba/(rbc*rab**3) + scbabc*ybc/(rab*rbc**3)
            dcos3(3,2) =  - (zba+zbc)/(rab*rbc) +
     $       scbabc*zba/(rbc*rab**3) + scbabc*zbc/(rab*rbc**3)

            dcos3(1,1) =  xbc/(rab*rbc) -(scbabc*xba)/(rbc*rab**3)
            dcos3(2,1) =  ybc/(rab*rbc) -(scbabc*yba)/(rbc*rab**3)
            dcos3(3,1) =  zbc/(rab*rbc) -(scbabc*zba)/(rbc*rab**3)
            dcos3(1,3) =  xba/(rab*rbc) -(scbabc*xbc)/(rab*rbc**3)
            dcos3(2,3) =  yba/(rab*rbc) -(scbabc*ybc)/(rab*rbc**3)
            dcos3(3,3) =  zba/(rab*rbc) -(scbabc*zbc)/(rab*rbc**3)
c
            cos14 = (xba*xbd + yba*ybd + zba*zbd)/(rab*rbd)
            scbabd = xba*xbd + yba*ybd + zba*zbd
c
            dcos4(1,2) = - (xba+xbd)/(rab*rbd) +
     $       scbabd*xba/(rbd*rab**3) + scbabd*xbd/(rab*rbd**3)
            dcos4(2,2) =  - (yba+ybd)/(rab*rbd) +
     $       scbabd*yba/(rbd*rab**3) + scbabd*ybd/(rab*rbd**3)
            dcos4(3,2) =  - (zba+zbd)/(rab*rbd) +
     $       scbabd*zba/(rbd*rab**3) + scbabd*zbd/(rab*rbd**3)

            dcos4(1,1) =  xbd/(rab*rbd) -(scbabd*xba)/(rbd*rab**3)
            dcos4(2,1) =  ybd/(rab*rbd) -(scbabd*yba)/(rbd*rab**3)
            dcos4(3,1) =  zbd/(rab*rbd) -(scbabd*zba)/(rbd*rab**3)

            dcos4(1,4) =  xba/(rab*rbd) -(scbabd*xbd)/(rab*rbd**3)
            dcos4(2,4) =  yba/(rab*rbd) -(scbabd*ybd)/(rab*rbd**3)
            dcos4(3,4) =  zba/(rab*rbd) -(scbabd*zbd)/(rab*rbd**3)

            cos15 = (xcd*xca + ycd*yca + zcd*zca)/(rcd*rac)
            sccdca = xcd*xca + ycd*yca + zcd*zca
c
            dcos5(1,3) = - (xca+xcd)/(rcd*rac) +
     $       sccdca*xcd/(rac*rcd**3) + sccdca*xca/(rcd*rac**3)
            dcos5(2,3) = - (yca+ycd)/(rcd*rac) +
     $       sccdca*ycd/(rac*rcd**3) + sccdca*yca/(rcd*rac**3)
            dcos5(3,3) = - (zca+zcd)/(rcd*rac) +
     $       sccdca*zcd/(rac*rcd**3) + sccdca*zca/(rcd*rac**3)

            dcos5(1,4) =  xca/(rcd*rac) -(sccdca*xcd)/(rac*rcd**3)
            dcos5(2,4) =  yca/(rcd*rac) -(sccdca*ycd)/(rac*rcd**3)
            dcos5(3,4) =  zca/(rcd*rac) -(sccdca*zcd)/(rac*rcd**3)

            dcos5(1,1) =  xcd/(rcd*rac) -(sccdca*xca)/(rcd*rac**3)
            dcos5(2,1) =  ycd/(rcd*rac) -(sccdca*yca)/(rcd*rac**3)
            dcos5(3,1) =  zcd/(rcd*rac) -(sccdca*zca)/(rcd*rac**3)

            cos16 = (xcd*xcb + ycd*ycb + zcd*zcb)/(rcd*rbc)
            sccdcb = xcd*xcb + ycd*ycb + zcd*zcb

            dcos6(1,3) = - (xcb+xcd)/(rcd*rbc) +
     $       sccdcb*xcd/(rbc*rcd**3) + sccdcb*xcb/(rcd*rbc**3)
            dcos6(2,3) =  - (ycb+ycd)/(rcd*rbc) +
     $       sccdcb*ycd/(rbc*rcd**3) + sccdcb*ycb/(rcd*rbc**3)
            dcos6(3,3) =  - (zcb+zcd)/(rcd*rbc) +
     $       sccdcb*zcd/(rbc*rcd**3) + sccdcb*zcb/(rcd*rbc**3)
     
            dcos6(1,4) =  xcb/(rcd*rbc) -(sccdcb*xcd)/(rbc*rcd**3)
            dcos6(2,4) =  ycb/(rcd*rbc) -(sccdcb*ycd)/(rbc*rcd**3)
            dcos6(3,4) =  zcb/(rcd*rbc) -(sccdcb*zcd)/(rbc*rcd**3)

            dcos6(1,2) =  xcd/(rcd*rbc) -(sccdcb*xcb)/(rcd*rbc**3)
            dcos6(2,2) =  ycd/(rcd*rbc) -(sccdcb*ycb)/(rcd*rbc**3)
            dcos6(3,2) =  zcd/(rcd*rbc) -(sccdcb*zcb)/(rcd*rbc**3)
c

            cos17 = (xcd*xad + ycd*yad + zcd*zad)/(rcd*rda) 
            sccdad = xcd*xad + ycd*yad + zcd*zad
c
            dcos7(1,4) =  - (xdc+xda)/(rcd*rda) +
     $       sccdad*xdc/(rda*rcd**3) + sccdad*xda/(rcd*rda**3)
            dcos7(2,4) =  - (ydc+yda)/(rcd*rda) +
     $       sccdad*ydc/(rda*rcd**3) + sccdad*yda/(rcd*rda**3)
            dcos7(3,4) =  - (zdc+zda)/(rcd*rda) +
     $       sccdad*zdc/(rda*rcd**3) + sccdad*zda/(rcd*rda**3)

            dcos7(1,3) =  xda/(rcd*rda) -(sccdad*xdc)/(rda*rcd**3)
            dcos7(2,3) =  yda/(rcd*rda) -(sccdad*ydc)/(rda*rcd**3)
            dcos7(3,3) =  zda/(rcd*rda) -(sccdad*zdc)/(rda*rcd**3)

            dcos7(1,1) =  xdc/(rcd*rda) -(sccdad*xda)/(rcd*rda**3)
            dcos7(2,1) =  ydc/(rcd*rda) -(sccdad*yda)/(rcd*rda**3)
            dcos7(3,1) =  zdc/(rcd*rda) -(sccdad*zda)/(rcd*rda**3)
            
            cos18 = (xdc*xdb + ydc*ydb + zdc*zdb)/(rcd*rbd)
            scdcdb = xdc*xdb + ydc*ydb + zdc*zdb

            dcos8(1,4) = - (xdc+xdb)/(rcd*rbd) +
     $       scdcdb*xdc/(rbd*rcd**3) + scdcdb*xdb/(rcd*rbd**3)
            dcos8(2,4) = - (ydc+ydb)/(rcd*rbd) +
     $       scdcdb*ydc/(rbd*rcd**3) + scdcdb*ydb/(rcd*rbd**3)
            dcos8(3,4) = - (zdc+zdb)/(rcd*rbd) +
     $       scdcdb*zdc/(rbd*rcd**3) + scdcdb*zdb/(rcd*rbd**3)

            dcos8(1,3) =  xdb/(rcd*rbd) -(scdcdb*xdc)/(rbd*rcd**3)
            dcos8(2,3) =  ydb/(rcd*rbd) -(scdcdb*ydc)/(rbd*rcd**3)
            dcos8(3,3) =  zdb/(rcd*rbd) -(scdcdb*zdc)/(rbd*rcd**3)

            dcos8(1,2) =  xdc/(rcd*rbd) -(scdcdb*xdb)/(rcd*rbd**3)
            dcos8(2,2) =  ydc/(rcd*rbd) -(scdcdb*ydb)/(rcd*rbd**3)
            dcos8(3,2) =  zdc/(rcd*rbd) -(scdcdb*zdb)/(rcd*rbd**3)

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
            dist1 = sqrt(dist12) 
            e11 = onab*olap1*olap1*oncd/dist1
            e12 = onab*olap12*olap12*oncd/dist12
c
            erep =  erep + e11*cvrep11 + e12*cvrep12 
c              
            dolap1ac(1,1) =  sasc*(ss11+ps11+sp11+pp11)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*xac/rac) 
     &       + sasc*cpa*csc*psac*dcos1(1,1) 
     &       + sasc*csa*cpc*spac*dcos5(1,1) 
     &       + sasc*cpa*cpc*2*dcos1(1,1)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(1,1) 
            dolap1ac(2,1) =  sasc*(ss11+ps11+sp11+pp11)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*yac/rac)
     &       + sasc*cpa*csc*psac*dcos1(2,1) 
     &       + sasc*csa*cpc*spac*dcos5(2,1) 
     &       + sasc*cpa*cpc*2*dcos1(2,1)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(2,1) 
            dolap1ac(3,1) =  sasc*(ss11+ps11+sp11+pp11)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*zac/rac)
     &       + sasc*cpa*csc*psac*dcos1(3,1) 
     &       + sasc*csa*cpc*spac*dcos5(3,1) 
     &       + sasc*cpa*cpc*2*dcos1(3,1)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(3,1) 
c
            dolap1ac2(1,1) =  sasc2*(ss11+ps11+sp11+pp11)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*xac/rac) 
     &       + sasc2*cpa*csc*psac*dcos1(1,1) 
     &       + sasc2*csa*cpc*spac*dcos5(1,1) 
     &       + sasc2*cpa*cpc*2*dcos1(1,1)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(1,1) 
            dolap1ac2(2,1) =  sasc2*(ss11+ps11+sp11+pp11)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*yac/rac)
     &       + sasc2*cpa*csc*psac*dcos1(2,1) 
     &       + sasc2*csa*cpc*spac*dcos5(2,1) 
     &       + sasc2*cpa*cpc*2*dcos1(2,1)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(2,1) 
            dolap1ac2(3,1) =  sasc2*(ss11+ps11+sp11+pp11)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*zac/rac)
     &       + sasc2*cpa*csc*psac*dcos1(3,1) 
     &       + sasc2*csa*cpc*spac*dcos5(3,1) 
     &       + sasc2*cpa*cpc*2*dcos1(3,1)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(3,1) 
c
            dolap1ac(1,3) =  - sasc*(ss11+ps11+sp11+pp11)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*xac/rac)
     &       + sasc*cpa*csc*psac*dcos1(1,3) 
     &       + sasc*csa*cpc*spac*dcos5(1,3) 
     &       + sasc*cpa*cpc*2*dcos1(1,3)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(1,3) 
            dolap1ac(2,3) =  - sasc*(ss11+ps11+sp11+pp11)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*yac/rac)
     &       + sasc*cpa*csc*psac*dcos1(2,3) 
     &       + sasc*csa*cpc*spac*dcos5(2,3) 
     &       + sasc*cpa*cpc*2*dcos1(2,3)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(2,3) 
            dolap1ac(3,3) =  - sasc*(ss11+ps11+sp11+pp11)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*zac/rac) 
     &       + sasc*cpa*csc*psac*dcos1(3,3) 
     &       + sasc*csa*cpc*spac*dcos5(3,3) 
     &       + sasc*cpa*cpc*2*dcos1(3,3)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(3,3) 
c
            dolap1ac2(1,3) =  - sasc2*(ss11+ps11+sp11+pp11)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*xac/rac)
     &       + sasc2*cpa*csc*psac*dcos1(1,3) 
     &       + sasc2*csa*cpc*spac*dcos5(1,3) 
     &       + sasc2*cpa*cpc*2*dcos1(1,3)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(1,3) 
            dolap1ac2(2,3) =  - sasc2*(ss11+ps11+sp11+pp11)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*yac/rac)
     &       + sasc2*cpa*csc*psac*dcos1(2,3) 
     &       + sasc2*csa*cpc*spac*dcos5(2,3) 
     &       + sasc2*cpa*cpc*2*dcos1(2,3)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(2,3) 
            dolap1ac2(3,3) =  - sasc2*(ss11+ps11+sp11+pp11)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jc)))*zac/rac) 
     &       + sasc2*cpa*csc*psac*dcos1(3,3) 
     &       + sasc2*csa*cpc*spac*dcos5(3,3) 
     &       + sasc2*cpa*cpc*2*dcos1(3,3)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(3,3) 
c
            dolap1ac(1,2) =  sasc*cpa*csc*psac*dcos1(1,2) 
     &       + sasc*csa*cpc*spac*dcos5(1,2) 
     &       + sasc*cpa*cpc*2*dcos1(1,2)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(1,2) 
            dolap1ac(2,2) = sasc*cpa*csc*psac*dcos1(2,2) 
     &       + sasc*csa*cpc*spac*dcos5(2,2) 
     &       + sasc*cpa*cpc*2*dcos1(2,2)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(2,2) 
            dolap1ac(3,2) =  sasc*cpa*csc*psac*dcos1(3,2) 
     &       + sasc*csa*cpc*spac*dcos5(3,2) 
     &       + sasc*cpa*cpc*2*dcos1(3,2)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(3,2) 
c
            dolap1ac2(1,2) =  sasc2*cpa*csc*psac*dcos1(1,2) 
     &       + sasc2*csa*cpc*spac*dcos5(1,2) 
     &       + sasc2*cpa*cpc*2*dcos1(1,2)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(1,2) 
            dolap1ac2(2,2) = sasc2*cpa*csc*psac*dcos1(2,2) 
     &       + sasc2*csa*cpc*spac*dcos5(2,2) 
     &       + sasc2*cpa*cpc*2*dcos1(2,2)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(2,2) 
            dolap1ac2(3,2) =  sasc2*cpa*csc*psac*dcos1(3,2) 
     &       + sasc2*csa*cpc*spac*dcos5(3,2) 
     &       + sasc2*cpa*cpc*2*dcos1(3,2)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(3,2) 
c
            dolap1ac(1,4) =  sasc*cpa*csc*psac*dcos1(1,4) 
     &       + sasc*csa*cpc*spac*dcos5(1,4) 
     &       + sasc*cpa*cpc*2*dcos1(1,4)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(1,4) 
            dolap1ac(2,4) =  sasc*cpa*csc*psac*dcos1(2,4) 
     &       + sasc*csa*cpc*spac*dcos5(2,4) 
     &       + sasc*cpa*cpc*2*dcos1(2,4)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(2,4) 
            dolap1ac(3,4) =  sasc*cpa*csc*psac*dcos1(3,4) 
     &       + sasc*csa*cpc*spac*dcos5(3,4) 
     &       + sasc*cpa*cpc*2*dcos1(3,4)*cos15 
     &       + sasc*cpa*cpc*2*cos11*dcos5(3,4) 
c
            dolap1ac2(1,4) =  sasc2*cpa*csc*psac*dcos1(1,4) 
     &       + sasc2*csa*cpc*spac*dcos5(1,4) 
     &       + sasc2*cpa*cpc*2*dcos1(1,4)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(1,4) 
            dolap1ac2(2,4) =  sasc2*cpa*csc*psac*dcos1(2,4) 
     &       + sasc2*csa*cpc*spac*dcos5(2,4) 
     &       + sasc2*cpa*cpc*2*dcos1(2,4)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(2,4) 
            dolap1ac2(3,4) =  sasc2*cpa*csc*psac*dcos1(3,4) 
     &       + sasc2*csa*cpc*spac*dcos5(3,4) 
     &       + sasc2*cpa*cpc*2*dcos1(3,4)*cos15 
     &       + sasc2*cpa*cpc*2*cos11*dcos5(3,4) 
c
            dolap1ad(1,1) =  sasd*(ss12+ps12+sp12+pp12)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*xad/rda) 
     &       + sasd*cpa*csd*psad*dcos2(1,1) 
     &       + sasd*csa*cpd*spad*dcos7(1,1) 
     &       + sasd*cpa*cpd*2*dcos2(1,1)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(1,1) 
            dolap1ad(2,1) =  sasd*(ss12+ps12+sp12+pp12)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*yad/rda) 
     &       + sasd*cpa*csd*psad*dcos2(2,1) 
     &       + sasd*csa*cpd*spad*dcos7(2,1) 
     &       + sasd*cpa*cpd*2*dcos2(2,1)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(2,1) 
            dolap1ad(3,1) =  sasd*(ss12+ps12+sp12+pp12)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*zad/rda) 
     &       + sasd*cpa*csd*psad*dcos2(3,1) 
     &       + sasd*csa*cpd*spad*dcos7(3,1) 
     &       + sasd*cpa*cpd*2*dcos2(3,1)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(3,1) 
c
            dolap1ad2(1,1) =  sasd2*(ss12+ps12+sp12+pp12)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*xad/rda) 
     &       + sasd2*cpa*csd*psad*dcos2(1,1) 
     &       + sasd2*csa*cpd*spad*dcos7(1,1) 
     &       + sasd2*cpa*cpd*2*dcos2(1,1)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(1,1) 
            dolap1ad2(2,1) =  sasd2*(ss12+ps12+sp12+pp12)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*yad/rda) 
     &       + sasd2*cpa*csd*psad*dcos2(2,1) 
     &       + sasd2*csa*cpd*spad*dcos7(2,1) 
     &       + sasd2*cpa*cpd*2*dcos2(2,1)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(2,1) 
            dolap1ad2(3,1) =  sasd2*(ss12+ps12+sp12+pp12)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*zad/rda) 
     &       + sasd2*cpa*csd*psad*dcos2(3,1) 
     &       + sasd2*csa*cpd*spad*dcos7(3,1) 
     &       + sasd2*cpa*cpd*2*dcos2(3,1)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(3,1) 
c
            dolap1ad(1,3) =  sasd*cpa*csd*psad*dcos2(1,3) 
     &       + sasd*csa*cpd*spad*dcos7(1,3) 
     &       + sasd*cpa*cpd*2*dcos2(1,3)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(1,3) 
            dolap1ad(2,3) =  sasd*cpa*csd*psad*dcos2(2,3) 
     &       + sasd*csa*cpd*spad*dcos7(2,3) 
     &       + sasd*cpa*cpd*2*dcos2(2,3)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(2,3) 
            dolap1ad(3,3) =  sasd*cpa*csd*psad*dcos2(3,3) 
     &       + sasd*csa*cpd*spad*dcos7(3,3) 
     &       + sasd*cpa*cpd*2*dcos2(3,3)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(3,3) 
c
            dolap1ad2(1,3) =  sasd2*cpa*csd*psad*dcos2(1,3) 
     &       + sasd2*csa*cpd*spad*dcos7(1,3) 
     &       + sasd2*cpa*cpd*2*dcos2(1,3)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(1,3) 
            dolap1ad2(2,3) =  sasd2*cpa*csd*psad*dcos2(2,3) 
     &       + sasd2*csa*cpd*spad*dcos7(2,3) 
     &       + sasd2*cpa*cpd*2*dcos2(2,3)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(2,3) 
            dolap1ad2(3,3) =  sasd2*cpa*csd*psad*dcos2(3,3) 
     &       + sasd2*csa*cpd*spad*dcos7(3,3) 
     &       + sasd2*cpa*cpd*2*dcos2(3,3)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(3,3) 
c
            dolap1ad(1,2) =  sasd*cpa*csd*psad*dcos2(1,2) 
     &       + sasd*csa*cpd*spad*dcos7(1,2) 
     &       + sasd*cpa*cpd*2*dcos2(1,2)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(1,2) 
            dolap1ad(2,2) =  sasd*cpa*csd*psad*dcos2(2,2) 
     &       + sasd*csa*cpd*spad*dcos7(2,2) 
     &       + sasd*cpa*cpd*2*dcos2(2,2)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(2,2) 
            dolap1ad(3,2) =  sasd*cpa*csd*psad*dcos2(3,2) 
     &       + sasd*csa*cpd*spad*dcos7(3,2) 
     &       + sasd*cpa*cpd*2*dcos2(3,2)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(3,2) 
c
            dolap1ad2(1,2) =  sasd2*cpa*csd*psad*dcos2(1,2) 
     &       + sasd2*csa*cpd*spad*dcos7(1,2) 
     &       + sasd2*cpa*cpd*2*dcos2(1,2)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(1,2) 
            dolap1ad2(2,2) =  sasd2*cpa*csd*psad*dcos2(2,2) 
     &       + sasd2*csa*cpd*spad*dcos7(2,2) 
     &       + sasd2*cpa*cpd*2*dcos2(2,2)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(2,2) 
            dolap1ad2(3,2) =  sasd2*cpa*csd*psad*dcos2(3,2) 
     &       + sasd2*csa*cpd*spad*dcos7(3,2) 
     &       + sasd2*cpa*cpd*2*dcos2(3,2)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(3,2) 
cc
            dolap1ad(1,4) =  - sasd*(ss12+ps12+sp12+pp12)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*xad/rda) 
     &       + sasd*cpa*csd*psad*dcos2(1,4) 
     &       + sasd*csa*cpd*spad*dcos7(1,4) 
     &       + sasd*cpa*cpd*2*dcos2(1,4)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(1,4) 
            dolap1ad(2,4) =  - sasd*(ss12+ps12+sp12+pp12)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*yad/rda) 
     &       + sasd*cpa*csd*psad*dcos2(2,4) 
     &       + sasd*csa*cpd*spad*dcos7(2,4) 
     &       + sasd*cpa*cpd*2*dcos2(2,4)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(2,4) 
            dolap1ad(3,4) =  - sasd*(ss12+ps12+sp12+pp12)*
     &       (alpha/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*zad/rda) 
     &       + sasd*cpa*csd*psad*dcos2(3,4) 
     &       + sasd*csa*cpd*spad*dcos7(3,4) 
     &       + sasd*cpa*cpd*2*dcos2(3,4)*cos17 
     &       + sasd*cpa*cpd*2*cos12*dcos7(3,4) 
c
            dolap1ad2(1,4) =  - sasd2*(ss12+ps12+sp12+pp12)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*xad/rda) 
     &       + sasd2*cpa*csd*psad*dcos2(1,4) 
     &       + sasd2*csa*cpd*spad*dcos7(1,4) 
     &       + sasd2*cpa*cpd*2*dcos2(1,4)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(1,4) 
            dolap1ad2(2,4) =  - sasd2*(ss12+ps12+sp12+pp12)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*yad/rda) 
     &       + sasd2*cpa*csd*psad*dcos2(2,4) 
     &       + sasd2*csa*cpd*spad*dcos7(2,4) 
     &       + sasd2*cpa*cpd*2*dcos2(2,4)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(2,4) 
            dolap1ad2(3,4) =  - sasd2*(ss12+ps12+sp12+pp12)*
     &       (alpha2/(4*sqrt(vdwrep(ia)*vdwrep(jd)))*zad/rda) 
     &       + sasd2*cpa*csd*psad*dcos2(3,4) 
     &       + sasd2*csa*cpd*spad*dcos7(3,4) 
     &       + sasd2*cpa*cpd*2*dcos2(3,4)*cos17 
     &       + sasd2*cpa*cpd*2*cos12*dcos7(3,4) 
c
            dolap1bc(1,1) =  sbsc*cpb*csc*psbc*dcos3(1,1) 
     &       + sbsc*csb*cpc*spbc*dcos6(1,1) 
     &       + sbsc*cpb*cpc*2*dcos3(1,1)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(1,1) 
            dolap1bc(2,1) =  sbsc*cpb*csc*psbc*dcos3(2,1) 
     &       + sbsc*csb*cpc*spbc*dcos6(2,1) 
     &       + sbsc*cpb*cpc*2*dcos3(2,1)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(2,1) 
            dolap1bc(3,1) =  sbsc*cpb*csc*psbc*dcos3(3,1) 
     &       + sbsc*csb*cpc*spbc*dcos6(3,1) 
     &       + sbsc*cpb*cpc*2*dcos3(3,1)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(3,1) 
c
            dolap1bc2(1,1) =  sbsc2*cpb*csc*psbc*dcos3(1,1) 
     &       + sbsc2*csb*cpc*spbc*dcos6(1,1) 
     &       + sbsc2*cpb*cpc*2*dcos3(1,1)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(1,1) 
            dolap1bc2(2,1) =  sbsc2*cpb*csc*psbc*dcos3(2,1) 
     &       + sbsc2*csb*cpc*spbc*dcos6(2,1) 
     &       + sbsc2*cpb*cpc*2*dcos3(2,1)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(2,1) 
            dolap1bc2(3,1) =  sbsc2*cpb*csc*psbc*dcos3(3,1) 
     &       + sbsc2*csb*cpc*spbc*dcos6(3,1) 
     &       + sbsc2*cpb*cpc*2*dcos3(3,1)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(3,1) 
c
            dolap1bc(1,3) =  - sbsc*(ss13+ps13+sp13+pp13)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*xbc/rbc)
     &       + sbsc*cpb*csc*psbc*dcos3(1,3) 
     &       + sbsc*csb*cpc*spbc*dcos6(1,3) 
     &       + sbsc*cpb*cpc*2*dcos3(1,3)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(1,3) 
            dolap1bc(2,3) =  - sbsc*(ss13+ps13+sp13+pp13)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*ybc/rbc)
     &       + sbsc*cpb*csc*psbc*dcos3(2,3) 
     &       + sbsc*csb*cpc*spbc*dcos6(2,3) 
     &       + sbsc*cpb*cpc*2*dcos3(2,3)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(2,3) 
            dolap1bc(3,3) =  - sbsc*(ss13+ps13+sp13+pp13)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*zbc/rbc)
     &       + sbsc*cpb*csc*psbc*dcos3(3,3) 
     &       + sbsc*csb*cpc*spbc*dcos6(3,3) 
     &       + sbsc*cpb*cpc*2*dcos3(3,3)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(3,3) 
c
            dolap1bc2(1,3) =  - sbsc2*(ss13+ps13+sp13+pp13)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*xbc/rbc)
     &       + sbsc2*cpb*csc*psbc*dcos3(1,3) 
     &       + sbsc2*csb*cpc*spbc*dcos6(1,3) 
     &       + sbsc2*cpb*cpc*2*dcos3(1,3)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(1,3) 
            dolap1bc2(2,3) =  - sbsc2*(ss13+ps13+sp13+pp13)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*ybc/rbc)
     &       + sbsc2*cpb*csc*psbc*dcos3(2,3) 
     &       + sbsc2*csb*cpc*spbc*dcos6(2,3) 
     &       + sbsc2*cpb*cpc*2*dcos3(2,3)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(2,3) 
            dolap1bc2(3,3) =  - sbsc2*(ss13+ps13+sp13+pp13)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*zbc/rbc)
     &       + sbsc2*cpb*csc*psbc*dcos3(3,3) 
     &       + sbsc2*csb*cpc*spbc*dcos6(3,3) 
     &       + sbsc2*cpb*cpc*2*dcos3(3,3)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(3,3) 
c
            dolap1bc(1,2) =  sbsc*(ss13+ps13+sp13+pp13)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*xbc/rbc)
     &       + sbsc*cpb*csc*psbc*dcos3(1,2) 
     &       + sbsc*csb*cpc*spbc*dcos6(1,2) 
     &       + sbsc*cpb*cpc*2*dcos3(1,2)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(1,2) 
            dolap1bc(2,2) =  sbsc*(ss13+ps13+sp13+pp13)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*ybc/rbc)
     &       + sbsc*cpb*csc*psbc*dcos3(2,2) 
     &       + sbsc*csb*cpc*spbc*dcos6(2,2) 
     &       + sbsc*cpb*cpc*2*dcos3(2,2)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(2,2) 
            dolap1bc(3,2) =  sbsc*(ss13+ps13+sp13+pp13)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*zbc/rbc)
     &       + sbsc*cpb*csc*psbc*dcos3(3,2) 
     &       + sbsc*csb*cpc*spbc*dcos6(3,2) 
     &       + sbsc*cpb*cpc*2*dcos3(3,2)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(3,2) 
c
            dolap1bc2(1,2) =  sbsc2*(ss13+ps13+sp13+pp13)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*xbc/rbc)
     &       + sbsc2*cpb*csc*psbc*dcos3(1,2) 
     &       + sbsc2*csb*cpc*spbc*dcos6(1,2) 
     &       + sbsc2*cpb*cpc*2*dcos3(1,2)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(1,2) 
            dolap1bc2(2,2) =  sbsc2*(ss13+ps13+sp13+pp13)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*ybc/rbc)
     &       + sbsc2*cpb*csc*psbc*dcos3(2,2) 
     &       + sbsc2*csb*cpc*spbc*dcos6(2,2) 
     &       + sbsc2*cpb*cpc*2*dcos3(2,2)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(2,2) 
            dolap1bc2(3,2) =  sbsc2*(ss13+ps13+sp13+pp13)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jc)))*zbc/rbc)
     &       + sbsc2*cpb*csc*psbc*dcos3(3,2) 
     &       + sbsc2*csb*cpc*spbc*dcos6(3,2) 
     &       + sbsc2*cpb*cpc*2*dcos3(3,2)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(3,2) 
c
            dolap1bc(1,4) =  sbsc*cpb*csc*psbc*dcos3(1,4) 
     &       + sbsc*csb*cpc*spbc*dcos6(1,4) 
     &       + sbsc*cpb*cpc*2*dcos3(1,4)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(1,4) 
            dolap1bc(2,4) =  sbsc*cpb*csc*psbc*dcos3(2,4) 
     &       + sbsc*csb*cpc*spbc*dcos6(2,4) 
     &       + sbsc*cpb*cpc*2*dcos3(2,4)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(2,4) 
            dolap1bc(3,4) =  sbsc*cpb*csc*psbc*dcos3(3,4) 
     &       + sbsc*csb*cpc*spbc*dcos6(3,4) 
     &       + sbsc*cpb*cpc*2*dcos3(3,4)*cos16 
     &       + sbsc*cpb*cpc*2*cos13*dcos6(3,4) 
c
            dolap1bc2(1,4) =  sbsc2*cpb*csc*psbc*dcos3(1,4) 
     &       + sbsc2*csb*cpc*spbc*dcos6(1,4) 
     &       + sbsc2*cpb*cpc*2*dcos3(1,4)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(1,4) 
            dolap1bc2(2,4) =  sbsc2*cpb*csc*psbc*dcos3(2,4) 
     &       + sbsc2*csb*cpc*spbc*dcos6(2,4) 
     &       + sbsc2*cpb*cpc*2*dcos3(2,4)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(2,4) 
            dolap1bc2(3,4) =  sbsc2*cpb*csc*psbc*dcos3(3,4) 
     &       + sbsc2*csb*cpc*spbc*dcos6(3,4) 
     &       + sbsc2*cpb*cpc*2*dcos3(3,4)*cos16 
     &       + sbsc2*cpb*cpc*2*cos13*dcos6(3,4) 
c
            dolap1bd(1,2) =  sbsd*(ss14+ps14+sp14+pp14)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*xbd/rbd)
     &       + sbsd*cpb*csd*psbd*dcos4(1,2) 
     &       + sbsd*csb*cpd*spbd*dcos8(1,2) 
     &       + sbsd*cpb*cpd*2*dcos8(1,2)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(1,2) 
            dolap1bd(2,2) =  sbsd*(ss14+ps14+sp14+pp14)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*ybd/rbd)
     &       + sbsd*cpb*csd*psbd*dcos4(2,2) 
     &       + sbsd*csb*cpd*spbd*dcos8(2,2) 
     &       + sbsd*cpb*cpd*2*dcos8(2,2)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(2,2) 
            dolap1bd(3,2) =  sbsd*(ss14+ps14+sp14+pp14)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*zbd/rbd)
     &       + sbsd*cpb*csd*psbd*dcos4(3,2) 
     &       + sbsd*csb*cpd*spbd*dcos8(3,2) 
     &       + sbsd*cpb*cpd*2*dcos8(3,2)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(3,2) 
c
            dolap1bd2(1,2) =  sbsd2*(ss14+ps14+sp14+pp14)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*xbd/rbd)
     &       + sbsd2*cpb*csd*psbd*dcos4(1,2) 
     &       + sbsd2*csb*cpd*spbd*dcos8(1,2) 
     &       + sbsd2*cpb*cpd*2*dcos8(1,2)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(1,2) 
            dolap1bd2(2,2) =  sbsd2*(ss14+ps14+sp14+pp14)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*ybd/rbd)
     &       + sbsd2*cpb*csd*psbd*dcos4(2,2) 
     &       + sbsd2*csb*cpd*spbd*dcos8(2,2) 
     &       + sbsd2*cpb*cpd*2*dcos8(2,2)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(2,2) 
            dolap1bd2(3,2) =  sbsd2*(ss14+ps14+sp14+pp14)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*zbd/rbd)
     &       + sbsd2*cpb*csd*psbd*dcos4(3,2) 
     &       + sbsd2*csb*cpd*spbd*dcos8(3,2) 
     &       + sbsd2*cpb*cpd*2*dcos8(3,2)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(3,2) 
c
            dolap1bd(1,4) =  - sbsd*(ss14+ps14+sp14+pp14)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*xbd/rbd)
     &       + sbsd*cpb*csd*psbd*dcos4(1,4) 
     &       + sbsd*csb*cpd*spbd*dcos8(1,4) 
     &       + sbsd*cpb*cpd*2*dcos8(1,4)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(1,4) 
            dolap1bd(2,4) =  - sbsd*(ss14+ps14+sp14+pp14)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*ybd/rbd)
     &       + sbsd*cpb*csd*psbd*dcos4(2,4) 
     &       + sbsd*csb*cpd*spbd*dcos8(2,4) 
     &       + sbsd*cpb*cpd*2*dcos8(2,4)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(2,4) 
            dolap1bd(3,4) = - sbsd*(ss14+ps14+sp14+pp14)*
     &       (alpha/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*zbd/rbd)
     &       + sbsd*cpb*csd*psbd*dcos4(3,4) 
     &       + sbsd*csb*cpd*spbd*dcos8(3,4) 
     &       + sbsd*cpb*cpd*2*dcos8(3,4)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(3,4) 
c
            dolap1bd2(1,4) =  - sbsd2*(ss14+ps14+sp14+pp14)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*xbd/rbd)
     &       + sbsd2*cpb*csd*psbd*dcos4(1,4) 
     &       + sbsd2*csb*cpd*spbd*dcos8(1,4) 
     &       + sbsd2*cpb*cpd*2*dcos8(1,4)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(1,4) 
            dolap1bd2(2,4) =  - sbsd2*(ss14+ps14+sp14+pp14)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*ybd/rbd)
     &       + sbsd2*cpb*csd*psbd*dcos4(2,4) 
     &       + sbsd2*csb*cpd*spbd*dcos8(2,4) 
     &       + sbsd2*cpb*cpd*2*dcos8(2,4)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(2,4) 
            dolap1bd2(3,4) = - sbsd2*(ss14+ps14+sp14+pp14)*
     &       (alpha2/(4*sqrt(vdwrep(ib)*vdwrep(jd)))*zbd/rbd)
     &       + sbsd2*cpb*csd*psbd*dcos4(3,4) 
     &       + sbsd2*csb*cpd*spbd*dcos8(3,4) 
     &       + sbsd2*cpb*cpd*2*dcos8(3,4)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(3,4) 
c
             dolap1bd(1,1) =  sbsd*cpb*csd*psbd*dcos4(1,1) 
     &       + sbsd*csb*cpd*spbd*dcos8(1,1) 
     &       + sbsd*cpb*cpd*2*dcos8(1,1)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(1,1) 
             dolap1bd(2,1) =  sbsd*cpb*csd*psbd*dcos4(2,1) 
     &       + sbsd*csb*cpd*spbd*dcos8(2,1) 
     &       + sbsd*cpb*cpd*2*dcos8(2,1)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(2,1) 
             dolap1bd(3,1) =  sbsd*cpb*csd*psbd*dcos4(3,1) 
     &       + sbsd*csb*cpd*spbd*dcos8(3,1) 
     &       + sbsd*cpb*cpd*2*dcos8(3,1)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(3,1) 
c
             dolap1bd2(1,1) =  sbsd2*cpb*csd*psbd*dcos4(1,1) 
     &       + sbsd2*csb*cpd*spbd*dcos8(1,1) 
     &       + sbsd2*cpb*cpd*2*dcos8(1,1)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(1,1) 
             dolap1bd2(2,1) =  sbsd2*cpb*csd*psbd*dcos4(2,1) 
     &       + sbsd2*csb*cpd*spbd*dcos8(2,1) 
     &       + sbsd2*cpb*cpd*2*dcos8(2,1)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(2,1) 
             dolap1bd2(3,1) =  sbsd2*cpb*csd*psbd*dcos4(3,1) 
     &       + sbsd2*csb*cpd*spbd*dcos8(3,1) 
     &       + sbsd2*cpb*cpd*2*dcos8(3,1)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(3,1) 
c
             dolap1bd(1,3) =  sbsd*cpb*csd*psbd*dcos4(1,3) 
     &       + sbsd*csb*cpd*spbd*dcos8(1,3) 
     &       + sbsd*cpb*cpd*2*dcos8(1,3)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(1,3) 
             dolap1bd(2,3) =  sbsd*cpb*csd*psbd*dcos4(2,3) 
     &       + sbsd*csb*cpd*spbd*dcos8(2,3) 
     &       + sbsd*cpb*cpd*2*dcos8(2,3)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(2,3) 
             dolap1bd(3,3) =  sbsd*cpb*csd*psbd*dcos4(3,3) 
     &       + sbsd*csb*cpd*spbd*dcos8(3,3) 
     &       + sbsd*cpb*cpd*2*dcos8(3,3)*cos14 
     &       + sbsd*cpb*cpd*2*cos18*dcos4(3,3) 
c
             dolap1bd2(1,3) =  sbsd2*cpb*csd*psbd*dcos4(1,3) 
     &       + sbsd2*csb*cpd*spbd*dcos8(1,3) 
     &       + sbsd2*cpb*cpd*2*dcos8(1,3)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(1,3) 
             dolap1bd2(2,3) =  sbsd2*cpb*csd*psbd*dcos4(2,3) 
     &       + sbsd2*csb*cpd*spbd*dcos8(2,3) 
     &       + sbsd2*cpb*cpd*2*dcos8(2,3)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(2,3) 
             dolap1bd2(3,3) =  sbsd2*cpb*csd*psbd*dcos4(3,3) 
     &       + sbsd2*csb*cpd*spbd*dcos8(3,3) 
     &       + sbsd2*cpb*cpd*2*dcos8(3,3)*cos14 
     &       + sbsd2*cpb*cpd*2*cos18*dcos4(3,3) 
c
              dolaptot(1,1) = dolap1ac(1,1) +
     $         dolap1ad(1,1) + dolap1bc(1,1) + dolap1bd(1,1)
              dolaptot(2,1) = dolap1ac(2,1) +
     $         dolap1ad(2,1) + dolap1bc(2,1) + dolap1bd(2,1)
              dolaptot(3,1) = dolap1ac(3,1) +
     $         dolap1ad(3,1) + dolap1bc(3,1) + dolap1bd(3,1)
c
              dolaptot2(1,1) = dolap1ac2(1,1) +
     $         dolap1ad2(1,1) + dolap1bc2(1,1) + dolap1bd2(1,1)
              dolaptot2(2,1) = dolap1ac2(2,1) +
     $         dolap1ad2(2,1) + dolap1bc2(2,1) + dolap1bd2(2,1)
              dolaptot2(3,1) = dolap1ac2(3,1) +
     $         dolap1ad2(3,1) + dolap1bc2(3,1) + dolap1bd2(3,1)
c
              dolaptot(1,2) = dolap1ac(1,2) +
     $         dolap1ad(1,2) + dolap1bc(1,2) + dolap1bd(1,2)
              dolaptot(2,2) = dolap1ac(2,2) +
     $         dolap1ad(2,2) + dolap1bc(2,2) + dolap1bd(2,2)
              dolaptot(3,2) = dolap1ac(3,2) +
     $         dolap1ad(3,2) + dolap1bc(3,2) + dolap1bd(3,2)
c
              dolaptot2(1,2) = dolap1ac2(1,2) +
     $         dolap1ad2(1,2) + dolap1bc2(1,2) + dolap1bd2(1,2)
              dolaptot2(2,2) = dolap1ac2(2,2) +
     $         dolap1ad2(2,2) + dolap1bc2(2,2) + dolap1bd2(2,2)
              dolaptot2(3,2) = dolap1ac2(3,2) +
     $         dolap1ad2(3,2) + dolap1bc2(3,2) + dolap1bd2(3,2)
c
              dolaptot(1,3) = dolap1ac(1,3) +
     $         dolap1ad(1,3) + dolap1bc(1,3) + dolap1bd(1,3)
              dolaptot(2,3) = dolap1ac(2,3) +
     $         dolap1ad(2,3) + dolap1bc(2,3) + dolap1bd(2,3)
              dolaptot(3,3) = dolap1ac(3,3) +
     $         dolap1ad(3,3) + dolap1bc(3,3) + dolap1bd(3,3)
c
              dolaptot2(1,3) = dolap1ac2(1,3) +
     $         dolap1ad2(1,3) + dolap1bc2(1,3) + dolap1bd2(1,3)
              dolaptot2(2,3) = dolap1ac2(2,3) +
     $         dolap1ad2(2,3) + dolap1bc2(2,3) + dolap1bd2(2,3)
              dolaptot2(3,3) = dolap1ac2(3,3) +
     $         dolap1ad2(3,3) + dolap1bc2(3,3) + dolap1bd2(3,3)
c
              dolaptot(1,4) = dolap1ac(1,4) +
     $         dolap1ad(1,4) + dolap1bc(1,4) + dolap1bd(1,4)
              dolaptot(2,4) = dolap1ac(2,4) +
     $         dolap1ad(2,4) + dolap1bc(2,4) + dolap1bd(2,4)
              dolaptot(3,4) = dolap1ac(3,4) +
     $         dolap1ad(3,4) + dolap1bc(3,4) + dolap1bd(3,4)
c
              dolaptot2(1,4) = dolap1ac2(1,4) +
     $         dolap1ad2(1,4) + dolap1bc2(1,4) + dolap1bd2(1,4)
              dolaptot2(2,4) = dolap1ac2(2,4) +
     $         dolap1ad2(2,4) + dolap1bc2(2,4) + dolap1bd2(2,4)
              dolaptot2(3,4) = dolap1ac2(3,4) +
     $         dolap1ad2(3,4) + dolap1bc2(3,4) + dolap1bd2(3,4)
c
              forcetemp(1) = cvrep11*onab*oncd*
     $         (2*dolaptot(1,1)*olap1*(1/dist1) +
     $         olap1**2*(-xdist1/(2*dist1**3)))
              forcetemp(2) = cvrep11*onab*oncd*
     $         (2*dolaptot(2,1)*olap1*(1/dist1) +
     $         olap1**2*(-ydist1/(2*dist1**3)))
              forcetemp(3) = cvrep11*onab*oncd*
     $         (2*dolaptot(3,1)*olap1*(1/dist1) +
     $         olap1**2*(-zdist1/(2*dist1**3)))
              forcetemp(1) = forcetemp(1) + cvrep12*onab*oncd*
     $         (2*dolaptot2(1,1)*olap12*(1/dist12) +
     $         olap12**2*(-xdist1/dist1**4))
              forcetemp(2) = forcetemp(2) + cvrep12*onab*oncd*
     $         (2*dolaptot2(2,1)*olap12*(1/dist12) +
     $         olap12**2*(-ydist1/dist1**4))
              forcetemp(3) = forcetemp(3) + cvrep12*onab*oncd*
     $         (2*dolaptot2(3,1)*olap12*(1/dist12) +
     $         olap12**2*(-zdist1/dist1**4))
              ialoc = loc(ia)
              derep(1,ialoc) = derep(1,ialoc) + forcetemp(1) 
              derep(2,ialoc) = derep(2,ialoc) + forcetemp(2) 
              derep(3,ialoc) = derep(3,ialoc) + forcetemp(3) 
c
              forcetemp(1) = cvrep11*onab*oncd*
     $         (2*dolaptot(1,2)*olap1*(1/dist1) +
     $         olap1**2*(-xdist1/(2*dist1**3)))
              forcetemp(2) = cvrep11*onab*oncd*
     $         (2*dolaptot(2,2)*olap1*(1/dist1) +
     $         olap1**2*(-ydist1/(2*dist1**3)))
              forcetemp(3) = cvrep11*onab*oncd*
     $         (2*dolaptot(3,2)*olap1*(1/dist1) +
     $         olap1**2*(-zdist1/(2*dist1**3)))
              forcetemp(1) = forcetemp(1) + cvrep12*onab*oncd*
     $         (2*dolaptot2(1,2)*olap12*(1/dist12) +
     $         olap12**2*(-xdist1/dist1**4))
              forcetemp(2) = forcetemp(2) + cvrep12*onab*oncd*
     $         (2*dolaptot2(2,2)*olap12*(1/dist12) +
     $         olap12**2*(-ydist1/dist1**4))
              forcetemp(3) = forcetemp(3) + cvrep12*onab*oncd*
     $         (2*dolaptot2(3,2)*olap12*(1/dist12) +
     $         olap12**2*(-zdist1/dist1**4))
              ibloc = loc(ib)
              derep(1,ibloc) = derep(1,ibloc) + forcetemp(1) 
              derep(2,ibloc) = derep(2,ibloc) + forcetemp(2) 
              derep(3,ibloc) = derep(3,ibloc) + forcetemp(3) 
c
              forcetemp(1) =  cvrep11*onab*oncd*
     $         (2*dolaptot(1,3)*olap1*(1/dist1) -
     $         olap1**2*(-xdist1/(2*dist1**3)))
              forcetemp(2) =  cvrep11*onab*oncd*
     $         (2*dolaptot(2,3)*olap1*(1/dist1) -
     $         olap1**2*(-ydist1/(2*dist1**3)))
              forcetemp(3) =  cvrep11*onab*oncd*
     $         (2*dolaptot(3,3)*olap1*(1/dist1) -
     $         olap1**2*(-zdist1/(2*dist1**3)))
              forcetemp(1) = forcetemp(1) + cvrep12*onab*oncd*
     $         (2*dolaptot2(1,3)*olap12*(1/dist12) -
     $         olap12**2*(-xdist1/dist1**4))
              forcetemp(2) = forcetemp(2) + cvrep12*onab*oncd*
     $         (2*dolaptot2(2,3)*olap12*(1/dist12) -
     $         olap12**2*(-ydist1/dist1**4))
              forcetemp(3) = forcetemp(3) + cvrep12*onab*oncd*
     $         (2*dolaptot2(3,3)*olap12*(1/dist12) -
     $         olap12**2*(-zdist1/dist1**4))
              jcloc = loc(jc)
              derep(1,jcloc) = derep(1,jcloc) + forcetemp(1) 
              derep(2,jcloc) = derep(2,jcloc) + forcetemp(2) 
              derep(3,jcloc) = derep(3,jcloc) + forcetemp(3) 
c
              forcetemp(1) =  cvrep11*onab*oncd*
     $         (2*dolaptot(1,4)*olap1*(1/dist1) -
     $         olap1**2*(-xdist1/(2*dist1**3)))
              forcetemp(2) =  cvrep11*onab*oncd*
     $         (2*dolaptot(2,4)*olap1*(1/dist1) -
     $         olap1**2*(-ydist1/(2*dist1**3)))
              forcetemp(3) = cvrep11*onab*oncd*
     $         (2*dolaptot(3,4)*olap1*(1/dist1) -
     $         olap1**2*(-zdist1/(2*dist1**3)))
              forcetemp(1) = forcetemp(1) + cvrep12*onab*oncd*
     $         (2*dolaptot2(1,4)*olap12*(1/dist12) -
     $         olap12**2*(-xdist1/dist1**4))
              forcetemp(2) = forcetemp(2) + cvrep12*onab*oncd*
     $         (2*dolaptot2(2,4)*olap12*(1/dist12) -
     $         olap12**2*(-ydist1/dist1**4))
              forcetemp(3) = forcetemp(3) + cvrep12*onab*oncd*
     $         (2*dolaptot2(3,4)*olap12*(1/dist12) -
     $         olap12**2*(-zdist1/dist1**4))
              jdloc = loc(jd)
              derep(1,jdloc) = derep(1,jdloc) + forcetemp(1) 
              derep(2,jdloc) = derep(2,jdloc) + forcetemp(2) 
              derep(3,jdloc) = derep(3,jdloc) + forcetemp(3) 
            end if
            nerep = nerep + 1
         end do
      end do   
c           
c     Compute Bond - Lone Pair repulsion
c
      do inl = 1, nbondlocnl
        i = bondglobnl(inl)
        if (inl.eq.0) cycle
        do ii = 1, nbondlplst(inl)
            jlp = bondlplst(ii,inl)
            k = lpatom(jlp)       
            call rotlp1(jlp,di,dix,diz)
c
c  correspondance between indexes and number (1,2,3,4,5)
c  1 -> ia, 2-> ib, 3 -> k, 4 -> ixlp(j), 5 -> izlp(j) 
c
            dcos21 = 0d0
            dcos22 = 0d0
            dcos23 = 0d0
            dcos24 = 0d0
            dolap2ak = 0d0
            dolap2bk = 0d0
            dolap2ak2 = 0d0
            dolap2bk2 = 0d0
            dolaptot = 0d0
            dolaptot2 = 0d0
            ddist2 = 0d0
            ia = ibnd(1,i)                            
            ib = ibnd(2,i)
            xa = x(ia)  
            ya = y(ia)
            za = z(ia)
            xb = x(ib)
            yb = y(ib)
            zb = z(ib)
            xab = xb - xa
            yab = yb - ya
            zab = zb - za
c            if (use_bounds)  call image (xab,yab,zab)
            xba = -xab
            yba = -yab
            zba = -zab
            rab2 = xab*xab + yab*yab + zab*zab
c
c             if (rab2.gt.(repcut*repcut)) cycle
c
            rab = sqrt(rab2)
            cia = rpole(1,ia) 
            cib = rpole(1,ib) 
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
             xka = -xak
             yka = -yak
             zka = -zak
             rak2 = xak*xak + yak*yak + zak*zak
             rak = sqrt(rak2)
             xbk = xk - xb
             ybk = yk - yb
             zbk = zk - zb
             if (use_bounds)  call image (xbk,ybk,zbk)
             xkb = -xbk
             ykb = -ybk
             zkb = -zbk
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
            
             cos21 =  (xka*xklp + yka*yklp + zka*zklp)/(rak*rklp)
             sckaklp = xka*xklp + yka*yklp + zka*zklp

             dcos21(1,1) = xklp/(rak*rklp) - 
     $      xka*sckaklp/(rklp*rak**3)
             dcos21(2,1) = yklp/(rak*rklp) - 
     $      yka*sckaklp/(rklp*rak**3)
             dcos21(3,1) = zklp/(rak*rklp) - 
     $      zka*sckaklp/(rklp*rak**3)
     
            dcos21(1,3) = - (xka+xklp)/(rak*rklp) + 
     $      xklp*sckaklp/(rak*rklp**3) + 
     $      xka*sckaklp/(rklp*rak**3)
            dcos21(2,3) = - (yka+yklp)/(rak*rklp) + 
     $      yklp*sckaklp/(rak*rklp**3) + 
     $      yka*sckaklp/(rklp*rak**3)
            dcos21(3,3) =  -(zka+zklp)/(rak*rklp) + 
     $      zklp*sckaklp/(rak*rklp**3) + 
     $      zka*sckaklp/(rklp*rak**3)
c
c         contribution of the rotation matrixes
c
             dcoslp1 = xka/(rak*rklp) - xklp*sckaklp/(rak*rklp**3) 
             dcoslp2 = yka/(rak*rklp) - yklp*sckaklp/(rak*rklp**3) 
             dcoslp3 = zka/(rak*rklp) - zklp*sckaklp/(rak*rklp**3) 
c
             dcos21(1,4) =  
     $    dix(1,1)*dcoslp1 + dix(1,2)*dcoslp2 + dix(1,3)*dcoslp3
             dcos21(2,4) =  
     $    dix(2,1)*dcoslp1 + dix(2,2)*dcoslp2 + dix(2,3)*dcoslp3
             dcos21(3,4) =  
     $    dix(3,1)*dcoslp1 + dix(3,2)*dcoslp2 + dix(3,3)*dcoslp3
c
             dcos21(1,5) =  
     $    diz(1,1)*dcoslp1 + diz(1,2)*dcoslp2 + diz(1,3)*dcoslp3
             dcos21(2,5) =  
     $    diz(2,1)*dcoslp1 + diz(2,2)*dcoslp2 + diz(2,3)*dcoslp3
             dcos21(3,5) =  
     $    diz(3,1)*dcoslp1 + diz(3,2)*dcoslp2 + diz(3,3)*dcoslp3
c
             dcos21(1,3) = dcos21(1,3) + 
     $    di(1,1)*dcoslp1 + di(1,2)*dcoslp2 + di(1,3)*dcoslp3
             dcos21(2,3) = dcos21(2,3) +  
     $    di(2,1)*dcoslp1 + di(2,2)*dcoslp2 + di(2,3)*dcoslp3
             dcos21(3,3) = dcos21(3,3) + 
     $    di(3,1)*dcoslp1 + di(3,2)*dcoslp2 + di(3,3)*dcoslp3
c
             cos22 =  (xkb*xklp + ykb*yklp + zkb*zklp)/(rbk*rklp)
             sckbklp = xkb*xklp + ykb*yklp + zkb*zklp
             
             dcos22(1,2) =  xklp/(rbk*rklp) - 
     $          xkb*sckbklp/(rklp*rbk**3)
             dcos22(2,2) =  yklp/(rbk*rklp) - 
     $          ykb*sckbklp/(rklp*rbk**3)
             dcos22(3,2) =  zklp/(rbk*rklp) - 
     $          zkb*sckbklp/(rklp*rbk**3)

             dcos22(1,3) = - (xkb+xklp)/(rbk*rklp) + 
     $       xklp*sckbklp/(rbk*rklp**3) + 
     $       xkb*sckbklp/(rklp*rbk**3)
             dcos22(2,3) = - (ykb+yklp)/(rbk*rklp) + 
     $       yklp*sckbklp/(rbk*rklp**3) + 
     $       ykb*sckbklp/(rklp*rbk**3)
             dcos22(3,3) =  - (zkb+zklp)/(rbk*rklp) + 
     $       zklp*sckbklp/(rbk*rklp**3) + 
     $       zkb*sckbklp/(rklp*rbk**3)
c
c         contribution of the rotation matrixes
c
             dcoslp1 = xkb/(rbk*rklp) - xklp*sckbklp/(rbk*rklp**3) 
             dcoslp2 = ykb/(rbk*rklp) - yklp*sckbklp/(rbk*rklp**3) 
             dcoslp3 = zkb/(rbk*rklp) - zklp*sckbklp/(rbk*rklp**3) 
c
             dcos22(1,4) =  
     $    dix(1,1)*dcoslp1 + dix(1,2)*dcoslp2 + dix(1,3)*dcoslp3
             dcos22(2,4) =  
     $    dix(2,1)*dcoslp1 + dix(2,2)*dcoslp2 + dix(2,3)*dcoslp3
             dcos22(3,4) =  
     $    dix(3,1)*dcoslp1 + dix(3,2)*dcoslp2 + dix(3,3)*dcoslp3
c
             dcos22(1,5) =  
     $    diz(1,1)*dcoslp1 + diz(1,2)*dcoslp2 + diz(1,3)*dcoslp3
             dcos22(2,5) =  
     $    diz(2,1)*dcoslp1 + diz(2,2)*dcoslp2 + diz(2,3)*dcoslp3
             dcos22(3,5) =  
     $    diz(3,1)*dcoslp1 + diz(3,2)*dcoslp2 + diz(3,3)*dcoslp3
c
             dcos22(1,3) = dcos22(1,3) + 
     $    di(1,1)*dcoslp1 + di(1,2)*dcoslp2 + di(1,3)*dcoslp3
             dcos22(2,3) = dcos22(2,3) +  
     $    di(2,1)*dcoslp1 + di(2,2)*dcoslp2 + di(2,3)*dcoslp3
             dcos22(3,3) = dcos22(3,3) + 
     $    di(3,1)*dcoslp1 + di(3,2)*dcoslp2 + di(3,3)*dcoslp3
           
              sp21 = csa*cpk*cos21*psak
              sp22 = csb*cpk*cos22*psbk
              
              cos23 = (xab*xak + yab*yak + zab*zak)/(rab*rak)  
              scabak = xab*xak + yab*yak + zab*zak
              dcos23(1,1) = - (xak+xab)/(rab*rak) +
     $       scabak*xab/(rak*rab**3) + scabak*xak/(rab*rak**3)
              dcos23(2,1) = - (yak+yab)/(rab*rak) +
     $       scabak*yab/(rak*rab**3) + scabak*yak/(rab*rak**3)
              dcos23(3,1) = - (zak+zab)/(rab*rak) +
     $       scabak*zab/(rak*rab**3) + scabak*zak/(rab*rak**3)
c     
              dcos23(1,2) =  xak/(rab*rak) -
     $       (scabak*xab)/(rak*rab**3)
              dcos23(2,2) =  yak/(rab*rak) -
     $       (scabak*yab)/(rak*rab**3)
              dcos23(3,2) =  zak/(rab*rak) -
     $       (scabak*zab)/(rak*rab**3)
c     
              dcos23(1,3) =  xab/(rab*rak) -(scabak*xak)/(rab*rak**3)
              dcos23(2,3) =  yab/(rab*rak) -(scabak*yak)/(rab*rak**3)
              dcos23(3,3) =  zab/(rab*rak) -(scabak*zak)/(rab*rak**3)

              cos24 =  (xba*xbk + yba*ybk + zba*zbk)/(rab*rbk)  
              scbabk = xba*xbk + yba*ybk + zba*zbk
c
              dcos24(1,2) =  - (xba+xbk)/(rab*rbk) +
     $         scbabk*xba/(rbk*rab**3) + scbabk*xbk/(rab*rbk**3)
              dcos24(2,2) =  - (yba+ybk)/(rab*rbk) +
     $         scbabk*yba/(rbk*rab**3) + scbabk*ybk/(rab*rbk**3)
              dcos24(3,2) =  - (zba+zbk)/(rab*rbk) +
     $         scbabk*zba/(rbk*rab**3) + scbabk*zbk/(rab*rbk**3)
c
              dcos24(1,1) =  xbk/(rab*rbk) -(scbabk*xba)/(rbk*rab**3)
              dcos24(2,1) =  ybk/(rab*rbk) -(scbabk*yba)/(rbk*rab**3)
              dcos24(3,1) =  zbk/(rab*rbk) -(scbabk*zba)/(rbk*rab**3)
c
              dcos24(1,3) =  xba/(rab*rbk) -(scbabk*xbk)/(rab*rbk**3)
              dcos24(2,3) =  yba/(rab*rbk) -(scbabk*ybk)/(rab*rbk**3)
              dcos24(3,3) =  zba/(rab*rbk) -(scbabk*zbk)/(rab*rbk**3)
            
              ps21 = cpa*csk*cos23*spak
              ps22 = cpb*csk*cos24*spbk
            
              pp21 = cpa*cpk*cos21*cos23*2 
              pp22 = cpb*cpk*cos22*cos24*2
            
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
              dolap2ak(1,1) = sask*(ss21+ps21+sp21+pp21)*
     $        (alpha/(4*sqrt(vdwrep(ia)*vdwreplp))*xak/rak)
     &        + sask*cpa*csk*spak*dcos23(1,1) 
     &        + sask*csa*cpk*psak*dcos21(1,1) 
     &        + sask*cpa*cpk*2*dcos21(1,1)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(1,1) 
              dolap2ak(2,1) = sask*(ss21+ps21+sp21+pp21)*
     $        (alpha/(4*sqrt(vdwrep(ia)*vdwreplp))*yak/rak)
     &        + sask*cpa*csk*spak*dcos23(2,1) 
     &        + sask*csa*cpk*psak*dcos21(2,1) 
     &        + sask*cpa*cpk*2*dcos21(2,1)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(2,1) 
              dolap2ak(3,1) = sask*(ss21+ps21+sp21+pp21)*
     $        (alpha/(4*sqrt(vdwrep(ia)*vdwreplp))*zak/rak)
     &        + sask*cpa*csk*spak*dcos23(3,1) 
     &        + sask*csa*cpk*psak*dcos21(3,1) 
     &        + sask*cpa*cpk*2*dcos21(3,1)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(3,1) 
c
              dolap2ak2(1,1) = sask2*(ss21+ps21+sp21+pp21)*
     $        (alpha2/(4*sqrt(vdwrep(ia)*vdwreplp))*xak/rak)
     &        + sask2*cpa*csk*spak*dcos23(1,1) 
     &        + sask2*csa*cpk*psak*dcos21(1,1) 
     &        + sask2*cpa*cpk*2*dcos21(1,1)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(1,1) 
              dolap2ak2(2,1) = sask2*(ss21+ps21+sp21+pp21)*
     $        (alpha2/(4*sqrt(vdwrep(ia)*vdwreplp))*yak/rak)
     &        + sask2*cpa*csk*spak*dcos23(2,1) 
     &        + sask2*csa*cpk*psak*dcos21(2,1) 
     &        + sask2*cpa*cpk*2*dcos21(2,1)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(2,1) 
              dolap2ak2(3,1) = sask2*(ss21+ps21+sp21+pp21)*
     $        (alpha2/(4*sqrt(vdwrep(ia)*vdwreplp))*zak/rak)
     &        + sask2*cpa*csk*spak*dcos23(3,1) 
     &        + sask2*csa*cpk*psak*dcos21(3,1) 
     &        + sask2*cpa*cpk*2*dcos21(3,1)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(3,1) 
c 
              dolap2ak(1,3) = -sask*(ss21+ps21+sp21+pp21)*
     $        (alpha/(4*sqrt(vdwrep(ia)*vdwreplp))*xak/rak)
     &        + sask*cpa*csk*spak*dcos23(1,3) 
     &        + sask*csa*cpk*psak*dcos21(1,3) 
     &        + sask*cpa*cpk*2*dcos21(1,3)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(1,3) 
              dolap2ak(2,3) = -sask*(ss21+ps21+sp21+pp21)*
     $        (alpha/(4*sqrt(vdwrep(ia)*vdwreplp))*yak/rak)
     &        + sask*cpa*csk*spak*dcos23(2,3) 
     &        + sask*csa*cpk*psak*dcos21(2,3) 
     &        + sask*cpa*cpk*2*dcos21(2,3)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(2,3) 
              dolap2ak(3,3) = -sask*(ss21+ps21+sp21+pp21)*
     $        (alpha/(4*sqrt(vdwrep(ia)*vdwreplp))*zak/rak)
     &        + sask*cpa*csk*spak*dcos23(3,3) 
     &        + sask*csa*cpk*psak*dcos21(3,3) 
     &        + sask*cpa*cpk*2*dcos21(3,3)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(3,3) 
c
              dolap2ak2(1,3) = -sask2*(ss21+ps21+sp21+pp21)*
     $        (alpha2/(4*sqrt(vdwrep(ia)*vdwreplp))*xak/rak)
     &        + sask2*cpa*csk*spak*dcos23(1,3) 
     &        + sask2*csa*cpk*psak*dcos21(1,3) 
     &        + sask2*cpa*cpk*2*dcos21(1,3)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(1,3) 
              dolap2ak2(2,3) = -sask2*(ss21+ps21+sp21+pp21)*
     $        (alpha2/(4*sqrt(vdwrep(ia)*vdwreplp))*yak/rak)
     &        + sask2*cpa*csk*spak*dcos23(2,3) 
     &        + sask2*csa*cpk*psak*dcos21(2,3) 
     &        + sask2*cpa*cpk*2*dcos21(2,3)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(2,3) 
              dolap2ak2(3,3) = -sask2*(ss21+ps21+sp21+pp21)*
     $        (alpha2/(4*sqrt(vdwrep(ia)*vdwreplp))*zak/rak)
     &        + sask2*cpa*csk*spak*dcos23(3,3) 
     &        + sask2*csa*cpk*psak*dcos21(3,3) 
     &        + sask2*cpa*cpk*2*dcos21(3,3)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(3,3) 
c
              dolap2ak(1,2) = 
     &        sask*cpa*csk*spak*dcos23(1,2) 
     &        + sask*csa*cpk*psak*dcos21(1,2) 
     &        + sask*cpa*cpk*2*dcos21(1,2)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(1,2) 
              dolap2ak(2,2) = 
     &        sask*cpa*csk*spak*dcos23(2,2) 
     &        + sask*csa*cpk*psak*dcos21(2,2) 
     &        + sask*cpa*cpk*2*dcos21(2,2)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(2,2) 
              dolap2ak(3,2) = 
     &        sask*cpa*csk*spak*dcos23(3,2) 
     &        + sask*csa*cpk*psak*dcos21(3,2) 
     &        + sask*cpa*cpk*2*dcos21(3,2)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(3,2) 
c
              dolap2ak2(1,2) = 
     &        sask2*cpa*csk*spak*dcos23(1,2) 
     &        + sask2*csa*cpk*psak*dcos21(1,2) 
     &        + sask2*cpa*cpk*2*dcos21(1,2)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(1,2) 
              dolap2ak2(2,2) = 
     &        sask2*cpa*csk*spak*dcos23(2,2) 
     &        + sask2*csa*cpk*psak*dcos21(2,2) 
     &        + sask2*cpa*cpk*2*dcos21(2,2)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(2,2) 
              dolap2ak2(3,2) = 
     &        sask2*cpa*csk*spak*dcos23(3,2) 
     &        + sask2*csa*cpk*psak*dcos21(3,2) 
     &        + sask2*cpa*cpk*2*dcos21(3,2)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(3,2) 
c
              dolap2ak(1,4) = 
     &        sask*cpa*csk*spak*dcos23(1,4) 
     &        + sask*csa*cpk*psak*dcos21(1,4) 
     &        + sask*cpa*cpk*2*dcos21(1,4)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(1,4) 
              dolap2ak(2,4) = 
     &        sask*cpa*csk*spak*dcos23(2,4) 
     &        + sask*csa*cpk*psak*dcos21(2,4) 
     &        + sask*cpa*cpk*2*dcos21(2,4)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(2,4) 
              dolap2ak(3,4) = 
     &        sask*cpa*csk*spak*dcos23(3,4) 
     &        + sask*csa*cpk*psak*dcos21(3,4) 
     &        + sask*cpa*cpk*2*dcos21(3,4)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(3,4) 
c
              dolap2ak2(1,4) = 
     &        sask2*cpa*csk*spak*dcos23(1,4) 
     &        + sask2*csa*cpk*psak*dcos21(1,4) 
     &        + sask2*cpa*cpk*2*dcos21(1,4)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(1,4) 
              dolap2ak2(2,4) = 
     &        sask2*cpa*csk*spak*dcos23(2,4) 
     &        + sask2*csa*cpk*psak*dcos21(2,4) 
     &        + sask2*cpa*cpk*2*dcos21(2,4)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(2,4) 
              dolap2ak2(3,4) = 
     &        sask2*cpa*csk*spak*dcos23(3,4) 
     &        + sask2*csa*cpk*psak*dcos21(3,4) 
     &        + sask2*cpa*cpk*2*dcos21(3,4)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(3,4) 
c
              dolap2ak(1,5) = 
     &        sask*cpa*csk*spak*dcos23(1,5) 
     &        + sask*csa*cpk*psak*dcos21(1,5) 
     &        + sask*cpa*cpk*2*dcos21(1,5)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(1,5) 
              dolap2ak(2,5) = 
     &        sask*cpa*csk*spak*dcos23(2,5) 
     &        + sask*csa*cpk*psak*dcos21(2,5) 
     &        + sask*cpa*cpk*2*dcos21(2,5)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(2,5) 
              dolap2ak(3,5) = 
     &        sask*cpa*csk*spak*dcos23(3,5) 
     &        + sask*csa*cpk*psak*dcos21(3,5) 
     &        + sask*cpa*cpk*2*dcos21(3,5)*cos23 
     &        + sask*cpa*cpk*2*cos21*dcos23(3,5) 
c
              dolap2ak2(1,5) = 
     &        sask2*cpa*csk*spak*dcos23(1,5) 
     &        + sask2*csa*cpk*psak*dcos21(1,5) 
     &        + sask2*cpa*cpk*2*dcos21(1,5)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(1,5) 
              dolap2ak2(2,5) = 
     &        sask2*cpa*csk*spak*dcos23(2,5) 
     &        + sask2*csa*cpk*psak*dcos21(2,5) 
     &        + sask2*cpa*cpk*2*dcos21(2,5)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(2,5) 
              dolap2ak2(3,5) = 
     &        sask2*cpa*csk*spak*dcos23(3,5) 
     &        + sask2*csa*cpk*psak*dcos21(3,5) 
     &        + sask2*cpa*cpk*2*dcos21(3,5)*cos23 
     &        + sask2*cpa*cpk*2*cos21*dcos23(3,5) 
c
c
              dolap2bk(1,2) = sbsk*(ss22+ps22+sp22+pp22)*
     $        (alpha/(4*sqrt(vdwrep(ib)*vdwreplp))*xbk/rbk)
     &        + sbsk*cpb*csk*spbk*dcos24(1,2) 
     &        + sbsk*csb*cpk*psbk*dcos22(1,2) 
     &        + sbsk*cpb*cpk*2*dcos22(1,2)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(1,2)
              dolap2bk(2,2) = sbsk*(ss22+ps22+sp22+pp22)*
     $        (alpha/(4*sqrt(vdwrep(ib)*vdwreplp))*ybk/rbk)
     &        + sbsk*cpb*csk*spbk*dcos24(2,2) 
     &        + sbsk*csb*cpk*psbk*dcos22(2,2) 
     &        + sbsk*cpb*cpk*2*dcos22(2,2)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(2,2)
              dolap2bk(3,2) = sbsk*(ss22+ps22+sp22+pp22)*
     $        (alpha/(4*sqrt(vdwrep(ib)*vdwreplp))*zbk/rbk)
     &        + sbsk*cpb*csk*spbk*dcos24(3,2) 
     &        + sbsk*csb*cpk*psbk*dcos22(3,2) 
     &        + sbsk*cpb*cpk*2*dcos22(3,2)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(3,2)
c
              dolap2bk2(1,2) = sbsk2*(ss22+ps22+sp22+pp22)*
     $        (alpha2/(4*sqrt(vdwrep(ib)*vdwreplp))*xbk/rbk)
     &        + sbsk2*cpb*csk*spbk*dcos24(1,2) 
     &        + sbsk2*csb*cpk*psbk*dcos22(1,2) 
     &        + sbsk2*cpb*cpk*2*dcos22(1,2)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(1,2)
              dolap2bk2(2,2) = sbsk2*(ss22+ps22+sp22+pp22)*
     $        (alpha2/(4*sqrt(vdwrep(ib)*vdwreplp))*ybk/rbk)
     &        + sbsk2*cpb*csk*spbk*dcos24(2,2) 
     &        + sbsk2*csb*cpk*psbk*dcos22(2,2) 
     &        + sbsk2*cpb*cpk*2*dcos22(2,2)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(2,2)
              dolap2bk2(3,2) = sbsk2*(ss22+ps22+sp22+pp22)*
     $        (alpha2/(4*sqrt(vdwrep(ib)*vdwreplp))*zbk/rbk)
     &        + sbsk2*cpb*csk*spbk*dcos24(3,2) 
     &        + sbsk2*csb*cpk*psbk*dcos22(3,2) 
     &        + sbsk2*cpb*cpk*2*dcos22(3,2)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(3,2)
c
              dolap2bk(1,3) = -sbsk*(ss22+ps22+sp22+pp22)*
     $        (alpha/(4*sqrt(vdwrep(ib)*vdwreplp))*xbk/rbk)
     &        + sbsk*cpb*csk*spbk*dcos24(1,3) 
     &        + sbsk*csb*cpk*psbk*dcos22(1,3) 
     &        + sbsk*cpb*cpk*2*dcos22(1,3)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(1,3)
              dolap2bk(2,3) = -sbsk*(ss22+ps22+sp22+pp22)*
     $        (alpha/(4*sqrt(vdwrep(ib)*vdwreplp))*ybk/rbk)
     &        + sbsk*cpb*csk*spbk*dcos24(2,3) 
     &        + sbsk*csb*cpk*psbk*dcos22(2,3) 
     &        + sbsk*cpb*cpk*2*dcos22(2,3)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(2,3)
              dolap2bk(3,3) = -sbsk*(ss22+ps22+sp22+pp22)*
     $        (alpha/(4*sqrt(vdwrep(ib)*vdwreplp))*zbk/rbk)
     &        + sbsk*cpb*csk*spbk*dcos24(3,3) 
     &        + sbsk*csb*cpk*psbk*dcos22(3,3) 
     &        + sbsk*cpb*cpk*2*dcos22(3,3)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(3,3)
c
              dolap2bk2(1,3) = -sbsk2*(ss22+ps22+sp22+pp22)*
     $        (alpha2/(4*sqrt(vdwrep(ib)*vdwreplp))*xbk/rbk)
     &        + sbsk2*cpb*csk*spbk*dcos24(1,3) 
     &        + sbsk2*csb*cpk*psbk*dcos22(1,3) 
     &        + sbsk2*cpb*cpk*2*dcos22(1,3)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(1,3)
              dolap2bk2(2,3) = -sbsk2*(ss22+ps22+sp22+pp22)*
     $        (alpha2/(4*sqrt(vdwrep(ib)*vdwreplp))*ybk/rbk)
     &        + sbsk2*cpb*csk*spbk*dcos24(2,3) 
     &        + sbsk2*csb*cpk*psbk*dcos22(2,3) 
     &        + sbsk2*cpb*cpk*2*dcos22(2,3)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(2,3)
              dolap2bk2(3,3) = -sbsk2*(ss22+ps22+sp22+pp22)*
     $        (alpha2/(4*sqrt(vdwrep(ib)*vdwreplp))*zbk/rbk)
     &        + sbsk2*cpb*csk*spbk*dcos24(3,3) 
     &        + sbsk2*csb*cpk*psbk*dcos22(3,3) 
     &        + sbsk2*cpb*cpk*2*dcos22(3,3)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(3,3)
c
              dolap2bk(1,1) = 
     &        sbsk*cpb*csk*spbk*dcos24(1,1) 
     &        + sbsk*csb*cpk*psbk*dcos22(1,1) 
     &        + sbsk*cpb*cpk*2*dcos22(1,1)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(1,1)
              dolap2bk(2,1) = 
     &        sbsk*cpb*csk*spbk*dcos24(2,1) 
     &        + sbsk*csb*cpk*psbk*dcos22(2,1) 
     &        + sbsk*cpb*cpk*2*dcos22(2,1)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(2,1)
              dolap2bk(3,1) = 
     &        sbsk*cpb*csk*spbk*dcos24(3,1) 
     &        + sbsk*csb*cpk*psbk*dcos22(3,1) 
     &        + sbsk*cpb*cpk*2*dcos22(3,1)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(3,1)
c
              dolap2bk2(1,1) = 
     &        sbsk2*cpb*csk*spbk*dcos24(1,1) 
     &        + sbsk2*csb*cpk*psbk*dcos22(1,1) 
     &        + sbsk2*cpb*cpk*2*dcos22(1,1)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(1,1)
              dolap2bk2(2,1) = 
     &        sbsk2*cpb*csk*spbk*dcos24(2,1) 
     &        + sbsk2*csb*cpk*psbk*dcos22(2,1) 
     &        + sbsk2*cpb*cpk*2*dcos22(2,1)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(2,1)
              dolap2bk2(3,1) = 
     &        sbsk2*cpb*csk*spbk*dcos24(3,1) 
     &        + sbsk2*csb*cpk*psbk*dcos22(3,1) 
     &        + sbsk2*cpb*cpk*2*dcos22(3,1)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(3,1)
c
              dolap2bk(1,4) = 
     &        sbsk*cpb*csk*spbk*dcos24(1,4) 
     &        + sbsk*csb*cpk*psbk*dcos22(1,4) 
     &        + sbsk*cpb*cpk*2*dcos22(1,4)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(1,4)
              dolap2bk(2,4) = 
     &        sbsk*cpb*csk*spbk*dcos24(2,4) 
     &        + sbsk*csb*cpk*psbk*dcos22(2,4) 
     &        + sbsk*cpb*cpk*2*dcos22(2,4)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(2,4)
              dolap2bk(3,4) = 
     &        sbsk*cpb*csk*spbk*dcos24(3,4) 
     &        + sbsk*csb*cpk*psbk*dcos22(3,4) 
     &        + sbsk*cpb*cpk*2*dcos22(3,4)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(3,4)
c
              dolap2bk2(1,4) = 
     &        sbsk2*cpb*csk*spbk*dcos24(1,4) 
     &        + sbsk2*csb*cpk*psbk*dcos22(1,4) 
     &        + sbsk2*cpb*cpk*2*dcos22(1,4)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(1,4)
              dolap2bk2(2,4) = 
     &        sbsk2*cpb*csk*spbk*dcos24(2,4) 
     &        + sbsk2*csb*cpk*psbk*dcos22(2,4) 
     &        + sbsk2*cpb*cpk*2*dcos22(2,4)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(2,4)
              dolap2bk2(3,4) = 
     &        sbsk2*cpb*csk*spbk*dcos24(3,4) 
     &        + sbsk2*csb*cpk*psbk*dcos22(3,4) 
     &        + sbsk2*cpb*cpk*2*dcos22(3,4)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(3,4)
c
              dolap2bk(1,5) = 
     &        sbsk*cpb*csk*spbk*dcos24(1,5) 
     &        + sbsk*csb*cpk*psbk*dcos22(1,5) 
     &        + sbsk*cpb*cpk*2*dcos22(1,5)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(1,5)
              dolap2bk(2,5) = 
     &        sbsk*cpb*csk*spbk*dcos24(2,5) 
     &        + sbsk*csb*cpk*psbk*dcos22(2,5) 
     &        + sbsk*cpb*cpk*2*dcos22(2,5)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(2,5)
              dolap2bk(3,5) = 
     &        sbsk*cpb*csk*spbk*dcos24(3,5) 
     &        + sbsk*csb*cpk*psbk*dcos22(3,5) 
     &        + sbsk*cpb*cpk*2*dcos22(3,5)*cos24 
     &        + sbsk*cpb*cpk*2*cos22*dcos24(3,5)
c
              dolap2bk2(1,5) = 
     &        sbsk2*cpb*csk*spbk*dcos24(1,5) 
     &        + sbsk2*csb*cpk*psbk*dcos22(1,5) 
     &        + sbsk2*cpb*cpk*2*dcos22(1,5)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(1,5)
              dolap2bk2(2,5) = 
     &        sbsk2*cpb*csk*spbk*dcos24(2,5) 
     &        + sbsk2*csb*cpk*psbk*dcos22(2,5) 
     &        + sbsk2*cpb*cpk*2*dcos22(2,5)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(2,5)
              dolap2bk2(3,5) = 
     &        sbsk2*cpb*csk*spbk*dcos24(3,5) 
     &        + sbsk2*csb*cpk*psbk*dcos22(3,5) 
     &        + sbsk2*cpb*cpk*2*dcos22(3,5)*cos24 
     &        + sbsk2*cpb*cpk*2*cos22*dcos24(3,5)
c
              dolaptot(1,1) = dolap2ak(1,1) + dolap2bk(1,1)
              dolaptot(2,1) = dolap2ak(2,1) + dolap2bk(2,1)
              dolaptot(3,1) = dolap2ak(3,1) + dolap2bk(3,1)
c
              dolaptot2(1,1) = dolap2ak2(1,1) + dolap2bk2(1,1)
              dolaptot2(2,1) = dolap2ak2(2,1) + dolap2bk2(2,1)
              dolaptot2(3,1) = dolap2ak2(3,1) + dolap2bk2(3,1)
c
              dolaptot(1,2) = dolap2ak(1,2) + dolap2bk(1,2)
              dolaptot(2,2) = dolap2ak(2,2) + dolap2bk(2,2)
              dolaptot(3,2) = dolap2ak(3,2) + dolap2bk(3,2)
c
              dolaptot2(1,2) = dolap2ak2(1,2) + dolap2bk2(1,2)
              dolaptot2(2,2) = dolap2ak2(2,2) + dolap2bk2(2,2)
              dolaptot2(3,2) = dolap2ak2(3,2) + dolap2bk2(3,2)
c
              dolaptot(1,3) = dolap2ak(1,3) + dolap2bk(1,3)
              dolaptot(2,3) = dolap2ak(2,3) + dolap2bk(2,3)
              dolaptot(3,3) = dolap2ak(3,3) + dolap2bk(3,3)
c
              dolaptot2(1,3) = dolap2ak2(1,3) + dolap2bk2(1,3)
              dolaptot2(2,3) = dolap2ak2(2,3) + dolap2bk2(2,3)
              dolaptot2(3,3) = dolap2ak2(3,3) + dolap2bk2(3,3)
c
              dolaptot(1,4) = dolap2ak(1,4) + dolap2bk(1,4)
              dolaptot(2,4) = dolap2ak(2,4) + dolap2bk(2,4)
              dolaptot(3,4) = dolap2ak(3,4) + dolap2bk(3,4)
c
              dolaptot2(1,4) = dolap2ak2(1,4) + dolap2bk2(1,4)
              dolaptot2(2,4) = dolap2ak2(2,4) + dolap2bk2(2,4)
              dolaptot2(3,4) = dolap2ak2(3,4) + dolap2bk2(3,4)
c
              dolaptot(1,5) = dolap2ak(1,5) + dolap2bk(1,5)
              dolaptot(2,5) = dolap2ak(2,5) + dolap2bk(2,5)
              dolaptot(3,5) = dolap2ak(3,5) + dolap2bk(3,5)
c
              dolaptot2(1,5) = dolap2ak2(1,5) + dolap2bk2(1,5)
              dolaptot2(2,5) = dolap2ak2(2,5) + dolap2bk2(2,5)
              dolaptot2(3,5) = dolap2ak2(3,5) + dolap2bk2(3,5)
c
c
              ddist2(1,1) =  0.5*xdist2/dist2 
              ddist2(2,1) =  0.5*ydist2/dist2 
              ddist2(3,1) =  0.5*zdist2/dist2 

              ddist2(1,2) =  0.5*xdist2/dist2 
              ddist2(2,2) =  0.5*ydist2/dist2 
              ddist2(3,2) =  0.5*zdist2/dist2 

              dlp1 = -xdist2/dist2 
              dlp2 = -ydist2/dist2 
              dlp3 = -zdist2/dist2 
c
              ddist2(1,3) = ddist2(1,3) +
     $         di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
              ddist2(2,3) = ddist2(2,3) +
     $         di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
              ddist2(3,3) = ddist2(3,3) +
     $         di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3
c
              ddist2(1,4) =  
     $         dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
              ddist2(2,4) = 
     $         dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
              ddist2(3,4) = 
     $         dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3
c
              ddist2(1,5) = 
     $         diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
              ddist2(2,5) = 
     $         diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
              ddist2(3,5) =
     $         diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3
c
              forcetemp(1) = cvrep21*onab*onlp*
     $          (2*dolaptot(1,1)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(1,1)/dist22)) 
              forcetemp(2) = cvrep21*onab*onlp*
     $          (2*dolaptot(2,1)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(2,1)/dist22)) 
              forcetemp(3) = cvrep21*onab*onlp*
     $          (2*dolaptot(3,1)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(3,1)/dist22)) 
              forcetemp(1) = forcetemp(1) + cvrep22*onab*onlp*
     $          (2*dolaptot2(1,1)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(1,1)/(dist2**3))
              forcetemp(2) = forcetemp(2) + cvrep22*onab*onlp*
     $          (2*dolaptot2(2,1)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(2,1)/(dist2**3))
              forcetemp(3) = forcetemp(3) + cvrep22*onab*onlp*
     $          (2*dolaptot2(3,1)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(3,1)/(dist2**3))
              ialoc = loc(ia)
              derep(1,ialoc) = derep(1,ialoc) + forcetemp(1) 
              derep(2,ialoc) = derep(2,ialoc) + forcetemp(2) 
              derep(3,ialoc) = derep(3,ialoc) + forcetemp(3) 
c
              forcetemp(1) = cvrep21*onab*onlp*
     $          (2*dolaptot(1,2)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(1,2)/dist22)) 
              forcetemp(2) = cvrep21*onab*onlp*
     $          (2*dolaptot(2,2)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(2,2)/dist22)) 
              forcetemp(3) = cvrep21*onab*onlp*
     $          (2*dolaptot(3,2)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(3,2)/dist22))
              forcetemp(1) = forcetemp(1) + cvrep22*onab*onlp*
     $          (2*dolaptot2(1,2)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(1,2)/(dist2**3))
              forcetemp(2) = forcetemp(2) + cvrep22*onab*onlp*
     $          (2*dolaptot2(2,2)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(2,2)/(dist2**3))
              forcetemp(3) = forcetemp(3) + cvrep22*onab*onlp*
     $          (2*dolaptot2(3,2)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(3,2)/(dist2**3))
              ibloc = loc(ib)
              derep(1,ibloc) = derep(1,ibloc) + forcetemp(1) 
              derep(2,ibloc) = derep(2,ibloc) + forcetemp(2) 
              derep(3,ibloc) = derep(3,ibloc) + forcetemp(3) 
c
              forcetemp(1) = cvrep21*onab*onlp*
     $          (2*dolaptot(1,3)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(1,3)/dist22)) 
              forcetemp(2) = cvrep21*onab*onlp*
     $          (2*dolaptot(2,3)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(2,3)/dist22)) 
              forcetemp(3) = cvrep21*onab*onlp*
     $          (2*dolaptot(3,3)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(3,3)/dist22)) 
              forcetemp(1) = forcetemp(1) + cvrep22*onab*onlp*
     $          (2*dolaptot2(1,3)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(1,3)/(dist2**3))
              forcetemp(2) = forcetemp(2) + cvrep22*onab*onlp*
     $          (2*dolaptot2(2,3)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(2,3)/(dist2**3))
              forcetemp(3) = forcetemp(3) + cvrep22*onab*onlp*
     $          (2*dolaptot2(3,3)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(3,3)/(dist2**3))
              kloc = loc(k)
              derep(1,kloc) = derep(1,kloc) + forcetemp(1) 
              derep(2,kloc) = derep(2,kloc) + forcetemp(2) 
              derep(3,kloc) = derep(3,kloc) + forcetemp(3) 
c
              forcetemp(1) = cvrep21*onab*onlp*
     $          (2*dolaptot(1,4)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(1,4)/dist22)) 
              forcetemp(2) = cvrep21*onab*onlp*
     $          (2*dolaptot(2,4)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(2,4)/dist22)) 
              forcetemp(3) = cvrep21*onab*onlp*
     $          (2*dolaptot(3,4)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(3,4)/dist22)) 
              forcetemp(1) = forcetemp(1) +cvrep22*onab*onlp*
     $          (2*dolaptot2(1,4)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(1,4)/(dist2**3))
              forcetemp(2) = forcetemp(2) +cvrep22*onab*onlp*
     $          (2*dolaptot2(2,4)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(2,4)/(dist2**3))
              forcetemp(3) = forcetemp(3) +cvrep22*onab*onlp*
     $          (2*dolaptot2(3,4)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(3,4)/(dist2**3))
c
              ixlploc = loc(ixlp(jlp))
              derep(1,ixlploc) = derep(1,ixlploc) + forcetemp(1)
              derep(2,ixlploc) = derep(2,ixlploc) + forcetemp(2)
              derep(3,ixlploc) = derep(3,ixlploc) + forcetemp(3)
c
              forcetemp(1) = cvrep21*onab*onlp*
     $          (2*dolaptot(1,5)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(1,5)/dist22))
              forcetemp(2) = cvrep21*onab*onlp*
     $          (2*dolaptot(2,5)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(2,5)/dist22))
              forcetemp(3) = cvrep21*onab*onlp*
     $          (2*dolaptot(3,5)*olap2*(1/dist2)
     $          -olap2**2*(ddist2(3,5)/dist22))
              forcetemp(1) = forcetemp(1) +cvrep22*onab*onlp*
     $          (2*dolaptot2(1,5)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(1,5)/(dist2**3))
              forcetemp(2) = forcetemp(2) +cvrep22*onab*onlp*
     $          (2*dolaptot2(2,5)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(2,5)/(dist2**3))
              forcetemp(3) = forcetemp(3) +cvrep22*onab*onlp*
     $          (2*dolaptot2(3,5)*olap22*(1/dist22)
     $          -olap22**2*2*ddist2(3,5)/(dist2**3))
              izlploc = loc(izlp(jlp))

              derep(1,izlploc) = derep(1,izlploc) + forcetemp(1)
              derep(2,izlploc) = derep(2,izlploc) + forcetemp(2)
              derep(3,izlploc) = derep(3,izlploc) + forcetemp(3)
c
              erep = erep + cvrep21*e21 + cvrep22*e22 
c
            end if
            nerep = nerep + 1
         end do
      end do  
c
c     lone pair / lone pair repulsion      
c
      do ilpnl = 1, nlplocnl
         ilp = lpglobnl(ilpnl)
         if (ilpnl.eq.0) cycle
         p = lpatom(ilp)   
         call rotlp1(ilp,di,dix,diz)
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xp = x(p)
         yp = y(p)
         zp = z(p)
         xip = xp - xi
         yip = yp - yi
         zip = zp - zi
         if (use_bounds)  call image (xip,yip,zip)
         xpi = -xip
         ypi = -yip
         zpi = -zip
         rip2 = xip*xip + yip*yip + zip*zip
         rip = sqrt(rip2)
         cilp = lpcharge(ilp)
         cip = rpole(1,p)
         do j = 1, nlplplst(ilpnl)
            jlp = lplplst(j,ilpnl)
c
c  correspondance between indexes and number (1,2,3,4,5)
c  1 -> p, 2-> k, 3 -> ixlp(i), 4 -> izlp(i), 5 -> ixlp(j), 6 -> izlp(j) 
c
            k = lpatom(jlp)
            call rotlp1(jlp,di2,dix2,diz2)
            dcos31 = 0d0
            dcos32 = 0d0
            dolaptot = 0d0
            ddist2 = 0d0
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
               xkj = -xjk
               ykj = -yjk
               zkj = -zjk
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
               xkp = -xpk
               ykp = -ypk
               zkp = -zpk
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

               cos31 = (xpi*xpk + ypi*ypk + zpi*zpk)/(rpk*rip)
               scpipk = xpi*xpk + ypi*ypk + zpi*zpk
c
               dcos31(1,1) =  -(xpk+xpi)/(rip*rpk) +
     $       scpipk*xpi/(rpk*rip**3) + scpipk*xpk/(rip*rpk**3)
               dcos31(2,1) =  -(ypk+ypi)/(rip*rpk) +
     $       scpipk*ypi/(rpk*rip**3) + scpipk*ypk/(rip*rpk**3)
               dcos31(3,1) =  -(zpk+zpi)/(rip*rpk) +
     $       scpipk*zpi/(rpk*rip**3) + scpipk*zpk/(rip*rpk**3)
c
               dcos31(1,2) =  xpi/(rip*rpk) -(scpipk*xpk)/(rip*rpk**3)
               dcos31(2,2) =  ypi/(rip*rpk) -(scpipk*ypk)/(rip*rpk**3)
               dcos31(3,2) =  zpi/(rip*rpk) -(scpipk*zpk)/(rip*rpk**3)
c
c         contribution of the rotation matrixes
c
               dcoslp1 = xpk/(rpk*rip) - xpi*scpipk/(rpk*rip**3) 
               dcoslp2 = ypk/(rpk*rip) - ypi*scpipk/(rpk*rip**3) 
               dcoslp3 = zpk/(rpk*rip) - zpi*scpipk/(rpk*rip**3) 
c
               dcos31(1,3) =  
     $    dix(1,1)*dcoslp1 + dix(1,2)*dcoslp2 + dix(1,3)*dcoslp3
               dcos31(2,3) =  
     $    dix(2,1)*dcoslp1 + dix(2,2)*dcoslp2 + dix(2,3)*dcoslp3
               dcos31(3,3) =  
     $    dix(3,1)*dcoslp1 + dix(3,2)*dcoslp2 + dix(3,3)*dcoslp3
c
               dcos31(1,4) =  
     $    diz(1,1)*dcoslp1 + diz(1,2)*dcoslp2 + diz(1,3)*dcoslp3
               dcos31(2,4) =  
     $    diz(2,1)*dcoslp1 + diz(2,2)*dcoslp2 + diz(2,3)*dcoslp3
               dcos31(3,4) =  
     $    diz(3,1)*dcoslp1 + diz(3,2)*dcoslp2 + diz(3,3)*dcoslp3
c
               dcos31(1,1) = dcos31(1,1) + 
     $    di(1,1)*dcoslp1 + di(1,2)*dcoslp2 + di(1,3)*dcoslp3
               dcos31(2,1) = dcos31(2,1) + 
     $    di(2,1)*dcoslp1 + di(2,2)*dcoslp2 + di(2,3)*dcoslp3
               dcos31(3,1) = dcos31(3,1) + 
     $    di(3,1)*dcoslp1 + di(3,2)*dcoslp2 + di(3,3)*dcoslp3
c
               cos32 = (xkj*xkp + ykj*ykp + zkj*zkp)/(rjk*rpk) 
               sckjkp = xkj*xkp + ykj*ykp + zkj*zkp
c
               dcos32(1,2) =  - (xkp+xkj)/(rjk*rpk) +
     $       sckjkp*xkj/(rpk*rjk**3) + sckjkp*xkp/(rjk*rpk**3)
               dcos32(2,2) =  - (ykp+ykj)/(rjk*rpk) +
     $       sckjkp*ykj/(rpk*rjk**3) + sckjkp*ykp/(rjk*rpk**3)
               dcos32(3,2) =  - (zkp+zkj)/(rjk*rpk) +
     $       sckjkp*zkj/(rpk*rjk**3) + sckjkp*zkp/(rjk*rpk**3)
c
               dcos32(1,1) =  xkj/(rjk*rpk) -(sckjkp*xkp)/(rjk*rpk**3)
               dcos32(2,1) =  ykj/(rjk*rpk) -(sckjkp*ykp)/(rjk*rpk**3)
               dcos32(3,1) =  zkj/(rjk*rpk) -(sckjkp*zkp)/(rjk*rpk**3)
c
c
c         contribution of the rotation matrixes
c
               dcoslp1 = xkp/(rpk*rjk) - xkj*sckjkp/(rpk*rjk**3) 
               dcoslp2 = ykp/(rpk*rjk) - ykj*sckjkp/(rpk*rjk**3) 
               dcoslp3 = zkp/(rpk*rjk) - zkj*sckjkp/(rpk*rjk**3) 
c
               dcos32(1,5) =  
     $    dix2(1,1)*dcoslp1 + dix2(1,2)*dcoslp2 + dix2(1,3)*dcoslp3
               dcos32(2,5) =  
     $    dix2(2,1)*dcoslp1 + dix2(2,2)*dcoslp2 + dix2(2,3)*dcoslp3
               dcos32(3,5) =  
     $    dix2(3,1)*dcoslp1 + dix2(3,2)*dcoslp2 + dix2(3,3)*dcoslp3
c
               dcos32(1,6) =  
     $    diz2(1,1)*dcoslp1 + diz2(1,2)*dcoslp2 + diz2(1,3)*dcoslp3
               dcos32(2,6) =  
     $    diz2(2,1)*dcoslp1 + diz2(2,2)*dcoslp2 + diz2(2,3)*dcoslp3
               dcos32(3,6) =  
     $    diz2(3,1)*dcoslp1 + diz2(3,2)*dcoslp2 + diz2(3,3)*dcoslp3
c
               dcos32(1,2) = dcos32(1,2) + 
     $    di2(1,1)*dcoslp1 + di2(1,2)*dcoslp2 + di2(1,3)*dcoslp3
               dcos32(2,2) = dcos32(2,2) + 
     $    di2(2,1)*dcoslp1 + di2(2,2)*dcoslp2 + di2(2,3)*dcoslp3
               dcos32(3,2) = dcos32(3,2) + 
     $    di2(3,1)*dcoslp1 + di2(3,2)*dcoslp2 + di2(3,3)*dcoslp3

               ss31 = csi*csj
               sp31 = csi*cpj*spij*cos32
               ps31 = cpi*csj*psij*cos31
               pp31 = cpi*cpj*2*cos31*cos32

               oni = cilp!2.0
               onj = cjlp!2.0
               olap3 = (ss31 + sp31 + ps31 + pp31)*sisj
               olap32 = (ss31 + sp31 + ps31 + pp31)*sisj2
c
               dolaptot(1,1) =  
     $         sisj*(ss31+sp31+ps31+pp31)*
     $         (alpha/(4*sqrt(vdwreplpi*vdwreplpj))*xpk/rpk)
     &       + sisj*csi*cpj*spij*dcos32(1,1) 
     &       + sisj*cpi*csj*psij*dcos31(1,1) 
     &       + sisj*cpi*cpj*2*dcos32(1,1)*cos31
     &       + sisj*cpi*cpj*2*cos32*dcos31(1,1)
               dolaptot(2,1) =  
     $         sisj*(ss31+sp31+ps31+pp31)*
     $         (alpha/(4*sqrt(vdwreplpi*vdwreplpj))*ypk/rpk)
     &       + sisj*csi*cpj*spij*dcos32(2,1) 
     &       + sisj*cpi*csj*psij*dcos31(2,1) 
     &       + sisj*cpi*cpj*2*dcos32(2,1)*cos31
     &       + sisj*cpi*cpj*2*cos32*dcos31(2,1)
               dolaptot(3,1) =  
     $         sisj*(ss31+sp31+ps31+pp31)*
     $         (alpha/(4*sqrt(vdwreplpi*vdwreplpj))*zpk/rpk)
     &       + sisj*csi*cpj*spij*dcos32(3,1) 
     &       + sisj*cpi*csj*psij*dcos31(3,1) 
     &       + sisj*cpi*cpj*2*dcos32(3,1)*cos31
     &       + sisj*cpi*cpj*2*cos32*dcos31(3,1)
c
               dolaptot2(1,1) =  
     $         sisj2*(ss31+sp31+ps31+pp31)*
     $         (alpha2/(4*sqrt(vdwreplpi*vdwreplpj))*xpk/rpk)
     &       + sisj2*csi*cpj*spij*dcos32(1,1) 
     &       + sisj2*cpi*csj*psij*dcos31(1,1) 
     &       + sisj2*cpi*cpj*2*dcos32(1,1)*cos31
     &       + sisj2*cpi*cpj*2*cos32*dcos31(1,1)
               dolaptot2(2,1) =  
     $         sisj2*(ss31+sp31+ps31+pp31)*
     $         (alpha2/(4*sqrt(vdwreplpi*vdwreplpj))*ypk/rpk)
     &       + sisj2*csi*cpj*spij*dcos32(2,1) 
     &       + sisj2*cpi*csj*psij*dcos31(2,1) 
     &       + sisj2*cpi*cpj*2*dcos32(2,1)*cos31
     &       + sisj2*cpi*cpj*2*cos32*dcos31(2,1)
               dolaptot2(3,1) =  
     $          sisj2*(ss31+sp31+ps31+pp31)*
     $         (alpha2/(4*sqrt(vdwreplpi*vdwreplpj))*zpk/rpk)
     &       + sisj2*csi*cpj*spij*dcos32(3,1) 
     &       + sisj2*cpi*csj*psij*dcos31(3,1) 
     &       + sisj2*cpi*cpj*2*dcos32(3,1)*cos31
     &       + sisj2*cpi*cpj*2*cos32*dcos31(3,1)
c
               dolaptot(1,2) =  
     $       - sisj*(ss31+sp31+ps31+pp31)*
     $         (alpha/(4*sqrt(vdwreplpi*vdwreplpj))*xpk/rpk)
     &       + sisj*csi*cpj*spij*dcos32(1,2) 
     &       + sisj*cpi*csj*psij*dcos31(1,2) 
     &       + sisj*cpi*cpj*2*dcos32(1,2)*cos31
     &       + sisj*cpi*cpj*2*cos32*dcos31(1,2)
               dolaptot(2,2) =  
     $       - sisj*(ss31+sp31+ps31+pp31)*
     $         (alpha/(4*sqrt(vdwreplpi*vdwreplpj))*ypk/rpk)
     &       + sisj*csi*cpj*spij*dcos32(2,2) 
     &       + sisj*cpi*csj*psij*dcos31(2,2) 
     &       + sisj*cpi*cpj*2*dcos32(2,2)*cos31
     &       + sisj*cpi*cpj*2*cos32*dcos31(2,2)
               dolaptot(3,2) =  
     $       - sisj*(ss31+sp31+ps31+pp31)*
     $         (alpha/(4*sqrt(vdwreplpi*vdwreplpj))*zpk/rpk)
     &       + sisj*csi*cpj*spij*dcos32(3,2) 
     &       + sisj*cpi*csj*psij*dcos31(3,2) 
     &       + sisj*cpi*cpj*2*dcos32(3,2)*cos31
     &       + sisj*cpi*cpj*2*cos32*dcos31(3,2)
c
               dolaptot2(1,2) =  
     $       - sisj2*(ss31+sp31+ps31+pp31)*
     $         (alpha2/(4*sqrt(vdwreplpi*vdwreplpj))*xpk/rpk)
     &       + sisj2*csi*cpj*spij*dcos32(1,2) 
     &       + sisj2*cpi*csj*psij*dcos31(1,2) 
     &       + sisj2*cpi*cpj*2*dcos32(1,2)*cos31
     &       + sisj2*cpi*cpj*2*cos32*dcos31(1,2)
               dolaptot2(2,2) =  
     $       - sisj2*(ss31+sp31+ps31+pp31)*
     $         (alpha2/(4*sqrt(vdwreplpi*vdwreplpj))*ypk/rpk)
     &       + sisj2*csi*cpj*spij*dcos32(2,2) 
     &       + sisj2*cpi*csj*psij*dcos31(2,2) 
     &       + sisj2*cpi*cpj*2*dcos32(2,2)*cos31
     &       + sisj2*cpi*cpj*2*cos32*dcos31(2,2)
               dolaptot2(3,2) =  
     $       - sisj2*(ss31+sp31+ps31+pp31)*
     $         (alpha2/(4*sqrt(vdwreplpi*vdwreplpj))*zpk/rpk)
     &       + sisj2*csi*cpj*spij*dcos32(3,2) 
     &       + sisj2*cpi*csj*psij*dcos31(3,2) 
     &       + sisj2*cpi*cpj*2*dcos32(3,2)*cos31
     &       + sisj2*cpi*cpj*2*cos32*dcos31(3,2)
c
               dolaptot(1,5) = 
     &         sisj*csi*cpj*spij*dcos32(1,5) 
     &       + sisj*cpi*cpj*2*dcos32(1,5)*cos31
               dolaptot(2,5) =  
     &         sisj*csi*cpj*spij*dcos32(2,5) 
     &       + sisj*cpi*cpj*2*dcos32(2,5)*cos31
               dolaptot(3,5) =  
     &         sisj*csi*cpj*spij*dcos32(3,5) 
     &       + sisj*cpi*cpj*2*dcos32(3,5)*cos31
c
               dolaptot2(1,5) = 
     &       + sisj2*csi*cpj*spij*dcos32(1,5) 
     &       + sisj2*cpi*cpj*2*dcos32(1,5)*cos31
               dolaptot2(2,5) =  
     &       + sisj2*csi*cpj*spij*dcos32(2,5) 
     &       + sisj2*cpi*cpj*2*dcos32(2,5)*cos31
               dolaptot2(3,5) =  
     &       + sisj2*csi*cpj*spij*dcos32(3,5) 
     &       + sisj2*cpi*cpj*2*dcos32(3,5)*cos31
c
               dolaptot(1,6) =  
     &       + sisj*csi*cpj*spij*dcos32(1,6) 
     &       + sisj*cpi*cpj*2*dcos32(1,6)*cos31
               dolaptot(2,6) =  
     &       + sisj*csi*cpj*spij*dcos32(2,6) 
     &       + sisj*cpi*cpj*2*dcos32(2,6)*cos31
               dolaptot(3,6) =  
     &       + sisj*csi*cpj*spij*dcos32(3,6) 
     &       + sisj*cpi*cpj*2*dcos32(3,6)*cos31
c
               dolaptot2(1,6) =  
     &       + sisj2*csi*cpj*spij*dcos32(1,6) 
     &       + sisj2*cpi*cpj*2*dcos32(1,6)*cos31
               dolaptot2(2,6) =  
     &       + sisj2*csi*cpj*spij*dcos32(2,6) 
     &       + sisj2*cpi*cpj*2*dcos32(2,6)*cos31
               dolaptot2(3,6) =  
     &       + sisj2*csi*cpj*spij*dcos32(3,6) 
     &       + sisj2*cpi*cpj*2*dcos32(3,6)*cos31
c
               dolaptot(1,3) =  
     &       + sisj*cpi*csj*psij*dcos31(1,3) 
     &       + sisj*cpi*cpj*2*cos32*dcos31(1,3)
               dolaptot(2,3) =  
     &       + sisj*cpi*csj*psij*dcos31(2,3) 
     &       + sisj*cpi*cpj*2*cos32*dcos31(2,3)
               dolaptot(3,3) =  
     &       + sisj*cpi*csj*psij*dcos31(3,3) 
     &       + sisj*cpi*cpj*2*cos32*dcos31(3,3)
c
               dolaptot2(1,3) =  
     &       + sisj2*cpi*csj*psij*dcos31(1,3) 
     &       + sisj2*cpi*cpj*2*cos32*dcos31(1,3)
               dolaptot2(2,3) =  
     &       + sisj2*cpi*csj*psij*dcos31(2,3) 
     &       + sisj2*cpi*cpj*2*cos32*dcos31(2,3)
               dolaptot2(3,3) =  
     &       + sisj2*cpi*csj*psij*dcos31(3,3) 
     &       + sisj2*cpi*cpj*2*cos32*dcos31(3,3)
c
               dolaptot(1,4) =  
     &       + sisj*cpi*csj*psij*dcos31(1,4) 
     &       + sisj*cpi*cpj*2*cos32*dcos31(1,4)
               dolaptot(2,4) =  
     &       + sisj*cpi*csj*psij*dcos31(2,4) 
     &       + sisj*cpi*cpj*2*cos32*dcos31(2,4)
               dolaptot(3,4) =  
     &       + sisj*cpi*csj*psij*dcos31(3,4) 
     &       + sisj*cpi*cpj*2*cos32*dcos31(3,4)
c
               dolaptot2(1,4) =  
     &       + sisj2*cpi*csj*psij*dcos31(1,4) 
     &       + sisj2*cpi*cpj*2*cos32*dcos31(1,4)
               dolaptot2(2,4) =  
     &       + sisj2*cpi*csj*psij*dcos31(2,4) 
     &       + sisj2*cpi*cpj*2*cos32*dcos31(2,4)
               dolaptot2(3,4) =  
     &       + sisj2*cpi*csj*psij*dcos31(3,4) 
     &       + sisj2*cpi*cpj*2*cos32*dcos31(3,4)
c
               dlp1 = -xij/rij
               dlp2 = -yij/rij
               dlp3 = -zij/rij
c
               ddist2(1,1) =  di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
               ddist2(2,1) =  di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
               ddist2(3,1) =  di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3
c
               ddist2(1,3) =  
     $          dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
               ddist2(2,3) =  
     $          dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
               ddist2(3,3) =  
     $          dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3
c
               ddist2(1,4) =  
     $          diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
               ddist2(2,4) =  
     $          diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
               ddist2(3,4) =  
     $          diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3
c
               ddist2(1,2) = 
     $          -di2(1,1)*dlp1 - di2(1,2)*dlp2 - di2(1,3)*dlp3
               ddist2(2,2) =  
     $          -di2(2,1)*dlp1 - di2(2,2)*dlp2 - di2(2,3)*dlp3
               ddist2(3,2) =  
     $          -di2(3,1)*dlp1 - di2(3,2)*dlp2 - di2(3,3)*dlp3
c
               ddist2(1,5) =  
     $          -dix2(1,1)*dlp1 - dix2(1,2)*dlp2 - dix2(1,3)*dlp3
               ddist2(2,5) =  
     $          -dix2(2,1)*dlp1 - dix2(2,2)*dlp2 - dix2(2,3)*dlp3
               ddist2(3,5) =  
     $          -dix2(3,1)*dlp1 - dix2(3,2)*dlp2 - dix2(3,3)*dlp3
c
               ddist2(1,6) =  
     $          -diz2(1,1)*dlp1 - diz2(1,2)*dlp2 - diz2(1,3)*dlp3
               ddist2(2,6) =  
     $          -diz2(2,1)*dlp1 - diz2(2,2)*dlp2 - diz2(2,3)*dlp3
               ddist2(3,6) =  
     $          -diz2(3,1)*dlp1 - diz2(3,2)*dlp2 - diz2(3,3)*dlp3
c
               forcetemp(1) = cvrep31*oni*onj*
     $           (2*dolaptot(1,1)*olap3*(1/rij)
     $           -olap3**2*(ddist2(1,1)/rij2)) 
               forcetemp(2) = cvrep31*oni*onj*
     $           (2*dolaptot(2,1)*olap3*(1/rij)
     $           -olap3**2*(ddist2(2,1)/rij2)) 
               forcetemp(3) = cvrep31*oni*onj*
     $           (2*dolaptot(3,1)*olap3*(1/rij)
     $           -olap3**2*(ddist2(3,1)/rij2)) 
               forcetemp(1) = forcetemp(1) + cvrep32*oni*onj*
     $           (2*dolaptot2(1,1)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(1,1)/(rij**3))
               forcetemp(2) = forcetemp(2) + cvrep32*oni*onj*
     $           (2*dolaptot2(2,1)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(2,1)/(rij**3))
               forcetemp(3) = forcetemp(3) + cvrep32*oni*onj*
     $           (2*dolaptot2(3,1)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(3,1)/(rij**3))
               ploc = loc(p)
               derep(1,ploc) = derep(1,ploc) + forcetemp(1) 
               derep(2,ploc) = derep(2,ploc) + forcetemp(2) 
               derep(3,ploc) = derep(3,ploc) + forcetemp(3) 
c
               forcetemp(1) = cvrep31*oni*onj*
     $           (2*dolaptot(1,3)*olap3*(1/rij)
     $           -olap3**2*(ddist2(1,3)/rij2)) 
               forcetemp(2) =  cvrep31*oni*onj*
     $           (2*dolaptot(2,3)*olap3*(1/rij)
     $           -olap3**2*(ddist2(2,3)/rij2)) 
               forcetemp(3) =  cvrep31*oni*onj*
     $           (2*dolaptot(3,3)*olap3*(1/rij)
     $           -olap3**2*(ddist2(3,3)/rij2)) 
               forcetemp(1) = forcetemp(1) + cvrep32*oni*onj*
     $           (2*dolaptot2(1,3)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(1,3)/(rij**3))
               forcetemp(2) = forcetemp(2) + cvrep32*oni*onj*
     $           (2*dolaptot2(2,3)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(2,3)/(rij**3))
               forcetemp(3) = forcetemp(3) + cvrep32*oni*onj*
     $           (2*dolaptot2(3,3)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(3,3)/(rij**3))
               ixlploc = loc(ixlp(ilp))
               derep(1,ixlploc) = derep(1,ixlploc) + forcetemp(1) 
               derep(2,ixlploc) = derep(2,ixlploc) + forcetemp(2) 
               derep(3,ixlploc) = derep(3,ixlploc) + forcetemp(3) 
c
               forcetemp(1) = cvrep31*oni*onj*
     $           (2*dolaptot(1,4)*olap3*(1/rij)
     $           -olap3**2*(ddist2(1,4)/rij2)) 
               forcetemp(2) = cvrep31*oni*onj*
     $           (2*dolaptot(2,4)*olap3*(1/rij)
     $           -olap3**2*(ddist2(2,4)/rij2)) 
               forcetemp(3) = cvrep31*oni*onj*
     $           (2*dolaptot(3,4)*olap3*(1/rij)
     $           -olap3**2*(ddist2(3,4)/rij2)) 
               forcetemp(1) = forcetemp(1) + cvrep32*oni*onj*
     $           (2*dolaptot2(1,4)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(1,4)/(rij**3))
               forcetemp(2) = forcetemp(2) + cvrep32*oni*onj*
     $           (2*dolaptot2(2,4)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(2,4)/(rij**3))
               forcetemp(3) = forcetemp(3) + cvrep32*oni*onj*
     $           (2*dolaptot2(3,4)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(3,4)/(rij**3))
               izlploc = loc(izlp(ilp))
               derep(1,izlploc) = derep(1,izlploc) + forcetemp(1) 
               derep(2,izlploc) = derep(2,izlploc) + forcetemp(2) 
               derep(3,izlploc) = derep(3,izlploc) + forcetemp(3) 
c
               forcetemp(1) = cvrep31*oni*onj*
     $           (2*dolaptot(1,2)*olap3*(1/rij)
     $           -olap3**2*(ddist2(1,2)/rij2)) 
               forcetemp(2) = cvrep31*oni*onj*
     $           (2*dolaptot(2,2)*olap3*(1/rij)
     $           -olap3**2*(ddist2(2,2)/rij2)) 
               forcetemp(3) = cvrep31*oni*onj*
     $           (2*dolaptot(3,2)*olap3*(1/rij)
     $           -olap3**2*(ddist2(3,2)/rij2))
               forcetemp(1) = forcetemp(1) + cvrep32*oni*onj*
     $           (2*dolaptot2(1,2)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(1,2)/(rij**3))
               forcetemp(2) = forcetemp(2) + cvrep32*oni*onj*
     $           (2*dolaptot2(2,2)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(2,2)/(rij**3))
               forcetemp(3) = forcetemp(3) + cvrep32*oni*onj*
     $           (2*dolaptot2(3,2)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(3,2)/(rij**3))
               kloc = loc(k)
               derep(1,kloc) = derep(1,kloc) + forcetemp(1) 
               derep(2,kloc) = derep(2,kloc) + forcetemp(2) 
               derep(3,kloc) = derep(3,kloc) + forcetemp(3) 
c
               forcetemp(1) = cvrep31*oni*onj*
     $           (2*dolaptot(1,5)*olap3*(1/rij)
     $           -olap3**2*(ddist2(1,5)/rij2)) 
               forcetemp(2) = cvrep31*oni*onj*
     $           (2*dolaptot(2,5)*olap3*(1/rij)
     $           -olap3**2*(ddist2(2,5)/rij2)) 
               forcetemp(3) = cvrep31*oni*onj*
     $           (2*dolaptot(3,5)*olap3*(1/rij)
     $           -olap3**2*(ddist2(3,5)/rij2)) 
               forcetemp(1) = forcetemp(1) + cvrep32*oni*onj*
     $           (2*dolaptot2(1,5)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(1,5)/(rij**3))
               forcetemp(2) = forcetemp(2) + cvrep32*oni*onj*
     $           (2*dolaptot2(2,5)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(2,5)/(rij**3))
               forcetemp(3) = forcetemp(3) + cvrep32*oni*onj*
     $           (2*dolaptot2(3,5)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(3,5)/(rij**3))
               ixlploc = loc(ixlp(jlp))
               derep(1,ixlploc) = derep(1,ixlploc) + forcetemp(1) 
               derep(2,ixlploc) = derep(2,ixlploc) + forcetemp(2) 
               derep(3,ixlploc) = derep(3,ixlploc) + forcetemp(3) 
c
               forcetemp(1) = cvrep31*oni*onj*
     $           (2*dolaptot(1,6)*olap3*(1/rij)
     $           -olap3**2*(ddist2(1,6)/rij2)) 
               forcetemp(2) = cvrep31*oni*onj*
     $           (2*dolaptot(2,6)*olap3*(1/rij)
     $           -olap3**2*(ddist2(2,6)/rij2)) 
               forcetemp(3) = cvrep31*oni*onj*
     $           (2*dolaptot(3,6)*olap3*(1/rij)
     $           -olap3**2*(ddist2(3,6)/rij2)) 
               forcetemp(1) = forcetemp(1) + cvrep32*oni*onj*
     $           (2*dolaptot2(1,6)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(1,6)/(rij**3))
               forcetemp(2) = forcetemp(2) + cvrep32*oni*onj*
     $           (2*dolaptot2(2,6)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(2,6)/(rij**3))
               forcetemp(3) = forcetemp(3) + cvrep32*oni*onj*
     $           (2*dolaptot2(3,6)*olap32*(1/rij2)
     $           -olap32**2*2*ddist2(3,6)/(rij**3))
               izlploc = loc(izlp(jlp))
               derep(1,izlploc) = derep(1,izlploc) + forcetemp(1) 
               derep(2,izlploc) = derep(2,izlploc) + forcetemp(2) 
               derep(3,izlploc) = derep(3,izlploc) + forcetemp(3) 
c
               e31 = oni*onj*olap3*olap3/rij
               e32 = oni*onj*olap32*olap32/rij2
               erep = erep + cvrep31*e31 + cvrep32*e32 
               nerep = nerep + 1
            end if 
         end do 
      end do     
c     
c
      deallocate (ddist2)
      deallocate (dolap2bk)
      deallocate (dolap2ak)
      deallocate (dolap2bk2)
      deallocate (dolap2ak2)
      deallocate (dolaptot2)
      deallocate (dolap1ac2)
      deallocate (dolap1ad2)
      deallocate (dolap1bc2)
      deallocate (dolap1bd2)
      deallocate (dolaptot)
      deallocate (dolap1ac)
      deallocate (dolap1ad)
      deallocate (dolap1bc)
      deallocate (dolap1bd)
      deallocate (dcos1)
      deallocate (dcos2)
      deallocate (dcos3)
      deallocate (dcos4)
      deallocate (dcos5)
      deallocate (dcos6)
      deallocate (dcos7)
      deallocate (dcos8)
c
      deallocate (dcos21)
      deallocate (dcos22)
      deallocate (dcos23)
      deallocate (dcos24)
      deallocate (dcos31)
      deallocate (dcos32)
      return
      end
