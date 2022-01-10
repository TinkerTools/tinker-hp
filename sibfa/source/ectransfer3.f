c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine ectransfer3  --  charge transfer energy & analysis  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "ectransfer3" calculates the charge transfer energy
c     and partitions the energy among the atoms
c
      subroutine ectransfer3
      implicit none
c
      call ectransfer3b
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ectransfer3b  --  charge transfer analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ectransfer3b" calculates the charge transfer energy
c     and partitions the energy among the atoms using neighbor lists
c
c
      subroutine ectransfer3b
      use action
      use analyz
      use atmtyp
      use atmlst
      use atoms
      use bound
      use cell
      use charge
      use chargetransfer
      use cutoff
      use domdec
      use energi
      use group
      use inform
      use inter
      use iounit
      use kct
      use molcul
      use mpole
      use neigh
      use potent
      use usage
      implicit none
      integer i,ii,kk,k1,k2,kpole,ilp,iacc,kaccept
      integer ilpnl,jj,knl,iacc1
      real*8 xlp,ylp,zlp,clp
      real*8 xa,ya,za
      real*8 xralp,yralp,zralp,ralp2,ralp
      real*8 xk1,yk1,zk1
      real*8 xk2,yk2,zk2
      real*8 xrk1k2,yrk1k2,zrk1k2,rk1k22,rk1k2
      real*8 xrak1,yrak1,zrak1,rak12,rak1
      real*8 xrak2,yrak2,zrak2,rak22,rak2
      real*8 rhoak1,rhoak2
      real*8 Vk2k1,rk1k23,rk1K25
      real*8 Vk1k2,ksih
      real*8 V,Vbis,Vbis1,Vbis2,rak13,rak23,rak15,rak25
      real*8 Vtemp,Vconv
      real*8 rvdwi,rvdwk1,rvdwk2
      real*8 cos1,cos2
      real*8 task2,mak2,tapk2,task1,mak1,tapk1
      real*8 e1,e2,e,deltae
      real*8 Cs,Cp
      real*8 conv,prop
      real*8 dkr,qkx,qky,qkr,scd,scq,fc,fd1,fd2,fq1,fq2
      real*8 vdwinc,ef0(3)
      real*8 ectloc
      real*8, allocatable :: Vtot1(:),Vtot2(:,:)
c
      allocate (Vtot1(nlp)) 
      allocate (Vtot2(2,nacceptbloc)) 
      Vtot1 = 0d0
      Vtot2 = 0d0
c
c     zero out the charge transfer interaction energy and partitioning
c
      V = 0.0d0
      Vtot1 = 0d0
      Vtot2 = 0d0
      nect = 0
      ect = 0.0d0
      e = 0.0d0
      aect = 0d0
      if (nlp .eq. 0)  return
      call rotlp
c
c     compute the potential of the whole system on the acceptors (except the local molecule)
c
      if (use_ctpot) then
        do iacc = 1, nacceptlocnl
          i = acceptglobnl(iacc)
          k1 = acceptor(1,i)
          k2 = acceptor(2,i)
          iacc1 = acceptloc(i)
          call potmole3a(k1,molcule(k1),0,Vtot2(1,iacc1),0,
     $      accpotlst(1,1,iacc),naccpotlst(1,iacc))
          call potmole3a(k2,molcule(k2),0,Vtot2(2,iacc1),0,
     $      accpotlst(1,2,iacc),naccpotlst(2,iacc))
        end do
c
c     send and receive potentials to/from the neighboring processes
c
        call commpotacc(Vtot2)
c
c     compute the potential of the whole system on the lp carriers (except the local molecule)
c
        do ii = 1, nlplocnl
          ilp = lpglobnl(ii)
          i = lpatom(ilp) 
          call potmole3a(i,molcule(i),0,Vtot1(ilp),0,lppotlst(1,ii),
     $      nlppotlst(ii))
        end do
c
c     send and receive potentials to/from the neighboring processes
c
        call commpotlp(Vtot1)
      end if
c
c     compute and partition the charge transfer interaction energy
c
c     Loop over the lone pairs
c
      do ii = 1, nlplocnl
         ilp = lpglobnl(ii)
         ectloc = 0.0d0
         i = lpatom(ilp)
         clp = lpcharge(ilp)
c
c     Get geometrical data about the LP and its carrier         
c
         xlp = rlonepair(1,ilp)
         ylp = rlonepair(2,ilp)
         zlp = rlonepair(3,ilp)
         xa = x(i)
         ya = y(i)
         za = z(i)
         xralp = xlp - xa
         yralp = ylp - ya
         zralp = zlp - za
c         if (use_bounds)  call image (xralp,yralp,zralp)
         ralp2 = xralp*xralp + yralp*yralp + zralp*zralp
         ralp  = sqrt(ralp2)
         V = Vtot1(ilp)
c
c     Loop over the acceptors
c
         do jj = 1, nlpacclst(ii)
            kk = lpacclst(jj,ii)
            kaccept = acceptloc(kk) 
            k1 = acceptor(1,kk)
            k2 = acceptor(2,kk)
            if (molcule(k1).ne.molcule(i)) then
c
c     Get geometrical data about the accepting(s) atom(s) 
c
              xk1 = x(k1)
              yk1 = y(k1)
              zk1 = z(k1)
              if (acceptor(2,kk).ne.0) then
                xk2 = x(k2)
                yk2 = y(k2)
                zk2 = z(k2)
                xrk1k2 = xk2 - xk1
                yrk1k2 = yk2 - yk1
                zrk1k2 = zk2 - zk1
                if (use_bounds)  call image (xrk1k2,yrk1k2,zrk1k2)
                rk1k22 = xrk1k2*xrk1k2 + yrk1k2*yrk1k2 + zrk1k2*zrk1k2
                rk1k2 = sqrt(rk1k22)
              end if
c
c     Get the interatomic distances rak1 and rak2
c
              xrak1 = xk1 - xa
              yrak1 = yk1 - ya
              zrak1 = zk1 - za
              if (use_bounds)  call image (xrak1,yrak1,zrak1)
              rak12 = xrak1*xrak1 + yrak1*yrak1 + zrak1*zrak1
c
              if (rak12.gt.ctransfercut*ctransfercut) cycle
c
c               
              rak1 = sqrt(rak12)
              if (acceptor(2,kk).ne.0) then
                xrak2 = xk2 - xa
                yrak2 = yk2 - ya
                zrak2 = zk2 - za
                if (use_bounds)  call image (xrak2,yrak2,zrak2)
                rak22 = xrak2*xrak2+ yrak2*yrak2 + zrak2*zrak2
                rak2 = sqrt(rak22)
              end if
c              
c    vdw radius is vdw of the carrier plus increment of the lp
c
              rvdwi = rvdwct1(i) + dincr_lpect(ilp)
c
c     Compute the interaction between acceptor molecule and the ii (lone pair) orbital on i
c
c
              cos1=(xralp*xrak1 + yralp*yrak1 + zralp*zrak1)/(ralp*rak1)
c
c
c     Compute the exponent term of the interaction
c
              if (acceptor(2,kk).ne.0) then
                rvdwk2 = rvdwct2(k2)
                rhoak2 = rak2/(4*sqrt(rvdwi*rvdwk2))
c
c     Compute the potential of the acceptor molecule on k2 (H) virtual orbital
c      (Approximation using Slater exponent ksih)
c
                ksih = 1.20d0
                kpole = ipole(k2)
                Vk2k1 = 0.375 * ksih + 0.625 * ksih * rpole(1,kpole)
                Vk2k1 = 1.88976*331.93359*Vk2k1
c
c     Potential of the whole system but the lp molecule on k2
c
                Vtemp = 0.0d0
                Vbis1 = Vtot2(2,kaccept)
                if (use_ctpot) call potmole3aloc(k2,molcule(i),Vtemp)
                Vbis1 = Vbis1 - Vtemp
                Vk2k1 = Vk2k1 + Vbis1
              end if
c
c     Compute the potential of the acceptor molecule on k1 (B) virtual orbital
c
              Vk1k2 = 0.0d0
              Vtemp = 0.0d0
              Vbis2 = Vtot2(1,kaccept)
              if (acceptor(1,kk).ne.0) then
                call potmoleint3a(k1,Vk1k2)
                if (use_ctpot) call potmole3aloc(k1,molcule(i),Vtemp)
                Vbis2 = Vbis2 - Vtemp
                Vk1k2 = Vk1k2 + Vbis2
              else
c
c     empiric values for cations, Zn
c                 
                Vk1k2 = 1497.102d0
                if (use_ctpot) call potmole3aloc(k1,molcule(i),Vtemp)
                Vbis2 = Vbis2 - Vtemp
              end if
              rvdwk1 = rvdwct2(k1)
              rhoak1 = rak1/(4*sqrt(rvdwi*rvdwk1))
c
c     Final assembling of all the terms of the numerator
c
c    condition à verifier pour identifier le couple accepteur
c
              task1 = tas(atomic(k1),atomic(i))
              tapk1 = tap(atomic(k1),atomic(i))
              mak1 = ma(atomic(k1),atomic(i))
              if (acceptor(2,kk).ne.0) then
                cos2=(xralp*xrak2 + yralp*yrak2+zralp*zrak2)/(ralp*rak2)
                task2 = tas(atomic(k2),atomic(i))
                tapk2 = tap(atomic(k2),atomic(i))
                mak2 = ma(atomic(k2),atomic(i))
              end if
c
              Cs = hybrid_lp(1,ilp)
              Cp = hybrid_lp(2,ilp)
c
              e1 = 0.0d0
              e2 = exp(-etaect*rhoak2)*(D2*Vk2k1-V)*(Cs*task2+
     $              Cp*mak2*tapk2*cos2)
c
c     Compute the denominator 
c
              Vbis1 = 0.0d0
              Vbis2 = 0.0d0
              if (acceptor(2,kk).eq.0) then
                Vbis1 = Vtot2(1,kaccept)
              else if (acceptor(2,kk).ne.0) then
                e1 = exp(-etaect*rhoak1)*(D1*Vk1k2-V)*(Cs*task1+
     $             Cp*mak1*tapk1*cos1)       
                Vbis1 = Vtot2(1,kaccept)
                Vbis2 = Vtot2(2,kaccept)
              end if
c
              Vbis1 = Vbis1/23.06
              Vbis2 = Vbis2/23.06
              Vconv = V/23.06
c
              if (acceptor(2,kk).ne.0) then
                deltae = ahct(i) + Vconv  - (aelct(k1) + Vbis2) 
              else
c
c  for the zinc
c
                deltae = ahct(k1) + Vconv - (17.99d0 + Vbis2) 
              end if
              e = -clp*(e1-e2)**2/deltae
c
              conv = 27.2/627.5442
              e  = e*conv
c
c     Facteur empirique, si cation : different
c
              if (acceptor(2,kk).ne.0) then
c                prop = 0.7712d0 
                prop = 0.6364d0 
              else
                prop = 0.61d0
              end if
              e = prop*e
              ect = ect + e
c            
c     increment the overall charge transfer energy component
c
              nect = nect + 1
              aect(loc(i)) = aect(loc(i)) + 0.5d0*e
              aect(loc(k1)) = aect(loc(k1)) + 0.25d0*e
              if (acceptor(2,kk).ne.0) then
                aect(loc(k2)) = aect(loc(k2)) + 0.25d0*e
              end if
c
c     increment the total intermolecular energy
c
              if (molcule(i) .ne. molcule(k1)) then
                einter = einter + e
              end if
            end if     
         end do
       end do
       return
       end
c
c     "potmole3a" calculates the electrostatic potential created by the whole
c     system except some molecule 
c
c     - rule = 0 : compute the whole potential on ia except the one created by imol1
c     - rule = 1 : compute the whole potential on ia except the one created by imol1 and imol2 
c
c
      subroutine potmole3a(ia,imol1,imol2,V,rule,list,nlist)
      use sizes
      use action
      use analyz
      use atmtyp
      use atmlst
      use atoms
      use bound
      use charge
      use chargetransfer
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use iounit
      use molcul
      use mplpot
      use mpole
      use neigh
      use polpot
      use polgrp
      use shunt
      use usage
      implicit none
      integer rule,imol1,imol2
      integer i,j,k,l,k1
      integer ii,kk,ia,im,ilp
      integer inl,kkk,kglob,kbis,kkpole
      integer kx,ky,kz,natom
      integer iipole,iglob
      integer list(maxelst),nlist
      real*8 e,ei,fgrp,V
      real*8 f,fm,fp,fi,fik
      real*8 xi,yi,zi
      real*8 r,r2,xr,yr,zr
      real*8 xrbis,yrbis,zrbis,rbis,r2bis
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck
      real*8 qix,qiy,qiz
      real*8 dir,fc,fd1,fd2,fq1,fq2,scd,scq
      real*8 scale7
      real*8 alpha1,alpha2,beta1,beta2,etaemtp,rvdw1,rvdw2
      real*8 fd,fq
      real*8 fbis
      real*8 rin,rout,s,ds
      logical usei,usek
      character*6 mode
c
 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')
c
c     check whether the molecules imol1 and imol2 exists
c
      if ((imol1.gt.nmol).or.(imol2.gt.nmol)) then
        write (iout,10)
   10   format (/,' POTMOLA  --  Unable to Find the molecule')
        call fatal
      end if
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      fbis = 331.93359 
c
      mode = 'EWALD'
      call switch (mode)
c
      rin = 1.0d0
      rout = 2.0d0
c      
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
      iipole = pollist(ia)
c
      do kkk = 1, nlist
        kkpole = list(kkk)
        kglob = ipole(kkpole)
        if (((rule.eq.0).and.(molcule(kglob).ne.imol1)).or.
     $    ((rule.eq.1).and.(molcule(kglob).ne.imol1).and.
     $    (molcule(kglob).ne.imol2))) then
          ci = rpole(1,kkpole)
          dix = rpole(2,kkpole)
          diy = rpole(3,kkpole)
          diz = rpole(4,kkpole)
          qixx = rpole(5,kkpole)
          qixy = rpole(6,kkpole)
          qixz = rpole(7,kkpole)
          qiyy = rpole(9,kkpole)
          qiyz = rpole(10,kkpole)
          qizz = rpole(13,kkpole)
          xr = x(ia) - x(kglob)
          yr = y(ia) - y(kglob)
          zr = z(ia) - z(kglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2.gt.(mpolectcut+rout)*(mpolectcut+rout)) cycle
          r = sqrt(r2)
c
c    get switchung function terms
c
          call switch_emtp(.false.,r,mpolectcut-rin,
     $       mpolectcut,s,ds)
          qix = qixx*xr + qixy*yr + qixz*zr
          qiy = qixy*xr + qiyy*yr + qiyz*zr
          qiz = qixz*xr + qiyz*yr + qizz*zr
          scd = dix*xr + diy*yr + diz*zr
          scq = qix*xr + qiy*yr + qiz*zr
c                
c     compute the energy contributions for this interaction
c
           rr1 = 1.0d0 / r
           rr3 = rr1 / r2
           rr5 = 3.0d0 * rr3 / r2
           fc  = ci*rr1
           fd = scd*rr3
           fq = scq*rr5/2.0d0
           e  = fc + fd + fq
c
c     apply the energy adjustments for scaled interactions
c
           fm = f 
           fp = f 
           e = fbis*e
c
c     increment the electrostatic energy of the molecule on the atom     
c
           V = V + s*e
          end if  
 100      continue
      end do
c
      return
      end
c
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine potmoleint3a  --  internal molecular electrostatic potential ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "potmoleint3a" calculates the electrostatic potential generated by a molecule on the central atom of an electron accepting bond (X-H inside the molecule)
c     The approximation made is based on : 
c     International Journal of Quantum CHemistry : Vol XXIX,101-118 (1986)
c
c
      subroutine potmoleint3a(ia,emol)
      use atmtyp
      use couple
      use mpole
      implicit none
      integer ia,i,k,kk
      real*8 emol,nb,nh
      real*8 fh,ksih
c
c     Initialization
c
      emol = 0.0d0
      nb = 0
      nh = 0
c
c     Value of the slater exponent for H
c
      ksih = 1.20d0
c
c     Loop over the atoms connected to ia
c
      do i = 1, n12(ia)
c
c     Check whether it's an hydrogen
c
        k = i12(i,ia)
        kk = ipole(k)
        nb = nb + 1
        if (atomic(k).eq.1) then
          nh = nh +1
          fh = 0.375 * ksih + 0.625 * ksih * rpole(1,kk)
          emol = emol + fh
        end if
      end do
c
      emol = emol / (nb * nh)
      emol = emol*331.93359*1.88976
c      
      return
      end
c
c     subroutine potmole3aloc : computes the electrostatic potential created by molecule
c     imol2 on ia
c
      subroutine potmole3aloc(ia,imol2,V)
      use sizes
      use action
      use analyz
      use atmtyp
      use atoms
      use bound
      use charge
      use chargetransfer
      use chgpot
      use couple
      use cutoff
      use energi
      use iounit
      use molcul
      use mplpot
      use mpole
      use polpot
      use polgrp
      use shunt
      use usage
      implicit none
      integer imol2
      integer i,j,k,l,k1
      integer ii,kk,ia,im,ilp
      integer ix,iy,iz
      integer kx,ky,kz,natom
      integer iglob,iloc
      real*8 e,ei,fgrp,V
      real*8 f,fm,fp,fi,fik
      real*8 xi,yi,zi
      real*8 r,r2,xr,yr,zr
      real*8 xrbis,yrbis,zrbis,rbis,r2bis
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck
      real*8 qix,qiy,qiz
      real*8 dir,fc,fd1,fd2,fq1,fq2,scd,scq
      real*8 scale7
      real*8 alpha1,alpha2,beta1,beta2,etaemtp,rvdw1,rvdw2
      real*8 fd,fq
      real*8 fbis
      real*8 rin,rout,s,ds
      logical usei,usek
      character*6 mode
      
c
c     check whether the molecules imol2 exists
c
      if (imol2.gt.nmol) then
        write (iout,10)
   10   format (/,' POTMOL3Aloc  --  Unable to Find the molecule')
        call fatal
      end if
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      fbis = 331.93359 
c
      mode = 'EWALD'
      call switch (mode)
      rin = 1.0d0
      rout = 2.0d0
c      
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
      do i = imol(1,imol2), imol(2,imol2)
        ii = pollist(i)
        iloc = loc(i)
        iglob = i
        ci = rpole(1,ii)
        dix = rpole(2,ii)
        diy = rpole(3,ii)
        diz = rpole(4,ii)
        qixx = rpole(5,ii)
        qixy = rpole(6,ii)
        qixz = rpole(7,ii)
        qiyy = rpole(9,ii)
        qiyz = rpole(10,ii)
        qizz = rpole(13,ii)
        xr = x(ia) - x(iglob)
        yr = y(ia) - y(iglob)
        zr = z(ia) - z(iglob)
        if (use_bounds)  call image (xr,yr,zr)
        r2 = xr*xr + yr* yr + zr*zr
        r = sqrt(r2)
        if (r2.gt.(mpolectcut+rout)*(mpolectcut+rout)) cycle
c
c    get switchung function terms
c
        call switch_emtp(.false.,r,mpolectcut-rin,mpolectcut,s,
     $   ds)
        qix = qixx*xr + qixy*yr + qixz*zr
        qiy = qixy*xr + qiyy*yr + qiyz*zr
        qiz = qixz*xr + qiyz*yr + qizz*zr
        scd = dix*xr + diy*yr + diz*zr
        scq = qix*xr + qiy*yr + qiz*zr
c                
c     compute the energy contributions for this interaction
c
        rr1 = 1.0d0 / r
        rr3 = rr1 / r2
        rr5 = 3.0d0 * rr3 / r2
        fc  = ci*rr1
        fd = scd*rr3
        fq = scq*rr5/2.0d0
        e  = fc + fd + fq
c
c     apply the energy adjustments for scaled interactions
c
        fm = f 
        fp = f 
        e = fbis*e
c
c     increment the electrostatic energy of the molecule on the atom     
c
        V = V + s*e
      end do
c
      return
      end
