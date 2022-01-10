c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine ectransfer1  --  charge transfer energy & derivatives  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "ectransfer1" calculates the charge transfer energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine ectransfer1
      implicit none
c
      call ectransfer1b
      return
      end
c
c
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine ectransfer1b  --  charge transfer derivs                    ##
c     ##                                                                         ##
c     #############################################################################
c
c
c     "ectransfer1b" calculates the charge transfer energy
c     and first derivatives
c
c
      subroutine ectransfer1b
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
      use inter
      use kct
      use molcul
      use mpole
      use neigh
      use potent
      use usage
      use virial
      implicit none
      integer i,j,ii,kk,k1,k2,kpole,kaccept,ilp,iacc,iglob
      integer ilpnl,jj,knl,iacc1,jglob
      real*8 xlp,ylp,zlp,clp
      real*8 xa,ya,za
      real*8 xralp,yralp,zralp,ralp2,ralp
      real*8 xk1,yk1,zk1
      real*8 xk2,yk2,zk2
      real*8 xrk1k2,yrk1k2,zrk1k2,rk1k22,rk1k2
      real*8 xrak1,yrak1,zrak1,rak12,rak1
      real*8 xrak2,yrak2,zrak2,rak22,rak2
      real*8 scalpk1,scalpk2
      real*8 rhoak1,rhoak2
      real*8 Cs,Cp
      real*8 Vk2k1,rk1k23,rk1K25,V,Vconv
      real*8 Vk1k2,ksih,Vtemp
      real*8 rak13,rak23,rak15,rak25
      real*8 cos1,cos2
      real*8 dcoslp1,dcoslp2,dcoslp3
      real*8 task2,mak2,tapk2,task1,mak1,tapk1
      real*8 e1,e2,e,deltae,ef0(3)
      real*8 rvdwi,rvdwk1,rvdwk2,Vbis1,Vbis2
      real*8 conv,prop
      real*8 di(3,3),dix(3,3),diz(3,3)
      real*8, allocatable :: Vtot1(:),Vtot2(:,:)
      real*8, allocatable :: dVtot1(:,:,:),dVtot1t(:,:,:)
      real*8, allocatable :: dVtot2(:,:,:,:),dVtot2t(:,:,:,:)
      real*8, allocatable :: dV(:,:),dVt(:,:),dcos1(:,:)
      real*8, allocatable :: dVbis1(:,:), dVbis1t(:,:)
      real*8, allocatable :: dVbis2(:,:), dVbis2t(:,:)
      real*8, allocatable :: dcos2(:,:),dnum(:,:)
      real*8, allocatable :: dectt(:,:)
c
      allocate (dectt(3,nbloc))
c
      allocate (Vtot1(nlp))
      allocate (Vtot2(2,nbloc))
c
      allocate (dV(3,nbloc))
      allocate (dVt(3,nbloc))
      allocate (dVtot1(3,nbloc,nlp))
      allocate (dVtot1t(3,nbloc,nlp))
      allocate (dVtot2(3,nbloc,nacceptbloc,2))
      allocate (dVtot2t(3,nbloc,nacceptbloc,2))
      
      allocate (dcos1(3,4))
      allocate (dVbis1(3,nbloc))
      allocate (dVbis1t(3,nbloc))
      allocate (dVbis2(3,nbloc))
      allocate (dVbis2t(3,nbloc))
      allocate (dcos2(3,4))
      allocate (dnum(3,nbloc))
c
      dectt = 0d0
      Vtot1 = 0d0
      Vtot2 = 0d0
      dVtot1 = 0d0
      dVtot1t = 0d0
      dVtot2 = 0d0
      dVtot2t = 0d0
      call aclear(3*nbloc,dVbis1)
      call aclear(3*nbloc,dVbis1t)
      call aclear(3*nbloc,dVbis2)
      call aclear(3*nbloc,dVbis2t)
c
      ect = 0.0d0
      dect = 0d0
      if (nlp .eq. 0)  return
c
      call rotlp
c
c     compute the potential of the whole system on the acceptors (except the local molecule)
c
      if (use_ctpot) then
        do iacc = 1, nacceptlocnl
          i = acceptglobnl(iacc)
          iacc1 = acceptloc(i)
          k1 = acceptor(1,i)
          k2 = acceptor(2,i)
          call potmole1a(k1,molcule(k1),0,Vtot2(1,iacc1),
     $    dVtot2(:,:,iacc1,1),dVtot2t(:,:,iacc1,1),0,
     $    accpotlst(1,1,iacc),naccpotlst(1,iacc))
          call potmole1a(k2,molcule(k2),0,Vtot2(2,iacc1),
     $    dVtot2(:,:,iacc1,2),dVtot2t(:,:,iacc1,2),0,
     $    accpotlst(1,2,iacc),naccpotlst(2,iacc))
        end do
c
c     send and receive potentials and derivatives to/from the neighboring processes
c
        call commpotacc(Vtot2)
        call commpotacc1(dVtot2,dVtot2t)
c
c     compute the potential of the whole system on the lp carriers (except the local molecule)
c
        do ii = 1, nlplocnl
          ilp = lpglobnl(ii)
          i = lpatom(ilp) 
          call potmole1a(i,molcule(i),0,Vtot1(ilp),dVtot1(:,:,ilp),
     $     dVtot1t(:,:,ilp),0,lppotlst(1,ii),nlppotlst(ii))
        end do
c
c     send and receive potentials and derivatives to/from the neighboring processes
c
        call commpotlp(Vtot1)
        call commpotlp1(dVtot1,dVtot1t)
      end if
c
c     compute the charge transfer interaction energy and derivatives
c
c     Loop over the lone pairs
c
      do ii = 1, nlplocnl
         ilp = lpglobnl(ii)
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
c
         V = Vtot1(ilp)
c
c     Precompute product of the local lp coordinates with derivatives of the
c      rotation matrixes
c
         call rotlp1(ilp,di,dix,diz)
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
c   for cos1 :
c   correspondance between indexes and number (1,2,3,4)
c   1->k1, 2->i, 3->ixlp(ii), 4-> izlp(ii)
c
c   for cos2 :
c   correspondance between indexes and number (1,2,3,4)
c   1->k2, 2->i, 3->ixlp(ii), 4-> izlp(ii)
c
c
              if (use_ctpot) then
                 call aclear(3*nbloc,dVbis1)
                 call aclear(3*nbloc,dVbis1t)
                 call aclear(3*nbloc,dVbis2)
                 call aclear(3*nbloc,dVbis2t)
              end if
              call aclear(3*4,dcos1)
              call aclear(3*4,dcos2)
              call aclear(3*nbloc,dnum)
c
c     Get geometrical data about the atoms forming the X-H electron accepting bond
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
              cos1=(xralp*xrak1 + yralp*yrak1 + zralp*zrak1)/(ralp*rak1)
c
              scalpk1 = xralp*xrak1 + yralp*yrak1 + zralp*zrak1
c
              dcos1(1,1) = dcos1(1,1) + xralp/(rak1*ralp) - 
     $          xrak1*scalpk1/(ralp*rak1**3)
              dcos1(2,1) = dcos1(2,1) + yralp/(rak1*ralp) - 
     $          yrak1*scalpk1/(ralp*rak1**3)
              dcos1(3,1) = dcos1(3,1) + zralp/(rak1*ralp) - 
     $          zrak1*scalpk1/(ralp*rak1**3)
c
              dcos1(1,2) = dcos1(1,2) - (xrak1+xralp)/(rak1*ralp) + 
     $          xralp*scalpk1/(rak1*ralp**3) + 
     $          xrak1*scalpk1/(ralp*rak1**3)
              dcos1(2,2) = dcos1(2,2) - (yrak1+yralp)/(rak1*ralp) + 
     $          yralp*scalpk1/(rak1*ralp**3) + 
     $          yrak1*scalpk1/(ralp*rak1**3)
              dcos1(3,2) = dcos1(3,2) - (zrak1+zralp)/(rak1*ralp) + 
     $          zralp*scalpk1/(rak1*ralp**3) + 
     $          zrak1*scalpk1/(ralp*rak1**3)
c
              dcoslp1 = xrak1/(rak1*ralp) -xralp*scalpk1/(rak1*ralp**3) 
              dcoslp2 = yrak1/(rak1*ralp) -yralp*scalpk1/(rak1*ralp**3) 
              dcoslp3 = zrak1/(rak1*ralp) -zralp*scalpk1/(rak1*ralp**3) 
c
c       contribution of the rotation matrixes
c
              dcos1(1,3) = dcos1(1,3) + dix(1,1)*dcoslp1 +
     $  dix(1,2)*dcoslp2 + dix(1,3)*dcoslp3
              dcos1(2,3) = dcos1(2,3) + dix(2,1)*dcoslp1 +
     $  dix(2,2)*dcoslp2 + dix(2,3)*dcoslp3
              dcos1(3,3) = dcos1(3,3) + dix(3,1)*dcoslp1 +
     $  dix(3,2)*dcoslp2 + dix(3,3)*dcoslp3
c     
              dcos1(1,4) = dcos1(1,4) + diz(1,1)*dcoslp1 +
     $  diz(1,2)*dcoslp2 + diz(1,3)*dcoslp3
              dcos1(2,4) = dcos1(2,4) + diz(2,1)*dcoslp1 +
     $  diz(2,2)*dcoslp2 + diz(2,3)*dcoslp3
              dcos1(3,4) = dcos1(3,4) + diz(3,1)*dcoslp1 +
     $  diz(3,2)*dcoslp2 + diz(3,3)*dcoslp3
c     
              dcos1(1,2) = dcos1(1,2) + di(1,1)*dcoslp1 + 
     $  di(1,2)*dcoslp2 + di(1,3)*dcoslp3
              dcos1(2,2) = dcos1(2,2) + di(2,1)*dcoslp1 +
     $  di(2,2)*dcoslp2 + di(2,3)*dcoslp3
              dcos1(3,2) = dcos1(3,2) + di(3,1)*dcoslp1 +
     $  di(3,2)*dcoslp2 + di(3,3)*dcoslp3
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
                if (use_ctpot) call potmole1aloc(k2,molcule(i),Vtemp,
     $            dVbis1,dVbis1t,dVtot2(:,:,kaccept,2),
     $           dVtot2t(:,:,kaccept,2))
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
                if (use_ctpot) call potmole1aloc(k1,molcule(i),Vtemp,
     $             dVbis2,dVbis2t,
     $             dVtot2(:,:,kaccept,1),dVtot2t(:,:,kaccept,1))
                Vbis2 = Vbis2 - Vtemp
                Vk1k2 = Vk1k2 + Vbis2
              else
c
c     empiric values for cations, Zn
c                 
                Vk1k2 = 1497.102d0
                if (use_ctpot) call potmole1aloc(k1,molcule(i),Vtemp,
     $            dVbis2,dVbis2t,
     $            dVtot2(:,:,kaccept,1),dVtot2t(:,:,kaccept,1))
                Vbis2 = Vbis2 - Vtemp
              end if
              rvdwk1 = rvdwct2(k1)
              rhoak1 = rak1/(4*sqrt(rvdwi*rvdwk1))
c
              task1 = tas(atomic(k1),atomic(i))
              tapk1 = tap(atomic(k1),atomic(i))
              mak1 = ma(atomic(k1),atomic(i))
c
              if (acceptor(2,kk).ne.0) then
c
              cos2=(xralp*xrak2 + yralp*yrak2+zralp*zrak2)/(ralp*rak2)
c
              scalpk2 = xralp*xrak2 + yralp*yrak2 + zralp*zrak2
c
              dcos2(1,1) = dcos2(1,1) + xralp/(rak2*ralp) - 
     $          xrak2*scalpk2/(ralp*rak2**3)
              dcos2(2,1) = dcos2(2,1) + yralp/(rak2*ralp) - 
     $          yrak2*scalpk2/(ralp*rak2**3)
              dcos2(3,1) = dcos2(3,1) + zralp/(rak2*ralp) - 
     $          zrak2*scalpk2/(ralp*rak2**3)
c
              dcos2(1,2) = dcos2(1,2) - (xrak2+xralp)/(rak2*ralp) + 
     $          xralp*scalpk2/(rak2*ralp**3) + 
     $          xrak2*scalpk2/(ralp*rak2**3)
              dcos2(2,2) = dcos2(2,2) - (yrak2+yralp)/(rak2*ralp) + 
     $          yralp*scalpk2/(rak2*ralp**3) + 
     $          yrak2*scalpk2/(ralp*rak2**3)
              dcos2(3,2) = dcos2(3,2) - (zrak2+zralp)/(rak2*ralp) + 
     $          zralp*scalpk2/(rak2*ralp**3) + 
     $          zrak2*scalpk2/(ralp*rak2**3)
c     
              dcoslp1 = xrak2/(rak2*ralp) - xralp*scalpk2/(rak2*ralp**3)
              dcoslp2 = yrak2/(rak2*ralp) - yralp*scalpk2/(rak2*ralp**3)
              dcoslp3 = zrak2/(rak2*ralp) - zralp*scalpk2/(rak2*ralp**3)
c
c       contribution of the rotation matrixes
c
              dcos2(1,3) = dcos2(1,3) + dix(1,1)*dcoslp1 +
     $  dix(1,2)*dcoslp2 + dix(1,3)*dcoslp3
              dcos2(2,3) = dcos2(2,3) + dix(2,1)*dcoslp1 +
     $  dix(2,2)*dcoslp2 + dix(2,3)*dcoslp3
              dcos2(3,3) = dcos2(3,3) + dix(3,1)*dcoslp1 +
     $  dix(3,2)*dcoslp2 + dix(3,3)*dcoslp3
c     
              dcos2(1,4) = dcos2(1,4) + diz(1,1)*dcoslp1 +
     $  diz(1,2)*dcoslp2 + diz(1,3)*dcoslp3
              dcos2(2,4) = dcos2(2,4) + diz(2,1)*dcoslp1 +
     $  diz(2,2)*dcoslp2 + diz(2,3)*dcoslp3
              dcos2(3,4) = dcos2(3,4) + diz(3,1)*dcoslp1 +
     $  diz(3,2)*dcoslp2 + diz(3,3)*dcoslp3
c     
              dcos2(1,2) = dcos2(1,2) + di(1,1)*dcoslp1 + 
     $  di(1,2)*dcoslp2 + di(1,3)*dcoslp3
              dcos2(2,2) = dcos2(2,2) + di(2,1)*dcoslp1 +
     $  di(2,2)*dcoslp2 + di(2,3)*dcoslp3
              dcos2(3,2) = dcos2(3,2) + di(3,1)*dcoslp1 +
     $  di(3,2)*dcoslp2 + di(3,3)*dcoslp3
c
                task2 = tas(atomic(k2),atomic(i))
                tapk2 = tap(atomic(k2),atomic(i))
                mak2 = ma(atomic(k2),atomic(i))
              end if
c
c     Final assembling of all the terms of the numerator
c
              Cs = hybrid_lp(1,ilp)
              Cp = hybrid_lp(2,ilp)
c
              e1 = 0.0d0
              e2 = exp(-etaect*rhoak2)*(D2*Vk2k1-V)*(Cs*task2+
     $             Cp*mak2*tapk2*cos2)
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
c            Derivative of the e2 term
c
              dnum(1,loc(i)) =  
     $         etaect*xrak2/(4*rak2*sqrt(rvdwi*rvdwk2))*
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*(Cs*task2 + 
     $         Cp*mak2*tapk2*cos2) + exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(1,2)
              dnum(2,loc(i)) = 
     $         etaect*yrak2/(4*rak2*sqrt(rvdwi*rvdwk2))*
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*(Cs*task2 +
     $         Cp*mak2*tapk2*cos2) + exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(2,2)
              dnum(3,loc(i)) = dnum(3,loc(i)) +
     $         etaect*zrak2/(4*rak2*sqrt(rvdwi*rvdwk2))*
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*(Cs*task2 +
     $         Cp*mak2*tapk2*cos2) + exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(3,2)
c
              dnum(1,loc(k2)) = dnum(1,loc(k2)) - 
     $         etaect*xrak2/(4*rak2*sqrt(rvdwi*rvdwk2))*
     $         exp(-etaect*rhoak2)*(D1*Vk2k1-V)*(Cs*task2 + 
     $         Cp*mak2*tapk2*cos2) + exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(1,1)
              dnum(2,loc(k2)) = dnum(2,loc(k2)) -
     $         etaect*yrak2/(4*rak2*sqrt(rvdwi*rvdwk2))*
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*(Cs*task2 + 
     $         Cp*mak2*tapk2*cos2) + exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(2,1) 
              dnum(3,loc(k2)) = dnum(3,loc(k2)) -
     $         etaect*zrak2/(4*rak2*sqrt(rvdwi*rvdwk2))*
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*(Cs*task2 + 
     $         Cp*mak2*tapk2*cos2) + exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(3,1)
c
              dnum(1,loc(ixlp(ilp))) = dnum(1,loc(ixlp(ilp))) + 
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(1,3)
              dnum(2,loc(ixlp(ilp))) = dnum(2,loc(ixlp(ilp))) + 
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(2,3)
              dnum(3,loc(ixlp(ilp))) = dnum(3,loc(ixlp(ilp))) + 
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(3,3)
c
              dnum(1,loc(izlp(ilp))) = dnum(1,loc(izlp(ilp))) + 
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(1,4)
              dnum(2,loc(izlp(ilp))) = dnum(2,loc(izlp(ilp))) + 
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(2,4)
              dnum(3,loc(izlp(ilp))) = dnum(3,loc(izlp(ilp))) + 
     $         exp(-etaect*rhoak2)*(D2*Vk2k1-V)*
     $         Cp*mak2*tapk2*dcos2(3,4)
c
              do  j = 1, nbloc
                  dnum(1,j) = dnum(1,j) + (D2*(dVbis1(1,j)+
     $              dVbis1t(1,j))-dVtot1(1,j,ilp)-dVtot1t(1,j,ilp))
     $               *exp(-etaect*rhoak2)*
     $               (Cs*task2 + Cp*mak2*tapk2*cos2)
                  dnum(2,j) = dnum(2,j) + (D2*(dVbis1(2,j)+
     $              dVbis1t(2,j))-dVtot1(2,j,ilp)-dVtot1t(2,j,ilp))
     $               *exp(-etaect*rhoak2)*
     $               (Cs*task2 + Cp*mak2*tapk2*cos2)
                  dnum(3,j) = dnum(3,j) + (D2*(dVbis1(3,j)+
     $              dVbis1t(3,j))-dVtot1(3,j,ilp)-dVtot1t(3,j,ilp))
     $               *exp(-etaect*rhoak2)*
     $               (Cs*task2 + Cp*mak2*tapk2*cos2)
              end do
c
c            Derivative of the e1 term
c
              dnum(1,loc(i)) = dnum(1,loc(i)) - 
     $         etaect*xrak1/(4*rak1*sqrt(rvdwi*rvdwk1))*
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*(Cs*task1 + 
     $         Cp*mak1*tapk1*cos1) - exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(1,2)
              dnum(2,loc(i)) = dnum(2,loc(i)) - 
     $         etaect*yrak1/(4*rak1*sqrt(rvdwi*rvdwk1))*
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*(Cs*task1 + 
     $         Cp*mak1*tapk1*cos1) - exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(2,2)
              dnum(3,loc(i)) = dnum(3,loc(i)) - 
     $         etaect*zrak1/(4*rak1*sqrt(rvdwi*rvdwk1))*
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*(Cs*task1 + 
     $         Cp*mak1*tapk1*cos1) - exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(3,2)
c
              dnum(1,loc(k1)) = dnum(1,loc(k1)) + 
     $         etaect*xrak1/(4*rak1*sqrt(rvdwi*rvdwk1))*
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*(Cs*task1 + 
     $         Cp*mak1*tapk1*cos1) - exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(1,1)
              dnum(2,loc(k1)) = dnum(2,loc(k1)) + 
     $         etaect*yrak1/(4*rak1*sqrt(rvdwi*rvdwk1))*
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*(Cs*task1 + 
     $         Cp*mak1*tapk1*cos1) - exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(2,1)
              dnum(3,loc(k1)) = dnum(3,loc(k1)) + 
     $         etaect*zrak1/(4*rak1*sqrt(rvdwi*rvdwk1))*
     $         exp(-etaect*rhoak1)*(D2*Vk1k2-V)*(Cs*task1 + 
     $         Cp*mak1*tapk1*cos1) - exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(3,1)
c
              dnum(1,loc(ixlp(ilp))) = dnum(1,loc(ixlp(ilp))) - 
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(1,3)
              dnum(2,loc(ixlp(ilp))) = dnum(2,loc(ixlp(ilp))) - 
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(2,3)
              dnum(3,loc(ixlp(ilp))) = dnum(3,loc(ixlp(ilp))) - 
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(3,3)
c
              dnum(1,loc(izlp(ilp))) = dnum(1,loc(izlp(ilp))) - 
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(1,4)
              dnum(2,loc(izlp(ilp))) = dnum(2,loc(izlp(ilp))) - 
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(2,4)
              dnum(3,loc(izlp(ilp))) = dnum(3,loc(izlp(ilp))) - 
     $         exp(-etaect*rhoak1)*(D1*Vk1k2-V)*
     $         Cp*mak1*tapk1*dcos1(3,4)
c
              do  j = 1, nbloc
                  dnum(1,j) = dnum(1,j) - (D1*(dVbis2(1,j)+
     $              dVbis2t(1,j))-dVtot1(1,j,ilp)-dVtot1t(1,j,ilp))
     $                *exp(-etaect*rhoak1)*
     $               (Cs*task1 + Cp*mak1*tapk1*cos1)
                  dnum(2,j) = dnum(2,j) - (D1*(dVbis2(2,j)+
     $              dVbis2t(2,j))-dVtot1(2,j,ilp)-dVtot1t(2,j,ilp))
     $                *exp(-etaect*rhoak1)*
     $               (Cs*task1 + Cp*mak1*tapk1*cos1)
                  dnum(3,j) = dnum(3,j) - (D1*(dVbis2(3,j)+
     $              dVbis2t(3,j))-dVtot1(3,j,ilp)-dVtot1t(3,j,ilp))
     $                *exp(-etaect*rhoak1)*
     $               (Cs*task1 + Cp*mak1*tapk1*cos1)
              end do
c
              do j = 1, nbloc
                dnum(1,j) = -2*clp*dnum(1,j)*(e2-e1) 
                dnum(2,j) = -2*clp*dnum(2,j)*(e2-e1) 
                dnum(3,j) = -2*clp*dnum(3,j)*(e2-e1) 
              end do
c
              Vbis1 = Vbis1/23.06
              Vbis2 = Vbis2/23.06
              Vconv = V/23.06
c
              if (acceptor(2,kk).ne.0) then
                deltae = ahct(i) + Vconv  - (aelct(k2) + Vbis2) 
              else
c
c  for the zinc
c
                deltae = ahct(k1) + Vconv - (17.99d0 + Vbis2) 
              end if
c
              e = - clp*(e1 - e2)**2/deltae
c
c     Facteur de Conversion
c
              conv = 27.2/627.5442
              if (acceptor(2,kk).ne.0) then
                prop = 0.6364d0 
              else
                prop = 0.61d0
              end if
              e = e*conv
              e = prop*e
              ect = ect + e
c
              do j = 1, nbloc
                dectt(1,j) = dectt(1,j) + (prop*conv)*(dnum(1,j)*deltae+
     $            clp*(e1-e2)**2*((dVtot1(1,j,ilp)+dVtot1t(1,j,ilp)
     $            -dVtot2(1,j,kaccept,2)-dVtot2t(1,j,kaccept,2))
     $            /23.06))/deltae**2
                dectt(2,j) = dectt(2,j) + (prop*conv)*(dnum(2,j)*deltae+
     $            clp*(e1-e2)**2*((dVtot1(2,j,ilp)+dVtot1t(2,j,ilp)
     $            -dVtot2(2,j,kaccept,2)-dVtot2t(2,j,kaccept,2))
     $            /23.06))/deltae**2
                dectt(3,j) = dectt(3,j) + (prop*conv)*(dnum(3,j)*deltae+
     $            clp*(e1-e2)**2*((dVtot1(3,j,ilp)+dVtot1t(3,j,ilp)
     $            -dVtot2(3,j,kaccept,2)-dVtot2t(3,j,kaccept,2))
     $            /23.06))/deltae**2
              end do
c            
c     increment the overall charge transfer energy component
c
              nect = nect + 1
c
c     increment the total intermolecular energy
c
              if (molcule(i) .ne. molcule(k1)) then
                einter = einter + e
              end if
            end if     
         end do
       end do
c
c      move forces to global arrays
c
       do i = 1, nbloc
         iglob = glob(i)
         dect(1,i) = dectt(1,i)
         dect(2,i) = dectt(2,i)
         dect(3,i) = dectt(3,i)
       end do
c
       deallocate (dectt)
       deallocate (dVtot1,dVtot1t)
       deallocate (dVtot2,dVtot2t)
       deallocate (Vtot1)
       deallocate (Vtot2)
       deallocate (dnum)
       deallocate (dV)
       deallocate (dVt)
       deallocate (dcos1)
       deallocate (dVbis1)
       deallocate (dVbis1t)
       deallocate (dVbis2)
       deallocate (dVbis2t)
       deallocate (dcos2)
       return
       end
c
c
c     "potmole1a" calculates the electrostatic potential created by the whole
c     system except some molecule on a site and its derivative
c
c     - rule = 0 : compute the whole potential on ia except the one created by imol1
c     - rule = 1 : compute the whole potential on ia except the one created by imol1 and imol2 
c
c
      subroutine potmole1a(ia,imol1,imol2,V,dV,dVt,rule,list,nlist)
      use action
      use analyz
      use atmtyp
      use atoms
      use atmlst
      use bound
      use charge
      use chargetransfer
      use chgpot
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use inter
      use iounit
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use usage
      implicit none
      integer rule,imol1,imol2
      integer i,j,k,l,k1,m,p,n1,n2,q1
      integer ii,kk,ia,im,ilp
      integer ialoc,iglob,iipole,iloc
      integer kkpole,kglob,inl,kkk,kbis
      integer in,ic,kn,kc
      integer natom
      integer list(maxelst),nlist
      real*8 e,ei,fgrp
      real*8 f,fm,fp,fi,fik
      real*8 xi,yi,zi
      real*8 rc4,rc5,rc6,rc7
      real*8 r,rb2,r2,xr,yr,zr
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz,qiyx,qiyy,qiyz,qizx,qizy,qizz
      real*8 ck
      real*8 qix,qiy,qiz
      real*8 scale3,scale5
      real*8 scale7
      real*8 d(3),qr(3,3),q2(3,3)
      real*8 scd1,scd2,scd3,sq1x,sq1y,sq1z,scq1,sq2x,sq2y,sq2z,scq2
      real*8 sq3x,sq3y,sq3z,scq3
      real*8 rin,rout,s,ds
      real*8 ftm2(3)
      real*8 fd1,fd2,scd,fq1,fq2,scq,fbis,fc
      real*8 dV(3,*),dVt(3,*),V
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
      ialoc = loc(ia)
      iipole = pollist(ia)
c
      do kkk = 1, nlist
         kkpole = list(kkk)
         kglob = ipole(kkpole)
         kbis = poleloc(kglob)
         k = loc(kglob)
         if (kbis.eq.0) then
           write(iout,1000)
           goto 100
         end if
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
          qiyx = rpole(8,kkpole)
          qiyy = rpole(9,kkpole)
          qiyz = rpole(10,kkpole)
          qizx = rpole(11,kkpole)
          qizy = rpole(12,kkpole)
          qizz = rpole(13,kkpole)
          xr = x(ia) - x(kglob)
          yr = y(ia) - y(kglob)
          zr = z(ia) - z(kglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr*yr + zr*zr
          if (r2.gt.(mpolectcut*mpolectcut)) cycle
          r = sqrt(r2)
c
c    get switchung function terms
c
          call switch_emtp(.true.,r,mpolectcut-rin,mpolectcut,s,ds)
c
          qix = qixx*xr + qixy*yr + qixz*zr
          qiy = qixy*xr + qiyy*yr + qiyz*zr
          qiz = qixz*xr + qiyz*yr + qizz*zr
          scd = dix*xr + diy*yr + diz*zr
          scq = qix*xr + qiy*yr + qiz*zr
c                
          rr1 = 1.0d0 / r
          rr3 = rr1 / r2
          rr5 = 3.0d0 * rr3 / r2
          rr7 = 5.0d0 * rr5 / r2
          fc  = - ci*rr3
          fd1 = rr3
          fd2 = -scd*rr5
          fq1 = rr5/2.0d0
          fq2 = -scq*rr7/2.0d0
          e  = ci*rr1 + scd*rr3 + scq*rr5/2.0d0
          ftm2(1) = fbis*(fd1*dix + fq1*2*qix + (fc + fd2 + fq2)*xr)
          ftm2(2) = fbis*(fd1*diy + fq1*2*qiy + (fc + fd2 + fq2)*yr)
          ftm2(3) = fbis*(fd1*diz + fq1*2*qiz + (fc + fd2 + fq2)*zr)
c
          e = fbis * e
c
          V = V + s*e
c
          dV(1,k) = dV(1,k) - s*ftm2(1) - ds*e*xr*rr1
          dV(2,k) = dV(2,k) - s*ftm2(2) - ds*e*yr*rr1
          dV(3,k) = dV(3,k) - s*ftm2(3) - ds*e*zr*rr1
          dV(1,ialoc) = dV(1,ialoc) + s*ftm2(1) + ds*e*xr*rr1
          dV(2,ialoc) = dV(2,ialoc) + s*ftm2(2) + ds*e*yr*rr1
          dV(3,ialoc) = dV(3,ialoc) + s*ftm2(3) + ds*e*zr*rr1
        end if
 100    continue
      end do
c
c     get the contribution due to the derivatives of the rotation matrixes
c
      call derrotpotmol(ia,imol1,imol2,dVt,rule,list,nlist)
c
      return
      end
c
c     subroutine derrotpotmol : computes the contribution to the forces (electrostatic part) due to the derivatives of the rotation matrixes 
c
      subroutine derrotpotmol(ia,imol1,imol2,dVt,rule,list,nlist)
      use atoms
      use atmlst
      use bound
      use chargetransfer
      use chgpot
      use cutoff
      use deriv
      use domdec
      use iounit
      use molcul
      use mpole
      use neigh
      use sizes
      use shunt
      implicit none
      integer rule,imol1,imol2
      integer i,j,ii,kz,kx,ky,ia,inl,iipole
      integer kxloc,kzloc,kyloc,kglob,kkpole,kbis,kkk
      integer ilp,im,natom
      integer l,m,n1,n2,p,k
      integer list(maxelst),nlist
      real*8 ck,xr,yr,zr,r2,rr3,rr5
      real*8 emol,fm,fbis,r, dx(3)
      real*8 d(3),di(3,3),q2(3,3),qr(3,3),qi(3,3,3)
      real*8 rmat(3,3),dri(3,3,3),driz(3,3,3),drix(3,3,3),driy(3,3,3)
      real*8 scqx, scqy, scqz
      real*8 rin,rout,s,ds
      real*8 dVt(3,*)
      logical doi,doix,doiz,doiy
      character*6 mode
c
 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')
c
      fbis = 331.93359 
c
      mode = 'EWALD'
      call switch (mode)
c
      rin = 1.0d0
      rout = 2.0d0
c
      iipole = pollist(ia)
c      
      do kkk = 1, nlist
         kkpole = list(kkk)
         kglob = ipole(kkpole)
         kbis = poleloc(kglob)
         k = loc(kglob)
         if (kbis.eq.0) then
           write(iout,1000)
           goto 100
         end if
        if (((rule.eq.0).and.(molcule(kglob).ne.imol1)).or.
     $    ((rule.eq.1).and.(molcule(kglob).ne.imol1).and.
     $    (molcule(kglob).ne.imol2))) then
        kz = zaxis(kkpole)
        if (kz.ne.0) kzloc = loc(kz)
        kx = xaxis(kkpole)
        if (kx.ne.0) kxloc = loc(kx)
        ky = yaxis(kkpole)
        if (ky.ne.0) kyloc = loc(ky)
        xr = x(ia) - x(kglob)
        yr = y(ia) - y(kglob)
        zr = z(ia) - z(kglob)
        if (use_bounds)  call image (xr,yr,zr)
        dx(1) = xr
        dx(2) = yr
        dx(3) = zr
        r2 = xr*xr + yr*yr + zr*zr
        if (r2.gt.(mpolectcut+rout)*(mpolectcut+rout)) cycle
c
        r = sqrt(r2)
c
c    get switchung function terms
c
        call switch_emtp(.true.,r,mpolectcut-rin,mpolectcut+rout,s,
     $   ds)
        rr3 = 1.0d0 / (r2*r)
        rr5 = 3.0d0 * rr3 / r2
c        
        d(1)    = pole(2,kkpole)
        d(2)    = pole(3,kkpole)
        d(3)    = pole(4,kkpole)
        q2(1,1)  = pole(5,kkpole)
        q2(1,2)  = pole(6,kkpole)
        q2(1,3)  = pole(7,kkpole)
        q2(2,1)  = pole(8,kkpole) 
        q2(2,2)  = pole(9,kkpole) 
        q2(2,3)  = pole(10,kkpole) 
        q2(3,1)  = pole(11,kkpole) 
        q2(3,2)  = pole(12,kkpole) 
        q2(3,3)  = pole(13,kkpole) 
        qr(1,1) = rpole(5,kkpole)
        qr(1,2) = rpole(6,kkpole)
        qr(1,3) = rpole(7,kkpole)
        qr(2,1) = rpole(8,kkpole) 
        qr(2,2) = rpole(9,kkpole) 
        qr(2,3) = rpole(10,kkpole) 
        qr(3,1) = rpole(11,kkpole) 
        qr(3,2) = rpole(12,kkpole) 
        qr(3,3) = rpole(13,kkpole)
c
c     Initialisation
c
        do l = 1, 3
          do m = 1, 3
              di(l,m) = 0.0d0
              do n1 = 1, 3
               qi(l,m,n1) = 0.0d0
              end do
          end do
        end do
c       
        call derrot(kglob,.true.,kkpole,kz,kx,ky,rmat,dri,driz,drix,
     $     driy)
c
        doi = .true.
        doix = .false.
        doiy = .false.
        doiz = .false.
        if (kz.ne.0) then
          doiz = .true.
        end if
        if (kx.ne.0) then
          doix = .true.
        end if
        if (ky.ne.0) then
          doiy = .true.
        end if
c
        fm = fbis
        if (doi) then
          call aclear(9,di)
          call aclear(27,qi)
          do l = 1, 3
            do m = 1, 3
              di(1,l) = di(1,l) + dri(1,l,m)*d(m)
              di(2,l) = di(2,l) + dri(2,l,m)*d(m)
              di(3,l) = di(3,l) + dri(3,l,m)*d(m)
              do n1 = 1, 3
                do n2 = 1, 3
                  qi(1,l,m) = qi(1,l,m) + 
     $                        q2(n1,n2)*(dri(1,l,n1)*rmat(m,n2) +
     $                        dri(1,m,n2)*rmat(l,n1))
                  qi(2,l,m) = qi(2,l,m) + 
     $                        q2(n1,n2)*(dri(2,l,n1)*rmat(m,n2) +
     $                        dri(2,m,n2)*rmat(l,n1))
                  qi(3,l,m) = qi(3,l,m) + 
     $                        q2(n1,n2)*(dri(3,l,n1)*rmat(m,n2) +
     $                        dri(3,m,n2)*rmat(l,n1))
                end do
              end do
            end do
          end do
c          
          scqx = 0.0d0
          scqy = 0.0d0
          scqz = 0.0d0
          do l = 1, 3
            do m = 1, 3
              scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
              scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
              scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
            end do
          end do
          dVt(1,k) = dVt(1,k) + fm*s*rr3*(di(1,1)*xr +
     $     di(1,2)*yr + di(1,3)*zr) + scqx
          dVt(2,k) = dVt(2,k) + fm*s*rr3*(di(2,1)*xr +
     $     di(2,2)*yr + di(2,3)*zr) + scqy     
          dVt(3,k) = dVt(3,k) + fm*s*rr3*(di(3,1)*xr +
     $     di(3,2)*yr + di(3,3)*zr) + scqz      
        end if
c
c       derivative with respect to the position of the iz defining atom
c
        if (doiz) then
          call aclear(9,di)
          call aclear(27,qi)
          do l = 1, 3
            do m = 1, 3
              di(1,l) = di(1,l) + driz(1,l,m)*d(m)
              di(2,l) = di(2,l) + driz(2,l,m)*d(m)
              di(3,l) = di(3,l) + driz(3,l,m)*d(m)
              do n1 = 1, 3
                do n2 = 1, 3
                  qi(1,l,m) = qi(1,l,m) + 
     $                        q2(n1,n2)*(driz(1,l,n1)*rmat(m,n2) +
     $                        driz(1,m,n2)*rmat(l,n1))
                  qi(2,l,m) = qi(2,l,m) + 
     $                        q2(n1,n2)*(driz(2,l,n1)*rmat(m,n2) +
     $                        driz(2,m,n2)*rmat(l,n1))
                  qi(3,l,m) = qi(3,l,m) + 
     $                        q2(n1,n2)*(driz(3,l,n1)*rmat(m,n2) +
     $                        driz(3,m,n2)*rmat(l,n1))
                end do
              end do
            end do
          end do
c          
          scqx = 0.0d0
          scqy = 0.0d0
          scqz = 0.0d0
          do l = 1, 3
            do m = 1, 3
              scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
              scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
              scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
            end do
          end do
c          
          dVt(1,kzloc) = dVt(1,kzloc) + fm*s*rr3*(di(1,1)*xr +
     $     di(1,2)*yr + di(1,3)*zr) + scqx    
          dVt(2,kzloc) = dVt(2,kzloc) + fm*s*rr3*(di(2,1)*xr +
     $     di(2,2)*yr + di(2,3)*zr) + scqy      
          dVt(3,kzloc) = dVt(3,kzloc) + fm*s*rr3*(di(3,1)*xr +
     $     di(3,2)*yr + di(3,3)*zr) + scqz      
        end if
c
c        derivative with respect to the position of the ix defining atom
c
        if (doix) then
          call aclear(9,di)
          call aclear(27,qi)
          do l = 1, 3
            do m = 1, 3
              di(1,l) = di(1,l) + drix(1,l,m)*d(m)
              di(2,l) = di(2,l) + drix(2,l,m)*d(m)
              di(3,l) = di(3,l) + drix(3,l,m)*d(m)
              do n1 = 1, 3
                do n2 = 1, 3
                  qi(1,l,m) = qi(1,l,m) + 
     $                        q2(n1,n2)*(drix(1,l,n1)*rmat(m,n2) +
     $                        drix(1,m,n2)*rmat(l,n1))
                  qi(2,l,m) = qi(2,l,m) + 
     $                        q2(n1,n2)*(drix(2,l,n1)*rmat(m,n2) +
     $                        drix(2,m,n2)*rmat(l,n1))
                  qi(3,l,m) = qi(3,l,m) + 
     $                        q2(n1,n2)*(drix(3,l,n1)*rmat(m,n2) +
     $                        drix(3,m,n2)*rmat(l,n1))
                end do
              end do
            end do
          end do
c          
          scqx = 0.0d0
          scqy = 0.0d0
          scqz = 0.0d0
          do l = 1, 3
            do m = 1, 3
              scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
              scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
              scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
            end do
          end do
c          
          dVt(1,kxloc) = dVt(1,kxloc) + fm*s*rr3*(di(1,1)*xr +
     $     di(1,2)*yr + di(1,3)*zr) + scqx      
          dVt(2,kxloc) = dVt(2,kxloc) + fm*s*rr3*(di(2,1)*xr +
     $     di(2,2)*yr + di(2,3)*zr) + scqy      
          dVt(3,kxloc) = dVt(3,kxloc) + fm*s*rr3*(di(3,1)*xr +
     $     di(3,2)*yr + di(3,3)*zr) + scqz      
        end if
c
c
c        derivative with respect to the position of the iy defining atom
c
         if (doiy) then
           call aclear(9,di)
           call aclear(27,qi)
           do l = 1, 3
             do m = 1, 3
               di(1,l) = di(1,l) + driy(1,l,m)*d(m)
               di(2,l) = di(2,l) + driy(2,l,m)*d(m)
               di(3,l) = di(3,l) + driy(3,l,m)*d(m)
               do n1 = 1, 3
                 do n2 = 1, 3
                   qi(1,l,m) = qi(1,l,m) + 
     $                         q2(n1,n2)*(driy(1,l,n1)*rmat(m,n2) +
     $                         driy(1,m,n2)*rmat(l,n1))
                   qi(2,l,m) = qi(2,l,m) + 
     $                         q2(n1,n2)*(driy(2,l,n1)*rmat(m,n2) +
     $                         driy(2,m,n2)*rmat(l,n1))
                   qi(3,l,m) = qi(3,l,m) + 
     $                         q2(n1,n2)*(driy(3,l,n1)*rmat(m,n2) +
     $                         driy(3,m,n2)*rmat(l,n1))
                 end do
               end do
             end do
           end do
c           
           scqx = 0.0d0
           scqy = 0.0d0
           scqz = 0.0d0
           do l = 1, 3
             do m = 1, 3
               scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
               scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
               scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
             end do
           end do
c           
           dVt(1,kyloc) = dVt(1,kyloc) + fm*s*rr3*(di(1,1)*xr +
     $      di(1,2)*yr + di(1,3)*zr) + scqx      
           dVt(2,kyloc) = dVt(2,kyloc) + fm*s*rr3*(di(2,1)*xr +
     $      di(2,2)*yr + di(2,3)*zr) + scqy      
           dVt(3,kyloc) = dVt(3,kyloc) + fm*s*rr3*(di(3,1)*xr +
     $      di(3,2)*yr + di(3,3)*zr) + scqy      
         end if
c
       end if
  100   continue
      end do
      return 
      end 
c
c     subroutine potmole1aloc : computes the electrostatic potential created by molecule 
c     imol2 on ia and derivatives
c
      subroutine potmole1aloc(ia,imol2,V,dV,dVt,dVtot,dVtott)
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
      use deriv
      use domdec
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use usage
      implicit none
      integer imol2
      integer i,j,k,l,k1,m,p,n1,n2,q1
      integer ii,kk,ia,im,ilp
      integer ialoc,iglob,iipole,iloc
      integer ix,iy,iz
      integer in,ic,kn,kc
      integer kx,ky,kz,natom
      real*8 e,ei,fgrp,V
      real*8 f,fm,fp,fi,fik
      real*8 xi,yi,zi
      real*8 rc4,rc5,rc6,rc7
      real*8 r,rb2,r2,xr,yr,zr
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz,qiyx,qiyy,qiyz,qizx,qizy,qizz
      real*8 ck
      real*8 qix,qiy,qiz
      real*8 scale3,scale5
      real*8 scale7
      real*8 d(3),qr(3,3),q2(3,3)
      real*8 scd1,scd2,scd3,sq1x,sq1y,sq1z,scq1,sq2x,sq2y,sq2z,scq2
      real*8 sq3x,sq3y,sq3z,scq3
      real*8 rin,rout,s,ds
      real*8 ftm2(3)
      real*8 fd1,fd2,scd,fq1,fq2,scq,fbis,fc
      real*8 dV(3,*),dVt(3,*),dVtot(3,*),dVtott(3,*)
      character*6 mode
c
      ialoc = loc(ia)
c
c     check whether the molecules imol1 and imol2 exists
c
      if (imol2.gt.nmol) then
        write (iout,10)
   10   format (/,' POTMOL1Aloc --  Unable to Find the molecule')
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
      rout = 1.0d0
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
          qiyx = rpole(8,ii)
          qiyy = rpole(9,ii)
          qiyz = rpole(10,ii)
          qizx = rpole(11,ii)
          qizy = rpole(12,ii)
          qizz = rpole(13,ii)
          xr = x(ia) - x(iglob)
          yr = y(ia) - y(iglob)
          zr = z(ia) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr*yr + zr*zr
          if (r2.gt.(mpolectcut+rout)*(mpolectcut+rout)) cycle
          r = sqrt(r2)
c
c    get switchung function terms
c
          call switch_emtp(.true.,r,mpolectcut-rin,mpolectcut+rout,s,
     $     ds)
c
          qix = qixx*xr + qixy*yr + qixz*zr
          qiy = qixy*xr + qiyy*yr + qiyz*zr
          qiz = qixz*xr + qiyz*yr + qizz*zr
          scd = dix*xr + diy*yr + diz*zr
          scq = qix*xr + qiy*yr + qiz*zr
c              
          rr1 = 1.0d0 / r
          rr3 = rr1 / r2
          rr5 = 3.0d0 * rr3 / r2
          rr7 = 5.0d0 * rr5 / r2
          fc  = - ci*rr3
          fd1 = rr3
          fd2 = -scd*rr5
          fq1 = rr5/2.0d0
          fq2 = -scq*rr7/2.0d0
          e  = ci*rr1 + scd*rr3 + scq*rr5/2.0d0
          ftm2(1) = fbis*(fd1*dix + fq1*2*qix + (fc + fd2 + fq2)*xr)
          ftm2(2) = fbis*(fd1*diy + fq1*2*qiy + (fc + fd2 + fq2)*yr)
          ftm2(3) = fbis*(fd1*diz + fq1*2*qiz + (fc + fd2 + fq2)*zr)
c
          e = fbis * e
c
          V = V + s*e
c
          dV(1,iloc) = dV(1,iloc) - s*ftm2(1) -ds*xr*e*rr1
          dV(2,iloc) = dV(2,iloc) - s*ftm2(2) -ds*yr*e*rr1
          dV(3,iloc) = dV(3,iloc) - s*ftm2(3) -ds*zr*e*rr1
          dV(1,ialoc) = dV(1,ialoc) + s*ftm2(1) + ds*xr*e*rr1
          dV(2,ialoc) = dV(2,ialoc) + s*ftm2(2) + ds*yr*e*rr1
          dV(3,ialoc) = dV(3,ialoc) + s*ftm2(3) + ds*zr*e*rr1
      end do
c
c     sum to get the actual force 
c
      do i = 1, nbloc
        dV(1,i) = dVtot(1,i) - dV(1,i)
        dV(2,i) = dVtot(2,i) - dV(2,i)
        dV(3,i) = dVtot(3,i) - dV(3,i)
      end do
c
c     get the contribution due to the derivatives of the rotation matrixes
c
      call derrotpotmolloc(ia,imol2,dVt,dVtott)
c
      return
      end
c
c     
      subroutine derrotpotmolloc(ia,imol2,dVt,dVtott)
      use atoms
      use atmlst
      use bound
      use chargetransfer
      use chgpot
      use cutoff
      use deriv
      use domdec
      use molcul
      use mpole
      use shunt
      implicit none
      integer rule,imol1,imol2
      integer i,j,ii,iz,ix,iy,ia
      integer ixloc,izloc,iyloc,iglob,iipole,iloc
      integer ilp,im,natom
      integer l,m,n1,n2,p,k
      real*8 ck,xr,yr,zr,r2,rr3,rr5
      real*8 emol,fm,fbis,r, dx(3)
      real*8 d(3),di(3,3),q2(3,3),qr(3,3),qi(3,3,3)
      real*8 rmat(3,3),dri(3,3,3),driz(3,3,3),drix(3,3,3),driy(3,3,3)
      real*8 scqx, scqy, scqz
      real*8 rin,rout,s,ds
      real*8 dVt(3,*),dVtott(3,*)
      logical doi,doix,doiz,doiy
      character*6 mode
c
      fbis = 331.93359 
      mode = 'EWALD'
      call switch (mode)
      rin = 1.0d0
      rout = 1.0d0
c
      do i = imol(1,imol2), imol(2,imol2)
        ii = pollist(i)
        iloc = loc(i)
        iglob = i
        iz = zaxis(ii)
        if (iz.ne.0) izloc = loc(iz)
        ix = xaxis(ii)
        if (ix.ne.0) ixloc = loc(ix)
        iy = yaxis(ii)
        if (iy.ne.0) iyloc = loc(iy)
        xr = x(ia) - x(iglob)
        yr = y(ia) - y(iglob)
        zr = z(ia) - z(iglob)
        if (use_bounds)  call image (xr,yr,zr)
        dx(1) = xr
        dx(2) = yr
        dx(3) = zr
        r2 = xr*xr + yr*yr + zr*zr
        if (r2.gt.(mpolectcut*mpolectcut)) cycle
        r = sqrt(r2)
c
c    get switchung function terms
c
        call switch_emtp(.true.,r,mpolectcut-rin,mpolectcut,s,
     $   ds)
        rr3 = 1.0d0 / (r2*r)
        rr5 = 3.0d0 * rr3 / r2
c        
        d(1)    = pole(2,ii)
        d(2)    = pole(3,ii)
        d(3)    = pole(4,ii)
        q2(1,1)  = pole(5,ii)
        q2(1,2)  = pole(6,ii)
        q2(1,3)  = pole(7,ii)
        q2(2,1)  = pole(8,ii) 
        q2(2,2)  = pole(9,ii) 
        q2(2,3)  = pole(10,ii) 
        q2(3,1)  = pole(11,ii) 
        q2(3,2)  = pole(12,ii) 
        q2(3,3)  = pole(13,ii) 
        qr(1,1) = rpole(5,ii)
        qr(1,2) = rpole(6,ii)
        qr(1,3) = rpole(7,ii)
        qr(2,1) = rpole(8,ii) 
        qr(2,2) = rpole(9,ii) 
        qr(2,3) = rpole(10,ii) 
        qr(3,1) = rpole(11,ii) 
        qr(3,2) = rpole(12,ii) 
        qr(3,3) = rpole(13,ii)
c
c     Initialisation
c
        do l = 1, 3
          do m = 1, 3
              di(l,m) = 0.0d0
              do n1 = 1, 3
               qi(l,m,n1) = 0.0d0
              end do
          end do
        end do
c       
        call derrot(iglob,.true.,ii,iz,ix,iy,rmat,dri,driz,drix,driy)
c
        doi = .true.
        doix = .false.
        doiy = .false.
        doiz = .false.
        if (iz.ne.0) then
          doiz = .true.
        end if
        if (ix.ne.0) then
          doix = .true.
        end if
        if (iy.ne.0) then
          doiy = .true.
        end if
c
        fm = fbis
        if (doi) then
          call aclear(9,di)
          call aclear(27,qi)
          do l = 1, 3
            do m = 1, 3
              di(1,l) = di(1,l) + dri(1,l,m)*d(m)
              di(2,l) = di(2,l) + dri(2,l,m)*d(m)
              di(3,l) = di(3,l) + dri(3,l,m)*d(m)
              do n1 = 1, 3
                do n2 = 1, 3
                  qi(1,l,m) = qi(1,l,m) + 
     $                        q2(n1,n2)*(dri(1,l,n1)*rmat(m,n2) +
     $                        dri(1,m,n2)*rmat(l,n1))
                  qi(2,l,m) = qi(2,l,m) + 
     $                        q2(n1,n2)*(dri(2,l,n1)*rmat(m,n2) +
     $                        dri(2,m,n2)*rmat(l,n1))
                  qi(3,l,m) = qi(3,l,m) + 
     $                        q2(n1,n2)*(dri(3,l,n1)*rmat(m,n2) +
     $                        dri(3,m,n2)*rmat(l,n1))
                end do
              end do
            end do
          end do
c          
          scqx = 0.0d0
          scqy = 0.0d0
          scqz = 0.0d0
          do l = 1, 3
            do m = 1, 3
              scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
              scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
              scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
            end do
          end do
          dVt(1,iloc) = dVt(1,iloc) + fm*s*rr3*(di(1,1)*xr +
     $     di(1,2)*yr + di(1,3)*zr) + scqx
          dVt(2,iloc) = dVt(2,iloc) + fm*s*rr3*(di(2,1)*xr +
     $     di(2,2)*yr + di(2,3)*zr) + scqy     
          dVt(3,iloc) = dVt(3,iloc) + fm*s*rr3*(di(3,1)*xr +
     $     di(3,2)*yr + di(3,3)*zr) + scqz      
        end if
c
c       derivative with respect to the position of the iz defining atom
c
        if (doiz) then
          call aclear(9,di)
          call aclear(27,qi)
          do l = 1, 3
            do m = 1, 3
              di(1,l) = di(1,l) + driz(1,l,m)*d(m)
              di(2,l) = di(2,l) + driz(2,l,m)*d(m)
              di(3,l) = di(3,l) + driz(3,l,m)*d(m)
              do n1 = 1, 3
                do n2 = 1, 3
                  qi(1,l,m) = qi(1,l,m) + 
     $                        q2(n1,n2)*(driz(1,l,n1)*rmat(m,n2) +
     $                        driz(1,m,n2)*rmat(l,n1))
                  qi(2,l,m) = qi(2,l,m) + 
     $                        q2(n1,n2)*(driz(2,l,n1)*rmat(m,n2) +
     $                        driz(2,m,n2)*rmat(l,n1))
                  qi(3,l,m) = qi(3,l,m) + 
     $                        q2(n1,n2)*(driz(3,l,n1)*rmat(m,n2) +
     $                        driz(3,m,n2)*rmat(l,n1))
                end do
              end do
            end do
          end do
c          
          scqx = 0.0d0
          scqy = 0.0d0
          scqz = 0.0d0
          do l = 1, 3
            do m = 1, 3
              scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
              scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
              scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
            end do
          end do
c          
          dVt(1,izloc) = dVt(1,izloc) + fm*s*rr3*(di(1,1)*xr +
     $     di(1,2)*yr + di(1,3)*zr) + scqx    
          dVt(2,izloc) = dVt(2,izloc) + fm*s*rr3*(di(2,1)*xr +
     $     di(2,2)*yr + di(2,3)*zr) + scqy      
          dVt(3,izloc) = dVt(3,izloc) + fm*s*rr3*(di(3,1)*xr +
     $     di(3,2)*yr + di(3,3)*zr) + scqz      
        end if
c
c        derivative with respect to the position of the ix defining atom
c
        if (doix) then
          call aclear(9,di)
          call aclear(27,qi)
          do l = 1, 3
            do m = 1, 3
              di(1,l) = di(1,l) + drix(1,l,m)*d(m)
              di(2,l) = di(2,l) + drix(2,l,m)*d(m)
              di(3,l) = di(3,l) + drix(3,l,m)*d(m)
              do n1 = 1, 3
                do n2 = 1, 3
                  qi(1,l,m) = qi(1,l,m) + 
     $                        q2(n1,n2)*(drix(1,l,n1)*rmat(m,n2) +
     $                        drix(1,m,n2)*rmat(l,n1))
                  qi(2,l,m) = qi(2,l,m) + 
     $                        q2(n1,n2)*(drix(2,l,n1)*rmat(m,n2) +
     $                        drix(2,m,n2)*rmat(l,n1))
                  qi(3,l,m) = qi(3,l,m) + 
     $                        q2(n1,n2)*(drix(3,l,n1)*rmat(m,n2) +
     $                        drix(3,m,n2)*rmat(l,n1))
                end do
              end do
            end do
          end do
c          
          scqx = 0.0d0
          scqy = 0.0d0
          scqz = 0.0d0
          do l = 1, 3
            do m = 1, 3
              scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
              scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
              scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
            end do
          end do
c          
          dVt(1,ixloc) = dVt(1,ixloc) + fm*s*rr3*(di(1,1)*xr +
     $     di(1,2)*yr + di(1,3)*zr) + scqx      
          dVt(2,ixloc) = dVt(2,ixloc) + fm*s*rr3*(di(2,1)*xr +
     $     di(2,2)*yr + di(2,3)*zr) + scqy      
          dVt(3,ixloc) = dVt(3,ixloc) + fm*s*rr3*(di(3,1)*xr +
     $     di(3,2)*yr + di(3,3)*zr) + scqz      
        end if
c
c
c        derivative with respect to the position of the iy defining atom
c
        if (doiy) then
          call aclear(9,di)
          call aclear(27,qi)
          do l = 1, 3
            do m = 1, 3
              di(1,l) = di(1,l) + driy(1,l,m)*d(m)
              di(2,l) = di(2,l) + driy(2,l,m)*d(m)
              di(3,l) = di(3,l) + driy(3,l,m)*d(m)
              do n1 = 1, 3
                do n2 = 1, 3
                  qi(1,l,m) = qi(1,l,m) + 
     $                        q2(n1,n2)*(driy(1,l,n1)*rmat(m,n2) +
     $                        driy(1,m,n2)*rmat(l,n1))
                  qi(2,l,m) = qi(2,l,m) + 
     $                        q2(n1,n2)*(driy(2,l,n1)*rmat(m,n2) +
     $                        driy(2,m,n2)*rmat(l,n1))
                  qi(3,l,m) = qi(3,l,m) + 
     $                        q2(n1,n2)*(driy(3,l,n1)*rmat(m,n2) +
     $                        driy(3,m,n2)*rmat(l,n1))
                end do
              end do
            end do
          end do
c          
          scqx = 0.0d0
          scqy = 0.0d0
          scqz = 0.0d0
          do l = 1, 3
            do m = 1, 3
              scqx = scqx + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(1,l,m)
              scqy = scqy + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(2,l,m)
              scqz = scqz + 0.5d0*fm*s*rr5*dx(l)*dx(m)*qi(3,l,m)
            end do
          end do
c          
          dVt(1,iyloc) = dVt(1,iyloc) + fm*s*rr3*(di(1,1)*xr +
     $     di(1,2)*yr + di(1,3)*zr) + scqx      
          dVt(2,iyloc) = dVt(2,iyloc) + fm*s*rr3*(di(2,1)*xr +
     $     di(2,2)*yr + di(2,3)*zr) + scqy      
          dVt(3,iyloc) = dVt(3,iyloc) + fm*s*rr3*(di(3,1)*xr +
     $     di(3,2)*yr + di(3,3)*zr) + scqz      
        end if
      end do
c
c     sum to get the actual forces
c
      do i = 1, nbloc
        dVt(1,i) = dVtott(1,i)-dVt(1,i)
        dVt(2,i) = dVtott(2,i)-dVt(2,i)
        dVt(3,i) = dVtott(3,i)-dVt(3,i)
      end do
      return 
      end 
