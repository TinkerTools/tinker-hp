c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine dispersion1  --  dispersion derivatives             ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "dispersion3" calculates the repulsion energy
c
c
      subroutine edispersion1
      use bound
      implicit none
      if (use_replica) then
        call edispersion1a
      else
        call edispersion1b
      end if
      return
      end
c
c     "dispersion1b" calculates the dispersion energy and gradients
c     
      subroutine edispersion1b
      use action
      use analyz
      use atmtyp
      use atoms
      use atmlst
      use bond
      use bound
      use charge
      use chargetransfer
      use cutoff
      use deriv
      use dispersion
      use domdec
      use energi
      use inter
      use iounit
      use molcul
      use mpole
      use neigh
      use shunt
      use usage
      use potent
      use repulsion
      implicit none
c
      integer i,j,k,iglob,jglob,ilp,jlp
      integer iloc,jloc,kloc,ixlploc,izlploc,lloc
      integer jj,ilpnl,inl,jnl
      real*8 xi,yi,zi,xj,yj,zj,rij2,rij
      real*8 xij,yij,zij
      real*8 rdispi,rdispj,mij,mij6,mij8,mij10
      real*8 bij,rdpi,rdpj,dij
      real*8 rij3,rij4,rij5
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
      real*8 clp,gdk,tkj
      integer l
      real*8 xl,yl,zl
      real*8 xkl,ykl,zkl,rkl2,rkl
      real*8 bkl,rdpl,dkl,rskl
      real*8 mkl,mkl6,mkl8,mkl10
      real*8 ck,cl,gdl,tkl
      real*8 rlp
      real*8 fac1,clp1,clp2
      real*8 ddij(3,2),dmij(3,2)
      real*8 drkj(3,4),ddkj(3,4),dmkj(3,4)
      real*8 drij(3,6),ddkl(3,6),dmkl(3,6)
      real*8 di(3,3),dix(3,3),diz(3,3)
      real*8 di2(3,3),dix2(3,3),diz2(3,3)
      real*8 dlp1,dlp2,dlp3
      real*8 forcetemp(3)
      real*8 xr,yr,zr
      real*8 e,de,taper,dtaper
      real*8 dedx,dedy,dedz
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

      nexdisp = 0
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
         iloc = loc(iglob)
         inl = locnl(iglob)
         if (inl.eq.0) cycle
         xi = x(iglob) 
         yi = y(iglob) 
         zi = z(iglob) 
         do j = 1, natatlst(inl)
            jglob = atatlst(j,inl)
            jloc = loc(jglob)
            if (molcule(iglob) .eq. molcule(jglob)) cycle
            ddij = 0d0
            dmij = 0d0
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
            rij = sqrt(rij2) 
c
c     use energy switching if near the cutoff distance
c
            taper = 1.0d0
            dtaper = 0.0d0
            if (rij2 .gt. cut2) then
               rij3 = rij2 * rij
               rij4 = rij2 * rij2
               rij5 = rij2 * rij3
               taper = c5*rij5 + c4*rij4 + c3*rij3
     &                    + c2*rij2 + c1*rij + c0
               dtaper = 5.0d0*c5*rij4 + 4.0d0*c4*rij3
     &                     + 3.0d0*c3*rij2 + 2.0d0*c2*rij + c1
            end if
c
            bij =  gorbrep(iglob)*gorbrep(jglob)
            rdpi = vdwdisp1(iglob)
            rdpj = vdwdisp1(jglob)
            dij = ((rdpi+rdpj)*bdmp/rij)-1
c
            ddij(1,1) = (rdpi+rdpj)*bdmp*xij/rij**3 
            ddij(2,1) = (rdpi+rdpj)*bdmp*yij/rij**3 
            ddij(3,1) = (rdpi+rdpj)*bdmp*zij/rij**3 
            ddij(1,2) = -(rdpi+rdpj)*bdmp*xij/rij**3 
            ddij(2,2) = -(rdpi+rdpj)*bdmp*yij/rij**3 
            ddij(3,2) = -(rdpi+rdpj)*bdmp*zij/rij**3 
c
            dmij(1,1) = -xij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(2,1) = -yij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(3,1) = -zij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(1,2) = xij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(2,2) = yij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(3,2) = zij/(2*sqrt(rdpi*rdpj)*rij) 
c
            rsij = (rdpi + rdpj)*bdmp
            expo6 = -admp6*dij
            mij = rij/(2*sqrt(rdpi*rdpj))
            mij6 = mij*mij*mij*mij*mij*mij 
            mij8 = mij6*mij*mij
            mij10 = mij8*mij*mij

            expo6 = admp6*dij 
c           si dump pas assez efficace , > 1
            if (rij > rsij) then
              dmp6 = bij
              dmp8 = bij
              dmp10 = bij
              disp6 = c6disp*dmp6/mij6
              disp8 = c8disp*dmp8/mij8
              disp10 = c10disp*dmp10/mij10
            else
              dmp6 = bij*exp(-expo6)
              dmp8 = bij*exp(-admp8*dij)
              dmp10 = bij*exp(-admp10*dij)
              disp6 = c6disp*dmp6/mij6
              disp8 = c8disp*dmp8/mij8
              disp10 = c10disp*dmp10/mij10
c
              forcetemp(1) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(1,1)/mij6
     $       - c8disp*admp8*dmp8*ddij(1,1)/mij8
     $       - c10disp*admp10*dmp10*ddij(1,1)/mij10)*discof
              forcetemp(2) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(2,1)/mij6
     $       - c8disp*admp8*dmp8*ddij(2,1)/mij8
     $       - c10disp*admp10*dmp10*ddij(2,1)/mij10)*discof
              forcetemp(3) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(3,1)/mij6
     $       - c8disp*admp8*dmp8*ddij(3,1)/mij8
     $       - c10disp*admp10*dmp10*ddij(3,1)/mij10)*discof

              dexdisp(1,iloc) = dexdisp(1,iloc)+taper*forcetemp(1)
              dexdisp(2,iloc) = dexdisp(2,iloc)+taper*forcetemp(2)
              dexdisp(3,iloc) = dexdisp(3,iloc)+taper*forcetemp(3)
c
              forcetemp(1) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(1,2)/mij6
     $       - c8disp*admp8*dmp8*ddij(1,2)/mij8
     $       - c10disp*admp10*dmp10*ddij(1,2)/mij10)*discof
              forcetemp(2) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(2,2)/mij6
     $       - c8disp*admp8*dmp8*ddij(2,2)/mij8
     $       - c10disp*admp10*dmp10*ddij(2,2)/mij10)*discof
              forcetemp(3) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(3,2)/mij6
     $       - c8disp*admp8*dmp8*ddij(3,2)/mij8
     $       - c10disp*admp10*dmp10*ddij(3,2)/mij10)*discof

              dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1)
              dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2)
              dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3)
            end if
c
            forcetemp(1) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(1,1)/(mij**7)
     $         - c8disp*dmp8*8*dmij(1,1)/(mij**9)
     $         - c10disp*dmp10*10*dmij(1,1)/(mij**11))*discof
            forcetemp(2) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(2,1)/(mij**7)
     $         - c8disp*dmp8*8*dmij(2,1)/(mij**9)
     $         - c10disp*dmp10*10*dmij(2,1)/(mij**11))*discof
            forcetemp(3) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(3,1)/(mij**7)
     $         - c8disp*dmp8*8*dmij(3,1)/(mij**9)
     $         - c10disp*dmp10*10*dmij(3,1)/(mij**11))*discof
c
            dexdisp(1,iloc) = dexdisp(1,iloc) + taper*forcetemp(1)
            dexdisp(2,iloc) = dexdisp(2,iloc) + taper*forcetemp(2)
            dexdisp(3,iloc) = dexdisp(3,iloc) + taper*forcetemp(3)
c
            forcetemp(1) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(1,2)/(mij**7)
     $         - c8disp*dmp8*8*dmij(1,2)/(mij**9)
     $         - c10disp*dmp10*10*dmij(1,2)/(mij**11))*discof
            forcetemp(2) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(2,2)/(mij**7)
     $         - c8disp*dmp8*8*dmij(2,2)/(mij**9)
     $         - c10disp*dmp10*10*dmij(2,2)/(mij**11))*discof
            forcetemp(3) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(3,2)/(mij**7)
     $         - c8disp*dmp8*8*dmij(3,2)/(mij**9)
     $         - c10disp*dmp10*10*dmij(3,2)/(mij**11))*discof
c
            dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
            dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
            dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
c           exchange disp              
c
            ci = rpole(1,iglob) 
            cj = rpole(1,jglob)
            gdi = 1 - ci/valemtp(atomic(iglob)) 
            gdj = 1 - cj/valemtp(atomic(jglob)) 
            tij = bij*gdi*gdj
            xdisp1 = tij*cxd*exp(-axd*mij)
c
            forcetemp(1) = - facdispij*xdisp1*axd*dmij(1,1) 
            forcetemp(2) = - facdispij*xdisp1*axd*dmij(2,1) 
            forcetemp(3) = - facdispij*xdisp1*axd*dmij(3,1) 

            dexdisp(1,iloc) =dexdisp(1,iloc)+forcetemp(1)
            dexdisp(2,iloc) =dexdisp(2,iloc)+forcetemp(2)
            dexdisp(3,iloc) =dexdisp(3,iloc)+forcetemp(3)
c
            forcetemp(1) = - facdispij*xdisp1*axd*dmij(1,2) 
            forcetemp(2) = - facdispij*xdisp1*axd*dmij(2,2) 
            forcetemp(3) = - facdispij*xdisp1*axd*dmij(3,2) 
c
            dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
            dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
            dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 

            e = facdispij*scdp*(disp6 + disp8 + disp10)*discof
     $       + xdisp1*facdispij
c
c           derivatives of the switching function
c 
            dedx = -e*dtaper*xij/rij
            dedy = -e*dtaper*yij/rij
            dedz = -e*dtaper*zij/rij
            dexdisp(1,iloc) = dexdisp(1,iloc) + dedx 
            dexdisp(2,iloc) = dexdisp(2,iloc) + dedy 
            dexdisp(3,iloc) = dexdisp(3,iloc) + dedz 
            dexdisp(1,jloc) = dexdisp(1,jloc) - dedx 
            dexdisp(2,jloc) = dexdisp(2,jloc) - dedy 
            dexdisp(3,jloc) = dexdisp(3,jloc) - dedz 
c
            e = e*taper
            exdisp = exdisp + e
            nexdisp = nexdisp + 1
         end do
      end do  
c
c    lp-atom dispersion
c
      do ilpnl = 1, nlplocnl
         ilp = lpglobnl(ilpnl)
         clp = lpcharge(ilp)!2.0d0
         if ((ilpnl).eq.0) cycle
         k = lpatom(ilp)
         xlp = rlonepair(1,ilp)
         ylp = rlonepair(2,ilp)
         zlp = rlonepair(3,ilp)
         call rotlp1(ilp,di,dix,diz)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do jj = 1, nlpatlst(ilpnl)
            j = lpatlst(jj,ilpnl)
            jloc = loc(j)
            if (molcule(k) .eq. molcule(j)) cycle
               drkj = 0d0
               ddkj = 0d0
               dmkj = 0d0
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
               rkj = sqrt(rkj2)
c
c     use energy switching if near the cutoff distance
c
               taper = 1.0d0
               dtaper = 0.0d0
               if (rkj2 .gt. cut2) then
                  rkj3 = rkj2 * rkj
                  rkj4 = rkj2 * rkj2
                  rkj5 = rkj2 * rkj3
                  taper = c5*rkj5 + c4*rkj4 + c3*rkj3
     &                       + c2*rkj2 + c1*rkj + c0
                  dtaper = 5.0d0*c5*rkj4 + 4.0d0*c4*rkj3
     &                        + 3.0d0*c3*rkj2 + 2.0d0*c2*rkj + c1
               end if
c
c
               drkj(1,1) = xkj/rkj
               drkj(2,1) = ykj/rkj
               drkj(3,1) = zkj/rkj
c
               dlp1 = -xkj/rkj
               dlp2 = -ykj/rkj
               dlp3 = -zkj/rkj
               drkj(1,2) =  
     $          di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
               drkj(2,2) =  
     $          di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
               drkj(3,2) =  
     $          di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3
      
               drkj(1,3) =  
     $          dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
               drkj(2,3) =  
     $          dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
               drkj(3,3) =  
     $          dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3
               drkj(1,4) =  
     $          diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
               drkj(2,4) =  
     $          diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
               drkj(3,4) =  
     $          diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3
c
               bkj = gorbrep(j)
               bkjx = gorbrep(k)*gorbrep(j)
c
c   vdw radius is the vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp2(k)! + dincr_lprep(ilp)
               rdpj = vdwdisp1(j)
               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
c
               ddkj(1,1) = -(rdpk+rdpj)*bdmp*drkj(1,1)/rkj**2 
               ddkj(2,1) = -(rdpk+rdpj)*bdmp*drkj(2,1)/rkj**2 
               ddkj(3,1) = -(rdpk+rdpj)*bdmp*drkj(3,1)/rkj**2 

               ddkj(1,2) = -(rdpk+rdpj)*bdmp*drkj(1,2)/rkj**2 
               ddkj(2,2) = -(rdpk+rdpj)*bdmp*drkj(2,2)/rkj**2 
               ddkj(3,2) = -(rdpk+rdpj)*bdmp*drkj(3,2)/rkj**2 

               ddkj(1,3) = -(rdpk+rdpj)*bdmp*drkj(1,3)/rkj**2 
               ddkj(2,3) = -(rdpk+rdpj)*bdmp*drkj(2,3)/rkj**2 
               ddkj(3,3) = -(rdpk+rdpj)*bdmp*drkj(3,3)/rkj**2 

               ddkj(1,4) = -(rdpk+rdpj)*bdmp*drkj(1,4)/rkj**2 
               ddkj(2,4) = -(rdpk+rdpj)*bdmp*drkj(2,4)/rkj**2 
               ddkj(3,4) = -(rdpk+rdpj)*bdmp*drkj(3,4)/rkj**2 
c
               dmkj(1,1) = drkj(1,1)/(2*sqrt(rdpk*rdpj))
               dmkj(2,1) = drkj(2,1)/(2*sqrt(rdpk*rdpj))
               dmkj(3,1) = drkj(3,1)/(2*sqrt(rdpk*rdpj))

               dmkj(1,2) = drkj(1,2)/(2*sqrt(rdpk*rdpj))
               dmkj(2,2) = drkj(2,2)/(2*sqrt(rdpk*rdpj))
               dmkj(3,2) = drkj(3,2)/(2*sqrt(rdpk*rdpj))

               dmkj(1,3) = drkj(1,3)/(2*sqrt(rdpk*rdpj))
               dmkj(2,3) = drkj(2,3)/(2*sqrt(rdpk*rdpj))
               dmkj(3,3) = drkj(3,3)/(2*sqrt(rdpk*rdpj))

               dmkj(1,4) = drkj(1,4)/(2*sqrt(rdpk*rdpj))
               dmkj(2,4) = drkj(2,4)/(2*sqrt(rdpk*rdpj))
               dmkj(3,4) = drkj(3,4)/(2*sqrt(rdpk*rdpj))
c
               rskj = (rdpk + rdpj)*bdmp
               expo6 = -admp6*dkj
               mkj = rkj/(2*sqrt(rdpk*rdpj))
               mkj6 = mkj*mkj*mkj*mkj*mkj*mkj 
               mkj8 = mkj6*mkj*mkj
               mkj10 = mkj8*mkj*mkj
               fac1 = -0.5*clp*colpa*facdispij*discof
               if (rkj > rskj) then
                 expo6=0.0d0
                 dmp6 = bkj*exp(expo6)
                 dmp8 = bkj*exp(-admp8*dkj)
                 dmp10 = bkj*exp(-admp10*dkj)
               else
                 dmp6 = bkj*exp(expo6)
                 dmp8 = bkj*exp(-admp8*dkj)
                 dmp10 = bkj*exp(-admp10*dkj)
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,1)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,1)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,1)/mkj6
c
                 jloc = loc(j)
                 dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
                 dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
                 dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,2)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,2)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,2)/mkj6
c
                 kloc = loc(k)
                 dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
                 dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2)
                 dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,3)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,3)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,3)/mkj6
c
                 ixlploc = loc(ixlp(ilp))
                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $             taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $             taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $             taper*forcetemp(3)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,4)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,4)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,4)/mkj6
c
                 izlploc = loc(izlp(ilp))
                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $             taper*forcetemp(1) 
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $             taper*forcetemp(2) 
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $             taper*forcetemp(3) 
               end if
               disp6 = fac1*c6disp*dmp6/mkj6
               disp8 = fac1*c8disp*dmp8/mkj8
               disp10 = fac1*c10disp*dmp10/mkj10
c
               forcetemp(1) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(1,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,1)/mkj8
     $          - dmp8*8*dmkj(1,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,1)/mkj10
     $            - dmp10*10*dmkj(1,1)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(2,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,1)/mkj8
     $          - dmp8*8*dmkj(2,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,1)/mkj10
     $            - dmp10*10*dmkj(2,1)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(3,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,1)/mkj8
     $          - dmp8*8*dmkj(3,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,1)/mkj10
     $            - dmp10*10*dmkj(3,1)/(mkj**11)))
c
               jloc = loc(j)
               dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1)
               dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2)
               dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3)
c
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(1,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,2)/mkj8
     $          - dmp8*8*dmkj(1,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,2)/mkj10
     $            - dmp10*10*dmkj(1,2)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(2,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,2)/mkj8
     $          - dmp8*8*dmkj(2,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,2)/mkj10
     $            - dmp10*10*dmkj(2,2)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(3,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,2)/mkj8
     $          - dmp8*8*dmkj(3,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,2)/mkj10
     $            - dmp10*10*dmkj(3,2)/(mkj**11)))
c
               kloc = loc(k)
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(1,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,3)/mkj8
     $          - dmp8*8*dmkj(1,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,3)/mkj10
     $            - dmp10*10*dmkj(1,3)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(2,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,3)/mkj8
     $          - dmp8*8*dmkj(2,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,3)/mkj10
     $            - dmp10*10*dmkj(2,3)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(3,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,3)/mkj8
     $          - dmp8*8*dmkj(3,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,3)/mkj10
     $            - dmp10*10*dmkj(3,3)/(mkj**11)))
c
               ixlploc = loc(ixlp(ilp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           taper*forcetemp(3) 
c
               forcetemp(1) =  fac1*(
     $          c6disp*(- dmp6*6*dmkj(1,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,4)/mkj8
     $          - dmp8*8*dmkj(1,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,4)/mkj10
     $            - dmp10*10*dmkj(1,4)/(mkj**11)))
               forcetemp(2) =  fac1*(
     $         c6disp* (- dmp6*6*dmkj(2,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,4)/mkj8
     $          - dmp8*8*dmkj(2,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,4)/mkj10
     $            - dmp10*10*dmkj(2,4)/(mkj**11)))
               forcetemp(3) =  fac1*(
     $          +c6disp*(- dmp6*6*dmkj(3,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,4)/mkj8
     $          - dmp8*8*dmkj(3,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,4)/mkj10
     $            - dmp10*10*dmkj(3,4)/(mkj**11)))
c
               izlploc = loc(izlp(ilp))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
c              bkj =  frep(atomic(k))*frep(atomic(j))
               bkjx = gorbrep(k)*gorbrep(j)
c               rdpk = rlp
               rdpj = vdwdisp1(j)
c               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
c
c              exchange disp
c
               ck = rpole(1,k)
               cj = rpole(1,j)
               gdk = 1 - ck/valemtp(atomic(k)) 
               gdj = 1 - cj/valemtp(atomic(j)) 
               tkj = bkjx*gdk*gdj
c
               xdisp2 = facdispij*0.5*clp*tkj*cxdla*exp(-axdla*mkj)
cc
               forcetemp(1) = -xdisp2*axdla*dmkj(1,1) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,1) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,1) 
               jloc = loc(j)
               dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
               dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
               dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
               forcetemp(1) = -xdisp2*axdla*dmkj(1,2) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,2) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,2) 
c
               kloc = loc(k)
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               forcetemp(1) =  -xdisp2*axdla*dmkj(1,3) 
               forcetemp(2) =  -xdisp2*axdla*dmkj(2,3) 
               forcetemp(3) =  -xdisp2*axdla*dmkj(3,3) 

               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           taper*forcetemp(3) 
c
               forcetemp(1) = -xdisp2*axdla*dmkj(1,4) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,4) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,4) 

               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $           taper*forcetemp(3) 
               e = (disp6 + disp8 + disp10) + xdisp2
c
c           derivatives of the switching function
c 
               dexdisp(1,jloc) = dexdisp(1,jloc) + e*dtaper*drkj(1,1) 
               dexdisp(2,jloc) = dexdisp(2,jloc) + e*dtaper*drkj(2,1) 
               dexdisp(3,jloc) = dexdisp(3,jloc) + e*dtaper*drkj(3,1) 
               dexdisp(1,kloc) = dexdisp(1,kloc) + e*dtaper*drkj(1,2) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + e*dtaper*drkj(2,2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + e*dtaper*drkj(3,2) 
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $            e*dtaper*drkj(1,3) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $            e*dtaper*drkj(2,3) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $            e*dtaper*drkj(3,3) 
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $            e*dtaper*drkj(1,4) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $            e*dtaper*drkj(2,4) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $            e*dtaper*drkj(3,4) 

c               
               e = e*taper
               exdisp = exdisp + e
         end do 
      end do  
c
c      lone pair / lone pair dispersion term 
c
      do ilpnl = 1, nlplocnl
         ilp = lpglobnl(ilpnl)
         if (ilpnl.eq.0) cycle
         call rotlp1(ilp,di,dix,diz)
         k = lpatom(ilp)    
         clp1 = lpcharge(ilp)
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, nlplplst(ilpnl)
            jlp = lplplst(j,ilpnl)
            call rotlp1(jlp,di2,dix2,diz2)
            drij = 0d0
            ddkl = 0d0
            dmkl = 0d0
            l = lpatom(jlp)
            clp2 = lpcharge(jlp)
            if (molcule(k) .eq. molcule(l)) cycle
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
c
c     use energy switching if near the cutoff distance
c
               taper = 1.0d0
               dtaper = 0.0d0
               if (rij2 .gt. cut2) then
                  rij3 = rij2 * rij
                  rij4 = rij2 * rij2
                  rij5 = rij2 * rij3
                  taper = c5*rij5 + c4*rij4 + c3*rij3
     &                       + c2*rij2 + c1*rij + c0
                  dtaper = 5.0d0*c5*rij4 + 4.0d0*c4*rij3
     &                        + 3.0d0*c3*rij2 + 2.0d0*c2*rij + c1
               end if

               bkl = 1  
c
c
c    vdw radius is vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp3(k)! + dincr_lprep(ilp) 
               rdpl = vdwdisp3(l)! + dincr_lprep(jlp) 
c
               dkl = ((rdpk+rdpl)*bdmp/rij)-1
               rskl = (rdpk + rdpl)*bdmp
               mkl = rij/(2*sqrt(rdpk*rdpl))
               mkl6 = mkl*mkl*mkl*mkl*mkl*mkl 
               mkl8 = mkl6*mkl*mkl
               mkl10 = mkl8*mkl*mkl
c
               dlp1 = -xij/rij
               dlp2 = -yij/rij
               dlp3 = -zij/rij

               drij(1,1) =  
     $          di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
               drij(2,1) =  
     $          di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
               drij(3,1) =  
     $          di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3

               drij(1,2) =  
     $          dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
               drij(2,2) =  
     $          dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
               drij(3,2) =  
     $          dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3

               drij(1,3) =  
     $          diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
               drij(2,3) =  
     $          diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
               drij(3,3) =  
     $          diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3

               drij(1,4) =  
     $          -di2(1,1)*dlp1 - di2(1,2)*dlp2 - di2(1,3)*dlp3
               drij(2,4) =  
     $          -di2(2,1)*dlp1 - di2(2,2)*dlp2 - di2(2,3)*dlp3
               drij(3,4) =  
     $          -di2(3,1)*dlp1 - di2(3,2)*dlp2 - di2(3,3)*dlp3

               drij(1,5) =  
     $          -dix2(1,1)*dlp1 - dix2(1,2)*dlp2 - dix2(1,3)*dlp3
               drij(2,5) =  
     $          -dix2(2,1)*dlp1 - dix2(2,2)*dlp2 - dix2(2,3)*dlp3
               drij(3,5) =  
     $          -dix2(3,1)*dlp1 - dix2(3,2)*dlp2 - dix2(3,3)*dlp3

               drij(1,6) =  
     $          -diz2(1,1)*dlp1 - diz2(1,2)*dlp2 - diz2(1,3)*dlp3
               drij(2,6) =  
     $          -diz2(2,1)*dlp1 - diz2(2,2)*dlp2 - diz2(2,3)*dlp3
               drij(3,6) =  
     $          -diz2(3,1)*dlp1 - diz2(3,2)*dlp2 - diz2(3,3)*dlp3
      
               ddkl(1,1) = -(rdpk+rdpl)*bdmp*drij(1,1)/rij**2 
               ddkl(2,1) = -(rdpk+rdpl)*bdmp*drij(2,1)/rij**2 
               ddkl(3,1) = -(rdpk+rdpl)*bdmp*drij(3,1)/rij**2 
               ddkl(1,2) = -(rdpk+rdpl)*bdmp*drij(1,2)/rij**2 
               ddkl(2,2) = -(rdpk+rdpl)*bdmp*drij(2,2)/rij**2 
               ddkl(3,2) = -(rdpk+rdpl)*bdmp*drij(3,2)/rij**2 
               ddkl(1,3) = -(rdpk+rdpl)*bdmp*drij(1,3)/rij**2 
               ddkl(2,3) = -(rdpk+rdpl)*bdmp*drij(2,3)/rij**2 
               ddkl(3,3) = -(rdpk+rdpl)*bdmp*drij(3,3)/rij**2 
               ddkl(1,4) = -(rdpk+rdpl)*bdmp*drij(1,4)/rij**2 
               ddkl(2,4) = -(rdpk+rdpl)*bdmp*drij(2,4)/rij**2 
               ddkl(3,4) = -(rdpk+rdpl)*bdmp*drij(3,4)/rij**2 
               ddkl(1,5) = -(rdpk+rdpl)*bdmp*drij(1,5)/rij**2 
               ddkl(2,5) = -(rdpk+rdpl)*bdmp*drij(2,5)/rij**2 
               ddkl(3,5) = -(rdpk+rdpl)*bdmp*drij(3,5)/rij**2 
               ddkl(1,6) = -(rdpk+rdpl)*bdmp*drij(1,6)/rij**2 
               ddkl(2,6) = -(rdpk+rdpl)*bdmp*drij(2,6)/rij**2 
               ddkl(3,6) = -(rdpk+rdpl)*bdmp*drij(3,6)/rij**2 
      
               dmkl(1,1) = drij(1,1)/(2*sqrt(rdpk*rdpl))
               dmkl(2,1) = drij(2,1)/(2*sqrt(rdpk*rdpl))
               dmkl(3,1) = drij(3,1)/(2*sqrt(rdpk*rdpl))
               dmkl(1,2) = drij(1,2)/(2*sqrt(rdpk*rdpl))
               dmkl(2,2) = drij(2,2)/(2*sqrt(rdpk*rdpl))
               dmkl(3,2) = drij(3,2)/(2*sqrt(rdpk*rdpl))
               dmkl(1,3) = drij(1,3)/(2*sqrt(rdpk*rdpl))
               dmkl(2,3) = drij(2,3)/(2*sqrt(rdpk*rdpl))
               dmkl(3,3) = drij(3,3)/(2*sqrt(rdpk*rdpl))
               dmkl(1,4) = drij(1,4)/(2*sqrt(rdpk*rdpl))
               dmkl(2,4) = drij(2,4)/(2*sqrt(rdpk*rdpl))
               dmkl(3,4) = drij(3,4)/(2*sqrt(rdpk*rdpl))
               dmkl(1,5) = drij(1,5)/(2*sqrt(rdpk*rdpl))
               dmkl(2,5) = drij(2,5)/(2*sqrt(rdpk*rdpl))
               dmkl(3,5) = drij(3,5)/(2*sqrt(rdpk*rdpl))
               dmkl(1,6) = drij(1,6)/(2*sqrt(rdpk*rdpl))
               dmkl(2,6) = drij(2,6)/(2*sqrt(rdpk*rdpl))
               dmkl(3,6) = drij(3,6)/(2*sqrt(rdpk*rdpl))
c

               expo6 = -admp6*dkl
               fac1 = -0.5*clp1*0.5*clp2*facdispij*colp*discof
               if (rij > rskl) then
                 expo6 = 0.0d0 
                 dmp6 = bkl*exp(expo6)
                 dmp8 = bkl*exp(-admp8*dkl)
                 dmp10 = bkl*exp(-admp10*dkl)
               else
                 dmp6 = bkl*exp(expo6)
                 dmp8 = bkl*exp(-admp8*dkl)
                 dmp10 = bkl*exp(-admp10*dkl)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,1)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,1)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,1)/mkl6
                 kloc = loc(k)
                 dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
                 dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
                 dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,2)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,2)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,2)/mkl6
                 ixlploc = loc(ixlp(ilp))
                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $            taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $            taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $            taper*forcetemp(3)
c
                 izlploc = loc(izlp(ilp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,3)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,3)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,3)/mkl6

                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $            taper*forcetemp(1)
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $            taper*forcetemp(2)
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $            taper*forcetemp(3)
c
                 lloc = loc(l)
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,4)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,4)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,4)/mkl6

                 dexdisp(1,lloc) = dexdisp(1,lloc) + taper*forcetemp(1) 
                 dexdisp(2,lloc) = dexdisp(2,lloc) + taper*forcetemp(2) 
                 dexdisp(3,lloc) = dexdisp(3,lloc) + taper*forcetemp(3) 
c
                 ixlploc = loc(ixlp(jlp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,5)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,5)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,5)/mkl6

                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $              taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $              taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $              taper*forcetemp(3)
c
                 izlploc = loc(izlp(jlp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,6)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,6)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,6)/mkl6

                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $             taper*forcetemp(1)
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $             taper*forcetemp(2)
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $             taper*forcetemp(3)
               end if
c
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(1,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,1)/mkl8
     $          - dmp8*8*dmkl(1,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,1)/mkl10
     $            - dmp10*10*dmkl(1,1)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,1)/mkl8
     $          - dmp8*8*dmkl(2,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,1)/mkl10
     $            - dmp10*10*dmkl(2,1)/(mkl**11)))
               forcetemp(3) = fac1*(
     $         c6disp*(- dmp6*6*dmkl(3,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,1)/mkl8
     $          - dmp8*8*dmkl(3,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,1)/mkl10
     $            - dmp10*10*dmkl(3,1)/(mkl**11)))
               kloc = loc(k)

               dexdisp(1,kloc) = dexdisp(1,kloc) +taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) +taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) +taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(ilp))
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(1,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,2)/mkl8
     $          - dmp8*8*dmkl(1,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,2)/mkl10
     $            - dmp10*10*dmkl(1,2)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(2,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,2)/mkl8
     $          - dmp8*8*dmkl(2,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,2)/mkl10
     $            - dmp10*10*dmkl(2,2)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(3,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,2)/mkl8
     $          - dmp8*8*dmkl(3,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,2)/mkl10
     $            - dmp10*10*dmkl(3,2)/(mkl**11)))

               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(ilp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,3)/mkl8
     $          - dmp8*8*dmkl(1,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,3)/mkl10
     $            - dmp10*10*dmkl(1,3)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,3)/mkl8
     $          - dmp8*8*dmkl(2,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,3)/mkl10
     $            - dmp10*10*dmkl(2,3)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,3)/mkl8
     $          - dmp8*8*dmkl(3,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,3)/mkl10
     $            - dmp10*10*dmkl(3,3)/(mkl**11)))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
               lloc = loc(l)
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,4)/mkl8
     $          - dmp8*8*dmkl(1,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,4)/mkl10
     $            - dmp10*10*dmkl(1,4)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,4)/mkl8
     $          - dmp8*8*dmkl(2,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,4)/mkl10
     $            - dmp10*10*dmkl(2,4)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,4)/mkl8
     $          - dmp8*8*dmkl(3,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,4)/mkl10
     $            - dmp10*10*dmkl(3,4)/(mkl**11)))
               dexdisp(1,lloc) = dexdisp(1,lloc) +taper*forcetemp(1) 
               dexdisp(2,lloc) = dexdisp(2,lloc) +taper*forcetemp(2) 
               dexdisp(3,lloc) = dexdisp(3,lloc) +taper*forcetemp(3) 

               ixlploc = loc(ixlp(jlp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,5)/mkl8
     $          - dmp8*8*dmkl(1,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,5)/mkl10
     $            - dmp10*10*dmkl(1,5)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,5)/mkl8
     $          - dmp8*8*dmkl(2,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,5)/mkl10
     $            - dmp10*10*dmkl(2,5)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,5)/mkl8
     $          - dmp8*8*dmkl(3,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,5)/mkl10
     $            - dmp10*10*dmkl(3,5)/(mkl**11)))
               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(jlp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,6)/mkl8
     $          - dmp8*8*dmkl(1,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,6)/mkl10
     $            - dmp10*10*dmkl(1,6)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,6)/mkl8
     $          - dmp8*8*dmkl(2,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,6)/mkl10
     $            - dmp10*10*dmkl(2,6)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,6)/mkl8
     $          - dmp8*8*dmkl(3,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,6)/mkl10
     $            - dmp10*10*dmkl(3,6)/(mkl**11)))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)

               disp6 = fac1*c6disp*dmp6/mkl6
               disp8 = fac1*c8disp*dmp8/mkl8
               disp10 = fac1*c10disp*dmp10/mkl10
c
c          Exchange dispersion                  
c
               bkl = gorbrep(k)*gorbrep(l)
               ck = rpole(1,k)
               cl = rpole(1,l)
               gdk = 1 - ck/valemtp(atomic(k))
               gdl = 1 - cl/valemtp(atomic(l))
               tkl = bkl*gdk*gdl
               xdisp3 = facdispij*0.5*clp1*0.5*clp2*tkl*cxdlp*
     $            exp(-axdlp*mkl) 
c
               kloc = loc(k)
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,1) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,1) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,1) 
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(ilp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,2) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,2) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,2) 
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           forcetemp(3) 
c
               izlploc = loc(izlp(ilp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,3) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,3) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,3) 

               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
               lloc = loc(l)
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,4) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,4) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,4) 
               dexdisp(1,lloc) = dexdisp(1,lloc) + taper*forcetemp(1) 
               dexdisp(2,lloc) = dexdisp(2,lloc) + taper*forcetemp(2) 
               dexdisp(3,lloc) = dexdisp(3,lloc) + taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(jlp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,5) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,5) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,5) 
               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(jlp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,6) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,6) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,6) 

               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $           taper*forcetemp(3) 
              
               e = disp6 + disp8 + disp10 + xdisp3
c
c           derivatives of the switching function
c 
               dexdisp(1,kloc) = dexdisp(1,kloc) + e*dtaper*drij(1,1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + e*dtaper*drij(2,1) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + e*dtaper*drij(3,1) 
               ixlploc = loc(ixlp(ilp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $          e*dtaper*drij(1,2) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $          e*dtaper*drij(2,2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $          e*dtaper*drij(3,2) 
               izlploc = loc(izlp(ilp))
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $          e*dtaper*drij(1,3) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $          e*dtaper*drij(2,3) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $          e*dtaper*drij(3,3) 

               dexdisp(1,lloc) = dexdisp(1,lloc) +e*dtaper*drij(1,4) 
               dexdisp(2,lloc) = dexdisp(2,lloc) +e*dtaper*drij(2,4) 
               dexdisp(3,lloc) = dexdisp(3,lloc) +e*dtaper*drij(3,4) 
               ixlploc = loc(ixlp(jlp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $          e*dtaper*drij(1,5) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $          e*dtaper*drij(2,5) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $          e*dtaper*drij(3,5) 
               izlploc = loc(izlp(jlp))
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $          e*dtaper*drij(1,6) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $          e*dtaper*drij(2,6) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $          e*dtaper*drij(3,6) 

               e = e*taper
               exdisp = exdisp + e
         end do
      end do
      einter = einter + exdisp
      return
      end
c
c     "dispersion1a" calculates the dispersion energy and gradients
c     
      subroutine edispersion1a
      use action
      use analyz
      use atmtyp
      use atoms
      use atmlst
      use bond
      use bound
      use boxes
      use cell
      use charge
      use chargetransfer
      use cutoff
      use deriv
      use dispersion
      use domdec
      use energi
      use inter
      use iounit
      use molcul
      use mpole
      use neigh
      use shunt
      use usage
      use potent
      use repulsion
      implicit none
c
      integer i,j,k,iglob,jglob,ilp,jlp,m
      integer iloc,jloc,kloc,ixlploc,izlploc,lloc
      integer jj,ilpnl,inl,jnl
      real*8 xi,yi,zi,xj,yj,zj,rij2,rij
      real*8 xij,yij,zij
      real*8 rdispi,rdispj,mij,mij6,mij8,mij10
      real*8 bij,rdpi,rdpj,dij
      real*8 rij3,rij4,rij5
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
      real*8 clp,gdk,tkj
      integer l
      real*8 xl,yl,zl
      real*8 xkl,ykl,zkl,rkl2,rkl
      real*8 bkl,rdpl,dkl,rskl
      real*8 mkl,mkl6,mkl8,mkl10
      real*8 ck,cl,gdl,tkl
      real*8 rlp
      real*8 fac1,clp1,clp2
      real*8 ddij(3,2),dmij(3,2)
      real*8 drkj(3,4),ddkj(3,4),dmkj(3,4)
      real*8 drij(3,6),ddkl(3,6),dmkl(3,6)
      real*8 di(3,3),dix(3,3),diz(3,3)
      real*8 di2(3,3),dix2(3,3),diz2(3,3)
      real*8 dlp1,dlp2,dlp3
      real*8 forcetemp(3)
      real*8 xr,yr,zr
      real*8 e,de,taper,dtaper
      real*8 dedx,dedy,dedz
      real*8 xmove,ymove,zmove
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

      nexdisp = 0
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
         iloc = loc(iglob)
         xi = x(iglob) 
         yi = y(iglob) 
         zi = z(iglob) 
         do j = 1, n
            jglob = j
            if (jglob.le.iglob) cycle
            jloc = loc(jglob)
            if (molcule(iglob) .eq. molcule(jglob)) cycle
            ddij = 0d0
            dmij = 0d0
            xj = x(jglob) 
            yj = y(jglob) 
            zj = z(jglob) 
            xij = xj - xi 
            yij = yj - yi 
            zij = zj - zi 
c            if (use_bounds)  call image (xij,yij,zij)
            rij2 = xij*xij + yij*yij + zij*zij
c
            if (rij2.gt.(dispcut*dispcut)) cycle
            rij = sqrt(rij2) 
c
c     use energy switching if near the cutoff distance
c
            taper = 1.0d0
            dtaper = 0.0d0
            if (rij2 .gt. cut2) then
               rij3 = rij2 * rij
               rij4 = rij2 * rij2
               rij5 = rij2 * rij3
               taper = c5*rij5 + c4*rij4 + c3*rij3
     &                    + c2*rij2 + c1*rij + c0
               dtaper = 5.0d0*c5*rij4 + 4.0d0*c4*rij3
     &                     + 3.0d0*c3*rij2 + 2.0d0*c2*rij + c1
            end if
c
            bij =  gorbrep(iglob)*gorbrep(jglob)
            rdpi = vdwdisp1(iglob)
            rdpj = vdwdisp1(jglob)
            dij = ((rdpi+rdpj)*bdmp/rij)-1
c
            ddij(1,1) = (rdpi+rdpj)*bdmp*xij/rij**3 
            ddij(2,1) = (rdpi+rdpj)*bdmp*yij/rij**3 
            ddij(3,1) = (rdpi+rdpj)*bdmp*zij/rij**3 
            ddij(1,2) = -(rdpi+rdpj)*bdmp*xij/rij**3 
            ddij(2,2) = -(rdpi+rdpj)*bdmp*yij/rij**3 
            ddij(3,2) = -(rdpi+rdpj)*bdmp*zij/rij**3 
c
            dmij(1,1) = -xij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(2,1) = -yij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(3,1) = -zij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(1,2) = xij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(2,2) = yij/(2*sqrt(rdpi*rdpj)*rij) 
            dmij(3,2) = zij/(2*sqrt(rdpi*rdpj)*rij) 
c
            rsij = (rdpi + rdpj)*bdmp
            expo6 = -admp6*dij
            mij = rij/(2*sqrt(rdpi*rdpj))
            mij6 = mij*mij*mij*mij*mij*mij 
            mij8 = mij6*mij*mij
            mij10 = mij8*mij*mij

            expo6 = admp6*dij 
c           si dump pas assez efficace , > 1
            if (rij > rsij) then
              dmp6 = bij
              dmp8 = bij
              dmp10 = bij
              disp6 = c6disp*dmp6/mij6
              disp8 = c8disp*dmp8/mij8
              disp10 = c10disp*dmp10/mij10
            else
              dmp6 = bij*exp(-expo6)
              dmp8 = bij*exp(-admp8*dij)
              dmp10 = bij*exp(-admp10*dij)
              disp6 = c6disp*dmp6/mij6
              disp8 = c8disp*dmp8/mij8
              disp10 = c10disp*dmp10/mij10
c
              forcetemp(1) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(1,1)/mij6
     $       - c8disp*admp8*dmp8*ddij(1,1)/mij8
     $       - c10disp*admp10*dmp10*ddij(1,1)/mij10)*discof
              forcetemp(2) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(2,1)/mij6
     $       - c8disp*admp8*dmp8*ddij(2,1)/mij8
     $       - c10disp*admp10*dmp10*ddij(2,1)/mij10)*discof
              forcetemp(3) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(3,1)/mij6
     $       - c8disp*admp8*dmp8*ddij(3,1)/mij8
     $       - c10disp*admp10*dmp10*ddij(3,1)/mij10)*discof

              dexdisp(1,iloc) = dexdisp(1,iloc)+taper*forcetemp(1)
              dexdisp(2,iloc) = dexdisp(2,iloc)+taper*forcetemp(2)
              dexdisp(3,iloc) = dexdisp(3,iloc)+taper*forcetemp(3)
c
              forcetemp(1) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(1,2)/mij6
     $       - c8disp*admp8*dmp8*ddij(1,2)/mij8
     $       - c10disp*admp10*dmp10*ddij(1,2)/mij10)*discof
              forcetemp(2) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(2,2)/mij6
     $       - c8disp*admp8*dmp8*ddij(2,2)/mij8
     $       - c10disp*admp10*dmp10*ddij(2,2)/mij10)*discof
              forcetemp(3) =facdispij*scdp*(
     $       - c6disp*admp6*dmp6*ddij(3,2)/mij6
     $       - c8disp*admp8*dmp8*ddij(3,2)/mij8
     $       - c10disp*admp10*dmp10*ddij(3,2)/mij10)*discof

              dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1)
              dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2)
              dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3)
            end if
c
            forcetemp(1) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(1,1)/(mij**7)
     $         - c8disp*dmp8*8*dmij(1,1)/(mij**9)
     $         - c10disp*dmp10*10*dmij(1,1)/(mij**11))*discof
            forcetemp(2) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(2,1)/(mij**7)
     $         - c8disp*dmp8*8*dmij(2,1)/(mij**9)
     $         - c10disp*dmp10*10*dmij(2,1)/(mij**11))*discof
            forcetemp(3) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(3,1)/(mij**7)
     $         - c8disp*dmp8*8*dmij(3,1)/(mij**9)
     $         - c10disp*dmp10*10*dmij(3,1)/(mij**11))*discof
c
            dexdisp(1,iloc) = dexdisp(1,iloc) + taper*forcetemp(1)
            dexdisp(2,iloc) = dexdisp(2,iloc) + taper*forcetemp(2)
            dexdisp(3,iloc) = dexdisp(3,iloc) + taper*forcetemp(3)
c
            forcetemp(1) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(1,2)/(mij**7)
     $         - c8disp*dmp8*8*dmij(1,2)/(mij**9)
     $         - c10disp*dmp10*10*dmij(1,2)/(mij**11))*discof
            forcetemp(2) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(2,2)/(mij**7)
     $         - c8disp*dmp8*8*dmij(2,2)/(mij**9)
     $         - c10disp*dmp10*10*dmij(2,2)/(mij**11))*discof
            forcetemp(3) = facdispij*scdp*(
     $       - c6disp*dmp6*6*dmij(3,2)/(mij**7)
     $         - c8disp*dmp8*8*dmij(3,2)/(mij**9)
     $         - c10disp*dmp10*10*dmij(3,2)/(mij**11))*discof
c
            dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
            dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
            dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
c           exchange disp              
c
            ci = rpole(1,iglob) 
            cj = rpole(1,jglob)
            gdi = 1 - ci/valemtp(atomic(iglob)) 
            gdj = 1 - cj/valemtp(atomic(jglob)) 
            tij = bij*gdi*gdj
            xdisp1 = tij*cxd*exp(-axd*mij)
c
            forcetemp(1) = - facdispij*xdisp1*axd*dmij(1,1) 
            forcetemp(2) = - facdispij*xdisp1*axd*dmij(2,1) 
            forcetemp(3) = - facdispij*xdisp1*axd*dmij(3,1) 

            dexdisp(1,iloc) =dexdisp(1,iloc)+forcetemp(1)
            dexdisp(2,iloc) =dexdisp(2,iloc)+forcetemp(2)
            dexdisp(3,iloc) =dexdisp(3,iloc)+forcetemp(3)
c
            forcetemp(1) = - facdispij*xdisp1*axd*dmij(1,2) 
            forcetemp(2) = - facdispij*xdisp1*axd*dmij(2,2) 
            forcetemp(3) = - facdispij*xdisp1*axd*dmij(3,2) 
c
            dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
            dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
            dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 

            e = facdispij*scdp*(disp6 + disp8 + disp10)*discof
     $       + xdisp1*facdispij
c
c           derivatives of the switching function
c 
            dedx = -e*dtaper*xij/rij
            dedy = -e*dtaper*yij/rij
            dedz = -e*dtaper*zij/rij
            dexdisp(1,iloc) = dexdisp(1,iloc) + dedx 
            dexdisp(2,iloc) = dexdisp(2,iloc) + dedy 
            dexdisp(3,iloc) = dexdisp(3,iloc) + dedz 
            dexdisp(1,jloc) = dexdisp(1,jloc) - dedx 
            dexdisp(2,jloc) = dexdisp(2,jloc) - dedy 
            dexdisp(3,jloc) = dexdisp(3,jloc) - dedz 
c
            e = e*taper
            exdisp = exdisp + e
            nexdisp = nexdisp + 1
         end do
      end do  
c
c    lp-atom dispersion
c
      do i = 1, nlploc
         ilp = lpglob(i)
         clp = lpcharge(ilp)!2.0d0
         k = lpatom(ilp)
         xlp = rlonepair(1,ilp)
         ylp = rlonepair(2,ilp)
         zlp = rlonepair(3,ilp)
         call rotlp1(ilp,di,dix,diz)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, n
            jloc = loc(j)
            if (molcule(k) .eq. molcule(j)) cycle
               drkj = 0d0
               ddkj = 0d0
               dmkj = 0d0
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
               rkj = sqrt(rkj2)
c
c     use energy switching if near the cutoff distance
c
               taper = 1.0d0
               dtaper = 0.0d0
               if (rkj2 .gt. cut2) then
                  rkj3 = rkj2 * rkj
                  rkj4 = rkj2 * rkj2
                  rkj5 = rkj2 * rkj3
                  taper = c5*rkj5 + c4*rkj4 + c3*rkj3
     &                       + c2*rkj2 + c1*rkj + c0
                  dtaper = 5.0d0*c5*rkj4 + 4.0d0*c4*rkj3
     &                        + 3.0d0*c3*rkj2 + 2.0d0*c2*rkj + c1
               end if
c
c
               drkj(1,1) = xkj/rkj
               drkj(2,1) = ykj/rkj
               drkj(3,1) = zkj/rkj
c
               dlp1 = -xkj/rkj
               dlp2 = -ykj/rkj
               dlp3 = -zkj/rkj
               drkj(1,2) =  
     $          di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
               drkj(2,2) =  
     $          di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
               drkj(3,2) =  
     $          di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3
      
               drkj(1,3) =  
     $          dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
               drkj(2,3) =  
     $          dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
               drkj(3,3) =  
     $          dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3
               drkj(1,4) =  
     $          diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
               drkj(2,4) =  
     $          diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
               drkj(3,4) =  
     $          diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3
c
               bkj = gorbrep(j)
               bkjx = gorbrep(k)*gorbrep(j)
c
c   vdw radius is the vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp2(k)! + dincr_lprep(ilp)
               rdpj = vdwdisp1(j)
               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
c
               ddkj(1,1) = -(rdpk+rdpj)*bdmp*drkj(1,1)/rkj**2 
               ddkj(2,1) = -(rdpk+rdpj)*bdmp*drkj(2,1)/rkj**2 
               ddkj(3,1) = -(rdpk+rdpj)*bdmp*drkj(3,1)/rkj**2 

               ddkj(1,2) = -(rdpk+rdpj)*bdmp*drkj(1,2)/rkj**2 
               ddkj(2,2) = -(rdpk+rdpj)*bdmp*drkj(2,2)/rkj**2 
               ddkj(3,2) = -(rdpk+rdpj)*bdmp*drkj(3,2)/rkj**2 

               ddkj(1,3) = -(rdpk+rdpj)*bdmp*drkj(1,3)/rkj**2 
               ddkj(2,3) = -(rdpk+rdpj)*bdmp*drkj(2,3)/rkj**2 
               ddkj(3,3) = -(rdpk+rdpj)*bdmp*drkj(3,3)/rkj**2 

               ddkj(1,4) = -(rdpk+rdpj)*bdmp*drkj(1,4)/rkj**2 
               ddkj(2,4) = -(rdpk+rdpj)*bdmp*drkj(2,4)/rkj**2 
               ddkj(3,4) = -(rdpk+rdpj)*bdmp*drkj(3,4)/rkj**2 
c
               dmkj(1,1) = drkj(1,1)/(2*sqrt(rdpk*rdpj))
               dmkj(2,1) = drkj(2,1)/(2*sqrt(rdpk*rdpj))
               dmkj(3,1) = drkj(3,1)/(2*sqrt(rdpk*rdpj))

               dmkj(1,2) = drkj(1,2)/(2*sqrt(rdpk*rdpj))
               dmkj(2,2) = drkj(2,2)/(2*sqrt(rdpk*rdpj))
               dmkj(3,2) = drkj(3,2)/(2*sqrt(rdpk*rdpj))

               dmkj(1,3) = drkj(1,3)/(2*sqrt(rdpk*rdpj))
               dmkj(2,3) = drkj(2,3)/(2*sqrt(rdpk*rdpj))
               dmkj(3,3) = drkj(3,3)/(2*sqrt(rdpk*rdpj))

               dmkj(1,4) = drkj(1,4)/(2*sqrt(rdpk*rdpj))
               dmkj(2,4) = drkj(2,4)/(2*sqrt(rdpk*rdpj))
               dmkj(3,4) = drkj(3,4)/(2*sqrt(rdpk*rdpj))
c
               rskj = (rdpk + rdpj)*bdmp
               expo6 = -admp6*dkj
               mkj = rkj/(2*sqrt(rdpk*rdpj))
               mkj6 = mkj*mkj*mkj*mkj*mkj*mkj 
               mkj8 = mkj6*mkj*mkj
               mkj10 = mkj8*mkj*mkj
               fac1 = -0.5*clp*colpa*facdispij*discof
               if (rkj > rskj) then
                 expo6=0.0d0
                 dmp6 = bkj*exp(expo6)
                 dmp8 = bkj*exp(-admp8*dkj)
                 dmp10 = bkj*exp(-admp10*dkj)
               else
                 dmp6 = bkj*exp(expo6)
                 dmp8 = bkj*exp(-admp8*dkj)
                 dmp10 = bkj*exp(-admp10*dkj)
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,1)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,1)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,1)/mkj6
c
                 jloc = loc(j)
                 dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
                 dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
                 dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,2)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,2)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,2)/mkj6
c
                 kloc = loc(k)
                 dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
                 dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2)
                 dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,3)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,3)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,3)/mkj6
c
                 ixlploc = loc(ixlp(ilp))
                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $             taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $             taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $             taper*forcetemp(3)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,4)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,4)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,4)/mkj6
c
                 izlploc = loc(izlp(ilp))
                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $             taper*forcetemp(1) 
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $             taper*forcetemp(2) 
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $             taper*forcetemp(3) 
               end if
               disp6 = fac1*c6disp*dmp6/mkj6
               disp8 = fac1*c8disp*dmp8/mkj8
               disp10 = fac1*c10disp*dmp10/mkj10
c
               forcetemp(1) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(1,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,1)/mkj8
     $          - dmp8*8*dmkj(1,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,1)/mkj10
     $            - dmp10*10*dmkj(1,1)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(2,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,1)/mkj8
     $          - dmp8*8*dmkj(2,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,1)/mkj10
     $            - dmp10*10*dmkj(2,1)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(3,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,1)/mkj8
     $          - dmp8*8*dmkj(3,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,1)/mkj10
     $            - dmp10*10*dmkj(3,1)/(mkj**11)))
c
               jloc = loc(j)
               dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1)
               dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2)
               dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3)
c
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(1,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,2)/mkj8
     $          - dmp8*8*dmkj(1,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,2)/mkj10
     $            - dmp10*10*dmkj(1,2)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(2,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,2)/mkj8
     $          - dmp8*8*dmkj(2,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,2)/mkj10
     $            - dmp10*10*dmkj(2,2)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(3,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,2)/mkj8
     $          - dmp8*8*dmkj(3,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,2)/mkj10
     $            - dmp10*10*dmkj(3,2)/(mkj**11)))
c
               kloc = loc(k)
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(1,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,3)/mkj8
     $          - dmp8*8*dmkj(1,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,3)/mkj10
     $            - dmp10*10*dmkj(1,3)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(2,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,3)/mkj8
     $          - dmp8*8*dmkj(2,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,3)/mkj10
     $            - dmp10*10*dmkj(2,3)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(3,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,3)/mkj8
     $          - dmp8*8*dmkj(3,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,3)/mkj10
     $            - dmp10*10*dmkj(3,3)/(mkj**11)))
c
               ixlploc = loc(ixlp(ilp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           taper*forcetemp(3) 
c
               forcetemp(1) =  fac1*(
     $          c6disp*(- dmp6*6*dmkj(1,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,4)/mkj8
     $          - dmp8*8*dmkj(1,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,4)/mkj10
     $            - dmp10*10*dmkj(1,4)/(mkj**11)))
               forcetemp(2) =  fac1*(
     $         c6disp* (- dmp6*6*dmkj(2,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,4)/mkj8
     $          - dmp8*8*dmkj(2,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,4)/mkj10
     $            - dmp10*10*dmkj(2,4)/(mkj**11)))
               forcetemp(3) =  fac1*(
     $          +c6disp*(- dmp6*6*dmkj(3,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,4)/mkj8
     $          - dmp8*8*dmkj(3,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,4)/mkj10
     $            - dmp10*10*dmkj(3,4)/(mkj**11)))
c
               izlploc = loc(izlp(ilp))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
c              bkj =  frep(atomic(k))*frep(atomic(j))
               bkjx = gorbrep(k)*gorbrep(j)
c               rdpk = rlp
               rdpj = vdwdisp1(j)
c               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
c
c              exchange disp
c
               ck = rpole(1,k)
               cj = rpole(1,j)
               gdk = 1 - ck/valemtp(atomic(k)) 
               gdj = 1 - cj/valemtp(atomic(j)) 
               tkj = bkjx*gdk*gdj
c
               xdisp2 = facdispij*0.5*clp*tkj*cxdla*exp(-axdla*mkj)
cc
               forcetemp(1) = -xdisp2*axdla*dmkj(1,1) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,1) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,1) 
               jloc = loc(j)
               dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
               dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
               dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
               forcetemp(1) = -xdisp2*axdla*dmkj(1,2) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,2) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,2) 
c
               kloc = loc(k)
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               forcetemp(1) =  -xdisp2*axdla*dmkj(1,3) 
               forcetemp(2) =  -xdisp2*axdla*dmkj(2,3) 
               forcetemp(3) =  -xdisp2*axdla*dmkj(3,3) 

               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           taper*forcetemp(3) 
c
               forcetemp(1) = -xdisp2*axdla*dmkj(1,4) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,4) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,4) 

               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $           taper*forcetemp(3) 
               e = (disp6 + disp8 + disp10) + xdisp2
c
c           derivatives of the switching function
c 
               dexdisp(1,jloc) = dexdisp(1,jloc) + e*dtaper*drkj(1,1) 
               dexdisp(2,jloc) = dexdisp(2,jloc) + e*dtaper*drkj(2,1) 
               dexdisp(3,jloc) = dexdisp(3,jloc) + e*dtaper*drkj(3,1) 
               dexdisp(1,kloc) = dexdisp(1,kloc) + e*dtaper*drkj(1,2) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + e*dtaper*drkj(2,2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + e*dtaper*drkj(3,2) 
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $            e*dtaper*drkj(1,3) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $            e*dtaper*drkj(2,3) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $            e*dtaper*drkj(3,3) 
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $            e*dtaper*drkj(1,4) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $            e*dtaper*drkj(2,4) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $            e*dtaper*drkj(3,4) 

c               
               e = e*taper
               exdisp = exdisp + e
         end do 
      end do  
c
c      lone pair / lone pair dispersion term 
c
      do i = 1, nlploc
         ilp = lpglob(i)
         call rotlp1(ilp,di,dix,diz)
         k = lpatom(ilp)    
         clp1 = lpcharge(ilp)
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, nlp
            jlp = j
            if (jlp.le.ilp) cycle
            call rotlp1(jlp,di2,dix2,diz2)
            drij = 0d0
            ddkl = 0d0
            dmkl = 0d0
            l = lpatom(jlp)
            clp2 = lpcharge(jlp)
            if (molcule(k) .eq. molcule(l)) cycle
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
c
c     use energy switching if near the cutoff distance
c
               taper = 1.0d0
               dtaper = 0.0d0
               if (rij2 .gt. cut2) then
                  rij3 = rij2 * rij
                  rij4 = rij2 * rij2
                  rij5 = rij2 * rij3
                  taper = c5*rij5 + c4*rij4 + c3*rij3
     &                       + c2*rij2 + c1*rij + c0
                  dtaper = 5.0d0*c5*rij4 + 4.0d0*c4*rij3
     &                        + 3.0d0*c3*rij2 + 2.0d0*c2*rij + c1
               end if

               bkl = 1  
c
c
c    vdw radius is vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp3(k)! + dincr_lprep(ilp) 
               rdpl = vdwdisp3(l)! + dincr_lprep(jlp) 
c
               dkl = ((rdpk+rdpl)*bdmp/rij)-1
               rskl = (rdpk + rdpl)*bdmp
               mkl = rij/(2*sqrt(rdpk*rdpl))
               mkl6 = mkl*mkl*mkl*mkl*mkl*mkl 
               mkl8 = mkl6*mkl*mkl
               mkl10 = mkl8*mkl*mkl
c
               dlp1 = -xij/rij
               dlp2 = -yij/rij
               dlp3 = -zij/rij

               drij(1,1) =  
     $          di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
               drij(2,1) =  
     $          di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
               drij(3,1) =  
     $          di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3

               drij(1,2) =  
     $          dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
               drij(2,2) =  
     $          dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
               drij(3,2) =  
     $          dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3

               drij(1,3) =  
     $          diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
               drij(2,3) =  
     $          diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
               drij(3,3) =  
     $          diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3

               drij(1,4) =  
     $          -di2(1,1)*dlp1 - di2(1,2)*dlp2 - di2(1,3)*dlp3
               drij(2,4) =  
     $          -di2(2,1)*dlp1 - di2(2,2)*dlp2 - di2(2,3)*dlp3
               drij(3,4) =  
     $          -di2(3,1)*dlp1 - di2(3,2)*dlp2 - di2(3,3)*dlp3

               drij(1,5) =  
     $          -dix2(1,1)*dlp1 - dix2(1,2)*dlp2 - dix2(1,3)*dlp3
               drij(2,5) =  
     $          -dix2(2,1)*dlp1 - dix2(2,2)*dlp2 - dix2(2,3)*dlp3
               drij(3,5) =  
     $          -dix2(3,1)*dlp1 - dix2(3,2)*dlp2 - dix2(3,3)*dlp3

               drij(1,6) =  
     $          -diz2(1,1)*dlp1 - diz2(1,2)*dlp2 - diz2(1,3)*dlp3
               drij(2,6) =  
     $          -diz2(2,1)*dlp1 - diz2(2,2)*dlp2 - diz2(2,3)*dlp3
               drij(3,6) =  
     $          -diz2(3,1)*dlp1 - diz2(3,2)*dlp2 - diz2(3,3)*dlp3
      
               ddkl(1,1) = -(rdpk+rdpl)*bdmp*drij(1,1)/rij**2 
               ddkl(2,1) = -(rdpk+rdpl)*bdmp*drij(2,1)/rij**2 
               ddkl(3,1) = -(rdpk+rdpl)*bdmp*drij(3,1)/rij**2 
               ddkl(1,2) = -(rdpk+rdpl)*bdmp*drij(1,2)/rij**2 
               ddkl(2,2) = -(rdpk+rdpl)*bdmp*drij(2,2)/rij**2 
               ddkl(3,2) = -(rdpk+rdpl)*bdmp*drij(3,2)/rij**2 
               ddkl(1,3) = -(rdpk+rdpl)*bdmp*drij(1,3)/rij**2 
               ddkl(2,3) = -(rdpk+rdpl)*bdmp*drij(2,3)/rij**2 
               ddkl(3,3) = -(rdpk+rdpl)*bdmp*drij(3,3)/rij**2 
               ddkl(1,4) = -(rdpk+rdpl)*bdmp*drij(1,4)/rij**2 
               ddkl(2,4) = -(rdpk+rdpl)*bdmp*drij(2,4)/rij**2 
               ddkl(3,4) = -(rdpk+rdpl)*bdmp*drij(3,4)/rij**2 
               ddkl(1,5) = -(rdpk+rdpl)*bdmp*drij(1,5)/rij**2 
               ddkl(2,5) = -(rdpk+rdpl)*bdmp*drij(2,5)/rij**2 
               ddkl(3,5) = -(rdpk+rdpl)*bdmp*drij(3,5)/rij**2 
               ddkl(1,6) = -(rdpk+rdpl)*bdmp*drij(1,6)/rij**2 
               ddkl(2,6) = -(rdpk+rdpl)*bdmp*drij(2,6)/rij**2 
               ddkl(3,6) = -(rdpk+rdpl)*bdmp*drij(3,6)/rij**2 
      
               dmkl(1,1) = drij(1,1)/(2*sqrt(rdpk*rdpl))
               dmkl(2,1) = drij(2,1)/(2*sqrt(rdpk*rdpl))
               dmkl(3,1) = drij(3,1)/(2*sqrt(rdpk*rdpl))
               dmkl(1,2) = drij(1,2)/(2*sqrt(rdpk*rdpl))
               dmkl(2,2) = drij(2,2)/(2*sqrt(rdpk*rdpl))
               dmkl(3,2) = drij(3,2)/(2*sqrt(rdpk*rdpl))
               dmkl(1,3) = drij(1,3)/(2*sqrt(rdpk*rdpl))
               dmkl(2,3) = drij(2,3)/(2*sqrt(rdpk*rdpl))
               dmkl(3,3) = drij(3,3)/(2*sqrt(rdpk*rdpl))
               dmkl(1,4) = drij(1,4)/(2*sqrt(rdpk*rdpl))
               dmkl(2,4) = drij(2,4)/(2*sqrt(rdpk*rdpl))
               dmkl(3,4) = drij(3,4)/(2*sqrt(rdpk*rdpl))
               dmkl(1,5) = drij(1,5)/(2*sqrt(rdpk*rdpl))
               dmkl(2,5) = drij(2,5)/(2*sqrt(rdpk*rdpl))
               dmkl(3,5) = drij(3,5)/(2*sqrt(rdpk*rdpl))
               dmkl(1,6) = drij(1,6)/(2*sqrt(rdpk*rdpl))
               dmkl(2,6) = drij(2,6)/(2*sqrt(rdpk*rdpl))
               dmkl(3,6) = drij(3,6)/(2*sqrt(rdpk*rdpl))
c

               expo6 = -admp6*dkl
               fac1 = -0.5*clp1*0.5*clp2*facdispij*colp*discof
               if (rij > rskl) then
                 expo6 = 0.0d0 
                 dmp6 = bkl*exp(expo6)
                 dmp8 = bkl*exp(-admp8*dkl)
                 dmp10 = bkl*exp(-admp10*dkl)
               else
                 dmp6 = bkl*exp(expo6)
                 dmp8 = bkl*exp(-admp8*dkl)
                 dmp10 = bkl*exp(-admp10*dkl)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,1)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,1)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,1)/mkl6
                 kloc = loc(k)
                 dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
                 dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
                 dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,2)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,2)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,2)/mkl6
                 ixlploc = loc(ixlp(ilp))
                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $            taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $            taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $            taper*forcetemp(3)
c
                 izlploc = loc(izlp(ilp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,3)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,3)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,3)/mkl6

                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $            taper*forcetemp(1)
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $            taper*forcetemp(2)
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $            taper*forcetemp(3)
c
                 lloc = loc(l)
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,4)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,4)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,4)/mkl6

                 dexdisp(1,lloc) = dexdisp(1,lloc) + taper*forcetemp(1) 
                 dexdisp(2,lloc) = dexdisp(2,lloc) + taper*forcetemp(2) 
                 dexdisp(3,lloc) = dexdisp(3,lloc) + taper*forcetemp(3) 
c
                 ixlploc = loc(ixlp(jlp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,5)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,5)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,5)/mkl6

                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $              taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $              taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $              taper*forcetemp(3)
c
                 izlploc = loc(izlp(jlp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,6)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,6)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,6)/mkl6

                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $             taper*forcetemp(1)
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $             taper*forcetemp(2)
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $             taper*forcetemp(3)
               end if
c
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(1,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,1)/mkl8
     $          - dmp8*8*dmkl(1,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,1)/mkl10
     $            - dmp10*10*dmkl(1,1)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,1)/mkl8
     $          - dmp8*8*dmkl(2,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,1)/mkl10
     $            - dmp10*10*dmkl(2,1)/(mkl**11)))
               forcetemp(3) = fac1*(
     $         c6disp*(- dmp6*6*dmkl(3,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,1)/mkl8
     $          - dmp8*8*dmkl(3,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,1)/mkl10
     $            - dmp10*10*dmkl(3,1)/(mkl**11)))
               kloc = loc(k)

               dexdisp(1,kloc) = dexdisp(1,kloc) +taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) +taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) +taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(ilp))
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(1,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,2)/mkl8
     $          - dmp8*8*dmkl(1,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,2)/mkl10
     $            - dmp10*10*dmkl(1,2)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(2,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,2)/mkl8
     $          - dmp8*8*dmkl(2,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,2)/mkl10
     $            - dmp10*10*dmkl(2,2)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(3,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,2)/mkl8
     $          - dmp8*8*dmkl(3,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,2)/mkl10
     $            - dmp10*10*dmkl(3,2)/(mkl**11)))

               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(ilp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,3)/mkl8
     $          - dmp8*8*dmkl(1,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,3)/mkl10
     $            - dmp10*10*dmkl(1,3)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,3)/mkl8
     $          - dmp8*8*dmkl(2,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,3)/mkl10
     $            - dmp10*10*dmkl(2,3)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,3)/mkl8
     $          - dmp8*8*dmkl(3,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,3)/mkl10
     $            - dmp10*10*dmkl(3,3)/(mkl**11)))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
               lloc = loc(l)
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,4)/mkl8
     $          - dmp8*8*dmkl(1,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,4)/mkl10
     $            - dmp10*10*dmkl(1,4)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,4)/mkl8
     $          - dmp8*8*dmkl(2,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,4)/mkl10
     $            - dmp10*10*dmkl(2,4)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,4)/mkl8
     $          - dmp8*8*dmkl(3,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,4)/mkl10
     $            - dmp10*10*dmkl(3,4)/(mkl**11)))
               dexdisp(1,lloc) = dexdisp(1,lloc) +taper*forcetemp(1) 
               dexdisp(2,lloc) = dexdisp(2,lloc) +taper*forcetemp(2) 
               dexdisp(3,lloc) = dexdisp(3,lloc) +taper*forcetemp(3) 

               ixlploc = loc(ixlp(jlp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,5)/mkl8
     $          - dmp8*8*dmkl(1,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,5)/mkl10
     $            - dmp10*10*dmkl(1,5)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,5)/mkl8
     $          - dmp8*8*dmkl(2,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,5)/mkl10
     $            - dmp10*10*dmkl(2,5)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,5)/mkl8
     $          - dmp8*8*dmkl(3,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,5)/mkl10
     $            - dmp10*10*dmkl(3,5)/(mkl**11)))
               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(jlp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,6)/mkl8
     $          - dmp8*8*dmkl(1,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,6)/mkl10
     $            - dmp10*10*dmkl(1,6)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,6)/mkl8
     $          - dmp8*8*dmkl(2,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,6)/mkl10
     $            - dmp10*10*dmkl(2,6)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,6)/mkl8
     $          - dmp8*8*dmkl(3,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,6)/mkl10
     $            - dmp10*10*dmkl(3,6)/(mkl**11)))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)

               disp6 = fac1*c6disp*dmp6/mkl6
               disp8 = fac1*c8disp*dmp8/mkl8
               disp10 = fac1*c10disp*dmp10/mkl10
c
c          Exchange dispersion                  
c
               bkl = gorbrep(k)*gorbrep(l)
               ck = rpole(1,k)
               cl = rpole(1,l)
               gdk = 1 - ck/valemtp(atomic(k))
               gdl = 1 - cl/valemtp(atomic(l))
               tkl = bkl*gdk*gdl
               xdisp3 = facdispij*0.5*clp1*0.5*clp2*tkl*cxdlp*
     $            exp(-axdlp*mkl) 
c
               kloc = loc(k)
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,1) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,1) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,1) 
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(ilp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,2) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,2) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,2) 
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           forcetemp(3) 
c
               izlploc = loc(izlp(ilp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,3) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,3) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,3) 

               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
               lloc = loc(l)
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,4) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,4) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,4) 
               dexdisp(1,lloc) = dexdisp(1,lloc) + taper*forcetemp(1) 
               dexdisp(2,lloc) = dexdisp(2,lloc) + taper*forcetemp(2) 
               dexdisp(3,lloc) = dexdisp(3,lloc) + taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(jlp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,5) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,5) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,5) 
               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(jlp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,6) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,6) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,6) 

               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $           taper*forcetemp(3) 
              
               e = disp6 + disp8 + disp10 + xdisp3
c
c           derivatives of the switching function
c 
               dexdisp(1,kloc) = dexdisp(1,kloc) + e*dtaper*drij(1,1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + e*dtaper*drij(2,1) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + e*dtaper*drij(3,1) 
               ixlploc = loc(ixlp(ilp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $          e*dtaper*drij(1,2) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $          e*dtaper*drij(2,2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $          e*dtaper*drij(3,2) 
               izlploc = loc(izlp(ilp))
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $          e*dtaper*drij(1,3) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $          e*dtaper*drij(2,3) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $          e*dtaper*drij(3,3) 

               dexdisp(1,lloc) = dexdisp(1,lloc) +e*dtaper*drij(1,4) 
               dexdisp(2,lloc) = dexdisp(2,lloc) +e*dtaper*drij(2,4) 
               dexdisp(3,lloc) = dexdisp(3,lloc) +e*dtaper*drij(3,4) 
               ixlploc = loc(ixlp(jlp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $          e*dtaper*drij(1,5) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $          e*dtaper*drij(2,5) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $          e*dtaper*drij(3,5) 
               izlploc = loc(izlp(jlp))
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $          e*dtaper*drij(1,6) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $          e*dtaper*drij(2,6) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $          e*dtaper*drij(3,6) 

               e = e*taper
               exdisp = exdisp + e
         end do
      end do
      einter = einter + exdisp
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy and forces with other unit cells
c
c
c     atomic dispersion computation
c
      do i = 1, nloc
         iglob = glob(i)
         iloc = loc(iglob)
         xi = x(iglob) 
         yi = y(iglob) 
         zi = z(iglob) 
         do j = 1, n
            jglob = j
c            if (jglob.le.iglob) cycle
            jloc = loc(jglob)
            do k = 2, ncell
c
c     set the distance to translate along each cell axis
c
              xmove = icell(1,k) * xbox
              ymove = icell(2,k) * ybox
              zmove = icell(3,k) * zbox
c              if (molcule(iglob) .eq. molcule(jglob)) cycle
              ddij = 0d0
              dmij = 0d0
              xj = x(jglob) 
              yj = y(jglob) 
              zj = z(jglob) 
              xij = xj - xi 
              yij = yj - yi 
              zij = zj - zi 
              xij = xij + xmove
              yij = yij + ymove
              zij = zij + zmove
c              if (use_bounds)  call image (xij,yij,zij)
              rij2 = xij*xij + yij*yij + zij*zij
c
              if (rij2.gt.(dispcut*dispcut)) cycle
              rij = sqrt(rij2) 
c
c     use energy switching if near the cutoff distance
c
              taper = 1.0d0
              dtaper = 0.0d0
              if (rij2 .gt. cut2) then
                 rij3 = rij2 * rij
                 rij4 = rij2 * rij2
                 rij5 = rij2 * rij3
                 taper = c5*rij5 + c4*rij4 + c3*rij3
     &                      + c2*rij2 + c1*rij + c0
                 dtaper = 5.0d0*c5*rij4 + 4.0d0*c4*rij3
     &                       + 3.0d0*c3*rij2 + 2.0d0*c2*rij + c1
              end if
c
              bij =  gorbrep(iglob)*gorbrep(jglob)
              rdpi = vdwdisp1(iglob)
              rdpj = vdwdisp1(jglob)
              dij = ((rdpi+rdpj)*bdmp/rij)-1
c
              ddij(1,1) = (rdpi+rdpj)*bdmp*xij/rij**3 
              ddij(2,1) = (rdpi+rdpj)*bdmp*yij/rij**3 
              ddij(3,1) = (rdpi+rdpj)*bdmp*zij/rij**3 
              ddij(1,2) = -(rdpi+rdpj)*bdmp*xij/rij**3 
              ddij(2,2) = -(rdpi+rdpj)*bdmp*yij/rij**3 
              ddij(3,2) = -(rdpi+rdpj)*bdmp*zij/rij**3 
c
              dmij(1,1) = -xij/(2*sqrt(rdpi*rdpj)*rij) 
              dmij(2,1) = -yij/(2*sqrt(rdpi*rdpj)*rij) 
              dmij(3,1) = -zij/(2*sqrt(rdpi*rdpj)*rij) 
              dmij(1,2) = xij/(2*sqrt(rdpi*rdpj)*rij) 
              dmij(2,2) = yij/(2*sqrt(rdpi*rdpj)*rij) 
              dmij(3,2) = zij/(2*sqrt(rdpi*rdpj)*rij) 
c
              rsij = (rdpi + rdpj)*bdmp
              expo6 = -admp6*dij
              mij = rij/(2*sqrt(rdpi*rdpj))
              mij6 = mij*mij*mij*mij*mij*mij 
              mij8 = mij6*mij*mij
              mij10 = mij8*mij*mij

              expo6 = admp6*dij 
c             si dump pas assez efficace , > 1
              if (rij > rsij) then
                dmp6 = bij
                dmp8 = bij
                dmp10 = bij
                disp6 = c6disp*dmp6/mij6
                disp8 = c8disp*dmp8/mij8
                disp10 = c10disp*dmp10/mij10
              else
                dmp6 = bij*exp(-expo6)
                dmp8 = bij*exp(-admp8*dij)
                dmp10 = bij*exp(-admp10*dij)
                disp6 = c6disp*dmp6/mij6
                disp8 = c8disp*dmp8/mij8
                disp10 = c10disp*dmp10/mij10
c
                forcetemp(1) =facdispij*scdp*(
     $         - c6disp*admp6*dmp6*ddij(1,1)/mij6
     $         - c8disp*admp8*dmp8*ddij(1,1)/mij8
     $         - c10disp*admp10*dmp10*ddij(1,1)/mij10)*discof
                forcetemp(2) =facdispij*scdp*(
     $         - c6disp*admp6*dmp6*ddij(2,1)/mij6
     $         - c8disp*admp8*dmp8*ddij(2,1)/mij8
     $         - c10disp*admp10*dmp10*ddij(2,1)/mij10)*discof
                forcetemp(3) =facdispij*scdp*(
     $         - c6disp*admp6*dmp6*ddij(3,1)/mij6
     $         - c8disp*admp8*dmp8*ddij(3,1)/mij8
     $         - c10disp*admp10*dmp10*ddij(3,1)/mij10)*discof

                dexdisp(1,iloc) = dexdisp(1,iloc)+taper*forcetemp(1)
                dexdisp(2,iloc) = dexdisp(2,iloc)+taper*forcetemp(2)
                dexdisp(3,iloc) = dexdisp(3,iloc)+taper*forcetemp(3)
c
                forcetemp(1) =facdispij*scdp*(
     $         - c6disp*admp6*dmp6*ddij(1,2)/mij6
     $         - c8disp*admp8*dmp8*ddij(1,2)/mij8
     $         - c10disp*admp10*dmp10*ddij(1,2)/mij10)*discof
                forcetemp(2) =facdispij*scdp*(
     $         - c6disp*admp6*dmp6*ddij(2,2)/mij6
     $         - c8disp*admp8*dmp8*ddij(2,2)/mij8
     $         - c10disp*admp10*dmp10*ddij(2,2)/mij10)*discof
                forcetemp(3) =facdispij*scdp*(
     $         - c6disp*admp6*dmp6*ddij(3,2)/mij6
     $         - c8disp*admp8*dmp8*ddij(3,2)/mij8
     $         - c10disp*admp10*dmp10*ddij(3,2)/mij10)*discof

                dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1)
                dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2)
                dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3)
              end if
c
              forcetemp(1) = facdispij*scdp*(
     $         - c6disp*dmp6*6*dmij(1,1)/(mij**7)
     $           - c8disp*dmp8*8*dmij(1,1)/(mij**9)
     $           - c10disp*dmp10*10*dmij(1,1)/(mij**11))*discof
              forcetemp(2) = facdispij*scdp*(
     $         - c6disp*dmp6*6*dmij(2,1)/(mij**7)
     $           - c8disp*dmp8*8*dmij(2,1)/(mij**9)
     $           - c10disp*dmp10*10*dmij(2,1)/(mij**11))*discof
              forcetemp(3) = facdispij*scdp*(
     $         - c6disp*dmp6*6*dmij(3,1)/(mij**7)
     $           - c8disp*dmp8*8*dmij(3,1)/(mij**9)
     $           - c10disp*dmp10*10*dmij(3,1)/(mij**11))*discof
c
              dexdisp(1,iloc) = dexdisp(1,iloc) + taper*forcetemp(1)
              dexdisp(2,iloc) = dexdisp(2,iloc) + taper*forcetemp(2)
              dexdisp(3,iloc) = dexdisp(3,iloc) + taper*forcetemp(3)
c
              forcetemp(1) = facdispij*scdp*(
     $         - c6disp*dmp6*6*dmij(1,2)/(mij**7)
     $           - c8disp*dmp8*8*dmij(1,2)/(mij**9)
     $           - c10disp*dmp10*10*dmij(1,2)/(mij**11))*discof
              forcetemp(2) = facdispij*scdp*(
     $         - c6disp*dmp6*6*dmij(2,2)/(mij**7)
     $           - c8disp*dmp8*8*dmij(2,2)/(mij**9)
     $           - c10disp*dmp10*10*dmij(2,2)/(mij**11))*discof
              forcetemp(3) = facdispij*scdp*(
     $         - c6disp*dmp6*6*dmij(3,2)/(mij**7)
     $           - c8disp*dmp8*8*dmij(3,2)/(mij**9)
     $           - c10disp*dmp10*10*dmij(3,2)/(mij**11))*discof
c
              dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
              dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
              dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
c             exchange disp              
c
              ci = rpole(1,iglob) 
              cj = rpole(1,jglob)
              gdi = 1 - ci/valemtp(atomic(iglob)) 
              gdj = 1 - cj/valemtp(atomic(jglob)) 
              tij = bij*gdi*gdj
              xdisp1 = tij*cxd*exp(-axd*mij)
c
              forcetemp(1) = - facdispij*xdisp1*axd*dmij(1,1) 
              forcetemp(2) = - facdispij*xdisp1*axd*dmij(2,1) 
              forcetemp(3) = - facdispij*xdisp1*axd*dmij(3,1) 

              dexdisp(1,iloc) =dexdisp(1,iloc)+forcetemp(1)
              dexdisp(2,iloc) =dexdisp(2,iloc)+forcetemp(2)
              dexdisp(3,iloc) =dexdisp(3,iloc)+forcetemp(3)
c
              forcetemp(1) = - facdispij*xdisp1*axd*dmij(1,2) 
              forcetemp(2) = - facdispij*xdisp1*axd*dmij(2,2) 
              forcetemp(3) = - facdispij*xdisp1*axd*dmij(3,2) 
c
              dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
              dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
              dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 

              e = facdispij*scdp*(disp6 + disp8 + disp10)*discof
     $         + xdisp1*facdispij
c
c             derivatives of the switching function
c 
              dedx = -e*dtaper*xij/rij
              dedy = -e*dtaper*yij/rij
              dedz = -e*dtaper*zij/rij
              dexdisp(1,iloc) = dexdisp(1,iloc) + dedx 
              dexdisp(2,iloc) = dexdisp(2,iloc) + dedy 
              dexdisp(3,iloc) = dexdisp(3,iloc) + dedz 
              dexdisp(1,jloc) = dexdisp(1,jloc) - dedx 
              dexdisp(2,jloc) = dexdisp(2,jloc) - dedy 
              dexdisp(3,jloc) = dexdisp(3,jloc) - dedz 
c
              e = e*taper
              exdisp = exdisp + e
              nexdisp = nexdisp + 1
            end do
         end do
      end do  
c
c    lp-atom dispersion
c
      do i = 1, nlploc
         ilp = lpglob(i)
         clp = lpcharge(ilp)!2.0d0
         k = lpatom(ilp)
         xlp = rlonepair(1,ilp)
         ylp = rlonepair(2,ilp)
         zlp = rlonepair(3,ilp)
         call rotlp1(ilp,di,dix,diz)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, n
            jloc = loc(j)
            do l = 2, ncell
c
c     set the distance to translate along each cell axis
c
               xmove = icell(1,l) * xbox
               ymove = icell(2,l) * ybox
               zmove = icell(3,l) * zbox
c               if (molcule(k) .eq. molcule(j)) cycle
               drkj = 0d0
               ddkj = 0d0
               dmkj = 0d0
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
c               if (use_bounds)  call image (xkj,ykj,zkj)
               rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
c
               if (rkj2.gt.(dispcut*dispcut)) cycle
               rkj = sqrt(rkj2)
c
c     use energy switching if near the cutoff distance
c
               taper = 1.0d0
               dtaper = 0.0d0
               if (rkj2 .gt. cut2) then
                  rkj3 = rkj2 * rkj
                  rkj4 = rkj2 * rkj2
                  rkj5 = rkj2 * rkj3
                  taper = c5*rkj5 + c4*rkj4 + c3*rkj3
     &                       + c2*rkj2 + c1*rkj + c0
                  dtaper = 5.0d0*c5*rkj4 + 4.0d0*c4*rkj3
     &                        + 3.0d0*c3*rkj2 + 2.0d0*c2*rkj + c1
               end if
c
c
               drkj(1,1) = xkj/rkj
               drkj(2,1) = ykj/rkj
               drkj(3,1) = zkj/rkj
c
               dlp1 = -xkj/rkj
               dlp2 = -ykj/rkj
               dlp3 = -zkj/rkj
               drkj(1,2) =  
     $          di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
               drkj(2,2) =  
     $          di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
               drkj(3,2) =  
     $          di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3
      
               drkj(1,3) =  
     $          dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
               drkj(2,3) =  
     $          dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
               drkj(3,3) =  
     $          dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3
               drkj(1,4) =  
     $          diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
               drkj(2,4) =  
     $          diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
               drkj(3,4) =  
     $          diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3
c
               bkj = gorbrep(j)
               bkjx = gorbrep(k)*gorbrep(j)
c
c   vdw radius is the vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp2(k)! + dincr_lprep(ilp)
               rdpj = vdwdisp1(j)
               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
c
               ddkj(1,1) = -(rdpk+rdpj)*bdmp*drkj(1,1)/rkj**2 
               ddkj(2,1) = -(rdpk+rdpj)*bdmp*drkj(2,1)/rkj**2 
               ddkj(3,1) = -(rdpk+rdpj)*bdmp*drkj(3,1)/rkj**2 

               ddkj(1,2) = -(rdpk+rdpj)*bdmp*drkj(1,2)/rkj**2 
               ddkj(2,2) = -(rdpk+rdpj)*bdmp*drkj(2,2)/rkj**2 
               ddkj(3,2) = -(rdpk+rdpj)*bdmp*drkj(3,2)/rkj**2 

               ddkj(1,3) = -(rdpk+rdpj)*bdmp*drkj(1,3)/rkj**2 
               ddkj(2,3) = -(rdpk+rdpj)*bdmp*drkj(2,3)/rkj**2 
               ddkj(3,3) = -(rdpk+rdpj)*bdmp*drkj(3,3)/rkj**2 

               ddkj(1,4) = -(rdpk+rdpj)*bdmp*drkj(1,4)/rkj**2 
               ddkj(2,4) = -(rdpk+rdpj)*bdmp*drkj(2,4)/rkj**2 
               ddkj(3,4) = -(rdpk+rdpj)*bdmp*drkj(3,4)/rkj**2 
c
               dmkj(1,1) = drkj(1,1)/(2*sqrt(rdpk*rdpj))
               dmkj(2,1) = drkj(2,1)/(2*sqrt(rdpk*rdpj))
               dmkj(3,1) = drkj(3,1)/(2*sqrt(rdpk*rdpj))

               dmkj(1,2) = drkj(1,2)/(2*sqrt(rdpk*rdpj))
               dmkj(2,2) = drkj(2,2)/(2*sqrt(rdpk*rdpj))
               dmkj(3,2) = drkj(3,2)/(2*sqrt(rdpk*rdpj))

               dmkj(1,3) = drkj(1,3)/(2*sqrt(rdpk*rdpj))
               dmkj(2,3) = drkj(2,3)/(2*sqrt(rdpk*rdpj))
               dmkj(3,3) = drkj(3,3)/(2*sqrt(rdpk*rdpj))

               dmkj(1,4) = drkj(1,4)/(2*sqrt(rdpk*rdpj))
               dmkj(2,4) = drkj(2,4)/(2*sqrt(rdpk*rdpj))
               dmkj(3,4) = drkj(3,4)/(2*sqrt(rdpk*rdpj))
c
               rskj = (rdpk + rdpj)*bdmp
               expo6 = -admp6*dkj
               mkj = rkj/(2*sqrt(rdpk*rdpj))
               mkj6 = mkj*mkj*mkj*mkj*mkj*mkj 
               mkj8 = mkj6*mkj*mkj
               mkj10 = mkj8*mkj*mkj
               fac1 = -0.5*clp*colpa*facdispij*discof
               if (rkj > rskj) then
                 expo6=0.0d0
                 dmp6 = bkj*exp(expo6)
                 dmp8 = bkj*exp(-admp8*dkj)
                 dmp10 = bkj*exp(-admp10*dkj)
               else
                 dmp6 = bkj*exp(expo6)
                 dmp8 = bkj*exp(-admp8*dkj)
                 dmp10 = bkj*exp(-admp10*dkj)
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,1)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,1)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,1)/mkj6
c
                 jloc = loc(j)
                 dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
                 dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
                 dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,2)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,2)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,2)/mkj6
c
                 kloc = loc(k)
                 dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
                 dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2)
                 dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,3)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,3)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,3)/mkj6
c
                 ixlploc = loc(ixlp(ilp))
                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $             taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $             taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $             taper*forcetemp(3)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkj(1,4)/mkj6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkj(2,4)/mkj6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkj(3,4)/mkj6
c
                 izlploc = loc(izlp(ilp))
                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $             taper*forcetemp(1) 
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $             taper*forcetemp(2) 
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $             taper*forcetemp(3) 
               end if
               disp6 = fac1*c6disp*dmp6/mkj6
               disp8 = fac1*c8disp*dmp8/mkj8
               disp10 = fac1*c10disp*dmp10/mkj10
c
               forcetemp(1) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(1,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,1)/mkj8
     $          - dmp8*8*dmkj(1,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,1)/mkj10
     $            - dmp10*10*dmkj(1,1)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(2,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,1)/mkj8
     $          - dmp8*8*dmkj(2,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,1)/mkj10
     $            - dmp10*10*dmkj(2,1)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(-dmp6*6*dmkj(3,1)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,1)/mkj8
     $          - dmp8*8*dmkj(3,1)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,1)/mkj10
     $            - dmp10*10*dmkj(3,1)/(mkj**11)))
c
               jloc = loc(j)
               dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1)
               dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2)
               dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3)
c
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(1,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,2)/mkj8
     $          - dmp8*8*dmkj(1,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,2)/mkj10
     $            - dmp10*10*dmkj(1,2)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(2,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,2)/mkj8
     $          - dmp8*8*dmkj(2,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,2)/mkj10
     $            - dmp10*10*dmkj(2,2)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(- dmp6*6*dmkj(3,2)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,2)/mkj8
     $          - dmp8*8*dmkj(3,2)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,2)/mkj10
     $            - dmp10*10*dmkj(3,2)/(mkj**11)))
c
               kloc = loc(k)
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(1,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,3)/mkj8
     $          - dmp8*8*dmkj(1,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,3)/mkj10
     $            - dmp10*10*dmkj(1,3)/(mkj**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(2,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,3)/mkj8
     $          - dmp8*8*dmkj(2,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,3)/mkj10
     $            - dmp10*10*dmkj(2,3)/(mkj**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkj(3,3)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,3)/mkj8
     $          - dmp8*8*dmkj(3,3)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,3)/mkj10
     $            - dmp10*10*dmkj(3,3)/(mkj**11)))
c
               ixlploc = loc(ixlp(ilp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           taper*forcetemp(3) 
c
               forcetemp(1) =  fac1*(
     $          c6disp*(- dmp6*6*dmkj(1,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(1,4)/mkj8
     $          - dmp8*8*dmkj(1,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(1,4)/mkj10
     $            - dmp10*10*dmkj(1,4)/(mkj**11)))
               forcetemp(2) =  fac1*(
     $         c6disp* (- dmp6*6*dmkj(2,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(2,4)/mkj8
     $          - dmp8*8*dmkj(2,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(2,4)/mkj10
     $            - dmp10*10*dmkj(2,4)/(mkj**11)))
               forcetemp(3) =  fac1*(
     $          +c6disp*(- dmp6*6*dmkj(3,4)/(mkj**7))
     $          + c8disp*(-admp8*dmp8*ddkj(3,4)/mkj8
     $          - dmp8*8*dmkj(3,4)/(mkj**9))
     $          + c10disp*(-admp10*dmp10*ddkj(3,4)/mkj10
     $            - dmp10*10*dmkj(3,4)/(mkj**11)))
c
               izlploc = loc(izlp(ilp))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
c              bkj =  frep(atomic(k))*frep(atomic(j))
               bkjx = gorbrep(k)*gorbrep(j)
c               rdpk = rlp
               rdpj = vdwdisp1(j)
c               dkj = ((rdpk+rdpj)*bdmp/rkj)-1
c
c              exchange disp
c
               ck = rpole(1,k)
               cj = rpole(1,j)
               gdk = 1 - ck/valemtp(atomic(k)) 
               gdj = 1 - cj/valemtp(atomic(j)) 
               tkj = bkjx*gdk*gdj
c
               xdisp2 = facdispij*0.5*clp*tkj*cxdla*exp(-axdla*mkj)
cc
               forcetemp(1) = -xdisp2*axdla*dmkj(1,1) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,1) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,1) 
               jloc = loc(j)
               dexdisp(1,jloc) = dexdisp(1,jloc) + taper*forcetemp(1) 
               dexdisp(2,jloc) = dexdisp(2,jloc) + taper*forcetemp(2) 
               dexdisp(3,jloc) = dexdisp(3,jloc) + taper*forcetemp(3) 
c
               forcetemp(1) = -xdisp2*axdla*dmkj(1,2) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,2) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,2) 
c
               kloc = loc(k)
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               forcetemp(1) =  -xdisp2*axdla*dmkj(1,3) 
               forcetemp(2) =  -xdisp2*axdla*dmkj(2,3) 
               forcetemp(3) =  -xdisp2*axdla*dmkj(3,3) 

               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           taper*forcetemp(3) 
c
               forcetemp(1) = -xdisp2*axdla*dmkj(1,4) 
               forcetemp(2) = -xdisp2*axdla*dmkj(2,4) 
               forcetemp(3) = -xdisp2*axdla*dmkj(3,4) 

               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $           taper*forcetemp(3) 
               e = (disp6 + disp8 + disp10) + xdisp2
c
c           derivatives of the switching function
c 
               dexdisp(1,jloc) = dexdisp(1,jloc) + e*dtaper*drkj(1,1) 
               dexdisp(2,jloc) = dexdisp(2,jloc) + e*dtaper*drkj(2,1) 
               dexdisp(3,jloc) = dexdisp(3,jloc) + e*dtaper*drkj(3,1) 
               dexdisp(1,kloc) = dexdisp(1,kloc) + e*dtaper*drkj(1,2) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + e*dtaper*drkj(2,2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + e*dtaper*drkj(3,2) 
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $            e*dtaper*drkj(1,3) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $            e*dtaper*drkj(2,3) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $            e*dtaper*drkj(3,3) 
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $            e*dtaper*drkj(1,4) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $            e*dtaper*drkj(2,4) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $            e*dtaper*drkj(3,4) 

c               
               e = e*taper
               exdisp = exdisp + e
             end do
         end do 
      end do  
c
c      lone pair / lone pair dispersion term 
c
      do i = 1, nlploc
         ilp = lpglob(i)
         call rotlp1(ilp,di,dix,diz)
         k = lpatom(ilp)    
         clp1 = lpcharge(ilp)
         xi = rlonepair(1,ilp)
         yi = rlonepair(2,ilp)
         zi = rlonepair(3,ilp)
         xk = x(k) 
         yk = y(k) 
         zk = z(k) 
         do j = 1, nlp
            jlp = j
            call rotlp1(jlp,di2,dix2,diz2)
            drij = 0d0
            ddkl = 0d0
            dmkl = 0d0
            l = lpatom(jlp)
            clp2 = lpcharge(jlp)
            do m = 2, ncell
c
c     set the distance to translate along each cell axis
c
               xmove = icell(1,m) * xbox
               ymove = icell(2,m) * ybox
               zmove = icell(3,m) * zbox
c            if (molcule(k) .eq. molcule(l)) cycle
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
c               if (use_bounds)  call image (xij,yij,zij)
               rij2 = xij*xij + yij*yij + zij*zij
c
               if (rij2.gt.(dispcut*dispcut)) cycle
c
               rij = sqrt(rij2)
c
c     use energy switching if near the cutoff distance
c
               taper = 1.0d0
               dtaper = 0.0d0
               if (rij2 .gt. cut2) then
                  rij3 = rij2 * rij
                  rij4 = rij2 * rij2
                  rij5 = rij2 * rij3
                  taper = c5*rij5 + c4*rij4 + c3*rij3
     &                       + c2*rij2 + c1*rij + c0
                  dtaper = 5.0d0*c5*rij4 + 4.0d0*c4*rij3
     &                        + 3.0d0*c3*rij2 + 2.0d0*c2*rij + c1
               end if

               bkl = 1  
c
c
c    vdw radius is vdw of the carrier! plus increment of the lp
c
               rdpk = vdwdisp3(k)! + dincr_lprep(ilp) 
               rdpl = vdwdisp3(l)! + dincr_lprep(jlp) 
c
               dkl = ((rdpk+rdpl)*bdmp/rij)-1
               rskl = (rdpk + rdpl)*bdmp
               mkl = rij/(2*sqrt(rdpk*rdpl))
               mkl6 = mkl*mkl*mkl*mkl*mkl*mkl 
               mkl8 = mkl6*mkl*mkl
               mkl10 = mkl8*mkl*mkl
c
               dlp1 = -xij/rij
               dlp2 = -yij/rij
               dlp3 = -zij/rij

               drij(1,1) =  
     $          di(1,1)*dlp1 + di(1,2)*dlp2 + di(1,3)*dlp3
               drij(2,1) =  
     $          di(2,1)*dlp1 + di(2,2)*dlp2 + di(2,3)*dlp3
               drij(3,1) =  
     $          di(3,1)*dlp1 + di(3,2)*dlp2 + di(3,3)*dlp3

               drij(1,2) =  
     $          dix(1,1)*dlp1 + dix(1,2)*dlp2 + dix(1,3)*dlp3
               drij(2,2) =  
     $          dix(2,1)*dlp1 + dix(2,2)*dlp2 + dix(2,3)*dlp3
               drij(3,2) =  
     $          dix(3,1)*dlp1 + dix(3,2)*dlp2 + dix(3,3)*dlp3

               drij(1,3) =  
     $          diz(1,1)*dlp1 + diz(1,2)*dlp2 + diz(1,3)*dlp3
               drij(2,3) =  
     $          diz(2,1)*dlp1 + diz(2,2)*dlp2 + diz(2,3)*dlp3
               drij(3,3) =  
     $          diz(3,1)*dlp1 + diz(3,2)*dlp2 + diz(3,3)*dlp3

               drij(1,4) =  
     $          -di2(1,1)*dlp1 - di2(1,2)*dlp2 - di2(1,3)*dlp3
               drij(2,4) =  
     $          -di2(2,1)*dlp1 - di2(2,2)*dlp2 - di2(2,3)*dlp3
               drij(3,4) =  
     $          -di2(3,1)*dlp1 - di2(3,2)*dlp2 - di2(3,3)*dlp3

               drij(1,5) =  
     $          -dix2(1,1)*dlp1 - dix2(1,2)*dlp2 - dix2(1,3)*dlp3
               drij(2,5) =  
     $          -dix2(2,1)*dlp1 - dix2(2,2)*dlp2 - dix2(2,3)*dlp3
               drij(3,5) =  
     $          -dix2(3,1)*dlp1 - dix2(3,2)*dlp2 - dix2(3,3)*dlp3

               drij(1,6) =  
     $          -diz2(1,1)*dlp1 - diz2(1,2)*dlp2 - diz2(1,3)*dlp3
               drij(2,6) =  
     $          -diz2(2,1)*dlp1 - diz2(2,2)*dlp2 - diz2(2,3)*dlp3
               drij(3,6) =  
     $          -diz2(3,1)*dlp1 - diz2(3,2)*dlp2 - diz2(3,3)*dlp3
      
               ddkl(1,1) = -(rdpk+rdpl)*bdmp*drij(1,1)/rij**2 
               ddkl(2,1) = -(rdpk+rdpl)*bdmp*drij(2,1)/rij**2 
               ddkl(3,1) = -(rdpk+rdpl)*bdmp*drij(3,1)/rij**2 
               ddkl(1,2) = -(rdpk+rdpl)*bdmp*drij(1,2)/rij**2 
               ddkl(2,2) = -(rdpk+rdpl)*bdmp*drij(2,2)/rij**2 
               ddkl(3,2) = -(rdpk+rdpl)*bdmp*drij(3,2)/rij**2 
               ddkl(1,3) = -(rdpk+rdpl)*bdmp*drij(1,3)/rij**2 
               ddkl(2,3) = -(rdpk+rdpl)*bdmp*drij(2,3)/rij**2 
               ddkl(3,3) = -(rdpk+rdpl)*bdmp*drij(3,3)/rij**2 
               ddkl(1,4) = -(rdpk+rdpl)*bdmp*drij(1,4)/rij**2 
               ddkl(2,4) = -(rdpk+rdpl)*bdmp*drij(2,4)/rij**2 
               ddkl(3,4) = -(rdpk+rdpl)*bdmp*drij(3,4)/rij**2 
               ddkl(1,5) = -(rdpk+rdpl)*bdmp*drij(1,5)/rij**2 
               ddkl(2,5) = -(rdpk+rdpl)*bdmp*drij(2,5)/rij**2 
               ddkl(3,5) = -(rdpk+rdpl)*bdmp*drij(3,5)/rij**2 
               ddkl(1,6) = -(rdpk+rdpl)*bdmp*drij(1,6)/rij**2 
               ddkl(2,6) = -(rdpk+rdpl)*bdmp*drij(2,6)/rij**2 
               ddkl(3,6) = -(rdpk+rdpl)*bdmp*drij(3,6)/rij**2 
      
               dmkl(1,1) = drij(1,1)/(2*sqrt(rdpk*rdpl))
               dmkl(2,1) = drij(2,1)/(2*sqrt(rdpk*rdpl))
               dmkl(3,1) = drij(3,1)/(2*sqrt(rdpk*rdpl))
               dmkl(1,2) = drij(1,2)/(2*sqrt(rdpk*rdpl))
               dmkl(2,2) = drij(2,2)/(2*sqrt(rdpk*rdpl))
               dmkl(3,2) = drij(3,2)/(2*sqrt(rdpk*rdpl))
               dmkl(1,3) = drij(1,3)/(2*sqrt(rdpk*rdpl))
               dmkl(2,3) = drij(2,3)/(2*sqrt(rdpk*rdpl))
               dmkl(3,3) = drij(3,3)/(2*sqrt(rdpk*rdpl))
               dmkl(1,4) = drij(1,4)/(2*sqrt(rdpk*rdpl))
               dmkl(2,4) = drij(2,4)/(2*sqrt(rdpk*rdpl))
               dmkl(3,4) = drij(3,4)/(2*sqrt(rdpk*rdpl))
               dmkl(1,5) = drij(1,5)/(2*sqrt(rdpk*rdpl))
               dmkl(2,5) = drij(2,5)/(2*sqrt(rdpk*rdpl))
               dmkl(3,5) = drij(3,5)/(2*sqrt(rdpk*rdpl))
               dmkl(1,6) = drij(1,6)/(2*sqrt(rdpk*rdpl))
               dmkl(2,6) = drij(2,6)/(2*sqrt(rdpk*rdpl))
               dmkl(3,6) = drij(3,6)/(2*sqrt(rdpk*rdpl))
c

               expo6 = -admp6*dkl
               fac1 = -0.5*clp1*0.5*clp2*facdispij*colp*discof
               if (rij > rskl) then
                 expo6 = 0.0d0 
                 dmp6 = bkl*exp(expo6)
                 dmp8 = bkl*exp(-admp8*dkl)
                 dmp10 = bkl*exp(-admp10*dkl)
               else
                 dmp6 = bkl*exp(expo6)
                 dmp8 = bkl*exp(-admp8*dkl)
                 dmp10 = bkl*exp(-admp10*dkl)
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,1)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,1)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,1)/mkl6
                 kloc = loc(k)
                 dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
                 dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
                 dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,2)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,2)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,2)/mkl6
                 ixlploc = loc(ixlp(ilp))
                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $            taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $            taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $            taper*forcetemp(3)
c
                 izlploc = loc(izlp(ilp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,3)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,3)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,3)/mkl6

                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $            taper*forcetemp(1)
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $            taper*forcetemp(2)
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $            taper*forcetemp(3)
c
                 lloc = loc(l)
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,4)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,4)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,4)/mkl6

                 dexdisp(1,lloc) = dexdisp(1,lloc) + taper*forcetemp(1) 
                 dexdisp(2,lloc) = dexdisp(2,lloc) + taper*forcetemp(2) 
                 dexdisp(3,lloc) = dexdisp(3,lloc) + taper*forcetemp(3) 
c
                 ixlploc = loc(ixlp(jlp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,5)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,5)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,5)/mkl6

                 dexdisp(1,ixlploc) = dexdisp(1,ixlploc) +
     $              taper*forcetemp(1)
                 dexdisp(2,ixlploc) = dexdisp(2,ixlploc) +
     $              taper*forcetemp(2)
                 dexdisp(3,ixlploc) = dexdisp(3,ixlploc) +
     $              taper*forcetemp(3)
c
                 izlploc = loc(izlp(jlp))
                 forcetemp(1) =-c6disp*fac1*admp6*dmp6*ddkl(1,6)/mkl6
                 forcetemp(2) =-c6disp*fac1*admp6*dmp6*ddkl(2,6)/mkl6
                 forcetemp(3) =-c6disp*fac1*admp6*dmp6*ddkl(3,6)/mkl6

                 dexdisp(1,izlploc) = dexdisp(1,izlploc) +
     $             taper*forcetemp(1)
                 dexdisp(2,izlploc) = dexdisp(2,izlploc) +
     $             taper*forcetemp(2)
                 dexdisp(3,izlploc) = dexdisp(3,izlploc) +
     $             taper*forcetemp(3)
               end if
c
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(1,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,1)/mkl8
     $          - dmp8*8*dmkl(1,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,1)/mkl10
     $            - dmp10*10*dmkl(1,1)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,1)/mkl8
     $          - dmp8*8*dmkl(2,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,1)/mkl10
     $            - dmp10*10*dmkl(2,1)/(mkl**11)))
               forcetemp(3) = fac1*(
     $         c6disp*(- dmp6*6*dmkl(3,1)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,1)/mkl8
     $          - dmp8*8*dmkl(3,1)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,1)/mkl10
     $            - dmp10*10*dmkl(3,1)/(mkl**11)))
               kloc = loc(k)

               dexdisp(1,kloc) = dexdisp(1,kloc) +taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) +taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) +taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(ilp))
               forcetemp(1) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(1,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,2)/mkl8
     $          - dmp8*8*dmkl(1,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,2)/mkl10
     $            - dmp10*10*dmkl(1,2)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(2,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,2)/mkl8
     $          - dmp8*8*dmkl(2,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,2)/mkl10
     $            - dmp10*10*dmkl(2,2)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          + c6disp*(- dmp6*6*dmkl(3,2)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,2)/mkl8
     $          - dmp8*8*dmkl(3,2)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,2)/mkl10
     $            - dmp10*10*dmkl(3,2)/(mkl**11)))

               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(ilp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,3)/mkl8
     $          - dmp8*8*dmkl(1,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,3)/mkl10
     $            - dmp10*10*dmkl(1,3)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,3)/mkl8
     $          - dmp8*8*dmkl(2,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,3)/mkl10
     $            - dmp10*10*dmkl(2,3)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,3)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,3)/mkl8
     $          - dmp8*8*dmkl(3,3)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,3)/mkl10
     $            - dmp10*10*dmkl(3,3)/(mkl**11)))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
               lloc = loc(l)
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,4)/mkl8
     $          - dmp8*8*dmkl(1,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,4)/mkl10
     $            - dmp10*10*dmkl(1,4)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,4)/mkl8
     $          - dmp8*8*dmkl(2,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,4)/mkl10
     $            - dmp10*10*dmkl(2,4)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,4)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,4)/mkl8
     $          - dmp8*8*dmkl(3,4)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,4)/mkl10
     $            - dmp10*10*dmkl(3,4)/(mkl**11)))
               dexdisp(1,lloc) = dexdisp(1,lloc) +taper*forcetemp(1) 
               dexdisp(2,lloc) = dexdisp(2,lloc) +taper*forcetemp(2) 
               dexdisp(3,lloc) = dexdisp(3,lloc) +taper*forcetemp(3) 

               ixlploc = loc(ixlp(jlp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,5)/mkl8
     $          - dmp8*8*dmkl(1,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,5)/mkl10
     $            - dmp10*10*dmkl(1,5)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,5)/mkl8
     $          - dmp8*8*dmkl(2,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,5)/mkl10
     $            - dmp10*10*dmkl(2,5)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,5)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,5)/mkl8
     $          - dmp8*8*dmkl(3,5)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,5)/mkl10
     $            - dmp10*10*dmkl(3,5)/(mkl**11)))
               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(jlp))
               forcetemp(1) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(1,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(1,6)/mkl8
     $          - dmp8*8*dmkl(1,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(1,6)/mkl10
     $            - dmp10*10*dmkl(1,6)/(mkl**11)))
               forcetemp(2) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(2,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(2,6)/mkl8
     $          - dmp8*8*dmkl(2,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(2,6)/mkl10
     $            - dmp10*10*dmkl(2,6)/(mkl**11)))
               forcetemp(3) = fac1*(
     $          c6disp*(- dmp6*6*dmkl(3,6)/(mkl**7))
     $          + c8disp*(-admp8*dmp8*ddkl(3,6)/mkl8
     $          - dmp8*8*dmkl(3,6)/(mkl**9))
     $          + c10disp*(-admp10*dmp10*ddkl(3,6)/mkl10
     $            - dmp10*10*dmkl(3,6)/(mkl**11)))
               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)

               disp6 = fac1*c6disp*dmp6/mkl6
               disp8 = fac1*c8disp*dmp8/mkl8
               disp10 = fac1*c10disp*dmp10/mkl10
c
c          Exchange dispersion                  
c
               bkl = gorbrep(k)*gorbrep(l)
               ck = rpole(1,k)
               cl = rpole(1,l)
               gdk = 1 - ck/valemtp(atomic(k))
               gdl = 1 - cl/valemtp(atomic(l))
               tkl = bkl*gdk*gdl
               xdisp3 = facdispij*0.5*clp1*0.5*clp2*tkl*cxdlp*
     $            exp(-axdlp*mkl) 
c
               kloc = loc(k)
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,1) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,1) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,1) 
               dexdisp(1,kloc) = dexdisp(1,kloc) + taper*forcetemp(1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + taper*forcetemp(2) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(ilp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,2) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,2) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,2) 
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $           forcetemp(1) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $           forcetemp(2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $           forcetemp(3) 
c
               izlploc = loc(izlp(ilp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,3) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,3) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,3) 

               dexdisp(1,izlploc) =dexdisp(1,izlploc)+taper*forcetemp(1)
               dexdisp(2,izlploc) =dexdisp(2,izlploc)+taper*forcetemp(2)
               dexdisp(3,izlploc) =dexdisp(3,izlploc)+taper*forcetemp(3)
c
               lloc = loc(l)
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,4) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,4) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,4) 
               dexdisp(1,lloc) = dexdisp(1,lloc) + taper*forcetemp(1) 
               dexdisp(2,lloc) = dexdisp(2,lloc) + taper*forcetemp(2) 
               dexdisp(3,lloc) = dexdisp(3,lloc) + taper*forcetemp(3) 
c
               ixlploc = loc(ixlp(jlp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,5) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,5) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,5) 
               dexdisp(1,ixlploc) =dexdisp(1,ixlploc)+taper*forcetemp(1)
               dexdisp(2,ixlploc) =dexdisp(2,ixlploc)+taper*forcetemp(2)
               dexdisp(3,ixlploc) =dexdisp(3,ixlploc)+taper*forcetemp(3)
c
               izlploc = loc(izlp(jlp))
               forcetemp(1) = -xdisp3*axdlp*dmkl(1,6) 
               forcetemp(2) = -xdisp3*axdlp*dmkl(2,6) 
               forcetemp(3) = -xdisp3*axdlp*dmkl(3,6) 

               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $           taper*forcetemp(1) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $           taper*forcetemp(2) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $           taper*forcetemp(3) 
              
               e = disp6 + disp8 + disp10 + xdisp3
c
c           derivatives of the switching function
c 
               dexdisp(1,kloc) = dexdisp(1,kloc) + e*dtaper*drij(1,1) 
               dexdisp(2,kloc) = dexdisp(2,kloc) + e*dtaper*drij(2,1) 
               dexdisp(3,kloc) = dexdisp(3,kloc) + e*dtaper*drij(3,1) 
               ixlploc = loc(ixlp(ilp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $          e*dtaper*drij(1,2) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $          e*dtaper*drij(2,2) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $          e*dtaper*drij(3,2) 
               izlploc = loc(izlp(ilp))
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $          e*dtaper*drij(1,3) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $          e*dtaper*drij(2,3) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $          e*dtaper*drij(3,3) 

               dexdisp(1,lloc) = dexdisp(1,lloc) +e*dtaper*drij(1,4) 
               dexdisp(2,lloc) = dexdisp(2,lloc) +e*dtaper*drij(2,4) 
               dexdisp(3,lloc) = dexdisp(3,lloc) +e*dtaper*drij(3,4) 
               ixlploc = loc(ixlp(jlp))
               dexdisp(1,ixlploc) = dexdisp(1,ixlploc) + 
     $          e*dtaper*drij(1,5) 
               dexdisp(2,ixlploc) = dexdisp(2,ixlploc) + 
     $          e*dtaper*drij(2,5) 
               dexdisp(3,ixlploc) = dexdisp(3,ixlploc) + 
     $          e*dtaper*drij(3,5) 
               izlploc = loc(izlp(jlp))
               dexdisp(1,izlploc) = dexdisp(1,izlploc) + 
     $          e*dtaper*drij(1,6) 
               dexdisp(2,izlploc) = dexdisp(2,izlploc) + 
     $          e*dtaper*drij(2,6) 
               dexdisp(3,izlploc) = dexdisp(3,izlploc) + 
     $          e*dtaper*drij(3,6) 

               e = e*taper
               exdisp = exdisp + e
             end do
         end do
      end do

      return
      end
