c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine epolar1tcg1  --  Ewald polarization energy            ##
c     ##                              and grad with tcg1                   ##
c     ##                                                                   ##
c     #######################################################################
c
c     "epolar1tcg2" calculates the induced dipole polarization gradients with the tcg2 method
c  
! Computes TCG2 dipoles with NO refinement (except the preconditioner if
! "precond" is TRUE)
! Allows for minimal number of fft and matvec
! # of matvec : 4 (efi, Tr0, T2r0, T3r0)
! # of ffts : 5 (efi, r0, Tr0, T2r0, T3r0)
      subroutine epolar1tcg2_2
      use atmlst
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use mpi
      implicit none
      logical :: precond
      integer nrhs

      integer :: betac, kk, i,
     $           j, iipole, irhs, k,
     $           ierr,iglob,ilocrec,kglob,kkpole,
     $           kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, term, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1






      real*8, allocatable, dimension(:,:,:,:) ::adte
      real*8, allocatable, dimension(:,:,:) :: denedt,denedtbis,
     $            mu_tcg2,mu_tcg2bis, efi,efibis,
     $            Tr0,Tr0rec,Tr0bis,Tr0recbis,T2r0,T2r0rec,T2r0bis,
     $           T2r0recbis,r0,r0bis, 
     $           T3r0,T3r0bis,T3r0rec,T3r0recbis,Tefi,Tefibis,
     $            oldr0,oldTr0, oldefi,ade, adme, 
     $            fphiTr0, fphit2r0, fphit3r0, fphiTE,fphir0,fphiE,
     $           adtebis
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu, 
     $            torq_tbis,torq_mubis,
     $            denedmu,denedmubis, denedr,denedrbis,cphi,
     $            arr_ded, arr_dedbis,arr_dep,arr_depbis,adtb,adtbbis,
     $            arr_dtr0, arr_dtr0bis,arr_dTTr0,arr_dttr0bis,
     $            arr_dTT2r0,arr_dtt2r0bis, 
     $            fphiarrdep, fphiarrded, fphiarrdtr0,
     $            fphiarrdttr0, fphiarrdtt2r0,adebis,admebis
      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi
      real*8:: time0
      real*8 :: xx(1)
      parameter (nrhs=2)

      if (rank.le.ndir-1) then
        allocate(adte(3,3,nrhs,npolebloc))
        allocate(denedt(3,3,npolebloc),efi(3,nrhs,max(1,npolebloc)))
        allocate(mu_tcg2(3,nrhs, npolebloc),Tr0(3,nrhs,npolebloc),
     $ T2r0(3,nrhs,npolebloc),Tr0rec(3,nrhs,npoleloc),
     $ T2r0rec(3,nrhs,npoleloc),r0(3,nrhs,npolebloc),
     $ T3r0(3,nrhs,npolebloc),T3r0rec(3,nrhs,npoleloc),
     $ Tefi(3,nrhs,npolebloc),ade(3,nrhs,npolebloc),
     $ adme(3,nrhs,npolebloc),oldr0(3,nrhs,npolebloc), 
     $ oldTr0(3,nrhs,npolebloc),oldefi(3,nrhs,npolebloc),
     $ cphi(10, npoleloc),torq_mu(3,nbloc),torq_t(3,nbloc),
     $ denedmu(3,npolebloc),denedr(3,npolebloc),arr_ded(3,npolebloc),
     $ arr_dep(3,npolebloc),adtb(3,npolebloc),arr_dtr0(3,npolebloc),
     $ arr_dTTr0(3,npolebloc),arr_dTT2r0(3,npolebloc))
        allocate (buffermpi(10,max(npoleloc,1)))
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        !1. prepare arrays...
        cphi = 0d0
        efi = 0d0
      else
        allocate(Tr0recbis(3,nrhs,npolerecloc),
     $ T2r0recbis(3,nrhs,npolerecloc),r0bis(3,nrhs,npolerecloc),
     $ T3r0recbis(3,nrhs,npolerecloc),Tefibis(3,nrhs,npolerecloc),
     $ efibis(3,nrhs,npolerecloc),Tr0bis(3,nrhs,npolerecloc),
     $ T2r0bis(3,nrhs,npolerecloc),T3r0bis(3,nrhs,npolerecloc),
     $ mu_tcg2bis(3,nrhs,npolerecloc),adebis(3,nlocrec),
     $ admebis(3,nlocrec),torq_mubis(3,nlocrec2),torq_tbis(3,nlocrec2),
     $ denedrbis(3,nlocrec),denedmubis(3,nlocrec),
     $ arr_dedbis(3,npolerecloc),arr_depbis(3,npolerecloc),
     $ adtbbis(3,nlocrec),arr_dtr0bis(3,npolerecloc),
     $ arr_dTTr0bis(3,npolerecloc),arr_dTT2r0bis(3,npolerecloc))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        allocate(adtebis(3,3,nlocrec))
        allocate(denedtbis(3,3,npolerecloc))
        allocate(fphitr0(20,2,max(1,npolerecloc)),
     $         fphit2r0(20,2,max(1,npolerecloc)),
     $         fphit3r0(20,2,max(1,npolerecloc)),
     $         fphitE(20,2,max(1,npolerecloc)),
     $         fphir0(20,2,max(1,npolerecloc)),
     $         fphiE(20,2,npolerecloc))
        allocate(fphiarrded(20,npolerecloc),
     $         fphiarrdep(20,npolerecloc), 
     $         fphiarrdtr0(20,npolerecloc),
     $         fphiarrdttr0(20,npolerecloc),
     $         fphiarrdtt2r0(20,npolerecloc))
      end if

      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))

      time0 = mpi_wtime()
      precond = tcgprec
      
      f = electric/dielec

      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,efi,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs, efi)
        call commfield(nrhs, efi)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            efi(j,1,i) = efi(j,1,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
            efi(j,2,i) = efi(j,2,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
          end do
        end do
        call commdirdir(nrhs,1,efi,reqrec,reqsend)
        call commdirdir(nrhs,2,efi,reqrec,reqsend)
        oldefi = efi
      else
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
      end if

      if (rank.le.ndir-1) then
        call commrecdirdip(nrhs,1,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (precond) then
          call diagvec(nrhs, efi, efi)
        end if
        oldr0 = oldefi
        r0 = efi

        call tmatxb_pme(nrhs, .true., r0, Tr0)
        call commfield(nrhs,Tr0)

        call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $    + Tr0rec(:,:,1:npoleloc) -term*r0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldTr0 = Tr0
        if (precond) then
          call diagvec(nrhs, Tr0 ,Tr0)
        end if
        Tefi = Tr0
      else
        call commrecdirdip(nrhs,0,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            efibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
          end do
        end if
        r0bis = efibis

        call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,fphir0)

        call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            Tr0bis(:,:,i) = polarity(iipole)*Tr0bis(:,:,i)
          end do
        end if

        fphiE = fphir0
        Tefibis = Tr0bis
      end if

      if (rank.le.ndir-1) then
        do irhs = 1,nrhs
           n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
           t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1

        call tmatxb_pme(nrhs, .true., Tr0, T2r0)
        call commfield(nrhs,T2r0)

        call commrecdirsolv(nrhs,0,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        T2r0(:,:,1:npoleloc) = T2r0(:,:,1:npoleloc) 
     $   + T2r0rec(:,:,1:npoleloc) - term*Tr0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T2r0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        if (precond) then
          call diagvec(nrhs, T2r0, T2r0)
        end if

        do k = 1,nrhs
           np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2(:,k,:) =  (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $     - g1(k)*t4(k)*Tr0(:,k,:)
        end do

        sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
        sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1

        call tmatxb_pme(nrhs, .true., T2r0, T3r0)
        call commfield(nrhs,T3r0)

        call commrecdirsolv(nrhs,0,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T3r0(:,:,1:npoleloc) = T3r0(:,:,1:npoleloc) 
     $   + T3r0rec(:,:,1:npoleloc) - term*T2r0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T3r0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
         if (precond) then
           call diagvec(nrhs, T3r0, T3r0)
         end if
      else
        n0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1

        call tmatxbrecipsave(Tr0,Tr0bis,nrhs,T2r0rec,T2r0recbis,fphiTr0)

        call commrecdirsolv(nrhs,1,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            T2r0bis(:,:,i) = polarity(iipole)*T2r0bis(:,:,i)
          end do
        end if
        fphiTE = fphiTr0

        do k = 1,nrhs
           np1(k) = 0d0 
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = 0d0 
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2bis(:,k,:) =  (g1(k)*t2(k) + t4(k))*r0bis(:,k,:) 
     $     - g1(k)*t4(k)*Tr0bis(:,k,:)
        end do
        sp0 = 0d0 
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = 0d0 
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1

        call tmatxbrecipsave(T2r0,T2r0bis,nrhs,T3r0rec,T3r0recbis,
     $                        fphit2r0)
        call commrecdirsolv(nrhs,1,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            T3r0bis(:,:,i) = polarity(iipole)*T3r0bis(:,:,i)
          end do
        end if
        call fftthatplz2(T3r0,T3r0bis,fphit3r0)
      end if

      a10 = t4(1) + g1(1)*t2(1)
      a1m1= -g1(1)*t4(1)
      a11 = 2d0*b1/t1(1) - 2d0*np1(1)*b2/t3(1) 
     $      - 2d0*np1(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1)*t1(1))
     $      + 2d0*t9(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + 2d0*np1(1)*sp0*g1(1)/(t1(1)*t1(1))
      a12 = -2d0*n0(1)*b1/(t1(1)*t1(1)) + 4d0*t1(1)*b2/t3(1) 
     $      - 2d0*n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 4d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 4d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a13 = -4d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0* n0(1)*b2/t3(1)
     $      + 2d0*n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a14 = 2d0*t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a20 = -g1(1)*t4(1)
      a21 = -n0(1)*b1/(t1(1)*t1(1)) + 2d0*t1(1)*b2/t3(1) 
     $      - n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 2d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a22 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a23 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a31 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*g1(1)*sp0/(t1(1)*t1(1))
      a41 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a32 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))


      if (rank.le.ndir-1) then
        arr_dEp  = mu_tcg2(:,1,:)
        arr_dEd = a10*efi(:,2,:) 
     $              + a1m1*Tefi(:,2,:)
     $              + a11*r0(:,1,:)
     $              + a12*Tr0(:,1,:)
     $              + a13*T2r0(:,1,:) 
     $              + a14*T3r0(:,1,:)
        arr_dTr0 = a20*efi(:,2,:)
     $              + a21*r0(:,1,:)
     $              + a22*Tr0(:,1,:)
     $              + a23*T2r0(:,1,:)
        arr_dTTr0 = a31*r0(:,1,:)
     $              + a32*Tr0(:,1,:)
        arr_dTT2r0 = a41*r0(:,1,:)
        ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
        ep = -.5d0*f*ep 

        denedr = 0d0
        denedmu = 0d0
        denedt = 0d0
        call scalderfieldzmat6(
     &                         arr_ded, arr_dep, arr_dtr0, r0(:,1,:),
     &                         arr_dTTr0,  Tr0(:,1,:), arr_dTT2r0,
     &                         T2r0(:,1,:), ade, adme, adte, adtb
     &                        )
        denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb
        denedmu=   adme(:,1,:)   + adme(:,2,:)    
        denedt =  adte(:,:,1,:) + adte(:,:,2,:)
        call torquetcg_dir(torq_mu,torq_t,denedmu,denedt)
        dep = 0d0

        do kk = 1, npolebloc
           kkpole = poleglob(kk)
           kglob = ipole(kkpole) 
           kloc = loc(kglob)
           do betac = 1, 3
              dep(betac,kloc) =  -.5d0*denedr(betac,kk)*f
           end do
        end do
        dep(:,:) = dep(:,:)-0.5d0*f*(torq_mu(:,:)+torq_t(:,:))

      else
        arr_dEpbis  = mu_tcg2bis(:,1,:)
        arr_dEdbis = a10*efibis(:,2,:) 
     $              + a1m1*Tefibis(:,2,:)
     $              + a11*r0bis(:,1,:)
     $              + a12*Tr0bis(:,1,:)
     $              + a13*T2r0bis(:,1,:) 
     $              + a14*T3r0bis(:,1,:)
        arr_dTr0bis = a20*efibis(:,2,:)
     $              + a21*r0bis(:,1,:)
     $              + a22*Tr0bis(:,1,:)
     $              + a23*T2r0bis(:,1,:)
        arr_dTTr0bis = a31*r0bis(:,1,:)
     $              + a32*Tr0bis(:,1,:)
        arr_dTT2r0bis = a41*r0bis(:,1,:)
        fphiarrdep = (t4(1)+g1(1)*t2(1))*fphir0(:,1,:) 
     $                  - g1(1)*t4(1)*fphiTr0(:,1,:)
        fphiarrded = a10*fphiE(:,2,:)
     $              + a1m1*fphiTe(:,2,:)
     $              + a11*fphir0(:,1,:)
     $              + a12*fphiTr0(:,1,:)
     $              + a13*fphiT2r0(:,1,:)
     $              + a14*fphiT3r0(:,1,:)
        fphiarrdTr0 = a20*fphie(:,2,:)
     $              + a21*fphir0(:,1,:)
     $              + a22*fphiTr0(:,1,:)
     $              + a23*fphiT2r0(:,1,:)
        fphiarrdTTr0  = a31*fphir0(:,1,:)
     $               + a32*fphiTr0(:,1,:)
        fphiarrdTT2r0  = a41*fphir0(:,1,:)

        call scalderfieldzmatrec6(
     &                            xx,arr_depbis,fphiarrdep,xx,
     &                            arr_dedbis,fphiarrded, xx,arr_dtr0bis,
     &                            fphiarrdtr0,xx,r0bis(:,1,:),
     &                            fphir0(:,1,:),xx,arr_dttr0bis,
     &                            fphiarrdttr0, xx,tr0bis(:,1,:),
     &                            fphitr0(:,1,:), xx,arr_dtt2r0bis,
     &                            fphiarrdtt2r0,xx, t2r0bis(:,1,:),
     &                            fphit2r0(:,1,:),adebis, admebis,
     &                            adtebis, adtbbis
     &                           )

        denedrbis = adebis + adtbbis
        denedmubis =  admebis
        denedtbis = adtebis
        call torquetcg_rec(torq_mubis,torq_tbis,denedmubis,denedtbis)

        deprec = 0d0

        do kk = 1, npolerecloc
           iipole = polerecglob(kk)
           iglob = ipole(iipole)
           ilocrec = locrec(iglob)
           do betac = 1, 3
              deprec(betac,ilocrec) =  -.5d0*f*denedrbis(betac,kk)
           end do
        end do
        deprec(:,:) = deprec(:,:)-0.5d0*f*(torq_mubis(:,:)+
     $   torq_tbis(:,:))
      end if
      return
      end

! Only guess is available (+ prec via 'precond')

      subroutine epolar1tcg2bis_2
      use atmlst
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use mpi
      implicit none

      logical :: precond
      integer :: betac, kk, i,
     $           j, iipole, irhs, k,
     $           ierr,iglob,ilocrec,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, term, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, 
     $          omega

      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1




      real*8, allocatable, dimension(:,:,:,:) ::adte
      real*8, allocatable, dimension(:,:,:) :: denedt, denedtbis,
     $            mu_tcg2,efi, efibis,
     $            tmpbig,tmpbigbis, oldr0, Tr0,Tr0bis, oldTr0, T2r0,
     $            T2r0bis, mu0,mu0bis, r0,r0bis,
     $            T3r0, T3r0bis,oldT3r0, 
     $            Tefi, Tefibis,oldTefi, ade, adme,
     $            adtebis, oldaefi,Tmu0rec,Tmu0recbis,
     $            oldefi,oldefibis,Tefirec,T3r0rec,T3r0recbis,
     $            fphiTr0, fphit2r0, fphit3r0, fphiTE, 
     $            fphir0, fphiE, fphiaE, fphimu0,
     $            tmpbig2,Tr0rec,Tr0recbis,
     $            T2r0rec,T2r0recbis,Tefirecbis,
     $            tmpbig2rec,tmpbig2recbis
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu, torq_tbis,
     $            torq_mubis,
     $            denedmu, denedmubis, denedr, denedrbis,cphi,
     $            arr_ded, arr_dedbis,arr_dep, arr_depbis, 
     $            adtb, adtbbis, 
     $            fphiarrded, fphiarrdep, fphiarrdtr0, 
     $            fphiarrdtt2r0, fphiatarrdr0,fphiarrdttr0,
     $            fphiarrA, fphiaTA,
     $            arr_dTr0, arr_dTr0bis,arr_dTTr0, arr_dTTr0bis,
     $            arr_dTT2r0, arr_dTT2r0bis, arrA, arrAbis, arrTA, 
     $            arrTAbis,arraTA,arraTAbis,admebis,adebis
      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi
      real*8 :: xx(1)
      parameter (nrhs=2)

      if (rank.le.ndir-1) then
        allocate(adte(3,3,nrhs, npolebloc))
        allocate(denedt(3,3,npolebloc),mu_tcg2(3,nrhs, npolebloc),
     $ efi(3,nrhs, npolebloc),tmpbig(3,nrhs,npolebloc),
     $ oldr0(3,nrhs,npolebloc),Tefirec(3,nrhs,npoleloc),
     $ Tr0(3,nrhs,npolebloc),oldTr0(3,nrhs,npolebloc),
     $ Tr0rec(3,nrhs,npoleloc),T2r0rec(3,nrhs,npoleloc),
     $ T2r0(3,nrhs,npolebloc),mu0(3,nrhs,npolebloc),
     $ r0(3,nrhs,npolebloc),T3r0(3,nrhs,npolebloc),
     $ T3r0rec(3,nrhs,npoleloc),oldT3r0(3,nrhs,npolebloc),
     $ Tefi(3,nrhs,npolebloc),oldTefi(3,nrhs,npolebloc),
     $ ade(3,nrhs,npolebloc),adme(3,nrhs,npolebloc),
     $ oldaefi(3,nrhs,npolebloc),oldefi(3,nrhs,npolebloc), 
     $ tmpbig2(3,nrhs,npolebloc),Tmu0rec(3,nrhs,npoleloc),
     $ tmpbig2rec(3,nrhs,npoleloc),cphi(10, npoleloc),
     $ torq_mu(3,nbloc),torq_t(3,nbloc), denedmu(3,npolebloc),
     $ denedr(3,npolebloc), arr_ded(3,npolebloc),
     $ arr_dep(3,npolebloc),adtb(3,npolebloc),arr_dtr0(3,npolebloc), 
     $ arr_dttr0(3,npolebloc),arr_dtt2r0(3,npolebloc),arrA(3,npolebloc),
     $ arrTA(3,npolebloc), arraTA(3,npolebloc))
       allocate (buffermpi(10,max(npoleloc,1)))
       allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
       !1. prepare arrays...
       cphi = 0d0
       efi = 0d0
       Tr0 = 0d0
      else
        allocate(adtebis(3,3,nlocrec),denedtbis(3,3,npolerecloc),
     $ efibis(3,nrhs,npolerecloc),tmpbigbis(3,nrhs,npolerecloc),
     $ Tefirecbis(3,nrhs,npolerecloc),Tr0bis(3,nrhs,npolerecloc),
     $ Tr0recbis(3,nrhs,npolerecloc),T2r0recbis(3,nrhs,npolerecloc),
     $ T2r0bis(3,nrhs,npolerecloc),mu0bis(3,nrhs,npolerecloc),
     $ r0bis(3,nrhs,npolerecloc),T3r0bis(3,nrhs,npolerecloc),
     $ T3r0recbis(3,nrhs,npolerecloc),Tefibis(3,nrhs,npolerecloc),
     $ admebis(3,nlocrec),oldefibis(3,nrhs,npolerecloc),
     $ Tmu0recbis(3,nrhs,npolerecloc),tmpbig2recbis(3,nrhs,npolerecloc),
     $ torq_tbis(3,nlocrec2),torq_mubis(3,nlocrec2),
     $ denedmubis(3,npolerecloc),denedrbis(3,npolerecloc),
     $ arr_dedbis(3,npolerecloc),arr_depbis(3,npolerecloc),
     $ adebis(3,nlocrec),adtbbis(3,npolerecloc),
     $ arr_dtr0bis(3,npolerecloc),arr_dttr0bis(3,npolerecloc),
     $ arr_dtt2r0bis(3,npolerecloc),arrAbis(3,npolerecloc),
     $ arrTAbis(3,npolerecloc), arraTAbis(3,npolerecloc))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        allocate(fphitr0(20,2,max(1,npolerecloc)),
     $           fphit2r0(20,2,max(1,npolerecloc)),
     $           fphit3r0(20,2,max(1,npolerecloc)),
     $           fphitE(20,2,max(1,npolerecloc)),
     $           fphir0(20,2,max(1,npolerecloc)),
     $           fphiE(20,2,npolerecloc),
     $           fphiaE(20,2,npolerecloc),
     $           fphimu0(20,2,npolerecloc))
        allocate(fphiarrded(20, npolerecloc), 
     $           fphiarrdep(20, npolerecloc),
     $           fphiarrdtr0(20,npolerecloc),
     $           fphiarrdtt2r0(20,npolerecloc), 
     $           fphiatarrdr0(20,npolerecloc),
     $           fphiarrdttr0(20,npolerecloc), 
     $           fphiarrA(20,npolerecloc), 
     $           fphiaTA(20,npolerecloc)) 
      end if
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      
      precond = tcgprec
      f = electric/dielec


      omega = tcgomega

      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,efi,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs, efi)
        call commfield(nrhs, efi)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            efi(j,1,i) = efi(j,1,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
            efi(j,2,i) = efi(j,2,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
          end do
        end do
        call commdirdir(nrhs,1,efi,reqrec,reqsend)
        call commdirdir(nrhs,2,efi,reqrec,reqsend)
        oldefi = efi
      else
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
      end if

      if (rank.le.ndir-1) then
        call commrecdirdip(nrhs,1,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (precond) then
          call diagvec(nrhs, efi, efi)
        end if
        call diagvec(nrhs, oldefi, mu0)
        call dbletmatxb_pme(nrhs, .true., efi, mu0, Tefi, tmpbig)
        call commfield(nrhs,Tefi)
        call commfield(nrhs,tmpbig)
        call commrecdirsolv(nrhs,0,Tmu0recbis,Tmu0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,0,Tefirecbis,Tefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        tmpbig(:,:,1:npoleloc) = tmpbig(:,:,1:npoleloc) + 
     $  Tmu0rec(:,:,1:npoleloc) - term*mu0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,tmpbig,reqrec,reqsend)
        call commdirdir(nrhs,1,tmpbig,reqrec,reqsend)
        call commdirdir(nrhs,2,tmpbig,reqrec,reqsend)

        call commrecdirdip(nrhs,1,tmpbigbis,tmpbig,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,tmpbigbis,tmpbig,buffermpimu,
     $ buffermpimu,req2rec,req2send)

        Tefi(:,:,1:npoleloc) = Tefi(:,:,1:npoleloc) + 
     $ Tefirec(:,:,1:npoleloc) - term*efi(:,:,1:npoleloc)

        call commdirdir(nrhs,0,Tefi,reqrec,reqsend)
        call commdirdir(nrhs,1,Tefi,reqrec,reqsend)
        call commdirdir(nrhs,2,Tefi,reqrec,reqsend)

        call commrecdirdip(nrhs,1,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        oldTefi = Tefi
        oldaefi =  mu0
        r0 = oldefi - tmpbig
        oldr0 = r0
        if (precond) then
          call diagvec2(nrhs, Tefi, r0, Tefi, r0)
        end if

        call tmatxb_pme(nrhs, .true., r0, Tr0)
        call commfield(nrhs,Tr0)
        call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) + 
     $     Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

        call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldTr0 = Tr0
        if (precond) then
          call diagvec(nrhs, Tr0, Tr0)
        end if
      else
        call commrecdirdip(nrhs,0,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldefibis = efibis
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            efibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
          end do
        end if
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          mu0bis(:,:,i) = polarity(iipole)*oldefibis(:,:,i)
        end do
        call tmatxbrecipsave(mu0, mu0bis, nrhs, Tmu0rec, Tmu0recbis,
     $                      fphiaE)
        call commrecdirsolv(nrhs,1,Tmu0recbis,Tmu0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call tmatxbrecipsave(efi, efibis, nrhs, Tefirec, Tefirecbis, 
     $                        fphiE)
        call commrecdirsolv(nrhs,1,Tefirecbis,Tefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call commrecdirdip(nrhs,0,tmpbigbis,tmpbig,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,tmpbigbis,tmpbig,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,0,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        r0bis = oldefibis - tmpbigbis
        fphimu0 = fphiae
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            r0bis(:,:,i) = polarity(iipole)*r0bis(:,:,i)
            Tefibis(:,:,i) = polarity(iipole)*Tefibis(:,:,i)
          end do
        end if
        call fftthatplz2(Tefi,Tefibis, fphiTE)

        call tmatxbrecipsave(r0, r0bis, nrhs, Tr0rec, Tr0recbis, 
     $                      fphir0)
        call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            Tr0bis(:,:,i) = polarity(iipole)*Tr0bis(:,:,i)
          end do
        end if
      end if



      if (rank.le.ndir-1) then
        do irhs = 1,nrhs
           n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
           t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
        sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)

        call tmatxb_pme(nrhs, .true., Tr0, T2r0)
        call commfield(nrhs,T2r0)

        call commrecdirsolv(nrhs,0,T2r0recbis,T2r0rec,buffermpi,
     $  buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $  buffermpi,reqrecdirrec,reqrecdirsend)

        T2r0(:,:,1:npoleloc) = T2r0(:,:,1:npoleloc) + 
     $ T2r0rec(:,:,1:npoleloc)- term*Tr0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T2r0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        if (precond) then
          call diagvec(nrhs, T2r0, T2r0)
        end if
        do k = 1,nrhs
           np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2(:,k,:) = mu0(:,k,:) + (g1(k)*t2(k) + t4(k))*r0(:,k,:)
     $     - g1(k)*t4(k)*Tr0(:,k,:)
        end do
        sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
        sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1

        call tmatxb_pme(nrhs, .true., T2r0, T3r0)
        call commfield(nrhs,T3r0)

        call commrecdirsolv(nrhs,0,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T3r0(:,:,1:npoleloc) = T3r0(:,:,1:npoleloc) 
     $   + T3r0rec(:,:,1:npoleloc) - term*T2r0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T3r0,reqrec,reqsend)
        oldT3r0 = T3r0

        call commrecdirdip(nrhs,1,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (precond) then
          call diagvec(nrhs, T3r0, T3r0)
        end if
      else
        n0 = 0d0
        t1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        sp0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
        sp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)

        call tmatxbrecipsave(Tr0, Tr0bis, nrhs, T2r0rec,T2r0recbis,
     $                      fphiTr0)
        call commrecdirsolv(nrhs,1,T2r0recbis,T2r0rec,buffermpi,
     $  buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $  buffermpi,reqrecdirrec,reqrecdirsend)

        call commrecdirdip(nrhs,0,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            T2r0bis(:,:,i) = polarity(iipole)*T2r0bis(:,:,i)
          end do
        end if 
        do k = 1,nrhs
           np1(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
        end do
        sp0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
        sp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1

        call tmatxbrecipsave(T2r0, T2r0bis, nrhs, T3r0rec,T3r0recbis,
     $                      fphit2r0)
        call commrecdirsolv(nrhs,1,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        call commrecdirdip(nrhs,0,T3r0bis,T3r0,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            T3r0bis(:,:,i) = polarity(iipole)*T3r0bis(:,:,i)
          end do
        end if
        call fftthatplz2(T3r0,T3r0bis, fphiT3r0)
      end if

      a10 = t4(1) + g1(1)*t2(1)
      a1m1= -g1(1)*t4(1)
      a11 = 2d0*b1/t1(1) - 2d0*np1(1)*b2/t3(1) 
     $      - 2d0*np1(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1)*t1(1))
     $      + 2d0*t9(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + 2d0*np1(1)*sp0*g1(1)/(t1(1)*t1(1))
      a12 = -2d0*n0(1)*b1/(t1(1)*t1(1)) + 4d0*t1(1)*b2/t3(1) 
     $      - 2d0*n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 4d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 4d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a13 = -4d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0* n0(1)*b2/t3(1)
     $      + 2d0*n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a14 = 2d0*t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a20 = -g1(1)*t4(1)
      a21 = -n0(1)*b1/(t1(1)*t1(1)) + 2d0*t1(1)*b2/t3(1) 
     $      - n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 2d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a22 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a23 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a31 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*g1(1)*sp0/(t1(1)*t1(1))
      a41 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a32 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))




      if (rank.le.ndir-1) then
        arr_dep =  a10*r0(:,1,:) 
     $              + a20*Tr0(:,1,:)
        arr_dEd = a10*efi(:,2,:) 
     $              + a1m1*Tefi(:,2,:)
     $              + a11*r0(:,1,:)
     $              + a12*Tr0(:,1,:)
     $              + a13*T2r0(:,1,:) 
     $              + a14*T3r0(:,1,:)
        arr_dTr0 = a20*efi(:,2,:)
     $              + a21*r0(:,1,:)
     $              + a22*Tr0(:,1,:)
     $              + a23*T2r0(:,1,:)
        arr_dTTr0 = a31*r0(:,1,:)
     $              + a32*Tr0(:,1,:)
        arr_dTT2r0 = a41*r0(:,1,:)
        ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
        ep = -.5d0*f*ep 
      else
        arr_depbis =  a10*r0bis(:,1,:) 
     $              + a20*Tr0bis(:,1,:)
        arr_dEdbis = a10*efibis(:,2,:) 
     $              + a1m1*Tefibis(:,2,:)
     $              + a11*r0bis(:,1,:)
     $              + a12*Tr0bis(:,1,:)
     $              + a13*T2r0bis(:,1,:) 
     $              + a14*T3r0bis(:,1,:)
        arr_dTr0bis = a20*efibis(:,2,:)
     $              + a21*r0bis(:,1,:)
     $              + a22*Tr0bis(:,1,:)
     $              + a23*T2r0bis(:,1,:)
        arr_dTTr0bis = a31*r0bis(:,1,:)
     $              + a32*Tr0bis(:,1,:)
        arr_dTT2r0bis = a41*r0bis(:,1,:)
        fphiarrdep = (t4(1)+g1(1)*t2(1))*fphir0(:,1,:) 
     $                 - g1(1)*t4(1)*fphiTr0(:,1,:)
        fphiarrdEd = a10*fphiE(:,2,:)
     $              + a1m1*fphiTe(:,2,:)
     $              + a11*fphir0(:,1,:)
     $              + a12*fphiTr0(:,1,:)
     $              + a13*fphiT2r0(:,1,:)
     $              + a14*fphiT3r0(:,1,:)
        fphiarrdTr0 = a20*fphie(:,2,:)
     $              + a21*fphir0(:,1,:)
     $              + a22*fphiTr0(:,1,:)
     $              + a23*fphiT2r0(:,1,:)
        fphiarrdTTr0 = a31*fphir0(:,1,:)
     $               + a32*fphiTr0(:,1,:)
        fphiarrdTT2r0 = a41*fphir0(:,1,:)
        ep = 0d0
      end if



      if (rank.le.ndir-1) then
        denedr = 0d0
        denedmu = 0d0
        denedt = 0d0
        arrA = arr_ded 
        tmpbig(:,1,:) = arrA
        tmpbig(:,2,:) = arrA
        call tmatxb_pme(nrhs, .true., tmpbig, tmpbig2)
        call commfield(nrhs,tmpbig2)
        call commrecdirsolv(nrhs,0,tmpbig2recbis,tmpbig2rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,tmpbig2recbis,tmpbig2rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        arrTA(:,1:npoleloc) = tmpbig2(:,1,1:npoleloc) + 
     $   tmpbig2rec(:,1,1:npoleloc) -term*arrA(:,1:npoleloc)

        call commdirdir(1,0,arrTA,reqrec,reqsend)
        call commdirdir(1,1,arrTA,reqrec,reqsend)
        call commdirdir(1,2,arrTA,reqrec,reqsend)

        call diagvec(1, arrTA, arraTA)

        call commrecdirdip(1,1,arrTAbis,arrTA,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(1,2,arrTAbis,arrTA,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        arr_ded = mu0(:,2,:) + arrA - arraTA
        arr_dep = arr_dep + mu0(:,1,:)

        call scalderfieldzmat8(
     &                         arr_ded, arr_dep, arr_dtr0,  r0(:,1,:),
     &                         arr_dTTr0,  Tr0(:,1,:), arr_dTT2r0,
     &                         T2r0(:,1,:), -arrA,  oldaefi(:,1,:),
     &                         ade, adme, adte, adtb
     &                        )

        denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb  
        denedmu=  adme(:,1,:)   + adme(:,2,:)    
        denedt =  adte(:,:,1,:) + adte(:,:,2,:)

        call torquetcg_dir(torq_mu,torq_t,denedmu,denedt)
        dep = 0d0

        do kk = 1, npolebloc
           kkpole = poleglob(kk)
           kglob = ipole(kkpole) 
           kloc = loc(kglob)
           do betac = 1, 3
              dep(betac,kloc) =  -.5d0*denedr(betac,kk)*f
           end do
        end do
        dep(:,:) = dep(:,:)-0.5d0*f*(torq_mu(:,:)+torq_t(:,:))
      else
        arrAbis = arr_dedbis
        tmpbigbis(:,1,:) = arrAbis
        tmpbigbis(:,2,:) = arrAbis
        fphiarrA = fphiarrded
        call tmatxbrecip(tmpbig, tmpbigbis, nrhs, tmpbig2rec,
     $   tmpbig2recbis)
        call commrecdirsolv(nrhs,1,tmpbig2recbis,tmpbig2rec,buffermpi,
     $  buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,tmpbig2recbis,tmpbig2rec,buffermpi,
     $   buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(1,0,arrTAbis,arrTA,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        call commrecdirdip(1,2,arrTAbis,arrTA,buffermpimu,
     $ buffermpimu,req2rec,req2send)
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          arraTabis(:,i) = polarity(iipole)*arrTAbis(:,i)
        end do
        call fftthatplz(arraTA, arraTAbis, fphiaTA)
        arr_dedbis = mu0bis(:,2,:) + arrAbis - arraTAbis
        fphiarrded = fphimu0(:,2,:) + fphiarrded - fphiaTA
        arr_depbis = arr_depbis + mu0bis(:,1,:)
        fphiarrdep = fphiarrdep + fphimu0(:,1,:)

        call scalderfieldzmatrec8(
     &                            xx,arr_dedbis,fphiarrded,xx,
     &                            arr_depbis,fphiarrdep, xx,arr_dtr0bis,
     &                            fphiarrdtr0,xx,r0bis(:,1,:),
     &                            fphir0(:,1,:),xx,arr_dttr0bis,
     &                            fphiarrdttr0, xx,tr0bis(:,1,:),
     &                            fphitr0(:,1,:), xx,arr_dtt2r0bis,
     &                            fphiarrdtt2r0,xx, t2r0bis(:,1,:),
     &                            fphit2r0(:,1,:),xx, - arrAbis,
     &                          - fphiarrA,xx,mu0bis(:,1,:),
     &                            fphimu0(:,1,:), adebis, admebis,
     &                            adtebis,adtbbis
     &                           )

        denedrbis = adebis + adtbbis
        denedmubis = admebis
        denedtbis = adtebis

        call torquetcg_rec(torq_mubis,torq_tbis,denedmubis,denedtbis)

        deprec = 0d0

        do kk = 1, npolerecloc
           iipole = polerecglob(kk)
           iglob = ipole(iipole)
           ilocrec = locrec(iglob)
           do betac = 1, 3
              deprec(betac,ilocrec) =  -.5d0*f*denedrbis(betac,kk)
           end do
        end do
        deprec(:,:)=deprec(:,:)-0.5d0*f*(torq_mubis(:,:)+torq_tbis(:,:))
      end if

      return
      end
! Obj : tcg2 with peek only, preconditioner available if 'precond' is
! TRUE
      subroutine epolar1tcg2ter_2
      use atmlst
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use mpi
      implicit none


      logical :: precond
      integer :: betac, kk, i,
     $           j, iipole, irhs, k,
     $           ierr,iglob,ilocrec,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, term, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, ap10a,
     $          ap11a, ap12a, ap11, ap12, ap13, ap14, ap2a0, ap21a,
     $          ap21, ap22, ap23, ap3a0, ap31, ap32, ap41, spp1, spp2,
     $          omega, k1, k2, k3
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1



      real*8 xx(1)

      real*8, allocatable, dimension(:,:,:,:) ::adte
      real*8, allocatable, dimension(:,:,:) :: denedt,denedtbis,mu_tcg2,
     $            mu_tcg2bis,efi,efibis,adtebis,
     $            tefi,tefibis, oldtefi, tefirec, tefirecbis,r0,
     $            r0bis,t3efi,T3efibis,T3efirec,T3efirecbis, 
     $            atefibis,at2efibis,
     $             oldtefibis,
     $            ade, adme, oldaefi,
     $            oldefi,oldefibis, atefi, at2efi, aefi,aefibis,
     $            Taefi,taefirec, T2aefi,T2aefirec, mupeek,mupeekbis,
     $            T2efi,
     $            Taefirecbis,T2efibis,T2aefibis,T2aefirecbis,
     $            T2efirec,T2efirecbis,fphiTaE, 
     $            Taefibis,Tr0rec,Tr0recbis,
     $            fphiT2ae,fphiTE, fphiT2E,
     $            fphiE, fphit3e, fphiaE,
     $            fphiate, fphiat2e,
     $            oldt2efibis
      real*8, allocatable, dimension(:,:) ::
     $            arr_dTr0,
     $            arr_dTr0bis,arr_dTTr0bis,arr_dTT2r0bis,
     $            arr_dTTr0, arr_dTT2r0,
     $            fphiarr_ded,
     $            fphiarr_dep, 
     $            fphiarrdtr0, fphiarrdttr0,
     $            fphiarrdtt2r0
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu,
     $            denedmu, denedr, cphi, 
     $            arr_ded, arr_dep, adtb,
     $            arr_dedbis, arr_depbis,torq_mubis,torq_tbis,
     $            denedmubis,denedrbis,
     $            adebis,admebis,adtbbis
      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi
      parameter (nrhs=2)


      if (rank.le.ndir-1) then
        allocate(adte(3,3,nrhs, npolebloc))
        allocate(denedt(3,3,npolebloc),mu_tcg2(3,nrhs,npolebloc),
     $ efi(3,nrhs, npolebloc),tefi(3,nrhs,npolebloc), 
     $ oldtefi(3,nrhs,npolebloc),Taefirec(3,nrhs,npoleloc),
     $ Tefirec(3,nrhs,npoleloc),T2efirec(3,nrhs,npoleloc),
     $ t2efi(3,nrhs,npolebloc),r0(3,nrhs,npolebloc),
     $ t3efi(3,nrhs,npolebloc),T3efirec(3,nrhs,npoleloc),
     $ arr_dTr0(3,npolebloc),arr_dTTr0(3,npolebloc),
     $ arr_dTT2r0(3,npolebloc),ade(3,nrhs,npolebloc), 
     $ adme(3,nrhs,npolebloc),oldaefi(3,nrhs,npolebloc),
     $ oldefi(3,nrhs,npolebloc),atefi(3,nrhs,npolebloc),
     $ at2efi(3,nrhs,npolebloc),aefi(3,nrhs,npolebloc),
     $ Taefi(3,nrhs,npolebloc),T2aefi(3,nrhs,npolebloc),
     $ T2aefirec(3,nrhs,npoleloc),Tr0rec(3,nrhs,npoleloc),
     $ mupeek(3,nrhs,npolebloc))
        allocate(cphi(10, npoleloc),torq_mu(3,nbloc),
     $ torq_t(3,nbloc),denedmu(3,npolebloc),denedr(3,npolebloc), 
     $ arr_ded(3,npolebloc),arr_dep(3,npolebloc),
     $ adtb(3,npolebloc))
        allocate (buffermpi(10,max(npoleloc,1)))
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        cphi = 0d0
        efi = 0d0
      else
        allocate(adtebis(3,3,npolerecloc))
        allocate(denedtbis(3,3,npolerecloc),
     $ mu_tcg2bis(3,nrhs,npolerecloc),efibis(3,nrhs,npolerecloc),
     $ atefibis(3,nrhs,npolerecloc),at2efibis(3,nrhs,npolerecloc),
     $ oldt2efibis(3,nrhs,npolerecloc),T2efirecbis(3,nrhs,npolerecloc),
     $ Taefirecbis(3,nrhs,npolerecloc),Tefibis(3,nrhs,npolerecloc),
     $ oldtefibis(3,nrhs,npolerecloc),T2efibis(3,nrhs,npolerecloc),
     $ r0bis(3,nrhs,npolerecloc),T3efibis(3,nrhs,npolerecloc),
     $ T3efirecbis(3,nrhs,npolerecloc),Tefirecbis(3,nrhs,npolerecloc),
     $ arr_dTr0bis(3,npolerecloc),arr_dTTr0bis(3,npolerecloc),
     $ arr_dTT2r0bis(3,npolerecloc),oldefibis(3,nrhs,npolerecloc),
     $ aefibis(3,nrhs,npolerecloc),Taefibis(3,nrhs,npolerecloc),
     $ T2aefibis(3,nrhs,npolerecloc),T2aefirecbis(3,nrhs,npolerecloc),
     $ Tr0recbis(3,nrhs,npolerecloc),mupeekbis(3,nrhs,npolerecloc))
        allocate(fphiTaE(20,nrhs,npolerecloc),
     $ fphit2ae(20,nrhs,npolerecloc),fphitE(20,nrhs,npolerecloc),
     $ fphit2E(20,nrhs,npolerecloc),fphit3E(20,nrhs,npolerecloc),
     $ fphiarr_ded(20,npolerecloc),fphiarr_dep(20,npolerecloc),
     $ fphiE(20,nrhs,npolerecloc),fphiaE(20,nrhs,npolerecloc),
     $ fphiatE(20,nrhs,npolerecloc),fphiat2E(20,nrhs,npolerecloc),
     $ fphiarrdtr0(20,npolerecloc),fphiarrdttr0(20,npolerecloc),
     $ fphiarrdtt2r0(20,npolerecloc))
        allocate(torq_mubis(3,nlocrec2),torq_tbis(3,nlocrec2),
     $ adebis(3,npolerecloc),adtbbis(3,npolerecloc),
     $ admebis(3,npolerecloc),arr_dedbis(3,npolerecloc),
     $ arr_depbis(3,npolerecloc),denedmubis(3,npolerecloc),
     $ denedrbis(3,nlocrec))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
      end if
      
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      
      f = electric/dielec
      precond = tcgprec

      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,efi,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs, efi)
        call commfield(nrhs, efi)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            efi(j,1,i) = efi(j,1,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
            efi(j,2,i) = efi(j,2,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
          end do
        end do
        call commdirdir(nrhs,1,efi,reqrec,reqsend)
        call commdirdir(nrhs,2,efi,reqrec,reqsend)
        oldefi = efi
      else
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
      end if

      omega = tcgomega

      if (rank.le.ndir-1) then
        call commrecdirdip(nrhs,1,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldefi = efi
        r0 = efi
        oldefi = efi
        if (precond) then
          call diagvec(nrhs, r0, r0)
        end if
        call diagvec(nrhs, efi, aefi)
        oldaefi = aefi

        if (.not.precond) then
           call dbletmatxb_pme(nrhs, .true., r0, aefi, tefi, taefi)
           call commfield(nrhs,Tefi)
           call commrecdirsolv(nrhs,0,Tefirecbis,Tefirec,buffermpi,
     $       buffermpi,reqrecdirrec,reqrecdirsend)
           call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $       buffermpi,reqrecdirrec,reqrecdirsend)

           Tefi(:,:,1:npoleloc) = Tefi(:,:,1:npoleloc) +
     $       Tefirec(:,:,1:npoleloc) - term*efi(:,:,1:npoleloc)

           call commdirdir(nrhs,0,Tefi,reqrec,reqsend)
           call commdirdir(nrhs,1,Tefi,reqrec,reqsend)
           call commdirdir(nrhs,2,Tefi,reqrec,reqsend)
           call commrecdirdip(nrhs,1,Tefibis,Tefi,buffermpimu,
     $      buffermpimu,req2rec,req2send)
           call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $      buffermpimu,req2rec,req2send)

           call commfield(nrhs,Taefi)
           call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,buffermpi,
     $       buffermpi,reqrecdirrec,reqrecdirsend)
           call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,buffermpi,
     $       buffermpi,reqrecdirrec,reqrecdirsend)

           Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) +
     $       Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)

           call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
           call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
           call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

           call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu,
     $      buffermpimu,req2rec,req2send)
           call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $      buffermpimu,req2rec,req2send)
        else
           call tmatxb_pme(nrhs, .true., r0, Tefi)
           call commfield(nrhs,Tefi)
           call commrecdirsolv(nrhs,0,Tefirecbis,Tefirec,buffermpi,
     $       buffermpi,reqrecdirrec,reqrecdirsend)
           call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $       buffermpi,reqrecdirrec,reqrecdirsend)
           Tefi(:,:,1:npoleloc) = Tefi(:,:,1:npoleloc) +
     $       Tefirec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)
           call commdirdir(nrhs,0,Tefi,reqrec,reqsend)
           call commdirdir(nrhs,1,Tefi,reqrec,reqsend)
           call commdirdir(nrhs,2,Tefi,reqrec,reqsend)
           call commrecdirdip(nrhs,1,Tefibis,Tefi,buffermpimu,
     $      buffermpimu,req2rec,req2send)
           call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $      buffermpimu,req2rec,req2send)
           Taefi = Tefi
        end if
      else
        call commrecdirdip(nrhs,0,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldefibis = efibis
        r0bis = efibis
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            r0bis(:,:,i) = polarity(iipole)*r0bis(:,:,i)
            aefibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
          end do
        end if

        if (.not.precond) then
          call tmatxbrecipsave(r0, r0bis, nrhs,Tefirec, Tefirecbis, 
     $                         fphie)
          call commrecdirsolv(nrhs,1,Tefirecbis,Tefirec,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $       buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirdip(nrhs,0,Tefibis,Tefi,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $      buffermpimu,req2rec,req2send)

          call tmatxbrecipsave(aefi, aefibis,nrhs,Taefirec,Taefirecbis, 
     $                           fphiaE)
          call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $     buffermpimu,req2rec,req2send)
        else
          call tmatxbrecipsave(r0,r0bis,nrhs,Tefirec,Tefirecbis,
     $    fphiE)
          call commrecdirsolv(nrhs,1,Tefirecbis,Tefirec,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirdip(nrhs,0,Tefibis,Tefi,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          Taefibis = Tefibis
          fphiae = fphiE
        end if
      end if


      if (rank.le.ndir-1) then
        oldtefi = tefi
        if (precond) then
          call diagvec(nrhs, tefi, tefi)
        end if
        do irhs = 1,nrhs
           n0(irhs) = sprod(3*npoleloc, oldefi(:,irhs,:),r0(:,irhs,:))
           t1(irhs) = sprod(3*npoleloc, oldefi(:,irhs,:),tefi(:,irhs,:))
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1

        call dbletmatxb_pme(nrhs, .true., Tefi, Taefi, T2efi, T2aefi)
        call commfield(nrhs,T2efi)

        call commrecdirsolv(nrhs,0,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        
        t2efi(:,:,1:npoleloc) = t2efi(:,:,1:npoleloc) + 
     $    t2efirec(:,:,1:npoleloc) - term*tefi(:,:,1:npoleloc)

        call commdirdir(nrhs,0,T2efi,reqrec,reqsend)
        call commdirdir(nrhs,1,T2efi,reqrec,reqsend)
        call commdirdir(nrhs,2,T2efi,reqrec,reqsend)
        call commrecdirdip(nrhs,1,t2efibis,t2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,t2efibis,t2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        call commfield(nrhs,T2aefi)
        call commrecdirsolv(nrhs,0,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        t2aefi(:,:,1:npoleloc) = t2aefi(:,:,1:npoleloc) + 
     $    t2aefirec(:,:,1:npoleloc) - term*taefi(:,:,1:npoleloc)

        call commdirdir(nrhs,0,T2aefi,reqrec,reqsend)
        call commdirdir(nrhs,1,T2aefi,reqrec,reqsend)
        call commdirdir(nrhs,2,T2aefi,reqrec,reqsend)

        call commrecdirdip(nrhs,1,t2aefibis,t2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,t2aefibis,t2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (precond) then
          call diagvec(nrhs, t2efi, t2efi)
        end if
        do k = 1,nrhs
           np1(k) = sprod(3*npoleloc, oldtefi(:,k,:), tefi(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = sprod(3*npoleloc, oldtefi(:,k,:), t2efi(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2(:,k,:) = (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $     - g1(k)*t4(k)*tefi(:,k,:)
        end do
        sp0 = sprod(3*npoleloc, r0(:,1,:), efi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = sprod(3*npoleloc, tefi(:,1,:), efi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1
      else
        oldtefibis = tefibis
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            Tefibis(:,:,i) = polarity(iipole)*Tefibis(:,:,i)
          end do
        end if
        n0 = 0d0
        t1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1

        call tmatxbrecipsave(tefi, tefibis, nrhs, t2efirec,
     $                      t2efirecbis,fphite)

        call commrecdirsolv(nrhs,1,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,t2efibis,t2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,t2efibis,t2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        call tmatxbrecipsave(taefi, taefibis, nrhs,t2aefirec,
     $                        t2aefirecbis,fphiTae)
        call commrecdirsolv(nrhs,1,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,t2aefibis,t2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,t2aefibis,t2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        call fftthatplz2(T2aefi,T2aefibis, fphit2ae)
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            t2efibis(:,:,i) = polarity(iipole)*t2efibis(:,:,i)
          end do
        end if
        do k = 1,nrhs
           np1(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2bis(:,k,:) = (g1(k)*t2(k) + t4(k))*r0bis(:,k,:) 
     $     - g1(k)*t4(k)*tefibis(:,k,:)
        end do
        sp0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1
      end if

      if (rank.le.ndir-1) then
        call tmatxb_pme(nrhs, .true., t2efi, t3efi)
        call commfield(nrhs,T3efi)
        call commrecdirsolv(nrhs,0,T3efirecbis,T3efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3efirecbis,T3efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T3efi(:,:,1:npoleloc) = T3efi(:,:,1:npoleloc) 
     $   + T3efirec(:,:,1:npoleloc) - term*T2efi(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T3efi,reqrec,reqsend)
        call commdirdir(nrhs,1,T3efi,reqrec,reqsend)
        call commdirdir(nrhs,2,T3efi,reqrec,reqsend)

        call commrecdirdip(nrhs,1,t3efibis,t3efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,t3efibis,t3efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (precond) then
          call diagvec(nrhs, t3efi, t3efi)
        end if

        if (.not.precond) then
          call diagvec3(nrhs, efi, tefi, t2efi, aefi, atefi, at2efi)
          mupeek = omega*aefi 
     $            - omega*(g1(1)*t2(1) + t4(1))*atefi 
     $            + omega*g1(1)*t4(1)*at2efi
        else
          mupeek = omega*r0
     $            - omega*(g1(1)*t2(1) + t4(1))*tefi 
     $            + omega*g1(1)*t4(1)*t2efi
        end if
        mu_tcg2 = mu_tcg2 + mupeek
        efi = r0
      else
        call tmatxbrecipsave(t2efi,t2efibis,nrhs,t3efirec,t3efirecbis,
     $                        fphit2e)
        call commrecdirsolv(nrhs,1,T3efirecbis,T3efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3efirecbis,T3efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,t3efibis,t3efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,t3efibis,t3efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            t3efibis(:,:,i) = polarity(iipole)*t3efibis(:,:,i)
          end do
        end if

        if (.not.precond) then
          do i = 1, npolerecloc
            iipole = polerecglob(i)
            aefibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
            atefibis(:,:,i) = polarity(iipole)*tefibis(:,:,i)
            at2efibis(:,:,i) = polarity(iipole)*t2efibis(:,:,i)
          end do
          call fftthatplz2(atefi,atefibis, fphiate)
          call fftthatplz2(at2efi,at2efibis, fphiat2e)
          mupeekbis = omega*aefibis 
     $            - omega*(g1(1)*t2(1) + t4(1))*atefibis
     $            + omega*g1(1)*t4(1)*at2efibis
        else
          mupeekbis = omega*r0bis
     $            - omega*(g1(1)*t2(1) + t4(1))*tefibis 
     $            + omega*g1(1)*t4(1)*t2efibis
        end if
        mu_tcg2bis = mu_tcg2bis + mupeekbis
        efibis = r0bis
        call fftthatplz2(t3efi,t3efibis, fphit3e)
      end if

      a10 = t4(1) + g1(1)*t2(1)
      a1m1= -g1(1)*t4(1)
      a11 = 2d0*b1/t1(1) - 2d0*np1(1)*b2/t3(1) 
     $      - 2d0*np1(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1)*t1(1))
     $      + 2d0*t9(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + 2d0*np1(1)*sp0*g1(1)/(t1(1)*t1(1))
      a12 = -2d0*n0(1)*b1/(t1(1)*t1(1)) + 4d0*t1(1)*b2/t3(1) 
     $      - 2d0*n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 4d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 4d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a13 = -4d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0* n0(1)*b2/t3(1)
     $      + 2d0*n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a14 = 2d0*t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a20 = -g1(1)*t4(1)
      a21 = -n0(1)*b1/(t1(1)*t1(1)) + 2d0*t1(1)*b2/t3(1) 
     $      - n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 2d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a22 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a23 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a31 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*g1(1)*sp0/(t1(1)*t1(1))
      a41 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a32 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))


      if (rank.le.ndir-1) then
        if (.not. precond) then
           spp1 = sprod(3*npoleloc, oldaefi(:,2,:), tefi(:,1,:))
           spp2 = sprod(3*npoleloc, oldaefi(:,2,:), t2efi(:,1,:))
        else if (precond) then
           spp1 = sprod(3*npoleloc, oldefi(:,2,:), tefi(:,1,:))
           spp2 = sprod(3*npoleloc, oldefi(:,2,:), t2efi(:,1,:))
        end if
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
      else
        if (.not. precond) then
           spp1 = 0d0 
           spp2 = 0d0 
        else if (precond) then
           spp1 = 0d0
           spp2 = 0d0
        end if
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
      end if


      k1 = omega*t4(1)*spp2 - omega*t2(1)*spp1
      k2 = omega*g1(1)*spp2 - omega*spp1
      k3 = -omega*g1(1)*spp1

      ap10a = omega
      ap11a = -omega*(t2(1)*g1(1) + t4(1))
      ap12a = omega*t4(1)*g1(1)
      ap2a0 = -omega*(g1(1)*t2(1) + t4(1))
      ap21a = omega*t4(1)*g1(1)
      ap3a0 = omega*g1(1)*t4(1)

      ap11 = k1*(-2d0*np1(1)/t3(1)
     $         - 2d0*np1(1)*np1(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t9(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(2d0/t1(1))
     $      +k3*(2d0*np1(1)/(t1(1)*t1(1)))
      ap12 = k1*(4d0*t1(1)/t3(1) 
     $        - 2d0*n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $        + 4d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $        - 2d0*t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(-2d0*n0(1)/(t1(1)*t1(1)))
     $      +k3*(-4d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap13 = k1*(-4d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1))
     $         - 2d0*n0(1)/t3(1))
     $      +k3*(2d0*n0(1)/(t1(1)*t1(1)))
      ap14 = k1*(2d0*t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap21 = k1*(2d0*t1(1)/t3(1) 
     $         - n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $         - t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(-n0(1)/(t1(1)*t1(1)))
     $      +k3*(-2d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap22 = k1*(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3*(n0(1)/(t1(1)*t1(1)))
      ap23 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap31 = k1*(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3*(n0(1)/(t1(1)*t1(1)))
      ap32 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap41 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))



      if (rank.le.ndir-1) then
        arr_dEp =(g1(1)*t2(1)+ t4(1))*efi(:,1,:)-g1(1)*t4(1)*tefi(:,1,:)
        arr_dEd = a10*efi(:,2,:)+ a1m1*Tefi(:,2,:)+ a11*efi(:,1,:)
     $              + a12*tefi(:,1,:)
     $              + a13*t2efi(:,1,:) 
     $              + a14*t3efi(:,1,:)
        arr_dTr0 = a20*efi(:,2,:)+ a21*efi(:,1,:)+ a22*tefi(:,1,:)
     $              + a23*t2efi(:,1,:)
        arr_dTTr0 = a31*efi(:,1,:)+ a32*tefi(:,1,:)
        arr_dTT2r0 = a41*efi(:,1,:)
        if (.not.precond) then
          arr_dep = arr_dep+ ap10a*aefi(:,1,:)+ ap11a*atefi(:,1,:) 
     $             + ap21a*at2efi(:,1,:)
          arr_ded = arr_ded 
     $              + ap10a*aefi(:,2,:)+ ap11a*Taefi(:,2,:)
     $              + ap12a*T2aefi(:,2,:)+ ap11*efi(:,1,:)
     $              + ap12*tefi(:,1,:)+ ap13*t2efi(:,1,:)
     $              + ap14*t3efi(:,1,:)
          arr_dTr0 = arr_dTr0 + ap2a0*aefi(:,2,:)+ ap21a*Taefi(:,2,:)
     $              + ap21*efi(:,1,:)+ap22*tefi(:,1,:)+ap23*t2efi(:,1,:)
          arr_dTTr0 = arr_dTTr0+ ap3a0*aefi(:,2,:)+ ap31*efi(:,1,:)
     $              + ap32*tefi(:,1,:)
          arr_dTT2r0 = arr_dTT2r0+ ap41*efi(:,1,:)
        else
          arr_dep = arr_dep + ap10a*efi(:,1,:)+ ap11a*tefi(:,1,:) 
     $             + ap21a*t2efi(:,1,:)
          arr_ded = arr_ded+ ap10a*efi(:,2,:)+ ap11a*Tefi(:,2,:)
     $              + ap12a*T2efi(:,2,:)+ ap11*efi(:,1,:)
     $              + ap12*tefi(:,1,:)+ ap13*t2efi(:,1,:)
     $              + ap14*t3efi(:,1,:)
          arr_dTr0 = arr_dTr0+ ap2a0*efi(:,2,:)
     $              + ap21a*Tefi(:,2,:)+ ap21*efi(:,1,:)
     $              + ap22*tefi(:,1,:)+ ap23*t2efi(:,1,:)
          arr_dTTr0 = arr_dTTr0+ ap3a0*efi(:,2,:)
     $              + ap31*efi(:,1,:)+ ap32*tefi(:,1,:)
          arr_dTT2r0 = arr_dTT2r0+ ap41*efi(:,1,:)
        end if

        ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
        ep = -.5d0*f*ep 

        denedr = 0d0
        denedmu = 0d0
        denedt = 0d0

        call scalderfieldzmat6(
     &                         arr_ded, arr_dep, arr_dtr0, efi(:,1,:),
     &                         arr_dTTr0, Tefi(:,1,:), arr_dTT2r0,
     &                         T2efi(:,1,:), ade, adme, adte, adtb
     &                        )
        denedr =  ade(:,1,:)    + ade(:,2,:)  + adtb     
        denedmu=  adme(:,1,:)   + adme(:,2,:)    
        denedt =  adte(:,:,1,:) + adte(:,:,2,:)

        call torquetcg_dir(torq_mu,torq_t,denedmu,denedt)
        dep = 0d0

        do kk = 1, npolebloc
           kkpole = poleglob(kk)
           kglob = ipole(kkpole) 
           kloc = loc(kglob)
           do betac = 1, 3
              dep(betac,kloc) =  -.5d0*denedr(betac,kk)*f
           end do
        end do
        dep(:,:) = dep(:,:)-0.5d0*f*(torq_mu(:,:)+torq_t(:,:))
      else
        arr_dEpbis =(g1(1)*t2(1)+ t4(1))*efibis(:,1,:) 
     $           - g1(1)*t4(1)*tefibis(:,1,:)
        arr_dEdbis = a10*efibis(:,2,:)+ a1m1*Tefibis(:,2,:)
     $              + a11*efibis(:,1,:)+ a12*tefibis(:,1,:)
     $              + a13*t2efibis(:,1,:)+ a14*t3efibis(:,1,:)
        arr_dTr0bis = a20*efibis(:,2,:)+ a21*efibis(:,1,:)
     $              + a22*tefibis(:,1,:)+ a23*t2efibis(:,1,:)
        arr_dTTr0bis = a31*efibis(:,1,:)+ a32*tefibis(:,1,:)
        arr_dTT2r0bis = a41*efibis(:,1,:)
        fphiarr_dep = (t4(1)+g1(1)*t2(1))*fphie(:,1,:) 
     $              - g1(1)*t4(1)*fphite(:,1,:)
        fphiarr_dEd = a10*fphiE(:,2,:)+ a1m1*fphiTe(:,2,:)
     $              + a11*fphie(:,1,:)+ a12*fphite(:,1,:)
     $              + a13*fphit2e(:,1,:)+ a14*fphit3e(:,1,:)
        fphiarrdTr0 = a20*fphie(:,2,:)+ a21*fphie(:,1,:)
     $              + a22*fphite(:,1,:)+ a23*fphit2e(:,1,:)
        fphiarrdTTr0 = a31*fphie(:,1,:)+ a32*fphite(:,1,:)
        fphiarrdTT2r0 = a41*fphie(:,1,:)
        if (.not.precond) then
          arr_depbis = arr_depbis+ ap10a*aefibis(:,1,:)
     $             + ap11a*atefibis(:,1,:) + ap21a*at2efibis(:,1,:)
          arr_dedbis = arr_dedbis+ ap10a*aefibis(:,2,:)
     $              + ap11a*Taefibis(:,2,:)+ ap12a*T2aefibis(:,2,:)
     $              + ap11*efibis(:,1,:)+ ap12*tefibis(:,1,:)
     $              + ap13*t2efibis(:,1,:)+ ap14*t3efibis(:,1,:)
          arr_dTr0bis = arr_dTr0bis+ ap2a0*aefibis(:,2,:)
     $              + ap21a*Taefibis(:,2,:)+ ap21*efibis(:,1,:)
     $              + ap22*tefibis(:,1,:)+ ap23*t2efibis(:,1,:)
          arr_dTTr0bis = arr_dTTr0bis+ ap3a0*aefibis(:,2,:)
     $              + ap31*efibis(:,1,:)+ ap32*tefibis(:,1,:)
          arr_dTT2r0bis = arr_dTT2r0bis+ ap41*efibis(:,1,:)
          fphiarr_dep = fphiarr_dep+ ap10a*fphiae(:,1,:)
     $              + ap11a*fphiate(:,1,:)+ ap21a*fphiat2e(:,1,:)
          fphiarr_ded = fphiarr_ded+ ap10a*fphiae(:,2,:)
     $              + ap11a*fphiTae(:,2,:)+ ap12a*fphiT2ae(:,2,:)
     $              + ap11*fphie(:,1,:)+ ap12*fphite(:,1,:)
     $              + ap13*fphit2e(:,1,:)+ ap14*fphit3e(:,1,:)
          fphiarrdTr0 = fphiarrdTr0 + ap2a0*fphiae(:,2,:)
     $              + ap21a*fphiTae(:,2,:)+ ap21 *fphie(:,1,:)
     $              + ap22 *fphite(:,1,:)+ ap23 *fphit2e(:,1,:)
          fphiarrdTTr0 = fphiarrdTTr0+ ap3a0*fphiae(:,2,:)
     $              + ap31 *fphie(:,1,:)+ ap32 *fphite(:,1,:)
          fphiarrdTT2r0 = fphiarrdTT2r0+ ap41*fphie(:,1,:)
        else
          arr_depbis = arr_depbis+ ap10a*efibis(:,1,:)
     $             + ap11a*tefibis(:,1,:)+ ap21a*t2efibis(:,1,:)
          arr_dedbis = arr_dedbis+ ap10a*efibis(:,2,:)
     $              + ap11a*Tefibis(:,2,:)+ ap12a*T2efibis(:,2,:)
     $              + ap11*efibis(:,1,:)+ ap12*tefibis(:,1,:)
     $              + ap13*t2efibis(:,1,:)+ ap14*t3efibis(:,1,:)
          arr_dTr0bis = arr_dTr0bis+ ap2a0*efibis(:,2,:)
     $              + ap21a*Tefibis(:,2,:)+ ap21*efibis(:,1,:)
     $              + ap22*tefibis(:,1,:)+ ap23*t2efibis(:,1,:)
          arr_dTTr0bis = arr_dTTr0bis+ ap3a0*efibis(:,2,:)
     $              + ap31*efibis(:,1,:)+ ap32*tefibis(:,1,:)
          arr_dTT2r0bis = arr_dTT2r0bis+ ap41*efibis(:,1,:)
          fphiarr_dep = fphiarr_dep+ ap10a*fphie(:,1,:)
     $              + ap11a*fphite(:,1,:)+ ap21a*fphit2e(:,1,:)
          fphiarr_ded = fphiarr_ded+ ap10a*fphie(:,2,:)
     $              + ap11a*fphiTe(:,2,:)+ ap12a*fphiT2e(:,2,:)
     $              + ap11 *fphie(:,1,:)+ ap12 *fphite(:,1,:)
     $              + ap13 *fphit2e(:,1,:)+ ap14 *fphit3e(:,1,:)
          fphiarrdTr0 = fphiarrdTr0+ ap2a0*fphie(:,2,:)
     $              + ap21a*fphiTe(:,2,:)+ ap21 *fphie(:,1,:)
     $              + ap22 *fphite(:,1,:)+ ap23 *fphit2e(:,1,:)
          fphiarrdTTr0 = fphiarrdTTr0+ ap3a0*fphie(:,2,:)
     $              + ap31 *fphie(:,1,:)+ ap32 *fphite(:,1,:)
          fphiarrdTT2r0 = fphiarrdTT2r0+ ap41*fphie(:,1,:)
        end if
        ep = 0d0

        call scalderfieldzmatrec6( 
     &                            xx,arr_depbis,fphiarr_dep,xx,
     &                            arr_dedbis, fphiarr_ded, xx,
     &                            arr_dtr0bis,fphiarrdtr0,
     &                            xx,efibis(:,1,:), fphie(:,1,:),xx,
     &                            arr_dttr0bis, fphiarrdttr0, xx,
     &                            tefibis(:,1,:), fphite(:,1,:), xx,
     &                            arr_dtt2r0bis, fphiarrdtt2r0,xx,
     &                            t2efibis(:,1,:), fphit2e(:,1,:),
     &                            adebis, admebis, adtebis, adtbbis
     &                           )

        denedrbis = adebis + adtbbis
        denedmubis = admebis
        denedtbis = adtebis

        call torquetcg_rec(torq_mubis,torq_tbis,denedmubis,denedtbis)

        deprec = 0d0

        do kk = 1, npolerecloc
           iipole = polerecglob(kk)
           iglob = ipole(iipole)
           ilocrec = locrec(iglob)
           do betac = 1, 3
              deprec(betac,ilocrec) =  -.5d0*f*denedrbis(betac,kk)
           end do
        end do
        deprec(:,:)=deprec(:,:)-0.5d0*f*(torq_mubis(:,:)+torq_tbis(:,:))
      end if
      return
      end


! TCG2 with guess AND peek-step, NO preconditioner

      subroutine epolar1tcg2quat_2
      use atmlst
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use mpi
      implicit none


      logical :: precond
      integer :: betac, kk, i,
     $           j, iipole, irhs, k,
     $           ierr,iglob,ilocrec,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, term, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, ap10a,
     $          ap11a, ap12a, ap11, ap12, ap13, ap14, ap2a0, ap21a,
     $          ap21, ap22, ap23, ap3a0, ap31, ap32, ap41, spp1, spp2,
     $          omega, k1, k2, k3
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1



      real*8 xx(1)

      real*8, allocatable, dimension(:,:,:,:) ::adte
      real*8, allocatable, dimension(:,:,:) :: denedt,denedtbis,mu_tcg2,
     $            mu_tcg2bis,efibis,efi,oldefibis,
     $            oldr0,oldr0bis,Tr0,Tr0rec,Tr0bis,Tr0recbis,
     $            oldTr0,oldTr0bis,T2r0,T2r0rec,T2r0bis,T2r0recbis,
     $            mu0,mu0bis,r0,r0bis,
     $            T3r0,T3r0bis,T3r0rec,T3r0recbis,
     $            Tefi,Tefibis,Tefirec,Tefirecbis,oldTefi,oldTefibis,
     $            ade,adme,adtebis,oldaefi,oldaefibis,Tarr_dr0,
     $            Tarr_dr0bis,Tarr_dr0rec,Tarr_dr0recbis,
     $            aTarr_dr0,aTarr_dr0bis,oldefi,
     $            ar0,ar0bis,atr0,aTr0bis, at2r0,aT2r0bis, aefi,aefibis,
     $            Taefi,Taefibis,Taefirec,Taefirecbis,T2aefi,T2aefibis,
     $            T2aefirec,T2aefirecbis,mupeek,mupeekbis,
     $            fphiTaE, fphiT2ae,fphir0,
     $            fphiE, fphiaE,fphimu0,
     $            fphiTr0, fphit2r0, fphit3r0, fphiTE, 
     $            fphiar0, fphiatr0, fphiat2r0,  
     $            fphiatarrdr0,tmpbig,tmpbigbis,fphitmp
      real*8, allocatable, dimension(:,:) ::
     $            arr_dTT2r0,arr_dTT2r0bis,
     $            arr_dTr0,arr_dTr0bis,arr_dTTr0,arr_dTTr0bis,
     $            fphiarrdr0,  fphiarr_ded,arr_dr0,arr_dr0bis,
     $            fphiarr_dep, 
     $            fphiarrdtr0, fphiarrdttr0,
     $            fphiarrdtt2r0,tmpfphi2
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu,
     $            torq_tbis,torq_mubis,
     $            denedmu, denedmubis,denedr,denedrbis, cphi, 
     $            arr_ded,arr_dedbis,arr_dep,arr_depbis,adebis,admebis,
     $            adtb,adtbbis,tmpsmall,tmpsmallbis

      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi
      parameter (nrhs=2)

      if (rank.le.ndir-1) then
        allocate(adte(3,3,nrhs, npolebloc))
        allocate(denedt(3,3,npolebloc),mu_tcg2(3,nrhs,npolebloc),
     $ efi(3,nrhs,npolebloc),oldr0(3,nrhs,npolebloc),
     $ Tr0(3,nrhs,npolebloc),Tr0rec(3,nrhs,npoleloc),
     $ oldTr0(3,nrhs,npolebloc),T2r0(3,nrhs,npolebloc),
     $ T2r0rec(3,nrhs,npoleloc),mu0(3,nrhs,npolebloc),
     $ r0(3,nrhs,npolebloc),T3r0(3,nrhs,npolebloc),
     $ T3r0rec(3,nrhs,npoleloc),arr_dTr0(3,npolebloc),
     $ arr_dTTr0(3,npolebloc),arr_dTT2r0(3,npolebloc),
     $ Tefi(3,nrhs,npolebloc),Tefirec(3,nrhs,npoleloc),
     $ oldTefi(3,nrhs,npolebloc),ade(3,nrhs,npolebloc), 
     $ adme(3,nrhs,npolebloc),oldaefi(3,nrhs,npolebloc),
     $ arr_dr0(3,npolebloc),Tarr_dr0(3,nrhs,npolebloc),
     $ Tarr_dr0rec(3,nrhs,npoleloc),atarr_dr0(3,nrhs,npolebloc),
     $ oldefi(3,nrhs,npolebloc),ar0(3,nrhs,npolebloc),
     $ atr0(3,nrhs,npolebloc),at2r0(3,nrhs,npolebloc),
     $ aefi(3,nrhs,npolebloc),Taefi(3,nrhs,npolebloc),
     $ Taefirec(3,nrhs,npoleloc),T2aefi(3,nrhs,npolebloc),
     $ T2aefirec(3,nrhs,npoleloc),mupeek(3,nrhs,npolebloc),
     $ tmpbig(3,nrhs,npolebloc))
        allocate(cphi(10,npoleloc),torq_mu(3,nbloc),
     $ torq_t(3,nbloc), denedmu(3,npolebloc),
     $ denedr(3,npolebloc),arr_ded(3,npolebloc),
     $ arr_dep(3,npolebloc),adtb(3,npolebloc),
     $ tmpsmall(3,npolebloc))
        allocate(buffermpi(10,max(npoleloc,1)))
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        cphi = 0d0
        efi = 0d0
      else
        allocate(adtebis(3,3,npolerecloc),denedtbis(3,3,npolerecloc),
     $ mu_tcg2bis(3,nrhs,npolerecloc),efibis(3,nrhs,npolerecloc),
     $ oldr0bis(3,nrhs,npolerecloc),Tr0bis(3,nrhs,npolerecloc),
     $ Tr0recbis(3,nrhs,npolerecloc),oldTr0bis(3,nrhs,npolerecloc),
     $ T2r0bis(3,nrhs,npolerecloc),T2r0recbis(3,nrhs,npolerecloc),
     $ mu0bis(3,nrhs,npolerecloc),r0bis(3,nrhs,npolerecloc),
     $ T3r0bis(3,nrhs,npolerecloc),T3r0recbis(3,nrhs,npolerecloc),
     $ arr_dTr0bis(3,npolerecloc),arr_dTTr0bis(3,npolerecloc),
     $ arr_dTT2r0bis(3,npolerecloc),Tefibis(3,nrhs,npolerecloc),
     $ Tefirecbis(3,nrhs,npolerecloc),oldTefibis(3,nrhs,npolerecloc),
     $ adebis(3,npolerecloc), admebis(3,npolerecloc),
     $ oldaefibis(3,nrhs,npolerecloc),arr_dr0bis(3,npolerecloc),
     $ Tarr_dr0bis(3,nrhs,npolerecloc),
     $ Tarr_dr0recbis(3,nrhs,npolerecloc),
     $ atarr_dr0bis(3,nrhs,npolerecloc),oldefibis(3,nrhs,npolerecloc),
     $ ar0bis(3,nrhs,npolerecloc),atr0bis(3,nrhs,npolerecloc),
     $ at2r0bis(3,nrhs,npolerecloc),aefibis(3,nrhs,npolerecloc),
     $ Taefibis(3,nrhs,npolerecloc),Taefirecbis(3,nrhs,npolerecloc),
     $ T2aefibis(3,nrhs,npolerecloc),T2aefirecbis(3,nrhs,npolerecloc),
     $ mupeekbis(3,nrhs,npolerecloc),tmpbigbis(3,nrhs,npolerecloc),
     $ fphitmp(20,nrhs,npolerecloc))
        allocate(fphiTaE(20,2,npolerecloc),fphit2ae(20,2,npolerecloc),
     $ fphitr0(20,2,npolerecloc),fphit2r0(20,2,npolerecloc),
     $ fphit3r0(20,2,npolerecloc),fphitE(20,2,npolerecloc),
     $ fphiarrdr0(20,npolerecloc),fphiatarrdr0(20,2,npolerecloc),
     $ fphiarr_ded(20,npolerecloc),fphiarr_dep(20,npolerecloc),
     $ fphir0(20,2,npolerecloc),fphiE(20,2,npolerecloc),
     $ fphiaE(20,2,npolerecloc),fphiarrdtr0(20,npolerecloc),
     $ fphiarrdttr0(20,npolerecloc),fphiarrdtt2r0(20,npolerecloc),
     $ fphiar0(20,2,npolerecloc),fphiaTr0(20,2,npolerecloc),
     $ fphiaT2r0(20,2,npolerecloc),fphimu0(20,2,npolerecloc),
     $ tmpfphi2(20,npolerecloc))
        allocate(torq_tbis(3,nlocrec2),torq_mubis(3,nlocrec2),
     $denedmubis(3,npolerecloc),denedrbis(3,npolerecloc),
     $arr_dedbis(3,npolerecloc),arr_depbis(3,npolerecloc),
     $adtbbis(3,npolerecloc),tmpsmallbis(3,npolerecloc))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
      end if
      

      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      
      f = electric/dielec
      precond = tcgprec

      omega = tcgomega

      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,efi,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs, efi)
        call commfield(nrhs, efi)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            efi(j,1,i) = efi(j,1,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
            efi(j,2,i) = efi(j,2,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
          end do
        end do
        call commdirdir(nrhs,1,efi,reqrec,reqsend)
        call commdirdir(nrhs,2,efi,reqrec,reqsend)
        oldefi = efi
      else
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
      end if


      if (rank.le.ndir-1) then
        call commrecdirdip(nrhs,1,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call diagvec(nrhs, efi, aefi)
        mu0 = aefi
        oldaefi = aefi
        call dbletmatxb_pme(nrhs, .true., aefi, efi, Taefi, Tefi)
        call commfield(nrhs,Tefi)
        call commfield(nrhs,Taefi)
        call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) +
     $    Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)

        call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
        call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
        call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

        call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirsolv(nrhs,0,Tefirecbis,Tefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Tefi(:,:,1:npoleloc) = Tefi(:,:,1:npoleloc) +
     $    Tefirec(:,:,1:npoleloc) - term*efi(:,:,1:npoleloc)
        call commdirdir(nrhs,0,Tefi,reqrec,reqsend)
        call commdirdir(nrhs,1,Tefi,reqrec,reqsend)
        call commdirdir(nrhs,2,Tefi,reqrec,reqsend)
        call commrecdirdip(nrhs,1,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        oldTefi = Tefi
        r0 = efi - Taefi
        oldr0 = r0

        call dbletmatxb_pme(nrhs, .true., Taefi, r0, T2aefi, Tr0)
        call commfield(nrhs,T2aefi)
        call commfield(nrhs,Tr0)
        call commrecdirsolv(nrhs,0,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T2aefi(:,:,1:npoleloc) = T2aefi(:,:,1:npoleloc) + 
     $    T2aefirec(:,:,1:npoleloc) - term*Taefi(:,:,1:npoleloc)
     
        call commdirdir(nrhs,0,T2aefi,reqrec,reqsend)
        call commdirdir(nrhs,1,T2aefi,reqrec,reqsend)
        call commdirdir(nrhs,2,T2aefi,reqrec,reqsend)

        call commrecdirdip(nrhs,1,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) + 
     $     Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

        call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,2,Tr0,reqrec,reqsend)
        call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldTr0 = Tr0

        call tmatxb_pme(nrhs, .true., Tr0, T2r0)
        call commfield(nrhs,T2r0)
        call commrecdirsolv(nrhs,0,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T2r0(:,:,1:npoleloc) = T2r0(:,:,1:npoleloc) +
     $    T2r0rec(:,:,1:npoleloc) - term*Tr0(:,:,1:npoleloc)

        call commdirdir(nrhs,0,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T2r0,reqrec,reqsend)
        call commrecdirdip(nrhs,1,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        call tmatxb_pme(nrhs, .true., T2r0, T3r0)
        call commfield(nrhs,T3r0)

        call commrecdirsolv(nrhs,0,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        T3r0(:,:,1:npoleloc) = T3r0(:,:,1:npoleloc) + 
     $     T3r0rec(:,:,1:npoleloc)- term*T2r0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T3r0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

      call diagvec3(nrhs,r0,Tr0,T2r0,ar0,aTr0,aT2r0)
      else
        call commrecdirdip(nrhs,0,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldefibis = efibis
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          aefibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
        end do
        mu0bis = aefibis
        oldaefibis = aefibis

        call tmatxbrecipsave(aefi,aefibis,nrhs,Taefirec,Taefirecbis,
     $      fphiaE)
        call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        fphimu0 = fphiae
        call tmatxbrecipsave(efi,efibis,nrhs,Tefirec,Tefirecbis,
     $                       fphiE)
        call commrecdirsolv(nrhs,1,Tefirecbis,Tefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tefirecbis,Tefirec,buffermpi,
     $  buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tefibis,Tefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldTefibis = Tefibis
        r0bis = efibis - Taefibis
        oldr0bis = r0bis

        call tmatxbrecipsave(Taefi,Taefibis,nrhs,T2aefirec,T2aefirecbis,
     $                       fphiTaE)
        call commrecdirsolv(nrhs,1,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis, 
     $                          fphir0)
        call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldTr0bis = Tr0bis

        call tmatxbrecipsave(Tr0,Tr0bis,nrhs,T2r0rec,T2r0recbis,
     $                          fphitr0)
        call commrecdirsolv(nrhs,1,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        call tmatxbrecipsave(T2r0,T2r0bis,nrhs,T3r0rec,T3r0recbis,
     $                          fphiT2r0)
        call commrecdirsolv(nrhs,1,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        call fftthatplz2(T2aefi,T2aefibis, fphiT2aE)
        call fftthatplz2(Tefi,Tefibis, fphiTE)
        call fftthatplz2(T3r0,T3r0bis, fphiT3r0)
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          ar0bis(:,:,i) = polarity(iipole)*r0bis(:,:,i)
          aTr0bis(:,:,i) = polarity(iipole)*Tr0bis(:,:,i)
          aT2r0bis(:,:,i) = polarity(iipole)*T2r0bis(:,:,i)
        end do
        call fftthatplz2(ar0,ar0bis, fphiar0)
        call fftthatplz2(aTr0,aTr0bis, fphiaTr0)
        call fftthatplz2(aT2r0,aT2r0bis, fphiaT2r0)
      end if

      if (rank.le.ndir-1) then
        do irhs = 1,nrhs
           n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
           t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        do k = 1,nrhs
           np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2(:,k,:) = mu0(:,k,:) + (g1(k)*t2(k)+ t4(k))*r0(:,k,:) 
     $     - g1(k)*t4(k)*Tr0(:,k,:)
        end do
        sp0 = sprod(3*npoleloc, r0(:,1,:), efi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = sprod(3*npoleloc, Tr0(:,1,:), efi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1
        mupeek = omega*ar0 
     $          - omega*(g1(1)*t2(1) + t4(1))*aTr0 
     $          + omega*g1(1)*t4(1)*at2r0
        mu_tcg2 = mu_tcg2 + mupeek
      else
        n0 = 0d0
        t1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        do k = 1,nrhs
           np1(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2bis(:,k,:) =  mu0bis(:,k,:) + (g1(k)*t2(k) + 
     $       t4(k))*r0bis(:,k,:)- g1(k)*t4(k)*Tr0bis(:,k,:)
        end do
        sp0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1
        mupeekbis = omega*ar0bis
     $          - omega*(g1(1)*t2(1) + t4(1))*aTr0bis
     $          + omega*g1(1)*t4(1)*at2r0bis
        mu_tcg2bis = mu_tcg2bis + mupeekbis
      end if


      !peek  no prec


      ! 2. build derivs arrays...
      ! 2.1 Prepare coefficients
      a10 = t4(1) + g1(1)*t2(1)
      a1m1= -g1(1)*t4(1)
      a11 = 2d0*b1/t1(1) - 2d0*np1(1)*b2/t3(1) 
     $      - 2d0*np1(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1)*t1(1))
     $      + 2d0*t9(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + 2d0*np1(1)*sp0*g1(1)/(t1(1)*t1(1))
      a12 = -2d0*n0(1)*b1/(t1(1)*t1(1)) + 4d0*t1(1)*b2/t3(1) 
     $      - 2d0*n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 4d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 4d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a13 = -4d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0* n0(1)*b2/t3(1)
     $      + 2d0*n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a14 = 2d0*t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a20 = -g1(1)*t4(1)
      a21 = -n0(1)*b1/(t1(1)*t1(1)) + 2d0*t1(1)*b2/t3(1) 
     $      - n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 2d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a22 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a23 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a31 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*g1(1)*sp0/(t1(1)*t1(1))
      a41 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a32 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))


      if (rank.le.ndir-1) then
        spp1 = sprod(3*npoleloc, oldaefi(:,2,:), Tr0(:,1,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        spp2 = sprod(3*npoleloc, oldaefi(:,2,:), T2r0(:,1,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
      else
        spp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        spp2 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
      end if

      k1 = omega*t4(1)*spp2 - omega*t2(1)*spp1
      k2 = omega*g1(1)*spp2 - omega*spp1
      k3 = -omega*g1(1)*spp1

      ap10a = omega
      ap11a = -omega*(t2(1)*g1(1) + t4(1))
      ap12a = omega*t4(1)*g1(1)
      ap2a0 = -omega*(g1(1)*t2(1) + t4(1))
      ap21a = omega*t4(1)*g1(1)
      ap3a0 = omega*g1(1)*t4(1)

      ap11 = k1*(-2d0*np1(1)/t3(1)
     $         - 2d0*np1(1)*np1(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t9(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(2d0/t1(1))
     $      +k3*(2d0*np1(1)/(t1(1)*t1(1)))
      ap12 = k1*(4d0*t1(1)/t3(1) 
     $        - 2d0*n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $        + 4d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $        - 2d0*t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(-2d0*n0(1)/(t1(1)*t1(1)))
     $      +k3*(-4d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap13 = k1*(-4d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1))
     $         - 2d0*n0(1)/t3(1))
     $      +k3*(2d0*n0(1)/(t1(1)*t1(1)))
      ap14 = k1*(2d0*t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap21 = k1*(2d0*t1(1)/t3(1) 
     $         - n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $         - t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(-n0(1)/(t1(1)*t1(1)))
     $      +k3*(-2d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap22 = k1*(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3*(n0(1)/(t1(1)*t1(1)))
      ap23 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap31 = k1*(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3*(n0(1)/(t1(1)*t1(1)))
      ap32 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap41 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))



      if (rank.le.ndir-1) then
        arr_dEp  = mu_tcg2(:,1,:)
        arr_dEd = a10*efi(:,2,:) 
     $              + a1m1*Tefi(:,2,:)
     $              + a11*r0(:,1,:)
     $              + a12*Tr0(:,1,:)
     $              + a13*T2r0(:,1,:) 
     $              + a14*T3r0(:,1,:)
        arr_dTr0 = a20*efi(:,2,:)
     $              + a21*r0(:,1,:)
     $              + a22*Tr0(:,1,:)
     $              + a23*T2r0(:,1,:)
        arr_dTTr0 = a31*r0(:,1,:)
     $              + a32*Tr0(:,1,:)
        arr_dTT2r0 = a41*r0(:,1,:)
        arr_ded = arr_ded 
     $            + ap10a*aefi(:,2,:)
     $            + ap11a*Taefi(:,2,:)
     $            + ap12a*T2aefi(:,2,:)
     $            + ap11*r0(:,1,:)
     $            + ap12*Tr0(:,1,:)
     $            + ap13*T2r0(:,1,:)
     $            + ap14*T3r0(:,1,:)
        arr_dTr0 = arr_dTr0 
     $            + ap2a0*aefi(:,2,:)
     $            + ap21a*Taefi(:,2,:)
     $            + ap21*r0(:,1,:)
     $            + ap22*Tr0(:,1,:)
     $            + ap23*T2r0(:,1,:)
        arr_dTTr0 = arr_dTTr0 
     $            + ap3a0*aefi(:,2,:)
     $            + ap31*r0(:,1,:)
     $            + ap32*Tr0(:,1,:)
        arr_dTT2r0 = arr_dTT2r0 
     $            + ap41*r0(:,1,:)
        ep = sprod(3*npoleloc, mu_tcg2(:,1,:), efi(:,2,:))
        ep = -.5d0*f*ep 
        denedr = 0d0
        denedmu = 0d0
        denedt = 0d0

        atarr_dr0 = 0d0
        arr_dr0 = arr_ded
        tmpbig(:,1,:) = arr_dr0 
        tmpbig(:,2,:) = arr_dr0 

        call tmatxb_pme(nrhs, .true., tmpbig, Tarr_dr0)
        call commfield(nrhs,Tarr_dr0)
        call commrecdirsolv(nrhs,0,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Tarr_dr0(:,:,1:npoleloc) = Tarr_dr0(:,:,1:npoleloc) + 
     $   Tarr_dr0rec(:,:,1:npoleloc) - term*tmpbig(:,:,1:npoleloc)
     
        call commdirdir(nrhs,0,Tarr_dr0,reqrec,reqsend)
        call commdirdir(nrhs,1,Tarr_dr0,reqrec,reqsend)
        call commdirdir(nrhs,2,Tarr_dr0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call diagvec(nrhs, Tarr_dr0, aTarr_dr0)
        arr_ded = arr_dr0 - aTarr_dr0(:,1,:) + mu0(:,2,:)


        call scalderfieldzmat8(
     &                         arr_ded, arr_dep, arr_dtr0,r0(:,1,:),
     &                         arr_dTTr0,Tr0(:,1,:), arr_dTT2r0,
     &                         T2r0(:,1,:), arr_dr0, - oldaefi(:,1,:),
     &                         ade, adme, adte, adtb
     &                        )

        denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb
        denedmu=   adme(:,1,:)   + adme(:,2,:)    
        denedt =  adte(:,:,1,:) + adte(:,:,2,:)

        !N. Compute dEne/dmu and dEne/dtheta (given E=.5*mu.field)

        ! Do the contraction dEne/dmuP * dmuP/dr
        call torquetcg_dir(torq_mu,torq_t,denedmu,denedt)
        dep = 0d0

        do kk = 1, npolebloc
           kkpole = poleglob(kk)
           kglob = ipole(kkpole) 
           kloc = loc(kglob)
           do betac = 1, 3
              dep(betac,kloc) =  -.5d0*denedr(betac,kk)*f
           end do
        end do
        dep(:,:) = dep(:,:)-0.5d0*f*(torq_mu(:,:)+torq_t(:,:))
      else
        arr_dEpbis  = mu_tcg2bis(:,1,:)
        arr_dEdbis = a10*efibis(:,2,:) 
     $              + a1m1*Tefibis(:,2,:)
     $              + a11*r0bis(:,1,:)
     $              + a12*Tr0bis(:,1,:)
     $              + a13*T2r0bis(:,1,:) 
     $              + a14*T3r0bis(:,1,:)
        arr_dTr0bis = a20*efibis(:,2,:)
     $              + a21*r0bis(:,1,:)
     $              + a22*Tr0bis(:,1,:)
     $              + a23*T2r0bis(:,1,:)
        arr_dTTr0bis = a31*r0bis(:,1,:)
     $              + a32*Tr0bis(:,1,:)
        arr_dTT2r0bis = a41*r0bis(:,1,:)
        fphiarr_dep = (t4(1)+g1(1)*t2(1))*fphir0(:,1,:) 
     $   - g1(1)*t4(1)*fphiTr0(:,1,:)
        fphiarr_dEd = a10*fphiE(:,2,:)
     $              + a1m1*fphiTe(:,2,:)
     $              + a11*fphir0(:,1,:)
     $              + a12*fphiTr0(:,1,:)
     $              + a13*fphiT2r0(:,1,:)
     $              + a14*fphiT3r0(:,1,:)
        fphiarrdTr0 = a20*fphie(:,2,:) 
     $              + a21*fphir0(:,1,:)
     $              + a22*fphiTr0(:,1,:)
     $              + a23*fphiT2r0(:,1,:)
        fphiarrdTTr0 = a31*fphir0(:,1,:)
     $               + a32*fphiTr0(:,1,:)
        fphiarrdTT2r0 = a41*fphir0(:,1,:)

        !guess
        fphiarr_dep = fphiarr_dep 
     $          + fphiaE(:,1,:)
        arr_dedbis = arr_dedbis
     $            + ap10a*aefibis(:,2,:)
     $            + ap11a*Taefibis(:,2,:)
     $            + ap12a*T2aefibis(:,2,:)
     $            + ap11*r0bis(:,1,:)
     $            + ap12*Tr0bis(:,1,:)
     $            + ap13*T2r0bis(:,1,:)
     $            + ap14*T3r0bis(:,1,:)
        arr_dTr0bis = arr_dTr0bis
     $            + ap2a0*aefibis(:,2,:)
     $            + ap21a*Taefibis(:,2,:)
     $            + ap21*r0bis(:,1,:)
     $            + ap22*Tr0bis(:,1,:)
     $            + ap23*T2r0bis(:,1,:)
        arr_dTTr0bis = arr_dTTr0bis
     $            + ap3a0*aefibis(:,2,:)
     $            + ap31*r0bis(:,1,:)
     $            + ap32*Tr0bis(:,1,:)
        arr_dTT2r0bis = arr_dTT2r0bis
     $            + ap41*r0bis(:,1,:)
        fphiarr_dep = fphiarr_dep 
     $            + ap10a*fphiar0(:,1,:)
     $            + ap11a*fphiaTr0(:,1,:)
     $            + ap21a*fphiaT2r0(:,1,:)
        fphiarr_ded = fphiarr_ded 
     $            + ap10a*fphiae(:,2,:)
     $            + ap11a*fphiTae(:,2,:)
     $            + ap12a*fphiT2ae(:,2,:)
     $            + ap11*fphir0(:,1,:)
     $            + ap12*fphiTr0(:,1,:)
     $            + ap13*fphiT2r0(:,1,:)
     $            + ap14*fphiT3r0(:,1,:)
        fphiarrdTr0 = fphiarrdTr0 
     $            + ap2a0*fphiae(:,2,:)
     $            + ap21a*fphiTae(:,2,:)
     $            + ap21 *fphir0(:,1,:)
     $            + ap22 *fphiTr0(:,1,:)
     $            + ap23 *fphiT2r0(:,1,:)
        fphiarrdTTr0 = fphiarrdTTr0 
     $            + ap3a0*fphiae(:,2,:)
     $            + ap31 *fphir0(:,1,:)
     $            + ap32 *fphiTr0(:,1,:)
        fphiarrdTT2r0 = fphiarrdTT2r0 
     $            + ap41*fphir0(:,1,:)
        ep = 0d0

        arr_dr0bis = arr_dedbis
        tmpbigbis(:,1,:) = arr_dr0bis 
        tmpbigbis(:,2,:) = arr_dr0bis 
        fphiarrdr0 = fphiarr_ded
        call tmatxbrecipsave(tmpbig,tmpbigbis,nrhs,Tarr_dr0rec,
     $        Tarr_dr0recbis,fphitmp)
        call commrecdirsolv(nrhs,1,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          aTarr_dr0bis(:,:,i) = polarity(iipole)*Tarr_dr0bis(:,:,i)
        end do
        call fftthatplz2(aTarr_dr0,aTarr_dr0bis, fphiaTarrdr0)
        arr_dedbis = arr_dr0bis - aTarr_dr0bis(:,1,:) + 
     $   mu0bis(:,2,:)
        fphiarr_ded = fphiarr_ded - fphiatarrdr0(:,1,:) + 
     $   fphimu0(:,2,:)

        call scalderfieldzmatrec8( 
     &                            xx,arr_depbis,fphiarr_dep,xx,
     &                            arr_dedbis, fphiarr_ded, xx,
     &                            arr_dtr0bis,fphiarrdtr0, xx,
     &                            r0bis(:,1,:), fphir0(:,1,:),xx,
     &                            arr_dttr0bis,fphiarrdttr0,xx,
     &                            tr0bis(:,1,:), fphitr0(:,1,:),xx,
     &                            arr_dtt2r0bis,fphiarrdtt2r0,xx,
     &                            T2r0bis(:,1,:), fphit2r0(:,1,:), xx,
     &                          - arr_dr0bis, - fphiarrdr0,xx,
     &                            oldaefibis(:,1,:), fphiaE(:,1,:),
     &                            adebis, admebis, adtebis, adtbbis
     &                           )

        denedrbis = adebis + adtbbis
        denedmubis = admebis
        denedtbis = adtebis
        call torquetcg_rec(torq_mubis,torq_tbis,denedmubis,denedtbis)

        deprec = 0d0

        do kk = 1, npolerecloc
           iipole = polerecglob(kk)
           iglob = ipole(iipole)
           ilocrec = locrec(iglob)
           do betac = 1, 3
              deprec(betac,ilocrec) =  -.5d0*f*denedrbis(betac,kk)
           end do
        end do
        deprec(:,:)=deprec(:,:)-0.5d0*f*(torq_mubis(:,:)+torq_tbis(:,:))
      end if

      return
      end

! TCG2 with fastgrad, fully optimized, all options 
! aboard
      subroutine epolar1tcg2cinq_2
      use atmlst
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use mpi
      implicit none
      logical :: precond, isguess, peek
      integer :: betac, kk, i,
     $           j, iipole, irhs, k,
     $           ierr,iglob,ilocrec,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, term, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, ap10a,
     $          ap11a, ap12a, ap11, ap12, ap13, ap14, ap2a0, ap21a,
     $          ap21, ap22, ap23, ap3a0, ap31, ap32, ap41, spp1, spp2,
     $          omega, k1, k2, k3
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1

      
      
      real*8 xx(1)

      real*8, allocatable, dimension(:,:,:,:) ::adte
      real*8, allocatable, dimension(:,:,:) :: denedt,denedtbis,mu_tcg2,
     $            mu_tcg2bis,efi,efibis,oldt2r0,oldt2r0bis,
     $            oldr0,oldr0bis,Tr0,Tr0bis,Tr0rec,Tr0recbis,
     $            oldTr0,oldTr0bis,T2r0,T2r0bis,T2r0rec,T2r0recbis,
     $            mu0,mu0bis,r0,r0bis,
     $            T3r0,T3r0bis,T3r0rec,T3r0recbis, 
     $            T2aefi, t2aefibis, t2aefirec, t2aefirecbis,
     $            Tefi,Tefibis,
     $            ade,adme,adtebis,oldaefi,oldaefibis,tarr_dr0,
     $            tarr_dr0bis,tarr_dr0rec,tarr_dr0recbis,
     $            atarr_dr0,atarr_dr0bis,oldefi,
     $            oldefibis,aefi,
     $            aefibis,Taefi,Taefibis,Taefirec,Taefirecbis,


     $            mupeek,mupeekbis,T2efi,T2efibis,T2efirec,
     $            T2efirecbis,fphiTaE,
     $            fphiTr0, fphit2r0, fphit3r0, fphiTE, fphiT2E,
     $            fphir0, 
     $            fphiE, fphiaE,
     $            fphiatarrdr0, fphimu0,
     $            tmpfphi2,tmpbig,tmpbigbis,fphitmp,
     $            fphioldefi,fphioldr0,fphioldtr0,fphioldt2r0
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu,
     $            torq_tbis,torq_mubis,arr_dr0,arr_dr0bis,
     $            arr_dTr0,arr_dTr0bis,arr_dTTr0,arr_dTTr0bis,
     $            arr_dTT2r0,arr_dTT2r0bis,fphiarrdr0,fphiarr_ded, 
     $            fphiarr_dep,fphiarrdtr0, fphiarrdttr0,fphiarrdtt2r0,
     $            denedmu,denedmubis,denedr,denedrbis,cphi,
     $            arr_ded,arr_dedbis,arr_dep,arr_depbis,adebis,admebis,
     $            adtb,adtbbis,tmpsmall,tmpsmallbis

      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi
      parameter (nrhs=2)

      if (rank.le.ndir-1) then
        allocate(adte(3,3,nrhs,npolebloc),denedt(3,3,npoleloc),
     $ efi(3,nrhs,npolebloc),oldr0(3,nrhs,npolebloc),
     $ oldt2r0(3,nrhs,npolebloc),Tr0(3,nrhs,npolebloc),
     $ Tr0rec(3,nrhs,npoleloc),oldTr0(3,nrhs,npolebloc),
     $ T3r0rec(3,nrhs,npoleloc),T2r0(3,nrhs,npolebloc),
     $ T2r0rec(3,nrhs,npoleloc),T2aefi(3,nrhs,npolebloc),
     $ T2aefirec(3,nrhs,npoleloc),mu0(3,nrhs,npolebloc),
     $ r0(3,nrhs,npolebloc),T3r0(3,nrhs,npolebloc),
     $ arr_dTr0(3,npolebloc),arr_dTTr0(3,npoleloc),
     $ arr_dTT2r0(3,npolebloc),Tefi(3,nrhs,npolebloc),
     $ ade(3,nrhs,npolebloc),adme(3,nrhs,npolebloc),
     $ oldaefi(3,nrhs,npolebloc),arr_dr0(3,npolebloc),
     $ Tarr_dr0(3,nrhs,npolebloc),Tarr_dr0rec(3,nrhs,npoleloc),
     $ atarr_dr0(3,nrhs,npolebloc),oldefi(3,nrhs,npolebloc),
     $ aefi(3,nrhs,npolebloc),Taefi(3,nrhs,npolebloc),
     $ Taefirec(3,nrhs,npoleloc),mupeek(3,nrhs,npolebloc),
     $ T2efi(3,nrhs,npolebloc),T2efirec(3,nrhs,npoleloc),
     $ tmpbig(3,nrhs,npolebloc),cphi(10, npoleloc),
     $ torq_mu(3,nbloc),torq_t(3,nbloc),denedmu(3,npolebloc),
     $ denedr(3,npolebloc),arr_ded(3,npolebloc),
     $ arr_dep(3,npolebloc),adtb(3,npolebloc),
     $ tmpsmall(3,npolebloc),mu_tcg2(3,nrhs,npolebloc))
        allocate (buffermpi(10,max(npoleloc,1)))
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        cphi = 0d0
        efi = 0d0
        Tr0 = 0d0
        Tarr_dr0 = 0d0
      else
        allocate(adtebis(3,3,npolerecloc),
     $ mu_tcg2bis(3,nrhs,npolerecloc),efibis(3,nrhs,npolerecloc),
     $ oldr0bis(3,nrhs,npolerecloc),oldt2r0bis(3,nrhs,npolerecloc),
     $ Tr0bis(3,nrhs,npolerecloc),T3r0bis(3,nrhs,npolerecloc),
     $ Tr0recbis(3,nrhs,npolerecloc),T3r0recbis(3,nrhs,npolerecloc),
     $ T2r0bis(3,nrhs,npolerecloc),T2r0recbis(3,nrhs,npolerecloc),
     $ T2aefibis(3,nrhs,npolerecloc),T2aefirecbis(3,nrhs,npolerecloc),
     $ mu0bis(3,nrhs,npolerecloc),r0bis(3,nrhs,npolerecloc),
     $ arr_dTr0bis(3,npolerecloc),arr_dTTr0bis(3,npolerecloc),
     $ arr_dTT2r0bis(3,npolerecloc),Tefibis(3,nrhs,npolerecloc),
     $ adebis(3,npolerecloc), admebis(3,npolerecloc),
     $ oldaefibis(3,nrhs,npolerecloc),arr_dr0bis(3,npolerecloc),
     $ Tarr_dr0bis(3,nrhs,npolerecloc),Taefibis(3,nrhs,npolerecloc),
     $ Tarr_dr0recbis(3,nrhs,npolerecloc),mupeekbis(3,nrhs,npolerecloc),
     $ atarr_dr0bis(3,nrhs,npolerecloc),Taefirecbis(3,nrhs,npolerecloc),
     $ oldefibis(3,nrhs,npolerecloc),aefibis(3,nrhs,npolerecloc),
     $ T2efibis(3,nrhs,npolerecloc),T2efirecbis(3,nrhs,npolerecloc),
     $ tmpbigbis(3,nrhs,npolerecloc),fphitmp(20,nrhs,npolerecloc))
        allocate(fphiTaE(20,2,npolerecloc),fphitr0(20,2,npolerecloc),
     $ fphit2r0(20,2,npolerecloc),fphit3r0(20,2,npolerecloc),
     $ fphitE(20,2,npolerecloc),fphit2E(20,2,npolerecloc),
     $ fphiarrdr0(20,npolerecloc),fphiatarrdr0(20,nrhs,npolerecloc),
     $ fphiarr_ded(20,npolerecloc),fphiarr_dep(20,npolerecloc),
     $ fphir0(20,2,npolerecloc),fphiE(20,2,npolerecloc),
     $ fphiaE(20,2,npolerecloc),fphiarrdtr0(20,npolerecloc),
     $ fphiarrdttr0(20,npolerecloc),fphiarrdtt2r0(20,npolerecloc),
     $ fphimu0(20,2,npolerecloc),tmpfphi2(20,2,npolerecloc),
     $ fphioldefi(20,2,npolerecloc),fphioldr0(20,2,npolerecloc),
     $ fphioldtr0(20,2,npolerecloc),fphioldt2r0(20,2,npolerecloc))
        allocate(torq_tbis(3,nlocrec2),torq_mubis(3,nlocrec2),
     $ denedmubis(3,npolerecloc),denedrbis(3,npolerecloc),
     $ arr_dedbis(3,npolerecloc),arr_depbis(3,npolerecloc),
     $ adtbbis(3,npolerecloc),tmpsmallbis(3,npolerecloc))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
      end if
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      
      f = electric/dielec

      precond = tcgprec
      isguess = tcgguess
      peek    = tcgpeek

      !1. prepare arrays...

      omega = tcgomega

      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,efi,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs, efi)
        call commfield(nrhs, efi)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            efi(j,1,i) = efi(j,1,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
            efi(j,2,i) = efi(j,2,i)- cphi(j+1,i)+ term*rpole(j+1,iipole)
          end do
        end do
        call commdirdir(nrhs,1,efi,reqrec,reqsend)
        call commdirdir(nrhs,2,efi,reqrec,reqsend)
        oldefi = efi
      else
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $ reqrecdirrec,reqrecdirsend)
      end if

      if (rank.le.ndir-1) then
        call commrecdirdip(nrhs,1,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call diagvec(nrhs, efi, efi)
        aefi = efi
        oldaefi = efi
        mu0 = efi
        call tmatxb_pme(nrhs, .true., aefi, Taefi)
        call commfield(nrhs,Taefi)
        call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) +
     $    Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)
        call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
        call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
        call commdirdir(nrhs,2,Taefi,reqrec,reqsend)
        call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        Tefi = Taefi
        call diagvec(nrhs, Tefi, Tefi)
        oldr0 = oldefi - Taefi
        call diagvec(nrhs, oldr0, r0)
        call dbletmatxb_pme(nrhs, .true.,Taefi,r0, T2aefi, Tr0)
        call commfield(nrhs,T2aefi)
        call commfield(nrhs,Tr0)
        call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) + 
     $     Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
        call commdirdir(nrhs,2,Tr0,reqrec,reqsend)
        call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirsolv(nrhs,0,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T2aefi(:,:,1:npoleloc) = T2aefi(:,:,1:npoleloc) 
     $         + T2aefirec(:,:,1:npoleloc) - term*Taefi(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T2aefi,reqrec,reqsend)
        call commdirdir(nrhs,1,T2aefi,reqrec,reqsend)
        call commdirdir(nrhs,2,T2aefi,reqrec,reqsend)

        call commrecdirdip(nrhs,1,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldTr0 = Tr0
        call diagvec(nrhs, Tr0, Tr0)
        call dbletmatxb_pme(nrhs, .true., Tr0, Tefi, T2r0, T2efi)
        call commfield(nrhs,T2r0)
        call commfield(nrhs,T2efi)
        call commrecdirsolv(nrhs,0,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T2efi(:,:,1:npoleloc) = T2efi(:,:,1:npoleloc) + 
     $    T2efirec(:,:,1:npoleloc) - term*Tefi(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T2efi,reqrec,reqsend)
        call commdirdir(nrhs,1,T2efi,reqrec,reqsend)
        call commdirdir(nrhs,2,T2efi,reqrec,reqsend)
        call commrecdirdip(nrhs,1,T2efibis,T2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2efibis,T2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call diagvec(nrhs, T2efi, T2efi)
        call commrecdirsolv(nrhs,0,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        T2r0(:,:,1:npoleloc) = T2r0(:,:,1:npoleloc) +
     $    T2r0rec(:,:,1:npoleloc) - term*Tr0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T2r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T2r0,reqrec,reqsend)
        call commrecdirdip(nrhs,1,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call diagvec(nrhs, T2r0, T2r0)

        call tmatxb_pme(nrhs, .true., T2r0, T3r0)
        call commfield(nrhs,T3r0)

        call commrecdirsolv(nrhs,0,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        T3r0(:,:,1:npoleloc) = T3r0(:,:,1:npoleloc) + 
     $     T3r0rec(:,:,1:npoleloc)- term*T2r0(:,:,1:npoleloc)
        call commdirdir(nrhs,0,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,1,T3r0,reqrec,reqsend)
        call commdirdir(nrhs,2,T3r0,reqrec,reqsend)

        call commrecdirdip(nrhs,1,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call diagvec(nrhs, T3r0, T3r0)
      else
        call commrecdirdip(nrhs,0,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldefibis = efibis
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          efibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
        end do
        aefibis = efibis
        oldaefibis = efibis
        mu0bis = efibis
        call tmatxbrecipsave(aefi,aefibis,nrhs,Taefirec,Taefirecbis,
     $    fphiaE)
        call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        fphiE = fphiae
        fphimu0 = fphiaE
        Tefibis = Taefibis
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          Tefibis(:,:,i) = polarity(iipole)*Tefibis(:,:,i)
        end do
        oldr0bis = oldefibis - Taefibis
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          r0bis(:,:,i) = polarity(iipole)*oldr0bis(:,:,i)
        end do
        call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis, 
     $                          fphir0)
        call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call tmatxbrecipsave(Taefi,oldefibis,nrhs,T2aefirec,
     $    T2aefirecbis,fphitae)
        call commrecdirsolv(nrhs,1,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2aefirecbis,T2aefirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2aefibis,T2aefi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldTr0bis = Tr0bis
        fphiTE = fphiTaE
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          Tr0bis(:,:,i) = polarity(iipole)*Tr0bis(:,:,i)
        end do
        call tmatxbrecipsave(Tefi,Tefibis,nrhs,T2efirec,T2efirecbis,
     $                       fphiTe)
        call commrecdirsolv(nrhs,1,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2efirecbis,T2efirec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T2efibis,T2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2efibis,T2efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          T2efibis(:,:,i) = polarity(iipole)*T2efibis(:,:,i)
        end do
        call tmatxbrecipsave(Tr0,Tr0bis,nrhs,T2r0rec,T2r0recbis,
     $                       fphitr0)
        call commrecdirsolv(nrhs,1,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T2r0recbis,T2r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T2r0bis,T2r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          T2r0bis(:,:,i) = polarity(iipole)*T2r0bis(:,:,i)
        end do
        call tmatxbrecipsave(T2r0,T2r0bis,nrhs,T3r0rec,T3r0recbis,
     $                    fphiT2r0)
        call commrecdirsolv(nrhs,1,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,T3r0recbis,T3r0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,T3r0bis,T3r0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          T3r0bis(:,:,i) = polarity(iipole)*T3r0bis(:,:,i)
        end do
        call fftthatplz2(T3r0,T3r0bis, fphit3r0)
        call fftthatplz2(T2efi,T2efibis, fphit2e)
      end if




      if (rank.le.ndir-1) then
        do irhs = 1,nrhs
           n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
           t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        do k = 1,nrhs
           np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $        COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2(:,k,:) = mu0(:,k,:)+ (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $     - g1(k)*t4(k)*Tr0(:,k,:)
        end do
        sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1

        mupeek = omega*r0 
     $           - omega*(g1(1)*t2(1) + t4(1))*Tr0 
     $           + omega*g1(1)*t4(1)*t2r0
        mu_tcg2 = mu_tcg2 + mupeek
      else
        do irhs = 1,nrhs
           n0(irhs) = 0d0 
           t1(irhs) = 0d0 
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        do k = 1,nrhs
           np1(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
           t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
           t9(k) = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $        COMM_TINKER,ierr)
           t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
           t3(k) = t1(k)*t8(k)
           t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
           g1(k)  = t10(k)/t3(k)
           mu_tcg2bis(:,k,:) = mu0bis(:,k,:) + (g1(k)*t2(k) + t4(k))*
     $     r0bis(:,k,:)- g1(k)*t4(k)*Tr0bis(:,k,:)
        end do
        sp0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        b1 = sp0 - g1(1)*sp1
        b2 = sp0*t2(1) - t4(1)*sp1
        mupeekbis = omega*r0bis
     $           - omega*(g1(1)*t2(1) + t4(1))*Tr0bis
     $           + omega*g1(1)*t4(1)*t2r0bis
        mu_tcg2bis = mu_tcg2bis + mupeekbis
      end if





      a10 = t4(1) + g1(1)*t2(1)
      a1m1= -g1(1)*t4(1)
      a11 = 2d0*b1/t1(1) - 2d0*np1(1)*b2/t3(1) 
     $      - 2d0*np1(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1)*t1(1))
     $      + 2d0*t9(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + 2d0*np1(1)*sp0*g1(1)/(t1(1)*t1(1))
      a12 = -2d0*n0(1)*b1/(t1(1)*t1(1)) + 4d0*t1(1)*b2/t3(1) 
     $      - 2d0*n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 4d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 4d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a13 = -4d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0* n0(1)*b2/t3(1)
     $      + 2d0*n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a14 = 2d0*t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a20 = -g1(1)*t4(1)
      a21 = -n0(1)*b1/(t1(1)*t1(1)) + 2d0*t1(1)*b2/t3(1) 
     $      - n0(1)*t9(1)*t10(1)*b2/(t1(1)*t3(1)*t3(1))
     $      + 2d0*t2(1)*np1(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - t8(1)*t10(1)*b2/(t3(1)*t3(1))
     $      - 2d0*n0(1)*np1(1)*sp0*g1(1)/(t1(1)*t1(1)*t1(1))
      a22 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*sp0*g1(1)/(t1(1)*t1(1))
      a23 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a31 = -n0(1)*b2/t3(1) - 2d0*t1(1)*t2(1)*t10(1)*b2/(t3(1)*t3(1))
     $      + n0(1)*g1(1)*sp0/(t1(1)*t1(1))
      a41 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))
      a32 = t1(1)*t4(1)*t10(1)*b2/(t3(1)*t3(1))


      if (rank.le.ndir-1) then
        spp1 = sprod(3*npoleloc, oldefi(:,2,:), Tr0(:,1,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        spp2 = sprod(3*npoleloc, oldefi(:,2,:), T2r0(:,1,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
      else
        spp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        spp2 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
      end if
     

      k1 = omega*t4(1)*spp2 - omega*t2(1)*spp1
      k2 = omega*g1(1)*spp2 - omega*spp1
      k3 = -omega*g1(1)*spp1

      ap10a = omega
      ap11a = -omega*(t2(1)*g1(1) + t4(1))
      ap12a = omega*t4(1)*g1(1)
      ap2a0 = -omega*(g1(1)*t2(1) + t4(1))
      ap21a = omega*t4(1)*g1(1)
      ap3a0 = omega*g1(1)*t4(1)

      ap11 = k1*(-2d0*np1(1)/t3(1)
     $         - 2d0*np1(1)*np1(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t9(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(2d0/t1(1))
     $      +k3*(2d0*np1(1)/(t1(1)*t1(1)))
      ap12 = k1*(4d0*t1(1)/t3(1) 
     $        - 2d0*n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $        + 4d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $        - 2d0*t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(-2d0*n0(1)/(t1(1)*t1(1)))
     $      +k3*(-4d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap13 = k1*(-4d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1))
     $         - 2d0*n0(1)/t3(1))
     $      +k3*(2d0*n0(1)/(t1(1)*t1(1)))
      ap14 = k1*(2d0*t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap21 = k1*(2d0*t1(1)/t3(1) 
     $         - n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $         - t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2*(-n0(1)/(t1(1)*t1(1)))
     $      +k3*(-2d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap22 = k1*(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3*(n0(1)/(t1(1)*t1(1)))
      ap23 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap31 = k1*(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3*(n0(1)/(t1(1)*t1(1)))
      ap32 = k1*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap41 = k*(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))

      ! 2.2 Assemble arrays
      
      if (rank.le.ndir-1) then
        arr_dEp  = mu_tcg2(:,1,:)
        arr_dEd = a10*efi(:,2,:) 
     $              + a1m1*Tefi(:,2,:)
     $              + a11*r0(:,1,:)
     $              + a12*Tr0(:,1,:)
     $              + a13*T2r0(:,1,:) 
     $              + a14*T3r0(:,1,:)
        arr_dTr0 = a20*efi(:,2,:) 
     $              + a21*r0(:,1,:)
     $              + a22*Tr0(:,1,:)
     $              + a23*T2r0(:,1,:)
        arr_dTTr0 = a31*r0(:,1,:)
     $              + a32*Tr0(:,1,:)
        arr_dTT2r0 = a41*r0(:,1,:)
        arr_ded = arr_ded 
     $            + ap10a*efi(:,2,:)
     $            + ap11a*Tefi(:,2,:)
     $            + ap12a*T2efi(:,2,:)
     $            + ap11*r0(:,1,:)
     $            + ap12*Tr0(:,1,:)
     $            + ap13*T2r0(:,1,:)
     $            + ap14*T3r0(:,1,:)
        arr_dTr0 = arr_dTr0 
     $            + ap2a0*efi(:,2,:)
     $            + ap21a*Tefi(:,2,:)
     $            + ap21*r0(:,1,:)
     $            + ap22*Tr0(:,1,:)
     $            + ap23*T2r0(:,1,:)
        arr_dTTr0 = arr_dTTr0 
     $            + ap3a0*efi(:,2,:)
     $            + ap31*r0(:,1,:)
     $            + ap32*Tr0(:,1,:)
        arr_dTT2r0 = arr_dTT2r0 
     $            + ap41*r0(:,1,:)
        ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
        ep = -.5d0*f*ep 
        denedr = 0d0
        denedmu = 0d0
        denedt = 0d0

        atarr_dr0 = 0d0
        arr_dr0 = arr_ded
        tmpbig(:,1,:) = arr_dr0
        tmpbig(:,2,:) = arr_dr0
        call tmatxb_pme(nrhs, .true., tmpbig, Tarr_dr0)
        call commfield(nrhs,Tarr_dr0)

        call commrecdirsolv(nrhs,0,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)

        Tarr_dr0(:,:,1:npoleloc) = Tarr_dr0(:,:,1:npoleloc) + 
     $   Tarr_dr0rec(:,:,1:npoleloc) - term*tmpbig(:,:,1:npoleloc)

        call commdirdir(nrhs,0,Tarr_dr0,reqrec,reqsend)
        call commdirdir(nrhs,1,Tarr_dr0,reqrec,reqsend)
        call commdirdir(nrhs,2,Tarr_dr0,reqrec,reqsend)
        call commrecdirdip(nrhs,1,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)

        call diagvec(nrhs, tarr_dr0, atarr_dr0)

        arr_ded = arr_dr0 - atarr_dr0(:,1,:) + mu0(:,2,:)

        call scalderfieldzmat8(
     &                         arr_ded, arr_dep, arr_dtr0, r0(:,1,:),
     &                         arr_dTTr0,  Tr0(:,1,:), arr_dTT2r0,
     &                         T2r0(:,1,:), arr_dr0,  - oldaefi(:,1,:),
     &                         ade, adme, adte, adtb
     &                        )

        denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb
        denedmu=   adme(:,1,:)   + adme(:,2,:)    
        denedt =  adte(:,:,1,:) + adte(:,:,2,:)

        !N. Compute dEne/dmu and dEne/dtheta (given E=.5*mu.field)
        call torquetcg_dir(torq_mu,torq_t,denedmu,denedt)

        dep = 0d0

        do kk = 1, npolebloc
           kkpole = poleglob(kk)
           kglob = ipole(kkpole) 
           kloc = loc(kglob)
           do betac = 1, 3
              dep(betac,kloc) =  -.5d0*denedr(betac,kk)*f
           end do
        end do
        dep(:,:) = dep(:,:)-0.5d0*f*(torq_mu(:,:)+torq_t(:,:))
      else
        arr_dEpbis  = mu_tcg2bis(:,1,:)
        arr_dEdbis = a10*efibis(:,2,:) 
     $              + a1m1*Tefibis(:,2,:)
     $              + a11*r0bis(:,1,:)
     $              + a12*Tr0bis(:,1,:)
     $              + a13*T2r0bis(:,1,:) 
     $              + a14*T3r0bis(:,1,:)
        arr_dTr0bis = a20*efibis(:,2,:) 
     $              + a21*r0bis(:,2,:)
     $              + a22*Tr0bis(:,1,:)
     $              + a23*T2r0bis(:,1,:)
        arr_dTTr0bis = a31*r0bis(:,1,:)
     $              + a32*Tr0bis(:,1,:)
        arr_dTT2r0bis = a41*r0bis(:,1,:)
        fphiarr_dep = (t4(1)+g1(1)*t2(1))*fphir0(:,1,:) 
     $   - g1(1)*t4(1)*fphiTr0(:,1,:)
        fphiarr_dEd = a10*fphiE(:,2,:)
     $              + a1m1*fphiTe(:,2,:)
     $              + a11*fphir0(:,1,:)
     $              + a12*fphiTr0(:,1,:)
     $              + a13*fphiT2r0(:,1,:)
     $              + a14*fphiT3r0(:,1,:)
        fphiarrdTr0 = a20*fphie(:,2,:) 
     $              + a21*fphir0(:,1,:)
     $              + a22*fphiTr0(:,1,:)
     $              + a23*fphiT2r0(:,1,:)
        fphiarrdTTr0 = a31*fphir0(:,1,:)
     $               + a32*fphiTr0(:,1,:)
        fphiarrdTT2r0 = a41*fphir0(:,1,:)

        ! guess adds...
        fphiarr_dep = fphiarr_dep 
     $              + fphiaE(:,1,:)
        arr_dedbis = arr_dedbis
     $            + ap10a*efibis(:,2,:)
     $            + ap11a*Tefibis(:,2,:)
     $            + ap12a*T2efibis(:,2,:)
     $            + ap11*r0bis(:,1,:)
     $            + ap12*Tr0bis(:,1,:)
     $            + ap13*T2r0bis(:,1,:)
     $            + ap14*T3r0bis(:,1,:)
        arr_dTr0bis = arr_dTr0bis 
     $            + ap2a0*efibis(:,2,:)
     $            + ap21a*Tefibis(:,2,:)
     $            + ap21*r0bis(:,1,:)
     $            + ap22*Tr0bis(:,1,:)
     $            + ap23*T2r0bis(:,1,:)
        arr_dTTr0bis = arr_dTTr0bis
     $            + ap3a0*efibis(:,2,:)
     $            + ap31*r0bis(:,1,:)
     $            + ap32*Tr0bis(:,1,:)
        arr_dTT2r0bis = arr_dTT2r0bis
     $            + ap41*r0bis(:,1,:)
        fphiarr_dep = fphiarr_dep
     $            + ap10a*fphir0(:,1,:)
     $            + ap11a*fphiTr0(:,1,:)
     $            + ap21a*fphiT2r0(:,1,:)
        fphiarr_ded = fphiarr_ded 
     $            + ap10a*fphie(:,2,:)
     $            + ap11a*fphiTe(:,2,:)
     $            + ap12a*fphiT2e(:,2,:)
     $            + ap11 *fphir0(:,1,:)
     $            + ap12 *fphiTr0(:,1,:)
     $            + ap13 *fphiT2r0(:,1,:)
     $            + ap14 *fphiT3r0(:,1,:)
        fphiarrdTr0 = fphiarrdTr0 
     $            + ap2a0*fphie(:,2,:)
     $            + ap21a*fphiTe(:,2,:)
     $            + ap21 *fphir0(:,1,:)
     $            + ap22 *fphiTr0(:,1,:)
     $            + ap23 *fphiT2r0(:,1,:)
        fphiarrdTTr0 = fphiarrdTTr0 
     $            + ap3a0*fphie(:,2,:)
     $            + ap31 *fphir0(:,1,:)
     $            + ap32 *fphiTr0(:,1,:)
        fphiarrdTT2r0 = fphiarrdTT2r0 
     $            + ap41*fphir0(:,1,:)
        ep = 0d0
        arr_dr0bis = arr_dedbis
        tmpbigbis(:,1,:) = arr_dr0bis
        tmpbigbis(:,2,:) = arr_dr0bis
        fphiarrdr0 = fphiarr_ded
        call tmatxbrecipsave(tmpbig,tmpbigbis,nrhs,Tarr_dr0rec,
     $                    Tarr_dr0recbis,fphitmp)
        call commrecdirsolv(nrhs,1,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirsolv(nrhs,2,Tarr_dr0recbis,Tarr_dr0rec,buffermpi,
     $    buffermpi,reqrecdirrec,reqrecdirsend)
        call commrecdirdip(nrhs,0,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,Tarr_dr0bis,Tarr_dr0,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          aTarr_dr0bis(:,:,i) = polarity(iipole)*Tarr_dr0bis(:,:,i)
        end do
        call fftthatplz2(atarr_dr0,atarr_dr0bis, fphiatarrdr0)
        arr_dedbis = arr_dr0bis - atarr_dr0bis(:,1,:)+mu0bis(:,2,:)
        fphiarr_ded = fphiarr_ded - fphiatarrdr0(:,1,:) + fphimu0(:,2,:)

        call scalderfieldzmatrec8(
     &                            xx,arr_depbis,fphiarr_dep,xx,
     &                            arr_dedbis,fphiarr_ded,xx,arr_dtr0bis,
     &                            fphiarrdtr0, xx,r0bis(:,1,:),
     &                            fphir0(:,1,:),xx, arr_dttr0bis,
     &                            fphiarrdttr0,xx, tr0bis(:,1,:),
     &                            fphitr0(:,1,:),xx, arr_dtt2r0bis,
     &                            fphiarrdtt2r0,xx, t2r0bis(:,1,:),
     &                            fphit2r0(:,1,:),xx, -arr_dr0bis,
     &                          - fphiarrdr0, xx, oldaefibis(:,1,:),
     &                            fphiaE(:,1,:), adebis, admebis,
     &                            adtebis, adtbbis
     &                           )
        denedrbis = adebis + adtbbis
        denedmubis =  admebis
        denedtbis = adtebis
        call torquetcg_rec(torq_mubis,torq_tbis,denedmubis,denedtbis)

        do kk = 1, npolerecloc
           iipole = polerecglob(kk)
           iglob = ipole(iipole)
           ilocrec = locrec(iglob)
           do betac = 1, 3
              deprec(betac,ilocrec) =  -.5d0*f*denedrbis(betac,kk)
           end do
        end do
        deprec(:,:)=deprec(:,:)-0.5d0*f*(torq_mubis(:,:)+torq_tbis(:,:))
       
      end if

      return
      end

