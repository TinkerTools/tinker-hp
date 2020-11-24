c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1tcg  --  induced dipole energy & gradient  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1tcg" calculates the induced dipole polarization energy and gradients with the tcg method
c
c
      subroutine epolar1tcg
      use domdec
      use iounit
      use polpot
      use potent
      implicit none
 1000 format(' illegal tcg order')
c
c     recompute the optimal value for omega (TCG-PEEK) if needed
c
      if (omegafitstep) call omegafit
c
c     choose the method for summing over polarization interactions
c
      if (use_polarshortreal) then
        if (tcgordershort.eq.1) then
          call epolar1tcg1shortreal
        else if (tcgordershort.eq.2) then
          if ((tcgprecshort).and.(tcgpeekshort).and.(tcgguessshort)) 
     $       then
            call epolar1tcg2cinqshortreal
          else if ((tcgguessshort.eqv..false.).and.
     $      (tcgpeekshort.eqv..false.)) then 
            call epolar1tcg2shortreal
          else if ((tcgguessshort).and.(tcgpeekshort.eqv..false.)) then
            call epolar1tcg2bisshortreal
          else if ((tcgpeekshort).and.(tcgguessshort.eqv..false.)) then
            call epolar1tcg2tershortreal
          else if ((tcgpeekshort).and.(tcgguessshort)) then
            call epolar1tcg2quatshortreal
          else
            if (rank.eq.0) write(iout,1000) 
            call fatal
          end if
        end if
      else
        if (tcgorder.eq.1) then
          if (use_pmecore) then
            call epolar1tcg1_2
          else
            call epolar1tcg1
          end if
        else if (tcgorder.eq.2) then
          if ((tcgprec).and.(tcgpeek).and.(tcgguess)) then
            if (use_pmecore) then
              call epolar1tcg2cinq_2
            else
              call epolar1tcg2cinq
            end if
          else if ((tcgguess.eqv..false.).and.(tcgpeek.eqv..false.)) 
     $ then 
            if (use_pmecore) then
              call epolar1tcg2_2
            else
              call epolar1tcg2
            end if
          else if ((tcgguess).and.(tcgpeek.eqv..false.)) then
            if (use_pmecore) then
              call epolar1tcg2bis_2
            else
              call epolar1tcg2bis
            end if
          else if ((tcgpeek).and.(tcgguess.eqv..false.)) then
            if (use_pmecore) then
              call epolar1tcg2ter_2
            else
              call epolar1tcg2ter
            end if
          else if ((tcgpeek).and.(tcgguess)) then
            if (use_pmecore) then
              call epolar1tcg2quat_2
            else
              call epolar1tcg2quat
            end if
          else
            if (rank.eq.0) write(iout,1000) 
            call fatal
          end if
        end if
      end if
      return
      end
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine epolar1tcg1  --  Ewald polarization energy             ##
c     ##                               and grad with tcg1                   ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "epolar1tcg1" calculates the induced dipole polarization gradients with the tcg1 method
c  
      subroutine epolar1tcg1
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
      integer  :: nrhs
      logical :: isguess, peek, precond
      integer ::  kk, irhs,
     $            i, iipole,j, ierr,iglob,ilocrec,
     $            kkpole,kglob,kloc,betac

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, sp0, sp1, term, omega, spp1, ap1a0, ap11a,
     $            ap11, ap12, ap2a0, ap21, e
      real*8 :: time0,time1
      real*8, dimension(2) :: n0, t4, t1, a10, a11, a12, a21



      real*8, allocatable, dimension(:,:,:,:) ::adte
      real*8, allocatable, dimension(:,:,:) :: dEnedt,efi,r0,r0bis,
     $         mu0,Tr0,mu_tcg1,adme,ade,adtebis, 
     $         aefi, ar0, atr0, tefi,  
     $         tmpbigbis,
     $         mupeek, taefi, tmpbig, 
     $         oldr0,oldr0bis, PTr0, oldTr0, oldefi, oldaefi, tmpbig2,
     $         Tmu0, efibis,aefibis,oldTr0bis,
     $         mu0bis,oldefibis,
     $         fphie, fphiae, fphir0, fphitr0,
     $         fphitae, fphite, fphiar0, fphiatr0,
     $         fphimu0,ar0bis,
     $         Tr0rec,Tr0recbis,Tr0bis,denedtbis,Taefirec,Taefirecbis,
     $         Tmu0rec,Tmu0bis,Tmu0recbis,
     $         aTr0bis,Taefibis,mupeekbis,Tefibis,
     $         oldaefibis,Tarrarec,Tarrarecbis

      real*8, allocatable, dimension(:,:) :: torq_mu, torq_t, grad_ene,
     $        dEnedmu, dEnedr, cphi, arr_dEp, 
     $         arr_depbis,arr_dEd,arr_dedbis, 
     $         adtb,adtbbis, adebis, admebis,
     $         denedrbis,denedmubis,torq_mubis,torq_tbis,
     $         arrA,arrTA,arraTA,fphiarrded,fphiarrdep,fphiarrA,fphiaTA,
     $         fphiarrdTr0,arr_dtr0,arr_dtr0bis,
     $         arrAbis,arrTAbis,arraTAbis


      real*8, allocatable, dimension(:,:,:) :: buffermpimu1,buffermpimu2
      real*8, allocatable, dimension(:,:) :: buffermpi1,buffermpi2

      parameter (nrhs=2)

      allocate(adte(3,3,nrhs,npolebloc))
      allocate(denedt(3,3,npolebloc), adtebis(3,3,nblocrec),
     $         efibis(3,nrhs,npolerecloc))
      allocate(efi(3,nrhs,max(1,npolebloc)), r0(3,nrhs,npolebloc),
     $         r0bis(3,nrhs,npolerecloc),
     $       aefibis(3,nrhs,npolerecloc),
     $       oldaefibis(3,nrhs,npolerecloc),
     $       mu0bis(3,nrhs,npolerecloc),
     $       ar0bis(3,nrhs,npolerecloc),
     $         mu0(3,nrhs,npolebloc), Tr0(3,nrhs,npolebloc),
     $         Tr0bis(3,nrhs,npolerecloc),
     $         aTr0bis(3,nrhs,npolerecloc),
     $         mu_tcg1(3,nrhs,npolebloc), ade(3,nrhs,npolebloc),
     $         adme(3,nrhs,npolebloc),
     $         Tmu0bis(3,nrhs,npolerecloc),
     $         Taefibis(3,nrhs,npolerecloc),
     $        Tefibis(3,nrhs,npolerecloc),
     $         aefi(3,nrhs,npolebloc),
     $         ar0(3,nrhs,npolebloc), atr0(3,nrhs,npolebloc),
     $         tefi(3,nrhs,npolebloc),
     $       arr_dTr0(3,npolebloc),arr_dTr0bis(3,npolerecloc),
     $        tmpbigbis(3,nrhs,npolerecloc),
     $         mupeek(3,nrhs,npolebloc), taefi(3,nrhs,npolebloc),
     $         tmpbig(3,nrhs,npolebloc), 
     $         oldr0(3,nrhs,npolebloc),oldr0bis(3,nrhs,npolerecloc),
     $         PTr0(3,nrhs,npoleloc), oldTr0(3,nrhs,npoleloc),
     $         oldefi(3,nrhs,npoleloc),oldTr0bis(3,nrhs,npolerecloc),
     $         Tmu0(3,nrhs,npolebloc),Tr0rec(3,nrhs,npoleloc),
     $         Tr0recbis(3,nrhs,npolerecloc),
     $          Taefirec(3,nrhs,npoleloc),
     $         Tmu0rec(3,nrhs,npoleloc),
     $         Tmu0recbis(3,nrhs,npolerecloc),
     $         Taefirecbis(3,nrhs,npolerecloc),
     $         TarrArec(3,nrhs,npolebloc),
     $         TarrArecbis(3,nrhs,npolerecloc),
     $         mupeekbis(3,nrhs,npolerecloc),
     $         oldefibis(3,nrhs,npolerecloc),
     $         oldaefi(3,nrhs,npolebloc), tmpbig2(3,nrhs,npolebloc),
     $         fphiE(20,2,npolerecloc), 
     $         fphiae(20,2,npolerecloc),
     $         fphir0(20,2,npolerecloc),
     $         fphiTr0(20,2,npolerecloc),
     $         fphiar0(20,2,npolerecloc), 
     $         fphitae(20,2,npolerecloc),
     $         fphiatr0(20,2,npolerecloc),
     $         fphite(20,2,npolerecloc),
     $         fphiarrded(20,npolerecloc),
     $         fphimu0(20,2,npolerecloc),
     $         fphiarrdep(20,npolerecloc),
     $         fphiarrdtr0(20,npolerecloc),
     $         fphiarrA(20,npolerecloc),
     $         fphiaTA(20,npolerecloc))
      if (allocated(fphirec)) deallocate(fphirec)
      allocate(fphirec(20, npolerecloc))
      if (allocated(cphirec)) deallocate(cphirec)
      allocate(cphirec(10, npolerecloc))
      allocate(cphi(10, npoleloc))
      allocate(torq_mu(3,nbloc), torq_t(3,nbloc),
     $         grad_ene(3,npolebloc), denedr(3,npolebloc),
     $         denedmu(3,npolebloc), arr_ded(3,npolebloc), 
     $         arr_dedbis(3,npolerecloc),
     $         arr_depbis(3,npolerecloc),
     $         arr_dep(3,npolebloc),adtbbis(3,nlocrec),
     $         adtb(3,npolebloc), adebis(3,nlocrec),
     $         admebis(3,nlocrec),
     $         torq_mubis(3,nblocrec),torq_tbis(3,nblocrec))
      allocate (arrA(3,npolebloc),arrTA(3,npolebloc),
     $           arraTA(3,npolebloc),arrAbis(3,npolerecloc),
     $           arrTAbis(3,npolerecloc),
     $           arraTAbis(3,npolerecloc))
      allocate(denedrbis(3,nlocrec),denedmubis(3,nlocrec),
     $         denedtbis(3,3,nlocrec))

      allocate (buffermpi1(10,max(npoleloc,1)))
      allocate (buffermpi2(10,max(npolerecloc,1)))
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      time0 = mpi_wtime()

      f = electric/dielec

      precond = tcgprec
      isguess = tcgguess
      peek    = tcgpeek

      omega = tcgomega

      !peek params
      efi = 0d0
      mu0 = 0d0
      mupeek = 0d0

      !1. Prepare electric field
      call efld0_recip(cphi)
      call commrecdirfields(0,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(1,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(2,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      time0 = mpi_wtime()
      call efld0_direct(nrhs, efi)
      call commfield(nrhs, efi)
      time1 = mpi_wtime()

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
          efi(j,1,i)  = efi(j,1,i) 
     $      - cphi(j+1,i) + term*rpole(j+1,iipole)
          efi(j,2,i)  = efi(j,2,i) 
     $      - cphi(j+1,i) + term*rpole(j+1,iipole)
        end do
      end do
      call commdirdir(nrhs,0,efi,reqrec,reqsend)
      call commdirdir(nrhs,1,efi,reqrec,reqsend)
      call commdirdir(nrhs,2,efi,reqrec,reqsend)

      oldefi = efi

      call commrecdirdip(nrhs,0,efibis,efi,buffermpimu1,
     $ buffermpimu2,req2rec,req2send)
      call commrecdirdip(nrhs,1,efibis,efi,buffermpimu1,
     $ buffermpimu2,req2rec,req2send)
      call commrecdirdip(nrhs,2,efibis,efi,buffermpimu1,
     $ buffermpimu2,req2rec,req2send)
      oldefibis = efibis

      if (.not. precond) then    
         if (.not. isguess) then 
            if (.not. peek) then !---
               r0 = efi
               r0bis = efibis
               oldr0 = r0
               time0 = mpi_wtime()
               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)
               time1 = mpi_wtime()

               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $            fphir0)

               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $           + Tr0rec(:,:,1:npoleloc) -term*r0(:,:,1:npoleloc)

               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
c               

               oldtr0 = Tr0
               fphie = fphir0
               time0 = mpi_wtime()
               call fftthatplz2(Tr0, Tr0bis, fphiTr0)
               time1 = mpi_wtime()
            else if (peek) then  !--P
               r0 = efi
               oldr0 = r0
               call diagvec(nrhs, efi, aefi)
               oldaefi = aefi
               call dbletmatxb_pme(nrhs,.true.,r0, aefi, Tr0, Taefi)
               call commfield(nrhs,Tr0)
               call commfield(nrhs,Taefi)

               call commrecdirdip(nrhs,0,r0bis,r0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,r0bis,r0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,r0bis,r0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               efibis = r0bis
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 aefibis(:,:,i) = polarity(iipole)*r0bis(:,:,i)
               end do

               call tmatxbrecipsave(r0, r0bis, nrhs, Tr0rec, Tr0recbis,
     $           fphir0)

               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc)- term*r0(:,:,1:npoleloc)

               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               fphiE = fphir0
               oldTr0 = Tr0
               call tmatxbrecipsave(aefi,aefibis,nrhs,Taefirec,
     $           Taefirecbis,fphiaE)

               call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)


               Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) 
     $          + Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)
               call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
               call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
               call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               call fftthatplz2(Tr0, Tr0bis, fphitr0)
               call fftthatplz2(Taefi, Taefibis, fphitae)
               call diagvec2(nrhs, r0, Tr0, ar0, aTr0)

               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 ar0bis(:,:,i)  = polarity(iipole)*r0bis(:,:,i)
                 aTr0bis(:,:,i)  = polarity(iipole)*Tr0bis(:,:,i)
               end do

               call fftthatplz2(ar0, ar0bis, fphiar0)
               call fftthatplz2(aTr0, aTr0bis, fphiatr0)
            end if
         else if (isguess) then  
            if (.not. peek) then !-G-
               call diagvec(nrhs, efi, mu0)
               aefi = mu0
               oldaefi = mu0
               call tmatxb_pme(nrhs, .true., mu0, Taefi)
               call commfield(nrhs,Taefi)

               call commrecdirdip(nrhs,0,aefibis,mu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,aefibis,mu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,aefibis,mu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               oldaefibis = aefibis
               mu0bis = aefibis

               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 if (polarity(iipole) .eq. 0d0) then
                    cycle
                 else
                   efibis(:,:,i) = aefibis(:,:,i)/polarity(iipole)
                 end if
               end do


               call tmatxbrecipsave(mu0,aefibis,nrhs,Taefirec,
     $            Taefirecbis,fphimu0)
               call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)

               Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) 
     $          + Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)
               call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
               call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
               call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               r0 = efi - Taefi 
               r0bis = efibis - Taefibis
               oldr0 = r0

               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)
               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $          fphir0)
               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)
               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)
               oldTr0 = Tr0 

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               call fftthatplz2(Tr0, Tr0bis, fphitr0)
               call fftthatplz2(efi, efibis, fphiE)

            else if (peek) then  !-GP
               call diagvec(nrhs, efi, aefi)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 aefibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
               end do
               mu0 = aefi
               mu0bis = aefibis
               oldaefi = aefi
               oldaefibis = aefibis
               call fftthatplz2(efi, efibis, fphiE)

               call tmatxb_pme(nrhs, .true., aefi, Taefi)
               call commfield(nrhs,Taefi)
               
               call tmatxbrecipsave(aefi,aefibis,nrhs,Taefirec,
     $          Taefirecbis,fphiae)
               fphimu0 = fphiae
               call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               
               Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) 
     $          + Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)

               call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
               call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
               call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call fftthatplz2(Taefi, Taefibis, fphiTaE)

               r0 = efi - Taefi
               r0bis = efibis - Taefibis
               oldr0 = r0
               oldr0bis = r0bis
               
               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)

               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $            fphir0)
               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)


               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               oldTr0bis = Tr0bis


               call fftthatplz2(Tr0, Tr0bis, fphitr0)

               call diagvec2(nrhs, r0, tr0, ar0, atr0)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 ar0bis(:,:,i)  = polarity(iipole)*r0bis(:,:,i)
                 aTr0bis(:,:,i)  = polarity(iipole)*Tr0bis(:,:,i)
               end do

               call fftthatplz2(ar0, ar0bis, fphiar0)
               call fftthatplz2(aTr0, aTr0bis, fphiaTr0)
            end if
         end if
      else if (precond) then
         if (.not. isguess) then 
            if (.not. peek) then !P--
               call diagvec(nrhs, efi, efi)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 efibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
               end do
               r0 = efi
               oldr0 = oldefi
               r0bis = efibis
               oldr0bis = oldefibis
               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)


               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $            fphir0)

               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0
               call diagvec(nrhs, Tr0, Tr0)

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               fphie = fphir0
               call fftthatplz2(Tr0, Tr0bis, fphiTr0)

            else if (peek) then  !P-P
               oldr0 = oldefi
               oldr0bis = oldefibis
               call diagvec(nrhs, efi, efi)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 efibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
               end do
               r0 = efi
               r0bis = efibis

               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)


               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $         fphir0)

               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)

               fphie = fphir0

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)


               oldtr0 = Tr0
               oldTr0bis = Tr0bis
               call diagvec(nrhs, Tr0,Tr0)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 Tr0bis(:,:,i) = Tr0bis(:,:,i)*polarity(iipole)
               end do

               Tefi = Tr0
               Tefibis = Tr0bis
               call fftthatplz2(Tr0, Tr0bis, fphitr0)
               fphiTE = fphiTr0
            end if
         else if (isguess) then
            if (.not. peek) then !PG-
               call diagvec2(nrhs, efi, oldefi, efi, mu0)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 efibis(:,:,i)  = efibis(:,:,i)*polarity(iipole)
                 mu0bis(:,:,i) = oldefi(:,:,i)*polarity(iipole)
               end do
               aefi = mu0
               oldaefi = mu0

               call tmatxb_pme(nrhs, .true., mu0, Tmu0)
               call commfield(nrhs,Tmu0)

               call commrecdirdip(nrhs,0,mu0bis,mu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,mu0bis,mu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,mu0bis,mu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               aefibis = mu0bis
               oldaefibis = aefibis

               call tmatxbrecipsave(mu0,mu0bis,nrhs,Tmu0rec,Tmu0recbis,
     $           fphimu0)
               call commrecdirsolv(nrhs,0,Tmu0recbis,Tmu0rec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tmu0recbis,Tmu0rec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)

               Tmu0(:,:,1:npoleloc) = Tmu0(:,:,1:npoleloc) 
     $          + Tmu0rec(:,:,1:npoleloc) - term*mu0(:,:,1:npoleloc)
               call commdirdir(nrhs,0,Tmu0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tmu0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tmu0,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Tmu0bis,Tmu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tmu0bis,Tmu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tmu0bis,Tmu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               fphiE = fphimu0 !"efi" est P.efi, et mu0 = P.efi !
               r0 = oldefi - Tmu0
               r0bis = oldefibis - Tmu0bis
               oldr0 = r0 
               call diagvec(nrhs, r0, r0)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 r0bis(:,:,i) = r0bis(:,:,i)*polarity(iipole)
               end do

               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)

               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $           fphir0)
               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)
               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0


               call diagvec(nrhs, Tr0, Tr0)

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               call fftthatplz2(Tr0, Tr0bis,fphitr0)
            else if (peek) then  !PGP
               call diagvec(nrhs, efi, efi)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 efibis(:,:,i) = polarity(iipole)*efibis(:,:,i) 
               end do
               mu0 = efi
               mu0bis = efibis
               aefi = mu0
               oldaefi = mu0
               aefibis = mu0bis
               oldaefibis = aefibis

               call tmatxb_pme(nrhs, .true., mu0, Tmu0)
               call commfield(nrhs,Tmu0)

               call tmatxbrecipsave(mu0,mu0bis,nrhs,Tmu0rec,Tmu0recbis,
     $           fphimu0)
               call commrecdirsolv(nrhs,0,Tmu0recbis,Tmu0rec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tmu0recbis,Tmu0rec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,
     $           buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
               fphiE = fphimu0

               Tmu0(:,:,1:npoleloc) = Tmu0(:,:,1:npoleloc) 
     $          + Tmu0rec(:,:,1:npoleloc) - term*mu0(:,:,1:npoleloc)
               call commdirdir(nrhs,0,Tmu0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tmu0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tmu0,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Tmu0bis,Tmu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tmu0bis,Tmu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tmu0bis,Tmu0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)


               oldr0 = oldefi - Tmu0
               r0bis = oldefibis - Tmu0bis
               oldr0bis = r0bis

               Tefi = Tmu0
               Tefibis = Tmu0bis

               call diagvec2(nrhs, oldr0, Tefi, r0, Tefi)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 r0bis(:,:,i) = polarity(iipole)*r0bis(:,:,i)
                 Tefibis(:,:,i) = polarity(iipole)*Tefibis(:,:,i)
               end do

               call fftthatplz2(Tefi, Tefibis, fphiTe)

               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)

               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $           fphir0)
               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi1,
     $           buffermpi2,reqrecdirrec,reqrecdirsend)

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu1,
     $          buffermpimu2,req2rec,req2send)

               oldTr0 = Tr0
               oldTr0bis = Tr0bis
               call diagvec(nrhs, Tr0, Tr0)
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 Tr0bis(:,:,i)  = polarity(iipole)*Tr0bis(:,:,i)
               end do


               call fftthatplz2(Tr0, Tr0bis, fphitr0)
            end if
         end if
      end if
       

      do irhs = 1,nrhs
         n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
         t1(irhs) = sprod(3*npoleloc,r0(:,irhs,:),oldTr0(:,irhs,:))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      t4 = n0/t1
      sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      do irhs = 1, nrhs
         mu_tcg1(:,irhs,:) =  mu0(:,irhs,:) + t4(irhs) * r0(:,irhs,:)
      end do

      if (peek .and. precond) then 
         mupeek(:,:,1:npolebloc) = omega*r0(:,:,1:npolebloc) 
     $     - omega*t4(1)*Tr0(:,:,1:npolebloc)
         mupeekbis(:,:,1:npolerecloc) = 
     $     omega*r0bis(:,:,1:npolerecloc) 
     $     - omega*t4(1)*Tr0bis(:,:,1:npolerecloc)
         mu_tcg1 = mu_tcg1 + mupeek
      else if (peek .and. .not. precond) then
         mupeek(:,:,1:npolebloc) = omega*ar0(:,:,1:npolebloc)
     $      - omega*t4(1)*aTr0(:,:,1:npolebloc)
         mupeekbis(:,:,1:npolerecloc) = omega*ar0bis(:,:,1:npolerecloc)
     $      - omega*t4(1)*aTr0bis(:,:,1:npolerecloc)
         mu_tcg1(:,:,1:npolebloc) = mu_tcg1(:,:,1:npolebloc) 
     $      + mupeek(:,:,1:npolebloc)
      end if
      e = sprod(3*npoleloc, oldefi(:,2,:), mu_tcg1(:,1,:))
      ep =  -.5d0*e*electric



      a10(1) = t4(1)
      a10(2) = 0d0
      a11(1) = 2d0*sp0/t1(1)
      a11(2) = t4(1)
      a12(1) = -2d0*sp0*n0(1)/(t1(1)*t1(1))
      a12(2) = 0d0
      a21 = - sp0*n0/(t1*t1)


      arr_dEd = a10(1)*efi(:,2,:) 
     $      + a11(1)*r0(:,1,:)
     $      + a12(1)*Tr0(:,1,:) 

      arr_dEdbis = a10(1)*efibis(:,2,:) 
     $      + a11(1)*r0bis(:,1,:) 
     $      + a12(1)*Tr0bis(:,1,:)

      arr_dEp = a11(2)*r0(:,1,:)
      arr_dEpbis = a11(2)*r0bis(:,1,:)
      arr_dTr0 = a21(1)*r0(:,1,:)
      arr_dTr0bis = a21(1)*r0bis(:,1,:)

      fphiarrded = a10(1)*fphie(:,2,:)
     $            + a11(1)*fphir0(:,1,:)
     $            + a12(1)*fphitr0(:,1,:)
      fphiarrdep = a11(2)*fphir0(:,1,:)
      fphiarrdtr0 = a21(1)*fphir0(:,1,:)
      time1 = mpi_wtime()
c      write(*,*) 'time build array = ',time1-time0

      if (peek .and. .not. precond) then
         spp1 = sprod(3*npoleloc, oldaefi(:,2,:), Tr0(:,1,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $      COMM_TINKER,ierr)

         ap1a0 = omega
         ap11a = -t4(1)*omega
         ap11 = -2d0*spp1*omega/t1(1)
         ap12 = 2d0*n0(1)*spp1*omega/(t1(1)*t1(1))
         ap2a0 = -t4(1)*omega
         ap21 = n0(1)*spp1*omega/(t1(1)*t1(1))

         arr_dEp = arr_dEp + mupeek(:,1,:)
         arr_dEpbis = arr_dEpbis + mupeekbis(:,1,:)
         arr_dEd = arr_dEd + ap1a0*aefi(:,2,:) 
     $             + ap11a*Taefi(:,2,:)
     $             + ap11*r0(:,1,:)
     $             + ap12*Tr0(:,1,:)
         arr_dEdbis = arr_dEdbis + ap1a0*aefibis(:,2,:) 
     $             + ap11a*Taefibis(:,2,:)
     $             + ap11*r0bis(:,1,:)
     $             + ap12*Tr0bis(:,1,:)
         arr_dTr0 = arr_dTr0 + ap21*r0(:,1,:)
     $                     + ap2a0*aefi(:,2,:)
         arr_dTr0bis = arr_dTr0bis + ap21*r0bis(:,1,:)
     $                     + ap2a0*aefibis(:,2,:)

         !To be tested
         fphiarrdep = fphiarrdep 
     $               + omega*fphiar0(:,1,:)
     $               - omega*t4(1)*fphiatr0(:,1,:)
         fphiarrded = fphiarrded 
     $               + ap1a0*fphiae(:,2,:)
     $               + ap11a*fphiTae(:,2,:)
     $               + ap11* fphir0(:,1,:)
     $               + ap12* fphiTr0(:,1,:)
         fphiarrdTr0 = fphiarrdTr0 
     $               + ap21 *fphir0(:,1,:)
     $               + ap2a0*fphiae(:,2,:)

      else if (peek .and. precond) then
         spp1 = sprod(3*npoleloc, oldefi(:,2,:), Tr0(:,1,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $      COMM_TINKER,ierr)

         ap1a0 = omega
         ap11a = -t4(1)*omega
         ap11 = -2d0*spp1*omega/t1(1)
         ap12 = 2d0*n0(1)*spp1*omega/(t1(1)*t1(1))
         ap2a0 = -t4(1)*omega
         ap21 = n0(1)*spp1*omega/(t1(1)*t1(1))

         arr_dEp = arr_dEp + mupeek(:,1,:)
         arr_dEpbis = arr_dEpbis + mupeekbis(:,1,:) 
         arr_dEd = arr_dEd + ap1a0*efi(:,2,:) 
     $             + ap11a*Tefi(:,2,:)
     $             + ap11*r0(:,1,:)
     $             + ap12*Tr0(:,1,:)
         arr_dEdbis = arr_dEdbis + ap1a0*efibis(:,2,:) 
     $             + ap11a*Tefibis(:,2,:)
     $             + ap11*r0bis(:,1,:)
     $             + ap12*Tr0bis(:,1,:)

         arr_dTr0 = arr_dTr0 
     $             + ap21*r0(:,1,:)
     $             + ap2a0*efi(:,2,:)

         arr_dTr0bis = arr_dTr0bis + ap21*r0bis(:,1,:)+
     $      ap2a0*efibis(:,2,:)

         !To be tested
         fphiarrdEp = fphiarrdEp 
     $             + omega*fphir0(:,1,:)
     $             - omega*t4(1)*fphitr0(:,1,:)
         fphiarrdEd = fphiarrdEd 
     $             + ap1a0*fphie(:,2,:)
     $             + ap11a*fphiTe(:,2,:)
     $             + ap11 *fphir0(:,1,:)
     $             + ap12 *fphiTr0(:,1,:)
         fphiarrdTr0 = fphiarrdTr0 
     $             + ap21 *fphir0(:,1,:)
     $             + ap2a0*fphie(:,2,:)


      end if

      time0 = mpi_wtime()
      denedr = 0d0
      denedmu = 0d0
      denedt = 0d0
      time1 = mpi_wtime()
c      write(*,*) 'time clear array  = ',time1-time0

      if (isguess) then
         arrA = arr_ded
         tmpbig(:,1,:) = arrA
         tmpbig(:,2,:) = arrA
         call commdirdir(nrhs,0,tmpbig,reqrec,reqsend)
         call commdirdir(nrhs,1,tmpbig,reqrec,reqsend)
         call commdirdir(nrhs,2,tmpbig,reqrec,reqsend)

         arrAbis = arr_dedbis
         tmpbigbis(:,1,:) = arrAbis
         tmpbigbis(:,2,:) = arrAbis
         fphiarrA = fphiarrded

         call tmatxb_pme(nrhs, .true., tmpbig, tmpbig2)
         call commfield(nrhs, tmpbig2)

         call tmatxbrecip(tmpbig,tmpbigbis,nrhs,TarrArec,TarrArecbis)

         call commrecdirsolv(nrhs,0,TarrArecbis,TarrArec,
     $     buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolv(nrhs,1,TarrArecbis,TarrArec,
     $     buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolv(nrhs,2,TarrArecbis,TarrArec,
     $     buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
         arrTA(:,1:npoleloc) = tmpbig2(:,1,1:npoleloc) 
     $    + TarrArec(:,1,1:npoleloc) - term*arrA(:,1:npoleloc)

         call commdirdir(1,0,arrTA,reqrec,reqsend)
         call commdirdir(1,1,arrTA,reqrec,reqsend)
         call commdirdir(1,2,arrTA,reqrec,reqsend)
c
         call commrecdirdip(1,0,arrTAbis,arrTA,buffermpimu1,
     $    buffermpimu2,req2rec,req2send)
         call commrecdirdip(1,1,arrTAbis,arrTA,buffermpimu1,
     $    buffermpimu2,req2rec,req2send)
         call commrecdirdip(1,2,arrTAbis,arrTA,buffermpimu1,
     $    buffermpimu2,req2rec,req2send)
         call diagvec(1, arrTA, arraTA)
         do i = 1, npolerecloc
           iipole = polerecglob(i)
           arraTAbis(:,i) = polarity(iipole)*arrTAbis(:,i)
         end do

         call fftthatplz(arraTa, arraTAbis, fphiaTA)

         arr_ded = mu0(:,2,:) + arrA - arraTA
         arr_dedbis = mu0bis(:,2,:) + arrAbis - arraTAbis
         fphiarrded= fphimu0(:,2,:) + fphiarrded- fphiaTA

         arr_dep = arr_dep + mu0(:,1,:)
         arr_depbis = arr_depbis + mu0bis(:,1,:)
         fphiarrdep= fphiarrdep+ fphimu0(:,1,:)

         time0 = mpi_wtime()
         call scalderfieldzmat3( arr_ded, arr_dep, 
     $                           arr_dtr0, r0(:,1,:),
     $                          -arrA, oldaefi(:,1,:),
     $                           ade, adme, adte, adtb)

         time1 = mpi_wtime()
         denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb
         denedmu=  adme(:,1,:)   + adme(:,2,:)    
         denedt =  adte(:,:,1,:) + adte(:,:,2,:)

         time0 = mpi_wtime()
         call scalderfieldzmatrec3(
     &                             arr_ded,arr_dedbis, fphiarrded,
     &                             arr_dep,arr_depbis,fphiarrdep,
     &                             arr_dtr0,arr_dtr0bis,fphiarrdtr0,
     &                             r0(:,1,:),r0bis(:,1,:),fphir0(:,1,:),
     &                           - arrA, - arrAbis, - fphiarrA,
     &                             oldaefi(:,1,:),oldaefibis(:,1,:),
     &                             fphimu0(:,1,:),adebis, admebis,
     &                             adtebis, adtbbis
     &                            )
         time1 = mpi_wtime()
c         write(*,*) 'time scalderrrec= ',time1-time0
         denedrbis = adebis + adtbbis
         denedmubis = admebis
         denedtbis = adtebis
      else if (.not. isguess) then

         time0 = mpi_wtime()
         call scalderfieldzmat1(
     &                          arr_ded, arr_dep, arr_dtr0,r0(:,1,:),
     &                          ade,adme,adte,adtb
     &                         )

         denedr =   ade(:,1,:)  + ade(:,2,:)  + adtb    
         denedmu=  adme(:,1,:)   + adme(:,2,:)    
         denedt =  adte(:,:,1,:) + adte(:,:,2,:)
         time1 = mpi_wtime()
c         write(*,*) 'time scalderreal = ',time1-time0

         time0 = mpi_wtime()
         call scalderfieldzmatrec1(
     &                             arr_ded,arr_dedbis, fphiarrded,
     &                             arr_dep,arr_depbis,fphiarrdep,
     &                             arr_dtr0,arr_dtr0bis,fphiarrdtr0,
     &                             r0(:,1,:),r0bis(:,1,:),fphir0(:,1,:),
     &                             adebis, admebis, adtebis, adtbbis
     &                            )

          denedrbis =  adebis + adtbbis      
          denedmubis=  admebis     
          denedtbis =  adtebis
         time1 = mpi_wtime()
c         write(*,*) 'time scalderrrec= ',time1-time0
      end if

      !N. Compute dEne/dmu and dEne/dtheta (given E=.5*mu.field)

      ! Do the contraction dEne/dmuP * dmuP/dr
      time0 = mpi_wtime()
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
      time1 = mpi_wtime()
c      write(*,*) 'contrac dir = ',time1 - time0
      time0 = mpi_wtime()

      ! Do the contraction dEne/dmuP * dmuP/dr
      call torquetcg_rec(torq_mubis,torq_tbis,denedmubis,denedtbis)

      deprec = 0d0

      do kk = 1, npolerecloc
         iipole = polerecglob(kk)
         iglob = ipole(iipole)
         ilocrec = locrec1(iglob)
         do betac = 1, 3
            deprec(betac,ilocrec) =  -.5d0*f*denedrbis(betac,kk)
         end do
      end do
      deprec(:,:) = deprec(:,:)-0.5d0*f*(torq_mubis(:,:)+torq_tbis(:,:))
      time1 = mpi_wtime()

c
      return
      end
