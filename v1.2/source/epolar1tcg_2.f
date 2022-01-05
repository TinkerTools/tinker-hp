c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine epolar1tcg1_2  --  Ewald polarization energy           ##
c     ##                                and grad with tcg1                  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "epolar1tcg1_2" calculates the induced dipole polarization gradients with the tcg1 method
c     whith separate cores for PME 
c  
      subroutine epolar1tcg1_2
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
      integer :: kk, betac, irhs,
     $            i, iipole, j, ierr,
     $            kkpole,kglob,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, sp0, sp1, term, omega, spp1, ap1a0, ap11a,
     $            ap11, ap12, ap2a0, ap21, e
      real*8, dimension(2) :: n0, t4, t1, a10, a11, a12, a21



      real*8, allocatable, dimension(:,:,:,:) :: adte
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


      real*8, allocatable, dimension(:,:,:) :: buffermpimu
      real*8, allocatable, dimension(:,:) :: buffermpi
      real*8 :: xx(1)

      parameter (nrhs=2)

      if (rank.le.ndir-1) then
        allocate(adte(3,3,nrhs,npolebloc))
        allocate(denedt(3,3,npolebloc),efi(3,nrhs,max(1,npolebloc)))
        allocate(r0(3,nrhs,npolebloc),mu0(3,nrhs,npolebloc),
     $ Tr0(3,nrhs,npolebloc),mu_tcg1(3,nrhs,npolebloc),
     $ ade(3,nrhs,npolebloc),adme(3,nrhs,npolebloc),
     $ aefi(3,nrhs,npolebloc),ar0(3,nrhs,npolebloc),
     $ atr0(3,nrhs,npolebloc),tefi(3,nrhs,npolebloc),
     $ mupeek(3,nrhs,npolebloc),taefi(3,nrhs,npolebloc),
     $ tmpbig(3,nrhs,npolebloc),oldr0(3,nrhs,npolebloc),
     $ PTr0(3,nrhs,npoleloc),oldTr0(3,nrhs,npoleloc),
     $ oldefi(3,nrhs,npoleloc),Tmu0(3,nrhs,npolebloc),
     $ Tr0rec(3,nrhs,npoleloc),Taefirec(3,nrhs,npoleloc),
     $ Tmu0rec(3,nrhs,npoleloc),TarrArec(3,nrhs,npolebloc),
     $ oldaefi(3,nrhs,npolebloc), tmpbig2(3,nrhs,npolebloc),
     $ cphi(10, npoleloc),torq_mu(3,nbloc), torq_t(3,nbloc),
     $ grad_ene(3,npolebloc),denedr(3,npolebloc),denedmu(3,npolebloc),
     $ arr_ded(3,npolebloc),adtb(3,npolebloc),arrA(3,npolebloc),
     $ arrTA(3,npolebloc),arraTA(3,npolebloc),arr_dTr0(3,npolebloc),
     $ arr_dep(3,npolebloc))
        allocate (buffermpi(10,max(npoleloc,1)))
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        !peek params
        efi = 0d0
        mu0 = 0d0
        mupeek = 0d0
      else
        allocate(r0bis(3,nrhs,npolerecloc),aefibis(3,nrhs,npolerecloc),
     $ oldaefibis(3,nrhs,npolerecloc),mu0bis(3,nrhs,npolerecloc),
     $ ar0bis(3,nrhs,npolerecloc),Tr0bis(3,nrhs,npolerecloc),
     $ aTr0bis(3,nrhs,npolerecloc),Tmu0bis(3,nrhs,npolerecloc),
     $ Taefibis(3,nrhs,npolerecloc),Tefibis(3,nrhs,npolerecloc),
     $ arr_dTr0bis(3,npolerecloc),tmpbigbis(3,nrhs,npolerecloc),
     $ oldr0bis(3,nrhs,npolerecloc),oldTr0bis(3,nrhs,npolerecloc),
     $ Tr0recbis(3,nrhs,npolerecloc),Tmu0recbis(3,nrhs,npolerecloc),
     $ Taefirecbis(3,nrhs,npolerecloc),TarrArecbis(3,nrhs,npolerecloc),
     $ mupeekbis(3,nrhs,npolerecloc),oldefibis(3,nrhs,npolerecloc),
     $ fphiE(20,2,npolerecloc),fphiae(20,2,npolerecloc),
     $ fphir0(20,2,npolerecloc),fphiTr0(20,2,npolerecloc),
     $ fphiar0(20,2,npolerecloc),fphitae(20,2,npolerecloc),
     $ fphiatr0(20,2,npolerecloc),fphite(20,2,npolerecloc),
     $ fphiarrded(20,npolerecloc),fphimu0(20,2,npolerecloc),
     $ fphiarrdep(20,npolerecloc),fphiarrdtr0(20,npolerecloc),
     $ fphiarrA(20,npolerecloc),fphiaTA(20,npolerecloc))
        if (allocated(fphirec)) deallocate(fphirec)
        allocate(fphirec(20,npolerecloc))
        if (allocated(cphirec)) deallocate(cphirec)
        allocate(cphirec(10,npolerecloc))
        allocate(arr_dedbis(3,npolerecloc),arr_depbis(3,npolerecloc),
     $ adtbbis(3,nlocrec),adebis(3,nlocrec),admebis(3,nlocrec),
     $ torq_mubis(3,nlocrec2),torq_tbis(3,nlocrec2))
        allocate (arrAbis(3,npolerecloc),arrTAbis(3,npolerecloc),
     $           arraTAbis(3,npolerecloc))
        allocate(denedrbis(3,nlocrec),denedmubis(3,nlocrec),
     $         denedtbis(3,3,nlocrec))
        allocate (buffermpi(10,max(npolerecloc,1)))
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        allocate(adtebis(3,3,nlocrec2),efibis(3,nrhs,npolerecloc))
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

      omega = tcgomega

      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,efi,reqrec,reqsend)
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        call efld0_direct(nrhs, efi)
        call commfield(nrhs, efi)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            efi(j,1,i)  = efi(j,1,i) 
     $        - cphi(j+1,i) + term*rpole(j+1,iipole)
            efi(j,2,i)  = efi(j,2,i) 
     $        - cphi(j+1,i) + term*rpole(j+1,iipole)
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
     $   reqrecdirrec,reqrecdirsend)
      end if

      if (rank.le.ndir-1) then
        call commrecdirdip(nrhs,1,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        if (.not. precond) then    
          if (.not. isguess) then 
            if (.not. peek) then !---
              r0 = efi
              oldr0 = r0

              call tmatxb_pme(nrhs, .true., r0, Tr0)
              call commfield(nrhs,Tr0)

              call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              
              Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) -term*r0(:,:,1:npoleloc)

              call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

              call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)

              oldtr0 = Tr0
            else if (peek) then  !--P
              r0 = efi
              oldr0 = r0
              call diagvec(nrhs, efi, aefi)
              oldaefi = aefi
              call dbletmatxb_pme(nrhs,.true.,r0, aefi, Tr0, Taefi)
              call commfield(nrhs,Tr0)
              call commfield(nrhs,Taefi)

              call commrecdirdip(nrhs,1,r0bis,r0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,r0bis,r0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              
              call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)

              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)

              Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $         + Tr0rec(:,:,1:npoleloc)- term*r0(:,:,1:npoleloc)

              call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

              call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)

              oldTr0 = Tr0
              call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

              Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) 
     $          + Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)
              call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
              call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
              call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

              call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)

              call diagvec2(nrhs, r0, Tr0, ar0, aTr0)
            end if
          else if (isguess) then  
            if (.not. peek) then !-G-
              call diagvec(nrhs, efi, mu0)
              aefi = mu0
              oldaefi = mu0
              call tmatxb_pme(nrhs, .true., mu0, Taefi)
              call commfield(nrhs,Taefi)
              call commrecdirdip(nrhs,1,aefibis,mu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,aefibis,mu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

              Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) 
     $         + Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)
              call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
              call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
              call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

              call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              r0 = efi - Taefi 
              oldr0 = r0
              call tmatxb_pme(nrhs, .true., r0, Tr0)
              call commfield(nrhs,Tr0)

              call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)

              Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $         + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)
              call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tr0,reqrec,reqsend)
              oldTr0 = Tr0 

              call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
            else if (peek) then  !-GP
              call diagvec(nrhs, efi, aefi)
              mu0 = aefi
              oldaefi = aefi
              call tmatxb_pme(nrhs, .true., aefi, Taefi)
              call commfield(nrhs,Taefi)
              call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              Taefi(:,:,1:npoleloc) = Taefi(:,:,1:npoleloc) 
     $         + Taefirec(:,:,1:npoleloc) - term*aefi(:,:,1:npoleloc)

              call commdirdir(nrhs,0,Taefi,reqrec,reqsend)
              call commdirdir(nrhs,1,Taefi,reqrec,reqsend)
              call commdirdir(nrhs,2,Taefi,reqrec,reqsend)

              call commrecdirdip(nrhs,1,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              r0 = efi - Taefi
              oldr0 = r0

              call tmatxb_pme(nrhs, .true., r0, Tr0)
              call commfield(nrhs,Tr0)

              call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $         + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)
              call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tr0,reqrec,reqsend)
              oldTr0 = Tr0

              call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
              call diagvec2(nrhs, r0, tr0, ar0, atr0)
            end if
          end if
        else if (precond) then
          if (.not. isguess) then 
            if (.not. peek) then !P--
              call diagvec(nrhs, efi, efi)
              r0 = efi
              oldr0 = oldefi
              call tmatxb_pme(nrhs, .true., r0, Tr0)
              call commfield(nrhs,Tr0)
              call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)

              Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $         + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

              call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

              oldTr0 = Tr0
              call diagvec(nrhs, Tr0, Tr0)

              call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
            else if (peek) then  !P-P
              oldr0 = oldefi
              call diagvec(nrhs, efi, efi)
              r0 = efi
              call tmatxb_pme(nrhs, .true., r0, Tr0)
              call commfield(nrhs,Tr0)
              call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)

              Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $         + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

              call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

              call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              oldtr0 = Tr0
              call diagvec(nrhs, Tr0,Tr0)
              Tefi = Tr0
            end if
          else if (isguess) then
            if (.not. peek) then !PG-
              call diagvec2(nrhs, efi, oldefi, efi, mu0)
              aefi = mu0
              oldaefi = mu0
              call tmatxb_pme(nrhs, .true., mu0, Tmu0)
              call commfield(nrhs,Tmu0)
              call commrecdirdip(nrhs,1,mu0bis,mu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,mu0bis,mu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirsolv(nrhs,0,Tmu0recbis,Tmu0rec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

              Tmu0(:,:,1:npoleloc) = Tmu0(:,:,1:npoleloc) 
     $         + Tmu0rec(:,:,1:npoleloc) - term*mu0(:,:,1:npoleloc)
              call commdirdir(nrhs,0,Tmu0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tmu0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tmu0,reqrec,reqsend)

              call commrecdirdip(nrhs,1,Tmu0bis,Tmu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tmu0bis,Tmu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              r0 = oldefi - Tmu0
              oldr0 = r0 
              call diagvec(nrhs, r0, r0)
              call tmatxb_pme(nrhs, .true., r0, Tr0)
              call commfield(nrhs,Tr0)

              call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)

              Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $         + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)
              call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
              call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

              oldTr0 = Tr0
              call diagvec(nrhs, Tr0, Tr0)

              call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
            else if (peek) then  !PGP
              call diagvec(nrhs, efi, efi)
               mu0 = efi
               aefi = mu0
               oldaefi = mu0
               call tmatxb_pme(nrhs, .true., mu0, Tmu0)
               call commfield(nrhs,Tmu0)

               call commrecdirsolv(nrhs,0,Tmu0recbis,Tmu0rec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

               Tmu0(:,:,1:npoleloc) = Tmu0(:,:,1:npoleloc) 
     $          + Tmu0rec(:,:,1:npoleloc) - term*mu0(:,:,1:npoleloc)
               call commdirdir(nrhs,0,Tmu0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tmu0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tmu0,reqrec,reqsend)

               call commrecdirdip(nrhs,1,Tmu0bis,Tmu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tmu0bis,Tmu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)

               oldr0 = oldefi - Tmu0
               Tefi = Tmu0
               call diagvec2(nrhs, oldr0, Tefi, r0, Tefi)
               call tmatxb_pme(nrhs, .true., r0, Tr0)
               call commfield(nrhs,Tr0)

               call commrecdirsolv(nrhs,0,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)

               Tr0(:,:,1:npoleloc) = Tr0(:,:,1:npoleloc) 
     $          + Tr0rec(:,:,1:npoleloc) - term*r0(:,:,1:npoleloc)

               call commdirdir(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdir(nrhs,2,Tr0,reqrec,reqsend)

               call commrecdirdip(nrhs,1,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)

               oldTr0 = Tr0
               call diagvec(nrhs, Tr0, Tr0)
            end if
          end if
        end if
      else
        call commrecdirdip(nrhs,0,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        call commrecdirdip(nrhs,2,efibis,efi,buffermpimu,
     $   buffermpimu,req2rec,req2send)
        oldefibis = efibis

        if (.not. precond) then    
          if (.not. isguess) then 
            if (.not. peek) then !---
              r0bis = efibis
              oldr0bis = r0bis
              call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $       fphir0)
              call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $      buffermpi,reqrecdirrec,reqrecdirsend)

              call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)

              fphie = fphir0
              call fftthatplz2(Tr0, Tr0bis, fphiTr0)
            else if (peek) then  !--P
              call commrecdirdip(nrhs,0,r0bis,r0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,r0bis,r0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              efibis = r0bis
              do i = 1, npolerecloc
                iipole = polerecglob(i)
                aefibis(:,:,i) = polarity(iipole)*r0bis(:,:,i)
              end do
              call tmatxbrecipsave(r0, r0bis, nrhs, Tr0rec, Tr0recbis,
     $           fphir0)
              call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)

              call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              fphiE = fphir0
              call tmatxbrecipsave(aefi,aefibis,nrhs,Taefirec,
     $           Taefirecbis,fphiaE)
              call commrecdirsolv(nrhs,0,Taefirecbis,Taefirec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

              call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)

              call fftthatplz2(Tr0, Tr0bis, fphitr0)
              call fftthatplz2(Taefi, Taefibis, fphitae)

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
              r0bis = efibis
              call commrecdirdip(nrhs,0,aefibis,mu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,aefibis,mu0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
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
     $           Taefirecbis,fphimu0)
              call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

              call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              r0bis = efibis - Taefibis
              call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $          fphir0)
              call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)

              call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)

              call fftthatplz2(Tr0, Tr0bis, fphitr0)
              call fftthatplz2(efi, efibis, fphiE)
            else if (peek) then  !-GP
              do i = 1, npolerecloc
                iipole = polerecglob(i)
                aefibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
              end do
              mu0bis = aefibis
              oldaefibis = aefibis
              call fftthatplz2(efi, efibis, fphiE)
              call tmatxbrecipsave(aefi,aefibis,nrhs,Taefirec,
     $         Taefirecbis,fphiae)
              fphimu0 = fphiae
              call commrecdirsolv(nrhs,1,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Taefirecbis,Taefirec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirdip(nrhs,0,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Taefibis,Taefi,buffermpimu,
     $         buffermpimu,req2rec,req2send)

              call fftthatplz2(Taefi, Taefibis, fphiTaE)

              r0bis = efibis - Taefibis
              oldr0bis = r0bis

              call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $           fphir0)
              call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)

              call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              oldTr0bis = Tr0bis
              call fftthatplz2(Tr0, Tr0bis, fphitr0)

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
              do i = 1, npolerecloc
                iipole = polerecglob(i)
                efibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
              end do
              r0bis = efibis
              oldr0bis = oldefibis

              call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $           fphir0)
              call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)

               fphie = fphir0
               call fftthatplz2(Tr0, Tr0bis, fphiTr0)
            else if (peek) then  !P-P
              oldr0bis = oldefibis
              do i = 1, npolerecloc
                iipole = polerecglob(i)
                efibis(:,:,i) = polarity(iipole)*efibis(:,:,i)
              end do
              r0bis = efibis
              call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $        fphir0)
              call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $          buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)

              fphie = fphir0

              call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)
              call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $         buffermpimu,req2rec,req2send)

               oldTr0bis = Tr0bis
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 Tr0bis(:,:,i) = Tr0bis(:,:,i)*polarity(iipole)
               end do
               Tefibis = Tr0bis
               call fftthatplz2(Tr0, Tr0bis, fphitr0)
               fphiTE = fphiTr0
            end if
          else if (isguess) then
            if (.not. peek) then !PG-
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 efibis(:,:,i)  = efibis(:,:,i)*polarity(iipole)
                 mu0bis(:,:,i) = oldefi(:,:,i)*polarity(iipole)
               end do
               call commrecdirdip(nrhs,0,mu0bis,mu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               call commrecdirdip(nrhs,2,mu0bis,mu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               aefibis = mu0bis
               oldaefibis = aefibis
               call tmatxbrecipsave(mu0,mu0bis,nrhs,Tmu0rec,Tmu0recbis,
     $           fphimu0)
               call commrecdirsolv(nrhs,1,Tmu0recbis,Tmu0rec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,
     $           buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
               call commrecdirdip(nrhs,0,Tmu0bis,Tmu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tmu0bis,Tmu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)

               fphiE = fphimu0 !"efi" est P.efi, et mu0 = P.efi !
               r0bis = oldefibis - Tmu0bis
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 r0bis(:,:,i) = r0bis(:,:,i)*polarity(iipole)
               end do
               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $           fphir0)
               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)
               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)

               call fftthatplz2(Tr0, Tr0bis,fphitr0)
            else if (peek) then  !PGP
              do i = 1, npolerecloc
                iipole = polerecglob(i)
                efibis(:,:,i) = polarity(iipole)*efibis(:,:,i) 
              end do
              mu0bis = efibis
              aefibis = mu0bis
              oldaefibis = aefibis
              call tmatxbrecipsave(mu0,mu0bis,nrhs,Tmu0rec,Tmu0recbis,
     $          fphimu0)
              call commrecdirsolv(nrhs,1,Tmu0recbis,Tmu0rec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              call commrecdirsolv(nrhs,2,Tmu0recbis,Tmu0rec,
     $          buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
              fphiE = fphimu0

               call commrecdirdip(nrhs,0,Tmu0bis,Tmu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tmu0bis,Tmu0,buffermpimu,
     $          buffermpimu,req2rec,req2send)

               r0bis = oldefibis - Tmu0bis
               oldr0bis = r0bis

               Tefibis = Tmu0bis
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 r0bis(:,:,i) = polarity(iipole)*r0bis(:,:,i)
                 Tefibis(:,:,i) = polarity(iipole)*Tefibis(:,:,i)
               end do
               call fftthatplz2(Tefi, Tefibis, fphiTe)
               call tmatxbrecipsave(r0,r0bis,nrhs,Tr0rec,Tr0recbis,
     $           fphir0)

               call commrecdirsolv(nrhs,1,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)
               call commrecdirsolv(nrhs,2,Tr0recbis,Tr0rec,buffermpi,
     $           buffermpi,reqrecdirrec,reqrecdirsend)

               call commrecdirdip(nrhs,0,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)
               call commrecdirdip(nrhs,2,Tr0bis,Tr0,buffermpimu,
     $          buffermpimu,req2rec,req2send)

               oldTr0bis = Tr0bis
               do i = 1, npolerecloc
                 iipole = polerecglob(i)
                 Tr0bis(:,:,i)  = polarity(iipole)*Tr0bis(:,:,i)
               end do

               call fftthatplz2(Tr0, Tr0bis, fphitr0)
            end if
          end if
        end if
      end if

      if (rank.le.ndir-1) then
        denedr = 0d0
        denedmu = 0d0
        denedt = 0d0
      end if


      if (rank.le.ndir-1) then
        do irhs = 1,nrhs
           n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
           t1(irhs) = sprod(3*npoleloc,r0(:,irhs,:),oldTr0(:,irhs,:))
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        do irhs = 1, nrhs
           mu_tcg1(:,irhs,:) =  mu0(:,irhs,:) + t4(irhs) * r0(:,irhs,:)
        end do

        if (peek .and. precond) then 
          mupeek(:,:,1:npolebloc) = omega*r0(:,:,1:npolebloc) 
     $      - omega*t4(1)*Tr0(:,:,1:npolebloc)
          mu_tcg1 = mu_tcg1 + mupeek
        else if (peek .and. .not. precond) then
          mupeek(:,:,1:npolebloc) = omega*ar0(:,:,1:npolebloc)
     $      - omega*t4(1)*aTr0(:,:,1:npolebloc)
          mu_tcg1(:,:,1:npolebloc) = mu_tcg1(:,:,1:npolebloc) 
     $      + mupeek(:,:,1:npolebloc)
        end if
        e = sprod(3*npoleloc, oldefi(:,2,:), mu_tcg1(:,1,:))
        ep =  -.5d0*e*electric

      ! Separated the fastgrad paper's terms a1i depending on the 
      ! which fields scale they are taken again.
      ! Ex : < a11(1)*r0, Ed' > 
      !  vs. < a11(2)*r0, Ep' >
        a10(1) = t4(1)
        a10(2) = 0d0
        a11(1) = 2d0*sp0/t1(1)
        a11(2) = t4(1)
        a12(1) = -2d0*sp0*n0(1)/(t1(1)*t1(1))
        a12(2) = 0d0
        a21 = - sp0*n0/(t1*t1)

        arr_dEd = a10(1)*efi(:,2,:) 
     $        + a11(1)*r0(:,1,:)
     $        + a12(1)*Tr0(:,1,:) 
        arr_dEp = a11(2)*r0(:,1,:)
        arr_dTr0 = a21(1)*r0(:,1,:)

      else
        n0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        t4 = n0/t1
        sp0 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)
        sp1 = 0d0
        call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $     COMM_TINKER,ierr)

        a10(1) = t4(1)
        a10(2) = 0d0
        a11(1) = 2d0*sp0/t1(1)
        a11(2) = t4(1)
        a12(1) = -2d0*sp0*n0(1)/(t1(1)*t1(1))
        a12(2) = 0d0
        a21 = - sp0*n0/(t1*t1)

        if (peek .and. precond) then 
          mupeekbis(:,:,1:npolerecloc) = omega*r0bis(:,:,1:npolerecloc) 
     $     - omega*t4(1)*Tr0bis(:,:,1:npolerecloc)
        else if (peek .and. .not. precond) then
          mupeekbis(:,:,1:npolerecloc) = omega*ar0bis(:,:,1:npolerecloc)
     $      - omega*t4(1)*aTr0bis(:,:,1:npolerecloc)
        end if

        arr_dEdbis = a10(1)*efibis(:,2,:) 
     $        + a11(1)*r0bis(:,1,:) 
     $        + a12(1)*Tr0bis(:,1,:)
        arr_dEpbis = a11(2)*r0bis(:,1,:)
        arr_dTr0bis = a21(1)*r0bis(:,1,:)
        fphiarrded = a10(1)*fphie(:,2,:)
     $              + a11(1)*fphir0(:,1,:)
     $              + a12(1)*fphitr0(:,1,:)
        fphiarrdep = a11(2)*fphir0(:,1,:)
        fphiarrdtr0 = a21(1)*fphir0(:,1,:)
      end if

      if (peek .and. .not. precond) then
        if (rank.le.ndir-1) then
           spp1 = sprod(3*npoleloc, oldaefi(:,2,:), Tr0(:,1,:))
           call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $        COMM_TINKER,ierr)
           ap1a0 = omega
           ap11a = -t4(1)*omega
           ap11 = -2d0*spp1*omega/t1(1)
           ap12 = 2d0*n0(1)*spp1*omega/(t1(1)*t1(1))
           ap2a0 = -t4(1)*omega
           ap21 = n0(1)*spp1*omega/(t1(1)*t1(1))

           arr_dEp = arr_dEp + mupeek(:,1,:)
           arr_dEd = arr_dEd + ap1a0*aefi(:,2,:) 
     $               + ap11a*Taefi(:,2,:)
     $               + ap11*r0(:,1,:)
     $               + ap12*Tr0(:,1,:)
           arr_dTr0 = arr_dTr0 + ap21*r0(:,1,:)
     $                     + ap2a0*aefi(:,2,:)
        else
           spp1 = 0d0
           call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $        COMM_TINKER,ierr)
           ap1a0 = omega
           ap11a = -t4(1)*omega
           ap11 = -2d0*spp1*omega/t1(1)
           ap12 = 2d0*n0(1)*spp1*omega/(t1(1)*t1(1))
           ap2a0 = -t4(1)*omega
           ap21 = n0(1)*spp1*omega/(t1(1)*t1(1))

           arr_dEpbis = arr_dEpbis + mupeekbis(:,1,:)
           arr_dEdbis = arr_dEdbis + ap1a0*aefibis(:,2,:) 
     $               + ap11a*Taefibis(:,2,:)
     $               + ap11*r0bis(:,1,:)
     $               + ap12*Tr0bis(:,1,:)
           arr_dTr0bis = arr_dTr0bis + ap21*r0bis(:,1,:)
     $                       + ap2a0*aefibis(:,2,:)
           fphiarrdep = fphiarrdep 
     $                 + omega*fphiar0(:,1,:)
     $                 - omega*t4(1)*fphiatr0(:,1,:)
           fphiarrded = fphiarrded 
     $                 + ap1a0*fphiae(:,2,:)
     $                 + ap11a*fphiTae(:,2,:)
     $                 + ap11* fphir0(:,1,:)
     $                 + ap12* fphiTr0(:,1,:)
           fphiarrdTr0 = fphiarrdTr0 
     $                 + ap21 *fphir0(:,1,:)
     $                 + ap2a0*fphiae(:,2,:)
        end if
      else if (peek .and. precond) then
        if (rank.le.ndir-1) then
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
          arr_dEd = arr_dEd + ap1a0*efi(:,2,:) 
     $             + ap11a*Tefi(:,2,:)
     $             + ap11*r0(:,1,:)
     $             + ap12*Tr0(:,1,:)
          arr_dTr0 = arr_dTr0 
     $             + ap21*r0(:,1,:)
     $             + ap2a0*efi(:,2,:)
        else
          spp1 = 0d0
          call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $       COMM_TINKER,ierr)
          ap1a0 = omega
          ap11a = -t4(1)*omega
          ap11 = -2d0*spp1*omega/t1(1)
          ap12 = 2d0*n0(1)*spp1*omega/(t1(1)*t1(1))
          ap2a0 = -t4(1)*omega
          ap21 = n0(1)*spp1*omega/(t1(1)*t1(1))
          arr_dEpbis = arr_dEpbis + mupeekbis(:,1,:) 
          arr_dEdbis = arr_dEdbis + ap1a0*efibis(:,2,:) 
     $             + ap11a*Tefibis(:,2,:)
     $             + ap11*r0bis(:,1,:)
     $             + ap12*Tr0bis(:,1,:)
         arr_dTr0bis = arr_dTr0bis + ap21*r0bis(:,1,:)+
     $      ap2a0*efibis(:,2,:)
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
      end if


      if (isguess) then
        if (rank.le.ndir-1) then
          arrA = arr_ded
          tmpbig(:,1,:) = arrA
          tmpbig(:,2,:) = arrA
          call commdirdir(nrhs,0,tmpbig,reqrec,reqsend)
          call commdirdir(nrhs,1,tmpbig,reqrec,reqsend)
          call commdirdir(nrhs,2,tmpbig,reqrec,reqsend)

          call tmatxb_pme(nrhs, .true., tmpbig, tmpbig2)
          call commfield(nrhs, tmpbig2)

          call commrecdirsolv(nrhs,0,TarrArecbis,TarrArec,
     $      buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,TarrArecbis,TarrArec,
     $      buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

          arrTA(:,1:npoleloc) = tmpbig2(:,1,1:npoleloc) 
     $     + TarrArec(:,1,1:npoleloc) - term*arrA(:,1:npoleloc)
          call commdirdir(1,0,arrTA,reqrec,reqsend)
          call commdirdir(1,1,arrTA,reqrec,reqsend)
          call commdirdir(1,2,arrTA,reqrec,reqsend)

          call commrecdirdip(1,1,arrTAbis,arrTA,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(1,2,arrTAbis,arrTA,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call diagvec(1, arrTA, arraTA)

          arr_ded = mu0(:,2,:) + arrA - arraTA
          arr_dep = arr_dep + mu0(:,1,:)

          call scalderfieldzmat3( arr_ded, arr_dep, 
     $                            arr_dtr0, r0(:,1,:),
     $                           -arrA, oldaefi(:,1,:),
     $                            ade, adme, adte, adtb)

          denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb
          denedmu=  adme(:,1,:)   + adme(:,2,:)    
          denedt =  adte(:,:,1,:) + adte(:,:,2,:)
        else
          arrAbis = arr_dedbis
          tmpbigbis(:,1,:) = arrAbis
          tmpbigbis(:,2,:) = arrAbis
          fphiarrA = fphiarrded
          call tmatxbrecip(tmpbig,tmpbigbis,nrhs,TarrArec,TarrArecbis)

          call commrecdirsolv(nrhs,1,TarrArecbis,TarrArec,
     $      buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,TarrArecbis,TarrArec,
     $     buffermpi,buffermpi,reqrecdirrec,reqrecdirsend)

          call commrecdirdip(1,0,arrTAbis,arrTA,buffermpimu,
     $     buffermpimu,req2rec,req2send)
          call commrecdirdip(1,2,arrTAbis,arrTA,buffermpimu,
     $     buffermpimu,req2rec,req2send)

          do i = 1, npolerecloc
            iipole = polerecglob(i)
            arraTAbis(:,i) = polarity(iipole)*arrTAbis(:,i)
          end do
          call fftthatplz(arraTa, arraTAbis, fphiaTA)

          arr_dedbis = mu0bis(:,2,:) + arrAbis - arraTAbis
          fphiarrded= fphimu0(:,2,:) + fphiarrded- fphiaTA
          arr_depbis = arr_depbis + mu0bis(:,1,:)
          fphiarrdep= fphiarrdep+ fphimu0(:,1,:)


        
          call scalderfieldzmatrec3(
     &                              xx,arr_dedbis,fphiarrded,
     &                              xx,arr_depbis,fphiarrdep,
     &                              xx,arr_dtr0bis,fphiarrdtr0,xx,
     &                              r0bis(:,1,:), fphir0(:,1,:),xx,
     &                            - arrAbis, - fphiarrA,xx,
     &                              oldaefibis(:,1,:), fphimu0(:,1,:),
     &                              adebis, admebis, adtebis, adtbbis
     &                             )
          

          denedrbis = adebis + adtbbis
          denedmubis = admebis
          denedtbis = adtebis
        end if
      else if (.not. isguess) then
        if (rank.le.ndir-1) then
          call scalderfieldzmat1(
     &                           arr_ded, arr_dep, arr_dtr0,r0(:,1,:),
     &                           ade,adme,adte,adtb
     &                          )
          denedr =   ade(:,1,:)  + ade(:,2,:)  + adtb    
          denedmu=  adme(:,1,:)   + adme(:,2,:)    
          denedt =  adte(:,:,1,:) + adte(:,:,2,:)
        else
          call scalderfieldzmatrec1(
     &                              xx,arr_dedbis, fphiarrded,xx,
     &                              arr_depbis,fphiarrdep, xx,
     &                              arr_dtr0bis,fphiarrdtr0,xx,
     &                              r0bis(:,1,:), fphir0(:,1,:),adebis,
     &                              admebis, adtebis, adtbbis
     &                             )

          denedrbis =  adebis + adtbbis      
          denedmubis=  admebis     
          denedtbis =  adtebis
        end if
      end if


      if (rank.le.ndir-1) then
        
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
        call torquetcg_rec(torq_mubis,torq_tbis,denedmubis,denedtbis)

        deprec = 0d0

        do kk = 1, npolerecloc
           do betac = 1, 3
              deprec(betac,kk) =  -.5d0*f*denedrbis(betac,kk)
           end do
        end do
        deprec(:,:) = deprec(:,:)-0.5d0*f*(torq_mubis(:,:)+
     $    torq_tbis(:,:))
      end if

      return
      end
