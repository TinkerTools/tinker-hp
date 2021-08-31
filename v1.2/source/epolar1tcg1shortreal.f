
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine epolar1tcg1shortreal  --  Ewald polarization energy     ##
c     ##                                       and grad with tcg1            ##
c     ##                                                                     ##
c     #########################################################################
c
c
c     "epolar1tcg1" calculates the induced dipole polarization gradients with the tcg1 method
c  
      subroutine epolar1tcg1shortreal
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
      integer :: kk, betac, irhs,ierr, kkpole,kglob,kloc
      integer, allocatable :: reqrec(:),reqsend(:)
      real*8 :: f, sprod, sp0, sp1, omega, spp1, ap1a0, ap11a,
     $            ap11, ap12, ap2a0, ap21, e
      real*8, dimension(2) :: n0, t4, t1, a10, a11, a12, a21
      real*8, allocatable, dimension(:,:,:,:) :: adte
      real*8, allocatable, dimension(:,:,:) :: dEnedt,efi,r0,
     $         mu0,Tr0,mu_tcg1,adme,ade,
     $         aefi, ar0, atr0, tefi,  
     $         mupeek, taefi, tmpbig, 
     $         oldr0, oldTr0, oldefi, oldaefi, tmpbig2,
     $         Tmu0

      real*8, allocatable, dimension(:,:) :: torq_mu, torq_t, grad_ene,
     $        dEnedmu, dEnedr, arr_dEp, arr_dEd, adtb, arrA, arrTA,
     $        arraTA,
     $        arr_dtr0

      parameter (nrhs=2)

      allocate(adte(3,3,nrhs,npolebloc))
      allocate(denedt(3,3,npolebloc))
      allocate(efi(3,nrhs,max(1,npolebloc)),r0(3,nrhs,max(1,npolebloc)),
     $         mu0(3,nrhs,max(1,npolebloc)), 
     $         Tr0(3,nrhs,max(1,npolebloc)),
     $         mu_tcg1(3,nrhs,max(1,npolebloc)), 
     $         ade(3,nrhs,max(1,npolebloc)),
     $         adme(3,nrhs,max(1,npolebloc)),
     $         aefi(3,nrhs,max(1,npolebloc)),
     $         ar0(3,nrhs,max(1,npolebloc)), 
     $         atr0(3,nrhs,max(1,npolebloc)),
     $         tefi(3,nrhs,max(1,npolebloc)),
     $         arr_dTr0(3,max(1,npolebloc)),
     $         mupeek(3,nrhs,max(1,npolebloc)), 
     $         taefi(3,nrhs,max(1,npolebloc)),
     $         tmpbig(3,nrhs,max(1,npolebloc)), 
     $         oldr0(3,nrhs,max(1,npolebloc)),
     $         oldTr0(3,nrhs,max(1,npoleloc)),
     $         oldefi(3,nrhs,max(1,npoleloc)),
     $         Tmu0(3,nrhs,max(1,npolebloc)),
     $         oldaefi(3,nrhs,max(1,npolebloc)), 
     $         tmpbig2(3,nrhs,max(1,npolebloc)))
      allocate(torq_mu(3,max(1,nbloc)), torq_t(3,max(1,nbloc)),
     $         grad_ene(3,max(1,npolebloc)),denedr(3,max(1,npolebloc)),
     $         denedmu(3,max(1,npolebloc)),arr_ded(3,max(1,npolebloc)), 
     $         arr_dep(3,max(1,npolebloc)),
     $         adtb(3,max(1,npolebloc)))
      allocate (arrA(3,max(1,npolebloc)),arrTA(3,max(1,npolebloc)),
     $           arraTA(3,max(1,npolebloc)))

      allocate (reqrec(nproc))
      allocate (reqsend(nproc))

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
      call efld0_shortreal(nrhs, efi)
      call commfieldshort(nrhs, efi)

      call commdirdirshort(nrhs,0,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,efi,reqrec,reqsend)

      oldefi = efi

      if (.not. precond) then    
         if (.not. isguess) then 
            if (.not. peek) then !---
               r0 = efi
               oldr0 = r0
               call tmatxb_shortreal(nrhs, .true., r0, Tr0)
               call commfieldshort(nrhs,Tr0)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

               oldtr0 = Tr0
            else if (peek) then  !--P
               r0 = efi
               oldr0 = r0
               call diagvec(nrhs, efi, aefi)
               oldaefi = aefi
               call dbletmatxb_shortreal(nrhs,.true.,r0, aefi, Tr0, 
     $         Taefi)
               call commfieldshort(nrhs,Tr0)
               call commfieldshort(nrhs,Taefi)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0
               call commdirdirshort(nrhs,0,Taefi,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Taefi,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Taefi,reqrec,reqsend)

               call diagvec2(nrhs, r0, Tr0, ar0, aTr0)
            end if
         else if (isguess) then  
            if (.not. peek) then !-G-
               call diagvec(nrhs, efi, mu0)
               aefi = mu0
               oldaefi = mu0
               call tmatxb_shortreal(nrhs, .true., mu0, Taefi)
               call commfieldshort(nrhs,Taefi)

               call commdirdirshort(nrhs,0,Taefi,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Taefi,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Taefi,reqrec,reqsend)

               r0 = efi - Taefi 
               oldr0 = r0

               call tmatxb_shortreal(nrhs, .true., r0, Tr0)
               call commfieldshort(nrhs,Tr0)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)
               oldTr0 = Tr0 
            else if (peek) then  !-GP
               call diagvec(nrhs, efi, aefi)
               mu0 = aefi
               oldaefi = aefi

               call tmatxb_shortreal(nrhs, .true., aefi, Taefi)
               call commfieldshort(nrhs,Taefi)
               
               call commdirdirshort(nrhs,0,Taefi,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Taefi,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Taefi,reqrec,reqsend)

               r0 = efi - Taefi
               oldr0 = r0
               
               call tmatxb_shortreal(nrhs, .true., r0, Tr0)
               call commfieldshort(nrhs,Tr0)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0
               call diagvec2(nrhs, r0, tr0, ar0, atr0)
            end if
         end if
      else if (precond) then
         if (.not. isguess) then 
            if (.not. peek) then !P--
               call diagvec(nrhs, efi, efi)
               r0 = efi
               oldr0 = oldefi
               call tmatxb_shortreal(nrhs, .true., r0, Tr0)
               call commfieldshort(nrhs,Tr0)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0
               call diagvec(nrhs, Tr0, Tr0)
            else if (peek) then  !P-P
               oldr0 = oldefi
               call diagvec(nrhs, efi, efi)
               r0 = efi

               call tmatxb_shortreal(nrhs, .true., r0, Tr0)
               call commfieldshort(nrhs,Tr0)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

               oldtr0 = Tr0
               call diagvec(nrhs, Tr0,Tr0)

               Tefi = Tr0
            end if
         else if (isguess) then
            if (.not. peek) then !PG-
               call diagvec2(nrhs, efi, oldefi, efi, mu0)
               aefi = mu0
               oldaefi = mu0

               call tmatxb_shortreal(nrhs, .true., mu0, Tmu0)
               call commfieldshort(nrhs,Tmu0)

               call commdirdirshort(nrhs,0,Tmu0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tmu0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tmu0,reqrec,reqsend)

               r0 = oldefi - Tmu0
               oldr0 = r0 
               call diagvec(nrhs, r0, r0)

               call tmatxb_shortreal(nrhs, .true., r0, Tr0)
               call commfieldshort(nrhs,Tr0)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0
               call diagvec(nrhs, Tr0, Tr0)
            else if (peek) then  !PGP
               call diagvec(nrhs, efi, efi)
               mu0 = efi
               aefi = mu0
               oldaefi = mu0

               call tmatxb_shortreal(nrhs, .true., mu0, Tmu0)
               call commfieldshort(nrhs,Tmu0)

               call commdirdirshort(nrhs,0,Tmu0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tmu0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tmu0,reqrec,reqsend)

               oldr0 = oldefi - Tmu0
               Tefi = Tmu0
               call diagvec2(nrhs, oldr0, Tefi, r0, Tefi)

               call tmatxb_shortreal(nrhs, .true., r0, Tr0)
               call commfieldshort(nrhs,Tr0)

               call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
               call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

               oldTr0 = Tr0
               call diagvec(nrhs, Tr0, Tr0)
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
     $      + a11(1)*r0(:,1,:)
     $      + a12(1)*Tr0(:,1,:) 

      arr_dEp = a11(2)*r0(:,1,:)
      arr_dTr0 = a21(1)*r0(:,1,:)


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
         arr_dEd = arr_dEd + ap1a0*aefi(:,2,:) 
     $             + ap11a*Taefi(:,2,:)
     $             + ap11*r0(:,1,:)
     $             + ap12*Tr0(:,1,:)
         arr_dTr0 = arr_dTr0 + ap21*r0(:,1,:)
     $                     + ap2a0*aefi(:,2,:)

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
         arr_dEd = arr_dEd + ap1a0*efi(:,2,:) 
     $             + ap11a*Tefi(:,2,:)
     $             + ap11*r0(:,1,:)
     $             + ap12*Tr0(:,1,:)

         arr_dTr0 = arr_dTr0 
     $             + ap21*r0(:,1,:)
     $             + ap2a0*efi(:,2,:)
      end if

      denedr = 0d0
      denedmu = 0d0
      denedt = 0d0

      if (isguess) then
         arrA = arr_ded
         tmpbig(:,1,:) = arrA
         tmpbig(:,2,:) = arrA
         call commdirdirshort(nrhs,0,tmpbig,reqrec,reqsend)
         call commdirdirshort(nrhs,1,tmpbig,reqrec,reqsend)
         call commdirdirshort(nrhs,2,tmpbig,reqrec,reqsend)

         call tmatxb_shortreal(nrhs, .true., tmpbig, tmpbig2)
         call commfieldshort(nrhs, tmpbig2)

         call commdirdirshort(1,0,arrTA,reqrec,reqsend)
         call commdirdirshort(1,1,arrTA,reqrec,reqsend)
         call commdirdirshort(1,2,arrTA,reqrec,reqsend)
c
         call diagvec(1, arrTA, arraTA)

         arr_ded = mu0(:,2,:) + arrA - arraTA

         arr_dep = arr_dep + mu0(:,1,:)

         call scalderfieldzmat3( arr_ded, arr_dep, 
     $                           arr_dtr0, r0(:,1,:),
     $                          -arrA, oldaefi(:,1,:),
     $                           ade, adme, adte, adtb)

         denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb
         denedmu=  adme(:,1,:)   + adme(:,2,:)    
         denedt =  adte(:,:,1,:) + adte(:,:,2,:)

      else if (.not. isguess) then

         call scalderfieldzmat1(
     &                          arr_ded, arr_dep, arr_dtr0,
     &                          r0(1:3,1,1:max(1,npolebloc)),
     &                          ade,adme,adte,adtb
     &                         )

         denedr =   ade(:,1,:)  + ade(:,2,:)  + adtb    
         denedmu=  adme(:,1,:)   + adme(:,2,:)    
         denedt =  adte(:,:,1,:) + adte(:,:,2,:)
      end if

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
      return
      end
