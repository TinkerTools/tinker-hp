c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine epolar1tcg1shortreal  --  Ewald polarization energy    ##
c     ##                                       and grad with tcg1           ##
c     ##                                                                    ##
c     ########################################################################
c
c     "epolar1tcg2shortreal" calculates the induced dipole polarization gradients with the tcg2 method
c  
! Computes TCG2 dipoles with NO refinement (except the preconditioner if
! "precond" is TRUE)
! Allows for minimal number of fft and matvec
! # of matvec : 4 (efi, Tr0, T2r0, T3r0)
! # of ffts : 5 (efi, r0, Tr0, T2r0, T3r0)
      subroutine epolar1tcg2shortreal
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

      integer :: betac, kk,
     $           irhs, k,
     $           ierr,kglob,kkpole,
     $           kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      real*8 :: f, sprod, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1




      real*8, allocatable, dimension(:,:,:,:) :: adte
      real*8, allocatable, dimension(:,:,:) :: denedt, mu_tcg2,efi,
     $            Tr0,T2r0, T3r0,Tefi, oldr0, oldTr0, oldefi, ade, adme,
     $            r0
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu, grad_ene,
     $            denedmu,denedr, arr_ded, arr_dep,adtb, arr_dtr0, 
     $            arr_dTTr0, arr_dTT2r0
      real*8:: time0,time1
      parameter (nrhs=2)

      allocate(adte(3,3,nrhs, npolebloc))
      allocate(denedt(3,3,npolebloc))
      allocate(mu_tcg2(3,nrhs, npolebloc), efi(3,nrhs, npolebloc),
     $         Tr0(3,nrhs,npolebloc),T2r0(3,nrhs,npolebloc),
     $         r0(3,nrhs,npolebloc),
     $         T3r0(3,nrhs,npolebloc),
     $         Tefi(3,nrhs,npolebloc),
     $         ade(3,nrhs,npolebloc),      adme(3,nrhs,npolebloc),
     $         oldr0(3,nrhs,npolebloc), oldTr0(3,nrhs,npolebloc),
     $         oldefi(3,nrhs,npolebloc))
      allocate(grad_ene(3,npolebloc), torq_mu(3,nbloc),
     $         torq_t(3,nbloc), 
     $         denedmu(3,npolebloc),
     $         denedr(3,npolebloc),   arr_ded(3,npolebloc),
     $         arr_dep(3,npolebloc),
     $         adtb(3,npolebloc),
     $         arr_dtr0(3,npolebloc),
     $         arr_dTTr0(3,npolebloc),
     $         arr_dTT2r0(3,npolebloc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))

      time0 = mpi_wtime()
      precond = tcgprec
      
      f = electric/dielec

      !1. prepare arrays...
      efi = 0d0


      !1.1 Electric field
      call efld0_shortreal(nrhs, efi)
      call commfieldshort(nrhs, efi)

      call commdirdirshort(nrhs,0,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,efi,reqrec,reqsend)

      !1.2 r0 
      oldefi = efi

      if (precond) call diagvec(nrhs, efi, efi)
      
      oldr0 = oldefi
      r0 = efi

      !1.3. Tr0 (+Tefi)
      call tmatxb_shortreal(nrhs, .true., r0, Tr0)
      call commfieldshort(nrhs,Tr0)

      call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

      oldTr0 = Tr0

      if (precond) call diagvec(nrhs, Tr0 ,Tr0)

      Tefi = Tr0
    
      do irhs = 1,nrhs
         n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
         t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      t4 = n0/t1

      !1.4. T2r0
      call tmatxb_shortreal(nrhs, .true., Tr0, T2r0)
      call commfieldshort(nrhs,T2r0)

      call commdirdirshort(nrhs,0,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2r0,reqrec,reqsend)

      if (precond) call diagvec(nrhs, T2r0, T2r0)

      do k = 1,nrhs
         np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
         t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
         t3(k) = t1(k)*t8(k)
         t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
         g1(k)  = t10(k)/t3(k)
         mu_tcg2(:,k,:) =  (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $   - g1(k)*t4(k)*Tr0(:,k,:)
      end do
      sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      b1 = sp0 - g1(1)*sp1
      b2 = sp0*t2(1) - t4(1)*sp1


      !1.5. T3r0
      call tmatxb_shortreal(nrhs, .true., T2r0, T3r0)
      call commfieldshort(nrhs,T3r0)

      call commdirdirshort(nrhs,0,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T3r0,reqrec,reqsend)

      if (precond) call diagvec(nrhs, T3r0, T3r0)
      
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


      ! 2.2 Assemble arrays
      arr_dEp  = mu_tcg2(:,1,:)

      arr_dEd = a10*efi(:,2,:) 
     $            + a1m1*Tefi(:,2,:)
     $            + a11*r0(:,1,:)
     $            + a12*Tr0(:,1,:)
     $            + a13*T2r0(:,1,:) 
     $            + a14*T3r0(:,1,:)

      !newones
      arr_dTr0 = a20*efi(:,2,:)
     $            + a21*r0(:,1,:)
     $            + a22*Tr0(:,1,:)
     $            + a23*T2r0(:,1,:)

      arr_dTTr0 = a31*r0(:,1,:)
     $            + a32*Tr0(:,1,:)

      arr_dTT2r0 = a41*r0(:,1,:)

      
      ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
      ep = -.5d0*f*ep 
      time1 = mpi_wtime()
c      write(*,*) 'time solver = ',time1-time0

      denedr = 0d0
      denedmu = 0d0
      denedt = 0d0

      time0 = mpi_wtime()
      call scalderfieldzmat6(
     &                       arr_ded, arr_dep, arr_dtr0, r0(:,1,:),
     &                       arr_dTTr0,  Tr0(:,1,:), arr_dTT2r0,
     &                       T2r0(:,1,:), ade, adme, adte, adtb
     &                      )

      time1 = mpi_wtime()
c         write(*,*) 'time scalderreal = ',time1-time0
      denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb
      denedmu=   adme(:,1,:)   + adme(:,2,:)    
      denedt =  adte(:,:,1,:) + adte(:,:,2,:)

      !N. Compute dEne/dmu and dEne/dtheta (given E=.5*mu.field)

      ! Do the contraction dEne/dmuP * dmuP/dr
      time0 = mpi_wtime()
      call torquetcg_dir(torq_mu,torq_t,denedmu,denedt)
c
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

      subroutine epolar1tcg2bisshortreal
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
      integer :: betac, kk,
     $           irhs, k,
     $           ierr,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: req2rec(:),req2send(:)
      real*8 :: f, sprod, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, 
     $          omega

      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1




      real*8, allocatable, dimension(:,:,:,:) :: adte
      real*8, allocatable, dimension(:,:,:) :: denedt,
     $            mu_tcg2,efi,
     $            tmpbig, oldr0, Tr0, oldTr0, T2r0,
     $            mu0,r0,
     $            T3r0,oldT3r0, 
     $            Tefi,oldTefi, ade, adme,
     $            oldaefi,
     $            oldefi,
     $            tmpbig2
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu,grad_ene,
     $            denedmu, denedr,
     $            arr_ded,arr_dep,adtb, 
     $            arr_dTr0, arr_dTTr0,
     $            arr_dTT2r0,arrA,arrTA, 
     $            arraTA
      real*8, allocatable, dimension(:,:,:) :: buffermpimu1,buffermpimu2
      real*8, allocatable, dimension(:,:) :: buffermpi1,buffermpi2
      parameter (nrhs=2)


      allocate(adte(3,3,nrhs, npolebloc))
      allocate(denedt(3,3,npolebloc))
      allocate(mu_tcg2(3,nrhs, npolebloc),
     $         efi(3,nrhs, npolebloc),
     $         tmpbig(3,nrhs,npolebloc),
     $         oldr0(3,nrhs,npolebloc),
     $         Tr0(3,nrhs,npolebloc),
     $         oldTr0(3,nrhs,npolebloc),
     $         T2r0(3,nrhs,npolebloc),
     $         mu0(3,nrhs,npolebloc),
     $         r0(3,nrhs,npolebloc),
     $         T3r0(3,nrhs,npolebloc),
     $         oldT3r0(3,nrhs,npolebloc),
     $         Tefi(3,nrhs,npolebloc),
     $         oldTefi(3,nrhs,npolebloc),
     $         ade(3,nrhs,npolebloc), 
     $         adme(3,nrhs,npolebloc),
     $         oldaefi(3,nrhs,npolebloc),
     $         oldefi(3,nrhs,npolebloc), tmpbig2(3,nrhs,npolebloc))
      allocate(grad_ene(3,npolebloc), torq_mu(3,nbloc),
     $         torq_t(3,nbloc), denedmu(3,npolebloc),
     $         denedr(3,npolebloc), arr_ded(3,npolebloc),
     $         arr_dep(3,npolebloc),
     $         adtb(3,npolebloc), 
     $         arr_dtr0(3,npolebloc), arr_dttr0(3,npolebloc),
     $         arr_dtt2r0(3,npolebloc), arrA(3,npolebloc),
     $         arrTA(3,npolebloc), arraTA(3,npolebloc))
      allocate (buffermpi1(10,max(npoleloc,1)))
      allocate (buffermpi2(10,max(npolerecloc,1)))
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      
      precond = tcgprec
      f = electric/dielec

      !1. prepare arrays...
      efi = 0d0
      Tr0 = 0d0

      omega = tcgomega


      !1.1 Electric field
      call efld0_shortreal(nrhs, efi)
      call commfieldshort(nrhs, efi)

      call commdirdirshort(nrhs,0,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,efi,reqrec,reqsend)

      oldefi = efi


      if (precond) then
        call diagvec(nrhs, efi, efi)
      end if

      !1.2 r0 
      !guess
      call diagvec(nrhs, oldefi, mu0)
      call dbletmatxb_shortreal(nrhs, .true., efi, mu0, Tefi, tmpbig)
      call commfieldshort(nrhs,Tefi)
      call commfieldshort(nrhs,tmpbig)

      call commdirdirshort(nrhs,0,tmpbig,reqrec,reqsend)
      call commdirdirshort(nrhs,1,tmpbig,reqrec,reqsend)
      call commdirdirshort(nrhs,2,tmpbig,reqrec,reqsend)

      call commdirdirshort(nrhs,0,Tefi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tefi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tefi,reqrec,reqsend)

      oldTefi = Tefi
      oldaefi =  mu0
      r0 = oldefi - tmpbig
      oldr0 = r0

      if (precond) then
        call diagvec2(nrhs, Tefi, r0, Tefi, r0)
      end if

      !1.3. Tr0 (+Tefi)
      call tmatxb_shortreal(nrhs, .true., r0, Tr0)
      call commfieldshort(nrhs,Tr0)

      call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

      oldTr0 = Tr0
      if (precond) then
        call diagvec(nrhs, Tr0, Tr0)
      end if


      do irhs = 1,nrhs
         n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
         t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
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

      !1.4. T2r0
      call tmatxb_shortreal(nrhs, .true., Tr0, T2r0)
      call commfieldshort(nrhs,T2r0)

      call commdirdirshort(nrhs,0,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2r0,reqrec,reqsend)

      if (precond) then
        call diagvec(nrhs, T2r0, T2r0)
      end if

      do k = 1,nrhs
         np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
         t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
         t3(k) = t1(k)*t8(k)
         t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
         g1(k)  = t10(k)/t3(k)
         mu_tcg2(:,k,:) = mu0(:,k,:) + (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $   - g1(k)*t4(k)*Tr0(:,k,:)
      end do
      sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      b1 = sp0 - g1(1)*sp1
      b2 = sp0*t2(1) - t4(1)*sp1


      !1.5. T3r0
      call tmatxb_shortreal(nrhs, .true., T2r0, T3r0)
      call commfieldshort(nrhs,T3r0)

      call commdirdirshort(nrhs,0,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T3r0,reqrec,reqsend)

      oldT3r0 = T3r0

      if (precond) then
        call diagvec(nrhs, T3r0, T3r0)
      end if

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

      ! 2.2 Assemble arrays
      arr_dep =  a10*r0(:,1,:) 
     $            + a20*Tr0(:,1,:)

      arr_dEd = a10*efi(:,2,:) 
     $            + a1m1*Tefi(:,2,:)
     $            + a11*r0(:,1,:)
     $            + a12*Tr0(:,1,:)
     $            + a13*T2r0(:,1,:) 
     $            + a14*T3r0(:,1,:)

      !newones
      arr_dTr0 = a20*efi(:,2,:)
     $            + a21*r0(:,1,:)
     $            + a22*Tr0(:,1,:)
     $            + a23*T2r0(:,1,:)
      arr_dTTr0 = a31*r0(:,1,:)
     $            + a32*Tr0(:,1,:)

      arr_dTT2r0 = a41*r0(:,1,:)

      ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
      ep = -.5d0*f*ep 

      denedr = 0d0
      denedmu = 0d0
      denedt = 0d0

      ! new guess computation. FANCYYYYY
      arrA = arr_ded 
      tmpbig(:,1,:) = arrA
      tmpbig(:,2,:) = arrA

      call tmatxb_shortreal(nrhs, .true., tmpbig, tmpbig2)
      call commfieldshort(nrhs,tmpbig2)

      arrTA = tmpbig2(:,1,:)

      call commdirdirshort(1,0,arrTA,reqrec,reqsend)
      call commdirdirshort(1,1,arrTA,reqrec,reqsend)
      call commdirdirshort(1,2,arrTA,reqrec,reqsend)

      call diagvec(1, arrTA, arraTA)

      arr_ded = mu0(:,2,:) + arrA - arraTA

      arr_dep = arr_dep + mu0(:,1,:)

      !Contraction
      call scalderfieldzmat8(
     &                       arr_ded, arr_dep, arr_dtr0,  r0(:,1,:),
     &                       arr_dTTr0,  Tr0(:,1,:), arr_dTT2r0,
     &                       T2r0(:,1,:), -arrA,     oldaefi(:,1,:),
     &                       ade, adme, adte, adtb
     &                      )

      denedr =  ade(:,1,:)    + ade(:,2,:)    + adtb  
      denedmu=  adme(:,1,:)   + adme(:,2,:)    
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

      return
      end



! Obj : tcg2 with peek only, preconditioner available if 'precond' is
! TRUE
      subroutine epolar1tcg2tershortreal
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
      integer :: betac, kk, 
     $           irhs, k,
     $           ierr,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      real*8 :: f, sprod, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, ap10a,
     $          ap11a, ap12a, ap11, ap12, ap13, ap14, ap2a0, ap21a,
     $          ap21, ap22, ap23, ap3a0, ap31, ap32, ap41, spp1, spp2,
     $          omega, k1, k2, k3
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1




      real*8, allocatable, dimension(:,:,:,:) :: adte
      real*8, allocatable, dimension(:,:,:) :: denedt,mu_tcg2,
     $            efi, tefi, oldtefi, r0, t3efi, ade, adme, oldaefi,
     $            oldefi, atefi, at2efi, aefi, Taefi, T2aefi, mupeek,
     $            T2efi
      real*8, allocatable, dimension(:,:) ::
     $            arr_dTr0, arr_dTTr0, arr_dTT2r0
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu, grad_ene,
     $            denedmu, denedr, arr_ded, arr_dep, adtb
      parameter (nrhs=2)

      allocate(adte(3,3,nrhs, npolebloc))
      allocate(denedt(3,3,npolebloc))
      allocate(mu_tcg2(3,nrhs,npolebloc),
     $         efi(3,nrhs, npolebloc),
     $         tefi(3,nrhs,npolebloc), oldtefi(3,nrhs,npolebloc), 
     $         t2efi(3,nrhs,npolebloc),
     $         r0(3,nrhs,npolebloc),
     $         t3efi(3,nrhs,npolebloc),
     $         arr_dTr0(3,npolebloc),
     $         arr_dTTr0(3,npolebloc), 
     $         arr_dTT2r0(3,npolebloc),
     $         ade(3,nrhs,npolebloc), adme(3,nrhs,npolebloc),
     $         oldaefi(3,nrhs,npolebloc),
     $         oldefi(3,nrhs,npolebloc), 
     $         atefi(3,nrhs,npolebloc), at2efi(3,nrhs,npolebloc),
     $         aefi(3,nrhs,npolebloc), Taefi(3,nrhs,npolebloc),
     $         T2aefi(3,nrhs,npolebloc),
     $         mupeek(3,nrhs,npolebloc))
      allocate(grad_ene(3,npolebloc), torq_mu(3,nbloc),
     $         torq_t(3,nbloc),
     $         denedmu(3,npolebloc),
     $         denedr(3,npolebloc), arr_ded(3,npolebloc),
     $         arr_dep(3,npolebloc), 
     $         adtb(3,npolebloc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      
      f = electric/dielec
      precond = tcgprec


      !1. prepare arrays...
      efi = 0d0
      omega = tcgomega


      !1.1 Electric field
      call efld0_shortreal(nrhs, efi)
      call commfieldshort(nrhs, efi)
      call commdirdirshort(nrhs,0,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,efi,reqrec,reqsend)
      oldefi = efi

      !1.2 r0 
      r0 = efi
      oldefi = efi
   
      ! If preconditioner is used, redefine every vector A as P^-1.A to
      ! simplify further notations.
      ! Keep a "non-preconditioned" copy (oldA) to compute scalars used
      ! in the algorithm.
      if (precond) call diagvec(nrhs, r0, r0)

      !1.3. Tr0  (= Tefi) + Taefi
      call diagvec(nrhs, efi, aefi)
      oldaefi = aefi

      if (.not. precond) then
         call dbletmatxb_shortreal(nrhs, .true., r0, aefi, tefi, taefi)
         call commfieldshort(nrhs,Tefi)
         call commdirdirshort(nrhs,0,Tefi,reqrec,reqsend)
         call commdirdirshort(nrhs,1,Tefi,reqrec,reqsend)
         call commdirdirshort(nrhs,2,Tefi,reqrec,reqsend)
         call commfieldshort(nrhs,Taefi)
         call commdirdirshort(nrhs,0,Taefi,reqrec,reqsend)
         call commdirdirshort(nrhs,1,Taefi,reqrec,reqsend)
         call commdirdirshort(nrhs,2,Taefi,reqrec,reqsend)
      else if (precond) then 
         call tmatxb_shortreal(nrhs, .true., r0, Tefi)
         call commfieldshort(nrhs,Tefi)
         call commdirdirshort(nrhs,0,Tefi,reqrec,reqsend)
         call commdirdirshort(nrhs,1,Tefi,reqrec,reqsend)
         call commdirdirshort(nrhs,2,Tefi,reqrec,reqsend)

         Taefi = Tefi
      end if

      oldtefi = tefi
      if (precond) call diagvec(nrhs, tefi, tefi)

      do irhs = 1,nrhs
         n0(irhs) = sprod(3*npoleloc, oldefi(:,irhs,:), r0(:,irhs,:))
         t1(irhs) = sprod(3*npoleloc, oldefi(:,irhs,:), tefi(:,irhs,:))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      t4 = n0/t1

      !1.4. T2r0
!      if (.not.precond) then
      call dbletmatxb_shortreal(nrhs, .true., Tefi, Taefi, 
     $                                          T2efi, T2aefi)
      call commfield(nrhs,T2efi)
      call commdirdirshort(nrhs,0,T2efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2efi,reqrec,reqsend)
      call commfieldshort(nrhs,T2aefi)
      call commdirdirshort(nrhs,0,T2aefi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2aefi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2aefi,reqrec,reqsend)

      if (precond) call diagvec(nrhs, t2efi, t2efi)

      do k = 1,nrhs
         np1(k) = sprod(3*npoleloc, oldtefi(:,k,:), tefi(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
         t9(k) = sprod(3*npoleloc, oldtefi(:,k,:), t2efi(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
         t3(k) = t1(k)*t8(k)
         t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
         g1(k)  = t10(k)/t3(k)
         mu_tcg2(:,k,:) = (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $   - g1(k)*t4(k)*tefi(:,k,:)
      end do
      sp0 = sprod(3*npoleloc, r0(:,1,:), efi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      sp1 = sprod(3*npoleloc, tefi(:,1,:), efi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      b1 = sp0 - g1(1)*sp1
      b2 = sp0*t2(1) - t4(1)*sp1


      !1.5. t3efi
      call tmatxb_shortreal(nrhs, .true., t2efi, t3efi)
      call commfieldshort(nrhs,T3efi)
      call commdirdirshort(nrhs,0,T3efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T3efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T3efi,reqrec,reqsend)

      if (precond) call diagvec(nrhs, t3efi, t3efi)
c
      ! Peek-step
      if (.not. precond) then
         call diagvec3(nrhs, efi, tefi, t2efi, aefi, atefi, at2efi)
         mupeek = omega*aefi 
     $            - omega*(g1(1)*t2(1) + t4(1))*atefi 
     $            + omega*g1(1)*t4(1)*at2efi
      else if (precond) then
         mupeek = omega*r0
     $            - omega*(g1(1)*t2(1) + t4(1))*tefi 
     $            + omega*g1(1)*t4(1)*t2efi
      end if
      mu_tcg2 = mu_tcg2 + mupeek

      efi = r0
      
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


      !peek terms
      if (.not. precond) then
         spp1 = sprod(3*npoleloc, oldaefi(:,2,:), tefi(:,1,:))
         spp2 = sprod(3*npoleloc, oldaefi(:,2,:), t2efi(:,1,:))
      else if (precond) then
         spp1 = sprod(3*npoleloc, oldefi(:,2,:), tefi(:,1,:))
         spp2 = sprod(3*npoleloc, oldefi(:,2,:), t2efi(:,1,:))
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)

      k1 = omega*t4(1)*spp2 - omega*t2(1)*spp1
      k2 = omega*g1(1)*spp2 - omega*spp1
      k3 = -omega*g1(1)*spp1

      ap10a = omega
      ap11a = -omega*(t2(1)*g1(1) + t4(1))
      ap12a = omega*t4(1)*g1(1)
      ap2a0 = -omega*(g1(1)*t2(1) + t4(1))
      ap21a = omega*t4(1)*g1(1)
      ap3a0 = omega*g1(1)*t4(1)

      ap11 = k1
     $       *(-2d0*np1(1)/t3(1)
     $         - 2d0*np1(1)*np1(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t9(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(2d0/t1(1))
     $      +k3
     $       *(2d0*np1(1)/(t1(1)*t1(1)))
      ap12 = k1
     $       *(4d0*t1(1)/t3(1) 
     $        - 2d0*n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $        + 4d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $        - 2d0*t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(-2d0*n0(1)/(t1(1)*t1(1)))
     $      +k3
     $       *(-4d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap13 = k1
     $       *(-4d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1))
     $         - 2d0*n0(1)/t3(1))
     $      +k3
     $       *(2d0*n0(1)/(t1(1)*t1(1)))
      ap14 = k1
     $       *(2d0*t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap21 = k1
     $       *(2d0*t1(1)/t3(1) 
     $         - n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $         - t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(-n0(1)/(t1(1)*t1(1)))
     $      +k3
     $       *(-2d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap22 = k1
     $       *(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3
     $       *(n0(1)/(t1(1)*t1(1)))
      ap23 = k1
     $       *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap31 = k1
     $       *(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3
     $       *(n0(1)/(t1(1)*t1(1)))
      ap32 = k1
     $       *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap41 = k1
     $        *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))

      ! 2.2 Assemble arrays
      arr_dEp =(g1(1)*t2(1) 
     $         + t4(1))*efi(:,1,:) 
     $         - g1(1)*t4(1)*tefi(:,1,:)

      arr_dEd = a10*efi(:,2,:) 
     $            + a1m1*Tefi(:,2,:)
     $            + a11*efi(:,1,:)
     $            + a12*tefi(:,1,:)
     $            + a13*t2efi(:,1,:) 
     $            + a14*t3efi(:,1,:)

      arr_dTr0 = a20*efi(:,2,:)
     $            + a21*efi(:,1,:)
     $            + a22*tefi(:,1,:)
     $            + a23*t2efi(:,1,:)
      arr_dTTr0 = a31*efi(:,1,:)
     $            + a32*tefi(:,1,:)
      arr_dTT2r0 = a41*efi(:,1,:)

      ! Peek-step arrays
      if (.not. precond) then
         arr_dep = arr_dep 
     $            + ap10a*aefi(:,1,:)
     $            + ap11a*atefi(:,1,:) 
     $            + ap21a*at2efi(:,1,:)

         arr_ded = arr_ded 
     $             + ap10a*aefi(:,2,:)
     $             + ap11a*Taefi(:,2,:)
     $             + ap12a*T2aefi(:,2,:)
     $             + ap11*efi(:,1,:)
     $             + ap12*tefi(:,1,:)
     $             + ap13*t2efi(:,1,:)
     $             + ap14*t3efi(:,1,:)

         arr_dTr0 = arr_dTr0 
     $             + ap2a0*aefi(:,2,:)
     $             + ap21a*Taefi(:,2,:)
     $             + ap21*efi(:,1,:)
     $             + ap22*tefi(:,1,:)
     $             + ap23*t2efi(:,1,:)

         arr_dTTr0 = arr_dTTr0 
     $             + ap3a0*aefi(:,2,:)
     $             + ap31*efi(:,1,:)
     $             + ap32*tefi(:,1,:)

         arr_dTT2r0 = arr_dTT2r0 
     $             + ap41*efi(:,1,:)


      else if (precond) then
         arr_dep = arr_dep 
     $            + ap10a*efi(:,1,:)
     $            + ap11a*tefi(:,1,:) 
     $            + ap21a*t2efi(:,1,:)

         arr_ded = arr_ded 
     $             + ap10a*efi(:,2,:)
     $             + ap11a*Tefi(:,2,:)
     $             + ap12a*T2efi(:,2,:)
     $             + ap11*efi(:,1,:)
     $             + ap12*tefi(:,1,:)
     $             + ap13*t2efi(:,1,:)
     $             + ap14*t3efi(:,1,:)

         arr_dTr0 = arr_dTr0 
     $             + ap2a0*efi(:,2,:)
     $             + ap21a*Tefi(:,2,:)
     $             + ap21*efi(:,1,:)
     $             + ap22*tefi(:,1,:)
     $             + ap23*t2efi(:,1,:)

         arr_dTTr0 = arr_dTTr0 
     $             + ap3a0*efi(:,2,:)
     $             + ap31*efi(:,1,:)
     $             + ap32*tefi(:,1,:)

         arr_dTT2r0 = arr_dTT2r0 
     $             + ap41*efi(:,1,:)
      end if

      ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
      ep = -.5d0*f*ep 
      denedr = 0d0
      denedmu = 0d0
      denedt = 0d0


      call scalderfieldzmat6(
     &                       arr_ded, arr_dep, arr_dtr0, efi(:,1,:),
     &                       arr_dTTr0, Tefi(:,1,:), arr_dTT2r0,
     &                       T2efi(:,1,:), ade, adme, adte, adtb
     &                      )

      denedr =  ade(:,1,:)    + ade(:,2,:)  + adtb     
      denedmu=  adme(:,1,:)   + adme(:,2,:)    
      denedt =  adte(:,:,1,:) + adte(:,:,2,:)

c
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


! TCG2 with guess AND peek-step, NO preconditioner

      subroutine epolar1tcg2quatshortreal
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
      integer :: betac, kk,
     $           irhs, k,
     $           ierr,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      real*8 :: f, sprod, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, ap10a,
     $          ap11a, ap12a, ap11, ap12, ap13, ap14, ap2a0, ap21a,
     $          ap21, ap22, ap23, ap3a0, ap31, ap32, ap41, spp1, spp2,
     $          omega, k1, k2, k3
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1




      real*8, allocatable, dimension(:,:,:,:) :: adte
      real*8, allocatable, dimension(:,:,:) :: denedt,mu_tcg2,
     $            efi, oldr0,Tr0, oldTr0,T2r0, mu0,r0, T3r0, Tefi,
     $            oldTefi, ade,adme,oldaefi,Tarr_dr0, aTarr_dr0,oldefi,
     $            ar0,atr0,at2r0,aefi, Taefi,T2aefi,mupeek, tmpbig
      real*8, allocatable, dimension(:,:) ::
     $            arr_dTT2r0, arr_dTr0,arr_dTTr0, arr_dr0
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu, grad_ene,
     $            denedmu, denedr, arr_ded, arr_dep, adtb, 
     $            tmpsmall
      parameter (nrhs=2)

      allocate(adte(3,3,nrhs, npolebloc))
      allocate(denedt(3,3,npolebloc)) 
      allocate(mu_tcg2(3,nrhs,npolebloc),
     $         efi(3,nrhs,npolebloc),
     $         oldr0(3,nrhs,npolebloc),
     $         Tr0(3,nrhs,npolebloc),
     $         oldTr0(3,nrhs,npolebloc),
     $         T2r0(3,nrhs,npolebloc),
     $         mu0(3,nrhs,npolebloc),
     $         r0(3,nrhs,npolebloc),
     $         T3r0(3,nrhs,npolebloc),
     $         arr_dTr0(3,npolebloc),
     $         arr_dTTr0(3,npolebloc),
     $         arr_dTT2r0(3,npolebloc),
     $         Tefi(3,nrhs,npolebloc),
     $         oldTefi(3,nrhs,npolebloc),
     $         ade(3,nrhs,npolebloc), adme(3,nrhs,npolebloc),
     $         oldaefi(3,nrhs,npolebloc),
     $         arr_dr0(3,npolebloc),
     $         Tarr_dr0(3,nrhs,npolebloc),
     $         atarr_dr0(3,nrhs,npolebloc),
     $         oldefi(3,nrhs,npolebloc),
     $         ar0(3,nrhs,npolebloc),
     $         atr0(3,nrhs,npolebloc),
     $         at2r0(3,nrhs,npolebloc),
     $         aefi(3,nrhs,npolebloc),
     $         Taefi(3,nrhs,npolebloc),
     $         T2aefi(3,nrhs,npolebloc),
     $         mupeek(3,nrhs,npolebloc),
     $         tmpbig(3,nrhs,npolebloc))
      allocate(grad_ene(3,npolebloc), torq_mu(3,nbloc),
     $         torq_t(3,nbloc), denedmu(3,npolebloc),
     $         denedr(3,npolebloc),arr_ded(3,npolebloc),
     $         arr_dep(3,npolebloc),
     $         adtb(3,npolebloc),
     $         tmpsmall(3,npolebloc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      
      f = electric/dielec
      precond = tcgprec

      !1. prepare arrays...
      efi = 0d0

      omega = tcgomega


      !1.1 Electric field
      call efld0_shortreal(nrhs, efi)
      call commfieldshort(nrhs, efi)
      call commdirdirshort(nrhs,0,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,efi,reqrec,reqsend)

      oldefi = efi

      call diagvec(nrhs, efi, aefi)

      mu0 = aefi
      oldaefi = aefi

      call dbletmatxb_shortreal(nrhs, .true., aefi, efi, Taefi, Tefi)
      call commfieldshort(nrhs,Tefi)
      call commfieldshort(nrhs,Taefi)
      call commdirdirshort(nrhs,0,Taefi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Taefi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Taefi,reqrec,reqsend)
      call commdirdirshort(nrhs,0,Tefi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tefi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tefi,reqrec,reqsend)

      oldTefi = Tefi
      r0 = efi - Taefi
      oldr0 = r0

      call dbletmatxb_shortreal(nrhs, .true., Taefi, r0, T2aefi, Tr0)
      call commfieldshort(nrhs,T2aefi)
      call commfieldshort(nrhs,Tr0)
      call commdirdirshort(nrhs,0,T2aefi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2aefi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2aefi,reqrec,reqsend)
      call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)

      oldTr0 = Tr0

      call tmatxb_shortreal(nrhs, .true., Tr0, T2r0)
      call commfieldshort(nrhs,T2r0)
      call commdirdirshort(nrhs,0,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2r0,reqrec,reqsend)

      call tmatxb_shortreal(nrhs, .true., T2r0, T3r0)
      call commfieldshort(nrhs,T3r0)
      call commdirdirshort(nrhs,0,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T3r0,reqrec,reqsend)

      call diagvec3(nrhs,r0,Tr0,T2r0,ar0,aTr0,aT2r0)
      do irhs = 1,nrhs
         n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
         t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      t4 = n0/t1

      do k = 1,nrhs
         np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
         t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
         t3(k) = t1(k)*t8(k)
         t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
         g1(k)  = t10(k)/t3(k)
         mu_tcg2(:,k,:) = mu0(:,k,:) + (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $   - g1(k)*t4(k)*Tr0(:,k,:)
      end do
      sp0 = sprod(3*npoleloc, r0(:,1,:), efi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      sp1 = sprod(3*npoleloc, Tr0(:,1,:), efi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      b1 = sp0 - g1(1)*sp1
      b2 = sp0*t2(1) - t4(1)*sp1

      !peek  no prec
      mupeek = omega*ar0 
     $        - omega*(g1(1)*t2(1) + t4(1))*aTr0 
     $        + omega*g1(1)*t4(1)*at2r0
      mu_tcg2 = mu_tcg2 + mupeek


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


      !peek noprec
      spp1 = sprod(3*npoleloc, oldaefi(:,2,:), Tr0(:,1,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      spp2 = sprod(3*npoleloc, oldaefi(:,2,:), T2r0(:,1,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)

      k1 = omega*t4(1)*spp2 - omega*t2(1)*spp1
      k2 = omega*g1(1)*spp2 - omega*spp1
      k3 = -omega*g1(1)*spp1

      ap10a = omega
      ap11a = -omega*(t2(1)*g1(1) + t4(1))
      ap12a = omega*t4(1)*g1(1)
      ap2a0 = -omega*(g1(1)*t2(1) + t4(1))
      ap21a = omega*t4(1)*g1(1)
      ap3a0 = omega*g1(1)*t4(1)

      ap11 = k1
     $       *(-2d0*np1(1)/t3(1)
     $         - 2d0*np1(1)*np1(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t9(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(2d0/t1(1))
     $      +k3
     $       *(2d0*np1(1)/(t1(1)*t1(1)))
      ap12 = k1
     $       *(4d0*t1(1)/t3(1) 
     $        - 2d0*n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $        + 4d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $        - 2d0*t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(-2d0*n0(1)/(t1(1)*t1(1)))
     $      +k3
     $       *(-4d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap13 = k1
     $       *(-4d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1))
     $         - 2d0*n0(1)/t3(1))
     $      +k3
     $       *(2d0*n0(1)/(t1(1)*t1(1)))
      ap14 = k1
     $       *(2d0*t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap21 = k1
     $       *(2d0*t1(1)/t3(1) 
     $         - n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $         - t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(-n0(1)/(t1(1)*t1(1)))
     $      +k3
     $       *(-2d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap22 = k1
     $       *(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3
     $       *(n0(1)/(t1(1)*t1(1)))
      ap23 = k1
     $       *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap31 = k1
     $       *(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3
     $       *(n0(1)/(t1(1)*t1(1)))
      ap32 = k1
     $       *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap41 = k1
     $        *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))

      ! 2.2 Assemble arrays
      arr_dEp  = mu_tcg2(:,1,:)

      arr_dEd = a10*efi(:,2,:) 
     $            + a1m1*Tefi(:,2,:)
     $            + a11*r0(:,1,:)
     $            + a12*Tr0(:,1,:)
     $            + a13*T2r0(:,1,:) 
     $            + a14*T3r0(:,1,:)

      arr_dTr0 = a20*efi(:,2,:)
     $            + a21*r0(:,1,:)
     $            + a22*Tr0(:,1,:)
     $            + a23*T2r0(:,1,:)

      arr_dTTr0 = a31*r0(:,1,:)
     $            + a32*Tr0(:,1,:)

      arr_dTT2r0 = a41*r0(:,1,:)


      !Peek noprec
      arr_ded = arr_ded 
     $          + ap10a*aefi(:,2,:)
     $          + ap11a*Taefi(:,2,:)
     $          + ap12a*T2aefi(:,2,:)
     $          + ap11*r0(:,1,:)
     $          + ap12*Tr0(:,1,:)
     $          + ap13*T2r0(:,1,:)
     $          + ap14*T3r0(:,1,:)

      arr_dTr0 = arr_dTr0 
     $          + ap2a0*aefi(:,2,:)
     $          + ap21a*Taefi(:,2,:)
     $          + ap21*r0(:,1,:)
     $          + ap22*Tr0(:,1,:)
     $          + ap23*T2r0(:,1,:)

      arr_dTTr0 = arr_dTTr0 
     $          + ap3a0*aefi(:,2,:)
     $          + ap31*r0(:,1,:)
     $          + ap32*Tr0(:,1,:)

      arr_dTT2r0 = arr_dTT2r0 
     $          + ap41*r0(:,1,:)

      ep = sprod(3*npoleloc, mu_tcg2(:,1,:), efi(:,2,:))
      ep = -.5d0*f*ep 

      denedr = 0d0
      denedmu = 0d0
      denedt = 0d0


      !GUESS
      ! There is some fundamental optimization to be done here :
      ! this matvec product should be avoidable through some smart
      ! combinations... or maybe not ?
      ! Either we tmatxb all arr_dr0, or we combine already tmatxb'd
      ! vectors composing arr_dr0. But for this second solution,
      ! we'll need tmatxb(T3r0) => a fourth tmatxb

      atarr_dr0 = 0d0


      arr_dr0 = arr_ded
      tmpbig(:,1,:) = arr_dr0 
      tmpbig(:,2,:) = arr_dr0 

      ! Compute big array
      call tmatxb_shortreal(nrhs, .true., tmpbig, Tarr_dr0)
      call commfieldshort(nrhs,Tarr_dr0)
      call commdirdirshort(nrhs,0,Tarr_dr0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tarr_dr0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tarr_dr0,reqrec,reqsend)

      call diagvec(nrhs, Tarr_dr0, aTarr_dr0)

      ! extra fft needed
      arr_ded = arr_dr0 - aTarr_dr0(:,1,:) + mu0(:,2,:)

      call scalderfieldzmat8(
     &                       arr_ded, arr_dep, arr_dtr0,  r0(:,1,:),
     &                       arr_dTTr0,  Tr0(:,1,:), arr_dTT2r0,
     &                       T2r0(:,1,:), arr_dr0, - oldaefi(:,1,:),
     &                       ade, adme, adte, adtb
     &                      )

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

      return
      end

! TCG2 with fastgrad, fully optimized, all options 
! aboard
      subroutine epolar1tcg2cinqshortreal
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
      integer :: betac, kk,
     $           irhs, k,
     $           ierr,nrhs,kglob,
     $           kkpole,kloc

      integer, allocatable :: reqrec(:),reqsend(:)
      real*8 :: f, sprod, sp1, sp0, b1, b2, a10, a1m1, a11, a12,
     $          a13, a14, a20, a21, a22, a23, a31, a41, a32, ap10a,
     $          ap11a, ap12a, ap11, ap12, ap13, ap14, ap2a0, ap21a,
     $          ap21, ap22, ap23, ap3a0, ap31, ap32, ap41, spp1, spp2,
     $          omega, k1, k2, k3
      real*8, dimension(2) :: n0, t1, t4, np1, t2, t9, t8, t3, t10, g1




      real*8, allocatable, dimension(:,:,:,:) :: adte
      real*8, allocatable, dimension(:,:,:) :: denedt, mu_tcg2,
     $            efi,oldt2r0, oldr0,Tr0, oldTr0,T2r0, mu0,r0,
     $            T3r0,T2aefi, Tefi,
     $            ade,adme,oldaefi,tarr_dr0, atarr_dr0,oldefi,
     $            aefi,Taefi, mupeek,T2efi,tmpbig
      real*8, allocatable, dimension(:,:) :: torq_t, torq_mu, grad_ene,
     $            arr_dr0, arr_dTr0,arr_dTTr0, arr_dTT2r0, denedmu,
     $            denedr, arr_ded,arr_dep, adtb,tmpsmall
      parameter (nrhs=2)

      allocate(adte(3,3,nrhs, npolebloc))
      allocate(denedt(3,3,npoleloc))
      allocate(mu_tcg2(3,nrhs,npolebloc),
     $         efi(3,nrhs,npolebloc),
     $         oldr0(3,nrhs,npolebloc),
     $         oldt2r0(3,nrhs,npolebloc),
     $         Tr0(3,nrhs,npolebloc),
     $         oldTr0(3,nrhs,npolebloc),
     $         T2r0(3,nrhs,npolebloc),
     $         T2aefi(3,nrhs,npolebloc),
     $         mu0(3,nrhs,npolebloc),
     $         r0(3,nrhs,npolebloc),
     $         T3r0(3,nrhs,npolebloc),
     $         arr_dTr0(3,npolebloc),
     $         arr_dTTr0(3,npoleloc),
     $         arr_dTT2r0(3,npolebloc),
     $         ade(3,nrhs,npolebloc), adme(3,nrhs,npolebloc),
     $         oldaefi(3,nrhs,npolebloc),
     $         arr_dr0(3,npolebloc),
     $         Tarr_dr0(3,nrhs,npolebloc),
     $         atarr_dr0(3,nrhs,npolebloc),
     $         oldefi(3,nrhs,npolebloc),
     $         aefi(3,nrhs,npolebloc),
     $         Taefi(3,nrhs,npolebloc),
     $         mupeek(3,nrhs,npolebloc),
     $         T2efi(3,nrhs,npolebloc),
     $         tmpbig(3,nrhs,npolebloc))
      allocate(grad_ene(3,npolebloc), torq_mu(3,nbloc),
     $         torq_t(3,nbloc), denedmu(3,npolebloc),
     $         denedr(3,npolebloc),arr_ded(3,npolebloc),
     $         arr_dep(3,npolebloc),
     $         adtb(3,npolebloc),
     $         tmpsmall(3,npolebloc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      
      f = electric/dielec

      precond = tcgprec
      isguess = tcgguess
      peek    = tcgpeek

      !1. prepare arrays...
      efi = 0d0
      Tr0 = 0d0
      Tarr_dr0 = 0d0

      omega = tcgomega


      !1.1 Electric field
      call efld0_shortreal(nrhs, efi)
      call commfieldshort(nrhs, efi)
      call commdirdirshort(nrhs,0,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,efi,reqrec,reqsend)

      oldefi = efi
      call diagvec(nrhs, efi, efi)
      aefi = efi
      oldaefi = efi
      mu0 = efi

      call tmatxb_shortreal(nrhs, .true., aefi, Taefi)
      call commfieldshort(nrhs,Taefi)
      call commdirdirshort(nrhs,0,Taefi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Taefi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Taefi,reqrec,reqsend)

      Tefi = Taefi

      call diagvec(nrhs, Tefi, Tefi)

      oldr0 = oldefi - Taefi
      call diagvec(nrhs, oldr0, r0)

      call dbletmatxb_shortreal(nrhs, .true.,Taefi,r0, T2aefi, Tr0)
      call commfieldshort(nrhs,T2aefi)
      call commfieldshort(nrhs,Tr0)
      call commdirdirshort(nrhs,0,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tr0,reqrec,reqsend)
      call commdirdirshort(nrhs,0,T2aefi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2aefi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2aefi,reqrec,reqsend)

      oldTr0 = Tr0

      call diagvec(nrhs, Tr0, Tr0)

      call dbletmatxb_shortreal(nrhs, .true., Tr0, Tefi, T2r0, T2efi)
      call commfieldshort(nrhs,T2r0)
      call commfieldshort(nrhs,T2efi)
      call commdirdirshort(nrhs,0,T2efi,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2efi,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2efi,reqrec,reqsend)

      call diagvec(nrhs, T2efi, T2efi)

      call commdirdirshort(nrhs,0,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T2r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T2r0,reqrec,reqsend)

      call diagvec(nrhs, T2r0, T2r0)

      !T3r0
      call tmatxb_shortreal(nrhs, .true., T2r0, T3r0)
      call commfieldshort(nrhs,T3r0)
      call commdirdirshort(nrhs,0,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,T3r0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,T3r0,reqrec,reqsend)

      call diagvec(nrhs, T3r0, T3r0)
c
      do irhs = 1,nrhs
         n0(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), r0(:,irhs,:))
         t1(irhs) = sprod(3*npoleloc, oldr0(:,irhs,:), Tr0(:,irhs,:))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,n0,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,t1,nrhs,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      t4 = n0/t1

      do k = 1,nrhs
         np1(k) = sprod(3*npoleloc, oldTr0(:,k,:), Tr0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,np1(k),1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
         t2(k) = n0(k)*np1(k)/(t1(k)*t1(k))
         t9(k) = sprod(3*npoleloc, oldTr0(:,k,:), T2r0(:,k,:))
         call MPI_ALLREDUCE(MPI_IN_PLACE,t9(k),1,MPI_REAL8,MPI_SUM,
     $      COMM_TINKER,ierr)
         t8(k) = t2(k)*np1(k) - t4(k)*t9(k)
         t3(k) = t1(k)*t8(k)
         t10(k) = t1(k)*t1(k) - n0(k)*np1(k)
         g1(k)  = t10(k)/t3(k)
         mu_tcg2(:,k,:) = mu0(:,k,:) + (g1(k)*t2(k) + t4(k))*r0(:,k,:) 
     $   - g1(k)*t4(k)*Tr0(:,k,:)
      end do
      sp0 = sprod(3*npoleloc, r0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp0,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      sp1 = sprod(3*npoleloc, Tr0(:,1,:), oldefi(:,2,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,sp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      b1 = sp0 - g1(1)*sp1
      b2 = sp0*t2(1) - t4(1)*sp1

      !peek
      mupeek = omega*r0 
     $         - omega*(g1(1)*t2(1) + t4(1))*Tr0 
     $         + omega*g1(1)*t4(1)*t2r0
      mu_tcg2 = mu_tcg2 + mupeek

      ! 2. build derivs arrays...
      ! 2.1 Prepare coefficients
      a10 = t4(1) + g1(1)*t2(1)
c      write(*,*) 'a10 = ',a10
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


      !peek terms
      spp1 = sprod(3*npoleloc, oldefi(:,2,:), Tr0(:,1,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,spp1,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)
      spp2 = sprod(3*npoleloc, oldefi(:,2,:), T2r0(:,1,:))
      call MPI_ALLREDUCE(MPI_IN_PLACE,spp2,1,MPI_REAL8,MPI_SUM,
     $   COMM_TINKER,ierr)

      k1 = omega*t4(1)*spp2 - omega*t2(1)*spp1
      k2 = omega*g1(1)*spp2 - omega*spp1
      k3 = -omega*g1(1)*spp1

      ap10a = omega
      ap11a = -omega*(t2(1)*g1(1) + t4(1))
      ap12a = omega*t4(1)*g1(1)
      ap2a0 = -omega*(g1(1)*t2(1) + t4(1))
      ap21a = omega*t4(1)*g1(1)
      ap3a0 = omega*g1(1)*t4(1)

      ap11 = k1
     $       *(-2d0*np1(1)/t3(1)
     $         - 2d0*np1(1)*np1(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t9(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(2d0/t1(1))
     $      +k3
     $       *(2d0*np1(1)/(t1(1)*t1(1)))
      ap12 = k1
     $       *(4d0*t1(1)/t3(1) 
     $        - 2d0*n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $        + 4d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $        - 2d0*t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(-2d0*n0(1)/(t1(1)*t1(1)))
     $      +k3
     $       *(-4d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap13 = k1
     $       *(-4d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1))
     $         - 2d0*n0(1)/t3(1))
     $      +k3
     $       *(2d0*n0(1)/(t1(1)*t1(1)))
      ap14 = k1
     $       *(2d0*t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap21 = k1
     $       *(2d0*t1(1)/t3(1) 
     $         - n0(1)*t9(1)*t10(1)/(t1(1)*t3(1)*t3(1))
     $         + 2d0*t2(1)*np1(1)*t10(1)/(t3(1)*t3(1))
     $         - t8(1)*t10(1)/(t3(1)*t3(1)))
     $      +k2
     $       *(-n0(1)/(t1(1)*t1(1)))
     $      +k3
     $       *(-2d0*n0(1)*np1(1)/(t1(1)*t1(1)*t1(1)))
      ap22 = k1
     $       *(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3
     $       *(n0(1)/(t1(1)*t1(1)))
      ap23 = k1
     $       *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap31 = k1
     $       *(-n0(1)/t3(1) - 2d0*t1(1)*t2(1)*t10(1)/(t3(1)*t3(1)))
     $      +k3
     $       *(n0(1)/(t1(1)*t1(1)))
      ap32 = k1
     $       *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))
      ap41 = k1
     $        *(t1(1)*t4(1)*t10(1)/(t3(1)*t3(1)))

      ! 2.2 Assemble arrays
      arr_dEp  = mu_tcg2(:,1,:)

      arr_dEd = a10*efi(:,2,:) 
     $            + a1m1*Tefi(:,2,:)
     $            + a11*r0(:,1,:)
     $            + a12*Tr0(:,1,:)
     $            + a13*T2r0(:,1,:) 
     $            + a14*T3r0(:,1,:)

      arr_dTr0 = a20*efi(:,2,:) 
     $            + a21*r0(:,1,:)
     $            + a22*Tr0(:,1,:)
     $            + a23*T2r0(:,1,:)

      arr_dTTr0 = a31*r0(:,1,:)
     $            + a32*Tr0(:,1,:)

      arr_dTT2r0 = a41*r0(:,1,:)

      ! peek terms (precond on ofc)
      arr_ded = arr_ded 
     $          + ap10a*efi(:,2,:)
     $          + ap11a*Tefi(:,2,:)
     $          + ap12a*T2efi(:,2,:)
     $          + ap11*r0(:,1,:)
     $          + ap12*Tr0(:,1,:)
     $          + ap13*T2r0(:,1,:)
     $          + ap14*T3r0(:,1,:)

      arr_dTr0 = arr_dTr0 
     $          + ap2a0*efi(:,2,:)
     $          + ap21a*Tefi(:,2,:)
     $          + ap21*r0(:,1,:)
     $          + ap22*Tr0(:,1,:)
     $          + ap23*T2r0(:,1,:)

      arr_dTTr0 = arr_dTTr0 
     $          + ap3a0*efi(:,2,:)
     $          + ap31*r0(:,1,:)
     $          + ap32*Tr0(:,1,:)

      arr_dTT2r0 = arr_dTT2r0 
     $          + ap41*r0(:,1,:)

      ep = sprod(3*npoleloc, mu_tcg2(:,1,:), oldefi(:,2,:))
      ep = -.5d0*f*ep 

      denedr = 0d0
      denedmu = 0d0
      denedt = 0d0


      !guess option
      ! There is some fundamental optimization to be done here :
      ! this matvec product should be avoidable through some smart
      ! combinations... or maybe not ?
      ! Either we tmatxb all arr_dr0, or we combine already tmatxb'd
      ! vectors composing arr_dr0. But for this second solution,
      ! we'll need tmatxb(T3r0) => a fourth tmatxb

      atarr_dr0 = 0d0


      arr_dr0 = arr_ded
      tmpbig(:,1,:) = arr_dr0
      tmpbig(:,2,:) = arr_dr0

      ! Compute big array
      call tmatxb_shortreal(nrhs, .true., tmpbig, Tarr_dr0)
      call commfieldshort(nrhs,Tarr_dr0)
      call commdirdirshort(nrhs,0,Tarr_dr0,reqrec,reqsend)
      call commdirdirshort(nrhs,1,Tarr_dr0,reqrec,reqsend)
      call commdirdirshort(nrhs,2,Tarr_dr0,reqrec,reqsend)

      call diagvec(nrhs, tarr_dr0, atarr_dr0)
      
      arr_ded = arr_dr0 - atarr_dr0(:,1,:) + mu0(:,2,:)

      call scalderfieldzmat8(
     &                       arr_ded, arr_dep, arr_dtr0,  r0(:,1,:),
     &                       arr_dTTr0,  Tr0(:,1,:), arr_dTT2r0,
     &                       T2r0(:,1,:), arr_dr0, - oldaefi(:,1,:),
     &                       ade, adme, adte, adtb
     &                      )

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

      return
      end

