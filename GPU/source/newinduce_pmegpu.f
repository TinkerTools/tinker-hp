c      
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
c     evaluate induced dipole moments and the polarization energy
c     using either a (preconditioned) conjugate gradient algorithm or
c     Jacobi iterations coupled with DIIS extrapolation.
c
c     literature reference:
c     "Scalable Evaluation of Polarization Energy and Associated Forces
c     in Polarizable Molecular Dynamics: II. Toward Massively Parallel
c     Computations Using Smooth Particle Mesh Ewald",L. Lagardere et al.,
c     J. Chem. Theory Comput., 2015, 11 (6), pp 2589�~@~S2599
c
#include "tinker_precision.h"
      subroutine newinduce_pmegpu
      use atoms  ,only: n,n_pe
      use atmlst
      use domdec
      use ewald
      use inform ,only: deb_Path
      use iounit
      use interfaces,only:inducepcg_pmegpu,tmatxb_p,
     &    efld0_directgpu_p,
     &    inducejac_pmegpu
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use units
      use utilcomm ,buffermpi=>buffermpi2d,buffermpimu=>buffermpimu1
      use utils
      use utilgpu
      use uprior
      use mpi
      use timestat
      implicit none
c
c     with separate cores for reciprocal part
c
      integer i, j, k, nrhs
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
c
c     MPI
c
      integer ierr, iipole, proc,tag
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer status(MPI_STATUS_SIZE)
c
      parameter (nrhs=2)
      real(t_p)  wtime0, wtime1, wtime2, udsum, upsum
      real(t_p)  term, xx(1)
      real(t_p), allocatable :: ef(:,:,:), mu(:,:,:),
     &                       murec(:,:,:)
      real(t_p), allocatable :: cphi(:,:)
c
      if (deb_Path) 
     &   write(*,'(2a,x)') 'newinduce_pmegpu'
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in newinduce: ',F14.5)
c
c     allocate some memory and clear the arrays:
c
      if (rank.le.ndir-1) then
         allocate (mu   (3,nrhs,max(1,npolebloc  )))
         allocate (murec(3,nrhs,max(1,npolerecloc)))
         allocate (ef   (3,nrhs,max(1,npolebloc  )))
         allocate (cphi (10,max(1,npoleloc)))
!$acc enter data create(mu,murec,ef,cphi) async
         j = max(npolerecloc,1)
         call prmem_request(cphirec,10,j,async=.true.)
         j = max(1,npoleloc)
         call prmem_request(buffermpi,10,j,async=.true.)
         call prmem_request(buffermpimu,3,nrhs,j,async=.true.)
         call set_to_zero2(mu,ef       ,3*nrhs*npolebloc,rec_queue)
         call set_to_zero1(murec,size(murec),rec_queue)
         call set_to_zero1(cphi,size(cphi),rec_queue)
         call set_to_zero1(cphirec,size(cphirec),rec_queue)
         call set_to_zero1(buffermpi,size(buffermpi),rec_queue)
         call set_to_zero1(buffermpimu,size(buffermpimu),rec_queue)
      else
         allocate (mu   (3,nrhs,max(1,npolebloc  )))
         allocate (murec(3,nrhs,max(1,npolerecloc)))
         allocate (cphi (10,max(1,npoleloc)))
!$acc enter data create(mu,murec,cphi) async
         j = max(1,npolerecloc)
         if (.not.use_mpole) then
            call prmem_request(cphirec,10,j,async=.true.)
            call prmem_request(fphirec,20,j,async=.true.)
         end if
         call prmem_request(buffermpi,10,j,async=.true.)
         call prmem_request(buffermpimu,3,nrhs,j,async=.true.)
         call set_to_zero1(mu           ,3*nrhs*npolebloc ,rec_queue)
         call set_to_zero2(cphirec,buffermpi,10*j,rec_queue)
         call set_to_zero1(murec        ,size(murec)      ,rec_queue)
         call set_to_zero1(cphi         ,size(cphi)       ,rec_queue)
         call set_to_zero1(fphirec      ,size(fphirec)    ,rec_queue)
         call set_to_zero1(buffermpimu  ,size(buffermpimu),rec_queue)
      end if
c
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c    compute the electric fields:
c
      wtime0 = mpi_wtime()
c
      if (rank.le.ndir-1) then
        call commdirdirgpu(nrhs,0,mu,reqrec,reqsend)
        call commrecdirfieldsgpu(0,cphirec,cphi,buffermpi,buffermpi,
     $       reqrecdirrec,reqrecdirsend)
      else
c
c    compute the reciprocal space contribution (fields)
c
        call timer_enter( timer_recdip )
        call efld0_recipgpu(cphi,ef)
!$acc wait
        call timer_exit ( timer_recdip,quiet_timers )
c
        call commrecdirfieldsgpu(1,cphirec,cphi,buffermpi,buffermpi,
     $       reqrecdirrec,reqrecdirsend)
        call commrecdirfieldsgpu(2,cphirec,cphi,buffermpi,buffermpi,
     $       reqrecdirrec,reqrecdirsend)
!$acc wait
      end if
c
c    The real space processes compute the real fields and  add them to the recip ones
c
      if (rank.le.ndir-1) then
        term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
c
        call timer_enter( timer_realdip )
        call efld0_directgpu_p(nrhs,ef)
        call timer_exit ( timer_realdip,quiet_timers )

        call commfieldgpu(nrhs,ef)
c
!$acc wait
        call commrecdirfieldsgpu(2,cphirec,cphi,buffermpi,buffermpi,
     $       reqrecdirrec,reqrecdirsend)
c
c       Add direct and reciprocal fields
c
        def_queue = rec_queue

!$acc parallel loop collapse(3) async(def_queue)
!$acc&         default(present)
        do i = 1,npoleloc
           do j = 1,nrhs
              do k = 1,3
                 iipole = poleglob(i)
                 ef(k,j,i) = ef(k,j,i) + term*rpole(k+1,iipole)
     &                                 - cphi(k+1,i)
              end do
           end do
        end do
      end if
      wtime1 = mpi_wtime()

c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      if (use_pred .and. nualt.eq.maxualt) then
         call pred_InitField(mu)
      end if

      if (rank.le.ndir-1) then
        if (.not.(use_pred.and.nualt.eq.maxualt)) then

        if (polgsf.eq.0) then
!$acc parallel loop collapse(3) async(def_queue)
!$acc&         default(present)
           do i = 1, npoleloc
             do k = 1, nrhs
               do j = 1, 3
                 iipole = poleglob(i)
                 mu(j,k,i) = polarity(iipole)*ef(j,k,i)
               end do
             end do
           end do
        else
!$acc parallel loop collapse(3) async(def_queue)
!$acc&         default(present)
           do i = 1, npoleloc
              do j = 1,nrhs
                 do k = 1, 3
                    iipole = poleglob(i)
                    if (btest(j,0)) then !j==1
                       mu(k,j,i) = uind(k,iipole)
                    else  !j/=1
                       mu(k,j,i) = uinp(k,iipole)
                    end if
                 end do
              end do
           end do
        end if

        end if
c
        call commdirdirgpu(nrhs,1,mu,reqrec,reqsend)
c
        call commrecdirdipgpu(nrhs,1,murec,mu,buffermpimu,buffermpimu,
     $       req2rec,req2send)
c
        call commdirdirgpu(nrhs,2,mu,reqrec,reqsend)
        call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $       req2rec,req2send)
      else
c
        call commrecdirdipgpu(nrhs,0,murec,mu,buffermpimu,buffermpimu,
     $       req2rec,req2send)
        call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $       req2rec,req2send)
      end if
c
c     now, call the proper solver.
c
      if (polalg.eq.1) then
        if (rank.le.ndir-1) then
          call inducepcg_pmegpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
        else
          call inducepcg_pmegpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
        end if
      else if (polalg.eq.2) then
        if (rank.le.ndir-1) then
          call inducejac_pmegpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
       else
          call inducejac_pmegpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
       end if
      else 
         if (rank.eq.0) write(iout,1000) 
         call fatal
      end if
      wtime2 = mpi_wtime()
      if (polprt.ge.1.and.rank.eq.0) then
        if (polprt.ge.2) then
          write(iout,1010) 'fields:  ', wtime1-wtime0
          write(iout,1010) 'dipoles: ', wtime2-wtime1
        end if
        write(iout,1020) wtime2 - wtime0
      end if
c
c     move the computed dipoles in the common block.
c
      def_queue = rec_queue
      if (rank.le.ndir-1) then
!$acc parallel loop collapse(2) async(def_queue)
!$acc&         default(present)
        do i = 1, npolebloc
          do j = 1, 3
            iipole = poleglob(i)
            uind(j,iipole) = mu(j,1,i)
            uinp(j,iipole) = mu(j,2,i)
          end do
        end do
      else
!$acc parallel loop collapse(2) async(def_queue)
!$acc&         default(present)
        do i = 1, npolerecloc
          do j = 1, 3
            iipole = polerecglob(i)
            uind(j,iipole) = murec(j,1,i)
            uinp(j,iipole) = murec(j,2,i)
          end do
        end do
      end if
c
c     update the lists of previous induced dipole values
c
      if (use_pred) call pred_SetAlt

      if (rank.le.ndir-1) then
!$acc exit data delete(mu,murec,ef,cphi) async
        deallocate (ef)
        deallocate (cphi)
        deallocate (mu)
      else
!$acc exit data delete(mu,murec,cphi) async
        deallocate (murec)
      end if
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqrecdirrec)
      deallocate (reqrecdirsend)

      end
c
      subroutine inducepcg_pmegpu(matvec,nrhs,precnd,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use inform    ,only: deb_Path,abort
      use interfaces,only: tmatxb_pmegpu
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use timestat
      use units
      use utils
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
     &             ,buffermpimu=>buffermpimu1
      use utilgpu
      use mpi
      implicit none
c
c     solves the polarizable dipoles linear equations by preconditioned
c     conjugate gradient. A diagonal preconditioner is used when precnd
c     is true, otherwise the standard conjugate gradient algorithm is
c     recovered by setting the preconditioner to one.
c
      integer  ,intent(in):: nrhs
      logical  ,intent(in)   :: precnd
      real(t_p),intent(in)   :: ef   (:,:,:)
      real(t_p),intent(inout):: mu   (:,:,:)
      real(t_p),intent(inout):: murec(:,:,:)
      procedure(tmatxb_pmegpu)::matvec

      integer i, iglob, it, j, k, iipole
      real(t_p) zero, one
      real(r_p) zerom, pt5
      real(r_p) resnrm
      parameter(zero=0.0,zerom=0,pt5=0.5,one=1.0)
      real(t_p) term
      real(t_p) ggold(2), alphacg(2), gg(2)
      real(r_p) gnorm(2), ggnew(2), ene(2)
      real(t_p),save:: ggold1=0.0,ggold2=0.0
      real(r_p),save:: ggnew1=0.0,ggnew2=0.0
      real(t_p),save:: ggdev1=0.0,ggdev2=0.0
      real(r_p),save:: gnorm1=0.0,gnorm2=0.0
      real(t_p),save:: gg1=0.0,gg2=0.0
      real(r_p),save:: ene1=0.0,ene2=0.0
      real(t_p),save:: alphacg1=0.0,alphacg2=0.0
      real(t_p), allocatable :: res(:,:,:), h(:,:,:),
     $           pp(:,:,:), zr(:,:,:), diag(:)
      real(t_p), allocatable :: dipfield(:,:,:),
     &                       dipfieldbis(:,:,:)
      real(8) time0,time1
      logical,save::first_in=.true.
      logical tinker_isnan_m
c
c     MPI
c
      integer ierr, tag, proc
      integer req1, req2, req3, req4
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer status(MPI_STATUS_SIZE)
c
 1000 format(' cgiter converged after ',I3,' iterations.',/,
     $       ' final energy        = ',2D14.7,/,
     $       ' final residual norm = ',2D14.7)
 1010 format(' energy and residual norm at iteration ',I3,':',4D12.2)
 1020 format(' Conjugate gradient solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Conjugate gradient solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
 1040 format(' Using a diagonal preconditioner.')
 1050 format(' Induce PCG solver not converged after',I6,' steps'
     &    ,/,' Residual norm ', d10.5)
c
c
c     allocate some memory and setup the preconditioner:
c
      call timer_enter ( timer_other )
      if (rank.le.ndir-1) then
         call prmem_request(buffermpimu,3,nrhs,max(npoleloc,1),
     &        async=.true.)
         allocate (dipfield   (3,nrhs,max(1,npoleloc   )))
         allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
         allocate (res(3,nrhs,max(1,npoleloc )))
         allocate (zr (3,nrhs,max(1,npoleloc )))
         allocate (h  (3,nrhs,max(1,npolebloc)))
         allocate (diag(npoleloc))
!$acc enter data create(dipfield,dipfieldbis,res,zr,h,diag) async
         call set_to_zero2(buffermpimu,dipfield,3*nrhs*npoleloc,
     &        rec_queue)
         call set_to_zero1(dipfieldbis,3*nrhs*npolerecloc,
     &        rec_queue)
      else 
         call prmem_request(buffermpimu,3,nrhs,max(npolerecloc,1),
     &        async=.true.)
         allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
         allocate (dipfield   (3,nrhs,max(1,npoleloc)))
         allocate (diag(npoleloc))
!$acc enter data create(diag,dipfield,dipfieldbis) async
         call set_to_zero2(buffermpimu,dipfieldbis,3*nrhs*npolerecloc,
     &        rec_queue)
         call set_to_zero1(dipfield,3*nrhs*npoleloc,rec_queue)
      end if
      allocate (pp(3,nrhs,max(1,npolebloc)))
!$acc enter data create(pp) async
      call prmem_request(buffermpi1,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpi2,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      call set_to_zero1(buffermpi1,3*nrhs*npoleloc,rec_queue)
      call set_to_zero1(buffermpi2,3*nrhs*npolerecloc,rec_queue)

      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqendrec(nproc))
      allocate (reqendsend(nproc))
      allocate (req2endsend(nproc))
      allocate (req2endrec(nproc))
c
      if (deb_Path) 
     &   write(*,'(3x,a)') 'inducepcg_pmegpu'
      if (first_in) then
!$acc enter data create(gg1,gg2) async
         first_in=.false.
      end if
c
      if (precnd) then
!$acc parallel loop present(poleglob,polarity,diag) async
         do i = 1, npoleloc
            iipole = poleglob(i)
            diag(i) = polarity(iipole)
         end do
         if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
!$acc parallel loop present(diag) async
         do i = 1,npoleloc
            diag(i) = 1.0_ti_p
         end do
      end if
c
c     initialize
c
      if (rank.le.ndir-1) then
         call set_to_zero2(pp,h,3*nrhs*npolebloc,rec_queue)
         call set_to_zero1(res,3*nrhs*npoleloc,rec_queue)
      end if
      call timer_exit ( timer_other,quiet_timers )
c
c     now, compute the initial direction
c
      ggold = 0.0_ti_p
      ggold1= 0.0_ti_p; ggold2= 0.0_ti_p;
c
c     MPI : begin reception
c
      if (rank.le.ndir-1) then
        call commdirdirgpu(nrhs,0,pp,reqrec,reqsend)
c
        call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $       buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     The PME processes compute the reciprocal matrix vector product
c
      else 
        call commrecdirdipgpu(nrhs,0,murec,pp,buffermpimu,
     $       buffermpimu,req2rec,req2send)
c
        call timer_enter( timer_recdip )
        call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
        call timer_exit ( timer_recdip,quiet_timers )

        call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $       buffermpi2,reqrecdirrec,reqrecdirsend)
        call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $       buffermpi2,reqrecdirrec,reqrecdirsend)
      end if
      if (rank.le.ndir-1) then
c
         call timer_enter( timer_realdip )
         call matvec(nrhs,.true.,mu,h)
         call timer_exit ( timer_realdip,quiet_timers )
!$acc wait
         call commfieldgpu(nrhs,h)
c        
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $        buffermpi2,reqrecdirrec,reqrecdirsend)
         term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
         call timer_enter( timer_other )
!$acc parallel loop collapse(3) default(present) async
         do k=1,npoleloc
            do j=1,nrhs
               do i=1,3
                  h(i,j,k)   = h(i,j,k)  + dipfield(i,j,k) 
     &                                   - term*mu(i,j,k)
                  res(i,j,k) = ef(i,j,k) - h(i,j,k)
                  zr(i,j,k)  = res(i,j,k)*diag(k)
                  pp(i,j,k)  = zr(i,j,k)
                  if (btest(j,0)) then  ! test de 1
                     ggold1  = ggold1 + res(i,j,k)*zr(i,j,k)
                  else
                     ggold2  = ggold2 + res(i,j,k)*zr(i,j,k)
                  end if
               end do
            end do
         end do
!$acc wait(rec_queue)
         ggold(1)=ggold1; ggold(2)=ggold2;

         call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_TPREC,
     $        MPI_SUM,comm_dir,req1,ierr)
         call timer_exit ( timer_other )
         call commdirdirgpu(nrhs,1,pp,reqrec,reqsend)
c
         call commrecdirdipgpu(nrhs,1,murec,pp,buffermpimu,
     $        buffermpimu,req2rec,req2send)
c
         call commdirdirgpu(nrhs,2,mu,reqrec,reqsend)
c
         call commrecdirdipgpu(nrhs,2,murec,pp,buffermpimu,
     $        buffermpimu,req2rec,req2send)
         call MPI_WAIT(req1,status,ierr)
         ggold1=ggold(1); ggold2=ggold(2);
c
      else 
         call commrecdirdipgpu(nrhs,2,murec,pp,buffermpimu,
     $        buffermpimu,req2rec,req2send)
      end if
c
c     now, start the main loop:
c
      do it = 1, politer
         !gg(1) = zero; gg(2)=zero;
!$acc serial async present(gg1,gg2)
         gg1   = zero; gg2  =zero;
!$acc end serial
         ggnew(1) = zerom; ggnew(2)= zerom;
         ggnew1   = zerom; ggnew2  = zerom;
c
c     MPI : begin reception
c
         if (rank.le.ndir-1) then
            call commdirdirgpu(nrhs,0,pp,reqrec,reqsend)
         else
            call commrecdirdipgpu(nrhs,0,murec,pp,buffermpimu,
     $           buffermpimu,req2rec,req2send)
         end if
c
         if (rank.le.ndir-1) then
            call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield
     $          ,buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
c
         else 
            call timer_enter( timer_recdip )
            call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
            call timer_exit ( timer_recdip,quiet_timers )
            call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield
     $          ,buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
            call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield
     $          ,buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
         end if
         if (rank.le.ndir-1) then
c
            call timer_enter( timer_realdip )
            call matvec(nrhs,.true.,pp,h)
            call timer_exit ( timer_realdip,quiet_timers )
            call commfieldgpu(nrhs,h)
c
            call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield
     $          ,buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
c
            call timer_enter( timer_other )
!$acc parallel loop collapse(3) present(gg1,gg2) async
!$acc&         default(present)
            do k=1,npoleloc
               do j=1,nrhs
                  do i=1,3
                  h(i,j,k) = h(i,j,k) + dipfield(i,j,k) - term*pp(i,j,k)
                     if (btest(j,0)) then
                        gg1 = gg1 + pp(i,j,k)*h(i,j,k)
                     else
                        gg2 = gg2 + pp(i,j,k)*h(i,j,k)
                     end if
                  end do
               end do
            end do

            if (ndir.gt.1) then
!$acc wait
!$acc host_data use_device(gg1,gg2)
            call MPI_ALLREDUCE(MPI_IN_PLACE,gg1,1,MPI_TPREC,MPI_SUM,
     $           comm_dir,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,gg2,1,MPI_TPREC,MPI_SUM,
     $           comm_dir,ierr)
!$acc end host_data
            end if

            ! This check might be unecessary
            ! the only for residu this to reach zero is for the field itself to be equal to zero
            ! TODO Find a proper place for that check
            if (deb_Path) then
!$acc update host(gg1,gg2) async
!$acc wait
               if (gg1.eq.zero) then
                  print*, 'Residu equals zero in polarisation solver'
                  goto 50
               end if
               ! GLITCHY exit only with direct group
               ! TODO Find a way to make reciproqual group to exit routine
               if (gg2.eq.zero) goto 50
            end if

            gg(1)=gg1; gg(2)=gg2;
            ggnew1 = zero; ggnew2 = zero
            ene1 = zerom; ene2 = zerom

!$acc parallel loop collapse(3) async
!$acc&   default(present) present(gg1,gg2)
            do k=1,npoleloc
               do j=1,nrhs
                  do i=1,3
                     if (btest(j,0)) then
                        mu (i,j,k) = mu (i,j,k) + (ggold1/gg1)*pp(i,j,k)
                        res(i,j,k) = res(i,j,k) - (ggold1/gg1)*h (i,j,k)
                        zr (i,j,k) = diag(k)* res(i,j,k)
                        ggnew1     = ggnew1 + res(i,j,k)*zr(i,j,k)
                        ene1 = ene1 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
                      else
                        mu (i,j,k) = mu (i,j,k) + (ggold2/gg2)*pp(i,j,k)
                        res(i,j,k) = res(i,j,k) - (ggold2/gg2)*h (i,j,k)
                        zr (i,j,k) = diag(k)* res(i,j,k)
                        ggnew2     = ggnew2 + res(i,j,k)*zr(i,j,k)
                        ene2 = ene2 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
                     end if
                  end do
               end do
            end do
!$acc wait(rec_queue)
            ene2 = -pt5*ene2; ene1 = -pt5*ene1
            ggnew(1)=ggnew1; ggnew(2)=ggnew2
              ene(1)=  ene1;   ene(2)=  ene2

            call timer_exit ( timer_other,quiet_timers )
         end if

         call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_RPREC,
     $                       MPI_SUM,COMM_TINKER,req3,ierr)
         call MPI_WAIT(req3,status,ierr)

         if (rank.le.ndir-1) then
            call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_RPREC,
     $                          MPI_SUM,comm_dir,req4,ierr)
            call MPI_WAIT(req4,status,ierr)

            call timer_enter ( timer_other )
            ggdev1 = ggnew(1)/ggold(1); ggdev2 = ggnew(2)/ggold(2)
!$acc parallel loop collapse(3) async
!$acc&         default(present)
            do k=1,npoleloc
               do j=1,nrhs
                  do i=1,3
                     if (btest(j,0)) then
                        pp(i,j,k) = zr(i,j,k)+ ggdev1*pp(i,j,k)
                     else
                        pp(i,j,k) = zr(i,j,k)+ ggdev2*pp(i,j,k)
                     end if
                  end do
               end do
            end do
            call timer_exit ( timer_other,quiet_timers )
c
            call commdirdirgpu(nrhs,1,pp,reqrec,reqsend)
c
            call commrecdirdipgpu(nrhs,1,murec,pp,buffermpimu,
     $           buffermpimu,req2rec,req2send)
c
            call commdirdirgpu(nrhs,2,pp,reqrec,reqsend)
            call commrecdirdipgpu(nrhs,2,murec,pp,buffermpimu,
     $           buffermpimu,req2rec,req2send)
         else 
            call commrecdirdipgpu(nrhs,2,murec,pp,buffermpimu,
     $           buffermpimu,req2rec,req2send)
         end if
c        
         ggold = ggnew
         ggold1=ggold(1); ggold2=ggold(2);
         resnrm = zero
         do k = 1, nrhs
            gnorm(k) = sqrt(ggnew(k)/real(3*npolar,r_p))
            resnrm   = max(resnrm,gnorm(k))
         end do
         if ((it.eq.politer.and.resnrm.gt.poleps)
     &      .or.tinker_isnan_m(gnorm(1))) then
            ! Not converged Abort call
            if (rank.eq.0) write(0,1050) it,resnrm
            abort = .true.
         end if
         if (resnrm.lt.poleps) then
            if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     $         (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
           goto 10
         end if
         if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $      it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
      end do
 10   continue
c
      if (rank.le.ndir-1) then
         call commdirdirgpu(nrhs,0,mu,reqendrec,reqendsend)
         call commdirdirgpu(nrhs,1,mu,reqendrec,reqendsend)
         call commrecdirdipgpu(nrhs,1,murec,mu,buffermpimu,
     $        buffermpimu,req2endrec,req2endsend)
c        
         call commdirdirgpu(nrhs,2,mu,reqendrec,reqendsend)
         call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu,
     $        buffermpimu,req2endrec,req2endsend)
      else 
         call commrecdirdipgpu(nrhs,0,murec,mu,buffermpimu,
     $        buffermpimu,req2endrec,req2endsend)
         call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu,
     $        buffermpimu,req2endrec,req2endsend)
      end if
c
c     finalize and debug printing:
c
      if (polprt.ge.3) then
         write(iout,1020)
!$acc wait
!$acc update host(glob,mu)
         do i = 1, npoleloc
            iglob = glob(i)
            write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
         end do
      end if
      if (polprt.ge.4) then
         write(iout,1021)
!$acc wait
!$acc update host(glob,mu)
         do i = 1, npoleloc
            iglob = glob(i)
            write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
         end do
      end if
50    continue
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqendsend)
      deallocate (reqendrec)
      deallocate (req2endsend)
      deallocate (req2endrec)

!$acc exit data delete(diag,pp) async
      deallocate (diag)
      deallocate (pp)
      if (rank.le.ndir-1) then
!$acc exit data delete(res,zr,dipfield,dipfieldbis,h) async
         deallocate (res)
         deallocate (zr)
         deallocate (dipfield)
         deallocate (h)
      else
!$acc exit data delete(dipfield,dipfieldbis) async
         deallocate (dipfieldbis)
      end if

      end
c
      subroutine inducejac_pmegpu(matvec,nrhs,dodiis,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use interfaces ,only: tmatxb_pmegpu
      use math
      use mpole
      use polar
      use polpot
      use potent
      use utils
      use utilgpu
      use timestat
      use mpi
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer  ,intent(in)    :: nrhs
      logical  ,intent(in)    :: dodiis
      real(t_p),intent(in)    :: ef   (:,:,:)
      real(t_p),intent(inout) :: mu   (:,:,:)
      real(t_p),intent(inout) :: murec(:,:,:)
      procedure(tmatxb_pmegpu):: matvec

      integer info
      integer i, j, k, iipole, it, ind, ndismx, nmat, lenb
      parameter (ndismx=25)
      real(t_p), allocatable :: munew(:,:,:), h(:,:,:)
      real(t_p), allocatable :: xdiis(:,:), ediis(:,:),
     $                     bmat(:,:), bloc(:,:), cex(:)
      integer ipiv(ndismx+1)
      real(t_p)  zero, one, rnorm(2), rr, xx(1)
      real(t_p) term
      real(8) time0,time1,time2,time3
      save    zero, one, xx
      data    zero/0.0_ti_p/, one/1.0_ti_p/, xx/0.0_ti_p/
      real(t_p), allocatable :: dipfield(:,:,:),
     &                      dipfieldbis(:,:,:)
c
c     MPI
c
      real(t_p), allocatable :: buffermpi1(:,:,:),
     &                         buffermpi2(:,:,:)
      real(t_p), allocatable :: buffermpimu(:,:,:)
      integer ierr, tag, proc
      integer reqnorm,reqdiis(2*ndismx+1)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer status(MPI_STATUS_SIZE)
c
 1000 format(' itsolv converged after ',I3,' iterations.',/,
     $       ' final residual norm = ',3D14.7)
 1010 format(' residual norm at iteration ',I3,':',3D12.2)
 1020 format(' Jacobi/DIIS solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
c
c
      if (rank.le.ndir-1) then
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0_ti_p
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0_ti_p
        allocate (h(3,nrhs,max(1,npolebloc)))
        h = 0.0_ti_p
        if (dodiis) then
          nmat = 1
          lenb = ndismx + 1
          allocate (xdiis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (ediis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (bmat(lenb,lenb))
          allocate (bloc(lenb,lenb))
          allocate (cex(lenb))
          bmat = 0.0_ti_p
        end if
      else
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0_ti_p
        allocate (buffermpimu(3,nrhs,max(1,npolerecloc)))
        buffermpimu = 0.0_ti_p
      end if
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0_ti_p
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0_ti_p
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0_ti_p
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0_ti_p
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
c     main loop:
c
      do it = 1, politer
        rnorm = 0.0_ti_p
c
c     MPI : begin reception
c
        if (rank.le.ndir-1) then
          call commdirdir(nrhs,0,mu,reqrec,reqsend)
        else 
          call commrecdirdip(nrhs,0,dipfieldbis,dipfield,
     $      buffermpimu,buffermpimu,req2rec,req2send)
        end if
c
        if (rank.le.ndir-1) then
           call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $     buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     The PME processes compute the reciprocal matrix vector product
c
c        else
        else 
          time0 = mpi_wtime()
          call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
          time1 = mpi_wtime()
          if (it.eq.1) timerecdip = timerecdip + time1-time0
          call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
        end if
c
        if (rank.le.ndir-1) then
c
c     jacobi step:
c
c    The real space processes extract the recip fields, compute the real fields
c    and add them
c
         term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
         time0 = mpi_wtime()
         call matvec(nrhs,.false.,mu,h)
         time1 = mpi_wtime()
         if (it.eq.1) timerealdip = timerealdip + time1-time0
         call commfield(nrhs,h)
c
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
c
          h(1:3,1:nrhs,1:npoleloc) = h(:,:,1:npoleloc)+
     $      dipfield(:,:,1:npoleloc)
     $     - term*mu(:,:,1:npoleloc)
c
          do i = 1, npoleloc
            iipole = poleglob(i)
            do k = 1, nrhs
              do j = 1, 3
                munew(j,k,i) = polarity(iipole)*(ef(j,k,i) - h(j,k,i))
              end do
            end do
          end do
c
          do k = 1, nrhs
            rnorm(k)=sum((munew(1:3,k,1:npoleloc)-
     $        mu(1:3,k,1:npoleloc))**2)
          end do
c
        end if
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_TPREC,
     $    MPI_SUM,COMM_TINKER,reqnorm,ierr)
        if (rank.le.ndir-1) then
          if (dodiis) then
            ind = 0
            do i = 1, npoleloc
              do k = 1, nrhs
                do j = 1, 3
                  ind = ind + 1
                  xdiis(ind,nmat) = munew(j,k,i)
                  ediis(ind,nmat) = munew(j,k,i) - mu(j,k,i)
                end do
              end do
            end do
c
c           Compute Pulay's Matrix and extrapolate
c
            call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $        bmat,nmat,reqdiis,comm_dir)
            do i = 1, 2*nmat-3
              call MPI_WAIT(reqdiis(i),status,ierr)
            end do
            bloc = bmat
            cex = 0.0_ti_p
            cex(1) = one
#ifdef _OPENACC
            call fatal_device("inducejac_pmegpu")
#else
            call M_gesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
            munew = 0.0_ti_p
            call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
          end if
        time2 = mpi_wtime()
c
          call commdirdir(nrhs,1,munew,reqrec,reqsend)
          call commrecdirdip(nrhs,1,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
c
          time0 = mpi_wtime()
          call commdirdir(nrhs,2,mu,reqrec,reqsend)
          call commrecdirdip(nrhs,2,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
          mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)
        else 
          call commrecdirdip(nrhs,2,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
        end if
c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/real(3*npolar,t_p))
c          rnorm(k) = sqrt(rnorm(k))
          rr = max(rnorm(k),rr)
        end do
        if (polprt.ge.2.and.rank.eq.0)
     $         write(6,1010) it, (rnorm(k), k = 1, nrhs)
        if (rr.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1000) it, (rnorm(k), k = 1, nrhs)
          goto 10
        end if
      end do
 10   continue
c
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iipole = poleglob(i)
          write(iout,1030) iipole, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iipole = poleglob(i)
          write(iout,1030) iipole, (mu(j,2,i), j = 1, 3)
        end do
      end if
c
c     free the memory.
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      deallocate (buffermpimu)
      if (rank.le.ndir-1) then
        deallocate (munew)
        deallocate (h)
        deallocate (dipfield)
        if (dodiis) then
          deallocate (xdiis)
          deallocate (ediis)
          deallocate (bmat)
          deallocate (bloc)
          deallocate (cex)
        end if
      else
        deallocate (dipfieldbis)
      end if
      return
      end

      subroutine mallocMpiGrid
      use domdec
      use fft
      !use pme
      use pme1
      implicit none
      integer i

      if (.not.allocated(qgridmpi).and.nrec_recep.ge.1) then
         allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
!$acc enter data create(qgridmpi)
      end if
      end subroutine
c
c     Compute the reciprocal space contribution to the electric field due to the permanent 
c     multipoles
c
      subroutine efld0_recipgpu(cphi,ef)
      use atmlst
      use bound
      use boxes
      use domdec
      use ewald
      use fft
      use inform     ,only: deb_Path
      use interfaces ,only: fphi_mpole_site_p
     &               ,grid_mpole_site_p,efld0_direct_cp
      use polar_temp ,only: cmp=>fphid,fmp=>fphip !Reuse module memory
      use math
      use mpole
      use pme
      use pme1
      use potent
      use timestat
      use utils
      use utilgpu
      use mpi
      implicit none
      integer ierr,gang_
      integer status(MPI_STATUS_SIZE),tag
      integer i,j,k,iglob,iipole,iloc
      integer k1,k2,k3,b1,b2,b3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer istat2,ied2,jstat2,jed2,kstat2,ked2
      integer qgridout_size
      integer qgridin_size
      integer,parameter::nrhs=2
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) cphi(10,max(npoleloc,1))
      real(t_p) ef(3,nrhs,npolebloc)
c     real(t_p),dimension(10,npolerecloc)::cmp,fmp
      integer, allocatable :: reqrec(:),reqsend(:)
      integer, allocatable :: reqbcastrec(:),reqbcastsend(:)
      integer nprocloc,commloc,rankloc,proc
      real(8) time0,time1,timem
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return

      if (deb_Path)
     &   write(*,'(3x,a)') 'efld0_recipgpu'
      call timer_enter( timer_efld0_recip )

      if (use_pmecore) then
         nprocloc = nrec
         commloc  = comm_rec
         rankloc  = rank_bis
      else
         nprocloc = nproc
         commloc  = COMM_TINKER
         rankloc  = rank
      end if

      allocate (reqrec(nprocloc))
      allocate (reqsend(nprocloc))
      allocate (reqbcastrec(nprocloc))
      allocate (reqbcastsend(nprocloc))

      call mallocMpiGrid
      call prmem_request(cmp,10,npolerecloc,async=.true.)
      call prmem_request(fmp,10,npolerecloc,async=.true.)

      istat2 = istart2(rankloc+1)
      ied2   = iend2  (rankloc+1)
      jstat2 = jstart2(rankloc+1)
      jed2   = jend2  (rankloc+1)
      kstat2 = kstart2(rankloc+1)
      ked2   = kend2  (rankloc+1)
      qgridin_size = 2*(nrec_send+1)*n1mpimax*n2mpimax*n3mpimax

!!$acc enter data create(cmp,fmp) async(rec_queue)
c
!!$acc data
!!$acc&     present(cphi)
!!$acc&     present(polerecglob,ipole,rpole,kstart2,kend2,
!!$acc&  jstart2,jend2,istart2,iend2,bsmod1,bsmod2,bsmod3,
!!$acc&  repart,qgridin_2d,qgridout_2d,use_bounds,qfac_2d,
!!$acc&  octahedron,use_pmecore,cphirec,fphirec,recip,cmp,fmp)
!!$acc&     async(rec_queue)
c
c     zero out the PME charge grid
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgridin_2d(1,1,1,1,1),
     &                  qgridin_size,rec_queue)
c
c     fill the pme grid, loop over the multipoles sites
c
!$acc parallel loop collapse(2) async(rec_queue)
!$acc&         present(polerecglob,rpole)
      do i = 1, npolerecloc
         do j = 1, 10
            iipole   = polerecglob(i)
            cmp(j,i) = rpole_scale(j)*rpole(rpole_ind_extract(j),iipole)
         end do
      end do
c
c     compute B-spline coefficients
c
      !call bspline_fill_sitegpu
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp_sitegpu(cmp,fmp)
c
c     assign PME grid
c
      call grid_mpole_site_p(fmp)
      call timer_exit ( timer_grid1,quiet_timers )
c
      !Exchange Grid Between reciprocal process
      call commGridFront( qgridin_2d,r_comm )
c
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         call start_dir_stream_cover
         call efld0_direct_cp(nrhs,ef)
      end if
#endif
c
      ! Wait for Grid communication
      call commGridFront( qgridin_2d,r_wait )
c
c     Perform 3-D FFT forward transform
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#else
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#endif
c
c     make the scalar summation over reciprocal lattice
c
      call timer_enter( timer_scalar )

      pterm   = (pi/aewald)**2
      volterm = pi * volbox
      nf1     = (nfft1+1) / 2
      nf2     = (nfft2+1) / 2
      nf3     = (nfft3+1) / 2
      qgridout_size= 2*ksize2(rankloc+1)*jsize2(rankloc+1)
     &                *isize2(rankloc+1)

!$acc serial async(rec_queue) present(qfac_2d)
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     &   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0_ti_p
      end if
!$acc end serial

!$acc parallel loop collapse(3) async(rec_queue)
!$acc&         present(qfac_2d,recip,octahedron,use_bounds)
      do k3 = kstat2,ked2
        do k2 = jstat2,jed2
          do k1 = istat2,ied2
             m1 = k1 - 1
             m2 = k2 - 1
             m3 = k3 - 1
             if (k1 .gt. nf1)  m1 = m1 - nfft1
             if (k2 .gt. nf2)  m2 = m2 - nfft2
             if (k3 .gt. nf3)  m3 = m3 - nfft3
             if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
             r1 = real(m1,t_p)
             r2 = real(m2,t_p)
             r3 = real(m3,t_p)
             h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
             h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
             h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
             hsq = h1*h1 + h2*h2 + h3*h3
             term = -pterm * hsq
             expterm = 0.0_ti_p
             if ((term .gt. -50.0_ti_p)) then
                denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
                expterm = exp(term) / denom
                if (.not. use_bounds) then
                   expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
                else if (octahedron) then
                   if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
                end if
                qfac_2d(k1-istat2+1,k2-jstat2+1,k3-kstat2+1) = expterm
             end if
 10          continue
          end do
        end do
      end do
c
c     account for the zeroth grid point for a finite system
c
!$acc serial async(rec_queue) default(present)
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     &   (kstart2(rankloc+1).eq.1)) then
         if (.not. use_bounds) then
            expterm = 0.5_ti_p * pi / xbox
            qfac_2d(1,1,1) = expterm
         end if
      end if
!$acc end serial
c
c     complete the transformation of the charge grid
c
      call convolution_product(qfac_2d,qgridout_size,qgridout_2d
     &                         ,rec_queue)
      call timer_exit ( timer_scalar,quiet_timers )
c
c     perform 3-D FFT backward transform
c
#ifdef _OPENACC
      call cufft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#else
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#endif
c
      call commGridBack( qgridin_2d,r_comm )
      ! Wait
      call commGridBack( qgridin_2d,r_wait )

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
c
c     get field
c
      call timer_enter( timer_grid2 )
      call fphi_mpole_site_p
      call fphi_to_cphi_sitegpu(fphirec,cphirec)

!$acc parallel loop async(rec_queue)
!$acc&         default(present)
      do i = 1, npolerecloc
         do j = 1, 10
         iipole = polerecglob(i)
         iglob  = ipole(iipole)
         if (.not.(use_pmecore).and.(repart(iglob).eq.rankloc)) then
            iloc         = poleloc(iipole)
            cphi(j,iloc) = cphirec(j,i)
         end if
         end do
      end do
      call timer_exit ( timer_grid2,quiet_timers )

!!$acc end data
c
!!$acc exit data delete(cmp,fmp) async(rec_queue)
      deallocate (reqbcastrec)
      deallocate (reqbcastsend)
      deallocate (reqrec)
      deallocate (reqsend)
      
      call timer_exit( timer_efld0_recip )
      end


c
c     Compute the reciprocal space contribution to the electric field due to the current  
c     value of the induced dipoles
c
      subroutine tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
      use atmlst
      use boxes
      use domdec
      use ewald
      use fft
      use math
      use mpole
      use inform    ,only: deb_Path
      use interfaces,only: fphi_uind_site2_p
     &              ,grid_uind_site_p,tmatxb_cp
      use pme
      use pme1
      use potent
      use polar_temp,only: fuind,fuinp,h
     &              ,fdip_phi1=>fphid,fdip_phi2=>fphip !Reuse module memory
      use timestat
      use mpi
      use utils
      use utilgpu,only: rec_queue,dir_queue,prmem_request
#ifdef _OPENACC
     &           ,start_dir_stream_cover,end_dir_stream_cover
#endif
      implicit none
      integer ierr,iglob
      integer status(MPI_STATUS_SIZE),tag
      integer nrhs,iipole
      integer i,j,k,i1,i2,iloc
      integer qgrid2_size
      real(t_p) time0,time1,timem
c     real(t_p) fuind(3), fuinp(3)
      real(t_p) term
      real(t_p),parameter ::zero=0.0_ti_p
      real(t_p),save:: a(3,3)
c     real(t_p) fdip_phi1(10), fdip_phi2(10), fdip_sum_phi(20)
      real(t_p),intent(inout):: dipfield(3,nrhs,npoleloc)
     &         ,dipfieldbis(3,nrhs,max(1,npolerecloc))
      real(t_p),intent(inout):: mu(3,nrhs,max(1,npolebloc))
     &         ,murec(3,nrhs,max(1,npolerecloc))
c     real(t_p),dimension(3,npolerecloc)::fuind,fuinp
c     real(t_p),dimension(10,npolerecloc)::fdip_phi1,fdip_phi2

      integer, allocatable :: reqbcastrec(:),reqbcastsend(:)
      integer, allocatable :: reqrec(:),reqsend(:)
      integer nprocloc,commloc,rankloc,proc
      logical,save::f_in=.true.

      interface
      subroutine fftTrC_tmat_alloc(nrhs,npolebloc
     &                      ,dodiag,argmu,argef)
      integer,intent(in)::nrhs,npolebloc
      logical,intent(in)::dodiag
      real(t_p),dimension(3,nrhs,npolebloc),target::argmu,argef
      end subroutine
      end interface
c
      if (deb_Path)
     &   write(*,'(4x,a)') 'tmatxbrecipgpu'
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return

      if (use_pmecore) then
         nprocloc = nrec
         commloc  = comm_rec
         rankloc  = rank_bis
      else
         nprocloc = nproc
         commloc  = COMM_TINKER
         rankloc  = rank
      end if
      call timer_enter( timer_tmatxrecipgpu )

      call mallocMpiGrid
      call prmem_request(fuind    , 3,npolerecloc,async=.false.)
      call prmem_request(fuinp    , 3,npolerecloc,async=.false.)
      call prmem_request(fdip_phi1,10,npolerecloc,async=.false.)
      call prmem_request(fdip_phi2,10,npolerecloc,async=.false.)

      allocate (reqbcastrec(nprocloc))
      allocate (reqbcastsend(nprocloc))
      allocate (reqrec(nprocloc))
      allocate (reqsend(nprocloc))

      if (f_in) then
!$acc enter data create(a)
      end if

!$acc data present(mu,murec)
!$acc&     present(fuind,fuinp,fdip_phi1,fdip_phi2
!$acc&    ,polerecglob,ipole,poleloc,repart
!$acc&    ,recip,a,ksize2,jsize2,isize2)
c
c     zero out the PME charge grid
c
      call timer_enter( timer_grid1 )
      call set_to_zero1(qgrid2in_2d(1,1,1,1,1),
     &     2*n1mpimax*n2mpimax*n3mpimax*(nrec_send+1),rec_queue)
c
c     fill the pme grid, loop over the induced dipoles sites
c
!$acc parallel loop collapse(2) async(rec_queue)
      do j = 1, 3
         do i = 1, 3
            if (i.eq.1) then
               a(1,j) = real(nfft1,t_p) * recip(j,1)
            else if (i.eq.2) then
               a(2,j) = real(nfft2,t_p) * recip(j,2)
            else
               a(3,j) = real(nfft3,t_p) * recip(j,3)
            end if
         end do
      end do
c
!$acc parallel loop collapse(2)
!$acc&         async(rec_queue)
      do i = 1, npolerecloc
         do k = 1, 3
            iipole = polerecglob(i)
            iglob  = ipole(iipole)
            iloc   = poleloc(iipole)
c         
c           Convert cartesian dipoles to fractional coordinates
c         
            if (repart(iglob).ne.rank) then
               fuind(k,i) = a(k,1)*murec(1,1,i)
     &                    + a(k,2)*murec(2,1,i)
     &                    + a(k,3)*murec(3,1,i)
               fuinp(k,i) = a(k,1)*murec(1,2,i) 
     &                    + a(k,2)*murec(2,2,i) 
     &                    + a(k,3)*murec(3,2,i)
            else
               fuind(k,i) = a(k,1)*mu(1,1,iloc) 
     &                    + a(k,2)*mu(2,1,iloc)
     &                    + a(k,3)*mu(3,1,iloc)
               fuinp(k,i) = a(k,1)*mu(1,2,iloc) 
     &                    + a(k,2)*mu(2,2,iloc)
     &                    + a(k,3)*mu(3,2,iloc)
            end if
         end do
      end do
c
c     assign PME grid
c
      call grid_uind_site_p(fuind,fuinp,qgrid2in_2d)
      call timer_exit ( timer_grid1,quiet_timers )
c
      !Exchange Grid Between reciprocal process
      call commGridFront( qgrid2in_2d,r_comm )
c
      ! Call to matvec to cover MPI communication
      ! Only enabled in parallel
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         call fftTrC_tmat_alloc(nrhs,npolebloc,.true.,mu,h)
c        call start_dir_stream_cover
c        call tmatxb_cp(nrhs,.true.,mu,h)
c        call incMatvecStep
      end if
#endif
c
      !Wait for commGrid and add accumulate result
      call commGridFront( qgrid2in_2d,r_wait )
c
#ifdef _OPENACC
      ! set an event to sync with rec_stream
      !if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
c
c     Perform 3-D FFT forward transform
c
#ifdef _OPENACC
      call cufft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#else
      call fft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#endif

      qgrid2_size = 2*isize2(rankloc+1)*jsize2(rankloc+1)
     &               *ksize2(rankloc+1)
#ifdef _OPENACC
      ! Increment Matrix-Vector product counter
      if (dir_queue.ne.rec_queue) call incMatvecStep
#endif
c
c     complete the transformation of the charge grid
c
      call timer_enter( timer_scalar )
      call convolution_product(qfac_2d,qgrid2_size,qgrid2out_2d
     &                        ,rec_queue)
      call timer_exit ( timer_scalar,quiet_timers )
c
c     perform 3-D FFT backward transform
c
#ifdef _OPENACC
      call cufft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#else
      call fft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     &     n3mpimax)
#endif
c
      ! Exchange Grid Between reciprocal process
      call commGridBack( qgrid2in_2d,r_comm )
c
      ! Compute second half of matvec
#ifdef _OPENACC
c     if (dir_queue.ne.rec_queue) then
c        call start_dir_stream_cover
c        call tmatxb_cp(nrhs,.true.,mu,h)
c        call resetMatvecStep
c     end if
#endif
c
      ! Wait for Backward Grid communication
      call commGridBack( qgrid2in_2d,r_wait )

#ifdef _OPENACC
      ! set an event to sync with rec_stream
c     if (dir_queue.ne.rec_queue) call end_dir_stream_cover
      if (dir_queue.ne.rec_queue) then
         call resetMatvecStep
         call fftTrC_tmat_free
      end if
#endif

c
c     get fields
c
      call timer_enter( timer_grid2 )
      call fphi_uind_site2_p(fdip_phi1,fdip_phi2)

!$acc parallel loop collapse(2)
!$acc&         async(rec_queue)
      do i = 1, npolerecloc
         do k = 1, 3
            iipole = polerecglob(i)
            iglob  = ipole(iipole)
            iloc   = poleloc(iipole)
            if (.not.(use_pmecore).and.(repart(iglob).eq.rankloc)) then
c
c     convert the dipole fields from fractional to Cartesian
c
               dipfield(k,1,iloc) = a(k,1)*fdip_phi1(2,i)
     &                            + a(k,2)*fdip_phi1(3,i)
     &                            + a(k,3)*fdip_phi1(4,i)
               dipfield(k,2,iloc) = a(k,1)*fdip_phi2(2,i)
     &                            + a(k,2)*fdip_phi2(3,i)
     &                            + a(k,3)*fdip_phi2(4,i)
            else
               dipfieldbis(k,1,i) = a(k,1)*fdip_phi1(2,i)
     &                            + a(k,2)*fdip_phi1(3,i)
     &                            + a(k,3)*fdip_phi1(4,i)
               dipfieldbis(k,2,i) = a(k,1)*fdip_phi2(2,i)
     &                            + a(k,2)*fdip_phi2(3,i)
     &                            + a(k,3)*fdip_phi2(4,i)
            end if
         end do
      end do
      call timer_exit ( timer_grid2,quiet_timers )

!$acc end data
c
      if (f_in) f_in=.false.
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (reqbcastsend)
      deallocate (reqbcastrec)
      
      call timer_exit( timer_tmatxrecipgpu )
      end

      subroutine print_arrayr(array,n,col,name)
      use tinheader
      implicit none
      real(t_p) array(*)
      character*12 name
      integer:: n,col,i=1

      print*,name
c     do while (i.le.n)
c        if (i+col-1<=n) then
c           print*,array(i:i+col-1)
c           i = i+col
c        else
c           print*,array(i:n)
c           i = n+1
c        end if
c     end do
      do i=1, n
         if (array(i).ne.0.0_ti_p)
     &      print*, i,array(i)
      end do
      print*,''
      end subroutine

      subroutine print_arraym(array,array2,array3,n,name)
      use tinheader
      implicit none
      real(t_p) array(*),array2(*),array3(*)
      character*12 name
      integer:: n,col,i=1

      print*,name
      do i=1, n
         print*, i
         print*, array(i),array2(i),array3(i)
      end do
      print*,''
      end subroutine

      subroutine print_arrayi(array,n,col,name)
      implicit none
      integer array(*)
      character*12 name
      integer:: n,col,i=1

      print*,name
      do while (i.le.n)
         if (i+col-1<=n) then
            print*,array(i:i+col-1)
            i = i+col
         else
            print*,array(i:n)
            i = n+1
         end if
      end do
      print*,''
      end subroutine
