c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################################################
c     ##                                                                                           ##
c     ##  subroutine baoab4sitepistonqtb  --  QTB-BAOAB/Langevin Piston NPT Langevin molecular dynamics step    ##
c     ##                                                                                           ##
c     ###############################################################################################
c
c
c     "baoabpiston" performs a single molecular dynamics time step
c     via the BAOAB recursion formula and using a Langevin Piston barostat
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c     
c      Constant pressure molecular dynamics simulation: The Langevin piston method
c      J. Chem. Phys. 103, 4613 (1995)
c      Scott E. Feller, Yuhong Zhang,Richard W. Pastor, Bernie Brooks
c
      subroutine baoab4sitepistonqtb (istep,dt)
      use adqtb
      use atmtyp
      use atoms
      use bath
      use cutoff
      use deriv
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use timestat
      use units
      use usage
      use neigh
      use mpi
      use qtb
      implicit none
      integer i,j,k,l,istep,iglob,ierr,iloc
      integer iloc1,iloc2,iloc3,ilocrec1,ilocrec2,ilocrec3
      real*8 dt,dt_2,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 oterm,hterm
      real*8 a1,a2,normal
      real*8 a1piston,a2piston,R
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
      real*8 scale,velfact
    
      real*8, allocatable :: mCvv_piston(:), Cvf_piston(:)
      real*8, allocatable :: dFDR_piston(:)
      double complex, allocatable :: s_in(:), s_out_v(:), s_out_r(:)
      integer*8 planv, planr,est,s
      character*120 Cvf_file
      character*120 mCvv_file

      INTERFACE
         function int_to_str(value) result(str)
           IMPLICIT NONE
           INTEGER, INTENT(in) :: value
           CHARACTER(:), ALLOCATABLE :: str
         end
      END INTERFACE
c
c     set time values and coefficients for BAOAB integration
c
      dt_2 = 0.5d0*dt
c      a1 = exp(-gamma*dt)
c      a2 = sqrt((1-a1**2)*boltzmann*kelvin)
      a1piston = exp(-gammapiston*dt)
      a2piston = sqrt((1-a1piston**2)*boltzmann*kelvin/masspiston)
c
c     set linear combination factors for TIP4P or Dang-Chang water
c
      oterm = 0.73612d0
      hterm = 0.13194d0
c
c     find quarter step velocities and half step positions via BAOAB recursion
c
      vextvol = vextvol + dt_2*aextvol
      
      do i = 1, nloc
         iglob = glob(i)
         if (atomic(iglob).eq.0) then
           v(:,iglob) = 0d0
         else if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + dt_2*a(j,iglob)
            end do
         end if
      end do
      
c
      if (use_rattle) call rattle2(dt)
c
      
      extvolold = extvol
      scale=exp(dt_2*vextvol) 
      velfact=sinh(dt_2*vextvol)/vextvol
      extvol = extvol*(scale**3) ! + vextvol*0.5*dt
      call rescale(istep)
      do i = 1, nloc
        iglob = glob(i)
        if (atomic(iglob).eq.0) then
           v(:,iglob) = 0d0
           xold(iglob) = x(iglob)
           yold(iglob) = y(iglob)
           zold(iglob) = z(iglob)
        else if (use(iglob)) then
          xold(iglob) = x(iglob)
          yold(iglob) = y(iglob)
          zold(iglob) = z(iglob)          
          x(iglob) = x(iglob)*scale + v(1,iglob)*velfact
          y(iglob) = y(iglob)*scale + v(2,iglob)*velfact
          z(iglob) = z(iglob)*scale + v(3,iglob)*velfact
          v(:,iglob)=v(:,iglob)/scale          
        end if
      end do
c
      if (use_rattle) call rattle(dt_2,xold,yold,zold)
      if (use_rattle) call rattle2(dt_2)
c
c      if (use_rattle) call rattle2(dt)
c
c     compute random part
c

      call qtbrandom(dt,istep)

    
      if (rank.eq.0) then
        R = normal()
        vextvol = a1piston*vextvol + a2piston*R
      end if
      if((rank.eq.0).and.(adaptive)) then
        k = mod(istep-1,nseg)+1
        vad_piston(k)=sqrt(a1piston)*vextvol+(a2piston*R/
     &              sqrt(masspiston))/2
        fad_piston(k)=a2piston*R*sqrt(masspiston)/dt
      endif


      call MPI_BCAST(vextvol,1,MPI_REAL8,0,COMM_BEAD,ierr)
c
      if (use_rattle) call rattle2(dt)
c
c     find full step positions via BAOAB recursion
c
      
      extvolold = extvol
      scale=exp(dt_2*vextvol) 
      velfact=sinh(dt_2*vextvol)/vextvol
      extvol = extvol*(scale**3) ! + vextvol*0.5*dt
      call rescale(istep)
      do i = 1, nloc
        iglob = glob(i)
        if (atomic(iglob).eq.0) then
           v(:,iglob) = 0d0
           xold(iglob) = x(iglob)
           yold(iglob) = y(iglob)
           zold(iglob) = z(iglob)
        else if (use(iglob)) then
          xold(iglob) = x(iglob)
          yold(iglob) = y(iglob)
          zold(iglob) = z(iglob)          
          x(iglob) = x(iglob)*scale + v(1,iglob)*velfact
          y(iglob) = y(iglob)*scale + v(2,iglob)*velfact
          z(iglob) = z(iglob)*scale + v(3,iglob)*velfact
          v(:,iglob)=v(:,iglob)/scale          
        end if
      end do
      
c
      if (use_rattle) call rattle(dt_2,xold,yold,zold)
      if (use_rattle) call rattle2(dt_2)

c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
      call reassignqtb(istep)
      
c
c     -> reciprocal space
c
      call reassignpme(.false.)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
c
c     communicate positions
c
      time0 = mpi_wtime()
      call commpos
      call commposrec
      time1 = mpi_wtime()
      timecommstep = timecommstep + time1 - time0
c
      allocate (derivs(3,nbloc))
      derivs = 0d0
c
c
c     move M center to linear combination of H-O-H coordinates
c
c      do i = 1, nlocnl
c         iglob = ineignl(i)
      do i = 1, nbloc
         iglob = glob(i)
         if (atomic(iglob) .eq. 0) then
            x(iglob) = oterm*x(iglob-3) + hterm*x(iglob-2) + 
     $           hterm*x(iglob-1)
            y(iglob) = oterm*y(iglob-3) + hterm*y(iglob-2) +
     $           hterm*y(iglob-1)
            z(iglob) = oterm*z(iglob-3) + hterm*z(iglob-2) + 
     $           hterm*z(iglob-1)
         end if
      end do
      call reinitnl(istep)
c
      time0 = mpi_wtime()
      call mechanicstep(istep)
      time1 = mpi_wtime()
c
      timeparam = timeparam + time1 - time0
c
      time0 = mpi_wtime()
      call allocstep
      time1 = mpi_wtime()
      timeclear = timeclear  + time1 - time0
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     transfer the forces on the M center to the H-O-H sites
c
      do i = 1, nlocnl
         iglob = ineignl(i)
         iloc = loc(iglob)
         if (atomic(iglob) .eq. 0) then
            iloc1 = loc(iglob-1)
            iloc2 = loc(iglob-2)
            iloc3 = loc(iglob-3)
            do j = 1, 3
               derivs(j,iloc3) = derivs(j,iloc3) + oterm*derivs(j,iloc)
               derivs(j,iloc2) = derivs(j,iloc2) + hterm*derivs(j,iloc)
               derivs(j,iloc1) = derivs(j,iloc1) + hterm*derivs(j,iloc)
            end do
         end if
      end do
      do i = 1, nlocrec
         iglob = globrec(i)
         iloc = locrec1(iglob)
         if (atomic(iglob) .eq. 0) then
            iloc1 = locrec1(iglob-1)
            iloc2 = locrec1(iglob-2)
            iloc3 = locrec1(iglob-3)
            do j = 1, 3
               decrec(j,iloc3) = decrec(j,iloc3) + oterm*decrec(j,iloc)
               decrec(j,iloc2) = decrec(j,iloc2) + hterm*decrec(j,iloc)
               decrec(j,iloc1) = decrec(j,iloc1) + hterm*decrec(j,iloc)
            end do
         end if
      end do
c
c     MPI : get total energy
c
      call reduceen(epot) 
c      if ((mod(istep,100).eq.0).and.(rank.eq.0)) then 
c        write(*,*) 'ec=',ec 
c        write(*,*) 'epot=', epot
c        write(*,*) 'esum=', esum
c        write(*,*) 'ea=', ea
c        write(*,*) 'eb=', eb
c        write(*,*) 'ev=', ev
c      endif
c
c     communicate forces
c
      call commforces(derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (atomic(iglob) .eq. 0) then
            do j = 1, 3
               v(j,iglob) = 0.0d0
               a(j,iglob) = 0.0d0
            end do
         else if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = v(j,iglob) + dt_2*a(j,iglob)
            end do
         end if
      end do
      
c
c
c     temperature and pressure estimation
c
      call temper (dt,eksum,ekin,temp)
c      call pressure2 (epot,temp)
      call pressure (dt,epot,ekin,temp,pres,stress,istep)
      
      aextvol = 3.d0*convert*(extvol*(pres-atmsph)/prescon 
     &  +gasconst*kelvin)/masspiston
      vextvol = vextvol + dt_2*aextvol
      
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c      call temper (dt,eksum,ekin,temp)
c      call pressure (dt,epot,ekin,temp,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + esum
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
c      call mdrest (istep)
      if ((mod(istep,nseg).eq.0).and.adaptive) then
                if (ir) then
                  call irspectra
                endif
                  compteur=compteur+1
                  write(*,*) compteur
            if(compteur .ge. skipseg) then
              call fftw(istep,dt)
              call adHnoise(dt)


              if(rank.eq.0) then
                allocate(s_in(3*nseg))
                allocate(s_out_v(3*nseg))
                allocate(s_out_r(3*nseg))

                allocate(mCvv_piston(0:nad-1))
                allocate(Cvf_piston(0:nad-1))
c               allocate(dFDR_piston(0:nad-1))

                s=3*nseg
                est=1
    
                do i=1,nad-1
                  mCvv_piston(i)=0d0
                  Cvf_piston(i)=0d0
c                 dFDR_piston(i)=0d0
                enddo
  
                if(compteur.le.startsavespec) then
                  mCvv_average_piston=0d0
                  Cvf_average_piston=0d0
c                 dFDR_average_piston=0d0
                endif

                Cvf_file='Cvf_piston.out'
                mCvv_file='mCvv_piston.out'

                open(168,file=Cvf_file)
                open(178, file=mCvv_file)

                call dfftw_plan_dft_1d(planv,s,s_in,s_out_v,1,est)
                call dfftw_plan_dft_1d(planr,s,s_in,s_out_r,1,est)
            
                s_in(:)=dcmplx(0d0,0d0)
                do k=1,nseg
                  s_in(k)=dcmplx(vad_piston(k),0d0)
                enddo
                call dfftw_execute(planv,s_in,s_out_v)
                s_in(:)=dcmplx(0d0,0d0)
                do k=1,nseg
                  s_in(k)=dcmplx(fad_piston(k),0d0)
                enddo
                call dfftw_execute(planr,s_in,s_out_r)
                mCvv_piston(:)=mCvv_piston(:)+masspiston*
     &                      abs(s_out_v(1:nad))**2/nseg
                Cvf_piston(:)=Cvf_piston(:)+real(s_out_v(1:nad)
     &                      *conjg(s_out_r(1:nad)/nseg))
              
                mCvv_average_piston=mCvv_average_piston+mCvv_piston
                Cvf_average_piston=Cvf_average_piston+Cvf_piston
c              dFDR_average_piston=dFDR_average_piston+dFDR_piston

                do i=0,nad-1
                  write(168,'('//int_to_str(1+typemax)//'E16.8)') 
     &                    domega*i,Cvf_average_piston(i)
     &                    /max(compteur-startsavespec+1,1)
                  write(178,'('//int_to_str(1+typemax)//'E16.8)')
     &                    domega*i,mCvv_average_piston(i)
     &                  /max(compteur-startsavespec+1,1)
c                 write(188,'('//int_to_str(1+typemax)//'E16.8)')
c     &               domega*l,dFDR_average_piston(:)
c     &              /max(compteur-startsavespec+1,1)
c                 write(170,'('//int_to_str(1+typemax)//'E16.8)') 
c     &              domega*l,gamma_piston(:)
c                 write(98,'('//int_to_str(1+typemax)//'E16.8)') 
c     &              domega*l,gamma_piston(:)
                enddo
                
              close(168)
              close(178)
              call dfftw_destroy_plan(planv)
              call dfftw_destroy_plan(planr)
              endif
            endif
      endif
      if ((mod(istep,nseg).eq.0)) call convolseg
      return
      end
