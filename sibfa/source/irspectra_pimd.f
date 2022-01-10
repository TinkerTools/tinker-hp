
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Tex_sizeas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine irspectra_pimd--initialize a dynamics trajectory   ##
!     ##                                                           ##
!     ###############################################################
!

c     subroutine irspectra_pimd: get the total dipole moment of the system (permanent+induced)

      subroutine irspectra_pimd(dt)
      use atmlst
      use atmtyp
      use atoms
      use bath
      use boxes, only: volbox
      use charge
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use math
      use mdstuf
      use moldyn
      use mpole
      use potent
      use qtb
      use adqtb
      use timestat
      use usage
      use units
      use mpi
      implicit none
      integer i,j,l,istep,ierr
      integer iglob,iipole,iichg
      real*8 oterm,hterm,q,beta,k_b,dt
      real*8, allocatable ::  vtot(:,:,:)
      real*8, allocatable :: Cmumu(:,:)
      real*8, allocatable ::  mudot(:,:)
      double complex, allocatable :: s_in(:), s_out(:)
      integer*8 plan, est,s

      oterm = 0.73612d0
      hterm = 0.13194d0

      s=3*nseg
      est=1
      k_b=boltzmann
      beta=1./(k_b*kelvin)
      domega = (2.*pi)/(3.*nseg*dt)
      nad=int(omegacut/domega)

      allocate(vtot(3,n,nseg))
      allocate(mudot(3,nseg))
     
      allocate(Cmumu(3,0:nad-1))

      Cmumu(:,:)=0d0
      mudot(:,:)=0d0

      if(compteur.le.startsavespec) then
        if (allocated(Cmumu_average)) deallocate (Cmumu_average)
         allocate(Cmumu_average(3,0:nad-1)) 
        Cmumu_average(:,:)=0d0
      endif

      open(12,file='IR_spectra_decomposed.out')
      open(13,file='IR_spectra.out')
      
      
      do i=1,n
          if (use(i) .and. (atomic(i).eq.0)) then
            do j=1,3
              vtot(1,i,:)=oterm*vad(1,i-3,:)+hterm*vad(1,i-2
     &                        ,:)+hterm*vad(1,i-1,:) 
              vtot(2,i,:)=oterm*vad(2,i-3,:)+hterm*vad(2,i-2
     &                        ,:)+hterm*vad(2,i-1,:) 
              vtot(3,i,:)=oterm*vad(3,i-3,:)+hterm*vad(3,i-2
     &                        ,:)+hterm*vad(3,i-1,:) 
           enddo
          else if (use(i)) then 
            vtot(1,i,:)=vad(1,i,:)
            vtot(2,i,:)=vad(2,i,:)
            vtot(3,i,:)=vad(3,i,:)
          endif
      enddo


        if (use_mpole) then
          do i= 1, npole
            iglob = ipole(i)
            q=rpole(1,i)
              do j=1,3
                mudot(j,:)=mudot(j,:)+q*(vtot(j,iglob,:))
              enddo
            enddo
        else if (use_charge) then
          do i = 1, nion
            iglob = iion(i)
            q = pchg(i)
              do j=1,3
                mudot(j,:)=mudot(j,:)+q*(vtot(j,iglob,:))
              enddo
          enddo
        endif
      

        
      allocate(s_in(3*nseg))   
      allocate(s_out(3*nseg))


      call dfftw_plan_dft_1d(plan,s,s_in,s_out,1,est)
        do j=1,3
          s_in(:)=dcmplx(0d0,0d0)
          s_in(1:nseg)=dcmplx(mudot(j,:),0d0)
          call dfftw_execute(plan,s_in,s_out)
          Cmumu(j,:)=Cmumu(j,:)+abs(s_out(1:nad))**2/nseg
        enddo

         Cmumu_average=Cmumu_average+Cmumu
        
        do i=0,nad-1
         write(12,'(4e20.5)')domega*i,Cmumu_average(:,i)
     &                /max(compteur-startsavespec+1,1)
         write(13,*) domega*i, ((pi*beta)/(3*lightspd*volbox))
     &                *(Cmumu_average(1,i)
     &                +Cmumu_average(2,i)
     &                +Cmumu_average(3,i))
     &                /max(compteur-startsavespec+1,1)
        enddo
      
      close(12)
      close(13)
      
      call dfftw_destroy_plan(plan) 
    


      return 
      end
