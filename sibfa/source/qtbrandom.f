
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Tex_sizeas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine qtbrandom--initialize a dynamics trajectory   ##
!     ##                                                           ##
!     ###############################################################
!

      subroutine qtbrandom(dt,istep)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use math
      use mdstuf
      use moldyn
      use qtb
      use adqtb
      use timestat
      use usage
      use units
      use mpi
      implicit none
      integer i,j,k,l, istep
      integer iglob
      real*8 a1,a2,dt
      real*8 normal
      a1 = exp(-gamma*dt)
      a2 = sqrt(2*gamma*dt)

!     create our random noise
!

      k = mod(istep-1,nseg)+1
      if (adaptive) then
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob) .AND. (atomic(iglob) /= 0)) then
            do j = 1, 3
                  vad(j,i,k) =sqrt(a1)*v(j,iglob)+
     $            (a2*rt(j,i,k)/sqrt(mass(iglob)))/2 
                 fad(j,i,k) = a2*rt(j,i,k)*sqrt(mass(iglob))/dt
            end do
          endif
        end do
      end if


      do i=1,nloc
            iglob=glob(i)
            if (use(iglob) .AND. (atomic(iglob) /= 0)) then
                  do j=1,3
                        v(j,iglob)=(a1*v(j,iglob)+
     $                  a2*rt(j,i,k)/sqrt(mass(iglob)))
                  enddo
            endif
      enddo

      end      
