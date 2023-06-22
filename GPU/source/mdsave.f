c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mdsave  --  save trajectory and restart files  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mdsave" writes molecular dynamics trajectory snapshots and
c     auxiliary files with velocity, force or induced dipole data;
c     also checks for user requested termination of a simulation
c
c
#include "tinker_macro.h"
      subroutine mdsave (istep,dt,epot)
      use atmtyp
      use atmlst
      use atoms   ,only: pbcunwrap,pbcWrap
      use atomsMirror
      use bound
      use boxes
      use domdec
      use dcdmod
      use files
      use group
      use inform
      use iounit
      use mdstuf
      use moldyn
      use mpole
      use output
      use polar
      use potent
      use replicas
      use titles
      use USampling ,only: US_step=>step,prtUS
     &              ,US_enable
      use units
      use mpi
      implicit none
      integer   i,j,k,istep,iloc,nprocloc,count
      integer   ixyz,iind,ivel,ifrc
      integer   iend1,idump,lext
      integer   freeunit,trimtext,moddump
      integer   iproc,ierr,req(nproc*nproc)
      integer   tagmpi,iglob,status(MPI_STATUS_SIZE)
      logical   exist
      real(r_p) dt,pico,wt,epot
      real(r_p),allocatable:: vtemp(:,:),postemp(:,:)
     &         ,atemp(:,:),aalttemp(:,:),uindtemp(:,:)
      integer  ,allocatable:: bufbegsave(:),globsave(:)
     &         ,pbcWrapTemp(:),bufbegpolesave(:),globpolesave(:)
      integer  ,parameter::PERIOD_INBOX_ATOMS=500
      character*7 ext
      character*240 endfile,xyzfile,velfile,frcfile,indfile,exten
      character*3 numberreps
c
      moddump = mod(istep,iwrite)
c
      if (mod(istep,PERIOD_INBOX_ATOMS).eq.0.and.use_bounds)
     &   call bounds
c
      ! Umbrella Sampling step follows
      if (US_enable) US_step = istep+1
      if (moddump.ne.0.and..not.f_mdsave)  return
c
c     wrap coordinates in unit cell
c
      !call molecule(.false.)   ! No need to be called
      if (use_bounds.and..not.f_mdsave) call bounds

!$acc update host(glob,v,a,x,y,z,aalt,epot) async
      if (pbcunwrap) then
!$acc update host(pbcWrap) async
      end if
      if (uindsave.and.use_polar) then
!$acc update host(uind) async
      end if

      if (nproc.eq.1) goto 34
c
c     Send positions,velocities and accelerations to the master
c
c
c      allocate temporary arrays
c
      if (rank.ne.0) then
        allocate (pbcWrapTemp(n))
        allocate (postemp(3,nloc))
        allocate (vtemp(3,nloc))
        allocate (atemp(3,nloc))
        allocate (aalttemp(3,nloc))
        if (uindsave .and. use_polar) then
          allocate (uindtemp(3,nloc))
          uindtemp = 0d0
        end if
        postemp  = 0_re_p
        vtemp    = 0_re_p
        atemp    = 0_re_p
        aalttemp = 0_re_p
        pbcWrapTemp = 0
      else
        allocate (pbcWrapTemp(n))
        allocate (postemp(3,n))
        allocate (vtemp(3,n))
        allocate (atemp(3,n))
        allocate (aalttemp(3,n))
        allocate (bufbegsave(nproc))
        allocate (globsave(n))
        if (uindsave.and.use_polar) then
          allocate (uindtemp(3,n))
          allocate (bufbegpolesave(nproc))
          allocate (globpolesave(n))
          uindtemp       = 0
          bufbegpolesave = 0
          globpolesave   = 0
        end if
        postemp  = 0_re_p
        vtemp    = 0_re_p
        atemp    = 0_re_p
        aalttemp = 0_re_p
        bufbegsave  = 0
        globsave    = 0
        pbcWrapTemp = 0
      end if
c
c     put the arrays to be sent in the right order
c
!$acc wait
      if (rank.ne.0) then
        do i = 1, nloc
          iglob         = glob(i)
          postemp (1,i) = x   (iglob)
          postemp (2,i) = y   (iglob)
          postemp (3,i) = z   (iglob)
          vtemp   (1,i) = v   (1,iglob)
          vtemp   (2,i) = v   (2,iglob)
          vtemp   (3,i) = v   (3,iglob)
          atemp   (1,i) = a   (1,iglob)
          atemp   (2,i) = a   (2,iglob)
          atemp   (3,i) = a   (3,iglob)
          aalttemp(1,i) = aalt(1,iglob)
          aalttemp(2,i) = aalt(2,iglob)
          aalttemp(3,i) = aalt(3,iglob)
          pbcWrapTemp(i)= pbcWrap(iglob)
        end do
        if (uindsave .and. use_polar) then
          do i = 1, npoleloc
            iglob = poleglob(i)
            uindtemp(1,i) = uind(1,iglob)
            uindtemp(2,i) = uind(2,iglob)
            uindtemp(3,i) = uind(3,iglob)
          end do
        end if
      end if
c
c     master get size of all domains
c
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_IRECV(domlen(iproc+1),1,MPI_INT,iproc,tagmpi,
     $     COMM_TINKER,req(tagmpi),ierr)
        end do
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_ISEND(domlen(rank+1),1,MPI_INT,0,tagmpi,
     $   COMM_TINKER,req(tagmpi),ierr)
        call MPI_WAIT(req(tagmpi),status,ierr)
      end if

      if (uindsave .and. use_polar) then
        if (rank.eq.0) then
          do iproc = 1, nproc-1
            tagmpi = iproc + 1
            call MPI_IRECV(domlenpole(iproc+1),1,MPI_INT,iproc,tagmpi,
     $       COMM_TINKER,req(tagmpi),ierr)
          end do
          do iproc = 1, nproc-1
            tagmpi = iproc + 1
            call MPI_WAIT(req(tagmpi),status,ierr)
          end do
        else
          tagmpi = rank + 1
          call MPI_ISEND(domlenpole(rank+1),1,MPI_INT,0,tagmpi,
     $     COMM_TINKER,req(tagmpi),ierr)
          call MPI_WAIT(req(tagmpi),status,ierr)
        end if
      end if
c
c     master get all the indexes
c
      if (rank.eq.0) then
        count = domlen(rank+1)
        do iproc = 1, nproc-1
          if (domlen(iproc+1).ne.0) then
            bufbegsave(iproc+1) = count + 1
            count = count + domlen(iproc+1)
          else
            bufbegsave(iproc+1) = 1
          end if
        end do
c
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_IRECV(globsave(bufbegsave(iproc+1)),domlen(iproc+1),
     $     MPI_INT,iproc,tagmpi,COMM_TINKER,req(tagmpi),ierr)
        end do
        do iproc = 1, nproc-1
          tagmpi = nproc*rank + iproc + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_ISEND(glob,domlen(rank+1),MPI_INT,0,
     $   tagmpi,COMM_TINKER,req(tagmpi),ierr)
        call MPI_WAIT(req(tagmpi),status,ierr)
      end if

      if (uindsave .and. use_polar) then
        if (rank.eq.0) then
          count = domlenpole(rank+1)
          do iproc = 1, nproc-1
            if (domlenpole(iproc+1).ne.0) then
              bufbegpolesave(iproc+1) = count + 1
              count = count + domlenpole(iproc+1)
            else
              bufbegpolesave(iproc+1) = 1
            end if
          end do
c
          do iproc = 1, nproc-1
            tagmpi = iproc + 1
            call MPI_IRECV(globpolesave(bufbegpolesave(iproc+1)),
     $       domlenpole(iproc+1),MPI_INT,iproc,tagmpi,COMM_TINKER,
     $       req(tagmpi),ierr)
          end do
          do iproc = 1, nproc-1
            tagmpi = nproc*rank + iproc + 1
            call MPI_WAIT(req(tagmpi),status,ierr)
          end do
        else
          tagmpi = rank + 1
          call MPI_ISEND(poleglob,domlenpole(rank+1),MPI_INT,0,
     $     tagmpi,COMM_TINKER,req(tagmpi),ierr)
          call MPI_WAIT(req(tagmpi),status,ierr)
        end if
      end if
c
c    send them to the master
c
c
c     positions
c
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_IRECV(postemp(1,bufbegsave(iproc+1)),
     $     3*domlen(iproc+1),
     $     MPI_RPREC,iproc,tagmpi,COMM_TINKER,req(tagmpi),ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_ISEND(postemp,3*nloc,MPI_RPREC,0,tagmpi,
     $  COMM_TINKER,req(tagmpi),ierr)
      end if
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end if
c
c     pbc Wrap State
c
      if (pbcunwrap) then
        if (rank.eq.0) then
          do iproc = 1, nproc-1
            tagmpi = iproc + 1
            call MPI_IRECV(pbcWrapTemp(bufbegsave(iproc+1))
     $       ,domlen(iproc+1),
     $       MPI_INT,iproc,tagmpi,COMM_TINKER,req(tagmpi),ierr)
          end do
        else
          tagmpi = rank + 1
          call MPI_ISEND(pbcWrapTemp,nloc,MPI_INT,0,tagmpi
     $    ,COMM_TINKER,req(tagmpi),ierr)
        end if
        if (rank.eq.0) then
          do iproc = 1, nproc-1
            tagmpi = iproc + 1
            call MPI_WAIT(req(tagmpi),status,ierr)
          end do
        else
          tagmpi = rank + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end if
      end if
c
c     velocities
c
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_IRECV(vtemp(1,bufbegsave(iproc+1)),3*domlen(iproc+1),
     $     MPI_RPREC,iproc,tagmpi,COMM_TINKER,req(tagmpi),ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_ISEND(vtemp,3*nloc,MPI_RPREC,0,tagmpi,COMM_TINKER,
     $  req(tagmpi),ierr)
      end if
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end if
c
c     accelerations
c
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_IRECV(atemp(1,bufbegsave(iproc+1)),3*domlen(iproc+1),
     $     MPI_RPREC,iproc,tagmpi,COMM_TINKER,req(tagmpi),ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_ISEND(atemp,3*nloc,MPI_RPREC,0,tagmpi,COMM_TINKER,
     $  req(tagmpi),ierr)
      end if
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end if
c
c     alternate accelerations
c
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_IRECV(aalttemp(1,bufbegsave(iproc+1)),
     $     3*domlen(iproc+1)
     $     ,MPI_RPREC,iproc,tagmpi,COMM_TINKER,req(tagmpi),ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_ISEND(aalttemp,3*nloc,MPI_RPREC,0,tagmpi,
     $   COMM_TINKER,req(tagmpi),ierr)
      end if
      if (rank.eq.0) then
        do iproc = 1, nproc-1
          tagmpi = iproc + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end do
      else
        tagmpi = rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end if
c
c
c     induced dipoles
c 
      if (uindsave .and. use_polar) then
        if (rank.eq.0) then
          do iproc = 1, nproc-1
            tagmpi = iproc + 1
            call MPI_IRECV(uindtemp(1,bufbegpolesave(iproc+1)),
     $       3*domlenpole(iproc+1),
     $       MPI_REAL8,iproc,tagmpi,COMM_TINKER,req(tagmpi),ierr)
          end do
        else
          tagmpi = rank + 1
          call MPI_ISEND(uindtemp,3*npoleloc,MPI_REAL8,0,tagmpi,
     $    COMM_TINKER,req(tagmpi),ierr)
        end if
        if (rank.eq.0) then
          do iproc = 1, nproc-1
            tagmpi = iproc + 1
            call MPI_WAIT(req(tagmpi),status,ierr)
          end do
        else
          tagmpi = rank + 1
          call MPI_WAIT(req(tagmpi),status,ierr)
        end if
      end if
c
c     put the array in global order to be written
c
      if (rank.eq.0) then
        if (use_pmecore) then
          nprocloc = ndir
        else
          nprocloc = nproc
        end if
        do iproc = 1, nprocloc-1
          do i = 1, domlen(iproc+1)
            if (domlen(iproc+1).ne.0) then
              iloc  = bufbegsave(iproc+1)+i-1
              iglob = globsave(iloc)
              x(iglob) = postemp(1,iloc)
              y(iglob) = postemp(2,iloc)
              z(iglob) = postemp(3,iloc)
              pbcWrap(iglob) = pbcWrapTemp(iloc)
              v(1,iglob) = vtemp(1,iloc)
              v(2,iglob) = vtemp(2,iloc)
              v(3,iglob) = vtemp(3,iloc)
              a(1,iglob) = atemp(1,iloc)
              a(2,iglob) = atemp(2,iloc)
              a(3,iglob) = atemp(3,iloc)
              aalt(1,iglob) = aalttemp(1,iloc)
              aalt(2,iglob) = aalttemp(2,iloc)
              aalt(3,iglob) = aalttemp(3,iloc)
            end if
          end do
          if (uindsave .and. use_polar) then
            do i = 1, domlenpole(iproc+1)
              if (domlenpole(iproc+1).ne.0) then
                iloc = bufbegpolesave(iproc+1)+i-1
                iglob = globpolesave(iloc)
                uind(1,iglob) = uindtemp(1,iloc)
                uind(2,iglob) = uindtemp(2,iloc)
                uind(3,iglob) = uindtemp(3,iloc)
              end if
            end do
          end if
c         if (frcsave) then
c           do i = 1, domlen(iproc+1)
c             if (domlen(iproc+1).ne.0) then
c               iloc = bufbegsave(iproc+1)+i-1
c               iglob = globsave(iloc)
c               derivstemp2(1,iglob) = derivstemp(1,iloc)
c               derivstemp2(2,iglob) = derivstemp(2,iloc)
c               derivstemp2(3,iglob) = derivstemp(3,iloc)
c             end if
c           end do
c         end if
        end do
      end if
      deallocate (postemp)
c
      deallocate (vtemp)
      deallocate (atemp)
      deallocate (aalttemp)
      if (uindsave .and. use_polar) deallocate (uindtemp)
      if (rank.ne.0) return
 34   continue
!$acc wait
c
c     get the sequence number of the current trajectory frame
c
      if (f_mdsave) then
         idump = merge(0,n_fwriten,n_fwriten.eq.-1)
         lext  = merge(1,        6,n_fwriten.eq.-1)
         call numeral (idump,ext,lext)
      else
         idump = nprior + istep/iwrite
         lext  = 3
         call numeral (idump,ext,lext)
      end if
c
c     print header for the instantaneous values at current step
c
      pico = real(istep,t_p) * dt
      write (iout,10)  istep
   10 format (/,' Instantaneous Values for Frame saved at',
     &           i10,' Dynamics Steps')
      write (iout,20)  pico
   20 format (/,' Current Time',8x,f15.4,' Picosecond')
      write (iout,30)  epot
   30 format (' Current Potential',3x,f15.4,' Kcal/mole')
      if (use_bounds) then
         write (iout,40)  xbox,ybox,zbox
   40    format (' Lattice Lengths',6x,3f14.6)
         write (iout,50)  alpha,beta,gamma
   50    format (' Lattice Angles',7x,3f14.6)
      end if
      write (iout,60)  idump
   60 format (' Frame Number',13x,i10)
c
c     update the information needed to restart the trajectory
c
      if (mod(istep,idumpdyn).eq.0) call prtdyn
c
c     save coordinates to an archive or numbered structure file, or dcd file
c
      if (dcdio) then
        call dcdio_write(istep,dt)
      else
        ixyz = freeunit ()
        if (archive) then
           xyzfile = filename(1:leng)
c
c          if multiple replicas, then number the traj outputs
c
           if (use_reps) then
             write(numberreps, '(i3.3)') rank_reploc
             exten='_reps'//numberreps
             xyzfile = filename(1:leng)//trim(exten)
           end if
           call suffix (xyzfile,'arc','old')

           inquire (file=xyzfile,exist=exist)
           if (exist) then
              call openend (ixyz,xyzfile)
           else
              open (unit=ixyz,file=xyzfile,status='new')
           end if
        else
           if (use_reps) then
             write(numberreps, '(i3.3)') rank_reploc
             exten='_reps'//numberreps
             xyzfile = filename(1:leng)//trim(exten)
     &                            //'.'//ext(1:lext)
           else
             xyzfile = filename(1:leng)//'.'//ext(1:lext)
           end if
           if (f_mdsave) then
              if (n_fwriten.eq.-1) then
              xyzfile = filename(1:leng)//'_err'
              else
              xyzfile = filename(1:leng)//'_'//ext(1:lext)
              end if
           else
             xyzfile = filename(1:leng)//'.'//ext(1:lext)
           end if
           call version (xyzfile,'new')
           open (unit=ixyz,file=xyzfile,status='new')
        end if
        call prtxyz (ixyz)
        close (unit=ixyz)
        write (iout,70)  xyzfile(1:trimtext(xyzfile))
   70   format (' Coordinate File',12x,a)
      end if
c
c     Save data for fred
c
      if (US_enable) call prtUS

c
c     save the velocity vector components at the current step
c
      if (velsave) then
         ivel = freeunit ()
         if (archive) then
            velfile = filename(1:leng)
            call suffix (velfile,'vel','old')
            inquire (file=velfile,exist=exist)
            if (exist) then
               call openend (ivel,velfile)
            else
               open (unit=ivel,file=velfile,status='new')
            end if
         else
            velfile = filename(1:leng)//'.'//ext(1:lext)//'v'
            call version (velfile,'new')
            open (unit=ivel,file=velfile,status='new')
         end if
         write (ivel,110)  n,title(1:ltitle)
  110    format (i6,2x,a)
         do i = 1, n
            write (ivel,120)  i,name(i),(v(j,i),j=1,3)
  120       format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
         end do
         close (unit=ivel)
         write (iout,130)  velfile(1:trimtext(velfile))
  130    format (' Velocity File',15x,a)
      end if
c
c     save the force vector components for the current step
c
      if (frcsave) then
         ifrc = freeunit ()
         if (archive) then
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc','old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               call openend (ifrc,frcfile)
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
         end if
         write (ifrc,140)  n,title(1:ltitle)
  140    format (i6,2x,a)
         do i = 1, n
            !TODO 1.2 Ask Louis for this part
            wt = mass(i) / convert
            write (ifrc,150)  i,name(i),(wt*a(j,i),j=1,3)
  150       format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
         end do
         close (unit=ifrc)
         write (iout,160)  frcfile(1:trimtext(frcfile))
  160    format (' Force Vector File',11x,a)
      end if
c
c     save the current induced dipole moment at each site
c
      if (uindsave .and. use_polar) then
         iind = freeunit ()
         if (archive) then
            indfile = filename(1:leng)
            call suffix (indfile,'uind','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
         else
            indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
         end if
         write (iind,170)  n,title(1:ltitle)
  170    format (i6,2x,a)
         do i = 1, npole
            if (polarity(i) .ne. 0.0_ti_p) then
               k = ipole(i)
               write (iind,180)  k,name(k),(debye*uind(j,i),j=1,3)
  180          format (i6,2x,a3,3f12.6)
            end if
         end do
         close (unit=iind)
         write (iout,190)  indfile(1:trimtext(indfile))
  190    format (' Induced Dipole File',10x,a)
      end if
c
c     test for requested termination of the dynamics calculation
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend1 = freeunit ()
            open (unit=iend1,file=endfile,status='old')
            close (unit=iend1,status='delete')
         end if
      end if
      if (exist) then
         write (iout,200)
  200    format (/,' MDSAVE  --  Dynamics Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
c
c     skip an extra line to keep the output formating neat
c
      moddump = mod(istep,iprint)
      if (verbose .and. moddump.ne.0) then
         write (iout,210)
  210    format ()
      end if
      end
