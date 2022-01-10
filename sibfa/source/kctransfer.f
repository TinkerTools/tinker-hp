c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine kctransfer  --  charge transfer parameter assignment  ##
c     ##                                                                   ##
c     #######################################################################
c
c     kctransfer reads a lone pair coordinates file, translates these coordinates in
c     the global frame and in cartesian coordinates, builds an array of electron acceptor
c     bonds and get charge transfer energy parameters
c
      subroutine kctransfer(init,istep)
      use atoms
      use atmtyp
      use atmlst
      use bond
      use domdec
      use keys
      use potent
      use chargetransfer
      use kct
      use couple
      use cutoff
      use neigh
      use repulsion
      use mpi
      implicit none
      integer i,j,next,k,iglob,iproc,istep,modnl
      integer l,m
      integer acceptcount,nacceptloc1
      integer iaccept,tagmpi
      integer status(MPI_STATUS_SIZE),ierr,req(nproc*nproc),count
      real*8 d
      logical header,init,docompute
      character*20 keyword
      character*120 record
      character*20 string
c
      if (init) then
c
c     defaults for charge transfer term (and parameters)
c
        use_ctransfer = .false.
        use_ctransferlist = .false.
        use_ctpot = .false.
c
c       search keywords for charge transfer commands
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           string = record(next:120)
           if (keyword(1:13) .eq. 'CTRANSFERTERM ') then
              use_ctransfer = .true.
              use_ctransferlist = .true.
           end if
           if (keyword(1:8) .eq. 'CTPOT ') then
              use_ctpot = .true.
           end if
        end do
c
c       hardcode some parameters
c
        call sibfaparam
c
        call readlp
c
        if (allocated(nbaccept)) deallocate (nbaccept)
        allocate (nbaccept(n))
        call acceptorlp

        if (allocated(lplen)) deallocate (lplen)
        allocate (lplen(nproc))
        if (allocated(acceptlen)) deallocate (acceptlen)
        allocate (acceptlen(nproc))
        if (allocated(bufbegacc)) deallocate (bufbegacc)
        allocate (bufbegacc(nproc))
        if (allocated(bufbeglp)) deallocate (bufbeglp)
        allocate (bufbeglp(nproc))
c
c       use atom class or type as index into ctransfer parameters
c
c         allocate global arrays
c
        if (allocated(chyb)) deallocate (chyb)
        allocate (chyb(2,n))
        if (allocated(ahct)) deallocate (ahct)
        allocate (ahct(n))
        if (allocated(aelct)) deallocate (aelct)
        allocate (aelct(n))
        k = 0
        do i = 1, n
           rvdwct1(i) = sibfact1(type(i))
           rvdwct2(i) = sibfact2(type(i))
           chyb(1,i)    = hybrid(1,type(i))
           chyb(2,i)    = hybrid(2,type(i))
           ahct(i) = ah(type(i)) 
           aelct(i) = ae(type(i)) 
        end do
c
c      get the correspondance between atoms and acceptors (the H atoms
c       because the heavy ones can belong to several acceptors)
c
        if (allocated(acceptlist)) deallocate (acceptlist)
        allocate (acceptlist(n))
        call iclear(n,acceptlist)
        do i = 1, n
          iaccept = nbaccept(i)
          if (atomic(i).ne.1) then
            do j = 1, n12(i)
              k = i12(j,i)
              if (atomic(k).eq.1) then
                iaccept = iaccept + 1
                acceptlist(k) = iaccept
              end if
            end do
          end if
          if (n12(i).eq.0) then
            iaccept = iaccept + 1
            acceptlist(i) = iaccept
          end if
        end do
      end if
c
c     get local lone pairs
c
      if (allocated(lpglob)) deallocate(lpglob)
      allocate (lpglob(8*nbloc))
      lpglob = 0
c
      nlploc = 0
      do i = 1, nloc
        iglob = glob(i)
        do j = 1, nilplst(iglob)
          nlploc = nlploc + 1
          lpglob(nlploc) = ilplst(j,iglob)
        end do
      end do
c      nlpbloc = nlploc
      lplen(rank+1) = nlploc
c      do iproc = 1, nbig_recep
c        do i = 1, domlen(pbig_recep(iproc)+1)
c          iglob = glob(bufbeg(pbig_recep(iproc)+1)+i-1)
c          do j = 1, nilplst(iglob)
c            nlpbloc = nlpbloc + 1
c            lpglob(nlpbloc) = ilplst(j,iglob)
c          end do
c        end do
c      end do
c
c     get local acceptors
c
      if (allocated(acceptglob)) deallocate(acceptglob)
      allocate (acceptglob(4*nbloc))
      if (allocated(acceptloc)) deallocate(acceptloc)
      allocate (acceptloc(n))
      call iclear(4*nbloc,acceptglob)
      call iclear(n,acceptloc)
      nacceptloc = 0
      do i = 1, nloc
        iglob = glob(i)
        acceptcount = nbaccept(iglob)
        nacceptloc1 = 0
        if (atomic(iglob).ne.1) then
          do j = 1, n12(iglob)
            k = i12(j,iglob)
c
c     Register the X-H atom pairs            
c
            if (atomic(k).eq.1) then
              nacceptloc = nacceptloc + 1
              nacceptloc1 = nacceptloc1 + 1
              acceptglob(nacceptloc) = acceptcount + nacceptloc1
              acceptloc(acceptcount+nacceptloc1) = nacceptloc
            end if
          end do
        end if
        if (n12(iglob).eq.0) then
          nacceptloc = nacceptloc + 1
          nacceptloc1 = nacceptloc1 + 1
          acceptglob(nacceptloc) = acceptcount + nacceptloc1
          acceptloc(acceptcount+nacceptloc1) = nacceptloc
        end if
      end do  
c
c     get the acceptlen values of the coupled real space processes
c

      acceptlen(rank+1) = nacceptloc
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_IRECV(acceptlen(pbig_recep(iproc)+1),1,MPI_INT,
     $   pbig_recep(iproc),tagmpi,MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(acceptlen(rank+1),1,MPI_INT,pbig_send(iproc),
     $   tagmpi,MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do
c
      bufbegacc(rank+1) = 1
      count = acceptlen(rank+1)
      do iproc = 1, nbig_recep
        if (acceptlen(pbig_recep(iproc)+1).ne.0) then
          bufbegacc(pbig_recep(iproc)+1) = count + 1
        else
          bufbegacc(pbig_recep(iproc)+1) = 1
        end if
        count = count + acceptlen(pbig_recep(iproc)+1)
      end do
      nacceptbloc = count
c
c     send and receive the accept indexes of the neighboring processes
c
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_IRECV(acceptglob(bufbegacc(pbig_recep(iproc)+1)),
     $   acceptlen(pbig_recep(iproc)+1),MPI_INT,pbig_recep(iproc),
     $   tagmpi,MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(acceptglob,acceptlen(rank+1),MPI_INT,
     $   pbig_send(iproc),tagmpi,MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do
c
      do iproc = 1, nbig_recep
        do i = 1, acceptlen(pbig_recep(iproc)+1)
          acceptloc(acceptglob(bufbegacc(pbig_recep(iproc)+1)+i-1)) =
     $      bufbegacc(pbig_recep(iproc)+1)+i-1
        end do
      end do
c
c     get the lplen values of the coupled real space processes
c
      lplen(rank+1) = nlploc
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_IRECV(lplen(pbig_recep(iproc)+1),1,MPI_INT,
     $   pbig_recep(iproc),tagmpi,MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(lplen(rank+1),1,MPI_INT,pbig_send(iproc),tagmpi,
     $   MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do
c
      bufbeglp(rank+1) = 1
      count = lplen(rank+1)
      do iproc = 1, nbig_recep
        if (lplen(pbig_recep(iproc)+1).ne.0) then
          bufbeglp(pbig_recep(iproc)+1) = count + 1
        else
          bufbeglp(pbig_recep(iproc)+1) = 1
        end if
        count = count + lplen(pbig_recep(iproc)+1)
      end do
      nlpbloc = count
c
c     send and receive the lp indexes of the neighboring processes
c
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_IRECV(lpglob(bufbeglp(pbig_recep(iproc)+1)),
     $   lplen(pbig_recep(iproc)+1),MPI_INT,pbig_recep(iproc),
     $   tagmpi,MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_ISEND(lpglob,lplen(rank+1),MPI_INT,
     $   pbig_send(iproc),tagmpi,MPI_COMM_WORLD,req(tagmpi),ierr)
      end do
      do iproc = 1, nbig_recep
        tagmpi = nproc*rank + pbig_recep(iproc) + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do
      do iproc = 1, nbig_send
        tagmpi = nproc*pbig_send(iproc) + rank + 1
        call MPI_WAIT(req(tagmpi),status,ierr)
      end do

c      do iproc = 1, nbig_recep
c        do i = 1, lplen(pbig_recep(iproc)+1)
c          lploc(lpglob(bufbeglp(pbig_recep(iproc)+1)+i-1)) =
c     $      bufbeglp(pbig_recep(iproc)+1)+i-1
c        end do
c      end do
      call rotlp
c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return
c
      if (allocated(lpglobnl)) deallocate(lpglobnl)
      allocate (lpglobnl(nlp))
      if (allocated(lplocnl)) deallocate(lplocnl)
      allocate (lplocnl(nlp))
c
      nlplocnl = 0
      do i = 1, nlocnl
        iglob = ineignl(i)
        call distprocpart(iglob,rank,d,.true.)
        if (repart(iglob).eq.rank) d = 0.0d0
        if (d*d.gt.(ctransferbuf2/4)) cycle
        do j = 1, nilplst(iglob)
          nlplocnl = nlplocnl + 1
          lpglobnl(nlplocnl) = ilplst(j,iglob)
          lplocnl(ilplst(j,iglob)) = nlplocnl
        end do
      end do
c
      if (allocated(acceptglobnl)) deallocate(acceptglobnl)
      allocate (acceptglobnl(naccept))
      if (allocated(acceptlocnl)) deallocate(acceptlocnl)
      allocate (acceptlocnl(naccept))
c
c     faire une neighbor list dediee pour interaction elec sur accepteurs et lp
c
      nacceptlocnl = 0
      do i = 1, nlocnl
        iglob = ineignl(i)
        call distprocpart(iglob,rank,d,.true.)
        if (repart(iglob).eq.rank) d = 0.0d0
        if (d*d.gt.(ctransferbuf2/4)) cycle
        if (atomic(iglob).ne.1) then
          do j = 1, n12(iglob)
           k = i12(j,iglob)
           call distprocpart(k,rank,d,.true.)
           if (repart(k).eq.rank) d = 0.0d0
           if (d*d.gt.(ctransferbuf2/4)) cycle
           if (atomic(k).eq.1) then
              iaccept = acceptlist(k)
              nacceptlocnl = nacceptlocnl + 1
              acceptlocnl(iaccept) = nacceptlocnl
              acceptglobnl(nacceptlocnl) = iaccept
            end if
          end do
        end if
        if (n12(iglob).eq.0) then
          nacceptlocnl = nacceptlocnl + 1
          acceptlocnl(iaccept) = nacceptlocnl
          acceptglobnl(nacceptlocnl) = iaccept
        end if
      end do
      return
      end
c
c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine readlp  -- get Lone Pair local    coordinates         ##
c     ##                                                                   ##
c     #######################################################################
c
c    Read a lone pair input file and stores the data in 3 arrays :
c     - lonepair(3,maxlp) : internal geometrical data
c     - lpcharge(maxlp) : partial charge of eache lone pair
c     - lpatom(maxlp) : index of the corresponding atom for each lone pair
c
c
      subroutine readlp
      use atoms
      use atmtyp
      use chargetransfer
      use files
      use iounit
      use kct
      use repulsion
      implicit none
      integer i,ilp,k
      integer freeunit
      integer next,size
      integer first,last
      integer nexttext,trimtext
      integer ltitle
      logical quit,abort,opened,exist
      character*120 lpfile
      character*120 record
      character*120 string
      character*120 title
c
      ilp = freeunit ()
c
c     initialize total number of lone pair in the system
      nlp = 0
c
c     open the input file if it has not already been done
      inquire (unit=ilp,opened=opened)
      if (.not. opened) then
         lpfile = filename(1:leng)//'.lp'
         inquire (file=lpfile,exist=exist)
         if (exist) then
            open (unit=ilp,file=lpfile,status='old')
            rewind (unit=ilp)
         else
            write (iout,10)
   10       format (/,' READLP  --  Unable to Find the LP',
     &                 ' Coordinates File')
            call fatal
         end if
      end if
c
c     read first line and return if already at end of file
c
      quit = .false.
      size = 0
      do while (size .eq. 0)
         read (ilp,20,err=70,end=70)  record
   20    format (a120)
         size = trimtext (record)
      end do
      quit = .true.
c
c     parse the title line to get the number of atoms
c
      i = 0
      next = 1
      call gettext (record,string,next)
      read (string,*,err=70,end=70)  nlp
c
c     allocate global arrays
c
      if (allocated(lpatom)) deallocate (lpatom)
      allocate (lpatom(nlp))
      if (allocated(acceptor)) deallocate (acceptor)
      allocate (acceptor(2,n))
      if (allocated(izlp)) deallocate (izlp)
      allocate (izlp(nlp))
      if (allocated(ixlp)) deallocate (ixlp)
      allocate (ixlp(nlp))
      if (allocated(lonepair)) deallocate (lonepair)
      allocate (lonepair(3,nlp))
      if (allocated(rlonepair)) deallocate (rlonepair)
      allocate (rlonepair(3,nlp))
      if (allocated(lpcharge)) deallocate (lpcharge)
      allocate (lpcharge(nlp))
      if (allocated(rvdwct1)) deallocate (rvdwct1)
      allocate (rvdwct1(n))
      if (allocated(rvdwct2)) deallocate (rvdwct2)
      allocate (rvdwct2(n))
      if (allocated(namelp)) deallocate (namelp)
      allocate (namelp(nlp))
      if (allocated(hybrid_lp)) deallocate (hybrid_lp)
      allocate (hybrid_lp(2,nlp))
      if (allocated(dincr_lpect)) deallocate (dincr_lpect)
      allocate (dincr_lpect(nlp))
      if (allocated(dincr_lprep)) deallocate (dincr_lprep)
      allocate (dincr_lprep(nlp))
      if (allocated(ilplst)) deallocate (ilplst)
      allocate (ilplst(8,n))
      if (allocated(nilplst)) deallocate (nilplst)
      allocate (nilplst(n))
c
c     extract the title and determine its length
c
      string = record(next:120)
      first = nexttext (string)
      last = trimtext (string)
      if (last .eq. 0) then
         title = ' '
         ltitle = 0
      else
         title = string(first:last)
         ltitle = trimtext (title)
      end if
c
c     check for too few or too many total atoms in the file
c
      if (nlp .le. 0) then
         write (iout,30)
   30    format (/,' READLP  --  The Coordinate File Does Not',
     &              ' Contain Any Lone Pair')
c         call fatal
      else if (nlp .gt. 8*n) then
         write (iout,40)  8*n
   40    format (/,' READLP  --  The Maximum of',i8,' Lone Pair',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     initialize coordinates, charges and atoms for each lp
c
      do i = 1, nlp
         lonepair(1,i) = 0.0d0
         lonepair(2,i) = 0.0d0
         lonepair(3,i) = 0.0d0
         ixlp(i)  = 0
         izlp(i)  = 0
         lpcharge(i) = 0.0d0
         lpatom(i) = 0.0d0
         namelp(i) = '    '
      end do
c
c     read the coordinates, charges and atoms for each lp
c
      do i = 1, nlp
         next = 1
         size = 0
         do while (size .eq. 0)
            read (ilp,50,err=70,end=70)  record
   50       format (a120)
            size = trimtext (record)
         end do
         read (record,*,err=70,end=70)  lpatom(i)
         call getword (record,namelp(i),next)
         string = record(next:120)
         read (string,*,err=60,end=60)  lonepair(1,i),lonepair(2,i),
     $     lonepair(3,i),ixlp(i),izlp(i),lpcharge(i),hybrid_lp(1,i),
     $     hybrid_lp(2,i),dincr_lprep(i),dincr_lpect(i)
   60    continue
      end do
      quit = .false.
   70 continue
      if (.not. opened)  close (unit=ilp)
c
c     an error occurred in reading the coordinate file
c
      if (quit) then
         write (iout,80)  i
   80    format (/,' READLP  --  Error in Coordinate File at LP',i6)
         call fatal
      end if
c
      call iclear(8*n,ilplst)
      call iclear(n,nilplst)
c
      do i = 1, nlp
        k = lpatom(i)
        nilplst(k) = nilplst(k) + 1
        ilplst(nilplst(k),k) = i
      end do
c
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotlp  --  rotate lone pairs to global frame    ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotlp" constructs the set of lone pairs in the global
c     frame in cartesian coordinates
c
c
      subroutine rotlp
      use atoms
      use atmlst
      use chargetransfer
      use domdec
      implicit none
      integer i,ii,k,l,ilp
      real*8 ri(3),riz(3),rix(3),a(3,3)
c
c
c     get the global cartesian coordinates for the lp
c
      do ilp = 1, nlpbloc
        i = lpglob(ilp)
        ii = lpatom(i)
        ri(1) = x(ii)
        ri(2) = y(ii)
        ri(3) = z(ii)
        rix(1) = x(ixlp(i))
        rix(2) = y(ixlp(i))
        rix(3) = z(ixlp(i))
        riz(1) = x(izlp(i))
        riz(2) = y(izlp(i))
        riz(3) = z(izlp(i))
        call rotmatlp(a,ri,rix,riz,1)
        do k = 1, 3
          rlonepair(k,i) = 0.0d0
          do l = 1, 3
            rlonepair(k,i) = rlonepair(k,i)+a(k,l)*lonepair(l,i)
          end do
          rlonepair(k,i) = ri(k) + rlonepair(k,i)
        end do
      end do
      return
      end
c
c     ###########################################################################
c     ##                                                                       ##
c     ##  subroutine acceptorlp  --  acceptor list for charge transfer energy  ##
c     ##                                                                       ##
c     ###########################################################################
c
c
c     "acceptorlp" sets up the array containing the indexes corresponding to the couples
c     X-H forming an electron accepting bond
c
c
      subroutine acceptorlp
      use atmtyp
      use atoms
      use chargetransfer
      use couple
      use kct
      implicit none
      integer i,j,k
      integer iaccept
c
c     Loop over the non hydrogen atoms and check if they are linked to an hydrogen one
c
      naccept = 0
      iaccept = 0
      do i = 1, n
        nbaccept(i) = iaccept
        if (atomic(i).ne.1) then
          do j = 1, n12(i)
            k = i12(j,i)
c
c     Register the X-H atom pairs            
c
            if (atomic(k).eq.1) then
              iaccept = iaccept + 1
              acceptor(1,iaccept) = i
              acceptor(2,iaccept) = k
            end if
          end do
        end if
        if (n12(i).eq.0) then
          iaccept = iaccept + 1
          acceptor(1,iaccept) = i
          acceptor(2,iaccept) = 0
        end if
      end do  
      naccept = iaccept
      return
      end
c
c     rotmatlp gives the rotation matrix to switch from local frame to global frame 
c     for lp coordinates stored with the  Z then X scheme
c
      subroutine rotmatlp(a,ri,rix,riz,i)
      implicit none
      real*8 random
      real*8 r,dot
      real*8 xi,yi,zi
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 dx3,dy3,dz3
      real*8 a(3,3),ri(3),rix(3),riz(3)
      real*8 xix,yix,zix,xiz,yiz,ziz,xiy,yiy,ziy
      integer i
c
c
c     get coordinates and frame definition for the multipole site
c
      xi = ri(1)
      yi = ri(2)
      zi = ri(3)
      xiz = riz(1)
      yiz = riz(2)
      ziz = riz(3)
      xix = rix(1)
      yix = rix(2)
      zix = rix(3)
c
c     use the identity matrix as the default rotation matrix
c
      a(1,1) = 1.0d0
      a(2,1) = 0.0d0
      a(3,1) = 0.0d0
      a(1,3) = 0.0d0
      a(2,3) = 0.0d0
      a(3,3) = 1.0d0
c
c  Z-Only scheme for local frame of the lone pair
c
      if (i.eq.0) then
        dx = xiz - xi
        dy = yiz - yi
        dz = ziz - zi
        r = sqrt(dx*dx + dy*dy + dz*dz)
        a(1,3) = dx / r
        a(2,3) = dy / r
        a(3,3) = dz / r
        dx = random ()
        dy = random ()
        dz = random ()
        dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
        dx = dx - dot*a(1,3)
        dy = dy - dot*a(2,3)
        dz = dz - dot*a(3,3)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        a(1,1) = dx / r
        a(2,1) = dy / r
        a(3,1) = dz / r
c
c  Z-then-X scheme for local frame of the lone pair
c
      else if (i.eq.1) then
        dx = xiz - xi
        dy = yiz - yi
        dz = ziz - zi
        r = sqrt(dx*dx + dy*dy + dz*dz)
        a(1,3) = dx / r
        a(2,3) = dy / r
        a(3,3) = dz / r
        dx = xix - xi
        dy = yix - yi
        dz = zix - zi
        dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
        dx = dx - dot*a(1,3)
        dy = dy - dot*a(2,3)
        dz = dz - dot*a(3,3)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        a(1,1) = dx / r
        a(2,1) = dy / r
        a(3,1) = dz / r
c
c  Bisector scheme for local frame of the lone pair
c
      else if (i.eq.2) then
         dx = xiz - xi
         dy = yiz - yi
         dz = ziz - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = xix - xi
         dy = yix - yi
         dz = zix - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c  Z-Bisect scheme for local frame of the lone pair
c
      else if (i.eq.3) then
         dx = xiz - xi
         dy = yiz - yi
         dz = ziz - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = xix - xi
         dy = yix - yi
         dz = zix - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = xiy - xi
         dy = yiy - yi
         dz = ziy - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx = dx / r
         dy = dy / r
         dz = dz / r
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     3-Fold method rotation matrix elements for z- and x-axes
c
        else if (i.eq.4) then
         dx = xiz - xi
         dy = yiz - yi
         dz = ziz - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = xix - xi
         dy = yix - yi
         dz = zix - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = xiy - xi
         dy = yiy - yi
         dz = ziy - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx3 = dx / r
         dy3 = dy / r
         dz3 = dz / r
         dx = dx1 + dx2 + dx3
         dy = dy1 + dy2 + dy3
         dz = dz1 + dz2 + dz3
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      end if
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
c
      return
      end
c
c
c     subroutine sibfaparam : hardcode some sibfa parameters (for ctransfer, repulsion and dispersion)
c
      subroutine sibfaparam
      use chargetransfer
      use repulsion
      use kct
      implicit none
      integer i
c
c       tas, tap 
c 
      call aclear(100*100,tas)
      call aclear(100*100,tap)
      tas(8,1) = 0.3488d0
      tap(8,1) = 0.5303d0
      tas(1,8) = 0.3488d0
      tap(1,8) = 0.5303d0
      tas(16,1) = 0.3488
      tap(16,1) = 0.5303 
      tas(34,1) = 0.3488
      tap(34,1) = 0.5303 
      tas(7,1) = 0.435
      tap(7,1) = 0.635 
      tas(9,1) = 0.25
      tap(9,1) = 0.25 
      tas(6,1) = 0.52
      tap(6,1) = 0.52 
      tas(6,7) = 0.5
      tap(6,7) = 0.5 
      tas(6,8) = 0.5
      tap(6,8) = 0.5d0
      tas(6,16) = 0.5d0
      tap(6,16) = 0.5 
      tas(8,8) = 0.5d0
      tap(8,8) = 0.5d0
      tas(16,16) = 0.5
      tap(16,16) = 0.5 
      tas(34,34) = 0.5
      tap(34,34) = 0.5 
      tas(7,7) = 0.5
      tap(7,7) = 0.5 
      tas(6,6) = 0.5
      tap(6,6) = 0.5 
      tas(8,7) = 0.415
      tap(8,7) = 0.48
      tas(9,6) = 0.5
      tap(9,6) = 0.5 
      tas(9,7) = 0.5
      tap(9,7) = 0.5 
      tas(9,8) = 0.5
      tap(9,8) = 0.5 
      tas(9,16) = 0.5
      tap(9,16) = 0.5 
      tas(9,9) = 0.5
      tap(9,9) = 0.5 
      tas(6,9) = 0.5
      tap(6,9) = 0.5 
      tas(7,9) = 0.5
      tap(7,9) = 0.5 
      tas(8,9) = 0.5
      tap(8,9) = 0.5 
      tas(16,9) = 0.5
      tap(16,9) = 0.5 
      tas(34,9) = 0.5
      tap(34,9) = 0.5 
c
c     ma
c       
      call aclear(100*100,ma)
c
      do i = 1, 100
        ma(24,i) = 1.0d0
        ma(i,24) = 1.0d0
      end do
c
      ma(8,1) = 1.03d0
      ma(1,8) = 1.03d0
      ma(34,1) = 1.03
      ma(34,8) = 1.49
      ma(34,7) = 1.35
      ma(34,16) = 1.49
      ma(34,7) = 1.35
      ma(34,6) = 1.2
      ma(34,3) = 0.29
      ma(34,11) = 0.22
      ma(34,12) = 0.29
      ma(34,19) = 0.23
      ma(34,20) = 0.6875
      ma(34,30) = 1.87
      ma(34,29) = 1.87
      ma(34,48) = 1.87
      ma(34,20) = 1.87
      ma(16,1) = 1.03
      ma(7,1) = 1.410
      ma(6,6) = 1.410
      ma(7,7) = 1.45
      ma(8,8) = 1.49d0
      ma(16,16) = 1.49
      ma(8,16) = 1.49
      ma(16,8) = 1.49
      ma(7,8) = 1.53
      ma(8,7) = 1.35
      ma(8,6) = 1.2
      ma(7,16) = 1.53
      ma(16,7) = 1.35
      ma(16,6) = 1.2
      ma(6,16) = 1.59
      ma(7,6) = 1.3
      ma(6,3) = 0.29
      ma(8,3) = 0.29
      ma(16,3) = 0.29
      ma(7,3) = 0.35
      ma(8,11) = 0.22
      ma(16,11) = 0.22
      ma(7,11) = 0.275
      ma(8,12) = 0.29
      ma(16,12) = 0.29
      ma(7,12) = 0.37
      ma(6,12) = 0.47
      ma(8,19) = 0.23
      ma(16,19) = 0.23
      ma(7,19) = 0.21
      ma(8,20) = 0.3
      ma(16,20) = 0.6875
      ma(7,20) = 0.34
      ma(8,30) = 0.66
      ma(7,30) = 0.53
      ma(8,48) = 0.66
      ma(7,48) = 0.53
      ma(6,48) = 0.46
      ma(6,1) = 1.23
      ma(6,7) = 1.51
      ma(6,8) = 1.59
      ma(6,11) = 0.32
      ma(6,12) = 0.46
      ma(6,19) = 0.19
      ma(6,20) = 0.38
      ma(6,30) = 0.46
c
c     forb
c       
      call aclear(100*100,forb)
c
      do i = 1, 100
        forb(24,i) = 1.0
        forb(i,24) = 1.0
      end do
c
      forb(34,1) = 1.03
      forb(34,8) = 1.49
      forb(34,7) = 1.35
      forb(34,16) = 1.49
      forb(34,7) = 1.35
      forb(34,6) = 1.2
      forb(34,3) = 0.29
      forb(34,11) = 0.22
      forb(34,12) = 0.29
      forb(34,19) = 0.23
      forb(34,20) = 0.6875
      forb(34,30) = 1.87
      forb(34,29) = 1.87
      forb(34,48) = 1.87
      forb(34,20) = 1.87
      forb(16,1) = 1.03
      forb(7,1) = 1.410
      forb(6,6) = 1.410
      forb(7,7) = 1.45
      forb(8,8) = 1.49d0
      forb(8,1) = 1.03d0
      forb(16,16) = 1.49
      forb(8,16) = 1.49
      forb(16,8) = 1.49
      forb(7,8) = 1.53
      forb(8,7) = 1.35
      forb(8,6) = 1.2
      forb(7,16) = 1.53
      forb(16,7) = 1.35
      forb(16,6) = 1.2
      forb(6,16) = 1.59
      forb(7,6) = 1.3
      forb(6,3) = 0.29
      forb(8,3) = 0.29
      forb(16,3) = 0.29
      forb(7,3) = 0.35
      forb(8,11) = 0.22
      forb(16,11) = 0.22
      forb(7,11) = 0.275
      forb(8,12) = 0.29
      forb(16,12) = 0.29
      forb(7,12) = 0.37
      forb(6,12) = 0.47
      forb(8,19) = 0.23
      forb(16,19) = 0.23
      forb(7,19) = 0.21
      forb(8,20) = 0.3
      forb(16,20) = 0.6875
      forb(7,20) = 0.34
      forb(8,30) = 0.66
      forb(7,30) = 0.53
      forb(8,48) = 0.66
      forb(7,48) = 0.53
      forb(6,48) = 0.46
      forb(6,1) = 1.23
      forb(6,7) = 1.51
      forb(6,8) = 1.59
      forb(6,11) = 0.32
      forb(6,12) = 0.46
      forb(6,19) = 0.19
      forb(6,20) = 0.38
      forb(6,30) = 0.46
c
      return
      end
c
c     subroutine rotlp1: computes the product of the derivatives of the
c     rotation matrixes with the non rotated lp vectors
c     
      subroutine rotlp1(ilp,di,dix,diz)
      use chargetransfer
      implicit none
      integer ilp,l,m
      real*8 di(3,3),dix(3,3),diz(3,3)
      real*8 rot(3,3),dri(3,3,3),drix(3,3,3),driz(3,3,3)
c
      call derrotlp(ilp,.true.,lpatom(ilp),izlp(ilp),
     $  ixlp(ilp),rot,dri,driz,drix)
c
      call aclear(9,di)
      call aclear(9,dix)
      call aclear(9,diz)
c
      do l = 1,3
        do m = 1,3
          di(1,l) = di(1,l) + dri(1,l,m)*lonepair(m,ilp)
          di(2,l) = di(2,l) + dri(2,l,m)*lonepair(m,ilp)
          di(3,l) = di(3,l) + dri(3,l,m)*lonepair(m,ilp)
          dix(1,l) = dix(1,l) + drix(1,l,m)*lonepair(m,ilp)
          dix(2,l) = dix(2,l) + drix(2,l,m)*lonepair(m,ilp)
          dix(3,l) = dix(3,l) + drix(3,l,m)*lonepair(m,ilp)
          diz(1,l) = diz(1,l) + driz(1,l,m)*lonepair(m,ilp)
          diz(2,l) = diz(2,l) + driz(2,l,m)*lonepair(m,ilp)
          diz(3,l) = diz(3,l) + driz(3,l,m)*lonepair(m,ilp)
        end do
      end do
c
      di(1,1) = di(1,1)+1
      di(2,2) = di(2,2)+1
      di(3,3) = di(3,3)+1
      return
      end
