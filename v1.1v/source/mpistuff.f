c
c    subroutine commdirdir : deal with communications of direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles 
c        - 2: wait for the communications to be done
c
      subroutine commdirdir(nrhs,rule,mu,reqrec,reqsend)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
      integer :: reqrec(nproc),reqsend(nproc)
      real*8 ef(3,nrhs,npolebloc)
      real*8 mu(3,nrhs,npolebloc)
 1000 format(' illegal rule in commdirdir.')
c
      if (rule.eq.0) then
c
c     MPI : begin reception
c
        do i = 1, n_recep1
          tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   MPI_COMM_WORLD,reqrec(i),ierr)
        end do
      else if (rule.eq.1) then
c
c     MPI : begin sending
c
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
          call MPI_ISEND(mu,3*nrhs*npoleloc,
     $     MPI_REAL8,p_send1(i),tag,MPI_COMM_WORLD,
     $     reqsend(i),ierr)
        end do
      else if (rule.eq.2) then
        do i = 1, n_recep1
           call MPI_WAIT(reqrec(i),status,ierr)
         end do
         do i = 1, n_send1
           call MPI_WAIT(reqsend(i),status,ierr)
         end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c    subroutine commdirdirshort : deal with communications of short range direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles 
c        - 2: wait for the communications to be done
c
cnon  subroutine commdirdirshort(nrhs,rule,mu,reqrec,reqsend)
cnon  use domdec
cnon  use iounit
cnon  use mpole
cnon  use mpi
cnon  implicit none
cnon  integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
cnon  integer :: reqrec(nproc),reqsend(nproc)
cnon  real*8 ef(3,nrhs,npolebloc)
cnon  real*8 mu(3,nrhs,npolebloc)
cnon0 format(' illegal rule in commdirdir.')
cnon
cnon  if (rule.eq.0) then
cnon
cnon  MPI : begin reception
cnon
cnon    do i = 1, n_recepshort1
cnon      tag = nproc*rank + p_recepshort1(i) + 1
cnon    call MPI_IRECV(mu(1,1,bufbegpole(p_recepshort1(i)+1)),
cnon $   3*nrhs*domlen(p_recepshort1(i)+1),MPI_REAL8,p_recepshort1(i),
cnon $   tag,MPI_COMM_WORLD,reqrec(i),ierr)
cnon    end do
cnon  else if (rule.eq.1) then
cnon
cnon  MPI : begin sending
cnon
cnon    do i = 1, n_sendshort1
cnon      tag = nproc*p_sendshort1(i) + rank + 1
cnon      call MPI_ISEND(mu,3*nrhs*npoleloc,
cnon $     MPI_REAL8,p_sendshort1(i),tag,MPI_COMM_WORLD,
cnon $     reqsend(i),ierr)
cnon    end do
cnon  else if (rule.eq.2) then
cnon    do i = 1, n_recepshort1
cnon       call MPI_WAIT(reqrec(i),status,ierr)
cnon     end do
cnon     do i = 1, n_sendshort1
cnon       call MPI_WAIT(reqsend(i),status,ierr)
cnon     end do
cnon  else
cnon    if (rank.eq.0) write(iout,1000) 
cnon    call fatal 
cnon  end if
cnon  return
cnon  end
c
c
c    subroutine commrecdir : deal with communications of  reciprocal fields
c
c    rule determines what to do:
c        - 0: start reception of reciprocal fields
c        - 1: start sending of reciprocal fields
c        - 2: wait for the communications to be done
c
      subroutine commrecdirfields(rule,efrec,ef,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: reqrecdirrec(nproc),reqrecdirsend(nproc)
      real*8 efrec(10,max(1,npolerecloc))
      real*8 ef(10,max(1,npoleloc))
      real*8 buffermpi1(10,max(npoleloc,1))
      real*8 buffermpi2(10,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdirfields.')
c
      if (rule.eq.0) then
c
c     Begin the reception of the reciprocal fields
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpi1(1,bufbeg1(proc+1)),
     $        10*buflen1(proc+1),
     $        MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdirrec(i),ierr)
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            do j = 0, buflen2(proc+1)-1
              call amove(10,efrec(1,polerecloc(
     $       buf2(bufbeg2(proc+1)+j))),buffermpi2(1,bufbeg2(proc+1)+j))
            end do
          end if
        end do
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpi2(1,bufbeg2(proc+1)),
     $       10*buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $       reqrecdirsend(i),ierr)
          end if
        end do
      else if (rule.eq.2) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirsend(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirrec(i),status,ierr)
            do j = 0, buflen1(proc+1)-1
              call amove(10,buffermpi1(1,bufbeg1(proc+1)+j),
     $          ef(1,poleloc(buf1(bufbeg1(proc+1)+j))))
            end do
          end if
        end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c    subroutine commrecdirsolv : deal with communications of reciprocal fields in the solver
c
c    rule determines what to do:
c        - 0: start reception of reciprocal fields
c        - 1: start sending of reciprocal fields
c        - 2: wait for the communications to be done
c
      subroutine commrecdirsolv(nrhs,rule,dipfieldbis,dipfield,
     $ buffermpi1,buffermpi2,reqrecdirrec,reqrecdirsend)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: reqrecdirrec(nproc),reqrecdirsend(nproc)
      real*8 dipfield(3,nrhs,max(1,npoleloc))
      real*8 dipfieldbis(3,nrhs,max(1,npolerecloc))
      real*8 buffermpi1(3,nrhs,max(npoleloc,1))
      real*8 buffermpi2(3,nrhs,max(1,npolerecloc))
 1000 format(' illegal rule in commrecdirdip.')
c
      if (rule.eq.0) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
           tag = nproc*rank + proc + 1
           call MPI_IRECV(buffermpi1(1,1,bufbeg1(proc+1)),
     $       3*nrhs*buflen1(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdirrec(i),ierr)
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,dipfieldbis(1,1,polerecloc(
     $          buf2(bufbeg2(proc+1)+j))),
     $          buffermpi2(1,1,bufbeg2(proc+1)+j))
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpi2(1,1,bufbeg2(proc+1)),
     $       3*nrhs*buflen2(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdirsend(i),ierr)
          end if
        end do
      else if (rule.eq.2) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirsend(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirrec(i),status,ierr)
            do j = 0, buflen1(proc+1)-1
              call amove(3*nrhs,buffermpi1(1,1,bufbeg1(proc+1)+j),
     $          dipfield(1,1,poleloc(buf1(bufbeg1(proc+1)+j))))
            end do
          end if
        end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
c
c    subroutine commrecdirdipsolv : deal with communications of dipoles for PME 
c
c    rule determines what to do:
c        - 0: start reception of dipoles
c        - 1: start sending of dipoles 
c        - 2: wait for the communications to be done
c
      subroutine commrecdirdip(nrhs,rule,diprec,dip,buffermpimu1,
     $ buffermpimu2,req2rec,req2send)
      use domdec
      use iounit
      use mpole
      use mpi
      implicit none
      integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i,j,proc
      integer :: req2rec(nproc),req2send(nproc)
      real*8 dip(3,nrhs,max(1,npoleloc))
      real*8 diprec(3,nrhs,max(1,npolerecloc))
c      real*8 buffermpimu1(3,nrhs,max(npoleloc,1))
c      real*8 buffermpimu2(3,nrhs,max(1,npolerecloc))
      real*8 buffermpimu1(3,nrhs,*)
      real*8 buffermpimu2(3,nrhs,*)
 1000 format(' illegal rule in commrecdirdip.')
c
      if (rule.eq.0) then
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),3*nrhs*
     $       buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $       req2rec(i),ierr)
          end if
        end do
      else if (rule.eq.1) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            do j = 0, buflen1(proc+1)-1
              call amove(3*nrhs,dip(1,1,poleloc(
     $          buf1(bufbeg1(proc+1)+j))),
     $          buffermpimu1(1,1,bufbeg1(proc+1)+j))
            end do
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),3*nrhs*
     $       buflen1(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $       req2send(i),ierr)
          end if
        end do
      else if (rule.eq.2) then
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(req2send(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(req2rec(i),status,ierr)
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $          diprec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
            end do
          end if
        end do
      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if
      return
      end
c
