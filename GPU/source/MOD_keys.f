c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module keys  --  contents of current keyword parameter file  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     nkey      number of nonblank lines in the keyword file
c     keyline   contents of each individual keyword file line
c
c
      module keys
      use sizes
      implicit none
      integer nkey
      character*240 keyline(maxkey)

      interface fetchkey
         module procedure fetchkeya
         module procedure fetchkeyl
         module procedure fetchkeyi
         module procedure fetchkeyi8
         module procedure fetchkeyr4
         module procedure fetchkeyr8
      end interface

      contains

      logical function existkey( keyword,index )
      implicit none
      character(*) ,intent(in ):: keyword
      integer      ,intent(out),optional::index
      character(64) keyw,record
      integer       i,siz

      existkey = .false.
      siz      = len(keyword)
      keyw     = merge(keyword,keyword//' ',keyword(siz:siz).eq.' ')
      call upcase( keyw )
      do i = 1, nkey
         record = keyline(i)(1:siz+1)
         call upcase(record)
         if (record.eq.keyw) then
            existkey = .true.
            if (present(index)) index=i
            return
         end if
      end do
      !print*,'test key::',trim(keyword), i,nkey,' no found'
      if (present(index)) index=nkey+1
      end function

      subroutine fetchkeya(keyword,keyvalue,default)
      implicit none
      character(*)  ,intent(In ):: keyword
      character(240),intent(Out):: keyvalue
      character(*)  ,optional,intent(in):: default
      integer         i,siz,next,temp
      character*240:: record,keyw

      call upcase(keyword)
      siz    = len(keyword)
      keyw   = keyword//" "

      do i = 1, nkey
         record = keyline(i)(1:siz+1)
         call upcase(record)
         if (record.eq.keyw) then
            keyvalue = keyline(i)(siz+1:240)
            return
         end if
      end do
      if (present(default)) keyvalue = default
      end subroutine
      subroutine fetchkeyl(keyword,keyvalue,default,forma)
      implicit none
      character(*) keyword
      logical      keyvalue
      logical     ,optional,intent(in):: default
      character(*),optional,intent(in):: forma
      integer         i,siz,next,temp
      character*240:: record,keyw

      call upcase(keyword)
      siz    = len(keyword)
      keyw   = keyword//" "

      do i = 1, nkey
         record = keyline(i)(1:siz+1)
         call upcase(record)
         if (record.eq.keyw) then
            record = keyline(i)(siz+1:240)
            next   = 1
            call gettext(record,keyw,next)
            if (present(forma)) then
               read (keyw,trim(forma),err=10,end=10) temp
            else
               read (keyw,*,err=10,end=10) temp
            end if
            keyvalue = merge(.true.,.false.,btest(temp,0))
            return
         end if
      end do
      if (present(default)) keyvalue = default
      return
 10   continue
      write(0,'(/,a,a,a,a)')
     &        "fetchkeyl:Issue accessing keyword ",keyword
     &       ,"    ",trim(keyw)
      call fatal
      end subroutine
      subroutine fetchkeyi(keyword,keyvalue,default,forma)
      implicit none
      character(*) keyword
      integer      keyvalue
      integer     ,optional,intent(in):: default
      character(*),optional,intent(in):: forma
      integer         i,siz,next
      character*240:: record,keyw

      call upcase(keyword)
      siz    = len(keyword)
      keyw   = keyword//" "

      do i = 1, nkey
         record = keyline(i)(1:siz+1)
         call upcase(record)
         if (record.eq.keyw) then
            record = keyline(i)(siz+1:240)
            next   = 1
            call gettext(record,keyw,next)
            if (present(forma)) then
               read (keyw,trim(forma),err=10,end=10) keyvalue
            else
               read (keyw,*,err=10,end=10) keyvalue
            end if
            return
         end if
      end do
      if (present(default)) keyvalue = default
      return
 10   continue
      write(0,'(/,a,a,a,a)') 
     &       "fetchkeyi:Issue accessing keyword ",keyword
     &      ,"    ",trim(keyw)
      call fatal
      end subroutine
      subroutine fetchkeyi8(keyword,keyvalue,default,forma)
      implicit none
      character(*) keyword
      integer(8)   keyvalue
      integer(8)  ,optional,intent(in):: default
      character(*),optional,intent(in):: forma
      integer         i,siz,next
      character*240:: record,keyw

      call upcase(keyword)
      siz    = len(keyword)
      keyw   = keyword//" "

      do i = 1, nkey
         record = keyline(i)(1:siz+1)
         call upcase(record)
         if (record.eq.keyw) then
            record = keyline(i)(siz+1:240)
            next   = 1
            call gettext(record,keyw,next)
            if (present(forma)) then
               read (keyw,trim(forma),err=10,end=10) keyvalue
            else
               read (keyw,*,err=10,end=10) keyvalue
            end if
            return
         end if
      end do
      if (present(default)) keyvalue = default
      return
 10   continue
      write(0,'(/,a,a,a,a)')
     &        "fetchkeyi8:Issue accessing keyword ",keyword
     &       ,"    ",trim(keyw)
      call fatal
      end subroutine
      subroutine fetchkeyr4(keyword,keyvalue,default,forma)
      implicit none
      character(*) keyword
      real(4)      keyvalue
      real(4)     ,optional,intent(in):: default
      character(*),optional,intent(in):: forma
      integer         i,siz,next
      character*240:: record,keyw

      call upcase(keyword)
      siz    = len(keyword)
      keyw   = keyword//" "

      do i = 1, nkey
         record = keyline(i)(1:siz+1)
         call upcase(record)
         if (record.eq.keyw) then
            record = keyline(i)(siz+1:240)
            next   = 1
            call gettext(record,keyw,next)
            if (present(forma)) then
               read (keyw,trim(forma),err=10,end=10) keyvalue
            else
               read (keyw,*,err=10,end=10) keyvalue
            end if
            return
         end if
      end do
      if (present(default)) keyvalue = default
      return
 10   continue
      write(0,'(/,a,a,a,a)')
     &        "fetchkeyr4:Issue accessing keyword ",keyword
     &       ,"    ",trim(keyw)
      call fatal
      end subroutine
      subroutine fetchkeyr8(keyword,keyvalue,default,forma)
      implicit none
      character(*) keyword
      real(8)      keyvalue
      real(8)     ,optional,intent(in):: default
      character(*),optional,intent(in):: forma
      integer         i,siz,next
      character*240:: record,keyw

      call upcase(keyword)
      siz    = len(keyword)
      keyw   = keyword//" "

      do i = 1, nkey
         record = keyline(i)(1:siz+1)
         call upcase(record)
         if (record.eq.keyw) then
            record = keyline(i)(siz+1:240)
            next   = 1
            call gettext(record,keyw,next)
            if (present(forma)) then
               read (keyw,trim(forma),err=10,end=10) keyvalue
            else
               read (keyw,*,err=10,end=10) keyvalue
            end if
            return
         end if
      end do
      if (present(default)) keyvalue = default
      return
 10   continue
      write(0,'(/,a,a,a,a)')
     &        "fetchkeyr8:Issue accessing keyword ",keyword
     &       ,"    ",trim(keyw)
      call fatal
      end subroutine
      end module
