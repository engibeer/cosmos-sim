      character*256 filename, accesstype 
      character*256 rewind
      parameter ( rewind = 'rewind ' )
      parameter ( filename='test_OPEN.dat ')
      iunit=13
      accesstype = 'SEQUENTIAL'
c      call Open_ASCII_File ( 13 , filename, accesstype )

      a=1.d9
      write (*,*)  a
      end

c     ------------------------------------------------------------
      subroutine Open_ASCII_File ( iunit , filename , accesstype )
c     ------------------------------------------------------------
c
c     filename is the string containing the name of the file to open
c              (should end with the space character)
c     accesstype is the string containing the type of file access 
c                (e.g., append)
c              (should end with the space character)
c

c
      integer iunit
      character*256 filename, accesstype
c
c      write (*,*)  filename
c      write (*,'(a)')  accesstype
      nfn  = index(filename, ' ') - 1
      nat  = index(accesstype, ' ') - 1
c      write (*,*)  filename(1:nfn)
c      write (*,*)  accesstype(1:nat)

c$$$      nat  = index(accesstype, ' ') - 1
c
c.... take into account different syntax on AIX 
c     (define or comment OS_AIX in a_def.h as needed)
c
c$$$#ifdef OS_AIX
c$$$      open ( iunit , file = filename(1:nfn), 
c$$$     &     position = accesstype(1:nat) )
c$$$#endif
c$$$#ifndef OS_AIX
c$$$      open ( iunit , file = filename(1:nfn), 
c$$$     &               access = accesstype(1:nat) )
c$$$#endif


      open ( iunit , file = filename(1:nfn), 
     &     access = accesstype(1:nat) )



c
      return
      end
