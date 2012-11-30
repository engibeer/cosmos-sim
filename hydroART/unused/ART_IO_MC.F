c     =============================================================
c
c      ART Version 3: ART_IO.f - routines dealing with Input/Output
c
c     =============================================================

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
# include "a_def.h"
c
      integer iunit
      character*256 filename, accesstype
c
      nfn  = index(filename, ' ') - 1
      nat  = index(accesstype, ' ') - 1
c
c.... take into account different syntax on AIX 
c     (define or comment OS_AIX in a_def.h as needed)
c
#ifdef OS_AIX
      open ( iunit , file = filename(1:nfn), 
     &     position = accesstype(1:nat) )
#endif
#ifndef OS_AIX
      open ( iunit , file = filename(1:nfn), 
     &               access = accesstype(1:nat) )
#endif
c
      return
      end
c
c     ------------------------------
      subroutine Save ( dt , iMode )
c     ------------------------------
c
c     iMode defines how the output file is named
c           = 0 , filename is just the jobname
c           = 1 , filename is jobname + global step number
c           = 2 , filename is jobname + expansion parameter       
c
#     include "a_def.h"
      include 'a_setup.h'
      include 'a_control.h'
      character*256 fname_gas, fname_part1, fname_part2, fname_part3
      character*256 fname_part4, fname_part5
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*5 chstep
      real*8 dt

      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1
c
      if ( iMode .eq. 0 ) then 
        fname_gas  = 
     &  path(1:lpath)//'/'//jobname1(1:ljob)//'.d '
        fname_part1 = 
     &  path(1:lpath)//'/PMcrd.DAT '
        fname_part2 = 
     &  path(1:lpath)//'/PMcrs0.DAT '
        fname_part3 = 
     &  path(1:lpath)//'/pt.dat '
        fname_part4 = 
     &  path(1:lpath)//'/stars.dat '
      endif
      if ( iMode .eq. 1 ) then
        is = istep 
        idummy = is / 10000
        is = is - idummy * 10000
        chstep = digits(idummy)
        ich = 1
        do i = 3 , 0 , -1
          idummy = is / 10**i        
          chstep = chstep(1:ich) // digits(idummy)
          ich = ich + 1
          is = is - idummy*10**i
        enddo
        fname_gas  = 
     &  path(1:lpath)//'/'//jobname1(1:ljob)//'_'//chstep//'.d '
        fname_part1 = 
     &  path(1:lpath)//'/PMcrd_'//chstep//'.DAT '
        fname_part2 = 
     &  path(1:lpath)//'/PMcrs0_'//chstep//'.DAT '
        fname_part3 = 
     &  path(1:lpath)//'/pt_'//chstep//'.dat '
        fname_part4 = 
     &  path(1:lpath)//'/stars_'//chstep//'.dat '
      endif
      if ( iMode .eq. 2 ) then 
          ic = int(aexpn*1000)
          i1 = ic / 1000 
          i2 = (ic - (ic/1000)*1000)/100 
          i3 = (ic - (ic/100)*100) / 10 
          i4 = (ic - (ic/10)*10)

          fname_gas = path(1:lpath)//'/'//jobname1(1:ljob)//'_a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.d '
          fname_part1 = path(1:lpath)//'/PMcrda'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.DAT '
          fname_part2 = path(1:lpath)//'/PMcrs0a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.DAT '
          fname_part3 = path(1:lpath)//'/pta'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
          fname_part4 = path(1:lpath)//'/stars_a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
      endif




      call Write_Gas_Binary ( dt , fname_gas )





#ifdef NBODY
      call Write_Particles_Binary ( fname_part1 , fname_part2 )
      call Write_Particle_Time    ( fname_part3 )
#endif
#ifdef STARFORM
      call Write_Stellar_Data ( fname_part4 )
#endif 



      return
      end
c
c     -----------------------------------
      subroutine Save_Part ( dt , iMode )
c     -----------------------------------
c
c     iMode defines how the output file is named
c           = 0 , filename is just the jobname
c           = 1 , filename is jobname + global step number
c           = 2 , filename is jobname + expansion parameter       
c
#     include "a_def.h"
      include 'a_setup.h'
      include 'a_control.h'
      character*256 fname_part1, fname_part2, fname_part3
      character*256 fname_part4, fname_part5
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*5 chstep
      real*8 dt

      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1
c
      if ( iMode .eq. 0 ) then 
        fname_part1 = 
     &  path(1:lpath)//'/PMcrd.DAT '
        fname_part2 = 
     &  path(1:lpath)//'/PMcrs0.DAT '
        fname_part3 = 
     &  path(1:lpath)//'/pt.dat '
        fname_part4 = 
     &  path(1:lpath)//'/stars.dat '
      endif
      if ( iMode .eq. 1 ) then
        is = istep 
        idummy = is / 10000
        is = is - idummy * 10000
        chstep = digits(idummy)
        ich = 1
        do i = 3 , 0 , -1
          idummy = is / 10**i        
          chstep = chstep(1:ich) // digits(idummy)
          ich = ich + 1
          is = is - idummy*10**i
        enddo
        fname_part1 = 
     &  path(1:lpath)//'/PMcrd_'//chstep//'.DAT '
        fname_part2 = 
     &  path(1:lpath)//'/PMcrs0_'//chstep//'.DAT '
        fname_part3 = 
     &  path(1:lpath)//'/pt_'//chstep//'.dat '
        fname_part4 = 
     &  path(1:lpath)//'/stars_'//chstep//'.dat '
      endif
      if ( iMode .eq. 2 ) then 
          ic = int(aexpn*1000)
          i1 = ic / 1000 
          i2 = (ic - (ic/1000)*1000)/100 
          i3 = (ic - (ic/100)*100) / 10 
          i4 = (ic - (ic/10)*10)

          fname_part1 = path(1:lpath)//'/PMcrda'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.DAT '
          fname_part2 = path(1:lpath)//'/PMcrs0a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.DAT '
          fname_part3 = path(1:lpath)//'/pta'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
          fname_part4 = path(1:lpath)//'/stars_a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
      endif
 
#ifdef NBODY
      call Write_Particles_Binary ( fname_part1 , fname_part2 )
      call Write_Particle_Time    ( fname_part3 )
#endif
#ifdef STARFORM
      call Write_Stellar_Data ( fname_part4 )
#endif 

      return
      end
c
c     ----------------------------
      subroutine Save_Check ( dt )
c     ----------------------------
      include 'a_setup.h'
      include 'a_control.h'

      real*8 dt

      do i = 1 , n_save
        if ( aexpn .gt. asave(i) .and. isave(i) .eq. 0 ) then 
          call Save ( dt , 2 ) 
c          call HaloFind ( 2 )
          isave(i) = 1
        endif
      enddo
      if ( mod(istep,nfsave) .eq. 0 ) then 
        call Save ( dt , 0 ) 
      endif
CEVERINO12052005
         call Save_InternalB( dt)
CEVERINO4/18/2005
c$$$c
c$$$c...  save DM and stellar particles every step (D. Nagai, 03/04/04)
c$$$c
c$$$      if ( mod(istep,2) .eq. 0 .and. aexpn .gt. 0.0625 ) then 
c$$$        call Save_Part ( dt , 2 ) 
c$$$        call HaloFind ( 3 )
c$$$      endif

      return
      end
c
c     ------------------------------------------
      subroutine Write_Gas_Binary ( dt , fname )
c     ------------------------------------------
c
c     purpose: writes files with control information and gas data
c              writes hydrodynamic variables for cells
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      character*256 fname
      real*8 b2a
      real*8 dt
c
      nfn = index ( fname , ' ' ) - 1
      call Get_MaxLevelNow ()
      open  ( 22 , file = fname(1:nfn) , form = 'unformatted' )
      write ( 22 ) jobname1
#ifdef COSMOLOGY
      adum = b2a ( (1.d0*t) )
#endif
#ifndef COSMOLOGY
      adum = 1.0 
#endif
      write ( 22 ) istep , t , dt , adum , ainit
      write ( 22 ) boxh, Om0, Oml0, Omb0, hubble
      write ( 22 ) nextra
      write ( 22 ) (extra(i), i=1,nextra)
      write ( 22 ) (lextra(i), i=1,nextra)
      write ( 22 ) MinLevel , MaxLevelNow
      write ( 22 ) ( tl(i)     , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( dtl(i)    , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( tlold(i)  , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( dtlold(i) , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( iSO(i)    , i = MinLevel , MaxLevelNow )
      write ( 22 ) ncell0
      write ( 22 ) (iOctCh(i), i=1,ncell0)
      write ( 22 ) ((hvar(ivar,i),ivar=1,nhvar), i=1,ncell0)
      write ( 22 ) ((var(ivar,i),ivar=2,3), i=1,ncell0)
      if ( MaxLevelNow .gt. MinLevel ) then 
        write ( 22 ) iOctFree , nOct
        do Level = MinLevel+1 , MaxLevelNow
          nLevel = iNOLL(Level)      
          write ( 22 ) Level , iNOLL(Level) , iHOLL(Level)
          iOct = iHOLL(Level)
          do ic1 = 1 , nLevel
            write(22) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct), 
     &                iOctLL2(iOct)
            iOct = iOctLL1(iOct)
          enddo
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel 
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idc = icell + ic2
              write(22) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar), 
     &                 (var(i,idc), i=2,3)
            enddo
          enddo
        enddo
      endif
      close (22)
      write(*,*) 'done writing...'
c
      return
      end
c
c     --------------------------------------------
      subroutine Read_Gas_Binary ( dt , FileName )
c     --------------------------------------------
c
c     purpose: reads files with control information and gas data
c              and hydrodynamic variables for cells
c
      include 'a_tree.h'
      include 'a_control.h'
      real*8 dt
      character*256 FileName
      character*256 jname
      real*8 b2a
c
      nfname = index ( FileName , ' ' ) - 1
      open  ( 19 , file = FileName(1:nfname) , form = 'unformatted' )
      write(*,*) 'reading tree data from file ',FileName(1:nfname)
      read ( 19 ) jname 
      if ( jname .ne. jobname1 ) then 
        write(*,*) '* warning : Read_Gas_Binary : jname .ne. jobname1'
        write(*,'(2(A256))') jname, jobname1
c        stop
      endif
      read ( 19 ) istep , t , dt, aexpn, ainit  
      write(*,20) 
     &      istep, t, dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      do Level = MinLevel , MaxLevel 
        aexp(Level) = aexpn
      enddo
c     if having trouble on linux, try inserting a dummy first in the read
c     list below
      read ( 19 ) boxh, Om0, Oml0, Omb0, hubble
      write(*,21) 
     &      boxh, Om0, Oml0, Omb0, hubble      
 21   format('Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f5.3,1x,
     &        ' Oml0= ',f5.3,1x,' Omb0= ',f7.5,1x,' h = ',f5.3)
      read ( 19 ) nextras
      if ( nextras .gt. nextra ) then 
        write(*,*) '* error: Read_Gas_Binary : nextras > nextra :',
     &             nextras, nextra
        write(*,*) '* (check nextra parameter in a_control.h) '
        stop
      endif
      read ( 19 ) (extra(i), i=1,nextras)
      read ( 19 ) (lextra(i), i=1,nextras)

      read ( 19 ) MinLev , MaxLevelNow
      if ( MaxLevelNow .gt. MaxLevel ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : MaxLevelNow > MaxLevel :',
     &        MaxLevelNow, MaxLevel 
         stop
      endif
      read ( 19 ) ( tl(i)     , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtl(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( tlold(i)  , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtlold(i) , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( iSO(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ncell
      if ( ncell .ne. ncell0 ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell not equal ncell0 :',
     &        ncell, ncell0
         stop        
      endif
      read ( 19 ) (iOctCh(i), i=1,ncell)
      read ( 19 ) ((hvar(ivar,i),ivar=1,nhvar), i=1,ncell)
      read ( 19 ) ((var(ivar,i),ivar=2,3), i=1,ncell)
      ntot = 0 
      if ( MaxLevelNow .gt. MinLevel ) then 
        read ( 19 ) iOctFree , nOct
        write(*,*) 'iOctFree =',iOctFree,' nOct =', nOct
        if ( ncell+nchild*nOct .gt. mcell ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell+nchild*nOct > mcell :',
     &        ncell+nchild*nOct, mcell 
         write(*,*)  
     &    '* (size of the input exceeds max. of cells (mcell))'
         write(*,*) '* increase mcell in a_setup.h and try again...'
         stop
        endif 

        do Lev = MinLevel+1 , MaxLevelNow
          read ( 19 ) Level , iNOLL(Level) , iHOLL(Level)
          iOct   = iHOLL(Level)
          nLevel = iNOLL(Level)
          nLevCells = nLevel*nchild
          if ( nLevCells .gt. nclmax ) then 
            write(*,*) 'error : L =',Lev,' nLevCell =',nLevCell,
     &        ' > nclmax =',nclmax,' set in a_setup.h'
            write(*,*) '=> increase nclmax and rerun.'
            close ( 19 )
            stop
          endif
          ntot = ntot + nLevel 
          write(*,*) 'reading tree (oct) data for level ', Lev
          do ic1 = 1 , nLevel
            read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct), 
     &                iOctLL2(iOct)
            iOct = iOctLL1(iOct)
          enddo
          write(*,*) 'reading cell data: ncells =', nLevel*nchild
          do ic1 = 1 , nLevel*nchild
            read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar), 
     &                 (var(i,idc), i=2,3)
          enddo
        enddo
      endif
      close (19)
      write(*,*) 'done reading',ncell0+nchild*ntot,' cells...'

c     set previous time moment for the energy cons. control
      ap0    = b2a ( (1.d0*(t-5.d-1*dt)) )      
c
      return
      end
c
c     ------------------------------------------------------
      subroutine Read_Gas_Binary_Adiabatic ( dt , FileName )
c     ------------------------------------------------------
c
c     purpose: reads files with control information and gas data
c              and hydrodynamic variables for cells
c
c     Only difference from Read_Gas_Binary () is 
c         nhvar - > nhvar-nchem
c
      include 'a_tree.h'
      include 'a_control.h'
      real*8 dt
c     Don't forget to change this to real*8
      character*256 FileName
      character*256 jname
      real*8 b2a
c
      nfname = index ( FileName , ' ' ) - 1
      open  ( 19 , file = FileName(1:nfname) , form = 'unformatted' )
      write(*,*) 'reading tree data from file ',FileName(1:nfname)
      read ( 19 ) jname 
      if ( jname .ne. jobname1 ) then 
        write(*,*) '* warning : Read_Gas_Binary : jname .ne. jobname1'
        write(*,'(2(A256))') jname, jobname1
c        stop
      endif
      read ( 19 ) istep , t , dt, aexpn, ainit  
c     quick fixes to match types 
c      t = t1
c      dt = dt1
      write(*,20) 
     &      istep, t, dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      do Level = MinLevel , MaxLevel 
        aexp(Level) = aexpn
      enddo
c     if having trouble on linux, try inserting a dummy first in the read
c     list below
      read ( 19 ) boxh, Om0, Oml0, Omb0, hubble
      write(*,21) 
     &      boxh, Om0, Oml0, Omb0, hubble      
 21   format('Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f5.3,1x,
     &        ' Oml0= ',f5.3,1x,' Omb0= ',f7.5,1x,' h = ',f5.3)
      read ( 19 ) nextras
      if ( nextras .gt. nextra ) then 
        write(*,*) '* error: Read_Gas_Binary : nextras > nextra :',
     &             nextras, nextra
        write(*,*) '* (check nextra parameter in a_control.h) '
        stop
      endif
      read ( 19 ) (extra(i), i=1,nextras)
      read ( 19 ) (lextra(i), i=1,nextras)

      read ( 19 ) MinLev , MaxLevelNow
      if ( MaxLevelNow .gt. MaxLevel ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : MaxLevelNow > MaxLevel :',
     &        MaxLevelNow, MaxLevel 
         stop
      endif
      read ( 19 ) ( tl(i)     , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtl(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( tlold(i)  , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtlold(i) , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( iSO(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ncell
      if ( ncell .ne. ncell0 ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell not equal ncell0 :',
     &        ncell, ncell0
         stop        
      endif
      read ( 19 ) (iOctCh(i), i=1,ncell)
      read ( 19 ) ((hvar(ivar,i),ivar=1,nhvar-nchem), i=1,ncell)
      read ( 19 ) ((var(ivar,i),ivar=2,3), i=1,ncell)
      ntot = 0 
      if ( MaxLevelNow .gt. MinLevel ) then 
        read ( 19 ) iOctFree , nOct
        write(*,*) 'iOctFree =',iOctFree,' nOct =', nOct
        if ( ncell+nchild*nOct .gt. mcell ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell+nchild*nOct > mcell :',
     &        ncell+nchild*nOct, mcell 
         write(*,*)  
     &    '* (size of the input exceeds max. of cells (mcell))'
         write(*,*) '* increase mcell in a_setup.h and try again...'
         stop
        endif 

        do Lev = MinLevel+1 , MaxLevelNow
          read ( 19 ) Level , iNOLL(Level) , iHOLL(Level)
          iOct   = iHOLL(Level)
          nLevel = iNOLL(Level)
          nLevCells = nLevel*nchild
          if ( nLevCells .gt. nclmax ) then 
            write(*,*) 'error : L =',Lev,' nLevCell =',nLevCell,
     &        ' > nclmax =',nclmax,' set in a_setup.h'
            write(*,*) '=> increase nclmax and rerun.'
            close ( 19 )
            stop
          endif
          ntot = ntot + nLevel 
          write(*,*) 'reading tree (oct) data for level ', Lev
          do ic1 = 1 , nLevel
            read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct), 
     &                iOctLL2(iOct)
            iOct = iOctLL1(iOct)
          enddo
          write(*,*) 'reading cell data: ncells =', nLevel*nchild
          do ic1 = 1 , nLevel*nchild
            read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar-nchem), 
     &                 (var(i,idc), i=2,3)
          enddo
        enddo
      endif
      close (19)
      write(*,*) 'done reading',ncell0+nchild*ntot,' cells...'

c     set previous time moment for the energy cons. control
      ap0    = b2a ( (1.d0*(t-5.d-1*dt)) )      
c
      return
      end

c     --------------------------------------------
      subroutine Read_Gas_Binary_old ( dt , FileName )
c     --------------------------------------------
c
c     purpose: reads files with control information and gas data
c              and hydrodynamic variables for cells
c
      include 'a_tree.h'
      include 'a_control.h'
      real*8 dt
c     Don't forget to change this to real*8
      real dt1, t1
      real dtl1(MinLevel:MaxLevel),tl1(MinLevel:MaxLevel)
      real dtlold1(MinLevel:MaxLevel),tlold1(MinLevel:MaxLevel)
      character*256 FileName
      character*256 jname
      real*8 b2a
c
      nfname = index ( FileName , ' ' ) - 1
      open  ( 19 , file = FileName(1:nfname) , form = 'unformatted' )
      write(*,*) 'reading tree data from file ',FileName(1:nfname)
      read ( 19 ) jname 
      if ( jname .ne. jobname1 ) then 
        write(*,*) '* warning : Read_Gas_Binary : jname .ne. jobname1'
        write(*,'(2(A256))') jname, jobname1
c        stop
      endif
      read ( 19 ) istep , t1 , dt1, aexpn, ainit  
c     quick fixes to match types 
      t = t1
      dt = dt1
      write(*,20) 
     &      istep, t, dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      do Level = MinLevel , MaxLevel 
        aexp(Level) = aexpn
      enddo
c     if having trouble on linux, try inserting a dummy first in the read
c     list below
      read ( 19 ) boxh, Om0, Oml0, Omb0, hubble
      write(*,21) 
     &      boxh, Om0, Oml0, Omb0, hubble      
 21   format('Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f5.3,1x,
     &        ' Oml0= ',f5.3,1x,' Omb0= ',f7.5,1x,' h = ',f5.3)
      read ( 19 ) nextras
      if ( nextras .gt. nextra ) then 
        write(*,*) '* error: Read_Gas_Binary : nextras > nextra :',
     &             nextras, nextra
        write(*,*) '* (check nextra parameter in a_control.h) '
        stop
      endif
      read ( 19 ) (extra(i), i=1,nextras)
      read ( 19 ) (lextra(i), i=1,nextras)

      read ( 19 ) MinLev , MaxLevelNow
      if ( MaxLevelNow .gt. MaxLevel ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : MaxLevelNow > MaxLevel :',
     &        MaxLevelNow, MaxLevel 
         stop
      endif
      read ( 19 ) ( tl1(i)     , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtl1(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( tlold1(i)  , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtlold1(i) , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( iSO(i)    , i = MinLevel , MaxLevelNow )
c
      do l = MinLevel, MaxLevelNow
        tl(l) = tl1(l)
        dtl(l) = dtl1(l)
        tlold(l) = tlold1(l)
        dtlold(l) = dtlold1(l)
      enddo
c

      read ( 19 ) ncell
      if ( ncell .ne. ncell0 ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell not equal ncell0 :',
     &        ncell, ncell0
         stop        
      endif
      read ( 19 ) (iOctCh(i), i=1,ncell)
      read ( 19 ) ((hvar(ivar,i),ivar=1,nhvar-nchem), i=1,ncell)
      read ( 19 ) ((var(ivar,i),ivar=2,3), i=1,ncell)
      ntot = 0 
      if ( MaxLevelNow .gt. MinLevel ) then 
        read ( 19 ) iOctFree , nOct
        write(*,*) 'iOctFree =',iOctFree,' nOct =', nOct
        if ( ncell+nchild*nOct .gt. mcell ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell+nchild*nOct > mcell :',
     &        ncell+nchild*nOct, mcell 
         write(*,*)  
     &    '* (size of the input exceeds max. of cells (mcell))'
         write(*,*) '* increase mcell in a_setup.h and try again...'
         stop
        endif 

        do Lev = MinLevel+1 , MaxLevelNow
          read ( 19 ) Level , iNOLL(Level) , iHOLL(Level)
          iOct   = iHOLL(Level)
          nLevel = iNOLL(Level)
          nLevCells = nLevel*nchild
          if ( nLevCells .gt. nclmax ) then 
            write(*,*) 'error : L =',Lev,' nLevCell =',nLevCell,
     &        ' > nclmax =',nclmax,' set in a_setup.h'
            write(*,*) '=> increase nclmax and rerun.'
            close ( 19 )
            stop
          endif
          ntot = ntot + nLevel 
          write(*,*) 'reading tree (oct) data for level ', Lev
          do ic1 = 1 , nLevel
            read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct), 
     &                iOctLL2(iOct)
            iOct = iOctLL1(iOct)
          enddo
          write(*,*) 'reading cell data: ncells =', nLevel*nchild
          do ic1 = 1 , nLevel*nchild
            read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar-nchem), 
     &                 (var(i,idc), i=2,3)
          enddo
        enddo
      endif
      close (19)
      write(*,*) 'done reading',ncell0+nchild*ntot,' cells...'

c     set previous time moment for the energy cons. control
      ap0    = b2a ( (1.d0*(t-5.d-1*dt)) )      
c
      return
      end
c
c     ----------------------------------------------------------
      subroutine Read_Particles_Binary ( FileName1 , FileName2 )
c     ----------------------------------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
      character*256 FileName1 , FileName2
c 
      ngrid = ng

      nfn1 = index ( FileName1 , ' ' ) - 1
      nfn2 = index ( FileName2 , ' ' ) - 1
c
      open ( 3 , file = FileName1(1:nfn1) , form = 'unformatted' )

c.... read control information and check whether it has proper structure

      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     &                  Ocurv,Omb0,extras
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRIDC,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &           1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &           1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &           1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &           1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &               ' Omegab_0=',F7.3,' Hubble=',f7.3)
      
c...  Working with new ICs (Daisuke added a line below) 
      if( nspecies .eq. 0 ) nspecies = 1
      if ( nspecies .gt. nspec ) then
        write (*,*)
     &      ' nspec in PARAMETER < nspecies in TAPE-FILE:'
         write (*,*) ' nspecies=',nspecies,' > ',nspec
         STOP
      endif
      If( npmax .lt. lspecies(nspecies) ) then
         write (*,*) ' Wrong number of particles: '
         write (*,*) ' should be =',lspecies(nspecies),' (lspecies)'
         write (*,*) ' but is set to ',npmax,' in a_setup.h...'
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo
         STOP
      Endif 

      nbyte  = nrecl * 4
      nacces = nbyte
#ifdef NBYTEWORD
      nacces = nrecl
#endif

      open ( 1 , file = FileName2(1:nfn2), access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3


      N_particles = lspecies(nspecies)
      np = N_particles 

      write(*,*) ' hello :',N_particles, np

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      do i = 1 , nspecies
        write(iOWork,*) 'ispecie =',i,lspecies(i)
      enddo
      write(iOWork,*) 'Read_Particles_Binary: reading np=',np,
     &                ' particles'
      close(iOWork)
#endif

      if ( np .eq. 0 .or. np .gt. npmax ) then 
        call Open_ASCII_File ( iOErr , errorfile , append ) 
        write(iOErr,*) '1 : error in Read_Particles_Binary: bad np =',
     &             np
        write(iOErr,*) 'stopping...'
        close ( iOErr )
        stop
      endif

      Npages      = (N_particles -1)/npage + 1
      N_in_last   = N_particles - npage*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',nspecies
      write (*,*) ' N_in_last=',N_in_last

      do irow = 1 , Npages         ! loop over particle pages
        In_page = npage
        if ( irow .eq. Npages ) In_page = N_in_last
c         write (*,*)' Read page=',IROW,' file=',ifile,' N=',In_page
        iL = npage * (irow-1)
        CALL GetRow(irow,1) ! read in a page of particles
        do in = 1 , In_page          ! Loop over particles
          ip = in + iL                     ! current particle number
          x(ip) = xpar(in)
          y(ip) = ypar(in)
          z(ip) = zpar(in)
          vx(ip) = vxx(in)
          vy(ip) = vyy(in)
          vz(ip) = vzz(in)
        enddo
      enddo

      close (1)
      close (3)

      do i = 1 , nspecies
        wpar(i) = wspecies(i)
        lsp(i)  = lspecies(i)
      enddo
      do i = 1 , nspecies
        if(i .eq. 1 ) then
          nsp(i,1) = 1
        else
           nsp(i,1) = lsp(i-1) + 1
        endif 
        nsp(i,2) = lsp(i)
      enddo 

C$OMP PARALLEL DO 
      do ic1 = 1 , np 
        pt(ic1) = t   ! set t to the current time moment
      enddo
c
c.... compute each particle's weight
c
      n_spec = nspecies
#ifdef STARFORM
      n_spec = nspecies-1
#endif
      do ic0 = 1 , n_spec 
        do ic1 = nsp(ic0,1) , nsp(ic0,2)
           pw(ic1) = wpar(ic0)
        enddo
      enddo 

#ifdef STARFORM
      if ( start ) then         
        if ( nspecies+1 .gt. nspec ) then
          write (*,*)
     &        ' nspec in PARAMETER < nspecies+1 for starformation:'
          write (*,*) ' nspecies+1=',nspecies+1,' > nspec =',nspec
          STOP
        endif
        nspecies = nspecies + 1
        nsp(nspecies,1) = lsp(nspecies-1) + 1
        nsp(nspecies,2) = lsp(nspecies-1)
        lsp(nspecies) = nsp(nspecies,2)
        lspecies(nspecies) = lsp(nspecies)


        write(*,*) 'in reading particles:'
        write(*,*) '    lsp      wpar         nsp1      nsp2 '
        do i = 1 , nspecies
           write(*,30) lsp(i), wpar(i), nsp(i,1), nsp(i,2)
        enddo
        write(*,*) ' '
      endif
#endif
      
 30   format (i10,2x,g11.5,2x,2(i8,2x))

      return
      end
c
c     ------------------------------------------------------------
      subroutine Write_Particles_Binary ( FileName1 , FileName2 )
c     ------------------------------------------------------------
c
c     purpose: writes control information and particles to the specified 
c              files (this routine is to be used when nspec = 1)
c
c     input  : FileName1 - C3CRD*; FileName2 - C3crs0*
c

      include 'a_tree.h'
      include 'a_control.h'
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))

      character*256 FileName1, FileName2

      nbyte  = nrecl * 4
      nacces = nbyte
#ifdef NBYTEWORD
      nacces = nrecl                                               
#endif

      nfn1 = index ( FileName1 , ' ' ) - 1
      nfn2 = index ( FileName2 , ' ' ) - 1

      open (4 , file = FileName1(1:nfn1),
     &           form = 'UNFORMATTED' , status = 'UNKNOWN')
 
      open (5 , file = FileName2(1:nfn2) , access = 'DIRECT',
     &	         status = 'UNKNOWN', recl = NACCES)

c.... write header and control data

      write (4) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                  ,Ocurv,Omb0,extras
      close (4)

      lspecies(nspecies) = lsp(nspecies)
      N_particles = lspecies(nspecies)   ! Total number of particles
      Npages = (N_particles -1)/npage +1
      N_in_last = N_particles - npage*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',nspecies
      write (*,*) ' N_in_last=',N_in_last

      do irow = 1 , Npages         ! Loop over particle pages
        In_page = npage
        If ( irow .eq. Npages ) In_page = N_in_last
c         write (*,*)' Write page=',IROW,' file=',ifile,' N=',In_page
        iL = npage * (irow-1)
        do in = 1 , In_page          ! Loop over particles
          ip = in + iL                     ! current particle number
	  xpar(in) = x(ip) 
	  ypar(in) = y(ip) 
	  zpar(in) = z(ip) 
	   vxx(in) = vx(ip)
	   vyy(in) = vy(ip)
	   vzz(in) = vz(ip)
         enddo
         call WriRow ( irow , 5 )
      enddo

      close (5)
      
      return
      end

c     ----------------------------------
      subroutine WriRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'a_tree.h'
      integer irow , ifile 
      write (ifile , rec = irow ) recdat
      return
      end
 
c     ----------------------------------
      subroutine GetRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'a_tree.h'
      read  ( ifile , rec = irow ) recdat
      return
      end

c     
c     ------------------------------------------
      subroutine Read_Particle_Time ( FileName )
c     ------------------------------------------
c
c     purpose: when the code starts or restarts
c              this routine takes care of array pt
c 
      include 'a_tree.h'
      include 'a_control.h'
      character*256 FileName 
c
      nfn = index ( FileName , ' ' ) - 1
      open ( 60 ,file = FileName(1:nfn),
     &       form = 'unformatted',status = 'old' )
      read (60) (pdt(ic1),ic1=1,np)
      close( 60 )
c
      return
      end

c     -------------------------------------------
      subroutine Write_Particle_Time ( FileName )
c     -------------------------------------------
c
c     purpose: write particle time moments in file FileName
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      character*256 FileName
      nfn = index ( FileName , ' ' ) - 1
#ifdef DEBUG
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'In Write_Particle_Time: np=',np
      close(iOWork)      
#endif
      open ( 60 , file = FileName(1:nfn) , form = 'unformatted' )
      write (60) (pdt(ic1),ic1=1,np)
      close ( 60 )
      return
      end
c     
c     ------------------------------------------
      subroutine Read_Stellar_Data ( FileName )
c     ------------------------------------------
c
c     purpose: when the code starts or restarts
c              this routine reads all data on stellar particles 
c 
      include 'a_tree.h'
      include 'a_control.h'
      character*256 FileName 
      real*8 tdum, adum 
      real*8 ws_old, ws_oldi
      common / s_old1 / ws_old, ws_oldi
c
      nfn = index ( FileName , ' ' ) - 1
      open ( 60 ,file = FileName(1:nfn),
     &       form = 'unformatted',status = 'old' )
      read(60) tdum, adum
      read(60) nstars

      if ( nstars .ne. nsp(nspecies,2) - nsp(nspecies,1) + 1 ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*) '1  : In Read_Stellar_Data: something is wrong:'
        write(13,*) 'nstars =',nstars,' is iconsistent with'
        write(13,*) 'nsp(nspecies,2)-nsp(nspecies,1)+1 =',
     &               nsp(nspecies,2) - nsp(nspecies,1) + 1
        write(13,*) 'stopping...'
        close(13)              
        stop
      endif
      if ( nstars .eq. 0 ) return
      if ( nstars .lt. 0 .or. nstars .gt. nstarmax ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*) '1  : In Read_Stellar_Data: bad nstars =',nstars
        write(13,*) 'must be 0 < nstars < nstarmax =',nstarmax        
        write(13,*) 'stopping...'
        close(13)              
        stop
      endif
      read (60) ws_old, ws_oldi
      read (60) (pw(ic1),ic1=nsp(nspecies,1),nsp(nspecies,2))   ! weights
      read (60) (pw0(ic1),ic1=1,nstars)      ! initial masses
      read (60) (tbirth(ic1),ic1=1,nstars)   ! birth times
#ifdef ENRICH
      read (60) (zstII(ic1),ic1=1,nstars)    ! metallicity 
#endif
#ifdef ENRICH_SNIa
      read (60) (zstIa(ic1),ic1=1,nstars)    ! metallicity 
#endif
      close( 60 )
c
      return
      end

c     -------------------------------------------
      subroutine Write_Stellar_Data ( FileName )
c     -------------------------------------------
c
c     purpose: write all stellar data to a file FileName 
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      real*8 ws_old, ws_oldi
      common / s_old1 / ws_old, ws_oldi

      character*256 FileName
      nfn = index ( FileName , ' ' ) - 1

      open ( 60 , file = FileName(1:nfn) , form = 'unformatted' )
      nstars = nsp(nspecies,2)-nsp(nspecies,1)+1
      write(60) t, aexp(0)
      write(60) nstars 
      write(60) ws_old, ws_oldi   ! current and init. tot. mass in stars
      write(60) (pw(ic1),ic1=nsp(nspecies,1),nsp(nspecies,2))   ! weights
      write(60) (pw0(ic1),ic1=1,nstars)      ! initial masses
      write(60) (tbirth(ic1),ic1=1,nstars)   ! birth times
#ifdef ENRICH
      write(60) (zstII(ic1),ic1=1,nstars)    ! metallicity 
#endif
#ifdef ENRICH_SNIa
      write(60) (zstIa(ic1),ic1=1,nstars)    ! metallicity 
#endif
      close ( 60 )
      return
      end
c--------------------------------------------------------------------
      SUBROUTINE Save_InternalB( dt)
c-------------------------------------------------------------------- 
c
c     This subroutine store all the information of a step in memory.
c     It is designed to work on  Altix superclusters.
c     Daniel Ceverino (2005)
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_buffer.h'
      real*8 b2a
      real*8 dt
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
      real*8 ws_old, ws_oldi
      common / s_old1 / ws_old, ws_oldi
#ifdef HYDRO
c
c     Gas data
c
      iB0(1) = istep
      iB0(2) = MaxLevelNow
      ainitbu = ainit  
      do i=MinLevel,MaxLevelNow
         ISObu(i) = iSO(i)
         tlbu(i) = tl(I)
         dtlbu(i) = dtl(i)
         tloldbu(i) = tlold(i)
         dtloldbu(i) = dtlold(i)
      enddo
      B0(1) = t
      B0(2) = dt
      B0(3) = adum
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO i=1,ncell0
         iB1(i) = iOctCh(i)
      enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(jcell,ivar)
      DO i=1,ncell0
         jcell = (i-1)*n0
         B1(jcell+1) = var(2,i)
         B1(jcell+2) = var(3,i)
         jcell = jcell + 2
         do ivar=1,nhvar
            B1(jcell+ivar) = hvar(ivar,i)
         enddo
      ENDDO
      if ( MaxLevelNow .gt. MinLevel ) then   
         iB0(3) = iOctFree
         iB0(4) = nOct
         do Level = MinLevel+1 , MaxLevelNow
            nLevel = iNOLL(Level)  
            iOct = iHOLL(Level)
            iB2(Level) = nLevel
            iB3(Level) = iOct
         enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(nLevel,iOct,ic1,jOct,i)
         DO Level = MinLevel+1 , MaxLevelNow
            nLevel = iNOLL(Level)            
            iOct = iHOLL(Level)
            do ic1 = 1 , nLevel
               jOct =  iOct*13 
               do i=1,3
                  iB4(Level, jOct+i)=iOctPs(i,iOct)
               enddo
               jOct = jOct+3
               do i=1, 6
                  iB4(Level, jOct+i)=iOctNb(i,iOct)
               enddo
               jOct = jOct + 6
               iB4(Level, jOct +1) = iOctPr(iOct)
               iB4(Level,jOct + 2) = iOctLv(iOct)
               iB4(Level,jOct + 3) = iOctLL1(iOct)
               iB4(Level,jOct + 4) = iOctLL2(iOct)
               iOct = iOctLL1(iOct)
            enddo
         ENDDO
         do Level = MinLevel+1 , MaxLevelNow
            nLevel = iNOLL(Level)
            call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell, ic2,idc,jcell  )            
            DO ic1 = 1 , nLevel 
               icell = iSelect(ic1)
               do ic2 = 0 , 7
                  idc = icell + ic2
                  jcell = (idc-1)*n0
                  iB5(idc) = iOctCh(idc)
                  do i=1,nhvar
                     B2(jcell+i) = hvar(i,idc)
                  enddo
                  jcell=jcell+nhvar
                  B2(jcell+1) = var(2,idc)
                  B2(jcell+2) = var(3,idc)
               enddo
            ENDDO
         enddo
      endif
#endif
#ifdef NBODY
c
c     Particles data
c
      lspecies(nspecies) = lsp(nspecies)
      N_particles = lsp(nspecies)
      AexpnB = Aexpn
      AstepB = Astep
      istepB = istep
      do i=1,100
         extrasB(i) = extras(i)
      enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO ip=1, N_particles
         Xbuff(ip) = X(ip)
         Ybuff(ip) = Y(ip)
         Zbuff(ip) = Z(ip)
         Vxbuff(ip) = Vx(ip)
         Vybuff(ip) = Vy(ip)
         Vzbuff(ip) = Vz(ip)
      ENDDO
c
c     Particles time
c
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO ip=1,np
         pdtbuff(ip) = pdt(ip)
      ENDDO
#endif
#ifdef STARFORM
c
c     stars data
c
      nstars = nsp(nspecies,2)-nsp(nspecies,1)+1
      ws_oldB = ws_old
      ws_oldiB= ws_oldi
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO i=1,nstars
         pwsB(i) = pw(nsp(nspecies,1)-1+i)
         pw0B(i) = pw0(i)
         tbirthB(i) = tbirth(i)
#ifdef ENRICH
         zstIIB(i) = zstII(i)
#endif
#ifdef ENRICH_SNIa
         zstIaB(i) = zstIa(i)
#endif
      ENDDO
#endif
	if (.not.internalB) internalB=.TRUE.
      return
      end
c--------------------------------------------------------------------
      SUBROUTINE Read_InternalB( dt)
c-------------------------------------------------------------------- 
c
c     This subroutine reads all the step information from memory.
c     It is designed to work on  Altix superclusters.
c     Daniel Ceverino (2005)
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_buffer.h'
      real*8 b2a
      real*8 dt
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
      real*8 ws_old, ws_oldi
      common / s_old1 / ws_old, ws_oldi
c
c     Gas Data
c
#ifdef  HYDRO
      istep = iB0(1)
      MaxLevelNow=iB0(2)
      ainit = ainitbu  
      do i=MinLevel,MaxLevelNow
         iSO(i)=iSObu(i)
         tl(i)=tlbu(i)
         dtl(i)=dtlbu(i)
         tlold(i)=tloldbu(i)
         dtlold(i)=dtloldbu(i)
      enddo
      t =B0(1) 
      dt=B0(2)
      adum=B0(3)
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO i=1,ncell0
         iOctCh(i) =iB1(i)
      enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(jcell,ivar)
      DO i=1,ncell0
         jcell = (i-1)*n0
         var(2,i)=B1(jcell+1)
         var(3,i)=B1(jcell+2)
         jcell = jcell + 2
         do ivar=1,nhvar
            hvar(ivar,i)=B1(jcell+ivar)
         enddo
      ENDDO
      if ( MaxLevelNow .gt. MinLevel ) then  
         iOctFree=iB0(3)
         nOct = iB0(4)
         do Level = MinLevel+1 , MaxLevelNow
            iNOLL(Level)=iB2(Level)
            iHOLL(Level)=iB3(Level)
         enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(nLevel,iOct,ic1,jOct,i)
         DO Level = MinLevel+1 , MaxLevelNow
            nLevel = iNOLL(Level)            
            iOct = iHOLL(Level)
            do ic1 = 1 , nLevel
               jOct = iOct*13 
               do i=1,3
                 iOctPs(i,iOct)= iB4(Level, jOct+i)
               enddo
               jOct = jOct+3
               do i=1, 6
                  iOctNb(i,iOct)=iB4(Level, jOct+i)
               enddo
               jOct = jOct + 6
               iOctPr(iOct)=iB4(Level, jOct +1)
               iOctLv(iOct) = iB4(Level,jOct + 2)
               iOctLL1(iOct)=iB4(Level,jOct + 3)
               iOctLL2(iOct)=iB4(Level,jOct + 4)
               iOct = iOctLL1(iOct)
            enddo
         ENDDO
         do Level = MinLevel+1 , MaxLevelNow
            nLevel = iNOLL(Level)
            call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell, ic2,idc,jcell  )            
            DO ic1 = 1 , nLevel 
               icell = iSelect(ic1)
               do ic2 = 0 , 7
                  idc = icell + ic2
                  jcell = (idc-1)*n0
                  iOctCh(idc)=iB5(idc)
                  do i=1,nhvar
                     hvar(i,idc)=B2(jcell+i)
                  enddo
                  jcell=jcell+nhvar
                  var(2,idc)=B2(jcell+1)
                  var(3,idc)=B2(jcell+2)
               enddo
            ENDDO
         enddo
      endif
#endif
#ifdef NBODY
c
c     Particles data
c
      Aexpn = AexpnB
      Astep = AstepB
      istep = istepB
      do i=1,100
         extras(i) = extrasB(i)
      enddo
      do i = 1 , nspecies
        wpar(i) = wspecies(i)
        lsp(i)  = lspecies(i)
      enddo
      do i = 1 , nspecies
        if(i .eq. 1 ) then
          nsp(i,1) = 1
        else
           nsp(i,1) = lsp(i-1) + 1
        endif 
        nsp(i,2) = lsp(i)
      enddo 
      N_particles = lsp(nspecies)
      np=N_particles
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO ip=1, N_particles
         X(ip) = Xbuff(ip)
         Y(ip) = Ybuff(ip)
         Z(ip) = Zbuff(ip)
         Vx(ip) = Vxbuff(ip)
         Vy(ip) = Vybuff(ip)
         Vz(ip) = Vzbuff(ip)
      ENDDO
c
c     Particles time
c
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO ip=1,np
         pdt(ip) = pdtbuff(ip)
      ENDDO
#endif
#ifdef STARFORM
c
c     stars data
c
      nstars = nsp(nspecies,2)-nsp(nspecies,1)+1
      ws_old = ws_oldB
      ws_oldi= ws_oldiB
C$OMP PARALLEL DO DEFAULT(SHARED)
      DO i=1,nstars
         pw(nsp(nspecies,1) -1+i) = pwsB(i)
         pw0(i) = pw0B(i)
         tbirth(i) = tbirthB(i)
#ifdef ENRICH
         zstII(i) = zstIIB(i)
#endif
#ifdef ENRICH_SNIa
         zstIa(i) = zstIaB(i)
#endif
      ENDDO
#endif
      return
      end     
c-------------------------------------------------------------------- 
                  
