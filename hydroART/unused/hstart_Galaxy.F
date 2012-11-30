c-------------------------------------------------------------------- 
c
c     Initial conditions for a galaxy simulation
c     from an already set files of particles coordinates
c
c     based on hstart.f ( Anatoly Klypin and Andrey Kravtsov (2001))
c     with the introduction of a initially refined mesh for the gas component.
c
c     Daniel Ceverino (2005)
c
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'hstart_Galaxy.h'
      character*256 fname,fname_part1, fname_part2, fname_part3
      real*8 dt

      call ReadControl () ! starting, reading control parameters ...'
      call Init_Parameters       ()
      call Init_Arrays           ()
      call Init_Tree             ()

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*)' Arrays initialized sucessfully'
      close(iOWork)
#endif

      CALL  init_Galaxy_Run(dt)

      write(*,*) 'writing hydro data...'
!      write (*,*)  'writting is off'
!      STOP
      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1
      fname =
     &     path(1:lpath)//'/'//jobname1(1:ljob)//'.d '
      CALL Write_Gas_Binary ( dt , fname )
      write(*,*) 'writing particles data...'   
      lpath = index(path    , ' ') - 1
      fname_part1 =
     &     path(1:lpath)//'/PMcrd.DAT '
      fname_part2 =
     &     path(1:lpath)//'/PMcrs0.DAT ' 

      write (*,*)  ' specie      weight (after baryons-Genesis)'
      do ispec = 1 , Nspecies
        wpar(ispec) = (1. - Omb0/Om0) * wpar(ispec)
        wspecies(ispec) =  wpar(ispec)
        write (*,*)   ispec , wpar(ispec) 
      enddo
      call Write_Particles_Binary ( fname_part1, fname_part2)
      write (*,*)  'writting particles time...'
      fname_part3 =path(1:lpath)//'/pt.dat ' 
      nfn3 = index ( fname_part3 , ' ' ) - 1
      open ( 60 , file = fname_part3(1:nfn3) , form = 'unformatted' )
      write(60) (pdt(ic1),ic1=1,np)
      close(60)
      END


C-------------------------------------------------------------------- 
      SUBROUTINE init_Galaxy_Run(dt)
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'hstart_Galaxy.h'
C
      real*8 tcorr
      real*8 ALPHA
      real*8 dt, dtold
      real*8 a2b
      double precision xmin, xmax, ymin, ymax, zmin, zmax
      double precision vxmin,vxmax,vymin,vymax,vzmin,vzmax
      Character  Hd*5,Tail*4
      character*256 FileName1 , FileName2

      Hd  ='PMcrs'
      Tail='.DAT'
     
!      read(*,'A256') Filename1  ! It should end with an empty space
!      read(*,'A256') Filename2  ! It should end with an empty space
      FileName1 ='PMcrdIC.DAT '
      FileName2='PMcrs0IC.DAT '
      call Read_Particles_Binary (FileName1 , FileName2)
      call LL_Construct ()     
      If(Nspecies.ne.0.and.2**(Nspecies-1).ne.Lblock)Then
         write (*,*) ' Wrong number of Levels for mass refinement (',
     &        Nspecies,') Block size =',Lblock,
     &        '       2**(Levels-1) must = Block size'
         Stop
      Endif 

      If(Nspecies.eq.0)Then
         np=lsp(1)
      Else
         np =lsp(Nspecies)
      EndIf
c      tcorr = 2./sqrt(Om0) * 0.5 ! *0.5 is to compensate for the NGRID change
!      tcorr = 1.0 
      tcorr =  2./sqrt(Om0)
      xmin = 1.e6
      xmax = -xmin
      ymin = 1.e6
      ymax = -ymin
      zmin = 1.e6
      zmax = -zmin
      vxmin =  1.e6
      vxmax = -vxmin
      vymin =  vxmin
      vymax = -vymin
      vzmin =  vxmin
      vzmax = -vzmin
      
      do ip = 1 , np
c        Xpt(ip) = 0.5 * (Xpt(ip) - 1.0) + 1.
c        Ypt(ip) = 0.5 * (Ypt(ip) - 1.0) + 1.
c        Zpt(ip) = 0.5 * (Zpt(ip) - 1.0) + 1.
        xmin = min(xmin,X(ip))
        xmax = max(xmax,X(ip))
        ymin = min(ymin,Y(ip))
        ymax = max(ymax,Y(ip))
        zmin = min(zmin,Z(ip))
        zmax = max(zmax,Z(ip))
        Vx(ip) = Vx(ip) * tcorr 
        Vy(ip) = Vy(ip) * tcorr
        Vz(ip) = Vz(ip) * tcorr
        vxmin   = min(vxmin,vx(ip) )
        vxmax   = max(vxmax,vx(ip))
        vymin   = min(vymin,vy(ip))
        vymax   = max(vymax,vy(ip))
        vzmin   = min(vzmin,vz(ip))
        vzmax   = max(vzmax,vz(ip))
c        pw(ip) = pw(ip)* 0.125 ! compensate for NGRID change        
      enddo
      write(*,*) 'rescaled particle coordinate range: tcorr= ',tcorr 
      write(*,*) 'xmin,max =',xmin,xmax
      write(*,*) 'ymin,max =',ymin,ymax
      write(*,*) 'zmin,max =',zmin,zmax
      write(*,*) 'vxmin,max   =',vxmin,vxmax
      write(*,*) 'vymin,max   =',vymin,vymax
      write(*,*) 'vzmin,max   =',vzmin,vzmax
     
c
c.... compute each particle's weight
c
      wsum=0.
      do ic0 = 1 , nspecies
        do ic1 = nsp(ic0,1) , nsp(ic0,2)
           pw(ic1) = wpar(ic0)
           wsum = wsum + pw(ic1)
        enddo
      enddo 
        write (*,*) 'N_particles= ',np, ' N_species= ', nspecies
        write (*,*)  ' Total mass[code units]= ', wsum
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
        write (iOwork,*) 'N_particles= ',np, ' N_species= ', nspecies
        write (iOwork,*)  ' Total mass[code units]= ', wsum
      close(iOWork)
#endif
c
c...Setting particles time:
c
           if ( aexpn .lt. astep ) then 
             write(*,*) '*error in Init_Run: aexpn < astep:',aexpn,astep
             stop
           endif
           t = a2b ( (1.d0*aexpn) ) 
           tprev = a2b ( (1.d0*(aexpn-astep)) )
           do ic1 = 1 , np
             pt(ic1) = t
             pdt(ic1) = t - tprev
           enddo
#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'aexpn =',aexpn,' aexp0 =',aexp0,
     &        ' astep =',astep
         WRITE (*,*)   ' t =', t , ' tprev= ', tprev
         close(iOWork)
#endif
c
c time
c
      ainit = aexp0 
      if ( aexpn .lt. astep ) then 
         write(*,*) '*error in Init_Run: aexpn < astep:',aexpn,astep
         stop
      endif
      AEXPV = AEXPN - ASTEP/2.
      fact  =  sqrt(Om0+Oml0*AEXPV**3+Ocurv*AEXPV)
      vfact =  sqrt(Om0+Oml0*AEXPN**3+Ocurv*AEXPN)
c
c.... correct for the difference in t0 units in N-body and Hydro
c     and for the difference in time position of particle and gas vels. 
c
      VCORR =  AEXPN*SQRT(AEXPN)/(AEXPV*SQRT(AEXPV))*vfact / fact
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'Velocity correction= ', VCORR
      close(iOWork)
#endif      
c
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'Set basic conversions units and SF init...'
      close(iOWork)
#endif  
      call SetUnits ()
#ifdef COOLING
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'Set Cooling...'
      close(iOWork)
#endif 
         call Set_Cooling ()
#endif
#ifdef REFINE
      const =  (1.0 - Omb0/Om0)
      constb = Omb0/Om0
      Ndiff=0. 

c$$$      trho(MinLevel) = const  ! Refinement thresholds
c$$$      trho(MinLevel+1) = 0.0001 * const
c$$$      trho(MinLevel+2) = 0.0001 * const
c$$$      tbmass(MinLevel) = 1.0*constb
c$$$      tbmass(MinLevel+1) = 1.0*constb
c$$$      tbmass(MinLevel+2) = 1.0*constb
c$$$      do Level = 3 , MaxLevel 
c$$$        trho(Level)  = const / 32  
c$$$        tbmass(Level) = constb / 16.0 
c$$$      enddo

c$$$      trho(MinLevel) = 6. * const
c$$$      trho(MinLevel+1) = 4. * const
c$$$      trho(MinLevel+2) = 4. * const
c$$$      tbmass(MinLevel) = 4. *constb
c$$$      tbmass(MinLevel+1) = 4. *constb
c$$$      tbmass(MinLevel+2) = 4. *constb
c$$$      do Level = 3 , MaxLevel 
c$$$        trho(Level)  =4 * const 
c$$$        tbmass(Level) = 4*constb  
c$$$      enddo

CEVERINO04302008: FOR NGRID=512
c$$$      trho(MinLevel) =   (wpar(1)*6.) /wsplit
c$$$      trho(MinLevel+1) =(wpar(1)*6.)/wsplit
c$$$      trho(MinLevel+2) = (wpar(1)*10.)/wsplit
c$$$      tbmass(MinLevel) = 10. *constb
c$$$      tbmass(MinLevel+1) = 10. *constb
c$$$      tbmass(MinLevel+2) = 10. *constb
c$$$      do Level = 3 , MaxLevel   
c$$$        trho(Level)  =wpar(1)/wsplit
c$$$        tbmass(Level) = 10.*constb  
c$$$      enddo
c$$$      do Level=MinLevel, MaxLevel
c$$$         write(*,*) 'L =',Level,' trho(L) =',trho(Level),
c$$$     &        ' trhog(L)=',tbmass(Level)
c$$$      enddo

CEVERINO03172009: FOR NGRID=1024
C      trho(MinLevel) =   (wpar(1)*6.) /wsplit
C      trho(MinLevel+1) =(wpar(1)*6.)/wsplit
C      trho(MinLevel+2) = (wpar(1)*6.)/wsplit
C      tbmass(MinLevel) = 10. *constb
C      tbmass(MinLevel+1) = 10. *constb
C      tbmass(MinLevel+2) = 10. *constb
C      tbmass(MinLevel+3) = 10. *constb
C      trho(MinLevel+3) = (wpar(1)*10.)/wsplit
C      do Level = 4 , MaxLevel   
C        trho(Level)  =wpar(1)/wsplit
C        tbmass(Level) = 10.*constb  
C      enddo
C      do Level=MinLevel, MaxLevel
C         write(*,*) 'L =',Level,' trho(L) =',trho(Level),
C     &        ' trhog(L)=',tbmass(Level)
C      enddo

      trho(MinLevel) =   (wpar(1)*6.) /wsplit
      trho(MinLevel+1) =(wpar(1)*6.)/wsplit
      trho(MinLevel+2) = (wpar(1)*6.)/wsplit
      trho(MinLevel+3) = (wpar(1)*6.)/wsplit
      trho(MinLevel+4) = (wpar(1)*10.)/wsplit
c
      tbmass(MinLevel) = 10. *constb
      tbmass(MinLevel+1) = 10. *constb
      tbmass(MinLevel+2) = 10. *constb
      tbmass(MinLevel+3) = 10. *constb
      tbmass(MinLevel+4) = 10. *constb
      do Level = 5 , MaxLevel
        trho(Level)  =wpar(1)/wsplit
        tbmass(Level) = 10.*constb
      enddo
      do Level=MinLevel, MaxLevel
         write(*,*) 'L =',Level,' trho(L) =',trho(Level),
     &        ' trhog(L)=',tbmass(Level)
      enddo


#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*)   'setting refinement thresholds...'
      do Level=MinLevel, MaxLevel
         write(iOWork,*) 'L =',Level,' trho(L) =',trho(Level),
     &        ' trhog(L)=',tbmass(Level)
      enddo
      close(iOWork)
#endif
c
c...  Refine meshes until the refinement converges. (Daisuke 4/3/03)
c
c      write(*,*)
c      write(*,*) ' Refining meshes immediately down to 2nd level'
      call Get_MaxLevelNow ()
      call Assign_Density_AND_GAS ( MinLevel, MaxLevelNow, 0 ,VCORR)
      call Modify ( MinLevel, 1, mtot )
      write(*,*) 'Level refined = ',MinLevel,' mtot = ',mtot
      Level = 1
      do while ( mtot .ne. nil ) 
         call Assign_Density ( Level, Level, 0)
         call Modify ( Level, 1, mtot )
         write(*,*) 'Level refined = ',Level,' mtot = ',mtot
         Level = Level + 1
      enddo
      call Assign_Density_AND_GAS ( MinLevel, Level-1, 0 ,VCORR)   
#endif  

c
c     BLACK MAGIC and setting pressure:
c
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'BLACK MAGIC and setting pressure... '
      close(iOWork)
#endif
      call Get_MaxLevelNow ()
      do L = MinLevel , MaxLevelNow         
         call Magic ( L )
         call TreeEoS ( L ) 
      enddo   
c
c Set level time:
c
      do Level = MinLevel , MaxLevel 
         aexp(Level) = aexpn
      enddo
c
c TIME STEP
c
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'choosing Time Step...'
      close(iOWork)
#endif

      call ChooseTimeStep ( dt )


#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*)' t= ', t , ' dt =', dt
      close(iOWork)
#endif
      do Level=MinLevel, MaxLevel
         tl(Level)=t
         dtl(Level) = dt / 2**iTimeBin(Level)
         tlold(Level) = t
         dtlold(Level) = dtl(Level)
         iSO (Level) =0
      enddo
      END

c
c     ----------------------------------------------------------
      subroutine Read_Particles_Binary ( FileName1 , FileName2 )
c     ----------------------------------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
c     ONLY DM particles.
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'hstart_Galaxy.h'
c
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
     &                  Ocurv,extras
      Omb0 = 0.022/(hubble)**2
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

      boxh = extras(100)
      N_particles = lspecies(nspecies)
      np = N_particles 

      write(*,*) ' hello :',N_particles, np
      write (*,*)  'hola: ', boxh
c 

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'Read_Particles_Binary: reading np=',np,
     &                ' particles'
      do i = 1 , nspecies
       write(iOWork,*) 'ispecie =',i,lspecies(i),' weight= ',wspecies(i)
      enddo
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

c... Write 
        write(*,*) 'in reading particles:'
        write(*,*) '    lsp      wpar         nsp1      nsp2 '
        do i = 1 , nspecies
           write(*,30) lsp(i), wpar(i), nsp(i,1), nsp(i,2)
        enddo
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*)  '    lsp      wpar         nsp1      nsp2 '
      do i = 1 , nspecies
         write(iOWork,30)  lsp(i), wpar(i), nsp(i,1), nsp(i,2)
      enddo
      close(iOWork)
#endif
 30   format (i10,2x,g11.5,2x,2(i8,2x))
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
c     ----------------------------------------------------------
      subroutine Read_Particles_Binary1 ( FileName1 , FileName2 )
c     ----------------------------------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
c     DM + stars
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
CEVERINO08312005: It was behind '.... compute each particle's weight' originally
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


      return
      end
c
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
c$$$#ifdef COSMOLOGY
c$$$      adum = b2a ( (1.d0*t) )
c$$$#endif
c$$$#ifndef COSMOLOGY
c$$$      adum = 1.0 
c$$$#endif
      write ( 22 ) istep , t , dt , aexpn , ainit
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
c     ------------------------------------------------------------
      subroutine Write_Particles_Binary1 ( FileName1 , FileName2 )
c     ------------------------------------------------------------
c
c
c     input  : FileName1 - C3CRD*; FileName2 - C3crs0*
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'hstart_Galaxy.h'

      character*15 FileName1
      character*16 FileName2

      nbyte  = nrecl * 4
      nacces = nbyte / nbyteword

      open (4 , file = FileName1,
     &           form = 'UNFORMATTED' , status = 'UNKNOWN')

      open (5 , file = FileName2 , access = 'DIRECT',
     &           status = 'UNKNOWN', recl = NACCES)

c.... write header and control data

      write (4) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRID,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                  ,Ocurv,Omb0,extras
      close (4)

      N_particles = lsp(nspecies)   ! Total number of particles
      Npages = (N_particles -1)/npage +1
      N_in_last = N_particles - npage*(Npages-1)
      write (*,*) 'N_particles =',N_particles
      write(*,*) ' Pages=',Npages,' Species=',nspecies
      write (*,*) ' N_in_last=',N_in_last

      do irow = 1 , Npages         ! Loop over particle pages
        In_page = npage
        If ( irow .eq. Npages ) In_page = N_in_last
c         write (*,*)' Write page=',IROW,' file=',ifile,' N=',In_page
        iL = npage * (irow-1)
        do in = 1 , In_page          ! Loop over particles
          ip = in + iL                     ! current particle number
          xpar(in) = X(ip)
          ypar(in) = Y(ip)
          zpar(in) = Z(ip)
           vxx(in) = Vx(ip)
           vyy(in) = Vy(ip)
           vzz(in) = Vz(ip)
         enddo
         call WriRow ( irow , 5 )
      enddo

      close (5)

      return
      end
c     ----------------------------------
      subroutine WriRow1 ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      integer irow , ifile 
      write (ifile , rec = irow ) recdat
      return
      end
c
c     ------------------------------------------------------------
      subroutine Write_Particles_Binary ( FileName1 , FileName2 )
c     ------------------------------------------------------------
c
c     purpose: writes control information and particles to the specified 
c              files
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
 
c
c     --------------------------------------------------------------
      double precision function INTEGRATE(FUNC,fp,np,a,b,dxinit,eps)
c     --------------------------------------------------------------
c
c     Quadrature using fifth order Runge-Kutta with adaptive step size.
c     Based on Press et al, Numerical Recipes in C, 2nd ed, pp 719-722.
c
c     Runge-Kutta driver with adaptive stepsize control.  Integrate starting
c     value y from a to b with accuracy eps, storing intermediate results in
c     global variables.  dxinit should be set as a guessed first stepsize.
c
c     Pass a vector of parameters of length np to FUNC in fp(np).
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC

      implicit none
      integer np 
      double precision a, b, eps, dxinit, FUNC, fp(np)
      external FUNC

      integer maxsteps
      parameter(maxsteps=100000000)

      double precision x, dx, dxnext, y, dydx, yscale
      integer  Nstep

      x     = a
      dx    = dxinit
      y     = 0.d0
      Nstep = 0

      do while ((x-b)*(b-a).lt.0.d0.and.Nstep.lt.maxsteps)
        Nstep = Nstep + 1
        dydx = FUNC(x,fp,np)
c
c       yscale is the scaling used to monitor accuracy.  This general-purpose
c       choice can be modified if need be.
c
        yscale = max(abs(y) + abs(dx*dydx), 1.d-20)
        if ((x+dx-b)*(x+dx-a).gt.0.d0)  ! If stepsize overshoots, decrease it.
     1    dx = b - x

        call RUNGE5VAR(y,dydx,x,dx,eps,yscale,dxnext,FUNC,fp,np)

        dx = dxnext
      end do

      if (Nstep.ge.maxsteps)
     1  write (*,*) 'WARNING: failed to converge in INTEGRATE.'

      INTEGRATE = y

      return
      end
c
c     -------------------------------------------------------------
      SUBROUTINE RUNGE5VAR(y,dydx,x,htry,eps,yscale,hnext,DERIVS,
     1                     fp,np)
c     -------------------------------------------------------------
c
c     Fifth-order Runge-Kutta step with monitoring of local truncation error
c     to ensure accuracy and adjust stepsize.  Input are the dependent
c     variable y and its derivative dydx at the starting value of the
c     independent variable x.  Also input are the stepsize to be attempted
c     htry, the required accuracy eps, and the value yscale, against which the
c     error is scaled.  On output, y and x are replaced by their new values.
c     hdid is the stepsize that was actually accomplished, and hnext is the
c     estimated next stepsize.  DERIVS is the user-supplied routine that
c     computes right-hand-side derivatives.  The argument fp is a vector 
c     of parameters (np parameters) to be passed to DERIVS 
c     (NOT integrated over).
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none
      integer np 
      double precision eps,hnext,htry,x,dydx,y,yscale,DERIVS,fp(np)
      external DERIVS

      double precision errmax,h,hold,htemp,xnew,yerr,ytemp

      double precision safety,pgrow,pshrink,errcon
      parameter (safety  =  0.9d0)
      parameter (pgrow   = -0.2d0)
      parameter (pshrink = -0.25d0)
      parameter (errcon  =  1.89d-4)

      h = htry                         ! Set stepsize to initial accuracy.
      errmax = 10.d0
      do while (errmax.gt.1.d0)
        call RUNGE(y,dydx,x,h,ytemp,yerr,DERIVS,fp,np)

        errmax = abs(yerr/yscale)/eps   ! Scale relative to required accuracy.
        if (errmax.gt.1.d0) then        ! Truncation error too large; reduce h
          htemp = safety*h*(errmax**pshrink)
          hold = h
          h = sign(max(abs(htemp),0.1d0*abs(h)),h)  ! No more than factor of 10
          xnew = x + h
          if (xnew.eq.x) then
            write (*,*) 'WARNING: ',
     1                  'Stepsize underflow in RUNGE5VAR().'
            h = hold
            errmax = 0.d0
          end if
        end if
      end do
c
c     Step succeeded.  Compute estimated size of next step.
c
      if (errmax.gt.errcon) then
        hnext = safety*h*(errmax**pgrow)
      else
        hnext = 5.d0 * h                ! No more than factor of 5 increase.
      end if
      x = x + h

      y = ytemp

      return
      end
c
c     ---------------------------------------------------
      SUBROUTINE RUNGE(y,dydx,x,h,yout,yerr,DERIVS,fp,np)
c     ---------------------------------------------------
c
c     Given values for a variable y and its derivative dydx known at x, use
c     the fifth-order Cash-Karp Runge-Kutta method to advance the solution
c     over an interval h and return the incremented variables as yout.  Also
c     return an estimate of the local truncation error in yout using the
c     embedded fourth order method.  The user supplies the routine
c     DERIVS(x,y,dydx), which returns derivatives dydx at x.
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none

      integer np
      double precision h,x,dydx,y,yerr,yout,DERIVS,fp(np)

      external DERIVS

      double precision ak3, ak4, ak5 ,ak6

      double precision a2,a3,a4,a5,a6
      double precision c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      parameter(a2  =    0.2d0)
      parameter(a3  =    0.3d0)
      parameter(a4  =    0.6d0)
      parameter(a5  =    1.d0)
      parameter(a6  =    0.875d0)
      parameter(c1  =   37.d0/378.d0)
      parameter(c3  =  250.d0/621.d0)
      parameter(c4  =  125.d0/594.d0)
      parameter(c6  =  512.d0/1771.d0)
      parameter(dc1 = c1 -  2825.d0/27648.d0)
      parameter(dc3 = c3 - 18575.d0/48384.d0)
      parameter(dc4 = c4 - 13525.d0/55296.d0)
      parameter(dc5 = -277.d0/14336.d0)
      parameter(dc6 = c6 -     0.25d0)

      ak3 = DERIVS(x+a3*h,fp,np)
      ak4 = DERIVS(x+a4*h,fp,np)
      ak5 = DERIVS(x+a5*h,fp,np)
      ak6 = DERIVS(x+a6*h,fp,np)
c
c     Estimate the fifth order value.
c
      yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6)
c
c     Estimate error as difference between fourth and fifth order
c
      yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

      return
      end
c
c     -------------------------------------
      FUNCTION zbrent(func,fp,np,x1,x2,tol)
c     -------------------------------------
c
      INTEGER ITMAX
      REAL zbrent
      integer np 
      real*8 fp(np)
      real*8 tol,x1,x2,EPS
      real*8 func
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b
      REAL c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a,fp,np)
      fb=func(b,fp,np)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then 
        write(*,*) 'root must be bracketed for zbrent'
        write(*,*) 'fa =',fa,'  fb =', fb
        pause
      endif
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,fp,np)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
c     -----------------
      FUNCTION RANDd(M)
c     -----------------
c     random number generator
c     initialize with arbitrary integer > 0
c     after that never change it
c     -----------------
      DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4
     +/453815927,2147483648.,2147483647,536870912,131072,256,16777216,
     +   4,16384,8388608,128/
      ML=M/K1*K1
      M1=(M-ML)*L1
      ML=M/K2*K2
      M2=(M-ML)*L2
      ML=M/K3*K3
      M3=(M-ML)*L3
      ML=M/K4*K4
      M4=(M-ML)*L4
      M5=KI-M
      IF(M1.GE.M5)M1=M1-KI-1
      ML=M+M1
      M5=KI-ML
      IF(M2.GE.M5)M2=M2-KI-1
      ML=ML+M2
      M5=KI-ML
      IF(M3.GE.M5)M3=M3-KI-1
      ML=ML+M3
      M5=KI-ML
      IF(M4.GE.M5)M4=M4-KI-1
      ML=ML+M4
      M5=KI-ML
      IF(LC.GE.M5)ML=ML-KI-1
      M=ML+LC
      RANDd=M/AM
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
