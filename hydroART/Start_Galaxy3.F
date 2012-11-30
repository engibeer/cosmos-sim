c
c
c     Start a simulation of a galaxy using HYDRO-ART code
c
c
c     Daniel Ceverino (2005)
c
c     ( Based on Start_Cosmology.F by  Andrey Kravtsov (November 2000))
c
c


c     --------------------------
      subroutine Init_Run ( dt )                
c     --------------------------       
c
      real*8 dt
c  
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'
#include "a_def.h"
      character*256 fname, fname_gas
      character*256 fname_part1, fname_part2, fname_part3
c
c
c
ccCEVERINO06182007
ccc
cc      IF (internalB) THEN
ccc
cc#ifdef WORKOUT
cc         call Open_ASCII_File ( iOWork , workfile , append )
cc           write(iOWork,*) 'Before reading from internal buffers'
cc         close(iOWork)
cc#endif
cc         call Read_InternalB ( dt )
cc         call SetUnits ()
cc#ifdef COOLING
cc         call Set_Cooling ()        
cc#endif
cc         call Get_MaxLevelNow ()
cc         do L = MinLevel , MaxLevelNow 
cc            call Magic ( L )
cc            call TreeEoS ( L ) 
cc         enddo
ccc
cc         call Init_Tree ()
ccc
cc         call LL_Construct ()  
ccc
cc#ifdef AGNFEEDBACK
cc        call Read_AGNs () 
cc#endif
ccc
cc          ELSE

#ifdef HYDRO
        ljob  = index(jobname1, ' ') - 1
        lpath = index(path    , ' ') - 1
        fname =
     &        path(1:lpath)//'/'//jobname1(1:ljob)//'.d '

        call Read_Gas_Binary ( dt , fname )
        call SetUnits ()
#ifdef COOLING
        call Set_Cooling ()        
#endif
        call Get_MaxLevelNow ()
        do L = MinLevel , MaxLevelNow         
          call Magic ( L )
          call TreeEoS ( L ) 
        enddo   
#endif
c
	call Init_Tree ()
c
#ifdef NBODY
        ljob  = index(jobname1, ' ') - 1
        lpath = index(path    , ' ') - 1
c
        fname_part1 =
     &                 path(1:lpath)//'/PMcrd.DAT '
        fname_part2 =
     &                 path(1:lpath)//'/PMcrs0.DAT '
        fname_part3 =
     &                 path(1:lpath)//'/pt.dat '

        call Read_Particles_Binary ( fname_part1 , fname_part2 )
        call Read_Particle_Time ( fname_part3 )  ! reads in dt for each particle
#ifdef STARFORM
        fname_part3 =
     &               path(1:lpath)//'/stars.dat '
      IF ( start ) THEN 
         write (*,*)  'No stars at the beginning'
      ELSE
        if ( irun .eq. 2) then
          call Read_Stellar_Data  ( fname_part3 )
        endif
      ENDIF
#endif
        call LL_Construct ()      
#endif


#ifdef AGNFEEDBACK
        call Read_AGNs () 
#endif
c
CCCCCCCCCCC        ENDIF
c
c...  Setting refinement thresholds
c
      const =  (1.0 - Omb0/Om0)
      constb = Omb0/Om0

c
#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'setting refinement thresholds...'
         close(iOWork)
#endif
c

         trho(MinLevel) =wpar(2)*2.0
         trho(MinLevel+1) =wpar(2)*2.0
         trho(MinLevel+2) =wpar(2)*2.0
         trho(MinLevel+3) =wpar(2)*2.0
         tbmass(MinLevel) =   1.*constb
         tbmass(MinLevel+1) = 1.*constb
         tbmass(MinLevel+2) = 1.*constb
         tbmass(MinLevel+3) = 1.*constb
         do Level = 4 , MaxLevel 
            trho(Level)  =  wpar(2)*2.0
            tbmass(Level) = constb/9./8.     ! constb/6.
         enddo
         do Level=MinLevel, MaxLevel
            write(*,*) 'L =',Level,' trho(L) =',trho(Level),
     &           ' trhog(L)=',tbmass(Level)
         enddo

c
c
c...  Refine meshes until the refinement converges. (Daisuke 4/3/03)
c
CEVERINO02152011: THIS IS NOT NEEDED ANYMORE
C      write(*,*)
C      write(*,*) ' Refining meshes immediately down to 2nd level'
C      if ( start ) then
C         call Get_MaxLevelNow ()
C         call Assign_Density ( MinLevel, MaxLevelNow )
C         call Modify ( MinLevel, 1, mtot )
C         write(*,*) 'Level = ',MinLevel,' mtot = ',mtot
C         Level = 1
C         do while ( mtot .ne. nil ) 
C            call Assign_Density ( Level, Level )
C            call Modify ( Level, 1, mtot )
C            write(*,*) 'Level = ',Level,' mtot = ',mtot
C            Level = Level + 1
C         enddo
C      endif
c
c...  Preparing arrays for saving moments
c
      write(*,*) 'preparing arrays'
      do i = 1 , n_save 
        if ( asave(i) .le. aexpn ) then 
          isave(i) = 1
        else
          isave(i) = 0 
        endif      
      enddo

      return
      end
