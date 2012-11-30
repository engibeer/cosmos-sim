c     ====================================================
c
c                Adaptive Refinement Tree (ART) 
c
c                N-body + hydrodynamics  solver
c                 for cosmological simulations
c
c         N-body : Version 3 - February  1997 - June 1998
c         Hydro  : Version 1 - September 1998 - 
c
c       Andrey Kravtsov, Anatoly Klypin, Alexei Khokhlov
c
c                         (c) 1996-2001
c
c     ====================================================
c
      program ART
c
c
#     include "a_def.h"
c
      include 'a_tree.h'
      include 'a_control.h'
      logical Continue_Run
c
      integer nMove
      integer mtot
      common / COURCELL / ivas, lvas
      real*8 b2a
      real*8 dt, dtold

c
c...  Initialize a seed for the random number generator
c
      mrand = 349884

CEVERINO12052005
c
c	Initialize flag for internal buffers to save previous step
c
	internalB = .FALSE.

c
c.... read/set in initial conditions, initialize variables and linked lists
c
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'starting, reading control parameters ...'
      close(iOWork)
#endif
c
      do i = 1 , ntiming
        CPU(i) = 0.
      enddo
c
      call ReadControl ()
c
      call Initialize ( dt )
c
      call Assign_Density ( MinLevel , MaxLevelNow, 0 )

c      Preaheating of Epreheat [keV/part]
c      call PREHEAT ( 0.3 )  ! 0.3keV per particle 
CEVERINO4/18/2004 it is in  ART_Heating.F
c
      if ( t .gt. 0. ) then 
        Continue_Run = .false.
        stop
      else
        Continue_Run = .true.         
      endif
c
      iStepExitCode = 0
      iStepp = 0
c

      DO WHILE ( Continue_Run )
c
c
c.... start timing the step
c


10      call Timing ( ntiming-1 , -1 )

        
        call Step ( dt , iStepExitCode )

        if ( iStepExitCode .eq. 0 ) then 
#ifdef WORKOUT
          call Open_ASCII_File ( iOWork , workfile , append )
          write(iOWork,*) 'choosing Time Step...'
          close(iOWork)
#endif
c
          dtold = dt
          call Timing ( 4 , -1 )
          call ChooseTimeStep ( dt )
          call Timing ( 4 ,  1 )
c
#ifdef WORKOUT
          call Open_ASCII_File ( iOWork , workfile , append )
          write(iOWork,*) 'old dt =',dtold, ' new dt =', dt
          close(iOWork)
#endif
        else
          if ( iStepExitCode .eq. -1 ) then
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append )
            write(iOWork,*) 'restarting with new dt...' 
            close(iOWork)
#endif
c           
            if ( istep .gt. 1 ) then 
              start = .false.

            endif
c
            call Timing( ntiming-1 , 1)

c
            call Initialize ( dtold )
c
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append )
            write(iOWork,*)'old dt =',dtold,
     &         ' (',dtold*t0*aexp(0)**2,' yrs)'
            write(iOWork,*)'new dt =',dt,
     &         ' (',dt*t0*aexp(0)**2,' yrs)'
            close(iOWork)
#endif
c
            go to 10 
c
          else
            call Open_ASCII_File ( iOErr , errorfile , sequent )
            write(iOErr,*) '1  *bad exit code from Step:',iStepExitCode
            close(iOErr)
            stop
          endif
        endif

c        call DataOutput ( tl(0) , dt )

        call Save_Check ( dt )

        if ( start ) start = .false. 
#ifdef COSMOLOGY
        a_exp = b2a ( dble(t) ) 
        if ( t .gt. 0.0 ) Continue_Run = .false.

        iStepp =iStepp +1
        if ( iStepp.ge.nsteprun ) Continue_Run = .false.
#endif
c
        call Write_Timing ()
c
      ENDDO                            ! main loop ends
c
      call Save ( dt , 0 ) 
c
      STOP
      END

c     -----------------------------
      subroutine Initialize ( dt )
c     -----------------------------
c
c
#     include "a_def.h"
c
      include 'a_tree.h'
      include 'a_control.h'

      real*8 b2a,dt

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'initializing arrays...'
      close(iOWork)
#endif
c
      call Init_Parameters       ()
      call Init_Arrays           ()
      call Init_Tree             ()
c
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'before Init_Run...'
      close(iOWork)
#endif

      call Init_Run ( dt )

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'after Init_Run...'
      close(iOWork)
#endif

      call Get_MaxLevelNow ()

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'MaxLevelNow =', MaxLevelNow
      close(iOWork)
#endif

#ifdef COSMOLOGY
      at = b2a ( dble(t) )
#endif
#ifndef COSMOLOGY
      at = 1.0
#endif
c
      aexpold(0) = at
      aexp(0)    = at
c
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'at =', at  
      close(iOWork)
#endif

c
      do i = 1 , ntiming 
        CPU(i) = 0.
      enddo
c
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'exiting Initialize...'
      close(iOWork)
#endif

      return
      end      

c     --------------------------
      subroutine Test_Density ()
c     -------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      include 'a_dnew.fh' 

      open ( 39 , file = 'test_dens.dat' )
      DO Level = MinLevel , MaxLevel 
        call Select_Cells ( Level , nLevel ) 
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          sum = 0.d0
          do ic2 = 0 , 7
            idc = icell + ic2
            write(39,*) Level, idc, var(1,idc), ref(idc)
          enddo
        enddo
      ENDDO
      close ( 39 )

      return
      end
