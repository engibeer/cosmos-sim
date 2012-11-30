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
      use density_alloc
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
c...  Initialize a seed for the random number generator
c
      mrand = 349884

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
CEVERINO04242011  
      allocate(vvar(NCPUS,nctot))
      allocate(vref(NCPUS,nctot))
c
      call Assign_Density ( MinLevel , MaxLevelNow )
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
c=======STG:
      iExitCode2 = 0
c==========
c

      DO WHILE ( Continue_Run )

c========STG: disabled this if-block to use Pedro's refinement scheme
CEVERINO10192007
c         if (aexpn.lt.aHR) then
c            MinL_Jeans= 9  
c            lmax_ref(1)=9
c            lmin_ref(2)=13 ! no gas ref
c            lmax_ref(2)=13
c         elseif (aexpn.ge.aHR.AND.aexpn.lt.0.160) then
c           if (MinL_Jeans.ne.10) write (*,*)  'Allow another level= ',10 ! check
c            MinL_Jeans=10
c            lmax_ref(1)=10
c            lmin_ref(2)=4 ! Gas ref 
c            lmax_ref(2)=10
c         elseif (aexpn.ge.0.160.AND.aexpn.lt.0.320) then
c           if (MinL_Jeans.ne.11) write (*,*)'Allow another level= ', 11 ! check
c            MinL_Jeans= 11
c            lmax_ref(1)=11
c            lmin_ref(2)=4 ! Gas ref 
c            lmax_ref(2)=11
c         elseif (aexpn.ge.0.320) then
c           if (MinL_Jeans.ne.12) write (*,*)'Allow another level= ', 12 ! check
c            MinL_Jeans= 12
c            lmax_ref(1)=12
c            lmin_ref(2)=4 ! Gas ref 
c            lmax_ref(2)=12
c         endif
c======================
c
c
c.... start timing the step
c
 10      call Timing ( ntiming-1 , -1 )

c======STG:
         iExitCode2 = 0
c=========
        
         call Step ( dt , iStepExitCode )
c======STG:
#ifdef WORKOUT
          call Open_ASCII_File ( iOWork , workfile , append )
          write(iOWork,*) 'iexitcode2 =', iExitCode2
          close(iOWork)
#endif    
c=========

c        if ( iStepExitCode .eq. 0 ) then 
c========STG:
         if ( iStepExitCode .eq. 0 .and. iExitCode2 .eq. 0 ) then 
c===========
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
c          if ( iStepExitCode .eq. -1 ) then
c=========STG:
           if ( iStepExitCode .eq. -1 .or. iExitCode2 .eq. -1) then
c===========
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
            call Timing ( ntiming-1 , 1 )
c 
 	    if (iExitCode2 .eq. -1) then
              dt = 0.3*dt
	    endif

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
        call Write_Timing () 
      ENDDO                            ! main loop ends
c
      call Save ( dt , 0 ) 
c



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

      real*8 b2a, dt

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

cc     --------------------------
c      subroutine Test_Density ()
cc     -------------------------- 
c      include 'a_tree.h'
c      include 'a_control.h'
c#     include "a_def.h"
c      include 'a_dnew.fh' 

c      open ( 39 , file = 'test_dens.dat' )
c      DO Level = MinLevel , MaxLevel 
c        call Select_Cells ( Level , nLevel ) 
c        do ic1 = 1 , nLevel
c          icell = iSelect(ic1)
c          sum = 0.d0
c          do ic2 = 0 , 7
c            idc = icell + ic2
c            write(39,*) Level, idc, var(1,idc), ref(idc)
c          enddo
c        enddo
c      ENDDO
c      close ( 39 )
c      return
c      end
