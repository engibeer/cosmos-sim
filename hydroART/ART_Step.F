C     =====================================================================
c                                                                         .
c          ART Version 3-1 : main hydrodynamic + N-body step              .
c                                                                         .
c            by Alexei Khokhlov and Andrey Kravtsov (1996-1999)           .
c                                                                         .   
c     =====================================================================
c
c     ----------------------------------
      subroutine Step ( dt , iExitCode )
c     ----------------------------------
c
c     iExitCode =  0 - OK
c     iExitCode = -1 - time step failure -> restart with new dt 
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
#ifdef CLOUDY_COOLING
      include 'a_clcool.h'
#endif

      real*8 dt, vas, rs, dtnow
      integer Max_Level
      integer iLStep(MinLevel:MaxLevel) , incStep(MinLevel:MaxLevel)
      integer iLCStep(MinLevel:MaxLevel)
#ifdef STARFORM 
      integer iSF(MinLevel:MaxLevel)
      integer nearest_int
      real logx
CEVERINO04072008
      real*8 DTDA, fp(2), INTEGRATE, age
      external DTDA
      real*8 tSFA
      integer iflagSF
c====== STG:
      real*8 dt_SF
c==========
#endif
      real*8 a2b, b2a
      external a2b, b2a
      real*8 v_tot, v_x, v_y, v_z, ei, e_kin, eratio
      real*8 epot_
      common / COURCELL / ivas, lvas
      integer Nbb(6)
#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append ) 
         write(iOWork,*) 'entering Step dt=',dt
         close(iOWork)
#endif

      iExitCode = nil 

c
c.... organize stepping order, assign time and timestep
c
      nmaxsteps = 2**iTimeBin(MaxLevel)
      nhalfsteps = nmaxsteps/2

#ifdef COSMOLOGY
      at = b2a ( dble(t) ) 
#endif
#ifndef COSMOLOGY
      at = 1.0
#endif
c
      gamma1 = gamma - 1.0 
c
      do Level = MinLevel , MaxLevel
        incStep(Level) = 2**iTimeBin(MaxLevel) / 2**iTimeBin(Level) 
        iLStep (Level) = incStep(Level)         
        iLCStep(Level) = 0 
        aexpold(Level) = aexp(Level)
        aexp   (Level) = at
        tl     (Level) = t 
        dtlold (Level) = dtl(Level)
        dtl    (Level) = dt / 2**iTimeBin(Level)
        if ( istep .eq. 0 ) then 
          iSO    (Level) = 0
          tlold  (Level) = tl(Level)
          dtlold (Level) = dtl(Level)
          aexpold(Level) = at
        endif
#ifdef STARFORM
        dt_SF = dtmin_SF / t0 / aexp(Level)**2
        dtratio = max(dt_SF / dtl(Level) , 1.d-30)
        id = max ( 0 , nearest_int(logx(dtratio,2.)) )
        id = min ( id , iTimeBin(Level) )
        iSF(Level) = 2**id     ! dt_SF(Level) = dtl(Level) * iSF(Level)
#endif
      enddo
      dtmin  = dtl(MaxLevel)
      dtime0 = dtl(MinLevel)

CEVERINO04072008
#ifdef STARFORM
c$$$c      iratioSF= int( (t+dt) / dt_SF) 
      fp(1)  = Om0
      fp(2)  = Oml0
      age =INTEGRATE(DTDA, fp, 2, 0.d0, aexp(0), 0.1d0* aexp(0), 1.d-8)/
     &     (100.d0 * hubble *gyr/mpc)*1.d9
       iratioSF= int( age / dtmin_SF) + 1
c$$$c      ratioSF=   age - 1.d0*iratioSF*dtmin_SF
       tSFA=  1.d0*iratioSF*dtmin_SF - age
      if (tSFA.gt.0.0.AND.tSFA.lt.dt*aexp(0)**2*t0) then 
         iflagSF=1
      else
         iflagSF=0
      endif

ccccc       iflagSF=1
c
c$$$      write (*,*)  'iflagSF= ', iflagSF, t+dt , age, 
c$$$     &     iratioSF ,
c$$$     &     iratioSF*dtmin_SF       , tSFA
#endif
c-------------------------------------------------------------------- 

#ifdef CLOUDY_COOLING

      rs = max ( 1.d0 / at - 1.0 , 0.d0 ) 
CEVERINO04302008      rs =  max ( rs, 8.d0)
      if (rs.lt.8.00) then      ! limit UVBackG to z=8-Flux for high densities
         ioptUV=1
      else
         ioptUV=0
      endif
c-------------------------------------------------------------------- 

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'setting Cloudy cooling rate table for z=',rs
         close(iOWork)
#endif

      call Set_Cooling_Rate_rs ( rs ) ! prepare tables for cooling 

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'done.'
         close(iOWork)
#endif

#endif

#ifndef DEBUG
      call Timing(10,-1)
      call Sync_Tree ()
      call Timing(10, 1)
#endif

#ifdef GRAVITY
      if ( istep .eq. 0 ) then 
        call Get_MaxLevelNow ()

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append ) 
         write(iOWork,*) 'before Assign_Density...'
         close(iOWork)
#endif

        call Timing(1,-1)
        call Assign_Density ( MinLevel , MaxLevelNow )
        call Timing(1,1)

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append ) 
         write(iOWork,*) 'before Poisson...'
         close(iOWork)
#endif

        call Solve_Poisson  ( MinLevel , MaxLevelNow , 0 )         
        call Timing(4,-1)
        call Copy_Potential  ( MinLevel , MaxLevelNow ) 
        call Timing(4,1)
      endif

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'after Poisson...'
      close(iOWork)
#endif
#endif



#ifdef AGNFEEDBACK
      call AGN_SETUP ( )  ! Run HF and set up AGN feedback parameters
      if ( lviewout ) then 
c         call DataOutput3 ( tl(0) , dt )
      endif
#endif

c
      MinModify =  1000
      MaxModify = -1000
      MaxLevelPrev = MaxLevelNow
 
      lstep = istep 
 
      do nstep = 1 , nmaxsteps

        do Level = MaxLevel , MinLevel , -1
          if (iLStep(Level) .eq. nstep ) then
           if ( Level .le. MaxLevelNow ) then 

#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,90)
     &         nstep,np,Level,tl(Level),dtl(Level),aexp(Level)
            close(iOWork)
 90         format(i5,2x,i8,1x,i2,2x,3(g14.6))
#endif

#ifdef HYDRO
#ifdef GRAVITY 
            call Timing(4,-1)
            call Interpolate_Potential ( Level , Level ) 
            call Timing(4,1)
#endif
c
c....       test if time step is still ok
c
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before TimeStepH'
            close(iOWork)
#endif


#ifdef AGNFEEDBACK
c           ViewDump at every some interval of small steps 
c           defined in the DataOutput3 (Daisuke, 7/24/03)
            if ( lviewout .and. Level .eq. MaxLevel .and. 
     &           nstep .eq. 1 ) then 
c              call DataOutput3 ( tl(0) , dt )
c              stop
            endif

c...        AGN feedback (Daisuke, 09/10/03)         
            call AGN_FEEDBACK ( Level )   ! inject energy 

c           ViewDump at every some interval of small steps 
c           defined in the DataOutput3  (Daisuke, 7/24/03)
            if ( lviewout .and. Level .eq. MaxLevel ) then 
c              call DataOutput3 ( tl(0) , dt )
            endif
#endif

            vas = zero
            call Timing(5,-1)
            call TimeStepH ( Level , vas ) 
            call Timing(5,1)

            dtnow = cfl * CellSize(Level) / vas

            if ( dtnow .lt. dtl(Level) ) then 
             dt = 0.5 * cfl * CellSize(MinLevel) / vas 
CEVERINO04222008            dt = 0.33 * cfl * CellSize(MinLevel) / vas 
#ifdef WORKOUT
              call Open_ASCII_File ( iOWork , workfile , append ) 
              write(iOWork,*) 'dt trouble: L =', Level, 'dtl =', 
     &          dtl(Level),' (',dtl(Level)*t0*aexp(Level)**2,' yrs)'
              write(iOWork,*) 'but need ',dtnow, 
     &         ' (',dtnow*t0*aexp(Level)**2,' yrs)'
              write(iOWork,*) 'restarting with new dt =',dt,
     &         ' (',dt*t0*aexp(0)**2,' yrs)'
              call Ps ( ivas , xc , yc , zc ) 
              Tcell = (hvar(7,ivas) - 1.0) * hvar(8,ivas)/hvar(1,ivas) *
     &             T_0 / aexp(lvas)**2
              v_x = hvar(3,ivas)/hvar(1,ivas) * v0 / aexp(lvas)
              v_y = hvar(4,ivas)/hvar(1,ivas) * v0 / aexp(lvas)
              v_z = hvar(5,ivas)/hvar(1,ivas) * v0 / aexp(lvas)
              write(iOWork,*)
     &          'Courant cell: ic L xc, yc, zc rho T v_x,y,z:'
              write(iOWork,91) 
     &                     ivas, lvas, xc, yc, zc, 
     &                     hvar(1,ivas),Tcell,
     &                     v_x, v_y, v_z
              close(iOWork)
 91           format(i8,1x,i2,1x,3(f10.6,1x),2(e11.6,1x),3(e11.3,1x))
#endif
              iExitCode = -1 
              return                 
            endif

#ifdef RadPressure
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before Radiation Pressure'
            close(iOWork)
#endif
            call  RadPre_Compu ( Level )
#endif

#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before StepH'
            close(iOWork)
#endif
            iSt = iSO(Level) 
            call StepH ( Level , dtl(Level) , iSt )
c======== STG: immediately restarts step if Magic gives exit code
            if (iExitCode2 .eq. -1) then
	       return
            endif
c============
            iSO(Level) = iSt
#endif


#ifdef NBODY
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before Move_Level'
            close(iOWork)
#endif
            
#ifdef STARFORM

CEVERINO04072008
            IF (iflagSF.EQ.1) THEN ! ALLOW SF
               if ( mod(iLCStep(Level),iSF(Level)) .eq. 0 ) then 
#ifdef WORKOUT
                  call Open_ASCII_File ( iOWork , workfile , append ) 
                  write(iOWork,*) 'before Star_Formation'
                  close(iOWork)
#endif
                  call Timing ( 9, -1 )
                  call Star_Formation ( Level , iSF(Level),  dt_SF ) ! CEVERINO04072008
                  call Timing ( 9, 1 )
               endif
            ENDIF
#endif

            call Timing(8,-1)
            call Move_Level ( Level )
            call Timing(8, 1)
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before LL_Update'
            close(iOWork)
#endif
            call Timing(8,-1)
            call LL_Update ( Level , MinModify , MaxModify )
            call Timing(8, 1)
#endif



#ifdef GRAVITY
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before Assign_Density'
            close(iOWork)
#endif
            call Timing(1,-1)
            call Assign_Density ( Level , Level )
            call Timing(1,1)
#endif        
            tlold(Level) = tl(Level)
            tl   (Level) = tl(Level) + dtl(Level)

#ifdef COSMOLOGY
            aexpold(Level) = aexp(Level)
            aexp(Level) = b2a ( dble(tl(Level)) )
#endif

#ifdef GRAVITY
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before Copy_Potential & Solve_Poisson'
            close(iOWork)
#endif
            call Timing(4,-1)
            call Copy_Potential ( Level , Level )
            call Timing(4, 1)
            call Solve_Poisson ( Level , Level , iLCStep(Level) ) 
#endif

#ifdef REFINE 
#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'before Modify'
            close(iOWork)
#endif
            call Timing(7,-1)
            call Modify ( Level, 0 , mtot )
            call Timing(7,1)
            call Get_MaxLevelNow ()
#endif

CEVERINO01292008
            Level1= min(Level+1,MaxLevelNow)
            Levelm1=max(0,Level-1)
CEVERINO02222012
            DO ilevel=Levelm1,Level1
               call Magic ( iLevel )
               call TreeEoS ( iLevel )
            ENDDO
c======== STG: immediately restarts step if Magic gives exit code
            if (iExitCode2 .eq. -1) then
	       return
            endif
c============


#ifdef WORKOUT
            call Open_ASCII_File ( iOWork , workfile , append ) 
            write(iOWork,*) 'after Modify. the end of nstep=', nstep
            close(iOWork)
#endif

CEVERINO01232008
C#ifdef WORKOUT
C      call Open_ASCII_File ( iOWork , workfile , append ) 
C      write(iOWork,*) 'Before Magic & EoS'
C      close(iOWork)
C#endif
C            call Magic ( Level )
C            call TreeEoS ( Level )
c-------------------------------------------------------------------- 

            iLStep(Level) = iLStep(Level) + incStep(Level)
            iLCStep(Level) = iLCStep(Level) + 1


           else ! level does not exist yet
             tlold(Level)  = tl(Level)
             dtlold(Level) = dtl(Level) 
             tl   (Level) = tl(Level) + dtl(Level)
#ifdef COSMOLOGY
             aexpold(Level) = aexp(Level)
             aexp(Level) = b2a ( dble(tl(Level)) )
#endif
             iLStep(Level) = iLStep(Level) + incStep(Level)
             iLCStep(Level) = iLCStep(Level) + 1            
           endif

          endif
        enddo ! end Level
 
c       ViewDump at every some interval of small steps 
c       defined in the DataOutput3  (Daisuke, 7/24/03)
c        if ( lviewout ) call DataOutput3 ( tl(0) , dt )
        
      enddo ! end nstep

      call Get_MaxLevelNow ()
c     
#ifdef REFINE 
#ifdef WORKOUT
        call Open_ASCII_File ( iOWork , workfile , append ) 
        write(iOWork,*) 'before final Density. MaxLevelNow =', 
     &                  MaxLevelNow
        close(iOWork)
#endif
        call Timing(1,-1)
        call Assign_Density ( MinLevel , MaxLevelNow )
        call Timing(1,1)
#ifdef WORKOUT
        call Open_ASCII_File ( iOWork , workfile , append ) 
        write(iOWork,*) 'before final Modify'
        close(iOWork)
#endif
        do Level = MinLevel , MaxLevelNow 
#ifdef WORKOUT
          call Open_ASCII_File ( iOWork , workfile , append ) 
          write(iOWork,*) 'Modifying level', Level
          close(iOWork)
#endif
          call Timing(7,-1)
          call Modify ( Level , 0 , mtot )
          call Timing(7,1)
        enddo
#endif
      call Get_MaxLevelNow ()

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'before AdvanceTime'
      close(iOWork)
#endif

c
c.... advance time variables, compute energies
c
      call AdvanceTime ( dt )

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'after Advance'
      close(iOWork)
#endif
c

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'finished Step, exiting with iExitCode =', 
     &                 iExitCode
      close(iOWork)
#endif
#ifdef timingStep
      write (*,*)  'timingStep'
      write (*,*)'Sweep , before Sw,  after Sw , 2nd part of Hstep , HS'
      write (*,60) CPU(12)/60. , CPU(13)/60. , CPU(14)/60. , CPU(15)/60. 
     & , CPU(6)/60.
 60   format( 5(1x,f7.3) )
      call Write_Timing () 
      stop
#endif
#ifdef timingDensity                                                               
      write (*,*)  'timingDensity'                                                            
      write (*,*)  'Density,  paral,  serial, LevelsLoop'                                     
      write (*,60) CPU(1)/60. , CPU(13)/60. , CPU(14)/60.                                     
     &     ,CPU(15)/60. , CPU(16)/60. , CPU(17)/60. , CPU(18)/60.                             
 60   format( 7(1x,f7.3) )                                                                    
      call Write_Timing ()                                                                    
      stop                                                                                    
#endif   
      return
      end          
c
c     ----------------------------------------
      subroutine StepH ( Level, dtLevel, iSt )
c     ----------------------------------------
c
c     hydrodynamical step on Level with timestep dtLevel
c
c     routine TreeEos specifies the equation of state for gas
c     and is to be supplied by user
c
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 dtLevel
      external TreeEos, Magic
c
      integer iorder(3,0:1)   ! sweep order (3D case only)
      data iorder / 2, 4, 6,  ! XYZ
     &              6, 4, 2 / ! ZYX

      call Timing ( 6 , -1 ) 
#ifdef timingStep
        call Timing (13, -1)
#endif
c
c.... Set vnw variables in leaves of Level-1 which will receive fluxes
c     and will accumulate vnw during integration of leaves of Level.
c
      if ( (iSt .eq. 0) .and. (Level .gt. MinLevel) ) then 
        call H_Old_to_New ( Level-1 , 1 )
      endif
c
      call H_Old2New_3b ( Level )
#ifdef timingStep
        call Timing (13, 1)
#endif
c
c.... save density and internal energy to estimate Gamma after the step
c
c      call Estimate_Effective_Gamma ( Level , 0 ) 
c
c.... Perform XYZ sweeps.
c
      do l = 1 , ndim   
c
c...... Set vnw variables in leaves of Level which are not
c       accumulating vnw during integration of leaves of Level+1
c       (e.g., leaves of Level with no split neighbors)
c
#ifdef timingStep
        call Timing (13, -1)
#endif
        call H_Old_to_New ( Level, 2 )
#ifdef timingStep
        call Timing (13, 1)
#endif
c
c...... Make a sweep in direction iDir
c

        iDir = iorder(l,iSt) 
#ifdef timingStep
        call Timing (12, -1)
#endif
        call Sweep ( iDir, dtLevel, Level )        
#ifdef timingStep
        call Timing (12, 1)
#endif

#ifdef timingStep
        call Timing (14, -1)
#endif
c
c.....  vnw --> hvar for all leaves with no split neighbors
c
        if ( l .ne. ndim ) call H_New_to_Old ( Level, 1 )
c
c.....  vnw --> hvar for all leaves
c
        if ( l .eq. ndim ) call H_New_to_Old ( Level, 0 )

c...... Make Magic, Compute EoS 

        if ( l .ne. ndim ) then 
          call Magic   ( Level )

          call TreeEoS ( Level )

        endif
#ifdef timingStep
        call Timing (14, 1)
#endif
      end do

      call Timing ( 6 , 1 ) 

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'Before Advance_IE'
      close(iOWork)
#endif

      call Timing ( 11 , -1 ) 
      call Advance_InternalEnergy ( Level )
      call Timing ( 11 , 1 ) 

#ifdef WORKOUT 
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'After Advance_IE'
      close(iOWork)
#endif

c
      call Timing ( 6 , -1 ) 
#ifdef timingStep
        call Timing (15, -1)
#endif

      if ( Level .lt. MaxLevelNow ) call SplitUpdate ( Level )

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'Before Magic '
      close(iOWork)
#endif

      call Magic ( Level )

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'After Magic'
      close(iOWork)
#endif

c
c      call Estimate_Effective_Gamma ( Level , 1 )
c      call Smooth_Gamma_eff ( Level ) 

      call TreeEoS ( Level )

      if ( Level .lt. MaxLevelNow ) call SplitUpdate ( Level )


#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'end of sweep'
      close(iOWork)
#endif
c
c.... change sweep direction for this level 
c

      iSt = abs(iSt - 1)
c
#ifdef timingStep
        call Timing (15, 1)
#endif
      call Timing ( 6 , 1 )

      return
      end

