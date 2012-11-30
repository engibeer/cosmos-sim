c     =====================================================================
c                                                                         
c               ART Version 3-1 : choose global time step              
c                                                                         
c              by Alexei Khokhlov and Andrey Kravtsov (1998)                
c                                                                         
c     =====================================================================
c
c     --------------------------------
      subroutine ChooseTimeStep ( dt ) 
c     --------------------------------
c
c     on input  dt is the current timestep
c     on output dt is the new timestep
c
c     cfl    - is the Courant number 
c     timinc - parameter controlling timestep increase
c
c     cfl and timinc are to be supplied by the user when run is initialized
c
      real*8 dt , dt1 , tdum, adum1, adum2
      real*8 a2b, b2a
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
C
CEVERINO09232005:
#ifdef DEBUGdcr
      call Open_ASCII_File ( iODebug , debugfile , append ) 
      write (iODebug,*) 'Debuging of ART_TimeStep'
      write (iODebug,*) '-----------------------------------------'
      write (iODebug,*) 'Ad initio, dt= ', dt
      close (iODebug)
#endif
c
      call TimeStep ( dt1 )
c
c      dt1 = min(dt1,0.125)
c
#ifdef COSMOLOGY
      adum1 = b2a ( (1.d0*t) ) 
      tdum  = t + dt1
      adum2 = b2a ( (1.d0*tdum) )           
      dda = min( (adum2 - adum1)/adum1 , 1.d-1 ) ! 
      dda = min( dda*adum1 , 3.d-3 )
      tdum = a2b ( (1.d0*(adum1+dda)) )
      dt1   = tdum - t
#endif
c
CEVERINO09232005:
#ifdef DEBUGdcr 
      call Open_ASCII_File ( iODebug , debugfile , append ) 
      write (iODebug,*) 'After Cosmology:'
      write (iODebug,*) 'dt= ', dt1 , t , tdum 
      write (iODebug,*) 'adum1= ',adum1, ' adum2= ', adum2,
     & ' dda= ', dda 
      close (iODebug)
#endif
      if ( istep .eq. 0 ) then 
        dt = dt1
      else
        if ( dt1 .lt. dt ) then
          dt = dt / 1.25
          if ( dt1 .lt. dt ) dt = dt1
          istep2 = 0
        else
          if ( istep2 .gt. ntc .and. dt1 .gt. timinc*dt ) then
            dt = dt * timinc
            istep2 = 0
          endif
        endif
      endif
CEVERINO09232005:
#ifdef DEBUGdcr 
      call Open_ASCII_File ( iODebug , debugfile , append ) 
      write (iODebug,*) 'Finally, before dt = min ( dt , 0.125d0 )  :'
      write (iODebug,*) 'dt= ', dt, istep
      close (iODebug)
#endif
      dt = min ( dt , 0.125d0 )  
c      
      return
      end
c
c     --------------------------
      subroutine TimeStep ( dt )
c     --------------------------
c
c     Output  :     dt - normalized timestep 
c
c     Courant condition from hydrodynamical variables
c
c     Note : loops are not parallel
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 dt, vas 
      integer Level 
c
      common / COURCELL / ivas, lvas
      common / COOLCELL / icool, lcool, cool_rate
c
      call Get_MaxLevelNow ()
      vas = zero
      ivas = 0 
      lvas = 0 
#ifdef HYDRO
      do Level = MinLevel , MaxLevelNow 
        call TimeStepH ( Level , vas ) 
      enddo 
      call Ps ( ivas , xc , yc , zc ) 
      Tcell = (gamma - 1.0) * hvar(8,ivas)/hvar(1,ivas) *
     &         T_0 / aexpn**2
      v_x = hvar(3,ivas)/hvar(1,ivas) * v0 / aexpn
      v_y = hvar(4,ivas)/hvar(1,ivas) * v0 / aexpn
      v_z = hvar(5,ivas)/hvar(1,ivas) * v0 / aexpn
      write(*,*) 'Courant cell: ',ivas, lvas, hvar(1,ivas),Tcell,
     &       v_x, v_y, v_z, hvar(7,ivas)
c
#endif
c
      if ( vas .le. zero ) then 
        write(*,*) 'TimeStep : error : vas = 0'
        stop
      endif
c
      dt = CellSize(MinLevel) / vas  ! required time step for level 0
      if ( cfl .eq. zero ) then
        write (*,*) 'ChooseTimeStep : error : cfl = 0'
        stop
      endif
      if ( timinc .eq. zero ) then
        write (*,*) 'ChooseTimeStep : error : timinc = 0'
        stop
      endif
c
      dt = cfl * dt
c
#ifdef COOLING 
c      dt_cool = 1.e6 
c      do Level = MinLevel , MaxLevelNow
c        call CoolingTimeStep ( Level , dt_cool )
c        write(*,*) 'L =',Level,' dt_cool =',dt_cool
c      enddo
c      icool = max(icool,1)
c      Tcell = (gamma - 1.0) * hvar(8,icool)/hvar(1,icool) * 
c     &        T_0 / aexpn**2
c      write(*,*) 'Cooling cell: ',icool,lcool,hvar(1,icool),Tcell,
c     &           hvar(7,icool)  
#endif

c
c      dt = min ( dt , dt_cool ) 
c
CEVERINO09232005:
#ifdef DEBUGdcr 
      call Open_ASCII_File ( iODebug , debugfile , append ) 
      write (iODebug,*) 'Courant cell: ',ivas, lvas, hvar(1,ivas),Tcell,
     &       v_x, v_y, v_z, hvar(7,ivas)
      write (iODebug,*) 'vas= ', vas, 'timinc =', timinc, ' cfl= ',cfl,
     & ' dt= ', dt 
      close (iODebug)
#endif
      return
      end
c
c     ------------------------------------
      subroutine TimeStepH ( Level , vas ) 
c     ------------------------------------
c
c     Courant condition from hydrodynamical variables
c
c     Note : loops are not parallel
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level , icell , idcell , ic
      real*8 vas1 , vas , rhor , as , vel
      common / COURCELL / ivas, lvas
#ifdef DEBUGdcr 
      real*8 dt , dt1 , tdum, adum1, adum2
      real*8 a2b, b2a
#endif
c
      vas1 = zero
c                
      IF ( Level .eq. MinLevel ) THEN 
c
ccC$OMP PARALLEL DO DEFAULT(PRIVATE) REDUCTION(MAX:vas1)
        do ic = 1 , ncell0 
          if ( iOctCh(ic) .eq. nil ) then  ! only leaves
            rhor = one / hvar(1,ic) 
c======STG:
#ifdef RadPressure
         as   = sqrt ( (hvar(6,ic)+RadPre(ic)) * gamma * rhor )
#endif
#ifndef RadPressure
         as   = sqrt ( hvar(6,ic) * gamma * rhor )
#endif
c==========
c            as   = sqrt ( hvar(6,ic) * gamma * rhor )
            as1  = sqrt ( (gamma-1.)*hvar(8,ic) * gamma * rhor )
            vel  = max ( abs ( hvar(3,ic) ),
     &                   abs ( hvar(4,ic) ),
     &                   abs ( hvar(5,ic) )  
     &                 ) * rhor + as
            if ( vel .gt. vas .and. vel .ge. vas1 ) then
	      ivas = ic
              lvas = 0  
            endif
              vas1 = max ( vas1 , vel )
          endif
        enddo
c
      ELSE
c
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            if ( iOctCh(idcell) .eq. nil ) then  ! only leaves
              rhor = one / hvar(1,idcell)
c======STG:
#ifdef RadPressure
         as   = sqrt ( (hvar(6,idcell)+RadPre(idcell)) * gamma * rhor )
#endif
#ifndef RadPressure
         as   = sqrt ( hvar(6,idcell) * gamma * rhor )
#endif
c==========
c              as   = sqrt ( hvar(6,idcell) * gamma * rhor )
              as1  = sqrt ( (gamma-1.)*hvar(8,idcell) * gamma * rhor )
              vel  = max ( abs ( hvar(3,idcell) ),
     &                     abs ( hvar(4,idcell) ),
     &                     abs ( hvar(5,idcell) )  
     &                   ) * rhor + as
              if ( vel .gt. vas .and. vel .ge. vas1 ) then 
                ivas = idcell
                lvas = Level 
              endif
              vas1 = max ( vas1 , vel )
            endif
          enddo
        enddo
c
      ENDIF
c
      vas = max ( vas, vas1 , 1.d0*small )     
c      write(*,*) ' Debug TimeStepH : ',ivas,vas,vas1,as,as1,
c     &     hvar(6,ivas),hvar(8,ivas),hvar(1,ivas)
c
CEVERINO09232005:
#ifdef DEBUGdcr 
      call Open_ASCII_File ( iODebug , debugfile , append ) 
      write (iODebug,*)  ' Debug TimeStepH : ',ivas,vas,vas1,as,as1,
     &     hvar(6,ivas),hvar(8,ivas),hvar(1,ivas)
              Tcell = (hvar(7,ivas) - 1.0) * hvar(8,ivas)/hvar(1,ivas) *
     &             T_0 / aexp(lvas)**2
              v_x = hvar(3,ivas)/hvar(1,ivas) * v0 / aexp(lvas)
              v_y = hvar(4,ivas)/hvar(1,ivas) * v0 / aexp(lvas)
              v_z = hvar(5,ivas)/hvar(1,ivas) * v0 / aexp(lvas)
c=======STG:
              call Ps ( ivas , xc , yc , zc ) 
c==============
              write(iODebug,*)
     &          'Courant cell: ic L xc, yc, zc rho T v_x,y,z:'
              write(iODebug,91) 
     &                     ivas, lvas, xc, yc, zc, 
     &                     hvar(1,ivas),Tcell,
     &                     v_x, v_y, v_z
 91           format(i8,1x,i2,1x,3(f10.6,1x),2(e11.6,1x),3(e11.3,1x))
              dtnow = cfl * CellSize(0) / vas 
              write (iODebug,*) 'dt_C= ', dtnow*t0*aexp(lvas)**2
c
      dt1=dtnow
      adum1 = b2a ( (1.d0*t) ) 
      tdum  = t + dt1
      adum2 = b2a ( (1.d0*tdum) )           
      dda = min( (adum2 - adum1)/adum1 , 1.d-1 ) ! 
      dda = min( dda*adum1 , 3.d-3 )
      tdum = a2b ( (1.d0*(adum1+dda)) )
      dt1   = tdum - t
c
      dtnow=dt1
              write (iODebug,*) 'dt_C= ', dtnow*t0*aexp(lvas)**2 
      close (iODebug)
#endif
      return
      end
c
