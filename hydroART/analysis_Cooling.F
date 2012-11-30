c--------------------------------------------------------------------
c--------------------------------------------------------------------
c     Code to test cooling tables from CLOUDY
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
c     by Daniel Ceverino (March 2006)
c
c     It should be compiled along with analysis.ART.F
c
c     It reads gas conditions from Gas_Cooling.dat:
c
c     T_Cooling.exe < Gas_Cooling.dat
c
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'
      include 'a_clcool.fh'
#     include "a_def.h"

      real*8  rs,Zdum
      character*256 filename, headerCo

      boxh = 6.0
      Om0 = 0.300
      Oml0= 0.700
      Omb0= 0.04490 
      hubble = 0.700
      gamma = 1.6666

      call Units ()
      call Set_Cooling ()
! choose tasks:
      read(*,*) iflag
      if (iflag.eq.1) write (*,*)  "Cooling Test:"
      if (iflag.eq.0) write (*,*)  "Compute equilibrium temperature"
!     Set gas conditions:
      read(*,*) headerCo
      read(*,*) rho             ! nH density in cm-3
      read(*,*) Zdum            ! Metallicity 
      read(*,*) rs              ! Redshift
      read(*,*) filename        ! OUtput file
c     
      aexpn = 1./(1+rs)
c     
      write(*,*) "Gas conditions for a Cooling Test:"
      write(*, '( "Hydrogen number density[cm-3]= ", F5.2)') rho             
      write(*, '( "Metallicity = ", F5.2)') Zdum  
      write(*,'( "Redshift= ", F5.2)') rs            
      write(*,'( " Output file= ", A30)') filename  
      write (*,*) "a= ",   aexpn
      
!     #ifdef CLOUDY_COOLING        
!     rs = max ( 1.d0 / aexpn - 1.0 , 0.d0 ) 
      call Set_Cooling_Rate_rs ( rs ) 
      write(*,*) 'setting Cloudy cooling rate table for z=',rs
!     #endif
      
      if (iflag.eq.1) then
         call Test_Cooling ( filename, rho, Zdum,rs )
      else if (iflag.eq.0) then
         call EqT (filename, rho, Zdum, rs)
      endif

      end
c--------------------------------------------------------------------
c     Subroutines taken from ART_Cooling.F
c
c     ------------------------------------
      subroutine Test_Cooling ( filename, rho, Zdum,rs )
c     ------------------------------------
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      real*8 Tdum, cr1, cr2, cr3, Zdum, rhogl, xion, rs
      character*256 filename 
c
      rhogl =log10(rho)
      rho2 = rho**2
      coolrate_transferCU = (v0**2 * rho0)/t0
      rhoCU = rho*aexpn**3 / (1.125e-5 * hubble**2*Om0)
      dummy1 = rhoCU**2 * aexpn
      eiFACTOR = rhoCU*aexpn**2/( T_0*(gamma-1) )
c
      lfn = index ( filename , ' ' ) - 1 
      open ( 20 , file = filename(1:lfn) ) 
      Tlmin = 1.d0 
      Tlmax = 9.d0
      dTlog=0.01 
      nt = int((Tlmax - Tlmin)/dTlog) + 1
      do it = 1 , nt  
         Tdum = 10.d0**(Tlmin + float(it-1)*dTlog)
         call Cooling_Rate ( rhogl, Tdum*1.d-4, Zdum, cr1 ) ! cr1 in erg/s/cm3
!         cr1 = ( (rho*aexpn**3/rho0)**2 * aexpn ) * cr1 ! in code units
c
c         cr1 = cr1*dummy1 ! Gamma net in code units
c
         ei = Tdum*eiFACTOR
         dtcool =abs( ei / ( cr1*dummy1) )
c
c!!!         cr1 = coolrate_transferCU * cr1*aexpn**-7 ! in ergs s-1 cm-3
c
          cr1 = cr1 / ( 1.d23*AL_SD) !      Gamma P net in erg s-1 cm3
c
          call  X_ion ( rhogl , Tdum , Zdum , rs , xion )
c
          write(20,90) Tdum, cr1, dtcool*aexpn**2*t0, xion
      enddo
 90   format (4(e14.6,1x))
      close ( 20 ) 

      return
      end   
c-------------------------------------------------------------------- 
c     Compute equilibrium temperature: cooling-heating balance
c
      SUBROUTINE EqT (filename, rho, Zdum, rs)
c
c
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      real*8 Teq, cr1,  Zdum, rhogl, xion, rs
      real*8 rhoglmin
      character*256 filename 

      parameter (rhoglMax=2.)
      parameter (drhogl=0.01)
      rhoglmin =log10(rho)
c
      lfn = index ( filename , ' ' ) - 1 
      open ( 20 , file = filename(1:lfn) ) 
c
      Tlmin = 1.d0 
      Tlmax = 9.d0
      dTlog=0.01 
      nt = int((Tlmax - Tlmin)/dTlog) + 1
      nrho=int ( (rhoglMax - rhoglmin)/drhogl) + 1
      do irho=1, nrho
         rhogl=rhoglmin + float(irho-1)*drhogl
         it=1
         cr1=100.

      DO WHILE (cr1.gt.0.)
         Teq = 10.d0**(Tlmax - float(it-1)*dTlog)
         call Cooling_Rate ( rhogl, Teq*1.d-4, Zdum, cr1 ) ! cr1 in erg/s/cm3
         it=it+1
         if ( it.ge.nt) then
            cr1=-100.
            Teq=0.
         endif
      ENDDO

      write (20,*)  rhogl , Teq 
      enddo
      end
