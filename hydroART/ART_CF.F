c     =====================================================================
c                                                                         
c     ART : source steps for some of the gas equations 
c                                                                         
c     by Andrey Kravtsov (2001)    
c
c
c     stellar heating implemented by Daniel Ceverino (2007)
c
c     =====================================================================
c

c     -------------------------------------------
      subroutine Advance_InternalEnergy ( Level )
c     -------------------------------------------
c
c     source step  for internal energy 
c
      integer Level
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      real*8  esmall, Zsolar, Zsolari
      parameter ( Zsolar = 0.02 , Zsolari = 1.d0 / Zsolar )
      parameter ( eta1 = 0.1, eta2 = 1.0e-3 ) 
      parameter ( dstep = 0.01 , esmall = 0. )
      parameter ( gammin = 1.01 )        
      real*8 b2a, t_e , dt_e, dtcool, dtcoolmin, dt3, ei1
      real*8 tbegin, tend, e_kin, ei, ei2, eratio, a1, a2, afact, acurr
      real*8 T_gas, ai, coolrate , coolcomp , rhog2, rhogi 
      real*8 dummy0, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6
      real*8 Tfact, Tfact1, Tfact2, rhogl, fact_nH
      real*8 e_pdV, e_cool, e_flux, e_small, e_min, Zdum 
CEVERINO04162007      real*8 GravGrad,PresGrad,absGravGrad, absPresGrad
      real*8 substractSN

      logical continue_cool

      integer Nbb(6)
#ifdef FEEDBACKenergy
      real*8 heatRate, heatRateSNII
      real*8 t_fb_code,t_fbi_code,cfbe,efact,TfactSN,heatfact
     &     ,e_ceil,heatfactSN, t_fbSN_code, tfbc
c----- STG changed T_ceil to 2e7:
      parameter ( T_ceil = 2.0d7 ) ! do not allow gas to be feedback-heated above this T [K]
c---------------------------
      parameter ( t_fbSN = 1.0d5 ) ! period of a supernova explosion in yrs
#endif
#ifdef FEEDBACK_SNIa
      real*8 heatRateSNIa, xd, dN_SNIadtdM
      real*8 t_SNIa_code, R_SNIa_fact, cfbeIa
#endif
#ifdef KineticFCONSTANT
      parameter (Sigma_KinFP=30 ) ! km/s
      integer iNbb(6)
#endif
      dt3 = dtl(Level) / 3.0
      tbegin = tl(Level)
      tend   = tl(Level) + dtl(Level)
#ifdef COSMOLOGY
      ai = 1.d0 / aexp(Level) 
      a1 = b2a ( dble(tl(Level)) )
      a2 = b2a ( dble(tl(Level) + dtl(Level)) )
      afact = (a2 - a1) / dtl(Level) ! stuff for linear interpolation ofa(t)
#endif
#ifndef COSMOLOGY
      a1 = 1.0
      a2 = 1.0 
      ai = 1.0
      afact = 0.
#endif

#ifdef SOLARMETALL
      Zdum=0.0
#endif
#ifndef METALCOOLING
      Zdum = -1000.0
#endif
#ifndef CLOUDY_COOLING
      rhogl = 0.0 
#endif
#ifdef CLOUDY_COOLING
c.... factor converting gas density in code units to n_H in cm^-3
      fact_nH = log10(1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3)
#endif
#ifdef STARFORM
      nstar1 = nsp(nspecies,1) - 1 
#endif
#ifdef FEEDBACKenergy
      t_fb_code = t_fb / t0 / aexp(Level)**2
      t_fbi_code= 1.0 / t_fb_code
#ifndef KineticF
      cfbe = C_fb * aexp(Level)**2 * CellVoli(Level) 
#endif
#ifdef KineticF
      cfbe = C_fb * aexp(Level)**2 * CellVoli(Level) *0.9 ! 0.1 FOR KINETIC FEEDBACK 
#endif
      TfactSN = T_0 / aexp(Level)**2
      efact = T_ceil / TfactSN / (gamma - 1.0) 
      heatfact = cfbe * t_fbi_code
c
      t_fbSN_code = t_fbSN / t0 / aexp(Level)**2
      t_fbSNi_code= 1.0 / t_fbSN_code
      heatfactSN  = cfbe *  t_fbSNi_code
ceverino10092007       heatfactSN  = cfbe / dtl(Level)

#ifdef FEEDBACK_SNIa
      t_SNIa_code = 1.d9 * t_SNIa / t0 / aexp(Level)**2
      R_SNIa_fact = RIaf * aexp(Level)**2
      cfbeIa = C_fbIa * CellVoli(Level) * aexp(Level)**2 
#endif
#endif

#ifdef KineticFCONSTANT
      Sigma_KinF= Sigma_KinFP/v0 *aexp(Level) ! code units
#endif
      gamma1 = gamma - 1.0
      Tfact = T_0 * gamma1 / 1.d4 ! T_0 must be in Kelvin
      dxi = 0.5 / CellSize(Level)
      e_min = T_min * a1**2 / T_0 / gamma1

      call Get_MaxLevelNow ()

CEVERINO04162007
cc      if (5785900.le.nsp(nspecies,2)) then
cc         icTrack= iFindCell ( MaxLevel , 
cc     &        x(5785900), y(5785900), z(5785900) )
cc      endif
c-------------------------------------------------------------------- 



      IF ( Level .eq. MinLevel ) THEN 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(t_e, dt_e, dt_e2, e_kin, ei, ei2, eratio, ei1, rhog2)
C$OMP+PRIVATE(acurr, ai, dummy0, dummy1, dummy2, dummy3, dummy4, dummy5)
C$OMP+PRIVATE(dummy6, rhogi, Tfact1, Tfact2, T_gas, coolrate, dtcool)
C$OMP+PRIVATE(e_small, Zdum, rhogl, Nbb)
C$OMP+PRIVATE(continue_cool)
C$OMP+PRIVATE(heatRate, idummy, istar, tp, tfbc,e_ceil)
C$OMP+PRIVATE(heatRateSNII, heatRateSNIa, xd, dN_SNIadtdM)
C$OMP+PRIVATE(substractSN)
C$OMP+PRIVATE(iNbb)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves 
c
c....       now take care of p*div V term and cooling/heating 
c
            t_e   = tbegin
            dummy0 = hvar(7,i1) - 1.0
            dummy5 = dummy0 * ref(i1)

            do while ( t_e .lt. tend )
              dt_e   = min ( dt3 , tend - t_e )
              hvar(8,i1) = hvar(8,i1) * (1.0 + dummy5 * dt_e)
              t_e = t_e + dt_e
            enddo ! end do while

c
c....       sync internal and total energy
c
            eratio = hvar(2,i1)
            e_kin = 0.5d0 * ( hvar(3,i1)**2 +
     &                        hvar(4,i1)**2 +
     &                        hvar(5,i1)**2   ) / hvar(1,i1)
            ei     = (dble(hvar(2,i1)) - e_kin)
            eratio = ( ei / eratio )
            if ( eratio .gt. eta2 ) then ! use e = E - rho * v**2 /2
              hvar(8,i1) = ei
            endif ! endif eratio
            hvar(7,i1) = gamma

#ifdef STARFORM
c
c...  Compute heating rates from SNII and SNIa 
c
#ifdef FEEDBACKenergy
            heatRate=0.0
            heatRateSNII=0.0
            e_ceil= efact * hvar(1,i1)
            if ( iCL(i1) .gt. nil ) then
               idummy = iCL(i1)  ! iparticle
               do while ( idummy .ne. nil )
                  if ( idummy .ge. nsp(nspecies,1) ) then ! a star
                     istar = idummy - nstar1
                     tp  = pt(idummy)                     
                     tfbc = tp - tbirth(istar) 
                     if ( tfbc .le. t_fb_code ) then ! young star
cc                        heatRateSNII= heatfact* pw0(istar)
cc                        heatRate=heatRate + heatRateSNII  
cccc                        if ( tfbc + pdt(idummy) .gt. t_fb_code ) then ! SN explosion
cccc                   heatRate=heatRate + 1.*cfbe/dtl(Level)* pw0(istar)
cccc                        endif

!                        if ( tfbc .ge. t_fb_code-t_fbSN_code) then ! SN explosion 
!                           heatRateSNII=heatfactSN*pw0(istar)
!                        else
                           heatRateSNII= heatfact* pw0(istar) ! continuous feedback
!                           if ( tfbc + pdt(idummy) .gt. t_fb_code ) then ! SN explosion
!                              heatRateSNII=heatfactSN*pw0(istar)
!                           endif
!                        endif
                        heatRate=heatRate + heatRateSNII  
                     endif
#ifdef FEEDBACK_SNIa
                     xd = t_SNIa_code / (tend - tbirth(istar))
                     dN_SNIadtdM = R_SNIa_fact * 
     &                    dexp(-xd*xd) * xd * dsqrt(xd)
                     heatRateSNIa=dN_SNIadtdM * cfbeIa* pw0(istar)
                     heatRate=heatRate + heatRateSNIa
#endif
                  endif
                  idummy         = iLL(idummy,1) ! next particle in the cell
               enddo            ! end do while
            endif
CEVERINO06222007
c            if (N_starsi1.gt.0) heatRate=Max(heatRate,  ! photoionization heating
c     & (1.d0*hvar(1,i1)) * hvar(1,i1)*a1*1.0d-1 * AL_SD )

!            if (heatRateSNII.gt.0.0)
!     &           heatRate= heatRate +
!     &           (1.d0*hvar(1,i1)) * hvar(1,i1)*a1 * AL_SD ! 10^-23 erg cm^3 s^-1 
c
#ifdef KineticFCONSTANT
            if (heatRateSNII.gt.0.0) then  ! CONSTANT KINETIC FEEDBACK
               call  NbAll0 ( i1 , iNbb )
       hvar(2,iNbb(1))=hvar(2,iNbb(1))+0.5*hvar(1,iNbb(1))*Sigma_kinF**2 
     &      -hvar(3,iNbb(1))*Sigma_kinF
       hvar(2,iNbb(2))=hvar(2,iNbb(2))+0.5*hvar(1,iNbb(2))*Sigma_kinF**2 
     &      +hvar(3,iNbb(2))*Sigma_kinF
       hvar(2,iNbb(3))=hvar(2,iNbb(3))+0.5*hvar(1,iNbb(3))*Sigma_kinF**2 
     &      -hvar(4,iNbb(3))*Sigma_kinF
       hvar(2,iNbb(4))=hvar(2,iNbb(4))+0.5*hvar(1,iNbb(4))*Sigma_kinF**2 
     &      +hvar(4,iNbb(3))*Sigma_kinF
       hvar(2,iNbb(5))=hvar(2,iNbb(5))+0.5*hvar(1,iNbb(5))*Sigma_kinF**2 
     &      -hvar(5,iNbb(5))*Sigma_kinF
       hvar(2,iNbb(6))=hvar(2,iNbb(6))+0.5*hvar(1,iNbb(6))*Sigma_kinF**2 
     &      +hvar(5,iNbb(6))*Sigma_kinF

       hvar(3, iNbb(1) )= hvar(3, iNbb(1)) - Sigma_kinF*hvar(1, iNbb(1)) 
       hvar(3, iNbb(2) )= hvar(3, iNbb(2)) + Sigma_kinF*hvar(1, iNbb(2))
       hvar(4, iNbb(3) )= hvar(4, iNbb(3)) - Sigma_kinF*hvar(1, iNbb(3))
       hvar(4, iNbb(4) )= hvar(4, iNbb(4)) + Sigma_kinF*hvar(1, iNbb(4))
       hvar(5, iNbb(5) )= hvar(5, iNbb(5)) - Sigma_kinF*hvar(1, iNbb(5))
       hvar(5, iNbb(6) )= hvar(5, iNbb(6)) + Sigma_kinF*hvar(1, iNbb(6))
            endif
#endif
#ifdef RadPressure
            RadPre(i1)=0.0 ! Reset Radiative Pressure in the leaves of Level
#endif
c       
#endif
#endif

#ifdef COOLING
            rhog2 = (1.d0*hvar(1,i1)) * hvar(1,i1)
            rhogi = 1.d0 / hvar(1,i1) 
#ifdef CLOUDY_COOLING
c....       rhog_code -> log10(n_H), [n_H]=[cm^-3]
            rhogl = log10(hvar(1,i1)) + fact_nH
#endif
            t_e   = tbegin
            acurr = a1
            Tfact1 = Tfact * rhogi
            e_small = max ( esmall * hvar(8,i1) , 
     &                      e_min * hvar(1,i1)) 
            continue_cool = .true. 

#ifdef METALCOOLING
#ifndef SOLARMETALL            
#ifdef ENRICH
            Zdum = hvar(izII,i1) 
#ifdef ENRICH_SNIa
            Zdum = Zdum + hvar(izIa,i1)
#endif
c
            Zdum = max ( Zdum , 1.d-10 )
            Zdum = log10(Zdum * rhogi * Zsolari)
CEVERINO04262007
C            Zdum= 0.0 ! solar metallicity
#endif
#endif
#endif      
            do while ( (t_e .lt. tend) .and. continue_cool ) 


              ai = 1.d0 / acurr
              dummy2 = ai * ai 
              Tfact2 = dummy2 * Tfact1
              T_gas = Tfact2 * hvar(8,i1)  
              call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 

              dummy1 = -rhog2 * acurr
              coolrate = coolrate * dummy1
CEVERINO04072008              dtcool = abs( hvar(8,i1) / coolrate ) 
              dtcool = abs( hvar(8,i1) / (coolrate+heatrate) ) 
              dtcool = dstep * dtcool
              dt_e   = min ( dtcool , tend - t_e )
c
              dt_e2  = dt_e * 0.5
#ifdef FEEDBACKenergy
              ei1 = max(hvar(8,i1)+(coolrate+heatrate) * dt_e2, e_small)
#endif
#ifndef FEEDBACKenergy
              ei1 = max(hvar(8,i1) + coolrate * dt_e2, e_small)
#endif
              T_gas = Tfact2 * ei1
              call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 
              coolrate = coolrate * dummy1
#ifdef FEEDBACKenergy
CEVERINO05022007
cc              if ((coolrate+heatRate) * dt_e.gt.efact * hvar(1,i1)) then
cc                 write (131,*) 'Cut in SN energy injection:', 
cc     &                (coolrate+heatRate) * dt_e  
cc              endif
C-------------------------------------------------------------------- 
c
CEVERINO09212007
cc              if (rhogl.gt.3.0) then
cc                 coolrate=0.0
cc              endif
c-------------------------------------------------------------------- 
                 dummy3 =min( (coolrate+heatRate) * dt_e, e_ceil )
ceverino07092007                 dummy3 =(coolrate+heatRate) * dt_e ! do not limit the energy release
                 hvar(8,i1) = hvar(8,i1) + dummy3
                 hvar(2,i1) = hvar(2,i1) + dummy3
#endif
#ifndef FEEDBACKenergy
              dummy3 =coolrate * dt_e
              hvar(8,i1) = hvar(8,i1) + dummy3
              hvar(2,i1) = hvar(2,i1) + dummy3
#endif
c
              if ( hvar(8,i1) .lt. e_small ) then 
                continue_cool = .false.
              endif
              hvar(8,i1) = max(dble(hvar(8,i1)), e_small)
              hvar(2,i1) = max(dble(hvar(2,i1)), e_small)
                      
              t_e = t_e + dt_e 
              acurr = a1 + afact * (t_e - tl(Level))
ceverino10042007
C              if (abs(dummy3/hvar(8,i1)).lt.1.e-8) then
c                 write (124,*)  i1,dummy3, hvar(8,i1)
C                 continue_cool = .false.
C              endif
           enddo                ! end do while

#ifdef FEEDBACKenergy
ceverino1062007
           if ( hvar(8,i1) .gt. e_ceil ) then  ! ceiling of temperature
              substractSN=  hvar(8,i1) -  e_ceil
              hvar(8,i1) = hvar(8,i1) - substractSN
              hvar(2,i1) = hvar(2,i1) - substractSN
           endif
#endif

c            endif
#endif
          endif
          
        enddo ! end i1

      ELSE

        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,iChild,i1, ici, ici0, Nbb, ei2)
C$OMP+PRIVATE(t_e, dt_e, dt_e2, e_kin, ei, eratio, ei1)
C$OMP+PRIVATE(rhog2, rhogi, rhogl)
C$OMP+PRIVATE(acurr, ai, dummy0, dummy1, dummy2, dummy3, dummy4, dummy5)
C$OMP+PRIVATE(Tfact1, Tfact2, T_gas, coolrate, coolcomp, dtcool, dummy6)
C$OMP+PRIVATE(e_flux, e_pdV, e_cool, e_small, dtcoolmin, Zdum)
C$OMP+PRIVATE(continue_cool)
C$OMP+PRIVATE(heatRate, idummy, istar, tp, tfbc,e_ceil)
C$OMP+PRIVATE(heatRateSNII, heatRateSNIa, xd, dN_SNIadtdM)
CEVERINOC$OMP+PRIVATE(aM, N_starsi1,N_sumInx)
CEVERINOC$OMP+PRIVATE(iHYDRO,GravGrad,PresGrad,absGravGrad, absPresGrad)
C$OMP+PRIVATE(substractSN)
C$OMP+PRIVATE(iNbb)
        do ic = 1 , nLevel
          icell = iSelect(ic) - 1
          ici0 = (ic-1)*nchild
          do iChild = 1 , 8 
            i1 = icell + iChild
            ici = ici0 + iChild
            if ( iOctCh(i1) .eq. nil ) then  ! only leaves 
c
c....         now take care of p*div V term and cooling/heating 
c
c            write(*,*) 'i1=',i1
              t_e   = tbegin
              dummy0 = hvar(7,i1) - 1.0
              dummy5 = dummy0 * ref(i1)

              do while ( t_e .lt. tend )
                dt_e   = min ( dt3 , tend - t_e )
                hvar(8,i1) = hvar(8,i1) * (1.0 + dummy5 * dt_e)
                t_e = t_e + dt_e
              enddo ! end do while

c
c....         sync internal and total energy
c
              eratio = hvar(2,i1)
              e_kin = 0.5d0 * ( hvar(3,i1)**2 +
     &                          hvar(4,i1)**2 +
     &                          hvar(5,i1)**2   ) / hvar(1,i1)
              ei     = (dble(hvar(2,i1)) - e_kin)
              eratio = ( ei / eratio )

              if ( eratio .gt. eta2 ) then ! use e = E - rho * v**2 /2
                hvar(8,i1) = ei
              endif ! endif eratio
              hvar(7,i1) = gamma

#ifdef STARFORM
c
c...  Compute heating rates from SNII and SNIa 
c
#ifdef FEEDBACKenergy
            heatRate=0.0
            heatRateSNII=0.0
            e_ceil= efact * hvar(1,i1)
CEVERINO04162007
ccc            if ( t.gt.0.1950000E-03.AND.Level.eq.1) then
cc                  if (5785900.le.nsp(nspecies,2)) then
cc               if(i1.eq.icTrack) then
cc                  aM=0.0
cc                  N_starsi1=0
cccc               N_sumInx=0
cc               endif
ccc                  endif
cc            endif
C-------------------------------------------------------------------- 
            if ( iCL(i1) .gt. nil ) then
               idummy = iCL(i1)  ! iparticle
               do while ( idummy .ne. nil )
                  if ( idummy .ge. nsp(nspecies,1) ) then ! a star
                     istar = idummy - nstar1
                     tp  = pt(idummy)                     
                     tfbc = tp - tbirth(istar) 
                     if ( tfbc .le. t_fb_code ) then ! young star
cc                        heatRateSNII= heatfact* pw0(istar) 
cc                        heatRate=heatRate + heatRateSNII
cccc                        if ( tfbc + pdt(idummy) .gt. t_fb_code ) then ! SN explosion
cccc                   heatRate=heatRate + 1.*cfbe/dtl(Level)* pw0(istar)
cccc                        endif

!                        if ( tfbc .ge. t_fb_code-t_fbSN_code) then ! SN explosion
!                           heatRateSNII=heatfactSN*pw0(istar)
!                        else
                           heatRateSNII= heatfact* pw0(istar) ! continuous feedback
!                           if ( tfbc + pdt(idummy) .gt. t_fb_code ) then ! SN explosion
!                              heatRateSNII=heatfactSN*pw0(istar)
!                           endif
!                        endif
                        heatRate=heatRate + heatRateSNII  

CEVERINO04162007
ccc                        if ( t.gt.0.1950000E-03.AND.Level.eq.1) then
cc                           if (5785900.le.nsp(nspecies,2)) then
cc                           if(i1.eq.icTrack) then
cc                              aM=aM +  pw0(istar)*aM0
cc                              N_starsi1=N_starsi1 + 1
cccc                     N_sumInx= N_sumInx+istar
ccc                           endif
cc                           endif
cc                        endif
C-------------------------------------------------------------------- 

                     endif
#ifdef FEEDBACK_SNIa
                     xd = t_SNIa_code / (tend - tbirth(istar))
                     dN_SNIadtdM = R_SNIa_fact * 
     &                    dexp(-xd*xd) * xd * dsqrt(xd)
                     heatRateSNIa=dN_SNIadtdM * cfbeIa* pw0(istar)
                     heatRate=heatRate + heatRateSNIa
#endif
                  endif
                  idummy         = iLL(idummy,1) ! next particle in the cell
               enddo            ! end do while
            endif
CEVERINO06222007
!            if (heatRateSNII.gt.0.0)
!     &           heatRate= heatRate +
c     &           heatRate=Max(heatRate, ! photoionization heating
!     &           (1.d0*hvar(1,i1)) * hvar(1,i1)*a1 * AL_SD ! 10^-23 erg cm^3 s^-1 
c
#ifdef KineticFCONSTANT
            if (heatRateSNII.gt.0.0) then ! CONSTANT KINETIC FEEDBACK
               call  NbAllL ( i1 , iNbb )
       hvar(2,iNbb(1))=hvar(2,iNbb(1))+0.5*hvar(1,iNbb(1))*Sigma_kinF**2 
     &      -hvar(3,iNbb(1))*Sigma_kinF
       hvar(2,iNbb(2))=hvar(2,iNbb(2))+0.5*hvar(1,iNbb(2))*Sigma_kinF**2 
     &      +hvar(3,iNbb(2))*Sigma_kinF
       hvar(2,iNbb(3))=hvar(2,iNbb(3))+0.5*hvar(1,iNbb(3))*Sigma_kinF**2 
     &      -hvar(4,iNbb(3))*Sigma_kinF
       hvar(2,iNbb(4))=hvar(2,iNbb(4))+0.5*hvar(1,iNbb(4))*Sigma_kinF**2 
     &      +hvar(4,iNbb(3))*Sigma_kinF
       hvar(2,iNbb(5))=hvar(2,iNbb(5))+0.5*hvar(1,iNbb(5))*Sigma_kinF**2 
     &      -hvar(5,iNbb(5))*Sigma_kinF
       hvar(2,iNbb(6))=hvar(2,iNbb(6))+0.5*hvar(1,iNbb(6))*Sigma_kinF**2 
     &      +hvar(5,iNbb(6))*Sigma_kinF

       hvar(3, iNbb(1) )= hvar(3, iNbb(1)) - Sigma_kinF*hvar(1, iNbb(1)) 
       hvar(3, iNbb(2) )= hvar(3, iNbb(2)) + Sigma_kinF*hvar(1, iNbb(2))
       hvar(4, iNbb(3) )= hvar(4, iNbb(3)) - Sigma_kinF*hvar(1, iNbb(3))
       hvar(4, iNbb(4) )= hvar(4, iNbb(4)) + Sigma_kinF*hvar(1, iNbb(4))
       hvar(5, iNbb(5) )= hvar(5, iNbb(5)) - Sigma_kinF*hvar(1, iNbb(5))
       hvar(5, iNbb(6) )= hvar(5, iNbb(6)) + Sigma_kinF*hvar(1, iNbb(6))
            endif
#endif
#ifdef RadPressure
            RadPre(i1)=0.0 ! Reset Radiative Pressure in the leaves of Level
#endif
c
#endif
#endif

#ifdef COOLING
              rhog2 = (1.d0*hvar(1,i1)) * hvar(1,i1)
              rhogi = 1.d0 / hvar(1,i1) 
#ifdef CLOUDY_COOLING
c....         rhog_code -> log10(n_H), [n_H]=[cm^-3]
              rhogl = log10(hvar(1,i1)) + fact_nH
#endif
              t_e   = tbegin
              acurr = a1
              Tfact1 = Tfact * rhogi 
              e_small = max ( esmall * hvar(8,i1) , 
     &                        e_min * hvar(1,i1)) 
              continue_cool = .true. 

#ifdef METALCOOLING
#ifndef SOLARMETALL            
#ifdef ENRICH
              Zdum = hvar(izII,i1) ! SNII metal density
#ifdef ENRICH_SNIa
              Zdum = Zdum + hvar(izIa,i1)  ! SNIa metal density 
#endif
c
              Zdum = max ( Zdum , 1.d-10 )
              Zdum = log10 ( Zdum * rhogi * Zsolari ) 
CEVERINO04262007
C            Zdum= 0.0 ! solar metallicity
#endif
#endif
#endif

!              write(132,*) i1, rhogi, Zdum, Zsolari

              do while ( (t_e .lt. tend) .and. continue_cool )
                ai = 1.d0 / acurr
                dummy2 = ai * ai 
                Tfact2 = dummy2 * Tfact1
                T_gas = Tfact2 * hvar(8,i1)  
                call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 

                dummy1 = -rhog2 * acurr
                coolrate = coolrate * dummy1
CEVERINO04072008                dtcool = abs(hvar(8,i1) / coolrate) 
                dtcool = abs(hvar(8,i1) / (coolrate+heatrate) ) 
                dtcool = dstep * dtcool
c
                dt_e   = min ( dtcool , tend - t_e )
c
                dt_e2  = dt_e * 0.5
#ifdef FEEDBACKenergy
              ei1 = max(hvar(8,i1)+(coolrate+heatrate) * dt_e2, e_small)
#endif
#ifndef FEEDBACKenergy
                ei1 = max(hvar(8,i1) + coolrate * dt_e2, e_small)
#endif
                T_gas = Tfact2 * ei1
                call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 
                coolrate = coolrate * dummy1
#ifdef FEEDBACKenergy   
cc              if ((coolrate+heatRate) * dt_e.gt.efact * hvar(1,i1)) then
cc                   write (131,*) 'Cut in SN energy injection:', i1, 
cc     &                  (coolrate+heatRate) * dt_e  
cc                endif
c
CEVERINO09212007
cc              if (rhogl.gt.3.0) then
cc                 coolrate=0.0
cc              endif
c-------------------------------------------------------------------- 
cc                if (T_gas.gt.10.0.AND.rhogl.ge.1.0 ) then
cc                 write (124,133) i1,  10**rhogl, T_gas*1e4,
cc     &                  coolrate, heatrate, coolrate+heatrate
cc 133               format (I10, 5(e15.8,4x))
cc                endif
ceverino10062007
cc                if ( (coolrate+heatRate) * dt_e.gt. e_ceil
cc     &               .AND.rhogl.ge.-1.0 ) then
cc                   write (123,*) '  limited energy release', i1,
cc     &                  (coolrate+heatRate) * dt_e,  e_ceil, 
cc     &                  rhogl, T_gas, coolrate, heatrate  
cc                endif
c-------------------------------------------------------------------- 
                dummy3 =min( (coolrate+heatRate) * dt_e, e_ceil )
ceverino07092007                   dummy3 = (coolrate+heatRate) * dt_e ! do not limit energy release

                hvar(8,i1) = hvar(8,i1) + dummy3
                hvar(2,i1) = hvar(2,i1) + dummy3
#endif
#ifndef FEEDBACKenergy
                dummy3 = coolrate * dt_e
                hvar(8,i1) = hvar(8,i1) + dummy3
                hvar(2,i1) = hvar(2,i1) + dummy3
#endif
CEVERINO04162007
cc                if ( t.gt.0.1950000E-03.AND.Level.eq.1) then
ccc                   if (3313.le.nsp(nspecies,2)) then
cc                      if(i1.eq.icTrack) then
cc                         write(128,132)  t
cc     &                        ,10**rhogl,T_gas*1.d4 , coolrate, heatRate
cc                      endif
ccc                   endif
cc                endif
c-------------------------------------------------------------------- 

                if ( hvar(8,i1) .lt. e_small ) then 
                  continue_cool = .false.
                endif
                hvar(8,i1) = max(dble(hvar(8,i1)), e_small)
                hvar(2,i1) = max(dble(hvar(2,i1)), e_small)

                t_e = t_e + dt_e 
                acurr = a1 + afact * (t_e - tl(Level))

ceverino10042007
!              if ( abs(dummy3/hvar(8,i1)).lt.1.e-7.AND.
!     &               abs(dummy3/hvar(2,i1)).lt.1.e-7) then
c                 write (124,*)  i1,dummy3/hvar(8,i1), dummy3, hvar(8,i1)
!                 continue_cool = .false.
!              endif

              enddo ! end do while
#ifdef FEEDBACKenergy
ceverino1062007
           if ( hvar(8,i1) .gt. e_ceil ) then  ! ceiling of temperature
              substractSN=  hvar(8,i1) -  e_ceil
              hvar(8,i1) = hvar(8,i1) - substractSN
              hvar(2,i1) = hvar(2,i1) - substractSN
           endif
#endif

CEVERINO04162007
cc                   if (5785900.le.nsp(nspecies,2)) then
cc                      if(i1.eq.icTrack) then
ccC     Hydrostatic equilibrium
cc                         absGravGrad= 0.0000
cc                         absPresGrad= 0.0000
cc                         do iHYDRO=1,3
cc                            ineigm = iNb ( i1 , iHYDRO*2-1 )
cc                            ineigp = iNb ( i1 , iHYDRO*2 )
ccc... grav. acceleration:
cc         GravGrad= ( var(2,ineigp)-var(2,ineigm) )/(2.0*CellSize(Level))
ccc... pressure gradient:
cc         PresGrad=(hvar(6,ineigp)-hvar(6,ineigm))/(2.0*CellSize(Level))
cc         PresGrad= PresGrad/ hvar(1,i1)
ccc... forces equilibrium min & max:
ccC         eqForcesMax =max(eqForcesMax, abs(PressGrad)/abs(GravGrad))
ccC         eqForcesMin =min(eqForcesMin, abs(PressGrad)/abs(GravGrad))
ccc
ccc... modulus:
cc         absGravGrad=absGravGrad + GravGrad**2
cc         absPresGrad=absPresGrad + PresGrad**2
ccc
ccc                            write (130,131) t,iHYDRO,GravGrad, PresGrad,
ccc     &       abs(PresGrad)-abs(GravGrad) , abs(PresGrad)/abs(GravGrad)
cc                             enddo          
ccc
cc                  write(129,127)  t, Level, i1, 
cc     &                        x(5785900), y(5785900), z(5785900)
ccc                         write(129,127)  t, 31.2, 2.7, 32.3
cc     &                        ,10**rhogl,T_gas*1.d4 , coolrate, heatRate
cc     &                        , aM
cc     &                        , N_starsi1
cc     &                        ,sqrt(absGravGrad), sqrt(absPresGrad)
cc                      endif
ccc                   endif
cc                   endif
c-------------------------------------------------------------------- 

c              endif
#endif

            endif ! end if iOctCh            
          enddo
        enddo
      ENDIF
      return
CEVERINO04162007
cc 127  format(e15.8,2x, 2(I9), 3(e15.8,2x), f9.3,2x,3(e15.8,2x),e11.4,2x 
cc     &     ,1(I9,2x), 2(e15.8,2x) )
cc 131  format (e15.8,2x,I9, 'GravGrad', g20.6,2x, 'PressGrad', g20.6, 2x, 
cc     & 'HydroS eq Diff '    , g20.6,2x, 'HydroS eq ratio ', g20.6)
cc 132  format (5(e15.8,2x))
c-------------------------------------------------------------------- 
      end
