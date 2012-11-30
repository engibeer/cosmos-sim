c
c....  routines handling the starformation processes 
c
c
c     -----------------------------------------------
      subroutine Star_Formation ( Level , istarform,  dt_SF )
c     -----------------------------------------------
c
c     form stars in eligible cells. The eligibility criteria, 
c     the minimum gas density rho_SF and maximum gas temperature T_SF
c     (in code units) 
c     must be be specified during the initialization 
c
c     the stars are formed in cells where hvar(1,ic) > rho_SF
c                                         hvar(6,ic) / hvar(1,ic) < T_SF
c                                         hvar(1,ic) * Om0 / Omb0 > deltavir
c
c     if level's step is smaller than a specified minimum time of SF
c     then SF is performed only during certain steps - this "frequency"
c     is specified by istarform (see ART_Step for details)

#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      parameter (G_cgs = 6.67e-8)
      parameter (eps_ff = 0.10 )
c
c.... do not allow stellar particles with mass > fraccell * mcell
c
      parameter ( fraccell = 0.667 ) 
c
c.... allow star formation only in the virialized objects 
c
      parameter ( deltavir = 200. ) 

c
CEVERINO09042007:  STOCHASTIC STAR FORMATION
c
      parameter (fstars=0.01)

      integer Level, istarform  
      real*8 n_H_min, dmst_min, dm_star, dum, dmcell, dmcell_tot
      real*8 t_ff, t_ffcu, rho_cgs, SF_period, dt_SF

c========= STG:
      logical formedstar

      formedstar = .false.
c============

      if ( Level .eq. MinLevel ) return
c
c.... mass factor 
c
c======== STG: new SF law
      fmass = aexp(Level)**(5.0-3.0*alpha_SF) * C_SFR * 
     &        dtl(Level) * istarform / 2.0**(3.0*Level)
c      fmass = aexp(Level)**2 * (0.05 * t0 * dtl(Level)) 
c     &         * istarform / 2.0**(3.0*Level)
c      fmass = (0.05 / dtl(Level)) / 2.0**(3.0*Level)
c==============

CEVERINO04072008
c      fmassMax=aexp(Level)**(5.0-3.0*alpha_SF) * C_SFR * 
c     &        dt_SF / 2.0**(3.0*Level)
c      fmass=max(fmass, fmassMax)
c      write (*,*)  'SF at L=',Level, fmass
c-------------------------------------------------------------------- 


      T_SF_max = T_SF / T_0 * aexp(Level)**2
      rho_SF_min = rho_SF_fact * aexp(Level)**3
c      rho_SF_min = rho_SF_fact
      dmst_min = dm_star_min / aM0 
      n_H_min = 10.0 *  ! <- n_crit in cm^-3
     &     rho_SF_fact / rho_SF * aexp(Level)**3

      cellfrac = 2.**(-3.0*Level) * fraccell 

c====== STG:
      fact_nH = 1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p)/aexp(Level)**3
c==========

c#ifdef DEBUGsf
c      call Open_ASCII_File ( 13 , 'debug ' , 'append ' ) 
c      write(13,*) 'Forming stars on Level=',Level
c      write(13,*) ' rho_SF_min=',rho_SF_min
c      write(13,*) 'dmst_min=',dmst_min
c      close(13)
c#endif
 
      IF ( Level .eq. MinLevel ) THEN
c
         do idc = 1 , ncell0
            if ( iOctCh(idc) .eq. nil ) then ! only leaves 
              if ( hvar(1,idc) .gt. rho_SF_min ) then ! only eligible cells
                Tcell = hvar(6,idc) / hvar(1,idc)
                dencell = hvar(1,idc) * Om0 / Omb0 
                if ( Tcell .lt. T_SF_max .and. 
     &               dencell .gt. deltavir ) then 
                  
CEVERINO09042007
                   coin=Randd(Nseed)
!               RandomTH=Min( sqrt(hvar(1,idc)/rho_SF_min)/10.0, 1.) ! for 10%SF at rho_SF_min
!                    RandomTH=Min( (hvar(1,idc)/rho_SF_min)/100.0,1.) ! for 1% SF at rho_SF_min
!                    RandomTH=Min(sqrt((hvar(1,idc)/rho_SF_min))/100.,1.) ! for 1% SF at rho_SF_min 1/2 law
                   RandomTH=Min(sqrt((hvar(1,idc)/rho_SF_min)/1000.),1.) ! for 3% SF at rho_SF_min 1/2 law
CEVERINO07012009
                    RandomTH=Min(RandomTH,0.2)  ! Max of 20%
                    if (coin.lt.RandomTH) then
c                   if (coin.lt.fstars) then

c
c....             compute the amount of stellar mass to create
c 
                      dmcell = hvar(1,idc) * cellfrac 
                      dum = fmass * hvar(1,idc)**alpha_SF
ccCEVERINO09262007
cc                      if ( hvar(1,idc)/rho_SF_min .gt. 100.0) then
cc                         dum =  dum * 100.
cc                      endif
c-------------------------------------------------------------------- 
CEVERINO09202007
cc                      coin=Randd(Nseed)
cc               RandomTH=Min( (hvar(1,idc)/rho_SF_min)**(1./3.)/10.0, 1.)
cc                      if (coin.lt.RandomTH) then
cc                         dum =  dum * 100.
cc                      endif
c-------------------------------------------------------------------- 
                      dm_star = max ( dum , dmst_min ) 
                      dm_star = min ( dm_star , dmcell ) 
CEVERINO10232006                  dm_star= min (dum, dmcell)

CEVERINO03102008                  IF (dm_star.GE.dmst_min)
                   call Create_Stellar_Particle (idc, Level, dm_star)

c                  endif
                   endif
                  
                endif
              endif
            endif 
         enddo
      ELSE
c
        call Select_Cells ( Level , nLevel ) 
        do ic1 = 1 , nLevel 
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idc = icell + ic2
            if ( iOctCh(idc) .eq. nil ) then ! only leaves 
              if ( hvar(1,idc) .gt. rho_SF_min ) then ! only eligible cells
                Tcell = hvar(6,idc) / hvar(1,idc)
                dencell = hvar(1,idc) * Om0 / Omb0 
                if ( Tcell .lt. T_SF_max .and. 
     &               dencell .gt. deltavir ) then 
c================ STG: removed random SF
cCEVERINO09042007
c                   coin=Randd(Nseed)
c!               RandomTH=Min( sqrt(hvar(1,idc)/rho_SF_min)/10.0, 1.) ! for 10%SF at rho_SF_min
c!                    RandomTH=Min( (hvar(1,idc)/rho_SF_min)/100.0,1.) ! for 1% SF at rho_SF_min
c!                    RandomTH=Min(sqrt((hvar(1,idc)/rho_SF_min))/100.,1.) ! for 1% SF at rho_SF_min 1/2 law
c                   RandomTH=Min(sqrt((hvar(1,idc)/rho_SF_min)/1000.),1.) ! for 3% SF at rho_SF_min 1/2 law    
cCEVERINO07012009
c                    RandomTH=Min(RandomTH,0.2) ! max of 20%
c                   if (coin.lt.RandomTH) then
c                   if (coin.lt.fstars) then

c=======================
c....             compute the amount of stellar mass to create
c                   
                      dmcell = hvar(1,idc) * cellfrac  ! max gas mass available in cell (2/3)m_g
c========= STG: new SF law
c                      dm_star = fmass * hvar(1,idc)**alpha_SF  ! SF law
                      dmcell_tot = hvar(1,idc) / 2.0**(3.0*Level)
	              rho_cgs = hvar(1,idc) * fact_nH * 1.67d-24
                      t_ff = sqrt(3.0 * 3.1416 / (32.0 * G_cgs 
     &                       * rho_cgs)) / 3.1536d7  ! free-fall time in yrs
		      t_ffcu = t_ff / (t0 * aexp(Level)**2)
                      SF_period = max(dt_SF,dtl(0)) ! SF occurs once per period (see ART_Step.F)
		  dm_star = (eps_ff * SF_period/t_ffcu) * dmcell_tot  ! eps_sf fraction of gas mass per t_ff
		      call Open_ASCII_File ( 13 , 'debug ' , 'append ' ) 
                      write(13,*) 'in SF loop - Level:', Level
		      write(13,*) idc, Tcell*T_0/aexp(Level)**2, 
     &                  rho_cgs/1.67d-24, 
     & 			t_ff, dtl(Level)*t0*aexp(Level)**2,
     &                  dm_star*aM0, 
     &                  dmcell_tot*aM0,
     &                  dm_star/dmcell_tot
                      write(13,*) 'timescales: ',
     &                  dt_SF*(t0*aexp(Level)**2),
     &                  dtl(0)*(t0*aexp(Level)**2),
     &                  SF_period*(t0*aexp(Level)**2)
                      close(13)
cc		      pause
c===================
CEVERINO09262007
cc                      if ( hvar(1,idc)/rho_SF_min .gt. 100.0) then
cc                         dm_star =  dm_star * 100.
cc                      endif
c-------------------------------------------------------------------- 
CEVERINO09202007
cc                       coin=Randd(Nseed)
cc               RandomTH=Min( (hvar(1,idc)/rho_SF_min)**(1./3.)/10.0, 1.)
cc                       if (coin.lt.RandomTH) then
cc                          dm_star =  dm_star * 100.
cc                       endif
c-------------------------------------------------------------------- 
c.....          choose particle mass between dmst_min and dmcell:
c                      dm_star = max ( dm_star , dmst_min ) 
c                      dm_star = min ( dm_star , dmcell ) 
CEVERINO10232006                  dm_star = fmass * hvar(1,idc)**alpha_SF
check
c                      if (dm_star.gt.dmcell) 
c    &                      write(22,23) t, Level,idc,dmcell,dm_star
c 23                   format(g16.6, I3, I15, g16.3, g16.3)

c                       write(22,23) t, Level, idc, 
c     &                     hvar(1,idc)/rho_SF_min, RandomTH, coin, 
c     &                     dmcell, dm_star,  Tcell/T_SF_max
c 23                    format(g16.6, I3, I15, F9.2, F9.3, F9.3,
c     &                      g16.3, g16.3, g16.3)
c--------------------------------------------------------------------  
CEVERINO03102008                 IF (dm_star.GE.dmst_min)
                      call Create_Stellar_Particle (idc, Level, dm_star)
                  
                      if ( CellSize(0) .lt. 0 ) then 
                         write(*,*) Level, idc, dmcell, dm_star,CellSize
                      endif
c========= STG:
		      formedstar = .true.
c=============
c                  endif
c                   endif
                  
                endif
              endif
            endif 
          enddo
        enddo
      ENDIF
c
#ifdef DEBUGsf
      if (formedstar .eqv. .true.) then
      	call Open_ASCII_File ( 13 , 'debug ' , append ) 
      	write(13,*) 'In Star_Formation...'
      	dtlev = dtl(Level) * t0 * aexp(Level)**2
      	dtsf = dtlev * istarform
      	write(13,*) 'step =', istep, 'L =',Level, ' t=',tl(Level),
     &      ' dt=',dtlev, ' dt_SF =', dtsf,' [yrs]' , 
     &      'istarform= ', istarform
      	write(13,*) 'nstars =',nsp(nspecies,2)-nsp(nspecies,1)+1
      	pwmin = 1.e20
      	pwmax = -pwmin
      	nstar1 = nsp(nspecies,1) - 1 
      	do ip = nsp(nspecies,1) , nsp(nspecies,2)
        	pwmin = min(pwmin,pw(ip))
        	pwmax = max(pwmax,pw(ip))
c        write(13,*) ip, x(ip), y(ip), z(ip), 
c     &    pw(ip)*aM0, tbirth(ip-nstar1)
      	enddo
      	pwmin = pwmin * aM0
      	pwmax = pwmax * aM0
      	write(13,*) 'Mstar_min,max =',pwmin,pwmax,' [Msun]'
      	close(13)
	formedstar = .false.
      endif
#endif      
c
c======= STG:
c	stop
c===========
      return
      end

c     ----------------------------------------------------------
      subroutine Create_Stellar_Particle ( ic , Level, dm_star )
c     ----------------------------------------------------------
c
c     create stellar particle with all the attributes and insert it in the linked
c     lists
c
c     stars is assumed to be the specie with the largest index
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_dnew.h'

      integer ic, Level  
      real xcell, ycell, zcell
      real*8 dm_star, rho_star

      integer lspecies(nspec)
      equivalence (lspecies(1),extras(11))
#ifdef RUNAWAY
      parameter ( V_RunS= 17.0 )  ! 17
      parameter ( coin_RunS=0.33) ! 0.33)
      v_fact=aexp(Level)/v0  
#endif

      if ( np + 1 .gt. npmax ) then 
        call Open_ASCII_File ( iOErr , errorfile , sequent ) 
        write(iOErr,*)'1  *error in Create_Stellar_Particle: np>npmax'
        write(iOErr,*)'   * np =',np+1,' npmax =',npmax
        write(iOErr,*)'   * increase npmax in a_setup.h and rerun'
        close (iOErr)
        stop
      endif

      np = np + 1
      lsp(nspecies) = lsp(nspecies) + 1
      lspecies(nspecies) = lsp(nspecies)
      nsp(nspecies,2) = nsp(nspecies,2) + 1
      pw(np) = dm_star

      call Ps ( ic , xcell, ycell, zcell )
      x(np) = xcell
      y(np) = ycell
      z(np) = zcell
      rhor = 1.0 / hvar(1,ic)
#ifdef RUNAWAY
      coin=Randd(Nseed)
      if (coin.lt.coin_RunS) then
         Vrun= -V_RunS*LOG(Randd(Nseed))*v_fact ! cu
!         Xram=sqrt(-2.0*LOG(Randd(Nseed))) ! Gauss
!         Yram=cos (pi*Randd(Nseed)/2.0)
!         Vrun=V_RunS*Xram*Yram*v_fact ! cu
         URUN=1.1
         do while (URUN.gt.1.00)
            UxRUN=2.*Randd(Nseed)-1.0
            UyRUN=2.*Randd(Nseed)-1.0
            UzRUN=2.*Randd(Nseed)-1.0
            URUN= UxRUN**2+UyRUN**2+UzRUN**2
         enddo
         URUN=sqrt(URUN)
         vx(np) = hvar(3,ic) * rhor + Vrun*UxRUN/URUN      
         vy(np) = hvar(4,ic) * rhor + Vrun*UyRUN/URUN  
         vz(np) = hvar(5,ic) * rhor + Vrun*UzRUN/URUN
CEVERINO04262007
cc         write (139,128) t, Vrun/ v_fact,
cc     &        hvar(3,ic) * rhor/ v_fact, hvar(4,ic) * rhor/ v_fact, 
cc     &        hvar(5,ic) * rhor/ v_fact,
cc     &        Vrun*UxRUN/URUN/ v_fact, Vrun*UyRUN/URUN/ v_fact,
cc     &        Vrun*UzRUN/URUN/ v_fact
cc 128     format ( 8(g16.3,2x) )
c-------------------------------------------------------------------- 
      else
         vx(np) = hvar(3,ic) * rhor 
         vy(np) = hvar(4,ic) * rhor
         vz(np) = hvar(5,ic) * rhor
      endif         
#endif
#ifndef RUNAWAY
      vx(np) = hvar(3,ic) * rhor 
      vy(np) = hvar(4,ic) * rhor
      vz(np) = hvar(5,ic) * rhor
#endif
      pt(np)  = tl(Level) 
      pdt(np) = dtl(Level)

      nstar = np - nsp(nspecies,1) + 1
      tbirth(nstar) = tl(Level)   ! mark the birth time
      pw0(nstar) = dm_star        ! save the initial mass of the particle

#ifdef ENRICH
      zstII(nstar) = hvar(izII,ic) / hvar(1,ic)
#endif
#ifdef ENRICH_SNIa
      zstIa(nstar) = hvar(izIa,ic) / hvar(1,ic)
#endif

      call LL_Insert ( ic , np )
      iPL(np) = Level

      rho_star = dm_star * 2.**(3.*Level)
      hvar(1,ic) = hvar(1,ic) - rho_star
      if ( hvar(1,ic) .le. 0.0 ) then 
        write(*,*) '* error in Create_Stellar_Particle: hvar(1,ic)<=0'
        write(*,*) '* ic, L, hvar(1,ic), dm_star =',
     &              ic,Level, hvar(1,ic), dm_star
        stop
      endif

      rhofact = rhor * hvar(1,ic)
      hvar(2,ic) = hvar(2,ic) * rhofact 
      hvar(3,ic) = hvar(3,ic) * rhofact 
      hvar(4,ic) = hvar(4,ic) * rhofact 
      hvar(5,ic) = hvar(5,ic) * rhofact 
      hvar(6,ic) = hvar(6,ic) * rhofact 
      hvar(8,ic) = hvar(8,ic) * rhofact 
#ifdef ENRICH
      hvar(izII,ic) = max ( 1.d-17, 
     &                      hvar(izII,ic) - zstII(nstar)*rho_star)
#endif
#ifdef ENRICH_SNIa
      hvar(izIa,ic) = max ( 1.d-17, 
     &                      hvar(izIa,ic) - zstIa(nstar)*rho_star)
#endif
      
      return
      end

c     ---------------------
      subroutine SF_Init ()
c     ---------------------
c
c     initialize starformation/feedback/enrichment related constants
c
#     include "a_def.h"
      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'

      real*8 dummy1, dummy2, dummy3, dummy4, dummy5, dumst
      real*8 fp(10), Aprime
      real*8 INTEGRATE, f_IMF, fm_IMF, f_IMF_plaw, fm_IMF_plaw
      real*8 fej_IMF, fej_IMF_plaw, f_SNIa_func
      real*8 aN_SNII, An_IMF
      external f_IMF, fm_IMF, f_IMF_plaw, fm_IMF_plaw
      external fej_IMF, fej_IMF_plaw, f_SNIa_func

#ifdef ENRICH
      if ( nhvar .lt. 9 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init : ENRICH is defined in a_def.h but'
          write(13,*)' nhvar =', nhvar,' < 9 in a_setup.h'
          write(13,*)'stopping...'
          close ( 13 ) 
          stop          
      endif
#endif      
#ifdef ENRICH_SNIa
#ifdef ENRICH
      if ( nhvar .lt. 10 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init: ENRICH and ENRICH_SNIa are defined in a_def.h'
          write(13,*)' but nhvar =', nhvar,' < 10 in a_setup.h'
          write(13,*)'stopping...'
          close ( 13 ) 
          stop          
      endif
#endif
#ifndef ENRICH
      if ( nhvar .lt. 9 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init : ENRICH_SNIa is defined in a_def.h'
          write(13,*)' but nhvar =', nhvar,' < 9 in a_setup.h'
          write(13,*)'stopping...'
          close ( 13 ) 
          stop          
      endif
#endif
#endif      
c
c.... below t0 must be in yrs, rho0 in Msun / Mpc^-3 (no h!) 
c
c
c.... starformation efficiency factor c_SFR
c
      C_SFR = eps_SF * 2.5*10.0**(6.0 - 16.0*alpha_SF) * 
     &        t0 * rho0**(alpha_SF-1.0)

      E_0 = 1.38d58 * r0**5 * Om0**2 / hubble  ! [ergs]
      write(*,*) 'E_0, r0, Om0, aM0, sun_mass =',
     &     E_0, r0, Om0, aM0, sun_mass  

      IF ( a_IMF .gt. 0.0 ) THEN  ! if power-law
        if ( a_IMF .eq. 2.0 ) then 
          Aprime = 1.0 / log(aM_stu/aM_stl)
        else
          Aprime = (2.0 - a_IMF) / 
     &             (aM_stu**(2.0-a_IMF) - aM_stl**(2.0-a_IMF)) 
        endif
        An_IMF = Aprime ! IMF normalization
        if ( a_IMF .eq. 1.0 ) then 
          Aprime = Aprime * log(aM_stu/aM_SNII)
        else
          Aprime = Aprime /
     &             (1.0 - a_IMF) * 
     &             (aM_stu**(1.0-a_IMF) - aM_SNII**(1.0-a_IMF)) 
        endif
#ifdef ENRICH
c
c....   initialize fraction of stellar mass returned as metals
c
        if ( a_IMF .eq. 2.0 ) then 
          dummy1 = 1.d0 / log(aM_stu/aM_stl)
        else
          dummy1 = (2.d0 - a_IMF) / 
     &             (aM_stu**(2.0-a_IMF) - aM_stl**(2.0-a_IMF)) 
        endif
        fp(1) = a_IMF
        dumst = 1.d-1 * (aM_stu - aM_SNII)
        dummy2 = 
     &    INTEGRATE(fej_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        fmass_met = dummy2 * dummy1 
#endif
        fp(1) = a_IMF
        dummy2 = 
     &    INTEGRATE(fm_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(f_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)

        aMSN_ave = dummy2 / dummy1
        aN_SNII = dummy1 * An_IMF
        
        dumst = 1.d-1 * (aM_stu - aM_stl)
        dummy2 = 
     &    INTEGRATE(fm_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(fm_IMF_plaw, fp, 1, aM_stl, aM_stu, dumst, 1.d-9)
        fmass_SN = dummy2 / dummy1 

#ifdef FEEDBACK_SNIa         
        dummy1 = 
     &   INTEGRATE(fm_IMF_plaw, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9)
        dummy2 = 
     &   INTEGRATE(f_IMF_plaw, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9)
        dummy3 = 
     &   INTEGRATE (f_SNIa_func, fp, 2, 1.d-2, 1.d3 , 1.d-1 , 1.d-9)
c
c....   to get SNIa rate in SN/yr: R_Ia = Mstar * 1.d-9 * RIaf * f_SNIa
c
        RIaf = 1.d-9 * t0 * aM0 * C_SNIa * dummy2 / dummy1 / dummy3 /
     &       t_SNIa
#endif
        
      ELSE ! Miller-Scalo IMF
        dumst = 1.d-1 * (aM_stu - aM_stl)
        dummy1 = 
     &    INTEGRATE ( fm_IMF, fp, 1, aM_stl, aM_stu, dumst, 1.d-9 )
        An_IMF = 1.d0 / dummy1

        dumst  = 1.d-1 * (aM_stu - aM_SNII)
        dummy2 = 
     &    INTEGRATE ( f_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9 )
        if ( dummy1 .le. 0.0 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init : something is wrong dummy1 <= 0 :', 
     &     dummy1
          close ( 13 ) 
          stop          
        endif
c
        Aprime = dummy2 / dummy1
c
        if ( Aprime .le. 0.0 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*)' 1  : SF_Init : something is wrong Aprime <= 0 :',
     &      Aprime
          close ( 13 ) 
          stop          
        endif 
#ifdef ENRICH
c
c....   initialize fraction of stellar mass returned as metals
c
        
        dumst = 1.d-1 * (aM_stu - aM_SNII)
        dummy2 = 
     &    INTEGRATE (fej_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)        
        fmass_met = dummy2 / dummy1         
#endif
        dummy2 = 
     &    INTEGRATE(fm_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(f_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)

        aMSN_ave = dummy2 / dummy1
        aN_SNII = dummy1 * An_IMF
 
        dumst = 1.d-1 * (aM_stu - aM_stl)
        dummy2 = 
     &    INTEGRATE(fm_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(fm_IMF, fp, 1, aM_stl, aM_stu, dumst, 1.d-9)
        fmass_SN = dummy2 / dummy1 
#ifdef FEEDBACK_SNIa
        dumst = 1.d-1 * ( aM_SNIa2 - aM_SNIa1 )
        dummy1 = 
     &    INTEGRATE ( fm_IMF, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9 )
        dummy2 = 
     &    INTEGRATE ( f_IMF, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9 )
        dummy3 = 
     &    INTEGRATE ( f_SNIa_func, fp, 2, 1.d-3, 1.d3 , 1.d-1 , 1.d-9)
c
c....   to get SNIa rate R_Ia =  Mstar * RIaf * f_SNIa(dt)
c


        RIaf = 1.d-9 * t0 * aM0 * 
     &         C_SNIa * dummy2 / dummy1 / dummy3 / t_SNIa
#endif
      ENDIF

#ifdef ENRICH
      write(*,*) 'fmass_met =',fmass_met
#endif
      write(*,*) 'fraction of mass in SNII:', fmass_SN
      write(*,*) 'average M_SN=', aMSN_ave
      write(*,*) 'N of SNII:', aN_SNII
      write(*,*) 'N of SNIa:', C_SNIa  * dummy2/dummy1

      C_fb = 1.d51 * E_51 * Aprime * aM0 / E_0       
      C_fbIa = 1.d51 * E_51 / E_0       

      rho_SF_fact = rho_SF * 8.9d4 / hubble**2 / Om0 / (1.d0 - Y_p)
      write(*,*) 'C_fb, C_fbIa, E_51, Aprime, aM0, E_0, rho_SF =',
     &      C_fb, C_fbIa, E_51, Aprime, aM0, E_0, rho_SF_fact

      return
      end

c     ----------------------------------------
      real*8 function f_SNIa_func ( xd , fp , np )
c     ----------------------------------------
c
c     function for integrating S R_SNIa(t) dt
c
c     xd = t_SNIa / t 
c
      integer np 
      real*8 xd, fp(np)
      
      f_SNIa_func = max ( dexp(-xd*xd) / dsqrt(xd) , 1.d-20 )

      return
      end

c     ----------------------------------------
      real*8 function f_SNIa ( xd , fp , np )
c     ----------------------------------------
c
c     function describing phenomenological shape of the 
c     R_SNIa(t) function
c
c     x = t_SNIa / t 
c
      integer np 
      real*8 xd, fp(np)
      
      f_SNIa = max ( dexp(-xd*xd) * xd**1.5 , 1.d-20 )

      return
      end

c     -----------------------------------------------
      real*8 function f_IMF_plaw ( amstar , fp , np )
c     ----------------------------------------------- 
c
c     unnormalized stellar IMF approximated by a single power law
c  
c     input: real*8 amstar - star mass in M_sun 
c            fp(1) - power index (e.g., 2.35 for the Salpeter IMF)

      integer np 
      real*8 amstar, fp(np)
c
      f_IMF_plaw = amstar**(-fp(1))
      
      return
      end

c     ------------------------------------------------
      real*8 function fm_IMF_plaw ( amstar , fp , np )
c     ------------------------------------------------ 
c
c     unnormalized stellar IMF approximated by a single power law
c         multiplied by mass
c  
c     input: real*8 amstar - star mass in M_sun 
c            fp(1) - power index (e.g., 2.35 for the Salpeter IMF)

      integer np 
      real*8 amstar, fp(np)
c
      fm_IMF_plaw = amstar**(1.d0 - fp(1))
      
      return
      end

c     ------------------------------------------------
      real*8 function fej_IMF_plaw ( amstar , fp , np )
c     ------------------------------------------------ 
c
c     unnormalized stellar IMF approximated by a single power law
c         multiplied by mass and by fraction of mass ejected as metals 
c  
c     input: real*8 amstar - star mass in M_sun 
c            fp(1) - power index (e.g., 2.35 for the Salpeter IMF)

      integer np 
      real*8 amstar, fp(np)
      real*8 fmet_ej
c
      fej_IMF_plaw = amstar**(1.d0 - fp(1)) * fmet_ej ( amstar )
      
c      write(*,*) amstar, fmet_ej ( amstar )

      return
      end

c     ------------------------------------------
      real*8 function f_IMF ( amstar , fp , np )
c     ------------------------------------------ 
c
c     unnormalized stellar IMF
c  
c     input: real*8 amstar - star mass in M_sun 
c            fp(np) - dummy array, is not used here

      integer np 
      real*8 amstar, fp(np)
      real*8 C_1, C_2
      parameter ( C_1 = 1.09 , C_2 = -1.02 )   ! parameters of ther M-S IMF
      real*8 am0_Ch, beta_Ch, alpha_Ch
      parameter ( am0_Ch = 716.4 , beta_Ch = 0.25, alpha_Ch = -3.3 )

c========= STG: switched from M-S to Chabrier:
c
c.... Miller-Scalo (1979, ApJS 41, 513, eq. 30, Table 7) IMF
c
c      f_IMF = dexp ( -C_1 * (dlog10(amstar) - C_2)**2 ) / amstar
c
c.... Chabrier, G. (2001, ApJ 554, 1274)
c
      f_IMF = dexp(-(am0_Ch/amstar)**beta_Ch) * amstar**alpha_Ch
c===========

      return
      end

c     ------------------------------------------
      real*8 function fm_IMF ( amstar , fp , np )
c     ------------------------------------------ 
c
c     unnormalized stellar IMF * m_star
c  
c     input: real*8 amstar - star mass in M_sun 
c            fp(np) - dummy array, is not used here

      integer np 
      real*8 amstar , fp(np), f_IMF

      fm_IMF = amstar * f_IMF ( amstar , fp , np )
      
      return
      end

c     ------------------------------------------
      real*8 function fej_IMF ( amstar , fp , np )
c     ------------------------------------------ 
c
c     unnormalized stellar IMF * m_star * fraction of mass ejected as metals 
c  
c     input: real*8 amstar - star mass in M_sun 
c            fp(np) - dummy array, is not used here

      integer np 
      real*8 amstar , fp(np), f_IMF, fmet_ej

      fej_IMF = amstar * f_IMF ( amstar , fp , np ) * fmet_ej ( amstar )
      
      return
      end

c     ----------------------------------
      real*8 function fmet_ej ( amstar )
c     ---------------------------------- 
c
c     for amstar >= M_SN, returns fraction of mass ejected as metals into ISM 
c
c     current approximation is based on the results of Woosley & Weaver 1995
c     and ignores metallicity dependence of the yields 
c
c     input : real*8 amstar in Msun
c
      real*8 amstar 
      
      fmet_ej = min ( 2.0d-1, max(1.0d-2*amstar - 6.d-2,1.d-20) )

      return
      end
