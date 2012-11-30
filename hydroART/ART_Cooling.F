
c     ---------------------------------------------------------
      subroutine Cooling_Rate ( rhogl, T_g , Z_met , coolrate )
c     ---------------------------------------------------------
c     input:  
c         rhogl - log10(n_H), where n_H is hydrogen number density in cm^-3
c         T_g  - gas temperature in units of 10^4 K
c         Z_met - log10([rho_metals/rho_gas/solar)
c         if ENRICH is not defined Z_met will be ignored 
c         if SD93_COOLING is defined, rhogl will be ignored 
c
c     output: coolrate - net output cooling/heating rate in code units
c        de/dt = e - coolrate * rho_g**2 / aexpn
c        where e is internal energy per unit volume and rho_g is 
c        gas density (e and rho_g are in code units)
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      real*8 rhogl, rho_g, T_g, Z_met, coolrate

#ifdef CLOUDY_COOLING
      include 'a_clcool.h' 
c
c.... use net cooling/heating rate tabulated using Cloudy
c
      real*8 Tlog
      real*8 ac, bc, ah, bh 

      Tlog = log10(T_g) + 4.d0 
      it1 = int((Tlog - tlmin)*dlti) + 1
      it2 = it1 + 1
      id1 = int((rhogl - dlmin)*dldi) + 1
      id2 = id1 + 1
#ifndef METALCOOLING 
      iz1 = 1
      iz2 = 1
#endif
#ifdef METALCOOLING
      iz1 = int((Z_met - zlmin)*dlzi) + 1
      iz2 = iz1 + 1
#endif
      it1 = max(it1,1) 
      it1 = min(it1,nlt)
      it2 = max(it2,1) 
      it2 = min(it2,nlt)
      id1 = max(id1,1) 
      id1 = min(id1,nld)
      id2 = max(id2,1) 
      id2 = min(id2,nld)
      iz1 = max(iz1,1) 
      iz1 = min(iz1,nlz)
      iz2 = max(iz2,1) 
      iz2 = min(iz2,nlz)

      td = tlmin + dlt * (it1 - 1) + dlt
      dd = dlmin + dld * (id1 - 1) + dld
      zd = zlmin + dlz * (iz1 - 1) + dlz
      t1 = (td - Tlog) * dlti
      d1 = 1.d0 - t1
      t2 = (dd - rhogl) * dldi
      d2 = 1.d0 - t2
      t3 = (zd - Z_met) * dlzi 
      d3 = 1.d0 - t3
c
c.... trilinear (CIC) interpolation 
c
CEVERINO04302008 
      IF ( rhogl.ge.-1..AND.ioptUV.eq.1) THEN ! ! limit UVBackG to z=8-Flux for high densities
         coolrate = t1*t2*t3 * coolcl(it1,id1,iz1,9) +
     &        d1*t2*t3 * coolcl(it2,id1,iz1,9) + 
     &        t1*d2*t3 * coolcl(it1,id2,iz1,9) + 
     &        d1*d2*t3 * coolcl(it2,id2,iz1,9) +
     &        t1*t2*d3 * coolcl(it1,id1,iz2,9) + 
     &        d1*t2*d3 * coolcl(it2,id1,iz2,9) + 
     &        t1*d2*d3 * coolcl(it1,id2,iz2,9) + 
     &        d1*d2*d3 * coolcl(it2,id2,iz2,9)   
      ELSE                      ! normal cooling
         coolrate = t1*t2*t3 * ccl_rs(it1,id1,iz1) +
     &        d1*t2*t3 * ccl_rs(it2,id1,iz1) + 
     &        t1*d2*t3 * ccl_rs(it1,id2,iz1) + 
     &        d1*d2*t3 * ccl_rs(it2,id2,iz1) +
     &        t1*t2*d3 * ccl_rs(it1,id1,iz2) + 
     &        d1*t2*d3 * ccl_rs(it2,id1,iz2) + 
     &        t1*d2*d3 * ccl_rs(it1,id2,iz2) + 
     &        d1*d2*d3 * ccl_rs(it2,id2,iz2)
      ENDIF


CEVERINO10022006
C      ELSE ! only cools, no heating.
c$$$      coolrate = t1*t2*t3 * ccl_rsC(it1,id1,iz1) +
c$$$     &           d1*t2*t3 * ccl_rsC(it2,id1,iz1) + 
c$$$     &           t1*d2*t3 * ccl_rsC(it1,id2,iz1) + 
c$$$     &           d1*d2*t3 * ccl_rsC(it2,id2,iz1) +
c$$$     &           t1*t2*d3 * ccl_rsC(it1,id1,iz2) + 
c$$$     &           d1*t2*d3 * ccl_rsC(it2,id1,iz2) + 
c$$$     &           t1*d2*d3 * ccl_rsC(it1,id2,iz2) + 
c$$$     &           d1*d2*d3 * ccl_rsC(it2,id2,iz2)         
C      ENDIF
      return
#endif

#ifdef SD93_COOLING
c
c.... use CIE cooling curves of Sutherland & Dopita 1993
c
      real*8 dummy, Z1
      real*8 csmall, Zsolar, Zsolari, Zpri, slope, zerop
      parameter ( ntab = 5000 , nmetal = 8 , csmall=1.d-10 )
      parameter ( Zsolar = 0.0199 , Zsolari = 1.d0/Zsolar )
      parameter ( Zpri = 1.d-5 * Zsolar )
      common / CLRATE / coolfun(0:ntab,nmetal), zval(nmetal)

      itab = int(100.0 * sqrt(sqrt(max(T_g,1.d-4))))
      itab = max( itab , 0 )
      itab = min( itab , ntab )
#ifndef METALCOOLING   
      coolrate = coolfun(itab,1) ! primordial composition 
#endif
#ifdef METALCOOLING 
c
c.... in the following, the small difference between primordial and Z/Zsun=-3
c     will be ignored 
c
      if( Z_met .gt. Zval(2) ) then
         Z1 = Z_met
      else
        coolrate = coolfun(itab,1) ! primordial composition
        return
      endif
c
c...  if Z > Fe/H=0.5, extrapolate linearly in log10(Z/Zsun)
c
      if ( Z1 .gt. zval(nmetal) ) then
        slope =  (coolfun(itab,nmetal) - coolfun(itab,nmetal-1)) /
     &           (zval(nmetal) - zval(nmetal-1))
        zerop = coolfun(itab,nmetal) - slope * zval(nmetal)
        coolrate = slope * Z1 + zerop
        return
      endif
c
c.... interpolate linearly in log(Z/Zsun)
c
      do j = 3 , nmetal
        if ( Z1 .le. zval(j) ) then
          slope = (coolfun(itab,j)-coolfun(itab,j-1)) / 
     &            (zval(j)-zval(j-1))
          zerop = coolfun(itab,j) - slope*zval(j)
          coolrate = slope * Z1 + zerop
          return
        endif
      enddo

      call Open_ASCII_File ( iOErr , errorfile , sequent ) 
      write(iOErr,*) '1 : error : in Cooling_Rate: something is wrong:'
      write(iOErr,*) 'T_g, Z_met =', T_g, Z_met
      close (iOErr) 
      stop
#endif
#endif

      return
      end

c     -------------------------
      subroutine Set_Cooling ()
c     -------------------------
c
c     tabulate cooling curve as a function of T in units of 10^4 K
c     table entry i corresponds to i = 100*T^1/4
c     AL_0 must be set in Set_Units prior to call to this routine
c      
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

#ifdef CLOUDY_COOLING
      include 'a_clcool.h' 
      real*8 cdum, hdum, ct_crit 

#ifdef SD93_COOLING
c.... check that other cooling tables are not defined 
      call Open_ASCII_File ( iOErr , errorfile , sequent ) 
      write(iOErr,*) '1 : error : in Set_Cooling:'
      write(iOErr,*)
     &    ' both CLOUDY_COOLING and SD93_COOLING are defined in a_def.h'
      write(iOErr,*)
     &     ' while only one must be defined. stopping... '
      close (iOErr) 
      stop
#endif
c
c.... use net cooling/heating rate tabulated using Cloudy
c
      open ( 40 , file = 'clcool.dat' )
      read(40,*) tlmin, tlmax, dlt, nlt
      read(40,*) dlmin, dlmax, dld, nld
      read(40,*) Zlmin, Zlmax, dlZ, nlz
      read(40,*) rsmin, rsmax, drs, nrs

      write(*,*) tlmin, tlmax, dlt, nlt
      write(*,*) dlmin, dlmax, dld, nld
      write(*,*) Zlmin, Zlmax, dlZ, nlz
      write(*,*) rsmin, rsmax, drs, nrs

      dlti = 1.d0 / dlt
      dldi = 1.d0 / dld
      dlzi = 1.d0 / dlz
      drsi = 1.d0 / drs 

      do irs = 1, nrs
        do ilz = 1 , nlz
          do ild = 1 , nld 

            ct_crit = 0.0
            do ilt = 1 , nlt
              read(40,*) d1, d2, d3, d4, d5, d6, d7, d8, d9, cdum, hdum
              if ( d1 .ge. 3.2 .and. ct_crit .eq. 0.0 ) then
                ct_crit = cdum  
              endif
              cdum = max ( cdum , smallrate )
              cdum = max ( cdum , ct_crit   )  ! fix for a trough in equilibrium H2 cooling curve
              hdum = max ( hdum , smallrate ) 
              coolcl(ilt,ild,ilz,irs) = (cdum - hdum) * 1.d23 * AL_SD
              coolclC(ilt,ild,ilz,irs) = (cdum) * 1.d23 * AL_SD
              coolclH(ilt,ild,ilz,irs) = (hdum) * 1.d23 * AL_SD              
              f_ion(ilt,ild,ilz,irs) = d6 / 10.d0**(d2)
 100          format(6(e12.6,1x))
            enddo
          enddo
        enddo
      enddo
      close ( 40 ) 

      return
#endif

#ifdef SD93_COOLING
c
c.... use CIE cooling curves of Sutherland & Dopita 1993
c
      real csmall, Zsolar
      parameter ( ntab = 5000, nmetal = 8  ) 
      parameter ( csmall = 1.0e-10, Zsolar = 0.0199 )
      common / CLRATE / coolfun(0:ntab,nmetal), zval(nmetal)
      real xd, temp(200), coolr(200,nmetal), ap(nmetal)
      character*79 Line
      real CINTR
      real*8 T4, cr_H2, Z_met

      open (40 , file = 'coolfunc.dat' )
      read (40,40) Line
 40   format(A)
      write (*,*) 'setting cooling rates...'
      write (*,40) Line
      read (40,*) nmet, (zval(j),j=1,nmet)
      if (nmet .ne. nmetal) then 
        call Open_ASCII_File ( iOErr , errorfile , sequent ) 
        write(iOErr,*) '1 : *error in Set_Cooling: nmet .ne. nmetal'
        write(iOErr,*) 'nmet =',nmet,' nmetal =',nmetal
        write(iOErr,*) 'stopping...'
        close (iOErr)
        stop
      endif

      n = 0
 20   read(40,*,end=30,err=30) aT,(ap(j),j=1,nmetal)
      n = n + 1
      temp(n) =10.d0**(aT-4.d0)
      do j = 1 , nmetal
        coolr(n,j) = (10.d0**(ap(j)+23.0) * AL_SD)
      enddo
      go to 20

 30   write (*,*) ' number of lines read for cooling=', n
      do j = 1 , nmetal 
        Z_met = 10.d0**(zval(j))
        do i = 0 , ntab
          xd = (float(i)/100.0)**4
          coolfun(i,j)  = CINTR ( xd, temp, coolr(1,j), n )
c
c....     H2 cooling at low T
c          
          T4 = xd 
          call CoolingRate_H2 ( T4 , Z_met , cr_H2 )
          cr_H2 = cr_H2 * 1.d23 * AL_SD
          if ( cr_H2 .gt. coolfun(i,j) .and. T4 .lt. 4.25 ) then 
            coolfun(i,j) = cr_H2
          endif
        enddo 
      enddo
c
c.... small ad hoc cooling rate below 10^4K 
c
c      do j = 1 , nmetal
c        do i = 0 , 99
c          dee = min( ((float(i)-100.) * 0.2)**2, 30. )
c          coolfun(i,j)=coolfun(100,j)/exp(dee) +
c     &                 5.0e-10 * coolfun(100,j)
c        enddo
c      enddo
      
      close (40)
      return
#endif
      
      call Open_ASCII_File ( iOErr , errorfile , sequent ) 
      write(iOErr,*) '1 : error : in Set_Cooling:'
      write(iOErr,*)
     &     ' routine is called but cool. rates are not defined'
      write(iOErr,*)
     &     ' CLOUDY_COOLING or SD93_COOLING must be defined in a_def.h'
      close (iOErr) 
      stop
      
      end
     
c     -------------------------------------------
      real function CINTR ( xd , temp , cool , n )
c     -------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
      integer n
      real xd, temp(n), cool(n)

      if ( xd .le. temp(1) ) then
        CINTR =xd/temp(1)*cool(1)
        return 
      endif
      if ( xd .ge. temp(n) ) then
        CINTR = sqrt(xd/temp(n)) * cool(n)
        return
      endif
      do i = 2 , n
        if ( xd .le. temp(i) ) then
          CINTR = (cool(i)-cool(i-1))/(temp(i)-temp(i-1))*(xd-temp(i-1))
     &             + cool(i-1) 
          return
        endif
      enddo

      call Open_ASCII_File ( iOErr , errorfile , sequent )
      write (iOErr,*) '1  * error in CINTR: xd=', xd, temp(1), temp(n)
      write (iOErr,*) 'stopping ...'
      close (iOErr)
      stop

      return
      end
c

c     ---------------------------------------------------
      subroutine CoolingRate_H2 ( T_g , Z_met, coolrate )
c     ---------------------------------------------------
c
c     approximation to cooling rate due to H2 
c     in the low-density regime 
c     [Galli & Palla 1998, A&A 335, 403 [eq. A7]]
c     in this regime (n_H <~ 10^2 cm^-3) cooling rate 
c     is independent of density. 
c        
c     T_g is gas temperature in units of 10^4 K      
c     Z_met - metallicity in solar units = [Fe/H] / [Fe/H]_solar
c     coolrate is in erg sm^3 s^-1
c
c     the rate assumes that for zero metallicity H2 reaches universal 
c     abundance of n_H2 / N_H = 1.e-3 to 1.e-1  [see Haiman and Oh 2002]
c     the assumed value is set in parameter H2_u below 
c
c     for non-zero metallicities the cooling rate is assumed to be boosted
c     with an hoc additive term which approximately mimicks 
c     results of Spaans & Norman 1997, ApJ 483, 87 
c     (see also Wada and Norman 2001, ApJ 547, 172)
c     
c     to get the cooling rate multiply coolrate by n_H^2
c
      real*8 T_g, Z_met, coolrate
      real*8 lt, lt2, tm
      real*8 H2_u  ! assumed universal fraction of H2
      parameter ( H2_u = 1.d-2 )
      real*8 sigt
      parameter ( sigt = 0.75 )

      tm = max(T_g, 1.3d-3) * 1.d4    ! no cooling below 13 Kelvin ...
      tm = min(tm, 1.d5)        ! fixes numerics
      lt = dlog10(tm)
      lt2 = lt * lt
c     low density limit from Galli and Palla (1998) 
c     multiplied by 1d23 
      coolrate = 
     &           H2_u * 
     &           10.d0**(-103.d0 + 97.59d0 * lt     - 
     &                           (48.05d0       -
     &                            10.8d0   * lt + 
     &                            0.9032d0 * lt2  ) * lt2 ) 
c
c.... add a boost due to metals to mimick Spaans cooling curve
c     [see Fig.1 in Wada & Norman 2001, ApJ 547, 172
c
      if ( T_g .lt. 5.d0 ) then 
        coolrate = coolrate + 
     &      2.5d-25* dexp(-5.d-1*(lt - 4.)**2/sigt**2) / sigt * 
     &      (1.d0 + 10.d0 * Z_met)
      endif

      return
      end      

c     ------------------------------------
      subroutine Test_Cooling ( filename )
c     ------------------------------------
      real*8 Tdum, cr1, cr2, cr3, Zdum, rhogl
      character*256 filename 

      lfn = index ( filename , ' ' ) - 1 
      open ( 20 , file = filename(1:lfn) ) 
      Tlmin = 1.d0 
      Tlmax = 9.d0
      dTl = 0.05 
      nt = int((Tlmax - Tlmin)/dTl) + 1
      do it = 1 , nt  
        Tdum = 10.d0**(Tlmin + float(it-1)*dTl)
        Zdum = 1.d-10
        call Cooling_Rate ( rhogl, Tdum*1.d-4, Zdum, cr1 )
        Zdum = 2.d-3
        call Cooling_Rate ( rhogl, Tdum*1.d-4, Zdum, cr2 )
        Zdum = 2.d-2
        call Cooling_Rate ( rhogl, Tdum*1.d-4, Zdum, cr3 )
        write(20,90) Tdum, cr1, cr2, cr3 
      enddo
 90   format (4(e11.6,1x))
      close ( 20 ) 

      return
      end

c     --------------------------------------------------
      subroutine X_ion ( rhogl , T_g , Z_met , rs , xion )
c     --------------------------------------------------
c
c     equilibrium ionized fraction for a given condition
c     calculated using Cloudy tables
c
c     input:  
c         rhogl - log10(n_H), where n_H is hydrogen number density in cm^-3
c         T_g  - gas temperature in K
c         Z_met - log10([rho_metals/rho_gas/solar)
c         if ENRICH is not defined Z_met will be ignored 
c
c     output: xion - ionized fraction = n_HII/n_H
c
#     include "a_def.h"
c
      include 'a_clcool.h' 
      include 'a_tree.h'
      include 'a_control.h'

      real*8 rs, rhogl, T_g, Z_met, xion 
      real*8 ac, bc, xion1, xion2

      Tlog = log10(T_g)  
      it1 = int((Tlog - tlmin)*dlti) + 1
      it2 = it1 + 1
      id1 = int((rhogl - dlmin)*dldi) + 1
      id2 = id1 + 1
#ifndef METALCOOLING 
      iz1 = 1
      iz2 = 1
#endif
#ifdef METALCOOLING
      iz1 = int((Z_met - zlmin)*dlzi) + 1
      iz2 = iz1 + 1
#endif
      it1 = max(it1,1) 
      it1 = min(it1,nlt)
      it2 = max(it2,1) 
      it2 = min(it2,nlt)
      id1 = max(id1,1) 
      id1 = min(id1,nld)
      id2 = max(id2,1) 
      id2 = min(id2,nld)
      iz1 = max(iz1,1) 
      iz1 = min(iz1,nlz)
      iz2 = max(iz2,1) 
      iz2 = min(iz2,nlz)

      td = tlmin + dlt * (it1 - 1) + dlt
      dd = dlmin + dld * (id1 - 1) + dld
      zd = zlmin + dlz * (iz1 - 1) + dlz
      t1 = (td - Tlog) * dlti
      d1 = 1.d0 - t1
      t2 = (dd - rhogl) * dldi
      d2 = 1.d0 - t2
      t3 = (zd - Z_met) * dlzi 
      d3 = 1.d0 - t3

      irs = int((rs - rsmin)*drsi) + 1
      irs1 = max(irs,1)
      irs1 = min(irs1,nrs)
      irs2 = min(irs+1,nrs)
      irs2 = max(irs+1,1)

      xion1 = t1*t2*t3 * f_ion(it1,id1,iz1,irs1) +
     &        d1*t2*t3 * f_ion(it2,id1,iz1,irs1) + 
     &        t1*d2*t3 * f_ion(it1,id2,iz1,irs1) + 
     &        d1*d2*t3 * f_ion(it2,id2,iz1,irs1) +
     &        t1*t2*d3 * f_ion(it1,id1,iz2,irs1) + 
     &        d1*t2*d3 * f_ion(it2,id1,iz2,irs1) + 
     &        t1*d2*d3 * f_ion(it1,id2,iz2,irs1) + 
     &        d1*d2*d3 * f_ion(it2,id2,iz2,irs1)

      xion2 = t1*t2*t3 * f_ion(it1,id1,iz1,irs2) +
     &        d1*t2*t3 * f_ion(it2,id1,iz1,irs2) + 
     &        t1*d2*t3 * f_ion(it1,id2,iz1,irs2) + 
     &        d1*d2*t3 * f_ion(it2,id2,iz1,irs2) +
     &        t1*t2*d3 * f_ion(it1,id1,iz2,irs2) + 
     &        d1*t2*d3 * f_ion(it2,id1,iz2,irs2) + 
     &        t1*d2*d3 * f_ion(it1,id2,iz2,irs2) + 
     &        d1*d2*d3 * f_ion(it2,id2,iz2,irs2)


      if ( irs1 .eq. irs2 ) then 
         xion = xion1
      else
        rs1 = rsmin + drs*(irs1-1)
        rs2 = rsmin + drs*(irs2-1)
        ac = (xion2 - xion1) /
     &         (rs2 - rs1)
        bc = xion1 - ac * rs1 
        xion = ac * rs + bc 
      endif

      return
      end

c     -------------------------------------
      subroutine Set_Cooling_Rate_rs ( rs ) 
c     -------------------------------------
c
c     prepare cooling rate table for a given redshift rs
c     from the Cloudy cooling rate table 
c
c     real*8 rs - redshift 
c     
c     this routine is only used if CLOUDY_COOLING is defined in a_def.h
c     it should be called only in the beginning of every step 
c
c
      include 'a_clcool.h' 
      real*8 rs, rsd
      real*8 ac, bc, ah, bh

      irs = int((rs - rsmin)*drsi) + 1
      irs1 = max(irs,1)
      irs1 = min(irs1,nrs)
      irs2 = min(irs+1,nrs)

      if ( irs1 .eq. irs2 ) then 
        do ilt = 1 , nlt
        do ild = 1 , nld
        do ilz = 1 , nlz 
          ccl_rs(ilt,ild,ilz) = coolcl(ilt,ild,ilz,irs1)
          ccl_rsC(ilt,ild,ilz) = coolclC(ilt,ild,ilz,irs1)
          ccl_rsH(ilt,ild,ilz) = coolclH(ilt,ild,ilz,irs1)
        enddo
        enddo
        enddo
      else
        rs1 = rsmin + drs*(irs1-1)
        rs2 = rsmin + drs*(irs2-1)

        do ilt = 1 , nlt
        do ild = 1 , nld
        do ilz = 1 , nlz 
          ac = (coolcl(ilt,ild,ilz,irs2) - coolcl(ilt,ild,ilz,irs1)) /
     &         (rs2 - rs1)
          bc = coolcl(ilt,ild,ilz,irs1) - ac * rs1 
          ccl_rs(ilt,ild,ilz) = ac * rs + bc 
c
          ac = (coolclC(ilt,ild,ilz,irs2) - coolclC(ilt,ild,ilz,irs1)) /
     &         (rs2 - rs1)
          bc = coolclC(ilt,ild,ilz,irs1) - ac * rs1 
          ccl_rsC(ilt,ild,ilz) = ac * rs + bc
c
          ac = (coolclH(ilt,ild,ilz,irs2) - coolclH(ilt,ild,ilz,irs1)) /
     &         (rs2 - rs1)
          bc = coolclH(ilt,ild,ilz,irs1) - ac * rs1 
          ccl_rsH(ilt,ild,ilz) = ac * rs + bc  
        enddo
        enddo
        enddo
      endif

      return
      end

