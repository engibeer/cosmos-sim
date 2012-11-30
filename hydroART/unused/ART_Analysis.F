c...  Analysis routines
c       Compute_Cold_Fraction
c       SFR_Rvir
c       SFR_r1_r2
c       AverageSF
c       Compute_Phase 
c       Output_Average_Profiles
c       Output_Average_Profiles_CSF
c       Output_Average_Profiles_Stars
c       Output_Fit_Profiles
c   
c     ----------------------------------------------------------
      subroutine Compute_Cold_Fraction ( xc, yc, zc, rvir, Tmax,
     &                                   totmass, coldmass )
c     ----------------------------------------------------------
c     compute fraction of gas with T < Tmax for an object centered
c     at (xc,yc,zc) and radius rvir
c     all quantities are in code units
c
c     input : xc, yc, zc, rvir, Tmax
c     output: mass of gas at T<Tmax, total mass of gas within rvir (masses are in code units)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      totmass = 0.
      coldmass = 0.

      call Get_MaxLevelNow ()

        DO Level = MinLevel , MaxLevelNow
          CellVolume = 1.0 * 2.0**(-3.0*Level)
          IF ( Level .eq. MinLevel ) THEN
          ELSE
            call Select_Cells ( Level , nLevel )
            do ic1 = 1 , nLevel
              icell = iSelect(ic1)
              do ic2 = 0 , 7
                ic = icell + ic2
                if ( iOctCh(ic) .eq. nil ) then
                  call Ps ( ic , Posx , Posy , Posz )
                  call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz)
                  rr = sqrt( (Posx-xc)**2 +
     &                       (Posy-yc)**2 +
     &                       (Posz-zc)**2 )
                  if ( rr .lt. rvir ) then
                    cellmass = hvar(1,ic) * CellVolume
                    totmass = totmass + cellmass
                    Tcell = hvar(6,ic) / hvar(1,ic)
                    if ( Tcell .lt. Tmax ) then
                      coldmass = coldmass + cellmass
c                      write(*,*) ic, cellmass, coldmass, totmass
                    endif
                  endif
                endif
              enddo
            enddo
          ENDIF
        ENDDO

      if ( totmass .eq. 0. ) then
        write(*,*) 'error: in Compute_Cold_Fraction: totmass = 0', totmass
        return
      endif
c
      return
      end
c
c     ------------------------------------------------------
      subroutine SFR_rvir ( xc, yc, zc, rvir ,
     &                      amin , amax , nabin , ibintype ,  
     &                      deltavir, fname_sfr ) 
c     ------------------------------------------------------
c
c     Computing SFR for stellar particles within < Rvir
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc , amin , amax , deltavir
      integer nabin, ibintype
      character*256 fname_sfr
c
      double precision xx, yy, zz 
      real sfr_bin(nabin),msi_bin(nabin),aexpn_bin(nabin)
      real t_bin(nabin),dt_bin(nabin)
      integer icnts_bin(nabin)

c
c...  prepare an array for SFR(z)
c

c      write(*,*) xc,yc,zc,rvir,amin,amax,nabin,ibintype,deltavir
c      stop

      ! Binning uniformly in the expansion factor
      da = (amax-amin) / nabin 
      tnow = age ( t )  ! age of the universe now
      
      do ibin = 1 , nabin
        sfr_bin(ibin) = 0. 
        msi_bin(ibin) = 0.
        icnts_bin(ibin) = 0
        ! expansion factor
        a1 = amin + da * (ibin-1) 
        a2 = amin + da * ibin 
        ! convert expansion factor to code time
        t1 = a2b ( a1*1.d0 )
        t2 = a2b ( a2*1.d0 )
        ! convert code time to time in Gyr
        t1 = age ( t1*1.d0 )
        t2 = age ( t2*1.d0 )
        aexpn_bin(ibin)  = ( a2 + a1 ) / 2.0
        dt_bin(ibin) = t2 - t1  ! time interval in Gyr
        t_bin(ibin) = tnow - ( t2 + t1 ) / 2.0  ! Look back time in Gyr
      enddo
      
c
c...  loop over stellar particles (only)
c
      write(*,*) 'analyzing stellar particles...' 
      nstars = nsp(nspecies,2)-nsp(nspecies,1)+1
      istar1 = nsp(nspecies,1) 
      
      do is = nsp(nspecies,2), nsp(nspecies,1), -1 
        xx = x(is)
        yy = y(is)
        zz = z(is)
c....   correct coordinate to insure periodic bundaries        
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz)
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
 
c...    count all the stellar particles
        if ( rp .lt. rvir ) then
          i = is - istar1 + 1
          ! this is when the star formed [Gyr]
c          tsf = age ( tbirth(i)*1.d0 )  
          ! age of the star [Gyr]
c          ts_age = tnow - age ( tbirth(i)*1.d0 )
c          write(*,*) is, tsf, ts_age 
          ! convert code time to expansion factor
          asf = b2a ( dble( tbirth(i) ) ) 
c          write(*,*) asf,tbirth(i)
          ibin = int( (asf - amin) /da ) + 1 
          if ( ibin .gt. 0 ) then 
            sfr_bin(ibin) = sfr_bin(ibin) + pw0(i)
            icnts_bin(ibin) = icnts_bin(ibin) + 1
          endif
        endif
      enddo

c
c.... computing SFR(<rvir) [Msun yr^-1]
c      
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2
      vol = 4.18879 * (rvir*r0)**3 / hubble**3 ! [Mpc^3]

      msi_bin(1) = sfr_bin(1)
      do ibin = 2 , nabin
        msi_bin(ibin) = sfr_bin(ibin) + msi_bin(ibin-1) 
      enddo

      do ibin = 1 , nabin
        sfr_bin(ibin) = sfr_bin(ibin) * aM0  ! [Msun]
        msi_bin(ibin) = msi_bin(ibin) * aM0  ! [Msun]
        dt_bin(ibin) = dt_bin(ibin) * 1.e+9  ! [yr]
        sfr_bin(ibin) = sfr_bin(ibin)/dt_bin(ibin) ! [Msun yr^-1]
c        write(*,*) ibin, icnts_bin(ibin), aexpn_bin(ibin), 
c     &       dt_bin(ibin), sfr_bin(ibin)
      enddo
      
c
c.... output SFR(<rvir) vs. aexpn to a file
c

      write(*,*) 'Writing data to a file'
      nsfr = index( fname_sfr , ' ' ) - 1
      open ( 12 , file = fname_sfr(1:nsfr) )
      write(12,9)  (1.0/aexpn-1.0), tnow
      write(12,10) xc , yc , zc  
      write(12,11) amin , amax , nabin , ibintype 
      write(12,12) rvir*r0, aMvir, deltavir
      write(12,13)
      write(12,14)
      write(12,15)
 9    format('# SFR(<rvir) at z=',f5.3,' t=',f8.5,' [Gyrs]')
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# amin = ',g12.6,' amax =',g12.6,
     &        ' nabin =',i5,' ibintype =',i1)
 12   format ('# Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 13   format ('# SFR(<rvir) vs. redshift :')
 14   format('#   ic  #part  aexpn   t   dt   SFR     Mi* ')
 15   format('#                     Gyr  yr  Msun/yr  Msun ')
      do i = 1 , nabin
         write(12,17) i,icnts_bin(i),aexpn_bin(i),t_bin(i),
     &        dt_bin(i),sfr_bin(i),msi_bin(i)
      enddo
 17   format(i3,i8,5(2x,g14.7))
      close ( 12 )
      write(*,*) 'done. writing...'

      return
      end
c
c     ------------------------------------------------------
      subroutine SFR_r1_r2 ( xc, yc, zc, r1, r2 ,
     &                      amin , amax , nabin , ibintype ,  
     &                      deltavir, fname_sfr ) 
c     ------------------------------------------------------
c
c     Computing SFR for stellar particles within r1 < r < r2
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc , r1, r2, amin , amax , deltavir
      integer nabin, ibintype
      character*256 fname_sfr
c
      double precision xx, yy, zz 
      real sfr_bin(nabin),msi_bin(nabin),aexpn_bin(nabin)
      real t_bin(nabin),dt_bin(nabin)
      integer icnts_bin(nabin)

c
c...  prepare an array for SFR(z)
c

c      write(*,*) xc,yc,zc,rvir,amin,amax,nabin,ibintype,deltavir
c      stop

      ! Binning uniformly in the expansion factor
      da = (amax-amin) / nabin 
      tnow = age ( t )  ! age of the universe now
      
      do ibin = 1 , nabin
        sfr_bin(ibin) = 0. 
        msi_bin(ibin) = 0.
        icnts_bin(ibin) = 0
        ! expansion factor
        a1 = amin + da * (ibin-1) 
        a2 = amin + da * ibin 
        ! convert expansion factor to code time
        t1 = a2b ( a1*1.d0 )
        t2 = a2b ( a2*1.d0 )
        ! convert code time to time in Gyr
        t1 = age ( t1*1.d0 )
        t2 = age ( t2*1.d0 )
        aexpn_bin(ibin)  = ( a2 + a1 ) / 2.0
        dt_bin(ibin) = t2 - t1  ! time interval in Gyr
        t_bin(ibin) = tnow - ( t2 + t1 ) / 2.0  ! Look back time in Gyr
      enddo
      
c
c...  loop over stellar particles (only)
c
      write(*,*) 'analyzing stellar particles...' 
      nstars = nsp(nspecies,2)-nsp(nspecies,1)+1
      istar1 = nsp(nspecies,1) 
      
      do is = nsp(nspecies,2), nsp(nspecies,1), -1 
        xx = x(is)
        yy = y(is)
        zz = z(is)
c....   correct coordinate to insure periodic bundaries        
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz)
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
 
c...    count all the stellar particles
        if ( rp .ge. r1 .and. rp .lt. r2 ) then
          i = is - istar1 + 1
          ! this is when the star formed [Gyr]
c          tsf = age ( tbirth(i)*1.d0 )  
          ! age of the star [Gyr]
c          ts_age = tnow - age ( tbirth(i)*1.d0 )
c          write(*,*) is, tsf, ts_age 
          ! convert code time to expansion factor
          asf = b2a ( dble( tbirth(i) ) ) 
c          write(*,*) asf,tbirth(i)
          ibin = int( (asf - amin) /da ) + 1 
          if ( ibin .gt. 0 ) then 
            sfr_bin(ibin) = sfr_bin(ibin) + pw0(i)
            icnts_bin(ibin) = icnts_bin(ibin) + 1
          endif
        endif
      enddo

c
c.... computing SFR(<rvir) [Msun yr^-1]
c      
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2
      vol = 4.18879 * (rvir*r0)**3 / hubble**3 ! [Mpc^3]

      msi_bin(1) = sfr_bin(1)
      do ibin = 2 , nabin
        msi_bin(ibin) = sfr_bin(ibin) + msi_bin(ibin-1) 
      enddo

      do ibin = 1 , nabin
        sfr_bin(ibin) = sfr_bin(ibin) * aM0  ! [Msun]
        msi_bin(ibin) = msi_bin(ibin) * aM0  ! [Msun]
        dt_bin(ibin) = dt_bin(ibin) * 1.e+9  ! [yr]
        sfr_bin(ibin) = sfr_bin(ibin)/dt_bin(ibin) ! [Msun yr^-1]
c        write(*,*) ibin, icnts_bin(ibin), aexpn_bin(ibin), 
c     &       dt_bin(ibin), sfr_bin(ibin)
      enddo
      
c
c.... output SFR(<rvir) vs. aexpn to a file
c

      write(*,*) 'Writing data to a file'
      nsfr = index( fname_sfr , ' ' ) - 1
      open ( 12 , file = fname_sfr(1:nsfr) )
      write(12,9)  (1.0/aexpn-1.0), tnow
      write(12,10) xc , yc , zc  
      write(12,11) amin , amax , nabin , ibintype 
      write(12,12) r1*r0, r2*r0, aMvir, deltavir
      write(12,13)
      write(12,14)
      write(12,15)
 9    format('# SFR(<rvir) at z=',f5.3,' t=',f8.5,' [Gyrs]')
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# amin = ',g12.6,' amax =',g12.6,
     &        ' nabin =',i5,' ibintype =',i1)
 12   format ('# R1, R2 =',f7.4, f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 13   format ('# SFR(<rvir) vs. redshift :')
 14   format('#   ic  #part  aexpn   t   dt   SFR     Mi* ')
 15   format('#                     Gyr  yr  Msun/yr  Msun ')
      do i = 1 , nabin
         write(12,17) i,icnts_bin(i),aexpn_bin(i),t_bin(i),
     &        dt_bin(i),sfr_bin(i),msi_bin(i)
      enddo
 17   format(i3,i8,5(2x,g14.7))
      close ( 12 )
      write(*,*) 'done. writing...'

      return
      end
c
c     ------------------------------------------------------
      subroutine AverageSF ( level, tsfave, fname )
c     ------------------------------------------------------
c
c     level - level on which to study SF law
c      
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'
      character*256 fname
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      real*8 xx, yy, zz 
      real*8 age
#ifdef PRESSUREFLOOR
      real*8 f_J, f_eJ
      parameter ( f_J = 10.  ) ! pressure floor parameter
#endif

#ifdef PRESSUREFLOOR
      DO L = MinL_Jeans , MaxLevelNow
c
c.... artificial pressure floor a-la Machacek et al. 2
        if ( L .ge. MinL_Jeans ) then
          f_eJ = 0.47746 * f_J * aexpn * 
     &           CellSize(MaxLevel)**2 * (gamma - 1.0)
        else
          f_eJ = 0.0
        endif
        if ( L .eq. MinLevel ) then
          do i1 = 1 , ncell0
            E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
            E_th = f_eJ * hvar(1,i1) * hvar(1,i1) 
            hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
          enddo 
        else
          call Select_Cells ( L , nLevel )
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              i1 = icell + ic2
              if ( iOctCh(i1) .eq. nil ) then
                E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
                E_th = f_eJ * hvar(1,i1) * hvar(1,i1)
                hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
              endif
            enddo
          enddo
        endif
      ENDDO
#endif
#ifndef PRESSUREFLOOR
      f_eJ = 0.0
#endif

      tsfaveb = tsfave * 1.e-9 ! yrs -> Gyrs

      do i = 1 , nctot
        ref(i) = 0.
      enddo

      nfpn = index ( fname , ' ' ) - 1
      open ( 60 , file = fname(1:nfpn) )
      nstars = nsp(nspecies,2)-nsp(nspecies,1)+1
      istar1 = nsp(nspecies,1) 

      tnow = age ( t ) 

c
c.... cell size in proper kpc
c
      rcell = 2.0**(-level) * r0/hubble * aexpn * 1000. 
      write(60,215) (1.0/aexpn-1.0), tnow 
      write(60,217) level, tsfave
      write(60,219) rcell 
      write(60,221)
      write(60,222)
 215  format('# average SFR vs. gas density at z=',
     & f5.3,' t=',f8.5,' [Gyrs]')
 217  format('# averaging level=',i2,
     &'; SFR average over tsfave=',e9.3,' [yrs]')
 219  format('# physical cell size=',e9.3,' [kpc] (proper)')
 221  format('#  ic      SFR3D      SFR2D        rho_g  ',
     &     'sig_g      Tg       Mcell')
 222  format('#         Msun/yr   Msun/yr/pc^2  Msun/pc^3  \ 
     &     Msun/pc^2   K        Msun')

      write(*,*) 'analyzing stellar particles...'
      do is = nsp(nspecies,2), nsp(nspecies,1), -1 
        xx = x(is)
        yy = y(is)
        zz = z(is) 
        ic = iFindcell0 ( level, xx, yy, zz )
c        write(*,*) is, xx, yy, zz, level, ic
        if ( ic .gt. 0 ) then 
          i = is - istar1 + 1
          ts = tnow - age ( tbirth(i)*1.d0 ) ! age of the star [Gyr]         
c          write(*,*) is, ic, ts, tnow, age (tbirth(i)*1.d0 ) , pw0(i) 
          if ( ts .le. tsfaveb ) then 
            ref(ic) = ref(ic) + pw0(i)
          else
            goto 190
          endif
        endif
      enddo

 190  write(*,*) 'done. writing...'
      sfrfact3d = aM0 / tsfave ! Msun yr^-1
      sfrfact2d = aM0 / tsfave / rcell**2 ! Msun yr^-1 kpc^-2
      rhofact = rho0 * 1.d-18 / aexpn**3  ! Msun pc^-3 
      sigfact = rhofact * ( 1000.* rcell ) ! Msun pc^-2 
      tfact = (gamma-1.0) * T_0 / aexpn**2 
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      rfact = r0 * 1.e6 * aexpn / hubble ! /h Mpc -> proper pc
      vfact = v0 / aexpn  ! code velocity -> peculiar in km/s
      amfact = aM0 ! code mass -> Msun
      CellVolume = 2.0**(-3.0*Level)

      call Select_Cells ( Level , nLevel )
      do ic1 = 1 , nLevel
        icell = iSelect(ic1)
        do ic2 = 0 , 7
          i = icell + ic2
          sfr3d = ref(i) * sfrfact3d ! Msun yr^-1
          sfr2d = ref(i) * sfrfact2d ! Msun yr^-1 kpc^-2
          rhocell = hvar(1,i) * rhofact ! Msun pc^-3
          sigcell = hvar(1,i) * sigfact ! Msun pc^-2
          if ( ref(i) .gt. 0.0 ) then 
            rhor = 1.d0 / hvar(1,i)
            dtemp = tfact * hvar(8,i) * rhor
            dmass = hvar(1,i)* CellVolume * aM0
            write(60,200) i, sfr3d, sfr2d, rhocell, sigcell, 
     &           dtemp, dmass
c            write(*,200) i, sfr, rhocell, dtemp, dmass
          endif
        enddo
      enddo
 200  format(i8,6(2x,e9.3))

      close ( 60 )

      return
      end
c
c     -----------------------------------------------------------
      subroutine Compute_Phase ( rhoth, plotname )
c     -----------------------------------------------------------
c
c     compute the phase of gas at cells selected randomly in log radius
c
c     character*256 plotname - full path and file name to output profile to 
c                               it must end with a space
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      real rhoth
      character*256 plotname 
c
      integer m, ic
      parameter ( npoints = 20000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real dg, press, temp, ent, zIa, zII
      common / POINTS / xpn, ypn, zpn
      integer Level, icell, iPack(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn

c     define parameters for cooling rate calculation
      real*8 dt_cool, dt_cool0, dt_cool1
      real*8  rhogl, T_gas, coolrate, fcool, Tfact, Zdum 
      real*8  fact_nH, fact_ex, fact_tau, Thomson, ex, tau
      parameter ( fcool = 1.d0, Thomson = 6.652e-25 )
      real*8 Zsolar, Zsolari
      parameter ( Zsolar = 0.0199 , Zsolari = 1.d0/Zsolar )
      common / COOLCELL / icool, lcool, cool_rate

      
      nfpn = index ( plotname , ' ' ) - 1

      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      ai = 1.0 / aexpn
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 

      m = 2714873

c...  set parameters for cooling rate
      dt_cool1 = 1.0e6
      Tfact = a2i * T_0 / 1.d4  ! T in code units -> T[10^4 Kelvin]
c      dt_cool2 = dt_cool * 2.0**(-Level)
#ifndef CLOUDY_COOLING
      rhogl = 0.0 
#endif
#ifdef CLOUDY_COOLING
c.... factor converting gas density in code units to n_H in cm^-3
      fact_nH = 1.12d-5*hubble*hubble*Om0*(1.d0 - Yp) * ai**3
#endif
      fact_ex = fact_nH**2 * ai / AL_SD * 1.0d-23 
     &     * (r0/hubble*mpc*1.d6)**3
      fact_tau = 1.12d-5*hubble*hubble*Om0*(1.d0-2.*Yp) 
     &     * (r0/hubble*mpc*1.d6) * a2i * Thomson

c
c.... Output to file
c
      open ( 12 , file = plotname(1:nfpn) )
      write(12,13)
      write(12,14)
      write(12,15)
 13   format ('# Columns:')
 14   format ('# Level delta_g  T_g  P_g  S_g  ZIa  ZII Mcell'
     &        '  Emissivity   tcool   tau')
 15   format ('# overdensities w.r.t. mean dens.',
     &        '  [Level]    [K]    [erg cm^-3]     [keV cm^2]  '
     &        '  [Zsun]    [Zsun]   [Msun]   [erg/s/cm^3]  '
     &        '  [Gyr]     [dimensionless]')
c
      DO Level = MinLevel+1 , MaxLevelNow
          CellVolume = 2.0**(-3.0*Level)
          call Select_Cells ( Level , nLevel ) ! nLevel
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              if ( iOctCh(ic) .eq. nil ) then ! look only at leaves 
                call Ps ( ic , Posx , Posy , Posz )
                diff_x = xc - Posx
                diff_y = yc - Posy
                diff_z = zc - Posz
                corr_x = zero
                corr_y = zero
                corr_z = zero

                if ( abs(diff_x) .gt. nf67 ) then
                  if ( diff_x .gt. 0. ) then
                    corr_x = ng
                  else
                    corr_x = -ng
                  endif
                endif
                if ( abs(diff_y) .gt. nf67 ) then
                  if ( diff_y .gt. 0. ) then
                    corr_y = ng
                  else
                    corr_y = -ng
                  endif
                endif
                if ( abs(diff_z) .gt. nf67 ) then
                  if ( diff_z .gt. 0. ) then
                    corr_z = ng
                  else
                    corr_z = -ng
                  endif
                endif

                rhocell = hvar(1,ic)  ! gas overdensity w.r.t mean mass density
                dum = RANDd ( m )

                if ( rhocell .gt. rhoth .and. dum .lt. 0.004 ) then
                  xcell = (Posx - 1.0) * r0
                  ycell = (Posy - 1.0) * r0
                  zcell = (Posz - 1.0) * r0
                  ! compute gas variables of the cell
                  press = hvar(6,ic) * P0 / a5 / 1.38066d-16 ! pressure [K cm^-3]

                  if ( hvar(1,ic) .gt. 0. ) then
                     rhogi = 1.d0 / hvar(1,ic)                     
                     T_code = gamma1 * hvar(8,ic) / hvar(1,ic)
                     temp = T_0 * T_code * a2i
                     ent  = S_0 * T_code / hvar(1,ic)**gamma1
                     zIa = hvar(izIa,ic) * Zsolari * rhogi
                     zII = hvar(izII,ic) * Zsolari * rhogi
                     dmass = hvar(1,ic) * CellVolume * aM0 ! code mass -> Msun
                     tau = hvar(1,ic) * 2.0**(-1.0*Level)
                     

#ifdef CLOUDY_COOLING
c....                rhog_code -> log10(n_H), [n_H]=[cm^-3]
                     rhogl = log10(hvar(1,ic)) + log10(fact_nH)
#endif
                     T_gas = Tfact * hvar(6,ic) / hvar(1,ic)
#ifdef METALCOOLING            
#ifdef ENRICH
                     Zdum = hvar(izII,ic) 
#ifdef ENRICH_SNIa
                     Zdum = Zdum + hvar(izIa,ic)
#endif
                     Zdum = max ( Zdum , 1.d-10 )
                     Zdum = log10(Zdum * Zsolari * rhogi )
#endif
#endif
                     call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate)
                     coolrate = coolrate * hvar(1,ic)**2 * aexpn
                     dt_cool0 = abs(fcool * hvar(8,ic) / coolrate)
                     dt_cool0 = min(dt_cool1,dt_cool0)

c...                 convert code units into physical units
                     ex = coolrate * fact_ex * CellVolume
                     dt_cool0 = dt_cool0 * t0 / 1.0e9 ! [Gyr]
                     tau = tau * fact_tau
                    
                     write(12,16) Level, rhocell, temp, press, ent, 
     &                    zIa, zII, dmass, ex, dt_cool0, tau
                  else
                     write(*,*) 'bad gas density in cell',ic,hvar(1,ic)
                  endif
                endif
              endif
            enddo ! ic2 
          enddo ! ic1
      ENDDO ! Level
 16   format(i2,1x,10(2x,g14.7))
      close ( 12 )

      write(*,*) 'done writing phases of gas. exiting...'
c
      return
      end
c
c     -----------------------------------------------------------
      subroutine Output_Average_Profiles ( xc , yc , zc , rvir, 
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fprofname )
c     -----------------------------------------------------------
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    deltavir : virial overdensity (with respect to mean)
c      integer ifindcenter (see below) 
c      logical lpointsout  - if .true, output sample points in points.dat
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc , rmin , rmax , deltavir
      integer nrbin, ibintype, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pent
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dgbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), dtibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real ndmp(0:nmaxbins),entbin(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) + 1

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 

      call Get_MaxLevelNow ()
c
c.... find a new center 
c
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin
c
c.... select random samples in the volume 
c
      do i = 1 , npoints
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m ) - pi
        phi   = pi2 * RANDd ( m )
        xpn(i) = rr * cos(theta) * cos(phi) + xc
        ypn(i) = rr * cos(theta) * sin(phi) + yc
        zpn(i) = rr * sin(theta) + zc
        if ( xpn(i) .lt. xmin ) xmin = xpn(i)
        if ( xpn(i) .gt. xmax ) xmax = xpn(i)
        if ( ypn(i) .lt. ymin ) ymin = ypn(i)
	if ( ypn(i) .gt. ymax ) ymax = ypn(i)
        if ( zpn(i) .lt. zmin ) zmin = zpn(i)
        if ( zpn(i) .gt. zmax ) zmax = zpn(i)
      enddo
c
      do i = 1 , npoints
        xx  = xpn(i)
        yy  = ypn(i)
        zz  = zpn(i)
        if ( xx .lt. 1.0 ) then 
          xx = xx + ng
        else 
          if ( xx .gt. narr ) then 
            xx = xx - ng
          endif
        endif
        if ( yy .lt. 1.0 ) then 
          yy = yy + ng
        else 
          if ( yy .gt. narr ) then 
            yy = yy - ng
          endif
        endif
        if ( zz .lt. 1.0 ) then 
          zz = zz + ng
        else 
          if ( zz .gt. narr ) then 
            zz = zz - ng
          endif
        endif
        
        idcell = iFindCell ( MaxLevel , xx , yy , zz )
c
c....   gas density
c
        pd(i) = hvar(1,idcell)
c
c....   pressure
c
        pp(i) = hvar(6,idcell)

        if ( hvar(1,idcell) .gt. 0. ) then
          T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
	  ptemp(i) = T_0 * T_code * a2i
          pent(i)  = S_0 * T_code / hvar(1,idcell)**gamma1
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo
c
c.... construct profiles
c
      do i = 0 , nmaxbins
        rbin(i) = 0.
        dgbin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        nbin(i) = 0
        ndmp(i) = 0.
      enddo

      do i = 1 , np
        xx = x(i)
        yy = y(i)
        zz = z(i)
c....   correct coordinate to insure periodic bundaries        
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          ddmbin(ibin) = ddmbin(ibin) + pw(i)
          ndmp(ibin) = ndmp(ibin) + 1
        endif
      enddo

      if ( lpointsout ) then 
        open ( 12 , file = 'points.dat' )
      endif
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin) = pbin(ibin) + pp(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          if ( lpointsout ) then 
            write(12,*) rp*r0, pd(ii), ptemp(ii),pent(ii)
          endif
        endif
      enddo
      if ( lpointsout ) then       
        close ( 12 )
      endif
c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
c          ddmbin(i) = ddmbin(i) / nbin(i)
          pbin(i) = pbin(i) / nbin(i)
          tbin(i) = tbin(i)/nbin(i)
          entbin(i) = entbin(i)/nbin(i)
	endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.0
      dti = 0.0
      ivirflag = 0 

      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
        else
          rl = 10.**(rlmin + float(i-1)*drl)
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
	if ( volr .gt. 0. ) then
          dmi = dmi + ddmbin(i) 
          dti = dti + ddmbin(i) + dgbin(i)*volr 
          ddmbin(i) = ddmbin(i)/volr 
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i) = dti / (4.18879 * rr**3)
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. deltavir .and.
     &         dtibin(i-1) .ge. deltavir ) then 
            ivir = i
            rvir = (deltavir * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
c            rvir = rr 
          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
      enddo
c
      open ( 12 , file = fprofname(1:nfpn) )
c
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype 
      write(12,12) ivir, rvir*r0, aMvir, deltavir
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r delta_g delta_dm deltai_dm',
     &        ' deltai_tot  T_g  S_g     N_dm_p' )
 16   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '             K  keV cm^2 ')
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) * r0
        endif
        rr = 10.**(rlmin + float(i)*drl)
	if ( nbin(i) .gt. 0 ) then
          write(12,17)
     &         rbin(i), rl*r0, rmid, rr*r0, 
     &         dgbin(i),ddmbin(i),ddmibin(i), dtibin(i), 
     &         tbin(i),entbin(i),int(ndmp(i))
	endif
      enddo
 17   format(10(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end
c
c     -----------------------------------------------------------
      subroutine Output_Average_Profiles_CSF ( xc , yc , zc , rvir, 
     &                            rcrit, rmin , rmax , nrbin , 
     &                            ibintype , deltavir, deltacrit, 
     &                            ifindcenter , lpointsout, 
     &                            l_force_center, fprofname )
c     -----------------------------------------------------------
c
c     compute the average profiles for simulations with cooling + star fomration
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    deltavir : virial overdensity (with respect to mean)
c      integer ifindcenter (see below) 
c      logical lpointsout  - if .true, output sample points in points.dat
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc , rmin , rmax
      real rvir, rcrit, deltavir, deltacrit
      integer nrbin, ibintype, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints),pzIa(npoints),pzII(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pent
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dtibin(0:nmaxbins)
      real dgbin(0:nmaxbins),dstbin(0:nmaxbins),ddmbin(0:nmaxbins)
      real dgibin(0:nmaxbins),dstibin(0:nmaxbins),ddmibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real zgIabin(0:nmaxbins),zgIIbin(0:nmaxbins)
      real zstIabin(0:nmaxbins),zstIIbin(0:nmaxbins)
      real ndmp(0:nmaxbins),nstp(0:nmaxbins),entbin(0:nmaxbins)
c     bins for cold gas fraction
      real mgi,mgi1,mgi2,mgi3,mgi4
      real mgbin(0:nmaxbins),mg1bin(0:nmaxbins),mg2bin(0:nmaxbins),
     +     mg3bin(0:nmaxbins),mg4bin(0:nmaxbins),mgi1bin(0:nmaxbins),
     +     mgi2bin(0:nmaxbins),mgi3bin(0:nmaxbins),mgi4bin(0:nmaxbins),
     +     mzgIabin(0:nmaxbins),mzgIIbin(0:nmaxbins),
     +     dzgIaibin(0:nmaxbins),dzgIIibin(0:nmaxbins),
     +     dzstIaibin(0:nmaxbins),dzstIIibin(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) + 1

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 

      call Get_MaxLevelNow ()
c
c.... find a new center 
c
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin
c
c.... select random samples in the volume 
c
      do i = 1 , npoints
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m ) - pi
        phi   = pi2 * RANDd ( m )
        xpn(i) = rr * cos(theta) * cos(phi) + xc
        ypn(i) = rr * cos(theta) * sin(phi) + yc
        zpn(i) = rr * sin(theta) + zc
        if ( xpn(i) .lt. xmin ) xmin = xpn(i)
        if ( xpn(i) .gt. xmax ) xmax = xpn(i)
        if ( ypn(i) .lt. ymin ) ymin = ypn(i)
	if ( ypn(i) .gt. ymax ) ymax = ypn(i)
        if ( zpn(i) .lt. zmin ) zmin = zpn(i)
        if ( zpn(i) .gt. zmax ) zmax = zpn(i)
      enddo
c
      do i = 1 , npoints
        xx  = xpn(i)
        yy  = ypn(i)
        zz  = zpn(i)
        if ( xx .lt. 1.0 ) then 
          xx = xx + ng
        else 
          if ( xx .gt. narr ) then 
            xx = xx - ng
          endif
        endif
        if ( yy .lt. 1.0 ) then 
          yy = yy + ng
        else 
          if ( yy .gt. narr ) then 
            yy = yy - ng
          endif
        endif
        if ( zz .lt. 1.0 ) then 
          zz = zz + ng
        else 
          if ( zz .gt. narr ) then 
            zz = zz - ng
          endif
        endif
        
        idcell = iFindCell ( MaxLevel , xx , yy , zz )
c
c....   gas density
c
        pd(i) = hvar(1,idcell)
c
c....   pressure
c
        pp(i) = hvar(6,idcell)

        if ( hvar(1,idcell) .gt. 0. ) then
          T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
	  ptemp(i) = T_0 * T_code * a2i
          pent(i)  = S_0 * T_code / hvar(1,idcell)**gamma1
          pzIa(i) = hvar(izIa,idcell) / 0.01989 / hvar(1,idcell)
          pzII(i) = hvar(izII,idcell) / 0.01989 / hvar(1,idcell)
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo
c
c.... construct profiles
c
      do i = 0 , nmaxbins
        rbin(i) = 0.
        dgbin(i) = 0.
        dstbin(i) = 0.
        ddmbin(i) = 0.
        dgibin(i) = 0.
        dstibin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        dzgIaibin(i) = 0.
        dzgIIibin(i) = 0.
        dzstIaibin(i) = 0.
        dzstIIibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        zgIabin(i) = 0.
        zgIIbin(i) = 0.
        zstIabin(i) = 0.
        zstIIbin(i) = 0.
        nbin(i) = 0
        ndmp(i) = 0

        mgbin(i) = 0.  
        mg1bin(i) = 0. 
        mg2bin(i) = 0.
        mg3bin(i) = 0.
        mg4bin(i) = 0.
        mzgIabin(i) = 0.
        mzgIIbin(i) = 0.
      enddo

c
c...  compute profiles of DM and Stars
c
      nstars = nsp(nspecies,2) - nsp(nspecies,1) + 1
      istar1 = nsp(nspecies,1)

      do i = 1 , np
        xx = x(i)
        yy = y(i)
        zz = z(i)
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
        ! compute DM and Star profiles
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          if ( i .lt. istar1 ) then ! DM particle
            ddmbin(ibin) = ddmbin(ibin) + pw(i)
            ndmp(ibin) = ndmp(ibin) + 1
          else ! stellar particle
            dstbin(ibin) = dstbin(ibin) + pw(i)
            zstIabin(ibin) = zstIabin(ibin) + zstIa(i-istar1+1)*pw(i)
            zstIIbin(ibin) = zstIIbin(ibin) + zstII(i-istar1+1)*pw(i)
            nstp(ibin) = nstp(ibin) + 1
          endif
        endif
      enddo

c...  compute gas profiles

      if ( lpointsout ) then 
        open ( 12 , file = 'points.dat' )
      endif
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin) = pbin(ibin) + pp(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          zgIabin(ibin) = zgIabin(ibin) + pzIa(ii) * pd(ii)
          zgIIbin(ibin) = zgIIbin(ibin) + pzII(ii) * pd(ii)
          if ( lpointsout ) then 
            write(12,*) rp*r0, pd(ii), ptemp(ii),pent(ii)
          endif
        endif
      enddo
      if ( lpointsout ) then       
        close ( 12 )
      endif

c
c...  compute cold gas fraction for 4 different temperature     
c

      ! temperature thresholds for computing cold gas fraction
      tg1 = 10000.
      tg2 = 20000.
      tg3 = 100000.
      tg4 = 1000000.

      call Get_MaxLevelNow ()

      DO Level = MinLevel , MaxLevelNow
         CellVolume = 1.0 * 2.0**(-3.0*Level)
         IF ( Level .eq. MinLevel ) THEN
         ELSE
            call Select_Cells ( Level , nLevel )
            do ic1 = 1 , nLevel
              icell = iSelect(ic1)
              do ic2 = 0 , 7
                ic = icell + ic2
                if ( iOctCh(ic) .eq. nil ) then
                  call Ps ( ic , Posx , Posy , Posz )
                  call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz )
                  rp = sqrt( (Posx-xc)**2 +
     &                       (Posy-yc)**2 +
     &                       (Posz-zc)**2 )
                  ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
                  cellmass = hvar(1,ic) * CellVolume
                  zgIamass = hvar(izIa,ic) * CellVolume
                  zgIImass = hvar(izII,ic) * CellVolume
                  if ( hvar(1,ic) .gt. 0. ) then
                     Tcell = gamma1 * hvar(8,ic) / hvar(1,ic)
                     Tcell = T_0 * Tcell * a2i
                  else
                     write(*,*) 'bad gas density in cell',ic,hvar(1,c)
                  endif
                  if ( ibin .le. nmaxbins ) then
                     mgbin(ibin) = mgbin(ibin) + cellmass
                     mzgIabin(ibin) = mzgIabin(ibin) + zgIamass
                     mzgIIbin(ibin) = mzgIIbin(ibin) + zgIImass
                     if ( Tcell .lt. tg1 ) 
     &                    mg1bin(ibin) = mg1bin(ibin) + cellmass
                     if ( Tcell .lt. tg2 ) 
     &                    mg2bin(ibin) = mg2bin(ibin) + cellmass
                     if ( Tcell .lt. tg3 ) 
     &                    mg3bin(ibin) = mg3bin(ibin) + cellmass
                     if ( Tcell .lt. tg4 ) 
     &                    mg4bin(ibin) = mg4bin(ibin) + cellmass
                  endif

c                  if ( Tcell .lt. tg3 .and. rp .gt. 0.50 
c     &                 .and. rp .lt. 1.6 ) then
c                      write(*,*) Tcell, rp, Posx, Posy, Posz
c                  endif

                endif
              enddo
            enddo
         ENDIF
      ENDDO

c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          ! compute average gas profiles 
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
          pbin(i) = pbin(i) / nbin(i)
          tbin(i) = tbin(i) / nbin(i)
          entbin(i) = entbin(i) / nbin(i)
          ! mass-weighted metallcity profiles of gas
          zgIabin(i) = zgIabin(i) / dgbin(i) / nbin(i)
          zgIIbin(i) = zgIIbin(i) / dgbin(i) / nbin(i)
	endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      ddmi = 0.0
      dsti = 0.0
      dgi  = 0.0
      dti  = 0.0
      mgi  = 1.0e-30
      mgi1 = 1.0e-30
      mgi2 = 1.0e-30
      mgi3 = 1.0e-30
      mgi4 = 1.0e-30
      dzgIai = 0.0
      dzgIIi = 0.0
      dzstIai = 0.0
      dzstIIi = 0.0

      dcritave = deltacrit * aexpn**3 / Om0
      ivirflag = 0 
   
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
        else
          rl = 10.**(rlmin + float(i-1)*drl)
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
	if ( volr .gt. 0. ) then
          ! mass of DM, star, gas, total in code units
          ddmi = ddmi + ddmbin(i) 
          dsti = dsti + dstbin(i)
          dgi = dgi + mgbin(i)
          dti = dti + ddmbin(i) + mgbin(i) + dstbin(i)
          ! mass of metals in gas and stars 
          dzgIai = dzgIai + mzgIabin(i)
          dzgIIi = dzgIIi + mzgIIbin(i)
          dzstIai = dzstIai + zstIabin(i)
          dzstIIi = dzstIIi + zstIIbin(i)
          ! mass-weighted stellar metallicity profiles of SNIa & II
          zstIabin(i) = zstIabin(i) / dstbin(i) / 0.01989
          zstIIbin(i) = zstIIbin(i) / dstbin(i) / 0.01989
          ! average overdensity of DM & star
          ddmbin(i) = ddmbin(i) / volr 
          dstbin(i) = dstbin(i) / volr
          ! cummulative overdensity profiles of DM, star, gas, total
          ddmibin(i) = ddmi / (4.18879 * rr**3)
          dstibin(i) = dsti / (4.18879 * rr**3)
          dgibin(i)  = dgi  / (4.18879 * rr**3)
          dtibin(i)  = dti  / (4.18879 * rr**3)
          dzgIaibin(i)  = dzgIai  / (4.18879 * rr**3)
          dzgIIibin(i)  = dzgIIi  / (4.18879 * rr**3)
          dzstIaibin(i) = dzstIai  / (4.18879 * rr**3)
          dzstIIibin(i) = dzstIIi  / (4.18879 * rr**3)
          ! compute cold gas mass fraction : Mgas(T<T*) / Mgas
          mgi  = mgi  + mgbin(i)
          mgi1 = mgi1 + mg1bin(i)
          mgi2 = mgi2 + mg2bin(i)
          mgi3 = mgi3 + mg3bin(i)
          mgi4 = mgi4 + mg4bin(i)
          mgi1bin(i) = mgi1 / mgi
          mgi2bin(i) = mgi2 / mgi
          mgi3bin(i) = mgi3 / mgi
          mgi4bin(i) = mgi4 / mgi
          ! compute virial radius
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. deltavir .and.
     &         dtibin(i-1) .ge. deltavir ) then 
            ivir = i
            rvir = (deltavir * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
c            rvir = rr 
          endif
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. dcritave .and.
     &         dtibin(i-1) .ge. dcritave ) then 
            icrit = i
            rcrit = (dcritave * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
      enddo

c
      open ( 12 , file = fprofname(1:nfpn) )
c
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      aMcrit = dcritave * 
     &        4.18879 * (rcrit*r0)**3 * rho0 / hubble**2 
      rr = 10.**(rlmin + float(ivir)*drl)
      aM_dm  = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
      aM_gas = dgibin(ivir)  * (4.18879 * rr**3) * aM0 * hubble
      aM_st  = dstibin(ivir) * (4.18879 * rr**3) * aM0 * hubble 
      aMvir2 = aM_dm + aM_gas + aM_st 
      
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype 
      write(12,12) ivir, rvir*r0, aMvir, deltavir
      write(12,22) icrit, rcrit*r0, aMcrit, deltacrit
      write(12,17) aM_dm, aM_gas, aM_st, aMvir2
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 22   format ('# icrit = ',i4,' Rcrit =',f7.4,
     &        ' [h^-1 Mpc], Mcrit =',g11.4,
     &  ' [h^-1 Msun] for Dcrit(overd w.r.t. rho_crit)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r  delta_g delta_st delta_dm  deltai_g',
     &        ' deltai_st deltai_dm  fg(T<10^4K) fg(T<2*10^4K)',
     &        ' fg(T<10^5K) fg(T<10^6K)   T_g  S_g  ZIa  ZII  N_dm_p' )
 16   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '             K  keV cm^2 ')
 17   format ('# Mdm = ',g11.4,' Mgas =',g11.4,
     &        ' Mstars =',g11.4,' Mtot = Mdm + Mgas + Mstars =',g11.4,
     &        ' [h^-1 Msun]')
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) * r0
        endif
        rr = 10.**(rlmin + float(i)*drl)
	if ( nbin(i) .gt. 0 ) then
          write(12,18)
     &         rbin(i), rl*r0, rmid, rr*r0, 
     &         dgbin(i),dstbin(i),ddmbin(i),dgibin(i),dstibin(i),
     &         ddmibin(i),mgi1bin(i),mgi2bin(i),mgi3bin(i),mgi4bin(i),
     &         tbin(i),entbin(i),
     &         zgIabin(i),zgIIbin(i),zstIabin(i),zstIIbin(i),
     &         dzgIaibin(i),dzgIIibin(i),dzstIaibin(i),dzstIIibin(i),
     &         int(ndmp(i))
	endif
      enddo
 18   format(24(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end
c
c     -----------------------------------------------------------
      subroutine Output_Average_Profiles_Stars (
     &                            xc , yc , zc , rvir, 
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, deltacrit, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fprofname )
c     -----------------------------------------------------------
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    deltavir : virial overdensity (with respect to mean)
c              deltacrit : virial overdensity with respect to the critical
c      integer ifindcenter (see below) 
c      logical lpointsout  - if .true, output sample points in points.dat
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc , rmin , rmax , deltavir, deltacrit
      integer nrbin, ibintype, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pent
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dgbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), dtibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real ndmp(0:nmaxbins),entbin(0:nmaxbins)
      real*8 amst(0:nmaxbins), amstot(0:nmaxbins), nst(0:nmaxbins)
      real*8 amgtot(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) + 1

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 

      call Get_MaxLevelNow ()
c
c.... find a new center 
c
      
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin
c
c.... select random samples in the volume 
c
      do i = 1 , npoints
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m ) - pi
        phi   = pi2 * RANDd ( m )
        xpn(i) = rr * cos(theta) * cos(phi) + xc
        ypn(i) = rr * cos(theta) * sin(phi) + yc
        zpn(i) = rr * sin(theta) + zc
        if ( xpn(i) .lt. xmin ) xmin = xpn(i)
        if ( xpn(i) .gt. xmax ) xmax = xpn(i)
        if ( ypn(i) .lt. ymin ) ymin = ypn(i)
	if ( ypn(i) .gt. ymax ) ymax = ypn(i)
        if ( zpn(i) .lt. zmin ) zmin = zpn(i)
        if ( zpn(i) .gt. zmax ) zmax = zpn(i)
      enddo
c
      do i = 1 , npoints
        xx  = xpn(i)
        yy  = ypn(i)
        zz  = zpn(i)
        if ( xx .lt. 1.0 ) then 
          xx = xx + ng
        else 
          if ( xx .gt. narr ) then 
            xx = xx - ng
          endif
        endif
        if ( yy .lt. 1.0 ) then 
          yy = yy + ng
        else 
          if ( yy .gt. narr ) then 
            yy = yy - ng
          endif
        endif
        if ( zz .lt. 1.0 ) then 
          zz = zz + ng
        else 
          if ( zz .gt. narr ) then 
            zz = zz - ng
          endif
        endif
        
        idcell = iFindCell ( MaxLevel , xx , yy , zz )
c
c....   gas density
c
        pd(i) = hvar(1,idcell)
c
c....   pressure
c
        pp(i) = hvar(6,idcell)

        if ( hvar(1,idcell) .gt. 0. ) then
          T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
	  ptemp(i) = T_0 * T_code * a2i
          pent(i)  = S_0 * T_code / hvar(1,idcell)**gamma1
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo
c
c.... construct profiles
c
      do i = 0 , nmaxbins
        rbin(i) = 0.
        dgbin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        nbin(i) = 0
        ndmp(i) = 0.
        nst(i) = 0
        amst(i) = 0
        amstot(i) = 0
        amgtot(i) = 0.
      enddo

      nstars = nsp(nspecies,2) - nsp(nspecies,1) + 1
      istar1 = nsp(nspecies,1)

      do i = 1 , np 
        xx = x(i)
        yy = y(i)
        zz = z(i)
        ! correct coordinate to insure periodic bundaries
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          if ( i .lt. istar1 ) then ! DM particle
            ddmbin(ibin) = ddmbin(ibin) + pw(i)
            ndmp(ibin) = ndmp(ibin) + 1
          else ! stellar particle
            amst(ibin) = amst(ibin) + pw(i)
            nst(ibin) = nst(ibin) + 1
          endif
        endif
      enddo

      if ( lpointsout ) then 
        open ( 12 , file = 'points.dat' )
      endif
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin) = pbin(ibin) + pp(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          if ( lpointsout ) then 
            write(12,*) rp*r0, pd(ii), ptemp(ii),pent(ii)
          endif
        endif
      enddo
      if ( lpointsout ) then       
        close ( 12 )
      endif
c
c.... prepare profiles
c
      do i = 0 , nbins
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
          pbin(i) = pbin(i) / nbin(i)
          tbin(i) = tbin(i)/nbin(i)
          entbin(i) = entbin(i)/nbin(i)
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.0
      dti = 0.0
      dcritave = deltacrit * aexpn**3 / Om0

      amst_tot = 0.0 
      amg_tot = 0.0 
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
        else
          rl = 10.**(rlmin + float(i-1)*drl)
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
	if ( volr .gt. 0. ) then
          dmi = dmi + ddmbin(i) 
          amst_tot = amst_tot + amst(i) 
          amstot(i) = amst_tot
          amg_tot = amg_tot + dgbin(i)*volr
          amgtot(i) = amg_tot
          dti = dti + ddmbin(i) + dgbin(i)*volr + amst(i)
          ddmbin(i) = ddmbin(i)/volr 
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i) = dti / (4.18879 * rr**3)
          if ( ibin .gt. 0 .and.
     &         dtibin(i) .lt. deltavir .and.
     &         dtibin(i-1) .ge. deltavir ) then 
            ivir = i
            rvir = (deltavir * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
c            rvir = rr 
          endif
          if ( ibin .gt. 0 .and.
     &         dtibin(i) .lt. dcritave .and.
     &         dtibin(i-1) .ge. dcritave ) then 
            icrit = i
            rcrit = (dcritave * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
      enddo
c
      open ( 12 , file = fprofname(1:nfpn) )
c
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      aMcrit = dcritave * 
     &        4.18879 * (rcrit*r0)**3 * rho0 / hubble**2 
      aM_gas = amgtot(ivir) * aM0 * hubble
      rr = 10.**(rlmin + float(ivir)*drl)
      aM_dm = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
      aM_st = amstot(ivir) * aM0 * hubble 
      aMvir2 = aM_dm + aM_gas + aM_st 


      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype 
      write(12,12) ivir, rvir*r0, aMvir, deltavir
      write(12,22) icrit, rcrit*r0, aMcrit, deltacrit
      write(12,17) aM_dm, aM_gas, aM_st, aMvir2
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g11.4,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 22   format ('# icrit = ',i4,' Rcrit =',f7.4,
     &        ' [h^-1 Mpc], Mcrit =',g11.4,
     &  ' [h^-1 Msun] for Dcrit(overd w.r.t. rho_crit)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r d_g d_* d_dm m_g m_* m_dm',
     &        ' T S N_dm N_*' )
 16   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '  [m]=[/h Msun]           K  keV cm^2 ')
 17   format ('# Mdm = ',g11.4,' Mgas =',g11.4,
     &        ' Mstars =',g11.4,' Mtot = Mdm + Mgas + Mstars =',g11.4,
     &        ' [h^-1 Msun]')
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) * r0
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        

        aMdm = ddmibin(i) * (4.18879 * rr**3) * aM0 * hubble
        rhost = amst(i) / volr 
          write(12,18)
     &         rbin(i), rl*r0, rmid, rr*r0, 
     &         dgbin(i),rhost, ddmbin(i), 
     &         amgtot(i)*aM0*hubble, 
     &         amstot(i)*aM0*hubble,
     &         aMdm,
     &         tbin(i),entbin(i),int(ndmp(i)),int(nst(i))
      enddo
 18   format(4(g11.4,1x),1x,8(e11.4,1x), 2(1x,i7))
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end
c
c     -----------------------------------------------------------
      subroutine Output_Fit_Profiles ( xc , yc , zc , rvir, 
     &                            rcrit, rmin , rmax , nrbin , 
     &                            ibintype , deltavir, deltacrit, 
     &                            ifindcenter , lpointsout, 
     &                            l_force_center, fprofname )
c     -----------------------------------------------------------
c
c     compute the average profiles for simulations with cooling + star fomration
c     for profile fitting
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    deltavir : virial overdensity (with respect to mean)
c      integer ifindcenter (see below) 
c      logical lpointsout  - if .true, output sample points in points.dat
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c

c
      real xc , yc , zc , rmin , rmax
      real rvir, rcrit, deltavir, deltacrit
      integer nrbin, ibintype, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints),pzIa(npoints),pzII(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pent
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dtibin(0:nmaxbins)
      real dgbin(0:nmaxbins),dstbin(0:nmaxbins),ddmbin(0:nmaxbins)
      real dgibin(0:nmaxbins),dstibin(0:nmaxbins),ddmibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real zIabin(0:nmaxbins),zIIbin(0:nmaxbins)
      real ndmp(0:nmaxbins),nstp(0:nmaxbins),entbin(0:nmaxbins)
c     bins for cold gas fraction
      real mgi,mgi1,mgi2,mgi3,mgi4
      real mgbin(0:nmaxbins),mg1bin(0:nmaxbins),mg2bin(0:nmaxbins),
     +     mg3bin(0:nmaxbins),mg4bin(0:nmaxbins),mgi1bin(0:nmaxbins),
     +     mgi2bin(0:nmaxbins),mgi3bin(0:nmaxbins),mgi4bin(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) + 1

      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a5i = 1.0 / a5
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 

      call Get_MaxLevelNow ()
c
c.... find a new center 
c
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin
c
c.... select random samples in the volume 
c
      do i = 1 , npoints
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m ) - pi
        phi   = pi2 * RANDd ( m )
        xpn(i) = rr * cos(theta) * cos(phi) + xc
        ypn(i) = rr * cos(theta) * sin(phi) + yc
        zpn(i) = rr * sin(theta) + zc
        if ( xpn(i) .lt. xmin ) xmin = xpn(i)
        if ( xpn(i) .gt. xmax ) xmax = xpn(i)
        if ( ypn(i) .lt. ymin ) ymin = ypn(i)
	if ( ypn(i) .gt. ymax ) ymax = ypn(i)
        if ( zpn(i) .lt. zmin ) zmin = zpn(i)
        if ( zpn(i) .gt. zmax ) zmax = zpn(i)
      enddo
c
      do i = 1 , npoints
        xx  = xpn(i)
        yy  = ypn(i)
        zz  = zpn(i)
        if ( xx .lt. 1.0 ) then 
          xx = xx + ng
        else 
          if ( xx .gt. narr ) then 
            xx = xx - ng
          endif
        endif
        if ( yy .lt. 1.0 ) then 
          yy = yy + ng
        else 
          if ( yy .gt. narr ) then 
            yy = yy - ng
          endif
        endif
        if ( zz .lt. 1.0 ) then 
          zz = zz + ng
        else 
          if ( zz .gt. narr ) then 
            zz = zz - ng
          endif
        endif
        
        idcell = iFindCell ( MaxLevel , xx , yy , zz )
c
c....   gas density
c
        pd(i) = hvar(1,idcell)
c
c....   pressure
c
        pp(i) = P0 * hvar(6,idcell) * a5i

        if ( hvar(1,idcell) .gt. 0. ) then
          T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
	  ptemp(i) = T_0 * T_code * a2i
          pent(i)  = S_0 * T_code / hvar(1,idcell)**gamma1
          pzIa(i) = hvar(izIa,idcell) / 0.01989 / hvar(1,idcell)
          pzII(i) = hvar(izII,idcell) / 0.01989 / hvar(1,idcell)
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo
c
c.... construct profiles
c
      do i = 0 , nmaxbins
        rbin(i) = 0.
        dgbin(i) = 0.
        dstbin(i) = 0.
        ddmbin(i) = 0.
        dgibin(i) = 0.
        dstibin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        zIabin(i) = 0.
        zIIbin(i) = 0.
        nbin(i) = 0
        ndmp(i) = 0

        mgbin(i) = 0.  
        mg1bin(i) = 0. 
        mg2bin(i) = 0.
        mg3bin(i) = 0.
        mg4bin(i) = 0.
      enddo

c
c...  compute profiles of DM and Stars
c
      nstars = nsp(nspecies,2) - nsp(nspecies,1) + 1
      istar1 = nsp(nspecies,1)

      do i = 1 , np
        xx = x(i)
        yy = y(i)
        zz = z(i)
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
        ! compute DM and Star profiles
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          if ( i .lt. istar1 ) then ! DM particle
            ddmbin(ibin) = ddmbin(ibin) + pw(i)
            ndmp(ibin) = ndmp(ibin) + 1
          else ! stellar particle
            dstbin(ibin) = dstbin(ibin) + pw(i)
            nstp(ibin) = nstp(ibin) + 1
          endif
        endif
      enddo

c...  compute gas profiles

      if ( lpointsout ) then 
        open ( 12 , file = 'points.dat' )
      endif
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin) = pbin(ibin) + pp(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          zIabin(ibin) = zIabin(ibin) + pzIa(ii) * pd(ii)
          zIIbin(ibin) = zIIbin(ibin) + pzII(ii) * pd(ii)
          if ( lpointsout ) then 
            write(12,*) rp*r0, pd(ii), ptemp(ii),pent(ii)
          endif
        endif
      enddo
      if ( lpointsout ) then       
        close ( 12 )
      endif

c
c...  compute cold gas fraction for 4 different temperature     
c

      ! temperature thresholds for computing cold gas fraction
      tg1 = 10000.
      tg2 = 20000.
      tg3 = 100000.
      tg4 = 1000000.

      call Get_MaxLevelNow ()

      DO Level = MinLevel , MaxLevelNow
         CellVolume = 1.0 * 2.0**(-3.0*Level)
         IF ( Level .eq. MinLevel ) THEN
         ELSE
            call Select_Cells ( Level , nLevel )
            do ic1 = 1 , nLevel
              icell = iSelect(ic1)
              do ic2 = 0 , 7
                ic = icell + ic2
                if ( iOctCh(ic) .eq. nil ) then
                  call Ps ( ic , Posx , Posy , Posz )
                  call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz )
                  rp = sqrt( (Posx-xc)**2 +
     &                       (Posy-yc)**2 +
     &                       (Posz-zc)**2 )
                  ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
                  cellmass = hvar(1,ic) * CellVolume
                  if ( hvar(1,ic) .gt. 0. ) then
                     Tcell = gamma1 * hvar(8,ic) / hvar(1,ic)
                     Tcell = T_0 * Tcell * a2i
                  else
                     write(*,*) 'bad gas density in cell',ic,hvar(1,c)
                  endif
                  if ( ibin .le. nmaxbins ) then
                     mgbin(ibin) = mgbin(ibin) + cellmass
                     if ( Tcell .lt. tg1 ) 
     &                    mg1bin(ibin) = mg1bin(ibin) + cellmass
                     if ( Tcell .lt. tg2 ) 
     &                    mg2bin(ibin) = mg2bin(ibin) + cellmass
                     if ( Tcell .lt. tg3 ) 
     &                    mg3bin(ibin) = mg3bin(ibin) + cellmass
                     if ( Tcell .lt. tg4 ) 
     &                    mg4bin(ibin) = mg4bin(ibin) + cellmass
                  endif

c                  if ( Tcell .lt. tg3 .and. rp .gt. 0.50 
c     &                 .and. rp .lt. 1.6 ) then
c                      write(*,*) Tcell, rp, Posx, Posy, Posz
c                  endif

                endif
              enddo
            enddo
         ENDIF
      ENDDO

c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          ! compute average gas profiles 
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
          pbin(i) = pbin(i) / nbin(i)
          tbin(i) = tbin(i) / nbin(i)
          entbin(i) = entbin(i) / nbin(i)
          zIabin(i) = zIabin(i) / dgbin(i) / nbin(i)
          zIIbin(i) = zIIbin(i) / dgbin(i) / nbin(i)
	endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      ddmi = 0.0
      dsti = 0.0
      dgi  = 0.0
      dti  = 0.0
      mgi  = 1.0e-30
      mgi1 = 1.0e-30
      mgi2 = 1.0e-30
      mgi3 = 1.0e-30
      mgi4 = 1.0e-30

      dcritave = deltacrit * aexpn**3 / Om0
      ivirflag = 0 
   
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
        else
          rl = 10.**(rlmin + float(i-1)*drl)
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
	if ( volr .gt. 0. ) then
          ! mass of DM, star, gas, total in code units
          ddmi = ddmi + ddmbin(i) 
          dsti = dsti + dstbin(i)
          dgi = dgi + mgbin(i)
          dti = dti + ddmbin(i) + mgbin(i) + dstbin(i)
          ! average overdensity of DM, star, gas
          ddmbin(i) = ddmbin(i) / volr 
          dstbin(i) = dstbin(i) / volr
          ! cummulative overdensity profiles of DM, star, gas, total
          ddmibin(i) = ddmi / (4.18879 * rr**3)
          dstibin(i) = dsti / (4.18879 * rr**3)
          dgibin(i)  = dgi  / (4.18879 * rr**3)
          dtibin(i)  = dti  / (4.18879 * rr**3)
          ! compute cold gas mass fraction : Mgas(T<T*) / Mgas
          mgi  = mgi  + mgbin(i)
          mgi1 = mgi1 + mg1bin(i)
          mgi2 = mgi2 + mg2bin(i)
          mgi3 = mgi3 + mg3bin(i)
          mgi4 = mgi4 + mg4bin(i)
          mgi1bin(i) = mgi1 / mgi
          mgi2bin(i) = mgi2 / mgi
          mgi3bin(i) = mgi3 / mgi
          mgi4bin(i) = mgi4 / mgi
          ! compute virial radius
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. deltavir .and.
     &         dtibin(i-1) .ge. deltavir ) then 
            ivir = i
            rvir = (deltavir * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
c            rvir = rr 
          endif
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. dcritave .and.
     &         dtibin(i-1) .ge. dcritave ) then 
            icrit = i
            rcrit = (dcritave * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
      enddo

c
      open ( 12 , file = fprofname(1:nfpn) )
c
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      aMcrit = dcritave * 
     &        4.18879 * (rcrit*r0)**3 * rho0 / hubble**2 
      rr = 10.**(rlmin + float(ivir)*drl)
      aM_dm  = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
      aM_gas = dgibin(ivir)  * (4.18879 * rr**3) * aM0 * hubble
      aM_st  = dstibin(ivir) * (4.18879 * rr**3) * aM0 * hubble 
      aMvir2 = aM_dm + aM_gas + aM_st 

      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype 
      write(12,12) ivir, rvir*r0, aMvir, deltavir
      write(12,22) icrit, rcrit*r0, aMcrit, deltacrit
      write(12,17) aM_dm, aM_gas, aM_st, aMvir2
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 22   format ('# icrit = ',i4,' Rcrit =',f7.4,
     &        ' [h^-1 Mpc], Mcrit =',g11.4,
     &  ' [h^-1 Msun] for Dcrit(overd w.r.t. rho_crit)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r  delta_g delta_st delta_dm  deltai_g',
     &        ' deltai_st deltai_dm  T_g  P_g S_g N_dm_p' )
 16   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '             K   erg cm^-3   keV cm^2 ')
 17   format ('# Mdm = ',g11.4,' Mgas =',g11.4,
     &        ' Mstars =',g11.4,' Mtot = Mdm + Mgas + Mstars =',g11.4,
     &        ' [h^-1 Msun]')
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) * r0
        endif
        rr = 10.**(rlmin + float(i)*drl)
	if ( nbin(i) .gt. 0 ) then
          write(12,18)
     &         rbin(i), rl*r0, rmid, rr*r0, 
     &         dgbin(i),dstbin(i),ddmbin(i),dgibin(i),dstibin(i),
     &         ddmibin(i),tbin(i),pbin(i),entbin(i),int(ndmp(i))
          write(*,*) entbin(i)
	endif
      enddo
 18   format(13(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end
