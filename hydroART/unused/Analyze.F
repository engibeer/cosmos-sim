c     ===============================================================
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
c     Last updated : 9/29/03 (D. Nagai : daisuke@oddjob.uchicago.edu)
c     
c     ==============================================================
c
      program ART
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'

      real*8 dt, rs
      character*5 fstep
      character*256 fname, hname, cldir, clname
      character*256 rs_dir, rs_fname, rs_spec_fname 
      integer lfrs
      logical lpointsout, l_force_center
c
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*5 chbox
c
      character*80 ctmp
      integer IARGC
      external IARGC

c
c.... read/set in initial conditions, initialize variables and linked lists
c
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'starting, reading control parameters ...'
      close(iOWork)

      call ReadControl ()
      call Initialize ( dt )


#ifdef CLOUDY_COOLING

      rs = max ( 1.d0 / aexpn - 1.0 , 0.d0 ) 

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'setting Cloudy cooling rate table for z=',rs
         close(iOWork)
#endif
      call Set_Cooling_Rate_rs ( rs ) ! prepare tables for cooling 
#endif

c
c...  set flags for analyses (1 - ON, 0- off)
c

      ! Output 3D & 2D data
      iflag_datacube          = 0 ! output data cube
      iflag_dataout_2dmap     = 0 ! output 2D map
      
      ! X-ray 
      iflag_pro_xray3d        = 0 ! 3D X-ray profiles
      iflag_pro_xray2d        = 0 ! 2D X-ray profiles
      iflag_pro_xray2d_smooth = 0 ! smoothed 2D X-ray profiles
      iflag_compute_xray      = 0 ! compute X-ray properties
      iflag_xspectra          = 0 ! compute X-ray spectra in annuli
      iflag_xrayvir           = 0 ! X-ray properties at virial radius

      ! SZ
      iflag_pro_sz3d          = 0 ! 3D SZ profiles
      iflag_pro_sz2d          = 0 ! 2D SZ profiles
      iflag_szvir             = 0 ! SZ properties at virial radius

      ! Halo properties
      iflag_pro_csf           = 0 ! compute average profiles
      iflag_pro_fit           = 0 ! profiles for fitting
      iflag_pro_HSE           = 0 ! profiles of HSE mass and the mach number
      iflag_gas_velocities    = 0 ! output data file with gas velocities
      iflag_sfr_rvir          = 0 ! compute SFR for stellar particles within < Rvir
      iflag_sfr_r1_r2         = 0 ! compute SFR for stellar particles within r1 < r < r2
      iflag_sf_law            = 0 ! compute the SFR per cell to study the SF-law
      iflag_phase             = 0 ! compute the phase of gas (density, temperature, entropy)
      iflag_HaloFind          = 1 ! output halo catalogs
      iflag_halo_properties   = 0 ! compute the halo properties (e.g., Mdm, Mgas, M*, Mcg)
      iflag_halo_five         = 0 ! compute the halo properties (e.g., SFR and profiles)
      iflag_stars_to_hop      = 0 ! stellar particle data for the HOP halo finder
      iflag_virial_radius     = 0 ! copute the virial radius 

c
c.... Compute observables (see ART_Analysis.F)
c
      rmin = 1.0e-3/r0 ! 1h^-1 kpc
      rmax = 10.0/r0   ! 10h^-1 Mpc
      nrbin = 80
      ibintype = 0 
      deltavir = 180.
      deltacrit = 500.
      ifindcenter = 3 
      lpointsout = .false.  
      l_force_center = .false.  
      nstep = 0 
      lpath = index(path    , ' ') - 1
      
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)      

c
c.... 2D projected data
c
c      call DataOutput1 ( tl(0) , dt )  

c
c.... Output 3D data cubes (see Start_Cosmology.F)
      if ( iflag_datacube .eq. 1 ) then
c
c      call DataOutput5( tl(0) , dt , 0.5 , 'L0.5Mpc' ) 
c      call DataOutput5( tl(0) , dt , 1.0 , 'L1Mpc' ) 
c       call DataOutput5( tl(0) , dt , 2.0 , 'L2Mpc' ) 
c       call DataOutput5( tl(0) , dt , 4.0 , 'L4Mpc' ) 
c       call DataOutput5( tl(0) , dt , 8.0 , 'L8Mpc' ) 

c....  3D data cubes in Binary (little-endian)
c       call DataOutput6( tl(0) , dt , 2.0 , 'L2' ) 
      endif


c
c....  Output 2D projected total density (little-endian)
c

      if ( iflag_dataout_2dmap .eq. 1 ) then
        iflag1 = 1  ! total surface density (DM+Gas+Stars)
        iflag2 = 1  ! surface mass density of DM only
        iflag3 = 0  ! X-ray surface brightness
        iflag4 = 0  ! emission weighted temperature
        iflag5 = 0  ! emission weighted metallicity
        iflag6 = 0  ! mass weighted temperature
        iflag7 = 0  ! mass weighted metallicity
        iflag8 = 1  ! thermal SZE (Compton-y)
        iflag9 = 1  ! kinetic SZE (b-parameter)
        iflag10 = 0 ! OVII  561eV [erg/s/cm^2]
        iflag11 = 0 ! OVII  568eV [erg/s/cm^2]
        iflag12 = 0 ! OVII  568eV [erg/s/cm^2]
        iflag13 = 0 ! OVII  665eV [erg/s/cm^2]
        iflag14 = 0 ! OVIII 653eV [erg/s/cm^2]

        if (iflag1 .eq. 1 )  call DataOutput_2Dmap( 1 )
        if (iflag2 .eq. 1 )  call DataOutput_2Dmap( 2 )
        if (iflag3 .eq. 1 )  call DataOutput_2Dmap( 3 )
        if (iflag4 .eq. 1 )  call DataOutput_2Dmap( 4 )
        if (iflag5 .eq. 1 )  call DataOutput_2Dmap( 5 )
        if (iflag6 .eq. 1 )  call DataOutput_2Dmap( 6 )
        if (iflag7 .eq. 1 )  call DataOutput_2Dmap( 7 )
        if (iflag8 .eq. 1 )  call DataOutput_2Dmap( 8 )
        if (iflag9 .eq. 1 )  call DataOutput_2Dmap( 9 )
        if (iflag10 .eq. 1 ) call DataOutput_2Dmap( 10 )
        if (iflag11 .eq. 1 ) call DataOutput_2Dmap( 11 )
        if (iflag12 .eq. 1 ) call DataOutput_2Dmap( 12 )
        if (iflag13 .eq. 1 ) call DataOutput_2Dmap( 13 )
        if (iflag14 .eq. 1 ) call DataOutput_2Dmap( 14 )

      endif

c
c...  Dump 2D slice for viewing with V visualization program
c
      
      ! DumpView centered on the main cluster
c      call ViewDump(MaxLevelNow,xc,1)
c      call ViewDump(MaxLevelNow,yc,2)
c      call ViewDump(MaxLevelNow,zc,3)

      ! location of cold sub-clump in CL6csf
c      xc = 16.592
c      yc = 22.213
c      zc = 128.744

      ! DumpView centered on the cold sub-clump
c      call ViewDump(MaxLevelNow,xc,4)
c      call ViewDump(MaxLevelNow,yc,5)
c      call ViewDump(MaxLevelNow,zc,6)
c      stop 


c
c...  compute X-ray observables (see ART_Analysis_Xray.F )
c

      ! compute the X-ray profiles
      ! idet  = 1   Chandra  [ACIS-s:0.5-7keV]
      !       = 2   BeppoSAX [MECS:2-10keV]
      if ( iflag_pro_xray3d .eq. 1 .or. iflag_pro_xray2d .eq. 1 .or. 
     &     iflag_pro_xray2d_smooth .eq. 1 ) then
        idet = 1 
        call Analyze_Xray ( iflag_pro_xray3d,iflag_pro_xray2d, 
     &        iflag_pro_xray2d_smooth, idet )           
        idet = 2
        call Analyze_Xray ( iflag_pro_xray3d,iflag_pro_xray2d, 
     &        iflag_pro_xray2d_smooth, idet ) 
      endif
      
      ! compute theoretical X-ray spectra in annuli
      if ( iflag_xspectra .eq. 1 ) then
         write(*,*) 'XSpectra'
         rs_dir = '/dworkin/home/daisuke/Raymond/DAT/'
         lfrs = index(rs_dir , ' ') - 1
         rs_spec_fname = rs_dir(1:lfrs)//'rs_spec_0.01-10.01_keV.dat'
         call XSpectra ( rs_spec_fname )
      endif

      ! compute X-ray observables
      if ( iflag_compute_xray .eq. 1 ) then
         write(*,*) 'Compute_Xray'
         call Compute_Xray ( xc, yc, zc, rcrit )
      endif

      ! compute X-ray properties at virial radius
      if ( iflag_xrayvir .eq. 1 ) then
         write(*,*) 'Output_XrayVir'
         call Output_XrayVir ( ifindcenter )
      endif

c
c...  compute SZ observables (see ART_Analysis_SZ.F )
c
      if ( iflag_pro_sz3d .eq. 1 .or. iflag_pro_sz2d .eq. 1 .or. 
     &     iflag_szvir .eq. 1 ) then
         call Analyze_SZ ( iflag_pro_sz3d, iflag_pro_sz2d, iflag_szvir )           
      endif

c
c.... Find the location of the potential minimum
c
      ifindcenter = 3  
      call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      write(iOut,*) 'Finding center : xc, yc, zc =', xc, yc, zc      

c
c...  compute average profiles for simulations with cooling + star formation
c
      if ( iflag_pro_csf .eq. 1 .or. iflag_sfr_rvir .eq. 1 .or. 
     &     iflag_sfr_r1_r2 .eq. 1 ) then 

c         xc = 17.8576
c         yc = 23.9472
c         zc = 125.586
c         l_force_center = .true.

         fname  = path(1:lpath)//'/pro'//'_csf_'//'a'//
     &        digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &        //'.dat '
         write(*,*) 'Output_Average_Profiles_CSF'
         call Output_Average_Profiles_CSF ( xc , yc , zc , rvir, rcrit,  
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, deltacrit, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fname ) 
       endif

c
c...  3D profiles for fitting
c
      if ( iflag_pro_fit .eq. 1 ) then 

         fname  = path(1:lpath)//'/pro'//'_fit_'//'a'//
     &        digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &        //'.dat '

         write(*,*) 'Output_Fit_Profiles'
         call Output_Fit_Profiles ( xc , yc , zc , rvir, rcrit,  
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, deltacrit, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fname ) 
      endif

c
c...  profiles of HSE mass and mach number
c
      if ( iflag_pro_HSE .eq. 1 ) then 

         rmin = 3.0e-3/r0 ! 3h^-1 kpc
         rmax = 10.0/r0   ! 10h^-1 Mpc
         vcut = 300.      ! 300km/s
         nrbin = 80
         deltavir = 180.
         deltacrit = 180.

         cldir  = 'CL6'
         clname = 'CL6a'
         nfn1 = index ( cldir , ' ' ) - 1 
         nfn2 = index ( clname , ' ' ) - 1
       
         hname = '/avalon/data/daisuke/Hydro/'//
     &     cldir(1:nfn1)//'/'//clname(1:nfn2)//'/HF/h_blist_a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'.dat '

         fname  = path(1:lpath)//'/pro'//'_HSE_'//'a'//
     &        digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &        //'.dat '

         write(*,*) 'Output_HSE_Profiles'
         call Output_HSE_Profiles ( xc , yc , zc , rmin , rmax , vcut, 
     &        deltavir, deltacrit, nrbin, l_force_center, fname, hname )
      endif      

c
c...  output data file for gas velocities
c
      if ( iflag_gas_velocities .eq. 1 ) then
         fname  = path(1:lpath)//'/vgas_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

         write(*,*) ' '
         write(*,*) 'Output datafile for gas velocities'
         call Output_Gas_Velocities ( fname )
      endif

c
c...  compute SFR for stellar particles within < Rvir
c
      if ( iflag_sfr_rvir .eq. 1 ) then
         fname  = path(1:lpath)//'/sfr_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

         amin = 0.06
         amax = 1.0
         nabin = 180

         write(*,*) ' '
         write(*,*) 'Computing SFR for stellar particles within rvir'
         call SFR_rvir ( xc, yc, zc, rvir,
     &        amin , amax , nabin , ibintype ,  
     &        deltavir, fname )
      endif

c
c...  compute SFR for stellar particles within r1 < r < r2
c
      if ( iflag_sfr_r1_r2 .eq. 1 ) then
         fname  = path(1:lpath)//'/sfr_0to60kpc_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

         r1 = 0.
         r2 = 60*1.e-3/r0  ! 60/h kpc comiving
         amin = 0.06
         amax = 1.0
         nabin = 180

         write(*,*) ' '
         write(*,*) 'Computing SFR for stellar particles within r1<r<r2'
         call SFR_r1_r2 ( xc, yc, zc, r1, r2, 
     &        amin , amax , nabin , ibintype ,  
     &        deltavir, fname )

         fname  = path(1:lpath)//'/sfr_60kpctorvir_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

         r1 = 60*1.e-3/r0  ! 60/h kpc comiving
         r2 = rvir         ! rvir

         write(*,*) ' '
         write(*,*) 'Computing SFR for stellar particles within r1<r<r2'
         call SFR_r1_r2 ( xc, yc, zc, r1, r2, 
     &        amin , amax , nabin , ibintype ,  
     &        deltavir, fname )

      endif

c
c...  study the SF-law 
c
      if ( iflag_sf_law .eq. 1 ) then
         fname  = path(1:lpath)//'/sf_law_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
         ! study the SF law on the Level [aiming for < 50/h kpc proper]
         if ( aexpn .le. 0.2 ) then 
            lev = 1
         else if ( aexpn .le. 0.5) then 
            lev = 2 
         else if ( aexpn .le. 0.7) then 
            lev = 3
         else 
            lev = 4
         endif
         tsfave = 5.0e8*aexpn   ! averaging over this many years
         write(*,*) ' '
         write(*,*) 'Creating data for studying the SF-law'
          call AverageSF ( lev , tsfave , fname)
      endif

c
c...  compute the phase of gas within < 10*rvir
c
      if ( iflag_phase .eq. 1 ) then
         fname  = path(1:lpath)//'/phase_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
         write(*,*) ' '
         write(*,*) 'Computing the phase of gas '
         delta = 0.0
         call Compute_Phase ( delta, fname )
c         call Compute_Phase ( xc, yc, zc, rmax, ifindcenter, 
c     &        l_force_center, fname )
      endif

c
c...  output halo catalog
c
      if ( iflag_HaloFind .eq. 1 ) then
         write(*,*) 'Output halo catalog'
         call HaloFind ( 3 )
      endif

c
c...  compute halo properties
c
      if ( iflag_halo_properties .eq. 1 ) then
         fname  = '/avalon/data/daisuke/Hydro/CL6/CL6csf1_sf2/HF'//
     &               '/hlist'//'_a'//
     &               digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &               //'.dat'
         write(*,*) ' '
         write(*,*) 'Compute Halo Properties'
         call Halo_Properties ( fname )
      endif

c
c...  compute SFR for stellar particles within < Rvir and 
c     radial profiles for the 5 most massive halos
c
      if ( iflag_halo_five .eq. 1 ) then
         fname  = '/avalon/data/daisuke/Hydro/CL6/CL6csf1_sf2_AGN/'//
     &               'HF/hlist'//'_'//
     &               digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &               //'.dat'
         write(*,*) ' '
         write(*,*) 'Compute Halo Properties of 5 most massive halos'
         call Halo_Properties_Massive ( fname )
      endif

c
c...  output the stellar particles for the HOP halo finder
c
      if ( iflag_stars_to_hop .eq. 1 ) then
         fname  = path(1:lpath)//'/st_'//'a'//
     &              digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &              //'_hop.in'
         write(*,*) ' '
         write(*,*) 'Output the input file for the HOP halo finder'
         call Output_Stars_to_HOP ( fname )

         fname  = path(1:lpath)//'/st_'//'a'//
     &              digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &              //'_hop.dat'
         write(*,*) ' '
         write(*,*) 'Output the stellar particles for HOP'
         call Write_Stars_ASCII2 ( fname, xc, yc, zc ) 
      endif

c
c...  compute virial radius
c
      if ( iflag_virial_radius .eq. 1 ) then
        ifindcenter = 3
        call Find_Halo_Center ( xc, yc, zc, ifindcenter )

        call Find_Halo_Radius ( xc, yc, zc, 180., 180.,
     &       i180, r180, aM180, aM180_dm, aM180_gas, aM180_st, 
     &       i180c, r180c, aM180c, aM180c_dm, aM180c_gas, aM180c_st, 
     &       ifindcenter )
  
        call Find_Halo_Radius ( xc, yc, zc, 200., 200., 
     &       i200, r200, aM200, aM200_dm, aM200_gas, aM200_st, 
     &       i200c, r200c, aM200c, aM200c_dm, aM200c_gas, aM200c_st, 
     &       ifindcenter )
  
        call Find_Halo_Radius ( xc, yc, zc, 500., 500., 
     &       i500, r500, aM500, aM500_dm, aM500_gas, aM500_st, 
     &       i500c, r500c, aM500c, aM500c_dm, aM500c_gas, aM500c_st, 
     &       ifindcenter )
  
        call Find_Halo_Radius ( xc, yc, zc, 2500., 2500., 
     &       i2500, r2500, aM2500, aM2500_dm, aM2500_gas, aM2500_st, 
     &       i2500c,r2500c,aM2500c, aM2500c_dm, aM2500c_gas, aM2500c_st, 
     &       ifindcenter )

        write(*,*)
        write(*,20) xc, yc, zc
        write(*,21) 
        write(*,23) 180., r180*r0, aM180, aM180_dm, aM180_gas, aM180_st
        write(*,23) 200., r200*r0, aM200, aM200_dm, aM200_gas, aM200_st
        write(*,23) 500., r500*r0, aM500, aM500_dm, aM500_gas, aM500_st
        write(*,23) 2500., r2500*r0, aM2500, aM2500_dm, aM2500_gas,
     &       aM2500_st

        write(*,*)
        write(*,22)
        write(*,23) 180., r180c*r0, aM180c, aM180c_dm, aM180c_gas, 
     &       aM180c_st
        write(*,23) 200., r200c*r0, aM200c, aM200c_dm, aM200c_gas, 
     &       aM200c_st
        write(*,23) 500., r500c*r0, aM500c, aM500c_dm, aM500c_gas, 
     &       aM500c_st
        write(*,23) 2500., r2500c*r0, aM2500c, aM2500c_dm, aM2500c_gas, 
     &       aM2500c_st

      endif

 20   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 21   format ('# Dvir(mean)  Rvir[1/hMpc]      Mvir       Mdm  ',
     &     '      Mgas        Mst  [1/h Msun]' )
 22   format ('# Dcrit(crit) Rcrit[1/h Mpc]    Mcrit      Mdm  ',
     &     '      Mgas        Mst  [1/h Msun]' )
 23   format (f9.2,6x,f7.3,1x,4(g11.4,1x))


      STOP
      END

c     -----------------------------
      subroutine Initialize ( dt )
c     -----------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
# include "a_def.h"
c
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
c
      do i = 1 , 10
        CPU(i) = 0.
      enddo
c
      return
      end      
