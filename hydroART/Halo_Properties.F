c     Halo_Properties ()
c     Halo_Properties_Massive ()
c     Output_GC_Cells ()
c     Identify_GCs ()
c     Estimate_GC_Mass ()
c     Output_HighDensity ()
c     DensityCenter_Profiles_HP()
c     DensityCenter_Profiles_HP2()


c     --------------------------------------
      subroutine Halo_Properties ( hc_name )
c     --------------------------------------
c
c     Output properties of all halos.
c
c     output : h_blist_aX.XXX.dat
c              h_bpro_aX.XXX.dat
c              h_bvpro_aX.XXX.dat
c
#     include "a_def.h"
      include 'a_analysis.h'      
      include 'a_hfind.h'
c
      character*5 fstep
      real*8 sum1, sum2, sum3, sum4
      character*256 fname_gas, hc_name, textline
      logical lpointsout, l_force_center
      integer imovstep
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*2 chbox
      real*8 pfact, rhofact 
      parameter ( nctrmax = 1000000 )
      integer ictr(nctrmax), indc(nctrmax)
      real pctr(nctrmax), dctr(nctrmax)
      common / CTR1 / ictr, indc
      common / CTR2 / pctr, dctr
c
#ifdef PRESSUREFLOOR
      real*8 f_J, f_eJ
      parameter ( f_J = 10.  ) ! pressure floor parameter
#endif

#ifdef PRESSUREFLOOR
      DO Level = MinL_Jeans , MaxLevelNow
c
c.... artificial pressure floor a-la Machacek et al. 2
        if ( Level .ge. MinL_Jeans ) then
          f_eJ = 0.47746 * f_J * aexpn * 
     &           CellSize(MaxLevel)**2 * (gamma - 1.0)
        else
          f_eJ = 0.0
        endif
        if ( Level .eq. MinLevel ) then
          do i1 = 1 , ncell0
            E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
            E_th = f_eJ * hvar(1,i1) * hvar(1,i1) 
            hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
          enddo 
        else
          call Select_Cells ( Level , nLevel )
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

c
c.... output profiles 
c

      rmin = 5.e-3/r0  ! min. radius for the profile (5/h kpc comoving) in code units
      rmax = 2.0/r0    ! max. radius for the profile (2000/h kpc comoving) in code units
      nrbin = 45 
      ibintype = 0 
      deltavir = 340
      deltacrit = 500
      ifindcenter = 3  ! 0/1/2/3 = gas density, DM density, Total density, potential min.
      lpointsout = .false.  
      l_force_center = .false.  

      lpath = index(path    , ' ') - 1

      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

      rhofact = rho0 * 1.d-18 / aexpn**3
      tfact = (gamma-1.0) * T_0 / aexpn**2 
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      rfact = r0 * 1.e3 * aexpn / hubble ! /h Mpc -> proper kpc
      vfact = v0 / aexpn  ! code velocity -> peculiar in km/s
      amfact = aM0 * hubble ! code mass -> Msun
      rmin = 3.e-2 / rfact 
      rmax = 1.5 / rfact 
      nrbin = 30 

      Tcold = 2.e5

c
c.... read DM halos 
c

      nfn = index ( hc_name , ' ' ) - 1
      open ( 20 , file = hc_name(1:nfn) ) 
c
c.... read header
c
      do i = 1 , 17 
        read(20,*) textline 
      enddo
c
c.... now read the actual halo list 
c
      nh = 0 
      do i = 1 , nhmax 
c       Reading outputs from the Andrey's HF
c        read(20,*,end=20) id1, id2, xhd, yhd, zhd, vxhd, vyhd, vzhd,
c     &     rhvd, amhd, npd, vmd, rvmd, rsd
        read(20,*,end=20) id1, id2, xhd, yhd, zhd, vxhd, vyhd, vzhd,
     &     rhd, rhvd, amhd, amvhd, npd, vmd, rvmd, rsd
        nh = nh + 1
        if ( nh .gt. nhmax ) then 
          write(*,*) '* error in Halo_Properties: nh>nhmax:',
     &    nh, nhmax
          write(*,*) 'increase nhmax and rerun. stopping...'
          stop
        endif
        ih1(nh) = id1
        ih2(nh) = id2
        xh(nh) = xhd
        yh(nh) = yhd
        zh(nh) = zhd
        vxh(nh) = vxhd
        vyh(nh) = vyhd
        vzh(nh) = vzhd
        rh(nh) = rhd
        rhvir(nh) = rhvd
        amh(nh) = amhd
        amhvir(nh) = amvhd
        nhp(nh) = npd 
        vhmax(nh) = vmd
        rhmax(nh) = rvmd
        rsh(nh) = rsd
      enddo
 20   close ( 20 )

c      
c.... now profiles around centers of DM halos 
c
      dout = 180
c      dout = 1.0 / rhofact ! outer density (total: baryons+DM) in Msun/pc^-3      
c      rmin = 5.e-3/r0  ! min. radius for the profile (5/h kpc comoving) in code units
c      rmax = 2.0/r0    ! max. radius for the profile (2000/h kpc comoving) in code units
      rmin = 5.0 / rfact    ! min. radius for the profile (5/h kpc proper)
      rmax = 2000.0 / rfact ! max. radius for the profile (2000/h kpc proper)
      nrbin = 108 
      Tcold = 2.e5

      fname_gas = path(1:lpath)//'/h_blist_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//
     &       '.dat '
      nfn = index( fname_gas , ' ' ) - 1 
      open ( 11 , file = fname_gas(1:nfn) ) 
      write(11,150) 
 150   format('# Halo list with properties determined from profiles')
      write(11,151) dout, Tcold
 151  format('# density (DM+bar.) defining the outer radius=',e9.3,
     &'[Msun pc^-3];  Tcold=',e9.3)
      write(11,152) rmin*rfact, rmax*rfact, nrbin 
 152  format('# rmin,max=',e9.3,1x,e9.3,' [kpc] (proper);'
     &  '  nrbin=',i4)
      write(11,156)
 156  format('#   x        y         z   min(rout,rmax)  rvir    Mg   
     & Mcg    M*     Mdm     <Zg>    <Z*>    <t*>    Mhalo ')
      write(11,158)
 158  format('# coord. in /h Mpc (comoving)      kpc (proper) 
     &all masses in /h Msun')

c
      fname_gas = path(1:lpath)//'/h_bpro_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//
     &       '.dat '
      nfn = index( fname_gas , ' ' ) - 1 
      open ( 12 , file = fname_gas(1:nfn) )
 
      fname_gas = path(1:lpath)//'/h_bvpro_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//
     &       '.dat '
      nfn = index( fname_gas , ' ' ) - 1 
      open ( 13 , file = fname_gas(1:nfn) ) 

c
      do i = 1 , nh 
        if ( amh(i) .gt. 2.e7 ) then 
        write(*,*) 'processing halo ',i
        ihc = i ! remember the index of current halo 
        Posx = xh(i) / r0 + 1.0 
        Posy = yh(i) / r0 + 1.0 
        Posz = zh(i) / r0 + 1.0 
        rmass = rh(i) * 1.e-3 / r0     ! comoving /h kpc -> code units
        rvir = rhvir(i) * 1.e-3 / r0   ! comoving /h kpc -> code units
        call DensityCenter_Profiles_HP 
     &    ( dc, dout, Tcold,  Posx, Posy, Posz,
     &      rmin , rmax , nrbin , rmass, rvir, 1 )
        endif
      enddo
      close ( 11 ) 
      close ( 12 ) 
      close ( 13 ) 

 220  return
      end

c     ----------------------------------------------
      subroutine Halo_Properties_Massive ( hc_name )
c     ----------------------------------------------
c
c     Output radial profiles and starformation rate within rvir
c     of the 5 most massive halos.
c
#     include "a_def.h"
      include 'a_analysis.h'      
      include 'a_hfind.h'
c
      character*5 fstep
      real*8 sum1, sum2, sum3, sum4
      character*256 fname, hc_name, textline
      logical lpointsout, l_force_center
      integer imovstep
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*2 chbox
      real*8 pfact, rhofact 
      parameter ( nctrmax = 1000000 )
      integer ictr(nctrmax), indc(nctrmax)
      real pctr(nctrmax), dctr(nctrmax)
      common / CTR1 / ictr, indc
      common / CTR2 / pctr, dctr
c
#ifdef PRESSUREFLOOR
      real*8 f_J, f_eJ
      parameter ( f_J = 10.  ) ! pressure floor parameter
#endif

#ifdef PRESSUREFLOOR
      DO Level = MinL_Jeans , MaxLevelNow
c
c.... artificial pressure floor a-la Machacek et al. 2
        if ( Level .ge. MinL_Jeans ) then
          f_eJ = 0.47746 * f_J * aexpn * 
     &           CellSize(MaxLevel)**2 * (gamma - 1.0)
        else
          f_eJ = 0.0
        endif
        if ( Level .eq. MinLevel ) then
          do i1 = 1 , ncell0
            E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
            E_th = f_eJ * hvar(1,i1) * hvar(1,i1) 
            hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
          enddo 
        else
          call Select_Cells ( Level , nLevel )
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

c
c.... output profiles 
c

      rmin = 1.e-3/r0   ! min. radius for the profile (5/h kpc comoving) in code units
      rmax = 10.0/r0    ! max. radius for the profile (2000/h kpc comoving) in code units
      nrbin = 80 
      ibintype = 0 
      deltavir = 340
      deltacrit = 500
      ifindcenter = 3  ! 0/1/2/3 = gas density, DM density, Total density, potential min.
      lpointsout = .false.  
      l_force_center = .true.  

      lpath = index(path    , ' ') - 1

      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

      rhofact = rho0 * 1.d-18 / aexpn**3
      tfact = (gamma-1.0) * T_0 / aexpn**2 
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      rfact = r0 * 1.e3 * aexpn / hubble ! /h Mpc -> proper kpc
      vfact = v0 / aexpn  ! code velocity -> peculiar in km/s
      amfact = aM0 * hubble ! code mass -> Msun
  
c...  parameters for the SFR calculation
      amin = 0.06
      amax = 1.0
      nabin = 180
        
c
c.... read DM halos 
c
      nfn = index ( hc_name , ' ' ) - 1
      open ( 20 , file = hc_name(1:nfn) ) 
c
c.... read header
c
      do i = 1 , 17 
        read(20,*) textline 
      enddo
c
c.... now read the actual halo list 
c
      nh = 0 
      do i = 1 , nhmax 
c       Reading outputs from the Andrey's HF
c        read(20,*,end=20) id1, id2, xhd, yhd, zhd, vxhd, vyhd, vzhd,
c     &     rhvd, amhd, npd, vmd, rvmd, rsd
        read(20,*,end=20) id1, id2, xhd, yhd, zhd, vxhd, vyhd, vzhd,
     &     rhd, rhvd, amhd, amvhd, npd, vmd, rvmd, rsd
        nh = nh + 1
        if ( nh .gt. nhmax ) then 
          write(*,*) '* error in Halo_Properties: nh>nhmax:',
     &    nh, nhmax
          write(*,*) 'increase nhmax and rerun. stopping...'
          stop
        endif
        
        xc = xhd / r0 + 1.0 
        yc = yhd / r0 + 1.0 
        zc = zhd / r0 + 1.0 
        rout = rhd / rfact ! proper kpc -> grid

        ! halo IDs of the five most massive halos
        ! CL6csf1_a1.001 (1,2,5,6,11)
        ! CL6csf1_sf2_a0.333 (1,2,3,6,5)
        ! CL6csf1_sf2_AGN_a0.334 (1,3,2,6,5)
        ihalo1 = 1
        ihalo2 = 2
        ihalo3 = 3
        ihalo4 = 6
        ihalo5 = 5

        if ( id1 .eq. ihalo1 .or. id1 .eq. ihalo2 .or. 
     &       id1 .eq. ihalo3 .or. id1 .eq. ihalo4 .or. 
     &       id1 .eq. ihalo5  ) then
          if ( id1 == ihalo1 )  icls = 1
          if ( id1 == ihalo2 )  icls = 2
          if ( id1 == ihalo3 )  icls = 3
          if ( id1 == ihalo4 )  icls = 4
          if ( id1 == ihalo5 )  icls = 5
          fname  = path(1:lpath)//'/pro'//'_csf_h'//digits(icls)//'_a'
     &         //digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &         //'.dat '
          
          write(*,*) 'Output_Average_Profiles_CSF '
          call Output_Average_Profiles_CSF ( xc , yc , zc , rvir, rcrit,  
     &                           rmin , rmax , nrbin , ibintype ,  
     &                           deltavir, deltacrit, ifindcenter , 
     &                           lpointsout, l_force_center,
     &                           fname ) 
          write(*,*) 'Done writing :',fname
  
          fname  = path(1:lpath)//'/sfr_h'//digits(icls)//'_a'
     &         //digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &         //'.dat '
          
          write(*,*) ' '
          write(*,*) 'Computing SFR for stellar particles within rout'
          call SFR_rvir ( xc, yc, zc, rout,
     &         amin , amax , nabin , ibintype ,  
     &         deltavir, fname )
          write(*,*) 'Done writing :',fname
        endif

      enddo

 20   close ( 20 )

 220  return
      end

c
c     -------------------------------------
      subroutine Output_GC_Cells ( rhoth )
c     -------------------------------------
c
c     rhoth - minimum threshold density to use in Msun/pc^3
c      
c
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'

      character*5 fstep
      real*8 sum1, sum2, sum3, sum4
      character*256 fname_gas
      logical lpointsout, l_force_center
      integer imovstep
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*2 chbox
      real*8 pfact, rhofact 
      parameter ( nctrmax = 1000000 )
      integer ictr(nctrmax), indc(nctrmax)
      real pctr(nctrmax), dctr(nctrmax)
      common / CTR1 / ictr, indc
      common / CTR2 / pctr, dctr

#ifdef PRESSUREFLOOR
      real*8 f_J, f_eJ
      parameter ( f_J = 10.  ) ! pressure floor parameter
#endif

#ifdef PRESSUREFLOOR
      DO Level = MinL_Jeans , MaxLevelNow
c
c.... artificial pressure floor a-la Machacek et al. 2
        if ( Level .ge. MinL_Jeans ) then
          f_eJ = 0.47746 * f_J * aexpn * 
     &           CellSize(MaxLevel)**2 * (gamma - 1.0)
        else
          f_eJ = 0.0
        endif
        if ( Level .eq. MinLevel ) then
          do i1 = 1 , ncell0
            E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
            E_th = f_eJ * hvar(1,i1) * hvar(1,i1) 
            hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
          enddo 
        else
          call Select_Cells ( Level , nLevel )
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

      do ic = 1 , nctot
        ind(ic) = 0
      enddo
c
c.... output profiles 
c
      xc = 32.95605469
      yc = 13.58300781 
      zc = 39.06152344
      rmin = 5.e-3/r0  ! min. radius for the profile
      rmax = 2.0/r0    ! max. radius for the profile

      nrbin = 45 
      ibintype = 0 
      deltavir = 340
      deltacrit = 500
      ifindcenter = 3 
      lpointsout = .false.  
      l_force_center = .false.  

      lpath = index(path    , ' ') - 1

      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

      fname_gas  = path(1:lpath)//'/pro'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

      call Output_Average_Profiles_Stars ( xc , yc , zc , rvir,  
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, deltacrit, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fname_gas ) 
c
c.... identify GC cells
c
      rhofact = rho0 * 1.d-18 / aexpn**3
      tfact = (gamma-1.0) * T_0 / aexpn**2 
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      rfact = r0 * 1.e6 * aexpn / hubble ! /h Mpc -> proper pc
      vfact = v0 / aexpn  ! code velocity -> peculiar in km/s
      amfact = aM0 * hubble ! code mass -> /h Msun
c
c.... efficiency of SF in proto-GC sites
c     be careful not too confuse with eps_SF in th rest of code
      e_SF = 0.5 
c.... crit. density for catastrophic s.f. in Msun/pc^3
      d_SF = 5.e3 
c.... min. GC mass to consider (in Msun) 
      gcm_min = 1.e4

      fname_gas = 'gcc__cloud_a'//
     &    digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//'.dat '
      nfn = index( fname_gas , ' ' ) - 1 
      open ( 12 , file = fname_gas(1:nfn) ) 

      fname_gas = 'ic_cloud_a'//
     &    digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//'.dat '
      nfn = index( fname_gas , ' ' ) - 1 
      open ( 13 , file = fname_gas(1:nfn) ) 

      write(12,18) 
 18   format('# properties of cells - candidate GC sites')
      write(12,21) boxh, Om0, Oml0, Omb0, hubble      
 21   format('# Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f5.3,1x,
     &        ' Oml0= ',f5.3,1x,' Omb0= ',f7.5,1x,' h = ',f5.3)
      write(12,20) istep, t, dt, aexpn, ainit
 20   format ('# istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      write(12,22) e_SF, d_SF, gcm_min
 22   format ('# e_SF =',e8.3,'; d_SF =',g15.8,' [Msun/pc^3]; ',
     &        'M_GC_min =',e8.3,' [Msun]')
      write(12,23)
      write(12,24)
c      write(12,25)
 23   format ('#rcell      M_GC      r_GC rho_cell Tcell   P_cell 
     &M_cell   Z_II     Z_Ia         x         y         z  
     &vx     vy     vz')
 24   format ('# pc       Msun        pc  Ms/pc^3    K     K cm^-3 
     &Msun    Zsun     Zsun                h^-1 Mpc                    
     &km/s')
 25   format ('#')

      write(13,19) 
 19   format('# indices of cells - candidate GC sites')
      write(13,21) boxh, Om0, Oml0, Omb0, hubble      
      write(13,20) istep, t, dt, aexpn, ainit
      write(13,22) e_SF, d_SF, gcm_min
      write(13,26)
      write(13,27)
c      write(13,25)
 26   format('# ic      L  rcell  rhocell  Pcell   Mcell')
 27   format('#             pc   Ms/pc^3  K cm^-3   Msun')
c
      nctr = 0 
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

                Posx = Posx + corr_x
                Posy = Posy + corr_y
                Posz = Posz + corr_z

c                rr = sqrt( (Posx-xc)**2 +
c     &                     (Posy-yc)**2 +
c     &                     (Posz-zc)**2 )
c                if ( rr .lt. rmax ) then
                  rhocell = hvar(1,ic) * rhofact 
                  if ( rhocell .gt. rhoth ) then 
                    rcell = 0.5 * CellSize(Level) * rfact 
c....               GC mass and radius      
                    rsf = sqrt(0.333333*rhocell/d_SF) * rcell
                    gcm = pi4 * e_SF * d_SF * rsf**3 
c                    if ( gcm .gt. gcm_min ) then 
                      gcr = rsf / e_SF 
                      nctr = nctr + 1
                      if ( nctr .gt. nctrmax ) then
                        write(*,*) 'nctr exceeds nctrmax =',nctrmax
                        write(*,*) 'increase nctrmax and rerun'
                        stop
                      endif
                      pcell = hvar(8,ic) * pfact ! pressure 
                      rhor = 1.d0 / hvar(1,ic)
                      dtemp = tfact * hvar(8,ic) * rhor
c....                 metallicity in units of solar 
                      dzII = 50.* hvar(izII,ic) * rhor
                      dzIa = 50.* hvar(izIa,ic) * rhor
                      dvcx = hvar(3,ic) * rhor * vfact 
                      dvcy = hvar(4,ic) * rhor * vfact
                      dvcz = hvar(5,ic) * rhor * vfact
                      dmass = hvar(1,ic)* CellVolume * amfact 
                      xcell = (Posx - xc) * r0
                      ycell = (Posy - yc) * r0
                      zcell = (Posz - zc) * r0
                      
c                      write(*,80) ic, Level, rcell, gcm, gcr, rsf, 
c     &                  rhocell, pcell, dmass, 
c     &                  dzII, dzIa, dvcx, dvcy, dvcz
                      write(13,79) ic, Level, rcell, 
     &                  rhocell, pcell, dmass
                      write(12,80) rcell, gcm, gcr, 
     &                  rhocell, dtemp, pcell, dmass, 
     &                  dzII, dzIa, xcell, ycell, zcell, 
     &                  dvcx, dvcy, dvcz
 79                   format(i8,1x,i2,1x,f6.1,1x,f6.1,2x,
     &                  2(e8.3,1x))
 80                   format(f6.1,2x,e11.5,1x,f6.1,1x,
     &                  f6.1,1x,5(e8.3,1x),1x,3(f9.6,1x),1x,3(f6.1,1x))
                      ictr(nctr) = ic
                      pctr(nctr) = pcell
                      dctr(nctr) = rhocell
c                    endif
                  endif
c                endif
              endif
            enddo
          enddo
      ENDDO
      close ( 12 ) 

c      return
c
      call indexx ( nctr , dctr , indc )

      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

      fname_gas = 'dpro_cloud_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//' '
      nfn = index( fname_gas , ' ' ) - 1 

      open ( 12 , file = fname_gas(1:nfn) ) 
      fname_gas = 'dc_cloud_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//' '
      nfn = index( fname_gas , ' ' ) - 1 

      open ( 11 , file = fname_gas(1:nfn) ) 
      fname_gas = 'vpro_cloud_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//' '
      nfn = index( fname_gas , ' ' ) - 1 

      open ( 13 , file = fname_gas(1:nfn) ) 

      dc = 10.0 / rhofact 
      dout = 1.0 / rhofact
      rmin = 30.0 / 1.e6 / r0 / aexpn   ! /h pc-> code units
      rmax = 3000.0 / 1.e6 / r0 / aexpn ! /h pc-> code units
      nrbin = 30 
      Pout = 1.e5 ! outer pressure in cm^-3 K
      Tout = 350. ! outer temperature in K
      Tcold = 1.e4

      do ic = nctr , 1 , -1 
        icell = ictr(indc(ic))
        write(*,*) 'processing center',ic
        if ( ind(icell) .eq. 0 ) then 
          call Ps ( icell , Posx , Posy , Posz )
          call DensityCenter_Profiles_HP2 
     &         ( Pout, Tout, Tcold,  Posx, Posy, Posz,
     &           rmin , rmax , nrbin , rmax, 0 )
        endif
      enddo
c
      close ( 11 ) 
      close ( 12 ) 
      close ( 13 ) 
c      
      return
      end
c
c     -----------------------------------------
      subroutine Identify_GCs ( rhoth , p_th, prefix, suffix  )
c     ----------------------------------------
c
c     rhoth - minimum threshold density to use in Msun/pc^3
c     p_th - pressure threshold for cloud identification
c      
c
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'

      character*256 prefix, suffix
      character*5 fstep
      real*8 sum1, sum2, sum3, sum4
      character*256 fname_gas
      logical lpointsout, l_force_center
      integer imovstep
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*2 chbox
      real*8 pfact, rhofact 
      parameter ( nctrmax = 1000000 )
      integer ictr(nctrmax), indc(nctrmax)
      real pctr(nctrmax), dctr(nctrmax)
      common / CTR1 / ictr, indc
      common / CTR2 / pctr, dctr
      parameter ( ncloud = 10000 )  ! max. # of clouds
      integer icmax(ncloud), ncc(ncloud)
      real*8 dcmax(ncloud), amc(ncloud), volc(ncloud), tempc(ncloud)
      real*8 zIac(ncloud), zIIc(ncloud), prc(ncloud), rhoc(ncloud)
      
      integer iNbb(6), iNC(27), nnc

#ifdef PRESSUREFLOOR
      real*8 f_J, f_eJ
      parameter ( f_J = 10.  ) ! pressure floor parameter
#endif

#ifdef PRESSUREFLOOR
      DO Level = MinL_Jeans , MaxLevelNow
c
c.... artificial pressure floor a-la Machacek et al. 2
        if ( Level .ge. MinL_Jeans ) then
          f_eJ = 0.47746 * f_J * aexpn * 
     &           CellSize(MaxLevel)**2 * (gamma - 1.0)
        else
          f_eJ = 0.0
        endif
        if ( Level .eq. MinLevel ) then
          do i1 = 1 , ncell0
            E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
            E_th = f_eJ * hvar(1,i1) * hvar(1,i1) 
            hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
          enddo 
        else
          call Select_Cells ( Level , nLevel )
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

      do ic = 1 , nctot
        ind(ic) = 0
      enddo
c
c.... identify the center of the main object
c
      xc = 32.95605469
      yc = 13.58300781 
      zc = 39.06152344
      rmin = 5.e-3/r0  ! min. radius for the profile
      rmax = 2.0/r0    ! max. radius for the profile

      nrbin = 45 
      ibintype = 0 
      deltavir = 340
      deltacrit = 500
      ifindcenter = 3 
      lpointsout = .false.  
      l_force_center = .false.  

      lpath = index(path    , ' ') - 1

      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

      fname_gas  = path(1:lpath)//'/pro'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

      call Output_Average_Profiles_Stars ( xc , yc , zc , rvir,  
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, deltacrit, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fname_gas ) 
c
c.... identify GC cells
c
      rhofact = rho0 * 1.d-18 / aexpn**3
      tfact = (gamma-1.0) * T_0 / aexpn**2 
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      rfact = r0 * 1.e6 * aexpn / hubble ! /h Mpc -> proper pc
      vfact = v0 / aexpn  ! code velocity -> peculiar in km/s
      amfact = aM0 * hubble ! code mass -> /h Msun
c
c.... now identify high-density clouds 
c
      nctr = 0 
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

                Posx = Posx + corr_x
                Posy = Posy + corr_y
                Posz = Posz + corr_z

c                rr = sqrt( (Posx-xc)**2 +
c     &                     (Posy-yc)**2 +
c     &                     (Posz-zc)**2 )
c                if ( rr .lt. rmax ) then
                  rhocell = hvar(1,ic) * rhofact 
                  pcell = hvar(8,ic) * pfact ! pressure 
c                  if ( rhocell .gt. rhoth .and. pcell .gt. p_th ) then 
c                  if ( pcell .gt. p_th ) then 
                  if ( rhocell .gt. rhoth ) then 
                    nctr = nctr + 1
                    if ( nctr .gt. nctrmax ) then
                      write(*,*) 'nctr exceeds nctrmax =',nctrmax
                      write(*,*) 'increase nctrmax and rerun'
                      stop
                    endif
                    ictr(nctr) = ic
                    pctr(nctr) = pcell
                    dctr(nctr) = rhocell
                  endif
              endif
            enddo
          enddo
      ENDDO
c
c.... sort selected cells by density 
c
      call indexx ( nctr , dctr , indc )
c
c.... mark cell clusters...
c     start with highest density cell, mark neighbors with the cluster # 
c     and so on. 
c
      ncl = 0 ! cluster counter
c
      do ic = nctr , 1 , -1 
        icell = ictr(indc(ic))
        write(*,*) 'processing center',ic
        if ( ind(icell) .eq. 0 ) then ! if not marked yet, start a new cluster
          ncl = ncl + 1       
          if ( ncl .gt. ncloud ) then 
            write(*,*) 
     &  '* in Identify_GC ncl exceeded max. allowed number ncloud...'
            write(*,*) 'increase ncloud and rerun. stopping...'
            stop
          endif
          icloud = ncl
          ind(icell) = icloud 
        else ! if marked mark neighbors with the same cluster number
          icloud = ind(icell) 
        endif
        Level = iLv(icell) 
        Size = CellSize(Level)
        
        call Find_Neighbors_and_Corners ( icell , Level , 
     &                                    Size  , iNC   , nnc )
        do ib = 1 , nnc
          if ( ind(iNC(ib)) .eq. 0 ) then 
            rhonb = hvar(1,iNC(ib)) * rhofact 
            pnb   = hvar(8,iNC(ib)) * pfact 
c            if ( rhonb .gt. rhoth .and. pnb .gt. p_th ) then 
c            if ( pnb .gt. p_th ) then 
            if ( rhonb .gt. rhoth ) then 
              ind(iNC(ib)) = icloud 
            endif
          endif
        enddo
c        call iNbAll ( icell , iNbb ) ! get neighbor indices 
c        do ib = 1 , 6
c          if ( ind(iNbb(ib)) .eq. 0 ) then 
c            rhonb = hvar(1,iNbb(ib)) * rhofact 
c            pnb   = hvar(8,iNbb(ib)) * pfact 
cc            if ( rhonb .gt. rhoth .and. pnb .gt. p_th ) then 
cc            if ( pnb .gt. p_th ) then 
c            if ( rhonb .gt. rhoth ) then 
c              ind(iNbb(ib)) = icloud 
c            endif
c          endif
c        enddo
      enddo
c
      write(*,*) 'identified ncl=',ncl,' clouds...'
c
c.... now analyze clouds
c
c....first determine the highest density cells of each cloud and its
c    mass, volume
c
      do i = 1 , ncloud 
        icmax(i) = 0 
        ncc(i) = 0
        dcmax(i) = -1.e6
        amc(i)   = 0.
        volc(i)  = 0.
        tempc(i) = 0.
        zIac(i)  = 0.
        zIIc(i)  = 0.
        prc(i)   = 0.
        rhoc(i)  = 0.
      enddo

      do ic = nctr, 1, -1 
        icell = ictr(indc(ic))
        if ( ind(icell) .gt. 0 ) then 
           if ( hvar(1,icell) .gt. dcmax(ind(icell)) ) then 
             icmax(ind(icell)) = icell 
             dcmax(ind(icell)) = hvar(1,icell)
           endif          
           CellVolume = 2.0**(-3.0*iLv(icell))
           dmass = hvar(1,icell)*CellVolume
           rhor = 1.d0 /hvar(1,icell)
           amc(ind(icell))  = amc(ind(icell)) + dmass
           volc(ind(icell)) = volc(ind(icell)) + CellVolume
           dtemp = hvar(8,icell) * rhor 
           tempc(ind(icell)) = tempc(ind(icell)) + dtemp*dmass 
           zIac(ind(icell)) = zIac(ind(icell)) + 
     &                        hvar(izIa,icell) * rhor * dmass 
           zIIc(ind(icell)) = zIIc(ind(icell)) + 
     &                        hvar(izII,icell) * rhor * dmass 
           prc(ind(icell))  = prc(ind(icell)) + hvar(8,icell)*dmass
           rhoc(ind(icell)) = rhoc(ind(icell)) + hvar(1,icell)*dmass
           ncc(ind(icell)) = ncc(ind(icell)) + 1
        endif
      enddo
c
c.... now output GCs and their parent cloud properties
c
c
c.... efficiency of SF in proto-GC sites
c     be careful not too confuse with eps_SF in th rest of code
      e_SF = 0.5 
c.... crit. density for catastrophic s.f. in Msun/pc^3
      d_SF = 5.e3 
c.... min. GC mass to consider (in Msun) 
      gcm_min = 1.e4
c
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

      npr = index ( prefix , ' ' ) - 1
      nsf = index ( suffix , ' ' ) - 1
c
c.... GC file 
c
      fname_gas = prefix(1:npr)//'gc_cloud_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//
     &       suffix(1:nsf)//' '
      nfn = index( fname_gas , ' ' ) - 1 
      open ( 12 , file = fname_gas(1:nfn) ) 
      write(12,18) 
 18   format('# properties of cells - candidate GC sites')
      write(12,21) boxh, Om0, Oml0, Omb0, hubble      
 21   format('# Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f5.3,1x,
     &        ' Oml0= ',f5.3,1x,' Omb0= ',f7.5,1x,' h = ',f5.3)
      write(12,20) istep, t, dt, aexpn, ainit
 20   format ('# istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      write(12,29) rhoth
 29   format('# rhoth =',e9.4,' Msun/pc^3;')
      write(12,22) e_SF, d_SF, gcm_min
 22   format ('# e_SF =',e8.3,'; d_SF =',g15.8,' [Msun/pc^3]; ',
     &        'M_GC_min =',e8.3,' [Msun]')
      write(12,23)
      write(12,24)
c      write(12,25)
 23   format ('#rcell      M_GC      r_GC rho_cell Tcell   P_cell 
     &M_cell   Z_II     Z_Ia         x         y         z  
     &vx     vy     vz')
 24   format ('# pc       Msun        pc  Ms/pc^3    K     K cm^-3 
     &Msun    Zsun     Zsun                h^-1 Mpc                    
     &km/s')
 25   format ('#')
c
c.... cloud file 
c
      fname_gas = prefix(1:npr)//'clouds_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//
     &       suffix(1:nsf)//' '
      nfn = index( fname_gas , ' ' ) - 1 
      open ( 13 , file = fname_gas(1:nfn) ) 
      write(13,19) 
 19   format('# high-density/pressure clouds - candidate GC sites')
      write(13,21) boxh, Om0, Oml0, Omb0, hubble      
      write(13,20) istep, t, dt, aexpn, ainit
      write(13,29) rhoth
      write(13,22) e_SF, d_SF, gcm_min
      write(13,26)
      write(13,27)
c      write(13,25)
 26   format('# ic rcl  Mcl  rhocl  Tcl  Pcl ZIa ZII rhomax pmax mmax ')
 27   format('#             pc   Ms/pc^3  K cm^-3   Msun')

      do i = 1 , ncl 
        ic = icmax(i) 
        if ( ic .gt. 0 ) then 
          Level = iLv(ic) 
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
c
          Posx = Posx + corr_x
          Posy = Posy + corr_y
          Posz = Posz + corr_z
c
          rcell = 0.5 * CellSize(Level) * rfact 
          CellVolume = CellSize(Level)**3
          rhocell = hvar(1,ic) * rhofact 
c....     GC mass and radius      
          rsf = sqrt(0.333333*rhocell/d_SF) * rcell
          gcm = pi4 * e_SF * d_SF * rsf**3 
          gcr = rsf / e_SF 
c         if ( gcm .gt. gcm_min ) then 
c         endif          
          pcell = hvar(8,ic) * pfact ! pressure 
          rhor = 1.d0 / hvar(1,ic)
          dtemp = tfact * hvar(8,ic) * rhor
c....     metallicity in units of solar 
          dzII = 50.* hvar(izII,ic) * rhor
          dzIa = 50.* hvar(izIa,ic) * rhor
          dvcx = hvar(3,ic) * rhor * vfact 
          dvcy = hvar(4,ic) * rhor * vfact
          dvcz = hvar(5,ic) * rhor * vfact
          dmass = hvar(1,ic)* CellVolume * amfact 
          xcell = (Posx - xc) * r0
          ycell = (Posy - yc) * r0
          zcell = (Posz - zc) * r0
c
c....     output globular clusters 
c
          write(12,80) rcell, gcm, gcr, 
     &      rhocell, dtemp, pcell, dmass, 
     &      dzII, dzIa, xcell, ycell, zcell, 
     &      dvcx, dvcy, dvcz, i 
 80       format(f6.1,2x,e11.5,1x,f6.1,1x,
     &           f6.1,1x,5(e8.3,1x),1x,3(f9.6,1x),1x,3(f6.1,1x),i4)
c
c....     output cloud properties 
c
          rcloud = (volc(i))**0.333333 * rfact 
          tcloud = tempc(i) / amc(i) * tfact 
          zIacl = zIac(i) / amc(i) / 0.019
          zIIcl = zIIc(i) / amc(i) / 0.019
          prcl  = prc(i) / amc(i) * pfact 
          rhocl = rhoc(i) / amc(i) * rhofact 
          write(13,81) i, rcloud, amc(i)*aM0, rhocl, tcloud, prcl, 
     &       zIacl, zIIcl, rhocell, pcell, dmass, ncc(i)
 81       format(i4,1x,10(e8.3,1x),i4)
        else
          write(*,*) 'Identify_GC: something is wrong... imax =',0
        endif
      enddo

      close ( 12 )
      close ( 13 ) 
c

      return
      end
c
c
c     --------------------------------------------------------
      subroutine Estimate_GC_Mass ( rcell, rhocell, gcm, gcr )
c     --------------------------------------------------------
c     
c     input:  rcell - cell size in physical pc 
c             rhocell - cell physical density in Msun/pc^3
c     output: gcm - estimate of the GC mass in Msun
c             gcr - "size" of the GC in physical pc
c     NOTE: no h in any of the above quantities!
c      
      real*8 pi, pi4, pi43
      parameter ( pi   = 3.14159265358978245      )
      parameter ( pi4  = 4.0 * pi                 )
      parameter ( pi43 = pi4 / 3.0                )

      real rcell, rhocell, gcm, gcr

c.... efficiency of SF in proto-GC sites
      e_SF = 0.5 ! be careful not too confuse with eps_SF in th rest of code
c.... critical density for catastrophic starformation in Msun/pc^3
      d_SF = 5.e3 
c.... GC mass and radius      
      rsf = sqrt(0.333333*rhocell/d_SF) * rcell
      gcm = pi4 * e_SF * d_SF * rsf**3 
      gcr = rsf / e_SF 
c
      return
      end
c
c     -----------------------------------------------------
      subroutine Output_HighDensity ( aboxh ,  rhoth, pth )
c     ----------------------------------------------------
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'

      character*5 fstep
      real svd(nsvx,nsvy,nsvz)
      real p_d1(nsvy,nsvz), p_d2(nsvy,nsvz), p_d3(nsvy,nsvz)
      real p_t1(nsvy,nsvz), p_t2(nsvy,nsvz), p_t3(nsvy,nsvz)
      real*8 sum1, sum2, sum3, sum4, vmax1, vmax2, em
      common / SVDUM / svd, p_d1, p_d2, p_d3, p_t1, p_t2, p_t3
      character*256 fname_gas
      logical lpointsout, l_force_center
      integer imovstep
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*2 chbox
      real*8 pfact, rhofact 
      parameter ( nctrmax = 100000 )
      integer ictr(nctrmax), indc(nctrmax)
      real pctr(nctrmax), dctr(nctrmax)

#ifdef PRESSUREFLOOR
      real*8 f_J, f_eJ
      parameter ( f_J = 2.5  ) ! pressure floor parameter
#endif

#ifdef PRESSUREFLOOR
      DO Level = MinL_Jeans , MaxLevelNow
c
c.... artificial pressure floor a-la Machacek et al. 2
        if ( Level .ge. MinL_Jeans ) then
          f_eJ = 0.47746 * f_J * aexpn * 
     &           CellSize(MaxLevel)**2 * (gamma - 1.0)
        else
          f_eJ = 0.0
        endif
        if ( Level .eq. MinLevel ) then
          do i1 = 1 , ncell0
            E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
            E_th = f_eJ * hvar(1,i1) * hvar(1,i1) 
            hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
          enddo 
        else
          call Select_Cells ( Level , nLevel )
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

      do ic = 1 , nctot
        ind(ic) = 0
      enddo
c
c.... output profiles 
c
      xc = 32.95605469
      yc = 13.58300781 
      zc = 39.06152344
      rmin = 5.e-3/r0
      rmax = aboxh/r0
      nrbin = 45 
      ibintype = 0 
      deltavir = 340
      deltacrit = 500
      ifindcenter = 3 
      lpointsout = .false.  
      l_force_center = .false.  

      lpath = index(path    , ' ') - 1

      fname_gas  = path(1:lpath)//'/pro'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'_'//chbox//'.dat '

      call Output_Average_Profiles_Stars ( xc , yc , zc , rvir,  
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, deltacrit, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fname_gas ) 
c
        rhofact = rho0 * 1.d-18 / aexpn**3
        tfact = (gamma-1.0) * T_0 / aexpn**2 
        pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
        nctr = 0 
c
        DO Level = MinLevel , MaxLevelNow
          CellVolume = 1.0 * 2.0**(-3.0*Level)
          IF ( Level .eq. MinLevel ) THEN
          ELSE
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

                  Posx = Posx + corr_x
                  Posy = Posy + corr_y
                  Posz = Posz + corr_z

c                  rr = sqrt( (Posx-xc)**2 +
c     &                       (Posy-yc)**2 +
c     &                       (Posz-zc)**2 )
c                  if ( rr .lt. rmax ) then
                    rhocell = hvar(1,ic) * rhofact 
                    pcell = hvar(8,ic) * pfact 
                    if ( rhocell .gt. rhoth ) then 
                      nctr = nctr + 1
                      if ( nctr .gt. nctrmax ) then
                        write(*,*) 'nctr exceeds nctrmax =',nctrmax
                        write(*,*) 'increase nctrmax and rerun'
                        stop
                      endif
                      ictr(nctr) = ic
                      pctr(nctr) = pcell
                      dctr(nctr) = rhocell
c                      write(22,'(i8,1x,i1,1x,3(f10.6,1x),2(g11.4,1x))') 
c     &                 ic, Level, Posx,Posy,Posz,rhocell,pcell
                    endif
c                  endif
                endif
              enddo
            enddo
          ENDIF
        ENDDO      
c
      call indexx ( nctr , dctr , indc )

      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

      fname_gas = 'dpro_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//' '
      nfn = index( fname_gas , ' ' ) - 1 

      open ( 12 , file = fname_gas(1:nfn) ) 
      fname_gas = 'dc_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//' '
      nfn = index( fname_gas , ' ' ) - 1 

      open ( 11 , file = fname_gas(1:nfn) ) 
      fname_gas = 'vpro_a'//
     &       digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)//' '
      nfn = index( fname_gas , ' ' ) - 1 

      open ( 13 , file = fname_gas(1:nfn) ) 

      dc = 10.0 / rhofact 
      dout = 1.0 / rhofact
      rmin = 30.0 / 1.e6 / r0 / aexpn ! /h pc-> code units
      rmax = 1500.0 / 1.e6 / r0 / aexpn ! /h pc-> code units
      nrbin = 30 
      Tcold = 1.e4

      do ic = nctr , 1 , -1 
        icell = ictr(indc(ic))
        write(*,*) 'processing center',ic
        if ( ind(icell) .eq. 0 ) then 
          call Ps ( icell , Posx , Posy , Posz )
          call DensityCenter_Profiles_HP 
     &         ( dc, dout, Tcold,  Posx, Posy, Posz,
     &           rmin , rmax , nrbin , rmax, 0 )
        endif
      enddo
c
      close ( 11 ) 
      close ( 12 ) 
      close ( 13 ) 
c      
      return
      end
c
c     -----------------------------------------------------------------
      subroutine DensityCenter_Profiles_HP ( dc, dout, Tcold, 
     &                            xc , yc , zc , rmin , rmax , 
     &                            nrbin , rmass, rvir, ictype )
c     -----------------------------------------------------------------
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    dc : central density in Msun/pc^3
c              dout : minimum density defining the outer radius of the cloud
c              Tcold - temperature (in K) to use to define cold gas 
c      integer: ictype - center type (0 - GC, 1 - halos) to determine what
c              to output   
c
c      count mass to dout or rmass whichever comes first 
c            rmass is in code units
c
c      itide = 1 if the halo is tidally truncated, 0 if not (rmass=rvir) 
c
#     include "a_def.h"
      include 'a_hfind.h'
c
      real xc , yc , zc , rmin , rmax , deltavir, deltacrit
      real rmass, rvir
      integer nrbin, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a, pfact
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 200000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real vcx(npoints), vcy(npoints), vcz(npoints), cw(npoints)
      real ptemp(npoints), pz(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pz
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dgbin(0:nmaxbins),dgcbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), dtibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real ndmp(0:nmaxbins),entbin(0:nmaxbins)
      real zbin(0:nmaxbins), zsbin(0:nmaxbins), zsbinIa(0:nmaxbins)
      real ppvx(0:nmaxbins), ppvy(0:nmaxbins), ppvz(0:nmaxbins)
      real vvrms(0:nmaxbins)
      real svx(0:nmaxbins), svy(0:nmaxbins), svz(0:nmaxbins)
      real svrms(0:nmaxbins)
      real gvx(0:nmaxbins), gvy(0:nmaxbins), gvz(0:nmaxbins)
      real gvrms(0:nmaxbins), gm(0:nmaxbins)
      real*8 amst(0:nmaxbins), amstot(0:nmaxbins), nst(0:nmaxbins)
      real*8 amgtot(0:nmaxbins), amgctot(0:nmaxbins)
      real*8 amdmtot(0:nmaxbins), tsbin(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(0:nmaxbins), ncbin(0:nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn
      real*8 rhofact, age
      real*8 zave, amgave, zsave, tave, amsave

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) 

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rhofact = rho0 * 1.d-18 / aexpn**3

      call Get_MaxLevelNow ()

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
c
c....   metallicity in units of solar 
c
          pz(i) = 50. * (hvar(izII,idcell)+hvar(izIa,idcell)) /
     &            hvar(1,idcell)
          rhor = 1.d0 / hvar(1,idcell)
          vcx(i) = hvar(3,idcell) * rhor 
          vcy(i) = hvar(4,idcell) * rhor 
          vcz(i) = hvar(5,idcell) * rhor 
          cw(i)  = hvar(1,idcell) * CellVol(iLv(idcell))
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
        dgcbin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        nbin(i) = 0
        ncbin(i) = 0
        ndmp(i) = 0.
        nst(i) = 0
        amst(i) = 0
        amstot(i) = 0
        amgtot(i) = 0.
        amgctot(i) = 0.
        amdmtot(i) = 0.
        zbin(i) = 0.
        zsbin(i) = 0.
        zsbinIa(i) = 0.
        tsbin(i) = 0.
        ppvx(i) = 0.
        ppvy(i) = 0.
        ppvz(i) = 0.
        vvrms(i) = 0.
        svx(i) = 0.
        svy(i) = 0.
        svz(i) = 0.
        svrms(i) = 0.
        gvx(i) = 0.
        gvy(i) = 0.
        gvz(i) = 0.
        gvrms(i) = 0.
        gm(i) = 0.
      enddo

      nstars = nsp(nspec,2) - nsp(nspec,1) + 1
      istar1 = nsp(nspec,1)

      do i = 1 , np 
        xx = x(i)
        yy = y(i)
        zz = z(i)
        diff_x = xc - xx
        diff_y = yc - yy
        diff_z = zc - zz
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

c....  correct coordinate to insure periodic bundaries
        xx = xx + corr_x
        yy = yy + corr_y
        zz = zz + corr_z

        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          if ( i .lt. istar1 ) then ! DM particle
            ddmbin(ibin) = ddmbin(ibin) + pw(i)
            ndmp(ibin) = ndmp(ibin) + 1
            ppvx(ibin) = ppvx(ibin) + vx(i)*pw(i) 
            ppvy(ibin) = ppvy(ibin) + vy(i)*pw(i)
            ppvz(ibin) = ppvz(ibin) + vz(i)*pw(i)
            vvrms(ibin) = vvrms(ibin) + 
     &                   (vx(i)**2 + vy(i)**2 + vz(i)**2)*pw(i)
          else ! stellar particle
            amst(ibin) = amst(ibin) + pw(i)
            nst(ibin) = nst(ibin) + 1
            zsbin(ibin) = zsbin(ibin) + pw(i)*zstII(i-istar1+1)
            zsbinIa(ibin) = zsbinIa(ibin) + pw(i)*zstIa(i-istar1+1)
            svx(ibin) = svx(ibin) + vx(i)*pw(i) 
            svy(ibin) = svy(ibin) + vy(i)*pw(i)
            svz(ibin) = svz(ibin) + vz(i)*pw(i)
            svrms(ibin) = svrms(ibin) + 
     &                   (vx(i)**2 + vy(i)**2 + vz(i)**2)*pw(i)
            tsbin(ibin) = tsbin(ibin) + tbirth(i-istar1+1)*pw(i)
          endif
        endif
      enddo

      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          if ( ptemp(ii) .le. Tcold ) then ! cold gas 
            ncbin(ibin) = ncbin(ibin) + 1
            dgcbin(ibin)  = dgcbin(ibin) + pd(ii)
          endif
          pbin(ibin) = pbin(ibin) + pp(ii)*cw(ii)
          zbin(ibin) = zbin(ibin) + pz(ii)*cw(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)*cw(ii) ! mass weighted
          gvx(ibin) = gvx(ibin) + vcx(ii)*cw(ii) 
          gvy(ibin) = gvy(ibin) + vcy(ii)*cw(ii)
          gvz(ibin) = gvz(ibin) + vcz(ii)*cw(ii)
          gvrms(ibin) = gvrms(ibin) + 
     &           (vcx(ii)**2 + vcy(ii)**2 + vcz(ii)**2)*cw(ii)
          gm(ibin) = gm(ibin) + cw(ii)
        endif
      enddo
c
c.... prepare profiles
c
      vfact = v0 / aexpn
      do i = 0 , nbins
        if ( nbin(i) .gt. 0 ) then 
          rnp = 1.d0 / gm(i)
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
          pbin(i) = pbin(i) * rnp 
          tbin(i) = tbin(i) * rnp
          zbin(i) = zbin(i) * rnp
          gvx(i) = gvx(i) * rnp 
          gvy(i) = gvy(i) * rnp 
          gvz(i) = gvz(i) * rnp 
          vmn = gvx(i)**2 + gvy(i)**2 + gvz(i)**2
          gvrms(i) = sqrt(abs(gvrms(i)*rnp - vmn)) * vfact
          gvx(i) = gvx(i) * vfact 
          gvy(i) = gvy(i) * vfact 
          gvz(i) = gvz(i) * vfact         
          dgcbin(i) = dgcbin(i) / nbin(i)
        endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.0
      dti = 0.0
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      ivir = 0
      irflag = 0 

      dlout = log10(dout)
c
c.... initialize variables for averages 
c
      zave  = 0.0 ! average metallicity of gas 
      zsave = 0.0 ! average metallicity of stars 
      tave  = 0.0 ! average stellar age
      amgave = 0.0 ! sum of the bin gas masses 
      amsave = 0.0 ! sume of the stellar masses in bins

      amst_tot = 0.0 
      amg_tot = 0.0 
      amgc_tot = 0.0 ! cold gas mass
      volr1 = 0.

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
          dmg = dgbin(i)*volr
          amg_tot = amg_tot + dmg 
          amgtot(i) = amg_tot
          amgc_tot = amgc_tot + dgcbin(i) * volr
          amgctot(i) = amgc_tot 
          amdmtot(i) = dmi
          dti = dti + ddmbin(i) + dgbin(i)*volr + amst(i)
          pbin(i) = pbin(i) * pfact 
c
c....     compute averages within rout
c
          if ( irflag .eq. 0 ) then ! if outer radius not reached yet
            zave   = zave + zbin(i) * dmg
            amgave = amgave + dmg 
            zsave  = zsave + (zsbin(i) + zsbinIa(i))
            tave   = tave + tsbin(i) 
            amsave = amsave + amst(i)
          endif
c
          ddmbin(i) = ddmbin(i)/volr 
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i) = dti / (4.18879 * rr**3)
c          dbi1 = dtibin(i-1)-ddmibin(i-1) ! baryons only
c          dbi2 = dtibin(i)-ddmibin(i)
          dbi1 = dtibin(i-1) ! total density 
          dbi2 = dtibin(i)

c
c...      compute the virial radius
c
          if ( (ibin .gt. 0 .and.
     &         dbi2 .lt. dout .and.
     &         dbi1 .ge. dout) ) then 
             rrl = log10(rr)
             rll = log10(rl) 
             rri = 1.d0 / (rrl - rll)
             dlbi1 = log10(dbi1)
             dlbi2 = log10(dbi2)
             rvir = 10.d0**((dlout * (rrl - rll) + rll*dlbi2 - 
     &            rrl*dlbi1) / (dlbi2 - dlbi1))
          endif

          if ( (ibin .gt. 0 .and.
     &         dbi2 .lt. dout .and.
     &         dbi1 .ge. dout) .and. irflag .eq. 0 ) then 
            ivir = i
            rrl = log10(rr)
            rll = log10(rl) 
            rri = 1.d0 / (rrl - rll)
            dlbi1 = log10(dbi1)
            dlbi2 = log10(dbi2)
            rout = 10.d0**((dlout * (rrl - rll) + rll*dlbi2 - 
     &             rrl*dlbi1) / (dlbi2 - dlbi1))
            rlout = log10(rout)
c
c....       gas mass
c
            if ( amgtot(i-1) .gt. 0 ) then 
              amgl1 = log10(amgtot(i-1))
            else
              amgl1 = -15.0
            endif
            if ( amgtot(i) .gt. 0 ) then 
              amgl2 = log10(amgtot(i))
            else
              amgl2 = -15.0
            endif
            ah = (amgl2 - amgl1) * rri
            bh = amgl1 - ah * rll 
            aM_gas = 10.d0**(ah*rlout + bh)
c
c....       cold gas mass
c
            if ( amgctot(i-1) .gt. 0 ) then 
              amgcl1 = log10(amgctot(i-1))
            else
              amgcl1 = -15.0
            endif
            if ( amgctot(i) .gt. 0 ) then 
              amgcl2 = log10(amgctot(i))
            else
              amgcl2 = -15.0
            endif
            ah = (amgcl2 - amgcl1) * rri
            bh = amgcl1 - ah * rll 
            aM_cgas = 10.d0**(ah*rlout + bh)
c
c....       stellar mass
c 
            if ( amstot(i-1) .gt. 0. ) then 
              amsl1 = log10(amstot(i-1))
            else
              amsl1 = -15.0 
            endif
            if ( amstot(i) .gt. 0. ) then 
              amsl2 = log10(amstot(i))
            else
              amsl2 = -15.0
            endif
            ah = (amsl2 - amsl1) * rri
            bh = amsl1 - ah * rll 
            aM_st = 10.d0**(ah*rlout + bh)
c
c....       DM mass
c
            if ( amdmtot(i-1) .gt. 0.0 ) then 
              amdl1 = log10(amdmtot(i-1))
            else
              amdl1 = -15.0 
            endif
            if ( amdmtot(i) .gt. 0.0 ) then 
              amdl2 = log10( amdmtot(i) )
            else
              amdl2 = -15.0 
            endif
            ah = (amdl2 - amdl1) * rri
            bh = amdl1 - ah * rll 
            aM_dm   = 10.d0**(ah*rlout + bh)
c            write(*,*)'Mdm:',amdmtot(i-1)*aM0*hubble,
c     &        amdmtot(i)*aM0*hubble,aM_dm*aM0*hubble

            aM_gas  = aM_gas * aM0 * hubble
            aM_cgas = aM_cgas * aM0 * hubble
            aM_st   = aM_st  * aM0 * hubble
            aM_dm   = aM_dm  * aM0 * hubble
            aMout   = aM_gas + aM_st 
            irflag  = 1 
          endif
c
c...  estimate mass if rmass is reached
c
          if ( (ibin .gt. 0 .and.
     &         rr .gt. rmass .and.
     &         rl .le. rmass) .and. irflag .eq. 0 ) then 
            ivir = i
            rrl = log10(rr)
            rll = log10(rl) 
            rri = 1.d0 / (rrl - rll)
            rout = rmass
            rlout = log10(rout)
c
c....       gas mass
c
            if ( amgtot(i-1) .gt. 0 ) then 
              amgl1 = log10(amgtot(i-1))
            else
              amgl1 = -15.0
            endif
            if ( amgtot(i) .gt. 0 ) then 
              amgl2 = log10(amgtot(i))
            else
              amgl2 = -15.0
            endif
            ah = (amgl2 - amgl1) * rri
            bh = amgl1 - ah * rll 
            aM_gas = 10.d0**(ah*rlout + bh)
c
c....       cold gas mass
c
            if ( amgctot(i-1) .gt. 0. ) then 
              amgcl1 = log10(amgctot(i-1))
            else
              amgcl1 = -15.0
            endif
            if ( amgctot(i) .gt. 0. ) then 
              amgcl2 = log10(amgctot(i))
            else
              amgcl2 = -15.0
            endif
            ah = (amgcl2 - amgcl1) * rri
            bh = amgcl1 - ah * rll 
            aM_cgas = 10.d0**(ah*rlout + bh)
c
c....       stellar mass
c 
            if ( amstot(i-1) .gt. 0. ) then 
              amsl1 = log10(amstot(i-1))
            else
              amsl1 = -15.0 
            endif
            if ( amstot(i) .gt. 0. ) then 
              amsl2 = log10(amstot(i))
            else
              amsl2 = -15.0
            endif
            ah = (amsl2 - amsl1) * rri
            bh = amsl1 - ah * rll 
            aM_st = 10.d0**(ah*rlout + bh)
c
c....       DM mass
c
            if ( amdmtot(i-1) .gt. 0.0 ) then 
              amdl1 = log10(amdmtot(i-1))
            else
              amdl1 = -15.0 
            endif
            if ( amdmtot(i) .gt. 0.0 ) then 
              amdl2 = log10( amdmtot(i) )
            else
              amdl2 = -15.0 
            endif
            ah = (amdl2 - amdl1) * rri
            bh = amdl1 - ah * rll 
            aM_dm   = 10.d0**(ah*rlout + bh)

            aM_gas  = aM_gas * aM0 * hubble
            aM_cgas = aM_cgas * aM0 * hubble
            aM_st   = aM_st  * aM0 * hubble
            aM_dm   = aM_dm  * aM0 * hubble
            aMout   = aM_gas + aM_st 
            irflag  = 1 
          endif

         
c          if ( ibin .gt. 0 .and.
c     &         tbin(i) .gt. 1.e3 .and. tbin(i-1) .lt. 1.e3 ) then 
c             irflag = 1
c             rout = rr
c             ivir = i
c          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
        volr1 = volr 
      enddo
      if ( ivir .eq. 0 ) then 
        rout = rmax
        ivir = nbins
        aM_gas = amgtot(ivir) * aM0 * hubble
        aM_cgas = amgctot(ivir) * aM0 * hubble
        rr = 10.**(rlmin + float(ivir)*drl)
        aM_dm = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
        aM_st = amstot(ivir) * aM0 * hubble 
        aMout = aM_gas + aM_st
      endif
c
      xcd=(xc-1.0)*r0
      ycd=(yc-1.0)*r0
      zcd=(zc-1.0)*r0
      rfact = r0 * 1.e3 * aexpn / hubble
      vfact = v0 / aexpn
      rdout  = rout * rfact   ! min(rt,rvir) in proper kpc
      rvdout = rvir * rfact   ! rvir in proper kpc 
      
c
c.... averages
c
      if ( amgave .gt. 0.0 ) then 
        zave = zave / amgave ! already solar
      endif
      if ( amsave .gt. 0.0 ) then 
        zsave = 50.0 * zsave / amsave ! -> solar 
        td = tave / amsave
        tave = age ( 1.d0*td ) ! code units -> Gyr
      endif
      
      if ( ictype .eq. 1 ) then ! DM halos
        write(11,10) xcd, ycd, zcd, rdout, rvdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc)
        write(12,10) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc)
        write(13,10) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc)
      endif
      if ( ictype .eq. 0 ) then ! GC
        write(11,11) xcd, ycd, zcd, rdout, rvdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave
        write(12,11) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave
        write(13,11) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm,
     &               zave, zsave, tave
      endif
 10   format (3(f9.6,1x),1x,2(f8.3,1x),1x,6(e9.3,1x),f6.3,1x,e9.3)
 11   format (3(f9.6,1x),1x,2(f8.3,1x),1x,6(e9.3,1x),f6.3)
c
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) 
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
        abtot = (amgtot(i)+ amstot(i))/(4.18879 * rr**3)
        if ( amst(i) .gt. 0.0 ) then ! correct for stellar metallicity
          rnp = 1.0 / amst(i)
          zbin(i) = (0.02*zbin(i)*dgbin(i)*volr) / volr /
     &               dgbin(i)*50.
          zsbin(i) = zsbin(i) * rnp * 50.
          zsbinIa(i) = zsbinIa(i) * rnp * 50.
          svx(i) = svx(i) * rnp
          svy(i) = svy(i) * rnp
          svz(i) = svz(i) * rnp
          vmn = svx(i)**2 + svy(i)**2 + svz(i)**2
          svrms(i) = sqrt(abs(svrms(i)*rnp - vmn)) * vfact
          svx(i) = svx(i) * vfact 
          svy(i) = svy(i) * vfact 
          svz(i) = svz(i) * vfact 
        endif
        if ( ddmbin(i) .gt. 0.0 ) then 
          dm = ddmbin(i) * volr 
          rnp = 1.0 / dm
          ppvx(i) = ppvx(i) * rnp
          ppvy(i) = ppvy(i) * rnp
          ppvz(i) = ppvz(i) * rnp
          vmn = ppvx(i)**2 + ppvy(i)**2 + ppvz(i)**2
          vvrms(i) = sqrt(abs(vvrms(i)*rnp - vmn)) * vfact
          ppvx(i) = ppvx(i) * vfact 
          ppvy(i) = ppvy(i) * vfact 
          ppvz(i) = ppvz(i) * vfact           
        endif
        aMdm = ddmibin(i) * (4.18879 * rr**3) * aM0 * hubble
        abtot = abtot*rhofact
        rhost = amst(i) / volr 
        hmassi = (amgtot(i)+amstot(i))*aM0 + aMdm/hubble
        rcirc = rr * r0/hubble * aexpn * 3.0856d24
        vcirc1 = 1.15206d8 * sqrt(hmassi/rcirc)
        hmassi = (amgtot(i)+amstot(i))*aM0
        vcirc2 = 1.15206d8 * sqrt(hmassi/rcirc)
        hmassi = (amstot(i))*aM0
        vcirc3 = 1.15206d8 * sqrt(hmassi/rcirc)
          write(12,18)
     &         rmid*rfact, rr*rfact, dtibin(i)*rhofact,
     &         dgbin(i)*rhofact,rhost*rhofact, ddmbin(i)*rhofact,
     &         amgtot(i)*aM0*hubble, 
     &         amgctot(i)*aM0*hubble, 
     &         amstot(i)*aM0*hubble,
     &         aMdm,tbin(i),zbin(i),zsbin(i),zsbinIa(i),
     &         int(ndmp(i)),int(nst(i))
          write(13,19)
     &         rmid*rfact, rr*rfact, 
     &         amgtot(i)*aM0*hubble,
     &         amstot(i)*aM0*hubble,
     &         aMdm, gvrms(i),svrms(i),vvrms(i), vcirc1, vcirc2, vcirc3,
     &         int(ndmp(i)),int(nst(i))
      enddo
 18   format(2(f8.3,1x),12(e9.3,1x),i6,1x,i6)
 19   format(2(f8.3,1x),1x,9(e9.3,1x),i6,1x,i6)
c
c
c.... now label cells falling within the radius of clump
c
      do i = 1 , npoints
        rr    = rout * RANDd ( m )
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
        ind(idcell) = 1
      enddo

      return
      end

c     -----------------------------------------------------------------
      subroutine DensityCenter_Profiles_HP2 ( Pout, Tout , Tcold, 
     &                            xc , yc , zc ,
     &                            rmin , rmax , nrbin , rmass, ictype )
c     -----------------------------------------------------------------
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    dc : central density in Msun/pc^3
c              dout : minimum density defining the outer radius of the cloud
c              Tcold - temperature (in K) to use to define cold gas 
c      integer: ictype - center type (0 - GC, 1 - halos) to determine what
c              to output   
c
c     count mass to Tout, Pout, rmass whichever comes first 
c            rmass is in code units
c
c     Tout and Pout are outer temperature and pressure
c
#     include "a_def.h"
      include 'a_hfind.h'
c
      real xc , yc , zc , rmin , rmax , deltavir, deltacrit
      integer nrbin, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a, pfact
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 200000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real vcx(npoints), vcy(npoints), vcz(npoints), cw(npoints)
      real ptemp(npoints), pz(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pz
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dgbin(0:nmaxbins),dgcbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), dtibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real ndmp(0:nmaxbins),entbin(0:nmaxbins)
      real zbin(0:nmaxbins), zsbin(0:nmaxbins), zsbinIa(0:nmaxbins)
      real ppvx(0:nmaxbins), ppvy(0:nmaxbins), ppvz(0:nmaxbins)
      real vvrms(0:nmaxbins)
      real svx(0:nmaxbins), svy(0:nmaxbins), svz(0:nmaxbins)
      real svrms(0:nmaxbins)
      real gvx(0:nmaxbins), gvy(0:nmaxbins), gvz(0:nmaxbins)
      real gvrms(0:nmaxbins), gm(0:nmaxbins)
      real*8 amst(0:nmaxbins), amstot(0:nmaxbins), nst(0:nmaxbins)
      real*8 amgtot(0:nmaxbins), amgctot(0:nmaxbins)
      real*8 amdmtot(0:nmaxbins), tsbin(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(0:nmaxbins), ncbin(0:nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn
      real*8 rhofact, age
      real*8 zave, amgave, zsave, tave, amsave
      logical lcond1, lcond2

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) 

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      ai = 1.0 / aexpn
      gamma1 = gamma - 1.0 
      rhofact = rho0 * 1.d-18 / aexpn**3
c.... factor converting gas density in code units to n_gas in cm^-3
      fact_nH = 1.12d-5*hubble*hubble*Om0* ai**3

      call Get_MaxLevelNow ()

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
c
c....   metallicity in units of solar 
c
          pz(i) = 50. * (hvar(izII,idcell)+hvar(izIa,idcell)) /
     &            hvar(1,idcell)
          rhor = 1.d0 / hvar(1,idcell)
          vcx(i) = hvar(3,idcell) * rhor 
          vcy(i) = hvar(4,idcell) * rhor 
          vcz(i) = hvar(5,idcell) * rhor 
          cw(i)  = hvar(1,idcell) * CellVol(iLv(idcell))
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
        dgcbin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        nbin(i) = 0
        ncbin(i) = 0
        ndmp(i) = 0.
        nst(i) = 0
        amst(i) = 0
        amstot(i) = 0
        amgtot(i) = 0.
        amgctot(i) = 0.
        amdmtot(i) = 0.
        zbin(i) = 0.
        zsbin(i) = 0.
        zsbinIa(i) = 0.
        tsbin(i) = 0.
        ppvx(i) = 0.
        ppvy(i) = 0.
        ppvz(i) = 0.
        vvrms(i) = 0.
        svx(i) = 0.
        svy(i) = 0.
        svz(i) = 0.
        svrms(i) = 0.
        gvx(i) = 0.
        gvy(i) = 0.
        gvz(i) = 0.
        gvrms(i) = 0.
        gm(i) = 0.
      enddo

      nstars = nsp(nspec,2) - nsp(nspec,1) + 1
      istar1 = nsp(nspec,1)
      write(*,*) ' hello : ',nspecies, nspec, nstars, istar1, np

      do i = 1 , np 
        xx = x(i)
        yy = y(i)
        zz = z(i)
        diff_x = xc - xx
        diff_y = yc - yy
        diff_z = zc - zz
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

c....  correct coordinate to insure periodic bundaries
        xx = xx + corr_x
        yy = yy + corr_y
        zz = zz + corr_z

        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          if ( i .lt. istar1 ) then ! DM particle
            ddmbin(ibin) = ddmbin(ibin) + pw(i)
            ndmp(ibin) = ndmp(ibin) + 1
            ppvx(ibin) = ppvx(ibin) + vx(i)*pw(i) 
            ppvy(ibin) = ppvy(ibin) + vy(i)*pw(i)
            ppvz(ibin) = ppvz(ibin) + vz(i)*pw(i)
            vvrms(ibin) = vvrms(ibin) + 
     &                   (vx(i)**2 + vy(i)**2 + vz(i)**2)*pw(i)
          else ! stellar particle
            amst(ibin) = amst(ibin) + pw(i)
            nst(ibin) = nst(ibin) + 1
            zsbin(ibin) = zsbin(ibin) + pw(i)*zstII(i-istar1+1)
            zsbinIa(ibin) = zsbinIa(ibin) + pw(i)*zstIa(i-istar1+1)
            svx(ibin) = svx(ibin) + vx(i)*pw(i) 
            svy(ibin) = svy(ibin) + vy(i)*pw(i)
            svz(ibin) = svz(ibin) + vz(i)*pw(i)
            svrms(ibin) = svrms(ibin) + 
     &                   (vx(i)**2 + vy(i)**2 + vz(i)**2)*pw(i)
            tsbin(ibin) = tsbin(ibin) + tbirth(i-istar1+1)*pw(i)
          endif
        endif
      enddo

      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          if ( ptemp(ii) .le. Tcold ) then ! cold gas 
            ncbin(ibin) = ncbin(ibin) + 1
            dgcbin(ibin)  = dgcbin(ibin) + pd(ii)
          endif
          pbin(ibin) = pbin(ibin) + pp(ii)*cw(ii)
          zbin(ibin) = zbin(ibin) + pz(ii)*cw(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)*cw(ii) ! mass weighted
          gvx(ibin) = gvx(ibin) + vcx(ii)*cw(ii) 
          gvy(ibin) = gvy(ibin) + vcy(ii)*cw(ii)
          gvz(ibin) = gvz(ibin) + vcz(ii)*cw(ii)
          gvrms(ibin) = gvrms(ibin) + 
     &           (vcx(ii)**2 + vcy(ii)**2 + vcz(ii)**2)*cw(ii)
          gm(ibin) = gm(ibin) + cw(ii)
        endif
      enddo
c
c.... prepare profiles
c
      vfact = v0 / aexpn
      do i = 0 , nbins
        if ( nbin(i) .gt. 0 ) then 
          rnp = 1.d0 / gm(i)
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
          pbin(i) = pbin(i) * rnp 
          tbin(i) = tbin(i) * rnp
          zbin(i) = zbin(i) * rnp
          gvx(i) = gvx(i) * rnp 
          gvy(i) = gvy(i) * rnp 
          gvz(i) = gvz(i) * rnp 
          vmn = gvx(i)**2 + gvy(i)**2 + gvz(i)**2
          gvrms(i) = sqrt(abs(gvrms(i)*rnp - vmn)) * vfact
          gvx(i) = gvx(i) * vfact 
          gvy(i) = gvy(i) * vfact 
          gvz(i) = gvz(i) * vfact         
          dgcbin(i) = dgcbin(i) / nbin(i)
        endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.0
      dti = 0.0
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      ivir = 0
      irflag = 0 

      dlout = log10(dout)
      dPlout = log10(Pout)
      dTlout = log10(Tout)
c
c.... initialize variables for averages 
c
      zave  = 0.0 ! average metallicity of gas 
      zsave = 0.0 ! average metallicity of stars 
      tave  = 0.0 ! average stellar age
      amgave = 0.0 ! sum of the bin gas masses 
      amsave = 0.0 ! sume of the stellar masses in bins

      amst_tot = 0.0 
      amg_tot = 0.0 
      amgc_tot = 0.0 ! cold gas mass
      volr1 = 0.

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
          dmg = dgbin(i)*volr
          amg_tot = amg_tot + dmg 
          amgtot(i) = amg_tot
          amgc_tot = amgc_tot + dgcbin(i) * volr
          amgctot(i) = amgc_tot 
          amdmtot(i) = dmi
          dti = dti + ddmbin(i) + dgbin(i)*volr + amst(i)
          pbin(i) = pbin(i) * pfact 
c
c....     compute averages within rout
c
          if ( irflag .eq. 0 ) then ! if outer radius not reached yet
            zave   = zave + zbin(i) * dmg
            amgave = amgave + dmg 
            zsave  = zsave + (zsbin(i) + zsbinIa(i))
            tave   = tave + tsbin(i) 
            amsave = amsave + amst(i)
          endif
c
          ddmbin(i) = ddmbin(i)/volr 
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i) = dti / (4.18879 * rr**3)
c          dbi1 = dtibin(i-1)-ddmibin(i-1) ! baryons only
c          dbi2 = dtibin(i)-ddmibin(i)
          dbi1 = dtibin(i-1) ! total density 
          dbi2 = dtibin(i)

          dn_H1 = dgbin(i-1) * fact_nH ! code density -> n cm^-3
          dT_1 = tbin(i-1) ! temperature in K
          dpn1  = dn_H1 * dT_1 ! nT in cm^-3 K
          dn_H2 = dgbin(i) * fact_nH ! code density -> n cm^-3
          dT_2 = tbin(i) ! temperature in K
          dpn2  = dn_H2 * dT_2 ! nT in cm^-3 K
c
          lcond1 = (ibin .gt. 0 .and.
     &         dpn2 .lt. Pout .and.
     &         dpn1 .ge. Pout) .and. irflag .eq. 0 
          lcond2 = (ibin .gt. 0 .and.
     &         dT_2 .gt. Tout .and.
     &         dT_1 .le. Tout) .and. irflag .eq. 0 
c
          write(*,*) 'pressure:', dpn1, dpn2, Pout, lcond1
c          write(*,*) 'temperature:', dT_1, dT_2, Tout, lcond2
          if ( lcond1 ) then 
            ivir = i
            rrl = log10(rr)
            rll = log10(rl) 
            rri = 1.d0 / (rrl - rll)
            dlbi1 = log10(dpn1)
            dlbi2 = log10(dpn2)
            rout = 10.d0**((dPlout * (rrl - rll) + rll*dlbi2 - 
     &             rrl*dlbi1) / (dlbi2 - dlbi1))
c            dlbi1 = log10(dT_1)
c            dlbi2 = log10(dT_2)
c            rout = 10.d0**((dTlout * (rrl - rll) + rll*dlbi2 - 
c     &             rrl*dlbi1) / (dlbi2 - dlbi1))
            rlout = log10(rout)
c
c....       stellar mass
c 
            if ( amstot(i-1) .gt. 0. ) then 
              amsl1 = log10(amstot(i-1))
            else
              amsl1 = -15.0 
            endif
            if ( amstot(i) .gt. 0. ) then 
              amsl2 = log10(amstot(i))
            else
              amsl2 = -15.0
            endif
            ah = (amsl2 - amsl1) * rri
            bh = amsl1 - ah * rll 
            aM_st = 10.d0**(ah*rlout + bh)
c
c....       DM mass
c
            if ( amdmtot(i-1) .gt. 0.0 ) then 
              amdl1 = log10(amdmtot(i-1))
            else
              amdl1 = -15.0 
            endif
            if ( amdmtot(i) .gt. 0.0 ) then 
              amdl2 = log10( amdmtot(i) )
            else
              amdl2 = -15.0 
            endif
            ah = (amdl2 - amdl1) * rri
            bh = amdl1 - ah * rll 
            aM_dm   = 10.d0**(ah*rlout + bh)
c            write(*,*)'Mdm:',amdmtot(i-1)*aM0*hubble,
c     &        amdmtot(i)*aM0*hubble,aM_dm*aM0*hubble

            aM_st   = aM_st  * aM0 * hubble
            aM_dm   = aM_dm  * aM0 * hubble
            irflag  = 1 
          endif
c
c...  estimate mass if rmass is reached
c
          if ( (ibin .gt. 0 .and.
     &         rr .gt. rmass .and.
     &         rl .le. rmass) .and. irflag .eq. 0 ) then 
            ivir = i
            rrl = log10(rr)
            rll = log10(rl) 
            rri = 1.d0 / (rrl - rll)
            rout = rmass
            rlout = log10(rout)
c
c....       stellar mass
c 
            if ( amstot(i-1) .gt. 0. ) then 
              amsl1 = log10(amstot(i-1))
            else
              amsl1 = -15.0 
            endif
            if ( amstot(i) .gt. 0. ) then 
              amsl2 = log10(amstot(i))
            else
              amsl2 = -15.0
            endif
            ah = (amsl2 - amsl1) * rri
            bh = amsl1 - ah * rll 
            aM_st = 10.d0**(ah*rlout + bh)
c
c....       DM mass
c
            if ( amdmtot(i-1) .gt. 0.0 ) then 
              amdl1 = log10(amdmtot(i-1))
            else
              amdl1 = -15.0 
            endif
            if ( amdmtot(i) .gt. 0.0 ) then 
              amdl2 = log10( amdmtot(i) )
            else
              amdl2 = -15.0 
            endif
            ah = (amdl2 - amdl1) * rri
            bh = amdl1 - ah * rll 
            aM_dm   = 10.d0**(ah*rlout + bh)

            aM_gas  = aM_gas * aM0 * hubble
            aM_cgas = aM_cgas * aM0 * hubble
            aM_st   = aM_st  * aM0 * hubble
            aM_dm   = aM_dm  * aM0 * hubble
            irflag  = 1 
          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
        volr1 = volr 
      enddo

c
      xcd=(xc-1.0)*r0
      ycd=(yc-1.0)*r0
      zcd=(zc-1.0)*r0
      rfact = r0 * 1.e3 * aexpn / hubble
      vfact = v0 / aexpn
      rdout = rout * rfact
      Tfact = gamma1 * T_0 * a2i
c
c.... now estimate gas masses by summing over the cells
c
      aM_gas = 0.0
      aM_cgas = 0.0
c
      DO L = MinLevel+1 , MaxLevelNow
        CellVolume = CellSize(L)**3
        call Select_Cells ( L , nLevel )
        do ic1 = 1 , nLevel
          ic0 = iSelect(ic1)
          do ic2 = 0 , 7
            ic = ic0 + ic2
            if ( iOctCh(ic) .eq. 0 ) then
               call Ps(ic, Posx, Posy, Posz)
               xx = Posx
               yy = Posy
               zz = Posz 
               diff_x = xc - xx
               diff_y = yc - yy
               diff_z = zc - zz
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
c
c....          correct coordinate to insure periodic bundaries
c
               xx = xx + corr_x
               yz = yy + corr_y
               zz = zz + corr_z

               rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
               Tcell = Tfact * hvar(8,ic)/hvar(1,ic)

               if ( rp .le. rout ) then 
                 write(*,*) 'ic=',ic, Tcell
               if ( Tcell .lt. Tcold ) then 
                 aM_gas = aM_gas + hvar(1,ic)* CellVolume 
               endif
               endif
            endif
          enddo
        enddo
      ENDDO      
      aM_gas = aM_gas * aM0 
      
c
c.... averages
c
      if ( amgave .gt. 0.0 ) then 
        zave = zave / amgave ! already solar
      endif
      if ( amsave .gt. 0.0 ) then 
        zsave = 50.0 * zsave / amsave ! -> solar 
        td = tave / amsave
        tave = age ( 1.d0*td ) ! code units -> Gyr
      endif
      
      if ( ictype .eq. 1 ) then ! DM halos
        write(11,10) xcd, ycd, zcd, rdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc), nhgc(ihc)
        write(12,10) xcd, ycd, zcd, rdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc), nhgc(ihc)
        write(13,10) xcd, ycd, zcd, rdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc), nhgc(ihc)
      endif
      if ( ictype .eq. 0 ) then ! GC
        write(11,11) xcd, ycd, zcd, rdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave
        write(12,11) xcd, ycd, zcd, rdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave
        write(13,11) xcd, ycd, zcd, rdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm,
     &               zave, zsave, tave
      endif
 10   format (3(f9.6,1x),1x,f8.3,1x,6(e9.3,1x),f6.3,1x,e9.3,1x,i5)
 11   format (3(f9.6,1x),1x,f8.3,1x,6(e9.3,1x),f6.3)
c
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) 
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
        abtot = (amgtot(i)+ amstot(i))/(4.18879 * rr**3)
        if ( amst(i) .gt. 0.0 ) then ! correct for stellar metallicity
          rnp = 1.0 / amst(i)
          zbin(i) = (0.02*zbin(i)*dgbin(i)*volr) / volr /
     &               dgbin(i)*50.
          zsbin(i) = zsbin(i) * rnp * 50.
          zsbinIa(i) = zsbinIa(i) * rnp * 50.
          svx(i) = svx(i) * rnp
          svy(i) = svy(i) * rnp
          svz(i) = svz(i) * rnp
          vmn = svx(i)**2 + svy(i)**2 + svz(i)**2
          svrms(i) = sqrt(abs(svrms(i)*rnp - vmn)) * vfact
          svx(i) = svx(i) * vfact 
          svy(i) = svy(i) * vfact 
          svz(i) = svz(i) * vfact 
        endif
        if ( ddmbin(i) .gt. 0.0 ) then 
          dm = ddmbin(i) * volr 
          rnp = 1.0 / dm
          ppvx(i) = ppvx(i) * rnp
          ppvy(i) = ppvy(i) * rnp
          ppvz(i) = ppvz(i) * rnp
          vmn = ppvx(i)**2 + ppvy(i)**2 + ppvz(i)**2
          vvrms(i) = sqrt(abs(vvrms(i)*rnp - vmn)) * vfact
          ppvx(i) = ppvx(i) * vfact 
          ppvy(i) = ppvy(i) * vfact 
          ppvz(i) = ppvz(i) * vfact           
        endif
        aMdm = ddmibin(i) * (4.18879 * rr**3) * aM0 * hubble
        abtot = abtot*rhofact
        rhost = amst(i) / volr 
        hmassi = (amgtot(i)+amstot(i))*aM0 + aMdm/hubble
        rcirc = rr * r0/hubble * aexpn * 3.0856d24
        vcirc1 = 1.15206d8 * sqrt(hmassi/rcirc)
        hmassi = (amgtot(i)+amstot(i))*aM0
        vcirc2 = 1.15206d8 * sqrt(hmassi/rcirc)
        hmassi = (amstot(i))*aM0
        vcirc3 = 1.15206d8 * sqrt(hmassi/rcirc)
          write(12,18)
     &         rmid*rfact, rr*rfact, dtibin(i)*rhofact,
     &         dgbin(i)*rhofact,rhost*rhofact, ddmbin(i)*rhofact,
     &         amgtot(i)*aM0*hubble, 
     &         amgctot(i)*aM0*hubble, 
     &         amstot(i)*aM0*hubble,
     &         aMdm,tbin(i),zbin(i),zsbin(i),zsbinIa(i),
     &         int(ndmp(i)),int(nst(i))
          write(13,19)
     &         rmid*rfact, rr*rfact, 
     &         amgtot(i)*aM0*hubble,
     &         amstot(i)*aM0*hubble,
     &         aMdm, gvrms(i),svrms(i),vvrms(i), vcirc1, vcirc2, vcirc3,
     &         int(ndmp(i)),int(nst(i))
      enddo
 18   format(2(f8.3,1x),12(e9.3,1x),i6,1x,i6)
 19   format(2(f8.3,1x),1x,9(e9.3,1x),i6,1x,i6)
c
c
c.... now label cells falling within the radius of clump
c
      do i = 1 , npoints
        rr    = rout * RANDd ( m )
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
        ind(idcell) = 1
      enddo

      return
      end
