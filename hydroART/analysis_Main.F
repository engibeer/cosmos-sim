!     program Main_schedule
!     list of flags in schedule.dat

      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      real*8 dt
      character*60 title
      character*256 FileNameG,  fname_part1, fname_part2, fname_part3
     &     , fprofname
      character*256 base, key, a , ext
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
      real*8 tmax,tmin
C-----------------Schedule of analysis (flags= 0 Off , 1 On)-------- 
      read(*,10) title
 10   format (A60)
      write (*,10) title

      read (*,*) iflagNbody
      if ( iflagNbody.eq.1) write (*,*)  'Read N-Body file'
      read (*,*) iflagSTARS
      if (iflagSTARS.eq.1) write (*,*)  'Read stellar file'
      read (*,*) iflagHYDRO
      if ( iflagHYDRO.eq.1) write (*,*)  'Read HYDRO file'
      read (*,*) iflagSetU
      if (iflagSetU.eq.1) write (*,*)  'Set units'
      
      read (*,*) iflag0
      if ( iflag0.eq.1) write (*,*)  'Find center using particles'
      read (*,*) iflag0a
      if (iflag0a.eq.1) write (*,*) "Position of minimum potential"
      read (*,*) iflag1
      if ( iflag1.eq.1) write (*,*) 'spherical profiles'
      read (*,*) iflag2
      if ( iflag2.eq.1) write (*,*)  'Find center using gas'
      read (*,*) iflag3
      if ( iflag3.eq.1) write (*,*)' spherical velocity profiles'
      read (*,*) iflag4
      if ( iflag4.eq.1) write (*,*)  'Check range of variables.'
      read (*,*) iflag4b
      if(iflag4b.eq.1)write (*,*) 'Check ranges in code unit'
      read (*,*) iflag5
      if ( iflag5.eq.1) write (*,*)  'slice of hydro data'
      read (*,*) iflag6
      if (iflag6.eq.1)write (*,*)'ASCII of first specie of DM or stars'
      read(*,*) iflag7
      if (iflag7.eq.1) write (*,*) 'Compute the angular momentum & 
     & new system of coordinates' 
      read(*,*) iflag8
      if (iflag8.eq.1) write (*,*)'Cylindrical profiles.'
      read(*,*) iflag9
      if (iflag9.eq.1) write (*,*)  'Mapping of the cooling time'
      read(*,*) iflag10
      if (iflag10.eq.1) write (*,*) 'HALO FINDER'
      read(*,*) iflag11
      if (iflag11.eq.1) write (*,*)  'Stellar analysis'
      read(*,*) iflag12
      if (iflag12.eq.1) write (*,*) 'Distribution functions' 
      read(*,*) iflag13
      if (iflag13.eq.1) write (*,*) 'Run a QSO line of sight'
      read(*,*) iflag14
      if (iflag14.eq.1) write (*,*)  'Young stars studies'
      read(*,*) iflag15
      if (iflag15.eq.1) write (*,*)  'Mock distributions'
C----- Pre- analysis -------------------------------------------------------------
      read (*, 10) title
      write (*, 10) title
! Read files:
      read(*,20) jobname1
      write (*,20) jobname1
      read(*,20) FileNameG
      read(*, 20) fname_part1
      read(*, 20) fname_part2
      read(*, 20) fname_part3
 20   format (A256)
 21   format (('("Jobname= ",A25)'))
#ifdef HYDRO
      if ( iflagHYDRO.eq.1) then
         write (*,11)    FileNameG
         call Read_Gas_Binary ( dt , FileNameG )
      else
         write (*,11)    FileNameG
         call Read_Gas_Binary_HEADER ( dt , FileNameG )
      endif
#endif
#ifdef NBODY
      if ( iflagNbody.eq.1) then
         write (*,11)    fname_part1
         write (*,11)    fname_part2         
         call Read_Particles_Binary(fname_part1,fname_part2)
      endif
#endif
#ifdef STARFORM
      if (iflagSTARS.eq.1) then
         write (*,11)    fname_part3
         call   Read_Stellar_Data ( fname_part3)
      endif
#endif
      read (*,*) gamma
      if ( iflagSetU.eq.1)      call units() ! Set physical units
 11   format(  'Reading  file: ', A30)
c-------------------------------------------------------------------- 
c     PRE-ANALYSIS:
c-------------------------------------------------------------------- 
!     Find center of the main Halo:
      read (*,*) ifindcenter     
      rmin=1.  ! in proper Kpc
      rmax=5.
      rmin  =rmin / ( r0 * 1000.0 / hubble * aexpn) ! in code units
      rmax=rmax/ ( r0 * 1000.0 / hubble * aexpn)
      if (iflag0.eq.1) then     !! For particles
         write (*,*)  'Center position and velocity of particles:'
         write (*,*) '-------------------------------------------------'
         xc=(boxh/2.)/r0 ! Initial center in code units.
         yc=(boxh/2.)/r0
         zc=(boxh/2.)/r0
         if ( iflag0a.eq.1) then
            call Find_Halo_Center ( xc, yc, zc, ifindcenter )
            write (*,*) 'Position of minimum potential in  proper Kpc:'
            write (*,*)  (xc-1.)*1000.0*aexpn*r0/hubble
            write (*,*)  (yc-1.)*1000.0*aexpn*r0/hubble
            write (*,*)  (zc-1.)*1000.0*aexpn*r0/hubble
         else
        write (*,*)'Initial center suplied by user, see analysis_Main.F'
            xc=(0.014557 +1.)/2. * boxh / r0  +1 ! From ifrit ! Be careful with the flip of axis in ifrit!!
            yc=(-0.0968247 +1.)/2. * boxh / r0 +1
            zc=(0.0997198  +1.)/2. * boxh / r0  +1
!        xc = (6.07216   / r0 ) +1. ! From Halo Finder
!        yc = (4.81862  / r0 ) +1.
!        zc = (5.38863  / r0 ) +1.
            write (*,*)  (xc-1.)*1000.0*aexpn*r0/hubble
            write (*,*)  (yc-1.)*1000.0*aexpn*r0/hubble
            write (*,*)  (zc-1.)*1000.0*aexpn*r0/hubble
       endif
        call ParticlesCenter (rmin, rmax, Xc , Yc , Zc, Vxc, Vyc, Vzc)
!         call STARcenter (rmin, rmax, Xc , Yc , Zc, Vxc, Vyc, Vzc)
         write (*,*) 'Center in proper Kpc:'
         write (*,*)   (xc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)   (yc-1.)*1000.0*aexpn*r0/hubble 
         write (*,*)   (zc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)  'Velocity of the center in Km/s'
         write (*,*)  v0*Vxc / aexpn
         write (*,*)  v0*Vyc / aexpn
         write (*,*)  v0*Vzc / aexpn
!
         write (*,*)'Center position using gas inside rmin'
         call GasCenter(rmin, Xc , Yc , Zc)
         write (*,*) 'Center in proper Kpc:'
         write (*,*)   (xc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)   (yc-1.)*1000.0*aexpn*r0/hubble 
         write (*,*)   (zc-1.)*1000.0*aexpn*r0/hubble
!
      endif
      if (iflag2.eq.1) then     ! Center using gas
         write (*,*)'Center position and velocity using gas'
         write (*,*) '-------------------------------------------------'
         xc=(boxh/2.)/r0 ! Initial center in code units.
         yc=(boxh/2.)/r0
         zc=(boxh/2.)/r0
         call Find_Halo_Center ( xc, yc, zc, 0 )
         write (*,*) 'Center in proper Kpc using gas density:'
         write (*,*)   (xc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)   (yc-1.)*1000.0*aexpn*r0/hubble 
         write (*,*)   (zc-1.)*1000.0*aexpn*r0/hubble
         npoints = 100000
         rmax = 10. ! in proper Kpc  
         rmin = 1.
         rmin  =rmin / ( r0 * 1000.0 / hubble * aexpn) ! in code units
         rmax=rmax/ ( r0 * 1000.0 / hubble * aexpn)
         call Vcenter(npoints,xc,yc,zc,rmax , rmin, vxc, vyc, vzc)
          write (*,*) 'Center in proper Kpc:'
         write (*,*)   (xc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)   (yc-1.)*1000.0*aexpn*r0/hubble  
         write (*,*)   (zc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)  'Velocity of the center in Km/s'
         write (*,*)  v0*Vxc / aexpn
         write (*,*)  v0*Vyc / aexpn
         write (*,*)  v0*Vzc / aexpn        
      endif
!     Find total angular momentum:
      if (iflag7.eq.1) then
         isys=2 ! isys=1 --> DM particles,  isys=2 --> stars
         Rs=5. ! 60 ! proper Kpc
         Rs =Rs / ( r0 * 1000.0 / hubble * aexpn) ! in code units
         r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
         v0P = v0 / aexpn       ! in km/s (proper units)
         call angularM(Vxc,Vyc,Vzc,xc,yc,zc,isys,Rs,aLx,aLy,aLz)
         write (*,*)  'Angular momentum in Kpc Km/s (proper):'
         write (*,*)  aLx*aM0*r0P*v0P
         write (*,*)  aLy*aM0*r0P*v0P
         write (*,*)  aLz*aM0*r0P*v0P
         aL = sqrt ( aLx**2 + aLy**2 + aLz**2)
         write (*,*)  'AM axis'
         write (*,*)  aLx/aL
         write (*,*)  aLy/aL
         write (*,*)  aLz/aL
         call NewAxes(aLx, aLy, aLz, AxX,AxY, AxZ)
c         write (*,*) ' Vector 1= ', AXX(1) , AxX(2), AxX(3)
c         write (*,*)  'Vector 2= ', AXY(1), AXY(2), AXY(3)
c         write (*,*)  'Vector 3= ', AxZ(1), AxZ(2), AxZ(3)
      endif
c-------------------------------------------------------------------- 
c    Subroutines of analysis
C-------------------------------------------------------------------- 
      write (*,*) '! Subroutines of analysis:'
      read (*,10) title
      write (*,10) title
      read (*,*) rmin           ! min. radius of the profile
      read (*,*) rmax           ! max. radius of the profile
      read (*,*) nrbin          ! number of bins
      read (*,20) fprofname !  full path and file name to output profile to  it must end with a space
      read (*,*) icomponent ! 0 / 1 / 2 / 3 --> all / gas /stars / dm
      base= fprofname
      write (a, 124) aexpn
      ext = '.dat'
      lbase = index(base, ' ') -1
      la = index(a, ' ' ) -1
      lext = index (ext, ' ') -1
      rmin  =rmin / ( r0 * 1000.0 / hubble * aexpn)
      rmax=rmax/ ( r0 * 1000.0 / hubble * aexpn)
 124  format (F5.3)
      if (iflag1.eq.1) then     ! Spherical profiles of density, temperature,...
         T_gas =1.0d8
         key ='SpGe'
         lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext)
         call ProfilesSp (xc,yc,zc,T_gas,rmin, rmax,nrbin,fprofname)
      endif
      if (iflag3.eq.1) then     ! Velocity profiles:
         if ( icomponent.eq.1) then ! Gas velocity profile
            T_gas =1.0d8
            key ='SpVg'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext)
            call ProfilesVGa(Vxc,Vyc,Vzc,xc,yc,zc,
     &      T_gas,rmin,rmax,nrbin,fprofname)
         endif
         if ( icomponent.eq.2) then ! Stellar velocity profiles
            key ='SpVs'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
         call ProfilesVst(Vxc,Vyc,Vzc,xc,yc,zc,
     &        rmin,rmax,nrbin,fprofname) 
         endif
      if (icomponent.eq.3) then ! DM velocity component
            key ='SpVd'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
             call ProfilesVdm  (Vxc, Vyc, Vzc,xc,yc,zc,
     &        rmin, rmax,nrbin,fprofname) 
          endif
          if (icomponent.eq.0) then ! All components
            T_gas =1.0d8
             key ='SpVg'
             lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext)
             call ProfilesVGa(Vxc,Vyc,Vzc,xc,yc,zc,
     &     T_gas,rmin,rmax,nrbin,fprofname)
             key ='SpVs'
             lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
            call ProfilesVst(Vxc,Vyc,Vzc,xc,yc,zc,
     &     rmin,rmax,nrbin,fprofname) 
            key ='SpVd'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
            call ProfilesVdm  (Vxc, Vyc, Vzc,xc,yc,zc,
     &        rmin, rmax,nrbin,fprofname) 
         endif
      endif
!     call Xvir (xc,yc,zc,rmin,rmax)
!     deltavir = 180. 
!     lpointsout = .false.  
!     l_force_center = .true.  
!     call Output_Average_Profiles ( xc , yc , zc , rvir,
!     &                            rmin , rmax , nrbin , ibintype ,
!     &                            deltavir, ifindcenter ,
!     &                            lpointsout, l_force_center,
!     &                            fprofname )
c
c
      if (iflag4.eq.1) call MinMax(fprofname)
      if (iflag4b.eq.1) call MinMaxCU(fprofname)
      if (iflag5.eq.1) then 
c$$$         read (*,*) kslice
c$$$         write (*,*)  'cutting a slice of the zero-level mesh with k= ',
c$$$     &        kslice
         call slicing()
      endif
c      if (iflag6.eq.1) call particles_plot()
      if (iflag6.eq.1) then
         write (*,*) " writing a ASCII file with "
         if  (icomponent.eq.2) then
            tAge=2.d12 ! stellar age
            write (*,340) tAge
            call Write_Stars_ASCII(fprofname, Vxc, Vyc, Vzc,xc,yc,zc, 
     & rmax,tAge)
         else if (icomponent.eq.3) then
            write (*,*)  "first DM specie."
            call Write_DM_ASCII(fprofname,Vxc, Vyc, Vzc, xc,yc,zc,rmax) 
         else
            write (*,*)  'icomponent is wrong'
         write(*,*) " 2 / 3 = stars / DM"
         endif
 340  format ( 'stars younger than ', E6.1, ' yrs')
      endif
      if (iflag8.eq.1) then
         if (icomponent.eq.1) then 
            write (*,*)  'Cylindrical profiles for gas'
            ZDmax = 1.       ! in Proper Kpc
            T_gas =1.0d8
            ZDmax=ZDmax/ ( r0 * 1000.0 / hubble * aexpn) 
            key ='CyVg'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext)
            call CylinPG(Vxc,Vyc,Vzc,xc,yc,zc,
     &     zDmax,T_gas,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ)
         endif
         if ( icomponent.eq.2) then
            write (*,*)  'Cylindrical profiles for stars:'
            tAge=8.d10           ! stellar age            
            ZDmax = 1.        ! in Proper Kpc
            ZDmax=ZDmax/ ( r0 * 1000.0 / hubble * aexpn)
            key ='CyVs'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
            call CylinPS(Vxc,Vyc,Vzc,xc,yc,zc,
     &           ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ,tage)
         endif
         if (icomponent.eq.3) then
            write (*,*)  'Cylindrical profiles for DM:'
            ZDmax = 1.0         ! in Proper Kpc
            ZDmax=ZDmax/ ( r0 * 1000.0 / hubble * aexpn)
            key ='CyVd'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
            call CylinPD(Vxc,Vyc,Vzc,xc,yc,zc,
     &           ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ)
         endif
         if ( icomponent.eq.0 ) then ! All components:
            tAge=1.d8           ! stellar age
            ZDmax = 1.0         ! in Proper Kpc
            ZDmax=ZDmax/ ( r0 * 1000.0 / hubble * aexpn) 
            T_gas =1.0d4
            key ='CyVg'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext)
            call CylinPG(Vxc,Vyc,Vzc,xc,yc,zc,
     &     zDmax,T_gas,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ)
            key ='CyVs'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
            call CylinPS(Vxc,Vyc,Vzc,xc,yc,zc,
     &           ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ,tage)
            key ='CyVd'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
            call CylinPD(Vxc,Vyc,Vzc,xc,yc,zc,
     &           ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ)
         endif
         if ( icomponent.eq.4 ) then  ! Surface brightness profiles 
            tAge=2.d10           ! stellar age            
            ZDmax = 1.        ! in Proper Kpc
            ZDmax=ZDmax/ ( r0 * 1000.0 / hubble * aexpn)
            key ='CySb'
            lkey = index(key, ' ' ) -1
      fprofname=base(1:lbase)//key(1:lkey)//'a'//a(1:la)//ext(1:lext) 
            call CylinPSb(Vxc,Vyc,Vzc,xc,yc,zc,
     &           ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ,tage)
         endif
      endif
      if (iflag9.eq.1) then
         write (*,*)  'Setting cooling'
         call Set_Cooling()
         write (*,*)  'Cooling time map'
         call  cooling_T_map()
         call Write_Gas_Binary ( dt , fprofname )
      endif
      if (iflag10.eq.1) then
         write (*,*)  'HALO FINDER'
!         call  HaloFind(3) !  iWrite = 0/1/2/3 (No/DM/DM+Gas+Stars/Both)
      endif
      if (iflag11.eq.1) then
         write (*,*) 'STELLAR ANALYSIS'
         tmax=-2.0
         tmin=-4.1
         call  rate(tmin,tmax,nrbin, fprofname)
      endif
      if(iflag12.eq.1) then
         write (*,*)  'Making distribution function of...'
         ivar=1
         Vmin=0.0000001
         Vmax=1000000
         if (ivar.eq.1) then
            write (*,*)  'Gas density.'
         endif
         call distriFunc(ivar,Vmin,Vmax,nrbin,fprofname)
      endif
      if (iflag13.eq.1) then
         write (*,*)  'Setting cooling'
         call Set_Cooling()
         write (*,*)  'Running a line of sight through the simulation'
         QSOn=0.000
         base= fprofname
         ext = '.dat'
         lbase = index(base, ' ') -1
         lext = index (ext, ' ') -1
         do j=0,6
            DO I=0,6
               QSOn=QSOn + 0.01
               write(a,125) QSOn
               la = index(a, ' ' ) -1
 125           format(F4.2)
               fprofname=base(1:lbase)//a(1:la)//ext(1:lext) 
         QSOb1 = -30.0 + I*10.0              ! Impact parameter in proper Kpc
         QSOb2 = -30.0 + j*10.0
         write(*,*) 'Shooting: ', fprofname,'(',QSOb1,QSOb2,')'
         QSOb1 = QSOb1  / ( r0 * 1000.0 / hubble * aexpn) ! in code units
         QSOb2 = QSOb2  / ( r0 * 1000.0 / hubble * aexpn) ! in code units
         IF (iflag7.EQ.0) THEN
            WRITE (*,*)  'The principal axis should be computed'
            write (*,*)  'Stopping...See schedule.dat'
         ENDIF
         QSODx = AXz(1)         ! Given Direction
         QSODy = AXz(2)
         QSODz = AXz(3)
c
         QSObx =QSOb1 * AXx(1) + QSOb2 * AXy(1) ! Given Impact parameter vector in code units
         QSOby =QSOb1 * AXx(2) + QSOb2 * AXy(2)
         QSObz =QSOb1 * AXx(3) + QSOb2 * AXy(3)
c
c         write (*,*)  'Using galactic center in code units:'
c         write (*,*)  xc, yc, zc
         write (*,*)  'Impact parameter vector in code units'
         write (*,*)  QSObx , QSOby, QSObz
         write (*,*)  'Direction:'
         write (*,*)  QSODx, QSODy,QSODz
c
         QSObx = xc + QSObx
         QSOby = yc + QSOby
         QSObz = zc + QSObz
c
!         QSOx0 = 64. ! Given starting point
!         QSOy0 = 64.
!         QSOz0 = 0. ! or...
         call Find_Starting_P(QSObx,QSOby,QSObz,QSODx,QSODy,QSODz,
     &        QSOx0,QSOy0,QSOz0,Ng) 
         write (*,*)  'Starting point:'
         write (*,*)  QSOx0, QSOy0, QSOz0
         call  QSOline(xc,yc,zc,Vxc,Vyc,Vzc,QSOb1, QSOb2,
     &        QSOx0,QSOy0,QSOz0,QSODx,QSODy,QSODz,fprofname)
         enddo
      ENDDO
      endif
      if (iflag14.eq.1) then
C         call Stellar_data(fprofname)
C         call YoungSt() 
         call SF_law(fprofname)
c         call Write_Gas_Binary ( dt , fprofname )
      endif
      if (iflag15.eq.1) then
         ioptMock=2             ! 1--> Density distribution 1/r^2
                                ! 2-->  Constant radial flow
         write (*,*)  ioptMock,xc,yc,zc,rmax,Vxc,Vyc,Vzc
         call Mock_distributions(ioptMock,xc,yc,zc,rmax,Vxc,Vyc,Vzc)
         call Write_Gas_Binary ( dt , fprofname )
      endif
      end

