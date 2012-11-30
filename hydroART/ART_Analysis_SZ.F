c...  Analysis routines (SZ effects)
c        Analyze_SZ
c        Output_SZ_3D_Profiles
c        Output_SZ_2D_Profiles
c        Output_SZvir
c
c     ---------------------------------------
      subroutine Analyze_SZ ( if1, if2, if3 )
c     ---------------------------------------
c     
c     Driver for computing SZ observables
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'

      character*256 fname
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      logical l_force_center

c
c.... Compute SZ observables 
c
      rmin = 1.0e-3/r0 ! 1h^-1 kpc
      rmax = 10.0/r0   ! 10h^-1 Mpc
      nrbin = 80
      ibintype = 0 
      deltavir = 180.
      deltacrit = 500.
      ifindcenter = 3 
      l_force_center = .false.  
      nstep = 0 
      lpath = index(path    , ' ') - 1
      
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)
c
c...  compute 3D SZ profiles by directly projecting down values of each cell
c
      if ( if1 .eq. 1 ) then
         fname  = path(1:lpath)//'/pro_sz3d'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

         write(*,*) 'Output_SZ_3D_Profiles'
         call Output_SZ_3D_Profiles ( 
     &                      xc , yc , zc , rvir, rcrit,  
     &                      rmin , rmax , nrbin , ibintype ,  
     &                      deltavir, deltacrit, ifindcenter , 
     &                      lpointsout, l_force_center,
     &                      fname, rs_fname ) 
      endif

c
c...  compute 2D SZ profiles by directly projecting down values of each cell
c     
      if ( if2 .eq. 1 ) then
         fname  = path(1:lpath)//'/pro_sz2d'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '

         write(*,*) 'Output_SZ_2D_Profiles'
         call Output_SZ_2D_Profiles ( xc, yc, zc, 
     &        deltavir, deltacrit, 
     &        rmin , rmax , nrbin , ibintype, 
     &        l_force_center, ifindcenter, fname )
      endif

c
c...  compute SZ properties at virial radius
c     
      if ( if3 .eq. 1 ) then
         fname  = path(1:lpath)//'/szvir'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
         write(*,*) 'Output_SZvir'
         call Output_SZvir ( ifindcenter, fname )
      endif

      return
      end

c
c     --------------------------------------------------------------
      subroutine Output_SZ_3D_Profiles ( 
     &                            xc , yc , zc , rvir, rcrit, 
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, deltacrit, ifindcenter, 
     &                            lpointsout, l_force_center,
     &                            fprofname, rs_fname )
c     --------------------------------------------------------------
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
      real xc , yc , zc , rvir, rcrit, rmin , rmax 
      real deltavir, deltacrit
      integer nrbin, ibintype, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname, rs_fname
c
      character*256 fname
c 
      real*8 a2b, b2a
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pe(npoints)
      real pzIa(npoints), pzII(npoints)
      real ptemp(npoints), pent(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pe, pzIa, pzII, 
     +     ptemp, pent
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dgibin(0:nmaxbins),dgbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins),dtibin(0:nmaxbins)
      real sztbin(0:nmaxbins),sztibin(0:nmaxbins),tmbin(0:nmaxbins)
      real tmibin(0:nmaxbins)
      real ndmp(0:nmaxbins)
      real szt, szk_x, szk_y, szk_z
      real tdum(8), entdum(8)
      integer nbin(nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz, exbin(0:nmaxbins)
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn
c
      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) + 1
c
      a5 = aexpn**5
      a4 = aexpn**4
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      a3i = 1.0 / a3
      a4i = 1.0 / a4
      a5i = 1.0 / a5
      gamma1 = gamma - 1.0 
c
      Ynow = 0.2496
      zsim = 1.0/aexpn - 1.0
      if ( zsim .le. 1.0e-20 ) zsim = 0.0
      
c...  factors for thermal and kinematic SZ effects
      szt_fact = 3.41262e-15 * r0 * T_0 * hubble * Om0 * a4i
      szk_fact = 6.75012e-11 * r0 * v0 * hubble * Om0 * a3i

      call Get_MaxLevelNow ()

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

c....   gas density
        pd(i) = hvar(1,idcell)

c....   internal energy
        pe(i) = hvar(8,idcell)

c....   metallicity
        pzIa(i) = hvar(izIa,idcell)
        pzII(i) = hvar(izII,idcell)

      enddo
c
c.... construct profiles
c
      do i = 0 , nmaxbins
        rbin(i) = 0.
        dgbin(i)  = 0.
        dgibin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i)  = 0.
        sztbin(i)  = 0.
        sztibin(i) = 0.
	tmbin(i)   = 0.
        tmibin(i)  = 0.
        nbin(i) = 0
        ndmp(i) = 0.
      enddo

      do i = 1 , np
        xx = x(i)
        yy = y(i)
        zz = z(i)
        ! correct coordinate to insure periodic bundaries
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          ddmbin(ibin) = ddmbin(ibin) + pw(i)
          ndmp(ibin) = ndmp(ibin) + 1
        endif
      enddo

c
c.... compute SZ observables
c
      if ( lpointsout ) then 
        open ( 12 , file = 'points.dat' )
      endif
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )

c...    compute temperature and metallicity of the cell
        if ( pd(ii) .gt. 0. ) then
           Tcell = gamma1 * pe(ii) / pd(ii)
           Tcell = T_0 * Tcell * a2i
           Zcell = pzIa(ii) + pzII(ii)
           Zcell = Zcell / 0.01989 / pd(ii)
        else
           write(*,*) 'bad gas density in cell', ii,pd(ii)
        endif

c...    compute TSZ & KSZ signals of the cell
        szt = szt_fact * gamma1 * pe(ii) 
        szk_x = szk_fact * hvar(3,ii) 
        szk_y = szk_fact * hvar(4,ii) 
        szk_z = szk_fact * hvar(5,ii) 
 
        if ( ibin .le. nmaxbins ) then    
          nbin(ibin)   = nbin(ibin)  + 1
 	  rbin(ibin)   = rbin(ibin)  + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
c...      add up values  & convert code units to physical units
          sztbin(ibin) = sztbin(ibin) + szt
          tmbin(ibin)  = tmbin(ibin) + Tcell * pd(ii)
          if ( lpointsout ) then 
            write(12,*) rp*r0, pd(ii), Tcell, exbin(ii)
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
          tmibin(i) = tmbin(i) / nbin(i)
          tmbin(i)  = tmbin(i) / dgbin(i)
          dgbin(i)  = dgbin(i) / nbin(i)
c          ddmbin(i) = ddmbin(i) / nbin(i)
          sztbin(i) = sztbin(i) / nbin(i)
	endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dgi = 0.0
      dmi = 0.0
      dti = 0.0
      szti = 0.0
      tmi = 0.0
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
          dgi = dgi + dgbin(i)*volr
          dmi = dmi + ddmbin(i) 
          dti = dti + ddmbin(i) + dgbin(i)*volr 
          ddmbin(i) = ddmbin(i)/volr 
          szti = szti + sztbin(i)*volr
          tmi  = tmi  + tmibin(i)*volr
          dgibin(i)  = dgi / (4.18879 * rr**3)
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i)  = dti / (4.18879 * rr**3)
          sztibin(i) = szti / (4.18879 * rr**3)
          tmibin(i)  = tmi / dgi
c          write(*,*) rr,dmi,dti,ddmbin(i),ddmibin(i),dtibin(i),volr
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
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
c
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2
      aMcrit = dcritave * 
     &        4.18879 * (rcrit*r0)**3 * rho0 / hubble**2 
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype
      write(12,12) ifindcenter, l_force_center
      write(12,13) 
      write(12,15) deltavir, ivir, rvir*r0, aMvir
      write(12,14) 
      write(12,15) deltacrit, icrit, rcrit*r0, aMcrit
      write(12,17)
      write(12,18)
      write(12,19)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# ifindcenter =',i2,' l_force_center=',l6)
 13   format ('# Dvir(mean)    ivir     Rvir [1/h Mpc]   ',
     &     ' Mvir [1/h Msun]')
 14   format ('# Dcrit(crit)   icrit    Rcrit [1/h Mpc]  ',
     &     ' Mcrit [1/h Msun]' )
 15   format (f9.2,6x,i4,6x,f7.4,4x,g13.5) 
 17   format ('# Columns:')
 18   format ('# <r> r_l r_m r_r delta_g deltai_g delta_dm deltai_dm',
     &        ' deltai_tot  szt  szti  tm  tmi  N_dm_p' )
 19   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '          erg/s/cm^2/Sr   Kelvin  Zsun  ')
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
          write(12,20)
     &         rbin(i), rl*r0, rmid, rr*r0, 
     &         dgbin(i),dgibin(i),ddmbin(i),ddmibin(i),dtibin(i), 
     &         sztbin(i),sztibin(i),tmbin(i),tmibin(i),
     &         int(ndmp(i))
	endif
      enddo
 20   format(14(2x,g14.7),1x,i7)
 21   format(8(2x,g14.7))
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end

c
c     -------------------------------------------------
      subroutine Output_SZ_2D_Profiles ( xc, yc, zc, 
     &        deltavir, deltacrit,  
     &        rmin , rmax , nrbin , ibintype, 
     &        l_force_center, ifindcenter, fprofname )
c     -------------------------------------------------
c     compute XSBP, Tew, Zew profiles for an object centered 
c     at (xc,yc,zc) and radius rvir
c     all quantities are in code units
c
c     input : xc, yc, zc, rvir, Tmax
c             integer ibintype : 0 - logarithmic, 1 - even (now only log)
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
c     output: mass of gas at T<Tmax, total mass of gas within rvir (masses are in code units)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real   xc, yc, zc, rcore, rmin, rmax
      double precision cnt_rate
      integer ibintype
      character*256 fprofname
c
      character*256 fname
      parameter ( nmaxbins = 1000 )
      real rxbin(0:nmaxbins),rybin(0:nmaxbins),rzbin(0:nmaxbins),
     +  szt_xbin(0:nmaxbins),szt_ybin(0:nmaxbins),
     +  szt_zbin(0:nmaxbins),szk_xbin(0:nmaxbins),
     +  szk_ybin(0:nmaxbins),szk_zbin(0:nmaxbins)
      real szti_xbin(0:nmaxbins),szti_ybin(0:nmaxbins),
     +  szti_zbin(0:nmaxbins),szki_xbin(0:nmaxbins),
     +  szki_ybin(0:nmaxbins),szki_zbin(0:nmaxbins)
      real nxbin(0:nmaxbins),nybin(0:nmaxbins),nzbin(0:nmaxbins)
      integer nfpn
c
      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 
      nbins = int((rlmax-rlmin)/drl) + 1
c
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      call Find_Halo_Radius ( xc, yc, zc, deltavir, deltacrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )

      write(*,*)
      write(*,20) xc, yc, zc
      write(*,21)
      write(*,23) deltavir, rvir*r0, aMvir, aMvir_dm, aMvir_gas, 
     &     aMvir_st
      write(*,22)
      write(*,23) deltacrit, rcrit*r0, aMcrit, aMcrit_dm, aMcrit_gas, 
     &     aMcrit_st

 20   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 21   format ('# Dvir(mean)  Rvir[1/hMpc]    Mvir       Mdm  ',
     &     '       Mgas        Mst  [1/h Msun]' )
 22   format ('# Dcrit(crit) Rvir[1/h Mpc]   Mvir       Mdm  ',
     &     '       Mgas        Mst  [1/h Msun]' )
 23   format (f9.2,5x,f7.3,3x,4(g12.4,1x))

c
      a5 = aexpn**5
      a4 = aexpn**4
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      a3i = 1.0 / a3
      a4i = 1.0 / a4
      a5i = 1.0 / a5
      gamma1 = gamma - 1.0 
c
      Ynow = 0.2496
      zsim = 1.0/aexpn - 1.0
      if ( zsim .le. 1.0e-20 ) zsim = 0.0
      
c...  factors for thermal and kinematic SZ effects
      szt_fact = 3.41262e-15 * r0 * T_0 * hubble * Om0 * a4i
      szk_fact = 6.75012e-11 * r0 * v0 * hubble * Om0 * a3i

      call Get_MaxLevelNow ()

c.... construct profiles
      do i = 0 , nmaxbins
        rxbin(i) = 0.
        szt_xbin(i) = 0.
	szk_xbin(i) = 0.
        szti_xbin(i) = 0.
        szki_xbin(i) = 0.
        nxbin(i) = 0.

        rybin(i) = 0.
        szt_ybin(i) = 0.
	szk_ybin(i) = 0.
        szti_ybin(i) = 0.
        szki_ybin(i) = 0.
        nybin(i) = 0.

        rzbin(i) = 0.
        szt_zbin(i) = 0.
	szk_zbin(i) = 0.
        szti_zbin(i) = 0.
        szki_zbin(i) = 0.
        nzbin(i) = 0.
      enddo

      szti_x = 0.
      szti_y = 0.
      szti_z = 0.
      szki_x = 0.
      szki_y = 0.
      szki_z = 0.

c...  Loop over pixel-by-pixel
      DO Level = MinLevel , MaxLevelNow
        CellVolume = 1.0 * 2.0**(-3.0*Level)
        CellArea = 1.0 * 2.0**(-2.0*Level)
        cs = CellSize(Level)
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
                rx = sqrt( (Posy-yc)**2 + (Posz-zc)**2 )
                ry = sqrt( (Posx-xc)**2 + (Posz-zc)**2 )
                rz = sqrt( (Posx-xc)**2 + (Posy-yc)**2 )
      
                ixbin = max ( int((log10(rx) - rlmin)/drl) + 1 , 0 )
                iybin = max ( int((log10(ry) - rlmin)/drl) + 1 , 0 )
                izbin = max ( int((log10(rz) - rlmin)/drl) + 1 , 0 )
      
c...            compute TSZ & KSZ signals of the cell
                szt = szt_fact * gamma1 * hvar(8,ic) * CellVolume
                szk_x = szk_fact * hvar(3,ic) * CellVolume
                szk_y = szk_fact * hvar(4,ic) * CellVolume
                szk_z = szk_fact * hvar(5,ic) * CellVolume
 
c...            projection along x-axis
                if ( ixbin .le. nrbin ) then    
                   nxbin(ixbin) = nxbin(ixbin) + CellVolume
                   rxbin(ixbin) = rxbin(ixbin) + rx * CellVolume
                   szt_xbin(ixbin) = szt_xbin(ixbin) + szt
                   szk_xbin(ixbin) = szk_xbin(ixbin) + szk_x
                endif

                if ( rx .le. rvir ) then
                   szti_x = szti_x + szt
                   szki_x = szki_x + szk_x
                endif
                
                if ( rx .le. rcrit ) then
                   szti_xc = szti_xc + szt
                   szki_xc = szki_xc + szk_x
                endif
      
c...            projection along y-axis
                if ( iybin .le. nrbin ) then    
                   nybin(iybin) = nybin(iybin) + CellVolume
                   rybin(iybin) = rybin(iybin) + ry * CellVolume
                   szt_ybin(iybin) = szt_ybin(iybin) + szt 
                   szk_ybin(iybin) = szk_ybin(iybin) + szk_y 
                endif
      
                if ( ry .le. rvir ) then
                   szti_y = szti_y + szt
                   szki_y = szki_y + szk_y
                endif

                if ( ry .le. rcrit ) then
                   szti_yc = szti_yc + szt
                   szki_yc = szki_yc + szk_y
                endif

c...            projection along z-axis
                if ( izbin .le. nrbin ) then    
                   nzbin(izbin) = nzbin(izbin) + CellVolume
                   rzbin(izbin) = rzbin(izbin) + rz * CellVolume
                   szt_zbin(izbin) = szt_zbin(izbin) + szt
                   szk_zbin(izbin) = szk_zbin(izbin) + szk_z
                endif
      
                if ( rz .le. rvir ) then
                   szti_z = szti_z + szt
                   szki_z = szki_z + szk_z
                endif

                if ( rz .le. rcrit ) then
                   szti_zc = szti_zc + szt
                   szki_zc = szki_zc + szk_z
                endif

              endif
            enddo
          enddo
        ENDIF
      ENDDO

c
c...  compute the integrated SZ profiles
c
      
      szti_xbin(0) = szti_xbin(0)
      szki_xbin(0) = szki_xbin(0)
      szti_ybin(0) = szti_ybin(0)
      szki_ybin(0) = szki_ybin(0)
      szti_zbin(0) = szti_zbin(0)
      szki_zbin(0) = szki_zbin(0)
      
      do i = 1, nrbin
         szti_xbin(i) = szti_xbin(i-1) + szt_xbin(i)
         szki_xbin(i) = szki_xbin(i-1) + szk_xbin(i)
         szti_ybin(i) = szti_ybin(i-1) + szt_ybin(i)
         szki_ybin(i) = szki_ybin(i-1) + szk_ybin(i)
         szti_zbin(i) = szti_zbin(i-1) + szt_zbin(i)
         szki_zbin(i) = szki_zbin(i-1) + szk_zbin(i)
      enddo

c
c.... prepare profiles
c     
      do i = 0 , nrbin
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlmin + float(i-1)*drl)
         endif
         rr = 10.**(rlmin + float(i)*drl)
         area = pi * (rr**2 - rl**2)  
         area2 = pi * rr**2

c...     x-axis
         if ( nxbin(i) .gt. 0. ) then
            rxbin(i) = rxbin(i) / nxbin(i) * r0
            szt_xbin(i) = szt_xbin(i) / area 
            szk_xbin(i) = szk_xbin(i) / area
            szti_xbin(i) = szti_xbin(i) / area2 
            szki_xbin(i) = szki_xbin(i) / area2
         endif
      
c...     y-axis
         if ( nybin(i) .gt. 0. ) then
            rybin(i) = rybin(i) / nybin(i) * r0
            szt_ybin(i) = szt_ybin(i) / area 
            szk_ybin(i) = szk_ybin(i) / area
            szti_ybin(i) = szti_ybin(i) / area2 
            szki_ybin(i) = szki_ybin(i) / area2
         endif
      
c...     z-axis
         if ( nzbin(i) .gt. 0. ) then
            rzbin(i) = rzbin(i) / nzbin(i) * r0
            szt_zbin(i) = szt_zbin(i) / area 
            szk_zbin(i) = szk_zbin(i) / area
            szti_zbin(i) = szti_zbin(i) / area2 
            szki_zbin(i) = szki_zbin(i) / area2
         endif            
      enddo


c
c.... prepare profiles
c     
      area_vir  = pi * rvir**2  
      area_crit = pi * rcrit**2

      szti_x = szti_x / area_vir
      szki_x = szki_x / area_vir
      szti_y = szti_y / area_vir
      szki_y = szki_y / area_vir
      szti_z = szti_z / area_vir
      szki_z = szki_z / area_vir
      
      szti_xc = szti_xc / area_crit
      szki_xc = szki_xc / area_crit
      szti_yc = szti_yc / area_crit
      szki_yc = szki_yc / area_crit
      szti_zc = szti_zc / area_crit
      szki_zc = szki_zc / area_crit

c...  x-axis
      if ( nxbin(i) .gt. 0. ) then
         rxbin(i) = rxbin(i) / nxbin(i) * r0
         szt_xbin(i) = szt_xbin(i) / area 
         szk_xbin(i) = szk_xbin(i) / area
         szti_xbin(i) = szti_xbin(i) / area2 
         szki_xbin(i) = szki_xbin(i) / area2
      endif
      
c...  y-axis
      if ( nybin(i) .gt. 0. ) then
         rybin(i) = rybin(i) / nybin(i) * r0
         szt_ybin(i) = szt_ybin(i) / area 
         szk_ybin(i) = szk_ybin(i) / area
         szti_ybin(i) = szti_ybin(i) / area2 
         szki_ybin(i) = szki_ybin(i) / area2
      endif
      
c...  z-axis
      if ( nzbin(i) .gt. 0. ) then
         rzbin(i) = rzbin(i) / nzbin(i) * r0
         szt_zbin(i) = szt_zbin(i) / area 
         szk_zbin(i) = szk_zbin(i) / area
         szti_zbin(i) = szti_zbin(i) / area2 
         szki_zbin(i) = szki_zbin(i) / area2
      endif            
      

c...  Writing data to the file
      open ( 12 , file = fprofname(1:nfpn) )
c
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype 
      write(12,12)
      write(12,14) deltavir, rvir*r0, szti_x, szki_x  
      write(12,14) deltavir, rvir*r0, szti_y, szki_y
      write(12,14) deltavir, rvir*r0, szti_z, szki_z
      write(12,13)
      write(12,14) deltacrit, rcrit*r0, szti_xc, szki_xc   
      write(12,14) deltacrit, rcrit*r0, szti_yc, szki_yc  
      write(12,14) deltacrit, rcrit*r0, szti_zc, szki_zc 
      write(12,15) 
      write(12,16)
      write(12,17)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# Dvir(mean)  Rvir [1/h Mpc]   <y>(<Rvir)   <b>(<Rvir)')
 13   format ('# Dcrit(crit) Rcrit [1/h Mpc]  <b>(<Rcrit)  <b>(<Rcrit)')
 14   format (f9.2,4x,f7.4,3x,g13.5,3x,g13.5)
 15   format ('# Columns:')
 16   format ('# r_l r_m r_r  r_x szt_x szk_x szti_x szki_x ',
     &     ' r_y szt_y szk_y szti_y szki_y  ',
     &     ' r_z szt_z szk_z szti_z szki_z ')
 17   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens. ')
      do i = 0 , nrbin
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) * r0
        endif
        rr = 10.**(rlmin + float(i)*drl)
        if ( nxbin(i) .gt. 0. .and. 
     +       nybin(i) .gt. 0. .and. 
     +       nzbin(i) .gt. 0.  ) then
          write(12,18)
     &       rl*r0, rmid, rr*r0, 
     &       rxbin(i),szt_xbin(i),szk_xbin(i),szti_xbin(i),szki_xbin(i),
     &       rybin(i),szt_ybin(i),szk_ybin(i),szti_ybin(i),szki_ybin(i),
     &       rzbin(i),szt_zbin(i),szk_zbin(i),szti_zbin(i),szki_zbin(i)
	endif
      enddo
 18   format(18(2x,g14.7))
c 18   format(7(2x,g14.7),1x,i7,4(2x,g14.7),1x,i7,4(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      
      return
      end
c
c     --------------------------------------------------
      subroutine Output_SZvir ( ifindcenter, fprofname )
c     --------------------------------------------------
c     compute mean tSZ & kSZ enclosed within an aperture defined by
c     variousvirial radius definitions.
c     all quantities are in code units
c
c     input : xc, yc, zc, rvir, Tmax
c             integer ibintype : 0 - logarithmic, 1 - even (now only log)
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
c     output: mass of gas at T<Tmax, total mass of gas within rvir (masses are in code units)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      character*256 fprofname
      integer nfpn
c
      a5 = aexpn**5
      a4 = aexpn**4
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      a3i = 1.0 / a3
      a4i = 1.0 / a4
      a5i = 1.0 / a5
      gamma1 = gamma - 1.0 
c
      Ynow = 0.2496
      zsim = 1.0/aexpn - 1.0
      if ( zsim .le. 1.0e-20 ) zsim = 0.0
      
c...  factors for thermal and kinematic SZ effects
      szt_fact = 3.41262e-15 * r0 * T_0 * hubble * Om0 * a4i
      szk_fact = 6.75012e-11 * r0 * v0 * hubble * Om0 * a3i
c
      call Get_MaxLevelNow ()
c
      call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      
      oma = Om0 / ( Om0 + Oml0 * aexpn**3 )
      xf = oma - 1.
      Deltacrit = ( 18.*pi**2 + 82.*xf - 39.*xf**2 )
      Deltavir = Deltacrit / (1.+xf)
      write(*,*) aexpn, Om0, Oml0, oma, xf, Deltavir, Deltacrit

      call Find_Halo_Radius ( xc, yc, zc, Deltavir, Deltacrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )

      call Find_Halo_Radius ( xc, yc, zc, 180., 180., 
     &     i180, r180, aM180, aM180_dm, aM180_gas, aM180_st, 
     &     i180c, r180c, aM180c, aM180c_dm, aM180c_gas, aM180c_st, 
     &     ifindcenter )

      call Find_Halo_Radius ( xc, yc, zc, 200., 200., 
     &     i200, r200, aM200, aM200_dm, aM200_gas, aM200_st, 
     &     i200c, r200c, aM200c, aM200c_dm, aM200c_gas, aM200c_st, 
     &     ifindcenter )

      call Find_Halo_Radius ( xc, yc, zc, 500., 500., 
     &     i500, r500, aM500, aM500_dm, aM500_gas, aM500_st, 
     &     i500c, r500c, aM500c, aM500c_dm, aM500c_gas, aM500c_st, 
     &     ifindcenter )

      call Find_Halo_Radius ( xc, yc, zc, 2500., 2500., 
     &     i2500, r2500, aM2500, aM2500_dm, aM2500_gas, aM2500_st, 
     &     i2500c, r2500c, aM2500c, aM2500c_dm, aM2500c_gas, aM2500c_st, 
     &     ifindcenter )

      call Find_Halo_Radius ( xc, yc, zc, 5000., 5000., 
     &     i5000, r5000, aM5000, aM5000_dm, aM5000_gas, aM5000_st, 
     &     i5000c, r5000c, aM5000c, aM5000c_dm, aM5000c_gas, aM5000c_st, 
     &     ifindcenter )
c
      area_vir  = pi * rvir**2
      area_crit = pi * rcrit**2
      area_180  = pi * r180**2  
      area_180c = pi * r180c**2
      area_200  = pi * r200**2  
      area_200c = pi * r200c**2
      area_500  = pi * r500**2  
      area_500c = pi * r500c**2
      area_2500  = pi * r2500**2  
      area_2500c = pi * r2500c**2
      area_5000  = pi * r5000**2  
      area_5000c = pi * r5000c**2
c
      szti_x = 0.
      szti_y = 0.
      szti_z = 0.
      szki_x = 0.
      szki_y = 0.
      szki_z = 0.

c...  Loop over pixel-by-pixel
      DO Level = MinLevel , MaxLevelNow
        CellVolume = 1.0 * 2.0**(-3.0*Level)
        CellArea = 1.0 * 2.0**(-2.0*Level)
        cs = CellSize(Level)
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
                rr = sqrt( (Posx-xc)**2 + (Posy-yc)**2 + (Posz-zc)**2 )
                rx = sqrt( (Posy-yc)**2 + (Posz-zc)**2 )
                ry = sqrt( (Posx-xc)**2 + (Posz-zc)**2 )
                rz = sqrt( (Posx-xc)**2 + (Posy-yc)**2 )
      
c                ixbin = max ( int((log10(rx) - rlmin)/drl) + 1 , 0 )
c                iybin = max ( int((log10(ry) - rlmin)/drl) + 1 , 0 )
c                izbin = max ( int((log10(rz) - rlmin)/drl) + 1 , 0 )
      
c...            compute TSZ & KSZ signals of the cell
                szt = szt_fact * gamma1 * hvar(8,ic) * CellVolume
                szk_x = szk_fact * hvar(3,ic) * CellVolume
                szk_y = szk_fact * hvar(4,ic) * CellVolume
                szk_z = szk_fact * hvar(5,ic) * CellVolume

c...            SZ signals from the cluster virial regions
                if ( rr .le. rvir ) then
                   sztvir_r = sztvir_r + szt / area_vir
                endif
                if ( rr .le. r180 ) then
                   szt180_r = szt180_r + szt / area_180
                endif
                if ( rr .le. r200 ) then
                   szt200_r = szt200_r + szt / area_200
                endif
                if ( rr .le. r500 ) then
                   szt500_r = szt500_r + szt / area_500
                endif
                if ( rr .le. r2500 ) then
                   szt2500_r = szt2500_r + szt / area_2500
                endif
                if ( rr .le. r5000 ) then
                   szt5000_r = szt5000_r + szt / area_5000
                endif
                
                if ( rr .le. rcrit ) then
                   sztcrit_r = sztcrit_r + szt / area_crit
                endif
                if ( rr .le. r180c ) then
                   szt180c_r = szt180c_r + szt / area_180c
                endif
                if ( rr .le. r200c ) then
                   szt200c_r = szt200c_r + szt / area_200c
                endif
                if ( rr .le. r500c ) then
                   szt500c_r = szt500c_r + szt / area_500c 
                endif
                if ( rr .le. r2500c ) then
                   szt2500c_r = szt2500c_r + szt / area_2500c
                endif       
                if ( rr .le.r5000c ) then
                   szt5000c_r = szt5000c_r + szt / area_5000c
                endif

c...            projection along x-axis
                if ( rx .le. rvir ) then
                   sztvir_x = sztvir_x + szt / area_vir
                   szkvir_x = szkvir_x + szk_x / area_vir
                endif
                if ( rx .le. r180 ) then
                   szt180_x = szt180_x + szt / area_180
                   szk180_x = szk180_x + szk_x / area_180
                endif
                if ( rx .le. r200 ) then
                   szt200_x = szt200_x + szt / area_200
                   szk200_x = szk200_x + szk_x / area_200
                endif
                if ( rx .le. r500 ) then
                   szt500_x = szt500_x + szt / area_500
                   szk500_x = szk500_x + szk_x / area_500
                endif
                if ( rx .le. r2500 ) then
                   szt2500_x = szt2500_x + szt / area_2500
                   szk2500_x = szk2500_x + szk_x / area_2500
                endif
                if ( rx .le. r5000 ) then
                   szt5000_x = szt5000_x + szt / area_5000
                   szk5000_x = szk5000_x + szk_x / area_5000
                endif

                if ( rx .le. rcrit ) then
                   sztcrit_x = sztcrit_x + szt / area_crit
                   szkcrit_x = szkcrit_x + szk_x / area_crit
                endif
                if ( rx .le. r180c ) then
                   szt180c_x = szt180c_x + szt / area_180c
                   szk180c_x = szk180c_x + szk_x / area_180c
                endif
                if ( rx .le. r200c ) then
                   szt200c_x = szt200c_x + szt / area_200c
                   szk200c_x = szk200c_x + szk_x / area_200c
                endif
                if ( rx .le. r500c ) then
                   szt500c_x = szt500c_x + szt / area_500c 
                   szk500c_x = szk500c_x + szk_x / area_500c
                endif
                if ( rx .le. r2500c ) then
                   szt2500c_x = szt2500c_x + szt / area_2500c
                   szk2500c_x = szk2500c_x + szk_x / area_2500c
                endif       
                if ( rx .le.r5000c ) then
                   szt5000c_x = szt5000c_x + szt / area_5000c
                   szk5000c_x = szk5000c_x + szk_x / area_5000c
                endif

c...            projection along y-axis
                if ( ry .le. rvir ) then
                   sztvir_y = sztvir_y + szt / area_vir
                   szkvir_y = szkvir_y + szk_y / area_vir
                endif
                if ( ry .le. r180 ) then
                   szt180_y = szt180_y + szt / area_180
                   szk180_y = szk180_y + szk_y / area_180
                endif
                if ( ry .le. r200 ) then
                   szt200_y = szt200_y + szt / area_200
                   szk200_y = szk200_y + szk_y / area_200
                endif
                if ( ry .le. r500 ) then
                   szt500_y = szt500_y + szt / area_500
                   szk500_y = szk500_y + szk_y / area_500
                endif
                if ( ry .le. r2500 ) then
                   szt2500_y = szt2500_y + szt / area_2500
                   szk2500_y = szk2500_y + szk_y / area_2500
                endif
                if ( ry .le. r5000 ) then
                   szt5000_y = szt5000_y + szt / area_5000
                   szk5000_y = szk5000_y + szk_y / area_5000
                endif
                
                if ( ry .le. rcrit ) then
                   sztcrit_y = sztcrit_y + szt / area_crit
                   szkcrit_y = szkcrit_y + szk_y / area_crit
                endif
                if ( ry .le. r180c ) then
                   szt180c_y = szt180c_y + szt / area_180c
                   szk180c_y = szk180c_y + szk_y / area_180c
                endif
                if ( ry .le. r200c ) then
                   szt200c_y = szt200c_y + szt / area_200c
                   szk200c_y = szk200c_y + szk_y / area_200c
                endif
                if ( ry .le. r500c ) then
                   szt500c_y = szt500c_y + szt / area_500c
                   szk500c_y = szk500c_y + szk_y / area_500c
                endif
                if ( ry .le. r2500c ) then
                   szt2500c_y = szt2500c_y + szt / area_2500c
                   szk2500c_y = szk2500c_y + szk_y / area_2500c
                endif
                if ( ry .le. r5000c ) then
                   szt5000c_y = szt5000c_y + szt / area_5000c
                   szk5000c_y = szk5000c_y + szk_y / area_5000c
                endif

c...            projection along z-axis
                if ( rz .le. rvir ) then
                   sztvir_z = sztvir_z + szt / area_vir
                   szkvir_z = szkvir_z + szk_z / area_vir
                endif
                if ( rz .le. r180 ) then
                   szt180_z = szt180_z + szt / area_180
                   szk180_z = szk180_z + szk_z / area_180
                endif
                if ( rz .le. r200 ) then
                   szt200_z = szt200_z + szt / area_200
                   szk200_z = szk200_z + szk_z / area_200
                endif
                if ( rz .le. r500 ) then
                   szt500_z = szt500_z + szt / area_500
                   szk500_z = szk500_z + szk_z / area_500
                endif
                if ( rz .le. r2500 ) then
                   szt2500_z = szt2500_z + szt / area_2500
                   szk2500_z = szk2500_z + szk_z / area_2500
                endif
                if ( rz .le. r5000 ) then
                   szt5000_z = szt5000_z + szt / area_5000
                   szk5000_z = szk5000_z + szk_z / area_5000
                endif
                
                if ( rz .le. rcrit ) then
                   sztcrit_z = sztcrit_z + szt / area_crit
                   szkcrit_z = szkcrit_z + szk_z / area_crit
                endif
                if ( rz .le. r180c ) then
                   szt180c_z = szt180c_z + szt / area_180c
                   szk180c_z = szk180c_z + szk_z / area_180c
                endif
                if ( rz .le. r200c ) then
                   szt200c_z = szt200c_z + szt / area_200c
                   szk200c_z = szk200c_z + szk_z / area_200c
                endif
                if ( rz .le. r500c ) then
                   szt500c_z = szt500c_z + szt / area_500c
                   szk500c_z = szk500c_z + szk_z / area_500c
                endif
                if ( rz .le. r2500c ) then
                   szt2500c_z = szt2500c_z + szt / area_2500c
                   szk2500c_z = szk2500c_z + szk_z / area_2500c
                endif
                if ( rz .le. r5000c ) then
                   szt5000c_z = szt5000c_z + szt / area_5000c
                   szk5000c_z = szk5000c_z + szk_z / area_5000c
                endif

              endif
            enddo
          enddo
        ENDIF
      ENDDO
            
c...  Writing data to the file
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,19) xc , yc , zc  
      write(12,20) ifindcenter
      write(12,21)
      write(12,23) Deltavir, rvir*r0, aMvir, aMvir_dm, aMvir_gas, 
     &     aMvir_st, sztvir_r, sztvir_x, szkvir_x, sztvir_y, 
     &     szkvir_y, sztvir_z, szkvir_z 
      write(12,23) 180., r180*r0, aM180, aM180_dm, aM180_gas, aM180_st, 
     &     szt180_r, szt180_x, szk180_x, szt180_y, szk180_y, 
     &     szt180_z, szk180_z 
      write(12,23) 200., r200*r0, aM200, aM200_dm, aM200_gas, aM200_st,
     &     szt200_r, szt200_x, szk200_x, szt200_y, szk200_y, 
     &     szt200_z, szk200_z 
      write(12,23) 500., r500*r0, aM500, aM500_dm, aM500_gas, aM500_st,
     &     szt500_r, szt500_x, szk500_x, szt500_y, szk500_y, 
     &     szt500_z, szk500_z 
      write(12,23) 2500., r2500*r0, aM2500, aM2500_dm, aM2500_gas, 
     &     aM2500_st, szt2500_r, szt2500_x, szk2500_x, 
     &     szt2500_y, szk2500_y, szt2500_z, szk2500_z 
      write(12,23) 5000., r5000*r0, aM5000, aM5000_dm, aM5000_gas, 
     &     aM5000_st, szt5000_r, szt5000_x, szk5000_x, 
     &     szt5000_y, szk5000_y, szt5000_z, szk5000_z 
      
      write(12,22)
      write(12,23) Deltacrit, rcrit*r0, aMcrit, aMcrit_dm, 
     &     aMcrit_gas, aMcrit_st, sztcrit_r, sztcrit_x, 
     &     szkcrit_x, sztcrit_y, szkcrit_y, sztcrit_z, szkcrit_z 
      write(12,23) 180., r180c*r0, aM180c, aM180c_dm, aM180c_gas, 
     &     aM180c_st, szt180c_r, szt180c_x, szk180c_x, 
     &     szt180c_y, szk180c_y, szt180c_z, szk180c_z 
      write(12,23) 200., r200c*r0, aM200c, aM200c_dm, aM200c_gas, 
     &     aM200c_st, szt200c_r, szt200c_x, szk200c_x, 
     &     szt200c_y, szk200c_y, szt200c_z, szk200c_z 
      write(12,23) 500., r500c*r0, aM500c, aM500c_dm, aM500c_gas, 
     &     aM500c_st, szt500c_r, szt500c_x, szk500c_x, 
     &     szt500c_y, szk500c_y, szt500c_z, szk500c_z 
      write(12,23) 2500., r2500c*r0, aM2500c, aM2500c_dm, aM2500c_gas, 
     &     aM2500c_st, szt2500c_r, szt2500c_x, szk2500c_x, 
     &     szt2500c_y, szk2500c_y, szt2500c_z, szk2500c_z 
      write(12,23) 5000., r5000c*r0, aM5000c, aM5000c_dm, aM5000c_gas, 
     &     aM5000c_st, szt5000c_r, szt5000c_x, szk5000c_x, 
     &     szt5000c_y, szk5000c_y, szt5000c_z, szk5000c_z 
 19   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 20   format ('# ifindcenter = ',i5)
 21   format ('# Dvir(mean)  Rvir [1/h Mpc]  Mvir        Mdm     ',
     &     '     Mgas       Mst [1/h Msun]      ', 
     &     '     <y>(<Rvir)      <b>(<Rvir)  for r,x,y,z projections')
 22   format ('# Dcrit(crit) Rcrit [1/h Mpc] Mcrit       Mdm     ', 
     &     '     Mgas       Mst [1/h Msun]      ',
     &     '     <y>(<Rcrit)     <b>(<Rcrit) for r,x,y,z projections')
 23   format (f9.2,4x,f7.4,3x,12(g13.5,1x))
      close ( 12 )
      write(*,*) 'done writing SZ data. exiting...'
    
      return
      end
