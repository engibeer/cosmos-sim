c...  Analysis routines
c        Analyze_Xray
c        Output_Xray_3D_Profiles
c        Output_Xray_2D_Profiles
c        Output_Xray_2D_Smooth_Profiles
c        Output_Xray_2D_Smooth2_Profiles
c        XSpectra
c        Output_XrayVir
c        Compute_Xray
c
c...  Raymond-Smith spectra 
c        Read_RSdet_Table
c        Read_RS_Table
c        RS()
c        Read_RS_Spec()
c        RS_Spec()
c
c     -----------------------------------------------
      subroutine Analyze_Xray ( if1, if2, if3, idet )
c     -----------------------------------------------
c     
c     Driver for computing X-ray observables
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'

      character*256 fname, rs_dir, rs_fname, rs_spec_fname, obs
      integer lfrs, lobs
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*5 chbox
      logical lpointsout, l_force_center

c
c...  define the cooling function 
c
      rs_dir = '/dworkin/home/daisuke/Raymond/DAT/'
      lfrs = index(rs_dir , ' ') - 1
      rs_spec_fname = rs_dir(1:lfrs)//'rs_spec_0.01-10.01_keV.dat'
     
      ! Chandra 
      if ( idet .eq. 1 ) then
         rs_fname  = rs_dir(1:lfrs)//
     &        'rs_det_chandra_acis-s-bi_0.5-7.0keV.dat'
         obs = 'chandra_0.5-7.0keV'
      endif

      ! BeppoSAX
      if ( idet .eq. 2 ) then
         rs_fname  = rs_dir(1:lfrs)//'rs_det_sax_mecs_2-10keV.dat'
         obs = 'sax_2-10keV'
      endif

      lobs = index(obs    , ' ') - 1

c
c.... Compute X-ray observables
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
c...  compute 3D X-ray profiles
c

      if ( if1 .eq. 1 .or. if5 .eq. 1 ) then

        fname  = path(1:lpath)//'/pro_xray3d'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'_'//obs(1:lobs)//'.dat '
      
        write(*,*) 'Output_Xray_3D_Profiles'
        call Output_Xray_3D_Profiles ( 
     &                      xc , yc , zc , rvir, rcrit,  
     &                      rmin , rmax , nrbin , ibintype ,  
     &                      deltavir, deltacrit, ifindcenter , 
     &                      lpointsout, l_force_center,
     &                      fname, rs_fname ) 

        write(*,11) xc, yc, zc
        write(*,12) rvir*r0
        write(*,13) rcrit*r0

      endif

 11   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 12   format ('# rvir = ', (g13.7,1x),' 1/h Mpc')
 13   format ('# rcrit = ', (g13.7,1x),' 1/h Mpc')

c
c...  compute 2D X-ray profiles by directly projecting down values of each cell
c
      
      if ( if2 .eq. 1 ) then
         nrbin = 80
         rmin  = 1.0e-3/r0        ! 1/h kpc
         rmax  = 10.0/r0          ! 10/h Mpc
         rcore = 6.0e-2/r0        ! core to exclude (< 60/h kpc)

         fname  = path(1:lpath)//'/pro_xray2d'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'_'//obs(1:lobs)//'.dat '

         write(*,*) 'Output_Xray_2D_Profiles'
         call Output_Xray_2D_Profiles ( xc, yc, zc, 
     &        deltavir, deltacrit, rcore, 
     &        cnt_rate, rmin , rmax , nrbin , ibintype, 
     &        l_force_center, ifindcenter, fname, rs_fname )
      endif

c
c...  compute 2D X-ray profiles by distributing photons with gaussian deviate
c
      if ( if3 .eq. 1 ) then
         nrbin = 80
         rmin  = 1.0e-3/r0       ! 1/h kpc
         rmax  = 10.0/r0         ! 10/h Mpc
         rcore = 6.0e-2/r0       ! core to exclude (< 60/h kpc)

         fname  = path(1:lpath)//'/pro_xray2d'//'_smooth_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'_'//obs(1:lobs)//'.dat '

         write(*,*) 'Output_Xray_2D_Smooth_Profiles'
         call Output_Xray_2D_Smooth_Profiles ( xc, yc, zc, 
     &        deltavir, deltacrit, rcore,  
     &        cnt_rate, rmin , rmax , nrbin , ibintype, 
     &        l_force_center, ifindcenter, fname, rs_fname )
      endif

      return
      end
c
c     --------------------------------------------------------------
      subroutine Output_Xray_3D_Profiles ( 
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
      real norm_rs
      double precision mp, mpc_2_cm, lambda, Lx, norm_xsbp
      double precision ex_fact, lx_fact
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
      real rbin(0:nmaxbins),dgbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), dtibin(0:nmaxbins)
      real xsbpbin(0:nmaxbins),tewbin(0:nmaxbins),zewbin(0:nmaxbins)
      real tew_irbin(0:nmaxbins),tmbin(0:nmaxbins),tm_irbin(0:nmaxbins)
      real ndmp(0:nmaxbins)
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

      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      a5i = 1.0 / a5
      gamma1 = gamma - 1.0 

         
      zsim = 1.0/aexpn - 1.0
      fy = (1.0-0.5*Y_p)*(1.0-Y_p)
      mp = 1.6726d-24 
      mpc_2_cm = mpc * 100000.0
      norm_rs = 10.**(-23)

      norm_xsbp = fy * (rho0/a3*sun_mass/mp)**2 / mpc_2_cm**6 
      norm_xsbp = norm_xsbp * norm_rs * r0 * mpc_2_cm / hubble
      norm_xsbp = norm_xsbp / (4.0*pi*(1+zsim)**4)*(4.0*(pi/21600)**2)

c...  Read in the Raymond-Smith cooling function from file
      irs = 1    ! 0 = RS_det(T,Z), 1 = RS_cgs(T,Z), 2 = RS(T,Z=0.3Msun),
      if ( irs .eq. 2 ) then
        call read_rs()
      else
        Ynow = 0.2496
        zsim = 1.0/aexpn - 1.0
        if ( zsim .le. 1.0e-20 ) zsim = 0.0
        ! Normalization for XSB [erg s^-1 cm^-2 arcmin^-2]
        fy = (1.0 - 0.5 * Ynow) * (1.0 - Ynow)
        ex_fact = 3.89535e-9 * fy * a5i * hubble**3 * Om0**2 * r0
        ex_fact = ex_fact / (4.0*pi*(1+zsim)**4) * (4.0*(pi/21600)**2)

        lx_fact = (4.0*pi*(1+zsim)**4)/(4.0*(pi/21600)**2)
        lx_fact = lx_fact * ( r0 * mpc_2_cm )**2 

        if ( irs .eq. 0 ) call Read_RSdet_Table ( 0, rs_fname )
        if ( irs .eq. 1 ) call Read_RSdet_Table ( 1, rs_fname )
        write(*,*) 'Finished reading the RS spectra...'
      endif

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
        dgbin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        xsbpbin(i) = 0.
	tewbin(i) = 0.
        zewbin(i) = 0.
        tew_irbin(i) = 0.
        exbin(i) = 0.
        tmbin(i) = 0.
        tm_irbin(i) = 0.
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
c.... compute X-ray observables
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
        
c...    compute emissivity
        if ( irs .eq. 2 ) then
           call lambda_rs( Tcell, lambda )   
           ex = pd(ii)**2 * lambda 
           ex = norm_xsbp * ex 
        else 
           Tlog = log10(Tcell)
           Zlog = log10(Zcell)
           ex = pd(ii)**2 * RS(Tlog,Zlog) 
           ex = ex_fact * ex
        endif

        if ( ibin .le. nmaxbins ) then    
          nbin(ibin)   = nbin(ibin)  + 1
 	  rbin(ibin)   = rbin(ibin)  + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)

c...      add up values  & convert code units to physical units
          exbin(ibin) = exbin(ibin) + ex
          xsbpbin(ibin) = xsbpbin(ibin) + ex
          tewbin(ibin) = tewbin(ibin) + Tcell * ex
	  zewbin(ibin) = zewbin(ibin) + Zcell * ex
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
c....   compute Tew(>r)
c
      tewir = 0.0
      exir = 0.0
      tmir = 0.0
      dgir = 0.0
      do i = nbins, 0, -1 
         if ( nbin(i) .gt. 0. ) then
            if ( exbin(i) .gt. 1.0e-30 ) then
               tewir = tewir + tewbin(i)
               exir  = exir  + exbin(i)
               tew_irbin(i) = tewir / exir
               tmir  = tmir + tmbin(i)
               dgir  = dgir  + dgbin(i)
               tm_irbin(i) = tmir / dgir
            endif
         endif
      enddo

c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
c          ddmbin(i) = ddmbin(i) / nbin(i)
          if ( exbin(i) .gt. 1.0e-30 ) then
            xsbpbin(i) = xsbpbin(i) / nbin(i) 
            tewbin(i) = tewbin(i) / exbin(i)
            zewbin(i) = zewbin(i) / exbin(i)
            tmbin(i)  = tmbin(i)  / dgbin(i) / nbin(i)
          endif
	endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.0
      dti = 0.0
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
          dmi = dmi + ddmbin(i) 
          dti = dti + ddmbin(i) + dgbin(i)*volr 
          ddmbin(i) = ddmbin(i)/volr 
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i) = dti / (4.18879 * rr**3)
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
 18   format ('# <r> r_l r_m r_r delta_g delta_dm deltai_dm',
     &        ' deltai_tot  XSBP  Tew  Zew     N_dm_p' )
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
     &         dgbin(i),ddmbin(i),ddmibin(i), dtibin(i), 
     &         xsbpbin(i),tewbin(i),zewbin(i),tew_irbin(i),
     &         tmbin(i),tm_irbin(i),int(ndmp(i))
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
c     ----------------------------------------------------------
      subroutine Output_Xray_2D_Profiles ( xc, yc, zc, 
     &        deltavir, deltacrit, rcore,  
     &        cnt_rate, rmin , rmax , nrbins , ibintype, 
     &        l_force_center, ifindcenter, fprofname, rs_fname )
c     ----------------------------------------------------------
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
      character*256 fprofname, rs_fname
c
      character*256 fname
      real norm_rs
      real tew_r, zew_r, tew_x, zew_x, tew_y, zew_y, tew_z, zew_z
      double precision mp, lambda, norm_xsbp, norm_lx, norm_cnts
      double precision mpc_2_cm, lx_r, lx_x, lx_y, lx_z
      double precision ex_r, ex_x, ex_y, ex_z
      double precision lx_rc, lx_xc, lx_yc, lx_zc
      double precision ex_rc, ex_xc, ex_yc, ex_zc
      double precision ex_fact, lx_fact
      parameter ( nmaxbins = 1000 )
      real rrbin(0:nmaxbins),rxbin(0:nmaxbins),rybin(0:nmaxbins),
     +     rzbin(0:nmaxbins),xsbp_rbin(0:nmaxbins),
     +     xsbp_xbin(0:nmaxbins),xsbp_ybin(0:nmaxbins),
     +     xsbp_zbin(0:nmaxbins),tew_rbin(0:nmaxbins),
     +     tew_xbin(0:nmaxbins),tew_ybin(0:nmaxbins),
     +     tew_zbin(0:nmaxbins),zew_rbin(0:nmaxbins),
     +     zew_xbin(0:nmaxbins),zew_ybin(0:nmaxbins),
     +     zew_zbin(0:nmaxbins)
      real tewi_rbin(0:nmaxbins),tewi_xbin(0:nmaxbins),
     +     tewi_ybin(0:nmaxbins),tewi_zbin(0:nmaxbins)
      real tewir_rbin(0:nmaxbins),tewir_xbin(0:nmaxbins),
     +     tewir_ybin(0:nmaxbins),tewir_zbin(0:nmaxbins)
      real ex_rbin(0:nmaxbins),ex_xbin(0:nmaxbins),
     +     ex_ybin(0:nmaxbins),ex_zbin(0:nmaxbins)
      real nrbin(0:nmaxbins),nxbin(0:nmaxbins),nybin(0:nmaxbins),
     +     nzbin(0:nmaxbins)
      integer nfpn
c
      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbins) 
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
      Mgas = 0.
      Lx   = 0.
      Tew  = 0.
      Zew  = 0.
      ex_sum = 0.
      coldmass = 0.
      cnt_rate = 0.
      Tmax = 100000.0
c
      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      a5i = 1.0 / a5
      gamma1 = gamma - 1.0 
c
c      zsim = 1.0/aexpn - 1.0
c      fy = (1.0-0.5*Y_p)*(1.0-Y_p)
c      mp = 1.6726d-24 
      mpc_2_cm = mpc * 100000.0
c      norm_rs = 10.**(-23)

c      norm_xsbp = fy * (rho0/a3*sun_mass/mp)**2 / mpc_2_cm**6 
c      norm_xsbp = norm_xsbp * norm_rs * r0 * mpc_2_cm / hubble 
c      norm_xsbp = norm_xsbp / (4.0*pi*(1+zsim)**4)*(4.0*(pi/21600)**2)
c
c      norm_lx = norm_xsbp * (4.0*pi*(1+zsim)**4)/(4.0*(pi/21600)**2)
c      norm_lx = norm_lx * r0**2 * mpc_2_cm**2 

c      norm_cnts = fy * (rho0/a3*sun_mass/mp)**2 / mpc_2_cm**6 
c      norm_cnts = norm_cnts * norm_rs * ( r0 * mpc_2_cm / hubble )**3 
      
c...  Read in the Raymond-Smith cooling function from file
      irs = 1    ! 0 = RS_det(T,Z), 1 = RS_cgs(T,Z), 2 = RS(T,Z=0.3Msun),
      if ( irs .eq. 2 ) then
        call read_rs()
      else
        Ynow = 0.2496
        zsim = 1.0/aexpn - 1.0
        if ( zsim .le. 1.0e-20 ) zsim = 0.0
        ! Normalization for XSB [erg s^-1 cm^-2 arcmin^-2]
        fy = (1.0 - 0.5 * Ynow) * (1.0 - Ynow)
        ex_fact = 3.89535e-9 * fy * a5i * hubble**3 * Om0**2 * r0
        ex_fact = ex_fact / (4.0*pi*(1+zsim)**4) * (4.0*(pi/21600)**2)

        lx_fact = (4.0*pi*(1+zsim)**4)/(4.0*(pi/21600)**2)
        lx_fact = lx_fact * ( r0 * mpc_2_cm )**2 

        norm_cnts = ( aexpn * r0 / hubble * mpc_2_cm )**2 

        if ( irs .eq. 0 ) call Read_RSdet_Table ( 0, rs_fname )
        if ( irs .eq. 1 ) call Read_RSdet_Table ( 1, rs_fname )
        write(*,*) 'Finished reading the RS spectra...'
      endif

      call Get_MaxLevelNow ()

c.... construct profiles
      do i = 0 , nmaxbins
        rrbin(i) = 0.
        xsbp_rbin(i) = 0.
	tew_rbin(i) = 0.
        zew_rbin(i) = 0.
        ex_rbin(i) = 0.
        nrbin(i) = 0.
         
        rxbin(i) = 0.
        xsbp_xbin(i) = 0.
	tew_xbin(i) = 0.
        zew_xbin(i) = 0.
        ex_xbin(i) = 0.
        nxbin(i) = 0.

        rybin(i) = 0.
        xsbp_ybin(i) = 0.
	tew_ybin(i) = 0.
        zew_ybin(i) = 0.
        ex_ybin(i) = 0.
        nybin(i) = 0.

        rzbin(i) = 0.
        xsbp_zbin(i) = 0.
	tew_zbin(i) = 0.
        zew_zbin(i) = 0.
        ex_zbin(i) = 0.
        nzbin(i) = 0.
      enddo

c...  Loop over pixel-by-pixel
      DO Level = MinLevel , MaxLevelNow
        CellVolume = 1.0 * 2.0**(-3.0*Level)
        CellArea = 1.0 * 2.0**(-2.0*Level)
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
      
                irbin = max ( int((log10(rr) - rlmin)/drl) + 1 , 0 )
                ixbin = max ( int((log10(rx) - rlmin)/drl) + 1 , 0 )
                iybin = max ( int((log10(ry) - rlmin)/drl) + 1 , 0 )
                izbin = max ( int((log10(rz) - rlmin)/drl) + 1 , 0 )
      
c...            compute temperature and metallicity of the cell
                if ( hvar(1,ic) .gt. 0. ) then
                   Tcell = gamma1 * hvar(8,ic) / hvar(1,ic)
                   Tcell = T_0 * Tcell * a2i
                   Zcell = hvar(izIa,ic) + hvar(izII,ic)
                   Zcell = Zcell / 0.01989 / hvar(1,ic)
                else
                   write(*,*) 'bad gas density in cell', ic,hvar(1,ic)
                endif
                
c...            compute emissivity of the cell
                if ( irs .eq. 2 ) then
                  call lambda_rs( Tcell, lambda )   
                  ex = hvar(1,ic)**2 * lambda * CellVolume
                else 
                  Tlog = log10(Tcell)
                  Zlog = log10(Zcell)
                  ex = hvar(1,ic)**2 * RS(Tlog,Zlog) * CellVolume
                  ex = ex_fact * ex
                endif

c...            exclude cells in 2 layers of transition regions in the cluster core
c...            if the neighboring pixels has T<10^5K

c                do ii = 1, 6
c                 inb1 = iNb(ic,ii)
c                 tnb1 = gamma1 * hvar(8,inb1) / hvar(1,inb1)
c                 tnb1 = tnb1 * T_0 * a2i
c                 if (tnb1 .le. Tmax) then
c                    ex = 0.0d0                    
c                 endif
c                 do jj = 1, 6
c                    inb2 = iNb(inb1,jj)
c                    tnb2 = gamma1 * hvar(8,inb2) / hvar(1,inb2)
c                    tnb2 = tnb2 * T_0 * a2i
c                    if (tnb2 .le. Tmax) then
c                       ex = 0.0d0                    
c                    endif
c                 enddo
c                enddo

c...            compute the net photon count rate in the entire box
                cnt_rate = cnt_rate + ex * norm_cnts / boltz / Tcell
c                write(*,*) ex, norm_cnts, boltz, Tcell, cnt_rate


c...            within a spherer of radius r
                if ( irbin .le. nmaxbins ) then    
                   nrbin(irbin)   = nrbin(irbin) + CellVolume
                   rrbin(irbin)   = rrbin(irbin) + rr * CellVolume
                   ex_rbin(irbin) = ex_rbin(irbin) + ex
                   xsbp_rbin(irbin) = xsbp_rbin(irbin) + ex
                   tew_rbin(irbin)  = tew_rbin(irbin)  + Tcell * ex
                   zew_rbin(irbin)  = zew_rbin(irbin)  + Zcell * ex
                endif

                if ( rr .ge. rcore .and. rr .le. rvir ) then
                   lx_r  = lx_r  + ex
                   tew_r = tew_r + Tcell * ex
                   zew_r = zew_r + Zcell * ex
                   ex_r  = ex_r  + ex
                endif
                
                if ( rr .ge. rcore .and. rr .le. rcrit ) then
                   lx_rc  = lx_rc  + ex
                   tew_rc = tew_rc + Tcell * ex
                   zew_rc = zew_rc + Zcell * ex
                   ex_rc  = ex_rc  + ex
                endif

c...            projection along x-axis
                if ( ixbin .le. nmaxbins ) then    
                   nxbin(ixbin)   = nxbin(ixbin) + CellVolume
                   rxbin(ixbin)   = rxbin(ixbin) + rx * CellVolume
                   ex_xbin(ixbin) = ex_xbin(ixbin) + ex
                   xsbp_xbin(ixbin) = xsbp_xbin(ixbin) + ex
                   tew_xbin(ixbin)  = tew_xbin(ixbin)  + Tcell * ex
                   zew_xbin(ixbin)  = zew_xbin(ixbin)  + Zcell * ex
                endif

                if ( rx .ge. rcore .and. rx .le. rvir ) then
                   lx_x  = lx_x  + ex
                   tew_x = tew_x + Tcell * ex
                   zew_x = zew_x + Zcell * ex
                   ex_x  = ex_x  + ex
                endif
                
                if ( rx .ge. rcore .and. rx .le. rcrit ) then
                   lx_xc  = lx_xc  + ex
                   tew_xc = tew_xc + Tcell * ex
                   zew_xc = zew_xc + Zcell * ex
                   ex_xc  = ex_xc  + ex
                endif
      
c...            projection along y-axis
                if ( iybin .le. nmaxbins ) then    
                   nybin(iybin)   = nybin(iybin) + CellVolume
                   rybin(iybin)   = rybin(iybin) + ry * CellVolume
                   ex_ybin(iybin) = ex_ybin(iybin) + ex
                   xsbp_ybin(iybin) = xsbp_ybin(iybin) + ex
                   tew_ybin(iybin)  = tew_ybin(iybin)  + Tcell * ex
                   zew_ybin(iybin)  = zew_ybin(iybin)  + Zcell * ex
                endif
      
                if ( ry .ge. rcore .and. ry .le. rvir ) then
                   lx_y  = lx_y  + ex
                   tew_y = tew_y + Tcell * ex
                   zew_y = zew_y + Zcell * ex
                   ex_y  = ex_y  + ex
                endif

                if ( ry .ge. rcore .and. ry .le. rcrit ) then
                   lx_yc  = lx_yc  + ex
                   tew_yc = tew_yc + Tcell * ex
                   zew_yc = zew_yc + Zcell * ex
                   ex_yc  = ex_yc  + ex
                endif

c...            projection along z-axis
                if ( izbin .le. nmaxbins ) then    
                   nzbin(izbin)   = nzbin(izbin) + CellVolume
                   rzbin(izbin)   = rzbin(izbin) + rz * CellVolume
                   ex_zbin(izbin) = ex_zbin(izbin) + ex
                   xsbp_zbin(izbin) = xsbp_zbin(izbin) + ex
                   tew_zbin(izbin)  = tew_zbin(izbin)  + Tcell * ex
                   zew_zbin(izbin)  = zew_zbin(izbin)  + Zcell * ex
                endif
      
                if ( rz .ge. rcore .and. rz .le. rvir ) then
                   lx_z  = lx_z  + ex
                   tew_z = tew_z + Tcell * ex
                   zew_z = zew_z + Zcell * ex
                   ex_z  = ex_z  + ex
                endif

                if ( rz .ge. rcore .and. rz .le. rcrit ) then
                   lx_zc  = lx_zc  + ex
                   tew_zc = tew_zc + Tcell * ex
                   zew_zc = zew_zc + Zcell * ex
                   ex_zc  = ex_zc  + ex
                endif

              endif
            enddo
          enddo
        ENDIF
      ENDDO

c
c....   compute Tew(>R)
c
      tewir_r = 0.0
      tewir_x = 0.0
      tewir_y = 0.0
      tewir_z = 0.0
      exir_r = 0.0
      exir_x = 0.0
      exir_y = 0.0
      exir_z = 0.0

      do i = nbins, 0, -1 

c...     r-axis
         if ( nrbin(i) .gt. 0. ) then
            if ( ex_rbin(i) .gt. 1.0e-30 ) then
               tewir_r = tewir_r + tew_rbin(i)
               exir_r  = exir_r  + ex_rbin(i)
               tewir_rbin(i) = tewir_r / exir_r
c               write(*,*) tew_xbin(i),ex_xbin(i),tewir_xbin(i)
            endif
         endif

c...     x-axis
         if ( nxbin(i) .gt. 0. ) then
            if ( ex_xbin(i) .gt. 1.0e-30 ) then
               tewir_x = tewir_x + tew_xbin(i)
               exir_x  = exir_x  + ex_xbin(i)
               tewir_xbin(i) = tewir_x / exir_x
c               write(*,*) tew_xbin(i),ex_xbin(i),tewir_xbin(i)
            endif
         endif
      
c...     y-axis
         if ( nybin(i) .gt. 0. ) then
            if ( ex_ybin(i) .gt. 1.0e-30 ) then
               tewir_y = tewir_y + tew_ybin(i)
               exir_y  = exir_y  + ex_ybin(i)
               tewir_ybin(i) = tewir_y / exir_y
            endif
         endif
      
c...     z-axis
         if ( nzbin(i) .gt. 0. ) then
            if ( ex_zbin(i) .gt. 1.0e-30 ) then
               tewir_z = tewir_z + tew_zbin(i)
               exir_z  = exir_z  + ex_zbin(i)
               tewir_zbin(i) = tewir_z / exir_z
            endif
         endif            
      enddo

c
c....   prepare profiles
c
      tewi_r = 0.0
      tewi_x = 0.0
      tewi_y = 0.0
      tewi_z = 0.0
      exi_r = 0.0
      exi_x = 0.0
      exi_y = 0.0
      exi_z = 0.0

      do i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlmin + float(i-1)*drl)
         endif
         rr = 10.**(rlmin + float(i)*drl)
         area = pi * (rr**2 - rl**2)  

c...     r-axis
         if ( nrbin(i) .gt. 0. ) then
            rrbin(i) = rrbin(i) / nrbin(i) * r0
            if ( ex_rbin(i) .gt. 1.0e-30 ) then
               if ( rrbin(i) .ge. rcore*r0 ) then
                 tewi_r = tewi_r + tew_rbin(i)
                 exi_r  = exi_r  + ex_rbin(i)
                 tewi_rbin(i) = tewi_r / exi_r
               endif
               xsbp_rbin(i) = xsbp_rbin(i) / area 
               tew_rbin(i) = tew_rbin(i) / ex_rbin(i)
               zew_rbin(i) = zew_rbin(i) / ex_rbin(i)
            endif
         endif
      
c...     x-axis
         if ( nxbin(i) .gt. 0. ) then
            rxbin(i) = rxbin(i) / nxbin(i) * r0
            if ( ex_xbin(i) .gt. 1.0e-30 ) then
               if ( rxbin(i) .ge. rcore*r0 ) then
                 tewi_x = tewi_x + tew_xbin(i)
                 exi_x  = exi_x  + ex_xbin(i)
                 tewi_xbin(i) = tewi_x / exi_x
               endif
               xsbp_xbin(i) = xsbp_xbin(i) / area 
               tew_xbin(i) = tew_xbin(i) / ex_xbin(i)
               zew_xbin(i) = zew_xbin(i) / ex_xbin(i)
            endif
         endif
      
c...     y-axis
         if ( nybin(i) .gt. 0. ) then
            rybin(i) = rybin(i) / nybin(i) * r0
            if ( ex_ybin(i) .gt. 1.0e-30 ) then
               if ( rybin(i) .ge. rcore*r0 ) then
                 tewi_y = tewi_y + tew_ybin(i)
                 exi_y  = exi_y  + ex_ybin(i)
                 tewi_ybin(i) = tewi_y / exi_y
               endif
               xsbp_ybin(i) = xsbp_ybin(i) / area 
               tew_ybin(i) = tew_ybin(i) / ex_ybin(i)
               zew_ybin(i) = zew_ybin(i) / ex_ybin(i)
            endif
         endif
      
c...     z-axis
         if ( nzbin(i) .gt. 0. ) then
            rzbin(i) = rzbin(i) / nzbin(i) * r0
            if ( ex_zbin(i) .gt. 1.0e-30 ) then
               if ( rzbin(i) .ge. rcore*r0 ) then
                 tewi_z = tewi_z + tew_zbin(i)
                 exi_z  = exi_z  + ex_zbin(i)
                 tewi_zbin(i) = tewi_z / exi_z
               endif
               xsbp_zbin(i) = xsbp_zbin(i) / area 
               tew_zbin(i) = tew_zbin(i) / ex_zbin(i)
               zew_zbin(i) = zew_zbin(i) / ex_zbin(i)
            endif
         endif            
      enddo


c...  compute emission-weighted temperature and metallicity      
      lx_r  = lx_r * lx_fact
      tew_r = tew_r / ex_r / (1.160485e+7)
      zew_r = zew_r / ex_r
      lx_x  = lx_x * lx_fact
      tew_x = tew_x / ex_x / (1.160485e+7)
      zew_x = zew_x / ex_x
      lx_y  = lx_y * lx_fact
      tew_y = tew_y / ex_y / (1.160485e+7)
      zew_y = zew_y / ex_y
      lx_z  = lx_z * lx_fact
      tew_z = tew_z / ex_z / (1.160485e+7)
      zew_z = zew_z / ex_z

      lx_rc  = lx_rc * lx_fact
      tew_rc = tew_rc / ex_rc / (1.160485e+7)
      zew_rc = zew_rc / ex_rc
      lx_xc  = lx_xc * lx_fact
      tew_xc = tew_xc / ex_xc / (1.160485e+7)
      zew_xc = zew_xc / ex_xc
      lx_yc  = lx_yc * lx_fact
      tew_yc = tew_yc / ex_yc / (1.160485e+7)
      zew_yc = zew_yc / ex_yc
      lx_zc  = lx_zc * lx_fact
      tew_zc = tew_zc / ex_zc / (1.160485e+7)
      zew_zc = zew_zc / ex_zc


c...  Writing data to the file
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
c
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbins , ibintype 
      write(12,19) rcore*r0*1000.
      write(12,12)
      write(12,14) deltavir, rvir*r0, lx_x, tew_x, zew_x  
      write(12,14) deltavir, rvir*r0, lx_y, tew_y, zew_y  
      write(12,14) deltavir, rvir*r0, lx_z, tew_z, zew_z
      write(12,13)
      write(12,14) deltacrit, rcrit*r0, lx_xc, tew_xc, zew_xc  
      write(12,14) deltacrit, rcrit*r0, lx_yc, tew_yc, zew_yc  
      write(12,14) deltacrit, rcrit*r0, lx_zc, tew_zc, zew_zc
      write(12,15) 
      write(12,16)
      write(12,17)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbins =',i5,' ibintype =',i1)
 19   format ('# rcore = ',f7.2,' [1/h kpc] is excluded') 
 12   format ('# Dvir(mean)  Rvir [1/h Mpc]  Lx [h^-2 erg/s] ',
     &     ' Tew [keV]  Zew [Zsun] ')
 13   format ('# Dcrit(crit) Rcrit [1/h Mpc] Lx [h^-2 erg/s] ',
     &     ' Tew [keV]  Zew [Zsun] ')
 14   format (f9.2,4x,f7.4,3x,g13.5,3x,g13.5,3x,g13.5)
 15   format ('# Columns:')
 16   format ('# r_l r_m r_r  rr xsbp_r tew_r zew_r nr  ',
     &     ' rx xsbp_x tew_x zew_x nx  ry xsbp_y ',
     &     ' tew_y zew_y ny  rz xsbp_z tew_z zew_z nz')
 17   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '     erg/s/cm^2/Sr   Kelvin  Zsun ')
      do i = 0 , nbins
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
     &         rl*r0, rmid, rr*r0,
     &         rrbin(i),xsbp_rbin(i),tew_rbin(i),zew_rbin(i),nrbin(i),
     &         rxbin(i),xsbp_xbin(i),tew_xbin(i),zew_xbin(i),nxbin(i),
     &         rybin(i),xsbp_ybin(i),tew_ybin(i),zew_ybin(i),nybin(i),
     &         rzbin(i),xsbp_zbin(i),tew_zbin(i),zew_zbin(i),nzbin(i),
     &         tewi_rbin(i), tewi_xbin(i),tewi_ybin(i),tewi_zbin(i),
     &         tewir_rbin(i), tewir_xbin(i),tewir_ybin(i),tewir_zbin(i)
	endif
      enddo
 18   format(31(2x,g14.7))
c 18   format(7(2x,g14.7),1x,i7,4(2x,g14.7),1x,i7,4(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      
      return
      end

c
c     ----------------------------------------------------------
      subroutine Output_Xray_2D_Smooth_Profiles ( xc, yc, zc, 
     &        deltavir, deltacrit, rcore,  
     &        cnt_rate, rmin , rmax , nrbin , ibintype, 
     &        l_force_center, ifindcenter, fprofname, rs_fname )
c     ----------------------------------------------------------
c     compute smooth XSBP, Tew, Zew profiles for an object centered 
c     at (xc,yc,zc) and radius rvir
c
c     Emission weighting scheme :
c       We simulate X-ray photons emitted from each cell and 
c       use them to compute smoothed X-ray radial profiles.
c
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
      character*256 fprofname, rs_fname
c
      character*256 fname
      real norm_rs, tew_x, zew_x, tew_y, zew_y, tew_z, zew_z
      double precision mp, texp, lambda, norm_xsbp, norm_lx, norm_cnts
      double precision mpc_2_cm, lx_x, lx_y, lx_z, ex_x, ex_y, ex_z
      double precision lx_xc, lx_yc, lx_zc, ex_xc, ex_yc, ex_zc
      double precision ex_fact, lx_fact
      parameter ( nmaxbins = 1000 )
      parameter ( nphotons = 50000000 )
      real rxbin(0:nmaxbins),rybin(0:nmaxbins),rzbin(0:nmaxbins),
     +  xsbp_xbin(0:nmaxbins),xsbp_ybin(0:nmaxbins),
     +  xsbp_zbin(0:nmaxbins),tew_xbin(0:nmaxbins),
     +  tew_ybin(0:nmaxbins),tew_zbin(0:nmaxbins),
     +  zew_xbin(0:nmaxbins),zew_ybin(0:nmaxbins),
     +  zew_zbin(0:nmaxbins)
      real tewi_xbin(0:nmaxbins),tewi_ybin(0:nmaxbins),
     +  tewi_zbin(0:nmaxbins)
      real ex_xbin(0:nmaxbins),ex_ybin(0:nmaxbins),ex_zbin(0:nmaxbins)
      real nxbin(0:nmaxbins),nybin(0:nmaxbins),nzbin(0:nmaxbins)
      integer nfpn, m
      logical continue_while
c
      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 
      nbins = int((rlmax-rlmin)/drl) + 1

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
 21   format ('# Dvir(mean)  Rvir[1/hMpc]     Mvir       Mdm  ',
     &     '       Mgas        Mst  [1/h Msun]' )
 22   format ('# Dcrit(crit) Rvir[1/h Mpc]    Mvir       Mdm  ',
     &     '       Mgas        Mst  [1/h Msun]' )
 23   format (f9.2,4x,f7.3,3x,4(g12.4,1x))

c
      Mgas = 0.
      Lx   = 0.
      Tew  = 0.
      Zew  = 0.
      ex_sum = 0.
      coldmass = 0.
      Tmax = 100000.0
c
c     seeds for gaussian deviate with zero mean and unit variance
c
      m = 28548535  
c
      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      a5i = 1.0 / a5
      gamma1 = gamma - 1.0 
c
c      zsim = 1.0/aexpn - 1.0
c      fy = (1.0-0.5*Y_p)*(1.0-Y_p)
c      mp = 1.6726d-24 
      mpc_2_cm = mpc * 100000.0
c      norm_rs = 10.**(-23)

c      norm_xsbp = fy * (rho0/a3*sun_mass/mp)**2 / mpc_2_cm**6 
c      norm_xsbp = norm_xsbp * norm_rs * r0 * mpc_2_cm / hubble 
c      norm_xsbp = norm_xsbp / (4.0*pi*(1+zsim)**4)*(4.0*(pi/21600)**2)

c      norm_lx = norm_xsbp * (4.0*pi*(1+zsim)**4)/(4.0*(pi/21600)**2)
c      norm_lx = norm_lx * r0**2 * mpc_2_cm**2 

c
c...  compute the exposure time [sec] needed to generate nphotons 
c
      texp = nphotons / cnt_rate
      norm_cnts = fy * (rho0/a3*sun_mass/mp)**2 / mpc_2_cm**6 
      norm_cnts = norm_cnts * norm_rs * ( r0 * mpc_2_cm / hubble )**3 
 
c...  Read in the Raymond-Smith cooling function from file
      irs = 1    ! 0 = RS_det(T,Z), 1 = RS_cgs(T,Z), 2 = RS(T,Z=0.3Msun),
      if ( irs .eq. 2 ) then
        call read_rs()
      else 
        Ynow = 0.2496
        zsim = 1.0/aexpn - 1.0
        if ( zsim .le. 1.0e-20 ) zsim = 0.0
        ! Normalization for XSB [erg s^-1 cm^-2 arcmin^-2]
        fy = (1.0 - 0.5 * Ynow) * (1.0 - Ynow)
        ex_fact = 3.89535e-9 * fy * a5i * hubble**3 * Om0**2 * r0
        ex_fact = ex_fact / (4.0*pi*(1+zsim)**4) * (4.0*(pi/21600)**2)

        lx_fact = (4.0*pi*(1+zsim)**4)/(4.0*(pi/21600)**2)
        lx_fact = lx_fact * ( r0 * mpc_2_cm )**2 

        norm_cnts = ( aexpn * r0 / hubble * mpc_2_cm )**2 

        if ( irs .eq. 0 ) call Read_RSdet_Table ( 0, rs_fname )
        if ( irs .eq. 1 ) call Read_RSdet_Table ( 1, rs_fname )
        write(*,*) 'Finished reading the RS spectra...'
      endif

      call Get_MaxLevelNow ()

c.... construct profiles
      do i = 0 , nmaxbins
        rxbin(i) = 0.
        xsbp_xbin(i) = 0.
	tew_xbin(i)  = 0.
        zew_xbin(i)  = 0.
        ex_xbin(i)   = 0.
        tewi_xbin(i) = 0.
        nxbin(i) = 0.

        rybin(i) = 0.
        xsbp_ybin(i) = 0.
	tew_ybin(i)  = 0.
        zew_ybin(i)  = 0.
        ex_ybin(i)   = 0.
        tewi_ybin(i) = 0.
        nybin(i) = 0.

        rzbin(i) = 0.
        xsbp_zbin(i) = 0.
	tew_zbin(i)  = 0.
        zew_zbin(i)  = 0.
        ex_zbin(i)   = 0.
        tewi_zbin(i) = 0.
        nzbin(i) = 0.
      enddo

c...  Loop over pixel-by-pixel
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

c...            compute temperature and metallicity of the cell
                if ( hvar(1,ic) .gt. 0. ) then
                   Tcell = gamma1 * hvar(8,ic) / hvar(1,ic)
                   Tcell = T_0 * Tcell * a2i
                   Zcell = hvar(izIa,ic) + hvar(izII,ic)
                   Zcell = Zcell / 0.01989 / hvar(1,ic)
                else
                   write(*,*) 'bad gas density in cell', ic,hvar(1,ic)
                endif
                
c...            compute emissivity of the cell
                if ( irs .eq. 2 ) then
                  call lambda_rs( Tcell, lambda )   
                  ex = hvar(1,ic)**2 * lambda * CellVolume
                else 
                  Tlog = log10(Tcell)
                  Zlog = log10(Zcell)
                  ex = hvar(1,ic)**2 * RS(Tlog,Zlog) * CellVolume
                  ex = ex_fact * ex
                endif

c...            exclude cells in 2 layers of transition regions in the cluster core
c...            if the neighboring pixels has T<10^5K

c                do ii = 1, 6
c                 inb1 = iNb(ic,ii)
c                 tnb1 = gamma1 * hvar(8,inb1) / hvar(1,inb1)
c                 tnb1 = tnb1 * T_0 * a2i
c                 if (tnb1 .le. Tmax) then
c                    ex = 0.0d0                    
c                 endif
c                 do jj = 1, 6
c                    inb2 = iNb(inb1,jj)
c                    tnb2 = gamma1 * hvar(8,inb2) / hvar(1,inb2)
c                    tnb2 = tnb2 * T_0 * a2i
c                    if (tnb2 .le. Tmax) then
c                       ex = 0.0d0                    
c                    endif
c                 enddo
c                enddo

c...            computes the number of photons emitted from this cell
                cnts = texp * ex * norm_cnts / boltz / Tcell
                ncnts = nearest_int( cnts )

c                write(*,*) texp, ex, norm_cnts, boltz, Tcell, cnts,ncnts

                ip = 1
                continue_while = .true.
                if ( ncnts .ge. 1 ) then
                   ex = ex / ncnts  ! [erg/s/photon]
                   fvol = CellVolume / ncnts
                endif

c...            loop over each photons
                do while ( continue_while )

c                  write(*,*) ncnts, cnts, fcnts, ip, continue_while
 
c...              perturb the positions with Gaussian random numbers
                  Posx = Posx + 0.05 * CellSize(Level) * gasdev( m ) 
                  Posy = Posy + 0.05 * CellSize(Level) * gasdev( m )
                  Posz = Posz + 0.05 * CellSize(Level) * gasdev( m )

c...              take care of periodic boundary condition
                  call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz )

                  rx = sqrt( (Posy-yc)**2 + (Posz-zc)**2 )
                  ry = sqrt( (Posx-xc)**2 + (Posz-zc)**2 )
                  rz = sqrt( (Posx-xc)**2 + (Posy-yc)**2 )
        
                  ixbin = max ( int((log10(rx) - rlmin)/drl) + 1 , 0 )
                  iybin = max ( int((log10(ry) - rlmin)/drl) + 1 , 0 )
                  izbin = max ( int((log10(rz) - rlmin)/drl) + 1 , 0 )

c...              projection along x-axis
                  if ( ixbin .le. nmaxbins ) then    
                     nxbin(ixbin) = nxbin(ixbin) + fvol
                     rxbin(ixbin) = rxbin(ixbin) + rx * fvol
                     ex_xbin(ixbin) = ex_xbin(ixbin) + ex
                     xsbp_xbin(ixbin) = xsbp_xbin(ixbin)+ ex 
                     tew_xbin(ixbin) = tew_xbin(ixbin) + Tcell * ex
                     zew_xbin(ixbin) = zew_xbin(ixbin) + Zcell * ex
                  endif
  
                  if ( rx .ge. rcore .and. rx .le. rvir ) then
                     lx_x  = lx_x  + ex 
                     tew_x = tew_x + Tcell * ex
                     zew_x = zew_x + Zcell * ex
                     ex_x  = ex_x  + ex
                  endif

                  if ( rx .ge. rcore .and. rx .le. rcrit ) then
                     lx_xc  = lx_xc  + ex
                     tew_xc = tew_xc + Tcell * ex
                     zew_xc = zew_xc + Zcell * ex
                     ex_xc  = ex_xc  + ex
                  endif

c...              projection along y-axis
                  if ( iybin .le. nmaxbins ) then    
                     nybin(iybin) = nybin(iybin) + fvol
                     rybin(iybin) = rybin(iybin) + ry * fvol
                     ex_ybin(iybin) = ex_ybin(iybin) + ex
                     xsbp_ybin(iybin) = xsbp_ybin(iybin) + ex 
                     tew_ybin(iybin) = tew_ybin(iybin) + Tcell * ex
                     zew_ybin(iybin) = zew_ybin(iybin) + Zcell * ex
                  endif
        
                  if ( ry .ge. rcore .and. ry .le. rvir ) then
                     lx_y  = lx_y  + ex
                     tew_y = tew_y + Tcell * ex
                     zew_y = zew_y + Zcell * ex
                     ex_y  = ex_y  + ex
                  endif

                  if ( ry .ge. rcore .and. ry .le. rcrit ) then
                     lx_yc  = lx_yc  + ex
                     tew_yc = tew_yc + Tcell * ex
                     zew_yc = zew_yc + Zcell * ex
                     ex_yc  = ex_yc  + ex
                  endif

  
c...              projection along z-axis
                  if ( izbin .le. nmaxbins ) then    
                     nzbin(izbin) = nzbin(izbin) + fvol
                     rzbin(izbin) = rzbin(izbin) + rz * fvol
                     ex_zbin(izbin) = ex_zbin(izbin) + ex
                     xsbp_zbin(izbin) = xsbp_zbin(izbin) + ex 
                     tew_zbin(izbin) = tew_zbin(izbin) + Tcell * ex
                     zew_zbin(izbin) = zew_zbin(izbin) + Zcell * ex
                  endif
        
                  if ( rz .ge. rcore .and. rz .le. rvir ) then
                     lx_z  = lx_z  + ex 
                     tew_z = tew_z + Tcell * ex
                     zew_z = zew_z + Zcell * ex
                     ex_z  = ex_z  + ex
                  endif

                  if ( rz .ge. rcore .and. rz .le. rcrit ) then
                     lx_zc  = lx_zc  + ex
                     tew_zc = tew_zc + Tcell * ex
                     zew_zc = zew_zc + Zcell * ex
                     ex_zc  = ex_zc  + ex
                  endif

c...              continue while-loop if ip < ncnts
                  if ( ip .lt. ncnts ) then
                     ip = ip + 1
                  else
                     continue_while = .false.
                  endif

                enddo  ! end do while
              endif 
            enddo
          enddo
        ENDIF
      ENDDO

c
c....   prepare profiles
c
      tewi_x = 0.0
      tewi_y = 0.0
      tewi_z = 0.0
      exi_x = 0.0
      exi_y = 0.0
      exi_z = 0.0

      do i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlmin + float(i-1)*drl)
         endif
         rr = 10.**(rlmin + float(i)*drl)
         area = pi * (rr**2 - rl**2)  

c...     x-axis
         if ( nxbin(i) .gt. 0. ) then
            rxbin(i) = rxbin(i) / nxbin(i) * r0
            if ( ex_xbin(i) .gt. 1.0e-30 ) then
               if ( rxbin(i) .ge. rcore*r0 ) then
                  tewi_x = tewi_x + tew_xbin(i)
                  exi_x  = exi_x  + ex_xbin(i)
                  tewi_xbin(i) = tewi_x / exi_x
               endif
               xsbp_xbin(i) = xsbp_xbin(i) / area * norm_xsbp 
               tew_xbin(i)  = tew_xbin(i)  / ex_xbin(i) 
               zew_xbin(i)  = zew_xbin(i)  / ex_xbin(i)
            endif
         endif
      
c...     y-axis
         if ( nybin(i) .gt. 0. ) then
            rybin(i) = rybin(i) / nybin(i) * r0
            if ( ex_ybin(i) .gt. 1.0e-30 ) then
               if ( rybin(i) .ge. rcore*r0 ) then
                  tewi_y = tewi_y + tew_ybin(i)
                  exi_y  = exi_y  + ex_ybin(i)
                  tewi_ybin(i) = tewi_y / exi_y
               endif
               xsbp_ybin(i) = xsbp_ybin(i) / area * norm_xsbp 
               tew_ybin(i)  = tew_ybin(i)  / ex_ybin(i) 
               zew_ybin(i)  = zew_ybin(i)  / ex_ybin(i) 
            endif
         endif
      
c...     z-axis
         if ( nzbin(i) .gt. 0. ) then
            rzbin(i) = rzbin(i) / nzbin(i) * r0
            if ( ex_zbin(i) .gt. 1.0e-30 ) then
               if ( rzbin(i) .ge. rcore*r0 ) then
                  tewi_z = tewi_z + tew_zbin(i)
                  exi_z  = exi_z  + ex_zbin(i)
                  tewi_zbin(i) = tewi_z / exi_z
               endif
               xsbp_zbin(i) = xsbp_zbin(i) / area * norm_xsbp
               tew_zbin(i)  = tew_zbin(i)  / ex_zbin(i) 
               zew_zbin(i)  = zew_zbin(i)  / ex_zbin(i) 
            endif
         endif            
      enddo


c...  compute emission-weighted temperature and metallicity
      lx_x  = lx_x * lx_fact
      tew_x = tew_x / ex_x / (1.160485e+7)
      zew_x = zew_x / ex_x
      lx_y  = lx_y * lx_fact
      tew_y = tew_y / ex_y / (1.160485e+7)
      zew_y = zew_y / ex_y
      lx_z  = lx_z * lx_fact
      tew_z = tew_z / ex_z / (1.160485e+7)
      zew_z = zew_z / ex_z

      lx_xc  = lx_xc * lx_fact
      tew_xc = tew_xc / ex_xc / (1.160485e+7)
      zew_xc = zew_xc / ex_xc
      lx_yc  = lx_yc * lx_fact
      tew_yc = tew_yc / ex_yc / (1.160485e+7)
      zew_yc = zew_yc / ex_yc
      lx_zc  = lx_zc * lx_fact
      tew_zc = tew_zc / ex_zc / (1.160485e+7)
      zew_zc = zew_zc / ex_zc


c...  Writing data to the file
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
c
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype
      write(12,19) rcore*r0*1000.
      write(12,12)
      write(12,14) deltavir, rvir*r0, lx_x, tew_x, zew_x  
      write(12,14) deltavir, rvir*r0, lx_y, tew_y, zew_y  
      write(12,14) deltavir, rvir*r0, lx_z, tew_z, zew_z
      write(12,13)
      write(12,14) deltacrit, rcrit*r0, lx_xc, tew_xc, zew_xc  
      write(12,14) deltacrit, rcrit*r0, lx_yc, tew_yc, zew_yc  
      write(12,14) deltacrit, rcrit*r0, lx_zc, tew_zc, zew_zc
      write(12,15) 
      write(12,16)
      write(12,17)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 19   format ('# rcore = ',f7.2,' [1/h kpc] is excluded') 
 12   format ('# Dvir(mean)  Rvir [1/h Mpc]  Lx [h^-2 erg/s] ',
     &     ' Tew [keV]  Zew [Zsun] ')
 13   format ('# Dcrit(crit) Rcrit [1/h Mpc] Lx [h^-2 erg/s] ',
     &     ' Tew [keV]  Zew [Zsun] ')
 14   format (f9.2,4x,f7.4,3x,g13.5,3x,g13.5,3x,g13.5)
 15   format ('# Columns:')
 16   format ('# r_l r_m r_r  r_x xsbp_x tew_x zew_x nx  r_y xsbp_y ',
     &     'tew_y zew_y ny  r_z xsbp_z tew_z zew_z nz')
 17   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '     erg/s/cm^2/Sr   Kelvin  Zsun ')
      do i = 0 , nbins
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
     &         rl*r0, rmid, rr*r0, 
     &         rxbin(i),xsbp_xbin(i),tew_xbin(i),zew_xbin(i),nxbin(i),
     &         rybin(i),xsbp_ybin(i),tew_ybin(i),zew_ybin(i),nybin(i),
     &         rzbin(i),xsbp_zbin(i),tew_zbin(i),zew_zbin(i),nzbin(i),
     &         tewi_xbin(i),tewi_ybin(i),tewi_zbin(i)
	endif
      enddo
 18   format(21(2x,g14.7))
c 18   format(7(2x,g14.7),1x,i7,4(2x,g14.7),1x,i7,4(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      
      return
      end

c
c     -----------------------------------------------------------
      subroutine Output_Xray_2D_Smooth_Profiles_old ( xc, yc, zc, 
     &        deltavir, deltacrit, rvir, rcrit,  
     &        cnt_rate, rmin , rmax , nrbin , ibintype, 
     &        l_force_center, ifindcenter, fprofname )
c     -----------------------------------------------------------
c     compute Mgas, Lx, Tew, Zew for an object centered 
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
      real   xc, yc, zc, rvir, rmin, rmax
      integer ibintype
      character*256 fprofname
c
      character*256 fname
      real   norm_rs, tew_x, zew_x, tew_y, zew_y, tew_z, zew_z
      double precision mp, mpc_2_cm, lambda, norm_xsbp, norm_lx
      double precision lx_x, lx_y, lx_z, ex_x, ex_y, ex_z
      double precision xx, yy, zz
      integer m
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      common / POINTS / xpn, ypn, zpn
      parameter ( nmaxbins = 1000 )
      real rxbin(0:nmaxbins),rybin(0:nmaxbins),rzbin(0:nmaxbins),
     +  xsbp_xbin(0:nmaxbins),xsbp_ybin(0:nmaxbins),
     +  xsbp_zbin(0:nmaxbins),tew_xbin(0:nmaxbins),
     +  tew_ybin(0:nmaxbins),tew_zbin(0:nmaxbins),
     +  zew_xbin(0:nmaxbins),zew_ybin(0:nmaxbins),
     +  zew_zbin(0:nmaxbins)
      real ex_xbin(0:nmaxbins),ex_ybin(0:nmaxbins),ex_zbin(0:nmaxbins)
      real nxbin(0:nmaxbins),nybin(0:nmaxbins),nzbin(0:nmaxbins)
      integer nfpn
c
      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 
      nbins = int((rlmax-rlmin)/drl) + 1

      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

c
      Mgas = 0.
      Lx   = 0.
      Tew  = 0.
      Zew  = 0.
      ex_sum = 0.
      coldmass = 0.
      Tmax = 100000.0
c
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
c
      zsim = 1.0/aexpn - 1.0
      fy = (1.0-0.5*Y_p)*(1.0-Y_p)
      mp = 1.6726d-24 
      mpc_2_cm = mpc * 100000.0
      norm_rs = 10.**(-23)

      norm_xsbp = fy * (rho0/a3*sun_mass/mp)**2 / mpc_2_cm**6 
      norm_xsbp = norm_xsbp * norm_rs * r0 * mpc_2_cm / hubble
      norm_xsbp = norm_xsbp / (4.0*pi*(1+zsim)**4)*(4.0*(pi/21600)**2)

      norm_lx = (4.0*pi*(1+zsim)**4)/(4.0*(pi/21600)**2)
      norm_lx = norm_lx * r0**2 * mpc_2_cm**2 

      
c...  Read in the Raymond-Smith cooling function from file
      call read_rs()
      
c.... construct profiles
      do i = 0 , nmaxbins
        rxbin(i) = 0.
        xsbp_xbin(i) = 0.
	tew_xbin(i) = 0.
        zew_ybin(i) = 0.
        ex_xbin(i) = 0.
        nxbin(i) = 0.

        rybin(i) = 0.
        xsbp_ybin(i) = 0.
	tew_ybin(i) = 0.
        zew_ybin(i) = 0.
        ex_ybin(i) = 0.
        nybin(i) = 0.

        rzbin(i) = 0.
        xsbp_zbin(i) = 0.
	tew_zbin(i) = 0.
        zew_zbin(i) = 0.
        ex_zbin(i) = 0.
        nzbin(i) = 0.
      enddo

      ! seed for a random number
      m = 28548535

c
c...    projection along x-axis
c      
      do i = 1 , npoints

c...    select random samples in the 2D plane
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m )
        zlos  = ng * RANDd ( m )
        xpn(i) = zlos
        ypn(i) = rr * cos(theta) + yc
        zpn(i) = rr * sin(theta) + zc
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
        
        icx = iFindCell ( MaxLevel , xx , yy , zz )
        rr = sqrt( (yy-yc)**2 + (zz-zc)**2 )
        ibin = max ( int((log10(rr) - rlmin)/drl) + 1 , 0 )
        
        if ( ibin .le. nmaxbins ) then    
           ! compute temperature and metallicity of the cell
           if ( hvar(1,icx) .gt. 0. ) then
              Tcell = gamma1 * hvar(8,icx) / hvar(1,icx)
              Tcell = T_0 * Tcell * a2i
              Zcell = hvar(izIa,icx) + hvar(izII,icx)
              Zcell = Zcell / 0.01989 / hvar(1,icx)
           else
              write(*,*) 'bad gas density in cell', icx,hvar(1,icx)
           endif
        
           ! compute emissivity
           call lambda_rs( Tcell, lambda )   
           ex = hvar(1,icx)**2 * lambda * CellVolume
c           write(*,*) Tcell
           
           ! assign values to arrays                                     
           nxbin(ibin)   = nxbin(ibin) + CellVolume
           rxbin(ibin)   = rxbin(ibin) + rr * CellVolume
           ex_xbin(ibin) = ex_xbin(ibin) + ex
           xsbp_xbin(ibin) = xsbp_xbin(ibin) + norm_xsbp * ex
           tew_xbin(ibin)  = tew_xbin(ibin)  + Tcell * ex
           zew_xbin(ibin)  = zew_xbin(ibin)  + Zcell * ex
        endif
       
        if ( rx .lt. rvir ) then
           lx_x  = lx_x  + norm_xsbp * ex
           tew_x = tew_x + Tcell * ex
           zew_x = zew_x + Zcell * ex
           ex_x  = ex_x  + ex
        endif
        
      enddo

c
c...    projection along y-axis
c
      do i = 1 , npoints

c...    select random samples in the 2D plane
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m )
        zlos  = ng * RANDd ( m )
        xpn(i) = rr * cos(theta) + xc
        ypn(i) = zlos
        zpn(i) = rr * sin(theta) + zc
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
        
        icy = iFindCell ( MaxLevel , xx, yy , zz )
        rr = sqrt( (xx-xc)**2 + (zz-zc)**2 )
        ibin = max ( int((log10(rr) - rlmin)/drl) + 1 , 0 )
        
        if ( ibin .le. nmaxbins ) then    
           ! compute temperature and metallicity of the cell
           if ( hvar(1,icy) .gt. 0. ) then
              Tcell = gamma1 * hvar(8,icy) / hvar(1,icy)
              Tcell = T_0 * Tcell * a2i
              Zcell = hvar(izIa,icy) + hvar(izII,icy)
              Zcell = Zcell / 0.01989 / hvar(1,icy)
           else
              write(*,*) 'bad gas density in cell', icy,hvar(1,icy)
           endif
        
           ! compute emissivity
           call lambda_rs( Tcell, lambda )   
           ex = hvar(1,icy)**2 * lambda * CellVolume

           ! assign values to arrays                                     
           nybin(ibin)   = nybin(ibin) + CellVolume
           rybin(ibin)   = rybin(ibin) + rr * CellVolume
           ex_ybin(ibin) = ex_ybin(ibin) + ex
           xsbp_ybin(ibin) = xsbp_ybin(ibin) + norm_xsbp * ex
           tew_ybin(ibin)  = tew_ybin(ibin)  + Tcell * ex
           zew_ybin(ibin)  = zew_ybin(ibin)  + Zcell * ex
        endif
        
        if ( ry .lt. rvir ) then
           lx_y  = lx_y  + norm_xsbp * ex
           tew_y = tew_y + Tcell * ex
           zew_y = zew_y + Zcell * ex
           ex_y  = ex_y  + ex
        endif
        
      enddo


c     
c...    projection along z-axis
c
      do i = 1 , npoints

c...    select random samples in the 2D plane
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m )
        zlos  = ng * RANDd ( m )
        xpn(i) = rr * sin(theta) + xc
        ypn(i) = rr * cos(theta) + yc
        zpn(i) = zlos
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
        
        icz = iFindCell ( MaxLevel , xx , yy, zz)
        rr = sqrt( (xx-xc)**2 + (yy-yc)**2 )
        ibin = max ( int((log10(rr) - rlmin)/drl) + 1 , 0 )
        
        if ( ibin .le. nmaxbins ) then    
           ! compute temperature and metallicity of the cell
           if ( hvar(1,icz) .gt. 0. ) then
              Tcell = gamma1 * hvar(8,icz) / hvar(1,icz)
              Tcell = T_0 * Tcell * a2i
              Zcell = hvar(izIa,icz) + hvar(izII,icz)
              Zcell = Zcell / 0.01989 / hvar(1,icz)
           else
              write(*,*) 'bad gas density in cell', icz,hvar(1,icz)
           endif
        
           ! compute emissivity
           call lambda_rs( Tcell, lambda )   
           ex = hvar(1,icz)**2 * lambda * CellVolume

           ! assign values to arrays                  
           nzbin(ibin)   = nzbin(ibin) + CellVolume
           rzbin(ibin)   = rzbin(ibin) + rr * CellVolume
           ex_zbin(ibin) = ex_zbin(ibin) + ex
           xsbp_zbin(ibin) = xsbp_zbin(ibin) + norm_xsbp * ex
           tew_zbin(ibin)  = tew_zbin(ibin)  + Tcell * ex
           zew_zbin(ibin)  = zew_zbin(ibin)  + Zcell * ex
        endif
        
        if ( rz .lt. rvir ) then
           lx_z  = lx_z  + norm_xsbp * ex
           tew_z = tew_z + Tcell * ex
           zew_z = zew_z + Zcell * ex
           ex_z  = ex_z  + ex
        endif

      enddo

c
c....   prepare profiles
c
      do i = 0 , nbins
c...     x-axis
         if ( nxbin(i) .gt. 0. ) then
            rxbin(i) = rxbin(i) / nxbin(i) * r0
            if ( ex_xbin(i) .gt. 1.0e-30 ) then
               xsbp_xbin(i) = xsbp_xbin(i) / nxbin(i)
               tew_xbin(i) = tew_xbin(i) / ex_xbin(i)
               zew_xbin(i) = zew_xbin(i) / ex_xbin(i)
            endif
         endif
      
c...     y-axis
         if ( nybin(i) .gt. 0. ) then
            rybin(i) = rybin(i) / nybin(i) * r0
            if ( ex_ybin(i) .gt. 1.0e-30 ) then
               xsbp_ybin(i) = xsbp_ybin(i) / nybin(i)
               tew_ybin(i) = tew_ybin(i) / ex_ybin(i)
               zew_ybin(i) = zew_ybin(i) / ex_ybin(i)
            endif
         endif
      
c...     z-axis
         if ( nzbin(i) .gt. 0. ) then
            rzbin(i) = rzbin(i) / nzbin(i) * r0
            if ( ex_zbin(i) .gt. 1.0e-30 ) then
               xsbp_zbin(i) = xsbp_zbin(i) / nzbin(i)
               tew_zbin(i) = tew_zbin(i) / ex_zbin(i)
               zew_zbin(i) = zew_zbin(i) / ex_zbin(i)
            endif
         endif            
      enddo


c...  compute emission-weighted temperature and metallicity
      lx_x  = lx_x * norm_lx
      tew_x = tew_x / ex_x / (1.160485e+7)
      zew_x = zew_x / ex_x
      lx_y  = lx_y * norm_lx
      tew_y = tew_y / ex_y / (1.160485e+7)
      zew_y = zew_y / ex_y
      lx_z  = lx_z * norm_lx
      tew_z = tew_z / ex_z / (1.160485e+7)
      zew_z = zew_z / ex_z


c...  Writing data to the file
      open ( 12 , file = fprofname(1:nfpn) )
c
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype 
      write(12,12) lx_x, tew_x, zew_x ,rvir*r0 
      write(12,13) lx_y, tew_y, zew_y ,rvir*r0 
      write(12,14) lx_z, tew_z, zew_z ,rvir*r0 
      write(12,15) 
      write(12,16)
      write(12,17)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# x : Lx =',g13.7,'[h^-2 erg/s]  Tew =  ',g13.7,
     &    '[keV]  Zew = ',g13.7,'[Zsun]  within R =',f7.4,' [1/h Mpc]')
 13   format ('# y : Lx =',g13.7,'[h^-2 erg/s]  Tew =  ',g13.7,
     &    '[keV]  Zew = ',g13.7,'[Zsun]  within R =',f7.4,' [1/h Mpc]')
 14   format ('# z : Lx =',g13.7,'[h^-2 erg/s]  Tew =  ',g13.7,
     &    '[keV]  Zew = ',g13.7,'[Zsun]  within R =',f7.4,' [1/h Mpc]')
 15   format ('# Columns:')
 16   format ('# r_l r_m r_r  r_x xsbp_x tew_x zew_x nx  xsbp_y ',
     &     'tew_y zew_y ny  r_z xsbp_z tew_z zew_z nz')
 17   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
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
        if ( nxbin(i) .gt. 0. .and. 
     +       nybin(i) .gt. 0. .and. 
     +       nzbin(i) .gt. 0.  ) then
          write(12,18)
     &         rl*r0, rmid, rr*r0, 
     &         rxbin(i),xsbp_xbin(i),tew_xbin(i),zew_xbin(i),nxbin(i),
     &         rybin(i),xsbp_ybin(i),tew_ybin(i),zew_ybin(i),nybin(i),
     &         rzbin(i),xsbp_zbin(i),tew_zbin(i),zew_zbin(i),nzbin(i)
	endif
      enddo
 18   format(18(2x,g14.7))
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      
      return
      end
c
c     -----------------------------------------------------------------
      subroutine XSpectra ( rs_spec_fname )
c     Output X-ray spectra in annuli (rl<r<rr) [erg/s/cm^2/keV]
c 
c          Use the smallest grid as a basis for uniform grid
c          Compute mass within an uniform cell in the projected plane
c          Compute the surface mass density [Msun Mpc^-2]
c
c     CAUTION : nsize<=1024 for X-ray analysis involving X-ray spectra 
c
c     -----------------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      parameter ( nxspec = 1100, nrsize = 100 )
      integer iflag
      real*8 ex
      real ex1(nxspec,nrsize), ex2(nxspec,nrsize),ex3(nxspec,nrsize) 
      common / XSPEC1 / ex1, ex2, ex3
      character*256 fname, fnh, rs_spec_fname
      character*5 chbox
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 

      real EeV1, EeV2, dE
      common / RS03 / Xpri, Xsol,EeV1,EeV2,nrs,dene,SOLA,defactor2
      common / RS04 / NET, NZT, NTT, dE, RSPEC
      
c      integer NRBIN
c      parameter ( NRBIN = 6 )
c      real rad(NRBIN)
c      DATA rad/0.,100.,200.,300.,400.,500./ ! [kpc]
     
c
c...  define parameters for constructing X-ray spectra
c
      
      ! radial bins
      NRBIN = 10
      rmin = 10.0 / 1000.0 / r0 * hubble   ! 10 kpc
      rmax = 500.0 / 1000.0 / r0 * hubble  ! 500 kpc
      rlmin = log10(rmin)
      rlmax = log10(rmax)
      dr = ( rmax - rmin ) / NRBIN
      drl = (rlmax - rlmin) / float(NRBIN) 

      do ir = 1 , NRBIN
         do ie = 1, NET
            ex1(ie,ir) = 0.0
            ex2(ie,ir) = 0.0
            ex3(ie,ir) = 0.0
         enddo
      enddo

      ! file name
      lpath = index(path    , ' ') - 1
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)

c
c...  Find the halo center
c
      deltavir  = 340.
      deltacrit = 500.
      ifindcenter = 3 

      call Find_Halo_Center ( xc, yc, zc, ifindcenter )
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
c...  set up var(1,*) and ref(1,*)
c
      call Get_MaxLevelNow() 
      if ( iflag .eq. 1 ) then 
         call Assign_Density ( MinLevel , MaxLevelNow , 0 ) 
      else if ( iflag .eq. 2 ) then
         call Assign_DMDensity ( MinLevel , MaxLevelNow , 0 ) 
      endif

c
c...  define useful variables
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

c...  normalization constants
      Ynow = 0.2496
      t_fact = gamma1 * T_0 * a2i
      ex_fact = 3.89535e-9 * (1.0 - 0.5 * Ynow) * (1.0 - Ynow) * 
     &   a5i * hubble**3 * Om0**2 * r0
      
      call Read_RS_Spec ( rs_spec_fname )
      write(*,*) 'Finished reading the RS spectra...'

c
c...  Loop over pixel-by-pixel
c

      lstart = MinLevel
      lend   = MaxLevelNow
      lstart = lstart + 1
      
      DO Level = lstart , lend
        CellVolume = 1.0 * 2.0**(-3.0*Level)
        CellArea = 1.0 * 2.0**(-2.0*Level)
        cs  = CellSize(Level)
        cs2 = CellSize(Level)/2.       
  
        IF ( Level .eq. MinLevel ) THEN
        ELSE

          nLevel = iNOLL(Level)
          write(*,*) 'Level = ',Level,' nLevel = ',nLevel
          call Select_Cells ( Level , nLevel )
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              if ( iOctCh(ic) .eq. nil ) then
                call Ps ( ic , Posx , Posy , Posz )
                call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz)
                rx = sqrt( (Posy-yc)**2 + (Posz-zc)**2 )
                ry = sqrt( (Posx-xc)**2 + (Posz-zc)**2 )
                rz = sqrt( (Posx-xc)**2 + (Posy-yc)**2 )

                ixbin = max ( int((log10(rx) - rlmin)/drl) + 1 , 0 )
                iybin = max ( int((log10(ry) - rlmin)/drl) + 1 , 0 )
                izbin = max ( int((log10(rz) - rlmin)/drl) + 1 , 0 )
      
c                ixbin = max ( int((rx-rmin)/dr) + 1 , 0 )
c                iybin = max ( int((ry-rmin)/dr) + 1 , 0 )
c                izbin = max ( int((rz-rmin)/dr) + 1 , 0 )
      
                if ( ixbin .le. NRBIN .or. iybin .le. NRBIN .or. 
     &               izbin .le. NRBIN ) then
                   do ie = 1, NET
                      Tcell = t_fact *  hvar(8,ic) / hvar(1,ic)
                      Zcell = (hvar(izIa,ic)+hvar(izII,ic)) / hvar(1,ic)
                      Zcell = Zcell / 0.01989  
c                      Zcell = 0.30
                      Tlog = log10(Tcell)
                      Zlog = log10(Zcell)
                      ex = ex_fact * hvar(1,ic)**2 * 
     &                     RS_Spec( ie, Tlog, Zlog ) * CellVolume
                      if ( ixbin .le. NRBIN ) then
                         ex1(ie,ixbin) = ex1(ie,ixbin) + ex
                      endif
                      if ( iybin .le. NRBIN ) then
                         ex2(ie,iybin) = ex2(ie,iybin) + ex
                      endif
                      if ( izbin .le. NRBIN ) then
                         ex3(ie,izbin) = ex3(ie,izbin) + ex
                      endif
                   enddo
                endif

              endif
            enddo !ic2
          enddo  !ic1

        ENDIF
      ENDDO

c
c...  convolve in the X-ray instrumental responses
c
c      do ir = 1 , NRBIN
c         do ie = 1, NET
c            ex(ie) = ex1(ie,ir)
c         enddo
c
c         q_d1(ix,iy) = rs_det( ex,  NET, EeV1, dE )
c         write(*,*) p_d1(ix,iy), q_d1(ix,iy)
c  
c         do ie = 1, NET
c            ex(ie) = ex2(ie,ix,iy) 
c            vex(ie) = vex2(ie,ix,iy)
c         enddo
c         p_d2(ix,iy) = rs_det( vex, NET, EeV1, dE )
c         q_d2(ix,iy) = rs_det( ex,  NET, EeV1, dE )
c  
c         do ie = 1, NET
c            ex(ie) = ex3(ie,ix,iy) 
c            vex(ie) = vex3(ie,ix,iy)
c         enddo  
c         p_d3(ix,iy) = rs_det( vex, NET, EeV1, dE )
c         q_d3(ix,iy) = rs_det( ex,  NET, EeV1, dE )
c
c      enddo
c
c      endif

c
c.... X-ray spectra in annuli (rl<r<rr) [erg/s/cm^2/keV]
c     
      do ir = 1 , NRBIN
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
c            rl = 10.**(rmin + float(i-1)*dr)
            rl = 10.**(rlmin + float(i-1)*drl)
         endif
c         rr = 10.**(rmin + float(i)*dr)
         rr = 10.**(rlmin + float(i)*drl)
         area = pi * (rr**2 - rl**2)  
         area2 = pi * rr**2
         
         write(*,*) rl, rr, area
         do ie = 1, NET
            ex1(ie,ir) = ex1(ie,ir) / area / (dE/1000.)
            ex2(ie,ir) = ex2(ie,ir) / area / (dE/1000.)
            ex3(ie,ir) = ex3(ie,ir) / area / (dE/1000.)
         enddo
      enddo

c
c...  Output ASCII data to the file (Total surface mass density)
c

      rmin = rmin * 1000.0 * r0 / hubble  ! kpc
      rmax = rmax * 1000.0 * r0 / hubble  ! kpc
      dr = dr * 1000.0 * r0 / hubble      ! kpc

      lpath = index(path    , ' ') - 1
      fnh = '/xspec'
      lfnh = index(fnh    , ' ') - 1

      fname = path(1:lpath)//fnh(1:lfnh)//'_a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_x.dat'
      nfname = index ( fname , ' ' ) - 1
      write(*,*) 'writing...'
      open ( 20 , file = fname(1:nfname) )
      write(20,30) aexpn
      write(20,31) NRBIN, rmin, rmax, dr
      write(20,31) NET, EeV1, EeV2, dE
      do ir = 1, NRBIN
         do ie = 1, NET
            write(20,32) ex1(ie,ir)
         enddo
      enddo
      close ( 20 )
      write(*,*) 'done.'

      fname = path(1:lpath)//fnh(1:lfnh)//'_a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_y.dat'
      nfname = index ( fname , ' ' ) - 1
      write(*,*) 'writing...'
      open ( 20 , file = fname(1:nfname) )
      write(20,30) aexpn
      write(20,31) NRBIN, rmin, rmax, dr
      write(20,31) NET, EeV1, EeV2, dE
      do ir = 1, NRBIN
         do ie = 1, NET
            write(20,32) ex2(ie,ir)
         enddo
      enddo
      close ( 20 )
      write(*,*) 'done.'

      fname = path(1:lpath)//fnh(1:lfnh)//'_a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_z.dat'
      nfname = index ( fname , ' ' ) - 1
      write(*,*) 'writing...'
      open ( 20 , file = fname(1:nfname) )
      write(20,30) aexpn
      write(20,31) NRBIN, rmin, rmax, dr
      write(20,31) NET, EeV1, EeV2, dE
      do ir = 1, NRBIN
         do ie = 1, NET
            write(20,32) ex3(ie,ir)
         enddo
      enddo
      close ( 20 )
      write(*,*) 'done.'

 30   format (f7.4)
 31   format (i5,2x,3(f10.3,1x))
 32   format (g13.7)

      return
      end

c
c     -----------------------------------------
      subroutine Output_XrayVir ( ifindcenter )
c     -----------------------------------------
c     compute the mass & emission weighted temperature and 
c     metallicity within an aperture defined by various virial 
c     radius definitions.
c     all quantities are in code units
c
c     output: mass of gas at T<Tmax, total mass of gas within rvir (masses are in code units)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      character*256 fname, rs_dir, rs_fname, obs
      integer nfpn, lfrs, lobs
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      real rA(12,6), xA(12,6), yA(12,6), zA(12,6), ahA(12,6)
      integer iA(12)
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
c
c...  define the cooling function 
c
      idet = 1
      rs_dir = '/dworkin/home/daisuke/Raymond/DAT/'
      lfrs = index(rs_dir , ' ') - 1
     
      ! Chandra 
      if ( idet .eq. 1 ) then
         rs_fname  = rs_dir(1:lfrs)//
     &        'rs_det_chandra_acis-s-bi_0.5-7.0keV.dat'
         obs = 'chandra_0.5-7.0keV'
      endif

      ! BeppoSAX
      if ( idet .eq. 2 ) then
         rs_fname  = rs_dir(1:lfrs)//'rs_det_sax_mecs_2-10keV.dat'
         obs = 'sax_2-10keV'
      endif

      ! ROSAT
      if ( idet .eq. 3 ) then
         rs_fname  = rs_dir(1:lfrs)//'rs_det_rosat_pspc_0.5-2.4keV.dat'
         obs = 'rosat_0.5-2.4keV'
      endif

      lobs = index(obs    , ' ') - 1

      
c...  Read in the Raymond-Smith cooling function from file
      irs = 0    ! 0 = RS_det(T,Z), 1 = RS_cgs(T,Z), 2 = RS(T,Z=0.3Msun),
      if ( irs .eq. 2 ) then
        call read_rs()
      else
        if ( irs .eq. 0 ) call Read_RSdet_Table ( 0, rs_fname )
        if ( irs .eq. 1 ) call Read_RSdet_Table ( 1, rs_fname )
        write(*,*) 'Finished reading the RS spectra...'
      endif

c
c...  Define the output file name
c
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)
      
      lpath = index(path    , ' ') - 1
      fname  = path(1:lpath)//'/xrayvir'//'_'//'a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'_'//obs(1:lobs)//'.dat '

c
c...  Initialize arrays
c
      do i = 1, 12
         do j = 1, 6
            rA(i,j) = 0.0
            xA(i,j) = 0.0
            yA(i,j) = 0.0
            zA(i,j) = 0.0
         enddo
      enddo

      do i = 1, 12
         do j = 1, 5
            ahA(i,j) = 0.0
         enddo
         iA(i) = 0
      enddo

c
c...  Find halo center and radius
c
      call Get_MaxLevelNow ()
      call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      
      oma = Om0 / ( Om0 + Oml0 * aexpn**3 )
      xf = oma - 1.
      ahA(7,1) = ( 18.*pi**2 + 82.*xf - 39.*xf**2 )
      ahA(1,1) = ahA(7,1) / (1.+xf)  
      call Find_Halo_Radius ( xc, yc, zc, ahA(1,1), ahA(7,1),
     &     iA(1), ahA(1,2), ahA(1,3), ahA(1,4), ahA(1,5), ahA(1,6), 
     &     iA(2), ahA(7,2), ahA(7,3), ahA(7,4), ahA(7,5), ahA(7,6), 
     &     ifindcenter )

      ahA(2,1) = 180.
      ahA(8,1) = 180.
      call Find_Halo_Radius ( xc, yc, zc, ahA(2,1), ahA(8,1),
     &     iA(3), ahA(2,2), ahA(2,3), ahA(2,4), ahA(2,5), ahA(2,6), 
     &     iA(4), ahA(8,2), ahA(8,3), ahA(8,4), ahA(8,5), ahA(8,6), 
     &     ifindcenter )

      ahA(3,1) = 200.
      ahA(9,1) = 200.
      call Find_Halo_Radius ( xc, yc, zc, ahA(3,1), ahA(9,1),
     &     iA(5), ahA(3,2), ahA(3,3), ahA(3,4), ahA(3,5), ahA(3,6), 
     &     iA(6), ahA(9,2), ahA(9,3), ahA(9,4), ahA(9,5), ahA(9,6), 
     &     ifindcenter )

      ahA(4,1) = 500.
      ahA(10,1) = 500.
      call Find_Halo_Radius ( xc, yc, zc, ahA(4,1), ahA(10,1),
     &     iA(7), ahA(4,2), ahA(4,3), ahA(4,4), ahA(4,5), ahA(4,6), 
     &     iA(8), ahA(10,2), ahA(10,3), ahA(10,4), ahA(10,5), ahA(10,6), 
     &     ifindcenter )

      ahA(5,1) = 2500.
      ahA(11,1) = 2500.
      call Find_Halo_Radius ( xc, yc, zc, ahA(5,1), ahA(11,1),
     &     iA(9), ahA(5,2), ahA(5,3), ahA(5,4), ahA(5,5), ahA(5,6), 
     &     iA(10), ahA(11,2), ahA(11,3), ahA(11,4), ahA(11,5),ahA(11,6), 
     &     ifindcenter )

      ahA(6,1) = 5000.
      ahA(12,1) = 5000.
      call Find_Halo_Radius ( xc, yc, zc, ahA(6,1), ahA(12,1),
     &     iA(11), ahA(6,2), ahA(6,3), ahA(6,4), ahA(6,5),ahA(6,6), 
     &     iA(12), ahA(12,2), ahA(12,3), ahA(12,4), ahA(12,5),ahA(12,6), 
     &     ifindcenter )


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
      
c...            compute temperature and metallicity of the cell
                if ( hvar(1,ic) .gt. 0. ) then
                   Tcell = gamma1 * hvar(8,ic) / hvar(1,ic)
                   Tcell = T_0 * Tcell * a2i
                   Zcell = hvar(izIa,ic) + hvar(izII,ic)
                   Zcell = Zcell / 0.01989 / hvar(1,ic)
                else
                   write(*,*) 'bad gas density in cell', ic,hvar(1,ic)
                endif
                
c...            compute gas mass and emissivity of the cell
                if ( irs .eq. 2 ) then
                  call lambda_rs( Tcell, lambda )   
                  ex = hvar(1,ic)**2 * lambda * CellVolume
                else 
                  Tlog = log10(Tcell)
                  Zlog = log10(Zcell)
                  gm = hvar(1,ic) * CellVolume
                  ex = hvar(1,ic)**2 * RS(Tlog,Zlog) * CellVolume
                endif

                tm = Tcell * gm
                zm = Zcell * gm
                tew = Tcell * ex
                zew = Zcell * ex
                                
c...            temperature and metallicity within the virial regions
                do i = 1, 12
                   ! radial
                   if ( rr .le. ahA(i,2) ) then
                      call Compute_Xray_Virial( i, rA, 
     &                  tm, zm, tew, zew, gm, ex )
                   endif
                   ! x-axis
                   if ( rx .le. ahA(i,2) ) then
                      call Compute_Xray_Virial( i, xA, 
     &                  tm, zm, tew, zew, gm, ex )
                   endif
                   ! y-axis
                   if ( ry .le. ahA(i,2) ) then
                      call Compute_Xray_Virial( i, yA, 
     &                  tm, zm, tew, zew, gm, ex )
                   endif
                   ! z-axis
                   if ( rz .le. ahA(i,2) ) then
                      call Compute_Xray_Virial( i, zA, 
     &                  tm, zm, tew, zew, gm, ex )
                   endif
                enddo

              endif
            enddo
          enddo
        ENDIF
      ENDDO

c...  mass and emission weighted temperature and metallicity

      do i = 1, 12
         rA(i,1) = rA(i,1) / rA(i,5) / 1.160485e+7  
         rA(i,2) = rA(i,2) / rA(i,5)                
         rA(i,3) = rA(i,3) / rA(i,6) / 1.160485e+7  
         rA(i,4) = rA(i,4) / rA(i,6)                

         xA(i,1) = xA(i,1) / xA(i,5) / 1.160485e+7  
         xA(i,2) = xA(i,2) / xA(i,5)                
         xA(i,3) = xA(i,3) / xA(i,6) / 1.160485e+7  
         xA(i,4) = xA(i,4) / xA(i,6)                

         yA(i,1) = yA(i,1) / yA(i,5) / 1.160485e+7  
         yA(i,2) = yA(i,2) / yA(i,5)                
         yA(i,3) = yA(i,3) / yA(i,6) / 1.160485e+7  
         yA(i,4) = yA(i,4) / yA(i,6)                

         zA(i,1) = zA(i,1) / zA(i,5) / 1.160485e+7  
         zA(i,2) = zA(i,2) / zA(i,5)                
         zA(i,3) = zA(i,3) / zA(i,6) / 1.160485e+7  
         zA(i,4) = zA(i,4) / zA(i,6)                

         ahA(i,2) = ahA(i,2) * r0
      enddo

c...  Writing data to the file
      nfpn = index ( fname , ' ' ) - 1
      open ( 12 , file = fname(1:nfpn) )
      write(12,19) xc , yc , zc  
      write(12,20) ifindcenter
      write(12,21)
      do i = 1, 6
         write(12,23) ahA(i,1), ahA(i,2), ahA(i,3), ahA(i,4), ahA(i,5), 
     &        ahA(i,6), rA(i,1), rA(i,2), rA(i,3), rA(i,4), 
     &        xA(i,1), xA(i,2), xA(i,3), xA(i,4),
     &        yA(i,1), yA(i,2), yA(i,3), yA(i,4),
     &        zA(i,1), zA(i,2), zA(i,3), zA(i,4)
      enddo

      write(12,22)
      do i = 7, 12 
         write(12,23) ahA(i,1), ahA(i,2), ahA(i,3), ahA(i,4), ahA(i,5), 
     &        ahA(i,6), rA(i,1), rA(i,2), rA(i,3), rA(i,4), 
     &        xA(i,1), xA(i,2), xA(i,3), xA(i,4),
     &        yA(i,1), yA(i,2), yA(i,3), yA(i,4),
     &        zA(i,1), zA(i,2), zA(i,3), zA(i,4)
      enddo
 19   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 20   format ('# ifindcenter = ',i5)
 21   format ('# Dvir(mean)  Rvir [1/h Mpc]  Mvir        Mdm     ',
     &     '     Mgas       Mst [1/h Msun]   ', 
     &     '     <tm>  <zm> <tew> <zew>(<Rvir) for r,x,y,z projections')
 22   format ('# Dcrit(crit) Rcrit [1/h Mpc] Mcrit       Mdm     ', 
     &     '     Mgas       Mst [1/h Msun]   ',
     &     '     <tm>  <zm> <tew> <zew>(<Rvir) for r,x,y,z projections')
 23   format (f9.3,4x,f7.4,3x,20(g13.5,1x))
      close ( 12 )
      write(*,*) 'done writing X-ray data. exiting...'
    
      return
      end

c
c     -----------------------------------------
      subroutine Compute_Xray_Virial ( 
     &     i, array, tm, zm, tew, zew, gm, ex ) 
c     -----------------------------------------
c
      real array(12,6)

      array(i,1) = array(i,1) + tm
      array(i,2) = array(i,2) + zm
      array(i,3) = array(i,3) + tew
      array(i,4) = array(i,4) + zew
      array(i,5) = array(i,5) + gm
      array(i,6) = array(i,6) + ex

      return
      end

c
c     ----------------------------------------------------------
      subroutine Compute_Xray ( xc, yc, zc, rvir )
c     ----------------------------------------------------------
c     compute Mgas, Lx, Tew, Zew for an object centered 
c     at (xc,yc,zc) and radius rvir
c     all quantities are in code units
c
c     input : xc, yc, zc, rvir, Tmax
c     output: mass of gas at T<Tmax, total mass of gas within rvir (masses are in code units)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real   xc, yc, zc, rvir
c
      character*256 fname
      real   Tcell, Tmax, Mgas, Tew, Zew, coldmass, norm_rs
      double precision mp, mpc_2_cm, lambda, Lx, norm_lx, ex_sum,ex_cell    
c
      Mgas = 0.
      Lx   = 0.
      Tew  = 0.
      Zew  = 0.
      ex_sum = 0.
      coldmass = 0.
      Tmax = 100000.0

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      
c...  Read in the Raymond-Smith cooling function from file
      call read_rs()
      call Get_MaxLevelNow ()

c...  Loop over pixel-by-pixel
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
c...                 compute gas mass of the cell
                     cellmass_gas = hvar(1,ic) * CellVolume
                     cellmass_gas = cellmass_gas * rho0 / a3 * r0**3
                     Mgas = Mgas + cellmass_gas  

c...                 compute temperature of the cell
                     if ( hvar(1,ic) .gt. 0. ) then
                        Tcell = gamma1 * hvar(8,ic)/hvar(1,ic)
                        Tcell = T_0 * Tcell * a2i
                        Zcell = hvar(izIa,ic) + hvar(izII,ic)
                        Zcell = Zcell / 0.01989 / hvar(1,ic)
                     else 
                        write(*,*) 'bad gas density in cell',
     &                       ic,hvar(1,ic)
                     endif

c...                 compute emissivity
                     call lambda_rs(Tcell,lambda)   
                     ex_cell = hvar(1,ic)**2 * lambda * CellVolume

c...                 add up values of the cell 
                     ex_sum  = ex_sum + ex_cell
                     Lx  = Lx  + ex_cell
                     Tew = Tew + Tcell * ex_cell
                     Zew = Zew + Zcell * ex_cell

                     if ( Tcell .lt. Tmax ) then
                        coldmass = coldmass + cellmass_gas 
                     endif
        
c                     write(*,'(2(e12.5))') Tcell, lambda
c                     write(*,'(i7,4(e12.5))') ic, ex_cell, Lx, Tew, Zew

                  endif
                endif
              enddo
            enddo
          ENDIF
        ENDDO

        
c...  compute emission-weighted temperature and metallicity
      Tew = Tew / ex_sum
      Zew = Zew / ex_sum

c...  convert code units to physical units            
      Tew = Tew * T_0 * a2i / (1.160485e+7)
          
      zsim = 1.0/aexpn - 1.0
      fy = (1.0-0.5*Y_p)*(1.0-Y_p)
      mp = 1.6726d-24 
      mpc_2_cm = mpc * 100000.0
      norm_rs = 10.**(-23)
      norm_lx = fy * (rho0/a3*sun_mass/mp)**2 / mpc_2_cm**3 * r0**3
c      norm_lx = norm_lx * (4.0*pi*(1+zsim)**4)/(4.0*(pi/21600)**2)
      norm_lx = norm_lx * norm_rs / hubble
      Lx = Lx * norm_lx 

      write(*,'("Tew " g10.4" [keV]")')  Tew
      write(*,'("Zew " g10.4" [Zsun]")') Zew
      write(*,'("Lx " g12.4 " [h^-2 erg s^-1]")') Lx

      write(*,'("Mcold[T<10^5K] " e" [Msun] ")') coldmass
      
      return
      end

c
c...  routines dealing with the RS code
c
c     ---------------------------------------------
      subroutine Read_RSdet_Table ( icgs, ftable )
c     ---------------------------------------------
c     
c     read in previously stored RS emissivity lookup table
c     input: ftable - table filename (must end with a space)
c
c     if icgs = 1 : the RS table in the CGS units [erg s^-1 cm^3]
c     if icgs = 0 : the RS table in the detector units [cnts s^-1 cm^5]
c
      parameter ( NZMAX = 100 , NTMAX = 1000 )
      integer NET, NZT, NTT, icgs
      real Zlmin, Zlmax, dlZ, Tlmin, Tlmax, dlT
      real Xpri, Xsol, EeV1, EeV2, dE, dene
      real*8 defactor2
      integer nrs
      real*8 RST(NTMAX,NZMAX)
      real SOLA(12)
      common / RS01 / RST
      common / RS02 / Zlmin, Zlmax, dlZ, Tlmin, Tlmax, dlT
      common / RS03 / Xpri,Xsol,EeV1,EeV2,nrs,dene,SOLA,defactor2
      common / RS04 / NET, NZT, NTT, dE, RSPEC

      character*256 ftable 
      character*20 dummy
      integer lname 

      lname = index(ftable, ' ') - 1
      open ( 30 , file = ftable(1:lname) )
      do i = 1 , 4
        read(30,*) dummy
      enddo
      read(30,*) NET, EeV1, EeV2, dE
      read(30,*) NZT, Zlmin, Zlmax, dlZ
      read(30,*) NTT, Tlmin, Tlmax, dlT

      if ( NZT .gt. NZMAX ) then 
       write(*,*) '* error in Read_RSdet_Table: NZT > NZMAX:',NZT, NZMAX
       write(*,*) '  increase NZMAX and rerun'       
       stop
      endif
      if ( NTT .gt. NTMAX ) then 
       write(*,*) '* error in Read_RSdet_Table: NTT > NTMAX:',NTT, NTMAX
       write(*,*) '  increase NTMAX and rerun'       
       stop
      endif
      if ( NZT .le. 1 ) then 
        write(*,*) '* error in Read_RSdet_Table: NZT <= 1:',NZT
        stop
      endif
      if ( NTT .le. 1 ) then 
        write(*,*) '* error in Read_RSdet_Table: NTT <= 1:',NTT
        stop
      endif
      
      do iz = 1, NZT
         do it = 1, NTT
           read(30,*) Zlog, Tlog, rs_det, rs_cgs
           if (icgs .eq. 1 ) RST(it,iz) = rs_cgs
           if (icgs .eq. 0 ) RST(it,iz) = rs_det
         enddo
      enddo
      close ( 30 )

      write(*,10) NET, EeV1, EeV2, dE 
      write(*,11) NZT, Zlmin, Zlmax, dlZ
      write(*,12) NTT, Tlmin, Tlmax, dlT
      if ( icgs. eq. 1 ) then 
        write(*,*) 'RS table in CGS units [erg s^-1 cm^3]'
      endif
      if ( icgs. eq. 2 ) then
        write(*,*) 'RS table in detector units [cnts s^-1 cm^5]'
      endif

 10   format('NET= ',i5,' EeV1 =',f8.3,' EeV2 = ',f8.3,' dE =',f8.3)
 11   format('NZT= ',i5,' Zlmin = ',f6.3,' Zlmax = ',f6.3,' dlZ =',f6.3)
 12   format('NTT= ',i5,' Tlmin = ',f6.3,' Zlmax = ',f6.3,' dlT =',f6.3)

      return
      end
c
c     -----------------------------------
      subroutine Read_RS_Table ( ftable )
c     -----------------------------------
c     
c     read in previously stored RS emissivity lookup table
c     input: ftable - table filename (must end with a space)
c
      parameter ( NZMAX = 100 , NTMAX = 1000 )
      integer NZT, NTT
      real Zlmin, Zlmax, Tlmin, Tlmax
      real Xpri, Xsol, EeV1, EeV2, dene
      real*8 defactor2
      integer nrs
      real*8 RST(NTMAX,NZMAX)
      real SOLA(12)
      common / RS01 / RST
      common / RS02 / Zlmin, Zlmax, dlZ, Tlmin, Tlmax, dlT
      common / RS03 / Xpri,Xsol,EeV1,EeV2,nrs,dene,SOLA,defactor2
      character*256 ftable 
      character*20 dummy
      integer lname 

      lname = index(ftable, ' ') - 1
      open ( 30 , file = ftable(1:lname) )
      do i = 1 , 12 
        read(30,*) dummy
      enddo
      read(30,*) NZT, NTT
      if ( NZT .gt. NZMAX ) then 
        write(*,*) '* error in Read_RS_Table: NZT > NZMAX:',NZT, NZMAX
        write(*,*) '  increase NZMAX and rerun'       
        stop
      endif
      if ( NTT .gt. NTMAX ) then 
        write(*,*) '* error in Read_RS_Table: NTT > NTMAX:',NTT, NTMAX
        write(*,*) '  increase NTMAX and rerun'       
        stop
      endif
      if ( NZT .le. 1 ) then 
        write(*,*) '* error in Read_RS_Table: NZT <= 1:',NZT
        stop
      endif
      if ( NTT .le. 1 ) then 
        write(*,*) '* error in Read_RS_Table: NTT <= 1:',NTT
        stop
      endif
      read(30,*) Zlmin, Zlmax, dlZ, Tlmin, Tlmax, dlT
      read(30,*) Xpri, Xsol, EeV1, EeV2, dene, nrs, defactor2
      read(30,*) (SOLA(ie), ie=1,12)
      read(30,*) ((RST(it,iz),it=1,NTT),iz=1,NZT)      
      close ( 30 )

      write(*,10) NET, EeV1, EeV2, dE 
      write(*,11) NZT, Zlmin, Zlmax, dlZ
      write(*,12) NTT, Tlmin, Tlmax, dlT

 10   format('NET= ',i5,' EeV1 =',f8.3,' EeV2 = ',f8.3,' dE =',f8.3)
 11   format('NZT= ',i5,' Zlmin = ',f6.3,' Zlmax = ',f6.3,' dlZ =',f6.3)
 12   format('NTT= ',i5,' Tlmin = ',f6.3,' Zlmax = ',f6.3,' dlT =',f6.3)

      return
      end
c
c     ----------------------------------
      real*8 function RS ( Tlog, Zlog )
c     ---------------------------------- 
c
c     uses CIC to interpolate 2D RST table for Zlog, Tlog
c     (log10 of Z in solar units and T in K)
c
c     WARNING: RST definition must be identical to that in Tabulate routine
c
      parameter ( NZMAX = 100 , NTMAX = 1000 )
      real Zlmin, Zlmax, Tlmin, Tlmax, dlZ, dlT
      real Zlog, Tlog
      real*8 RST(NTMAX,NZMAX)
      real   SOLA(12) 
      real*8 tt, zt, td, zd
c
      real*8 defactor2
      common / RS01 / RST
      common / RS02 / Zlmin, Zlmax, dlZ, Tlmin, Tlmax, dlT
      common / RS03 / Xpri, Xsol,EeV1,EeV2,nrs,dene,SOLA,defactor2
      common / RS04 / NET, NZT, NTT, dE, RSPEC
      common/PARAMS/ NJ(12),ABUNJ(12),ABUND,BINMIN,BINSYZ,NBIN
      common/com/cnc(12,30),PTOT(12,220),abin(1000),bin(1000)

      ! no X-ray emission from gas with T<Tmin
      if ( Tlog .lt. Tlmin ) then 
        RS = 0.
        return
      endif 

      izb1 = (Zlog - Zlmin) / dlZ + 1
      itb1 = (Tlog - Tlmin) / dlT + 1
      izb2 = izb1 + 1
      itb2 = itb1 + 1

      IF ( izb1.ge.1 .and. itb1.ge.1 .and. 
     &     Zlog .lt. Zlmax           .and. 
     &     Tlog .lt. Tlmax ) THEN 
        Zl1  = Zlmin + dlZ * (izb1 - 1)
        Tl1  = Tlmin + dlT * (itb1 - 1)
        Zl2  = Zl1 + dlZ 
        Tl2  = Tl1 + dlT
        tt  = (Tl2 - Tlog)/dlT
        td  = 1.0 - tt
        zt  = (Zl2 - Zlog)/dlZ
        zd  = 1.0 - zt 
       
        RS = tt * zt * RST(itb1,izb1) + 
     &       td * zt * RST(itb2,izb1) + 
     &       tt * zd * RST(itb1,izb2) + 
     &       td * zd * RST(itb2,izb2)
c        write(*,*) 'rs :',Zl1,Tl1,itb1,izb1,RST(itb1,izb1)
      ELSE ! if outside the range
        if ( izb1 .lt. 1 ) then
           izb1 = 1
           Zlog = Zlmin
        endif
        if ( itb1 .lt. 1 ) then
           itb1 = 1
           Tlog = Tlmin
        endif
        if ( Zlog .gt. Zlmax ) then
           izb1 = (Zlmax - Zlmin) / dlZ + 1
           Zlog = Zlmax
        endif
        if ( Tlog .gt. Tlmax ) then
           itb1 = (Tlmax - Tlmin) / dlT + 1
           Tlog = Tlmax
        endif
        izb2 = izb1 + 1
        itb2 = itb1 + 1
        Zl1  = Zlmin + dlZ * (izb1 - 1)
        Tl1  = Tlmin + dlT * (itb1 - 1)
        Zl2  = Zl1 + dlZ 
        Tl2  = Tl1 + dlT
        tt  = (Tl2 - Tlog)/dlT
        td  = 1.0 - tt
        zt  = (Zl2 - Zlog)/dlZ
        zd  = 1.0 - zt 
       
        RS = tt * zt * RST(itb1,izb1) + 
     &       td * zt * RST(itb2,izb1) + 
     &       tt * zd * RST(itb1,izb2) + 
     &       td * zd * RST(itb2,izb2)
        if ( RS .lt. 0.0 ) then
           write(*,*) 'error : RS < 0.0' 
           write(*,*) tt, td, zt, zd, Tlog, Zlog
           write(*,*) Zl2, Zlog, dlZ
           write(*,*) itb1, izb1, itb2, izb2
           write(*,*) RST(itb1,izb1),RST(itb2,izb1)
           write(*,*) RST(itb1,izb2),RST(itb2,izb2)
           write(*,*) 'RS = ',RS
        endif
      ENDIF

      return
      end
c
c     -----------------------------------
      subroutine Read_RS_Spec ( ftable )
c     -----------------------------------
c     
c     read in previously stored RS emissivity lookup table
c     input: ftable - table filename (must end with a space)
c     
c     units : n_e * n_H * 1.0e-23 * RS_Spec [erg s^-1 cm^-3] 
c
      parameter ( NZMAX = 100 , NTMAX = 1000, NEMAX = 1000 )
      integer NET, NZT, NTT
      real Zlmin, Zlmax, Tlmin, Tlmax
      real Xpri, Xsol, EeV1, EeV2, dene
      real*8 defactor2
      integer nrs
      real*8 RST(NTMAX,NZMAX),RSPEC(1000,100,20)
      real SOLA(12)
      common / RS01 / RST
      common / RS02 / Zlmin, Zlmax, dlZ, Tlmin, Tlmax, dlT
      common / RS03 / Xpri,Xsol,EeV1,EeV2,nrs,dene,SOLA,defactor2
      common / RS04 / NET, NZT, NTT, dE, RSPEC
      character*256 ftable 
      character*20 dummy
      integer lname 

      lname = index(ftable, ' ') - 1
      open ( 30 , file = ftable(1:lname) )
      do i = 1 , 13 
        read(30,*) dummy
      enddo
      read(30,*) NET, EeV1, EeV2, dE
      read(30,*) NZT, Zlmin, Zlmax, dlZ
      read(30,*) NTT, Tlmin, Tlmax, dlT

      if ( NZT .gt. NZMAX ) then 
        write(*,*) '* error in Read_RS_Table: NZT > NZMAX:',NZT, NZMAX
        write(*,*) '  increase NZMAX and rerun'       
        stop
      endif
      if ( NTT .gt. NTMAX ) then 
        write(*,*) '* error in Read_RS_Table: NTT > NTMAX:',NTT, NTMAX
        write(*,*) '  increase NTMAX and rerun'       
        stop
      endif
      if ( NZT .le. 1 ) then 
        write(*,*) '* error in Read_RS_Table: NZT <= 1:',NZT
        stop
      endif
      if ( NTT .le. 1 ) then 
        write(*,*) '* error in Read_RS_Table: NTT <= 1:',NTT
        stop
      endif
      read(30,*) Xpri, Xsol, EeV1, EeV2, dene, nrs, defactor2
      read(30,*) (SOLA(ie), ie=1,12)
      do it = 1 , NTT 
         do iz = 1 , NZT 
            read(30,*) (RSPEC(ii,it,iz),ii=1,NET)
         enddo
      enddo
      close ( 30 )

      write(*,10) NET, EeV1, EeV2, dE 
      write(*,11) NZT, Zlmin, Zlmax, dlZ
      write(*,12) NTT, Tlmin, Tlmax, dlT

 10   format('NET= ',i5,' EeV1 =',f8.3,' EeV2 = ',f8.3,' dE =',f8.3)
 11   format('NZT= ',i5,' Zlmin = ',f6.3,' Zlmax = ',f6.3,' dlZ =',f6.3)
 12   format('NTT= ',i5,' Tlmin = ',f6.3,' Zlmax = ',f6.3,' dlT =',f6.3)
      return
      end

c
c     ------------------------------------------
      real*8 function RS_Spec ( ie, Tlog, Zlog )
c     ------------------------------------------ 
c
c     uses CIC to interpolate 2D RS spectra at energy bin, ie, 
c     for Zlog, Tlog (log10 of Z in solar units and T in K)
c
c     WARNING: RST definition must be identical to that in Tabulate routine
c
      parameter ( NZMAX = 100 , NTMAX = 1000, NEMAX = 1000 )

      real Zlmin, Zlmax, Tlmin, Tlmax, dlZ, dlT
      real EeV1, EeV2, dE
      real Zlog, Tlog
      real*8 RST(NTMAX,NZMAX),RSPEC(1000,100,20)
      real   SOLA(12) 
      real*8 tt, zt, td, zd
c
      real*8 defactor2
      common / RS01 / RST
      common / RS02 / Zlmin, Zlmax, dlZ, Tlmin, Tlmax, dlT
      common / RS03 / Xpri, Xsol,EeV1,EeV2,nrs,dene,SOLA,defactor2
      common / RS04 / NET, NZT, NTT, dE, RSPEC
      common/PARAMS/ NJ(12),ABUNJ(12),ABUND,BINMIN,BINSYZ,NBIN
      common/com/cnc(12,30),PTOT(12,220),abin(1000),bin(1000)

      ! no X-ray emission from gas with T<Tmin
      if ( Tlog .lt. Tlmin ) then 
        RS = 0.
        return
      endif 

      izb1 = (Zlog - Zlmin) / dlZ + 1
      itb1 = (Tlog - Tlmin) / dlT + 1
      izb2 = izb1 + 1
      itb2 = itb1 + 1

      IF ( izb1.ge.1 .and. itb1.ge.1 .and. 
     &     Zlog .lt. Zlmax           .and. 
     &     Tlog .lt. Tlmax ) THEN 
        Zl1  = Zlmin + dlZ * (izb1 - 1)
        Tl1  = Tlmin + dlT * (itb1 - 1)
        Zl2  = Zl1 + dlZ 
        Tl2  = Tl1 + dlT
        tt  = (Tl2 - Tlog)/dlT
        td  = 1.0 - tt
        zt  = (Zl2 - Zlog)/dlZ
        zd  = 1.0 - zt 
 
        RS_Spec = tt * zt * RSPEC(ie,itb1,izb1) + 
     &            td * zt * RSPEC(ie,itb2,izb1) + 
     &            tt * zd * RSPEC(ie,itb1,izb2) + 
     &            td * zd * RSPEC(ie,itb2,izb2)
      ELSE ! use the upper/lower limit, if outside the range
        if ( izb1 .lt. 1 ) izb1 = 1
        if ( itb1 .lt. 1 ) itb1 = 1
        if ( Zlog .gt. Zlmax ) izb1 = (Zlmax - Zlmin) / dlZ + 1
        if ( Tlog .gt. Tlmax ) itb1 = (Tlmax - Tlmin) / dlT + 1
        izb2 = izb1 + 1
        itb2 = itb1 + 1
        Zl1  = Zlmin + dlZ * (izb1 - 1)
        Tl1  = Tlmin + dlT * (itb1 - 1)
        Zl2  = Zl1 + dlZ 
        Tl2  = Tl1 + dlT
        tt  = (Tl2 - Tlog)/dlT
        td  = 1.0 - tt
        zt  = (Zl2 - Zlog)/dlZ
        zd  = 1.0 - zt 
       
        RS_Spec = tt * zt * RSPEC(ie,itb1,izb1) + 
     &            td * zt * RSPEC(ie,itb2,izb1) + 
     &            tt * zd * RSPEC(ie,itb1,izb2) + 
     &            td * zd * RSPEC(ie,itb2,izb2)
      ENDIF

      return
      end
