c     Routines for analysis of hydrostatic equilibirum
c
c     Output_HSE_Profiles
c     Output_HSE_Profiles1
c     Read_Halo_Catalog
c     Mark_Halo_Particles
c     Mark_Halo_Cells
c     Output_HSE_Profiles2
c     Output_HSE_Profiles3
c     Output_HSE_Vir
c     
c     --------------------------------------------------------
      subroutine Output_HSE_Profiles ( xc, yc, zc, rmin, rmax, 
     &                        vcut, dvir, dcrit, nrbin, 
     &                        l_force_center, fname, hname )
c     --------------------------------------------------------
c
c     1. compute the profiles of hydrostatic mass and mach number
c     2. compute the velocity histogram 
c
c     input:  
c
c      integer ifindcenter (see below) 
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fname - full path and file name to output profile to 
c                               it must end with a space
c                   hname - full path and file name of the halo catalog
c     
      include 'a_hfind.h'
c
      real xc, yc, zc, rmin, rmax
      real vcut, dvir, dcrit
      integer nrbin
      logical l_force_center 
      character*256 fname, hname
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      real rvir, rcrit
      integer ifindcenter
      integer m, mbin, Level, icell, nfpn
      double precision amHSEfact, ptfact, xx, yy, zz
      character*256 cldir, clname
c
      integer cnthalo(30000)
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints),pzIa(npoints), pzII(npoints)
      real cs(npoints), vxp(npoints), vyp(npoints), vzp(npoints)
      integer icmark(npoints)
      common / POINTS01 / xpn, ypn, zpn, pdm, pd, pp, 
     +     ptemp, pent, pzIa, pzII, cs, vxp, vyp, vzp, icmark
      parameter ( nmaxbins = 1000 )
      integer nbin(0:nmaxbins),n2bin(0:nmaxbins)
      real rbin(0:nmaxbins),dtibin(0:nmaxbins)
      real dgbin(0:nmaxbins),dstbin(0:nmaxbins),ddmbin(0:nmaxbins)
      real dgibin(0:nmaxbins),dstibin(0:nmaxbins),ddmibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),entbin(0:nmaxbins)
      real zIabin(0:nmaxbins),zIIbin(0:nmaxbins),mgbin(0:nmaxbins)
      real amHSEbin1(0:nmaxbins),amHSEbin2(0:nmaxbins),
     +     amHSEbin3(0:nmaxbins)
      real ppvx(0:nmaxbins),ppvy(0:nmaxbins),ppvz(0:nmaxbins)
      real pvx2(0:nmaxbins),pvy2(0:nmaxbins),pvz2(0:nmaxbins)
      real ptbin(0:nmaxbins),vrms1(0:nmaxbins),csbin(0:nmaxbins)
      real dg2bin(0:nmaxbins),p2bin(0:nmaxbins)
      real qvx(0:nmaxbins),qvy(0:nmaxbins),qvz(0:nmaxbins)
      real qvx2(0:nmaxbins),qvy2(0:nmaxbins),qvz2(0:nmaxbins)
      real pt2bin(0:nmaxbins),vrms2(0:nmaxbins),cs2bin(0:nmaxbins) 
      real pSGbin(0:nmaxbins),ptSGbin(0:nmaxbins),pt2SGbin(0:nmaxbins)
      real ndmp(0:nmaxbins),nstp(0:nmaxbins)


c     arrays for interpolations (polynomicals/cubic spline)
      real rsp(nmaxbins),psp(nmaxbins),psp2nd(nmaxbins)
      real ptsp(nmaxbins),ptsp2nd(nmaxbins)
      real pt2sp(nmaxbins),pt2sp2nd(nmaxbins)
      real dgsp(nmaxbins),dgsp2nd(nmaxbins)
      real rpol(nmaxbins),ppol(nmaxbins),dgpol(nmaxbins),cof(nmaxbins)
      real pSG(nmaxbins),ptSG(nmaxbins),pt2SG(nmaxbins),cofSG(nmaxbins)

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 
      nbins = int((rlmax-rlmin)/drl) + 1
      np1 = nsp(1,2)   ! only the first specie 

      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a5i = 1.0 / a5
      a3i = 1.0 / a3
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      vfact  = v0 / aexpn  ! [km/s]
      ptfact = rho0 * a3i * sun_mass / (mpc*1.0e5)**3 *
     &     (vfact*1.0e5)**2  ! turbulent pressure [erg cm^-3]
      amHSEfact = -1.0 * r0 * P0 / grav_c / rho0 / aexpn *
     &     ((mpc*1.0e5)**2/sun_mass)**2  ! [1/h Msun]

c
c.... calculate basic properties of the main halo 
c
      ifindcenter = 3
      if (  .not. l_force_center ) then 
         write(*,*) ' Find_Halo_Center...'
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      write(*,*) ' Find_Halo_Radius...'
      call Find_Halo_Radius ( xc, yc, zc, dvir, dcrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )

c      write(*,*) xc, yc, zc
c      write(*,*) rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st

    
c
c...  read halo catalog
c
      write(*,*) ' Read_Halo_Catalog...'
      call Read_Halo_Catalog( hname )
        
c
c...  mark all cells associated with subclumps in the main cluster
c    
      write(*,*) 'Mark cells in subhalos...'
      ammin = 0.0     ! [1/h Msun]
      ammax = 1.0e15  ! [1/h Msun]
      do ih = 1, 30000
         cnthalo(ih) = nil
      enddo 
c      call Mark_Halo_Cells( ammin, ammax )
      do ic = 1, nctot
         cnthalo(ind(ic)) = cnthalo(ind(ic)) + 1
      enddo
      do ih = 1, nhalo
         write(*,*) 'halo : ',ih, cnthalo(ih)
      enddo

c
c...  compute the center-of-mass, gas mass profiles and 
c     gas mass-weighted pecular velocity of the cluster 
c
      call Get_MaxLevelNow ()

      write(*,*) ' Computing the center-of-mass position/velocities...'
      do i = 0 , nmaxbins
        mgbin(i)  = 0.
      enddo

      call Compute_Velocities ( xc, yc, zc, vxc, vyc, vzc, rvir, 1 )
      call Compute_Velocities ( xc, yc, zc, vxc, vyc, vzc, rvir, 2 )
      call Compute_Velocities ( xc, yc, zc, vxc, vyc, vzc, rvir, 3 )
    
c
c.... select random samples in the volume 
c
      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin

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
c...  mark all cells associated with subclumps in the main cluster
c 
c      write(*,*) 'Mark cells in subhalos...'
c      ammin = 0.0     ! [1/h Msun]
c      ammax = 1.0e15  ! [1/h Msun]
c      call Mark_Halo_Cells2( ammin, ammax )
c      icnt1 = 0
c      icnt2 = 0
c      do i = 1 , npoints
c         if ( icmark(i) .eq. nil ) then
c            icnt1 = icnt1 + 1
c         else
c            icnt2 = icnt2 + 1
c         endif
c      enddo
c      write(*,*) 'icnts : ',icnt1, icnt2, icnt2/(icnt1+icnt2)

c
c...  compute observables @ the location of each test particle
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
c...    mark cell in sub-halos
c
        icmark(i) = ind(idcell)

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
c         gas momentum
          vxp(i) = hvar(3,idcell)
          vyp(i) = hvar(4,idcell)
          vzp(i) = hvar(5,idcell)
c         sound speed [km/s]
          cs(i) = gamma * boltz * ptemp(i) / wmu / p_mass
          cs(i) = sqrt(cs(i)) / (1.0e+5)  ! km/s
c          vtot(i) = (vxp-vxhalo)**2 + (vyp-vyhalo)**2 + (vzp-vzhalo)**2
c          vtot(i) = vfact * sqrt(vtot(i)) ! km/s
c          ppt(i)  = ptfact * pd(i) * vtot(i)**2
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo

c
c.... construct profiles
c
      do i = 0 , nmaxbins
        nbin(i)  = 0
        ndmp(i)  = 0

        rbin(i)    = 0.
        dgbin(i)   = 0.
        dstbin(i)  = 0.
        ddmbin(i)  = 0.
        dgibin(i)  = 0.
        dstibin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i)  = 0.
        amHSEbin1(i) = 0.
        amHSEbin2(i) = 0.
        amHSEbin3(i) = 0.

        ppvx(i) = 0.
        ppvy(i) = 0.
        ppvz(i) = 0.
        pvx2(i) = 0.
        pvy2(i) = 0.
        pvz2(i) = 0.

        pbin(i)    = 0.
        ptbin(i)   = 0.
        tbin(i)    = 0.
        entbin(i)  = 0.
        vrms1(i)   = 0.
        csbin(i)   = 0.
        zIabin(i)  = 0.
        zIIbin(i)  = 0.

        qvx(i) = 0.
        qvy(i) = 0.
        qvz(i) = 0.
        qvx2(i) = 0.
        qvy2(i) = 0.
        qvz2(i) = 0.

        p2bin(i)  = 0.
        pt2bin(i) = 0.
        vrms2(i)  = 0.
        cs2bin(i) = 0.

        pSGbin(i)  = 0.
        ptSGbin(i) = 0.
        pt2SGbin(i) = 0.
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

c
c...  compute gas profiles
c
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp 
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin)   = pbin(ibin)  + pp(ii)
          tbin(ibin)   = tbin(ibin)  + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          zIabin(ibin) = zIabin(ibin) + pzIa(ii) * pd(ii)
          zIIbin(ibin) = zIIbin(ibin) + pzII(ii) * pd(ii)
          ppvx(ibin) = ppvx(ibin) + vxp(ii)
          ppvy(ibin) = ppvy(ibin) + vyp(ii)
          ppvz(ibin) = ppvz(ibin) + vzp(ii)
          pvx2(ibin) = pvx2(ibin) + vxp(ii)**2 / pd(ii)
          pvy2(ibin) = pvy2(ibin) + vyp(ii)**2 / pd(ii)
          pvz2(ibin) = pvz2(ibin) + vzp(ii)**2 / pd(ii)
          csbin(ibin) = csbin(ibin) + cs(ii) * pd(ii)

          if ( icmark(ii) .eq. nil ) then
            n2bin(ibin)  = n2bin(ibin) + 1
            dg2bin(ibin) = dg2bin(ibin) + pd(ii)
            p2bin(ibin)  = p2bin(ibin) + pp(ii)
            qvx(ibin)  = qvx(ibin) + vxp(ii)
            qvy(ibin)  = qvy(ibin) + vyp(ii)
            qvz(ibin)  = qvz(ibin) + vzp(ii)
            qvx2(ibin) = qvx2(ibin) + vxp(ii)**2 / pd(ii)
            qvy2(ibin) = qvy2(ibin) + vyp(ii)**2 / pd(ii)
            qvz2(ibin) = qvz2(ibin) + vzp(ii)**2 / pd(ii)
            cs2bin(ibin) = cs2bin(ibin) + cs(ii) * pd(ii)
          endif
          
        endif
      enddo

c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          ! compute average gas profiles
          rnp = 1.0/nbin(i)
          rdp = 1.0/dgbin(i)
          rbin(i)   = rbin(i) * rnp
          dgbin(i)  = dgbin(i) * rnp
          tbin(i)   = tbin(i) * rnp
          pbin(i)   = pbin(i) * rnp
          entbin(i) = entbin(i) * rnp
          zIabin(i) = zIabin(i) * rdp
          zIIbin(i) = zIIbin(i) * rdp
          rhov   = ( ppvx(i)**2 + ppvy(i)**2 + ppvz(i)**2 ) * rnp * rnp
          rhov2  = ( pvx2(i) + pvy2(i) + pvz2(i) ) * rnp
          rhos2  = rhov2 - rhov / dgbin(i) 
          ptbin(i) = ptfact * rhos2
          vrms1(i) = vfact * sqrt ( rhos2 / dgbin(i) )
          csbin(i) = csbin(i) * rdp
	endif

        if ( n2bin(i) .gt. 0 ) then
           rnp = 1.0/n2bin(i)
           rdp = 1.0/dg2bin(i)
           p2bin(i)  = p2bin(i) * rnp
           dg2bin(i) = dg2bin(i) * rnp
           rhov   = ( qvx(i)**2 + qvy(i)**2 + qvz(i)**2 ) * rnp * rnp
           rhov2  = ( qvx2(i) + qvy2(i) + qvz2(i) ) * rnp
           rhos2  = rhov2 - rhov / dg2bin(i) 
           pt2bin(i) = ptfact * rhos2
           vrms2(i) = vfact * sqrt ( rhos2 / dg2bin(i) )
           cs2bin(i) = cs2bin(i) * rdp
        endif
      enddo

c
c...  spline interpolate pressure & gas density
c
      do i = 0 , nbins
        rsp(i+1)  = rbin(i)
        psp(i+1)  = pbin(i) / P0 / a5i
        ptsp(i+1) = p2bin(i) / P0 / a5i
        pt2sp(i+1) = (p2bin(i)+pt2bin(i)) / P0 / a5i
        dgsp(i+1) = dgbin(i)
c        write(*,*) pbin(i),ptbin(i),ptbin(i)/pbin(i)
      enddo

      ! interpolate thermal pressure (with clumps)
      yp1 = 0.0
      ypn = 0.0
      call spline(rsp,psp,nbins+1,yp1,ypn,psp2nd)

      ! interpolate thermal pressure (no clump)
      call spline(rsp,ptsp,nbins+1,yp1,ypn,ptsp2nd)

      ! interpolate thermal + turbulent pressure (no clump)
      call spline(rsp,pt2sp,nbins+1,yp1,ypn,pt2sp2nd)

      ! interpolate gas density
      call spline(rsp,dgsp,nbins+1,yp1,ypn,dgsp2nd)

c
c...  set up arrays for the S.G. smoothing filter
c     Array needs to be uniformly space in x
c
      do i = 0 , nbins
        rr = 10.**(rlmin + float(i)*drl)
        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure (with clumps)
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal pressure (no clump)
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent (no clump) 
        pSG(i+1)  = log10(pr)
        ptSG(i+1) = log10(ptr)
        pt2SG(i+1) = log10(pt2r)
c        write(*,*) rr, pSG(i+1),ptSG(i+1)
      enddo
      
c
c...  compute HSE mass profiles 
c
      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      ddmi = 0.0
      dsti = 0.0
      dgi  = 0.0
      dti  = 0.0

      dcritave = dcrit * aexpn**3 / Om0
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
          ! compute virial radius
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. dvir .and.
     &         dtibin(i-1) .ge. dvir ) then 
            ivir = i
            rvir = (dvir * (rr - rl) + rl*dtibin(i) - 
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
        
        ! compute the HSE mass profiles
        rr = 10.**(rlmin + float(i)*drl)
        if ( rr .le. rbin(0) ) then
           rr = rbin(0)
        else if ( rr .ge. rbin(nbins) ) then
           rr = rbin(nbins)
        endif

        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal + turbulent + coherent
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent 
        call splint(rsp,dgsp,dgsp2nd,nbins+1,rr,dgr) ! gas density

        ! 1 : derivative using the Savitzky-Golay smoothing filter
        !  m=2 smooth (similar to f.d. with 4 pts on each side)
        !  m=4 ratty 

        ! thermal pressure (with clumps)
        mbin = 4
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pSGbin(i)  = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr )
        pSGbin(i)  = 10.**pSGbin(i) * P0 * a5i
        val1 = pSGbin(i) / P0 / a5i
        write(*,*) i, rr, val1, pr, (val1-pr)/pr
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr1 = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr1 = pr/rr * dpdr1 / drl

        ! thermal pressure (no clump)
        mbin = 4
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        ptSGbin(i) = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr )
        ptSGbin(i) = 10.**ptSGbin(i) * P0 * a5i
        val2 = ptSGbin(i) / P0 / a5i
        write(*,*) i, rr, val2, ptr, (val2-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr2 = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr2 = ptr/rr * dpdr2 / drl

        ! thermal + turbulent pressure (no clump)
        mbin = 4
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pt2SGbin(i) = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr )
        pt2SGbin(i) = 10.**pt2SGbin(i) * P0 * a5i
        val3 = pt2SGbin(i) / P0 / a5i
        write(*,*) i, rr, val3, ptr, (val3-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr3 = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr) 
        dpdr3 = pt2r/rr * dpdr3 / drl

        if ( dpdr1 .ge. 0.0 ) then
           dpdr1 = -1.0e-30
        endif
        if ( dpdr2 .ge. 0.0 ) then
           dpdr2 = -1.0e-30
        endif
        if ( dpdr3 .ge. 0.0 ) then
           dpdr3 = -1.0e-30
        endif
        if ( dgr .ge. 1.0e-30 ) then
           amHSEbin1(i) = amHSEfact * rr**2 / dgr * dpdr1
           amHSEbin2(i) = amHSEfact * rr**2 / dgr * dpdr2
           amHSEbin3(i) = amHSEfact * rr**2 / dgr * dpdr3
        else
           amHSEbin1(i) = 1.0e-30
           amHSEbin2(i) = 1.0e-30
           amHSEbin3(i) = 1.0e-30
        endif
                  
      enddo

c
c...  output profiles
c
      nfpn = index ( fname , ' ' ) - 1
      open ( 12 , file = fname(1:nfpn) )
c
      aMvir = dvir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      aMcrit = dcritave * 
     &        4.18879 * (rcrit*r0)**3 * rho0 / hubble**2 
      rr = 10.**(rlmin + float(ivir)*drl)
      aM_dm  = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
      aM_gas = dgibin(ivir)  * (4.18879 * rr**3) * aM0 * hubble
      aM_st  = dstibin(ivir) * (4.18879 * rr**3) * aM0 * hubble 
      aMvir2 = aM_dm + aM_gas + aM_st 

      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin  
      write(12,12) ivir, rvir*r0, aMvir, dvir
      write(12,22) icrit, rcrit*r0, aMcrit, dcrit
      write(12,17) aM_dm, aM_gas, aM_st, aMvir2
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,' nrbin =',i5)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 22   format ('# icrit = ',i4,' Rcrit =',f7.4,
     &        ' [h^-1 Mpc], Mcrit =',g11.4,
     &  ' [h^-1 Msun] for Dcrit(overd w.r.t. rho_crit)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r  delta_g delta_st delta_dm  deltai_g',
     &        ' deltai_st deltai_dm  massHSE1 massHSE2 massHSE3 ',
     &        ' Tg  Pg  Sg Pturb  vrms  cs  Pturb2  vrms2  cs2 ',
     &        ' PgSG  PtSG  Pt2SG  N_dm_p' )
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
     &         rbin(i)*r0, rl*r0, rmid, rr*r0, 
     &         dgbin(i),dstbin(i),ddmbin(i),dgibin(i),dstibin(i),
     &         ddmibin(i),amHSEbin1(i),amHSEbin2(i),amHSEbin3(i),
     &         tbin(i),entbin(i),
     &         pbin(i),ptbin(i),vrms1(i),csbin(i),
     &         p2bin(i),pt2bin(i),vrms2(i),cs2bin(i),
     &         pSGbin(i),ptSGbin(i),pt2SGbin(i),int(ndmp(i))
	endif
      enddo
 18   format(26(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'


      write(*,*) ' Output_Velocity_Histogram'
      call Output_Velocity_Histogram ( )

      return
      end
c
c     ----------------------------------------------------------
      subroutine Compute_Velocities ( xc, yc, zc, vxc, vyc, vzc, 
     &     rvir, imode )
c     ----------------------------------------------------------
c     
      include 'a_hfind.h'
c
      real xc, yc, zc, vxc, vyc, vzc, rvir
      integer imode 
c
      parameter ( nmaxbins = 1000 )
      real mgbin(0:nmaxbins),vrad(0:nmaxbins),vrot(0:nmaxbins)
      real mJ(0:nmaxbins)
      real mJx(0:nmaxbins),mJy(0:nmaxbins),mJz(0:nmaxbins)
      real ex(0:nmaxbins),ey(0:nmaxbins),ez(0:nmaxbins)
      real vr, Jx, Jy, Jz, Jtot, rJx, rJy, rJz, norm

      call Get_MaxLevelNow ()

      if ( imode .eq. 0 ) then
        xcm = 0.0
        ycm = 0.0
        zcm = 0.0
        vxcm = 0.0
        vycm = 0.0
        vzcm = 0.0
        amcm = 0.0        
        do i = 0 , nmaxbins
            mgbin(i)  = 0.
        enddo
      endif


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
                  xx = Posx 
                  yy = Posy 
                  zz = Posz 
                  vgx = hvar(3,ic) / hvar(1,ic) 
                  vgy = hvar(4,ic) / hvar(1,ic) 
                  vgz = hvar(5,ic) / hvar(1,ic) 
c                  ! use the center of mass as a halo center 
c                  if ( imode .eq. 1 ) then
c                     call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
c                     rr = sqrt( (xx-xc)**2 + (yy-yc)**2 + (zz-zc)**2 )
c                     if ( rr .le. rvir ) then
c                        xcm = xcm + xx * cellmass
c                        ycm = ycm + yy * cellmass
c                        zcm = zcm + zz * cellmass
c                        vxcm = vxcm + vgx * cellmass 
c                        vycm = vycm + vgy * cellmass
c                        vzcm = vzcm + vgz * cellmass
c                        amcm = amcm + cellmass
c                     endif
c                  else
c                     xx = xx - xcm
c                     yy = yy - ycm
c                     zz = zz - zcm
c                     vgx = vgx - vxcm 
c                     vgy = vgy - vycm
c                     vgz = vgz - vzcm
c                  endif

                  ! use the densest cell as a halo center
                  if ( imode .eq. 1 ) then
                     call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )  
                     rr = sqrt( (xx-xc)**2 + (yy-yc)**2 + (zz-zc)**2 )
                     if ( rr .le. rvir ) then
                       xx = xx - xc
                       yy = yy - yc
                       zz = zz - zc
                       vgx = vgx - vxc 
                       vgy = vgy - vyc
                       vgz = vgz - vzc
                       ! radius from the C.O.M.
                       rr = sqrt( xx**2 + yy**2 + zz**2 )
                       vr = ( xx*vgx + yy*vgy + zz*vgz ) / rr
                       ibin = max ( int((log10(rr)-rlmin)/drl) + 1, 0 )
                       ! total angular momentum of each cell
                       Jx = yy*vgz - zz*vgy
                       Jy = zz*vgx - xx*vgz
                       Jz = xx*vgy - yy*vgx
                       cellmass = hvar(1,ic) * CellVolume
                     endif
                  endif

                  ! compute angular momentum of each shell
                  if ( imode .eq. 2 .and. ibin .le. nmaxbins ) then
                     ! gas mass profiles
                     mgbin(ibin) = mgbin(ibin) + cellmass
                     ! average radial velocity of the shell
                     vrad(ibin) = vrad(ibin) + vr * cellmass  
                     ! average angular momentum of each shell
                     mJx(ibin) = mJx(ibin) + Jx * cellmass
                     mJy(ibin) = mJy(ibin) + Jy * cellmass
                     mJz(ibin) = mJz(ibin) + Jz * cellmass
                  endif
  
                  if ( imode .eq. 3 .and. ibin .le. nmaxbins ) then
                     ! angular momentum of each cell
                     Jtot = ex(ibin)*Jx + ey(ibin)*Jy + ez(ibin)*Jz
                     ! projections of the total J on the average J
                     Jx = ex(ibin)*Jtot
                     Jy = ey(ibin)*Jtot
                     Jz = ez(ibin)*Jtot   
                     ! r x J = r x ( r x v ) 
                     rJx = yy*Jz - zz*Jy
                     rJy = zz*Jx - xx*Jz
                     rJz = xx*Jy - yy*Jx
                     norm = -1.0*(Jx**2+Jy**2+Jz**2)/
     &                    (rJx**2+rJy**2+rJz**2)
                     ! rotation velocity
                     vrotx = rJx*norm
                     vroty = rJy*norm
                     vrotz = rJz*norm
                     ! subtruct rotation
                     vgx = vgx - vrotx
                     vgy = vgy - vroty
                     vgz = vgz - vrotz
                     ! radial velocity
                     vradx = xx/rr * vrad(ibin)
                     vrady = yy/rr * vrad(ibin)
                     vradz = zz/rr * vrad(ibin)
                     ! subtract net radial motions
                     vgx = vgx - vradx
                     vgy = vgy - vrady
                     vgz = vgz - vradz

c                     v2 = vgx**2 + vgy**2 + vgz**2
c                     rp = sqrt(-1./norm)
                     vrot(ibin) = sqrt(vrotx**2+vroty**2+vrotz**2)
                  endif

                endif
              enddo
            enddo
         ENDIF
      ENDDO

      do ii = 1, nmaxbins   
         if ( imode .eq. 1 ) then
           xcm = xcm / amcm
           ycm = ycm / amcm
           zcm = zcm / amcm
           vxcm = vxcm / amcm
           vycm = vycm / amcm
           vzcm = vzcm / amcm
         endif
         if ( imode .eq. 2 ) then
           ! compute the radial velocity
           vrad(ii) = vrad(ii) / mgbin(ii) 
           ! total angular momentum of the shell
           mJ(ii) = sqrt( mJx(ii)**2 + mJy(ii)**2 + mJz(ii)**2 )
           ex(ii) = mJx(ii) / mJ(ii)
           ey(ii) = mJy(ii) / mJ(ii)
           ez(ii) = mJz(ii) / mJ(ii)
         endif
         if ( imode .eq. 3 ) then
         endif
      enddo

      return
      end
c
c     ---------------------------------------------------------
      subroutine Output_HSE_Profiles1 ( xc, yc, zc, rmin, rmax, 
     &                        vcut, dvir, dcrit, nrbin, 
     &                        l_force_center, fname, hname )
c     ---------------------------------------------------------
c
c     1. compute the profiles of hydrostatic mass and mach number
c     2. compute the velocity histogram 
c
c     input:  
c
c      integer ifindcenter (see below) 
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fname - full path and file name to output profile to 
c                               it must end with a space
c                   hname - full path and file name of the halo catalog
c     
      include 'a_hfind.h'
c
      real xc, yc, zc, rmin, rmax
      real vcut, dvir, dcrit
      integer nrbin
      logical l_force_center 
      character*256 fname, hname
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      real rvir, rcrit
      integer ifindcenter
      integer m, mbin, Level, icell, nfpn
      double precision amHSEfact, ptfact, xx, yy, zz
      character*256 cldir, clname
c
      integer cnthalo(30000)
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints),pzIa(npoints), pzII(npoints)
      real cs(npoints), vxp(npoints), vyp(npoints), vzp(npoints)
      integer icmark(npoints)
      common / POINTS01 / xpn, ypn, zpn, pdm, pd, pp, 
     +     ptemp, pent, pzIa, pzII, cs, vxp, vyp, vzp, icmark
      parameter ( nmaxbins = 1000 )
      integer nbin(0:nmaxbins),n2bin(0:nmaxbins)
      real rbin(0:nmaxbins),dtibin(0:nmaxbins)
      real dgbin(0:nmaxbins),dstbin(0:nmaxbins),ddmbin(0:nmaxbins)
      real dgibin(0:nmaxbins),dstibin(0:nmaxbins),ddmibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),entbin(0:nmaxbins)
      real zIabin(0:nmaxbins),zIIbin(0:nmaxbins),mgbin(0:nmaxbins)
      real amHSEbin1(0:nmaxbins),amHSEbin2(0:nmaxbins),
     +     amHSEbin3(0:nmaxbins)
      real ppvx(0:nmaxbins),ppvy(0:nmaxbins),ppvz(0:nmaxbins)
      real pvx2(0:nmaxbins),pvy2(0:nmaxbins),pvz2(0:nmaxbins)
      real ptbin(0:nmaxbins),vrms1(0:nmaxbins),csbin(0:nmaxbins)
      real dg2bin(0:nmaxbins),p2bin(0:nmaxbins)
      real qvx(0:nmaxbins),qvy(0:nmaxbins),qvz(0:nmaxbins)
      real qvx2(0:nmaxbins),qvy2(0:nmaxbins),qvz2(0:nmaxbins)
      real pt2bin(0:nmaxbins),vrms2(0:nmaxbins),cs2bin(0:nmaxbins) 
      real pSGbin(0:nmaxbins),ptSGbin(0:nmaxbins),pt2SGbin(0:nmaxbins)
      real ndmp(0:nmaxbins),nstp(0:nmaxbins)
      
c     arrays for interpolations (polynomicals/cubic spline)
      real rsp(nmaxbins),psp(nmaxbins),psp2nd(nmaxbins)
      real ptsp(nmaxbins),ptsp2nd(nmaxbins)
      real pt2sp(nmaxbins),pt2sp2nd(nmaxbins)
      real dgsp(nmaxbins),dgsp2nd(nmaxbins)
      real rpol(nmaxbins),ppol(nmaxbins),dgpol(nmaxbins),cof(nmaxbins)
      real pSG(nmaxbins),ptSG(nmaxbins),pt2SG(nmaxbins),cofSG(nmaxbins)

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 
      nbins = int((rlmax-rlmin)/drl) + 1
      np1 = nsp(1,2)   ! only the first specie 

      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a5i = 1.0 / a5
      a3i = 1.0 / a3
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      vfact  = v0 / aexpn  ! [km/s]
      ptfact = rho0 * a3i * sun_mass / (mpc*1.0e5)**3 *
     &     (vfact*1.0e5)**2  ! turbulent pressure [erg cm^-3]
      amHSEfact = -1.0 * r0 * P0 / grav_c / rho0 / aexpn *
     &     ((mpc*1.0e5)**2/sun_mass)**2  ! [1/h Msun]

c
c.... calculate basic properties of the main halo 
c
      ifindcenter = 3
      if (  .not. l_force_center ) then 
         write(*,*) ' Find_Halo_Center...'
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      write(*,*) ' Find_Halo_Radius...'
      call Find_Halo_Radius ( xc, yc, zc, dvir, dcrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )

c      write(*,*) xc, yc, zc
c      write(*,*) rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st

c
c...  compute the gas mass profiles and 
c     gas mass-weighted pecular velocity of the cluster 
c
      call Get_MaxLevelNow ()

      write(*,*) ' Computing gas mass profiles etc...'
      amgw = 0.0
      vxhalo = 0.0
      vyhalo = 0.0
      vzhalo = 0.0
      do i = 0 , nmaxbins
        mgbin(i)  = 0.
      enddo
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
                  if ( ibin .le. nmaxbins ) then
                     mgbin(ibin) = mgbin(ibin) + cellmass
c                     write(*,*) 'mg :',ic, ibin, mgbin(ibin), cellmass
                  endif
                  if ( rp .le. rvir ) then
                     amgw = amgw + cellmass
                     vxhalo = vxhalo + hvar(3,ic) * CellVolume
                     vyhalo = vyhalo + hvar(4,ic) * CellVolume
                     vzhalo = vzhalo + hvar(5,ic) * CellVolume
c                     write(*,*) 'v :',rp, rvir, amgw, vxhalo, vyhalo, 
c     &                    vzhalo
                  endif
                endif
              enddo
            enddo
         ENDIF
      ENDDO

      ! the gas mass-weighted cluster peculiar velocity [code units]
      vxhalo = vxhalo / amgw
      vyhalo = vyhalo / amgw
      vzhalo = vzhalo / amgw
      write(*,*) 'vxh, vyh, vzh = ',vxhalo, vyhalo, vzhalo
      
c
c...  read halo catalog
c
      write(*,*) ' Read_Halo_Catalog...'
      call Read_Halo_Catalog( hname )
        
c
c...  mark all cells associated with subclumps in the main cluster
c    
      write(*,*) 'Mark cells in subhalos...'
      ammin = 0.0     ! [1/h Msun]
      ammax = 1.0e15  ! [1/h Msun]
      do ih = 1, 30000
         cnthalo(ih) = nil
      enddo 
      call Mark_Halo_Cells( ammin, ammax )
      do ic = 1, nctot
         cnthalo(ind(ic)) = cnthalo(ind(ic)) + 1
      enddo
      do ih = 1, nhalo
         write(*,*) 'halo : ',ih, cnthalo(ih)
      enddo
      
c
c.... select random samples in the volume 
c
      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin

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
c...  mark all cells associated with subclumps in the main cluster
c 
c      write(*,*) 'Mark cells in subhalos...'
c      ammin = 0.0     ! [1/h Msun]
c      ammax = 1.0e15  ! [1/h Msun]
c      call Mark_Halo_Cells2( ammin, ammax )
c      icnt1 = 0
c      icnt2 = 0
c      do i = 1 , npoints
c         if ( icmark(i) .eq. nil ) then
c            icnt1 = icnt1 + 1
c         else
c            icnt2 = icnt2 + 1
c         endif
c      enddo
c      write(*,*) 'icnts : ',icnt1, icnt2, icnt2/(icnt1+icnt2)

c
c...  compute observables @ the location of each test particle
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
c...    mark cell in sub-halos
c
        icmark(i) = ind(idcell)

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
c         gas momentum
          vxp(i) = hvar(3,idcell)
          vyp(i) = hvar(4,idcell)
          vzp(i) = hvar(5,idcell)
c         sound speed [km/s]
          cs(i) = gamma * boltz * ptemp(i) / wmu / p_mass
          cs(i) = sqrt(cs(i)) / (1.0e+5)  ! km/s
c          vtot(i) = (vxp-vxhalo)**2 + (vyp-vyhalo)**2 + (vzp-vzhalo)**2
c          vtot(i) = vfact * sqrt(vtot(i)) ! km/s
c          ppt(i)  = ptfact * pd(i) * vtot(i)**2
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo

c
c.... construct profiles
c
      do i = 0 , nmaxbins
        nbin(i)  = 0
        ndmp(i)  = 0

        rbin(i)    = 0.
        dgbin(i)   = 0.
        dstbin(i)  = 0.
        ddmbin(i)  = 0.
        dgibin(i)  = 0.
        dstibin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i)  = 0.
        amHSEbin1(i) = 0.
        amHSEbin2(i) = 0.
        amHSEbin3(i) = 0.

        ppvx(i) = 0.
        ppvy(i) = 0.
        ppvz(i) = 0.
        pvx2(i) = 0.
        pvy2(i) = 0.
        pvz2(i) = 0.

        pbin(i)    = 0.
        ptbin(i)   = 0.
        tbin(i)    = 0.
        entbin(i)  = 0.
        vrms1(i)   = 0.
        csbin(i)   = 0.
        zIabin(i)  = 0.
        zIIbin(i)  = 0.

        qvx(i) = 0.
        qvy(i) = 0.
        qvz(i) = 0.
        qvx2(i) = 0.
        qvy2(i) = 0.
        qvz2(i) = 0.

        p2bin(i)  = 0.
        pt2bin(i) = 0.
        vrms2(i)  = 0.
        cs2bin(i) = 0.

        pSGbin(i)  = 0.
        ptSGbin(i) = 0.
        pt2SGbin(i) = 0.
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

c
c...  compute gas profiles
c
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp 
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin)   = pbin(ibin)  + pp(ii)
          tbin(ibin)   = tbin(ibin)  + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          zIabin(ibin) = zIabin(ibin) + pzIa(ii) * pd(ii)
          zIIbin(ibin) = zIIbin(ibin) + pzII(ii) * pd(ii)
          ppvx(ibin) = ppvx(ibin) + vxp(ii)
          ppvy(ibin) = ppvy(ibin) + vyp(ii)
          ppvz(ibin) = ppvz(ibin) + vzp(ii)
          pvx2(ibin) = pvx2(ibin) + vxp(ii)**2 / pd(ii)
          pvy2(ibin) = pvy2(ibin) + vyp(ii)**2 / pd(ii)
          pvz2(ibin) = pvz2(ibin) + vzp(ii)**2 / pd(ii)
          csbin(ibin) = csbin(ibin) + cs(ii) * pd(ii)

          if ( icmark(ii) .eq. nil ) then
            n2bin(ibin)  = n2bin(ibin) + 1
            dg2bin(ibin) = dg2bin(ibin) + pd(ii)
            p2bin(ibin)  = p2bin(ibin) + pp(ii)
            qvx(ibin)  = qvx(ibin) + vxp(ii)
            qvy(ibin)  = qvy(ibin) + vyp(ii)
            qvz(ibin)  = qvz(ibin) + vzp(ii)
            qvx2(ibin) = qvx2(ibin) + vxp(ii)**2 / pd(ii)
            qvy2(ibin) = qvy2(ibin) + vyp(ii)**2 / pd(ii)
            qvz2(ibin) = qvz2(ibin) + vzp(ii)**2 / pd(ii)
            cs2bin(ibin) = cs2bin(ibin) + cs(ii) * pd(ii)
          endif
          
        endif
      enddo

c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          ! compute average gas profiles
          rnp = 1.0/nbin(i)
          rdp = 1.0/dgbin(i)
          rbin(i)   = rbin(i) * rnp
          dgbin(i)  = dgbin(i) * rnp
          tbin(i)   = tbin(i) * rnp
          pbin(i)   = pbin(i) * rnp
          entbin(i) = entbin(i) * rnp
          zIabin(i) = zIabin(i) * rdp
          zIIbin(i) = zIIbin(i) * rdp
          rhov   = ( ppvx(i)**2 + ppvy(i)**2 + ppvz(i)**2 ) * rnp * rnp
          rhov2  = ( pvx2(i) + pvy2(i) + pvz2(i) ) * rnp
          rhos2  = rhov2 - rhov / dgbin(i) 
          ptbin(i) = ptfact * rhos2
          vrms1(i) = vfact * sqrt ( rhos2 / dgbin(i) )
          csbin(i) = csbin(i) * rdp
	endif

        if ( n2bin(i) .gt. 0 ) then
           rnp = 1.0/n2bin(i)
           rdp = 1.0/dg2bin(i)
           p2bin(i)  = p2bin(i) * rnp
           dg2bin(i) = dg2bin(i) * rnp
           rhov   = ( qvx(i)**2 + qvy(i)**2 + qvz(i)**2 ) * rnp * rnp
           rhov2  = ( qvx2(i) + qvy2(i) + qvz2(i) ) * rnp
           rhos2  = rhov2 - rhov / dg2bin(i) 
           pt2bin(i) = ptfact * rhos2
           vrms2(i) = vfact * sqrt ( rhos2 / dg2bin(i) )
           cs2bin(i) = cs2bin(i) * rdp
        endif
      enddo

c
c...  spline interpolate pressure & gas density
c
      do i = 0 , nbins
        rsp(i+1)  = rbin(i)
        psp(i+1)  = pbin(i) / P0 / a5i
        ptsp(i+1) = p2bin(i) / P0 / a5i
        pt2sp(i+1) = (p2bin(i)+pt2bin(i)) / P0 / a5i
        dgsp(i+1) = dgbin(i)
c        write(*,*) pbin(i),ptbin(i),ptbin(i)/pbin(i)
      enddo

      ! interpolate thermal pressure (with clumps)
      yp1 = 0.0
      ypn = 0.0
      call spline(rsp,psp,nbins+1,yp1,ypn,psp2nd)

      ! interpolate thermal pressure (no clump)
      call spline(rsp,ptsp,nbins+1,yp1,ypn,ptsp2nd)

      ! interpolate thermal + turbulent pressure (no clump)
      call spline(rsp,pt2sp,nbins+1,yp1,ypn,pt2sp2nd)

      ! interpolate gas density
      call spline(rsp,dgsp,nbins+1,yp1,ypn,dgsp2nd)

c
c...  set up arrays for the S.G. smoothing filter
c     Array needs to be uniformly space in x
c
      do i = 0 , nbins
        rr = 10.**(rlmin + float(i)*drl)
        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure (with clumps)
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal pressure (no clump)
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent (no clump) 
        pSG(i+1)  = log10(pr)
        ptSG(i+1) = log10(ptr)
        pt2SG(i+1) = log10(pt2r)
c        write(*,*) rr, pSG(i+1),ptSG(i+1)
      enddo
      
c
c...  compute HSE mass profiles 
c
      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      ddmi = 0.0
      dsti = 0.0
      dgi  = 0.0
      dti  = 0.0

      dcritave = dcrit * aexpn**3 / Om0
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
          ! compute virial radius
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. dvir .and.
     &         dtibin(i-1) .ge. dvir ) then 
            ivir = i
            rvir = (dvir * (rr - rl) + rl*dtibin(i) - 
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
        
        ! compute the HSE mass profiles
        rr = 10.**(rlmin + float(i)*drl)
        if ( rr .le. rbin(0) ) then
           rr = rbin(0)
        else if ( rr .ge. rbin(nbins) ) then
           rr = rbin(nbins)
        endif

        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal + turbulent + coherent
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent 
        call splint(rsp,dgsp,dgsp2nd,nbins+1,rr,dgr) ! gas density

        ! 1 : derivative using the Savitzky-Golay smoothing filter
        !  m=2 smooth (similar to f.d. with 4 pts on each side)
        !  m=4 ratty 

        ! thermal pressure (with clumps)
        mbin = 4
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pSGbin(i)  = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr )
        pSGbin(i)  = 10.**pSGbin(i) * P0 * a5i
        val1 = pSGbin(i) / P0 / a5i
        write(*,*) i, rr, val1, pr, (val1-pr)/pr
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr1 = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr1 = pr/rr * dpdr1 / drl

        ! thermal pressure (no clump)
        mbin = 4
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        ptSGbin(i) = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr )
        ptSGbin(i) = 10.**ptSGbin(i) * P0 * a5i
        val2 = ptSGbin(i) / P0 / a5i
        write(*,*) i, rr, val2, ptr, (val2-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr2 = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr2 = ptr/rr * dpdr2 / drl

        ! thermal + turbulent pressure (no clump)
        mbin = 4
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pt2SGbin(i) = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr )
        pt2SGbin(i) = 10.**pt2SGbin(i) * P0 * a5i
        val3 = pt2SGbin(i) / P0 / a5i
        write(*,*) i, rr, val3, ptr, (val3-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr3 = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr) 
        dpdr3 = pt2r/rr * dpdr3 / drl

        if ( dpdr1 .ge. 0.0 ) then
           dpdr1 = -1.0e-30
        endif
        if ( dpdr2 .ge. 0.0 ) then
           dpdr2 = -1.0e-30
        endif
        if ( dpdr3 .ge. 0.0 ) then
           dpdr3 = -1.0e-30
        endif
        if ( dgr .ge. 1.0e-30 ) then
           amHSEbin1(i) = amHSEfact * rr**2 / dgr * dpdr1
           amHSEbin2(i) = amHSEfact * rr**2 / dgr * dpdr2
           amHSEbin3(i) = amHSEfact * rr**2 / dgr * dpdr3
        else
           amHSEbin1(i) = 1.0e-30
           amHSEbin2(i) = 1.0e-30
           amHSEbin3(i) = 1.0e-30
        endif
                  
      enddo

c
c...  output profiles
c
      nfpn = index ( fname , ' ' ) - 1
      open ( 12 , file = fname(1:nfpn) )
c
      aMvir = dvir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      aMcrit = dcritave * 
     &        4.18879 * (rcrit*r0)**3 * rho0 / hubble**2 
      rr = 10.**(rlmin + float(ivir)*drl)
      aM_dm  = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
      aM_gas = dgibin(ivir)  * (4.18879 * rr**3) * aM0 * hubble
      aM_st  = dstibin(ivir) * (4.18879 * rr**3) * aM0 * hubble 
      aMvir2 = aM_dm + aM_gas + aM_st 

      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin  
      write(12,12) ivir, rvir*r0, aMvir, dvir
      write(12,22) icrit, rcrit*r0, aMcrit, dcrit
      write(12,17) aM_dm, aM_gas, aM_st, aMvir2
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,' nrbin =',i5)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 22   format ('# icrit = ',i4,' Rcrit =',f7.4,
     &        ' [h^-1 Mpc], Mcrit =',g11.4,
     &  ' [h^-1 Msun] for Dcrit(overd w.r.t. rho_crit)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r  delta_g delta_st delta_dm  deltai_g',
     &        ' deltai_st deltai_dm  massHSE1 massHSE2 massHSE3 ',
     &        ' Tg  Pg  Sg Pturb  vrms  cs  Pturb2  vrms2  cs2 ',
     &        ' PgSG  PtSG  Pt2SG  N_dm_p' )
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
     &         rbin(i)*r0, rl*r0, rmid, rr*r0, 
     &         dgbin(i),dstbin(i),ddmbin(i),dgibin(i),dstibin(i),
     &         ddmibin(i),amHSEbin1(i),amHSEbin2(i),amHSEbin3(i),
     &         tbin(i),entbin(i),
     &         pbin(i),ptbin(i),vrms1(i),csbin(i),
     &         p2bin(i),pt2bin(i),vrms2(i),cs2bin(i),
     &         pSGbin(i),ptSGbin(i),pt2SGbin(i),int(ndmp(i))
	endif
      enddo
 18   format(26(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'


      write(*,*) ' Output_Velocity_Histogram'
      call Output_Velocity_Histogram ( )

      return
      end
c
c     -------------------------------------
      subroutine Read_Halo_Catalog( hname )
c     -------------------------------------
c
c     purpose : read in halo catalog
c 
      include 'a_hfind.h'

      character*256 hname, textline 
      
      
      nfn = index ( hname , ' ' ) - 1 
      write(*,*) ' Reading : ',hname(1:nfn)
      open ( 20 , file = hname(1:nfn) )
      do ic1 = 1 , 5
        read(20,*) textline
      enddo

      nh = 0 
      rfact = r0 * 1.e3 * aexpn / hubble ! code -> proper kpc
      do ic1 = 1 , nhmax 
        read(20,*,end=10) xhd, yhd, zhd, rvd, rv, amg, amcg, amst, 
     &        amdm, zg, zst, tst, amhalo
        nh = nh + 1
        if ( nh .gt. nhmax ) then 
          write(*,*) '* error in Halo_Properties: nh>nhmax:',
     &    nh, nhmax
          write(*,*) 'increase nhmax and rerun. stopping...'
          stop
        endif
        xh(nh) = xhd/r0 + 1.0 
        yh(nh) = yhd/r0 + 1.0 
        zh(nh) = zhd/r0 + 1.0 
c       For isolated halo : rvir(nh) = rv  & amv(nh) = amv
c       For sub-halos     : rvir(nh) = rvd & am(nh)  = amd
        rh(nh)    = rvd / rfact  ! proper kpc -> code unit
        rhvir(nh) = rv  / rfact  ! proper kpc -> code unit
        amh(nh)   = amg + amst + amdm
c        write(*,*) nh, xh(nh), yh(nh), zh(nh), rh(nh), rhvir(nh),amh(nh)
      enddo
 10   close(20)

      nhalo = nh
      write(*,*) 'Read in ',nhalo,' halos'

      return
      end
c
c     ----------------------------------------------
      subroutine Mark_Halo_Particles( ammin, ammax )
c     ----------------------------------------------
c
c     purpose : mark cells associated with subclumps in the main halo
c 
      include 'a_hfind.h'
      
c
c...  Initialize arrays and construct linked list
c
      write(*,*) 'Initialize arrays for linked list...'
      do ic1 = 1 , nll
        do ic2 = 1 , nll
          do ic3 = 1 , nll 
            iCL2(ic3,ic2,ic1) = nil 
          enddo
        enddo
      enddo

      do ic1 = 1 , np1max 
        iLL2(ic1) = nil 
        iSp(ic1)  = nil 
      enddo

      write(*,*) 'Constructing Linked List...'
      call LL_Construct_HF ()

c
c...  Figure out ID of the most massive halo
c
      imaxhalo = 0
      amaxhalo = -1.0e30
      do ic1 = 1 , nhalo  
c         write(*,*) ic1, amh(ic1)
         if ( amh(ic1) .ge. amaxhalo ) then
            imaxhalo = ic1
            amaxhalo = amh(ic1)
         endif
      enddo

      if ( imaxhalo .eq. 0 ) then
         write(*,*) 'imaxhalo = ',imaxhalo
         write(*,*) ' Did not find the most massive halo...'
      else
         write(*,*) ' ih, mh = ',imaxhalo, amaxhalo
      endif

c
c...  Mark cells in all halos except in the most most massive halo
c
c      ih = 0 
      do ic1 = 1 , nhalo  
        icnt = 0
        if ( amh(ic1) .gt. ammin .and. amh(ic1) .lt. ammax .and. 
     &        ic1 .ne. imaxhalo ) then
          rhalo = 1.05 * rh(ic1)
          xhalo = xh(ic1) 
          yhalo = yh(ic1) 
          zhalo = zh(ic1) 
          write(*,*) rh(ic1), xhalo, yhalo, zhalo

c....     convert from original grid to chaining mesh          
          xd    = xhalo - rhalo
          yd    = yhalo - rhalo
          zd    = zhalo - rhalo
          imin  = int(float(nll) * (xd - 1) / ng + 1)
          jmin  = int(float(nll) * (yd - 1) / ng + 1)
          kmin  = int(float(nll) * (zd - 1) / ng + 1)
          xd    = xhalo + rhalo
          yd    = yhalo + rhalo
          zd    = zhalo + rhalo
          imax  = int(float(nll) * (xd - 1) / ng + 1)
          jmax  = int(float(nll) * (yd - 1) / ng + 1)
          kmax  = int(float(nll) * (zd - 1) / ng + 1)

c....     sweep over neighbors checking periodic boundary conditions
          iph = 0
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idm = iCL2(ic,jc,kc)  ! read LL head
c                write(*,*) ' iCL2 :',ic,jc,kc,idm, iCL2(ic,jc,kc)
                do while ( idm .ne. nil )
                  xd = x(idm)
                  yd = y(idm)
                  zd = z(idm)
                  vxd = vx(idm)
                  vyd = vy(idm)
                  vzd = vz(idm)
c....             enforce periodic boundary conditions
                  diff_x = xd - xhalo                  
                  diff_y = yd - yhalo 
                  diff_z = zd - zhalo 
                  corr_x = zero 
                  corr_y = zero 
                  corr_z = zero 
                  
                  if ( abs(diff_x) .gt. ng/2 ) then 
                    if ( diff_x .gt. zero ) then 
                      corr_x = -ng 
                    else
                      corr_x =  ng
                    endif
                  endif
                  if ( abs(diff_y) .gt. ng/2 ) then 
                    if ( diff_y .gt. zero ) then 
                      corr_y = -ng 
                    else
                      corr_y =  ng
                    endif
                  endif
                  if ( abs(diff_z) .gt. ng/2 ) then 
                    if ( diff_z .gt. zero ) then 
                      corr_z = -ng 
                    else
                      corr_z =  ng
                    endif
                  endif                 
c
c....             compute distance particle - halo center squared
c                 corr_x , corr_y , corr_z - take care of periodicity
                  rr = sqrt( (xd - xhalo + corr_x)**2 + 
     &                       (yd - yhalo + corr_y)**2 + 
     &                       (zd - zhalo + corr_z)**2   )
                  if ( rr .le. rh(ic1) ) then
                     iSp(idm) = 1
                     icnt = icnt + 1
c                     write(*,*) ic1, icnt, xd, yd, zd, rh(ic1), iSp(idm)
                  endif
                  idm = iLL2(idm) ! next particle from linked list
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! end i
       
        endif
      enddo ! end ic1

      return
      end

c
c     ------------------------------------------
      subroutine Mark_Halo_Cells( ammin, ammax )
c     ------------------------------------------
c
c     purpose : mark cells associated with subclumps in the main halo
c 
      include 'a_hfind.h'
      
      do ic = 1, nctot
         ind(ic) = nil
      enddo
c
c...  Figure out ID of the most massive halo
c
      imaxhalo = 0
      amaxhalo = -1.0e30
      do ih = 1 , nhalo  
         if ( amh(ih) .ge. amaxhalo ) then
            imaxhalo = ih
            amaxhalo = amh(ih)
            xch = xh(ih)
            ych = yh(ih)
            zch = zh(ih)
            rch = 1.5 * rh(ih) 
         endif
      enddo

      if ( imaxhalo .eq. 0 ) then
         write(*,*) 'imaxhalo = ',imaxhalo
         write(*,*) ' Did not find the most massive halo...'
      else
         write(*,*) ' ih, mh = ',imaxhalo, amaxhalo
      endif


      DO Level = MinLevel, MaxLevelNow
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
                ! loop around halos
                do ih = 1 , nhalo  
                   if ( amh(ih) .gt. ammin .and. 
     &                  amh(ih) .lt. ammax .and. 
     &                  ih .ne. imaxhalo ) then
                      call PeriodicBoundary( 
     &                     xh(ih), yh(ih), zh(ih), 
     &                     Posx, Posy, Posz )
                      rr = sqrt( (Posx-xh(ih))**2 + 
     &                     (Posy-yh(ih))**2 + 
     &                     (Posz-zh(ih))**2 )
                      if ( rr .le. rh(ih) ) then
                         ind(ic) = ih
                         goto 199
                      endif
                   endif
                enddo ! ih
 199          endif 
            enddo !ic2
          enddo  !ic1

        ENDIF
      ENDDO
      
      return
      end
c
c     -------------------------------------------
      subroutine Mark_Halo_Cells2( ammin, ammax )
c     -------------------------------------------
c
c     purpose : mark cells associated with subclumps in the main halo
c 
      include 'a_hfind.h'

      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints),pzIa(npoints), pzII(npoints)
      real cs(npoints), vxp(npoints), vyp(npoints), vzp(npoints)
      integer icmark(npoints)
      common / POINTS01 / xpn, ypn, zpn, pdm, pd, pp, 
     +     ptemp, pent, pzIa, pzII, cs, vxp, vyp, vzp, icmark
      double precision xx, yy, zz

c
c...  Figure out ID of the most massive halo
c
      imaxhalo = 0
      amaxhalo = -1.0e30
      do ih = 1 , nhalo  
         if ( amh(ih) .ge. amaxhalo ) then
            imaxhalo = ih
            amaxhalo = amh(ih)
            xch = xh(ih)
            ych = yh(ih)
            zch = zh(ih)
            rch = 1.5 * rh(ih) 
         endif
      enddo

      if ( imaxhalo .eq. 0 ) then
         write(*,*) 'imaxhalo = ',imaxhalo
         write(*,*) ' Did not find the most massive halo...'
      else
         write(*,*) ' ih, mh = ',imaxhalo, amaxhalo
      endif

c
c.... now label cells falling within the radius of clump
c      
      do ii = 1, npoints
         icmark(ii) = nil
      enddo

      do ii = 1, npoints
         xx  = xpn(ii)
         yy  = ypn(ii)
         zz  = zpn(ii)
         do ih = 1 , nhalo
            
            if ( amh(ih) .gt. ammin .and. 
     &           amh(ih) .lt. ammax .and. 
     &           ih .ne. imaxhalo ) then
               call PeriodicBoundaryD( xh(ih), yh(ih), zh(ih), 
     &              xx, yy, zz)
               rp = sqrt((xx-xh(ih))**2+(yy-yh(ih))**2+(zz-zh(ih))**2)
               write(*,*) rp,rh(ih),amh(ih),xh(ih),yh(ih),zh(ih)
               if ( rp .le. rh(ih) ) then
                  icmark(ii) = ih
               endif
            endif
         enddo ! ih
      enddo ! npoints

      return
      end
c
c     -------------------------------------------
      subroutine Mark_Halo_Cells3( ammin, ammax )
c     -------------------------------------------
c
c     purpose : mark cells associated with subclumps in the main halo
c 
      include 'a_hfind.h'

      parameter ( npoints = 6000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      common / MARK01 / xpn, ypn, zpn
      double precision xx, yy, zz

      do ic = 1, nctot
         ind(ic) = nil
      enddo
c
c...  Figure out ID of the most massive halo
c
      imaxhalo = 0
      amaxhalo = -1.0e30
      do ih = 1 , nhalo  
         if ( amh(ih) .ge. amaxhalo ) then
            imaxhalo = ih
            amaxhalo = amh(ih)
            xch = xh(ih)
            ych = yh(ih)
            zch = zh(ih)
            rch = 1.5 * rh(ih) 
         endif
      enddo

      if ( imaxhalo .eq. 0 ) then
         write(*,*) 'imaxhalo = ',imaxhalo
         write(*,*) ' Did not find the most massive halo...'
      else
         write(*,*) ' ih, mh = ',imaxhalo, amaxhalo
      endif
c
c.... now label cells falling within the radius of clump
c
      do ih = 1 , nhalo  
         if ( amh(ih) .gt. ammin .and. 
     &        amh(ih) .lt. ammax .and. 
     &        ih .ne. imaxhalo ) then

            ! use 431500 pts for 1x10^15Msun halo
            ! use 216000 pts for 1x10^14Msun halo
            ! use 108000 pts for 1x10^13Msun halo
            ! use 13600 pts  for 1x10^10Msun halo
            if ( amh(ih) .ge. 3.0e15 ) then
               write(*,*) ' amh(ih) = ',amh(ih),' > 3.0e15 '
               stop
            else
               npts = npoints * (amh(ih)/3.0e15)**0.5
            endif
            write(*,*) ih, amh(ih), npts, rh(ih), xh(ih), yh(ih), zh(ih)

            do i = 1 , npts
               rr    = rh(ih) * RANDd ( mrand )
               theta = pi2 * RANDd ( mrand ) - pi
               phi   = pi2 * RANDd ( mrand )
               xpn(i) = rr * cos(theta) * cos(phi) + xh(ih)
               ypn(i) = rr * cos(theta) * sin(phi) + yh(ih)
               zpn(i) = rr * sin(theta) + zh(ih)
c               xpn(i) = xh(ih) + 2.0*rh(ih) * ( RANDd ( mrand ) - 0.5 )
c               ypn(i) = yh(ih) + 2.0*rh(ih) * ( RANDd ( mrand ) - 0.5 )
c               zpn(i) = zh(ih) + 2.0*rh(ih) * ( RANDd ( mrand ) - 0.5 )      
            enddo
c     
            do i = 1 , npts
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
        
c               rp = (xx-xh(ih))**2 + (yy-yh(ih))**2 + (zz-zh(ih))**2
c               rp = sqrt(rp)
c               if ( rp .le. rh(ih) ) then
                  idcell = iFindCell ( MaxLevel , xx , yy , zz )
                  ind(idcell) = ih
c               endif
            enddo

         endif
      enddo  ! ih
      
      return
      end
c
c     --------------------------------------------------------
      subroutine Output_HSE_Profiles2 ( xc, yc, zc, rmin, rmax, 
     &                        vcut, deltavir, deltacrit, nrbin, 
     &                        l_force_center, fname )
c     --------------------------------------------------------
c
c     compute the profiles of hydrostatic mass and mach number
c
c     input:  
c
c      integer ifindcenter (see below) 
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc, yc, zc, rmin, rmax
      real vcut, deltavir, deltacrit
      integer nrbin
      logical l_force_center 
      character*256 fname 
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      real rvir, rcrit
      integer ifindcenter
      integer m, mbin, Level, icell, nfpn
      double precision amHSEfact, ptfact, xx, yy, zz
c
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints),pzIa(npoints), pzII(npoints)
      real cs(npoints), vtot(npoints), ppt(npoints)
      common / POINTS01 / xpn, ypn, zpn, pdm, pd, pp, 
     +     ptemp, pent, pzIa, pzII, vtot, cs, vxp, vyp, vzp
      parameter ( nmaxbins = 1000 )
      integer nbin(0:nmaxbins),n2bin(0:nmaxbins)
      real rbin(0:nmaxbins),dtibin(0:nmaxbins)
      real dgbin(0:nmaxbins),dstbin(0:nmaxbins),ddmbin(0:nmaxbins)
      real dgibin(0:nmaxbins),dstibin(0:nmaxbins),ddmibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),entbin(0:nmaxbins)
      real zIabin(0:nmaxbins),zIIbin(0:nmaxbins),mgbin(0:nmaxbins)
      real amHSEbin1(0:nmaxbins),amHSEbin2(0:nmaxbins),
     +     amHSEbin3(0:nmaxbins)
      real ptbin(0:nmaxbins),vrms(0:nmaxbins),csbin(0:nmaxbins)
      real dg2bin(0:nmaxbins)
      real pt2bin(0:nmaxbins),vrms2(0:nmaxbins),cs2bin(0:nmaxbins) 
      real pSGbin(0:nmaxbins),ptSGbin(0:nmaxbins),pt2SGbin(0:nmaxbins)
      real ndmp(0:nmaxbins),nstp(0:nmaxbins)
      
c     arrays for interpolations (polynomicals/cubic spline)
      real rsp(nmaxbins),psp(nmaxbins),psp2nd(nmaxbins)
      real ptsp(nmaxbins),ptsp2nd(nmaxbins)
      real pt2sp(nmaxbins),pt2sp2nd(nmaxbins)
      real dgsp(nmaxbins),dgsp2nd(nmaxbins)
      real rpol(nmaxbins),ppol(nmaxbins),dgpol(nmaxbins),cof(nmaxbins)
      real pSG(nmaxbins),ptSG(nmaxbins),pt2SG(nmaxbins),cofSG(nmaxbins)

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 
      nbins = int((rlmax-rlmin)/drl) + 1

      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a5i = 1.0 / a5
      a3i = 1.0 / a3
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      vfact  = v0 / aexpn  ! [km/s]
      ptfact = rho0 * a3i * sun_mass / (mpc*1.0e5)**3 *
     &     (1.0e5)**2  ! turbulent pressure [erg cm^-3]
      amHSEfact = -1.0 * r0 * P0 / grav_c / rho0 / aexpn *
     &     ((mpc*1.0e5)**2/sun_mass)**2  ! [1/h Msun]
      call Get_MaxLevelNow ()

c
c.... find a new center 
c
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      ifindcenter = 3 
      call Find_Halo_Radius ( xc, yc, zc, deltavir, deltacrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )

c
c...  compute the gas mass-weighted pecular velocity of the cluster & 
c     gas mass profile
c
      call Get_MaxLevelNow ()

      amgw = 0.0
      vxhalo = 0.0
      vyhalo = 0.0
      vzhalo = 0.0
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
                  if ( ibin .le. nmaxbins ) then
                     mgbin(ibin) = mgbin(ibin) + cellmass
                  endif
                  if ( rp .le. rvir ) then
                     amgw = amgw + cellmass
                     vxhalo = vxhalo + hvar(3,ic)/hvar(1,ic) * cellmass
                     vyhalo = vyhalo + hvar(4,ic)/hvar(1,ic) * cellmass
                     vzhalo = vzhalo + hvar(5,ic)/hvar(1,ic) * cellmass
                  endif
                endif
              enddo
            enddo
         ENDIF
      ENDDO

      ! the gas mass-weighted cluster peculiar velocity [code units]
      vxhalo = vxhalo / amgw
      vyhalo = vyhalo / amgw
      vzhalo = vzhalo / amgw
c      write(*,*) 'vxhalo, vyhalo, vzhalo = ',vxhalo, vyhalo, vzhalo
      
c
c.... select random samples in the volume 
c
      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin

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
c...  compute observables @ the location of each test particle
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
          vxp = hvar(3,idcell)/hvar(1,idcell)
          vyp = hvar(4,idcell)/hvar(1,idcell)
          vzp = hvar(5,idcell)/hvar(1,idcell)
c         turbulent gas flow [km/s]
          vtot(i) = (vxp-vxhalo)**2 + (vyp-vyhalo)**2 + (vzp-vzhalo)**2
          vtot(i) = vfact * sqrt(vtot(i)) ! km/s
          ppt(i)  = ptfact * pd(i) * vtot(i)**2
c         sound speed [km/s]
          cs(i) = gamma * boltz * ptemp(i) / wmu / p_mass
          cs(i) = sqrt(cs(i)) / (1.0e+5)  ! km/s
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo

c
c.... construct profiles
c
      do i = 0 , nmaxbins
        nbin(i)  = 0
        ndmp(i)  = 0

        rbin(i)    = 0.
        dgbin(i)   = 0.
        dstbin(i)  = 0.
        ddmbin(i)  = 0.
        dgibin(i)  = 0.
        dstibin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i)  = 0.
        amHSEbin1(i) = 0.
        amHSEbin2(i) = 0.
        amHSEbin3(i) = 0.

        pbin(i)    = 0.
        ptbin(i)   = 0.
        tbin(i)    = 0.
        entbin(i)  = 0.
        vrms(i)    = 0.
        csbin(i)   = 0.
        zIabin(i)  = 0.
        zIIbin(i)  = 0.

        pt2bin(i)   = 0.
        vrms2(i)    = 0.
        cs2bin(i)   = 0.

        pSGbin(i)  = 0.
        ptSGbin(i) = 0.
        pt2SGbin(i) = 0.
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

c
c...  compute gas profiles
c
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp 
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin)   = pbin(ibin)  + pp(ii)
          ptbin(ibin)  = ptbin(ibin) + ppt(ii)
          tbin(ibin)   = tbin(ibin)  + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          zIabin(ibin) = zIabin(ibin) + pzIa(ii) * pd(ii)
          zIIbin(ibin) = zIIbin(ibin) + pzII(ii) * pd(ii)
          vrms(ibin) = vrms(ibin) + vtot(ii) * pd(ii)
          csbin(ibin) = csbin(ibin) + cs(ii) * pd(ii)
          if ( vtot(ii) .lt. vcut ) then
            n2bin(ibin)  = n2bin(ibin) + 1
            dg2bin(ibin) = dg2bin(ibin) + pd(ii)
            pt2bin(ibin) = pt2bin(ibin) + ppt(ii)
            vrms2(ibin) = vrms2(ibin) + vtot(ii) * pd(ii)
            cs2bin(ibin) = cs2bin(ibin) + cs(ii) * pd(ii)
          endif
        endif
      enddo

c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          ! compute average gas profiles
          rnp = 1.0/nbin(i)
          rdp = 1.0/dgbin(i)
          rbin(i)   = rbin(i) * rnp
          dgbin(i)  = dgbin(i) * rnp
          tbin(i)   = tbin(i) * rnp
          pbin(i)   = pbin(i) * rnp
          ptbin(i)  = ptbin(i) * rnp
          entbin(i) = entbin(i) * rnp
          zIabin(i) = zIabin(i) * rdp
          zIIbin(i) = zIIbin(i) * rdp
          vrms(i)  = vrms(i) * rdp
          csbin(i) = csbin(i) * rdp
	endif

        if ( n2bin(i) .gt. 0 ) then
           rnp = 1.0/n2bin(i)
           rdp = 1.0/dg2bin(i)
           dg2bin(i) = dg2bin(i) * rnp
           pt2bin(i) = pt2bin(i) * rnp
           vrms2(i)  = vrms2(i) * rdp
           cs2bin(i) = cs2bin(i) * rdp
        endif
      enddo

c
c...  spline interpolate pressure & gas density
c
      do i = 0 , nbins
        rsp(i+1)  = rbin(i)
        psp(i+1)  = pbin(i) / P0 / a5i
        ptsp(i+1) = (pbin(i)+ptbin(i)) / P0 / a5i
        pt2sp(i+1) = (pbin(i)+pt2bin(i)) / P0 / a5i
        dgsp(i+1) = dgbin(i)
        write(*,*) pbin(i),ptbin(i),pt2bin(i)
      enddo

      ! interpolate thermal pressure
      yp1 = 0.0
      ypn = 0.0
      call spline(rsp,psp,nbins+1,yp1,ypn,psp2nd)

      ! interpolate thermal + turbulent pressure + coherent
      call spline(rsp,ptsp,nbins+1,yp1,ypn,ptsp2nd)

      ! interpolate thermal + turbulent pressure
      call spline(rsp,pt2sp,nbins+1,yp1,ypn,pt2sp2nd)

      ! interpolate gas density
      call spline(rsp,dgsp,nbins+1,yp1,ypn,dgsp2nd)

c
c...  set up arrays for the S.G. smoothing filter
c     Array needs to be uniformly space in x
c
      do i = 0 , nbins
        rr = 10.**(rlmin + float(i)*drl)
        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal + turbulent + coherent
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent  
        pSG(i+1)  = log10(pr)
        ptSG(i+1) = log10(ptr)
        pt2SG(i+1) = log10(pt2r)
        write(*,*) rr, pSG(i+1),ptSG(i+1)
      enddo
      
c
c...  compute HSE mass profiles 
c
      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      ddmi = 0.0
      dsti = 0.0
      dgi  = 0.0
      dti  = 0.0

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
        
        ! compute the HSE mass profiles
        rr = 10.**(rlmin + float(i)*drl)
        if ( rr .le. rbin(0) ) then
           rr = rbin(0)
        else if ( rr .ge. rbin(nbins) ) then
           rr = rbin(nbins)
        endif

        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal + turbulent + coherent
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent 
        call splint(rsp,dgsp,dgsp2nd,nbins+1,rr,dgr) ! gas density

        ! 1 : derivative using the Savitzky-Golay smoothing filter
        !  m=2 smooth (similar to f.d. with 4 pts on each side)
        !  m=4 ratty 

        ! thermal pressure
        mbin = 8
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pSGbin(i)  = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr )
        pSGbin(i)  = 10.**pSGbin(i) * P0 * a5i
        val1 = pSGbin(i) / P0 / a5i
        write(*,*) i, rr, val1, pr, (val1-pr)/pr
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr1 = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr1 = pr/rr * dpdr1 / drl

        ! thermal + turbulent pressure + bulk motion
        mbin = 8
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        ptSGbin(i) = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr )
        ptSGbin(i) = 10.**ptSGbin(i) * P0 * a5i
        val2 = ptSGbin(i) / P0 / a5i
        write(*,*) i, rr, val2, ptr, (val2-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr2 = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr2 = ptr/rr * dpdr2 / drl

        ! thermal + turbulent pressure 
        mbin = 8
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pt2SGbin(i) = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr )
        pt2SGbin(i) = 10.**pt2SGbin(i) * P0 * a5i
        val3 = pt2SGbin(i) / P0 / a5i
        write(*,*) i, rr, val3, ptr, (val3-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr3 = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr) 
        dpdr3 = pt2r/rr * dpdr3 / drl

        ! 2 : finite-differencing (plus-minus 4 points)
        !     4pts works the best (only small bumps and wiggles)
        !     2pts are not bad (but more sensitive to kinks)
        !     1pt is very ratty
c        r1 = 10.**(rlmin + float(i-4)*drl)
c        r2 = 10.**(rlmin + float(i+4)*drl)
c        call splint(rsp,psp,psp2nd,nbins+1,r1,p1)
c        call splint(rsp,psp,psp2nd,nbins+1,r2,p2)
c        dpdr3 = pr/rr * 
c     &       (log10(p2)-log10(p1))/(log10(r2)-log10(r1))

        ! 3 : derivative of the polynomical interpolation (in linear)
        !     CAUTION : in polcof, mbin <= 15 (mbin~8-10 should work fine)
        !     mbin=3-7 (ratty, similar to f.d. with 1p)
        !     mbin=15  (out of question, way off scale)
c        mbin = 1
c        mstart = max(i-mbin,0)
c        mend   = min(i+mbin,nbins)
c        do j = mstart , mend
c           k = j-mstart+1
c           rpol(k)  = rbin(j)
c           ppol(k)  = pbin(j) / P0 / a5i
c           dgpol(k) = dgbin(j)
c        enddo
c        mbin = mend - mstart + 1

c        call polint( rpol, ppol, mbin, rr, val1, dy )
c        write(*,*) 'e1: ',rr, val1, pr, (val1-pr)/pr, dy
c        call polcof( rpol, ppol, mbin, cof )
c        val2  = poleval( rr, mbin, cof )
c        write(*,*) 'e2: ',rr, val2, pr, (val2-pr)/pr
c        dpdr2 = polderiv( rr, mbin, cof )

        if ( dpdr1 .ge. 0.0 ) then
           dpdr1 = -1.0e-30
        endif
        if ( dpdr2 .ge. 0.0 ) then
           dpdr2 = -1.0e-30
        endif
        if ( dpdr3 .ge. 0.0 ) then
           dpdr3 = -1.0e-30
        endif
        if ( dgr .ge. 1.0e-30 ) then
           amHSEbin1(i) = amHSEfact * rr**2 / dgr * dpdr1
           amHSEbin2(i) = amHSEfact * rr**2 / dgr * dpdr2
           amHSEbin3(i) = amHSEfact * rr**2 / dgr * dpdr3
        else
           amHSEbin1(i) = 1.0e-30
           amHSEbin2(i) = 1.0e-30
           amHSEbin3(i) = 1.0e-30
        endif
                  
      enddo

c
c...  output profiles
c
      nfpn = index ( fname , ' ' ) - 1
      open ( 12 , file = fname(1:nfpn) )
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
      write(12,11) rmin , rmax , nrbin  
      write(12,12) ivir, rvir*r0, aMvir, deltavir
      write(12,22) icrit, rcrit*r0, aMcrit, deltacrit
      write(12,17) aM_dm, aM_gas, aM_st, aMvir2
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,' nrbin =',i5)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 22   format ('# icrit = ',i4,' Rcrit =',f7.4,
     &        ' [h^-1 Mpc], Mcrit =',g11.4,
     &  ' [h^-1 Msun] for Dcrit(overd w.r.t. rho_crit)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r  delta_g delta_st delta_dm  deltai_g',
     &        ' deltai_st deltai_dm  massHSE1 massHSE2 massHSE3 ',
     &        ' Tg  Pg  Sg Pturb  vrms  cs  Pturb2  vrms2  cs2 ',
     &        ' PgSG  PtSG  Pt2SG  N_dm_p' )
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
     &         rbin(i)*r0, rl*r0, rmid, rr*r0, 
     &         dgbin(i),dstbin(i),ddmbin(i),dgibin(i),dstibin(i),
     &         ddmibin(i),amHSEbin1(i),amHSEbin2(i),amHSEbin3(i),
     &         tbin(i),pbin(i),entbin(i),
     &         ptbin(i),vrms(i),csbin(i),
     &         pt2bin(i),vrms2(i),cs2bin(i),
     &         pSGbin(i),ptSGbin(i),pt2SGbin(i),int(ndmp(i))
	endif
      enddo
 18   format(25(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end
c
c     --------------------------------------------------------
      subroutine Output_HSE_Profiles3 ( xc, yc, zc, rmin, rmax, 
     &                        vcut, deltavir, deltacrit, nrbin, 
     &                        l_force_center, fname )
c     --------------------------------------------------------
c
c     compute the profiles of hydrostatic mass and mach number
c
c     input:  
c
c      integer ifindcenter (see below) 
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc, yc, zc, rmin, rmax
      real vcut, deltavir, deltacrit
      integer nrbin
      logical l_force_center 
      character*256 fname 
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      real rvir, rcrit
      integer ifindcenter
      integer m, mbin, Level, icell, nfpn
      double precision amHSEfact, ptfact, xx, yy, zz
c
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints),pzIa(npoints), pzII(npoints)
      real cs(npoints), vxp(npoints), vyp(npoints), vzp(npoints)
      common / POINTS01 / xpn, ypn, zpn, pdm, pd, pp, 
     +     ptemp, pent, pzIa, pzII, cs, vxp, vyp, vzp
      parameter ( nmaxbins = 1000 )
      integer nbin(0:nmaxbins),n2bin(0:nmaxbins)
      real rbin(0:nmaxbins),dtibin(0:nmaxbins)
      real dgbin(0:nmaxbins),dstbin(0:nmaxbins),ddmbin(0:nmaxbins)
      real dgibin(0:nmaxbins),dstibin(0:nmaxbins),ddmibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),entbin(0:nmaxbins)
      real zIabin(0:nmaxbins),zIIbin(0:nmaxbins),mgbin(0:nmaxbins)
      real amHSEbin1(0:nmaxbins),amHSEbin2(0:nmaxbins),
     +     amHSEbin3(0:nmaxbins)
      real pvx(0:nmaxbins),pvy(0:nmaxbins),pvz(0:nmaxbins)
      real pvx2(0:nmaxbins),pvy2(0:nmaxbins),pvz2(0:nmaxbins)
      real ptbin(0:nmaxbins),vrms(0:nmaxbins),csbin(0:nmaxbins)
      real dg2bin(0:nmaxbins)
      real qvx(0:nmaxbins),qvy(0:nmaxbins),qvz(0:nmaxbins)
      real qvx2(0:nmaxbins),qvy2(0:nmaxbins),qvz2(0:nmaxbins)
      real pt2bin(0:nmaxbins),vrms2(0:nmaxbins),cs2bin(0:nmaxbins) 
      real pSGbin(0:nmaxbins),ptSGbin(0:nmaxbins),pt2SGbin(0:nmaxbins)
      real ndmp(0:nmaxbins),nstp(0:nmaxbins)
      
c     arrays for interpolations (polynomicals/cubic spline)
      real rsp(nmaxbins),psp(nmaxbins),psp2nd(nmaxbins)
      real ptsp(nmaxbins),ptsp2nd(nmaxbins)
      real pt2sp(nmaxbins),pt2sp2nd(nmaxbins)
      real dgsp(nmaxbins),dgsp2nd(nmaxbins)
      real rpol(nmaxbins),ppol(nmaxbins),dgpol(nmaxbins),cof(nmaxbins)
      real pSG(nmaxbins),ptSG(nmaxbins),pt2SG(nmaxbins),cofSG(nmaxbins)

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 
      nbins = int((rlmax-rlmin)/drl) + 1

      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      a5i = 1.0 / a5
      a3i = 1.0 / a3
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      vfact  = v0 / aexpn  ! [km/s]
      ptfact = rho0 * a3i * sun_mass / (mpc*1.0e5)**3 *
     &     (vfact*1.0e5)**2  ! turbulent pressure [erg cm^-3]
      amHSEfact = -1.0 * r0 * P0 / grav_c / rho0 / aexpn *
     &     ((mpc*1.0e5)**2/sun_mass)**2  ! [1/h Msun]
      call Get_MaxLevelNow ()

c
c.... find a new center 
c
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      ifindcenter = 3 
      call Find_Halo_Radius ( xc, yc, zc, deltavir, deltacrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )

c
c...  compute the gas mass-weighted pecular velocity of the cluster & 
c     gas mass profile
c
      call Get_MaxLevelNow ()

      amgw = 0.0
      vxhalo = 0.0
      vyhalo = 0.0
      vzhalo = 0.0
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
                  if ( ibin .le. nmaxbins ) then
                     mgbin(ibin) = mgbin(ibin) + cellmass
                  endif
                  if ( rp .le. rvir ) then
                     amgw = amgw + cellmass
                     vxhalo = vxhalo + hvar(3,ic)/hvar(1,ic) * cellmass
                     vyhalo = vyhalo + hvar(4,ic)/hvar(1,ic) * cellmass
                     vzhalo = vzhalo + hvar(5,ic)/hvar(1,ic) * cellmass
                  endif
                endif
              enddo
            enddo
         ENDIF
      ENDDO

      ! the gas mass-weighted cluster peculiar velocity [code units]
      vxhalo = vxhalo / amgw
      vyhalo = vyhalo / amgw
      vzhalo = vzhalo / amgw
c      write(*,*) 'vxhalo, vyhalo, vzhalo = ',vxhalo, vyhalo, vzhalo
      
c
c.... select random samples in the volume 
c
      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin

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
c...  compute observables @ the location of each test particle
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
c         gas momentum
          vxp(i) = hvar(3,idcell)
          vyp(i) = hvar(4,idcell)
          vzp(i) = hvar(5,idcell)
c         sound speed [km/s]
          cs(i) = gamma * boltz * ptemp(i) / wmu / p_mass
          cs(i) = sqrt(cs(i)) / (1.0e+5)  ! km/s
c          vtot(i) = (vxp-vxhalo)**2 + (vyp-vyhalo)**2 + (vzp-vzhalo)**2
c          vtot(i) = vfact * sqrt(vtot(i)) ! km/s
c          ppt(i)  = ptfact * pd(i) * vtot(i)**2
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo

c
c.... construct profiles
c
      do i = 0 , nmaxbins
        nbin(i)  = 0
        ndmp(i)  = 0

        rbin(i)    = 0.
        dgbin(i)   = 0.
        dstbin(i)  = 0.
        ddmbin(i)  = 0.
        dgibin(i)  = 0.
        dstibin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i)  = 0.
        amHSEbin1(i) = 0.
        amHSEbin2(i) = 0.
        amHSEbin3(i) = 0.

        pvx(i) = 0.
        pvy(i) = 0.
        pvz(i) = 0.
        pvx2(i) = 0.
        pvy2(i) = 0.
        pvz2(i) = 0.

        pbin(i)    = 0.
        ptbin(i)   = 0.
        tbin(i)    = 0.
        entbin(i)  = 0.
        vrms(i)    = 0.
        csbin(i)   = 0.
        zIabin(i)  = 0.
        zIIbin(i)  = 0.

        qvx(i) = 0.
        qvy(i) = 0.
        qvz(i) = 0.
        qvx2(i) = 0.
        qvy2(i) = 0.
        qvz2(i) = 0.

        pt2bin(i)   = 0.
        vrms2(i)    = 0.
        cs2bin(i)   = 0.

        pSGbin(i)  = 0.
        ptSGbin(i) = 0.
        pt2SGbin(i) = 0.
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

c
c...  compute gas profiles
c
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp 
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin)   = pbin(ibin)  + pp(ii)
          tbin(ibin)   = tbin(ibin)  + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          zIabin(ibin) = zIabin(ibin) + pzIa(ii) * pd(ii)
          zIIbin(ibin) = zIIbin(ibin) + pzII(ii) * pd(ii)
          pvx(ibin)  = pvx(ibin) + vxp(ii)
          pvy(ibin)  = pvy(ibin) + vyp(ii)
          pvz(ibin)  = pvz(ibin) + vzp(ii)
          pvx2(ibin) = pvx2(ibin) + vxp(ii)**2 / pd(ii)
          pvy2(ibin) = pvy2(ibin) + vyp(ii)**2 / pd(ii)
          pvz2(ibin) = pvz2(ibin) + vzp(ii)**2 / pd(ii)
          vrms(ibin) = vrms(ibin) + ( vxp(ii)**2 + 
     +                                vyp(ii)**2 + 
     +                                vzp(ii)**2 ) / pd(ii)
          csbin(ibin) = csbin(ibin) + cs(ii) * pd(ii)

          vtot = ( (vxp(ii)/pd(ii)-vxhalo)**2 + 
     +             (vyp(ii)/pd(ii)-vyhalo)**2 + 
     +             (vzp(ii)/pd(ii)-vzhalo)**2 ) 
          vtot = vfact * sqrt(vtot) 
          if ( vtot .lt. vcut ) then
            n2bin(ibin)  = n2bin(ibin) + 1
            dg2bin(ibin) = dg2bin(ibin) + pd(ii)
            qvx(ibin)  = qvx(ibin) + vxp(ii)
            qvy(ibin)  = qvy(ibin) + vyp(ii)
            qvz(ibin)  = qvz(ibin) + vzp(ii)
            qvx2(ibin) = qvx2(ibin) + vxp(ii)**2 / pd(ii)
            qvy2(ibin) = qvy2(ibin) + vyp(ii)**2 / pd(ii)
            qvz2(ibin) = qvz2(ibin) + vzp(ii)**2 / pd(ii)
            vrms2(ibin) = vrms2(ibin) + ( vxp(ii)**2 + 
     +                                   vyp(ii)**2 + 
     +                                   vzp(ii)**2 ) / pd(ii)
            cs2bin(ibin) = cs2bin(ibin) + cs(ii) * pd(ii)
          endif
          
        endif
      enddo

c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          ! compute average gas profiles
          rnp = 1.0/nbin(i)
          rdp = 1.0/dgbin(i)
          rbin(i)   = rbin(i) * rnp
          dgbin(i)  = dgbin(i) * rnp
          tbin(i)   = tbin(i) * rnp
          pbin(i)   = pbin(i) * rnp
          entbin(i) = entbin(i) * rnp
          zIabin(i) = zIabin(i) * rdp
          zIIbin(i) = zIIbin(i) * rdp
          rhov   = sqrt( pvx(i)**2 + pvy(i)**2 + pvz(i)**2 ) * rnp
          rhov2  = sqrt( pvx2(i)**2 + pvy2(i)**2 + pvz2(i)**2 ) * rnp
          pvx(i) = pvx(i) * rdp ! mass-weighted vx
          pvy(i) = pvy(i) * rdp ! mass-weighted vy
          pvz(i) = pvz(i) * rdp ! mass-weighted vz
          vmn2 = pvx(i)**2 + pvy(i)**2 + pvz(i)**2
          vmn  = sqrt(vmn2)
          vrms(i) = vfact * sqrt(abs(vrms(i)*rdp - vmn2))
          ptbin(i) = ptfact *( rhov2 - 2.*vmn*rhov + dgbin(i)*vmn2 )
          csbin(i) = csbin(i) * rdp
	endif

        if ( n2bin(i) .gt. 0 ) then
           rnp = 1.0/n2bin(i)
           rdp = 1.0/dg2bin(i)
           dg2bin(i) = dg2bin(i) * rnp
           rhov   = sqrt( qvx(i)**2 + qvy(i)**2 + qvz(i)**2 ) * rnp
           rhov2  = sqrt( qvx2(i)**2 + qvy2(i)**2 + qvz2(i)**2 ) * rnp
           qvx(i) = qvx(i) * rdp ! mass-weighted vx
           qvy(i) = qvy(i) * rdp ! mass-weighted vy
           qvz(i) = qvz(i) * rdp ! mass-weighted vz
           vmn2 = qvx(i)**2 + qvy(i)**2 + qvz(i)**2
           vmn  = sqrt(vmn2)
           vrms2(i) = vfact * sqrt(abs(vrms2(i)*rdp - vmn2))
           pt2bin(i) = ptfact *( rhov2 - 2.*vmn*rhov + dg2bin(i)*vmn2 )
           cs2bin(i) = cs2bin(i) * rdp
        endif
      enddo

c
c...  spline interpolate pressure & gas density
c
      do i = 0 , nbins
        rsp(i+1)  = rbin(i)
        psp(i+1)  = pbin(i) / P0 / a5i
        ptsp(i+1) = (pbin(i)+ptbin(i)) / P0 / a5i
        pt2sp(i+1) = (pbin(i)+pt2bin(i)) / P0 / a5i
        dgsp(i+1) = dgbin(i)
        write(*,*) pbin(i),ptbin(i),ptbin(i)/pbin(i)
      enddo

      ! interpolate thermal pressure
      yp1 = 0.0
      ypn = 0.0
      call spline(rsp,psp,nbins+1,yp1,ypn,psp2nd)

      ! interpolate thermal + turbulent pressure + coherent
      call spline(rsp,ptsp,nbins+1,yp1,ypn,ptsp2nd)

      ! interpolate thermal + turbulent pressure
      call spline(rsp,pt2sp,nbins+1,yp1,ypn,pt2sp2nd)

      ! interpolate gas density
      call spline(rsp,dgsp,nbins+1,yp1,ypn,dgsp2nd)

c
c...  set up arrays for the S.G. smoothing filter
c     Array needs to be uniformly space in x
c
      do i = 0 , nbins
        rr = 10.**(rlmin + float(i)*drl)
        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal + turbulent + coherent
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent  
        pSG(i+1)  = log10(pr)
        ptSG(i+1) = log10(ptr)
        pt2SG(i+1) = log10(pt2r)
        write(*,*) rr, pSG(i+1),ptSG(i+1)
      enddo
      
c
c...  compute HSE mass profiles 
c
      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      ddmi = 0.0
      dsti = 0.0
      dgi  = 0.0
      dti  = 0.0

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
        
        ! compute the HSE mass profiles
        rr = 10.**(rlmin + float(i)*drl)
        if ( rr .le. rbin(0) ) then
           rr = rbin(0)
        else if ( rr .ge. rbin(nbins) ) then
           rr = rbin(nbins)
        endif

        call splint(rsp,psp,psp2nd,nbins+1,rr,pr)    ! thermal pressure
        call splint(rsp,ptsp,ptsp2nd,nbins+1,rr,ptr) ! thermal + turbulent + coherent
        call splint(rsp,pt2sp,pt2sp2nd,nbins+1,rr,pt2r) ! thermal + turbulent 
        call splint(rsp,dgsp,dgsp2nd,nbins+1,rr,dgr) ! gas density

        ! 1 : derivative using the Savitzky-Golay smoothing filter
        !  m=2 smooth (similar to f.d. with 4 pts on each side)
        !  m=4 ratty 

        ! thermal pressure
        mbin = 8
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pSGbin(i)  = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr )
        pSGbin(i)  = 10.**pSGbin(i) * P0 * a5i
        val1 = pSGbin(i) / P0 / a5i
        write(*,*) i, rr, val1, pr, (val1-pr)/pr
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr1 = SGeval( i+1, pSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr1 = pr/rr * dpdr1 / drl

        ! thermal + turbulent pressure + bulk motion
        mbin = 4
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        ptSGbin(i) = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr )
        ptSGbin(i) = 10.**ptSGbin(i) * P0 * a5i
        val2 = ptSGbin(i) / P0 / a5i
        write(*,*) i, rr, val2, ptr, (val2-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr2 = SGeval( i+1, ptSG, cofSG, nbins+1, nc, nl, nr) 
        dpdr2 = ptr/rr * dpdr2 / drl

        ! thermal + turbulent pressure 
        mbin = 8
        if ( i .lt. mbin ) then
          nl = i
          nr = mbin
        else if ( i .gt. (nbins-mbin) ) then
          nl = mbin
          nr = nbins-i
        else
          nl = mbin
          nr = nl
        endif
        nc = nl+nr+1
        call savgol( cofSG, nc, nl, nr, 0, 2 )
        pt2SGbin(i) = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr )
        pt2SGbin(i) = 10.**pt2SGbin(i) * P0 * a5i
        val3 = pt2SGbin(i) / P0 / a5i
        write(*,*) i, rr, val3, ptr, (val3-ptr)/ptr 
        call savgol( cofSG, nc, nl, nr, 1, 2 )
        dpdr3 = SGeval( i+1, pt2SG, cofSG, nbins+1, nc, nl, nr) 
        dpdr3 = pt2r/rr * dpdr3 / drl

        ! 2 : finite-differencing (plus-minus 4 points)
        !     4pts works the best (only small bumps and wiggles)
        !     2pts are not bad (but more sensitive to kinks)
        !     1pt is very ratty
c        r1 = 10.**(rlmin + float(i-4)*drl)
c        r2 = 10.**(rlmin + float(i+4)*drl)
c        call splint(rsp,psp,psp2nd,nbins+1,r1,p1)
c        call splint(rsp,psp,psp2nd,nbins+1,r2,p2)
c        dpdr3 = pr/rr * 
c     &       (log10(p2)-log10(p1))/(log10(r2)-log10(r1))

        ! 3 : derivative of the polynomical interpolation (in linear)
        !     CAUTION : in polcof, mbin <= 15 (mbin~8-10 should work fine)
        !     mbin=3-7 (ratty, similar to f.d. with 1p)
        !     mbin=15  (out of question, way off scale)
c        mbin = 1
c        mstart = max(i-mbin,0)
c        mend   = min(i+mbin,nbins)
c        do j = mstart , mend
c           k = j-mstart+1
c           rpol(k)  = rbin(j)
c           ppol(k)  = pbin(j) / P0 / a5i
c           dgpol(k) = dgbin(j)
c        enddo
c        mbin = mend - mstart + 1

c        call polint( rpol, ppol, mbin, rr, val1, dy )
c        write(*,*) 'e1: ',rr, val1, pr, (val1-pr)/pr, dy
c        call polcof( rpol, ppol, mbin, cof )
c        val2  = poleval( rr, mbin, cof )
c        write(*,*) 'e2: ',rr, val2, pr, (val2-pr)/pr
c        dpdr2 = polderiv( rr, mbin, cof )

        if ( dpdr1 .ge. 0.0 ) then
           dpdr1 = -1.0e-30
        endif
        if ( dpdr2 .ge. 0.0 ) then
           dpdr2 = -1.0e-30
        endif
        if ( dpdr3 .ge. 0.0 ) then
           dpdr3 = -1.0e-30
        endif
        if ( dgr .ge. 1.0e-30 ) then
           amHSEbin1(i) = amHSEfact * rr**2 / dgr * dpdr1
           amHSEbin2(i) = amHSEfact * rr**2 / dgr * dpdr2
           amHSEbin3(i) = amHSEfact * rr**2 / dgr * dpdr3
        else
           amHSEbin1(i) = 1.0e-30
           amHSEbin2(i) = 1.0e-30
           amHSEbin3(i) = 1.0e-30
        endif
                  
      enddo

c
c...  output profiles
c
      nfpn = index ( fname , ' ' ) - 1
      open ( 12 , file = fname(1:nfpn) )
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
      write(12,11) rmin , rmax , nrbin  
      write(12,12) ivir, rvir*r0, aMvir, deltavir
      write(12,22) icrit, rcrit*r0, aMcrit, deltacrit
      write(12,17) aM_dm, aM_gas, aM_st, aMvir2
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,' nrbin =',i5)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 22   format ('# icrit = ',i4,' Rcrit =',f7.4,
     &        ' [h^-1 Mpc], Mcrit =',g11.4,
     &  ' [h^-1 Msun] for Dcrit(overd w.r.t. rho_crit)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r  delta_g delta_st delta_dm  deltai_g',
     &        ' deltai_st deltai_dm  massHSE1 massHSE2 massHSE3 ',
     &        ' Tg  Pg  Sg Pturb  vrms  cs  Pturb2  vrms2  cs2 ',
     &        ' PgSG  PtSG  Pt2SG  N_dm_p' )
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
     &         rbin(i)*r0, rl*r0, rmid, rr*r0, 
     &         dgbin(i),dstbin(i),ddmbin(i),dgibin(i),dstibin(i),
     &         ddmibin(i),amHSEbin1(i),amHSEbin2(i),amHSEbin3(i),
     &         tbin(i),pbin(i),entbin(i),
     &         ptbin(i),vrms(i),csbin(i),
     &         pt2bin(i),vrms2(i),cs2bin(i),
     &         pSGbin(i),ptSGbin(i),pt2SGbin(i),int(ndmp(i))
	endif
      enddo
 18   format(25(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end
c
c     ----------------------------------------------------------
      real function SGeval ( ii, ySG, cofSG, nsize, nc, nl, nr )
c     ----------------------------------------------------------
c
c     evaluate values for the Savitzky-Golay smoothing filters
c
      integer nsize, ii, nc, nl, nr, k, kk
      real val
      real ySG(nsize), cofSG(nc)
c      parameter ( nmaxbins = 1000 )
c      real xSG(nmaxbins),ySG(nmaxbins),cofSG(nmaxbins)
c      common / SG01 / xSG, ySG, cofSG
       
      val = 0.0
      do k = -nl, nr
         kk = mod((nc-k),nc)+1 ! wrap-around order
         val = val + cofSG(kk) * ySG(ii+k)
c         write(*,*) ii, k, kk, cofSG(kk), ySG(ii+k), val
      enddo
      SGeval = val
      
      return
      end
c
c     -----------------------------------
      real function polderiv( x, m, cof ) 
c     -----------------------------------
c
c     evaluate derivative of the polynomial
c
      integer m
      real x, val, cof(m)

      val = 0.0
      do i = 2, m
         val = val + float(i-1) * cof(i) * x**(i-2)
      enddo
      polderiv = val
      return
      end
c
c     ----------------------------------
      real function poleval( x, m, cof ) 
c     ----------------------------------
c
c     evaluate derivative of the polynomial
c
      integer m
      real x, val, cof(m)

      val = 0.0
      do i = 1, m
         val = val + cof(i) * x**(i-1)
      enddo
      poleval = val
      return
      end

c
c     -----------------------------------------
      subroutine Output_HSEVir ( ifindcenter )
c     -----------------------------------------
c
c     compute the profiles of hydrostatic mass and mach number
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
c...  Define the output file name
c
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)
      
      lpath = index(path    , ' ') - 1
      fname  = path(1:lpath)//'/mHSEvir'//'_'//'a'//
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
c     ----------------------------------------
      subroutine Output_Velocity_Histogram ( )
c     ----------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      character*256 fname, fname2, fname4, fname8
      integer nfpn
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 

c
c...  Define the output file name
c
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)
      
      lpath = index(path, ' ') - 1
      fname = path(1:lpath)//'/vhist'//'_'//'a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_rvir.dat '
      fname2 = path(1:lpath)//'/vhist'//'_'//'a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_r2vir.dat '
      fname4 = path(1:lpath)//'/vhist'//'_'//'a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_r4vir.dat '
      fname8 = path(1:lpath)//'/vhist'//'_'//'a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_r8vir.dat '

c
c...  Find halo center and radius
c
      call Get_MaxLevelNow ()
      ifindcenter = 3
      call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      
c      oma = Om0 / ( Om0 + Oml0 * aexpn**3 )
c      xf = oma - 1.
c      dcrit = ( 18.*pi**2 + 82.*xf - 39.*xf**2 )
c      dvir  = dcrit / (1.+xf)
      dvir = 180.
      dcrit = 180.
      call Find_Halo_Radius ( xc, yc, zc, dvir, dcrit,
     &       ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &       icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &       ifindcenter )
  
      rvir2 = rvir / 2.0
      rvir4 = rvir / 4.0
      rvir8 = rvir / 8.0

      call Compute_Velocity_Histogram ( xc, yc, zc, rvir, fname )
      call Compute_Velocity_Histogram ( xc, yc, zc, rvir2, fname2 )
      call Compute_Velocity_Histogram ( xc, yc, zc, rvir4, fname4 )
      call Compute_Velocity_Histogram ( xc, yc, zc, rvir8, fname8 )

      return
      end
c
c     -----------------------------------------------------------------
      subroutine Compute_Velocity_Histogram ( xc, yc, zc, rout, fname )
c     -----------------------------------------------------------------
c
c     compute the profiles of hydrostatic mass and mach number
c
      include 'a_tree.h'
      include 'a_control.h'
c
      character*256 fname
      integer nfpn
      integer cnthalo(30000)
      real xc, yc, zc, rout
      real fgx(500), fgy(500), fgz(500)
      real fgx2(500), fgy2(500), fgz2(500)

c
c...  read halo catalog & mark cells in subhalos
c
c      write(*,*) ' Read_Halo_Catalog...'
c      call Read_Halo_Catalog( hname )
c      write(*,*) 'Mark cells in subhalos...'
c      ammin = 0.0     ! [1/h Msun]
c      ammax = 1.0e15  ! [1/h Msun]
c      do ih = 1, 30000
c         cnthalo(ih) = nil
c      enddo 
c      call Mark_Halo_Cells( ammin, ammax )
c      do ic = 1, nctot
c         cnthalo(ind(ic)) = cnthalo(ind(ic)) + 1
c      enddo
c      do ih = 1, nhalo
c         write(*,*) 'halo : ',ih, cnthalo(ih)
c      enddo

c
c...  computing histogram
c
      vfact  = v0 / aexpn  ! [km/s]
      nv   = 301
      vmax = 3000.0
      vmin = -1.0 * vmax
      dv = (vmax-vmin)/(nv-1)

      do i = 1, nv
         fgx(i) = 0.
         fgy(i) = 0.
         fgz(i) = 0.
         fgx2(i) = 0.
         fgy2(i) = 0.
         fgz2(i) = 0.
      enddo
      
      fgxt = 0.
      fgyt = 0.
      fgzt = 0.
      fgxt2 = 0.
      fgyt2 = 0.
      fgzt2 = 0.

      amgw = 0.0
      vxhalo = 0.0
      vyhalo = 0.0
      vzhalo = 0.0          
      amgw2 = 0.0
      vxhalo2 = 0.0
      vyhalo2 = 0.0
      vzhalo2 = 0.0      
      
c
c...  Loop over pixel-by-pixel
c
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
                cellmass = hvar(1,ic) * CellVolume
                
c               gas velocity [km/s]
                if ( hvar(1,ic) .gt. 0. ) then
                  vxp = vfact * hvar(3,ic) / hvar(1,ic)
                  vyp = vfact * hvar(4,ic) / hvar(1,ic)
                  vzp = vfact * hvar(5,ic) / hvar(1,ic)
                endif

                ixbin = ( vxp + vmax ) / dv + 1.5
                iybin = ( vyp + vmax ) / dv + 1.5
                izbin = ( vzp + vmax ) / dv + 1.5
 
                if ( rr .le. rout ) then
                 if ( ixbin .ge. 1 .or. ixbin .le. nv ) then
                     fgx(ixbin) = fgx(ixbin) + cellmass
                     fgxt = fgxt + cellmass
                  endif
                  if ( iybin .ge. 1 .or. iybin .le. nv ) then
                     fgy(iybin) = fgy(iybin) + cellmass
                     fgyt = fgyt + cellmass
                  endif
                  if ( izbin .ge. 1 .or. izbin .le. nv ) then
                     fgz(izbin) = fgz(izbin) + cellmass
                     fgzt = fgzt + cellmass
                  endif

                  amgw = amgw + cellmass
                  vxhalo = vxhalo + hvar(3,ic)/hvar(1,ic) * cellmass
                  vyhalo = vyhalo + hvar(4,ic)/hvar(1,ic) * cellmass
                  vzhalo = vzhalo + hvar(5,ic)/hvar(1,ic) * cellmass
                endif

                if ( rr .le. rout .and. ind(ic) .eq. nil ) then
                  if ( ixbin .ge. 1 .or. ixbin .le. nv ) then
                     fgx2(ixbin) = fgx2(ixbin) + cellmass
                     fgxt2 = fgxt2 + cellmass
                  endif
                  if ( iybin .ge. 1 .or. iybin .le. nv ) then
                     fgy2(iybin) = fgy2(iybin) + cellmass
                     fgyt2 = fgyt2 + cellmass
                  endif
                  if ( izbin .ge. 1 .or. izbin .le. nv ) then
                     fgz2(izbin) = fgz2(izbin) + cellmass
                     fgzt2 = fgzt2 + cellmass
                  endif

                  amgw2 = amgw2 + cellmass
                  vxhalo2 = vxhalo2 + hvar(3,ic)/hvar(1,ic) * cellmass
                  vyhalo2 = vyhalo2 + hvar(4,ic)/hvar(1,ic) * cellmass
                  vzhalo2 = vzhalo2 + hvar(5,ic)/hvar(1,ic) * cellmass
                endif

              endif
            enddo
          enddo
        ENDIF
      ENDDO

c
c...  normalize the histogram
c
      do i = 1, nv
         fgx(i) = fgx(i) / fgxt
         fgy(i) = fgy(i) / fgyt
         fgz(i) = fgz(i) / fgzt
         fgx2(i) = fgx2(i) / fgxt2
         fgy2(i) = fgy2(i) / fgyt2
         fgz2(i) = fgz2(i) / fgzt2
      enddo
c
c...  the gas mass-weighted cluster peculiar velocity [code units]
c
      vxhalo = vfact * vxhalo / amgw
      vyhalo = vfact * vyhalo / amgw
      vzhalo = vfact * vzhalo / amgw
      vxhalo2 = vfact * vxhalo2 / amgw2
      vyhalo2 = vfact * vyhalo2 / amgw2
      vzhalo2 = vfact * vzhalo2 / amgw2
      

c...  Writing data to the file
      nfpn = index ( fname , ' ' ) - 1
      open ( 12 , file = fname(1:nfpn) )
      write(12,11) xc , yc , zc 
      write(12,12) rout, vxhalo, vyhalo, vzhalo
      write(12,12) rout, vxhalo2, vyhalo2, vzhalo2
      write(12,13)
      do i = 1, nv
         vel = vmin + float(i-1) * dv 
         write(12,14) vel, fgx(i), fgy(i), fgz(i), 
     &        fgx2(i), fgy2(i), fgz2(i)
      enddo

 11   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 12   format (4(g13.5,1x))
 13   format ('#  velocity [km/s]       fx        fy       fz ' )
 14   format (7(g13.5,1x))
      close ( 12 )
      write(*,*) 'done writing velocity histogram...'
    
      return
      end
c
c     ----------------------------------------------
      subroutine Output_Gas_Velocities ( fprofname )
c     ----------------------------------------------
c
c     input:  
c
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      character*256 fprofname 
      integer nfpn

c
c.... calculate basic properties of the main halo 
c
      ifindcenter = 3
      dvir = 180.
      dcrit = 180.
      write(*,*) ' Find_Halo_Center...'
      call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      
      write(*,*) ' Find_Halo_Radius...'
      call Find_Halo_Radius ( xc, yc, zc, dvir, dcrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )

c      write(*,*) xc, yc, zc, rvir, rcrit, ivir, icrit

      vfact  = v0 / aexpn  ! [km/s]
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
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
                  rr = sqrt( (Posx-xc)**2+(Posy-yc)**2+(Posz-zc)**2 )
                  if ( rr .le. 1.5*rvir ) then
                     rhogas = hvar(1,ic)  
                     cellmass = hvar(1,ic) * CellVolume
                     T_code = gamma1 * hvar(8,ic) / hvar(1,ic)
                     temp = T_0 * T_code * a2i
                     vgx = vfact * hvar(3,ic) / hvar(1,ic)
                     vgy = vfact * hvar(4,ic) / hvar(1,ic)
                     vgz = vfact * hvar(5,ic) / hvar(1,ic)                  
                     write(12,10) Posx, Posy, Posz, vgx, vgy, vgz,
     &                    rhogas, cellmass, temp
                  endif
                endif
              enddo
            enddo
         ENDIF
      ENDDO
      close ( 12 )

 10   format (9(g13.4,1x))

      write(*,*) 'done writing the gas velocities. exiting...'
c
      return
      end
