CEVERINO02082006
c     Useful routines 
c       Find_Halo_Center
c       Find_Halo_Radius
c       PeriodicBoundary
c       PeriodicBoundaryD
c       age

c
c     -------------------------------------------------------
      subroutine Find_Halo_Center ( xc, yc, zc, ifindcenter )
c     Finding a new center of the halo
c     -------------------------------------------------------
c
c     ifindcenter = 0 -- gas density
c                 = 1 -- DM density
c                 = 2 -- total density
c                 = 3 -- potential
c                >= 4 -- no new center
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc 
      integer ifindcenter
c

c      write(*,*) 'Find_Halo_Center...'
      call Get_MaxLevelNow ()
c
c .... find a new center 
c
      if ( ifindcenter .le. 3 ) then

        rhogmax = -1.e6  ! max. gas density 
        rhotmax = -1.e6  ! max. total density
        phimin  = 1.e6   ! minimum of potential
        rhodmmax = -1.e6 
        irhogmax = -1
        irhotmax = -1 
        iphimin  = -1 
        irhodmmax = -1
        
        rmax = 10.0/r0
  
        DO Level = MinLevel , MaxLevelNow  
          IF ( Level .eq. MinLevel ) THEN
          ELSE 
            call Select_Cells ( Level , nLevel ) 
            do ic1 = 1 , nLevel
              icell = iSelect(ic1)
              do ic2 = 0 , 7
                ic = icell + ic2
                if ( iOctCh(ic) .eq. nil ) then 
                  call Ps ( ic , Posx , Posy , Posz )
c                  call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz)
c                  rr = sqrt( (Posx-xc)**2 +
c       &                     (Posy-yc)**2 +
c       &                     (Posz-zc)**2 )
c                  if ( rr .lt. rmax ) then                   
                    dcell1 = var(1,ic)*2.0**(3.*Level) + 1.0 ! tot. density
                    dcell2 = var(2,ic)  ! potential 
                    dcell3 = hvar(1,ic) ! gas density
                    dcell4 = ref(ic) ! mass of DM + stars
c                    write(*,*) dcell1, dcell2, dcell3, dcell4
                    if ( dcell1 .gt. rhotmax ) then 
                      rhotmax  = dcell1 
                      irhotmax = ic 
                    endif
                    if ( dcell2 .lt. phimin ) then 
                      phimin  = dcell2
                      iphimin = ic 
                    endif
                    if ( dcell3 .gt. rhogmax ) then
                      rhogmax  = dcell3
                      irhogmax = ic
                    endif
                    if ( dcell4 .gt. rhodmmax ) then 
                      rhodmmax = dcell4
                      irhodmmax = ic
                    endif
c                  endif
                endif
              enddo
            enddo
          ENDIF
        ENDDO      
c       
        if ( ifindcenter .eq. 0 ) then 
          if ( irhogmax .eq. -1 ) then 
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find gas density maximum'
            write(*,*) '   will use input center '
          else
            call Ps ( irhogmax , xc1 , yc1 , zc1 )
            write(*,*) 'input center x, y, z:', xc, yc, zc 
            write(*,*) 'identified gas density maximum rhogmax=',rhogmax
            write(*,*) 'at x, y, z =', xc1, yc1, zc1            
            write(*,*) 'will use it to construct profiles' 
            xc = xc1
            yc = yc1 
            zc = zc1
          endif
        endif
        if ( ifindcenter .eq. 1 ) then
          if ( irhodmmax .eq. -1 ) then
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find DM density maximum'
            write(*,*) '   will use input center '
          else
            call Ps ( irhodmmax , xc1 , yc1 , zc1 )
            write(*,*) 'input center x, y, z:', xc, yc, zc
            write(*,*)'identified DM density maximum rhodmmax=',rhodmmax
            write(*,*) 'at x, y, z =', xc1, yc1, zc1
            write(*,*) 'will use it to construct profiles'
            xc = xc1
            yc = yc1
            zc = zc1
          endif
        endif
        if ( ifindcenter .eq. 2 ) then  ! total density
          if ( irhotmax .eq. -1 ) then 
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find total density maximum'
            write(*,*) '   will use input center '
          else
            call Ps ( irhotmax , xc1 , yc1 , zc1 )
            write(*,*) 'input center x, y, z:', xc, yc, zc 
            write(*,*) 'identified tot density maximum rhotmax=',rhotmax
            write(*,*) 'at x, y, z =', xc1, yc1, zc1            
            write(*,*) 'will use it to construct profiles' 
            xc = xc1
            yc = yc1 
            zc = zc1
          endif
        endif
        if ( ifindcenter .eq. 3 ) then  ! potential 
          if ( iphimin .eq. -1 ) then 
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find potential minimum'
            write(*,*) '   will use input center '
          else
            call Ps ( iphimin , xc1 , yc1 , zc1 )
c            write(*,*) 'input center x, y, z:', xc, yc, zc 
c            write(*,*) 'identified minimum of the potential=',phimin
c            write(*,*) 'at x, y, z =', xc1, yc1, zc1            
c            write(*,*) 'will use it to construct profiles' 
            xc = xc1
            yc = yc1 
            zc = zc1
          endif
        endif
        
      endif

      return
      end
c
c     --------------------------------------------------------------
      subroutine Find_Halo_Radius ( xc, yc, zc, deltavir, deltacrit, 
     &     ivir, rvir, aMvir, aMvir_dm, aMvir_gas, aMvir_st, 
     &     icrit, rcrit, aMcrit, aMcrit_dm, aMcrit_gas, aMcrit_st, 
     &     ifindcenter )
c     Finding radius of the halo
c     --------------------------------------------------------------
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
      real zIabin(0:nmaxbins),zIIbin(0:nmaxbins)
      real ndmp(0:nmaxbins),nstp(0:nmaxbins),entbin(0:nmaxbins)
      real mgi
      real mgbin(0:nmaxbins)
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

      rmin = 1.0e-3/r0 ! 1h^-1 kpc
      rmax = 10.0/r0   ! 10h^-1 Mpc
      nrbin = 80

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
c ...  find a new center
c
      call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      
      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin
c
c .... select random samples in the volume 
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
c ....   gas density
c
        pd(i) = hvar(1,idcell)
c
c ....   pressure
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
c  . ... construct profiles
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
        nbin(i) = 0
        ndmp(i) = 0
        mgbin(i) = 0.  
      enddo

c
c  ...  compute profiles of DM and Stars
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

c  ...  compute gas profiles

      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
        endif
      enddo

c
c  ...  compute gas mass 
c

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
                  if ( ibin .le. nmaxbins ) then
                     mgbin(ibin) = mgbin(ibin) + cellmass
                  endif
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
	endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      ddmi = 0.0
      dsti = 0.0
      dgi  = 0.0
      dti  = 0.0
      mgi  = 1.0e-30

c      dcritave = deltacrit * aexpn**3 / Om0 
      dcritave = deltacrit * ( Om0 + Oml0 * aexpn**3 + 
     &     (1.0-Om0-Oml0) * aexpn ) / Om0  
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

c      aMvir = deltavir * 
c     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
c      aMcrit = dcritave * 
c     &        4.18879 * (rcrit*r0)**3 * rho0 / hubble**2 
      
      rr = 10.**(rlmin + float(ivir)*drl)
      aMvir_dm  = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
      aMvir_gas = dgibin(ivir)  * (4.18879 * rr**3) * aM0 * hubble
      aMvir_st  = dstibin(ivir) * (4.18879 * rr**3) * aM0 * hubble 
      aMvir = aMvir_dm + aMvir_gas + aMvir_st 

      rr = 10.**(rlmin + float(icrit)*drl)
      aMcrit_dm  = ddmibin(icrit) * (4.18879 * rr**3) * aM0 * hubble
      aMcrit_gas = dgibin(icrit)  * (4.18879 * rr**3) * aM0 * hubble
      aMcrit_st  = dstibin(icrit) * (4.18879 * rr**3) * aM0 * hubble 
      aMcrit = aMcrit_dm + aMcrit_gas + aMcrit_st 

      return
      end
c
c     -----------------------------------------------------------
      subroutine PeriodicBoundary ( xc, yc, zc, Posx, Posy, Posz)
c     Take care of periodic boundary conditions for analysis
c     -----------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc, yc, zc, Posx, Posy, Posz
c      double precision xc, yc, zc, Posx, Posy, Posz
c
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

      return
      end
c
c     -----------------------------------------------------------
      subroutine PeriodicBoundaryD ( xc, yc, zc, xx, yy, zz)
c     Take care of periodic boundary conditions for analysis
c     -----------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc, yc, zc
      double precision xx, yy, zz
c
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
      
      xx = xx + corr_x
      yy = yy + corr_y
      zz = zz + corr_z

      return
      end

c     ---------------------------
      real*8 function age ( td )
c     ---------------------------
c
c     returns age of the Universe at td (time in code units)
c
c     uses Oleg's formula for flat LCDM 
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 td, b2a

      as = b2a ( td ) ! convert code time to expansion factor
      ff = Om0/(1.-Om0)/as**3 
c
c  .... calculate age of the universe in Gyrs 
c
      age = 9.779/hubble * 2./3./
     &      sqrt(1.-Om0)*log((1.+sqrt(1.+ff))/sqrt(ff))

      return
      end
