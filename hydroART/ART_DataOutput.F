c     Output_Stars_to_HOP
c     DataOutput_2Dmap
c     Extract_Particles()
c     Extract_Volume()
c
c     --------------------------------------------
      subroutine Output_Stars_to_HOP ( fprofname )
c     --------------------------------------------
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

      nfpn = index ( fprofname , ' ' ) - 1

      nstars = nsp(nspecies,2) - nsp(nspecies,1) + 1
      istar1 = nsp(nspecies,1)

      open ( 12 , file = fprofname(1:nfpn) )

      write(12,*) nstars
      do i = 1 , np 
        if ( i .gt. istar1 ) then ! stellar particle
          write(12,*) x(i), y(i), z(i), pw(i)
        endif
      enddo
      close ( 12 )
      write(*,*) 'done writing the stellar particles. exiting...'
c
      return
      end
c
c     ---------------------------------------------------------------
      subroutine DataOutput_2Dmap ( iflag )
c     1. Output average profiles : Output_Average_Profiles()
c     2. Output projected 2D map 
c          1 = total surface density (DM+Gas+Stars)
c          2 = surface mass density of DM only
c          3 = X-ray surface brightness
c          4 = emission weighted temperature
c          5 = emission weighted metallicity
c          6 = mass weighted temperature
c          7 = mass weighted metallicity
c          8 = thermal SZE (Compton-y)
c          9 = kinetic SZE (b-parameter)
c          10 = OVII  561eV [erg/s/cm^2]
c          11 = OVII  568eV [erg/s/cm^2]
c          12 = OVII  568eV [erg/s/cm^2]
c          13 = OVII  665eV [erg/s/cm^2]
c          14 = OVIII 653eV [erg/s/cm^2]
c 
c          Use the smallest grid as a basis for uniform grid
c          Compute mass within an uniform cell in the projected plane
c          Compute the surface mass density [Msun Mpc^-2]
c
c     CAUTION : nsize<=1024 for X-ray analysis involving X-ray spectra 
c
c     ---------------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer nspec
      parameter ( nsize = 1024 )
      integer iflag
      real p_d1(nsize,nsize), p_d2(nsize,nsize), p_d3(nsize,nsize)
      real q_d1(nsize,nsize), q_d2(nsize,nsize), q_d3(nsize,nsize)
      common / MAP1 / p_d1, p_d2, p_d3
      common / MAP2 / q_d1, q_d2, q_d3
c      real p_t1(nsize,nsize), p_t2(nsize,nsize), p_t3(nsize,nsize)
      real*8 sum1, sum2, sum3, sum4, vmax1, vmax2, em
      double precision emiss
      character*256 fname, fnh
      logical lpointsout, l_force_center
      integer imovstep
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*5 chbox
      character*256 fname

      real EeV1, EeV2, dE
      common / RS03 / Xpri, Xsol,EeV1,EeV2,nrs,dene,SOLA,defactor2
      common / RS04 / NET, NZT, NTT, dE, RSPEC

c
c.... output profiles 
c
      xc = 16.53711
      yc = 22.15430
      zc = 128.7441 
      rmin = 1.e-3/r0
      rmax = 10.0/r0
      nrbin = 45 
      ibintype = 0 
      deltavir = 180.
      deltacrit = 180.
      ifindcenter = 3 
      lpointsout = .false.  
      l_force_center = .false.  

      lpath = index(path    , ' ') - 1

        ic = int(aexpn*1000)
        i1 = ic / 1000 
        i2 = (ic - (ic/1000)*1000)/100 
        i3 = (ic - (ic/100)*100) / 10 
        i4 = (ic - (ic/10)*10)

      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif
      
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

c...  Normalization of gas variables
      t_fact = gamma1 * T_0 * a2i
      
c...  X-ray emission [RS in units of n_e n_H 10^-23 erg s^-1 cm^-3]
      if ( iflag .ge. 3 .and. iflag .le. 5 ) then
        Ynow = 0.2496
        zsim = 1.0/aexpn - 1.0
        if ( zsim .le. 1.0e-20 ) zsim = 0.0
        ! Normalization for XSB [erg s^-1 cm^-2 arcmin^-2]
        fy = (1.0 - 0.5 * Ynow) * (1.0 - Ynow)
        ex_fact = 3.89535e-9 * fy * a5i * hubble**3 * Om0**2 * r0
        ex_fact = ex_fact / (4.0*pi*(1+zsim)**4) * (4.0*(pi/21600)**2)

c        fname  = '/dworkin/home/daisuke/ARTsim/Data/RS/' //
c     &       'rs_spec_0.01-10.01_keV.dat'
c        call Read_RS_Spec ( fname )

        fname  = '/dworkin/home/daisuke/Raymond/DAT/'//
     &       'rs_det_chandra_acis-s-bi_0.3_7.0keV.dat'        
        if ( iflag .eq. 3 ) call Read_RSdet_Table ( 1, fname )
        if ( iflag .eq. 4 .or. iflag .eq. 5 ) 
     &       call Read_RSdet_Table ( 0, fname )
        write(*,*) 'Finished reading the RS spectra...'
      endif

c...  factors for thermal and kinematic SZ effects
      szt_fact = 3.41262e-15 * r0 * T_0 * hubble * Om0 * a4i
      szk_fact = 6.75012e-11 * r0 * v0 * hubble * Om0 * a3i

c...  Read in the Oxygen emission line emissivities [erg cm^3 s^-1]
      imodel = 2 ! 0 = simulation, 1 = 0.2 Zsun, 2 = 0.02 * (rho/rho_b)^0.3
      if ( iflag .eq. 10 )  call load_oxygen(1)
      if ( iflag .eq. 11 )  call load_oxygen(2)
      if ( iflag .eq. 12 )  call load_oxygen(3)
      if ( iflag .eq. 13 )  call load_oxygen(4)
      if ( iflag .eq. 14 )  call load_oxygen(5)

c      Tcell = 100000.0
c      call lambda_oxygen( Tcell, emiss )
c      write(*,*) Tcell, emiss

c      Tcell = 680397.2 
c      call lambda_oxygen( Tcell, emiss )
c      write(*,*) Tcell, emiss

c      Tcell = 700000.0
c      call lambda_oxygen( Tcell, emiss )
c      write(*,*) Tcell, emiss

c      Tcell = 1000000.0
c      call lambda_oxygen( Tcell, emiss )
c      write(*,*) Tcell, emiss
c      stop
      
      oxgen_fact = 2.99944e+14 * hubble**3 * Om0**2 * r0 * a5i
      oxgen_fact = oxgen_fact / 4.0 / pi * aexpn**4

c
c.... define a box
c
      icenter = int(nsize/2.+small)
      csMaxLev  = CellSize(MaxLevelNow) 
      csMaxLev2 = csMaxLev/2.
      dx = csMaxLev * (r0/hubble)  ! pixel size [Mpc]
      sum4 = 0.
      
c     initialize 2d maps
      do iy = 1 , nsize
        do ix = 1 , nsize
          p_d1(ix,iy) = 0.0
          p_d2(ix,iy) = 0.0
          p_d3(ix,iy) = 0.0
          q_d1(ix,iy) = 0.0
          q_d2(ix,iy) = 0.0
          q_d3(ix,iy) = 0.0
        enddo
      enddo

c
c...  Loop over pixel-by-pixel
c

      lstart = MinLevel
      lend   = MaxLevelNow
      if ( iflag .ge. 3 .and. iflag .le. 7 .or. iflag .ge. 10 ) then
        lstart = lstart + 1
      endif
      

      DO Level = lstart , lend
        CellVolume = 1.0 * 2.0**(-3.0*Level)
        CellArea = 1.0 * 2.0**(-2.0*Level)
        cs  = CellSize(Level)
        cs2 = CellSize(Level)/2.       
  
        IF ( Level .eq. MinLevel ) THEN

          write(*,*) 'Level = ',Level,' nLevel = ',ncell0
          do ic = 1 , ncell0
            if ( iOctCh(ic) .eq. nil ) then
              call Ps ( ic , Posx , Posy , Posz )
              call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz)
              
              ixmin = icenter + int((Posx-xc-cs2)/csMaxLev+0.5) 
              ixmax = icenter + int((Posx-xc+cs2)/csMaxLev+0.5)     
              iymin = icenter + int((Posy-yc-cs2)/csMaxLev+0.5) 
              iymax = icenter + int((Posy-yc+cs2)/csMaxLev+0.5) 
              izmin = icenter + int((Posz-zc-cs2)/csMaxLev+0.5)
              izmax = icenter + int((Posz-zc+cs2)/csMaxLev+0.5)
              
              idx = ixmax - ixmin
              idy = iymax - iymin
              idz = izmax - izmin
              
              if ( idx .ne. idy ) then
                 write(*,*) 'STOP : idx=',idx,' not equal to idy=',idy
                 write(*,*) ' ixmax=',ixmax,' ixmin=',ixmin
                 write(*,*) ' iymax=',iymax,' iymin=',iymin  
                 write(*,*) Posx,Posy,xc,yc,cs2
                 stop 
              endif
              if ( idx .ne. idz ) then
                 write(*,*) 'STOP : idx=',idx,' not equal to idz=',idz
                 write(*,*) ' ixmax=',ixmax,' ixmin=',ixmin
                 write(*,*) ' izmax=',izmax,' izmin=',izmin
                 write(*,*) Posx,Posz,xc,zc,cs2
                 stop
              endif
              if ( idy .ne. idz ) then
                 write(*,*) 'STOP : idy=',idy,' not equal to idz=',idz
                 write(*,*) ' iymax=',iymax,' iymin=',iymin
                 write(*,*) ' izmax=',izmax,' izmin=',izmin
                 write(*,*) Posy,Posz,yc,zc,cs2
                 stop
              endif
             
              fcell = cs / idx / idy   

c...          compute the cell variables
              if ( iflag .eq. 1 ) then
                 ! compute total mass in each cell [Msun]
                 ! m = (var(1,ic) + cellvolume) * aM0 is correct, not * (times)!
                 ax = (var(1,ic) + CellVolume) * aM0
              else if ( iflag .eq. 2 ) then
                 ax = ref(ic) * aM0
              else if ( iflag .eq. 6 ) then
                 Tcell = t_fact *  hvar(8,ic) / hvar(1,ic)
                 bx = hvar(1,ic) * fcell
                 ax = Tcell * bx
              else if ( iflag .eq. 7 ) then
                 Zcell = (hvar(izIa,ic)+hvar(izII,ic)) / hvar(1,ic)
                 Zcell = Zcell / 0.01989
                 bx = hvar(1,ic) * fcell
                 ax = Zcell * bx
              else if ( iflag .eq. 8 ) then
                 ax = szt_fact * gamma1 * hvar(8,ic) * fcell
              else if ( iflag .eq. 9 ) then
                 ax = szk_fact * hvar(3,ic) * fcell 
                 ay = szk_fact * hvar(4,ic) * fcell 
                 az = szk_fact * hvar(5,ic) * fcell 
              endif
              
              if ( iflag .ne. 9 ) then
                 ay = ax
                 az = ax
                 by = bx
                 bz = bx
              endif

              sum4 = sum4 + ax

c...          projection along x-axis   
              if ( iymin .gt. nsize .or. iymax .lt. 1 .or. 
     &             izmin .gt. nsize .or. izmax .lt. 1 ) then
              else
                do iz = izmin, izmax-1
                   do iy = iymin, iymax-1
                      if ( iy .ge. 1 .and. iy .le. nsize .and.
     &                     iz .ge. 1 .and. iz .le. nsize ) then
                         p_d1(iy,iz) = p_d1(iy,iz) + ax
                         q_d1(iy,iz) = q_d1(iy,iz) + bx
                      endif
                   enddo
                enddo
              endif
              
c...          projection along y-axis
              if ( izmin .gt. nsize .or. izmax .lt. 1 .or. 
     &             ixmin .gt. nsize .or. ixmax .lt. 1 ) then
              else
                 do ix = ixmin, ixmax-1
                    do iz = izmin, izmax-1
                       if ( iz .ge. 1 .and. iz .le. nsize .and.
     &                      ix .ge. 1 .and. ix .le. nsize ) then
                          p_d2(iz,ix) = p_d2(iz,ix) + ay
                          q_d2(iz,ix) = q_d2(iz,ix) + by
                       endif
                    enddo
                 enddo
              endif

c...          projection along z-axis
              if ( ixmin .gt. nsize .or. ixmax .lt. 1 .or. 
     &             iymin .gt. nsize .or. iymax .lt. 1 ) then
              else
                 do iy = iymin, iymax-1
                    do ix = ixmin, ixmax-1
                       if ( ix .ge. 1 .and. ix .le. nsize .and.
     &                      iy .ge. 1 .and. iy .le. nsize ) then
                         p_d3(ix,iy) = p_d3(ix,iy) + az 
                         q_d3(ix,iy) = q_d3(ix,iy) + bz
                       endif
                    enddo
                 enddo
              endif
              
            endif   
          enddo  ! ic
        

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
                ixmin = icenter + int((Posx-xc-cs2)/csMaxLev+0.5) 
                ixmax = icenter + int((Posx-xc+cs2)/csMaxLev+0.5)     
                iymin = icenter + int((Posy-yc-cs2)/csMaxLev+0.5) 
                iymax = icenter + int((Posy-yc+cs2)/csMaxLev+0.5) 
                izmin = icenter + int((Posz-zc-cs2)/csMaxLev+0.5)
                izmax = icenter + int((Posz-zc+cs2)/csMaxLev+0.5)
                
                idx = ixmax - ixmin
                idy = iymax - iymin
                idz = izmax - izmin

                if ( idx .ne. idy ) then
                   write(*,*) 'STOP : idx=',idx,' not equal to idy=',idy
                   write(*,*) ' ixmax=',ixmax,' ixmin=',ixmin
                   write(*,*) ' iymax=',iymax,' iymin=',iymin  
                   write(*,*) Posx,Posy,xc,yc,cs2
                   stop 
                endif
                if ( idx .ne. idz ) then
                   write(*,*) 'STOP : idx=',idx,' not equal to idz=',idz
                   write(*,*) ' ixmax=',ixmax,' ixmin=',ixmin
                   write(*,*) ' izmax=',izmax,' izmin=',izmin
                   write(*,*) Posx,Posz,xc,zc,cs2
                   stop
                endif
                if ( idy .ne. idz ) then
                   write(*,*) 'STOP : idy=',idy,' not equal to idz=',idz
                   write(*,*) ' iymax=',iymax,' iymin=',iymin
                   write(*,*) ' izmax=',izmax,' izmin=',izmin
                   write(*,*) Posy,Posz,yc,zc,cs2
                   stop
                endif

                fcell = cs / idx / idy  
                
c...            compute the cell variables
                if ( iflag .eq. 1 ) then
                   ! compute total mass in each cell [Msun]
                   ! m = (var(1,ic) + cellvolume) * aM0 is correct, not * (times)!
                   ax = (var(1,ic) + CellVolume) * aM0
                else if ( iflag .eq. 2 ) then
                   ax = ref(ic) * aM0
                else if ( iflag .ge. 3 .and. iflag .le. 5 ) then
                  Tcell = t_fact *  hvar(8,ic) / hvar(1,ic)
                  Zcell = (hvar(izIa,ic)+hvar(izII,ic)) / hvar(1,ic)
                  Zcell = Zcell / 0.01989  
                  Tlog = log10(Tcell)
                  Zlog = log10(Zcell)
                  bx = ex_fact * hvar(1,ic)**2 * RS(Tlog,Zlog) * fcell
                  if ( iflag .eq. 4 ) ax = Tcell * bx 
                  if ( iflag .eq. 5 ) ax = Zcell * bx
                else if ( iflag .eq. 6 ) then
                   Tcell = t_fact *  hvar(8,ic) / hvar(1,ic)
                   bx = hvar(1,ic) * fcell
                   ax = Tcell * bx
                else if ( iflag .eq. 7 ) then
                   Zcell = (hvar(izIa,ic)+hvar(izII,ic)) / hvar(1,ic)
                   Zcell = Zcell / 0.01989
                   bx = hvar(1,ic) * fcell
                   ax = Zcell * bx
                else if ( iflag .eq. 8 ) then
                   ax = szt_fact * gamma1 * hvar(8,ic) * fcell 
                else if ( iflag .eq. 9 ) then
                   ax = szk_fact * hvar(3,ic) * fcell 
                   ay = szk_fact * hvar(4,ic) * fcell 
                   az = szk_fact * hvar(5,ic) * fcell 
                else if ( iflag .ge. 10 .and. iflag .le. 14 ) then
                   Tcell = t_fact * hvar(8,ic) / hvar(1,ic)
                   Zcell = hvar(izIa,ic) + hvar(izII,ic)
                   if (imodel.eq.0) Zcell = Zcell / 0.01989 / hvar(1,ic)
                   if (imodel.eq.1) Zcell = 0.2
                   if (imodel.eq.2) Zcell = 0.02 * 
     &                  (hvar(1,ic)*Om0/Omb0)**(0.3)
                   call lambda_oxygen( Tcell, emiss )
                   ax = oxgen_fact * hvar(1,ic)**2 * emiss * Zcell 
     &                  * fcell
                   if ( ax .lt. -1.0e-20) then
                      write(*,*) oxgen_fact,hvar(1,ic),Tcell,
     &                     Zcell,emiss,ax
                   endif
                endif

                if ( iflag .ne. 9 ) then
                 ay = ax
                 az = ax
                 by = bx
                 bz = bx
                endif
                
c...            projection along x-axis   
                if ( iymin .gt. nsize .or. iymax .lt. 1 .or. 
     &               izmin .gt. nsize .or. izmax .lt. 1 ) then
                else
                   do iz = izmin, izmax-1
                      do iy = iymin, iymax-1
                         if ( iy .ge. 1 .and. iy .le. nsize .and.
     &                        iz .ge. 1 .and. iz .le. nsize ) then
                            p_d1(iy,iz) = p_d1(iy,iz) + ax
                            q_d1(iy,iz) = q_d1(iy,iz) + bx
                         endif
                      enddo
                   enddo
                endif

c...            projection along y-axis
                if ( izmin .gt. nsize .or. izmax .lt. 1 .or. 
     &               ixmin .gt. nsize .or. ixmax .lt. 1 ) then
                else
                   do ix = ixmin, ixmax-1
                      do iz = izmin, izmax-1
                         if ( iz .ge. 1 .and. iz .le. nsize .and.
     &                        ix .ge. 1 .and. ix .le. nsize ) then
                            p_d2(iz,ix) = p_d2(iz,ix) + ay
                            q_d2(iz,ix) = q_d2(iz,ix) + by
                         endif
                       enddo
                   enddo
                endif
              
c...            projection along z-axis
                if ( ixmin .gt. nsize .or. ixmax .lt. 1 .or. 
     &               iymin .gt. nsize .or. iymax .lt. 1 ) then
                else
                   do iy = iymin, iymax-1
                      do ix = ixmin, ixmax-1
                         if ( ix .ge. 1 .and. ix .le. nsize .and.
     &                        iy .ge. 1 .and. iy .le. nsize ) then
                            p_d3(ix,iy) = p_d3(ix,iy) + az 
                            q_d3(ix,iy) = q_d3(ix,iy) + bz
                         endif
                      enddo
                   enddo
                endif

              endif
            enddo !ic2
          enddo  !ic1

        ENDIF
      ENDDO

c
c...  compute the projected quantities 
c
      pix_area = dx**2  ! pixel area in [Mpc^2]
      do iy = 1 , nsize
        do ix = 1 , nsize
          if ( iflag .eq. 1 .or. iflag .eq. 2 ) then
            p_d1(ix,iy) = p_d1(ix,iy) / pix_area ! surface mass density [Msun Mpc^-2]
            p_d2(ix,iy) = p_d2(ix,iy) / pix_area
            p_d3(ix,iy) = p_d3(ix,iy) / pix_area
          else if ( iflag .eq. 3 ) then
            p_d1(ix,iy) = q_d1(ix,iy) 
            p_d2(ix,iy) = q_d2(ix,iy)
            p_d3(ix,iy) = q_d3(ix,iy) 
          else if ( iflag .ge. 4 .and. iflag .le. 7 ) then
            p_d1(ix,iy) = p_d1(ix,iy) / q_d1(ix,iy) 
            p_d2(ix,iy) = p_d2(ix,iy) / q_d2(ix,iy)
            p_d3(ix,iy) = p_d3(ix,iy) / q_d3(ix,iy)
          endif
        enddo
      enddo

c
c...  computing mass within the projected map [Msun]
c

      if ( iflag .ge. 1 .and. iflag .le. 2 ) then
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        do iy = 1 , nsize
          do ix = 1 , nsize
            sum1 = sum1 + p_d1(ix,iy) 
            sum2 = sum2 + p_d2(ix,iy) 
            sum3 = sum3 + p_d3(ix,iy) 
          enddo
        enddo
        write(*,10) dx*nsize,sum1,sum2,sum3
      endif
 10   format ('Lbox [Mpc] = ',g9.3,
     &     ' (Sx,Sy,Sz) = ', 3(g13.5,1x))

c
c...  computing mean, minimum, maxmum of the maps
c

      if ( iflag .gt. 2 ) then
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        smin1 = 1.0e30
        smin2 = smin1
        smin3 = smin1
        smax1 = -1.0e30
        smax2 = smax1
        smax3 = smax1
        do iy = 1 , nsize
          do ix = 1 , nsize
            sum1 = sum1 + p_d1(ix,iy) 
            sum2 = sum2 + p_d2(ix,iy)
            sum3 = sum3 + p_d3(ix,iy)
            if ( smin1 .gt. p_d1(ix,iy) ) smin1 = p_d1(ix,iy)
            if ( smin2 .gt. p_d2(ix,iy) ) smin2 = p_d2(ix,iy)
            if ( smin3 .gt. p_d3(ix,iy) ) smin3 = p_d3(ix,iy)
            if ( smax1 .lt. p_d1(ix,iy) ) smax1 = p_d1(ix,iy)
            if ( smax2 .lt. p_d2(ix,iy) ) smax2 = p_d2(ix,iy)
            if ( smax3 .lt. p_d3(ix,iy) ) smax3 = p_d3(ix,iy)
          enddo
        enddo
        write(*,13) sum1/nsize**2, smin1, smax1
        write(*,14) sum2/nsize**2, smin2, smax2
        write(*,15) sum3/nsize**2, smin3, smax3
      endif
 13   format ('x : mean = ',(g13.5,1x),'min = ',(g13.5,1x),
     &     ' max = ',(g13.5,1x))
 14   format ('y : mean = ',(g13.5,1x),'min = ',(g13.5,1x),
     &     ' max = ',(g13.5,1x))
 15   format ('z : mean = ',(g13.5,1x),'min = ',(g13.5,1x),
     &     ' max = ',(g13.5,1x))
  
c
c...  Output ASCII image data to the file 
c

      lpath = index(path    , ' ') - 1  
      if ( iflag .eq. 1 ) fnh = '/sigma_tot'
      if ( iflag .eq. 2 ) fnh = '/sigma_dm'
      if ( iflag .eq. 3 ) fnh = '/xsb'
      if ( iflag .eq. 4 ) fnh = '/tew'
      if ( iflag .eq. 5 ) fnh = '/zew'
      if ( iflag .eq. 6 ) fnh = '/tm'
      if ( iflag .eq. 7 ) fnh = '/zm'
      if ( iflag .eq. 8 ) fnh = '/tsz'
      if ( iflag .eq. 9 ) fnh = '/ksz'
      if ( iflag .eq. 10 ) fnh = '/OVII_561eV_ZDen'
      if ( iflag .eq. 11 ) fnh = '/OVII_568eV_ZDen'
      if ( iflag .eq. 12 ) fnh = '/OVII_574eV_ZDen'
      if ( iflag .eq. 13 ) fnh = '/OVII_665eV_ZDen'
      if ( iflag .eq. 14 ) fnh = '/OVIII_653eV_ZDen'
      lfnh = index(fnh    , ' ') - 1
      
      fname = path(1:lpath)//fnh(1:lfnh)//'_a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_x'
      nfname = index ( fname , ' ' ) - 1
      write(*,*) 'writing...'
      open ( 20 , file = fname(1:nfname) )
      write(20,*) aexpn, nsize, nsize
      write(20,*) dx, dx
      write(20,*) p_d1
      close ( 20 )
      write(*,*) 'done.'
      
      fname = path(1:lpath)//fnh(1:lfnh)//'_a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_y'
      nfname = index ( fname , ' ' ) - 1
      write(*,*) 'writing...'
      open ( 20 , file = fname(1:nfname) )
      write(20,*) aexpn, nsize, nsize
      write(20,*) dx, dx
      write(20,*) p_d2
      close ( 20 )
      write(*,*) 'done.'
      
      fname = path(1:lpath)//fnh(1:lfnh)//'_a'//
     &     digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &     //'_z'
      nfname = index ( fname , ' ' ) - 1
      write(*,*) 'writing...'
      open ( 20 , file = fname(1:nfname) )
      write(20,*) aexpn, nsize, nsize
      write(20,*) dx, dx
      write(20,*) p_d3
      close ( 20 )
      write(*,*) 'done.'

      return
      end
