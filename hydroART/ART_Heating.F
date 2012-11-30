c     ========================================================================
c     
c           Routines for AGN feedback & Preheating
c
c            Daisuke Nagai (October 2003)
c
c     ========================================================================
c
c     AGN_SETUP        : activate/reset AGN feedback parameters
c       Initialize_AGN   : initialize basic parameters for AGN feedback 
c       Setup_Jets       : set up bipolar mehcnical jets from AGN
c       E61_AGN          : return total energy of the AGN feedback [10^61 erg] 
c       t7_AGN           : return the duty cycle of AGN [10^7 yrs]
c       fke_AGN          : returns the fraction of energy that goes into 
c                          the kinetic energy of the AGN feedback
c       Check_Active_AGN : check if there are active AGNs in this halo
c       iActivate_AGN    : compute a probability of forming the AGN per unit 
c                          timestep of the simulation given the AGN duty cycle
c       Compute_Galaxy_Properties : compute the properties of AGN hosting galaxies
c
c     AGN_FEEDBACK     : perform AGN feedback 
c       AGN_UPDATE_POS   : find a new location of the AGN using potential min. 
c       AGN_UPDATE_POS2  : update the positions of AGNs using previous positions and velocities 
c       Add_Energy       : add thermal energy to hvar(2,ic) & hvar(6,ic)
c       Window           : two-dimensional Gaussian window (for velocity) 
c       Read_AGNs ()     : read the AGN data to the log file.
c       Write_AGNs ()    : write the AGN data to the log file.
c       Find_Neighbors_and_Corners : finds neighbors and adjacent corners for icell
c
c     PREHEAT : uniform preaheating with a specified energy/particle
c
c     ========================================================================
c
c     ------------------------
      subroutine AGN_SETUP ( )
c     ------------------------
c     
c     initialize the AGN feedback parameters :
c
c     (1) activate AGNs only once on each big step
c     (2) reset AGNs if iAGNstep == istep 
c
c     Most of CPU time is now spent on computing galaxy properites 
c     (see Compute_Galaxy_Properties), rather than halo finding.
c
#     include "a_def.h"
      include 'a_hfind.h'
      include 'a_agn.h'

c     parameters for computing galaxy properties
      parameter ( nmin = 5 )       ! Min. # of particles
      parameter ( aMmin = 5.0e10 ) ! Min. halo [1/h Msun]
      parameter ( rout = 10. )     ! 10/h kpc physical
      
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) ' '
      write(iOWork,*) 'before AGN_SETUP....'
      close(iOWork)
#endif 

#ifdef COSMOLOGY
      aexpn  = b2a ( (1.d0*t) ) 
#endif

c
c...  Halo Finding
c     
      call HaloFind ( 0 )
      
c
c...  Compute AGN feedback parameters in the halo with M > Mmin
c

      iAGN = nAGN ! nAGN needs to be initialized to 0 when the code was first started

      do ih = 1 , nhalo 
         if ( nhp(ih) .gt. nmin ) then 

            amhalo = amhvir(ih) * pmmsun ! halo mass [1/h Msun]
            if ( amhalo .gt. aMmin ) then

               ! check if there are active AGNs in this halo already
               call Check_Active_AGN ( iFind, jAGN,  
     &              xh(ih), yh(ih), zh(ih), rhvir(ih) )

               if ( iAGNstep .lt. istep ) then
                 if ( iFind .eq. 0 .and. iActivate_AGN () .eq. 1 ) then               
                    write(*,*) ' Activate AGN' 
                    iAGN = iAGN + 1
                    nAGN = iAGN
                    call Initialize_AGN ( iAGN, ih, amhalo, rout)
                 endif 
               else 
                 if ( iFind .eq. 1 .and. 
     &                istepAGN(jAGN) .eq. istep ) then
                    write(*,*) ' Reset AGN '
                    call Initialize_AGN ( jAGN, ih, amhalo, rout )
                 endif 
               endif 

            endif ! M > Mmin
         endif ! nhp < nmin
      enddo ! ih

      ! set up AGN only once on each big step
      if ( iAGNstep .lt. istep ) then
         iAGNstep = istep  
         call Write_AGNs ( ) 
      endif

      write(*,*) ' '
      write(*,15)  nAGN, aexpn

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,15)  nAGN, aexpn
      write(iOWork,*) 'after AGN_SETUP....'
      close(iOWork)
#endif 
 15   format (' # of active AGNs = ',(i5,1x),'at aexpn = ',(f9.4,1x))

 100  return
      end

c
c     ----------------------------------------------------
      subroutine Initialize_AGN ( iAGN, ih, amhalo, rout )
c     ----------------------------------------------------
c
c     initialize AGN parameters
c
#     include "a_def.h"
      include 'a_hfind.h'
      include 'a_agn.h'

      iAGN_ON(iAGN)  = 1 ! turning the AGN ON
      istepAGN(iAGN) = istep
      aexpnAGN(iAGN) = aexpn 
      xAGN(iAGN) = xh(ih)
      yAGN(iAGN) = yh(ih)
      zAGN(iAGN) = zh(ih)
      call Compute_Galaxy_Properties ( 
     &     xh(ih), yh(ih), zh(ih), rout, vc, 
     &     amtot, amst, amgas, amcold )
      rhAGN(iAGN)     = rhvir(ih)
      amhAGN(iAGN)    = amhalo
      amtotAGN(iAGN)  = amtot
      amstAGN(iAGN)   = amst
      amgasAGN(iAGN)  = amgas
      amcoldAGN(iAGN) = amcold
      
c     vxAGN(iAGN) = zero      
c     vyAGN(iAGN) = zero
c     vzAGN(iAGN) = zero
c     LevAGN_old(iAGN) = iLv( xAGN(iAGN), yAGN(iAGN), zAGN(iAGN) )         
      
      call Setup_Jets ( iAGN )

      E61(iAGN)  = E61_AGN ( iAGN, amhalo, amgas, amcold )
      t7 (iAGN)  = t7_AGN  ( amtot )
      fke(iAGN)  = fke_AGN ( amtot )
      tAGN_now(iAGN) = zero ! initialize time
      eAGN_now(iAGN) = zero ! initialize energy            

      write(*,*) ' '
      write(*,10) iAGN
      write(*,11) xAGN(iAGN),yAGN(iAGN),zAGN(iAGN)
      write(*,12) amhAGN(iAGN), amtotAGN(iAGN), amstAGN(iAGN), 
     &     amgasAGN(iAGN), amcoldAGN(iAGN) 
      write(*,13) rout, vc 
      write(*,14) E61(iAGN), t7(iAGN), fke(iAGN)

 10   format (' Activate the AGN : #',i5)
 11   format (' (xc,yc,zc) = ', 3(g13.7,1x))        
 12   format (' Mh, mtot, mst, mgas, mcold = ',5(g13.5,1x))
 13   format (' Vc( r < ',(f6.2,1x),'kpc) = ',(g13.7,1x),' km s^-1')
 14   format (' E61, t7, fke = ',3(f11.5,1x))

      return
      end
c
c     -----------------------------
      subroutine Setup_Jets( iAGN )
c     -----------------------------
c
c     set up bipolar mechanical jets
c
      include 'a_tree.h'
      include 'a_agn.h'

      integer nnc
      real xc, yc, zc
      dimension iNC(27)    ! array for neighbors and corners

      drhox = 0.0
      drhoy = 0.0
      drhoz = 0.0
      sum1 = 0.0
      sum2 = 0.0
      
      idcell = iFindCell( MaxLevel, xAGN(iAGN), yAGN(iAGN), zAGN(iAGN) )
      call Ps ( idcell, xc, yc, zc )
      LevAGN = iLv(idcell)
      Size = CellSize(LevAGN)
c      write(*,*) idcell, xc, yc, zc, LevAGN, Size
      call Find_Neighbors_and_Corners ( idcell ,
     &     LevAGN, Size , iNC , nnc   )
      
      do ii = 1 , nnc
         iAcc = iNC(ii)
         call Ps ( iAcc, Posx, Posy, Posz )
         ix = (Posx-xc)/Size
         iy = (Posy-yc)/Size
         iz = (Posz-zc)/Size
         ! write(*,*) ii, iAcc, ix, iy, iz
      
         if ( ix .ne. 0 .and. iy .eq. 0 .and. iz .eq. 0 ) 
     &       drhox = drhox + hvar(1,iAcc)
         if ( ix .eq. 0 .and. iy .ne. 0 .and. iz .eq. 0 )
     &       drhoy = drhoy + hvar(1,iAcc)
         if ( ix .eq. 0 .and. iy .eq. 0 .and. iz .ne. 0 ) 
     &       drhoz = drhoz + hvar(1,iAcc)
         
         if ( drhox .lt. drhoy ) then
            if ( drhox .lt. drhoz ) then
               iJet(iAGN) = 1
            else 
               iJet(iAGN) = 3
            endif
         else
            if ( drhoy .lt. drhoz ) then
               iJet(iAGN) = 2
            else 
               iJet(iAGN) = 3
            endif
         endif
      
      enddo             
      
c     
c...  Determine the 2D plane from which jets are emanating
c     
      do ii = 1 , nnc
         iAcc = iNC(ii)
         call Ps ( iAcc, Posx, Posy, Posz )
         ix = (Posx-xc) / Size
         iy = (Posy-yc) / Size
         iz = (Posz-zc) / Size
         ! write(*,*) ii, iAcc, ix, iy, iz
      
         if ( iJet(iAGN) .eq. 1 ) then
            if ( ix .gt. 0 ) sum1 = sum1 + hvar(1,iAcc)
            if ( ix .lt. 0 ) sum2 = sum2 + hvar(1,iAcc)
         endif
      
         if ( iJet(iAGN) .eq. 2 ) then
            if ( iy .gt. 0 ) sum1 = sum1 + hvar(1,iAcc)
            if ( iy .lt. 0 ) sum2 = sum2 + hvar(1,iAcc)
         endif
      
         
         if ( iJet(iAGN) .eq. 3 ) then
            if ( iz .gt. 0 ) sum1 = sum1 + hvar(1,iAcc)
            if ( iz .lt. 0 ) sum2 = sum2 + hvar(1,iAcc)
         endif
      
      enddo             
      
      if ( sum1 .lt. sum2 ) iJet(iAGN) = iJet(iAGN) + 3
c      write(*,*) drhox, drhoy, drhoz, sum1, sum2
c      write(*,*) 'iJet(iAGN) =',iJet(iAGN)
      
      return
      end
c
c     -----------------------------------------------------
      real*8 function E61_AGN ( id, amhalo, amgas, amcold )
c     -----------------------------------------------------
c
c     returns the total energy of the AGN feedback 
c     in units of 10^61 erg for a given cold gas mass 
c
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_agn.h'

      integer  i, LevAGN, idcell, nnc
      real*8   Size
      dimension iNC(27)  ! array for neighbors and corners
      parameter ( eps_acc = 0.01, eps_rad = 0.1 ) 

c      E61_AGN = 1.787e-3 * ( eps_acc / 0.01 ) * 
c     &     (eps_rad / 0.1 ) * ( amcold / 1.0e8 ) 

c     Run 1 & 2 : the original simple model 
c      E61_AGN = 0.1 * ( amhalo / 1.0e13 ) 

c     Run 3 : the 2nd model (6/3/04) with Mgas(<10/h kpc comoving)
c      E61_AGN = 0.1 * ( amgas / 1.0e10 )

c     Run 4 : the 3rd model (6/21/04) with Mgas(<10/h kpc physical)
      E61_AGN = 0.005* ( amgas / 1.0e10 )

c     Run 5 : the 4th model (7/8/04) with Mgas(<10/h kpc physical)
      E61_AGN = 0.005* ( amgas / 1.0e11 )
                                                                                                                
c...  Limit E61 so that vjet < vCFL
      idcell = iFindCell( MaxLevel, xAGN(id), yAGN(id), zAGN(id) ) 
      call Ps ( idcell, xc, yc, zc )
      LevAGN  = iLv( idcell )
      Size = CellSize(LevAGN)
      anow = aexp(LevAGN)
      ameff = 0.
      write(*,*) MaxLevel, id, xAGN(id), yAGN(id), zAGN(id), idcell
      write(*,*) LevAGN, Size, anow
      call Find_Neighbors_and_Corners ( idcell ,
     &        LevAGN, Size , iNC , nnc   )
      wi2max = 0.0
      do ii = 1 , nnc
         iAcc  = iNC(ii)
         call Ps ( iAcc, Posx, Posy, Posz )
         ix = (Posx-xc)/Size
         iy = (Posy-yc)/Size
         iz = (Posz-zc)/Size
         wi = Window( iJet(id), ix, iy, iz )
         wi2 = wi**2
         wi2max = max ( wi2, wi2max )
         amcell = hvar(1,iAcc) * Size**3 * wi2
         ameff = ameff + amcell
c         write(*,*) iAcc, ix, iy, iz, wi, ameff 
      enddo 
      
c     Model 1 : (Daisuke : sometime in March 2004) -> Limit too much! 
c      vmax = 1.2e4 * anow
c      if ( anow .le. 0.417 ) vmax = 5000.

c     Model 2 : last updated : 6/3/04 
c      vmax = 30000. (Daisuke : 4/25/04) -> Try to overshoot!  
      vmax = 20000. 

c     Limit the max. energy to inject
      E61max = 2.761e-7 * ameff / hubble
     &     * Om0 * r0**3 * vmax**2 / wi2max
c      E61max = 115.466 * ameff * anow**2 / 
c     &     hubble * Om0 * r0**3 / wi2max
      write(*,*) ' ameff, E61max = ',ameff, E61max 
      if ( E61_AGN .ge. E61max ) then
         E61_AGN = E61max
      endif

c     Save the maximum velocity [km/s]
      vmaxAGN(id) = sqrt( E61_AGN / 2.761e-7 / ameff * hubble
     &     / Om0 / r0**3 * wi2max )

      return
      end
c
c     --------------------------------
      real*8 function t7_AGN ( amass )
c     --------------------------------
c
c     returns the duty cycle of AGN in units of 10^7 yrs.
c     set t7_AGN very small for instantaneous energy injection
c

      t7_AGN = 5.0
      t7_AGN = 0.0001

      return
      end
c
c     ---------------------------------
      real*8 function fke_AGN ( amass )
c     ---------------------------------
c
c     returns the fraction of energy that goes into
c     the kinetic energy of the AGN feedback.
c

      fke_AGN = 1.0

      return
      end
c
c     -------------------------------------------------------------
      subroutine Check_Active_AGN ( iFind, jAGN, xc, yc, zc, rvir ) 
c     -------------------------------------------------------------
c
c     Check if there are active AGNs in this halo
c
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_agn.h'

      ! Initialize iFind and jAGN
      iFind = 0
      jAGN  = 0

      ! Check if there are active AGNs in this halo
      DO i = 1, nAGN 
         Posx = xAGN(i)
         Posy = yAGN(i)
         Posz = zAGN(i)
         call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz )
         rr = sqrt( (Posx-xc)**2 + (Posy-yc)**2 + (Posz-zc)**2 )
         if ( rr .le. rvir ) then
            if ( iAGN_ON(i) .eq. 1 ) then
               iFind = 1
               jAGN  = i 
            endif
         endif
      ENDDO

      return
      end
c
c     ----------------------------------
      integer function iActivate_AGN ( ) 
c     ----------------------------------
c
c     Activate AGN in this halo?
c
c     Compute a probability of forming the AGN 
c     per unit timestep of the simulation given
c     the AGN duty cycle of 1Gyr.
c
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_agn.h'

c      parameter ( tcycle = 0.0001 ) 
      parameter ( tcycle = 1.0 ) 

      iActivate_AGN = 0

      tsim = dtl(MinLevel) * t0 * aexp(MinLevel)**2 / 1.0e9  ! timestep in Gyr
      prob = tsim / tcycle
      dice = RANDd ( mrand )
      if ( dice .lt. prob ) then
         iActivate_AGN = 1          
      endif

      return
      end

c
c     --------------------------------------------------
      subroutine Compute_Galaxy_Properties ( xc, yc, zc, 
     &     rc, vc, amtot, amst, amgas, amcold )
c     --------------------------------------------------
c
c     input : xc, yc, zc, rr
c     output: circular velocity [km s^-1] and mass [Msun] 
c     within r < rr [kpc]
c
c     Look at L >= 5
c
      include 'a_tree.h'
      include 'a_control.h'

      double precision xx, yy, zz
      parameter ( Tcold = 10000.0 ) 
c
      amtot = 0.0
      amst  = 0.0
      amgas = 0.0
      ampar = 0.0
      amcold = 0.0
c
      rout = rc / 1000.                   ! convert radius from kpc to Mpc
c      rout = rout / r0 / aexpn * hubble   ! radius in grid units for rout [Mpc physical]
      rout = rout / r0 / aexpn            ! radius in grid units for rout [1/h Mpc physical]
c      rout = rout / r0                    ! radius in grid units for rout [1/h Mpc comoving]

      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0
c
c...  compute mass of gas and cold gas (T<Tcold)
c

      write(*,*) ' Computing mass of gas and cold gas'
      call Get_MaxLevelNow ()
      DO Level = MaxLevelNow-5 , MaxLevelNow
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
     &                     (Posy-yc)**2 +
     &                     (Posz-zc)**2 )
                if ( rr .lt. rout ) then
                  cellmass = hvar(1,ic) * CellVolume
                  amgas = amgas + cellmass
                  Tcell = gamma1 * hvar(8,ic) / hvar(1,ic)
                  Tcell = T_0 * Tcell * a2i
                  if ( Tcell .lt. Tcold ) then
                    amcold = amcold + cellmass
c                    write(*,*) ic, cellmass, amgas, amcold
                  endif
                endif
              endif
            enddo
          enddo
        ENDIF
      ENDDO

c
c...  compute mass of DM + Stars
c
      nstars = nsp(nspecies,2) - nsp(nspecies,1) + 1
      istar1 = nsp(nspecies,1)

      write(*,*) ' Computing DM+Stars'
      do i = 1 , np
        xx = x(i)
        yy = y(i)
        zz = z(i)
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
        ! compute DM and Star profiles
        rr = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        if ( rr .le. rout ) then
          ampar = ampar + pw(i)
          if ( i .ge. istar1 ) then ! stellar particle
             amst = amst + pw(i)
          endif
        endif
      enddo

      amtot  = ampar + amgas
      amst   = amst   * aM0 / hubble   ! stellar mass [Msun]
      amgas  = amgas  * aM0 / hubble   ! gas mass [Msun]
      amcold = amcold * aM0 / hubble   ! cold gas mass [Msun]
      amtot  = amtot  * aM0 / hubble   ! total mass [Msun]
c      r50 = 0.050 / r0 / aexpn * hubble  ! 50/h kpc physical in grid units
      r50 = 0.050 / r0                    ! 50/h kpc comoving in grid units
      vc  = 293.302 * ( amtot / 1.0e12 )**(0.5) * ( rout / r50 )**(-0.5)  
c      write(*,*) amtot, amst, amgas, amcold, vc

      return
      end
c
c     --------------------------------------------------------
      subroutine AGN_UPDATE_POS( ifindcenter, xx, yy, zz, rr )
c     --------------------------------------------------------
c
c     find a new location of the AGN within r < rr 
c     centered on xx, yy, zz.  Look for only Level >= 6.
c
c     input:  
c
c      integer ifindcenter (see below) 
c                  max. gas density (ifindcenter=0)   
c                  max. DM density      (=1)
c                  max. total density   (=2)
c                  min. grav. potential (=3)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer ifindcenter, Level, icell
      double precision xx, yy, zz
      real*8 rr

      xc = xx
      yc = yy
      zc = zz
      rfind = rr

c
c...  identify the cell with the highest density
c
      rhogmax  = -1.e6  ! max. gas density
      rhodmmax = -1.e6  ! max. DM density
      rhotmax  = -1.e6  ! max. total density
      phimin   = 1.e6   ! minimum of potential
      rhodmmax = -1.e6 
      irhogmax = -1
      irhotmax = -1 
      iphimin  = -1 
      irhodmmax = -1

c
c...  Search for the cell with L >= 6
c
      call Get_MaxLevelNow ()

      DO Level = MinLevel , MaxLevelNow  
        IF ( Level .lt. 6 ) THEN
        ELSE 
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel
	    icell = iSelect(ic1)
            do ic2 = 0 , 7
	      ic = icell + ic2
              if ( iOctCh(ic) .eq. nil ) then 
                call Ps ( ic , Posx , Posy , Posz )
                call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz )
                rr = sqrt( (Posx-xc)**2 +
     &                     (Posy-yc)**2 +
     &                     (Posz-zc)**2 )
                if ( rr .le. rfind ) then
                  dcell1 = var(1,ic)*2.0**(3.*Level) + 1.0 ! tot. density
                  dcell2 = var(2,ic)  ! potential 
                  dcell3 = hvar(1,ic) ! gas density
                  dcell4 = ref(ic) ! mass of DM + stars 
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
                endif
              endif
	    enddo
          enddo
        ENDIF
      ENDDO      

      if ( ifindcenter .eq. 0 ) then 
        if ( irhogmax .eq. -1 ) then 
          write(*,*) 'error in Output_Profiles:'
          write(*,*) '   failed to find gas density maximum'
          write(*,*) '   will use input center '
        else
          call Ps ( irhogmax , xc1 , yc1 , zc1 )
c          write(*,*) 'input center x, y, z:', xc, yc, zc 
c          write(*,*) 'new center   x, y, z:', xc1,yc1,zc1
          xc = xc1
          yc = yc1 
          zc = zc1
        endif
      endif
      if ( ifindcenter .eq. 1 ) then
        if ( irhodmmax .eq. -1 ) then
          write(*,*) 'error in Output_Profiles:'
          write(*,*) '   failed to find DM density maximum'
          write(*,*) '   will use input center '
        else
          call Ps ( irhodmmax , xc1 , yc1 , zc1 )
c          write(*,*) 'input center x, y, z:', xc, yc, zc
c          write(*,*) 'new center   x, y, z:', xc1,yc1,zc1
          xc = xc1
          yc = yc1
          zc = zc1
        endif
      endif
      if ( ifindcenter .eq. 2 ) then  ! total density
        if ( irhotmax .eq. -1 ) then 
          write(*,*) 'error in Output_Profiles:'
          write(*,*) '   failed to find total density maximum'
          write(*,*) '   will use input center '
        else
          call Ps ( irhotmax , xc1 , yc1 , zc1 )
c          write(*,*) 'input center x, y, z:', xc, yc, zc 
c          write(*,*) 'new center   x, y, z:', xc1,yc1,zc1
          xc = xc1
          yc = yc1 
          zc = zc1
        endif
      endif
      if ( ifindcenter .eq. 3 ) then  ! potential 
        if ( iphimin .eq. -1 ) then 
          write(*,*) 'error in Output_Profiles:'
          write(*,*) '   failed to find potential minimum'
          write(*,*) '   will use input center '
        else
          call Ps ( iphimin , xc1 , yc1 , zc1 )
c          write(*,*) 'input center x, y, z:', xc, yc, zc 
c          write(*,*) 'new center   x, y, z:', xc1,yc1,zc1
          xc = xc1
          yc = yc1 
          zc = zc1
        endif
      endif

      xx = xc
      yy = yc
      zz = zc
      rr = rfind

      return
      end
c
c     -----------------------------------
      subroutine AGN_UPDATE_POS2 ( Level )
c     -----------------------------------
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_agn.h'

c
c     update the positions of the AGNs based on 
c     their previous positions and velocities with
c     the step size used to advance the hydro variables
c
      integer  i, LevAGN, idcell
      real*8   dx, dy, dz

      ! new positions of the AGNs
      DO i = 1, nAGN 
         if ( Level .eq. LevAGN_old(i) ) then
          dtl_step = dtl(LevAGN_old(i))
          dx = vxAGN(i) * dtl_step
          dy = vyAGN(i) * dtl_step
          dz = vzAGN(i) * dtl_step
          xAGN(i) = xAGN(i) + dx
          yAGN(i) = yAGN(i) + dy
          zAGN(i) = zAGN(i) + dz
        endif
      ENDDO

      return
      end

c     --------------------------------------
      subroutine AGN_FEEDBACK ( Level )
c     --------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_agn.h'
c
c     AGN feedback via thermal energy injection and bipolar mechanical jets 
c     with energy E61 [10^61 erg] for the period of t7 [10^7 yrs].  
c     The fraction fke of energy goes into the mechanical feedback of 18 cells
c     centered on the cell with the highest density in such a way that the
c     the resulting velocity is approximately characterized by 2D Gaussian.
c
c     In each energy injection step, do not allow gas to be AGN-heated
c     above T_ceil [K].
c
c     Last updated : 12/10/03 (D. Nagai)
c

      integer  i, LevAGN, idcell, nnc
c      real*8   Tcell
c      real*8   einj, einj_cgs, dpjet, Size
      dimension iNC(27)                 ! array for neighbors and corners

c
c...  AGN energy injection
c

      dt7yrs = dtl(8) * t0 * aexp(8)**2 / (1.0e7)
c      write(*,*) ' HELLOOO dt7yrs = ',dt7yrs


#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) ' '
      write(iOWork,*) 'before AGN_Feedback...' 
      close(iOWork)
#endif 

      DO i = 1, nAGN 

       IF ( iAGN_ON(i)  .eq. 1 ) then
        IF ( tAGN_now(i) .le. t7(i) ) then

          ! Update positions of AGN 
          call AGN_UPDATE_POS( 3, xAGN(i), yAGN(i), zAGN(i), rhAGN(i) )
          idcell = iFindCell( MaxLevel, xAGN(i), yAGN(i), zAGN(i) ) 
          call Ps ( idcell, xc, yc, zc )
          LevAGN  = iLv( idcell )
          Size = CellSize(LevAGN)
          eps  = 0.1 * Size
 
          if ( Level .eq. LevAGN ) then

            anow = aexp(Level)
            a5 = anow**5
            a3 = anow**3
            a2 = anow**2
            ai = 1.0 / anow
            a2i = 1.0 / a2
            a3i = 1.0 / a3

c...        compute the total energy to inject in this step
            dt7yrs = dtl(Level) * t0 * anow**2 / (1.0e7)
            tAGN_now(i) = tAGN_now(i) + dt7yrs
            if ( tAGN_now(i) .le. t7(i) ) then 
               deinj = E61(i) / t7(i) * dt7yrs
            else
               deinj = E61(i) / t7(i) * ( t7(i)-tAGN_now(i)+dt7yrs )  
            endif
            
c...        inject momentum and energy in neighbors and corners          
            call Find_Neighbors_and_Corners ( idcell ,
     &           LevAGN, Size , iNC , nnc   )
         
c...        compute the effective mass
            ameff = 0.
            do ii = 1 , nnc
              iAcc  = iNC(ii)
              call Ps ( iAcc, Posx, Posy, Posz )
              ix = (Posx-xc)/Size
              iy = (Posy-yc)/Size
              iz = (Posz-zc)/Size
              ! write(*,*) idcell, ii, iAcc, ix, iy, iz
              wi = Window( iJet(i), ix, iy, iz )
              wi2 = wi**2
              amcell = hvar(1,iAcc) * Size**3 * wi2
              ameff = ameff + amcell
            enddo 

c...        compute the energy and momentum to inject 
            vjet_max = 0.0 
            do ii = 1 , nnc

              iAcc  = iNC(ii)
              call Ps ( iAcc, Posx, Posy, Posz )
              ix = (Posx-xc)/Size
              iy = (Posy-yc)/Size
              iz = (Posz-zc)/Size
              wi = Window( iJet(i), ix, iy, iz )
              wi2 = wi**2
              
c...          compute the total energy to inject
              einj = 724.64 * deinj / ameff * a2 * hubble / 
     &             r0**5 / Om0**2 * hvar(1,iAcc) * wi2
              einj_cgs = 1.3802e-3 * einj * r0**5 / hubble * Om0**2 *
     &             2.0**(-3*Level) * a2i
              eAGN_now(i)  = eAGN_now(i) + einj_cgs
              
c...          compute the momentum of jets
              dpjet = 38.066 * anow / Om0 * hvar(1,iAcc) * wi *
     &             ( deinj * fke(i) * hubble / ameff / r0**5 )** 0.5        
              vjet_cgs = v0 / anow * dpjet / hvar(1,iAcc)
              vjet_est = 1903.1 * ( deinj * fke(i) * hubble / 
     &             ameff / Om0 / r0**3 )** 0.5 * wi
              vjet_max = max ( vjet_max, vjet_cgs )
              
              einj_KE = einj_KE + 2.761e-7 * vjet_cgs**2 * 
     &             hvar(1,iAcc) * 2.0**(-3*Level) / hubble * Om0 * r0**3
              einj_KE2 = einj_KE2 + fke(i) * einj_cgs
              write(*,*) ' vjet, eKE (cgs,est) = ',vjet_cgs, vjet_est,
     &             einj_KE,einj_KE2

c              e_kin = 0.5d0 * ( hvar(3,iAcc)**2 +
c     &                          hvar(4,iAcc)**2 +
c     &                          hvar(5,iAcc)**2   ) / hvar(1,iAcc)
c              eth1 = hvar(8,iAcc)
c              eth2 = (dble(hvar(2,iAcc)) - e_kin)
 
c...          inject energy and momentum
              if ( iJet(i) .eq. 1 .or. iJet(i) .eq. 4 ) then                
                if ( ix .gt. 0 ) hvar(3,iAcc) = hvar(3,iAcc) + dpjet
                if ( ix .lt. 0 ) hvar(3,iAcc) = hvar(3,iAcc) - dpjet
                if ( ix .gt. 0 .or. ix .lt. 0 ) then
                   call Add_Energy ( iAcc, fke(i), einj )
                endif
              endif
              if ( iJet(i) .eq. 2 .or. iJet(i) .eq. 5 ) then
                if ( iy .gt. 0 ) hvar(4,iAcc) = hvar(4,iAcc) + dpjet
                if ( iy .lt. 0 ) hvar(4,iAcc) = hvar(4,iAcc) - dpjet
                if ( iy .gt. 0 .or. iy .lt. 0 ) then
                   call Add_Energy ( iAcc, fke(i), einj )
                endif
              endif
              if ( iJet(i) .eq. 3 .or. iJet(i) .eq. 6 ) then
                if ( iz .gt. 0 ) hvar(5,iAcc) = hvar(5,iAcc) + dpjet
                if ( iz .lt. 0 ) hvar(5,iAcc) = hvar(5,iAcc) - dpjet
                if ( iz .gt. 0 .or. iz .lt. 0 ) then
                   call Add_Energy ( iAcc, fke(i), einj )
                endif
              endif
               
c              e_kin = 0.5d0 * ( hvar(3,iAcc)**2 +
c     &                          hvar(4,iAcc)**2 +
c     &                          hvar(5,iAcc)**2   ) / hvar(1,iAcc)
c              eth3 = hvar(8,iAcc)
c              eth4 = (dble(hvar(2,iAcc)) - e_kin)
c              write(*,*) ii, iAcc, eth1, eth2, eth3, eth4

            enddo 

            ! Save Level and velocities of the AGN cell
c            LevAGN_old(i) = LevAGN
c            vxAGN(i) = hvar(3,idcell) / hvar(1,idcell)
c            vyAGN(i) = hvar(4,idcell) / hvar(1,idcell)
c            vzAGN(i) = hvar(5,idcell) / hvar(1,idcell)
  

            write(*,*) ' '
            write(*,14) i, anow, dt7yrs
            write(*,15) LevAGN, xAGN(i), yAGN(i), zAGN(i)
            write(*,16) E61(i), t7(i), fke(i), eAGN_now(i), tAGN_now(i)
            do ii = 1 , nnc
              iAcc  = iNC(ii)
              Tcell = (hvar(7,iAcc)-1.0) * hvar(8,iAcc) / 
     &          hvar(1,iAcc) * T_0 * a2i
              write(*,*) ii, iAcc, hvar(1,iAcc), hvar(8,iAcc), Tcell,
     &             hvar(3,iAcc),hvar(4,iAcc),hvar(5,iAcc)
            enddo

          endif ! Level = LevAGN 

        ELSE
   
           ! Update positions of AGN 
           call AGN_UPDATE_POS( 3, xAGN(i), yAGN(i), zAGN(i), rhAGN(i) )
           idcell = iFindCell( MaxLevel, xAGN(i), yAGN(i), zAGN(i) ) 
           call Ps ( idcell, xc, yc, zc )
           LevAGN  = iLv( idcell )
           Size = CellSize(LevAGN)
           eps  = 0.1 * Size  
           anow = aexp(Level)
           a2 = anow**2
           a2i = 1.0 / a2
           
           write(*,*) ' '
           write(*,14) i, anow, dt7yrs
           write(*,15) LevAGN, xAGN(i), yAGN(i), zAGN(i)
           write(*,16) E61(i), t7(i), fke(i), eAGN_now(i), tAGN_now(i)

           call Find_Neighbors_and_Corners ( idcell ,
     &          LevAGN, Size , iNC , nnc )
 
           do ii = 1 , nnc
              iAcc  = iNC(ii)
              Tcell = (hvar(7,iAcc)-1.0) * hvar(8,iAcc) / 
     &          hvar(1,iAcc) * T_0 * a2i
              write(*,*) ii, iAcc, hvar(1,iAcc), hvar(8,iAcc), Tcell, 
     &             hvar(3,iAcc),hvar(4,iAcc),hvar(5,iAcc)
           enddo

           iAGN_ON(i) = 0
           write(*,*) ' '
           write(*,*) ' AGN #',i,' is turned off at t =',
     &          tAGN_now(i)*1.0e7,' yrs'
           
        ENDIF  ! if ( tAGN_now(i) .le. t7(i) )
         
       ENDIF ! if (  iAGN_ON(i) == 1 )
 
      ENDDO

 14   format (' iAGN, aexpn, dt7yrs = ',i5,1x,2(f12.5,1x))
 15   format (' LevAGN, xAGN(i), yAGN(i), zAGN(i) = ',
     &     (i3,1x),3(g12.6,1x))
 16   format (' E61, t7, fke, eAGN_now, tAGN_now = ',5(f12.6,1x)) 
 17   format (' Rhog, Tcell, vjet_max, = ',3(g13.7,1x))

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'after AGN_Feedback...' 
      close(iOWork)
#endif 
       
      return
      end

c     -------------------------------------------
      subroutine Add_Energy ( idcell, fke, einj )
c     -------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
c
      integer  idcell
      real*8   fke
      
      hvar(2,idcell) = hvar(2,idcell) + einj
      if ( fke .lt. 1.00 ) then
         hvar(8,idcell) = hvar(8,idcell) + einj * ( 1.0 - fke )
         hvar(6,idcell) = dble((hvar(7,idcell)-1.0) *
     &        hvar(8,idcell))
c         write(*,*) einj, fke, hvar(8,idcell), hvar(6,idcell)
      endif

      return
      end      
c
c     --------------------------------------------------------
      real function Window ( iaxis, ix, iy, iz )
c     --------------------------------------------------------
c     
c     Return the weights of each cell based on the 2D 
c     Gaussian normalized to the unity.
c
c     Two dimensional Gaussian :
c       sigma = 0.5  -> f1 = 0.46866, f2 = 0.10800, f3 = 0.02484
c       sigma = 1.0  -> f1 = 0.19538, f2 = 0.12333, f3 = 0.05843
c       sigma = 1.5  -> f1 = 0.14629, f2 = 0.11810, f3 = 0.09533
c       sigma = 2.0  -> f1 = 0.13037, f2 = 0.11535, f3 = 0.10206
c
      integer  iaxis, ix, iy, iz, ii, jj, kk
      real*8   dpjet, einj
      real     f1, f2, f3, frac
c      parameter ( f1 = 0.46866, f2 = 0.10800 , f3 = 0.02484 )
c      parameter ( f1 = 0.19538, f2 =0.12333 , f3 = 0.077825 )
      parameter ( f1 = 0.14629, f2 = 0.11810, f3 = 0.09533 )
c      parameter ( f1 = 0.13037, f2 = 0.11535, f3 = 0.10206 )

      ii = ix
      jj = iy
      kk = iz
      ipos = -1

      if ( iaxis .gt. 3 ) then
         iaxis = iaxis - 3
      endif

      if ( ix .lt. 0 ) ii = -1 * ii
      if ( iy .lt. 0 ) jj = -1 * jj
      if ( iz .lt. 0 ) kk = -1 * kk

      if ( iaxis .eq. 1 .and. ix .gt. 0 ) ipos = jj + kk
      if ( iaxis .eq. 2 .and. iy .gt. 0 ) ipos = ii + kk
      if ( iaxis .eq. 3 .and. iz .gt. 0 ) ipos = ii + jj
      if ( iaxis .eq. 1 .and. ix .lt. 0 ) ipos = jj + kk
      if ( iaxis .eq. 2 .and. iy .lt. 0 ) ipos = ii + kk
      if ( iaxis .eq. 3 .and. iz .lt. 0 ) ipos = ii + jj
      
      if ( ipos .eq. -1 ) Window = 0.0 
      if ( ipos .eq. 0 ) Window = f1 / 2.0
      if ( ipos .eq. 1 ) Window = f2 / 2.0
      if ( ipos .eq. 2 ) Window = f3 / 2.0
      ! write(*,*) ix, iy, iz, ii, jj, ipos, Window

      return
      end     
c
c     --------------------------------------------------------
      real function Window2 ( iaxis, ix, iy, iz )
c     --------------------------------------------------------
c     
c     Return the weights of each cell based on the 2D 
c     Gaussian normalized to the unity.
c
c     Two dimensional Gaussian :
c       sigma = 0.5  -> f1 = 0.46866, f2 = 0.10800, f3 = 0.02484
c       sigma = 1.0  -> f1 = 0.19538, f2 = 0.12333, f3 = 0.05843
c       sigma = 1.5  -> f1 = 0.14629, f2 = 0.11810, f3 = 0.09533
c       sigma = 2.0  -> f1 = 0.13037, f2 = 0.11535, f3 = 0.10206
c
      integer  iaxis, ix, iy, iz, ii, jj, kk
      real*8   dpjet, einj
      real     f1, f2, f3, frac
c      parameter ( f1 = 0.46866, f2 = 0.10800 , f3 = 0.02484 )
c      parameter ( f1 = 0.19538, f2 =0.12333 , f3 = 0.077825 )
      parameter ( f1 = 0.14629, f2 = 0.11810, f3 = 0.09533 )
c      parameter ( f1 = 0.13037, f2 = 0.11535, f3 = 0.10206 )

      ii = ix
      jj = iy
      kk = iz
      ipos = -1

      if ( ix .lt. 0 ) ii = -1 * ii
      if ( iy .lt. 0 ) jj = -1 * jj
      if ( iz .lt. 0 ) kk = -1 * kk

      if ( iaxis .eq. 1 .and. ix .ge. 0 ) ipos = jj + kk
      if ( iaxis .eq. 2 .and. iy .ge. 0 ) ipos = ii + kk
      if ( iaxis .eq. 3 .and. iz .ge. 0 ) ipos = ii + jj
      if ( iaxis .eq. 4 .and. ix .le. 0 ) ipos = jj + kk
      if ( iaxis .eq. 5 .and. iy .le. 0 ) ipos = ii + kk
      if ( iaxis .eq. 6 .and. iz .le. 0 ) ipos = ii + jj
      
      if ( ipos .eq. -1 ) Window2 = 0.0 
      if ( ipos .eq. 0 ) Window2 = f1 / 2.0
      if ( ipos .eq. 1 ) Window2 = f2 / 2.0
      if ( ipos .eq. 2 ) Window2 = f3 / 2.0
      ! write(*,*) ix, iy, iz, ii, jj, ipos, Window2

      return
      end     
c
c     -----------------------
      subroutine Read_AGNs ()
c     -----------------------
c     
c     Read the AGN data from the log file.
c
      include 'a_setup.h'
      include 'a_control.h'
      include 'a_agn.h'
      character*256 textline
c
      open  ( 12 , file = agnfile , form = 'formatted' )

c....   read header
        read(12,*) nAGN, iAGNstep
        read(12,*) textline
  
c...    read data
        DO i = 1, nAGN 
           read(12,*)
     &          iAGN, iAGN_ON(i), istepAGN(i), iJet(i),
     &          aexpnAGN(i), xAGN(i), yAGN(i), zAGN(i),
     &          E61(i), t7(i), fke(i), 
     &          tAGN_now(i), eAGN_now(i),
     &          rhAGN(i), amhAGN(i), amtotAGN(i),
     &          amstAGN(i), amgasAGN(i), amcoldAGN(i),
     &          vmaxAGN(i)
        ENDDO

      close(12)

      return
      end
c
c     ------------------------
      subroutine Write_AGNs ()
c     ------------------------
c     
c     Write the AGN data to the log file.
c
      include 'a_setup.h'
      include 'a_control.h'
      include 'a_agn.h'
c
      call Open_ASCII_File ( 12 , agnfile , sequent )

        write(12,'(i5,3x,i5,1x)') nAGN, iAGNstep
        write(12,103)
        
        DO i = 1, nAGN 
           write(12,104)
     &          i, iAGN_ON(i), istepAGN(i), iJet(i),
     &          aexpnAGN(i), xAGN(i), yAGN(i), zAGN(i),
     &          E61(i), t7(i), fke(i), 
     &          tAGN_now(i), eAGN_now(i),
     &          rhAGN(i), amhAGN(i), amtotAGN(i),
     &          amstAGN(i), amgasAGN(i), amcoldAGN(i),
     &          vmaxAGN(i)
        ENDDO

      close(12)

 103  format(' iAGN ',1x,'ON/OFF',2x,'istep',2x,'iJet ',1x,
     &       ' aexpn      x         y         z  ',2x,
     &       '    E61      t7       fke  ',1x,
     &       '  tnow  ',1x,'  Enow ', 1x,
     &       '  rh[grid]      Mh         Mtot   ',1x,
     &       '     Mst          Mg        Mcold [/h Msun]',1x,
     &       '  vmax [km/s] ' )
 104  format (i5,4x,i2,2x,i5,4x,i2,3x,f7.4,2x,3(f8.4,2x),
     &     5(f8.5,1x),7(g11.4,1x))

      return
      end
c
c     -----------------------------------------------------------
      subroutine Find_Neighbors_and_Corners ( icell , Level , 
     &                                        Size  , iNC   , nnc )
c     -----------------------------------------------------------
c
c     purpose: finds neighbors and adjacent corners for icell - 
c              only those of the icell level are taken
c   
c     this version is for "parallel" use - some computations are manually
c     unrolled to optimize the computations and memory accesses
c
c     input  : icell - cell for which routine finds neighbors and corners
c              Level - iLv(icell); MUST be > MinLevel
c              Size  - CellSize(level)
c
c     output : iNC   - indices of neighbors and corners
c              nnc   - number of neighbors and corners (size of iNC)     
c     
c
      include 'a_tree.h'
c
      dimension iNC(27)                 ! array for neighbors and corners
      dimension iENb(8,3)               ! external neighbors for given child
      dimension iChDir(6,4)             ! which children of ext.neighb to take
      dimension iNbDir(6,3)             ! vector describing a given direction
      dimension iDir(3,3)               ! specifies the order of Directions 
      integer   i1 , i2 , i3 
      integer   icell, Level , nnc
      real      Size
c
      data iENb    /  1, 2, 1, 2, 1, 2, 1, 2,       ! "x"
     &                3, 3, 4, 4, 3, 3, 4, 4,       ! "y"
     &                5, 5, 5, 5, 6, 6, 6, 6  /     ! "z"
c
      data iChDir  /  2, 1, 3, 1, 5, 1,             
     &                4, 3, 4, 2, 6, 2,
     &                6, 5, 7, 5, 7, 3,
     &                8, 7, 8, 6, 8, 4  /
c      
      data iNbDir  / -1, 1, 0, 0, 0, 0, 
     &                0, 0,-1, 1, 0, 0,
     &                0, 0, 0, 0,-1, 1  /
c
      data iDir    /  1, 2, 3,
     &                2, 3, 1,
     &                3, 1, 2  /
c
      icn   = mod ( icell - ncell0 - 1 , nchild ) + 1 ! get child number
      iPar  = iPr(icell)           ! get parent
      iCh1  = iCh(iPar,1)
      iFlag = nil                  ! zero flag
      nnc   = nil                  ! zero counter 

      iNC(1) = iCh1
      iNC(2) = iCh1 + 1
      iNC(3) = iCh1 + 2
      iNC(4) = iCh1 + 3
      iNC(5) = iCh1 + 4
      iNC(6) = iCh1 + 5
      iNC(7) = iCh1 + 6
      iNC(8) = iCh1 + 7
            
      nnc = nnc + 8

      iENb1 = iENb(icn,1)
      iENb2 = iENb(icn,2)
      iENb3 = iENb(icn,3)
    
c....   first direction

        iDir1 = iENb1 
        iDir2 = iENb2
        iDir3 = iENb3
        iNb1  = iNb(iPar,iDir1)             ! current external neighbor of iPar
        iNb2  = iNb(iNb1,iDir2)
        iNb3  = iNb(iNb2,iDir3)
        iOC = iOctCh(iNb1)
        if ( iOC .gt. nil ) then ! if neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0
          do ic2 = 1 , 4
            nnc      = nnc + 1
            iNC(nnc) = iChsh +iChDir(iDir1,ic2)
          enddo
        endif

c....   iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * (k + 1) + 1
        iOC = iOctCh(iNb2)
        if ( iOC .gt. nil ) then ! if neighbors neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0
          i1       = - iNbDir(iDir1,1) - iNbDir(iDir2,1) + 1
          i2       = - iNbDir(iDir1,2) - iNbDir(iDir2,2) + 4
          i3       = - iNbDir(iDir1,3) - iNbDir(iDir2,3) + 1
          nCh      = i1/2 + i2 + 2 * i3 
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh -4
        endif

        if ( iFlag .eq. nil ) then
           iOC =iOctCh(iNb3)
          if ( (iOC  .gt.   nil) 
     &                        .and. 
     &         (iLv(iNb3) + 1 .eq. Level) ) then
                iFlag = 1
                i1 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,1)  
     &               - iNbDir(iDir3,1) 
                i2 = - iNbDir(iDir1,2) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,2)
                i3 = - iNbDir(iDir1,3) 
     &               - iNbDir(iDir2,3)  
     &               - iNbDir(iDir3,3)  
              nCh = (i1 + 1)/2 + i2 + 2 * i3 + 4
              nnc      = nnc + 1
              iNC(nnc) = ( iOC - 1 ) * nchild + nCh + ncell0
          endif
        endif

c....   second direction    

        iDir1 = iENb2 
        iDir2 = iENb3
        iDir3 = iENb1
        iNb1  = iNb(iPar,iDir1)             ! current external neighbor of iPar
        iNb2  = iNb(iNb1,iDir2)
        iNb3  = iNb(iNb2,iDir3)
        iOC   = iOctCh(iNb1)
        if ( iOC .gt. nil ) then ! if neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0
          do ic2 = 1 , 4
            nnc      = nnc + 1
            iNC(nnc) = iChsh +iChDir(iDir1,ic2)
          enddo
        endif

c....   iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * (k + 1) + 1
        iOC = iOctCh(iNb2)
        if ( iOC .gt. nil ) then ! if neighbors neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0 
          i1 = - iNbDir(iDir1,1) - iNbDir(iDir2,1) + 1
          i2 = - iNbDir(iDir1,2) - iNbDir(iDir2,2) + 4
          i3 = - iNbDir(iDir1,3) - iNbDir(iDir2,3)  
          nCh      = (i1 + 1)/2 + i2 + 2 * i3 
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh -1
        endif

        if ( iFlag .eq. nil ) then
          iOC =iOctCh(iNb3) 
          if ( (iOC  .gt.   nil) 
     &                        .and. 
     &         (iLv(iNb3) + 1 .eq. Level) ) then            
                iFlag = 1
                i1 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,3)  
                i2 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,3)  
                i3 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,3)  
              nCh = (i1 + 1)/2 + i2 + 2 * i3 + 4
              nnc      = nnc + 1
              iNC(nnc) = ( iOC - 1 ) * nchild + nCh + ncell0
          endif
        endif

c....   third direction
     
        iDir1 = iENb3 
        iDir2 = iENb1
        iDir3 = iENb2
        iNb1  = iNb(iPar,iDir1)             ! current external neighbor of iPar
        iNb2  = iNb(iNb1,iDir2)
        iNb3  = iNb(iNb2,iDir3)
        iOC   = iOctCh(iNb1)
        if ( iOC .gt. nil ) then ! if neighbor is split
           iChsh = ( iOC - 1 ) * nchild  + ncell0
           do ic2 = 1 , 4
            nnc      = nnc + 1
            iNC(nnc) = iChsh +iChDir(iDir1,ic2)
          enddo
        endif

c....   iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * (k + 1) + 1
        iOC = iOctCh(iNb2)
        if ( iOC .gt. nil ) then ! if neighbors neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0           
          i1       = - iNbDir(iDir1,1) - iNbDir(iDir2,1) + 1
          i2       = - iNbDir(iDir1,2) - iNbDir(iDir2,2) + 5 
          i3       = - iNbDir(iDir1,3) - iNbDir(iDir2,3) 
          nCh      = i1/2 + i2 + 2 * i3
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh -2
        endif

        if ( iFlag .eq. nil ) then
          iOC =iOctCh(iNb3) 
          if ( (iOC  .gt.   nil) 
     &                        .and. 
     &         (iLv(iNb3) + 1 .eq. Level) ) then
              iFlag = 1
                i1 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,1)  
     &               - iNbDir(iDir3,1)  
                i2 = - iNbDir(iDir1,2) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,2)  
                i3 = - iNbDir(iDir1,3) 
     &               - iNbDir(iDir2,3)  
     &               - iNbDir(iDir3,3)  
              nCh = (i1 + 1)/2 + i2 + 2 * i3 + 4
              nnc      = nnc + 1
              iNC(nnc) = ( iOC - 1 ) * nchild + nCh + ncell0
          endif
        endif

      return
      end
c
c     ----------------------------
      subroutine PREHEAT ( epart )
c     ----------------------------
      include 'a_tree.h'
      include 'a_control.h'

c
c.... preheat by (epart) keV per particle
c
      DO Level = MinLevel , MaxLevelNow
        heatconst = aexpn**2 * 63.85 / r0**2 / Om0 * epart
        write(*,*) 'heatconst = ',heatconst
        IF ( Level .eq. MinLevel ) THEN
c
C$OMP PARALLEL DO DEFAULT(SHARED)
          do ic1 = 1 , ncell0
            hvar(2,ic1) = hvar(2,ic1) + heatconst * hvar(1,ic1)
            hvar(8,ic1) = hvar(8,ic1) + heatconst * hvar(1,ic1)
          enddo
        ELSE
          call Select_Cells ( Level , nLevel )
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(ic)
          do ic1 = 1 , nLevel
            ic = iSelect(ic1)
            hvar(2,ic)   = hvar(2,ic)   + heatconst * hvar(1,ic)
            hvar(2,ic+1) = hvar(2,ic+1) + heatconst * hvar(1,ic+1)
            hvar(2,ic+2) = hvar(2,ic+2) + heatconst * hvar(1,ic+2)
            hvar(2,ic+3) = hvar(2,ic+3) + heatconst * hvar(1,ic+3)
            hvar(2,ic+4) = hvar(2,ic+4) + heatconst * hvar(1,ic+4)
            hvar(2,ic+5) = hvar(2,ic+5) + heatconst * hvar(1,ic+5)
            hvar(2,ic+6) = hvar(2,ic+6) + heatconst * hvar(1,ic+6)
            hvar(2,ic+7) = hvar(2,ic+7) + heatconst * hvar(1,ic+7)
            hvar(8,ic)   = hvar(8,ic)   + heatconst * hvar(1,ic)
            hvar(8,ic+1) = hvar(8,ic+1) + heatconst * hvar(1,ic+1)
            hvar(8,ic+2) = hvar(8,ic+2) + heatconst * hvar(1,ic+2)
            hvar(8,ic+3) = hvar(8,ic+3) + heatconst * hvar(1,ic+3)
            hvar(8,ic+4) = hvar(8,ic+4) + heatconst * hvar(1,ic+4)
            hvar(8,ic+5) = hvar(8,ic+5) + heatconst * hvar(1,ic+5)
            hvar(8,ic+6) = hvar(8,ic+6) + heatconst * hvar(1,ic+6)
            hvar(8,ic+7) = hvar(8,ic+7) + heatconst * hvar(1,ic+7)
          enddo
        ENDIF
      ENDDO

      return
      end



****************************************************************************
**** OLD AGNFEEDBACK code (distribute KE uniformly in neighboring cells ****
****************************************************************************

c     ------------------------------------
      subroutine AGN_FEEDBACK_KE ( Level )
c     ------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_agn.h'
c
c     AGN feedback via thermal energy injection and bipolar mechanical jets 
c     with energy E61 [10^61 erg] for the period of t7 [10^7 yrs].  The fraction 
c     fke of energy goes into the mechanical feedback of 18 cells centered on 
c     the cell with the highest density.  The kinetic energy is distributed on
c     each cell by 2D Gaussian.
c
c     In each energy injection step, do not allow gas to be AGN-heated
c     above T_ceil [K].
c
      integer  i, LevAGN, idcell, nnc
      real*8   Tcell
      real*8   einj, einj_cgs, dpjet, Size
      dimension iNC(27)                 ! array for neighbors and corners

c
c...  AGN energy injection
c

      dt7yrs = dtl(8) * t0 * aexp(8)**2 / (1.0e7)
      write(*,*) ' HELLOOO dt7yrs = ',dt7yrs


#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) ' '
      write(iOWork,*) 'before AGN_Feedback...' 
      close(iOWork)
#endif 

      DO i = 1, nAGN 

        ! write(*,*) tAGN_now(i),t7,Level,LevAGN,iAGN_ON(i)

       IF ( iAGN_ON(i)  .eq. 1 ) then
        IF ( tAGN_now(i) .le. t7(i) ) then

          ! Update positions of AGN 
          call AGN_UPDATE_POS( 3, xAGN(i), yAGN(i), zAGN(i), rhAGN(i) )
          idcell = iFindCell( MaxLevel, xAGN(i), yAGN(i), zAGN(i) ) 
          call Ps ( idcell, xc, yc, zc )
          LevAGN  = iLv( idcell )
          Size = CellSize(LevAGN)
          eps  = 0.1 * Size
 
          if ( Level .eq. LevAGN ) then

            anow = aexp(Level)
            a5 = anow**5
            a3 = anow**3
            a2 = anow**2
            ai = 1.0 / anow
            a2i = 1.0 / a2
            a3i = 1.0 / a3

c...        compute the energy to inject in this step
            dt7yrs = dtl(Level) * t0 * anow**2 / (1.0e7)
            tAGN_now(i) = tAGN_now(i) + dt7yrs
            if ( tAGN_now(i) .le. t7(i) ) then 
               deinj = E61(i) / t7(i) * dt7yrs
            else
               deinj = E61(i) / t7(i) * ( t7(i)-tAGN_now(i)+dt7yrs )  
            endif
            einj = 724.64 * deinj * a2 * hubble / r0**5 / 
     &           2.0**(-3*Level) / Om0**2
            
            einj_cgs = 1.3802e-3 * einj * r0**5 / hubble * Om0**2 *
     &           2.0**(-3*Level) * a2i
            eAGN_now(i)  = eAGN_now(i) + einj_cgs 

c...        inject momentum and energy in neighbors and corners          
            call Find_Neighbors_and_Corners ( idcell ,
     &           LevAGN, Size , iNC , nnc   )
          
            vjet_max = 0.0 
            do ii = 1 , nnc

              iAcc  = iNC(ii)
              call Ps ( iAcc, Posx, Posy, Posz )
              ix = (Posx-xc)/Size
              iy = (Posy-yc)/Size
              iz = (Posz-zc)/Size
              ! write(*,*) idcell, ii, iAcc, ix, iy, iz
              
c...          compute the momentum of jets
c              dpjet = 38.066 * anow / Om0 *
c     &             ( hvar(1,iAcc) * E61(i) * dt7yrs * fke(i) * hubble * 
c     &               2.0**(3*Level) / r0**5 / t7(i) )** 0.5
              dpjet = 38.066 * anow / Om0 *
     &             ( hvar(1,iAcc) * deinj * fke(i) * hubble * 
     &               2.0**(3*Level) / r0**5 )** 0.5
                              
              if ( iJet(i) .eq. 1 .or. iJet(i) .eq. 4 ) then
                 v1 = v0 / anow * hvar(3,iAcc) / hvar(1,iAcc)
              endif
              if ( iJet(i) .eq. 2 .or. iJet(i) .eq. 5 ) then
                 v1 = v0 / anow * hvar(4,iAcc) / hvar(1,iAcc)
              endif
              if ( iJet(i) .eq. 3 .or. iJet(i) .eq. 6 ) then
                 v1 = v0 / anow * hvar(5,iAcc) / hvar(1,iAcc)
              endif
              
              if ( iJet(i) .eq. 1 ) then                
                call Window_KE( dpjet, einj, iJet(i), ix, iy, iz )
                if ( ix .gt. 0 ) hvar(3,iAcc) = hvar(3,iAcc) + dpjet
                if ( ix .eq. 0 ) hvar(3,iAcc) = hvar(3,iAcc) - dpjet
                if ( ix .ge. 0 ) call Add_Energy ( iAcc, fke(i), einj )
                if ( ix .ge. 0 ) then
                   vjet_cgs = v0 / anow * dpjet / hvar(1,iAcc)
                   vjet_max = max ( vjet_max, vjet_cgs )
                endif
              endif
              if ( iJet(i) .eq. 2 ) then
                call Window_KE( dpjet, einj, iJet(i), ix, iy, iz )
                if ( iy .gt. 0 ) hvar(4,iAcc) = hvar(4,iAcc) + dpjet
                if ( iy .eq. 0 ) hvar(4,iAcc) = hvar(4,iAcc) - dpjet
                if ( iy .ge. 0 ) call Add_Energy ( iAcc, fke(i), einj )
                if ( iy .ge. 0 ) then
                   vjet_cgs = v0 / anow * dpjet / hvar(1,iAcc)
                   vjet_max = max ( vjet_max, vjet_cgs )
                endif
              endif
              if ( iJet(i) .eq. 3 ) then
                call Window_KE( dpjet, einj, iJet(i), ix, iy, iz )
                if ( iz .gt. 0 ) hvar(5,iAcc) = hvar(5,iAcc) + dpjet
                if ( iz .eq. 0 ) hvar(5,iAcc) = hvar(5,iAcc) - dpjet
                if ( iz .ge. 0 ) call Add_Energy ( iAcc, fke(i), einj )
                if ( iz .ge. 0 ) then
                   vjet_cgs = v0 / anow * dpjet / hvar(1,iAcc)
                   vjet_max = max ( vjet_max, vjet_cgs )
                endif
             endif
              if ( iJet(i) .eq. 4 ) then
                call Window_KE( dpjet, einj, iJet(i), ix, iy, iz )
                if ( ix .eq. 0 ) hvar(3,iAcc) = hvar(3,iAcc) + dpjet
                if ( ix .lt. 0 ) hvar(3,iAcc) = hvar(3,iAcc) - dpjet
                if ( ix .le. 0 ) call Add_Energy ( iAcc, fke(i), einj )
                if ( ix .le. 0 ) then
                   vjet_cgs = v0 / anow * dpjet / hvar(1,iAcc)
                   vjet_max = max ( vjet_max, vjet_cgs )
                endif
             endif
              if ( iJet(i) .eq. 5 ) then
                call Window_KE( dpjet, einj, iJet(i), ix, iy, iz )
                if ( iy .eq. 0 ) hvar(4,iAcc) = hvar(4,iAcc) + dpjet
                if ( iy .lt. 0 ) hvar(4,iAcc) = hvar(4,iAcc) - dpjet
                if ( iy .le. 0 ) call Add_Energy ( iAcc, fke(i), einj )
                if ( iy .le. 0 ) then
                   vjet_cgs = v0 / anow * dpjet / hvar(1,iAcc)
                   vjet_max = max ( vjet_max, vjet_cgs )
                endif
              endif
              if ( iJet(i) .eq. 6 ) then
                call Window_KE( dpjet, einj, iJet(i), ix, iy, iz )
                if ( iz .eq. 0 ) hvar(5,iAcc) = hvar(5,iAcc) + dpjet
                if ( iz .lt. 0 ) hvar(5,iAcc) = hvar(5,iAcc) - dpjet
                if ( iz .le. 0 ) call Add_Energy ( iAcc, fke(i), einj )
                if ( iz .le. 0 ) then
                   vjet_cgs = v0 / anow * dpjet / hvar(1,iAcc)
                   vjet_max = max ( vjet_max, vjet_cgs )
                endif
             endif

              if ( iJet(i) .eq. 1 .or. iJet(i) .eq. 4 ) then
                 v2 = v0 / anow * hvar(3,iAcc) / hvar(1,iAcc)
              endif
              if ( iJet(i) .eq. 2 .or. iJet(i) .eq. 5 ) then
                 v2 = v0 / anow * hvar(4,iAcc) / hvar(1,iAcc)
              endif
              if ( iJet(i) .eq. 3 .or. iJet(i) .eq. 6 ) then
                 v2 = v0 / anow * hvar(5,iAcc) / hvar(1,iAcc)
              endif
              
c              write(*,*) ' rhog = ', hvar(1,iAcc)
c              write(*,*) ' v1, v2, vjet =',v1, v2, v2-v1,' [km/s]'
               
            enddo 

            ! Save Level and velocities of the AGN cell
c            LevAGN_old(i) = LevAGN
c            vxAGN(i) = hvar(3,idcell) / hvar(1,idcell)
c            vyAGN(i) = hvar(4,idcell) / hvar(1,idcell)
c            vzAGN(i) = hvar(5,idcell) / hvar(1,idcell)
  
            Tcell = (hvar(7,idcell)-1.0) * hvar(8,idcell) / 
     &           hvar(1,idcell) * T_0 * a2i

            write(*,*) ' '
            write(*,14) i, anow, dt7yrs
            write(*,15) LevAGN, xAGN(i), yAGN(i), zAGN(i)
            write(*,16) E61(i), t7(i), fke(i), eAGN_now(i), tAGN_now(i)
            write(*,17) hvar(1,idcell), Tcell, vjet_max

          endif ! Level = LevAGN 

        ELSE

           iAGN_ON(i) = 0
           write(*,*) ' '
           write(*,*) ' AGN is turned off at t =',
     &          tAGN_now(i)*1.0e7,' yrs'
           
        ENDIF  ! if ( tAGN_now(i) .le. t7(i) )
         
       ENDIF ! if (  iAGN_ON(i) == 1 )
 
      ENDDO

 14   format (' iAGN, aexpn, dt7yrs = ',i5,1x,2(f12.5,1x))
 15   format (' LevAGN, xAGN(i), yAGN(i), zAGN(i) = ',
     &     (i3,1x),3(g12.6,1x))
 16   format (' E61, t7, fke, eAGN_now, tAGN_now = ',5(f12.6,1x)) 
 17   format (' Rhog, Tcell, vjet_max, = ',3(g13.7,1x))

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      write(iOWork,*) 'after AGN_Feedback...' 
      close(iOWork)
#endif 
       
      return
      end
c
c
c     ---------------------------------------------------------
      subroutine Window_KE( dpjet, einj, iaxis, ix, iy, iz )
c     ---------------------------------------------------------
c     
c     Compute the fraction of momentum and energy according to
c     the two dimensional Gaussian normalized to the unity.
c
c     Two dimensional Gaussian :
c       sigma = 1.0  -> f1 = 0.19538, f2 = 0.12333, f3 = 0.05843
c       sigma = 0.5  -> f1 = 0.46866, f2 = 0.10800, f3 = 0.02484
c
      integer  iaxis, ix, iy, iz, ii, jj, kk
      real*8   dpjet, einj
      real     f1, f2, f3, frac
c      parameter ( f1 = 0.46866, f2=0.10800 , f3=0.02484 )
      parameter ( f1 = 0.19538, f2=0.12333 , f3=0.05843 )

      ii = ix
      jj = iy
      kk = iz

      if ( iaxis .eq. 1 .or. iaxis .eq. 4 ) then
         if ( jj .lt. 0 ) jj = -1 * jj
         if ( kk .lt. 0 ) kk = -1 * kk
         ipos = jj + kk
      endif

      if ( iaxis .eq. 2 .or. iaxis .eq. 5 ) then
         if ( ii .lt. 0 ) ii = -1 * ii
         if ( kk .lt. 0 ) kk = -1 * kk
         ipos = ii + kk
      endif

      if ( iaxis .eq. 3 .or. iaxis .eq. 6 ) then
         if ( ii .lt. 0 ) ii = -1 * ii
         if ( jj .lt. 0 ) jj = -1 * jj
         ipos = ii + jj
      endif

      if ( ipos .eq. 0 ) frac = f1 / 2.0
      if ( ipos .eq. 1 ) frac = f2 / 2.0
      if ( ipos .eq. 2 ) frac = f3 / 2.0
      ! write(*,*) ix, iy, iz, ii, jj, ipos, frac

      dpjet = dpjet * frac**0.5
      einj  = einj  * frac

      return
      end     
