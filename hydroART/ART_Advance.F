c     =====================================================================
c                                                                         
c       ART Version 3-1: Advancing control variables to the next time step  
c                                                                         
c             by Andrey Kravtsov and Anatoly Klypin (1997-1998)             
c                                                                         
c     =====================================================================

c     -----------------------------
      subroutine AdvanceTime ( dt )
c     -----------------------------
c
c     advances time variables, 
c     computes kinetic and potential energies
c     and energy conservation error
c
      real*8 dt
      double precision e_gas, eth_gas, ekin_gas, rho_gas, volcell
c
      double precision b2a
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      real*8 epot_
      real xmin,xmax,ymin,ymax,zmin,zmax
      real vxmin,vxmax,vymin,vymax,vzmin,vzmax
      real dxmin,dxmax,dymin,dymax,dzmin,dzmax
      real*8 DTDA, fp(2), INTEGRATE, ae, age, dtyrs
      real*8 LevelVol(MinLevel:MaxLevel), resvol 
      external DTDA
#ifdef STARFORM
      real*8 wstot, wstoti, asmin, asmax
      real*8 ws_old, ws_oldi, dws, dwsi, sfr
      common / s_old1 / ws_old, ws_oldi
      real*8 tdum, xsmin, xsmax, ysmin, ysmax, zsmin, zsmax
      real*8 vxsmin, vxsmax, vysmin, vysmax, vzsmin, vzsmax
#endif 

c
#ifdef COSMOLOGY
      aexpn0 = b2a ( (1.d0*t) ) 
#endif
#ifndef COSMOLOGY
      aexpn0 = 1.0
#endif

      t = t + dt
      
#ifdef COSMOLOGY
      aexpn  = b2a ( (1.d0*t) ) 
      ap1    = b2a ( (1.d0*(t-5.d-1*dt)) )
      astep  = aexpn - aexpn0
      ae = aexpn
      fp(1)  = Om0
      fp(2)  = Oml0

      age = INTEGRATE(DTDA, fp, 2, 0.d0, ae, 0.1d0*ae, 1.d-8)/
     &       (100.d0 * hubble *gyr/mpc)
#endif
#ifndef COSMOLOGY
      aexpn = 1.0 
      ap1 = 1.
      astep = 0.
      age = 0.0
#endif

      istep  = istep + 1
      istep2 = istep2 + 1
c
      dtyrs = dt * t0 * aexpn**2 
c
#ifdef HYDRO
      do Level = MinLevel , MaxLevel 
        LevelVol(Level) = 0.0 
      enddo

        e_gas    = zero 
        eth_gas  = zero
        ekin_gas = zero
        rho_gas  = zero
        am_gas_tot = Omb0 / Om0 * float(ncell0) 
        do Level = MinLevel , MaxLevelNow
          if ( Level .eq. MinLevel ) then 
            do idcell = 1 , ncell0
              if ( iOctCh(idcell) .eq. 0 ) then 
                ethermal = hvar(6,idcell)/(gamma-1) 
                ekinetic = ( hvar(3,idcell)**2 +
     &                       hvar(4,idcell)**2 + 
     &                       hvar(5,idcell)**2 ) * 
     &                       onehalf / hvar(1,idcell)  
                eth_gas  = eth_gas + ethermal
                ekin_gas = ekin_gas + ekinetic 
                e_gas    = e_gas + ethermal + ekinetic
                rho_gas  = rho_gas +  hvar(1,idcell) 
              endif
            enddo
          else
            nLevel = iNOLL(Level)
            call Select_Cells ( Level , nLevel ) 
            volcell = CellSize(Level)**3
            do i1 = 1 , nLevel
              icell = iSelect(i1)
              do i2 = 0 , 7
                idcell = icell + i2
                if ( iOctCh(idcell) .eq. 0 ) then 
                  LevelVol(Level) = LevelVol(Level) + volcell 
                  ethermal = hvar(6,idcell)/(gamma-1) * 
     &                       volcell
                  ekinetic = ( hvar(3,idcell)**2 +
     &                         hvar(4,idcell)**2 + 
     &                         hvar(5,idcell)**2 ) * 
     &                         onehalf / hvar(1,idcell) *
     &                       volcell
                  eth_gas  = eth_gas + ethermal
                  ekin_gas = ekin_gas + ekinetic 
                  e_gas    = e_gas + ethermal + ekinetic
ceverino09062007                  rho_gas  = rho_gas + hvar(1,idcell) * vol
              rho_gas  = rho_gas + hvar(1,idcell) * volcell
                endif
              enddo
            enddo
          endif          
        enddo
#endif

ceverino09062007
        write(*,*) istep, 'rho_gas=',  rho_gas

#ifdef NBODY
        do ic1 = 1 , nspecies
          ekin_(ic1) = zero 
          epot_      = zero 
        enddo 
        do ic0 = 1 , nspecies
          do ic1 = nsp(ic0,1) , nsp(ic0,2)
            ekin_(ic0) = ekin_(ic0) + wpar(ic0) * 
     &                   (vx(ic1)**2 + vy(ic1)**2 + vz(ic1)**2)
            epot_ = epot_ + pw(ic1)*pdummy(ic1)
          enddo      
          write(*,'(''iSpecie ='',i3,'' wSpec ='',g13.6,
     &              '' Np ='',i8,'' Epot ='',g15.5)') 
     &        ic0,wpar(ic0),(nsp(ic0,2)-nsp(ic0,1)+1),epot_*0.5
       enddo
c
        ENkin = zero     
        do ic0 = 1 , nspec 
          ekin_(ic0) = ekin_(ic0) / 2.0 / ap1**2 
          ENkin      = ENkin + ekin_(ic0)
        enddo
        ENpot = epot_ / 2.0


#ifdef STARFORM

        if ( istep .eq. 1 ) then 
           ws_old  = 0.0 
           ws_oldi = 0.0
        endif
c
c....   output info about stellar particles
c
        simin =  1.d30
        simax = -1.d30
        smmin =  1.d30
        smmax = -1.d30
        asmin =  1.d30
        asmax = -1.d30
        xsmin =  1.d30
        xsmax = -1.d30 
        ysmin =  1.d30
        ysmax = -1.d30 
        zsmin =  1.d30
        zsmax = -1.d30 
        vxsmin =  1.d30
        vxsmax = -1.d30 
        vysmin =  1.d30
        vysmax = -1.d30 
        vzsmin =  1.d30
        vzsmax = -1.d30 
        wstot  = 0.
        wstoti = 0.
#ifdef ENRICH
        z2min =  1.d30
        z2max = -1.d30
#endif
#ifdef ENRICH_SNIa
        z1min =  1.d30
        z1max = -1.d30
#endif
        nstar1 = nsp(nspecies,1) - 1 
        do is = nsp(nspecies,1) , nsp(nspecies,2)
          istar = is - nstar1
          tdum  = t - tbirth(istar)
          asmin = min ( asmin , tdum )
          asmax = max ( asmax , tdum )
          smmin = min ( smmin , pw(is) )
          smmax = max ( smmax , pw(is) )
          wstoti = wstoti + pw0(istar)
          wstot  = wstot  + pw(is)
          simin = min ( simin , pw0(istar) )
          simax = max ( simax , pw0(istar) )
          xsmin = min ( xsmin , x(is) ) 
          xsmax = max ( xsmax , x(is) ) 
          ysmin = min ( ysmin , y(is) ) 
          ysmax = max ( ysmax , y(is) ) 
          zsmin = min ( zsmin , z(is) ) 
          zsmax = max ( zsmax , z(is) ) 
          vxsmin = min ( vxsmin , vx(is) ) 
          vxsmax = max ( vxsmax , vx(is) ) 
          vysmin = min ( vysmin , vy(is) ) 
          vysmax = max ( vysmax , vy(is) ) 
          vzsmin = min ( vzsmin , vz(is) ) 
          vzsmax = max ( vzsmax , vz(is) ) 
#ifdef ENRICH
          z2min = min ( z2min , zstII(istar) )
          z2max = max ( z2max , zstII(istar) ) 
#endif
#ifdef ENRICH_SNIa
          z1min = min ( z1min , zstIa(istar) )
          z1max = max ( z1max , zstIa(istar) ) 
#endif
c
c....   add stellar mass to gas masses to monitor conservation of mass
c
ceverino09062007          rho_gas = rho_gas + pw(is)
           rho_gas = rho_gas + dble( pw(is) )
        enddo

ceverino09062007
         write(*,*) 'bary_t=' ,  rho_gas

c
c....   compute refined volume 
c       below I assume that starformation can only occur in a subset of
c       the total volume, defined by the refined zero level cells
c       if SF can occur anywhere in the box, resvol must be set to boxh**3
c
        resvol = 0.
        do Level = MinLevel+1 , MaxLevel 
          resvol = resvol + LevelVol(Level) 
        enddo
        resvol = resvol * boxh**3 / ncell0 

#ifdef COSMOLOGY
c        ae = b2a ( asmin ) 
c        fp(1)  = Om0
c        fp(2)  = Oml0

c        asmin = INTEGRATE(DTDA, fp, 2, 0.d0, ae, 0.1d0*ae, 1.d-8)/
c     &          (100.d0 * hubble *gyr/mpc)
c        ae = b2a ( asmax )
c        asmax = INTEGRATE(DTDA, fp, 2, 0.d0, ae, 0.1d0*ae, 1.d-8)/
c     &          (100.d0 * hubble *gyr/mpc)
#endif
#ifndef COSMOLOGY
        asmin = asmin * t0
        asmax = asmax * t0
#endif
        dws  = wstot  - ws_old 
        dwsi = wstoti - ws_oldi 
        write(*,89) nsp(nspecies,2) - nsp(nspecies,1) + 1
        write(*,87) wstoti*aM0, dwsi*aM0
        write(*,88) wstot*aM0, dws*aM0
        write(*,90) simin*aM0, simax*aM0
        write(*,91) smmin*aM0, smmax*aM0
        write(*,92) asmin, asmax
#ifdef ENRICH
        write(*,93) z2min/0.0199, z2max/0.0199
#endif
#ifdef ENRICH_SNIa
        write(*,94) z1min/0.0199, z1max/0.0199
#endif
        write(*,95)  xsmin, xsmax 
        write(*,96)  ysmin, ysmax 
        write(*,97)  zsmin, zsmax 
        write(*,98)  vxsmin*v0/aexpn, vxsmax*v0/aexpn
        write(*,99)  vysmin*v0/aexpn, vysmax*v0/aexpn 
        write(*,100) vzsmin*v0/aexpn, vzsmax*v0/aexpn 
 87     format('Tot. init. mass of *, change [/hMsun] =',
     &       g13.6,1x,g13.6)
 88     format('Tot. mass in *, change [/hMsun] =',
     &       g13.6,1x,g13.6)
 89     format('Number of stellar particles    =',i8)
 90     format('min/max initial mass [/h Msun] =',2(g13.6,1x))
 91     format('min/max current mass [/h Msun] =',2(g13.6,1x))
 92     format('min/max age  [yrs]             =',2(g13.6,1x))
 93     format('min/max SNII metal. [Z/Zsun]   =',2(g13.6,1x))
 94     format('min/max SNIa metal. [Z/Zsun]   =',2(g13.6,1x))
 95     format('xmin/max [grid units]          =',2(g13.6,1x))
 96     format('ymin/max [grid units]          =',2(g13.6,1x))
 97     format('zmin/max [grid units]          =',2(g13.6,1x))
 98     format('vxmin/max [km/s]               =',2(g13.6,1x))
 99     format('vymin/max [km/s]               =',2(g13.6,1x))
 100    format('vzmin/max [km/s]               =',2(g13.6,1x))

        if ( istep .eq. 1 ) then 
          call Open_ASCII_File ( 12 , starfile , append ) 
          write(12,101)
          close ( 12 ) 
 101      format('istep',1x,'     t     ',1x,'      dt ',1x,
     &           '     aexpn  ',1x,'t [Gyrs] dt [yrs]   ',
     &           '   N*   ',
     &       '    M*         dM*        Mi*       dMi* [/h Msun]',
     &       '    SFR [Msun/yr/Mpc^3]' )
        endif
c
        if ( resvol .gt. 0.0 ) then 
          sfr = dwsi * aM0 * hubble * hubble / dtyrs / resvol 
        else
          sfr = 0.0 
        endif

        call Open_ASCII_File ( 12 , starfile , append ) 
        write(12,102) istep, t, dt, aexpn, age, dtyrs,
     &      nsp(nspecies,2) - nsp(nspecies,1) + 1, 
     &      wstot*aM0, dws*aM0, wstoti*aM0, dwsi*aM0, sfr
        close ( 12 ) 
 102    format (i5,1x,e13.6,1x,e11.6,1x,f8.6,1x,f6.3,1x,e11.6,
     &          i9,1x,5(e11.6,1x))
c
        ws_old = wstot
        ws_oldi = wstoti 
c
#endif
#endif
c
c.... energy conservation control
c
        Ekin1 = Ekin
        Ekin  = ENkin
        ncells = noct * nchild + ncell0
        IF ( istep .eq. 1 ) THEN
          au0   = aexpn0 * ENpot
          aeu0  = au0 + aexpn0 * Ekin 
          TINTG = 0.
c
          call Open_ASCII_File ( 12 , 'run.log ' , append ) 
          write(*,'(''istep'',1x,''     t     '',1x,''  dt '',1x,
     &              ''  aexpn  '',1x,'' MaxL '',1x,''  ncells '',1x,
     &              ''  Ekin_p  '',1x,''  Epot_p  '',1x,'' error '',1x,
     &              ''   rho_g  '',1x,''  eth_g  '',1x,
     &              ''  ekin_g  '',1x,''  etot_g  '')') 
          write(*,'(i5,1x,e13.6,1x,e11.6,1x,f8.6,1x,
     &              i2,1x,i8,1x,7(g13.7,2x))')
     &              istep, t, dt, aexpn, MaxLevelNow, ncells, 
     &              Enkin, ENpot,error,
     &              rho_gas, eth_gas, ekin_gas, e_gas
          write(12,'(''istep'',1x,''     t     '',1x,''  dt '',1x,
     &              ''  aexpn  '',1x,'' MaxL '',1x,''  ncells '',1x,
     &              ''  Ekin_p  '',1x,''  Epot_p  '',1x,'' error '',1x,
     &              ''   rho_g  '',1x,''  eth_g  '',1x,
     &              ''  ekin_g  '',1x,''  etot_g  '')') 
          write(12,'(i5,1x,e13.6,1x,e11.6,1x,f8.6,1x,
     &              i2,1x,i8,1x,7(g14.7,2x))')
     &              istep, t, dt, aexpn, MaxLevelNow, ncells, 
     &              Enkin, ENpot,error,
     &              rho_gas, eth_gas, ekin_gas, e_gas
         close ( 12 )
c
        ELSE
          TINTG = TINTG + 2.*(aexpn0 - ap0)*EKIN1
          Ekinp = (Ekin * (ap1-ap0-astep*0.5)+Ekin1*0.5*astep)/(ap1-ap0)
          error = (aexpn0*(Ekinp+ENpot) - aeu0 + TINTG) /
     &            (aexpn0*ENpot - au0)

          call Get_MaxLevelNow ()
          call Open_ASCII_File ( 12 , 'run.log ' , append ) 
          write(*,'(i5,1x,e13.6,1x,e11.6,1x,f8.6,1x,
     &              i2,1x,i8,1x,7(g13.7,2x))' )
     &              istep, t, dt, aexpn, MaxLevelNow, ncells, 
     &              Enkin, ENpot, error,
     &              rho_gas, eth_gas, ekin_gas, e_gas
          write(12,'(i5,1x,e13.6,1x,e11.6,1x,f8.6,1x,
     &              i2,1x,i8,1x,7(g13.7,2x))' )
     &              istep, t, dt, aexpn, MaxLevelNow, ncells, 
     &              Enkin, ENpot, error,
     &              rho_gas, eth_gas, ekin_gas, e_gas
          close ( 12 )
        ENDIF

      ap0 = ap1
c
c.... check if format conforms with the number of species
c     namely check number of e12.4
c
 40   format (i5 , 1x , f10.4 , 1x , e12.4 , e12.4 , 2(1x,e12.4) )
 50   format (i5 , 1x ,  f7.5 , 1x , i2 , 1x ,  i8 , 4(1x,E13.5) )

c
c.... write the AGN feedback parameters
c
#ifdef AGNFEEDBACK     
      call Write_AGNs ( ) 
#endif
c
c
c
#ifdef DEBUG
      call Open_ASCII_File ( iODebug , debugfile , append ) 
      call Particles_minmax ( iODebug )
      close(iODebug)
#endif

#ifdef DEBUG
        call Open_ASCII_File ( iODebug , debugfile , append ) 
        rhomin = 1.e20
        rhomax = -rhomin
        phimin = 1.e20
        phimax = -phimin
        do Level = MinLevel , MaxLevelNow
          if ( Level .eq. MinLevel ) then 
            do idcell = 1 , ncell0
                phimin = min(phimin,var(2,idcell))
                phimax = max(phimax,var(2,idcell))
                rhomin = min(rhomin,var(1,idcell))
                rhomax = max(rhomax,var(1,idcell))
            enddo
          else
            nLevel = iNOLL(Level)
            call Select_Cells ( Level , nLevel ) 
            do i1 = 1 , nLevel
              icell = iSelect(i1)
              do i2 = 0 , 7
                idcell = icell + i2
                if ( iOctCh(idcell) .eq. 0 ) then 
                  phimin = min(phimin,var(2,idcell))
                  phimax = max(phimax,var(2,idcell))
                  rhomin = min(rhomin,var(1,idcell))
                  rhomax = max(rhomax,var(1,idcell))
                endif
              enddo
            enddo
          endif          
          write(iODebug,*) 'L=',Level,' phi rho min,max=',
     &       phimin, phimax, rhomin,rhomax
        enddo
        close(iODebug)
#endif

      return
      end

c     -------------------------------------
      subroutine Particles_minmax ( iunit )
c     -------------------------------------
      include 'a_tree.h'
      include 'a_control.h'

        xmin = 1.e6
        xmax = -xmin
        ymin = 1.e6
        ymax = -xmin
        zmin = 1.e6
        zmax = -xmin
        vxmin = 1.e6
        vxmax = -xmin
        vymin = 1.e6
        vymax = -xmin
        vzmin = 1.e6
        vzmax = -xmin
        do ip = 1 , np
          if ( x(ip) .lt. xmin ) then
            xmin = x(ip)
            dxmin = ip
          endif
          if ( y(ip) .lt. ymin ) then
            ymin = y(ip)
            dymin = ip
          endif
          if ( z(ip) .lt. zmin ) then
            zmin = z(ip)
            dzmin = ip
          endif
          if ( x(ip) .gt. xmax ) then
            xmax = x(ip)
            dxmax = ip
          endif
          if ( y(ip) .gt. ymax ) then
            ymax = y(ip)
            dymax = ip
          endif
          if ( z(ip) .gt. zmax ) then
            zmax = z(ip)
            dzmax = ip
          endif
        enddo
        write(iunit,*) 'xmin,max=',xmin,xmax,dxmin,dxmax
        write(iunit,*) 'ymin,max=',ymin,ymax,dymin,dymax
        write(iunit,*) 'zmin,max=',zmin,zmax,dzmin,dzmax

        do ip = 1 , np
          if ( vx(ip) .lt. vxmin ) then
            vxmin = vx(ip)
            dxmin = ip
          endif
          if ( vy(ip) .lt. vymin ) then
            vymin = vy(ip)
            dymin = ip
          endif
          if ( vz(ip) .lt. vzmin ) then
            vzmin = vz(ip)
            dzmin = ip
          endif
          if ( vx(ip) .gt. vxmax ) then
            vxmax = vx(ip)
            dxmax = ip
          endif
          if ( vy(ip) .gt. vymax ) then
            vymax = vy(ip)
            dymax = ip
          endif
          if ( vz(ip) .gt. vzmax ) then
            vzmax = vz(ip)
            dzmax = ip
          endif
        enddo
        write(iunit,*) 'vxmin,max=',vxmin,vxmax,dxmin,dxmax
        write(iunit,*) 'vymin,max=',vymin,vymax,dymin,dymax
        write(iunit,*) 'vzmin,max=',vzmin,vzmax,dzmin,dzmax

      return
      end
