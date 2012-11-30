c
c     Set up initial conditions for Zeldovich pancake collapse test
c

c     ----------------------                                                 
      subroutine Init_Run ()                
c     ----------------------           
c
      real*8 dt
      common /initrun_COM/ xc, astart, aend,
     &                     dens0, pre0, across, alambda, gamma
      common / CENTER / x0, y0, z0
      real*8 a2b, b2a, qsolve, dgrowthdt
      real   zbrent 
      real*8 fp(5), x1, x2
      external qsolve
      real*8  v_tot
      real*8 eta1
      parameter ( eta1 = 1.d-3 )
      real*8 growth
      external growth
c
#include "a_def.h"
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
c
      dens0  = 1.0
      astart = 1.0 / 51.0 
      across = 1./6. 
      aend   = 1.0 
      alambda = float(ng)
      gamma   = 5.0/3.0 
      radin(1) = gamma
      Om0  = 1.0
      Omb0 = 0.1
      Oml0 = 0.0
      hubble = 0.5
      boxh   = 5.0 

      write (*,'(1x,''Om0           '',1pe12.4)') Om0
      write (*,'(1x,''Omb0          '',1pe12.4)') Omb0
      write (*,'(1x,''Oml0          '',1pe12.4)') Oml0
      write (*,'(1x,''hubble        '',1pe12.4)') hubble
      write (*,'(1x,''Lbox [/h Mpc] '',1pe12.4)') boxh
      write (*,'(1x,''CFL     '',1pe12.4)') cfl
      write (*,'(1x,''RHO0    '',1pe12.4)') dens0
      write (*,'(1x,''P0      '',1pe12.4)') pre0
      write (*,'(1x,''lambda  '',1pe12.4)') alambda
      write (*,'(1x,''across  '',1pe12.4)') across
      write (*,'(1x,''astart   '',1pe12.4)') astart
      write (*,'(1x,''aend    '',1pe12.4)') aend
      write (*,'(1x,''gamma   '',1pe12.4)') gamma
         tbmass(MinLevel) = 0.2
         tbmass(1) = 0.6
c
      IF ( start ) THEN 
c
         write(*,*) 'Setting up initial conditions...'
c
         aexpn = astart
         t = a2b ( (1.d0*astart) ) 
         istep = 0

         f_gas = Omb0/Om0 
         call SetUnits ()

         a_th = 1.0 / (1.e3*(Omb0*hubble**2)**0.4)
         
         if ( aexpn .lt. a_th ) then 
           TinitK = 2.726 / aexpn
         else
           TinitK = 2.726 / a_th * (a_th/aexpn)**2
         endif
         Tinit = TinitK * aexpn**2 / (0.31*T_0/0.6) ! 0.31 is for neutral gas
         pre0 = Tinit * dens0 * (Omb0/Om0)

         write(*,*) 'Tinit =',TinitK
         pause

         do L = MinLevel , MaxLevel
           aexp(L) = astart
           tl(L) = t
         enddo
c
         m_tot = 1
         xc = (1.*ng/2) + 0.5
         yc = (1.*ng/2) + 0.5
         zc = (1.*ng/2) + 0.5
         dcross = growth(dble(across)) 
         ak = 2.0 * pi / alambda
         ampl  = 1.d0 / dcross / ak 
         dinit = growth(dble(astart))
         dampl = dinit * ampl
         damplk = dampl * ak
         fp(2) = ak
         fp(3) = dampl         
         x1 = 0. 
         x2 = dble(ng) 
         fp(4) = x1
         fp(5) = x2
         dddt  = dgrowthdt ( dble(astart) ) 
         dddta = dddt * ampl
c
c....    set up gas 
c
         do ic1 = 1 , ncell0 
           hvar(4,ic1) = 0.
           hvar(5,ic1) = 0.
           call iCellijk ( ic1 , i , j , k )
           xcell = float(i) + 0.5
           fp(1) = xcell - 1.
           aq    = zbrent(qsolve, fp, 5 , x1, x2, 1.d-6)
           hvar(1,ic1) = f_gas * dens0 / (1.d0 + damplk * cos(ak*aq))
           hvar(3,ic1) = dddta * sin(ak*aq) * hvar(1,ic1)
           v_tot = dble(hvar(3,ic1))**2
           hvar(8,ic1) = dble(pre0)/dble((gamma-1)) 
           hvar(2,ic1) = hvar(8,ic1) + 
     &                   0.5d0 * v_tot/dble(hvar(1,ic1))
         enddo            
         
         call Magic ( 0 )
         call TreeEoS ( 0 ) 
         call ChooseTimeStep ( dt )
         call DataOutput ( t , dt )
c
c....    set up particles
c
         qfact = float(ng) / float(nrow)
         xcons = dampl
         anow2 = b2a ( dble(t-0.5*dt) )
         vcons = ampl * dgrowthdt ( dble(anow2) )          
         ip = 0 
c         
         pweight = (1.0 - f_gas) * float(ncell0)/float(np)
         if ( nspec .gt. 1 ) then 
           write(*,*) 'error: setting up Pancake IC: nspec > 1 :',nspec
           write(*,*) 'set to nspec = 1 in a_setup.h and rerun'
           stop
         endif
         wpar(1) = pweight
         lsp(1)  = np
         lspecies(1) = np
         wspecies(1) = pweight
         nspecies = 1

         nsp(1,1) = 1
         nsp(1,2) = lsp(1)

         do i = 1 , nrow
           qi = qfact * float(i-1) + 1.5 
           dx  = xcons * sin(ak*(qi-1.)) 
           dvx = vcons * sin(ak*(qi-1.))
           do j = 1 , nrow 
             qj = qfact * float(j-1) + 1.5
             do k = 1 , nrow 
               qk = qfact * float(k-1) + 1.5
               ip  = ip + 1
               x(ip)  = qi + dx
               y(ip)  = qj 
               z(ip)  = qk 
               vx(ip) = dvx 
               vy(ip) = 0.
               vz(ip) = 0.
               pt(ip) = t
               pdt(ip) = dt
               pw(ip) = pweight
             enddo
           enddo
         enddo
         call LL_Construct ()
      ELSE
        call Init_Tree ()

#ifdef HYDRO
          call Read_Gas_Binary ( 'ls1.d', 5 , dt )
          call SetUnits ()        
          call Get_MaxLevelNow ()
          do L = MinLevel , MaxLevelNow         
            call Magic ( L )
            call TreeEoS ( L ) 
          enddo   
#endif
#ifdef NBODY
          call Read_Particles_Binary ()
          call Read_Particle_Time ()  ! reads in dt for each particle
          call LL_Construct ()      

           xrefmin = 129.
           xrefmax = 1.
           yrefmin = 129.
           yrefmax = 1.
           zrefmin = 129.
           zrefmax = 1.
          do ic1 = nsp(1,1) , nsp(1,2)
            if ( x(ic1) .lt. xrefmin ) xrefmin = x(ic1)
            if ( x(ic1) .gt. xrefmax ) xrefmax = x(ic1)
            if ( y(ic1) .lt. yrefmin ) yrefmin = y(ic1)
            if ( y(ic1) .gt. yrefmax ) yrefmax = y(ic1)
            if ( z(ic1) .lt. zrefmin ) zrefmin = z(ic1)
            if ( z(ic1) .gt. zrefmax ) zrefmax = z(ic1)
          enddo
          write(*,*) 'xrefmin,max =',xrefmin,xrefmax
          write(*,*) 'yrefmin,max =',yrefmin,yrefmax
          write(*,*) 'zrefmin,max =',zrefmin,zrefmax


C$OMP PARALLEL DO 
          do ic1 = 1 , np 
            pt(ic1) = t   ! set t to the current time moment
          enddo
c
c.... compute each particle's weight
c
         pweight = (1.0 - f_gas) * float(ncell0)/float(np)
         if ( nspec .gt. 1 ) then 
           write(*,*) 'error: setting up Pancake IC: nspec > 1 :',nspec
           write(*,*) 'set to nspec = 1 in a_setup.h and rerun'
           stop
         endif
         wpar(1) = pweight
         lsp(1)  = np

         nsp(1,1) = 1
         nsp(1,2) = lsp(1)

          do ic0 = 1 , nspec 
            do ic1 = nsp(ic0,1) , nsp(ic0,2)
              pw(ic1) = wpar(ic0)
            enddo
          enddo
#endif

      ENDIF
      
c
      call ViewDump ( 0 , float(ng/2) + 0.5 , 3 )
      call DataOutput( t , dt )

      end
c  
c     -----------------------
      real function a2t ( ad )
c     -----------------------
c     convert expansion factor ad to cosmic time t
c     here only for Omega_0 = 1 universe
c
      real*8 ad 
c
      a2t = ad**1.5
c
      return
      end
c
c     --------------------------------
      real*8 function dgrowthdt ( ae ) 
c     --------------------------------
c
      real*8 ae, dae, twothirds
      parameter ( dae = 1.d-5 ) 
      parameter ( twothirds = 2.d0/3.d0 )     
      real a2t
c
      dgrowthdt = 2.d0 * ae**1.5
c
      return
      end
c
c     --------------------------------
      real*8 function dgrowthdt2 ( ae ) 
c     --------------------------------
c
      real*8 ae, dae
      parameter ( dae = 1.d-6 ) 
      real*8 growth, a2b

      dgrowthdt2 = (growth(ae+dae) - growth(ae-dae))/
     &            (a2b(ae+dae) - a2b(ae-dae))
      return
      end
c
c     ---------------------------------------
      real*8 function qsolve ( xq , fp , np )
c     ---------------------------------------
c
c     solve for the lagrangian coordinates knowing eulerian
c     assuming : 
c                fp(1) = x (Eulerian coordinate)
c                fp(2) = k (wavenumber = 2 pi / lambda)
c                fp(3) = D+(a) * A (growth factor at a times amplitude)
c                fp(4) = xmin of the cube
c                fp(5) = xmax of the cube     
c
      integer np 
      real*8 xq, fp(np)
c
      qsolve = fp(1) - xq - fp(3)*sin(fp(2)*xq) 
c      if ( qsolve .lt. fp(4) ) qsolve = qsolve + (fp(5)-fp(4))
c      if ( qsolve .gt. fp(5) ) qsolve = qsolve - (fp(5)-fp(4))
c      
      return
      end
c
c     -----------------------------
      real*8 function growth ( ae )
c     -----------------------------
c
c     Growth function of linear perturbations
c     using the accurate fitting formula of 
c     Caroll, Press & Turner 1992, ARA&A, 30, 499
c
c     WARNING: now works only for Omega_h = 0.0 (no HDM) universe
c
      real*8 ae
      real*8 Omega, Omegal
      real*8 g_CPT
c
      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'
c
      if ( (Om0 .eq. 1.0) .and. (Oml0 .eq. 0.0) ) then 
        growth = ae
      else
        call Omegas ( dble(Om0) , dble(Oml0) , Omega , Omegal , ae )
        growth = ae * g_CPT ( Omega , Omegal ) / 
     &                  g_CPT ( dble(Om0) , dble(Oml0) )
      endif

      return
      end
c
c     ----------------------------------------
      real*8 function g_CPT ( Omega , Omegal )
c     ----------------------------------------
c
c     approximation function for the growth factor from 
c     Caroll, Press & Turner 1992, ARA&A, 30, 499
c
c
      real*8 Omega , Omegal
      real*8 pow
      parameter ( pow = 4.d0 / 7.d0 )
c
      g_CPT = 2.5d0 * Omega / ( 
     &                          Omega**pow - 
     &                          Omegal     + 
     &                          (1.d0+0.5d0*Omega)*(1.d0+Omegal/70.0)
     &                        )
c
      return
      end
c
c     --------------------------------------------------------------
      subroutine Omegas ( Omega0 , Omegal0 , Omega , Omegal , ae )
c     --------------------------------------------------------------
c
      real*8 Omega0 , Omegal0 , Omega , Omegal , ae    
      real*8 dum 
c
      dum = ae + Omega0 * (1.d0-ae) + Omegal0 * (ae**3-ae)
      Omega  = Omega0 / dum
      Omegal = ae**3 * Omegal0 / dum
c
      return
      end
c
c     ---------------------
      subroutine setphys ()
c     ---------------------
      include 'a_numbers.h'
      return
      end
c
c     ---------------------------
      subroutine tremain ( trem )   
c     ---------------------------
c
      INCLUDE 'a_numbers.h'
      INCLUDE 'a_setup.h'
      INCLUDE 'a_control.h'
c
      trem = tend - t
c
      return
      end
c
c     ----------------------------------
      subroutine DataOutput ( tt , dtt )
c     ----------------------------------
c
      common /initrun_COM/ xc, astart, aend,
     &                     dens0, pre0, across, alambda, gamma
c
      real*8 growth, dgrowthdt, qsolve
      external qsolve      
      real*8 fp(5), fp2(5), x1, x2, b2a
      real*8 aq, aq1, aq2, phia, phia1, phia2, df, trfi2, differ
      real*8 INTEGRATE, rho_func
      external rho_func, INTEGRATE
c
      include 'a_tree.h'
      include 'a_control.h'
c
      j = ng/2
      k = ng/2
      yj = float(j) + 0.5
      zk = float(j) + 0.5
      xc = float(ng/2) + 0.5
      yc = xc - 0.01
      zc = xc - 0.01
c 
      anow =  b2a ( dble(tt) ) 
      anow2 = b2a ( dble(tt-0.5*dtt) )
      dcross = growth(dble(across)) 
      ak = 2.*pi / alambda
      ampl  = 1.d0 / dcross / ak 
      dinit = growth(dble(anow))
      dampl = dinit * ampl
      damplk = dampl * ak
      fp(2) = ak
      fp(3) = dampl   
      x1 = 0. 
      x2 = dble(ng) 
      fp(4) = x1
      fp(5) = x2
      dddt = dgrowthdt ( dble(anow) ) 
      dddta = dddt * ampl
      c1 = 6.d0 / anow
      c12 = 6.d0 / anow2
      c2 = dens0 / ak**2
      c3 = growth(dble(anow))*ampl*ak
      fp(1) = xc-1.
      aqc    = zbrent(qsolve,fp, 5,x1,x2,1.d-6)

         dstart = growth(dble(astart))
         dstampl = dstart * ampl
         dstamplk = dstampl * ak
         fp2(2) = ak
         fp2(3) = dampl         
         x1 = 0. 
         x2 = dble(ng) 
         fp2(4) = x1
         fp2(5) = x2
         dddtst  = dgrowthdt ( dble(astart) ) 
         dddtast = dddtst * ampl

      open ( unit = 25 , file = 'row.dat' ) 
      open ( unit = 26 , file = 'rowan.dat' ) 
      open ( unit = 27 , file = 'pxv.dat' ) 
      write(25,'(1x , ''t = ''      , e14.5, 3x,
     &                 ''dt = ''     , e10.5)' ) t, dt
      write(26,'(1x , ''t = ''      , e14.5, 3x,
     &                 ''dt = ''     , e10.5)' ) t, dt
      dx = .25
      ib = int((x2-x1)/dx)
      xcell = 1. - 0.5 * dx
      gamma = hvar(7,1)
      do ic1 = 1 , ib-1
        xcell = xcell + dx
        fp(1)  = xcell - 1
        fp2(1) = fp(1)
        aq    = zbrent( qsolve, fp, 5, x1, x2, 1.d-6 )
        
        if ( ic1 .gt. 1 ) then 
           fp(1) = xcell - dx - 1
           aq1   = zbrent(qsolve,fp, 5,x1,x2,1.d-6)
           fp(1) = xcell + dx - 1
           aq2   = zbrent(qsolve,fp, 5,x1,x2,1.d-6)
        endif
        rhoa  = dens0 / (1.d0 + damplk * cos(ak*aq))
        rhoai = dens0 / (1.d0 + dstamplk * cos(ak*aq)) 
        ea   = pre0  / (gamma-1.) * (rhoa/rhoai)**gamma
        vxa   = dddta * sin(ak*aq)
        acca  = c1 *(dens0*aq - (xcell-1.))
        aqk = aq*ak
        phia = c1*(c2*(0.5d0*aqk**2 + c3*(aqk*sin(aqk)+cos(aqk)-1.)) - 
     &             5.d-1*(xcell-1)**2)
        if ( ic1 .gt. 1 ) then 
        aqk = aq1*ak
        phia1= c1*(c2*(0.5d0*aqk**2 + c3*(aqk*sin(aqk)+cos(aqk)-1.)) - 
     &             5.d-1*(xcell-1.-dx)**2)
        aqk = aq2*ak
        phia2= c1*(c2*(0.5d0*aqk**2 + c3*(aqk*sin(aqk)+cos(aqk)-1.)) - 
     &             5.d-1*(xcell-1.+dx)**2)
        df = (phia1 + phia2 - 2.d0 * phia)/dx**2
        trfi2 = 6.d0 * anow * dx**2 * ( rhoa - 1. )
        differ = df - trfi2
        endif
        write(26,'(2(2x,e10.5),5(3x, g13.7))' ) 
     &      xcell, aq, rhoa, vxa, acca, phia, ea

      enddo
      close ( 26 )

c
c....    set up gas 
c
      do ic1 = 1 , ncell0 
        call iCellijk ( ic1 , i , jd , kd )
        xcell = float(i) + 0.5
        fp(1) = xcell - 1.
        aq    = zbrent(qsolve, fp, 5 , x1, x2, 1.d-6)
        gacc(ic1) = 0.1 / (1.d0 + damplk * cos(ak*aq))
        gacc(ic1) = pre0/(gamma-1.0)/gacc(ic1)**gamma
        gacc(ic1) = gacc(ic1) * hvar(1,ic1)**gamma
      enddo            

      j = ng/2
      k = ng/2

      do Level = MinLevel , MaxLevel
        IF ( Level .eq. MinLevel ) THEN 
            Posy = float(j) + 0.5
            Posz = float(k) + 0.5
          do i = 1 , ng
            Posx = float(i) + 0.5
            icell = iCellIndex ( i , j , k ) 
            if ( iOctCh(icell) .eq. 0 ) then 
              rd = sqrt((Posy-yc)**2+(Posz-zc)**2)
              fp(1) = Posx - 1
              aq    = zbrent(qsolve,fp, 5,x1,x2,1.d-6)
              rhoa  = dens0 / (1.d0 + damplk * cos(ak*aq))
              ecell1 = hvar(8,icell)
              rhor   = 1.0 / hvar(1,icell)
              vtot  = 0.5d0 *(hvar(3,icell)**2 + hvar(4,icell)**2 + 
     &                        hvar(5,icell)**2) * rhor              
              ecell2 = hvar(2,icell) - vtot
              vxa   = dddta * sin(ak*aq)
              acca = c1*(dens0*aq-(Posx-1.))
              aqk = aq*ak
        phia = c1*(c2*(0.5d0*aqk**2 + c3*(aqk*sin(aqk)+cos(aqk)-1.)) - 
     &             5.d-1*(Posx-1)**2)
              rhocell = hvar(1,icell) 
              write(25,'(1x,i2,2x,i7,3x,e10.5,15(3x, g13.7))' ) 
     &            Level, icell, Posx, hvar(1,icell), 
     &            hvar(3,icell), hvar(4,icell), hvar(5,icell), 
     &            hvar(6,icell), rhoa, vxa, acca, gacc(icell),
     &            var(2,icell), var(3,icell), phia, rhocell, 
     &            ecell1, ecell2
            endif
          enddo
        ELSE
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          cs = 0.51*CellSize(Level)
          cs3 = CellSize(Level)**3
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idcell = icell + ic2
              if ( iOctCh(idcell) .eq. 0 ) then 
                call Ps ( idcell, Posx, Posy, Posz )
                dcy = abs(Posy-yc) 
                dcz = abs(Posz-zc)
                if ( dcy .le. cs .and. dcz .le. cs ) then 
                  fp(1) = Posx - 1
                  aq    = zbrent(qsolve,fp, 5,x1,x2,1.d-6)
                  rhoa  = dens0 / (1.d0 + damplk * cos(ak*aq))
                  vxa   = dddta * sin(ak*aq)
                  acca = c1*(dens0*aq-(Posx-1.))
                  aqk = aq*ak
        phia = c1*(c2*(0.5d0*aqk**2 + c3*(aqk*sin(aqk)+cos(aqk)-1.)) - 
     &             5.d-1*(Posx-1)**2)
                  rhocell = var(1,idcell)/cs3 + 1.
              ecell1 = hvar(8,idcell)
              rhor   = 1.0 / hvar(1,idcell)
              vtot  = 0.5d0 *(hvar(3,idcell)**2 + hvar(4,idcell)**2 + 
     &                        hvar(5,idcell)**2) * rhor              
              ecell2 = hvar(2,idcell) - vtot
                  write(25,'(1x,i2,2x,i7,3x,e10.5,15(3x, g13.7))' ) 
     &             Level, idcell, Posx, hvar(1,idcell), 
     &             hvar(3,idcell), hvar(4,idcell), hvar(5,idcell), 
     &             hvar(6,idcell), rhoa, vxa, acca, gacc(idcell), 
     &             var(2,idcell), var(3,idcell), phia, rhocell,
     &             ecell1, ecell2
                endif
              endif
            enddo
          enddo
        ENDIF
      enddo
c
      close ( unit = 25 )
c
      ip = 1
      qfact = float(ng) / float(nrow)
      xcons = dampl
      anow2 = b2a ( dble(tt-0.5*dtt) )
      vcons = ampl * dgrowthdt ( dble(anow2) )          
c         
      do i = 1 , nrow
        qi =  qfact * float(i-1) + 1.5 
        dx =  qi + xcons * sin(ak*(qi-1.)) 
        dvx = vcons * sin(ak*(qi-1.))
        
        write(27,'(8(1x,g13.7))') x(ip),y(ip),z(ip),
     &                              vx(ip),vy(ip),vz(ip), dx, dvx 
        ip = ip + npage
      enddo
      close ( unit = 27 )
c
      return
      end

c     ---------------------------------------
      subroutine ViewDump ( Lev , zc , iDir )
c     ---------------------------------------
c
      character*256 plotname
      character*5 fstep
      
      parameter ( nvars = 10 ) 
      character*10 varnames(nvars)
      real*8 ve
      real CellPos(3)
      integer iOD(3,2)
      data iOD / 2 , 1 , 1 , 
     &           3 , 3 , 2   / 
c
      real*8 a2b, b2a, qsolve, dgrowthdt
      real   zbrent 
      real*8 fp(5), x1, x2
      external qsolve
      real*8  v_tot
      real*8 growth
      external growth
      common /initrun_COM/ xc, astart, aend,
     &                     dens0, pre0, across, alambda, gamma

      include 'a_tree.h'
      include 'a_control.h'
c
      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1
      call GetName ( imoviestep , fstep )
      if ( iDir .eq. 1 ) then 
      plotname  = 
     & path(1:lpath)//'/'//jobname1(1:ljob)//'_x.v'
      endif
      if ( iDir .eq. 2 ) then 
      plotname  = 
     & path(1:lpath)//'/'//jobname1(1:ljob)//'_y.v'
      endif
      if ( iDir .eq. 3 ) then 
      plotname  = 
     & path(1:lpath)//'/'//jobname1(1:ljob)//'_'//fstep//'_z.v'
      endif

c
      nlayers = 0
      open ( 25 , file = plotname )
      write (unit=25,fmt='(a10,i10,1pe12.4,i10)')
     &            jobname1,imoviestep,aexp(Lev), nlayers
      write (unit=25,fmt='(''zones         1 1 1'')')
c
      nfuncs=8
      varnames(1)='Rho_g'
      varnames(2)='Rho_dm'
      varnames(3)='phi'
      varnames(4)='E'
      varnames(5)='e'
      varnames(6)='Ee'
      varnames(7)='vtot'
      varnames(8)='T'

      
      write ( unit = 25 , fmt = '(i3,20(1x,a10))' )
     &               nfuncs,(varnames(i),i=1,nfuncs)
      ifuncs = 0 
      write ( unit = 25 , fmt = '(i3)' ) ifuncs
c
      call Get_MaxLevelNow ()
c
c.... figure out number of cells first 
c      
      numcell = 0 
      do Level = MinLevel , Lev
        cs = CellSize(Level)/2.          ! halfsize of a cell
        IF ( Level .eq. MinLevel ) THEN 
          do ic = 1 , ncell0
            if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then 
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
              if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                numcell = numcell + 1
              endif
            endif
          enddo ! ic
        ELSE
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then     
                call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
                if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                  numcell = numcell + 1
                endif
              endif
            enddo ! ic2
          enddo ! ic1            
        ENDIF
      end do
c
c.... write cell corners
c
c      write(*,*) 'numcell =',numcell
      write (25,*)  numcell*4

      do Level = MinLevel , Lev
        cs = CellSize(Level)/2.          ! halfsize of a cell
        IF ( Level .eq. MinLevel ) THEN 
          do ic = 1 , ncell0
            if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then 
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
              if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                xi = CellPos(iOD(iDir,1))  ! X-coord of cell center
                yi = CellPos(iOD(iDir,2))  ! Y-coord of cell center
                write (unit=25,fmt='(1p4e12.4)') xi-cs,yi-cs,0.,0.
                write (unit=25,fmt='(1p4e12.4)') xi-cs,yi+cs,0.,0.
                write (unit=25,fmt='(1p4e12.4)') xi+cs,yi+cs,0.,0.
                write (unit=25,fmt='(1p4e12.4)') xi+cs,yi-cs,0.,0.
              endif
            endif
          enddo ! ic
        ELSE
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev )
     +             then     
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
                if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                  xi = CellPos(iOD(iDir,1))  ! X-coord of cell center
                  yi = CellPos(iOD(iDir,2))  ! Y-coord of cell center
                  write (unit=25,fmt='(1p4e12.4)') xi-cs,yi-cs,0.,0.
                  write (unit=25,fmt='(1p4e12.4)') xi-cs,yi+cs,0.,0.
                  write (unit=25,fmt='(1p4e12.4)') xi+cs,yi+cs,0.,0.
                  write (unit=25,fmt='(1p4e12.4)') xi+cs,yi-cs,0.,0.
                endif
              endif
            enddo ! ic2
          enddo ! ic1            
        ENDIF
      end do
c
c.... now write variables 
c      
      nzone = 1
      write (*,*)  numcell
      write (25,*)  numcell
      numcell = 0      

      Tfactor = T_0 / aexpn**2
      do Level = MinLevel , Lev
        cs = CellSize(Level)/2.          ! halfsize of a cell
        cs3 = CellSize(Level)**3
        IF ( Level .eq. MinLevel ) THEN 
          do ic = 1 , ncell0
            if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then 
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
              if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                numcell = numcell + 1
                np1 = numcell * 4 - 3
                np2 = np1 + 1
                np3 = np2 + 1
                np4 = np3 + 1
                v_x = hvar(3,ic)/hvar(1,ic)
                v_y = hvar(4,ic)/hvar(1,ic)
                v_z = hvar(5,ic)/hvar(1,ic)
                v_tot = sqrt(v_x**2+v_y**2+v_z**2)
                rho_dm = (var(1,ic) - hvar(1,ic) + 1.) 
                rho_g  = hvar(1,ic) 
                vsound = sqrt(hvar(7,ic)*hvar(6,ic)/hvar(1,ic))
                ve = (1.d0*hvar(2,ic)-0.5d0*rho_g*v_tot**2)
                vmach  = v_tot / vsound
                T_gas = Tfactor * hvar(6,ic) / hvar(1,ic)
                write (unit=25,fmt='(4i10,i5,1p8e14.6)')
     &                np1,np2,np3,np4,nzone,
     &             rho_g, gacc(ic), var(3,ic), hvar(2,ic), hvar(8,ic),
     &                ve, v_tot, T_gas
              endif
            endif
          enddo ! ic
        ELSE
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then     
                call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
                if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                  numcell = numcell + 1
                  np1 = numcell * 4 - 3
                  np2 = np1 + 1
                  np3 = np2 + 1
                  np4 = np3 + 1
                  v_x = hvar(3,ic)/hvar(1,ic)
                  v_y = hvar(4,ic)/hvar(1,ic)
                  v_z = hvar(5,ic)/hvar(1,ic)
                  v_tot = sqrt(v_x**2+v_y**2+v_z**2)
                  rho_dm = (var(1,ic) - hvar(1,ic) + 1.) 
                  rho_g  = hvar(1,ic) 
                  vsound = sqrt(hvar(7,ic)*hvar(6,ic)/hvar(1,ic))
                  ve = (1.d0*hvar(2,ic)-0.5d0*rho_g*v_tot**2)
                  vmach  = v_tot / vsound
                  T_gas = Tfactor * hvar(6,ic) / hvar(1,ic)
                  write (unit=25,fmt='(4i10,i5,1p8e14.6)')
     &                 np1,np2,np3,np4,nzone,
     &                 rho_g, rho_dm, var(2,ic), hvar(2,ic), hvar(8,ic),
     &                 ve, v_tot, T_gas

                endif
              endif
            enddo ! ic2
          enddo ! ic1            
        ENDIF
      end do
      write (unit=25,fmt='(''interfaces   0'')')
      close ( 25 ) 
c
      write(*,*) 'Saved data for view...'
c
      return
      end
c     --------------------------------
      subroutine GetName (ic , fname )
c     --------------------------------
      integer ic
      character*5 fname
      character n1(10) , n2(10) , n3(10) , n4(10), n5(10)
      data n1 / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      data n2 / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      data n3 / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      data n4 / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      data n5 / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 

      i1 = ic / 10000 + 1
      i2 = (ic - (ic/10000)*10000)/1000 + 1
      i3 = (ic - (ic/1000)*1000) / 100 + 1
      i4 = (ic - (ic/100)*100)/10 + 1
      i5 = (ic - (ic/10)*10) + 1
      fname = n1(i1)//n2(i2)//n3(i3)//n4(i4)//n5(i5)

      return
      end
c
