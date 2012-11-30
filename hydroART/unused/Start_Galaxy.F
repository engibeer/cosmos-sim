c     ========================================================================
c     
c           Start a cosmological simulation 
c
c            Andrey Kravtsov (November 2000)
c
c     ========================================================================


c     --------------------------
      subroutine Init_Run ( dt )                
c     --------------------------       
c
      real*8 a2b
      common /initrun_COM/ xc, fpert, eps
c  
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'
#     include "a_def.h"

      real*8 dt
      
      real xmin,xmax,ymin,ymax,zmin,zmax
      real*8 tempc, coolrate
      character*256 fname, fname_gas
      character*256 fname_part1, fname_part2, fname_part3
c
c
      IF ( start ) THEN 

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append ) 
         write(iOWork,*) 'starting for ICs...'
         close(iOWork)
#endif
          
c
c....    set up particles
c
         xc = float(ng/2) + 1.
         yc = xc
         zc = xc

#ifdef NBODY

           ljob  = index(jobname1, ' ') - 1
           lpath = index(path    , ' ') - 1
c
           fname_part1 = 
     &                    path(1:lpath)//'/PMcrd.DAT '
           fname_part2 = 
     &                    path(1:lpath)//'/PMcrs0.DAT '
           fname_part3 = 
     &                    path(1:lpath)//'/pt.dat '

           call Read_Particles_Binary ( fname_part1 , fname_part2 ) 
           call LL_Construct ()      
 
           xrefmin = 1.
           xrefmax = ng+1.
           yrefmin = 1.0
           yrefmax = ng+1.
           zrefmin = 1.0
           zrefmax = ng+1.
           write(*,*) 'xrefmin,max =',xrefmin,xrefmax
           write(*,*) 'yrefmin,max =',yrefmin,yrefmax
           write(*,*) 'zrefmin,max =',zrefmin,zrefmax

           if ( aexpn .lt. astep ) then 
             call Open_ASCII_File ( iOErr , errorfile , sequent ) 
             write(iOErr,*) '1  *error in Init_Run: aexpn < astep:',
     &          aexpn,astep
             close ( iOErr )
             stop
           endif

           t = a2b ( (1.d0*aexpn) ) -0.5 
           tprev = a2b ( (1.d0*(aexpn-astep)) ) -0.5
           do ic1 = 1 , np
             pt(ic1) = t
             pdt(ic1) = t - tprev
           enddo
c          
c.... compute each particle's weight
c
           Vconst = 2./sqrt(Om0) 
           fgas = Omb0 / Om0       

           wsum = 0.
           do ic0 = 1 , nspecies-1
             
             do ic1 = nsp(ic0,1) , nsp(ic0,2)
              pw(ic1) = (1.0 - fgas) * wpar(ic0)
c               x(ic1) = 0.5 * (x(ic1) - 1.0) + 1.0
c               y(ic1) = 0.5 * (y(ic1) - 1.0) + 1.0
c               z(ic1) = 0.5 * (z(ic1) - 1.0) + 1.0
               vx(ic1) = Vconst * vx(ic1)
               vy(ic1) = Vconst * vy(ic1)
               vz(ic1) = Vconst * vz(ic1)
             enddo
             wsum = wsum + wpar(ic0)*float(nsp(ic0,2)-nsp(ic0,1)+1)
           enddo
           wsum = wsum / ( 1. - Omb0 )
#ifdef WORKOUT
           call Open_ASCII_File ( iOWork , workfile , append ) 
           write(iOWork,*) 'isp =',ic0,'wsum =',wsum,' ncell0 =',ncell0
           close(iOWork)
#endif
#endif

         ainit = aexp0 
         do Level = MinLevel , MaxLevel 
           aexp(Level) = aexpn
         enddo

         boxh   = extras(100) 

         write(*,*) 'initializing for box size =', boxh
c
         call SetUnits ()
#ifdef COOLING 
         call Set_Cooling ()
         call Test_Cooling ( 'testcool.dat ' )
#endif
c
         const = wpar(1)
c
#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'setting refinement thresholds...'
         close(iOWork)
#endif
c
        trho(MinLevel) = 0.01 * const ! 2.0
        tbmass(MinLevel) = 1000.0 
        trho(MinLevel+1) = 0.
        trho(MinLevel+1) = 0.
        tbmass(MinLevel+1) = 1000.0 
        do Level = 3 , MaxLevel 
          trho(Level)  = const * 4.0
          tbmass(Level) = 1000.0
          write(*,*) 'L =',Level,' trho(L) =',trho(Level),
     &               ' trhog(L)=',tbmass(Level)
        enddo
 
        do ic = 1 , mcell
          hvar(1,ic) = 1.e-15
          hvar(2,ic) = 1.e-15
          hvar(3,ic) = 1.e-15
          hvar(4,ic) = 1.e-15 
          hvar(5,ic) = 1.e-15
          hvar(6,ic) = 1.e-15
          hvar(7,ic) = 1.e-15
          hvar(8,ic) = 1.e-15
#ifdef ADVECT_SPECIES
          do j = 9 , nhvar
            hvar(j,ic) = 1.e-15           
          enddo
#endif          
        enddo

        Ndiff = 0
        do iL = MinLevel , MaxLevel 
          call Get_MaxLevelNow ()
          call Assign_Density ( MinLevel , MaxLevelNow , 0 )
          do Level = MinLevel , MaxLevelNow 
            call Modify ( Level , 0 , mtot )
          enddo
          ncell = noct * nchild + ncell0
          call Get_MaxLevelNow ()
          write(*,*) 'MaxLevelNow = ',MaxLevelNow,' Number of Cells=',
     &          iNOLL(MaxLevelNow)*nchild
        enddo

        Ndiff = 4
        call Assign_Density ( MinLevel , MaxLevelNow , 0 ) 
        call Solve_Poisson ( MinLevel , MaxLevelNow , 0 ) 
         
        phimin = 1.e9 

        IF ( MaxLevelNow .gt. MinLevel ) THEN 
          nLevel = iNOLL(MaxLevelNow)
          call Select_Cells ( MaxLevelNow , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idcell = icell + ic2
              if ( var(2,idcell) .lt. phimin ) then                  
                phimin = var(2,idcell)
                imin = idcell
              endif
            enddo
          enddo
        ELSE
          do ic = 1 , ncell0 
            if ( var(2,ic) .lt. phimin ) then                  
              phimin = var(2,ic)
              imin = ic
            endif
          enddo
        ENDIF

        write(*,*) 'imin, phimin =',imin, phimin
        call Ps ( imin , xphimin, yphimin, zphimin ) 

        istep = 0 
c
c....  initialize gas variables
c
       write(*,*) 'initializing gas variables...'
       write(*,*) 'BoxSize =',boxh,' [/h Mpc]; r0 =',r0,' /h Mpc'
       write(*,*) 'aexpn =', aexpn


       a3 = aexpn**3 
       a2 = aexpn**2 
       
       rho0i = a3 / rho0 
       v0i = aexpn / v0 
       e0i = aexpn**5 / P0

       rmin = 0.5 * CellSize ( MaxLevelNow ) * r0 * 1000.0 / hubble 

       do ic = 1 , ncell0 
         call Ps ( ic , xd , yd , zd ) 
         xd = (xd - zphimin)
         yd = (yd - yphimin)
         zd = (zd - zphimin)
         
         xd = xd * r0 * 1000.0 / hubble 
         yd = yd * r0 * 1000.0 / hubble  
         zd = zd * r0 * 1000.0 / hubble 

         call Init_Gas ( xd, yd, zd, rmin, gamma, boxh/hubble, fgas,
     &                   aMbox, energy, dens_gas, px, py, pz )
         hvar(1,ic) = dens_gas * rho0i
         ei = energy * e0i 
         hvar(2,ic) = ei + 0.5 * hvar(1,ic) * v0i * 
     &                sqrt(px**2 + py**2 + pz**2)
         hvar(3,ic) = px * v0i * hvar(1,ic)
         hvar(4,ic) = py * v0i * hvar(1,ic)
         hvar(5,ic) = pz * v0i * hvar(1,ic)
         hvar(6,ic) = (gamma - 1) * ei 
         hvar(7,ic) = gamma
         hvar(8,ic) = ei 
c         write(*,*) 'x,y,z =',zd,yd,zd
c         write(*,*) ' energy, px,py,pz=',energy,px,py,pz
c         pause
       enddo

       pxmin = 1.e20
       pxmax = -1.e20
       pymin = 1.e20
       pymax = -1.e20
       pzmin = 1.e20
       pzmax = -1.e20
       DO Level = MinLevel+1 , MaxLevelNow 
         nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              call Ps ( ic , xd , yd , zd ) 
              xd = (xd - zphimin) * r0 * 1000.0 / hubble
              yd = (yd - yphimin) * r0 * 1000.0 / hubble 
              zd = (zd - zphimin) * r0 * 1000.0 / hubble 
              call Init_Gas ( xd, yd, zd, rmin, gamma, 
     &                       boxh/hubble, fgas, aMbox, 
     &                       energy, dens_gas, px, py, pz )
              hvar(1,ic) = dens_gas * rho0i
              ei = energy * e0i 
              hvar(2,ic) = ei + 0.5 * hvar(1,ic) * v0i * 
     &                     sqrt(px**2 + py**2 + pz**2)
              hvar(3,ic) = px * v0i * hvar(1,ic)
              hvar(4,ic) = py * v0i * hvar(1,ic)
              hvar(5,ic) = pz * v0i * hvar(1,ic)
              hvar(6,ic) = (gamma - 1) * ei 
              hvar(7,ic) = gamma
              hvar(8,ic) = ei 
              pxmin = min(pxmin,px)
              pxmax = max(pxmax,px)
              pymin = min(pymin,py)
              pymax = max(pymax,py)
              pzmin = min(pzmin,pz)
              pzmax = max(pzmax,pz)
            enddo
          enddo
       ENDDO

       write(*,*) 'pxmin,max =',pxmin,pxmax
       write(*,*) 'pymin,max =',pymin,pymax
       write(*,*) 'pzmin,max =',pzmin,pzmax

       write(*,*) 'Mass in the box =',aMbox,' [Msun]'
c
c....  re-define mass scale 
c
       C_fb = C_fb / aM0
       RIaf = RIaf / aM0
       aM0  = aMbox / ncell0 
       C_fb = C_fb * aM0
       RIaf = RIaf * aM0

       DO Level = MinLevel , MaxLevelNow-1 
         call SplitUpdate ( Level ) 
       ENDDO
c
         call Magic ( 0 )
#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'setting pressure... gamma =',gamma
         close(iOWork)
#endif

         call TreeEoS ( 0 ) 

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'choosing time step...'
         close(iOWork)
#endif

         call ChooseTimeStep ( dt )
         dt = 1.e-5 !!!!! temporary, to start w/o restarts in runs with feedback
c
        write(*,*) 'zphimin =', zphimin
        
        call ViewDump ( MaxLevelNow , zphimin , 3 )
c        write(*,*) 'xphimin =', xphimin
c        call ViewDump ( MaxLevelNow , xphimin , 1 )
c
        call Save ( dt , 0 ) 
        start = .false. 
c
      ELSE

#ifdef HYDRO
          ljob  = index(jobname1, ' ') - 1
          lpath = index(path    , ' ') - 1
          fname  =
     &          path(1:lpath)//'/'//jobname1(1:ljob)//'.d '

          call Read_Gas_Binary ( dt , fname )
          call SetUnits ()
#ifdef COOLING
          call Set_Cooling ()        
#endif
          call Get_MaxLevelNow ()
          do L = MinLevel , MaxLevelNow         
            call Magic ( L )
            call TreeEoS ( L ) 
          enddo   
#endif
c
c....  re-define mass scale 
c
          C_fb = C_fb / aM0
          RIaf = RIaf / aM0
          aM0 = 4.286924E+09 / ncell0 
          C_fb = C_fb * aM0
          RIaf = RIaf * aM0
          
c
	  call Init_Tree ()
c
#ifdef NBODY
          ljob  = index(jobname1, ' ') - 1
          lpath = index(path    , ' ') - 1
c
          fname_part1 =
     &                   path(1:lpath)//'/PMcrd.DAT '
          fname_part2 =
     &                   path(1:lpath)//'/PMcrs0.DAT '
          fname_part3 =
     &                   path(1:lpath)//'/pt.dat '

          call Read_Particles_Binary ( fname_part1 , fname_part2 )
          call Read_Particle_Time ( fname_part3 )  ! reads in dt for each particle

#ifdef STARFORM
          fname_part3 =
     &                   path(1:lpath)//'/stars.dat '
          call Read_Stellar_Data  ( fname_part3 )
#endif
          call LL_Construct ()      

#endif

      ENDIF
c
c      const =  (1.0 - Omb0/Om0)
c      constb = Omb0/Om0

        const = wpar(1)

        trho(MinLevel) = 2.0 * const
        tbmass(MinLevel) = 1000.0 
        trho(MinLevel+1) = 0.
        trho(MinLevel+2) = 0.
        tbmass(MinLevel+1) = 1000.0 
        do Level = 3 , MaxLevel 
          trho(Level)  = 4.0 * const 
          tbmass(Level) = 1000.0
          write(*,*) 'L =',Level,' trho(L) =',trho(Level),
     &               ' trhog(L)=',tbmass(Level)
        enddo

      do i = 1 , n_save 
        if ( asave(i) .le. aexpn ) then 
          isave(i) = 1
        else
          isave(i) = 0 
        endif      
      enddo

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
      common /initrun_COM/ xc, fpert, eps
c
      return
      end
c
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
      SAVE idumpp
      data iOD / 2 , 1 , 1 , 
     &           3 , 3 , 2   / 
      DATA idumpp/0/
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1
      idumpp = idumpp +1
c      call GetName ( istep , fstep )
      call GetName ( idumpp , fstep )
      if ( iDir .eq. 1 ) then 
      plotname  = 
     & path(1:lpath)//'/'//jobname1(1:ljob)//'_'//fstep//'_x.v'
      endif
      if ( iDir .eq. 2 ) then 
      plotname  = 
     & path(1:lpath)//'/'//jobname1(1:ljob)//'_'//fstep//'_y.v'
      endif
      if ( iDir .eq. 3 ) then 
      plotname  = 
     & path(1:lpath)//'/'//jobname1(1:ljob)//'_'//fstep//'_z.v'
      endif

c
      nlayers = 0
      open ( 25 , file = plotname )
      write (25,25)
     &            jobname1,imoviestep,aexp(Lev), nlayers
      write (25,26)
 25   format(a10,i10,1pe12.4,i10)
 26   format('zones         1 1 1')
c
      nfuncs=10
      varnames(1)='Rho_g'
      varnames(2)='T_g'
      varnames(3)='v'
      varnames(4)='phi'
      varnames(5)='Rho'
      varnames(6)='Lev'
      varnames(7)='vx'
      varnames(8)='vy'
      varnames(9)='ZII'
      varnames(10)='ZIa'
      
      write ( 25 , 27 )
     &               nfuncs,(varnames(i),i=1,nfuncs)
 27   format(i3,6(1x,a10))
      ifuncs = 0 
      write ( 25 , 28) ifuncs
 28   format(i3)
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
                write (25,29) xi-cs,yi-cs,0.,0.
                write (25,29) xi-cs,yi+cs,0.,0.
                write (25,29) xi+cs,yi+cs,0.,0.
                write (25,29) xi+cs,yi-cs,0.,0.
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
                  write (25,29) xi-cs,yi-cs,0.,0.
                  write (25,29) xi-cs,yi+cs,0.,0.
                  write (25,29) xi+cs,yi+cs,0.,0.
                  write (25,29) xi+cs,yi-cs,0.,0.
                endif
              endif
            enddo ! ic2
          enddo ! ic1            
        ENDIF
      end do
 29   format(1p4e13.6)
c
c.... now write variables 
c      
      nzone = 1
      write (*,*)  numcell
      write (25,*)  numcell
      numcell = 0      
      a3 = aexpn**3
      a2 = aexpn**2
      Tcons =  T_0 / a2

      Scons = S_0
      vfact = v0 / aexpn 

      do Level = MinLevel , Lev
        fjfact = 2.**Level*sqrt(2.0*pi**2/3.*2.5/aexpn)
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
                T_g = Tcons * hvar(6,ic) / hvar(1,ic)
                v_tot = sqrt(hvar(3,ic)**2 +
     &                       hvar(4,ic)**2 + 
     &                       hvar(5,ic)**2 ) /
     &                                         hvar(1,ic) * vfact

      ZII = 0.0 
      ZIa = 0.0
#ifdef ENRICH
                ZII = hvar(izII,ic) / hvar(1,ic) / 0.01989
#endif
#ifdef ENRICH_SNIa 
                ZIa = hvar(izIa,ic) / hvar(1,ic) / 0.01989
#endif
              write(25,30)
     &                np1,np2,np3,np4,nzone,
     &                hvar(1,ic), T_g, v_tot, var(2,ic),
     &                var(1,ic)/cs3,Level,
     &                hvar(3,ic)/ hvar(1,ic) * vfact,
     &                hvar(4,ic)/ hvar(1,ic) * vfact,
     &                ZII, ZIa
              endif
            endif
          enddo ! ic
 30       format(4(i7,1x),i5,1x,1p5e14.6,i3,4e14.5)
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
                  T_g = Tcons * hvar(6,ic) / hvar(1,ic)
                  v_tot = sqrt(hvar(3,ic)**2 +
     &                         hvar(4,ic)**2 +
     &                         hvar(5,ic)**2) /
     &                                          hvar(1,ic) * vfact
      ZII = 0.0 
      ZIa = 0.0
#ifdef ENRICH
                ZII = hvar(izII,ic) / hvar(1,ic) / 0.01989
#endif
#ifdef ENRICH_SNIa 
                ZIa = hvar(izIa,ic) / hvar(1,ic) / 0.01989
#endif
              write(25,30)
     &                np1,np2,np3,np4,nzone, 
     &                hvar(1,ic), T_g, v_tot, var(2,ic),
     &                var(1,ic)/cs3,Level,
     &                hvar(3,ic)/ hvar(1,ic) * vfact,
     &                hvar(4,ic)/ hvar(1,ic) * vfact,
     &                ZII, ZIa
                endif
              endif
            enddo ! ic2
          enddo ! ic1            
        ENDIF
      end do
      write (25,31)
 31   format('interfaces   0')
      close ( 25 ) 
c
      write(*,*) 'Saved data for view...'
c
      return
      end
c

c
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
c
C----------------------------------------------  
      SUBROUTINE Init_Gas(x1,y1,z1, rmin, gamma, boxsize, fgas, aMbox, 
     *     energy,dens_gas,px,py,pz)           
C----------------------------------------------  
c
c     input:  x1, y1, z1 - coordinates of the cell with respect to the center in kpc (no h!)
c             rmin  - minimum radius in kpc 
c             gamma - gas polytropic index 
c             boxsize in Mpc (no h!)
c             fgas  - gas fraction in the halo
c     output: 
c             aMbox - total mass in the computational box
c                     computed only on the first call 
c             energy - internal energy per unit volume [ergs cm^-3]
c             dens_gas - gas density in Msun Mpc^-3
c             px, py, pz - velocity components in km/s 
c
c      parameters.h should match with /param/
      IMPLICIT REAL*4 (A-H,O-Z)  
      Common /param/xout,xtrunc,Overdens,CON, 
     *     OM0,aMass_vir,alpha,
     *     h,Rvir,Rsnfw,V_0,rho_1,rho_0,e_0,
     *     x0,p0,energytrunc,energyout,fout 
      DATA Index/1/
      SAVE Index,dens_gastrunc,dens2trunc
      external Mmx,F,AINTG,Mmx_ext  

      If(Index.eq.1)then
         Index =2
c                        Initialize parameters
      OM0 = .3 ! Omega matter at z=0
      Overdens = 340.  
      Con = 15.  ! Halo Concentration     
      aMass_vir = 1.4e8  !Msun   
      
      Xout = Con                ! DM truncation
      Xtrunc = 2.*Xout          !Gas truncation radius   
       
      fout = F(xout)            ! DM Mass ( units aMass_vir/F(c))           
      
  
      h = .7                    !hubble/100 
      Rvir     = 442.97/h * (h*aMass_vir/1.e+11/(Om0*Overdens))**0.3333 ! kpc    
      Rsnfw = Rvir/Con          !kpc
      Fcon = F(con)             ! Halo Virial Mass ( units aMass_vir/F(c))  
      aM0      = (aMass_vir)/Fcon ! in units of Msun, massscale       
      t_0       = 4.693e+11/sqrt(aM0/Rsnfw**3) ! in   yrs   
c     V_0= 9.767e+8 *Rsnfw/(t_0)   ! in units of km/sec       
      V_0 = sqrt(aM0/Rsnfw)*2.047e-3 !km/s    
      rho_1 = 5.81e-33*(aM0 /Rsnfw**3 ) !gm/cm3  
      rho_0 = 7.95e7*( aM0/Rsnfw**3) ! Msun/Mpc3 
      e_0 = 2.55e-28*(aM0**2/Rsnfw**4 ) ! g cm-1 s-2   
   
      alpha = fgas 
c     alpha: dens(gas)/dens(dm)  
c     gamma: gas index 
    
      x0 = .02  ! Point of Match between analytical and numerical solution 
      xinf = 500. 
      p0 = log(x0)*.5 + 
     *     AINTG(Mmx,x0,xout) + fout*AINTG(Mmx_ext,xout,xinf)   
      energytrunc = fout*AINTG(Mmx_ext,xtrunc,xinf) 
      energyout = fout*AINTG(Mmx_ext,xout,xinf)
      dens_gastrunc=alpha*rho_0/max(xtrunc*(1.+xtrunc)**2,1.e-10)
      dens2trunc=alpha*rho_1/max(xtrunc*(1. + xtrunc)**2,1.e-10) 
         write (*,*) ' Vir mass      =', aMass_vir 
         write (*,*) ' Vir radius    =', Rvir 
         write (*,*) ' Concentrat    =', Con
         write (*,*) ' Gamma         =', gamma
         write (*,*) ' Trunc gas(kpc)=', Xtrunc*Rsnfw
         write (*,*) ' dens_gastrunc =', dens_gastrunc
         Rout = Xtrunc * Rsnfw * 1.e-3
         rho_const = fgas * rho_0 / Xtrunc / (1.0 + Xtrunc)**2
         aMbox = rho_const * (boxsize - 33.510322 * (Rout)**3) +
     &           (1.-fgas) * aMass_vir + 
     &           fgas * aMass_vir * F(Xtrunc)/Fcon 
         write (*,*) ' Mbox [Msun]   =', aMbox         
      Endif  
     
      px = 0. 
      py = 0.
      pz = 0.   

      xd = x1 / Rsnfw
      yd = y1 / Rsnfw
      zd = z1 / Rsnfw
      x = max ( Sqrt(xd**2 + yd**2 +zd**2 ) , rmin / Rsnfw )
      
      If(x.lt.xtrunc)then  
         
         If(x.le.x0)then  
            energy = p0 -  log(x)*.5   
            dens_gas = alpha*rho_0/max(x*(1. +x)**2,1.e-10) 
            dens2 = alpha*rho_1/max(x*(1. + x)**2,1.e-10)
         Else             
            If(x.lt.xout)then  
               energy = AINTG(Mmx,x,xout) +  energyout     
            else
               energy =   energyout -fout*AINTG(Mmx_ext,xout,x)
            Endif
            dens_gas = alpha*rho_0/max(x*(1. + x)**2,1.e-10)
            dens2 = alpha*rho_1/max(x*(1. + x)**2,1.e-10) 
             
            
         Endif
         
      Else           
         energy = energytrunc 
         dens_gas = dens_gastrunc
         dens2 =   dens2trunc 
      Endif
       
      
      energy = energy*alpha*e_0 /(gamma-1.) !energy density (g cm-1 s-2)      
c      dens_gas = alpha*rho_0/max(x*(1. + x)**2,1.e-10) !density (Msun/Mpd3)       
c      dens2 = alpha*rho_1/max(x*(1. + x)**2,1.e-10) ! density (g/ cm3)  
      
      r = x*Rsnfw  !r kpc
      px = V_0*px  ! velocity km/s  
      py = V_0*py
      pz = V_0*pz
            


      Return
      End

C---------------------------------------------- 
      SUBROUTINE Init_Gas_Disk (x1,y1,z1,rmin,gamma,   
     *    energy,dens_gas,px,py,pz)           
C---------------------------------------------- 

c
c     input:  x1, y1, z1 - coordinates of the cell with respect to the center in kpc (no h!)
c             rmin  - minimum radius in kpc 
c             gamma - gas polytropic index 
c     output: 
c             energy - internal energy per unit volume [ergs cm^-3]
c             dens_gas - gas density in Msun Mpc^-3
c             px, py, pz - velocity components in km/s 
c
c      parameters.h should match with /param/
 
      IMPLICIT REAL*4 (A-H,O-Z)  
 
      Common /param/xout,xtrunc,Overdens,CON, 
     *     OM0,aMass_vir,alpha,
     *     h,Rvir,Rsnfw,V_0,rho_1,rho_0,e_0,
     *     x0,p0,energytrunc,energyout,fout,
     *     aMassd,rd,z0,aMg,alpha_d,RMvirMd,
     *      fcon,rdtrunc,zdtrunc,e0d,Tgasd0,
     *     p0d
      DATA Index/1/
      SAVE Index,dens_gastrunc,dens2trunc,
     *     rhodg0,rhodg2,rdtrunc2,vcdtrunc,
     *     fcdtrunc,Tgasext0
      external Mmx,F,AINTG,Mmx_ext,Mmd 
      external  bessi0,bessk0,bessi1,bessk1,
     *     Zdensity,Vcd_ext  

      If(Index.eq.1)then
         Index =2 
c                        Initialize parameters
      OM0 = .3 ! Omega matter at z=0
      Overdens = 340.  
      Con = 15.  ! Halo Concentration     
      aMass_vir = 1.4e8  !Msun   
       
      Xout = Con                ! DM truncation
      Xtrunc = 2.*Xout          !Gas truncation radius   
       
      fout = F(xout)            ! DM Mass ( units aMass_vir/F(c))           
      
  
      h = .7                    !hubble/100 
      Rvir     = 442.97/h * (h*aMass_vir/1.e+11/(Om0*Overdens))**0.3333 ! kpc    
      Rsnfw = Rvir/Con          !kpc
      Fcon = F(con)             ! Halo Virial Mass ( units aMass_vir/F(c))  
      aM0      = (aMass_vir)/Fcon ! in units of Msun, massscale       
      t_0       = 4.693e+11/sqrt(aM0/Rsnfw**3) ! in   yrs   
c     V_0= 9.767e+8 *Rsnfw/(t_0)   ! in units of km/sec       
      V_0 = sqrt(aM0/Rsnfw)*2.047e-3 !km/s    
      rho_1 = 5.81e-33*(aM0 /Rsnfw**3 ) !gm/cm3  
      rho_0 = 7.95e7*( aM0/Rsnfw**3) ! Msun/Mpc3 
      e_0 = 2.55e-28*(aM0**2/Rsnfw**4 ) ! g cm-1 s-2   
      amu = .6 ! mean molecular weight 
      Temp0 = 1.212e-8*amu! 7.269e-9 ! kelvin 
   
      alpha =.01 
c      gamma = 1.666   
c     alpha: dens(gas)/dens(dm)  
c     gamma: gas index 
 
        
c****************************************************************
C    Disk Parameters   
      aMassd = 4.2e10           ! Msun
      rd = 3.2                  !   kpc 
      z0 = .325                 !  kpc
      alpha_d = 10.     
      aMg = aMassd/alpha_d  
      RMvirMd =aMass_vir/aMassd  
      V0d = 2.1087e-3*Sqrt(aMassd/rd) !km/s 
      rhodg0 = 7.95e7*aMassd/(z0*rd**2)/alpha_d !Msun/Mpc3  
      rhodg2 = 4.95e-33*aMassd/(z0*rd**2)/alpha_d !g/cm3  
      Tgasd0 = 5.08214e-4*(z0*aMassd/rd**2*
     *     (1.+1./alpha_d)*amu)   
      e0d = Tgasd0*rhodg2*8.254e7/amu
 
      rdtrunc = 5. !rd units 
      rdtrunc2 = 3.*rdtrunc
      zdtrunc = 3.*z0/rd !rd units 
      vcdtrunc = Vcd_ext(rdtrunc)
      fcdtrunc = F(rdtrunc*rd/Rsnfw)   
      Tgasext0 = (1.+1./alpha_d)*alpha_d 

c****************************************************************


    
      x0 = .02  ! Point of Match between analytical and numerical solution 
      xinf = 500. 
 

      energytrunc = fout*AINTG(Mmx_ext,xtrunc,xinf) 
      energyout = fout*AINTG(Mmx_ext,xout,xinf) 
      dens_gastrunc=alpha*rho_0/max(xtrunc*(1.+xtrunc)**2,1.e-10)
      dens2trunc=alpha*rho_1/max(xtrunc*(1. + xtrunc)**2,1.e-10) 

       
      
      p0 = log(x0)*.5 +   
     *     AINTG(Mmx,x0,xout) + fout*AINTG(Mmx_ext,xout,xinf)   
      
        
      p0d = fout/RMvirMd*(x0*(Rsnfw/rd)  - log(x0) )* 
     *     (Rsnfw/rd)**2 +  
     *     AINTG(Mmd,x0,xout)*fout/RMvirMd +  
     *     energyout/RMvirMd       
      
      

      
         write (*,*) ' Vir mass      =', aMass_vir 
         write (*,*) ' Vir radius    =', Rvir 
         write (*,*) ' Concentrat    =', Con
         write (*,*) ' Gamma         =', gamma
         write (*,*) ' Trunc gas(kpc)=', Xtrunc*Rsnfw
         write (*,*) ' dens_gastrunc =', dens_gastrunc
         write (*,*) ' Temp_0 =',Temp0 
         write (*,*) ' V_0 =', V_0
         write (*,*) ' Md = ',aMassd
         write (*,*) ' rd(kpc) = ',rd
         write (*,*) ' z0(kpc)= ',z0
         write (*,*) ' V0d = ',V0d 
         write (*,*) 'rhodg0(Msun/Mpc3) = ',rhodg0 
         write (*,*) 'rhodg0(Msun/pc3)= ',rhodg0*1.e-18
         write (*,*) 'rhodg2(g/cm3)= ',rhodg2
         write (*,*) 'Tgasd0 (K)= ',Tgasd0
         write (*,*) 'e0d= ',e0d


      Endif  
     
 

      px = 0. 
      py = 0.
      pz = 0.   


      xd2 = x1 / Rsnfw 
      yd2 = y1 / Rsnfw
      zd2 = z1 / Rsnfw
 
      x = max (Sqrt(xd2**2 + yd2**2 +zd2**2),rmin/Rsnfw) ! Spherical radius

       
                                  
      

c**************************************************************      
      If(x.lt.xtrunc)then   
           
         If(x.le.x0)then   
            energy = p0 -  log(x)*.5   + 
     *           p0d -   
     *          fout/RMvirMd*(log(x) - x*(Rsnfw/rd))*(Rsnfw/rd)**2


            dens_gas = alpha*rho_0/max(x*(1. +x)**2,1.e-10) 
            dens2 = alpha*rho_1/max(x*(1. + x)**2,1.e-10)
         Else              
            If(x.lt.xout)then  
               energy = AINTG(Mmx,x,xout)+
     *              energyout*(1.+1./RMvirMd)+
     *              AINTG(Mmd,x,xout)*fout/RMvirMd
            else 
               energy =   energyout*(1.+1./RMvirMd) -
     *              (1.+1./RMvirMd)*fout*AINTG(Mmx_ext,xout,x)   
            Endif   
            dens_gas = alpha*rho_0/max(x*(1. + x)**2,1.e-10) 
            dens2 = alpha*rho_1/max(x*(1. + x)**2,1.e-10) 
              
         Endif
          
      Else           
         energy = energytrunc*(1.+1./RMvirMd)    
         dens_gas = dens_gastrunc 
         dens2 =   dens2trunc 
      Endif

      energy = energy*alpha*(e_0)/(gamma-1.)

c************************************************************** 
C     Disk Gas initialization

       
 
      xd = xd2*Rsnfw/rd  
      yd = yd2*Rsnfw/rd 
      zd = zd2*Rsnfw/rd  
      rc = Sqrt(xd**2 + yd**2) ! Cylindrical radius (Rd units) 
       
        

      If(rc.le.rdtrunc.and.abs(zd).le.zdtrunc)then 

         
 
      Vphi2 = (Vcd(rc) - 3.*z0/rd*
     *     exp(-rc)*rc)*(1.+1./alpha_d)   
      

      Vphi = Sqrt(Vphi2)        !Rotation velocity in disk units 

      vxd = yd*(-1.)/max(rc,1.e-8)*Vphi !Cartesian Gas Velocities 
      vyd = xd/max(rc,1.e-8)*Vphi  
      vzd = 0. 
      
      
      densdiskgas = rhodg0*Exp(-rc)*Zdensity(zd) 
      densdiskgas2 = rhodg2*Exp(-rc)*Zdensity(zd)
      Tgasd = Tgasd0*Exp(-rc)     
      egasd = e0d*Exp(-2.*rc)/(gamma -1.)       
      

      energy = energy + egasd
      dens_gas = dens_gas + densdiskgas
      dens2 = dens2 + densdiskgas2 
         
         px = px + vxd*V0d      ! velocity km/s  
         py = py + vyd*V0d
         pz = pz + vzd*V0d

       

      else If(rc.le.rdtrunc2.and.abs(zd).le.zdtrunc)then
         
c     vcdtrunc stellar: component total contribution(until rtrunc)
c     Vcd_ext(rc)/alpha: gas component contribution at rc
c     F(rc*rd/Rsnfw)/(rc*rd/Rsnfw): halo contribution at rc
   
         y  = rc*rd/Rsnfw   
         Vphi2 = vcdtrunc +  
     *        RMvirMd/Fcon*F(y)/rc+Vcd_ext(rc)/alpha_d -  
     *      3.*z0/rd*exp(-rc)*rc/Tgasext0 
  

         

         Vphi = Sqrt(Vphi2)     !Rotation velocity in disk units 
         vxd = yd*(-1.)/max(rc,1.e-8)*Vphi !Cartesian Gas Velocities 
         vyd = xd/max(rc,1.e-8)*Vphi  
         vzd = 0. 
         
         densdiskgas = rhodg0*Exp(-rc)*Zdensity(zd) 
         densdiskgas2 = rhodg2*Exp(-rc)*Zdensity(zd)
         Tgasd = Tgasd0*Exp(-rc)/Tgasext0    
         egasd = e0d*Exp(-2.*rc)/(gamma -1.)/Tgasext0        
         

         energy = energy + egasd
         dens_gas = dens_gas + densdiskgas
         dens2 = dens2 + densdiskgas2 
         
         px = px + vxd*V0d      ! velocity km/s  
         py = py + vyd*V0d
         pz = pz + vzd*V0d
         
         if(dens_gas.lt.0.)write(*,*)dens_gas

      endif 
       
      
      

c**************************************************************

        
       
c      dens_gas = dens_gas + densdisk_gas !Total density (Msun/Mpd3)
c      dens2 =dens2 + densdisk_gas2  ! Total density (g/ cm3)
c      energy = energy*alpha*(e_0)/(gamma-1.)+
c     *     egas_d               !Total energy density (g cm-1 s-2)      

    

      r = x*Rsnfw  !r kpc

            
   
     

      Return
      End
C----------------------------------------------
      SUBROUTINE TRAPZD(FUNCC,A,B,S,N)
C--------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
        SAVE IT
        EXTERNAL FUNCC
      IF (N.EQ.1) THEN
        S=0.5d0*(B-A)*(FUNCC(A)+FUNCC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5D0*DEL
        SUM=0.0D0
        DO 11 J=1,IT
          SUM=SUM+FUNCC(X)
          X=X+DEL
11      CONTINUE
        S=0.5D0*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
   
C----------------------------------- Simpson integration
      REAL*4 FUNCTION AINTG(FUNC,A,B)
C--------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      PARAMETER (EPS=2.0d-5, JMAX=22) 
      EXTERNAL FUNC
c      REAL*4 AINTG 
      OST=-1.d30 
      OS  = -1.d30
      ST   =0. 
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,ST,J)
        AINTG=(4.0d0*ST-OST)/3.0d0
        IF (ABS(AINTG-OS).Le.EPS*ABS(OS).and.J.GT.6) RETURN
        OS=AINTG
        OST=ST
11    CONTINUE
      WRITE (16,*)'Integration did not converge'
      RETURN 
      END
C---------------------------------------
      REAL*4 Function Mmx(x) 
C---------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
 
         Mmx =(LOG(1.+x) -x/(1.+x))/x**3/(1.+x)**2
      Return
      End 
C---------------------------------------
      REAL*4 Function Mmx_ext(x) 
C---------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
c      REAL*4  x
      Mmx_ext =(1./x**3/(1.+x)**2 ) 
      Return
      End 
C--------------------------------------- 
      REAL*4 Function F(x)   ! x =r/r_s,  F= mass(r)/4pi r_s^3
C---------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
         y =max(x,1.e-10) 
         F =Max(LOG(1.+y) -y/(1.+y),1.e-16)
      Return
      End  
C--------------------------------------- 
C--------------------------------------------------------------
      REAL*4 function Vcd(x)!Total Vcircular**2  
C--------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      Common /param/xout,xtrunc,Overdens,CON, 
     *     OM0,aMass_vir,alpha,
     *     h,Rvir,Rsnfw,V_0,rho_1,rho_0,e_0,
     *     x0,p0,energytrunc,energyout,fout,
     *     aMassd,rd,z0,aMg,alpha_d,RMvirMd,
     *     fcon,rdtrunc,zdtrunc,e0d,Tgasd0,
     *     p0d
C    Mvir/Md= 1.26d12/4.d10 =31.5 

      
      External  bessi0,bessk0,bessi1,bessk1,F 
 
     

      y=.5*x 
      Vcd = 2.*y*y*( bessi0(y) *bessk0(y)-bessi1(y) *bessk1(y))  +  
     *     RMvirMd/Fcon*(F(x*Rd/Rsnfw)/x)   
     
      
      Return     
      End

C--------------------------------------------------------------
      REAL*4 function Vcd_ext(x)!Extended Disk Vcircular**2  
C--------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      Common /param/xout,xtrunc,Overdens,CON, 
     *     OM0,aMass_vir,alpha,
     *     h,Rvir,Rsnfw,V_0,rho_1,rho_0,e_0,
     *     x0,p0,energytrunc,energyout,fout,
     *     aMassd,rd,z0,aMg,alpha_d,RMvirMd,
     *     fcon,rdtrunc,zdtrunc,e0d,Tgasd0,
     *     p0d

      
      External  bessi0,bessk0,bessi1,bessk1,F 
 
     

      y=.5*x 
      Vcd_ext = 2.*y*y*( bessi0(y) *bessk0(y)-bessi1(y) *bessk1(y))   

      
      
      Return     
      End

C---------------------------------------------------------------
      REAL*4 Function bessi0(x) ! Modified Bessel I_o Function  
C----------------------------------------------------------- 
      IMPLICIT REAL*4(A-H,O-Z)
   
c      Real*4 x                  !bessi0, 
C     Returns the modified Bessel Function  I_o(x) for any real x
c      Real ax
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y ! Accumulate polinomials in double precision
      Save  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      Data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *     1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      Data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *     0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     *     0.2635537d-1,-0.1647633d-1,0.392377d-2/ 
      
      if (abs(x).lt.3.75)then
         y=(x/3.75)**2
         bessi0 = p1 + y*(p2+y*(p3 + y*(p4 + y*(p5 + y*(p6 +y*p7)))))
      else
         ax = abs(x)
         y = 3.75/ax
         bessi0 = (exp(ax)/sqrt(ax))*(q1 + y*(q2 +y*(q3 + y*(q4 + 
     *        y*(q5 +y*(q6 + y*(q7 + y*(q8 + y*q9)))))))) 
      endif

      Return
      End 


C------------------------------------------------------------
      REAL*4 Function bessk0(x) ! Modified Bessel K_o Function 
C----------------------------------------------------------- 
 
C     Returns the modified Bessel Function  K_o(x) for positive real x  
      IMPLICIT REAL*4(A-H,O-Z)
c      Real*4 x                  
      external bessi0  
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      Save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7 
      Data p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *     0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/ 
      Data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *     -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/ 
      
      If(x.le.2.0)then  
         y = x*x/4.0
         bessk0 = (-log(x/2.0)* bessi0(x)) + (p1 + y*(p2 +y*(p3 + 
     *        y*(p4 + y*(p5 + y*(p6 + y*p7))))))
      Else 
         y = (2.0/x)
         bessk0 = ( exp(-x)/sqrt(x))*(q1 + y*(q2 + y*(q3 +
     *        y*(q4 + y*(q5 + y*(q6 + y*q7))))))
      Endif
      Return 
      End
C------------------------------------------------------------
      REAL*4 Function bessi1(x) ! Modified Bessel I_1 Function 
C----------------------------------------------------------- 
      IMPLICIT REAL*4(A-H,O-Z)
c      Real*4 x   
C     Returns the modified Bessel Function  I_o(x) for any real x 
c      Real ax
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y ! Accumulate polinomials in double precision
      Save  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      Data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *     0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      Data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *     -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     *     -0.2895312d-1,0.1787654d-1,-0.420059d-2/
      If(abs(x).lt.3.75)then
         y = (x/3.75)**2  
         bessi1 = x*(p1 + y*(p2 +y*(p3+ y*(p4+ y*(p5+ y*(p6 + y*p7)))))) 
      Else
         ax = abs(x) 
         y = 3.75/ax
         bessi1 = (exp(ax)/sqrt(ax))*(q1 + y*(q2 + y*(q3 + y*(q4 +
     *        y*(q5 + y*(q6 + y*(q7 + y*(q8 + y*q9))))))))
         if(x.lt.0.)bessi1= -bessi1
      Endif
      Return 
      End

C------------------------------------------------------------
      REAL*4 Function bessk1(x) ! Modified Bessel K_1 Function 
C----------------------------------------------------------- 
      IMPLICIT REAL*4(A-H,O-Z)
c      Real*4 x   
c      Real*4 bessi1 
      external bessi1 
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      Save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7 
      Data p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *     -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      Data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *     0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      If(x.le.2.0)then
         y = x*x/4.0
         bessk1 = (log(x/2.0)*bessi1(x)) + (1.0/x)*(p1 + y*(p2+
     *        y*(p3+y*(p4+y*(p5 + y*(p6 + y*p7)))))) 
      Else
         y = 2.0/x  
         bessk1 = ( exp(-x)/sqrt(x))*(q1 + y*(q2+y*(q3 +
     *        y*(q4+y*(q5 + y*(q6 + y*q7))))))
      Endif
      Return
      End
 
C------------------------------------------------------------
      REAL*4 Function Zdensity(x)
C------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      
      Zdensity = 4./(Exp(2.*x)+ 2. +Exp(-2.*x))
       
      return
      end
C
C--------------------------------------------------------------
      REAL*4 Function Mmd(x) !Disk contribution for the halo velocities.
C-----------------------------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
      
      Common /param/xout,xtrunc,Overdens,CON, 
     *     OM0,aMass_vir,alpha,
     *     h,Rvir,Rsnfw,V_0,rho_1,rho_0,e_0,
     *     x0,p0,energytrunc,energyout,fout,
     *     aMassd,rd,z0,aMg,alpha_d,RMvirMd,
     *     fcon,rdtrunc,zdtrunc,e0d,Tgasd0,
     *     p0d
          
      External d

      
      Mmd = d(x*Rsnfw/(Rd))/x**3/(1.d0+x)**2         
        
      
      
      Return  
      End 
C--------------------------------------------------------------------------------- 
      REAL*4 Function d(x) ! This is the dimensionless disk mass. 
C----------------------------------------------------------- 
      IMPLICIT REAL*4(A-H,O-Z)
      Common /param/xout,xtrunc,Overdens,CON, 
     *     OM0,aMass_vir,alpha,
     *     h,Rvir,Rsnfw,V_0,rho_1,rho_0,e_0,
     *     x0,p0,energytrunc,energyout,fout,
     *     aMassd,rd,z0,aMg,alpha_d,RMvirMd, 
     *     fcon,rdtrunc,zdtrunc,e0d,Tgasd0,
     *     p0d
      
 
      Rt = rdtrunc   ! disk truncation  
      If(x.le.RT)then  
         d = (1.d0 - EXP(-X)*(1.d0 + X))/        
     *        (1.d0 - EXP(-RT)*(1.d0 + RT)) 
      Else
         d =1.   
      Endif 
c     The normalization is made with respect the finite size of our disk
c     (1.d0 - 1.d0/EXP(RT)*(1.d0 + RT))* zcoef   not to the infinite disk 
c     Zcoeff is cancelled out  because of the normalization factor (March 2 01)
c      d = (1.d0 - 1.d0/EXP(X)*(1.d0 + X))*zcoef   this the case of an infinite disk
c     zcoef : mass integrated in the verticla direction at radius: r 
      Return   
      End   
C--------------------------------------------------------------
 
