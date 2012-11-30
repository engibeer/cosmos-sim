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
      real*8 dt
CEVERINO09232005      real*8 a2bcc
      REAL*8 a2b
      common /initrun_COM/ xc, fpert, eps
c  
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'
#include "a_def.h"
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
c           fname_part1 = 
c     &                    path(1:lpath)//'/IC/PMcrd_a.DAT '
c           fname_part2 = 
c     &                    path(1:lpath)//'/IC/PMcrs0_a.DAT '

           fname_part1 =
     &                    path(1:lpath)//'/PMcrd.DAT '
           fname_part2 =
     &                    path(1:lpath)//'/PMcrs0.DAT '

           call Read_Particles_Binary ( fname_part1 , fname_part2 ) 
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
        if ( x(ic1).lt.0. .or. y(ic1).lt.0. .or. z(ic1).lt. 0.) then 
               write(*,*) ic1, x(ic1),y(ic1),z(ic1)
             endif
           enddo
           write(*,*) 'xrefmin,max =',xrefmin,xrefmax
           write(*,*) 'yrefmin,max =',yrefmin,yrefmax
           write(*,*) 'zrefmin,max =',zrefmin,zrefmax

           if ( aexpn .lt. astep ) then 
             write(*,*) '*error in Init_Run: aexpn < astep:',aexpn,astep
             stop
           endif
           t = a2b ( (1.d0*aexpn) ) 
           tprev = a2b ( (1.d0*(aexpn-astep)) )
           do ic1 = 1 , np
             pt(ic1) = t
             pdt(ic1) = t - tprev
           enddo
CEVERINO09182005
           WRITE (*,*)   'a= ', aexpn ,' t =', t
c          
c.... compute each particle's weight
c
           wsum = 0.
           do ic0 = 1 , nspec 
             do ic1 = nsp(ic0,1) , nsp(ic0,2)
               pw(ic1) = wpar(ic0)
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

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'aexpn =',aexpn,' aexp0 =',aexp0,
     &        ' astep =',astep
         close(iOWork)
#endif

         ainit = aexp0 
         do Level = MinLevel , MaxLevel 
           aexp(Level) = aexpn
         enddo

         fname_gas  = 
     &     path(1:lpath)//'/tr_ic.dat'

         open ( 50 , file = fname_gas, form = 'unformatted' )
         read(50) boxh
c         boxh = 80.0
         read(50) ainit, astep
         read(50) ncells

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'reading IC grid: ainit =',ainit,
     &                   ' ncells =',ncells
         if ( ncells .ne. ncell0 ) then 
           write(iOWork,*) 'bad ICs: ncells .ne. ncell0'
         endif
         close(iOWork)
#endif
         call SetUnits ()

#ifdef COOLING
         call Set_Cooling ()
#endif
         read(50) (hvar(1,ic),ic=1,ncell0)
         read(50) (hvar(3,ic),ic=1,ncell0)
         read(50) (hvar(4,ic),ic=1,ncell0)
         read(50) (hvar(5,ic),ic=1,ncell0)
         read(50) (hvar(2,ic),ic=1,ncell0)
         read(50) (hvar(8,ic),ic=1,ncell0)
         close(50)
c
         rhomin = 1.e20
         rhomax = -1.e20
         tgmin = 1.e20
         tgmax = -1.e20
         do ic = 1 , ncell0 
           hvar(7,ic) = gamma
           tgas = T_0 * 
     &         (gamma-1.) * hvar(8,ic) / hvar(1,ic) / ainit / ainit
           rhomin = min(rhomin,hvar(1,ic))
           rhomax = max(rhomax,hvar(1,ic))
           tgmin = min(tgmin,tgas)
           tgmax = max(tgmax,tgas)
         enddo
         write(*,*) 'rhogmin,max=',rhomin,rhomax
         write(*,*) 'Tgmin,max=',tgmin,tgmax
c         
         istep = 0 
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

c
         call ChooseTimeStep ( dt )  
      
#ifdef AGNFEEDBACK
       nAGN = 0  ! Initialize the AGN feedback
#endif

      ELSE

#ifdef HYDRO
        ljob  = index(jobname1, ' ') - 1
        lpath = index(path    , ' ') - 1
        fname =
     &        path(1:lpath)//'/'//jobname1(1:ljob)//'.d '

        call Read_Gas_Binary ( dt , fname )

c       For analysis of old outputs
c        if ( irun .eq. 0 ) then
c          call Read_Gas_Binary_old ( dt , fname )
c        else if ( irun .eq. 1 ) then
c          call Read_Gas_Binary_Adiabatic ( dt , fname )
c        else if ( irun .eq. 2 ) then
c          call Read_Gas_Binary ( dt , fname )
c        endif

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
	call Init_Tree ()
c
#ifdef NBODY
        ljob  = index(jobname1, ' ') - 1
        lpath = index(path    , ' ') - 1
c
        fname_part1 =
     &                 path(1:lpath)//'/PMcrd.DAT '
        fname_part2 =
     &                 path(1:lpath)//'/PMcrs0.DAT '
        fname_part3 =
     &                 path(1:lpath)//'/pt.dat '

        call Read_Particles_Binary ( fname_part1 , fname_part2 )
        call Read_Particle_Time ( fname_part3 )  ! reads in dt for each particle
#ifdef STARFORM
        fname_part3 =
     &               path(1:lpath)//'/stars.dat '

        if ( irun .eq. 2 ) then
          call Read_Stellar_Data  ( fname_part3 )
        endif
c        call Write_Stars_ASCII ( 'stdump.dat ' ) 
#endif
        call LL_Construct ()      
#endif


#ifdef AGNFEEDBACK
        call Read_AGNs () 
#endif

      ENDIF

c
c...  Setting refinement thresholds
c
      const =  (1.0 - Omb0/Om0)
      constb = Omb0/Om0
c
#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'setting refinement thresholds...'
         close(iOWork)
#endif
c
c$$$      trho(MinLevel) = const
c$$$      trho(MinLevel+1) = 0.0001 * const
c$$$      trho(MinLevel+2) = 0.0001 * const
c$$$      tbmass(MinLevel) = 1.0*constb
c$$$      tbmass(MinLevel+1) = 1.0*constb
c$$$      tbmass(MinLevel+2) = 1.0*constb
c$$$      do Level = 3 , MaxLevel 
c$$$        trho(Level)  = const / 32  
c$$$        tbmass(Level) = constb / 16.0 
c$$$        write(iOut,*) 'L =',Level,' trho(L) =',trho(Level),
c$$$     &             ' trhog(L)=',tbmass(Level)
c$$$      enddo

c
         trho(MinLevel) = wpar(1)*2.
         trho(MinLevel+1) =  wpar(1)*2.
         trho(MinLevel+2) = wpar(1)*2.
      tbmass(MinLevel) = 1.*constb
      tbmass(MinLevel+1) = 1.*constb
      tbmass(MinLevel+2) = 1.*constb
      do Level=3,4
         trho(Level) = wpar(1)*2.
         tbmass(Level) = constb/3. 
      enddo
      do Level = 5 , MaxLevel 
        trho(Level)  =wpar(1)*2.
        tbmass(Level) = constb/10.  
      enddo
      do Level=MinLevel, MaxLevel
         write(*,*) 'L =',Level,' trho(L) =',trho(Level),
     &        ' trhog(L)=',tbmass(Level)
      enddo

CEVERINO42205
c$$$      trho(MinLevel) = 6. * const
c$$$      trho(MinLevel+1) = 4. * const
c$$$      trho(MinLevel+2) = 4. * const
c$$$      tbmass(MinLevel) = 4. *constb
c$$$      tbmass(MinLevel+1) = 4. *constb
c$$$      tbmass(MinLevel+2) = 4. *constb
c$$$
c$$$      do Level = 3 , MaxLevel 
c$$$        trho(Level)  =4 * const 
c$$$        tbmass(Level) = 4*constb  
c$$$        write(iOut,*) 'L =',Level,' trho(L) =',trho(Level),
c$$$     &             ' trhog(L)=',tbmass(Level)
c$$$      enddo

c
c...  Refine meshes until the refinement converges. (Daisuke 4/3/03)
c
      write(*,*)
      write(*,*) ' Refining meshes immediately down to 2nd level'
      if ( start ) then
         call Get_MaxLevelNow ()
         call Assign_Density ( MinLevel, MaxLevelNow, 0 )
         call Modify ( MinLevel, 1, mtot )
         write(*,*) 'Level = ',MinLevel,' mtot = ',mtot
         Level = 1
         do while ( mtot .ne. nil ) 
            call Assign_Density ( Level, Level, 0 )
            call Modify ( Level, 1, mtot )
            write(*,*) 'Level = ',Level,' mtot = ',mtot
            Level = Level + 1           
         enddo
      endif
        
c
c...  Preparting arrays for saving moments
c
      write(*,*) 'preparing arrays'
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
c     ---------------------------------------
      subroutine ViewDump ( Lev , zc , iDir )
c     ---------------------------------------
c
      character*256 plotname
      character*5 fstep
      parameter ( nvars = 12 )
      character*10 varnames(nvars)
      real*8 ve
      real CellPos(3)
      integer iOD(3,2)
      SAVE idumpp1, idumpp2, idumpp3
      data iOD / 2 , 1 , 1 , 
     &           3 , 3 , 2   / 
      DATA idumpp1/0/
      DATA idumpp2/0/
      DATA idumpp3/0/
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1
c      call GetName ( istep , fstep )
      if ( iDir .eq. 1 ) then
        idumpp1 = idumpp1 + 1
        call GetName ( idumpp1, fstep )
        plotname  = 
     &    path(1:lpath)//'/'//jobname1(1:ljob)//'_'//fstep//'_x.v'
      endif
      if ( iDir .eq. 2 ) then
         idumpp2 = idumpp2 + 1
         call GetName ( idumpp2, fstep )
         plotname  = 
     &     path(1:lpath)//'/'//jobname1(1:ljob)//'_'//fstep//'_y.v'
      endif
      if ( iDir .eq. 3 ) then
         idumpp3 = idumpp3 + 1
         call GetName ( idumpp3, fstep )
         plotname  = 
     &     path(1:lpath)//'/'//jobname1(1:ljob)//'_'//fstep//'_z.v'
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
      nfuncs=12
      varnames(1)='Rho_g'
      varnames(2)='T_g'
      varnames(3)='P_g'
      varnames(4)='S_g'
      varnames(5)='v'
      varnames(6)='phi'
      varnames(7)='Rho'
      varnames(8)='Lev'
      varnames(9)='vx'
      varnames(10)='vy'
      varnames(11)='ZII'
      varnames(12)='ZIa'
      
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
      write (25,*)  numcell
      numcell = 0     
      a5 = aexpn**5
      a3 = aexpn**3
      a2 = aexpn**2
      Tcons =  0.666666 * T_0 / a2

      Scons = S_0
      vfact = v0 / aexpn 
      sfact = 21.5953 * (0.666666) * r0**2 * (Om0/hubble)**(1/3) 

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
                T_g = Tcons * hvar(8,ic) / hvar(1,ic)
                P_g = P0 * hvar(6,ic) / a5
                S_g = sfact * hvar(8,ic) / hvar(1,ic)**(5./3.)
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
     &                hvar(1,ic), T_g, P_g, S_g, v_tot, 
     &                var(2,ic), var(1,ic)/cs3, Level,
     &                hvar(3,ic)/ hvar(1,ic) * vfact,
     &                hvar(4,ic)/ hvar(1,ic) * vfact,
     &                ZII, ZIa
              endif
            endif
          enddo ! ic
 30       format(4(i7,1x),i5,1x,1p7e14.6,i3,4e14.5)
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
                  T_g = Tcons * hvar(8,ic) / hvar(1,ic)
                  P_g = P0 * hvar(6,ic) / a5
                  S_g = sfact * hvar(8,ic) / hvar(1,ic)**(5./3.)
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
     &                hvar(1,ic), T_g, P_g, S_g, v_tot, 
     &                var(2,ic), var(1,ic)/cs3, Level,
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
      write(*,*) 'Saved',numcell,' cells for view...'
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
c$$$c
c$$$c     ----------------------------------
c$$$      subroutine DataOutput ( tt , dtt )
c$$$c     ---------------------------------- 
c$$$c
c$$$c     ViewDump slicing through the max gas density or
c$$$c     the min potential at every small step
c$$$c
c$$$#     include "a_def.h"
c$$$      include 'a_tree.h'
c$$$      include 'a_control.h'
c$$$c
c$$$      real*8 tt, dtt 
c$$$      character*256 plotname, pname, sname
c$$$      character*5 fstep  
c$$$      logical lpointsout,l_force_center
c$$$      integer Nbb(6)
c$$$      real*8 n_H, nH_Nb, T_Nb, p_Nb(6), v_x, v_y, v_z
c$$$c
c$$$
c$$$      if ( mod(istep,1) .eq. 0 ) then
c$$$c
c$$$           rhomax = -1.e6
c$$$           phimin = 1.e6
c$$$	   do Level = 1 , MaxLevelNow
c$$$             nLevel = iNOLL(Level)
c$$$             call Select_Cells ( Level , nLevel ) 
c$$$             do i1 = 1 , nLevel
c$$$               icell = iSelect(i1)
c$$$               do i2 = 0 , 7
c$$$                 idcell = icell + i2
c$$$                 if ( iOctCh(idcell) .eq. nil ) then   ! only leaves
c$$$                   if ( hvar(1,idcell) .gt. rhomax ) then 
c$$$                     call Ps (idcell,Posx,Posy,Posz)
c$$$                     if ( Posx .lt. 40. .and. Posy .gt. 10.) then 
c$$$                     rhomax = hvar(1,idcell)
c$$$                     irhomax = idcell
c$$$                     endif
c$$$                   endif
c$$$                   if ( var(2,idcell) .lt. phimin ) then 
c$$$                     call Ps (idcell,Posx,Posy,Posz)
c$$$                     if ( Posx .lt. 40. .and. Posy .gt. 10.) then 
c$$$                     phimin = var(2,idcell)
c$$$                     iphimin = idcell
c$$$                     endif
c$$$                   endif
c$$$                 endif
c$$$               enddo
c$$$             enddo
c$$$            enddo
c$$$            irhomax = max(irhomax,1)
c$$$            iphimin = max(iphimin,1)
c$$$            call Ps ( irhomax , xc0 , yc0 , zc0 )
c$$$            call Ps ( iphimin , xc , yc , zc )
c$$$            write(iOut,*) 'Output Data:'
c$$$            write(iOut,*) 'irhomax =',irhomax,' rhog=',hvar(1,irhomax)
c$$$            write(iOut,*) 'xc, yc, zc =',xc0,yc0,zc0
c$$$            write(iOut,*) 'iphimin =',iphimin,' rhog=',hvar(1,iphimin)
c$$$            write(iOut,*) 'xc, yc, zc =', xc, yc, zc
c$$$            
c$$$            t_gas = (gamma - 1.0) * hvar(8,irhomax)/hvar(1,irhomax)
c$$$     &         * T_0/aexpn**2
c$$$             write(iOut,*) 'maxrho:',irhomax,iLv(irhomax),
c$$$     &         hvar(1,irhomax),t_gas,hvar(8,irhomax)
c$$$c            irhomax = ivas 
c$$$
c$$$             rmin = 2.5e-5/r0
c$$$             rmax = 1.0/r0
c$$$             nrbin = 80 
c$$$             ibintype = 0 
c$$$             deltavir = 340
c$$$             ifindcenter = 3 
c$$$             lpointsout = .false.  
c$$$             l_force_center = .false.  
c$$$             nstep = 0 
c$$$             lpath = index(path    , ' ') - 1
c$$$             call GetName ( istep , fstep )
c$$$
c$$$             plotname  =
c$$$     &         path(1:lpath)//'/'//'pro'//'_'//fstep//'.dat '
c$$$
c$$$              call Output_Average_Profiles ( xc , yc , zc , rvir, 
c$$$     &                             rmin , rmax , nrbin , ibintype ,  
c$$$     &                             deltavir, ifindcenter , 
c$$$     &                             lpointsout, l_force_center,
c$$$     &                             plotname ) 
c$$$
c$$$              call ViewDump(MaxLevelNow,zc,3)
c$$$              call ViewDump(MaxLevelNow,xc,1)
c$$$              call ViewDump(MaxLevelNow,yc,2)
c$$$
c$$$       endif
c$$$
c$$$       return
c$$$       end
c$$$c
c$$$c     ----------------------------------
c$$$      subroutine DataOutput1 ( tt , dtt )
c$$$c     ----------------------------------
c$$$c
c$$$c     Extracting 3D boxs interpolated onto the uniform grids,
c$$$c     profiles and ViewDump
c$$$c
c$$$      common /initrun_COM/ xc, fpert, eps
c$$$c
c$$$      include 'a_tree.h'
c$$$      include 'a_control.h'
c$$$      include 'a_analysis.h'
c$$$
c$$$      character*5 fstep
c$$$      real svd(nsvx,nsvy,nsvz)
c$$$      real p_d1(nsvy,nsvz), p_d2(nsvy,nsvz), p_d3(nsvy,nsvz)
c$$$      real p_t1(nsvy,nsvz), p_t2(nsvy,nsvz), p_t3(nsvy,nsvz)
c$$$      real*8 sum1, sum2, sum3, sum4, vmax1, vmax2, em
c$$$      common / SVDUM / svd, p_d1, p_d2, p_d3, p_t1, p_t2, p_t3
c$$$      character*256 fname_gas
c$$$      logical lpointsout, l_force_center
c$$$      integer imovstep
c$$$c
c$$$c.... output profiles 
c$$$c
c$$$      xc = 16.53711
c$$$      yc = 22.15430
c$$$      zc = 128.7441 
c$$$      rmin = 1.e-3/r0
c$$$      rmax = 10.0/r0
c$$$      nrbin = 45 
c$$$      ibintype = 0 
c$$$      deltavir = 340
c$$$      ifindcenter = 3 
c$$$      lpointsout = .false.  
c$$$      l_force_center = .false.  
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$
c$$$      open ( 22 , file = 'nstep' )
c$$$      read(22,*) imovstep       
c$$$      close ( 22 )
c$$$
c$$$      call GetName ( imovstep , fstep )
c$$$
c$$$      fname_gas  =
c$$$     &         path(1:lpath)//'/'//'pro'//'_'//fstep//'.dat '
c$$$
c$$$      call Output_Average_Profiles ( xc , yc , zc , rvir,  
c$$$     &                            rmin , rmax , nrbin , ibintype ,  
c$$$     &                            deltavir, ifindcenter , 
c$$$     &                            lpointsout, l_force_center,
c$$$     &                            fname_gas ) 
c$$$c
c$$$c.... extract boxes 
c$$$c
c$$$      Posx = xc
c$$$      Posy = yc
c$$$      Posz = zc
c$$$      sizebox2 = 4./2.
c$$$      dx = 2.* sizebox2 / r0 / float(nsvx)
c$$$      dy = 2.* sizebox2 / r0 / float(nsvy)
c$$$      dz = 2.* sizebox2 / r0 / float(nsvz)
c$$$      xgmin = Posx - sizebox2/r0
c$$$      xgmax = Posx + sizebox2/r0
c$$$      ygmin = Posy - sizebox2/r0
c$$$      ygmax = Posy + sizebox2/r0
c$$$      zgmin = Posz - sizebox2/r0
c$$$      zgmax = Posz + sizebox2/r0
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax , 
c$$$     &                      ygmin , ygmax , 
c$$$     &                      zgmin , zgmax , 1 , 1 ) ! gas density 
c$$$
c$$$      rhofact = rho0 / aexpn**3 ! Msun / Mpc^3
c$$$      write(*,*) 'rhofact =',rhofact
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            svd(ix,iy,iz) = sv(ix,iy,iz)
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) * rhofact 
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      write(*,*) 'extracting pressure...'      
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 6 )  ! pressure 
c$$$
c$$$      tfact = T_0 / aexpn**2 
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * tfact       
c$$$            svd(ix,iy,iz) = svd(ix,iy,iz) * rhofact
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = 1 , nsvx
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nsvx
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L4Mpc_yz_full.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      ib = nsvx / 2
c$$$      dslice = 1.0 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L4Mpc_yz_d1Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      dslice = 0.5 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L4Mpc_yz_d0.5Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      dslice = 0.1 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L4Mpc_yz_d0.1Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$c
c$$$c.... now 2/h Mpx box
c$$$c
c$$$c
c$$$      sizebox2 = 2./2.
c$$$      dx = 2.* sizebox2 / r0 / float(nsvx)
c$$$      dy = 2.* sizebox2 / r0 / float(nsvy)
c$$$      dz = 2.* sizebox2 / r0 / float(nsvz)
c$$$      xgmin = Posx - sizebox2/r0
c$$$      xgmax = Posx + sizebox2/r0
c$$$      ygmin = Posy - sizebox2/r0
c$$$      ygmax = Posy + sizebox2/r0
c$$$      zgmin = Posz - sizebox2/r0
c$$$      zgmax = Posz + sizebox2/r0
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax , 
c$$$     &                      ygmin , ygmax , 
c$$$     &                      zgmin , zgmax , 1 , 1 ) ! gas density 
c$$$
c$$$      rhofact = rho0 / aexpn**3 ! Msun / Mpc^3
c$$$      write(*,*) 'rhofact =',rhofact
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            svd(ix,iy,iz) = sv(ix,iy,iz)
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) * rhofact 
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      write(*,*) 'extracting pressure...'      
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 6 )  ! pressure 
c$$$
c$$$      tfact = T_0 / aexpn**2 
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * tfact       
c$$$            svd(ix,iy,iz) = svd(ix,iy,iz) * rhofact
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = 1 , nsvx
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nsvx
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L2Mpc_yz_full.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      ib = nsvx / 2
c$$$      dslice = 1.0 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L2Mpc_yz_d1Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      dslice = 0.5 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L2Mpc_yz_d0.5Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      dslice = 0.1 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L2Mpc_yz_d0.1Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$c
c$$$c.... now 8/h Mpx box
c$$$c
c$$$c
c$$$      sizebox2 = 8./2.
c$$$      dx = 2.* sizebox2 / r0 / float(nsvx)
c$$$      dy = 2.* sizebox2 / r0 / float(nsvy)
c$$$      dz = 2.* sizebox2 / r0 / float(nsvz)
c$$$      xgmin = Posx - sizebox2/r0
c$$$      xgmax = Posx + sizebox2/r0
c$$$      ygmin = Posy - sizebox2/r0
c$$$      ygmax = Posy + sizebox2/r0
c$$$      zgmin = Posz - sizebox2/r0
c$$$      zgmax = Posz + sizebox2/r0
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax , 
c$$$     &                      ygmin , ygmax , 
c$$$     &                      zgmin , zgmax , 1 , 1 ) ! gas density 
c$$$
c$$$      rhofact = rho0 / aexpn**3 ! Msun / Mpc^3
c$$$      write(*,*) 'rhofact =',rhofact
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            svd(ix,iy,iz) = sv(ix,iy,iz)
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) * rhofact 
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      write(*,*) 'extracting pressure...'      
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 6 )  ! pressure 
c$$$
c$$$      tfact = T_0 / aexpn**2 
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * tfact       
c$$$            svd(ix,iy,iz) = svd(ix,iy,iz) * rhofact
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = 1 , nsvx
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nsvx
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L8Mpc_yz_full.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      ib = nsvx / 2
c$$$      dslice = 1.0 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L8Mpc_yz_d1Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      dslice = 0.5 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L8Mpc_yz_d0.5Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      dslice = 0.1 / r0  ! slice thickness -> code units
c$$$      nrange = int( dslice / dx ) + 1
c$$$      nr1 = ib - nrange/2
c$$$      nr2 = ib + nrange/2
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          sum1 = 0.0
c$$$          sum2 = 0.0
c$$$          sum3 = 0.0
c$$$          sum4 = 0.0
c$$$          vmax1 = -1.e6
c$$$          vmax2 = -1.e6
c$$$          do ix = nr1 , nr2
c$$$            sum1 = sum1 + svd(ix,iy,iz)
c$$$            vmax1 = max(vmax1,1.d0*svd(ix,iy,iz))
c$$$            vmax2 = max(vmax2,1.d0*sv(ix,iy,iz))
c$$$            em = (1.d0*svd(ix,iy,iz))**2 * sqrt(sv(ix,iy,iz))
c$$$            sum2 = sum2 + sv(ix,iy,iz) * em
c$$$            sum4 = sum4 + svd(ix,iy,iz) * em
c$$$            sum3 = sum3 + em
c$$$          enddo
c$$$          sum1 = sum1 / nrange
c$$$          sum4 = sum4 / sum3
c$$$          p_d1(iy,iz) = sum1
c$$$          p_d2(iy,iz) = sum4
c$$$          p_d3(iy,iz) = vmax1
c$$$          sum2 = sum2 / sum3
c$$$          p_t1(iy,iz) = sum2
c$$$          p_t2(iy,iz) = vmax2
c$$$        enddo
c$$$      enddo
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/proj'//'_'//fstep
c$$$     &            //'_L8Mpc_yz_d0.1Mpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      
c$$$      write(20,*) p_d1
c$$$      write(20,*) p_d2
c$$$      write(20,*) p_d3
c$$$      write(20,*) p_t1
c$$$      write(20,*) p_t2
c$$$
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      if ( lviewout .and. mod(istep,1) .eq. 0) then 
c$$$        zd = zc
c$$$        call ViewDump ( MaxLevelNow , zd , 3 )
c$$$c       call ViewDump ( MaxLevelNow , Posx , 1 )
c$$$c       call ViewDump ( MaxLevelNow , Posy , 2 )
c$$$      endif
c$$$c
c$$$      open ( 22 , file = 'nstep' )
c$$$      imovstep = imovstep + 1
c$$$      write(22,*) imovstep
c$$$      close ( 22 ) 
c$$$c
c$$$      return
c$$$      end
c$$$c
c$$$c     -----------------------------------
c$$$      subroutine DataOutput2 ( tt , dtt )
c$$$c     -----------------------------------
c$$$c
c$$$c     Extracting 3D boxs interpolated onto the uniform grids
c$$$c
c$$$      common /initrun_COM/ xc, fpert, eps
c$$$c
c$$$      include 'a_tree.h'
c$$$      include 'a_control.h'
c$$$      include 'a_analysis.h'
c$$$
c$$$      character*5 fstep
c$$$      real svd(nsvx,nsvy,nsvz), svd1(nsvx,nsvy,nsvz)
c$$$      real svd2(nsvx,nsvy,nsvz), svd3(nsvx,nsvy,nsvz)
c$$$      common / SVDUM / svd, svd1, svd2, svd3
c$$$      character*256 fname_gas
c$$$      real*8 rhogl, Tgas, Zmet, rs, xion
c$$$      logical lpointsout, l_force_center
c$$$      integer imovstep
c$$$      character digits(0:9)
c$$$      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
c$$$c
c$$$c.... output profiles 
c$$$c
c$$$      xc = 32.96973
c$$$      yc = 13.57129      
c$$$      zc = 39.05566
c$$$      rmin = 1.e-3/r0
c$$$      rmax = 10.0/r0
c$$$      nrbin = 45 
c$$$      ibintype = 0 
c$$$      deltavir = 340
c$$$      ifindcenter = 3 
c$$$      lpointsout = .false.  
c$$$      l_force_center = .false.  
c$$$
c$$$      ic = int(aexpn*1000)
c$$$      i1 = ic / 1000 
c$$$      i2 = (ic - (ic/1000)*1000)/100 
c$$$      i3 = (ic - (ic/100)*100) / 10 
c$$$      i4 = (ic - (ic/10)*10)
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$      fname_gas = path(1:lpath)//'/pro'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'.dat '
c$$$
c$$$      call Output_Average_Profiles ( xc , yc , zc , rvir,  
c$$$     &                            rmin , rmax , nrbin , ibintype ,  
c$$$     &                            deltavir, ifindcenter , 
c$$$     &                            lpointsout, l_force_center,
c$$$     &                            fname_gas ) 
c$$$c
c$$$c.... extract boxes 
c$$$c
c$$$      Posx = xc
c$$$      Posy = yc
c$$$      Posz = zc
c$$$      sizebox2 = 0.1875/2.
c$$$      dx = 2.* sizebox2 / r0 / float(nsvx)
c$$$      dy = 2.* sizebox2 / r0 / float(nsvy)
c$$$      dz = 2.* sizebox2 / r0 / float(nsvz)
c$$$      xgmin = Posx - sizebox2/r0
c$$$      xgmax = Posx + sizebox2/r0
c$$$      ygmin = Posy - sizebox2/r0
c$$$      ygmax = Posy + sizebox2/r0
c$$$      zgmin = Posz - sizebox2/r0
c$$$      zgmax = Posz + sizebox2/r0
c$$$
c$$$      tfact = 0.666666 * T_0 / aexpn**2 
c$$$      rcons = 1.125e-5 * hubble * hubble * Om0 / aexpn**3
c$$$      vcons = v0 / aexpn
c$$$      rhofact = rho0 / aexpn**3 ! Msun / Mpc^3
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax , 
c$$$     &                      ygmin , ygmax , 
c$$$     &                      zgmin , zgmax , 1 , 1 ) ! gas density 
c$$$
c$$$      write(*,*) 'rhofact =',rhofact
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            svd(ix,iy,iz) = sv(ix,iy,iz)
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) * rcons
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      fname_gas = path(1:lpath)//'/n_H'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      write(*,*) 'extracting pressure...'      
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 8 )  ! pressure 
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * tfact       
c$$$            svd1(ix,iy,iz) = sv(ix,iy,iz)
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$
c$$$      fname_gas = path(1:lpath)//'/tg'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      write(*,*) 'extracting velocity components...'      
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 3 )  ! v_x
c$$$
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vcons
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      fname_gas = path(1:lpath)//'/v_x'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 4 )  ! v_y
c$$$
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vcons
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      fname_gas = path(1:lpath)//'/v_y'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 5 )  ! v_z
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vcons
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      fname_gas = path(1:lpath)//'/v_z'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      write(*,*) 'extracting ZIa...'
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , izIa )  ! Z_Ia
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) / 0.01989
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      fname_gas = path(1:lpath)//'//zIa'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      write(*,*) 'extracting ZII...'
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , izII )  ! Z_II
c$$$
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) / 0.01989            
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      fname_gas = path(1:lpath)//'//zII'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      write(*,*) 'computing ionized fractions...'
c$$$      rs = 1.d0 / aexpn - 1.0 
c$$$
c$$$      open ( 33 , file = 'test_fion' ) 
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            Zmet = log10(sv(ix,iy,iz))
c$$$            rhogl = log10(svd(ix,iy,iz) * rcons)
c$$$            Tgas = svd1(ix,iy,iz)
c$$$            call X_ion ( rhogl, Tgas, Zmet, rs, xion )
c$$$            sv(ix,iy,iz) = xion 
c$$$            if ( rhogl .ge. 2 ) then 
c$$$              write(33,93) rhogl, Tgas, Zmet, xion
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$ 93   format(4(e12.6,1x))
c$$$      close ( 33 ) 
c$$$      
c$$$      fname_gas = path(1:lpath)//'//fion'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname), form = 'unformatted' )
c$$$      write(20) aexpn, nsvx, nsvy, nsvz
c$$$      write(20) dx * r0 , dy * r0 , dz * r0
c$$$      write(20) sv
c$$$      close ( 20 )
c$$$      fname_gas = path(1:lpath)//'//fion'//'_a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'L187.5kpc_ascii.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      write(20,*) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      return
c$$$      end
c$$$c
c$$$c     -----------------------------------
c$$$      subroutine DataOutput3 ( tt , dtt )
c$$$c     -----------------------------------
c$$$c
c$$$c     ViewDump for debugging the AGN feedback
c$$$c     The slice is through the pixel of the AGN energy injection
c$$$c     
c$$$      common /initrun_COM/ xc, fpert, eps
c$$$c
c$$$      include 'a_tree.h'
c$$$      include 'a_control.h'
c$$$      include 'a_analysis.h'
c$$$      include 'a_agn.h'
c$$$
c$$$      character*5 fstep
c$$$      character*256 fname_gas
c$$$      logical lpointsout, l_force_center
c$$$      integer imovstep, inum
c$$$      SAVE idumpp
c$$$      DATA idumpp/0/
c$$$
c$$$c
c$$$c.... output profiles 
c$$$c
c$$$      xc = xAGN(1)
c$$$      yc = yAGN(1)
c$$$      zc = zAGN(1) 
c$$$      rmin = 1.e-3/r0
c$$$      rmax = 10.0/r0
c$$$      nrbin = 80 
c$$$      ibintype = 0 
c$$$      deltavir = 340
c$$$      deltacrit = 500
c$$$      ifindcenter = 3 
c$$$      lpointsout = .false.  
c$$$c      l_force_center = .true.  
c$$$      l_force_center = .false. 
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$
c$$$      nstep = 128
c$$$      nstep = 16
c$$$      idumpp = idumpp +1
c$$$      if ( lviewout .and. mod(idumpp,nstep) .eq. 0) then
c$$$        inum = idumpp / nstep
c$$$        call GetName ( inum , fstep )
c$$$        fname_gas  =
c$$$     &           path(1:lpath)//'/'//'pro'//'_'//fstep//'.dat '
c$$$        call Output_Average_Profiles_CSF ( xc , yc , zc , rvir,  
c$$$     &                              rcrit, rmin , rmax , nrbin , 
c$$$     &                              ibintype, deltavir, deltacrit, 
c$$$     &                              ifindcenter, lpointsout, 
c$$$     &                              l_force_center, fname_gas )
c$$$      endif
c$$$
c$$$c
c$$$c...  ViewDump at every 1 small steps
c$$$c      
c$$$      call Get_MaxLevelNow()
c$$$      if ( lviewout .and. mod(idumpp,nstep) .eq. 0) then 
c$$$        write(*,*) 'Slicing through the location of the AGN'
c$$$        write(*,10) xc, yc, zc
c$$$ 10     format ('  ( xc, yc, zc ) = ',3(g13.7,1x))
c$$$        cs = CellSize(MaxLevelNow) 
c$$$c        call ViewDump ( MaxLevelNow , xc+cs , 1 )
c$$$        call ViewDump ( MaxLevelNow , xc    , 1 )
c$$$c        call ViewDump ( MaxLevelNow , xc-cs , 1 )        
c$$$c        call ViewDump ( MaxLevelNow , yc , 2 )
c$$$c        call ViewDump ( MaxLevelNow , zc , 3 )
c$$$      endif
c$$$c
c$$$      return
c$$$      end
c$$$c
c$$$c     --------------------------------------------------
c$$$      subroutine DataOutput5 ( tt , dtt , aboxh, chbox )
c$$$c     --------------------------------------------------
c$$$c
c$$$      common /initrun_COM/ xc, fpert, eps
c$$$c
c$$$      include 'a_tree.h'
c$$$      include 'a_control.h'
c$$$      include 'a_analysis.h'
c$$$
c$$$      character*5 fstep
c$$$      real svd(nsvx,nsvy,nsvz)
c$$$      real p_d1(nsvy,nsvz), p_d2(nsvy,nsvz), p_d3(nsvy,nsvz)
c$$$      real p_t1(nsvy,nsvz), p_t2(nsvy,nsvz), p_t3(nsvy,nsvz)
c$$$      real*8 sum1, sum2, sum3, sum4, vmax1, vmax2, em
c$$$      common / SVDUM / svd, p_d1, p_d2, p_d3, p_t1, p_t2, p_t3
c$$$      character*256 fname_gas
c$$$      logical lpointsout, l_force_center
c$$$      integer imovstep
c$$$      character digits(0:9)
c$$$      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
c$$$      character*2 chbox
c$$$
c$$$#ifdef PRESSUREFLOOR
c$$$      real*8 f_J, f_eJ
c$$$      parameter ( f_J = 10.  ) ! pressure floor parameter
c$$$#endif
c$$$
c$$$#ifdef PRESSUREFLOOR
c$$$      DO Level = MinL_Jeans , MaxLevelNow
c$$$c
c$$$c.... artificial pressure floor a-la Machacek et al. 2
c$$$        if ( Level .ge. MinL_Jeans ) then
c$$$          f_eJ = 0.47746 * f_J * aexpn * 
c$$$     &           CellSize(MaxLevel)**2 * (gamma - 1.0)
c$$$        else
c$$$          f_eJ = 0.0
c$$$        endif
c$$$        E_min = 300.0 / T_0 * aexpn**2 * hvar(1,i1) * 1.5
c$$$        if ( Level .eq. MinLevel ) then
c$$$          do i1 = 1 , ncell0
c$$$            E_th = f_eJ * hvar(1,i1) * hvar(1,i1) 
c$$$            hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
c$$$          enddo 
c$$$        else
c$$$          call Select_Cells ( Level , nLevel )
c$$$          do ic1 = 1 , nLevel
c$$$            icell = iSelect(ic1)
c$$$            do ic2 = 0 , 7
c$$$              i1 = icell + ic2
c$$$              if ( iOctCh(i1) .eq. nil ) then
c$$$                E_th = f_eJ * hvar(1,i1) * hvar(1,i1)
c$$$                hvar(8,i1) = max ( E_min , hvar(8,i1)-E_th )
c$$$              endif
c$$$            enddo
c$$$          enddo
c$$$        endif
c$$$      ENDDO
c$$$#endif
c$$$#ifndef PRESSUREFLOOR
c$$$      f_eJ = 0.0
c$$$#endif
c$$$
c$$$c
c$$$c.... output profiles 
c$$$c
c$$$      xc = 16.53711
c$$$      yc = 22.15430
c$$$      zc = 128.7441 
c$$$      rmin = 1.e-3/r0
c$$$      rmax = 10.0/r0
c$$$      nrbin = 45 
c$$$      ibintype = 0 
c$$$      deltavir = 340
c$$$      deltacrit = 500
c$$$      ifindcenter = 3 
c$$$      lpointsout = .false.  
c$$$      l_force_center = .false.  
c$$$
c$$$      lpath = index(path    , ' ') - 1
c$$$
c$$$        ic = int(aexpn*1000)
c$$$        i1 = ic / 1000 
c$$$        i2 = (ic - (ic/1000)*1000)/100 
c$$$        i3 = (ic - (ic/100)*100) / 10 
c$$$        i4 = (ic - (ic/10)*10)
c$$$
c$$$      fname_gas  = path(1:lpath)//'/pro'//'_'//'a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'_'//chbox//'.dat '
c$$$
c$$$
c$$$      call Output_Average_Profiles_Stars ( xc , yc , zc , rvir,  
c$$$     &                            rmin , rmax , nrbin , ibintype ,  
c$$$     &                            deltavir, deltacrit, ifindcenter , 
c$$$     &                            lpointsout, l_force_center,
c$$$     &                            fname_gas ) 
c$$$
c$$$c      return
c$$$
c$$$c
c$$$c.... extract boxes 
c$$$c
c$$$      Posx = xc
c$$$      Posy = yc
c$$$      Posz = zc
c$$$      sizebox2 = aboxh/2.
c$$$      
c$$$      dx = 2.* sizebox2 / r0 / float(nsvx)
c$$$      dy = 2.* sizebox2 / r0 / float(nsvy)
c$$$      dz = 2.* sizebox2 / r0 / float(nsvz)
c$$$      xgmin = Posx - sizebox2/r0
c$$$      xgmax = Posx + sizebox2/r0
c$$$      ygmin = Posy - sizebox2/r0
c$$$      ygmax = Posy + sizebox2/r0
c$$$      zgmin = Posz - sizebox2/r0
c$$$      zgmax = Posz + sizebox2/r0
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax , 
c$$$     &                      ygmin , ygmax , 
c$$$     &                      zgmin , zgmax , 1 , 1 ) ! gas density 
c$$$
c$$$      rhofact = rho0 / aexpn**3 ! Msun / Mpc^3
c$$$      write(*,*) 'rhofact =',rhofact
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            svd(ix,iy,iz) = sv(ix,iy,iz)
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) * rhofact 
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$
c$$$
c$$$      write(*,*) 'extracting pressure...'      
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 8 )  ! pressure 
c$$$
c$$$      tfact = 2./3. * T_0 / aexpn**2 
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * tfact       
c$$$            svd(ix,iy,iz) = svd(ix,iy,iz) * rhofact
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      fname_gas  = path(1:lpath)//'/rhog'//'_'//'a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'_'//chbox//'.dat '
c$$$
c$$$      
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      write(20,*) svd
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      fname_gas  = path(1:lpath)//'/tg'//'_'//'a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'_'//chbox//'.dat '
c$$$
c$$$      
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      write(20,*) sv
c$$$      close ( 20 )
c$$$      write(*,*) 'done.'
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 3 )  ! vx
c$$$
c$$$       
c$$$      vfact = v0 / aexpn * rhofact
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vfact       
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$
c$$$      fname_gas  = path(1:lpath)//'/v'//'_'//'a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'_'//chbox//'.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      write(20,*) sv
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 4 )  ! vx
c$$$
c$$$       
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vfact       
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      write(20,*) sv
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 5 )  ! vx
c$$$
c$$$       
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vfact       
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      write(20,*) sv
c$$$
c$$$      close (20)
c$$$
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 9 )  ! ZII
c$$$
c$$$       
c$$$      vfact = rhofact / 0.0198
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vfact 
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$
c$$$      fname_gas  = path(1:lpath)//'/z'//'_'//'a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'_'//chbox//'.dat '
c$$$      nfname = index ( fname_gas , ' ' ) - 1
c$$$      write(*,*) 'writing...'
c$$$      open ( 20 , file = fname_gas(1:nfname) )
c$$$      write(20,*) aexpn, nsvx, nsvy, nsvz
c$$$      write(20,*) dx * r0 , dy * r0 , dz * r0
c$$$      write(20,*) sv
c$$$
c$$$      call Extract_Volume ( xgmin , xgmax ,
c$$$     &                      ygmin , ygmax ,
c$$$     &                      zgmin , zgmax , 1 , 10 )  ! ZIa
c$$$
c$$$       
c$$$      do iz = 1 , nsvz
c$$$        do iy = 1 , nsvy
c$$$          do ix = 1 , nsvx
c$$$            sv(ix,iy,iz) = sv(ix,iy,iz) / svd(ix,iy,iz) * vfact      
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      write(20,*) sv
c$$$      close(20)
c$$$
c$$$      fname_gas  = path(1:lpath)//'/st'//'_'//'a'//
c$$$     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$     &            //'.dat '
c$$$      call Write_Stars_ASCII ( fname_gas, xc,yc,zc )
c$$$
c$$$c      fname_gas  = path(1:lpath)//'/dm'//'_'//'a'//
c$$$c     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$c     &            //'.dat '
c$$$
c$$$c      call Write_DM_ASCII ( fname_gas, xc,yc,zc )
c$$$
c$$$
c$$$      return
c$$$      end
c
c     ---------------------------------------------
      subroutine Extract_Particles ( fname, 
     &                               xgmin, xgmax,
     &                               ygmin, ygmax,
     &                               zgmin, zgmax )
c     ---------------------------------------------
c
c     extract particles from the rectangular volume 
c     and write positions, velocities, and masses into ASCII file fname
c
c     fname should end with a space and not contain any other space characters
c
      include 'a_tree.h'
      include 'a_control.h'
      character*256 fname 
      integer lfn 

      lfn = index ( fname , ' ' ) - 1

      open ( 20 , file = fname(1:lfn) )
      write(20,*) 'xmin, xmax, ymin, ymax, zmin, zmax [/h Mpc]'
      write(20,*) (xgmin-1.)*r0, (xgmax-1.)*r0, 
     &            (ygmin-1.)*r0, (ygmax-1.)*r0,
     &            (zgmin-1.)*r0, (zgmax-1.)*r0
      write(20,*) 'aexpn =',aexpn
      vfact = v0 / aexpn
      do i = 1 , nspec
         wpart = am0 * wpar(i)
         do ip = nsp(i,1) , nsp(i,2)
          if ( x(ip) .ge. xgmin .and. x(ip) .le. xgmax ) then
          if ( y(ip) .ge. ygmin .and. y(ip) .le. ygmax ) then
          if ( z(ip) .ge. zgmin .and. z(ip) .le. zgmax ) then
            xd = (x(ip) - 1.0)*r0 
            yd = (y(ip) - 1.0)*r0
            zd = (z(ip) - 1.0)*r0
            vxd = vx(ip) * vfact
            vyd = vy(ip) * vfact
            vzd = vz(ip) * vfact
            write(20,90) xd, yd, zd, vxd, vyd, vzd, wpart
          endif
          endif
          endif
        enddo
      enddo
 90   format ( 3(f9.6,1x),3(f9.2,1x),g11.5)
      close ( 20 )

      return
      end

c     --------------------------------------------
      subroutine Extract_Volume ( xgmin , xgmax , 
     &                            ygmin , ygmax , 
     &                            zgmin , zgmax , 
     &                            ivartype , ivar )
c     ---------------------------------------------
c
c     extract variable ivar
c          [ivartype = 0 -> from var(ivar,*)
c           ivartype = 1 -> hvar(ivar,*)]
c     on a uniform rectangular grid sv covering sub-volume,
c     defined by 
c           xgmin, xgmax, ygmin, ygmax, zgmin, zgmax 
c
c     The grid size [nsvx, nsvy, nsvz] and array sv are 
c     specified in a_analysis.h
c
c     it is recommended that edges of the subvolume are aligned with 
c     zero level cells [i.e., xmin = 1.0, 2.0, 3.0, ...]
c     and size of the grid nsvx, nsvy, nsvz are powers of 2  
c
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_analysis.h'
c
      real*8 xd, yd, zd, xx, yy, zz 
c
      
      if ( xgmin .gt. xgmax ) then 
        write(*,*) 'error in Extract_Volume: xgmin > xgmax :',
     &             xgmin, xgmax
        write(*,*) 'exiting..'
        return
      endif
      if ( ygmin .gt. ygmax ) then
        write(*,*) 'error in Extract_Volume: ygmin > ygmax :',
     &             ygmin, ygmax
        write(*,*) 'exiting..'
        return
      endif
      if ( zgmin .gt. zgmax ) then
        write(*,*) 'error in Extract_Volume: zgmin > zgmax :',
     &             zgmin, zgmax
        write(*,*) 'exiting..'
        return
      endif

      if ( ivartype .eq. 0 .or. ivartype .eq. 1 ) then 
        if ( ivartype .eq. 0 ) then 
          if ( ivar .lt. 1 .or. ivar .gt. nvar ) then 
            write(*,*) 'error in Extract_Volume: wrong ivar:',ivar
            write(*,*) 'for ivartype=0 must be 1 <= ivar <= nvar =',
     &                  nvar
            write(*,*) 'exiting...'
            return
          endif 
        else
          if ( ivar .lt. 1 .or. ivar .gt. nhvar ) then 
            write(*,*) 'error in Extract_Volume: wrong ivar:',ivar
            write(*,*) 'for ivartype=1, must be 1 <= ivar <= nhvar =',
     &                  nhvar
            write(*,*) 'exiting...'
            return
          endif          
        endif
      else
        write(*,*) 'error in Extract_Volume: wrong ivarype:',
     &        ivartype
        write(*,*) 'but must be 0 or 1, exiting...'
        return
      endif
c
      dxsv = (xgmax - xgmin) / float(nsvx)
      dysv = (ygmax - ygmin) / float(nsvy)
      dzsv = (zgmax - zgmin) / float(nsvz)
      dxsv2 = 0.5 * dxsv
      dysv2 = 0.5 * dysv
      dzsv2 = 0.5 * dzsv
      dxmax = max ( dxsv , dysv , dzsv ) 
c
c.... effective level of the grid 
c
      iGL = nearest_int ( -log10(dxmax) / log10(2.0) )

      write(*,*) 'extracting volume...', 
     &       xgmin,xgmax,ygmin,ygmax,zgmin,zgmax
      write(*,*) 'effective resolution level ',iGL
      write(*,*) '=> comoving cell size of',r0/2**iGL

      zd = zgmin - dxsv2

      do iz = 1 , nsvz
        zd = zd + dxsv
	yd = ygmin - dysv2
        zz = zd
        if ( zz .lt. 1.0 ) then
          zz = zz + ng
        else
          if ( zz .gt. narr ) then
            zz = zz - ng
          endif
        endif

        do iy = 1 , nsvy
	  yd = yd + dysv
          yy = yd
          if ( yy .lt. 1.0 ) then
            yy = yy + ng
          else
            if ( yy .gt. narr ) then
              yy = yy - ng
            endif
          endif

          xd = xgmin - dxsv2
          do ix = 1 , nsvx
	    xd = xd + dxsv
            xx = xd
            if ( xx .lt. 1.0 ) then
              xx = xx + ng
            else
              if ( xx .gt. narr ) then
                xx = xx - ng
              endif
            endif

            icell = iFindCell ( iGL , xx, yy, zz )
            if ( ivartype .eq. 0 ) then 
               cellvolume = CellVol(iLv(icell))
               if ( ivar .eq. 1 ) then
                  sv(ix,iy,iz) = var(ivar,icell) + cellvolume
                  sv(ix,iy,iz) = sv(ix,iy,iz) / cellvolume ! total density
               else 
                  sv(ix,iy,iz) = ref(icell) / cellvolume ! DM density
               endif
            else
              sv(ix,iy,iz) = hvar(ivar,icell)
            endif
          enddo
	enddo
      enddo

      return
      end


