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
#     include "a_def.h"
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'

      real*8 dt
      real xmin,xmax,ymin,ymax,zmin,zmax
      real*8 tempc, coolrate
      character*256 fname, fname_gas
      character*256 fname_part1, fname_part2, fname_part3
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

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append ) 
         write(iOWork,*) 'reading particles from files'
         write(iOWork,'(A)') fname_part1
         write(iOWork,'(A)') fname_part2
         write(iOWork,'(A)') path
         close(iOWork)
#endif

         call Read_Particles_Binary ( fname_part1 , fname_part2 ) 

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append ) 
         write(iOWork,*) 'constructing linked list...'
         close(iOWork)
#endif
           call LL_Construct ()      

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append ) 
         write(iOWork,*) 'done.'
         close(iOWork)
#endif

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
#ifdef WORKOUT
c          
c.... compute each particle's weight
c
#ifdef STARFORM
           wsum = 0.
           do ic0 = 1 , nspecies-1
             do ic1 = nsp(ic0,1) , nsp(ic0,2)
               pw(ic1) = wpar(ic0)
             enddo
             wsum = wsum + wpar(ic0)*float(nsp(ic0,2)-nsp(ic0,1)+1)
           enddo
           wsum = wsum / ( 1. - Omb0 )
           call Open_ASCII_File ( iOWork , workfile , append ) 
           write(iOWork,*) 'isp =',ic0,'wsum =',wsum,' ncell0 =',ncell0
           close(iOWork)
#endif
#ifndef STARFORM
           wsum = 0.
           do ic0 = 1 , nspecies 
             do ic1 = nsp(ic0,1) , nsp(ic0,2)
               pw(ic1) = wpar(ic0)
             enddo
             wsum = wsum + wpar(ic0)*float(nsp(ic0,2)-nsp(ic0,1)+1)
           enddo
           wsum = wsum / ( 1. - Omb0 )
           call Open_ASCII_File ( iOWork , workfile , append ) 
           write(iOWork,*) 'isp =',ic0,'wsum =',wsum,' ncell0 =',ncell0
           close(iOWork)
#endif
#endif
#endif

#ifdef WORKOUT
           call Open_ASCII_File ( iOWork , workfile , append )
           write(iOWork,*) 'aexpn =',aexpn,' aexp0 =',aexp0,' astep =',astep
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
         read(50) ainit, astep
         read(50) ncells

#ifdef WORKOUT
         call Open_ASCII_File ( iOWork , workfile , append )
         write(iOWork,*) 'reading IC grid: ainit =',ainit,' ncells =',ncells
         if ( ncells .ne. ncell0 ) then 
           write(iOWork,*) 'bad ICs: ncells .ne. ncell0'
         endif
         close(iOWork)
#endif
         call SetUnits ()
         call Set_Cooling ()

         read(50) (hvar(1,ic),ic=1,ncell0)
         read(50) (hvar(3,ic),ic=1,ncell0)
         read(50) (hvar(4,ic),ic=1,ncell0)
         read(50) (hvar(5,ic),ic=1,ncell0)
         read(50) (hvar(2,ic),ic=1,ncell0)
         read(50) (hvar(8,ic),ic=1,ncell0)
         close(50)
c
         do ic = 1 , ncell0 
           hvar(7,ic) = gamma
         enddo
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

         call ChooseTimeStep ( dt )
c       
      ELSE

#ifdef HYDRO
          ljob  = index(jobname1, ' ') - 1
          lpath = index(path    , ' ') - 1
          fname  =
     &          path(1:lpath)//'/'//jobname1(1:ljob)//'.d '

          call Read_Gas_Binary ( dt , fname )
          call SetUnits ()
          call Set_Cooling ()        
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
     &                   path(1:lpath)//'/PMcrd.DAT '
          fname_part2 =
     &                   path(1:lpath)//'/PMcrs0.DAT '
          fname_part3 =
     &                   path(1:lpath)//'/pt.dat '

          call Read_Particles_Binary ( fname_part1 , fname_part2 )
          call Read_Particle_Time ( fname_part3 )  ! reads in dt for each particle
          call LL_Construct ()      

C$OMP PARALLEL DO 
          do ic1 = 1 , np 
            pt(ic1) = t   ! set t to the current time moment
          enddo
c
c.... compute or read each particle's weight
c
#ifdef STARFORM 
          do ic0 = 1 , nspecies-1
            do ic1 = nsp(ic0,1) , nsp(ic0,2)
              pw(ic1) = wpar(ic0)
            enddo
          enddo
c
c....     read weights of stellar particles 
c
          nstars = nsp(nspecies,2) - nsp(nspecies,1) 
          if ( nstars .gt. 0 ) then 
            open ( 20 , file = 'pw.dat' )
            read(20) nstf
            if ( nstf .ne. nstars ) then 
              write(*,*) '* error in Start: mismatch in number of stars'
              write(*,*) '* nstf .ne. nstars:',
     &          nstf, nstars
              stop
            endif
            read(20) (pw(i),i=nsp(nspecies,1),nsp(nspecies,2))
            close ( 20 )
          endif
#endif
#ifndef STARFORM 
          do ic0 = 1 , nspecies
            do ic1 = nsp(ic0,1) , nsp(ic0,2)
              pw(ic1) = wpar(ic0)
            enddo
          enddo
#endif
#endif

      ENDIF
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
      trho(MinLevel) = 2.0*wpar(1)
      trho(MinLevel+1) = 2.0*wpar(1)
      tbmass(MinLevel) = constb
      tbmass(MinLevel+1) = constb
      do Level = 2 , MaxLevel 
        trho(Level)  = 2.0*wpar(1)
        tbmass(Level) = constb / 8.0 
        write(*,*) 'L =',Level,' trho(L) =',trho(Level),
     &             ' trhog(L)=',tbmass(Level)
        tphi(Level) = -6000.0
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
      parameter ( nvars = 8 )
      character*10 varnames(nvars)
      real*8 ve
      real CellPos(3)
      integer iOD(3,2)
      data iOD / 2 , 1 , 1 , 
     &           3 , 3 , 2   / 
c
      include 'a_tree.h'
      include 'a_control.h'
c
      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1
      call GetName ( istep , fstep )
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
      write (unit=25,fmt='(a10,i10,1pe12.4,i10)')
     &            jobname1,imoviestep,aexp(Lev), nlayers
      write (unit=25,fmt='(''zones         1 1 1'')')
c
      nfuncs=8
      varnames(1)='Rho_g'
      varnames(2)='T'
      varnames(3)='v'
      varnames(4)='phi'
      varnames(5)='vx'
      varnames(6)='vy'
      varnames(7)='vz'
      varnames(8)='gamma'

      
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
      a3 = aexpn**3
      a2 = aexpn**2
      Tcons =  T_0 * (gamma - 1.0) / a2
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
                v_x = vfact * hvar(3,ic)/hvar(1,ic)
                v_y = vfact * hvar(4,ic)/hvar(1,ic)
                v_z = vfact * hvar(5,ic)/hvar(1,ic)
c                v_tot = sqrt(v_x**2+v_y**2+v_z**2)
cc                v_tot = v_tot * vcons / aexpn
c                g_x  = gacc(iOD(iDir,1),ic)
c                g_y  = gacc(iOD(iDir,2),ic)
c                g_z  = gacc(iDir,ic)
c                gtot = sqrt(g_x**2+g_y**2+g_z**2)
c                rho_dm = (var(1,ic) - hvar(1,ic) + 1.) 
         v_tot = sqrt(hvar(3,ic)**2+hvar(4,ic)**2+hvar(5,ic)**2)/
     &           hvar(1,ic) * vfact 
                E_1 = v_tot*0.5 + hvar(8,ic)
                e1 = hvar(2,ic) - v_tot*0.5
                rho_g  = hvar(1,ic) 
c                T_code = hvar(6,ic)/hvar(1,ic)
                T_code = hvar(8,ic)/hvar(1,ic)

                T_g = Tcons * T_code 
                S_g = Scons * T_code / hvar(1,ic)**(2./3.)
              if ( iDir .eq. 1 ) vdum = v_x
              if ( iDir .eq. 2 ) vdum = v_y
              if ( iDir .eq. 3 ) vdum = v_z
              fj = 1.0/(fjfact * sqrt(hvar(8,ic)/hvar(1,ic)**2))
              v_x = accc(ic,1)
              v_y = accc(ic,2)
              v_z = accc(ic,3)
                write (unit=25,fmt='(4i10,i5,1p9e14.6)')
     &                np1,np2,np3,np4,nzone,
     &                rho_g, T_g,v_tot, var(2,ic), 
     &                v_x, v_y, v_z, hvar(7,ic)
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
                  v_x = vfact * hvar(3,ic)/hvar(1,ic)
                  v_y = vfact * hvar(4,ic)/hvar(1,ic)
                  v_z = vfact * hvar(5,ic)/hvar(1,ic)
c                  v_tot = sqrt(v_x**2+v_y**2+v_z**2)
cc                  v_tot = v_tot * vcons / aexpn
c                  g_x  = gacc(iOD(iDir,1),ic)
c                  g_y  = gacc(iOD(iDir,2),ic)
c                  g_z  = gacc(iDir,ic)
c                  gtot = sqrt(g_x**2+g_y**2+g_z**2)
c                  rho_dm = var(1,ic) + cs3 - hvar(1,ic)*cs3
                  rho_g  = hvar(1,ic) 
c                T_code = hvar(6,ic)/hvar(1,ic)
                T_code = hvar(8,ic)/hvar(1,ic)

                T_g = Tcons * T_code
                S_g = Scons * T_code / hvar(1,ic)**(2./3.)

         v_tot = sqrt(hvar(3,ic)**2+hvar(4,ic)**2+hvar(5,ic)**2)/
     &           hvar(1,ic) * vfact
                E_1 = v_tot*0.5 + hvar(8,ic)
                E1 = hvar(2,ic) - v_tot*0.5
              if ( iDir .eq. 1 ) vdum = v_x
              if ( iDir .eq. 2 ) vdum = v_y
              if ( iDir .eq. 3 ) vdum = v_z
              fj = 1.0/(fjfact * sqrt(hvar(8,ic)/hvar(1,ic)**2))
              v_x = accc(ic,1)
              v_y = accc(ic,2)
              v_z = accc(ic,3)

                write (unit=25,fmt='(4i10,i5,1p9e14.6)')
     &                np1,np2,np3,np4,nzone,
     &                rho_g, T_g,v_tot, var(2,ic),
     &                v_x, v_y, v_z, hvar(7,ic) 

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
