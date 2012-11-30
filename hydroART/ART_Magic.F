c     =====================================================================
c                                                                         .
c          ART Version 1-0 : magic on hydro variables                     .
c                                                                         .
c      by Alexei Khokhlov and Andrey Kravtsov (1996-1999)                 .
c                                                                         .   
c     =====================================================================
c
c     --------------------------
      subroutine Magic ( Level ) 
c     --------------------------
c
c     rhohdr, enghdr - are user-supplied floor values for density 
c                      and energy, read in in ART_Init.f
c
      include 'a_tree.h'
      include 'a_control.h'
c
# include "a_def.h"
c
      integer Level 
      integer i1, i2, icell, idcell, nLevel, Nbb(6), i3
      real*8  ekinetic, v_x, v_y, v_z, E_th, T_minc
CEVERINO02222012
      real*8 T_maxc
CEVERINO01192012
      parameter (VeloMax=1000.) ! km/s
      real*8 VeloMax_code_Mo
c
      T_minc = T_min * aexp(Level)**2 / T_0 / (gamma - 1.0)
c
CEVERINO02222012
      T_maxc = 1d8 * aexp(Level)**2 / T_0 / (gamma - 1.0)
CEVERINO01192012
      VeloMax_code=VeloMax/v0 *aexp(Level) 


      IF ( Level .eq. MinLevel ) THEN

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(v_x,v_y,v_z,v_tot,ekinetic,E_th,rhoav,Nbb,i3)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then   ! only leaves
c==========STG:
             if ( hvar(1,i1) .le. 0.0 ) then
#ifdef WORKOUT
                call Open_ASCII_File ( iOWork , workfile , append )
                write(iOWork,*) '>>>>>>>>>>> negative density!'
                write(iOWork,*) '-----> density:', 
     &               istep, i1, hvar(1,i1), 
     &               hvar(2,i1), hvar(8,i1), hvar(6,i1),
     &               hvar(3,i1), hvar(4,i1), hvar(5,i1),
     &               iExitCode2
                close(iOWork)
#endif
                iExitCode2 = -1
                call Open_ASCII_File ( 98 , 'dens.log' , append )
                write(98,*) 'negative density:', 
     &               istep, i1, hvar(1,i1), 
     &               hvar(2,i1), hvar(8,i1), hvar(6,i1),
     &               hvar(3,i1), hvar(4,i1), hvar(5,i1),
     &               iExitCode2
                close(98)
c                Continue_Run = .false.
             endif
c==============
            if ( hvar(1,i1) .lt. rhohdr ) then   
              rhoav = 0.
              call NbAll0 ( i1 , Nbb )
              do i3 = 1 , 6
                rhoav = rhoav + hvar(1,Nbb(i3)) 
              enddo
              hvar(1,i1) = max ( 0.1666666*rhoav , rhohdr )
              v_x   = hvar(3,i1)
              v_y   = hvar(4,i1)
              v_z   = hvar(5,i1)
c========STG:
CEVERINO01192012: CUT VELOCITIES LARGER THAN VeloMax KM/S
              VeloMax_code_Mo=VeloMax_code*hvar(1,i1)
              v_x   = min (v_x , VeloMax_code_Mo) 
              v_x   = max (v_x , -VeloMax_code_Mo)
              v_y   = min (v_y , VeloMax_code_Mo)    
              v_y   = max (v_y , -VeloMax_code_Mo)         
              v_z   = min (v_z , VeloMax_code_Mo)          
              v_z   = max (v_z , -VeloMax_code_Mo)         
              hvar(3,i1)=v_x                           
              hvar(4,i1)=v_y                           
              hvar(5,i1)=v_z                           
              ekinetic =                                   
     &             0.5d0*(v_x**2 + v_y**2 + v_z**2) / hvar(1,i1)
              E_th  = T_minc * hvar(1,i1)
              hvar(8,i1) = max(hvar(8,i1) , sngl(E_th) )
CEVERINO02222012
              E_th  = T_maxc * hvar(1,i1)
              hvar(8,i1) = min(hvar(8,i1) , sngl(E_th) )
              hvar(2,i1) =  hvar(8,i1) +  sngl(ekinetic)    
c===============              
            else 
              v_x   = hvar(3,i1)
              v_y   = hvar(4,i1)
              v_z   = hvar(5,i1)
CEVERINO01192012: CUT VELOCITIES LARGER THAN VeloMax KM/S
              VeloMax_code_Mo=VeloMax_code*hvar(1,i1)
              v_x   = min (v_x , VeloMax_code_Mo) 
              v_x   = max (v_x , -VeloMax_code_Mo)
              v_y   = min (v_y , VeloMax_code_Mo)    
              v_y   = max (v_y , -VeloMax_code_Mo)         
              v_z   = min (v_z , VeloMax_code_Mo)          
              v_z   = max (v_z , -VeloMax_code_Mo)         
              hvar(3,i1)=v_x                           
              hvar(4,i1)=v_y                           
              hvar(5,i1)=v_z
              ekinetic =                                   
     &             0.5d0*(v_x**2 + v_y**2 + v_z**2) / hvar(1,i1)       
              E_th  = T_minc * hvar(1,i1)
              hvar(8,i1) = max(hvar(8,i1) , sngl(E_th) )
CEVERINO02222012
              E_th  = T_maxc * hvar(1,i1)
              hvar(8,i1) = min(hvar(8,i1) , sngl(E_th) )
c--------------------------------------------------------------------
              hvar(2,i1) =  hvar(8,i1) +  sngl(ekinetic)
c-------------------------------------------------------------------- 
c===============
           endif
#ifdef ADVECT_SPECIES
            do j = 9 , nhvar 
              hvar(j,i1) = max ( 1.e-15 , hvar(j,i1) )
            enddo
#endif            
          endif
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,idcell,i2)
C$OMP+PRIVATE(v_x,v_y,v_z,v_tot,ekinetic,E_th,rhoav,Nbb,i3)
C$OMP+PRIVATE(VeloMax_code_Mo)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            if ( iOctCh(idcell) .eq. nil ) then   ! only leaves
c==========STG:
               if ( hvar(1,idcell) .le. 0.0 ) then
#ifdef WORKOUT
                  call Open_ASCII_File ( iOWork , workfile , append )
                  write(iOWork,*) '>>>>>>>>>>>>>>> negative density!'
                  write(iOWork,*) '-----> density:', 
     &                 idcell, hvar(1,idcell), 
     &                 hvar(2,idcell), hvar(8,idcell), hvar(6,idcell),
     &                 hvar(3,idcell), hvar(4,idcell), hvar(5,idcell),
     &                 istep
                  close(iOWork)
#endif
                  iExitCode2 = -1
                  call Open_ASCII_File ( 98 , 'dens.log' , append )
                  write(98,*) 'negative density:', 
     &               istep, i1, hvar(1,i1), 
     &               hvar(2,i1), hvar(8,i1), hvar(6,i1),
     &               hvar(3,i1), hvar(4,i1), hvar(5,i1),
     &               iExitCode2
                  close(98)
              endif
c==============
              if ( hvar(1,idcell) .lt. rhohdr ) then 
                rhoav = 0.
                call NbAllL ( idcell , Nbb )
                do i3 = 1 , 6
                  rhoav = rhoav + hvar(1,Nbb(i3))
                enddo
CEVERINO02212012
check
                write (128,*) level, idcell, hvar(1,idcell),
     &                        0.16666*rhoav,rhohdr
c-------------------------------------------------------------------- 
                hvar(1,idcell) = max ( 0.1666666*rhoav , rhohdr )
                v_x   = hvar(3,idcell)
                v_y   = hvar(4,idcell)
                v_z   = hvar(5,idcell)
CEVERINO01192012: CUT VELOCITIES LARGER THAN VeloMax KM/S
                VeloMax_code_Mo=VeloMax_code*hvar(1,idcell)
                v_x   = min (v_x , VeloMax_code_Mo) 
                v_x   = max (v_x , -VeloMax_code_Mo)
                v_y   = min (v_y , VeloMax_code_Mo)    
                v_y   = max (v_y , -VeloMax_code_Mo)         
                v_z   = min (v_z , VeloMax_code_Mo)          
                v_z   = max (v_z , -VeloMax_code_Mo)         
                hvar(3,idcell)=v_x                           
                hvar(4,idcell)=v_y                           
                hvar(5,idcell)=v_z                           
                ekinetic =                                   
     &             0.5d0*(v_x**2 + v_y**2 + v_z**2) / hvar(1,idcell)
                E_th  = T_minc * hvar(1,idcell)
                hvar(8,idcell) = max(hvar(8,idcell) , sngl(E_th) )
CEVERINO02222012
                E_th  = T_maxc * hvar(1,idcell)
                hvar(8,idcell) = min(hvar(8,idcell) , sngl(E_th) )
c-------------------------------------------------------------------- 
                hvar(2,idcell) =  hvar(8,idcell) +  sngl(ekinetic)    
              else 
                v_x   = hvar(3,idcell)
                v_y   = hvar(4,idcell)
                v_z   = hvar(5,idcell)
CEVERINO01192012: CUT VELOCITIES LARGER THAN VeloMax KM/S
                VeloMax_code_Mo=VeloMax_code*hvar(1,idcell)
                v_x   = min (v_x , VeloMax_code_Mo) 
                v_x   = max (v_x , -VeloMax_code_Mo)
                v_y   = min (v_y , VeloMax_code_Mo)    
                v_y   = max (v_y , -VeloMax_code_Mo)         
                v_z   = min (v_z , VeloMax_code_Mo)          
                v_z   = max (v_z , -VeloMax_code_Mo)         
                hvar(3,idcell)=v_x                           
                hvar(4,idcell)=v_y                           
                hvar(5,idcell)=v_z
                ekinetic =                                   
     &             0.5d0*(v_x**2 + v_y**2 + v_z**2) / hvar(1,idcell)       
                E_th  = T_minc * hvar(1,idcell)
                hvar(8,idcell) = max(hvar(8,idcell) , sngl(E_th) )
CEVERINO02222012
                E_th  = T_maxc * hvar(1,idcell)
                hvar(8,idcell) = min(hvar(8,idcell) , sngl(E_th) )
c--------------------------------------------------------------------
                hvar(2,idcell) =  hvar(8,idcell) +  sngl(ekinetic)
c-------------------------------------------------------------------- 
             endif
#ifdef ADVECT_SPECIES
              do j = 9 , nhvar 
                hvar(j,idcell) = max ( 1.e-15 , hvar(j,idcell) )
              enddo
#endif            
            endif
          enddo
        enddo
      ENDIF
c
      return
      end

c     ---------------------------
      subroutine FILTER ( Level )
c     ----------------------------
c
      include 'a_tree.h' 
      include 'a_control.h'
      parameter ( w0 = 100.0, w1 = 1.0, w01 = 1./106. )
      real*8 wsum1,wsum2,wsum3,wsum4,wsum5,wsum6 
      integer Nbb(6)

      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i2,Nbb,i3)
C$OMP+PRIVATE(wsum1,wsum2,wsum3,wsum4,wsum5,wsum6)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then   ! only leaves
            call NbAll0 ( i1 , Nbb ) 
            wsum1 = w0 * hvar(1,i1)
            wsum2 = w0 * hvar(2,i1)
            wsum3 = w0 * hvar(3,i1)
            wsum4 = w0 * hvar(4,i1)
            wsum5 = w0 * hvar(5,i1)
            wsum6 = w0 * hvar(8,i1)
            do i2 = 1 , 6
              wsum1 = wsum1 + w1 * hvar(1,Nbb(i2))
              wsum2 = wsum2 + w1 * hvar(2,Nbb(i2))
              wsum3 = wsum3 + w1 * hvar(3,Nbb(i2))
              wsum4 = wsum4 + w1 * hvar(4,Nbb(i2))
              wsum5 = wsum5 + w1 * hvar(5,Nbb(i2))
              wsum6 = wsum6 + w1 * hvar(8,Nbb(i2))
            enddo
            hvar(1,i1) = wsum1 * w01
            hvar(2,i1) = wsum2 * w01
            hvar(3,i1) = wsum3 * w01
            hvar(4,i1) = wsum4 * w01
            hvar(5,i1) = wsum5 * w01
            hvar(8,i1) = wsum6 * w01
          endif
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel )
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,idcell,i2,Nbb,i3)
C$OMP+PRIVATE(wsum1,wsum2,wsum3,wsum4,wsum5,wsum6)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            if ( iOctCh(idcell) .eq. nil ) then   ! only leaves
              call NbAllL ( idcell , Nbb ) 
              wsum1 = w0 * hvar(1,idcell)
              wsum2 = w0 * hvar(2,idcell)
              wsum3 = w0 * hvar(3,idcell)
              wsum4 = w0 * hvar(4,idcell)
              wsum5 = w0 * hvar(5,idcell)
              wsum6 = w0 * hvar(8,idcell)
              do i3 = 1 , 6
                wsum1 = wsum1 + w1 * hvar(1,Nbb(i3))
                wsum2 = wsum2 + w1 * hvar(2,Nbb(i3))
                wsum3 = wsum3 + w1 * hvar(3,Nbb(i3))
                wsum4 = wsum4 + w1 * hvar(4,Nbb(i3))
                wsum5 = wsum5 + w1 * hvar(5,Nbb(i3))
                wsum6 = wsum6 + w1 * hvar(8,Nbb(i3))
              enddo
              hvar(1,idcell) = wsum1 * w01
              hvar(2,idcell) = wsum2 * w01
              hvar(3,idcell) = wsum3 * w01
              hvar(4,idcell) = wsum4 * w01
              hvar(5,idcell) = wsum5 * w01
              hvar(8,idcell) = wsum6 * w01
            endif
          enddo
        enddo
      ENDIF

      return
      end

c     --------------------------------------
      subroutine Smooth_Potential  ( Level )
c     --------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
      parameter ( w0 = 1.0, w1 = 1.0, w01 = 1./(w0 + 6.*w1) )
      real*8 wsum1,wsum2,wsum3,wsum4,wsum5,wsum6
      integer Nbb(6)

      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i2,Nbb,i3)
C$OMP+PRIVATE(wsum1,wsum2,wsum3,wsum4,wsum5,wsum6)
        do i1 = 1 , ncell0
          call NbAll0 ( i1 , Nbb )
          wsum1 = w0 * var(2,i1)
          do i2 = 1 , 6
            wsum1 = wsum1 + w1 * var(2,Nbb(i2))
          enddo
          vnw(1,i1) = wsum1 * w01
        enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
        do i1 = 1 , ncell0
          var(2,i1) = vnw(1,i1)
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel )
       do i = 1 , 1 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,idcell,i2,Nbb,i3)
C$OMP+PRIVATE(wsum1,wsum2,wsum3,wsum4,wsum5,wsum6)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
              call NbAllL ( idcell , Nbb )
              wsum1 = w0 * var(2,idcell)
              do i3 = 1 , 6
                wsum1 = wsum1 + w1 * var(2,Nbb(i3))
              enddo
              vnw(1,idcell) = wsum1 * w01
          enddo
        enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,idcell,i2)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            var(2,idcell) = vnw(1,idcell)
          enddo
        enddo
       enddo
      ENDIF

      return
      end
