c     =====================================================================
c                                                                         
c     ART : source steps for some of the gas equations 
c                                                                         
c     by Andrey Kravtsov (2001)    
c                                                                         
c     =====================================================================
c

c     -------------------------------------------
      subroutine Advance_InternalEnergy ( Level )
c     -------------------------------------------
c
c     source step  for internal energy 
c
      integer Level
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      real*8  esmall, Zsolar, Zsolari
      parameter ( Zsolar = 0.02 , Zsolari = 1.d0 / Zsolar )
      parameter ( eta1 = 0.1, eta2 = 1.0e-3 ) 
      parameter ( dstep = 0.01 , esmall = 0. )
      parameter ( gammin = 1.01 )        
      real*8 b2a, t_e , dt_e, dtcool, dtcoolmin, dt3, ei1
      real*8 tbegin, tend, e_kin, ei, ei2, eratio, a1, a2, afact, acurr
      real*8 T_gas, ai, coolrate , coolcomp , rhog2, rhogi 
      real*8 dummy0, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6
      real*8 Tfact, Tfact1, Tfact2, rhogl, fact_nH
      real*8 e_pdV, e_cool, e_flux, e_small, e_min, Zdum 

      logical continue_cool

      integer Nbb(6)

      dt3 = dtl(Level) / 3.0
      tbegin = tl(Level)
      tend   = tl(Level) + dtl(Level)
#ifdef COSMOLOGY
      ai = 1.d0 / aexp(Level) 
      a1 = b2a ( dble(tl(Level)) )
      a2 = b2a ( dble(tl(Level) + dtl(Level)) )
      afact = (a2 - a1) / dtl(Level) ! stuff for linear interpolation ofa(t)
#endif
#ifndef COSMOLOGY
      a1 = 1.0
      a2 = 1.0 
      ai = 1.0
      afact = 0.
#endif

#ifndef METALCOOLING
      Zdum = -1000.0
#endif
#ifndef CLOUDY_COOLING
      rhogl = 0.0 
#endif
#ifdef CLOUDY_COOLING
c.... factor converting gas density in code units to n_H in cm^-3
      fact_nH = log10(1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3)
#endif

      gamma1 = gamma - 1.0
      Tfact = T_0 * gamma1 / 1.d4 ! T_0 must be in Kelvin
      dxi = 0.5 / CellSize(Level)
      e_min = T_min * a1**2 / T_0 / gamma1


      call Get_MaxLevelNow ()

      IF ( Level .eq. MinLevel ) THEN 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(t_e, dt_e, dt_e2, e_kin, ei, ei2, eratio, ei1, rhog2)
C$OMP+PRIVATE(acurr, ai, dummy0, dummy1, dummy2, dummy3, dummy4, dummy5)
C$OMP+PRIVATE(dummy6, rhogi, Tfact1, Tfact2, T_gas, coolrate, dtcool)
C$OMP+PRIVATE(e_small, Zdum, rhogl, Nbb)
C$OMP+PRIVATE(continue_cool)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves 
c
c....       now take care of p*div V term and cooling/heating 
c
            t_e   = tbegin
            dummy0 = hvar(7,i1) - 1.0
            dummy5 = dummy0 * ref(i1)

            do while ( t_e .lt. tend )
              dt_e   = min ( dt3 , tend - t_e )
              hvar(8,i1) = hvar(8,i1) * (1.0 + dummy5 * dt_e)
              t_e = t_e + dt_e
            enddo ! end do while

c
c....       sync internal and total energy
c
            eratio = hvar(2,i1)
            e_kin = 0.5d0 * ( hvar(3,i1)**2 +
     &                        hvar(4,i1)**2 +
     &                        hvar(5,i1)**2   ) / hvar(1,i1)
            ei     = (dble(hvar(2,i1)) - e_kin)
            eratio = ( ei / eratio )
            if ( eratio .gt. eta2 ) then ! use e = E - rho * v**2 /2
              hvar(8,i1) = ei
            endif ! endif eratio
            hvar(7,i1) = gamma


#ifdef COOLING
            rhog2 = (1.d0*hvar(1,i1)) * hvar(1,i1)
            rhogi = 1.d0 / hvar(1,i1) 
#ifdef CLOUDY_COOLING
c....       rhog_code -> log10(n_H), [n_H]=[cm^-3]
            rhogl = log10(hvar(1,i1)) + fact_nH
#endif
            t_e   = tbegin
            acurr = a1
            Tfact1 = Tfact * rhogi
            e_small = max ( esmall * hvar(8,i1) , 
     &                      e_min * hvar(1,i1)) 
            continue_cool = .true. 

#ifdef METALCOOLING            
#ifdef ENRICH
            Zdum = hvar(izII,i1) 
#ifdef ENRICH_SNIa
            Zdum = Zdum + hvar(izIa,i1)
#endif
            Zdum = max ( Zdum , 1.d-10 )
            Zdum = log10(Zdum * rhogi * Zsolari)
#endif
#endif

            do while ( (t_e .lt. tend) .and. continue_cool ) 


              ai = 1.d0 / acurr
              dummy2 = ai * ai 
              Tfact2 = dummy2 * Tfact1
              T_gas = Tfact2 * hvar(8,i1)  
              call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 

              dummy1 = -rhog2 * acurr
              coolrate = coolrate * dummy1
              dtcool = abs( hvar(8,i1) / coolrate ) 
              dtcool = dstep * dtcool
              dt_e   = min ( dtcool , tend - t_e )
c
              dt_e2  = dt_e * 0.5
              ei1 = max(hvar(8,i1) + coolrate * dt_e2, e_small)
              T_gas = Tfact2 * ei1
              call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 
              coolrate = coolrate * dummy1
              dummy3 = coolrate * dt_e
              hvar(8,i1) = hvar(8,i1) + dummy3
              hvar(2,i1) = hvar(2,i1) + dummy3
              if ( hvar(8,i1) .lt. e_small ) then 
                continue_cool = .false.
              endif
              hvar(8,i1) = max(dble(hvar(8,i1)), e_small)
              hvar(2,i1) = max(dble(hvar(2,i1)), e_small)
                      
              t_e = t_e + dt_e 
              acurr = a1 + afact * (t_e - tl(Level))
            enddo ! end do while

c            endif
#endif
          endif
          
        enddo ! end i1

      ELSE

        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,iChild,i1, ici, ici0, Nbb, ei2)
C$OMP+PRIVATE(t_e, dt_e, dt_e2, e_kin, ei, eratio, ei1)
C$OMP+PRIVATE(rhog2, rhogi, rhogl)
C$OMP+PRIVATE(acurr, ai, dummy0, dummy1, dummy2, dummy3, dummy4, dummy5)
C$OMP+PRIVATE(Tfact1, Tfact2, T_gas, coolrate, coolcomp, dtcool, dummy6)
C$OMP+PRIVATE(e_flux, e_pdV, e_cool, e_small, dtcoolmin, Zdum)
C$OMP+PRIVATE(continue_cool)
        do ic = 1 , nLevel
          icell = iSelect(ic) - 1
          ici0 = (ic-1)*nchild
          do iChild = 1 , 8 
            i1 = icell + iChild
            ici = ici0 + iChild
            if ( iOctCh(i1) .eq. nil ) then  ! only leaves 
c
c....         now take care of p*div V term and cooling/heating 
c
c            write(*,*) 'i1=',i1
              t_e   = tbegin
              dummy0 = hvar(7,i1) - 1.0
              dummy5 = dummy0 * ref(i1)

              do while ( t_e .lt. tend )
                dt_e   = min ( dt3 , tend - t_e )
                hvar(8,i1) = hvar(8,i1) * (1.0 + dummy5 * dt_e)
                t_e = t_e + dt_e
              enddo ! end do while

c
c....         sync internal and total energy
c
              eratio = hvar(2,i1)
              e_kin = 0.5d0 * ( hvar(3,i1)**2 +
     &                          hvar(4,i1)**2 +
     &                          hvar(5,i1)**2   ) / hvar(1,i1)
              ei     = (dble(hvar(2,i1)) - e_kin)
              eratio = ( ei / eratio )

              if ( eratio .gt. eta2 ) then ! use e = E - rho * v**2 /2
                hvar(8,i1) = ei
              endif ! endif eratio
              hvar(7,i1) = gamma


#ifdef COOLING
              rhog2 = (1.d0*hvar(1,i1)) * hvar(1,i1)
              rhogi = 1.d0 / hvar(1,i1) 
#ifdef CLOUDY_COOLING
c....         rhog_code -> log10(n_H), [n_H]=[cm^-3]
              rhogl = log10(hvar(1,i1)) + fact_nH
#endif
              t_e   = tbegin
              acurr = a1
              Tfact1 = Tfact * rhogi 
              e_small = max ( esmall * hvar(8,i1) , 
     &                        e_min * hvar(1,i1)) 
              continue_cool = .true. 

#ifdef METALCOOLING            
#ifdef ENRICH
              Zdum = hvar(izII,i1) ! SNII metal density
#ifdef ENRICH_SNIa
              Zdum = Zdum + hvar(izIa,i1)  ! SNIa metal density 
#endif
              Zdum = max ( Zdum , 1.d-10 )
              Zdum = log10 ( Zdum * rhogi * Zsolari ) 
#endif
#endif

c              write(*,*) i1, rhogi, Zdum, Zsolari

              do while ( (t_e .lt. tend) .and. continue_cool )
                ai = 1.d0 / acurr
                dummy2 = ai * ai 
                Tfact2 = dummy2 * Tfact1
                T_gas = Tfact2 * hvar(8,i1)  
                call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 

                dummy1 = -rhog2 * acurr
                coolrate = coolrate * dummy1
                dtcool = abs(hvar(8,i1) / coolrate) 
                dtcool = dstep * dtcool
c
                dt_e   = min ( dtcool , tend - t_e )
c
                dt_e2  = dt_e * 0.5
                ei1 = max(hvar(8,i1) + coolrate * dt_e2, e_small)
                T_gas = Tfact2 * ei1
                call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 
                coolrate = coolrate * dummy1
                dummy3 = coolrate * dt_e
                hvar(8,i1) = hvar(8,i1) + dummy3
                hvar(2,i1) = hvar(2,i1) + dummy3
                if ( hvar(8,i1) .lt. e_small ) then 
                  continue_cool = .false.
                endif
                hvar(8,i1) = max(dble(hvar(8,i1)), e_small)
                hvar(2,i1) = max(dble(hvar(2,i1)), e_small)

CEVERINO02222006
c                if ( Level .eq. MaxLevelNow ) then 
c                   ei     = (dble(hvar(2,i1)) - e_kin) * Tfact2*1.d4
c
c                   write(123,*) i1,
c     &            t_e*t0*acurr**2, tend*t0*acurr**2 ,dtcool*t0*acurr**2,
c     &                    rhogl,Zdum,coolrate,T_gas*1.d4,ei
c                endif

                t_e = t_e + dt_e 
                acurr = a1 + afact * (t_e - tl(Level))

              enddo ! end do while

c              endif
#endif

            endif ! end if iOctCh            
          enddo
        enddo
      ENDIF
      return
      end
c
c     --------------------------------------------
      subroutine GravWork ( iDir, dtLevel, Level )
c     --------------------------------------------
c
c     post hoc gravitation correction to gas energy and momenta
c     to be applied after StepH, MoveLevel, Poisson if no grav. correction
c     is done in StepH itself
c
c     Input:   Level      - level of leaves for which to compute grav. work
c              dtLevel    - timestep
c              iDir       - direction of gravity
c
c     it is assumed that new gravity is calculated at aexpn
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer iD , iDir1 , Level , nLevel , idcell , icell
      integer i1 , i2 , iNbl , iNbr , Nbb(6)
      real dtLevel , dt 
      real*8 acc, phi_lo , phi_ro , phi_ln , phi_rn , vnew
      real*8 factor1, factor2, rhor, rhor1, rho1, rho2
      real*8 b2a
      logical iLevelFlag
      integer iNbDir(3), iNbCh(8,3)
      data    iNbDir / 1 , 3 , 5 / 
      data    iNbCh  / 2 , 1 , 4 , 3 , 6 , 5 , 8 , 7 , 
     &                 3 , 4 , 1 , 2 , 7 , 8 , 5 , 6 , 
     &                 5 , 6 , 7 , 8 , 1 , 2 , 3 , 4   /
c
      real vasa 
      common / VASA1 / vasa
c
      iDir1  = iDir
      
      alev = aexp(Level)
      
c      factor1 = 2.5d-1 * dtLevel * alev**2 / CellSize(Level)
      factor1 = 2.5d-1 * dtLevel * aexp(Level)**2 / CellSize(Level)
      factor2 = 2.5d-1 * dtLevel * aexpold(Level)**2 / CellSize(Level)
      dtla2 = 1.0 / dtLevel / alev**2
c                
      IF ( Level .eq. MinLevel ) THEN 
c
c....   it is assumed here that there can't be cells of Level < MinLevel
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(Nbb,rho1,rho2,accx,accy,accz,rhor,rhor1,vnew)
          do i1 = 1 , ncell0
            if ( iOctCh(i1) .eq. nil ) then  ! only leaves 
              call NbAll0 ( i1 , Nbb )
c
              rho2 = factor2 * hvar(1,i1)
c              rho1 = factor1 * vnw(1,i1)
              rho1 = rho2
              accx = rho2 * (var(3,Nbb(2)) - var(3,Nbb(1))) +
     &               rho1 * (var(2,Nbb(2)) - var(2,Nbb(1)))
              accy = rho2 * (var(3,Nbb(4)) - var(3,Nbb(3))) +
     &               rho1 * (var(2,Nbb(4)) - var(2,Nbb(3)))
              accz = rho2 * (var(3,Nbb(6)) - var(3,Nbb(5))) +
     &               rho1 * (var(2,Nbb(6)) - var(2,Nbb(5)))

              if ( hvar(1,i1) .le. rhogrv ) then 
                accx = 0.
                accy = 0.
                accz = 0.
              endif

              rhor  = 1.0 / hvar(1,i1) 
              rhor1 = rhor * dtla2
c              gacc(1,i1) = accx * rhor1
c              gacc(2,i1) = accy * rhor1
c              gacc(3,i1) = accz * rhor1

              rhor = rhor * 0.5
              vnew = dble(hvar(3,i1)) - accx
              hvar(2,i1) = dble(hvar(2,i1)) + 
     &                    (vnew**2 - dble(hvar(3,i1))**2) * rhor
              hvar(3,i1) = vnew

              vnew = dble(hvar(4,i1)) - accy
              hvar(2,i1) = dble(hvar(2,i1)) +
     &                    (vnew**2 - dble(hvar(4,i1))**2) * rhor
              hvar(4,i1) = vnew

              vnew = dble(hvar(5,i1)) - accz
              hvar(2,i1) = dble(hvar(2,i1)) +
     &                    (vnew**2 - dble(hvar(5,i1))**2) * rhor
              hvar(5,i1) = vnew

            endif
          end do
c
      ELSE
c
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idcell,iNbl,iNbr,iChN)
C$OMP+PRIVATE(phi_ln,phi_lo,phi_rn,phi_ro)
C$OMP+PRIVATE(Nbb,rho1,rho2,accx,accy,accz,rhor,rhor1,vnew)
        do i1 = 1 , nLevel
          icell = iSelect(i1) - 1
          do i2 = 1 , 8 
            idcell = icell + i2
            if ( iOctCh(idcell) .eq. nil ) then  ! only leaves
              call NbAllL ( idcell , Nbb )  
              rho2 = factor2 * hvar(1,idcell)
c              rho1 = factor1 * vnw(1,idcell)
              rho1 = rho2

              iNbl  = Nbb(1) 
              iNbr  = Nbb(2)
              if ( iLv(iNbl).lt.Level  ) then 
                iChN = iNbCh ( i2 , 1 )
                phi_ln = Pyramide ( iNbl , iChN , 2 )
                phi_lo = Pyramide ( iNbl , iChN , 3 ) 
              else
                phi_ln = var(2,iNbl)
                phi_lo = var(3,iNbl)
              endif
              if ( iLv(iNbr).lt.Level ) then 
                iChN = iNbCh ( i2 , 1 )
                phi_rn = Pyramide ( iNbr , iChN , 2 )
                phi_ro = Pyramide ( iNbr , iChN , 3 ) 
              else
                phi_rn = var(2,iNbr)
                phi_ro = var(3,iNbr)
              endif
c
              accx = rho2 * (phi_ro - phi_lo) +
     &               rho1 * (phi_rn - phi_ln)

              iNbl  = Nbb(3)
              iNbr  = Nbb(4)
              if ( iLv(iNbl).lt.Level  ) then
                iChN = iNbCh ( i2 , 2 )
                phi_ln = Pyramide ( iNbl , iChN , 2 )
                phi_lo = Pyramide ( iNbl , iChN , 3 )
              else
                phi_ln = var(2,iNbl)
                phi_lo = var(3,iNbl)
              endif
              if ( iLv(iNbr).lt.Level ) then
                iChN = iNbCh ( i2 , 2 )
                phi_rn = Pyramide ( iNbr , iChN , 2 )
                phi_ro = Pyramide ( iNbr , iChN , 3 )
              else
                phi_rn = var(2,iNbr)
                phi_ro = var(3,iNbr)
              endif
c
              accy = rho2 * (phi_ro - phi_lo) +
     &               rho1 * (phi_rn - phi_ln)

              iNbl  = Nbb(5)
              iNbr  = Nbb(6)
              if ( iLv(iNbl).lt.Level  ) then
                iChN = iNbCh ( i2 , 3 )
                phi_ln = Pyramide ( iNbl , iChN , 2 )
                phi_lo = Pyramide ( iNbl , iChN , 3 )
              else
                phi_ln = var(2,iNbl)
                phi_lo = var(3,iNbl)
              endif
              if ( iLv(iNbr).lt.Level ) then
                iChN = iNbCh ( i2 , 3 )
                phi_rn = Pyramide ( iNbr , iChN , 2 )
                phi_ro = Pyramide ( iNbr , iChN , 3 )
              else
                phi_rn = var(2,iNbr)
                phi_ro = var(3,iNbr)
              endif
c
              accz = rho2 * (phi_ro - phi_lo) +
     &               rho1 * (phi_rn - phi_ln)

              if ( hvar(1,idcell) .le. rhogrv ) then 
                accx = 0.
                accy = 0.
                accz = 0.
              endif

              rhor  = 1.0 / hvar(1,idcell)
              rhor1 = rhor * dtla2
              rhor  = rhor * 0.5

c              gacc(1,idcell) = accx * rhor1
c              gacc(2,idcell) = accy * rhor1
c              gacc(3,idcell) = accz * rhor1

              vnew = dble(hvar(3,idcell)) - accx
              hvar(2,idcell) = dble(hvar(2,idcell)) + 
     &                        (vnew**2-dble(hvar(3,idcell))**2) * rhor
              hvar(3,idcell) = vnew

              vnew = dble(hvar(4,idcell)) - accy
              hvar(2,idcell) = dble(hvar(2,idcell)) +
     &                        (vnew**2-dble(hvar(4,idcell))**2) * rhor
              hvar(4,idcell) = vnew

              vnew = dble(hvar(5,idcell)) - accz
              hvar(2,idcell) = dble(hvar(2,idcell)) +
     &                        (vnew**2-dble(hvar(5,idcell))**2) * rhor
              hvar(5,idcell) = vnew

            endif
          enddo
        enddo

      ENDIF
c
      return
      end
c
