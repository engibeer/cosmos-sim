c     ---------------------------------------
      subroutine FluxH ( dtx, dtx2, g, v, c, f )
c     ---------------------------------------

c     Calculates Euler fluxes across the 0-1 interface 
c
c          -1         0         1         2
c     -!----*----!----*----!----*----!----*----!-

c     Input      : dt               - time step 
c                  dx               - cell size
c                  dtx  = dt / dx 
c                  dtx2 = 0.5*dtx
c                  g(n,2)           - cell-centered accelerations for cells 0,1
c                                   = -a**2 * dphi/dx * 0.5 * dt
c                  c(n,0:1)         - slope interpolation correction factors
c                  v(n,nvar,-1:2)   - zone averaged values 
c                                     (.,1,.)     - density
c                                     (.,2,.)     - pressure
c                                     (.,3,.)     - longitudinal velocity
c                                     (.,4,.)     - transversal  velocity
c                                     (.,5,.)     - transversal  velocity
c                                     (.,6,.)     - gamma = P/E + 1
c                                     (.,7,.)     - Gamma = adiabatic index
c                                     (.,>7,.)    - chem. vars.

c     Output     : f(n,5)           - fluxes
c                                     (.,1)       - mass
c                                     (.,2)       - longitudinal momentum
c                                     (.,3)       - transversal momentum
c                                     (.,4)       - transversal momentum
c                                     (.,5)       - total energy density
c                                     (.,>5)      - species fluxes
c
#     include "a_def.h"

      include 'a_tree.h'
c
      real*8 dtx, dtx2
      real*8 v(nhvarA,-1:2)
      real*8 g(2), c(0:1)
      real*8 f(nhvarA-1)
c
      parameter ( gammin = 1.01, gammax = 10., smalll = 1.e-20 )
      real*8 dv(6,0:1)
      real*8  stl(5), str(5), sta(4)
c
      parameter ( epsP = 0.3, om1 = 0.75, om2 = 10. )
c
      real*8 dp0, dp1, dp, dvel, wj, r01, fl, dv0, dv1, dv2
      real*8 dv0a, dv1a, dv2a, dv11, dv20, dlq0, dlq1, dv00, dv01
      real*8 rhor_l, u_l, a_l, c_l, c2_l, cp_l, cm_l, x_l, rhow_l
      real*8 pw_l, uw_l, vw_l, ww_l, gamw_l, b_l, x0_l, rho0_l
      real*8 v0_l, w0_l, gam0_l, xx1, xx2, b0_l, v_l, w_l, p_l, gam_l
      real*8 rhor_r, u_r, a_r, c_r, c2_r, cp_r, cm_r, x_r, rhow_r
      real*8 pw_r, uw_r, vw_r, ww_r, gamw_r, b_r, x0_r, rho0_r
      real*8 v0_r, w0_r, gam0_r, xx3, xx4, b0_r, v_r, w_r, p_r, gam_r
      real*8 rho, vu, pre, gam, xup_r, xup_l, vv, vw, fmass, predtx
      real*8 pvu, vudtx, u
      real*8 fu, fv, fw, fe
c
c.... Flattening
c
      dp0 = v(2,1) - v(2,0)
      dp1 = v(2,2) - v(2,-1)
      dp  = abs ( dp0 ) / min ( v(2,1), v(2,0) )
      dvel= v(3,1) - v(3,0)
      wj  = ( onehalf + sign ( onehalf, dp - epsP ) ) *
     &      ( onehalf - sign ( onehalf, dvel ) )
      r01 = sign ( abs ( dp0 ) / ( abs ( dp1 ) + 1.e-30 ),
     &             dp0 * dp1 )
c      fl  = 1. -  wj * max ( zero, 
c     &                       min ( one, sngl(om2 * ( r01 - om1 )) ) )
      fl = 1.    ! ???
      stl(1) = fl * c(0)
      stl(2) = fl * c(1)
c
c.... Slopes
c
      do j = 1 , 6
        dv0  = v(j,0) - v(j,-1)
        dv1  = v(j,1) - v(j,0)
        dv2  = v(j,2) - v(j,1)
        dv0a = 2 * abs ( dv0 )
        dv1a = 2 * abs ( dv1 )
        dv2a = 2 * abs ( dv2 )
        dv11 = v(j,1) - v(j,-1) 
        dv20 = v(j,2) - v(j,0) 
        dlq0 = ( onehalf + sign ( onehalf, dv1*dv0 ) )
     &       * min ( dv0a , dv1a )
        dlq1 = ( onehalf + sign ( onehalf, dv1*dv2 ) ) 
     &       * min ( dv1a , dv2a )
        dv00 = sign ( min ( 0.5 * abs ( dv11 ), dlq0 ), dv11 )
     &       * stl(1)
        dv01 = sign ( min ( 0.5 * abs ( dv20 ), dlq1 ), dv20 )
     &       * stl(2)
        dv(j,0) = sign ( min ( abs ( dv00 ), dv1a, dv0a ), dv11 )
        dv(j,1) = sign ( min ( abs ( dv01 ), dv2a, dv1a ), dv20 )
      end do
c
c.... Left states
c          
      rhor_l   = 1. / v(1,0)
      u_l      = v(3,0) 
      a_l      = sqrt ( v(7,0) * v(2,0) * rhor_l )
      c_l      = a_l * v(1,0)
      c2_l     = c_l**2
      cp_l     = u_l + a_l
      cm_l     = u_l - a_l
      x_l      =  0.5 * ( 1. - dtx * max ( sngl(cp_l) , zero ) )
      rhow_l   = max ( smalll , sngl(v(1,0)     + x_l * dv(1,0)) )
      pw_l     = max ( smalll , sngl(v(2,0)     + x_l * dv(2,0)) )
      uw_l     = v(3,0)     + x_l * dv(3,0)
      vw_l     = v(4,0)     + x_l * dv(4,0)
      ww_l     = v(5,0)     + x_l * dv(5,0)
      gamw_l   = v(6,0)     + x_l * dv(6,0)
      b_l      = ( sign ( onehalf , cm_l ) + onehalf ) 
     &         * dtx2 * rhor_l * ( dv(2,0) / c_l - dv(3,0) )
      x0_l     = 0.5 * ( 1. - dtx * v(3,0) )
      rho0_l   = v(1,0) + x0_l * dv(1,0)
      v0_l     = v(4,0) + x0_l * dv(4,0)
      w0_l     = v(5,0) + x0_l * dv(5,0)
      gam0_l   = v(6,0) + x0_l * dv(6,0)
      if ( v(3,0) .gt. 0. ) then
        xx1 = 1.
      else
        xx1 = 0.
      end if
      xx2      = 1. - xx1
      b0_l     = xx1 * a_l * dtx2
     &               * (   dv(1,0) / ( rhow_l * rho0_l )
     &                   - dv(2,0) / c2_l                ) 
      v_l      = xx1 * v0_l + xx2 * vw_l
      w_l      = xx1 * w0_l + xx2 * ww_l
      dv(1,0)  = v_l
      dv(2,0)  = w_l
      stl(1)   = max ( smalll , 
     &                 sngl(rhow_l / ( 1. - ( b0_l + b_l ) * rhow_l )) )
      stl(2)   = uw_l - b_l * c_l
! + g(1) ! gravity correction
      p_l      = max ( smalll, sngl(pw_l + b_l * c2_l) )
      stl(3)   = p_l
      stl(4)   = v(7,0)
      gam_l    = xx1 * (   gam0_l
     &                   + 2. * ( 1. - v(6,0) / v(7,0) )
     &                        * ( v(6,0) - 1. )
     &                        * ( p_l - v(2,0) ) 
     &                        / ( p_l + v(2,0) )
     &                 ) 
     &         + xx2 * gamw_l
      stl(5)   = max ( gammin, min ( gammax , sngl(gam_l) ) )
c
c.... Right states
c
      rhor_r   = 1. / v(1,1)
      u_r      = v(3,1)
      a_r      = sqrt ( v(7,1) * v(2,1) * rhor_r )
      c_r      = a_r * v(1,1)
      c2_r     = c_r**2
      cm_r     = u_r - a_r
      cp_r     = u_r + a_r
      x_r      = -0.5 * ( 1. + dtx * min ( sngl(cm_r) , zero ) )
      rhow_r   = max ( smalll , sngl(v(1,1)     + x_r * dv(1,1)) )
      pw_r     = max ( smalll , sngl(v(2,1)     + x_r * dv(2,1)) )
      uw_r     = v(3,1)     + x_r * dv(3,1)
      vw_r     = v(4,1)     + x_r * dv(4,1)
      ww_r     = v(5,1)     + x_r * dv(5,1)
      gamw_r   = v(6,1)     + x_r * dv(6,1)
      b_r      = ( sign ( onehalf , cp_r ) - onehalf )
     &         * dtx2 * rhor_r * ( dv(3,1) + dv(2,1) / c_r )
      x0_r     =  - 0.5 * ( 1. + dtx*v(3,1) )
      rho0_r   = v(1,1) + x0_r * dv(1,1)
      v0_r     = v(4,1) + x0_r * dv(4,1)
      w0_r     = v(5,1) + x0_r * dv(5,1)
      gam0_r   = v(6,1) + x0_r * dv(6,1)
      if ( v(3,1) .lt. 0. ) then
        xx3    =  1.
      else
        xx3    =  0.
      end if
      xx4      = 1. - xx3
      b0_r     = xx3 * a_r * dtx2
     &               * (   dv(2,1) / c2_r
     &                   - dv(1,1) / ( rhow_r * rho0_r ) )
      v_r      = xx3 * v0_r + xx4 * vw_r
      w_r      = xx3 * w0_r + xx4 * ww_r
      dv(1,1)  = v_r
      dv(2,1)  = w_r
      str(1)   = max ( smalll,
     &                 sngl(rhow_r / ( 1. - ( b0_r + b_r ) * rhow_r ) ))
c     This is different from : uw_r + b_r * c_r, in Collela & Glaz (a typo?)
      str(2)   = uw_r - b_r * c_r
c  + g(2)  ! gravity correction (AVK)
      p_r      = max ( smalll, sngl(pw_r + b_r * c2_r) )
      str(3)   = p_r
      str(4)   = v(7,1)
      gam_r    = xx3 * (   gam0_r
     &                     + 2. * ( 1. - v(6,1) / v(7,1) )
     &                          * ( v(6,1) - 1. )
     &                          * ( p_r - v(2,1) ) 
     &                          / ( p_r + v(2,1) )
     &                   )
     &           + xx4 * gamw_r
      str(5)   = max ( gammin, min ( gammax , sngl(gam_r) ) )
c
c.... Solve Riemann problem
c
      call Riemann ( stl, str, sta )
c
c.... Fluxes of hydro variables
c
        v_l    = dv(1,0) 
        w_l    = dv(2,0) 
        v_r    = dv(1,1) 
        w_r    = dv(2,1) 
        rho    = sta(1)
        vu     = sta(2)
        pre    = sta(3)
        gam    = sta(4)
        xup_r  = 0.5 - sign ( onehalf , vu )
        xup_l  = 1. - xup_r
        vv     = xup_r * v_r + xup_l * v_l 
        vw     = xup_r * w_r + xup_l * w_l 
        vudtx  = vu * dtx
        fmass  = rho * vudtx
        predtx = pre * dtx
        fu     = fmass * vu + predtx
        fv     = fmass * vv
        fw     = fmass * vw
        pvu    = vu  * predtx / ( gam - 1. )
        fe     = pvu * gam 
     &         + 0.5 * fmass * ( vu**2 + vv**2 + vw**2 )
        f(1) = fmass 
        f(2) = fu 
        f(3) = fv 
        f(4) = fw 
        f(5) = fe 
        f(6) = pvu  
        f(7) = vu

#ifdef ADVECT_SPECIES
c
c.... Fluxes of advected variables
c
        do j = 8 , nhvarA-1
          u    = sta(2) 
          x0   = 0.5 + sign ( onehalf , u )
          x1   = 1. - x0
          dv0  = v(j,0) - v(j,-1)
          dv1  = v(j,1) - v(j,0)
          dv2  = v(j,2) - v(j,1)
          dv11 = v(j,1) - v(j,-1) 
          dv20 = v(j,2) - v(j,0) 
          dlq0 = ( 1. + sign ( 1.d0, dv1*dv0 ) )
     &           * min ( abs ( dv0 ) , abs ( dv1 ) )
          dlq1 = ( 1. + sign ( 1.d0, dv1*dv2 ) ) 
     &           * min ( abs ( dv1 ) , abs ( dv2 ) )
          fl   = v(j,0)
     &           + ( 0.5 - dtx2 * u ) 
     &           * sign ( min ( 0.5 * abs ( dv11 ), dlq0 ), dv11 )
          fr   = v(j,1)
     &            - ( 0.5 + dtx2 * u )
     &           * sign ( min ( 0.5 * abs ( dv20 ), dlq1 ), dv20 )
          f(j) = vudtx * ( fl * x0 + fr * x1 )
          
        end do
#endif

      return
      end
