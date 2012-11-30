c
c     ------------------------------------
      subroutine Riemann ( stl, str, sta )
c     ------------------------------------
c
c     Solves the Riemann problem at x/t=0
c     
c     by A.M.Khokhlov,   27 June 1994
c
c     Input     : nm           - vector length
c                 n            - number of arguments
c                 str(n,5)     - right states    
c                 stl(n,5)     - left states
c                                  (.,1) - density
c                                  (.,2) - velocity
c                                  (.,3) - pressure
c                                  (.,4) - Gamma = true adiabatic index
c                                  (.,5) - gamma = P/E + 1
c     
c     Output    : sta(n,4)     - state at x/t=ug
c                                  (.,1) - density
c                                  (.,2) - velocity
c                                  (.,3) - pressure
c                                  (.,4) - gamma = P/E + 1
c
c     
c     Tested for:                        1  <  rho    <    10**10
c                                 -5*10**9  <  u      <  5*10**9
c                                        1  <  P      <    10**30
c                                     1.01  <  Gamma  <  10
c                max ( 1.01, 0.5 * Gamma )  <  gamma  <  1.5 * Gamma
c
      include 'a_numbers.h'
      include 'a_setup.h'
c
      real*8 str(5), stl(5), sta(4)
c
c.... Internal variables/arrays
c
      real*8 p_l, p_r, p0, bgam_l, rgam_l, gam_l, xl2, xl3
      real*8 al, bl, cl, xx_l, cc_l, q_l, bgam_r, rgam_r, gam_r
      real*8 xr2, xr3, ar, br, cr, xx_r, cc_r, q_r, ul, ur
      real*8 p1, ul_0, ur_0, p_0, p_1, xxl, w2l, ul1, xxr, w2r, ur1
      real*8 p2, devi, dev, u, rho_s, u_s, p_s, bgam_s, gam_s
      real*8 a_s, b_s, c_s, w_s, rho, gam, xx1, xx2, xx3, xx4, xx5
      real*8 a1, a2, a3, a23, fs, f

c      dimension p_0 , p_1
c      dimension ul_0, ur_0
c      dimension al, bl, cl
c      dimension ar, br, cr
c
      real*8 small_R, eps
      parameter ( small_R   = 1.e-20 ,   ! small number for the solver
     &            eps     = 1.e-6    ,   ! convergence criterion
     &            maxit   = 50         ) ! maximum number of iterations
c
c.... Initial guess for the secant iteration scheme
c
        p_l      = stl(3)                   
        p_r      = str(3)                   
        p0       = 0.5 * ( p_l + p_r )
        bgam_l   = stl(4)
        rgam_l   = 1. / bgam_l
        gam_l    = stl(5)
        xl2      = 2. * gam_l * rgam_l - 1.
        cl       = xl2 * p_l
        xl3      = 0.5 * ( gam_l - 1. ) * xl2
        bl       = xl3 * p_l
        al       = gam_l - xl3
c        xx_l     = 0.5 * ( bgam_l - 1. ) * rgam_l
        q_l      = sqrt ( p_l * stl(1) * 
     &                    (( 1. + bgam_l ) * p0 / p_l + bgam_l - 1.) )
        bgam_r   = str(4)
        rgam_r   = 1. / bgam_r
        gam_r    = str(5)
        xr2      = 2. * gam_r * rgam_r - 1.
        cr       = xr2 * p_r
        xr3      = 0.5 * ( gam_r - 1. ) * xr2
        br       = xr3 * p_r
        ar       = gam_r - xr3
c        xx_r     = 0.5 * ( bgam_r - 1. ) * rgam_r
        q_r      = sqrt ( p_r * str(1) * 
     &                   (( 1. + bgam_r ) * p0 / p_r + bgam_r - 1.) )
        ul       = stl(2)
        ur       = str(2)
        p1       = (((ul - ur)*q_l*0.707106781d0 + p_l)*q_r + p_r*q_l)
     &           / ( q_l + q_r )
        ul_0     = ul     
        ur_0     = ur     
        p_0      = p0      
        p_1      = max ( small_R , p1 )
c
c.... Riemann solver - Secant iterations for pressure
c
      iter = 0
 1    dev = 0.
      if ( iter .le. 0 ) then
          xxl = ( al * p_1 + bl ) / ( p_1 + cl )
          w2l = 1./sqrt(max(small_R, xxl * stl(1) * (p_1 + stl(3))))
          ul1 = stl(2) + ( stl(3) - p_1 ) * w2l 
          xxr = ( ar * p_1  + br ) / ( p_1 + cr )
          w2r = 1./sqrt(max(small_R, xxr * str(1) * (p_1 + str(3))))
          ur1 = str(2) + ( p_1 - str(3) ) * w2r
          p2  = max ( small_R , 1.0000001 * p_1 - ( ur1 - ul1 ) 
     &                                  * abs ( p_1 - p_0 ) 
     &                                  / (   abs ( ur1 - ur_0 ) 
     &                                      + abs ( ul1 - ul_0 ) 
     &                                      + small_R ) )
          p_0   = p_1
          p_1   = p2
          ul_0  = ul1
          ur_0  = ur1
          devi     = abs ( p2 - p_0 ) / ( p2 + p_0 ) 
          sta(1) = devi 
          dev      = max ( dev , devi )
      else
        if ( sta(1) .gt. eps ) then
          xxl = ( al * p_1  + bl ) / ( p_1 + cl )
          w2l = 1./sqrt(max(small_R, xxl * stl(1) * (p_1 + stl(3))))
          ul1 = stl(2) + ( stl(3) - p_1 ) * w2l 
          xxr = ( ar * p_1  + br ) / ( p_1 + cr )
          w2r = 1./sqrt(max(small_R, xxr * str(1) * ( p_1 + str(3))))
          ur1 = str(2) + ( p_1 - str(3) ) * w2r
          p2  = max ( small_R , 1.0000001 * p_1 - ( ur1 - ul1 ) 
     &                                  * abs ( p_1 - p_0 ) 
     &                                  / (   abs ( ur1 - ur_0 ) 
     &                                      + abs ( ul1 - ul_0 ) 
     &                                      + small_R ) )
          p_0    = p_1
          p_1    = p2
          ul_0   = ul1
          ur_0   = ur1
          devi   = abs ( p2 - p_1 ) / ( p2 + p_1 ) 
          sta(1) = devi 
          dev    = max ( dev , devi )
        endif
      end if
      iter = iter + 1 
      if ( iter .le. maxit .and. dev .gt. eps ) go to 1
c
c..... Bad news !!!     
c
      if ( dev .gt. eps ) then
        write(*,'(1x,''Riemann_l solver iteration failure'')')
        stop
      end if
c
c.... State at x/t=0     
c
        u      = 0.5 * ( ul_0 + ur_0 )        
        ind_r  = int ( 0.9 - sign ( onehalf , u ) )
        rho_s  = ind_r * (str(1) - stl(1)) + stl(1)
        u_s    = ind_r * (str(2) - stl(2)) + stl(2)
        p_s    = ind_r * (str(3) - stl(3)) + stl(3)
        bgam_s = ind_r * (str(4) - stl(4)) + stl(4)
        gam_s  = ind_r * (str(5) - stl(5)) + stl(5)
        a_s    = ind_r * (ar - al) + al
        b_s    = ind_r * (br - bl) + bl
        c_s    = ind_r * (cr - cl) + cl
        w_s    = ( a_s * p_1  + b_s ) / ( p_1 + c_s ) 
        w_s    = max ( small_R, w_s * rho_s * ( p_1 + p_s ) )
        rho    = max ( small_R, 
     &                 rho_s / ( 1. - rho_s * ( p_1 - p_s ) / w_s ) )
        gam    = gam_s + 2. * ( gam_s - 1. ) * ( 1. - gam_s / bgam_s )
     &         / ( p_1 + p_s ) * ( p_1 - p_s )
        indd   = 2 * ind_r - 1
        if ( (p_1-p_s) .gt. 0 ) then 
          a2 = indd * u_s + sqrt ( w_s ) / rho_s
          a3 = 1.0e-10
        else
          xx4    = indd * u_s + sqrt ( bgam_s * p_s/rho_s )
          xx5    = indd * u +
     &             sqrt ( 0.5 * (stl(4)+str(4)) * p_1 / rho )
          a2 = xx5
          a3 = xx4 - xx5
        endif

        a23 = a2*a3

        if ( a23 .gt. zero ) then
c          fs = 0.
          sta(1) = rho
          sta(2) = u
          sta(3) = p_1
          sta(4) = gam
        else
          fs = min ( 1.d0 , -a2/a3 )
          sta(1) = rho + fs * (rho_s - rho)
          sta(2) = u   + fs * (u_s - u)
          sta(3) = p_1 + fs * (p_s - p_1)
          sta(4) = gam + fs * (gam_s - gam)
        endif

c        fs     = max ( zero , -a2/a3 )
c        fs     = min ( one , fs)
c        f      = 1. - fs  ! - sign is crucial
c        sta(1) = rho * f + rho_s * fs
c        sta(2) = u   * f + u_s   * fs
c        sta(3) = p_1 * f + p_s   * fs
c        sta(4) = gam * f + gam_s * fs
c
      return
      end




