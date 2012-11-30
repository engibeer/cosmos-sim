c     ------------------------------------------------------
      subroutine FU ( dtx, dtx2, dtx22, index, g, v, c, f )
c     ------------------------------------------------------

c.... Computes total flux due to all physical processes

c     Input    :        dt    - time step
c                       dx    - cell size
c                       dtx   = dt / dx 
c                       dtx2  = 0.5 * dtx      
c                       dtx22 = 2. * dt / dx**2
c                       index - leaves left and right to interface
c                       g     - accelerations left and right to interface
c                       v     - variables for four leaves around interface
c                               1 - density
c                               2 - total energy
c                               3 - longitudinal momentum
c                             4,5 - transversal momenta
c                               6 - pressure
c                               7 - true adiabatic index

c     Output   :        f     - (flux * dt/dx)  across interface
c                        1    - density
c                        2-4  - momenta
c                        5    - total energy
c                        6    - internal energy (transport step flux)
c                        7    - dv/dx for div V in the source step
c                       v     - modified variables
c
c
#     include "a_def.h"
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'

      real*8 dtx, dtx2, dtx22

      integer index(2) 
      real*8 v(nhvarA,-1:2)         ! input for Interface Solvers 
      real*8 g(2)                 
      real*8 c(0:1)                 
      real*8 f(nhvarA-1)
      real gameng
#ifdef RadPressure
      real*8 RadP(-1:2)
c-------STG:
      common / RADPREVAR / RadP
C$OMP THREADPRIVATE (/RADPREVAR/)
c-----------
#endif

      double precision press, rhor, v1, v2, v3, ekine, eng
      real gammin, gammax
      double precision gam1
      parameter ( gammin = 1.01, gammax = 10. )
      parameter ( gam1 = 1.d-6 / ( gammax - 1. ) ) ! smallest pressure
#ifdef PRESSUREFLOOR
      real*8 f_eJ, f_eJ2 
      common / PFLOOR1 / f_eJ, f_eJ2
#endif

c      gam1 = enghdr / ( gammax - 1. )
c
      do k = -1 , 2
        rhor  = 1. / v(1,k)
        v1    = v(3,k) * rhor 
        v2    = v(4,k) * rhor       
        v3    = v(5,k) * rhor       
c        ekine = v(2,k) - 0.5 * ( v1*v1 + v2*v2 + v3*v3 ) * v(1,k)
c        ekine = v(2,k)
c        eng   = max ( enghdr , ekine )
#ifdef PRESSUREFLOOR

CEVERINO03242008        v(2,k) = max ( v(6,k) , f_eJ * v(1,k) * v(1,k) )

        v(2,k) = v(6,k) + f_eJ2 * v(1,k) ! Add turbulent pressure
#ifdef RadPressure
        v(2,k)=v(2,k) + RadP(k)   ! Add radiation pressure

c---------------STG:  
cccC$OMP CRITICAL
c        if (index(1) == 50009711) then
c           write(*,*) '--->',index(1), k, RadP(k)!, f_eJ*v(1,k)*v(1,k)
c        endif
c
c        if (f_eJ * v(1,k) * v(1,k) > v(2,k)) then
c           write(*,*) 'floor =', f_eJ * v(1,k) * v(1,k)
c        endif
ccC$OMP END CRITICAL
c---------------
#endif
        v(2,k) = max ( v(2,k), f_eJ * v(1,k) * v(1,k) )
  

CEVERINO04032008        v(2,k) = v(2,k) + f_eJ * v(1,k) * v(1,k) ! Add pressure floor
#endif
#ifndef PRESSUREFLOOR 
        v(2,k) = v(6,k) 
#endif
        v(3,k) = v1 
        v(4,k) = v2
        v(5,k) = v3
c        v(6,k) = press / eng + 1. ! effective gamma
        v(6,k) = v(7,k)  ! effective gamma 
        v(7,k) = gamma  ! true gamma 
      end do
c
      f(1) = 0.
      f(2) = 0.
      f(3) = 0.
      f(4) = 0.
      f(5) = 0.
      f(6) = 0.
      f(7) = 0.

#ifdef ADVECT_SPECIES
        do i = 8 , nhvarA-1
          f(i) = 0.
        enddo
#endif
c                

#ifdef HYDRO
      call FluxH ( dtx, dtx2, g, v, c, f )
#endif
c
      return
      end
c
c     -----------------------------------------------------------
      subroutine FD ( iDir, dtx2, index, Nbb1, Nbb2, g, v, f  )
c     -----------------------------------------------------------
c
c.... Adds Lapidus diffusion flux to the hydrodynamical flux computed by FU
c
c     Input    :        iDir  - sweep direction, 0 - turns off diffusion
c                       dtx2  = 0.5 * dt / dx, where
c                               dt    - time step
c                               dx    - cell size
c                       index  - leaves left and right to interface
c                       Nbb1,2 - neighbors of index(1) and index(2)
c                       g      - accelerations left and right to interface
c                       v      - variables for four leaves around interface
c                               1 - density
c                               2 - total energy
c                               3 - longitudinal momentum
c                             4,5 - transversal momenta
c                               6 - pressure
c                               7 - true adiabatic index
c
c     Output   :        f     - (flux * dt/dx)  across interface
c                               with added diffusion flux
c
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 dtx2
      dimension index(2)
      real*8 v(nhvarA,-1:2), g(2), f(nhvarA-1)
      real*8 diffk, diff, dvisc
      integer   iCh00, iCh10,
     &          iCh01, iCh11
c
      parameter ( diff1 = 0.1 )    ! controlls amount of diffusion
      parameter ( dviscmax = 0.1 ) 
      parameter ( drhomax = 0.2 )
c
      integer jpermute(7,6)
      data jpermute / 1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7 /
      integer j3, j4, j5 
      integer Nbb1(6), Nbb2(6)      
c
      j3 = jpermute (3,iDir)
      j4 = jpermute (4,iDir)
      j5 = jpermute (5,iDir)
c
      diffk = dtx2 * diff1 
c
c.... Compute velocity divergence
c
c     velocity difference in the direction of hydrodynamical sweep iDir
c
      gvisc = 2. * ( v(3,0) - v(3,1) )
c
c.... viscosity to smooth density gradients  
c
#ifdef DENSGRADSMOOTH
      xx = drhomax * max( v(1,1) , v(1,0)) / min( v(1,1) , v(1,0) )
      dvisc = max ( 0.0 , dviscmax * (xx - 1.0) / (xx + 1.0) )
#endif
   
c     
c     velocity difference in the directions orthogonal to iDir
c
      do j = 1 , ndim
        j1 = 2 * j
        j0 = j1 - 1
        jv = 2 + j
        if ( j0 .ne. iDir .and. j1 .ne. iDir ) then
          iCh00 = Nbb1 ( j0 )
          iCh01 = Nbb1 ( j1 )
          iCh10 = Nbb2 ( j0 )
          iCh11 = Nbb2 ( j1 )          
          v00 = hvar(jv,iCh00) / hvar(1,iCh00)
          v01 = hvar(jv,iCh01) / hvar(1,iCh01)
          v10 = hvar(jv,iCh10) / hvar(1,iCh10)
          v11 = hvar(jv,iCh11) / hvar(1,iCh11)
          gvisc = gvisc + v10 - v11 + v00 - v01 
        end if
      end do
c
c.... Add Lapidus diffusion flux to hydrodynamical fluxes
c
      i0   = index(1)
      i1   = index(2)
      diff = diffk * max ( 0.0 , gvisc )

#ifdef DENSGRADSMOOTH
      diff = max ( diff , dvisc )
#endif

      f(1) = f(1) + diff * ( hvar(1 ,i0) - hvar(1 ,i1) )
      f(2) = f(2) + diff * ( hvar(j3,i0) - hvar(j3,i1) )
      f(3) = f(3) + diff * ( hvar(j4,i0) - hvar(j4,i1) )
      f(4) = f(4) + diff * ( hvar(j5,i0) - hvar(j5,i1) )
      f(5) = f(5) + diff * ( hvar(2 ,i0) - hvar(2 ,i1) )
      f(6) = f(6) + diff * ( hvar(8 ,i0) - hvar(8 ,i1) )

#ifdef ADVECT_SPECIES
      do j = 8 , nhvarA-1
        f(j) = f(j) + diff * ( hvar(j+1,i0) - hvar(j+1,i1) )
      enddo
#endif
c
      return
      end
