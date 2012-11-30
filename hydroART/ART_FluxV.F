c     ------------------------------------------
      subroutine FluxV ( dtx22, g, v, c, f )
c     ------------------------------------------

c     Calculates viscous fluxes across the 0-1 interface
c     assuming incompressible flow.
c
c          -1         0         1         2
c     -!----*----!----*----!----*----!----*----!-


c     Input      : n                - number of arguments
c                  dt               - time step
c                  dx               - cell size, same for all cells
c                  dtx22            = 2. * dt / dx**2
c                  g(n,2)           - accelerations for cells 0,1
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
c     kinematic viscosity can be implemented to be different for cells 0,1
c
      include 'a_tree.h'
      parameter ( visc_k = 0.05 ) ! kinematic viscosity coefficient
c
      real*8 dtx22
      real*8 g(2), v(nhvar,-1:2), c(0:1)
      real*8 f(nhvar-1)
c
c      dtx22  = 2. * dt / dx**2
c
        rho0 =  v(1,0)
        rho1 =  v(1,1)
        d0 = visc_k * rho0
        d1 = visc_k * rho1
        visc = d0 * d1 / ( d0 + d1 )
        f(2) = f(2) + dtx * visc * ( v(3,0) - v(3,1) )
        f(3) = f(3) + dtx * visc * ( v(4,0) - v(4,1) )
        f(4) = f(4) + dtx * visc * ( v(5,0) - v(5,1) )
c
      return
      end
