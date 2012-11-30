c     --------------------------------------------------------------------
c
c     a_hfind.h - definitions for hfind.f (Andrey Kravtsov, december 1996)
c
c     --------------------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'

c...  constants
      parameter ( rhoaver = 1.0 )! average density of particles
      parameter ( akpc   = 3.0856e21 )        ! kiloparsec in cm
c      parameter ( h      = 0.05      )       ! resolution in grid units
c      parameter ( wpar = ((1.0*ng))**3/ncold ) ! particle weight
c      parameter ( rinit   = 100.0    )       ! initial halo radius in kpc
					      
      parameter ( nll     = 128      )        ! # of chain mesh cells in 1D
      parameter ( np1max  = 20000000 )        ! max # of the first specie
      parameter ( nhmax   = 300000   )        ! max # of halos
					      
c.... grid arrays 			      
      common / GRID1  / iCL2(nll,nll,nll)     ! particle linked list headers
      common / GRID2  / iCH(nll,nll,nll)      ! halo linked list headers
c      common / GRID3  / rho(ng,ng,ng)         ! density
					      
c.... particle arrays			      
      common / HPART01 / np1                  ! # of the first specie 
      common / HPART02 / pot(np1max)          ! potential at particle location
      common / HPART03 / pw2(np1max)          ! particle mass
      common / HPART04 / dnb2(np1max)         ! distance to the Nth clos.nb
      common / HPART05 / ind2(np1max)         ! index array for sorting     
      common / HPART06 / iSp(np1max)          ! 0 - free, 1 - belongs to a halo
      common / HPART07 / iLL2(np1max)         ! particle-mesh linked list 

c
c.... halo arrays 
c
      common / HALO01 /  xh(nhmax)            ! halo x-coordinates
      common / HALO02 /  yh(nhmax)            ! halo y-coordinates
      common / HALO03 /  zh(nhmax)            ! halo z-coordinates
      common / HALO04 / vxh(nhmax)            ! halo x-velocities
      common / HALO05 / vyh(nhmax)            ! halo y-velocities
      common / HALO06 / vzh(nhmax)            ! halo z-velocities
      common / HALO07 /  rh(nhmax)            ! halo radius up to Deltavir
      common / HALO08 / rhvir(nhmax)          ! virial radius of the halo
      common / HALO09 / rsh(nhmax)            ! scale radius r_s
      common / HALO10 / rst(nhmax)            ! stop radius (profile flattens)
      common / HALO11 / nhp(nhmax)            ! # of parts inside rh
      common / HALO12 / amh(nhmax)            ! halo mass 
      common / HALO13 / amhvir(nhmax)         ! halo virial mass 
      common / HALO14 / iHP(nhmax)            ! halo-particle LL headers
      common / HALO15 / iLH(nhmax)            ! particle-halo linked list
      common / HALO16 / iSh(nhmax)            ! 1 - radius reached 0 - not
      common / HALO17 / aMc(nhmax)            ! core mass
      common / HALO18 / DeltaMin , Deltavir   ! halo overdensity
      common / HALO19 / nhalo, ihc            ! current # of haloes 
      common / HALO20 / vhmax(nhmax)          ! velocity at rhmax
      common / HALO21 / rhmax(nhmax)          ! radius with max circular velocity
      common / HALO22 / ih1(nhmax)            ! halo ID
      common / HALO23 / ih2(nhmax)            ! halo ID

c
c...  halo profile arrays
c      
      parameter ( nbmax    = 200 )         ! size of the profile arrays
      real pn   (nil:nbmax,nhmax)          ! mass  profile
      real pnt  (nil:nbmax,nhmax)          ! tot mass within r 
      real na   (nil:nbmax,nhmax)          ! auxiliary arrays:
      real pvx  (nil:nbmax,nhmax)          ! x-velocity profile 
      real pvy  (nil:nbmax,nhmax)          ! y-velocity profile 
      real pvz  (nil:nbmax,nhmax)          ! z-velocity profile 
      real vcirc(nil:nbmax,nhmax)          ! circular velocity profile
      real vrms (nil:nbmax,nhmax)          ! 3D velocity dispersion profile

      common / PROFILES / pn, pnt, na, pvx, pvy, pvz, vcirc, vrms

c
c.... conversion factors
c
      common / FACTORS / box0 , box , pmmsun , fb,  
     &                   rMpc2g , rkpc2g , rg2Mpc , rg2kpc , vg2kms,
     &                   rg2pMpc, rg2pkpc

c
c...  globular cluster arrays
c
      parameter ( ngcmax = 100000 )
      real amgc(ngcmax)
      real*8 xgc(ngcmax), ygc(ngcmax), zgc(ngcmax), rgc(ngcmax)
      real vxgc(ngcmax), vygc(ngcmax), vzgc(ngcmax)
      real gcZII(ngcmax), gcZIa(ngcmax)
      integer igch(ngcmax)
      common / GC1 / xgc, ygc, zgc, rgc, vxgc, vygc, vzgc
      common / GC2 / gcZII, gcZIa, amgc

c
c...  Halo properties 
c
      integer nhgc(nhmax)
      real zhgc(nhmax), agcmax(nhmax), amtgc(nhmax)
      common / HH1 / igch, nhgc, zhgc, agcmax, amtgc

