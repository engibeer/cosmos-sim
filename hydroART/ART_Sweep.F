c
c     Parallel version of Sweep 
c
c     Andrey Kravtsov (2000-2002)
c


      subroutine iSwNb0 ( iC , idir , iR1 , iR2 , iL1  )

      include 'a_tree.h'

      integer iC, idir, iR1 , iR2 , iL1 , iL2
      parameter ( ng1 = ng - 1)
      dimension lnbr(0:8,6), kiot(0:8,6), ishift(6,3)
      data lnbr / 0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 5, 6, 7, 8, 1, 2, 3, 4,      ! idir = 5,6
     &            0, 5, 6, 7, 8, 1, 2, 3, 4/      ! idir = 5,6
      data kiot / 1, 0, 1, 0, 1, 0, 1, 0, 1,      ! idir = 1
     &            1, 1, 0, 1, 0, 1, 0, 1, 0,      ! idir = 2
     &            1, 0, 0, 1, 1, 0, 0, 1, 1,      ! idir = 3
     &            1, 1, 1, 0, 0, 1, 1, 0, 0,      ! idir = 4
     &            1, 0, 0, 0, 0, 1, 1, 1, 1,      ! idir = 5
     &            1, 1, 1, 1, 1, 0, 0, 0, 0/      ! idir = 6
      data ishift / -1, 1,  0, 0,  0, 0,
     &               0, 0, -1, 1,  0, 0,
     &               0, 0,  0, 0, -1, 1  /

        i = ( iC - 1 ) / ng2
        j = ( iC - 1 - i*ng2  ) / ng
        k =   iC - i*ng2 - j*ng
        i1 = i + ishift(idir,1)
        j1 = j + ishift(idir,2)
        k1 = k + ishift(idir,3)
        if ( i1 .gt. ng1 ) i1 = 0              ! periodic box check
        if ( j1 .gt. ng1 ) j1 = 0
        if ( k1 .gt. ng ) k1 = 1
        if ( i1 .lt.  0 ) i1 = ng1
        if ( j1 .lt.  0 ) j1 = ng1
        if ( k1 .lt.  1 ) k1 = ng
        i2 = i - ishift(idir,1)
        j2 = j - ishift(idir,2)
        k2 = k - ishift(idir,3)
        if ( i2 .gt. ng1 ) i2 = 0              ! periodic box check
        if ( j2 .gt. ng1 ) j2 = 0
        if ( k2 .gt. ng ) k2 = 1
        if ( i2 .lt.  0 ) i2 = ng1
        if ( j2 .lt.  0 ) j2 = ng1
        if ( k2 .lt.  1 ) k2 = ng

        iR1  = (i1*ng + j1)*ng + k1
        iL1  = (i2*ng + j2)*ng + k2

        i1 = i1 + ishift(idir,1)
        j1 = j1 + ishift(idir,2)
        k1 = k1 + ishift(idir,3)
        if ( i1 .gt. ng1 ) i1 = 0              ! periodic box check
        if ( j1 .gt. ng1 ) j1 = 0
        if ( k1 .gt. ng ) k1 = 1
        if ( i1 .lt.  0 ) i1 = ng1
        if ( j1 .lt.  0 ) j1 = ng1
        if ( k1 .lt.  1 ) k1 = ng
        iR2 = (i1*ng + j1)*ng + k1

      return
      end

      subroutine iSwNbL ( iC , idir , iR1 , iR2 , iL1 , iL2 , 
     &                                iLvR1 , iLvR2 , iLvL1  )
c
c
c      iL2 | iL1 | iC | iR1 | iR2 
c
      include 'a_tree.h'
      integer iC, idir, iR1, iR2, iL1, iL2, iLvR1, iLvR2, iLvL1
      dimension lnbr(0:8,6), kiot(0:8,6), nbnb(8,6)
      data lnbr / 0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 5, 6, 7, 8, 1, 2, 3, 4,      ! idir = 5,6
     &            0, 5, 6, 7, 8, 1, 2, 3, 4/      ! idir = 5,6
      data kiot / 1, 0, 1, 0, 1, 0, 1, 0, 1,      ! idir = 1
     &            1, 1, 0, 1, 0, 1, 0, 1, 0,      ! idir = 2
     &            1, 0, 0, 1, 1, 0, 0, 1, 1,      ! idir = 3
     &            1, 1, 1, 0, 0, 1, 1, 0, 0,      ! idir = 4
     &            1, 0, 0, 0, 0, 1, 1, 1, 1,      ! idir = 5
     &            1, 1, 1, 1, 1, 0, 0, 0, 0/      ! idir = 6
      data nbnb / 0, 1, 0, 3, 0, 5, 0, 7, 
     &            2, 0, 4, 0, 6, 0, 8, 0,
     &            0, 0, 1, 2, 0, 0, 5, 6, 
     &            3, 4, 0, 0, 7, 8, 0, 0,
     &            0, 0, 0, 0, 1, 2, 3, 4, 
     &            5, 6, 7, 8, 0, 0, 0, 0 /  
      integer iopp(6)
      data iopp / 2, 1, 4, 3, 6, 5 /

        iodir = iopp(idir)
        iC_   = iC - ncell0
        iO    = ishft ( iC_ + nshift , -ndim )
        iO_   = iO - 1                           ! cells oct
        icloc = iC_ - ishft( iO_ , ndim ) ! cells local # ( 1 - 8 )
        inl   = lnbr(icloc,idir)         ! neighbours local number
        iON   = iOctNb(idir,iO)             ! cells octs neighbouring cell
        iONC  = iOctCh(iON) - 1                ! neighbouring oct

      if ( kiot(icloc,idir) .eq. 0 ) Then          ! 0/1 - out/in key
        if ( iONC .eq. -1 ) then
          iR1 = iON
          iLvR1 = iLv(iR1)
          iR2 = iNb(iR1,idir)  ! work on this later!
          iLvR2 = iLv(iR2)
        else
          idum = ishft(iONC,ndim) + ncell0
          iR1 = inl  + idum
          iR2 = nbnb(inl,idir) + idum
          iLvR1 = iOctLv(iONC+1)
          iLvR2 = iLvR1
        endif
      else
        iR1   =  inl  + ishft ( iO_, ndim ) + ncell0
        iLvR1 = iOctLv(iO)
        if ( iONC .eq. -1 ) then
          iR2 = iON
          iLvR2 = iLv(iON)
        else            
        iR2 = icloc + ishft(iONC,ndim) + ncell0
c         iR2 = iNb(iR1,idir)
        iLvR2 = iOctLv(iONC+1)
        endif          
      endif

      inl  = lnbr(icloc,iodir)
      iON   = iOctNb(iodir,iO)             ! cells octs neighbouring cell
      iONC  = iOctCh(iON) - 1                ! neighbouring oct
      if ( kiot(icloc,iodir) .eq. 0 ) then          ! 0/1 - out/in key
      if ( iONC .eq. -1 ) then
        iL1 = iON
        iLvL1 = iLv(iON)
        iL2 = iNb(iL1,iodir)  ! work on this later!
      else
        idum = ishft(iONC,ndim) + ncell0
        iL1 = inl  + idum
        iLvL1 = iOctLv(iONC+1)
        iL2 = nbnb(inl,iodir) + idum
      endif
      else
        iL1   =  inl  + ishft ( iO_, ndim ) + ncell0
        iLvL1 = iOctLv(iO)
        if ( iONC .eq. -1 ) then
          iL2 = iON
        else
          iL2 = icloc + ishft(iONC,ndim) + ncell0
        endif
      endif

      return
      end

c     ------------------------------------
      subroutine Sweep ( iDir, dt, Level )
c     ------------------------------------
c     
c     Input     :    iDir            - direction of sweep ( 1 - 6 )
c                    dt              - time step
c                    Level           - level of sweep in the tree
c
#     include "a_def.h"
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 dt, dtx, dtx2, dtx22
      real*8 v(nhvarA,-1:2) ! input for Interface Solvers 
      real*8 g(2)                 
      real*8 c(2) 
      integer index(2)
      real*8 ff(nhvarA-1), f(nhvarA-1,nclmax)      
      common / Sweep_Com1 / f
      integer icell(2)           
      integer Nbb1(6), Nbb2(6)
c
      integer iP, iR1, iR2, iL1, iL2
      integer iLvR1, iLvR2, iLvL1
      integer iKyR1, j3, j4, j5
c
      integer jpermute(7,6)
      data jpermute / 1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7 /
c
      dimension aspect(2), asp(3)
      data asp / 0.5, 0.25, 0.125 /
      integer iopposite(6)
      data iopposite / 2, 1, 4, 3, 6, 5 /
      real*8 c112, c122
      parameter ( c112 = 1. / 1.25, c122 = 1. / 1.5 )
c
      integer iopp(6)
      data iopp / 2, 1, 4, 3, 6, 5 /
c
      integer iNbCh(8,3)
      data    iNbCh  / 2 , 1 , 4 , 3 , 6 , 5 , 8 , 7 , 
     &                 3 , 4 , 1 , 2 , 7 , 8 , 5 , 6 , 
     &                 5 , 6 , 7 , 8 , 1 , 2 , 3 , 4   /
      dimension iPyr(nchild,3)   ! interpolation pyramid vertices 
c
      data iPyr / 1, 2, 1, 2, 1, 2, 1, 2,   
     &            3, 3, 4, 4, 3, 3, 4, 4,
     &            5, 5, 5, 5, 6, 6, 6, 6  /

      integer iONB(nchild,3)
      
      data iONB / 2, 1, 4, 3, 6, 5, 8, 7,
     &            3, 4, 1, 2, 7, 8, 5, 6, 
     &            5, 6, 7, 8, 1, 2, 3, 4  / 
c
      real*8 rho_new,gravadd, rhov_new, dx, dxi, dxi2
      real*8 a2half, a2half2, ahalf, a2b, b2a
      external a2b, b2a
c
#ifdef TimeSweep
124	format('Level= ',I3,' iDir= ',I3,' NOct= ',I9,' HyTime= ',F8.3)
126	format('Level= ',I3,' iDir= ',I3,' NOct= ',I9,' GravTm= ',F8.3)
#endif

#ifdef PRESSUREFLOOR
      real*8 f_J, f_eJ
      real*8 f_J2, f_eJ2
      parameter ( f_J = 10.  ) ! pressure floor parameter
      parameter ( f_J2= 0.0 ) ! 3.047d4   )  ! turbulent pressure. 1.25 in code units. CEVERINO03242008
      common / PFLOOR1 / f_eJ, f_eJ2
#endif
#ifdef RadPressure
      real*8 RadP(-1:2)
c--------STG:
      common / RADPREVAR / RadP
C$OMP THREADPRIVATE (/RADPREVAR/)
c------------
#endif
#ifdef PRESSUREFLOOR
c
c.... artificial pressure floor a-la Machacek et al. 2
      if ( Level .ge. MinL_Jeans ) then
CEVERINO03272008        f_eJ = 0.47746 * f_J * aexp(Level) * CellSize(MaxLevel)**2 
         f_eJ = 0.47746 * f_J * aexp(Level) * CellSize(Level)**2 
         f_eJ2=       f_J2 * aexp(Level)**3 * CellSize(Level) 
CEVERINO04072008
      elseif ( Level .ge. MinL_Jeans-2.AND.Level.lt.MinL_Jeans) then 
         f_eJ = 0.0
         f_eJ2=       f_J2 * aexp(Level)**3 * CellSize(Level)          
      else
         f_eJ = 0.0
         f_eJ2= 0.0
      endif
#endif
#ifndef PRESSUREFLOOR
      f_eJ = 0.0
      f_eJ2= 0.0
#endif

c----STG:
c      write(*,*) 'f_eJ=',f_eJ, MinL_Jeans
c-----------------------

c
      dx    = CellSize(Level)   
      dxi   = 1.0 / dx
      dxi2  = 0.5 / dx
      dtx   = dt * dxi
      dtx2  = 0.5*dtx
      dtx22 = 2. * dt / dx**2
c
c.... a_1/2 =  a(tl+dtl/2) / (2*dx)
c
c
#ifdef COSMOLOGY
      ahalf   = b2a ( dble(tl(Level) + 0.5*dt) )
      ahalf   = ahalf * ahalf 
      a2half  = -0.5d0  * ahalf * dtx
      a2half2 = -0.25d0 * ahalf * dtx
#endif
#ifndef COSMOLOGY
      a2half  = -0.5d0 * dtx
      a2half2 = -0.25d0 * dtx
#endif
c              
      j3 = jpermute(3,iDir)
      j4 = jpermute(4,iDir)
      j5 = jpermute(5,iDir)
      iodir = iopp(idir)
      iDir2 = iDir / 2
c
      nif = 0
c
      IF ( Level .eq. MinLevel ) THEN
#ifdef TimeSweep
	call WallClockTime ( WallTime )
#endif
C.... OpenMP
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iP,j,iopp,iR1,iR2,iL1,iL2,iLvL1,iLvR1,iLvR2)
C$OMP+PRIVATE(v,icell,index,Nbb1,Nbb2,c,g,ff)
C$OMP+PRIVATE(ic0,ic1,i2)
C$OMP+FIRSTPRIVATE(nif)
        do iP = 1 , ncell0
          f(1,iP) = 0.
          f(2,iP) = 0.
          f(3,iP) = 0.
          f(4,iP) = 0.
          f(5,iP) = 0.          
          f(6,iP) = 0.          
          f(7,iP) = 0.
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            f(j,iP) = 0.
          enddo
#endif
          gacc(iP) = 0.
c          g(1) = 0.
c          g(2) = 0.

          if ( iOctCh(iP) .eq. nil ) then   ! only leaves
            call iSwNb0 ( iP , idir , iR1 , iR2 , iL1  )
c              if ( (hvar(1,iR1)-hvar(1,iP))*(hvar(1,iP)-hvar(1,iL1)) 
c     &             .lt. 0. ) then 
c                dd1 = min(abs(hvar(1,iR1)-hvar(1,iL1)),
c     &                    abs(2.0*(hvar(1,iP)-hvar(1,iL1))),
c     &                    abs(2.0*(hvar(1,iR1)-hvar(1,iP)))) * 
c     &                    sign(1.0, (hvar(1,iR1)-hvar(1,iL1))) / 
c     &                    hvar(1,iP)
c              else
c                dd1 = 0.
c              endif
c              if ( (hvar(1,iR2)-hvar(1,iR1))*(hvar(1,iR1)-hvar(1,iP)) 
c     &             .lt. 0. ) then 
c                dd2 = min(abs(hvar(1,iR2)-hvar(1,iP)),
c     &                    abs(2.0*(hvar(1,iR1)-hvar(1,iP))),
c     &                    abs(2.0*(hvar(1,iR2)-hvar(1,iR1)))) * 
c     &                    sign(1, (hvar(1,iR2)-hvar(1,iP))) / 
c     &                    hvar(1,iR1)
c              else
c                dd2 = 0.
c              endif
c              g(1) = a2half2 * (var(3,iR1) - var(3,iL1) + 
c     &               0.083333333d0 * dd1 * 
c     &                (var(3,iR1) - 2.0 *var(3,iP) + var(3,iL1))  )
c              g(2) = a2half2 * (var(3,iR2) - var(3,iP) + 
c     &               0.083333333d0 * dd2 * 
c     &                (var(3,iR2) - 2.0 *var(3,iR1) + var(3,iP))  )
c             gacc(iP) = 2.d0 * g(1)
#ifdef GRAVITY
            gacc(iP) = a2half * (var(3,iR1) - var(3,iL1)) 
#endif
c            g(2) = a2half2 * (var(3,iR2) - var(3,iP))
c            g(1) = gacc(iP) * 0.5d0
c
            if ( iOctCh(iR1) .eq. 0 ) then ! here we need to consider only =|= case
              v(1,-1) = hvar(1,iL1)
              v(1, 0) = hvar(1,iP )
              v(1, 1) = hvar(1,iR1)
              v(1, 2) = hvar(1,iR2)
              v(2,-1) = hvar(8,iL1)
              v(2, 0) = hvar(8,iP )
              v(2, 1) = hvar(8,iR1)
              v(2, 2) = hvar(8,iR2)
              v(3,-1) = hvar(j3,iL1)
              v(3, 0) = hvar(j3,iP )
              v(3, 1) = hvar(j3,iR1)
              v(3, 2) = hvar(j3,iR2)
              v(4,-1) = hvar(j4,iL1)
              v(4, 0) = hvar(j4,iP )
              v(4, 1) = hvar(j4,iR1)
              v(4, 2) = hvar(j4,iR2)
              v(5,-1) = hvar(j5,iL1)
              v(5, 0) = hvar(j5,iP )
              v(5, 1) = hvar(j5,iR1)
              v(5, 2) = hvar(j5,iR2)
              v(6,-1) = hvar(6,iL1)
              v(6, 0) = hvar(6,iP )
              v(6, 1) = hvar(6,iR1)
              v(6, 2) = hvar(6,iR2)
              v(7,-1) = hvar(7,iL1)
              v(7, 0) = hvar(7,iP )
              v(7, 1) = hvar(7,iR1)
              v(7, 2) = hvar(7,iR2)              
#ifdef ADVECT_SPECIES
              do j = 8 , nhvarA-1
                v(j,-1) = hvar(j+1,iL1)
                v(j, 0) = hvar(j+1,iP)
                v(j, 1) = hvar(j+1,iR1)
                v(j, 2) = hvar(j+1,iR2)
              enddo
#endif
#ifdef RadPressure
              RadP(-1)=RadPre(iL1)
              RadP(0)= RadPre(iP)
              RadP(1)= RadPre(iR1)
              RadP(2)= RadPre(iR2)
#endif
c
              index(1)      = iP
              index(2)      = iR1
              c(1)     = 1.0
              c(2)     = 1.0
              call FU ( dtx, dtx2, dtx22, index, g, v, c, ff )


c              if ( iP .eq. 4096 .or. iR1 .eq. 4096 ) then
c                write(*,*) '4096:', iP, iR1, dtx, j4, 
c     &          hvar(2,iP), hvar(2,iR1), ff(5), var(2,iP),var(2,iR1)
c	do jj = 1 , 7
c		 write(*,*) jj,v(jj,-1),v(jj,0),v(jj,1),v(jj,2) 
c                enddo  
c              endif


#ifdef LAPIDUS 
              call NbAll0 ( iP , Nbb1 )
              call NbAll0 ( iR1 , Nbb2 )
              call FD ( iDir, dtx2, index, Nbb1, Nbb2, g, v, ff )
#endif
              f(1,iP) = -ff(1)
              f(2,iP) = -ff(2)
              f(3,iP) = -ff(3)
              f(4,iP) = -ff(4)
              f(5,iP) = -ff(5)
              f(6,iP) = -ff(6)
              f(7,iP) = -ff(7) * dxi
#ifdef ADVECT_SPECIES
              do j = 8 , nhvarA-1
                f(j,iP) = -ff(j)
              enddo
#endif
              vnw(1 ,iR1) = vnw(1 ,iR1) + ff(1)
              vnw(2 ,iR1) = vnw(2 ,iR1) + ff(5)
              vnw(j3,iR1) = vnw(j3,iR1) + ff(2)
              vnw(j4,iR1) = vnw(j4,iR1) + ff(3)
              vnw(j5,iR1) = vnw(j5,iR1) + ff(4)
              vnw(6 ,iR1) = vnw(6 ,iR1) + ff(6)
              ref(iR1)    = ref(iR1)    + ff(7) * dxi
#ifdef ADVECT_SPECIES
              do j = 8 , nhvarA-1
                vnw(j-1,iR1) = vnw(j-1,iR1) + ff(j)
              enddo
#endif
              
            endif
          endif
c
         enddo

#ifdef TimeSweep
	call WallClockTime ( WallTime )
	write (27,124) Level, iDir, ncell0, WallTime
#endif
C.... OpenMP
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iP,j,rho_new,rhov_new,gravadd)
        do iP = 1 , ncell0
          if ( iOctCh(iP) .eq. nil ) then 
            rho_new    = vnw(1 ,iP) + f(1,iP)
            vnw(j4,iP) = vnw(j4,iP) + f(3,iP)
            vnw(j5,iP) = vnw(j5,iP) + f(4,iP)                     
            rhov_new   = vnw(j3,iP) + f(2,iP)            
            gravadd    = rho_new * gacc(iP)
            vnw(1 ,iP) = rho_new 
c            vnw(2 ,iP) = vnw(2 ,iP) + f(5,iP) 
c     &              +    0.5d0 * (rhov_new2**2 - 
c     &                            rhov_new1**2) / rho_new
            vnw(2 ,iP) = vnw(2 ,iP) + f(5,iP) + gacc(iP) * 
     &                   ( rhov_new + 0.5d0 * gravadd )
            vnw(j3,iP) = rhov_new   + gravadd
            vnw(6 ,iP) = vnw(6 ,iP) + f(6,iP)
            ref(iP)    = ref(iP)    + f(7,iP)
#ifdef ADVECT_SPECIES
            do j = 8 , nhvarA-1
              vnw(j-1,iP) = vnw(j-1,iP) + f(j,iP)
            enddo
#endif
          endif
        enddo
#ifdef TimeSweep
        call WallClockTime ( WallTime )
        write (27,126) Level, iDir, ncell0, WallTime
#endif
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
#ifdef TimeSweep
        call WallClockTime ( WallTime )

#endif
C.... OpenMP
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i1,j,ic,indc)
C$OMP+PRIVATE(iP,iT,iopp,iR1,iR2,iL1,iL2,iLvL1,iLvR1,iLvR2)
C$OMP+PRIVATE(v,icell,index,c,aspect,g,ff,iNb1,iNb2,iNb3)
C$OMP+PRIVATE(i2,i3,Nbb1,Nbb2,iChN,phi_b,phi_l,phi_r)
        do i1 = 1 , nLevel
          ic = iSelect(i1) - 1
          indc = (i1-1) * nchild 
          do i2 = 1 , 8
            iP = ic   + i2
               
            iT = indc + i2
            f(1,iT) = 0.
            f(2,iT) = 0.
            f(3,iT) = 0.
            f(4,iT) = 0.
            f(5,iT) = 0.
            f(6,iT) = 0.
            f(7,iT) = 0.
#ifdef ADVECT_SPECIES
            do j = 8 , nhvarA-1
              f(j,iT) = 0.0
            enddo
#endif

            gacc(iT) = 0.
c            g(1) = 0.
c            g(2) = 0.
            if ( iOctCh(iP) .eq. nil ) then   ! only leaves
              call iSwNbL ( iP , iDir , iR1 , iR2 , iL1 , iL2 , 
     &                      iLvR1 , iLvR2 , iLvL1 )
              iChN  = iNbCh ( i2 , iDir2 )
              if ( iLvR1 .eq. Level .and. iLvL1 .eq. Level ) then
                phi_l = var(3,iL1)
                phi_r = var(3,iR1)
              else
                if ( iLvL1 .lt. Level ) then
                  phi_l = Pyramide ( iL1 , iChN , 2 )
c                  phi_r = Pyramide ( iPr(iP) , iONb(i2,iDir2) , 2 )
                  phi_r = var(2,iR1)
                else
                  phi_l = var(2,iL1) 
c                  phi_l = Pyramide ( iPr(iP) , iONb(i2,iDir2) , 2 )
                  phi_r = Pyramide ( iR1 , iChN , 2 )
                endif
              endif
#ifdef GRAVITY
              gacc(iT) = a2half  * (phi_r - phi_l)
#endif

c              g(1) = gacc(iT) * 0.5d0 

              if ( iOctCh(iR1) .eq. 0 ) then 
                if ( iLvR1 .eq. Level ) then  ! = | = case
                  v(1,-1) = hvar(1,iL1)
                  v(2,-1) = hvar(8,iL1)
                  v(3,-1) = hvar(j3,iL1)
                  v(4,-1) = hvar(j4,iL1)
                  v(5,-1) = hvar(j5,iL1)
                  v(6,-1) = hvar(6,iL1)
                  v(7,-1) = hvar(7,iL1)
                  v(1, 2) = hvar(1, iR2)
                  v(2, 2) = hvar(8, iR2)
                  v(3, 2) = hvar(j3,iR2)
                  v(4, 2) = hvar(j4,iR2)
                  v(5, 2) = hvar(j5,iR2)
                  v(6, 2) = hvar(6, iR2)
                  v(7, 2) = hvar(7, iR2)
                  v(1, 0) = hvar(1,iP )
                  v(1, 1) = hvar(1,iR1)
                  v(2, 0) = hvar(8,iP ) 
                  v(2, 1) = hvar(8,iR1)
                  v(3, 0) = hvar(j3,iP )
                  v(3, 1) = hvar(j3,iR1)
                  v(4, 0) = hvar(j4,iP )
                  v(4, 1) = hvar(j4,iR1)
                  v(5, 0) = hvar(j5,iP )
                  v(5, 1) = hvar(j5,iR1)
                  v(6, 0) = hvar(6,iP )
                  v(6, 1) = hvar(6,iR1)
                  v(7, 0) = hvar(7,iP )
                  v(7, 1) = hvar(7,iR1)
#ifdef ADVECT_SPECIES
                  do j = 8 , nhvarA-1
                    v(j,-1) = hvar(j+1,iL1)
                    v(j, 0) = hvar(j+1,iP)
                    v(j, 1) = hvar(j+1,iR1)
                    v(j, 2) = hvar(j+1,iR2)
                  enddo
#endif
#ifdef RadPressure
              RadP(-1)=RadPre(iL1)
              RadP(0)=RadPre(iP)
              RadP(1)=RadPre(iR1)
              RadP(2)=RadPre(iR2)
#endif
c                  
                  index(1) = iP
                  index(2) = iR1

                  if ( iLvL1 .eq. Level ) then 
                    c(1) = 1.0
                  else
                    c(1)  = c112
                  endif
                  if ( iLvR2 .eq. Level ) then 
                    c(2) = 1.0 
                  else
                    c(2)  = c112
                  endif

c                  if (index(1).eq.482579 .or. index(2).eq.482579) then 
c                    idummy = 1
c                    write(*,*)  hvar(9,iP),(v(8,jjj),jjj=-1,2)
c                    pause
c                  endif

c----------STG:
cccC$OMP CRITICAL
c              if (iP == 50009711 ) then
c                write(*,*) iP,RadPre(iP),RadP(0),'1'
c              endif
cccC$OMP END CRITICAL
c--------------

                  call FU ( dtx, dtx2, dtx22, index, g, v, c, ff)
#ifdef LAPIDUS
                  call NbAllL ( iP , Nbb1 )
                  call NbAllL ( iR1 , Nbb2 )
                  call FD ( iDir, dtx2, index, Nbb1, Nbb2, g, v, ff )
#endif

                  f(1,iT) = -ff(1)
                  f(2,iT) = -ff(2)
                  f(3,iT) = -ff(3)
                  f(4,iT) = -ff(4)
                  f(5,iT) = -ff(5)
                  f(6,iT) = -ff(6)
                  f(7,iT) = -ff(7) * dxi
#ifdef ADVECT_SPECIES
                  do j = 8 , nhvarA-1
                    f(j,iT) = -ff(j)
                  enddo
#endif
                  vnw(1 ,iR1) = vnw(1 ,iR1) + ff(1)
                  vnw(2 ,iR1) = vnw(2 ,iR1) + ff(5)
                  vnw(j3,iR1) = vnw(j3,iR1) + ff(2)
                  vnw(j4,iR1) = vnw(j4,iR1) + ff(3)
                  vnw(j5,iR1) = vnw(j5,iR1) + ff(4)
                  vnw(6 ,iR1) = vnw(6 ,iR1) + ff(6)
                  ref(iR1)    = ref(iR1)    + ff(7) * dxi
#ifdef ADVECT_SPECIES
                  do j = 8 , nhvarA-1
                    vnw(j-1,iR1) = vnw(j-1,iR1) + ff(j)
                  enddo
#endif
                else  ! iLvR1 > Level => iLvL1 = Level
                  v(1,-1) = hvar(1,iL1)
                  v(1, 0) = hvar(1,iP )
                  v(2,-1) = hvar(8,iL1)
                  v(2, 0) = hvar(8,iP )
                  v(3,-1) = hvar(j3,iL1)
                  v(3, 0) = hvar(j3,iP )
                  v(4,-1) = hvar(j4,iL1)
                  v(4, 0) = hvar(j4,iP )
                  v(5,-1) = hvar(j5,iL1)
                  v(5, 0) = hvar(j5,iP )
                  v(6,-1) = hvar(6,iL1)
                  v(6, 0) = hvar(6,iP )
                  v(7,-1) = hvar(7,iL1)
                  v(7, 0) = hvar(7,iP )
                  v(1,1) = hvar(1,iR1)
                  v(2,1) = hvar(8,iR1)
                  v(3,1) = hvar(j3,iR1)
                  v(4,1) = hvar(j4,iR1)
                  v(5,1) = hvar(j5,iR1)
                  v(6,1) = hvar(6,iR1)
                  v(7,1) = hvar(7,iR1)
                  v(1,2) = hvar(1,iR2)
                  v(2,2) = hvar(8,iR2)
                  v(3,2) = hvar(j3,iR2)
                  v(4,2) = hvar(j4,iR2)
                  v(5,2) = hvar(j5,iR2)
                  v(6,2) = hvar(6,iR2)
                  v(7,2) = hvar(7,iR2)
#ifdef ADVECT_SPECIES
                  do j = 8 , nhvarA-1
                    v(j,-1) = hvar(j+1,iL1)
                    v(j, 0) = hvar(j+1,iP)
                    v(j, 1) = hvar(j+1,iR1)
                    v(j, 2) = hvar(j+1,iR2)
                  enddo
#endif
#ifdef RadPressure
              RadP(-1)=RadPre(iL1)
              RadP(0)=RadPre(iP)
              RadP(1)=RadPre(iR1)
              RadP(2)=RadPre(iR2)
#endif
                  index(1) = iP
                  index(2) = iR1
                  c(1)      = c112
                  c(2)      = c122
                  
c                  phi_r = Pyramide ( iR1 , i2 , 2 )
c                  g(2) = a2half2  * (phi_r - var(3,iP)) 

c                  if (index(1).eq.482579 .or. index(2).eq.482579) then 
c                    idummy = 1
c                    write(*,*)  hvar(9,iP),(v(8,jjj),jjj=-1,2)
c                    pause
c                  endif


c----------STG:
cccC$OMP CRITICAL
c              if (iP == 50009711 ) then
c                write(*,*) RadPre(iP),RadP(0),'2'
c              endif
cccC$OMP END CRITICAL
c--------------

                  call FU ( dtx, dtx2, dtx22, index, g, v, c, ff)

#ifdef LAPIDUS
                  call NbAllL ( iP , Nbb1 )
                  if ( iLvR1 .eq. MinLevel ) then 
                    call NbAll0 ( iR1 , Nbb2 )
                  else
                    call NbAllL ( iR1 , Nbb2 )
                  endif
                  call FD ( iDir, dtx2, index, Nbb1, Nbb2, g, v, ff )
#endif

                  f(1,iT) = -ff(1)
                  f(2,iT) = -ff(2)
                  f(3,iT) = -ff(3)
                  f(4,iT) = -ff(4)
                  f(5,iT) = -ff(5)
                  f(6,iT) = -ff(6)
                  f(7,iT) = -ff(7) * dxi
#ifdef ADVECT_SPECIES
                  do j = 8 , nhvarA-1
                    f(j,iT) = -ff(j)
                  enddo
#endif
                  vnw(1 ,iR1) = vnw(1 ,iR1) + ff(1)*asp(ndim)
                  vnw(2 ,iR1) = vnw(2 ,iR1) + ff(5)*asp(ndim)
                  vnw(j3,iR1) = vnw(j3,iR1) + ff(2)*asp(ndim)
                  vnw(j4,iR1) = vnw(j4,iR1) + ff(3)*asp(ndim)
                  vnw(j5,iR1) = vnw(j5,iR1) + ff(4)*asp(ndim)
                  vnw(6, iR1) = vnw(6, iR1) + ff(6)*asp(ndim)
                  ref(iR1)    = ref(iR1)    + ff(7)*asp(ndim) * dxi2
#ifdef ADVECT_SPECIES
                  do j = 8 , nhvarA-1
                    vnw(j-1,iR1) = vnw(j-1,iR1) + ff(j)*asp(ndim)
                  enddo
#endif
                endif
              endif
              if ( iLvL1 .eq. Level-1 ) then ! => iLvR1 = Level
                v(1,-1) = hvar(1,iL2)
                v(2,-1) = hvar(8,iL2)
                v(3,-1) = hvar(j3,iL2)
                v(4,-1) = hvar(j4,iL2)
                v(5,-1) = hvar(j5,iL2)
                v(6,-1) = hvar(6,iL2)
                v(7,-1) = hvar(7,iL2)
                v(1,0) = hvar(1,iL1)
                v(2,0) = hvar(8,iL1)
                v(3,0) = hvar(j3,iL1)
                v(4,0) = hvar(j4,iL1)
                v(5,0) = hvar(j5,iL1)
                v(6,0) = hvar(6,iL1)
                v(7,0) = hvar(7,iL1)

                v(1, 1) = hvar(1,iP )
                v(1, 2) = hvar(1,iR1)
                v(2, 1) = hvar(8,iP )
                v(2, 2) = hvar(8,iR1)
                v(3, 1) = hvar(j3,iP )
                v(3, 2) = hvar(j3,iR1)
                v(4, 1) = hvar(j4,iP )
                v(4, 2) = hvar(j4,iR1)
                v(5, 1) = hvar(j5,iP )
                v(5, 2) = hvar(j5,iR1)
                v(6, 1) = hvar(6,iP )
                v(6, 2) = hvar(6,iR1)
                v(7,-1) = hvar(7,iL2)
                v(7, 0) = hvar(7,iL1)
                v(7, 1) = hvar(7,iP )
                v(7, 2) = hvar(7,iR1)
#ifdef ADVECT_SPECIES
                do j = 8 , nhvarA-1
                  v(j,-1) = hvar(j+1,iL2)
                  v(j, 0) = hvar(j+1,iL1)
                  v(j, 1) = hvar(j+1,iP)
                  v(j, 2) = hvar(j+1,iR1)
                enddo
#endif
#ifdef RadPressure
              RadP(-1)=RadPre(iL1)
              RadP(0)=RadPre(iP)
              RadP(1)=RadPre(iR1)
              RadP(2)=RadPre(iR2)
#endif
                index(1) = iL1
                index(2) = iP
                c(1) = c122
                c(2) = c112

c                g(2) = g(1)
c                phi_l = Pyramide ( iL1 , i2 , 2 )
c                g(1) = a2half2  * (var(3,iP) - phi_l) 

c                  if (index(1).eq.482579 .or. index(2).eq.482579) then 
c                    idummy = 1
c                    write(*,*)  hvar(9,iP),(v(8,jjj),jjj=-1,2)
c                    pause
c                  endif

c----------STG:
cccC$OMP CRITICAL
c              if (iP == 50009711 ) then
c                write(*,*) RadPre(iP),RadP(0),'3'
c              endif
cccC$OMP END CRITICAL
c--------------

                call FU ( dtx, dtx2, dtx22, index, g, v, c, ff )

#ifdef LAPIDUS
                if ( iLvL1 .eq. MinLevel ) then 
                  call NbAll0 ( iL1 , Nbb1 )
                else
                  call NbAllL ( iL1 , Nbb1 )
                endif
                call NbAllL ( iP , Nbb2 )
                call FD ( iDir, dtx2, index, Nbb1, Nbb2, g, v, ff )
#endif

                f(1,iT) = f(1,iT) + ff(1)
                f(2,iT) = f(2,iT) + ff(2)
                f(3,iT) = f(3,iT) + ff(3)
                f(4,iT) = f(4,iT) + ff(4)
                f(5,iT) = f(5,iT) + ff(5)
                f(6,iT) = f(6,iT) + ff(6)
                f(7,iT) = f(7,iT) + ff(7) * dxi 
#ifdef ADVECT_SPECIES
                do j = 8 , nhvarA-1
                  f(j,iT) = f(j,iT) + ff(j)
                enddo
#endif
                vnw(1 ,iL1) = vnw(1 ,iL1) - ff(1)*asp(ndim)
                vnw(2 ,iL1) = vnw(2 ,iL1) - ff(5)*asp(ndim)
                vnw(j3,iL1) = vnw(j3,iL1) - ff(2)*asp(ndim)
                vnw(j4,iL1) = vnw(j4,iL1) - ff(3)*asp(ndim)
                vnw(j5,iL1) = vnw(j5,iL1) - ff(4)*asp(ndim)
                vnw(6, iL1) = vnw(6, iL1) - ff(6)*asp(ndim)
                ref(iL1)    = ref(iL1)    - ff(7)*asp(ndim) * dxi2
#ifdef ADVECT_SPECIES
                do j = 8 , nhvarA-1
                  vnw(j-1,iL1) = vnw(j-1,iL1) - ff(j)*asp(ndim)
                enddo
#endif
              endif
            endif
          enddo
       enddo
#ifdef TimeSweep
        call WallClockTime ( WallTime )
        write (27,124) Level, iDir, nLevel, WallTime
#endif
C.... OpenMP
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i1,j,ic,iT)
C$OMP+PRIVATE(rho_new,rhov_new,gravadd)
        do i1 = 1 , nLevel
          ic = iSelect(i1)
          iT = (i1-1) * nchild + 1
                       ! 1          
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif

          ic = ic + 1  ! 2
          iT = iT + 1
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif

          ic = ic + 1 ! 3
          iT = iT + 1
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif

          ic = ic + 1 ! 4
          iT = iT + 1
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif

          ic = ic + 1 ! 5
          iT = iT + 1
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif

          ic = ic + 1 ! 6 
          iT = iT + 1
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif

          ic = ic + 1 ! 7 
          iT = iT + 1
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif

          ic = ic + 1 ! 8 
          iT = iT + 1
          rho_new    = vnw(1 ,ic) + f(1,iT)
          vnw(j4,ic) = vnw(j4,ic) + f(3,iT)
          vnw(j5,ic) = vnw(j5,ic) + f(4,iT)          
          rhov_new   = vnw(j3,ic) + f(2,iT)
          gravadd    = rho_new * gacc(iT)
          vnw(1 ,ic) = rho_new 
          vnw(2 ,ic) = vnw(2 ,ic) + f(5,iT) + 
     &                 gacc(iT) * (rhov_new + 0.5d0*gravadd)
          vnw(j3,ic) = rhov_new + gravadd
          vnw(6 ,ic) = vnw(6 ,ic) + f(6,iT)
          ref(ic)    = ref(ic)    + f(7,iT)
#ifdef ADVECT_SPECIES
          do j = 8 , nhvarA-1
            vnw(j-1,ic) = vnw(j-1,ic) + f(j,iT)
          enddo
#endif 
       enddo
#ifdef TimeSweep
        call WallClockTime ( WallTime )
        write (27,126) Level, iDir, nLevel, WallTime
#endif
      ENDIF
c
      return
      end


