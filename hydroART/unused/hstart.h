c
c     header file for the hydro initial conditions code 
c
c
ceverino08232005: commented out INDEX(NF67) for FFT due to incompatibilities with the intrinsic function index
c
      IMPLICIT REAL*8 (D)
      parameter (NROW=512,NGRID=128, NPAGE=NROW**2, NMAX=NGRID/2)
c      parameter (NRECL= NPAGE*6, NARR =MAX(2*NROW+1,NGRID+1))
c      parameter (NF67=MAX(NROW,NGRID/2))
      PARAMETER (NRECL= NPAGE*6, NARR =2*NROW+1)
      PARAMETER (NF67=NROW)
      PARAMETER (Nmaxpart = 6000000)
      PARAMETER (Lblock=2,Nblocks=NROW/Lblock)
      PARAMETER (NSPEC = NROW/2)   ! No waves shorter than Ny
      PARAMETER (marr  = NROW+1)   ! Arrays for FFT
      PARAMETER (mf67  = NROW/2)
c.... defines length of direct access record
c     Linux/absoft; SP2; Exemplar = 1
c     PC, O2K = 4
c
      parameter ( nbyteword = 1          )  

      COMMON / CONTROL/ AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                  ,Ocurv,Omb0,extras(100)

      real               ns      
      COMMON / TRUNCOM / Om,Omb,Omc,Omnu,Par(6),ns,qqscaleb,
     &                   QSCALE, SCALEL

      CHARACTER*45      HEADER
      COMMON / HEADDR/  HEADER

      COMMON /FOURAR/Zf(NARR),Yf(NARR)
      COMMON /F67COM/
     +                 IBC,      IP,       ISL,     L1,     N2,
     +                 N3,       N4,        N7,
     +                 SI(NF67), 
ceverino08232005   INDEX(NF67)

      common / SPEC01 / nsp(10,2)
      COMMON / ROW /	XPAR(NPAGE),YPAR(NPAGE),ZPAR(NPAGE),
     +			VX(NPAGE),VY(NPAGE),VZ(NPAGE)
      COMMON / BWEIG/ iWeight(Nmaxpart)

      DIMENSION         RECDAT(NRECL)  ,wspecies(10),lspecies(10)
      EQUIVALENCE    (RECDAT(1),XPAR(1))
     +                               ,(wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))

      COMMON / MASK/	iMask(Nblocks,Nblocks,Nblocks) ! Refinement Mask

      common / HYDRO /  rho0(ngrid, ngrid, ngrid),
     &                   vx0(ngrid, ngrid, ngrid),
     &                   vy0(ngrid, ngrid, ngrid),
     &                   vz0(ngrid, ngrid, ngrid),
     &                  eng0(ngrid, ngrid, ngrid),
     &                   ei0(ngrid, ngrid, ngrid)


      COMMON / BUFF /	XPt(Nmaxpart),YPt(Nmaxpart),ZPt(Nmaxpart),
     &			VXt(Nmaxpart),VYt(Nmaxpart),VZt(Nmaxpart), 
     &                   PW(NMaxpart)

