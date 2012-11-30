c     ===================================================================
c
c      ART Version 3: ART_Auxiliary.f - auxiliary routines & functions
c
c     ===================================================================
c   
c     Routines for timing :
c       Timing
c       WallClockTime
c       seconds
c       Write_Timing
c     
c
c     Routines from the Numerical Recipe : 
c      ran0
c      DADT
c      ran0
c      ran3
c      DTDA
c      INTEGRATE
c      RUNGE5VAR
c      RUNGE
c      ZBRENT
c      RANDD
c      GASDEV
c      NEAREST_INT
c      LOGX
c      SPLINE
c      SPLINT
c      CHEBFT
c      CHEBEV
c      CHDER
c      CHINT
c      DFRIDR
c      POLCOF
c      POLINT
c      SAVGOL
c      LUDCMP
c      LUBKSB
c
c     --------------------------------------
      subroutine Timing ( ielement , isign )
c     --------------------------------------
c
c     purpose: to time a particular block of the code;
c              timings for a block are stored in CPU(ielement)
c 
c     input  : ielement - number of the block to time:
c              
c                     1 - Density
c                     2 - FFT
c                     3 - FMG_Relax
c                     4 - FFT + FMG_Relax
c                     5 - Move
c                     6 - Move Min Level
c                     7 - LL_Update
c                     8 - Modify 
c                     9 - Advance
c                    10 - Total time of a step
c
c              isign    - specificator taking care of details of 
c                         the timing routine; for AIX (timef) should be:
c                         isign = -1 before routine, 1 - after
c                         for dtime:
c                         isign = 0 before routine,  1 - after
c
c     uses   : seconds
c
c     
      integer ielement 
      include 'a_tree.h'
      include 'a_control.h'

      CPU(ielement) = CPU(ielement) + float(isign) * seconds()

      return
      end

c     -------------------------------------
      subroutine WallClockTime ( WallTime )
c     -------------------------------------
c
c     wrapper for HPCWALL routine
c
c     call this routine before and after a block of code to be timed 
c     WallTime after the 2nd call will be the wall clock time 
c     spent in the block
c
      implicit none
c
      real WallTime 
      real*8 WBegin, WEnd
      common / HPCTIME / WBegin
c     
      call HPCWALL ( WEnd )
c
      WallTime = WEnd - WBegin
      WBegin = WEnd
c
      return
      end

c     ------------------------
      real function seconds ()
c     ------------------------
c
c     purpose: returns elapsed time in seconds
c     uses   : xlf90 utility function timef()
c              (in this case - to be compiled with xlf90 or with -lxlf90)
c              or dtime (Exemplar) or etime (Power Challenge)
c      real*8 timef
      real*8 dummy
      real tarray(2)
      real*8 WCTime

c.... for apollon
c      dummy   = timef ()
c      seconds = dummy / 1000

c.... for exemplar
c      dummy = dtime(tarray) 
c      seconds = dtime(tarray)

c.... for power challenge
c      dummy = etime(tarray)
c      seconds = etime(tarray)
c
c.... general (uses gettimeofday in hpcwall.c)
c
      call HPCWALL ( WCTime )
      seconds = sngl(WCTime) 

      return
      end

c     -------------------------
      subroutine Write_Timing ()
c     -------------------------

      include 'a_tree.h'
      include 'a_control.h'

      call Timing ( ntiming-1 , 1 ) ! end timing step 

      CPU(ntiming) = CPU(1)+CPU(2)+CPU(3)+CPU(4)+CPU(5)+
     &               CPU(6)+CPU(7)+CPU(8)+CPU(9)+CPU(11)

      do i = 1 , ntiming
        CPU(i) = CPU(i) / 60.0  ! sec -> min
      enddo

      call Open_ASCII_File ( 40 , timingfile , append )

      ncells = noct * nchild + ncell0
      if ( istep .eq. 1 ) then
        write(40,58) HEADER
        write(40,59) 
      endif
      write(40,60) istep , t, aexpn , MaxLevelNow , ncells , 
     &             CPU(1), CPU(2) , CPU(3) , CPU(4) , CPU(5) , 
     &             CPU(6), CPU(7) , CPU(8) , CPU(9), CPU(11),
     &             CPU(ntiming), CPU(ntiming-1)

 58   format(A)
 59   format('istep','       t   ','      aexpn ',' MaxL ',
     &       'ncells ','  Dens ','    Pot ','   Relax ',
     &       '  CP+IPot ',
     &       ' TSH ', '   StepH ', '  Modify ', ' Move+LL ',
     &       ' SF ', '     IE  ',
     &       '   total   tstep  [all times in min]')
 60   format (i5,1x,e13.6,1x,f8.6,1x,i2,1x, i8,11(1x,f7.3),1x,f8.3 )

      close ( 40 )

      do i = 1 , ntiming
        CPU(i) = 0.
      enddo
c
      return
      end
c     -------------------
      function ran0(idum)     
c     -------------------
c
c     random number generator from NR
c
c     probably not to be used for serious purposes
c

      dimension v(97)
      data iff /0/

      if(idum.lt.0. or. iff.eq.0) then
         iff=1
         iseed=abs(idum)
         idum=1
         do j=1,97
            dum=ran3(iseed)
         enddo
         do j=1,97
            v(j)=ran3(iseed)
         enddo
         y=ran3(iseed)
      endif
      j=1+int(97.*y)
      if(j.gt.97 .or. j.lt.1) pause
      y=v(j)
      ran0=y
      v(j)=ran3(iseed)
      return
      end

      real function ran3(d)
      integer *4 k,i,j,l,m
      data i/646143019/,j/1048729165/,l/0/,m/1073741823/
      if (l.gt.0) goto 10
      fact=4./256.**4
      l=l+1
   10 continue
      k=iand(i+j,m)
      i=j
      j=k
    4 ran3=k*fact
      return
      end
c
c     ------------------------------------------
      double precision function DTDA(a, fp, np )
c     ------------------------------------------
c
c     Find dt/da given a, from the Friedmann Equation.  This is exact for
c     any isotropic-metric cosmology consistent with General Relativity.
c     Here, "t" is understood to be in units of the inverse Hubble constant
c     (i.e. "t" = H0*t).
c
c     Definitions for parameters are as in Peebles 1993, eqn (5.53).
c
c     a - expansion factor
c     fp(1) - Om0 - matter density wrt critical at z=0
c     fp(2) - Oml0 
c
      implicit none
      
      integer np 
      double precision a, fp(np), omegam, omegal

      double precision temp

      omegam = fp(1) 
      omegal = fp(2)

      temp = omegam + omegal*a**3 + (1.d0-omegam-omegal)*a

      if (temp.gt.0.d0) then
        DTDA = sqrt(a/temp)
      else
        DTDA = 0.d0
      end if

      return
      end

c     --------------------------------------------------------------
      double precision function INTEGRATE(FUNC,fp,np,a,b,dxinit,eps)
c     --------------------------------------------------------------
c
c     Quadrature using fifth order Runge-Kutta with adaptive step size.
c     Based on Press et al, Numerical Recipes in C, 2nd ed, pp 719-722.
c
c     Runge-Kutta driver with adaptive stepsize control.  Integrate starting
c     value y from a to b with accuracy eps, storing intermediate results in
c     global variables.  dxinit should be set as a guessed first stepsize.
c
c     Pass a vector of parameters of length np to FUNC in fp(np).
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC

      implicit none
      integer np 
      double precision a, b, eps, dxinit, FUNC, fp(np)
      external FUNC

      integer maxsteps
      parameter(maxsteps=100000000)

      double precision x, dx, dxnext, y, dydx, yscale
      integer  Nstep

      x     = a
      dx    = dxinit
      y     = 0.d0
      Nstep = 0

      do while ((x-b)*(b-a).lt.0.d0.and.Nstep.lt.maxsteps)
        Nstep = Nstep + 1
        dydx = FUNC(x,fp,np)
c
c       yscale is the scaling used to monitor accuracy.  This general-purpose
c       choice can be modified if need be.
c
        yscale = max(abs(y) + abs(dx*dydx), 1.d-20)
        if ((x+dx-b)*(x+dx-a).gt.0.d0)  ! If stepsize overshoots, decrease it.
     1    dx = b - x

        call RUNGE5VAR(y,dydx,x,dx,eps,yscale,dxnext,FUNC,fp,np)

        dx = dxnext
      end do

      if (Nstep.ge.maxsteps)
     1  write (*,*) 'WARNING: failed to converge in INTEGRATE.'

      INTEGRATE = y

      return
      end
c
c     -------------------------------------------------------------
      SUBROUTINE RUNGE5VAR(y,dydx,x,htry,eps,yscale,hnext,DERIVS,
     1                     fp,np)
c     -------------------------------------------------------------
c
c     Fifth-order Runge-Kutta step with monitoring of local truncation error
c     to ensure accuracy and adjust stepsize.  Input are the dependent
c     variable y and its derivative dydx at the starting value of the
c     independent variable x.  Also input are the stepsize to be attempted
c     htry, the required accuracy eps, and the value yscale, against which the
c     error is scaled.  On output, y and x are replaced by their new values.
c     hdid is the stepsize that was actually accomplished, and hnext is the
c     estimated next stepsize.  DERIVS is the user-supplied routine that
c     computes right-hand-side derivatives.  The argument fp is a vector 
c     of parameters (np parameters) to be passed to DERIVS 
c     (NOT integrated over).
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none
      integer np 
      double precision eps,hnext,htry,x,dydx,y,yscale,DERIVS,fp(np)
      external DERIVS

      double precision errmax,h,hold,htemp,xnew,yerr,ytemp

      double precision safety,pgrow,pshrink,errcon
      parameter (safety  =  0.9d0)
      parameter (pgrow   = -0.2d0)
      parameter (pshrink = -0.25d0)
      parameter (errcon  =  1.89d-4)

      h = htry                         ! Set stepsize to initial accuracy.
      errmax = 10.d0
      do while (errmax.gt.1.d0)
        call RUNGE(y,dydx,x,h,ytemp,yerr,DERIVS,fp,np)

        errmax = abs(yerr/yscale)/eps   ! Scale relative to required accuracy.
        if (errmax.gt.1.d0) then        ! Truncation error too large; reduce h
          htemp = safety*h*(errmax**pshrink)
          hold = h
          h = sign(max(abs(htemp),0.1d0*abs(h)),h)  ! No more than factor of 10
          xnew = x + h
          if (xnew.eq.x) then
            write (*,*) 'WARNING: ',
     1                  'Stepsize underflow in RUNGE5VAR().'
            h = hold
            errmax = 0.d0
          end if
        end if
      end do
c
c     Step succeeded.  Compute estimated size of next step.
c
      if (errmax.gt.errcon) then
        hnext = safety*h*(errmax**pgrow)
      else
        hnext = 5.d0 * h                ! No more than factor of 5 increase.
      end if
      x = x + h

      y = ytemp

      return
      end
c
c     ---------------------------------------------------
      SUBROUTINE RUNGE(y,dydx,x,h,yout,yerr,DERIVS,fp,np)
c     ---------------------------------------------------
c
c     Given values for a variable y and its derivative dydx known at x, use
c     the fifth-order Cash-Karp Runge-Kutta method to advance the solution
c     over an interval h and return the incremented variables as yout.  Also
c     return an estimate of the local truncation error in yout using the
c     embedded fourth order method.  The user supplies the routine
c     DERIVS(x,y,dydx), which returns derivatives dydx at x.
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none

      integer np
      double precision h,x,dydx,y,yerr,yout,DERIVS,fp(np)

      external DERIVS

      double precision ak3, ak4, ak5 ,ak6

      double precision a2,a3,a4,a5,a6
      double precision c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      parameter(a2  =    0.2d0)
      parameter(a3  =    0.3d0)
      parameter(a4  =    0.6d0)
      parameter(a5  =    1.d0)
      parameter(a6  =    0.875d0)
      parameter(c1  =   37.d0/378.d0)
      parameter(c3  =  250.d0/621.d0)
      parameter(c4  =  125.d0/594.d0)
      parameter(c6  =  512.d0/1771.d0)
      parameter(dc1 = c1 -  2825.d0/27648.d0)
      parameter(dc3 = c3 - 18575.d0/48384.d0)
      parameter(dc4 = c4 - 13525.d0/55296.d0)
      parameter(dc5 = -277.d0/14336.d0)
      parameter(dc6 = c6 -     0.25d0)

      ak3 = DERIVS(x+a3*h,fp,np)
      ak4 = DERIVS(x+a4*h,fp,np)
      ak5 = DERIVS(x+a5*h,fp,np)
      ak6 = DERIVS(x+a6*h,fp,np)
c
c     Estimate the fifth order value.
c
      yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6)
c
c     Estimate error as difference between fourth and fifth order
c
      yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

      return
      end
c
c     -------------------------------------
      FUNCTION zbrent(func,fp,np,x1,x2,tol)
c     -------------------------------------
c
      INTEGER ITMAX
      REAL zbrent
      integer np 
      real*8 fp(np)
      real*8 tol,x1,x2,EPS
      real*8 func
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b
      REAL c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a,fp,np)
      fb=func(b,fp,np)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then 
        write(*,*) 'root must be bracketed for zbrent'
        write(*,*) 'fa =',fa,'  fb =', fb
        pause
      endif
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,fp,np)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
c
c     -----------------
      FUNCTION RANDd(M)
c     -----------------
c     random number generator
c     initialize with arbitrary integer > 0
c     after that never change it
c     -----------------
      DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4
     +/453815927,2147483648.,2147483647,536870912,131072,256,16777216,
     +   4,16384,8388608,128/
      ML=M/K1*K1
      M1=(M-ML)*L1
      ML=M/K2*K2
      M2=(M-ML)*L2
      ML=M/K3*K3
      M3=(M-ML)*L3
      ML=M/K4*K4
      M4=(M-ML)*L4
      M5=KI-M
      IF(M1.GE.M5)M1=M1-KI-1
      ML=M+M1
      M5=KI-ML
      IF(M2.GE.M5)M2=M2-KI-1
      ML=ML+M2
      M5=KI-ML
      IF(M3.GE.M5)M3=M3-KI-1
      ML=ML+M3
      M5=KI-ML
      IF(M4.GE.M5)M4=M4-KI-1
      ML=ML+M4
      M5=KI-ML
      IF(LC.GE.M5)ML=ML-KI-1
      M=ML+LC
      RANDd=M/AM
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
c
c     -----------------
      FUNCTION gasdev(IDUM)
c     -----------------
c     return normally distributed deviate with 
c     zero mean and unit variance
c     initialize with arbitrary integer > 0
c     after that never change it
c     -----------------
      SAVE ISET,GSET
      DATA ISET/0/
      IF (IDUM.LT.0) ISET=0
      IF (ISET.EQ.0) THEN
1       V1=2.*RANDd(IDUM)-1.
        V2=2.*RANDd(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END
c
c     ----------------------------------
      integer function nearest_int ( x )
c     ----------------------------------
c     find nearest integer of the real x
c
      real x, frac 
      integer ix 

      ix = int(x)
      frac = x - ix
      
      if ( frac .lt. 0.5 ) then 
        nearest_int = ix 
      else
        nearest_int = ix + 1
      endif

      return
      end
c
c     --------------------------------
      real function logx ( xx , base )
c     --------------------------------
c
c     compute logarithm of base 'base'
c
      real xx, base
c
      logx = log(xx)/log(base)
c
      return
      end


c     -----------------------------------------------------------
c     SPLINE : Cubic spline interpolation 
c     -----------------------------------------------------------

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

c     -----------------------------------------------------------
c     SPLINT 
c     -----------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

c     -----------------------------------------------------------
c     CHEBFT : Chebyshev approximations 
c     -----------------------------------------------------------
      SUBROUTINE CHEBFT(A,B,C,N,FUNC)
      REAL*8 PI
      PARAMETER (NMAX=50, PI=3.141592653589793D0)
      REAL*8 SUM
      DIMENSION C(N),F(NMAX)
      BMA=0.5*(B-A)
      BPA=0.5*(B+A)
      DO 11 K=1,N
        Y=COS(PI*(K-0.5)/N)
        F(K)=FUNC(Y*BMA+BPA)
11    CONTINUE
      FAC=2./N
      DO 13 J=1,N
        SUM=0.D0
        DO 12 K=1,N
          SUM=SUM+F(K)*COS((PI*(J-1))*((K-0.5D0)/N))
12      CONTINUE
        C(J)=FAC*SUM
13    CONTINUE
      RETURN
      END
 

c     -----------------------------------------------------------
c     CHEBEV : Chebyshev evaluation
c     -----------------------------------------------------------
      FUNCTION CHEBEV(A,B,C,M,X)
      DIMENSION C(M)
      IF ((X-A)*(X-B).GT.0.) PAUSE 'X not in range.'
      D=0.
      DD=0.
      Y=(2.*X-A-B)/(B-A)
      Y2=2.*Y
      DO 11 J=M,2,-1
        SV=D
        D=Y2*D-DD+C(J)
        DD=SV
11    CONTINUE
      CHEBEV=Y*D-DD+0.5*C(1)
      RETURN
      END
 
c     -----------------------------------------------------------
c     CHDER
c     -----------------------------------------------------------
      SUBROUTINE CHDER(A,B,C,CDER,N)
      DIMENSION C(N),CDER(N)
      CDER(N)=0.
      CDER(N-1)=2*(N-1)*C(N)
      IF(N.GE.3)THEN
        DO 11 J=N-2,1,-1
          CDER(J)=CDER(J+2)+2*J*C(J+1)
11      CONTINUE
      ENDIF
      CON=2./(B-A)
      DO 12 J=1,N
        CDER(J)=CDER(J)*CON
12    CONTINUE
      RETURN
      END

c     -----------------------------------------------------------
c     CHINT
c     -----------------------------------------------------------
      SUBROUTINE CHINT(A,B,C,CINT,N)
      DIMENSION C(N),CINT(N)
      CON=0.25*(B-A)
      SUM=0.
      FAC=1.
      DO 11 J=2,N-1
        CINT(J)=CON*(C(J-1)-C(J+1))/(J-1)
        SUM=SUM+FAC*CINT(J)
        FAC=-FAC
11    CONTINUE
      CINT(N)=CON*C(N-1)/(N-1)
      SUM=SUM+FAC*CINT(N)
      CINT(1)=2.*SUM
      RETURN
      END

c     -----------------------------------------------------------
c     DFRIDR : derivative of func at x
c     -----------------------------------------------------------
      FUNCTION dfridr(func,x,h,err)
      INTEGER NTAB
      REAL dfridr,err,h,x,func,CON,CON2,BIG,SAFE
      PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.E30,NTAB=10,SAFE=2.)
      EXTERNAL func
C     UUSES func
      INTEGER i,j
      REAL errt,fac,hh,a(NTAB,NTAB)
      if(h.eq.0.) pause 'h must be nonzero in dfridr'
      hh=h
      a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)
      err=BIG
      do 12 i=2,NTAB
        hh=hh/CON
        a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
        fac=CON2
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dfridr=a(j,i)
          endif
11      continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12    continue
      return
      END
c     ----------------------------------------------------------------
c     POLCOF : return the coefficients of the interpolating polinomial
c     ----------------------------------------------------------------
      SUBROUTINE POLCOF(XA,YA,N,COF)
      PARAMETER (NMAX=15)
      DIMENSION XA(N),YA(N),COF(N),X(NMAX),Y(NMAX)
      DO 11 J=1,N
        X(J)=XA(J)
        Y(J)=YA(J)
11    CONTINUE
      DO 14 J=1,N
        CALL POLINT(X,Y,N+1-J,0.,COF(J),DY)
        XMIN=1.E38
        K=0
        DO 12 I=1,N+1-J
          IF (ABS(X(I)).LT.XMIN)THEN
            XMIN=ABS(X(I))
            K=I
          ENDIF
          IF(X(I).NE.0.)Y(I)=(Y(I)-COF(J))/X(I)
12      CONTINUE
        IF (K.LT.N+1-J) THEN
          DO 13 I=K+1,N+1-J
            Y(I-1)=Y(I)
            X(I-1)=X(I)
13        CONTINUE
        ENDIF
14    CONTINUE
      RETURN
      END
c     ----------------------------------------------------------------
c     POLINT : polinomial interpolation
c     ----------------------------------------------------------------
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10) 
      REAL X,Y,DY
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
c     ----------------------------------------------------------------
c     SAVGOL : Savitzky-Golay fitted polynomial
c     ----------------------------------------------------------------
      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      INTEGER ld,m,nl,np,nr,MMAX
      REAL c(np)
      PARAMETER (MMAX=6)
CU    USES lubksb,ludcmp
      INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)
      REAL d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
      if(np.lt.nl+nr+
     *1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m)
     *pause 'bad args in savgol'
      do 14 ipj=0,2*m
        sum=0.
        if(ipj.eq.0)sum=1.
        do 11 k=1,nr
          sum=sum+float(k)**ipj
11      continue
        do 12 k=1,nl
          sum=sum+float(-k)**ipj
12      continue
        mm=min(ipj,2*m-ipj)
        do 13 imj=-mm,mm,2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
13      continue
14    continue
      call ludcmp(a,m+1,MMAX+1,indx,d)
      do 15 j=1,m+1
        b(j)=0.
15    continue
      b(ld+1)=1.
      call lubksb(a,m+1,MMAX+1,indx,b)
      do 16 kk=1,np
        c(kk)=0.
16    continue
      do 18 k=-nl,nr
        sum=b(1)
        fac=1.
        do 17 mm=1,m
          fac=fac*k
          sum=sum+b(mm+1)*fac
17      continue
        kk=mod(np-k,np)+1
        c(kk)=sum
18    continue
      return
      END
c     ----------------------------------------------------------------
c     LUDCMP : LU Decomposition
c     ----------------------------------------------------------------
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
c     ----------------------------------------------------------------
c     LUBKSB : backsubstituion
c     ----------------------------------------------------------------
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
 
