c     ===================================================================
c
c      ART Version 3: ART_Auxiliary.f - auxiliary routines & functions
c
c     ===================================================================


c     --------------------------------------
      subroutine Timing ( ielement , jsign )
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
c              jsign    - specificator taking care of details of 
c                         the timing routine; for AIX (timef) should be:
c                         jsign = -1 before routine, 1 - after
c                         for dtime:
c                         jsign = 0 before routine,  1 - after
c
c     uses   : seconds
c
c     
      include 'a_control.h'
      integer ielement 
      REAL         last_second
      COMMON /ls/  last_second

      If(jsign.eq.1)Then
        CPU(ielement) = CPU(ielement) + ( seconds() - last_second )
      Else
        last_second = seconds()
      EndIf
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

      dummy = 0.0

c.... for IBM SP, intel compiler
c      dummy   = timef ()
c      seconds = dummy / 1000

c.... for exemplar
c      dummy = dtime(tarray) 
c      seconds = dtime(tarray)

c.... for hitachi
cccc      dummy = dtime(tarray)
c      seconds = second()

c.... for power challenge, pgi compiler
       dummy = etime(tarray)
       seconds = etime(tarray)

      return
      end 

c     -------------------------
      subroutine WriteTiming ()
c     -------------------------

      include 'a_tree.h'
      include 'a_control.h'

      CPU(10) = CPU(1)+CPU(2)+CPU(3)+CPU(5)+CPU(7)+CPU(8)

      if ( istep .eq. 1 ) then
        write(40,*) HEADER
        write(40,*) ' STEP ',' A ' , ' DEPTH ', ' Ncell ',
     &              ' DENS ' , ' POTENT ', ' RELAX ' , ' GRAVTOT ',
     &              ' MOVE ', ' LL_UPDATE ', 
     &              ' MODIFY ' , ' TOTAL '
        write(40,60) istep , aexpn , MaxLevelNow , ncell , 
     &               CPU(1) , CPU(2) , CPU(3) , CPU(4) , CPU(5) , 
     &               CPU(7) , CPU(8) , CPU(10)
      else
        write(40,60) istep , aexpn , MaxLevelNow , ncell , 
     &               CPU(1) , CPU(2) , CPU(3) , CPU(4) , CPU(5) , 
     &               CPU(7) , CPU(8) , CPU(10)
      endif

 60   format (i5 , 1x ,  f7.5 , 1x , i2 , 1x ,  i10 , 8(1x,f12.3) )

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

