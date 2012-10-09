c     ===================================================================
c
c      ART Version 3: ART_Init.f - initialization of variables and arrays
c
c     ===================================================================

c     -----------------------------
      subroutine Init_Parameters ()
c     -----------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))

c.... define cube corners

      xmin = 1.0
      ymin = 1.0 
      zmin = 1.0
      xmax = float(ng) + 1.0
      ymax = float(ng) + 1.0
      zmax = float(ng) + 1.0

      niter = 20
      niter1= 30
      niter2= 10
c.... initialize parameters for pyramidal interpolation
c      wa   = sqrt(11.)/(sqrt(11.) + 3.*sqrt(3.))   ! 0.38
c      wbcd = sqrt(3.)/(sqrt(11.) + 3.*sqrt(3.))

      wa   = 0.25
      wbcd = (1. - wa)/3

      MaxLevelNow = nil
      
c.... define small refinement box

      fsmall = 0.15 
      lbox   = 0.15 * ng
      if ( (0.15*ng-lbox) .lt. 0.5 ) lbox = lbox + 1
      if ( mod(lbox,2) .eq. 0 ) lbox = lbox + 1

      ic = ng / 2 + 1
      imin = ic - lbox/2
      imax = ic + lbox/2
      jmin = ic - lbox/2
      jmax = ic + lbox/2
      kmin = ic - lbox/2
      kmax = ic + lbox/2

c.... initialize particle weights if not MPI
c      wpar(1) = (wp1*float(ncell0)) / (wp1*ncold + wp2*(np-ncold))
c      wpar(2) = (wp2*float(ncell0)) / (wp1*ncold + wp2*(np-ncold))

      IF(mpi_yes .NE. 1234567 ) THEN      
       Do i =1,Nspecies
          wpar(i) = wspecies(i)
          lsp(i)     = lspecies(i)
       Enddo 
        Do i =1,Nspecies
           If(i.eq.1)Then
              nsp(i,1) = 1
           Else
              nsp(i,1) = lsp(i-1) +1
           Endif 
           nsp(i,2) = lsp(i)
        Enddo 
        Do i =1,Nspecies
           Do ip=nsp(i,1),nsp(i,2)
              Wpart(ip) =wpar(i)
           EndDo
        EndDo
        wpar_min = wpar(1)*1.0001  ! increased for test .LE.
      ENDIF

c.... initialize level-dependent density threshold

      nsplit(MinLevel) = 0
      nsplit(1) = 0
      nsplit(2) = 4

      trho(MinLevel) = wpar_min * nsplit(MinLevel) - 1.0  
      do Level = MinLevel + 1 , MaxLevel
        if ( Level .gt. 2 ) then 
c          nsplit(Level) = 3
c        else
          nsplit(Level) = 4
        endif
        if ( Level .ge. MaxLevel-1 )nsplit(Level) =10000 ! do not allow Maxlevel
        trho(Level) = wpar_min * nsplit(Level)
      enddo
      write(13,*) 'INIT: wpar_min =',wpar_min,' niter = ', niter
      write(13,100)  (nsplit(Level), Level =MinLevel,MaxLevel)
 100  FORMAT('INIT:  nsplit (level) = ',20I6)

c.... initialize cell sizes

      do Level = MinLevel , MaxLevel
        CellSize(Level) = 1.0 / 2**Level        
      enddo

      do Level = MinLevel , MaxLevel 
        Aexp    (Level) = aexpn
        iTimeBin(Level) = Level 
        dAexp   (Level) = astep / 2.0**iTimeBin(Level) 
      enddo

      damin = astep / 2.0**iTimeBin(MaxLevel)

      return
      end
      

c     -------------------------
      subroutine Init_Arrays ()
c     -------------------------
c
      include 'a_tree.h'

c*poption parallel
C$omp  parallel do default(shared)
C$omp& private( i , j  )
      do i = 1 , nctot
        iCL(i)    = nil 
        iOctCh(i) = nil 
        do j = 1 , nvar
          var(j,i) = zero
        enddo
      enddo
cC$omp end parallel do

      do Level = MinLevel+1 , MaxLevel
        iHOLL(Level) = nil
        iNOLL(Level) = nil
      enddo

c*poption parallel
C$omp  parallel do default(shared)
C$omp& private( j)
      do j = 1 , moct
        iOctLL1(j) = nil
        iOctLL2(j) = nil 
      enddo
cC$omp end parallel do
      return
      end

c     -----------------------
      subroutine Init_Tree ()  
c     -----------------------
c
c     initializes tree arrays & variables
c 
      include 'a_tree.h'

c*poption parallel
C$omp  parallel do default(shared)
C$omp& private(i)
      do i = 1 , moct
        iOctPr(i) = i+1          ! establish linked list of Octs #>0
        iOctLv(i) = iFreeLevel   ! oct type = free
      end do
cC$omp end parallel do
      iOctPr(moct) = 0
      iOctFree     = 1
      nOct         = 0

      return
      end

