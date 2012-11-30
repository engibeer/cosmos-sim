c     =====================================================================
c                                                                         .
c          ART Version 3-1 : Gravitational Solver                         .
c                                                                         .
c         by Andrey Kravtsov and Anatoly Klypin (1995-99)                 .
c                                                                         .   
c     =====================================================================


c     --------------------------------------------------
      subroutine Solve_Poisson ( MinModify , MaxModify , ist )
c     --------------------------------------------------
c 
c     purpose: to organize routines which solve Poisson equations
c
      include 'a_tree.h'
      include 'a_control.h'
      
      integer MinModify , MaxModify 

      if ( (MinModify .lt. MinLevel)
     &                .or. 
     &     (MinModify .gt. MaxLevel) ) return
         
      if ( (MaxModify .lt. MinLevel)
     &                .or. 
     &     (MaxModify .gt. MaxLevel) ) return

      if ( MinModify .eq. MinLevel ) then

        call Timing ( 2 , -1 )
         call Potent ()
        call Timing ( 2 ,  1 )
        call Timing ( 3 , -1 )
         call Relax ( MinModify+1 , MaxModify , ist )
        call Timing ( 3 ,  1 )
      else
        call Timing ( 3 , -1 )
         call Relax ( MinModify , MaxModify , ist )
        call Timing ( 3 ,  1 )
      endif
c
      return
      end
c
c     ---------------------------------------------------
      subroutine Copy_Potential ( MinModify , MaxModify ) 
c     ---------------------------------------------------
c
c     copy current potential in var(2,*) to var(3,*)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level 
c
      DO Level = MinModify , MaxModify 
        IF ( Level .eq. MinLevel ) THEN
c
C$OMP PARALLEL DO DEFAULT(SHARED)
          do ic1 = 1 , ncell0
            var(3,ic1) = var(2,ic1)
          enddo
        ELSE
          call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell)
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            var(3,icell)   = var(2,icell)
            var(3,icell+1) = var(2,icell+1) 
            var(3,icell+2) = var(2,icell+2)
            var(3,icell+3) = var(2,icell+3)
            var(3,icell+4) = var(2,icell+4)
            var(3,icell+5) = var(2,icell+5)
            var(3,icell+6) = var(2,icell+6)
            var(3,icell+7) = var(2,icell+7)
          enddo
        ENDIF
      ENDDO
c
      return
      end
c
c     ----------------------------------------------------------
      subroutine Interpolate_Potential ( MinModify , MaxModify ) 
c     ----------------------------------------------------------
c
c     interpolate potential at time t_n in var(2,*) to t_n+1/2 in var(3,*)
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level 
c
      DO Level = MinModify , MaxModify 
        dtdt2 = 0.5 * dtl(Level) / dtlold(Level)
        IF ( Level .eq. MinLevel ) THEN
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(ic1)
          do ic1 = 1 , ncell0
            var(3,ic1) = var(2,ic1) + (var(2,ic1) - var(3,ic1)) * dtdt2
            rhophi(1,ic1) = 0.0
          enddo
        ELSE
          call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(ic1,icell,iT)
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            iT = (ic1-1) * nchild + 1

            var(3,icell)   = var(2,icell)   + 
     &                      (var(2,icell)   - var(3,icell)  ) * dtdt2
            var(3,icell+1) = var(2,icell+1) + 
     &                      (var(2,icell+1) - var(3,icell+1)) * dtdt2 
            var(3,icell+2) = var(2,icell+2) + 
     &                      (var(2,icell+2) - var(3,icell+2)) * dtdt2
            var(3,icell+3) = var(2,icell+3) + 
     &                      (var(2,icell+3) - var(3,icell+3)) * dtdt2
            var(3,icell+4) = var(2,icell+4) + 
     &                      (var(2,icell+4) - var(3,icell+4)) * dtdt2
            var(3,icell+5) = var(2,icell+5) + 
     &                      (var(2,icell+5) - var(3,icell+5)) * dtdt2
            var(3,icell+6) = var(2,icell+6) + 
     &                      (var(2,icell+6) - var(3,icell+6)) * dtdt2
            var(3,icell+7) = var(2,icell+7) + 
     &                      (var(2,icell+7) - var(3,icell+7)) * dtdt2
            rhophi(1,iT)   = 0.
            rhophi(1,iT+1) = 0.
            rhophi(1,iT+2) = 0.
            rhophi(1,iT+3) = 0.
            rhophi(1,iT+4) = 0.
            rhophi(1,iT+5) = 0.
            rhophi(1,iT+6) = 0.
            rhophi(1,iT+7) = 0.
          enddo
        ENDIF
      ENDDO
c
      return
      end

c     ------------------------------------------
      subroutine Relax ( MinModify , MaxModify , ist ) 
c     ------------------------------------------
c
c     purpose: main driver for relaxation iterations
c     input  : MinModify , MaxModify - levels to work on
c     note   : MinModify should be > MinLevel 
c              to ensure this is caller's responsibility
c
      include 'a_tree.h' 
      include 'a_control.h'
c
      iOctMax = moct
      do while ( iOctLv(iOctMax) .eq. iFreeLevel ) 
        iOctMax = iOctMax - 1
      enddo
      iOctMax = iOctMax + 1
c
      do Level = MinModify , MaxModify
         niter = 10



        if ( (Level.gt.MinLevel) .and. ist .eq. 0 ) then 
          call Prolongate ( Level-1 )
          niter = 60
        endif
        call Smooth ( Level )
        if ( Level .gt. MinLevel+1 ) then 
          call Restrict ( Level )
        endif
      enddo
      return
      end

c     -------------------------------
      subroutine Prolongate ( Level ) 
c     -------------------------------
c
c     makes prolongation sweep through Level
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer nLevel , Level 
      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iOC, iChild_1)
        do icell = 1 , ncell0
          iOC = iOctCh(icell)
          if ( iOC .gt. nil ) then ! if cell is split
            iChild_1 = ( iOC - 1 ) * nchild + 1 + ncell0
            var(2,iChild_1  ) = Pyramide ( icell , 1 , 2 ) 
            var(2,iChild_1+1) = Pyramide ( icell , 2 , 2 ) 
            var(2,iChild_1+2) = Pyramide ( icell , 3 , 2 )
            var(2,iChild_1+3) = Pyramide ( icell , 4 , 2 )
            var(2,iChild_1+4) = Pyramide ( icell , 5 , 2 )
            var(2,iChild_1+5) = Pyramide ( icell , 6 , 2 )
            var(2,iChild_1+6) = Pyramide ( icell , 7 , 2 )
            var(2,iChild_1+7) = Pyramide ( icell , 8 , 2 )
          endif
        enddo
      ELSE

        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell, iCh1,iCh2,iCh3,iCh4,iCh5,iCh6,iCh7,iCh8)
C$OMP+PRIVATE(idcell,iChild_1)
        do j = 1 , nLevel
          icell = iSelect(j)
          iCh1 = icell
          iCh2 = icell + 1
          iCh3 = icell + 2
          iCh4 = icell + 3
          iCh5 = icell + 4
          iCh6 = icell + 5
          iCh7 = icell + 6
          iCh8 = icell + 7

c.... 1st child
          idcell   = iCh1
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 )
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

c.... 2nd child
          idcell   = iCh2
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 ) 
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

c.... 3rd child
          idcell   = iCh3
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 ) 
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

c.... 4th child
          idcell   = iCh4
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 ) 
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

c.... 5th child
          idcell   = iCh5
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 ) 
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

c.... 6th child
          idcell   = iCh6
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 ) 
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

c.... 7th child
          idcell   = iCh7
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 ) 
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

c.... 8th child
          idcell   = iCh8
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var(2,iChild_1   ) = Pyramide ( idcell , 1 , 2 ) 
            var(2,iChild_1 +1) = Pyramide ( idcell , 2 , 2 ) 
            var(2,iChild_1 +2) = Pyramide ( idcell , 3 , 2 ) 
            var(2,iChild_1 +3) = Pyramide ( idcell , 4 , 2 ) 
            var(2,iChild_1 +4) = Pyramide ( idcell , 5 , 2 ) 
            var(2,iChild_1 +5) = Pyramide ( idcell , 6 , 2 ) 
            var(2,iChild_1 +6) = Pyramide ( idcell , 7 , 2 ) 
            var(2,iChild_1 +7) = Pyramide ( idcell , 8 , 2 ) 
          endif

       enddo

      ENDIF

      return
      end
c
c     ---------------------------------
      subroutine Prolongation ( iCell )  
c     ---------------------------------
c
c     purpose: prolongates density from parent (iCell) 
c     and its neighbors to the iCell's children
c     algorithm: "pyramidal" interpolation is used
c
      include 'a_tree.h'
      include 'a_control.h'
c
      dimension iPyr(nchild,3)   ! interpolation pyramid vertices 
c
      data iPyr / 1, 2, 1, 2, 1, 2, 1, 2,   
     &            3, 3, 4, 4, 3, 3, 4, 4,
     &            5, 5, 5, 5, 6, 6, 6, 6  /
c
      pdens = wa * var(2,iCell)      ! density contribution from parent

      do ic1 = 1 , nchild
        iChild = iCh(iCell,ic1)

        iNb1   = iNb(iCell,iPyr(ic1,1))
        iNb2   = iNb(iCell,iPyr(ic1,2))
        iNb3   = iNb(iCell,iPyr(ic1,3))

        var(2,iChild) = pdens + 
     &                  wbcd  * ( var(2,iNb1) +
     &                            var(2,iNb2) +         
     &                            var(2,iNb3) )

      enddo

      return
      end

c     -----------------------------
      subroutine Restrict ( Level ) 
c     -----------------------------
c
c     makes restriction sweep through Level
c
c     must be Level > MinLevel 
c
      include 'a_tree.h'
      real*8 sum, rfw
      parameter ( rfw = 1.d0 / nchild )
c
      nLevel = iNOLL(Level) ! cells are already selected in smooth
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell, ic2, idc, sum)
      do ic1 = 1 , nLevel
        icell = iSelect(ic1)
        sum = 0.d0
        do ic2 = 0 , 7
          idc = icell + ic2
          sum = sum + var(2,idc)          
        enddo
        var(2,iPr(icell)) = sum * rfw 
      enddo
      
      return
      end
c
c     ---------------------------
      subroutine Smooth ( Level )   
c     ---------------------------
c
c     purpose: smoothes Level > 0
c     input  : Level  - level to be smoothed
c
      include 'a_tree.h'
      include 'a_control.h'

#     include "a_def.h"

c      parameter ( nrb = (nctot-ncell0) )
      parameter ( nrb = nclmax+1 ) 

      real*8 wsor, wsor6, six, rhoJ
      real*8 Size, trfi1, trfi2
      real*8 phi, phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8
      real*8 phicell, enorm, differ, sum
      integer iLevNb(nrb,3)
      common / GravCom1 / iLevNb

      dimension iENb(nchild,3)   ! external neighbors
      dimension iINb(nchild,3)   ! internal neighbors
      dimension iNbC(nchild,3)   ! neighbors' child (used in interpolation)

      data iENb / 1, 2, 1, 2, 1, 2, 1, 2,   
     &            3, 3, 4, 4, 3, 3, 4, 4,
     &            5, 5, 5, 5, 6, 6, 6, 6  /

      data iINb / 2, 1, 2, 1, 2, 1, 2, 1,   
     &            4, 4, 3, 3, 4, 4, 3, 3,
     &            6, 6, 6, 6, 5, 5, 5, 5  /
c
      data iNbC / 2, 1, 4, 3, 6, 5, 8, 7,   
     &            3, 4, 1, 2, 7, 8, 5, 6,
     &            5, 6, 7, 8, 1, 2, 3, 4  /

      Size  = dble(CellSize(Level))     ! cell size on this level      
      trfi1 = aexp(Level) / Size
c      trfi2 = aexp(Level) * Size * Size * wsor

      IF ( Level .eq. MinLevel+1 ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
        do j = 1 , ncell0
          ind(j) = 0 
        enddo
      ELSE
        call Select_Cells ( Level-1 , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell)
        do j = 1 , nLevel
          icell = iSelect(j)
          ind(icell)   = 0
          ind(icell+1) = 0 
          ind(icell+2) = 0 
          ind(icell+3) = 0 
          ind(icell+4) = 0 
          ind(icell+5) = 0 
          ind(icell+6) = 0 
          ind(icell+7) = 0 
        enddo
      ENDIF

      call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell)
      do j = 1 , nLevel
        icell = iSelect(j)
        ind(icell)   = 0
        ind(icell+1) = 0 
        ind(icell+2) = 0 
        ind(icell+3) = 0 
        ind(icell+4) = 0 
        ind(icell+5) = 0 
        ind(icell+6) = 0 
        ind(icell+7) = 0 
      enddo
c
c.... red
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,jcell,idcell)
      do j = 1 , nLevel 
        icell = iSelect(j)
        jcell = j*4 - 3
        rhophi(1,jcell) = var(1,icell)
        rhophi(2,jcell) = var(2,icell)
        ind(icell) = jcell

        idcell = icell + 3
        jcell  = jcell + 1
        rhophi(1,jcell) = var(1,idcell)
        rhophi(2,jcell) = var(2,idcell)
        ind(idcell) = jcell

        idcell = icell + 5
        jcell  = jcell + 1
        rhophi(1,jcell) = var(1,idcell)
        rhophi(2,jcell) = var(2,idcell)
        ind(idcell) = jcell

        idcell = icell + 6
        jcell  = jcell + 1
        rhophi(1,jcell) = var(1,idcell)
        rhophi(2,jcell) = var(2,idcell)
        ind(idcell) = jcell
      enddo
c
c.... black
c
      iblack = nLevel*4 - 3

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,idcell,jcell)
      do j = 1 , nLevel 
        icell = iSelect(j) 
        idcell = icell + 1
        jcell = iblack + j*4 
        rhophi(1,jcell) = var(1,idcell)
        rhophi(2,jcell) = var(2,idcell)
        ind(idcell) = jcell

        idcell = icell + 2
        jcell  = jcell + 1
        rhophi(1,jcell) = var(1,idcell)
        rhophi(2,jcell) = var(2,idcell)
        ind(idcell) = jcell

        idcell = icell + 4
        jcell  = jcell + 1
        rhophi(1,jcell) = var(1,idcell)
        rhophi(2,jcell) = var(2,idcell)
        ind(idcell) = jcell

        idcell = icell + 7
        jcell  = jcell + 1
        rhophi(1,jcell) = var(1,idcell)
        rhophi(2,jcell) = var(2,idcell)
        ind(idcell) = jcell
      enddo
c
c.... Assign Neighbors 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,jcell,idcell)
      do j = 1 , nLevel
        icell = iSelect(j)
        jcell = j*4 - 3
c
c.... 1st sibling
c        
        iLevNb(jcell,1) = iNb(icell,1)
        iLevNb(jcell,2) = iNb(icell,3)
        iLevNb(jcell,3) = iNb(icell,5)              
c   
c.... 4th sibling
c
        idcell = icell + 3
        jcell  = jcell + 1
        iLevNb(jcell,1) = iNb(idcell,2)
        iLevNb(jcell,2) = iNb(idcell,4)
        iLevNb(jcell,3) = iNb(idcell,5)
c
c.... 6th sibling
c
        idcell = icell + 5
        jcell  = jcell + 1
        iLevNb(jcell,1) = iNb(idcell,2)
        iLevNb(jcell,2) = iNb(idcell,3)
        iLevNb(jcell,3) = iNb(idcell,6)
c
c.... 7th sibling
c
        idcell = icell + 6
        jcell  = jcell + 1
        iLevNb(jcell,1) = iNb(idcell,1)
        iLevNb(jcell,2) = iNb(idcell,4)
        iLevNb(jcell,3) = iNb(idcell,6)
c
      enddo ! end j

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,jcell,idcell)
      do j = 1 , nLevel
        icell = iSelect(j) + 1
        jcell = iblack + j*4
c
c.... 2nd sibling
c        
        iLevNb(jcell,1) = iNb(icell,2)
        iLevNb(jcell,2) = iNb(icell,3)
        iLevNb(jcell,3) = iNb(icell,5)              
c
c.... 3rd sibling
c
        idcell = icell + 1
        jcell  = jcell + 1
        iLevNb(jcell,1) = iNb(idcell,1)
        iLevNb(jcell,2) = iNb(idcell,4)
        iLevNb(jcell,3) = iNb(idcell,5)
c
c.... 5th sibling
c
        idcell = icell + 3
        jcell  = jcell + 1
        iLevNb(jcell,1) = iNb(idcell,1)
        iLevNb(jcell,2) = iNb(idcell,3)
        iLevNb(jcell,3) = iNb(idcell,6)
c
c.... 8th sibling
c
        idcell = icell + 6
        jcell = jcell + 1
        iLevNb(jcell,1) = iNb(idcell,2)
        iLevNb(jcell,2) = iNb(idcell,4)
        iLevNb(jcell,3) = iNb(idcell,6)
c
      enddo ! end j
c
c.... pre-compute potential on the border
c     the following loop is to be executed SERIALLY
c
c.... red
c
      ibc = nLevel*8
      do j = 1 , nLevel
        jcell = j*4 - 3
c
c....   1st sibling
c
        if ( ibc+12 .gt. nclmax ) then 
          call Open_ASCII_File ( iOErr , errorfile , sequent ) 
          write(iOErr,*) '1  * error: in Smooth:'
          write(iOErr,*) 'L =',Level,' ibc=',ibc
          write(iOErr,*) ' ibc+12 > nclmax =',nclmax
          write(iOErr,*) ' increase nclmax in a_setup.h and rerun.'
          write(iOErr,*) ' stopping...'
          close(iOErr)
          stop        
        endif

        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,2,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,3,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,5,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif
c
c....   4th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,3,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,2,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,8,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif
c
c....   6th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,5,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,8,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,2,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif
c
c....   7th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,8,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,5,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,3,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif
c
      enddo
c
c.... black
c
      do j = 1 , nLevel
        jcell = iblack + j*4 

        if ( ibc+12 .gt. nclmax ) then 
          call Open_ASCII_File ( iOErr , errorfile , sequent ) 
          write(iOErr,*) '1  * error: in Smooth:'
          write(iOErr,*) 'L =',Level,' ibc=',ibc
          write(iOErr,*) ' ibc+12 > nclmax =',nclmax
          write(iOErr,*) ' increase nclmax in a_setup.h and rerun.'
          write(iOErr,*) ' stopping...'
          close(iOErr)
          stop        
        endif

        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,1,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,4,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,6,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif

        jcell = jcell + 1
        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,4,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,1,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,7,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif
c
        jcell = jcell + 1
        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,6,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,7,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,1,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif

        jcell = jcell + 1
        inb1 = iLevNb(jcell,1)
        inb2 = iLevNb(jcell,2)
        inb3 = iLevNb(jcell,3)
        if (ind(inb1) .eq. 0 ) then
          ibc = ibc + 1          
          rhophi(2,ibc) = Pyramide(inb1,7,2)
          iLevNb(jcell,1) = ibc
        else
          iLevNb(jcell,1) = ind(inb1)
        endif
        if (ind(inb2) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb2,6,2)
          iLevNb(jcell,2) = ibc
        else
          iLevNb(jcell,2) = ind(inb2)
        endif
        if (ind(inb3) .eq. 0) then
          ibc = ibc + 1
          rhophi(2,ibc) = Pyramide(inb3,4,2)
          iLevNb(jcell,3) = ibc
        else
          iLevNb(jcell,3) = ind(inb3)
        endif
c
      enddo
 
      phimin = 1.d-2 * phimin
CEVERINO11022005      niter = 60 ! it was 15  
      wsor  = 1.d0
      rhoJ  = 0.95d0
      Size2 = Size * Size
      maxit = 60
      minit = 3
      iter  = nil 
      six = 6.d0 

      DO iter = 1 , niter

        wsor6 = wsor / 6.d0  
        trfi2 = wsor / aexp(Level) / Size 
c
c.....  split sweep in 2 (odd/even) parts
c          
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j4,jcell1,jcell2,jdcell,inb1,inb2,inb3)
C$OMP+PRIVATE(phi,phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8)
        do j = 1 , nLevel
          j4 = j * 4   
          jcell1 = j4 - 3
          jcell2 = j4 + iblack
c
c....     1st sibling
c              
          inb1 = iLevNb(jcell1,1)
          inb2 = iLevNb(jcell1,2)
          inb3 = iLevNb(jcell1,3)              
          phi  = rhophi(2,inb1)
          phi  = phi + rhophi(2,inb2)
          phi  = phi + rhophi(2,inb3)

          phi2  = rhophi(2,jcell2)
          phi3  = rhophi(2,jcell2+1)
          phi5  = rhophi(2,jcell2+2)
          phi1  = rhophi(2,jcell1)

          phi1 =  phi1 + 
     &           wsor6 * ( 
     &                           (phi + phi2 + phi3 + phi5) - 
     &                           six*phi1
     &                         ) - trfi2 * rhophi(1,jcell1)
c   
c....     4th sibling
c
          jdcell = jcell1 + 1
          inb1   = iLevNb(jdcell,1)
          inb2   = iLevNb(jdcell,2)
          inb3   = iLevNb(jdcell,3)
          phi    = rhophi(2,inb1)
          phi    = phi + rhophi(2,inb2)
          phi    = phi + rhophi(2,inb3)

          phi8  = rhophi(2,jcell2+3)
          phi4  = rhophi(2,jdcell)

          phi4 = phi4 + 
     &            wsor6 * ( 
     &                           (phi + phi3 + phi2 + phi8) - 
     &                            six*phi4
     &                          ) - trfi2 * rhophi(1,jdcell)
c
c....     6th sibling
c
          jdcell = jcell1 + 2
          inb1   = iLevNb(jdcell,1)
          inb2   = iLevNb(jdcell,2)
          inb3   = iLevNb(jdcell,3)
          phi    = rhophi(2,inb1)
          phi    = phi + rhophi(2,inb2)
          phi    = phi + rhophi(2,inb3)

          phi6 = rhophi(2,jdcell)
          phi6 = phi6 + 
     &           wsor6 * ( 
     &                           (phi + phi5 + phi8 + phi2) - 
     &                            six* phi6
     &                          ) - trfi2 * rhophi(1,jdcell)
c
c....     7th sibling
c
          jdcell = jcell1 + 3
          inb1   = iLevNb(jdcell,1)
          inb2   = iLevNb(jdcell,2)
          inb3   = iLevNb(jdcell,3)
          phi    = rhophi(2,inb1)
          phi    = phi + rhophi(2,inb2)
          phi    = phi + rhophi(2,inb3)

          phi7 = rhophi(2,jdcell)
          phi7 =  phi7 + 
     &           wsor6 * ( 
     &                           (phi + phi8 + phi5 + phi3) - 
     &                            six* phi7
     &                    ) - trfi2 * rhophi(1,jdcell)
c
          rhophi(2,jcell1)   = phi1
          rhophi(2,jcell1+1) = phi4
          rhophi(2,jcell1+2) = phi6
          rhophi(2,jcell1+3) = phi7          
c
        enddo ! end j
        


        if ( iter .eq. 1 ) then 
          wsor = 1.d0 / (1.d0 - 5.d-1 * rhoJ**2)
        else
          wsor = 1.d0 / (1.d0 - 2.5d-1 * rhoJ**2 * wsor)
        endif

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j4,jcell1,jcell2,jdcell,inb1,inb2,inb3)
C$OMP+PRIVATE(phi,phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8)
        do j = 1 , nLevel
          j4 = j * 4   
          jcell1 = j4 + iblack
          jcell2 = j4 - 3

          phi1  = rhophi(2,jcell2)
c          phi2  = rhophi(2,jcell1)
c          phi3  = rhophi(2,jcell1+1)
          phi4  = rhophi(2,jcell2+1)
c          phi5  = rhophi(2,jcell1+2)
          phi6  = rhophi(2,jcell2+2)
          phi7  = rhophi(2,jcell2+3)
c          phi8  = rhophi(2,jcell1+3)
c
c....     2nd sibling
c         
          inb1   = iLevNb(jcell1,1)
          inb2   = iLevNb(jcell1,2)
          inb3   = iLevNb(jcell1,3)
          phi    = rhophi(2,inb1)
          phi    = phi + rhophi(2,inb2)
          phi    = phi + rhophi(2,inb3)

          rhophi(2,jcell1) = rhophi(2,jcell1) + 
     &            wsor6 * ( 
     &                           (phi + phi1 + phi4 + phi6) - 
     &                           six * rhophi(2,jcell1)
     &                          ) - trfi2 * rhophi(1,jcell1)
c
c....     3rd sibling
c
          jdcell = jcell1 + 1
          inb1   = iLevNb(jdcell,1)
          inb2   = iLevNb(jdcell,2)
          inb3   = iLevNb(jdcell,3)
          phi    = rhophi(2,inb1)
          phi    = phi + rhophi(2,inb2)
          phi    = phi + rhophi(2,inb3)

          rhophi(2,jdcell) = rhophi(2,jdcell) + 
     &            wsor6 * ( 
     &                           (phi + phi4 + phi1 + phi7) - 
     &                            six * rhophi(2,jdcell)
     &                          ) - trfi2 * rhophi(1,jdcell)
c
c....     5th sibling
c
          jdcell = jdcell + 1
          inb1   = iLevNb(jdcell,1)
          inb2   = iLevNb(jdcell,2)
          inb3   = iLevNb(jdcell,3)
          phi    = rhophi(2,inb1)
          phi    = phi + rhophi(2,inb2)
          phi    = phi + rhophi(2,inb3)

          rhophi(2,jdcell) = rhophi(2,jdcell) + 
     &           wsor6 * ( 
     &                           (phi + phi6 + phi7 + phi1) - 
     &                            six*  rhophi(2,jdcell)
     &                          ) - trfi2 * rhophi(1,jdcell)
c
c....     8th sibling
c
          jdcell = jdcell + 1
          inb1   = iLevNb(jdcell,1)
          inb2   = iLevNb(jdcell,2)
          inb3   = iLevNb(jdcell,3)
          phi    = rhophi(2,inb1)
          phi    = phi + rhophi(2,inb2)
          phi    = phi + rhophi(2,inb3)

          rhophi(2,jdcell) = rhophi(2,jdcell) +
     &            wsor6 * ( 
     &                           (phi + phi7 + phi6 + phi4) - 
     &                            six* rhophi(2,jdcell)
     &                          ) - trfi2 * rhophi(1,jdcell)
c
        enddo ! end j
        
        if ( iter .eq. 1 ) then 
          wsor = 1.d0 / (1.d0 - 5.d-1 * rhoJ**2)
        else
          wsor = 1.d0 / (1.d0 - 2.5d-1 * rhoJ**2 * wsor)
        endif

      ENDDO ! end DO iter
c
c.... write back to var
c      
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j4,jcell1,jcell2,icell)
      do j = 1 , nLevel 
        j4 = j*4
        jcell1 = j4 - 3
        jcell2 = j4 + iblack
        icell  = iSelect(j)
        var(2,icell)   = rhophi(2,jcell1)
        var(2,icell+1) = rhophi(2,jcell2)
        var(2,icell+2) = rhophi(2,jcell2+1)
        var(2,icell+3) = rhophi(2,jcell1+1)
        var(2,icell+4) = rhophi(2,jcell2+2)
        var(2,icell+5) = rhophi(2,jcell1+2)
        var(2,icell+6) = rhophi(2,jcell1+3)
        var(2,icell+7) = rhophi(2,jcell2+3)
      enddo

      return
      end
c
c     -------------------------------------------
      function Pyramide ( iCell , iChild , ivar )  
c     -------------------------------------------
c
c     purpose: interpolates on a boundary using pyramidal interpolation
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer iCell, iChild, ivar
c
      dimension iPyr(nchild,3)   ! interpolation pyramid vertices 
c
      data iPyr / 1, 2, 1, 2, 1, 2, 1, 2,   
     &            3, 3, 4, 4, 3, 3, 4, 4,
     &            5, 5, 5, 5, 6, 6, 6, 6  /
c
        iNb1   = iNb(iCell,iPyr(iChild,1))
        iNb2   = iNb(iCell,iPyr(iChild,2))
        iNb3   = iNb(iCell,iPyr(iChild,3))

        Pyramide = 0.25 * ( var(ivar,iCell) + 
     &                      var(ivar,iNb1) +
     &                      var(ivar,iNb2) +         
     &                      var(ivar,iNb3) )
c
      return
      end
c
c     --------------------------------------------
      function Pyramideh ( iCell , iChild , ivar )  
c     --------------------------------------------
c
c     purpose: interpolates on a boundary using pyramidal interpolation
c              operates on array hvar
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer iCell, iChild, ivar
c
      dimension iPyr(nchild,3)   ! interpolation pyramid vertices 
c
      data iPyr / 1, 2, 1, 2, 1, 2, 1, 2,   
     &            3, 3, 4, 4, 3, 3, 4, 4,
     &            5, 5, 5, 5, 6, 6, 6, 6  /
c
        iNb1   = iNb(iCell,iPyr(iChild,1))
        iNb2   = iNb(iCell,iPyr(iChild,2))
        iNb3   = iNb(iCell,iPyr(iChild,3))

       Pyramideh = 0.25 * ( hvar(ivar,iCell) + 
     &                       hvar(ivar,iNb1) +
     &                       hvar(ivar,iNb2) +         
     &                       hvar(ivar,iNb3) )
c
      return
      end
c
c     ---------------------------------------
      function Average_Int ( icell , icell0 )
c     ---------------------------------------
c
c     function interpolating on the boundaries.
c     algorithm: 
c               a boundary cell is weight-averaged  
c               with higher-level neighbor. weights
c               were found empirically by minimizing
c               self-action.
c     input:
c           icell  - index of lower-level neighb.
c           icell0 - index of a boundary cell 
c
      include 'a_tree.h'

      parameter ( wbig   = 0.8 ) ! weights
      parameter ( wsmall = 0.2 )

      Average_Int = wbig * var(2,icell) + wsmall * var(2,icell0)

      return
      end

c     --------------------
      subroutine Potent ()
c     --------------------
c
c	    Find potential on Grid FI=var(i,1):	DENSITY    ->	POTENTIAL
c
c		   O 1		    ^ - Fourier component
c		   |
c	     1	   |-4	 1	^      ^	2Pi
c	     O-----O-----O     Fi    =	Rho	/ (2cos(---  (i-1))+
c		   |	   i,j		i,j	Ngrid
c		   |
c		   O 1			  2Pi
c				       2cos(---  (j-1))-4)
c		       ^			Ngrid
c		       Fi	= 1 (?)
c			 11
c		   2
c		NABLA  Fi = 3/2  /A * (Rho - <Rho>) ;
c		   X
c			      <Rho> = 1
c
      include 'a_tree.h'
      include 'a_control.h'

      dimension       greenc(ng)
      dimension Zf(narr) , Yf(narr)
      dimension si(nf67) , index(nf67)

      integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

c			Set coeffitient in Poisson eq.
c			  (2*ng)**3 - from FFT

      trfi= 6.0 / aexp(MinLevel) / (2.*ng)**3
 
c.... set green function components

      p16 = 2.0 * pi / ng

      ngrid2 = ng/2 + 2
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(xx)
      do i = 1 , ng
	xx        = p16 * (i - 1) 
	greenc(i) = 2. * cos(xx)
	if (i .ge. ngrid2) greenc(i) = -greenc(i)
      enddo

C			Ngrid = 2**IQ

		 iq = int(alog(float(ng))/alog(2.)+0.5)
C$OMP PARALLEL DO DEFAULT(SHARED)
      do ic1 = 1 , ncell0
        var(2,ic1) = var(1,ic1) 
      enddo

      Ndex = 1
      ib1  = 3
C					  ALONG X-DIRECTION
      call setf67(ib1,iq,ibc,ip,isl,l1,n2,n3,n4,n7,si,index,Zf,Yf)
100   continue

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j,i,icell,zf,yf,ip,isl,l1,n2,n4)
      do k = 1 , ng
	    do j = 1 , ng
		  do i = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
                    zf(i) = var(2,icell)
		  enddo
 
		  call four67(ib1,iq,ip,isl,l1,n2,n7,
     &                 si,index,Zf,Yf)
 
		  do i = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
	            var(2,icell) = yf(i)
		  enddo
	    enddo
C					  ALONG Y-DIRECTION

	    do i = 1 , ng
		  do j = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
                    zf(j) = var(2,icell)
		  enddo
		  call four67(ib1,iq,ip,isl,l1,n2,n7,si,index,
     &                        Zf,Yf)
		  do j = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
		    var(2,icell) = yf(j)
		  enddo
	    enddo
      enddo
c					  EXIT IF IT IS THE SECOND LOOP
      if (Ndex .EQ. 2) return

C					  ALONG Z-DIRECTION
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i,k,icell,ib1,a1,a2,a3,zf,yf,ip,isl,l1,n2,n4)
      do j = 1 , ng
	  ib1 = 3
	  do i = 1 , ng
	    do k = 1 , ng
              icell = ((i-1)*ng + (j-1))*ng + k
              zf(k) = var(2,icell)
            enddo

	    call four67(ib1,iq,ip,isl,l1,n2,n7,
     &                  si,index,Zf,Yf)

            do k = 1 , ng
              icell        = ((i-1)*ng + (j-1))*ng + k
              var(2,icell) = yf(k)
            enddo
	    enddo
C					  BACK IN Z
	    a3  = greenc(j) - 6.0
	    ib1 = 4
	    do i = 1 , ng
		  a2 = greenc(i) + a3
		  do k = 1 , ng
	            a1 = a2 + greenc(k)
		    if (abs(a1) .lt. 1.e-4) a1 = 1.0
                    icell = ((i-1)*ng + (j-1))*ng + k
		    zf(k) = var(2,icell) * trfi / a1
		  enddo

	      call four67(ib1,iq,ip,isl,l1,n2,n7,
     &                    si,index,Zf,Yf)

		  do k = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
		    var(2,icell) = yf(k)
		  enddo
	    enddo
      enddo

      ib1  = 4 
      Ndex = 2
      goto 100

      end
