c     --------------------
c        poisson solver
c     --------------------

c     --------------------------------------------------
      subroutine Solve_Poisson ( MinModify , MaxModify , iProlongFlag )
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
         call Relax ( MinModify+1 , MaxModify , iProlongFlag )
        call Timing ( 3 ,  1 )
      else
        call Timing ( 3 , -1 )

         call Relax ( MinModify , MaxModify , iProlongFlag )
         call Timing ( 3 ,  1 )
      endif

      return
      end

c     ------------------------------------------
      subroutine Relax ( MinModify , MaxModify , iProlongFlag ) 
c     ------------------------------------------
c
c     purpose: main driver for relaxation iterations
c     input  : MinModify , MaxModify - levels to work on
c     note   : MinModify should be > MinLevel 
c              to ensure this is caller's responsibility
c
      include 'a_tree.h'
c
      iOctMax = moct
      do while ( iOctLv(iOctMax) .eq. iFreeLevel ) 
        iOctMax = iOctMax - 1
      enddo
      iOctMax = iOctMax + 1

      do Level = MinModify , MaxModify
         if ( iProlongFlag .eq. 1 ) then 
            call Prolongate ( Level-1 )
         endif
        call Smooth ( Level )
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

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( icell,var1,iOC,iChild_1)
        do icell = 1 , ncell0
c          iChild_1 = iCh(icell,1)
           iOC = iOctCh(icell)
           if ( iOC .le. nil ) then
             iChild_1 = nil
           else
             iChild_1 = ( iOC - 1 ) * nchild + 1 + ncell0
           end if
           if ( iChild_1 .gt. nil ) then ! if cell is split
            var1 = pot(icell)
            pot(iChild_1  )   = var1  
            pot(iChild_1+1) = var1 
            pot(iChild_1+2) = var1
            pot(iChild_1+3) = var1
            pot(iChild_1+4) = var1
            pot(iChild_1+5) = var1
            pot(iChild_1+6) = var1
            pot(iChild_1+7) = var1
          endif
        enddo
      ELSE

        call Select_Cells ( Level , nLevel ) 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( j,icell,idcell,pdens,iChild_1)
C$OMP+PRIVATE (iCh1,iCh2,iCh3,iCh4,iCh5,iCh6,iCh7,iCh8,var1)
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
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
          endif

c.... 2nd child
          idcell   = iCh2
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
          endif

c.... 3rd child
          idcell   = iCh3
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
          endif

c.... 4th child
          idcell   = iCh4
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
          endif

c.... 5th child
          idcell   = iCh5
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
          endif

c.... 6th child
          idcell   = iCh6
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
          endif

c.... 7th child
          idcell   = iCh7
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
          endif

c.... 8th child
          idcell   = iCh8
          iChild_1 = iCh(idcell,1)
          if ( iChild_1 .ne. nil ) then
            var1 = pot(idcell)
            pot(iChild_1   )   = var1
            pot(iChild_1 +1) = var1
            pot(iChild_1 +2) = var1
            pot(iChild_1 +3) = var1
            pot(iChild_1 +4) = var1
            pot(iChild_1 +5) = var1
            pot(iChild_1 +6) = var1
            pot(iChild_1 +7) = var1
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
c     algorythm: "pyramidal" interpolation is used
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
      pdens = wa * pot(iCell)      ! density contribution from parent

      do ic1 = 1 , nchild
        iChild = iCh(iCell,ic1)

        iNb1   = iNb(iCell,iPyr(ic1,1))
        iNb2   = iNb(iCell,iPyr(ic1,2))
        iNb3   = iNb(iCell,iPyr(ic1,3))

        pot(iChild) = pdens + 
     &                  wbcd  * ( pot(iNb1) +
     &                            pot(iNb2) +         
     &                            pot(iNb3) )

      enddo

      return
      end

c     -----------------------------
      subroutine Restrict ( Level ) 
c     -----------------------------
c
c     makes restriction sweep through Level
c
      include 'a_tree.h'

      if ( Level .le. MinLevel ) then
         print*, 'Restrict Error: cannot restrict level <= 0'
         stop
      endif

      nLevel = iNOLL(Level)         ! get level boundary for index array
      icell  = iGet_Cell_Index ( iHOLL(Level) )
         

      do j = 1 , nLevel
        icell = icell + nchild      
        sum   = zero
        do k = 1 , nchild
          idcell = icell - k
          sum    = sum   + pot(idcell)
        enddo
        iParent        = iPr(idcell)
        pot(iParent) = sum / nchild
        iOct           = iGet_Oct_Number ( idcell     ) 
        icell          = iGet_Cell_Index ( iOctLL1(iOct) )
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
      include 'a_mpi.h'
      include 'a_control.h'

c       common /REF03/  iLevNb(3,nctot-ncell0)
       real*4 wsor, wsor6, rhoJ
      real*4 Size,  trfi2
      real*4 phi, phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8


      logical ContinueSweep
      
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

      call Select_Cells ( Level , nLevel ) 

      Size  = (CellSize(Level))     ! cell size on this level      


C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) , SCHEDULE(DYNAMIC,2000) 
C$OMP+PRIVATE ( j,icell,idcell,jcell,jdcell )
          do j = 1 , nLevel
            icell = iSelect(j)
            jcell = (j-1)*8 + 1
c
c.... 1st sibling
              iLevNb(1,jcell) = iNb(icell,1)
              iLevNb(2,jcell) = iNb(icell,3)
              iLevNb(3,jcell) = iNb(icell,5)              
c.... 2nd sibling
              idcell = icell + 1
              jdcell = jcell + 1
              iLevNb(1,jdcell) = iNb(idcell,2)
              iLevNb(2,jdcell) = iNb(idcell,3)
              iLevNb(3,jdcell) = iNb(idcell,5)
c              write(13,3) j,icell,jcell,jdcell,
c     &         iLevNb(1,jcell),iLevNb(2,jcell),iLevNb(3,jcell),
c     &          iLevNb(1,jdcell),iLevNb(2,jdcell),iLevNb(3,jdcell) 
c 3            format(i7,3i8,3x,3i8,3x,3i8)
c.... 3rd sibling
              idcell = icell + 2
              jdcell = jcell + 2
              iLevNb(1,jdcell) = iNb(idcell,1)
              iLevNb(2,jdcell) = iNb(idcell,4)
              iLevNb(3,jdcell) = iNb(idcell,5)
c.... 4th sibling
              idcell = icell + 3
              jdcell = jcell + 3
              iLevNb(1,jdcell) = iNb(idcell,2)
              iLevNb(2,jdcell) = iNb(idcell,4)
              iLevNb(3,jdcell) = iNb(idcell,5)
c.... 5th sibling
              idcell = icell + 4
              jdcell = jcell + 4
              iLevNb(1,jdcell) = iNb(idcell,1)
              iLevNb(2,jdcell) = iNb(idcell,3)
              iLevNb(3,jdcell) = iNb(idcell,6)
c.... 6th sibling
              idcell = icell + 5
              jdcell = jcell + 5
              iLevNb(1,jdcell) = iNb(idcell,2)
              iLevNb(2,jdcell) = iNb(idcell,3)
              iLevNb(3,jdcell) = iNb(idcell,6)
c.... 7th sibling
              idcell = icell + 6
              jdcell = jcell + 6
              iLevNb(1,jdcell) = iNb(idcell,1)
              iLevNb(2,jdcell) = iNb(idcell,4)
              iLevNb(3,jdcell) = iNb(idcell,6)
c.... 8th sibling
              idcell = icell + 7
              jdcell = jcell + 7
              iLevNb(1,jdcell) = iNb(idcell,2)
              iLevNb(2,jdcell) = iNb(idcell,4)
              iLevNb(3,jdcell) = iNb(idcell,6)
          enddo ! end j
c
c.... pre-compute potential on the border
c     the following loop is to be executed SERIALLY
c
      if ( iOctMax .gt. 0 ) then          
        ibc = ncell0 + (iOctMax-1)*nchild
      endif
      if ( nctot - ibc .lt. 24 ) then 
        write(13,*) 'error: Smooth: not enough cells:'
        write(13,*) 'L =',Level,' ibc =',ibc,' nctot =',nctot
        write(13,*) 'increase mcell and rerun'
        write(name,'(a,i4.4,a)')'ERROR_RUN',irank+1,'.DAT'
        accesstype = 'append'
        call Open_ASCII_FILE(19,name,accesstype)
        rewind 19
        write(19,*) 
     +  ' 5  : ART_Gravity.f, error: Smooth: not enough cells:'
        write(19,*) 'L =',Level,' ibc =',ibc,' nctot =',nctot
        write(19,*) 'increase mcell and rerun'
        close(19) 
        stop

      endif
c
      do j = 1 , nLevel
        jcell = (j-1)*8 + 1
c
c....   1st sibling
c
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,2)
          iLevNb(1,jcell) = ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,3)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,5)
          iLevNb(3,jcell) = ibc
        endif
c
c....   2nd sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,1)
          iLevNb(1,jcell) = ibc
c          If(j.le.10)write(13,*)'     ',inb1,pot(ibc),ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,4)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,6)
          iLevNb(3,jcell) = ibc
        endif
c
c....   3rd sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,4)
          iLevNb(1,jcell) = ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,1)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,7)
          iLevNb(3,jcell) = ibc
        endif
c
c....   4th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,3)
          iLevNb(1,jcell) = ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,2)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,8)
          iLevNb(3,jcell) = ibc
        endif
c
c....   5th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,6)
          iLevNb(1,jcell) = ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,7)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,1)
          iLevNb(3,jcell) = ibc
        endif
c
c....   6th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,5)
          iLevNb(1,jcell) = ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,8)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,2)
          iLevNb(3,jcell) = ibc
        endif
c
c....   7th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,8)
          iLevNb(1,jcell) = ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,5)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,3)
          iLevNb(3,jcell) = ibc
        endif
c
c....   8th sibling
c
        jcell = jcell + 1
        inb1 = iLevNb(1,jcell)
        inb2 = iLevNb(2,jcell)
        inb3 = iLevNb(3,jcell)
        if (iLv(inb1) .ne. Level) then
          ibc = ibc + 1          
          pot(ibc) = Pyramide(inb1,7)
          iLevNb(1,jcell) = ibc
        endif
        if (iLv(inb2) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb2,6)
          iLevNb(2,jcell) = ibc
        endif
        if (iLv(inb3) .ne. Level) then
          ibc = ibc + 1
          pot(ibc) = Pyramide(inb3,4)
          iLevNb(3,jcell) = ibc
        endif
c
        if ( nctot - ibc .lt. 24 ) then 
          write(13,*) 'error: Smooth: not enough cells:'
          write(13,*) 'L =',Level,' nLevel =',nLevel,' j =',j
          write(13,*) ' ibc =',ibc,' nctot =',nctot
          write(13,*) 'increase mcell and rerun'
        write(name,'(a,i4.4,a)')'ERROR_RUN',irank+1,'.DAT'
          accesstype = 'append'
          call Open_ASCII_FILE(19,name,accesstype)
          rewind 19
          write(19,*) 
     +  ' 5  : ART_Gravity.f, error: Smooth: not enough cells:'
          write(19,*) 'L =',Level,' ibc =',ibc,' nctot =',nctot
          write(19,*) ' ibc =',ibc,' nctot =',nctot
          write(19,*) 'increase mcell and rerun'
          close(19) 
          stop
        endif
      enddo

C............................................................
      wsor = 1.d0
      rhoJ = 0.9995d0
      Size2 = Size * Size

      DO iter = 1 , niter

        wsor6 = wsor / 6.d0  
        trfi2 = 0.25 * Om0 * wsor  / aexp(Level)/Size
c
c.....    split sweep in 2 (odd/even) parts

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(DYNAMIC,2000) 
C$OMP+PRIVATE ( j,icell,idcell,phi,inb1,inb2,inb3)
C$OMP+PRIVATE ( phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8)
C$OMP+PRIVATE ( jcell,jdcell )
          do j = 1 , nLevel
            icell = iSelect(j)
              phi1 = pot(icell)
              phi2 = pot(icell+1)
              phi3 = pot(icell+2)
              phi4 = pot(icell+3)
              phi5 = pot(icell+4)
              phi6 = pot(icell+5)
              phi7 = pot(icell+6)
              phi8 = pot(icell+7)
              jcell = (j-1)*8 + 1
c
c.... 2nd sibling
c
              idcell = icell + 1
              jdcell = jcell + 1
              inb1 = iLevNb(1,jdcell)
              inb2 = iLevNb(2,jdcell)
              inb3 = iLevNb(3,jdcell)
              phi = pot(inb1) + pot(inb2) + pot(inb3)
              pot(idcell) = pot(idcell) + 
     &            wsor6 * ( 
     &                           (phi + phi1 + phi4 + phi6) - 
     &                           6.d0* pot(idcell)
     &                          ) - trfi2 * var(1,idcell)
c         If(j.le.10)write(13,2)iter,j,idcell, 
c     &           phi,pot(idcell)
 2       format(i3,2i8,3g11.5,g11.5)
c
c.... 3rd sibling
c
              idcell = icell + 2
              jdcell = jcell + 2
              inb1 = iLevNb(1,jdcell)
              inb2 = iLevNb(2,jdcell)
              inb3 = iLevNb(3,jdcell)
              phi =  pot(inb1)+ pot(inb2)+ pot(inb3)
              pot(idcell) = pot(idcell) + 
     &            wsor6 * ( 
     &                           (phi + phi4 + phi1 + phi7) - 
     &                            6.* pot(idcell)
     &                          ) - trfi2 * var(1,idcell)
c
c.... 5th sibling
c
              idcell = icell + 4
              jdcell = jcell + 4
              inb1 = iLevNb(1,jdcell)
              inb2 = iLevNb(2,jdcell)
              inb3 = iLevNb(3,jdcell)
              phi =  pot(inb1)+ pot(inb2)+ pot(inb3)
              pot(idcell) =  pot(idcell) + 
     &           wsor6 * ( 
     &                           (phi + phi6 + phi7 + phi1) - 
     &                            6.*  pot(idcell)
     &                          ) - trfi2 * var(1,idcell)
c
c.... 8th sibling
c
              idcell = icell + 7
              jdcell = jcell + 7
              inb1 = iLevNb(1,jdcell)
              inb2 = iLevNb(2,jdcell)
              inb3 = iLevNb(3,jdcell)
              phi = pot(inb1) + pot(inb2) + pot(inb3)
              pot(idcell) = pot(idcell) +
     &            wsor6 * ( 
     &                           (phi + phi7 + phi6 + phi4) - 
     &                            6.* pot(idcell)
     &                          ) - trfi2 * var(1,idcell)
c
          enddo ! end j
        
        if ( iter .eq. 1 ) then 
          wsor = 1.d0 / (1.d0 - 5.d-1 * rhoJ**2)
        else
          wsor = 1.d0 / (1.d0 - 2.5d-1 * rhoJ**2 * wsor)
        endif

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED)  ! , SCHEDULE(DYNAMIC,32)  
C$OMP+PRIVATE ( j,icell,idcell,phi,inb1,inb2,inb3)
C$OMP+PRIVATE ( phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8 )
C$OMP+PRIVATE ( jcell,jdcell)
          do j = 1 , nLevel
            icell = iSelect(j)
              phi1 = pot(icell)
              phi2 = pot(icell+1)
              phi3 = pot(icell+2)
              phi4 = pot(icell+3)
              phi5 = pot(icell+4)
              phi6 = pot(icell+5)
              phi7 = pot(icell+6)
              phi8 = pot(icell+7)
              jcell = (j-1)*8 + 1
c
c.... 1st sibling
c
              inb1 = iLevNb(1,jcell)
              inb2 = iLevNb(2,jcell)
              inb3 = iLevNb(3,jcell)              
              phi = pot(inb1) + pot(inb2) + pot(inb3)
              pot(icell) =  pot(icell) + 
     &           wsor6 * ( 
     &                           (phi + phi2 + phi3 + phi5) - 
     &                           6.*pot(icell)
     &                         ) - trfi2 * var(1,icell)
c   
c.... 4th sibling
c
              idcell = icell + 3
              jdcell = jcell + 3
              inb1 = iLevNb(1,jdcell)
              inb2 = iLevNb(2,jdcell)
              inb3 = iLevNb(3,jdcell)
              phi =  pot(inb1)+ pot(inb2) + pot(inb3)
              pot(idcell) =   pot(idcell) + 
     &            wsor6 * ( 
     &                           (phi + phi3 + phi2 + phi8) - 
     &                            6.*  pot(idcell)
     &                          ) - trfi2 * var(1,idcell)
c
c.... 6th sibling
c
              idcell = icell + 5
              jdcell = jcell + 5
              inb1 = iLevNb(1,jdcell)
              inb2 = iLevNb(2,jdcell)
              inb3 = iLevNb(3,jdcell)
              phi =  pot(inb1) + pot(inb2) + pot(inb3)
              pot(idcell) = pot(idcell) + 
     &           wsor6 * ( 
     &                           (phi + phi5 + phi8 + phi2) - 
     &                            6.* pot(idcell)
     &                          ) - trfi2 * var(1,idcell)
c
c.... 7th sibling
c
              idcell = icell + 6
              jdcell = jcell + 6
              inb1 = iLevNb(1,jdcell)
              inb2 = iLevNb(2,jdcell)
              inb3 = iLevNb(3,jdcell)
              phi = pot(inb1)+ pot(inb2)+ pot(inb3)
              pot(idcell) =  pot(idcell) + 
     &           wsor6 * ( 
     &                           (phi + phi8 + phi5 + phi3) - 
     &                            6.*  pot(idcell)
     &                          ) - trfi2 * var(1,idcell)
c
          enddo ! end j
        
        if ( iter .eq. 1 ) then 
          wsor = 1.d0 / (1.d0 - 5.d-1 * rhoJ**2)
        else
          wsor = 1.d0 / (1.d0 - 2.5d-1 * rhoJ**2 * wsor)
        endif


      ENDDO ! end DO WHILE
c      write (13,*) '               out'
c      write (13,*) ' Level=',Level,nlevel,' NPyr=',nPyr,' NiNb=',niNb
      return
      end
c
c
c     ------------------------------------
      function Pyramide ( iCell , iChild )  
c     ------------------------------------
c
c     purpose: interpolates on a boundary using pyramidal interpolation
c
      include 'a_control.h'
      include 'a_tree.h'
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

        Pyramide = 0.25 * ( pot(iCell) + 
     &                      pot(iNb1) +
     &                      pot(iNb2) +         
     &                      pot(iNb3) )
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

      Average_Int = wbig * pot(icell) + wsmall * pot(icell0)

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
      dimension si(nf67) , indx(nf67)

      integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

      double precision dumb, dumb1(ng)

c			Set coeffitient in Poisson eq.
c			  (2*ng)**3 - from FFT

      trfi=1.5 * Om0 / aexpn / (2.*ng)**3
 
c.... set green function components

      P16 = 2.0 * pi

      ngrid2 = ng/2 + 2
      do i = 1 , ng
        xx        = p16 * (i - 1.) / ng
        greenc(i) = 2. * cos(xx)
        if (i .ge. ngrid2) greenc(i) = -greenc(i)
      enddo

C			Ngrid = 2**IQ

		 iq = int(alog(float(ng))/alog(2.)+0.5)

      dumb = 0.
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1) REDUCTION(+:dumb)
      do ic1 = 1 , ncell0
        pot(ic1) = var(1,ic1) 
        dumb = dumb + var(1,ic1)
      enddo

      Ndex = 1
      ib1  = 3
C					  ALONG X-DIRECTION
      call setf67(ib1,iq,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
100   continue

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,icell,zf,yf,ip,isl,l1,n2,n4 )
      do k = 1 , ng
	    do j = 1 , ng
		  do i = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
                    zf(i) = pot(icell)
		  enddo
 
		  call four67(ib1,iq,ip,isl,l1,n2,n7,
     &                 si,indx,Zf,Yf)
 
		  do i = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
	            pot(icell) = yf(i)
		  enddo
	    enddo
C					  ALONG Y-DIRECTION

	    do i = 1 , ng
		  do j = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
                    zf(j) = pot(icell)
		  enddo
		  call four67(ib1,iq,ip,isl,l1,n2,n7,si,indx,
     &                        Zf,Yf)
		  do j = 1 , ng
                    icell = ((i-1)*ng + (j-1))*ng + k
		    pot(icell) = yf(j)
		  enddo
	    enddo
      enddo
c					  EXIT IF IT IS THE SECOND LOOP
      if (Ndex .EQ. 2) return

C					  ALONG Z-DIRECTION
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (j,i,k,icell,ib1,a1,a2,a3,zf,yf,ip,isl,l1,n2,n4 )
      do j = 1 , ng
        ib1 = 3
        do i = 1 , ng
          do k = 1 , ng
            icell = ((i-1)*ng + (j-1))*ng + k
            zf(k) = pot(icell)
          enddo

	       call four67(ib1,iq,ip,isl,l1,n2,n7,
     &                  si,indx,Zf,Yf)

          do k = 1 , ng
            icell        = ((i-1)*ng + (j-1))*ng + k
            pot(icell) = yf(k)
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
            zf(k) = pot(icell) * trfi / a1
          enddo

          call four67(ib1,iq,ip,isl,l1,n2,n7,
     &                    si,indx,Zf,Yf)

          do k = 1 , ng
            icell = ((i-1)*ng + (j-1))*ng + k
            pot(icell) = yf(k)
            dumb1(k) = var(1,icell)
		    enddo
	     enddo
      enddo

      ib1  = 4 
      Ndex = 2
      goto 100

      end
