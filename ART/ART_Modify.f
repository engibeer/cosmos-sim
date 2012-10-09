c     =====================================================================
c                                                                         .
c              ART Version 3:  Refinement/Derefinement                    .
c                                                                         .
c                      by Andrey Kravtsov (1997)                          .
c                                                                         .   
c     =====================================================================

c     ----------------------------------------------------------------
      subroutine Modify ( mtot , MinModify , MaxModify , MaxL , iJob )
c     ----------------------------------------------------------------
c
c     purpose: Makes refinement/derefinement on a tree
c
c     input  : MinModify , MaxModify - min & max level to process
c              MaxL - refinement are not allowed beyond this level
c              iJob - kind of job to be done
c                     1 - refinement only
c                     0 - refinement & derefinement
c                    -1 - derefinement only
c
c     output : mtot - total number of refined/derefined cells
c          
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_mpi.h'
      integer MinModify , MaxModify , MaxL 
      integer mtot         ! total # of cell marked to split 
      integer nmark        ! # of marked cells

      mtot  = nil 
      nmark = nil 

      call Mark_Cells ( nmark , MinModify , MaxModify , MaxL )

      if ( iJob .eq. 1 ) then
        call Refine ( mtot , MinModify , MaxModify , MaxL )
      else
        if ( iJob .eq. 0 ) then
          call DeRefine ( mtot , MinModify , MaxModify )          
          call Refine   ( mtot , MinModify , MaxModify , MaxL )
        else
          if ( iJob .eq. -1 ) then
          if ( MinModify .gt. MinLevel ) then
            call DeRefine ( mtot , MinModify , MaxModify )     
          endif
          else
            write(13,*) 'Modify error: bad iJob:',iJob
            name = 'ERROR_RUN'
            accesstype = 'append'
            call Open_ASCII_FILE(19,name,accesstype)
            rewind 19
            write(19,*) 
     +      ' 3  : ART_Modify.f,  Modify error: bad iJob:',iJob
            close(19)
            call MPI_Abort(MPI_COMM_WORLD,ierr1,ierr2)
            stop
          endif
        endif
      endif

      return
      end

c     --------------------------------------------------------------
      subroutine Mark_Cells ( nmark , MinModify , MaxModify , MaxL )
c     --------------------------------------------------------------
c
c     Computes weights for refinement/derefinement
c     according to distribution of DM particles
c
c.... get refinement/derefinement indicators for specified levels
c
      do Level = MaxModify , MinModify , -1
        call Zero_Indicators ( Level )
      enddo

      do Level = min(MaxModify,MaxL) , MinModify , -1
        call Smooth_Refinement_Indicators ( Level )
      enddo 

      return
      end
c     ------------------------------------
      subroutine Zero_Indicators ( Level ) 
c     ------------------------------------
c
c     purpose: to prepare working space for 
c              marking. var(1,*) does not have to be zeroed 
c              because it is used to construct original indicators
c
      include 'a_control.h'
      include 'a_tree.h'

      if ( Level .eq. MinLevel ) then 
        wconst = wsplit / (trho(MinLevel)+1.0)
      else
        wconst = wsplit / trho(Level)
      endif

      nmark  = nil 

      IF ( Level .eq. MinLevel ) THEN 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1)
        do ic1 = 1 , ncell0 
          var(1,ic1) = wconst * var(2,ic1)
          var(2,ic1) = zero 
        enddo

      ELSE

        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1,icell,idcell,ic2)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            var(1,idcell) = wconst * var(2,idcell)
            var(2,idcell) = zero
          enddo
        enddo

      ENDIF

      return
      end

c     ------------------------------------
      subroutine Zero_Indicators_fixed ( Level ) 
c     ------------------------------------
c
c     purpose: to prepare working space for 
c              marking. var(1,*) does not have to be zeroed 
c              because it is used to construct original indicators
c
      include 'a_tree.h'
      include 'a_control.h'
      logical TestRefinementZone

      wconst = wsplit / trho(Level)
      nmark  = nil 

      IF ( Level .eq. MinLevel ) THEN 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1 ) 
        do ic1 = 1 , ncell0 
          if ( TestRefinementZone(ic1) ) then 
            var(1,ic1) = wsplit
          else
            var(1,ic1) = zero
          endif
          var(2,ic1) = zero 
        enddo

      ELSE

        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1,icell,idcell,ic2)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            if ( TestRefinementZone(idcell) ) then 
              var(1,idcell) = wconst * var(2,idcell)
            else
              var(1,idcell) = zero
            endif
            var(2,idcell) = zero 
          enddo
        enddo

      ENDIF

      return
      end

c     ------------------------------------------
      logical function TestRefinementZone ( iC )
c     ------------------------------------------
c
c     test if cell iC is located within the refinement zone
c
      include 'a_tree.h'
      integer iC , i , j , k
      real Posx, Posy, Posz

      call Ps ( iC , Posx , Posy , Posz )
      i = int(Posx)
      j = int(Posy)
      k = int(Posz)
      TestRefinementZone = .true.
      if ( (i.lt.imin) .or. (i.gt.imax) ) then 
        TestRefinementZone = .false.
        return
      endif
      if ( (j.lt.jmin) .or. (j.gt.jmax) ) then 
        TestRefinementZone = .false.
        return
      endif
      if ( (k.lt.kmin) .or. (k.gt.kmax) ) then
        TestRefinementZone = .false.
        return
      endif

      return
      end

c     -------------------------------------------------
      subroutine Smooth_Refinement_Indicators ( Level )
c     -------------------------------------------------
c
c     purpose     : smoothes refinement indices on a given Level
c     input       : Level - level to smooth
c
      include 'a_tree.h'
      include 'a_control.h'

      parameter ( nPass   = 6 )  ! how many passes to make
      dimension nSmooth(nPass)

      data nSmooth / 1 , 2 , 2 , 1 , 2 , 2 /  ! smoothing mode for every pass

      IF ( Level .eq. MinLevel ) THEN 
        do ic0 = 1 , nPass
          n1 = mod(ic0+1,2) + 1    ! index of working array column 
          n2 = mod(ic0+2,2) + 1    ! index of auxiliary array column 
          if ( nSmooth(ic0) .eq. 1 ) then 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1,iDir)
            do ic1 = 1 , ncell0
              if ( var(n1,ic1) .ge. wsplit ) then 
                var(n2,ic1) = var(n1,ic1)
                do iDir = 1 , neighb 
                  var(n2,iNb(ic1,iDir)) = wsplit
                enddo               
              endif
            enddo ! end ic1

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1)
            do ic1 = 1 , ncell0
              var(n1,ic1) = zero 
            enddo

          else 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1,iDir,Neighbor)
            do ic1 = 1 , ncell0
              if ( var(n1,ic1) .ge. wsplit ) then 
                var(n2,ic1) = var(n1,ic1)
                do iDir = 1 , neighb 
                  Neighbor = iNb(ic1,iDir)
                  var(n2,Neighbor) = var(n2,Neighbor) + 
     &                               0.5 * wsplit + 0.05
                enddo               
              endif              
            enddo ! end ic1

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1)
            do ic1 = 1 , ncell0
              var(n1,ic1) = zero 
              if ( var(n2,ic1) .lt. wsplit ) var(n2,ic1) = zero 
            enddo

          endif
        enddo ! end ic0
      ELSE                 ! if Level .ne. MinLevel

        nLevel = iNOLL(Level)       
        call Select_Cells ( Level , nLevel ) 
        do ic0 = 1 , nPass
          n1 = mod(ic0+1,2) + 1    ! index of working array column 
          n2 = mod(ic0+2,2) + 1    ! index of auxiliary array column 
          IF ( nSmooth(ic0) .eq. 1 ) THEN

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( j , icell , k , idcell , iDir )
            do j = 1 , nLevel
              iCell = iSelect(j)
              do k = 0 , nchild-1
                idcell   = icell + k 
                if ( var(n1,idcell) .ge. wsplit ) then 
                  var(n2,idcell) = var(n1,idcell)
                  do iDir = 1 , neighb 
                    var(n2,iNb(idcell,iDir)) = wsplit
                  enddo               
                endif
              enddo ! end k
            enddo ! end j

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (  j , icell , k , idcell )
            do j = 1 , nLevel 
              iCell = iSelect(j)
              do k = 0 , nchild-1
                idcell = icell + k 
                var(n1,idcell) = zero
              enddo
            enddo ! end j

          ELSE

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( j , icell , k , idcell , iDir , Neighbor ) 
            do j = 1 , nLevel
              iCell = iSelect(j)
              do k = 0 , nchild-1
                idcell = icell + k 
                if ( var(n1,idcell) .ge. wsplit ) then 
                  var(n2,idcell) = var(n1,idcell)
                  do iDir = 1 , neighb 
                    Neighbor = iNb(idcell,iDir)
                    var(n2,Neighbor) = var(n2,Neighbor) + 
     &                                 0.5 * wsplit + 0.05
                  enddo               
                endif              
              enddo ! end k
            enddo ! end j

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (  j , icell , k , idcell )
            do j = 1 , nLevel 
              iCell = iSelect(j)
              do k = 0 , nchild-1
                idcell = icell + k 
                var(n1,idcell) = zero
                if ( var(n2,idcell) .lt. wsplit ) var(n2,idcell) = zero 
              enddo
            enddo ! end j

          ENDIF
        enddo ! end ic0
      ENDIF

      return
      end

c     ---------------------------------------------------------
      subroutine Refine ( mtot , MinModify , MaxModify , MaxL )
c     ---------------------------------------------------------
c
c     purpose: Makes refinement/derefinement on a tree
c     input  : MinModify , MaxModify - min & max level to process
c     output : mtot - total number of refined/derefined cells
c
c     uses   : Zero_Final_Weights, 
c              DM_Indicators, 
c              Split_Cells, 
c     var(*,1) - is used to mark cells for splitting & joining
c
      integer MinModify , MaxModify , MaxL
      integer mtot         ! total # of cell marked to split 

c.... split all marked cells       
      call Split ( mtot , MinModify , MaxModify , MaxL )

      return
      end

c     --------------------------------------------------------
      subroutine Split ( mtot , MinModify , MaxModify , MaxL )
c     --------------------------------------------------------
c
c     purpose: splits cells marked to split
c
      include 'a_tree.h'
      include 'a_control.h'
c      include 'a_sfc.h'

      integer idummy , idcell , mtot
      integer MinMod, MaxMod
      common / MinMaxMod / MinMod, MaxMod

      MinMod = MinModify
      MaxMod = MaxModify

c.... split marked cells


      do Level = min(MaxModify,MaxL) , MinModify , -1
        IF ( Level .eq. MinLevel ) THEN
c           write (13,*) ' filling curve'
         do ic1 = 1 , ncell0 
c          do icell = 1 , ncell0                        !! filling curve
c             ic1 = indx_sort(icell)                  !!
            if ( iOctCh(ic1) .eq. nil ) then 
              if (var(1,ic1)  .ge. wsplit ) then          
                if ( MaxL .ne. MinLevel ) then 
                  mtot = mtot + 1
                  call Force_Split_Cell ( ic1 , ires )
                endif
              endif
            endif
          enddo          
        ELSE
          nLevel = iNOLL(Level)        ! get level boundary for index array
          icell  = iGet_Cell_Index ( iHOLL(Level) ) 
          do j = 1 , nLevel
            icell = icell + nchild    ! index array contains only 1st child
            do k = 1 , nchild
              idcell = icell - k 
              if ( iOctCh(idcell) .eq. nil ) then    ! if leaf
                if ( var(1,idcell) .ge. wsplit ) then ! if marked to split
                  if ( Level .lt. MaxL ) then
                    mtot = mtot + 1
                    call Force_Split_Cell ( idcell , ires )
                  endif
                endif
              endif
            enddo
            iOct  = iGet_Oct_Number ( idcell     )
            icell = iGet_Cell_Index ( iOctLL1(iOct) )
          enddo
        ENDIF
      enddo

      MinModify = MinMod
      MaxModify = MaxMod

      return
      end

c     --------------------------------------------------
      subroutine DeRefine ( mtot , MinModify , MaxModify )
c     --------------------------------------------------
c
c     purpose: Makes derefinement on a tree
c     input  : MinModify , MaxModify - min & max level to process
c     output : mtot - total number of refined/derefined cells
c     uses   : Join
c     var(*,1) - is used to mark cells for splitting & joining
c
      integer mtot         ! total # of cell marked to split 

c.... join all marked cells 

      do Level = MaxModify , MinModify , -1
        call Join ( mtot , Level )
      enddo

      return
      end

c     --------------------------------
      subroutine Join ( mtot , Level )
c     --------------------------------
c
c     purpose: splits cells marked to split
c
      include 'a_tree.h'
      include 'a_control.h'

      integer idummy

c.... join marked cells

      IF ( Level .eq. MinLevel ) THEN
        do j = 1 , ncell0
          idcell = j
          if ( iOctCh(idcell) .gt. nil ) then 
            if ( var(1,idcell) .lt. wjoin ) then 
              nleaf = nil
              iCh1  = iCh(idcell,1) + nchild
              do ic1 = 1 , nchild
                iChild = iCh1 - ic1
                if ( (iOctCh(iChild) .eq. nil)
     &                          .and. 
     &               (var(1,iChild) .lt. wjoin) ) then
                  nleaf = nleaf + 1
                endif
              enddo
           
              if ( nleaf .eq. nchild ) then 
                idummy = iJoinCell ( idcell ) 
                if ( idummy .eq. nil ) mtot = mtot + 1
              endif         
            endif
          endif
        enddo
      ELSE
        nLevel = iNOLL(Level)        ! get level boundary for index array
        icell  = iGet_Cell_Index ( iHOLL(Level) )
        do j = 1 , nLevel
          icell = icell + nchild        ! index array contains only 1st child
          do k = 1 , nchild
            idcell = icell - k
            if ( iOctCh(idcell) .gt. nil ) then 
              if ( var(1,idcell) .lt. wjoin ) then 
                nleaf = nil
                iCh1  = iCh(idcell,1) + nchild
                do ic1 = 1 , nchild
                  iChild = iCh1 - ic1
                  if ( (iOctCh(iChild) .eq. nil)
     &                          .and.
     &                 (var(1,iChild) .lt. wjoin) ) then
                    nleaf = nleaf + 1
                  endif
                enddo
                if ( nleaf .eq. nchild ) then 
                  idummy = iJoinCell ( idcell ) 
                  if ( idummy .eq. nil ) mtot = mtot + 1
                endif
              endif
            endif
          enddo
          iOct  = iGet_Oct_Number ( idcell     )
          icell = iGet_Cell_Index ( iOctLL1(iOct) )
        enddo

      ENDIF

      return
      end
      
c     --------------------------------------------
      subroutine Force_Split_Cell ( icell , ires )
c     --------------------------------------------
c
c     purpose: splits cell splitting its neighbors if required
c
c     input  : icell - index of cell to split
c     output : ires  - success or not
c     note   : to assure that icell is leaf - is caller's responsibility
c
      include 'a_tree.h'
      include 'a_mpi.h'
      integer iCell , ires

 620  idummy0 = icell
      ires =  iSplitCell ( idummy0 )

      if ( ires .eq. -4 ) then   ! if neighbor rule is violated

 630     continue

         Level = iLv(idummy0)
         do ic2 = 1 , neighb           ! sweep through the neighbors
           inbdum = iNb(idummy0,ic2)
           if (iLv(inbdum).lt. Level) then
             idummy1 = inbdum
             ires = iSplitCell(idummy1)
             if (ires .eq. -4) then
               idummy0 = idummy1
               goto 630           
             endif

           endif
         enddo

         goto 620

       else 
         if ( (ires .ne. 0) .and. (ires .ne. -1) ) then
           write(13,*)'* Force Split Cell error: ires - ',ires
           close (13)
           call MPI_Abort(MPI_COMM_WORLD,ierr1,ierr2)
           CALL mpi_finalize(ierr)
           STOP
         end if         
       endif
            
      return
      end

c     -------------------------
      function iSplitCell ( i )            
c     -------------------------
c
c     Purpose: splits cell i
c
c     Returns:          0  - success
c                      -1  - already split
c                      -2  - run out of free Octs
c                      -4  - +/- 1 rule is violated
c
      include 'a_tree.h'
      dimension iENb(8,3)               ! external neighbors for given child
      data iENb    / 1, 2, 1, 2, 1, 2, 1, 2,       ! "x"
     &               3, 3, 4, 4, 3, 3, 4, 4,       ! "y"
     &               5, 5, 5, 5, 6, 6, 6, 6  /     ! "z"
      integer MinMod, MaxMod
      common / MinMaxMod / MinMod, MaxMod

      if ( iOctCh(i) .gt. 0 ) then
        iSplitCell = -1                   ! already split
        return
      end if
      NewOct = iGetOct()
      if ( NewOct .eq. nil ) then
        iSplitCell = -2                   ! no Octs left
           write(13,*)'* Force Split Cell error: ires - ',ires                    |     &               5, 5, 5, 5, 6, 6, 6, 6  /     ! "z"
           close (13)                                                             |        iSplitCell = -2                   ! no Octs left
           call MPI_Abort(MPI_COMM_WORLD,ierr1,ierr2)                             |           write(13,*)'* Force Split Cell error: ires - ',ires                   $
           CALL mpi_finalize(ierr)
        return
      end if
      Level = iLv(i)
c.... check for possible neighbor rule violations
      if ( Level .gt. MinLevel ) then 
        iChN = mod ( i - ncell0 - 1 , nchild ) + 1
        do ic2 = 1 , 3        ! over external neighbors
          iDir   = iENb(iChN,ic2)
          idummy = iNb (i,iDir)
          if ((iOctCh(idummy).eq.nil).and.(iLv(idummy).lt.Level)) then
            iSplitCell = -4
            return          
          endif
        enddo                   
      endif
      iOctCh(i) = NewOct               ! pointer from i to NewOct 
      iOctPr(NewOct) = i               ! parent of NewOct
      iOctLv(NewOct) = Level + 1       ! level of NewOct

      if  ( MinMod .gt. Level+1 ) MinMod = Level + 1
      if  ( MaxMod .lt. Level+1 ) MaxMod = Level + 1

      do k = 1 , nchild                ! cells of NewOct = leaves
        iOctCh(ncell0+(NewOct-1)*nchild+k) = 0  
      end do
      if ( i .le. ncell0 ) then 
        ic = (i - 1)/ng2 + 1
        jc = (i - (ic-1)*ng2 - 1)/ng + 1
        kc =  i - (ic-1)*ng2 - (jc-1)*ng
        id = ishft ( 1 , MaxLevel-1 )
        iOctPs(NewOct,1) = ic * isize + id
        iOctPs(NewOct,2) = jc * isize + id
        iOctPs(NewOct,3) = kc * isize + id
        do k = 1 , nneigh                   ! NewOct-neighbours = i-neighbours
          iOctNb(NewOct,k) = iNb(i,k)
        end do
      else
        iO = ishft ( i + nbshift , - ndim )
        id = ishft ( 1 , MaxLevel - iOctLv(iO) - 1 )  
        j  = mod ( i - ncell0 - 1 , nchild ) + 1
        do k = 1 , ndim
          iposik = iOctPs(iO,k) + sign ( id , idelta(j,k) ) 
          iOctPs(NewOct,k) = iposik         ! NewOct-coordinates
        end do
        do k = 1 , nneigh                   ! NewOct-neighbours = i-neighbours
          iOctNb(NewOct,k) = iNb(i,k)
        end do
      end if
c
      call iKy_Insert ( NewOct )        ! update cell linked list
      call LL_Split   (   i    )        ! split particle linked list
c
      do j = 1 , nchild 
        iChild = ( NewOct - 1 ) * nchild + j + ncell0
        pot(iChild) = Pyramide ( i , j )
      enddo
c
c.... successful split
c
      iSplitCell = 0
c
      return
      end



