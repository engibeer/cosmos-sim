c     =====================================================================
c                                                                         
c          ART Version 3-1 : routines copying old values of variables
c     at time t to the storage for new variables at time t+dt and vice versa
c                                                                         
c                         Andrey Kravtsov (1996-1999)           
c                                                                         
c     =====================================================================
c
c
c
c     --------------------------------------
      subroutine H_Old_to_New ( Level, Key )
c     --------------------------------------

c     Moves physical variables from old (var) to new (var) arrays.
c     Input: Level - level of leaves to move
c            Key   - 0   moves all leaves
c                    1   moves leaves that do     have split neighbours
c                    2   moves leaves that do not have split neighbours
c 
      include 'a_numbers.h'
c
      integer Level , Key
c
      external H_Old2New_1,
     &         H_Old2New_2,
     &         H_Old2New_3
c
      if ( Key .eq. 0 ) then
        call H_Old2New_1 ( Level )
      else if ( Key .eq. 1 ) then
        call H_Old2New_2 ( Level ) 
      else if ( Key .eq. 2 ) then
        call H_Old2New_3 ( Level )
      end if
c
      return
      end
c
c     --------------------------------
      subroutine H_Old2New_1 ( Level )
c     --------------------------------
c
#     include "a_def.h"
      include 'a_tree.h'

      integer Level
      integer i1, i2, j, idcell, icell, nLevel
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j)
        do i1 = 1 , ncell0 
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves
            vnw(1,i1) = hvar(1,i1)
            vnw(2,i1) = hvar(2,i1)
            vnw(3,i1) = hvar(3,i1)
            vnw(4,i1) = hvar(4,i1)
            vnw(5,i1) = hvar(5,i1)
            vnw(6,i1) = hvar(8,i1)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,i1) = hvar(j,i1)
            enddo
#endif
          endif
        enddo        
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,j)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
            vnw(1,icell) = hvar(1,icell)
            vnw(2,icell) = hvar(2,icell)
            vnw(3,icell) = hvar(3,icell)
            vnw(4,icell) = hvar(4,icell)
            vnw(5,icell) = hvar(5,icell)
            vnw(6,icell) = hvar(8,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,icell) = hvar(j,icell)
            enddo
#endif
          endif
        enddo           
      ENDIF
c
      return
      end
c
c     ----------------------------------
      subroutine H_Old2New_2 ( Level )
c     ----------------------------------
c
c     hvar --> vnw for leaves with split neighbors, 
c     in addition zeroes ref(*) of these cells 
c
#     include "a_def.h"
      include 'a_tree.h'
      integer Level
      integer i1, i2, j, idcell, icell, nLevel, Neighbors(6)
      logical split
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j,split,Neighbors)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then   ! only leaves
            split = .false. 
            call NbAll0 ( i1 , Neighbors )
            if ( iOctCh(Neighbors(1)) .gt. nil ) then 
              split = .true. 
              goto 10
            endif
            if ( iOctCh(Neighbors(2)) .gt. nil ) then 
              split = .true. 
              goto 10
            endif
            if ( iOctCh(Neighbors(3)) .gt. nil ) then 
              split = .true. 
              goto 10
            endif
            if ( iOctCh(Neighbors(4)) .gt. nil ) then 
              split = .true. 
              goto 10
            endif
            if ( iOctCh(Neighbors(5)) .gt. nil ) then 
              split = .true. 
              goto 10
            endif
            if ( iOctCh(Neighbors(6)) .gt. nil ) then 
              split = .true. 
            endif
 10         if ( split ) then 
              vnw(1,i1) = hvar(1,i1)
              vnw(2,i1) = hvar(2,i1)
              vnw(3,i1) = hvar(3,i1)
              vnw(4,i1) = hvar(4,i1)
              vnw(5,i1) = hvar(5,i1)
              vnw(6,i1) = hvar(8,i1)
              ref(i1)   = 0.0
#ifdef ADVECT_SPECIES
              do j = 9 , nhvarA
                vnw(j-2,i1) = hvar(j,i1)
              enddo
#endif
            endif
          endif
        enddo        
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,j,i2,idcell,split,Neighbors)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            if ( iOctCh(idcell) .eq. nil ) then 
              split = .false. 
              call NbAllL ( idcell , Neighbors ) 
              if ( iOctCh(Neighbors(1)) .gt. nil ) then 
                split = .true. 
                goto 20
              endif
              if ( iOctCh(Neighbors(2)) .gt. nil ) then 
                split = .true. 
                goto 20
              endif
              if ( iOctCh(Neighbors(3)) .gt. nil ) then 
                split = .true. 
                goto 20
              endif
              if ( iOctCh(Neighbors(4)) .gt. nil ) then 
                split = .true. 
                goto 20
              endif
              if ( iOctCh(Neighbors(5)) .gt. nil ) then 
                split = .true. 
                goto 20
              endif
              if ( iOctCh(Neighbors(6)) .gt. nil ) then 
                split = .true. 
              endif
 20           if ( split ) then 
                vnw(1,idcell) = hvar(1,idcell)
                vnw(2,idcell) = hvar(2,idcell)
                vnw(3,idcell) = hvar(3,idcell)
                vnw(4,idcell) = hvar(4,idcell)
                vnw(5,idcell) = hvar(5,idcell)
                vnw(6,idcell) = hvar(8,idcell)
                ref(idcell) = 0.0
#ifdef ADVECT_SPECIES
                do j = 9 , nhvarA
                  vnw(j-2,idcell) = hvar(j,idcell)
                enddo
#endif
              endif
            endif 
          enddo  ! i2
        enddo ! i1          
      ENDIF
c
      return
      end
c
c     --------------------------------
      subroutine H_Old2New_3 ( Level )
c     --------------------------------
c
#     include "a_def.h"
      include 'a_tree.h'
      integer Level
      integer i1, i2, idcell, icell, nLevel, Neighbors(6)
      logical leaf
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j,Neighbors)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves
            call NbAll0 ( i1 , Neighbors )
            if ( iOctCh(Neighbors(1)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(2)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(3)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(4)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(5)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(6)) .gt. nil ) then 
              goto 30 
            endif
            vnw(1,i1) = hvar(1,i1)
            vnw(2,i1) = hvar(2,i1)
            vnw(3,i1) = hvar(3,i1)
            vnw(4,i1) = hvar(4,i1)
            vnw(5,i1) = hvar(5,i1)
            vnw(6,i1) = hvar(8,i1)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              vnw(j-2,i1) = hvar(j,i1)
            enddo
#endif
 30         continue
          endif
        enddo        
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idcell,j,Neighbors)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            if ( iOctCh(idcell) .eq. nil ) then 
              call NbAllL ( idcell , Neighbors )
              if ( iOctCh(Neighbors(1)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(2)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(3)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(4)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(5)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(6)) .gt. nil ) then 
                goto 40 
              endif
              vnw(1,idcell) = hvar(1,idcell)
              vnw(2,idcell) = hvar(2,idcell)
              vnw(3,idcell) = hvar(3,idcell)
              vnw(4,idcell) = hvar(4,idcell)
              vnw(5,idcell) = hvar(5,idcell)
              vnw(6,idcell) = hvar(8,idcell)
#ifdef ADVECT_SPECIES
              do j = 9 , nhvarA
                vnw(j-2,idcell) = hvar(j,idcell)
              enddo
#endif
 40           continue
            endif
          enddo
        enddo           
      ENDIF
c
      return
      end
c
c     ---------------------------------
      subroutine H_Old2New_3b ( Level )
c     ---------------------------------
c
c     zeroes ref
c
      include 'a_tree.h'
      integer Level
      integer i1, i2, idcell, icell, nLevel, Neighbors(6)
      logical leaf
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(Neighbors)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves
            call NbAll0 ( i1 , Neighbors )
            if ( iOctCh(Neighbors(1)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(2)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(3)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(4)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(5)) .gt. nil ) then 
              goto 30
            endif
            if ( iOctCh(Neighbors(6)) .gt. nil ) then 
              goto 30 
            endif
            ref(i1) = 0.
            vnw(6,i1) = hvar(8,i1)
 30         continue
          endif
        enddo        
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idcell,Neighbors)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            if ( iOctCh(idcell) .eq. nil ) then 
              call NbAllL ( idcell , Neighbors )
              if ( iOctCh(Neighbors(1)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(2)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(3)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(4)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(5)) .gt. nil ) then 
                goto 40
              endif
              if ( iOctCh(Neighbors(6)) .gt. nil ) then 
                goto 40 
              endif
              ref(idcell) = 0.0
              vnw(6,idcell) = hvar(8,idcell)
 40           continue
            endif
          enddo
        enddo           
      ENDIF
c
      return
      end
c
c
c     --------------------------------------
      subroutine H_New_to_Old ( Level, Key )
c     --------------------------------------

c     Moves physical variables from new (vnw) to old (var) arrays.
c     Input: Level - level of leaves to move
c            Key   - if .ne. 0 , moves only leaves with no split neighbours

      INCLUDE 'a_numbers.h'
c
      if ( Key .eq. 0 ) then
        call H_New2Old_1 ( Level )
      else
        call H_New2Old_2 ( Level )
      end if
      return
      end
c
c     --------------------------------
      subroutine H_New2Old_1 ( Level )
c     --------------------------------
c
c     new --> old variables (except for internal energy)
c
c
#     include "a_def.h"
      include 'a_tree.h'
      integer Level
      integer i1, i2, idcell, icell, nLevel
      parameter ( rmin = 1.e-30 )
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(j)
        do i1 = 1 , ncell0 
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves
c            hvar(1,i1) = max ( rmin , vnw(1,i1) )
            hvar(1,i1) = vnw(1,i1)
            hvar(2,i1) = vnw(2,i1)
            hvar(3,i1) = vnw(3,i1)
            hvar(4,i1) = vnw(4,i1)
            hvar(5,i1) = vnw(5,i1)
            hvar(8,i1) = vnw(6,i1)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,i1) = vnw(j-2,i1)
            enddo
#endif
          endif
        enddo        
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,j)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
          icell = icell + 1
          if ( iOctCh(icell) .eq. nil ) then   ! only leaves
c            hvar(1,icell) = max ( rmin , vnw(1,icell) )
            hvar(1,icell) = vnw(1,icell)
            hvar(2,icell) = vnw(2,icell)
            hvar(3,icell) = vnw(3,icell)
            hvar(4,icell) = vnw(4,icell)
            hvar(5,icell) = vnw(5,icell)
            hvar(8,icell) = vnw(6,icell)
#ifdef ADVECT_SPECIES
            do j = 9 , nhvarA
              hvar(j,icell) = vnw(j-2,icell)
            enddo
#endif
          endif
        enddo           
      ENDIF
c
      return
      end
c
c     --------------------------------
      subroutine H_New2Old_2 ( Level )
c     --------------------------------
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      integer Level
      integer i1, i2, idcell, icell, nLevel, Neighbors(6)
      logical leaf
      parameter ( rmin = 1.e-30 )
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(leaf,Neighbors,j)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves
            leaf = .true. 
            call NbAll0 ( i1 , Neighbors ) 
            if ( iOctCh(Neighbors(1)) .gt. nil ) then 
              leaf = .false. 
              goto 50
            endif
            if ( iOctCh(Neighbors(2)) .gt. nil ) then 
              leaf = .false. 
              goto 50
            endif
            if ( iOctCh(Neighbors(3)) .gt. nil ) then 
              leaf = .false. 
              goto 50
            endif
            if ( iOctCh(Neighbors(4)) .gt. nil ) then 
              leaf = .false. 
              goto 50
            endif
            if ( iOctCh(Neighbors(5)) .gt. nil ) then 
              leaf = .false. 
              goto 50
            endif
            if ( iOctCh(Neighbors(6)) .gt. nil ) then 
              leaf = .false. 
            endif
 50         if ( leaf ) then               
c              hvar(1,i1) = max ( rmin , vnw(1,i1) )
              hvar(1,i1) = vnw(1,i1)
              hvar(2,i1) = vnw(2,i1)
              hvar(3,i1) = vnw(3,i1)
              hvar(4,i1) = vnw(4,i1)
              hvar(5,i1) = vnw(5,i1)
              hvar(8,i1) = vnw(6,i1)
#ifdef ADVECT_SPECIES
              do j = 9 , nhvarA
                hvar(j,i1) = vnw(j-2,i1)
              enddo
#endif
            endif
          endif
        enddo        
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idcell,leaf,Neighbors,j)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            
            if ( iOctCh(idcell) .eq. nil ) then 
              leaf = .true. 
              call NbAllL ( idcell , Neighbors ) 

              if ( iOctCh(Neighbors(1)) .gt. nil ) then 
                leaf = .false. 
                goto 60
              endif
              if ( iOctCh(Neighbors(2)) .gt. nil ) then 
                leaf = .false. 
                goto 60
              endif
              if ( iOctCh(Neighbors(3)) .gt. nil ) then 
                leaf = .false.
                goto 60
              endif
              if ( iOctCh(Neighbors(4)) .gt. nil ) then 
                leaf = .false. 
                goto 60
              endif
              if ( iOctCh(Neighbors(5)) .gt. nil ) then 
                leaf = .false. 
                goto 60
              endif
              if ( iOctCh(Neighbors(6)) .gt. nil ) then 
                leaf = .false. 
              endif
 60           if ( leaf ) then 
c                hvar(1,idcell) = max ( rmin , vnw(1,idcell) )
                hvar(1,idcell) = vnw(1,idcell)
                hvar(2,idcell) = vnw(2,idcell) 
                hvar(3,idcell) = vnw(3,idcell)
                hvar(4,idcell) = vnw(4,idcell)
                hvar(5,idcell) = vnw(5,idcell)
                hvar(8,idcell) = vnw(6,idcell)
#ifdef ADVECT_SPECIES
                do j = 9 , nhvarA
                  hvar(j,idcell) = vnw(j-2,idcell)
                enddo
#endif
              endif
            endif
          enddo
        enddo           
      ENDIF
c
      return
      end
c
c
