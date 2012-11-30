c     =====================================================================
c                                                                         .
c                ART Version 3: Density assignment routines               .
c                                                                         .
c               by Andrey Kravtsov and Anatoly Klypin (1997)              .
c                                                                         .
c     =====================================================================

c     ------------------------------------------------------
      subroutine Assign_Density ( Level1 , Level2 , NOFlag )
c     ------------------------------------------------------
c
c     purpose: assigns density to cells on levels Level1-Level2
c      Level1 <= Level2
c
c     NOFlag = 0 : assign old density in hvar(1,*)
c     NOFlag = 1 : assign new density in vnw(1,*)
c
      integer Level, Level1, Level2, NOFlag
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      do Level = Level1, Level2
        call Initialize_Density ( Level ) 
      enddo
c
#ifdef HYDRO
      do Level = Level1 , Level2
        call Assign_GasDensity ( Level , NOFlag )
      enddo
#endif
c
#ifdef NBODY
      call Assign_ParticleDensity ( Level1 , Level2 )
#endif
c
      return
      end
c
c     ---------------------------------------
      subroutine Initialize_Density ( Level ) 
c     ---------------------------------------
c
      integer Level 
      double precision Size3
      include 'a_tree.h'
      include 'a_control.h'

      IF ( Level .eq. MinLevel ) THEN
c
C$OMP PARALLEL DO DEFAULT(SHARED)
        do ic1 = 1 , ncell0
          var(1,ic1) = -1.0
          ref(ic1) = 0.
        enddo
      ELSE
c
        Size3 = -CellSize(Level)**3
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          var(1,icell)   = Size3 
          var(1,icell+1) = Size3
          var(1,icell+2) = Size3
          var(1,icell+3) = Size3
          var(1,icell+4) = Size3
          var(1,icell+5) = Size3
          var(1,icell+6) = Size3
          var(1,icell+7) = Size3
          ref(icell)   = 0.
          ref(icell+1) = 0.
          ref(icell+2) = 0.
          ref(icell+3) = 0.
          ref(icell+4) = 0.
          ref(icell+5) = 0.
          ref(icell+6) = 0.
          ref(icell+7) = 0.
        enddo
      ENDIF
      
      return
      end
c
c     -----------------------------------------------
      subroutine Assign_GasDensity ( Level , NOFlag )
c     -----------------------------------------------
c
c     purpose: initializes density on a given Level with gas density
c
      integer Level, NOFlag  
      double precision Size3
      include 'a_tree.h'
      include 'a_control.h'
c
      IF ( Level .eq. MinLevel ) THEN
c
        IF ( NOFlag .eq. 0 ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
          do ic1 = 1 , ncell0
            var(1,ic1) = var(1,ic1) + hvar(1,ic1)
          enddo
        ELSE
C$OMP PARALLEL DO DEFAULT(SHARED)
          do ic1 = 1 , ncell0
            var(1,ic1) = var(1,ic1) + vnw(1,ic1)
          enddo          
        ENDIF
c
      ELSE
c
        Size3 = CellSize(Level)**3
        call Select_Cells ( Level , nLevel ) 
        IF ( NOFlag .eq. 0 ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell)
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            var(1,icell)   = var(1,icell)   + hvar(1,icell)   * Size3
            var(1,icell+1) = var(1,icell+1) + hvar(1,icell+1) * Size3
            var(1,icell+2) = var(1,icell+2) + hvar(1,icell+2) * Size3
            var(1,icell+3) = var(1,icell+3) + hvar(1,icell+3) * Size3 
            var(1,icell+4) = var(1,icell+4) + hvar(1,icell+4) * Size3 
            var(1,icell+5) = var(1,icell+5) + hvar(1,icell+5) * Size3
            var(1,icell+6) = var(1,icell+6) + hvar(1,icell+6) * Size3
            var(1,icell+7) = var(1,icell+7) + hvar(1,icell+7) * Size3 
          enddo
        ELSE
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell)
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            var(1,icell)   = var(1,icell)   + vnw(1,icell)   * Size3
            var(1,icell+1) = var(1,icell+1) + vnw(1,icell+1) * Size3
            var(1,icell+2) = var(1,icell+2) + vnw(1,icell+2) * Size3
            var(1,icell+3) = var(1,icell+3) + vnw(1,icell+3) * Size3 
            var(1,icell+4) = var(1,icell+4) + vnw(1,icell+4) * Size3 
            var(1,icell+5) = var(1,icell+5) + vnw(1,icell+5) * Size3
            var(1,icell+6) = var(1,icell+6) + vnw(1,icell+6) * Size3
            var(1,icell+7) = var(1,icell+7) + vnw(1,icell+7) * Size3 
          enddo
        ENDIF
      ENDIF
c
      return
      end
c
c     -----------------------------------------------------
      subroutine Assign_ParticleDensity ( Level1 , Level2 )
c     -----------------------------------------------------
c
c     purpose: assigns density to all existing cells     
c
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      COMMON /DensBuff/ iBuffer(LBuffer,Nproc),wBuffer(LBuffer,Nproc), 
     &                  rBuffer(LBuffer,Nproc)
      DIMENSION         InBuffer(Nproc)
      integer Level, Level1, Level2
      
      if ((Level1.lt.MinLevel).or.(Level1.gt.MaxLevel)) return          
      if ((Level2.lt.MinLevel).or.(Level2.gt.MaxLevel)) return          
      

      if ( Level1 .eq. MinLevel ) then
        Ldummy = Level1 + 1
        call Assign_ParticleDensity0 ()
      else
        Ldummy = Level1
      endif

      

      maxbuf =0
      do L = MaxLevel , Ldummy , -1 
        call Select_Cells ( L , nLevel )  ! all cells on this level

      If(nLevel.ne.0)Then
        LChunk = min(nLevel,MLChunk) ! length of a chunk
        nChunks= (nLevel-1)/LChunk +1        ! number of chunks        
        Lrecord= (nLevel)/(Nproc*nChunks)  
        Lastrec= nLevel -(nChunks*Nproc-1)*Lrecord ! last may be shorter
        If(Lastrec.le.0)Then
           write (*,*) ' ERR: nLevel=',nLevel,' Chunks=',nChunks,
     .                 ' Length=',LChunk,Lrecord,Lastrec 
           STOP
        EndIf

        Do iChunk = 1 , nChunks    ! Loop over all chunks

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(ic,ioffset,icount,ic0,icell,Lrec0)
          Do iproc = 1 , Nproc
            ioffset= ((iproc-1)+(iChunk-1)*Nproc)*Lrecord
            icount =0
            Lrec0  =Lrecord

            if(iproc.eq.Nproc.and.iChunk.eq.nChunks)Lrec0=Lastrec
            do ic0 = ioffset+1 , ioffset+Lrec0
             icell = iSelect(ic0) 
             call Dens_Child(icell, 
     &              L,Ldummy,Level2,iproc,icount)
#ifdef DEBUG
        call Open_ASCII_File ( iOWork , workfile, append ) 
        write(iOWork,*) iproc, ic0, icell, L, Ldummy, Level2,icount
        close(iOWork)
#endif
            enddo    ! end ic0
            InBuffer(iproc) =icount
          Enddo         ! end iproc

          do iproc = 1 , Nproc
            icount = InBuffer(iproc)
            if ( icount .gt. LBuffer ) then
              call Open_ASCII_File ( iOErr , errorfile , sequent )
              write (iOErr,*) '1  * Error in Dens_Child:',icount,iproc
              write (iOErr,*) 'stopping...'
              close(iOErr)
              stop
            endif
          enddo

          Do iproc =1,Nproc
             maxbuf =max(maxbuf,InBuffer(iproc))
             Do iB =1,InBuffer(iproc)
                iAcc =iBuffer(iB,iproc)
                var(1,iAcc) = var(1,iAcc) +wBuffer(iB,iproc)
                ref(iAcc)   = ref(iAcc) + rBuffer(iB,iproc)
             enddo   
          enddo 
        enddo                      ! end iChunk
        EndIf
      enddo                     ! end L
#ifdef DEBUG
      write (*,11)InBuffer,maxbuf,L,nChunks     
 11   format(9i7,3i3)
      write (*,*) ' ic0=',ic0,' must be =',nLevel
#endif
      return
      end
c
c     ---------------------------------------------------
      subroutine Dens_Child (icell,Level,Level2,MaxModify,iproc,icount)
c     ---------------------------------------------------
c
c     purpose: find density contributions by all particles in
c              cell 'icell'; store them in (i/w)Buffer
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      COMMON /DensBuff/ iBuffer(LBuffer,Nproc),wBuffer(LBuffer,Nproc),
     &                  rBuffer(LBuffer,Nproc)
      dimension iNC(27)            ! array for neighbors and corners

        do ic1 = 0 , nchild-1
          idcell = icell + ic1        
          if ( iOctCh(idcell)  .eq.      nil) then     ! work only with leaves 
            iHead  = iCL(idcell)     ! first particle in the cell
          if (iHead           .ne.      nil ) then    ! which are non-empty
            idummy1 = idcell
            
            do ic2 = Level , Level2 , -1 ! assign on leaf level and above
              if ( ic2 .le. MaxModify ) then ! only on Level1-Level2
              Size = CellSize(ic2)
              rSize = 1.0 / Size
              call Find_Neighbors_and_Corners ( idummy1 ,
     &                                              ic2 , 
     &                                             Size , 
     &                                              iNC , 
     &                                              nnc   )
              call Ps ( idummy1 , x0,y0,z0)
c              write (*,*) ' ic2=',ic2,nnc,' x=',x0,y0,z0
              do ic3 = 1 , nnc
                iAcc  = iNC(ic3)             
                if ( iAcc .ne. idummy1 ) then 
                call Ps ( iAcc , Posx,Posy,Posz )
                diff_x = x0 - Posx
                diff_y = y0 - Posy
                diff_z = z0 - Posz
       
                if ( abs(diff_x) .gt. nf67 ) then   ! boundary effects
                  if ( diff_x .gt. zero ) then      ! nf67 = ng / 2
                    Posx = Posx + ng
                  else
                    Posx = Posx - ng
                  endif
                endif 
                if ( abs(diff_y) .gt. nf67 ) then              
                  if ( diff_y .gt. zero ) then
                    Posy = Posy + ng
                  else
                    Posy = Posy - ng
                  endif
                endif 
                if ( abs(diff_z) .gt. nf67 ) then              
                  if ( diff_z .gt. zero ) then
                    Posz = Posz + ng
                  else
                    Posz = Posz - ng
                  endif
                endif 
                else
                  Posx = x0
                  Posy = y0
                  Posz = z0
                endif
                idummy2 = iHead               ! loop over particles
                do while ( idummy2 .ne. nil )

                  t1 = abs(x(idummy2) - Posx) * rSize 
                  if ( t1 .lt. 1. ) then  
                  t2 = abs(y(idummy2) - Posy) * rSize
                  if ( t2 .lt. 1. ) then 
                  t3 = abs(z(idummy2) - Posz) * rSize
                  if ( t3 .lt. 1. ) then 
                    icount  = icount + 1

                    iBuffer(icount,iproc) = iAcc
                    wBuffer(icount,iproc) = (1.0-t3)*(1.0-t1)*(1.0-t2)* 
     &                                      pw(idummy2)
#ifdef STARFORM 
                    if ( idummy2 .le. nsp(1,2) .or.
     &                   idummy2 .ge. nsp(nspecies,1) ) then 
#endif
#ifndef STARFORM 
                    if ( idummy2 .le. nsp(1,2) ) then 
#endif
                    rBuffer(icount,iproc) = (1.0-t3)*(1.0-t1)*(1.0-t2)*
     &                                      pw(idummy2)
                    else
                       rBuffer(icount,iproc) = 0.
                    endif
                  endif               
                  endif
                  endif
                  idummy2 = iLL(idummy2,1)
                enddo  ! end do while
       
              enddo            !end ic3
              endif
              idummy1 = iPr(idummy1)
            enddo  ! end ic2
          endif
          endif
       enddo                  ! end ic1
      return
      end
c
c     -----------------------------------------------------------
      subroutine Find_Neighbors_and_Corners ( icell , Level , 
     &                                        Size  , iNC   , nnc )
c     -----------------------------------------------------------
c
c     purpose: finds neighbors and adjacent corners for icell - 
c              only those of the icell level are taken
c   
c     this version is for "parallel" use - some computations are manually
c     unrolled to optimize the computations and memory accesses
c
c     input  : icell - cell for which routine finds neighbors and corners
c              Level - iLv(icell); MUST be > MinLevel
c              Size  - CellSize(level)
c
c     output : iNC   - indices of neighbors and corners
c              nnc   - number of neighbors and corners (size of iNC)     
c     
c
      include 'a_tree.h'
c
      dimension iNC(27)                 ! array for neighbors and corners
      dimension iENb(8,3)               ! external neighbors for given child
      dimension iChDir(6,4)             ! which children of ext.neighb to take
      dimension iNbDir(6,3)             ! vector describing a given direction
      dimension iDir(3,3)               ! specifies the order of Directions 
      integer   i1 , i2 , i3 
      integer   icell, Level , nnc
      real      Size
c
      data iENb    /  1, 2, 1, 2, 1, 2, 1, 2,       ! "x"
     &                3, 3, 4, 4, 3, 3, 4, 4,       ! "y"
     &                5, 5, 5, 5, 6, 6, 6, 6  /     ! "z"
c
      data iChDir  /  2, 1, 3, 1, 5, 1,             
     &                4, 3, 4, 2, 6, 2,
     &                6, 5, 7, 5, 7, 3,
     &                8, 7, 8, 6, 8, 4  /
c      
      data iNbDir  / -1, 1, 0, 0, 0, 0, 
     &                0, 0,-1, 1, 0, 0,
     &                0, 0, 0, 0,-1, 1  /
c
      data iDir    /  1, 2, 3,
     &                2, 3, 1,
     &                3, 1, 2  /
c
      icn   = mod ( icell - ncell0 - 1 , nchild ) + 1 ! get child number
      iPar  = iPr(icell)           ! get parent
      iCh1  = iCh(iPar,1)
      iFlag = nil                  ! zero flag
      nnc   = nil                  ! zero counter 

      iNC(1) = iCh1
      iNC(2) = iCh1 + 1
      iNC(3) = iCh1 + 2
      iNC(4) = iCh1 + 3
      iNC(5) = iCh1 + 4
      iNC(6) = iCh1 + 5
      iNC(7) = iCh1 + 6
      iNC(8) = iCh1 + 7
            
      nnc = nnc + 8

      iENb1 = iENb(icn,1)
      iENb2 = iENb(icn,2)
      iENb3 = iENb(icn,3)
    
c....   first direction

        iDir1 = iENb1 
        iDir2 = iENb2
        iDir3 = iENb3
        iNb1  = iNb(iPar,iDir1)             ! current external neighbor of iPar
        iNb2  = iNb(iNb1,iDir2)
        iNb3  = iNb(iNb2,iDir3)
        iOC = iOctCh(iNb1)
        if ( iOC .gt. nil ) then ! if neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0
          do ic2 = 1 , 4
            nnc      = nnc + 1
            iNC(nnc) = iChsh +iChDir(iDir1,ic2)
          enddo
        endif

c....   iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * (k + 1) + 1
        iOC = iOctCh(iNb2)
        if ( iOC .gt. nil ) then ! if neighbors neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0
          i1       = - iNbDir(iDir1,1) - iNbDir(iDir2,1) + 1
          i2       = - iNbDir(iDir1,2) - iNbDir(iDir2,2) + 4
          i3       = - iNbDir(iDir1,3) - iNbDir(iDir2,3) + 1
          nCh      = i1/2 + i2 + 2 * i3 
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh -4
        endif

        if ( iFlag .eq. nil ) then
           iOC =iOctCh(iNb3)
          if ( (iOC  .gt.   nil) 
     &                        .and. 
     &         (iLv(iNb3) + 1 .eq. Level) ) then
                iFlag = 1
                i1 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,1)  
     &               - iNbDir(iDir3,1) 
                i2 = - iNbDir(iDir1,2) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,2)
                i3 = - iNbDir(iDir1,3) 
     &               - iNbDir(iDir2,3)  
     &               - iNbDir(iDir3,3)  
              nCh = (i1 + 1)/2 + i2 + 2 * i3 + 4
              nnc      = nnc + 1
              iNC(nnc) = ( iOC - 1 ) * nchild + nCh + ncell0
          endif
        endif

c....   second direction    

        iDir1 = iENb2 
        iDir2 = iENb3
        iDir3 = iENb1
        iNb1  = iNb(iPar,iDir1)             ! current external neighbor of iPar
        iNb2  = iNb(iNb1,iDir2)
        iNb3  = iNb(iNb2,iDir3)
        iOC   = iOctCh(iNb1)
        if ( iOC .gt. nil ) then ! if neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0
          do ic2 = 1 , 4
            nnc      = nnc + 1
            iNC(nnc) = iChsh +iChDir(iDir1,ic2)
          enddo
        endif

c....   iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * (k + 1) + 1
        iOC = iOctCh(iNb2)
        if ( iOC .gt. nil ) then ! if neighbors neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0 
          i1 = - iNbDir(iDir1,1) - iNbDir(iDir2,1) + 1
          i2 = - iNbDir(iDir1,2) - iNbDir(iDir2,2) + 4
          i3 = - iNbDir(iDir1,3) - iNbDir(iDir2,3)  
          nCh      = (i1 + 1)/2 + i2 + 2 * i3 
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh -1
        endif

        if ( iFlag .eq. nil ) then
          iOC =iOctCh(iNb3) 
          if ( (iOC  .gt.   nil) 
     &                        .and. 
     &         (iLv(iNb3) + 1 .eq. Level) ) then            
                iFlag = 1
                i1 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,3)  
                i2 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,3)  
                i3 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,3)  
              nCh = (i1 + 1)/2 + i2 + 2 * i3 + 4
              nnc      = nnc + 1
              iNC(nnc) = ( iOC - 1 ) * nchild + nCh + ncell0
          endif
        endif

c....   third direction
     
        iDir1 = iENb3 
        iDir2 = iENb1
        iDir3 = iENb2
        iNb1  = iNb(iPar,iDir1)             ! current external neighbor of iPar
        iNb2  = iNb(iNb1,iDir2)
        iNb3  = iNb(iNb2,iDir3)
        iOC   = iOctCh(iNb1)
        if ( iOC .gt. nil ) then ! if neighbor is split
           iChsh = ( iOC - 1 ) * nchild  + ncell0
           do ic2 = 1 , 4
            nnc      = nnc + 1
            iNC(nnc) = iChsh +iChDir(iDir1,ic2)
          enddo
        endif

c....   iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * (k + 1) + 1
        iOC = iOctCh(iNb2)
        if ( iOC .gt. nil ) then ! if neighbors neighbor is split
          iChsh = ( iOC - 1 ) * nchild  + ncell0           
          i1       = - iNbDir(iDir1,1) - iNbDir(iDir2,1) + 1
          i2       = - iNbDir(iDir1,2) - iNbDir(iDir2,2) + 5 
          i3       = - iNbDir(iDir1,3) - iNbDir(iDir2,3) 
          nCh      = i1/2 + i2 + 2 * i3
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh
          nnc      = nnc + 1
          iNC(nnc) = iChsh + nCh -2
        endif

        if ( iFlag .eq. nil ) then
          iOC =iOctCh(iNb3) 
          if ( (iOC  .gt.   nil) 
     &                        .and. 
     &         (iLv(iNb3) + 1 .eq. Level) ) then
              iFlag = 1
                i1 = - iNbDir(iDir1,1) 
     &               - iNbDir(iDir2,1)  
     &               - iNbDir(iDir3,1)  
                i2 = - iNbDir(iDir1,2) 
     &               - iNbDir(iDir2,2)  
     &               - iNbDir(iDir3,2)  
                i3 = - iNbDir(iDir1,3) 
     &               - iNbDir(iDir2,3)  
     &               - iNbDir(iDir3,3)  
              nCh = (i1 + 1)/2 + i2 + 2 * i3 + 4
              nnc      = nnc + 1
              iNC(nnc) = ( iOC - 1 ) * nchild + nCh + ncell0
          endif
        endif

      return
      end
c
c     -------------------------------------
      subroutine Assign_ParticleDensity0 ()
c     -------------------------------------
c
c     assigns density on the root mesh
c     algorithm: cloud-in-cell 
c
      include 'a_tree.h'

c.... Loop over particles

      do in = 1 , nsp(1,2) ! 1st specie first to avoid checking for ref

	xx = x(in)
        yy = y(in)
        zz = z(in)
        i  = int(xx - 0.5)      ! 0.5 - to make my mesh consistent with CIC
        j  = int(yy - 0.5)
	k  = int(zz - 0.5)

        i1  = i + 1
        j1  = j + 1
        k1  = k + 1

        d1  = (xx - i) - 0.5
        d2  = (yy - j) - 0.5
        d3  = (zz - k) - 0.5
        t1  = 1.0 - d1
        t2  = 1.0 - d2
        t3  = 1.0 - d3

        t2w = t2 * pw(in)
        d2w = d2 * pw(in)

        if (i  .lt.  1) i  = ng
        if (j  .lt.  1) j  = ng
        if (k  .lt.  1) k  = ng
        if (i1 .gt. ng) i1 = 1
        if (j1 .gt. ng) j1 = 1
        if (k1 .gt. ng) k1 = 1

c....   get indices of cell vertices

        ind1 = iCellIndex (  i ,  j ,  k ) 
        ind2 = iCellIndex ( i1 ,  j ,  k )  
        ind3 = iCellIndex (  i , j1 ,  k )  
        ind4 = iCellIndex ( i1 , j1 ,  k )  
        ind5 = iCellIndex (  i ,  j , k1 )  
        ind6 = iCellIndex ( i1 ,  j , k1 )  
        ind7 = iCellIndex (  i , j1 , k1 )  
        ind8 = iCellIndex ( i1 , j1 , k1 )   

        var(1,ind1) = var(1,ind1) + t3 * t1 * t2w
        var(1,ind2) = var(1,ind2) + t3 * d1 * t2w
        var(1,ind3) = var(1,ind3) + t3 * t1 * d2w
        var(1,ind4) = var(1,ind4) + t3 * d1 * d2w
        var(1,ind5) = var(1,ind5) + d3 * t1 * t2w
        var(1,ind6) = var(1,ind6) + d3 * d1 * t2w
        var(1,ind7) = var(1,ind7) + d3 * t1 * d2w
        var(1,ind8) = var(1,ind8) + d3 * d1 * d2w

        ref(ind1) = ref(ind1) + t3 * t1 * t2w
        ref(ind2) = ref(ind2) + t3 * d1 * t2w
        ref(ind3) = ref(ind3) + t3 * t1 * d2w
        ref(ind4) = ref(ind4) + t3 * d1 * d2w
        ref(ind5) = ref(ind5) + d3 * t1 * t2w
        ref(ind6) = ref(ind6) + d3 * d1 * t2w
        ref(ind7) = ref(ind7) + d3 * t1 * d2w
        ref(ind8) = ref(ind8) + d3 * d1 * d2w
      enddo

      do in = nsp(1,2)+1 , np

	xx = x(in)
        yy = y(in)
        zz = z(in)
        i  = int(xx - 0.5)      ! 0.5 - to make my mesh consistent with CIC
        j  = int(yy - 0.5)
	k  = int(zz - 0.5)

        i1  = i + 1
        j1  = j + 1
        k1  = k + 1

        d1  = (xx - i) - 0.5
        d2  = (yy - j) - 0.5
        d3  = (zz - k) - 0.5
        t1  = 1.0 - d1
        t2  = 1.0 - d2
        t3  = 1.0 - d3

        t2w = t2 * pw(in)
        d2w = d2 * pw(in)

        if (i  .lt.  1) i  = ng
        if (j  .lt.  1) j  = ng
        if (k  .lt.  1) k  = ng
        if (i1 .gt. ng) i1 = 1
        if (j1 .gt. ng) j1 = 1
        if (k1 .gt. ng) k1 = 1

c....   get indices of cell vertices

        ind1 = iCellIndex (  i ,  j ,  k ) 
        ind2 = iCellIndex ( i1 ,  j ,  k )  
        ind3 = iCellIndex (  i , j1 ,  k )  
        ind4 = iCellIndex ( i1 , j1 ,  k )  
        ind5 = iCellIndex (  i ,  j , k1 )  
        ind6 = iCellIndex ( i1 ,  j , k1 )  
        ind7 = iCellIndex (  i , j1 , k1 )  
        ind8 = iCellIndex ( i1 , j1 , k1 )   

        var(1,ind1) = var(1,ind1) + t3 * t1 * t2w
        var(1,ind2) = var(1,ind2) + t3 * d1 * t2w
        var(1,ind3) = var(1,ind3) + t3 * t1 * d2w
        var(1,ind4) = var(1,ind4) + t3 * d1 * d2w
        var(1,ind5) = var(1,ind5) + d3 * t1 * t2w
        var(1,ind6) = var(1,ind6) + d3 * d1 * t2w
        var(1,ind7) = var(1,ind7) + d3 * t1 * d2w
        var(1,ind8) = var(1,ind8) + d3 * d1 * d2w
      enddo

      return
      end




