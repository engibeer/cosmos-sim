c     =====================================================================
c                                                                         .
c                ART Version 3: Density assignment routines               .
c                                                                         .
c               by Andrey Kravtsov and Anatoly Klypin (1997)              .
c                                                                         .
c     =====================================================================

c     ---------------------------------
      subroutine Zero_Density ( Level )
c     ---------------------------------
c
c     purpose: zeroes density storage on a given Level
c
      include 'a_tree.h'

      integer Level 

      IF ( Level .eq. MinLevel ) THEN

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1)
        do ic1 = 1 , ncell0
          var(1,ic1) = -1.0
          var(2,ic1) = zero
        enddo

      ELSE

      Size  = (CellSize(Level))     ! cell size on this level   
      Size3 = - Size**3
        call Select_Cells ( Level , nLevel ) 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1 , icell )
       do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          var(2,icell)     = zero
          var(2,icell+1)   =  zero
          var(2,icell+2)   =  zero
          var(2,icell+3)   =  zero
          var(2,icell+4)   =  zero
          var(2,icell+5)   =  zero
          var(2,icell+6)   =  zero
          var(2,icell+7)   =  zero
          var(1,icell)   = Size3
          var(1,icell+1) = Size3
          var(1,icell+2) = Size3
          var(1,icell+3) = Size3
          var(1,icell+4) = Size3
          var(1,icell+5) = Size3
          var(1,icell+6) = Size3
          var(1,icell+7) = Size3
        enddo
      ENDIF

      return
      end

c     ---------------------------------------------------
      subroutine Assign_Density ( MinModify , MaxModify)
c     ---------------------------------------------------
c
c     purpose: assigns density to all existing cells     
c
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_mpi.h'

      PARAMETER (Lbound = 50000 )
      DIMENSION         InBuffer(Nproc),Mbound(0:Lbound)
      integer MinModify , MaxModify
      
      if ((MinModify.lt.MinLevel).or.(MinModify.gt.MaxLevel)) return          
      if ((MaxModify.lt.MinLevel).or.(MaxModify.gt.MaxLevel)) return

      do Level = MaxModify , MinModify , -1
        call Zero_Density ( Level )
      enddo

c      write (13,*) '   --- MinModify/Max=',MinModify , MaxModify

      if ( MinModify .eq. MinLevel ) then
        Level2 = MinModify + 1
        call Assign_Density0 ()
c        If(irank.eq.0)Then
c        do k=1,65
c              write(51,*) k
c           do j=1,65
c              write(51,'(i4,2x,16f8.4)')
c     &          j,(var(1,i+(j-1)*ng+(k-1)*ng**2),i=1,16)
c          enddo
c        enddo
c        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
c        stop
c        endIf
      else
        Level2 = MinModify
      endif

      do Level = MaxLevel , Level2 , -1 
        Call Select_Cells (Level , nLevel )  ! all cells on this level
      If(nLevel.ne.0)Then
c                                     find block of cells for load balancing
       If(MinModify.eq.MaxModify.and.Level.le.6)Then
           LChunk = min(nLevel,MLChunk) ! length of a chunk
        Else
           LChunk = min(nLevel,MLChunk) ! length of a chunk
        Endif 
c        write (13,*) ' LChunk =',LChunk,' nLevel=',nLevel,
c     &       MinModify,Level,MaxModify

        maxbuf =0
        minbuf   =10000000
        Mbound(0) =1
        idummy = 0
        nks   = 0
        Do ic =1,nLevel
           ic0 = iSelect(ic)
           Do i=0,nchild-1
              maxbuf =maxbuf +iNp(ic0+i)
              idummy =idummy +iNp(ic0+i)
           EndDo
           If(idummy.ge.LChunk)Then
              nks = nks +1
              Mbound(nks) = ic
c              write (13,'(" ic=",i9," Np=",i7," nCh=",i4,i6)')
c     &            ic,idummy,nks,LChunk
              idummy = 0
           EndIf
        EndDo
c           If(nks.gt.Lbound)Then
c              write (13,*) ' Error: Buffer for chunks ',
c     &                      ' Mbound is too small=',nks,Lbound
c              write (13,*) '       Increase Lbound'
c              name = 'ERROR_RUN'
c              accesstype = 'append'
c              call Open_ASCII_FILE(19,name,accesstype)
c              rewind 19
c              write (19,*) 
c     &              ' 8  :ART_Density0.f,  Error: Buffer for chunks',
c     &                      ' LBuffer is too small=',
c     &                       nks,Lbound,' Increase Lbound'
c              close(19) 
c              call MPI_Abort(MPI_COMM_WORLD,ierr1,ierr2)
c              Stop
c           EndIf

c                                                           find number of Chunks
        if(idummy.ne.0)nks = nks +1
        Mbound(nks) = nLevel
        nChunks= (nks-1)/Nproc +1
        
c           write (13,200)maxbuf,nChunks,(Mbound(i),i=1,nChunks)
c 200       format(' Npart=',i7,' chunks=',i4,' bounds=',/(10i7))
        maxbuf =0
c        If(MinModify.le.1.and.Level.eq.MaxModify.and.Level.le.6)Then
c           LChunk = min(nLevel,MLChunk/5) ! length of a chunk
c        Else
c           LChunk = min(nLevel,MLChunk) ! length of a chunk
c        Endif 
c        nChunks= (nLevel-1)/LChunk +1        ! number of chunks        
c        Lrecord= (nLevel)/(Nproc*nChunks)  
c        Lastrec= nLevel -(nChunks*Nproc-1)*Lrecord ! last may be shorter
c        If(Lastrec.le.0)Then
c           write (13,*) ' ERR: nLevel=',nLevel,' Chunks=',nChunks,
c     .                 ' Length=',LChunk,Lrecord
c           STOP
c        EndIf
c        write (13,300)iChunk,nChunks,nlevel,Level
c 300    format(' chunk=',i3,' nChunks=',i2,
c     &                 ' nlevel=',i9,' Level=',i3)

        Do iChunk =1,nChunks    ! Loop over all chunks
           LProc = Nproc
           If(iChunk.eq.nChunks)LProc =nks-(nChunks-1)*Nproc
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED), !SCHEDULE(DYNAMIC) 
C$OMP+PRIVATE ( iproc, ic, jCellEnd,jCellSt,Mend,icount,ic0,icell)
          Do iproc =1,LProc
             icount = 0
             jCellSt = (iChunk-1)*Nproc+(iproc-1)
             jCellEnd =jCellSt +1
             Mend =Mbound(jCellEnd)-1
          if(iChunk.eq.nChunks.and.iproc.eq.LProc)Mend=Mbound(jCellEnd)
            do ic0 = Mbound(jCellSt),Mend
             icell = iSelect(ic0) +nchild
c            write (13,'("   ip=",i6," Bounds=",2i7," chunk=",2i8
c     &        )')iproc,Mbound(iChunk-1),Mbound(iChunk)-1,
c     &        iChunk,ic0

             call Dens_Child(icell, 
     .              Level,Level2,MaxModify,iproc,icount)
            enddo                  ! end ic
            InBuffer(iproc) =icount
          Enddo         ! end iproc

          Do iproc =1,LProc
             maxbuf =max(maxbuf,InBuffer(iproc))
             minbuf =min(minbuf,InBuffer(iproc))
             Do iB =1,InBuffer(iproc)
                iAcc =iBuffer(iB,iproc)
                var(1,iAcc) = var(1,iAcc) +wBuffer(iB,iproc)
                var(2,iAcc)   = var(2,iAcc) + rBuffer(iB,iproc)
             enddo   
          enddo 
        enddo                      ! end iChunk
c        idummy = 0
c        Do ic =1,nLevel
c           ic0 = iSelect(ic)
c           Do i=0,nchild-1
c              idummy =idummy +iNp(ic0+i)
c           EndDo
c        EndDo
c           write (13,400) Level,nLevel,nChunks,
c     .         LChunk,Lproc,minbuf,maxbuf,
c     .           Float(maxbuf)/MAX(10,minbuf)
c 400       format(' Lev=',i3,' nLevel=',i9,' Chunks=',i4,
c     .         ' Length=',i6,i6,' Buffer=',i8,i9,2x,f7.1)

        EndIf
      enddo                     ! end Level

      return
      end
c     ---------------------------------------------------
      subroutine Dens_Child (icell,Level,Level2,MaxModify,iproc,
     &               icount)
c     ---------------------------------------------------
c
c     purpose: find density contributions by all particles in
c              cell 'icell'; store them in (i/w)Buffer
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'
c      DIMENSION iBuffer(LBuffer,Nproc),
c     &                  wBuffer(LBuffer,Nproc),
c     &                  rBuffer(LBuffer,Nproc)
c      COMMON /DensBuffa/ iBuffer(LBuffer,Nproc)
c      COMMON /DensBuffb/ wBuffer(LBuffer,Nproc)
c      COMMON /DensBuffc/ rBuffer(LBuffer,Nproc)
      common /sumparticles/ nnpart     !!!
      dimension iNC(27)            ! array for neighbors and corners
      double precision diff_x,diff_y,diff_z,
     &                          corr_x,corr_y,corr_z
      Logical OutBuffer
      OutBuffer = .false.
 
      IF(mpi_yes .NE. 1234567 ) THEN
         lsp1 = lsp(1)
      ELSE
         lsp1 = n_refin 
      ENDIF

        do ic1 = 1 , nchild 
          idcell = icell - ic1 
          iHead  = iCL(idcell)     ! first particle in the cell
c          write (13,*) '    icell=',icell,idcell,ic1,iHead

          if ( (iOctCh(idcell)  .eq.      nil)      ! work only with leaves 
     &                         .and.                ! which are non-empty
     &         (iHead           .ne.      nil) ) then 
            idummy1 = idcell
            do ic2 = Level , Level2 , -1 ! assign on leaf level and above
              if ( ic2 .le. MaxModify ) then
              Size = CellSize(ic2)
              call Find_Neighbors_and_Corners ( idummy1 ,
     &                                              ic2 , 
     &                                             Size , 
     &                                              iNC , 
     &                                              nnc   )
              call Ps ( idummy1 , x0,y0,z0)
              do ic3 = 1 , nnc
                iAcc  = iNC(ic3)             
                call Ps ( iAcc , Posx,Posy,Posz )
                diff_x = x0 - Posx
                diff_y = y0 - Posy
                diff_z = z0 - Posz
       
                corr_x = zero
                corr_y = zero
                corr_z = zero        
                if ( abs(diff_x) .gt. nf67 ) then   ! boundary effects
                  if ( diff_x .gt. zero ) then      ! nf67 = ng / 2
                    corr_x = -ng
                  else
                    corr_x = ng
                  endif
                endif 
                if ( abs(diff_y) .gt. nf67 ) then              
                  if ( diff_y .gt. zero ) then
                    corr_y = -ng
                  else
                    corr_y = ng
                  endif
                endif 
                if ( abs(diff_z) .gt. nf67 ) then              
                  if ( diff_z .gt. zero ) then
                    corr_z = -ng
                  else
                    corr_z = ng
                  endif
                endif 
                idummy2 = iHead               ! loop over particles
                do while ( idummy2 .ne. nil )
                  t1 = 1.0 - abs(x(idummy2) + corr_x - Posx) / Size 
                  t2 = 1.0 - abs(y(idummy2) + corr_y - Posy) / Size
                  t3 = 1.0 - abs(z(idummy2) + corr_z - Posz) / Size
                  if ((t1.gt.0.).and.(t2.gt.0.).and.(t3.gt.0.)) then
c                    ispecie     = iWhichSpecie (idummy2)
c                    t3t1t2w     = t3 * t1 * t2 * wpar(ispecie)
                    t3t1t2w     = t3 * t1 * t2 * Wpart(idummy2)
                    icount =icount +1
                  if(icount.gt.LBuffer)OutBuffer =.true.
c                  if(OutBuffer.or.icount.gt.LBuffer-100)Then
c                    write (13,*) ' Error: Buffer for density',
c     &                  ' LBuffer is too small=',icount,' cell=',icell
c                    write (13,*) '   Particle=',idummy2,
c     &               x(idummy2),y(idummy2),z(idummy2)
c                    write (13,*) '   Proc=',iproc,' Weight=',
c     &               Wpart(idummy2),' First Specie=',lsp1
c                    EndIf
c                    If(OutBuffer) Stop

                    iBuffer(icount,iproc) = iAcc
                    wBuffer(icount,iproc) = t3t1t2w
c                    if ( ispecie .eq. 1 ) then 
                    IF(mpi_yes .NE. 1234567 ) THEN
                      if ( idummy2 .le. lsp1 ) then  
                        rBuffer(icount,iproc) = t3t1t2w
                      else
                        rBuffer(icount,iproc) = 0.
                      endif
                    ELSE
                      if ( idummy2 .le. lsp1 .and.
     &                     wpart( idummy2 ) .le. wpar_min ) then  
                        rBuffer(icount,iproc) = t3t1t2w
                      else
                        rBuffer(icount,iproc) = 0.
                      endif
                    ENDIF
                 endif               
                  idummy2 = iLL(idummy2,1)
                enddo  ! end do while
       
              enddo            !end ic3
              endif
              idummy1 = iPr(idummy1)
            enddo  ! end ic2
          endif
       enddo                  ! end ic1
             if(OutBuffer)Then
                 write (13,*) ' Error: Buffer for density',
     &                      ' LBuffer is too small=',icount,icell
                 write(name,'(a,i4.4,a)')'ERROR_RUN',irank+1,'.DAT'
                 accesstype = 'append'
                 call Open_ASCII_FILE(19,name,accesstype)
                 rewind 19
                 write (19,*) 
     &              ' 7  :ART_Density0.f,  Error: Buffer for density',
     &                      ' LBuffer is too small=',
     &                       icount,iproc
                 close(19)
                 close(13)
                 CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
                 Stop
             EndIf
      return
      end
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

c     -----------------------------------------------------------
      subroutine Find_Neighbors_and_Corners_o ( icell , Level , 
     &                                        Size  , iNC   , nnc )
c     -----------------------------------------------------------
c
c     purpose: finds neighbors and adjacent corners for icell - 
c              only those of the icell level are taken
c   
c     this version is for "serial" use - it is more compact and easier
c     to read `
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
      dimension i(3)
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

      do ic1 = 1 , 3      

        iDir1 = iENb(icn,ic1) 
        iDir2 = iENb(icn,iDir(ic1,2))
        iDir3 = iENb(icn,iDir(ic1,3))
        iNb1  = iNb(iPar,iDir1)             ! current external neighbor of iPar
        iNb2  = iNb(iNb1,iDir2)
        iNb3  = iNb(iNb2,iDir3)

        if ( iOctCh(iNb1) .gt. nil ) then ! if neighbor is split
          do ic2 = 1 , 4
            nnc      = nnc + 1
            iNC(nnc) = iCh(iNb1,iChDir(iDir1,ic2))
          enddo
        endif

c....   iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * (k + 1) + 1

        if ( iOctCh(iNb2) .gt. nil ) then ! if neighbors neighbor is split
          do ic2 = 1 , 3
            i(ic2) = - iNbDir(iDir1,ic2) - iNbDir(iDir2,ic2)  
          enddo          
          idm      = iDir(ic1,3)
          i(idm)   = i(idm) + 1
          nCh      = (i(1) + 1)/2 + i(2) + 2 * i(3) + 4
          nnc      = nnc + 1
          iNC(nnc) = iCh(iNb2,nCh)
          i(idm)   = i(idm) - 2
          nCh      = (i(1) + 1)/2 + i(2) + 2 * i(3) + 4
          nnc      = nnc + 1
          iNC(nnc) = iCh(iNb2,nCh)
        endif

        if ( iFlag .eq. nil ) then
          if ( (iOctCh(iNb3)  .gt.   nil) 
     &                        .and. 
     &         (iLv(iNb3) + 1 .eq. Level) ) then
              iFlag = 1
              do ic2 = 1 , 3
                i(ic2) = - iNbDir(iDir1,ic2) 
     &                   - iNbDir(iDir2,ic2)  
     &                   - iNbDir(iDir3,ic2)  
              enddo
              nCh = (i(1) + 1)/2 + i(2) + 2 * i(3) + 4
              nnc      = nnc + 1
              iNC(nnc) = iCh(iNb3,nCh)
          endif
        endif

      enddo

      return
      end

c     -----------------------------
      subroutine Assign_Density0 ()
c     -----------------------------
c
c     assigns density on the root mesh
c     algorithm: cloud-in-cell 
c
      include 'a_tree.h'
      include 'a_control.h'
      double precision xx , yy , zz
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))

c.... Loop over particles

      IF(mpi_yes .NE. 1234567 ) THEN
         lsp1 = lsp(1)
      ELSE
         lsp1 = n_refin 
      ENDIF
      idtot = 0 
      IF(mpi_yes .NE. 1234567 ) n_all = lspecies(Nspecies)  
      do in = 1 , n_all

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

c        ispecie = iWhichSpecie ( in )
c        t2w = t2 * wpar(ispecie)
c        d2w = d2 * wpar(ispecie)
        t2w = t2 *  Wpart(in)
        d2w = d2 *  Wpart(in)

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
c        if ( ispecie .eq. 1 ) then 
        if ( in .le. lsp1 .and. wpart( in ) .le. wpar_min ) then 
          var(2,ind1) = var(2,ind1) + t3 * t1 * t2w
          var(2,ind2) = var(2,ind2) + t3 * d1 * t2w
          var(2,ind3) = var(2,ind3) + t3 * t1 * d2w
          var(2,ind4) = var(2,ind4) + t3 * d1 * d2w
          var(2,ind5) = var(2,ind5) + d3 * t1 * t2w
          var(2,ind6) = var(2,ind6) + d3 * d1 * t2w
          var(2,ind7) = var(2,ind7) + d3 * t1 * d2w
          var(2,ind8) = var(2,ind8) + d3 * d1 * d2w
        endif

       
      enddo

        return
      end




