c     ------------------------------------
      subroutine Sweep22 ( iDir, dt, Level )
c     ------------------------------------
c     
c     Input     :    iDir            - direction of sweep ( 1 - 6 )
c                    dt              - time step
c                    Level           - level of sweep in the tree
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 v(nhvar,-1:2) ! input for Interface Solvers 
      real*8 g(2)                 
      real*8 c(2) 
      real*8 index(2)
      real*8 ff(5), f(5,iSwChunk), fdum(5)
      integer icell(2, iSwChunk)           
c
      integer iP, iR1, iR2, iL1, iL2
      integer iLvR1, iLvR2, iLvL1
      integer iKyR1, j3, j4, j5
c
      integer jpermute(7,6)
      data jpermute / 1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7 /
c
      dimension aspect(2,iSwChunk), asp(3)
      data asp / 0.5, 0.25, 0.125 /
      integer iopposite(6)
      data iopposite / 2, 1, 4, 3, 6, 5 /
      real*8 c112, c122
      parameter ( c112 = 1. / 1.25, c122 = 1. / 1.5 )
c
c      common / Sweep_Com1 / v
c      common / Sweep_Com2 / g
c      common / Sweep_Com3 / c
c      common / Sweep_Com4 / f
c      common / Sweep_Com5 / icell
c      common / Sweep_Com6 / nif
c
      dimension lnbr(0:8,6), kiot(0:8,6), nbnb(8,6)
      data lnbr / 0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 5, 6, 7, 8, 1, 2, 3, 4,      ! idir = 5,6
     &            0, 5, 6, 7, 8, 1, 2, 3, 4/      ! idir = 5,6
      data kiot / 1, 0, 1, 0, 1, 0, 1, 0, 1,      ! idir = 1
     &            1, 1, 0, 1, 0, 1, 0, 1, 0,      ! idir = 2
     &            1, 0, 0, 1, 1, 0, 0, 1, 1,      ! idir = 3
     &            1, 1, 1, 0, 0, 1, 1, 0, 0,      ! idir = 4
     &            1, 0, 0, 0, 0, 1, 1, 1, 1,      ! idir = 5
     &            1, 1, 1, 1, 1, 0, 0, 0, 0/      ! idir = 6
      data nbnb / 0, 1, 0, 3, 0, 5, 0, 7, 
     &            2, 0, 4, 0, 6, 0, 8, 0,
     &            0, 0, 1, 2, 0, 0, 5, 6, 
     &            3, 4, 0, 0, 7, 8, 0, 0,
     &            0, 0, 0, 0, 1, 2, 3, 4, 
     &            5, 6, 7, 8, 0, 0, 0, 0 /  
      integer iopp(6)
      data iopp / 2, 1, 4, 3, 6, 5 /


c      call Timing(1,0)
c
      dx = CellSize(Level)   
      dtx = dt / dx 
      dtx2 = 0.5*dtx
c              
      j3 = jpermute(3,iDir)
      j4 = jpermute(4,iDir)
      j5 = jpermute(5,iDir)
      iodir = iopp(idir)
c
      nif = 0
c
      IF ( Level .eq. MinLevel ) THEN
C.... OpenMP
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iP,iopp,iR1,iR2,iL1,iL2,iLvL1,iLvR1,iLvR2)
C$OMP+PRIVATE(v,icell,index,c,g,ff,f)
C$OMP+PRIVATE(ic0,ic1,i2)
C$OMP+FIRSTPRIVATE(nif)
        do i1 = 1 , ncell0
          iP = i1
          if ( iOctCh(iP) .eq. nil ) then   ! only leaves
c             call iSwNb0 ( iP , idir , iR1 , iR2 , iL1  )
             iLvL1 = 0
             iLvR1 = 0 
             iLvR2 = 0 
            iR1   = iNb ( iP  , iDir )
            iR2   = iNb ( iR1 , iDir )
            iL1   = iNb ( iP  , iodir )
            iL2   = iNb ( iL1 , iodir )
            iLvL1 = iLv ( iL1 )
            iLvR1 = iLv ( iR1 )
            iLvR2 = iLv ( iR2 )
            if ( iOctCh(iR1) .eq. 0 ) then ! here we need to consider only =|= case
              nif = nif + 1
              v(1,-1) = hvar(1,iL1)
              v(1, 0) = hvar(1,iP )
              v(1, 1) = hvar(1,iR1)
              v(1, 2) = hvar(1,iR2)
              v(2,-1) = hvar(2,iL1)
              v(2, 0) = hvar(2,iP )
              v(2, 1) = hvar(2,iR1)
              v(2, 2) = hvar(2,iR2)
              v(3,-1) = hvar(j3,iL1)
              v(3, 0) = hvar(j3,iP )
              v(3, 1) = hvar(j3,iR1)
              v(3, 2) = hvar(j3,iR2)
              v(4,-1) = hvar(j4,iL1)
              v(4, 0) = hvar(j4,iP )
              v(4, 1) = hvar(j4,iR1)
              v(4, 2) = hvar(j4,iR2)
              v(5,-1) = hvar(j5,iL1)
              v(5, 0) = hvar(j5,iP )
              v(5, 1) = hvar(j5,iR1)
              v(5, 2) = hvar(j5,iR2)
              v(6,-1) = hvar(6,iL1)
              v(6, 0) = hvar(6,iP )
              v(6, 1) = hvar(6,iR1)
              v(6, 2) = hvar(6,iR2)
              v(7,-1) = hvar(7,iL1)
              v(7, 0) = hvar(7,iP )
              v(7, 1) = hvar(7,iR1)
              v(7, 2) = hvar(7,iR2)              
c
              icell(1, nif) = iP
              icell(2, nif) = iR1
              index(1)      = iP
              index(2)      = iR1
              c(1)     = 1.0
c              if ( iLvL1 .lt. Level ) c(1) = c112
              c(2)     = 1.0
c              if ( iLvR2 .lt. Level ) c(2) = c112
              call FU ( dt, dx, dtx, dtx2, index, g, v, c, ff )
              if ( lapidus ) then 
                call FD ( iDir, dt, dx, index, g, v, ff )
              endif
              f(1,nif) = ff(1)
              f(2,nif) = ff(2)
              f(3,nif) = ff(3)
              f(4,nif) = ff(4)
              f(5,nif) = ff(5)
            endif
          endif
c
c           ! the loop below must be guarded in parallel execution
          if ( (nif .eq. iSwChunk).or.(i1 .eq. ncell0) ) then 
            do i2 = 1 , nif 
              ic0 = icell(1,i2)
c$omp critical 
              vnw(1 ,ic0) = vnw(1 ,ic0) - f(1,i2)    
              vnw(2 ,ic0) = vnw(2 ,ic0) - f(5,i2)    
              vnw(j3,ic0) = vnw(j3,ic0) - f(2,i2)    
              vnw(j4,ic0) = vnw(j4,ic0) - f(3,i2)    
              vnw(j5,ic0) = vnw(j5,ic0) - f(4,i2)    
c$omp end critical
              ic1 = icell(2,i2)
C$OMP CRITICAL            
              vnw(1 ,ic1) = vnw(1 ,ic1) + f(1,i2)
              vnw(2 ,ic1) = vnw(2 ,ic1) + f(5,i2)
              vnw(j3,ic1) = vnw(j3,ic1) + f(2,i2)
              vnw(j4,ic1) = vnw(j4,ic1) + f(3,i2)
              vnw(j5,ic1) = vnw(j5,ic1) + f(4,i2)
C$OMP END CRITICAL
            enddo
            nif = 0 
          endif       
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C.... OpenMP
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(ic,idcell)
C$OMP+PRIVATE(iP,iopp,iR1,iR2,iL1,iL2,iLvL1,iLvR1,iLvR2)
C$OMP+PRIVATE(iKyR1,v,icell,index,c,aspect,g,ff,f)
C$OMP+PRIVATE(ic0,ic1,i2,i3,fdum)
C$OMP+FIRSTPRIVATE(nif)
        do i1 = 1 , nLevel
          ic = iSelect(i1)
          do i2 = 0 , 7
            idcell = ic + i2
            iP = idcell
            if ( iOctCh(iP) .eq. nil ) then   ! only leaves
              iR1   = iNb ( iP  , iDir )
              iR2   = iNb ( iR1 , iDir )
              iL1   = iNb ( iP  , iodir )
              iL2   = iNb ( iL1 , iodir )
              iLvL1 = iLv ( iL1 )
              iLvR1 = iLv ( iR1 )
              iLvR2 = iLv ( iR2 )
c              call iSwNbL ( iP , iDir , iR1 , iR2 , iL1 , iL2 , 
c     &                      iLvR1 , iLvR2 , iLvL1 )

              if ( iOctCh(iR1) .eq. 0 ) then  
                iKyR1 = 0   ! leaf
              else
                iKyR1 = 1   ! split
              endif            
              if ( iKyR1 .eq. 0 ) then 
                if ( iLvR1 .eq. Level ) then  ! = | = case
                  nif = nif + 1
                  v(1,-1) = hvar(1,iL1)
                  v(1, 0) = hvar(1,iP )
                  v(1, 1) = hvar(1,iR1)
                  v(1, 2) = hvar(1,iR2)
                  v(2,-1) = hvar(2,iL1)
                  v(2, 0) = hvar(2,iP )
                  v(2, 1) = hvar(2,iR1)
                  v(2, 2) = hvar(2,iR2)
                  v(3,-1) = hvar(j3,iL1)
                  v(3, 0) = hvar(j3,iP )
                  v(3, 1) = hvar(j3,iR1)
                  v(3, 2) = hvar(j3,iR2)
                  v(4,-1) = hvar(j4,iL1)
                  v(4, 0) = hvar(j4,iP )
                  v(4, 1) = hvar(j4,iR1)
                  v(4, 2) = hvar(j4,iR2)
                  v(5,-1) = hvar(j5,iL1)
                  v(5, 0) = hvar(j5,iP )
                  v(5, 1) = hvar(j5,iR1)
                  v(5, 2) = hvar(j5,iR2)
                  v(6,-1) = hvar(6,iL1)
                  v(6, 0) = hvar(6,iP )
                  v(6, 1) = hvar(6,iR1)
                  v(6, 2) = hvar(6,iR2)
                  v(7,-1) = hvar(7,iL1)
                  v(7, 0) = hvar(7,iP )
                  v(7, 1) = hvar(7,iR1)
                  v(7, 2) = hvar(7,iR2)              
c
                  icell(1, nif) = iP
                  icell(2, nif) = iR1
                  c(1)     = 1.0
                  if ( iLvL1 .lt. Level ) c(1) = c112
                  c(2)     = 1.0
                  if ( iLvR2 .lt. Level ) c(2) = c112
                  aspect(1, nif) = 1.
                  aspect(2, nif) = 1.
                  call FU ( dt, dx, dtx, dtx2 , index , g, v, c, ff)
                  if ( lapidus ) then 
                    call FD( iDir, dt, dx, index, g, v, ff)
                  endif
                  f(1,nif) = ff(1)
                  f(2,nif) = ff(2)
                  f(3,nif) = ff(3)
                  f(4,nif) = ff(4)
                  f(5,nif) = ff(5)
                else
                  nif = nif + 1
                  v(1,-1) = hvar(1,iL1)
                  v(1, 0) = hvar(1,iP )
                  v(1, 1) = hvar(1,iR1)
                  v(1, 2) = hvar(1,iR2)
                  v(2,-1) = hvar(2,iL1)
                  v(2, 0) = hvar(2,iP )
                  v(2, 1) = hvar(2,iR1)
                  v(2, 2) = hvar(2,iR2)
                  v(3,-1) = hvar(j3,iL1)
                  v(3, 0) = hvar(j3,iP )
                  v(3, 1) = hvar(j3,iR1)
                  v(3, 2) = hvar(j3,iR2)
                  v(4,-1) = hvar(j4,iL1)
                  v(4, 0) = hvar(j4,iP )
                  v(4, 1) = hvar(j4,iR1)
                  v(4, 2) = hvar(j4,iR2)
                  v(5,-1) = hvar(j5,iL1)
                  v(5, 0) = hvar(j5,iP )
                  v(5, 1) = hvar(j5,iR1)
                  v(5, 2) = hvar(j5,iR2)
                  v(6,-1) = hvar(6,iL1)
                  v(6, 0) = hvar(6,iP )
                  v(6, 1) = hvar(6,iR1)
                  v(6, 2) = hvar(6,iR2)
                  v(7,-1) = hvar(7,iL1)
                  v(7, 0) = hvar(7,iP )
                  v(7, 1) = hvar(7,iR1)
                  v(7, 2) = hvar(7,iR2)              
                  icell(1,nif)  = iP
                  icell(2,nif)  = iR1
                  index(1) = iP
                  index(2) = iR1
                  aspect(1,nif) = 1.
                  aspect(2,nif) = asp(ndim)
                  c(1)      = c112
                  c(2)      = c122
                  call FU ( dt, dx, dtx, dtx2 , index , g, v, c, ff)
                  if ( lapidus ) then 
                    call FD( iDir, dt, dx, index, g, v, ff)
                  endif
                  f(1,nif) = ff(1)
                  f(2,nif) = ff(2)
                  f(3,nif) = ff(3)
                  f(4,nif) = ff(4)
                  f(5,nif) = ff(5)
                endif
              endif
              if ( nif .eq. iSwChunk ) then 
c               ! the loop below must be guarded in parallel execution
                do i3 = 1 , nif 
                  ic0 = icell(1,i3)
                  fdum(1) = aspect(1,i3) * f(1,i3)
                  fdum(2) = aspect(1,i3) * f(2,i3)
                  fdum(3) = aspect(1,i3) * f(3,i3)
                  fdum(4) = aspect(1,i3) * f(4,i3)
                  fdum(5) = aspect(1,i3) * f(5,i3)
C$OMP CRITICAL
                  vnw(1 ,ic0) = vnw(1 ,ic0) - fdum(1)    
                  vnw(2 ,ic0) = vnw(2 ,ic0) - fdum(5)    
                  vnw(j3,ic0) = vnw(j3,ic0) - fdum(2)    
                  vnw(j4,ic0) = vnw(j4,ic0) - fdum(3)    
                  vnw(j5,ic0) = vnw(j5,ic0) - fdum(4)    
C$OMP END CRITICAL
                  ic1 = icell(2,i3)
                  fdum(1) = aspect(2,i3) * f(1,i3)
                  fdum(2) = aspect(2,i3) * f(2,i3)
                  fdum(3) = aspect(2,i3) * f(3,i3)
                  fdum(4) = aspect(2,i3) * f(4,i3)
                  fdum(5) = aspect(2,i3) * f(5,i3)
C$OMP CRITICAL
                  vnw(1 ,ic1) = vnw(1 ,ic1) + fdum(1)
                  vnw(2 ,ic1) = vnw(2 ,ic1) + fdum(5)
                  vnw(j3,ic1) = vnw(j3,ic1) + fdum(2)
                  vnw(j4,ic1) = vnw(j4,ic1) + fdum(3)
                  vnw(j5,ic1) = vnw(j5,ic1) + fdum(4)
C$OMP END CRITICAL
                enddo
                nif = 0 
              endif       
              if ( iLvL1 .eq. Level-1 ) then 
                nif = nif + 1
                v(1,-1) = hvar(1,iL2)
                v(1, 0) = hvar(1,iL1)
                v(1, 1) = hvar(1,iP )
                v(1, 2) = hvar(1,iR1)
                v(2,-1) = hvar(2,iL2)
                v(2, 0) = hvar(2,iL1)
                v(2, 1) = hvar(2,iP )
                v(2, 2) = hvar(2,iR1)
                v(3,-1) = hvar(j3,iL2)
                v(3, 0) = hvar(j3,iL1)
                v(3, 1) = hvar(j3,iP )
                v(3, 2) = hvar(j3,iR1)
                v(4,-1) = hvar(j4,iL2)
                v(4, 0) = hvar(j4,iL1)
                v(4, 1) = hvar(j4,iP )
                v(4, 2) = hvar(j4,iR1)
                v(5,-1) = hvar(j5,iL2)
                v(5, 0) = hvar(j5,iL1)
                v(5, 1) = hvar(j5,iP )
                v(5, 2) = hvar(j5,iR1)
                v(6,-1) = hvar(6,iL2)
                v(6, 0) = hvar(6,iL1)
                v(6, 1) = hvar(6,iP )
                v(6, 2) = hvar(6,iR1)
                v(7,-1) = hvar(7,iL2)
                v(7, 0) = hvar(7,iL1)
                v(7, 1) = hvar(7,iP )
                v(7, 2) = hvar(7,iR1)
                icell (1, nif) = iL1
                icell (2, nif) = iP
                index(1) = iL1
                index(2) = iP
                aspect(1, nif) = asp(ndim)
                aspect(2, nif) = 1.
                c(1) = c122
                c(2) = c112
                call FU ( dt, dx, dtx , dtx2 , index , g, v, c, ff )
                if ( lapidus ) then 
                  call FD ( iDir, dt, dx, index, g, v, ff )
                endif
                f(1,nif) = ff(1)
                f(2,nif) = ff(2)
                f(3,nif) = ff(3)
                f(4,nif) = ff(4)
                f(5,nif) = ff(5)
              endif
            endif

            if (         (nif .eq. iSwChunk)
     &                        .or. 
     &           ((i1.eq. nLevel) .and. (i2.eq.7)) ) then 
c             ! the loop below must be guarded in parallel execution
                do i3 = 1 , nif 
                  ic0 = icell(1,i3)
                  fdum(1) = aspect(1,i3) * f(1,i3)
                  fdum(2) = aspect(1,i3) * f(2,i3)
                  fdum(3) = aspect(1,i3) * f(3,i3)
                  fdum(4) = aspect(1,i3) * f(4,i3)
                  fdum(5) = aspect(1,i3) * f(5,i3)
C$OMP CRITICAL
                  vnw(1 ,ic0) = vnw(1 ,ic0) - fdum(1)    
                  vnw(2 ,ic0) = vnw(2 ,ic0) - fdum(5)    
                  vnw(j3,ic0) = vnw(j3,ic0) - fdum(2)    
                  vnw(j4,ic0) = vnw(j4,ic0) - fdum(3)    
                  vnw(j5,ic0) = vnw(j5,ic0) - fdum(4)    
C$OMP END CRITICAL
                  ic1 = icell(2,i3)
                  fdum(1) = aspect(2,i3) * f(1,i3)
                  fdum(2) = aspect(2,i3) * f(2,i3)
                  fdum(3) = aspect(2,i3) * f(3,i3)
                  fdum(4) = aspect(2,i3) * f(4,i3)
                  fdum(5) = aspect(2,i3) * f(5,i3)
C$OMP CRITICAL
                  vnw(1 ,ic1) = vnw(1 ,ic1) + fdum(1)
                  vnw(2 ,ic1) = vnw(2 ,ic1) + fdum(5)
                  vnw(j3,ic1) = vnw(j3,ic1) + fdum(2)
                  vnw(j4,ic1) = vnw(j4,ic1) + fdum(3)
                  vnw(j5,ic1) = vnw(j5,ic1) + fdum(4)
C$OMP END CRITICAL
                enddo
              nif = 0 
            endif       
          enddo
        enddo
      ENDIF
c
c      call Timing ( 1 , 0 )
c
      return
      end

c     ------------------------------------
      subroutine Sweep_oldest ( iDir, dt, Level )
c     ------------------------------------
c     
c     Input     :    iDir            - direction of sweep ( 1 - 6 )
c                    dt              - time step
c                    Level           - level of sweep in the tree
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 v(nhvar,-1:2) ! input for Interface Solvers 
      real*8 g(2)                 
      real*8 c(2) 
      real*8 index(2)
      real*8 ff(5), f(5,iSwChunk), fdum(5)
      integer icell(2, iSwChunk)           
c
      integer iP, iR1, iR2, iL1, iL2
      integer iLvR1, iLvR2, iLvL1
      integer iKyR1, iopp, j3, j4, j5
c
      integer jpermute(7,6)
      data jpermute / 1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7 /
c
      dimension aspect(2,iSwChunk), asp(3)
      data asp / 0.5, 0.25, 0.125 /
      integer iopposite(6)
      data iopposite / 2, 1, 4, 3, 6, 5 /
      real*8 c112, c122
      parameter ( c112 = 1. / 1.25, c122 = 1. / 1.5 )
c
c      common / Sweep_Com1 / v
c      common / Sweep_Com2 / g
c      common / Sweep_Com3 / c
c      common / Sweep_Com4 / f
c      common / Sweep_Com5 / icell
c      common / Sweep_Com6 / nif
c
      call Timing(1,0)
c
      dx = CellSize(Level)   
      dtx = dt / dx 
      dtx2 = 0.5*dtx
c              
      j3 = jpermute(3,iDir)
      j4 = jpermute(4,iDir)
      j5 = jpermute(5,iDir)
c
      nif = 0
c
      IF ( Level .eq. MinLevel ) THEN
C.... OpenMP
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP PRIVATE(iP,iopp,iR1,iR2,iL1,iL2,iLvL1,iLvR1,iLvR2)
ccC$OMP PRIVATE(iKyR1,v,icell,index,c,g,ff,f)
ccC$OMP PRIVATE(ic0,ic1,i2)
ccC$OMP FIRSTPRIVATE(nif)
        do i1 = 1 , ncell0
          iP = i1
          if ( iOctCh(iP) .eq. nil ) then   ! only leaves
            iopp  = iopposite(iDir)
            iR1   = iNb ( iP  , iDir )
            iR2   = iNb ( iR1 , iDir )
            iL1   = iNb ( iP  , iopp )
            iL2   = iNb ( iL1 , iopp )
            iLvL1 = iLv ( iL1 )
            iLvR1 = iLv ( iR1 )
            iLvR2 = iLv ( iR2 )
c            call iSwNb0 ( iP, iDir, iR1_2, iR2_2, iL1_2 )
c            if ( iR1 .ne. iR1_2 ) then
c            write(*,*) 'iR1 =',iR1,' iR1_2 =', iR1_2, 'iDir =',iDir
c            endif 
c            if ( iR2 .ne. iR2_2 ) then
c              write(*,*) 'iR2 =',iR2,' iR2_2 =', iR2_2, 'iDir =',iDir
c            endif 
c            if ( iL1 .ne. iL1_2 ) then
c              write(*,*) 'iL1 =',iL1,' iL1_2 =', iL1_2, 'iOpp =',iopp
c            endif 
            if ( iOctCh(iR1) .eq. 0 ) then  
              iKyR1 = 0   ! leaf
            else
              iKyR1 = 1   ! split
            endif            
            if ( iKyR1 .eq. 0 ) then ! here we need to consider only =|= case
              nif = nif + 1
              v(1,-1) = hvar(1,iL1)
              v(1, 0) = hvar(1,iP )
              v(1, 1) = hvar(1,iR1)
              v(1, 2) = hvar(1,iR2)
              v(2,-1) = hvar(2,iL1)
              v(2, 0) = hvar(2,iP )
              v(2, 1) = hvar(2,iR1)
              v(2, 2) = hvar(2,iR2)
              v(3,-1) = hvar(j3,iL1)
              v(3, 0) = hvar(j3,iP )
              v(3, 1) = hvar(j3,iR1)
              v(3, 2) = hvar(j3,iR2)
              v(4,-1) = hvar(j4,iL1)
              v(4, 0) = hvar(j4,iP )
              v(4, 1) = hvar(j4,iR1)
              v(4, 2) = hvar(j4,iR2)
              v(5,-1) = hvar(j5,iL1)
              v(5, 0) = hvar(j5,iP )
              v(5, 1) = hvar(j5,iR1)
              v(5, 2) = hvar(j5,iR2)
              v(6,-1) = hvar(6,iL1)
              v(6, 0) = hvar(6,iP )
              v(6, 1) = hvar(6,iR1)
              v(6, 2) = hvar(6,iR2)
              v(7,-1) = hvar(7,iL1)
              v(7, 0) = hvar(7,iP )
              v(7, 1) = hvar(7,iR1)
              v(7, 2) = hvar(7,iR2)              
c
              icell(1, nif) = iP
              icell(2, nif) = iR1
              index(1)      = iP
              index(2)      = iR1
              c(1)     = 1.0
              if ( iLvL1 .lt. Level ) c(1) = c112
              c(2)     = 1.0
              if ( iLvR2 .lt. Level ) c(2) = c112
              call FU ( dt, dx, dtx , dtx2 , index, g, v, c, ff )
              if ( lapidus ) then 
                call FD ( iDir, dt, dx, index, g, v, ff )
              endif
              f(1,nif) = ff(1)
              f(2,nif) = ff(2)
              f(3,nif) = ff(3)
              f(4,nif) = ff(4)
              f(5,nif) = ff(5)
            endif
          endif
c
c           ! the loop below must be guarded in parallel execution
          if ( (nif .eq. iSwChunk).or.(i1 .eq. ncell0) ) then 
            do i2 = 1 , nif 
              ic0 = icell(1,i2)
cc$omp critical 
              vnw(1 ,ic0) = vnw(1 ,ic0) - f(1,i2)    
              vnw(2 ,ic0) = vnw(2 ,ic0) - f(5,i2)    
              vnw(j3,ic0) = vnw(j3,ic0) - f(2,i2)    
              vnw(j4,ic0) = vnw(j4,ic0) - f(3,i2)    
              vnw(j5,ic0) = vnw(j5,ic0) - f(4,i2)    
cc$omp end critical
              ic1 = icell(2,i2)
cC$OMP CRITICAL            
              vnw(1 ,ic1) = vnw(1 ,ic1) + f(1,i2)
              vnw(2 ,ic1) = vnw(2 ,ic1) + f(5,i2)
              vnw(j3,ic1) = vnw(j3,ic1) + f(2,i2)
              vnw(j4,ic1) = vnw(j4,ic1) + f(3,i2)
              vnw(j5,ic1) = vnw(j5,ic1) + f(4,i2)
cC$OMP END CRITICAL
            enddo
            nif = 0 
          endif       
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C.... OpenMP
cC$OMP PARALLEL DO DEFAULT(SHARED)
cC$OMP+PRIVATE(ic,i2,idcell)
cC$OMP+PRIVATE(iP,iopp,iR1,iR2,iL1,iL2,iLvL1,iLvR1,iLvR2)
cC$OMP+PRIVATE(iKyR1,nif,v,icell,index,c,g,ff,f)
cC$OMP+PRIVATE(ic0,ic1,i2,i3,fdum)
cC$OMP+FIRSTPRIVATE(nif)
        do i1 = 1 , nLevel
          ic = iSelect(i1)
          do i2 = 0 , 7
            idcell = ic + i2
c            if ( idcell .eq. 3079209) then 
c              idummiest = 1
c            endif
            iP = idcell
            if ( iOctCh(iP) .eq. nil ) then   ! only leaves
              iopp  = iopposite(iDir)
              iR1   = iNb ( iP  , iDir )
              iR2   = iNb ( iR1 , iDir )
              iL1   = iNb ( iP  , iopp )
              iL2   = iNb ( iL1 , iopp )
              iLvL1 = iLv ( iL1 )
              iLvR1 = iLv ( iR1 )
              iLvR2 = iLv ( iR2 )
c            call iSwNbL ( iP, iDir, iR1_2, iR2_2, iL1_2, iL2_2 ,
c     &                    iLvR1_2, iLvR2_2, iLvL1_2)
c            if ( iLvL1 .ne. iLvL1_2 ) then
c           write(*,*) iP,' iL1 =',iL1,' iLvL1,_2 =',iLvL1,iLvL1_2
c            endif 
c            if ( iLvR1 .ne. iLvR1_2 ) then
c           write(*,*) iP,' iR1 =',iR1,' iLvR1,_2 =',iLvR1,iLvR1_2
c            endif 
c            if ( iLvR2 .ne. iLvR2_2 ) then
c           write(*,*) iP,' iR2 =',iR2,' iLvR2,_2 =',iLvR2,iLvR2_2
c            endif 

              if ( iOctCh(iR1) .eq. 0 ) then  
                iKyR1 = 0   ! leaf
              else
                iKyR1 = 1   ! split
              endif            
              if ( iKyR1 .eq. 0 ) then 
                if ( iLvR1 .eq. Level ) then  ! = | = case
                  nif = nif + 1
                  v(1,-1) = hvar(1,iL1)
                  v(1, 0) = hvar(1,iP )
                  v(1, 1) = hvar(1,iR1)
                  v(1, 2) = hvar(1,iR2)
                  v(2,-1) = hvar(2,iL1)
                  v(2, 0) = hvar(2,iP )
                  v(2, 1) = hvar(2,iR1)
                  v(2, 2) = hvar(2,iR2)
                  v(3,-1) = hvar(j3,iL1)
                  v(3, 0) = hvar(j3,iP )
                  v(3, 1) = hvar(j3,iR1)
                  v(3, 2) = hvar(j3,iR2)
                  v(4,-1) = hvar(j4,iL1)
                  v(4, 0) = hvar(j4,iP )
                  v(4, 1) = hvar(j4,iR1)
                  v(4, 2) = hvar(j4,iR2)
                  v(5,-1) = hvar(j5,iL1)
                  v(5, 0) = hvar(j5,iP )
                  v(5, 1) = hvar(j5,iR1)
                  v(5, 2) = hvar(j5,iR2)
                  v(6,-1) = hvar(6,iL1)
                  v(6, 0) = hvar(6,iP )
                  v(6, 1) = hvar(6,iR1)
                  v(6, 2) = hvar(6,iR2)
                  v(7,-1) = hvar(7,iL1)
                  v(7, 0) = hvar(7,iP )
                  v(7, 1) = hvar(7,iR1)
                  v(7, 2) = hvar(7,iR2)              
c
                  icell(1, nif) = iP
                  icell(2, nif) = iR1
                  c(1)     = 1.0
                  if ( iLvL1 .lt. Level ) c(1) = c112
                  c(2)     = 1.0
                  if ( iLvR2 .lt. Level ) c(2) = c112
                  aspect(1, nif) = 1.
                  aspect(2, nif) = 1.
                  call FU ( dt, dx , dtx, dtx2,index , g, v, c, ff)
                  if ( lapidus ) then 
                    call FD( iDir, dt, dx, index, g, v, ff)
                  endif
                  f(1,nif) = ff(1)
                  f(2,nif) = ff(2)
                  f(3,nif) = ff(3)
                  f(4,nif) = ff(4)
                  f(5,nif) = ff(5)
                else
                  nif = nif + 1
                  v(1,-1) = hvar(1,iL1)
                  v(1, 0) = hvar(1,iP )
                  v(1, 1) = hvar(1,iR1)
                  v(1, 2) = hvar(1,iR2)
                  v(2,-1) = hvar(2,iL1)
                  v(2, 0) = hvar(2,iP )
                  v(2, 1) = hvar(2,iR1)
                  v(2, 2) = hvar(2,iR2)
                  v(3,-1) = hvar(j3,iL1)
                  v(3, 0) = hvar(j3,iP )
                  v(3, 1) = hvar(j3,iR1)
                  v(3, 2) = hvar(j3,iR2)
                  v(4,-1) = hvar(j4,iL1)
                  v(4, 0) = hvar(j4,iP )
                  v(4, 1) = hvar(j4,iR1)
                  v(4, 2) = hvar(j4,iR2)
                  v(5,-1) = hvar(j5,iL1)
                  v(5, 0) = hvar(j5,iP )
                  v(5, 1) = hvar(j5,iR1)
                  v(5, 2) = hvar(j5,iR2)
                  v(6,-1) = hvar(6,iL1)
                  v(6, 0) = hvar(6,iP )
                  v(6, 1) = hvar(6,iR1)
                  v(6, 2) = hvar(6,iR2)
                  v(7,-1) = hvar(7,iL1)
                  v(7, 0) = hvar(7,iP )
                  v(7, 1) = hvar(7,iR1)
                  v(7, 2) = hvar(7,iR2)              
                  icell(1,nif)  = iP
                  icell(2,nif)  = iR1
                  index(1) = iP
                  index(2) = iR1
                  aspect(1,nif) = 1.
                  aspect(2,nif) = asp(ndim)
                  c(1)      = c112
                  c(2)      = c122
                  call FU ( dt, dx, dtx, dtx2, index , g, v, c, ff)
                  if ( lapidus ) then 
                    call FD( iDir, dt, dx, index, g, v, ff)
                  endif
                  f(1,nif) = ff(1)
                  f(2,nif) = ff(2)
                  f(3,nif) = ff(3)
                  f(4,nif) = ff(4)
                  f(5,nif) = ff(5)
                endif
              endif
              if ( nif .eq. iSwChunk ) then 
c               ! the loop below must be guarded in parallel execution
                do i3 = 1 , nif 
                  ic0 = icell(1,i3)
                  fdum(1) = aspect(1,i3) * f(1,i3)
                  fdum(2) = aspect(1,i3) * f(2,i3)
                  fdum(3) = aspect(1,i3) * f(3,i3)
                  fdum(4) = aspect(1,i3) * f(4,i3)
                  fdum(5) = aspect(1,i3) * f(5,i3)
cC$OMP CRITICAL
                  vnw(1 ,ic0) = vnw(1 ,ic0) - fdum(1)    
                  vnw(2 ,ic0) = vnw(2 ,ic0) - fdum(5)    
                  vnw(j3,ic0) = vnw(j3,ic0) - fdum(2)    
                  vnw(j4,ic0) = vnw(j4,ic0) - fdum(3)    
                  vnw(j5,ic0) = vnw(j5,ic0) - fdum(4)    
cC$OMP END CRITICAL
                  ic1 = icell(2,i3)
                  fdum(1) = aspect(2,i3) * f(1,i3)
                  fdum(2) = aspect(2,i3) * f(2,i3)
                  fdum(3) = aspect(2,i3) * f(3,i3)
                  fdum(4) = aspect(2,i3) * f(4,i3)
                  fdum(5) = aspect(2,i3) * f(5,i3)
cC$OMP CRITICAL
                  vnw(1 ,ic1) = vnw(1 ,ic1) + fdum(1)
                  vnw(2 ,ic1) = vnw(2 ,ic1) + fdum(5)
                  vnw(j3,ic1) = vnw(j3,ic1) + fdum(2)
                  vnw(j4,ic1) = vnw(j4,ic1) + fdum(3)
                  vnw(j5,ic1) = vnw(j5,ic1) + fdum(4)
cC$OMP END CRITICAL
                enddo
                nif = 0 
              endif       
              if ( iLvL1 .eq. Level-1 ) then 
                nif = nif + 1
                v(1,-1) = hvar(1,iL2)
                v(1, 0) = hvar(1,iL1)
                v(1, 1) = hvar(1,iP )
                v(1, 2) = hvar(1,iR1)
                v(2,-1) = hvar(2,iL2)
                v(2, 0) = hvar(2,iL1)
                v(2, 1) = hvar(2,iP )
                v(2, 2) = hvar(2,iR1)
                v(3,-1) = hvar(j3,iL2)
                v(3, 0) = hvar(j3,iL1)
                v(3, 1) = hvar(j3,iP )
                v(3, 2) = hvar(j3,iR1)
                v(4,-1) = hvar(j4,iL2)
                v(4, 0) = hvar(j4,iL1)
                v(4, 1) = hvar(j4,iP )
                v(4, 2) = hvar(j4,iR1)
                v(5,-1) = hvar(j5,iL2)
                v(5, 0) = hvar(j5,iL1)
                v(5, 1) = hvar(j5,iP )
                v(5, 2) = hvar(j5,iR1)
                v(6,-1) = hvar(6,iL2)
                v(6, 0) = hvar(6,iL1)
                v(6, 1) = hvar(6,iP )
                v(6, 2) = hvar(6,iR1)
                v(7,-1) = hvar(7,iL2)
                v(7, 0) = hvar(7,iL1)
                v(7, 1) = hvar(7,iP )
                v(7, 2) = hvar(7,iR1)
                icell (1, nif) = iL1
                icell (2, nif) = iP
                index(1) = iL1
                index(2) = iP
                aspect(1, nif) = asp(ndim)
                aspect(2, nif) = 1.
                c(1) = c122
                c(2) = c112
                call FU ( dt, dx, dtx, dtx2, index , g, v, c, ff )
                if ( lapidus ) then 
                  call FD ( iDir, dt, dx, index, g, v, ff )
                endif
                f(1,nif) = ff(1)
                f(2,nif) = ff(2)
                f(3,nif) = ff(3)
                f(4,nif) = ff(4)
                f(5,nif) = ff(5)
              endif
            endif

            if (         (nif .eq. iSwChunk)
     &                        .or. 
     &           ((i1.eq. nLevel) .and. (i2.eq.7)) ) then 
c             ! the loop below must be guarded in parallel execution
                do i3 = 1 , nif 
                  ic0 = icell(1,i3)
                  fdum(1) = aspect(1,i3) * f(1,i3)
                  fdum(2) = aspect(1,i3) * f(2,i3)
                  fdum(3) = aspect(1,i3) * f(3,i3)
                  fdum(4) = aspect(1,i3) * f(4,i3)
                  fdum(5) = aspect(1,i3) * f(5,i3)
cC$OMP CRITICAL
                  vnw(1 ,ic0) = vnw(1 ,ic0) - fdum(1)    
                  vnw(2 ,ic0) = vnw(2 ,ic0) - fdum(5)    
                  vnw(j3,ic0) = vnw(j3,ic0) - fdum(2)    
                  vnw(j4,ic0) = vnw(j4,ic0) - fdum(3)    
                  vnw(j5,ic0) = vnw(j5,ic0) - fdum(4)    
cC$OMP END CRITICAL
                  ic1 = icell(2,i3)
                  fdum(1) = aspect(2,i3) * f(1,i3)
                  fdum(2) = aspect(2,i3) * f(2,i3)
                  fdum(3) = aspect(2,i3) * f(3,i3)
                  fdum(4) = aspect(2,i3) * f(4,i3)
                  fdum(5) = aspect(2,i3) * f(5,i3)
cC$OMP CRITICAL
                  vnw(1 ,ic1) = vnw(1 ,ic1) + fdum(1)
                  vnw(2 ,ic1) = vnw(2 ,ic1) + fdum(5)
                  vnw(j3,ic1) = vnw(j3,ic1) + fdum(2)
                  vnw(j4,ic1) = vnw(j4,ic1) + fdum(3)
                  vnw(j5,ic1) = vnw(j5,ic1) + fdum(4)
cC$OMP END CRITICAL
                enddo
              nif = 0 
            endif       
          enddo
        enddo
      ENDIF
c
      call Timing ( 1 , 0 )
c
      return
      end

c     ------------------------------------
      subroutine Sweep4 ( iDir, dt, Level )
c     ------------------------------------
c     
c     Input     :    iDir            - direction of sweep ( 1 - 6 )
c                    dt              - time step
c                    Level           - level of sweep in the tree
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 vd2(nhvar,-1:2,iSwChunk) ! input for Interface Solvers 
      real*8 g(2,iSwChunk)                 
      real*8 c(2,iSwChunk)                 
      real*8 f(5,iSwChunk)
      integer icell(2,iSwChunk)           
c
      integer iP, iR1, iR2, iL1, iL2
      integer iLvR1, iLvR2, iLvL1
      integer iKyR1, iopp, j3, j4, j5
c
      integer jpermute(7,6)
      data jpermute / 1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 3, 4, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 4, 3, 5, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7,
     &                1, 2, 5, 4, 3, 6, 7 /
c
      dimension aspect(2,iSwChunk), asp(3)
      data asp / 0.5, 0.25, 0.125 /
      integer iopposite(6)
      data iopposite / 2, 1, 4, 3, 6, 5 /
      real*8 c112, c122
      parameter ( c112 = 1. / 1.25, c122 = 1. / 1.5 )
c
      common / Sweep_Com_1 / vd2
      common / Sweep_Com_2 / g
      common / Sweep_Com_3 / c
      common / Sweep_Com_4 / f
      common / Sweep_Com_5 / icell
      common / Sweep_Com_6 / nif
c
c      call Timing(1,0)
c
      dx = CellSize(Level)   
      dtx = dt / dx
      dtx2 = 0.5*dtx
c              
      j3 = jpermute(3,iDir)
      j4 = jpermute(4,iDir)
      j5 = jpermute(5,iDir)
c
      nif = 0
c
      IF ( Level .eq. MinLevel ) THEN
        do i1 = 1 , ncell0
          iP = i1
          if ( iOctCh(iP) .eq. nil ) then   ! only leaves
            iopp  = iopposite(iDir)
            iR1   = iNb ( iP  , iDir )
            iR2   = iNb ( iR1 , iDir )
            iL1   = iNb ( iP  , iopp )
            iL2   = iNb ( iL1 , iopp )
            iLvL1 = iLv ( iL1 )
            iLvR1 = iLv ( iR1 )
            iLvR2 = iLv ( iR2 )
            if ( iOctCh(iR1) .eq. 0 ) then  
              iKyR1 = 0   ! leaf
            else
              iKyR1 = 1   ! split
            endif            
            if ( iKyR1 .eq. 0 ) then ! here we need to consider only =|= case
              nif = nif + 1
              vd2(1,-1, nif) = hvar(1,iL1)
              vd2(1, 0, nif) = hvar(1,iP )
              vd2(1, 1, nif) = hvar(1,iR1)
              vd2(1, 2, nif) = hvar(1,iR2)
              vd2(2,-1, nif) = hvar(2,iL1)
              vd2(2, 0, nif) = hvar(2,iP )
              vd2(2, 1, nif) = hvar(2,iR1)
              vd2(2, 2, nif) = hvar(2,iR2)
              vd2(3,-1, nif) = hvar(j3,iL1)
              vd2(3, 0, nif) = hvar(j3,iP )
              vd2(3, 1, nif) = hvar(j3,iR1)
              vd2(3, 2, nif) = hvar(j3,iR2)
              vd2(4,-1, nif) = hvar(j4,iL1)
              vd2(4, 0, nif) = hvar(j4,iP )
              vd2(4, 1, nif) = hvar(j4,iR1)
              vd2(4, 2, nif) = hvar(j4,iR2)
              vd2(5,-1, nif) = hvar(j5,iL1)
              vd2(5, 0, nif) = hvar(j5,iP )
              vd2(5, 1, nif) = hvar(j5,iR1)
              vd2(5, 2, nif) = hvar(j5,iR2)
              vd2(6,-1, nif) = hvar(6,iL1)
              vd2(6, 0, nif) = hvar(6,iP )
              vd2(6, 1, nif) = hvar(6,iR1)
              vd2(6, 2, nif) = hvar(6,iR2)
              vd2(7,-1, nif) = hvar(7,iL1)
              vd2(7, 0, nif) = hvar(7,iP )
              vd2(7, 1, nif) = hvar(7,iR1)
              vd2(7, 2, nif) = hvar(7,iR2)              
c
              icell(1, nif) = iP
              icell(2, nif) = iR1
              c(1, nif)     = 1.0
              if ( iLvL1 .lt. Level ) c(1, nif) = c112
              c(2, nif)     = 1.0
              if ( iLvR2 .lt. Level ) c(2, nif) = c112
            endif
          endif
c
          if ( (nif .eq. iSwChunk).or.(i1 .eq. ncell0) ) then 
            call FU_p ( dt, dx, dtx, dtx2 )
            if ( lapidus ) call FD_p ( iDir, dt, dx )
c           ! the loop below must be guarded in parallel execution
            do i2 = 1 , nif 
              ic0 = icell(1,i2)
              vnw(1 ,ic0) = vnw(1 ,ic0) - f(1,i2)    
              vnw(2 ,ic0) = vnw(2 ,ic0) - f(5,i2)    
              vnw(j3,ic0) = vnw(j3,ic0) - f(2,i2)    
              vnw(j4,ic0) = vnw(j4,ic0) - f(3,i2)    
              vnw(j5,ic0) = vnw(j5,ic0) - f(4,i2)    
              ic1 = icell(2,i2)
              vnw(1 ,ic1) = vnw(1 ,ic1) + f(1,i2)
              vnw(2 ,ic1) = vnw(2 ,ic1) + f(5,i2)
              vnw(j3,ic1) = vnw(j3,ic1) + f(2,i2)
              vnw(j4,ic1) = vnw(j4,ic1) + f(3,i2)
              vnw(j5,ic1) = vnw(j5,ic1) + f(4,i2)
            enddo
            nif = 0 
          endif       
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
        do i1 = 1 , nLevel
          ic = iSelect(i1)
          do i2 = 0 , 7
            idcell = ic + i2
            iP = idcell
            if ( iOctCh(iP) .eq. nil ) then   ! only leaves
              iopp  = iopposite(iDir)
              iR1   = iNb ( iP  , iDir )
              iR2   = iNb ( iR1 , iDir )
              iL1   = iNb ( iP  , iopp )
              iL2   = iNb ( iL1 , iopp )
              iLvL1 = iLv ( iL1 )
              iLvR1 = iLv ( iR1 )
              iLvR2 = iLv ( iR2 )
              if ( iOctCh(iR1) .eq. 0 ) then  
                iKyR1 = 0   ! leaf
              else
                iKyR1 = 1   ! split
              endif            
              if ( iKyR1 .eq. 0 ) then 
                if ( iLvR1 .eq. Level ) then  ! = | = case
                  nif = nif + 1
                  vd2(1,-1, nif) = hvar(1,iL1)
                  vd2(1, 0, nif) = hvar(1,iP )
                  vd2(1, 1, nif) = hvar(1,iR1)
                  vd2(1, 2, nif) = hvar(1,iR2)
                  vd2(2,-1, nif) = hvar(2,iL1)
                  vd2(2, 0, nif) = hvar(2,iP )
                  vd2(2, 1, nif) = hvar(2,iR1)
                  vd2(2, 2, nif) = hvar(2,iR2)
                  vd2(3,-1, nif) = hvar(j3,iL1)
                  vd2(3, 0, nif) = hvar(j3,iP )
                  vd2(3, 1, nif) = hvar(j3,iR1)
                  vd2(3, 2, nif) = hvar(j3,iR2)
                  vd2(4,-1, nif) = hvar(j4,iL1)
                  vd2(4, 0, nif) = hvar(j4,iP )
                  vd2(4, 1, nif) = hvar(j4,iR1)
                  vd2(4, 2, nif) = hvar(j4,iR2)
                  vd2(5,-1, nif) = hvar(j5,iL1)
                  vd2(5, 0, nif) = hvar(j5,iP )
                  vd2(5, 1, nif) = hvar(j5,iR1)
                  vd2(5, 2, nif) = hvar(j5,iR2)
                  vd2(6,-1, nif) = hvar(6,iL1)
                  vd2(6, 0, nif) = hvar(6,iP )
                  vd2(6, 1, nif) = hvar(6,iR1)
                  vd2(6, 2, nif) = hvar(6,iR2)
                  vd2(7,-1, nif) = hvar(7,iL1)
                  vd2(7, 0, nif) = hvar(7,iP )
                  vd2(7, 1, nif) = hvar(7,iR1)
                  vd2(7, 2, nif) = hvar(7,iR2)              
c
                  icell(1, nif) = iP
                  icell(2, nif) = iR1
                  c(1, nif)     = 1.0
                  if ( iLvL1 .lt. Level ) c(1, nif) = c112
                  c(2, nif)     = 1.0
                  if ( iLvR2 .lt. Level ) c(2, nif) = c112
                  aspect(1, nif) = 1.
                  aspect(2, nif) = 1.
                else
                  nif = nif + 1
                  vd2(1,-1, nif) = hvar(1,iL1)
                  vd2(1, 0, nif) = hvar(1,iP )
                  vd2(1, 1, nif) = hvar(1,iR1)
                  vd2(1, 2, nif) = hvar(1,iR2)
                  vd2(2,-1, nif) = hvar(2,iL1)
                  vd2(2, 0, nif) = hvar(2,iP )
                  vd2(2, 1, nif) = hvar(2,iR1)
                  vd2(2, 2, nif) = hvar(2,iR2)
                  vd2(3,-1, nif) = hvar(j3,iL1)
                  vd2(3, 0, nif) = hvar(j3,iP )
                  vd2(3, 1, nif) = hvar(j3,iR1)
                  vd2(3, 2, nif) = hvar(j3,iR2)
                  vd2(4,-1, nif) = hvar(j4,iL1)
                  vd2(4, 0, nif) = hvar(j4,iP )
                  vd2(4, 1, nif) = hvar(j4,iR1)
                  vd2(4, 2, nif) = hvar(j4,iR2)
                  vd2(5,-1, nif) = hvar(j5,iL1)
                  vd2(5, 0, nif) = hvar(j5,iP )
                  vd2(5, 1, nif) = hvar(j5,iR1)
                  vd2(5, 2, nif) = hvar(j5,iR2)
                  vd2(6,-1, nif) = hvar(6,iL1)
                  vd2(6, 0, nif) = hvar(6,iP )
                  vd2(6, 1, nif) = hvar(6,iR1)
                  vd2(6, 2, nif) = hvar(6,iR2)
                  vd2(7,-1, nif) = hvar(7,iL1)
                  vd2(7, 0, nif) = hvar(7,iP )
                  vd2(7, 1, nif) = hvar(7,iR1)
                  vd2(7, 2, nif) = hvar(7,iR2)              
                  icell(1,nif)  = iP
                  icell(2,nif)  = iR1
                  aspect(1,nif) = 1.
                  aspect(2,nif) = asp(ndim)
                  c(1,nif)      = c112
                  c(2,nif)      = c122
                endif
              endif
              if ( nif .eq. iSwChunk ) then 
                call FU_p ( dt, dx, dtx, dtx2 )
                if ( lapidus ) call FD_p ( iDir, dt, dx )
c               ! the loop below must be guarded in parallel execution
                do i3 = 1 , nif 
                  ic0 = icell(1,i3)
                  vnw(1 ,ic0) = vnw(1 ,ic0) - aspect(1,i3) * f(1,i3)    
                  vnw(2 ,ic0) = vnw(2 ,ic0) - aspect(1,i3) * f(5,i3)    
                  vnw(j3,ic0) = vnw(j3,ic0) - aspect(1,i3) * f(2,i3)    
                  vnw(j4,ic0) = vnw(j4,ic0) - aspect(1,i3) * f(3,i3)    
                  vnw(j5,ic0) = vnw(j5,ic0) - aspect(1,i3) * f(4,i3)    
                  ic1 = icell(2,i3)
                  vnw(1 ,ic1) = vnw(1 ,ic1) + aspect(2,i3) * f(1,i3)
                  vnw(2 ,ic1) = vnw(2 ,ic1) + aspect(2,i3) * f(5,i3)
                  vnw(j3,ic1) = vnw(j3,ic1) + aspect(2,i3) * f(2,i3)
                  vnw(j4,ic1) = vnw(j4,ic1) + aspect(2,i3) * f(3,i3)
                  vnw(j5,ic1) = vnw(j5,ic1) + aspect(2,i3) * f(4,i3)
                enddo
                nif = 0 
              endif       
              if ( iLvL1 .eq. Level-1 ) then 
                nif = nif + 1
                vd2(1,-1, nif) = hvar(1,iL2)
                vd2(1, 0, nif) = hvar(1,iL1)
                vd2(1, 1, nif) = hvar(1,iP )
                vd2(1, 2, nif) = hvar(1,iR1)
                vd2(2,-1, nif) = hvar(2,iL2)
                vd2(2, 0, nif) = hvar(2,iL1)
                vd2(2, 1, nif) = hvar(2,iP )
                vd2(2, 2, nif) = hvar(2,iR1)
                vd2(3,-1, nif) = hvar(j3,iL2)
                vd2(3, 0, nif) = hvar(j3,iL1)
                vd2(3, 1, nif) = hvar(j3,iP )
                vd2(3, 2, nif) = hvar(j3,iR1)
                vd2(4,-1, nif) = hvar(j4,iL2)
                vd2(4, 0, nif) = hvar(j4,iL1)
                vd2(4, 1, nif) = hvar(j4,iP )
                vd2(4, 2, nif) = hvar(j4,iR1)
                vd2(5,-1, nif) = hvar(j5,iL2)
                vd2(5, 0, nif) = hvar(j5,iL1)
                vd2(5, 1, nif) = hvar(j5,iP )
                vd2(5, 2, nif) = hvar(j5,iR1)
                vd2(6,-1, nif) = hvar(6,iL2)
                vd2(6, 0, nif) = hvar(6,iL1)
                vd2(6, 1, nif) = hvar(6,iP )
                vd2(6, 2, nif) = hvar(6,iR1)
                vd2(7,-1, nif) = hvar(7,iL2)
                vd2(7, 0, nif) = hvar(7,iL1)
                vd2(7, 1, nif) = hvar(7,iP )
                vd2(7, 2, nif) = hvar(7,iR1)
                icell (1, nif) = iL1
                icell (2, nif) = iP
                aspect(1, nif) = asp(ndim)
                aspect(2, nif) = 1.
                c(1, nif) = c122
                c(2, nif) = c112
              endif
            endif
            if (         (nif .eq. iSwChunk)
     &                        .or. 
     &           ((i1.eq. nLevel) .and. (i2.eq.7)) ) then 
              call FU_p ( dt, dx, dtx , dtx2 )
              if ( lapidus ) call FD_p ( iDir, dt, dx )
c             ! the loop below must be guarded in parallel execution
              do i3 = 1 , nif 
                ic0 = icell(1,i3)
                vnw(1 ,ic0) = vnw(1 ,ic0) - aspect(1,i3) * f(1,i3)    
                vnw(2 ,ic0) = vnw(2 ,ic0) - aspect(1,i3) * f(5,i3)    
                vnw(j3,ic0) = vnw(j3,ic0) - aspect(1,i3) * f(2,i3)    
                vnw(j4,ic0) = vnw(j4,ic0) - aspect(1,i3) * f(3,i3)    
                vnw(j5,ic0) = vnw(j5,ic0) - aspect(1,i3) * f(4,i3)    
                ic1 = icell(2,i3)
                vnw(1 ,ic1) = vnw(1 ,ic1) + aspect(2,i3) * f(1,i3)
                vnw(2 ,ic1) = vnw(2 ,ic1) + aspect(2,i3) * f(5,i3)
                vnw(j3,ic1) = vnw(j3,ic1) + aspect(2,i3) * f(2,i3)
                vnw(j4,ic1) = vnw(j4,ic1) + aspect(2,i3) * f(3,i3)
                vnw(j5,ic1) = vnw(j5,ic1) + aspect(2,i3) * f(4,i3)
              enddo
              nif = 0 
            endif       
          enddo
        enddo
      ENDIF
c
c      call Timing ( 1 , 0 )
c
      return
      end
