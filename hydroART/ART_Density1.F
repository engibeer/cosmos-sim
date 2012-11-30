c     ------------------------------------------------------
      subroutine Assign_Density ( Level1 , Level2 )
c     ------------------------------------------------------
c
c     purpose: assigns density to cells on levels Level1-Level2
c      Level1 <= Level2
c
c
      integer Level, Level1, Level2
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
        call Assign_GasDensity ( Level )
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
      use density_alloc
      integer Level 
      double precision Size3
      include 'a_tree.h'
      include 'a_control.h'

      IF ( Level .eq. MinLevel ) THEN
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i_TH)
        do ic1 = 1 , ncell0
          var(1,ic1) = -1.0
          ref(ic1) = 0.
          do i_TH=1,NCPUS
             vref(i_TH,ic1) = 0.
             vvar(i_TH,ic1) = 0.
          enddo
        enddo
      ELSE
c
        Size3 = -CellSize(Level)**3
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i,i_TH)
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
          do i=icell, icell+7
             do i_TH=1,NCPUS
                vref(i_TH,i) = 0.
                vvar(i_TH,i) = 0.
             enddo
          enddo
       enddo
      ENDIF
      
      return
      end
c
c     -----------------------------------------------
      subroutine Assign_GasDensity ( Level )
c     -----------------------------------------------
c
c     purpose: initializes density on a given Level with gas density
c
      integer Level  
      double precision Size3
      include 'a_tree.h'
      include 'a_control.h'
c
      IF ( Level .eq. MinLevel ) THEN
c
C$OMP PARALLEL DO DEFAULT(SHARED)
          do ic1 = 1 , ncell0
            var(1,ic1) = var(1,ic1) + hvar(1,ic1)
          enddo
c
      ELSE
c
        Size3 = CellSize(Level)**3
        call Select_Cells ( Level , nLevel ) 
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
      use density_alloc
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      include 'a_dnew.h' 
      integer Level, Level1, Level2
      double precision xx, yy, zz, xdum, ydum, zdum
      real*8 xx1, yy1, zz1, d1,d2,d3,t1,t2,t3,t2w,d2w
      

      if ((Level1.lt.MinLevel).or.(Level1.gt.MaxLevel)) return          
      if ((Level2.lt.MinLevel).or.(Level2.gt.MaxLevel)) return          
      if ( Level1.eq.MinLevel ) then
        Ldummy = Level1 + 1
        call Assign_ParticleDensity0 ()
      else
        Ldummy = Level1
      endif
      DO Level = Ldummy , Level2 

        Sizei  = CellSizei(Level)
        hSize  = CellSize(Level) * 0.5
        Size   = CellSize(Level)
C---SERIAL LOOP: 
        npch = 0
        DO ip = 1 , np
c
c....     preselect NpCH particles that can contribute density to this Level
c           
           if ( iPL(ip) .ge. Level ) then 
              npch = npch + 1         
              ipc(npch) = ip
           endif 
        ENDDO
CEVERINO01202009
c        NCdensity=0
c-------------------------------------------------------------------- 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i,iPar,id,xx,yy,zz,wp,Posx,Posy,Posz,xdum,ydum,zdum)
C$OMP+PRIVATE(iac1,iac2,iac3,iac4,iac5,iac6,iac7,iac8,xx1,yy1,zz1)
C$OMP+PRIVATE(iCIC,jCIC,kCIC,d1,d2,d3,t1,t2,t3,t2w,d2w,dpcScalar)
C$OMP+PRIVATE(iREF,npch_TH,ip_first,ip_last)
        DO i_TH = 1 , NCPUS
           npch_TH=int(1.0*npch/NCPUS)
           if (i_TH.eq.NCPUS) then
              npch_TH=npch- (NCPUS-1)*npch_TH
              ip_first=npch-npch_TH+1
              ip_last=npch
           else
              ip_first= (i_TH-1) * npch_TH + 1
              ip_last=  ip_first + npch_TH - 1
           endif
c
           do i=ip_first,ip_last
c
               id = ipc(i)
               xx = x(id)
               yy = y(id)
               zz = z(id)
               wp = pw(id)

#ifdef STARFORM
                  if ( id .le. nsp(1,2) .or. 
     &                 id .ge. nsp(nspecies,1) ) then 
#endif
#ifndef STARFORM 
                  if ( id .le. nsp(1,2) ) then 
#endif
                     iREF=1
                  else
                     iREF=0
                  endif
c     
c....         find neighbors to which this particle can contribute
c              
               xdum = xx - hSize              
               ydum = yy - hSize
               zdum = zz - hSize
c
CEVERINO01202009: check periodic boundaries is not needed
!               if ( xdum .le. 1.0 ) xdum = xdum + yn
!               if ( ydum .le. 1.0 ) ydum = ydum + yn
!               if ( zdum .le. 1.0 ) zdum = zdum + yn
CEVERINO01202009
c               IF (iac1.eq.iFindCell0( Level ,xdum,ydum,zdum )) then
c                  if (iac1.ne.0) then
cC$OMP ATOMIC
c                     NCdensity=NCdensity+1
c                  endif
c                ENDIF
c--------------------------------------------------------------------
               iac1 = iFindCell0 ( Level , xdum, ydum, zdum ) 

               xdum = xx + hSize
               ydum = yy - hSize
               zdum = zz - hSize
!               if ( xdum .ge.  xn ) xdum = xdum - yn
!               if ( ydum .le. 1.0 ) ydum = ydum + yn
!               if ( zdum .le. 1.0 ) zdum = zdum + yn
               iac2 = iFindCell0 ( Level , xdum, ydum, zdum ) 

               xdum = xx - hSize
               ydum = yy + hSize
               zdum = zz - hSize
!              if ( xdum .le. 1.0 ) xdum = xdum + yn
!              if ( ydum .ge. xn  ) ydum = ydum - yn
!              if ( zdum .le. 1.0 ) zdum = zdum + yn
               iac3 = iFindCell0 ( Level , xdum, ydum, zdum ) 
               
               xdum = xx + hSize
               ydum = yy + hSize
               zdum = zz - hSize
!               if ( xdum .ge. xn  ) xdum = xdum - yn
!               if ( ydum .ge. xn  ) ydum = ydum - yn
!               if ( zdum .le. 1.0 ) zdum = zdum + yn
               iac4 = iFindCell0 ( Level , xdum, ydum, zdum ) 
               
               xdum = xx - hSize
               ydum = yy - hSize
               zdum = zz + hSize
!               if ( xdum .le. 1.0 ) xdum = xdum + yn
!               if ( ydum .le. 1.0 ) ydum = ydum + yn
!               if ( zdum .ge. xn  ) zdum = zdum - yn
               iac5 = iFindCell0 ( Level , xdum, ydum, zdum ) 

               xdum = xx + hSize
               ydum = yy - hSize
               zdum = zz + hSize
!               if ( xdum .ge. xn  ) xdum = xdum - yn
!               if ( ydum .le. 1.0 ) ydum = ydum + yn
!               if ( zdum .ge. xn  ) zdum = zdum - yn
               iac6 = iFindCell0 ( Level , xdum, ydum, zdum ) 

               xdum = xx - hSize
               ydum = yy + hSize
               zdum = zz + hSize
!               if ( xdum .le. 1.0 ) xdum = xdum + yn
!               if ( ydum .ge. xn  ) ydum = ydum - yn
!               if ( zdum .ge. xn  ) zdum = zdum - yn
               iac7 = iFindCell0 ( Level , xdum, ydum, zdum ) 

               xdum = xx + hSize
               ydum = yy + hSize
               zdum = zz + hSize
!               if ( xdum .ge. xn  ) xdum = xdum - yn
!               if ( ydum .ge. xn  ) ydum = ydum - yn
!               if ( zdum .ge. xn  ) zdum = zdum - yn
               iac8 = iFindCell0 ( Level , xdum, ydum, zdum ) 
c
c
               xx1=xx*Sizei
               yy1=yy*Sizei
               zz1=zz*Sizei
               iCIC= int(xx1-0.500) 
               jCIC= int(yy1-0.500)
               kCIC= int(zz1-0.500)
               d1=xx1-iCIC - 0.500
               d2=yy1-jCIC - 0.500
               d3=zz1-kCIC - 0.500
               t1  = 1.0 - d1
               t2  = 1.0 - d2
               t3  = 1.0 - d3
               t2w = t2 * wp
               d2w = d2 * wp
c
c   8 contributions
c
               if ( iac1 .gt. nil ) then
                  dpcScalar= t3 * t1 * t2w
                  vvar(i_TH,iac1)= vvar(i_TH,iac1)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac1)= vref(i_TH,iac1)+dpcScalar
                  endif
               endif


               if ( iac2 .gt. nil ) then 
                  dpcScalar= t3 * d1 * t2w
                  vvar(i_TH,iac2)= vvar(i_TH,iac2)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac2)= vref(i_TH,iac2)+dpcScalar
                  endif
               endif

               if ( iac3 .gt. nil ) then 
                  dpcScalar=  t3 * t1 * d2w
                  vvar(i_TH,iac3)= vvar(i_TH,iac3)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac3)= vref(i_TH,iac3)+dpcScalar
                  endif
               endif

               if ( iac4 .gt. nil ) then 
                  dpcScalar= t3 * d1 * d2w
                  vvar(i_TH,iac4)= vvar(i_TH,iac4)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac4)= vref(i_TH,iac4)+dpcScalar
                  endif
               endif

               if ( iac5 .gt. nil ) then 
                  dpcScalar= d3 * t1 * t2w
                  vvar(i_TH,iac5)= vvar(i_TH,iac5)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac5)= vref(i_TH,iac5)+dpcScalar
                  endif
               endif

               if ( iac6 .gt. nil ) then 
                  dpcScalar= d3 * d1 * t2w
                  vvar(i_TH,iac6)= vvar(i_TH,iac6)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac6)= vref(i_TH,iac6)+dpcScalar
                  endif
               endif

               if ( iac7 .gt. nil ) then 
                  dpcScalar= d3 * t1 * d2w
                  vvar(i_TH,iac7)= vvar(i_TH,iac7)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac7)= vref(i_TH,iac7)+dpcScalar
                  endif
               endif

               if ( iac8 .gt. nil ) then 
                  dpcScalar= d3 * d1 * d2w              
                  vvar(i_TH,iac8)= vvar(i_TH,iac8)+dpcScalar
                  if (iREF.eq.1) then
                     vref(i_TH,iac8)= vref(i_TH,iac8)+dpcScalar
                  endif
               endif
            enddo 
         ENDDO    
! parallel update of global arrays:
         call Select_Cells ( Level , nLevel )
C$OMP PARALLEL DO  DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idcell, i_TH)
         DO i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idcell = icell + i2
               do i_TH=1,NCPUS
                  var(1,idcell)=var(1,idcell)+vvar(i_TH,idcell)
                  ref(idcell)=ref(idcell)+vref(i_TH,idcell)
               enddo
            enddo
         ENDDO
      ENDDO
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
c
