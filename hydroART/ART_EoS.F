C     =====================================================================
c                                                                         .
c          ART Version 1-0 : Equation of state for the gas                .
c                                                                         .
c         by Alexei Khokhlov and Andrey Kravtsov (1996-1999)              .
c                                                                         .   
c     =====================================================================
c
c     ----------------------------
      subroutine TreeEoS ( Level )
c     ----------------------------
c
c     equation of state for (currently ideal) gas
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 eta1
      integer Level 
      integer i1, i2, idcell, icell, nLevel
c
      gamma1 = gamma - 1.

      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
        do i1 = 1 , ncell0 
          if ( iOctCh(i1) .eq. nil ) then  ! only leaves    
             hvar(6,i1) = (hvar(7,i1)-1.0) * hvar(8,i1) 
          endif
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idc)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idc = icell + i2
            if ( iOctCh(idc) .eq. nil ) then   ! only leaves
c  
               hvar(6,idc) = (hvar(7,idc)-1.0)*hvar(8,idc)
c
            endif
          enddo
        enddo
      ENDIF
c
      return
      end
c
c     --------------------------------------
      subroutine Smooth_Gamma_eff  ( Level )
c     --------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
      parameter ( w0 = 1.0, w1 = 1.0, w01 = 1./(w0 + 6.*w1) )
      real*8 wsum1,wsum2,wsum3,wsum4,wsum5,wsum6
      integer Nbb(6)

      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i2,Nbb,i3)
C$OMP+PRIVATE(wsum1)
        do i1 = 1 , ncell0
          call NbAll0 ( i1 , Nbb )
          wsum1 = w0 * hvar(7,i1)
          do i2 = 1 , 6
            wsum1 = wsum1 + w1 * hvar(7,Nbb(i2))
          enddo
          vnw(1,i1) = wsum1 * w01
        enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
        do i1 = 1 , ncell0
          hvar(7,i1) = min(gamma,max(1.01,vnw(1,i1)))
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel )
       do i = 1 , 1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,idcell,i2,Nbb,i3)
C$OMP+PRIVATE(wsum1)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
              call NbAllL ( idcell , Nbb )
              wsum1 = w0 * hvar(7,idcell)
              do i3 = 1 , 6
                wsum1 = wsum1 + w1 * hvar(7,Nbb(i3))
              enddo
              vnw(1,idcell) = wsum1 * w01
          enddo
        enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,idcell,i2)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            hvar(7,idcell) = min ( gamma ,
     &                             max(1.01,vnw(1,idcell)))
          enddo
        enddo
       enddo
      ENDIF

      return
      end

c
c     ---------------------------------------------------------
      subroutine Estimate_Effective_Gamma ( Level , iWorkType )
c     ---------------------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
c
      if ( iWorkType .eq. 0 ) then
        call EEG_Work0 ( Level ) 
        return 
      endif
      if ( iWorkType .eq. 1 ) then
        call EEG_Work1 ( Level )  
        return
      endif

      return
      end
c
c     ------------------------------
      subroutine EEG_Work0 ( Level )
c     ------------------------------
      include 'a_tree.h'
      include 'a_control.h'
c
      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then 
            rhophi(1,i1) = hvar(1,i1) 
            rhophi(2,i1) = hvar(8,i1)
          endif 
        enddo 
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel )
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(ic,indc,i2)
C$OMP+PRIVATE(iP,iT)
        do i1 = 1 , nLevel
          ic = iSelect(i1) - 1
          indc = (i1-1) * nchild
          do i2 = 1 , nchild  
            iP = ic   + i2
            iT = indc + i2
            if ( iOctCh(iP) .eq. nil ) then 
              rhophi(1,iT) = hvar(1,iP)
              rhophi(2,iT) = hvar(8,iP)
            endif
          enddo
        enddo 
      ENDIF 
c
      return
      end
c
c     ------------------------------
      subroutine EEG_Work1 ( Level )
c     ------------------------------
      include 'a_tree.h'
      include 'a_control.h'
c
      parameter ( gammin = 1.01 ) 
c
      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(g_eff)
        do i1 = 1 , ncell0
          if ( iOctCh(i1) .eq. nil ) then
            g_eff = (rhophi(2,i1) - hvar(8,i1))  * hvar(1,i1) / 
     &              ((rhophi(1,i1) - hvar(1,i1)) * hvar(8,i1))  
            g_eff = max ( gammin , g_eff )
            hvar(7,i1) = min ( gamma  , g_eff )
          endif
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel )
cC$OMP PARALLEL DO DEFAULT(SHARED)
cC$OMP+PRIVATE(ic,indc,i2)
cC$OMP+PRIVATE(iP,iT,g_eff)
       if ( Level .eq. 1 ) then 
        open ( 13 , file = 'g_eff.dat' )
       endif
        do i1 = 1 , nLevel
          ic = iSelect(i1) - 1
          indc = (i1-1) * nchild
          do i2 = 1 , nchild
            iP = ic   + i2
            iT = indc + i2
            if ( iOctCh(iP) .eq. nil ) then
              g_eff =  (rhophi(2,iT) - hvar(8,iP)) * hvar(1,iP) /
     &                ((rhophi(1,iT) - hvar(1,iP)) * hvar(8,iP))
           if ( Level .eq. 1 ) then 
           write(13,'(5(g12.6,1x))') rhophi(1,iT), hvar(1,iP),
     &        rhophi(2,iT), hvar(8,iP), g_eff
           endif

              g_eff = max ( gammin , g_eff )
              hvar(7,iP) = min ( gamma  , g_eff )
            endif
          enddo
        enddo
        if ( Level .eq. 1 ) then 
	  close(13)
          stop
        endif
      ENDIF
c
      return
      end

