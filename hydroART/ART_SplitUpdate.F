c     =====================================================================
c                                                                         
c          ART Version 3-1 : update split cells by averaging variables 
c                            of their children             
c                                                                         
c            by Alexei Khokhlov and Andrey Kravtsov (1996-1999)           
c                                                                         .  
c     =====================================================================

c     --------------------------------
      subroutine SplitUpdate ( Level ) 
c     --------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
      integer Level 
      integer i1, i2, iv, jch, icell, idcell, iCh1, iChild, nLevel
c
      real    factor
      parameter ( factor = 0.5**ndim )
c
      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iCh1,iv,jch,iChild)
        do i1 = 1 , ncell0
          iCh1 = iOctCh(i1)
          if ( iCh1 .gt. nil ) then   ! only split cells
            iCh1 = iCh1*nchild + mbshift
            do iv = 1 , nhvar   ! 1st child separate to avoid zeroing 
              hvar(iv,i1) = hvar(iv,iCh1)
            enddo            
            do jch = 1 , nchild-1
              iChild = iCh1 + jch
              do iv = 1 , nhvar 
                hvar(iv,i1) = hvar(iv,i1) + hvar(iv,iChild)
              enddo
            enddo
           do iv = 1 , nhvar
             hvar(iv,i1) = factor * hvar(iv,i1)
           enddo
          endif
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idcell,iCh1,iv,jch,iChild)
        do i1 = 1 , nLevel
          icell = iSelect(i1)
          do i2 = 0 , 7
            idcell = icell + i2
            iCh1 = iOctCh(idcell) 
            if ( iCh1 .gt. nil ) then   ! only split cells
              iCh1 = iCh1*nchild + mbshift
              do iv = 1 , nhvar   ! 1st child separate to avoid zeroing 
                hvar(iv,idcell) = hvar(iv,iCh1)
              enddo
              do jch = 1 , nchild-1
                iChild = iCh1 + jch
                do iv = 1 , nhvar 
                  hvar(iv,idcell) = hvar(iv,idcell) + hvar(iv,iChild)
                enddo
              enddo
              do iv = 1 , nhvar
                hvar(iv,idcell) = factor * hvar(iv,idcell)
              enddo              
            endif
          enddo
        enddo
      ENDIF
c
      return
      end
