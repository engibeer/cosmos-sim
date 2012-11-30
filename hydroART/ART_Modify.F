c     =====================================================================
c                                                                         
c              HART Version 1:  Refinement/Derefinement                   
c                                                                         
c          by Alexei Khokhlov and Andrey Kravtsov (1997-1999)           
c                                                                            
c     =====================================================================
c
c     -----------------------------------------
      subroutine Modify ( Level , iJob , mtot )
c     -----------------------------------------
c
c     purpose: Makes refinement/derefinement for Level
c
c     input  : int Level - level to process
c              int iJob  - kind of job to be done:
c                             1 - refinement only
c                             0 - refinement & derefinement
c                            -1 - derefinement only
c
c     output : int mtot - total number of refined/derefined cells
c          
      include 'a_tree.h'
c
      integer Level 
      integer mtot         ! total # of cell marked to split 
      integer nmark        ! # of marked cells
c
      mtot  = nil 
      nmark = nil 
c
      if ( iJob .eq. 1 ) then
          call Refine ( Level , mtot )
      else
        if ( iJob .eq. 0 ) then
          call Refine   ( Level , mtot )
          call DeRefine ( Level , mtot )          
        else
            call Open_ASCII_File ( iOErr , errorfile , sequent )
            write (iOErr,*) '1  *Modify error: bad iJob:',iJob
            close ( iOErr )
            stop
        endif
      endif      

      if ( iNOLL(Level+1)*nchild .ge. nclmax ) then 
        call Open_ASCII_File ( iOErr , errorfile , sequent ) 
        write(iOErr,*) '1  *error: in Modify: after refinement:'
        write(iOErr,*) 'L =',Level+1,' Ncells(L)=',iNOLL(Level+1)*nchild
        write(iOErr,*) ' Ncells > nclmax =',nclmax
        write(iOErr,*) ' increase nclmax in a_setup.h and rerun.'
        close(iOErr)
        stop
      endif
c
      return
      end
c
c     ----------------------------------
      subroutine Refine ( Level , mtot )
c     ----------------------------------
c
c     purpose: performs refinement of cells of level Level
c     input  : Level - level to process
c     output : mtot  - total number of refined/derefined cells
c
c
c      include 'a_setup.h'
      include 'a_tree.h'
c
      integer Level
      integer mtot         ! total # of cell marked to split 
c
c.... mark cells for refinement 
c
      call Mark_Cells_to_Refine ( nmark , Level )
c
c.... split all marked cells       
c
      if ( Level .lt. MaxLevel ) then 
        call Split ( Level , mtot )
      endif
c
      return
      end
c
c     -------------------------------------------------
      subroutine Mark_Cells_to_Refine ( nmark , Level )
c     -------------------------------------------------
c
c     Computes weights for refinement/derefinement
c     according to distribution of DM particles
c
c.... get refinement/derefinement indicators for specified levels
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level , nmark      
c
      call Zero_Indicators            ( Level )

      call Mark_Refinement_Indicators ( Level )

      Ndiffusion = Ndiff

      do iDiff = 1 , NDiffusion
        call Diffuse_Indicators       ( Level ) 
      enddo      
c
      call CheckRefinementMask ( Level )
c
      return
      end
c
c     ------------------------------------
      subroutine Zero_Indicators ( Level ) 
c     ------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
        do ic1 = 1 , ncell0 
          vnw(1,ic1) = zero
          vnw(2,ic1) = zero
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(icell,ic2,idcell)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            vnw(1,idcell) = zero
            vnw(2,idcell) = zero
          enddo
        enddo
      ENDIF
c
      return
      end
c
c     -----------------------------------------------
      subroutine Mark_Refinement_Indicators ( Level ) 
c     -----------------------------------------------
c
c     purpose: to prepare working space for 
c              marking. var(1,*) does not have to be zeroed 
c              because it is used to construct original indicators
c
      include 'a_tree.h'
      include 'a_control.h'
c
      parameter ( epsP = 10. , epsR = 0.6 , epsS = 10. )
      parameter ( epsdrho = 1.0 )
c
      parameter ( fl_jeans = 4.0 ) ! # of cells per Jeans wavelength
c
      integer iNbC(nchild,3), iNC(27), Neighbors(nneigh)
c
      data iNbC / 2, 1, 4, 3, 6, 5, 8, 7,   
     &            3, 4, 1, 2, 7, 8, 5, 6,
     &            5, 6, 7, 8, 1, 2, 3, 4  /
c
      logical lrefok(n_ref)
c
      do iref = 1 , n_ref
        lrefok(iref) = .false. 
        if ( lmin_ref(iref) .le. Level .and. 
     &       lmax_ref(iref) .gt. Level ) lrefok(iref) = .true. 
      enddo

      nmark  = nil       
c
      IF ( Level .eq. MinLevel ) THEN 
c
c.... total mass
c
        if ( w_ref(1) .gt. 0. .and. lrefok(1) ) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rhocell)
          do ic1 = 1 , ncell0 
c            rhocell = (var(1,ic1) + 1. - hvar(1,ic1))/trho(Level)
             rhocell = ref(ic1)/trho(MinLevel)
             ref(ic1) = 0.
c            if ( rhocell .gt. trho(Level) ) then 
              vnw(1,ic1) = max ( min(rhocell,w_ref(1)) , 0. )
c            endif
          enddo
        endif
        if ( w_ref(2) .gt. 0. .and. lrefok(2) ) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dmbar)
          do ic1 = 1 , ncell0 
            dmbar = hvar(1,ic1) / tbmass(MinLevel)
c            if ( dmbar .gt. 1.0 ) then 
              vnw(1,ic1) = max( vnw(1,ic1) , min(dmbar,w_ref(2)) )
c            endif
          enddo
        endif
        if ( w_ref(8) .gt. 0. .and. lrefok(8) ) then
          fw_jeans = fl_jeans /
     &               sqrt(2.0*pi**2/3.0*(gamma/(gamma-1.0))/aexpn) 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(w_jeans)
          do ic1 = 1 , ncell0
             w_jeans = fw_jeans * hvar(1,ic1) * sqrt(1.0 / hvar(8,ic1))
             vnw(1,ic1) = max( vnw(1,ic1) , min(w_jeans,w_ref(8)) )
          enddo
        endif
c
c....   potential
c
        if ( w_ref(9) .gt. 0. .and. lrefok(9) ) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dmbar)
          do ic1 = 1 , ncell0
            dmbar = max ( var(2,ic1) / tphi(MinLevel) , 0.0 )
            dmbar = min ( dmbar , w_ref(9) )
c            if ( dmbar .gt. 1.0 ) then
              vnw(1,ic1) = max ( vnw(1,ic1) , dmbar )
c            endif
          enddo
        endif

c
        IF ( w_ref(3).gt.0. .or. w_ref(4).gt.0. .or. w_ref(5).gt.0. .or.
     &       w_ref(6).gt.0. .or. w_ref(7).gt.0. ) THEN 
        do j = 1 , 3   ! analyze three spatial dimensions
           j0 = 2*j - 1
           j1 = j0 + 1
           jv = 2 + j
C.... OpenMP
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iNb0,iNb1,rhoC,rhoC0,rhoC1,preC,preC1,preC0)
C$OMP+PRIVATE(entC,entC0,entC1,dp,dv,Sind,drho,Cind)
C$OMP+PRIVATE(Rhoind0,Rhoind1,Rhoind,Preind0,Preind1,Preind)
C$OMP+PRIVATE(Entind0,Entind1,Entind)
c
           do ic1 = 1 , ncell0 
c
             iNb0  = iNb(ic1,j0)
             iNb1  = iNb(ic1,j1)
             rhoC  = hvar(1,ic1)
             rhoC0 = hvar(1,iNb0)
             rhoC1 = hvar(1,iNb1)
             preC  = hvar(6,ic1)
             preC0 = hvar(6,iNb0)
             preC1 = hvar(6,iNb1)
             entC  = preC  / rhoC**gamma
             entC0 = preC0 / rhoC0**gamma
             entC1 = preC1 / rhoC1**gamma
c
             dp = abs ( preC0 - preC1 ) / min ( preC0, preC1 )
             dv = abs(hvar(jv,iNb0) / rhoC0 - hvar(jv,iNb1) / rhoC1)
c             dp = abs ( preC1 - preC ) / min ( preC, preC1 )
c             dv = abs(hvar(jv,ic1) / rhoC - hvar(jv,iNb1) / rhoC1)
             
c
c....        shock indicator
c
             if ( (dp .gt. epsP) .and. 
     &            (dv .gt. 0.  ) .and. 
     &            (rhoC .gt. floor_ref(3)) .and. lrefok(3) ) 
     &       then
               Sind = w_ref(3)
             else
               Sind = 0.
             end if
c
c....        contact discontinuity indicator
c             
             drho = abs ( rhoC0 - rhoC1 ) / min ( rhoC0, rhoC1 )
c             drho = abs ( rhoC - rhoC1 ) / min ( rhoC, rhoC1 )
             if ( (drho .gt. epsR) .and. 
     &            (dp   .lt. epsP) .and.
     &            (rhoC .gt. floor_ref(4)) .and. lrefok(4) ) 
     &       then
               Cind = w_ref(4)
             else
               Cind = 0.
             end if
c
c....        density gradient indicator 
c
             Rhoind0 = abs ( rhoC - rhoC0 ) /
     &                 min ( rhoC , rhoC0 ) 
             Rhoind1 = abs ( rhoC - rhoC1 ) /
     &                 min ( rhoC , rhoC1 )
             drho = max ( drho , Rhoind0 , Rhoind1 ) 
             if ( drho .gt. epsdrho .and. 
     &            rhoC .gt. floor_ref(5) .and. lrefok(5) ) then
               Rhoind  = w_ref(5) 
             else
               Rhoind  = 0.0 
             end if
c
c....        pressure gradient indicator 
c
             Preind0 = abs ( preC - preC0 ) / ( preC + preC0 ) 
             Preind1 = abs ( preC - preC1 ) / ( preC + preC1 ) 
             if ( rhoC .gt. floor_ref(6) .and. lrefok(6) ) then
               Preind  = w_ref(6) * max ( Preind0 , Preind1 )
             else
               Preind  = 0.
             end if
c
c....        entropy gradient indicator 
c
             Entind0 = abs ( entC - entC0 ) / ( entC + entC0 ) 
             Entind1 = abs ( entC - entC1 ) / ( entC + entC1 ) 
             dS = abs( ( entC0 - entC1 ) / min ( entC0, entC1 ))

             if ( rhoC .gt. floor_ref(7) .and. lrefok(7) .and. 
     &            dS .gt. epsS ) then 
               Entind  = w_ref(7) * max ( Entind0 , Entind1 )
             else
               Entind  = 0.
             end if
c
             vnw(1,ic1) = max (vnw(1,ic1), Sind, Cind, Rhoind, 
     &                                     Preind, Entind      )
           enddo ! end ic1
        enddo ! end j           
        ENDIF
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C
        vc = CellSize(Level)**3
        if ( w_ref(1) .gt. 0. .and. lrefok(1) ) then 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,ic2,idcell,rhocell)
           do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idcell = icell + ic2
c              rhocell = (var(1,idcell) + vc - hvar(1,idcell)*vc) / 
c     &                   trho(Level)
               rhocell = ref(idcell)/trho(Level)
               ref(idcell) = 0.
c              if ( rhocell .gt. trho(Level) ) then 
                vnw(1,idcell) = max ( min(rhocell,w_ref(1)) , 0. )
c              endif
            enddo
          enddo
        endif
        
        if ( w_ref(2) .gt. 0. .and. lrefok(2) ) then 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,ic2,idcell,dmbar)
           do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7 
              idcell = icell + ic2
              dmbar  = hvar(1,idcell) * vc / tbmass(Level)
c              if ( dmbar .gt. 1.0 ) then 
                vnw(1,idcell) = max( vnw(1,idcell),
     &                               min(w_ref(2),dmbar) )
c              endif
            enddo
          enddo
        endif
c
        if ( w_ref(8) .gt. 0. .and. lrefok(8) ) then
          fw_jeans = CellSize(Level) * fl_jeans /  
     &               sqrt(2.0*pi**2/3.0*(gamma/(gamma-1.0))/aexp(Level))  
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,ic2,idc,w_jeans)
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idc = icell + ic2
              w_jeans = fw_jeans *
     &                  hvar(1,idc) * sqrt(1.0 / hvar(8,idc))
              vnw(1,idc) = max( vnw(1,idc) , min(w_jeans,w_ref(8)) )
            enddo
          enddo
        endif
c
c....   potential
c
        if ( w_ref(9) .gt. 0. .and. lrefok(9) ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,ic2,idcell,dmbar)
           do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idcell = icell + ic2
              dmbar = max ( var(2,idcell) / tphi(Level) , 0.0 )
              dmbar = min ( dmbar , w_ref(9) )
              vnw(1,idcell) = max ( vnw(1,idcell) , dmbar )
            enddo
          enddo
        endif

        IF ( w_ref(3).gt.0. .or. w_ref(4).gt.0. .or. w_ref(5).gt.0. .or.
     &       w_ref(6).gt.0. .or. w_ref(7).gt.0. ) THEN 
        do j = 1 , 3   ! analyze three spatial dimensions
          j0 = 2*j - 1
          j1 = j0 + 1
          jv = 2 + j
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,ic2,idcell,iNb0,iNb1,rhoC,rhoC0,rhoC1)
C$OMP+PRIVATE(preC,preC1,preC0,entC,entC0,entC1,dp,dv,drho)
C$OMP+PRIVATE(Rhoind0,Rhoind1,Rhoind,Preind0,Preind1,Preind)
C$OMP+PRIVATE(Entind0,Entind1,Entind,Sind,Cind)
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idcell = icell + ic2
c
              iNb0  = iNb(idcell,j0)
              iNb1  = iNb(idcell,j1)
c
              if ( iLv(iNb0).eq.Level .and. iLv(iNb1).eq.Level ) then
                 if ( iLv(iNb(iNb0,j0)) .eq. Level .and.
     &                iLv(iNb(iNb0,j0)) .eq. Level ) then 
c
                 rhoC  = hvar(1,idcell)
                 preC  = hvar(6,idcell)
                 entC  = preC  / rhoC**gamma
                 rhoC0 = hvar(1,iNb0)
                 preC0 = hvar(6,iNb0)
                 entC0 = preC0 / rhoC0**gamma
                 rhoC1 = hvar(1,iNb1)
                 preC1 = hvar(6,iNb1)
                 entC1 = preC1 / rhoC1**gamma

              dp = abs ( preC0 - preC1 ) / min ( preC0, preC1 )
              dv = abs(hvar(jv,iNb0) / rhoC0 - hvar(jv,iNb1) / rhoC1)
c              dp = abs ( preC - preC1 ) / min ( preC, preC1 )
c              dv = abs(hvar(jv,idcell) / rhoC - hvar(jv,iNb1) / rhoC1)
c
c....         shock indicator
c
              if ( (dp .gt. epsP) .and. 
     &             (dv .gt. 0.  ) .and. 
     &             (rhoC .gt. floor_ref(3)) .and. 
     &             lrefok(3) ) 
     &        then
                Sind = w_ref(3) 
              else
                Sind = 0.
              end if
              

c....         contact discontinuity indicator
c             
              drho = abs ( rhoC0 - rhoC1 ) / min ( rhoC0, rhoC1 )
c              drho = abs ( rhoC - rhoC1 ) / min ( rhoC, rhoC1 )
              if ( (drho .gt. epsR) .and. 
     &             (dp   .lt. epsP) .and.
     &             (rhoC .gt. floor_ref(4)) .and. lrefok(4) ) 
     &        then
                Cind = w_ref(4) 
              else
                Cind = 0.
              end if
c
c....        density gradient indicator
c
             Rhoind0 = abs ( rhoC - rhoC0 ) /
     &                 min ( rhoC , rhoC0 )
             Rhoind1 = abs ( rhoC - rhoC1 ) /
     &                 min ( rhoC , rhoC1 )
             drho = max ( drho , Rhoind0 , Rhoind1 )
             if ( drho .gt. epsdrho .and. 
     &            rhoC .gt. floor_ref(5) .and. lrefok(5) ) then
               Rhoind  = w_ref(5) 
             else
               Rhoind  = 0.0
             end if
c
c....         pressure gradient indicator 
c
              Preind0 = abs ( preC - preC0 ) / ( preC + preC0 ) 
              Preind1 = abs ( preC - preC1 ) / ( preC + preC1 ) 
              if ( rhoC .gt. floor_ref(6) .and. lrefok(6) ) then
                Preind  = w_ref(6) * max ( Preind0 , Preind1 )
              else
                Preind  = 0.
              end if
c
c....         entropy gradient indicator 
c 
              Entind0 = abs ( entC - entC0 ) / ( entC + entC0 ) 
              Entind1 = abs ( entC - entC1 ) / ( entC + entC1 ) 
              dS = abs( ( entC0 - entC1 ) / min ( entC0, entC1 ))
              if ( rhoC .gt. floor_ref(7) .and. dS .gt. epsS .and.
     &             lrefok(7) ) then 
                Entind  = w_ref(7) * max ( Entind0 , Entind1 )
              else
                Entind  = 0.
              end if
c
              vnw(1,idcell) = max ( vnw(1,idcell), Sind, Cind,  
     &                              Rhoind, Preind, Entind      )
              endif
              endif
            enddo ! end ic2
          enddo ! end ic1
        enddo ! end j 
        ENDIF 

      ENDIF

      return
      end
c
c     ---------------------------------------
      subroutine Diffuse_Indicators ( Level ) 
c     ---------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level
c
      parameter ( react = 0.1 )
      parameter ( cdf = 0.15 ) 

      integer Neighbors(6)
c
c.... Reaction 
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(xx)
        do ic1 = 1 , ncell0 
          xx = vnw(1,ic1)
          if ( xx .ge. wsplit ) xx = min ( 1.0 , xx + react )
          vnw(1,ic1) = xx
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE(icell,ic2,idcell,xx)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            xx = vnw(1,idcell)
            if ( xx .ge. wsplit ) xx = min ( 1.0 , xx + react )
            vnw(1,idcell) = xx
          enddo
        enddo
      ENDIF
c
c.... Diffusion
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED) 
        do ic1 = 1 , ncell0 
          vnw(2,ic1) = vnw(1,ic1)
        enddo
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE(ic2,Neighbors,sum,vcell)
        do ic1 = 1 , ncell0 
          call NbAll0 ( ic1 , Neighbors ) 
          sum = 0. 
          vcell = vnw(1,ic1) 
          do ic2 = 1 , nneigh 
            sum = sum + (vnw(1,Neighbors(ic2)) - vcell)
          enddo
          vnw(2,ic1) = vnw(2,ic1) + cdf*sum  
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE(icell,ic2,idcell)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            vnw(2,idcell) = vnw(1,idcell)            
          enddo
        enddo
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE(icell,ic2,idcell,ic3,sum,vcell,Neighbors)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            call NbAllL ( idcell , Neighbors ) 
            sum = 0.
            vcell = vnw(1,idcell)
            do ic3 = 1 , nneigh              
              if ( iLv(Neighbors(ic3)) .eq. Level ) then 
                sum = sum + (vnw(1,Neighbors(ic3)) - vcell)
              endif
            enddo
            vnw(2,idcell) = vnw(2,idcell) + cdf*sum
          enddo
        enddo
      ENDIF
c
c.... Map back
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED) 
        do ic1 = 1 , ncell0 
          vnw(1,ic1) = min(1.0,vnw(2,ic1))
          vnw(2,ic1) = 0.
        enddo
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE(icell,ic2,idcell)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idcell = icell + ic2
            vnw(1,idcell) = min(1.0,vnw(2,idcell))
            vnw(2,idcell) = 0.
          enddo
        enddo
      ENDIF
c
      return
      end
c
c     ----------------------------------------
      subroutine CheckRefinementMask ( Level ) 
c     ----------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level
      real xcell, ycell, zcell
c
      IF ( Level .eq. MinLevel ) THEN 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(xcell,ycell,zcell)
        do ic1 = 1 , ncell0 
          call Ps ( ic1 , xcell, ycell, zcell )
          if ( (xcell .lt. xrefmin) .or. (xcell .gt. xrefmax) .or.
     &         (ycell .lt. yrefmin) .or. (ycell .gt. yrefmax) .or.
     &         (zcell .lt. zrefmin) .or. (zcell .gt. zrefmax) ) then
             vnw(1,ic1) = 0.
          endif
        enddo
      ELSE
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell, ic2, idc, xcell,ycell,zcell)
        do ic0 = 1 , nLevel
          icell = iSelect(ic0)
          do ic2 = 0 , 7
            idc = icell + ic2
            call Ps ( idc , xcell, ycell, zcell )
            if ( (xcell .lt. xrefmin) .or. (xcell .gt. xrefmax) .or.
     &           (ycell .lt. yrefmin) .or. (ycell .gt. yrefmax) .or.
     &           (zcell .lt. zrefmin) .or. (zcell .gt. zrefmax) ) then
              vnw(1,idc) = 0.
            endif
          enddo
        enddo

      ENDIF
c
      return
      end
c
c     ---------------------------------
      subroutine Split ( Level , mtot )
c     ---------------------------------
c
c     purpose: splits cells marked to split
c
c     input  : Level - level to process
c     output : mtot  - # of cells just split
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      integer mtot, Level
      integer idcell 
      real*8 e_kin, e_ip
      real*8 whvar(nhvar), wvar1, wvar2, wvar3 
      dimension iPyr(nchild,3)   ! interpolation pyramid vertices 
c
      data iPyr / 1, 2, 1, 2, 1, 2, 1, 2,   
     &            3, 3, 4, 4, 3, 3, 4, 4,
     &            5, 5, 5, 5, 6, 6, 6, 6  /
c
c.... Warning! The loops below are to be executed SERIALLY
c
      IF ( Level .eq. MinLevel ) THEN
        do ic1 = 1 , ncell0
          if ( vnw(1,ic1) .gt. wsplit ) then
            ires = iSplitCell ( ic1 ) 
            if ( ires .eq. nil ) then 
              mtot = mtot + 1
              iOC = iOctCh(ic1)     
              v_p   = hvar(3,ic1)**2 +
     &                hvar(4,ic1)**2 +
     &                hvar(5,ic1)**2        
              e_kin = 0.5d0 * v_p / hvar(1,ic1)
              v_p = sqrt(v_p) 
              e_ip  = hvar(2,ic1) - e_kin 
              wvar1 = 0.0
              wvar2 = 0.0
              wvar3 = 0.0
              do j = 1 , nhvar
                whvar(j) = 0.0
              enddo
              iChild1 = ( iOC - 1 ) * nchild + ncell0
              iChild8 = iChild1 + 8
              do j = 1 , nchild 
                iChild = iChild1 + j 
c
                iNb1   = iNb(ic1,iPyr(j,1))
                iNb2   = iNb(ic1,iPyr(j,2))
                iNb3   = iNb(ic1,iPyr(j,3))

                e_kin = 0.5d0*( hvar(3,iNb1)**2 +
     &                          hvar(4,iNb1)**2 + 
     &                          hvar(5,iNb1)**2 ) / hvar(1,iNb1)
                e_inb1  = hvar(2,iNb1) - e_kin 
                e_kin = 0.5d0*( hvar(3,iNb2)**2 +
     &                          hvar(4,iNb2)**2 + 
     &                          hvar(5,iNb2)**2 ) / hvar(1,iNb2)
                e_inb2  = hvar(2,iNb2) - e_kin 
                e_kin = 0.5d0*( hvar(3,iNb3)**2 +
     &                          hvar(4,iNb3)**2 + 
     &                          hvar(5,iNb3)**2 ) / hvar(1,iNb3)
                e_inb3  = hvar(2,iNb3) - e_kin 
                e_ip = 0.25 * (e_ip + e_inb1 + e_inb2 + e_inb3)
                var(1,iChild) = 0.25 * ( var(1,ic1) + 
     &                                   var(1,iNb1) +
     &                                   var(1,iNb2) +         
     &                                   var(1,iNb3) )
                wvar1 = wvar1 + var(1,iChild)
                var(2,iChild) = 0.25 * ( var(2,ic1) + 
     &                                   var(2,iNb1) +
     &                                   var(2,iNb2) +         
     &                                   var(2,iNb3) )
                wvar2 = wvar2 + var(2,iChild)
                var(3,iChild) = 0.25 * ( var(3,ic1) + 
     &                                   var(3,iNb1) +
     &                                   var(3,iNb2) +         
     &                                   var(3,iNb3) )
                wvar3 = wvar3 + var(3,iChild)
                hvar(1,iChild) = 0.25 * ( hvar(1,ic1) + 
     &                                    hvar(1,iNb1) +
     &                                    hvar(1,iNb2) +         
     &                                    hvar(1,iNb3) )
                whvar(1) = whvar(1) + hvar(1,iChild)
                hvar(3,iChild) = 0.25 * ( hvar(3,ic1) + 
     &                                       hvar(3,iNb1) +
     &                                       hvar(3,iNb2) +         
     &                                       hvar(3,iNb3) )
                hvar(4,iChild) = 0.25 * ( hvar(4,ic1) +
     &                                       hvar(4,iNb1) +
     &                                       hvar(4,iNb2) +
     &                                       hvar(4,iNb3) )
                hvar(5,iChild) = 0.25 * ( hvar(5,ic1) +
     &                                       hvar(5,iNb1) +
     &                                       hvar(5,iNb2) +
     &                                       hvar(5,iNb3) )
                hvar(6,iChild) = 0.25 * ( hvar(6,ic1) +
     &                                       hvar(6,iNb1) +
     &                                       hvar(6,iNb2) +
     &                                       hvar(6,iNb3) )
                hvar(8,iChild) = 0.25 * ( hvar(8,ic1) +
     &                                       hvar(8,iNb1) +
     &                                       hvar(8,iNb2) +
     &                                       hvar(8,iNb3) )
                hvar(7,iChild) = hvar(7,ic1)
                v_ch = hvar(3,iChild)**2 +
     &                 hvar(4,iChild)**2 +
     &                 hvar(5,iChild)**2
                e_kin = 0.5d0 * v_ch / hvar(1,iChild) 
                whvar(3) = whvar(3) + sqrt(v_ch) 
                hvar(2,iChild) = e_ip + e_kin
                whvar(2) = whvar(2) + hvar(2,iChild)
                whvar(6) = whvar(6) + hvar(6,iChild)
                whvar(8) = whvar(8) + hvar(8,iChild)
#ifdef ADVECT_SPECIES
                do jj = 9 , nhvar
                  hvar(jj,iChild) = 0.25 * ( hvar(jj,ic1)  +
     &                                       hvar(jj,iNb1) +
     &                                       hvar(jj,iNb2) +
     &                                       hvar(jj,iNb3) )
                  whvar(jj) = whvar(jj) + hvar(jj,iChild)
                enddo
#endif
              enddo
c....         enforce conservation
              wvar1 = 8.0d0 * var(1,ic1) / wvar1
              wvar2 = 8.0d0 * var(2,ic1) / wvar2
              wvar3 = 8.0d0 * var(3,ic1) / wvar3
              whvar(1) = 8.0d0 * hvar(1,ic1) / whvar(1)
              whvar(2) = 8.0d0 * hvar(2,ic1) / whvar(2)
              whvar(3) = 8.0d0 * v_p / whvar(3)
              whvar(6) = 8.0d0 * hvar(6,ic1) / whvar(6)
              whvar(8) = 8.0d0 * hvar(8,ic1) / whvar(8)
#ifdef ADVECT_SPECIES
              do jj = 9 , nhvar
                whvar(jj) = 8.0d0 * hvar(jj,ic1) / whvar(jj)
              enddo
#endif

              do j = iChild1+1 , iChild8
                var(1,j) = wvar1 * var(1,j)
                var(2,j) = wvar2 * var(2,j)
                var(3,j) = wvar3 * var(3,j)
                hvar(1,j) = whvar(1) * hvar(1,j)
                hvar(2,j) = whvar(2) * hvar(2,j)
                hvar(3,j) = whvar(3) * hvar(3,j)
                hvar(4,j) = whvar(3) * hvar(4,j)
                hvar(5,j) = whvar(3) * hvar(5,j)
                hvar(6,j) = whvar(6) * hvar(6,j)
                hvar(8,j) = whvar(8) * hvar(8,j)
#ifdef ADVECT_SPECIES
                do jj = 9 , nhvar
                  hvar(jj,j) = whvar(jj) * hvar(jj,j)
                enddo
#endif
              enddo
            else
              if ( ires .eq. -2 ) then 
                call Open_ASCII_File ( iOErr , errorfile , sequent )
                write(iOErr,*) '1  *error in Split: ran out of cells'
                close ( iOErr )
                stop
              endif               
            endif
          endif
        enddo          
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idc = icell + ic2
            if ( vnw(1,idc) .gt. wsplit ) then
              ires = iSplitCell ( idc )
              if ( ires .eq. nil ) then 
                mtot = mtot + 1
                iOC = iOctCh(idc)  
                v_p   = hvar(3,idc)**2 +
     &                  hvar(4,idc)**2 +
     &                  hvar(5,idc)**2
                e_kin = 0.5d0 * v_p / hvar(1,idc)
                v_p = sqrt(v_p)
                e_ip  = hvar(2,idc) - e_kin 
                wvar1 = 0.0
                wvar2 = 0.0
                wvar3 = 0.0
                do j = 1 , nhvar
                  whvar(j) = 0.0
                enddo

                iChild1 = ( iOC - 1 ) * nchild + ncell0 
                iChild8 = iChild1 + 8 
                do j = 1 , nchild
                  iChild = iChild1 + j  
                  iNb1   = iNb(idc,iPyr(j,1))
                  iNb2   = iNb(idc,iPyr(j,2))
                  iNb3   = iNb(idc,iPyr(j,3))
c                  iLv1 = iLv(iNb1)
c                  iLv2 = iLv(iNb2)
c                  iLv3 = iLv(iNb3)
                  e_kin = 0.5d0*( hvar(3,iNb1)**2 +
     &                            hvar(4,iNb1)**2 + 
     &                            hvar(5,iNb1)**2 ) / hvar(1,iNb1)
                  e_inb1  = hvar(2,iNb1) - e_kin 
                  e_kin = 0.5d0*( hvar(3,iNb2)**2 +
     &                            hvar(4,iNb2)**2 + 
     &                            hvar(5,iNb2)**2 ) / hvar(1,iNb2)
                  e_inb2  = hvar(2,iNb2) - e_kin 
                  e_kin = 0.5d0*( hvar(3,iNb3)**2 +
     &                            hvar(4,iNb3)**2 + 
     &                            hvar(5,iNb3)**2 ) / hvar(1,iNb3)
                  e_inb3  = hvar(2,iNb3) - e_kin 
                  e_ip = 0.25 * (e_ip + e_inb1 + e_inb2 + e_inb3)

                  var(1,iChild) = 0.25 * ( var(1,idc) + 
     &                                     var(1,iNb1) + 
     &                                     var(1,iNb2) + 
     &                                     var(1,iNb3) ) 
c                  wvar1 = wvar1 + var(1,iChild)
                  var(2,iChild) = 0.25 * ( var(2,idc) + 
     &                                     var(2,iNb1) +
     &                                     var(2,iNb2) +         
     &                                     var(2,iNb3) )
c                  wvar2 = wvar2 + var(2,iChild)
                  var(3,iChild) = 0.25 * ( var(3,idc) + 
     &                                     var(3,iNb1) +
     &                                     var(3,iNb2) +         
     &                                     var(3,iNb3) )
c                  wvar3 = wvar3 + var(3,iChild) 
                  hvar(1,iChild) = 0.25 * ( hvar(1,idc) + 
     &                                      hvar(1,iNb1) +
     &                                      hvar(1,iNb2) +         
     &                                      hvar(1,iNb3) )
                  whvar(1) = whvar(1) + hvar(1,iChild) 
                  hvar(3,iChild) = 0.25 * ( hvar(3,idc) +
     &                                      hvar(3,iNb1) +
     &                                      hvar(3,iNb2) +
     &                                      hvar(3,iNb3) )
                  hvar(4,iChild) = 0.25 * ( hvar(4,idc) +
     &                                      hvar(4,iNb1) +
     &                                      hvar(4,iNb2) +
     &                                      hvar(4,iNb3) )
                  hvar(5,iChild) = 0.25 * ( hvar(5,idc) +
     &                                      hvar(5,iNb1) +
     &                                      hvar(5,iNb2) +
     &                                      hvar(5,iNb3) )
                  hvar(6,iChild) = 0.25 * ( hvar(6,idc) +
     &                                      hvar(6,iNb1) +
     &                                      hvar(6,iNb2) +
     &                                      hvar(6,iNb3) )
                  hvar(8,iChild) = 0.25 * ( hvar(8,idc) +
     &                                      hvar(8,iNb1) +
     &                                      hvar(8,iNb2) +
     &                                      hvar(8,iNb3) )

                  hvar(7,iChild) = hvar(7,idc)
                  v_ch = hvar(3,iChild)**2 +
     &                   hvar(4,iChild)**2 +
     &                   hvar(5,iChild)**2
                  e_kin = 0.5d0 * v_ch / hvar(1,iChild)
                  whvar(3) = whvar(3) + sqrt(v_ch)
                  hvar(2,iChild) = e_ip + e_kin
                  whvar(2) = whvar(2) + hvar(2,iChild)
                  whvar(6) = whvar(6) + hvar(6,iChild)
                  whvar(8) = whvar(8) + hvar(8,iChild)
#ifdef ADVECT_SPECIES
               do jj = 9 , nhvar
                    hvar(jj,iChild) = 0.25 * ( hvar(jj,idc)  +
     &                                         hvar(jj,iNb1) +
     &                                         hvar(jj,iNb2) +
     &                                         hvar(jj,iNb3) )
                    whvar(jj) = whvar(jj) + hvar(jj,iChild)
 44              enddo
#endif
                enddo
c....           enforce conservation
c                wvar1 = 8.0d0 * var(1,idc) / wvar1
c                wvar2 = 8.0d0 * var(2,idc) / wvar2
c                wvar3 = 8.0d0 * var(3,idc) / wvar3
                whvar(1) = 8.0d0 * hvar(1,idc) / whvar(1)
                whvar(2) = 8.0d0 * hvar(2,idc) / whvar(2)
                whvar(3) = 8.0d0 * v_p / whvar(3)
                whvar(6) = 8.0d0 * hvar(6,idc) / whvar(6)
                whvar(8) = 8.0d0 * hvar(8,idc) / whvar(8)
                  
#ifdef ADVECT_SPECIES
                do jj = 9 , nhvar
                  if ( whvar(jj) .gt. 0.0 ) then 
                    whvar(jj) = 8.0d0 * hvar(jj,idc) / whvar(jj)
                  endif
                enddo
#endif
                do j = iChild1+1 , iChild8
c                  var(1,j) = wvar1 * var(1,j) 
c                  var(2,j) = wvar2 * var(2,j)
c                  var(3,j) = wvar3 * var(3,j)
                  hvar(1,j) = whvar(1) * hvar(1,j)
                  hvar(2,j) = whvar(2) * hvar(2,j)
                  hvar(3,j) = whvar(3) * hvar(3,j)
                  hvar(4,j) = whvar(3) * hvar(4,j)
                  hvar(5,j) = whvar(3) * hvar(5,j)
                  hvar(6,j) = whvar(6) * hvar(6,j)
                  hvar(8,j) = whvar(8) * hvar(8,j)
#ifdef ADVECT_SPECIES
                  do jj = 9 , nhvar
                    if ( whvar(jj) .gt. 0.0 ) then 
                      hvar(jj,j) = whvar(jj) * hvar(jj,j)
                    endif
                  enddo
#endif
                enddo

              else
                if ( ires .eq. -2 ) then 
                  call Open_ASCII_File ( iOErr , errorfile , sequent )
                  write(iOErr,*) '1  *error in Split: ran out of cells'
                  close ( iOErr )
                  stop
                endif               
              endif
            endif            
          enddo ! end ic2
        enddo ! end ic1
      ENDIF
c     
      return
      end
c
c     ------------------------------------
      subroutine DeRefine ( Level , mtot )
c     ------------------------------------
c
c     purpose: Makes derefinement on Level
c     input  : Level - level to process
c     output : mtot  - total number of refined/derefined cells
c     uses   : Join
c     vnw(*,2) - is used to mark cells for splitting & joining
c
      integer mtot         ! total # of cell marked to split 
      integer Level 
c
c.... mark cells for joining 
c
      call Mark_Cells_to_Join ( Level , mtot ) 

      call Join ( Level , mtot )

      return
      end
c
c     ----------------------------------------------
      subroutine Mark_Cells_to_Join ( Level , mtot )
c     ----------------------------------------------
c
c     purpose: marks cells for joining
c
c     input  : Level - level to process
c
c     warning: designed to work only in 3D
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer mtot , Level 
      logical unrefine
      integer idummy, iNbs(nneigh)
c
c.... mark cells according to their refinement indicator
c
c     exclude : cells that have split children
c               cells that have children with split neighbors
c               cells that do have both neighbors split and 
c                     not marked for joining 
c
      mtot = nil
c
      IF ( Level .eq. MinLevel ) THEN
        do j = 1 , ncell0
          iOC = iOctCh(j)
          if ( iOC .gt. nil ) then         ! if split
            unrefine = (vnw(1,j).lt.wjoin)
            if ( unrefine ) then      ! if join criterion is satisfied
              iCh0 = ( iOC - 1 ) * nchild + ncell0
              do ic1 = 1 , nchild
                iChild = iCh0 + ic1
c....           do not mark if a child is split
                if ( iOctCh(iChild) .gt. nil) then 
                  unrefine = .false.                  
                endif
c....           do not mark if a child's neighbor is split                
                do ic2 = 1 , nneigh
                  if ( iOctCh(iNb(iChild,ic2)) .gt. nil ) then 
                    unrefine = .false. 
                  endif
                enddo
              enddo
              if ( unrefine ) then 
                vnw(1,j) = -1.0
                mtot = mtot + 1
              endif         
            endif
          endif
        enddo ! end j 
c
        do j = 1 , ncell0
          if ( vnw(1,j) .lt. -0.5 ) then 
            call NbAll0 ( j , iNbs )
            unrefine = .true. 
            if ( (iOctCh(iNbs(1)) .gt. nil ) .and. 
     &           (iOctCh(iNbs(2)) .gt. nil ) .and.
     &           ( vnw(1,iNbs(1)) .gt.-0.5 ) .and.
     &           ( vnw(1,iNbs(2)) .gt.-0.5 )        ) then 
              unrefine = .false. 
            endif
            if ( unrefine ) then 
              if ( (iOctCh(iNbs(3)) .gt. nil ) .and. 
     &             (iOctCh(iNbs(4)) .gt. nil ) .and.
     &             ( vnw(1,iNbs(3)) .gt.-0.5 ) .and.
     &             ( vnw(1,iNbs(4)) .gt.-0.5 )        ) then 
                unrefine = .false. 
              endif
            endif
            if ( unrefine ) then 
              if ( (iOctCh(iNbs(5)) .gt. nil ) .and. 
     &             (iOctCh(iNbs(6)) .gt. nil ) .and.
     &             ( vnw(1,iNbs(5)) .gt.-0.5 ) .and.
     &             ( vnw(1,iNbs(6)) .gt.-0.5 )        ) then 
                unrefine = .false. 
              endif
            endif
            if ( .not. unrefine ) then
              vnw(1,j) = 0.
              mtot = mtot - 1
            endif
          endif
        enddo ! end j 
      ELSE
        nLevel = iNOLL(Level)
        call Select_Cells ( Level , nLevel ) 
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idc = icell + ic2
            iOC = iOctCh(idc)
            if ( iOC .gt. nil ) then         ! if split
              unrefine = (vnw(1,idc).lt.wjoin)

              if ( unrefine ) then    ! if join criterion is satisfied
                iCh0 = ( iOC - 1 ) * nchild + ncell0
                unrefine = .true. 
                do ic3 = 1 , nchild
                  iChild = iCh0 + ic3
c....             do not mark if a child is split
                  if ( iOctCh(iChild) .gt. nil) then 
                    unrefine = .false.                  
                  endif
c....             do not mark if a child's neighbor is split
                  call NbAllL(iChild,iNbs)
                  do ic4 = 1 , nneigh
                    if ( iOctCh(iNbs(ic4)) .gt. nil ) then 
                      unrefine = .false. 
                    endif
                  enddo
                enddo ! end ic3
                if ( unrefine ) then 
                  vnw(1,idc) = -1.0
                  mtot = mtot + 1
                endif         
              endif
            endif
          enddo ! end ic2
        enddo ! end ic1
c
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idc = icell + ic2
            if ( vnw(1,idc) .lt. -0.5 ) then 
              call NbAllL ( idc , iNbs )
              unrefine = .true. 
              if ( iLv(iNbs(1)) .eq. Level 
     &                        .and.
     &             iLv(iNbs(2)) .eq. Level ) then 
                if ( (iOctCh(iNbs(1)) .gt. nil ) .and. 
     &               (iOctCh(iNbs(2)) .gt. nil ) .and.
     &               ( vnw(1,iNbs(1)) .gt.-0.5 ) .and.
     &               ( vnw(1,iNbs(2)) .gt.-0.5 )        ) then 
                  unrefine = .false. 
                endif
              endif
              if ( iLv(iNbs(3)) .eq. Level 
     &                        .and.
     &             iLv(iNbs(4)) .eq. Level ) then 
                if ( unrefine ) then 
                  if ( (iOctCh(iNbs(3)) .gt. nil ) .and. 
     &                 (iOctCh(iNbs(4)) .gt. nil ) .and.
     &                 ( vnw(1,iNbs(3)) .gt.-0.5 ) .and.
     &                 ( vnw(1,iNbs(4)) .gt.-0.5 )        ) then 
                    unrefine = .false. 
                  endif
                endif
              endif
              if ( iLv(iNbs(5)) .eq. Level 
     &                        .and.
     &             iLv(iNbs(6)) .eq. Level ) then 
                if ( unrefine ) then 
                  if ( (iOctCh(iNbs(5)) .gt. nil ) .and. 
     &                 (iOctCh(iNbs(6)) .gt. nil ) .and.
     &                 ( vnw(1,iNbs(5)) .gt.-0.5 ) .and.
     &                 ( vnw(1,iNbs(6)) .gt.-0.5 )        ) then 
                    unrefine = .false. 
                  endif
                endif
              endif
              if ( .not. unrefine ) then 
                vnw(1,idc) = 0.0
                mtot = mtot - 1
              endif
            endif
          enddo ! end ic2
        enddo ! end ic1
c
      ENDIF
c
      return
      end
c
c     --------------------------------
      subroutine Join ( Level , mtot )
c     --------------------------------
c
c     purpose: splits cells marked to split
c
c     input  : Level - level to process
c     output : mtot  - total number of cells just split
c
      include 'a_tree.h'
      include 'a_control.h'
c
      integer Level
c
c.... join marked cells
c
      IF ( Level .eq. MinLevel ) THEN
        do j = 1 , ncell0
          if ( vnw(1,j) .lt. -0.5 ) then 
            ires = iJoinCell ( j )
            if ( ires .ne. 0 ) then
              write (*,'(''*Join: cannot join '',i10,2x,i2)') j, Level
              write (*,'(''*Join: ires        '',i10)') ires
              stop
            endif
            mtot = mtot + 1
          endif
        enddo
      ELSE
        nLevel = iNOLL(Level)      
        call Select_Cells ( Level , nLevel ) 
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do ic2 = 0 , 7
            idc = icell + ic2
            if ( vnw(1,idc) .lt. -0.5 ) then 
              ires = iJoinCell ( idc )
              if ( ires .ne. 0 ) then
                write (*,'(''*Join: cannot join '',i10,2x,i2)') idc, Level
                write (*,'(''*Join: ires        '',i10)') ires
                stop
              endif
              mtot = mtot + 1
            endif            
          enddo ! end ic2
        enddo ! end ic1
      ENDIF
c
      return
      end
c      
