c     =====================================================================
c                                                   
c            Particle routines                      
c                                                   
c     =====================================================================
      
c     ------------------
      subroutine Move ()
c     ------------------
c
c     purpose: advances particles on one largest time step (=astep)
c     uses   : Move_Level , Move_MinLevel
c
      include 'a_tree.h'
      include 'a_control.h'
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))

      Common /AuxTime/t1,t2,t3,t4,t5

      dimension iLStep(MinLevel:MaxLevel) , incStep(MinLevel:MaxLevel)
      double precision tpot, tden
      
      nmaxsteps = 2**iTimeBin(MaxLevelNow)
      tpot =0.
      tden =0.
      Max_Level = MaxLevelNow
      do Level = MinLevel , Max_Level
        incStep(Level) = 2**iTimeBin(MaxLevelNow) / 2**iTimeBin(Level) 
        iLStep (Level) = incStep(Level) 
        !iLStep (Level) = incStep(Level)/2+1 
        Aexp   (Level) = aexpn 
      enddo

      dummy =seconds()      
      IF(mpi_yes .NE. 1234567 ) n_all = lspecies(Nspecies) 

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1)
      do ic1 = 1 , n_all
        ipt(ic1) = 0. 
      enddo
      t1 =t1+seconds()

      call Timing ( 1 , -1 )
       call Assign_Density ( MinLevel , MaxLevelNow )
      call Timing ( 1 ,  1 )
      iProlongFlag = 1

      niter = niter1       ! more iterations at the start
      call Solve_Poisson  ( MinLevel , MaxLevelNow , iProlongFlag )    

      dummy =seconds()
      do ic1 = 1 , ncell0
         tpot = tpot + pot(ic1)
         tden = tden + var(1,ic1)
      enddo

      t1 =t1+seconds()

      MinModify =  1000
      MaxModify = -1000
      iProlongFlag = 0

      do nstep = 1 , nmaxsteps 
         If(nstep.eq.1)Then
            niter = niter1
         Else
            niter = niter2
         EndIf
        do Level = Max_Level , MinLevel , -1
!        do Level = MinLevel , Max_Level 
          if ( iLStep(Level) .eq. nstep ) then 
            iLStep(Level) = iLStep(Level) + incStep(Level)
            MinWorkLevel  = MinModify 
            MaxWorkLevel  = MaxModify  
            write(13,'(a,6i8)')' nstep/Level/Modify=',nstep,Level,
     &                      MinModify,MaxModify 
            call Timing ( 1 , -1 )
             call Assign_Density ( MinModify , MaxModify       )
            call Timing ( 1 ,  1 )

            call Solve_Poisson (MinWorkLevel,MaxWorkLevel,iProlongFlag)         
            call Timing ( 5 , -1 )
             call Move_Level      ( Level                        )
            call Timing ( 5 ,  1 )
            call Timing ( 7 , -1 )
              call LL_Update      ( Level , MinModify , MaxModify   )
            call Timing ( 7 ,  1 )

            call Timing ( 8 , -1 )
             call Modify ( mtot, MinModify, MaxModify, Max_Level, 0)
            call Timing ( 8 ,  1 )
          endif
        enddo
      enddo   

      CPU(4) = CPU(2) + CPU(3)

      return
      end

c     -------------------------------
      subroutine Move_Level ( Level )
c     -------------------------------
c
c     purpose: moves particles on a given Level
c     input  : Level  - level to move
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'

      parameter ( wbig   = 0.8 )
      parameter ( wsmall = 0.2 )

      integer   Level 
      dimension iPack(8)                     ! pack of cells to be used in 
                                             ! acceleration interpolation
      dimension g_x(8)   , g_y(8) , g_z(8)
      dimension f_p(3)   , f_n(3)            ! potential in p&n directions
      dimension iODir(6)                     ! opposite directions
      dimension iIP(8,3)                     ! set of pointers to the leftmost
                                             ! cells - to be used as kernels
      double precision xx , yy , zz , vvx , vvy , vvz
      double precision diff_x,diff_y,diff_z,
     &                          corr_x,corr_y,corr_z,
     &                          xdum,ydum,zdum
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /

      data iODir / 2 , 1 , 4 , 3 , 6 , 5 / 
     
 
      MinMove =  10000
      MaxMove = -10000
      iFlag = 0

      if ( Level .eq. MinLevel ) then 
        call Move_Min_Level ()
        return
      endif

      Size      = CellSize(Level)
      a_level   = Aexp    (Level)
      da_level  = dAexp   (Level)
      a_next    = a_level + da_level
      a_limit   = a_next - 0.5 * damin
      a_rel     = a_next - aexpn

      pconst = - sqrt( a_level / (Om0+Oml0*a_level**3) )*0.5
      ahalf  = a_level + da_level / 2.0
      xconst = 1.0 / sqrt (ahalf * (Om0+Oml0*ahalf**3)) / ahalf

      call Select_Cells ( Level , nLevel )
      write(80+Level,135) '--- Level=',Level,nLevel,Size,a_level,
     &      da_level,pconst,xconst
 135  format(a,i3,i9,5g12.4)
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1,ic2,icell,idcell,iParent,iHead,idummy,iap,ap,dap)
C$OMP+PRIVATE (xx,yy,zz,vvx,vvy,vvz,gx,gy,gz,pdum,Level1,Size,iPar)
C$OMP+PRIVATE (ilc,icn,iCounter,iDim,idcell1)
C$OMP+PRIVATE (diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE (xdum,ydum,zdum,d1,d2,d3,t1,t2,t3)
C$OMP+PRIVATE (t2t1,t2d1,d2t1,d2d1,inb2,inb4,inb24,iPack)
C$OMP+PRIVATE ( ic3,iWorkLevel,iDir,iNb_p,iNb_n,f_p,f_n,g_x,g_y,g_z)
C$OMP+PRIVATE (pconst1,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1)
C$OMP+PRIVATE (d3t2t1,d3t2d1,d3d2t1,d3d2d1,delta_x,idcold)
C$OMP+PRIVATE (Posx,Posy,Posz,ifl,jfl,kfl)
      do ic1 = 1 , nLevel
        icell = iSelect(ic1)
        do ic2 = 0 , nchild-1
          idcell  = icell + ic2

          iParent = iPr(idcell)
          iHead   = iCL(idcell)

          if ( iHead .gt. nil ) then
            idummy = iHead
            do while ( idummy .ne. nil )

              iap = ipt(idummy)        
              ap  = aexpn + damin * iap 
              dap = pt(idummy) 
              if ( ap .lt. a_limit ) then 
                xx  = x (idummy) 
                yy  = y (idummy)
                zz  = z (idummy) 
                vvx = vx(idummy)
                vvy = vy(idummy)
                vvz = vz(idummy)

                Level1 = Level
                idcell1= idcell
 2553           iPar   = idcell1
                ilc    = idcell1
                Size   = CellSize(Level1)
                call Ps ( iPar , Posx,Posy,Posz )
                if ( xx .eq. Posx ) xx = xx + dsmall
                if ( yy .eq. Posy ) yy = yy + dsmall
                if ( zz .eq. Posz ) zz = zz + dsmall
                ifl    = sign ( one ,  REAL(xx - Posx) )
                jfl    = sign ( one ,  REAL(yy - Posy) )
                kfl    = sign ( one , REAL(zz - Posz) )
                icn = (ifl + 1)/2 + (jfl + 1) + 2 * (kfl + 1) + 1
                iCounter = nil

                do iDim = 1 , nDim
                  if ( iIP(icn,iDim) .eq. nil ) goto 2555
                  ilc = iNb(ilc,iIP(icn,iDim))              
                enddo

                
 2555           call Ps ( ilc , Posx,Posy,Posz )
                diff_x = Posx - xx
                diff_y = Posy - yy
                diff_z = Posz - zz

                corr_x  = zero
                corr_y  = zero
                corr_z  = zero
        
                if ( abs(diff_x) .gt. nf67 ) then              
                  if ( diff_x .gt. 0. ) then 
                    corr_x = ng
                  else
                    corr_x = - ng
                  endif
                endif 
                if ( abs(diff_y) .gt. nf67 ) then              
                  if ( diff_y .gt. 0. ) then
                    corr_y = ng
                  else
                    corr_y = -ng
                  endif
                endif 
                if ( abs(diff_z) .gt. nf67 ) then              
                  if ( diff_z .gt. 0. ) then
                    corr_z = ng
                  else
                    corr_z = -ng
                  endif
                endif 
c....           correct coordinate to insure periodic bundaries
                xdum = xx + corr_x
                ydum = yy + corr_y
                zdum = zz + corr_z
c....           compute weights
                d1   = abs(xdum - Posx) / Size
                d2   = abs(ydum - Posy) / Size
                d3   = abs(zdum - Posz) / Size
                t1   = 1.0 - d1
                t2   = 1.0 - d2
                t3   = 1.0 - d3
                t2t1 = t2 * t1
                t2d1 = t2 * d1
                d2t1 = d2 * t1
                d2d1 = d2 * d1

                inb2  = iNb(ilc,2)
                inb4  = iNb(ilc,4)
                inb24 = iNb(inb2,4)

                iPack(1) = ilc
                iPack(2) = inb2
                iPack(3) = inb4
                iPack(4) = inb24
                iPack(5) = iNb(ilc,6)
                iPack(6) = iNb(inb2,6)
                iPack(7) = iNb(inb4,6)
                iPack(8) = iNb(inb24,6)

                do ic3 = 1 , nchild 
                  if ( iLv(iPack(ic3)) .ne. Level1 ) then 
                    idcell1 = iPr(idcell1)
                    Level1 = Level1 - 1
                    goto 2553
                  endif
                enddo

                iWorkLevel = iLv(ilc)
 
                do ic3 = 1 , nchild
                  iDir = -1
                  do iDim = 1 , ndim
                    iDir  = iDir + 2
                    iNb_p = iNb(iPack(ic3),iDir+1)   ! neighbor in positive dir
                    iNb_n = iNb(iPack(ic3),iDir  )
                    if (iLv(iNb_p) .lt. iWorkLevel ) then
                      idcell1 = iPr(idcell1)                           !!
                      Level1 = Level1 - 1                           !!
                      goto 2553                               !!
c                      f_p(iDim) = wbig   * pot(iNb_p) +   !!
c     &                            wsmall * pot(iPack(ic3))    !!
                    else
                      f_p(iDim) = pot(iNb_p)
                    endif
                    if (iLv(iNb_n) .lt. iWorkLevel ) then
                      idcell1 = iPr(idcell1)
                      Level1 = Level1 - 1
                      goto 2553
c                      f_n(iDim) = wbig   * pot(iNb_n) +  !!
c     &                            wsmall * pot(iPack(ic3))   !!
                    else
                      f_n(iDim) = pot(iNb_n)
                    endif
                  enddo
                  g_x(ic3) = f_p(1) - f_n(1)
                  g_y(ic3) = f_p(2) - f_n(2)
                  g_z(ic3) = f_p(3) - f_n(3)
                enddo
c
c....           interpolate accelerations to the particle location
c
                pconst1 = pconst * (( a_level - ap ) + 
     &                    (da_level + dap)*0.5 )/ CellSize(iWorkLevel)

                pt3     = pconst1 * t3
                pd3     = pconst1 * d3
                t3t2t1  = pt3 * t2t1
                t3t2d1  = pt3 * t2d1
                t3d2t1  = pt3 * d2t1
                t3d2d1  = pt3 * d2d1
                d3t2t1  = pd3 * t2t1
                d3t2d1  = pd3 * t2d1
                d3d2t1  = pd3 * d2t1
                d3d2d1  = pd3 * d2d1

                gx = t3t2t1*g_x(1) + t3t2d1*g_x(2) + 
     &               t3d2t1*g_x(3) + t3d2d1*g_x(4) +
     &               d3t2t1*g_x(5) + d3t2d1*g_x(6) + 
     &               d3d2t1*g_x(7) + d3d2d1*g_x(8)
 
                gy = t3t2t1*g_y(1) + t3t2d1*g_y(2) + 
     &               t3d2t1*g_y(3) + t3d2d1*g_y(4) +
     & 	             d3t2t1*g_y(5) + d3t2d1*g_y(6) + 
     &               d3d2t1*g_y(7) + d3d2d1*g_y(8)

                gz = t3t2t1*g_z(1) + t3t2d1*g_z(2) + 
     &               t3d2t1*g_z(3) + t3d2d1*g_z(4) +
     &	             d3t2t1*g_z(5) + d3t2d1*g_z(6) + 
     &               d3d2t1*g_z(7) + d3d2d1*g_z(8)
                
      pdummy(idummy) = ( t3t2t1 * pot(iPack(1)) +
     &                   t3t2d1 * pot(iPack(2)) + 
     &                   t3d2t1 * pot(iPack(3)) + 
     &                   t3d2d1 * pot(iPack(4)) +
     &                   d3t2t1 * pot(iPack(5)) + 
     &                   d3t2d1 * pot(iPack(6)) + 
     &                   d3d2t1 * pot(iPack(7)) + 
     &                   d3d2d1 * pot(iPack(8)) ) / pconst1
!      if(xx < 35. .and.  abs(zz-33) < 2.)
!     & write(80+Level,1100)xx,yy,zz,gz,gy,gz,pconst1,
!     &                  CellSize(iWorkLevel), 
!     &                  pdummy(idummy),Level,Level1,idummy,
!     &                  i_par(idummy)

! 1100 format(3f9.4,3g12.4,2g13.4,g13.5,2i4,2i10)

                vvx   = vvx + gx 
                vvy   = vvy + gy 
                vvz   = vvz + gz

                delta_x = xconst * (a_next - ap)

	        xx    = xx  + vvx * delta_x
	        yy    = yy  + vvy * delta_x
	        zz    = zz  + vvz * delta_x

                if ( xx .lt. 1. ) xx = xx + yn    ! imposing periodic boundary
                if ( xx .ge. xn ) xx = xx - yn    ! conditions
                if ( yy .lt. 1. ) yy = yy + yn
                if ( yy .ge. xn ) yy = yy - yn
                if ( zz .lt. 1. ) zz = zz + yn
                if ( zz .ge. xn ) zz = zz - yn
                if(xx .ge. xn ) xx = xn-1.d-5
                if(yy .ge. xn ) yy = xn-1.d-5
                if(zz .ge. xn ) zz = xn-1.d-5


                if ( xx .lt. 1. .or.  xx .ge. xn)Then
                   write(13,*) ' error in coords: new=',xx,yy,zz
                   write(13,*) '                : old=',
     &                       x(idummy),y(idummy),z(idummy)
                   write(13,*) '          particle=',idummy
                   write(13,*) '          pconst=',pconst,
     &               CellSize(iWorkLevel),iWorkLevel,da_level,dap
                   write(13,*) '          xconst=',xconst,pconst1,
     &                     ' da=',(a_next - ap)
                   write(13,*) '          veloc =',vvx,vvy,vvz
                   write(13,*) '          accel =',gx,gy,gz
                   iFlag = 1
                endif
                   
                 x(idummy) =  xx
                 y(idummy) =  yy
                 z(idummy) =  zz
                vx(idummy) = vvx
                vy(idummy) = vvy
                vz(idummy) = vvz

                ap         = a_rel
                ipt(idummy) = nint( ap / damin )
                pt(idummy) =  da_level 
              endif

              idummy     = iLL(idummy,1)

 1553       continue
           enddo  ! end do while
         endif
        enddo
      enddo

      Aexp (Level) = a_level + da_level 
c                If(iFlag.eq.1)Then
c                            close(13)
c        CALL MPI_Abort(MPI_COMM_WORLD,ier1,ier2)
c                   STOP
c                endif



      return
      end
      
c     ----------------------------
      subroutine Move_Min_Level ()
c     ----------------------------
c
c     purpose: moves particles on the root level
c     
c     
      include 'a_tree.h'
      include 'a_control.h'

      parameter ( wbig   = 0.8 )
      parameter ( wsmall = 0.2 )

      dimension iPack(8)                     ! pack of cells to be used in 
                                             ! acceleration interpolation
      dimension g_x(8), g_y(8), g_z(8)
      dimension iIP(8,3)                     ! set of pointers to the leftmost
                                             ! cells - to be used as kernels
      double precision xx , yy , zz , vvx , vvy , vvz
      double precision diff_x,diff_y,diff_z,
     &                          corr_x,corr_y,corr_z,
     &                          xdum,ydum,zdum

      data iIP / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &           3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &           5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /

      Size    = 1.0
      pconst1 = 1.0 / Size

      a_level   = Aexp    (MinLevel)
      da_level  = dAexp   (MinLevel)
      a_next    = a_level + da_level
      a_limit   = a_next - 0.5 * damin
      a_rel     = a_next - aexpn

      pconst = - sqrt( a_level / (Om0+Oml0*a_level**3)) * 0.5
      ahalf  = a_level + da_level / 2.0
      xconst = 1.0 / sqrt (ahalf * (Om0+Oml0*ahalf**3)) / ahalf

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1,idummy,iap,ap,dap,pdum)
C$OMP+PRIVATE(xx,yy,zz,vvx,vvy,vvz,gx,gy,gz)
C$OMP+PRIVATE(iPar,ilc,icn,iCounter,iDim)
C$OMP+PRIVATE(diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE(xdum,ydum,zdum,i,j,k,pos1,pos2,pos3)
C$OMP+PRIVATE(d1,d2,d3,t1,t2,t3,t2t1,t2d1,d2t1,d2d1)
C$OMP+PRIVATE(inb2,inb4,inb24,iPack,ic3,g_x,g_y,g_z)
C$OMP+PRIVATE(pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1)
C$OMP+PRIVATE(d3t2t1,d3t2d1,d3d2t1,d3d2d1,pconst1,delta_x)
C$OMP+PRIVATE(Posx,Posy,Posz,ifl,jfl,kfl)
      do ic1 = 1 , ncell0
        if ( iCL(ic1) .gt. nil ) then
          idummy = iCL(ic1)
          do while ( idummy .ne. nil )
            iap = ipt(idummy)
            ap  = aexpn + damin * iap 
            dap = pt(idummy)

            if ( ap .lt. a_limit ) then 
              xx  = x(idummy) 
              yy  = y(idummy)
              zz  = z(idummy) 
              vvx = vx(idummy)
              vvy = vy(idummy)
              vvz = vz(idummy)
            
              iPar = ic1
              ilc  = iPar
              call Ps ( iPar , Posx,Posy,Posz )
              if ( xx .eq. Posx ) xx = xx + dsmall
              if ( yy .eq. Posy ) yy = yy + dsmall
              if ( zz .eq. Posz ) zz = zz + dsmall
              ifl = sign ( one , REAL(xx - Posx) )
              jfl = sign ( one , REAL(yy - Posy) )
              kfl = sign ( one , REAL(zz - Posz) )
              icn = (ifl + 1)/2 + (jfl + 1) + 2 * (kfl + 1) + 1
              iCounter = nil

              do iDim = 1 , nDim
                if ( iIP(icn,iDim) .eq. nil ) goto 2556
                ilc = iNb(ilc,iIP(icn,iDim))             
              enddo

 2556         call Ps ( ilc , Posx,Posy,Posz )
              diff_x = Posx - xx
              diff_y = Posy - yy
              diff_z = Posz - zz

              corr_x = zero
              corr_y = zero
              corr_z = zero        
              if ( abs(diff_x) .gt. nf67 ) then              
                if ( diff_x .gt. 0. ) then
                  corr_x = ng
                else
                  corr_x = - ng
                endif
              endif 
              if ( abs(diff_y) .gt. nf67 ) then              
                if ( diff_y .gt. 0. ) then
                  corr_y = ng
                else
                  corr_y = -ng
                endif
              endif 
              if ( abs(diff_z) .gt. nf67 ) then              
                if ( diff_z .gt. 0. ) then
                  corr_z = ng
                else
                  corr_z = -ng
                endif
              endif 
c....         correcting position to account for periodic boundaries
              xdum = xx + corr_x
              ydum = yy + corr_y
              zdum = zz + corr_z

c....         compute position of lefmost corner cell (ilc)
              i = (ilc - 1)/ng2 + 1
              j = (ilc - (i-1)*ng2 - 1)/ng + 1
              k =  ilc - (i-1)*ng2 - (j-1)*ng

              pos1 = i + 0.5
              pos2 = j + 0.5
              pos3 = k + 0.5

c....         compute CIC coefficients

              d1   = abs(xdum - pos1)  ! if CellSize(MinLevel).ne.1.0 
              d2   = abs(ydum - pos2)  ! MUST divide each of these by 
              d3   = abs(zdum - pos3)  ! CellSize(MinLevel)
              t1   = 1.0 - d1
              t2   = 1.0 - d2
              t3   = 1.0 - d3
              t2t1 = t2 * t1
              t2d1 = t2 * d1
              d2t1 = d2 * t1
              d2d1 = d2 * d1
            
              inb2  = iNb(ilc,2)
              inb4  = iNb(ilc,4)
              inb24 = iNb(iNb(ilc,2),4)

              iPack(1) = ilc
              iPack(2) = inb2
              iPack(3) = inb4
              iPack(4) = inb24
              iPack(5) = iNb(ilc,6)
              iPack(6) = iNb(inb2,6)
              iPack(7) = iNb(inb4,6)
              iPack(8) = iNb(inb24,6)

              do ic3 = 1 , nchild
                g_x(ic3) = pot(iNb(iPack(ic3),2)) - 
     &                     pot(iNb(iPack(ic3),1))
                g_y(ic3) = pot(iNb(iPack(ic3),4)) - 
     &                     pot(iNb(iPack(ic3),3))
                g_z(ic3) = pot(iNb(iPack(ic3),6)) - 
     &                     pot(iNb(iPack(ic3),5))
              enddo

c....         interpolate accelerations to the particle location

              pconst1 = pconst * ((a_level - ap) + (da_level + dap)*0.5)
              pt3     =  t3 * pconst1
              pd3     =  d3 * pconst1
              t3t2t1  = pt3 * t2t1
              t3t2d1  = pt3 * t2d1
              t3d2t1  = pt3 * d2t1
              t3d2d1  = pt3 * d2d1
              d3t2t1  = pd3 * t2t1
              d3t2d1  = pd3 * t2d1
              d3d2t1  = pd3 * d2t1
              d3d2d1  = pd3 * d2d1

              gx = t3t2t1*g_x(1) + t3t2d1*g_x(2) + 
     &             t3d2t1*g_x(3) + t3d2d1*g_x(4) +
     &             d3t2t1*g_x(5) + d3t2d1*g_x(6) + 
     &             d3d2t1*g_x(7) + d3d2d1*g_x(8)

              gy = t3t2t1*g_y(1) + t3t2d1*g_y(2) + 
     &             t3d2t1*g_y(3) + t3d2d1*g_y(4) +
     &	           d3t2t1*g_y(5) + d3t2d1*g_y(6) + 
     &             d3d2t1*g_y(7) + d3d2d1*g_y(8)

              gz = t3t2t1*g_z(1) + t3t2d1*g_z(2) + 
     &             t3d2t1*g_z(3) + t3d2d1*g_z(4) +
     &	           d3t2t1*g_z(5) + d3t2d1*g_z(6) + 
     &             d3d2t1*g_z(7) + d3d2d1*g_z(8)

c.... get potential energy at the particle location
      pdummy(idummy) = ( t3t2t1 * pot(iPack(1)) +
     &                   t3t2d1 * pot(iPack(2)) + 
     &                   t3d2t1 * pot(iPack(3)) + 
     &                   t3d2d1 * pot(iPack(4)) +
     &                   d3t2t1 * pot(iPack(5)) + 
     &                   d3t2d1 * pot(iPack(6)) + 
     &                   d3d2t1 * pot(iPack(7)) + 
     &                   d3d2d1 * pot(iPack(8)) ) / pconst1

              delta_x = (a_next - ap) * xconst

              vvx   = vvx + gx
              vvy   = vvy + gy
              vvz   = vvz + gz
	      xx    = xx  + vvx * delta_x
	      yy    = yy  + vvy * delta_x
	      zz    = zz  + vvz * delta_x

              if ( xx .lt. 1. ) xx = xx + yn     ! imposing periodic boundary
              if ( xx .ge. xn ) xx = xx - yn     ! conditions
              if ( yy .lt. 1. ) yy = yy + yn
              if ( yy .ge. xn ) yy = yy - yn
              if ( zz .lt. 1. ) zz = zz + yn
              if ( zz .ge. xn ) zz = zz - yn

                if(xx .ge. xn ) xx = xn-1.d-7
                if(yy .ge. xn ) yy = xn-1.d-7
                if(zz .ge. xn ) zz = xn-1.d-7

              x(idummy)  = xx
              y(idummy)  = yy
              z(idummy)  = zz
              vx(idummy) = vvx
              vy(idummy) = vvy
              vz(idummy) = vvz
              ap         = a_rel
              ipt(idummy) = nint( ap / damin )
              pt(idummy) =  da_level
            endif
            idummy         = iLL(idummy,1)

          enddo  ! end do while
        endif
      enddo

      Aexp (MinLevel) = a_next

      return
      end

c     ------------------------------------------------------
      subroutine LL_Update ( Level , MinModify , MaxModify )
c     ------------------------------------------------------
c
c     purpose: to update particle linked list
c              after particle have been moved
c              (this is necessary mainly when
c               code is run in parallel)
c     input  : Level - level to update  
c     output : MinModify , MaxModify - range of affected levels
c
      include 'a_tree.h'
      include 'a_control.h'

      
      integer Level , MinModify , MaxModify
      double precision xx , yy , zz 

      MinModify = Level 
      MaxModify = Level 

      IF ( Level .eq. MinLevel ) THEN

        do icell = 1 , ncell0
          if ( iCL(icell) .gt. nil ) then
            idummy = iCL(icell)

            do while ( idummy .ne. nil )
              xx        = x(idummy) 
              yy        = y(idummy)
              zz        = z(idummy)
              Next_Cell = iFindCell ( MaxLevel , xx , yy , zz )
              idummy2   = iLL(idummy,1)
c              if(Next_Cell.lt.0)Then                      !!!!!!
c                 write (13,*) ' error !!!'
c              endif 

              if ( Next_Cell .ne. icell ) then
                Next_Level = iLv(Next_Cell)
                if ( Next_Level.lt.MinModify ) MinModify = Next_Level
                if ( Next_Level.gt.MaxModify ) MaxModify = Next_Level
                call LL_Delete ( icell     , idummy )
                call LL_Insert ( Next_Cell , idummy )
              endif
              idummy = idummy2
            enddo ! end do while
          endif
        enddo ! end icell

      ELSE

        nLevel = iNOLL(Level)
        icell  = iGet_Cell_Index ( iHOLL(Level) ) 
c        write(13,*) ' nLevel=',nLevel,' icell=',icell,Level
        do ic0 = 1 , nLevel
          icell = icell + nchild 
          do ic1 = 1 , nchild
            idcell = icell - ic1
            if ( iCL(idcell) .gt. nil ) then
              idummy = iCL(idcell)
              do while ( idummy .ne. nil )
                xx        = x(idummy)
                yy        = y(idummy)
                zz        = z(idummy)
c                If(idcell.gt.nctot.or.idcell.lt.0)
c     &           write(13,'(" idcell:",3i10,3g11.3)')idummy,idcell,
c     &                        MaxLevel , xx , yy , zz
c                If(xx.lt.1.00.or.xx.gt.513.0)
c     &           write(13,'(" xx:",3i10,3g11.3)')idummy,idcell,
c     &                        MaxLevel , xx , yy , zz
c                If(Level.eq.2)
c     &           write(13,'("   :",3i10,3g11.3)')idummy,idcell,
c     &                        MaxLevel , xx , yy , zz

                Next_Cell = iFindCell ( MaxLevel , xx , yy , zz )
                idummy2   = iLL(idummy,1)

                if ( Next_Cell .ne. idcell ) then
                  Next_Level = iLv(Next_Cell)
                  if ( Next_Level.lt.MinModify ) MinModify = Next_Level
                  if ( Next_Level.gt.MaxModify ) MaxModify = Next_Level
                  call LL_Delete ( idcell    , idummy )
                  call LL_Insert ( Next_Cell , idummy )
                endif
                idummy = idummy2
              enddo ! end do while
            endif
          enddo ! end do ic1
          iOct  = iGet_Oct_Number ( idcell     )
          icell = iGet_Cell_Index ( iOctLL1(iOct) )            
        enddo

      ENDIF

      return
      end

c     --------------------------
      subroutine LL_Construct ()
c     --------------------------
c
c     constructs  doubly linked list of particles
c     for the root mesh (only when the code starts)
c     see description of algorithms in  
c       T.H.Cormen , C.E.Leiserson , R.L.Rivest
c       "Introduction to algorithms" pp. 204-208
c
      include 'a_tree.h'
      include 'a_control.h'
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))

      Do i=1,mcell
        iNP(i) = 0                         ! initialize particle counter
      EndDo
      IF(mpi_yes .NE. 1234567 ) n_all = lspecies(Nspecies)
!       write(13,*) ' LL_con: nall=',n_all
      do ic1 = 1 , n_all

        ip    = int(x(ic1))
        jp    = int(y(ic1))
        kp    = int(z(ic1))
        icell = iCellIndex ( ip , jp , kp )
!        if(ip<1.or.ip>ng)write(13,'(a,2i8,g12.4)')' ip=',ip,ic1,x(ic1)
!        if(jp<1.or.jp>ng)write(13,'(a,2i8,g12.4)')' jp=',jp,ic1,y(ic1)
!        if(kp<1.or.kp>ng)write(13,'(a,2i8,g12.4)')' kp=',kp,ic1,z(ic1)
        call LL_Insert ( icell , ic1 )     ! add partilces
      enddo
      ip = 0
      Do i=1,ncell0
         ip =ip + iNP(i)
      EndDo
!      write (13,*) ' LL_Construct: Nparticles=',ip,n_all
      return
      end

c     --------------------------------------
      subroutine LL_Insert ( icell , ipart )
c     --------------------------------------
c
c     inserts particle ipart into linked list for cell icell
c
c     input: icell , ipart
c
c     iLL(icell,1) - pointer to the next list entry
c     iLL(icell,2) - pointer to the previous list entry
c
      include 'a_tree.h'
      
      iLL(ipart,1) = iCL(icell)
      if ( iCL(icell) .ne. nil ) then
        iLL(iCL(icell),2) = ipart
      endif
      iCL(icell)   = ipart
      iLL(ipart,2) = nil
      iNp(icell) = iNp(icell) + 1
      return
      end

c     --------------------------------------
      subroutine LL_Delete ( icell , ipart )
c     --------------------------------------
c
c     deletes particle ipart from linked list for cell icell
c     input: icell , ipart
c
      include 'a_tree.h'

      if ( iLL(ipart,2) .ne. nil ) then
        iLL(iLL(ipart,2),1) = iLL(ipart,1)
      else
        iCL(icell) = iLL(ipart,1)
      endif
      if ( iLL(ipart,1) .ne. nil ) then
        iLL(iLL(ipart,1),2) = iLL(ipart,2)
      endif
      iNp(icell) = iNp(icell) - 1
      return
      end

c     -----------------------------
      subroutine LL_Split ( icell )
c     -----------------------------
c
c     constructs linked lists for icell's children
c     cancels linked list for icell
c
      include 'a_tree.h'

      dimension iChildren(nchild)
      double precision xx , yy , zz


      idummy = iCL(icell)
      if ( idummy .eq. nil ) return

      do ic1 = 1 , nchild
        iChildren(ic1) = iCh(icell,ic1)
      enddo

      call Ps ( icell , Posx,Posy,Posz )
      do while ( idummy .ne. nil )  ! construct linked lists for children
        xx = x(idummy)
        yy = y(idummy)
        zz = z(idummy)
        if ( xx .eq. Posx ) xx = xx + dsmall
        if ( yy .eq. Posy ) yy = yy + dsmall
        if ( zz .eq. Posz ) zz = zz + dsmall
        ifl = sign ( one , REAL(xx - Posx) )
        jfl = sign ( one , REAL(yy - Posy) )
        kfl = sign ( one , REAL(zz - Posz) )
        iChild  = (ifl + 1)/2 + (jfl + 1) + 2 * (kfl + 1) + 1
        idummy2 = iLL(idummy,1)
        call LL_Insert ( iChildren(iChild) , idummy )
        idummy = idummy2
      enddo
      iCL(icell) = nil       ! cancel linked list for LL
      iNP(icell) = nil
      return
      end

c     ----------------------------
      subroutine LL_Join ( icell )
c     ----------------------------
c
c     constructs linked list for icell
c     joining linked lists of its former children
c     cancels linked lists for children
c
      include 'a_tree.h'
      include 'a_control.h'

      dimension iNonEmpty(nchild)  ! array for non-empty children
      integer ne 
c.... error if icell is a leaf
!      if ( iOctCh(icell) .eq. nil ) then
!        write(13,*) '*LL_Join error: cell is a leaf'
!        name = 'ERROR_RUN'
!        accesstype = 'append'
!        call Open_ASCII_FILE(19,name,accesstype)
!        rewind 19
!        write(19,*) 
!     +  ' 2  : ART_Particles.f,  *LL_Join error: cell is a leaf'
!        close(19) 
!        stop
!      endif
c.... find non-empty kids 
      ne = nil
      idummy = iCh(icell,1) - 1
      do ic1 = 1 , nchild
        idummy = idummy + 1
        if ( iCL(idummy) .gt. nil ) then
          ne = ne + 1
          iNonEmpty(ne) = idummy
        endif
      enddo
c.... re-arrange linked lists
      iNP(icell) =0                    ! initialize number of particles
      if ( ne .gt. nil ) then 
        iCL(icell) = iCL(iNonEmpty(1))
        do ic1 = 1 , ne - 1
          idummy1 = iCL(iNonEmpty(ic1))
          do while ( idummy1 .ne. nil ) 
            idummy2 = idummy1 
            idummy1 = iLL(idummy1,1)
          enddo
          iLL(idummy2,1) = iCL(iNonEmpty(ic1+1))
          iLL(iLL(idummy2,1),2) = idummy2
        enddo
        do ic1 = 1 , ne
          iCL(iNonEmpty(ic1)) = nil      ! empty header
          iNP(icell) =iNP(icell) + iNP( iNonEmpty(ic1)) ! add particles
          iNP( iNonEmpty(ic1)) = 0      ! empy particle count
        enddo
      endif
      return
      end
c
c     ------------------------------------
c      integer function iWhichSpecie ( ip )      not used, should be not used 
c     ------------------------------------       in the MPI version
c
c      include 'a_tree.h'
c
c      integer ip, isd, is
c
c      isd = 0
c      is  = 1000
c      do while ( is .gt. 0 )  
c        isd = isd + 1
c        is  = (ip-1) / lsp(isd)
c      enddo
c
c      iWhichSpecie = isd
c
c      return
c      end
