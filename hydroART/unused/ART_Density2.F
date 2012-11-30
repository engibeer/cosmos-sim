c
c     Density assignment routines   
c-------------------------------------------------------------------- 
c
c     based on ART_Density1.F
c
c     modifications by Daniel Ceverino (2005)
c     for galaxy initial conditions. ( hstart_Galaxy.F)
c                  Assignment of Gas properties from DM density field.
c
c-------------------------------------------------------------------- 
c     ------------------------------------------------------
      subroutine Assign_Density_AND_GAS( Level1 , Level2 , NOFlag,vcorr)
c     ------------------------------------------------------
c
c     purpose: assigns density to cells on levels Level1-Level2
c      Level1 <= Level2
c
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
        call Initialize_Density_AND_GAS ( Level ) 
      enddo
c
      call Assign_ParticleDensity_AND_GAS ( Level1 , Level2 )
      do Level = Level1 , Level2
         call Assign_GAS_Final (Level,  vcorr )
      enddo
#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append )
      write(iOWork,*) 'Density and gas assignment from level ', Level1,
     $ ' to level ', Level2 
      close(iOWork)
#endif
c
c      call Derefine_Gas ()
      return
      end
c
c     ---------------------------------------
      subroutine Initialize_Density_AND_GAS ( Level ) 
c     ---------------------------------------
c
      integer Level 
c      double precision Size3
      include 'a_tree.h'
      include 'a_control.h'

      IF ( Level .eq. MinLevel ) THEN
c
C$OMP PARALLEL DO DEFAULT(SHARED)
        do ic1 = 1 , ncell0
          var(1,ic1) = -1.0
          ref(ic1) = 0.

          hvar(1,ic1)=0.

          hvar(3,ic1)=0.
          hvar(4,ic1)=0.
          hvar(5,ic1)=0.

          hvar(2,ic1)=0.
          hvar(6,ic1)=0.
          hvar(7,ic1)=0.
          hvar(8,ic1)=0.

          hvar(izII,ic1) = 0.    
          hvar(izIa,ic1) = 0.              
        enddo
      ELSE
c
        Size3 = -CellSize(Level)**3
        call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell)
        do ic1 = 1 , nLevel
          icell = iSelect(ic1)
          do iadd=0,7
             var(1,icell+iadd) =  Size3 
             ref (icell+iadd) = 0.

             hvar(1, icell+iadd) = 0.

             hvar(3, icell+iadd) = 0.             
             hvar(4, icell+iadd) = 0.   
             hvar(5, icell+iadd) = 0.   

          hvar(2, icell+iadd)=0.
          hvar(6, icell+iadd)=0.
          hvar(7, icell+iadd)=0.
          hvar(8, icell+iadd)=0.

          hvar(izII, icell+iadd)=0.
          hvar(izIa, icell+iadd ) = 0.    
          enddo
       enddo
      ENDIF
      
      return
      end
c
c
c     -----------------------------------------------------
      subroutine Assign_ParticleDensity_AND_GAS ( Level1 , Level2 )
c     -----------------------------------------------------
c
c     purpose: assigns density to all existing cells     
c
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      include 'a_dnew.fh' 
      integer Level, Level1, Level2
      double precision xx, yy, zz, xdum, ydum, zdum
      real*8 t1, t2, t3, diff_x, diff_y, diff_z


      if ((Level1.lt.MinLevel).or.(Level1.gt.MaxLevel)) return          
      if ((Level2.lt.MinLevel).or.(Level2.gt.MaxLevel)) return          
      
      if ( Level1 .eq. MinLevel ) then
        Ldummy = Level1 + 1
        call Assign_ParticleDensity_AND_GAS0 ()
      else
        Ldummy = Level1
      endif

      DO Level = Ldummy , Level2 

        Sizei  = CellSizei(Level)
        hSize  = CellSize(Level) * 0.5 
        npch = 0
        DO ip = 1 , np
c
c....     preselect Np_chunk particles that can contribute density to this Level
c           
           if ( iPL(ip) .ge. Level ) then 
              npch = npch + 1             
              ipc(npch) = ip
           endif 
           if ( npch .eq. Np_chunk .or. ip .eq. np ) then
              write (*,*)  'Level= ', Level, 'ip= ', ip, 'npch= ', npch
              write (*,*)  '---------************************________'
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i)
              do i = 1 , npch 
                 dpc(i,1) = zero
                 dpc(i,2) = zero
                 dpc(i,3) = zero
                 dpc(i,4) = zero
                 dpc(i,5) = zero
                 dpc(i,6) = zero
                 dpc(i,7) = zero
                 dpc(i,8) = zero
                 iac(i,1) = nil
                 iac(i,2) = nil
                 iac(i,3) = nil
                 iac(i,4) = nil
                 iac(i,5) = nil
                 iac(i,6) = nil
                 iac(i,7) = nil
                 iac(i,8) = nil
              enddo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(i,iPar,id,xx,yy,zz,wp,Posx,Posy,Posz,xdum,ydum,zdum)
C$OMP+PRIVATE(ilc,iAcc,diff_x,diff_y,diff_z,t1,t2,t3)
              do i = 1 , npch
                 id = ipc(i)
                 xx = x(id)
                 yy = y(id)
                 zz = z(id)
                 wp = pw(id)
c     
c....         find neighbors to which this particle can contribute
c              
                 xdum = xx - hSize              
                 ydum = yy - hSize
                 zdum = zz - hSize
                 if ( xdum .le. 1.0 ) xdum = xdum + yn
                 if ( ydum .le. 1.0 ) ydum = ydum + yn
                 if ( zdum .le. 1.0 ) zdum = zdum + yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,1) = ilc
                 
                 xdum = xx + hSize
                 ydum = yy - hSize
                 zdum = zz - hSize
                 if ( xdum .ge.  xn ) xdum = xdum - yn
                 if ( ydum .le. 1.0 ) ydum = ydum + yn
                 if ( zdum .le. 1.0 ) zdum = zdum + yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,2) = ilc

                 xdum = xx - hSize
                 ydum = yy + hSize
                 zdum = zz - hSize
                 if ( xdum .le. 1.0 ) xdum = xdum + yn
                 if ( ydum .ge. xn  ) ydum = ydum - yn
                 if ( zdum .le. 1.0 ) zdum = zdum + yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,3) = ilc
                 
                 xdum = xx + hSize
                 ydum = yy + hSize
                 zdum = zz - hSize
                 if ( xdum .ge. xn  ) xdum = xdum - yn
                 if ( ydum .ge. xn  ) ydum = ydum - yn
                 if ( zdum .le. 1.0 ) zdum = zdum + yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,4) = ilc

                 xdum = xx - hSize
                 ydum = yy - hSize
                 zdum = zz + hSize
                 if ( xdum .le. 1.0 ) xdum = xdum + yn
                 if ( ydum .le. 1.0 ) ydum = ydum + yn
                 if ( zdum .ge. xn  ) zdum = zdum - yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,5) = ilc

                 xdum = xx + hSize
                 ydum = yy - hSize
                 zdum = zz + hSize
                 if ( xdum .ge. xn  ) xdum = xdum - yn
                 if ( ydum .le. 1.0 ) ydum = ydum + yn
                 if ( zdum .ge. xn  ) zdum = zdum - yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,6) = ilc
                 
                 xdum = xx - hSize
                 ydum = yy + hSize
                 zdum = zz + hSize
                 if ( xdum .le. 1.0 ) xdum = xdum + yn
                 if ( ydum .ge. xn  ) ydum = ydum - yn
                 if ( zdum .ge. xn  ) zdum = zdum - yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,7) = ilc

                 xdum = xx + hSize
                 ydum = yy + hSize
                 zdum = zz + hSize
                 if ( xdum .ge. xn  ) xdum = xdum - yn
                 if ( ydum .ge. xn  ) ydum = ydum - yn
                 if ( zdum .ge. xn  ) zdum = zdum - yn
                 ilc = iFindCell0 ( Level , xdum, ydum, zdum ) 
                 iac(i,8) = ilc
                
                 do iAcc = 1 , 8 
                    if ( iac(i,iAcc) .gt. nil ) then ! if iLv(iac) = Level 
                       call Ps ( iac(i,iAcc), Posx, Posy, Posz ) 
                       diff_x = xx - Posx
                       diff_y = yy - Posy
                       diff_z = zz - Posz
                       if ( abs(diff_x) .gt. nf67 ) then ! boundary effects
                          if ( diff_x .gt. zero ) then ! nf67 = ng / 2
                             diff_x = diff_x - ng
                          else
                             diff_x = diff_x + ng
                          endif
                       endif 
                       if ( abs(diff_y) .gt. nf67 ) then              
                          if ( diff_y .gt. zero ) then
                             diff_y = diff_y - ng
                          else
                             diff_y = diff_y + ng
                          endif
                       endif 
                       if ( abs(diff_z) .gt. nf67 ) then              
                          if ( diff_z .gt. zero ) then
                             diff_z = diff_z - ng
                          else
                             diff_z = diff_z + ng
                          endif
                       endif 
                       t1 = abs(diff_x) * Sizei
                       t2 = abs(diff_y) * Sizei
                       t3 = abs(diff_z) * Sizei
                       dpc(i,iAcc) = (1.0-t3)*(1.0-t1)*(1.0-t2)* 
     &                      wp
                    endif
                 enddo
              enddo

              do i = 1 , npch
                 do iAcc = 1 , 8 
                    ilc = iac(i,iAcc) 
                    if ( ilc .gt. nil ) then 
                       var(1,ilc) = var(1,ilc) + dpc(i,iAcc)
                       hvar(1,ilc) = hvar(1,ilc) + dpc(i,iAcc)
                       i_pch = ipc(i)
                  hvar(3,ilc) =hvar(3,ilc) + dpc(i,iAcc)*Vx(i_pch)
                  hvar(4,ilc) =hvar(4,ilc) + dpc(i,iAcc)*Vy(i_pch)
                  hvar(5,ilc) =hvar(5,ilc) + dpc(i,iAcc)*Vz(i_pch)
#ifdef STARFORM 
                       if ( ipc(i) .le. nsp(1,2) .or. 
     &                      ipc(i) .ge. nsp(nspecies,1) ) then 
#endif
#ifndef STARFORM 
                       if ( ipc(i) .le. nsp(1,2) ) then 
#endif
                          ref(ilc) = ref(ilc) + dpc(i,iAcc) ! DM Mass inside that cell.
                       endif
                    endif
                 enddo
              enddo
         npch = 0
      endif
 10   continue
      ENDDO        
      
      ENDDO
      return
      end
c
c     -------------------------------------
      subroutine Assign_ParticleDensity_AND_GAS0 ()
c     -------------------------------------
c
c     assigns density on the root mesh
c     algorithm: cloud-in-cell 
c
      include 'a_tree.h'
 

c.... Loop over particles

      do in = 1 ,  nsp(1,2) ! 1st specie first to avoid checking for ref

         xx = x(in)
         yy = y(in)
         zz = z(in)

         Vxp = Vx(in) 
         Vyp = Vy(in) 
         Vzp =Vz(in) 
c          Vzp =SNGL ( Vz(in) )

         i  = int(xx - 0.5)     ! 0.5 - to make my mesh consistent with CIC
         j  = int(yy - 0.5)
         k  = int(zz - 0.5)
         
         i1  = i + 1
         j1  = j + 1
         k1  = k + 1
         
         d1  = (xx - i) -0.5
         d2  = (yy - j) -0.5
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

c.... get indices of cell vertices

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

        hvar(1,ind1) =  hvar(1,ind1) +t3 * t1 * t2w
         hvar(1,ind2) = hvar(1,ind2) + t3 * d1 * t2w
         hvar(1,ind3) = hvar(1,ind3) + t3 * t1 * d2w
         hvar(1,ind4) =  hvar(1,ind4) + t3 * d1 * d2w
         hvar(1,ind5) = hvar(1,ind5) + d3 * t1 * t2w
         hvar(1,ind6) = hvar(1,ind6) +  d3 * d1 * t2w
         hvar(1,ind7) = hvar(1,ind7) +  d3 * t1 * d2w
         hvar(1,ind8) = hvar(1,ind8) + d3 * d1 * d2w     

         hvar(3,ind1) =  hvar(3,ind1) +t3 * t1 * t2w*Vxp
         hvar(3,ind2) = hvar(3,ind2) + t3 * d1 * t2w*Vxp
         hvar(3,ind3) = hvar(3,ind3) + t3 * t1 * d2w*Vxp
         hvar(3,ind4) =  hvar(3,ind4) + t3 * d1 * d2w*Vxp
         hvar(3,ind5) = hvar(3,ind5) + d3 * t1 * t2w*Vxp
         hvar(3,ind6) = hvar(3,ind6) +  d3 * d1 * t2w*Vxp
         hvar(3,ind7) = hvar(3,ind7) +  d3 * t1 * d2w*Vxp
         hvar(3,ind8) = hvar(3,ind8) + d3 * d1 * d2w*Vxp

         hvar(4,ind1) =  hvar(4,ind1) +t3 * t1 * t2w*Vyp
         hvar(4,ind2) = hvar(4,ind2) + t3 * d1 * t2w*Vyp
         hvar(4,ind3) = hvar(4,ind3) + t3 * t1 * d2w*Vyp
         hvar(4,ind4) =  hvar(4,ind4) + t3 * d1 * d2w*Vyp
         hvar(4,ind5) = hvar(4,ind5) + d3 * t1 * t2w*Vyp
         hvar(4,ind6) = hvar(4,ind6) +  d3 * d1 * t2w*Vyp
         hvar(4,ind7) = hvar(4,ind7) +  d3 * t1 * d2w*Vyp
         hvar(4,ind8) = hvar(4,ind8) + d3 * d1 * d2w*Vyp

         hvar(5,ind1) =  hvar(5,ind1) +t3 * t1 * t2w*Vzp
         hvar(5,ind2) = hvar(5,ind2) + t3 * d1 * t2w*Vzp
         hvar(5,ind3) = hvar(5,ind3) + t3 * t1 * d2w*Vzp
         hvar(5,ind4) =  hvar(5,ind4) + t3 * d1 * d2w*Vzp
         hvar(5,ind5) = hvar(5,ind5) + d3 * t1 * t2w*Vzp
         hvar(5,ind6) = hvar(5,ind6) +  d3 * d1 * t2w*Vzp
         hvar(5,ind7) = hvar(5,ind7) +  d3 * t1 * d2w*Vzp
         hvar(5,ind8) = hvar(5,ind8) + d3 * d1 * d2w*Vzp
         
      enddo

      do in = nsp(1,2)+1 , np

         xx = x(in)
         yy = y(in)
         zz = z(in)
         Vxp = Vx(in) 
         Vyp = Vy(in) 
         Vzp =Vz(in) 
         i  = int(xx-0.5)       ! 0.5 - to make my mesh consistent with CIC
         j  = int(yy-0.5)
         k  = int(zz-0.5)

        i1  = i + 1
        j1  = j + 1
        k1  = k + 1

        d1  = (xx - i)  -0.5
        d2  = (yy - j)  -0.5
        d3  = (zz - k) -0.5
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

        hvar(1,ind1) =  hvar(1,ind1) +t3 * t1 * t2w
         hvar(1,ind2) = hvar(1,ind2) + t3 * d1 * t2w
         hvar(1,ind3) = hvar(1,ind3) + t3 * t1 * d2w
         hvar(1,ind4) =  hvar(1,ind4) + t3 * d1 * d2w
         hvar(1,ind5) = hvar(1,ind5) + d3 * t1 * t2w
         hvar(1,ind6) = hvar(1,ind6) +  d3 * d1 * t2w
         hvar(1,ind7) = hvar(1,ind7) +  d3 * t1 * d2w
         hvar(1,ind8) = hvar(1,ind8) + d3 * d1 * d2w       

         hvar(3,ind1) =  hvar(3,ind1) +t3 * t1 * t2w*Vxp
         hvar(3,ind2) = hvar(3,ind2) + t3 * d1 * t2w*Vxp
         hvar(3,ind3) = hvar(3,ind3) + t3 * t1 * d2w*Vxp
         hvar(3,ind4) =  hvar(3,ind4) + t3 * d1 * d2w*Vxp
         hvar(3,ind5) = hvar(3,ind5) + d3 * t1 * t2w*Vxp
         hvar(3,ind6) = hvar(3,ind6) +  d3 * d1 * t2w*Vxp
         hvar(3,ind7) = hvar(3,ind7) +  d3 * t1 * d2w*Vxp
         hvar(3,ind8) = hvar(3,ind8) + d3 * d1 * d2w*Vxp

         hvar(4,ind1) =  hvar(4,ind1) +t3 * t1 * t2w*Vyp
         hvar(4,ind2) = hvar(4,ind2) + t3 * d1 * t2w*Vyp
         hvar(4,ind3) = hvar(4,ind3) + t3 * t1 * d2w*Vyp
         hvar(4,ind4) =  hvar(4,ind4) + t3 * d1 * d2w*Vyp
         hvar(4,ind5) = hvar(4,ind5) + d3 * t1 * t2w*Vyp
         hvar(4,ind6) = hvar(4,ind6) +  d3 * d1 * t2w*Vyp
         hvar(4,ind7) = hvar(4,ind7) +  d3 * t1 * d2w*Vyp
         hvar(4,ind8) = hvar(4,ind8) + d3 * d1 * d2w*Vyp

         hvar(5,ind1) =  hvar(5,ind1) +t3 * t1 * t2w*Vzp
         hvar(5,ind2) = hvar(5,ind2) + t3 * d1 * t2w*Vzp
         hvar(5,ind3) = hvar(5,ind3) + t3 * t1 * d2w*Vzp
         hvar(5,ind4) =  hvar(5,ind4) + t3 * d1 * d2w*Vzp
         hvar(5,ind5) = hvar(5,ind5) + d3 * t1 * t2w*Vzp
         hvar(5,ind6) = hvar(5,ind6) +  d3 * d1 * t2w*Vzp
         hvar(5,ind7) = hvar(5,ind7) +  d3 * t1 * d2w*Vzp
         hvar(5,ind8) = hvar(5,ind8) + d3 * d1 * d2w*Vzp
         
      enddo



      return
      end
c
c
c     -----------------------------------------------
      subroutine Assign_GAS_Final (Level ,vcorr)
c     -----------------------------------------------
c
c     NOFlag = 0 : assign old density in hvar(1,*) ( default)
c
c     purpose: 
c           Finish the assignment of gas fields
c           initializes density on a given Level with gas density
c
      integer Level, NOFlag  
      include 'a_tree.h'
      include 'a_control.h'
c
      b_f = Omb0 / Om0
      gamma1 = gamma -1.
      sum    = 0.
      sum2   = 0.
      rhomax = -1.e6
      rhomin = -rhomax
      vxmin =  -rhomax
      vxmax = -vxmin
      vymin =  vxmin
      vymax = -vymin
      vzmin =  vxmin
      vzmax = -vzmin
      e_min = rhomin
      e_max = -e_min
c   
c...  initial temperature in code units
c          
               a_th = 1.0 / (1.e3*(Omb0*hubble**2)**0.4)

               if ( aexpn .lt. a_th ) then 
                  TinitK = 2.726 / aexpn
               else
                  TinitK = 2.726 / a_th * (a_th/aexpn)**2
               endif     
               Tinit = aexpn**2 * TinitK / T_0 / gamma1 ! specific internal energy in code Units
c
      IF ( Level .eq. MinLevel ) THEN
c
!C$OMP PARALLEL DO DEFAULT(SHARED)
            do ic1 = 1 , ncell0
               rho00 =  hvar(1,ic1)
               vx0 = hvar(3,ic1)
               vy0= hvar(4,ic1)
               vz0 = hvar(5,ic1)

               sum2 = sum2 + rho00
               rho00 =max (rho00*b_f ,0.)            
               sum  = sum  + rho00

               vx0 = vx0*b_f*vcorr
               vy0 =  vy0*b_f*vcorr
               vz0 = vz0*b_f*vcorr

               var(1,ic1) = var(1,ic1) + rho00

               eps         = Tinit * rho00 ! internal energy density.
               Ekin        = (vx0**2+vy0**2+vz0**2)/ rho00 / 2. ! Kinetic energy density
              

               rhomin = min(rho00,rhomin)
               rhomax = max(rho00,rhomax)
               vxmin   = min(vxmin,vx0)
               vxmax   = max(vxmax,vx0)
               vymin   = min(vymin,vy0)
               vymax   = max(vymax,vy0)
               vzmin   = min(vzmin,vz0)
               vzmax   = max(vzmax,vz0)
               e_min   = min(e_min,eps)
               e_max   = max(e_max,eps)

               hvar(1,ic1)=rho00
               hvar(2,ic1)=Ekin+eps
               hvar(3,ic1)=vx0
               hvar(4,ic1)=vy0
               hvar(5,ic1)=vz0
               hvar(8,ic1)= eps
               hvar(7,ic1) = gamma

c               hvar(izII,ic1) =0.
            enddo
c
         ELSE
c
         Size3 = CellSize(Level)**3
         call Select_Cells ( Level , nLevel ) 
!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP+PRIVATE(icell)
            do ic1 = 1 , nLevel
               icell = iSelect(ic1)
               do iadd=0,7
                  rho00 =  hvar(1,icell+iadd)
                  vx0 = hvar(3,icell+iadd)
                  vy0= hvar(4,icell+iadd)
                  vz0 = hvar(5,icell+iadd)

                  sum2 = sum2 + rho00

                  rho00=max(rho00*b_f /Size3 , 0.)
                  sum  = sum  + rho00*Size3
                  
                  var(1,icell+iadd)=  var(1,icell+iadd)+rho00*Size3
                  var(2,icell+iadd)=0.
                  var(3,icell+iadd)=0.                  
                  
                  vx0=vx0*b_f*vcorr / Size3
                  vy0=vy0*b_f*vcorr / Size3
                  vz0=vz0*b_f*vcorr / Size3

                  eps         = Tinit * rho00 ! internal energy density
                  Ekin        = (vx0**2+vy0**2+ vz0**2)/ rho00 / 2. ! kinetic energy density


               rhomin = min(rho00,rhomin)
               rhomax = max(rho00,rhomax)
               vxmin   = min(vxmin,vx0)
               vxmax   = max(vxmax,vx0)
               vymin   = min(vymin,vy0)
               vymax   = max(vymax,vy0)
               vzmin   = min(vzmin,vz0)
               vzmax   = max(vzmax,vz0)
               e_min   = min(e_min,eps)
               e_max   = max(e_max,eps)

               hvar(1,icell+iadd) =rho00
               hvar(2,icell+iadd) = Ekin + eps
               hvar(3,icell+iadd) =vx0
               hvar(4,icell+iadd) =vy0
               hvar(5,icell+iadd) =vz0
               hvar(8,icell+iadd)= eps
               hvar(7,icell+iadd) = gamma     

c               hvar(izII,icell+iadd) = 0.               
            enddo
         enddo
      ENDIF
c     range of variables:

      write (*,*)  'Level= ', Level
      write (*,*)  '---------------------'
      write(*,*) ' a_th =', a_th, Omb0, hubble
      write(*,*) 'Tinit =',Tinit, T_0
      write(*,*) 'rho_min,max =',rhomin,rhomax
      write(*,*) 'ei_min,max  =',e_min,e_max
      write(*,*) 'vxmin,max   =',vxmin,vxmax
      write(*,*) 'vymin,max   =',vymin,vymax
      write(*,*) 'vzmin,max   =',vzmin,vzmax

      write (*,*)  'Total mass= ', sum2
      write (*,*)  'Mass in gas= ', sum
      write (*,*)  'Baryonic fraction= ', sum/sum2
      write (*,*)  'compared with universal= ',  b_f 
c     
      return
      end
c$$$c
c$$$c-------------------------------------------------------------------- 
c$$$      SUBROUTINE Derefine_Gas ()
c$$$c
c$$$c     If a given specie occupied a cell in a higher refined level,
c$$$c     for example, a zero-mass-resolution particle sitting in the Level=1,
c$$$c     the gas properties of that cell are taken from the parent cell, 
c$$$c     as well as the gas properties of the childs in the same parent cell.
c$$$c
c$$$c-------------------------------------------------------------------- 
c$$$      integer iChild(nchild)
c$$$      double precision xx, yy, zz
c$$$
c$$$      NderefGas=0.
c$$$      call Get_MaxLevelNow ()
c$$$      do ispec=2,nspec
c$$$         do ip=nsp(ispec,1), nsp(ispec,2) ! particles of ispec
c$$$            Lp = MaxLevel +1 -ispec
c$$$            xx = X(ip)
c$$$            yy = Y(ip)
c$$$            zz = Z(ip)
c$$$            idcell = iFindCell0 (Lp+1,xx,yy,zz)
c$$$            if (idcell.gt.0) then ! Derefine_Gas
c$$$               iParent = iPr(idcell)
c$$$               hvar1 = hvar(1,iParent)
c$$$               hvar2 = hvar(2,iParent)
c$$$               hvar3 = hvar(3,iParent)
c$$$               hvar4 = hvar(4,iParent)        
c$$$               hvar5 = hvar(5,iParent)
c$$$               hvar7 = hvar(7,iParent)
c$$$               hvar8 = hvar(8,iParent)
c$$$               call  iChAll ( iParent , iChild )
c$$$               do iCh=1,nchild
c$$$                  NderefGas = NderefGas + 1       
c$$$                  ichildren= iChild(iCh)
c$$$                  hvar(1,ichildren) = hvar1
c$$$                  hvar(2,ichildren) = hvar2
c$$$                  hvar(3,ichildren) = hvar3
c$$$                  hvar(4,ichildren) = hvar4       
c$$$                  hvar(5,ichildren) = hvar5
c$$$                  hvar(7,ichildren) = hvar7
c$$$                  hvar(8,ichildren) = hvar8    
c$$$               enddo
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$      write (*,*)  'NderefGas= ', NderefGas
c$$$      return
c$$$      end
