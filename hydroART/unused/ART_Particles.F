c     =====================================================================
c                                                                         .
c                 ART Version 3:  Particle routines                       .
c                                                                         .
c                      by Andrey Kravtsov (1997)                          .
c                                                                         .
c     =====================================================================

c     -------------------------------
      subroutine Move_Level ( Level )
c     -------------------------------
c
c     purpose: moves particles on a given Level
c     input  : Level  - level to move
c
#     include "a_def.h"
      include 'a_tree.h'
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
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /

      data iODir / 2 , 1 , 4 , 3 , 6 , 5 / 
     
      double precision xx , yy , zz , vvx , vvy , vvz
      double precision diff_x,diff_y,diff_z, delta_x
      double precision pconst1, a2cs(MinLevel:MaxLevel)

#ifdef FEEDBACK 
      real*8 erel, edum, efact, T_fact, T_ceil, cfbe
      real*8 t_fb_code, t_fbi_code, dum1, dum2
c------STG changed T_ceil to 2e7 K for speedup:
      parameter ( T_ceil = 2.0d7 ) ! do not allow gas to be feedback-heated above this T [K]
c--------------------------
#endif
#ifdef FEEDBACK_SNIa
      real*8 t_SNIa_code, R_SNIa_fact, f_SNIa, dN_SNIa, cfbeIa, xd
#endif
#ifdef ENRICH_SNIa
      real*8 ejM_SNIa_code
#endif
#ifdef STELLARMASSLOSS
      real*8 rhor, dmloss, T0_ml_code, rhofact, e_old      
#endif

      MinMove =  10000
      MaxMove = -10000

      do L = MinLevel , MaxLevel 
        a2cs(L) = -5.d-1*aexp(L)*aexp(L) / CellSize(L)
      enddo

      if ( Level .eq. MinLevel ) then 
        call Move_Min_Level ()
        return
      endif

      Size      = CellSize(Level)
      t_level   = tl    (Level)
      dt_level  = dtl   (Level)
      t_next    = t_level + dt_level
      t_limit   = t_next - 0.5 * dtmin
      t_rel     = t_next - t
      pconst    = 0.5 * aexp(Level)**2

#ifdef STARFORM
      nstar1 = nsp(nspecies,1) - 1 
#endif

#ifdef FEEDBACK 
      cfbe = C_fb * aexp(Level)**2 * CellVoli(Level)
      t_fb_code = t_fb / t0 / aexp(Level)**2
      t_fbi_code = 1.0 / t_fb_code
      T_fact = T_0 / aexp(Level)**2
      efact = T_ceil / T_fact / (gamma - 1.0) 
#ifdef FEEDBACK_SNIa
      t_SNIa_code = 1.d9 * t_SNIa / t0 / aexp(Level)**2
      R_SNIa_fact = RIaf * aexp(Level)**2
      cfbeIa = C_fbIa * CellVoli(Level) * aexp(Level)**2 
#ifdef ENRICH_SNIa
      ejM_SNIa_code = ejM_SNIa / aM0 * CellVoli(Level)
#endif
#endif
#endif
#ifdef STELLARMASSLOSS
      T0_ml_code = T0_ml * 1.d6 / t0 / aexp(Level)**2
#endif

      call Select_Cells ( Level , nLevel )

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1,ic2,icell,idcell,iParent,iHead,idummy)
C$OMP+PRIVATE (xx,yy,zz,vvx,vvy,vvz,gx,gy,gz,pdum,Level1,Size,iPar)
C$OMP+PRIVATE (ilc,icn,iDim,idcell1)
C$OMP+PRIVATE (diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE (xdum,ydum,zdum,d1,d2,d3,t1,t2,t3)
C$OMP+PRIVATE (t2t1,t2d1,d2t1,d2d1,inb2,inb4,inb24,iPack)
C$OMP+PRIVATE ( ic3,iWorkLevel,iDir,iNb_p,iNb_n,f_p,f_n,g_x,g_y,g_z)
C$OMP+PRIVATE (tp,dtp,pconst1,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1)
C$OMP+PRIVATE (d3t2t1,d3t2d1,d3d2t1,d3d2d1,delta_x,idcold)
C$OMP+PRIVATE (Posx,Posy,Posz,ifl,jfl,kfl)
C$OMP+PRIVATE (tfbc,dte,efrac,erel,edum,istar,fmetej)
C$OMP+PRIVATE (rhor,dmloss,rhofact,e_old,xd,dN_SNIa)
C$OMP+PRIVATE (CellVolume,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8)
      do ic1 = 1 , nLevel
        icell = iSelect(ic1)
        do ic2 = 0 , nchild-1
          idcell  = icell + ic2

          iParent = iPr(idcell)
          iHead   = iCL(idcell)

          if ( iHead .gt. nil ) then
            idummy = iHead
#ifdef FEEDBACK
            edum = T_fact * hvar(6,idcell) / hvar(1,idcell)
#endif
            do while ( idummy .ne. nil )

              tp  = pt(idummy)
              dtp = pdt(idummy)
              if ( tp .lt. t_limit ) then 
                xx  = x (idummy) 
                yy  = y (idummy)
                zz  = z (idummy) 
                vvx = vx(idummy)
                vvy = vy(idummy)
                vvz = vz(idummy)

#ifdef GRAVITY 
                Level1 = Level
                idcell1= idcell
 2553           iPar   = idcell1
                ilc    = idcell1
                Size   = CellSizei(Level1)
                call Ps ( iPar , Posx , Posy , Posz )
                if ( xx .eq. Posx ) xx = xx + dsmall
                if ( yy .eq. Posy ) yy = yy + dsmall
                if ( zz .eq. Posz ) zz = zz + dsmall
                ifl    = sign ( one , sngl(xx - Posx) )
                jfl    = sign ( one , sngl(yy - Posy) )
                kfl    = sign ( one , sngl(zz - Posz) )
                icn = (ifl + 1)/2 + jfl + 2 * kfl + 4  ! child number 

                do iDim = 1 , nDim
                  if ( iIP(icn,iDim) .eq. nil ) goto 2555
                  ilc = iNb(ilc,iIP(icn,iDim))              
                enddo

 2555           call Ps ( ilc, Posx, Posy, Posz )
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
                    corr_x = -ng
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
c
c....           compute weights
c
                d1   = abs(xdum - Posx) * Size
                d2   = abs(ydum - Posy) * Size
                d3   = abs(zdum - Posz) * Size
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
                      f_p(iDim) = wbig   * var(2,iNb_p) + 
     &                            wsmall * var(2,iPack(ic3))
                    else
                      f_p(iDim) = var(2,iNb_p)
                    endif
                    if (iLv(iNb_n) .lt. iWorkLevel ) then
                      f_n(iDim) = wbig   * var(2,iNb_n) + 
     &                            wsmall * var(2,iPack(ic3))
                    else
                      f_n(iDim) = var(2,iNb_n)
                    endif
                  enddo
                  g_x(ic3) = f_p(1) - f_n(1)
                  g_y(ic3) = f_p(2) - f_n(2)
                  g_z(ic3) = f_p(3) - f_n(3)
                 
                enddo
c
c....           interpolate accelerations to the particle location
c
                pconst1 = a2cs(iWorkLevel) * 
     &                    ( t_level - tp + (dt_level + dtp)*0.5 ) 

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

                vvx   = vvx + gx 
                vvy   = vvy + gy 
                vvz   = vvz + gz

c....           get potential energy at the particle location
                
      pdummy(idummy) = ( t3t2t1 * var(2,iPack(1)) +
     &                   t3t2d1 * var(2,iPack(2)) + 
     &                   t3d2t1 * var(2,iPack(3)) + 
     &                   t3d2d1 * var(2,iPack(4)) +
     &                   d3t2t1 * var(2,iPack(5)) + 
     &                   d3t2d1 * var(2,iPack(6)) + 
     &                   d3d2t1 * var(2,iPack(7)) + 
     &                   d3d2d1 * var(2,iPack(8)) ) / pconst1

cc...           get total density at the particle location

                t3t2t1  = t3 * t2t1
                t3t2d1  = t3 * t2d1
                t3d2t1  = t3 * d2t1
                t3d2d1  = t3 * d2d1
                d3t2t1  = d3 * t2t1
                d3t2d1  = d3 * t2d1
                d3d2t1  = d3 * d2t1
                d3d2d1  = d3 * d2d1

                CellVolume = CellSize(Level1)**3
                rho1 = var(1,iPack(1)) / CellVolume + 1.
                rho2 = var(1,iPack(2)) / CellVolume + 1.
                rho3 = var(1,iPack(3)) / CellVolume + 1.
                rho4 = var(1,iPack(4)) / CellVolume + 1.
                rho5 = var(1,iPack(5)) / CellVolume + 1.
                rho6 = var(1,iPack(6)) / CellVolume + 1.
                rho7 = var(1,iPack(7)) / CellVolume + 1.
                rho8 = var(1,iPack(8)) / CellVolume + 1.
                if ( Level1 .ne. iLv(iPack(1)) .or. 
     &               Level1 .ne. iLv(iPack(2)) .or. 
     &               Level1 .ne. iLv(iPack(3)) .or. 
     &               Level1 .ne. iLv(iPack(4)) .or. 
     &               Level1 .ne. iLv(iPack(5)) .or. 
     &               Level1 .ne. iLv(iPack(6)) .or. 
     &               Level1 .ne. iLv(iPack(7)) .or. 
     &               Level1 .ne. iLv(iPack(8)) ) then
                  write(*,*) ' Move_Level : Level1 .ne. iLv(iPack(i)) ',
     &                 Level1,iLv(iPack(1)),iLv(iPack(2)),iLv(iPack(3)),
     &                 iLv(iPack(4)),iLv(iPack(5)),iLv(iPack(6)),
     &                 iLv(iPack(7)),iLv(iPack(8))
                endif
c              write(*,*) 'Level1, CellVolume =',Level1,CellVolume
             
      ddummy(idummy) = ( t3t2t1 * rho1 + t3t2d1 * rho2 + 
     &                   t3d2t1 * rho3 + t3d2d1 * rho4 +
     &                   d3t2t1 * rho5 + d3t2d1 * rho6 + 
     &                   d3d2t1 * rho7 + d3d2d1 * rho8 )
c               write(*,*) 'components :'
c               write(*,*) t3t2t1, t3t2d1, t3d2t1, d3d2d1
c               write(*,*) d3t2t1, d3t2t1, d3t2d1, d3d2d1
c               write(*,*) rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8
c               write(*,*) 'idum, ddum =',idummy,ddummy(idummy)

#endif

                delta_x = t_next - tp

#ifdef STARFORM 
c
c....           release feedback energy from the recently (t-tbirth)<t_fb
c               born stellar particles into parent cell 
c
                if ( idummy .ge. nsp(nspecies,1) ) then 
                  istar = idummy - nstar1
                  tfbc = tp - tbirth(istar) 
#ifdef FEEDBACK
                  if ( tfbc .le. t_fb_code ) then 
                    dte = min ( delta_x , 
     &                         (tbirth(istar) + t_fb_code - tp) )
                    efrac = dte * t_fbi_code   ! delta_x / t_fb_code
#ifdef ENRICH
                    fmetej = efrac * fmass_met * pw0(istar)
                    hvar(izII,idcell) = hvar(izII,idcell) + 
     &                               fmetej * CellVoli(Level)
#endif
                    erel = efrac * cfbe * pw0(istar)
c
c....               don't allow to explode in an already hot bubble
c                   and limit energy release 
c
                    erel = min ( erel , 
     &                             efact * hvar(1,idcell) )
                    if ( edum .lt. T_ceil ) then  
                      hvar(2,idcell) = hvar(2,idcell) + erel 
                      hvar(8,idcell) = hvar(8,idcell) + erel 
                      hvar(6,idcell) = dble((hvar(7,idcell)-1.0) *
     &                                       hvar(8,idcell))
                    endif
                  endif                  
#endif
#ifdef FEEDBACK_SNIa
                  xd = t_SNIa_code / (t_next - tbirth(istar))
                  dN_SNIa = R_SNIa_fact * 
     &                   dexp(-xd*xd) * xd * dsqrt(xd) * 
     &                   delta_x * pw0(istar)
                  erel = dN_SNIa * cfbeIa
                  erel = min ( erel , 
     &                         efact * hvar(1,idcell) )
                  if ( edum .lt. T_ceil ) then 
                    hvar(2,idcell) = hvar(2,idcell) + erel 
                    hvar(8,idcell) = hvar(8,idcell) + erel 
                    hvar(6,idcell) = dble((hvar(7,idcell)-1.0) *
     &                                       hvar(8,idcell))
                  endif
#ifdef ENRICH_SNIa
                  hvar(izIa,idcell) = hvar(izIa,idcell) + 
     &                                dN_SNIa * ejM_SNIa_code
#endif
#endif
#ifdef STELLARMASSLOSS
                  dmloss = min ( pw0(istar) * delta_x * c0_ml / 
     &                           (tfbc + T0_ml_code) ,
     &                           1.d-1 * pw(idummy) )
                  pw(idummy) = pw(idummy) - dmloss 
                  dmloss = dmloss * CellVoli(Level) ! m -> rho
c
c....             account for the momentum change 

                  rhor = 1.d0 / hvar(1,idcell) 
                  e_old = dble(hvar(2,idcell)) - 
     &                    0.5d0 * ( hvar(3,idcell)**2 +
     &                              hvar(4,idcell)**2 +
     &                              hvar(5,idcell)**2   ) * rhor
                  hvar(1,idcell) = hvar(1,idcell) + dmloss 
                  rhofact = rhor * hvar(1,idcell)
                  hvar(3,idcell) = hvar(3,idcell) + dmloss * vvx 
                  hvar(4,idcell) = hvar(4,idcell) + dmloss * vvy
                  hvar(5,idcell) = hvar(5,idcell) + dmloss * vvz
                  hvar(2,idcell) = e_old  + 0.5d0 *  
     &                    (hvar(3,idcell)**2 +
     &                     hvar(4,idcell)**2 +
     &                     hvar(5,idcell)**2   ) / hvar(1,idcell)
                  hvar(6,idcell) = hvar(6,idcell) * rhofact 
                  hvar(8,idcell) = hvar(8,idcell) * rhofact
#ifdef ENRICH
c
c....             mass lost by the stellar population has metallicity
c                 with with the stars were born
c
                  hvar(izII,idcell) = hvar(izII,idcell) + 
     &                                zstII(istar) * dmloss 
#ifdef ENRICH_SNIa
                  hvar(izIa,idcell) = hvar(izIa,idcell) + 
     &                                zstIa(istar) * dmloss 
#endif
#endif                 
 
#endif
                endif
#endif

	        xx    = xx  + vvx * delta_x
	        yy    = yy  + vvy * delta_x
	        zz    = zz  + vvz * delta_x

                if ( xx .le. 1. ) xx = xx + yn ! imposing periodic 
                if ( xx .ge. xn ) xx = xx - yn ! boundary conditions
                if ( yy .le. 1. ) yy = yy + yn
                if ( yy .ge. xn ) yy = yy - yn
                if ( zz .le. 1. ) zz = zz + yn
                if ( zz .ge. xn ) zz = zz - yn

                 x (idummy) =  xx
                 y (idummy) =  yy
                 z (idummy) =  zz
                vx (idummy) = vvx
                vy (idummy) = vvy
                vz (idummy) = vvz
c                acx(idummy) = gx
c                acy(idummy) = gy
c                acz(idummy) = gz

                pt (idummy) = t_next
                pdt(idummy) = dt_level

              endif

              idummy     = iLL(idummy,1)

 1553       continue
           enddo  ! end do while
         endif
        enddo
      enddo

      return
      end
      
c     ----------------------------
      subroutine Move_Min_Level ()
c     ----------------------------
c
c     purpose: moves particles on the root level
c     
c     
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      dimension iPack(8)                     ! pack of cells to be used in 
                                             ! acceleration interpolation
      dimension g_x(8), g_y(8), g_z(8)
      dimension iIP(8,3)                     ! set of pointers to the leftmost
                                             ! cells - to be used as kernels
      data iIP / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &           3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &           5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /

      double precision xx , yy , zz , vvx , vvy , vvz, a2

      real*8 delta_x

#ifdef FEEDBACK 
      real*8 erel, edum, T_fact, T_ceil, efact
c------STG changed T_ceil to 2e7 K for speedup:
      parameter ( T_ceil = 2.0d7 ) ! do not allow gas to be feedback-heated above this T [K]
                                  ! T_ceil here should be the same as in Move_Level
c----------------------------
      real*8 t_fb_code, t_fbi_code

#endif
#ifdef FEEDBACK_SNIa
      real*8 t_SNIa_code, R_SNIa_fact, f_SNIa, dN_SNIa, cfbeIa, xd
#endif
#ifdef ENRICH_SNIa
      real*8 ejM_SNIa_code
#endif
#ifdef STELLARMASSLOSS
      real*8 rhor, dmloss, T0_ml_code, rhofact, e_old
#endif

      Size    = 1.0
      pconst1 = 1.0 / Size

      a2 = -5.d-1 * aexp(MinLevel)**2
      t_level   = tl    (MinLevel)
      dt_level  = dtl   (MinLevel)
      t_next    = t_level + dt_level
      t_limit   = t_next - 0.5 * dtmin
      t_rel     = t_next - t

#ifdef STARFORM
      nstar1 = nsp(nspecies,1) - 1 
#endif

#ifdef FEEDBACK 
      cfbe = C_fb * aexp(MinLevel)**2 
      t_fb_code = t_fb / t0 / aexp(MinLevel)**2
      t_fbi_code = 1.0 / t_fb_code
      T_fact = T_0 / aexp(MinLevel)**2
      efact = T_ceil / T_fact / (gamma - 1.0) 
#ifdef FEEDBACK_SNIa
      t_SNIa_code = 1.d9 * t_SNIa / t0 / aexp(MinLevel)**2
      R_SNIa_fact = RIaf * aexp(MinLevel)**2
      cfbeIa = C_fbIa * aexp(MinLevel)**2 
#ifdef ENRICH_SNIa
      ejM_SNIa_code = ejM_SNIa / aM0
#endif
#endif
#endif
#ifdef STELLARMASSLOSS
      T0_ml_code = T0_ml * 1.d6 / t0 / aexp(MinLevel)**2
#endif


C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE(ic1,idummy)
C$OMP+PRIVATE(xx,yy,zz,vvx,vvy,vvz,gx,gy,gz)
C$OMP+PRIVATE(iPar,ilc,icn,iDim)
C$OMP+PRIVATE(diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE(xdum,ydum,zdum,i,j,k,pos1,pos2,pos3)
C$OMP+PRIVATE(d1,d2,d3,t1,t2,t3,t2t1,t2d1,d2t1,d2d1)
C$OMP+PRIVATE(inb2,inb4,inb24,iPack,ic3,g_x,g_y,g_z)
C$OMP+PRIVATE(tp,dtp,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1)
C$OMP+PRIVATE(d3t2t1,d3t2d1,d3d2t1,d3d2d1,pconst1,delta_x)
C$OMP+PRIVATE(Posx,Posy,Posz,ifl,jfl,kfl)
C$OMP+PRIVATE(tfbc,dte,efrac,erel,edum,istar,fmetej)
C$OMP+PRIVATE(rhor,dmloss,rhofact,e_old,xd,dN_SNIa)
C$OMP+PRIVATE(CellVolume,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8)
      do ic1 = 1 , ncell0
        if ( iCL(ic1) .gt. nil ) then
          idummy = iCL(ic1)
          do while ( idummy .ne. nil )
            tp  = pt(idummy)
            dtp = pdt(idummy)
            if ( tp .lt. t_limit ) then 
              xx  = x(idummy) 
              yy  = y(idummy)
              zz  = z(idummy) 
              vvx = vx(idummy)
              vvy = vy(idummy)
              vvz = vz(idummy)

#ifdef GRAVITY             
              iPar = ic1
              ilc  = iPar
              call Ps ( iPar, Posx, Posy, Posz )
              if ( xx .eq. Posx ) xx = xx + dsmall
              if ( yy .eq. Posy ) yy = yy + dsmall
              if ( zz .eq. Posz ) zz = zz + dsmall
              ifl = sign ( one , (sngl(xx) - Posx) )
              jfl = sign ( one , (sngl(yy) - Posy) )
              kfl = sign ( one , (sngl(zz) - Posz) )
              icn = (ifl + 1)/2 + jfl + 2 * kfl + 4

              do iDim = 1 , nDim
                if ( iIP(icn,iDim) .eq. nil ) goto 2556
                ilc = iNb(ilc,iIP(icn,iDim))             
              enddo

 2556         call Ps ( ilc, Posx, Posy, Posz )
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
c
c....         compute position of lefmost corner cell (ilc)
c
              i = (ilc - 1)/ng2 + 1
              j = (ilc - (i-1)*ng2 - 1)/ng + 1
              k =  ilc - (i-1)*ng2 - (j-1)*ng

              pos1 = 1.0 * i + 0.5
              pos2 = 1.0 * j + 0.5
              pos3 = 1.0 * k + 0.5
c
c....         compute CIC coefficients
c
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
                g_x(ic3) = var(2,iNb(iPack(ic3),2)) - 
     &                     var(2,iNb(iPack(ic3),1))
                g_y(ic3) = var(2,iNb(iPack(ic3),4)) - 
     &                     var(2,iNb(iPack(ic3),3))
                g_z(ic3) = var(2,iNb(iPack(ic3),6)) - 
     &                     var(2,iNb(iPack(ic3),5))

              enddo

c....         interpolate accelerations to the particle location

              pconst1 = a2 * ((t_level - tp) + (dt_level + dtp)*0.5)
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

              vvx   = vvx + gx
              vvy   = vvy + gy
              vvz   = vvz + gz

c....         get potential energy at the particle location

      pdummy(idummy) = ( t3t2t1 * var(2,iPack(1)) +
     &                   t3t2d1 * var(2,iPack(2)) + 
     &                   t3d2t1 * var(2,iPack(3)) + 
     &                   t3d2d1 * var(2,iPack(4)) +
     &                   d3t2t1 * var(2,iPack(5)) + 
     &                   d3t2d1 * var(2,iPack(6)) + 
     &                   d3d2t1 * var(2,iPack(7)) + 
     &                   d3d2d1 * var(2,iPack(8)) ) / pconst1

c...          get total density at the particle location

              t3t2t1  = t3 * t2t1
              t3t2d1  = t3 * t2d1
              t3d2t1  = t3 * d2t1
              t3d2d1  = t3 * d2d1
              d3t2t1  = d3 * t2t1
              d3t2d1  = d3 * t2d1
              d3d2t1  = d3 * d2t1
              d3d2d1  = d3 * d2d1

CEVERINO12102005              CellVolume = CellSize(Level1)**3
              CellVolume = CellSize(MinLevel)**3
              rho1 = var(1,iPack(1)) / CellVolume + 1.
              rho2 = var(1,iPack(2)) / CellVolume + 1.
              rho3 = var(1,iPack(3)) / CellVolume + 1.
              rho4 = var(1,iPack(4)) / CellVolume + 1.
              rho5 = var(1,iPack(5)) / CellVolume + 1.
              rho6 = var(1,iPack(6)) / CellVolume + 1.
              rho7 = var(1,iPack(7)) / CellVolume + 1.
              rho8 = var(1,iPack(8)) / CellVolume + 1.
      
      ddummy(idummy) = ( t3t2t1 * rho1 + t3t2d1 * rho2 + 
     &                   t3d2t1 * rho3 + t3d2d1 * rho4 +
     &                   d3t2t1 * rho5 + d3t2d1 * rho6 + 
     &                   d3d2t1 * rho7 + d3d2d1 * rho8 )

#endif

              delta_x = (t_next - tp) 

#ifdef STARFORM
c
c....           release feedback energy from the recently (t-tbirth)<t_fb
c               born stellar particles into parent cell 
c
                if ( idummy .ge. nsp(nspecies,1) ) then 
                  istar = idummy - nstar1
                  tfbc = tp - tbirth(istar) 
#ifdef FEEDBACK
                  edum = T_fact * hvar(6,ic1) / hvar(1,ic1)
                  if ( tfbc .le. t_fb_code ) then 
                    dte = min ( delta_x , 
     &                         (tbirth(istar) + t_fb_code - tp) )
                    efrac = dte * t_fbi_code   ! delta_x / t_fb_code
                    erel = efrac * cfbe * pw0(istar) 

#ifdef ENRICH
                    fmetej = efrac * fmass_met * pw0(istar)
                    hvar(izII,ic1) = hvar(izII,ic1) + fmetej  ! update SNII rho_metals
#endif
c
c....               limit increase in energy and don't allow to release E in an 
c                   already hot bubble
c
                    erel = min ( erel , 
     &                           efact * hvar(1,ic1) )

                    if ( edum .lt. T_ceil ) then 
                      hvar(8,ic1) = hvar(8,ic1) + erel
                      hvar(2,ic1) = hvar(2,ic1) + erel
                      hvar(6,ic1) = dble((hvar(7,ic1)-1.0)*hvar(8,ic1))
                    endif
                  endif                  
#endif
#ifdef FEEDBACK_SNIa
                  xd = t_SNIa_code / (t_next - tbirth(istar))
                  dN_SNIa = R_SNIa_fact * 
     &                   dexp(-xd*xd) * xd * dsqrt(xd) * 
     &                   delta_x * pw0(istar)
                  erel = dN_SNIa * cfbeIa

                  erel = min ( erel , 
     &                         efact * hvar(1,ic1) )
                  if ( edum .lt. T_ceil ) then 
                    hvar(2,ic1) = hvar(2,ic1) + erel 
                    hvar(8,ic1) = hvar(8,ic1) + erel 
                    hvar(6,ic1) = dble((hvar(7,ic1)-1.0) *
     &                                       hvar(8,ic1))
                  endif
#ifdef ENRICH_SNIa
                  hvar(izIa,ic1) = hvar(izIa,ic1) + 
     &                                dN_SNIa * ejM_SNIa_code
#endif
#endif
#ifdef STELLARMASSLOSS
                    dmloss = min ( pw0(istar) * delta_x * c0_ml / 
     &                             (tfbc + T0_ml_code) ,
     &                             1.d-1 * pw(idummy) )
                    pw(idummy) = pw(idummy) - dmloss 
c
c....               account for the momentum change 
c
                    rhor = 1.0 / hvar(1,ic1) 
                    e_old = dble(hvar(2,ic1)) - 
     &                      0.5d0 * ( hvar(3,ic1)**2 +
     &                                hvar(4,ic1)**2 +
     &                                hvar(5,ic1)**2   ) * rhor
                    hvar(1,ic1) = hvar(1,ic1) + dmloss 
                    rhofact = rhor * hvar(1,ic1)
                    hvar(3,ic1) = hvar(3,ic1) + dmloss * vvx 
                    hvar(4,ic1) = hvar(4,ic1) + dmloss * vvy
                    hvar(5,ic1) = hvar(5,ic1) + dmloss * vvz
                    hvar(2,ic1) = e_old * rhofact + 
     &                      0.5d0 * ( hvar(3,ic1)**2 +
     &                                hvar(4,ic1)**2 +
     &                                hvar(5,ic1)**2   ) / hvar(1,ic1)
                    hvar(6,ic1) = hvar(6,ic1) * rhofact 
                    hvar(8,ic1) = hvar(8,ic1) * rhofact                     
#ifdef ENRICH
c
c....             mass lost by the stellar population has metallicity
c                 with with the stars were born
c
                  hvar(izII,ic1) = hvar(izII,ic1) + 
     &                             zstII(istar) * dmloss 
#ifdef ENRICH_SNIa
                  hvar(izIa,ic1) = hvar(izIa,ic1) + 
     &                             zstIa(istar) * dmloss 
#endif
#endif                 
#endif
                endif
#endif

   	       xx   = xx  + vvx * delta_x
	       yy   = yy  + vvy * delta_x
	       zz   = zz  + vvz * delta_x

              if ( xx .le. 1. ) xx = xx + yn ! imposing periodic 
              if ( xx .ge. xn ) xx = xx - yn ! boundary conditions
              if ( yy .le. 1. ) yy = yy + yn
              if ( yy .ge. xn ) yy = yy - yn
              if ( zz .le. 1. ) zz = zz + yn
              if ( zz .ge. xn ) zz = zz - yn

              x  (idummy) = xx
              y  (idummy) = yy
              z  (idummy) = zz
              vx (idummy) = vvx
              vy (idummy) = vvy
              vz (idummy) = vvz
c              acx(idummy) = gx
c              acy(idummy) = gy
c              acz(idummy) = gz
              pt (idummy) = t_next
              pdt(idummy) = dt_level

            endif
            idummy         = iLL(idummy,1)

          enddo  ! end do while
        endif
      enddo

      return
      end

c     ---------------------------------------
      subroutine Assign_Dummy_Level ( Level )
c     ---------------------------------------
c
c     purpose: assign the local potential and overdensity
c              to the dummy variables on a given Level
c     input  : Level  - level to move
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      integer   Level 
      dimension iPack(8)                     ! pack of cells to be used in 
      dimension iODir(6)                     ! opposite directions
      dimension iIP(8,3)                     ! set of pointers to the leftmost
                                             ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /

      data iODir / 2 , 1 , 4 , 3 , 6 , 5 / 
     
      double precision xx , yy , zz , vvx , vvy , vvz
      double precision diff_x,diff_y,diff_z, delta_x
      double precision pconst1, a2cs(MinLevel:MaxLevel)

      do L = MinLevel , MaxLevel 
        a2cs(L) = -5.d-1*aexp(L)*aexp(L) / CellSize(L)
      enddo

      if ( Level .eq. MinLevel ) then 
        call Assign_Dummy_Min_Level ()
        return
      endif

      Size      = CellSize(Level)
      t_level   = tl    (Level)
      dt_level  = dtl   (Level)
      t_next    = t_level + dt_level
      t_limit   = t_next - 0.5 * dtmin
      t_rel     = t_next - t
      pconst    = 0.5 * aexp(Level)**2

      call Select_Cells ( Level , nLevel )

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (ic1,ic2,icell,idcell,iParent,iHead,idummy)
C$OMP+PRIVATE (xx,yy,zz,vvx,vvy,vvz,gx,gy,gz,pdum,Level1,Size,iPar)
C$OMP+PRIVATE (ilc,icn,iDim,idcell1)
C$OMP+PRIVATE (diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE (xdum,ydum,zdum,d1,d2,d3,t1,t2,t3)
C$OMP+PRIVATE (t2t1,t2d1,d2t1,d2d1,inb2,inb4,inb24,iPack)
C$OMP+PRIVATE ( ic3,iWorkLevel,iDir,iNb_p,iNb_n,f_p,f_n,g_x,g_y,g_z)
C$OMP+PRIVATE (tp,dtp,pconst1,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1)
C$OMP+PRIVATE (d3t2t1,d3t2d1,d3d2t1,d3d2d1,delta_x,idcold)
C$OMP+PRIVATE (Posx,Posy,Posz,ifl,jfl,kfl)
C$OMP+PRIVATE (tfbc,dte,efrac,erel,edum,istar,fmetej)
C$OMP+PRIVATE (rhor,dmloss,rhofact,e_old,xd,dN_SNIa)
C$OMP+PRIVATE (CellVolume,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8)
      do ic1 = 1 , nLevel
        icell = iSelect(ic1)
        do ic2 = 0 , nchild-1
          idcell  = icell + ic2

          iParent = iPr(idcell)
          iHead   = iCL(idcell)

          if ( iHead .gt. nil ) then
            idummy = iHead
            do while ( idummy .ne. nil )

              tp  = pt(idummy)
              dtp = pdt(idummy)
              if ( tp .lt. t_limit ) then 
                xx  = x (idummy) 
                yy  = y (idummy)
                zz  = z (idummy) 
                vvx = vx(idummy)
                vvy = vy(idummy)
                vvz = vz(idummy)

#ifdef GRAVITY 
                Level1 = Level
                idcell1= idcell
 2653           iPar   = idcell1
                ilc    = idcell1
                Size   = CellSizei(Level1)
                call Ps ( iPar , Posx , Posy , Posz )
                if ( xx .eq. Posx ) xx = xx + dsmall
                if ( yy .eq. Posy ) yy = yy + dsmall
                if ( zz .eq. Posz ) zz = zz + dsmall
                ifl    = sign ( one , sngl(xx - Posx) )
                jfl    = sign ( one , sngl(yy - Posy) )
                kfl    = sign ( one , sngl(zz - Posz) )
                icn = (ifl + 1)/2 + jfl + 2 * kfl + 4  ! child number 

                do iDim = 1 , nDim
                  if ( iIP(icn,iDim) .eq. nil ) goto 2655
                  ilc = iNb(ilc,iIP(icn,iDim))              
                enddo

 2655           call Ps ( ilc, Posx, Posy, Posz )
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
                    corr_x = -ng
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
c
c....           compute weights
c
                d1   = abs(xdum - Posx) * Size
                d2   = abs(ydum - Posy) * Size
                d3   = abs(zdum - Posz) * Size
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
                    goto 2653
                  endif
                enddo

                iWorkLevel = iLv(ilc)

c 
c....           interpolate accelerations to the particle location
c
                pconst1 = a2cs(iWorkLevel) * 
     &                    ( t_level - tp + (dt_level + dtp)*0.5 ) 

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

c....           get potential energy at the particle location
                
      pdummy(idummy) = ( t3t2t1 * var(2,iPack(1)) +
     &                   t3t2d1 * var(2,iPack(2)) + 
     &                   t3d2t1 * var(2,iPack(3)) + 
     &                   t3d2d1 * var(2,iPack(4)) +
     &                   d3t2t1 * var(2,iPack(5)) + 
     &                   d3t2d1 * var(2,iPack(6)) + 
     &                   d3d2t1 * var(2,iPack(7)) + 
     &                   d3d2d1 * var(2,iPack(8)) ) / pconst1

c...            get total density at the particle location

                t3t2t1  = t3 * t2t1
                t3t2d1  = t3 * t2d1
                t3d2t1  = t3 * d2t1
                t3d2d1  = t3 * d2d1
                d3t2t1  = d3 * t2t1
                d3t2d1  = d3 * t2d1
                d3d2t1  = d3 * d2t1
                d3d2d1  = d3 * d2d1
  
                CellVolume = CellSize(Level1)**3
                rho1 = var(1,iPack(1)) / CellVolume + 1.
                rho2 = var(1,iPack(2)) / CellVolume + 1.
                rho3 = var(1,iPack(3)) / CellVolume + 1.
                rho4 = var(1,iPack(4)) / CellVolume + 1.
                rho5 = var(1,iPack(5)) / CellVolume + 1.
                rho6 = var(1,iPack(6)) / CellVolume + 1.
                rho7 = var(1,iPack(7)) / CellVolume + 1.
                rho8 = var(1,iPack(8)) / CellVolume + 1.
                if ( Level1 .ne. iLv(iPack(1)) .or. 
     &               Level1 .ne. iLv(iPack(2)) .or. 
     &               Level1 .ne. iLv(iPack(3)) .or. 
     &               Level1 .ne. iLv(iPack(4)) .or. 
     &               Level1 .ne. iLv(iPack(5)) .or. 
     &               Level1 .ne. iLv(iPack(6)) .or. 
     &               Level1 .ne. iLv(iPack(7)) .or. 
     &               Level1 .ne. iLv(iPack(8)) ) then
                  write(*,*) ' Move_Level : Level1 .ne. iLv(iPack(i)) ',
     &                 Level1,iLv(iPack(1)),iLv(iPack(2)),iLv(iPack(3)),
     &                 iLv(iPack(4)),iLv(iPack(5)),iLv(iPack(6)),
     &                 iLv(iPack(7)),iLv(iPack(8))
                endif
c              write(*,*) 'Level1, CellVolume =',Level1,CellVolume
             
      ddummy(idummy) = ( t3t2t1 * rho1 + t3t2d1 * rho2 + 
     &                   t3d2t1 * rho3 + t3d2d1 * rho4 +
     &                   d3t2t1 * rho5 + d3t2d1 * rho6 + 
     &                   d3d2t1 * rho7 + d3d2d1 * rho8 )
c               if ( idummy .eq. 44953 ) then
c                 write(*,*) 'components :'
c                 write(*,*) t3t2t1, t3t2d1, t3d2t1, d3d2d1
c                 write(*,*) d3t2t1, d3t2t1, d3t2d1, d3d2d1
c               write(*,*) rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8
c                 write(*,*) 'idum, ddum =',idummy,ddummy(idummy)
c               endif

#endif

              endif

              idummy     = iLL(idummy,1)

 1653       continue
           enddo  ! end do while
         endif
        enddo
      enddo

      return
      end
      
c     ------------------------------------
      subroutine Assign_Dummy_Min_Level ()
c     ------------------------------------
c
c     purpose: assign the local potential and overdensity
c              to the dummy variables on a minimum Level
c     
c     
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      dimension iPack(8)                     ! pack of cells to be used in 
      dimension iIP(8,3)                     ! set of pointers to the leftmost
                                             ! cells - to be used as kernels
      data iIP / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &           3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &           5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /

      double precision xx , yy , zz , vvx , vvy , vvz, a2


      Size    = 1.0
      pconst1 = 1.0 / Size

      a2 = -5.d-1 * aexp(MinLevel)**2
      t_level   = tl    (MinLevel)
      dt_level  = dtl   (MinLevel)
      t_next    = t_level + dt_level
      t_limit   = t_next - 0.5 * dtmin
      t_rel     = t_next - t

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE(ic1,idummy)
C$OMP+PRIVATE(xx,yy,zz,vvx,vvy,vvz,gx,gy,gz)
C$OMP+PRIVATE(iPar,ilc,icn,iDim)
C$OMP+PRIVATE(diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE(xdum,ydum,zdum,i,j,k,pos1,pos2,pos3)
C$OMP+PRIVATE(d1,d2,d3,t1,t2,t3,t2t1,t2d1,d2t1,d2d1)
C$OMP+PRIVATE(inb2,inb4,inb24,iPack,ic3,g_x,g_y,g_z)
C$OMP+PRIVATE(tp,dtp,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1)
C$OMP+PRIVATE(d3t2t1,d3t2d1,d3d2t1,d3d2d1,pconst1,delta_x)
C$OMP+PRIVATE(Posx,Posy,Posz,ifl,jfl,kfl)
C$OMP+PRIVATE(tfbc,dte,efrac,erel,edum,istar,fmetej)
C$OMP+PRIVATE(rhor,dmloss,rhofact,e_old,xd,dN_SNIa)
C$OMP+PRIVATE(CellVolume,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8)
      do ic1 = 1 , ncell0
        if ( iCL(ic1) .gt. nil ) then
          idummy = iCL(ic1)
          do while ( idummy .ne. nil )
            tp  = pt(idummy)
            dtp = pdt(idummy)
            if ( tp .lt. t_limit ) then 
              xx  = x(idummy) 
              yy  = y(idummy)
              zz  = z(idummy) 
              vvx = vx(idummy)
              vvy = vy(idummy)
              vvz = vz(idummy)

#ifdef GRAVITY             
              iPar = ic1
              ilc  = iPar
              call Ps ( iPar, Posx, Posy, Posz )
              if ( xx .eq. Posx ) xx = xx + dsmall
              if ( yy .eq. Posy ) yy = yy + dsmall
              if ( zz .eq. Posz ) zz = zz + dsmall
              ifl = sign ( one , (sngl(xx) - Posx) )
              jfl = sign ( one , (sngl(yy) - Posy) )
              kfl = sign ( one , (sngl(zz) - Posz) )
              icn = (ifl + 1)/2 + jfl + 2 * kfl + 4

              do iDim = 1 , nDim
                if ( iIP(icn,iDim) .eq. nil ) goto 2656
                ilc = iNb(ilc,iIP(icn,iDim))             
              enddo

 2656         call Ps ( ilc, Posx, Posy, Posz )
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
c
c....         compute position of lefmost corner cell (ilc)
c
              i = (ilc - 1)/ng2 + 1
              j = (ilc - (i-1)*ng2 - 1)/ng + 1
              k =  ilc - (i-1)*ng2 - (j-1)*ng

              pos1 = 1.0 * i + 0.5
              pos2 = 1.0 * j + 0.5
              pos3 = 1.0 * k + 0.5
c
c....         compute CIC coefficients
c
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

c....         interpolate accelerations to the particle location

              pconst1 = a2 * ((t_level - tp) + (dt_level + dtp)*0.5)
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

c....         get potential energy at the particle location

      pdummy(idummy) = ( t3t2t1 * var(2,iPack(1)) +
     &                   t3t2d1 * var(2,iPack(2)) + 
     &                   t3d2t1 * var(2,iPack(3)) + 
     &                   t3d2d1 * var(2,iPack(4)) +
     &                   d3t2t1 * var(2,iPack(5)) + 
     &                   d3t2d1 * var(2,iPack(6)) + 
     &                   d3d2t1 * var(2,iPack(7)) + 
     &                   d3d2d1 * var(2,iPack(8)) ) / pconst1

c...          get total density at the particle location

              t3t2t1  = t3 * t2t1
              t3t2d1  = t3 * t2d1
              t3d2t1  = t3 * d2t1
              t3d2d1  = t3 * d2d1
              d3t2t1  = d3 * t2t1
              d3t2d1  = d3 * t2d1
              d3d2t1  = d3 * d2t1
              d3d2d1  = d3 * d2d1

              CellVolume = CellSize(Level1)**3
              rho1 = var(1,iPack(1)) / CellVolume + 1.
              rho2 = var(1,iPack(2)) / CellVolume + 1.
              rho3 = var(1,iPack(3)) / CellVolume + 1.
              rho4 = var(1,iPack(4)) / CellVolume + 1.
              rho5 = var(1,iPack(5)) / CellVolume + 1.
              rho6 = var(1,iPack(6)) / CellVolume + 1.
              rho7 = var(1,iPack(7)) / CellVolume + 1.
              rho8 = var(1,iPack(8)) / CellVolume + 1.
      
      ddummy(idummy) = ( t3t2t1 * rho1 + t3t2d1 * rho2 + 
     &                   t3d2t1 * rho3 + t3d2d1 * rho4 +
     &                   d3t2t1 * rho5 + d3t2d1 * rho6 + 
     &                   d3d2t1 * rho7 + d3d2d1 * rho8 )

#endif

            endif
            idummy         = iLL(idummy,1)

          enddo  ! end do while
        endif
      enddo

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
      include 'a_dnew.fh'
      
      integer Level , MinModify , MaxModify
      double precision xx, yy, zz

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
              if ( Next_Cell .ne. icell ) then
                Next_Level = iLv(Next_Cell)
                if ( Next_Level.lt.MinModify ) MinModify = Next_Level
                if ( Next_Level.gt.MaxModify ) MaxModify = Next_Level
                call LL_Delete ( icell     , idummy )
                call LL_Insert ( Next_Cell , idummy )
                iPL(idummy) = Next_Level 
              endif
              idummy = idummy2
            enddo ! end do while
          endif
        enddo ! end icell

      ELSE

        nLevel = iNOLL(Level)
        icell  = iGet_Cell_Index ( iHOLL(Level) ) 
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
                Next_Cell = iFindCell ( MaxLevel , xx , yy , zz )
                idummy2   = iLL(idummy,1)
                if ( Next_Cell .ne. idcell ) then
                    Next_Level = iLv(Next_Cell)
                   if ( Next_Level.lt.MinModify ) MinModify = Next_Level
                   if ( Next_Level.gt.MaxModify ) MaxModify = Next_Level
                    call LL_Delete ( idcell    , idummy )
                    call LL_Insert ( Next_Cell , idummy )
                    iPL(idummy) = Next_Level 
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
      include 'a_dnew.fh'

      integer iChildren(nchild)
      real    xx, yy, zz

      idummy = iCL(icell)
      if ( idummy .eq. nil ) return

      do ic1 = 1 , nchild
        iChildren(ic1) = iCh(icell,ic1)
      enddo

      call Ps ( icell , Posx,Posy,Posz )
      do while ( idummy .ne. nil )  ! construct linked lists for children
        iPL(idummy) = iPL(idummy) + 1
        xx = sngl(x(idummy))
        yy = sngl(y(idummy))
        zz = sngl(z(idummy))
        if ( xx .eq. Posx ) xx = xx + dsmall
        if ( yy .eq. Posy ) yy = yy + dsmall
        if ( zz .eq. Posz ) zz = zz + dsmall
        ifl = sign ( one , (xx - Posx) )
        jfl = sign ( one , (yy - Posy) )
        kfl = sign ( one , (zz - Posz) )
        iChild  = (ifl + 1)/2 + (jfl + 1) + 2 * (kfl + 1) + 1
        idummy2 = iLL(idummy,1)
        call LL_Insert ( iChildren(iChild) , idummy )
        idummy = idummy2
      enddo
      iCL(icell) = nil       ! cancel linked list for LL
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
      include 'a_dnew.fh'

      dimension iNonEmpty(nchild)  ! array for non-empty children
      integer ne 
c.... error if icell is a leaf
      if ( iOctCh(icell) .eq. nil ) then
        call Open_ASCII_File ( iOErr , errorfile , sequent )
        write(iOErr,*) '1  *LL_Join error: cell is a leaf'
        close(iOErr)
        stop
      endif
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
      if ( ne .gt. nil ) then 
        iCL(icell) = iCL(iNonEmpty(1))
        do ic1 = 1 , ne - 1
          idummy1 = iCL(iNonEmpty(ic1))
          do while ( idummy1 .ne. nil ) 
            iPL(idummy1) = iPL(idummy1) - 1 
            idummy2 = idummy1 
            idummy1 = iLL(idummy1,1)
          enddo
          iLL(idummy2,1) = iCL(iNonEmpty(ic1+1))
          iLL(iLL(idummy2,1),2) = idummy2
        enddo
        do ic1 = 1 , ne
          iCL(iNonEmpty(ic1)) = nil
        enddo
      endif
      return
      end

c     --------------------------
      subroutine LL_Construct ()
c     --------------------------
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_dnew.fh'
      double precision xx, yy, zz
c
      do ip = 1 , np 
        xx = x(ip)
        yy = y(ip)
        zz = z(ip)
        idc = iFindCell ( MaxLevel , xx , yy , zz )
        if ( (idc .le. nil) .or. (idc .gt. mcell) ) then 
          call Open_ASCII_File ( iOErr , errorfile , sequent ) 
          write(iOErr,*) '1 *LL_Create error: bad idc =',idc, xx, yy, zz
          close(iOErr)
          stop
        endif
        call LL_Insert ( idc , ip )
        iPL(ip) = iLv(idc)
      enddo
c      
      return
      end
c
c     ------------------------------------
      integer function iWhichSpecie ( ip ) 
c     ------------------------------------
c
      include 'a_tree.h'
c
      integer ip, isd, is
c
      isd = 0
      is  = 1000
      do while ( is .gt. 0 )  
        isd = isd + 1
        is  = (ip-1) / lsp(isd)
      enddo
c
      iWhichSpecie = isd
c
      return
      end



