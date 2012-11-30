c-------------------------------------------------------------------- 
      SUBROUTINE QSOline(xc,yc,zc,Vxc,Vyc,Vzc,
     & b,x0,y0,z0,Dx,Dy,Dz,filename)
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"
      character*90 filename
      real*8 xx,yy,zz,Tline
      real*8 d1, d0
      real*8 rs, rhogl, T_g, Z_met, xHII, xHI
! Units:
      d0 = dble(zero)
      d1 = dble(one)
      rs = max ( 1.d0 / aexpn - 1.0 , 0.d0 ) 
      aNg1 = float(Ng-1)
      gamma1 = gamma - 1.0 
      v0P = v0 / aexpn ! in km/s (proper units)
      Zsolar = 0.0199 
      a3 = aexpn**3
      a2 = aexpn**2
      ai = 1.0 / aexpn
      a2i = 1.0 / a2
      r0P = r0 / hubble * aexpn ! in Mpc (Physical units)
      r0Pk= r0*1000.0 / hubble * aexpn ! in Kpc (Physical units)
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      fact_nH = log10(1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3) ! factor converting gas density in code units to n_H in cm^-3
! Write:
      nfpn = index ( filename , ' ' ) - 1
      open ( 12 , file = filename(1:nfpn) )
      write (12,6) (xc-1.0)*r0P , (yc-1.0)*r0P, (zc-1.0)*r0P,
     & Vxc*v0P ,Vyc*v0P , Vzc*v0P 
      write (12,7)  b*r0Pk 
      write (12,8) (x0-1.0)*r0P, (y0-1.0)*r0P, (z0-1.0)*r0P
     &     , Dx,Dy,  Dz
      write (12,9) 
!     initialize cell volumes
      do iLevel=MinLevel,MaxLevel
         CellVol(iLevel)= 2.0**(-3.0*iLevel) 
      enddo
! Choice 1: Constant steps:
      call Get_MaxLevelNow () 
      iLevel=MaxLevelNow
      Tline= 2.0**(-iLevel) / 10 ! 4. 
      write (*,*)  'Choice 1: Constant steps:', Tline
      i=0
      xx= x0 ! initial values
      yy= y0
      zz= z0
      ICold = -10
      DO WHILE (xx.le.aNg1.AND.yy.le.aNg1.AND.zz.le.aNg1)

         iC = iFindCell ( MaxLevel , xx , yy , zz )
         if (IC.ne.ICold) then
            iL = iLv ( iC )
            call  Ps ( iC , Posx,Posy,Posz )
            rho = hvar(1,iC)
            Ztotal =( hvar(izII,iC) + hvar(izIa,iC) )/rho/Zsolar
            vpx = hvar(3,iC)/rho
            vpy = hvar(4,iC)/rho
            vpz = hvar(5,iC)/rho
            if ( rho .gt. 0. ) then
               T_code = gamma1 * hvar(8,iC)/rho
               ptemp = T_0 * T_code * a2i
               pent  = S_0 * T_code / rho**gamma1
            else
               write(*,*) 'bad gas density in cell',iC,rho
            endif
c
            T_g=dble(ptemp)
            rhogl=dble( log10(rho) + fact_nH)
            Z_met =dble( log10(Ztotal))
            call  X_ion ( rhogl , T_g , Z_met , rs , xHII )
            xHI =max( d0, d1 - xHII)
            
c
            write(12,17) 
     &           iL,IC, (Posx-1.0)*r0P, (Posy-1.0)*r0P ,(Posz-1.0)*r0P,
     &           vpx*v0P, vpy*v0P, vpz*v0P ,CellVol(iL)*r0Pk**3,
     &           rho*rho0Cm , Ztotal, ptemp, xHI
            ICold = IC
         endif

! Next step:
         i=i + 1
         xx = x0 + Dx*i*Tline
         yy = y0 + Dy*i*Tline
         zz = z0 + Dz*i*Tline
      enddo
 6    format ('Galaxy center in Mpc:', 3(TR3,F14.6)
     & , TR5, 'Peculiar velocity center in km/s:', 3(TR3,F10.3))
 7    format ('Impact parameter in Kpc:', F14.6 )
 8    format ('los: r=r0+Dt: r0:'3(TR3,F14.6)'  D:',3(TR3,F14.6))
 9    format('Level',TR5,'Cell #',TR5, 'X(Mpc)',TR5,'Y(Mpc)'
     & ,TR5,'Z(Mpc)'
     & ,TR5,'Vx(Km/s)',TR5,'Vy(Km/s)',TR5,'Vz(Km/s)'
     & ,TR3,'V(Kpc3)',TR3,'NH(cm-3)',TR3,'Met(solar)',TR3,'T(K)',
     & ,TR3,'HI.Frac')
 17   format (I3,2x,I9,3(TR3,F9.6),3(TR3,F10.3),4(TR3,g7.2),TR3,E10.2 )
      write (*,*)  'Output file: ', filename
      end
