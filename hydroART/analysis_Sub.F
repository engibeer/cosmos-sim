c   Subroutines for analysis
C
CEVERINO07062005
C
c-------------------------------------------------------------------- 
      SUBROUTINE ProfilesSp(xc,yc,zc,T_gas,rmin, rmax,nrbin,fprofname)
c-------------------------------------------------------------------- 
c     Make spherical profiles for gas properties and all components 
c     based on Output_Average_profile in ART_Analysis.F
c
c     input in code units.
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmin, rmax,rlogmin,rlmax
      integer nrbin
      character*256 fprofname 

      integer m
      parameter ( npoints = 1000000 )
      double precision xx , yy, zz
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins)
      real dgbin(0:nmaxbins), dgibin(0:nmaxbins)
      real dstbin(0:nmaxbins) , dstibin(0:nmaxbins), nstbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), ndmp(0:nmaxbins)
      real dtibin(0:nmaxbins),  bfibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins), entbin(0:nmaxbins)
      real gm(0:nmaxbins), gmc(0:nmaxbins)
      integer nbin(0:nmaxbins)
      a3 = aexpn**3
      a2 = aexpn**2
      ai = 1.0 / aexpn
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      deltavir = 180.           ! with respect to the critical density
!     deltavir = 340. ! with respect to the background density, see Xvir()
      
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
#ifdef STARFORM
      nDMparticles=lsp(nspec-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspec)
#endif
!     initialize cell volumes
      do iLevel=MinLevel,MaxLevel
         CellVol(iLevel)= 2.0**(-3.0*iLevel) 
      enddo
! Initialize profiles: 
      do i = 0 , nmaxbins
         rbin(i) = 0. ! averaged radius inside shell
         dgbin(i) = 0. ! gas density
         ddmbin(i) = 0. ! DM density
         pbin(i) = 0. ! pressure
         tbin(i) = 0. ! temperature
         entbin(i) = 0. ! entropy
         nbin(i) = 0 ! number of virtual points inside shell 
         ndmp(i) = 0. ! number of DM particles inside shell
         
         bfibin(i)=0            ! baryonic fraction

         dgibin(i)=0.  ! Gas mass inside radius r
         ddmibin(i) = 0. ! DM mass inside radius r
         dtibin(i) = 0. ! Total mass inside radius r

         dstbin(i)=0. ! Stellar density 
         dstibin(i)=0. ! Stellar mass inside radius r

         gm(i)=0. ! Gas mass proved inside bin
         gmc(i)=0. ! Cold Gas mass proved inside bin
      enddo
!     Sample of random points inside a sphere of rmax:
      m = 28548535              !  seed
      call Get_MaxLevelNow ()
      DO i = 1 , npoints
         rr    = rmax * RANDd ( m )
         theta = pi2 * RANDd ( m ) - pi
         phi   = pi2 * RANDd ( m )
         xx = rr * cos(theta) * cos(phi) + xc
         yy = rr * cos(theta) * sin(phi) + yc
         zz = rr * sin(theta) + zc
!     Check periodic boundaries is not necesary.
         idcell = iFindCell ( MaxLevel , xx , yy , zz )
         pd= hvar(1,idcell)
         pp =  hvar(6,idcell)
         cw = pd * CellVol(iLv(idcell)) ! Cell mass

         if ( hvar(1,idcell) .gt. 0. ) then
            T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
            ptemp = T_0 * T_code * a2i
            pent  = S_0 * T_code / hvar(1,idcell)**gamma1
         else
            write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
         endif

         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            nbin(ibin)   = nbin(ibin) + 1
            rbin(ibin)   = rbin(ibin) + rp
            dgbin(ibin)  = dgbin(ibin) + pd
            pbin(ibin) = pbin(ibin) + pp*cw
            tbin(ibin) = tbin(ibin) + ptemp*cw
            entbin(ibin) = entbin(ibin) + pent*cw
            gm(ibin) = gm(ibin) + cw
            IF (ptemp.lt.T_gas) then
               gmc(ibin) = gmc(ibin) + cw
            ENDIF
         endif
      ENDDO
      DO ip=1,nDMparticles ! DM
         xx = x(ip)
         yy = y (ip)
         zz = z (ip)
         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            ndmp(ibin) = ndmp(ibin) + 1
            ddmbin(ibin) = ddmbin(ibin) + pw(ip)
         endif
      ENDDO
#ifdef STARFORM
      DO ip=lsp(nspec-1)+1 , lsp(nspec)  ! Stars
         xx = x(ip)
         yy = y (ip)
         zz = z (ip)
         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            nstbin(ibin) = nstbin(ibin) + 1
            dstbin(ibin) = dstbin(ibin) + pw(ip)
         endif
      ENDDO
#endif
!     Normalize profiles: Mass weighted profiles for pressure, Temperature and entropy
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rnp = 1.d0 / gm(i)
            rbin(i) = rbin(i) / nbin(i) *r0P
            dgbin(i) = dgbin(i) / nbin(i)
            pbin(i) = pbin(i) * rnp 
            tbin(i) = tbin(i) *rnp
            entbin(i) = entbin(i) *rnp
         endif
      enddo
      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.                  ! DM mass inside radius
      dsi =0.                   ! stellar mass inside radius
      dgi =0.                   ! Gas mass inside radius
      dti = 0.                  ! Total mass inside radius
      DO i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlogmin + float(i-1)*drl)
         endif
         rr = 10.**(rlogmin + float(i)*drl)
         volr = 4.18879 * (rr**3 - rl**3)    
         if ( volr .gt. 0. ) then
            dmi = dmi + ddmbin(i) 
            dsi = dsi + dStbin(i)
            dgi = dgi + dgbin(i)*volr
            dti = dti + ddmbin(i) + dgbin(i)*volr + dStbin(i)

            ddmbin(i) = ddmbin(i)/ volr
            dStbin(i)  = dStbin(i) / volr
 
!            ddmibin(i) = dmi  / (4.18879 * rr**3) 
            ddmibin(i) = dmi 
            dStibin(i) = dsi  
            dgibin(i) = dgi  
            dtibin(i) = dti  / (4.18879 * rr**3)

           bfibin(i) = (dgi*(gmc(i)/gm(i))+dsi) / dti  
!	bfibin(i) = dsi / dti   ! stellar to total mass ratio       
            if ( i .gt. 0 .and.
     &           dtibin(i) .lt. deltavir .and.
     &           dtibin(i-1) .ge. deltavir ) then 
               ivir = i
               rvir = (deltavir * (rr - rl) + rl*dtibin(i) - 
     &              rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
c     rvir = rr 
            endif
            if ( i .gt. 0 .and.
     &           ddmibin(i) .lt. deltavir .and.
     &           ddmibin(i-1) .ge. deltavir ) then
            ivirDM = i
            rvirDM = (deltavir * (rr - rl) + rl*ddmibin(i) -
     &           rr*ddmibin(i-1)) / (ddmibin(i) - ddmibin(i-1))
c     rvir = rr 
         endif  
      else
         write(*,*) 'bad volr =',volr,i, rl, rr
      endif
      ENDDO
      aMvir = deltavir * 
     &     4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      aMvirDM = deltavir *
     &     4.18879 * (rvirDM*r0)**3 * rho0 / hubble**2
      write(12,12) ivirDM, rvirDM*r0P, aMvirDM/hubble, deltavir
      
c$$$!     call Xvir(deltavir,xc,yc,zc,rmin,rmax,Rvir,aMvir,Vcmax)
c$$$!     Bin of confidence:
c$$$      do i = 1 , nbins          ! The first bin is excluded
c$$$         if (ndmp(i).GT.200) then 
c$$$            ibinconfidence=i
c$$$            goto 302
c$$$         endif
c$$$      enddo
c$$$ 302  continue
c$$$      
c$$$      
c$$$! NFW fit: Between ibinconfidence and ivirDM. Compute concentration
c$$$      c=10.                     ! initial c
c$$$      sq = 0.
c$$$      do i = ibinconfidence , ivirDM
c$$$         rr = 10.**(rlogmin + float(i)*drl)
c$$$         rrvir = rr / rvirDM
c$$$         sq = sq + ( ddmbin(i) - fitNFW(deltavir,rrvir,c) )**2
c$$$      enddo
c$$$      i_LSF=1
c$$$ 301  continue
c$$$      i_LSF = i_LSF + 1
c$$$      sqold=sq
c$$$      sq=0.
c$$$      c = c + 1.
c$$$      DO i = ibinconfidence , ivirDM
c$$$         rr = 10.**(rlogmin + float(i)*drl)
c$$$         rrvir = rr / rvirDM
c$$$         sq = sq + ( ddmbin(i) - fitNFW(deltavir,rrvir,c) )**2
c$$$      enddo
c$$$      
c$$$      if (sq.lt.sqold)  GOTO 301
c$$$      c = c -1.
!     Write outputfile:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0  , (yc-1.)*r0 , (zc-1.)*r0    
      write(12,11) rmin*r0P, rmax*r0P , nrbin 
      write(12,12) ivir, rvir*r0P, aMvir/hubble, deltavir
      write (12,13) c
      write(12,14) aexpn , 1/aexpn -1. , T_gas
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[Mpc h-1 comoving]= ', 3(F8.3,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5)
 12   format ('# ivir = ',i4,' Rvir =',f8.1,
     &     ' [Kpc], Mvir =',g13.6,
     &     ' [Msun] for Dvir(mean)=',f7.2)
 13   format ('NFW fit: concentration= ', f3.0)
 14   format
     &('Distances in Kpc ,densities in cm-3 and M in M_Sun,
     &  for a= ',f6.3, ' and redshift= ', f6.2, 
     &' Baryonic fraction of cold gas', E7.1, ' K ')
 15   format ('# <r>',14x,' r_l',14x,' r_m',14x,' r_r',14x,
     &     ' Rho_g ',3x, 'Rho_St ',8x, 'Rho_dm',10x,
     &     'M_insideR_Gas',3x,  'M_insideR_Stars',3x,
     &     ' M_insideR_DM ',3x,' M_insideR_tot',4x,
     &     'T_g',10x,'S_g',8x,
     &     'N_dm_Pibin',2x,'Baryonic_F  [Uni]' )
 16   format('---------------------')
!     16   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
!     &        '             K  keV cm^2 ')
      NdmPibin=0.
      do i = 0 , nbins
         NdmPibin =NdmPibin + ndmp(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin 
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(12,17)
     &           rbin(i), rl*r0P  , 
     &         rmid*r0P , rr*r0P  , 
     &           dgbin(i)*rho0Cm ,dStbin(i)*rho0Cm, ddmbin(i)*rho0Cm, 
     &           dgibin(i)*aM0   ,
     &           dstibin(i)*aM0    ,
     &           ddmibin(i)*aM0 ,
     &           dtibin(i)*4.18879*(rr*r0)**3*rho0/hubble**3 ,  
     &           tbin(i),entbin(i),int(NdmPibin), bfibin(i)/(Omb0/Om0)
!    dgibin(i)*4.18879*(rr*r0)**3*rho0/hubble**3
         endif
      enddo
 17   format(13(2x,g14.7),1x,i7,1x,g14.7)
      close ( 12 )
      write(*,*) 'done writing spherical profiles. exiting...'
      return
      end
c-------------------------------------------------------------------- 
      SUBROUTINE CheckVirP(xc,yc,zc,rmin, rmax,nrbin,fprofname)
c-------------------------------------------------------------------- 
c
c     Check the distribution of virtual points used to make gas profiles
c
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmin, rmax,rlogmin,rlmax
      integer nrbin
      character*256 fprofname 

      integer m
      parameter ( npoints = 1000000 )
      double precision xx , yy, zz
      parameter ( nmaxbins = 1000 )


      DO i = 1 , npoints
         rr    = rmax * RANDd ( m )
         theta = pi2 * RANDd ( m ) - pi
         phi   = pi2 * RANDd ( m )
         xx = rr * cos(theta) * cos(phi) + xc
         yy = rr * cos(theta) * sin(phi) + yc
         zz = rr * sin(theta) + zc

      ENDDO
      END


c$$$c-------------------------------------------------------------------- 
c$$$      SUBROUTINE ProfilesStars (xc,yc,zc,rmin, rmax,nrbin,fprofname)
c$$$c-------------------------------------------------------------------- 
c$$$c     Make profiles for stars
c$$$c
c$$$c     input in code units.
c$$$c-------------------------------------------------------------------- 
c$$$      include 'a_tree.h'
c$$$      include 'a_control.h'
c$$$#    include "a_def.h"
c$$$
c$$$      double precision xx , yy, zz
c$$$      real xc,yc,zc,rmin, rmax,rlogmin,rlmax
c$$$      parameter ( nmaxbins = 1000 )
c$$$      integer nstbin(0:nmaxbins)
c$$$      real dstbin(0:nmaxbins)
c$$$      real dstibin(0:nmaxbins)
c$$$      real rbin(0:nmaxbins)
c$$$
c$$$      a3 = aexpn**3
c$$$      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
c$$$      r0P = r0 *1000.0 * hubble * aexpn ! in Kpc (Physical units)
c$$$      rlogmin = log10(rmin)
c$$$      rlmax = log10(rmax)
c$$$      drl = (rlmax - rlogmin) / float(nrbin) 
c$$$      nbins = int((rlmax-rlogmin)/drl) + 1
c$$$! Initialize profiles: 
c$$$      do i = 0 , nmaxbins
c$$$         rbin(i)=0.
c$$$         nstbin(i)=0
c$$$         dstbin(i)=0.
c$$$         dstibin(i) =0.
c$$$      enddo
c$$$      DO ip=lsp(nspec-1)+1 , lsp(nspec)
c$$$         xx = x(ip)
c$$$         yy = y (ip)
c$$$         zz = z (ip)
c$$$         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
c$$$         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
c$$$         if ( ibin .le. nmaxbins ) then
c$$$            nstbin(ibin) = nstbin(ibin) + 1
c$$$            dstbin(ibin) = dstbin(ibin) + pw(ip)
c$$$            rbin(ibin) = rbin(ibin) + rp
c$$$         endif
c$$$      ENDDO
c$$$!     Normalize profiles:
c$$$      do i = 0 , nbins
c$$$         if ( nbin(i) .gt. 0 ) then
c$$$            rbin(i) = rbin(i) / nstbin(i) *r0P
c$$$         enddo
c$$$         dsti =0.               ! Mass in stars inside radius
c$$$      DO i = 0 , nbins
c$$$         if ( i .eq. 0 ) then 
c$$$            rl = 0.0
c$$$         else
c$$$            rl = 10.**(rlogmin + float(i-1)*drl)
c$$$         endif
c$$$         rr = 10.**(rlogmin + float(i)*drl)
c$$$         volr = 4.18879 * (rr**3 - rl**3)        
c$$$         if ( volr .gt. 0. ) then
c$$$            dsti = dsti +  dstbin(i)
c$$$            dstibin(i) = dsti
c$$$            dstbin(i) =    dstbin(i) / volr
c$$$         else
c$$$            write(*,*) 'bad volr =',volr,i, rl, rr
c$$$         endif
c$$$      ENDDO
c$$$!     Write outputfile:
c$$$      nfpn = index ( fprofname , ' ' ) - 1
c$$$      open ( 12 , file = fprofname(1:nfpn) )
c$$$      write(12,10) (xc-1.)*r0  , (yc-1.)*r0 , (zc-1.)*r0    
c$$$      write(12,11) rmin*r0P, rmax*r0P , nrbin 
c$$$      write(12,14) aexpn , 1/aexpn -1.
c$$$      write(12,15)
c$$$      write(12,16)
c$$$ 10   format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
c$$$ 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
c$$$     &     ' nrbin =',i5)
c$$$ 14   format
c$$$     &('Distances in Kpc ,densities in Msun/pc and M in M_Sun,
c$$$     &  for a= ',f6.3, ' and redshift= ', f6.2)
c$$$ 15   format ('# <r>',14x,' r_l',14x,' r_m',14x,' r_r',14x,
c$$$     &     ' Rho_stars ',3x,' M_insideRR_ST ',3x,
c$$$     &     'N_St_Pibin' )
c$$$ 16   format('---------------------')
c$$$      NstPibin=0.
c$$$      do i = 0 , nbins
c$$$         NstPibin =NstPibin + nstbin(i)
c$$$         if ( i .eq. 0 ) then 
c$$$            rl = 0.0  
c$$$            rmid = 0.5*10.**rlogmin 
c$$$         else
c$$$            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
c$$$            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
c$$$         endif
c$$$         rr = 10.**(rlogmin + float(i)*drl) ! rightside
c$$$         if ( nbin(i) .gt. 0 ) then
c$$$            write(12,17)
c$$$     &           rbin(i), rl*r0P  , 
c$$$     &         rmid*r0P , rr*r0P  , 
c$$$     &         dstbin(i)*rhoC , dstibin(i),
c$$$     &         NstPibin
c$$$         endif
c$$$      enddo
c$$$ 17   format (6(2x,g14.7),1x,i7)
c$$$      close ( 12 )
c$$$      write(*,*) 'done writing profiles for stars. exiting...'
c$$$      return
c$$$      end
c-------------------------------------------------------------------- 
c------------------------------------------------------------------ 
      SUBROUTINE Vcenter (npoints,xc,yc,zc, rmax, rmin, vcx, vcy, vcz)
c-------------------------------------------------------------------- 
c     Mean velocity inside spheres centered in xc,yc,zc
c 
c-------------------------------------------------------------------- 
       include 'a_tree.h'
      include 'a_control.h'

      integer npoints
      real rmax, rmin
      real xc,yc,zc, vcx, vcy, vcz
      double precision xx , yy, zz
      real*8 x8, y8, z8, v8,v8x, v8y, v8z
      call Get_MaxLevelNow ()
!     initialize cell volumes
      do iLevel=MinLevel,MaxLevel
         CellVol(iLevel)= 2.0**(-3.0*iLevel) 
      enddo
!     Sample of random points inside a sphere of rmax:
      Rselect = rmax
      Ncount=0
 300  CONTINUE
      v8x=0.
      v8y =0.
      v8z =0.
      x8 = 0.
      y8 = 0.
      z8 = 0.
      amass=0.
      m = 28548535              !  seed
      DO i = 1 , npoints
         rr    = Rselect * RANDd ( m )
         theta = pi2 * RANDd ( m ) - pi
         phi   = pi2 * RANDd ( m )
         xx = rr * cos(theta) * cos(phi) + xc
         yy = rr * cos(theta) * sin(phi) + yc
         zz = rr * sin(theta) + zc 
!     Check periodic boundaires is not necesary.
         idcell = iFindCell ( MaxLevel , xx , yy , zz )
         pd = hvar(1,idcell)
         vpx = hvar(3,idcell)/hvar(1,idcell)
         vpy = hvar(4,idcell)/hvar(1,idcell)
         vpz = hvar(5,idcell)/hvar(1,idcell)
          cw = pd * CellVol(iLv(idcell)) ! Cell mass
         if ( vpx.eq.0.) write(*,*) 'bad cell', idcell, vpx
         v8x = v8x + vpx*cw
         v8y = v8y + vpy*cw
         v8z = v8z + vpz*cw
         x8 = x8 + xx*cw
         y8 = y8 + yy*cw
         z8 = z8 + zz*cw
         amass = amass +cw
      ENDDO
      vcx = v8x / amass
      vcy = v8y / amass
      vcz = v8z / amass
      xc = x8 / amass
      yc = y8 / amass
      zc = z8 / amass
c      write (*,12)  Rselect*1000.0*aexpn*r0*hubble ,
c     &  (xc-1.)*1000.0*aexpn*r0*hubble,
c     &  (yc-1.)*1000.0*aexpn*r0*hubble,  
c     &   (zc-1.)*1000.0*aexpn*r0*hubble,
c     $ vcx * v0 , vcy * v0 , vcz * v0 
      Ncount = Ncount + 1
      Rselect = Rselect / 1.1
      if ( Ncount.lt.100.AND.Rselect.gt.rmin) GOTO 300
      return
 12   format ( 7(f10.3))
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE ProfilesVGa (vcx, vcy, vcz, xc,yc,zc,
     &     T_gas, rmin, rmax,nrbin,fprofname)
c-------------------------------------------------------------------- 
c     Make velocity profiles for gas in spherical coordinates
c
c     input in code units. 
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'

      real vcx, vcy, vcz,xc,yc,zc,rmin, rmax
      integer nrbin
      character*256 fprofname 
      integer m
      parameter ( npoints = 1000000 )
      double precision xx , yy, zz
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),vbin(0:nmaxbins) 
      real vrbin(0:nmaxbins), vthbin(0:nmaxbins), vphbin(0:nmaxbins)
      integer nbin(0:nmaxbins)
      real vxbin(0:nmaxbins) , vybin(0:nmaxbins), vzbin(0:nmaxbins)
      real rmsVxbin(0:nmaxbins), vx2Sum(0:nmaxbins)
      real rmsVybin(0:nmaxbins), vy2Sum(0:nmaxbins)
      real rmsVzbin(0:nmaxbins), vz2Sum(0:nmaxbins)
      real rmsVbin(0:nmaxbins), v2Sum(0:nmaxbins)
      real rmsVrbin(0:nmaxbins), vr2Sum(0:nmaxbins)
      real gm(0:nmaxbins)
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      v0P = v0 / aexpn ! in km/s (proper units)
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
!     initialize cell volumes
      do iLevel=MinLevel,MaxLevel
         CellVol(iLevel)= 2.0**(-3.0*iLevel) 
      enddo
! Initialize profiles:
      do i = 0 , nmaxbins
         rbin(i) = 0.
         nbin(i) = 0
         vbin(i)=0.
         vrbin(i)=0.
         vthbin(i)=0.
         vphbin(i)=0.
         
         vxbin(i)=0.
         vybin(i)=0.
         vzbin(i)=0.
         
         vx2sum(i)=0.
         rmsVxbin(i)=0.
         vy2sum(i)=0.
         rmsVybin(i)=0.
         vz2sum(i)=0.
         rmsVzbin(i)=0.

         v2sum(i) = 0.
         rmsVbin(i) =0.
         vr2sum(i) = 0.
         rmsVrbin(i) =0.

         gm(i)=0.
      enddo
!      call Vcenter(npoints,xc,yc,zc, rmax, rmin, vcx, vcy, vcz)
!      write (*,*) xc, yc, zc, vcx, vcy, vcz

! Sample of random points inside a sphere of rmax:
      m = 28548535 !  seed
      rlogmin = log10(rmin)
      call Get_MaxLevelNow ()
      DO i = 1 , npoints
         rr    = rmax * RANDd ( m )
         theta = pi2 * RANDd ( m ) - pi
         phi   = pi2 * RANDd ( m )
         costheta = cos(theta)
         sintheta = sin(theta)
         cosphi =  cos(phi)
         sinphi = sin(phi)
         xx = rr * costheta * cosphi + xc
         yy = rr * costheta * sinphi + yc
         zz = rr * sintheta + zc
!     Check periodic boundaires is not necesary.
         idcell = iFindCell ( MaxLevel , xx , yy , zz )
         if ( hvar(1,idcell) .gt. 0. ) then
            T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
            ptemp = T_0 * T_code * a2i
         else
            write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
         endif
         IF (ptemp.lt.T_gas) then
            cw =hvar(1,idcell) * CellVol(iLv(idcell)) ! Cell mass
            vpx= hvar(3,idcell)/hvar(1,idcell) - vcx
            vpy =  hvar(4,idcell)/hvar(1,idcell) - vcy
            vpz =  hvar(5,idcell)/hvar(1,idcell) - vcz
            vpx2 = vpx**2
            vpy2 = vpy**2
            vpz2 = vpz**2
            vp2 = vpx**2 + vpy**2 + vpz**2
            vp = sqrt (vp2)
            vpr =vpx*cosphi*costheta + vpy*sinphi*costheta +vpz*sintheta
            vpth=vpx*cosphi*sintheta +vpy*sinphi*sintheta - vpz*costheta
            vpph =-vpx*sinphi + vpy*cosphi
            vpr2 = vpr**2
        
            rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
            ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
            if ( ibin .le. nmaxbins ) then
               rbin(ibin) = rbin(ibin) + rp
               nbin(ibin) = nbin(ibin) + 1
               vbin(ibin)= vbin(ibin) + vp*cw
               vrbin(ibin)= vrbin(ibin) + vpr*cw
               vthbin(ibin)= vthbin(ibin) + vpth*cw
               vphbin(ibin) = vphbin(ibin) + vpph*cw
               
               vxbin(ibin) = vxbin(ibin) + vpx*cw
               vybin(ibin) = vybin(ibin) + vpy*cw
               vzbin(ibin) = vzbin(ibin) + vpz*cw

               vx2Sum(ibin) =  vx2Sum(ibin) + Vpx2*cw
               vy2sum(ibin) = vy2sum(ibin) + Vpy2*cw
               vz2sum(ibin) =   vz2sum(ibin) + Vpz2*cw
               v2sum(ibin) = v2sum(ibin) + Vp2*cw
               vr2sum(ibin) =  vr2sum(ibin) + vpr2*cw
            
               gm(ibin)= gm(ibin) + cw
            endif
         ENDIF
      ENDDO
! rms velocities without normalization:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rmsVxbin(i) = vx2Sum(i) - vxbin(i)**2 / gm(i) ! nbin(i) 
            rmsVybin(i) = vy2Sum(i) - vybin(i)**2 / gm(i)
            rmsVzbin(i) = vz2Sum(i) - vzbin(i)**2 / gm(i)
            rmsVbin(i)  = v2Sum(i)   -  vbin(i)**2 / gm(i)
            rmsVrbin(i) = vr2Sum(i) - vrbin(i)**2 / gm(i)
         endif
      enddo
!     Normalize profiles:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rbin(i) = rbin(i)/  nbin(i) *r0P
            vbin(i)= vbin(i) / gm(i)
            vrbin(i)= vrbin(i) / gm(i)
            vthbin(i)= vthbin(i) / gm(i)
            vphbin(i) = vphbin(i) / gm(i)

            vxbin(i) = vxbin(i) / gm(i)
            vybin(i) = vybin(i) / gm(i)
            vzbin(i) = vzbin(i) / gm(i)

            if (rmsVxbin(i).lt.0.)  then
               if (abs(rmsVxbin(i)).lt.100.)then
                  rmsVxbin(i)=0.0
               else
                  write (*,*)  'rmsVxbin(i)'
                  write (*,*)  rmsVxbin(i)
                  STOP
               endif
            endif
            if (rmsVybin(i).lt.0.)  then
               if (abs(rmsVybin(i)).lt.100.) then
                  rmsVybin(i)=0.0
               else
                  write (*,*)  'rmsVybin(i)'
                  write (*,*)  rmsVybin(i)
                  STOP
               endif
            endif
            if (rmsVzbin(i).lt.0.)  then
               if (abs(rmsVzbin(i)).lt.100.) then 
                  rmsVzbin(i)=0.0
               else
                  write (*,*)  'rmsVzbin(i)'
                  write (*,*)  rmsVzbin(i)
                  STOP
               endif
            endif
            if (rmsVbin(i).lt.0.)  then
               if (abs(rmsVbin(i)).lt.100.) then
                  rmsVbin(i)=0.0
               else
                  write (*,*)  'rmsVbin(i)'
                  write (*,*)  rmsVbin(i)
                  STOP
               endif
            endif
            if (rmsVrbin(i).lt.0.)  then
               if (abs(rmsVrbin(i)).lt.100.) then
                  rmsVrbin(i)=0.0
               else
                  write (*,*) ' rmsVrbin(i)'
                  write (*,*)  rmsVrbin(i)
                  STOP
               endif
            endif

            rmsVxbin(i) = sqrt( rmsVxbin(i) / gm(i) )
            rmsVybin(i) = sqrt( rmsVybin(i) / gm(i) )
            rmsVzbin(i) = sqrt( rmsVzbin(i) / gm(i) )
            rmsVbin(i)   = sqrt( rmsVbin(i) / gm(i) )
            rmsVrbin(i) = sqrt( rmsVrbin(i) / gm(i) )
         endif
      enddo
!     Write output:
         nfpn = index ( fprofname , ' ' ) - 1
         open ( 12 , file = fprofname(1:nfpn) )
         write(12,10) (xc-1.)*r0 , (yc-1.)*r0 , (zc-1.)*r0   
         write(12,11) vcx*v0 , vcy*v0, vcz*v0
         write(12,12) rmin*r0P, rmax*r0P , nrbin 
         write (12,13)
         write(12,14) T_gas     
         write(12,15)
         write(12,16)
 10      format ('# (xc,yc,zc)[/h Mpc comoving] = ', 3(g13.7,1x))
 11      format ('# (Vcx,Vcy,Vcz)[Km/s comoving] = ', 3(g13.7,1x))
 12      format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5)
 13      format ('GAS velocity mass weighted profiles ')
 14      format (' with temperature bellow ', E7.1)
 15      format ('# <r>',14x,'r_l',14x,' r_m',14x,'r_r',14x,
     &        'v',7x,'vr',7x,'vtheta',7x,'vphi',7x,
     &        'vx',7x,' vy',7x,'vz',7x,
     &        'rmsVx',5x,'rmsVy ',5x,'rmsVz',
     &        5x,'rmsV',5x,'rmsVr',2x,'NVPibin')
 16      format('---------------------------------------')
         Npibin=0.
         do i = 0 , nbins
            Npibin =Npibin + nbin(i)
            if ( i .eq. 0 ) then 
               rl = 0.0  
               rmid = 0.5*10.**rlogmin
            else
               rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
               rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
            endif
            rr = 10.**(rlogmin + float(i)*drl) ! rightside
            if ( nbin(i) .gt. 0 ) then
               write(12,17)
     &              rbin(i), rl*r0P, 
     &              rmid*r0P, rr*r0P, 
     &           vbin(i)*v0P, vrbin(i)*v0P ,vthbin(i)*v0P,vphbin(i)*v0P,
     &              vxbin(i)*v0P , vybin(i)*v0P , vzbin(i)*v0P,
     &              rmsVxbin(i)*v0P , rmsVybin(i)*v0P , rmsVzbin(i)*v0P,
     &              rmsVbin(i)*v0P , rmsVrbin(i)*v0P,
     &              npibin
            endif
         enddo
 17      format (4(2x,g14.7), 12(2x,f8.3),  1x, I6)
         close ( 12 )
         write(*,*) 'done writing velocity profiles for gas. exiting...'
         return
         END
c-------------------------------------------------------------------- 
      SUBROUTINE ProfilesVdm (vcx, vcy, vcz,xc,yc,zc
     &     ,rmin, rmax,nrbin,fprofname)
c-------------------------------------------------------------------- 
c     Make velocity profiles for DARK MATTER 
c--------------------------------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real vcx, vcy, vcz,xc,yc,zc,rmin, rmax
      integer nrbin
      character*256 fprofname 
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),vbin(0:nmaxbins) 
      real vrbin(0:nmaxbins), vthbin(0:nmaxbins), vphbin(0:nmaxbins)
      integer nbin(0:nmaxbins)
      real vxbin(0:nmaxbins) , vybin(0:nmaxbins), vzbin(0:nmaxbins)
      real rmsVxbin(0:nmaxbins),rmsVybin(0:nmaxbins),
     &     rmsVzbin(0:nmaxbins)
      real rmsVbin(0:nmaxbins) 
      real vr2sum(0:nmaxbins) , rmsVrbin(0:nmaxbins)

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      v0P = v0 / aexpn ! in km/s (proper units)
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
#ifdef STARFORM
      nDMparticles=lsp(nspec-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspec)
#endif
!     Initialize profiles:
      do i = 0 , nmaxbins
         rbin(i) = 0.
         nbin(i) = 0
         vbin(i)=0.
         vrbin(i)=0.
         vthbin(i)=0.
         vphbin(i)=0.
         
         vxbin(i)=0.
         vybin(i)=0.
         vzbin(i)=0.

         rmsVxbin(i)=0.
         rmsVybin(i)=0.
         rmsVzbin(i)=0.

         rmsVbin(i)=0.
         rmsVrbin(i)=0.

         vr2sum(i)=0.
      enddo
     
!      call DMcenter (rmin, rmax,Xc , Yc , Zc, Vcx, Vcy, Vcz)
!      write (*,*) xc, yc, zc, vcx, vcy, vcz
      DO ip=1,nDMparticles
         xp = x(ip) -xc
         yp = y (ip) - yc
         zp = z (ip) -zc
         rp = sqrt((xp)**2+(yp)**2+(zp)**2)
         vpx = vx(ip) - vcx
         vpy = vy(ip) - vcy
         vpz = vz(ip) - vcz
         vp = sqrt (vpx**2 + vpy**2 + vpz**2)
         vpr = ( vpx*xp + vpy*yp + vpz*zp)/rp
         vpr2=vpr**2
         if ( xp.eq.0.) then
            PhiP=0.
         else
            PhiP=atan(Yp/Xp)
            if (Xp.gt.0.) then
               if(Yp.lt.0.) then
                  PhiP =PhiP + 2.*pi
               endif
            else
               PhiP=   PhiP + pi
            endif   
         endif
         thetaP = asin (zp/rp)
         vpth=vpx*cos(PhiP)*sin(thetaP) + vpy*sin(PhiP)*Sin(thetaP)
     $        -vpz*cos(thetaP)
         vpph=-vpx*sin(PhiP) + vpy*cos(PhiP)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            rbin(ibin) = rbin(ibin) + rp
            nbin(ibin) = nbin(ibin) + 1
            vbin(ibin)= vbin(ibin) + vp
            vrbin(ibin)= vrbin(ibin) + vpr
            vthbin(ibin)= vthbin(ibin) + vpth
            vphbin(ibin) = vphbin(ibin) + vpph
            
            vxbin(ibin) = vxbin(ibin) + vpx
            vybin(ibin) = vybin(ibin) + vpy
            vzbin(ibin) = vzbin(ibin) + vpz

            vr2sum(ibin) =    vr2sum(ibin) + vpr2
         endif
      ENDDO
! rms velocities without normalization:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rmsVrbin(i) = vr2Sum(i) - vrbin(i)**2 / nbin(i)
         endif
      enddo
! Normalize profiles:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rbin(i) = rbin(i)/  nbin(i) *r0P
            vbin(i)= vbin(i) / nbin(i)
            vrbin(i)= vrbin(i) / nbin(i)
            vthbin(i)= vthbin(i) / nbin(i)
            vphbin(i) = vphbin(i) / nbin(i)
            
            vxbin(i) = vxbin(i) / nbin(i)
            vybin(i) = vybin(i) / nbin(i)
            vzbin(i) = vzbin(i) / nbin(i)

           rmsVrbin(i) = sqrt( rmsVrbin(i) / nbin(i) )
        endif
      enddo
!     rms velocity:  OLD
      DO i = 1 ,nDMparticles

         rp = sqrt((x(i)-xc)**2+(y(i)-yc)**2+(z(i)-zc)**2)

         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            rmsVxbin(ibin)= rmsVxbin(ibin) +( Vx(i)-Vcx -vxbin(ibin))**2 
            rmsVybin(ibin)= rmsVybin(ibin) +( Vy(i)-Vcy -vybin(ibin))**2
            rmsVzbin(ibin)= rmsVzbin(ibin) + (Vz(i)-Vcz -vzbin(ibin))**2  
         endif
      ENDDO
      do i = 0 , nbins
         rmsVbin(i)= sqrt ( (rmsVxbin(i)+rmsVybin(i)+rmsVzbin(i))
     &        / nbin(i) )
         rmsVxbin(i) =sqrt( rmsVxbin(i)/ nbin(i))
         rmsVybin(i) =sqrt( rmsVybin(i) /  nbin(i))
         rmsVzbin(i) =sqrt( rmsVzbin(i) / nbin(i))
      enddo
!     Write output:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0, (yc-1.)*r0 , (zc-1.)*r0   
      write(12,11) vcx*v0 , vcy*v0, vcz*v0
      write(12,12) rmin*r0P, rmax*r0P , nrbin 
      write (12,13) aexpn , 1/aexpn -1.
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[/h Mpc comoving] = ', 3(g13.7,1x))
 11   format ('# (Vcx,Vcy,Vcz)[Km/s comoving] = ', 3(g13.7,1x))
 12   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5)
 13   format ('Distances in proper Kpc, proper velocities in Km/s
     &     for a= ',f6.3, ' and redshift= ', f6.2)
 14   format ('#DARK MATTER velocity profiles: Columns:')
 15   format ('# <r>',14x,'r_l',14x,' r_m',14x,'r_r',14x,
     &     'v',7x,'vr',7x,'vtheta',7x,'vphi',7x,
     &     'vx',7x,' vy',7x,'vz',7x,
     &     'rmsVx',5x,'rmsVy ',5x,'rmsVz',5x,'rmsV',5x,'rmsVr', 
     &     1x,'NdmPibin')
 16   format('----------------------')
      NdmPibin=0.
      do i = 0 , nbins
         NdmPibin =NdmPibin + nbin(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(12,17)
     &           rbin(i), rl*r0P, 
     &           rmid*r0P, rr*r0P, 
     &         vbin(i)*v0P,vrbin(i)*v0P , vthbin(i)*v0P , vphbin(i)*v0P,
     &           vxbin(i)*v0P , vybin(i)*v0P , vzbin(i)*v0P,
     &           rmsVxbin(i)*v0P , rmsVybin(i)*v0P , rmsVzbin(i)*v0P,
     &           rmsVbin(i)*v0P, rmsVrbin(i)*v0P  , ndmPibin
         endif
      enddo
 17   format (4(2x,g14.7), 12(2x,f8.3), 1x, I6)
      close ( 12 )
      write(*,*) 'done writing velocity profiles for DM. exiting...'
      return
      END
c--------------------------------------------------------------------
c--------------------------------------------------------------------
c-------------------------------------------------------------------- 
      SUBROUTINE ProfilesVSt(vcx, vcy, vcz,xc,yc,zc
     &     ,rmin, rmax,nrbin,fprofname)
c-------------------------------------------------------------------- 
c     Make velocity profiles for the stellar component 
c--------------------------------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real vcx, vcy, vcz,xc,yc,zc,rmin, rmax
      integer nrbin
      character*256 fprofname 
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),vbin(0:nmaxbins) 
      real vrbin(0:nmaxbins), vthbin(0:nmaxbins), vphbin(0:nmaxbins)
      integer nbin(0:nmaxbins)
      real vxbin(0:nmaxbins) , vybin(0:nmaxbins), vzbin(0:nmaxbins)
      real rmsVxbin(0:nmaxbins),rmsVybin(0:nmaxbins),
     &     rmsVzbin(0:nmaxbins)
      real rmsVbin(0:nmaxbins)
      real vr2sum(0:nmaxbins) , rmsVrbin(0:nmaxbins)
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      v0P = v0 / aexpn ! in km/s (proper units)
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
!     Initialize profiles:
      do i = 0 , nmaxbins
         rbin(i) = 0.
         nbin(i) = 0
         vbin(i)=0.
         vrbin(i)=0.
         vthbin(i)=0.
         vphbin(i)=0.
         
         vxbin(i)=0.
         vybin(i)=0.
         vzbin(i)=0.

         rmsVxbin(i)=0.
         rmsVybin(i)=0.
         rmsVzbin(i)=0.

         rmsVbin(i)=0.
         
         vr2sum(i)=0.
      enddo
      DO ip=nsp(nspec,1),nsp(nspec,2)  ! STARS
         xp = x(ip) -xc
         yp = y (ip) - yc
         zp = z (ip) -zc
         rp = sqrt((xp)**2+(yp)**2+(zp)**2)
         vpx = vx(ip) - vcx
         vpy = vy(ip) - vcy
         vpz = vz(ip) - vcz
         vp = sqrt (vpx**2 + vpy**2 + vpz**2)
         vpr = ( vpx*xp + vpy*yp + vpz*zp)/rp
         vpr2=vpr**2
         if ( xp.eq.0.) then
            PhiP=0.
         else
            PhiP=atan(Yp/Xp)
            if (Xp.gt.0.) then
               if(Yp.lt.0.) then
                  PhiP =PhiP + 2.*pi
               endif
            else
               PhiP=   PhiP + pi
            endif   
         endif
         thetaP = asin (zp/rp)
         vpth=vpx*cos(PhiP)*sin(thetaP) + vpy*sin(PhiP)*Sin(thetaP)
     $        -vpz*cos(thetaP)
         vpph=-vpx*sin(PhiP) + vpy*cos(PhiP)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            rbin(ibin) = rbin(ibin) + rp
            nbin(ibin) = nbin(ibin) + 1
            vbin(ibin)= vbin(ibin) + vp
            vrbin(ibin)= vrbin(ibin) + vpr
            vthbin(ibin)= vthbin(ibin) + vpth
            vphbin(ibin) = vphbin(ibin) + vpph
            
            vxbin(ibin) = vxbin(ibin) + vpx
            vybin(ibin) = vybin(ibin) + vpy
            vzbin(ibin) = vzbin(ibin) + vpz

            vr2sum(ibin) =    vr2sum(ibin) + vpr2
         endif
      ENDDO
! rms velocities without normalization:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rmsVrbin(i) = vr2Sum(i) - vrbin(i)**2 / nbin(i)
         endif
      enddo
! Normalize profiles:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rbin(i) = rbin(i)/  nbin(i) *r0P
            vbin(i)= vbin(i) / nbin(i)
            vrbin(i)= vrbin(i) / nbin(i)
            vthbin(i)= vthbin(i) / nbin(i)
            vphbin(i) = vphbin(i) / nbin(i)
            
            vxbin(i) = vxbin(i) / nbin(i)
            vybin(i) = vybin(i) / nbin(i)
            vzbin(i) = vzbin(i) / nbin(i)

            rmsVrbin(i) = sqrt( rmsVrbin(i) / nbin(i) )
         endif
      enddo
!     rms velocity: OLD
      DO i = nsp(nspec,1),nsp(nspec,2)  ! STARS
         rp = sqrt((x(i)-xc)**2+(y(i)-yc)**2+(z(i)-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            rmsVxbin(ibin)= rmsVxbin(ibin) +( Vx(i)-Vcx -vxbin(ibin))**2 
            rmsVybin(ibin)= rmsVybin(ibin) +( Vy(i)-Vcy -vybin(ibin))**2
            rmsVzbin(ibin)= rmsVzbin(ibin) + (Vz(i)-Vcz -vzbin(ibin))**2
         endif
      ENDDO
      do i = 0 , nbins
         rmsVbin(i)= sqrt ( (rmsVxbin(i)+rmsVybin(i)+rmsVzbin(i))
     &        / nbin(i) )
         rmsVxbin(i) =sqrt( rmsVxbin(i)/ nbin(i))
         rmsVybin(i) =sqrt( rmsVybin(i) /  nbin(i))
         rmsVzbin(i) =sqrt( rmsVzbin(i) / nbin(i))
      enddo
!     Write output:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0, (yc-1.)*r0 , (zc-1.)*r0   
      write(12,11) vcx*v0 , vcy*v0, vcz*v0
      write(12,12) rmin*r0P, rmax*r0P , nrbin 
      write (12,13) aexpn , 1/aexpn -1.
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[/h Mpc comoving] = ', 3(g13.7,1x))
 11   format ('# (Vcx,Vcy,Vcz)[Km/s comoving] = ', 3(g13.7,1x))
 12   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5)
 13   format ('Distances in proper Kpc, proper velocities in Km/s
     &     for a= ',f6.3, ' and redshift= ', f6.2)
 14   format ('# STELLAR VELOCITY PROFILES: Columns:')
 15   format ('# <r>',14x,'r_l',14x,' r_m',14x,'r_r',14x,
     &     'v',7x,'vr',7x,'vtheta',7x,'vphi',7x,
     &     'vx',7x,' vy',7x,'vz',7x,
     &     'rmsVx',5x,'rmsVy ',5x,'rmsVz',5x,'rmsV' ,5x,'rmsVr', 
     &     1x,'NPibin')
 16   format('----------------------')
      NdmPibin=0.
      do i = 0 , nbins
         NdmPibin =NdmPibin + nbin(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(12,17)
     &           rbin(i), rl*r0P, 
     &           rmid*r0P, rr*r0P, 
     &       vbin(i)*v0P, vrbin(i)*v0P , vthbin(i)*v0P , vphbin(i)*v0P ,
     &           vxbin(i)*v0P , vybin(i)*v0P , vzbin(i)*v0P,
     &           rmsVxbin(i)*v0P , rmsVybin(i)*v0P , rmsVzbin(i)*v0P,
     &           rmsVbin(i)*v0P,  rmsVrbin(i)*v0P ,ndmPibin
         endif
      enddo
 17   format (4(2x,g14.7), 12(2x,f8.3), 1x, I6)
      close ( 12 )
      write(*,*)'done writing of velocity profiles for stars. Exiting..'
      return
      END
c--------------------------------------------------------------------
c--------------------------------------------------------------------

      SUBROUTINE ParticlesCenter (rmin, rmax,
     &     Xdmc , Ydmc , Zdmc, Vxdmc, Vydmc, Vzdmc)
c-------------------------------------------------------------------- 
c     Find center of the particles distribution (DM+stars) 
c     based on the iterative method in PMgalaxy.f
c
c  input in code units
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real Xdmc , Ydmc , Zdmc, Vxdmc, Vydmc, Vzdmc

      real*8 x8, y8, z8,v8x, v8y, v8z

c#ifdef STARFORM
c      nDMparticles=lsp(nspec-1)
c#endif
c#ifndef STARFORM
c      nDMparticles=lsp(nspec)
c#endif
      Rselect=rmax
      Ncount = 0
 300  CONTINUE
      R2 = Rselect**2
      nn=0
      x8 =0.
      y8 = 0.
      z8 = 0.
      v8x=0.
      v8y =0.
      v8z =0.
      wdummy=0.
      do ip= 1, np
         rp2 =(x(ip)-Xdmc)**2 + (y(ip)-Ydmc )**2 + (z(ip) - Zdmc)**2
         if (rp2.lt.R2) then
            nn = nn + 1
            x8 = x8 + x(ip)*pw(ip)
            y8 = y8 + y(ip)*pw(ip)
            z8 = z8 + z(ip)*pw(ip)
            v8x = v8x + vx(ip)*pw(ip)
            v8y = v8y + vy(ip)*pw(ip)
            v8z = v8z + vz(ip)*pw(ip)
            wdummy=wdummy + pw(ip)
         endif
      enddo
      If(nn.eq.0)write (*,*)' Error: no particles found. R=',Rselect
      Xdmc = x8/wdummy ! nn
      Ydmc = y8 /wdummy !  nn
      Zdmc = z8 / wdummy ! nn
      Vxdmc = v8x / wdummy ! nn
      Vydmc = v8y / wdummy ! nn
      Vzdmc = v8z / wdummy ! nn
c      write (*,12) Ncount,  Rselect*1000.0*aexpn*r0*hubble ,
c     &  (xdmc-1.)*1000.0*aexpn*r0*hubble,
c     &  (ydmc-1.)*1000.0*aexpn*r0*hubble,  
c     &   (zdmc-1.)*1000.0*aexpn*r0*hubble,
c     $ vxdmc* v0 /aexpn, vydmc * v0 /aexpn , vzdmc * v0 / aexpn,
c     & nn 
      Rselect = Rselect / 1.1
      Ncount = Ncount + 1
      if ( Ncount.lt.100.AND.Rselect.gt.rmin.AND.nn.gt.20) goto 300
 12   format (I3, 7(f10.3), I6)
      return
      end
c--------------------------------------------------------------------
c--------------------------------------------------------------------
c--------------------------------------------------------------------
      SUBROUTINE STARcenter (rmin, rmax,
     &     Xdmc , Ydmc , Zdmc, Vxdmc, Vydmc, Vzdmc)
c-------------------------------------------------------------------- 
c     Find center of the stellar distribution 
c     based on the iterative method in PMgalaxy.f
c
c  input in code units
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real Xdmc , Ydmc , Zdmc, Vxdmc, Vydmc, Vzdmc
      real*8 x8, y8, z8,v8x, v8y, v8z
      Rselect=rmax
      Ncount = 0
 300  CONTINUE
      R2 = Rselect**2
      nn=0
      x8 =0.
      y8 = 0.
      z8 = 0.
      v8x=0.
      v8y =0.
      v8z =0.
      do ip= nsp(nspec,1),nsp(nspec,2)
         rp2 =(x(ip)-Xdmc)**2 + (y(ip)-Ydmc )**2 + (z(ip) - Zdmc)**2
         if (rp2.lt.R2) then
            nn = nn + 1
            x8 = x8 + x(ip)
            y8 = y8 + y(ip)
            z8 = z8 + z(ip)
            v8x = v8x + vx(ip)
            v8y = v8y + vy(ip)
            v8z = v8z + vz(ip)
         endif
      enddo
      If(nn.eq.0)write (*,*)' Error: no particles found. R=',Rselect
      Xdmc = x8/nn
      Ydmc = y8 / nn
      Zdmc = z8 / nn
      Vxdmc = v8x / nn
      Vydmc = v8y / nn
      Vzdmc = v8z / nn
c      write (*,12) Ncount,  Rselect*1000.0*aexpn*r0*hubble ,
c     &  (xdmc-1.)*1000.0*aexpn*r0*hubble,
c     &  (ydmc-1.)*1000.0*aexpn*r0*hubble,  
c     &   (zdmc-1.)*1000.0*aexpn*r0*hubble,
c     $ vxdmc* v0 , vydmc * v0 , vzdmc * v0,
c     & nn 
      Rselect = Rselect / 1.1
      Ncount = Ncount + 1
      if ( Ncount.lt.100.AND.Rselect.gt.rmin.AND.nn.gt.20) goto 300
 12   format (I3, 7(f10.3), I6)
      return
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE GasCenter (rmin, Xc0 , Yc0 ,Zc0)
c-------------------------------------------------------------------- 
c     Find center of the gas as the maximum density/min potential inside a
c     sphere rmin around a proposed center (Xc0,Yc0,Zc0)
c
c  input in code units
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real Xc0 , Yc0 , Zc0

      rmin2=rmin**2
      rhomax=-1.e20
      phimin= 1.e20
      idcenter=0

      call Get_MaxLevelNow ()
      DO Level =  MaxLevelNow-1, MaxLevelNow 
         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel ) 
         do i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idc = icell + i2
               if ( iOctCh(idc) .eq. nil ) then ! only leaves
                  call  Ps ( idc , Posx,Posy,Posz )
                  r2=(Posx-Xc0)**2+(Posy-Yc0)**2+(Posz-Zc0)**2
                  if (r2.lt.rmin2) then ! cell inside rmax
                     rho00 =  hvar(1,idc)
                     phi =  var(2,idc)
c$$$                     if (rh00.gt.rhomax) then
c$$$                        rhomax =  rho00 
c$$$                        idcenter=idc
                     if (phi.lt.phimin) then
                        phimin=phi
                        idcenter=idc
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
      call  Ps ( idcenter , Posx,Posy,Posz )
      write (*,*)  'cell center id: ', idcenter
      Xc0=Posx
      Yc0=Posy
      Zc0=Posz
      end
c--------------------------------------------------------------------
c-------------------------------------------------------------------- 
      SUBROUTINE Xvir(deltavir,xc,yc,zc,rmin,rmax,Rvir,aMvir,Vcmax)
c-------------------------------------------------------------------- 
c     Compute Rvir, Mvir,Vcmax
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'

      real xc,yc,zc,rmin,rmax
      double precision xx, yy, zz, dM
	real Rvir, aMvir, Vcmax
      rhodmbackground= Om0 ! in units of rhoCr (code units)
!      deltavir =340.
	Vcmax=0.
      dr=(1.*hubble) /(r0*1000.0)

      Rselect = rmin ! initial

 300  dM=0.
      Rselect =  Rselect + dr
      if ( Rselect.gt.rmax) then
         write (*,*)  ' Rvir exceeds rmax'
         stop
      endif
      R2 = Rselect**2

      DO ip=1,np
         xx = x(ip)
         yy = y(ip)
         zz = z(ip)
         rp2 = (xx-xc)**2 + (yy-yc)**2 + (zz-zc)**2
         if ( rp2.lt.R2) then
            dM = dM + pw(ip) 
         endif
      ENDDO
	Vc = dM / Rselect
	if (Vc.gt.Vcmax) Vcmax = Vc 
      overD = (dM/ (4.18879*Rselect**3)) /  rhodmbackground ! in code units
      if ( overD.gt.deltavir) GOTO 300
      Rvir = Rselect*r0*1000./hubble
      aMvir = deltavir*rhodmbackground* 
     &        4.18879 * (Rvir*r0)**3 * rho0 / hubble**3 
	Vcmax =1.15206d8*sqrt( (Vcmax*aM0*hubble) / (r0*aexpn*3.0856d24)) 
	return
        end
c-------------------------------------------------------------------- 
C-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      function fitNFW(deltavir,x,c)
c-------------------------------------------------------------------- 
c  NFW profile fit.
c     rho in units of the critical density (rho0) and
c     radius in units of virial radius.
c-------------------------------------------------------------------- 
      fitNFW = (deltavir/3.) * (c**3/(log(1+c) - c/(1+c))) /
     &  ( x*c*(1+c *x)**2 )

      end
c-------------------------------------------------------------------- 
      SUBROUTINE MinMax(filename)
c
c     check range of all quantities
c
c-------------------------------------------------------------------- 
      character*90 filename
      include 'a_tree.h'
      include 'a_control.h'

      ai=1./aexpn
      ai2 =1./aexpn**2
      ai3=1./aexpn**3
      ai5 =1./aexpn**5
c.... factor converting gas density in code units to n_H in cm^-3
      fact_nH = 1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3
      rho0C =  rho0 * ai3 / 1.E+18     ! in M_sun/pc3 (proper units)
      rhomax= -1.e20
      rhomin =  1.e20

      Emax = -1.e20
      Emin =  1.e20

      Vxmax =  -1.e20
      Vxmin = 1.e20
      Vymax =  -1.e20
      Vymin = 1.e20
      Vzmax =   -1.e20
      Vzmin = 1.e20

      Prmax =   -1.e20
      Prmin =  1.e20

      Gmax =  -1.e20
      Gmin =  1.e20

      Eimax =  -1.e20
      Eimin =  1.e20

      rhoTmax =  -1.e20
      rhoTmin =   1.e20

      phiNmax =  -1.e20
      phiNmin = 1.e20

      phiOmax =  -1.e20
      phiOmin =  1.e20

      Tgmax =  -1.e20
      Tgmin =  1.e20

      Sgmax =  -1.e20
      Sgmin = 1.e20

      call Get_MaxLevelNow ()
      DO i1 = 1 , ncell0
         iCh1 = iOctCh(i1)
         if ( iOctCh(i1) .eq. nil ) then ! only leaves

            rho00 =  hvar(1,i1)
            Et = hvar(2,i1)
            vx0 =  hvar(3,i1) / rho00
            vy0 = hvar(4,i1) / rho00
            vz0 = hvar(5,i1) / rho00
            Pre = hvar(6,i1)
            gammaC = hvar(7,i1)
            Ei =  hvar(8,i1)

            rho00T =  var(1,i1)
            phiN = var(2,i1)
            phio = var(3,i1)

            rhomax = max (rhomax, rho00 )
            rhomin = min (rhomin,  rho00 )
            Emax = max (Emax , Et )
            Emin = min (Emin , Et )
            Vxmax = max (Vxmax , vx0 )
            Vxmin = min ( Vxmin , vx0 )
            Vymax = max (Vymax, vy0 )
            Vymin = min (Vymin , vy0 )
            Vzmax = max (Vzmax, vz0 )
            Vzmin = min (Vzmin, vz0 )
            Prmax = max (Prmax, Pre )
            Prmin = min (Prmin , Pre )
            Gmax = max (Gmax , gammaC )
            Gmin = min (Gmin, gammaC )
            Eimax = max (Eimax, Ei)
            Eimin = min (Eimin, Ei )
            
            rhoTmax = max (rhoTmax, rho00T )
            rhoTmin = min (rhoTmin , rho00T )
            phiNmax = max (phiNmax , phiN )
            phiNmin = min (phiNmin , phiN )
            phiOmax = max (phiOmax, phio )
            phiOmin = min ( phiOmin, phio )

            Tg=(gammaC -1.)*Ei / rho00
            Tgmax =max(Tgmax, Tg)
            Tgmin =min(Tgmin, Tg)

            Sg = Pre / rho00**gammaC
            Sgmax = max ( Sgmax, Sg)
            Sgmin = min ( Sgmin, Sg)
          endif
       ENDDO
       DO Level =  MinLevel+1, MaxLevelNow 
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do i1 = 1 , nLevel
             icell = iSelect(i1)
             do i2 = 0 , 7
                idc = icell + i2
                if ( iOctCh(idc) .eq. nil ) then ! only leaves

                   rho00 =  hvar(1,idc)
                   Et = hvar(2,idc)
                   vx0 =  hvar(3,idc) / rho00
                   vy0 = hvar(4,idc) / rho00
                   vz0 = hvar(5,idc) / rho00 
                   Pre = hvar(6,idc)
                   gammaC = hvar(7,idc)
                   Ei =  hvar(8,idc)

                   rho00T =  var(1,idc)
                   phiN = var(2,idc)
                   phio = var(3,idc)

                   rhomax = max (rhomax, rho00 )
                   rhomin = min (rhomin,  rho00 )
                   Emax = max (Emax , Et )
                   Emin = min (Emin , Et )
                   Vxmax = max (Vxmax , vx0 )
                   Vxmin = min ( Vxmin , vx0 )
                   Vymax = max (Vymax, vy0 )
                   Vymin = min (Vymin , vy0 )
                   Vzmax = max (Vzmax, vz0 )
                   Vzmin = min (Vzmin, vz0 )
                   Prmax = max (Prmax, Pre )
                   Prmin = min (Prmin , Pre )
                   Gmax = max (Gmax , gammaC )
                   Gmin = min (Gmin, gammaC )
                   Eimax = max (Eimax, Ei )
                   Eimin = min (Eimin, Ei )
            
                   rhoTmax = max (rhoTmax, rho00T )
                   rhoTmin = min (rhoTmin , rho00T )
                   phiNmax = max (phiNmax , phiN )
                   phiNmin = min (phiNmin , phiN )
                   phiOmax = max (phiOmax, phio )
                   phiOmin = min ( phiOmin, phio )

                   Tg=(gammaC-1.)*Ei / rho00
                   Tgmax =max(Tgmax, Tg)
                   Tgmin =min(Tgmin, Tg)

                   Sg = Pre / rho00**gammaC
                   Sgmax = max ( Sgmax, Sg)
                   Sgmin = min ( Sgmin, Sg)
                endif
             enddo
          enddo
       ENDDO
! Write:
            nfpn = index ( filename , ' ' ) - 1
      open ( 12 , file = filename(1:nfpn) )
       write (12,*)  'Check units:'
      write (12,'(1x,''r0 [h^{-1} Mpc]    = '',1pe12.4)') r0
      write (12,'(1x,''t0 [yrs]           = '',1pe12.4)') t0
      write (12,'(1x,''rho0 [M_sun/Mpc^3] = '',1pe12.4)') rho0
      write (12,'(1x,''v0 [km/s]          = '',1pe12.4)') v0
      write (12,'(1x,''T0 [K]             = '',1pe12.4)') T_0
      write (12,'(1x,''P0 [g/cm/s^2]      = '',1pe12.4)') P0
      write (12,'(1x,''S0 [keV cm^2]      = '',1pe12.4)') S_0
      write (12,*)  'Check range of parameters:'
      write (12,*)  'in physical and proper units'
      write(12,*) 'rhogmin,max [M_sun/pc^3]= ',rhomin*rho0C , 
     &     rhomax*rho0C 
       write (12,*)  'Emin,max= ', Emin*P0*ai5, Emax*P0*ai5
       write (12,*)  'Vxmin,max= ' ,Vxmin*v0*ai, Vxmax*v0*ai
       write (12,*)  'Vymin,max= ', Vymin*v0*ai, Vymax*v0*ai
       write (12,*)  'Vzmin,max= ', Vzmin*v0*ai, Vzmax*v0*ai
       write (12,*)  'Prmin,max= ', Prmin*P0*ai5, Prmax*P0*ai5
       write (12,*)  'Gmin,max= ', Gmin, Gmax
       write (12,*)  'Eimin,max= ', Eimin*P0*ai5, Eimax*P0*ai5
       write (12,*)  'Tgmin, max= ', Tgmin*T_0*ai2,Tgmax*T_0*ai2
       write (12,*)  'Sgmin, max= ', Sgmin*S_0 , Sgmax*S_0
      write (12,*)  'In code units...'
       write (12,*)  'rhoTmin,max= ', rhoTmin , rhoTmax
       write (12,*)  'phiNmin,max= ', phiNmin , phiNmax
       write (12,*)  'phiOmin,max= ', phiOmin , phiOmax
       write (12,*)  'rhoG (cm-3) = ', rhomin*fact_nH  , rhomax*fact_nH 
 

       END
c-------------------------------------------------------------------- 
      SUBROUTINE MinMaxCU(filename)
c
c     check range of all quantities in code units
c
c-------------------------------------------------------------------- 
      character*90 filename
      include 'a_tree.h'
      include 'a_control.h'

      rhomax= -1.e20
      rhomin =  1.e20

      Emax = -1.e20
      Emin =  1.e20

      Vxmax =  -1.e20
      Vxmin = 1.e20
      Vymax =  -1.e20
      Vymin = 1.e20
      Vzmax =   -1.e20
      Vzmin = 1.e20

      Prmax =   -1.e20
      Prmin =  1.e20

      Gmax =  -1.e20
      Gmin =  1.e20

      Eimax =  -1.e20
      Eimin =  1.e20

      rhoTmax =  -1.e20
      rhoTmin =   1.e20

      phiNmax =  -1.e20
      phiNmin = 1.e20

      phiOmax =  -1.e20
      phiOmin =  1.e20

      Tgmax =  -1.e20
      Tgmin =  1.e20

      Sgmax =  -1.e20
      Sgmin = 1.e20

      ZIImax=-1.e35
      ZIImin=-ZIImax

      ZIamax=-1.e35
      ZIamin=-ZIImax

      call Get_MaxLevelNow ()
      DO i1 = 1 , ncell0
         iCh1 = iOctCh(i1)
         if ( iOctCh(i1) .eq. nil ) then ! only leaves

            rho00 =  hvar(1,i1)
            Et = hvar(2,i1)
            vx0 =  hvar(3,i1) 
            vy0 = hvar(4,i1) 
            vz0 = hvar(5,i1) 
            Pre = hvar(6,i1)
            gammaC = hvar(7,i1)
            Ei =  hvar(8,i1)
            ZII = hvar(9,i1)
            Zia = hvar(10,i1)

            rho00T =  var(1,i1)
            phiN = var(2,i1)
            phio = var(3,i1)

            rhomax = max (rhomax, rho00 )
            rhomin = min (rhomin,  rho00 )
            Emax = max (Emax , Et )
            Emin = min (Emin , Et )
            Vxmax = max (Vxmax , vx0 )
            Vxmin = min ( Vxmin , vx0 )
            Vymax = max (Vymax, vy0 )
            Vymin = min (Vymin , vy0 )
            Vzmax = max (Vzmax, vz0 )
            Vzmin = min (Vzmin, vz0 )
            Prmax = max (Prmax, Pre )
            Prmin = min (Prmin , Pre )
            Gmax = max (Gmax , gammaC )
            Gmin = min (Gmin, gammaC )
            Eimax = max (Eimax, Ei)
            Eimin = min (Eimin, Ei )
            ZIImax = max ( ZIImax, ZII)
            ZIImin = min (ZIImin, ZII)
            ZIamax = max ( ZIamax, ZIa)
            ZIamin = min (ZIamin, ZIa)
            
            rhoTmax = max (rhoTmax, rho00T )
            rhoTmin = min (rhoTmin , rho00T )
            phiNmax = max (phiNmax , phiN )
            phiNmin = min (phiNmin , phiN )
            phiOmax = max (phiOmax, phio )
            phiOmin = min ( phiOmin, phio )

            Tg=(gammaC -1.)*Ei / rho00
            Tgmax =max(Tgmax, Tg)
            Tgmin =min(Tgmin, Tg)

            Sg = Pre / rho00**gammaC
            Sgmax = max ( Sgmax, Sg)
            Sgmin = min ( Sgmin, Sg)
          endif
       ENDDO
       nref = 0
       DO Level =MinLevel+1, MaxLevelNow 
         
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do i1 = 1 , nLevel
             icell = iSelect(i1)
             do i2 = 0 , 7
                idc = icell + i2
                if ( iOctCh(idc) .eq. nil ) then ! only leaves
                   rho00 =  hvar(1,idc)
                   Et = hvar(2,idc)
                   vx0 =  hvar(3,idc) 
                   vy0 = hvar(4,idc) 
                   vz0 = hvar(5,idc)  
                   Pre = hvar(6,idc)
                   gammaC = hvar(7,idc)
                   Ei =  hvar(8,idc)
                   ZII = hvar(9,idc)
                   ZIa = hvar(10,idc)
                   if (ZIa.eq.1.) write (*,*)  ZIa, idc

                   rho00T =  var(1,idc)
                   phiN = var(2,idc)
                   phio = var(3,idc)
! To check refinement conditions:
!                   refi = rho00*3.0517*10.**-5 / (1.49*10.**-2) 
!                   if ( refi.gt.0.2) then
!                      write (*,*)  idc, rho00, refi
!                      nref = nref + 1
!                   endif

                   rhomax = max (rhomax, rho00 )
                   rhomin = min (rhomin,  rho00 )
                   Emax = max (Emax , Et )
                   Emin = min (Emin , Et )
                   Vxmax = max (Vxmax , vx0 )
                   Vxmin = min ( Vxmin , vx0 )
                   Vymax = max (Vymax, vy0 )
                   Vymin = min (Vymin , vy0 )
                   Vzmax = max (Vzmax, vz0 )
                   Vzmin = min (Vzmin, vz0 )
                   Prmax = max (Prmax, Pre )
                   Prmin = min (Prmin , Pre )
                   Gmax = max (Gmax , gammaC )
                   Gmin = min (Gmin, gammaC )
                   Eimax = max (Eimax, Ei )
                   Eimin = min (Eimin, Ei )
                   ZIImax = max (ZIImax, ZII)
                   ZIImin = min (ZIImin, ZII)
                   ZIamax = max (ZIamax, ZIa)
                   ZIamin = min (ZIamin, ZIa)
            
                   rhoTmax = max (rhoTmax, rho00T )
                   rhoTmin = min (rhoTmin , rho00T )
                   phiNmax = max (phiNmax , phiN )
                   phiNmin = min (phiNmin , phiN )
                   phiOmax = max (phiOmax, phio )
                   phiOmin = min ( phiOmin, phio )

                   Tg=(gammaC-1.)*Ei / rho00
                   Tgmax =max(Tgmax, Tg)
                   Tgmin =min(Tgmin, Tg)

                   Sg = Pre / rho00**gammaC
                   Sgmax = max ( Sgmax, Sg)
                   Sgmin = min ( Sgmin, Sg)
                endif
             enddo
          enddo
       ENDDO
! Write:
      nfpn = index ( filename , ' ' ) - 1
      open ( 12 , file = filename(1:nfpn) )

      write (12,*)  'Check range of parameters:'
      write (12,*)  'in code units'
       write(12,*) 'rhogmin,max= ',rhomin, rhomax
       write (12,*)  'Emin,max= ', Emin, Emax
       write (12,*)  'Pxmin,max= ' ,Vxmin, Vxmax
       write (12,*)  'Pymin,max= ', Vymin, Vymax
       write (12,*)  'Pzmin,max= ', Vzmin, Vzmax
       write (12,*)  'Prmin,max= ', Prmin, Prmax
       write (12,*)  'Gmin,max= ', Gmin, Gmax
       write (12,*)  'Eimin,max= ', Eimin, Eimax
       write (12,*)  'Tgmin, max= ', Tgmin,Tgmax
       write (12,*)  'Sgmin, max= ', Sgmin , Sgmax
        write (12,*) 'ZIImin, ZIImax= ', ZIImin, ZIImax      
        write (12,*) 'ZIamin, ZIamax= ', ZIamin, ZIamax  
      write (12,*)  'In code units...'
       write (12,*)  'rhoTmin,max= ', rhoTmin , rhoTmax
       write (12,*)  'phiNmin,max= ', phiNmin , phiNmax
       write (12,*)  'phiOmin,max= ', phiOmin , phiOmax
 

       END
c-------------------------------------------------------------------- 
      SUBROUTINE  slicing()
c-------------------------------------------------------------------- 
c
c     Visualize a slice of the simulation at zero level
c     with k=cte.
c
c---------------------------------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
c
c$$$      parameter (N=ng)
c$$$      real slice(N*N)
c$$$      parameter (kslice= 64)
c$$$      EQUIVALENCE (slice(1) , var(2,(kslice-1)*N*N+1 ) )
c$$$! Plot:
c$$$      sum=0.
c$$$      Potmin=1000.
c$$$      Potmax=-1000.
c$$$      open ( 1 , file = 'slice0.dat' )
c$$$      do j=1,N
c$$$         do i=1,N
c$$$            write(1,*) j , i , slice(i+(j-1)*N)
c$$$            sum = sum +  slice(i+(j-1)*N)
c$$$            Potmin= min (Potmin, slice(i+(j-1)*N))
c$$$            Potmax= max(Potmax,slice(i+(j-1)*N))
c$$$         enddo
c$$$      enddo
c$$$      write (*,*)  sum, Potmin, Potmax
      END
c$$$c-------------------------------------------------------------------- 
c$$$      SUBROUTINE particles_plot()
c$$$c
c$$$c     Plot in sm a selected area of the simulation.
c$$$c
c$$$c     Only for DM particles.
c$$$c
c$$$c-------------------------------------------------------------------- 
c$$$      character Option
c$$$      parameter (NXY=20)
c$$$      real x(NXY),y(NXY)
c$$$      integer sm_device
c$$$c     
c$$$      do 1 i=1,NXY
c$$$         x(i) = i
c$$$         y(i) = i*i
c$$$ 1    continue
c$$$
c$$$      IFlag =sm_device('x11 -bg black -fg white ')
c$$$      call sm_graphics
c$$$
c$$$      call sm_limits(0.,1.,0.,1.)
c$$$      call sm_ctype('red')
c$$$ 15   Do i=1,20
c$$$      call sm_box(1,2,0,0)
c$$$      call sm_ptype(40.,1)
c$$$      call sm_points(x,y,NXY)
c$$$      call sm_xlabel('X axis')
c$$$      call sm_ylabel('Not x axis')
c$$$
c$$$      write (*,*)  'Do it again? y/n'
c$$$      read (*,*) Option
c$$$      if (Option.eq.'y'.or.Option.eq.'Y') then 
c$$$         call sm_erase
c$$$         call sm_ctype('yellow')
c$$$      else
c$$$         stop
c$$$      endif
c$$$      ENDDO
c$$$
c$$$      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE angularM(vcx,vcy,vcz,xc,yc,zc,isys,Rs,aLx,aLy,aLz)
c-------------------------------------------------------------------- 
c
c     Calculate the total angular momentum inside a 
c     sphere Rs for a given center. 
c     isys specifies the system:
c
c     isys=1 --> DM particles
c     isys=2 --> stars
c--------------------------------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      double precision sumx, sumy, sumz
      sumx=0.0
      sumy=0.0
      sumz=0.0
      nn=0
      Rs2 = Rs**2
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      IF (isys.eq.1) then       ! DM particles:
#ifdef STARFORM
      nDMparticles=lsp(nspecies-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspecies)
#endif
         DO ip=1,nDMparticles
            xp = x(ip) -xc
            yp = y(ip) - yc
            zp = z(ip) -zc
            rp2 = xp**2 + yp**2 + zp**2
            vpx =  vx(ip) - vcx
            vpy =  vy(ip) - vcy
            vpz =  vz(ip) - vcz
            IF ( rp2.lt.Rs2) then
               nn= nn+1
               sumx = sumx + pw(ip)*( yp*vpz - zp*vpy)
               sumy = sumy + pw(ip)*( zp*vpx - xp*vpz) 
               sumz = sumz + pw(ip)*( xp*vpy - yp*vpx)
            ENDIF
         ENDDO
         write (*,*)  'The number of particles inside ',Rs,' is ',nn
         aLx = sumx
         aLy = sumy
         aLz = sumz
      ELSEIF (isys.eq.2) then   ! Stars
check         write (*,*)  nsp(nspec,1),nsp(nspec,2)
         DO ip=nsp(nspecies,1),nsp(nspecies,2)
            xp = x(ip) -xc
            yp = y(ip) - yc
            zp = z(ip) -zc
            rp2 = xp**2 + yp**2 + zp**2
            vpx =  vx(ip) - vcx
            vpy =  vy(ip) - vcy
            vpz =  vz(ip) - vcz
            IF ( rp2.lt.Rs2) then
               nn = nn+1
               sumx = sumx + pw(ip)*( yp*vpz - zp*vpy)
               sumy = sumy + pw(ip)*( zp*vpx - xp*vpz) 
               sumz = sumz + pw(ip)*( xp*vpy - yp*vpx)
            ENDIF
         ENDDO
         write (*,*)  'The number of particles inside ',Rs*r0P,' is ',nn
            aLx = sumx
            aLy = sumy
            aLz = sumz
      ENDIF
      END
C-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE NewAxes(aLx,aLy,aLz,v1,v2,v3)
c-------------------------------------------------------------------- 
c
c     Find new system of coordinates 
c     in which aM L=(aLx,aLy,aLz) is parallel to V3
c
c     based on PMgalaxy.f
c
c     input in code units
c-------------------------------------------------------------------- 
      real Lx, Ly, Lz, L , Ld
      real v1(3) , v2(3), v3(3)

! Find new system of coordinates: ex,ey,ez
      L = sqrt( aLx**2 +aLy**2 + aLz**2)
      Lx = aLx/L
      Ly = aLy/L
      Lz = aLz/L
      Ld =sqrt(Lx**2 +Lz**2)
         ex1 =  Lz/Ld
         ey1 = 0.
         ez1 = -Lx/Ld
         ex2 = -Lx*abs(Ly)/Ld
         ey2 =  Ld
         ez2 = -Lz*abs(Ly)/Ld
c         If(Ly.lt.0.) ey2 = -Ld
c         If(Ly.lt.0.) then
c            ey2 = -Ld
c            ex1=-ex1
c            ey1=-ey1
c            ez1=-ez1
c         endif
         IF (Ly.lt.0.) then
            ex2=-ex2
            ez2=-ez2
         ENDIF
            
      ex3 =Lx
      ey3 =Ly
      ez3 =Lz
c                     Check if the vectors are (1) unit (2) orthogonal
      prod1 = ex1*Lx +ey1*Ly +ez1*Lz
      prod2 = ex2*Lx +ey2*Ly +ez2*Lz
      prod3 = ex1*ex2 +ey1*ey2 +ez1*ez2
      a1       = Lx**2 +Ly**2 +Lz**2
      a2       = ex1**2 +ey1**2 +ez1**2
      a3       = ex2**2 +ey2**2 +ez2**2
      write (*,*)  'New system of coordinates: Vector 1, 2 ,3' 
      write (*,*) ' Vector 1=',ex1,ey1,ez1 
      write (*,*) ' Vector 2=',ex2,ey2,ez2 
      write (*,*) ' Vector 3=',ex3,ey3,ez3 
      write (*,*) ' Vector lengths=', a1,a2,a3,' must be 1'
      write (*,*) ' Vector products=', prod1,prod2,prod3,' must be 0'
      if(abs(a1-1.).gt.1.e-2)stop
      if(abs(a2-1.).gt.1.e-2)stop
      if(abs(a2-1.).gt.1.e-2)stop
      if(abs(prod1).gt.1.5e-1)stop
      if(abs(prod2).gt.1.5e-1)stop
      if(abs(prod3).gt.1.5e-1)stop
      v1(1) = ex1
      v1(2) = ey1
      v1(3) = ez1
      v2(1) = ex2
      v2(2) = ey2
      v2(3) = ez2
      v3(1) = ex3
      v3(2) = ey3
      v3(3) = ez3
      return
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE CylinPS(Vcx,Vcy,Vcz,xc,yc,zc,
     &        ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ,tage)
c-------------------------------------------------------------------- 
c     Profiles in cylindrical coordinates,
c     where v1,v2,v3 define its axis
c     and Vcx,Vcy,Vcz,xc,yc,zc, the velocity and position of its center
c     ZDmax defines the maximum height of the sampled disk.
c     in code unis
c-------------------------------------------------------------------- 
c
c     Velocity and density profiles of the stellar component
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
      
       real vcx, vcy, vcz,xc,yc,zc,rmin, rmax,tage
      integer nrbin
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles

      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),vbin(0:nmaxbins) 
      real vrbin(0:nmaxbins), vphbin(0:nmaxbins)
      integer nbin(0:nmaxbins)
      real vxbin(0:nmaxbins) , vybin(0:nmaxbins), vzbin(0:nmaxbins)
      real rmsVxbin(0:nmaxbins), vx2Sum(0:nmaxbins)
      real rmsVybin(0:nmaxbins), vy2Sum(0:nmaxbins)
      real rmsVzbin(0:nmaxbins), vz2Sum(0:nmaxbins)
      real rmsVbin(0:nmaxbins), v2Sum(0:nmaxbins)
      real rmsVrbin(0:nmaxbins), vr2Sum(0:nmaxbins)
      real rmsVphbin(0:nmaxbins), vph2sum(0:nmaxbins)

      real dmbin(0:nmaxbins)
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      Sigma0P =  rho0 * r0 / a2 / 1.E+12     ! in M_sun/pc2 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      v0P = v0 / aexpn ! in km/s (proper units)
      tfact = t0*aexpn**2
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
!     Initialize profiles:
      do i = 0 , nmaxbins
         rbin(i) = 0.
         nbin(i) = 0
         vbin(i)=0.
         vrbin(i)=0.
         vphbin(i)=0.
         
         vxbin(i)=0.
         vybin(i)=0.
         vzbin(i)=0.

         vx2sum(i)=0.
         rmsVxbin(i)=0.
         vy2sum(i)=0.
         rmsVybin(i)=0.
         vz2sum(i)=0.
         rmsVzbin(i)=0.

         v2sum(i) = 0.
         rmsVbin(i) =0.
         vr2sum(i) = 0.
         rmsVrbin(i) =0.
         vph2sum(i) =0.
         rmsVphbin(i) =0.

         dmbin(I)=0.
      enddo
      Nstars=0
      DO ip=nsp(nspecies,1),nsp(nspecies,2) ! STARS
         istar =ip - nsp(nspecies,1) + 1
         if ( (t-tbirth(istar))*tfact.lt.tage) then
            Nstars=Nstars+1
            xp = x(ip) -xc
            yp = y(ip) - yc
            zp = z(ip) -zc
            rp3D2 = xp**2 + yp**2 + zp**2 
            Xpp =xp*AxX(1) +yp*AxX(2) +zp*AxX(3)
            Ypp =xp*AxY(1) +yp*AxY(2) +zp*AxY(3)
            Zpp =xp*AxZ(1) +yp*AxZ(2) +zp*AxZ(3)
            rp3D2p = xpp**2 + ypp**2 + zpp**2 
         if (abs( rp3D2 -  rp3D2p).gt.0.01) then
            write (*,*) 'error in the rotation of coordinates. Rp3D2= ',
     &           rp3D2, rp3D2p  
            write (*,*)  'it will not stop'
!            STOP
         endif
         IF ( ABS(ZPP).LT.ZDmax) then
            Xp = Xpp
            Yp = Ypp
            Zp = Zpp
            vpx = vx(ip) - vcx
            vpy = vy(ip) - vcy
            vpz = vz(ip) - vcz
            Vp3D2 = vpx**2 + vpy**2 + vpz**2
            Vppx =vpx*AxX(1) + vpy*AxX(2) + vpz*AxX(3)
            Vppy =vpx*AxY(1) + vpy*AxY(2) + vpz*AxY(3)
            Vppz =vpx*AxZ(1) + vpy*AxZ(2) + vpz*AxZ(3)         
            vp3D2p = vppx**2 + vppy**2 + vppz**2
         if (abs( vp3D2 -  vp3D2p).gt.0.01) then
            write (*,*) 'error in the rotation of coordinates. Vp3D2= ',
     &              vp3D2, vp3D2p  
            write (*,*)  'it will not stop'
!            STOP
         endif
         Vpx=Vppx
         Vpx2 =  Vpx**2
         Vpy=Vppy
         Vpy2 =  Vpy**2
         Vpz=Vppz 
         Vpz2 = Vpz**2
         rp =  sqrt(xp**2+yp**2 )
         vp2 = vpx2 + vpy2 + vpz2
         vp = sqrt ( vp2 )
         if (rp.eq.0.) then
            vpr=0.
            vpr2=0.
            vpph=0.
            vpph2 =0.
         else
            vpr =(Vpx *xp + vpy*yp ) / rp
            vpr2 = vpr**2
            vpph = ( -Vpx *Yp + Vpy*Xp) /rp
            vpph2 = vpph**2
         endif
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            rbin(ibin) = rbin(ibin) + rp
            nbin(ibin) = nbin(ibin) + 1
            vbin(ibin)= vbin(ibin) + vp
            vrbin(ibin)= vrbin(ibin) + vpr
            vphbin(ibin) = vphbin(ibin) + vpph
            vxbin(ibin) = vxbin(ibin) + vpx
            vybin(ibin) = vybin(ibin) + vpy
            vzbin(ibin) = vzbin(ibin) + vpz
            dmbin(ibin) = dmbin(ibin) + pw(ip)
            vx2Sum(ibin) =  vx2Sum(ibin) + Vpx2
            vy2sum(ibin) = vy2sum(ibin) + Vpy2
            vz2sum(ibin) =   vz2sum(ibin) + Vpz2
            v2sum(ibin) = v2sum(ibin) + Vp2
            vr2sum(ibin) =  vr2sum(ibin) + vpr2
            vph2sum(ibin) =  vph2sum(ibin) + vpph2
         endif
      ENDIF
      endif
      ENDDO
         
!     rms velocities without normalization:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rmsVxbin(i) = vx2Sum(i) - vxbin(i)**2 / nbin(i) 
            rmsVybin(i) = vy2Sum(i) - vybin(i)**2 / nbin(i)
            rmsVzbin(i) = vz2Sum(i) - vzbin(i)**2 / nbin(i)
            rmsVbin(i)  = v2Sum(i)   -  vbin(i)**2 / nbin(i)
            rmsVrbin(i) = vr2Sum(i) - vrbin(i)**2 / nbin(i)
            rmsVphbin(i) = vph2sum(i) -vphbin(i)**2  / nbin(i)
         endif
      enddo
! Normalize profiles:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rbin(i) = rbin(i)/  nbin(i) 
            vbin(i)= vbin(i) / nbin(i)
            vrbin(i)= vrbin(i) / nbin(i)
            vphbin(i) = vphbin(i) / nbin(i)
            
            vxbin(i) = vxbin(i) / nbin(i)
            vybin(i) = vybin(i) / nbin(i)
            vzbin(i) = vzbin(i) / nbin(i)

            if (rmsVxbin(i).lt.0.) then
               if (  nbin(i).eq.1) then
                  rmsVxbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVybin(i).lt.0.) then
               if  (  nbin(i).eq.1) then
                  rmsVybin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVzbin(i).lt.0.) then
               if  (  nbin(i).eq.1) then
                  rmsVzbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVbin(i).lt.0.)  then
               if  (  nbin(i).eq.1) then
                  rmsVbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVrbin(i).lt.0.)then
               if  (  nbin(i).eq.1) then
                  rmsVrbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVphbin(i).lt.0.)then
               if  (  nbin(i).eq.1) then
                  rmsVphbin(i)=0.
               else
                  STOP
               endif
            endif
            rmsVxbin(i) = sqrt( rmsVxbin(i) / nbin(i) )
            rmsVybin(i) = sqrt( rmsVybin(i) / nbin(i) )
            rmsVzbin(i) = sqrt( rmsVzbin(i) / nbin(i) )
            rmsVbin(i)   = sqrt( rmsVbin(i) / nbin(i) )
            rmsVrbin(i) = sqrt( rmsVrbin(i) / nbin(i) )
            rmsVphbin(i)= sqrt(  rmsVphbin(i) / nbin(i))
         endif
      enddo
! Density profile:
      DO i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlogmin + float(i-1)*drl)
         endif
         rr = 10.**(rlogmin + float(i)*drl)
         surf = pi * (rr**2 - rl**2)          
         if ( surf .gt. 0. ) then
            dmbin(i) = dmbin(i)/ surf
         else
            write(*,*) 'bad surf =',surf,i, rl, rr
         endif
      ENDDO
!     Write output:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0,(yc-1.)*r0,(zc-1.)*r0,vcx*v0,vcy*v0,vcz*v0
      write(12,11) AxZ(1) , AxZ(2) , AxZ(3)
      write(12,12) rmin*r0P, rmax*r0P , nrbin , ZDmax*r0P
      write (12,13) aexpn , 1/aexpn -1.
      write(12,14) tage
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[Mpc h-1] = ', 3(g13.7,1x),
     &     '# (Vcx,Vcy,Vcz)[Km/s] = ', 3(g13.7,1x))
 11   format ('axis direction= ', 3(g13.7,1x))
 12   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5,' ZDmax= ',g12.6)
 13   format ('Distances in proper Kpc, proper velocities in Km/s and densities in Msun/pc2
     &     for a= ',f6.3, ' and redshift= ', f6.2)
 14   format ('# CYLINDRICAL PROFILES FOR STARS younger than', 
     & E8.1, ' yrs : Columns:')
 15   format ('# <r>',14x,'r_l',14x,' r_m',14x,'r_r',14x,
     &     'v',7x,'vr',7x,'vphi',3x,'rhoStars',11x,
     &     'vx',7x,' vy',7x,'vz',7x,
     &     'rmsVx',5x,'rmsVy ',6x,'rmsVz',4x,
     &     'rmsV',6x,'rmsVr',6x, 'rmsVph',2x,'NPibin')
 16   format('----------------------')
      Npibin=0.
      do i = 0 , nbins
         Npibin =Npibin + nbin(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(12,17)
     &           rbin(i)*r0P, rl*r0P, 
     &           rmid*r0P, rr*r0P, 
     &       vbin(i)*v0P, vrbin(i)*v0P ,vphbin(i)*v0P, dmbin(i)*Sigma0P,
     &           vxbin(i)*v0P , vybin(i)*v0P , vzbin(i)*v0P,
     &           rmsVxbin(i)*v0P , rmsVybin(i)*v0P , rmsVzbin(i)*v0P,
     &           rmsVbin(i)*v0P, rmsVrbin(i)*v0P, rmsVphbin(i)*v0P,
     &           npibin
         endif
      enddo
 17   format (4(2x,g14.7), 3(2x,f8.3), 1x,g14.7, 9(2x,f8.3),  1x, I6)
      close ( 12 )
      write(*,*)'done writing cylindrical profiles for stars. exiting..'
      return 
      END
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE CylinPG(Vcx,Vcy,Vcz,xc,yc,zc,
     &     zDmax,T_gas,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ)
c-------------------------------------------------------------------- 
c     Profiles in cylindrical coordinates,
c     where v1,v2,v3 define its axis
c     and Vcx,Vcy,Vcz,xc,yc,zc, the velocity and position of its center.
c     ZDmax defines the maximum height of the sampled disk.
C     T_gas defines the maximum tempetature of the sampled gas.
c     In code unis.
c-------------------------------------------------------------------- 
c
c     Velocity and density profile of the gas component
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
      
       real vcx, vcy, vcz,xc,yc,zc,rmin, rmax, zDmax
      integer nrbin
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
      double precision xx,yy,zz
      integer m
      parameter ( npoints = 10000000 )
      parameter ( nmaxbins = 1000 )
      real*8 rbin(0:nmaxbins),vbin(0:nmaxbins) 
      real*8 vrbin(0:nmaxbins), vphbin(0:nmaxbins)
      integer nbin(0:nmaxbins)
      real*8 vxbin(0:nmaxbins) , vybin(0:nmaxbins), vzbin(0:nmaxbins)
      real*8 rmsVxbin(0:nmaxbins), vx2Sum(0:nmaxbins)
      real*8 rmsVybin(0:nmaxbins), vy2Sum(0:nmaxbins)
      real*8 rmsVzbin(0:nmaxbins), vz2Sum(0:nmaxbins)
      real*8 rmsVbin(0:nmaxbins), v2Sum(0:nmaxbins)
      real*8 rmsVrbin(0:nmaxbins), vr2Sum(0:nmaxbins)
      real*8 rmsVphbin(0:nmaxbins), vph2sum(0:nmaxbins)
      
      real*8 dmbin(0:nmaxbins)
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      Sigma0P =  rho0 * r0 / a2 / 1.E+12 ! in M_sun/pc2 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      v0P = v0 / aexpn ! in km/s (proper units)
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
!     Initialize profiles:
      do i = 0 , nmaxbins
         rbin(i) = 0.
         nbin(i) = 0
         vbin(i)=0.
         vrbin(i)=0.
         vphbin(i)=0.
         
         vxbin(i)=0.
         vybin(i)=0.
         vzbin(i)=0.

         vx2sum(i)=0.
         rmsVxbin(i)=0.
         vy2sum(i)=0.
         rmsVybin(i)=0.
         vz2sum(i)=0.
         rmsVzbin(i)=0.

         v2sum(i) = 0.
         rmsVbin(i) =0.
         vr2sum(i) = 0.
         rmsVrbin(i) =0.
         vph2sum(i) =0.
         rmsVphbin(i) =0.

         dmbin(I)=0.
      enddo
! Sample of random points in a disk of size 2rmax and
! width 2Zmax
      DO ip=1, npoints
         rr =  rmax * RANDd ( m )
         phi =  pi2 * RANDd ( m )
         zp =zDmax*(2.* RANDd ( m ) - 1. )
         xp = rr * cos(phi)
         yp = rr * sin(phi)
         xx = Xc + Xp*AxX(1) + Yp*AxY(1) + Zp*AxZ(1)
         yy = Yc + Xp*AxX(2) + Yp*AxY(2) + Zp*AxZ(2)
         zz = Zc + Xp*AxX(3) + Yp*AxY(3) + Zp*AxZ(3)
!     Check periodic boundaries is not necesary.
         idcell = iFindCell ( MaxLevel , xx , yy , zz )
         pd= hvar(1,idcell)
         if ( pd .gt. 0. ) then
            T_code= gamma1 * hvar(8,idcell)/hvar(1,idcell)
            ptemp = T_0 * T_code * a2i
         else
            write(*,*)'bad gas density in cell',idcell,pd
         endif           
         IF (ptemp.lt.T_gas) then
   
            vpx= hvar(3,idcell)/pd - vcx
            vpy= hvar(4,idcell)/pd - vcy
            vpz= hvar(5,idcell)/pd - vcz

            Vp3D2 = vpx**2 + vpy**2 + vpz**2
            Vppx =vpx*AxX(1) + vpy*AxX(2) + vpz*AxX(3)
            Vppy =vpx*AxY(1) + vpy*AxY(2) + vpz*AxY(3)
            Vppz =vpx*AxZ(1) + vpy*AxZ(2) + vpz*AxZ(3)         
            vp3D2p = vppx**2 + vppy**2 + vppz**2
            if (abs( vp3D2 -  vp3D2p).gt.0.1) then
             write (*,*)'error in the rotation of coordinates. Vp3D2= ',
     &              vp3D2, vp3D2p  
               stop
            endif
            Vpx=Vppx
            Vpy=Vppy
            Vpz=Vppz
            Vpx2 = Vppx**2 
            Vpy2 = Vppy**2
            Vpz2 = Vppz**2
            rp=rr
            vp2 =  Vpx2 + Vpy2 + Vpz2 
            vp = sqrt ( vp2 )
            if (rp.eq.0.) then
               vpr = 0.
               vpr2 = 0.
               vpph=0.
               vpph2 =0.
            else
               vpr =(Vpx *xp + vpy*yp ) / rp 
               vpr2 = vpr**2
               vpph = ( -Vpx *Yp + Vpy*Xp) /rp
               vpph2 =  vpph **2
            endif
            ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
            if ( ibin .le. nmaxbins ) then  ! Density weighted profiles.
               rbin(ibin) = rbin(ibin) + rp
               nbin(ibin) = nbin(ibin) + 1
               vbin(ibin)= vbin(ibin) + vp*pd
               vrbin(ibin)= vrbin(ibin) + vpr*pd
               vphbin(ibin) = vphbin(ibin) + vpph*pd
               vxbin(ibin) = vxbin(ibin) + vpx*pd
               vybin(ibin) = vybin(ibin) + vpy*pd
               vzbin(ibin) = vzbin(ibin) + vpz*pd
               dmbin(ibin) = dmbin(ibin) + pd
               vx2Sum(ibin) =  vx2Sum(ibin) + Vpx2*pd
               vy2sum(ibin) = vy2sum(ibin) + Vpy2*pd
               vz2sum(ibin) =   vz2sum(ibin) + Vpz2*pd
               v2sum(ibin) = v2sum(ibin) + Vp2*pd
               vr2sum(ibin) =  vr2sum(ibin) + vpr2*pd
               vph2sum(ibin) =  vph2sum(ibin) + vpph2*pd
            endif
         ENDIF
      ENDDO
! rms velocities without normalization:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rmsVxbin(i) = vx2Sum(i) - vxbin(i)**2 /  dmbin(i)
            rmsVybin(i) = vy2Sum(i) - vybin(i)**2 /  dmbin(i)
            rmsVzbin(i) = vz2Sum(i) - vzbin(i)**2 /  dmbin(i)
            rmsVbin(i)  = v2Sum(i)   -  vbin(i)**2 /  dmbin(i)
            rmsVrbin(i) = vr2Sum(i) - vrbin(i)**2 /  dmbin(i)
            rmsVphbin(i) = vph2sum(i) -vphbin(i)**2  /  dmbin(i)
check:
c            write (*,*)  rmsVbin(i) ,  v2Sum(i), vbin(i), dmbin(i)
         endif
      enddo
c$$$! Density profile:
c$$$      DO i = 0 , nbins
c$$$         if ( i .eq. 0 ) then 
c$$$            rl = 0.0
c$$$         else
c$$$            rl = 10.**(rlogmin + float(i-1)*drl)
c$$$         endif
c$$$         rr = 10.**(rlogmin + float(i)*drl)
c$$$         surf = pi * (rr**2 - rl**2)          
c$$$         if ( surf .gt. 0. ) then
c$$$            dmbin(i) = dmbin(i)/ surf
c$$$         else
c$$$            write(*,*) 'bad surf =',surf,i, rl, rr
c$$$         endif
c$$$      ENDDO
! Normalize profiles (Density weighted):
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rbin(i) = rbin(i)/  nbin(i) 
            vbin(i)= vbin(i) /  dmbin(i)
            vrbin(i)= vrbin(i) / dmbin(i)
            vphbin(i) = vphbin(i) / dmbin(i)            
            vxbin(i) = vxbin(i) / dmbin(i)
            vybin(i) = vybin(i) / dmbin(i)
            vzbin(i) = vzbin(i) / dmbin(i)

            if (rmsVxbin(i).lt.0.)  then
               if (abs(rmsVxbin(i)).lt.100.)then
                  rmsVxbin(i)=0.0
               else
                  write (*,*)  'rmsVxbin(i)'
                  write (*,*)  rmsVxbin(i)
                  STOP
               endif
            endif
            if (rmsVybin(i).lt.0.)  then
               if (abs(rmsVybin(i)).lt.100.) then
                  rmsVybin(i)=0.0
               else
                  write (*,*)  'rmsVybin(i)'
                  write (*,*)  rmsVybin(i)
                  STOP
               endif
            endif
            if (rmsVzbin(i).lt.0.)  then
               if (abs(rmsVzbin(i)).lt.100.) then 
                  rmsVzbin(i)=0.0
               else
                  write (*,*)  'rmsVzbin(i)'
                  write (*,*)  rmsVzbin(i)
                  STOP
               endif
            endif
            if (rmsVbin(i).lt.0.)  then
               if (abs(rmsVbin(i)).lt.100.) then
                  rmsVbin(i)=0.0
               else
                  write (*,*)  'rmsVbin(i)'
                  write (*,*)  rmsVbin(i)
                  STOP
               endif
            endif
            if (rmsVrbin(i).lt.0.)  then
               if (abs(rmsVrbin(i)).lt.100.) then
                  rmsVrbin(i)=0.0
               else
                  write (*,*) ' rmsVrbin(i)'
                  write (*,*)  rmsVrbin(i)
                  STOP
               endif
            endif
            if (rmsVphbin(i).lt.0.)  then
               if (abs(rmsVphbin(i)).lt.100.) then
                  rmsVphbin(i)=0.0
               else
                  write (*,*)  'rmsVphbin(i)'
                  write (*,*)  rmsVphbin(i)
                  STOP
               endif
            endif

            rmsVxbin(i) = sqrt( rmsVxbin(i) / dmbin(i) )
            rmsVybin(i) = sqrt( rmsVybin(i) / dmbin(i) )
            rmsVzbin(i) = sqrt( rmsVzbin(i) / dmbin(i) )
            rmsVbin(i)   = sqrt( rmsVbin(i) / dmbin(i) )
            rmsVrbin(i) = sqrt( rmsVrbin(i) / dmbin(i) )
            rmsVphbin(i)= sqrt(  rmsVphbin(i) / dmbin(i) )

            dmbin(i) = dmbin(i) / nbin(i)
         endif
      enddo
!     Write output:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0,(yc-1.)*r0,(zc-1.)*r0,vcx*v0,vcy*v0,vcz*v0
      write(12,11) AxZ(1) , AxZ(2) , AxZ(3)
      write(12,12) rmin*r0P, rmax*r0P , nrbin , ZDmax*r0P
      write (12,13) aexpn , 1/aexpn -1.
      write(12,14) T_gas
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[Mpc h-1] = ', 3(g13.7,1x),
     &     '# (Vcx,Vcy,Vcz)[Km/s] = ', 3(g13.7,1x))
 11   format ('axis direction= ', 3(g13.7,1x))
 12   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5,' ZDmax= ',g12.6)
 13   format ('Distances in proper Kpc, proper velocities in Km/s and densities in Msun/pc2
     &     for a= ',f6.3, ' and redshift= ', f6.2)
 14   format ('# CYLINDRICAL PROFILES FOR GAS with temperature bellow '
     &, E7.1,' :  Columns:  Density weighted.')
 15   format ('# <r>',14x,'r_l',14x,' r_m',14x,'r_r',14x,
     &     'v',7x,'vr',7x,'vphi',11x,'rhoGas',7x,
     &     'vx',7x,' vy',7x,'vz',7x,
     &     'rmsVx',4x,'rmsVy',6x,'rmsVz',6x,
     &     'rmsV',6x,'rmsVr',4x, 'rmsVph',2x,'NVPibin')
 16   format('----------------------')
      Npibin=0.
      do i = 0 , nbins
         Npibin =Npibin + nbin(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(12,17)
     &           rbin(i)*r0P, rl*r0P, 
     &           rmid*r0P, rr*r0P, 
     &       vbin(i)*v0P, vrbin(i)*v0P, vphbin(i)*v0P, dmbin(i)*Sigma0P,
     &           vxbin(i)*v0P , vybin(i)*v0P , vzbin(i)*v0P,
     &           rmsVxbin(i)*v0P , rmsVybin(i)*v0P , rmsVzbin(i)*v0P,
     &           rmsVbin(i)*v0P, rmsVrbin(i)*v0P, rmsVphbin(i)*v0P,
     &           npibin
         endif
      enddo
 17   format (4(2x,g14.7), 3(2x,f8.3), 1x,g14.3, 9(2x,f8.3),  1x, I7)
      close ( 12 )
      write(*,*) 'done writing cylindrical profiles for gas. exiting...'
      return
      END
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE CylinPD(Vcx,Vcy,Vcz,xc,yc,zc,
     &        ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ)
c-------------------------------------------------------------------- 
c     Profiles in cylindrical coordinates,
c     where v1,v2,v3 define its axis
c     and Vcx,Vcy,Vcz,xc,yc,zc, the velocity and position of its center
c     ZDmax defines the maximum height of the sampled disk.
c     in code unis
c-------------------------------------------------------------------- 
c
c     Velocity and density profiles of the DM component
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
      
       real vcx, vcy, vcz,xc,yc,zc,rmin, rmax
      integer nrbin
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles

      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),vbin(0:nmaxbins) 
      real vrbin(0:nmaxbins), vphbin(0:nmaxbins)
      integer nbin(0:nmaxbins)
      real vxbin(0:nmaxbins) , vybin(0:nmaxbins), vzbin(0:nmaxbins)
      real rmsVxbin(0:nmaxbins), vx2Sum(0:nmaxbins)
      real rmsVybin(0:nmaxbins), vy2Sum(0:nmaxbins)
      real rmsVzbin(0:nmaxbins), vz2Sum(0:nmaxbins)
      real rmsVbin(0:nmaxbins), v2Sum(0:nmaxbins)
      real rmsVrbin(0:nmaxbins), vr2Sum(0:nmaxbins)
      real rmsVphbin(0:nmaxbins), vph2sum(0:nmaxbins)

      real dmbin(0:nmaxbins)
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      Sigma0P =  rho0 * r0 / a2 / 1.E+12     ! in M_sun/pc2 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      v0P = v0 / aexpn ! in km/s (proper units)
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
#ifdef STARFORM
      nDMparticles=lsp(nspec-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspec)
#endif
!     Initialize profiles:
      do i = 0 , nmaxbins
         rbin(i) = 0.
         nbin(i) = 0
         vbin(i)=0.
         vrbin(i)=0.
         vphbin(i)=0.
         
         vxbin(i)=0.
         vybin(i)=0.
         vzbin(i)=0.

         vx2sum(i)=0.
         rmsVxbin(i)=0.
         vy2sum(i)=0.
         rmsVybin(i)=0.
         vz2sum(i)=0.
         rmsVzbin(i)=0.

         v2sum(i) = 0.
         rmsVbin(i) =0.
         vr2sum(i) = 0.
         rmsVrbin(i) =0.
         vph2sum(i) =0.
         rmsVphbin(i) =0.

         dmbin(I)=0.
      enddo
      DO ip=1,nDMparticles  ! DM
         xp = x(ip) -xc
         yp = y(ip) - yc
         zp = z(ip) -zc
         rp3D2 = xp**2 + yp**2 + zp**2 
         Xpp =xp*AxX(1) +yp*AxX(2) +zp*AxX(3)
         Ypp =xp*AxY(1) +yp*AxY(2) +zp*AxY(3)
         Zpp =xp*AxZ(1) +yp*AxZ(2) +zp*AxZ(3)
         rp3D2p = xpp**2 + ypp**2 + zpp**2 
         if (abs( rp3D2 -  rp3D2p).gt.0.01) then
            write (*,*) 'error in the rotation of coordinates. Rp3D2= ',
     &           rp3D2, rp3D2p  
            stop
         else
            Xp = Xpp
            Yp = Ypp
            Zp = Zpp
         endif
            vpx = vx(ip) - vcx
            vpy = vy(ip) - vcy
            vpz = vz(ip) - vcz
            Vp3D2 = vpx**2 + vpy**2 + vpz**2
            Vppx =vpx*AxX(1) + vpy*AxX(2) + vpz*AxX(3)
            Vppy =vpx*AxY(1) + vpy*AxY(2) + vpz*AxY(3)
            Vppz =vpx*AxZ(1) + vpy*AxZ(2) + vpz*AxZ(3)         
            vp3D2p = vppx**2 + vppy**2 + vppz**2
            if (abs( vp3D2 -  vp3D2p).gt.0.01) then
            write (*,*) 'error in the rotation of coordinates. Vp3D2= ',
     &           vp3D2, vp3D2p  
            stop
         else
            Vpx=Vppx
            Vpx2 =  Vpx**2
            Vpy=Vppy
            Vpy2 =  Vpy**2
            Vpz=Vppz 
            Vpz2 = Vpz**2
         endif
         IF ( ABS(ZP).LT.ZDmax) then
            rp =  sqrt(xp**2+yp**2 )
            vp2 = vpx**2 + vpy**2 + vpz**2
            vp = sqrt ( vp2 )
            if (rp.eq.0.) then
               vpr=0.
               vpr2=0.
               vpph=0.
               vpph2 =0.
            else
               vpr =(Vpx *xp + vpy*yp ) / rp
               vpr2 = vpr**2
               vpph = ( -Vpx *Yp + Vpy*Xp) /rp
               vpph2 = vpph**2
            endif
            ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
            if ( ibin .le. nmaxbins ) then
               rbin(ibin) = rbin(ibin) + rp
               nbin(ibin) = nbin(ibin) + 1
               vbin(ibin)= vbin(ibin) + vp
               vrbin(ibin)= vrbin(ibin) + vpr
               vphbin(ibin) = vphbin(ibin) + vpph
               vxbin(ibin) = vxbin(ibin) + vpx
               vybin(ibin) = vybin(ibin) + vpy
               vzbin(ibin) = vzbin(ibin) + vpz
               dmbin(ibin) = dmbin(ibin) + pw(ip)
               vx2Sum(ibin) =  vx2Sum(ibin) + Vpx2
               vy2sum(ibin) = vy2sum(ibin) + Vpy2
               vz2sum(ibin) =   vz2sum(ibin) + Vpz2
               v2sum(ibin) = v2sum(ibin) + Vp2
               vr2sum(ibin) =  vr2sum(ibin) + vpr2
               vph2sum(ibin) =  vph2sum(ibin) + vpph2
            endif
         ENDIF
      ENDDO
         
!     rms velocities without normalization:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rmsVxbin(i) = vx2Sum(i) - vxbin(i)**2 / nbin(i) 
            rmsVybin(i) = vy2Sum(i) - vybin(i)**2 / nbin(i)
            rmsVzbin(i) = vz2Sum(i) - vzbin(i)**2 / nbin(i)
            rmsVbin(i)  = v2Sum(i)   -  vbin(i)**2 / nbin(i)
            rmsVrbin(i) = vr2Sum(i) - vrbin(i)**2 / nbin(i)
            rmsVphbin(i) = vph2sum(i) -vphbin(i)**2  / nbin(i)
         endif
      enddo
! Normalize profiles:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rbin(i) = rbin(i)/  nbin(i) 
            vbin(i)= vbin(i) / nbin(i)
            vrbin(i)= vrbin(i) / nbin(i)
            vphbin(i) = vphbin(i) / nbin(i)
            
            vxbin(i) = vxbin(i) / nbin(i)
            vybin(i) = vybin(i) / nbin(i)
            vzbin(i) = vzbin(i) / nbin(i)

            if (rmsVxbin(i).lt.0.) then
               if (  nbin(i).eq.1) then
                  rmsVxbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVybin(i).lt.0.) then
               if  (  nbin(i).eq.1) then
                  rmsVybin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVzbin(i).lt.0.) then
               if  (  nbin(i).eq.1) then
                  rmsVzbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVbin(i).lt.0.)  then
               if  (  nbin(i).eq.1) then
                  rmsVbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVrbin(i).lt.0.)then
               if  (  nbin(i).eq.1) then
                  rmsVrbin(i)=0.
               else
                  STOP
               endif
            endif
            if (rmsVphbin(i).lt.0.)then
               if  (  nbin(i).eq.1) then
                  rmsVphbin(i)=0.
               else
                  STOP
               endif
            endif

            rmsVxbin(i) = sqrt( rmsVxbin(i) / nbin(i) )
            rmsVybin(i) = sqrt( rmsVybin(i) / nbin(i) )
            rmsVzbin(i) = sqrt( rmsVzbin(i) / nbin(i) )
            rmsVbin(i)   = sqrt( rmsVbin(i) / nbin(i) )
            rmsVrbin(i) = sqrt( rmsVrbin(i) / nbin(i) )
            rmsVphbin(i)= sqrt(  rmsVphbin(i) / nbin(i))
         endif
      enddo
! Density profile:
      DO i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlogmin + float(i-1)*drl)
         endif
         rr = 10.**(rlogmin + float(i)*drl)
         surf = pi * (rr**2 - rl**2)          
         if ( surf .gt. 0. ) then
            dmbin(i) = dmbin(i)/ surf
         else
            write(*,*) 'bad surf =',surf,i, rl, rr
         endif
      ENDDO
!     Write output:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0,(yc-1.)*r0,(zc-1.)*r0,vcx*v0,vcy*v0,vcz*v0
      write(12,11) AxZ(1) , AxZ(2) , AxZ(3)
      write(12,12) rmin*r0P, rmax*r0P , nrbin , ZDmax*r0P
      write (12,13) aexpn , 1/aexpn -1.
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[Mpc h-1] = ', 3(g13.7,1x),
     &     '# (Vcx,Vcy,Vcz)[Km/s] = ', 3(g13.7,1x))
 11   format ('axis direction= ', 3(g13.7,1x))
 12   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5,' ZDmax= ',g12.6)
 13   format ('Distances in proper Kpc, proper velocities in Km/s and densities in Msun/pc2
     &     for a= ',f6.3, ' and redshift= ', f6.2)
 14   format ('# CYLINDRICAL PROFILES FOR DM: Columns:')
 15   format ('# <r>',14x,'r_l',14x,' r_m',14x,'r_r',14x,
     &     'v',7x,'vr',7x,'vphi',3x,'rhoDM',11x,
     &     'vx',7x,' vy',7x,'vz',7x,
     &     'rmsVx',5x,'rmsVy ',6x,'rmsVz',
     &     'rmsV',6x,'rmsVr',6x, 'rmsVph',2x,'NPibin')
 16   format('----------------------')
      Npibin=0.
      do i = 0 , nbins
         Npibin =Npibin + nbin(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(12,17)
     &           rbin(i)*r0P, rl*r0P, 
     &           rmid*r0P, rr*r0P, 
     &       vbin(i)*v0P, vrbin(i)*v0P , vphbin(i)*v0P,dmbin(i)*Sigma0P,
     &           vxbin(i)*v0P , vybin(i)*v0P , vzbin(i)*v0P,
     &           rmsVxbin(i)*v0P , rmsVybin(i)*v0P , rmsVzbin(i)*v0P,
     &           rmsVbin(i)*v0P, rmsVrbin(i)*v0P, rmsVphbin(i)*v0P,
     &           npibin
         endif
      enddo
 17   format (4(2x,g14.7), 3(2x,f8.3), 1x,g14.7, 9(2x,f8.3),  1x, I6)
      close ( 12 )
      write(*,*) 'done writing cylindrical profiles for DM. exiting...'
      return
      END


c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE CylinPSb(Vcx,Vcy,Vcz,xc,yc,zc,
     &        ZDmax,rmin,rmax,nrbin,fprofname,AxX,AxY,AxZ,tage)
c-------------------------------------------------------------------- 
c     Profiles in cylindrical coordinates,
c     where v1,v2,v3 define its axis
c     and Vcx,Vcy,Vcz,xc,yc,zc, the velocity and position of its center
c     ZDmax defines the maximum height of the sampled disk.
c     in code unis
c-------------------------------------------------------------------- 
c
c     surface brightness profile of the stellar component
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
      
       real vcx, vcy, vcz,xc,yc,zc,rmin, rmax,tage
      integer nrbin
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
      real*8 ageS, ageU, age, tyoungGyr

      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins)
      real surfBrU(0:nmaxbins),surfBrB(0:nmaxbins),surfBrV(0:nmaxbins)
     &    ,surfBrK(0:nmaxbins),surfBrbo(0:nmaxbins)
      real dmbin(0:nmaxbins)
      integer nbin(0:nmaxbins)

      character dummy
      real FlUp(nstarmax),FlBp(nstarmax),FlVp(nstarmax)
      real FlKp(nstarmax),Flbop(nstarmax)
      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      Dangular=8 ! Angular-size redshift relation for z=1 in Kpc/arcsec
      Dlum=6607.1 ! Luminosity distance at z=1 in Mpc
      Dlum=Dlum*1.d6 ! in pc
!      Dangular=1./20609 ! Angular-size relation at 10 pc in pc/arcsec 
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
!      r0Pp=r0P*1000.0 ! in pc (Physical units)
      r0arc=r0P/Dangular ! in arcsec 
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
! Read SSPs luminosities:
      open(13,file ="/home/zevs5/nayra/SSPs/SSP_Luminos_a0.500.dat")
      read(13,*) dummy
      read(13,*) dummy
      read(13,*) Nssps
      DO i=Nssps,1,-1
         read(13,100) ID,FlUp(i),FlBp(i),FlVp(i),FlKp(i),Flbop(i)
      ENDDO
100    format ( I7 , 2x, 5(1pE12.6,2x) )
!     Initialize profiles:
      do i = 0 , nmaxbins
         rbin(i) = 0.
         nbin(i) = 0
         surfBrU(i)=0.
         surfBrB(i)=0.
         surfBrV(i)=0.
         surfBrK(i)=0.
         surfBrbo(i)=0.
         dmbin(I)=0.
      enddo
      Nstars=0
      ageU = age(t)
      tyoungGyr=dble(tage)/ 1.e9 ! in Gyr
      DO ip=nsp(nspecies,1),nsp(nspecies,2) ! STARS
         istar =ip - nsp(nspecies,1) + 1
         ageS = ageU - age(dble( tbirth(istar)))
         IF (ageS.le.tyoungGyr) THEN
            Nstars=Nstars+1
            xp = x(ip) -xc
            yp = y(ip) -yc
            zp = z(ip) -zc
            rp3D2 = xp**2 + yp**2 + zp**2 
            Xpp =xp*AxX(1) +yp*AxX(2) +zp*AxX(3)
            Ypp =xp*AxY(1) +yp*AxY(2) +zp*AxY(3)
            Zpp =xp*AxZ(1) +yp*AxZ(2) +zp*AxZ(3)
            rp3D2p = xpp**2 + ypp**2 + zpp**2 
            if (abs( rp3D2 -  rp3D2p).gt.0.01) then
            write (*,*) 'error in the rotation of coordinates. Rp3D2= ',
     &           rp3D2, rp3D2p  
                STOP
             endif
             IF ( ABS(ZPP).LT.ZDmax) then
                rp =  sqrt(xpp**2+ypp**2)
                ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
                if ( ibin .le. nmaxbins ) then
                   nbin(ibin) = nbin(ibin) + 1
                   rbin(ibin) = rbin(ibin) + rp
                   surfBrU(ibin)=surfBrU(ibin) + FlUp(istar)
                   surfBrB(ibin)=surfBrB(ibin) + FlBp(istar)
                   surfBrV(ibin)=surfBrV(ibin) + FlVp(istar)
                   surfBrK(ibin)=surfBrK(ibin) + FlKp(istar)
                   surfBrbo(ibin)=surfBrbo(ibin) + Flbop(istar)
                   dmbin(ibin) = dmbin(ibin) + pw(ip)
                endif
             ENDIF
          ENDIF
       ENDDO
! Normalize profiles:
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rbin(i) = rbin(i)/  nbin(i) 
         endif
      enddo
! make profiles:
      DO i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlogmin + float(i-1)*drl)
         endif
         rr = 10.**(rlogmin + float(i)*drl)
         surf = pi * (rr**2 - rl**2)          
         if ( surf .gt. 0. ) then
            dmbin(i) = dmbin(i)/surf
            surfBrU(i)= -2.5*log10(surfBrU(i)) ! absolute mag
            surfBrU(i)=surfBrU(i) -5.0 + 5*log10(Dlum) ! aparent mag           
            surfBrU(i)=surfBrU(i)+2.5*log10(surf*r0arc**2) ! Surface brightness in mag/arcsec2

            surfBrB(i)= -2.5*log10(surfBrB(i)) ! absolute mag
            surfBrB(i)=surfBrB(i) -5.0 + 5*log10(Dlum) ! aparent mag           
            surfBrB(i)=surfBrB(i)+2.5*log10(surf*r0arc**2) ! Surface brightness in mag/arcsec2

            surfBrV(i)= -2.5*log10(surfBrV(i)) ! absolute mag
            surfBrV(i)=surfBrV(i) -5.0 + 5*log10(Dlum) ! aparent mag           
            surfBrV(i)=surfBrV(i)+2.5*log10(surf*r0arc**2) ! Surface brightness in mag/arcsec2

            surfBrK(i)= -2.5*log10(surfBrK(i)) ! absolute mag
            surfBrK(i)=surfBrK(i) -5.0 + 5*log10(Dlum) ! aparent mag           
            surfBrK(i)=surfBrK(i)+2.5*log10(surf*r0arc**2) ! Surface brightness in mag/arcsec2

            surfBrbo(i)= -2.5*log10(surfBrbo(i)) ! absolute mag
            surfBrbo(i)=surfBrbo(i) -5.0 + 5*log10(Dlum) ! aparent mag           
            surfBrbo(i)=surfBrbo(i)+2.5*log10(surf*r0arc**2) ! Surface brightness in mag/arcsec2

!            surfBr(i)= -2.5*log10(surfBr(i)) ! absolute mag
!            surfBr(i)=surfBr(i)+2.5*log10(surf*r0arc**2) ! Surface brightness in mag/arcsec2
         else
            write(*,*) 'bad surf =',surf,i, rl, rr
         endif
      ENDDO
!     Write output:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0,(yc-1.)*r0,(zc-1.)*r0,vcx*v0,vcy*v0,vcz*v0
      write(12,11) AxZ(1) , AxZ(2) , AxZ(3)
      write(12,12) rmin*r0P, rmax*r0P , nrbin , ZDmax*r0P
      write (12,13) aexpn , 1/aexpn -1.
      write(12,14) tage
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[Mpc h-1] = ', 3(g13.7,1x),
     &     '# (Vcx,Vcy,Vcz)[Km/s] = ', 3(g13.7,1x))
 11   format ('axis direction= ', 3(g13.7,1x))
 12   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5,' ZDmax= ',g12.6)
 13   format ('Distances in proper Kpc, proper velocities in Km/s and densities in Msun/pc2
     &     for a= ',f6.3, ' and redshift= ', f6.2)
 14   format ('# CYLINDRICAL PROFILES FOR STARS younger than', 
     & E8.1, ' yrs : Columns:')
 15   format ('# <r>',14x,'r_l',14x,' r_m',14x,'r_r',14x,
     &     'Surf.BrightnessU',2x,'Surf.BrightnessB',2x,
     &     'Surf.BrightnessV',2x,'Surf.BrightnessK',2x,
     &     'Surf.Brightnessbo',2x,'NPibin')
 16   format('----------------------')
      Npibin=0.
      do i = 0 , nbins
         Npibin =Npibin + nbin(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(12,17)
     &           rbin(i)*r0P, rl*r0P, 
     &           rmid*r0P, rr*r0P, 
     &           surfBrU(i),surfBrB(i),surfBrV(i),surfBrK(i),
     &           surfBrbo(i),
     &           npibin
         endif
      enddo
 17   format (4(2x,g14.7),5(4x,g14.3),2x,I6)
      close ( 12 )
      write(*,*)'done writing surface brightness profiles for stars. 
     &     exiting..'
      return 
      END
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE cooling_T_map()
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      integer Level
      real*8 Zsolar, Zsolari
      parameter ( Zsolar = 0.02 , Zsolari = 1.d0 / Zsolar )
      real rhoCell, ei
      real*8 rhog2, rhogi , rhogl, fact_nH, Tfact2, Tfact1, Tfact
      real*8 dummy1, dummy2, ai, acurr
      real*8 T_gas, coolrate, dtcool, Zdum,xion,rs
      d1 = dble(one)
      acurr = aexpn  
      rs = max ( 1.d0 / acurr - 1.0 , 0.d0 ) 
      write (*,*)  'setting Cloudy cooling rate table for z=',rs
      call Set_Cooling_Rate_rs ( rs ) 
      ai = 1.d0 / acurr
      dummy2 = ai * ai 
      a2=acurr * acurr
      r0P = r0 *1000.0 / hubble * acurr ! in Kpc (Physical units)
      r0cm=r0P*1000.0*3.0857*1.d18 ! in cm
      gamma1 = gamma -1.0
      Tfact = T_0 * gamma1 / 1.d4 ! T_0 must be in Kelvin
c.... factor converting gas density in code units to n_H in cm^-3
      fact_nH = log10(1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3)
      call Get_MaxLevelNow ()
      DO Level=MinLevel,MaxLevelNow
         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif
         
         nch = 128/2**Level
         if(nch .lt. 1) nch = 1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell, iChild, i1)
C$OMP+PRIVATE(rhoCell, ei, rhogl, rhogi, rhog2, T_gas, Zdum)
C$OMP+PRIVATE(Tfact1, Tfact2)
C$OMP+PRIVATE(coolrate, dummy1, dtcool,xion)
C$OMP+SCHEDULE(dynamic,nch)
         DO ic = 1,nLevel
            if(Level .eq. MinLevel) then
               icell = ic - 1
            else
               icell = iSelect(ic) - 1
            endif
            do iChild = 1,maxChild 
               i1 = icell + iChild
               IF ( iOctCh(i1) .eq. nil ) then ! only leaves 
                  
                  rhoCell= hvar(1,i1)
                  ei = hvar(8,i1)  
                  rhogl = log10(rhoCell) + fact_nH
                  rhogi = 1.d0 / rhoCell 
                  rhog2 = (1.d0*rhoCell ) * rhoCell
                  Tfact1 = Tfact * rhogi 
                  Tfact2 = dummy2 * Tfact1
                  T_gas = Tfact2 * ei ! in units of 10^4 K
                  Zdum = hvar(izII,i1) 
#ifdef ENRICH_SNIa
                  Zdum = Zdum + hvar(izIa,i1)
#endif
                  Zdum = max ( Zdum , 1.d-10 )
                  Zdum = log10(Zdum * rhogi * Zsolari)
                  call Cooling_Rate ( rhogl, T_gas , Zdum , coolrate ) 
                  dummy1 = -rhog2 * acurr
                  coolrate = coolrate * dummy1
!                  dtcool = abs( ei / coolrate ) 
                  dtcool =  ei / coolrate  
c
                  call  X_ion ( rhogl , T_gas*1.d4 , Zdum , rs , xion )
                  xHI=max ( zero , sngl(d1 - xion))
                  var(2,i1)=
     & log10(max(1.e-10, xHI*r0cm/2.**Level)) +rhogl !  neutral gas column density.
                  
c
c                  hvar(10,i1)=1. ! dtcool *t0*acurr**2 / 1.d9 ! in Gyr
!                  if ( dtcool.lt.0.) then
!                     var(2,i1)=log10(-dtcool *t0*a2 ) ! in log( yr)
c                  var(2,i1)=log10((1.00-xion)*10**rhogl*r0cm/2**Level) !  neutral gas column density.
!                  else
!                     var(2,i1)=-1.
!                  endif
c                  hvar(10,i1)= (hvar(1,i1)/50)*100.
c                     if ( i1.eq.13803147) then
c                        write (*,*)  rhogl , T_gas*1.d4 , Zdum ,rs ,xion 
c                     endif
c     
               ENDIF
            enddo 
         ENDDO                  !nLevel
      ENDDO                     ! Level


      end
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_ASCII(fprofname, Vxc, Vyc, Vzc,
     &     xc,yc,zc, rmax,tage)
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax,tage
      real*8 Xph,Yph,Zph
      character*256 fprofname 

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      tfact = t0*aexpn**2
      rmax2= rmax**2

      nfn  = index(fprofname, ' ') - 1
      open ( 33, file=fprofname(1:nfn) )

      write(33,300) t*tfact, aexpn
      write(33,301) (xc-1.)*r0 ,(yc-1.)*r0 ,(zc-1.)*r0
      write(33,302)

      Nstars=0
      do is=nsp(nspecies,1),nsp(nspecies,2) ! stars
         istar =is - nsp(nspecies,1) + 1
         Xph =  x(is)-xc
         Yph =  y(is)-yc
         Zph =  z(is)-zc

         IF ( Xph**2 + Yph**2 + Zph**2 .lt. rmax2) THEN
            if ( (t-tbirth(istar))*tfact.lt.tage) then
               Nstars=Nstars+1
               write(33,334) istar,
     &              (Xph)*xfact , 
     &              (Yph)*xfact ,
     &              (Zph)*xfact ,
     &              vx(is)*vfact,vy(is)*vfact,vz(is)*vfact,
     &              pw(is)*aM0 , 
     &              (t-tbirth(istar))*tfact
            endif
         ENDIF
      enddo
      write(33,*)"Nstars inside ",rmax*xfact," Kpc :",Nstars
      write(33,340) tage
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format (" X Y Z [Kpc] ,Vx  Vy  Vz [ km/s]  M[Msun] age[yr]")
 334  format(I7 , 2x, 8(e14.6,1x))
 340  format ( ' younger than ', E6.1, ' yrs')
      close ( 33 ) 
      end
c-------------------------------------------------------------------- 
      SUBROUTINE Write_DM_ASCII(fprofname,Vxc, Vyc, Vzc, xc,yc,zc, rmax)
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax
      real*8 Xph,Yph,Zph,Vxph,Vyph,Vzph
      character*256 fprofname 

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      tfact = t0*aexpn**2
      rmax2= rmax**2

      nfn  = index(fprofname, ' ') - 1
      open ( 33, file=fprofname(1:nfn) )

      write(33,300) t*tfact, aexpn
      write(33,301) (xc-1.)*r0 ,(yc-1.)*r0 ,(zc-1.)*r0
      write(33,302)
      Ndm=0
      do is=nsp(1,1),nsp(1,2) ! first DM specie
         Xph =  x(is)-xc
         Yph =  y(is)-yc
         Zph =  z(is)-zc
         Vxph= vx(is)-Vxc
         Vyph=vy(is) -vyc
         Vzph=vz(is) -vzc
         IF ( Xph**2 + Yph**2 + Zph**2 .lt. rmax2) THEN
            Ndm=Ndm+1
            write(33,334) is,
     &              Xph*xfact , 
     &              Yph*xfact ,
     &              Zph*xfact ,
     &              vxph*vfact,vyph*vfact,vzph*vfact
         ENDIF
      enddo
      write(33,*)"#DM particles inside ",rmax*xfact," Kpc :",Ndm
      write(33,340) pw(nsp(nspecies-1,1))*aM0
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format (" X Y Z [Kpc] ,Vx  Vy  Vz [ km/s]  M[Msun] age[yr]")
 334  format(I7 , 2x, 8(e14.6,1x))
 340  format ( ' with mass equal to ', E6.2, ' Msun')
      close ( 33 ) 
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
c     STELLAR ANALYSIS
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE rate(tmin,tmax,nbin, fprofname)
c-------------------------------------------------------------------- 
c     This subroutine compute #stars per each age bin
c
c     input in code units.
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      real*8 tmin,tmax,dbin,age
      real*8 aMrate
      character*256 fprofname
      parameter ( nmaxbins = 10000 )
      integer nstbin(0:nmaxbins)
      real*8 Agebin(0:nmaxbins)
      real aMassbin(0:nmaxbins), aMassaVbin(0:nmaxbins)
      real zstIIbin(0:nmaxbins)
      real aMratebin(0:nmaxbins)
!      aexpn=0.11736155
      a2 = aexpn**2
      Zsolar = 0.0199 
      aM0L=log10(aM0)
      dbin=(tmax-tmin) / float (nbin)
      nbins = int((tmax-tmin)/dbin) + 1
      write (*,*)  dbin, tmin, tmax, nbin, nbins
      ! Initialize profiles: 
      do i = 0 , nmaxbins
         nstbin(I)=0
         Agebin(i)=0.0
         aMassbin(i)=0.0
         aMassaVbin(i)=0.0
          zstIIbin(i)=0.0

          aMratebin(i)=0.0
      enddo
      nstars= nsp(nspecies,2) - nsp(nspecies,1) + 1
      write (*,*)  '# stars ', nstars
      DO i=1,Nstars
         age = tbirth(i)
         aMass = pw0(i)
         Zdum= zstII(i)
         Zdum = log10(Zdum / Zsolar)


         ibin = max ( int((age - tmin)/dbin + 1) , 0 )
         if ( ibin .le. nmaxbins ) then
            nstbin( ibin) = nstbin( ibin) + 1
            Agebin( ibin) = Agebin( ibin) + age
            aMassbin(ibin) = aMassbin(ibin) + aMass      
            zstIIbin(ibin) =  zstIIbin(ibin) + zdum
         endif
      ENDDO
      aMrate=0.000                  ! Gas consumption rate or effective star formation rate
!     Normalize profiles:
                  Nstarcheck = 0
            do i = 0 , nbins
               if ( nstbin(i) .gt. 0 ) then
                  Nstarcheck = Nstarcheck + nstbin(i)
                  Agebin(i)=Agebin(i)/nstbin(i)
                  aMrate=aMrate + aMassbin(i) 
                  aMassaVbin(i)=    aMassbin(i) / nstbin(i)
                  zstIIbin(i) =      zstIIbin(i)  / nstbin(i)

                  aMratebin(i) = aMrate !  / dbin
               endif
            enddo
            if (Nstarcheck.ne.nstars)
     &           write (*,*)  'wrong',Nstarcheck,nstars
!     Write outputfile:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,*)tmin*t0*a2/1.d9, tmax*t0*a2/1.d9, aexpn
      write(12,*) dbin*t0*a2
      write(12,*)
      write(12,15)

      do i = 0 , nbins
         if ( nstbin(i) .gt. 0 ) then
         write(12,14) i, nstbin(i) , agebin(i),  ! *t0*a2/1.d9,  
     &        log10(aMassbin(i)) +aM0L, log10(aMassaVbin(i))+aM0L,
!     &         aMratebin(i)*aM0  / (t0*a2/1.d9),
!     &          zstIIbin(i)
     &  log10 (aMratebin(i)*aM0), ! cumulative mass converted into stars.
     & aMratebin(i)/aMrate ! cumulative mass converted into stars (normalized)
      endif
      enddo
      close ( 12 )
      write(*,*) 'done writing stellar analysis. exiting...'
 14   format ( I3 , TR5, I6, TR5, F6.3, 4( TR5, F6.3) ) 
 15   format ( "bin", TR2, "#stars inside bin", TR2, "Mean age", TR2,
     & "Mass inside bin", TR2, "Average mass", TR2, "cumulative mass",
     & "cumulative mass (normalized)" )
      return
      end
c-------------------------------------------------------------------- 
      SUBROUTINE YoungSt()
c
c
c
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"
      real*8 xx, yy, zz

      do i=1,mcell
         var(2,i)=0
      enddo

      nstars = nsp(nspecies,2)-nsp(nspecies,1)+1
      tyoung = tbirth(nstars)   ! youngest star
      call Get_MaxLevelNow () 
      do  is=nsp(nspecies,1),nsp(nspecies,2) ! stars
         istar =is - nsp(nspecies,1) + 1
         IF (tbirth(istar).eq.tyoung) THEN
            xx = x(is)
            yy = y(is)
            zz = z(is)
            iC = iFindCell ( MaxLevel , xx , yy , zz )
            var(2,iC)=
     &           hvar(1,Ic)*(pw0(istar)*aM0/10000)*(100*aexpn/v0)**2 ! pressure feedback
     &           - hvar(6,Ic)
            hvar(6,iC)= 
     &           hvar(1,iC)*(pw0(istar)*aM0/10000)*(100*aexpn/v0)**2 ! pressure feedback
         ENDIF
      enddo
      end
c-------------------------------------------------------------------- 
      SUBROUTINE SF_law(filename)
c
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      character*256 filename
      real*8 ageS, ageU, age,tyoungGyr 
      real*8 xx, yy, zz
      tyoung = 3.e7 ! in yrs
      tyoungGyr =  tyoung / 1.e9   ! in Gyr
      call Get_MaxLevelNow () 
      Level=2
      CellSize(Level)=1.0 / 2**Level
      CellArea = CellSize(Level)**2
!     units:
      a3 = aexpn**3
      rho0C =  rho0 / a3 / 1.E+18 ! in M_sun/pc3 (proper units)
      r0Pk= r0*1000.0 / hubble * aexpn ! in Kpc (Physical units)
      r0Pp=r0Pk*1000.0          ! in pc (Physical units)
      CellAreaP = CellArea*r0Pk**2
      CellSizeP = CellSize(Level)*r0Pp
      do i=1,mcell
         var(2,i)=0.
      enddo

      ageU = age(t)
      DO  is=nsp(nspecies,1),nsp(nspecies,2) ! stars
         istar =is - nsp(nspecies,1) + 1
         ageS = ageU - age(dble( tbirth(istar)))
         IF (ageS.le.tyoungGyr) THEN
            xx = x(is)
            yy = y(is)
            zz = z(is)
            iC = iFindCell ( Level , xx , yy , zz )
            var(2,iC)= var(2,iC) + pw0(istar)
         ENDIF
      ENDDO
! rescale & write:
      nfpn = index ( filename , ' ' ) - 1
      open ( 12 , file = filename(1:nfpn) )
      open ( 12 , file = filename(1:nfpn) )
      write(12,*) 'IC ,sigma( M_sun/pc2) ,SFR(M_sun/yr/Kpc^2)'
      nLevel = iNOLL(Level)
      call Select_Cells ( Level , nLevel ) 
      do i1 = 1 , nLevel
         icell = iSelect(i1)
         do i2 = 0 , 7
            idc = icell + i2
            var(2,idc)= var(2,idc)*aM0/tyoung/CellAreaP ! SFR in M_sun/yr/Kpc^2
            if ( var(2,idc).gt.0.0) then
               sigma= hvar(1,idc)*rho0C*CellSizeP ! Gas Surface density in M_sun/pc2
               write (12,*)  idc , sigma , var(2,idc)
            endif
         enddo
      enddo
check:
      T_ceil= 1.0d4
      T_fact = T_0 / aexpn**2
      efact = T_ceil / T_fact / (gamma - 1.0) 
      write (*,*)  efact
      do i=1,mcell
         hvar(8,i)=efact
      enddo
      END
c-------------------------------------------------------------------- 
      SUBROUTINE  Stellar_data(filename)
c
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      character*256 filename
      real*8 ageS, ageU, age
      Zsolar = 0.0199 
      a2 = aexpn**2
      nfpn = index ( filename , ' ' ) - 1
      open ( 12 , file = filename(1:nfpn) )
      write (12,*)  'Age(Gyr)  '  ,' Initial mass(M_sun) ',
     &     'Metallicity(Mass fraction in metals): SNII, SNIa, total'
      
      nstars= nsp(nspecies,2) - nsp(nspecies,1) + 1
      ageU = age(t)
      write (*,*)  '# stars ',nstars, ' at t= ', ageU
      DO i=1,Nstars
         ageS = ageU - age(dble( tbirth(i)))
         aMass = pw0(i)
         Zdum= zstII(i) + zstIa(i)         
         write(12,13) ageS, aMass*aM0,zstII(i),zstIa(i),Zdum 
      ENDDO
 13   format (5(e10.3,3x))
      END



c-------------------------------------------------------------------- 
      SUBROUTINE distriFunc(ivar,Vmin,Vmax,nbin,filename)
c-------------------------------------------------------------------- 
c
c     This subroutine compute the distribution function of a given variable 
c
c                        ivar=1 --> Gas density
c
c     in the whole box
c
c-------------------------------------------------------------------- 
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      character*90 filename
      parameter ( nmaxbins = 1000000 )
      integer NValuebin(0:nmaxbins)
      real Valuebin(0:nmaxbins)
      Vlogmin = log10(Vmin)
      Vlmax = log10(Vmax)
      dVL= (Vlmax - Vlogmin) / nbin
      ai=1./aexpn
      a3 = aexpn**3
c.... factor converting gas density in code units to n_H in cm^-3
      fact_nH = 1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3
      fact_nHl=log10( fact_nH)
      rho0C =  rho0 / a3 / 1.E+18 ! in M_sun/pc3 (proper units)
      rho0Cl = log10(rho0C)
      ng3L = log10( 1.*ng**3)
! Initialize profiles: 
      do i = 0 , nmaxbins
         NValuebin(i)=0
         Valuebin(i)=0.0
      enddo
      N_leaves=0                       ! leaves count
      call Get_MaxLevelNow ()
c$$$      DO i1 = 1 , ncell0
c$$$         iCh1 = iOctCh(i1)
c$$$         if ( iOctCh(i1) .eq. nil ) then ! only leaves
c$$$            if (ivar.eq.1) then ! Gas Density
c$$$c               Value=log10(hvar(1,i1)) + fact_nHl ! Density in cm-3
c$$$                Value= log10(hvar(1,i1)) + rho0Cl ! M_sun/pc3
c$$$            endif
c$$$            N_leaves=N_leaves+1
c$$$            ibin = max ( int((Value - Vlogmin)/dVl) + 1 , 0 )
c$$$            if ( ibin .le. nmaxbins ) then
c$$$               NValuebin(ibin)=NValuebin(ibin)+1
c$$$               Valuebin(ibin) =    Valuebin(ibin) + Value
c$$$            ELSE
c$$$               WRITE (*,*)  'out of limits: ', ibin , value
c$$$            endif
c$$$         endif
c$$$      ENDDO
       DO Level =  MinLevel+1, MaxLevelNow 
c      DO Level =  8, MaxLevelNow 
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do i1 = 1 , nLevel
             icell = iSelect(i1)
             do i2 = 0 , 7
                idc = icell + i2
                if ( iOctCh(idc) .eq. nil ) then ! only leaves
                   if (ivar.eq.1) then ! Gas Density
                      Value =  log10(hvar(1,idc)) + fact_nHl ! Density in cm-3
c                      Value= log10(hvar(1,idc)) + rho0Cl ! M_sun/pc3
                      if (Value.lt. Vlogmin) write (*,*)  Value
                   endif
                   N_leaves=N_leaves+1
                   ibin =max( int( (Value - Vlogmin)/dVl) + 1 , 0)
                   if ( ibin .le. nmaxbins ) then
                      NValuebin(ibin)=NValuebin(ibin)+1
                      Valuebin(ibin) =    Valuebin(ibin) + Value
                   else
                      WRITE (*,*)  'out of limits: ', ibin , value
                   endif
                endif
             enddo
          enddo
       ENDDO
!     Normalize profiles:
       Ncounts=0
       do i = 0 , nmaxbins
          Ncounts=Ncounts +  NValuebin(i)
       enddo
! Write:
       nfpn = index ( filename , ' ' ) - 1
       open ( 12 , file = filename(1:nfpn) )
             write (12,*)  'Limits: ', Vlogmin, Vlmax, dvL, nbin
       IF (Ncounts.eq.N_leaves) then
          write (*,*)  'found ', Ncounts , ' leaves'
          write (*,*)  'from a total number of ',N_leaves,' leaves'
          do i = 0 , nmaxbins
             if (NValuebin(i).gt.0) then
                Valuebin(i) =    Valuebin(i) / NValuebin(i) 
                pdf = log10( 1.*NValuebin(i)/Ncounts) -dvL
                write (12,*)  i , Valuebin(i) , NValuebin(i), pdf
             endif
          enddo
       ELSE
          write (*,*)  'We are missing some leaves'
          write (*,*)  'found ', Ncounts , ' leaves'
          write (*,*)  'from a total number of ',N_leaves,' leaves'
       ENDIF
       write (*,*)  'written file: ', filename
      END
c-------------------------------------------------------------------- 
      SUBROUTINE QSOline(xc,yc,zc,Vxc,Vyc,Vzc,
     & bx,bz,x0,y0,z0,Dx,Dy,Dz,filename)
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
      write (12,6)aexpn, (xc-1.0)*r0P , (yc-1.0)*r0P, (zc-1.0)*r0P,
     & Vxc*v0P ,Vyc*v0P , Vzc*v0P 
      write (12,7)  bx*r0Pk ,  bz*r0Pk
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
c      write (*,*)  'Choice 1: Constant steps:', Tline
      i=0
      xx= x0 ! initial values
      yy= y0
      zz= z0
      ICold = -10
      DO WHILE (xx.le.aNg1.AND.yy.le.aNg1.AND.zz.le.aNg1
     &     .AND.xx.ge.0.00.AND.yy.ge.0.00.AND.zz.ge.0.00)

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
 6    format ('aexpn= ', F6.3,' Galaxy center in Mpc:', 3(TR1,F14.6)
     & , TR5, 'Peculiar velocity center in km/s:', 3(TR1,F10.3))
 7    format ('Impact parameter in Kpc: bx= ', F14.6,' bz= ',F14.6 )
 8    format ('los: r=r0+Dt: r0:'3(TR3,F14.6)'  D:',3(TR3,F14.6))
 9    format('Level',TR5,'Cell #',TR5, 'X(Mpc)',TR5,'Y(Mpc)'
     & ,TR5,'Z(Mpc)'
     & ,TR5,'Vx(Km/s)',TR5,'Vy(Km/s)',TR5,'Vz(Km/s)'
     & ,TR3,'V(Kpc3)',TR3,'NH(cm-3)',TR3,'Met(solar)',TR3,'T(K)'
     & ,TR3,'HI.Frac')
 17   format (I3,2x,I9,3(TR3,F9.6),3(TR3,F10.3),4(TR3,g7.2),TR3,E10.2 )
c      write (*,*)  'Output file: ', filename
      end
c--------------------------------------------------------------------
      SUBROUTINE Find_Starting_P( bx, by, bz, Dx, Dy, Dz,
     &         x0, y0, z0,Ng) 
c-------------------------------------------------------------------- 
c
!     Find smallest t : r=r0 + Dt
      aNg1 = float(Ng-1)

      Tline=Max(-1.e6, Min (1.e6 , (  bx - 1.) /  Dx ))
!      WRITE (*,*)Tline  
      TlineOLD=Tline
c
      Tline=Max(-1.e6, Min (1.e6 , (  bx - aNg1) /  Dx ))
!      WRITE (*,*)aNg1, Tline
      IF (abs(Tline).LT.abs(TlineOLD)) TlineOLD=Tline
c
      Tline=Max(-1.e6,Min (1.e6 , (  by - 1. ) /  Dy ))
!      WRITE (*,*)Tline  
      IF (abs(Tline).LT.abs(TlineOLD)) TlineOLD=Tline
c
      Tline=Max(-1.e6,Min (1.e6 , (  by - aNg1 ) /  Dy ))
!      WRITE (*,*)aNg1, Tline 
      IF (abs(Tline).LT.abs(TlineOLD)) TlineOLD=Tline
c
      Tline=Max(-1.e6,Min (1.e6 , (  bz - 1.) /  Dz ))
!      WRITE (*,*)Tline  
      IF (abs(Tline).LT.abs(TlineOLD)) TlineOLD=Tline
c
      Tline=Max(-1.e6,Min (1.e6 , (  bz - aNg1) /  Dz ))     
!      WRITE (*,*)aNg1, Tline  
      IF (abs(Tline).LT.abs(TlineOLD)) TlineOLD=Tline
c
!      write (*,*)  'smallest t: ', TlineOLD
c
      x0 = bx - Dx * TlineOLD
      y0 = by - Dy * TlineOLD
      z0 = bz - Dz * TlineOLD
      IF (x0.LT.1.OR.x0.GT.aNg1) write (*,*)  'wrong',x0
      IF (y0.LT.1.OR.y0.GT.aNg1) write (*,*)  'wrong',y0
      IF (z0.LT.1.OR.z0.GT.aNg1) write (*,*)  'wrong',z0

      IF ( TlineOLD.lt.0. ) THEN
         write (*,*)  'Changed sign of direction'
         Dx= -Dx
         Dy= -Dy
         Dz= -Dz
      ENDIF
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE Mock_distributions(iopt,xc,yc,zc,rmax,Vxc,Vyc,Vzc)
c-------------------------------------------------------------------- 
c
c     This subroutine generate Mock distributions of densities
c     and velocities using a given mesh and a center.
c
c     iopt--> 1 --> Density distribution 1/r^2
c
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax,Vxc,Vyc,Vzc
      
      rmax2=rmax**2
      rcellOld=0.0
      i_alpha=3  !  Density distribution 1/r^i_alpha

      call Get_MaxLevelNow ()
      DO Level =  MinLevel+1, MaxLevelNow 
         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel ) 
         do i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idc = icell + i2
               if ( iOctCh(idc) .eq. nil ) then ! only leaves
                  
                  call  Ps ( idc , Posx,Posy,Posz )
                  rcell2 = (Posx - xc)**2 + (Posy-yc)**2+(Posz-zc)**2
                  IF (rcell2.lt.rmax2) THEN
                        
                     rcell=sqrt(rcell2)
                     if (rcell.gt.0.0) then
                        IF (IOPT.EQ.1) THEN !  Density distribution 1/r^i_alpha
                           hvar(1,idc)= 1.0/rcell**i_alpha
                          
                        ELSEIF (IOPT.EQ.2) THEN ! Constant radial flow
                           vxCell=(Posx-Xc)/rcell + Vxc
                           vyCell=(Posy-Yc)/rcell + Vyc
                           vzCell=(Posz-Zc)/rcell + Vzc

                            hvar(3,idc)=vxCell*hvar(1,idc)
                            hvar(4,idc)=vyCell*hvar(1,idc)
                            hvar(5,idc)=vzCell*hvar(1,idc)
                           
                        ENDIF
                     else
                        if (rcellOld.le.0.0) STOP
                        IF (IOPT.EQ.1) THEN !  Density distribution 1/r^i_alpha
                           hvar(1,idc)= 1.0/rcellOld**i_alpha
                        ELSEIF (IOPT.EQ.2) THEN ! Constant radial flow
                            hvar(3,idc)=Vxc*hvar(1,idc)
                            hvar(4,idc)=Vyc*hvar(1,idc)
                            hvar(5,idc)=Vzc*hvar(1,idc)

                        ENDIF
                     endif 
                     rcellOld=rcell
                  ELSE
                     hvar(1,idc)=0.0
                     hvar(3,idc)=0.0
                     hvar(4,idc)=0.0
                     hvar(5,idc)=0.0
                  ENDIF
               endif
            enddo
         enddo
      ENDDO
      END
