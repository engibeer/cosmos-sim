c     =====================================================================
c                                                                         
c     ART : Radiation pressure 
c                                                                         
c     by Daniel Ceverino, Jerusalem, January 2010    
c
c     =====================================================================
c

c     -------------------------------------------
      subroutine RadPre_Compu ( Level )
c     -------------------------------------------
c 
c     compute RadPre values 
c
      integer Level
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      real*8 t_fb_code,t_rp_code,cfbe,heatfact, tfbc

c      Parameter (ColumnTh=21.00) ! minimum column density for radiative feedback (in cm^-2)
      Parameter (ColumnTh=0.00) ! minimum column density for radiative feedback (in cm^-2)
      Parameter (GammaPFactor= 230.0) ! 8 ! factor times standard stellar feedback
      Parameter (t_rp= 5.0e6)  ! duration of RP feedback in yrs
      Parameter (boost = 11.0) ! due to tau_IR or due to efficiency of SF
      real*8 RadPSum

      logical followstar

      ai = 1.d0 / aexp(Level) 

 
c.... factor converting gas density in code units to n_H in cm^-3
      fact_nH = log10(1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3)

      nstar1 = nsp(nspecies,1) - 1 

      t_fb_code = t_fb / t0 / aexp(Level)**2
      t_rp_code = t_rp / t0 / aexp(Level)**2   ! duration of RP feedback
      t_fbi_code= 1.0 / t_fb_code
      cfbe = C_fb * aexp(Level)**2 
      heatfact = cfbe * t_fbi_code

      RadPreF=v0/3.e5 / aexp(Level) / ( CellSize(Level)/2)**2

      SizePl=log10(CellSize(Level)*r0/hubble*aexp(Level)*3.0857e24) ! in log cm

      r0Pk= (r0*1000.0 / hubble) * aexp(Level) ! in Kpc (Physical units)
      v0P = v0 / aexp(Level)        ! in km/s (proper units)
      P0cgs= 1.2643d-19/(aexp(Level)**5)/1.380e-16  ! P/k

c      star_mass =0
c      star_age =0
c      radpre_max = 1.0e20/P0cgs   !5.0e7/P0cgs
  
 
c------------------------
      followstar = .false.
c
c      if (followstar) then
c         idstar1 = 9542647
c         idstar2 = 9538452
c         idstar3 = 9542312
c         idstar4 = 9539143
c      endif
c---------------------------

      call Get_MaxLevelNow ()

      IF ( Level .eq. MinLevel ) THEN 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(idummy,istar,tp,tfbc,rhogl,heatRateSNII)
         do i1 = 1 , ncell0
            if ( iOctCh(i1) .eq. nil ) then ! only leaves 

               if ( iCL(i1) .gt. nil ) then
                  idummy = iCL(i1) ! iparticle
                  do while ( idummy .ne. nil )
                     if ( idummy .ge. nsp(nspecies,1) ) then ! a star
                        istar=idummy - nstar1
                        tp  = pt(idummy)                     
                        tfbc = tp - tbirth(istar) ! stellar age
c===========STG: changed time cutoff:
                        if ( tfbc .le. t_rp_code ) then !  star younger than t_rp
c============       
                           rhogl = log10(hvar(1,i1)) + fact_nH ! rhog_code -> log10(n_H), [n_H]=[cm^-3]
c                           if (rhogl+sizePl.gt.ColumnTh) then ! column density limit for radiation feedback
                           if (rhogl .gt. -2.0) then ! volume density limit to avoid too large soundspeed
                              heatRateSNII= heatfact* pw0(istar) ! SN feedback
                              RadPre(i1) = RadPre(i1) + 
     &                             GammaPFactor*heatRateSNII*RadPreF ! radiation pressure in that cell
c            if (RadPSum.gt.0.00) write(25,*) i1, rhogl+sizePl, rhogl, 
c     &          sizePl,RadPre(i1)
                           endif
                        endif
                     endif
                     idummy         = iLL(idummy,1) ! next particle in the cell
                  enddo
               endif
            endif
         enddo

      ELSE

         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel ) 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,ici0,iChild,i1,ici)
C$OMP+PRIVATE(idummy,istar,tp,tfbc,rhogl,heatRateSNII) 
C$OMP+PRIVATE(Posx,Posy,Posz,star_mass) 
         do ic = 1 , nLevel
            icell = iSelect(ic) - 1
            ici0 = (ic-1)*nchild
            do iChild = 1 , 8 
               i1 = icell + iChild
               ici = ici0 + iChild

               if ( iOctCh(i1) .eq. nil ) then ! only leaves
                  if ( iCL(i1) .gt. nil ) then
                     idummy = iCL(i1) ! iparticle
                     do while ( idummy .ne. nil )
                        if ( idummy .ge. nsp(nspecies,1) ) then ! a star
                           istar=idummy - nstar1
                           tp  = pt(idummy)                     
                           tfbc = tp - tbirth(istar) 
c=============== STG changed time cutoff:
                           if ( tfbc .le. t_rp_code ) then !  star younger than t_rp
c================STG:
               if ( tfbc < 2.e4/t0/aexp(Level)**2 ) then
                  delay = 0.1
               else if ( tfbc < 4.e4/t0/aexp(Level)**2 ) then
                  delay = 0.2
               else if ( tfbc < 6.e4/t0/aexp(Level)**2 ) then
                  delay = 0.5
               else 
                  delay = 1.0
               endif
c======================================
               rhogl = log10(hvar(1,i1)) + fact_nH ! rhog_code -> log10(n_H), [n_H]=[cm^-3]
c              if (rhogl+sizePl.gt.ColumnTh) then ! column density limit for radiation feedback
               if (rhogl .gt. -4.0) then ! volume density limit to avoid too large soundspeed
                  smooth = 0.5*( 1 + erf((rhogl + 1.0)/(1.414*0.5)) )
                  heatRateSNII= heatfact* pw0(istar) ! SN feedback
                  RadPre(i1) = RadPre(i1) + 
     &          smooth*delay*boost*GammaPFactor*heatRateSNII*RadPreF ! radiation pressure in that cell
c--------------------------------
c                 RadPre(i1) = min( RadPre(i1) , radpre_max )  ! cap on radpressure to avoid large temp
c--------------------------------
c             write(25,26) iLv(ic), ic, istar, rhogl+sizePl, rhogl, 
c     &          sizePl,RadPre(ic), heatRateSNII, RadPreF, pw0(istar)*aM0

               endif
c-----------------------------------
c                              star_count = star_count + 1 
c                              mass_tot = mass_tot + pw0(istar)*aM0
c
c---------------------------------- only rad pre around particles of interest:
c$$$       if (i1 .ne. iFindCell(10,x(idstar1),y(idstar1),z(idstar1)) .and.
c$$$     &      i1 .ne. iFindCell(10,x(idstar2),y(idstar2),z(idstar2)) .and.
c$$$     &      i1 .ne. iFindCell(10,x(idstar3),y(idstar3),z(idstar3)) .and.
c$$$     &      i1 .ne. iFindCell(10,x(idstar4),y(idstar4),z(idstar4))) then
c$$$                    RadPre(i1) = 0
c$$$       endif
c-------------------------------- add rad_pre to neighbors:
c                 do idir=1,6
c                    RadPre(iNb(i1,idir)) = RadPre(i1)/4.0 
c                 enddo 
c--------------------------------
c---------------------------------- this fixes the star particle to the gas parcel:
c$$$                              if (idummy == idstar1 .or. 
c$$$     &                             idummy == idstar2 .or.
c$$$     &                             idummy == idstar3 .or.
c$$$     &                             idummy == idstar4) then
c$$$                                 vx(idummy) = hvar(3,i1) / hvar(1,i1)
c$$$                                 vy(idummy) = hvar(4,i1) / hvar(1,i1)
c$$$                                 vz(idummy) = hvar(5,i1) / hvar(1,i1) 
c$$$                              endif
c-------------------------------------------------------------------

c-----------properties of all young stars in cold cells within a specified radius:
c
c$$$                              call Ps (i1,Posx,Posy,Posz)
c$$$                              vx0 =    hvar(3,i1) / hvar(1,i1)
c$$$                              vy0 =    hvar(4,i1) / hvar(1,i1)
c$$$                              vz0 =    hvar(5,i1) / hvar(1,i1)
c$$$c
c$$$                  if (Level == 10 .and. nstep == 1 .or.
c$$$     &                Level == 9  .and. nstep == 2   ) then
c$$$                     if (hvar(6,i1)*P0cgs/(10.**rhogl) < 10000.0 .and. 
c$$$     &                   tfbc*t0*aexp(Level)**2 < 5.e6 .and.
c$$$     &                  (Posx*r0Pk-1659.37)**2+(Posy*r0Pk-1543.71)**2 
c$$$     &                   +(Posz*r0Pk-1794.99)**2 < 10.**2) then
c$$$                        write(28,26) i1, ic, iLv(i1), 
c$$$     &                       Posx*r0Pk , Posy*r0Pk , Posz*r0Pk,
c$$$     &                       rhogl , hvar(6,i1)*P0cgs ,
c$$$     &                       RadPre(i1)*P0cgs ,
c$$$     &                    pw0(istar)*aM0, (tfbc*t0*aexp(Level)**2)/1.e6
c$$$                        write(30,'(i10,11g16.6)') 
c$$$     &                       idummy,
c$$$     &                       vx0*v0P ,vy0*v0P, vz0*v0P, 
c$$$     &                     pw0(istar)*aM0,(tfbc*t0*aexp(Level)**2)/1.e6,
c$$$     &                     x(idummy)*r0Pk,y(idummy)*r0Pk,z(idummy)*r0Pk,
c$$$     &                     vx(idummy)*v0P,vy(idummy)*v0P,vz(idummy)*v0P
c$$$                     endif
c$$$                  endif
c
c--------------------------------------------------------------------------                             
                           endif ! end age condition

                        endif ! end star condition
                        idummy         = iLL(idummy,1) ! next particle in the cell
                     enddo
                  endif
c----------------------------- properties of all cells at given levels:
c
c$$$                  if (Level == 10 .and. nstep == 1 .or.
c$$$     &                Level == 9  .and. nstep == 2   ) then
c$$$                     write(27,26) i1, ic, iLv(i1), 
c$$$     &                    Posx*r0Pk , Posy*r0Pk , Posz*r0Pk,
c$$$     &                    rhogl , hvar(6,i1)*P0cgs ,
c$$$     &                    RadPre(i1)*P0cgs ,
c$$$     &                    vx0*v0P ,vy0*v0P, vz0*v0P, 
c$$$     &                    star_mass, star_age/1.e6,
c$$$     &                    star_vx, star_vy, star_vz 
c$$$                  endif
c
c--------------------------------------------------------------
c                  write(50,*) i1, star_count, mass_tot
c                  star_count =0
c                  mass_tot =0
c-------------------------------------
               endif
            enddo
         enddo
      ENDIF
 26    format ( 3(I9), 3(g16.6), f10.4, 7(g16.6))
c      write (*,*)  '----------------- done L=',Level, 'step=',nstep



c----------------------------------------------------------------------------------

      if ( followstar ) then

      write (*,*)  '----------------- done L=',Level, 'step=',nstep


c- print all properties of a single gas parcel containing a star across cells in time: 

         open(40,file='a.dat',access='append')
         open(41,file='b.dat',access='append')
         open(42,file='c.dat',access='append')
         open(43,file='d.dat',access='append')

         idstar = idstar1   ! idstar is the same as idummy above

         idstarcell  = 
     &        iFindCell(10,x(idstar),y(idstar),z(idstar)) 

       if ( Level == iLv( idstarcell ) ) then

         star_mass = pw0(idstar-nstar1)*aM0
         tp  = pt(idstar)                     
         tfbc = tp - tbirth(idstar-nstar1) 
         star_age= tfbc*t0*aexp(Level)**2
         star_x = x(idstar)*r0Pk
         star_y = y(idstar)*r0Pk
         star_z = z(idstar)*r0Pk
         star_vx = vx(idstar)*v0P
         star_vy = vy(idstar)*v0P
         star_vz = vz(idstar)*v0P
         
         
         i1 = idstarcell
         
         write(*,'(a,2i10,3f12.3)') 'star:', 
     &        idstar,  idstarcell, 
     &        star_vx, star_vy, star_vz
         
         call Ps (i1,Posx,Posy,Posz)
         vx_g =    hvar(3,i1) / hvar(1,i1)
         vy_g =    hvar(4,i1) / hvar(1,i1)
         vz_g =    hvar(5,i1) / hvar(1,i1)
         
         rhogl = log10(hvar(1,i1)) + fact_nH
         
         write(*,'(a,i10,2g16.6,es12.4,f10.4,6f12.3,f14.6,
     &        i5,f12.6)') 'cell :', i1, 
     &        star_mass, RadPre(i1)*P0cgs, 
     &        hvar(6,i1)*P0cgs,
     &        rhogl,
     &        vx_g*v0P, vy_g*v0P, vz_g*v0P,
     &        Posx*r0Pk , Posy*r0Pk , Posz*r0Pk, 
     &        star_age/1.e6 , iLv(i1), 
     &        CellSize(iLv(i1))*r0Pk*1000.
         do idir=1,6
            vx0 =hvar(3,iNb(i1,idir)) / hvar(1,iNb(i1,idir))
            vy0 =hvar(4,iNb(i1,idir)) / hvar(1,iNb(i1,idir))
            vz0 =hvar(5,iNb(i1,idir)) / hvar(1,iNb(i1,idir))
            write(*,'(i2,i10,2es12.4,2f10.3,i4)')
     &           idir, iNb(i1,idir),
     &           RadPre(iNb(i1,idir))*P0cgs, 
     &           hvar(6,iNb(i1,idir))*P0cgs, 
     &           log10(hvar(1,iNb(i1,idir))) + fact_nH,
     &           sqrt((vx0*v0P)**2+(vy0*v0P)**2+(vz0*v0P)**2),
     &           iLv(iNb(i1,idir))
         enddo
         
         write(40,'(2i10,i5,2g16.6,es12.4,f10.4,6f12.3,f14.6,i5)')
     &        idstar, i1, iLv(i1), 
     &        star_mass, RadPre(i1)*P0cgs, 
     &        hvar(6,i1)*P0cgs, rhogl,
     &        vx_g*v0P, vy_g*v0P, vz_g*v0P,
     &        Posx*r0Pk , Posy*r0Pk , Posz*r0Pk, 
     &        star_age/1.e6, nstep 
         
       endif
         
c     ---------- another star:
         
         idstar = idstar2   ! idstar is the same as idummy above

         idstarcell  = 
     &        iFindCell(10,x(idstar),y(idstar),z(idstar)) 

       if ( Level == iLv( idstarcell ) ) then

         star_mass = pw0(idstar-nstar1)*aM0
         tp  = pt(idstar)                     
         tfbc = tp - tbirth(idstar-nstar1) 
         star_age= tfbc*t0*aexp(Level)**2
         star_x = x(idstar)*r0Pk
         star_y = y(idstar)*r0Pk
         star_z = z(idstar)*r0Pk
         star_vx = vx(idstar)*v0P
         star_vy = vy(idstar)*v0P
         star_vz = vz(idstar)*v0P
         
         
         i1 = idstarcell
         
c         write(*,'(a,2i10,3f12.3)') 'star:', 
c     &        idstar,  idstarcell, 
c     &        star_vx, star_vy, star_vz
         
         call Ps (i1,Posx,Posy,Posz)
         vx_g =    hvar(3,i1) / hvar(1,i1)
         vy_g =    hvar(4,i1) / hvar(1,i1)
         vz_g =    hvar(5,i1) / hvar(1,i1)
         
         rhogl = log10(hvar(1,i1)) + fact_nH
         

         
         write(41,'(2i10,i5,2g16.6,es12.4,f10.4,6f12.3,f14.6,i5)')
     &        idstar, i1, iLv(i1), 
     &        star_mass, RadPre(i1)*P0cgs, 
     &        hvar(6,i1)*P0cgs, rhogl,
     &        vx_g*v0P, vy_g*v0P, vz_g*v0P,
     &        Posx*r0Pk , Posy*r0Pk , Posz*r0Pk, 
     &        star_age/1.e6, nstep 
    
       endif
         
c     ---------- another star:

         idstar = idstar3   ! idstar is the same as idummy above

         idstarcell  = 
     &        iFindCell(10,x(idstar),y(idstar),z(idstar)) 

       if ( Level == iLv( idstarcell ) ) then

         star_mass = pw0(idstar-nstar1)*aM0
         tp  = pt(idstar)                     
         tfbc = tp - tbirth(idstar-nstar1) 
         star_age= tfbc*t0*aexp(Level)**2
         star_x = x(idstar)*r0Pk
         star_y = y(idstar)*r0Pk
         star_z = z(idstar)*r0Pk
         star_vx = vx(idstar)*v0P
         star_vy = vy(idstar)*v0P
         star_vz = vz(idstar)*v0P
         
         
         i1 = idstarcell
         
c         write(*,'(a,2i10,3f12.3)') 'star:', 
c     &        idstar,  idstarcell, 
c     &        star_vx, star_vy, star_vz
         
         call Ps (i1,Posx,Posy,Posz)
         vx_g =    hvar(3,i1) / hvar(1,i1)
         vy_g =    hvar(4,i1) / hvar(1,i1)
         vz_g =    hvar(5,i1) / hvar(1,i1)
         
         rhogl = log10(hvar(1,i1)) + fact_nH
         

         
         write(42,'(2i10,i5,2g16.6,es12.4,f10.4,6f12.3,f14.6,i5)')
     &        idstar, i1, iLv(i1), 
     &        star_mass, RadPre(i1)*P0cgs, 
     &        hvar(6,i1)*P0cgs, rhogl,
     &        vx_g*v0P, vy_g*v0P, vz_g*v0P,
     &        Posx*r0Pk , Posy*r0Pk , Posz*r0Pk, 
     &        star_age/1.e6, nstep 

       endif
         
c     ---------- another star:

         idstar = idstar4   ! idstar is the same as idummy above

         idstarcell  = 
     &        iFindCell(10,x(idstar),y(idstar),z(idstar)) 

       if ( Level == iLv( idstarcell ) ) then

         star_mass = pw0(idstar-nstar1)*aM0
         tp  = pt(idstar)                     
         tfbc = tp - tbirth(idstar-nstar1) 
         star_age= tfbc*t0*aexp(Level)**2
         star_x = x(idstar)*r0Pk
         star_y = y(idstar)*r0Pk
         star_z = z(idstar)*r0Pk
         star_vx = vx(idstar)*v0P
         star_vy = vy(idstar)*v0P
         star_vz = vz(idstar)*v0P
         
         
         i1 = idstarcell
         
c         write(*,'(a,2i10,3f12.3)') 'star:', 
c     &        idstar,  idstarcell, 
c     &        star_vx, star_vy, star_vz
         
         call Ps (i1,Posx,Posy,Posz)
         vx_g =    hvar(3,i1) / hvar(1,i1)
         vy_g =    hvar(4,i1) / hvar(1,i1)
         vz_g =    hvar(5,i1) / hvar(1,i1)
         
         rhogl = log10(hvar(1,i1)) + fact_nH
         

         
         write(43,'(2i10,i5,2g16.6,es12.4,f10.4,6f12.3,f14.6,i5)')
     &        idstar, i1, iLv(i1), 
     &        star_mass, RadPre(i1)*P0cgs, 
     &        hvar(6,i1)*P0cgs, rhogl,
     &        vx_g*v0P, vy_g*v0P, vz_g*v0P,
     &        Posx*r0Pk , Posy*r0Pk , Posz*r0Pk, 
     &        star_age/1.e6, nstep 

       endif

         close(40)
         close(41)
         close(42)
         close(43)
c--------------------------------------------------------------------------
         
      
      endif

c$$$      if (Level < 10) then
c$$$         stop
c$$$      endif


      end

   

    



 


