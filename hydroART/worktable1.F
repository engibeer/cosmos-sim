
      SUBROUTINE DelayFeedback()
      logical continue_heat
      real*8 dtheat, theat
C$OMP+PRIVATE ( continue_heat, dtheat, theat)

#ifdef FEEDBACK
c
c             introducing a delay in the energy release of t_fb years ( D. Ceverino, 2006)
c
                  if ( tfbc .le. t_fb_code ) then 
                     continue_heat= .TRUE.
                     dtheat= tbirth(istar) + t_fb_code -tp ! 1 step
                     dte = dtheat 
                     theat = tp +  dtheat ! first step
                     DO WHILE ( continue_heat) 
                        IF ( theat .gt. t_next) THEN
!                           dte = dtheat -( theat-t_next)
                           continue_heat= .FALSE.
                        ELSE    ! For a delay in energy release
!                           efrac =  dte * t_fbi_code ! delta_x / t_fb_code
                           efrac = 1. ! For a delay in energy release
#ifdef ENRICH
                           fmetej = efrac * fmass_met * pw0(istar)
                           hvar(izII,idcell) = hvar(izII,idcell) + 
     &                          fmetej * CellVoli(Level)
#endif
                           erel = efrac * cfbe * pw0(istar)
c
c....               don't allow to explode in an already hot bubble
c                   and limit energy release 
c
                           erel = min ( erel , 
     &                          efact * hvar(1,idcell) )
                           if ( edum .lt. T_ceil ) then  
                              hvar(2,idcell) = hvar(2,idcell) + erel 
                              hvar(8,idcell) = hvar(8,idcell) + erel 
                            hvar(6,idcell) = dble((hvar(7,idcell)-1.0)*
     &                             hvar(8,idcell))
CEVERINO01032006 
C                      hvar(10,idcell)=(hvar(1,idcell)/50.)*erel  ! Mark sf regions.
                           endif
                        ENDIF
                        theat =  theat + dtheat
              IF (theat -dtheat/2. .gt. tbirth(istar) + t_fb_code) 
     &                       continue_heat= .FALSE.  ! end of feedback
                     ENDDO
                  endif
#endif                     
