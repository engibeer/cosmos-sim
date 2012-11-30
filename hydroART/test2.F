#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      include 'hstart_Galaxy.h'


      open (3, file='PMcrdIC.DAT', form = 'unformatted' )
      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     &                  Ocurv,extras
c
      Omb0 = 0.022/(hubble)**2
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRIDC,NRECL,Om0,Oml0,Omb0,hubble
 100    format (1X,'Header=>',A45,/
     &           1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &           1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &           1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &           1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &               ' Omegab_0=',F7.3,' Hubble=',f7.3)
c

      end
