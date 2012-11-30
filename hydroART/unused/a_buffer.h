c    include file for the previous step
c
c    Daniel Ceverino 2005
c

c
c..... Gas data
c
      parameter ( n0 = 2 + nhvar)
      parameter ( ncell0n0 = n0 *ncell0 )
      parameter (moct13 = 13*moct)
      parameter (nctotn0 = n0*nctot)

       integer iB0(4),            ! istep, MaxLevelNow, ainit, iOctFree, nOct
     & iSObu(MinLevel:MaxLevel),
     & iB1(ncell0),               ! iOctCh(ncell0)
     & iB2(MinLevel:MaxLevel+1),  ! iNOLL
     & iB3(MinLevel:MaxLevel+1),  ! iHOLL
     & iB4(MinLevel+1:MaxLevel,moct13),  ! Pointers: iOctPs(3) , iOctNb(6), iOctPr, iOctLv, iOctLL1, iOctLL2
     & iB5(0:nctot) ! iOctCh
       real ainitbu , 
     &                    B1(ncell0n0) ,  ! hydro data of Level 0
     &                    B2(nctotn0)        ! hydro data of L>0
      real*8                 B0(3),                           ! t , dt, adum
     &                   tlbu(MinLevel:MaxLevel)      ,  
     &                   tloldbu(MinLevel:MaxLevel)   ,   
     &                   dtlbu(MinLevel:MaxLevel)     ,   
     &                   dtloldbu(MinLevel:MaxLevel)  

        common /Bu00/ iSObu
        common /Bu0/   iB0, ainitbu
        common /Bu1/   B0 , 
     &                   tlbu      ,  
     &                   tloldbu   ,   
     &                   dtlbu     ,   
     &                   dtloldbu
        common /Bu2/  iB1, iB2 , IB3
        common /Bu3/ B1
        common /Bu4/ iB4
        common /Bu5/ iB5
        common /Bu6/ B2
c
c..........   Particles data
c
           real xbuff(npmax) , ybuff(npmax), zbuff(npmax),
     &          vxbuff(npmax), vybuff(npmax) , vzbuff(npmax)

         real extrasB(100), AexpnB, AstepB
         integer istepB
         common /Bu7/ xbuff, ybuff, zbuff, vxbuff, vybuff, vzbuff
         common /Bu8/ extrasB, AexpnB, AstepB, istepB
c
c.....   Particles time
c
         real pdtbuff(npmax)
         common /Bu9/ pdtbuff
c
c....   Stars data
c
        integer nstarsB
        real*8  ws_oldB , ws_oldiB
        real pwsB(nstarmax),
     &          pw0B(nstarmax),
     &          tbirthB(nstarmax), 
     &          zstIIB(nstarmax), 
     &          zstIaB(nstarmax)
        common /Bu10/  ws_oldB , ws_oldiB
        common /Bu11/ pwsB, pw0B, tbirthB, zstIIB, zstIaB
        common /Bu12/ nstarsB
