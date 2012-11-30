c     ---------------------------
c     Parameters for AGN feedback
c     ---------------------------
 
c 
c     .............AGN feedback............... 
c 
      integer nmaxAGN
      parameter ( nmaxAGN = 10000 )
      integer  nAGN, iAGNstep
      integer  iAGN_ON(nmaxAGN), iJet(nmaxAGN), LevAGN_old(nmaxAGN)
      integer  istepAGN(nmaxAGN)
      double precision xAGN(nmaxAGN), yAGN(nmaxAGN), zAGN(nmaxAGN)
      real*8  vxAGN(nmaxAGN), vyAGN(nmaxAGN), vzAGN(nmaxAGN)
      real*8  amhAGN(nmaxAGN), rhAGN(nmaxAGN), aexpnAGN(nmaxAGN)
      real*8  amtotAGN(nmaxAGN),amstAGN(nmaxAGN)
      real*8  amgasAGN(nmaxAGN),amcoldAGN(nmaxAGN)
      real*8  E61(nmaxAGN), t7(nmaxAGN), fke(nmaxAGN)
      real*8  tAGN_now(nmaxAGN), eAGN_now(nmaxAGN)
      real*8  vmaxAGN(nmaxAGN)
 
      common /AGN01/ nAGN, iAGNstep
      common /AGN02/ iAGN_ON, iJET, LevAGN_old, istepAGN
      common /AGN03/ xAGN, yAGN, zAGN
      common /AGN04/ vxAGN, vyAGN, vzAGN 
      common /AGN05/ amhAGN, rhAGN, aexpnAGN
      common /AGN06/ amtotAGN, amstAGN, amgasAGN, amcoldAGN
      common /AGN07/ E61, t7, fke
      common /AGN08/ tAGN_now, eAGN_now, vmaxAGN
