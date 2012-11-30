c
c     header file of the initial conditions for a galaxy simulation
c
c adapted from hstart.h
c
c      Daniel Ceverino (2005)

      PARAMETER (Lblock=16,Nblocks=NROW/Lblock) ! For mass refinements

      DIMENSION         wspecies(10),lspecies(10)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
