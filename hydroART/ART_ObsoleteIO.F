c     -----------------------------------
      subroutine Read_Particles_Binary_old ()
c     -----------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
      include 'a_tree.h'
      include 'a_control.h'
 
      character  Hd*5,Tail*4,Nm*40

      Hd    = 'C3crs'
      Tail  = '.DAT'

      ngrid = ng

      open ( 2 , file = 'Result.DAT' , status = 'unknown' )
      open ( 3 , file = 'C3CRD.DAT ', form = 'unformatted', 
     &            status = 'old'                           )

c.... read control information and check whether it has proper structure

      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,hubble
      write (2,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,hubble
 100  format (1X,'Header=>',A45,/
     &            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3)

      if ( nrow .ne. nrowc ) then
        write (*,*) ' NROW in PARAMETER and in TAPE-FILE are different'
        stop
      endif
      if ( ngrid .ne. ngridc ) then
        write (*,*)
     &           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      endif

c.... open particle file on disk
c     for Exemplar, SP: nacces = nbyte; for PC: nacces = nrecl 

      nbyte  = nrecl * 4
      nacces = nbyte/4 
 
      open ( 1 , file = 'C3crs0.DAT', access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3

      ip = 0

      do irow = 1 , nrow
        call GetRow ( irow , 1 )
	do in = 1 , npage
          ip    = ip + 1
	   x(ip) = xpar(in)
	   y(ip) = ypar(in)
	   z(ip) = zpar(in)
	  vx(ip) = vxx(in)
	  vy(ip) = vyy(in)
	  vz(ip) = vzz(in)
         enddo
      enddo

      close (1)
      close (2)
      close (3)

      return
      end
c     -----------------------------------
      subroutine Read_Particles_Binary_old2 ()
c     -----------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
c 
      character  Hd*5,Tail*4,Nm*40

      Hd    = 'PMcrs'
      Tail  = '.DAT'

      ngrid = ng

      open ( 2 , file = 'Result.DAT' , status = 'unknown' )
      open ( 3  , file = 'SB2/PMcrda0.988.DAT ', form = 'unformatted', 
     &            status = 'old'                           )

c.... read control information and check whether it has proper structure

      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                  ,Ocurv,  extras
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,Omb0,hubble
      write (2,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Omb_0=',f7.3,
     &            1x,' H0=',F7.3)

      if ( nrowreal .ne. nrowc ) then
        write (*,*) ' NROW in PARAMETER and in TAPE-FILE are different'
        stop
      endif
      if ( ngrid .ne. ngridc ) then
        write (*,*)
     &           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      endif

c.... open particle file on disk
c     for Exemplar, SP: nacces = nbyte; for PC: nacces = nrecl 

      nbyte  = nrecl * 4
      nacces = nbyte
 
      open ( 1 , file = 'SB2/PMcrs0a0.988.DAT', access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3

      ip = 0
      nrows = nrow**3/nrowreal**2

      do irow = 1 , nrows
        call GetRow ( irow , 1 )
	do in = 1 , npage
          ip    = ip + 1
	   x(ip) = xpar(in)
	   y(ip) = ypar(in)
	   z(ip) = zpar(in)
	  vx(ip) = vxx(in)
	  vy(ip) = vyy(in)
	  vz(ip) = vzz(in)
         enddo
      enddo

      close (1)
      close (2)
      close (3)

      return
      end
c     -----------------------------------
      subroutine Read_Particles_Binary_old22 ()
c     -----------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
c 
      character  Hd*5,Tail*4,Nm*40

      Hd    = 'PMcrs'
      Tail  = '.DAT'

      ngrid = ng

      open ( 2 , file = 'Result.DAT' , status = 'unknown' )
      open ( 3  , file = 'SB2/PMcrd_i.DAT ', form = 'unformatted', 
     &            status = 'old'                           )

c.... read control information and check whether it has proper structure

      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                  ,Ocurv, extras
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,Omb0,hubble
      write (2,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Omb_0=',f7.3,
     &            1x,' H0=',F7.3)

      if ( nrowreal .ne. nrowc ) then
        write (*,*) ' NROW in PARAMETER and in TAPE-FILE are different'
        stop
      endif
      if ( ngrid .ne. ngridc ) then
        write (*,*)
     &           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      endif

c.... open particle file on disk
c     for Exemplar, SP: nacces = nbyte; for PC: nacces = nrecl 

      nbyte  = nrecl * 4
      nacces = nbyte
 
      open ( 1 , file = 'SB2/PMcrs0_i.DAT', access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3

      ip = 0
      nrows = nrow**3/nrowreal**2

      do irow = 1 , nrows
        call GetRow ( irow , 1 )
	do in = 1 , npage
          ip    = ip + 1
	   x(ip) = xpar(in)
	   y(ip) = ypar(in)
	   z(ip) = zpar(in)
	  vx(ip) = vxx(in)
	  vy(ip) = vyy(in)
	  vz(ip) = vzz(in)
         enddo
      enddo

      close (1)
      close (2)
      close (3)

      return
      end
c
c     ------------------------------------------------------------
      subroutine Write_Particles_Binary1_old2 ( FileName1 , FileName2 )
c     ------------------------------------------------------------
c
c     purpose: writes control information and particles to the specified 
c              files (this routine is to be used when nspec = 1)
c
c     input  : FileName1 - C3CRD*; FileName2 - C3crs0*
c

      include 'a_tree.h'
      include 'a_control.h'

      character*15 FileName1 
      character*16 FileName2

c     for Exemplar, SP: nacces = nbyte; for PC: nacces = nrecl 

      nbyte  = nrecl * 4
      nacces = nbyte

      open (4 , file = FileName1,
     &           form = 'UNFORMATTED' , status = 'UNKNOWN')
 
      open (5 , file = FileName2 , access = 'DIRECT',
     &	         status = 'UNKNOWN', recl = NACCES)

c.... write header and control data

      write (4) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                  ,Ocurv,extras
      close (4)

      ip = 0
      nrows = nrow**3/nrowreal**2
      do irow = 1 , nrows
	do in = 1 , npage
          ip   =ip+1
	  xpar(in) = x(ip) 
	  ypar(in) = y(ip) 
	  zpar(in) = z(ip) 
	   vxx(in)  =VX(ip)
	   vyy(IN)  =VY(ip)
	   vzz(IN)  =VZ(ip)
         EndDo
         CALL WRIROW(IROW,5)
      EndDo

      close (5)

      return
      end
c
c     -------------------------------------------------------------
      subroutine Write_Particles_Binary2 ( FileName1 , FileName2 , 
     &                                     FileName3 , FileName4  )
c     -------------------------------------------------------------
c
c     purpose: writes control information and particles to the specified 
c              files (this routine is to be used when nspec = 3)
c
c     input: FileName1 - C3CRD*;  FileName2 - C3crs0*
c            FileName3 - C3crs1*; FileName4 - C3crs2*
c
      include 'a_tree.h'
      include 'a_control.h'

      character*15 FileName1 
      character*16 FileName2
      character*16 FileName3
      character*16 FileName4

c     for Exemplar, SP: nacces = nbyte; for PC: nacces = nrecl 

      nbyte  = nrecl * 4
      nacces = nbyte / 4

      open (4 , file = FileName1,
     &          form = 'unformatted' , status = 'unknown')
 
      open (5 , file = FileName2 , access = 'direct',
     &	        status = 'unknown', recl = nacces     )
      open (6 , file = FileName3 , access = 'direct',
     &	        status = 'unknown', recl = nacces     )
      open (7 , file = FileName4 , access = 'direct',
     &	        status = 'unknown', recl = nacces     )

      write (4) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
      close (4)
      ip = 0
      nrows = nrow**3/nrowreal**2
      do ic1 = 1 , 3
        do irow = 1 , nrow
	  do in = 1 , npage
            ip       = ip + 1
	    xpar(in) = x(ip) 
	    ypar(in) = y(ip) 
	    zpar(in) = z(ip) 
	     vxx(in) = vx(ip)
	     vyy(in) = vy(ip)
	     vzz(in) = vz(ip)
          enddo
          call WriRow ( irow , 4+ic1 )
        enddo
      enddo

      close (5)
      close (6)
      close (7)

      return
      end
