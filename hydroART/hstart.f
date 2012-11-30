c
c          Set up ICs for gas from the already set files with 
c                      particle coordinates
c  Anatoly Klypin and Andrey Kravtsov (2001)
c
      include 'hstart.h'
C
      COMMON / KINENRG/ SKINE,SX,SY,SZ,SX2,SY2,SZ2
      real*8 ALPHA
      Character  Hd*5,Tail*4
      DATA PI / 3.1415926535 /

      Hd  ='PMcrs'
      Tail='.DAT'
     
      h100 = 0.70 
      Omb0 = 0.022/(h100)**2
      SCALEL = 10.0/h100  ! Boxh/h in Mpc


      call Read_Particles_Binary ()

c      call Write_Particles_ASCII ('PMcrd_a.DAT      ',
c     &                             'PMcrs0_a.DAT      ')

      write (*,*) 'Nspecies =',Nspecies
      write (*,*) 'aexpn =',AEXPN,' astep =',ASTEP
      AEXP0 = AEXPN
      AEXPV = AEXPN - ASTEP/2.

      NspecM   = Nspecies - 1 
      Nmin_lev = 1  ! Nspecies
c      Nspecies = 1
      write (*,*) ' Max resolution Level=',NspecM,
     &               ' Min Level=',Nmin_lev
      write (16,*) ' Max resolution Level=',NspecM,
     &               ' Min Level=',Nmin_lev
      If(Nspecies.ne.0.and.2**(Nspecies-1).ne.Lblock)Then
         write (*,*) ' Wrong number of Levels for mass refinement (',
     &                  Nspecies,') Block size =',Lblock,
     &                  '       2**(Levels-1) must = Block size'
c         Stop
      Endif 
      If(Nspecies.eq.0)Then
         Nparticles =lspecies(1)
      Else
         Nparticles =lspecies(Nspecies)
      EndIf

c      tcorr = 2./sqrt(Om0) * 0.5 ! *0.5 is to compensate for the NGRID change
CEVERINO09272005      tcorr = 1.0 
      tcorr =  2./sqrt(Om0)
      xmin = 1.e6
      xmax = -xmin
      ymin = 1.e6
      ymax = -ymin
      zmin = 1.e6
      zmax = -zmin
      vxmin =  1.e6
      vxmax = -vxmin
      vymin =  vxmin
      vymax = -vymin
      vzmin =  vxmin
      vzmax = -vzmin  
      do ip = 1 , Nparticles 
c        Xpt(ip) = 0.5 * (Xpt(ip) - 1.0) + 1.
c        Ypt(ip) = 0.5 * (Ypt(ip) - 1.0) + 1.
c        Zpt(ip) = 0.5 * (Zpt(ip) - 1.0) + 1.
        xmin = min(xmin,Xpt(ip))
        xmax = max(xmax,Xpt(ip))
        ymin = min(ymin,Ypt(ip))
        ymax = max(ymax,Ypt(ip))
        zmin = min(zmin,Zpt(ip))
        zmax = max(zmax,Zpt(ip))
        Vxt(ip) = Vxt(ip) * tcorr 
        Vyt(ip) = Vyt(ip) * tcorr
        Vzt(ip) = Vzt(ip) * tcorr
        vxmin   = min(vxmin,vxt(ip) )
        vxmax   = max(vxmax,vxt(ip))
        vymin   = min(vymin,vyt(ip))
        vymax   = max(vymax,vyt(ip))
        vzmin   = min(vzmin,vzt(ip))
        vzmax   = max(vzmax,vzt(ip))
c        pw(ip) = pw(ip)* 0.125 ! compensate for NGRID change        
      enddo
      write(*,*) 'rescaled particle coordinate range: tcorr= ',tcorr 
      write(*,*) 'xmin,max =',xmin,xmax
      write(*,*) 'ymin,max =',ymin,ymax
      write(*,*) 'zmin,max =',zmin,zmax
      write(*,*) 'vxmin,max   =',vxmin,vxmax
      write(*,*) 'vymin,max   =',vymin,vymax
      write(*,*) 'vzmin,max   =',vzmin,vzmax

c
      fact  =  sqrt(Om0+Oml0*AEXPV**3+Ocurv*AEXPV)
      vfact =  sqrt(Om0+Oml0*AEXPN**3+Ocurv*AEXPN)
c
c.... correct for the difference in t0 units in N-body and Hydro
c     and for the difference in time position of particle and gas vels. 
c
      VCORR =  AEXPN*SQRT(AEXPN)/(AEXPV*SQRT(AEXPV))*vfact / fact

c
c.... Now set up and write hydro variables on the (ngrid,ngrid,ngrid) grid
c

      write(*,*) 'setting hydro variables on the grid...'
      write(*,*) 'setting hydro variables...'
      call SetHydro2 ( VCORR )
      write(*,*) 'writing hydro data...'
      call WriteHydro (SCALEL*h100)
c
      write(*,*) 'writing particles...'      
      do i = 1 , Nspecies
        wspecies(i) = (1. - Omb0/Om0) * wspecies(i)
      enddo

      call Write_Particles_Binary ('PMcrd.DAT      ', 
     &                             'PMcrs0.DAT      ')

c      call Write_Particles_ASCII ('PMcrd_a.DAT      ',
c     &                             'PMcrs0_a.DAT      ')

      Do ic1 =1,Nparticles
         XPt(ic1) =astep
      Enddo 
      open ( 60 , file = 'pt.dat' , form = 'unformatted' )
      write(60) (XPt(ic1),ic1=1,Nparticles)
      close(60)


      STOP
      END
c
c
c     -----------------------------------
      subroutine Read_Particles_Binary ()
c     -----------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
      include 'hstart.h'
c
c      real    wspecies(nspec)
c      integer lspecies(nspec)
c      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
c 
c      ngrid = ng

      open ( 2,file = 'Result.DAT')
      open ( 3 ,file ='PMcrd_i.DAT', form = 'unformatted')

c.... read control information and check whether it has proper structure

      read      (3) HEADER, 
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                   ,Ocurv,extras
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
     &            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I8,/
     &            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &               ' Omegab_0=',F7.3,' Hubble=',f7.3)
      if ( NGRIDC .ne. NGRID ) then 
        write(*,*) 'NGRIDC .ne. NGRID :',NGRIDC, NGRID
        write(*,*) 'hope this is ok...'
      endif
      write(*,*) 'Nspecies as read from particle file=',nspecies
      if(nspecies .eq. 0 ) nspecies = 1
      If( Nmaxpart .lt. lspecies(nspecies) ) then
      write (*,*) ' Wrong number of particles: '
      write (*,*) ' must be at least=',lspecies(nspecies),' (lspecies)'
      write (*,*) ' but is set to ',Nmaxpart,' in a_setup.h...'
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo
         STOP
      Endif 
      write(*,*) 'numbers of particles of different species:'
      do ispec = 1 , nspecies
        write(*,*) ispec, lspecies(ispec)
      enddo
      
      nbyte  = nrecl * 4
      nacces = nbyte / nbyteword
 
      open ( 1 , file = 'PMcrs0_i.DAT', access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3

      N_particles = lspecies(nspecies)   ! Total number of particles
      Npages      = (N_particles -1)/npage + 1
      N_in_last   = N_particles - npage*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) 'N_particles =',N_particles, ' N_in_last=',N_in_last
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo

      do irow = 1 , Npages         ! loop over particle pages
        In_page = npage
        if ( irow .eq. Npages ) In_page = N_in_last
         write (*,*)' Read page=',IROW,' file=',ifile,' N=',In_page
        iL = npage * (irow-1)
        CALL GetRow(irow,1) ! read in a page of particles
        do in = 1 , In_page          ! Loop over particles
          ip = in + iL                     ! current particle number
          XPt(ip) = xpar(in)
          YPt(ip) = ypar(in)
          ZPt(ip) = zpar(in)
          VXt(ip) = vx(in)
          VYt(ip) = vy(in)
          VZt(ip) = vz(in)
          if ( XPt(ip) .gt. 257. .or.
     &         YPt(ip) .gt. 257. .or.
     &         ZPt(ip) .gt. 257. ) then 
          write(*,*) ip, XPt(ip), YPt(ip), ZPt(ip)
          write(*,*) In_page, iL
          pause
          endif
        enddo
      enddo

      close (1)
      close (2)
      close (3)

      do i = 1 , NSpecies
        if(i .eq. 1 ) then
          nb1 = 1
        else
          nb1 = lspecies(i-1) + 1
        endif 
        nb2 = lspecies(i)
        do ip = nb1 , nb2
          pw(ip) = wspecies(i)
        enddo
      enddo 

      return
      end

c     ------------------------------------------------------------
      subroutine Write_Particles_Binary ( FileName1 , FileName2 )
c     ------------------------------------------------------------
c
c
c     input  : FileName1 - C3CRD*; FileName2 - C3crs0*
c

      include 'hstart.h'

      character*15 FileName1
      character*16 FileName2

      nbyte  = nrecl * 4
      nacces = nbyte / nbyteword

      open (4 , file = FileName1,
     &           form = 'UNFORMATTED' , status = 'UNKNOWN')

      open (5 , file = FileName2 , access = 'DIRECT',
     &           status = 'UNKNOWN', recl = NACCES)

c.... write header and control data

      write (4) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRID,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                  ,Ocurv,Omb0,extras
      close (4)

      N_particles = lspecies(nspecies)   ! Total number of particles
      Npages = (N_particles -1)/npage +1
      N_in_last = N_particles - npage*(Npages-1)
      write (*,*) 'N_particles =',N_particles
      write(*,*) ' Pages=',Npages,' Species=',nspecies
      write (*,*) ' N_in_last=',N_in_last

      do irow = 1 , Npages         ! Loop over particle pages
        In_page = npage
        If ( irow .eq. Npages ) In_page = N_in_last
c         write (*,*)' Write page=',IROW,' file=',ifile,' N=',In_page
        iL = npage * (irow-1)
        do in = 1 , In_page          ! Loop over particles
          ip = in + iL                     ! current particle number
          xpar(in) = Xpt(ip)
          ypar(in) = Ypt(ip)
          zpar(in) = Zpt(ip)
           vx(in) = Vxt(ip)
           vy(in) = Vyt(ip)
           vz(in) = Vzt(ip)
         enddo
         call WriRow ( irow , 5 )
      enddo

      close (5)

      return
      end
c
c     ------------------------------------------------------------
      subroutine Write_Particles_ASCII ( FileName1 , FileName2 )
c     ------------------------------------------------------------
c
c     purpose: writes control information and particles to the specified 
c              files (this routine is to be used when nspec = 1)
c
c     input  : FileName1 - C3CRD*; FileName2 - C3crs0*
c

      include 'hstart.h'

      character*15 FileName1 
      character*16 FileName2

      nbyte  = nrecl * 4
      nacces = nbyte / nbyteword

      open (4 , file = FileName1,
     &          status = 'UNKNOWN')
 
      open (5 , file = FileName2 , 
     &	        status = 'UNKNOWN')

c.... write header and control data

      write (4,*)
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                  ,Ocurv,Omb0,extras
      close (4)

      N_particles = lspecies(nspecies)   ! Total number of particles
      Npages = (N_particles -1)/npage +1
      N_in_last = N_particles - npage*(Npages-1)
      write (*,*) 'N_particles =',N_particles
      write(*,*) ' Pages=',Npages,' Species=',nspecies
      write (*,*) ' N_in_last=',N_in_last

      if ( nspecies .eq. 1 ) then 
        np0 = N_particles
      else
        np0 = lspecies(2)
      endif

      frac = 0.1 
      m = 38233663
      do i = 1 , np0 
        fd = RANDd ( m )
        if ( fd .lt. frac ) then 
          write(5,190) Xpt(i),Ypt(i),Zpt(i),Vxt(i),Vyt(i),Vzt(i)
 190      format(3(f9.5,2x),1x,3(e10.4,2x))
        endif
      enddo
      close (5)

      return
      end

c     ----------------------------------
      subroutine WriRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'hstart.h'
      integer irow , ifile 
      write (ifile , rec = irow ) recdat
      return
      end
 

c     ------------------------
      subroutine WriteHydro ( boxh )
c     ------------------------
c
      include 'hstart.h'
      real abuff(ngrid**3)
      common / TMPGRID / abuff
c
      open ( 50 , file = 'tr_ic.dat', form ='unformatted' )
      write(50) boxh
      write(50) aexpn, astep
      write(50) ngrid**3
c
c.... gas density
c
      ic = 0 
      do i = 1 , ngrid 
        do j = 1 , ngrid 
          do k = 1 , ngrid 
            ic = ic + 1
            abuff(ic) = rho0(i,j,k)
          enddo
        enddo
      enddo
      write(50) (abuff(ic),ic=1,ngrid**3)
c
c.... components of the gas momentum
c
      ic = 0 
      do i = 1 , ngrid 
        do j = 1 , ngrid 
          do k = 1 , ngrid 
            ic = ic + 1
            abuff(ic) = vx0(i,j,k)
          enddo
        enddo
      enddo
      write(50) (abuff(ic),ic=1,ngrid**3)
      ic = 0 
      do i = 1 , ngrid 
        do j = 1 , ngrid 
          do k = 1 , ngrid 
            ic = ic + 1
            abuff(ic) = vy0(i,j,k)
          enddo
        enddo
      enddo
      write(50) (abuff(ic),ic=1,ngrid**3)
      ic = 0 
      do i = 1 , ngrid 
        do j = 1 , ngrid 
          do k = 1 , ngrid 
            ic = ic + 1
            abuff(ic) = vz0(i,j,k)
          enddo
        enddo
      enddo
      write(50) (abuff(ic),ic=1,ngrid**3)
c
c...  gas energy
c
      ic = 0 
      do i = 1 , ngrid 
        do j = 1 , ngrid 
          do k = 1 , ngrid 
            ic = ic + 1
            abuff(ic) = eng0(i,j,k)
          enddo
        enddo
      enddo
      write(50) (abuff(ic),ic=1,ngrid**3)
      ic = 0
      do i = 1 , ngrid
        do j = 1 , ngrid
          do k = 1 , ngrid
            ic = ic + 1
            abuff(ic) = ei0(i,j,k)
          enddo
        enddo
      enddo
      write(50) (abuff(ic),ic=1,ngrid**3)

c
      close ( 50 )
c
      return
      end
c
c     -----------------------------
      subroutine SetHydro2 ( vcorr )
c     -----------------------------
      include 'hstart.h'

      do k = 1 , ngrid
        do j = 1 , ngrid 
          do i = 1 , ngrid 
            rho0(i,j,k) = 0.
            eng0(i,j,k) = 0.
             vx0(i,j,k) = 0.
             vy0(i,j,k) = 0.
             vz0(i,j,k) = 0.
             ei0(i,j,k) = 0.
          enddo
        enddo
      enddo
           
        do ippp = 1 , lspecies(NSpecies)
          x  = XPt(ippp)
          y  = YPt(ippp)
          z  = ZPt(ippp)
          v1 = VXt(ippp) * vcorr
          v2 = VYt(ippp) * vcorr
          v3 = VZt(ippp) * vcorr
          i  = int(x)
          j  = int(y)
          k  = int(z)

          s1  = x - float(i)
          s2  = y - float(j)
          s3  = z - float(k)
          t1  = 1.0 - s1
          t2  = 1.0 - s2
          t3  = 1.0 - s3
          t2w = t2 * pw(ippp)
          s2w = s2 * pw(ippp)

          i1 = i + 1
          if ( i1 .gt. ngrid ) i1 = 1
          j1 = j + 1
          if ( j1 .gt. ngrid ) j1 = 1
          k1 = k + 1
          if ( k1 .gt. ngrid ) k1 = 1

          t3t1t2 = t3*t1*t2w
          t3s1t2 = t3*s1*t2w
          t3t1s2 = t3*t1*s2w
          t3s1s2 = t3*s1*s2w         
          s3t1t2 = s3*t1*t2w
          s3s1t2 = s3*s1*t2w
          s3t1s2 = s3*t1*s2w
          s3s1s2 = s3*s1*s2w
         
          rho0(i ,j ,k ) = rho0(i ,j ,k ) + t3t1t2
          rho0(i1,j ,k ) = rho0(i1,j ,k ) + t3s1t2
          rho0(i ,j1,k ) = rho0(i ,j1,k ) + t3t1s2
          rho0(i1,j1,k ) = rho0(i1,j1,k ) + t3s1s2
          rho0(i ,j ,k1) = rho0(i ,j ,k1) + s3t1t2
          rho0(i1,j ,k1) = rho0(i1,j ,k1) + s3s1t2
          rho0(i ,j1,k1) = rho0(i ,j1,k1) + s3t1s2
          rho0(i1,j1,k1) = rho0(i1,j1,k1) + s3s1s2
         
          vx0(i ,j ,k ) = vx0(i ,j ,k ) + t3t1t2*v1
          vx0(i1,j ,k ) = vx0(i1,j ,k ) + t3s1t2*v1
          vx0(i ,j1,k ) = vx0(i ,j1,k ) + t3t1s2*v1
          vx0(i1,j1,k ) = vx0(i1,j1,k ) + t3s1s2*v1
          vx0(i ,j ,k1) = vx0(i ,j ,k1) + s3t1t2*v1
          vx0(i1,j ,k1) = vx0(i1,j ,k1) + s3s1t2*v1
          vx0(i ,j1,k1) = vx0(i ,j1,k1) + s3t1s2*v1
          vx0(i1,j1,k1) = vx0(i1,j1,k1) + s3s1s2*v1
         
          vy0(i ,j ,k ) = vy0(i ,j ,k ) + t3t1t2*v2
          vy0(i1,j ,k ) = vy0(i1,j ,k ) + t3s1t2*v2
          vy0(i ,j1,k ) = vy0(i ,j1,k ) + t3t1s2*v2
          vy0(i1,j1,k ) = vy0(i1,j1,k ) + t3s1s2*v2
          vy0(i ,j ,k1) = vy0(i ,j ,k1) + s3t1t2*v2
          vy0(i1,j ,k1) = vy0(i1,j ,k1) + s3s1t2*v2
          vy0(i ,j1,k1) = vy0(i ,j1,k1) + s3t1s2*v2
          vy0(i1,j1,k1) = vy0(i1,j1,k1) + s3s1s2*v2
         
          vz0(i ,j ,k ) = vz0(i ,j ,k ) + t3t1t2*v3
          vz0(i1,j ,k ) = vz0(i1,j ,k ) + t3s1t2*v3
          vz0(i ,j1,k ) = vz0(i ,j1,k ) + t3t1s2*v3
          vz0(i1,j1,k ) = vz0(i1,j1,k ) + t3s1s2*v3
          vz0(i ,j ,k1) = vz0(i ,j ,k1) + s3t1t2*v3
          vz0(i1,j ,k1) = vz0(i1,j ,k1) + s3s1t2*v3
          vz0(i ,j1,k1) = vz0(i ,j1,k1) + s3t1s2*v3
          vz0(i1,j1,k1) = vz0(i1,j1,k1) + s3s1s2*v3
         
        enddo

c      call FILTER

      sum    = 0.
      sum2   = 0.
      disp   = 0.
      rhomax = -1.e6
      rhomin = -rhomax
      x0     = SCALEL*hubble/ngrid ! cell size in units of 100 kpc
      amu    = 0.6 ! mean molecular weight of the primordial gas
      gamma  = 5./3.      
      write(*,*) 'SCALEL =',SCALEL, ngrid
c   
c...  initial temperature in code units
c
      a_th = 1.0 / (1.e3*(Omb0*hubble**2)**0.4)
         
      write(*,*) '3: a_th =', a_th, Omb0, hubble
      if ( aexpn .lt. a_th ) then 
        TinitK = 2.726 / aexpn
      else
        TinitK = 2.726 / a_th * (a_th/aexpn)**2
      endif

      Tinit  = aexpn**2/amu/Om0/x0**2 * (TinitK/3.03e5)/(gamma - 1.)      
      write(*,*) 'Tinit =',Tinit, amu, x0, gamma

      rho_min = 1.e9
      rho_max = -rho_min
      vxmin = rho_min
      vxmax = -vxmin
      vymin = rho_min
      vymax = -vymin
      vzmin = rho_min
      vzmax = -vzmin
      e_min = rho_min
      e_max = -e_min

      f_b = Omb0/Om0
      do k = 1 , ngrid
        do j = 1 , ngrid
          do i = 1 , ngrid
            sum2 = sum2 + rho0(i,j,k)
            disp = disp + rho0(i,j,k)**2
            rhomin = min(rho0(i,j,k),rhomin)
            rhomax = max(rho0(i,j,k),rhomax)
            rho0(i,j,k) = max(rho0(i,j,k)*f_b,0.)
            vx0(i,j,k)  = vx0(i,j,k) * f_b               
            vy0(i,j,k)  = vy0(i,j,k) * f_b
            vz0(i,j,k)  = vz0(i,j,k) * f_b      
            eps         = Tinit * rho0(i,j,k) 
            Ekin        = (vx0(i,j,k)**2+vy0(i,j,k)**2+vz0(i,j,k)**2)/
     &                    rho0(i,j,k) / 2.
            eng0(i,j,k) = Ekin + eps
            ei0(i,j,k)  = eps
            sum  = sum  + rho0(i,j,k)
            rho_min = min(rho_min,rho0(i,j,k))
            rho_max = max(rho_max,rho0(i,j,k))
            e_min   = min(e_min,eps)
            e_max   = max(e_max,eps)
            vxmin   = min(vxmin,vx0(i,j,k))
            vxmax   = max(vxmax,vx0(i,j,k))
            vymin   = min(vymin,vy0(i,j,k))
            vymax   = max(vymax,vy0(i,j,k))
            vzmin   = min(vzmin,vz0(i,j,k))
            vzmax   = max(vzmax,vz0(i,j,k))
          enddo
        enddo
      enddo

      write(*,*) 'rho_min,max =',rho_min,rho_max
      write(*,*) 'ei_min,max  =',e_min,e_max
      write(*,*) 'vxmin,max   =',vxmin,vxmax
      write(*,*) 'vymin,max   =',vymin,vymax
      write(*,*) 'vzmin,max   =',vzmin,vzmax

      cells = ngrid * ngrid * ngrid
      disp  = sqrt(disp/cells-(sum2/cells)**2)

      write (*,*) ' Mass of gas =',sum
      write (*,*) ' Total mass =',sum2,' Ncells=',cells,' Sigma=',disp
      write (*,*) ' rhomax =',rhomax,' rhomin =',rhomin

      return
      end

c
C--------------------------------------------------
C     sqrt(Power spectrum)
C     k = (2pi/L) wk
      REAL  FUNCTION TRUNF(wk)
C--------------------------------------------------
      include 'hstart.h'
      real k
      IF (WK.GE.FLOAT(NSPEC)) THEN
	       TRUNF =0.
	return
      ENDIF
      If(Par(6).ne.0.)Then   ! Holtzman approx
        k = QSCALE*wk
        sk= sqrt(k)
        TRUNF = k**(ns/2.) /
     .         (1.+sk*(Par(2)
     .            +sk*(Par(3)
     .            +sk*(Par(4)
     .            +sk*(Par(5) )))) )**(Par(6))
      Else                   ! BBKS + Sugiyama approx
c        Gammaeff =Om*hsmall/exp(Omb*(1.+sqrt(hsmall/0.5)/Om))
c        Q = wk/hsmall/Gammaeff
        k = QSCALE*wk
        Q = k*qqscaleb
        TRUNF = k**(ns/2.)* LOG(1.+Par(1)*Q)/(Par(1)*Q)/
     .          sqrt(sqrt(1.+Q*(Par(2)+
     .                  Q*(Par(3)**2+
     .                  Q*(Par(4)**3+
     .                  Q*(Par(5)**4) )))))
      EndIf
      RETURN
      END

C*********************************************************************
C			  INITIALIZE CONSTANTS:
C			      Scalel    =length of square in MPC
C			      AEXPN =expansion factor (current)
C			      ASTEP	  =increment of expansion factor: AEXPN =AEXP_0+N*ASTEP
C                     PARTW	=weight of a particle: total 'mass' per cell must be unity
c                                          even for LCDM or OCDM 
C			      AMPLT	=amplitude of perturbations inside the simulation box
C			      NBYTE	=record length in bytes
C
      SUBROUTINE InitValues(NBYTE)
      include 'hstart.h'
      COMMON / FERMI  / Vscale
      Character            Answer*1
      Real                 INPUT
      DATA PI		         /3.1415926535/

       Write (*,*) 'Would you like to use a model provided in cdm.fit '
       Write (*,*) 'or your own model?' 
       Write (*,'(A,$)') ' Enter Yes for cdm.fit;  No for your model ='
       Read  (*,*)  Answer
       Answer = 'Y'
       HEADER='N=128x256 L=20h-1CDMsig=0.239 z=30-----------'
       write (*,*) '------ Enter Header for the run up to 45 characters'
       write (*,*) '       Example:'
       write (*,'(A)') HEADER                                   
       read  (*,'(A)') HEADER
       write (*,'(A)') HEADER
       AEXPN =INPUT(' Initial expansion parameter (0-1)=')
       ASTEP =INPUT(' Step of the expansion parameter  =')
       AMPLT =INPUT(' Amplitude of density fluctuations=')
       ns    =INPUT(' Slope of the Power spectrum     n=')
       SCALEL=INPUT(' Box size (Mpc/h)                 =')
       Ocurv =INPUT(' Omega_curvature   at present     =')
       Nseed =INPUT(' Random seed (Integer  1-2^31)    =') 
       If(Answer.eq.'Y' .or. Answer.eq.'y')Then
          CALL MODEL(hsmall)
          hubble=hsmall
          Om0   =Om
          Oml0  =1. -Om0 -Ocurv
       Else
       write (*,*) ' You use your cosmological model.'
       write (*,*) ' Be sure you provide routine TRUNF'
       hubble=INPUT(' The Hubble constant (h=H/100)    =')
       Om0   =INPUT(' Omega_0 matter    at present     =')
       Oml0  =INPUT(' Omega_lambda      at present     =')
       EndIf
       SCALEL=SCALEL/hubble   ! scale it to real megaparsecs      
       NBYTE = NPAGE*6*4
       Nspecies =INPUT(' Number of mass species   = ')
             W     = (FLOAT(NGRID)/FLOAT(NROW))**3
            PARTW = W 
      If(Nspecies.eq.0)Then  ! old constant mass resolution
              write(*,*) ' You use PM code with one particle set '
       Else    ! multiple mass resolution
              write(*,*) 
              write(*,*) ' You use multiple mass resolution'
       endif 
       ISTEP = 0
       TINTG = 0.
       AU0   = 0.
       AEU0  = 0.
       EKIN  = 0.
       EKIN1 = 0.
       EKIN2 = 0.
       NROWC = NROW
       NGRIDC= NGRID
	     QSCALE = 2.*PI/SCALEL
	     QS     = hubble**(-2)
        write (*,*) ' Qscale=',QSCALE,SCALEL,ns
        write(*,*) 'AMPLT =',AMPLT

      RETURN
      END
C________________________________________Read parameters of the model
C                                 Om       = Omega_0
C                                 Omb     = Omega_baryon
C                                 Omnu  = Omega_neutrino
C                                 hsmall  = hubble constant
C                                 Par      = fitting parameters
C                                 Par(6) = 0  --- bbks-style+Hu&Sugiyama
C                                 Par(6) ne 0 --- holtzman-style
      SUBROUTINE MODEL(hsmall)
C---------------------------------------
      include 'hstart.h'
      Real     INPUT
      Character             Header1*79

      OPEN(2,file='cdm.fit')
c      Line1 =INPUT(' Enter Line Number in cdm.fit     =')
c      Read(2,'(A)') Header1
c      If(Line1.gt.2)Then
c         Do i=1,Line1-2
c            Read (2,*) a
c         EndDo
c      EndIf
      Read (2,*) Om,Omb,Omc,Gamma,hsmall,Par
      CLOSE (2) 
      
      theta = 2.726/2.7  ! = T_cmb/2.7
      Ob0   =Omb/Om      ! Hu&Sugiyama fitting
      Omh2  =Om*hsmall**2
      a1    =(46.9*Omh2)**0.670*(1.+(32.1*Omh2)**(-0.532))
      a2    =(12.0*Omh2)**0.424*(1.+(45.*Omh2)**(-0.582))
      alpha =1./a1**Ob0/a2**(Ob0**3)
      qqscaleb = theta**2/(1.-Ob0)**0.60/sqrt(alpha)/Omh2
      
      qqscaleb = 1.d0 / (Gamma * hsmall )

      Write (*,20)Om,Omb,Omc,Gamma,hsmall,Par
c      Write (16,20)Om,Omb,Omc,Omnu,hsmall,Par
 20   Format(' Model: Om_0=',F5.3,' Om_baryon=',F5.3,
     .       ' Om_cold=',F5.3,' Gamma=',F5.3,
     .       ' hsmall=',F4.2,/8x,'Parameters=',(6G9.3))
      Return
      End
C------------------------------------------------
C---------------------------------------------------
C---------------------------------------------------
C                                  Read  current data from disk/tape,
C                                  Open files
C                                  Nrecl is the number of values in a record
C                                  Npage is the number of particles in a page
      SUBROUTINE RDTAPE
C---------------------------------------------------
      include 'hstart.h'
      Character  Hd*5,Tail*4,Nm*40
      Hd  ='PMcrs'
      Tail='.DAT'
C                                     Open file on a tape
      OPEN(UNIT=9,FILE='PMcrd.DAT',
     +                FORM='UNFORMATTED', STATUS='UNKNOWN')
C                                 Read control information
C                                 and see whether it has proper
C                                 structure
      READ  (9,err=10,end=10) HEADER,
     +                  AEXPN,AEXP0,AMPLT0,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                   ,Ocurv,extras
      Ocurv =0.
      WRITE (*,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT0,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv
      WRITE (16,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT0,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv
100   FORMAT(1X,'Header=>',A45,/
     +            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3)
      IF(NROW.NE.NROWC) THEN
         WRITE (*,*)
     +            ' NROW in PARAMETER and in TAPE-FILE are different'
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (*,*)
     +           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      ENDIF
C                                         Open work file on disk
c     Intel PC/Linux: NACCES = NBYTE
 10   NBYTE = NRECL*4
      NACCES= NBYTE

         Nm =Hd//Char(48)//Tail
         iun=21
         write (*,*) ' File>',ifile,' unit=',iun,' Name=',Nm
         OPEN(UNIT=iun,FILE=Nm,ACCESS='DIRECT',
     +	               STATUS='UNKNOWN',RECL=NACCES)

 
      REWIND 9
      RETURN
      END
C---------------------------------------------
C                       Write current data to  disk/tape
C
      SUBROUTINE WRTAPE
C----------------------------------------------
      include 'hstart.h'
C                                       write header and control data
      WRITE (9) HEADER,
     +           AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +           TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +           NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     +           Ocurv,extras
      REWIND 9
      RETURN
      END
cC--------------------------------------------------
cC                             Set Weights of particles for
cC                             fast access
c      SUBROUTINE SetWeight
cC--------------------------------------------------
c      INCLUDE 'PMparameters.h'
c      If(Nspecies.eq.0)Then  ! old  constant weights
c         N_particles =lspecies(1)
c         Do i=1,N_particles
c            iWeight(i) =PARTW
c         EndDo
c      Else
c         N_particles =lspecies(Nspecies)
c         jstart =1
c         Do j=1,Nspecies
c            jend =lspecies(j)
c            Do k=jstart ,jend
c               iWeight(k) =wspecies(j)
c            EndDo
c            jstart =jend
c         EndDo
c      EndIf
c      write (*,*) ' Set Weights for ',N_particles,' particles'
c      RETURN
c      END
c     ----------------------------------
      subroutine GetRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'hstart.h'
      read  ( ifile , rec = irow ) recdat
      return
      end

C------------------------------------------------
C				                                       random number generator
      FUNCTION RANDd(M)
C------------------------------------------------
      DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4
     +	/453815927,2147483648.,2147483647,536870912,131072,256,
     +	 16777216,4,16384,8388608,128/
      ML=M/K1*K1
      M1=(M-ML)*L1
      ML=M/K2*K2
      M2=(M-ML)*L2
      ML=M/K3*K3
      M3=(M-ML)*L3
      ML=M/K4*K4
      M4=(M-ML)*L4
      M5=KI-M
      IF(M1.GE.M5)M1=M1-KI-1
      ML=M+M1
      M5=KI-ML
      IF(M2.GE.M5)M2=M2-KI-1
      ML=ML+M2
      M5=KI-ML
      IF(M3.GE.M5)M3=M3-KI-1
      ML=ML+M3
      M5=KI-ML
      IF(M4.GE.M5)M4=M4-KI-1
      ML=ML+M4
      M5=KI-ML
      IF(LC.GE.M5)ML=ML-KI-1
      M=ML+LC
      RANDd=M/AM
      RETURN
      END
C--------------------------------------
C				normal random numbers
      FUNCTION GAUSS(M)
C--------------------------------------
      X=0.
      DO  I=1,5
         X=X+RANDd(M)
      EndDo
      X2   =1.5491933*(X-2.5)
      GAUSS=X2*(1.-0.01*(3.-X2*X2))
      RETURN
      END

C---------------------------------- Read in variables      
      REAL FUNCTION INPUT(text)
C------------------------------------------------
      Character text*(*)
          write (*,'(A,$)')text
          read (*,*) x
          INPUT =x
      Return
      End
C--------------------------------------   Fourier Transform
      SUBROUTINE SETF67(JBC1,JQ1,jbc,jp,jsl,ll1,k2,k3,k4,k7,
     &                  Qi,jndx,Uf,Vf)
c     -----------------------------------------------------
      include 'hstart.h'

      dimension Uf(marr) , Vf(marr)
      dimension QI(mf67) , jNDX(mf67)

      PI=DATAN(1.D0)*4.
       JBC=JBC1
       JQ=JQ1
       IF(JBC.LT.3) GO TO 101
       JQ=JQ-1
  101      K3=2**JQ
       K7=K3/2
       N5=K3/4
       I=1
       JNDX(I)=N5
       QI(I)=0.5*SQRT(2.)
       K=I
       I=I+1
  102     IL=I
       IF(I.EQ.K7) GO TO 104
  103     K1=JNDX(K)/2
       JNDX(I)=K1
       QI(I)=SIN(PI*FLOAT(K1)/FLOAT(K3))
       K1=K7-K1
       I=I+1
       JNDX(I)=K1
       QI(I)=SIN(PI*FLOAT(K1)/FLOAT(K3))
       K=K+1
       I=I+1
       IF(K.EQ.IL) GO TO 102
       GO TO 103
  104     RETURN
       END
C-----------------------------------------
       SUBROUTINE TFOLD(IS,L,ZZZ,k2)
       include 'hstart.h'


       DIMENSION ZZZ(MARR)
       IH2=K2/2-1
       DO 100 I=IS,IH2
       I1=I+L
       I2=K2-I+L
       A=ZZZ(I1)
       ZZZ(I1)=A-ZZZ(I2)
       ZZZ(I2)=A+ZZZ(I2)
  100     CONTINUE
       RETURN
       END
C------------------------------------
       SUBROUTINE NEG(I1,I3,I2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       include 'hstart.h'
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)
 
       DO 100 K=I1,I3,I2
       Vf(K)=-Vf(K)
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE REVNEG(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       include 'hstart.h'
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       DO 100 I=1,K7
       J=K3+1+I
       K=K4+1-I
       A=Vf(J)
       Vf(J)=-Vf(K)
       Vf(K)=-A
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE ZEERO(L,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       include 'hstart.h'
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       DO 100 I=1,K2
       LI=L+I
       Uf(LI-1)=0.0
  100     CONTINUE
       RETURN
       END
C------------------------------------------
       SUBROUTINE TFOLD1(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
      include 'hstart.h'
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       II=K2-1
       DO 100 I=1,II
       I1=JSL+I
       I2=LL1-I
       A=Uf(I1)
       Uf(I1)=A+Uf(I2)
       Uf(I2)=A-Uf(I2)
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       include 'hstart.h'
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       JS1=K2
       I=1
       J5=JSL+K2
       IS1=JSL
       IC1=LL1
       JS1=JS1/2
       IF(JS1.NE.1) GO TO 200
       K1=JNDX(I)
       SN=QI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*(Uf(IC1)-Uf(IS1))
       Vf(K1+1)=Uf(IC0)+ODD1
       K3MK1=K3-K1
       Vf(K3MK1+1)=Uf(IC0)-ODD1
       IF(JBC.LT.3) GO TO 110
       ODD2=SN*(Uf(IC1)+Uf(IS1))
       K3PK1=K3+K1
       Vf(K3PK1+1)=Uf(IS0)+ODD2
       K4MK1=K4-K1
       Vf(K4MK1+1)=-Uf(IS0)+ODD2
  110     RETURN
  200     SN=QI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  210     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*(Uf(IC1)-Uf(IS1))
       ODD2=SN*(Uf(IC1)+Uf(IS1))
       Uf(IC1)=Uf(IC0)-ODD1
       Uf(IS1)=-Uf(IS0)+ODD2
       Uf(IC0)=Uf(IC0)+ODD1
       Uf(IS0)=Uf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 210
       I=I+1
  300     IS1=JSL
       IC1=LL1
       JS1=JS1/2
       IF(JS1.EQ.1) GO TO 400
  310     SN=QI(I)
       I=I+1
       CS=QI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  320     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=CS*Uf(IC1)-SN*Uf(IS1)
       ODD2=SN*Uf(IC1)+CS*Uf(IS1)
       Uf(IC1)=Uf(IC0)-ODD1
       Uf(IC0)=Uf(IC0)+ODD1
       Uf(IS1)=-Uf(IS0)+ODD2
       Uf(IS0)=Uf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 320
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  330     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*Uf(IC1)-CS*Uf(IS1)
       ODD2=CS*Uf(IC1)+SN*Uf(IS1)
       Uf(IC1)=Uf(IC0)-ODD1
       Uf(IC0)=Uf(IC0)+ODD1
       Uf(IS1)=-Uf(IS0)+ODD2
       Uf(IS0)=Uf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 330
       I=I+1
       IF(IS1.EQ.J5) GO TO 300
       GO TO 310
  400     K1=JNDX(I)
       SN=QI(I)
       I=I+1
       CS=QI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=CS*Uf(IC1)-SN*Uf(IS1)
       Vf(K1+1)=Uf(IC0)+ODD1
       K3MK1=K3-K1
       Vf(K3MK1+1)=Uf(IC0)-ODD1
       IF(JBC.LT.3) GO TO 410
       ODD2=SN*Uf(IC1)+CS*Uf(IS1)
       K3PK1=K3+K1
       Vf(K3PK1+1)=Uf(IS0)+ODD2
       K4MK1=K4-K1
       Vf(K4MK1+1)=-Uf(IS0)+ODD2
  410     IS1=IS1+1
       IC1=IC1+1
       K1=JNDX(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*Uf(IC1)-CS*Uf(IS1)
       Vf(K1+1)=Uf(IC0)+ODD1
       K3MK1=K3-K1
       Vf(K3MK1+1)=Uf(IC0)-ODD1
       IF(JBC.LT.3) GO TO 420
       ODD2=CS*Uf(IC1)+SN*Uf(IS1)
       K3PK1=K3+K1
       Vf(K3PK1+1)=Uf(IS0)+ODD2
       K4MK1=K4-K1
       Vf(K4MK1+1)=-Uf(IS0)+ODD2
  420     IS1=IS1+1
       IC1=IC1+1
       I=I+1
       IF(IS1.NE.J5) GO TO 400
       RETURN
       END

c      ------------------------------------------------------
       SUBROUTINE FOUR67(JBC1,JQ1,jp1,jsl1,ll11,k21,k71,
     &                   Qi,jndx,Uf,Vf)
c      ------------------------------------------------------
       include 'hstart.h'
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       JBC=JBC1
       JQ=JQ1
       jp = jp1
       jsl = jsl1
       ll1 = ll11
       k2 = k21
       k7 = k71
       A5=0.5*SQRT(2.0)
       K4=2**JQ
       K3=K4
       GO TO (103,103,101,102),JBC
  101     Uf(1)=Uf(1)/2.0
       Uf(K3+1)=Uf(1)
       K2=K3

       CALL TFOLD(0,1,Uf,k2)

       K3=K3/2
       JQ=JQ-1
       GO TO 103
  102     K3=K3/2
       JQ=JQ-1
  103     N5=K3/4
       K7=K3/2
       N11=3*K7
       K31=K3+1
       GO TO(300,400,500,600),JBC
  300     Uf(K31)=0.0
       Uf(1)=0.0
       K2=K3
       DO 301 I=2,JQ

       CALL TFOLD(1,1,Uf,k2)

  301     K2=K2/2
       Vf(K7+1)=Uf(2)
       JF=N5
       JSL=1
       DO 302 JP=2,JQ
       LL1=K2+1

       CALL ZEERO(1,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       I1=3*JF+1
       I2=4*JF
       I3=I1+(K2/2-1)*I2

       CALL NEG(I1,I3,I2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       K2=K2+K2
  302     JF=JF/2
       RETURN
  400     Uf(1)=Uf(1)/2.0
       Uf(K31)=Uf(K31)/2.0
       K2=K3
       DO 401 I=2,JQ

       CALL TFOLD(0,K31-K2,Uf,k2)

  401     K2=K2/2
       LL1=K31-K2
       A=Uf(LL1)+Uf(LL1+2)
       Vf(1)=-A-Uf(LL1+1)
       Vf(K31)=-A+Uf(LL1+1)
       Vf(K7+1)=Uf(LL1)-Uf(LL1+2)
       DO 402 JP=2,JQ
       JSL=K31-K2
       LL1=JSL-K2
       idum = jsl

       CALL ZEERO(idum,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

  402     K2=K2+K2

       CALL NEG(1,K31,2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       RETURN
  500     K2=K3
       L2=K4
       DO 501 JP=2,JQ

       CALL TFOLD(1,1,Uf,k2)
       CALL TFOLD(0,L2-K2+1,Uf,k2)

  501     K2=K2/2
       LL1=L2-K2+1
       A=Uf(LL1)+Uf(LL1+2)
       Vf(K7+1)=2.0*(-Uf(LL1)+Uf(LL1+2))
       Vf(1)=2.0*(A+Uf(LL1+1))
       Vf(K31)=2.0*(A-Uf(LL1+1))
       Vf(N11+1)=2.0*Uf(2)
       DO 502 JP=2,JQ
       Uf(K2+1)=2.0*Uf(K2+1)
       JSL=K2+1

       CALL TFOLD1(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       LL1=LL1-K2
       Uf(LL1)=-2.0*Uf(LL1)

       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

  502     K2=K2+K2
       Vf(1)=Vf(1)*A5
       Vf(K31)=Vf(K31)*A5
       RETURN
  600     Uf(1)=Uf(1)*A5
       Uf(K31)=Uf(K31)*A5
       K2=K3
       DO 601 JP=2,JQ

       CALL TFOLD(0,K31-K2,Uf,k2)
       CALL TFOLD(1,K31,Uf,k2)

  601     K2=K2/2
       LL1=K31-K2
       A=Uf(LL1)+Uf(LL1+2)
       Vf(1)=2.0*(A+Uf(LL1+1))
       Vf(K31)=2.0*(A-Uf(LL1+1))
       Vf(K7+1)=2.0*(-Uf(LL1)+Uf(LL1+2))
       Vf(N11+1)=2.0*Uf(K3+2)
       DO 602 JP=2,JQ
       JSL=K31+K2
       Uf(JSL)=2.0*Uf(JSL)

       CALL TFOLD1(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       LL1=LL1-K2
       Uf(LL1)=-2.0*Uf(LL1)

       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

  602     K2=K2+K2

       CALL REVNEG(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       idum = k3

       CALL NEG(2,idum,2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       K2=K4

       CALL TFOLD(1,1,Vf,k2)

       Vf(K4+1)=0.0

       RETURN
       END

