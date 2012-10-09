C-------------------------------------------------
C        A.Klypin  1997,2009 (NMSU)
C                            Bound-Density-Maxima code:
C                            MPI version.
C                            Find  maxima of mass in spheres of radius Rsearch1 
C                            for particles in a periodic cube     of size     Box
C                           Coordinates: 0-NGRID,  Box =NGRID
C       All variables and constants are scaled to the Hubble constant of your
C       model. Only FINAL values are rescaled to H=100km/s/Mpc
C       More detailed description of the algorithm and the code can be
C       found in file BDM_describe
C--------------------------------------------------   
C                     MPI:   (Nx,Ny,Nz) mpi tasks
C                             each task reads all the particle data and
C                             selects only partiocles, which are in its domain
C                             plus buffer of width of Rhalo_maximum
C                             Computational box split in equal-volume 
C                             parallelipipeds. 
C--------------------------------------------------   
C       Arrays Label and Lst provide quick access to neighbor particles
C       Label(i,j,k) = number of the first particle in
c                      the cell (i,j,k), if Label .ne. 0
C                    = 0 - empty cell
C       Cell         = size of a cell for Label
C       Lst(i)       = number of the next particle in the cell,
C                      where i-th was found.
C                      If Lst=0 - last point in the cell was found.
C     
      Module Structures

C       Nc - number of centers, Np - number of dm particles
C       Nrad - number of shells for halos
C
C                       Configuration parameters
      Integer, PARAMETER  :: NROW =4096,       ! maximum number of particles in 1D
     &                                          NGRID =256   ,       ! zero-level mesh in 1D
     &                                          Nrad     =   90       ! Number of shells for halo profile                   
      Real, PARAMETER        ::
     &                                          FracSearch = 15.     ! Factor of Rmax for small halos

      PARAMETER (NPAGE = NROW**2)    ! number of particles in a record
      PARAMETER (NMAX  =NGRID/2)
      PARAMETER (NRECL= NPAGE*6)  ! number of words in a record
      PARAMETER (NARR  =MAX(2*NROW+1,NGRID+1))   !      need in FFT
      PARAMETER (NF67    =MAX(NROW,NGRID/2))
      Real, PARAMETER :: MaxMem =8.

      Integer*4 :: Np               ! number of particles on this node
      Integer*8 lspecies(10)
      Real ::           AEXPN,ASTEP,PARTW, Om0,Oml0,hubble,hsmall,
     &                    wspecies(10),
     &                    Xscale,Vscale,Dscale,aMassOne
      Integer ::      Nspecies,Nseed,ISTEP,Nx,Ny,Nz,jStep
      CHARACTER*45           HEADER

      Real,         ALLOCATABLE,   DIMENSION(:) :: Xp,Yp,Zp,VXp,VYp,Vzp
      Integer*8, ALLOCATABLE, DIMENSION(:) :: id_part
      Character*80 :: fNames

       Real        ::  Box,xL,xR,yL,yR,zL,zR,dRbuff

      PARAMETER (Nrad_lim =Nrad)

C                       Parameters to tune halo finder
      PARAMETER (dClose  =    75. )         ! Only one seed within distance
                                           ! dRdubl/dClose will be kept; dRdubl is distance
                                           ! to remove duplicates.
	                                   ! This is used in CATALOG.
                                           ! It is maximum of the radius of the first bin
                                           ! Rad(1) and  dRdubl/dClose

      PARAMETER (nFirst  =    15 )          ! minimum number of particles inside 
                                           ! radius for chosing first bin for overdensity
                                           ! and mass profiles
      PARAMETER (mFirst  =    30)          ! minimum number of particles in side 
                                           ! radius for chosing first bin for radial velocity
					   ! profile Vrad. This is for <Vradial>=0.5Vcirc test,
                                           ! which used to truncate outer radius of halos

      PARAMETER (dRmaxrot =    2.)         ! defines minimum concentration for a halo when
                                           ! deciding on escape velocity.
                                           ! Radius of maximum circular velocity is
                                           ! limited to be smaller than Rhalo/dRmaxrot.
                                           ! Concentration is larger than dRmaxrot*2.15

      PARAMETER (nmin     = 30 )           ! minimum number of particles, which defines
                                           ! radius of a sphere used for bulk velocity of a halo
 
      PARAMETER (iCentr   = 2 )           ! minimum bin number for  bulk velocity of a halo:
                                           ! Radial bin for the velocity is
                                           ! the largest of iCentr and radius containing nmin
 
c                       Optimization parameters
      PARAMETER (Nm      =    -5)          ! define size of linker mesh
      PARAMETER (Nlinker =    250)         ! 
      PARAMETER (Nb      =    Nlinker-Nm)  ! define size of linker mesh
      PARAMETER (mDynamic = 15000 )     ! dynamical scheduling for OpenMP
      PARAMETER (kDynamic = 50 )     ! these parameters define length of
                                     ! chunk of halos to assign for 
                                     ! one OMP thread
                                     ! Integer  4 <->8
      Integer*4,  ALLOCATABLE :: Lst(:)
      Integer*4,  ALLOCATABLE :: Label(:,:,:)
      Integer*4, ALLOCATABLE, DIMENSION(:) :: iRadius,Lab

      Real,  ALLOCATABLE, DIMENSION(:,:) :: Wxc,Wyc,Wzc,
     &                                        Vrmc,Rrc,Vrada
      Integer,  ALLOCATABLE,DIMENSION(:,:) :: Nbc,NbcM,
     &                                            NbcG,NbcS
      Real, ALLOCATABLE, DIMENSION(:) :: Xm,Ym,Zm,Amc,Rmc,Vrm,RmaxV

      Integer :: Nmx,Nbx,Nmy,Nby,Nmz,Nbz
      Real     ::  RadSr1,RadSr2,aNvir,Cell,
     +             Rad(Nrad),Rad2(Nrad),weightSmall,
     &            TotalMemory

      End Module Structures

!-------------------------------------------------------------------
      Program BDM

      use Structures
      include 'mpif.h'

      Real             INPUT,Mhalo
      Real             Overdens(Nrad),Mass(Nrad),MnMass(Nrad)
      COMMON /OVER/ Ovlim,Ovcnt
      DIMENSION      Radsc(Nrad)
      Character ::  fName1*80,fName2*80

      write (*,*)  ' <=== IN MPI_ART ==>'                                                                                  
      TotalMemory  = 0.5   ! give initial memory in Gb
      !irank =41
      mpisize =0
                                                                                                                           
      CALL mpi_init(ierr)                                 
      CALL mpi_comm_size(MPI_COMM_WORLD, mpisize, ierr )    
      CALL mpi_comm_rank(MPI_COMM_WORLD, irank, ierr )           
      node = irank +1                                                                                                      
      write (*,*)  '  size,node=', mpisize,node                                                                                 
                                                               ! Read data  -----------------------------
      CaLL ReadParticles(node,mpisize) ! Open files and read data for each node
          write (13,*) '     time after SetWeights =',seconds()
                                                            ! Set scaling factors and parameters

      Ovlim   = OverdenVir()        ! set virial overdensity
      Ovcnt   = Ovlim   ! Min. Center Overdensity for Halos 
                                   ! scale factor to get number of virial particles
      aNvir  = 4.1888*Ovlim*(NROW/float(NGRID))**3/weightSmall
      Amassl  = AmassOne*5.  ! Minimum halo mass in Msun/h 

      If(node==1)Then
      Rsearch1=INPUT('Enter comoving search radius(Mpc/h)       => ')
         write(13,*)
      Rsearch2=INPUT('Enter smaller radius(Mpc/h) of final halos=> ')
         write(13,*)
      EndIf
      CALL MPI_BCAST(Rsearch1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Rsearch2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      Rminhalo=0.0001   ! min.radius for halos(Mpc/h) 
                             
      If(node==1)Then
      Toohot  =INPUT('Enter rejection velocity limit (V/Vescape)=> ')
         write(13,*)
      EndIf
      CALL MPI_BCAST(Toohot,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      If(node==1)Then
      dRdubl  =INPUT('Distance to check for Velocity duplicates => ')
         write(13,*)
      EndIf
      CALL MPI_BCAST(dRdubl,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      fVdubl  =0.15    ! Define duplicates if  (v1-v2)/Vrms   <   

       Boxx =Box    !  store box size in Mpc
C     ------------------------------------------ Constants ---------
           Omcold= Om0 
           Omhot = 0.
           Amassl= Amassl/hsmall     ! Scale to real masses
           Radius= Rsearch1/hsmall  ! Comoving Search radius in real Mpc
           Rminhalo=Rminhalo/hsmall
           dRdubl= dRdubl/hsmall
C                         Dscale is a factor used to convert particles to mass 
C                                    =  mass of a particle for Omega=1 model
C                         Xscale is a factor to convert PM coordinates into Mpc
C                                    = Box / Ngrid = length of a PM cell in Mpc
C                        Vscale is a factor to convert PM pomenta to V (km/sec)
C                                    = the Hubble velocity at the distance 
C                                         of 1 PM cell
C                the factor AEXPN is needed to go from momenta to pec.velocity
           Radius=Radius/Xscale   ! comoving radius in cell units
           RadSr1 =Rsearch1/hsmall/Xscale ! max Min search radii
           RadSr2 =Rsearch2/hsmall/Xscale           
           Box   =NGRID           ! all internal variables are in Grid units
           Riter =0.0001          ! Error for iteration of radius in cell units
           Rmax  =Radius          ! Treat maxima Identical  if dR<Rmax
           Cell   = max(Radius,(xR-xl)/Nlinker)
           write(13,*) ' R search=', RadSr1,RadSr2,Radius
           write (13,*) ' Cell size=',Cell,
     +               ' Search Radius in cell units=',Radius
         DtoNscale = Om0*Dscale ! scale dark matter particles
         ANlimit =aNvir*weightSmall ! factor to calculate
                                      ! central overdensity
                                      ! It will be 1 times larger than
                                      ! virial overdensity
           R0    =0.004
           Rnrad =R0*((Nrad-1)/4.+1.)**2 
           write (13,'(/" The first and the last bins have radii:",
     &                   2f8.4," Mpc/h, comoving")')  R0,Rnrad
           write (13,*)'   If you want to change shells, enter a factor'
           write (13,'(A,$)') '    or 1 to keep current radii : '
           If(irank==0)Then
               read  (*,*) factor
           EndIf
      CALL MPI_BCAST(factor,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
           If(abs(factor-1.).ge.1.e-3.and.factor.gt.0.)Then
              R0    =0.004 * factor
              Rnrad =R0*((Nrad-1)/4.+1.)**2           
              write (13,*) '   The first and the last radii for halos'
              write (13,*) '   will be',R0,Rnrad,' Mpc/h, comoving'
           EndIf
           Do ir=1,Nrad                            ! Radii for halo shells
              Rad(ir) =R0/hsmall*((ir-1)/4.+1.)**2/Xscale ! comoving
              Rad2(ir)=Rad(ir)**2
           EndDo
           write (13,'(3(a,g12.3))') 
     +                ' Final Radii(Mpc/h): 1st=',Rad(1)*Xscale*hsmall,
     +                ' last=',Rad(Nrad)*Xscale*hsmall,' Cell=',Cell

      Call LimList                       ! find limits for linker-list
          !  call mpi_barrier(MPI_COMM_WORLD,ierr)
          !  CALL mpi_finalize(ierr)                                                                                            
          !  STOP  
      ALLOCATE (Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
      ALLOCATE (Lst(Np))
         Call Memory((Nbx-Nmx+1)*(Nby-Nmy+1)*(Nbz-Nmz+1))
         Call Memory(Np)

      Call List                                ! make linker List, Label

         write (13,*)  'List made. time=',seconds()
      CALL ReadCent(Ncentr,Boxx,0,irank)     ! count initial centers
      Nc = Ncentr
      ALLOCATE (Xm(Nc),Ym(Nc),Zm(Nc),Amc(Nc),Rmc(Nc))
      Call Memory(5*Nc)
      CALL ReadCent(Ncentr,Boxx,1,irank)     ! set initial centers
          write (13,'(2(a,i11),a,f9.2)') ' Nparticles=',Np,
     &                        ' Ncenters=',Ncentr,' time=',seconds()

C----------------------------- Find maxima
      ALLOCATE (Lab(Ncentr))
      Call Memory(Ncenter)
      CALL Pairs(Ncentr,Riter)       ! find maxima 
             write(13,*)     '  Time after Pairs=',seconds()
      DEALLOCATE (Lab)
      Call Memory(-Ncenter)
      write (13,16)  Amassl*hsmall,ANlimit,Rminhalo*hsmall,
     &          dRdubl/dClose*hsmall*1.e3,Rsearch1,Rsearch2,Toohot,Cell
      If(node==1)
     & write (2,16)  Amassl*hsmall,ANlimit,Rminhalo*hsmall,
     &          dRdubl/dClose*hsmall*1.e3,Rsearch1,Rsearch2,Toohot,Cell
!      write (20,16) Amassl*hsmall,ANlimit,Rminhalo*hsmall,
!     &          dRdubl/dClose*hsmall*1.e3,Rsearch1,Rsearch2,Toohot,Cell
 16   format(5x,'Small halos with Mass(M_sun/h) less than',
     .       g10.3,'(N_eff<',g11.3,') were removed',/
     .       5x,'Minimum Radius of Halo(Mpc/h)=',f7.4,
     &          ' Minimum distance beteen halos(kpc/h)=',f7.4,/
     .       5x,'Comov.Search radius(Mpc/h)=',f7.4,
     .         ' Second radius(Mpc/h)=',f7.4,/
     .       5x,'Particles with V>',f5.1,'V_escape were removed',
     .       ' Cell (in grids) for linker list=',f8.3)
           write (13,*)  'Pairs made'
      CALL SmallRem(ANlimit,Ncentr)  ! Remove small halos
           write (13,*)  'Small Rem made. Ncentr=',Ncentr
!           call mpi_barrier(MPI_COMM_WORLD,ierr)
!           CALL mpi_finalize(ierr)                                                                                            
!            STOP  
C                           Final touch: if Rsearch1 and Rsearch2 are different
C                           then  improve position of halos by shrinking
C                           the search Radius to Rsearch2 in Nit steps

      Rsearch2 = Rsearch2/4.          ! reduce the small radius 
c      RadSr2 =Rsearch2/hsmall/Xscale  ! for more accurate search of center
      ALLOCATE (Lab(Ncentr))
      Call Memory(Ncentr)
      IF(ABS(Rsearch1-Rsearch2).gt.0.001*Rsearch1)Then
         Nit =2
         hNit =(Rsearch1/Rsearch2)**(1./Nit)
         !If(iDebug.eq.1)write (13,*) ' Nit =',Nit,' Factor:',hNit
        Do ii =1,Nit
            n1 = 0
            r1  =0.
            Do i=1,Ncentr
      	     If(Rmc(i).gt.RadSr2)Then
                    n1     = n1 +1
                    Rmc(i) = max(Rmc(i)/hNit,RadSr2)
                    r1     = r1 +Rmc(i)
                 Endif
             Enddo 
!          write (13,*)'Iterate:',
!     &               ' Mean Comoving Radius for halos(Mpc/h)=',
!     &                   r1/(max(n1,1))*Xscale*hsmall,' n=',n1
          CALL Pairs(Ncentr,Riter) 
        EndDo
        Radius=Rsearch1
      EndIf

       write(13,*) '   ---- do catalog --'
      dRdubl =dRdubl/Xscale                ! scale to cell units
      CALL CATALOG(Ncentr,dRdubl)      ! Remove duplicates
      DEALLOCATE (Lab)
      Call Memory(Ncentr)
      
       write(13,*) '   ---- end catalog --'
      DEALLOCATE (Label)
         Call Memory(-(Nbx-Nmx+1)*(Nby-Nmy+1)*(Nbz-Nmz+1))
      Celloptimal = (Rad(Nrad)/FracSearch)
      CellOld     = Cell
      Cell = max(Celloptimal,(xR-xL)/float(Nlinker))
         ff   = CellOld/Cell      ! scaling factor to new limits
         Nmx  = (Nmx+5)*ff-5
         Nmy  = (Nmy+5)*ff-5
         Nmz  = (Nmz+5)*ff-5
         Nbx  = (Nbx-5)*ff+5
         Nby  = (Nby-5)*ff+5
         Nbz  = (Nbz-5)*ff+5
      
      ALLOCATE (Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz), Stat =iStat)
         Call Memory((Nbx-Nmx+1)*(Nby-Nmy+1)*(Nbz-Nmz+1))
      Call List              ! make linker List, Label

       write (13,20) Np,Ncentr,Box*Xscale*hsmall,Ovlim,Cell
       If(node==1)
     &  write (2,20) Np,Ncentr,Box*Xscale*hsmall,Ovlim,Cell
!       write (20,20)Np,Ncentr,Box*Xscale*hsmall,Ovlim,Cell
20         Format(5x,'Statistics for ',i12,' points. Maxima=',i9,/
     .         5x,'Box(Mpc)=',f6.1,
     .            5x,'Overdensity Limit=',f6.1,
     .             ' Cell for linker list(grids)=',f8.3)
          write (13,*) ' time  before Spairs(sec)=',seconds()
      ALLOCATE (Wxc(Nrad,Ncentr),Wyc(Nrad,Ncentr),Wzc(Nrad,Ncentr))
      Call Memory(3*Nrad*Ncentr)
      ALLOCATE (Vrmc(Nrad,Ncentr),Rrc(Nrad,Ncentr),Nbc(Nrad,Ncentr))
      Call Memory(3*Nrad*Ncentr)

      ALLOCATE (Vrada(Nrad,Ncentr)) !,NbcM(Nrad,Ncentr),NbcG(Nrad,Ncentr),
!     &            NbcS(Nrad,Ncentr),Stat=iStat)
      CALL SPairs(Ncentr)  ! Accumulate  statistics 
          write (13,*) ' time for Spairs (sec)=',seconds()
!            call mpi_barrier(MPI_COMM_WORLD,ierr)
!           CALL mpi_finalize(ierr)                                                                                            
!            STOP  
      Amassl  = AmassOne*8.  ! Minimum halo mass in Msun/h 
C-----------------------------      Auxiliary Settings
      nfalse =0
      Do ir=1,Nrad                  ! all scaled to hubble, not h^(-1)
         Radsc(ir) =Rad(ir) *Xscale ! scale radii to comoving Mpc
         MnMass(ir)=Om0*4.188*2.746e+11*hubble**2*Radsc(ir)**3 
c                                   mass at the mean density
      EndDo
      write (13,12) (Radsc(i)*hsmall*1000.,i=1,Nrad)
      If(node==1)
     &  write (2,12) (Radsc(i)*hsmall*1000.,i=1,Nrad)
!      write (20,12) (Radsc(i)*hsmall*1000.,i=1,Nrad)
 12   Format(5x,'Radii of shells in comoving kpc/h:',/20(4x,10f8.2/))
C---------------------------- Process Halos : Unbind, Dublicates 
      ALLOCATE (Vrm(Ncentr),RmaxV(Ncentr),iRadius(Ncentr))
      Call Memory(3*Ncentr)
      aa =1.          !! send this to wake up all MPI tasks
      CALL MPI_BCAST(aa,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL CLEAN(Ncentr,Radsc,MnMass,Toohot,Rminhalo)
      write (13,*)  '  CLEAN: Ncenters=',Ncentr
      write (13,*) ' time after CLEAN(sec)=',seconds()
      aa =1.          !! send this to wake up all MPI tasks
      CALL MPI_BCAST(aa,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      IF(fVdubl .gt.0. .and. dRdubl.gt.0.)Then    !  Remove VelDublicates
      ALLOCATE (Lab(Ncentr)) 
      Call Memory(Ncentr)
         CALL RemVelDublA(Ncentr,dRdubl,fVdubl,Xscale,Vscale)
      DEALLOCATE (Lab) 
      aa =1.          !! send this to wake up all MPI tasks
      CALL MPI_BCAST(aa,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      Call Memory(-Ncentr)
         write (13,25) Ncentr,dRdubl*Xscale*hsmall ,fVdubl
         If(node==1)
     &    write (2,25) Ncentr,dRdubl*Xscale*hsmall ,fVdubl
!         write (20,25)Ncentr,dRdubl*Xscale*hsmall ,fVdubl
 25       Format(5x,'Number of centers=',i8,
     &                   2x,'Velocity Duplicates: dR(Mpc/h)<',f6.3,
     &                   2x,'dV/Vrms<',f6.3)
       ENDIF

       write(13,*) ' Step=',jStep,ISTEP
       write(fName1,'(a,i3.3,a)')'CATALOGS/CatshortA.',ISTEP,'.DAT' 
C------------------------------- Print Results
      iCount = 0
      write(13,'(a,1P,6g12.3)') ' Boundaries =',xL,xR,yL,yR,zL,zR

      Do inode = 1,Nx*Ny*Nz
       Close(2)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      If(node==inode)Then
         Open(2,file=TRIM(fName1),position='append')
!      write ( 2,*)' ---- New node',node
      !write (20,13)
 13   Format(6x,'Coordinates Mpc',12x,'Velocity km/s',7x,
     & 'Mass     Radius Vrms(3D) Vcirc Npart(<Rad) Rmax(kpch)',
     & ' Conctr Bins')
      ! write (20,14)
 14   Format('R(kpc/h)',5x,' N(<R)    M/Msun   Overdens   Vrms  Vcirc',
     &                  ' Dens(Msun/Mpc^3) Nsmall Nlrg Nbig Vrad(km/s)')
      Do i=1,Ncentr          ! Scale all to physical units
              If(Xm(i).ge.xL+dRbuff.and.Xm(i).lt.xR-dRbuff
     & .and.Ym(i).ge.yL+dRbuff.and.Ym(i).lt.yR-dRbuff
     & .and.Zm(i).ge.zL+dRbuff.and.Zm(i).le.zR-dRbuff)Then         
         Xc    =(Xm(i)-xshift)  *Xscale
         Yc    =(Ym(i)-yshift)  *Xscale
         Zc    =(Zm(i)-zshift)  *Xscale
         irad  =iRadius(i)
         Mhalo =Amc(i)
         Rhalo =Rmc(i)
         RmaxH =RmaxV(i)
         Vcirc =Vrm(i)

!         iS      =NbcS(irad,i)   ! small-size particles
!         iM      =NbcM(irad,i)   ! Medium-size particles
!         iG      =NbcG(irad,i)   ! Giant  particles
         Nnp   =INT(Mhalo/(Dscale*Om0))
         Do ir =1,Nrad
            Summass     = Om0*Nbc(ir,i)
            Mass(ir)    = Dscale*Summass
            Overdens(ir)= Mass(ir)/MnMass(ir)
         EndDo
         im =Nbc(Nrad,i)
             Do ir =Nrad-1,1,-1
                if(Nbc(ir,i).ne.im)Then
                   iLast =ir
                   goto 220
                EndIf 
             EndDo
 220         iLast = min(Nrad,iLast+1)
       If(irad.eq.0.or.Rhalo.le.Rminhalo
     .      .or.Mhalo.le.Amassl)Then
              nfalse=nfalse+1
!            irr =Max(irad,2)         ! don't take V from the first bin
!            VXc =Wxc(irr,i)   *Vscale
!            VYc =Wyc(irr,i)   *Vscale
!            VZc =Wzc(irr,i)   *Vscale
!            Nnp =INT(Mhalo/(Dscale*Om0))
!            write(30,45) Xc*hsmall,Yc*hsmall,Zc*hsmall,
!     &                      VXc,VYc,VZc,Mhalo*hsmall,
!     &                      Rhalo*hsmall*1000.,Vrmc(irr,i),
!     &                      Vcirc,Nnp,
!     &                      RmaxH*hsmall*1000.,Rhalo/RmaxH*2.15,iLast
       EndIf
         If((irad.ne.0.and.Rhalo.gt.Rminhalo)
     .                .and.Mhalo.gt.Amassl
     .                )Then
            irr = Max(irad,2)
            irr2 = 2  !Max(irad,2)         ! don't take V from the first bin
 10         if(Nbc(irr2,i).lt.nmin.and.irr2.lt.iCentr)Then
               irr2 =irr2+1
               goto 10
            EndIf 

            VXc =Wxc(irr2,i)   *Vscale
            VYc =Wyc(irr2,i)   *Vscale
            VZc =Wzc(irr2,i)   *Vscale
            Nnp =INT(Mhalo/(Dscale*Om0))

!             write (20,45) Xc*hsmall,Yc*hsmall,Zc*hsmall,
!     &                      VXc,VYc,VZc,Mhalo*hsmall,
!     &                      Rhalo*hsmall*1000.,Vrmc(irr,i),
!     &                      Vcirc,Nnp,
!     &                      RmaxH*hsmall*1000.,Rhalo/RmaxH*2.15,iLast

             write ( 2,46) Xc*hsmall,Yc*hsmall,Zc*hsmall,
     &                      VXc,VYc,VZc,Mhalo*hsmall,
     &                      Rhalo*hsmall*1000.,Vrmc(irr,i),
     &                      Vcirc,Nnp,
     &                      RmaxH*hsmall*1000.,Rhalo/RmaxH*2.15


!            Do ir=1,iLast
!               If(ir.eq.1)Then
!                  Vvc =Vrmc(ir,i)
!                  Volume =4.188*Radsc(ir)**3
!                  DenDm=Mass(ir)/Volume
!                  If(Nbc(ir,i).le.2)Then ! too few particles
!                      Radmean =Radsc(ir)/2.**0.3333
!                  Else         ! take real mean radius
!                      Radmean  =Rrc(ir,i)*Xscale/Nbc(ir,i)
!                  EndIf
!                  Vcr    =sqrt(Mass(1)/Radsc(1)/AEXPN)*6.58e-5
!               Else              ! remove contribution from inner rad
!                 Vvc =sqrt(MAX((Vrmc(ir,i)**2*Nbc(ir,i) -
!     .              Vrmc(ir-1,i)**2*Nbc(ir-1,i))/
!     .              max(Nbc(ir,i) -Nbc(ir-1,i),1),1.e-10))
!c                  Vvc =(Vrmc(ir,i)*Nbc(ir,i) -  !radial velocity
!c     .              Vrmc(ir-1,i)*Nbc(ir-1,i))/
!c     .              max(Nbc(ir,i) -Nbc(ir-1,i),1)
!                  Volume =4.188*(Radsc(ir)**3- Radsc(ir-1)**3) 
!                  DenDm  =(Mass(ir)-Mass(ir-1)) /Volume
!                  If(Nbc(ir,i) -Nbc(ir-1,i).le.2)Then ! too few particles
!                     Radmean=((Radsc(ir)**3+Radsc(ir-1)**3)/2.)**0.333
!                  Else         ! take real mean radius
!                     Radmean=(Rrc(ir,i)-Rrc(ir-1,i))*Xscale/
!     .                        (Nbc(ir,i)-Nbc(ir-1,i))
!                  EndIf
!                  Vcr    =sqrt(Mass(ir)/Radsc(ir)/AEXPN)*6.58e-5
!                EndIf
!               write (20,40) Radmean*hsmall*1000.,
!     +                       Nbc(ir,i),Mass(ir)*hsmall,Overdens(ir),
!     +                       Vvc,Vcr,DenDm/hsmall**2
!!     +                       ,NbcS(ir,i),NbcM(ir,i),NbcG(ir,i),
!     &                      ,Vrada(ir,i)*Vscale
!            EndDo  ! end loop ir
45          Format(F9.4,2F9.4,3F8.1,g11.3,f8.2,2f7.1,I9,f8.1,f8.2,i4)
46          Format(F9.4,2F9.4,3F8.1,g11.3,f8.2,2F7.1,i9,f8.1,f8.2,i3)
47         Format(F8.3,2F8.3,g11.3,f8.2,i9,i8,i7,i5)
40          Format(g11.4,I9,g11.3,g10.3,2f7.1,1P,g10.3,g11.3)
41          Format(f7.4,I7,I6,g11.3,g10.3,2f7.1,f6.3,g10.3,4G10.3)
         Endif
         EndIf    ! end boundaries test
      EndDo    !  end loop Ncentr
      write (13,50) Ncentr-nfalse,nfalse
      EndIf
      EndDo   !  inode
c      write (2,50) Ncentr,nfalse
 50   Format('---------- Final Catalog has ',i7,' halos ------------',
     .      /'                    ',i7,' were too puffy or too small') 
      close (2)
      close (13)
      !close (20)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      CALL mpi_finalize(ierr) 

      Stop
      End Program BDM
!-------------------------------------------------------------------------------- 
!                          Define configuration of files
!                          Read control infor from the first file
!                          allocate arrays
      Subroutine ReadParticles(node,mpisize)
      Use Structures
      include 'mpif.h'
      character*80 :: Name,fName1,fName2

      write(fName1,'(a,i4.4,a,i4.4,a)')'CATALOGS/outputA.',node,
     &                '.',jStep,'.dat'
      open(13,file=TRIM(fName1),buffered='NO')
      write(13,*) ' My rank=',node,mpisize

      If(node==1)Then
         write(*,'(a,$)')' Enter snapshot number ='
         read(*,*)jStep
       EndIf
      CALL MPI_BCAST(jStep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

                                     ! Open files
      write(Name,'(a,i4.4,a,i4.4,a)') 'PMss/PMss.',node,'.',jStep,'.DAT'
      open(1,file=Name,form='unformatted')
      If(node==1)Then      
         write(fName1,'(a,i3.3,a)')'CATALOGS/CatshortA.',jStep,'.DAT' 
         !write(fName2,'(a,i3.3,a)')'CatalogB.',node,'.DAT' 
         Open( 2,file=TRIM(fName1),position='append')! short list of halos
         !Open(20,file=TRIM(fName2),Status='UNKNOWN') ! catalog of halos
       EndIf
        Read(1)HEADER
        Read(1)AEXPN,ASTEP,ISTEP,NROWC,NGRIDC,Nspecies,  
     &                            Nseed,Om0,Oml0,hubble,Box
        Read(1) k,Nx,Ny,Nz,dRbuff
        Read(1) xL,xR,yL,yR,zL,zR
        xL = xL -1.
        xR = xR -1.
        yL = yL -1.
        yR = yR -1.
        zL = zL -1.
        zR = zR -1.
                                               ! Set scales
       weightSmall = 8.
       wspecies(1) = (float(NGRID)/float(NROW))**3*weightSmall
           hsmall = hubble            ! Hubble/100.
           Box     = Box /hsmall       ! Scale to real Mpc          
           Xscale = Box/Ngrid        ! Scale for comoving coordinates
           Vscale = 100.*hsmall*Xscale/AEXPN ! Scale for velocities
           Dscale = 2.746e+11*hsmall**2*(Box/NROW)**3 *weightSmall ! mass scale
           aMassOne= Om0*Dscale*hsmall

      write(13,'(a)')HEADER
      write(13,'(4x,a,f9.4,5x,a,2i5)')'aexpn =',AEXPN,' Step=',ISTEP
      write (13,*) '   Box size (Mpc/h)      =  ',Box
      write (13,*) '   Ngrid      =',NGRID
      write (13,*) '   Nrow       =',NROW
      write(13,*) '    Boundaries =',xL*Xscale,xR*Xscale
      write(13,*) '    Boundaries =',yL*Xscale,yR*Xscale
      write(13,*) '    Boundaries =',zL*Xscale,zR*Xscale
      write(13,*) '    Buffer zone=',dRbuff*Xscale
      If(node==1)
     &    WRITE (2,100) HEADER, AEXPN,ASTEP,ISTEP,
     +               NROW,NGRID,Om0,Oml0,hubble,dRbuff*Xscale
100      FORMAT(1X,'Header=>',A45,/
     +         ' A=',F8.3,' Step=',F8.3,/
     +         ' I =',I4,' Nrow=',I4,' Ngrid=',I4,/' Omega_0=',f6.3,
     +         ' Omega_L=',f6.3,' h=',f6.3,' Buffer zone (Mpch)=',f6.3)

      if(ISTEP .ne. jStep)Stop ' -- wrong time step --'
      if(k .ne. node)Stop ' -- wrong node --'
      IF(mpisize .NE. Nx*Ny*Nz) THEN                                                                                        
           write(13,*) 'mpisize .NE. n_nodes',mpisize,Nx*Ny*Nz 
            CALL mpi_finalize(ierr)                                                                                            
            STOP  
      EndIf
      write(13,'(4x,a,3i5)')'Configuration: in each direction=',
     &                           Nx,Ny,Nz
      write (13,'(6x,"Mass of smallest particle",/9x,
     +    "in units M_sun/h is   =",3x,g10.3)') 
     +      aMassOne
      write(13,'(6x,"first species=",2g11.3)')wspecies(1),weightSmall
      If(Nspecies.gt.1)Then
      write (13,'(6x,"Mass of largest particle",/9x,
     +    "in units M_sun/h is   =",3x,g10.3)') 
     +   Om0*Dscale*wspecies(Nspecies)*hsmall/(NGRID**3/FLOAT(NROW)**3)
      EndIf

        read(1)Np       ! total number of particles in this file
        write(13,*) '     Number of particles in this node   =',Np
        ALLOCATE(Xp(Np),Yp(Np),Zp(Np),
     &                   VXp(Np),VYp(Np),VZp(Np),id_part(Np))
        Call Memory(8*Np)
        ioffset =0
      Do                                              ! read coordinates and velocities
         read(1,iostat=ierr)mm         ! number of particles in this record        
         If(ierr /= 0) exit
         write(13,*) '     Number of particles in this record =',mm
         read(1)                  (Xp(i),Yp(i),Zp(i),              
     &                               VXp(i),VYp(i),VZp(i),      
     &                              id_part(i),i=1+ioffset,mm+ioffset)
          ioffset =ioffset + mm
      EndDo

      DEALLOCATE(id_part)
        Call Memory(-2*Np)
      If(ioffset.ne.Np)Then
            write(13,*) '  Wrong Number of particles:',Np,
     &                        ' should=',ioffset
            CALL mpi_finalize(ierr)                                                                                            
            STOP  
      EndIf        
           lspecies(1) = Np
          Nspecies    = 1        
      write(13,*) ' Done Reading Points. N=',Np
 
!$OMP PARALLEL DO DEFAULT(SHARED) 
!$OMP+ PRIVATE (i)                 
!$OMP+ REDUCTION(MAX:xmax,ymax,zmax) 
!$OMP+ REDUCTION(MIN:xmin,ymin,zmin) 
      Do i   =1,Np 
            Xp(i)  =Xp(i)-1.
            Yp(i)  =Yp(i)-1.
            Zp(i)  =Zp(i)-1.
!      if(Xp(i).le.xzero)write (13,'("X:",3f11.6,i12)')Xp(i),Yp(i),Zp(i),i
!      if(Yp(i).le.xzero)write (13,'("Y:",3f11.6,i12)')Xp(i),Yp(i),Zp(i),i
!      if(Zp(i).le.xzero)write (13,'("Z:",3f11.6,i12)')Xp(i),Yp(i),Zp(i),i
!            If(Xp(i).gt.NGRID)Xp(i)  =Xp(i)-NGRID
!            If(Yp(i).gt.NGRID)Yp(i)  =Yp(i)-NGRID
!            If(Zp(i).gt.NGRID)Zp(i)  =Zp(i)-NGRID
!            If(Xp(i).le.0.)Xp(i)  =Xp(i)+NGRID
!            If(Yp(i).le.0.)Yp(i)  =Yp(i)+NGRID
!            If(Zp(i).le.0.)Zp(i)  =Zp(i)+NGRID
            If(Xp(i).ge.NGRID)Xp(i)  =Xp(i)-1.e-4
            If(Yp(i).ge.NGRID)Yp(i)  =Yp(i)-1.e-4
            If(Zp(i).ge.NGRID)Zp(i)  =Zp(i)-1.e-4
!      if(Xp(i).ge.xbig)
!     &            write(50,'(" X:",3f11.6,i12)')Xp(i),Yp(i),Zp(i),i
!      if(Yp(i).ge.xbig)
!     &            write(51,'(" Y:",3f11.6,i12)')Xp(i),Yp(i),Zp(i),i
!      if(Zp(i).ge.xbig)
!     &            write(52,'(" Z:",3f11.6,i12)')Xp(i),Yp(i),Zp(i),i
               xmin = min(Xp(i),xmin)
               xmax =max(Xp(i),xmax)
               ymin = min(Yp(i),ymin)
               ymax =max(Yp(i),ymax)
               zmin = min(Zp(i),zmin)
               zmax =max(Zp(i),zmax)
      ENDDO 

             call mpi_barrier(MPI_COMM_WORLD,ierr)
          ! CALL mpi_finalize(ierr)                                                                                            
          !  STOP  


      end Subroutine ReadParticles
C--------------------------------------------------------------
C                        virial overdensity for cosmological model
C                        at different expansion parameter AEXPN
      Function OverdenVir()
C--------------------------------------------------------------
      use Structures

      xx =-(1.-Om0)*AEXPN**3/(Om0+(1.-Om0)*AEXPN**3)
      Ovlim =(178.+82.*xx-39.*xx**2)/(1.+xx)
      write (13,*)  '      Overdensity Delta =',Ovlim
                              ! FOR DARK ENERGY
      ISUGRA = 0   ! 1- for RP
      If(ISUGRA.eq.1)Then
         Qlambda  = 3    ! log(lambda/Gev) for RP
         a = -2.119e-2*Qlambda-0.259
         b = -1.833e-2*Qlambda+0.975
         c = -6.975e-3*Qlambda+9.771e-2
         z = 1./AEXPN-1.
         Omz = 1.-(1.-Om0)*AEXPN**(a+b*abs(z)**c)
         write (13,*) ' a,b,c=',a,b,c
         write (13,*) ' power =',a+b*abs(z)**c
         a = -1.45e-2*Qlambda+0.186
         b =  -1.1e-2*Qlambda +0.22
         Ovlim =178.*Omz**(a+b*Omz-1)
         write (13,*) ' OverdensLimit=',Ovlim
         write (13,*) ' Omega_matter(z) =',Omz
         write (13,*) ' Omega_matter(0) =',Om0
      EndIf

      OverdenVir = Ovlim
c      Ovlim   =INPUT('Enter Overdensity Threshold for Halos     => ')

 
      RETURN
      END
C--------------------------------------------------------------
C                                       Find limits for coordinates
      SUBROUTINE LimList
C--------------------------------------------------------------
      use Structures
      Real*8 :: xmin,xmax,ymin,ymax,zmin,zmax

      write(13,*) '         LimList: N       =',Np
      write(13,*) '         Cell                  =',Cell

      first =weightSmall*1.1   ! add small factor to avoid rounding errors
      
      xmin = 1.e10
      ymin = 1.e10
      zmin = 1.e10
      xmax = -1.e10
      ymax = -1.e10
      zmax = -1.e10
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i) REDUCTION(MIN:xmin,ymin,zmin) 
C$OMP+REDUCTION(MAX:xmax,ymax,zmax)
         Do i=1,Np
!            If(iWeight(i).lt.first)Then
               xmin =min(Xp(i),xmin)
               xmax =max(Xp(i),xmax)
               ymin =min(Yp(i),ymin)
               ymax =max(Yp(i),ymax)
               zmin =min(Zp(i),zmin)
               zmax =max(Zp(i),zmax)
!            EndIf
         EndDo          

         Nmx=INT(xmin/Cell)-5
         Nmy=INT(ymin/Cell)-5
         Nmz=INT(zmin/Cell)-5
         Nbx=INT(xmax/Cell)+5
         Nby=INT(ymax/Cell)+5
         Nbz=INT(zmax/Cell)+5

      write(13,*) '       Particles:',Np
      write(13,*) '               x:',xmin,xmax
      write(13,*) '               y:',ymin,ymax
      write(13,*) '               z:',zmin,zmax
      write(13,*) ' Label limits  x:',Nmx,Nbx
      write(13,*) '               y:',Nmy,Nby
      write(13,*) '               z:',Nmz,Nbz
      write(13,'(a,g12.4)')' Memory required (Mb)=',
     &                    4.*(Nbz-Nmz)*(Nby-Nmy)*(Nbx-Nmx)/1024**2

      Return
      End
C--------------------------------------------------------------
C                          Make linker lists of particles in each cell
      SUBROUTINE List
C--------------------------------------------------------------
      use Structures
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i)
             Do i=1,Np
                Lst(i)=-1
             EndDo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i,j,k)
             Do k=Nmz,Nbz
             Do j=Nmy,Nby
             Do i=Nmx,Nbx
                Label(i,j,k)=0
             EndDo
             EndDo
             EndDo  
      Do jp=1,Np
         i=INT(Xp(jp)/Cell)
         j=INT(Yp(jp)/Cell)
         k=INT(Zp(jp)/Cell)
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
         Lst(jp)      =Label(i,j,k)
         Label(i,j,k) =jp
      EndDo

      write(13,*) '     Linker list is done '

      Return
      End
C--------------------------------------------------------------
C                                 Remove dublicates 
      SUBROUTINE RemVelDublA(Ncentr,dRdubl,fVdubl)
C--------------------------------------------------------------
      use Structures

      Real      Mhalo,Mhalo2
      double precision q1
      Logical  inside

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i)
             Do i=1,Ncentr
                Lst(i)=-1
             EndDo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i,j,k)
             Do k=Nmz,Nbz
             Do j=Nmy,Nby
             Do i=Nmx,Nbx
                Label(i,j,k)=0
             EndDo
             EndDo
             EndDo  
      Do jp=1,Ncentr
         i=INT(Xm(jp)/Cell)
         j=INT(Ym(jp)/Cell)
         k=INT(Zm(jp)/Cell)
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
         Lst(jp)      =Label(i,j,k)
         Label(i,j,k) =jp
      EndDo

      dR2 = dRdubl**2
      fV2 = fVdubl**2
      Vol_2 =Rad(2)**3-Rad(1)**3
      nrem= 0


C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i,Mhalo,irad) REDUCTION(+:nrem)
      Do i=1,Ncentr
        Mhalo  =Amc(i)
        irad   =Max(iRadius(i),2)         
        If(Mhalo.le.1.e-7 .or. irad.eq.0)Then
          Lab(i) =-1
          nrem = nrem +1
        Else
          Lab(i) =1
        EndIf
c        if(Lab(i).eq.1)
      EndDo
      write(13,*) ' RemVelDubl:',Ncentr,Xscale,Vscale,nrem
      nrem =0
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i,jp,Mhalo,Mhalo2,irad,jrad,Rh_1,Rh_2,Rhh,Vrh2,dd,  
C$OMP+dv,VV,irr,q1,Xc,Yc,Zc,i1,j1,k1,i2,j2,k2,ii,jj,kk,inside)
C$OMP+REDUCTION(+:nrem) SCHEDULE(DYNAMIC,mDynamic)
      Do i=1,Ncentr
c        If(i/1000*1000.eq.i)write (13,*)  ' RmDubl i=',i
        If(Lab(i).ne.-1)Then         ! do not take empty centers
        Mhalo  =Amc(i)
        irad   =Max(iRadius(i),2)
        Rh_1   =Rad(irad)
        Vrh2   =(Vrmc(irad,i)/Vscale)**2

        Xc =Xm(i)
        Yc =Ym(i)
        Zc =Zm(i)
c                      limits for Label
         i1=INT((Xc-dRdubl)/Cell)
         j1=INT((Yc-dRdubl)/Cell)
         k1=INT((Zc-dRdubl)/Cell)
            i1=MIN(MAX(Nmx,i1),Nbx)
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
         i2=INT((Xc+dRdubl)/Cell)
         j2=INT((Yc+dRdubl)/Cell)
         k2=INT((Zc+dRdubl)/Cell)
            i2=MIN(MAX(Nmx,i2),Nbx)
            j2=MIN(MAX(Nmy,j2),Nby)
            k2=MIN(MAX(Nmz,k2),Nbz)
C                                        Look for neibhours in nearest cells
         Do kk =k1,k2
         Do jj =j1,j2
         Do ii =i1,i2
            jp =Label(ii,jj,kk)   !  Dark matter Halo 
 10         If(jp.ne.0)Then

          If(i.ne.jp.and.Lab(jp).ne.-1)Then ! do not compare with empty centers
          dd =sqrt((Xm(jp)-Xc)**2 + (Ym(jp)-Yc)**2+ (Zm(jp)-Zc)**2)     
          Mhalo2 =Amc(jp)
          jrad=Max(iRadius(jp),2)

          Rh_2 =Rad(jrad)
          Rhh   =MAX(Rh_1,Rh_2,1.e-8)

          inside = (dd.lt.Rhh).and.
     &               (Mhalo2.gt.Mhalo/1.5.and.Mhalo2.lt.Mhalo*1.5)
          inside = inside.or.
     &               (dd.lt.Rhh/10.).and.
     &               (Mhalo2.gt.Mhalo/3.0.and.Mhalo2.lt.Mhalo*3.0)
          If(inside)Then    ! close distance pair
!             jrad=Max(iRadius(jp),2)
!             dv =(Wxc(jrad,jp)-Wxc(irad,i))**2 +
!     .           (Wyc(jrad,jp)-Wyc(irad,i))**2 +      
!     .           (Wzc(jrad,jp)-Wzc(irad,i))**2
!             VV =MAX(Vrh2,(Vrmc(jrad,jp)/Vscale)**2) ! max Vrms of the two
c             If(dv.lt.(6.-5.*dd/dR2)*fV2*VV)Then   ! close in velocity
c             If(dv.lt.(2.-dd/(dRdubl*Rhh))*fV2*VV)Then ! close in velocity
!             If(dv.lt.fV2*VV.or.dd.lt.dR2/100.)Then ! close in velocity
!             If(dv.lt.fV2*VV)Then ! close in velocity
                 irr = min(10,Nrad)             ! costruct weight to find
                q1 = (Mhalo/Mhalo2)*     ! which halo is larger
     &               (Vrmc(irr,i)/Vrmc(irr,jp))**2*
     &              (Nbc(3,i)/real(max(Nbc(3,jp),1)))*
     &              Nbc(irr,i)/real(max(Nbc(irr,jp),1))
                if(abs(q1-1.).lt.1.e-5)q1=jp/float(i)
                nrem = nrem +1

31           format(f10.4,f9.3,f11.2,2i14,3f10.3,2g11.3,i4,2i5,2i3)               
 30             format(3i6,f10.4,2f9.3,2g12.4,4f10.3,3x,3f9.4)
 35             format(' MM ',2i8,f10.4,4i6,g11.3,6f10.3,2g11.3,i4)                
                If(q1.lt.1.)Then
!!$OMP critical
!                   If( Xm(i)*Xscale*0.7> 58. .and.Xm(i)*Xscale*0.7<62.)
!     &  write(50,30)i,jp,Lab(jp),sqrt(dd)*Xscale*0.7,sqrt(dv),sqrt(VV),
!     &          Mhalo,Mhalo2,Vrmc(irr,i),Vrmc(irr,jp),Rh_1,Rh_2,
!     &          Xm(i)*Xscale*0.7,Ym(i)*Xscale*0.7,Zm(i)*Xscale*0.7  
!!$OMP end critical
                   Lab(i)=0     ! kill dublicate
                EndIf 
  !           EndIf                 ! end dv     check
          EndIf                    ! end dd<dR2 + Mhalo2>Mhalo/1.5 check
          EndIf                    ! end i=jp     check
          jp = Lst(jp)
          GoTo 10
          EndIf                    ! end jp neq 0
        EndDo                      ! end j loop
        EndDo                      ! end j loop
        EndDo                      ! end j loop
      EndIf                        ! end Lab(i)=-1 check
      EndDo                        ! end i loop

      m  =0
      meanmass =0
      write (13,*)  '  packing Halos in remVelDubl '
      !Lab = 1

      Do i=1,Ncentr
         If(Lab(i).eq.1)Then
            m = m +1
            Xm(m)      =Xm(i)
            Ym(m)      =Ym(i)
            Zm(m)      =Zm(i)
            Amc(m)     =Amc(i)   
            Rmc(m)     =Rmc(i)  
            RmaxV(m)   =RmaxV(i)
            Vrm(m)     =Vrm(i) 
            iRadius(m) =iRadius(i)

            Do ir=1,Nrad  !  Store Results 
               Wxc(ir,m)  =Wxc(ir,i) 
               Wyc(ir,m)  =Wyc(ir,i) 
               Wzc(ir,m)  =Wzc(ir,i) 
               Nbc(ir,m)  =Nbc(ir,i) 
!               NbcS(ir,m)  =NbcS(ir,i) 
!               NbcM(ir,m)  =NbcM(ir,i) 
!               NbcG(ir,m)  =NbcG(ir,i) 
               Rrc(ir,m)  =Rrc(ir,i)
               Vrmc(ir,m) =Vrmc(ir,i)
               Vrada(ir,m) =Vrada(ir,i)
            EndDo
            meanmass =meanmass +Nbc(Nrad,m)
         EndIf
      EndDo
      xmeanmass = FLOAT(meanmass)/FLOAT(m)      
      write(13,15) m,Ncentr,xmeanmass
      write(20,15) m,Ncentr,xmeanmass
      write(13,20) dRdubl*Xscale*hubble,fVdubl
 15   format(5x,'Vel.Dublicates: New number of objects=',i9,' Old=',i9,
     .       ' mean N_particles=',g11.3)
 20   format('     Distance Difference=',f8.4,/
     .       '     Vrelative/Vrms <    ',f8.4)
      Ncentr =m

      Return
      End
C--------------------------------------------------------------
C                                Find all neibhours for a center
C                                Xc,Yc,Zc - center; a0 -its weight
C
      SUBROUTINE Neib(Xnew,Ynew,Znew,Amnew,Xc,Yc,Zc,Radius)
C--------------------------------------------------------------
      use Structures
      Real*8  Xw,Yw,Zw,aMw,ww
C                                       Initiate counters
        Xnew =0.
        Ynew =0.
        Znew =0.
  	Amnew=0.
        Xw =0.
        Yw =0.
        Zw =0.
  	aMw=0.
        Radius2=Radius**2
        weig1 = 1.01*weightSmall
c                      limits for Label
         i1=INT((Xc-Radius)/Cell)
         j1=INT((Yc-Radius)/Cell)
         k1=INT((Zc-Radius)/Cell)
            i1=MIN(MAX(Nmx,i1),Nbx)
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
         i2=INT((Xc+Radius)/Cell)
         j2=INT((Yc+Radius)/Cell)
         k2=INT((Zc+Radius)/Cell)
            i2=MIN(MAX(Nmx,i2),Nbx)
            j2=MIN(MAX(Nmy,j2),Nby)
            k2=MIN(MAX(Nmz,k2),Nbz)
C                                        Look for neibhours
         nn =0
         Do k=k1,k2
         Do j=j1,j2
         Do i=i1,i2
            jp =Label(i,j,k)
 10         If(jp.ne.0)Then
               dd =(Xc-Xp(jp))**2+(Yc-Yp(jp))**2+(Zc-Zp(jp))**2
!               If(dd.lt.Radius2.and.iWeight(jp).le.weig1)Then
               If(dd.lt.Radius2)Then
                  nn =nn +1
                  ww =weightSmall  !iWeight(jp)
	             Xw =Xw +Xp(jp)*ww
	             Yw =Yw +Yp(jp)*ww
	             Zw =Zw +Zp(jp)*ww
	             aMw=aMw+ww
               EndIf  
                  jp =Lst(jp)
               GoTo10
            EndIf
         EndDo
         EndDo
         EndDo
         	If(aMw.GT.0.)Then
	            Xnew =Xw /aMw
	            Ynew =Yw /aMw
	            Znew =Zw /aMw
            EndIf
            Amnew = aMw
      Return
      End
C--------------------------------------------------------------
C                                Find all neibhours for a center
C                                Xc,Yc,Zc - center; 
C                               
      SUBROUTINE SNeib(Vx1,Vy1,Vz1,Nob1,Xc,Yc,Zc)
C--------------------------------------------------------------
      use Structures

      DIMENSION  Vx1(Nrad),Vy1(Nrad),Vz1(Nrad),Nob1(Nrad)
      Real*8     Wx1(Nrad),Wy1(Nrad),Wz1(Nrad)
      Integer*8 N, Ncp,i,jp,nn 

      iFlag = 0
         Radmax =Rad(Nrad)/FracSearch 
         d0      = Rad(1)
 50      dmax = Radmax**2
         Ninside = aNvir*Radmax**3


         Do i=1,Nrad
              Nob1(i) =0
              Wx1(i)  =0.
              Wy1(i)  =0.
              Wz1(i)  =0.
         EndDo
c                      limits for Label
         i1=INT((Xc-Radmax)/Cell)
         j1=INT((Yc-Radmax)/Cell)
         k1=INT((Zc-Radmax)/Cell)
            i1=MIN(MAX(Nmx,i1),Nbx)
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
         i2=INT((Xc+Radmax)/Cell)
         j2=INT((Yc+Radmax)/Cell)
         k2=INT((Zc+Radmax)/Cell)
            i2=MIN(MAX(Nmx,i2),Nbx)
            j2=MIN(MAX(Nmy,j2),Nby)
            k2=MIN(MAX(Nmz,k2),Nbz)
            nps =0
C                                        Look for neibhours
         Do k=k1,k2
         Do j=j1,j2
         Do i=i1,i2
            jp =Label(i,j,k)   !  Dark matter
 10         If(jp.ne.0)Then
               dd =(Xc-Xp(jp))**2+(Yc-Yp(jp))**2+(Zc-Zp(jp))**2
               If(dd.lt.dmax)Then
                  d1 =sqrt(dd)
                    ir =INT(max(0.,4.*(sqrt(d1/d0)-1.) +1)) +1
                    ir = min(max(1,ir),Nrad)
                    wu =weightSmall  !iWeight(jp)
                    iww = INT(wu/weightSmall+0.05)
                    Nob1(ir)= Nob1(ir)+ iww
c             write(30,'(20x,i9,i5,i4,6g11.3)')jp,iww,ir,wu,dd,Rad2(ir)
c                   If(iww.eq.iw1)NobS(ir)= NobS(ir)+ 1
c                    If(iww.eq.iw8)NobM(ir)= NobM(ir)+ 1
c                    If(iww.eq.iw64)NobG(ir)= NobG(ir)+ 1
c                    Rob1(ir)= Rob1(ir)+ d1*iww
c                    Vxx(ir) = Vxx(ir) +
c     +                            iww*(VXp(jp)**2+VYp(jp)**2+VZp(jp)**2)
                    Wx1(ir) = Wx1(ir) +iww*VXp(jp)  
                    Wy1(ir) = Wy1(ir) +iww*VYp(jp)
                    Wz1(ir) = Wz1(ir) +iww*VZp(jp)
               EndIf
               jp =Lst(jp)
               GoTo10
            EndIf
         EndDo
         EndDo
         EndDo
         Do ir =2,Nrad
            Nob1(ir)= Nob1(ir)   + Nob1(ir-1)
            Wx1(ir) = Wx1(ir)      + Wx1(ir-1)
            Wy1(ir) = Wy1(ir)      + Wy1(ir-1) 
            Wz1(ir) = Wz1(ir)      + Wz1(ir-1) 
         EndDo
 
         Do ir=1,Nrad
           Vx1(ir)  =0.
           Vy1(ir)  =0.
           Vz1(ir)  =0.
           If(Nob1(ir).GT.0)Then
               Nobj    = Nob1(ir)
               Vx1(ir) = Wx1(ir)/Nobj
               Vy1(ir) = Wy1(ir)/Nobj
               Vz1(ir) = Wz1(ir)/Nobj
           EndIf
         EndDo
!         If(Nob1(Nrad).gt.Ninside.and.iFlag.eq.0)Then
!             Radmax =Rad(Nrad)
!             Iflag = 1
!             goto50
!          EndIf 
      Return
      End
C-------------------------------------------------------------
C                     Find positions of 'Ncentr' spheres, which maximize number of particles
C                    inside 'radius Radius'. 
C                    Iterations are repeated until 1) mass stops increasing   OR
C                                                       2) the sphere stops moving: dR is less than Riter 
C                    Looks for all particles inside Radius for periodical set 
      SUBROUTINE Pairs(Ncentr,Riter)
C--------------------------------------------------------------
      Use Structures

      write (13,*) '  Rsearch =',RadSr1,RadSr2
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i)
      Do i=1,Ncentr
      	Lab(i) =1
         Amc(i) =0
         If(Rmc(i).lt.1.e-8)write (13,*) ' Error in Rmc: zero',Rmc(i),i
         If(Rmc(i).le.0.5*RadSr2)
     &   write (13,'(" Error Rmc too small=",g11.3," min=",g11.3,i8)')
     &              Rmc(i),RadSr2,i
         If(Rmc(i).gt.RadSr1)write (13,*) ' Error in Rmc: big ',Rmc(i),i
c         Rmc(i) =Radius
      EndDo
      Niter =0
      iDynamic = mDynamic
10    Iter  =0
      Niter =Niter +1
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,Xc,Yc,Zc,a0,in,Radius)
C$OMP+PRIVATE (Xnew,Ynew,Znew,Amnew,Dr)
C$OMP+SCHEDULE(DYNAMIC,iDynamic) REDUCTION(+:Iter)
      Do i=1,  Ncentr
         Xc =Xm(i)
         Yc =Ym(i)
         Zc =Zm(i)
         a0 =Amc(i)
         in =Lab(i)
         Radius =Rmc(i)
         If(Radius.lt.0.999*RadSR2)in =-Niter
         If(in.gt.0)Then         ! If not converged, keep iterating
            Call Neib(Xnew,Ynew,Znew,Amnew,Xc,Yc,Zc,Radius) ! new center of mass
            Dr =MAX(ABS(Xnew-Xc),ABS(Ynew-Yc),ABS(Znew-Zc))
            IF(Dr.LT.Riter.or.Amnew.lt.Amc(i))
     .                  Lab(i)=-Niter     ! iterations converged
            Iter =Iter +1
!	         IF(Xnew.lt.0.)Xnew =Xnew+Box
!	         IF(Ynew.lt.0.)Ynew =Ynew+Box
!	         IF(Znew.lt.0.)Znew =Znew+Box
!	         IF(Xnew.gt.Box)Xnew =Xnew-Box
!	         IF(Ynew.gt.Box)Ynew =Ynew-Box
!	         IF(Znew.gt.Box)Znew =Znew-Box
            Xm(i)  = Xnew
            Ym(i)  = Ynew
            Zm(i)  = Znew
            Amc(i) = Amnew
!            write (13,'(i6,6g12.5)') i,
!     &           Xm(i),Ym(i),Zm(i),Amc(i),Radius 
         EndIf
      EndDo

      Write (13,*) ' Iteration=',Niter,' Centers left=',Iter,iDynamic
      iDynamic =max(min(INT(Iter/float(Ncentr)*mDynamic),
     &                                     mDynamic),kDynamic)
      If(Iter.GE.1 .and. Niter.LT.25)Goto10

      Return
      End

C-------------------------------- Update Statistics of pairs 
C                   
      SUBROUTINE SPairs(Ncentr)
C---------------------------------------------------
      Use Structures

      Integer*8 N 
      DIMENSION  Vx1(Nrad),Vy1(Nrad),Vz1(Nrad),Vxx(Nrad),Nob1(Nrad),
     .             Rob1(Nrad),NobM(Nrad),NobG(Nrad),NobS(Nrad)

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,ir)
      Do i=1,Ncentr
            Do ir=1,Nrad
               Wxc(ir,i) = 0.
               Wyc(ir,i) = 0.
               Wzc(ir,i) = 0.
               Nbc(ir,i) = 0
!               NbcS(ir,i) = 0
!               NbcM(ir,i) = 0
!               NbcG(ir,i) = 0
               Rrc(ir,i) = 0.
               Vrmc(ir,i)= 0.
            EndDo
      EndDo
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,ir,Xc,Yc,Zc,Vx1,Vy1,Vz1,Nob1)
C$OMP+SCHEDULE (DYNAMIC,kDynamic)
       Do i=1,Ncentr
!          if(mod(i,10000).eq.0)write (13,*)  ' Sneib halo=',i
         Xc =Xm(i)
         Yc =Ym(i)
         Zc =Zm(i)
         Call SNeib(Vx1,Vy1,Vz1,Nob1,Xc,Yc,Zc)
         Do ir=1,Nrad  !  Store Results 
             Wxc(ir,i) =Vx1(ir)
             Wyc(ir,i) =Vy1(ir)
             Wzc(ir,i) =Vz1(ir)
             Nbc(ir,i) =Nob1(ir)
         EndDo
      EndDo

      Return
      End

C-------------------------------- Remove small halos
C           
C           
      SUBROUTINE SmallRem(ANlimit,Ncentr)
C--------------------------------------------------------------
      Use Structures

      New =0
      Do i=1,Ncentr
         If(Amc(i).GT.ANlimit*Rmc(i)**3)Then
          New =New +1
          Xm(New) =Xm(i)
          Ym(New) =Ym(i)
          Zm(New) =Zm(i)
          Amc(New)=Amc(i)
          Rmc(New)=Rmc(i)
         EndIf
      EndDo
      Ncentr =New
      write (13,*) ' small centers removed. New Ncentr=',Ncentr
      Return
      End
C-------------------------------- Get Catalog of Halos
C                        Remove close centers, keep only the largest
C               
      SUBROUTINE CATALOG(Ncentr,dRdubl)
C--------------------------------------------------------------
      Use Structures

      Real       Mhalo
      double precision q1
      Logical    inside

      dR2 =(max(Rad(1),dRdubl/dClose))**2
      Do i=1,Ncentr
        Lst(i)=-1
      EndDo
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i,j,k)
       Do k=Nmz,Nbz
       Do j=Nmy,Nby
       Do i=Nmx,Nbx
          Label(i,j,k)=0
       EndDo
       EndDo
       EndDo
  
      Do jp=1,Ncentr
         i=INT(Xm(jp)/Cell)
         j=INT(Ym(jp)/Cell)
         k=INT(Zm(jp)/Cell)
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
         Lst(jp)      =Label(i,j,k)
         Label(i,j,k) =jp
      EndDo
      write (13,'(a,i10,a,1P,g12.4)') ' CATALOG:',Ncentr,
     &           ' MinDist=',sqrt(dR2)*Xscale*hsmall
      nrem =0
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i,Mhalo,irad)
C$OMP+REDUCTION(+:nrem)
      Do i=1,Ncentr
        If(Amc(i).le.1.e-7)Then
          Lab(i) =-1
          nrem =nrem +1
        Else
          Lab(i) =1
        EndIf
      EndDo

      write (13,*) ' CATALOG:',Ncentr,' Removed Empty=',nrem

      nrem =0
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i,jp,Mhalo,dd,  
C$OMP+q1,Xc,Yc,Zc,i1,j1,k1,i2,j2,k2,ii,jj,kk,inside)
C$OMP+REDUCTION(+:nrem) SCHEDULE(DYNAMIC,500)
      Do i=1,Ncentr
c        If(i/10000*10000.eq.i)write (13,*)  ' RmClose i=',i
        If(Lab(i).ne.-1)Then         ! do not take empty centers
        Mhalo  =Amc(i)

        Xc =Xm(i)
        Yc =Ym(i)
        Zc =Zm(i)
c                      limits for Label
         i1=INT((Xc-dRdubl)/Cell)
         j1=INT((Yc-dRdubl)/Cell)
         k1=INT((Zc-dRdubl)/Cell)
            i1=MIN(MAX(Nmx,i1),Nbx)
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
         i2=INT((Xc+dRdubl)/Cell)
         j2=INT((Yc+dRdubl)/Cell)
         k2=INT((Zc+dRdubl)/Cell)
            i2=MIN(MAX(Nmx,i2),Nbx)
            j2=MIN(MAX(Nmy,j2),Nby)
            k2=MIN(MAX(Nmz,k2),Nbz)
C                                        Look for neibhours in nearest cells
         Do kk =k1,k2
         Do jj =j1,j2
         Do ii =i1,i2
            jp =Label(ii,jj,kk)   !  Dark matter Halo 
 10         If(jp.ne.0)Then

          If(i.ne.jp.and.Lab(jp).ne.-1)Then ! do not compare with empty centers
          dd =(Xm(jp)-Xc)**2 + (Ym(jp)-Yc)**2+ (Zm(jp)-Zc)**2     
          inside = (dd.lt.dR2)
           If(inside)Then    ! close distance pair
             q1   =Mhalo/Amc(jp)
             if(abs(q1-1.d0).lt.1.d-5.and.i.lt.jp)q1=0.9
             if(abs(q1-1.d0).lt.1.d-5.and.i.gt.jp)q1=1.1
             If(q1.lt.1.)Then
                Lab(i)=0     ! kill dublicate
                nrem  = nrem +1
             EndIf 
          EndIf                    ! end dd<dR2 
          EndIf                    ! end i=jp     check
          jp = Lst(jp)
          GoTo 10
          EndIf                    ! end jp neq 0
        EndDo                      ! end j loop
        EndDo                      ! end j loop
        EndDo                      ! end j loop
      EndIf                        ! end Lab(i)=-1 check
      EndDo                        ! end i loop

      write (13,*)' CATALOG: done removed halos=',nrem 
C                                 Remove maxima with Label =0
      New =0
      Do i=1,Ncentr
        If(Lab(i).eq.1)Then
          New =New +1
          Xm(New) =Xm(i)
          Ym(New) =Ym(i)
          Zm(New) =Zm(i)
          Amc(New)=Amc(i)
          Rmc(New)=Rmc(i)
      !      write (13,'(i6,3f9.3,g12.3,f8.2)')
      !&     New,Xm(i),Ym(i),Zm(i),Amc(i),Rmc(i)  
         EndIf 
      EndDo
      Ncentr =New
      write (13,*) ' dublicates removed. New Ncentr=',Ncentr

!      nrem =0
!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP+PRIVATE (i,jp,Mhalo,Xc,Yc,Zc,dd)  
!C$OMP+REDUCTION(+:nrem) 
!      Do i=1,Ncentr
!c        If(i/10000*10000.eq.i)write (13,*)  ' RmClose i=',i
!        Mhalo  =Amc(i)

!        Xc =Xm(i)
!        Yc =Ym(i)
!        Zc =Zm(i)
!      Do jp=1,Ncentr
!c        If(i/10000*10000.eq.i)write (13,*)  ' RmClose i=',i
!         if(jp.ne.i)Then
!          dd       = (Xm(jp)-Xc)**2 + (Ym(jp)-Yc)**2+ (Zm(jp)-Zc)**2     
!          If(dd.lt.dR2)Then
!             nrem =nrem+1
!             write(13,'(2i9,1P,13g12.4)')i,jp,sqrt(dd)*Xscale*hsmall,
!     &              Mhalo,Amc(jp),Xc,Xm(jp),Yc,Ym(jp),Zc,Zm(jp)
!          EndIf
            
!        EndIf
!        EndDo
!        EndDo



      Return
      End
C---------------------------------------------------------
C                              Select particles
c                              as initial seeds of centers
c                              Need to run it twice
c                      iFlag = 0 => count centers
c                            = 1 => get centers
      SUBROUTINE ReadCent(Ncentr,Boxx,iFlag,irank)
C--------------------------------------------------------------
      Use Structures
      Include 'mpif.h'

      Real INPUT,Iweigh_first,Iweigh_last,Iweigh_min
      Logical      Ldiff
      SAVE Nline,jneib

!      Xccc =49. ; Yccc =59. ; Zccc =106.
!      Rccc =2.
!      Xccc =INPUT(' Center x=')
!         write (13,*) Xccc
!      Yccc =INPUT(' Center y=')
!         write (13,*) Yccc
!      Zccc =INPUT(' Center z=')
!         write (13,*) Zccc
!      Rccc =INPUT(' Center r=')
!         write (13,*) Rccc
!         hsmall= hubble            ! Hubble/100.
!       Xccc =Xccc/Boxx*NGRID/hsmall
!      Yccc =Yccc/Boxx*NGRID/hsmall
!      Zccc =Zccc/Boxx*NGRID/hsmall
!      Rccc =Rccc/Boxx*NGRID/hsmall
!      write (13,'(a,4f8.2)')' Center for seeds: coords,Rad(Cells)=',
!     &       Xccc,Yccc,Zccc,Rccc 
!      write (13,'(a,f8.2,i5)')' Boxx,Ngrid=',Boxx,NGRID
      node = irank +1                                                                                                      
      Ncentr =0
      Iweigh_first =weightSmall
      Iweigh_last  =weightSmall  !iWeight(N)
      Iweigh_min  = weightSmall !Iweigh_first

        n1        = 0  
        nhalf    = 0
        n2        = 0
      If(RadSr1.lt.1.1*RadSr2.or.ABS(RadSr2).lt.1.e-8)Then
                                   ! there is only one radius
         Ldiff =.false.
         Radius = RadSr1
         R1        = RadSr1
         R2        = RadSr1
         Rhalf    = RadSr1
      Else
        Ldiff =.true.
        R1       = 1.1*RadSr2
        R2       = 0.9*RadSr1
        Rhalf   =(RadSr1+RadSr2)/2.
        ialpha =2
        RRs  = (RadSr1-RadSr2)
      Endif 

      If(iFlag.eq.0)Then
         If(irank ==0)Then
      Nline=INPUT(' Enter Number of particles for initial seeds=> ')
         write (13,*)
      jneib=INPUT(' Enter Number of neibhours for a seed       => ')
         write (13,*)
         EndIf
      CALL MPI_BCAST(Nline,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(jneib,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      EndIf

      If(Nline.eq.0)Goto 50
         fr    = Np/FLOAT(Nline)
         xrand =1.00001
         Nstep =MAX(INT(fr+1.),1)
         fr    =1./fr
         Nseed = 121071
         write (13,*) ' Every ',Nstep,' particle is used',
     &                ' as initial center',' Fraction=',fr
         Nplast = lspecies(1)
         !Do i=1,N        ! take centers if weight is less than Iweigh_min
         !  If( iWeight(i).le.Iweigh_min)Then
         !     Nplast =Nplast +1
         !  Endif 
         !Enddo 
         write (13,*) ' Total pool for particles centers=',Nplast
          Do i=1,Nplast
             If(fr.lt.0.999)xrand =RANDd(Nseed) 
               IF(xrand.lt.fr)Then 
!                 IF(ABS(Xp(i)-Xccc).lt.Rccc.and. 
!     .                ABS(Yp(i)-Yccc).lt.Rccc.and.
!     .                ABS(Zp(i)-Zccc).lt.Rccc)THEN
                  Ncentr           =Ncentr+1
                  IF(Ldiff)Then           ! different radii
                        Rrrr   =RadSr2+
     &                            RRs*(MAX(RANDd(Nseed),1.e-8))**ialpha
                        If(iFlag.eq.1) Rmc(Ncentr) = Rrrr

                        If(Rrrr.lt.R1)n1 = n1+1
                        If(Rrrr.lt.R2.and.Rrrr.ge.R1)nhalf = nhalf+1
                        If(Rrrr.gt.R2)n2 = n2+1
                  Else
                        If(iFlag.eq.1) Rmc(Ncentr) =Radius
                  Endif 
                        If(iFlag.eq.1)Then 
                           Amc(Ncentr) =1.
                           Xm(Ncentr)   =Xp(i)
                           Ym(Ncentr)   =Yp(i)
                           Zm(Ncentr)   =Zp(i)
                        EndIf
                  j=Ncentr

!                 EndIf      ! endif for small box test
               EndIf
         EndDo
         write (13,*) ' Ncenters (random fraction)=',Ncentr
      If(jneib.eq.0)goto 150               ! No more points
 50   Ncnt_max = 0         ! select a fraction of particles from
                           ! each cell
      Ncnt_lrg = 0
      Ncnt_pmax = 0
      Ncnt_p        = 0
      Ncnt_p1      = 0
!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP+PRIVATE (i,j,k,jp,Ncnt,jstep,jcount,Ncount) 
      Do k=Nmz,Nbz
         If(mod(k,10)==0)Write(13,*) '       k=',k,' Ncenter=',Ncentr
      Do j=Nmy,Nby
      Do i=Nmx,Nbx
         jp =Label(i,j,k)
         If(jp.ne.0)Then
           Ncnt  =0
           Do while (jp.ne.0)  ! count number of particles in the cell
              Ncnt = Ncnt+1
              jp=Lst(jp) 
           End Do
             !    step for selecting centers  
           jstep = jneib*min(100,max(INT(Ncnt/500.),1))
             !
           jp        = Label(i,j,k)
           jcount    = 0
           Ncount    = 0
 10        jcount    = jcount+1
              If(jp.ne.0.and.jcount.lt.jstep)Then
                     jp=Lst(jp) ! find if there are enough neib
                     goto 10
              endif 
               jcount    =0
           If(jp.ne.0)Then 
!           If(iWeight(jp).le.Iweigh_min)Then !  take this particle 
              Ncount =Ncount +1
!                IF(ABS(Xp(jp)-Xccc).lt.Rccc.and. 
!     .                 ABS(Yp(jp)-Yccc).lt.Rccc.and.
!     .                 ABS(Zp(jp)-Zccc).lt.Rccc)THEN
              If(iFlag == 0)Then
!$OMP atomic
                Ncentr=Ncentr+1
              Else
!$OMP critical
                Ncentr=Ncentr+1
                Xm(Ncentr)  =Xp(jp)
                Ym(Ncentr)  =Yp(jp)
                Zm(Ncentr)  =Zp(jp)
!$OMP end critical
                 EndIf     ! iFlag ==0

!                  EndIf           ! endif for small-box test
c                 If(Ncount.lt.300) Goto 10
                 Goto 10
!           EndIf    ! end (iWeight< min) condition
           EndIf    ! end (jp>0)         condition
           EndIf    ! end (jp>0)         condition
           If(Ncount.gt.0)Then
             Ncnt_max = max(Ncnt_max,Ncnt)
             Ncnt_pmax = max(Ncnt_pmax,Ncount)
             Ncnt_p =  Ncnt_p +Ncount
             Ncnt_p1 =  Ncnt_p1 +1
             If(Ncnt.gt.1e4)Ncnt_lrg = Ncnt_lrg+1
           EndIf 
      EndDo
      EndDo
      EndDo
      If(iFlag ==1)Then
         Do i=1,Ncentr
                  IF(Ldiff)Then           ! different radii
                        Rrrr   =RadSr2+
     &                            RRs*(MAX(RANDd(Nseed),1.e-8))**ialpha
                        Rmc(i) = Rrrr
                  Else
                        Rmc(i) =Radius
                  Endif 
                           Amc(i) =1.
         EndDo
      EndIf
150   write(13,*) ' Ncentr (total)=',Ncentr  !,' Max allowed=',Nc
      write(13,*) ' Max number of particles in a cell   =',Ncnt_max
      write(13,*) ' N cells with more than 1e4 particles=',Ncnt_lrg
      write(13,*) ' Max number of selected    particles =',Ncnt_pmax
      write(13,*) ' Average number of selected particles=',
     &                                Ncnt_p/max(Ncnt_p1,1)
      !If(Ncentr.gt.Nc)Then
      !   write (13,*) ' Too many centers'
      !   STOP ' Too many centers'
      !Endif 
         If(Ncentr.gt.0.and.iFlag.eq.1)Then
            nempty = 0
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE (i) REDUCTION(+:n1,n2,nhalf,nempty) 
            Do i=1,Ncentr
                   If(Rmc(i).lt.RadSr2)nempty = nempty+1
                   If(Rmc(i).lt.R1)n1 = n1+1
                   If(Rmc(i).lt.R2.and.Rmc(i).ge.R1)nhalf = nhalf+1
                   If(Rmc(i).gt.R2)n2 = n2+1
             EndDo
             write(13,*) '     Empty centers =',nempty
            write (13,200) RadSR2,R1,n1,R1,R2,nhalf,R2,RadSR1,n2
            If(node==1)
     &        write (2,200) RadSR2,R1,n1,R1,R2,nhalf,R2,RadSR1,n2
            !write (20,200) RadSR2,R1,n1,R1,R2,nhalf,R2,RadSR1,n2

 200   format(' Number of centers with radii between ',3x,2f7.4,
     &                ' was ',i8,
     &              /18x,' with radii between ',3x,2f7.4,' was ',i8,
     &              /18x,' with radii between ',3x,2f7.4,' was ',i8)
         Endif   
      Return
      End

C---------------------------------------------------
C                                Remove Unbound particles, find radius and mass
C                                No removal of unbound particles If:
C                                       Toohot =< 0 or Toohot =>5
C                             Vscale =100.*hsmall*Xscale/AEXPN 
C                              Xscale=Box/Ngrid 
C                              Dscale=2.75e+11*hsmall**2*(Box/NROW)**3 
      SUBROUTINE CLEAN(Ncentr,Radsc,MnMass,Toohot,Rminhalo)  
C--------------------------------------------------------------
      Use Structures

      dimension Mass(Nrad),MnMass(Nrad),Overdens(Nrad),Radsc(Nrad),
     &          Vrad(Nrad)
      Real          Mass,MnMass,Mhalo

C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,ir,iter,Xc,Yc,Zc,Summass,Mass,Overdens)
C$OMP+PRIVATE (irad,Mhalo,Rhalo,Vmaxrot,Rmaxrot,Vcirc,Vcc)
C$OMP+PRIVATE (NobG,Rob1,Vxx,ifirst,iflag,Vrad )
C$OMP+SCHEDULE (DYNAMIC,kDynamic)
      Do i=1,Ncentr          
         if(mod(i,50000).eq.0)write(13,*)' Clean halo=',i,Ncentr
         Xc =Xm(i)  *Xscale !          Scale all to physical units: comoving
         Yc =Ym(i)  *Xscale
         Zc =Zm(i)  *Xscale
         iflag =0
         ifirst=0
         Do ir =1,Nrad
            Summass     = Om0*Nbc(ir,i)
            If(Nbc(ir,i).gt.nFirst.and.iflag.eq.0)Then
               ifirst = ir            ! find first significant bin
               iflag  = 1
            EndIf 
            Mass(ir)    = Dscale*Summass
            Overdens(ir)= Mass(ir)/MnMass(ir)
            Vrad(ir) =0.
         EndDo  

C                                     find size, mass, rad ... for a halo
         Call Decide(Mass,Radsc,Overdens,Nrad,irad,Mhalo,Rhalo,
     +               Vrad,Vmaxrot,Rmaxrot,ifirst,dRmaxrot )
         Vcirc=sqrt(Vmaxrot/AEXPN)*6.58e-5 ! =10^-5*sqrt(G 2e33/3.08e24) 

         If(Toohot.gt.0. .and. Toohot.lt.5.)Then
         Do iter =0,3      !  Remove High Energy particles from halos
c            If(i==9724)write (35,'(i7,3f8.4,2i5,4g12.4)')
c     &           i,Xm(i),Ym(i),Zm(i),Nrad,irad,Vcirc,Mhalo,Rhalo 
c         If(irad .lt.1 .or. Rhalo.le.Rminhalo)Then
c           iic = iic +1
c            write (30,'(i7,i4,2g11.3,3f9.3,3x,10i6)')
c     &         i,irad,Rhalo,Mhalo,Xc,Yc,Zc,(Nbc(ir,i),ir=1,10) 
c         else 
c            write (32,'(i7,2i4,2g11.3,3f9.3,3x,10i6)')
c     &         i,irad,iter,Rhalo,Mhalo,Xc,Yc,Zc,(Nbc(ir,i),ir=1,10) 
c         EndIf 
         If(irad .ge.1 .and. Rhalo.gt.Rminhalo)Then
           Vcc =Vcirc/Vscale          ! back to dimensionless units
           Vmaxrot=Vcc*Toohot !*1.03**(3-iter) ! Increase the Rotatioanl Vel.

           Call RemEng(i,irad,Vmaxrot,Rmaxrot/Xscale)  ! Remove Unbound Particles
           iflag =0
           ifirst=0
           Do ir =1,Nrad
            If(Nbc(ir,i).gt.nFirst.and.iflag.eq.0)Then
               ifirst = ir            ! find first significant bin
               iflag  = 1
            EndIf 
            If(Nbc(ir,i).gt.mFirst.and.ir.gt.2)then
               Vrad(ir)    = abs(Vrada(ir,i))*Vscale/
     &           (sqrt(Mass(ir)/Rad(ir)/Xscale/AEXPN)*6.58e-5)
            else
               Vrad(ir) =0.
            EndIf
              Vrmc(ir,i)=Vrmc(ir,i) *Vscale
             Summass     =Om0*Nbc(ir,i)
             Mass(ir)    = Dscale*Summass
             Overdens(ir)=Mass(ir)/MnMass(ir)
           EndDo  

          Call Decide(Mass,Radsc,Overdens,Nrad,irad,Mhalo,Rhalo,
     +                Vrad, Vmaxrot,Rmaxrot,ifirst,dRmaxrot)
          Vcirc=sqrt(Vmaxrot/AEXPN)*6.58e-5 !=10^-5*sqrt(2e33/3.08e24)

         EndIf
         EndDo       ! end iter energy
         EndIf
         Amc(i)     =Mhalo ! store results
         Rmc(i)     = Rhalo 
         RmaxV(i)   = Rmaxrot
         Vrm(i)     =Vcirc
         iRadius(i) =irad

       Enddo   ! end i (Ncentr)
       write(13,*) ' Clean:',Ncentr,Xscale
!      Do i=1,Ncentr          
!         Xc =Xm(i)  *Xscale*0.7 !          Scale all to physical units: comoving
!         Yc =Ym(i)  *Xscale*0.7
!         Zc =Zm(i)  *Xscale*0.7
!             write (51,'(2i6,3f8.3,5g11.3)')i,iRadius(i),
!     &     Xm(i)*xscale,Ym(i)*xscale,Zm(i)*xscale,
!     &     Amc(i), Rmc(i), Vrm(i)
!      EndDo 

       Return
       End
C------------------------------------------------------------- 
C                                       Find  halo radius and mass
      SUBROUTINE Decide(Mass,Rad,Overdens,Nrad,irad,Mhalo,Rhalo,
     +                  Vrad,Vmaxrot,Rmaxrot,ifirst,dRmaxrotA )
C--------------------------------------------------------------
      COMMON /OVER/ Ovlim,Ovcnt
      Real   Rad(Nrad),Overdens(Nrad),Mass(Nrad),Mhalo,
     &       Vrad(Nrad)
      DATA  slope/0.4999/
         Vmaxrot =0
         Rmaxrot =0. 
         Mhalo=0.
         Rhalo=0.   
      If(ifirst.eq.0)Then
         irad =0         ! too few particles
         RETURN
      EndIf

      If(MAX(Overdens(ifirst),Overdens(2)).LT.Ovcnt)Then
         irad =0         ! max overdensity is less then Ovcnt
         RETURN
      EndIf
      iov  =0
      Do ir=Nrad,ifirst,-1
         If(Overdens(ir).gt.Ovlim .and. iov.eq.0)Then
            iov  =ir
         EndIf
      EndDo
      If(iov.eq.0)Then
         irad =0         ! max overdensity is less then Ovlim
         RETURN
      EndIf
      iovv =0
      Do ir=ifirst,iov
         irleft =max(ir-1,1)
         irr1   =min(ir+1,Nrad)
         irr2   =min(ir+2,Nrad)
         Overright =(2.*Overdens(irr1)+Overdens(irr2))/3.
         Overleft  =(2.*Overdens(ir)+Overdens(irleft))/3.
         If(Overright.gt.0.97*Overleft.or.Vrad(ir).gt.0.5)Then !stop:rising profile
!         If(Overright.gt.1.2*Overleft.or.Vrad(ir).gt.0.5)Then !stop:rising profile
            If(ir.eq.1)Then ! softer conditions for the first bin
               If(Overright.gt.1.2*Overleft)Then
                   iovv =ir
                   goto 20
               EndIf
            else 
                iovv =ir
               goto 20
            EndIf
         EndIf
      EndDo
 20   If(iovv.ne.0)iov =iovv
         irad =iov
            If(irad.lt.ifirst)write (13,*)' err: irad=',irad,ifirst


 30      irleft  =max(irad-2,1)  ! check if overdens is steep
         irright =min(irad+1,Nrad)
         gradlf  =Overdens(irleft) * (Rad(irleft))**slope
         gradrt  =Overdens(irright)* (Rad(irright))**slope
         If(gradlf.lt.gradrt)Then  ! not a steep gradient
            irad =irad -1        ! reduce radius
            iovv =irad
            If(irad.lt.ifirst)Then   ! zero radius -> quit
               irad =0
               RETURN
            Else
               goto 30             ! iterate
            EndIf
         EndIf
c            If(irad.lt.ifirst)write (13,*)' err2: irad=',irad,ifirst
         If(irad.lt.Nrad)Then    ! interpolate to Overlim
            If(Ovlim.gt.Overdens(irad+1) .and.
     .         Overdens(irad).gt.Ovlim)Then       
              Mhalo=(Mass(irad)*(Ovlim-Overdens(irad+1))+
     .             Mass(irad+1)*(Overdens(irad)-Ovlim)) /
     .             (Overdens(irad)-Overdens(irad+1))
              Rhalo=(Rad(irad)*(Ovlim-Overdens(irad+1))+
     .             Rad(irad+1)*(Overdens(irad)-Ovlim)) /
     .             (Overdens(irad)-Overdens(irad+1))
            Else                 ! Overdens is > Overlim 
              Mhalo =Mass(irad)
              Rhalo =Rad(irad)
            EndIf
         Else                 ! Overdens is > Overlim 
            Mhalo =Mass(irad)
            Rhalo =Rad(irad)
         EndIf

         imax  = 0
         Do ir =1,irad  ! find max rotational velocity
            V2    = Mass(ir)/Rad(ir)
            If (V2.gt.Vmaxrot)Then
               Vmaxrot =V2               
               Rmaxrot =Rad(ir)
               imax = ir
            EndIf
         EndDo
                                ! improve Rmaxrot value by finding
                              ! radius containing F(1)/F(2.15)=0.415 fraction
                              ! of mass at Vmax. Then Rmax=2.15*R
         imax = max(imax,1)
         aMrs = 0.415*Mass(imax)
c         write (13,*)  ' == imax=',imax,' mass=',aMrs
c         write (13,'(1P,10G11.3)') Mass
         Do ir=1,irad
            If(Mass(ir).gt.aMrs)Then
               if(ir.eq.1)Then
                  aM0 =0.
                  r0     =0.
               Else
                  aM0 =Mass(ir-1)
                  r0     =Rad(ir-1)
               EndIf
               Rmaxrot2 =2.15*(r0+(Rad(ir)-r0)*(aMrs-aM0)/
     &                          max(Mass(ir)-aM0,1.e-3))
               Rmaxrot =min(Rmaxrot,Rmaxrot2,Rhalo/dRmaxrotA)

               return
            EndIf 
          EndDo 
 50       write (13,'(4i4,1P,10G11.3)') ir,irad,ifirst,iovv,
     &             Mass(ir),Ovlim,   Rmaxrot,
     &                Rmaxrot2,Rad(ir),r0
          write (13,'(1P,10g11.3)')Overdens  

      Return
      End
C------------------------------------------------
C                                  Remove unbound particles for i-th halo
C                                  irad   = radius of the halo
C                                 Vmaxrot = max rotational Velocity
C                                 Rmaxrot = radius at max rot.velocity
      SUBROUTINE RemEng(i,irad,Vmaxrot,Rmaxrot)
C-----------------------------------------------
      Use Structures

      DIMENSION VXxc(Nrad), VYyc(Nrad), VZzc(Nrad),
     &             Vescape2(Nrad)            
      DIMENSION  Vx1(Nrad),Vy1(Nrad),Vz1(Nrad),Vxx(Nrad),Nob1(Nrad),
     .      Rob1(Nrad),NobM(Nrad),NobG(Nrad),NobS(Nrad),Varad(Nrad)

         Xc  =Xm(i)  ! store center of mass and velocity      
         Yc  =Ym(i)     
         Zc  =Zm(i)     
            irr = 1  !Max(irad,2)         ! don't take V from the first bin
 10         if(Nbc(irr,i).lt.nmin.and.irr.lt.iCentr)Then
               irr =irr+1
               goto 10
            EndIf 
            VXc =Wxc(irr,i)
            VYc =Wyc(irr,i)
            VZc =Wzc(irr,i)

         Do ir =1,Nrad ! store mean velocity
            VXxc(ir) =Wxc(ir,i) 
            VYyc(ir) =Wyc(ir,i) 
            VZzc(ir) =Wzc(ir,i) 
                                         ! escape velocity for NFW profile
            Vescape2(ir) = (2.15*Vmaxrot)**2 * 
     &           max(log(1.0+2.*Rad(ir)/Rmaxrot)/(Rad(ir)/Rmaxrot),
     &                   0.25) 
               Wxc(ir,i) = 0.
               Wyc(ir,i) = 0.
               Wzc(ir,i) = 0.
               Nbc(ir,i) = 0
!               NbcS(ir,i) = 0
!               NbcM(ir,i) = 0
!               NbcG(ir,i) = 0
               Rrc(ir,i) = 0.
               Vrmc(ir,i)= 0.
               Vrada(ir,i)= 0.
         EndDo
            Call SrejNeib(VXxc,VYyc,VZzc,Xc,Yc,Zc,Vescape2,
     &           Vx1,Vy1,Vz1,Nob1,NobS,NobM,NobG,Rob1,Vxx
     &           ,VXc,VYc,VZc,Varad)
            Do ir=1,Nrad  !  Store Results 
               Wxc(ir,i) =Vx1(ir)
               Wyc(ir,i) =Vy1(ir)
               Wzc(ir,i) =Vz1(ir)
               Nbc(ir,i) =Nob1(ir)
!               NbcS(ir,i) =NobS(ir)
!               NbcM(ir,i) =NobM(ir)
!               NbcG(ir,i) =NobG(ir)
               Rrc(ir,i)  =Rob1(ir)
               Vrmc(ir,i) =Vxx(ir)
               Vrada(ir,i)=Varad(ir)
            EndDo
      Return
      End
C--------------------------------------------------------------
C                                Accumilate statistics only
C                                if relative velocity is not large
C                                Xc,Yc,Zc - center; 
C                               
      SUBROUTINE SrejNeib(VXxc,VYyc,VZzc,Xc,Yc,Zc,Vescape2,
     &          Vx1,Vy1,Vz1,Nob1,NobS,NobM,NobG,Rob1,Vxx
     &           ,VXc,VYc,VZc,Varad)
C--------------------------------------------------------------
      Use Structures

      Integer*8 jp
      DIMENSION VXxc(Nrad), VYyc(Nrad), VZzc(Nrad),Vescape2(Nrad)
      DIMENSION  Vx1(Nrad),Vy1(Nrad),Vz1(Nrad),Vxx(Nrad),Nob1(Nrad),
     .       Rob1(Nrad),NobM(Nrad),NobG(Nrad),NobS(Nrad),Varad(Nrad)

      iFlag = 0
         Radmax =Rad(Nrad)/FracSearch
         d0      = Rad(1)

 50      dmax = Radmax**2
         Ninside = aNvir*Radmax**3
         iw1 = 1
         iw8 = 8
         iw64 = 64

         Do i=1,Nrad
              Nob1(i) =0
!              NobS(i) =0
!              NobM(i) =0
!              NobG(i) =0
              Rob1(i) =0.
              Vx1(i)  =0.
              Vy1(i)  =0.
              Vz1(i)  =0.
              Vxx(i)  =0.
!              Varad(i) =0.
         EndDo
c                                           limits for Label
         i1=INT((Xc-Radmax)/Cell)
         j1=INT((Yc-Radmax)/Cell)
         k1=INT((Zc-Radmax)/Cell)
            i1=MIN(MAX(Nmx,i1),Nbx)
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
         i2=INT((Xc+Radmax)/Cell)
         j2=INT((Yc+Radmax)/Cell)
         k2=INT((Zc+Radmax)/Cell)
            i2=MIN(MAX(Nmx,i2),Nbx)
            j2=MIN(MAX(Nmy,j2),Nby)
            k2=MIN(MAX(Nmz,k2),Nbz)
C                                        Look for neibhours
         Do k=k1,k2
         Do j=j1,j2
         Do i=i1,i2
           jp =Label(i,j,k)       !  Dark matter
             Do while (jp.ne.0)
                dd =(Xc-Xp(jp))**2+(Yc-Yp(jp))**2+(Zc-Zp(jp))**2
               If(dd.lt.dmax)Then
                  d1 =sqrt(dd)
                    iov =INT(   max(0.,4.*(sqrt(d1/d0)-1.) +1)) +1
                    iov = min(max(1,iov),Nrad)
c                    dv =  (VXxc(iov)-VXp(jp))**2+
c     +                     (VYyc(iov)-VYp(jp))**2+
c     +                     (VZzc(iov)-VZp(jp))**2
                    dv =  (VXc-VXp(jp))**2+
     +                    (VYc-VYp(jp))**2+
     +                    (VZc-VZp(jp))**2
                If(dv .lt. Vescape2(iov))Then
                    ww =weightSmall   !iWeight(jp)
                    iww =INT(ww/weightSmall+0.005)
                    Nob1(iov)= Nob1(iov)+ iww
!                    If(iww.eq.iw1)NobS(iov)= NobS(iov)+ 1
!                    If(iww.eq.iw8)NobM(iov)= NobM(iov)+ 1
!                    If(iww.eq.iw64)NobG(iov)= NobG(iov)+ 1
                    Rob1(iov)= Rob1(iov)+ d1*iww
                    Vxx(iov) = Vxx(iov) +
     +                          iww*(VXp(jp)**2+VYp(jp)**2+VZp(jp)**2)
                    Vx1(iov) = Vx1(iov) +iww*VXp(jp)  
                    Vy1(iov) = Vy1(iov) +iww*VYp(jp)
                    Vz1(iov) = Vz1(iov) +iww*VZp(jp)
                    Varad(iov) =Varad(iov) + iww*
     +                   ( (VXp(jp)-VXc)*(Xp(jp)-Xc)+
     +                     (VYp(jp)-VYc)*(Yp(jp)-Yc)+
     +                     (VZp(jp)-VZc)*(Zp(jp)-Zc) )
     &                   /max(d1,1.e-8)
c               write (13,'(i7,i3,2i5,g11.4,2x,3g11.3,2x,6g11.4)') 
c     &            jp,iov,Nob1(iov),iww,Varad(ir),
c     &             Xp(jp),Yp(jp),Zp(jp),VXp(jp),VYp(jp),VZp(jp),
c     &              VXc,VYc,VZc
                EndIf     !  dv
             EndIf        ! dd
           jp =Lst(jp)
            EndDo    !  jp 
         EndDo
         EndDo
         EndDo

         Do ir =1,Nrad
     	      Varad(ir) = Varad(ir)/max(Nob1(ir),1)
         EndDo
         Do ir =2,Nrad-1
            If(Nob1(ir).eq.0)Varad(ir)=(Varad(ir-1)+Varad(ir+1))/2.
         EndDo

         Do ir =2,Nrad  ! make intergal quantities
            Nob1(ir)= Nob1(ir)   + Nob1(ir-1)
!            NobS(ir)= NobS(ir)   + NobS(ir-1)
!            NobM(ir)= NobM(ir)   + NobM(ir-1)
!            NobG(ir)= NobG(ir)   + NobG(ir-1)
            Rob1(ir)= Rob1(ir)   + Rob1(ir-1)
            Vxx(ir) = Vxx(ir)    + Vxx(ir-1)
            Vx1(ir) = Vx1(ir)    + Vx1(ir-1)
            Vy1(ir) = Vy1(ir)    + Vy1(ir-1) 
            Vz1(ir) = Vz1(ir)    + Vz1(ir-1) 
         EndDo

        Do ir=1,Nrad
           Vmn =0.
     	      If(Nob1(ir).GT.0)Then
               Nobj =Nob1(ir)
               Vx1(ir) = Vx1(ir)/Nobj
               Vy1(ir) = Vy1(ir)/Nobj
               Vz1(ir) = Vz1(ir)/Nobj
              Vmn     = Vx1(ir)**2+Vy1(ir)**2+Vz1(ir)**2
               Vxx(ir) = SQRT(ABS(Vxx(ir)/Nobj -Vmn))
	         EndIf
         EndDo
         If(Nob1(Nrad).gt.Ninside.and.iFlag.eq.0)Then
             Radmax =Rad(Nrad)
             Iflag = 1
             goto50
          EndIf 

      Return
      End
C---------------------------------------------------
C                     update memory usage      
C                     stop if too much memory requested
      SUBROUTINE Memory(Nelements)
C---------------------------------------------------
      USE STRUCTURES
      PARAMETER (Gb    =1024.**3)
      include 'mpif.h'

      TotalMemory =TotalMemory + Nelements*4./Gb
          write(13,*)' Memory =',TotalMemory
          If(TotalMemory>MaxMem)Then
              CALL mpi_finalize(ierr)                                                                                            
               Stop ' not enough memory'
          EndIf

      End SUBROUTINE Memory
C---------------------------------------------------
C                                  Read  current data from disk/tape,
C                                  Open files
C                                  Nrecl is the number of values in a record
C                                  Npage is the number of particles in a page
      SUBROUTINE RDdata
C---------------------------------------------------

      USE STRUCTURES

      Character ::   fname*120
      Open(4,file ='PMcrd.DAT',form ='UNFORMATTED',status ='UNKNOWN')

      READ  (4,err=20,end=20) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                   ,Ocurv,extras
      WRITE (13,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv
!      WRITE (16,100) HEADER,
!     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
!     +                  EKIN,EKIN1,EKIN2,
!     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
!     +                  Ocurv
100   FORMAT(1X,'Header=>',A45,/
     +          1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I9,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3)
      IF(NROW.NE.NROWC) THEN
         WRITE (13,*)
     +            ' NROW in PARAMETER and in TAPE-FILE are different'
         write (13,*)  ' NROW,NGRID (PMparamters.h)=',NROW,NGRID
         write (13,*)  ' NROW,NGRID (PMcrd.DAT)    =',NROWC,NGRIDC
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (13,*)
     +           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (13,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      ENDIF
      write(13,*) ' Number of particles =',lspecies(Nspecies)

      CLOSE (4)
      NACCES= NRECL*4 / nbyteword 
      Jfiles =1
      If(lspecies(Nspecies)> 1024**3)Jfiles=8
         Do i =1,Jfiles
            write(fname,'(a,i1,a)')'PMcrs',i-1,'.DAT'
            OPEN(UNIT=20+i,FILE=TRIM(fname),ACCESS='DIRECT', 
     &	               FORM='unformatted',STATUS='UNKNOWN',RECL=NACCES)
         EndDo

      RETURN
 20   write (13,*) ' Error reading the header file: Abort'
      stop
      END SUBROUTINE RDdata

C---------------------------------- Read in variables      
      REAL FUNCTION INPUT(text)
C------------------------------------------------
      Character text*(*)
          write (13,'(A,$)')text
          read (*,*) x
          INPUT =x
      Return
      End  FUNCTION INPUT
c     ------------------------
      real function seconds ()
c     ------------------------
c
c     purpose: returns elapsed time in seconds
      Integer, SAVE :: first=0,rate=0
      Real   , SAVE :: t0=0.
      !real*8 dummy
      !real tarray(2)
      If(first==0)Then
         CALL SYSTEM_CLOCK(i,rate)
         first =1
         t0    =float(i)/float(rate)
         seconds = 0.
      Else
!        seconds = etime(tarray)
         CALL SYSTEM_CLOCK(i)
         seconds = float(i)/float(rate)-t0
      EndIf

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

