!-------------------------------------------------------------------------------- 
!                  Read Files
!                  
!                  (nx,ny,nz) - number of boxes in each direction
!                  Rhalo  = width of buffer in Mpch arounf each region
Module SetArrs

  Integer, PARAMETER    ::                         &
                         NROW        = 4096,             &
                         NGRID         = 256,               &
                         Nrecord       = 20e6,             &   ! max number of particles in buffer
                         nbyteword  = 4,                    &
                         NPAGE       = NROW**2,     & ! # particles in a record
                         NRECL        = NPAGE*6     
  Integer*8,PARAMETER ::Np=2048_8**3  ! max number of particles 
  Real ::           AEXPN,ASTEP,PARTW,  & 
                        Om0,Oml0,hubble
  Integer ::      Nspecies,Nseed,ISTEP,Nx,Ny,Nz
  CHARACTER*45           HEADER

  Real,          DIMENSION(Nrecord) :: Xp,Yp,Zp,VXp,VYp,Vzp
  Integer*8, DIMENSION(Nrecord) :: id_part
  Character*80 :: fNames

  Real                ::  Box,Rad,xc,yc,zc,xxL,xxR,yyL,yyR,zzL,zzR,Xscale, &
                              Vscale
  Integer*8       ::  iCount
  Integer, Allocatable           :: ListFiles(:,:,:)

Contains
!-------------------------------------------------------------------------------- 
!                          Define configuration of files
!                          Read control infor from the first file
!                          allocate arrays
Subroutine ReadConfig
  character*80 :: Name

  write(*,'(a,$)')' Enter snapshot number ='
  read(*,*)jStep
                         ! read header information from the first file
  write(Name,'(a,i4.4,a)') 'PMss.0001.',jStep,'.DAT'
  open(1,file=Name,form='unformatted')

        Read(1)HEADER
        Read(1)AEXPN,ASTEP,ISTEP,NROWC,NGRIDC,Nspecies,  &
                                 Nseed,Om0,Oml0,hubble,Box
        Read(1) k,Nx,Ny,Nz,dR
  write(*,'(a)')HEADER
  write(*,'(a,f9.4,5x,a,2i5)')' aexpn=',AEXPN,' Step=',ISTEP,jStep
  if(ISTEP .ne. jStep)Stop ' -- wrong time step --'
  Allocate (ListFiles(Nx,Ny,Nz))
  write(*,'(a,3i5)')' Configuration of files: in each direction=',Nx,Ny,Nz
  write(*,'(a,$)')' Enter list (triplets) of nodes to read (zero for all files) ='
  Xscale= Box/Ngrid        ! Scale for comoving coordinates
  Vscale= Box*100./Ngrid/AEXPN        ! Scale for comoving coordinates
  ListFiles = 0
  Nlist =0

  Do                                 ! read list of files
     read(*,*,iostat=ierr) kx, ky, kz
     if(ierr /= 0)exit
     If(kx == 0)Then
        Nlist = 0
        exit
     EndIf
     if(kx<1.or.kx>Nx)Stop ' wrong node in X direction'
     if(ky<1.or.ky>Ny)Stop ' wrong node in Y direction'
     if(kz<1.or.kz>Nz)Stop ' wrong node in Z direction'
     Nlist =Nlist +1
     ListFiles(kx,ky,kz) = 1
  EndDo  ! end reading file list
  If(Nlist == 0)Then    !  if list is empty, take all files
     Nlist = Nx*Ny*Nz
     ListFiles = 1
  EndIf
  write (*,*) ' Number of files to read=',Nlist
  Close (1)

  Do i=1,Nx        ! open files
     Do j=1,Ny
        Do k=1,Nz
           node = i+(j-1)*Nx +(k-1)*Nx*Ny
           If(ListFiles(i,j,k) ==1)Then
                  write(Name,'(a,i4.4,a,i4.4,a)') 'PMss.',node,'.',jStep,'.DAT'
                  open(100+node,file=Name,form='unformatted')
           EndIf
        EndDo ! k
     EndDo ! j
  EndDo ! i
  write(Name,'(a,i4.4,a)') 'Part.',jStep,'.DAT'
  open(1,file=Name,form='formatted')
end Subroutine ReadConfig
end Module SetArrs
!-------------------------------------------------------------------------------- 
!
Program ReadFiles
  use SetArrs
  Logical :: inside
 
  Call ReadConfig                     ! set parameters: list of nodes to read
                                                    !                     names of files; open files
   write(*,'(/a,$)')' Enter Radius and center ='
   read(*,*)Rad,xc,yc,zc
   xxL = (xc-Rad)/Xscale+1. ; xxR = (xc+Rad)/Xscale +1.  ! limits in internal units
   yyL = (yc-Rad)/Xscale+1. ; yyR = (yc+Rad)/Xscale +1.
   zzL = (zc-Rad)/Xscale+1. ; zzR = (zc+Rad)/Xscale +1.
   xxL =-10.
   xxR = 300.
   yyL =-10.
   yyR = 300.

   write(*,'(/a,$)')' Enter fraction of particles to write ='
   read(*,*)fraction
   iSeed = 102191 ! seed for random numbers
  Do i=1,Nx        ! read files, select and dump particles
     Do j=1,Ny
        Do k=1,Nz
           node = i+(j-1)*Nx +(k-1)*Nx*Ny
           If(ListFiles(i,j,k) == 1)Then
                  write(*,'(a,3x,4i4)')' file=',i,j,k,node
                  Call ReadHeader(node,xL,xR,yL,yR,zl,zR)
                  inside = (min(xR,xxR)>max(xL,xxL))
                  inside = (min(yR,yyR)>max(yL,yyL)).and.inside
                  inside = (min(zR,zzR)>max(zL,zzL)).and.inside
                  If(inside)Call ReadParticles(node,iSeed,fraction)
                  close(100+node)
           EndIf
        EndDo ! k
     EndDo ! j
  EndDo ! i

  Stop
end Program ReadFiles
!---------------------------------------------------------
!                                        Read particles
      SUBROUTINE ReadHeader(node,xL,xR,yL,yR,zl,zR)
!--------------------------------------------------------------
       use SetArrs
       write(*,*) '   Reading file=',node
        Read(100+node)HEADER
        Read(100+node)AEXPN,ASTEP,ISTEP,NROWC,NGRIDC,Nspecies,  &
                                 Nseed,Om0,Oml0,hubble,Box
        Read(100+node) k,Nx,Ny,Nz,dR
        Read(100+node) xL,xR,yL,yR,zl,zR
        write(*,'(a,3x,6f8.2)') '          limits=',xL,xR,yL,yR,zl,zR
end SUBROUTINE ReadHeader
!---------------------------------------------------------
!                                        Read particles
      SUBROUTINE ReadParticles(node,iSeed,fraction)
!--------------------------------------------------------------
       use SetArrs
       Logical :: inside,takeit

       takeit = .true.
   read(100+node)nn        ! total number of particles in this file
        write(*,*) '     Number of particles in this node   =',nn
   Do 
      read(100+node,iostat=ierr)mm         ! number of particles in this record        
      If(ierr /= 0) exit
        write(*,*) '     Number of particles in this record =',mm
      read(100+node)(Xp(i),Yp(i),Zp(i),              &
                                    VXp(i),VYp(i),VZp(i),      &
                                    id_part(i),i=1,mm)
      Do i=1,mm
         inside = (xxL <= Xp(i)).and. (Xp(i)<xxR) .and.        &
                        (yyL <= Yp(i)).and. (Yp(i)<yyR) .and.      &
                        (zzL <= Zp(i)).and. (Zp(i)<zzR)
         If(inside)Then
            If(fraction<0.999)takeit =(Randd(iSeed)<fraction)
            If(takeit) &
            write(1,'(i10,3f10.4,1P,3g12.4)')                        &
                      id_part(i),(Xp(i)-1.)*Xscale,(Yp(i)-1.)*Xscale,(Zp(i)-1.)*Xscale, &
                                       VXp(i)*Vscale,VYp(i)*Vscale,VZp(i)*Vscale
        endif ! inside
      EndDo
   EndDo

end SUBROUTINE ReadParticles
!------------------------------------------------
!				                                       random number generator
      FUNCTION RANDd(M)
!------------------------------------------------
      DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4  &
      	/453815927,2147483648.,2147483647,536870912,131072,256, &
     	 16777216,4,16384,8388608,128/
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
      END FUNCTION Randd


