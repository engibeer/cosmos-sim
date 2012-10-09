!-------------------------------------------------------------------------------- 
!                  Split simulation box in smaller chunks 
!                  for BDM
!                  (nx,ny,nz) - number of boxes in each direction
!                  Rhalo  = width of buffer in Mpch arounf each region
Module SetArrs

Integer, PARAMETER    ::                          &
                         NROW        = 4096,            &
                         NGRID         = 256,              &
                         nx                 = 5,                  &  !number of regions in X
                         ny                 = 5,                  &
                         nz                 = 5,                  &
                         Nfiles            = nx*ny*nz,    &   !  number of files
                         Nrecord       = 170e6,         &   ! max number of particles in buffer
                         Lchunk         = 20e6,           &   ! max number of particles per record
                         nbyteword  = 4,                   &
                         NPAGE       = NROW**2,    & ! # particles in a record
                         NRECL        = NPAGE*6     
      Integer*8,PARAMETER ::Np=2048_8**3  ! max number of particles 
  COMMON / ROW /	        XPAR(NPAGE),YPAR(NPAGE),ZPAR(NPAGE), &
     			VX(NPAGE),VY(NPAGE),VZ(NPAGE)
  COMMON / CONTROL/      AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,  & 
                        TINTG,EK,EK1,EK2,AU0,AEU0,     &
                        NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5, &
                        Ocurv,extras(100)
  COMMON / HEADDR/       HEADER
  CHARACTER*45           HEADER

 Real*4,     PARAMETER        ::  Rhalo=2.2 
 Real,         DIMENSION(Np)  :: Xp,Yp,Zp,VXp,VYp,Vzp
 Integer*8, DIMENSION(Nrecord,Nfiles) :: id_part
 Character*80 :: fNames
 Real*4 :: bndry(2,3,Nfiles)

  Real                ::  Box,RECDAT(NRECL),wspecies(10)
  Integer*8       ::  lspecies(10)
  Integer*8       ::  jCount(Nfiles),iCount(Nfiles)


  EQUIVALENCE           (wspecies(1),extras(1)), &
                       (lspecies(1),extras(11))
!$OMP THREADPRIVATE(/ROW/)
Contains
!--------------------------------------------------------------------------------------------
!
!                            Open Files, read control information
!
Subroutine OpenPM
      Character ::   fname*120
      Open(4,file ='PMcrd.DAT',form ='UNFORMATTED',status ='UNKNOWN')

      READ  (4,err=20,end=20) HEADER, &
                      AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0, &
                       NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5 &
                        ,Ocurv,extras
      WRITE (*,100) HEADER, &
                       AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       EKIN,EKIN1,EKIN2, &
                       NROWC,NGRID,NRECL,Om0,Oml0,hubble, &
                       Ocurv
100   FORMAT(1X,'Header=>',A45,/ &
               1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/ &
                 1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/ &
                 1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I9,/ &
                 1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/ &
                 1x,' Omega_curvature=',F7.3)
      IF(NROW.NE.NROWC) THEN
         WRITE (*,*) &
                 ' NROW in PARAMETER and in TAPE-FILE are different'
         write (*,*)  ' NROW,NGRID (PMparamters.h)=',NROW,NGRID
         write (*,*)  ' NROW,NGRID (PMcrd.DAT)    =',NROWC,NGRIDC
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (*,*) &
                ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      ENDIF
      write(*,*) ' Number of particles =',lspecies(Nspecies)
      Box = extras(100)

      CLOSE (4)
      NACCES= NRECL*4 / nbyteword 
      Jfiles =1
      If(lspecies(Nspecies)> 1024**3)Jfiles=8
         Do i =1,Jfiles
            write(fname,'(a,i1,a)')'PMcrs',i-1,'.DAT'
            OPEN(UNIT=20+i,FILE=TRIM(fname),ACCESS='DIRECT',  &
     	               FORM='unformatted',STATUS='UNKNOWN',RECL=NACCES)
         EndDo

      RETURN
 20   write (*,*) ' Error reading the header file: Abort'
      stop
end Subroutine OpenPM

end Module SetArrs
!-------------------------------------------------------------------------------- 
!
Program SplitBox
  use SetArrs
  Character :: Line*80
  integer*8 :: Lsmall,N_particles,JPAGE,iL,IN,N
 
  dR  =Rhalo ! boundary width in Mpch units
  Call OpenPM
           Xscale= Box/Ngrid        ! Scale for comoving coordinates
           dR = dR/Xscale
  Call ReadPntP(N)
  Call Split(N,dR)
  Call WriteFiles


  Stop
end Program SplitBox
!---------------------------------------------------------
!                                        Read particles
      SUBROUTINE WriteFiles
!--------------------------------------------------------------
       use SetArrs

       dB = NGRID
!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE (node,i,j,nn,mm,Nchunks,xL,yL,zL,xS,yS,zS,i0,i1) 
Do node=1,Nfiles/2
   nn = iCount(node)
   Nchunks =(nn-1)/Lchunk+1 
    xL = bndry(1,1,node)   ! boundaries of the box
    yL = bndry(1,2,node) 
    zL = bndry(1,3,node) 
    xS = xL -NGRID   ! boundaries of the box
    yS= yL -NGRID 
    zS = zL -NGRID 
    !write(*,'(a,3f9.3)')' Bounadries=',xL,yL,zL
   write(100+node)nn
   Do j=1,Nchunks
      mm = Lchunk
      if(j==Nchunks)mm =nn-(Nchunks-1)*Lchunk
      i0= 1+(j-1)*Lchunk    !  start
      i1= i0+mm-1               ! finish
      write(100+node)mm
      write(100+node)(mod(Xp(id_part(i,node))-xS,dB)+xL,  &
                                     mod(Yp(id_part(i,node))-yS,dB)+yL,  &
                                     mod(Zp(id_part(i,node))-zS,dB)+zL  , &
                                    VXp(id_part(i,node)),VYp(id_part(i,node)),VZp(id_part(i,node)), &
                                id_part(i,node),i=i0,i1)

      EndDo
EndDo   ! node

!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE (node,i,j,nn,mm,Nchunks,xL,yL,zL,xS,yS,zS,i0,i1) 
Do node=Nfiles/2+1, Nfiles
   nn = iCount(node)
   Nchunks =(nn-1)/Lchunk+1 
    xL = bndry(1,1,node)   ! boundaries of the box
    yL = bndry(1,2,node) 
    zL = bndry(1,3,node) 
    xS = xL -NGRID   ! boundaries of the box
    yS= yL -NGRID 
    zS = zL -NGRID 
    !write(*,'(a,3f9.3)')' Bounadries=',xL,yL,zL
   write(100+node)nn
   Do j=1,Nchunks
      mm = Lchunk
      if(j==Nchunks)mm =nn-(Nchunks-1)*Lchunk
      i0= 1+(j-1)*Lchunk    !  start
      i1= i0+mm-1               ! finish
      write(100+node)mm
      write(100+node)(mod(Xp(id_part(i,node))-xS,dB)+xL,  &
                                     mod(Yp(id_part(i,node))-yS,dB)+yL,  &
                                     mod(Zp(id_part(i,node))-zS,dB)+zL  , &
                                    VXp(id_part(i,node)),VYp(id_part(i,node)),VZp(id_part(i,node)), &
                                id_part(i,node),i=i0,i1)

      EndDo
EndDo   ! node


100 format(/3f10.3,i10)

      end SUBROUTINE WriteFiles
!---------------------------------------------------------
!                                        Read particles
      SUBROUTINE Split(N,dR)
!--------------------------------------------------------------
       use SetArrs

      Integer*8 :: N,i,N_particles,Lsmall,Jpage,iL,IN 
      Real*8 :: xmin,xmax,ymin,ymax,zmin,zmax
      Logical :: inx,iny,inz

      dx = float(NGRID)/float(nx)
      dy = float(NGRID)/float(ny)
      dz = float(NGRID)/float(nz)
      xbig = 1.+NGRID
      xzero =1.
                                     ! define boundaries of regions
      kk =0

      Do k=1,nz
            bz1 = 1.-dR+ (k-1.)*dz
            bz2 = 1.+dR+ k        *dz
      Do j=1,ny
            by1 = 1.-dR+ (j-1.)*dy
            by2 = 1.+dR+ j        *dy
      Do i=1,nx
           bx1 = 1.-dR+ (i-1.)*dx
           bx2 = 1.+dR+ i     *dx
            kk =kk +1
            bndry(1,1,kk) =bx1
            bndry(2,1,kk) =bx2
            bndry(1,2,kk) =by1
            bndry(2,2,kk) =by2
            bndry(1,3,kk) =bz1
            bndry(2,3,kk) =bz2
         EndDo
       EndDo
     EndDo
     If(kk/= Nfiles)Stop ' Error in Nfiles'
     write(*,*) ' File         x             y             z'
     Do k=1,Nfiles
        write(*,'(i3,5x,3(2x,2f8.2))')k,((bndry(i,j,k),i=1,2),j=1,3)
      EndDo
                                      ! Open files and write control information
     Do k=1,Nfiles
        write(fNames,'(a,i4.4,a,i4.4,a)') 'PMss.',k,'.',ISTEP,'.DAT'
        open(100+k,file=fNames,form='unformatted') 
        write(100+k)HEADER
        write(100+k)AEXPN,ASTEP,ISTEP,NROWC,NGRIDC,Nspecies,  &
                                 Nseed,Om0,Oml0,hubble,Box
        write(100+k) k,Nx,Ny,Nz,dR
        write(100+k) ((bndry(i,j,k),i=1,2),j=1,3)
!        write(100+k,*)HEADER
!        write(100+k,100)AEXPN,ASTEP,ISTEP,NROWC,NGRIDC,Nspecies,  &
!                                 Nseed,Om0,Oml0,hubble,Box
!        write(100+k,'(4i5)') k,Nx,Ny,Nz
!        write(100+k,'(6f8.3)') ((bndry(i,j,k),i=1,2),j=1,3)
       EndDo
100 format(f9.4,f9.5,i6,2i6,i4,i10,4f8.3)


      iCount = 0
!                              Find particles for each region.
!                              Once buffer is full, dump to a file
!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE (node,i,xx,yy,zz,xL,xR,yL,yR,zL,zR,dx,dy,dz,inx,iny,inz) 
Do node=1,Nfiles
    iCount(node) =0
    xL = bndry(1,1,node)    ! boundaries of the box
    xR = bndry(2,1,node)
    yL = bndry(1,2,node)
    yR = bndry(2,2,node)
    zL = bndry(1,3,node)
    zR = bndry(2,3,node)

    dx = 0.                            ! need for periodical conditions
    if(xL<1.0)dx =-NGRID
    if(xR>xbig)dx =NGRID
    dy = 0.
    if(yL<1.0)dy =-NGRID
    if(yR>xbig)dy =NGRID
    dz = 0.
    if(zL<1.0)dz =-NGRID
    if(zR>xbig)dz =NGRID
    Do i=1,N
       xx = Xp(i)
       inx = (xx>=xL.and.xx<xR)
       If(.not.inx.and.abs(dx)>1.)Then
          xx = xx + dx
          inx = (xx>=xL.and.xx<xR)
       EndIf
       If(inx)Then
          yy = Yp(i)
          iny = (yy>=yL.and.yy<yR)
          If(.not.iny.and.abs(dy)>1.)Then
             yy = yy + dy
             iny = (yy>=yL.and.yy<yR)
          EndIf
          If(iny)Then
             zz = Zp(i)
             inz = (zz>=zL.and.zz<zR)
             If(.not.inz.and.abs(dz)>1.)Then
                zz = zz + dz
                inz = (zz>=zL.and.zz<zR)
             EndIf
             If(inz)Then
                iCount(node) = iCount(node) +1
               If( iCount(node)<Nrecord)  id_part(iCount(node),node) = i
             EndIf    ! inz
          EndIf     ! iny
       EndIf     ! inx
     EndDo    ! i
   EndDo      ! node
       write(*,*) ' End Split. Number of particles on nodes='
       write(*,'(5i10)') (iCount(i),i=1,Nfiles)
       Do node=1,Nfiles
                If(iCount(node)>Nrecord)Then
                   write(*,'(a,i10)') ' Too many particles in node =',node, &
                                    ' Increase Nrecord=',Nrecord, &
                                    ' Found = ',iCount(node)
                   Stop
                EndIf
       EndDo
end SUBROUTINE Split

!---------------------------------------------------------
!                                        Read particles
      SUBROUTINE ReadPntP(N)    
!--------------------------------------------------------------
       use SetArrs

      Integer*8 :: N,i,N_particles,Lsmall,Jpage,iL,IN 
      Real*8 :: xmin,xmax,ymin,ymax,zmin,zmax,xbig


      ! lspecies(1) = lspecies(1)/2_8   !!!
      ! write (*,*) ' Temp: N=',lspecies(1) 
      N =0
      Jpage  = NPAGE

      xmin =1.e10
      xmax =-1.e10
      ymin =1.e10
      ymax =-1.e10
      zmin =1.e10
      zmax =-1.e10
      xzero = 1.
      xbig    = float(NGRID)+1.d0-1.d-4
      write (*,*) ' X small/big=',xzero,xbig
      Nsmall =0
      Lsmall =lspecies(1)

         N_particles =lspecies(Nspecies)   ! Total number of particles
         Npages = (N_particles -1)/Jpage +1
         N_in_last=N_particles -Jpage*(Npages-1)

      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) ' N_in_last=',N_in_last
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP  PRIVATE (IROW,In_page,iL,IN,j) 
      DO  IROW=1, Npages         ! Loop over particle pages
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
            iL = JPAGE*(IROW-1)
            if(mod(IROW,10).eq.0)write (*,*)' Read page=',IROW,' Np=',iL
!				       Loop over particles
          If(lspecies(Nspecies)>1024**3)Then
            j = (IROW-1)/64 +1
            if(j>8.or.j<1)write(*,*) ' Error in file number: ',j
            read(20+j,REC=mod(IROW-1,64)+1)xpar,ypar,zpar,vx,vy,vz        
          Else   
            READ  (21,REC=IROW) XPAR,YPAR,ZPAR,VX,VY,VZ
          EndIf
	   DO   IN=1,In_page
                Xp(iL+IN)  =XPAR(IN)
                Yp(iL+IN)  =YPAR(IN)
                Zp(iL+IN)  =ZPAR(IN)
                VXp(iL+IN)  =VX(IN)
                VYp(iL+IN)  =VY(IN)
                VZp(iL+IN)  =VZ(IN)
           EndDo      ! IN
      EndDo           ! IROW
      N      = lspecies(Nspecies)
      write(*,*) ' Done Reading Points. N=',N
 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)                 &
!$OMP  REDUCTION(MAX:xmax,ymax,zmax) &
!$OMP REDUCTION(MIN:xmin,ymin,zmin) 
      Do i   =1,N 
         If(Xp(i).gt.xbig)Xp(i)=Xp(i)-1.e-4
         If(Yp(i).gt.xbig)Yp(i)=Yp(i)-1.e-4
         If(Zp(i).gt.xbig)Zp(i)=Zp(i)-1.e-4
            If(i.le.Lsmall)Then
               xmin = min(Xp(i),xmin)
               xmax =max(Xp(i),xmax)
               ymin = min(Yp(i),ymin)
               ymax =max(Yp(i),ymax)
               zmin = min(Zp(i),zmin)
               zmax =max(Zp(i),zmax)
            EndIf  
      ENDDO 


      write(*,*) ' Small Particles:',Lsmall
      write(*,*) '               x:',xmin,xmax
      write(*,*) '               y:',ymin,ymax
      write(*,*) '               z:',zmin,zmax
      write (*,*)  ' Nmaxpart=',Nmaxpart,Np
      !do i=1,1000
      !write (*,'(i10," Coords=",6g13.6)')i,Xp(i),Yp(i),Zp(i),
      !&                                             VXp(i),VYp(i),VZp(i)
      !EndDo 

      Return
      End  SUBROUTINE ReadPntP
!-------------------------------------------------------------------- 
      SUBROUTINE Node_to_IJK(node,i2,j2,k2)
!                             give {ijk} for node
!-------------------------------------------------------------------- 
        use SetArrs

         k2 =INT((node-1)/(nx*ny))+1
         jj = node-(k2-1)*nx*ny
         j2 =INT((jj-1)/nx) +1
         i2 = jj-(j2-1)*nx
      Return
      End SUBROUTINE Node_to_IJK

