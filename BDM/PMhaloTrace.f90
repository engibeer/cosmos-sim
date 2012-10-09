!-------------------------------------------------
!    Trace Major progenitors of halos
!     
MODULE   setHalos
PUBLIC      

Integer, parameter  :: Ncats =174,  &        ! total number of  snapshots to analyze
                       Nsnaps=6,       &       !  number of buffered snapshots
                       Nmax  = 13e6,     &       ! max #  halos in a snapshot
                       Lmax   =  70e6, &       ! max number of halos in merging tree
                       Nmx      = 0,         &      ! dimensions of linker list
                       Nbx       = 400,     &
                       Nmy      = 0,         &      ! dimensions of linker list
                       Nby       = 400,     &
                       Nmz      = 0,         &      ! dimensions of linker list
                       Nbz       = 400,     &
                       Nbnd  = 50                  ! max #  of bound particles particles  
Real, parameter :: Box=250.,   &
                                 Omega0 = 0.27, &  ! Omega_matter
                                 hubble = 0.70, &
                                 Rmax   = 0.5,      &     !  Max distance to search for progenotors in linker-list
                                 Wmin   = 2.0,    &       ! Minimum weight for accepting matches of halos
                                 dmax   = 0.5,    &       ! max distance for search of progenitors
                                 d2max  = dmax**2
Integer, DIMENSION(Ncats)             :: Nhalos 
Real,    DIMENSION(Ncats)             :: Aexpn 
Integer :: iGood(Ncats,Lmax),     &             ! pointers to halos in the same merging brunch
                                                                          ! iGood(:,k) = progenitors of halo 'k'
                 iGlobal(Nsnaps,Nmax),  &           ! iGlobal(i,k) is  pointer of halo 'k' at current moment to
                                                                          !                  a position in 'iGood' array at snapshot i
                 jNotFound(Nmax),        &             ! indexes of halos, which are not detected at next moment
                 Lst(Nmax),                     &              ! part of linker-list
                 Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz), &    
                 jTree(Lmax),                           &        !  pointer to position in Tree(:,:)
                  Nselected,                          &         ! number of selected halos
                  Last                                                  ! number of halos in the list of good matches
Logical    :: Selected(Lmax)                           ! List of selected halos for final Tree dump
Real*4     :: wGood(Ncats,Lmax)                  ! 'Goodness' of matching halos. Maximum =4.
Real*4     :: Cell,dCorrect
Character*10, DIMENSION(Ncats)  :: iNames

TYPE :: haloData
  Integer    :: num                                                         ! halo number
  Real         ::   x,y,z,vx,vy,vz,aM,Rm,rmsV,circV    ! data for current halo
  Integer    :: sel                                                            ! number of selected particles
  Integer*8    :: bound(Nbnd)                                      ! list of bound particles
End TYPE haloData

TYPE :: haloRecord
  Integer    :: nh                                                              ! halo number in catalog 
  Integer    :: Iso                                                             ! index of isolation: 0 - distinct, /=0  pointer for parent
  Real         :: x,y,z,vx,vy,vz,aM,Rm,circV
End TYPE haloRecord

TYPE :: haloSnap
  TYPE(haloData) :: halo(Nmax)     ! data for current halo
End TYPE haloSnap

Type(haloSnap) ::         Snapshot(Nsnaps) 
Type(haloRecord) :: Catalog(Nmax),Tree(Ncats-Nsnaps,Lmax)
Contains
!---------------------------------------------------------
!                         Read 'iC' halos catalog including particles IDs
!                         Store in  buffer 'iCput'
      SUBROUTINE ReadHalos(iC,iCput)
!--------------------------------------------------------------
Character*120   :: listname, Line, Txt*5, Txt2*50
Integer*8            :: nbd(Nbnd)
Type(haloData) :: halo2
                              !                                                                                 Open file
      listname='HaloLists/HaloListA.'//TRIM(iNames(iC))//'.DAT'
      Open( iC+100,file=listname,Status='UNKNOWN',form='unformatted')! short list of halos
   read(iC+100) aa
           !write(*,'(a5,f8.3,a)')' A=',aa
   Aexpn(iC) = aa
           !write(*,*)' iC=',iC,Aexpn(iC)
   Read(iC+100)Nmax_s,Ncenter,Rout,dLog,Risolate,Unbound
           !write(*,110) Nmax_s,Ncenter,Rout,dLog,Risolate,Unbound
 110  format(1x,'Max number of selected particles=',i5, &
                  /1x,'Number of halos=',i12,&
                  /1x,' Outer radius in virial radius= ',f6.2,' dLog =',f6.3,  &
                           ' Isolation radius (Rvir units)=',f6.2, &
                         /1x,'Unbound factor (Vcirc units)= ',f6.2)
    
      Do kh =1,Nmax
           !read(iC+100,*,iostat=iStat) i
           !read(iC+100,*,iostat=iStat)j,halo1,nSelect
           !if(mod(kh,50000)==0)write(*,*) ' Reading ',iC,kh
        read(iC+100,iostat=iStat)j,x,y,z,vx,vy,vz,aM,Rm,rmsV,circV,nSelect
        If(iStat.ne.0)Exit
        If(kh /= j) Then
           write(*,*) ' Error: iC=',iC,' kh =',kh,' halo=',j,' iStat=',iStat
           Exit
        EndIf
        If(nSelect /= 0)Then
        backspace (iC+100)
        read(iC+100,iostat=iStat)j,x,y,z,vx,vy,vz,aM,Rm,rmsV,circV,nSelect,(nbd(k),k=1,nSelect)
                halo2%num =j
                halo2%x = x ; halo2%y = y ; halo2%z = z 
                halo2%vx = vx ; halo2%vy = vy ; halo2%vz = vz 
                halo2%aM = aM ; halo2%Rm = Rm 
                halo2%rmsV = rmsV ; halo2%circV = circV
                halo2%sel = nSelect
                Do k =1,nSelect
                    halo2%bound(k) =nbd(k)
                EndDo
          !write(*,'(2i8,3x,i6,10g12.4,3x,i3)')iC,i,Snapshot(iC)%halo(i)%num, &
          !                         Snapshot(iC)%halo(i), &
          !                         Snapshot(iC)%halo(i)%sel
        Else
                halo2%num =j
                halo2%x = x ; halo2%y = y ; halo2%z = z 
                halo2%vx = vx ; halo2%vy = vy ; halo2%vz = vz 
                halo2%aM = aM ; halo2%Rm = Rm 
                halo2%rmsV = rmsV ; halo2%circV = circV
                halo2%sel = nSelect
                halo2%bound =0
        EndIf
        Snapshot(iCput)%halo(j) = halo2

      EndDo
10    Nhalos(iC) =Min(j,Nmax)

      write(*,'(a,i8,a,f7.4,a,i3)') ' Read Halos=',Nhalos(iC),' a=',Aexpn(iC),' Buffer=',iCput
      close(100+iC)
      Return

      End SUBROUTINE ReadHalos
!---------------------------------------------------------
!                         Read halo catalog without id's
!                         Store in 'Catalog'
      SUBROUTINE LoadSnap(iC)
!--------------------------------------------------------------
Character*120   :: listname, Line, Txt*5, Txt2*50
Integer*8            :: nbd(Nbnd)
Type(haloData) :: halo2
                              !                                                                                 Open file
      listname='HaloLists/HaloListA.'//TRIM(iNames(iC))//'.DAT'
      Open( iC+100,file=listname,Status='UNKNOWN',form='unformatted')! short list of halos
   read(iC+100) aa
           !write(*,'(a5,f8.3,a)')' A=',aa
   Aexpn(iC) = aa
           !write(*,*)' iC=',iC,Aexpn(iC)
   Read(iC+100)Nmax_s,Ncenter,Rout,dLog,Risolate,Unbound
           !write(*,110) Nmax_s,Ncenter,Rout,dLog,Risolate,Unbound
 110  format(1x,'Max number of selected particles=',i5, &
                  /1x,'Number of halos=',i12,&
                  /1x,' Outer radius in virial radius= ',f6.2,' dLog =',f6.3,  &
                           ' Isolation radius (Rvir units)=',f6.2, &
                         /1x,'Unbound factor (Vcirc units)= ',f6.2)
    
      Do kh =1,Nmax
        read(iC+100,iostat=iStat)j,x,y,z,vx,vy,vz,aM,Rm,rmsV,circV
        If(iStat.ne.0)Exit
                Catalog(j)%nh =j
                Catalog(j)%x   = x    ; Catalog(j)%y   = y   ; Catalog(j)%z   = z 
                Catalog(j)%vx = vx ; Catalog(j)%vy = vy ; Catalog(j)%vz = vz 
                Catalog(j)%aM = aM ; Catalog(j)%Rm = Rm 
                Catalog(j)%circV = circV
      EndDo
10    Nhalos(iC) =Min(j,Nmax)

      Catalog(:)%Iso = -1
      write(*,'(a,i8,a,f7.4,a,i3)') ' Read Halos=',Nhalos(iC),' a=',Aexpn(iC)
      close(100+iC)
      Return

      End SUBROUTINE LoadSnap
!--------------------------------------------------------------
!                                Make list of  distinct halos
!    Input:            Ncenter - numer of halos, {Xc ...} = halo data     
!                           Cell = size of linker-list cell
!                           Radmax = maximum halo radius    
!    Output:         jNotFound(i) = 0 for distinct halo
!                          
      SUBROUTINE HaloDistinct(Ncenter,iNext,iFlag)
!--------------------------------------------------------------
Real, parameter :: Radmax = 2.2   ! maximum virial radius

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
        Do i=1,Ncenter
           jNotFound(i) =0
        EndDo

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ic,i,j,k,x,y,z,R0,i1,j1,k1,i2,j2,k2,jp, & 
!$OMP                 dd,i3,j3,k3,aM0,xp,yp,zp,R0p,aM0p,dx,dy,dz)
        Do ic=1,Ncenter
           If(iFlag ==0) Then
             x = Snapshot(iNext)%halo(ic)%x
             y = Snapshot(iNext)%halo(ic)%y
             z = Snapshot(iNext)%halo(ic)%z
             R0 =Snapshot(iNext)%halo(ic)%Rm/1000.
             aM0 =Snapshot(iNext)%halo(ic)%aM
           Else
             x = Catalog(ic)%x
             y = Catalog(ic)%y
             z = Catalog(ic)%z
             R0 =Catalog(ic)%Rm/1000.
             aM0 =Catalog(ic)%aM
           EndIf

           i1=INT((x-Radmax)/Cell+100)-100
           j1=INT((y-Radmax)/Cell+100)-100
           k1=INT((z-Radmax)/Cell+100)-100
            i1=MIN(MAX(Nmx,i1),Nbx)
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
           i2=INT((x+Radmax)/Cell)
           j2=INT((y+Radmax)/Cell)
           k2=INT((z+Radmax)/Cell)
           i2=MIN(MAX(Nmx,i2),Nbx)
           j2=MIN(MAX(Nmy,j2),Nby)
           k2=MIN(MAX(Nmz,k2),Nbz)
!                                        Look for neibhours
         Do k3=k1,k2
            k =mod(k3+50*(Nbz+1),Nbz+1)  ! periodical boundaries
         Do j3=j1,j2
            j =mod(j3+50*(Nby+1),Nby+1)  ! periodical boundaries
         Do i3=i1,i2
            i =mod(i3+50*(Nbz+1),Nbz+1)  ! periodical boundaries
            jp =Label(i,j,k)   
           Do while (jp.ne.0)
           If(iFlag ==0) Then
             xp = Snapshot(iNext)%halo(jp)%x
             yp = Snapshot(iNext)%halo(jp)%y
             zp = Snapshot(iNext)%halo(jp)%z
             R0p =Snapshot(iNext)%halo(jp)%Rm/1000.
             aM0p =Snapshot(iNext)%halo(jp)%aM
           Else
             xp = Catalog(jp)%x
             yp = Catalog(jp)%y
             zp = Catalog(jp)%z
             R0p =Catalog(jp)%Rm/1000.
             aM0p =Catalog(jp)%aM
           EndIf
               If(aM0p.ge.aM0)Then
                dx =  xp - x ; dy =  yp - y ; dz =  zp - z 

             dx = min(dx,abs(dx-Box),abs(dx+Box))                   ! use periodical boundaries
             dy = min(abs(dy),abs(dy+Box),abs(dy-Box)) 
             dz = min(abs(dz),abs(dz+Box),abs(dz-Box)) 
                 dd =(x-xp)**2+(y-yp)**2+(z-zp)**2
                  If(dd< (Max(R0p,R0))**2) Then
                       If(aM0p.gt.aM0.or.jp>ic)jNotFound(ic) =jp  ! not distinct
                   EndIf        ! dd<Rmc
               EndIf          ! Rmc => R0
               jp =Lst(jp)
           End Do ! jp /=0
         EndDo   ! i
         EndDo   ! j
         EndDo   ! k

         EndDo  ! ic

         nn=0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i) REDUCTION(+:nn)
        Do i=1,Ncenter
           If(jNotFound(i).eq.0)nn=nn+1
        EndDo
        write(*,*) ' Number of distinct halos =',nn
      Return
      End SUBROUTINE HaloDistinct

!--------------------------------------------------------------
!              Select halos for final merging trees
!               Only selected halos are dumped to
!               files for subsequent analysis
      SUBROUTINE SelectHalos
!--------------------------------------------------------------
        Nselected = 0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
       Do i=1,Lmax
          Selected(i) = .false.
          jTree(i)        = 0
       EndDo
!$OMP PARALLEL DO DEFAULT(SHARED)
       Do i =1,Nsnaps
           CALL ReadHalos(i,i)                    ! read halo catalogs
        EndDo
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i,nFound)
       Do i=1,Last
          nFound =0
          Do j=1,Ncats-Nsnaps
              If(iGood(j,i)>0)nFound =nFound +1
           EndDo
          If(nFound>10)Selected(i) =.true.
       EndDo
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (i,nFound)
!       Do i=1,Nhalos(1)
!          If(Selected(i))Then
!              nFound =0
!              Do j=1,Nsnaps
!                 If(iGood(1,i)>0)Then
!                    If(Snapshot(j)%halo(iGood(j,i))%circV>70.)Nfound =1 
!                    If(abs(Snapshot(1)%halo(iGood(1,i))%x-21.8)<2. .and.   &
!                    abs(Snapshot(1)%halo(iGood(1,i))%y-13.9)<2. .and.   &
!                    abs(Snapshot(1)%halo(iGood(1,i))%z-9.5)< 2. .and.   &
!                    abs(Snapshot(1)%halo(iGood(1,i))%circV >70. ))Selected(i) = .true.
!                !EndIf
!              !EndDo
!              !If(nFound==0)Selected(i) =.false.
!           !EndIf
!       EndDo
       j = 0
       Do i=1,Last          ! fill array of pointers to Tree
          If(Selected(i))Then
             j = j+1 ; jTree(i) = j
          EndIf
       EndDo
       Nselected =  j
       write(*,*) '  Number of selected halos =',nSelected
      Return
      End SUBROUTINE SelectHalos
!--------------------------------------------------------------
!              find major progenitors of selected halos
!              results are placed in Tree(iC,jj)
!                        iC - snapshot moment
!                        jj   - selected halo number
      SUBROUTINE Progenitors
!--------------------------------------------------------------
Character*120 :: Fname
                   ! start with already loaded snapshots
        Do iC=1,Nsnaps
           CALL  ListHalos(iC,Nhalos(iC),0)
           CALL HaloDistinct(Nhalos(iC),iC,0)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (j,jj,ih)
        Do j =1,Last
              If(Selected(j).and.iGood(iC,j)>0)Then
                 ih =  iGood(iC,j)           ! halo number in catalog iC
                 jj  = jTree(j)                  ! position in Tree
                 Tree(iC,jj)%nh  = ih
                 Tree(iC,jj)%Iso  = jNotFound(ih)
                 Tree(iC,jj)%x      = Snapshot(iC)%halo(ih)%x
                 Tree(iC,jj)%y      = Snapshot(iC)%halo(ih)%y
                 Tree(iC,jj)%z      = Snapshot(iC)%halo(ih)%z
                 Tree(iC,jj)%vx    = Snapshot(iC)%halo(ih)%vx
                 Tree(iC,jj)%vy    = Snapshot(iC)%halo(ih)%vy
                 Tree(iC,jj)%vz    = Snapshot(iC)%halo(ih)%vz
                 Tree(iC,jj)%aM   = Snapshot(iC)%halo(ih)%aM
                 Tree(iC,jj)%rM   = Snapshot(iC)%halo(ih)%Rm
                 Tree(iC,jj)%circV= Snapshot(iC)%halo(ih)%circV
              EndIf   ! Selected
        EndDo
            write(*,*) ' Done Tree for moment =',Aexpn(iC),iC
        EndDo
                   !  read new snapshots and add them to the Tree
        Do iC=Nsnaps+1,Ncats-Nsnaps
           CALL LoadSnap(iC)
           CALL ListHalos(1,Nhalos(iC),1)
           CALL HaloDistinct(Nhalos(iC),iC,1)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (j,jj,ih)
           Do j =1,Last
              If(Selected(j).and.iGood(iC,j)>0)Then
                 ih =  iGood(iC,j)           ! halo number in catalog
                 jj  = jTree(j)                  ! position in Tree
                 Tree(iC,jj)%nh  = ih
                 Tree(iC,jj)%Iso  =  jNotFound(ih)
                 Tree(iC,jj)%x      = Catalog(ih)%x
                 Tree(iC,jj)%y      = Catalog(ih)%y
                 Tree(iC,jj)%z      = Catalog(ih)%z
                 Tree(iC,jj)%vx    = Catalog(ih)%vx
                 Tree(iC,jj)%vy    = Catalog(ih)%vy
                 Tree(iC,jj)%vz    = Catalog(ih)%vz
                 Tree(iC,jj)%aM   = Catalog(ih)%aM
                 Tree(iC,jj)%circV= Catalog(ih)%circV
                 Tree(iC,jj)%Rm= Catalog(ih)%Rm
              EndIf   ! Selected
        EndDo
           ! write(*,*) ' Done Tree for moment =',Aexpn(iC),iC
        EndDo

        inFile =Nselected   /100
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (jFile,Fname,i,j,iC)
        Do jFile =1,100
           !jFile =1
           write(Fname,'(a,i3.3,a)') 'HaloTraceASCII.',jFile,'.dat'
           open(10+jFile,file=TRIM(Fname))
           write(Fname,'(a,i3.3,a)') 'HaloTrace.',jFile,'.dat'
           open(500+jFile,file=TRIM(Fname),form='unformatted')
           Do i =1,inFile
              j = i + (jFile-1)*InFile
              write(10+jFile,'(/a,i10)') ' ----  Halo =',j
              write(500+jFile) 
              write(500+jFile) j
              Do iC =1,Ncats-Nsnaps
                 If(Tree(iC,j)%nh>0) Then
                 write(10+jFile,'(f9.4,2i10,3f9.4,3f9.2,1pg12.4,0p2f8.2)')Aexpn(iC), Tree(iC,j)
                 write(500+jFile)Aexpn(iC), Tree(iC,j)
                 EndIf
              EndDo
           EndDo
           close(10+jFile) ; close(500+jFile)
      EndDo
        inFile =Nhalos(1)   /10
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (jFile,Fname,i,j,iC)
        Do jFile =1,10
           !jFile =1
           write(Fname,'(a,i3.3,a)') 'HaloTraceZ0.',jFile,'.dat'
           open(10+jFile,file=TRIM(Fname))
           Do i =1,inFile
              j = i + (jFile-1)*InFile
              write(10+jFile,'(/a,i10)') ' ----  Halo =',j
              Do iC =1,Ncats-Nsnaps
                 If(Tree(iC,j)%nh>0) &
                 write(10+jFile,'(f9.4,2i9,3f9.4,3f9.2,1pg12.4,0p2f8.2)')Aexpn(iC), Tree(iC,j)
              EndDo
           EndDo
           close(10+jFile) 
      EndDo
 
      Return
      End SUBROUTINE Progenitors
!---------------------------------------------------------
!              compare particles of halo 'ih' in Catalog iFirst with all halos in catalog iNext
!                  jh = identified halo in iG
!                  iCross = total number of identical particles
!                  Good = goodness of match: weighted
      SUBROUTINE IndPart(ih,iC,iFirst,iNext,jh,iCross,weight,iFlag)
!--------------------------------------------------------------
Integer  :: ListMax(100),jCross(100)
Real      :: wMax(100)

      Nc = Snapshot(iFirst)%halo(ih)%sel
      xh  = Snapshot(iFirst)%halo(ih)%x
      yh  = Snapshot(iFirst)%halo(ih)%y
      zh  = Snapshot(iFirst)%halo(ih)%z
      vxh  = Snapshot(iFirst)%halo(ih)%vx
      vyh  = Snapshot(iFirst)%halo(ih)%vy
      vzh  = Snapshot(iFirst)%halo(ih)%vz
      Vrms2  = Snapshot(iFirst)%halo(ih)%rmsV**2
      Vcirc   =Snapshot(iFirst)%halo(ih)%circV
      aM   =Snapshot(iFirst)%halo(ih)%aM
      Rmh2  =(1.e-3*Snapshot(iFirst)%halo(ih)%Rm)**2
      !ind =0
      !If(abs(xh-21.8)<1 .and.abs(yh-13.9)<1. .and. abs(zh-9.5)<1..and.aM>3.e10) Then
      !   ind = 1
       !write(*,'(a,i8,3x,3f9.4,f8.1,g12.3,f7.4)') ' halo =',ih,xh,yh,zh,Vcirc,Snapshot(iFirst)%halo(ih)%aM,Aexpn(iC)
      !EndIf
      xh = xh - vxh*dCorrect                ! extrapolate coordinates from moment iC to iG
      yh = yh - vyh*dCorrect
      zh = zh - vzh*dCorrect

      ListMax = 0
      wMax     =0.
      nCandidate = 0

              N0x = INT((xh-Rmax)/Cell+100)-101 ; N1x = (xh+Rmax)/Cell+1  ! limits for linker list
              N0y = INT((yh-Rmax)/Cell+100)-101 ; N1y = (yh+Rmax)/Cell+1 
              N0z = INT((zh-Rmax)/Cell+100)-101 ; N1z = (zh+Rmax)/Cell+1 
              !write(*,*) ' Cell=',Cell,' Rmax=',Rmax
              !write(*,*) ' Limits =',N0x,N1x
              !write(*,*) ' Limits =',N0y,N1y
              !write(*,*) ' Limits =',N0z,N1z
     Do k0 =N0z,N1z
        k1 =mod(k0+50*(Nbz+1),Nbz+1)  ! periodical boundaries
     Do j0 =N0y,N1y
        j1 =mod(j0+50*(Nby+1),Nby+1)  ! periodical boundaries
     Do i0 =N0x,N1x
        i1 =mod(i0+50*(Nbx+1),Nbx+1)  ! periodical boundaries
        j = Label(i1,j1,k1)
        Do while (j.ne.0)                         ! loop over all particles in cell
             Ng =Snapshot(iNext)%halo(j)%sel
             xj  = Snapshot(iNext)%halo(j)%x
             yj  = Snapshot(iNext)%halo(j)%y
             zj  = Snapshot(iNext)%halo(j)%z
             dx =  xj - xh ; dy =  yj - yh ; dz =  zj - zh 
             !  write(*,'(a,i8,i3,3x,8f9.4)') ' halo2=      ',j,Ng,xj,yj,zj

             dx = min(dx,abs(dx-Box),abs(dx+Box))                   ! use periodical boundaries
             dy = min(abs(dy),abs(dy+Box),abs(dy-Box)) 
             dz = min(abs(dz),abs(dz+Box),abs(dz-Box)) 
              dd = dx**2 +dy**2 +dz**2
              nCross = 0
             If(dd < d2max)Then                                                      ! check only halos with small differences
                If(Ng.ne.0.and.Nc.ne.0)Then
                  Do i=1,Nc                                                                   ! loop through all particles of first halo
                    ii =Snapshot(iFirst)%halo(ih)%bound(i) 
                    Do k=1,Ng                                                               ! loop through all particles of second halo
                       kk =Snapshot(iNext)%halo(j)%bound(k)
                       if(kk==ii)nCross =nCross +1
                    EndDo
                  EndDo
                 EndIf   ! Ng /= 0
               dv = ((Snapshot(iNext)%halo(j)%vx -  vxh)**2 + &
                                (Snapshot(iNext)%halo(j)%vy - vyh)**2 + &
                                (Snapshot(iNext)%halo(j)%vz - vzh)**2)
               Wdist =min(0.05*Rmh2/dd,1.5)
               If(nCross > 25)Then
                   w_part = 1.5
               Else
                  w_part = float(nCross)/(Ng+Nc)*2.  
               EndIf
               Vrat = Snapshot(iNext)%halo(j)%circV/Vcirc
               Vrat = max(Vrat,1./Vrat)
               wv    = min(-log10(Vrat-1.),1.)
               ratM = Snapshot(iNext)%halo(j)%aM/aM
               ratM = max(ratM,1./ratM)
               wm    = min(-1.5*log10(ratM-1.),1.)

                 w =  w_part                                                                + &      ! number of common particles
                         Wdist                                                                 + &      ! distance
                         min(0.25*Vrms2/max(dv,100.),1.5)                          + &      ! velocity
                         wv + & ! Vcirc
                         wm  !  Mass
      !If(ind == 1.and. w >1.) &
      ! write(*,'(a,i8,3i4,3x,f9.4,2f9.2,3f9.4,3f8.2,g12.4)') '     j=',j,nCross,Nc,Ng,w,sqrt(dv),1000.*sqrt(dd), &
       !     Snapshot(iNext)%halo(j)%x, Snapshot(iNext)%halo(j)%y, Snapshot(iNext)%halo(j)%z, & 
       !     Snapshot(iNext)%halo(j)%vx, Snapshot(iNext)%halo(j)%vy, Snapshot(iNext)%halo(j)%vz, & 
       !     Snapshot(iNext)%halo(j)%aM
              Else
                w =0.
             EndIf      ! dd < d2max
             If(w > wMin)Then           ! this is a candidate for halo matching
                nCandidate = nCandidate +1
                if(nCandidate > 100) Stop ' Too many candidates for matching halos. Increase length of ListMax'
                ListMax(nCandidate) = j 
                jCross(nCandidate)  = nCross
                wMax(nCandidate)    = w
             EndIf
             j = Lst(j)
      EndDo     ! j - loop inside cell
    EndDo       ! kc
    EndDo       ! jc 
    EndDo       ! ic 

    If(nCandidate > 0)Then      !  there are candidates. Find the best match
       jmx = Transfer(MAXLOC(wMax),0)   ! position of maximum
       jh    =  ListMax(jmx)
       weight  = wMax(jmx)
       iCross  = jCross(jmx)
    Else
       jh    =  0
       weight  = 0
       iCross  = 0
    EndIf

         If(iFlag==0.and.ih<200.and.jh/=0)Then
               ! write(1,'(i6,3x,i3,T76,f9.3)') ih,Nc,Snapshot(iFirst)%halo(ih)%circV
          xj  = Snapshot(iNext)%halo(jh)%x
          yj  = Snapshot(iNext)%halo(jh)%y
          zj  = Snapshot(iNext)%halo(jh)%z
          dx =  xj - xh ; dy =  yj - yh ; dz =  zj - zh ; 
          dx = 1.e3*min(abs(dx),abs(dx+Box),abs(dx-Box))  ! periodical boundaries
          dy = 1.e3*min(abs(dy),abs(dy+Box),abs(dy-Box)) 
          dz = 1.e3*min(abs(dz),abs(dz+Box),abs(dz-Box)) 
          dd = sqrt(dx**2 +dy**2 +dz**2)
          dv = sqrt( (Snapshot(iNext)%halo(jh)%vx - &
                              Snapshot(iFirst)%halo(ih)%vx)**2 + &
                             (Snapshot(iNext)%halo(jh)%vy - &
                              Snapshot(iFirst)%halo(ih)%vy)**2 + &
                              (Snapshot(iNext)%halo(jh)%vz - &
                              Snapshot(iFirst)%halo(ih)%vz)**2 )
                write(1,10) ih,jh, Snapshot(iNext)%halo(jh)%sel,iCross,&
                     dx,dy,dz,dd, &
                     Snapshot(iNext)%halo(jh)%circV, &
                     dv, weight
     EndIf
10    format(9x,3i6,3x,i5,T40,4f8.2,T76,2f8.2,3x,1p3g11.3)

      End SUBROUTINE IndPart

!-------------------------------------------------------------------------
!                    Make linker list of Halos
      SUBROUTINE ListHalos(iNext,Ncenter,iFlag)
!------------------------------------------------------------------------
!write(*,*) ' Moment=',iNames(iNext),' a=',Aexpn(iNext)
!write(*,*) ' Ncenter=',Ncenter
!write(*,*) ' Limits=',Nmx,Nbx
!write(*,*) ' Limits=',Nmy,Nby
!write(*,*) ' Limits=',Nmz,Nbz
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
             Do i=1,Ncenter
                Lst(i)= -1
             EndDo

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,k)
       Do k=Nmz,Nbz
       Do j=Nmy,Nby
       Do i=Nmx,Nbx
          Label(i,j,k)=0
       EndDo
       EndDo
       EndDo
       !write(*,*) ' Ncenter=',Ncenter,' Cell=',Cell
      If(iFlag ==0)Then                       ! use Snapshot
         Do jp=1,Ncenter
            i=INT(Snapshot(iNext)%halo(jp)%x/Cell)
            j=INT(Snapshot(iNext)%halo(jp)%y/Cell)
            k=INT(Snapshot(iNext)%halo(jp)%z/Cell)
            i=MIN(MAX(Nmx,i),Nbx)
            j=MIN(MAX(Nmy,j),Nby)
            k=MIN(MAX(Nmz,k),Nbz)
            Lst(jp)      =Label(i,j,k)
            Label(i,j,k) =jp
         EndDo
      Else                                          ! use Catalog
         Do jp=1,Ncenter
            i=INT(Catalog(jp)%x/Cell)
            j=INT(Catalog(jp)%y/Cell)
            k=INT(Catalog(jp)%z/Cell)
            i=MIN(MAX(Nmx,i),Nbx)
            j=MIN(MAX(Nmy,j),Nby)
            k=MIN(MAX(Nmz,k),Nbz)
            Lst(jp)      =Label(i,j,k)
            Label(i,j,k) =jp
         EndDo
      EndIf

       write(*,*) ' Done Linker List'
      Return
      End SUBROUTINE ListHalos
!-------------------------------------------------------------------------
!                    Initialize arrays and structures
      SUBROUTINE Initialize
!------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED)
   Do i=1,Lmax                                   ! initialize tables 
      iGood(:,i)   =0
      wGood(:,i) =0.
   EndDo
!$OMP PARALLEL DO DEFAULT(SHARED)
   Do i=1,Nmax                                   ! initialize tables 
      Catalog(i)%nh =0
      Catalog(i)%Iso = -1
      Catalog(i)%x   = 0.    ; Catalog(i)%y   = 0.   ; Catalog(i)%z   = 0. 
      Catalog(i)%vx = 0.    ; Catalog(i)%vy = 0. ; Catalog(i)%vz = 0. 
      Catalog(i)%aM = 0.
      Catalog(i)%circV = 0.
      Do j=1,Nsnaps
          iGlobal(j,i) = i
       EndDo
  EndDo
Do iSnap =1,Nsnaps
!$OMP PARALLEL DO DEFAULT(SHARED)
Do i=1,Nmax
                Snapshot(iSnap)%halo(i)%num = 0
                 Snapshot(iSnap)%halo(i)%x = 0. 
                 Snapshot(iSnap)%halo(i)%y = 0.
                 Snapshot(iSnap)%halo(i)%z = 0.
                 Snapshot(iSnap)%halo(i)%vx = 0.
                 Snapshot(iSnap)%halo(i)%vy = 0.
                 Snapshot(iSnap)%halo(i)%vz = 0.
                 Snapshot(iSnap)%halo(i)%aM = 0.
                 Snapshot(iSnap)%halo(i)%Rm = 0.
                 Snapshot(iSnap)%halo(i)%rmsV = 0.
                 Snapshot(iSnap)%halo(i)%circV = 0.
                 Snapshot(iSnap)%halo(i)%sel = 0
                 Snapshot(iSnap)%halo(i)%bound(:) = 0
EndDo
EndDo

Do iSnap =1,Ncats-Nsnaps
!$OMP PARALLEL DO DEFAULT(SHARED)
Do i=1,Nmax
                Tree(iSnap,i)%nh = 0
EndDo
EndDo

      End SUBROUTINE INITIALIZE
!-------------------------------------------------------------------------
!                    read first snapshots
      SUBROUTINE START
!------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED)
       Do i =1,Nsnaps
           CALL ReadHalos(i,i)                    ! read halo catalogs
        EndDo
        write(*,*) ' finished reading'
!$OMP PARALLEL DO DEFAULT(SHARED)
   Do  ih=1,Nhalos(1)   ! first column initialized to all halos  
      iGood(1,ih) = ih
      wGood(1,ih) = 4.
   EndDo 
   
      End SUBROUTINE START

End MODULE SetHalos

!-------------------------------------------------
!            Find halos at different moments, which
!            have common particles
PROGRAM HaloCrossIdentity

use     SetHalos
Character*120                                     :: file1,listname,txt*8     
!     ------------------------------------------        Read data
      Open(1,file='HaloCrossList.dat')

      Call CPU_time(t0); call system_clock(iT0,iTc,ib)
   Do iC=1,Ncats
      write(*,'(a,$)')' Enter snapshot='
      read(*,'(a)')iNames(iC)
   endDo
   Cell = Box/(Nbx-Nmx+1)          ! cell size for linker list

   Ncatalogs = Ncats
   nc =0
   nm =0

   CALL Initialize                             ! initialize arrays and structures
   write(*,*)
   Call CPU_time(t1); call system_clock(iT2,iTc,ib) ; write(*,*) '          time for Init =',t1-t0,float(iT2-iT0)/iTc

   CALL START                               ! Read first snapshots
   Last = Nhalos(1)   ! last halo in the list
   Call CPU_time(t2); call system_clock(iT3,iTc,ib) ; write(*,*) '          time for START =',t2-t1,float(iT3-iT2)/iTc
                              ! ---------------------------------------------   Main loop 
Do iC =1, Ncats-Nsnaps
   Call CPU_time(t2);call system_clock(iT2,iTc,ib) 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
   Do i=1,Nhalos(iC)
      jNotFound(i) = 0
   EndDo
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
   Do i=1,Nmax
      iGlobal(:,i)     = 0
   EndDo

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i) 
    Do i =1,Last                  ! find pointers of halos to merging tree entries
       Do j=1,Nsnaps
          If(iGood(iC+j-1,i).ne.0)iGlobal(j,iGood(iC+j-1,i)) =i 
       EndDo
    Enddo
   iFirst = mod(iC-1,Nsnaps)+1                                                ! pointer to first Snapshot 
    Do i=1,Nhalos(iC)
       If(iGlobal(1,i) == 0)Then 
         If( Snapshot(iFirst)%halo(i)%circV>40.)Then          ! If halo 'i' does not have tree entry, then
           Last = Last +1                       ! add  halo to the list of progenitors: create a new line
           if(Last>Lmax)Stop ' Not enough space for halo progenitors:Lmax'
              iGood(iC,Last) = i
              wGood(iC,Last) = 4.
              iGlobal(1,i) =Last 
              !write(1,'(a,6i6)') '       halo=',i,iGlobal(i),iGood(iC,Last)
              Else
               jNotFound(i) = -1  
         EndIf
       EndIf
    EndDo       
   !write(1,'(3a,i10,a,36i8)') ' Catalog=',TRIM(iNames(iC)),' Halos=',Last,' Pointers :',(iGlobal(1,i),i=1,min(Nhalos(iC),20))
   !write(1,'(29x,a,36i8)') ' Pointers2:',(iGlobal(2,i),i=1,min(Nhalos(iC+1),20))

   Do jCompare = 1,Nsnaps-1                                                    ! loop over next  snapshots
      iG = iC + jCompare                                                                ! absolute index of next snapshot
      iNext  = mod(iFirst +jCompare-1,Nsnaps)+1                    ! pointer to comparison Snapshot
      Call ListHalos(iNext,Nhalos(iG),0)                                     ! make Linker List
                       !write(*,*)  '  Comparison of ',iC,iFirst,' with ',iG,iNext
      da = Aexpn(iC) - Aexpn(iG)
      dCorrect = 0.01*da/sqrt(Omega0*Aexpn(iC)+(1.-Omega0)*Aexpn(iC)**4)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ih,j,nn,weight,jj,xh,yh,zh,aM) 
   Do ih =1,Nhalos(iC)                                        ! find halos for iNext moment 
      If(jNotFound(ih) == 0)Then                        ! consider only halos not identified yet                     
        CALL IndPart(ih,iC,iFirst,iNext,j,nn,weight,1)
        !if(mod(ih,500000)==0)write(*,*) ' Halo=',ih,j,weight
        If(j/= 0)Then
           If(iGood(iG,iGlobal(1,ih))==0.and.iGlobal(jCompare+1,j)==0)Then ! if halo not taken, take it
               iGood(iG,iGlobal(1,ih)) = j
               wGood(iG,iGlobal(1,ih)) = weight
               iGlobal(jCompare+1,j) = iGlobal(1,ih)
               jNotFound(ih) =1
      xh  = Snapshot(iFirst)%halo(ih)%x
      yh  = Snapshot(iFirst)%halo(ih)%y
      zh  = Snapshot(iFirst)%halo(ih)%z               
      aM  = Snapshot(iFirst)%halo(ih)%aM             
      !If(abs(xh-21.8)<1. .and.abs(yh-13.9)<1. .and. abs(zh-9.5)<1. .and.aM>3.e10)&
      !         write(*,'(a,2i8,a,i8,f6.2,a,2i9)')'    0   ih,j=',ih,j,' Good=',iGood(iG,iGlobal(1,ih)),wGood(iG,iGlobal(1,ih)), &
      !                                                                ' current=',iGlobal(1,ih)
           Else
              !if(j<12)write(1,'(a,2i8,a,i8,f6.2,a,2i9)')'       ih,j=',ih,j,' Good=',iGood(iG,iGlobal(jCompare+1,j)), &
              !                                                                          wGood(iG,iGlobal(jCompare+1,j)),     &
              !                                                                         ' current=',iGlobal(1,ih) ,iGlobal(jCompare+1,j)           
             If(wGood(iG,iGlobal(jCompare+1,j))<weight)Then     ! If halo was taken, but had smaller weight, release the 

              jj = iGlobal(jCompare+1,j)
              iGood(iG,jj) = 0 ; wGood(iG,jj) = 0.
               iGood(iG,iGlobal(1,ih)) = j
               wGood(iG,iGlobal(1,ih)) = weight 
               iGlobal(jCompare+1,j) = iGlobal(1,ih)
               jNotFound(ih) =1               

      !If(abs(xh-21.8)<1. .and.abs(yh-13.9)<1. .and. abs(zh-9.5)<1..and.aM>3.e10)&
      !         write(*,'(a,2i8,a,i8,f6.2,a,2i9)')'    1   ih,j=',ih,j,' Good=',iGood(iG,iGlobal(1,ih)),wGood(iG,iGlobal(1,ih)), &
      !                                                                ' current=',iGlobal(1,ih)
           Else
                iGood(iG,iGlobal(1,ih)) = 0
                wGood(iG,iGlobal(1,ih)) = 0.
            EndIf    ! wGood
           EndIf     ! iGood
      
        Else
           iGood(iG,iGlobal(1,ih)) = 0
           wGood(iG,iGlobal(1,ih)) = 0.
           jNotFound(ih) = 0 
        EndIf  !  j/=0
     EndIf     ! jNotFound ==0
   EndDo    ! ih

   EndDo    ! jCompare

   Call CPU_time(t3) ; call system_clock(iT3,iT0,ib)
   write(*,'(3(a,i10),a,2f8.2)') ' Current Halos=',Nhalos(iC),' Total Halos=',Last,' iC=',iC,  &
                    '   time for step =',t3-t2,float(iT3-iT2)/iT0
   CALL ReadHalos(iC+Nsnaps,iFirst)
EndDo

!                                          Collect information on each halo in the tree
!                                          Dump the tree
write(*,*) '  Total Number of halos in the tree=',Last
Nempty = 0
Do i =1,Last
   Non_zero =0
   Do j =1,Ncats-2
      If(iGood(j,i).ne.0)Non_zero=Non_zero +1
   EndDo
   If(Non_zero>3)Nempty =Nempty+1
EndDo
write(*,*) '  Number of halos found more than 3 times = ',Nempty
write(*,*) ' -------------------------  Map of progenitors ------------------ '
write(*,'(6x,40f9.3)') (Aexpn(iC),iC=1,Min(Ncats-2,20))
write(*,'(6x,40a9)') (TRIM(iNames(iC)),iC=1,Min(Ncats-2,20))
CALL ReadHalos(1,1)
Do j=1,600
If( Snapshot(1)%halo(j)%circV >80.)   write (*,'(i8,2x,40i9)') j,(iGood(iC,j),iC=1,Min(Ncats-2,20))
EndDo
Do j=1,Nhalos(1)
If( Snapshot(1)%halo(j)%circV >800.)   write (*,'(i8,2x,40i9)') j,(iGood(iC,j),iC=1,Min(Ncats-2,20))
EndDo

CALL SelectHalos
CALL Progenitors

!Do ih=1,300
!If(abs(Snapshot(1)%halo(ih)%x-37.3)<0.7 .and. abs(Snapshot(1)%halo(ih)%y-5.1)<0.7) &
!If( Snapshot(1)%halo(ih)%circV>70. .and. Snapshot(1)%halo(ih)%vz<0.) &
!              Call MajorProgenitor(ih)
!EndDo

stop
end Program HaloCrossIdentity
