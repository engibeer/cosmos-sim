c     ====================================================================
c                                    
c                ART Version 4: Find Large Particles
c               Anatoly Klypin, Stefan Gottloeber  (2002)                  
c
c     For current node defined by iN0(i,j,k):
c          - create meshes used for large particles
c          - loop over all nodes iN1 { 
c                            find Boundary(iN0,iN1);
c                            find Small particles, add them to send buffer    
c                                              } end loop over nodes      
c          - complete send buffer 
c          - send/receive Small particles
c          - append Small particles to x(),vx(),pt(),wpart(),i_par()
c          - repeat all the steps of the loop for Large particles
c     ====================================================================
c-------------------------------------------------------------------- 
      SUBROUTINE BoundNode(iN0,iN1,Nbound)
c-------------------------------------------------------------------- 
c
c     purpose :   find boundaries of two nodes in 3D
c     input   :   iN0=(i0,j0,k0)  - position of the node, which creates particles
c                 iN1 =(i1,j1,k1) - node, which receives particles
c     output  : Nbound - number of boundaries: 0 -no boundaries
c                                              -1 = the same node
c                                          there are  maximum 8 boundaries
c               isBound(8,6) - boundaries in units of zero-level cells
c               iBound(8,6) - boundaries in units of lowest resolution mass mesh
c                           - first index is the boundary number
c                           - second index is the coordinate of the boundary 
c                  i(s)Bound(i,j): j=1: x_left,  j=2: x_right,  
c                                    j=3: y_left,  j=4: y_right,   j=5: z_left,  j=6: z_right
c                mBound()  - length of mesh for largest particles 
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      DIMENSION iN0(3),iN1(3),jBound(3),kBound(3,2,2)

      Nbound = 0
      If((iN0(1).eq.iN1(1)).and.      ! the same node -> return
     &  (iN0(2).eq.iN1(2)).and.
     &  (iN0(3).eq.iN1(3)))Then
         Nbound =-1
         Return
      EndIf 
      iScale =2**(N_lev_lg-1)/N_sub_lg

      Do idir =1,3
         If(idir.eq.1)Then                 ! find boundaries of iN0 (creating) node
                                                     !               and  iN1 (receiving) node
           i0Left       =  l_divX(iN0(1))
           i0Right     =  l_divX(iN0(1)+1)
           i1Left       =  l_divX(iN1(1))
           i1Right     =  l_divX(iN1(1)+1)

         Else if (idir.eq.2)Then
           i0Left       =  l_divY(iN0(1),iN0(2))
           i0Right     =  l_divY(iN0(1),iN0(2)+1)
           i1Left       =  l_divY(iN1(1),iN1(2))
           i1Right     =  l_divY(iN1(1),iN1(2)+1)
         Else
           i0Left       =  l_divZ(iN0(1),iN0(2),iN0(3))
           i0Right     =  l_divZ(iN0(1),iN0(2),iN0(3)+1)
           i1Left       =  l_divZ(iN1(1),iN1(2),iN1(3))
           i1Right     =  l_divZ(iN1(1),iN1(2),iN1(3)+1)
         EndIf
         jBound(idir) = 0
            If(iDebug_lg)write (13,50) i0Left,i0Right,idir,mBound(idir)
 50         format(' i0Left/Right=',2i4,' idir=',i3, ' MeshSize=',i3)
            Do iperiod = -1,1   ! periodical boundaries
               iLeft   = max(i0Left,i1Left - L_buff_lg+iperiod*ng) 
               iRight  = min(i0Right,i1Right + L_buff_lg+iperiod*ng)
               If(iLeft .lt. iRight)Then   ! there is a boundary
c               write (*,10) i1Left,i1Right,iLeft, iRight,iperiod
c 10            format('    i1Left/Right=',2i4,' Bound=',4i4)
                  jBound(idir)     = jBound(idir) +1
                  kBound(idir,jBound(idir),1) = iLeft
                  kBound(idir,jBound(idir),2) = iRight
               EndIf
            EndDo    ! end 1D test along idir direction
         If(jBound(idir).eq.0)then   ! no common boundaries 
            Nbound = 0
            iBound(1,1) = 0
            Return
         EndIf
      EndDo           ! end all directions tests

      If(iDebug_lg)write (13,*) '       iScale=', iScale

       Nbound = 0     
           i0Left       =  l_divX(iN0(1))   ! low left corner of iN0 node
           j0Left       =  l_divY(iN0(1),iN0(2))
           k0Left       =  l_divZ(iN0(1),iN0(2),iN0(3))

      Do i =1, jBound(1)
      Do j =1, jBound(2)
      Do k =1, jBound(3)
       Nbound =        Nbound +1
       isBound(Nbound,1) =kBound(1,i,1)  !  x boundaries
       isBound(Nbound,2) =kBound(1,i,2)
       isBound(Nbound,3) =kBound(2,j,1)  !  y boundaries
       isBound(Nbound,4) =kBound(2,j,2)
       isBound(Nbound,5) =kBound(3,k,1)  !  z boundaries
       isBound(Nbound,6) =kBound(3,k,2)
         iBound(Nbound,1) =(isBound(Nbound,1)-i0Left)/iScale+1
         iBound(Nbound,2) =(isBound(Nbound,2)-i0Left-1)/iScale+1
         iBound(Nbound,3) =(isBound(Nbound,3)-j0Left)/iScale+1
         iBound(Nbound,4) =(isBound(Nbound,4)-j0Left-1)/iScale+1
         iBound(Nbound,5) =(isBound(Nbound,5)-k0Left)/iScale+1
         iBound(Nbound,6) =(isBound(Nbound,6)-k0Left-1)/iScale+1
         if(iDebug_lg)
     & write (13,30) Nbound,mBound,(iBound(Nbound,jj),jj=1,6)
 30    format(' Nbound=',i4,' MeshSize=',3i4,3x,
     &             ' Boundary: x=',2i4,' y=',2i4,' z=',2i4)
      EndDo         
      EndDo         
      EndDo         

      Return
      End
c-------------------------------------------------------------------- 
      SUBROUTINE MakeMapNode(Nbound)
c-------------------------------------------------------------------- 
c          
c     purpose :   construct map of mass levels for Large Particles
c     input   :   Nbound - number of boundaries: 0 (no boundaries)
c                 iBound(8,6) - boundaries in units of lowest resolution mass mesh 
c     output  : MeshLG(i,j,k) - map of levels 
c               = 1 for first level (highest mass resolution)
c               = 2 second level: 8 times larger cells
c               = 3 third level
c               = 4 forth level (lowest mass resolution) = N_lev_lg 
c               MeshLG has the same size as the lowest mass resolution meshes
c               
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
 
      If(N_lev_lg.eq.4)Then
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k)
         Do k= 1,N_m4_lg    ! initialize the map
         Do j= 1,N_m4_lg
         Do i= 1,N_m4_lg
            MeshLG(i,j,k) = N_lev_lg
         EndDo
         EndDo
         EndDo
      Else
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k)
         Do k= 1,N_m3_lg    ! initialize the map
         Do j= 1,N_m3_lg
         Do i= 1,N_m3_lg
            MeshLG(i,j,k) = N_lev_lg
         EndDo
         EndDo
         EndDo
      EndIf

      If(Nbound.le.0) Return
c      DO ib = 1,Nbound
c      write(13,*) 'ib =',ib 
c      write(13,*)  iBound(ib,5),iBound(ib,6),
c     +    iBound(ib,3),iBound(ib,4),
c     +     iBound(ib,1),iBound(ib,2)
c      ENDDO

         Do ib = 1,Nbound  ! mark highest resolution cells
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k)
            Do k =iBound(ib,5),iBound(ib,6)
            Do j =iBound(ib,3),iBound(ib,4)
            Do i =iBound(ib,1),iBound(ib,2)
               MeshLG(i,j,k) = 1
            EndDo
            EndDo
            EndDo
         EndDo   ! end loop over all boundaries
         mboundX = mBound(1)   ! length of the mesh
         mboundY = mBound(2)
         mboundZ = mBound(3)


         Do k= 1,mboundZ   ! add extra shell of high-res cells 
         Do j= 1,mboundY
         Do i= 1,mboundX
            If(MeshLG(i,j,k).eq.N_lev_lg)Then ! unmarked cell
             kp =min(k+1, mboundZ) 
             jp =min(j+1, mboundY) 
             ip =min(i+1, mboundX)
             km = max(k-1,1)
             jm = max(j-1,1)
             im = max(i-1,1)
             Do kk=km,kp  ! mark the cell if a neighbor is 1 level down
             Do jj=jm,jp
             Do ii=im,ip
                If(MeshLG(ii,jj,kk).eq.1)MeshLG(i,j,k)=-2
             EndDo
             EndDo
             EndDo
            EndIf     ! end Mesh = 0
         EndDo
         EndDo
         EndDo

         Do k= 1,mboundZ   ! add extra shell of high-res cells 
         Do j= 1,mboundY
         Do i= 1,mboundX
            If(MeshLG(i,j,k).eq.N_lev_lg)Then ! unmarked cell
             kp =min(k+1, mboundZ) 
             jp =min(j+1, mboundY) 
             ip =min(i+1, mboundX)
             km = max(k-1,1)
             jm = max(j-1,1)
             im = max(i-1,1)
             Do kk=km,kp  ! mark the cell if a neighbor is 1 level down
             Do jj=jm,jp
             Do ii=im,ip
                If(MeshLG(ii,jj,kk).eq.-2)MeshLG(i,j,k)=2
             EndDo
             EndDo
             EndDo
            EndIf     ! end Mesh = 0
         EndDo
         EndDo
         EndDo

!         Do k= 1,mboundZ   ! add extra shell of high-res cells 
!         Do j= 1,mboundY
!         Do i= 1,mboundX
!            If(MeshLG(i,j,k).eq.N_lev_lg)Then ! unmarked cell
!             kp =min(k+1, mboundZ) 
!             jp =min(j+1, mboundY) 
!             ip =min(i+1, mboundX)
!             km = max(k-1,1)
!             jm = max(j-1,1)
!             im = max(i-1,1)
!             Do kk=km,kp  ! mark the cell if a neighbor is 1 level down
!             Do jj=jm,jp
!             Do ii=im,ip
!                If(MeshLG(ii,jj,kk).eq.2)MeshLG(i,j,k)=-2
!             EndDo
!             EndDo
!             EndDo
!            EndIf     ! end Mesh = 0
!         EndDo
!         EndDo
!         EndDo


         Do k= 1,mboundZ   ! add extra shell of high-res cells 
         Do j= 1,mboundY
         Do i= 1,mboundX
                If(MeshLG(i,j,k).eq.-2)MeshLG(i,j,k)=2
         EndDo
         EndDo
         EndDo

         
         Do m =2,N_lev_lg   ! mark neighbors of marked cells
         Do k= 1,mboundZ   
         Do j= 1,mboundY
         Do i= 1,mboundX
            If(MeshLG(i,j,k).eq.N_lev_lg)Then ! unmarked cell
             kp =min(k+1, mboundZ) 
             jp =min(j+1, mboundY) 
             ip =min(i+1, mboundX)
             km = max(k-1,1)
             jm = max(j-1,1)
             im = max(i-1,1)
             Do kk=km,kp  ! mark the cell if a neighbor is 1 level down
             Do jj=jm,jp
             Do ii=im,ip
                If(MeshLG(ii,jj,kk).eq.m-1)MeshLG(i,j,k)=m
             EndDo
             EndDo
             EndDo
            EndIf     ! end Mesh = 0
         EndDo
         EndDo
         EndDo
         EndDo          

!         m1 =0
!         m2 =0
!         m3 =0
!         Do k= 1,mboundZ   ! add extra shell of high-res cells 
!         Do j= 1,mboundY
!         Do i= 1,mboundX
!                If(MeshLG(i,j,k).eq.1)m1 =m1 +1
!                If(MeshLG(i,j,k).eq.2)m2 =m2 +1
!                If(MeshLG(i,j,k).eq.3)m3 =m3 +1
!         EndDo
!         EndDo
!         EndDo
!         write(13,*) ' Mesh =',mBoundX,mboundY,mboundZ
!         write(13,*) ' n cells=',m1,m2,m3
c         If(iDebug_lg)Then
c         Write(13,*) '---- Map: length in Z=',mboundZ
c         Do k=1,mboundZ
c            write(13,*) '      k =',k
c            Do j=1,mboundY
c               write(13,10) (MeshLG(i,j,k),i=1,mboundX )
c            EndDo 
c         EndDo
c 10      format(16i3)
c         EndIf   ! end debug
      Return
      End
c------------------------------------------------------------- 
      SUBROUTINE TestMeshes
c------------------------------------------------------------- 
c          
c     purpose :   test grids with coordinates, velocities,
c                        and masses of large particles 
c     input   :   
c                 
c     output  : 
      include 'a_control.h'
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
       REAL*8  sum1(8),sum2(8),sum3(8),sum4(8)

      Do i=1,8          ! initialize test counters
         sum1(i) =0.
         sum2(i) =0.
         sum3(i) =0.
         sum4(i) =0.
      EndDo          

c             test whether the mass, momentum, coords  are preserved
            Do k=1,N_m1_lg
            Do j=1,N_m1_lg
            Do i=1,N_m1_lg
               sum1(8) = sum1(8) +ReshLG1(8,i,j,k)
               Do m =1,6
                  sum1(m) =sum1(m) +ReshLG1(m,i,j,k)*ReshLG1(8,i,j,k)
               EndDo 
            EndDo
            EndDo
            EndDo
            Do k=1,N_m2_lg
            Do j=1,N_m2_lg
            Do i=1,N_m2_lg
               sum2(8) = sum2(8) +ReshLG2(8,i,j,k)
               Do m =1,6
                  sum2(m) =sum2(m) +ReshLG2(m,i,j,k)*ReshLG2(8,i,j,k)
               EndDo 
            EndDo
            EndDo
            EndDo
            Do k=1,N_m3_lg
            Do j=1,N_m3_lg
            Do i=1,N_m3_lg
               sum3(8) = sum3(8) +ReshLG3(8,i,j,k) 
               Do m =1,6
                  sum3(m) =sum3(m) +ReshLG3(m,i,j,k)*ReshLG3(8,i,j,k) 
               EndDo 
            EndDo
            EndDo
            EndDo
            Do k=1,N_m4_lg
            Do j=1,N_m4_lg
            Do i=1,N_m4_lg
               sum4(8) = sum4(8) +ReshLG4(8,i,j,k)
               Do m =1,6
                  sum4(m) =sum4(m) +ReshLG4(m,i,j,k)*ReshLG4(8,i,j,k)
               EndDo 
            EndDo
            EndDo
            EndDo
            write (13,20) 1,(sum1(i)/sum1(8),i=1,6),sum1(8)
            write (13,20) 2,(sum2(i)/sum2(8),i=1,6),sum2(8)
            write (13,20) 3,(sum3(i)/sum3(8),i=1,6),sum3(8)
            write (13,20) 4,(sum4(i)/sum4(8),i=1,6),sum4(8)
 20         format(' Level ',i2,1P,' <x>=',3g12.5,' <v>=',3g12.5,
     &                '  sum mass=',g12.5) 
            Do k=1,8
               write (13,*)' --------------- k=',k,'  X: meshes 1 and 2'
               Do j=1,8
                  write (13,30) (ReshLG1(1,i,j,k),i=1,8),
     &                             (ReshLG2(1,i,j,k),i=1,8)
               EndDo 
            EndDo 
            Do k=1,8
               write (13,*)' --------------- k=',k,'     X: 3 4'
               Do j=1,8
                  write (13,30) (ReshLG3(1,i,j,k),i=1,8),
     &                             (ReshLG4(1,i,j,k),i=1,8)
               EndDo 
            EndDo 
            Do k=1,8
               write (13,*) ' ---------------- k=',k,'     M: 1 2'
               Do j=1,8
                  write (13,30) (ReshLG1(8,i,j,k),i=1,8),
     &                             (ReshLG2(8,i,j,k),i=1,8)
               EndDo 
            EndDo 
            Do k=1,8
               write (13,*) ' ---------------- k=',k,'     M: 3 4'
               Do j=1,8
                  write (13,30) (ReshLG3(8,i,j,k),i=1,8),
     &                             (ReshLG4(8,i,j,k),i=1,8)
               EndDo 
            EndDo 

 30               format(8f8.2,5x,8f8.2) 
      Return
      End
c-------------------------------------------------------------------- 
      SUBROUTINE MakeMeshes(iN0,np_node)
c-------------------------------------------------------------------- 
c          
c     purpose :   construct grids with coordinates, velocities,
c                        and masses of large particles 
c     input   :   iN0 - node coordinates
c                 np_node - number of primary particles in this node    
c     output  : ReshLG#(8,i,j,k) = grids at level #  
c               total number of grids is N_lev_lg 
c                ReshLG#(1,...) = x, (2,... =y, z, vx, vy,vz, pt, mass
c                mBound(3) - size of mesh for particles 
c               # = 4 forth level (lowest mass resolution) = N_lev_lg 
c     1-level grid first collects Sum_i(x_i*Wpart_i).          
c     each cell on 2,3,4-levels give sums of 8 contributions
c                from lower level grid
c     At the end every grid is re-normalized to give <x>,<V>,mass
c                by deviding each cell by its mass
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      PARAMETER  ( xmax_b = ng +1.-5.e-5)
      DIMENSION iN0(3)

      iScale =2**(N_lev_lg-1)/N_sub_lg
         iLeft       =  l_divX(iN0(1))               ! x boundary
         jLeft       =  l_divY(iN0(1),iN0(2))               ! y
         kLeft       =  l_divZ(iN0(1),iN0(2),iN0(3))   ! z
         iRight      =  l_divX(iN0(1)+1)               ! x boundary
         jRight      =  l_divY(iN0(1),iN0(2)+1)               ! y
         kRight      =  l_divZ(iN0(1),iN0(2),iN0(3)+1)   ! z
         xLeft       =  iLeft
         yLeft       =  jLeft
         zLeft       =  kLeft

         mBound(1) = (iRight-iLeft-1)/iScale+1
         mBound(2) = (jRight-jLeft-1)/iScale+1
         mBound(3) = (kRight-kLeft-1)/iScale+1

c         If(mBound(idir)*iScale.ne.(i0Right-i0Left))Then
cc            write (13,*) ' Node =',irank+1
c            write (13,*) ' Error: incompatible mesh size: mesh=',
c     &              i0Right-i0Left,'   N_sub_lg=',N_sub_lg, iScale    
c            name = 'ERROR_RUN'
c            accesstype = 'append'
c            call Open_ASCII_FILE(19,name,accesstype)
cc            rewind 19
c            j123=123
c            write(19,*) j123, ' Error: incompatible mesh size: mesh=',
c     &              i0Right-i0Left,'   N_sub_lg=',N_sub_lg, iScale    
c            close(19) 
c            Stop
c         EndIf
         If(iDebug_lg)write (13,50) iLeft,iRight,mBound
 50      format(' iLeft/Right=',2i3,' Mesh Size=',3i3)

c                              Initialize 1-level grid
         mboundX =mBound(1)*2**(N_lev_lg-1)
         mboundY =mBound(2)*2**(N_lev_lg-1)
         mboundZ =mBound(3)*2**(N_lev_lg-1)
         if(max(mboundX,mboundY,mboundZ).gt.N_m1_lg) then
            write (13,*)' Error MakeMeshes: too small mesh size',
     &                      ' N_m1_lg=',N_m1_lg,
     &                      ' requested:',mboundX,mboundY,mboundZ
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
            accesstype = 'append'
            call Open_ASCII_FILE(19,name,accesstype)
            rewind 19
            j123=123
            write(19,*) j123,' Error MakeMeshes: too small mesh size',
     &                      ' N_m1_lg=',N_m1_lg,
     &                      ' requested:',mboundX,mboundY,mboundZ 
            close(19)
            CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
            STOP
         endif

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (m,i,j,k)
         Do m =1,8
            Do k=1,mboundZ
            Do j=1,mboundY
            Do i=1,mboundX
               ReshLG1(m,i,j,k) =0.
            EndDo
            EndDo
            EndDo
      EndDo
c                           1-level grid - get info from particles      
      Do i= 1,np_node    ! loop over primary particles only
         i_x = INT((x(i)-xLeft)*N_sub_lg)+1
         i_y = INT((y(i)-yLeft)*N_sub_lg)+1
         i_z = INT((z(i)-zLeft)*N_sub_lg)+1
            If(i_x.lt.1.or.i_y.lt.1.or.i_z.lt.1)Then
               write (13,*) ' Error in left boundary: ind=',
     &              i_x,i_y,i_z,' part=',i
               write (13,*) '      coords=',
     &              x(i),y(i),z(i)
            EndIf
            If(i_x.gt.mboundX.or.
     &         i_y.gt.mboundY.or.
     &         i_z.gt.mboundZ)Then
               write (13,*) ' Error in right boundary: ind=',
     &              i_x,i_y,i_z
               write (13,*) '      coords=',
     &              x(i),y(i),z(i)
               close(13)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
               STOP
            EndIf
            If(i_x.gt.N_m1_lg.or.
     &         i_y.gt.N_m1_lg.or.
     &         i_z.gt.N_m1_lg)Then
               write (13,*) ' Error in i_: ind=',
     &              i_x,i_y,i_z
               write (13,*) '      coords=',
     &              x(i),y(i),z(i)
               write (13,*) '      LeftBounds=',
     &              xLeft,yLeft,zLeft

               close(13)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
               STOP
            EndIf

         ReshLG1(1,i_x,i_y,i_z)=ReshLG1(1,i_x,i_y,i_z)+
     &                                          x(i)*Wpart(i)
         ReshLG1(2,i_x,i_y,i_z)=ReshLG1(2,i_x,i_y,i_z)+
     &                                          y(i)*Wpart(i)
         ReshLG1(3,i_x,i_y,i_z)=ReshLG1(3,i_x,i_y,i_z)+
     &                                          z(i)*Wpart(i)
         ReshLG1(4,i_x,i_y,i_z)=ReshLG1(4,i_x,i_y,i_z)+
     &                                          vx(i)*Wpart(i)
         ReshLG1(5,i_x,i_y,i_z)=ReshLG1(5,i_x,i_y,i_z)+
     &                                          vy(i)*Wpart(i)
         ReshLG1(6,i_x,i_y,i_z)=ReshLG1(6,i_x,i_y,i_z)+
     &                                          vz(i)*Wpart(i)
         ReshLG1(7,i_x,i_y,i_z)=ReshLG1(7,i_x,i_y,i_z)+
     &                                         pt(i)*Wpart(i)
         ReshLG1(8,i_x,i_y,i_z)=ReshLG1(8,i_x,i_y,i_z)+Wpart(i)
      EndDo    ! loop over particles

c                                              Second grid
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k,i2,j2,k2,m)
               Do m=1,8
            Do k=1,N_m2_lg
               k2 =2*k-1
            Do j=1,N_m2_lg
               j2 =2*j-1
            Do i=1,N_m2_lg
               i2 =2*i-1
                  ReshLG2(m,i,j,k) =
     &             ReshLG1(m,i2,j2,k2)+ReshLG1(m,i2+1,j2,k2)+  
     &             ReshLG1(m,i2,j2+1,k2)+ReshLG1(m,i2+1,j2+1,k2)+  
     &             ReshLG1(m,i2,j2,k2+1)+ReshLG1(m,i2+1 ,j2,k2+1)+  
     &             ReshLG1(m,i2,j2+1,k2+1)+ReshLG1(m,i2+1,j2+1,k2+1)  
               EndDo
            EndDo
            EndDo
            EndDo
c                                              Third grid
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k,i2,j2,k2,m)
                Do m=1,8
            Do k=1,N_m3_lg
               k2 =2*k-1
            Do j=1,N_m3_lg
               j2 =2*j-1
            Do i=1,N_m3_lg
               i2 =2*i-1
                  ReshLG3(m,i,j,k) =
     &             ReshLG2(m,i2,j2,k2)+ReshLG2(m,i2+1,j2,k2)+  
     &             ReshLG2(m,i2,j2+1,k2)+ReshLG2(m,i2+1,j2+1,k2)+  
     &             ReshLG2(m,i2,j2,k2+1)+ReshLG2(m,i2+1 ,j2,k2+1)+  
     &             ReshLG2(m,i2,j2+1,k2+1)+ReshLG2(m,i2+1,j2+1,k2+1)  
               EndDo
            EndDo
            EndDo
            EndDo
            If(N_lev_lg.eq.4)Then
c                                              Forth grid
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k,i2,j2,k2,m)
                Do m=1,8
            Do k=1,N_m4_lg
               k2 =2*k-1
            Do j=1,N_m4_lg
               j2 =2*j-1
            Do i=1,N_m4_lg
               i2 =2*i-1
                  ReshLG4(m,i,j,k) =
     &             ReshLG3(m,i2,j2,k2)+ReshLG3(m,i2+1,j2,k2)+  
     &             ReshLG3(m,i2,j2+1,k2)+ReshLG3(m,i2+1,j2+1,k2)+  
     &             ReshLG3(m,i2,j2,k2+1)+ReshLG3(m,i2+1 ,j2,k2+1)+  
     &             ReshLG3(m,i2,j2+1,k2+1)+ReshLG3(m,i2+1,j2+1,k2+1)  
               EndDo
            EndDo
            EndDo
            EndDo
         EndIf
c                                        Normalize all grids
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k,m)
      Do m=1,7
            Do k=1,N_m1_lg
            Do j=1,N_m1_lg
            Do i=1,N_m1_lg
               If(ReshLG1(8,i,j,k).gt.1.e-10)
     &          ReshLG1(m,i,j,k) =ReshLG1(m,i,j,k) /
     &                                      ReshLG1(8,i,j,k)
            EndDo
            EndDo
            EndDo
            Do k=1,N_m2_lg
            Do j=1,N_m2_lg
            Do i=1,N_m2_lg
              If(ReshLG2(8,i,j,k).gt.1.e-10)
     &          ReshLG2(m,i,j,k) =ReshLG2(m,i,j,k) /
     &                                      ReshLG2(8,i,j,k)
            EndDo
            EndDo
            EndDo
            Do k=1,N_m3_lg
            Do j=1,N_m3_lg
            Do i=1,N_m3_lg
              If(ReshLG3(8,i,j,k).gt.1.e-10)
     &          ReshLG3(m,i,j,k) =ReshLG3(m,i,j,k) /
     &                                      ReshLG3(8,i,j,k)
            EndDo
            EndDo
            EndDo
            If(N_lev_lg.eq.4)Then
            Do k=1,N_m4_lg
            Do j=1,N_m4_lg
            Do i=1,N_m4_lg
              If(ReshLG4(8,i,j,k).gt.1.e-10)
     &          ReshLG4(m,i,j,k) =ReshLG4(m,i,j,k) /
     &                                      ReshLG4(8,i,j,k)
            EndDo
            EndDo
            EndDo
            EndIf   ! end 4-th level
         EndDo       ! end loop over m=1,6
c                             do not allow coordinates large than ng+1
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (i,j,k,m)
      Do m=1,3
            Do k=1,N_m1_lg
            Do j=1,N_m1_lg
            Do i=1,N_m1_lg
               ReshLG1(m,i,j,k) =min(ReshLG1(m,i,j,k),xmax_b)
            EndDo
            EndDo
            EndDo
            Do k=1,N_m2_lg
            Do j=1,N_m2_lg
            Do i=1,N_m2_lg
               ReshLG2(m,i,j,k) =min(ReshLG2(m,i,j,k),xmax_b)
            EndDo
            EndDo
            EndDo
            Do k=1,N_m3_lg
            Do j=1,N_m3_lg
            Do i=1,N_m3_lg
               ReshLG3(m,i,j,k) =min(ReshLG3(m,i,j,k),xmax_b)
            EndDo
            EndDo
            EndDo
            If(N_lev_lg.eq.4)Then
            Do k=1,N_m4_lg
            Do j=1,N_m4_lg
            Do i=1,N_m4_lg
               ReshLG4(m,i,j,k) =min(ReshLG4(m,i,j,k),xmax_b)
            EndDo
            EndDo
            EndDo
            EndIf   ! end 4-th level
         EndDo       ! end loop over m=1,3
      Return 
      End
c-------------------------------------------------------------------- 
      SUBROUTINE MakeFake(iN0,nn)
c-------------------------------------------------------------------- 
c          
c     purpose :   create fake particles for tests
c     input:         iN0 - node         
c               nn - number of particles
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      DIMENSION iN0(3)
      
      nn = 0
         iLeft       =  l_divX(iN0(1))               ! x boundary
         jLeft       =  l_divY(iN0(1),iN0(2))               ! y
         kLeft       =  l_divZ(iN0(1),iN0(2),iN0(3))   ! z
         iRight      =  l_divX(iN0(1)+1)               ! x boundary
         jRight      =  l_divY(iN0(1),iN0(2)+1)               ! y
         kRight      =  l_divZ(iN0(1),iN0(2),iN0(3)+1)   ! z

         ii =0
         jj = 0
         kk = 0
      Do k=kLeft,kRight-1
      Do j=jLeft,jRight-1
      Do i=iLeft,iRight-1
c         do kk =0,1
c         do jj =0,1
c         do ii =0,1
         nn = nn +1
            x(nn) = i+ii/2.+0.0001
            y(nn) = j+jj/2.+0.0001
            z(nn) = k+kk/2.+0.0001
            vx(nn) = i+0.0001
            vy(nn) = j+0.0001
            vz(nn) = k+0.0001
            Wpart(nn) = 1.
            i_par(nn) = nn
            pt(nn) = 0.001
c         EndDo 
c         EndDo 
c         EndDo 
      EndDo
      EndDo
      EndDo

      n_all = nn
      n_recv_prim = nn
      write (13,*) ' Number of fake particles =',n_all
      Return
      End
c-------------------------------------------------------------------- 
      SUBROUTINE Find_Small(iN0,Nbound,np_node,
     &                                           nn,ncount)
c-------------------------------------------------------------------- 
c
c  purpose: find small (primary) particles, which node iN0 has (will send)
c                  for node iN1. These are buffer particles. They are
c                  used to increase the volume around the node for
c                  refinement.
c     input: iN0(3) - coordinates of sending node
c            Nbound - number of boundaries
c            isBound(8,6) - boundaries (units of zero-level mesh)
c            np_node - number of primary particles in iN0 node
c     output:writes the small particles at the end of *_se buffer
c            ncount - number of small particles to be send by iN0 node
c            nn   - current number of particles in the buffer *_se
c    conditions for taking a particle:
c            isBound(*,1)=< X(i) < isBound(*,2), 
c                                  Y: (*,3-4), Z: (*,5-6)
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      DIMENSION iN0(3)
      Logical   inside,in_x,in_y,in_z

      ncount =0
      If(Nbound.le.0)RETURN  ! no boundaries or the same node 
         
      Do i= 1,np_node    ! loop over primary particles only
         i_x = x(i)              ! index of the particle in global mesh
         i_y = y(i)
         i_z = z(i)
         inside =.false. 
         Do jb =1,Nbound ! check if particle is inside boundaries
            If(.NOT.inside)Then     ! do it only until first boundary
               in_x = (i_x.ge.isBound(jb,1).and.i_x.lt.isBound(jb,2))
               in_y = (i_y.ge.isBound(jb,3).and.i_y.lt.isBound(jb,4))
               in_z = (i_z.ge.isBound(jb,5).and.i_z.lt.isBound(jb,6))
               If(in_x.and.in_y.and.in_z)Then ! take the particle, place in buffer
                  ncount = ncount +1
                  nn     = nn     +1 
                  inside = .true. 
                  if(nn.gt.max_send) then 
                    write (13,*)' Error Find_Small:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                    call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=123
                    write(19,*) j123,' Error Find_Small:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif
  
                  x_se(nn)   =  x(i)    
                  y_se(nn)   =  y(i)    
                  z_se(nn)   =  z(i)    
                  vx_se(nn)  =  vx(i)   
                  vy_se(nn)  =  vy(i)    
                  vz_se(nn)  =  vz(i)   
                  pt_se(nn)  =  pt(i)  
                  wpar_se(nn)=  wpart(i)   
                  ip_se(nn)  =  i_par(i)   
               EndIf
            EndIf             ! end .not.inside condition 
         EndDo              ! end boundary loop
      EndDo                     ! end particles loop
      RETURN
      END
c-------------------------------------------------------------------- 
      SUBROUTINE Find_Large3(iN0,Nbound,nn,ncount)
c-------------------------------------------------------------------- 
c
c  purpose: find large particles, which node iN0 creates (will send)
c                  for node iN1
c     input: iN0(3) - coordinates of sending node
c            nn   - current number of particles in the buffer *_se
c     output:writes the large particles at the end of *_se buffer
c             ncount - number of particles produced by this node 
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      DIMENSION iN0(3),imBound(8,6)
      Logical   inside,in_x,in_y,in_z
      REAL*8          summW,wcount1,wcount2,wcount3,wcount4

      If(Nbound.eq.-1)Then  ! the same node
         ncount =0
         RETURN 
      EndIf 
      summW =0.
      If(Nbound.ne.0)Then   ! scale boundaries to ReshG1 mesh
         Do idir =1,3 
            If(idir.eq.1)Then
                i0Left       =  l_divX(iN0(1))
              else if (idir.eq.2) Then
                i0Left       =  l_divY(iN0(1),iN0(2))
              else
                i0Left       =  l_divZ(iN0(1),iN0(2),iN0(3))
            endif
            Do jb =1,Nbound
               imBound(jb,idir*2-1) = (isBound(jb,idir*2-1)-i0Left )
     &                                                  *N_sub_lg+1
               imBound(jb,idir*2)   = (isBound(jb,idir*2)-i0Left)
     &                                                  *N_sub_lg+1
            EndDo     ! end boundaries
         EndDo     ! end directions
         If(iDebug_lg)Then
            write (13,10) Nbound,((imBound(jb,i),i=1,6),jb=1,Nbound)
            write (13,20) ((isBound(jb,i),i=1,6),jb=1,Nbound)
 10      format(' Scale boundaries: Nbound=',i3,/
     &                (6i4))
 20      format(' Old boundaries: ', (6i4))
         EndIf         ! end debug mode
      EndIf            ! end boundary scaling

      ncount =0
      ncount1 =0
      ncount2 =0
      ncount3 =0
      wcount1 =0
      wcount2 =0
      wcount3 =0
         Do k= 1,mBound(3)    ! check every cell of the map
         Do j= 1,mBound(2)
         Do i= 1,mBound(1)
            If(MeshLG(i,j,k).eq.1)Then       ! switch for different meshes
               Do kk =k*4-3,k*4
               Do jj =j*4-3,j*4
               Do ii =i*4-3,i*4
                  If( ReshLG1(8,ii,jj,kk).gt.small_lg)Then ! nonempty cell 
                  inside =.false. 
                  Do jb =1,Nbound              
                    in_x = (ii.ge.imBound(jb,1).and.ii.lt.imBound(jb,2))
                    in_y = (jj.ge.imBound(jb,3).and.jj.lt.imBound(jb,4))
                    in_z = (kk.ge.imBound(jb,5).and.kk.lt.imBound(jb,6))
                    If(in_x.and.in_y.and.in_z)inside =.true. 
                  EndDo 
                  If(.NOT.inside)Then     ! take only cells which are outside
                                                        ! of boundaries. Inside cells are
                                                        ! presented by primary small particles 
                  ncount = ncount +1
                  ncount1= ncount1+1
                  nn     = nn     +1 
                  if(nn.gt.max_send)then
                    write (13,*)' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                    call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=123
                    write(19,*) j123,' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    write(19,*) j123,' mBound array=',
     &                 mBound(1),mBound(2),mBound(3)
                    write(19,*) j123,' estimate:',
     &                 mBound(1)*mBound(2)*mBound(3)*4**3 
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif

                  x_se(nn)   =  ReshLG1(1,ii,jj,kk)    
                  y_se(nn)   =  ReshLG1(2,ii,jj,kk) 
                  z_se(nn)   =  ReshLG1(3,ii,jj,kk) 
                  vx_se(nn)  =  ReshLG1(4,ii,jj,kk) 
                  vy_se(nn)  =  ReshLG1(5,ii,jj,kk) 
                  vz_se(nn)  =  ReshLG1(6,ii,jj,kk) 
                  pt_se(nn)  =  ReshLG1(7,ii,jj,kk) 
                  wpar_se(nn)=  ReshLG1(8,ii,jj,kk)
                  ip_se(nn)  =  0
                  summW = summW +ReshLG1(8,ii,jj,kk)
                  wcount1 = wcount1 +ReshLG1(8,ii,jj,kk)
                  EndIf       ! end not-inside boundaries
                  EndIf       ! end non-empty cell test
               EndDo 
               EndDo 
               EndDo 

            else If(MeshLG(i,j,k).eq.2) Then ! 2nd level
               Do kk =k*2-1,k*2
               Do jj =j*2-1,j*2
               Do ii =i*2-1,i*2
                  If( ReshLG2(8,ii,jj,kk).gt.small_lg)Then ! nonempty cell
                  ncount = ncount +1
                  ncount2= ncount2+1
                  nn     = nn     +1 
                  if(nn.gt.max_send)then
                    write (13,*)' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                    call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=123
                    write(19,*) j123,' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif

                  x_se(nn)   =  ReshLG2(1,ii,jj,kk)    
                  y_se(nn)   =  ReshLG2(2,ii,jj,kk) 
                  z_se(nn)   =  ReshLG2(3,ii,jj,kk) 
                  vx_se(nn)  =  ReshLG2(4,ii,jj,kk) 
                  vy_se(nn)  =  ReshLG2(5,ii,jj,kk) 
                  vz_se(nn)  =  ReshLG2(6,ii,jj,kk) 
                  pt_se(nn)  =  ReshLG2(7,ii,jj,kk) 
                  wpar_se(nn)=  ReshLG2(8,ii,jj,kk)
                  ip_se(nn)  =  0
                  summW = summW +ReshLG2(8,ii,jj,kk)
                  wcount2 = wcount2 +ReshLG2(8,ii,jj,kk)
                  EndIf 
               EndDo 
               EndDo 
               EndDo 
               
            else                             ! 3th level
               If( ReshLG3(8,i,j,k).gt.small_lg)Then ! nonempty cell
                  ncount = ncount +1
                  ncount3= ncount3+1
                  nn     = nn     +1 
                   if(nn.gt.max_send)then 
                    write (13,*)' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                    call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=122
                    write(19,*) j123,' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    write(19,*) j123,' mBound array=',
     &                 mBound(1),mBound(2),mBound(3)
                    write(19,*) j123,' estimate:',
     &                 mBound(1)*mBound(2)*mBound(3)
                    write(19,*) j123,' counts=',
     &                  ncount1,ncount2,ncount3                             
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif

                  x_se(nn)   =  ReshLG3(1,i,j,k)    
                  y_se(nn)   =  ReshLG3(2,i,j,k) 
                  z_se(nn)   =  ReshLG3(3,i,j,k) 
                  vx_se(nn)  =  ReshLG3(4,i,j,k) 
                  vy_se(nn)  =  ReshLG3(5,i,j,k) 
                  vz_se(nn)  =  ReshLG3(6,i,j,k) 
                  pt_se(nn)  =  ReshLG3(7,i,j,k) 
                  wpar_se(nn)=  ReshLG3(8,i,j,k)
                  ip_se(nn)  =  0
                  summW = summW +ReshLG3(8,i,j,k)
                  wcount3 = wcount3 +ReshLG3(8,i,j,k)
               EndIf 
            EndIf 
         EndDo
         EndDo
         EndDo

         If(iDebug_lg)write (13,*)  ' Nparticle created=',
     &    ncount1,ncount2,ncount3
         If(iDebug_lg)write (13,*)  ' Total weight =',summW,
     &    wcount1,wcount2,wcount3
         RETURN 
         End
c-------------------------------------------------------------------- 
      SUBROUTINE Find_Large4(iN0,Nbound,nn,ncount)
c-------------------------------------------------------------------- 
c
c  purpose: find large particles, which node iN0 creates (will send)
c                  for node iN1
c     input: iN0(3) - coordinates of sending node
c            nn   - current number of particles in the buffer *_se
c     output:writes the large particles at the end of *_se buffer
c             ncount - number of particles produced by this node 
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      DIMENSION iN0(3),imBound(8,6)
      Logical   inside,in_x,in_y,in_z
      REAL*8          summW,wcount1,wcount2,wcount3,wcount4

      If(Nbound.eq.-1)Then  ! the same node
         ncount =0
         RETURN 
      EndIf 
      summW =0.
      If(Nbound.ne.0)Then   ! scale boundaries to ReshG1 mesh
         Do idir =1,3                                        
            If(idir.eq.1)Then
                i0Left       =  l_divX(iN0(1))
              else if (idir.eq.2) Then
                i0Left       =  l_divY(iN0(1),iN0(2))
              else
                i0Left       =  l_divZ(iN0(1),iN0(2),iN0(3))
            endif
            Do jb =1,Nbound
               imBound(jb,idir*2-1) = (isBound(jb,idir*2-1)-i0Left )
     &                                                  *N_sub_lg+1
               imBound(jb,idir*2)   = (isBound(jb,idir*2)-i0Left)
     &                                                  *N_sub_lg+1
            EndDo     ! end boundaries
         EndDo     ! end directions
         If(iDebug_lg)Then
            write (13,10) Nbound,((imBound(jb,i),i=1,6),jb=1,Nbound)
            write (13,20) ((isBound(jb,i),i=1,6),jb=1,Nbound)
 10      format(' Scale boundaries: Nbound=',i3,/
     &                (6i4))
 20      format(' Old boundaries: ', (6i4))
         EndIf         ! end debug mode
      EndIf            ! end boundary scaling

      ncount =0
      ncount1 =0
      ncount2 =0
      ncount3 =0
      ncount4 =0
      wcount1 =0
      wcount2 =0
      wcount3 =0
      wcount4 =0
         Do k= 1,mBound(3)    ! check every cell of the map
         Do j= 1,mBound(2)
         Do i= 1,mBound(1)
            If(MeshLG(i,j,k).eq.1)Then       ! switch for different meshes
               Do kk =k*8-7,k*8
               Do jj =j*8-7,j*8
               Do ii =i*8-7,i*8
                  If( ReshLG1(8,ii,jj,kk).gt.small_lg)Then ! nonempty cell 
                  inside =.false. 
                  Do jb =1,Nbound              
                    in_x = (ii.ge.imBound(jb,1).and.ii.lt.imBound(jb,2))
                    in_y = (jj.ge.imBound(jb,3).and.jj.lt.imBound(jb,4))
                    in_z = (kk.ge.imBound(jb,5).and.kk.lt.imBound(jb,6))
                    If(in_x.and.in_y.and.in_z)inside =.true. 
                  EndDo 
                  If(.NOT.inside)Then     ! take only cells which are outside
                                                        ! of boundaries. Inside cells are
                                                        ! presented by primary small particles 
                  ncount = ncount +1
                  ncount1= ncount1+1
                  nn     = nn     +1 
                  if(nn.gt.max_send)then
                    write (13,*)' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                   call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=125
                    write(19,*) j123,' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif

                  x_se(nn)   =  ReshLG1(1,ii,jj,kk)    
                  y_se(nn)   =  ReshLG1(2,ii,jj,kk) 
                  z_se(nn)   =  ReshLG1(3,ii,jj,kk) 
                  vx_se(nn)  =  ReshLG1(4,ii,jj,kk) 
                  vy_se(nn)  =  ReshLG1(5,ii,jj,kk) 
                  vz_se(nn)  =  ReshLG1(6,ii,jj,kk) 
                  pt_se(nn)  =  ReshLG1(7,ii,jj,kk) 
                  wpar_se(nn)=  ReshLG1(8,ii,jj,kk)
                  ip_se(nn)  =  0
                  summW = summW +ReshLG1(8,ii,jj,kk)
                  wcount1 = wcount1 +ReshLG1(8,ii,jj,kk)
                  EndIf       ! end not-inside boundaries
                  EndIf       ! end non-empty cell test
               EndDo 
               EndDo 
               EndDo 

            else If(MeshLG(i,j,k).eq.2) Then ! 2nd level
               Do kk =k*4-3,k*4
               Do jj =j*4-3,j*4
               Do ii =i*4-3,i*4
                  If( ReshLG2(8,ii,jj,kk).gt.small_lg)Then ! nonempty cell
                  ncount = ncount +1
                  ncount2= ncount2+1
                  nn     = nn     +1 
                  if(nn.gt.max_send)then
                    write (13,*)' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                    call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=126
                    write(19,*) j123,' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif

                  x_se(nn)   =  ReshLG2(1,ii,jj,kk)    
                  y_se(nn)   =  ReshLG2(2,ii,jj,kk) 
                  z_se(nn)   =  ReshLG2(3,ii,jj,kk) 
                  vx_se(nn)  =  ReshLG2(4,ii,jj,kk) 
                  vy_se(nn)  =  ReshLG2(5,ii,jj,kk) 
                  vz_se(nn)  =  ReshLG2(6,ii,jj,kk) 
                  pt_se(nn)  =  ReshLG2(7,ii,jj,kk) 
                  wpar_se(nn)=  ReshLG2(8,ii,jj,kk)
                  ip_se(nn)  =  0
                  summW = summW +ReshLG2(8,ii,jj,kk)
                  wcount2 = wcount2 +ReshLG2(8,ii,jj,kk)
                  EndIf 
               EndDo 
               EndDo 
               EndDo 
               
            else If(MeshLG(i,j,k).eq.3) Then ! 3rd level
               Do kk =k*2-1,k*2
               Do jj =j*2-1,j*2
               Do ii =i*2-1,i*2
               If( ReshLG3(8,ii,jj,kk).gt.small_lg)Then ! nonempty cell
                  ncount = ncount +1
                  ncount3= ncount3+1
                  nn     = nn     +1 
                  if(nn.gt.max_send)then
                    write (13,*)' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                    call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=127
                    write(19,*) j123,' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif
                  x_se(nn)   =  ReshLG3(1,ii,jj,kk)    
                  y_se(nn)   =  ReshLG3(2,ii,jj,kk) 
                  z_se(nn)   =  ReshLG3(3,ii,jj,kk) 
                  vx_se(nn)  =  ReshLG3(4,ii,jj,kk) 
                  vy_se(nn)  =  ReshLG3(5,ii,jj,kk) 
                  vz_se(nn)  =  ReshLG3(6,ii,jj,kk) 
                  pt_se(nn)  =  ReshLG3(7,ii,jj,kk) 
                  wpar_se(nn)=  ReshLG3(8,ii,jj,kk)
                  ip_se(nn)  =  0
                  summW = summW +ReshLG3(8,ii,jj,kk)
                  wcount3 = wcount3 +ReshLG3(8,ii,jj,kk)
                  EndIf 
               EndDo 
               EndDo 
               EndDo 
 
            else                             ! 4th level
               If( ReshLG4(8,i,j,k).gt.small_lg)Then ! nonempty cell
                  ncount = ncount +1
                  ncount4= ncount4+1
                  nn     = nn     +1 
                   if(nn.gt.max_send)then 
                    write (13,*)' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    name = 'ERROR_RUN'
            write(name,'(a,i4.4,a)') 'ERROR_RUN',irank+1,'.DAT'
                    accesstype = 'append'
                    call Open_ASCII_FILE(19,name,accesstype)
                    rewind 19
                    j123=128
                    write(19,*) j123,' Error Find_Large:',
     &              ' buffer _se: ',nn,'  maximum is ',max_send  
                    close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2) 
                    STOP
                  endif

                  x_se(nn)   =  ReshLG4(1,i,j,k)    
                  y_se(nn)   =  ReshLG4(2,i,j,k) 
                  z_se(nn)   =  ReshLG4(3,i,j,k) 
                  vx_se(nn)  =  ReshLG4(4,i,j,k) 
                  vy_se(nn)  =  ReshLG4(5,i,j,k) 
                  vz_se(nn)  =  ReshLG4(6,i,j,k) 
                  pt_se(nn)  =  ReshLG4(7,i,j,k) 
                  wpar_se(nn)=  ReshLG4(8,i,j,k)
                  ip_se(nn)  =  0
                  summW = summW +ReshLG4(8,i,j,k)
                  wcount4 = wcount4 +ReshLG4(8,i,j,k)
        If(iDebug_lg.and.(x_se(nn)<1..or.x_se(nn)>129.))
     &  WRITE(13,'(" Large4:wrong X:",3g12.4,4i12)')
     &           x_se(nn),y_se(nn),z_se(nn),
     &         nn,k,node

               EndIf 
            EndIf 
         EndDo
         EndDo
         EndDo

         If(iDebug_lg)write (13,*)  ' Nparticle created=',
     &    ncount1,ncount2,ncount3,ncount4
         If(iDebug_lg)write (13,*)  ' Total weight =',summW,
     &    wcount1,wcount2,wcount3,wcount4
         RETURN 
         End

c-------------------------------------------------------------------- 
      SUBROUTINE Test_Part(node)
c-------------------------------------------------------------------- 
c
c  purpose: check coords and vels of particles
c
c     input: 

      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      Real*8  Smass,smass0,smass1,smass2
c      Real*4  xmin(7,3),xmax(7,3)
      Real*8  xmin(7,3),xmax(7,3),wwpart

      np_node = np_sb(node)
      Smass =0.
      smass0 =0.
      smass1 =0.
      smass2 =0.
      Do i=1,7
         Do j=1,3
            xmin(i,j) =1.e+10
            xmax(i,j) =-1.e+10
         EndDo 
      EndDo 

      do i=1,n_all
         Smass =Smass + wpart(i)
         wwpart = wpart(i)
         If(x(i).ge.xn.or.x(i).lt.1.d-0)
     &     write(13,*)i,x(i),y(i),z(i),xn
         If(y(i).ge.xn.or.y(i).lt.1.d-0)
     &     write(13,*)i,x(i),y(i),z(i),xn
         If(z(i).ge.xn.or.z(i).lt.1.d-0)
     &     write(13,*)i,x(i),y(i),z(i),xn

         If(i.le.np_node)Then
            smass0=smass0 + wpart(i)
            xmin(1,1) =min(x(i),xmin(1,1))
            xmin(2,1) =min(y(i),xmin(2,1))
            xmin(3,1) =min(z(i),xmin(3,1))
            xmin(4,1) =min(vx(i),xmin(4,1))
            xmin(5,1) =min(vy(i),xmin(5,1))
            xmin(6,1) =min(vz(i),xmin(6,1))
            xmax(1,1) =max(x(i),xmax(1,1))
            xmax(2,1) =max(y(i),xmax(2,1))
            xmax(3,1) =max(z(i),xmax(3,1))
            xmax(4,1) =max(vx(i),xmax(4,1))
            xmax(5,1) =max(vy(i),xmax(5,1))
            xmax(6,1) =max(vz(i),xmax(6,1))
c            xmin(7,1) =min(xmin(7,1),wpart(i)
c            xmax(7,1) =max(xmax(7,1),wpart(i))
            xmin(7,1) =min(xmin(7,1),wwpart)
            xmax(7,1) =max(xmax(7,1),wwpart)
         Else If(i.le.n_refin)Then
            smass1=smass1 + wpart(i)
            xmin(1,2) =min(x(i),xmin(1,2))
            xmin(2,2) =min(y(i),xmin(2,2))
            xmin(3,2) =min(z(i),xmin(3,2))
            xmin(4,2) =min(vx(i),xmin(4,2))
            xmin(5,2) =min(vy(i),xmin(5,2))
            xmin(6,2) =min(vz(i),xmin(6,2))
            xmax(1,2) =max(x(i),xmax(1,2))
            xmax(2,2) =max(y(i),xmax(2,2))
            xmax(3,2) =max(z(i),xmax(3,2))
            xmax(4,2) =max(vx(i),xmax(4,2))
            xmax(5,2) =max(vy(i),xmax(5,2))
            xmax(6,2) =max(vz(i),xmax(6,2))
c            xmin(7,2) =min(xmin(7,2),wpart(i))
c            xmax(7,2) =max(xmax(7,2),wpart(i))
            xmin(7,2) =min(xmin(7,2),wwpart)
            xmax(7,2) =max(xmax(7,2),wwpart)
         Else
            smass2=smass2 + wpart(i)
            xmin(1,3) =min(x(i),xmin(1,3))
            xmin(2,3) =min(y(i),xmin(2,3))
            xmin(3,3) =min(z(i),xmin(3,3))
            xmin(4,3) =min(vx(i),xmin(4,3))
            xmin(5,3) =min(vy(i),xmin(5,3))
            xmin(6,3) =min(vz(i),xmin(6,3))
            xmax(1,3) =max(x(i),xmax(1,3))
            xmax(2,3) =max(y(i),xmax(2,3))
            xmax(3,3) =max(z(i),xmax(3,3))
            xmax(4,3) =max(vx(i),xmax(4,3))
            xmax(5,3) =max(vy(i),xmax(5,3))
            xmax(6,3) =max(vz(i),xmax(6,3))
            xmin(7,3) =min(xmin(7,3),wwpart)
            xmax(7,3) =max(xmax(7,3),wwpart)
         EndIf
      EndDo 
         write(13,100)n_all,n_refin,np_node,Smass,node
 100     format(' Test: N_all=',i10,' N_small=',i10,
     &           ' N_primary=',i7,' Total weight=',g14.6,i4)
         write(13,110)node,smass0,(xmin(i,1),i=1,7),
     &                                (xmax(i,1),i=1,7)
         write(13,120)node,smass1,(xmin(i,2),i=1,7),
     &                                (xmax(i,2),i=1,7)
         write(13,130)node,smass2,(xmin(i,3),i=1,7),
     &                                (xmax(i,3),i=1,7)
 110     format(' Primary: mass =',i4,g12.5,' Min=',7g10.3,/
     &    40x,'Max=',7g10.3)
 120     format(' Buffer : mass =',i4,g12.5,' Min=',7g10.3,/
     &    40x,'Max=',7g10.3)
 130     format(' Large  : mass =',i4,g12.5,' Min=',7g10.3,/
     &    40x,'Max=',7g10.3)
         RETURN 
         End
c-------------------------------------------------------------------- 
      SUBROUTINE Send_Small()
ccc      SUBROUTINE Send_Small(node)
c-------------------------------------------------------------------- 
c
c  purpose: gathers and sends small particles
c
c     input: iN0(3) - coordinates of sending node
c           
c     output: sends particles, sets n_refin =number of particles
c                                                                 (primaries + small)
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      double precision xx,yy,zz,xLf,xRt,yLf,yRt,zLf,zRt,
     &             xLf1,xRt1,yLf1,yRt1,zLf1,zRt1       
      INTEGER iN0(3),iN1(3)
      Logical Lout,Loutx,Louty,Loutz

      node = irank + 1
       write(13,*) ' SendSmall node=',node                                                                                      
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!        CALL mpi_finalize(ierr)
!        STOP  

      CALL Node_to_IJK(node,iN0)
c      Node_my = iN0(1) +(iN0(2)-1)*n_div
c     &                             +(iN0(3)-1)*n_div**2
      np_node = np_sb(node)
      nn = 0
!       write(13,*) ' SendSmall node=',node       
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!        CALL mpi_finalize(ierr)                   
!        STOP   
      Do kn =1,n_divz      !  Loop over other nodes
         iN1(3) = kn
      Do jn =1,n_divy
         iN1(2) = jn
      Do in =1,n_divx
         iN1(1) = in
         node1 = in +(jn-1)*n_divx+(kn-1)*n_divx*n_divy
         sdispl(node1)    =  nn 
         CALL BoundNode(iN0,iN1,Nbound)
         CALL Find_Small(iN0,Nbound,np_node,nn,ncount)
         sendcount(node1) = ncount
         If(iDebug_lg)Then
            write (13,20) Nbound,iN0,iN1,nn,ncount
 20      format(/' =======   Nbound=',i3,' Nodes=',2(3i3,2x)/
     &               5x,' global counter=',i8,' node counter=',i6)
         If(Nbound.gt.0)Then
            write(13,*) '     iBound: large particles boundaries'
            write(13,10)((iBound(i,j),j=1,6),i,i=1,Nbound)
            write(13,*) '     isBound: small particles boundaries'
            write(13,10)((isBound(i,j),j=1,6),i,i=1,Nbound)
         EndIf
 10      format(3(3x,2i4),' boundary=',i3)
         EndIf           ! end debug mode  
      EndDo 
      EndDo 
      EndDo 

!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

      CALL Send_Receive()

c                       construct working arrays for ART from primary 
c                      (belonging to  the node) and secondary (sent) data
      i_rec = 0
      DO k = 1, n_nodes
        i_rec = i_rec + recvcount(k)
      ENDDO
      If(iDebug_lg)WRITE(13,*)i_rec, ' primaries added to ',
     &     np_sb(node)       
      n_recv_prim = i_rec               ! received primaries 
c                 calculate  boundaries of current node +/- width of buffer
c                 for small particles, which define refinement
c                 The width of the buffer is dBuffer      
           xLf      =  l_divX(iN0(1))                       -dBuffer
           xRt     =  l_divX(iN0(1)+1)                   +dBuffer
           yLf      =  l_divY(iN0(1),iN0(2))             -dBuffer
           yRt     =  l_divY(iN0(1),iN0(2)+1)         +dBuffer
           zLf      =  l_divZ(iN0(1),iN0(2),iN0(3))   -dBuffer
           zRt     =  l_divZ(iN0(1),iN0(2),iN0(3)+1) +dBuffer
           xLf1     = xLf +ng
           yLf1     = yLf +ng
           zLf1     = zLf +ng
           xRt1     = xRt -ng
           yRt1     = yRt -ng
           zRt1     = zRt -ng
      If(iDebug_lg)WRITE(13,*)' boundaries for refinement particles'
      If(iDebug_lg)WRITE(13,'(4f8.2)')xRt1,xLf,xRt,xLf1
      If(iDebug_lg)WRITE(13,'(4f8.2)')yRt1,yLf,yRt,yLf1
      If(iDebug_lg)WRITE(13,'(4f8.2)')zRt1,zLf,zRt,zLf1

c                Put particles, which are  inside the buffer
c                these and primary particles define refinement
      icount = 0
      DO k = 1, i_rec
         xx = x_re(k)
         yy = y_re(k)
         zz = z_re(k)
         Loutx=(xx.gt.xRt1.and.xx.lt.xLf).or.(xx.gt.xRt.and.xx.lt.xLf1)
         Louty=(yy.gt.yRt1.and.yy.lt.yLf).or.(yy.gt.yRt.and.yy.lt.yLf1)
         Loutz=(zz.gt.zRt1.and.zz.lt.zLf).or.(zz.gt.zRt.and.zz.lt.zLf1)
         Lout = Loutx.or.Louty.or.Loutz
         If(.not.Lout)Then
!      If(iDebug_lg)WRITE(13,'(" in:",3g12.4,4i3)')xx,yy,zz,
!     &                                           Loutx,Louty,Loutz,Lout
           icount = icount +1
           i = np_sb(node) + icount
           x(i)   =  xx
           y(i)   =  yy
           z(i)   =  zz
           vx(i)  =  vx_re(k)
           vy(i)  =  vy_re(k)
           vz(i)  =  vz_re(k)
           pt(i)  =  pt_re(k)
           wpart(i) =  wpar_re(k)
           i_par(i) =  ip_re(k)
        EndIf 
      ENDDO
       n_refin = icount + np_sb(node) 
      If(iDebug_lg)WRITE(13,*)' Particles for possible refinement: ',
     +              n_refin

c                Add particles, which are  outside the buffer
c                these particles do not define refinements
      DO k = 1, i_rec
         xx = x_re(k)
         yy = y_re(k)
         zz = z_re(k)
         Loutx=(xx.gt.xRt1.and.xx.lt.xLf).or.(xx.gt.xRt.and.xx.lt.xLf1)
         Louty=(yy.gt.yRt1.and.yy.lt.yLf).or.(yy.gt.yRt.and.yy.lt.yLf1)
         Loutz=(zz.gt.zRt1.and.zz.lt.zLf).or.(zz.gt.zRt.and.zz.lt.zLf1)
         Lout = Loutx.or.Louty.or.Loutz
         If(Lout)Then
!      If(iDebug_lg)WRITE(13,'(" out:",3g12.4,4i3)')xx,yy,zz,
!     &                                           Loutx,Louty,Loutz,Lout

           icount = icount +1
           i = np_sb(node) + icount
           x(i)   =  xx
           y(i)   =  yy
           z(i)   =  zz
           vx(i)  =  vx_re(k)
           vy(i)  =  vy_re(k)
           vz(i)  =  vz_re(k)
           pt(i)  =  pt_re(k)
           wpart(i) =  wpar_re(k)
           i_par(i) =  ip_re(k)
           If(iDebug_lg.and.(x(i)<1..or.x(i)>129.))
     &          WRITE(13,'(" wrong X:",3g12.4,4i12)')xx,yy,zz,i,node
        EndIf 
      ENDDO
       n_all = icount + np_sb(node) 
c                     n_all is the current number of particles in the node
c

      RETURN 
      End
c-------------------------------------------------------------------- 
      SUBROUTINE Send_Large()
ccc      SUBROUTINE Send_Large(node)
c-------------------------------------------------------------------- 
c
c  purpose: creates and sends large particles,
c                  which node iN0 creates 
c                  for all other nodes
c     input: iN0(3) - coordinates of sending node
c           
c     output: sends particles, sets n_all
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'
      INTEGER iN0(3),iN1(3)

      node = irank + 1 
      CALL Node_to_IJK(node,iN0)
c      Node_my = iN0(1) +(iN0(2)-1)*n_div
c     &                             +(iN0(3)-1)*n_div**2
      CALL MakeMeshes(iN0,np_sb(node))

      nn = 0                    !  initiate global counter
      Do kn =1,n_divz      !  Loop over other nodes
         iN1(3) = kn
      Do jn =1,n_divy
         iN1(2) = jn
      Do in =1,n_divx
         iN1(1) = in
         node1 = in +(jn-1)*n_divx+(kn-1)*n_divx*n_divy
         sdispl(node1)    =  nn 
         CALL BoundNode(iN0,iN1,Nbound)
         CALL MakeMapNode(Nbound)
c         CALL TESTMESHES
         If(N_lev_lg.eq.4)Then
            CALL Find_Large4(iN0,Nbound,nn,ncount)
         Else
            CALL Find_Large3(iN0,Nbound,nn,ncount)
         EndIf
         sendcount(node1) = ncount
         If(iDebug_lg)Then         
            write (13,20) Nbound,iN0,iN1,nn,ncount
 20         format(/' =======   Nbound=',i3,' Nodes=',2(3i3,2x)/
     &                  5x,' global counter=',i8,' node counter=',i6)
            If(Nbound.gt.0)Then
               write(13,*) '     iBound: large particles boundaries'
               write(13,10)((iBound(i,j),j=1,6),i,i=1,Nbound)
               write(13,*) '     isBound: small particles boundaries'
               write(13,10)((isBound(i,j),j=1,6),i,i=1,Nbound)
            EndIf
         EndIf                ! end debug
      EndDo 
      EndDo 
      EndDo 
 10      format(3(3x,2i4),' boundary=',i3)
      CALL Send_Receive()
c                       construct working arrays for ART from primary 
c                      (belonging to  the node) and secondary (sent) data

      i_rec = 0
      DO k = 1, n_nodes
        i_rec = i_rec + recvcount(k)
      ENDDO
      If(iDebug_lg)WRITE(13,*)i_rec, ' primaries added to ',
     &     np_sb(node)       
      n_recv_prim = i_rec               ! received primaries 
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (k,i)
       DO k = 1, i_rec
        i = n_all + k
        x(i)   =  x_re(k)
        y(i)   =  y_re(k)
        z(i)   =  z_re(k)
        vx(i)  =  vx_re(k)
        vy(i)  =  vy_re(k)
        vz(i)  =  vz_re(k)
        pt(i)  =  pt_re(k)
        wpart(i) =  wpar_re(k)
        i_par(i) =  ip_re(k)
        If(iDebug_lg.and.(x(i)<1..or.x(i)>129.))
     &  WRITE(13,'(" Large:wrong X:",3g12.4,4i12)')x(i),y(i),z(i),
     &         i,k,node
      ENDDO  
      n_all = n_all+i_rec     ! redefine n_all

      If(iDebug_lg)WRITE(13,*)' Total Particles: ',
     +              n_all,' refine particles=',n_refin,
     +             ' prime=',np_sb(node)


      RETURN 
      End

c-------------------------------------------------------------------- 
      SUBROUTINE Send_Receive()
c-------------------------------------------------------------------- 
c
c  purpose: sends and reives data for particles
c                  
c     input: 
c           
c     output: 
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'


c               send integers with the lengths of all arrays which will be
c               sent (sendcount) and received (recvcount) to/from  all nodes  
               If(iDebug_lg)write (13,*)  ' Start sending counters '
      CALL mpi_allgather(sendcount, n_nodes, MPI_INTEGER,
     +                   sendcount_all, n_nodes, MPI_INTEGER,
     +                   MPI_COMM_WORLD, ierr)

      rdispl_new = 0
      DO i = 1, n_nodes
        rdispl(i)    = rdispl_new 
        recvcount(i) = sendcount_all(irank+1+n_nodes*(i-1))
        rdispl_new   = rdispl(i) + recvcount(i)
      ENDDO
         If(iDebug_lg)Then
            WRITE(13,'(''sendcount: '',100I8)') sendcount
            WRITE(13,'(''recvcount: '',100I8)') recvcount
            WRITE(13,'(''sdispl   : '',100I8)') sdispl
            WRITE(13,'(''rdispl   : '',100I8)') rdispl
         EndIf 
c               receive data from other nodes
               If(iDebug_lg)write (13,*)  ' Start sending data '
      t1 = MPI_Wtime()

      CALL mpi_alltoallv(x_se, sendcount, sdispl, MPI_REAL8,
     +                   x_re, recvcount, rdispl, MPI_REAL8,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(y_se, sendcount, sdispl, MPI_REAL8,
     +                   y_re, recvcount, rdispl, MPI_REAL8,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(z_se, sendcount, sdispl, MPI_REAL8,
     +                   z_re, recvcount, rdispl, MPI_REAL8,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(vx_se, sendcount, sdispl, MPI_REAL8,
     +                   vx_re, recvcount, rdispl, MPI_REAL8,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(vy_se, sendcount, sdispl, MPI_REAL8,
     +                   vy_re, recvcount, rdispl, MPI_REAL8,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(vz_se, sendcount, sdispl, MPI_REAL8,
     +                   vz_re, recvcount, rdispl, MPI_REAL8,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(pt_se, sendcount, sdispl, MPI_REAL,
     +                   pt_re, recvcount, rdispl, MPI_REAL,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(wpar_se, sendcount, sdispl, MPI_REAL,
     +                   wpar_re, recvcount, rdispl, MPI_REAL,
     +                  MPI_COMM_WORLD,ierr)
      CALL mpi_alltoallv(ip_se, sendcount, sdispl, MPI_INTEGER,
     +                   ip_re, recvcount, rdispl, MPI_INTEGER,
     +                  MPI_COMM_WORLD,ierr)
      t2 = MPI_Wtime()
      If(iDebug_lg)WRITE(13,*) ' transfer time total:',t2-t1 

      RETURN
      END
c-------------------------------------------------------------------- 
      SUBROUTINE Node_to_IJK(node,iN0)
c                             give {ijk} for node
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_mpi.h'
      DIMENSION iN0(3)
         k =INT((node-1)/(n_divx*n_divy))+1
         jj = node-(k-1)*n_divx*n_divy
         j =INT((jj-1)/n_divx) +1
         i = jj-(j-1)*n_divx
         iN0(1) =i
         iN0(2) =j
         iN0(3) =k
      Return
      End
c-------------------------------------------------------------------- 
c                           
      SUBROUTINE LoadBalance2()
c                   Uses: previos step timing time_load =w_node...
c                             boundary displacements are
c                             not more than iStepB
c                    Output: new l_div's
c                    incrementally moves boundaries
c                    towards the equal distribution of the load
c        recommendations:
c                    use iStepB =2. Code will chose values of
c                    displacement  iStepB or less.
c                    Code converges on partition, which is very
c                     close to instantenious equal-load partition
c                     Unfortunately 'equal-load' often has disbalance
c                     of factor 1.3 away from real equal load.
c                     This happens because  boundaries are
c                      moved in steps, not continiously. 
c                     If disbalance is 1.3 or larger, try LoadBalance1,
c                      which is based on random search of minimum
c                      load.
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      REAL time_load_balance_all(n_nodes)
      DIMENSION w_node(n_divx,n_divy,n_divz),
     &                iL_x( n_divx +1),iL_y(n_divx +1, n_divy +1),
     &                iL_z(n_divx +1,n_divy +1,n_divz +1)
      EQUIVALENCE (time_load_balance_all(1),w_node(1,1,1))
      REAL*8 wx(n_divx+1),wy(n_divy+1),wz(n_divz+1)
      REAL*8 ww(ng+1),wf(ng+1)
      DIMENSION  Ldx(n_divx+1),Ldy(n_divy+1),Ldz(n_divz+1)
      REAL*8 wwx,wwy,wwz,wc,wcc

      time_load_balance = end_time - start_time 
      mpi_send = 1



      CALL mpi_gather(time_load_balance,mpi_send,MPI_REAL,
     +              time_load_balance_all,mpi_send,MPI_REAL,
     +              root,MPI_COMM_WORLD,ierr)

      IF (irank .EQ. root) THEN
        if ( istep .eq. 1 ) then
          write(41,*) HEADER
c          write(41,662) 
        endif
       !write(13,'(4(4F7.2,2x))')w_node
 
        xx_min = 1.E9
        xx_max = -1.E9
        w_aver = 0.
        DO i = 1, n_nodes
          xx_min = min(xx_min,time_load_balance_all(i))
          xx_max = max(xx_max,time_load_balance_all(i))
          w_aver = w_aver + time_load_balance_all(i)
        ENDDO
          w_aver = w_aver /n_nodes     ! average time per node
        ratio = xx_max/max(xx_min,1.e-10)
        
         If(n_nodes.eq.8)Then
           write(41,62) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes),
     +          (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),j=2,n_divy),i=1,n_divx),
     &          (((l_divZ(i,j,k),k=2,n_divz),j=1,n_divy),i=1,n_divx)

        else  If(n_nodes.eq.27)Then
           write(41,65) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes),
     +          (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),i=1,n_divx),j=2,n_divy),
     &          (((l_divZ(i,j,k),k=2,n_divz),j=1,n_divy),i=1,n_divx)

        else  If(n_nodes.eq.32)Then
           write(41,66) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes),
     +          (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),j=2,n_divy),i=1,n_divx),
     &          (((l_divZ(i,j,k),k=2,n_divz),j=1,n_divy),i=1,n_divx)

        else  If(n_nodes.eq.64)Then
           write(41,67) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes)
           write(41,670) (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),j=2,n_divy),i=1,n_divx),
     &          (((l_divZ(i,j,k),k=1,n_divz),j=1,n_divy),i=1,n_divx)
        else
           write(41,68) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes)
           write(41,680) (l_divX(i),i=2,n_divx)
           write(41,682) ((l_divY(i,j),j=2,n_divy),i=1,n_divx)
c           write(41,684) 
c     &          (((l_divZ(i,j,k),k=1,n_divz),j=1,n_divy),i=1,n_divx)
     
        EndIf 
 62   format (i5,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,4(2f5.2,3x),/T30,'X:',i4,' Y:',2i4,' Z:',4i4)     ! nodes=8=2**3
 65   format (i5,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,3(3f5.2,3x),2(/T52,3(3f5.2,3x)),
     &     /T30,'X:',2i4,' Y:',6i4,' Z:',3(2x,6i4))     ! nodes=27=3**3
 66   format (i5,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,4(4f5.2,3x),/T52,4(4f5.2,3x),
     &     /T30,'X:',3i4,T50,' Y:',4(3i5,3x),     ! nodes=32=4x4x2
     &    4(/T52,4(4i5,3x))) 

 67   format (i5 ,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,4(4f5.2,3x),3(/T52,4(4f5.2,3x)))   ! nodes=64
 670  format (30x,' X:',3i4,T50,'Y:',3(4i5,3x),
     &    4(/T52,4(4i5,3x)))                             ! nodes=64
 680  format (30x,'X:',100i4)
 682  format (30x,'Y:',100(4i4))
 684  format (30x,'Z:',100i4)
 68   format (i5 ,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,5(5f5.2,1x),5(/T52,5(5f5.2,1x)))   ! generic

      !!!return  !!!!!!
 
         wx(n_divx+1) =0.
         wy(n_divy+1) =0.
         wz(n_divz+1) =0.
 
      wwx = 0.       ! x-direction -------------------------------
      Do i=1,n_divx     ! find projections on x- axis
         wx(i) = 0.          ! sum over all nodes 
         Do k=1,n_divz
         Do j=1,n_divy
            wx(i) = wx(i) + w_node(i,j,k)
         EndDo 
         EndDo 
         iLeft =l_divX(i)    ! find width of i-nodes
         iRight=l_divX(i+1)
         aver = wx(i)/MAX(iRight-iLeft,1) ! average weight
                                                   ! of a cell in i-th strip
         Do j=iLeft ,iRight-1           ! set weights for ng elements
            ww(j) =aver
         EndDo 
         wwx = wwx +wx(i)
      EndDo
c      write (*,*)  ' Sum of weights=',wwx
c      write (*,'(16f8.3)') (ww(i)/wwx*ng,i=1,ng+1)
      Do i=2,ng-1
         wf(i) = (ww(i-1)+ww(i)+ww(i+1))/3.
      EndDo 
      wf(1) = (ww(ng)+ww(1)+ww(2))/3.
      wf(ng) = (ww(ng)+ww(1)+ww(ng-1))/3.
      Do i=1,ng
         ww(i) = wf(i)
      EndDo 



      wideal = wwx/n_divx                  ! x- direction: set l_divx
      wc          = 0.
      wcc        = 0.
      node_c  = 1
      node_w = 0
      Do i=1,ng
         wc          = wc          +ww(i)
        wcc         = wcc         +ww(i)
         node_w = node_w +1
         If(wc.gt.wideal.or.wc.gt.wideal-ww(i+1)/2.)Then 
            node_c = node_c +1
            ldx(node_c) = i +1
c               write(*,'(" X:   weight=",F7.2," ideal=",F7.2,
c     &                      " left=",F7.2)')wc,wideal,wwx-wcc
            wc       = 0
            node_w = 0
         EndIf 
      EndDo
c                       Ldx is ideal position of boundaries
c                       Move l_divX toward Ldx, but limit
c                       displacements by iStepB         
c      write(*,*) 'OldivX=',l_divX
      Do i=2,n_divx
         idispl = MIN(iStepB,abs(l_divX(i)-ldx(i)))
         l_divX(i) =l_divX(i) + SIGN(idispl,ldx(i)-l_divX(i))
      EndDo 
      l_divX(1) = 1                          ! set boundary values for l_div'x
      l_divX(n_divx+1) = ng+1
c      write(*,*) 'l_divX=',l_divX
c      write(*,*) 'Ideal =',ldx

      Do i=1,n_divx     ! y-direction ----------------------
           wwy = 0.        ! find projections on y- axis, {i}-th strip      
         Do j=1,n_divy  
            wy(j) =0.        ! sum nodes with given {i}-coordinates
            Do k=1,n_divz
               wy(j) = wy(j) +w_node(i,j,k)
            EndDo 
            jLeft  =l_divY(i,j)
            jRight=l_divY(i,j+1)
            aver = wy(j)/MAX(jRight-jLeft,1)
            Do k=jLeft ,jRight-1
               ww(k) =aver
            EndDo 
            wwy = wwy +wy(j)
         EndDo     ! end j
c         write (*,*)  ' Sum of weights=',wwy,' strip=',i
c         write (*,'(16f8.3)') (ww(j)/wwy*ng,j=1,ng)
      Do ii=2,ng-1             ! smooth map of weights
         wf(ii) = (ww(ii-1)+ww(ii)+ww(ii+1))/3.
      EndDo 
      wf(1)   = (ww(ng)+ww(1)+ww(2))/3.
      wf(ng) = (ww(ng)+ww(1)+ww(ng-1))/3.
      Do ii=1,ng
         ww(ii) = wf(ii)
      EndDo 
c      write (*,*)  
c         write (*,'(16f8.3)') (ww(j)/wwy*ng,j=1,ng)
c      write (*,*)  

         wideal = wwy/n_divy                  ! set l_divy for i-th strip
         wc          = 0.
         wcc         = 0.
         node_c  = 1
         node_w = 0
         Do j=1,ng
            wc          = wc           +ww(j)
           wcc         = wcc         +ww(j)
            node_w = node_w +1
            If(wc.gt.wideal.or.wc.gt.wideal-ww(j+1)/2.)Then
               node_c = node_c +1
               ldy(node_c) = j +1
c               write(*,'(" Y:   weight=",F7.2," ideal=",F7.2,
c     &                      " left=",F7.2)')wc,wideal,wwy-wcc
               wc       = 0
               node_w = 0
            EndIf 
         EndDo
c                       Ldy is ideal position of boundaries
c                       Move l_divY toward Ldy, but limit
c                       displacements by iStepB         
c         write(*,*) 'OldivY=',(l_divY(i,j),j=1,n_divy+1)
         Do j=2,n_divy
            idispl = MIN(iStepB,abs(l_divY(i,j)-ldy(j)))
            l_divY(i,j) =l_divY(i,j) + SIGN(idispl,ldy(j)-l_divY(i,j))
         EndDo 
         l_divY(i,1) = 1                          ! set boundary values for l_div'y
         l_divY(i,n_divy+1) = ng+1
c         write(*,*) 'l_divY=',(l_divY(i,j),j=1,n_divy+1)
c         write(*,*) 'Ideal =',ldy
      EndDo              ! end i --  setup of l_divY   

      Do i=1,n_divx     ! z-direction ----------------------
      Do j=1,n_divy  
        wwz = 0.     ! find projections on z- axis, {i,j}-th strip      
                           ! sum nodes with given {i,j}-coordinates
         Do k=1,n_divz
            wz(k) = w_node(i,j,k)
            kLeft =l_divZ(i,j,k)
            kRight=l_divZ(i,j,k+1)
            aver = wz(k)/MAX(kRight-kLeft,1)
            Do kk=kLeft ,kRight-1
               ww(kk) =aver
            EndDo 
            wwz = wwz +wz(k)
         EndDo
c         write (13,*)  ' Sum of weights=',wwz,' strip=',i,j
c         write (13,'(16f8.3)') (ww(k)/wwz*ng,k=1,ng)
         Do ii=2,ng-1           ! smooth map of weights
            wf(ii) = (ww(ii-1)+ww(ii)+ww(ii+1))/3.
         EndDo 
          wf(1)  = (ww(ng)+ww(1)+ww(2))/3.
          wf(ng) = (ww(ng)+ww(1)+ww(ng-1))/3.
         Do ii=1,ng
             ww(ii) = wf(ii)
         EndDo 

         wideal = wwz/n_divz                  ! set l_divy for i-th strip
         wc          = 0.
         wcc         = 0.
         node_c  = 1
         node_w = 0
         Do k=1,ng
           wc          = wc           +ww(k)
           wcc         = wcc         +ww(k)
            node_w = node_w +1
            If(wc.gt.wideal.or.wc.gt.wideal-ww(k+1)/2.)Then
               node_c = node_c +1
               ldz(node_c) = k +1
c              write(13,'(" Z:   weight=",F7.2," ideal=",F7.2,
c     &          " left=",F7.2,2i4)')wc,wideal,wwz-wcc,node_c,k
               wc       = 0
               node_w = 0
            EndIf 
         EndDo
c         write (13,*)  ' i,j =',i,j
c         write(13,*) 'OldivZ=',(l_divZ(i,j,k),k=1,n_divz+1)
         Do k=2,n_divz
            idispl = MIN(iStepB,abs(l_divZ(i,j,k)-ldz(k)))
            l_divZ(i,j,k) =l_divZ(i,j,k) + 
     &                          SIGN(idispl,ldz(k)-l_divZ(i,j,k))

         EndDo 
         l_divZ(i,j,1) = 1                          ! set boundary values for l_div'z
         l_divZ(i,j,n_divz+1) = ng+1
c         write(13,*) 'l_divZ=',(l_divZ(i,j,k),k=1,n_divz+1)
c         write(13,*) 'Ideal =',ldz
      EndDo              ! end j -- direction   
      EndDo              ! end i -- direction  -- setup of l_divZ 



c                                     --- Adjust Y-boundaries to avoid too small nodes
            Do i   = 2,n_divx     ! loop over all divisions in this direction
               If(l_divX(i).lt.l_divX(i-1)+min_cell_numb)
     &             l_divX(i) =l_divX(i-1)+min_cell_numb 
            EndDo                   
            Do i   = 1,n_divx     ! loop over all divisions in this direction
            Do j   = 2,n_divy     
               If(l_divY(i,j).lt.l_divY(i,j-1)+min_cell_numb
     &             .and.l_divY(i,j).lt.ng+1-min_cell_numb)
     &             l_divY(i,j) =l_divY(i,j-1)+min_cell_numb 
            EndDo
            EndDo 
            Do i   = 1,n_divx     ! loop over all divisions in this direction
            Do j   = 1,n_divy     ! loop over all divisions in this direction
            Do k   = 2,n_divz     ! loop over all divisions in this direction
               If(l_divZ(i,j,k).lt.l_divZ(i,j,k-1)+min_cell_numb)
     &             l_divZ(i,j,k) =l_divZ(i,j,k-1)+min_cell_numb 
            EndDo
            EndDo
            EndDo

      EndIf             ! end root out put
      !! return !!!!
      is =  n_divX +1 
      CALL MPI_BCAST(l_divX,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      is =  n_divX*(n_divY +1) 
      CALL MPI_BCAST(l_divY,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      is =  n_divX*n_divY* (n_divZ +1)
      CALL MPI_BCAST(l_divZ,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

      WRITE(13,'(''new l_divXY:'',T15,2x,2I4,"  Y:",5(2x,3i4))') 
     &    (l_divx(i),i=2,n_divx),
     &    ((l_divy(i,j),i=1,n_divx),j=2,n_divy)          
!      WRITE(13,'(''new l_div Z:'',T15,4(2x,3i4),4(/T15,4(2x,3I4)))')
!     &           (((l_divz(i,j,k),
!     &             i=1,n_divx),j=1,n_divy),k=2,n_divz)

c      close (13)
c      CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
c      STOP

      RETURN 
      END
c-------------------------------------------------------------------- 
c                        Set boundaries using many random displacements
c                         Uses:     time_load_balance_all
c                         Output: new l_div's
      SUBROUTINE LoadBalance1()
c                Recommendations:
c                        Code has slow convergence. iStepB is
c                        limited to 1. For larger displacements
c                        code does not reach satisfactory balance
c                        in test runs. 
c                        Use this code in case when LoadBalance2
c                        gives load disbalance (max/average) larger
c                        than 1.3. During transition disbalance may
c                        temporary increase before it decreases. It takes
c                        about 10 iterations before improvement is
c                        observed. Balance may oscillate.
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      PARAMETER (n_divt =n_divx*n_divy*n_divz-1) 
                                                         !  total number of variable boundaries 
c      PARAMETER (n_divL =MIN(3**n_divt,60000) ) ! number of combinations
      PARAMETER (n_divL =10000 ) ! number of combinations

      REAL time_load_balance_all(n_nodes)
      DIMENSION iComb(n_divt,n_divL)
      DIMENSION w_node(n_divx,n_divy,n_divz),
     &                iL_x( n_divx +1),iL_y(n_divx +1, n_divy +1),
     &                iL_z(n_divx +1,n_divy +1,n_divz +1)
      SAVE iComb
      EQUIVALENCE (time_load_balance_all(1),w_node(1,1,1))

      time_load_balance = end_time - start_time 
      mpi_send = 1

c      write (13,*)  ' node=',irank+1,' time=',time_load_balance
c      write(13,*)' size: n_divt=',n_divt,' n_divL=',n_divL
c      write (13,*) ' In Balance 1. nodes=',n_nodes


      CALL mpi_gather(time_load_balance,mpi_send,MPI_REAL,
     +              time_load_balance_all,mpi_send,MPI_REAL,
     +              root,MPI_COMM_WORLD,ierr)

      IF (irank .EQ. root) THEN
        if ( istep .eq. 1 ) then
          write(41,*) HEADER
c          write(41,662) 
        endif

c      write(13,*)' size: n_divt=',n_divt,' n_divL=',n_divL
!      write(13,'(8(5F7.2,2x))')w_node
!      write (13,*) ' In Balance 1'


c      call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
c      CALL mpi_finalize(ierr)
c      stop

        xx_min = 1.E9
        xx_max = -1.E9
        w_aver = 0.
        DO i = 1, n_nodes
          xx_min = min(xx_min,time_load_balance_all(i))
          xx_max = max(xx_max,time_load_balance_all(i))
          w_aver = w_aver + time_load_balance_all(i)
        ENDDO
          w_aver = w_aver /n_nodes     ! average time per node
        ratio = xx_max/max(xx_min,1.e-10)
        If(n_nodes.eq.8)Then
           write(41,62) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes),
     +          (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),j=2,n_divy),i=1,n_divx),
     &          (((l_divZ(i,j,k),k=2,n_divz),j=1,n_divy),i=1,n_divx)

        else  If(n_nodes.eq.27)Then
           write(41,65) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes),
     +          (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),i=1,n_divx),j=2,n_divy),
     &          (((l_divZ(i,j,k),k=2,n_divz),j=1,n_divy),i=1,n_divx)

        else  If(n_nodes.eq.32)Then
c           write (*,*)  ' inside format. node=',irank+1
           write(41,66) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes),
     +          (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),j=2,n_divy),i=1,n_divx),
     &          (((l_divZ(i,j,k),k=2,n_divz),j=1,n_divy),i=1,n_divx)

        else  If(n_nodes.eq.64)Then
           write(41,67) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes)
           write(41,670) (l_divX(i),i=2,n_divx),
     +          ((l_divY(i,j),j=2,n_divy),i=1,n_divx),
     &          (((l_divZ(i,j,k),k=1,n_divz),j=1,n_divy),i=1,n_divx)
        else
           write(41,68) istep , xx_max, xx_min, w_aver,
     +          (time_load_balance_all(i)/xx_max, i = 1,n_nodes)
           write(41,680) (l_divX(i),i=2,n_divx)
           write(41,682) ((l_divY(i,j),j=2,n_divy),i=1,n_divx)
c           write(41,684) 
c     &          (((l_divZ(i,j,k),k=1,n_divz),j=1,n_divy),i=1,n_divx)
     
        EndIf 
 62   format (i5,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,4(2f5.2,3x),/T30,'X:',i4,' Y:',2i4,' Z:',4i4)     ! nodes=8=2**3
 65   format (i5,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,3(3f5.2,3x),2(/T52,3(3f5.2,3x)),
     &     /T30,'X:',2i4,' Y:',6i4,' Z:',3(2x,6i4))     ! nodes=27=3**3
 66   format (i5,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,4(4f5.2,3x),/T52,4(4f5.2,3x),
     &     /T30,'X:',3i4,T50,' Y:',4(3i5,3x),     ! nodes=32=4x4x2
     &    4(/T52,4(4i5,3x))) 

 67   format (i5 ,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,4(4f5.2,3x),3(/T52,4(4f5.2,3x)))   ! nodes=64
 670  format (30x,' X:',3i4,T50,'Y:',3(4i5,3x),
     &    4(/T52,4(4i5,3x)))                             ! nodes=64
 680  format (30x,'X:',100i4)
 682  format (30x,'Y:',100(4i4))
 684  format (30x,'Z:',100i4)
 68   format (i5 ,' load: max=',F8.2,' min=',F8.2,' aver=',F8.2,
     &    T52,5(5f5.2,3x),5(/T52,5(5f5.2,3x)))   ! generic


      iStepP = 1             ! step to move node boundaries
      w_aver =0.
      Do i=1,n_divx+1 ! left and right boundaries do not move
         iL_x(i) =0          ! set all of them to zero; change what is needed later
         Do j=1,n_divy+1
            iL_y(i,j) =0
            Do k=1,n_divz+1
               iL_z(i,j,k) =0
            EndDo 
         EndDo 
      EndDo 
c   
c                    Set all possible combinations of boundaries
c                    This can be done only once for the run                          
      CALL SetComb(iComb)
c      call mpi_barrier(MPI_COMM_WORLD,ierr)
c      write (13,*) ' In Balance 1: setComb done'
c      close (13)
c      CALL mpi_finalize(ierr)
c      stop

c   
c      open(80,file='/data1/aklypin/Box20/outB.dat')
      Niter_bal =n_divL
      If(n_divt.lt.8) Niter_bal =3**n_divt
               bal_min = 1.e20  ! minimum of current balance
               qbal_min   = 0.       ! estimate of quality of balance
      Do iter =1, Niter_bal   ! make iterations to find min balance
c        write (80,*) ' In Balance 1: iter=',iter,bal_min
                                         ! copy realization to iL_* arrays
         ibnd = 0                    ! current boundary number
            Do i   = 2,n_divx     ! loop over all divisions in this direction
               ibnd = ibnd +1
                  iL_x(i) = iComb(ibnd,iter)
             EndDo                   
            Do i   = 1,n_divx
            Do j   = 2,n_divy     ! loop over all divisions in this direction
               ibnd = ibnd +1
                  iL_y(i,j) = iComb(ibnd,iter)
             EndDo
             EndDo
            Do i   = 1,n_divx     ! loop over all divisions in this direction
            Do j   = 1,n_divy     ! loop over all divisions in this direction
            Do k   = 2,n_divz     ! loop over all divisions in this direction
               ibnd = ibnd +1
                  iL_z(i,j,k) = iComb(ibnd,iter)
             EndDo
             EndDo
             EndDo
         CALL TryBalance(w_node,iL_x,iL_y, iL_z,balance,qbal)
         If(balance.le.1.01*bal_min)Then    ! store minimum of max time
            If(balance.le.0.99*bal_min.or.qbal.lt.qbal_min)then
              kBal = iter
              bal_min =balance
              qbal_min = qbal
            EndIf 
          EndIf 
      EndDo                        ! end iter
               qbal_min = sqrt(qbal_min/n_nodes -w_aver**2)
c       write (13,*) ' new balance=',bal_min,kBal
c       write (14,*) ' new balance =',bal_min,' iter=',kBal,qbal_min
c                                     --- Set new boundaries to min.balance
         ibnd = 0                    ! current boundary number
            Do i   = 2,n_divx     ! loop over all divisions in this direction
               ibnd = ibnd +1
                  iL_x(i)     = iComb(ibnd,kBal)
                  l_divX(i) = l_divX(i) + iL_x(i)*iStepP 
             EndDo                   
            Do i   = 1,n_divx     ! loop over all divisions in this direction
            Do j   = 2,n_divy     !
               ibnd = ibnd +1
                  iL_y(i,j)     = iComb(ibnd,kBal)
                  l_divY(i,j) = l_divY(i,j) + iL_y(i,j)*iStepP 
             EndDo
             EndDo 
            Do i   = 1,n_divx     ! loop over all divisions in this direction
            Do j   = 1,n_divy     ! loop over all divisions in this direction
            Do k   = 2,n_divz     ! loop over all divisions in this direction
               ibnd = ibnd +1
                  iL_z(i,j,k) = iComb(ibnd,kBal)
                  l_divZ(i,j,k) = l_divZ(i,j,k) + iL_z(i,j,k)*iStepP
             EndDo
             EndDo
             EndDo

c                                     --- Adjust boundaries to avoid too small nodes
            Do i   = 2,n_divx     ! loop over all divisions in this direction
               If(l_divX(i).lt.l_divX(i-1)+min_cell_numb)
     &             l_divX(i) =l_divX(i-1)+min_cell_numb 
            EndDo                   
            Do i   = 1,n_divx     ! loop over all divisions in this direction
            Do j   = 2,n_divy     
               If(l_divY(i,j).lt.l_divY(i,j-1)+min_cell_numb
     &             .and.l_divY(i,j).lt.ng+1-min_cell_numb)
     &             l_divY(i,j) =l_divY(i,j-1)+min_cell_numb 
            EndDo
            EndDo 
            Do i   = 1,n_divx     ! loop over all divisions in this direction
            Do j   = 1,n_divy     ! loop over all divisions in this direction
            Do k   = 2,n_divz     ! loop over all divisions in this direction
               If(l_divZ(i,j,k).lt.l_divZ(i,j,k-1)+min_cell_numb)
     &             l_divZ(i,j,k) =l_divZ(i,j,k-1)+min_cell_numb 
            EndDo
            EndDo
            EndDo
      EndIf             ! end root out put
      is =  n_divX +1 
      CALL MPI_BCAST(l_divX,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      is =  n_divX*(n_divY +1) 
      CALL MPI_BCAST(l_divY,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      is =  n_divX*n_divY* (n_divZ +1)
      CALL MPI_BCAST(l_divZ,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

c      WRITE(13,'(''new l_divXY:'',T15,2x,3I4,"  Y:",32(2x,3i4))') 
c     &    (l_divx(i),i=2,n_divx),
c     &    ((l_divy(i,j),j=2,n_divy),i=1,n_divx)          
c      WRITE(13,'(''new l_div Z:'',T15,4(2x,4i4),4(/T15,4(2x,4I4)))')
c     &           (((l_divz(i,j,k),
c     &             i=1,n_divx),j=1,n_divy),k=2,n_divz)

c      close (13)
c      CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
c      STOP

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
c-------------------------------------------------------------------- 
c                                                            Estimate load balance 
      SUBROUTINE TryBalance(w_node,iL_x,iL_y, iL_z,balance,qbal)
c               Estimate load balance using old boundaries l_div
c                and new displacements of boundaries iL_div
c               w_node = load for each node
c               Use linear approximation.
c               balance = estimate of maximum load of all nodes
c               qbal = sum of squares of estimated times
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      Common /NewBalance/  w_NewBal(n_divx,n_divy,n_divz)
       DIMENSION w_node(n_divx,n_divy,n_divz),
     &                iL_x( n_divx +1),iL_y(n_divx +1, n_divy +1),
     &                iL_z(n_divx +1,n_divy +1,n_divz +1),iTake(-1:1) 
 
      balance = -1.e12 
      qbal       = 0.
      iStepP =1     ! step to move node boundaries 

      Do i=1,n_divX                    ! loop over all nodes
            iTake(-1) =   -iStepP*min(0,iL_x(i))
            iTake(0) = l_divX(i+1)- l_divX(i)+ iStepP*
     &                      (min(0,iL_x(i+1))-max(0,iL_x(i)))
             iTake(1) =    iStepP*max(0,iL_x(i+1))
           If(iTake(-1).lt.0.)write (13,*)' Error xTake=',iTake(-1),
     &        i,iL_x(i)
           If(iTake(0).lt.0.)write (13,*)' Error xTake=',iTake(0)
           If(iTake(1).lt.0.)write (13,*)' Error xTake=',iTake(1)
           iMinn = 0            ! find which nearest nodes to test
           if(iL_x(i).eq.-1)iMinn = -1
           iMaxn = 0
           if(iL_x(i+1).eq.1)iMaxn = 1
 
      Do j=1,n_divY                       

           Do k=1,n_divZ             ! {i,j,k} - is the current node  
           w_new = 0.                   ! counter for estaimate of new load
           Do ib=iMinn,iMaxn       ! go over all neighbors of current node
                ii =i+ib
c                if(ii.lt.1)stop ' error balance: ii<1'
c                if(ii.gt.n_divx)stop 'error balance: ii>n_divx'
           Do jj=1,n_divY
                jLeft = max(l_divY(ii,jj),l_divY(i,j)
     &                                         +iStepP*iL_y(i,j))
                jRight= min(l_divY(ii,jj+1),l_divY(i,j+1)
     &                                         +iStepP*iL_y(i,j+1))
                jTake = max(jRight -jLeft,0)

           Do kk=1,n_divZ                 
                kLeft = max(l_divZ(ii,jj,kk),l_divZ(i,j,k)
     &                                         +iStepP*iL_z(i,j,k))
                kRight= min(l_divZ(ii,jj,kk+1),l_divZ(i,j,k+1)
     &                                         +iStepP*iL_z(i,j,k+1))
                kTake = max(kRight -kLeft,0)

                V_new = iTake(ib)*jTake*kTake
                V_old =(l_divX(ii+1)-l_divX(ii)) *
     &                    (l_divY(ii,jj+1)-l_divY(ii,jj)) *
     &                    (l_divZ(ii,jj,kk+1)-l_divZ(ii,jj,kk)) 
                if(V_new.gt.V_old)Then
                 write(13,*)' Error TryBalance: V_new =',V_new,
     &          '  V_old=', V_old
                 write (13,*)  '    test    node=',ib,jb,kk
                 write (13,*)  '    current node=',i,j,k
                EndIf 
              
                w_new = w_new +w_node(ii,jj,kk)*
     &                         V_new /V_old
c                write (*,50) i,j,k,ii,jj,kk,w_node(ii,jj,kk),
c     &             V_new /V_old
c 50             format(' node=',3i3,' takes from=',3i3,
c     &             ' with time=',f8.2,' fraction=',f7.3)
             EndDo                  ! kb - loop in z-direction
             EndDo                  ! jb
             EndDo                  ! kb
c                 write (*,60) w_new,w_node(i,j,k),i,j,k
c 60             format(15x,'new time=',f8.3,'old time=',f8.3,3i3)
c             w_NewBal(i,j,k) = w_new
             balance = max(balance,w_new)
             qbal       = qbal + w_new**4
c             qq          = qq   + w_new
c             balance = min(balance,w_new)
      EndDo 
      EndDo 
      EndDo 

      dummy   = balance
      balance = sqrt(sqrt(qbal/(n_divX*n_divY*n_divZ)))
      qbal    = dummy
     
c      write(*,*) qq,qbal
      RETURN 
      End
c-------------------------------------------------------------------- 
      SUBROUTINE SetComb(iComb)
c-------------------------------------------------------------------- 
c     purpose :   create all combinations of boundaries for n_divt
c
      include 'a_tree.h'
      include 'a_mpi.h'

      PARAMETER (mcomb3 =(n_divz-1)*n_divy*n_divx) ! 
      PARAMETER (mcomb2 =(n_divy-1)*n_divx)             ! 
      PARAMETER (n_comx = 3**(n_divx-1) ) 
      PARAMETER (n_divt =n_divx*n_divy*n_divz-1) 
                                                         !  total number of variable boundaries 
c      PARAMETER (n_divL =MIN(3**n_divt,60000) ) ! number of combinations
      PARAMETER (n_divL =10000 ) ! number of combinations

      PARAMETER (nxy      = n_divx*n_divy) !  
      DIMENSION iCom1d(n_divx-1,n_comx)
      DIMENSION iCom2d(n_divy-1,n_divx)
      DIMENSION iCom3d(n_divz-1,n_divy,n_divx)
      DIMENSION kCom2d(mcomb2)
      DIMENSION kCom3d(mcomb3)
      DIMENSION iComb(n_divt,n_divL)
      SAVE iCom1d,iCom2d
      External iCos, iCos2, iSin, iConst, nConst, iZero,
     &               nCos, nCos2, nSin

      equivalence (kCom2d(1),iCom2d(1,1))
      equivalence (kCom3d(1),iCom3d(1,1,1))
 

      nseed = 121071
      If(n_divt.lt.8)Then   !  generate all possible combinations
         iloop = 1                        ! number contigious regions
         Nloop = 3**(n_divt-1)   ! length of a contigious region
         Do i= 1,n_divt
            do k=1,iloop
            ioffset = 3*Nloop*(k-1)  ! offset for current region
            Do j=1,Nloop
               iComb(i,j+ioffset            ) = -1
               iComb(i,j+ioffset +  Nloop) = 0
               iComb(i,j+ioffset +2*Nloop) = 1
c               n = n +3
            EndDo 
            EndDo 
            iloop =iloop*3
            Nloop =Nloop/3
         EndDo 
         Return
      EndIf 

c      write (*,*) ' n_divx=',n_divx,' n_divy=',n_divy
c      write (*,*)  ' n_comx =',n_comx,n_divt,mcomb2,mcomb3
      Do j=1,n_divL          ! initialize the matrix of combinations
         Do i=1,n_divt
            iComb(i,j) = 0
         EndDo 
      EndDo 
      CALL SetXcom(iCom1d)     ! all combinations in x-direction
          do i=1,n_divt
            iComb(i,1) =-1
            iComb(i,2) =0
            iComb(i,3) =1
         EndDo 
      iter = 4
      Do j=1,n_comx                    ! copy x-combinations to main matrix 
         Do jj   = 0,59
         Do i=1,n_divx-1
            iComb(i,iter+jj)=iCom1d(i,j)
         EndDo 
         EndDo 
         CALL SetZcom(iComb(nxy,iter  ),iZero,     iZero,    iZero)
         CALL SetZcom(iComb(nxy,iter+1),nConst,nConst,nConst)
         CALL SetZcom(iComb(nxy,iter+2),iConst,iConst, iConst)
         CALL SetZcom(iComb(nxy,iter+3),iCos,iConst, iConst)

         CALL SetZcom(iComb(nxy,iter+4),iCos,    iCos  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+5),iCos,    iSin  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+6),iCos,    nSin  ,iSin  )
         CALL SetZcom(iComb(nxy,iter+7),iCos,    nCos  ,nSin  )
         CALL SetZcom(iComb(nxy,iter+8),iCos,    nCos  ,nCos  )
         CALL SetZcom(iComb(nxy,iter+9),iCos,    iCos2 ,nCos2 )
         CALL SetZcom(iComb(nxy,iter+10),iCos,    iConst ,nCos  )
         CALL SetZcom(iComb(nxy,iter+11),iCos,    nCos2 ,iConst )

         CALL SetZcom(iComb(nxy,iter+12),iSin,    iCos  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+13),iSin,    iSin  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+14),iSin,    nSin  ,iSin  )
         CALL SetZcom(iComb(nxy,iter+15),iSin,    nCos  ,nSin  )
         CALL SetZcom(iComb(nxy,iter+16),iSin,    nCos  ,nCos  )
         CALL SetZcom(iComb(nxy,iter+17),iSin,    iCos2 ,nCos2 )
         CALL SetZcom(iComb(nxy,iter+18),iSin,    iConst ,nCos  )
         CALL SetZcom(iComb(nxy,iter+19),iSin,    nCos2 ,iConst )

         CALL SetZcom(iComb(nxy,iter+20),nCos,    iCos  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+21),nCos,    iSin  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+22),nCos,    nSin  ,iSin  )
         CALL SetZcom(iComb(nxy,iter+23),nCos,    nCos  ,nSin  )
         CALL SetZcom(iComb(nxy,iter+24),nCos,    nCos  ,nCos  )
         CALL SetZcom(iComb(nxy,iter+25),nCos,    iCos2 ,nCos2 )
         CALL SetZcom(iComb(nxy,iter+26),nCos,    iConst ,nCos  )
         CALL SetZcom(iComb(nxy,iter+27),nCos,    nCos2 ,iConst )

         CALL SetZcom(iComb(nxy,iter+28),nSin,    iCos  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+29),nSin,    iSin  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+30),nSin,    nSin  ,iSin  )
         CALL SetZcom(iComb(nxy,iter+31),nSin,    nCos  ,nSin  )
         CALL SetZcom(iComb(nxy,iter+32),nSin,    nCos  ,nCos  )
         CALL SetZcom(iComb(nxy,iter+33),nSin,    iCos2 ,nCos2 )
         CALL SetZcom(iComb(nxy,iter+34),nSin,    iConst ,nCos  )
         CALL SetZcom(iComb(nxy,iter+35),nSin,    nCos2 ,iConst )


         CALL SetZcom(iComb(nxy,iter+36),iCos2,    iCos  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+37),iCos2,    iSin  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+38),iCos2,    nSin  ,iSin  )
         CALL SetZcom(iComb(nxy,iter+39),iCos2,    nCos  ,nSin  )
         CALL SetZcom(iComb(nxy,iter+40),iCos2,    nCos  ,nCos  )
         CALL SetZcom(iComb(nxy,iter+41),iCos2,    iCos2 ,nCos2 )
         CALL SetZcom(iComb(nxy,iter+42),iCos2,    iConst ,nCos  )
         CALL SetZcom(iComb(nxy,iter+43),iCos2,    nCos2 ,iConst )

         CALL SetZcom(iComb(nxy,iter+44),nCos2,    iCos  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+45),nCos2,    iSin  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+46),nCos2,    nSin  ,iSin  )
         CALL SetZcom(iComb(nxy,iter+47),nCos2,    nCos  ,nSin  )
         CALL SetZcom(iComb(nxy,iter+48),nCos2,    nCos  ,nCos  )
         CALL SetZcom(iComb(nxy,iter+49),nCos2,    iCos2 ,nCos2 )
         CALL SetZcom(iComb(nxy,iter+50),nCos2,    iConst ,nCos  )
         CALL SetZcom(iComb(nxy,iter+51),nCos2,    nCos2 ,iConst )

         CALL SetZcom(iComb(nxy,iter+52),iConst,    iCos  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+53),iConst,    iSin  ,iCos  )
         CALL SetZcom(iComb(nxy,iter+54),iConst,    nSin  ,iSin  )
         CALL SetZcom(iComb(nxy,iter+55),iConst,    nCos  ,nSin  )
         CALL SetZcom(iComb(nxy,iter+56),iConst,    nCos  ,nCos  )
         CALL SetZcom(iComb(nxy,iter+57),iConst,    iCos2 ,nCos2 )
         CALL SetZcom(iComb(nxy,iter+58),iConst,    iConst ,nCos  )
         CALL SetZcom(iComb(nxy,iter+59),iConst,    nCos2 ,iConst )

         Do jj=0,5
            i = iter+10*jj 
            CALL SetYcom(iComb(n_divx,i  ),iZero,iConst  )
            CALL SetYcom(iComb(n_divx,i+1),iConst,iSin   )
            CALL SetYcom(iComb(n_divx,i+2),nConst,iSin  )
            CALL SetYcom(iComb(n_divx,i+3),iSin   ,iCos )
            CALL SetYcom(iComb(n_divx,i+4),iSin   ,iSin )
            CALL SetYcom(iComb(n_divx,i+5),iSin   ,nCos)
            CALL SetYcom(iComb(n_divx,i+6),iCos   ,iConst  )
            CALL SetYcom(iComb(n_divx,i+7),iCos   ,iCos )
            CALL SetYcom(iComb(n_divx,i+8),iCos   ,nSin)
            CALL SetYcom(iComb(n_divx,i+9),iConst ,iCos)
         EndDo 
         iter = iter +60               ! advance the offset
      EndDo 
      iter = iter -1                     ! make adjusment to get current iter
c      write (*,*)  ' current iteration counter =',iter
c     Do jj =1,iter
c        write(*,'(i7,2x,64i2)') jj,(iComb(i,jj),i=1,n_divt)
c     EndDo 
c                         Fill up the rest of the table:
c                        copy as many replicas of x-table as possible
      j = 0           
 300  If(iter.lt.n_divL)Then
            iter = iter +1
            j = j +1
            If(j.gt.n_comx)j =1
            Do i=1,n_divx-1                  ! copy the x-table
               iComb(i,iter)=iCom1d(i,j)
            EndDo 
            Do i =n_divx,n_divt            ! fill the rest with random numbers
               iComb(i,iter) =INT(3.*RANDd(Nseed))-1
            EndDo
             goto 300
         EndIf 
      write (*,*) ' Set Combinations for N   =',n_divt,' boundaries'
      write (*,*) ' Number of Combinations   =',n_divL,' =',iter

c      write (*,*)  ' current iteration counter =',iter
c      Do jj =1,iter
c         write(*,'(i7,2x,64i2)') jj,(iComb(i,jj),i=1,n_divt)
c      EndDo 
c
      RETURN 
      End
c-------------------------------------------------------------------- 
c      make a vector of displacements for z-component
c      the displacement  is = iF1(x)*iF2(y)*iF3(z)
c-------------------------------------------------------------------- 
      SUBROUTINE SetZcom(iCom3d,iF1,iF2,iF3)
      include 'a_tree.h'
      include 'a_mpi.h'
      DIMENSION iCom3d(n_divz-1,n_divy,n_divx)
      external iF1,iF2,iF3
       Do i=1,n_divx
         iC  = iF1(i)
      Do j=1,n_divy
         jC = iC*iF2(j)
         Do k=1,n_divz-1
                     iCom3d(k,j,i) = jC*iF3(k)
         EndDo 
      EndDo 
      EndDo 
      End
c-------------------------------------------------------------------- 
c      make a vector of displacements for y-component
c      the displacement  is = iF1(x)*iF2(y)
c-------------------------------------------------------------------- 
      SUBROUTINE SetYcom(iCom2d,iF1,iF2)
      include 'a_tree.h'
      include 'a_mpi.h'
      DIMENSION iCom2d(n_divy-1,n_divx)
      external iF1,iF2
       Do i=1,n_divx
         iC  = iF1(i)
         Do j=1,n_divy-1
                     iCom2d(j,i) = iC*iF2(j)
         EndDo 
      EndDo 
      End

c-------------------------------------------------------------------- 
c      make a matrix of all displacements for x-component
c
c-------------------------------------------------------------------- 
      SUBROUTINE SetXcom(iCom1d)
      include 'a_tree.h'
      include 'a_mpi.h'
      PARAMETER (n_comx = 3**(n_divx-1) ) ! number of all x-combinations 
      DIMENSION iCom1d(n_divx-1,n_comx)

       iloop = 1                        ! number contigious regions
      Nloop = 3**(n_divx-2)   ! length of a contigious region
      Do i= 1,n_divx-1
         do k=1,iloop
         ioffset = 3*Nloop*(k-1)  ! offset for current region
         Do j=1,Nloop
            iCom1d(i,j+ioffset            ) = -1
            iCom1d(i,j+ioffset +  Nloop) = 0
            iCom1d(i,j+ioffset +2*Nloop) = 1
         EndDo 
         EndDo 
         iloop =iloop*3
         Nloop =Nloop/3
      EndDo 

      return
      End

c-------------------------------------------------------------------- 
      Integer Function iSin(i)
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_mpi.h'
          iSin = Min(Max(INT(2.*sin((i-1.)/(n_divx-1)*pi)),-1),1)
       RETURN 
       End
c-------------------------------------------------------------------- 
      Integer Function nSin(i)
c-------------------------------------------------------------------- 
         nSin = -iSin(i) 
       RETURN 
       End
c-------------------------------------------------------------------- 
      Integer Function nCos(i)
c-------------------------------------------------------------------- 
         nCos = -iCos(i) 
       RETURN 
       End
c-------------------------------------------------------------------- 
      Integer Function nCos2(i)
c-------------------------------------------------------------------- 
         nCos2 = -iCos2(i) 
       RETURN 
       End
c-------------------------------------------------------------------- 
      Integer Function iCos(i)
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_mpi.h'
         iCos = Min(Max(INT(2.*cos((i-1.)/(n_divx-1)*pi)),-1),1)
      RETURN 
      End
c-------------------------------------------------------------------- 
      Integer Function iConst(i)
c-------------------------------------------------------------------- 
      iConst = 1
      RETURN 
      End
c-------------------------------------------------------------------- 
      Integer Function nConst(i)
c-------------------------------------------------------------------- 
      nConst = -1
      RETURN 
      End
c-------------------------------------------------------------------- 
      Integer Function iZero(i)
c-------------------------------------------------------------------- 
      iZero = 0
      RETURN 
      End
c-------------------------------------------------------------------- 
      Integer Function iCos2(i)
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_mpi.h'
         iCos2 = Min(Max(INT(2.*cos(2.*(i-1.)/(n_divx-1)*pi)),-1),1)
      RETURN 
      End


