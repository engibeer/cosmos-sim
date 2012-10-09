c.................. Routines  FOR MPI VERSION .....................

c     -----------------------------------
      subroutine Redistribute_Primaries() 
c     -----------------------------------
c
c  purpose: redistribute primaries after each step 
c  i.e. send primaries which left the node to the corresponding node
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'

!      COMMON /redist/ indx(np),no_se(np)

      node = irank + 1

c      WRITE(13,*) ' Redistribute Primaries: '

      DO i = 1, max_send
        ip_se(i) = 0
      ENDDO
      DO i = 1, np
        no_se(i) = 0
      ENDDO

      i_send = 0
      i_new = 0
      i_store =0
      DO i = 1, np_sb(node)
            DO  n = 1,n_divX          ! find node to which the particles belongs
              IF(x(i) .GE. l_divX(n)  .AND.  
     +           x(i)  .LT. l_divX(n+1) ) THEN
                id1 = n          ! detected in n-th node in x-direction
                GOTO 222
              ENDIF 
            ENDDO
            write (13,*) ' Redistribute: Error find node'  
 222        CONTINUE
            DO  n = 1,n_divY 
              IF(y(i) .GE. l_divY(id1,n)  .AND.  
     +           y(i)  .LT. l_divY(id1,n+1) ) THEN
                id2 = n          ! detected in n
                GOTO 223
              ENDIF 
            ENDDO
            write (13,*) ' Redistribute: Error find node'  
 223        CONTINUE

            DO  n = 1,n_divZ 
              IF(z(i) .GE. l_divZ(id1,id2,n)  .AND.  
     +           z(i)  .LT. l_divZ(id1,id2,n+1) ) THEN
                id3 = n          ! detected in n
                GOTO 224
              ENDIF 
            ENDDO
            write (13,*) ' Redistribute: Error finding node=',id,
     &            x(i),y(i),z(i)  
 224        CONTINUE


         node_send  = id1+(id2-1)*n_divX+(id3-1)* n_divX*n_divY
         IF(node_send .NE. node)Then
            i_send = i_send + 1
            i_store = np_sb(node) + i_send  ! store temporary 
            no_se(i_store) = node_send ! i_store_th particles must be send
                                       ! to node # node_send
            x(i_store)   =  x(i)       ! store temporary at the end
            y(i_store)   =  y(i)       ! of the primaries
            z(i_store)   =  z(i)    
            vx(i_store)  =  vx(i)   
            vy(i_store)  =  vy(i)    
            vz(i_store)  =  vz(i)   
            pt(i_store)  =  pt(i)  
            wpart(i_store)=  wpart(i)   
            i_par(i_store)  =  i_par(i)              
         Else 
                                       ! particle belongs to the given node
            i_new = i_new + 1     
            x(i_new)   =  x(i)      ! filling up all arrays with "new"   
            y(i_new)   =  y(i)      ! primaries
            z(i_new)   =  z(i)    
            vx(i_new)  =  vx(i)   
            vy(i_new)  =  vy(i)    
            vz(i_new)  =  vz(i)   
            pt(i_new)  =  pt(i)  
            wpart(i_new)=  wpart(i)   
            i_par(i_new)  =  i_par(i)              
         EndIf 
      ENDDO
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        write(13,*) ' node=',node,' Npart=',np_sb(node)
        write(13,*) ' i_send=',i_send,' np=',np
        write(13,*) ' i_store=',i_store,' i_new=',i_new
c        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
c        CALL mpi_finalize(ierr)
c        STOP

      IF(np_sb(node)+i_send .GT. np) THEN
         write(name,'(a,i4.4,a)')'ERROR_RUN',node,'.DAT'
         accesstype = 'append'
         j123=456
         call Open_ASCII_FILE(19,name,accesstype)
         rewind 19
         write(19,*) j123, 
     +             ' maximum number of particles exeeded in ', 
     +             'Redistribute_Primaries'
          write(19,*) 'np_sb(node)+i_send .GT. np'
          write(19,*)  np_sb(node)+i_send,np
          close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
        STOP
      ENDIF
c      write(13,*) i_send,' primaries must be sent to neighboring nodes'
      IF(i_send .GT. max_send-2) THEN
         write(name,'(a,i4.4,a)')'ERROR_RUN',node,'.DAT'   
         accesstype = 'append'
         j123=456
         call Open_ASCII_FILE(19,name,accesstype)
         rewind 19
         write(19,*) j123, 
     +    ' Error in Redistribute Primaries'
         write(19,*) 'i_send .GT. max_send-2'
         write(19,*)  i_send,max_send
         close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
        STOP
      ENDIF

      DO knode = 1,n_nodes
        sdispl(knode)  = 0
        sendcount(knode) = 0
      ENDDO

      IF(i_send .EQ. 0) GOTO 555              
      
      norder = np
      CALL  Indexx_i (norder,no_se,indxx)
      knodemin = no_se(indxx(np-i_send+1)) ! lowest node number to which
                                                ! particles are sent
      IF(iDebug) write(13,*) 'knodemin = ',knodemin
      is  = 0
      iss = 0
      knode = knodemin


      DO ii = 1,i_send               ! fill send buffers
        i = indxx(np-i_send+ii)    ! particle number to be sent to
        k =  no_se(i)                    ! node number 
        IF(k .GT. knode)  THEN           ! new node 
          sendcount(knode) = iss         ! the previous node
          DO k_fill = knode+1,k
            sdispl(k_fill)  = is 
          ENDDO   
          knode = k
          iss = 0
        ENDIF
        is = is + 1
        iss = iss + 1
c        write(13,*) ii,i,knode,x(i),y(i),z(i)
!        If(iDebug)
!        write(13,1301) i_par(i),knode,x(i),y(i),z(i),
!     +                        vx(i),vy(i),vz(i)
        x_se(is)   =  x(i)    
        y_se(is)   =  y(i)    
        z_se(is)   =  z(i)    
        vx_se(is)  =  vx(i)   
        vy_se(is)  =  vy(i)    
        vz_se(is)  =  vz(i)   
        pt_se(is)  =  pt(i)  
        wpar_se(is)=  wpart(i)   
        ip_se(is)  =  i_par(i)              
      ENDDO    ! ii

      sendcount(knode) = iss
      IF(knode .LT. n_nodes) THEN
        DO k_fill = knode+1,n_nodes
          sdispl(k_fill)  = is   
        ENDDO
      ENDIF

 555  i_all_send = 0 
      DO i = 1, n_nodes
        i_all_send = i_all_send +  sendcount(i)
      ENDDO       
      IF(i_all_send .NE. i_send ) THEN
         write(name,'(a,i4.4,a)')'ERROR_RUN',node,'.DAT'   
         accesstype = 'append'
         j123=456
         call Open_ASCII_FILE(19,name,accesstype)
         rewind 19
         write(19,*) j123,
     +    ' Error in Redistribute Primaries'
         write(19,*) 'i_all_send .NE. i_send'
         write(19,*)  i_all_send, i_send
         close(19)
               CALL mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
        STOP
      ENDIF

c  send integers with the lengths of all arrays which will be
c  sent (sendcount) and received (recvcount) to/from  all nodes  

      CALL mpi_allgather(sendcount, n_nodes, MPI_INTEGER,
     +                   sendcount_all, n_nodes, MPI_INTEGER,
     +                   MPI_COMM_WORLD, ierr)

      rdispl_new = 0
      DO i = 1, n_nodes
        rdispl(i)    = rdispl_new 
        recvcount(i) = sendcount_all(irank+1+n_nodes*(i-1))
        rdispl_new   = rdispl(i) + recvcount(i)
      ENDDO

c      WRITE(13,'(''sendcount: '',100I8)') sendcount
c      WRITE(13,'(''recvcount: '',100I8)') recvcount
c      WRITE(13,'(''sdispl   : '',100I8)') sdispl
c      WRITE(13,'(''rdispl   : '',100I8)') rdispl
c      WRITE(13,*) ' '

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
      WRITE(13,*) ' transfer time total:',t2-t1 


c Add new primaries to old ones 
C 
      i_rec = 0
      DO k = 1, n_nodes
        i_rec = i_rec + recvcount(k)
      ENDDO
c      WRITE(13,*) i_send, ' of ', np_sb(node),' primaries sent '
c      WRITE(13,*) i_rec, ' primaries added to ', i_new       
      n_recv_prim = i_rec               ! received primaries 

      DO k = 1, i_rec
        i = i_new + k
        x(i)   =  x_re(k)    ! now adding primaries sent from other nodes
        y(i)   =  y_re(k)
        z(i)   =  z_re(k)
        vx(i)  =  vx_re(k)
        vy(i)  =  vy_re(k)
        vz(i)  =  vz_re(k)
        pt(i)  =  pt_re(k)
        wpart(i) =  wpar_re(k)
        i_par(i) =  ip_re(k)
!        If(iDebug)
!         WRITE(13,1302) i,i_par(i),x(i),y(i),z(i),
!     +                          vx(i),vy(i),vz(i)
      ENDDO  
 1301 FORMAT( I8,' --> ',I4,3F8.3,3X,3g12.4)
 1302 FORMAT(2I8,' <-- ',3F8.3,3X,3g12.4)

      n_par_in_box = 0
      DO i = 1, n_nodes
        n_par_in_box = n_par_in_box + np_sb(i)
      ENDDO 
c      WRITE(13,*) 'np_sb = ', np_sb
      WRITE(13,*) 'n_par_in_box = ', n_par_in_box

      n_new = i_new + n_recv_prim 
      np_sb(node) = n_new
      WRITE(13,'(" new # of particles : ",T25,16i9,(/T25,16i9))')
     &             np_sb(node)

c  exchange of information about the number of particles      

      mpi_send = 1
      CALL mpi_allgather(n_new, mpi_send, MPI_INTEGER,
     +                   np_sb, mpi_send, MPI_INTEGER,
     +                   MPI_COMM_WORLD, ierr)
      WRITE(13,*) 'new np_sb = ', np_sb

      n_par_in_box = 0
      DO i = 1, n_nodes
        n_par_in_box = n_par_in_box + np_sb(i)
      ENDDO
      WRITE(13,*) 'n_par_in_box = ', n_par_in_box

      IF(n_par_in_box .NE. n_total_part ) THEN
        write(13,*) ' ERROR: particles lost, n_par_in_box .NE. np',
     +         n_par_in_box, n_total_part
      ENDIF

      return
      end


c     -----------------------------------
      SUBROUTINE indexx_i(n,arr,indxx)
c     -----------------------------------
c
c  purpose: order particles by node number
C    Indexing  Sect. 8.4. of Numerical Recipis
C    Changed to integer field

      INTEGER n,indxx(n),M,NSTACK
      INTEGER arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      INTEGER a

      do 11 j=1,n
        indxx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indxx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indxx(i)).le.a)goto 2
            indxx(i+1)=indxx(i)
12        continue
          i=0
2         indxx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indxx(k)
        indxx(k)=indxx(l+1)
        indxx(l+1)=itemp
        if(arr(indxx(l+1)).gt.arr(indxx(ir)))then
          itemp=indxx(l+1)
          indxx(l+1)=indxx(ir)
          indxx(ir)=itemp
        endif
        if(arr(indxx(l)).gt.arr(indxx(ir)))then
          itemp=indxx(l)
          indxx(l)=indxx(ir)
          indxx(ir)=itemp
        endif
        if(arr(indxx(l+1)).gt.arr(indxx(l)))then
          itemp=indxx(l+1)
          indxx(l+1)=indxx(l)
          indxx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indxx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indxx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indxx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indxx(i)
        indxx(i)=indxx(j)
        indxx(j)=itemp
        goto 3
5       indxx(l)=indxx(j)
        indxx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.


!c     -----------------------------------
!      subroutine Test () 
!c     -----------------------------------
!c
!c  purpose: test which particles are on each node

!c  TEST EXAMPLE: the node sends all particles to all nodes  !!!!!!!

!      include 'a_tree.h'
!      include 'a_mpi.h'
!      include 'a_control.h'

!      REAL*8 xmin,xmax,ymin,ymax,zmin,zmax
!      REAL*8 vxmin,vxmax,vymin,vymax,vzmin,vzmax
!      INTEGER indx(np)

!c      write(13,*) '   test output begin *******************************'

!      write(13,*) 'n_all, n_recv_larg, n_recv_prim,n_refin'
!      write(13,*) n_all, n_recv_larg, n_recv_prim,n_refin
!      node = irank + 1


!       xmin = 1.E9
!       xmax = -1.E9
!       ymin = 1.E9
!       ymax = -1.E9
!       zmin = 1.E9
!       zmax = -1.E9
!      vxmin = 1.E9
!      vxmax = -1.E9
!      vymin = 1.E9
!      vymax = -1.E9
!      vzmin = 1.E9
!      vzmax = -1.E9

!ccc      write(*,*) 'here we are',n_all
!      DO i = 1, n_all
!           IF(i_par(i) .EQ. 100) i100 = i
!           IF(i_par(i) .EQ. 1000) i1000 = i
!           IF(i_par(i) .EQ. 10000) i10000 = i
!           IF(i_par(i) .EQ. 30000) i30000 = i
!           xmin = min(xmin,x( i))
!           xmax = max(xmax,x( i))
!           ymin = min(ymin,y( i))
!           ymax = max(ymax,y( i))
!           zmin = min(zmin,z( i))
!           zmax = max(zmax,z( i))
!          vxmin = min(vxmin,vx( i))
!          vxmax = max(vxmax,vx( i))
!          vymin = min(vymin,vy( i))
!          vymax = max(vymax,vy( i))
!          vzmin = min(vzmin,vz( i))
!          vzmax = max(vzmax,vz( i))
!      ENDDO

!      write(13,*) n_all,' particles will be tested: '
!      write(13,*) ' xmin,  xmax: ',xmin,xmax
!      write(13,*) ' ymin,  ymax: ',ymin,ymax
!      write(13,*) ' zmin,  zmax: ',zmin,zmax
!      write(13,*) 'vxmin, vxmax: ',vxmin,vxmax
!      write(13,*) 'vymin, vymax: ',vymin,vymax
!      write(13,*) 'vzmin, vzmax: ',vzmin,vzmax

!c      write(13,*) 'i  x y z vx vy vz ' 
!c      ii = i100
!c      write(13,*)  istep,ii, 'coord ',x(ii), y(ii),z(ii)
!c      write(13,*)  istep,ii, 'veloc ',vx(ii),vy(ii),vz(ii)
!c      ii = i1000
!c      write(13,*)  istep,ii, 'coord ',x(ii), y(ii),z(ii)
!c      write(13,*)  istep,ii, 'veloc ',vx(ii),vy(ii),vz(ii)
!c      ii = i10000
!c      write(13,*)  istep,ii, 'coord ',x(ii), y(ii),z(ii)
!c      write(13,*)  istep,ii, 'veloc ',vx(ii),vy(ii),vz(ii)
!c      ii = i30000
!c      write(13,*)  istep,ii, 'coord ',x(ii), y(ii),z(ii)
!c      write(13,*)  istep,ii, 'veloc ',vx(ii),vy(ii),vz(ii)

c      write(13,*) '   test output end *********************************'

!      return
!      end

c     -----------------------------------
      subroutine Test_coord () 
c     -----------------------------------
c
c  purpose: test position of particles which are shifted after the
c           second  integration step

      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'

c      PARAMETER (n_test = 32)
      PARAMETER (n_test = 12)
      INTEGER   i_test (n_test)

c      DATA i_test / 117,  5370,  6772,  8646, 10302, 10866, 11201,
c     +     11323, 12150, 12350, 13446, 13874, 14037, 14328, 17001,
c     +     19092, 19155, 19363, 20296, 21483, 22552, 24120, 24606,
c     +     24681, 27316, 27399, 27817, 28021, 28841, 30356, 31239,
c     +     32568 /
      

      DATA i_test / 13468, 16889, 16892, 18553, 20923, 20924, 21988,
     +      22948, 30304, 30326, 31195,  32533 /       
      

      DO i = 1, n_all
        DO k = 1, n_test
          IF( i_par(i) .EQ. i_test(k) ) THEN
              write(1717,100) i_par(i),x(i),y(i),z(i),vx(i),vy(i),vz(i) 
            GOTO 1234
          ENDIF
        ENDDO
 1234   CONTINUE
      ENDDO
 100  FORMAT(I8,3F8.2,4X,3F8.2)


      return
      end








