c     ------------------------
      subroutine Read_ART_MPI_Inp ()
c     ------------------------
c
c     purpose: read parameters from ART_MPI.input and distribute them
c 
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'
      CHARACTER*120    tmp1,tmp2,tmp3
      
      IF ( irank .eq. root ) then

         OPEN(9,file='ART_MPI.input', STATUS = 'OLD')
         READ(9,'(A)') name    ! directory for output
         READ(9,*) que1, que2   ! maximum time in queue h min 
         queue_time = ( 60.*que1 + que2 ) * 60.   ! seconds
         write(13,*) ' time to run(secs)=',queue_time
         READ(9,*) mstep       ! number of time steps to run
         READ(9,*) iostep      ! dump step 
         READ(9,*) n_save      ! 
         IF (n_save .gt. 0) THEN
           DO i = 1, n_save
             READ(9,*) asave(i)
           ENDDO
         ENDIF
         CLOSE(9)  
      Endif
!      write(81,*) ' done init reading '  
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

c... distribute input data

      is = 120    ! must be the length of name in a_control.h
      CALL MPI_BCAST(name,is,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
      is = 1
      CALL MPI_BCAST(queue_time,is,MPI_REAL,root,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mstep,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iostep,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(n_save,is,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      IF (n_save .gt. 0) THEN
        DO i = 1, n_save
          CALL MPI_BCAST(asave(i),is,MPI_REAL,root,MPI_COMM_WORLD,ierr)
        ENDDO
      Endif   

!      write(81,*) ' done init reading 2 n_save=',n_save  
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP


      IF ( irank .eq. 1 ) write(81,*) TRIM(name),
     +         '   mstep=',mstep,' iostep=',iostep,' n_save=',n_save


      write(name1,'(a)')'/RUN'
      write(directory,'(a)')TRIM(name)
      write(name,'(a,i4.4,a)') 'LOGS/ART_MPI.',irank+1,'.output'

      accesstype = 'append'
C  error file  ALWAYS write(19,........    !!!
C  output data ALWAYS write(13,........    !!!
      write(81,*) ' done init reading 2 Open file 13'  
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP


      call Open_ASCII_FILE(13,name,accesstype)
      write(81,*) ' done init reading 3 Opened file 13. Write'  
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP


      write (13,*) ' Start '
      write (13,*) 'mstep  = ',mstep
      write (13,*) 'iostep = ',iostep
      write (13,*) 'n_save = ',n_save
      write(13,*) ' name=',TRIM(name)
      write(13,*) ' directory=',TRIM(directory)
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

      return
      end

c     ----------------------
      subroutine Open_Log ()
c     ----------------------
c
c     purpose: open files used to log the current run
c    
c     the files are opened to be appended, not rewritten
c

      include 'a_tree.h'
      include 'a_control.h'
      include 'a_mpi.h'

      write(name,'(a,i4.4,a)') 'LOGS/timing.',irank+1,'.log'
      accesstype = 'append'
      call Open_ASCII_FILE(40,name,accesstype)

      write(name,'(a)') 'LOGS/load_balance.log'
      accesstype = 'append'
      if (irank .eq. root) 
     +          call Open_ASCII_FILE(41,name,accesstype)

      write(name,'(a,i4.4,a)') 'LOGS/run.',irank+1,'.log'
      accesstype = 'append'
      call Open_ASCII_FILE(50,name,accesstype)

      return
      end

c     -----------------------
      subroutine Close_Log ()
c     -----------------------
c
      close (40)
      close (41)
      close (50)
      return
      end

c     -----------------------
      subroutine Init_Control_Files ()  
c     -----------------------
c
c     initializes control_files for queue management
c 
 
      include 'a_control.h'
      include 'a_numbers.h'


      return
      end
      
c     -----------------------
      subroutine check_remaining_time () 
c     -----------------------
c
c     check remaining time 
c 

      include 'a_control.h'

      delta_time = seconds() - time1    ! time of last step 
      time1      = seconds() - time_start
      rest       = queue_time - time1
      write (13,*)  ' time remaining=',rest,delta_time 
      IF(rest .LT. delta_time*1.2) THEN
        Stop_run= .true. 
        write (13,777)  istep, delta_time, rest
 777    FORMAT(' Stop_run = .true.: istep =',I4,
     +                  ' delta_time [s] = ', F10.1,
     +                  ' rest [s] = ', F10.1 )
      ENDIF

      write (13,*) ' Stop_run=',Stop_run  
      return
      end

c     ------------------------------------------------------------
      subroutine Open_ASCII_File ( iunit , name , accesstype )
c     ------------------------------------------------------------
c
c     accesstype is the string containing the type of file access
c                (e.g., append)
c

      include 'a_control.h'

c construct here the different MPI files 
 
      name1 = TRIM(directory)//'/'//TRIM(name)
c.... open takes into account different syntax on AIX or others
c.... position='append' - for AIX; access = 'append' for others


c      open ( iunit , file = name1(in1:jn1) ,
c     &               position = accesstype(ia1:ja1) )
c      write(13,*) ' open file=',iunit,in1,jn1
c      write(13,*) ' open file=',id1,jd1
c      write(13,*) ' open file=',directory(id1:jd1)
c      write(13,*) ' open file=',iunit,name1(in1:jn1)
      open ( iunit , file = TRIM(name1) ,
     &               access = accesstype, STATUS='unknown' )
      
      return
      end

c     -----------------------------------
      subroutine Read_Control()
c     -----------------------------------
c
c     purpose: reads MPI control file PMPICrd.DAT
c

      CHARACTER*120   tmp1
 
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'

      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))

      node = irank + 1
      i_chunk = 50
!      Do !!!!!!
      write(tmp1,'(a,i4.4,a)') '/RUN/PMPICrd.',node,'.DAT'
c      write(13,*)   directory(id1:jd1)//tmp1(it11:jt11)
      open ( 3,file=TRIM(directory)//TRIM(tmp1),
     +        form = 'unformatted', STATUS = 'OLD')

c.... read control information and check whether it has proper structure
      write (81,*) '   Read Control file'
      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecs,Nseed,Om0,Oml0,hubble,Wp5
     &                   ,Ocurv,extras
      write (81,*) '   Read Control one'
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

      write (13,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nseed,Om0,Oml0,hubble,
     +                  Nspecs

      If(extras(100).gt.0.)Then
         Box =extras(100)
         write(13,*)'  Box size =',Box,' Mpc/h'
      endif
      if ( Nspecs.ne. Nspecies ) then
        write (13,*)
     &      ' Nspecies in PARAMETER and in TAPE-FILE are different:'
         write (13,*) ' Nspecies=',Nspecies,' not equal ',Nspecs
         STOP
      endif
      write (81,*) '   Read Control one'
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

      l_all = 0
      DO i = 1, Nspecs
        n_part = lspecies(i)-l_all
        write(13,*)  n_part,' particles, weight = ',wspecies(i)
        l_all = l_all + n_part
      ENDDO
      wpar_min = wspecies(1)*1.0001  ! increased for test .LE.  
      READ (3) n_divvX,n_divvY,n_divvZ, n_nodess
      IF(n_divvX.NE.n_divX) STOP 'incorrect ndiv or nndivX'
      IF(n_divvY.NE.n_divY) STOP 'incorrect ndiv or nndivY'
      IF(n_divvZ.NE.n_divZ) STOP 'incorrect ndiv or nndivZ'
      IF(n_nodes.NE.n_nodess) STOP 'incorrect n_nodes or n_nodess'
      READ(3)   l_divX,l_divY,l_divZ
c      WRITE(13,'(''l_div: '',100I8)') l_divX,l_divY,l_divZ
      READ(3)   (np_sb(inode) ,inode = 1,n_nodes)
      n_total_part = l_all 
      write(13,*) 'total number of particles:',n_total_part
      l_check = 0
      DO i = 1,n_nodes
        l_check = l_check + np_sb(i)
      ENDDO
      IF( l_check .NE. n_total_part) THEN
         write(name,'(a,i4.4,a)') 'ERROR_RUN.',node,'.DAT'
         accesstype = 'append'
         j123=999
         call Open_ASCII_FILE(19,name,accesstype)
         rewind 19
         write(19,*) j123, 
     +           'wrong particle numbers '
        write(19,*) 'l_check .NE. n_total_part',l_check,n_total_part
        write(19,'(20I9)') np_sb
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
        CALL mpi_finalize(ierr)
        STOP
      ENDIF
      write(13,*) 'number of particles on all nodes:',l_check
!      write(13,'(16I9)') np_sb
 
      close (3)

      return

100   FORMAT(1X,'Header=>',A45,/
     +   1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.6,' Step=',F8.3,/
     +   1X,' I =',I4,' WEIGHT=',F8.3,/
     +   1X,'   TINTG =',G10.3,'  EKIN =',G10.3,' EKIN1 =',G10.3,
     +          ' EKIN2 =',G10.3,/
     +   1X, '   AU0  = ',G10.3,' AEU0 = ',G10.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nseed=',I8,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,
     +            1x,' number of species = ',I10)
      end


c     -----------------------------------
      subroutine Read_Particles ()
c     -----------------------------------
c
c     purpose: read  particle coordinates, momenta and time steps
c     of  file # node
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'

      CHARACTER*120   tmp1
      REAL*8 xmin,xmax,ymin,ymax,zmin,zmax

      node = irank + 1

      write(tmp1,'(a,i4.4,a)') '/RUN/PMPICrs.',node,'.DAT'

      write(13,*) ' Open file= '
      write(13,*)   TRIM(directory)//TRIM(tmp1)
      OPEN(9,file = TRIM(directory)//TRIM(tmp1),
     +                    FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')

      write(13,*) ' opened coords file'


      READ (9) nb,ngg,istep
      IF(ng.NE.ngg) write (13,*) ' wrong ng or ngg'
      IF(ng.NE.ngg) write (13,*) ' node=',node,' ng=',ng,' ngg=',ngg
c      IF(ng.NE.ngg) call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
c      IF(ng.NE.ngg) STOP 'wrong ng or ngg'
      IF(nb.GT.np) Then
        write (13,*) ' wrong np=',np
        write (13,300) node,nb,ng,istep
        CALL mpi_finalize(ierr)
        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
        STOP 'increase number of particles np in a_setup.h'
      EndIf
      READ(9) l_div_x1,l_div_x2,l_div_y1,l_div_y2,
     +               l_div_z1,l_dix_z2
      write(13,300) node,nb,ng,istep
 300  FORMAT('node = ',I3, ' nb = ',I8, '  ng = ',I4, ' istep = ',I4)
      write(13,*)   ' borders of subbox:  l_diff (x1,x2 y1,y2 z1,z2)'
      write(13,*)   l_div_x1,l_div_x2,'   ',l_div_y1,l_div_y2,'   ',
     +               l_div_z1,l_dix_z2
      write(13,*) ' Start reading coords Np=',nb
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!       CALL mpi_finalize(ierr)
!       Stop


      READ(9) ( x( i), i = 1,nb)
      write(13,*) ' x is done'
      READ(9) ( y( i), i = 1,nb)
      write(13,*) ' y is done'

      READ(9) ( z( i), i = 1,nb)
      write(13,*) ' z is done'

      write(13,*) ' Start reading vels Np=',nb
      
      READ(9) (vx( i), i = 1,nb)
      READ(9) (vy( i), i = 1,nb)
      READ(9) (vz( i), i = 1,nb)
      write(13,*) ' vs are done'

      READ(9) (pt( i), i = 1,nb)
      write(13,*) ' pt is done'

      READ(9) (wpart( i), i = 1,nb)
      write(13,*) ' wpar is done'

      READ(9) (i_par( i), i = 1,nb)
      write(13,*) ' i_par is done'

!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!       CALL mpi_finalize(ierr)
!       Stop

      CLOSE(9)

c some testing
         xmin = 1.E9
         xmax = -1.E9
         ymin = 1.E9
         ymax = -1.E9
         zmin = 1.E9
         zmax = -1.E9
         ptmin = 1.E9
         ptmax = -1.E9
         wparmin = 1.E9
         wparmax = -1.E9
         i_par_min = 1024**3+1
         i_par_max = -1
        DO i = 1, nb
           xmin = min(xmin,x( i))
           xmax = max(xmax ,x( i))
           ymin = min(ymin,y( i))
           ymax = max(ymax ,y( i))
           zmin = min(zmin,z( i))
           zmax = max(zmax ,z( i))
           ptmin = min(ptmin,pt( i))
           ptmax = max(ptmax ,pt( i))
           i_par_min = min(i_par_min,i_par( i))
           i_par_max = max(i_par_max ,i_par( i))
           wparmin = min(wparmin,wpart( i))
           wparmax = max(wparmax ,wpart( i))
       If(abs(vx(i)).gt.1.d5)
     &     write(13,*) ' velocity ERROR X=',vx(i),i
       If(abs(vy(i)).gt.1.d5)
     &     write(13,*) ' velocity ERROR Y=',vy(i),i
       If(abs(vz(i)).gt.1.d5)
     &     write(13,*) ' velocity ERROR Y=',vz(i),i


        ENDDO
        IF(xmin  .LT.  l_div_x1 .OR. xmax .GE. l_div_x2 .OR. 
     +     ymin  .LT.  l_div_y1 .OR. ymax .GE. l_div_y2 .OR. 
     +     zmin  .LT.  l_div_z1 .OR. zmax .GE. l_dix_z2 ) THEN
             write(13,*) ' coordinate ERROR '
             write(13,*)   
     +           ' borders of subbox:  l_diff (x1,x2 y1,y2 z1,z2)'
             write(13,*)   
     +           l_div_x1,l_div_x2,'   ',l_div_y1,l_div_y2,'   ',
     +               l_div_z1,l_dix_z2
             write(13,*)  'xmin/max   ymin/max    zmin/max'
             write(13,*)   xmin,xmax,'    ', ymin,ymax,'       ', 
     +                         zmin,zmax
         ENDIF
      write(13,*) nb,' particles read'
             write(13,*)  'xmin/max   ymin/max    zmin/max'
             write(13,*)   xmin,xmax,'    ', ymin,ymax,'       ', 
     +                         zmin,zmax
             write(13,*)  ' i_par_min i_par_max Npar=1,2,3'
             write(13,'(7i11)') i_par_min,i_par_max, 
     &                  i_par(1),i_par(2),i_par(3)

       call mpi_barrier(MPI_COMM_WORLD,ierr)
c      CALL mpi_finalize(ierr)
c      Stop
       
      return
      end


c     ----------------------------
      subroutine Save_Check ()
c     ----------------------------
c     Andrey's output routine, modified
c
      include 'a_control.h'

c   n_save > 0  : save the time moments read from ART_MPI.input
c                 use aexpn for filename
c   n_save = 0  : save all steps, use  aexpn for file name
c   n_save < 0  : save every |n_save| step, use step number for file name

      if (n_save .gt. 0) then
        do i = 1 , n_save
          if (( aexpn .ge. (asave(i) - 0.5*astep)) .and. 
     &        ( aexpn .lt. (asave(i) + 0.5*astep) ))   
     &        call Save ( 2 ) 
        enddo
      endif
      if (n_save .eq. 0) call Save ( 1 )
      if (n_save .lt. 0) then
        nn_save = - n_save
        if ( mod(istep,nn_save) .eq. 0 ) then 
          call Save ( 1 ) 
        endif
      endif

      if ( mod(istep,iostep) .eq. 0 ) then 
        call Save ( 0 ) 
      endif

      if(aexpn .ge. 1.0000) then 
        name = 'STOP_RUN'
        accesstype = 'asis'
        CALL Open_ASCII_FILE(19,name,accesstype)
        rewind 19
        write(19,'(A)') ' 1   STOP'
        close(19)
        Stop_Run = .true.
      endif

      return
      end
c

c     ------------------------------
      subroutine Save ( iMode )
c     ------------------------------
c     Andrey's output routine, modified
c
c     iMode defines how the output file is named
c        = 0 , filenames are just the standard PMPICrd.DAT, PMPICrs.DAT
c        = 1 , filename is jobname + global step number
c        = 2 , filename is jobname + expansion parameter       
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'

      node = irank +1
      if ( iMode .eq. 0 ) then 
        write(name1,'(a,i4.4,a)')'/RUN/PMPICrd.',node,'.DAT'
        write(name2,'(a,i4.4,a)')'/RUN/PMPICrs.',node,'.DAT'
      endif

      if ( iMode .eq. 1 ) then
        write(name1,10)'/SNAP/PMPICrd.',node,'.',iStep,'.DAT'
        write(name2,10)'/SNAP/PMPICrs.',node,'.',iStep,'.DAT'
 10     format(a,i4.4,a,i4.4,a)
      EndIf
      if ( iMode .eq. 2 ) then
        ii =min(int(10000.*aexpn),9999)
        write(name1,10)'/SNAP/PMPICrd.',node,'.',ii,'.DAT'
        write(name2,10)'/SNAP/PMPICrs.',node,'.',ii,'.DAT'
      EndIf

      call Write_MPI ( name1 , name2 )

      return
      end


c     -----------------------------------
      subroutine Write_MPI (name1,name2)
c     -----------------------------------
c
c     purpose: write particle coordinates, momenta and time steps
c     of the given node
c
      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_control.h'

      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))
      INTEGER ijk(3)

c write control data

      open ( 3,file = TRIM(directory)//TRIM(name1),
     +         form = 'unformatted', STATUS = 'UNKNOWN')

      write (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecs,Nseed,Om0,Oml0,hubble,Wp5
     &                   ,Ocurv,extras
      write (3)  n_divX,n_divY,n_divZ, n_nodes
      write (3)  l_divX,l_divY,l_divZ
      write (3)  (np_sb(inode) ,inode = 1,n_nodes)
      close(3)

c write particles

      node = irank + 1
c      write(13,*) 'irank = ',irank,' node = ',node
c      write(13,*) '   writing ',name2(in21:jn21)
c      write(13,*) ' node = ',node,' nb = ',np_sb(node),' istep = ',istep
      OPEN(3,file = TRIM(directory)//TRIM(name2),
     +       FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')
      write(3) np_sb(node) ,ng,istep  
      CALL Node_to_IJK(node,ijk)
c      write(3) (l_div(k,ijk(k)),l_div(k,ijk(k)+1),k = 1,3) 
      write(3) l_divX(ijk(1)),l_divX(ijk(1)+1),
     &            l_divY(ijk(1),ijk(2)),l_divY(ijk(1),ijk(2)+1),
     &            l_divZ(ijk(1),ijk(2),ijk(3)),
     &            l_divZ(ijk(1),ijk(2),ijk(3)+1)
                                           ! six borders of the subbox
      write(3) ( x(i)   , i = 1,np_sb(node))
      write(3) ( y(i)   , i = 1,np_sb(node))
      write(3) ( z(i)   , i = 1,np_sb(node))
      write(3) (vx(i)   , i = 1,np_sb(node))
      write(3) (vy(i)   , i = 1,np_sb(node))
      write(3) (vz(i)   , i = 1,np_sb(node))
      write(3) (pt(i)   , i = 1,np_sb(node))
      write(3) (Wpart(i), i = 1,np_sb(node))
      write(3) (i_par(i), i = 1,np_sb(node))
      CLOSE(3)
 
      return
      end

