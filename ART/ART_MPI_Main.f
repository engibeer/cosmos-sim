C  ART_MPI_Main.f
c
c  shell program around  
c     ====================================================
c
c          Adaptive Refinement Tree (ART) N-body solver
c
c                 Version 3 - February 1997
c
c       Andrey Kravtsov, Anatoly Klypin, Alexei Khokhlov
c
c     ====================================================

      include 'a_tree.h'
      include 'a_mpi.h'
      include 'a_lg.h'
      include 'a_control.h'

      Common /AuxTime/t1,t2,t3,t4,t5

      integer nMove
      integer mtot

      CHARACTER*120    tmp1,tmp2,tmp3

      mpi_yes = 1234567


      write (*,*)  ' <=== IN MPI_ART ==>'

      CALL mpi_init(ierr)
      CALL mpi_comm_size(MPI_COMM_WORLD, mpisize, ierr )
      CALL mpi_comm_rank(MPI_COMM_WORLD, irank, ierr )

      write (*,*)  ' initiated MPI '
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP
      node = irank +1
      Do i = 1,1
       write(tmp1,'(a,i1.1,a,i3.3,a)')'/home1/aklypin/Box120c/Accel.',
     &                                 i,'.',node,'.dat'
       open(80+i,file=tmp1)
      EndDo

      IF(mpisize .NE. n_nodes) THEN
        write(*,*) 'mpisize .NE. n_nodes',mpisize,n_nodes
       call mpi_barrier(MPI_COMM_WORLD,ierr)
!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
        CALL mpi_finalize(ierr)
        STOP
      ENDIF
        LengthC = 3*nctot+moct*(nneigh+ndim)
        LengthB = 8*N_m1_lg**3
        iStepB      = 1    !  max step to move node boundaries

      CALL Read_ART_MPI_Inp ()
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1)
        do ic1 = 1 , mcell
          var(1,ic1) = -1.0
          var(2,ic1) = -1.0
          pot(ic1)   = zero
        enddo
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

c.... read in initial conditions, initialize variables and linked lists
      call Read_Control()
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

      call Read_Particles()
      write (13,*) ' Done Read_Particles' 
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

      call Init_Control_Files    ()
      call Init_Parameters       () 
      dBuffer = 0.25             ! set width of refinement buffer manually
      Stop_Run = .false. 
      time_start = seconds()
      time1=0.
      node = irank + 1

      write(13,*)' Initialization is done node=',node
      close (13)
c.... main loop over mstep (read from input) integration steps
      DO ijkl = 1, 1  !!!!mstep
c            gradually increase density threshold for
c            splitting the zero level. This is to
c            save memory.
!        If(iStep.eq.16.or.iStep.eq.32.or.iStep.eq.64.or.iStep.eq.96)Then
!         If(ijkl == 1 )Then
!       If(iStep.eq.72.or.iStep.eq.86.or.iStep.eq.120)Then
!            write(*,*) ' --- Increase step:',astep,iStep
!            astep = astep/1.5
!            do Level = MinLevel , MaxLevel 
!               dAexp(Level) = astep / 2.0**iTimeBin(Level) 
!           enddo
!           damin = astep / 2.0**iTimeBin(MaxLevel)
!         EndIf
!       If(iStep.eq.230)Then
!            !write(*,*) ' --- Decrease step:',astep,iStep
!            astep = astep/2.
!            do Level = MinLevel , MaxLevel 
!               dAexp   (Level) = astep / 2.0**iTimeBin(Level) 
!           enddo
!           damin = astep / 2.0**iTimeBin(MaxLevel)
!         EndIf

         dsplit =0.
         If(aexpn.lt.0.08)Then
            nsplit(MinLevel) = 2
            nsplit(MinLevel+1) = 2
            nsplit(MinLevel+2) = 0
            nsplit(MinLevel+3) = 2
            nsplit(MinLevel+4) = 2
            dsplit = 0.6
         Else
         If(aexpn.lt.0.1)Then
            nsplit(MinLevel) = 2
            nsplit(MinLevel+1) = 2 
            nsplit(MinLevel+2) = 0
            nsplit(MinLevel+3) = 2
            nsplit(MinLevel+4) = 2
            dsplit = 0.8
         Else
         If(aexpn.lt.0.12)Then
            nsplit(MinLevel) = 2
            nsplit(MinLevel+1) = 2 
            nsplit(MinLevel+2) = 1
            nsplit(MinLevel+3) = 2
            nsplit(MinLevel+4) = 2
            dsplit = 0.
         Else
         If(aexpn.lt.0.40)Then
            nsplit(MinLevel) = 2
            nsplit(MinLevel+1) = 2 
            nsplit(MinLevel+2) = 3
            nsplit(MinLevel+3) = 4
            nsplit(MinLevel+4) = 4
         Else
            If(aexpn.lt.0.5)Then
               nsplit(MinLevel) = 2
               nsplit(MinLevel+1) = 3
               nsplit(MinLevel+2) = 4
               nsplit(MinLevel+3) = 5
               nsplit(MinLevel+4) = 5
            Else
               nsplit(MinLevel) = 4
               nsplit(MinLevel+1) = 5
               nsplit(MinLevel+2) = 5
               nsplit(MinLevel+3) = 5
               nsplit(MinLevel+4) = 5
            EndIf
         EndIf
         EndIf
         EndIf
         EndIf
         IF(ISTEP.gt.5) iStepB      = 1
         IF(ISTEP.gt.30) iStepB     = 1

         trho(MinLevel)=wpar_min*max(REAL(nsplit(MinLevel)),0.001)-1.0
         trho(MinLevel+1)=wpar_min*max(REAL(nsplit(MinLevel+1)),0.001)
         trho(MinLevel+2)=wpar_min*max(REAL(nsplit(MinLevel+2)),0.001)
         trho(MinLevel+3)=wpar_min*max(REAL(nsplit(MinLevel+3))+dsplit,
     &                     0.001)
         trho(MinLevel+4)=wpar_min*max(REAL(nsplit(MinLevel+4)),
     &                     0.001)

        write(name,'(a,i4.4,a)') 'LOGS/ART_MPI.',node,'.output'
        accesstype = 'append'
        call Open_ASCII_FILE(13,name,accesstype)
        write(13,'("nsplit=",6i6)')(nsplit(i),i=0,4)
        write(13,'("rsplit=",6f8.4)')(trho(i),i=0,4)
!        do i=1,np_sb(node)
!          if(z(i) > 21. .and. z(i) < 24.)
!     &      write(13,'(i9,3f9.4,4g12.4)') i,x(i),y(i),z(i)
!!     &                 vx(i),vy(i),vz(i),wpart(i)
!        enddo
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

c                distribute small and large particles
        CALL Send_Small ()
!        write(13,*) ' Done Send_small '
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

        CALL Send_Large ()
        write(13,*) ' Done Send_large '
        Call Test_Part(node)
        !If(iDebug)call Test ()
!        If(ijkl ==2)Then
!        do i=1,n_all
!!          if(y(i) > 38. .and. y(i) < 46.)
!           write(13,'(i9,3f9.4,4g12.4)') i,x(i),y(i),z(i),
!     &                 vx(i),vy(i),vz(i),wpart(i)
!        enddo
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP
!        EndIf
        start_time = MPI_WTIME()
        t1 =0.
        t2 =0.
        t3 =0.
        dummy =seconds()  
!        Do i=1,n_all
!         if(abs(vx(i))>1.d5.or.abs(vy(i))>1.d5.or.abs(vz(i))>1.d5)
!     &      write(13,'(a,6g14.5,i10,2i10)')
!     &      ' Vy error:',x(i),y(i),z(i),vx(i),vy(i),vz(i),i,
!     &      np_sb(node),n_refin
!        EndDo
c                 init arrays
        call Init_Arrays           ()
        call Init_Tree             ()
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

        call LL_Construct          ()  ! construct linked list of particles
        call Open_Log              ()

c                rebuild refinements if code is restarted     
        do ic1 = 1 , 10 
          CPU(ic1) = 0.0
        enddo
        dummy =seconds()
        mtot  = 1
        ncell = noct * nchild + ncell0

        if ( istep .ge. 0 ) then 
          do while ( (mtot .ne. nil) )
            call Get_MaxLevelNow ()
            ncell = noct * nchild + ncell0
            call Timing ( 1 , -1 )
            call Assign_Density ( MaxLevelNow , MaxLevelNow )
            call Timing ( 1 ,  1 )
            write(13,*) ' ncells=',ncell ,' MaxLevel=', MaxLevelNow

            call Timing ( 8 , -1 )
            call Modify ( mtot, MaxLevelNow, MaxLevelNow, MaxLevel, 1)
            call Timing ( 8 ,  1 )
            ncell = noct * nchild + ncell0


          enddo
        endif
            write(13,*) ' Density is done.'

        t2 =t2 + seconds()
        call Get_MaxLevelNow ()
        do Level = MinLevel , MaxLevelNow 
          call LL_Update ( Level , MinModify , MaxModify )
        enddo
        t3 =t3 +seconds()

        call Move ()
            write(13,*) ' move is done. '
c                        advance time variables, compute energies
!           call Timing ( 10 ,  -1 )
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

        call Advance
           call Timing ( 10 ,  1 )
            write(13,*) ' advance is done. '

c                        timing for load balance
        end_time = MPI_WTIME()
c        If(ijkl.le.5)Then
           call LoadBalance2
c        Else
c           call LoadBalance1
c        EndIf
!            write(13,*) ' balance is done. '
!        call mpi_barrier(MPI_COMM_WORLD,ierr)
!c        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

c                        redistribution of primary particles

        call Redistribute_Primaries()
c                        write timing 
        call WriteTiming ()
c                        write output, if necessary
        call Save_Check ()
        call check_remaining_time ()

        if(aexpn .ge. 1.0000) then 
          name = 'STOP_RUN'
          accesstype = 'append'
          call Open_ASCII_FILE(19,name,accesstype)
          rewind 19
          write(19,'(A)') ' 1   STOP'
          close(19)
          Stop_Run = .true.
        endif
        If(iDebug)call Test ()       
        call Close_Log ()

        IF ( Stop_Run ) Then
              iStop=1
        Else
              iStop=0
        EndIf
        write (13,*) ' End of main loop: Stop_run=',Stop_Run
        jStop = 0
        CALL MPI_ALLREDUCE(iStop,jStop,1,MPI_INTEGER,MPI_MAX,
     &                     MPI_COMM_WORLD,ierr)  
        IF ( jStop.eq.1) GOTO 999
        close (13)
!        call mpi_barrier(MPI_COMM_WORLD,ierr)
!c        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP

      ENDDO 

 999  Continue
      CALL Save(0)
      CALL mpi_finalize(ierr)
      stop
      END











