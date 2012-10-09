c     =====================================================================
c                                                                         .
c       ART Version 3: Advancing control variables to the next time step  .
c                                                                         .
c                by Andrey Kravtsov and Anatoly Klypin (1997)             .
c                                                                         .
c     =====================================================================

c     ---------------------
      subroutine Advance ()
c     ---------------------
c
c     advances time variables, 
c     computes kinetic and potential energies
c     and energy conservation error
c
      include 'a_tree.h'
      include 'a_control.h'

      double precision epotb,dspl
      double precision ekinb(nspecies)                   ! kinetic energy 
      double precision Ndisplace(0:50)

      ahalf = aexpn + 0.5*astep
      dconst = astep/ sqrt (ahalf * (Om0+Oml0*ahalf**3)) /
     &     ahalf  ! constant to convert velocity to
                                ! displacement in units of zero-level cell
      dspl    = 0.
      ndspl   = 0
      dspl_max= 0.
      Ndisplace=0.
!      write(13,'(a,3g14.5)') ' dconst=',dconst, astep, aexpn
!      write(13,*) Om0,Oml0,ahalf
        ekinb = 0.
        epotb = 0.


!      do i=1,n_all
!         write(13,'(3g14.5,3g14.5,i10)')vx(i),vy(i),vz(i),
!     &             x(i),y(i),z(i),i
!      enddo

      IF(mpi_yes .NE. 1234567 ) THEN  
        do ic0 = 1 , nspecies 
          do ic1 = nsp(ic0,1) , nsp(ic0,2)
            v2   = vx(ic1)**2 + vy(ic1)**2 + vz(ic1)**2
            ekinb(ic0) = ekinb(ic0) + wpar(ic0) * v2
            epotb = epotb + wpar(ic0)*pdummy(ic1)
            vv = dconst*sqrt(v2)
            iv = MAX(MIN(INT(10.*vv),50),0)            ! units of 1/10 of a cell
!            if(iv<0)Then
!               write(13,'(a,2i10,6g12.4)')' Advance:',
!     &          ic1,ic0,x(ic1),y(ic1),y(ic1),vx(ic1),vy(ic1),vz(ic1)
!             iv =50
!            endif
            Ndisplace(iv) =Ndisplace(iv) +1.
            dspl = dspl + vv
            ndspl = ndspl +1
            dspl_max =max(vv,dspl_max)
          enddo
        enddo
      ELSE
        ic0 = 1 
          do ic1 = 1, n_all
            v2   = vx(ic1)**2 + vy(ic1)**2 + vz(ic1)**2
            ekinb(ic0) = ekinb(ic0) + Wpart(ic1) * v2
            epotb = epotb + Wpart(ic1)*pdummy(ic1)
            vv = dconst*sqrt(v2)
            iv = MIN(INT(10.*vv),50)            ! units of 1/10 of a cell
            Ndisplace(iv) =Ndisplace(iv) +1.
            dspl = dspl + vv
            ndspl = ndspl +1
            dspl_max =max(vv,dspl_max)
          enddo
      ENDIF

      nn =0
      iv   =0
       do i=50,0,-1
         nn =nn +Ndisplace(i)
         if(nn.gt.100)Then
            iv =i
            goto 5
         endif 
      enddo 
! 5    write (13,*)ndspl
!      write (13,*)Ndisplace
!       call mpi_barrier(MPI_COMM_WORLD,ierr)
!!        call mpi_abort(MPI_COMM_WORLD,ierr1,ierr2)
!        CALL mpi_finalize(ierr)
!        STOP


 5    write (13,100) ndspl,dspl/ndspl,dspl_max,iv,nn
 100  format(i8,' mean displacement=',g11.3,' max=',g11.3,
     &         ' 100_fast=',i10,' nfast=',i10,' units: 1 cell')
      write (13,110) Ndisplace
 110  format(5x,' N with given displacement: 1/10 of cell',/6(10g12.3/))
      enkin = zero
      do ic0 = 1 , nspecies 
        ekinb(ic0) = ekinb(ic0) / 2.0 / (aexpn + 0.5*astep)**2 
        enkin      = enkin + ekinb(ic0)
      enddo
        ENpot = epotb / 2.0
 
      istep = istep + 1
      aexpn = aexpn + astep
c
c.... energy conservation control
c
      IF ( istep .eq. 1 ) THEN
        Ekin1 = Ekin
        Ekin2 = 0.
        Ekin  = ENkin
        au0   = aexp0 * ENpot
        aeu0  = aexp0 * ENpot + aexp0 * (Ekin + Ekin1) / 2.0
        TINTG = 0.

        write(50,*) HEADER
        write(13,40) istep , aexpn , Enkin , ENpot , au0 , aeu0
        write(50, *) ' istep ' , ' aexpn ' , 
     &               ' Ekin ' , ' Ekin_s ' , ' Ekin_l ' , 'PE'
        write(50,40) istep , aexpn , Enkin , 
     &               (Ekinb(ic0),ic0=1,nspecies) , Enpot
        write(50, *) ' istep ' , ' aexpn ' , ' MaxLevel ' , ' ncell ',
     &               ' Ekin ' , ' Ekin1 ' , ' Ekin2 ' , 
     &               '  PE  ' , ' TINTG '

      ELSE
        Ekin2 = Ekin1
        Ekin1 = Ekin
        Ekin  = ENkin
        TINTG = TINTG + astep * (Ekin1 + (Ekin - 2.*Ekin1 + Ekin2)/24.)
        error = ((aexpn - astep) * ((Ekin  + Ekin1) / 2.0 + ENpot) - 
     &            aeu0  + TINTG) / ((aexpn - astep) * ENpot)

        ncell = noct * nchild + ncell0

        call Get_MaxLevelNow ()

        write  (13,50)  istep , aexpn , MaxLevelNow , ncell , 
     &                  Enkin ,  (Ekinb(ic0),ic0=1,nspecies) , 
     &                  Enpot , error , TINTG
        write (50,50)  istep , aexpn , MaxLevelNow , ncell , 
     &                  Enkin ,  (Ekinb(ic0),ic0=1,nspecies) ,
     &                  Enpot , error , TINTG

      ENDIF

c.... check if format conforms with the number of species
c     namely check number of e12.4

40    format (i5 , 1x , f10.4 , 1x , 4(1x,e12.4) )
50    format (i5 , 1x ,  f7.5 , 1x , i2 , 1x ,  i10 , 6(1x,E13.5) )


      return
      end

