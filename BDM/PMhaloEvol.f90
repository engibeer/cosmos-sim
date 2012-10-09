!---------------------------------------------------------------------------------
!     Read HaloTrace files
!     Extract information on halo evolution
!--------------------------------------------------------------------------------
MODULE   setHalos
PUBLIC      

Integer, parameter  :: Ncats =170          !  Maximum number of  snapshots to analyze

Real, parameter :: Box=250.,   &
                                 Omega0 = 0.27, &  ! Omega_matter
                                 hubble = 0.70
TYPE :: hal
  Real*4         ::   aexpn                      ! expansion parameter
  Integer*4    ::   InCat                        ! pointer to halo in Catshort
  Integer*4    ::   Parent                      ! if not 0, pointer to parent halo in Catshort
  Real         ::   x,y,z,vx,vy,vz        !  coordinates (Mpch units), peculiar velocities(km/s)
  Real         ::  aM,Rm,circV            ! Mass(Msunh), Vir radius(kpch), max circ Velocity
End TYPE hal

TYPE  :: haloEvolution
   Integer :: N                                   ! number of moments for which halo was traced
   Type(hal) :: h(Ncats)
End TYPE haloEvolution

Type(hal) ::     halo   
Contains
!---------------------------------------------------------
!                         Read 'iC' halo catalog 
!
      SUBROUTINE ReadHalos(iC)
!--------------------------------------------------------------
Character*120   :: listname, Line, Txt*5, Txt2*50
Type(hal)                      :: halo2
Type(haloEvolution)     :: History
                              !                                                                                 Open file
      write(listname,'(a,i3.3,a)') 'HaloTree/HaloTrace.',iC,'.dat'
      Open( iC+100,file=TRIM(listname),Status='UNKNOWN',form='unformatted')

   read(iC+100)      ! first line is empty
   nHalos = 0
   Do
      read(iC+100,iostat=iStat) ih
      If(iStat.ne.0)exit
      Do
         read(iC+100,iostat=iStat) halo2                  ! there is a bug in reading these data
         backspace(iC+100)                                       ! read it twice to get around the bug
         read(iC+100,iostat=iStat) halo2%aexpn   !     read aexpn again
         If(iStat.ne.0)exit
         !write(*,'(i9,i8,f9.4,2i9,3f9.4,3f9.2,1pg12.4,0p2f8.2)')ih,nHalos, halo2
      EndDo
      nHalos = Nhalos +1
      !if(nHalos> 5)stop
   EndDo

      write(*,'(a,i8,a,f7.4,a,i3)') ' Read Halos=',nHalos
      close(100+iC)
      Return

      End SUBROUTINE ReadHalos
    End module setHalos
!-------------------------------------------------
!            
!            
PROGRAM HaloCrossIdentity

use     SetHalos
Character*120                                     :: file1,listname,txt*8     
!     ------------------------------------------        Read data
      Open(1,file='HaloCrossList.dat')

      Call CPU_time(t0); call system_clock(iT0,iTc,ib)
 
   write(*,*)
   Call CPU_time(t1); call system_clock(iT2,iTc,ib) ; write(*,*) '          time for Init =',t1-t0,float(iT2-iT0)/iTc

                              ! ---------------------------------------------   Main loop 
Do iC =1, 100
   Call ReadHalos(iC)
   Call CPU_time(t1); call system_clock(iT2,iTc,ib) ; write(*,*) '          time  =',t1-t0,float(iT2-iT0)/iTc
EndDo

stop
end Program HaloCrossIdentity
