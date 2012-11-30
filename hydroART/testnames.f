c$$$
c$$$      character*256 filename
c$$$      character*256 base, key, a , ext
c$$$
c$$$      aexpn = 0.285
c$$$
c$$$      base='profile'
c$$$      key ='CyVg'
c$$$      write (a, 'F5.3') aexpn
c$$$      ext = '.dat'
c$$$
c$$$c      ljob  = index(jobname1, ' ') - 1
c$$$c      lpath = index(path    , ' ') - 1
c$$$c          fname_gas = path(1:lpath)//'/'//jobname1(1:ljob)//'_a'//
c$$$c     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
c$$$c     &            //'.d '
c$$$
c$$$      lbase = index(base, ' ') -1
c$$$      lkey = index(key, ' ' ) -1
c$$$      la = index(a, ' ' ) -1
c$$$
c$$$      filename=base(1:lbase) // key(1:lkey) // 'a' // a(1:la) // '.dat'
c$$$      write (*,10) filename
c$$$ 10   format(A30)

      character*256 file1, file2
      character*256 inputPATH, outputPATH,  b_a_file, name, momentC
      real moment, moment0
      bar_angle_OLD =0.  ! Initial values
      I_shift=0.
! Global files:

! MODEL C:
!      open (101, file= ! put executable in scratch
!     $ '/callisto3/danielcv/RESULTS/MW_LARGE/bar_angle.G2.dat')
! NEW MODELS:
!      open (101, file= ! put executable in scratch
!     $ '/callisto3/danielcv/RESULTS/NoBULGE2low/bar_angle.G3.dat')
! ZEVS:
!      open (101, file= ! put executable in HOME
!     $     '/zevs5/danielcv/RESULTS/NoBULGE1many/bar_angle.G1.dat')
! BASSI:
      inputPATH='/scratch/scratchdirs/aklypin/Khb/FILES/FILES1/'
      outputPATH='/scratch/scratchdirs/aklypin/RESONANCES/Khb/'
      b_a_file='RESULTS/bar_angle.dat'
      lb_a_file = index ( b_a_file, ' ') -1
      loutputPATH = index(outputPATH, ' ') -1
     
      name=outputPATH(1:loutputPATH)//b_a_file(1:lb_a_file)
      lname = index(name, ' ') -1
c      open (101, file=name(1:lname) )
      write (*,*)  'Writing in...'
      write (*,*)   name(1:lname)


      moment0 =0.1600           ! First file. BASSI
      Nrange=1                  ! Number of files
      Ndelta = 0.0001           ! increment

         moment=moment0+Ndelta*(i-1) ! ZEVS
         write (*,*) moment

!  BASSI:
         write(momentC, 'F6.4') moment
         lmomentC = index (momentC, ' ') -1

         linputPATH = index(inputPATH, ' ') -1
         File1=inputPATH(1:linputPATH)//"PMcrda"
     & //momentC(1:lmomentC)//".DAT"
         iFile1= index(File1, ' ') -1

         File2=inputPATH(1:linputPATH)//"PMcrs0a"
     & //momentC(1:lmomentC)//".DAT"
         iFile2= index(File2, ' ') -1

         write(*,*) File1
         write(*,*) File2

      end
      
