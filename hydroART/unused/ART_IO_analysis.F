c     
c     ---------------------------
      subroutine DumpParticles ()
c     ---------------------------
      include 'a_tree.h'

      open ( 45 , file = 'particles.dump' )
      do ip = 1 , np 
        write(45,'(i8,6(2x,g13.6))') ip, x(ip),y(ip),z(ip),
     &                                   vx(ip),vy(ip),vz(ip)
      enddo
      close ( 45 )
      
      return
      end
c
c     ----------------------------------------------------------
      subroutine Read_Particles_ASCII ( FileName1 , FileName2 )
c     ----------------------------------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c     the same as Read_Particles_Binary but reads files in ASCII format
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
      character*256 FileName1 , FileName2
c 
      ngrid = ng

      nfn1 = index ( FileName1 , ' ' ) - 1
      nfn2 = index ( FileName2 , ' ' ) - 1
c
      open ( 3 , file = FileName1(1:nfn1) , status = 'old' )

c.... read control information and check whether it has proper structure
c
c     set the header here 
c
      HEADER = 'GF RUN w/ Cooling and SF'
      read      (3,*) 
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                   ,Ocurv,Omb0,extras
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRIDC,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &               ' Omegab_0=',F7.3,' Hubble=',f7.3)
      if( nspecies .eq. 0 ) nspecies = 1
      if ( nspecies .gt. nspec ) then
        write (*,*)
     &      ' nspec in PARAMETER < nspecies in TAPE-FILE:'
         write (*,*) ' nspecies=',nspecies,' > ',nspec
         STOP
      endif
      If( npmax .lt. lspecies(nspecies) ) then
         write (*,*) ' Wrong number of particles: '
         write (*,*) ' should be =',lspecies(nspecies),' (lspecies)'
         write (*,*) ' but is set to ',npmax,' in a_setup.h...'
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo
         STOP
      Endif 

      nbyte  = nrecl * 4
      nacces = nbyte
#ifdef NBYTEWORD
      nacces = nrecl
#endif

      open ( 1 , file = FileName2(1:nfn2), 
     &	         status = 'old' )
 
      rewind 3


      N_particles = lspecies(nspecies)
      np = N_particles 

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      do i = 1 , nspecies
        write(iOWork,*) 'ispecie =',i,lspecies(i)
      enddo
      write(iOWork,*) 'Read_Particles_Binary: reading np=',np,
     &                ' particles'
      close(iOWork)
#endif

      if ( np .eq. 0 .or. np .gt. npmax ) then 
        call Open_ASCII_File ( iOErr , errorfile , append ) 
        write(iOErr,*) '1 : error in Read_Particles_Binary: bad np =',
     &             np
        write(iOErr,*) 'stopping...'
        close ( iOErr )
        stop
      endif

      Npages      = (N_particles -1)/npage + 1
      N_in_last   = N_particles - npage*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',nspecies
      write (*,*) ' N_in_last=',N_in_last

      do i = 1 , N_particles
        read(1,*) x(i), y(i), z(i), vx(i), vy(i), vz(i)
      enddo

      close (1)
      close (3)

      do i = 1 , nspecies
        wpar(i) = wspecies(i)
        lsp(i)  = lspecies(i)
      enddo
      do i = 1 , nspecies
        if(i .eq. 1 ) then
          nsp(i,1) = 1
        else
           nsp(i,1) = lsp(i-1) + 1
        endif 
        nsp(i,2) = lsp(i)
      enddo 

C$OMP PARALLEL DO 
      do ic1 = 1 , np 
        pt(ic1) = t   ! set t to the current time moment
      enddo
c
c.... compute each particle's weight
c
      n_spec = nspecies-1
c#ifdef STARFORM
c      n_spec = nspecies-1
c#endif
      do ic0 = 1 , n_spec 
        do ic1 = nsp(ic0,1) , nsp(ic0,2)
           pw(ic1) = wpar(ic0)
        enddo
      enddo 

#ifdef STARFORM
      if ( start ) then 
        if ( nspecies .gt. nspec ) then
          write (*,*)
     &        ' nspec in PARAMETER < nspecies+1 for starformation:'
          write (*,*) ' nspecies+1=',nspecies+1,' > nspec =',nspec
          STOP
        endif
c        nspecies = nspecies + 1
        nsp(nspecies,1) = lsp(nspecies-1) + 1
        nsp(nspecies,2) = lsp(nspecies-1)
        lsp(nspecies) = nsp(nspecies,2)
        lspecies(nspecies) = lsp(nspecies)
      endif
      write(*,*) 'in reading particles:',nsp(nspecies,1),nsp(nspecies,2)
#endif

      return
      end
c
c     -----------------------------------------------------
      subroutine Write_Stars_ASCII ( FileName, xc, yc, zc )
c     ----------------------------------------------------- 
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      
      character*256 FileName
      call Open_ASCII_File ( 33 , FileName , sequent ) 
      nstars = nsp(nspecies,2) - nsp(nspecies,1) + 1
      istar1 = nsp(nspecies,1) 
      vfact = v0 / aexpn 
      write(33,*) t, aexpn
      write(33,*) nstars
      write(33,*) xc,yc,zc
      do i = 0 , nstars-1 
        is = istar1 + i 
c        ts = t0 * (t - tbirth(i+1)) * aexpn**2 
        write(33,333) x(is),y(is),z(is),
     &    vx(is)*vfact,vy(is)*vfact,vz(is)*vfact,
     &   pw(is)*aM0, tbirth(i+1)
      enddo
 333  format(8(e12.6,1x))
      close ( 33 ) 
      
      return
      end
c
c     ------------------------------------------------------
      subroutine Write_Stars_ASCII2 ( FileName, xc, yc, zc )
c     ------------------------------------------------------ 
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      
      integer i, is, istar1, nstars
      real*8 amsti, amstnow, ts, zstarII, zstarIa

      character*256 FileName
      call Open_ASCII_File ( 33 , FileName , sequent ) 
      nstars = nsp(nspecies,2) - nsp(nspecies,1) + 1
      istar1 = nsp(nspecies,1) 
      vfact = v0 / aexpn 
      tnow = age ( t ) 
      write(33,*) aexpn
      write(33,*) nstars
      write(33,*) xc,yc,zc
      do i = 0 , nstars-1 
        is = istar1 + i 
        amsti   = pw0(i+1) * aM0
        amstnow = pw(is)   * aM0
        ts = tnow - age ( tbirth(i+1)*1.d0 ) 
        zstarII = pw(i+1)*zstII(i+1)
        zstarIa = pw(i+1)*zstIa(i+1)
        write(33,333) i, x(is),y(is),z(is),
     &    vx(is)*vfact,vy(is)*vfact,vz(is)*vfact,
     &    amsti, amstnow, ts, zstarII, zstarIa
      enddo
 333  format((i7,1x),11(e12.6,1x))
      close ( 33 ) 
      
      return
      end
c
c     -----------------------------------------------------
      subroutine Write_DM_ASCII ( FileName , xc , yc , zc )
c     ----------------------------------------------------- 
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
c.... write the first specie of DM
c      
      character*256 FileName
      call Open_ASCII_File ( 33 , FileName , sequent ) 
      nstars = nsp(1,2) - nsp(1,1) + 1
      istar1 = nsp(1,1) 
      vfact = v0 / aexpn 
      write(33,*) t, aexpn
      write(33,*) nstars
      write(33,*) xc,yc,zc
      do is = istar1 , nstars
        write(33,334) x(is),y(is),z(is),
     &    vx(is)*vfact,vy(is)*vfact,vz(is)*vfact,
     &   pw(is)*aM0
      enddo
 334  format(7(e12.6,1x))
      close ( 33 ) 
      
      return
      end
