c     ------------------------------------------------------------
c
c     HALO FINDER (fall 2002)
c 
c     Finds virial radius correctly. Removes unbound particles.
c
c     Escape velocity at a given halo centric radius is computed via
c     analytical expression for a NFW halo
c
c     Andrey Kravtsov 
c
c     Last Updated : 09/11/03 (Daisuke Nagai)
c       Added a routine Recompute_Halo_Mass()  - 03/09/03
c       HF is implemented on the HART code     - 09/11/03
c
c     Check hfind.inp a_hfind.h?
c     For isolated halo, you can comment out RemoveSmall() & Cannibalism()
c
c     ------------------------------------------------------------


c     ------------------------------
      subroutine HaloFind ( iWrite ) 
c     ------------------------------
c
c     Output Halo Catalog?    iWrite = 0/1/2/3 (No/DM/DM+Gas+Stars/Both)
c
c     aDeltaMin : minimum peak overdensity for halo selection
c     aDeltavir : virial overdensity w.r.t mean mass density of the universe
c     rfind     : search radius [kpc]  - exclude all particles within rfind
c     ismooth   : ismooth = 1 : use the smoothed density field
c                 ismooth = 0 : use the total density field 
c

      include 'a_hfind.h'
      character*40 FileName
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
CEVERINO052206      character*5 aname
      character*256 fname1, fname2, fname3, fname4, fname5, fname6
      character *256 datpath, respath
      integer ldpath, lrpath
      common / abin / alpha
      data alpha / 1.5 /      
      character*80 ctmp
      character*5 chaexp

      parameter( aDeltaMin = 5000., aDeltavir = 180., rfind = 100. )
      parameter( ismooth = 1 )

      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)
      chaexp = digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)

      DeltaMin = aDeltaMin
      Deltavir = aDeltavir
   
c      write(*,*) ' input parameters :'
c      write(*,*) aexpn, chaexp, DeltaMin, Deltavir, rfind, iIN

c
c.... read control data
c
      open ( 15 , file = 
     &     './hfind.inp', 
     &     status = 'old' )
      read  ( unit = 15 , fmt = '(10x,    a256)')  datpath
      read  ( unit = 15 , fmt = '(10x,    a256)')  respath
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  aboxh
      read  ( unit = 15 , fmt = '(70x,     i10)')  ibtype
      read  ( unit = 15 , fmt = '(70x,     i10)')  imark
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  rmin
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  rmax
      read  ( unit = 15 , fmt = '(70x,     i10)')  nbins
      read  ( unit = 15 , fmt = '(70x,     i10)')  miter
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  rvel
      read  ( unit = 15 , fmt = '(70x,     i10)')  nmin      
      read  ( unit = 15 , fmt = '(70x,     i10)')  nmax
      close ( 15 ) 
c
c.... write it out for control
c
      write  ( unit = * , fmt ='(/'' data path:'',/A)')  datpath
      write  ( unit = * , fmt ='(/'' results path:'',/A)') respath
      write  ( unit = * , fmt = '( 1x,''Lbox:   '', 1pe10.2)') aboxh
      write  ( unit = * , fmt = '( 1x,''ibtype:   '', i2,4x,
     &               ''imark:   '', i2,4x, 
     &               ''miter:  '', i2)')  ibtype, imark, miter
      write  ( unit = * , fmt = '( 1x,''rmin:   '', 1pe10.2,4x,
     &               ''rmax:   '', 1pe10.2,4x, 
     &               ''nbins:  '', i4)')  rmin, rmax, nbins
      write  ( unit = * , fmt = '( 1x,''rvel:   '', 1pe10.2,4x,
     &               ''nmin:  '', i4)') rvel, nmin

      ldpath = index( datpath , ' ' ) - 1 
      fname1 = datpath(1:ldpath)//'/PMcrd.DAT '
      fname2 = datpath(1:ldpath)//'/PMcrs0.DAT '
      fname3 = datpath(1:ldpath)//'/rho_smooth_a'//chaexp//'.dat '

      call InitArrays    ()
c      call Read_Particles_Binary_HF ( iIN, fname1, fname2, fname3) 
      call Initialize_HF ( ismooth ) 
      call Conversions   ( aboxh )

      ! [comoving 1/h kpc] -> grid units
      rhalo     = rfind * rkpc2g 
      rmin      = rmin * rkpc2g 
      rmax      = rmax * rkpc2g 
c      ibtype = 2
c      imark = 0
c      miter = 5
c      rvel = 1.0 
c      nmin = 5

      if ( ibtype .eq. 2 ) then 
        nbins = int(2.*(rmax/rmin-1.0)**(2./3.)-1.0) 
      endif
      if ( nbins .gt. nbmax ) then 
        write(*,*) '* error: nbins =',nbins,' > nbmax in header: ',nbmax
        write(*,*) '* increase nbmax or decrease nbins and rerun...'
        stop
      else
        write(*,*) 'will construct profiles with nbins=',nbins
      endif

      write(*,*) 'Constructing Linked List...'
      call LL_Construct_HF  ()
      
c      do i = 1 , np1
c        if ( (dnb2(ind2(i)) .gt. 100.)
c     &       .and.(dnb2(ind2(i)).lt.50000.)) then
c        write(*,*) i,dnb2(ind2(i))
c        endif
c      enddo

      write(*,*) 'Finding Halos...'
      call FindHaloes ( rhalo )
      call RemoveSmall   ( 5 )
c      write(*,*) 'Iterating halos...'
c      call IterateHaloCM ( rhalo , rminshift )      
      
c      write(*,*) 'Eating...'
c      call Cannibalism   ()
      call HL_Construct  ()

      write(*,*) 'Computing halo properties...'
      call Compute_Halo_Properties_DM ( ibtype, rmin, rmax, rvel, 
     &                               nbins, imark, miter )

c      write(*,*) 'Recomputing halo mass without binning profiles...'
c      call Recompute_Halo_Mass ( nmin, nmax )

      call Remove_Velocity_Duplicates ()
c      call Cannibalism200()
c      call ComputeVmax ()
      if ( iWrite .eq. 1 .or. iWrite .eq. 3 ) then
         lrpath = index( respath , ' ' ) - 1 
         fname1 = respath(1:lrpath)//'/hlist_a'//chaexp//'.dat '
         fname2 = respath(1:lrpath)//'/hpro_a'//chaexp//'.dat '
         fname3 = respath(1:lrpath)//'/hp_a'//chaexp//'.dat '
         fname4 = respath(1:lrpath)//'/hprs_a'//chaexp//'.dat '
         fname5 = respath(1:lrpath)//'/p_a'//chaexp//'.dat '
         fname6 = respath(1:lrpath)//'/hpn_a'//chaexp//'.dat '
         
         write(*,*) 'Writing halo catalogs...'
         call Write_Halo_Catalogs ( fname1, fname2 ,
     &        ibtype , rmin , rmax , nbins , nmin ,
     &        rvel, rhalo )
         write(*,*) 'Writing halo particles to the file...'
         call Write_Halo_Particles ( fname3, fname4, 
     &                              ibtype, rmin, rmax, rvel, 
     &                              nbins, imark, nmin )
         write(*,*) 'Converting hp_*** to hpn***'
         call Write_Halo_Particles_Convert ( fname3, fname6 )
         call Write_Density_HF ( ismooth )
         call Write_Particles_HF ( fname5 )
      endif

      
      if ( iWrite .eq. 2 .or. iWrite .eq. 3 ) then
         lrpath = index( respath , ' ' ) - 1 
         fname1 = respath(1:lrpath)//'/h_blist_a'//chaexp//'.dat '
         fname2 = respath(1:lrpath)//'/h_bpro_a'//chaexp//'.dat '
         fname3 = respath(1:lrpath)//'/h_bvpro_a'//chaexp//'.dat '

         ! want to use different rmin, rmax, nbins here?
c         Deltavir = 340.
c         rmin = 5.0 * rkpc2g     ! min. radius 
c         rmax = 2000.0 * rkpc2g  ! max. radius 
c         nbins = 45 
         write(*,*) 'Computing halo properties with DM, Gas and Stars..'
         call Compute_Halo_Properties ( fname1, fname2, fname3,
     &        rmin, rmax, nbins, nmin, 1 )
      endif

      return
      end
c
c     ------------------------------------
      subroutine Initialize_HF ( ismooth ) 
c     ------------------------------------
c
c     Prepare density fields for halo finding
c         ismooth = 0   total density
c         ismooth = 1   smoothed density 
c
      include 'a_hfind.h'
c      real    wspecies(nspec)
c      integer lspecies(nspec)
c      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))

c
c...  use only the first species
c
      np1 = nsp(1,2)            ! only the first specie
      write(*,*) ' Use only the first species for HF : np1 =',np1

c  
c...  correct DM particle weight by baryon fraction
c
      fb = Omb0 / Om0
      do ic1 = 1, np1
         pw2(ic1) = wpar(1) / (1.0-fb)
      enddo

c
c...  assign density to dnb2 (with baryon correction x 1/(1-fb) in pw(ic) )
c
      if ( ismooth .eq. 1 ) then
         write(*,*) ' smoothing density field...'
         nsmooth = 24  ! # of part. to use for computing smoothed density 
         call smoothden_fortran ( nsmooth, ng )
         write(*,*) ' after SmoothDen() routine in ART_HaloFind.F '
      else 
         write(*,*) ' assigning the local total density to ddummy...'
         write (*,*)'CEVERINO: the ismooth flag is down. It is always 1'
         STOP
c$$$         call Get_MaxLevelNow ()
c$$$         call Assign_Density ( MinLevel , MaxLevelNow, 0 )
c$$$         do Level = MaxLevelNow , MinLevel , -1
c$$$            if ( Level .le. MaxLevelNow ) then
c$$$c              write(*,*) 'Level = ',Level
c$$$               call Assign_Dummy_Level ( Level )
c$$$            endif
c$$$         enddo
c$$$         call Assign_Density_HF()
      endif

      write(*,*) 'Sorting particles...'
      call SortParticles ()

      return
      end
c
c     ------------------------------------------------------------
      subroutine Compute_Halo_Properties ( fname1, fname2, fname3,
     &                           rmin, rmax, nbins, nmin, ictype )
c     ------------------------------------------------------------
c
c     Output properties of all halos.
c
c     rmin, rmax  ! min/max radius for the profiling in code units
c
c     output : h_blist_aX.XXX.dat
c              h_bpro_aX.XXX.dat
c              h_bvpro_aX.XXX.dat
c
      include "a_hfind.h"
c
      character*256 fname1, fname2, fname3
      parameter ( Tcold = 2.e5) ! Kelvin

      ! grid units -> kpc (proper)
      rfact = r0 * 1.e3 * aexpn / hubble 

c
c.... profiles around centers of DM halos 
c
      nfn = index( fname1 , ' ' ) - 1 
      open ( 11 , file = fname1(1:nfn) ) 
      write(11,150) 
 150   format('# Halo list with properties determined from profiles')
      write(11,151) Deltavir, Tcold
 151  format('# overdensity (DM+bar.) defining the outer radius=',e9.3,
     &';  Tcold=',e9.3)
      write(11,152) rmin*rfact, rmax*rfact, nbins 
 152  format('# rmin,max=',e9.3,1x,e9.3,' [1/h kpc] (proper);'
     &  '  nrbin=',i4)
      write(11,156)
 156  format('#   x        y         z   min(rout,rmax)  rvir    Mg   
     & Mcg    M*     Mdm     <Zg>    <Z*>    <t*>    Mhalo ')
      write(11,158)
 158  format('# coord. in /h Mpc (comoving)      kpc (proper) 
     &all masses in /h Msun')
c
      nfn = index( fname2 , ' ' ) - 1 
      open ( 12 , file = fname2(1:nfn) )

      nfn = index( fname3 , ' ' ) - 1 
      open ( 13 , file = fname3(1:nfn) ) 
c
      do i = 1 , nhalo 
        if ( nhp(i) .gt. nmin ) then 
          write(*,*) 'processing halo ',i
          ihc = i ! remember the index of current halo
          call DensityCenter_Profiles 
     &      ( Deltavir, Tcold, xh(i), yh(i), zh(i), rh(i),
     &        rmin , rmax , nbins , ictype )
c          write(*,*) i, xh(i),yh(i),zh(i),rh(i),rhvir(i)
        endif
      enddo
      close ( 11 ) 
      close ( 12 ) 
      close ( 13 ) 

 220  return
      end
c
c     ---------------------------------------------------------
      subroutine DensityCenter_Profiles ( dout , Tcold , 
     &                            xc , yc , zc , rmass , 
     &                            rmin , rmax , nrbin , ictype )
c     ---------------------------------------------------------
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    dout : minimum density defining the outer radius of the cloud
c              Tcold - temperature (in K) to use to define cold gas 
c      integer: ictype - center type (0 - GC, 1 - halos) to determine what
c              to output   
c
c      count mass to dout or rmass whichever comes first 
c            rmass is in code units
c
c      itide = 1 if the halo is tidally truncated, 0 if not (rmass=rvir) 
c
#     include "a_def.h"
      include 'a_hfind.h'
c
      real xc , yc , zc , rmin , rmax 
      real rmass, rvir
      integer nrbin, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a, pfact
      real a2t, anow, tnow
      parameter ( npoints = 200000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real vcx(npoints), vcy(npoints), vcz(npoints), cw(npoints)
      real ptemp(npoints), pz(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pz
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dgbin(0:nmaxbins),dgcbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), dtibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real ndmp(0:nmaxbins),entbin(0:nmaxbins)
      real zbin(0:nmaxbins), zsbin(0:nmaxbins), zsbinIa(0:nmaxbins)
      real ppvx(0:nmaxbins), ppvy(0:nmaxbins), ppvz(0:nmaxbins)
      real vvrms(0:nmaxbins)
      real svx(0:nmaxbins), svy(0:nmaxbins), svz(0:nmaxbins)
      real svrms(0:nmaxbins)
      real gvx(0:nmaxbins), gvy(0:nmaxbins), gvz(0:nmaxbins)
      real gvrms(0:nmaxbins), gm(0:nmaxbins)
      real*8 amst(0:nmaxbins), amstot(0:nmaxbins), nst(0:nmaxbins)
      real*8 amgtot(0:nmaxbins), amgctot(0:nmaxbins)
      real*8 amdmtot(0:nmaxbins), tsbin(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(0:nmaxbins), ncbin(0:nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn
      real*8 rhofact, age
      real*8 zave, amgave, zsave, tave, amsave

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) 

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rhofact = rho0 * 1.d-18 / aexpn**3

      call Get_MaxLevelNow ()

      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin
c
c.... select random samples in the volume 
c
      do i = 1 , npoints
        rr    = rmax * RANDd ( mrand )
        theta = pi2 * RANDd ( mrand ) - pi
        phi   = pi2 * RANDd ( mrand )
        xpn(i) = rr * cos(theta) * cos(phi) + xc
        ypn(i) = rr * cos(theta) * sin(phi) + yc
        zpn(i) = rr * sin(theta) + zc
        if ( xpn(i) .lt. xmin ) xmin = xpn(i)
        if ( xpn(i) .gt. xmax ) xmax = xpn(i)
        if ( ypn(i) .lt. ymin ) ymin = ypn(i)
	if ( ypn(i) .gt. ymax ) ymax = ypn(i)
        if ( zpn(i) .lt. zmin ) zmin = zpn(i)
        if ( zpn(i) .gt. zmax ) zmax = zpn(i)
      enddo
c
      do i = 1 , npoints
        xx  = xpn(i)
        yy  = ypn(i)
        zz  = zpn(i)
        if ( xx .lt. 1.0 ) then 
          xx = xx + ng
        else 
          if ( xx .gt. narr ) then 
            xx = xx - ng
          endif
        endif
        if ( yy .lt. 1.0 ) then 
          yy = yy + ng
        else 
          if ( yy .gt. narr ) then 
            yy = yy - ng
          endif
        endif
        if ( zz .lt. 1.0 ) then 
          zz = zz + ng
        else 
          if ( zz .gt. narr ) then 
            zz = zz - ng
          endif
        endif
        
        idcell = iFindCell ( MaxLevel , xx , yy , zz )
c
c....   gas density
c
        pd(i) = hvar(1,idcell)
c
c....   pressure
c
        pp(i) = hvar(6,idcell)

        if ( hvar(1,idcell) .gt. 0. ) then
          T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
	  ptemp(i) = T_0 * T_code * a2i
c
c....   metallicity in units of solar 
c
          pz(i) = 50. * (hvar(izII,idcell)+hvar(izIa,idcell)) /
     &            hvar(1,idcell)
          rhor = 1.d0 / hvar(1,idcell)
          vcx(i) = hvar(3,idcell) * rhor 
          vcy(i) = hvar(4,idcell) * rhor 
          vcz(i) = hvar(5,idcell) * rhor 
          cw(i)  = hvar(1,idcell) * CellVol(iLv(idcell))
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo
c
c.... construct profiles
c
      do i = 0 , nmaxbins
        rbin(i) = 0.
        dgbin(i) = 0.
        dgcbin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        nbin(i) = 0
        ncbin(i) = 0
        ndmp(i) = 0.
        nst(i) = 0
        amst(i) = 0
        amstot(i) = 0
        amgtot(i) = 0.
        amgctot(i) = 0.
        amdmtot(i) = 0.
        zbin(i) = 0.
        zsbin(i) = 0.
        zsbinIa(i) = 0.
        tsbin(i) = 0.
        ppvx(i) = 0.
        ppvy(i) = 0.
        ppvz(i) = 0.
        vvrms(i) = 0.
        svx(i) = 0.
        svy(i) = 0.
        svz(i) = 0.
        svrms(i) = 0.
        gvx(i) = 0.
        gvy(i) = 0.
        gvz(i) = 0.
        gvrms(i) = 0.
        gm(i) = 0.
      enddo

      nstars = nsp(nspec,2) - nsp(nspec,1) + 1
      istar1 = nsp(nspec,1)

      ! Daiuke added the lines below (1/31/04)
#ifndef STARFORM
      istar1 = nsp(nspec-1,2) + 1
#endif
      write(*,*) nspec, nsp(nspec,2), nsp(nspec,1)
      write(*,*) ' from HF : ',nstars, nspec, istar1, np 

      do i = 1 , np 
        xx = x(i)
        yy = y(i)
        zz = z(i)
        diff_x = xc - xx
        diff_y = yc - yy
        diff_z = zc - zz
        corr_x = zero
        corr_y = zero
        corr_z = zero
        
        if ( abs(diff_x) .gt. nf67 ) then              
          if ( diff_x .gt. 0. ) then 
            corr_x = ng
          else
            corr_x = -ng
          endif
        endif 
        if ( abs(diff_y) .gt. nf67 ) then              
          if ( diff_y .gt. 0. ) then
            corr_y = ng
          else
            corr_y = -ng
          endif
        endif 
        if ( abs(diff_z) .gt. nf67 ) then              
          if ( diff_z .gt. 0. ) then
            corr_z = ng
          else
            corr_z = -ng
          endif
        endif 

c....  correct coordinate to insure periodic bundaries
        xx = xx + corr_x
        yy = yy + corr_y
        zz = zz + corr_z
        
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          if ( i .lt. istar1 ) then ! DM particle
            ddmbin(ibin) = ddmbin(ibin) + pw(i)
            ndmp(ibin) = ndmp(ibin) + 1
            ppvx(ibin) = ppvx(ibin) + vx(i)*pw(i) 
            ppvy(ibin) = ppvy(ibin) + vy(i)*pw(i)
            ppvz(ibin) = ppvz(ibin) + vz(i)*pw(i)
            vvrms(ibin) = vvrms(ibin) + 
     &                   (vx(i)**2 + vy(i)**2 + vz(i)**2)*pw(i)
          else ! stellar particle
            amst(ibin) = amst(ibin) + pw(i)
            nst(ibin) = nst(ibin) + 1
            zsbin(ibin) = zsbin(ibin) + pw(i)*zstII(i-istar1+1)
            zsbinIa(ibin) = zsbinIa(ibin) + pw(i)*zstIa(i-istar1+1)
            svx(ibin) = svx(ibin) + vx(i)*pw(i) 
            svy(ibin) = svy(ibin) + vy(i)*pw(i)
            svz(ibin) = svz(ibin) + vz(i)*pw(i)
            svrms(ibin) = svrms(ibin) + 
     &                   (vx(i)**2 + vy(i)**2 + vz(i)**2)*pw(i)
            tsbin(ibin) = tsbin(ibin) + tbirth(i-istar1+1)*pw(i)
          endif
        endif
      enddo

      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          if ( ptemp(ii) .le. Tcold ) then ! cold gas 
            ncbin(ibin) = ncbin(ibin) + 1
            dgcbin(ibin)  = dgcbin(ibin) + pd(ii)
          endif
          pbin(ibin) = pbin(ibin) + pp(ii)*cw(ii)
          zbin(ibin) = zbin(ibin) + pz(ii)*cw(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)*cw(ii) ! mass weighted
          gvx(ibin) = gvx(ibin) + vcx(ii)*cw(ii) 
          gvy(ibin) = gvy(ibin) + vcy(ii)*cw(ii)
          gvz(ibin) = gvz(ibin) + vcz(ii)*cw(ii)
          gvrms(ibin) = gvrms(ibin) + 
     &           (vcx(ii)**2 + vcy(ii)**2 + vcz(ii)**2)*cw(ii)
          gm(ibin) = gm(ibin) + cw(ii)
        endif
      enddo
c
c.... prepare profiles
c
      vfact = v0 / aexpn
      do i = 0 , nbins
        if ( nbin(i) .gt. 0 ) then 
          rnp = 1.d0 / gm(i)
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
          pbin(i) = pbin(i) * rnp 
          tbin(i) = tbin(i) * rnp
          zbin(i) = zbin(i) * rnp
          gvx(i) = gvx(i) * rnp 
          gvy(i) = gvy(i) * rnp 
          gvz(i) = gvz(i) * rnp 
          vmn = gvx(i)**2 + gvy(i)**2 + gvz(i)**2
          gvrms(i) = sqrt(abs(gvrms(i)*rnp - vmn)) * vfact
          gvx(i) = gvx(i) * vfact 
          gvy(i) = gvy(i) * vfact 
          gvz(i) = gvz(i) * vfact         
          dgcbin(i) = dgcbin(i) / nbin(i)
        endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.0
      dti = 0.0
      pfact = (gamma-1.0) * P0 / aexpn**5 / 1.38066d-16
      ivir = 0
      irflag = 0 

      dlout = log10(dout)
c
c.... initialize variables for averages 
c
      zave  = 0.0 ! average metallicity of gas 
      zsave = 0.0 ! average metallicity of stars 
      tave  = 0.0 ! average stellar age
      amgave = 0.0 ! sum of the bin gas masses 
      amsave = 0.0 ! sume of the stellar masses in bins

      amst_tot = 0.0 
      amg_tot = 0.0 
      amgc_tot = 0.0 ! cold gas mass
      volr1 = 0.

      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
        else
          rl = 10.**(rlmin + float(i-1)*drl)
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
	if ( volr .gt. 0. ) then
          dmi = dmi + ddmbin(i) 
          amst_tot = amst_tot + amst(i) 
          amstot(i) = amst_tot
          dmg = dgbin(i)*volr
          amg_tot = amg_tot + dmg 
          amgtot(i) = amg_tot
          amgc_tot = amgc_tot + dgcbin(i) * volr
          amgctot(i) = amgc_tot 
          amdmtot(i) = dmi
          dti = dti + ddmbin(i) + dgbin(i)*volr + amst(i)
          pbin(i) = pbin(i) * pfact 
c
c....     compute averages within rout
c
          if ( irflag .eq. 0 ) then ! if outer radius not reached yet
            zave   = zave + zbin(i) * dmg
            amgave = amgave + dmg 
            zsave  = zsave + (zsbin(i) + zsbinIa(i))
            tave   = tave + tsbin(i) 
            amsave = amsave + amst(i)
          endif
c
          ddmbin(i) = ddmbin(i)/volr 
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i) = dti / (4.18879 * rr**3)
c          dbi1 = dtibin(i-1)-ddmibin(i-1) ! baryons only
c          dbi2 = dtibin(i)-ddmibin(i)
          dbi1 = dtibin(i-1) ! total density 
          dbi2 = dtibin(i)

c
c...      compute the virial radius
c
          if ( (dbi2 .lt. dout .and.
     &          dbi1 .ge. dout) ) then 
             rrl = log10(rr)
             rll = log10(rl) 
             rri = 1.d0 / (rrl - rll)
             dlbi1 = log10(dbi1)
             dlbi2 = log10(dbi2)
             rvir = 10.d0**((dlout * (rrl - rll) + rll*dlbi2 - 
     &            rrl*dlbi1) / (dlbi2 - dlbi1))
          endif

          if ( (dbi2 .lt. dout .and.
     &          dbi1 .ge. dout) .and. irflag .eq. 0 ) then 
            ivir = i
            rrl = log10(rr)
            rll = log10(rl) 
            rri = 1.d0 / (rrl - rll)
            dlbi1 = log10(dbi1)
            dlbi2 = log10(dbi2)
            rout = 10.d0**((dlout * (rrl - rll) + rll*dlbi2 - 
     &             rrl*dlbi1) / (dlbi2 - dlbi1))
            rlout = log10(rout)
c
c....       gas mass
c
            if ( amgtot(i-1) .gt. 0 ) then 
              amgl1 = log10(amgtot(i-1))
            else
              amgl1 = -15.0
            endif
            if ( amgtot(i) .gt. 0 ) then 
              amgl2 = log10(amgtot(i))
            else
              amgl2 = -15.0
            endif
            ah = (amgl2 - amgl1) * rri
            bh = amgl1 - ah * rll 
            aM_gas = 10.d0**(ah*rlout + bh)
c
c....       cold gas mass
c
            if ( amgctot(i-1) .gt. 0 ) then 
              amgcl1 = log10(amgctot(i-1))
            else
              amgcl1 = -15.0
            endif
            if ( amgctot(i) .gt. 0 ) then 
              amgcl2 = log10(amgctot(i))
            else
              amgcl2 = -15.0
            endif
            ah = (amgcl2 - amgcl1) * rri
            bh = amgcl1 - ah * rll 
            aM_cgas = 10.d0**(ah*rlout + bh)
c
c....       stellar mass
c 
            if ( amstot(i-1) .gt. 0. ) then 
              amsl1 = log10(amstot(i-1))
            else
              amsl1 = -15.0 
            endif
            if ( amstot(i) .gt. 0. ) then 
              amsl2 = log10(amstot(i))
            else
              amsl2 = -15.0
            endif
            ah = (amsl2 - amsl1) * rri
            bh = amsl1 - ah * rll 
            aM_st = 10.d0**(ah*rlout + bh)
c
c....       DM mass
c
            if ( amdmtot(i-1) .gt. 0.0 ) then 
              amdl1 = log10(amdmtot(i-1))
            else
              amdl1 = -15.0 
            endif
            if ( amdmtot(i) .gt. 0.0 ) then 
              amdl2 = log10( amdmtot(i) )
            else
              amdl2 = -15.0 
            endif
            ah = (amdl2 - amdl1) * rri
            bh = amdl1 - ah * rll 
            aM_dm   = 10.d0**(ah*rlout + bh)
c            write(*,*)'Mdm:',amdmtot(i-1)*aM0*hubble,
c     &        amdmtot(i)*aM0*hubble,aM_dm*aM0*hubble

            aM_gas  = aM_gas * aM0 * hubble
            aM_cgas = aM_cgas * aM0 * hubble
            aM_st   = aM_st  * aM0 * hubble
            aM_dm   = aM_dm  * aM0 * hubble
            aMout   = aM_gas + aM_st 
            irflag  = 1 
          endif
c
c...  estimate mass if rmass is reached
c
          
          if ( (rr .gt. rmass .and.
     &          rl .le. rmass) .and. irflag .eq. 0 ) then 
            ivir = i
            rrl = log10(rr)
            rll = log10(rl) 
            rri = 1.d0 / (rrl - rll)
            rout = rmass
            rlout = log10(rout)
c
c....       gas mass
c
            if ( amgtot(i-1) .gt. 0 ) then 
              amgl1 = log10(amgtot(i-1))
            else
              amgl1 = -15.0
            endif
            if ( amgtot(i) .gt. 0 ) then 
              amgl2 = log10(amgtot(i))
            else
              amgl2 = -15.0
            endif
            ah = (amgl2 - amgl1) * rri
            bh = amgl1 - ah * rll 
            aM_gas = 10.d0**(ah*rlout + bh)
c
c....       cold gas mass
c
            if ( amgctot(i-1) .gt. 0. ) then 
              amgcl1 = log10(amgctot(i-1))
            else
              amgcl1 = -15.0
            endif
            if ( amgctot(i) .gt. 0. ) then 
              amgcl2 = log10(amgctot(i))
            else
              amgcl2 = -15.0
            endif
            ah = (amgcl2 - amgcl1) * rri
            bh = amgcl1 - ah * rll 
            aM_cgas = 10.d0**(ah*rlout + bh)
c
c....       stellar mass
c 
            if ( amstot(i-1) .gt. 0. ) then 
              amsl1 = log10(amstot(i-1))
            else
              amsl1 = -15.0 
            endif
            if ( amstot(i) .gt. 0. ) then 
              amsl2 = log10(amstot(i))
            else
              amsl2 = -15.0
            endif
            ah = (amsl2 - amsl1) * rri
            bh = amsl1 - ah * rll 
            aM_st = 10.d0**(ah*rlout + bh)
c
c....       DM mass
c
            if ( amdmtot(i-1) .gt. 0.0 ) then 
              amdl1 = log10(amdmtot(i-1))
            else
              amdl1 = -15.0 
            endif
            if ( amdmtot(i) .gt. 0.0 ) then 
              amdl2 = log10( amdmtot(i) )
            else
              amdl2 = -15.0 
            endif
            ah = (amdl2 - amdl1) * rri
            bh = amdl1 - ah * rll 
            aM_dm   = 10.d0**(ah*rlout + bh)

            aM_gas  = aM_gas * aM0 * hubble
            aM_cgas = aM_cgas * aM0 * hubble
            aM_st   = aM_st  * aM0 * hubble
            aM_dm   = aM_dm  * aM0 * hubble
            aMout   = aM_gas + aM_st 
            irflag  = 1 
          endif

         
c          if ( ibin .gt. 0 .and.
c     &         tbin(i) .gt. 1.e3 .and. tbin(i-1) .lt. 1.e3 ) then 
c             irflag = 1
c             rout = rr
c             ivir = i
c          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
        volr1 = volr 
      enddo
      if ( ivir .eq. 0 ) then 
        rout = rmax
        ivir = nbins
        aM_gas = amgtot(ivir) * aM0 * hubble
        aM_cgas = amgctot(ivir) * aM0 * hubble
        rr = 10.**(rlmin + float(ivir)*drl)
        aM_dm = ddmibin(ivir) * (4.18879 * rr**3) * aM0 * hubble
        aM_st = amstot(ivir) * aM0 * hubble 
        aMout = aM_gas + aM_st
        write(*,*) ' error : ',ivir, rmass, rout, rvir, irflag
        write(*,*) ' error : ',xc, yc, zc, rp, ibin
      endif
c
      xcd=(xc-1.0)*r0
      ycd=(yc-1.0)*r0
      zcd=(zc-1.0)*r0
      rfact = r0 * 1.e3 * aexpn / hubble
      vfact = v0 / aexpn
      rdout  = rout * rfact   ! min(rt,rvir) in proper kpc
      rvdout = rvir * rfact   ! rvir in proper kpc 
      
c
c.... averages
c
      if ( amgave .gt. 0.0 ) then 
        zave = zave / amgave ! already solar
      endif
      if ( amsave .gt. 0.0 ) then 
        zsave = 50.0 * zsave / amsave ! -> solar 
        td = tave / amsave
        tave = age ( 1.d0*td ) ! code units -> Gyr
      endif
      
      if ( ictype .eq. 1 ) then ! DM halos
        write(11,10) xcd, ycd, zcd, rdout, rvdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc)*pmmsun
        write(12,10) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc)*pmmsun
        write(13,10) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave, amh(ihc)*pmmsun
      endif
      if ( ictype .eq. 0 ) then ! GC
        write(11,11) xcd, ycd, zcd, rdout, rvdout, 
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave
        write(12,11) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm, 
     &               zave, zsave, tave
        write(13,11) xcd, ycd, zcd, rdout, rvdout,  
     &               aM_gas, aM_cgas, aM_st, aM_dm,
     &               zave, zsave, tave
      endif
 10   format (3(f9.5,1x),1x,2(f8.3,1x),1x,6(e9.3,1x),f6.3,1x,e9.3)
 11   format (3(f9.5,1x),1x,2(f8.3,1x),1x,6(e9.3,1x),f6.3)
c
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) 
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
        abtot = (amgtot(i)+ amstot(i))/(4.18879 * rr**3)
        if ( amst(i) .gt. 0.0 ) then ! correct for stellar metallicity
          rnp = 1.0 / amst(i)
          zbin(i) = (0.02*zbin(i)*dgbin(i)*volr) / volr /
     &               dgbin(i)*50.
          zsbin(i) = zsbin(i) * rnp * 50.
          zsbinIa(i) = zsbinIa(i) * rnp * 50.
          svx(i) = svx(i) * rnp
          svy(i) = svy(i) * rnp
          svz(i) = svz(i) * rnp
          vmn = svx(i)**2 + svy(i)**2 + svz(i)**2
          svrms(i) = sqrt(abs(svrms(i)*rnp - vmn)) * vfact
          svx(i) = svx(i) * vfact 
          svy(i) = svy(i) * vfact 
          svz(i) = svz(i) * vfact 
        endif
        if ( ddmbin(i) .gt. 0.0 ) then 
          dm = ddmbin(i) * volr 
          rnp = 1.0 / dm
          ppvx(i) = ppvx(i) * rnp
          ppvy(i) = ppvy(i) * rnp
          ppvz(i) = ppvz(i) * rnp
          vmn = ppvx(i)**2 + ppvy(i)**2 + ppvz(i)**2
          vvrms(i) = sqrt(abs(vvrms(i)*rnp - vmn)) * vfact
          ppvx(i) = ppvx(i) * vfact 
          ppvy(i) = ppvy(i) * vfact 
          ppvz(i) = ppvz(i) * vfact           
        endif
        aMdm = ddmibin(i) * (4.18879 * rr**3) * aM0 * hubble
        abtot = abtot*rhofact
        rhost = amst(i) / volr 
        hmassi = (amgtot(i)+amstot(i))*aM0 + aMdm/hubble
        rcirc = rr * r0/hubble * aexpn * 3.0856d24
        vcirc1 = 1.15206d8 * sqrt(hmassi/rcirc)
        hmassi = (amgtot(i)+amstot(i))*aM0
        vcirc2 = 1.15206d8 * sqrt(hmassi/rcirc)
        hmassi = (amstot(i))*aM0
        vcirc3 = 1.15206d8 * sqrt(hmassi/rcirc)
          write(12,18)
     &         rmid*rfact, rr*rfact, dtibin(i)*rhofact,
     &         dgbin(i)*rhofact,rhost*rhofact, ddmbin(i)*rhofact,
     &         amgtot(i)*aM0*hubble, 
     &         amgctot(i)*aM0*hubble, 
     &         amstot(i)*aM0*hubble,
     &         aMdm,tbin(i),zbin(i),zsbin(i),zsbinIa(i),
     &         int(ndmp(i)),int(nst(i))
          write(13,19)
     &         rmid*rfact, rr*rfact, 
     &         amgtot(i)*aM0*hubble,
     &         amstot(i)*aM0*hubble,
     &         aMdm, gvrms(i),svrms(i),vvrms(i), vcirc1, vcirc2, vcirc3,
     &         int(ndmp(i)),int(nst(i))
      enddo
 18   format(2(f8.3,1x),12(e9.3,1x),i6,1x,i6)
 19   format(2(f8.3,1x),1x,9(e9.3,1x),i6,1x,i6)
c
c
c.... now label cells falling within the radius of clump
c
      do i = 1 , npoints
        rr    = rout * RANDd ( mrand )
        theta = pi2 * RANDd ( mrand ) - pi
        phi   = pi2 * RANDd ( mrand )
        xpn(i) = rr * cos(theta) * cos(phi) + xc
        ypn(i) = rr * cos(theta) * sin(phi) + yc
        zpn(i) = rr * sin(theta) + zc
        if ( xpn(i) .lt. xmin ) xmin = xpn(i)
        if ( xpn(i) .gt. xmax ) xmax = xpn(i)
        if ( ypn(i) .lt. ymin ) ymin = ypn(i)
	if ( ypn(i) .gt. ymax ) ymax = ypn(i)
        if ( zpn(i) .lt. zmin ) zmin = zpn(i)
        if ( zpn(i) .gt. zmax ) zmax = zpn(i)
      enddo
c
      do i = 1 , npoints
        xx  = xpn(i)
        yy  = ypn(i)
        zz  = zpn(i)
        if ( xx .lt. 1.0 ) then 
          xx = xx + ng
        else 
          if ( xx .gt. narr ) then 
            xx = xx - ng
          endif
        endif
        if ( yy .lt. 1.0 ) then 
          yy = yy + ng
        else 
          if ( yy .gt. narr ) then 
            yy = yy - ng
          endif
        endif
        if ( zz .lt. 1.0 ) then 
          zz = zz + ng
        else 
          if ( zz .gt. narr ) then 
            zz = zz - ng
          endif
        endif
        
        idcell = iFindCell ( MaxLevel , xx , yy , zz )
        ind(idcell) = 1
      enddo

      return
      end
c
c     -----------------------------------------------------------------
      subroutine FindRBin ( rd, ibtype, rmind, rmaxd, drd, ibin ) 
c     -----------------------------------------------------------------
c
c     given radius rd, and bin type ibtype find the corresponding
c     bin number ibin
c
c
      integer ibtype, ibin
      real rd, rmind, rmaxd, drd
      common / abin / alpha
c
      if ( rd .le. rmind ) then 
        ibin = 0 
      else
        if ( ibtype .eq. 0 ) then 
          ibin = int((log10(rd)-rmind)/drd)+1
        elseif ( ibtype .eq. 1 ) then ! linear 
          ibin = int((rd-rmind)/drd) + 1
        elseif ( ibtype .eq. 2 ) then ! custom
          ibin = int(2.*(rd/rmind)**(1./alpha)-1)
        else
          write(*,*) '* error in FindRBin unknown bin type =',ibtype
          return
        endif
      endif

      return
      end
c
c     -----------------------------------------------------------------
      subroutine FindIBin (ibin, ibtype, rmind, rmaxd, 
     &                     drd, rd1, rm, rd2, volr, voli  ) 
c     -----------------------------------------------------------------
c     given a bin number and binning info, return left, middle, and
c     right edges of the bin: rd1, rm, rd2, as well as differential volume
c     of the shell corresponding to the radial bin and integral volume
c     voli = V(<rr)
c
      include 'a_hfind.h'
      common / abin / alpha

      IF ( ibin .eq. 0 ) then 
        if ( ibtype .eq. 0 ) then 
          rd2 = 10.0**(rmind)
        else
          rd2 = rmind
        endif
        rd1 = 0
        rm  = rd2*0.5
        voli = pi43 * rd2**3 * rhoaver
        volr = voli 
      ELSE 
        if ( ibtype .eq. 0 ) then ! log10 bins
          rd1 = 10.0**(rmind + float(ibin-1)*drd)
          rm  = rd1 * 10.0**(drd*0.5)
          rd2 = 10.0**(rmind + float(ibin)*drd)
c          write(*,*) ibin, ibtype, rmind,rmaxd,drd,rd1,rm,rd2
c          pause
        elseif ( ibtype .eq. 1 ) then 
          rd1 = rmind + float(ibin-1)*drd
          rm  = rd1 + 0.5*drd 
          rd2 = rmind + float(ibin)*drd
        elseif ( ibtype .eq. 2 ) then ! custom 
          rd1 = rmind * ((ibin-1.0)/2. + 1)**alpha 
          rm  = rmind * ((ibin-0.5)/2. + 1)**alpha
          rd2 = rmind * (float(ibin)/2. + 1)**alpha
        else
          write(*,*) '* error in FindIBin unknown bin type =',ibtype
          return
        endif
        voli = pi43 * rd2**3 * rhoaver
        volr = pi43 * (rd2**3 - rd1**3) * rhoaver
      ENDIF
c
      return
      end

c     ------------------------------------------------------------------
      subroutine Compute_Halo_Properties_DM( ibtype, rpmin, rpmax, rvel, 
     &                                    npbin,imark, miter )
c     ------------------------------------------------------------------
c
c     ibtype = 0 - log10 bins, 1 - even bins
c     rpmin,max - min/max. radius for profile construction in grid units
c     imark = 1 - mark particles as they assigned to a halo, so that they
c                 they are not included in subsequent halo profiles
c                (this will prevent detection of most of the subhalos but
c                 will greatly speed-up finding of isolated halos)
c     imark = 0 - don't
c     rvel = radius within which to measure mean halo velocity 
c     miter = number of iterations for unbound particles 
c

      include 'a_hfind.h'
      real toohot 
      common / TOOHOT1 / toohot
c      parameter ( toohot = 2. )
      parameter ( dstopslope = -0.5 ) ! stop if density profile becomes flatter than this slope
      integer ibtype, npbin, miter
      dimension rr(1000) , vr(1000) , er(1000), dir(1000)
      real di(0:nbmax), dd(0:nbmax), ddi(0:nbmax)
      real*8 hmassi 

      toohot = 2.0

      if     ( ibtype .eq. 0 ) then 
        rmax   = rpmax
        rmin   = rpmin
        rlmax  = log10(rmax)
        rlmin  = log10(rmin)
        drl    = (rlmax - rlmin)/float(npbin)
        nbins  = int((rlmax - rlmin)/drl) + 1
        rmind  = rlmin
        rmaxd  = rlmax
        drd = drl 
      elseif ( ibtype .eq. 1 ) then 
        rmax   = rpmax
        rmin   = rpmin
        drbin  = (rmax - rmin)/float(npbin) ! bin width in grid units
        nbins = npbin
        rmind  = rmin
        rmaxd  = rmax
        drd = drbin
      elseif ( ibtype .eq. 2 ) then 
        rmax  = rpmax
        rmin  = rpmin 
        rmaxd = rpmax
        rmind = rpmin
        nbins = npbin
        drd = (rpmax - rpmin)/float(npbin)
      else
        write(*,*) 'Compute_Halo_Properties: unknown bin type:',ibtype
        write(*,*) 'exiting...'
        return
      endif
      
      if ( ibtype .eq. 0 ) then 
      write(*,*) 'ibtype =',ibtype
      write(*,*) 'rmin =',rmind,' rmax=',rmaxd,' dr =',drd
      do ic2 = 0 , nbins
        call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                  rd1, rm, rd2, volr, voli ) 
        write(*,87) ic2, rd1*rg2kpc, rm*rg2kpc, rd2*rg2kpc, volr, voli
      enddo
      endif

      if ( ibtype .eq. 1 ) then   
      write(*,*) 'ibtype =',ibtype
      write(*,*) 'rmin =',rmind,' rmax=',rmaxd,' dr =',drd
      do ic2 = 0 , nbins
        call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                  rd1, rm, rd2, volr, voli ) 
        write(*,87) ic2, rd1*rg2kpc, rm*rg2kpc, rd2*rg2kpc, volr, voli
      enddo
 87   format(i4,1x,3(f8.3,1x),1x,2(e11.6,1x)) 
      endif
      if ( ibtype .eq. 2 ) then 
c      write(*,*) 'ibtype =',ibtype
c      write(*,*) 'rmin =',rmind,' rmax=',rmaxd,' dr =',drd
      do ic2 = 0 , nbins
        call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                  rd1, rm, rd2, volr, voli ) 
        call FindRBin ( rm , ibtype, rmind, rmaxd, drd,
     &                  ibin ) 
c        write(*,88) ic2, ibin, 
c     &     rd1*rg2kpc, rm*rg2kpc, rd2*rg2kpc, volr, voli
      enddo
      endif
 88   format(i4,1x,i5,3(f9.3,1x),1x,2(e11.6,1x))     
c      pause

      dvirlog = log10(Deltavir)

      if ( nbins .gt. nbmax ) then 
        write(*,*) 'error : Compute_Halo_Properties : nbins > nbmax'
        write(*,*) 'nbins =',nbins,' nbmax =',nbmax
        stop
      endif

      rhalo  = rmax
      rhalo2 = rmax * rmax

      do ic1 = 1 , np1max 
        iSp(ic1) = nil 
      enddo

      nr = 5 ! # of points to use for linear regression
c      nr2 = max(nr/2,1)
      nr2 = 1
      if ( nr .lt. 2 ) then 
        write(*,*) 'number of bins for linear regression too small!!!'
        write(*,*) 'nr = nbins/4 =',nr
        write(*,*) 'stopping...'
        stop
      endif

      open ( 20 , file ='test.dat' )

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iter,nhpold,ic2,ic3,di,dd,ddi,ptot,ivir,xhalo,yhalo,zhalo)
C$OMP+PRIVATE(vxhalo,vyhalo,vzhalo,vxhnew,vyhnew,vzhnew,vmassnew)
C$OMP+PRIVATE(wdummy,ndummy,rhalo,rdvel,xd,yd,zd)
C$OMP+PRIVATE(imin,jmin,kmin,imax,jmax,kmax)
C$OMP+PRIVATE(i,ic,j,jc,k,kc,idm)
C$OMP+PRIVATE(vxd,vyd,vzd,diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE(rd,dv,vesc,resc,ibin,irmaxflag)
C$OMP+PRIVATE(rnp,vmn,rd1,rd2,rm,volr,voli,rld1,rld2,dil1,dil2,rvir)
C$OMP+PRIVATE(rhdum,ntot,ph,rvmax,vmax,imaxflag,irsflag,nrg,er,rr,dir)
C$OMP+PRIVATE(vr,a,ae,b,be,csq,q,hmassi,rcirc,rcircd,pdum,ptotd)
C$OMP+PRIVATE(ph1,ph2,irflag,pl1,pl2,ah,bh)
C$OMP+PRIVATE(rdkpc1,rdkpc2,od,odi)
      do ic1 = 1 , nhalo 
        
c        write(*,*) 'processing halo',ic1,'...',nhalo-ic1,' left'
c
c....   iterate to remove unbound particles 
c
        do iter = 0 , miter
c          write(*,*) 'iter=',iter,' rm=',rhmax(ic1),' vm=',vhmax(ic1),
c     &     ' rh=',rh(ic1)*rg2pkpc,' rvir=',rhvir(ic1)*rg2pkpc
          do ic2 = nil , nbins
            na   (ic2,ic1) = nil 
            pn   (ic2,ic1) = zero
            pnt  (ic2,ic1) = zero
            pvx  (ic2,ic1) = zero
            pvy  (ic2,ic1) = zero
            pvz  (ic2,ic1) = zero
            vcirc(ic2,ic1) = zero
            vrms (ic2,ic1) = zero
            di(ic2) = zero
            dd(ic2) = zero 
          enddo
          ptot = 0.
          ivir = 0 
          xhalo = xh(ic1) 
          yhalo = yh(ic1) 
          zhalo = zh(ic1) 
          
          vxhalo = vxh(ic1)
          vyhalo = vyh(ic1)
          vzhalo = vzh(ic1)
          vxhnew = 0.
          vyhnew = 0.
          vzhnew = 0.
          vmassnew = 0.0 
          wdummy = 0.0
          ndummy = nil                           ! # of halo particles
c
c....     define the radius for vhalo estimate
c
          if ( iter .eq. 0 ) then 
            rhalo = rmax
            rdvel = rh(ic1)/5.0
          else
            rhalo = min(4.*rh(ic1),rhvir(ic1))
            rdvel = rhmax(ic1) * rvel 
          endif
             
c....     convert from original grid to chaining mesh          
          xd    = xhalo - rhalo
          yd    = yhalo - rhalo
          zd    = zhalo - rhalo
          imin  = int(float(nll) * (xd - 1) / ng + 1)
          jmin  = int(float(nll) * (yd - 1) / ng + 1)
          kmin  = int(float(nll) * (zd - 1) / ng + 1)
          xd    = xhalo + rhalo
          yd    = yhalo + rhalo
          zd    = zhalo + rhalo
          imax  = int(float(nll) * (xd - 1) / ng + 1)
          jmax  = int(float(nll) * (yd - 1) / ng + 1)
          kmax  = int(float(nll) * (zd - 1) / ng + 1)
c....     sweep over neighbors checking periodic boundary conditions
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idm = iCL2(ic,jc,kc)                 ! read LL head
                do while ( idm .ne. nil )
                  if ( iSp(idm) .eq. nil ) then 
                  xd = x(idm)
                  yd = y(idm)
                  zd = z(idm)
                  vxd = vx(idm)
                  vyd = vy(idm)
                  vzd = vz(idm)
c....             enforce periodic boundary conditions
                  diff_x = xd - xhalo                  
                  diff_y = yd - yhalo 
                  diff_z = zd - zhalo 
                  corr_x = zero 
                  corr_y = zero 
                  corr_z = zero 

                  if ( abs(diff_x) .gt. ng/2 ) then 
                    if ( diff_x .gt. zero ) then 
                      corr_x = -ng 
                    else
                      corr_x =  ng
                    endif
                  endif
                  if ( abs(diff_y) .gt. ng/2 ) then 
                    if ( diff_y .gt. zero ) then 
                      corr_y = -ng 
                    else
                      corr_y =  ng
                    endif
                  endif
                  if ( abs(diff_z) .gt. ng/2 ) then 
                    if ( diff_z .gt. zero ) then 
                      corr_z = -ng 
                    else
                      corr_z =  ng
                    endif
                  endif

c....             compute distance particle - halo center squared
c                 corr_x , corr_y , corr_z - take care of periodicity
                  rd  = sqrt( (xd - xhalo + corr_x)**2 + 
     &                        (yd - yhalo + corr_y)**2 + 
     &                        (zd - zhalo + corr_z)**2   )   
                  rd = max(rd,1.e-7)
c                  write(*,*) 'rd=',rd,xd,yd,zd,xhalo,yhalo,zhalo
                  if ( iter .eq. 0 ) then 
                    dv  = 0.       
                    vesc = 1.e10
                  else
                    dv = sqrt( (vxd - vxhalo)**2 +
     &                         (vyd - vyhalo)**2 + 
     &                         (vzd - vzhalo)**2   ) * vg2kms                    
                    resc = rd * rg2pkpc
                    vesc = ve ( resc , rhmax(ic1) , vhmax(ic1) )
                    vesc = vesc * toohot * 1.15**(miter-iter)
                  endif
c
c....             count particle if it falls inside halo radius and v < vesc 
c                   
                  if ( rd .lt. rmax .and. dv .lt. vesc ) then             
                    call FindRBin(rd, ibtype, rmind, rmaxd, drd, ibin )          
                    ndummy = ndummy + 1
                    wdummy = wdummy + pw2(idm) 
                    if ( rd .lt. min(rdvel,rh(ic1)) ) then 
                      vxhnew = vxhnew + pw2(idm)*vx(idm)
                      vyhnew = vyhnew + pw2(idm)*vy(idm)
                      vzhnew = vzhnew + pw2(idm)*vz(idm)
                      vmassnew = vmassnew + pw2(idm)
                    endif
                    na (ibin,ic1) =  na(ibin,ic1) + 1
                    pn (ibin,ic1) =  pn(ibin,ic1) + pw2(idm) 
                    pvx(ibin,ic1) = pvx(ibin,ic1) + vx(idm)*pw2(idm) 
                    pvy(ibin,ic1) = pvy(ibin,ic1) + vy(idm)*pw2(idm)
                    pvz(ibin,ic1) = pvz(ibin,ic1) + vz(idm)*pw2(idm)
                    vrms(ibin,ic1) = vrms(ibin,ic1) + 
     &                                 (vx(idm)**2 + 
     &                                  vy(idm)**2 + 
     &                                  vz(idm)**2)*pw2(idm)
                    if ( imark .eq. 1 ) then 
                      if ( iter .eq. miter .and. 
     &                   rd .lt. rh(ic1) ) then 
                        iSp(idm) = 1
                      endif
                    endif
                  endif ! ( rd .le. rhalo2 )
                  endif
                  idm = iLL2(idm) ! next particle from linked list
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! end i
c
c....     compute halo velocity 
c          
          if ( vmassnew .gt. 0.0 ) then 
            vxh(ic1) = vxhnew / vmassnew 
            vyh(ic1) = vyhnew / vmassnew 
            vzh(ic1) = vzhnew / vmassnew 
          else
            vxh(ic1) = 0.
            vyh(ic1) = 0.
            vzh(ic1) = 0.
          endif
c
c....     compute velocity components, and rms of radial shells 
c
          ptot = 0.
          irmaxflag = 0 
          if ( iter .eq. 0 ) then 
            rh(ic1) = rmax          
            rhvir(ic1) = rmax
          endif

          do ic2 = 0 , nbins 
            if ( na(ic2,ic1) .gt. 0 ) then 
              rnp = 1.0/pn(ic2,ic1)
              ptot = ptot + pn(ic2,ic1)
              pvx(ic2,ic1) = pvx(ic2,ic1) * rnp
              pvy(ic2,ic1) = pvy(ic2,ic1) * rnp
              pvz(ic2,ic1) = pvz(ic2,ic1) * rnp         
              vmn = pvx(ic2,ic1)**2 + pvy(ic2,ic1)**2 + pvz(ic2,ic1)**2
              vrms(ic2,ic1) = sqrt(abs(vrms(ic2,ic1)*rnp - vmn))
            endif
            pnt(ic2,ic1) = ptot
            call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd2, volr, voli ) 
            di(ic2) = ptot / voli            
            dd(ic2) = pn(ic2,ic1) / volr
c
c....       estimate the virial radius 
c                      
            if ( ic2 .gt. 0 .and. ivir .eq. 0 .and. 
     &         di(ic2) .lt. Deltavir .and.
     &         di(ic2-1) .ge. Deltavir .and. iter .eq. 0 ) then 
              ivir = ic2
              rld1 = log10(rd1)
              rld2 = log10(rd2)
              if ( di(ic2-1) .gt. 0.0 ) then 
                dil1 = log10(di(ic2-1))
              else
                dil1 = -3.0 
              endif
              if ( di(ic2) .gt. 0.0 ) then 
                dil2 = log10(di(ic2))
              else
                dil2 = -3.0 
              endif
              rvir = 10.0**(
     &                      (dvirlog * (rld2 - rld1) + rld1*dil2 - 
     &                       rld2*dil1) / (dil2 - dil1)
     &                     ) ! interpolate linearly in log r - log delta
c              if ( ic1 .eq. 2467 ) then 
c                write(*,*) 'op-pa:', rld1,rld2,dil1,dil2,dvirlog,rvir
c                write(*,*) rd1*rg2kpc,rd2*rg2kpc,10.**dil1,10.**dil2,
c     &                      rvir*rg2kpc
c                pause
c              endif
              rhvir(ic1) = rvir 
              amhvir(ic1) = pi43 * rvir**3 * rhoaver * Deltavir
            endif
            if ( iter .eq. 0 ) then 
              ddi(ic2) = di(ic2)
              rh(ic1) = rhvir(ic1)
            elseif (       (ddi(ic2) .gt. 0.0)
     &               .and. (ddi(ic2)-di(ic2)) .gt. di(ic2) 
     &               .and. irmaxflag .eq. 0 ) then 
              irmaxflag = 1 
              rh(ic1) = rd2
            endif
c          write(*,*) rd2, ddi(ic2), di(ic2), rh(ic1)
          enddo
c
c....     compute halo's extent
c
          rh(ic1) = min(rh(ic1),rhvir(ic1))
          if ( iter .eq. 0 ) then 
            rhdum = rhvir(ic1)
          else
            rhdum = rh(ic1)
          endif

c
c....     compute components of halo's velocity 
c
          vxhalo = 0.
          vyhalo = 0.
          vzhalo = 0.
          ntot = 0 
          ph = 0.
          rd = zero
          do ic2 = 0 , nbins
            call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd, volr, voli ) 
c            if ( rhdum .le. rd ) then 
            if ( rd .le. rhdum ) then
              ph = ph + pn(ic2,ic1)
            endif
          enddo
          if ( ph .gt. 0 ) then 
            nhp(ic1) = int(ph/pw2(1))
            amh(ic1) = ph
          else
            nhp(ic1) = 0
            amh(ic1) = 0.
            amhvir(ic1) = 0.
c            goto 20
          endif
c
c....     construct circular velocity profile 
c 
          rd   = zero
          vmax = -1000.
          imaxflag = 0 
          irmaxflag = 0 
          irsflag = 0
          if ( iter .eq. 0 ) then 
            rsh(ic1) = rh(ic1)
          else
            rsh(ic1) = rhmax(ic1)/2.2
          endif

          do ic2 = 0 , nbins 
            call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd, volr, voli ) 
            ptot = pnt(ic2,ic1)
            ptotd = ptot - pn(ic2,ic1)
            
            if ( irmaxflag .eq. 0 .and. rd .le. rh(ic1) ) then 
              if ( ic2 .lt. nbins-nr ) then 
                nrg = 0
                do ic3 = 0 , nr-1
                  nrg = nrg + 1
                  call FindIBin ( ic2+ic3 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rr(nrg), volr, voli ) 
                  er(nrg) = 1.0
                  rr(nrg) = log10(rr(nrg) * rg2pkpc)
                  if ( dd(ic2+ic3) .gt. 0. ) then 
                    dir(nrg) = log10(dd(ic2+ic3))
                  else
                    dir(nrg) = -3.
                  endif
                enddo
                call fit ( rr, dir, nrg+1, er, 0, a, b, ae, be, csq, q )
c                if ( ic1 .eq. 48 ) then 
c                  do i = 1 , nr
c                    write(*,1000) rd*rg2kpc, rr(i), dir(i), b, be, 
c     &                        rh(ic1)*rg2kpc
c 1000               format(6(e9.3,1x))
c                  enddo
c                endif
                if ( b-be .ge. dstopslope ) then
c                  write(*,*) 10.**rr(1),10.**rr(2),10.**rr(3),
c     &                 10.** rr(4),10.**rr(5)
c                  write(*,*) dd(1),dd(2),dd(3),dd(4),dd(5)
c                  write(*,*) rd*rg2kpc,rhvir(ic1)*rg2kpc,b,be,b-be
c                  write(*,*)
                  rh(ic1) = rd
                  irmaxflag = 1
                endif
                if ( irsflag .eq. 0 ) then 
                  if ( b-be .le. -2.0 ) then 
                    rsh(ic1) = 10.0**(rr(nr2))
                    irsflag = 1
                  endif
                endif  
              endif 
            endif
c
c....  now work on vmax
c
            hmassi     = ptot * pmmsun
            rcirc      = rd * rg2pkpc * akpc
            vcirc(ic2,ic1) = 1.d-5 * 
     &           dsqrt(grav_c * hmassi * sun_mass/rcirc)
            if ( imaxflag .eq. 0 ) then 
              if ( ic2 .lt. nbins-nr ) then 
                nrg = 0
                do ic3 = 0 , nr-1
                  nrg = nrg + 1
                  call FindIBin ( ic2+ic3 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rr(nrg), volr, voli ) 
                  er(nrg) = 1.0
                  ptotd   = ptotd + pn(ic2+ic3,ic1)
                  hmassi  = ptotd * pmmsun
                  rr(nrg) = rr(nrg) * rg2pkpc ! r -> (physical) kpc
                  rcirc   = rr(nrg) * akpc ! r -> cm            
                  vr(nrg) = 1.d-5 * 
     &                 dsqrt(grav_c * hmassi * sun_mass/rcirc) 
                  if ( vr(nrg) .gt. vmax ) then
                    vmax = vr(nrg)
                    rvmax = rr(nrg) 
                  endif
                enddo
                call fit ( rr, vr, nr, er, 0, a, b, ae, be, csq, q )
                if ( rd .le. rh(ic1) ) then
                  rhmax(ic1) = rvmax
                  vhmax(ic1) = vmax
                  if ( b-be .lt. 0 ) then 
                    imaxflag = 1
                  endif
                else                  
                  rcirc = min(2.16*rsh(ic1),rh(ic1)*rg2pkpc)
                  rcircd = rcirc/rg2pkpc
                  ic3 = 0 
                  pdum = 0.0 
                  do while ( pdum .eq. 0.0 .and. ic3 .le. nbins ) 
                    call FindIBin ( ic3 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd2, volr, voli ) 
                    if ( rcircd .gt. rd1 .and. rcircd .le. rd2 ) then
                      pdum = pnt(ic3,ic1)
                      rcirc = rd2 * rg2pkpc
                    endif                        
c                    if ( ic1 .eq. 5 ) then 
c                      write(*,'(i6,2x,6(e9.3,1x))') 
c     &                   ic3,rcirc,rcircd,rd1,rd2,pdum,pnt(ic3,ic1)
c                    endif
                    ic3 = ic3 + 1
                  enddo        
c                  if ( ic1 .eq. 5 ) then 
c                    pause
c                  endif
                  rhmax(ic1) = rcirc
                  rcirc = rcirc * akpc
                  hmassi = pdum*pmmsun
                  vhmax(ic1) = 1.d-5 * 
     &                 dsqrt(grav_c * hmassi * sun_mass/rcirc)
                  imaxflag = 1                   
                endif
              endif 
            endif
c            if ( ic1 .eq. 1 ) then 
c                do i = 1 , nr
c                  write(*,999)ic2, rd*rg2pkpc, rr(i), vr(i), b, be, 
c     &                rvmax, vmax,
c     &                rh(ic1)*rg2kpc,rhmax(ic1),2.16*rsh(ic1),
c     &                vhmax(ic1)
c 999              format(i5,1x,11(e9.3,1x))
c                enddo
c                pause
c            endif
c            if ( ic2 .eq. nbins .and. imaxflag .eq. 0 ) then 
c              rhmax(ic1) = rcirc  
c              vhmax(ic1) = vcirc(ic2,ic1) 
c            endif
          enddo ! ic2


c          write(*,*) ic1, rh(ic1)*rg2kpc, rhvir(ic1)*rg2kpc

c          if ( ic1 .eq. 48 ) then 
c            pause
c          endif
c
c....     compute halo radius, mass and the # of particles it contains
c
          rhdum = min(rh(ic1),rhvir(ic1))         
          if ( iter .gt. 0 ) then 
            rd = 0.
            ph1 = 0.      
            ph2 = 0.
            irflag = 0 
            do ic2 = 0 , nbins
              ph2 = ph2 + pn(ic2,ic1)
              call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                        rd1, rm, rd2, volr, voli ) 
              if ( rd1 .le. rhdum .and. rd2 .ge. rhdum ) then 
                if ( ic2 .eq. 0 ) then 
                  ph = pn(ic2,ic1)
                else
                  rld1 = log10(rd1)
                  rld2 = log10(rd2)
                  if ( ph1 .gt. 0.0 ) then 
                    pl1 = log10(ph1)
                  else
                    pl1 = -15.0 
                  endif
                  if ( ph2 .gt. 0.0 ) then 
                    pl2 = log10(ph2)
                  else
                    pl2 = -15.0
                  endif
                  ah = (pl2 - pl1)/(rld2 - rld1) 
                  bh = pl1 - ah * rld1
                  ph = 10.0**(ah*log10(rhdum) + bh) ! interpolate linearly in log r - log M
                endif
              endif            
              ph1 = ph2 
            enddo
            nhp(ic1) = int(ph/pw2(1))
            amh(ic1) = ph
          endif
c
          if ( nhp(ic1) .eq. 0 ) goto 20

c            write(*,96)ic1,iter,xhalo*rg2Mpc,yhalo*rg2Mpc,zhalo*rg2Mpc,
c     &              vxh(ic1)*vg2kms, vyh(ic1)*vg2kms, vzh(ic1)*vg2kms,
c     &                vhmax(ic1), rhmax(ic1), rsh(ic1),rh(ic1)*rg2kpc, 
c     &              rhvir(ic1)*rg2kpc, nhp(ic1)

          if ( iter .gt. -1000 ) then 

            write(20,96)ic1,iter,xhalo*rg2Mpc,yhalo*rg2Mpc,zhalo*rg2Mpc,
     &              vxh(ic1)*vg2kms, vyh(ic1)*vg2kms, vzh(ic1)*vg2kms,
     &              vhmax(ic1), rhmax(ic1), rsh(ic1), rh(ic1)*rg2kpc, 
     &              rhvir(ic1)*rg2kpc, nhp(ic1)
 96         format(i5,1x,i1,1x,3(1x,f8.4),1x,3(1x,f8.2),2x,f10.2,2x,
     &             4(f7.2,1x),1x,i8)

            do ic2 = 0 , nbins
              call FindIBin ( ic2, ibtype, rmind, rmaxd, drd,
     &                        rd1, rm, rd2, volr, voli ) 
              rdkpc1 = rd1 * rg2kpc
              rdkpc2 = rd2 * rg2kpc
              od  = pn(ic2,ic1)  / volr
              odi = pnt(ic2,ic1) / voli 
              write(20,97) rdkpc1, rdkpc2, odi, ddi(ic2), od, 
     &          vcirc(ic2,ic1),
     &          vrms(ic2,ic1)*vg2kms, ve(rdkpc2,rhmax(ic1),vhmax(ic1)),
     &          pnt(ic2,ic1)*pmmsun, pnt(ic2,ic1)
             enddo

 97          format(2(2x,g11.6),1x,3(2x,g13.6),
     &              1x,3(2x,g11.6),3x,g13.6,3x,g13.6)
           endif

        enddo ! iter

 20     continue

      enddo ! end ic1

      close ( 20 )

c      open ( 21 , file = 'pisp.dat' , form='unformatted' )
c      write(21) np1 
c      write(21) (iSp(ic1),ic1=1,np1)      
c      close ( 21 ) 

      return
      end

c     --------------------------------------------
      subroutine Recompute_Halo_Mass( nmin, nmax )
c     --------------------------------------------
c
c     rvel = radius within which to measure mean halo velocity 
c     nmin = recompute mass for halos with nhp > nmin 
c     nmax = recompute mass for halos with nhp < nmax
c
      include 'a_hfind.h'
      parameter ( nbin = 1000000 )
      real ra(nbin),weight(nbin)
      integer indra(nbin), nmin, nmax
      
c
c...  Loop over halos
c

c      ih = 0 
      do ic1 = 1 , nhalo 
        if ( nhp(ic1) .gt. nmin .and. nhp(ic1) .lt. nmax ) then 
          rhalo  = rhvir(ic1)
c         Look at a large enough sub-volume with 2 x rhalo on sides
          rhalo  = 2.0 * rhalo
          xhalo = xh(ic1) 
          yhalo = yh(ic1) 
          zhalo = zh(ic1) 

c....     convert from original grid to chaining mesh          
          xd    = xhalo - rhalo
          yd    = yhalo - rhalo
          zd    = zhalo - rhalo
          imin  = int(float(nll) * (xd - 1) / ng + 1)
          jmin  = int(float(nll) * (yd - 1) / ng + 1)
          kmin  = int(float(nll) * (zd - 1) / ng + 1)
          xd    = xhalo + rhalo
          yd    = yhalo + rhalo
          zd    = zhalo + rhalo
          imax  = int(float(nll) * (xd - 1) / ng + 1)
          jmax  = int(float(nll) * (yd - 1) / ng + 1)
          kmax  = int(float(nll) * (zd - 1) / ng + 1)

c....     sweep over neighbors checking periodic boundary conditions
          iph = 0
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idm = iCL2(ic,jc,kc)  ! read LL head
                do while ( idm .ne. nil )
                  xd = x(idm)
                  yd = y(idm)
                  zd = z(idm)
                  vxd = vx(idm)
                  vyd = vy(idm)
                  vzd = vz(idm)
c....             enforce periodic boundary conditions
                  diff_x = xd - xhalo                  
                  diff_y = yd - yhalo 
                  diff_z = zd - zhalo 
                  corr_x = zero 
                  corr_y = zero 
                  corr_z = zero 
                  
                  if ( abs(diff_x) .gt. ng/2 ) then 
                    if ( diff_x .gt. zero ) then 
                      corr_x = -ng 
                    else
                      corr_x =  ng
                    endif
                  endif
                  if ( abs(diff_y) .gt. ng/2 ) then 
                    if ( diff_y .gt. zero ) then 
                      corr_y = -ng 
                    else
                      corr_y =  ng
                    endif
                  endif
                  if ( abs(diff_z) .gt. ng/2 ) then 
                    if ( diff_z .gt. zero ) then 
                      corr_z = -ng 
                    else
                      corr_z =  ng
                    endif
                  endif
                  
                  
c....             compute distance particle - halo center squared
c                 corr_x , corr_y , corr_z - take care of periodicity
                  iph = iph + 1 
                  ra(iph) = sqrt( (xd - xhalo + corr_x)**2 + 
     &                        (yd - yhalo + corr_y)**2 + 
     &                        (zd - zhalo + corr_z)**2   )   
                  weight(iph) = pw2(idm)
c                  write(*,*) iph, ra(iph), weight(iph)

                  idm = iLL2(idm) ! next particle from linked list
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! end i
       
c...      Sort the array by the distance from the halo center
          call indexx ( iph , ra , indra )

c...      Count particle if it falls inside halo radius and delta < Delta           
c          do i = 1, iph
c             write(*,*) i,ra(i),indra(i),ra(indra(i)),weight(iph)
c          enddo
c          exit

          ir = 1
          amhalo = weight(indra(ir))
          overd = 1.0e+30
          do while ( overd .gt. Deltavir )  
             ir = ir + 1
             amhalo = amhalo + weight(indra(ir))
             rhalo = ra(indra(ir))
             vol = pi43 * rhalo**3
             overd = amhalo / vol / rhoaver
c             write(*,*) ir, indra(ir), ra(indra(ir)), overd
          enddo

          ! prepare for interpolation
          r1 = ra(indra(ir-1)) 
          r2 = ra(indra(ir))   
          vol1 = pi43 * r1**3
          vol2 = pi43 * r2**3
          overd1 = ( amhalo - weight(indra(ir-1)) ) / vol1 / rhoaver
          overd2 = amhalo / vol2 / rhoaver

          ! take log10 
          r1 = log10( r1 )
          r2 = log10( r2 )
          overd1 = log10( overd1 )
          overd2 = log10( overd2 )
          dvirlog = log10(Deltavir)

          ! interpolate linearly in log r - log delta
c          rhalo = 10.0**( (dvirlog * (r2-r1) + r1*overd2 - 
c     &                     r2*overd1) / (overd2 - overd1) ) 
          rhalo = (r2-r1)/(overd2-overd1)*(dvirlog-overd2) + r2
          rhalo = 10.**(rhalo)
 
          ! compute Rvir and Mvir in code units
          write(*,*) ic1, rhvir(ic1), amhvir(ic1)
          rhvir(ic1) = rhalo 
          amhvir(ic1) = pi43 * rhvir(ic1)**3 * rhoaver * Deltavir
c          write(*,*) ic1, rhvir(ic1), amhvir(ic1), 10.**r1, rhalo, 
c     &         10.**r2, 10.**overd1, 10.**overd2
          write(*,*) ic1, rhvir(ic1), amhvir(ic1)

        endif
      enddo ! end ic1

      return
      end
c
c     ------------------------------------------------------------------
      subroutine Write_Halo_Catalogs ( fname1, fname2, 
     &                               ibtype, rpmin, rpmax, npbin, nmin,
     &                               rvel, rhalo )
c     ------------------------------------------------------------------
      include 'a_hfind.h'
      real toohot 
      common / TOOHOT1 / toohot
      real dvdup
      common / VDUP / dvdup 
      integer ibtype, npbin
      character*256 fname1, fname2

      nfn1 = index ( fname1 , ' ' ) - 1
      nfn2 = index ( fname2 , ' ' ) - 1 

      if     ( ibtype .eq. 0 ) then 
        rmax   = rpmax
        rmin   = rpmin
        rmaxd  = log10(rmax)
        rmind  = log10(rmin)
        drd    = (rmaxd - rmind)/float(npbin)
        nbins  = int((rmaxd - rmind)/drd) + 1
      elseif ( ibtype .eq. 1 ) then 
        rmaxd  = rpmax
        rmind  = rpmin        
        drd    = (rpmax - rpmin)/float(npbin)
        nbins = npbin
      elseif ( ibtype .eq. 2 ) then 
        rmaxd = rpmax
        rmind = rpmin
        nbins = npbin
      else
        write(*,*) 'Compute_Halo_Properties: unknown bin type:',ibtype
        write(*,*) 'exiting...'
        return
      endif
      if ( nbins .gt. nbmax ) then 
        write(*,*) 'error : Compute_Halo_Properties : nbins > nbmax'
        write(*,*) 'nbins =',nbins,' nbmax =',nbmax
        stop
      endif

      open ( 33 , file = fname1(1:nfn1) )
      open ( 34 , file = fname2(1:nfn2) )
      write (33,105) HEADER,
     &               AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &               EKIN,
     &               NROWC,NGRIDC,Om0,Oml0,Omb0,hubble
      write (34,105) HEADER,
     &               AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &               EKIN,
     &               NROWC,NGRIDC,Om0,Oml0,Omb0,hubble

 105  format (' ',A45,/
     &        ' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &        ' I =',I4,' WEIGHT=',F8.3,' Ekin=',E12.3,/
     &        ' Nrow=',I4,' Ngrid=',I4,
     &        '  Omega_0=',F7.3,' OmLam_0=',F7.4,
     &        '  Omegab_0=',F7.3,' Hubble=',f7.3)
      pm1 = pw2(1) * pmmsun
      write(33,106) box0,rg2Mpc, vg2kms
      write(33,107) pmmsun,pm1
      write(33,120) DeltaMin
      write(33,121) Deltavir
      write(33,122) rhalo*rg2kpc
      write(33,108) rpmin*rg2kpc, rpmax*rg2kpc, npbin+1, ibtype
      write(33,109) rvel, toohot,dvdup
      write(33,110) nmin
      write(33,111) 
      write(33,112) 
      write(33,113) 
      write(33,114) 
      write(33,115) 
      write(34,106) box0,rg2Mpc, vg2kms
      write(34,107) pmmsun,pm1
      write(34,120) DeltaMin
      write(34,121) Deltavir
      write(34,122) rhalo*rg2kpc
      write(34,108) rpmin*rg2kpc, rpmax*rg2kpc, npbin+1, ibtype
      write(34,109) rvel, toohot, dvdup
      write(34,110) nmin
      write(34,116) 
      write(34,112) 
      write(34,113) 
      write(34,114) 
      write(34,115) 
      write(34,117) 
      write(34,118) 
      write(34,119) 

 120  format(' DeltaMin=',e10.3,' [min. overd. for density minima]')
 121  format(' Deltavir=',e10.3,' [virial overdensity]')
 122  format(' rhalo=',e9.2,' [/h kpc comoving=search radius]')
 106  format(' Lbox =',f7.2,' [/h Mpc comoving];',
     &       ' rg2Mpc=',e9.4,' [grid units->comoving /h Mpc];',
     &       ' vg2kms=',e9.4,' [grid units-> km/s peculiar velocity]')
 107  format(' M0 =',e9.4,' [/h Msun mass unit];',
     &       ' M(1)=',e9.4,' [/h Msun = mass of the 1st specie]')
 108  format(' rmin=',f8.3,'; rmax=',f8.3,' [/h kpc comoving];'
     &       ' nbins=',i4,' ibtype=',i2)
 109  format(' rvel=',f5.3,' [r/rhmax]; toohot=',f6.2,' dvdup=',f6.2)
 110  format(' min. number of particles =',i6)
 111  format(' column description:')
 112  format(' index x y z [/h Mpc comoving] vx vy vz [km/s peculiar]')
 113  format(' min(rt,rvir,rmax) min(rvir,rmax) [/h kpc comoving]')
 114  format(' M(<rh) M(<rvir) [/h Msun] Np(<rh)')
 115  format(' vmax [km/s] rmax r_s [/h kpc proper]')
 116  format(' column description for halo header:')
 117  format(' column description for halo profiles:')
 118  format(' rm rr [/h kpc comoving] d(<r) d(r) Vc(r) Vrms(r) [km/s]')
 119  format(' m(<r) [/h Msun] Np(r)=# of particles in the shell')
c
      ih = 0
c      write(*,*) 'Write_Halo_Catalog : ',nhalo, nhp(1)
      do ic1 = 1 , nhalo
        if ( nhp(ic1) .gt. nmin ) then 
        ih = ih + 1 
c        write(*,*) 'halo #',ic1,' nhp =',nhp(ic1)
        rkpc   = rh(ic1) * rg2kpc
        amassh = pmmsun * amh(ic1)
        amvir = pmmsun * amhvir(ic1)
        xmpc  = (xh(ic1)-1.) * rg2Mpc
        ympc  = (yh(ic1)-1.) * rg2Mpc
        zmpc  = (zh(ic1)-1.) * rg2Mpc
        vxkms = vxh(ic1) * vg2kms
        vykms = vyh(ic1) * vg2kms
        vzkms = vzh(ic1) * vg2kms
        rvir = rhvir(ic1) * rg2kpc

        write(33,94)  ih,ic1, xmpc ,  ympc ,  zmpc ,
     &                vxkms, vykms, vzkms , 
     &                rkpc, rvir, amassh, amvir, nhp(ic1), 
     &                vhmax(ic1), rhmax(ic1), rsh(ic1)
        write(34,94)  ih,ic1,xmpc ,  ympc ,  zmpc ,
     &                vxkms, vykms, vzkms , 
     &                rkpc, rvir, amassh, amvir, nhp(ic1), 
     &                vhmax(ic1), rhmax(ic1), rsh(ic1)
        
        rd = zero
        do ic2 = 0 , nbins
          call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                    rd1, rdm, rd2, volr, voli ) 
          rdkpc1 = rdm * rg2kpc
          rdkpc2 = rd2 * rg2kpc
          odi = pnt(ic2,ic1)/voli
          od = pn(ic2,ic1)/volr
          write(34,95) rdkpc1, rdkpc2, odi, od, vcirc(ic2,ic1),
     &       vrms(ic2,ic1)*vg2kms, pnt(ic2,ic1)*pmmsun, int(na(ic2,ic1))
        enddo
        endif
      enddo

 94   format(2(i6,1x),1x,3(1x,f7.3),1x,3(1x,f8.2),1x,2(1x,f7.2),
     &        1x,2(1x,e10.3),2x,i8,1x,f9.3,2(1x,f8.2))
 95   format(2(f8.3,1x),2x,2(e13.6,1x),1x,
     &       2(1x,f9.3),2x,e13.6,2x,i7)

      close ( 33 ) 
      close ( 34 ) 

      return
      end
c
c     ------------------------------------------------------------
      subroutine Write_Halo_Particles( fname, fname2, 
     &                                 ibtype, rpmin, rpmax, rvel, 
     &                                 npbin, imark, nmin )
c     ------------------------------------------------------------
c
c     ibtype = 0 - log10 bins, 1 - even bins
c     rpmin,max - min/max. radius for profile construction in grid units
c     imark = 1 - mark particles as they assigned to a halo, so that they
c                 they are not included in subsequent halo profiles
c                (this will prevent detection of most of the subhalos but
c                 will greatly speed-up finding of isolated halos)
c     imark = 0 - don't
c     rvel = radius within which to measure mean halo velocity 
c     miter = number of iterations for unbound particles 
c
      include 'a_hfind.h'
      character*256 fname,fname2
      real toohot ! should be initialized in Compute_Halo_Properties
      common / TOOHOT1 / toohot
c      parameter ( toohot = 2. )
      integer ibtype, npbin, miter

      nhd = 0 
      do ic1 = 1 , nhalo 
        if ( nhp(ic1) .gt. nmin ) then 
          nhd = nhd + 1
        endif
      enddo

      nfn1 = index ( fname , ' ' ) - 1
      nfn2 = index ( fname2 , ' ' ) - 1
      open ( 20 , file = fname(1:nfn1), form='unformatted' )
c      open ( 22 , file = fname2(1:nfn2), form='unformatted' )
      write(20) aexpn
      write(20) nhd, np1
c      write(22) aexpn
c      write(22) nhd, np1

      do ic1 = 1 , np1max 
        iSp(ic1) = nil 
      enddo

      ih = 0 
      do ic1 = 1 , nhalo 
        if ( nhp(ic1) .gt. nmin ) then 
          rhalo  = rh(ic1)
          rhmaxd = rhmax(ic1)/rg2pkpc
          ih = ih + 1
          write(20) ih
c          write(22) ih
          xhalo = xh(ic1) 
          yhalo = yh(ic1) 
          zhalo = zh(ic1) 

          vxhalo = vxh(ic1)
          vyhalo = vyh(ic1)
          vzhalo = vzh(ic1)
c....     convert from original grid to chaining mesh          
          xd    = xhalo - rhalo
          yd    = yhalo - rhalo
          zd    = zhalo - rhalo
          imin  = int(float(nll) * (xd - 1) / ng + 1)
          jmin  = int(float(nll) * (yd - 1) / ng + 1)
          kmin  = int(float(nll) * (zd - 1) / ng + 1)
          xd    = xhalo + rhalo
          yd    = yhalo + rhalo
          zd    = zhalo + rhalo
          imax  = int(float(nll) * (xd - 1) / ng + 1)
          jmax  = int(float(nll) * (yd - 1) / ng + 1)
          kmax  = int(float(nll) * (zd - 1) / ng + 1)
c....     sweep over neighbors checking periodic boundary conditions
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idm = iCL2(ic,jc,kc)                 ! read LL head
                do while ( idm .ne. nil )
                  if ( iSp(idm) .eq. nil ) then 
                  xd = x(idm)
                  yd = y(idm)
                  zd = z(idm)
                  vxd = vx(idm)
                  vyd = vy(idm)
                  vzd = vz(idm)
c....             enforce periodic boundary conditions
                  diff_x = xd - xhalo                  
                  diff_y = yd - yhalo 
                  diff_z = zd - zhalo 
                  corr_x = zero 
                  corr_y = zero 
                  corr_z = zero 

                  if ( abs(diff_x) .gt. ng/2 ) then 
                    if ( diff_x .gt. zero ) then 
                      corr_x = -ng 
                    else
                      corr_x =  ng
                    endif
                  endif
                  if ( abs(diff_y) .gt. ng/2 ) then 
                    if ( diff_y .gt. zero ) then 
                      corr_y = -ng 
                    else
                      corr_y =  ng
                    endif
                  endif
                  if ( abs(diff_z) .gt. ng/2 ) then 
                    if ( diff_z .gt. zero ) then 
                      corr_z = -ng 
                    else
                      corr_z =  ng
                    endif
                  endif

c....             compute distance particle - halo center squared
c                 corr_x , corr_y , corr_z - take care of periodicity
                  rd  = sqrt( (xd - xhalo + corr_x)**2 + 
     &                        (yd - yhalo + corr_y)**2 + 
     &                        (zd - zhalo + corr_z)**2   )   
                  rd = max(rd,1.e-7)
                  dv = sqrt( (vxd - vxhalo)**2 +
     &                       (vyd - vyhalo)**2 + 
     &                       (vzd - vzhalo)**2   ) * vg2kms                    
                  resc = rd * rg2pkpc
                  vesc = ve ( resc , rhmax(ic1) , vhmax(ic1) )
                  vesc = vesc * toohot 
c
c....             count particle if it falls inside halo radius and v < vesc 
c                   
                  if ( rd .lt. rhalo .and. dv .lt. vesc ) then             
                    write(20) idm
c                    if ( ih .eq. 1 ) then 
c                      write(*,*) idm,x(idm),y(idm),z(idm)
c                    endif
                  endif
c                  if ( rd .lt. rhmaxd .and. dv .lt. vesc ) then  
c                    write(22) idm
c                  endif
                  endif 
                  idm = iLL2(idm) ! next particle from linked list
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! end i
          write(20) nil
c         write(22) nil
        endif
      enddo ! end ic1
c
      close ( 20 ) 
c      close ( 22 ) 

      return
      end
c
c     ------------------------------------------------------------
      subroutine Write_Halo_Particles_Convert ( fname, fname2 )
c     ------------------------------------------------------------
c
      include 'a_hfind.h'
      character*256 fname, fname2

      nfn1 = index ( fname , ' ' ) - 1
      nfn2 = index ( fname2 , ' ' ) - 1
      open ( 21 , file = fname(1:nfn1), form='unformatted' )
      open ( 22 , file = fname2(1:nfn2), form='unformatted' )
      read(21) aexpnn
      read(21) nhd, npp
      write(22) aexpnn
      write(22) nhd, npp

      ih = 0 
      do ic1 = 1 , nhd
        read(21) ih
        idm = 1
        nhp1 = 0 
        do while ( idm .ne. 0 ) 
          read(21) idm
          if ( idm .ne. 0 ) then 
            nhp1 = nhp1 + 1
            ind2(nhp1) = idm 
          endif
        enddo
        write(22) ih, nhp1, (ind2(i),i=1,nhp1)
c        write(*,*) ih, nhp1, (ind2(i),i=1,nhp1)
      enddo  ! end ic1
c
      close ( 21 ) 
      close ( 22 ) 

      return
      end
c
c     ------------------------
      subroutine InitArrays ()
c     ------------------------
c     purpose: to initialize working arrays
c     -------------------------------------
      include 'a_hfind.h'

      do ic1 = 1 , nll
        do ic2 = 1 , nll
          do ic3 = 1 , nll 
            iCL2(ic3,ic2,ic1) = nil 
          enddo
        enddo
      enddo

      do ic1 = 1 , np1max 
        iLL2(ic1) = nil 
        iSp(ic1)  = nil 
      enddo

      do ic1 = 1 , nhmax 
        iSh(ic1) = nil 
        nhp(ic1) = nil 
        rhvir(ic1) = 0.0
        rh(ic1) = 0.0
        rsh(ic1) = 0.0 
        amh(ic1) = 0.0
        amhvir(ic1) = 0.0
      enddo

      return
      end

c     --------------------------------
      subroutine Conversions ( boxh0 ) 
c     --------------------------------
c     purpose: sets up various grid <-> physical 
c              conversion factors
c     input  : boxh0 - box size in h^{-1} Mpc at z=0
c     -------------------------------------------------------
      include 'a_hfind.h'
      real boxh0  
c
c.... particle mass in h^{-1} solar masses (CDM,LCDM)
c     divide by hubble to get physical mass
c
      fb = Omb0 / Om0 ! 1/(1-fb) to correct for baryon mass in hydro runs
      pmmsun = 2.760368e11 * Om0 * 
     &         (boxh0)**3 / (ng**3) ! must be in /h Msun !
      pm1 = pw2(1) * pmmsun ! mass of the smallest particles
      box0 = boxh0

c
c.... grid scale <-> physical scale conversions
c     using boxh will give units in h^{-1} 
c
      rMpc2g  = float(ng)/boxh0    ! conversion /h Mpc -> grid units
      rkpc2g  = rMpc2g / 1000.0    ! conversion /h kpc -> grid units
      rg2Mpc  = boxh0 / ng         ! conversion grid units -> /h Mpc
      rg2kpc  = 1000. * boxh0 / ng ! conversion grid units -> /h kpc
      rg2pMpc = boxh0 * aexpn / ng ! conversion grid units -> proper /h Mpc
      rg2pkpc = 1000. * boxh0 * aexpn / ng ! conversion grid units -> proper /h kpc
c
c.... grid velocity <-> km/s = (x0*Ho)/aexpn
c
c      vg2kms = 100.0*boxh0 / (aexpn*ngridc) ! ngridc instead of ng (for N-body sim)
      vg2kms = 50.0*boxh0 / (aexpn*ng) * Om0**0.5 ! (for Hydro sim)
      write(*,*) 'Lbox =',boxh0,' [/h Mpc]'
      write(*,*) 'M0 =',pmmsun,' [/h Msun]'
      write(*,*) 'M(ispec=1) =',pm1,' [/h Msun]'
      write(*,*) 'pw2(1) =',pw2(1)
      write(*,*) 'rkpc2g =',rkpc2g
      write(*,*) 'vg2kms =',vg2kms
      return
      end

c     -----------------------------
      subroutine LL_Construct_HF ()
c     -----------------------------
c     constructs particle linked list 
c     see description of algorithms in  
c       T.H.Cormen , C.E.Leiserson , R.L.Rivest
c       "Introduction to algorithms" pp. 204-208
c     ------------------------------------------
      include 'a_hfind.h'
      
c.... construct linked list
      do ic1 = 1 , np1
        ip            = int(float(nll)*(x(ic1)-1.0)/float(ng)+1)
        jp            = int(float(nll)*(y(ic1)-1.0)/float(ng)+1)
        kp            = int(float(nll)*(z(ic1)-1.0)/float(ng)+1)
        iLL2(ic1)      = iCL2(ip,jp,kp)
        iCL2(ip,jp,kp) = ic1        
      enddo
      return
      end

c     -------------------------------
      subroutine Assign_Density_HF ()
c     -------------------------------
c     assigns potential at the locations 
c     of the first mass species
c     ----------------------------
      include 'a_hfind.h'

      dmin = 1.0e20
      dmax = -1.0e20
      do ic1 = 1 , np1 
        dmin = min( ddummy(ic1), dmin )
        dmax = max( ddummy(ic1), dmax )
        dnb2(ic1) = ddummy(ic1)
        if ( dnb2(ic1) .le. 1.0e-10 ) then
           call Open_ASCII_File ( iOErr , errorfile , sequent ) 
           write(iOErr,*)'1  *error in Assign_Density_HF :'
           write(iOErr,*)'   * dnb2(ic1) < 1.0e-10 '
           write(iOErr,*)'   * check ddummy() array '
           close (iOErr)
        endif
      enddo
      write(*,*) 'dmin, dmax =',dmin, dmax

      return
      end
c
c     ----------------------------------------
      subroutine Write_Density_HF ( ismooth )
c     ----------------------------------------
c     write density field @ particle locations 
c     ----------------------------------------
      include 'a_hfind.h'
      character*256 fname
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      
      ic = int(aexpn*1000)
      i1 = ic / 1000 
      i2 = (ic - (ic/1000)*1000)/100 
      i3 = (ic - (ic/100)*100) / 10 
      i4 = (ic - (ic/10)*10)
      if ( ismooth .eq. 1 ) then
         fname  = 'rho_smooth_a'//
     &        digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &        //'.den'
      else
         fname  = 'rhotot_a'//
     &        digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &        //'.den'
      endif

      nfn = index ( fname , ' ' ) - 1 
      write(*,*) ' '
      write(*,*) 'writing total density : ',fname(1:nfn)
      
      open(21, file = fname(1:nfn),form='unformatted')
      write(21) np1
      write(21) aexpn
      write(21) (dnb2(i),i=1,np1)
      close (21)

      write(*,*) 'np1, aexpn =',np1,aexpn

      return
      end
c
c     ---------------------------------------
      subroutine Write_Particles_HF ( fname )
c     ---------------------------------------
c     write particle locations 
c     ---------------------------------------
      include 'a_hfind.h'
      character*256 fname
c
      nfn = index ( fname , ' ' ) - 1 
      write(*,*) ' '
      write(*,*) 'writing particle locations...'

      open(21, file = fname(1:nfn),form='unformatted')
      write(21) np1
      write(21) aexpn
      write(21) (x(i),i=1,np1)
      write(21) (y(i),i=1,np1)
      write(21) (z(i),i=1,np1)
      close (21)

      write(*,*) 'np1, aexpn =',np1,aexpn

      return
      end
c
c     ---------------------------
      subroutine SortParticles ()
c     ---------------------------
c     sort particles according to their local density/potential
c     ind2(np1) has index of particle with largest local density/
c     smallest local potential
c
c     uses: indexx subroutine from NR
c     -----------------------------------------------------
      include 'a_hfind.h'

      call indexx ( np1 , dnb2 , ind2 )

      return
      end

c     -------------------------------
      subroutine FindHaloes ( rhalo )
c     -------------------------------
c     purpose  : finds initial halo centers 
c
c     algorithm: particle with the smallest local potential (in dnb2)
c                determines the approximate location of the first 
c                halo center. Particles located inside a sphere
c                of radius rhalo centered at the halo center are assigned 
c                to the same halo and are excluded from  the list of 
c                particles used to identify halos (iSp is set to 1). 
c                The radius rhalo is an adjustable parameter; its values
c                should be of order of spatial resolution of a simulation. 
c                The procedure repeats for the particle with the smallest 
c                local potential in the list of remaining particles. The peaks 
c                are identified until there are no particles in the list.
c     ----------------------------------------------------------------------
      include 'a_hfind.h'

      nhalo = nil       
c
c.... the particles are sorted so that ind2(np1) has index of particle 
c     with largest local density - this determines the inverse order 
c     in the loop below

      rhalo2  = rhalo**2

      do ic1 = np1 , 1 , -1  
        if (mod(ic1,50000).eq.nil) then 
          write(*,*) ic1,' processed delta=',dnb2(ind2(ic1))
        endif

c        if ( dnb2(ind2(ic1)) .gt. DeltaMin ) then
c           write(*,*) dnb2(ind2(ic1)), iSP(ind2(ic1))
c        endif
c....   determine center of the next halo
        if ( iSp(ind2(ic1)) .eq. nil .and. 
     &       dnb2(ind2(ic1)) .gt. DeltaMin ) then
c          write(*,*) ic1, dnb2(ind2(ic1)), rhalo 
          iCenter      = ind2(ic1)     
          dcen         = dnb2(ind2(ic1))
          xhalo        = x(iCenter)
          yhalo        = y(iCenter)
          zhalo        = z(iCenter)
c....     convert from original grid to chaining mesh          
          ihalo        = int(float(nll) * (xhalo - 1) / ng + 1)
          jhalo        = int(float(nll) * (yhalo - 1) / ng + 1)
          khalo        = int(float(nll) * (zhalo - 1) / ng + 1)
          xd           = xhalo - rhalo
          yd           = yhalo - rhalo
          zd           = zhalo - rhalo
          imin         = int(float(nll) * (xd - 1) / ng + 1)
          jmin         = int(float(nll) * (yd - 1) / ng + 1)
          kmin         = int(float(nll) * (zd - 1) / ng + 1)
          xd           = xhalo + rhalo
          yd           = yhalo + rhalo
          zd           = zhalo + rhalo
          imax         = int(float(nll) * (xd - 1) / ng + 1)
          jmax         = int(float(nll) * (yd - 1) / ng + 1)
          kmax         = int(float(nll) * (zd - 1) / ng + 1)
c....     sweep over neighbors checking periodic boundary conditions
          xcm = 0.
          ycm = 0. 
          zcm = 0. 
          ddmax = -1.e6
          wdummy = nil                           ! # of halo particles
          n10    = nil 
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idummy = iCL2(ic,jc,kc)                 ! read LL head
                do while ( idummy .ne. nil )
c                  write(*,*) ic1, idummy, iLL2(idummy)
                  if ( iSp(idummy) .eq. nil ) then
                    xd = x(idummy)
                    yd = y(idummy)
                    zd = z(idummy)
                    diff_x = xd - xhalo 
                    diff_y = yd - yhalo 
                    diff_z = zd - zhalo 
                    corr_x = zero 
                    corr_y = zero 
                    corr_z = zero 
                    if ( abs(diff_x) .gt. ng/2 ) then 
                      if ( diff_x .gt. zero ) then 
                        corr_x = -ng 
                      else
                        corr_x = ng
                      endif
                    endif
                    if ( abs(diff_y) .gt. ng/2 ) then 
                      if ( diff_y .gt. zero ) then 
                        corr_y = -ng 
                      else
                        corr_y = ng
                      endif
                    endif
                    if ( abs(diff_z) .gt. ng/2 ) then 
                      if ( diff_z .gt. zero ) then 
                        corr_z = -ng 
                      else
                        corr_z = ng
                      endif
                    endif

                    rd = (xd - xhalo + corr_x)**2 + 
     &                   (yd - yhalo + corr_y)**2 + 
     &                   (zd - zhalo + corr_z)**2 
c                    write(*,*) rd, xd, yd, zd, xhalo, yhalo, zhalo

                    if ( rd .le. rhalo2 ) then
                      wdummy = wdummy + pw2(idummy)
                      xcm = xcm + x(idummy)*pw2(idummy)
                      ycm = ycm + y(idummy)*pw2(idummy)
                      zcm = zcm + z(idummy)*pw2(idummy)
                      ddmax = max(ddmax,dnb2(idummy))
                      iSp(idummy) = 1                      
                    endif
                  endif
                  idummy = iLL2(idummy)
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! and i
          vol   = 4.0 * pi * rhalo**3 / 3.0
          overd = wdummy / vol / rhoaver

c          write(*,*)ic1, nhalo, wdummy, int(wdummy/pw2(1)),overd, rhoaver
          if ( (wdummy.gt.0.) .and. dcen .ge. ddmax ) then
            nhalo = nhalo + 1
            if ( nhalo .gt. nhmax ) then 
              write(*,*) '* error: nhmax is too small: nh =',nhalo
              write(*,*) '* while more halos are found. stopping...'
              stop
            endif
c            xh(nhalo) = xcm / ndummy 
c            yh(nhalo) = ycm / ndummy
c            zh(nhalo) = zcm / ndummy
            xh(nhalo)  = xhalo
            yh(nhalo)  = yhalo
            zh(nhalo)  = zhalo
            rh(nhalo)  = rhalo 
            nhp(nhalo) = int(wdummy/pw2(1))
            amh(nhalo) = wdummy
          endif
        endif
      enddo

      write(*,*) ' nhalo =', nhalo 

      return
      end

c     ------------------------------------------
      subroutine IterateHaloCM ( r , rminshift )
c     ------------------------------------------
c     purpose: finds centers of mass for haloes iteratively by 
c              finding the center of mass of all particles inside r and 
c              displacing the center of the sphere to the center of mass. 
c              The procedure is iterated until convergence (displacement 
c              is less than rminshift or if number of particles inside r 
c              starts to decrease).
c     -----------------------------------------------------
      include 'a_hfind.h'
      
      real r , rminshift 
      real xt , yt , zt , xnew , ynew , znew , riter
      
      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
         rh(ic1) = r
        nhp(ic1) = 0 
      enddo

      riter = r
      miter = nil 
      nstop = nil 
 7    iter  = nil 
      miter = miter + 1
      do ic1 = 1 , nhalo 
        xt = xh(ic1) 
        yt = yh(ic1) 
        zt = zh(ic1) 
        nt = nhp(ic1)
        it = iSh(ic1)
        
        if ( it .eq. nil ) then 
          call CM ( xt , yt , zt , xnew , ynew , znew , 
     &              vxnew , vynew , vznew , nnew , riter )
          rshift = max (abs(xnew-xt) , abs(ynew-yt) , abs(znew-zt))
          if ((rshift.lt.rminshift) .or. (nnew.le.nt)) then 
            iSh(ic1) = 1
          else
            iter = iter + 1
            if ( xnew .lt. 1.0 ) xnew = xnew + ng 
            if ( ynew .lt. 1.0 ) ynew = ynew + ng 
            if ( znew .lt. 1.0 ) znew = znew + ng 
            if ( xnew .gt. xn  ) xnew = xnew - ng 
            if ( ynew .gt. xn  ) ynew = ynew - ng 
            if ( znew .gt. xn  ) znew = znew - ng 
            xh(ic1)  = xnew 
            yh(ic1)  = ynew 
            zh(ic1)  = znew 
            if ( nnew .gt. 0 ) then 
              vxh(ic1) = vxnew 
              vyh(ic1) = vynew 
              vzh(ic1) = vznew 
              nhp(ic1) = nnew       
            endif
          endif
        endif
      enddo
      write(*,*) 'iterating halo mass centers ', iter, 'left'
      if ((iter .gt. 0).and.(miter.lt.100)) goto 7 
      
      return
      end

c     -------------------------------------------------------------
      subroutine CM ( xold , yold , zold , 
     &                xnew , ynew , znew , vxnew , vynew , vznew , 
     &                nnew , r )
c     -------------------------------------------------------------
      include 'a_hfind.h'
      real xold , yold , zold , xnew , ynew , znew , r
      
      r2 = r**2
      
      xnew   = 0.0
      ynew   = 0.0
      znew   = 0.0 
      vxnew  = 0.0 
      vynew  = 0.0 
      vznew  = 0.0 

      wdummy = 0.0

c.... convert from original grid to chaining mesh          
      xd    = xold - r
      yd    = yold - r
      zd    = zold - r
      imin  = int(float(nll) * (xd - 1) / ng + 1)
      jmin  = int(float(nll) * (yd - 1) / ng + 1)
      kmin  = int(float(nll) * (zd - 1) / ng + 1)
      xd    = xold + r
      yd    = yold + r
      zd    = zold + r
      imax  = int(float(nll) * (xd - 1) / ng + 1)
      jmax  = int(float(nll) * (yd - 1) / ng + 1)
      kmax  = int(float(nll) * (zd - 1) / ng + 1)
c.... sweep over neighbors checking periodic boundary conditions
      ndummy = nil                           ! # of halo particles
      do i = imin , imax
        ic = i 
        if ( ic .lt. 1   ) ic = ic + nll
        if ( ic .gt. nll ) ic = ic - nll
        do j = jmin , jmax 
          jc = j 
          if ( jc .lt. 1   ) jc = jc + nll
          if ( jc .gt. nll ) jc = jc - nll
          do k = kmin , kmax               
            kc = k 
            if ( kc .lt. 1   ) kc = kc + nll
            if ( kc .gt. nll ) kc = kc - nll 
            idummy = iCL2(ic,jc,kc)                 ! read LL head
            do while ( idummy .ne. nil )
              xd = x(idummy)
              yd = y(idummy)
              zd = z(idummy)
              diff_x = xd - xold
              diff_y = yd - yold
              diff_z = zd - zold
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ng/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ng 
                else
                  corr_x = ng
                endif
              endif
              if ( abs(diff_y) .gt. ng/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ng 
                else
                  corr_y = ng
                endif
              endif
              if ( abs(diff_z) .gt. ng/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ng 
                else
                  corr_z = ng
                endif
              endif

              rd = (xd - xold + corr_x)**2 + 
     &             (yd - yold + corr_y)**2 + 
     &             (zd - zold + corr_z)**2 

              if ( rd .le. r2 ) then
                 wdummy = wdummy + pw2(idummy)
                 xnew   = xnew  +  x(idummy)*pw2(idummy)
                 ynew   = ynew  +  y(idummy)*pw2(idummy)
                 znew   = znew  +  z(idummy)*pw2(idummy)
                 vxnew  = vxnew + vx(idummy)*pw2(idummy)
                 vynew  = vynew + vy(idummy)*pw2(idummy)
                 vznew  = vznew + vz(idummy)*pw2(idummy)
              endif
              idummy = iLL2(idummy)
            enddo  ! end do while 
          enddo  ! end k 
        enddo  ! end j 
      enddo  ! and i            
      
      if (wdummy .eq. 0.0 ) then 
        xnew  =  xold
        ynew  =  yold
        znew  =  zold
        nnew  = 0
      else
       xnew  =  xnew / wdummy 
       ynew  =  ynew / wdummy 
       znew  =  znew / wdummy 
       vxnew = vxnew / wdummy
       vynew = vynew / wdummy
       vznew = vznew / wdummy
       nnew  = int(wdummy/pw2(1))
      endif
      

      return
      end

c     -------------------------------
      subroutine RemoveSmall ( nmin )
c     -------------------------------
c     purpose: removes small (with nhp < nmin) halos
c     input  : nmin - keep only halos with nhp > nmin 
c     ---------------------------------------------
      include 'a_hfind.h'
      integer nmin 

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( nhp(ic1) .gt. nmin ) then 
          nnew      = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      write(*,*) nhalo-nnew, ' small haloes were removed'

      nhalo = nnew

      return
      end

c     ----------------------------------------
      subroutine Remove_Velocity_Duplicates ()
c     ----------------------------------------
c     neighboring haloes eat each other 
c     ---------------------------------
      include 'a_hfind.h'
      common / VDUP / dvdup 
      dvdup = 0.1 
c
      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
      enddo 

      do ic1 = 1 , nhalo-1
c        write(*,*) 'processed',ic1,' halos'
        if ( iSh(ic1) .eq. nil ) then
          vxhd1 = vxh(ic1)
          vyhd1 = vyh(ic1)
          vzhd1 = vzh(ic1)
          vhdumi = 1.0 / sqrt(vxhd1**2 + vyhd1**2 + vzhd1**2)
          do ic2 = ic1+1 , nhalo 
            if ( iSh(ic2) .eq. nil ) then 
              vxhd2 = vxh(ic2)
              vyhd2 = vyh(ic2)
              vzhd2 = vzh(ic2)
              diff_x = xh(ic2) - xh(ic1) 
              diff_y = yh(ic2) - yh(ic1) 
              diff_z = zh(ic2) - zh(ic1)
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ng/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ng 
                else
                  corr_x = ng
                endif
              endif
              if ( abs(diff_y) .gt. ng/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ng 
                else
                  corr_y = ng
                endif
              endif
              if ( abs(diff_z) .gt. ng/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ng 
                else
                  corr_z = ng
                endif
              endif
              rd = (xh(ic2) - xh(ic1) + corr_x)**2 +
     &             (yh(ic2) - yh(ic1) + corr_y)**2 +
     &             (zh(ic2) - zh(ic1) + corr_z)**2 
              rhmaxd = min(rhmax(ic1),rhmax(ic2))
              
              dv = sqrt((vxhd1-vxhd2)**2 + 
     &                  (vyhd1-vyhd2)**2 + 
     &                  (vzhd1-vzhd2)**2 ) * vhdumi 
              if ( sqrt(rd)*rg2pkpc .lt. rhmaxd .and. 
     &             dv .lt. dvdup ) then
                if ( vhmax(ic2) .le. vhmax(ic1) ) then
                  iSh(ic2) = 1 
                else
                  iSh(ic1) = 1 
                  go to 5 
                endif
              endif
            endif
          enddo
        endif
 5      continue 
      enddo

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then 
          nnew = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      nhalo = nnew 

      write(*,*) 'Removing vel. duplicates: ',nhalo, ' haloes survived'

      return
      end

c     -------------------------
      subroutine Cannibalism ()
c     -------------------------
c     neighboring haloes eat each other 
c     ---------------------------------
      include 'a_hfind.h'

      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
      enddo 

      do ic1 = 1 , nhalo-1
c        write(*,*) 'processed',ic1,' halos'
        if ( iSh(ic1) .eq. nil ) then
          do ic2 = ic1+1 , nhalo 
            if ( iSh(ic2) .eq. nil ) then 
              diff_x = xh(ic2) - xh(ic1) 
              diff_y = yh(ic2) - yh(ic1) 
              diff_z = zh(ic2) - zh(ic1)
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ng/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ng 
                else
                  corr_x = ng
                endif
              endif
              if ( abs(diff_y) .gt. ng/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ng 
                else
                  corr_y = ng
                endif
              endif
              if ( abs(diff_z) .gt. ng/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ng 
                else
                  corr_z = ng
                endif
              endif
              rd = (xh(ic2) - xh(ic1) + corr_x)**2 +
     &             (yh(ic2) - yh(ic1) + corr_y)**2 +
     &             (zh(ic2) - zh(ic1) + corr_z)**2 
              rmax = max(rh(ic1),rh(ic2))
              if ( sqrt(rd) .lt. rmax ) then
                if ( nhp(ic2) .le. nhp(ic1) ) then
                  iSh(ic2) = 1 
                else
                  iSh(ic1) = 1 
                  go to 5 
                endif
              endif
            endif
          enddo
        endif
 5      continue 
      enddo

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then 
          nnew = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      nhalo = nnew 

      write(*,*) 'cannibalism: ',nhalo, ' haloes survived'

      return
      end

c     ----------------------------
      subroutine Cannibalism200 ()
c     ----------------------------
c     neighboring haloes eat each other after growing 
c     -----------------------------------------------
      include 'a_hfind.h'

      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
      enddo 

      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then
          do ic2 = 1 , nhalo 
            if ( (iSh(ic2) .eq. nil)
     &                 .and. 
     &           ( ic1 .ne. ic2   ) ) then 
              diff_x = xh(ic2) - xh(ic1) 
              diff_y = yh(ic2) - yh(ic1) 
              diff_z = zh(ic2) - zh(ic1)
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ng/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ng 
                else
                  corr_x = ng
                endif
              endif
              if ( abs(diff_y) .gt. ng/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ng 
                else
                  corr_y = ng
                endif
              endif
              if ( abs(diff_z) .gt. ng/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ng 
                else
                  corr_z = ng
                endif
              endif
              rd = (xh(ic2) - xh(ic1) + corr_x)**2 +
     &             (yh(ic2) - yh(ic1) + corr_y)**2 +
     &             (zh(ic2) - zh(ic1) + corr_z)**2 
              if ( rd .lt. rhvir(ic1)**2 ) then
                if ( amh(ic2) .le. amh(ic1) ) then
                  iSh(ic2) = 1 
                else
                  iSh(ic1) = 1 
                  go to 5 
                endif
              endif
            endif
          enddo
        endif
 5      continue 
      enddo

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then 
          nnew = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      nhalo = nnew 

      write(*,*) 'cannibalism: ',nhalo, ' haloes survived'

      close (40) 

      return
      end

c     -----------------------------------
      subroutine HOld2New ( nold , nnew )
c     -----------------------------------
c     copy halo variables from nold to nnew location
c     -----------------------------------------------
      include 'a_hfind.h'
c
      xh(nnew)  = xh(nold) 
      yh(nnew)  = yh(nold) 
      zh(nnew)  = zh(nold) 
      vxh(nnew) = vxh(nold) 
      vyh(nnew) = vyh(nold) 
      vzh(nnew) = vzh(nold) 
      rh(nnew)  = rh(nold)
      rsh(nnew)  = rsh(nold)
      rhvir(nnew) = rhvir(nold)
      amh(nnew) = amh(nold)
      amhvir(nnew) = amhvir(nold)
      rst(nnew)  = rst(nold)
      nhp(nnew) = nhp(nold) 
      vhmax(nnew) = vhmax(nold) 
      rhmax(nnew) = rhmax(nold) 
      rhvir(nnew) = rhvir(nold) 
      iHP(nnew) = iHP(nold) 
      iLH(nnew) = iLH(nold) 
      do i = 0 , nbmax
        pn(i,nnew) = pn(i,nold)
        pnt(i,nnew) = pnt(i,nold)
        na(i,nnew) = na(i,nold)
        pvx(i,nnew) = pvx(i,nold)
        pvy(i,nnew) = pvy(i,nold)
        pvz(i,nnew) = pvz(i,nold)
        vcirc(i,nnew) = vcirc(i,nold)
        vrms(i,nnew) = vrms(i,nold)
      enddo
c
      return
      end

c     -------------------------------
      function ve ( r , rmax , vmax )
c     -------------------------------
c
c     computes escape velocity for a NFW profile
c     c is approximated roughly
c
c.... compute escape velocity
c
      x = r/rmax
      ve = 2.15 * vmax 
     &          * sqrt(
     &                  log(1.0+2.0*x) / x
     &                )       
      return
      end

c     --------------------
      function con ( vmass )
c     --------------------
c     
c     concentration for mass vmass
c
      dimension cm(9)
c
c     approximation for c(M):
c
c     M < 10**   11.0  11.5  12.0  12.5  13.0  13.5  14.0  14.5  15.0
      data cm /  20.0, 17.8, 16.0, 14.3, 12.3, 11.0,  9.6,  8.3,  7.2 /
c
c.... approximate c
c
      if ( vmass .le. 0.5e11 ) then
        con  = 20.0
      else
        im = int(log10(vmass)/0.5) - 20
        if ( im .gt. 9 ) then
          con = 7.0
        else
          con = cm(im)
        endif
      endif

      return
      end

c     --------------------------
      subroutine HL_Construct ()
c     --------------------------
c     constructs halo linked list 
c     see description of algorithms in  
c       T.H.Cormen , C.E.Leiserson , R.L.Rivest
c       "Introduction to algorithms" pp. 204-208
c     ------------------------------------------
      include 'a_hfind.h'
      
c.... construct halo linked list

      do ic1 = 1 , nhalo 
        ih            = int(float(nll) * (xh(ic1) - 1) / ng + 1)
        jh            = int(float(nll) * (yh(ic1) - 1) / ng + 1)
        kh            = int(float(nll) * (zh(ic1) - 1) / ng + 1)
        iLH(ic1)      = iCH(ih,jh,kh)
        iCH(ih,jh,kh) = ic1
      enddo

      return
      end

c     ------------------------------------------------------------
      subroutine Read_Particles_Binary_HF ( isim,fname1,fname2,fname3 )
c     ------------------------------------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
      include 'a_hfind.h'
      character*256 fname1, fname2, fname3
c
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
c 
c      ngrid = ng

      nfn1 = index ( fname1 , ' ' ) - 1 
      nfn2 = index ( fname2 , ' ' ) - 1
      nfn3 = index ( fname3 , ' ' ) - 1

      open ( 3 ,file =fname1(1:nfn1), form = 'unformatted')

c.... read control information and check whether it has proper structure
c.... to read N-body ART, get rid of Omb0 read below
      Omb0 = 0.
      if ( isim .eq. 0 ) then
         read      (3) HEADER, 
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,
     &                  hubble,Wp5
     &                   ,Ocurv,extras
      else 
         read      (3) HEADER, 
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,
     &                  hubble,Wp5
     &                   ,Ocurv,Omb0,extras
      endif

      fb = Omb0/Om0

      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &               ' Omegab_0=',F7.3,' Hubble=',f7.3)
      if ( NGRIDC .ne. NGRID ) then 
        write(*,*) 'NGRIDC .ne. NGRID :',NGRIDC, NGRID
        write(*,*) 'hope this is ok...'
      endif
      if(nspecies .eq. 0 ) nspecies = 1
      If( np1max .lt. lspecies(nspecies) ) then
      write (*,*) ' Wrong number of particles: '
      write (*,*) ' must be at least=',lspecies(nspecies),' (lspecies)'
      write (*,*) ' but is set to ',np1max,' in a_hfind.h...'
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo
         STOP
      Endif 
      
      nbyte  = nrecl * 4
      nacces = nbyte / nbyteword
 
      open ( 1 , file = fname2(1:nfn2), access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3

      N_particles = lspecies(nspecies)   ! Total number of particles

      Npages      = (N_particles -1)/npage + 1
      N_in_last   = N_particles - npage*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) 'N_particles =',N_particles, ' N_in_last=',N_in_last
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec), wspecies(ispec)
         enddo

      do i = 1 , NSpecies
        if(i .eq. 1 ) then
          nb1 = 1
        else
          nb1 = lspecies(i-1) + 1
        endif 
        nb2 = lspecies(i)
        do ip = nb1 , nb2
c         correct for baryon fraction (Daisuke, 4/25/03)
          pw2(ip) = wspecies(i) / (1.0-fb)
c          write(*,*) wspecies(i),pw2(ip),fb
        enddo
      enddo 
      np1 = lspecies(1) ! only the first specie
      open (21,file = fname3(1:nfn3),
     &         form='unformatted' )
      read(21) npd
      read(21) aexpnd
      write(*,*) 'density file: np=',npd,' aexpn =',aexpnd
      if ( npd .ne. np1 ) then 
        write(*,*) '# of particles in density and particle files'
        write(*,*) ' do not match:'
        write(*,*) 'npd =',npd
        write(*,*) 'np1 =',np1
        write(*,*) 'aexpnd =',aexpnd
        stop
      endif
      
      read(21) (dnb2(i),i=1,np1)
      close (21)

c      xbox = 2.3507 * 256./25.0 + 1.0 
c      ybox = 2.6745 * 256./25.0 + 1.0 
c      zbox = 13.5273 * 256./25.0 + 1.0 
c      rsel = 2.0* 256 / 25.0 

      ip = 0 
      
      do irow = 1 , Npages         ! loop over particle pages
        In_page = npage
        if ( irow .eq. Npages ) In_page = N_in_last
c         write (*,*)' Read page=',IROW,' file=',ifile,' N=',In_page
        iL = npage * (irow-1)
        CALL GetRow_HF(irow,1) ! read in a page of particles
        do in = 1 , In_page          ! Loop over particles
          ipp = in + iL                     ! current particle number
c          rd = sqrt((xpar(in)-xbox)**2 +
c     &              (ypar(in)-ybox)**2 +
c     &              (zpar(in)-zbox)**2 )
c          if ( rd .le. rsel .and. pw2(ipp) .le. wspecies(1)) then 
          if ( ipp .le. np1 ) then 
            ip = ip + 1
            x(ip) = xpar(in)
            y(ip) = ypar(in)
            z(ip) = zpar(in)
            vx(ip) = vxx(in)
            vy(ip) = vyy(in)
            vz(ip) = vzz(in)
            pw2(ip) = pw2(ipp)
            dnb2(ip) = dnb2(ipp)
          endif
        enddo
      enddo
      np1 = ip 
      write(*,*) 'will use np1=',np1,' particles...'

      close (1)
      close (3)

      return
      end
c


c     -------------------------------
      subroutine ReadParticles_Binary_old
c     -------------------------------
      include 'a_hfind.h'

      Character  Hd*5,Tail*4,Nm*40
      Hd  ='C3crs'
      Tail='.DAT'
C                                     Open file on a tape

c      OPEN(16,FILE='Result.DAT',STATUS='UNKNOWN',POSITION='APPEND')
      OPEN(UNIT=9,FILE='../DAT/CDM2.5/C3CRDa0.999.DAT',
     +                FORM='UNFORMATTED', STATUS='UNKNOWN')
C                                 Read control information
C                                 and see whether it has proper
C                                 structure
      READ  (9) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
      WRITE (*,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble
c      WRITE (16,100) HEADER,
c     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
c     +                  EKIN,EKIN1,EKIN2,
c     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble
100   FORMAT(1X,'Header=>',A45,/
     +            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3)
      IF(NROW.NE.NROWC) THEN
         WRITE (*,*)
     +            ' NROW in PARAMETER and in INPUT file are different'
c         STOP
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (*,*)
     +           ' NGRID in PARAMETER and in INPUT file are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      ENDIF
C                              Open working file on disk
      NBYTE = NRECL*4
c
c.... for DEC alpha nacces=nrecl
c
      nacces = nrecl
c      NACCES= NBYTE

         OPEN(UNIT=1,FILE='../DAT/CDM2.5/C3crs0a0.999.DAT',
     &        ACCESS='DIRECT', STATUS='UNKNOWN',RECL=NACCES)
 
      REWIND 9
      ip   = 0
      Do IROW =1,NROW
         CALL GETROW_HF(IROW,-19)
	 Do  IN=1,NPAGE
            ip      = ip+1
	    x(ip)  = XPAR(IN)
	    y(ip)  = YPAR(IN)
	    z(ip)  = ZPAR(IN)
	    vx(ip) = VXX(IN)
	    vy(ip) = VYY(IN)
	    vz(ip) = VZZ(IN)
         EndDo
      EndDo
      RETURN
      END

C*********************************************************************
C                       Read current PAGE from disk
C
C                       NRECL - length of ROW block in words
 
      SUBROUTINE GETROW_HF(IROW,Ifile)
      INCLUDE 'a_hfind.h'
         READ  (Ifile,REC=IROW) RECDAT
      RETURN
      END


c     -----------------------------
      real function rinput ( text )
c     -----------------------------
      character text*(*)
       
      write (*,'(A,$)') text
      read (*,*) x
      rinput = x

      return
      end

c     -------------------------------------
      character*40 function tinput ( text )
c     -------------------------------------
      character text*(*)
      character*40 textline
       
      write (*,'(A,$)') text
      read (*, 90 ) textline
      tinput = textline

 90   format (A)
      return
      end

c     -----------------------------
      SUBROUTINE indexx(n,arr,indx)
c     -----------------------------
c     (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
c     ----------------------------------------------------------
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
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

      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
      INTEGER mwt,ndata
      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
CU    USES gammq
      INTEGER i
      REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      sx=0.
      sy=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
        ss=0.
        do 11 i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.
      if(mwt.eq.0) then
        do 15 i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        q=1.
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do 16 i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
16      continue
        q=gammq(0.5*(ndata-2),0.5*chi2)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      FUNCTION gammq(a,x)
      REAL a,gammq,x
CU    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.

