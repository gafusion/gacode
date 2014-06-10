    MODULE Nfreya_routines 
      USE nrtype,                                       ONLY : DP,SP,I4B,I2B
                  
      USE neutral_beams,                                ONLY : atwb => atw_beam,bntot , &
                                                               no_injectors,ebeam,ebkev,&
                                                               bion,vbeam,bneut,        &
                                                               n_izpt,time_now,         &
                                                               time_dep_beam,bcur,fbcur,&
                                                               neg_ion_source,nap,nbdep,&
                                                               nfreya_plot_unit,        &
                                                               vx_izpt,vy_izpt,vz_izpt, &
                                                               pitch_a,x_izpt, y_izpt,  &
                                                               r_izpt,z_izpt,iborb,     &
                                                               pbeam,fap,fwall,forb,    &
                                                               hibr,hdep,angmpf,        &
                                                               npart_pseudo,            &
                                                               npart_all_beamlines,     &
                                                               no_physical_injectors,   &
                                                               nsample_izpt,            &
                                                               fb00,fb01,fb10,fb11,     &
                                                               wb00,wb01,wb10,wb11,fber

      USE common_constants,                             ONLY : zeroc,PI,Rad_per_Deg,    &
                                                               kevperg,Proton_Mass,     &
                                                               izero,im32icm3,m2cm

      USE ions_gcnmp,                                   ONLY : nion

      USE MPI_data,                                     ONLY : myid,master,numprocs

      USE error_handler

      USE io_gcnmp,                                     ONLY : ncrt,nlog

      USE ran_num_gen,                                  ONLY : random12,ranorm

      USE zonal_data,                                   ONLY : psivol,nw,nh,nwh,nh2,    &
                                                               n_zone_cntr,r_zone_cntr, &
                                                               z_zone_cntr,mf


      CONTAINS 

 
        SUBROUTINE freyas (mb, injctr_id, codeid, psi, r,                             &
                        rin, rmax, z, zax, zmin, zmax, iexcit,                        &
                        time,time0)
! ----------------------------------------------------------------------
!
!  this subroutine calculates the particle and energy sources
!     due to neutral beam injection.
!
!  the input quantities are:
!
!  mb              Number of neutral beam injectors
!  anglev(ib) (through cangv,sangv)      Vertical angle (degrees) between optical axis
!                    and horizontal plane; a positive value indicates
!                    particles move upward
!  angleh(ib)  (through cangh,sangh)    Horizontal angle (degrees) between optical axis and
!                    vertical plane passing through pivot point and
!                    toroidal axis; a zero value denotes perpendicular
!                    injection, while a positive value indicates par-
!                    ticles move in the co-current direction
!  nsourc          Number of sources per beamline.
!                    If 1, source is centered on beamline axis.
!                    If nsourc = 2, distinguish between the beamline
!                     axis and the source centerline (optical axis).
!                     The two sources are assumed to be mirror images
!                     through the beamline axis.
!                    In either case, the exit grid plane is perpendicula
!                    to the beamline axis, and contains the source
!                    exit grid center(s).
!                    If nsourc = 2, the alignment of the sources w.r.t.
!                    the beamline axis is specified through bhofset,
!                    bvofset, and bleni (described further below).
!  bvofset(ib)     Vertical offset from beamline axis to center
!                    of each source (cm; used only for nsourc = 2)
!  bhofset(ib)     Horizontal offset from beamline axis to center
!                    of each source (cm; used only for nsourc = 2)
!  bleni(ib)       Length along source centerline (source optical axis)
!                    source to intersection point with the beamline axis
!  bcur(ib)        Total current (a) in ion beam (used only if bptor
!                    is zero)
!  bptor(ib)       Total power (w) through aperture into torus; when
!                    nonzero, bptor takes precedence over bcur
!  nbshape(ib)      Beam shape
!                    'circ': circular
!                    'rect': rectangular
!                    'rect-lps'  HSJ
!  bheigh(ib)      Height of source (cm)
!  bwidth(ib)      Width of source (cm); diameter for
!                     circular source.
!  bhfoc(ib)       Horizontal focal length of source (cm)
!  bvfoc(ib)       Vertical focal length of source (cm)
!  bhdiv(ib)       Horizontal divergence of source (degrees)
!  bvdiv(ib)       Vertical divergence of source (degrees)
!  ebkev(ib)       Maximum particle energy in source (keV)
!  fbcur(ie,ib)    Fraction of current at energy ebkeV/ie
!  iborb           Flag for modeling orbit effects on beam-
!                    generated fast ions
!                    3, use mcgo for fast ion slowing down,initial
!                       conditions for the fast ions are generated
!                       using freya.
!                    2, model prompt loss with STAMBAUGH model
!                    1, model orbit effects
!                    0, do not model orbit effects
!  npart           Number of particles followed into plasma
!                    (suggest 10000)
!  npskip          Ratio of number of particles followed into plasma
!                    to number of source particles (suggest 1)
!  naptr           Total number of apertures encountered by a particle
!                    as is moves from the source into the plasma chamber
!                    Maximum is specified by parameter nap ( = 10).
!                    First set of apertures encountered by the particle
!                    are assumed centered on the source axis, and subsequent
!                    apertures are centered on the beamline axis; the
!                    distinction is made through ashape.
!  nashape(iap,ib)  Aperture shape.
!                   Prefix 's-' indicates source axis centered.
!                   Prefix 'b-' indicates beamline axis centered.
!                     's-circ'          'b-circ'
!                     's-rect'          'b-rect'
!                     's-vert'          'b-vert'
!                     's-horiz'         'b-horiz'
!                                       'b-d3d'
!                   (circ = circular aperture, rect=rectagular,
!                    vert = limits vertical height of source particles,
!                    horiz = limits horizontal height of source particles,
!                    d3d= special DIII-D polygonal aperture)
!                    'circ': circular
!                    'rect': rectangular
!  aheigh(iap,ib)  Height of aperture (cm)
!  awidth(iap,ib)  Width  of aperture (cm); diameter for circular aperture
!  alen(iap,ib)    Length from source to aperture for 's-type' aperatur
!                    and from exit grid plane along beamline axis for
!                    'b-type' apertures.
!  blenp(ib)       Distance along beamline axis from source exit
!                    plane to the fiducial "pivot" point.
!  rpivot(ib)      Radial position of pivot (cm)
!  zpivot(ib)      Axial position of pivot (cm)
!     sfrac1(ib)      fraction of source current per beamline coming
!                     from upper source if nsourc =1 sfrac1 is not used
!     iexcit          hexnb switch
!                     0, (default) do not account for atomic excitation
!                        in neutral stopping calculation
!                    >0, use hexnb to account for atomic excitation
!     inubpat        switch controlling 2-d beam deposition calculation
!                    0, (default) bypass 2-d deposition
!                    1, determine beam deposition on (r,z) grid, defined
!                       by npat, and output neutral beam deposition
!                       information to file 'beamdep'
!     npat           dimensions of (r,z) grid for neutral beam
!                    deposition calculation.  used if inubpat = 1.0
!                    npat(1) = number of new 'r' elements
!                    npat(2) = number of new 'z' elements
!                    default:  npat(1) = mi, npat(2) = mj
!     maxp           max number of pseudo particle birth points to be saved.
!                    determined by size of vx_izpt,etc.
!     pbeam(ie,ib)   non normalized power(see sub beam_prop).
!                    normalization to bptor is done here
!
!
!     atw(k)          atomic mass of ion species k
!     codeid          flux surface geometry
!                       "onedee"     : elliptical
!                       anything else: nonelliptical
!     elong           elongation (height/width) of elliptical
!                       cross-section plasma
!     ibion           ion index of beam species
!     mf              number of flux zones
!     mi              number of radial mesh points for nonelliptical
!                       plasma
!     mj              number of axial mesh points for nonelliptical
!                       plasma
!     nion            number of ion species
!
!  b1ins(i)        ratio of R*B to Rax*Bax along horizontal
!                    chord inside and through the magnetic axis vs.
!                    a uniform mesh (rinsid(i)) in major radius;
!                    b1ins is needed only if iborb .gt. 0
!  b1ots(i)        ratio of R*B to Rax*Bax along horizontal
!                    chord outside and through the magnetic axis vs.
!                    a uniform mesh (rotsid(i)) in major radius;
!                    b1ots is needed only if iborb = 1
!  b2ins(i)        ratio of B to Btor along horizontal
!                    chord inside and through the magnetic axis vs.
!                    a uniform mesh (rinsid(i)) in major radius;
!                    b2ins is needed only if iborb = 1
!  b2ots(i)        ratio of B to Btor along horizontal
!                    chord outside and through the magnetic axis vs.
!                    a uniform mesh (rotsid(i)) in major radius;
!                    b2ots is needed only if iborb = 1
!     pinsid(i)       poloidal magnetic flux (G-cm2) along horizontal
!                       chord inside and through the magnetic axis vs.
!                       a uniform mesh (rinsid(i)) in major radius;
!                       pinsid is needed only if iborb .gt. 0
!     potsid(i)       poloidal magnetic flux (G-cm2) along horizontal
!                       chord outside and through the magnetic axis vs.
!                       a uniform mesh (rotsid(i)) in major radius;
!                       potsid(1) and potsid(mf) are needed if
!                       codeid .ne. 'onedee'; the entire potsid array is
!                       needed if iborb = 1
!     psi(i,j)        poloidal magnetic flux (G-cm2) at mesh point i,j
!                       (needed only if codid .ne. 'onedee')
!     csgn            sign of the scalar product of plasma current and
!                       toroidal unit vector phi, j dot phi, for the
!                       right-handed triad (R,phi,Z) where Z points
!                       upward in the lab (DIII-D) frame.
!     psivol(i)       volume (cm3) of flux zone i; depending upon
!                       whether codeid = 'onedee' or not, psivol is
!                       chosen such that either r or SQRT (psi-psiax)
!                       varies a constant amount from one flux surface
!                       to the next
!     rinsid(i)       major radius (cm) of uniform mesh along horizontal
!                       chord inside and through the magnetic axis;
!                       rinsid is needed only if iborb = 1
!     rotsid(i)       major radius (cm) of uniform mesh along horizontal
!                       chord outside and through the magnetic axis;
!                       rotsid(1) and rotsid(mf) are needed if
!                       codeid = 'onedee'; the entire rotsid array is
!                       needed if iborb = 1
!     r(i)            radial mesh point i for nonelliptical plasma (cm)
!     rin             major radius of inside of vacuum vessel (cm)
!                     for 2D equilibria this number must coincide with
!                     rmhdgrid(1).
!     rmax            maximum radial position of plasma (cm)
!     z(j)            axial mesh point j for nonelliptical plasma (cm)
!     zax             axial position of magnetic axis (cm) for
!                       elliptical plasma
!     zmin            minimum axial position of plasma (cm)
!     zmax            maximum axial position of plasma (cm)
!     zne(i)          electron density in zone i (cm-3)
!     zni(i,k)        density of ion species k in zone i (cm-3)
!     zte(i)          electron temperature in zone i (keV)
!     zzi(i,k)        charge number of ion species k in zone i
!  the output quantities are:
!
!     bion(ie,ib)     intensity of ion beam (particles/s)
!     bneut(ie,ib)    intensity of neutral beam (particles/s)
!     pbeam(ie,ib)     beam power to aperture (w)
!     fap(ie,ib)      fraction of beam stopped by aperture
!     fwall(ie,ib)    fraction of beam incident on wall (shinethrough)
!     forb(ie,ib)     fraction of beam lost on orbits
!     ftrapfi(i,ie,ib)fraction of trapped fast ions in each zone
!     fb11(ie,ib)     fraction of ions passing and axis-encircling
!     fb10(ie,ib)     fraction of ions passing and not encircling
!     fb01(ie,ib)     fraction of ions trapped and axis-encircling
!     fb00(ie,ib)     fraction of ions trapped and not encircling
!     fber(ie,ib)     fraction of ions trapped for which error was detected
!     hibrz(i,ie,ib)  normalized hot ion birth rate
!     hdepz(i,ie,ib)  normalized hot ion deposition rate
!     hicmz(i,ie,ib,imd) hot ion creation mode (i.e. ionization or cx)
!     x_izpt(ii,ie,jb)        x coordinate of birth point
!     y_izpt(ii,ie,jb)        y coordinate of birth point
!     z_izpt(ii,ie,jb)        z coordinate of birth point
!     vx_izpt(ii,ie,jb)          x component of birth velocity
!     vy_izpt(ii,ie,jb)          y component of birth velocity
!     vz_izpt(ii,ie,jb)          z component of birth velocity
!     pitch_a(ii,ie,jb)          pitch angle of ion         
!     nsample_izpt(ie,jb)  # points in above arrays 
!     nmbrz(1:mfmi,jb)   counts ions born in zone i, assumes minimum of 1 at least
!
!     wb11(ie,ib)     orbit width of fb11 ions (cm)
!     wb10(ie,ib)     orbit width of fb10 ions (cm)
!     wb01(ie,ib)     orbit width of fb01 ions (cm)
!     wb00(ie,ib)     orbit width of fb00 ions (cm)
!     zetaz(i,ie,ib)  average pitch angle cosine of deposited hot ions
!     angmpz(i,ie,ib) average toroidal angular momentum of a single ion
!                     born in zone i.  no orbit smearing is done.  for a
!                     pencil beam, angmpz will be constant across all
!                     zones.  the orbit smearing of the toroidal
!                     momentum is accomplished (in postnub) by
!                     multiplying angmpf by the smeared deposition rate
!                     of fast ions. units are g cm2/sec
! ----------------------------------------------------------------------
!
  ! moved from argument list to module input:
  ! elong,mi,mj 



      USE nf_param,                                       ONLY : kcmp1,kbe,ke,kb



      USE neutral_beams,                               ONLY :  drpat,dzpat,nrpat,nzpat,    &
                                                               norb,idebug,mi,mj,          &
                                                               mim1,mjm1,elong,            &
                                                               vbeam,cangv ,cangh ,sangv,  &
                                                               sangh ,thetp , npskip,      &
                                                               thetpp ,costp ,sintp ,      &
                                                               costpp ,sintpp,nsourc,      &
                                                               iatype,kt,npulse,pbeamOn,   &
                                                               pbeamOff,beam_on,            &
                                                               beam_end,source2_phase,     &
                                                               nap,nashape,bvofset,        &
                                                               bhofset,bleni,              &
                                                               ftrapfi,hibrz,hdepz,angmpz, &
                                                               zetaz,hicmz,olossc,nbshape, &
                                                               bheigh,bwidth,bhfoc,        &
                                                               bvfoc,bhdiv,bvdiv,sfrac1,   &
                                                               naptr,aheigh,awidth,alen,   &
                                                               blenp,rpivot,zpivot,de_tk,  &
                                                               ne_tk,ds_tk,inubpat,        &
                                                               nmbrz,bptor,nfreya_vb
                                                               
      USE Plasma_properties ,                           ONLY : mhd_dat

      USE zonal_data,                                   ONLY : mf,mfm1,rotsid,rinsid,      &
                                                               potsid,pinsid,              &
                                                               zangrot,b1ins,b1ots,        &
                                                               b2ins,b2ots

      USE xsct,                                         ONLY : sgxn,sgxnloc,sgxnmi,hxfrac

      IMPLICIT NONE

! argument list input: 
      REAL(DP)   psi(nw,nh),r(nw),z(nh)
      REAL(DP)   rin,rmax,zax,zmin,zmax,time,time0
      INTEGER(I4B) iexcit,mb,injctr_id
      CHARACTER(*) codeid

! local vars:   
!      INTEGER,PARAMETER ::  nxx2 = nw*2, nyx2 = nh*2
      REAL(DP)   one,one_millionth,volume,rmajor,elongi,dr,dz,        &
                 dri,dzi
      REAL(DP)   fin(mf),fot(mf),zin(mf),zot(mf),wt(mfm1),zetorb(mfm1)
      REAL(DP)   drin,drot,drutp,drini,droti,drutpi
!     REAL(DP)  rzpat(nxx2,nyx2,ke,kb)
!      REAL(DP)   rzpat(2*nw,2*nh,ke,kb)
      REAL(DP), ALLOCATABLE,DIMENSION(:,:,:,:) ::  rzpat
      REAL(DP)   x0,y0,z0,vx0,vy0,vz0,pzone,rzone,rpos,xpos,ypos,zpos,tenter,&
                 texit,smax,vplane,zetai,vtroid,vrad,angmtm,wid,xnorm,xloss1,&
                 xloss2,bptorx
      INTEGER(I4B) i,ib,ie,iskip,ic,insrc,j,n_izpt,npls,nparx,newpar,ipar,npar,  &
                   mlost,nout,izone,isourc,csgn,ipass,iaxis,ier,izp,nwx2,      &
                   nhx2,nwh,ib_start,ib_end,maxp



!
! ----------------------------------------------------------------------
! general FREYA initialization
! ----------------------------------------------------------------------
!

      maxp          = SIZE(vx_izpt,1) ! also for other *_izpt arrays
      npskip        = 1
      one           = 1.0_DP
      one_millionth = 1.0e-6_DP
      nout = ncrt
      csgn = 1._DP
      nwx2 = 2*nw ; nhx2 = 2*nh ; nwh = nw*nh
      IF(mhd_dat%tot_cur .LT. zeroc) csgn = -1._DP
             


!  set up some data for optional (r,z) deposition calculation
!
      IF (nfreya_vb  .AND. myid .EQ. master) &
             WRITE (ncrt, '('' starting FREYA,time = '',1pe18.8)')time

!  zero orbit output parameters
!  if output is desired, these need to be set accordingly
!  NOTE norb .ne. 0 causes printout in orbit_12. This
!  will need attention if it is activate due to
!  multiple cpu use. HSJ
      norb = 0 ! otherwise set to iounit to be written to.
      idebug(:) = zeroc



!
!     initialize flux surface quantities
!
!     mfm1   = mf-1
      rmajor = rotsid(1)
      drin   = (rinsid(1)-rinsid(mf))/mfm1
      drot   = (rotsid(mf)-rotsid(1))/mfm1
      drutp  = SQRT (potsid(mf)-potsid(1))/mfm1
      drin   =  MAX (drin , one_millionth)
      drot   =  MAX (drot , one_millionth)
      drutp  =  MAX (drutp, one_millionth)
      elong  =  MAX (elong, one_millionth)
      drini  = 1.0_DP / drin
      droti  = 1.0_DP / drot
      drutpi = 1.0_DP / drutp
      elongi = 1.0_DP / elong
!

      mim1  = mi - 1
      mjm1  = mj - 1
      dr    = (r(mi)-r(1))/mim1
      dz    = (z(mj)-z(1))/mjm1
      dri   = 1.0_DP / dr
      dzi   = 1.0_DP / dz


!     calculate total plasma volume
!
      volume = zeroc
      DO 30 i=1,mfm1
   30 volume = volume + psivol(i) ! cm^3


!
! ----------------------------------------------------------------------
!     begin loop over beams
! ----------------------------------------------------------------------
!
      iskip = 1 + (npart_all_beamlines-1)/maxp
      ic    = izero



      IF(injctr_id == izero)Then
         ib_start = 1
         ib_end = mb
      ELSE
         ib_start = injctr_id
         ib_end   = injctr_id
      ENDIF

      DO 200 ib = ib_start,ib_end     !  note assumes do loop is executed at least once
                                      !  P_Nfreya uses time_dep_beam = -1, set in sub nfreya_init

         IF(time_dep_beam .EQ. 1) THEN 
!           note that we assume that all sources of beam ib are on
!           or off together. (ie source2_phase .ne. 0.0  is not
!           implemented at this time)
!           The logic in sub source dictated that Freya is called.
!           We now have to decide what beam(s) are actually on.
!           Skip those beams that are currently off:
            IF(nfreya_vb)PRINT *,'beam_on , source2_p,ib', &
                           beam_on(ib),source2_phase(ib),ib
            IF(time .LT. beam_on(ib)+source2_phase(ib))   go to 200
!            if(time .gt. beam_end(ib)+source2_phase(ib)) go to 200
!              skip beamlines/sources  that aren't currently on:
!               do nsrc = 1, nsourc
!                   insrc = nsrc
                   insrc = 1              ! source phase not implemented
                   DO npls =1,kt          ! parameter (kt = 200) 
                      npulse(ib) = npls
                      IF( time .GE. pbeamOn(npls,insrc,ib) .AND. &
                        time .LE. pbeamOff(npls,insrc,ib)) go to 202
                   ENDDO
!               enddo
               npulse(ib) = 0
               IF(nfreya_vb)PRINT *,'skipping beam #',ib,' at time',time
               go to 200  !source is not currently on,skip this beam

           ENDIF



202        CONTINUE



!
!       Determine aperture designators
!
 

        DO i=1,nap ! nap is parameter

           ! old, pre April 2011 
           !          IF (nashape(i,ib) .EQ. 's-circ' ) iatype(i,ib) = 1
           !          IF (nashape(i,ib) .EQ. 's-rect' ) iatype(i,ib) = 2
           !          IF (nashape(i,ib) .EQ. 's-vert' ) iatype(i,ib) = 3
           !          IF (nashape(i,ib) .EQ. 's-horiz') iatype(i,ib) = 4
           !          IF (nashape(i,ib) .EQ. 'b-circ' ) iatype(i,ib) = 5
           !          IF (nashape(i,ib) .EQ. 'b-rect' ) iatype(i,ib) = 6
           !          IF (nashape(i,ib) .EQ. 'b-vert' ) iatype(i,ib) = 7
           !          IF (nashape(i,ib) .EQ. 'b-horiz') iatype(i,ib) = 8
           !          IF (nashape(i,ib) .EQ. 'b-d3d'  ) iatype(i,ib) = 9
           !          END DO

           ! new after april 2011 (onetwo v5.4)( corrections made by Bob Harvey) HSJ
           !BH110315          IF (nashape(i,ib) .EQ. 's-circ' ) iatype(i,ib) = 1
           !BH110315          IF (nashape(i,ib) .EQ. 's-rect' ) iatype(i,ib) = 2
           !BH110315          IF (nashape(i,ib) .EQ. 's-vert' ) iatype(i,ib) = 3
           !BH110315          IF (nashape(i,ib) .EQ. 's-horiz') iatype(i,ib) = 4
          IF (nsourc.eq.1) then
             if(nashape(i,ib).eq.'s-circ')  iatype(i,ib)=5 ! changed ashape to nashape HSJ 4/1/2011
             if(nashape(i,ib).eq.'s-rect')  iatype(i,ib)=6
             if(nashape(i,ib).eq.'s-vert')  iatype(i,ib)=7
             if(nashape(i,ib).eq.'s-horiz') iatype(i,ib)=8
          ELSE IF (nsourc.gt.1) then   !     iatype .le. 4 are source-centered apertures HSJ
             if(nashape(i,ib).eq.'s-circ')  iatype(i,ib)=1
             if(nashape(i,ib).eq.'s-rect')  iatype(i,ib)=2
             if(nashape(i,ib).eq.'s-vert')  iatype(i,ib)=3
             if(nashape(i,ib).eq.'s-horiz') iatype(i,ib)=4
          ENDIF                        !     iatype .ge. 5 are beamline-centered apertures HSJ
          IF (nashape(i,ib) .EQ. 'b-circ' ) iatype(i,ib) = 5
          IF (nashape(i,ib) .EQ. 'b-rect' ) iatype(i,ib) = 6
          IF (nashape(i,ib) .EQ. 'b-vert' ) iatype(i,ib) = 7
          IF (nashape(i,ib) .EQ. 'b-horiz') iatype(i,ib) = 8
          IF (nashape(i,ib) .EQ. 'b-d3d'  ) iatype(i,ib) = 9

            
        END DO



!
!       Some angles for subroutine ROTATE

        thetp(ib)  = &
            ATAN2 (bvofset(ib), SQRT (bleni(ib)**2-bvofset(ib)**2))
        costp(ib)  = COS (thetp(ib))
            sintp(ib)  = SIN (thetp(ib))
        thetpp(ib) = &
            ATAN2 (bhofset(ib), SQRT (bleni(ib)**2-bvofset(ib)**2))
        costpp(ib) = COS (thetpp(ib))
        sintpp(ib) = SIN (thetpp(ib))




!
! ----------------------------------------------------------------------
! begin loop over beam energy components
! ----------------------------------------------------------------------

        nmbrz(:,ib) = izero
               ! --- nmbrz(i,ib) counts ions born in zone i,injector ib
               ! --- it accumulates all three energy componenets
               ! --- minimum of 1 is forced below

        DO 201 ie=1,3
           n_izpt       = izero
           nsample_izpt(ie,ib) = izero
           fap(ie,ib)   = zeroc
           fwall(ie,ib) = zeroc
           forb(ie,ib)  = zeroc
           fb11(ie,ib)  = zeroc
           fb10(ie,ib)  = zeroc
           fb01(ie,ib)  = zeroc
           fb00(ie,ib)  = zeroc
           fber(ie,ib)  = zeroc
           wb11(ie,ib)  = zeroc
           wb10(ie,ib)  = zeroc
           wb01(ie,ib)  = zeroc
           wb00(ie,ib)  = zeroc

           DO  i=1,mfm1
              ftrapfi(i,ie,ib) = zeroc
              hibrz(i,ie,ib)   = zeroc
              hdepz(i,ie,ib)   = zeroc
              angmpz(i,ie,ib)  = zeroc
              zetaz(i,ie,ib)   = zeroc
              hicmz(i,ie,ib,1) = zeroc
              hicmz(i,ie,ib,2) = zeroc
              hicmz(i,ie,ib,3) = zeroc
              olossc(i,ie,ib)  = zeroc   
           ENDDO


           ! ----------------------------------------------------------------------
           ! begin loop over particles
           ! ----------------------------------------------------------------------
               
               !!npar = (bneut(ie,ib)/bntot)*npart            ! #neutrals to launch for this
                                                              ! energy component,beamline
               npar = npart_pseudo(ie,ib)

               IF (npar .EQ. 0)  go to 201
               nparx  = 0
               newpar = 0

 
 
               DO 180 ipar = 1,npar                        !particle loop starts here

 
                  IF (MOD (ipar-1,npskip) .EQ. 0)  newpar = 1

                  IF (newpar .EQ. 0)  go to 120
                  !
                  !  generate neutral particle at beam source
                  !
                  ! NOTE currently based on beam no need to change to beam
                  ! id of some sort when it becomes available in namelist input
                  IF(injctr_id == 0 .AND. ( ib .GE. 3 .AND. &
                       ib .LE. 10))THEN   ! single cpu case

                      CALL sorspt_150bm(nbshape,bheigh,bwidth,               &
                                        bhfoc,bvfoc,bhdiv,bvdiv,ib,ie,       &
                                        isourc,ke,nsourc,sfrac1,vbeam,x0,    &
                                        y0,z0,vx0,vy0,vz0)

                  ELSE IF( (3 .LE. injctr_id) .AND. (injctr_id .le.10))THEN
                       ! multiple cpu case

                       CALL sorspt_150bm(nbshape,bheigh,bwidth,              &
                                         bhfoc,bvfoc,bhdiv,bvdiv,ib,ie,      &
                                         isourc,ke,nsourc,sfrac1,vbeam,x0,   &
                                         y0,z0,vx0,vy0,vz0)
 
                  ELSE 
                     ! both single and multiple cpus

                     CALL sorspt(nbshape,bheigh,bwidth,bhfoc,bvfoc,bhdiv,bvdiv,ib,ie, &
                       isourc,ke,nsourc,sfrac1,vbeam,x0,y0,z0,vx0,vy0,vz0)
 
                  !
                  !  transform coordinates and advance particle to pivot point
                  !
                  ENDIF
 
                  CALL rotate(naptr,iatype,aheigh,awidth,alen,bhofset,bvofset,cangv, &
                       cangh,ib,isourc,costp,sintp,costpp,sintpp,blenp, &
                       nsourc,sangv,sangh,rpivot,zpivot,mlost,x0,y0,z0,vx0, &
                       vy0,vz0)
 

                  !  skip injection if particle is lost at aperture
                  !
120               IF (mlost .NE. 0)  go to 160

                  !
                  !  inject particle into plasma, i.e., follow particle from pivot
                  !     point into or through the plasma
                  !

 
                  CALL inject( codeid, de_tk,drutpi,droti,dri,ds_tk,dzi,       &
                       elongi,ib,ie,kcmp1,kbe,ke,nw,mfm1,mim1,mjm1,ne_tk,newpar, &
                       nout,potsid(1),psi,r,rmajor,rin,rmax,                        &
                       x0,y0,z0,vx0,vy0,vz0,vbeam,z,zangrot,zax,zmin, &
                       zmax,izone,pzone,rzone,rpos,xpos,ypos,zpos, &
                       tenter,smax,texit)

                  !
                  !  skip birth data if:  particle missed plasma
                  !

                  IF (izone .GE. mf)  go to 170
                  !
                  !  accumulate hot ion birth rate and creation mode
                  !

                  hibrz(izone,ie,ib)   = hibrz(izone,ie,ib) + 1.0
                  hicmz(izone,ie,ib,1) = hicmz(izone,ie,ib,1) + sgxnloc(1)
                  hicmz(izone,ie,ib,2) = hicmz(izone,ie,ib,2) + sgxnloc(2)
                  hicmz(izone,ie,ib,3) = hicmz(izone,ie,ib,3) + sgxnloc(3)

                  !  calculate and accumulate pitch angle at birth point
                  !  calculate angular momentum deposited by each monte carlo ion
                  !
                  !  NOTE AS OF 10/18/94 THE LINES MARKED WITH ! rs WERE ADDED. SOME OF
                  !  THESE CHANGES WILL AFFECT THE NEUTRAL BEAM RESULTS (PRESUMABLY ONLY
                  !  SLIGHTLY) EVEN WHEN THE STAMBAUGH ORBIT LOSS MODEL IS NOT USED. ... HSJ
                  !    should replace zetai this with  full formula HSJ
                  !
                  !***  zetai  = csgn*(xpos*vy0-ypos*vx0)/(rpos*vbeam(ie,ib))   ! ORIGINAL
                  vplane = SQRT (vx0**2+vy0**2)                           ! rs
                  zetai  = csgn*(xpos*vy0-ypos*vx0)/(rpos*vplane)         ! rs
                  zetai  = MIN (zetai,  one)
                  zetai  = MAX (zetai, -one)
                  !***  vtroid = csgn*zetai*vbeam(ie,ib)                        ! ORIGINAL
                  vtroid = csgn*zetai*vplane                              ! rs
                  vrad   = -SQRT (1.0_DP - zetai**2) * vplane             ! rs
                  angmtm = rpos*vtroid*atwb*1.673e-24  ! g cm2/sec
                  angmpz(izone,ie,ib) = angmpz(izone,ie,ib) + angmtm
                  IF (iborb .NE. 1)                                    &  ! rs
                       zetaz(izone,ie,ib) =  zetaz(izone,ie,ib) + zetai
                  !      t_zetaz(izone,ie,ib) =  zetaz(izone,ie,ib)
                  !
                  !  save occasional birth point for subsequent plotting
                  !
                  ic = ic + 1 ! ic counts particles over all beams,energies
                  !IF (MOD (ic-1, iskip) .EQ. 0 ) THEN
                  IF(n_izpt .LT. SIZE(vx_izpt,1))THEN
                     n_izpt       = n_izpt + 1
                     n_izpt=MIN(n_izpt,maxp)            ! dont let the arrays overflow 
                     x_izpt(n_izpt,ie,ib)    = xpos              
                     y_izpt(n_izpt,ie,ib)    = ypos
                     z_izpt(n_izpt,ie,ib)    = zpos
                     r_izpt(n_izpt,ie,ib)    = SQRT(xpos**2 + ypos**2)
                     vx_izpt(n_izpt,ie,ib)   = vx0
                     vy_izpt(n_izpt,ie,ib)   = vy0
                     vz_izpt(n_izpt,ie,ib)   = vz0
                     pitch_a(n_izpt,ie,ib)   = zetai
                  END IF

                  !
                  !  calculate (r,z) grid location
                  !
                  IF (inubpat .GT. 0 ) THEN
                     IF(.NOT. ALLOCATED(rzpat))THEN 
                        ALLOCATE(rzpat(nwx2,nhx2,ke,kb))
                        rzpat(:,:,:,:) = zeroc
                     ENDIF
                     i = (rpos-r(1))/drpat + 1.0
                     j = (zpos-z(1))/dzpat + 1.0
                     rzpat(i,j,ie,ib) = rzpat(i,j,ie,ib) + 1.0_DP
                  END IF
                  !

150               IF (iborb .EQ. 1) THEN                                     
                     nparx  = nparx + 1

                     CALL orbit_12(atwb,b1ins,b1ots,b2ins,b2ots, codeid,ic,idebug, &
                          iskip, &
                          izone,mf,norb,pinsid,potsid,pzone,rinsid,rotsid, &
                          rzone,rpos,vbeam(ie,ib),zetai,zpos,fin,fot,zin,zot, &
                          ipass,iaxis,ier,izp,wid,wt,zetorb)

                     IF (ier .NE. 0) THEN
                        fber(ie,ib) = fber(ie,ib) + 1.0

                     ELSE IF (izp .GT. mfm1) THEN
                        go to 175
                     ELSE
                        nmbrz(izp,ib) = nmbrz(izp,ib) + 1
                        IF      (ipass .EQ. 1 .AND. iaxis .EQ. 1) THEN
                           fb11(ie,ib) = fb11(ie,ib) + 1.0
                           wb11(ie,ib) = wb11(ie,ib) + wid

                        ELSE IF (ipass .EQ. 1 .AND. iaxis .EQ. 0) THEN
                           fb10(ie,ib) = fb10(ie,ib) + 1.0
                           wb10(ie,ib) = wb10(ie,ib) + wid

                        ELSE IF (ipass .EQ. 0 .AND. iaxis .EQ. 1) THEN
                           fb01(ie,ib) = fb01(ie,ib) + 1.0
                           wb01(ie,ib) = wb01(ie,ib) + wid
                           ftrapfi(izp,ie,ib) = ftrapfi(izp,ie,ib) + 1.0

                        ELSE IF (ipass .EQ. 0 .AND. iaxis .EQ. 0) THEN
                           fb00(ie,ib) = fb00(ie,ib) + 1.0
                           wb00(ie,ib) = wb00(ie,ib) + wid
                           ftrapfi(izp,ie,ib) = ftrapfi(izp,ie,ib)+1.0

                        END IF
                     END IF
                     !
                     !  accumulate hot ion deposition rate and average pitch angle cosine
                     !
                     DO i=1,mfm1
                        hdepz(i,ie,ib) = hdepz(i,ie,ib) + wt(i)
                        zetaz(i,ie,ib) = zetaz(i,ie,ib) + wt(i)*zetorb(i)
                     END DO
                  END IF                !iborb=1
                  go to 180
                  !
                  !  accumulate particles that are not deposited in plasma
                  !
160               fap(ie,ib)   = fap(ie,ib) + 1.0
                  go to 180
170               fwall(ie,ib) = fwall(ie,ib) + 1.0
                  go to 180
175               forb(ie,ib)  = forb(ie,ib) + 1.0


                  ! ----------------------------------------------------------------------
                  ! end loop over particles
                  ! ----------------------------------------------------------------------
                 
180               newpar = 0
                  nsample_izpt(ie,ib) = n_izpt


                  !
                  !  normalize average pitch angle cosine.  normalize momentum and birth
                  !     mode in each shell to a single particle.
                  !
                  DO 182 i=1,mfm1
                     xnorm = hibrz(i,ie,ib)
                     IF (xnorm .NE. zeroc) THEN
                        angmpz(i,ie,ib)  = angmpz(i,ie,ib)/xnorm
                        hicmz(i,ie,ib,1) = hicmz(i,ie,ib,1)/xnorm
                        hicmz(i,ie,ib,2) = hicmz(i,ie,ib,2)/xnorm
                        hicmz(i,ie,ib,3) = hicmz(i,ie,ib,3)/xnorm
                     END IF
                     IF (iborb .NE. 0)  xnorm = hdepz(i,ie,ib)
                     IF (xnorm .NE. zeroc) xnorm = 1.0/xnorm
                     !***      angmpz(i,ie,ib) = xnorm*angmpz(i,ie,ib)
182                  zetaz(i,ie,ib) = xnorm*zetaz(i,ie,ib)
                     !
                     !  get the fraction of trapped ions in each zone.  if orbit effects are
                     !     not turned on or if itrapfi = 0 then this effect is not included.
                     !
                     DO i=1,mfm1
                        nmbrz(i,ib) = MAX0 (nmbrz(i,ib),1)
                        ftrapfi(i,ie,ib) = ftrapfi(i,ie,ib)/nmbrz(i,ib)
                     END DO
                     !
                     !  normalize loss fractions and hot ion birth rate
                     !
                     fap(ie,ib)   = fap(ie,ib)/npar
                     fwall(ie,ib) = fwall(ie,ib)/npar
                     forb(ie,ib)  = forb(ie,ib)/npar
                     xloss1       = fap(ie,ib) + fwall(ie,ib)
                     xloss2       = fap(ie,ib) + fwall(ie,ib) + forb(ie,ib)
                     !***  if (xloss1 .ge. 1.0)  go to 200
                     IF (xloss1 .GE. 1.0)  go to 201
                     !
                     DO i=1,mfm1
                        hibrz(i,ie,ib) = &
                             hibrz(i,ie,ib)*volume/((1.0-xloss1)*npar*psivol(i))
                        IF      ( iborb .EQ. 0  ) THEN
                           hdepz(i,ie,ib) = hibrz(i,ie,ib)
                        ELSE IF (xloss2 .LT. 1.0) THEN
                           hdepz(i,ie,ib) = &
                                hdepz(i,ie,ib)*volume/((1.0-xloss2)*npar*psivol(i))
                        END IF

                     END DO
                     !
                     !  normalize orbit widths and fractions
                     !
                     IF (iborb  .EQ. 1) THEN
                        IF (fb11(ie,ib) .NE. zeroc)  wb11(ie,ib) = wb11(ie,ib)/fb11(ie,ib)
                        IF (fb10(ie,ib) .NE. zeroc)  wb10(ie,ib) = wb10(ie,ib)/fb10(ie,ib)
                        IF (fb01(ie,ib) .NE. zeroc)  wb01(ie,ib) = wb01(ie,ib)/fb01(ie,ib)
                        IF (fb00(ie,ib) .NE. zeroc)  wb00(ie,ib) = wb00(ie,ib)/fb00(ie,ib)
                        fb11(ie,ib) = fb11(ie,ib)/nparx
                        fb10(ie,ib) = fb10(ie,ib)/nparx
                        fb01(ie,ib) = fb01(ie,ib)/nparx
                        fb00(ie,ib) = fb00(ie,ib)/nparx
                        fber(ie,ib) = fber(ie,ib)/nparx
                     END IF
                     !
              IF (nfreya_vb  ) &
                  WRITE (ncrt, '("process",i3, " did   injection for beam # ", &
                           i2," energy index  ",i2, " time  = ",1pe14.6,"npar  =",i9)') myid,ib,ie,time,npar

  201          CONTINUE                ! end loop over energy components

 

  200 CONTINUE                ! end loop over beams






!-----------------------------------------------------------------------
! -- renormalize currents and powers to bptor
!-----------------------------------------------------------------------

      DO ib= ib_start,ib_end
        IF (bptor(ib) .GT. zeroc) THEN
          bptorx = zeroc
          DO  ie=1,3
             bptorx = bptorx + (1.0_DP-fap(ie,ib))*pbeam(ie,ib)
          ENDDO
          IF (bptorx .GT. zeroc) THEN
            xnorm = bptor(ib)/bptorx              ! normalize to  bptor (power into torus)
            bcur(ib) = xnorm*bcur(ib)
            DO  ie=1,3                                        
              bion (ie,ib)  = xnorm*bion(ie,ib)
              bneut(ie,ib)  = xnorm*bneut(ie,ib)
              pbeam (ie,ib) = xnorm*pbeam(ie,ib)  ! note pbeam is pwr to apperture
            ENDDO                                 ! bptor is power in torus
          END IF
        END IF
      END DO


!
!     calculate neutral beam density on (r,z) grid
!
      IF (inubpat .EQ. 1)      &
        CALL nbdep2d (psi,nw,nwx2,nh,nhx2, nwh,mi, mj, r, z, potsid,rzpat, nrpat, &
                      nzpat, ke, mb, vbeam, hxfrac)

 !go to 1010
       ie = 3
       DO j=ib_start,ib_end
         WRITE(ncrt,FMT='("beam no",i5,"      beam power into torus ",1pe12.4)')j,bptor(j) 
         WRITE(ncrt,FMT='("beam energy kev ")')
         WRITE(ncrt,1)(ebeam(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("intensity of ion beam #/sec",i5)')
         WRITE(ncrt,1)(bion(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("intensity of neutral beam #/sec ")')
         WRITE(ncrt,1)(bneut(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("beam power to aperature,w")')
         WRITE(ncrt,1)(pbeam(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction stopped by apperature")')
         WRITE(ncrt,1)(fap(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction incident on wall")')
         WRITE(ncrt,1)(fwall(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction of beam lost on orbits")')
         WRITE(ncrt,1)(forb(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction of ions passing and axis-encircling")')
         WRITE(ncrt,1)(fb11(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction of ions passing and not  axis-encircling")')
         WRITE(ncrt,1)(fb10(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction of ions trapped and axis-encircling")')
         WRITE(ncrt,1)(fb01(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction of ions trapped and not encircling")')
         WRITE(ncrt,1)(fb00(i,j),i= 1,ie)
         WRITE(ncrt,FMT='("fraction of ions trapped for which error was detected")')
         WRITE(ncrt,1)(fber(i,j),i= 1,ie)
       ENDDO
1      FORMAT(3(5x,1pe12.4))
1010      IF (nfreya_vb  )  WRITE (ncrt, '('' done with FREYA'')') 
 

      RETURN

      END SUBROUTINE freyas

      SUBROUTINE  beam_prop
!----------------------------------------------------------------------
! calculate beam power; account for neutralizer efficiency
! NOTE: scaling to bptor (namelist input power) is done in freya after
!       the deposition calculations
! INPUT (from neutral_beams module)
!  atw_beam,no_injectors,ebkev,neg_ion_source,fbcur,bcur
!
!
! OUTPUT (to neutral_beams module) see description on freya:
!  atwb
!  bntot
!  ebeam
!  bion
!  bneut
!  pbeam
!  vbeam
!  npart_pseudo(ie,ib) no pseudo particles by energy, and injectors

! LOCAL
! ebx
! beff
! ebev
!  --------------------------------------------------------------------

      USE nf_param,                                     ONLY :  ke,kb

      USE Plasma_properties,                            ONLY : neut_beam

      USE common_constants,                             ONLY : joupkev

      USE nub,                                          ONLY : bfr_neutrlz
     
      USE neutral_beams,                                ONLY : bptor

      IMPLICIT NONE

      REAL(DP)  ebx,beff,ebev,tot_pow_inj,frac,fracp
      INTEGER ib,ie,e_err

      bntot = zeroc ; e_err = izero ;            tot_pow_inj = zeroc

      !---------------------------------------------------------------------------
      ! -- In  P_Nfreya we want to use the input beam current fractions
      ! -- given after the neutralizer (bfr_neutrlz = .FALSE. )
      !-------------------------------------------------------------HSJ-9/28/11---
      DO  ib=1,no_injectors
         DO  ie=1,3
            IF(ebkev(ib) .lt. 1._DP)e_err = e_err + 1 ! expect > 1 KEV
            ebeam(ie,ib) = ebkev(ib)/ie
            ebx          = ebeam(ie,ib)/atwb
            IF (neg_ion_source(ib) .GT. 0) THEN
               beff = 0.98           ! arbitrary efficiency HSJ
            ELSE
               CALL logint (ebx, beff)
            END IF
            ebev          = 1.0e3*ebx
            vbeam(ie,ib)  = 1.384e6_DP * SQRT (ebev)
            IF(bfr_neutrlz)THEN ! standard before neutralizer input of fbcur
                bion(ie,ib)   = 0.625e19_DP*fbcur(ie,ib)*bcur(ib)     ! # singly charged ions/sec
                                                                  ! entering the neutralizer
                bneut(ie,ib)  = ie*beff*bion(ie,ib)
            ELSE !  fbcur is given after neutralizer 
                bion(ie,ib)   = 0.625e19_DP*fbcur(ie,ib)*bcur(ib)/beff     ! # singly charged ions/sec
                                                                  ! entering the neutralizer
                bneut(ie,ib)  = ie*bion(ie,ib)
            ENDIF
!           pbeam(ie,ib)  = ebeam(ie,ib)*bneut(ie,ib)/0.625e16_DP ! this is power to aperture
            pbeam(ie,ib)   = ebeam(ie,ib)*bneut(ie,ib)*joupkev 
            bntot         = bntot + bneut(ie,ib)                  ! renormalization to bptor 
                                                                  ! is done in freyas
         ENDDO
            tot_pow_inj   = tot_pow_inj + bptor(ib)
      ENDDO


      IF(e_err .GT. izero)THEN
         WRITE(ncrt,FMT='("ERROR:, beam energies not set ")')
         lerrno = iomaxerr + 357_I4B
         CALL terminate(lerrno,nlog)
      ENDIF

      ALLOCATE(npart_pseudo(ke,no_injectors))
      DO  ib=1,no_injectors
         DO  ie=1,3
 !           npart_pseudo(ie,ib)   =  (bneut(ie,ib)/bntot) * npart_all_beamlines
             frac                  =  fbcur(ie,ib)/(fbcur(1,ib)+fbcur(2,ib)+fbcur(3,ib))
             fracp                 =  bptor(ib)/tot_pow_inj
             npart_pseudo(ie,ib)   =  npart_all_beamlines*frac*fracp
 !            npart_pseudo(ie,ib)   =  MAX(npart_pseudo(ie,ib),100)
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE  beam_prop
 

      SUBROUTINE getsgxn (e, iz, ns, ktk, ib, ie, nbins, &
                          debin, kcmp1,kbe,nout,sgxnout)
! ----------------------------------------------------------------------
! this subroutine locates n*sigma in the lookup table, sgxn
! ----------------------------------------------------------------------
!
! --- input through argument list:
!          debin          - energy bin width     (keV/amu)
!          e(ns)          - energy to be located (keV/amu)
!          ib             - beam line index
!          ie             - beam energy group index
!          iz(ns)         - flux zone index
!          nbins          - total number of energy bins
!          ns             - total number of points to be evaluated
!          sgxn(i,iz,j,k) - array of neutral stopping data (module xsct)
!                           i  - data type index
!                                =1, fraction of interactions producing
!                                    electrons;
!                                =2, fraction of interactions producing
!                                    neutral species 1;
!                                =3, fraction of interactions producing
!                                    neutral species 2;
!                                =4, inverse mean free path (cm**-1)
!                           iz - flux zone index
!                           j  - beam/energy index, j = 3*(ib-1)+ie
!                           k  - energy bin index
!  NOTE this is actually sgxntab(1:ktk),not sgxntab(i) HSJ

! --- output through argument list:
!          if ns == 1 then 
!            sgxnout(1:4) is returned
!          if ns > 1 then
!           sgxnout(1:ns) is returned loaded with index i=4
!
!          sgxnout(i)     - local neutral stopping data (cm**-1)
!                           i - data type index (see above)
!                           note:  iftns>1, only the inverse mean free
!                           path is evaluated (i = 4) and ns data points
!                           are passed.  this facilitates faster
!                           execution when evaluating the max imfp along
!                           the collisionless neutral trajectory (used
!                           in rotating discharges only).  otherwise
!                           more detailed information is passed in the
!                           first four elements only (ns = 1).
! ----------------------------------------------------------------------
!
!   Created:    Glenn Sager?               Date:  ??-???-????
!
!   changes from original version:
!      1) The code will now terminate if the maximum rotational energy
!         bin is exceeded. Detailed Error message provide to both the
!         OUTONE file (nout) and the standard output device (ncrt = 6).
!      2) Fixed some math in array index calculations.
!      3) Today is the 50th Aniversary of Hiroshima. NEVER AGAIN!
!
!                                  D.F. Finkenthal 6-AUG-95
!
! ----------------------------------------------------------------------
!!
      USE xsct,                                         ONLY : sgxn
!      USE MPI_data,                   ONLY : myid ! 888889999
      IMPLICIT NONE
!
      REAL(DP)   e(*),sgxnout(*)
      REAL(DP) debin,emax
      INTEGER(I4B),SAVE :: imax ,istart 
      INTEGER(I4B) i,imaxa(200),ncrt,ind,ns,ibin,ktk,ib,ie,nbins, &
                   kcmp1,kbe,nout, iz(*)
!
      DATA istart /0/
      ncrt = 6
      IF(istart == 0)imaxa(:) = 0
!
      ind = 3 * (ib - 1) + ie
      imax = 0 !!hsj
      IF (ns .GT. 1) THEN
        DO i=1,ns  ! get ns values of sgxnout wiht type 4 index
          ibin       = e(i) / debin + 1.0
          imaxa(i)   = ibin
          ibin       = MIN0 (ibin, nbins)
          sgxnout(i) = sgxn (4, iz(i), ind, ibin)
          imax = MAX(imax,imaxa(i)) !! HSJ
        END DO
!        imax = MAXVAL (imaxa) !! HSJ
      ELSE

        ibin       = e(1) / debin + 1.0
        imax       = MAX0 (ibin, imax )
        ibin       = MIN0 (ibin, nbins)
        sgxnout(1) = sgxn (1, iz(1), ind, ibin)
        sgxnout(2) = sgxn (2, iz(1), ind, ibin)
        sgxnout(3) = sgxn (3, iz(1), ind, ibin)
        sgxnout(4) = sgxn (4, iz(1), ind, ibin)
      END IF
!



      IF (imax .GT. nbins) THEN
!        emax = amaxaf (e, 1, ns)
         emax = e(1)
         DO i=2,ns
            emax = MAX(emax,e(i))
         ENDDO
        WRITE (nout, 110)  nbins * debin, emax,myid
        WRITE (ncrt, 110)  nbins * debin, emax,myid
!        WRITE(940+myid,110)nbins * debin, emax,myid ! 8888889999
!        CALL STOP ('subroutine GETSGXN: max bin energy exceeded', 51)
          lerrno = iomaxerr + 228_I4B
         CALL terminate(lerrno,nlog)
      END IF
!
 110  FORMAT (/ &
       ' ERROR in GETSGXN: maximum rotational energy bin exceeded'     / &
       ' The highest calculated energy bin is now ', f10.2, ' keV/amu' / &
       ' but GETSGXN tried to look up a value of  ', f10.2, ' keV/amu' / &
       ' The FE_TK input parameter must be increased to overcome', &
       ' this problem!'                                                / &
       ' ONETWO will be terminated to avoid further problems.',i5)
!
      istart = istart + 1

      RETURN
!
      END SUBROUTINE getsgxn




      FUNCTION imax (f, mfm1)
! ----------------------------------------------------------------------
!   This function computes the index of the maximum element of the vector
!   f excluding the first and last elements.
! ----------------------------------------------------------------------

      REAL(DP)  f(*),fmax
      INTEGER(I4B) imax,mfm1,i
      fmax = f(2)
      imax = 2
      DO i=3,mfm1
        IF (f(i) .GE. fmax) THEN
          fmax = f(i)
          imax = i
        END IF
      END DO
      RETURN

      END FUNCTION imax


      SUBROUTINE inject (codeid, debin,drutpi,droti,dri,ds1,dzi, &
                         elongi,ib,ie,kcmp1,kbe,ke,ki,mfm1,mim1,mjm1,nebin, &
                         newpar,nout,psiax,psi,r,rmajor,rin,rmax, &
                         x0,y0,z0,vx0,vy0,vz0,vbeam,z, &
                         zangrot,zax,zmin,zmax,izone,pzone,rzone,rpos, &
                         xpos,ypos,zpos,tenter,smax,texit)
! ----------------------------------------------------------------------
! --  this subroutine follows the particle from the pivot point into,
! --  through, or around the plasma.  note:  bilinear interpolation is
! --  used to obtain psi as a function of track length along a
! --  collisionless neutral trajectory (rotating discharges only).
! --  bicubic spline interpolation was tested, but found to provide no
! --  appreciable increase in accuracy.
! ----------------------------------------------------------------------
!
!
      USE xsct,                                         ONLY : sgxn,sgxnloc,sgxnmi

      IMPLICIT NONE


      REAL(DP)  cvec(200)

!
      INTEGER(I4B),PARAMETER :: ktk = 100
      CHARACTER(*)    codeid
      INTEGER(I4B) ib,ie,kcmp1,kbe,ke,ki,mfm1,mim1,mjm1,nebin,newpar
      REAL(DP)     psi(ki,*),                                            &
                   r(*), vbeam(ke,*), z(*), zangrot(*), e1(ktk),         &
                   sgxntab(ktk),debin,drutpi,droti,ds1,dri,delt,dzi

      REAL(DP)     x0,y0,z0,vx0,vy0,vz0,zmin,zmax,tenter,texit,          &
                   x1,y1,z1,r1,rmajor,zax,elongi,psiax,rin,rmax,pzone,   &
                   rzone,xmassp,xpos,ypos,zpos,rpos,smax,dt1,pi,usq,     &
                   vdotu,vrel1,area1,area2,area3,area4,p1,xnorm,tt,      &
                   smin_step,dfac,rzone2,psix,smin_type,smin_time,       &
                   dum,tstep,ptest,vrel2,eova(1)   ! eova(1) required for getsgxn
      INTEGER(i4B) i1(ktk),ir,iz,nout,izone,n1,ne1,i,j,izone1,index,izoned(1)
      REAL(SP) seed0 


      DATA       seed0    /0.0_SP/

      DATA      cvec &
           /  1.0_DP,  2.0_DP,  3.0_DP,  4.0_DP,  5.0_DP,  6.0_DP,  7.0_DP,  8.0_DP,  9.0_DP, 10.0_DP, &
             11.0_DP, 12.0_DP, 13.0_DP, 14.0_DP, 15.0_DP, 16.0_DP, 17.0_DP, 18.0_DP, 19.0_DP, 20.0_DP, &
             21.0_DP, 22.0_DP, 23.0_DP, 24.0_DP, 25.0_DP, 26.0_DP, 27.0_DP, 28.0_DP, 29.0_DP, 30.0_DP, &
             31.0_DP, 32.0_DP, 33.0_DP, 34.0_DP, 35.0_DP, 36.0_DP, 37.0_DP, 38.0_DP, 39.0_DP, 40.0_DP, &
             41.0_DP, 42.0_DP, 43.0_DP, 44.0_DP, 45.0_DP, 46.0_DP, 47.0_DP, 48.0_DP, 49.0_DP, 50.0_DP, &
             51.0_DP, 52.0_DP, 53.0_DP, 54.0_DP, 55.0_DP, 56.0_DP, 57.0_DP, 58.0_DP, 59.0_DP, 60.0_DP, &
             61.0_DP, 62.0_DP, 63.0_DP, 64.0_DP, 65.0_DP, 66.0_DP, 67.0_DP, 68.0_DP, 69.0_DP, 70.0_DP, &
             71.0_DP, 72.0_DP, 73.0_DP, 74.0_DP, 75.0_DP, 76.0_DP, 77.0_DP, 78.0_DP, 79.0_DP, 80.0_DP, &
             81.0_DP, 82.0_DP, 83.0_DP, 84.0_DP, 85.0_DP, 86.0_DP, 87.0_DP, 88.0_DP, 89.0_DP, 90.0_DP, &
             91.0_DP, 92.0_DP, 93.0_DP, 94.0_DP, 95.0_DP, 96.0_DP, 97.0_DP, 98.0_DP, 99.0_DP,100.0_DP, &
            101.0_DP,102.0_DP,103.0_DP,104.0_DP,105.0_DP,106.0_DP,107.0_DP,108.0_DP,109.0_DP,110.0_DP, &
            111.0_DP,112.0_DP,113.0_DP,114.0_DP,115.0_DP,116.0_DP,117.0_DP,118.0_DP,119.0_DP,120.0_DP, &
            121.0_DP,122.0_DP,123.0_DP,124.0_DP,125.0_DP,126.0_DP,127.0_DP,128.0_DP,129.0_DP,130.0_DP, &
            131.0_DP,132.0_DP,133.0_DP,134.0_DP,135.0_DP,136.0_DP,137.0_DP,138.0_DP,139.0_DP,140.0_DP, &
            141.0_DP,142.0_DP,143.0_DP,144.0_DP,145.0_DP,146.0_DP,147.0_DP,148.0_DP,149.0_DP,150.0_DP, &
            151.0_DP,152.0_DP,153.0_DP,154.0_DP,155.0_DP,156.0_DP,157.0_DP,158.0_DP,159.0_DP,160.0_DP, &
            161.0_DP,162.0_DP,163.0_DP,164.0_DP,165.0_DP,166.0_DP,167.0_DP,168.0_DP,169.0_DP,170.0_DP, &
            171.0_DP,172.0_DP,173.0_DP,174.0_DP,175.0_DP,176.0_DP,177.0_DP,178.0_DP,179.0_DP,180.0_DP, &
            181.0_DP,182.0_DP,183.0_DP,184.0_DP,185.0_DP,186.0_DP,187.0_DP,188.0_DP,189.0_DP,190.0_DP, &
            191.0_DP,192.0_DP,193.0_DP,194.0_DP,195.0_DP,196.0_DP,197.0_DP,198.0_DP,199.0_DP,200.0_DP/

!


!     the following assumes that data from previous
!     particle is saved. this ok only if passed through
!     argument list. HSJ
      xmassp = Proton_mass*1000._DP ! xmassp in grams
      IF (newpar .EQ. izero )  go to 100
!
! calculate times for particle to enter and exit toroidal box surrounding plasma
!
      CALL timtor (rin,rmax,x0,y0,z0,vx0,vy0,vz0,zmin,zmax,tenter,texit)

      IF (tenter .LE. -1.0e10_DP)  go to 140
!
! advance particle to edge of box
!
      x0 = x0 + vx0*tenter
      y0 = y0 + vy0*tenter
      z0 = z0 + vz0*tenter
!
 
 
      IF (nebin .NE. izero) THEN
!
! ----------------------------------------------------------------------
! follow collisionless neutral trajectory to obtain minimum mean free path.
! this is required to account for toroidal rotation.
!
! ALSO:
!    Only include those energy groups that fall in a valid flux zone,
!    (izone1 .le. mfm1) for the search in subroutine GETSGXN. The array
!    of energy values e1 is now ne1 elements long rather than n1.
!                                            Daniel Finkenthal    9-6-95
! ----------------------------------------------------------------------
!
      smax = vbeam(ie,ib) * (texit-tenter)
      n1   = 2.0 + smax/ds1
      IF (n1 .GT. ktk) THEN
        n1  = ktk
        ds1 = smax / FLOAT (n1-1)
        WRITE (nout, 200) ds1
      END IF
      dt1 = smax / (FLOAT (n1-1) * vbeam(ie,ib))
      ne1 = 0
      IF (codeid .EQ. 'onedee') THEN
        DO i=1,n1
          delt   = (cvec(i) - 1.0_DP) * dt1
          x1     = x0 + delt*vx0
          y1     = y0 + delt*vy0
          z1     = z0 + delt*vz0
          r1     = SQRT (x1**2 + y1**2)
          ir     = (r1-r(1))*dri + 1.0_DP
          iz     = (z1-z(1))*dzi + 1.0_DP
          ir     = MIN0 (ir,mim1)
          iz     = MIN0 (iz,mjm1)
          IF (ir .LE. 0 .OR. iz .LE. 0)go to 10  !HSJ 1/28/2000
          p1     = SQRT ((r1-rmajor)**2+(elongi*(z1-zax))**2)
          izone1 = p1*droti + 1.0
          IF (izone1 .LE. mfm1) THEN
            ne1     = ne1 + 1
            i1(ne1) = izone1
            usq     = ( x1**2 + y1**2 ) * zangrot(i1(ne1))**2
            vdotu   = (x1*vy0 - y1*vx0) * zangrot(i1(ne1))
            vrel1   = vbeam(ie,ib)**2 + usq - 2.0*vdotu
            e1(ne1) = 0.5 * xmassp * vrel1 * kevperg
            e1(ne1) = ABS (e1(ne1))
          END IF
        END DO
      ELSE
        DO i=1,n1
          delt  = (cvec(i) - 1.0) * dt1
          x1    = x0 + delt*vx0
          y1    = y0 + delt*vy0
          z1    = z0 + delt*vz0
          r1    = SQRT (x1**2+y1**2)
          ir    = (r1-r(1))*dri + 1.0
          iz    = (z1-z(1))*dzi + 1.0
          ir    = MIN0 (ir,mim1)
          iz    = MIN0 (iz,mjm1)
          IF (ir .LE. 0 .OR. iz .LE. 0)go to 10  !HSJ 1/28/2000
          area1 = (r1-r(ir))*(z1-z(iz))
          area2 = (r(ir+1)-r1)*(z1-z(iz))
          area3 = (r(ir+1)-r1)*(z(iz+1)-z1)
          area4 = (r1-r(ir))*(z(iz+1)-z1)
          p1    = (area3*psi(ir,iz)   + area4*psi(ir+1,iz) &
                +  area2*psi(ir,iz+1) + area1*psi(ir+1,iz+1))*dri*dzi
          p1     =  MAX (p1,psiax)
          izone1 = SQRT (p1-psiax)*drutpi + 1.0
          IF (izone1 .LE. mfm1) THEN
            ne1     = ne1 + 1
            i1(ne1) = izone1
            usq     = (x1**2+y1**2)*zangrot(i1(ne1))**2
            vdotu   = (x1*vy0-y1*vx0)*zangrot(i1(ne1))
            vrel1   = vbeam(ie,ib)**2 + usq - 2.0*vdotu
            e1(ne1) = 0.5 * xmassp * vrel1 * kevperg
            e1(ne1) = ABS (e1(ne1))
          END IF
        END DO
      END IF
 
 10   CALL getsgxn (e1, i1, ne1, ktk, ib, ie, nebin, debin, kcmp1,kbe, &
                     nout,sgxntab)


!      xnorm         = amaxaf (sgxntab, 1, ne1)
       xnorm = sgxntab(1)
       DO i =2,ne1  ! MAXVAL not restricted to i <= ne1
            IF(sgxntab(i) .GT. xnorm)xnorm  = sgxntab(i) 
       ENDDO

       sgxnmi(ie,ib) = 1.0_DP / xnorm
      END IF
!
! ----------------------------------------------------------------------
! inject neutral into plasma
! ----------------------------------------------------------------------
!
!  set coordinates and time for entering box
!
  100 CONTINUE

      xpos      = x0
      ypos      = y0
      zpos      = z0
      tt        = tenter
      izone     = mfm1 + 1                         ! initially neutral is outside the plasma
      smin_step = 0.1_DP                           ! 0.1 cm min step or
      smin_step = MIN (smin_step, smax/1000.0_DP)  ! make scale-independent
      smin_time = smin_step / SQRT (vx0**2 + vy0**2 + vz0**2)
!
!  follow particle into plasma


  110 dfac  = -LOG (RANDOM12 (seed0))
!
!     if neutral is not yet in the plasma (izone ge mf) then take steps
!     of minimum size 1 mm until we enter the plasma. we could find the
!     exact plasma boundary but that would be overkill. with 1mm step
!     size we certainly are within any physics scales we could resolve
!     near the plasma edge. this avoids wasting a lot of steps outside
!     the plasma and will be a significant savings if the cross sections
!     are large. dfac is 1/exponentially distributed, so for large steps
!     don't modify tstep. --------------------------- 27 Oct 95 ---- HSJ
!
      IF (izone .LE. mfm1) THEN                     !  inside the plasma
        tstep =      dfac * sgxnmi(ie,ib) / vbeam(ie,ib)
      ELSE
        tstep = MAX (dfac * sgxnmi(ie,ib) / vbeam(ie,ib), &
                     smin_time)                     ! outside the plasma
      END IF
!
      tt    = tt  + tstep
      IF (tt .GE. texit)  go to 140
      xpos = xpos + vx0*tstep
      ypos = ypos + vy0*tstep
      zpos = zpos + vz0*tstep
      rpos = SQRT (xpos**2 + ypos**2)
!
!  determine zone in which particle collides for 'onedee' geometry
!
      IF (codeid .EQ. 'onedee') THEN
        rzone2 = (rpos-rmajor)**2 + (elongi*(zpos-zax))**2
        rzone  = SQRT (rzone2)
        izone  = rzone*droti + 1.0
      ELSE
!
!  determine zone in which particle collides for general geometry;
!     use bilinear interpolation away from magnetic axis,
!     and biquadratic interpolation near the axis.
!
        i     = (rpos-r(1))*dri+1.0
        j     = (zpos-z(1))*dzi+1.0
        i     = MIN0 (i,mim1)
        j     = MIN0 (j,mjm1)
        psix  = MIN  (psi(i,j),psi(i+1,j),psi(i,j+1),psi(i+1,j+1))
        ptest = (psix-psiax)*(drutpi/mfm1)**2
        IF (ptest .GE. 0.02) THEN
          area1 = (rpos-r(i))*(zpos-z(j))
          area2 = (r(i+1)-rpos)*(zpos-z(j))
          area3 = (r(i+1)-rpos)*(z(j+1)-zpos)
          area4 = (rpos-r(i))*(z(j+1)-zpos)
          pzone = (area3*psi(i,j) + area4*psi(i+1,j) &
                + area1*psi(i+1,j+1) + area2*psi(i,j+1))*dri*dzi
        ELSE
          CALL pfit(psi(i-1,j-1),r(i-1),z(j-1),rpos,zpos,ki,pzone,dum, &
                    dum)
        END IF
        pzone =  MAX (pzone,psiax)
        izone = SQRT (pzone-psiax)*drutpi + 1.0
      END IF
!
!     if particle has psuedo-collision, continue following particle.
!     if particle has real collision, return.
!
      IF (izone .GT. mfm1)  go to 110 ! the particle is inside the box..
!                                     ..but still outside the plasma
      IF (nebin .NE. izero   ) THEN
        usq   = (rpos*zangrot(izone))**2
        vdotu = (xpos*vy0-ypos*vx0)*zangrot(izone)
        vrel2 = vbeam(ie,ib)**2 + usq - 2.0*vdotu
        eova(1)  = 0.5_DP * xmassp*vrel2*kevperg
        eova(1)  = ABS (eova(1))
        izoned(1) = izone ! for ifort and declaration of getsgxn
        CALL getsgxn (eova, izoned, 1, ktk, ib, ie, nebin, &
                      debin, kcmp1,kbe,nout,sgxnloc)

      ELSE
        index      = ke*(ib-1) + ie
        sgxnloc(1) = sgxn(1,izone,index,1)  ! fraction of reactions producing electrons;
        sgxnloc(2) = sgxn(2,izone,index,1)  ! fraction of reactions producing species 1 ion
        sgxnloc(3) = sgxn(3,izone,index,1)  ! fraction of reactions producing species 2 ion
        sgxnloc(4) = sgxn(4,izone,index,1)  ! total inverse mean free path
      END IF
!     if (RANF   (     ) .gt. sgxnloc(4) * sgxnmi(ie,ib))  go to 110
      IF (RANDOM12 (seed0) .GT. sgxnloc(4) * sgxnmi(ie,ib))  go to 110
      RETURN
!
!  set flag to indicate that particle hit wall
!
  140 izone = mfm1 + 1
  200 FORMAT (' WARNING from subroutine INJECT:'                  / &
              '         maximum number of grid elements exceeded' / &
              '         increasing ds1 to ', e10.3, ' cm')
      RETURN
!
      END SUBROUTINE inject


      SUBROUTINE logint (x, y)
! ----------------------------------------------------------------------
! interpolates y(x) quadratically and logarithmically
! ----------------------------------------------------------------------
!

      IMPLICIT NONE

      REAL(DP)   xdat(15), ydat(15)

      REAL(DP)   ylogp,ylog0,ylogm,x,y,d0,d0m,dp0,dpm,dmp,facm,fac0,facp, &
                 ylog,dm,dpq,d0p,dm0
      INTEGER(I4B) i0,mdat,mdatm,im,ip
!
      DATA xdat /4.0_DP,6.0_DP,8.0_DP,10.0_DP,20.0_DP,30.0_DP,40.0_DP,  &
                60.0_DP,80.0_DP,100.0_DP, &
                200.0_DP,300.0_DP,400.0_DP,600.0_DP,800.0_DP/

      DATA ydat /8.95e-01_DP,8.75e-01_DP,8.70e-01_DP,8.65e-01_DP,8.20e-01_DP, &
                7.25e-01_DP,6.25e-01_DP,4.40e-01_DP,2.90e-01_DP,1.90e-01_DP,2.40e-02_DP, &
                5.25e-03_DP,1.20e-03_DP,1.60e-04_DP,5.40e-05_DP/
!
      mdat  = 15
      mdatm = mdat - 1
!
      DO i0=2,mdatm
        IF (xdat(i0) .GE. x)  go to 11
      END DO
!
   11 IF (i0 .GT. mdatm) i0 = mdatm
      im = i0-1
      ip = i0+1
      ylogp = LOG (ydat(ip))
      ylog0 = LOG (ydat(i0))
      ylogm = LOG (ydat(im))
      dm = x-xdat(im)
      d0 = x-xdat(i0)
      dpq = x-xdat(ip)
      d0m = xdat(i0)-xdat(im)
      dp0 = xdat(ip)-xdat(i0)
      dpm = xdat(ip)-xdat(im)
      dm0 = -d0m
      d0p = -dp0
      dmp = -dpm
      facm = d0*dpq/(dm0*dmp)
      fac0 = dm*dpq/(d0m*d0p)
      facp = dm*d0/(dpm*dp0)
      ylog = facm*ylogm+fac0*ylog0+facp*ylogp
      y = EXP (ylog)
!
      RETURN
!
      END SUBROUTINE logint




      SUBROUTINE nbdep2d (psi, nw,nwx2,nh,nhx2,nwh,mi, mj, r, z, potsid, & 
                          rzpat, nrpat,nzpat, me, mb, vbeam, hxfrac)
! ----------------------------------------------------------------------
! calculates neutral beam deposition on (r,z) grid.  grid size
! determined by nrpat and nzpat, the number of equally spaced
! elements in the r and z axes.  outputs the 3rd excited state
! component, rzhex, to file 'beamdep' for postprocessing.
! ----------------------------------------------------------------------
!
!
       USE nf_param,                           ONLY : ke,kb

       USE xsct,                               ONLY : sgxn

       USE zonal_data,                         ONLY : mf,mfm1

       USE replace_imsl,                       ONLY : my_ibcieu

      IMPLICIT NONE



!       INTEGER,PARAMETER ::kzm1 = kz-1
!  argument list:
      REAL(DP)     hxfrac(ke,kb), vbeam(ke,kb), &
                   potsid(*),r(nw),z(nh),       &
                   psi(nw,nh),rzpat(nwx2,nhx2,ke,kb)
      INTEGER(I4B) mi,mj,me,mb,nzpat,nrpat,nw,nh,nwx2,nhx2,nwh


! local storage:                  
      REAL(DP)     bpar1(4),bpar2(4),dr,dz,                 &
                   psikp(nw*nh),psipat(nwx2,nhx2),rr(nwx2), &
                   rznub(nwx2,nhx2,ke,kb),rzsig(nwh), &
                   zz(nh2)

!      REAL(DP)     capsig(kz),cspln1(kz-1,3),pc(kz),

      REAL(DP), ALLOCATABLE,DIMENSION(:)           ::  splnwk,capsig,pc

      REAL(DP), ALLOCATABLE,DIMENSION(:,:)         ::  cspln1

      INTEGER(I4B) i,ib,ie,j,inc(nh2),isupp(4,ke,kb),igrid,ndum,ier,  &
                   ifail,ind,n,kzm1

      INTEGER(I2B)  get_next_io_unit,nout
      CHARACTER*4   string
      CHARACTER *36 filename

      LOGICAL    exist
!
!     set (r,z) grid modification option
!
      IF(.NOT. ALLOCATED(capsig))ALLOCATE(capsig(mf))
       capsig(:) = 0.0_DP
      IF(.NOT. ALLOCATED(pc))ALLOCATE(pc(mf))
      pc(:)= 0.0_DP
      IF(.NOT. ALLOCATED(cspln1))ALLOCATE(cspln1(mfm1,3))
      cspln1(:,:)= 0.0_DP
      IF(.NOT. ALLOCATED(splnwk))ALLOCATE(splnwk(3*(nh-1)+nh))
      splnwk(:)= 0.0_DP
      kzm1 = mfm1
      igrid = 0
      IF (nrpat .NE. mi .OR. nzpat .NE. mj)  igrid = 1 ! per Jeff Moller
!
!     zero out arrays
!
      DO i=1,4
        bpar1(i) = 0.0_DP
        bpar2(i) = 0.0_DP
      END DO
!
!     set up integer array to allow vectorization further on
!
      ndum = MAX (nrpat, nzpat)
      IF (ndum .GT. 2*nh)THEN
!        CALL STOP ('subroutine NBDEP2D: Jeff Moller bailout', 68)
         WRITE(ncrt,FMT='("error by myid,ndum,2*nh =",3(i5,x))')myid,ndum,2*nh
         lerrno = iomaxerr + 227_I4B
         CALL terminate(lerrno,nlog)
      ENDIF
      DO i=1,ndum
        inc(i) = i-1
      END DO
!
!     get psi on zone centers
!
!      mfm1 = mf-1
      DO i=1,mfm1
        pc(i) = 0.5 * (potsid(i) + potsid(i+1))
      END DO
!
!     if user defined new (r,z) grid, interpolate psi(i,j) onto it
!
      IF (igrid .EQ. 1) THEN
        dr      = (r(mi)-r(1)) / (nrpat-1)
        rr(1)   = r(1)
        DO i=2,nrpat
          rr(i) = rr(1)+inc(i)*dr
        END DO
        dz      = (z(mj)-z(1))/(nzpat-1)
        zz(1)   = z(1)
        DO i=2,nzpat
          zz(i) = zz(1)+inc(i)*dz
        END DO
!        CALL ibcieu (psi, nw, r, mi, z, mj, rr, nrpat, zz, nzpat, &
!                     psipat, 2*nw, splnwk, ier)
        CALL my_ibcieu (psi, nw, r, mi, z, mj, rr, nrpat, zz, nzpat, &
                     psipat, 2*nw, splnwk, ier)
      END IF
!
!     begin loop over beam and beam energy
!
      DO   ib=1,mb
        DO ie=1,me
!
!         get an estimate of the support of rzpat
!
          CALL support (rzpat, nrpat, nzpat, ie, ib,nw,nwx2,nh,nhx2, isupp, ifail)
!
!         get spline fits to macroscopic neutral beam attenuation cross
!         sections as a function of psi
!
          ind = ke*(ib-1) + ie
          DO i=1,mfm1
            capsig(i) = sgxn(4,i,ind,1)
          END DO
!
          CALL icsicu1 (pc, capsig, mfm1, bpar1, cspln1, kzm1, ier)
!
!         loop over (r,z) points
!
          n = izero ! added 2/1/11 HSJ
          DO   i=isupp(1,ie,ib),isupp(2,ie,ib)
            DO j=isupp(3,ie,ib),isupp(4,ie,ib)
              IF (rzpat(i,j,ie,ib) .NE. 0) THEN
                n = n+1
                IF (igrid .EQ. 0) THEN
                  psikp(n) = psi   (i,j)
                ELSE
                  psikp(n) = psipat(i,j)
                END IF
              END IF
            END DO
          END DO
!
          CALL icsevu1 (pc, capsig, mfm1, cspln1, kzm1, &
                        psikp, rzsig, n,ier)
          n = 0
          DO   i=isupp(1,ie,ib),isupp(2,ie,ib)
            DO j=isupp(3,ie,ib),isupp(4,ie,ib)
              IF (rzpat(i,j,ie,ib) .NE. 0.0_DP) THEN
                n = n + 1
                rznub(i,j,ie,ib) = rzpat(i,j,ie,ib) &
                                 / (vbeam(ie,ib) * rzsig(n))
              END IF
            END DO
          END DO
        END DO
      END DO


!
!     output subset of rzhex to file 'beamdep'
!
!      INQUIRE (unit = nout, exist = exist)
       WRITE(string,FMT='(i3)' )myid
       filename = 'beamdep'//ADJUSTL(string)
       INQUIRE (FILE =filename, EXIST  = exist)
        nbdep = get_next_io_unit()
        nout = nbdep
      OPEN  (unit = nout, file = filename, status = 'REPLACE')
      WRITE (nout, 1000)  nrpat, nzpat, me, mb
      IF (igrid .EQ. 0) THEN
        WRITE (nout, 1010)   r(1),  r(mi)   ,  z(1),  z(mj)
      ELSE
        WRITE (nout, 1010)  rr(1), rr(nrpat), zz(1), zz(nzpat)
      END IF



!
      DO   ib=1,mb
        DO ie=1,me
          WRITE (nout, 1000) (isupp(i,ie,ib),i=1,4)
          WRITE (nout, 1010) hxfrac(ie,ib)
          WRITE (nout, 1010) &
          ((rznub(i,j,ie,ib), i=isupp(1,ie,ib),isupp(2,ie,ib)), &
                              j=isupp(3,ie,ib),isupp(4,ie,ib))
        END DO
      END DO
!
 1000 FORMAT (4(5x,i4))
 1010 FORMAT (7(2x,e16.7))
!      CALL giveupus(nout)
      CLOSE (unit = nout)
      RETURN
!
      END SUBROUTINE nbdep2d


      SUBROUTINE NF_nubplt
!-----------------------------------------------------------------------
!-- call after postnub (hibrz ==> hibr,etc)
!-----------------------------------------------------------------------
      USE nf_param,                                     ONLY : ke

      USE Plasma_properties ,                           ONLY : profile,neut_beam

      USE tport_mhd_grid_data,                          ONLY : rplasmax,rplasmin,  &
                                                               hibr,hdep,angmpf

      USE grid_class,                                   ONLY : nj,rho_grid

      USE zonal_data,                                   ONLY : rotsid,mfm1,max_ctr_pts

      USE P_nfreya_interface,                           ONLY : beam_data,d_fast_ion

      USE nub,                                          ONLY : hdepsmth,fidiff_on, &
                                                               bfr_neutrlz



      IMPLICIT NONE
      INTEGER(I4B) i,ie,j,jb,hdep_int,last_point
      REAL(DP),DIMENSION(ke) ::  adiff_ap,adiff_0p,adiff_xpinp,adiff_xpoutp
      hdep_int = INT(hdepsmth)

!     IF (beam_on(1) .GT. timmax)  RETURN

      ! The _izpt and pitch_a  are loaded on master for the total
      ! number of computational injectors used ( not jsut the number of physical
      ! injectors.) Below we only write out the data for first set of injector beamlines
      ! and neglect the psuedo bemaline output .
      ! the plotting program only uses a single n_izpt. So until we fix the plot program
      ! we will simply take a common value:
      ! NOTE we have neut_beam%nsample_izpt(1:3,1:no_physical_injectors) loaded but arent using it.
      ! to simplify plotting use the following, where all arrays are dumped at size n_izpt
      ! even if output doesnt exist (just copy last point) HSJ
      !n_izpt = MINVAL(nsample_izpt)
      n_izpt  = MAXVAL(nsample_izpt)
      DO jb   = 1,no_physical_injectors
         DO ie = 1,ke
            last_point = nsample_izpt(ie,jb)
            IF(last_point == 0)THEN
               last_point = MAX(last_point,1)
               x_izpt(last_point,ie,jb)    = MAXVAL(r_zone_cntr)
               y_izpt(last_point,ie,jb)    = MAXVAL(r_zone_cntr)
               z_izpt(last_point,ie,jb)    = MAXVAL(z_zone_cntr)
               r_izpt(last_point,ie,jb)    = MAXVAL(r_zone_cntr)
               vx_izpt(last_point,ie,jb)   = zeroc
               vy_izpt(last_point,ie,jb)   = zeroc
               vy_izpt(last_point,ie,jb)   = zeroc
               pitch_a(last_point,ie,jb)   = zeroc
            ENDIF
            IF( last_point  .LT. n_izpt)THEN
               DO i = last_point, n_izpt
                  x_izpt(i,ie,jb)  = x_izpt(last_point,ie,jb)
                  y_izpt(i,ie,jb)  = y_izpt(last_point,ie,jb)
                  z_izpt(i,ie,jb)  = z_izpt(last_point,ie,jb)
                  r_izpt(i,ie,jb)  = z_izpt(last_point,ie,jb)
                  vx_izpt(i,ie,jb) = vx_izpt(last_point,ie,jb)
                  vy_izpt(i,ie,jb) = vy_izpt(last_point,ie,jb)
                  vz_izpt(i,ie,jb) = vy_izpt(last_point,ie,jb)
                  pitch_a(i,ie,jb) = pitch_a(last_point,ie,jb)
               ENDDO
            ENDIF
          ENDDO
      ENDDO



      WRITE (nfreya_plot_unit, '(a)') '****continue****'
      WRITE (nfreya_plot_unit,  8005)  time_now,profile%ene%data(1)*im32icm3, &
                                       profile%te%data(1),rplasmax,rplasmin,n_izpt
      WRITE (nfreya_plot_unit,   8030)((nsample_izpt(ie,jb),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,   8021) (((x_izpt(i,ie,jb),   i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,  8021) (((y_izpt(i,ie,jb),   i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,  8021) (((z_izpt(i,ie,jb),   i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,  8021) (((r_izpt(i,ie,jb),   i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,  8021) (((vx_izpt(i,ie,jb),  i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,  8021) (((vy_izpt(i,ie,jb),  i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,  8021) (((vz_izpt(i,ie,jb),  i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)
      WRITE (nfreya_plot_unit,  8021) (((pitch_a(i,ie,jb),  i=1,n_izpt),ie = 1,ke),jb = 1,no_physical_injectors)

!----------------------------------------------------------------------------------------
! there is an iborb clause (eg iborb > 1) at this location  in  P_NF_nubplt 
! that is not printed out here because it is assumed that iborb <= 1 here.
!----------------------------------------------------------------------------------------



       ! ----------------------------------------------------------------------------------------
       ! fast ion diffusion arrays wont be set if option was not used so just put
       ! zeros in plotfile to avoid having to carry logic into plot program
       ! (plot program checks for valid data using d_fast_ion%fidif_on)
       ! ----------------------------------------------------------------------------------------
          IF(.NOT. ASSOCIATED(d_fast_ion%fi_den))THEN
             ALLOCATE(d_fast_ion%fi_den(nj,ke,no_physical_injectors), &
                      d_fast_ion%fi_d(nj,ke,no_physical_injectors))
            d_fast_ion%fi_d(:,:,:)   = zeroc
            d_fast_ion%fi_den(:,:,:) = zeroc
          ENDIF

      DO jb=1,no_physical_injectors
 
        WRITE (nfreya_plot_unit, 8021) (neut_beam%bneut(i,jb),i=1,3)
        WRITE (nfreya_plot_unit, 8021) (neut_beam%ebeam(i,jb),i=1,3)
        WRITE (nfreya_plot_unit, 8021) (neut_beam%pbeam(i,jb),i=1,3)
        WRITE (nfreya_plot_unit, 8021) (neut_beam%fap(i,jb),i=1,3)
        WRITE (nfreya_plot_unit, 8021) (neut_beam%fwall(i,jb),i=1,3)
        WRITE (nfreya_plot_unit, 8021) (neut_beam%forb(i,jb),i=1,3)
        DO ie=1,3
          WRITE (nfreya_plot_unit, 8021) (neut_beam%hibr  (j,ie,jb), j=1,nj)
          WRITE (nfreya_plot_unit, 8021) (neut_beam%hdep  (j,ie,jb), j=1,nj)
          WRITE (nfreya_plot_unit, 8021) (neut_beam%angmpf(j,ie,jb), j=1,nj)
          WRITE (nfreya_plot_unit, 8021) (neut_beam%sb(j,ie,jb), j=1,nj)
          WRITE (nfreya_plot_unit, 8021) (d_fast_ion%fi_den(j,ie,jb),j=1,nj)
          WRITE (nfreya_plot_unit, 8021) (d_fast_ion%fi_d(j,ie,jb),j=1,nj)
        END DO
        ! some (not all) new  output for NF_nubplt
        WRITE(nfreya_plot_unit, 8030)(neut_beam%nmbrz(i,jb),i=1,mfm1)
        WRITE(nfreya_plot_unit, 8021)(fb00(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(fb01(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(fb10(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(fb11(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(fber(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(wb00(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(wb01(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(wb10(ie,jb),ie=1,3)
        WRITE(nfreya_plot_unit, 8021)(wb11(ie,jb),ie=1,3)
      END DO
        ! new output for NF_nubplt
      WRITE(nfreya_plot_unit, 8021)(rotsid(i),i = 1,mf) ! rotsid is in cm ??




      WRITE (nfreya_plot_unit, 8021)  (rho_grid%data (j)*m2cm, j=1,nj)
      WRITE (nfreya_plot_unit, 8101)   max_ctr_pts,d_fast_ion%fidif_on,hdep_int
      WRITE (nfreya_plot_unit, 8107)   bfr_neutrlz

      adiff_0p(:)       = d_fast_ion%adiff_0l
      adiff_ap(:)       = d_fast_ion%adiff_al
      adiff_xpinp(:)    = d_fast_ion%adiff_xpinl
      adiff_xpoutp(:)   = d_fast_ion%adiff_xpoutl
      WRITE (nfreya_plot_unit, 8021)  adiff_ap,adiff_0p,adiff_xpinp,adiff_xpoutp

      DO i=1,mf-1

        WRITE (nfreya_plot_unit, 8100)   n_zone_cntr(i)

        WRITE (nfreya_plot_unit, 8120)  (r_zone_cntr(j,i),z_zone_cntr(j,i), j=1,n_zone_cntr(i))
      END DO

      DO i=1,no_physical_injectors
         WRITE(nfreya_plot_unit, 8150)beam_data%nlco(i)
         WRITE(nfreya_plot_unit, 8155)beam_data%beam_id(i)(1:LEN_TRIM(beam_data%beam_id(i)))
      ENDDO

      RETURN 

 8100 FORMAT (i6)
 8101 FORMAT (i6,2x,i2,2x,i4)
 8107 FORMAT (2x,L)
 8120 FORMAT (6e12.5)
 8150 FORMAT (L)
 8155 FORMAT (a)
 8005 FORMAT (5e12.3,2i10)
 8021 FORMAT (6e14.6)
 8030 FORMAT (6(i7,x))

      END SUBROUTINE NF_nubplt





      SUBROUTINE process_nubplt
!------------------------------------------------------------------------
! -- handle file structure for nubplt plotting
!------------------------------------------------------------------------
    USE nrtype,                     ONLY : I4B,DP

    USE common_constants,           ONLY : izero

    USE error_handler,              ONLY : load_errno,lerrno,iomaxerr,terminate

    USE io_gcnmp,                   ONLY : ncrt,nlog


      IMPLICIT NONE
          INTEGER(I4B)task,error

           task = 1      ! write mhd related data
           error = izero 
           CALL setup_plot_file(task,error)          ! nfreya_init.f90
           IF(error .ne. izero)THEN
              lerrno = iomaxerr + 242_I4B
              CALL terminate(lerrno,nlog)
           ENDIF

           task = 3      ! write nfreya calculated data (calls  nubplt)
           error = izero ! temp
           CALL setup_plot_file(task,error)
           IF(error .ne. izero)THEN
              lerrno = iomaxerr + 242_I4B
              CALL terminate(lerrno,nlog)
           ENDIF

           task = 2      ! write closing clause
           error = izero 
           CALL setup_plot_file(task,error)
           IF(error .ne. izero)THEN
              lerrno = iomaxerr + 242_I4B
              CALL terminate(lerrno,nlog)
           ENDIF

           task = -1    ! close the file
           CALL setup_plot_file(task,error)
           IF(error .ne. izero)THEN
              lerrno = iomaxerr + 243_I4B
              CALL terminate(lerrno,nlog)
           ENDIF

          RETURN

      END SUBROUTINE process_nubplt


      SUBROUTINE orbit_12 (atwb,b1ins,b1ots,b2ins,b2ots, codeid, &
                        ic,idebug,iskip,izi,mf,norb, &
                        pinsid,potsid,pzone,rinsid,rotsid,rzone,ri,v, &
                        zetai,zi,finsid,fotsid,zinsid,zotsid, &
                        ipass,iaxis,ier,izp,wid,wt,zeta)
!------------------------------------------------------------------------
!  This subroutine calculates the following orbit parameters for a
!     fast ion:
!       ipass: 1, particle is passing
!              0, particle is trapped
!       iaxis: 1, particle circles the axis
!              0, particle does not circle the axis
!       ier  : 1, routine failed to obtain orbit
!              0, routine obtained reasonable orbit
!       izp  : outermost zone of orbit
!       wid  : width of orbit
!       wt   : weight giving approximate fraction of time spent
!              by fast ion in each zone
!       zeta : approximate pitch angle cosine in each zone along
!              orbit
!--------------------------------------------------------------------------


      IMPLICIT NONE
      CHARACTER(*)  codeid 
      CHARACTER  flagm*1, flagp*1
      REAL(DP)      b1ins(*), b1ots(*), b2ins(*), b2ots(*)
      INTEGER(I4B)  idebug(*),mfm1,mf,norb,ier,izzm,izp,izi,           &
                    i,ii,ipass,iaxis,nin,not,im,ic,iskip,izm,ip,       &
                    isp,ism
      REAL(DP)      pinsid(*), potsid(*), rinsid(*), rotsid(*)
      REAL(DP)      finsid(*), fotsid(*), zinsid(*), zotsid(*)
      REAL(DP)      wt(*), zeta(*),v,rzone,zetai,zi,wtp
      REAL(DP)      rin(3), rot(3),rmajor,psiax,rm,rp,zetam,zetap,wid, &
                    xmass,atwb,psii,pzone,psiref,pang,c1,c2,c3,rtip,    &
                    eoverc,xmassp,zero,drini,droti,drutpi,b1i,b2i,     &
                    b1tip,psitip,pangt,zetax2,ri,pangx,psip,b2p,xzp,   &
                    psim,b2m,xzm,wtm,wtx,zdif,zetax
      DATA eoverc / 1.600e-20_DP /

!
!     initialize parameters
!

      xmassp = Proton_mass*1000._DP ! xmassp in grams
      zero   = zeroc
      mfm1   = mf-1
      drini  = mfm1/(rinsid(1)-rinsid(mf))
      droti  = mfm1/(rotsid(mf)-rotsid(1))
      drutpi = mfm1 / SQRT (potsid(mf)-potsid(1))

      flagm  = ' '
      flagp  = ' '
      ier    = 0
      izm    = izi
      izp    = izi
      DO 10 i=1,3
      rin(i) = 0.0_DP
   10 rot(i) = 0.0_DP
      rm     = 0.0_DP
      rp     = 0.0_DP
      zetam  = 0.0_DP
      zetap  = 0.0_DP
      wid    = 0.0_DP
      DO 20 i=1,mfm1
      wt(i)  = 0.0_DP
   20 zeta(i) = 0.0_DP
      wt(izi) = 1.0
      zeta(izi) = zetai
!
! calculate ion mass times c/e; initialize rmajor and psiax
!
      xmass  = atwb*xmassp/eoverc
      rmajor = rotsid(1)
      psiax  = potsid(1)
!
!  calculate poloidal magnetic flux at the fast ion birth point
!
      IF (codeid .EQ. 'onedee') THEN
        psii = yinter(droti,mfm1,zeroc,rzone,potsid)
      ELSE
        psii = pzone
      END IF
 
!
!  approximate ratios of magnetic fields at the fast ion birth point
!
      IF (ri .GE. rmajor) THEN
        b1i = yinter(droti,mfm1,rmajor,ri,b1ots)
        b2i = yinter(droti,mfm1,rmajor,ri,b2ots)
      ELSE
        b1i = yinter(drini,mfm1,rmajor,ri,b1ins)
        b2i = yinter(drini,mfm1,rmajor,ri,b2ins)
      END IF
!
!  calculate constants for orbit
!
      psiref = xmass*v*ri
      pang   = psiref*zetai/b2i - psii
      c1     = (1.0-zetai**2)*ri/b1i
      c2     = ri/psiref
      c3     = pang
!
!  determine whether fast ion is passing or trapped
!
      IF (ABS (zetai) .EQ. 1.0) THEN
        ipass = 1
      ELSE IF (zetai .EQ. 0.0_DP) THEN
        ipass = 0
      ELSE
        rtip = (1.0-zetai**2)*ri/b1i
        IF (rtip .GE. rmajor) THEN
          b1tip  = yinter(droti,mfm1,rmajor,rtip,b1ots)
          rtip   = rtip*b1tip
          psitip = yinter(droti,mfm1,rmajor,rtip,potsid)
        ELSE
 
          b1tip  = yinter(drini,mfm1,rmajor,rtip,b1ins)
          rtip   = rtip*b1tip
          psitip = yinter(drini,mfm1,rmajor,rtip,pinsid)
        END IF
        pangt = -psitip
        IF (pang .GT. pangt) THEN
          ipass = 1
        ELSE
          ipass = 0
        END IF
      END IF
!
!  determine whether ion circles the magnetic axis or not
!
      iaxis  = 1
      zetax2 = 1.0_DP - (1.0_DP-zetai**2)*ri/(rmajor*b1i)
      IF (zetax2 .LE. 0.0_DP)  go to 80
      zetax  = SQRT (zetax2)
      IF (ipass .EQ. 0 .OR. zetai .LT. 0.0_DP) zetax = -zetax
      pangx = xmass*v*zetax*rmajor - psiax
      IF (ipass .EQ. 1 .AND. pang .LT. pangx)  go to 90
      IF (ipass .EQ. 0 .AND. pang .GT. pangx)  go to 90
   80 iaxis = 0
!
!  find roots on inside and outside of axis
!
   90 DO i=1,mf
        finsid(i) = rinsid(i)**2 - c1*rinsid(i)*b1ins(i) &
                    - (c2*b2ins(i)*(c3+pinsid(i)))**2
        fotsid(i) = rotsid(i)**2 - c1*rotsid(i)*b1ots(i) &
                    - (c2*b2ots(i)*(c3+potsid(i)))**2
      END DO
!
      nin = 0
      not = 0
!
      DO i=1,mfm1
        IF (finsid(i) .EQ. 0.0_DP .OR. finsid(i)*finsid(i+1) .LT. 0.0_DP) THEN
          nin = nin+1
          IF(nin .LE. 3) &
           rin(nin) = (rinsid(i+1)*finsid(i)-rinsid(i)*finsid(i+1)) &
                    /(finsid(i)-finsid(i+1))
        END IF
        IF (fotsid(i) .EQ. 0.0_DP .OR. fotsid(i)*fotsid(i+1) .LT. 0.0_DP) THEN
          not = not+1
          IF(not .LE. 3) &
           rot(not) = (rotsid(i+1)*fotsid(i)-rotsid(i)*fotsid(i+1)) &
                    /(fotsid(i)-fotsid(i+1))
        END IF
      END DO
!
!  match roots to orbit extrema
!
      IF (nin .GT. 3 .OR. not .GT. 3) THEN
        flagp = '*'
        go to 340
      END IF
      IF (ipass .EQ. 0 .AND. iaxis .EQ. 0) THEN
        IF (not .LT. 2)  go to 210
        rp = rot(2)
        rm = rot(1)
      ELSE IF (ipass .EQ. 1 .AND. iaxis .EQ. 1) THEN
        IF (zetai .GE. 0.0_DP) THEN
          IF (nin .EQ. 1) THEN
            IF (not .EQ. 0)  go to 210
            rp = rot(1)
            rm = rin(1)
          ELSE IF (nin .EQ. 2) THEN
            IF (not .LE. 1)  go to 210
            rp = rot(2)
            rm = rin(2)
            im = (rmajor-rm)*drini + 2.0
            IF (finsid(im-2) .LT. 0.0_DP) THEN
              rm = rfine(b1ins,b2ins,c1,c2,c3,finsid,im,1,mfm1,pinsid, &
                         rinsid,rmajor)
            END IF
          ELSE IF (nin .EQ. 0) THEN
            IF (not .LE. 1)  go to 210
            rp = rot(2)
            im = imax(finsid,mfm1)
            IF (finsid(im+1) .GT. finsid(im-1)) im = im+1
            rm = rfine(b1ins,b2ins,c1,c2,c3,finsid,im,1,mfm1,pinsid, &
                       rinsid,rmajor)
            IF (rm .EQ. 0.0_DP) THEN
              flagm = '*'
              go to 340
            END IF
          ELSE
            IF (not .EQ. 0)  go to 210
            rp = rot(1)
            rm = rin(3)
          END IF
        ELSE
          IF (nin .EQ. 0)  go to 210
          rp = rin(1)
          rm = rot(1)
          ip = (rmajor-rp)*drini + 1.0
          IF (finsid(ip+2) .LT. 0.0_DP) THEN
            rp = rfine(b1ins,b2ins,c1,c2,c3,finsid,ip,-1,mfm1,pinsid, &
                       rinsid,rmajor)
          END IF
        END IF
      ELSE IF (ipass .EQ. 1 .AND. iaxis .EQ. 0) THEN
        IF (zetai .GE. 0.0_DP) THEN
          IF (not .LT. 2)  go to 210
          rp = rot(2)
          rm = rot(1)
        ELSE
          IF (nin .LT. 2)  go to 210
          rp = rin(2)
          rm = rin(1)
        END IF
      ELSE IF (ipass .EQ. 0 .AND. iaxis .EQ. 1) THEN
        IF (not .EQ. 0)  go to 210
        rp = rot(1)
        rm = rin(1)
      END IF
      go to 220
  210 izp = mf
      go to 350
!
!  calculate poloidal magnetic flux and index of outermost zone of orbit
!
  220 IF (rp .GE. rmajor) THEN
        psip = yinter(droti,mfm1,rmajor,rp,potsid)
        b2p = yinter(droti,mfm1,rmajor,rp,b2ots)
      ELSE
        psip = yinter(drini,mfm1,rmajor,rp,pinsid)
        b2p = yinter(drini,mfm1,rmajor,rp,b2ins)
      END IF
      IF (codeid .EQ. 'onedee') THEN
        xzp = ABS (rp-rmajor)*droti + 1.0
      ELSE
        xzp = SQRT (psip-psiax)*drutpi + 1.0
      END IF
      izp = xzp
      izp = MIN0 (izp,mf)
      zetap = (pang+psip)*b2p/(xmass*v*rp)
      isp = 1
      IF (zetai .LT. 0. .AND. ipass .EQ. 1) isp = -1
      IF (isp*zetap .LT. -0.01 .OR. izp .LT. izi-1) THEN
        flagp = '*'
        go to 340
      END IF
!
!  determine whether ion is confined or lost; return if it is lost
!
      IF (izp .GT. mfm1)  go to 350
!        if(rm .eq. 0)then
!           print *,ipass,iaxis
!           print *,'rin =',rin
!           print *,'rot =',rot
!           print *,'zetai,nin =',zetai,nin
!        endif
      IF(rm .LE. zeroc) go to 350  !ctr injection lost orbit ??
!
!  calculate poloidal magnetic flux and index of innermost zone of orbit
!
      IF (rm .LT. rmajor) THEN
        psim = yinter(drini,mfm1,rmajor,rm,pinsid)
        b2m = yinter(drini,mfm1,rmajor,rm,b2ins)
      ELSE
        psim = yinter(droti,mfm1,rmajor,rm,potsid)
        b2m = yinter(droti,mfm1,rmajor,rm,b2ots)
      END IF
      IF (codeid .EQ. 'onedee') THEN
        xzm = ABS (rmajor-rm)*droti + 1.0
      ELSE
        xzm = SQRT (psim-psiax)*drutpi + 1.0
      END IF
      izm = xzm
      zetam = (pang+psim)*b2m/(xmass*v*rm)
      ism = -1
      IF (zetai .GE. 0. .AND. ipass .EQ. 1) ism = 1
      IF (ism*zetam .LT. -0.01 .OR. izm .GT. izi+1) THEN
        flagm = '*'
        go to 340
      END IF
!
!  calculate orbit width
!
      wid = ABS (rp-rmajor) - ABS (rm-rmajor)
      wid = MAX (wid, zero)
!
!  calculate weights by zone; assume ion spends equal time in each
!     zone it fully traverses and proportionately less time in the
!     innermost and outermost zones of the orbit
!
      IF (izp .EQ. izm)  go to 350
      wtm = 1.0 - (xzm - izm)
      wtp = xzp - izp
      wtx = 1.0/(wtm+wtp+izp-izm-1.0)
      wt(izm) = wtm*wtx
      wt(izp) = wtp*wtx
      IF (izp .NE. izm+1) THEN
      DO 310 i=izm+1,izp-1
  310 wt(i) = wtx
      END IF
!
!  calculate pitch angle cosines by zone; assume cosine varies
!     linearly with zone index
!
      zdif      = zetap - zetam
      zeta(izm) = zetam + 0.5 * wt(izm)*zdif
      zeta(izp) = zetap - 0.5 * wt(izp)*zdif
!
      IF (izp .NE. izm+1) THEN
        DO i=izm+1,izp-1
          zeta(i) = zeta(i-1) + wt(i)*zdif
        END DO
      END IF
!
      go to 350
!
!  set error flag if error was detected
!
  340 izp = izi
      ier = 1
!
!  print out selected results
!
  350 IF (norb .EQ. 0)  RETURN
      IF (  ic .EQ. 1)  WRITE (norb, 1010)
      IF (MOD (ic-1,iskip) .NE. 0)  RETURN
      WRITE (norb, 1020)  ic, ipass, iaxis, ism, isp, nin, not, ri, zi, &
                          zetai, rm, flagm, zetam, rp, flagp, zetap, &
                          wid, izi, izm, izp
      IF (ic .NE. idebug(1) .AND. ic .NE. idebug(2) .AND. &
          ic .NE. idebug(3) .AND. ic .NE. idebug(4) .AND. &
                                  ic .NE. idebug(5))  RETURN
      DO i=1,mf
        zinsid(i) = (pang+pinsid(i))*b2ins(i)/(xmass*v*rinsid(i))
        zotsid(i) = (pang+potsid(i))*b2ots(i)/(xmass*v*rotsid(i))
      END DO
      DO ii=mf,11,-10
        WRITE (norb,2010) (i, i=ii,ii-9,-1)
        WRITE (norb,2020) (rinsid(i), i=ii,ii-9,-1)
        WRITE (norb,2030) (zinsid(i), i=ii,ii-9,-1)
        WRITE (norb,2040) (finsid(i), i=ii,ii-9,-1)
      END DO
      DO ii=1,mfm1,10
        WRITE (norb,2010) (i, i=ii,ii+9)
        WRITE (norb,2020) (rotsid(i), i=ii,ii+9)
        WRITE (norb,2030) (zotsid(i), i=ii,ii+9)
        WRITE (norb,2040) (fotsid(i), i=ii,ii+9)
      END DO
      RETURN
!
 1010 FORMAT (/' orbit parameters' // &
        4x,'ic  ip  ix ism isp nin not',7x,'ri',7x,'zi',4x,'zetai', &
        7x,'rm',4x,'zetam',7x,'rp',4x,'zetap',5x,'wid','  izi izm izp')
 1020 FORMAT (i6,6i4,2f9.2,f9.4,2(f9.2,a1,f8.4),f8.2,1x,3i4)
 2010 FORMAT (/ ' i   :',10i12)
 2020 FORMAT (  ' r   :',10f12.2)
 2030 FORMAT (  ' zeta:',10f12.4)
 2040 FORMAT (  ' f   :',1p10e12.4)
!
      END SUBROUTINE orbit_12

      SUBROUTINE pfit (p, x, y, xv, yv, nx, f, dfdx, dfdy)
!---------------------------------------------------------------------
!     bi-quadratic interpolator
!
!     input
!     1.  p     - effectivly a 4x4 matrix of function values
!                 but for generality the first dimension of
!                 p is nx.
!     2.  x     - associates a grid to the first dimension of p
!     3.  y     - associates a grid to the second dimension of p
!     4.  xv    - location at which bi-quadratic function is evaluated
!     5.  yv    - location at which bi-quadratic function is evaluated
!     6.  nx    - first dimension of p in calling program
!
!     output
!     1.  f     - interpolated value of p at (xv,yv)
!     2.  dfdx  - interpolated value of x-partial of p at (xv,yv)
!     3.  dfdy  - interpolated value of y-partial of p at (xv,yv)
!
! ----------------------------------------------------------------------
!

      IMPLICIT NONE
!
      REAL(DP)  p(nx,*),x(*),y(*),xv,yv,f,dfdx,dfdy,        &
                a1,a2,a3,a4,b1,b2,b3,b4,d2fdx2,d2fdy2
      REAL(DP)  cx(4),cy(4),cxp(4),cyp(4), c2xp(4),c2yp(4)
      INTEGER(I4B) nx,i,j
!
      a1    = (x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))
      a2    = (x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))
      a3    = (x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))
      a4    = (x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))
!
      cx(1) = (xv-x(2))*(xv-x(3))*(xv-x(4))/a1
      cx(2) = (xv-x(1))*(xv-x(3))*(xv-x(4))/a2
      cx(3) = (xv-x(1))*(xv-x(2))*(xv-x(4))/a3
      cx(4) = (xv-x(1))*(xv-x(2))*(xv-x(3))/a4
!
      cxp(1) = ((xv-x(3))*(xv-x(4)) &
         +      (xv-x(2))*(xv-x(3)) &
         +      (xv-x(2))*(xv-x(4)))/a1
      cxp(2) = ((xv-x(3))*(xv-x(4)) &
         +      (xv-x(1))*(xv-x(3)) &
         +      (xv-x(1))*(xv-x(4)))/a2
      cxp(3) = ((xv-x(1))*(xv-x(2)) &
         +      (xv-x(1))*(xv-x(4)) &
         +      (xv-x(2))*(xv-x(4)))/a3
      cxp(4) = ((xv-x(1))*(xv-x(2)) &
         +      (xv-x(1))*(xv-x(3)) &
         +      (xv-x(2))*(xv-x(3)))/a4
!
      c2xp(1) = 2.0 * (3.0*xv-x(2)-x(3)-x(4))/a1
      c2xp(2) = 2.0 * (3.0*xv-x(1)-x(3)-x(4))/a2
      c2xp(3) = 2.0 * (3.0*xv-x(1)-x(2)-x(4))/a3
      c2xp(4) = 2.0 * (3.0*xv-x(1)-x(2)-x(3))/a4
!
      b1      = (y(1)-y(2))*(y(1)-y(3))*(y(1)-y(4))
      b2      = (y(2)-y(1))*(y(2)-y(3))*(y(2)-y(4))
      b3      = (y(3)-y(1))*(y(3)-y(2))*(y(3)-y(4))
      b4      = (y(4)-y(1))*(y(4)-y(2))*(y(4)-y(3))
!
      cy(1)   = (yv-y(2))*(yv-y(3))*(yv-y(4))/b1
      cy(2)   = (yv-y(1))*(yv-y(3))*(yv-y(4))/b2
      cy(3)   = (yv-y(1))*(yv-y(2))*(yv-y(4))/b3
      cy(4)   = (yv-y(1))*(yv-y(2))*(yv-y(3))/b4
!
      cyp(1) = ((yv-y(3))*(yv-y(4)) &
         +      (yv-y(2))*(yv-y(3)) &
         +      (yv-y(2))*(yv-y(4)))/b1
      cyp(2) = ((yv-y(3))*(yv-y(4)) &
         +      (yv-y(1))*(yv-y(3)) &
         +      (yv-y(1))*(yv-y(4)))/b2
      cyp(3) = ((yv-y(1))*(yv-y(2)) &
         +      (yv-y(1))*(yv-y(4)) &
         +      (yv-y(2))*(yv-y(4)))/b3
      cyp(4) = ((yv-y(1))*(yv-y(2)) &
         +      (yv-y(1))*(yv-y(3)) &
         +      (yv-y(2))*(yv-y(3)))/b4
      c2yp(1) = 2.0_DP * (3.0*yv-y(2)-y(3)-y(4))/b1
      c2yp(2) = 2.0_DP * (3.0*yv-y(1)-y(3)-y(4))/b2
      c2yp(3) = 2.0_DP * (3.0*yv-y(1)-y(2)-y(4))/b3
      c2yp(4) = 2.0_DP * (3.0*yv-y(1)-y(2)-y(3))/b4
!
      f      = zeroc
      dfdx   = zeroc
      dfdy   = zeroc
      d2fdx2 = zeroc
      d2fdy2 = zeroc
!
      DO 10 i=1,4
      DO 10 j=1,4
        f      = f + cx(i)*cy(j)*p(i,j)
        dfdx   = dfdx+cxp(i)*cy(j)*p(i,j)
        dfdy   = dfdy+cx(i)*cyp(j)*p(i,j)
        d2fdx2 = d2fdx2+c2xp(i)*cy(j)*p(i,j)
        d2fdy2 = d2fdy2+c2yp(j)*cx(i)*p(i,j)
   10 CONTINUE
      IF (nx .GT. izero)  RETURN
      dfdx = dfdx/d2fdx2
      dfdy = dfdy/d2fdy2
      RETURN
!
      END SUBROUTINE pfit



      SUBROUTINE rotate (naptr,iatype,aheigh,awidth,alen,bhofset, &
                         bvofset,cangv,cangh,ib,isourc,costp,sintp, &
                         costpp,sintpp,blenp,nsourc,sangv,sangh, &
                         rpivot,zpivot,mlost,x0,y0,z0,vx0,vy0,vz0)
!--------------------------------------------------------------------------------
!  this subroutine advances a particle from source to the pivot point,
!  and transforms coordinates.
!
!     translation(s):  Particle to each aperture aligned along the source
!                      axis.  Return with mlost = 1 if particle hits aperture.
!     rotations:       about y-axis so z-axis is aligned vertically;
!                      about new z-axis so x-axis parallel to beamline axis.
!     translation:     coordinates to intersection of beamline axis with
!                      source exit plane.
!     translation:     particles to each of beamline axis centered apertures,
!                      checking for mlost = 1 condition.
!     translation:     particle and then coordinate axis to pivot point.
!     rotation:        x-axis through pivot point and torus center.
!     translation:     origin to torus center.
!---------------------------------------------------------------------------------------

      IMPLICIT NONE
      REAL,SAVE :: iap
!      INTEGER, PARAMETER ::                                nap = 10 ! part of module
      REAL(DP)  aheigh(nap,*), awidth(nap,*), &
                alen(nap,*), bhofset(*), bvofset(*), blenp(*), cangv(*), &
                cangh(*), sangv(*), sangh(*), rpivot(*), zpivot(*), &
                costp(*), sintp(*), costpp(*),  sintpp(*)
      REAL(DP) alen1,alen2,tpvt,x0,y0,z0,vx0,vy0,vz0,zcos,zsin,temp
      INTEGER(I4B) nsourc,naptr,iatype(nap,*),mlost,i,ib,isourc
!
!      IF (nsourc .EQ. 1)  go to 19 added Bob H corection 4/1/2011 HSJ
!
!BH110314:  IF (nsourc .EQ. 1)  go to 19
!BH110314:  Define iap for nsourc=1 case
      IF(nsourc.eq.1)  THEN
         iap=0
         go to 19
      ENDIF
!
!     Move particle to each of source-axis centered apertures,
!     and test for particle passage.
!
      alen1 = zeroc
      alen2 = zeroc
      DO  i=1,naptr
         !
         !     iatype .le. 4 are source-centered apertures
         !
         IF (iatype(i,ib) .GT. 4)  go to 11   ! assumes orderd apperatures ??? HSJ
         alen1 = alen(i,ib)-alen2             ! for DIII-D fist two apertures are iatype =2
         alen2 = alen(i,ib)
         tpvt = -alen1/vx0
         x0 = x0+vx0*tpvt
         y0 = y0+vy0*tpvt
         z0 = z0+vz0*tpvt
         !
         mlost = 0
         go to (12,13,14,15),  iatype(i,ib)
12       IF ((y0**2+z0**2) .LE. (0.5 * awidth(i,ib))**2)  go to 10
         go to 16
13       IF (ABS (y0) .LE. 0.5 * awidth(i,ib) .AND. &
              ABS (z0) .LE. 0.5 * aheigh(i,ib))  go to 10
         go to 16
14       IF (ABS (z0) .LE. 0.5 * aheigh(i,ib))  go to 10
         go to 16
15       IF (ABS (y0) .LE. 0.5 * awidth(i,ib))  go to 10
16       mlost = 1
         RETURN
   10    iap = i
      ENDDO


   11 CONTINUE
!
!     Source center is taken to lie at vertical distance
!     bvofset and horiz. distance bhofset from beam line axis.
!     Source centerline intersects beam line axis at distance
!     bleni from source.  Coords. are changed from source
!     coordinates to those with x-axis along beamline axis
!     directed in +R direction, z in the vertical direction, and
!     origin in source exit plane.
!
!     Rotate about y-axis of source coordinate system so x-axis
!     is vertical and in exit plane defined by beam axis
!     and source centers.
      temp = x0*costp(ib)-isourc*z0*sintp(ib)
      z0 = isourc*x0*sintp(ib)+z0*costp(ib)
      x0 = temp
      temp = vx0*costp(ib)-isourc*vz0*sintp(ib)
      vz0 = isourc*vx0*sintp(ib)+vz0*costp(ib)
      vx0 = temp
!
!     Rotate about z-axis of source system so y-axis lies
!     in source exit plane.
      temp = x0*costpp(ib)-isourc*y0*sintpp(ib)
      y0 = isourc*x0*sintpp(ib)+y0*costpp(ib)
      x0 = temp
      temp = vx0*costpp(ib)-isourc*vy0*sintpp(ib)
      vy0 = isourc*vx0*sintpp(ib)+vy0*costpp(ib)
      vx0 = temp
!
!     Translate coordinate axes to beamline axis.
      x0 = x0
      y0 = y0+isourc*bhofset(ib)
      z0 = z0+isourc*bvofset(ib)
!
!     Translate particle to beamline centered apertures and set
!     mlost = 1 if particle lost.
!
   19 CONTINUE
      IF (iap .EQ. naptr)  go to 29            ! no more apertures to contend with HSJ
      alen1 = 0.0_DP
      alen2 = ABS (x0)
      DO 20 i=iap+1,naptr                      ! loop over beamline ctr aperatures
      alen1 = alen(i,ib)-alen2
      IF (alen1 .LT. 0.0_DP)THEN
        IF(myid == master)WRITE(ncrt,FMT='("Error in sub rotate ")')
        lerrno =  iomaxerr + 226_I4B
        CALL terminate(lerrno,nlog)
      ENDIF
      alen2 = alen(i,ib)
      tpvt = -alen1/vx0                        ! vx0 is always negative 
      x0 = x0+vx0*tpvt                         ! here x0 on rhs is previous aperature 
      y0 = y0+vy0*tpvt                         ! and we advance the particle 
      z0 = z0+vz0*tpvt                         ! to next aperature HSJ
      mlost = 0
      go to (22,23,24,25,26),  iatype(i,ib)-4
   22 IF ((y0**2+z0**2) .LE. (0.5 * awidth(i,ib))**2)  go to 20
      go to 27
   23 IF (ABS (y0) .LE. 0.5 * awidth(i,ib) .AND. &
          ABS (z0) .LE. 0.5 * aheigh(i,ib))  go to 20
      go to 27
   24 IF (ABS (z0) .LE. 0.5 * aheigh(i,ib))  go to 20
      go to 27
   25 IF (ABS (y0) .LE. 0.5 * awidth(i,ib))  go to 20
      go to 27
!     DIII-D  special case:  iatype(i,ib) = 9
   26 IF (ABS (z0) .GT. 24.75)  go to 27
      IF (ABS (y0) .LE. 13.45)  go to 20
      IF (ABS (y0) .GT. 21.7)  go to 27
      IF (ABS (z0) .GT. 24.75*((ABS (y0)-13.45)/(-8.25)+1.))  go to 27
      go to 20
   27 mlost = 1
      RETURN
   20 CONTINUE
   29 CONTINUE
!
!  translate particle and then coordinate axes to pivot point
!
      tpvt = -(blenp(ib)+x0)/vx0
      x0 = 0.0_DP
      y0 = y0+vy0*tpvt
      z0 = z0+vz0*tpvt
!
      IF (sangv(ib) .EQ. 0.0_DP)  go to 30
      zcos = cangv(ib)
      zsin = sangv(ib)
      temp = x0*zcos+z0*zsin
       z0 = -x0*zsin+z0*zcos
      x0 = temp
      temp = vx0*zcos+vz0*zsin
      vz0 = -vx0*zsin+vz0*zcos
      vx0 = temp
   30 CONTINUE
!
      IF (sangh(ib) .EQ. 0.0_DP)  go to 40
      zcos = cangh(ib)
      zsin = sangh(ib)
      temp = x0*zcos+y0*zsin
      y0 = -x0*zsin+y0*zcos
      x0 = temp
      temp = vx0*zcos+vy0*zsin
      vy0 = -vx0*zsin+vy0*zcos
      vx0 = temp
   40 CONTINUE
!
      x0 = x0 + rpivot(ib)
      z0 = z0 + zpivot(ib)
      RETURN
!
      END SUBROUTINE ROTATE


      SUBROUTINE sorspt (nbshape, bheigh, bwidth, bhfoc, bvfoc, &
                         bhdiv, bvdiv, ib, ie, isourc, ke, &
                         nsourc, sfrac1, vbeam, &
                         x0, y0, z0, vx0, vy0, vz0)
!--------------------------------------------------------------------------------
! --- generates a particle at the injector surface with
! --- coordinates and velocities x0,y0,z0,vx0,vy0,vz0
! --- These coordinates are attached to the source center, with the
! --- x-direction perpendicular to the source along the source centerline
!--------------------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER*8 nbshape(*)
      REAL(DP)    bheigh(*), bwidth(*), bhfoc(*), bvfoc(*),    &
                  bhdiv(*), bvdiv(*), sfrac1(*), vbeam(ke,*),  &
                  rt2,x0,y0,z0,vx0,vy0,vz0,zsqu,vdx,     &
                  vdy,vdz,thz,thy,vsqrt
      REAL(SP) seed0
      INTEGER(I4B) isourc,ib,ie,ke,nsourc

      DATA        rt2    /1.414213562_DP/
      DATA        seed0  /0.0_SP/
!
      x0     = 0.0_DP
      isourc = 1
      IF (nsourc .EQ. 1)  go to 10
!
!     two sources
!     sfrac1 is fraction of source current coming from source 1
!     (upper or rightmost source, assuming positive offsets)
!

      IF (random12 (seed0) .GT. sfrac1(ib))  isourc = -1
   10 IF ( nbshape(ib) .NE. 'circ'    )  go to 20

   12 y0   = random12 (seed0) - 0.5

      z0   = random12 (seed0) - 0.5
      zsqu = y0**2 + z0**2                     ! particles must be in  disk of  radius r =1,
                                               ! disk is in x = 0 plane 
      IF (zsqu .GT. 0.25)  go to 12            ! bwidth diameter of  circular source,cm
      y0   = y0*bwidth(ib)                     ! This scales calulation to r = bwidth/2 HSJ
      z0   = z0*bwidth(ib)
      go to 30

   20 y0 = bwidth(ib) * (random12 (seed0) - 0.5) ! (0,y0,z0) is birth point of ion on source

      z0 = bheigh(ib) * (random12 (seed0) - 0.5) ! surface (source is rectangular here) HSJ
!
!  Special coding for DIII-D Long Pulse Sources  
!
   30 IF (nbshape(ib) .NE. 'rect-lps')  go to 32 !For DIII-d nbshape is always 'rect-lps'
      IF (   ABS (z0) .GT.  12.0     )  go to 31 ! hard wired 12.0 is not universal SO DIII-D specific HSJ
!
!  particle on central 2 modules
!
      vdx = -1.0_DP
      vdy =  zeroc
      vdz =  zeroc
      go to 33
!
!  particle on upper or lower modules
!
   31 vdx = -1.0_DP
      vdy =  zeroc
      vdz = -z0 / bvfoc(ib)
      go to 33
!
   32 vdx = -1.0_DP
      vdy = -y0 / bhfoc(ib)
      vdz = -z0 / bvfoc(ib)
!
   33 vsqrt = 1.0 / SQRT (vdx**2+vdy**2+vdz**2)
      vx0   = vbeam(ie,ib)*vdx*vsqrt
      vy0   = vbeam(ie,ib)*vdy*vsqrt
      vz0   = vbeam(ie,ib)*vdz*vsqrt
!
      thz   = ranorm ( ) * bvdiv(ib) / rt2   ! bhdiv,bvdiv in degrees
      thy   = ranorm ( ) * bhdiv(ib) / rt2
      vz0   =  vz0 + thz * Rad_per_deg * vx0
      vy0   =  vy0 + thy * Rad_per_deg * vx0 ! this yields v > vbeam for central sources HSJ ??
      RETURN
!
      END SUBROUTINE sorspt


      SUBROUTINE sorspt_150bm (nbshape, bheigh, bwidth, bhfoc, bvfoc, &
                         bhdiv, bvdiv, ib, ie, isourc, ke, &
                         nsourc, sfrac1, vbeam, &
                         x0, y0, z0, vx0, vy0, vz0)
!--------------------------------------------------------------------------------
! --- generates a particle at the injector surface with
! --- coordinates and velocities x0,y0,z0,vx0,vy0,vz0
! --- These coordinates are attached to the source center, with the
! --- x-direction perpendicular to the source along the source centerline
!--------------------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER*8 nbshape(*)
      REAL(DP)    bheigh(*), bwidth(*), bhfoc(*), bvfoc(*),    &
                  bhdiv(*), bvdiv(*), sfrac1(*), vbeam(ke,*),  &
                  rt2,x0,y0,z0,vx0,vy0,vz0,zsqu,vdx,     &
                  vdy,vdz,thz,thy,vsqrt,height
      REAL(SP) seed0,zr
      INTEGER(I4B) isourc,ib,ie,ke,nsourc

      DATA        rt2    /1.414213562_DP/
      DATA        seed0  /0.0_SP/
!
      x0     = zeroc
      isourc = 1

      IF (nsourc .EQ. 1)  go to 10  ! we use nsourc=2 and control with sfrac instead  for DIII-DHSJ
!
!     two sources
!     sfrac1 is fraction of source current coming from source 1
!     (upper or rightmost source, assuming positive offsets)
!

      IF (random12 (seed0) .GT. sfrac1(ib))  isourc = -1 ! sfrac1 =0 or 1 for DIII-D beamlets
   10 IF ( nbshape(ib) .NE. 'circ'    )  go to 20        ! using rect-lps for DIII-D beamlets

   12 y0   = random12 (seed0) - 0.5

      z0   = random12 (seed0) - 0.5
      zsqu = y0**2 + z0**2                     ! particles must be in  disk of  radius r =1,
                                               ! disk is in x = 0 plane 
      IF (zsqu .GT. 0.25)  go to 12            ! bwidth diameter of  circular source,cm
      y0   = y0*bwidth(ib)                     ! This scales calulation to r = bwidth/2 HSJ
      z0   = z0*bwidth(ib)
      go to 30

      !---------------------------------------------------------------
      ! rect-lps sources pick upper ,mid upper mid lower or lower 
      ! (UP,MU,ML,LO)  based on beam (ib) index
      ! (0,y0,z0) is birth point of ion on source
      !---------------------------------------------------------------
   20 y0 = bwidth(ib) * (random12 (seed0) - 0.5) 

  
      HEIGHT = bheigh(ib)        !!!!! NOTE : bheigh is 12 for 150 source see beamlet table
      zr = height * random12 (seed0)  ! zr in (0.,bheigh)
      IF( ib == 3 .OR. ib == 7) THEN  ! UP (LT or RT) z0 in (12.,24.)
           z0 = zr + 12._DP   

      ELSE IF( ib == 4 .OR. ib == 8  )THEN ! MU (LT or RT) z0 in (0.,12.)
           z0 = zr  ! z0 in (0.,bheigh)

      ELSE IF( ib == 5  .OR. ib == 9 )THEN ! ML(LT or RT) z0 in (0.,-12.)
           z0 = -zr

      ELSE IF( ib == 6 .OR. ib == 10  )THEN ! LO (LT or RT) z0 in (-12.,-24.)
           z0 = -(zr+12._DP)

      ELSE  ! standard beamline case
           z0 = height * (random12 (seed0) -0.5_DP) ! z0 in (-bheigh,bheigh)

      ENDIF

     z0 = zr ! no offset ech source defined separately by nubeam input

!  Special coding for DIII-D Long Pulse Sources  
!
   30 CONTINUE
      IF (nbshape(ib) .NE. 'rect-lps')  go to 32      ! For DIII-d nbshape is always 'rect-lps'
      IF (  ABS (z0) .GT.  12.0 )       go to 31      ! hard wired 12.0 , DIII-D specific HSJ
!
!  particle on central 2 modules
!
      vdx = -1.0_DP
      vdy =  zeroc
      vdz =  zeroc
      go to 33
!
!  particle on upper or lower modules
!
   31 vdx = -1.0_DP
      vdy =  zeroc
      vdz = -z0 / bvfoc(ib)                            ! vz0 will point down for UP and up for LO
      go to 33
!
   32 vdx = -1.0_DP
      vdy = -y0 / bhfoc(ib)
      vdz = -z0 / bvfoc(ib)
!
   33 vsqrt = 1.0 / SQRT (vdx**2+vdy**2+vdz**2)
      vx0   = vbeam(ie,ib)*vdx*vsqrt
      vy0   = vbeam(ie,ib)*vdy*vsqrt
      vz0   = vbeam(ie,ib)*vdz*vsqrt
!
      thz   = ranorm ( ) * bvdiv(ib) / rt2   ! bhdiv,bvdiv in degrees
      thy   = ranorm ( ) * bhdiv(ib) / rt2
      vz0   =  vz0 + thz * Rad_per_deg * vx0
      vy0   =  vy0 + thy * Rad_per_deg * vx0 ! this does not conserve energy (vbeam**2) HSJ ??
      RETURN
!
      END SUBROUTINE sorspt_150bm





      FUNCTION rfine (b1ins, b2ins, c1, c2, c3, finsid, i, ileft, &
                          mfm1, pinsid, rinsid, rmajor)
!-----------------------------------------------------------------------
!  This function uses a fine mesh to calculate:
!   a. the leftmost root of finsid within the interval from rinsid(i)
!      to rinsid(i-1) when ileft = 1 or
!   b. the rightmost root of finsid within the interval from rinsid(i+1)
!      to rinsid(i) when ileft = -1.
!-------------------------------------------------------------------------
      IMPLICIT NONE

      REAL(DP)  b1ins(*), b2ins(*), finsid(*), pinsid(*), rinsid(*)
      REAL(DP)  r1,f1,drin,r2,psi2,b12,b22,f2,rfine,c1,c2,c3,rmajor, &
                drini

      INTEGER(I4B) i,ileft,j,mfm1

      r1 = rinsid(i)
      f1 = finsid(i)
      drin = rinsid(i-1)-rinsid(i)
      drini = 1.0/drin
      DO 10 j=1,10
      r2 = r1 + ileft*0.1*drin
      psi2 = yinter(drini,mfm1,rmajor,r2,pinsid)
      b12 = yinter(drini,mfm1,rmajor,r2,b1ins)
      b22 = yinter(drini,mfm1,rmajor,r2,b2ins)
      f2 = r2**2 - c1*r2*b12 - (c2*b22*(c3+psi2))**2
      IF (f2 .GT. zeroc)  go to 20
      r1 = r2
      f1 = f2
   10 CONTINUE
      rfine = 0.0_DP
      RETURN
   20 rfine = (r1*f2-r2*f1)/(f2-f1)
      RETURN

      END FUNCTION rfine

      SUBROUTINE support (a, ni, nj, nk, nl, nw,nwx2,nh,nhx2,isupp, ifail)
! ----------------------------------------------------------------------
! determines crude 2d support of four-dimensional array
! a(i,j,k,l)   array searched for support approximation
! isupp(4,k,l) output.  with k,l fixed, contains (in order),
!   i1,i2,j1,j2, the indices of a corresponding to the
!   vertices of the smallest rectangle supporting a.
! ----------------------------------------------------------------------
      USE nf_param,                                      ONLY : ke,kb


      IMPLICIT NONE

!      PARAMETER (nxx2 = nw*2, nyx2 = nh*2)
!     argument list:
      INTEGER(I4B) isupp(4,ke,kb),ni,nj,nk,nl,nw,nh,nwx2,nhx2,ifail
      REAL(DP)     a(nwx2,nhx2,ke,kb)


!     local storage
      INTEGER(I4B) i,j
!
      DO 10 i=1,ni
      DO 10 j=1,nj
        IF (a(i,j,nk,nl) .NE. 0)  go to 12
   10 CONTINUE
      go to 99
   12 CONTINUE
      isupp(1,nk,nl) = i
      DO 20 i=ni,1,-1
      DO 20 j=1,nj
      IF (a(i,j,nk,nl) .NE. 0.0)  go to 22
   20 CONTINUE
      go to 99
   22 CONTINUE
      isupp(2,nk,nl) = i
      DO 30 j=1,nj
      DO 30 i=1,ni
      IF (a(i,j,nk,nl) .NE. 0.0_DP)  go to 32
   30 CONTINUE
      go to 99
   32 CONTINUE
      isupp(3,nk,nl) = j
      DO 40 j=nj,1,-1
      DO 40 i=1,ni
      IF (a(i,j,nk,nl) .NE. 0.0_DP)  go to 42
   40 CONTINUE
      go to 99
   42 CONTINUE
      isupp(4,nk,nl) = j
      ifail = 0
      RETURN
!
   99 ifail = 1
      RETURN
!
      END SUBROUTINE support



      SUBROUTINE timtor (rin, rmax, x0, y0, z0, vx0, vy0, vz0, &
                        zmin, zmax, tenter, texit)
! ----------------------------------------------------------------------
!  This subroutine calculates the times for a particle to enter and
!     exit a toroidal box surrounding the plasma starting from the
!     point (x0,y0,z0) and moving with velocity (vx0,vy0,vz0).
! ----------------------------------------------------------------------
!
      IMPLICIT NONE 

      REAL(DP)   edge(4), timsol(6),rin,rmax,x0,y0,z0,vx0,vy0,vz0,&
                 zmin,zmax,tenter,texit,aa,bb,cc,arg,sqr,tt,zz,xx,&
                 yy,rr
                 
      INTEGER(i4B) isol,iin,i
!
!  specify edges of toroidal box surrounding plasma
!
      edge(1) = rin
      edge(2) = rmax
      edge(3) = zmin
      edge(4) = zmax

!
!  find times for particle to intersect inside and outside of box
!
      isol = izero
      aa   = vx0**2+vy0**2
      IF (aa .EQ. 0.0_DP)  go to 50
      bb   = 2.0 * (vx0*x0+vy0*y0)
      DO 40 i=1,2
      cc   = x0**2+y0**2-edge(i)**2
      arg  = bb**2-4.0 * aa*cc
      IF (arg) 40, 30, 10
   10 sqr  = SQRT (arg)
      tt   = (-bb-sqr)/(2.0*aa)
      zz   = z0+vz0*tt
      IF (zz .LT. zmin)  go to 20
      IF (zz .GT. zmax)  go to 20
      isol = isol+1
      timsol(isol) = tt
   20 tt   = (-bb+sqr)/(2.0*aa)
      zz   = z0+vz0*tt
      IF (zz .LT. zmin)  go to 40
      IF (zz .GT. zmax)  go to 40
      isol = isol+1
      timsol(isol) = tt
      go to 40
   30 tt   = -bb/(2.0*aa)
      zz   = z0+vz0*tt
      IF (zz .LT. zmin)  go to 40
      IF (zz .GT. zmax)  go to 40
      isol = isol+1
      timsol(isol) = tt
   40 CONTINUE
!
!  find times for particle to intersect top and bottom of box
!
   50 IF (vz0 .EQ. 0.0_DP)  go to 70
      DO 60 i=3,4
      tt = (edge(i)-z0)/vz0
      xx = x0+vx0*tt
      yy = y0+vy0*tt
      rr = SQRT (xx**2+yy**2)
      IF (rr .LT. rin)  go to 60
      IF (rr .GT. rmax)  go to 60
      isol = isol+1
      timsol(isol) = tt
   60 CONTINUE
   70 CONTINUE
!
!  return if particle misses box
!
      tenter = -1.0e10
      IF (isol .EQ. 0)  RETURN
!
!  calculate times to enter and exit box
!
      tenter = timsol(1)
      iin    = 1
      DO 80 i=2,isol
      IF (timsol(i) .GE. tenter)  go to 80
      iin = i
      tenter = timsol(i)
   80 CONTINUE
      texit = 1.0e10
      DO 90 i=1,isol
      IF (i .EQ. iin)  go to 90
      IF (timsol(i) .LT. texit) texit = timsol(i)
   90 CONTINUE
      RETURN
!
      END SUBROUTINE timtor



      FUNCTION yinter (dxi, nxm1, x0, x, ytab)
! ----------------------------------------------------------------------
! this function performs a fast linear interpolation (or extrapolation)
! for y from x given ytab vs xtab, where xtab starts at x0 and has
! uniform spacing 1/dxi
! ----------------------------------------------------------------------

      IMPLICIT NONE
      REAL(DP) yinter,dxi,x0,x,ytab(*),xx,wt
      INTEGER(I4B) nxm1,i


      xx     = ABS (x-x0)*dxi + 1.0_DP
      i      = xx
      i      = MIN0 (i, nxm1)
      wt     = xx - i
      yinter = (1.0_DP-wt)*ytab(i) + wt*ytab(i+1)
      RETURN

      END FUNCTION yinter

    END MODULE Nfreya_routines


