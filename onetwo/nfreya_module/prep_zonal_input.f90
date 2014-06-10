   MODULE prep_zones

      CONTAINS

      SUBROUTINE prenub  (icall_prenub)
! ----------------------------------------------------------------------
!
! this subroutine calculates certain flux surface information needed by FREYA
! if orbit effects are to be considered, then calculate the
! magnetic flux (pinsid and potsid) vs. major radius (rinsid and rotsid)
! through the magnetic axis.  if orbit effects are not being
! modeled, then simply store the flux and radius values at the
! magnetic axis and limiter in the appropriate arrays.  all flux
! and radius values are in CGS units here.
! --- output is on zonal mf grid---------------------------------------------
!          icall_prenub is icremented by 1 used to set iterative density
!          psif(mf)
!          psivol(mfm1)     volume of flux zones
!          zangrot(mfm1)        angrot on mfm1 grid
!          zne(mfm1)
!          zni(mfm1,i)
!          zte(mfm1)
!          zzi(mfm1,i)
!
! 
!          b1ins(mf)
!          b2ins(mf)
!          b1ots(mf)
!          b2ots(mf)
!          pinsid(mf)
!          potsid(mf)
!          rinsid(mf)
!          rotsid(mf)
!
!   
!          zti(mfm1)
!
! 
!          p(nw,nh)   for 1-1/2d case p is changed to gauss-cm**2 on output
!                     p(nw,nh) is changed back to kgauss-cm**2 in postnub
!
! --------------------------------------------------HSJ 1/14/2011----------
!

      USE nrtype,                    ONLY : DP,I4B,I2B

      USE error_handler

      USE grid_class,                ONLY : nj,psir_grid,reverse

      USE MPI_data,                  ONLY : myid,master,numprocs

      USE io_gcnmp,                  ONLY : nlog,ncrt
    
      USE neutral_beams,             ONLY : iborb,nbion,fdbeam,nameb,ibion,         &
                                            enbeam_tot,fast_ion_target,elong,       &
                                            nfreya_plot_unit
                                            ! NOTE: fd_beam is defined in statefile read
                                            ! fdbeam is namelist input (which is used here)
                                            ! nameb(1:nbion) defined in statefile and 
                                            ! in namelist input

      USE ions_gcnmp,                ONLY : h_index ,t_index,he_index,                      &
                                            dt_index,d_index,fi_index

      USE vector_class,              ONLY : get_values

      USE common_constants,          ONLY : vs2kgcm2,zeroc,izero,M2cm,im32icm3,M32cm3,      &
                                            Tm2gcm

      USE Plasma_properties ,        ONLY : mhd_dat,dischg,profile

      USE tension_spline,            ONLY : spline_intrp

      USE contour,                   ONLY : cntour,cntour1,cntour2,getrmaj,                  &
                                            fixedcntour,volcalc,limiter_check

      USE zonal_data,                ONLY :      wnoperm,zenbeam,zenbeamold,                 &
                                                 zone_volume,mf,mfm1,                        &
                                                 zone_area,potsid,psif,                      &
                                                 zone_area,rinsid,psivol,             &
                                                 rotsid,pinsid,b1ins,b2ins,                  &
                                                 b1ots,b2ots,fpsio,fpsii,zne,                &
                                                 n_zone_cntr,r_zone_cntr,z_zone_cntr,        &
                                                 zte,zti,zangrot,zni,zzi,nw,nh,nwh,nh2,      &
                                                 max_ctr_pts
                                                 
   
      USE tport_mhd_grid_data,       ONLY :      psir,  psival,rmhdgrid, zmhdgrid,psivolp,   &
                                                 psi_g_cm2,en,xlimiter,ylimiter,rplasbdry,   &
                                                 zplasbdry, &
                                                 bpmag,rcontour,zcontour,fpsi,ene,te,ti,     &
                                                 angrot,xmagn1,ymagn1,rcmin,rcmax,zcmin,     &
                                                 zcmax,zrcmin,zrcmax,rzcmin,rzcmax,volume,   &
                                                 xmin_lim,xmax_lim,ymin_lim,ymax_lim,        &
                                                 psiax,psilim,ncontour,                      &
                                                 rplasmin,rplasmax,zplasmin,zplasmax

      USE ions_gcnmp,                ONLY :      namep,namei,nion,                           &
                                                 namen,nprim,nimp,nneu,z,zsq,zeff,           &
                                                 name_size,fd_thermal

      USE replace_imsl,              ONLY :      my_ibcccu,my_dbcevl1

      IMPLICIT  NONE  
 
      REAL(DP),DIMENSION(:,:,:),ALLOCATABLE   :: cspln
      REAL(DP) pds(6)

      INTEGER,PARAMETER  ::  n2cspln = 2
      INTEGER(I4B) icall_prenub,nwork,nplasbdry,     &
                   nlimiter,npts,iz
 

 

      LOGICAL get_scoef

      REAL(DP) drutp,dsp,d2sp,sdl,sdr,cconst,csgn

      INTEGER(I4B) i,j,asize ,nmx,sier,iconvg,iauto,ierr,iautoc,  &
                   isigncur,npsi,ib,maxpts,task,error

      REAL(DP),ALLOCATABLE,DIMENSION(:) :: work,cparam,zdum

      REAL(DP) bperr,arcl, taxis,tlim,a,bincp,delta_psi,          &
               psi_psif,ptrace,dang,drx,dry,drot,drin, brelax,    &
               absdifa,absdifb

      REAL(DP),ALLOCATABLE,DIMENSION(:) :: enbeams
     INTERFACE
         ! utils.f90:
         SUBROUTINE to_lower_case(string) 
           USE nrtype,            ONLY : I4B
           IMPLICIT NONE
           INTEGER(I4B) l
           CHARACTER*(*), INTENT (INOUT) :: string
         END SUBROUTINE to_lower_case

     END INTERFACE


     error = izero


! ----------------------------------------------------------------------
! calculate psif, an array of flux zone boundaries used in FREYA.
!     for a 1-1/2-D run the zones are of equal width in the square
!     root of the flux.  psif is in units of 
!     kgauss-cm**2 for 1-1/2-D runs.
!     NOTE: square root makes flux zones near axis too small to get
!           accurate contours. square of flux makes zones near axis
!           too large. Therefore the best choice still is linear in psi.
!           the dropoff in the density and increase in volume near the
!           plasma boundary approximately lead to constant deposition.
!     NOTE: the square root spacing is also used in FREYA so if it
!           is changed here it must be changed in FREYA also!
!     ALSO: The same is true for subroutine INJECT.
! ----------------------------------------------------------------------
!

          IF(.NOT. ALLOCATED(psif))ALLOCATE(psif(mf))
          IF(.NOT.  ALLOCATED(psir))ALLOCATE(psir(nj))
          psir(:)    =  get_values(psir_grid)
          psir(nj)   =  mhd_dat%psibdry      ! some difference here
          psir(:)    =  psir(:) * vs2kgcm2   ! volt sec  to kgauss cm^2
          drutp = SQRT (psir(nj)-psir(1))/mfm1
          DO  i=1,mf
            psif(i) = psir(1) + ((i-1)*drutp)**2
          ENDDO
! optional linear spacing:
!         dpsi = (psir(nj)-psir(1))/mfm1
!         DO i=1,mf
!            psif(i) = psir(1)+(i-1)*dpsi
!         ENDDO
          psif(mf)   = psir(nj)  ! avoid roundoff psif(1) = axis,psif(mf) = edge

!
! ----------------------------------------------------------------------
!  calculate volumes of flux zones for nonelliptical cross sections:
!           flux is in units of kgauss-cm**2; dimensions are in cm;
! ----------------------------------------------------------------------
!
          IF(ALLOCATED(work))THEN
             asize = SIZE(work)
             IF(asize .NE. mhd_dat%npsi)DEALLOCATE(work)
          ENDIF
          IF(.NOT. ALLOCATED(work))ALLOCATE(work(mhd_dat%npsi))

 

          asize = 0 
          IF(ALLOCATED(psival))DEALLOCATE(psival)
          ALLOCATE(psival(mhd_dat%npsi))
          IF(ALLOCATED(psivolp))DEALLOCATE(psivolp)
          ALLOCATE(psivolp(mhd_dat%npsi))
          npsi = mhd_dat%npsi
          psival(:)  = get_VALUES(mhd_dat%psivalnpsi)

          psival(:)  = psival(:)* vs2kgcm2   ! volt sec  to kgauss cm^2
          psivolp(:) = get_VALUES(dischg%psivolpnpsi)
          psivolp(:) = psivolp(:)*M32cm3
          ! NOTE psival(1) = plasma edge, psival(npsi) plasma center
          ! NOTE psivolp(1) = plasma edge, psivolp(npsi) plasma center
          IF(ALLOCATED(zdum))DEALLOCATE(zdum)
          ALLOCATE(zdum(mf))

 
          CALL reverse (npsi,psival)
          CALL reverse ( npsi,psivolp)  
         ! psival,psivolp,psif  are now numbered from axis to edge


          get_scoef = .TRUE.                   ! determine cparam on call to spline_intrp
          IF(ALLOCATED(cparam))DEALLOCATE(cparam)
          ALLOCATE(cparam(mhd_dat%npsi))
          ib =2 ! supply second derivatives
          sdl  =zeroc ; sdr = zeroc
          ! set cparam, psif(1),zdum(1) are not accessed
          CALL spline_intrp(npsi,psival,psivolp,cparam,&
                             get_scoef,ib,sdl,sdr, psif(1),zdum(1),dsp,d2sp) 
          get_scoef = .FALSE.                ! cparam is now set
          DO j=1,mf
              CALL spline_intrp(npsi,psival,psivolp,cparam,&
                   get_scoef,ib,sdl,sdr, psif(j),zdum(j),dsp,d2sp)
          ENDDO


          IF(ALLOCATED(psivol))DEALLOCATE(psivol)
          ALLOCATE(psivol(mfm1))
          DO i=1,mfm1
             psivol(i) = zdum(i+1)-zdum(i)  ! zdum(1) =0,zdum(mf) = total plasma volume
          ENDDO

!
! ----------------------------------------------------------------------
! psi_g_cm2 - poloidal flux array as is demanded by Nfreya code
! convert flux values to gauss-cm**2 in array p 
! mhd_dat%psi is in volt sec /rad 
! ----------------------------------------------------------------------
!

          IF(ALLOCATED(psi_g_cm2))THEN
             nw = SIZE(psi_g_cm2,1)
             nh = SIZE(psi_g_cm2,2)
             IF(nw .NE.dischg%nr_mhd .OR. nh .NE. dischg%nz_mhd)DEALLOCATE(psi_g_cm2)
          ENDIF
          IF( .NOT. ALLOCATED(psi_g_cm2))ALLOCATE(psi_g_cm2(dischg%nr_mhd,dischg%nz_mhd))

          psi_g_cm2(:,:) = mhd_dat%psi(:,:)*1.e3_DP * vs2kgcm2           ! gauss-cm**2
 
        

 
          psiax  = 1.0e3_DP*psif( 1)                        ! pisf is in kg-cm2
          psilim = 1.0e3_DP*psif(mf)
          absdifa  = DABS(mhd_dat%psiaxis*1.e3_DP * vs2kgcm2  - psiax)
          absdifb  = DABS(mhd_dat%psibdry*1.e3_DP * vs2kgcm2  - psilim)
          IF( absdifa .GT.   1.e-8_DP .OR. absdifb  .GT. 1.e-8_DP )THEN

             IF(myid == master)THEN
                 WRITE(ncrt,FMT='( "sub prenub: absdifa ,absdifb =",2(1pe12.4,x))')absdifa,absdifb
                 WRITE(ncrt,FMT='( "sub prenub: Error in psi axis/ psi lim values")')
                 WRITE(ncrt,FMT='("psiax,psilim = ",2(1pe24.15,x))'),psiax,psilim
                 WRITE(ncrt,FMT='("mhd_dat%(psiax,psilim) = ",2(1pe24.15,x))') &
                    mhd_dat%psiaxis*1.e3_DP * vs2kgcm2,mhd_dat%psibdry*1.e3_DP * vs2kgcm2
             ENDIF
             lerrno = iomaxerr + 205_I4B
             CALL terminate(lerrno,nlog)
          ENDIF


!
! ----------------------------------------------------------------------
! get bicubic representation of psi ( = psi_g_cm2):
! (reverse the sign of psi so that CNTOUR will work)
! ----------------------------------------------------------------------
!
          IF(ALLOCATED(rmhdgrid))THEN
             nw = SIZE(rmhdgrid)
             IF(nw .NE. dischg%nr_mhd )DEALLOCATE(rmhdgrid)
          ENDIF
          IF( .NOT. ALLOCATED(rmhdgrid))ALLOCATE(rmhdgrid(dischg%nr_mhd))
          IF(ALLOCATED(zmhdgrid))THEN
             nh = SIZE(zmhdgrid)
             IF(nh .NE. dischg%nz_mhd )DEALLOCATE(zmhdgrid)
          ENDIF
          IF( .NOT. ALLOCATED(zmhdgrid))ALLOCATE(zmhdgrid(dischg%nz_mhd))

          nwh = dischg%nr_mhd * dischg%nz_mhd
          nw = dischg%nr_mhd ; nh = dischg%nz_mhd
          nh2 = 2*nh
          nmx     = MAX(nw,nh)
          nwork   = 2*nwh+2*nmx        ! use this if nw > nh

          IF(ALLOCATED(cspln))DEALLOCATE(cspln)
          ALLOCATE(cspln(n2cspln,nw,nh2))
          IF(ALLOCATED(wnoperm))DEALLOCATE(wnoperm)
          ALLOCATE(wnoperm(nwork))
          rmhdgrid(:)  =    get_values( dischg%rmhdgrid)* M2cm 
          zmhdgrid(:)  =    get_values( dischg%zmhdgrid)* M2cm 
 
          cconst = -1.0_DP
          psi_g_cm2(:,:) = psi_g_cm2(:,:)*cconst
!          CALL ibcccu  (psi_g_cm2,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,wnoperm,sier)
          CALL my_ibcccu  (psi_g_cm2,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,wnoperm,sier)
          psi_g_cm2(:,:) = psi_g_cm2(:,:)*cconst


!-----------------------------------------------------------------------
! -- load limiter contour from statefile
!-----------------------------------------------------------------------
          nlimiter  =  dischg%nlimiter
          IF(ALLOCATED(xlimiter))DEALLOCATE(xlimiter)
          IF(ALLOCATED(ylimiter))DEALLOCATE(ylimiter)
          ALLOCATE(xlimiter(nlimiter+2),ylimiter(nlimiter+2))
          xlimiter(1:SIZE(xlimiter)) = zeroc
          ylimiter(1:SIZE(ylimiter)) = zeroc
          XLIMITER(NLIMITER+1)  = 1.E100_dp  ! THIS IS FOR SUB limiter_check
          yLIMITER(NLIMITER+1)  = 1.E100_dp
          XLIMITER(NLIMITER+2)  = -1.E100_dp
          yLIMITER(NLIMITER+2)  = -1.E100_dp
          xlimiter(1:nlimiter)  = get_values(dischg%rlimiter)
          ylimiter(1:nlimiter)  = get_values(dischg%zlimiter)
          xlimiter(1:nlimiter)           = xlimiter(1:nlimiter)* M2cm
          ylimiter(1:nlimiter)           = ylimiter(1:nlimiter)* M2cm
          xmin_lim = MINVAL(xlimiter(1:nlimiter)) ; xmax_lim = MAXVAL(xlimiter(1:nlimiter))
          ymin_lim = MINVAL(ylimiter(1:nlimiter)) ; ymax_lim = MAXVAL(ylimiter(1:nlimiter))


!-----------------------------------------------------------------------
! -- set some arrays to store contour point information
!-----------------------------------------------------------------------
          nplasbdry   =  dischg%nplasbdry 
          max_ctr_pts = 3*nplasbdry
          maxpts      = max_ctr_pts

          IF(ALLOCATED(bpmag))DEALLOCATE(bpmag)
          IF(ALLOCATED(zplasbdry))DEALLOCATE(zplasbdry) ! psilim contour
          IF(ALLOCATED(rplasbdry))DEALLOCATE(rplasbdry)
          ALLOCATE(bpmag(maxpts)) ; bpmag(:) = zeroc
          ALLOCATE(zplasbdry(maxpts)) ; zplasbdry(1:maxpts) = zeroc
          ALLOCATE(rplasbdry(maxpts)) ; rplasbdry(1:maxpts) = zeroc
          rplasbdry(1:nplasbdry) = get_values(dischg%rplasbdry )
          zplasbdry(1:nplasbdry) = get_values(dischg%zplasbdry )
          rplasbdry(1:nplasbdry)           = rplasbdry(1:nplasbdry)*M2cm
          zplasbdry(1:nplasbdry)           = zplasbdry(1:nplasbdry)*M2cm

          IF(ALLOCATED(zcontour))DEALLOCATE(zcontour) !interior psi  contour r,z,points
          IF(ALLOCATED(rcontour))DEALLOCATE(rcontour)
          ALLOCATE(zcontour(maxpts))
          ALLOCATE(rcontour(maxpts))
          IF(ALLOCATED(zone_volume))DEALLOCATE(zone_volume)
          IF(ALLOCATED(zone_area))DEALLOCATE(zone_area)
          ALLOCATE(zone_volume(mf))
          ALLOCATE(zone_area(mf))

          IF(ALLOCATED(n_zone_cntr))DEALLOCATE(n_zone_cntr)
          IF(ALLOCATED(r_zone_cntr))DEALLOCATE(r_zone_cntr)
          IF(ALLOCATED(z_zone_cntr))DEALLOCATE(z_zone_cntr)  
          ALLOCATE(n_zone_cntr(mf))
          n_zone_cntr(:) = izero
          ALLOCATE(r_zone_cntr(maxpts,mf))
          ALLOCATE(z_zone_cntr(maxpts,mf))

! ----------------------------------------------------------------------
! generate the psi contours corresponding to the FREYA grid. these are
! copied into the neutral beam plotting file in subroutine NUBPLT.
! ----------------------------------------------------------------------
!

          xmagn1 =  dischg%rma * M2cm
          ymagn1 =  dischg%zma * M2cm


          bperr     = 0.05_DP
          iconvg    = izero
          arcl      = 2.0_DP    ! 2 cm arclength increment
          taxis     = 5.0_DP
          tlim      = 30.0_DP
          a         = (tlim-taxis)/(1000.0_DP*(psif(1)-psif(mf)))
          bincp     = taxis
          delta_psi = 1000.0_DP * (-psif(mf-1)+psif(mf))
          iz        = izero

          DO j=1,mfm1

             i = mf-j+1
             psi_psif = psif(i)

             ptrace = -psi_psif * 1000.0_DP
!
! --- prevent searching for plasma boundary by moving psilim in slightly:
!

             IF (j .EQ. 1)  ptrace = ptrace +0.01_DP*delta_psi
!
             iauto  = izero
             iconvg = izero
             dang   = a*(ptrace+1000.0_DP*psif(1))+bincp
             drx    = zeroc
             dry    = zeroc
             IF (j .GT. 1 )THEN 

                 CALL cntour (xmagn1,ymagn1,ptrace,rcmin,rcmax,zcmin,         &
                          zcmax,zrcmin,zrcmax,rzcmin,rzcmax,dang,arcl,        &
                          bperr,drx,dry,xmin_lim,xmax_lim,ymin_lim,ymax_lim,  &
                          iauto,iautoc,rcontour,zcontour,ncontour,rmhdgrid,nw,&
                          zmhdgrid,nh,cspln,n2cspln,nh2,ncrt,max_ctr_pts,     &
                          ierr,bpmag,iconvg,delta_psi)

                 IF (ierr .NE. 0)THEN
                    IF(myid == master)     &
                      WRITE(ncrt,FMT='("subroutine PRENUB: returned CNTOUR error")')
                    lerrno = 213  + iomaxerr
                    CALL terminate(lerrno,nlog)
                 ENDIF

              ELSE

                 !get rcmin,rcmax,zcmin,zcmax,rzcmin,rzcmax,zrcmin,zrcmax
                 !bpmag from given plasma boundary:

                 CALL fixedcntour (rplasbdry,zplasbdry,nplasbdry, &
                              rcontour,zcontour,ncontour,         &
                              rcmin,rcmax,zcmin,zcmax,            &
                              rzcmin,rzcmax,zrcmin,zrcmax,        &
                              rmhdgrid,zmhdgrid,nw,nh,            &
                              bpmag,cspln,n2cspln,nh2,pds)

                 CALL limiter_check(rcmin,rcmax,zcmin,zcmax,       &
                                    xlimiter,ylimiter,nlimiter)      ! contour_module.f90


                 rplasmin = rcmin
                 rplasmax = rcmax
                 zplasmin = zcmin
                 zplasmax = zcmax
             ENDIF

             CALL volcalc (rcontour,zcontour,ncontour,xmagn1,ymagn1,zone_volume(j),zone_area(j))
!            save contours for use  in plotting
             iz = iz + 1  
             n_zone_cntr(iz) = ncontour 
             r_zone_cntr(1:n_zone_cntr(iz),iz) = rcontour(1:n_zone_cntr(iz))
             z_zone_cntr(1:n_zone_cntr(iz),iz) = zcontour(1:n_zone_cntr(iz))
          END DO

          zone_volume(mf)   = zeroc     ! zone_volume(mf) = volume at magnetic axis =0.0
          volume            = zeroc
          DO j=1,mfm1
            psivol(mfm1-j+1) = zone_volume(j)-zone_volume(j+1)
            volume = volume+psivol(mfm1-j+1)
          END DO


!
! ----------------------------------------------------------------------
! get psi and grad psi for these grids. IBCCCU sets up bicubic
! spline array cspln, for psi values given in array p:
! ----------------------------------------------------------------------
!
          IF(ALLOCATED(potsid))DEALLOCATE(potsid)
          ALLOCATE(potsid(mf))

          IF(ALLOCATED(pinsid))DEALLOCATE(pinsid)
          ALLOCATE(pinsid(mf))

          IF(ALLOCATED(b1ins))DEALLOCATE(b1ins)
          ALLOCATE(b1ins(mf))

          IF(ALLOCATED(b1ots))DEALLOCATE(b1ots)
          ALLOCATE(b1ots(mf))

          IF(ALLOCATED(b2ins))DEALLOCATE(b2ins)
          ALLOCATE(b2ins(mf))

          IF(ALLOCATED(b2ots))DEALLOCATE(b2ots)
          ALLOCATE(b2ots(mf))

          IF(ALLOCATED(rinsid))DEALLOCATE(rinsid)
          ALLOCATE(rinsid(mf))

          IF(ALLOCATED(rotsid))DEALLOCATE(rotsid)
          ALLOCATE(rotsid(mf))

          potsid(1)  = psiax
          potsid(mf) = psilim

          IF (iborb .NE. 0) THEN
!              CALL ibcccu (psi_g_cm2,rmhdgrid,nw,zmhdgrid,nh,cspln,nw, &
!                           wnoperm,sier)
              CALL my_ibcccu (psi_g_cm2,rmhdgrid,nw,zmhdgrid,nh,cspln,nw, &
                           wnoperm,sier)
!------------------------------------------------------------------------------
! --- first get inside and outside major radius at elevation of magnetic axis
!------------------------------------------------------------------------------
              isigncur = -1
              IF (psiax .GT. psilim)  isigncur = 1
              npts = 3
              IF(ALLOCATED(work))DEALLOCATE(work)
              ALLOCATE(work(2*(npts+1)))
              CALL getrmaj (ymagn1, npts, ierr, work(1), work(npts+1), &
                            isigncur, psilim,cspln,n2cspln,nw,nh,nh2,pds)
              rinsid(mf) = work(1)
              rotsid(mf) = work(npts)
!------------------------------------------------------------------------------
! --- next get inboard and outboard (uniform) spatial grids
! --- at magnetic axis elevation and dpsi/dR
!------------------------------------------------------------------------------
              drot = (rotsid(mf)-xmagn1)/mfm1
              drin = (xmagn1-rinsid(mf))/mfm1
              DO  i=1,mfm1
                rotsid(i) = xmagn1+(i-1)*drot
                rinsid(i) = xmagn1+(1-i)*drin
              ENDDO
              DO  i=1,mf
                CALL my_dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,nw, &
                             rotsid(i),ymagn1,pds,sier,2)
                potsid(i) = pds(1)     ! psi     at (rotsid(i),ymagn1)
                b2ots(i) = pds(2)      ! dpsi/dr at (rotsid(i),ymagn1)
                CALL my_dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,nw, &
                             rinsid(i),ymagn1,pds,sier,2)
                pinsid(i) = pds(1)
                b2ins(i) = pds(2) ! note this is just dpsi/dR at this point
              ENDDO
              b2ins(1) = zeroc
              b2ots(1) = zeroc

!------------------------------------------------------------------------------
! --- given f(psi) ( = fpsi) on the psival grid,get f(psi) on the rotsid
! --- ( = fpsio) and rinsid (=fpsii) grids. The MHD grids are numbered
! --- from the plasma edge inward so we have to temporarily reverse them:
!------------------------------------------------------------------------------
              IF(ALLOCATED(fpsi))DEALLOCATE(fpsi)
              ALLOCATE(fpsi(npsi))
              fpsi(:) = get_values(mhd_dat%fpsinpsi)

              IF(ALLOCATED(fpsii))DEALLOCATE(fpsii)
              ALLOCATE(fpsii(mf))
              IF(ALLOCATED(fpsio))DEALLOCATE(fpsio)
              ALLOCATE(fpsio(mf))

              DO  i=1,npsi
                 psival(i) = 1.0e3_DP*psival(i)
              ENDDO

!             psival(1) = axis, psival(npsi) = edge, in gauss cm^2

              CALL reverse(npsi,fpsi)
              fpsi(:) = fpsi(:)/Tm2gcm             ! fpsi in gauss cm
!             fpsi(1) = axis ,fpsi(npsi) = edge = R0*Bt0  

!              CALL intrp (1,1,psival,fpsi,npsi,potsid,fpsio,mf)

               get_scoef = .TRUE.       ! true on first call to spline_intrp
               IF(ALLOCATED(cparam))DEALLOCATE(cparam)
               ALLOCATE(cparam(mhd_dat%npsi))
               ib =2 ! supply second derivatives,sdl,sdr
               sdl  =zeroc ; sdr = zeroc
               CALL spline_intrp(npsi,psival,fpsi,cparam,&
                                  get_scoef,ib,sdl,sdr, potsid(1),fpsio(1),dsp,d2sp)
               get_scoef = .FALSE.    ! cparam is now set
               DO j=1,mf
                   CALL spline_intrp(npsi,psival,fpsi,cparam,&
                   get_scoef,ib,sdl,sdr, potsid(j),fpsio(j),dsp,d2sp)
               ENDDO

!               CALL intrp (1,1,psival,fpsi,npsi,pinsid,fpsii,mf)
               get_scoef = .TRUE.       ! true on first call to spline_intrp
               IF(ALLOCATED(cparam))DEALLOCATE(cparam)
               ALLOCATE(cparam(mhd_dat%npsi))
               ib =2 ! supply second derivatives
               sdl  =zeroc ; sdr = zeroc
               CALL spline_intrp(npsi,psival,fpsi,cparam,&
                                 get_scoef,ib,sdl,sdr, pinsid(1),fpsii(1),dsp,d2sp)
               get_scoef = .FALSE.                ! cparam is now set
               DO j=1,mf
                   CALL spline_intrp(npsi,psival,fpsi,cparam,&
                                     get_scoef,ib,sdl,sdr, pinsid(j),fpsii(j),dsp,d2sp)
               ENDDO


              DO  i=1,npsi
                 psival(i) = 1.0e-3_DP*psival(i)          ! restore psival to kg cm^2
              ENDDO

              CALL reverse(npsi,psival)

              CALL reverse(npsi,fpsi)
!
! --- next define the magnetic field ratios required in orbit:
!
              DO i=1,mf
                 b1ins(i) = SQRT (fpsii(i)**2+b2ins(i)**2)
                 b1ots(i) = SQRT (fpsio(i)**2+b2ots(i)**2)
                 b2ins(i) = b1ins(i) / ABS (fpsii(i))
                 b2ots(i) = b1ots(i) / ABS (fpsio(i))
              ENDDO
              DO i=mf,1,-1
                 b1ins(i) = b1ins(i)/b1ins(1)
                 b1ots(i) = b1ots(i)/b1ots(1)
              ENDDO
          END IF ! iborb = 1 section 



!
! ----------------------------------------------------------------------
! convert densities, temperatures, and angular rotation speed from
!   transport mesh points to FREYA zones for 1-d and 1-1/2d cases:
! ----------------------------------------------------------------------
!
       IF(ALLOCATED(ene))DEALLOCATE(ene)
       ALLOCATE(ene(nj))
       IF(ALLOCATED(zne))DEALLOCATE(zne)
       ALLOCATE(zne(mf))
       ene(:)      = get_values(profile%ene)*im32icm3


       CALL tozone(psir,ene,nj,psif,zne,mfm1) 

       IF(ALLOCATED(te))DEALLOCATE(te)
       ALLOCATE(te(nj))
       te(:) = get_values(profile%te)
       IF(ALLOCATED(zte))DEALLOCATE(zte)
       ALLOCATE(zte(mf))
       CALL tozone(psir,te,nj,psif,zte,mfm1)

       IF(ALLOCATED(ti))DEALLOCATE(ti)
       ALLOCATE(ti(nj))
       ti(:) = get_values(profile%ti)
       IF(ALLOCATED(zti))DEALLOCATE(zti)
       ALLOCATE(zti(mf))
       CALL tozone(psir,ti,nj,psif,zti,mfm1)

 
       IF(ALLOCATED(angrot))DEALLOCATE(angrot)
       ALLOCATE(angrot(nj))
       IF(ALLOCATED(zangrot))DEALLOCATE(zangrot)
       ALLOCATE(zangrot(mf))
       angrot(:) = get_values(profile%angrot)
       CALL tozone(psir,angrot,nj,psif,zangrot,mfm1)
       zangrot(mf) = zeroc  ! only  1:mfm1 is set above
       ! set zangrot(mf) = 0 arbitrarily (not used)
      IF(ALLOCATED(en))DEALLOCATE(en)
      ALLOCATE(en(nj,nion))
      IF(ALLOCATED(zni))DEALLOCATE(zni)
      ALLOCATE(zni(mfm1,nion))
      IF(ALLOCATED(zzi))DEALLOCATE(zzi)
      ALLOCATE(zzi(mfm1,nion+1)) ! nion+1,not nion see below

      DO i =1,nion
         en(:,i) = get_values(profile%en(i))*im32icm3
         ! z(:,i),zsq(:,i) were  set in statefile read
      ENDDO


      IF(ALLOCATED(work))DEALLOCATE(work)
      ALLOCATE(work(nj))
      DO  i=1,nion
         CALL tozone(psir,en(1,i),nj,psif,zni(1,i),mfm1)
         work(:) = z(:,i) ! z is pointer need interface or this
         CALL tozone(psir,work, nj,psif,zzi(1,i),mfm1)
      ENDDO

!-------------------------------------------------------------------------------------
! -- Assume that the stored fast ion density presents a target with
! -- similar cross sections as the corresponding thermal density
! -- add beam density to thermal ion density
! -- this mimicks the effect of the beam seeing circulating fast ions as part
! -- of the target for injected neutrals. Obviously this involves lots of assumptions
! -- zenbeam may be zero
!-------------------------------------------------------------------------------------



      IF (fast_ion_target .GT. 0) THEN
         IF(.NOT. ALLOCATED(zenbeamold))ALLOCATE(zenbeamold(mfm1))
         IF(ALLOCATED(zenbeam))DEALLOCATE(zenbeam)
         ALLOCATE(zenbeam(mfm1))
         IF(ALLOCATED(enbeams))DEALLOCATE(enbeams)          ! local to prenub
         ALLOCATE(enbeams(nj))
         enbeams(:) = enbeam_tot(:)*im32icm3   ! enbeam_tot is defined in statefile input
         CALL tozone(psir,enbeams,nj,psif,zenbeam,mfm1)

!
!        add the beam density to the thermal density
!        that will be used by FREYA
!
         brelax = 0.5_DP
         IF (icall_prenub .GT. 1) THEN
            DO i=1,mfm1
               zenbeam(i)    = brelax*zenbeam(i)+(1.-brelax)*zenbeamold(i)
               zenbeamold(i) = zenbeam(i)
            END DO
         ELSE
           zenbeamold(i) = zenbeam(i) ! on first call zenbeam is laoded from
                                      ! statefile enbeam_tot
                                      ! if beams are not on then  enbeam_tot =0
         END IF


         IF (ibion .GT. 0) THEN         !logic assumes one beam species
           zni(:,ibion) =  zenbeam(:) + zni(:,ibion)
         ELSE ! if dt mixture
            zni(:,d_index) = zni(:,d_index) + fdbeam * zenbeam(:)
            zni(:,t_index) = zni(:,t_index) + (1._DP - fdbeam)*zenbeam(:)

         END IF
      END IF

! ----------------------------------------------------------------------
! -- Put the Zeff array in an extra row of the zzi flux-zone array 
! ----------------------------------------------------------------------
      ! zeff(1:nj) direct from statefile
      CALL tozone (psir,zeff,nj,psif,zzi(1,nion+1),mfm1)



! ----------------------------------------------------------------------
! transfrom to FREYA (lab) coordinate system.
! psi is assumed to correspond to plasma current in the positive
! (FREYA) toroidal system, so carry the sign in var:csgn.
! this needs to be checked out HSJ ?? 
! ----------------------------------------------------------------------
         csgn = 1._DP
         IF(mhd_dat%tot_cur .LT. zeroc)csgn = -1._DP

         zangrot(:) = csgn*zangrot(:)

         elong = dischg%kappa

    IF(ALLOCATED(work))DEALLOCATE(work)
    RETURN
!
    END SUBROUTINE prenub 



      SUBROUTINE tozone (x, y, n, xz, yz, nz)
! ----------------------------------------------------------------------
! convert from point values to zone values
! INPUT
!   x(1..n)      abscissa values
!   y(1..n)      ordinate
!        n       # pys
!   xz(1..nz)    abscissa values of zones
!    nz          # zones 
! OUTPUT
!    yz(1..nz)   nz          zonal vlaues
! ----------------------------------------------------------------------
!
      USE nrtype,                             ONLY : DP,I4B

      USE error_handler

      USE io_gcnmp,                           ONLY : nlog,ncrt

      IMPLICIT  NONE 
    
      INTEGER(I4B)  ip,iz,nz,ipm1,n
      REAL(DP) tol ,yzlast,del,s

      REAL(DP)  x(*), y(*), xz(*), yz(*),reldif
      DATA       tol/1.0e-16_DP/
!
      ip     = 1
      iz     = 1
      yzlast = y(n)
!
! ----------------------------------------------------------------------
! interpolate to new grid that is compute yz(xz(iz))
! ----------------------------------------------------------------------
!
 2000 IF (iz .GT. nz)  go to 2210
 2100 IF (ip .GT. n )  go to 9200
      del    = reldif(x(ip),xz(iz))
      IF (del .GT. tol)  go to 2110
      yz(iz) = y(ip)
      iz     = iz + 1
      go to 2000
 2110 IF (x(ip) .GT. xz(iz))  go to 2120
      ip   = ip + 1
      go to 2100
 2120 ipm1   = ip - 1
      s      = (y(ip)-y(ipm1))/(x(ip)-x(ipm1))
      yz(iz) = y(ipm1)+(xz(iz)-x(ipm1))*s
      iz     = iz + 1
      go to 2000
!
! ----------------------------------------------------------------------
! convert to zone
! ----------------------------------------------------------------------
!
 2210 DO iz=1,nz-1
        yz(iz) = 0.5 * (yz(iz)+yz(iz+1))
      END DO
      yz(nz) = 0.5 * (yz(nz)+yzlast)
      RETURN
!
! ----------------------------------------------------------------------
! fatal errors
! ----------------------------------------------------------------------
!
 9200 WRITE (nlog, 8000)
      WRITE (ncrt, 8000)
 8000 FORMAT (' FATAL ERROR in subroutine TOZONE' //,                &
        ' xz(iz) or xz(izp1) is not contained in interval x(1),x(n)' /)
          lerrno = 203  + iomaxerr
          CALL terminate(lerrno,nlog)
      END SUBROUTINE tozone
      





   END MODULE prep_zones
