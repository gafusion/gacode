    MODULE prep_tport_outpt

      USE nrtype,                                ONLY : DP,I4B
      
      USE error_handler

      USE io_gcnmp,                              ONLY : ncrt,nlog

      USE nf_param,                              ONLY : ke,kcm
 
      USE neutral_beams,                         ONLY : no_physical_injectors,  &
                                                        atw_beam
 
      USE zonal_data,                            ONLY : mf,mfm1,psif,rho_zone

      USE prep_zones,                            ONLY : tozone

      USE grid_class,                            ONLY : nj,r,rho_grid

      USE common_constants,                      ONLY : izero,zeroc,Pi,M2cm,Proton_mass

      USE Plasma_properties ,                    ONLY : mhd_dat,neut_beam,dischg

      USE tension_spline,                        ONLY : spline_intrp

      USE tport_mhd_grid_data,                   ONLY : psir,volume

      USE vector_class,                          ONLY : get_values

      USE P_nfreya_interface,                    ONLY : d_fast_ion

      USE fast_ion_diffusion,                    ONLY : fi_steady_state, get_d_fast_ion, &
                                                        fi_refactor_c

      USE nub,                                   ONLY : hdepsmth,fidiff_on

      CONTAINS


      SUBROUTINE set_profiles
!----------------------------------------------------------------------
! --  this subroutine does the smoothing for the FREYA output
! --  and moves output onto the transport grid of size nj_out
!----------------------------------------------------------------------
!
 


      IMPLICIT NONE


      REAL(DP) hdep_lcl(nj),hibr_lcl(nj)
!     REAL(DP),ALLOCATABLE,DIMENSION(:) :: hcap_zone
      REAL(DP),ALLOCATABLE,DIMENSION(:) :: spbrt
      REAL(DP),ALLOCATABLE,DIMENSION(:) :: cparam
      REAL(DP) sdl,sdr,dsp,d2sp,csgn,vhdep,vhibr,volfac,slope,drrj,constv ,&
               xbir,xdep,xloss,vbeam,bpflav,stbeamt,stoloss,volume_m3,     &
               hdepfctrn,hdepfctrd,hibrfctrn,hibrfctrd,dummy


      INTEGER(I4B) i,ie,j,jb,jj,k,nj_out,ib,nsmooth,npass,nhdep,           &
                   ij,nj_outsm

      INTEGER(I4B)  ideg(4), nuse(4), nupdate(4)

      LOGICAL get_scoef,inv_weight,did_smooth


       nj_out = nj          ! hdep_lcl and hibr_lcl dims also
       neut_beam%nj_beam = nj_out
       IF(.NOT. fidiff_on)d_fast_ion%fidif_on =0 ! fidiff_on, set in nfreya namelist,
                                                 ! overides setting of d_fast_ion%fidif_on
       inv_weight = .TRUE.
       inv_weight = .FALSE.
 print *,'hdepsmth ,fidiff_on,d_fast_ion%fidif_on =',hdepsmth ,fidiff_on,d_fast_ion%fidif_on
      !----------------------------------------------------------------
      ! --- zero the FREYA calculated total beam torque:
      !-----------------------------------------------------------------
      ALLOCATE(spbrt(nj_out))
      spbrt(:) = zeroc
      !spbolt(:) = zeroc
      
      !----------------------------------------------------------------
      ! --- transport rho_grid (in meters) is obtained from statefile:
      !-----------------------------------------------------------------    
            IF(.NOT. ALLOCATED(r))ALLOCATE(r(nj_out))
            
            r(1:nj_out)        =  get_values(rho_grid)*M2cm   ! r in cm



!
! -------------------------------------------------------------------------------
! given r on psir grid get rho_zone on zone centers in cm . Note that  psif grid
! defines zone boundaries. tozone accounts for this
! obtain rho  values (cm)   (rho_zone) corresponding to psi values used in FREYA
! -------------------------------------------------------------------------------
!

       IF(.NOT. ALLOCATED(rho_zone))ALLOCATE(rho_zone(mf))
       CALL tozone (psir, r, nj_out,psif,rho_zone,mf)
 !PRINT *,'rho zone =',rho_zone(1:mf) ! 88889999
 !PRINT *,'r  =',r(1:nj) ! 88889999

!
! ----------------------------------------------------------------------
! transform sign of angmpz to transport coordinate system
! ----------------------------------------------------------------------
      csgn = 1._DP
      IF(mhd_dat%tot_cur .LT. zeroc)csgn = -1._DP
      neut_beam%angmpz(:,:,:)  = csgn*neut_beam%angmpz(:,:,:)       ! angmpz in KG M^2/SEC from load_plasma_prop
 


! -----------------------------------------------------------------------------
! convert hibrz, hdepz, zetaz, ftrapfi and angmpz from FREYA grid to transport grid
! -----------------------------------------------------------------------------


       ALLOCATE(neut_beam%hdep(nj_out,ke,no_physical_injectors), &
                                neut_beam%hibr(nj_out,ke,no_physical_injectors))
       ALLOCATE(neut_beam%zeta(nj_out,ke,no_physical_injectors))
       ALLOCATE(neut_beam%ftrapfit(nj_out,ke,no_physical_injectors))
       ALLOCATE(neut_beam%angmpf(nj_out,ke,no_physical_injectors)) 
       ALLOCATE(neut_beam%hicm(nj_out,ke,no_physical_injectors,kcm))
       neut_beam%nbeams = no_physical_injectors


         ALLOCATE(d_fast_ion%fi_d(nj_out,ke,no_physical_injectors),                  &
                  d_fast_ion%fi_den(nj_out,ke,no_physical_injectors))
                  d_fast_ion%fi_den(:,:,:) = zeroc
                  d_fast_ion%fi_d(:,:,:)   = zeroc

           IF(ALLOCATED(cparam))DEALLOCATE(cparam)
           ALLOCATE(cparam(mf))



  beamline_loop : DO  jb = 1,no_physical_injectors    ! LOOP OVER BEAMLINES

    energy_loop : DO  ie = 1,3   
     
           get_scoef = .TRUE.                        ! determine cparam on call to spline_intrp
            hibr_lcl(1:mfm1) = neut_beam%hibrz(1:mfm1,ie,jb)
            hdep_lcl(1:mfm1) = neut_beam%hdepz(1:mfm1,ie,jb)

  go to 1111
           ib =2 ! supply second derivatives
           sdl  =zeroc ; sdr = zeroc
           ! set cparam 
           CALL spline_intrp(mfm1,rho_zone,hibr_lcl,cparam,&
                             get_scoef,ib,sdl,sdr, r(1),dummy,dsp,d2sp) 
           get_scoef = .FALSE.                ! cparam is now set
           DO i = 1,nj_out
              CALL spline_intrp(mfm1,rho_zone,hibr_lcl,cparam,&
                   get_scoef,ib,sdl,sdr,r(i),neut_beam%hibr(i,ie,jb),dsp,d2sp)
           ENDDO
 
           ib =2 ! supply second derivatives
           sdl  =zeroc ; sdr = zeroc
           ! set cparam 
           get_scoef = .TRUE.
           CALL spline_intrp(mfm1,rho_zone,hdep_lcl,cparam,&
                             get_scoef,ib,sdl,sdr, r(1),dummy,dsp,d2sp) 
           get_scoef = .FALSE.                ! cparam is now set
           DO i = 1,nj_out
              CALL spline_intrp(mfm1,rho_zone,hdep_lcl,cparam,&
                   get_scoef,ib,sdl,sdr,r(i),neut_beam%hdep(i,ie,jb),dsp,d2sp)
           ENDDO
 1111       CONTINUE
        !call newgrid (rho_zone,hibr_lcl,mfm1,r,neut_beam%hibr(1,ie,jb),nj) ! convert from zone to rho grid
         call newgrid (rho_zone,hibr_lcl,mfm1,r,hdep_lcl,nj) ! convert from zone to rho grid
            neut_beam%hibr(1:nj,ie,jb) =  hdep_lcl(1:nj)
            hdep_lcl(1:mfm1) = neut_beam%hdepz(1:mfm1,ie,jb)
            !call newgrid (rho_zone,hdep_lcl,mfm1,r, neut_beam%hdep(1,ie,jb),nj) ! convert from zone to rho grid
            call newgrid (rho_zone,hdep_lcl,mfm1,r,hibr_lcl,nj) ! convert from zone to rho gri

            neut_beam%hdep(1:nj,ie,jb) = hibr_lcl(1:nj)


! ----------------------------------------------------------------------
! --- do some smoothing of deposition profiles hibr and hdep HSJ
! --- NOTE: smoothing can lead to non physical profiles ! as of 9/22/87 the
! --- smoothing can be turned on or off by input switch hdepsmth.
! --- see more hdepsmoothing further below
! ----------------------------------------------------------------------

            did_smooth = .FALSE.
 otrsmooth: IF(hdepsmth .GE. 0 )THEN



 smoothif:  IF (hdepsmth .GT. 10 .AND. d_fast_ion%fidif_on ==  0) THEN

            npass = hdepsmth - 10
            nhdep = 6  !note that this introduces a grid dependence
            DO ij=1,npass
                 neut_beam%hibr(1,ie,jb)=( 2._DP*neut_beam%hibr(1,ie,jb)+ &
                                     neut_beam%hibr(2,ie,jb)+ &
                                     neut_beam%hibr(3,ie,jb)+ &
                                     neut_beam%hibr(4,ie,jb)+ &
                                     neut_beam%hibr(5,ie,jb))/nhdep
                 neut_beam%hibr(2,ie,jb)=( neut_beam%hibr(1,ie,jb)+ &
                                     2._DP*neut_beam%hibr(2,ie,jb)+ &
                                     neut_beam%hibr(3,ie,jb)+ &
                                     neut_beam%hibr(4,ie,jb)+ &
                                     neut_beam%hibr(5,ie,jb))/nhdep
                neut_beam%hdep(1,ie,jb)=( 2._DP*neut_beam%hdep(1,ie,jb)+ &
                                    neut_beam%hdep(2,ie,jb)+ &
                                    neut_beam%hdep(3,ie,jb)+ &
                                    neut_beam%hdep(4,ie,jb)+ &
                                    neut_beam%hdep(5,ie,jb))/nhdep
                neut_beam%hdep(2,ie,jb)=(neut_beam%hdep(1,ie,jb)+ &
                                     2._DP*neut_beam%hdep(2,ie,jb)+ &
                                    neut_beam%hdep(3,ie,jb)+ &
                                    neut_beam%hdep(4,ie,jb)+ &
                                    neut_beam%hdep(5,ie,jb))/nhdep
                 DO j=3,nj_out-2
                    neut_beam%hibr(j,ie,jb)=(neut_beam%hibr(j-2,ie,jb)+ &
                                     neut_beam%hibr(j-1,ie,jb)+ &
                                     2._DP*neut_beam%hibr(j  ,ie,jb)+ &
                                     neut_beam%hibr(j+1,ie,jb)+ &
                                     neut_beam%hibr(j+2,ie,jb))/nhdep
                   neut_beam%hdep(j,ie,jb)=(neut_beam%hdep(j-2,ie,jb)+ &
                                    neut_beam%hdep(j-1,ie,jb)+ &
                                     2._DP*neut_beam%hdep(j  ,ie,jb)+ &
                                    neut_beam%hdep(j+1,ie,jb)+ &
                                    neut_beam%hdep(j+2,ie,jb))/nhdep
                 END DO
                 neut_beam%hibr(nj_out-1,ie,jb)= (neut_beam%hibr(nj_out-4,ie,jb)+ &
                                     neut_beam%hibr(nj_out-3,ie,jb)+ &
                                     neut_beam%hibr(nj_out-2,ie,jb)+ &
                                     2._DP*neut_beam%hibr(nj_out-1,ie,jb)+ &
                                     neut_beam%hibr(nj_out,ie,jb))/nhdep
                 neut_beam%hibr(nj_out,ie,jb)=(neut_beam%hibr(nj_out-4,ie,jb)+ &
                                     neut_beam%hibr(nj_out-3,ie,jb)+ &
                                     neut_beam%hibr(nj_out-2,ie,jb)+ &
                                     neut_beam%hibr(nj_out-1,ie,jb)+ &
                                     2._DP*neut_beam%hibr(nj_out,ie,jb))/nhdep
                neut_beam%hdep(nj_out-1,ie,jb)=(neut_beam%hdep(nj_out-4,ie,jb)+ &
                                    neut_beam%hdep(nj_out-3,ie,jb)+ &
                                    neut_beam%hdep(nj_out-2,ie,jb)+ &
                                     2._DP*neut_beam%hdep(nj_out-1,ie,jb)+ &
                                    neut_beam%hdep(nj_out,ie,jb))/nhdep
                neut_beam%hdep(nj_out,ie,jb)=(neut_beam%hdep(nj_out-4,ie,jb)+ &
                                    neut_beam%hdep(nj_out-3,ie,jb)+ &
                                    neut_beam%hdep(nj_out-2,ie,jb)+ &
                                    neut_beam%hdep(nj_out-1,ie,jb)+ &
                                     2._DP*neut_beam%hdep(nj_out,ie,jb))/nhdep
              END DO
             did_smooth = .TRUE.


      ELSE IF(hdepsmth .GT. 10 .AND. d_fast_ion%fidif_on == 1) THEN   smoothif ! fi_d weighted smoothing
            !get the fast ion diffusion coefficient based on nubeam parameterization model
            CALL get_d_fast_ion(d_fast_ion,ie,jb) 
            IF(inv_weight)d_fast_ion%fi_d(:,ie,jb) = 1._DP/d_fast_ion%fi_d(:,ie,jb)

            npass = hdepsmth - 10
            nj_outsm = nj_out-2


           DO ij=1,npass

            hibr_lcl(1:nj_out) = neut_beam%hibr(1:nj_out,ie,jb)
            hdep_lcl(1:nj_out) = neut_beam%hdep(1:nj_out,ie,jb)
                 ! grid point 1
                    hibrfctrn  =                                                   &
                       d_fast_ion%fi_d(2,ie,jb)*hibr_lcl(2)                      + &
                                        d_fast_ion%fi_d(3,ie,jb)*hibr_lcl(3)
                    hibrfctrd = d_fast_ion%fi_d(3,ie,jb)+ d_fast_ion%fi_d(2,ie,jb)
                 neut_beam%hibr(1,ie,jb)=(hibr_lcl(1)                              &
                                                 + hibrfctrn/hibrfctrd)/2._DP

                    hdepfctrn  =                                                   &
                       d_fast_ion%fi_d(2,ie,jb)*hdep_lcl(2)                      + &
                                    d_fast_ion%fi_d(3,ie,jb)*hdep_lcl(3)
                    hdepfctrd = hibrfctrd                                                  
                 neut_beam%hdep(1,ie,jb)=(hdep_lcl(1) + hdepfctrn/hdepfctrd)/2._DP


                 ! grid point 2
                    hibrfctrn  =  hibrfctrn + d_fast_ion%fi_d(4,ie,jb)*hibr_lcl(4)
                    hibrfctrd =  hibrfctrd  + d_fast_ion%fi_d(4,ie,jb)
                 neut_beam%hibr(2,ie,jb)= (hibr_lcl(1) + hibrfctrn/hibrfctrd)/2._DP
                    hdepfctrn  =  hdepfctrn + d_fast_ion%fi_d(4,ie,jb)*hdep_lcl(4)
                    hdepfctrd =   hibrfctrd
                 neut_beam%hdep(2,ie,jb)=(hdep_lcl(1) + hdepfctrn/hdepfctrd)/2._DP



                 ! grid point nj_out ( 3 pt avg )
                    hibrfctrn  =                                                     &
                       d_fast_ion%fi_d(nj_out-2,ie,jb)*hibr_lcl(nj_out-2)          + &
                       d_fast_ion%fi_d(nj_out-1,ie,jb)*hibr_lcl(nj_out-1)          + &
                       d_fast_ion%fi_d(nj_out,ie,jb)*hibr_lcl(nj_out)
                    hibrfctrd =                                                      &
                       d_fast_ion%fi_d(nj_out-2,ie,jb)                             + & 
                                  d_fast_ion%fi_d(nj_out-1,ie,jb)                  + &
                                      d_fast_ion%fi_d(nj_out,ie,jb)
                 neut_beam%hibr(nj_out,ie,jb) = hibrfctrn/hibrfctrd
                    hdepfctrn  =                                                     &
                       d_fast_ion%fi_d(nj_out-2,ie,jb)*hdep_lcl(nj_out-2)          + &
                       d_fast_ion%fi_d(nj_out-1,ie,jb)*hdep_lcl(nj_out-1)          + &
                             d_fast_ion%fi_d(nj_out,ie,jb)*hdep_lcl(nj_out)
                    hdepfctrd =   hibrfctrd 
                    neut_beam%hdep(nj_out,ie,jb) = hdepfctrn/hdepfctrd


                 ! grid point nj_out-1 ( 4 pt avg )
                    hibrfctrn  =                                                     &
                       d_fast_ion%fi_d(nj_out-3,ie,jb)*hibr_lcl(nj_out-3)         +  &
                       d_fast_ion%fi_d(nj_out-2,ie,jb)*hibr_lcl(nj_out-2)         +  &
                       d_fast_ion%fi_d(nj_out-1,ie,jb)*hibr_lcl(nj_out-1)         +  &
                       d_fast_ion%fi_d(nj_out,ie,jb)*hibr_lcl(nj_out)
                    hibrfctrd =                                                      &
                                d_fast_ion%fi_d(nj_out-3,ie,jb)                   +  &
                                d_fast_ion%fi_d(nj_out-2,ie,jb)                   +  &
                                d_fast_ion%fi_d(nj_out-1,ie,jb)                   +  &
                                d_fast_ion%fi_d(nj_out,ie,jb)  
                 neut_beam%hibr(nj_out-1,ie,jb) = hibrfctrn/hibrfctrd

                    hdepfctrn  =                                                     &
                       d_fast_ion%fi_d(nj_out-3,ie,jb)*hdep_lcl(nj_out-3)         +  &
                       d_fast_ion%fi_d(nj_out-2,ie,jb)*hdep_lcl(nj_out-2)         +  &
                       d_fast_ion%fi_d(nj_out-1,ie,jb)*hdep_lcl(nj_out-1)         +  &
                       d_fast_ion%fi_d(nj_out,ie,jb)*hdep_lcl(nj_out)
                    hdepfctrd =  hibrfctrd                                                 
                 neut_beam%hdep(nj_out-1,ie,jb) = hdepfctrn/hdepfctrd


               DO j=3,nj_outsm   ! 5 point average  
                    hibrfctrn  =                                                      &
                       d_fast_ion%fi_d(j-2,ie,jb)*hibr_lcl(j-2)                    +  &
                       d_fast_ion%fi_d(j-1,ie,jb)*hibr_lcl(j-1)                    +  &
                       d_fast_ion%fi_d(j,ie,jb)*hibr_lcl(j)                        +  &
                       d_fast_ion%fi_d(j+1,ie,jb)*hibr_lcl(j+1)                    +  &
                       d_fast_ion%fi_d(j+2,ie,jb)*hibr_lcl(j+2)   
                    hibrfctrd =                                                       &
                       d_fast_ion%fi_d(j-2,ie,jb)+ d_fast_ion%fi_d(j-1,ie,jb)       + &
                       d_fast_ion%fi_d(j,ie,jb)  + d_fast_ion%fi_d(j+1,ie,jb)       + &
                       d_fast_ion%fi_d(j+2,ie,jb)
                 neut_beam%hibr(j,ie,jb)= hibrfctrn/hibrfctrd

                    hdepfctrn  =                                                      &
                       d_fast_ion%fi_d(j-2,ie,jb)*hdep_lcl(j-2)                    +  &
                       d_fast_ion%fi_d(j-1,ie,jb)*hdep_lcl(j-1)                    +  &
                       d_fast_ion%fi_d(j,ie,jb)*hdep_lcl(j)                        +  &
                       d_fast_ion%fi_d(j+1,ie,jb)*hdep_lcl(j+1)                    +  &
                       d_fast_ion%fi_d(j+2,ie,jb)*hdep_lcl(j+2)
                     hdepfctrd = hibrfctrd 
                 neut_beam%hdep(j,ie,jb)=  hdepfctrn/hdepfctrd

               END DO

            END DO
            did_smooth = .TRUE.
      ELSEIF(d_fast_ion%fidif_on == 1 .AND.  hdepsmth == 0) THEN smoothif

            ! ----------------------------------------------------------------------
            ! Apply full fast ion diffusion  model 
            ! here we solve the steady state diffusion equation for the fast ion density
            ! using the (non smoothed) source rate neut_beam%sb. Should add loss
            ! term to this equation as well. Note:  not complete because hibr,hdep,etc
            ! are not changed to reflect the calculated fast ion density,d_fast_ion%fi_den
            ! Hence information is never passed to Onetwo 888888999999
            ! ----------------------------------------------------------------------

            d_fast_ion%fi_na = 1.e18_DP
            fi_refactor_c = .TRUE.
            CALL get_d_fast_ion(d_fast_ion,ie,jb) ! neut_beam%sb is not yet defined!!!  
            CALL fi_steady_state(neut_beam%sb(:,ie,jb),d_fast_ion%fi_d(:,ie,jb),    &
                              d_fast_ion%fi_na,d_fast_ion%fi_den(:,ie,jb))
            did_smooth = .TRUE.

      ELSE  smoothif
         ! remaining option is 0 .LT. hdepsmth .LT. 10 with d_fast_ion%fidif_on =0,1
         ! currently no option for this case (if option is added P_NF_nubplt eeds to
         ! be made aware of it as well
            did_smooth = .FALSE.
            hdepsmth = -1.

      END IF  smoothif
          IF(did_smooth)THEN
              !----------------------------------------------------------------------------
              ! renormalize the smoothed profiles to give the plasma volume when integrated
              !----------------------------------------------------------------------------
              vhdep = zeroc
              vhibr = zeroc
              volfac = 4._DP*Pi*Pi*dischg%rmajor*M2cm ! cm 
              CALL trapv (r,neut_beam%hdep(1,ie,jb), mhd_dat%hcap%data, nj_out, vhdep)
              CALL trapv (r, neut_beam%hibr(1,ie,jb), mhd_dat%hcap%data , nj_out, vhibr)
              vhdep=vhdep*volfac
              vhibr=vhibr*volfac

              ! it is possible for some beam/energy components to be zero
              ! because freya assignes zero particles to such a channel.
              ! (this could happend for example of the beam powers were greatly different)
              ! In this case renormalization should be skipped because the hdep,hibr are
              ! identically zero:
              IF(vhdep*vhibr .NE. 0.0)THEN
                 DO j=1,nj_out
                    neut_beam%hdep(j,ie,jb)=neut_beam%hdep(j,ie,jb)*(volume/vhdep)
                    neut_beam%hibr(j,ie,jb)=neut_beam%hibr(j,ie,jb)*(volume/vhibr)
                 END DO
              ENDIF
         ENDIF

     ENDIF  otrsmooth
      !------------------------------------------------------
      ! -- end smoothing
      !------------------------------------------------------



           ib =2 ! supply second derivatives
           sdl  =zeroc ; sdr = zeroc
           ! set cparam 
           get_scoef = .TRUE.
           CALL spline_intrp(mfm1,rho_zone,neut_beam%zetaz(:,ie,jb),cparam,&
                             get_scoef,ib,sdl,sdr, r(1),neut_beam%zeta(1,1,1),dsp,d2sp) 
           get_scoef = .FALSE.                ! cparam is now set
           DO i = 1,nj_out
              CALL spline_intrp(mfm1,rho_zone,neut_beam%zetaz(:,ie,jb),cparam,&
                   get_scoef,ib,sdl,sdr,r(i),neut_beam%zeta(i,ie,jb),dsp,d2sp)
           ENDDO


           ib =2 ! supply second derivatives
           sdl  =zeroc ; sdr = zeroc
           ! set cparam 
           get_scoef = .TRUE.
           CALL spline_intrp(mfm1,rho_zone,neut_beam%ftrapfi(:,ie,jb),cparam,&
                             get_scoef,ib,sdl,sdr, r(1),neut_beam%ftrapfi(1,1,1),dsp,d2sp) 
           get_scoef = .FALSE.                ! cparam is now set
           DO i = 1,nj_out
              CALL spline_intrp(mfm1,rho_zone,neut_beam%ftrapfi(:,ie,jb),cparam,&
                   get_scoef,ib,sdl,sdr,r(i),neut_beam%ftrapfit(i,ie,jb),dsp,d2sp)
           ENDDO



           ib =2 ! supply second derivatives
           sdl  =zeroc ; sdr = zeroc
           ! set cparam 
           get_scoef = .TRUE.
           CALL spline_intrp(mfm1,rho_zone,neut_beam%angmpz(:,ie,jb),cparam,&
                             get_scoef,ib,sdl,sdr, r(1),neut_beam%angmpz(1,1,1),dsp,d2sp) 
           get_scoef = .FALSE.                ! cparam is now set
           DO i = 1,nj_out
              CALL spline_intrp(mfm1,rho_zone,neut_beam%angmpz(:,ie,jb),cparam,&
                   get_scoef,ib,sdl,sdr,r(i),neut_beam%angmpf(i,ie,jb),dsp,d2sp)
           ENDDO

           DO k=1,3
              ib =2 ! supply second derivatives
              sdl  =zeroc ; sdr = zeroc
              ! set cparam 
              get_scoef = .TRUE.
              CALL spline_intrp(mfm1,rho_zone,neut_beam%hicmz(:,ie,jb,k),cparam,&
                             get_scoef,ib,sdl,sdr, r(1),neut_beam%hicmz(1,1,1,k),dsp,d2sp) 
              get_scoef = .FALSE.                ! cparam is now set
              DO i = 1,nj_out
                 CALL spline_intrp(mfm1,rho_zone,neut_beam%hicmz(:,ie,jb,k),cparam,&
                      get_scoef,ib,sdl,sdr,r(i),neut_beam%hicm(i,ie,jb,k),dsp,d2sp)
              ENDDO
           ENDDO

!
! --- negative zeta is possible. Undo the zeros set in newgrid if necessary:
!
      IF (neut_beam%zetaz(mfm1,ie,jb) .LT. zeroc) THEN
        DO j=nj_out,1,-1
          jj = j
          IF (r(jj) .LT. rho_zone(mfm1))  go to 510
        END DO
        WRITE  (ncrt, 600)
  600   FORMAT (' subroutine zet_profiles detected error in zetaz intrp.' /)
        lerrno = 244 + iomaxerr
        CALL terminate(lerrno,nlog)

  510   jj = jj + 1
        slope = (neut_beam%zetaz(mfm1,ie,jb)-neut_beam%zetaz(mfm1-1,ie,jb)) &
              / (rho_zone(mfm1)-rho_zone(mfm1-1))
        DO j=jj,nj_out
          drrj = r(jj)-rho_zone(mfm1)
          neut_beam%zeta(j,ie,jb) = neut_beam%zetaz(mfm1,ie,jb)+drrj*slope
        END DO
      END IF
!

 
    ENDDO  energy_loop
 
  ENDDO    beamline_loop


!
! ----------------------------------------------------------------------
! renormalize hot ion birth rate and deposition rate
! ----------------------------------------------------------------------
!
      constv = 4._DP*Pi*Pi*dischg%rmajor*M2cm                    ! cm 
      DO  jb=1,no_physical_injectors
         DO  ie=1,3
            CALL trapv (r,neut_beam%hibr(1,ie,jb),mhd_dat%hcap%data,nj_out,xbir)
            CALL trapv (r,neut_beam%hdep(1,ie,jb),mhd_dat%hcap%data,nj_out,xdep)
            IF (xbir .NE. zeroc) xbir = volume/(xbir*constv)
            IF (xdep .NE. zeroc) xdep = volume/(xdep*constv)
            DO  j=1,nj_out
               neut_beam%hibr(j,ie,jb) =neut_beam%hibr(j,ie,jb)*xbir
               neut_beam%hdep(j,ie,jb) =neut_beam%hdep(j,ie,jb)*xdep
            ENDDO
          ENDDO
      ENDDO
 


!
! ----------------------------------------------------------------------
! calculate particle, energy, and parallel momentum sources in plasma
!    due to neutral beam(s)
!      sb  is in particles/m3-s
!      qb  is in w/m3
!      spb is in kg/(m2-s2)  parallel momentum rate density
!      spbr (toroidal angular momentum source) is in (nt m/m3) = kg/(m-sec**2)
!      spbolr (toroidal torque density from orbit loss current, kg/(m-sec)   ! rs
!      pb0 flux average average initial parallel momentum per ion is in kg m/sec


! ----------------------------------------------------------------------
!
      ALLOCATE(neut_beam%sb(nj_out,ke,no_physical_injectors),neut_beam%qb(nj_out,ke,no_physical_injectors))
      ALLOCATE(neut_beam%sbpure(nj_out,ke,no_physical_injectors),neut_beam%spb(nj_out,ke,no_physical_injectors))
      ALLOCATE(neut_beam%spbr(nj_out,ke,no_physical_injectors),neut_beam%pb0(nj_out,ke,no_physical_injectors))
      volume_m3 = volume*1.e-6
      DO  jb=1,no_physical_injectors
         DO ie=1,3
            xloss = neut_beam%fap(ie,jb) + neut_beam%fwall(ie,jb) + neut_beam%forb(ie,jb)
            vbeam = 1.384e4_DP * SQRT (1.0e3_DP*neut_beam%ebeam(ie,jb)/atw_beam) ! m/sec
            DO  j=1,nj_out
               neut_beam%qb(j,ie,jb)   = (1.0-xloss)*neut_beam%pbeam(ie,jb)*neut_beam%hdep(j,ie,jb)/volume_m3 ! w/m3
               neut_beam%sb(j,ie,jb)   = 0.625e16_DP*neut_beam%qb(j,ie,jb)/neut_beam%ebeam(ie,jb)     ! #/(m3 sec)
               !
               ! --- sbpure(j,ie,jb) is set to the beam deposition sb here. sbpure is
               ! --- saved so that the modification due to time and charge exchange
               ! --- effects which are folded into sb in slow2 are not introduced.
               !
               neut_beam%sbpure(j,ie,jb)= neut_beam%sb(j,ie,jb)            !   #/(m3 sec)              
               neut_beam%spb(j,ie,jb)   = atw_beam*Proton_Mass*vbeam* &
                            neut_beam%zeta(j,ie,jb)*neut_beam%sb(j,ie,jb)    ! Kg/(m2-s2)
               neut_beam%spbr(j,ie,jb)  = neut_beam%angmpf(j,ie,jb)*neut_beam%sb(j,ie,jb)                   ! Kg/(m sec^2)
               spbrt(j)       = spbrt(j) + neut_beam%spbr(j,ie,jb)
               bpflav = zeroc                                              
               !***  if (j .ne. 1)                                              ! rsHSJ
               !*** .  bpflav = rbp(j)/(fcap(j)*mhd_dat%hcap%data(j)*gcap(j)*r(j)) / 1000.0  ! rsHSJ
               !bpflav = bprmaj(j)                           ! bprmaj is in kgauss
               !spbolr(j,ie,jb) = -oloss(j,ie,jb)*bpflav*mhd_dat%rcap%data(j)*drho(j)   ! rsHSJ
               !spbolr(j,ie,jb) = spbolr(j,ie,jb)*100.0                    ! rsHSJ
               !
               ! --- dv = (4*pi**2) * R0 * mhd_dat%hcap%data * r * dr = constv * hcap * r * dr ! HSJ
               ! --- assume oloss(amps),bpflav is flux averaged bp in (kgauss)  
               ! --- rcap(j) is (time-evolved) flux-averaged major radius (cm), 
               ! --- drho has units  (cm),rbp is in gauss-cm                    
               ! --- spbolr (gm-cm/sec**2)                                      
               ! --- accumulate total orbit loss current torque                 
               !       for routine source                                      
               !
               !spbolt(j)    = spbolt(j)+spbolr(j,ie,jb)                   ! rs
               neut_beam%pb0(j,ie,jb) = neut_beam%angmpf(j,ie,jb)*mhd_dat%rcap%data(j)/mhd_dat%r2capi%data(j)                 !kg m/sec
            ENDDO
         ENDDO
      ENDDO

!-------------------------------------------------------------------------
! --- calculate the beam and orbit loss current torque in nt-m:  ! rs
! --- this is the prompt beam torque as determined directly from FREYA
!-------------------------------------------------------------------------
      CALL trapv (r,spbrt,mhd_dat%hcap%data,nj_out,stbeamt)     ! r in cm, spbrt in  Kg/(m sec^2)  
      stbeamt = stbeamt*1.0e-06*constv                          ! constv in cm ,stbeamt in nt m

!      CALL trapv (r,spbolt,mhd_dat%hcap%data,nj_out,stoloss)                      ! rs
!      stoloss = stoloss*1.0E-06                                              ! rs
!
! ----------------------------------------------------------------------
! convert flux values to kgauss-cm**2
! ----------------------------------------------------------------------
!
!      IF (codeid .EQ. 'onedee')  go to 2320
!      DO 2310 i=1,nw
!      DO 2310 j=1,nh
! 2310   p(i,j) = 1.0e-3*p(i,j)
! 2320 CONTINUE

!
#ifdef prepout
#endif
      RETURN
      END  SUBROUTINE set_profiles



      SUBROUTINE newgrid (xold, yold, nold, xnew, ynew, nnew)
! ----------------------------------------------------------------------
! convert to new grid
! ----------------------------------------------------------------------

      IMPLICIT  NONE

      REAL(DP)   xold(*), yold(*), xnew(*), ynew(*)
      REAL(DP) tol,sex,del,s
      INTEGER(I4B) iold,inew,nnew,nold,ipm1
      DATA       tol/1.0e-20_DP/

      iold = 1
      inew = 1
      sex  = (yold(nold)-yold(nold-1))/(xold(nold)-xold(nold-1))

! ----------------------------------------------------------------------
! interpolate to new grid that is compute ynew(xnew(inew))
! ----------------------------------------------------------------------

 2000 IF (inew .GT. nnew)  RETURN
 2100 IF (iold .GT. nold)  go to 2200
      IF (xnew(inew) .GT. xold(1))  go to 2105
      ynew(inew) = yold(1)
      inew       = inew+1
      go to 2000

 2105 del  = reldif(xold(iold),xnew(inew))
      IF (del .GT. tol)  go to 2110
      ynew(inew) = yold(iold)
      inew = inew+1
      go to 2000
 2110 IF (xold(iold) .GT. xnew(inew))  go to 2120
      iold = iold+1
      go to 2100
 2120 ipm1 = iold-1
      s    = (yold(iold)-yold(ipm1))/(xold(iold)-xold(ipm1))
      ynew(inew) = yold(ipm1)+(xnew(inew)-xold(ipm1))*s
      inew = inew+1
      go to 2000
 2200 ynew(inew) = yold(nold)+sex*(xnew(inew)-xold(nold))
      IF (ynew(inew) .LT.  0.0_DP) ynew(inew) = 0.0_DP
      inew       = inew + 1
      go to 2000

      END SUBROUTINE newgrid



      REAL(DP)  FUNCTION reldif (x1, x2)

      IMPLICIT  NONE
      REAL(DP) x1,x2,del,xabs

      del    = ABS (x1-x2)
      xabs   = ABS (x1   )
      IF (xabs .GT. 1.0e-20_DP)  del = del / xabs
      reldif = del
      RETURN

      END FUNCTION reldif



    END MODULE prep_tport_outpt
