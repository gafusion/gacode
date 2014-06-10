!-------------------------------------------------------------------------------
! --  multimode model version 7.1 ,2/29/2012,v7.1.11,10/?/2012
! --  variable storage
!-----------------------------------------------------------------------HSJ-----
    MODULE mmm_in
      USE nrtype,                                           ONLY : DP,I4B

      USE modmmm7_1,                                        ONLY : MAXNOPT

      USE solcon_gcnmp,                                     ONLY : ntot,select_solver,&
                                                                   itran_max


         IMPLICIT NONE
         PRIVATE

         REAL(DP), PARAMETER :: &
                                zepslon = 2._DP**(-53) ! As defined by IEEE 754-2008 standard

         INTEGER, PARAMETER  :: &
         NPMAX  = 200, &! Maximum data points allowed in arrays
         hfIn   = 34,  &! Input file handle
         hfOut  = 35,  &! Output file handle
         hfDebug= 36    ! Diagnostic output

        !.. Constants for missing input dectection
         REAL(DP),PUBLIC, PARAMETER  :: BADREAL = -1E10_DP
         INTEGER, PUBLIC, PARAMETER  :: BADINT = -1000000


        !------- Input Variables --------------------------------------------------
        ! The content of the input variables follows their corresponding arguments
        ! of the mmm7_1 subroutine. See modmmm7_1.f90 for more details.
        !--------------------------------------------------------------------------
         INTEGER,PUBLIC :: lprint, npoints, nerr

         REAL(DP),PUBLIC, DIMENSION(3)       :: &
              mmm_cmodel                                   ! Internal model weights

         REAL(DP),PUBLIC, DIMENSION(MAXNOPT) :: &
              mmm_cW20, mmm_cDBM, mmm_cETG                 ! To be passed to cswitch

         INTEGER, PUBLIC,DIMENSION(MAXNOPT)  :: &
              mmm_lW20, mmm_lDBM, mmm_lETG                 ! To be passed to lswitch



         !---------------------------------------------------------------------------
         !.. Input profile variables of the first kind
         ! -- **_ze ==> zone_edge , from j= 1 (mag axis) to j = npoints (plasma edge)
         !---------------------------------------------------------------------------
        ! REAL(DP), DIMENSION(NPMAX) :: &
              !rmin   ,&! half-width of the mgnetic surface [m]
              !rmaj   ,&! major radius to geometric center of the magnetic surface [m]
              !elong  ,&! local elongation
              !ne     ,&! electron density [m^-3]
              !nh     ,&! thermal hydrogenic ion densitie [m^-3]
              !nz     ,&! impurity ion densitie [m^-3]
              !nf     ,&! density from fast (non-thermal) ions [m^-3]
              !zeff   ,&! Z_eff
              !te     ,&! T_e (electron temperature) [keV]
              !ti     ,&! T_i (temperature of thermal ions) [keV]
              !q      ,&! mgnetic q-value
              !btor   ,&! ( R B_tor ) / rmaj(jz)  [Tesla]
              !zimp   ,&! mean charge of impurities
              !aimp   ,&! mean atomic mass of impurities
              !ahyd   ,&! mean atomic mass of hydrogen ions
              !aimass ,&! mean atomic mass of thermal ions
              !wexbs  ,&! ExB shearing rate in [rad/s]
              !mmm_gne_ze    ,&! -R ( d n_e / d r ) / n_e
              !mmm_gni_ze    ,&! -R ( d n_i / d r ) / n_i
              !mmm_gnh_ze    ,&! -R ( d n_h / d r ) / n_h
              !mmm_gnz_ze    ,&! -R ( d Z n_Z / d r ) / ( Z n_Z )
              !mmm_gte_ze    ,&! -R ( d T_e / d r ) / T_e zone_edge
              !mmm_gti_ze    ,&! -R ( d T_i / d r ) / T_i
              !gq     ,&!  R ( d q   / d r ) / q
              !gvtor  ,&! Normalized toroidal velocity gradient (r/v_tor)*(dv_tor/dr)
              !vtor   ,&! Toroidal velocity [m/s]
              !gvpol  ,&! Normalized poloidal velocity gradient (r/v_pol)*(dv_pol/dr)
              !vpol   ,&! Poloidal velocity [m/s]
              !gvpar  ,&! Normalized parallel velocity gradient (r/v_par)*(dv_par/dr)
              !vpar,   &! Parallel velocity [m/s]
              !gelong   ! dkappa/dr with r = rmin/R
         REAL(DP),PUBLIC ::     mmm_rmaj0
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gne_ze    ! -R ( d n_e / d r ) / n_e
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gni_ze    ! -R ( d n_i / d r ) / n_i
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gnh_ze    ! -R ( d n_h / d r ) / n_h
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gnz_ze    ! -R ( d Z n_Z / d r ) / ( Z n_Z )
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gte_ze    ! -R ( d T_e / d r ) / T_e 
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gti_ze    ! -R ( d T_i / d r ) / T_i
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gq_ze     !  R ( d q   / d r ) / q
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gvtor_ze  !  R ( d vtor  / d r ) / q
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gelong_ze  !  R ( d vtor  / d r ) / q
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_rmaj_ctr_ze   ! major radius to geometric center of the magnetic surface [m]

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_rmin_ze   ! half-width of the magnetic surface [m]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_elong_ze  ! local elongation

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ne_ze     ! electron density [m^-3]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_nh_ze     ! thermal hydrogenic ion densitie [m^-3]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_nz_ze     ! impurity ion densitie [m^-3]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ni_ze     ! total thermal ion density [m^-3]

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_nf_ze     ! density from fast (non-thermal) ions [m^-3]

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_zeff_ze   ! Z_eff
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_te_ze     ! T_e (electron temperature) [keV]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ti_ze     ! T_i (temperature of thermal ions) [keV]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_q_ze      ! magnetic q-value

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_btor_ze   ! ( R B_tor ) / rmaj(jz)  [Tesla]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_bpol_ze   ! Bpol [Tesla]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_zimp_ze   ! mean charge of impurities
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_aimp_ze   ! mean atomic mass of impurities
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ahyd_ze   ! mean atomic mass of hydrogen ions
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_aimass_ze ! mean atomic mass of thermal ions  
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_wexbs_ze  ! ExB shearing rate in [rad/s]

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vpol_ze   !  Vpoloidal,m/sec
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gvpol_ze  !  Normalized poloidal velocity gradient (r/v_pol)*(dv_pol/dr)
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vpar_ze   !  Vparallel,m/sec
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gvpar_ze  !  Normalized parallel velocity gradient (r/v_par)*(dv_par/dr)

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vtor_ze   !  Vtoroidal,m/sec


         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_er_ze     !  radial eledtric field V/m

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_theta_ze  !  approx = Bp0/Bt0

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vdia_ze   !  diamagentic velocity


         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vexb_ze   !  
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_tor_ang_speed_ze  ! Er/(RBp)


         !---------------------------------------------------------------------------------
         ! zone ctr quantities 
         !---------------------------------------------------------------------------------
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gne_zct    ! -R ( d n_e / d r ) / n_e
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gni_zct    ! -R ( d n_i / d r ) / n_i
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gnh_zct    ! -R ( d n_h / d r ) / n_h
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gnz_zct    ! -R ( d Z n_Z / d r ) / ( Z n_Z )
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gte_zct    ! -R ( d T_e / d r ) / T_e 
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gti_zct    ! -R ( d T_i / d r ) / T_i
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gq_zct     !  R ( d q   / d r ) / q
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gvtor_zct  !  R ( d vtor  / d r ) / q
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gelong_zct  !  R ( d vtor  / d r ) / q
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_rmaj_ctr_zct   ! major radius to geometric center of the magnetic surface [m]

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_rmin_zct   ! half-width of the magnetic surface [m]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_elong_zct  ! local elongation

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ne_zct     ! electron density [m^-3]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_nh_zct     ! thermal hydrogenic ion densitie [m^-3]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_nz_zct     ! impurity ion densitie [m^-3]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ni_zct     ! total thermal ion density [m^-3]

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_nf_zct     ! density from fast (non-thermal) ions [m^-3]


         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_te_zct     ! T_e (electron temperature) [keV]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ti_zct     ! T_i (temperature of thermal ions) [keV]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_q_zct      ! magnetic q-value
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_zeff_zct    ! Z_eff
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_btor_zct   ! ( R B_tor ) / rmaj(jz)  [Tesla]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_bpol_zct   ! Bpol [Tesla]
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_zimp_zct   ! mean charge of impurities
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_aimp_zct   ! mean atomic mass of impurities
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_ahyd_zct   ! mean atomic mass of hydrogen ions
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_aimass_zct ! mean atomic mass of thermal ions  
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_wexbs_zct  ! ExB shearing rate in [rad/s]

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vpol_zct   !  Vpoloidal,m/sec
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gvpol_zct  !  Normalized poloidal velocity gradient (r/v_pol)*(dv_pol/dr)
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vpar_zct   !  Vparallel,m/sec
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_gvpar_zct  !  Normalized parallel velocity gradient (r/v_par)*(dv_par/dr)

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vtor_zct   !  Vtoroidal,m/sec


         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_er_zct     !  radial eledtric field V/m

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_theta_zct  !  approx = Bp0/Bt0

         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vdia_zct   !  diamagentic velocity


         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_vexb_zct   !  
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_tor_ang_speed_zct  ! Er/(RBp)
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_drhodrmin_zct ! dr/drmin
         REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:) ::  mmm_drhodrmin_ze ! dr/drmin
         !---------------------------------------------------------------------------
         !.. Input variables of the second kind
         !---------------------------------------------------------------------------
         REAL(DP) :: &
              k_rminor ,&! Minor radius of plasma [m]
              k_rmajor ,&! Major radius at plasma axis [m]
              k_elong  ,&! Elongtion
              k_denmin ,&! Minimal electron density [m^-3]
              k_temin  ,&! Minimal electron temperature [keV]
              k_btor   ,&! Toroidal mgnetic field [Tesla]
              k_amassh ,&! Mean atomic mass of hydrogenic ions
              k_amassz   ! Mean atomic mass of impurities



         !---------------------------------------------------------------------------
         !
         !.. Parameters for generting profiles
         !   See documenttion for more explantion
         !---------------------------------------------------------------------------
         REAL(DP) :: &
              denhaxis, denhedge, denhexp, &! Ion density
              denzaxis, denzedge, denzexp, &! Average impurity density
              chrzaxis, chrzedge, chrzexp, &! Average impurity charge
              denfaxis, denfedge, denfexp, &! Fast ion density
              chrfaxis, chrfedge, chrfexp, &! Average charge of super-thermal ions
              teaxis,   teedge,   teexp,   &! Electron temperature
              tiaxis,   tiedge,   tiexp,   &! Ion temperature
              qaxis,    qedge,    qexp,    &! Mgnetic q
              wexbmax,  xwexbinn, xwexbout  ! ExB flow shear rate


         !---------------------------------------------------------------------------
         !------- Output Variables ------------------------------------------
         ! The content of the output variables follows their corresponding arguments
         ! of the mmm7_1 subroutine. See modmmm7_1.f90 for more details.
         !---------------------------------------------------------------------------
             !REAL(DP), DIMENSION(NPMAX) :: &
             !     xti, xdi, xte, xdz, xvt, xvp, gammaDBM, omegaDBM,      &
             !     xtiW20, xdiW20, xteW20, xtiDBM, xdiDBM, xteDBM, xteETG
             REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE ::                         &     
                  mmm_xti_ze,  mmm_xdi_ze,  mmm_xte_ze,  mmm_xdz_ze,      &
                  mmm_xvt_ze,  mmm_xvp_ze, mmm_gammaDBM_ze,               &
                  mmm_omegaDBM_ze, mmm_xtiW20_ze,  mmm_xdiW20_ze,         &
                  mmm_xteW20_ze, mmm_xtiDBM_ze,mmm_xdiDBM_ze,             &
                  mmm_xteDBM_ze,  mmm_xteETG_ze
             REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE ::                         &     
                  mmm_xti_zct,  mmm_xdi_zct,  mmm_xte_zct,  mmm_xdz_zct,      &
                  mmm_xvt_zct,  mmm_xvp_zct, mmm_gammaDBM_zct,               &
                  mmm_omegaDBM_zct, mmm_xtiW20_zct,  mmm_xdiW20_zct,         &
                  mmm_xteW20_zct, mmm_xtiDBM_zct,mmm_xdiDBM_zct,             &
                  mmm_xteDBM_zct,  mmm_xteETG_zct

             !REAL(DP), DIMENSION(4,NPMAX) :: &
             !     gammaW20, omegaW20
             REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE ::  mmm_gammaW20_ze,         &
                                                      mmm_omegaW20_ze
             REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE ::  mmm_gammaW20_out_ze,    &
                                                      mmm_omegaW20_out_ze
             REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE ::  mmm_gammaW20_zct,         &
                                                      mmm_omegaW20_zct
             REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE ::  mmm_gammaW20_out_zct,    &
                                                      mmm_omegaW20_out_zct

             !REAL(DP), DIMENSION(6,NPMAX) :: &
             !     vconv, vflux
            REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE :: mmm_vconv_ze, mmm_vflux_ze
            REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE :: mmm_vconv_out_ze, mmm_vflux_out_ze
            REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE :: mmm_vconv_zct, mmm_vflux_zct
            REAL(DP),PUBLIC, DIMENSION(:,:),ALLOCATABLE :: mmm_vconv_out_zct, mmm_vflux_out_zct
         !---------------------------------------------------------------------------
         !------- Local Variables ------------------------------------------
         !---------------------------------------------------------------------------
             !.. Internal variables
             !REAL(DP), DIMENSION(NPMAX) :: &
             !     zchrgz, zchrgf, zxb, zdensf, zdensi, zgrdnf
             REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE ::   &
                  mmm_zchrgz_ze, mmm_zchrgf_ze, mmm_zxb_ze, mmm_zdensf_ze, &
                  mmm_zdensi_ze, mmm_zgrdnf_ze
             REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE ::   &
                  mmm_zchrgz_zct, mmm_zchrgf_zct, mmm_zxb_zct, mmm_zdensf_zct, &
                  mmm_zdensi_zct, mmm_zgrdnf_zct

             REAL(DP) :: &
                  zdx, zxfact

             REAL(DP),PUBLIC :: mmm_cswitch(MAXNOPT,4)

             INTEGER,PUBLIC  :: mmm_lswitch(MAXNOPT,4)

             INTEGER  :: jr

             INTEGER  :: input_kind

             INTEGER  :: errcount

             CHARACTER(len=255) :: strbuf
    
             LOGICAL  :: mmm_initialized



         !---------------------------------------------------------------
         ! -- for gcnmp
         !---------------------------------------------------------------
             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:)           ::   &
                  mmm_p_output,mmm_e_output,mmm_dcoef_zc,              &
                  mmm_fluxipt_zc

             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:)             ::   &
                  mmm_tot_angmtm_zc,mmm_receive_packed,                &
                  mmm_send_packed,mmm_fluxeth_zc,                      &
                  mmm_fluxith_zc,mmm_typ_val
 



             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:,:)         ::   &  ! used for nwt_pred_cor solver
                  mmm_vp,mmm_vm,mmm_d                                     ! allocated only if solver is selected
                                                                          ! (see load_mmm_flux_mat)

             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:,:)         ::   &  ! used for nwt_pred_cor solver
                  mmm_p_flux_zc,mmm_e_flux_zc,mmm_m_flux_zc                        ! allocated only if solver is selected
                                                                          ! (see collect_mmm_output)

             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:)           ::   &  ! used for nwt_pred_cor solver
                  mmm_val_pert_ze,  mmm_val_base_ze                       ! allocated only if solver is selected
             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:)           ::   &  ! used for nwt_pred_cor solver
                  mmm_val_pert_zct,  mmm_val_base_zct                       ! allocated only if solver is selected
                                                         
             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:)             ::   &  ! used for concistency of interface only
                  mmm_angrot_ze,mmm_gangrot_ze
             REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:)             ::   &  ! used for concistency of interface only
                  mmm_angrot_zct,mmm_gangrot_zct


             REAL(DP),PUBLIC ::  mmm_d_mult,mmm_v_mult,mmm_scfctr

             INTEGER, PUBLIC, PARAMETER  :: nspecies_mmm = 4 !(fast ion,hydrogenic ion,electrons,impurity ion)

             INTEGER(I4B),PUBLIC :: mmm_max_pack
             LOGICAL, PUBLIC     :: recall_mmm
             LOGICAL, PUBLIC     :: use_mmm_pert,mmm_unstable,use_mmm_zct,mmm_loaded_ze

             PUBLIC  mmm_descrip, mmm_allocate
             PUBLIC  set_mmm_vars,set_mmm_vars_zct,mmm_load_statefile_vectors
             PUBLIC  mmm_single_grid_point,mmm_single_grid_point_zct,mmm_zct2ze
     


    CONTAINS

         SUBROUTINE mmm_descrip
         !--------------------------------------------------------------------------
         ! -- write values of multimode switchdes to log file
         !--------------------------------------------------------------------------
              USE io_gcnmp,                  ONLY : ncrt,nlog

              IMPLICIT  NONE
              INTEGER(I4B) j


                 WRITE(nlog,FMT = '(2x,35("_"),"multimode  transport",35("_"))')
                 
                 WRITE(nlog,FMT='("  Weiland, Dribm, Etg weights = ",3(1pe8.1,2x))')mmm_cmodel
                 WRITE(nlog,FMT='("    Weiland adjustable switches:")')
                 WRITE(nlog,FMT='("       E x B shear multiplier: ",1pe12.4)')mmm_cswitch(1,1)
                 WRITE(nlog,FMT='("       Momentum Pinch scale factor: ",1pe12.4)')mmm_cswitch(2,1)

                 WRITE(nlog,FMT='("       Lower bound electron thermal diffusivity: ",1pe12.4)')mmm_cswitch(3,1)
                 WRITE(nlog,FMT='("       Upper bound electron thermal diffusivity: ",1pe12.4)')mmm_cswitch(4,1)
                 WRITE(nlog,FMT='("       Lower bound ion thermal diffusivity: ",1pe12.4)')mmm_cswitch(5,1)
                 WRITE(nlog,FMT='("       Upper bound ion thermal diffusivity: ",1pe12.4)')mmm_cswitch(6,1)

                 WRITE(nlog,FMT='("  DRIBM adjustable switches:")')
                 WRITE(nlog,FMT='("     E x B shear multiplier: ",1pe12.4)')mmm_cswitch(1,2)



                 WRITE(nlog,FMT='("  ETG   adjustable switches:")')
                 WRITE(nlog,FMT='(" PUBLIC    Scaling factor for electrostatic regime: ",1pe12.4)')mmm_cswitch(1,3)
                 WRITE(nlog,FMT='("     Scaling factor for electromagnetic regime: ",1pe12.4)')mmm_cswitch(2,3)
                 WRITE(nlog,FMT='("     ETG threshold selector",i3)')mmm_lswitch(1,3)
                 WRITE(nlog,FMT = '(2x,80("-"),//,//)')
              RETURN

         END SUBROUTINE mmm_descrip



         SUBROUTINE mmm_allocate(njin)
         !--------------------------------------------------------------------------
         ! -- ALLOCATE mmm storage arrays
         !--------------------------------------------------------------------------
              USE grid_class,                      ONLY : nj

              USE solcon_gcnmp,                    ONLY :                    &
                                                          n_slave_grid_tot,  &
                                                          t_slave_compt_tot, &
                                                          select_solver

              USE ions_gcnmp,                      ONLY : nion

              USE MPI_data,                        ONLY : master,myid,numprocs
     
              USE common_constants,                ONLY : izero,zeroc

              IMPLICIT NONE 
              INTEGER(I4B) njin


              mmm_initialized     = .FALSE.
              mmm_max_pack      = izero         ! counts # of items to pack and send
                                                  ! in master/slave routines


              IF(ALLOCATED(mmm_p_output))DEALLOCATE(mmm_p_output)
              ALLOCATE(mmm_p_output(njin,3))
              mmm_p_output= zeroc

              IF(ALLOCATED(mmm_e_output))DEALLOCATE(mmm_e_output)
              ALLOCATE(mmm_e_output(njin,3))
              mmm_e_output = zeroc

              IF(ALLOCATED(mmm_gte_ze))DEALLOCATE(mmm_gte_ze)
              ALLOCATE(mmm_gte_ze(njin))
              mmm_gte_ze = zeroc

              IF(ALLOCATED(mmm_gti_ze))DEALLOCATE(mmm_gti_ze)
              ALLOCATE(mmm_gti_ze(njin))
              mmm_gti_ze = zeroc

              IF(ALLOCATED(mmm_gne_ze))DEALLOCATE(mmm_gne_ze)
              ALLOCATE(mmm_gne_ze(njin))
              mmm_gne_ze = zeroc

              IF(ALLOCATED(mmm_gnh_ze))DEALLOCATE(mmm_gnh_ze)
              ALLOCATE(mmm_gnh_ze(njin))
              mmm_gnh_ze = zeroc

              IF(ALLOCATED(mmm_gnz_ze))DEALLOCATE(mmm_gnz_ze)
              ALLOCATE(mmm_gnz_ze(njin))
              mmm_gnz_ze = zeroc

              IF(ALLOCATED(mmm_gni_ze))DEALLOCATE(mmm_gni_ze)
              ALLOCATE(mmm_gni_ze(njin))
              mmm_gni_ze = zeroc

              IF(ALLOCATED(mmm_gq_ze))DEALLOCATE(mmm_gq_ze)
              ALLOCATE(mmm_gq_ze(njin))
              mmm_gq_ze = zeroc

              IF(ALLOCATED(mmm_rmaj_ctr_ze))DEALLOCATE(mmm_rmaj_ctr_ze)
              ALLOCATE(mmm_rmaj_ctr_ze(njin))
              mmm_rmaj_ctr_ze = zeroc

              IF(ALLOCATED(mmm_rmin_ze))DEALLOCATE(mmm_rmin_ze)
              ALLOCATE(mmm_rmin_ze(njin))
              mmm_rmin_ze = zeroc

              IF(ALLOCATED(mmm_elong_ze))DEALLOCATE(mmm_elong_ze)
              ALLOCATE(mmm_elong_ze(njin))
              mmm_elong_ze = zeroc

              IF(ALLOCATED(mmm_ne_ze))DEALLOCATE(mmm_ne_ze)
              ALLOCATE(mmm_ne_ze(njin))
              mmm_ne_ze = zeroc

              IF(ALLOCATED(mmm_nh_ze))DEALLOCATE(mmm_nh_ze)
              ALLOCATE(mmm_nh_ze(njin))
              mmm_nh_ze = zeroc

              IF(ALLOCATED(mmm_ni_ze))DEALLOCATE(mmm_ni_ze)
              ALLOCATE(mmm_ni_ze(njin))
              mmm_ni_ze = zeroc

              IF(ALLOCATED(mmm_nz_ze))DEALLOCATE(mmm_nz_ze)
              ALLOCATE(mmm_nz_ze(njin))
              mmm_nz_ze = zeroc

              IF(ALLOCATED(mmm_nf_ze))DEALLOCATE(mmm_nf_ze)
              ALLOCATE(mmm_nf_ze(njin))
              mmm_nf_ze = zeroc

              IF(ALLOCATED(mmm_te_ze))DEALLOCATE(mmm_te_ze)
              ALLOCATE(mmm_te_ze(njin))
              mmm_te_ze = zeroc

              IF(ALLOCATED(mmm_ti_ze))DEALLOCATE(mmm_ti_ze)
              ALLOCATE(mmm_ti_ze(njin))
              mmm_ti_ze = zeroc

              IF(ALLOCATED(mmm_zeff_ze))DEALLOCATE(mmm_zeff_ze)
              ALLOCATE(mmm_zeff_ze(njin))
              mmm_zeff_ze = zeroc

              IF(ALLOCATED(mmm_q_ze))DEALLOCATE(mmm_q_ze)
              ALLOCATE(mmm_q_ze(njin))
              mmm_q_ze = zeroc

              IF(ALLOCATED(mmm_btor_ze))DEALLOCATE(mmm_btor_ze)
              ALLOCATE(mmm_btor_ze(njin))
              mmm_btor_ze = zeroc

              IF(ALLOCATED(mmm_bpol_ze))DEALLOCATE(mmm_bpol_ze)
              ALLOCATE(mmm_bpol_ze(njin))
              mmm_bpol_ze = zeroc

              IF(ALLOCATED(mmm_zimp_ze))DEALLOCATE(mmm_zimp_ze)
              ALLOCATE(mmm_zimp_ze(njin))
              mmm_zimp_ze = zeroc

              IF(ALLOCATED(mmm_aimp_ze))DEALLOCATE(mmm_aimp_ze)
              ALLOCATE(mmm_aimp_ze(njin))
              mmm_aimp_ze = zeroc

              IF(ALLOCATED(mmm_ahyd_ze))DEALLOCATE(mmm_ahyd_ze)
              ALLOCATE(mmm_ahyd_ze(njin))
              mmm_ahyd_ze = zeroc

              IF(ALLOCATED(mmm_aimass_ze))DEALLOCATE(mmm_aimass_ze)
              ALLOCATE(mmm_aimass_ze(njin))
              mmm_aimass_ze  = zeroc


              IF(ALLOCATED(mmm_wexbs_ze))DEALLOCATE(mmm_wexbs_ze)
              ALLOCATE(mmm_wexbs_ze(njin))
              mmm_wexbs_ze = zeroc


              IF(ALLOCATED(mmm_vpol_ze))DEALLOCATE(mmm_vpol_ze)
              ALLOCATE(mmm_vpol_ze(njin)) 
              mmm_vpol_ze = zeroc

              IF(ALLOCATED(mmm_gvpol_ze))DEALLOCATE(mmm_gvpol_ze)
              ALLOCATE(mmm_gvpol_ze(njin))
              mmm_gvpol_ze = zeroc

              IF(ALLOCATED(mmm_vpar_ze))DEALLOCATE(mmm_vpar_ze)
              ALLOCATE(mmm_vpar_ze(njin))
              mmm_vpar_ze= zeroc

              IF(ALLOCATED(mmm_gvpar_ze))DEALLOCATE(mmm_gvpar_ze)
              ALLOCATE(mmm_gvpar_ze(njin))
              mmm_gvpar_ze= zeroc

              IF(ALLOCATED(mmm_vtor_ze))DEALLOCATE(mmm_vtor_ze)
              ALLOCATE(mmm_vtor_ze(njin))
              mmm_vtor_ze= zeroc

              IF(ALLOCATED(mmm_gvtor_ze))DEALLOCATE(mmm_gvtor_ze)
              ALLOCATE(mmm_gvtor_ze(njin))
              mmm_gvtor_ze= zeroc

              IF(ALLOCATED(mmm_gelong_ze))DEALLOCATE(mmm_gelong_ze)
              ALLOCATE(mmm_gelong_ze(njin))
              mmm_gelong_ze= zeroc

              IF(ALLOCATED(mmm_er_ze))DEALLOCATE(mmm_er_ze)
              ALLOCATE(mmm_er_ze(njin))
              mmm_er_ze = zeroc

              IF(ALLOCATED(mmm_theta_ze))DEALLOCATE(mmm_theta_ze)
              ALLOCATE(mmm_theta_ze(njin))
              mmm_theta_ze = zeroc


              IF(ALLOCATED(mmm_vdia_ze))DEALLOCATE(mmm_vdia_ze)
              ALLOCATE(mmm_vdia_ze(njin))
              mmm_vdia_ze= zeroc


              IF(ALLOCATED(mmm_vexb_ze))DEALLOCATE(mmm_vexb_ze)
              ALLOCATE(mmm_vexb_ze(njin))
              mmm_vexb_ze = zeroc

              IF(ALLOCATED(mmm_tor_ang_speed_ze))DEALLOCATE(mmm_tor_ang_speed_ze)
              ALLOCATE(mmm_tor_ang_speed_ze(njin))
              mmm_tor_ang_speed_ze = zeroc

              IF(ALLOCATED(mmm_angrot_ze))DEALLOCATE(mmm_angrot_ze)
              ALLOCATE(mmm_angrot_ze(njin))
              mmm_angrot_ze = zeroc


              IF(ALLOCATED(mmm_gangrot_ze))DEALLOCATE(mmm_gangrot_ze)
              ALLOCATE(mmm_gangrot_ze(njin))
              mmm_gangrot_ze = zeroc

              IF(ALLOCATED(mmm_drhodrmin_ze))DEALLOCATE(mmm_drhodrmin_ze)
              ALLOCATE(mmm_drhodrmin_ze(njin))
              mmm_drhodrmin_ze = zeroc

              !IF(select_solver == 'nwt_pred_cor')THEN
              ! must be defined for message passing even for maccormack solver
                 IF(ALLOCATED(mmm_val_pert_ze))DEALLOCATE(mmm_val_pert_ze)
                 ALLOCATE(mmm_val_pert_ze(-itran_max:itran_max,0:njin-1))
                 mmm_val_pert_ze = zeroc
                 IF(ALLOCATED(mmm_val_base_ze))DEALLOCATE(mmm_val_base_ze)
                 ALLOCATE(mmm_val_base_ze(-itran_max:itran_max,0:njin-1))
                 mmm_val_base_ze = zeroc

                 IF(ALLOCATED(mmm_typ_val))DEALLOCATE(mmm_typ_val)
                 ALLOCATE(mmm_typ_val(-itran_max : itran_max))
              !ENDIF


              IF(ALLOCATED(mmm_gte_zct))DEALLOCATE(mmm_gte_zct)
              ALLOCATE(mmm_gte_zct(njin-1))
              mmm_gte_zct = zeroc

              IF(ALLOCATED(mmm_gti_zct))DEALLOCATE(mmm_gti_zct)
              ALLOCATE(mmm_gti_zct(njin-1))
              mmm_gti_zct = zeroc

              IF(ALLOCATED(mmm_gne_zct))DEALLOCATE(mmm_gne_zct)
              ALLOCATE(mmm_gne_zct(njin-1))
              mmm_gne_zct = zeroc

              IF(ALLOCATED(mmm_gnh_zct))DEALLOCATE(mmm_gnh_zct)
              ALLOCATE(mmm_gnh_zct(njin-1))
              mmm_gnh_zct = zeroc

              IF(ALLOCATED(mmm_gnz_zct))DEALLOCATE(mmm_gnz_zct)
              ALLOCATE(mmm_gnz_zct(njin-1))
              mmm_gnz_zct = zeroc

              IF(ALLOCATED(mmm_gni_zct))DEALLOCATE(mmm_gni_zct)
              ALLOCATE(mmm_gni_zct(njin-1))
              mmm_gni_zct = zeroc

              IF(ALLOCATED(mmm_gq_zct))DEALLOCATE(mmm_gq_zct)
              ALLOCATE(mmm_gq_zct(njin-1))
              mmm_gq_zct = zeroc

              IF(ALLOCATED(mmm_rmaj_ctr_zct))DEALLOCATE(mmm_rmaj_ctr_zct)
              ALLOCATE(mmm_rmaj_ctr_zct(njin-1))
              mmm_rmaj_ctr_zct = zeroc

              IF(ALLOCATED(mmm_drhodrmin_zct))DEALLOCATE(mmm_drhodrmin_zct)
              ALLOCATE(mmm_drhodrmin_zct(njin-1))
              mmm_drhodrmin_zct = zeroc

              IF(ALLOCATED(mmm_rmin_zct))DEALLOCATE(mmm_rmin_zct)
              ALLOCATE(mmm_rmin_zct(njin-1))
              mmm_rmin_zct = zeroc

              IF(ALLOCATED(mmm_elong_zct))DEALLOCATE(mmm_elong_zct)
              ALLOCATE(mmm_elong_zct(njin-1))
              mmm_elong_zct = zeroc

              IF(ALLOCATED(mmm_zeff_zct))DEALLOCATE(mmm_zeff_zct)
              ALLOCATE(mmm_zeff_zct(njin-1))
              mmm_zeff_zct = zeroc

              IF(ALLOCATED(mmm_ne_zct))DEALLOCATE(mmm_ne_zct)
              ALLOCATE(mmm_ne_zct(njin-1))
              mmm_ne_zct = zeroc

              IF(ALLOCATED(mmm_nh_zct))DEALLOCATE(mmm_nh_zct)
              ALLOCATE(mmm_nh_zct(njin-1))
              mmm_nh_zct = zeroc

              IF(ALLOCATED(mmm_ni_zct))DEALLOCATE(mmm_ni_zct)
              ALLOCATE(mmm_ni_zct(njin-1))
              mmm_ni_zct = zeroc

              IF(ALLOCATED(mmm_nz_zct))DEALLOCATE(mmm_nz_zct)
              ALLOCATE(mmm_nz_zct(njin-1))
              mmm_nz_zct = zeroc

              IF(ALLOCATED(mmm_nf_zct))DEALLOCATE(mmm_nf_zct)
              ALLOCATE(mmm_nf_zct(njin-1))
              mmm_nf_zct = zeroc

              IF(ALLOCATED(mmm_te_zct))DEALLOCATE(mmm_te_zct)
              ALLOCATE(mmm_te_zct(njin-1))
              mmm_te_zct = zeroc

              IF(ALLOCATED(mmm_ti_zct))DEALLOCATE(mmm_ti_zct)
              ALLOCATE(mmm_ti_zct(njin-1))
              mmm_ti_zct = zeroc


              IF(ALLOCATED(mmm_q_zct))DEALLOCATE(mmm_q_zct)
              ALLOCATE(mmm_q_zct(njin-1))
              mmm_q_zct = zeroc

              IF(ALLOCATED(mmm_btor_zct))DEALLOCATE(mmm_btor_zct)
              ALLOCATE(mmm_btor_zct(njin-1))
              mmm_btor_zct = zeroc

              IF(ALLOCATED(mmm_bpol_zct))DEALLOCATE(mmm_bpol_zct)
              ALLOCATE(mmm_bpol_zct(njin-1))
              mmm_bpol_zct = zeroc

              IF(ALLOCATED(mmm_zimp_zct))DEALLOCATE(mmm_zimp_zct)
              ALLOCATE(mmm_zimp_zct(njin-1))
              mmm_zimp_zct = zeroc

              IF(ALLOCATED(mmm_aimp_zct))DEALLOCATE(mmm_aimp_zct)
              ALLOCATE(mmm_aimp_zct(njin-1))
              mmm_aimp_zct = zeroc

              IF(ALLOCATED(mmm_ahyd_zct))DEALLOCATE(mmm_ahyd_zct)
              ALLOCATE(mmm_ahyd_zct(njin-1))
              mmm_ahyd_zct = zeroc

              IF(ALLOCATED(mmm_aimass_zct))DEALLOCATE(mmm_aimass_zct)
              ALLOCATE(mmm_aimass_zct(njin-1))
              mmm_aimass_zct  = zeroc


              IF(ALLOCATED(mmm_wexbs_zct))DEALLOCATE(mmm_wexbs_zct)
              ALLOCATE(mmm_wexbs_zct(njin-1))
              mmm_wexbs_zct = zeroc


              IF(ALLOCATED(mmm_vpol_zct))DEALLOCATE(mmm_vpol_zct)
              ALLOCATE(mmm_vpol_zct(njin-1)) 
              mmm_vpol_zct = zeroc

              IF(ALLOCATED(mmm_gvpol_zct))DEALLOCATE(mmm_gvpol_zct)
              ALLOCATE(mmm_gvpol_zct(njin-1))
              mmm_gvpol_zct = zeroc

              IF(ALLOCATED(mmm_vpar_zct))DEALLOCATE(mmm_vpar_zct)
              ALLOCATE(mmm_vpar_zct(njin-1))
              mmm_vpar_zct= zeroc

              IF(ALLOCATED(mmm_gvpar_zct))DEALLOCATE(mmm_gvpar_zct)
              ALLOCATE(mmm_gvpar_zct(njin-1))
              mmm_gvpar_zct= zeroc

              IF(ALLOCATED(mmm_vtor_zct))DEALLOCATE(mmm_vtor_zct)
              ALLOCATE(mmm_vtor_zct(njin-1))
              mmm_vtor_zct= zeroc

              IF(ALLOCATED(mmm_gvtor_zct))DEALLOCATE(mmm_gvtor_zct)
              ALLOCATE(mmm_gvtor_zct(njin-1))
              mmm_gvtor_zct= zeroc

              IF(ALLOCATED(mmm_gelong_zct))DEALLOCATE(mmm_gelong_zct)
              ALLOCATE(mmm_gelong_zct(njin-1))
              mmm_gelong_zct= zeroc

              IF(ALLOCATED(mmm_er_zct))DEALLOCATE(mmm_er_zct)
              ALLOCATE(mmm_er_zct(njin-1))
              mmm_er_zct = zeroc

              IF(ALLOCATED(mmm_theta_zct))DEALLOCATE(mmm_theta_zct)
              ALLOCATE(mmm_theta_zct(njin-1))
              mmm_theta_zct = zeroc


              IF(ALLOCATED(mmm_vdia_zct))DEALLOCATE(mmm_vdia_zct)
              ALLOCATE(mmm_vdia_zct(njin-1))
              mmm_vdia_zct= zeroc


              IF(ALLOCATED(mmm_vexb_zct))DEALLOCATE(mmm_vexb_zct)
              ALLOCATE(mmm_vexb_zct(njin-1))
              mmm_vexb_zct = zeroc

              IF(ALLOCATED(mmm_tor_ang_speed_zct))DEALLOCATE(mmm_tor_ang_speed_zct)
              ALLOCATE(mmm_tor_ang_speed_zct(njin-1))
              mmm_tor_ang_speed_zct = zeroc

              IF(ALLOCATED(mmm_angrot_zct))DEALLOCATE(mmm_angrot_zct)
              ALLOCATE(mmm_angrot_zct(njin-1))
              mmm_angrot_zct = zeroc


              IF(ALLOCATED(mmm_gangrot_zct))DEALLOCATE(mmm_gangrot_zct)
              ALLOCATE(mmm_gangrot_zct(njin-1))
              mmm_gangrot_zct = zeroc


              !IF(select_solver == 'nwt_pred_cor')THEN
              ! must be defined for message passing even for maccormack solver
                 IF(ALLOCATED(mmm_val_pert_zct))DEALLOCATE(mmm_val_pert_zct)
                 ALLOCATE(mmm_val_pert_zct(-itran_max:itran_max,0:njin-1))
                 mmm_val_pert_zct = zeroc
                 IF(ALLOCATED(mmm_val_base_zct))DEALLOCATE(mmm_val_base_zct)
                 ALLOCATE(mmm_val_base_zct(-itran_max:itran_max,0:njin-1))
                 mmm_val_base_zct = zeroc
              !ENDIF




! OUTPUT:
              IF(ALLOCATED(mmm_xti_ze))DEALLOCATE( mmm_xti_ze)
              ALLOCATE(mmm_xti_ze(njin))
              mmm_xti_ze = zeroc

              IF(ALLOCATED(mmm_xdi_ze))DEALLOCATE(mmm_xdi_ze)
              ALLOCATE(mmm_xdi_ze(njin))
              mmm_xdi_ze = zeroc

              IF(ALLOCATED(mmm_xte_ze))DEALLOCATE( mmm_xte_ze)
              ALLOCATE(mmm_xte_ze(njin))
              mmm_xte_ze = zeroc

              IF(ALLOCATED(mmm_xdz_ze))DEALLOCATE(mmm_xdz_ze)
              ALLOCATE(mmm_xdz_ze(njin))
              mmm_xdz_ze = zeroc

              IF(ALLOCATED(mmm_xvt_ze))DEALLOCATE( mmm_xvt_ze)
              ALLOCATE(mmm_xvt_ze(njin))
              mmm_xvt_ze = zeroc

              IF(ALLOCATED(mmm_xvp_ze))DEALLOCATE( mmm_xvp_ze)
              ALLOCATE(mmm_xvp_ze(njin))
              mmm_xvp_ze = zeroc

              IF(ALLOCATED(mmm_gammaDBM_ze))DEALLOCATE(mmm_gammaDBM_ze)
              ALLOCATE(mmm_gammaDBM_ze(njin)) 
              mmm_gammaDBM_ze = zeroc

              IF(ALLOCATED(mmm_omegaDBM_ze))DEALLOCATE(mmm_omegaDBM_ze)   
              ALLOCATE(mmm_omegaDBM_ze(njin)) 
              mmm_omegaDBM_ze = zeroc
 
              IF(ALLOCATED(mmm_xtiW20_ze))DEALLOCATE(mmm_xtiW20_ze) 
              ALLOCATE(mmm_xtiW20_ze(njin))
              mmm_xtiW20_ze = zeroc

              IF(ALLOCATED(mmm_xdiW20_ze))DEALLOCATE(mmm_xdiW20_ze)
              ALLOCATE(mmm_xdiW20_ze(njin))
              mmm_xdiW20_ze= zeroc
 
              IF(ALLOCATED(mmm_xteW20_ze))DEALLOCATE( mmm_xteW20_ze)
              ALLOCATE(mmm_xteW20_ze(njin))
              mmm_xteW20_ze = zeroc

              IF(ALLOCATED(mmm_xtiDBM_ze))DEALLOCATE(mmm_xtiDBM_ze)
              ALLOCATE(mmm_xtiDBM_ze(njin)) 
              mmm_xtiDBM_ze = zeroc

              IF(ALLOCATED(mmm_xdiDBM_ze))DEALLOCATE(mmm_xdiDBM_ze)
              ALLOCATE(mmm_xdiDBM_ze(njin)) 
              mmm_xdiDBM_ze = zeroc

              IF(ALLOCATED(mmm_xteDBM_ze))DEALLOCATE( mmm_xteDBM_ze)
              ALLOCATE(mmm_xteDBM_ze(njin))
              mmm_xteDBM_ze = zeroc

              IF(ALLOCATED(mmm_xteETG_ze))DEALLOCATE( mmm_xteETG_ze) 
              ALLOCATE(mmm_xteETG_ze(njin)) 
              mmm_xteETG_ze = zeroc
         
              IF(ALLOCATED(mmm_gammaW20_ze))DEALLOCATE(mmm_gammaW20_ze) 
              ALLOCATE(mmm_gammaW20_ze(4,njin))           
              mmm_gammaW20_ze  = zeroc
         
              IF(ALLOCATED(mmm_omegaW20_ze))DEALLOCATE(mmm_omegaW20_ze) 
              ALLOCATE(mmm_omegaW20_ze(4,njin)) 
              mmm_omegaW20_ze = zeroc

              IF(ALLOCATED(mmm_vconv_ze))DEALLOCATE(mmm_vconv_ze) 
              ALLOCATE(mmm_vconv_ze(6,njin)) 
              mmm_vconv_ze = zeroc


              IF(ALLOCATED(mmm_vflux_ze))DEALLOCATE(mmm_vflux_ze) 
              ALLOCATE(mmm_vflux_ze(4,njin))
              mmm_vflux_ze = zeroc

              IF(ALLOCATED(mmm_zchrgz_ze))DEALLOCATE(mmm_zchrgz_ze)
              ALLOCATE(mmm_zchrgz_ze(njin))
              mmm_zchrgz_ze = zeroc

              IF(ALLOCATED(mmm_zchrgf_ze))DEALLOCATE( mmm_zchrgf_ze)
              ALLOCATE(mmm_zchrgf_ze(njin))
              mmm_zchrgf_ze = zeroc

              IF(ALLOCATED(mmm_zxb_ze))DEALLOCATE(mmm_zxb_ze) 
              ALLOCATE(mmm_zxb_ze(njin))
              mmm_zxb_ze = zeroc

              IF(ALLOCATED(mmm_zdensf_ze))DEALLOCATE( mmm_zdensf_ze)
              ALLOCATE(mmm_zdensf_ze(njin))
              mmm_zdensf_ze = zeroc

              IF(ALLOCATED(mmm_zdensi_ze))DEALLOCATE(mmm_zdensi_ze) 
              ALLOCATE(mmm_zdensi_ze(njin))
              mmm_zdensi_ze = zeroc

              IF(ALLOCATED(mmm_zgrdnf_ze))DEALLOCATE(mmm_zgrdnf_ze)
              ALLOCATE(mmm_zgrdnf_ze(njin))
              mmm_zgrdnf_ze = zeroc



              IF(ALLOCATED(mmm_xti_zct))DEALLOCATE( mmm_xti_zct)
              ALLOCATE(mmm_xti_zct(njin-1))
              mmm_xti_zct = zeroc

              IF(ALLOCATED(mmm_xdi_zct))DEALLOCATE(mmm_xdi_zct)
              ALLOCATE(mmm_xdi_zct(njin-1))
              mmm_xdi_zct = zeroc

              IF(ALLOCATED(mmm_xte_zct))DEALLOCATE( mmm_xte_zct)
              ALLOCATE(mmm_xte_zct(njin-1))
              mmm_xte_zct = zeroc

              IF(ALLOCATED(mmm_xdz_zct))DEALLOCATE(mmm_xdz_zct)
              ALLOCATE(mmm_xdz_zct(njin-1))
              mmm_xdz_zct = zeroc

              IF(ALLOCATED(mmm_xvt_zct))DEALLOCATE( mmm_xvt_zct)
              ALLOCATE(mmm_xvt_zct(njin-1))
              mmm_xvt_zct = zeroc

              IF(ALLOCATED(mmm_xvp_zct))DEALLOCATE( mmm_xvp_zct)
              ALLOCATE(mmm_xvp_zct(njin-1))
              mmm_xvp_zct = zeroc

              IF(ALLOCATED(mmm_gammaDBM_zct))DEALLOCATE(mmm_gammaDBM_zct)
              ALLOCATE(mmm_gammaDBM_zct(njin-1)) 
              mmm_gammaDBM_zct = zeroc

              IF(ALLOCATED(mmm_omegaDBM_zct))DEALLOCATE(mmm_omegaDBM_zct)   
              ALLOCATE(mmm_omegaDBM_zct(njin-1)) 
              mmm_omegaDBM_zct = zeroc
 
              IF(ALLOCATED(mmm_xtiW20_zct))DEALLOCATE(mmm_xtiW20_zct) 
              ALLOCATE(mmm_xtiW20_zct(njin-1))
              mmm_xtiW20_zct = zeroc

              IF(ALLOCATED(mmm_xdiW20_zct))DEALLOCATE(mmm_xdiW20_zct)
              ALLOCATE(mmm_xdiW20_zct(njin-1))
              mmm_xdiW20_zct= zeroc
 
              IF(ALLOCATED(mmm_xteW20_zct))DEALLOCATE( mmm_xteW20_zct)
              ALLOCATE(mmm_xteW20_zct(njin-1))
              mmm_xteW20_zct = zeroc

              IF(ALLOCATED(mmm_xtiDBM_zct))DEALLOCATE(mmm_xtiDBM_zct)
              ALLOCATE(mmm_xtiDBM_zct(njin-1)) 
              mmm_xtiDBM_zct = zeroc

              IF(ALLOCATED(mmm_xdiDBM_zct))DEALLOCATE(mmm_xdiDBM_zct)
              ALLOCATE(mmm_xdiDBM_zct(njin-1)) 
              mmm_xdiDBM_zct = zeroc

              IF(ALLOCATED(mmm_xteDBM_zct))DEALLOCATE( mmm_xteDBM_zct)
              ALLOCATE(mmm_xteDBM_zct(njin-1))
              mmm_xteDBM_zct = zeroc

              IF(ALLOCATED(mmm_xteETG_zct))DEALLOCATE( mmm_xteETG_zct) 
              ALLOCATE(mmm_xteETG_zct(njin-1)) 
              mmm_xteETG_zct = zeroc
         
              IF(ALLOCATED(mmm_gammaW20_zct))DEALLOCATE(mmm_gammaW20_zct) 
              ALLOCATE(mmm_gammaW20_zct(4,njin-1))           
              mmm_gammaW20_zct  = zeroc
         
              IF(ALLOCATED(mmm_omegaW20_zct))DEALLOCATE(mmm_omegaW20_zct) 
              ALLOCATE(mmm_omegaW20_zct(4,njin-1)) 
              mmm_omegaW20_zct = zeroc

              IF(ALLOCATED(mmm_vconv_zct))DEALLOCATE(mmm_vconv_zct) 
              ALLOCATE(mmm_vconv_zct(6,njin-1)) 
              mmm_vconv_zct = zeroc


              IF(ALLOCATED(mmm_vflux_zct))DEALLOCATE(mmm_vflux_zct) 
              ALLOCATE(mmm_vflux_zct(4,njin-1))
              mmm_vflux_zct = zeroc

              IF(ALLOCATED(mmm_zchrgz_zct))DEALLOCATE(mmm_zchrgz_zct)
              ALLOCATE(mmm_zchrgz_zct(njin-1))
              mmm_zchrgz_zct = zeroc

              IF(ALLOCATED(mmm_zchrgf_zct))DEALLOCATE( mmm_zchrgf_zct)
              ALLOCATE(mmm_zchrgf_zct(njin-1))
              mmm_zchrgf_zct = zeroc

              IF(ALLOCATED(mmm_zxb_zct))DEALLOCATE(mmm_zxb_zct) 
              ALLOCATE(mmm_zxb_zct(njin-1))
              mmm_zxb_zct = zeroc

              IF(ALLOCATED(mmm_zdensf_zct))DEALLOCATE( mmm_zdensf_zct)
              ALLOCATE(mmm_zdensf_zct(njin-1))
              mmm_zdensf_zct = zeroc

              IF(ALLOCATED(mmm_zdensi_zct))DEALLOCATE(mmm_zdensi_zct) 
              ALLOCATE(mmm_zdensi_zct(njin-1))
              mmm_zdensi_zct = zeroc

              IF(ALLOCATED(mmm_zgrdnf_zct))DEALLOCATE(mmm_zgrdnf_zct)
              ALLOCATE(mmm_zgrdnf_zct(njin-1))
              mmm_zgrdnf_zct = zeroc





!------------------------------------------------------------------
!------------------------------------------------------------------
              IF(ALLOCATED(mmm_dcoef_zc))DEALLOCATE(mmm_dcoef_zc)
              ALLOCATE(mmm_dcoef_zc(ntot+1,njin-1))  ! zone ctr grid, ntot + electrons
              mmm_dcoef_zc = zeroc

              IF(ALLOCATED(mmm_tot_angmtm_zc))DEALLOCATE(mmm_tot_angmtm_zc)
              ALLOCATE(mmm_tot_angmtm_zc(njin-1))
              mmm_tot_angmtm_zc = zeroc

              IF(ALLOCATED(mmm_fluxeth_zc))DEALLOCATE(mmm_fluxeth_zc)
             ALLOCATE(mmm_fluxeth_zc(njin-1))
             mmm_fluxeth_zc = zeroc

              IF(ALLOCATED(mmm_fluxith_zc))DEALLOCATE(mmm_fluxith_zc)
              ALLOCATE(mmm_fluxith_zc(njin-1))
              mmm_fluxith_zc = zeroc

             IF(ALLOCATED(mmm_fluxipt_zc))DEALLOCATE(mmm_fluxipt_zc)
             ALLOCATE(mmm_fluxipt_zc(njin-1,nion))
             mmm_fluxipt_zc = zeroc

! mpi related:
              !------------------------------------------------------------------
              !In sequence , MUST BE COORDINATED WITH mmm_pack_data and 
              !mmm_store_packed
              !------------------------------------------------------------------- 
              mmm_max_pack = 1   ! grid pt number
              mmm_max_pack = mmm_max_pack + 11       ! DRB and Weiland growth rates and freq
              mmm_max_pack = mmm_max_pack + 4        ! parrticle,energy fluxes
              mmm_max_pack = mmm_max_pack + 6        ! convective velocities
              mmm_max_pack = mmm_max_pack + 1        ! t_slave_compt
              mmm_max_pack = mmm_max_pack + 11       ! -itran_max,itran_max,mmm_val_pert
              mmm_max_pack = mmm_max_pack + 11       ! -itran_max,itran_max,mmm_base_pert
              mmm_max_pack = mmm_max_pack + 1        ! xti
              mmm_max_pack = mmm_max_pack + 1        ! xdi
              mmm_max_pack = mmm_max_pack + 1        ! xte
              mmm_max_pack = mmm_max_pack + 1        ! xdz
              mmm_max_pack = mmm_max_pack + 1        ! xvt
              mmm_max_pack = mmm_max_pack + 1        ! xvp
              mmm_max_pack = mmm_max_pack + 1        ! xtiW20
              mmm_max_pack = mmm_max_pack + 1        ! xdiW20
              mmm_max_pack = mmm_max_pack + 1        ! xteW20
              mmm_max_pack = mmm_max_pack + 1        ! xtiDBM
              mmm_max_pack = mmm_max_pack + 1        ! xdiDBM
              mmm_max_pack = mmm_max_pack + 1        ! xteDBM
              mmm_max_pack = mmm_max_pack + 1        ! xteETG

              
              mmm_max_pack = mmm_max_pack + 3*nspecies_mmm*(2*itran_max+1) ! mmm_pflux,mmm_eflux,mmm_m_flux        


!              mmm_max_pack = mmm_max_pack + ntot+1   ! mmm_dcoef_zc
              IF(ALLOCATED(mmm_send_packed))DEALLOCATE(mmm_send_packed)
              IF(ALLOCATED(t_slave_compt_tot))DEALLOCATE(t_slave_compt_tot)
              IF(ALLOCATED(n_slave_grid_tot))DEALLOCATE(n_slave_grid_tot)
              IF(myid .NE. master)THEN
                 ALLOCATE(mmm_send_packed(mmm_max_pack))
              ELSE
                 ALLOCATE(n_slave_grid_tot(0:numprocs-1),t_slave_compt_tot(0:numprocs-1))
                 t_slave_compt_tot(:)    = zeroc 
                 n_slave_grid_tot(:)     = izero
              ENDIF

              IF(ALLOCATED(mmm_receive_packed))DEALLOCATE(mmm_receive_packed)
              IF(myid .EQ. master)ALLOCATE(mmm_receive_packed(mmm_max_pack))


              RETURN

         END SUBROUTINE mmm_allocate

 
   SUBROUTINE  mmm_flux_processing(njoldin,njnew,array,k,ksp,kit,nsp,out)
!------------------------------------------------------------------------------------------
! -- prepare mmm  flux values for output (used if use_mmm_flux =0)
! -- note: array can be mmm_p_flux,mmm_e_flux or mmm_m_flux because of
! -- assumed structure.
! -- we take zone center values of  array(1:njoldin-1,nsp,0)  and convert to zone
! -- edge values  of size(1:njnew) in array out(1:njnew):
! -- INPUT: 
! --   array(k,ksp,kit)
! --   njoldin,njnew
! -- OUTPUT:
! --   out(1:njnew)
!---------------------------------------------------------------------------HSJ-----------
      USE nrtype,                                             ONLY : DP,I4B

      USE Grid_class,                                         ONLY : rnew_zc,rold_zc, &
                                                                     rold,r

      USE solcon_gcnmp,                                       ONLY : use_mmm_flux


      IMPLICIT NONE 
    
      INTEGER(I4B) j,k,ksp,kit,nsp,njnewm1,njnew, &
                   jold,njoldinm1,njoldin,njold,njoldm1
!      REAL(DP),DIMENSION(k,ksp, kit:-kit)  :: array
      REAL(DP),DIMENSION(k,ksp, -kit:kit)  :: array  ! 999,
      REAL(DP),DIMENSION(njnew)      :: out   ! temp array
      REAL(DP) df,drold,drx

       njnewm1 = njnew-1  ; njoldinm1 = njoldin -1
       IF(ALLOCATED(rold_zc))DEALLOCATE(rold_zc)
       IF(ALLOCATED(rnew_zc))DEALLOCATE(rnew_zc)
       ALLOCATE(rold_zc(njoldinm1))
       ALLOCATE(rnew_zc(njnewm1))


       DO j=1,njnewm1
           rnew_zc(j) = (r(j+1)+r(j))*0.5_DP
       ENDDO


      njold = SIZE(rold)
      njoldm1 = njold-1

       DO j=1,njoldm1
           rold_zc(j) = (rold(j+1)+rold(j))*0.5_DP
       ENDDO

          DO j=1,njnew-1
             DO jold =1,njoldm1-1
                  IF(rnew_zc(j) .LT. rold_zc(1))THEN
                        !take flux as zero at rho =0 so no extrapolation necessary:
                        df =  array(1,nsp,0)      ! j=1,first zc pt,nsp = species, 0 => pertrb =0
                        drold = rold_zc(1)
                        drx =  rnew_zc(j) 
                        out(j) = (df/drold)*drx 
                  ELSEIF(rnew_zc(j)  .GE. rold_zc(jold).AND. rnew_zc(j) .LE. rold_zc(jold+1))THEN
                       drx = rnew_zc(j) - rold_zc(jold)
                       df  = array(jold+1,nsp,0)- array(jold,nsp,0) ! [1 to njoldm1] 
                       drold = rold_zc(jold+1)- rold_zc(jold)                   ! [1 to njoldm1] 
                       out(j) = array(jold,nsp,0)  + drx*df/drold               ! [1 to njnew-1]
  
                  ELSEIF(rnew_zc(j) .GT. rold_zc(njoldm1))THEN 
                    ! value requires extrapolation since edge flux value is not known in this version of code
                    df = array(njoldm1,nsp,0)- array(njoldm1-1,nsp,0)
                    drold = rold_zc(njoldm1)- rold_zc(njoldm1-1)
                    drx =  rnew_zc(j) -  rold_zc(njoldm1)
                    out(j) = array(njoldm1,nsp,0)  + (df/drold)*drx

                 ENDIF
               ENDDO

           ENDDO

   END SUBROUTINE  mmm_flux_processing


         SUBROUTINE mmm_init
         !--------------------------------------------------------------
         ! -- Set those quantities that dont change in time here:
         ! -- dischg% ==> Set in plasma_contours
         !--------------------------------------------------------------
               USE Plasma_properties,              ONLY : dischg

               USE Vector_class,                   ONLY : get_values

               USE grid_class,                     ONLY : roa,r,nj

               USE curden_terms,                   ONLY : q

               USE plasma_contours,                ONLY : bp_outbd,btor_outbd

               USE common_constants,               ONLY : izero,zeroc

               IMPLICIT NONE
               INTEGER(i4B) j


               IF(SIZE(mmm_rmin_ze) .NE. dischg%rmin_half_width%size)THEN
                  DEALLOCATE(mmm_rmin_ze,mmm_rmaj_ctr_ze,mmm_elong_ze)
                  ALLOCATE(mmm_rmin_ze(dischg%rmin_half_width%size),   &
                  mmm_rmaj_ctr_ze(dischg%rmin_half_width%size),        &
                  mmm_elong_ze(dischg%rmin_half_width%size))
               ENDIF

               IF(SIZE(mmm_rmin_zct) .NE. dischg%rmin_half_width%size-1)THEN
                  DEALLOCATE(mmm_rmin_zct,mmm_rmaj_ctr_zct,mmm_elong_zct)
                  ALLOCATE(mmm_rmin_zct(dischg%rmin_half_width%size-1),   &
                  mmm_rmaj_ctr_zct(dischg%rmin_half_width%size-1),        &
                  mmm_elong_zct(dischg%rmin_half_width%size-1))
               ENDIF


               mmm_rmin_ze(:)      = get_values(dischg%rmin_half_width)
               mmm_rmaj_ctr_ze(:)  = get_values(dischg%rmaj_geom) ! rmajor to center of flux surface
               mmm_elong_ze(:)     = get_values(dischg%elongxnj)

               DO j=1,dischg%rmin_half_width%size-1
                  mmm_rmin_zct(j)      = 0.5_DP*(mmm_rmin_ze(j+1) + mmm_rmin_ze(j))
                  mmm_elong_zct(j)     = 0.5_DP*(mmm_elong_ze(j+1) + mmm_elong_ze(j))
                  mmm_rmaj_ctr_zct(j)  = 0.5_DP*(mmm_rmaj_ctr_ze(j+1) + mmm_rmaj_ctr_ze(j))
                  mmm_theta_zct        = (bp_outbd(j+1)+bp_outbd(j+1))/ &
                                         (btor_outbd(j+1)+btor_outbd(j))
               ENDDO

               mmm_rmaj0           = dischg%rma
               !mmm_theta_ze(:)    = r(nj)*roa(:)/(ABS(q(:)*dischg%rmajor)) ! = Bp0/Bt0
               mmm_theta_ze(:)     = bp_outbd(:)/btor_outbd(:)


 
               mmm_drhodrmin_ze(1) = zeroc
               DO j=1,nj-1
                  IF(j .GT. 1) &
                  mmm_drhodrmin_ze(j)    =   (r(j+1)-r(j-1))/(mmm_rmin_ze(j+1)-mmm_rmin_ze(j-1))
                  mmm_drhodrmin_zct(j)    =  (r(j+1)-r(j))/(mmm_rmin_ze(j+1)-mmm_rmin_ze(j))
               ENDDO


               mmm_initialized     = .TRUE.

               RETURN

         END SUBROUTINE mmm_init

         SUBROUTINE  set_mmm_vars
         !--------------------------------------------------------------------------
         ! -- setup input variables for mmm
         ! -- NOTE that these are defined on the outboard side in terms of half width 
         ! -- of each flux surface, as defined by dischg%rmin_half_width
         ! -- set variables at zone edges,** _ze
         !--------------------------------------------------------------------------
           USE grid_class,                           ONLY : nj,fcap,gcap,hcap,       &
                                                            use_compact_schemes,     &
                                                            psir_grid,r

           USE common_constants,                     ONLY : izero,zeroc,joupkev,     &
                                                            Electron_Charge

           USE error_handler,                        ONLY : iomaxerr, lerrno,        &
                                                            terminate
           USE io_gcnmp,                             ONLY : nlog, ncrt
      
           USE dep_var,                              ONLY : rbp,te,ti,angrot,en,     &
                                                            ene,vphi,vpol_nclass,    &
                                                            vpol,vpar_nclass,        &
                                                            er_tot_nclass

           USE ions_gcnmp,                           ONLY : nion,nprim,zeff,z,       &
                                                            atw
  
           USE neutral_beams,                        ONLY : nbion,enbeam_tot

           USE curden_terms,                         ONLY : q,bp

           USE Plasma_properties,                    ONLY : mhd_dat,dischg,profile

           USE solcon_gcnmp,                         ONLY : use_mmm_vpol_input,      &
                                                            use_mmm_vpar_input,      &
                                                            use_mmm_er_input,        &
                                                            use_forcebal

           USE Vector_class,                         ONLY : get_values

           USE plasma_contours,                      ONLY : bp_outbd,rmaj_outbd,btor_outbd
           USE MPI_data,              ONLY : master,myid,numprocs ! debug only ! 888889999

           IMPLICIT NONE
           
           INTEGER(I4B) jr
           REAL(DP), ALLOCATABLE,DIMENSION(:) :: s_den
           REAL(DP) alpha_dia,bt,bpsq,fc,eps, cfctr,grad_press,grad_ti,grad_den, &
                    K1,bps,dpsi,btot,grad_denni,grad_te,grad_angrot,grad_ne
!------------ temp settings for test of compact deriv method below
           LOGICAL    done
           INTEGER(I4B)inpt_dpl, inpt_dpr,nrhs,zctr_deriv, msp ,order
           REAL(DP) dpl,dpr
           REAL(DP) dpdg(nj)
              inpt_dpl = 1  ; dpl = zeroc             ! all profiles have dudr =0 at r =0
              inpt_dpr = izero                        ! no derivative condition at r = a 
              nrhs = 1 ! nrhs > 1 not implemented 


             !  use sixth order uniform grid spacing if possible:

             zctr_deriv = 0  ! get deriv on zone faces
             msp = 1                                                      ! note rmin_fg is non uniform grid
             order =6      
             done = .false.
!---------- end temp settings
               IF( .NOT. mmm_initialized)CALL mmm_init

               IF(ALLOCATED(s_den))DEALLOCATE(s_den)
               ALLOCATE(s_den(SIZE(mmm_ne_ze)))
               
               alpha_dia = 1._DP ! multiplier ,assign later if wanted
               npoints = nj
    
               mmm_gte_ze(1)   = zeroc
               mmm_gti_ze(1)   = zeroc
               mmm_gne_ze(1)   = zeroc
               mmm_gnh_ze(1)   = zeroc
               mmm_gni_ze(1)   = zeroc
               mmm_gnz_ze(1)   = zeroc
               mmm_gq_ze(1)    = zeroc      !  R ( d q   / d r ) / q

               mmm_ne_ze(:)    = ene(:)
               mmm_nf_ze(:)    = enbeam_tot(:)
               mmm_nh_ze(:)    = zeroc
               mmm_ni_ze(:)    = zeroc
               DO jr = 1, nprim
                  mmm_nh_ze(:)    = mmm_nh_ze(:) +  en(:,jr)
               ENDDO

               mmm_te_ze(:)    = te(:)
               mmm_ti_ze(:)    = ti(:)
               mmm_zeff_ze(:)  = zeff(:)
               mmm_q_ze(:)     = ABS(q(:))
               !mmm_btor_ze(:) = get_values(mhd_dat%fcap)
               !mmm_btor_ze(:) = mmm_btor_ze(:)*dischg%rmajor*mhd_dat%btor/rmaj_outbd(:)
               mmm_btor_ze(:)  = btor_outbd(:)
               mmm_bpol_ze(:)  = bp_outbd(:)
               mmm_zimp_ze(:)  = zeroc
               s_den(:)        = zeroc
               mmm_aimp_ze(:)  = zeroc
               DO jr = nprim+1,nion
                  mmm_zimp_ze(:)  = mmm_zimp_ze(:) + en(:,jr)*z(:,jr)
                  s_den(:)       = s_den(:) + en(:,jr)
                  mmm_aimp_ze(:)  = mmm_aimp_ze(:) + en(:,jr)*atw(jr)
               ENDDO
               mmm_nz_ze(:)       = s_den(:)
               mmm_zimp_ze(:)     = mmm_zimp_ze(:)/s_den(:)
               mmm_aimp_ze(:)     = mmm_aimp_ze(:)/s_den(:)
               mmm_ahyd_ze(:)     = zeroc
               s_den(:)           = zeroc

               DO jr = 1,nprim
                  mmm_ahyd_ze(:)  = mmm_ahyd_ze(:) + en(:,jr)*atw(jr)
                  s_den(:)        = s_den(:) + en(:,jr)
               ENDDO
               mmm_aimass_ze(:)   = mmm_ahyd_ze(:)
               mmm_ahyd_ze(:)     = mmm_ahyd_ze(:)/s_den(:)
               
               DO jr =nprim+1,nion
                  mmm_aimass_ze(:)= mmm_aimass_ze(:) + en(:,jr)*atw(jr)
                  s_den(:)        = s_den(:) + en(:,jr)
               ENDDO
               mmm_aimass_ze(:)   = mmm_aimass_ze(:)/s_den(:)

               mmm_ni_ze(:)       = mmm_nz_ze(:) + mmm_nh_ze(:)  ! total thermal ions density

          cfctr              = -joupkev/Electron_Charge              ! Electron_Charge is negative
          Bt                 = mhd_dat%btor
          grad_press         = zeroc ; grad_ti = zeroc ; grad_den = zeroc
          vphi(:)            = dischg%rmajor*angrot(:)
          mmm_angrot_ze(:)   = angrot(:)
          mmm_vtor_ze(:)     = vphi(:)


!write(940+myid,FMT='("use_mmm_vpol_input,use_mmm_vpar_input,  =",x,i5,x,i5)')use_mmm_vpol_input,use_mmm_vpar_input
!write(940+myid,FMT='("use_forcebal =",x,i5)')use_forcebal


          ! need to set here because we look forward (j+1) in Do loop"
          IF(use_mmm_vpol_input == 0 )THEN

                IF(use_forcebal == 1)THEN ! Nclass formulation
                     mmm_vpol_ze(:)  = vpol_nclass(:)
                ELSE
                   DO jr=2,npoints-1 ! Y.B.KIM formulation 
                       grad_ti    = ( mmm_ti_ze(jr+1)-mmm_ti_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) ) 
                       eps        = mmm_rmin_ze(jr)/mmm_rmaj_ctr_ze(jr)
                       fc         = 1.00_DP +(0.46_DP*eps -1.46_DP)*SQRT(eps) ! trapped partcle farc
                       K1         = 0.8839_DP*fc/(0.3477_DP+0.4058_DP*fc)
                       Bpsq       = mmm_bpol_ze(jr)*mmm_bpol_ze(jr) 
                       mmm_vpol_ze(jr) = alpha_dia*K1*cfctr*grad_ti*Bt/(Bt*Bt + Bpsq) ! m/sec (1,nj)
                   ENDDO

                     mmm_vpol_ze(1) = zeroc
                     grad_ti    = ( mmm_ti_ze(npoints)-mmm_ti_ze(npoints-1) )/ &
                                       ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) ) 
                     eps        = mmm_rmin_ze(npoints)/mmm_rmaj_ctr_ze(npoints)
                     fc         = 1.00_DP +(0.46_DP*eps -1.46_DP)*SQRT(eps) ! trapped partcle farc
                     K1         = 0.8839_DP*fc/(0.3477_DP+0.4058_DP*fc)
                     Bpsq       = mmm_bpol_ze(npoints)*mmm_bpol_ze(npoints) 
                     mmm_vpol_ze(npoints) = alpha_dia*K1*cfctr*grad_ti*Bt/(Bt*Bt + Bpsq) ! m/sec (1,nj)
                ENDIF

          ELSE           ! use vpol or vpol_nclass,from statefile
                         ! in this case these values dont change as the profiles evolve in time
               IF(use_forcebal == 1)THEN
                 !profile%vpol_nclass (from statefile) was put  on current grid in refactor_profiles
                 !Note that use_forcebal here means use profile%vpol_nclass%data
                 ! instead of profile%vpol%data but the two may infact be the same
                 mmm_vpol_ze(:) = profile%vpol_nclass%data(:)
               ELSE
                 ! put profile%vpol (from state file) on current grid
                 mmm_vpol_ze(:) = profile%vpol%data(:)

               ENDIF

          ENDIF

          IF(use_mmm_vpar_input == 0 )THEN

                IF(use_forcebal == 1)THEN ! Nclass formulation
                     mmm_vpar_ze(:)  = vpar_nclass(:)
                ELSE
                   IF(myid == master) &
                   write(ncrt,FMT='("error, sub set_mmm_vars mmm_vpar_ze not available")')
                   lerrno = iomaxerr+99
                   CALL terminate(lerrno,nlog)  
                     !assume vtor = vpar *btor/B to get vpar ! this is crap
                     mmm_vpar_ze(:) = mmm_vtor_ze(:)*SQRT(bp_outbd(:)**2 +  &
                                                          btor_outbd(:)**2)/btor_outbd(:)

                     !mmm_vpol_ze(1) = zeroc
                     grad_ti    = ( mmm_ti_ze(npoints)-mmm_ti_ze(npoints-1) )/ &
                                       ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) ) 
                     eps        = mmm_rmin_ze(npoints)/mmm_rmaj_ctr_ze(npoints)
                     fc         = 1.00_DP +(0.46_DP*eps -1.46_DP)*SQRT(eps) ! trapped partcle farc
                     K1         = 0.8839_DP*fc/(0.3477_DP+0.4058_DP*fc)
                     Bpsq       = mmm_bpol_ze(npoints)*mmm_bpol_ze(npoints) 
                    ! mmm_vpol_ze(npoints) = alpha_dia*K1*cfctr*grad_ti*Bt/(Bt*Bt + Bpsq) ! m/sec (1,nj)
                ENDIF
          ELSE           ! use vpar from statefile
                         ! in this case these values dont change as the 
                         ! profiles evolve in time
               IF(use_forcebal == 1)THEN
                 !profile%vpol_nclass (from statefile) was put  on current 
                 !grid in refactor_profiles
                 !Note that use_forcebal here means use profile%vpol_nclass%data
                 ! instead of profile%vpol%data but the two may infact be the same
                 mmm_vpar_ze(:) = profile%vpar_nclass%data(:)
               ELSE
                 ! put profile%vpar (from state file) on current grid
                 mmm_vpar_ze(:) = profile%vpar%data(:)

               ENDIF
          ENDIF



          IF(use_mmm_er_input == 0 )THEN

                IF(use_forcebal == 1)THEN ! Nclass formulation
                     mmm_er_ze(:)    = er_tot_nclass(:)   ! total radial electric field
                ELSE
                   IF(myid == master) &
                   write(ncrt,FMT='("error, sub set_mmm_vars  1  mmm_er_ze not available")')
                   lerrno = iomaxerr+99
                   CALL terminate(lerrno,nlog)  

                ENDIF

          ELSE           ! use Er from statefile
                         ! in this case these values dont change as the 
                         ! profiles evolve in time
              IF(myid == master) &
                ! write(ncrt,FMT='("error, sub set_mmm_vars 2  Er not available")')
                ! lerrno = iomaxerr+99
                ! CALL terminate(lerrno,nlog)  
                 ! THIS STOP EXISTS BECAUSE profile%er_nclass%data and profile%er%data need to be introduced

               mmm_er_ze(:)    = er_tot_nclass(:)   ! total radial electric field 
                                                    ! not currently taken from statfile
          ENDIF

          grad_ti = zeroc

    cmpct:     IF(use_compact_schemes == izero )THEN
                    ! FINITE DIFF FORM OF GRADIENTS.  mmm_rmin_ze MAY BE NON UNIFORM GRID
                    ! sO ACCURACY IS O(H) (NOT O(H^2) - SHOULD IMPLEMENT COMPACT SCHEMES HERE) hsj

                    DO jr=1,npoints 

                       IF(jr .GT. 1 .AND. jr .LT. npoints)THEN
                          mmm_gte_ze(jr) = - mmm_rmaj_ctr_ze(jr)/mmm_te_ze(jr)* &
                            (mmm_te_ze(jr+1)- mmm_te_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )
                            grad_ti        = ( mmm_ti_ze(jr+1)-mmm_ti_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) ) 
                          mmm_gti_ze(jr) = - mmm_rmaj_ctr_ze(jr)/ mmm_ti_ze(jr)* grad_ti
             
                        !-------------------test compact deriv scheme (temp setup for ti only)
                        done = .true. ! scheme is off
                        IF(.not. done)THEN
                           done = .true.
                           inpt_dpl = 1  ; dpl = zeroc             ! all profiles have dudr =0 at r =0
                           inpt_dpr = izero                        ! no derivative condition at r = a 
                           nrhs = 1 ! nrhs > 1 not implemented 
                           !  use sixth order uniform grid spacing if possible:
                           zctr_deriv = 0  ! get deriv on zone faces
                           msp = 1         ! note rmin_ze is non uniform grid
                           order =4        ! even if r is uniform , order =2,4
                           CALL compact_deriv(mmm_ti_ze,mmm_rmin_ze,nj,order,zctr_deriv,msp,dpl,inpt_dpl,dpr, &
                           inpt_dpr,nrhs,dpdg)
                            mmm_gti_ze(:) = - mmm_rmaj_ctr_ze(:)/ mmm_ti_ze(:)*dpdg(:)
                           !CALL midpt_interp(drho_drmin_fg(0),nj,nrhs,drho_drmin_zc(0))
                        ENDIF
             !---------------
                          mmm_gne_ze(jr) = - mmm_rmaj_ctr_ze(jr)/mmm_ne_ze(jr)* &
                            ( mmm_ne_ze(jr+1)-mmm_ne_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )
                          grad_den = ( mmm_nh_ze(jr+1) - mmm_nh_ze(jr-1))/mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1)
                          grad_denni =( mmm_ni_ze(jr+1) - mmm_ni_ze(jr-1))/mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1)
                          mmm_gnh_ze(jr) = - mmm_rmaj_ctr_ze(jr)/mmm_nh_ze(jr)* grad_den

                          mmm_gni_ze(jr) = - mmm_rmaj_ctr_ze(jr)/mmm_ni_ze(jr)* grad_denni

                          mmm_gnz_ze(jr) = - mmm_rmaj_ctr_ze(jr)/mmm_nz_ze(jr)* &
                            ( mmm_nz_ze(jr+1)-mmm_nz_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )
                          mmm_gq_ze(jr) =   mmm_rmaj_ctr_ze(jr)/mmm_q_ze(jr)* &
                            ( mmm_q_ze(jr+1)-mmm_q_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )


                          mmm_gvtor_ze(jr) =  mmm_rmaj_ctr_ze(jr)/mmm_vtor_ze(jr)* &
                            ( mmm_vtor_ze(jr+1)-mmm_vtor_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )

                          mmm_gvpol_ze(jr) =  mmm_rmaj_ctr_ze(jr)/mmm_vpol_ze(jr)* &
                            ( mmm_vpol_ze(jr+1)-mmm_vpol_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )

                          mmm_gvpar_ze(jr) =  mmm_rmaj_ctr_ze(jr)/mmm_vpar_ze(jr)* &
                            ( mmm_vpar_ze(jr+1)-mmm_vpar_ze(jr-1) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )

           
                          mmm_gelong_ze(jr)    =  (mmm_elong_ze(jr+1)- mmm_elong_ze(jr-1))/ &
                                                   (mmm_rmin_ze(jr+1)/mmm_rmaj_ctr_ze(jr+1) &
                                                  - mmm_rmin_ze(jr-1)/mmm_rmaj_ctr_ze(jr-1))

                          mmm_gangrot_ze(jr)   = (mmm_rmaj_ctr_ze(jr)/mmm_angrot_ze(jr))*      &
                                                  (mmm_angrot_ze(jr+1)- mmm_angrot_ze(jr-1))/  &
                                                  ( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr-1) )
                         grad_press = (grad_den*mmm_ti_ze(jr) + grad_ti* mmm_nh_ze(jr))*joupkev

                       ENDIF

 

                       Bpsq       = mmm_bpol_ze(jr)*mmm_bpol_ze(jr)                 ! Bp at (R,zma)

                       IF(jr == npoints)THEN
                             grad_ti = (mmm_ti_ze(jr)-mmm_ti_ze(jr-1))/( mmm_rmin_ze(jr)-mmm_rmin_ze(jr-1))
                             grad_den  =  (mmm_nh_ze(jr) - mmm_nh_ze(jr-1))/(mmm_rmin_ze(jr)-mmm_rmin_ze(jr-1))
                             grad_denni  =  (mmm_ni_ze(jr) - mmm_ni_ze(jr-1))/(mmm_rmin_ze(jr)-mmm_rmin_ze(jr-1))
                             grad_press = (grad_den*mmm_ti_ze(jr) + grad_ti* mmm_nh_ze(jr))*joupkev

                             mmm_gte_ze(npoints) = - mmm_rmaj_ctr_ze(npoints)/mmm_te_ze(npoints)* &
                                  ( mmm_te_ze(npoints)-mmm_te_ze(npoints-1) )/&
                                  ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )
                             mmm_gti_ze(npoints) = - mmm_rmaj_ctr_ze(npoints)/mmm_ti_ze(npoints)* &
                                  grad_ti
                             mmm_gne_ze(npoints) = - mmm_rmaj_ctr_ze(npoints)/mmm_ne_ze(npoints)* &
                                  ( mmm_ne_ze(npoints)-mmm_ne_ze(npoints-1) )/&
                                  ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )
                             mmm_gnh_ze(npoints) = - mmm_rmaj_ctr_ze(npoints)/mmm_nh_ze(npoints)* &
                                  grad_den
                             mmm_gni_ze(npoints) = - mmm_rmaj_ctr_ze(npoints)/mmm_ni_ze(npoints)* &
                                  grad_denni

                             mmm_gnz_ze(npoints) = - mmm_rmaj_ctr_ze(npoints)/mmm_nz_ze(npoints)* &
                                  ( mmm_nz_ze(npoints) - mmm_nz_ze(npoints-1) )/&
                                  ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )
                             mmm_gq_ze(npoints) =  mmm_rmaj_ctr_ze(npoints)/mmm_q_ze(npoints)* &
                                  ( mmm_q_ze(npoints) - mmm_q_ze(npoints-1) )/&
                                  ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )

                             mmm_gvtor_ze(npoints) =  mmm_rmaj_ctr_ze(npoints)/mmm_vtor_ze(npoints)* &
                                  ( mmm_vtor_ze(npoints) - mmm_vtor_ze(npoints-1) )/&
                                  ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )

                             mmm_gvpol_ze(npoints) =  mmm_rmaj_ctr_ze(npoints)/mmm_vpol_ze(npoints)* &
                                  ( mmm_vpol_ze(npoints) - mmm_vpol_ze(npoints-1) )/&
                                  ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )
                             mmm_gvpar_ze(npoints) =  mmm_rmaj_ctr_ze(npoints)/mmm_vpar_ze(npoints)* &
                                  ( mmm_vpar_ze(npoints) - mmm_vpar_ze(npoints-1) )/&
                                  ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )

                             mmm_gelong_ze(npoints)    =  (mmm_elong_ze(npoints)- mmm_elong_ze(npoints-1))/ &
                                                   (mmm_rmin_ze(npoints)/mmm_rmaj_ctr_ze(npoints) &
                                                  - mmm_rmin_ze(npoints-1)/mmm_rmaj_ctr_ze(npoints-1))

                             mmm_gangrot_ze(npoints)    =  (mmm_rmaj_ctr_ze(npoints)/mmm_angrot_ze(npoints))* & 
                                                            (mmm_angrot_ze(npoints)- mmm_angrot_ze(npoints-1))/ &
                                                            ( mmm_rmin_ze(npoints)-mmm_rmin_ze(npoints-1) )
                       ENDIF



                        bps = mmm_bpol_ze(jr)
                       ! radial electric field ,use primary ions,z=1, both grad_press and mm_nh_ze  
                       ! are for total primary ions ( Z=1)
                       IF(use_forcebal .NE. 1) &
                                               mmm_er_ze(jr)   = mmm_vtor_ze(jr)*bps       &
                                              - mmm_vpol_ze(jr)*mmm_btor_ze(jr) &
                                              + grad_press/ABS(Electron_Charge *mmm_nh_ze(jr))

                       mmm_vdia_ze(jr) = alpha_dia*grad_press/(-Electron_Charge*mmm_nh_ze(jr)*mmm_btor_ze(jr)) ! m/sec (1,nj)


                       mmm_vexb_ze(jr) = -mmm_theta_ze(jr)*vphi(jr) - mmm_vdia_ze(jr)    &
                                         + 2._DP* mmm_vpol_ze(jr)/hcap(jr)       ! m/sec (1,nj)  ( Vexb =- Er/Bt )
 
                       ! calc d/dpsi(Er/(R Bpol)): We use Bpol(R,zma) on outboard side:
                       ! psi on the nj grid is given by
                         IF(jr .NE. 1) mmm_tor_ang_speed_ze(jr) = mmm_er_ze(jr)/(rmaj_outbd(jr)*bps)


                    END DO
                           
 
                    !Get psi deriv of mmm_tor_ang_speed_ze and finally mmm_wexbs_ze(:)
                    mmm_tor_ang_speed_ze(1) = mmm_tor_ang_speed_ze(2)           ! punt here
                    DO jr = 2,npoints -1
                       dpsi         = mmm_tor_ang_speed_ze(jr+1) - mmm_tor_ang_speed_ze(jr-1)
                       dpsi         = dpsi/(psir_grid%data(jr+1)  - psir_grid%data(jr-1))
                       bps          = mmm_bpol_ze(jr)
                       btot         = SQRT(mmm_btor_ze(jr)**2 + bps**2)
                       mmm_wexbs_ze(jr) = dpsi*rmaj_outbd(jr)*rmaj_outbd(jr)*bps*bps/btot
                    ENDDO
                    mmm_wexbs_ze(1)  = mmm_wexbs_ze(2)*(rmaj_outbd(1)/rmaj_outbd(2))**2
                    dpsi         = mmm_tor_ang_speed_ze(npoints) - mmm_tor_ang_speed_ze(npoints-1)
                    dpsi         = dpsi/(psir_grid%data(npoints) - psir_grid%data(npoints-1))
                    bps          = mmm_bpol_ze(npoints)
                    btot         = SQRT(mmm_btor_ze(npoints)**2 + bps**2)
                    mmm_wexbs_ze(npoints) = dpsi*rmaj_outbd(npoints)*rmaj_outbd(npoints)*bps*bps/btot


                    !Take care of on axis values:
                             grad_ti     = (mmm_ti_ze(2)-mmm_ti_ze(1))/( mmm_rmin_ze(2)-mmm_rmin_ze(1))
                             grad_den    = (mmm_nh_ze(2) - mmm_nh_ze(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))
                             grad_denni  = (mmm_ni_ze(2) - mmm_ni_ze(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))
                             grad_press  = (grad_den*mmm_ti_ze(1) + grad_ti* mmm_nh_ze(1))*joupkev
                             grad_te     = zeroc
                             grad_ti     = zeroc
                             grad_den    = zeroc
                             grad_press  = zeroc
                             grad_angrot = zeroc
                             grad_ne     = zeroc
                             grad_denni  = zeroc
                             mmm_gte_ze(1) = - mmm_rmaj_ctr_ze(1)/mmm_te_ze(1)* &
                                  grad_te
                             mmm_gti_ze(1) = - mmm_rmaj_ctr_ze(1)/mmm_ti_ze(1)* &
                                  grad_ti
                             mmm_gne_ze(1) = - mmm_rmaj_ctr_ze(1)/mmm_ne_ze(1)* &
                                  grad_ne
                             mmm_gnh_ze(1) = - mmm_rmaj_ctr_ze(1)/mmm_nh_ze(1)* &
                                  grad_den
                             mmm_gni_ze(1) = - mmm_rmaj_ctr_ze(1)/mmm_ni_ze(1)* &
                                  grad_denni

                             mmm_gnz_ze(1) = - mmm_rmaj_ctr_ze(1)/mmm_nz_ze(1)* &
                                  ( mmm_nz_ze(2) - mmm_nz_ze(1) )/&
                                  ( mmm_rmin_ze(2)-mmm_rmin_ze(1) )
                             mmm_gq_ze(1) =  mmm_rmaj_ctr_ze(1)/mmm_q_ze(1)* &
                                  ( mmm_q_ze(2) - mmm_q_ze(1) )/&
                                  ( mmm_rmin_ze(2)-mmm_rmin_ze(1) )

                             mmm_gvtor_ze(1) =  mmm_rmaj_ctr_ze(1)/mmm_vtor_ze(1)* &
                                  ( mmm_vtor_ze(2) - mmm_vtor_ze(1) )/&
                                  ( mmm_rmin_ze(2)-mmm_rmin_ze(1) )

                             !mmm_gvpol_ze(1) =  mmm_rmaj_ctr_ze(1)/mmm_vpol_ze(1)* &
                             !     ( mmm_vpol_ze(2) - mmm_vpol_ze(1) )/&
                             !     ( mmm_rmin_ze(2)-mmm_rmin_ze(1) )
                             mmm_gvpol_ze(1) =  mmm_gvpol_ze(2) ! need limit here

 
                             mmm_gvpar_ze(1) =  mmm_rmaj_ctr_ze(1)/mmm_vpar_ze(1)* &
                                  ( mmm_vpar_ze(2) - mmm_vpar_ze(1) )/&
                                  ( mmm_rmin_ze(2)-mmm_rmin_ze(1) )

                             mmm_gelong_ze(1)   =  (mmm_elong_ze(2)- mmm_elong_ze(1))/ &
                                                   (mmm_rmin_ze(2)/mmm_rmaj_ctr_ze(2) &
                                                  - mmm_rmin_ze(1)/mmm_rmaj_ctr_ze(1))

                             mmm_gangrot_ze(1)   =   mmm_rmaj_ctr_ze(1)/mmm_angrot_ze(1)* &
                                  grad_angrot

                 ELSE   cmpct
                    WRITE(ncrt,FMT='("   Error, sub set_mmm_vars,",/, &
                                     "      use _compact schemes not implemented")')
                    lerrno = iomaxerr+187
                    CALL terminate(lerrno,nlog)  



                 ENDIF  cmpct


    RETURN


    WRITE(969,FMT='("mmm_rmaj_ctr_ze(jr),mmm_rmin_ze(jr),r(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_rmaj_ctr_ze(jr),mmm_rmin_ze(jr),r(jr)
    ENDDO

    WRITE(969,FMT='("mmm_elong_ze(jr)(jr),mmm_q_ze(jr),bp_bpol_ze(jr)")')
    DO jr =1,npoints
       WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_elong_ze(jr),mmm_q_ze(jr),mmm_bpol_ze(jr)
    ENDDO

    WRITE(969,FMT='("mmm_theta_ze(jr),mmm_btor_ze(jr),mmm_bpol_ze(jr)")')
    DO jr =1,npoints
       WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_theta_ze(jr),mmm_btor_ze(jr),mmm_bpol_ze(jr)
    ENDDO
    WRITE(969,FMT='("mmm_te_ze(jr),mmm_ti_ze(jr),mmm_zeff_ze(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_te_ze(jr),mmm_ti_ze(jr),mmm_zeff_ze(jr)
    ENDDO

 
    WRITE(969,FMT='("mmm_aimp_ze(jr),mmm_zimp_ze(jr),mmm_nz_ze(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_aimp_ze(jr),mmm_zimp_ze(jr),mmm_nz_ze(jr)
    ENDDO

    WRITE(969,FMT='("mmm_ahyd_ze(jr),mmm_aimass_ze(jr),mmm_nz_ze(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_ahyd_ze(jr),mmm_aimass_ze(jr),mmm_nz_ze(jr)
    ENDDO

     WRITE(969,FMT='("mmm_gte_ze(jr),mmm_gti_ze(jr),mmm_gne_ze(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gte_ze(jr),mmm_gti_ze(jr),mmm_gne_ze(jr)
    ENDDO

    WRITE(969,FMT='("mmm_gnh_ze(jr),mmm_gnz_ze(jr),mmm_gq_ze(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gnh_ze(jr),mmm_gnz_ze(jr),mmm_gq_ze(jr)
    ENDDO
  
    WRITE(969,FMT='("mmm_gni_ze(:),mmm_gnh_ze(:),mmm_gnz_ze(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gni_ze(jr),mmm_gnh_ze(jr),mmm_gnz_ze(jr)
    ENDDO

    WRITE(969,FMT='("mmm_gvtor_ze(:),mmm_vtor_ze(:),mmm_gnz_ze(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gvtor_ze(jr),mmm_vtor_ze(jr),mmm_gnz_ze(jr)
    ENDDO

    WRITE(969,FMT='("mmm_vpol_ze(:),mmm_gvpol_ze(:),mmm_vdia_ze(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_vpol_ze(jr),mmm_gvpol_ze(jr),mmm_vdia_ze(jr)
    ENDDO

    WRITE(969,FMT='("mmm_er_ze(:),mmm_vexb_ze(:),mmm_vpar_ze(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_er_ze(jr),mmm_vexb_ze(jr),mmm_vpar_ze(jr)
    ENDDO

    WRITE(969,FMT='(" mmm_tor_ang_speed_ze(jr),mmm_wexbs_ze(jr),mmm_vpar_ze(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_tor_ang_speed_ze(jr),mmm_wexbs_ze(jr),mmm_vpar_ze(jr)
    ENDDO


!   CALL debug_stop('set_mmm_vars testing line 748 mmm_in.f90')
 

    RETURN


    END SUBROUTINE set_mmm_vars


         SUBROUTINE  set_mmm_vars_zct
         !--------------------------------------------------------------------------
         ! -- setup input variables for mmm
         ! -- NOTE that these are defined on the outboard side in terms of half width 
         ! -- of each flux surface, as defined by dischg%rmin_half_width
         ! -- set variables at zone ceneters,**_zct
         !--------------------------------------------------------------------------
           USE grid_class,                           ONLY : nj,fcap,gcap,hcap,       &
                                                            use_compact_schemes,     &
                                                            psir_grid,r

           USE common_constants,                     ONLY : izero,zeroc,joupkev,     &
                                                            Electron_Charge

           USE error_handler,                        ONLY : iomaxerr, lerrno,        &
                                                            terminate
           USE io_gcnmp,                             ONLY : nlog, ncrt
      
           USE dep_var,                              ONLY : rbp,te,ti,angrot,en,     &
                                                            ene,vphi,vpol_nclass,    &
                                                            vpol,vpar_nclass,        &
                                                            er_tot_nclass

           USE ions_gcnmp,                           ONLY : nion,nprim,zeff,z,       &
                                                            atw
  
           USE neutral_beams,                        ONLY : nbion,enbeam_tot

           USE curden_terms,                         ONLY : q,bp

           USE Plasma_properties,                    ONLY : mhd_dat,dischg,profile

           USE solcon_gcnmp,                         ONLY : use_mmm_vpol_input,      &
                                                            use_mmm_vpar_input,      &
                                                            use_mmm_er_input,      &
                                                            use_forcebal

           USE Vector_class,                         ONLY : get_values

           USE plasma_contours,                      ONLY : bp_outbd,rmaj_outbd,btor_outbd

           USE MPI_data,                             ONLY : master,myid,numprocs

           IMPLICIT NONE
           
           INTEGER(I4B) jr,j
           REAL(DP), ALLOCATABLE,DIMENSION(:) :: s_den_zct,s_den_ze
           REAL(DP) alpha_dia,bt,bpsq,fc,eps, cfctr,grad_press,grad_ti,grad_den, &
                    K1,bps,dpsi,btot,grad_denni,grad_te,grad_angrot,grad_ne,     &
                    rmaj_zct


               IF( .NOT. mmm_initialized)CALL mmm_init

               IF(ALLOCATED(s_den_zct))DEALLOCATE(s_den_zct)
               ALLOCATE(s_den_zct(SIZE(mmm_ne_ze-1)))

               IF(ALLOCATED(s_den_ze))DEALLOCATE(s_den_ze)
               ALLOCATE(s_den_ze(SIZE(mmm_ne_ze-1)))

               Bt                   = mhd_dat%btor                         
               cfctr                = -joupkev/Electron_Charge    ! Electron_Charge is negative
               alpha_dia            = 1._DP             ! multiplier ,assign later if wanted
               vphi(:)              = dischg%rmajor*angrot(:)
               mmm_vtor_ze(:)       = vphi(:)



               npoints = nj

               mmm_rmin_ze(:)       = get_values(dischg%rmin_half_width)
               mmm_rmaj_ctr_ze(:)   = get_values(dischg%rmaj_geom) ! rmajor to center of flux surface



               mmm_nh_zct(:)        = zeroc            ! total primary  ion density
               mmm_nh_ze(:)         = zeroc 
               mmm_nz_ze(:)         = zeroc
               mmm_ni_zct(:)        = zeroc            ! total thermal  ion density

               ! pick up edge values for zone edge grid parameters:
               mmm_ne_ze(nj)        = ene(nj)
               mmm_nf_ze(nj)        = enbeam_tot(nj)
               DO jr = 1, nprim
                     mmm_nh_ze(nj)  = mmm_nh_ze(nj)  + en(nj,jr)
               ENDDO
               mmm_ni_ze(nj)        = mmm_nh_ze(nj)
               DO jr = nprim+1,nion
                  mmm_ni_ze(nj)     = mmm_ni_ze(nj)  + en(nj,jr)
                  mmm_nz_ze(nj)     = mmm_nz_ze(nj)  + en(nj,jr)
               ENDDO
               mmm_q_ze(nj)         = ABS(q(nj))


               DO j=1,nj-1
                  mmm_rmaj_ctr_zct(j)= 0.5_DP*(mmm_rmaj_ctr_ze(j+1)+mmm_rmaj_ctr_ze(j))
                  mmm_te_zct(j)    = 0.5_DP*(te(j+1)+ te(j))
                  mmm_ti_zct(j)    = 0.5_DP*(ti(j+1)+ ti(j))
                  mmm_zeff_zct(j)  = 0.5_DP*(zeff(j+1)+ zeff(j))
                  mmm_q_zct(j)     = 0.5_DP*(q(j+1)+ q(j))
                  mmm_q_zct(j)     = ABS(mmm_q_zct(j))
                  mmm_q_ze(j)      = ABS(q(j))
                  mmm_ne_zct(j)    = 0.5_DP*(ene(j+1)+ ene(j))
                  mmm_ne_ze(j)     = ene(j)
                  mmm_nf_zct(j)    = 0.5_DP*(enbeam_tot(j+1)+ enbeam_tot(j))
                  mmm_nf_ze(j)     = enbeam_tot(j)
                  mmm_btor_zct(j)  = 0.5_DP*(btor_outbd(j+1)+ btor_outbd(j))
                  mmm_bpol_zct(j)  = 0.5_DP*(bp_outbd(j+1)+ bp_outbd(j))
                  mmm_angrot_zct(j)= 0.5_DP*(angrot(j+1)+angrot(j))
                  mmm_vtor_zct(j)  = 0.5_DP*(vphi(j+1)+vphi(j))


                  DO jr = 1, nprim
                     mmm_nh_zct(j) = mmm_nh_zct(j) + 0.5_dp*(  en(j+1,jr) +en(j,jr))
                     mmm_nh_ze(j)  = mmm_nh_ze(j)  + en(j,jr)
                  ENDDO

                  mmm_zimp_zct(j)    = zeroc
                  s_den_zct(j)       = zeroc ; s_den_ze(j) = zeroc
                  mmm_aimp_zct(j)    = zeroc
                  DO jr = nprim+1,nion
                    mmm_zimp_zct(j)  = mmm_zimp_zct(j) + 0.5_dp*(en(j+1,jr)*z(j+1,jr) &
                                                                 +en(j,jr)*z(j,jr))
                    s_den_zct(j)     = s_den_zct(j) + 0.5_DP*(en(j+1,jr)+en(j,jr))
                    s_den_ze(j)      = s_den_ze(j) + en(j,jr)
                    mmm_aimp_zct(j)  = mmm_aimp_zct(j) + 0.5_DP*(en(j+1,jr) +en(j,jr))*atw(jr)
                 ENDDO
                 mmm_nz_zct(j)       = s_den_zct(j)
                 mmm_nz_ze(j)        = s_den_ze(j)
                 mmm_ni_zct(j)       = mmm_nh_zct(j) + mmm_nz_zct(j)
                 mmm_ni_ze(j)        = mmm_nh_ze(j)  + mmm_nz_ze(j)
                 mmm_zimp_zct(j)     = mmm_zimp_zct(j)/s_den_zct(j)
                 mmm_aimp_zct(j)     = mmm_aimp_zct(j)/s_den_zct(j)


                 ! get mean mass of  hydrogenic,thermal and impurity ions:
                 mmm_ahyd_zct(j)     = zeroc
                 s_den_zct(j)        = zeroc
                 DO jr = 1,nprim
                    mmm_ahyd_zct(j)  = mmm_ahyd_zct(j) + 0.5+DP*(en(j+1,jr)+ en(j,jr))*atw(jr)
                    s_den_zct(j)     = s_den_zct(j) + 0.5_DP*(en(j+1,jr)+en(j,jr))
                 ENDDO
                 mmm_aimass_zct(j)   = mmm_ahyd_zct(j)           ! hydrogenic  part of total
                 mmm_ahyd_zct(j)     = mmm_ahyd_zct(j)/s_den_zct(j)  ! mean hydrogenic  part
                
                 DO jr =nprim+1,nion    ! add  impurity part
                    mmm_aimass_zct(j)= mmm_aimass_zct(j) + 0.5_DP*(en(j+1,jr)+en(j,jr))*atw(jr)
                    s_den_zct(j)         = s_den_zct(j) + 0.5_DP*(en(j+1,jr)+en(j,jr))
                 ENDDO
                 mmm_aimass_zct(j)   = mmm_aimass_zct(j)/s_den_zct(j) ! mean mass of thermal ions

              ENDDO




               mmm_gte_zct(1)   = -(mmm_rmaj_ctr_zct(1)/mmm_te_zct(1))* &
                                   (te(2) - te(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))
               mmm_gti_zct(1)   = -(mmm_rmaj_ctr_zct(1)/mmm_ti_zct(1))* &
                                   (ti(2) - ti(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))
               mmm_gne_zct(1)   = -(mmm_rmaj_ctr_zct(1)/mmm_ne_zct(1))* &
                                   (ene(2) - ene(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))


               mmm_gnh_zct(1)   = -(mmm_rmaj_ctr_zct(1)/mmm_nh_zct(1))* &
                                   (mmm_nh_ze(2)-mmm_nh_ze(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))
               mmm_gni_zct(1)   = -(mmm_rmaj_ctr_zct(1)/mmm_ni_zct(1))* &
                                   (mmm_ni_ze(2)-mmm_ni_ze(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))
               mmm_gnz_zct(1)   = -(mmm_rmaj_ctr_zct(1)/mmm_nz_zct(1))* &
                                   (mmm_nz_ze(2)-mmm_nz_ze(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))

               !  R ( d q   / d r ) / q
               mmm_gq_zct(1)    =  (mmm_rmaj_ctr_zct(1)/mmm_q_zct(1))* &
                                    (mmm_q_ze(2)-mmm_q_ze(1))/(mmm_rmin_ze(2)-mmm_rmin_ze(1))
               grad_press         = zeroc ; grad_ti = zeroc ; grad_den = zeroc



          ! need to set here because we look forward (j+1) in Do loop"
          IF(use_mmm_vpol_input == 0 )THEN

                IF(use_forcebal == 1)THEN ! Nclass formulation
                   mmm_vpol_ze(:)      = vpol_nclass(:)
                   DO j=1,nj-1
                     mmm_vpol_zct(j)   = 0.5_DP*(vpol_nclass(j+1)+vpol_nclass(j))
                   ENDDO
                ELSE
                   DO jr=1,npoints-1 ! Y.B.KIM formualtion 
                       grad_ti    = ( ti(jr+1)-ti(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) ) 
                       eps        = mmm_rmin_zct(jr)/mmm_rmaj_ctr_zct(jr)
                       fc         = 1.00_DP +(0.46_DP*eps -1.46_DP)*SQRT(eps) ! trapped partcle farc
                       K1         = 0.8839_DP*fc/(0.3477_DP+0.4058_DP*fc)
                       Bpsq       = mmm_bpol_zct(jr)*mmm_bpol_zct(jr) 
                       mmm_vpol_zct(jr) = alpha_dia*K1*cfctr*grad_ti*Bt/(Bt*Bt + Bpsq) ! m/sec (1,nj)
                   ENDDO
 
                ENDIF

          ELSE           ! use vpol or vpol_nclass,vpar from statefile
                         ! in this case these values dont change as the profiles evolve in time
               DO j=1,npoints-1
                  IF(use_forcebal == 1)THEN
                     !profile%vpol_nclass (from statefile) was put  on current grid in refactor_profiles

                     mmm_vpol_zct(j) = 0.5_DP*(profile%vpol_nclass%data(j+1) + profile%vpol_nclass%data(j+1))
                  ELSE
                     ! put profile%vpol (from state file) on current grid
                     mmm_vpol_zct(j) = 0.5_DP*(profile%vpol%data(j+1)+ profile%vpol%data(j))
                  ENDIF
               ENDDO
          ENDIF


          IF(use_mmm_vpar_input == 0 )THEN

                IF(use_forcebal == 1)THEN     ! Nclass formulation
                   mmm_vpar_ze(:)      = vpar_nclass(:)
                   DO j=1,nj-1
                     mmm_vpar_zct(j)   = 0.5_DP*(vpar_nclass(j+1)+vpar_nclass(j))
                   ENDDO
                ELSE

                   IF(myid == master) &
                   write(ncrt,FMT='("error, sub set_mmm_vars_zct vpar not available")')
                   lerrno = iomaxerr+99
                   CALL terminate(lerrno,nlog)  

                   DO jr=1,npoints-1 ! Y.B.KIM formualtion 
                       grad_ti    = ( ti(jr+1)-ti(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) ) 
                       eps        = mmm_rmin_zct(jr)/mmm_rmaj_ctr_zct(jr)
                       fc         = 1.00_DP +(0.46_DP*eps -1.46_DP)*SQRT(eps) ! trapped partcle farc
                       K1         = 0.8839_DP*fc/(0.3477_DP+0.4058_DP*fc)
                       Bpsq       = mmm_bpol_zct(jr)*mmm_bpol_zct(jr) 
                      ! mmm_vpol_zct(jr) = alpha_dia*K1*cfctr*grad_ti*Bt/(Bt*Bt + Bpsq) ! m/sec (1,nj)
                       !mmm_vpar_zct(jr) = mmm_vtor_zct(jr)*0.5_DP*(SQRT(bp_outbd(jr+1)**2 +        & ! this is crap!
                       !                                   btor_outbd(jr+1)**2)/btor_outbd(jr+1) +  &
                       !                                   SQRT(bp_outbd(jr)**2 +        &
                       !                                   btor_outbd(jr)**2)/btor_outbd(jr))
                   ENDDO
 
                ENDIF

          ELSE           ! use vpar or vpar_nclass,from statefile
                         ! in this case these values dont change as the profiles evolve in time
               DO j=1,npoints-1
                  IF(use_forcebal == 1)THEN
                     !profile%vpar_nclass (from statefile) was put  on current grid in refactor_profiles
                     mmm_vpar_zct(j) = 0.5_DP*(profile%vpar_nclass%data(j+1) + profile%vpar_nclass%data(j))
                  ELSE
                     mmm_vpar_zct(j) = 0.5_DP*(profile%vpar%data(j+1)+ profile%vpar%data(j))
                  ENDIF
               ENDDO
          ENDIF



          IF(use_mmm_er_input == 0 )THEN

                IF(use_forcebal == 1)THEN ! Nclass formulation
                   mmm_er_ze(:)        = er_tot_nclass(:)   ! total radial electric field
                   DO j=1,nj-1
                     mmm_er_zct(j)     = 0.5_DP*(er_tot_nclass(j+1)+er_tot_nclass(j))   ! total radial electric field
                   ENDDO
                ELSE 
                   IF(myid == master) &
                   write(ncrt,FMT='("error, sub set_mmm_vars_zct mmm_er_ze not available")')
                   lerrno = iomaxerr+99
                   CALL terminate(lerrno,nlog)  
                ENDIF

          ELSE           ! use er ,from statefile
                         ! in this case these values dont change as the profiles evolve in time

              IF(myid == master) &
                 write(ncrt,FMT='("error, sub set_mmm_vars_zct Er not available")')
                 lerrno = iomaxerr+99
                 CALL terminate(lerrno,nlog)  
                 ! THIS STOP EXISTS BECAUSE profile%er_nclass%data and profile%er%data need to be introduced
               DO j=1,npoints-1
                  IF(use_forcebal == 1)THEN
                     ! need profile%er_nclass%data here, not yet defined
                     mmm_er_zct(j)     = 0.5_DP*(er_tot_nclass(j+1)+er_tot_nclass(j))   ! total radial electric field
                  ELSE 
                      ! need profile%er%data here, not yet defined
                     mmm_er_zct(j)     = 0.5_DP*(er_tot_nclass(j+1)+er_tot_nclass(j))   ! total radial electric field
                  ENDIF
               ENDDO
          ENDIF





          grad_ti = zeroc

    cmpct:     IF(use_compact_schemes == izero )THEN
                    ! FINITE DIFF FORM OF GRADIENTS.  mmm_rmin_zct MAY BE NON UNIFORM GRID
                    ! sO ACCURACY IS O(H) (NOT O(H^2) - SHOULD IMPLEMENT COMPACT SCHEMES HERE) hsj

                    DO jr=1,npoints-1 

                       IF(jr .GT. 1 .AND. jr .LT. npoints)THEN
                          mmm_gte_zct(jr) = - mmm_rmaj_ctr_zct(jr)/mmm_te_zct(jr)* &
                            (te(jr+1)- te(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )
                            grad_ti        = (ti(jr+1)-ti(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) ) 
                          mmm_gti_zct(jr) = - mmm_rmaj_ctr_zct(jr)/ mmm_ti_zct(jr)* grad_ti
             
 
                          mmm_gne_zct(jr) = - mmm_rmaj_ctr_zct(jr)/mmm_ne_zct(jr)* &
                            ( mmm_ne_ze(jr+1)-mmm_ne_ze(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )
                          grad_den = ( mmm_nh_ze(jr+1) - mmm_nh_ze(jr))/mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr)
                          grad_denni =( mmm_ni_ze(jr+1) - mmm_ni_ze(jr))/mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr)
                          mmm_gnh_zct(jr) = - mmm_rmaj_ctr_zct(jr)/mmm_nh_zct(jr)* grad_den

                          mmm_gni_zct(jr) = - mmm_rmaj_ctr_zct(jr)/mmm_ni_zct(jr)* grad_denni

                          mmm_gnz_zct(jr) = - mmm_rmaj_ctr_zct(jr)/mmm_nz_zct(jr)* &
                            ( mmm_nz_ze(jr+1)-mmm_nz_ze(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )
                          mmm_gq_zct(jr) =   mmm_rmaj_ctr_zct(jr)/mmm_q_zct(jr)* &
                            ( mmm_q_ze(jr+1)-mmm_q_ze(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )


                          mmm_gvtor_zct(jr) =  mmm_rmaj_ctr_zct(jr)/mmm_vtor_zct(jr)* &
                            ( mmm_vtor_ze(jr+1)-mmm_vtor_ze(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )

                          mmm_gvpol_zct(jr) =  mmm_rmaj_ctr_zct(jr)/mmm_vpol_zct(jr)* &
                            ( mmm_vpol_ze(jr+1)-mmm_vpol_ze(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )

                          mmm_gvpar_zct(jr) =  mmm_rmaj_ctr_zct(jr)/mmm_vpar_zct(jr)* &
                            ( mmm_vpar_ze(jr+1)-mmm_vpar_ze(jr) )/( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )

           
                          mmm_gelong_zct(jr)    =  (mmm_elong_ze(jr+1)- mmm_elong_ze(jr))/ &
                                                   (mmm_rmin_ze(jr+1)/mmm_rmaj_ctr_ze(jr+1) &
                                                  - mmm_rmin_ze(jr)/mmm_rmaj_ctr_ze(jr))

                          mmm_gangrot_zct(jr)   = (mmm_rmaj_ctr_zct(jr)/mmm_angrot_zct(jr))*      &
                                                  (angrot(jr+1)- angrot(jr))/  &
                                                  ( mmm_rmin_ze(jr+1)-mmm_rmin_ze(jr) )
                         grad_press = (grad_den*mmm_ti_zct(jr) + grad_ti* mmm_nh_zct(jr))*joupkev

                       ENDIF
 



                       Bpsq       = mmm_bpol_zct(jr)*mmm_bpol_zct(jr)                 ! Bp at (R,zma)

                       

                        bps = mmm_bpol_zct(jr)
                       ! radial electric field ,use primary ions,z=1, both grad_press and mm_nh_zct  
                       ! are for total primary ions ( Z=1)
                       IF(use_forcebal .NE. 1) &
                                               mmm_er_zct(jr)   = mmm_vtor_zct(jr)*bps       &
                                              - mmm_vpol_zct(jr)*mmm_btor_zct(jr) &
                                              + grad_press/ABS(Electron_Charge *mmm_nh_zct(jr))

                       mmm_vdia_zct(jr) = alpha_dia*grad_press/(-Electron_Charge*mmm_nh_zct(jr)*mmm_btor_zct(jr)) ! m/sec (1,nj)


                       mmm_vexb_zct(jr) = -mmm_theta_zct(jr)*vphi(jr) - mmm_vdia_zct(jr)    &
                                         + 2._DP* mmm_vpol_zct(jr)/hcap(jr)       ! m/sec (1,nj)  ( Vexb =- Er/Bt )
 
                       ! calc d/dpsi(Er/(R Bpol)): We use Bpol(R,zma) on outboard side:
                       ! psi on the nj grid is given by
                         mmm_tor_ang_speed_zct(jr) = mmm_er_zct(jr)/ &
                              (0.5_DP*(rmaj_outbd(jr+1) + rmaj_outbd(jr)) *bps)
                    END DO


                    !Get psi deriv of mmm_tor_ang_speed_zct and finally mmm_wexbs_zct(:)
                    mmm_tor_ang_speed_zct(1) = mmm_tor_ang_speed_zct(2)           ! punt here
                    DO jr = 1,npoints -1
                       dpsi         = mmm_tor_ang_speed_ze(jr+1) - mmm_tor_ang_speed_ze(jr)
                       dpsi         = dpsi/(psir_grid%data(jr+1)  - psir_grid%data(jr))
                       bps          = mmm_bpol_zct(jr)
                       btot         = SQRT(mmm_btor_zct(jr)**2 + bps**2)
                       rmaj_zct     = (rmaj_outbd(jr) +  rmaj_outbd(jr+1))/2._DP
                       mmm_wexbs_zct(jr) = dpsi*rmaj_zct*rmaj_zct*bps*bps/btot
                    ENDDO
 

                 ELSE   cmpct
                    WRITE(ncrt,FMT='("   Error, sub set_mmm_vars_zct,",/, &
                                     "      use _compact schemes not implemented")')
                    lerrno = iomaxerr+187
                    CALL terminate(lerrno,nlog)  


 
                 ENDIF  cmpct


    RETURN


    WRITE(969,FMT='("mmm_rmaj_ctr_zct(jr),mmm_rmin_zct(jr),r(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_rmaj_ctr_zct(jr),mmm_rmin_zct(jr),r(jr)
    ENDDO

    WRITE(969,FMT='("mmm_elong_zct(jr)(jr),mmm_q_zct(jr),bp_bpol_zct(jr)")')
    DO jr =1,npoints
       WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_elong_zct(jr),mmm_q_zct(jr),mmm_bpol_zct(jr)
    ENDDO

    WRITE(969,FMT='("mmm_theta_zct(jr),mmm_btor_zct(jr),mmm_bpol_zct(jr)")')
    DO jr =1,npoints
       WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_theta_zct(jr),mmm_btor_zct(jr),mmm_bpol_zct(jr)
    ENDDO
    WRITE(969,FMT='("mmm_te_zct(jr),mmm_ti_zct(jr),mmm_zeff_zct(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_te_zct(jr),mmm_ti_zct(jr),mmm_zeff_zct(jr)
    ENDDO

 
    WRITE(969,FMT='("mmm_aimp_zct(jr),mmm_zimp_zct(jr),mmm_nz_zct(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_aimp_zct(jr),mmm_zimp_zct(jr),mmm_nz_zct(jr)
    ENDDO

    WRITE(969,FMT='("mmm_ahyd_zct(jr),mmm_aimass_zct(jr),mmm_nz_zct(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_ahyd_zct(jr),mmm_aimass_zct(jr),mmm_nz_zct(jr)
    ENDDO

     WRITE(969,FMT='("mmm_gte_zct(jr),mmm_gti_zct(jr),mmm_gne_zct(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gte_zct(jr),mmm_gti_zct(jr),mmm_gne_zct(jr)
    ENDDO

    WRITE(969,FMT='("mmm_gnh_zct(jr),mmm_gnz_zct(jr),mmm_gq_zct(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gnh_zct(jr),mmm_gnz_zct(jr),mmm_gq_zct(jr)
    ENDDO
  
    WRITE(969,FMT='("mmm_gni_zct(:),mmm_gnh_zct(:),mmm_gnz_zct(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gni_zct(jr),mmm_gnh_zct(jr),mmm_gnz_zct(jr)
    ENDDO

    WRITE(969,FMT='("mmm_gvtor_zct(:),mmm_vtor_zct(:),mmm_gnz_zct(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_gvtor_zct(jr),mmm_vtor_zct(jr),mmm_gnz_zct(jr)
    ENDDO

    WRITE(969,FMT='("mmm_vpol_zct(:),mmm_gvpol_zct(:),mmm_vdia_zct(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_vpol_zct(jr),mmm_gvpol_zct(jr),mmm_vdia_zct(jr)
    ENDDO

    WRITE(969,FMT='("mmm_er_zct(:),mmm_vexb_zct(:),mmm_vpar_zct(:)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_er_zct(jr),mmm_vexb_zct(jr),mmm_vpar_zct(jr)
    ENDDO

    WRITE(969,FMT='(" mmm_tor_ang_speed_zct(jr),mmm_wexbs_zct(jr),mmm_vpar_zct(jr)")')
    DO jr =1,npoints
      WRITE(969,FMT='(i5,2x,3(1pe12.4,2x))')jr, mmm_tor_ang_speed_zct(jr),mmm_wexbs_zct(jr),mmm_vpar_zct(jr)
    ENDDO


!   CALL debug_stop('set_mmm_vars_zct testing line 748 mmm_in.f90')
 

    RETURN


    END SUBROUTINE set_mmm_vars_zct

 


    SUBROUTINE mmm_single_grid_point_zct(jm,pertrb)    
      !--------------------------------------------------------------------------
      ! single grid point interface to mmm
      ! jm is zone edge point,convert to zone center here
      !--------------------------------------------------------------------------

      USE nrtype,                                              ONLY : DP,I4B,I2B

      USE grid_class,                                          ONLY : nj,r2capi,r

      USE modmmm7_1,                                           ONLY : mmm7_1

      USE plasma_properties,                                   ONLY : dischg

      USE dep_var,                                             ONLY : ene,en,angrot,ti,te

      USE io_gcnmp,                                            ONLY : nlog, ncrt

      USE ions_gcnmp,                                          ONLY : nprim,nion,atw,nimp

      USE error_handler,                                       ONLY : iomaxerr, lerrno,terminate 

      USE common_constants,                                    ONLY : izero,zeroc,Proton_Mass,   &
           joupkev

      USE modmmm7_1,                                           ONLY : mmm7_1

      IMPLICIT NONE

      INTEGER(I4B) i,j,k

      INTEGER(I2B) nprout

      REAL(DP) ena,mass,sqrteta,scfctr,drmin,denhgrdm,denhgrdp,dtegrdm_sv,       &
           dtigrdm_sv,dwgrdm,dwgrdp,dvtgrdm,dvtgrdp,d2grd,dengradp,              &
           stepsize,dtidrho,dtidrho_new,dtidrmin,dtedrho,dtedrho_new,dtedrmin

      REAL(DP),DIMENSION(1) :: rminlcl,rmajlcl,elonglcl,nelcl,nhlcl,nzlcl,nflcl, &
           zefflcl,telcl,tilcl,qlcl,btorlcl,zimplcl,aimplcl, &
           ahydlcl,aimasslcl,wexbslcl,gnelcl,gnilcl,gnhlcl,  &
           gnzlcl,gtelcl,gtilcl,gqlcl,vtorlcl,vparlcl,       &
           vpollcl,gvtorlcl,gelonglcl,gvpollcl,gvparlcl,     &
           xtilcl,xdilcl,xdzlcl,xtelcl,xvtlcl,xvplcl,        &
           xtiW20lcl, xdiW20lcl,xteW20lcl, xtiDBMlcl,        &
           xdiDBMlcl,xteETGlcl,xteDBMlc,gammaDBMlcl,         &
           omegaDBMlcl,xteDBMlcl

      REAL(DP),DIMENSION(4,1) ::   gammaW20lcl,omegaW20lcl,vfluxlcl

      REAL(DP),DIMENSION(6,1) ::   vconvlcl 

      INTEGER,INTENT(IN) :: jm,pertrb 
      INTEGER(I4B) lcl_size,ok
      LOGICAL testing_mmm

      !testing_mmm = .true. 
      testing_mmm =.false.
      IF(jm .lt. 1 .OR. jm .GE. nj)THEN
         lerrno = iomaxerr + 189
         CALL terminate(lerrno,nlog)
      ENDIF


      nprout = nlog 
      lprint = 1  ! > =6 much output
      sqrteta = 1.e-2_DP*mmm_scfctr ! for 1.e-n perturbation  will be in digit no n
      !-----------------------------------------------------------------------------
      ! -- perturbation values:
      !-----------------------------------------------------------------------------

      slsv: IF(select_solver == 'nwt_pred_cor')THEN
         drmin = mmm_rmin_ze(jm+1)- mmm_rmin_ze(jm)
         scfctr = mmm_scfctr*drmin/mmm_rmaj_ctr_zct(jm)
         IF(pertrb == izero )THEN
            lcl_size = 0 
            IF(ALLOCATED(mmm_val_base_zct)) lcl_size = SIZE(mmm_val_base_zct,2)
            IF(lcl_size .NE. nj-1)THEN
               DEALLOCATE(mmm_val_base_zct,stat = ok)   ! non fatal dealoc
               ALLOCATE(mmm_val_base_zct(-itran_max:itran_max,nj-1))
            ENDIF
            lcl_size = 0 
            IF(ALLOCATED(mmm_val_pert_zct)) lcl_size = SIZE(mmm_val_pert_zct,2)
            IF(lcl_size .NE. nj-1)THEN
               DEALLOCATE(mmm_val_pert_zct,stat = ok)   ! non fatal dealoc
               ALLOCATE(mmm_val_pert_zct(-itran_max:itran_max,nj-1))
            ENDIF

            !note that val_pert is size of perturbation (not perturbed value)
            mmm_val_pert_zct(0,jm)  = zeroc  ! index 0 ==> base values,val_pert =0

            !mmm_val_base_zct(-5,jm) = (mmm_angrot_zct(jm)/mmm_rmaj_ctr_zct(jm))*mmm_gangrot_zct(jm)
            mmm_val_base_zct(-5,jm)  = (mmm_angrot_ze(jm+1)-mmm_angrot_ze(jm))/(r(jm+1)-r(jm))
            mmm_typ_val(-5) = ABS((mmm_angrot_zct(nj-1)-mmm_angrot_zct(1))/(mmm_rmin_zct(nj-1)-mmm_rmin_zct(1)))

            mmm_val_base_zct(-4,jm)  =  zeroc  
            mmm_typ_val(-4)          =  zeroc

            !mmm_val_base_zct(-3,jm) =  (mmm_ti_zct(jm)/mmm_rmaj_ctr_zct(jm))*mmm_gti_zct(jm)
            mmm_val_base_zct(-3,jm)  =  (mmm_ti_ze(jm+1)-mmm_ti_ze(jm))/(r(jm+1)-r(jm))
!           mmm_typ_val(-3)          =  ABS((mmm_ti_zct(nj-1)- mmm_ti_zct(1))/(mmm_rmin_zct(nj-1)-mmm_rmin_zct(1)))
            mmm_typ_val(-3)          =  ABS((mmm_ti_zct(nj-1)- mmm_ti_zct(1))/(r(nj)-r(1)))

            !mmm_val_base_zct(-2,jm) =  (mmm_te_zct(jm)/mmm_rmaj_ctr_zct(jm))*mmm_gte_zct(jm)
            mmm_val_base_zct(-2,jm)  =  (mmm_te_ze(jm+1)-mmm_te_ze(jm))/(r(jm+1)-r(jm))
            mmm_typ_val(-2)          =  ABS((mmm_te_zct(nj-1)- mmm_te_zct(1))/(r(nj)-r(1)))

            !mmm_val_base_zct(-1,jm) = (mmm_nh_zct(jm)/mmm_rmaj_ctr_zct(jm))*mmm_gnh_zct(jm)                mmm_val_base_zct(-1,jm)  = (mmm_nh_ze(jm+1)-mmm_nh_ze(jm))/(r(jm+1)-r(jm))
            mmm_typ_val(-1)          = ABS((mmm_nh_zct(nj-1)- mmm_nh_zct(1))/(mmm_rmin_zct(nj-1)-mmm_rmin_zct(1)))



            mmm_val_base_zct( 0,jm)  =  zeroc               ;  mmm_typ_val(0)  = zeroc
            mmm_val_base_zct( 1,jm)  =  mmm_nh_zct(jm)      ;  mmm_typ_val(1)  = mmm_nh_zct(nj/2)
            mmm_val_base_zct( 2,jm)  =  mmm_te_zct(jm)      ;  mmm_typ_val(2)  = mmm_te_zct(nj/2)
            mmm_val_base_zct( 3,jm)  =  mmm_ti_zct(jm)      ;  mmm_typ_val(3)  = mmm_ti_zct(nj/2)
            mmm_val_base_zct( 4,jm)  =  zeroc               ;  mmm_typ_val(4)  = zeroc
            mmm_val_base_zct( 5,jm)  =  mmm_angrot_zct(jm)  ;  mmm_typ_val(5)  = ABS(mmm_angrot_zct(nj/2))

            ! note step size on gradients is derived from abs(pertrb):
            !stepsize                = sqrteta*MAX(ABS(mmm_val_base(ABS(pertrb),jm)),mmm_typ_val(ABS(pertrb))) &
            !                                       *SIGN(1.D0,mmm_val_base(ABS(pertrb),jm))

            !-----------------------------------------------------------------------------
            ! -- perturb one of the simulation mode variables
            ! -- the magnitude of the pertubation is done based on the current
            ! -- spatial rate of change of the dependent variable.(du/dr)
            ! -- note that pertrb .ne. 0 is active only if dependent variable 
            ! -- is in simulation mode.
            !-----------------------------------------------------------------------------
         ELSE IF (pertrb == 1  )THEN
            mmm_val_pert_zct(pertrb,jm)= sqrteta*MAX(ABS(mmm_val_base_zct(pertrb,jm)),&
                    mmm_typ_val(pertrb))*SIGN(1.D0,mmm_val_base_zct(pertrb,jm))
            !perturb variable given to mmm7_1, restore to original value below:
            mmm_nh_zct(jm) = mmm_val_base_zct(pertrb,jm) + mmm_val_pert_zct(pertrb,jm)


         ELSE IF (pertrb == 2  )THEN
            mmm_val_pert_zct(pertrb,jm)   = sqrteta*MAX(ABS(mmm_val_base_zct(pertrb,jm)),mmm_typ_val(pertrb)) &
                 *SIGN(1.D0,mmm_val_base_zct(pertrb,jm))
            !perturb variable given to mmm7_1, restore to original value below:
            mmm_te_zct(jm) = mmm_val_base_zct(pertrb,jm) + mmm_val_pert_zct(pertrb,jm)



         ELSE IF (pertrb == 3  )THEN
            mmm_val_pert_zct(pertrb,jm) = sqrteta*MAX(ABS(mmm_val_base_zct(pertrb,jm)),mmm_typ_val(pertrb)) &
                 *SIGN(1.D0,mmm_val_base_zct(pertrb,jm))
            !perturb variable given to mmm7_1, restore to original value below:
            mmm_ti_zct(jm) = mmm_val_base_zct( pertrb,jm) + mmm_val_pert_zct(pertrb,jm)



         ELSE IF (pertrb == 4  )THEN 


         ELSE IF (pertrb == 5  )THEN ! perturb omega (toroidal rot speed)
               mmm_val_pert_zct(pertrb,jm)= sqrteta*MAX(ABS(mmm_val_base_zct(pertrb,jm)),&
                    mmm_typ_val(pertrb))*SIGN(1.D0,mmm_val_base_zct(pertrb,jm))

            !perturb variable given to mmm7_1, restore to original value below:
            ! vtor = R*omega:
            mmm_vtor_zct(jm) =  dischg%rmajor*(mmm_val_base_zct(pertrb,jm) + mmm_val_pert_zct(pertrb,jm))
            !mmm_vexb_zct(jm) = -mmm_theta_zct(jm)*mmm_vtor_zct(jm) - mmm_vdia_zct(jm)    &
            !                           + 2._DP* mmm_vpol_zct(jm)/hcap(jm)       ! m/sec (1,nj)  ( Vexb =- Er/Bt )
            !mmm_wexbs_zct(jm) = ??
            ! vpar not changed

         ELSE IF (pertrb == -1 )THEN  
            denhgrdm                   =  mmm_gnh_zct(jm)*mmm_nh_zct(jm)/mmm_rmaj_ctr_zct(jm) ! base grad value at jm
            denhgrdp                   =  mmm_gnh_zct(jm+1)*mmm_nh_zct(jm+1)/mmm_rmaj_ctr_zct(jm+1)
            d2grd                      = (denhgrdp-denhgrdm)/drmin
            mmm_val_pert_zct(pertrb,jm) = d2grd*mmm_scfctr
            dengradp                   = denhgrdm + mmm_val_pert_zct(pertrb,jm)
            mmm_gnh_zct(jm)             = -dengradp*mmm_rmaj_ctr_zct(jm)/mmm_nh_zct(jm)



         ELSE IF (pertrb == -2 )THEN 
            mmm_val_pert_zct(pertrb,jm)   = sqrteta*MAX(ABS(mmm_val_base_zct(pertrb,jm)),mmm_typ_val(pertrb)) &
                 *SIGN(1.D0,mmm_val_base_zct(pertrb,jm))
              dtedrho_new =  mmm_val_base_zct(pertrb,jm) + mmm_val_pert_zct(pertrb,jm) ! this is perturbed dte/drho
              dtedrmin = dtedrho_new*mmm_drhodrmin_zct(jm)
              dtegrdm_sv = mmm_gte_zct(jm) ! enables restore below 
              ! this is what  mmm_7.1 wants as input:
               mmm_gte_zct(jm)              = -(mmm_rmaj_ctr_zct(jm)/mmm_te_zct(jm))*dtedrmin



         ELSE IF (pertrb == -3 )THEN
            mmm_val_pert_zct(pertrb,jm)   = sqrteta*MAX(ABS(mmm_val_base_zct(pertrb,jm)),mmm_typ_val(pertrb)) &
                 *SIGN(1.D0,mmm_val_base_zct(pertrb,jm)) 
              dtidrho_new =  mmm_val_base_zct(pertrb,jm) + mmm_val_pert_zct(pertrb,jm) ! this is perturbed dti/drho
              dtidrmin = dtidrho_new*mmm_drhodrmin_zct(jm)

               dtigrdm_sv = mmm_gti_zct(jm) ! enables restore below 
              ! this is what  mmm_7.1 wants as input:            
               mmm_gti_zct(jm)              = -(mmm_rmaj_ctr_zct(jm)/mmm_ti_zct(jm))*dtidrmin
 

         ELSE IF (pertrb == -4 )THEN

         ELSE IF (pertrb == -5 )THEN
            dwgrdm =  mmm_gangrot_zct(jm)*mmm_angrot_zct(jm)/mmm_rmaj_ctr_zct(jm) ! base grad value at jm
            dwgrdp =  mmm_gangrot_zct(jm+1)*mmm_angrot_zct(jm+1)/mmm_rmaj_ctr_zct(jm+1)
            d2grd = (dwgrdp-dwgrdm)/drmin
            mmm_val_pert_zct(pertrb,jm) = d2grd*mmm_scfctr
            dwgrdp                     = dwgrdm + mmm_val_pert_zct(pertrb,jm)
            mmm_gangrot_zct(jm)         = dwgrdp*mmm_rmaj_ctr_zct(jm)/mmm_angrot_zct(jm)

            dvtgrdm =  mmm_gvtor_zct(jm)*mmm_vtor_zct(jm)/mmm_rmaj_ctr_zct(jm) ! base grad value at jm
            dvtgrdp =  mmm_gvtor_zct(jm+1)*mmm_vtor_zct(jm+1)/mmm_rmaj_ctr_zct(jm+1)
            d2grd = (dvtgrdp-dvtgrdm)/drmin
            dvtgrdp                    = dvtgrdm + d2grd*mmm_scfctr
            mmm_gvtor_zct(jm)           = dvtgrdp*mmm_rmaj_ctr_zct(jm)/mmm_angrot_zct(jm)

         ELSE IF (pertrb .NE. 0)THEN
            lerrno = iomaxerr + 149
            CALL terminate(lerrno,nlog)
         ENDIF

      ENDIF  slsv


      !---------------------------------------------------------------------------
      ! mmm7_1  this call is setup to do a single grid point only:
      ! because of argument associatin we have to go through this crap:  HSJ
      !---------------------------------------------------------------------------
      IF(pertrb == izero)mmm_unstable = .FALSE.
      npoints = 1
      rminlcl(:)       = mmm_rmin_zct(jm)        ;  rmajlcl(:)    = mmm_rmaj_ctr_zct(jm)
      elonglcl(:)      = mmm_elong_zct(jm)       ;  nelcl(:)      = mmm_ne_zct(jm)  
      nhlcl(:)         = mmm_nh_zct(jm)
      nzlcl(:)         = mmm_nz_zct(jm)          ;  nflcl(:)      = mmm_nf_zct(jm)
      zefflcl(:)       = mmm_zeff_zct(jm)        ;  telcl(:)      = mmm_te_zct(jm)  
      qlcl(:)          = mmm_q_zct(jm)           ;  btorlcl(:)    = mmm_btor_zct(jm)     
      zimplcl(:)       = mmm_zimp_zct(jm)        ;  aimplcl       = mmm_aimp_zct(jm)  
      aimasslcl(:)     = mmm_aimass_zct(jm)      ;  wexbslcl(:)   = mmm_wexbs_zct(jm)
      gnelcl(:)        = mmm_gne_zct(jm)         ;  gnilcl(:)     = mmm_gni_zct(jm) 
      gnzlcl(:)        = mmm_gnz_zct(jm)         ;  gtelcl(:)     = mmm_gte_zct(jm)  
      gtilcl(:)        = mmm_gti_zct(jm)         ;  gnhlcl(:)     = mmm_gnh_zct(jm)   
      gqlcl(:)         = mmm_gq_zct(jm)          ;  ahydlcl       = mmm_ahyd_zct(jm)
      vtorlcl(:)       = mmm_vtor_zct(jm)        ;  vpollcl(:)    = mmm_vpol_zct(jm)   
      vparlcl(:)       = mmm_vpar_zct(jm)        ;  gvtorlcl(:)   = mmm_gvtor_zct(jm)  
      gvpollcl(:)      = mmm_gvpol_zct(jm)       ;  gvparlcl      = mmm_gvpar_zct(jm)  
      gelonglcl        = mmm_gelong_zct(jm)      ;  tilcl(:)      = mmm_ti_zct(jm) 
      !OUTPUTS:
      xtilcl(:)        = zeroc                  ;  xdilcl(:)     = zeroc     
      xdzlcl(:)        = zeroc                  ;  xtelcl(:)     = zeroc 
      xvtlcl(:)        = zeroc                  ;  xvplcl(:)     = zeroc 
      xtiW20lcl(:)     = zeroc                  ;  xdiW20lcl(:)  = zeroc 
      xteW20lcl(:)     = zeroc                  ;  xtiDBMlcl(:)  = zeroc
      xdiDBMlcl(:)     = zeroc                  ;  xteETGlcl(:)  =  zeroc
      xteDBMlcl(:)     = zeroc                  
      gammaDBMlcl(:)   = zeroc                                            
      omegaDBMlcl(:)   = zeroc
      gammaW20lcl(:,:) = zeroc
      omegaW20lcl(:,:) = zeroc
      vconvlcl(:,:)    = zeroc  
      vfluxlcl(:,:)    = zeroc  
      IF(.NOT.  testing_mmm)THEN                                                  
         CALL mmm7_1 (                                                                      &  
              rmin = rminlcl,rmaj = rmajlcl,rmaj0 = mmm_rmaj0, elong = elonglcl,   &
              ne = nelcl, nh=nhlcl,nz=nzlcl,nf=nflcl,                              &
              zeff=zefflcl,te=telcl,ti=tilcl,q=qlcl,                               &
              btor=btorlcl,zimp=zimplcl,  aimp=aimplcl,   ahyd=ahydlcl,            &
              aimass=aimasslcl, wexbs=wexbslcl,                                    &
              gne=gnelcl,   gni=gnilcl,    gnh=gnhlcl,                             &
              gnz=gnzlcl,   gte=gtelcl,    gti=gtilcl,                             &
              gq=gqlcl,                                                            &
              vtor=vtorlcl,    vpol=vpollcl,    vpar=vparlcl,                      &
              gvtor=gvtorlcl,  gvpol=gvpollcl,  gvpar=gvparlcl,                    &   
              gelong = gelonglcl,                                                  & ! outputs
              xti=xtilcl, xdi= xdilcl, xte=xtelcl,xdz=xdzlcl,                      &
              xvt=xvtlcl,    xvp=xvplcl,                                           &
              xtiW20=xtiW20lcl, xdiW20=xdiW20lcl, xteW20 = xteW20lcl,              &
              xtiDBM=xtiDBMlcl,                                                    &
              xdiDBM=xdiDBMlcl, xteDBM=xteDBMlcl, xteETG=xteETGlcl,                &
              gammaW20=gammaW20lcl, omegaW20=omegaW20lcl,                          &
              gammaDBM=gammaDBMlcl,                                                &
              omegaDBM=omegaDBMlcl, vconv=vconvlcl,  vflux=vfluxlcl,               &
                                !.. Feature controls
              npoints=npoints, lprint=lprint, nprout=nprout, nerr=nerr,            &
              cmodel= mmm_cmodel, cswitch=mmm_cswitch,                             & 
              lswitch=mmm_lswitch )

         !effective ion thermal diffusivity, M^2/sec
         mmm_xti_zct(jm)          = xtilcl(1)                    ;  mmm_xdi_zct(jm)     = xdilcl(1)
         mmm_xdz_zct(jm)          = xdzlcl(1)                    ;  mmm_xte_zct(jm)     = xtelcl(1)
         mmm_xvt_zct(jm)          = xvtlcl(1)                    ;  mmm_xvp_zct(jm)     = xvplcl(1) 
         mmm_xtiW20_zct(jm)       = xtiW20lcl(1)                 ;  mmm_xdiW20_zct(jm)  = xdiW20lcl(1)  
         mmm_xteW20_zct(jm)       = xteW20lcl(1)                 ;  mmm_xtiDBM_zct(jm)  = xtiDBMlcl(1)
         mmm_xdiDBM_zct(jm)       = xdiDBMlcl(1)                 ;  mmm_xteETG_zct(jm)  = xteETGlcl(1)
         mmm_xteDBM_zct(jm)       = xteDBMlcl(1)                 
         mmm_gammaDBM_zct(jm)     = gammaDBMlcl(1)                                             
         mmm_omegaDBM_zct(jm)     = omegaDBMlcl(1) 

         mmm_gammaW20_zct(:,jm) = gammaW20lcl(:,1)
         mmm_omegaW20_zct(:,jm) = omegaW20lcl(:,1)
         mmm_vconv_zct(:,jm)    = vconvlcl(:,1)
         ! note vconvlcl(1:3,:) is returned as zero by  mmm7_1
         mmm_vflux_zct(:,jm)    = vfluxlcl(:,1)
         !---------------------------------------------------------------------------------------
         ! -- apparently the energy vflux output needs to be multiplied by density and temperature
         ! -- WhaT appears to be put out by mmm7_1 is JUST VFLUX = XT(E,I)*(-rDT/DR)(1/(R*T))
         ! --------------------------------------------------------------------------------------
         mmm_vflux_zct(1,jm) =  mmm_vflux_zct(1,jm)*mmm_ti_zct(jm)*mmm_nh_zct(jm)*joupkev ! W/m^2)
         mmm_vflux_zct(3,jm) =  mmm_vflux_zct(3,jm)*mmm_te_zct(jm)*mmm_ne_zct(jm)*joupkev ! W/m^2)
      ELSE    ! testing_mmm set some simple values
         !effective ion thermal diffusivity, M^2/sec
         mmm_xti_zct(jm)          = 1._DP                        ;  mmm_xdi_zct(jm)     = 1._DP
         mmm_xdz_zct(jm)          = zeroc                        ;  mmm_xte_zct(jm)     = 1._DP
         mmm_xvt_zct(jm)          = zeroc                        ;  mmm_xvp_zct(jm)     = zeroc 
         mmm_xtiW20_zct(jm)       = 1._DP                        ;  mmm_xdiW20_zct(jm)  = 1._DP 
         mmm_xteW20_zct(jm)       = 1._DP                        ;  mmm_xtiDBM_zct(jm)  = zeroc
         mmm_xdiDBM_zct(jm)       = zeroc                        ;  mmm_xteETG_zct(jm)  = zeroc
         mmm_xteDBM_zct(jm)       = zeroc                
         mmm_gammaDBM_zct(jm)     = zeroc                                            
         mmm_omegaDBM_zct(jm)     = zeroc

         mmm_gammaW20_zct(:,jm) = zeroc
         mmm_omegaW20_zct(:,jm) = zeroc
         mmm_vconv_zct(:,jm)    = zeroc
         !mmm_vflux_zct(:,jm)   = zeroc
         dtidrho = (ti(jm+1)- ti(jm) )/( r(jm+1)-r(jm) )
         dtedrho = (te(jm+1)- te(jm) )/( r(jm+1)-r(jm) )


         ! flux,ion thermal ,w/m^2
         mmm_vflux_zct(1,jm)       =  mmm_xti_zct(jm)*joupkev*1.e20* mmm_gti_zct(jm)*mmm_ti_zct(jm)/mmm_rmaj_ctr_zct(jm)
         !mmm_vflux_zct(1,jm)       = -mmm_xti_zct(jm)*joupkev*1.0e20*dtidrho   ! 8888888899999

         ! flux,hydrogenic,particle
         mmm_vflux_zct(2,jm)       =  mmm_xdi_zct(jm)*mmm_gnh_zct(jm)*mmm_nh_zct(jm)/mmm_rmaj_ctr_zct(jm)


         ! flux,electron  thermal ,w/m^2
         mmm_vflux_zct(3,jm)       =  mmm_xte_zct(jm)*joupkev*1.5e20*mmm_gte_zct(jm)*mmm_te_zct(jm)/mmm_rmaj_ctr_zct(jm) ! 8888888899999
         !mmm_vflux_zct(3,jm)       = - mmm_xte_zct(jm)*joupkev*1.5e20*dtedrho



         ! flux,impurities,particle
         mmm_vflux_zct(4,jm)       =  mmm_xdz_zct(jm)*mmm_gnz_zct(jm)*mmm_nz_zct(jm)/mmm_rmaj_ctr_zct(jm)

      ENDIF


      !------------------------------------------------------------------------------
      ! --  check for unstable growth rates. This info is used to determine if
      ! --  any variables need to be perturbed on this grid point, at this time step.
      !------------------------------------------------------------------------------
      IF(pertrb == izero)THEN
         IF(mmm_gammaDBM_zct(jm) .Gt. zeroc)mmm_unstable = .TRUE.
         DO  j=1,4
            IF(mmm_gammaW20_zct(j,jm) .GT. zeroc)mmm_unstable = .TRUE.
         ENDDO
      ENDIF




      ! restore the dependent variable to its base value
      ! the perturbations are recorded in mmm_val_pert
      IF(select_solver == 'nwt_pred_core')THEN
         IF(pertrb == 1) mmm_nh_zct(jm) = mmm_val_base_zct( 1,jm)
         IF(pertrb == 2) mmm_te_zct(jm) = mmm_val_base_zct( 2,jm)
         IF(pertrb == 3) mmm_ti_zct(jm) = mmm_val_base_zct( 3,jm)

         IF(pertrb == 5)THEN
            mmm_angrot_zct(jm) = mmm_val_base_zct(5,jm)
            mmm_vtor_zct(jm)   = dischg%rmajor*mmm_val_base_zct(5,jm)
            ! mmm_vexb_zct(jm)  = -mmm_theta_zct(jm)*mmm_vtor_zct(jm) - mmm_vdia_zct(jm)    &
            !                            + 2._DP* mmm_vpol_zct(jm)/hcap(jm) ! m/sec (1,nj)  ( Vexb =- Er/Bt )


         ENDIF

         IF(pertrb == -1) mmm_gnh_zct(jm) = -denhgrdm*mmm_rmaj_ctr_zct(jm)/mmm_nh_zct(jm)
         IF(pertrb == -2) mmm_gte_zct(jm) = dtegrdm_sv
         IF(pertrb == -3) mmm_gti_zct(jm) = dtigrdm_sv
!         IF(pertrb == -1) mmm_gnh_zct(jm) = mmm_val_base_zct( -1,jm)
!         IF(pertrb == -2) mmm_gte_zct(jm) = mmm_val_base_zct( -2,jm)
!         IF(pertrb == -3) mmm_gti_zct(jm) = mmm_val_base_zct( -3,jm)

         IF(pertrb == -5)THEN
            mmm_gangrot_zct(jm)           =  dwgrdm*mmm_rmaj_ctr_zct(jm)/mmm_angrot_zct(jm)
            mmm_gvtor_zct(jm)             =  dvtgrdm*mmm_rmaj_ctr_zct(jm)/mmm_angrot_zct(jm)
         ENDIF

      ENDIF
      mmm_loaded_ze = .FALSE.

      IF(nerr .ne. 0)THEN ! nerr is global to this routine
         lerrno = 57
         CALL terminate(lerrno,nlog)
      ENDIF

      RETURN

    END SUBROUTINE mmm_single_grid_point_zct
       



    SUBROUTINE mmm_single_grid_point(jm,pertrb)    
      !--------------------------------------------------------------------------
      ! single grid point interface to mmm
      !--------------------------------------------------------------------------

      USE nrtype,                                              ONLY : DP,I4B,I2B

      USE grid_class,                                          ONLY : nj,r2capi,r

      USE modmmm7_1,                                           ONLY : mmm7_1

      USE plasma_properties,                                   ONLY : dischg

      USE dep_var,                                             ONLY : ene,en,angrot

      USE io_gcnmp,                                            ONLY : nlog, ncrt

      USE ions_gcnmp,                                          ONLY : nprim,nion,atw,nimp

      USE error_handler,                                       ONLY : iomaxerr, lerrno,terminate 

      USE common_constants,                                    ONLY : izero,zeroc,Proton_Mass,   &
           joupkev

      USE modmmm7_1,                                           ONLY : mmm7_1

      IMPLICIT NONE

      INTEGER(I4B) i,j,k

      INTEGER(I2B) nprout

      REAL(DP) ena,mass,sqrteta,scfctr,drmin,denhgrdm,denhgrdp,dtegrdm,dtegrdp, &
           dtigrdm,dtigrdp,dwgrdm,dwgrdp,dvtgrdm,dvtgrdp,d2grd,dengradp,    &
           stepsize,dtidrho,dtedrho_new,dtidrho_new,dtedrmin,dtidrmin,      &
           dtegrdm_sv,dtigrdm_sv

      REAL(DP),DIMENSION(1) :: rminlcl,rmajlcl,elonglcl,nelcl,nhlcl,nzlcl,nflcl, &
           zefflcl,telcl,tilcl,qlcl,btorlcl,zimplcl,aimplcl, &
           ahydlcl,aimasslcl,wexbslcl,gnelcl,gnilcl,gnhlcl,  &
           gnzlcl,gtelcl,gtilcl,gqlcl,vtorlcl,vparlcl,       &
           vpollcl,gvtorlcl,gelonglcl,gvpollcl,gvparlcl,     &
           xtilcl,xdilcl,xdzlcl,xtelcl,xvtlcl,xvplcl,        &
           xtiW20lcl, xdiW20lcl,xteW20lcl, xtiDBMlcl,        &
           xdiDBMlcl,xteETGlcl,xteDBMlc,gammaDBMlcl,         &
           omegaDBMlcl,xteDBMlcl

      REAL(DP),DIMENSION(4,1) ::   gammaW20lcl,omegaW20lcl,vfluxlcl

      REAL(DP),DIMENSION(6,1) ::   vconvlcl 

      INTEGER,INTENT(IN) :: jm,pertrb 
      INTEGER(I4B) lcl_size,ok
      LOGICAL testing_mmm

      testing_mmm = .true. 
      testing_mmm =.false.
      IF(jm .lt. 1 .OR. jm .GE. nj)THEN
         lerrno = iomaxerr + 189
         CALL terminate(lerrno,nlog)
      ENDIF


      nprout = nlog 
      lprint = 1  ! > =6 much output
      sqrteta = 1.e-2_DP*mmm_scfctr ! for 1.e-n perturbation  will be in digit no n
      !-----------------------------------------------------------------------------
      ! -- perturbation values:
      !-----------------------------------------------------------------------------

      slsv: IF(select_solver == 'nwt_pred_cor')THEN
         drmin = mmm_rmin_ze(jm+1)-mmm_rmin_ze(jm)
         scfctr = mmm_scfctr*drmin/mmm_rmaj_ctr_ze(jm)
         IF(pertrb == izero )THEN
            lcl_size = 0 
            IF(ALLOCATED(mmm_val_base_ze)) lcl_size = SIZE(mmm_val_base_ze,2)
            IF(lcl_size .NE. nj-1)THEN
               DEALLOCATE(mmm_val_base_ze,stat = ok)   ! non fatal dealoc
               ALLOCATE(mmm_val_base_ze(-itran_max:itran_max,nj-1))
            ENDIF
            lcl_size = 0 
            IF(ALLOCATED(mmm_val_pert_ze)) lcl_size = SIZE(mmm_val_pert_ze,2)
            IF(lcl_size .NE. nj-1)THEN
               DEALLOCATE(mmm_val_pert_ze,stat = ok)   ! non fatal dealoc
               ALLOCATE(mmm_val_pert_ze(-itran_max:itran_max,nj-1))
            ENDIF

             IF(jm .GT. 1)THEN
                !note that val_pert is size of perturbation (not perturbed value)
                mmm_val_pert_ze(0,jm)   = zeroc  ! index 0 ==> base values,val_pert =0

                !mmm_val_base_ze(-5,jm)  = (mmm_angrot_ze(jm)/mmm_rmaj_ctr_ze(jm))*mmm_gangrot_ze(jm)

                mmm_val_base_ze(-5,jm)  = (mmm_angrot_ze(jm+1) - mmm_angrot_ze(jm-1))/(r(jm+1)-r(jm-1))

                !mmm_typ_val(-5) = ABS((mmm_angrot_ze(nj)-mmm_angrot_ze(1))/(mmm_rmin_ze(nj)-mmm_rmin_ze(1)))
                mmm_typ_val(-5) = ABS((mmm_angrot_ze(nj)-mmm_angrot_ze(1))/(r(nj)-r(1)))

                mmm_val_base_ze(-4,jm)   =  zeroc  
                mmm_typ_val(-4)          =  zeroc

                !mmm_val_base_ze(-3,jm)  =  (mmm_ti_ze(jm)/mmm_rmaj_ctr_ze(jm))*mmm_gti_ze(jm)
                mmm_val_base_ze(-3,jm)   =  (mmm_ti_ze(jm+1)- mmm_ti_ze(jm-1))/(r(jm+1)-r(jm-1))
                mmm_typ_val(-3)          =  ABS((mmm_ti_ze(nj)- mmm_ti_ze(1))/(r(nj)-r(1)))

               ! mmm_val_base_ze(-2,jm)  =  (mmm_te_ze(jm)/mmm_rmaj_ctr_ze(jm))*mmm_gte_ze(jm)
                 mmm_val_base_ze(-2,jm)  = (mmm_te_ze(jm+1)- mmm_te_ze(jm-1))/(r(jm+1)-r(jm-1))
                mmm_typ_val(-2)          =  ABS((mmm_te_ze(nj)- mmm_te_ze(1))/(r(nj)-r(1)))

                !mmm_val_base_ze(-1,jm)  = (mmm_nh_ze(jm)/mmm_rmaj_ctr_ze(jm))*mmm_gnh_ze(jm) 
                mmm_val_base_ze(-1,jm)   = (mmm_nh_ze(jm+1)- mmm_nh_ze(jm-1))/(r(jm+1)-r(jm-1))
                mmm_typ_val(-1)          = ABS((mmm_nh_ze(nj)- mmm_nh_ze(1))/(r(nj)- r(1)))

            ELSEIF(jm ==1)THEN
               mmm_val_base_ze(-5,jm)   = zeroc
               mmm_val_base_ze(-4,jm)   = zeroc
               mmm_val_base_ze(-3,jm)   = zeroc
               mmm_val_base_ze(-2,jm)   = zeroc
               mmm_val_base_ze(-1,jm)   = zeroc
            ENDIF

            mmm_val_base_ze( 0,jm)  =  zeroc              ;  mmm_typ_val(0)  = zeroc
            mmm_val_base_ze( 1,jm)  =  mmm_nh_ze(jm)      ;  mmm_typ_val(1)  = mmm_nh_ze(nj/2)
            mmm_val_base_ze( 2,jm)  =  mmm_te_ze(jm)      ;  mmm_typ_val(2)  = mmm_te_ze(nj/2)
            mmm_val_base_ze( 3,jm)  =  mmm_ti_ze(jm)      ;  mmm_typ_val(3)  = mmm_ti_ze(nj/2)
            mmm_val_base_ze( 4,jm)  =  zeroc              ;  mmm_typ_val(4)  = zeroc
            mmm_val_base_ze( 5,jm)  =  mmm_angrot_ze(jm)  ;  mmm_typ_val(5)  = ABS(mmm_angrot_ze(nj/2))

            ! note step size on gradients is derived from abs(pertrb):
            !stepsize                = sqrteta*MAX(ABS(mmm_val_base(ABS(pertrb),jm)),mmm_typ_val(ABS(pertrb))) &
            !                                       *SIGN(1.D0,mmm_val_base(ABS(pertrb),jm))

            !-----------------------------------------------------------------------------
            ! -- perturb one of the simulation mode variables
            ! -- the magnitude of the pertubation is done based on the current
            ! -- spatial rate of change of the dependent variable.(du/dr)
            ! -- note that pertrb .ne. 0 is active only if dependent variable 
            ! -- is in simulation mode.
            !-----------------------------------------------------------------------------
         ELSE IF      (pertrb == 1  )THEN

               mmm_val_pert_ze(pertrb,jm)= sqrteta*MAX(ABS(mmm_val_base_ze(pertrb,jm)),&
                    mmm_typ_val(pertrb))*SIGN(1.D0,mmm_val_base_ze(pertrb,jm))
 
            !perturb variable given to mmm7_1, restore to original value below:
            mmm_nh_ze(jm) = mmm_val_base_ze(pertrb,jm) + mmm_val_pert_ze(pertrb,jm)


         ELSE IF (pertrb == 2  )THEN
            !perturb variable given to mmm7_1, restore to original value below:
            mmm_val_pert_ze(pertrb,jm)   = sqrteta*MAX(ABS(mmm_val_base_ze(pertrb,jm)),mmm_typ_val(pertrb)) &
                 *SIGN(1.D0,mmm_val_base_ze(pertrb,jm))
            mmm_te_ze(jm) = mmm_val_base_ze(pertrb,jm) + mmm_val_pert_ze(pertrb,jm)



         ELSE IF (pertrb == 3  )THEN
            mmm_val_pert_ze(pertrb,jm)= sqrteta*MAX(ABS(mmm_val_base_ze(pertrb,jm)),&
                 mmm_typ_val(pertrb))*SIGN(1.D0,mmm_val_base_ze(pertrb,jm))
            !perturb variable given to mmm7_1, restore to original value below:
            mmm_ti_ze(jm) = mmm_val_base_ze( pertrb,jm) + mmm_val_pert_ze(pertrb,jm)



         ELSE IF (pertrb == 4  )THEN 


         ELSE IF (pertrb == 5  )THEN ! perturb omega (toroidal rot speed)
                    ! mag axis
               mmm_val_pert_ze(pertrb,jm)= sqrteta*MAX(ABS(mmm_val_base_ze(pertrb,jm)),&
                    mmm_typ_val(pertrb))*SIGN(1.D0,mmm_val_base_ze(pertrb,jm))
 
            !perturb variable given to mmm7_1, restore to original value below:
            ! vtor = R*omega:
            mmm_vtor_ze(jm) =  dischg%rmajor*(mmm_val_base_ze(pertrb,jm) + mmm_val_pert_ze(pertrb,jm))
            !mmm_vexb_ze(jm) = -mmm_theta_ze(jm)*mmm_vtor_ze(jm) - mmm_vdia_ze(jm)    &
            !                           + 2._DP* mmm_vpol_ze(jm)/hcap(jm)       ! m/sec (1,nj)  ( Vexb =- Er/Bt )
            !mmm_wexbs_ze(jm) = ??
            ! vpar not changed

         ELSE IF (pertrb == -1 )THEN  
            mmm_val_pert_ze(pertrb,jm)= sqrteta*MAX(ABS(mmm_val_base_ze(pertrb,jm)),&
                    mmm_typ_val(pertrb))*SIGN(1.D0,mmm_val_base_ze(pertrb,jm))
            denhgrdm                   =  mmm_gnh_ze(jm)*mmm_nh_ze(jm)/mmm_rmaj_ctr_ze(jm) ! base grad value at jm
            denhgrdp                   =  mmm_gnh_ze(jm+1)*mmm_nh_ze(jm+1)/mmm_rmaj_ctr_ze(jm+1)
            d2grd                      = (denhgrdp-denhgrdm)/drmin
            mmm_val_pert_ze(pertrb,jm) = d2grd*mmm_scfctr
            dengradp                   = denhgrdm + mmm_val_pert_ze(pertrb,jm)


            mmm_gnh_ze(jm)             = -dengradp*mmm_rmaj_ctr_ze(jm)/mmm_nh_ze(jm)
         ELSE IF (pertrb == -2 )THEN 
            mmm_val_pert_ze(pertrb,jm)   = sqrteta*MAX(ABS(mmm_val_base_ze(pertrb,jm)),mmm_typ_val(pertrb)) &
                 *SIGN(1.D0,mmm_val_base_ze(pertrb,jm))
               dtedrho_new =  mmm_val_base_ze(pertrb,jm) + mmm_val_pert_ze(pertrb,jm) ! this is perturbed dte/drho
              dtedrmin = dtedrho_new*mmm_drhodrmin_ze(jm)
              dtegrdm_sv = mmm_gte_zct(jm) ! enables restore below 
              ! this is what  mmm_7.1 wants as input:
            mmm_gte_ze(jm)              = -(mmm_rmaj_ctr_ze(jm)/mmm_te_ze(jm))*dtedrmin



         ELSE IF (pertrb == -3 )THEN
            mmm_val_pert_ze(pertrb,jm)   = sqrteta*MAX(ABS(mmm_val_base_ze(pertrb,jm)),mmm_typ_val(pertrb)) &
                 *SIGN(1.D0,mmm_val_base_ze(pertrb,jm)) 
              dtidrho_new =  mmm_val_base_ze(pertrb,jm) + mmm_val_pert_ze(pertrb,jm) ! this is perturbed dti/drho
              dtidrmin = dtidrho_new*mmm_drhodrmin_ze(jm)

               dtigrdm_sv = mmm_gti_ze(jm) ! enables restore below 
              ! this is what  mmm_7.1 wants as input:            
               mmm_gti_ze(jm)              = -(mmm_rmaj_ctr_ze(jm)/mmm_ti_ze(jm))*dtidrmin
 


         ELSE IF (pertrb == -4 )THEN

         ELSE IF (pertrb == -5 )THEN
            dwgrdm =  mmm_gangrot_ze(jm)*mmm_angrot_ze(jm)/mmm_rmaj_ctr_ze(jm) ! base grad value at jm
            dwgrdp =  mmm_gangrot_ze(jm+1)*mmm_angrot_ze(jm+1)/mmm_rmaj_ctr_ze(jm+1)
            d2grd = (dwgrdp-dwgrdm)/drmin
            mmm_val_pert_ze(pertrb,jm) = d2grd*mmm_scfctr
            dwgrdp                     = dwgrdm + mmm_val_pert_ze(pertrb,jm)
            mmm_gangrot_ze(jm)         = dwgrdp*mmm_rmaj_ctr_ze(jm)/mmm_angrot_ze(jm)

            dvtgrdm =  mmm_gvtor_ze(jm)*mmm_vtor_ze(jm)/mmm_rmaj_ctr_ze(jm) ! base grad value at jm
            dvtgrdp =  mmm_gvtor_ze(jm+1)*mmm_vtor_ze(jm+1)/mmm_rmaj_ctr_ze(jm+1)
            d2grd = (dvtgrdp-dvtgrdm)/drmin
            dvtgrdp                    = dvtgrdm + d2grd*mmm_scfctr
            mmm_gvtor_ze(jm)           = dvtgrdp*mmm_rmaj_ctr_ze(jm)/mmm_angrot_ze(jm)

         ELSE IF (pertrb .NE. 0)THEN
            lerrno = iomaxerr + 149
            CALL terminate(lerrno,nlog)
         ENDIF

      ENDIF  slsv


      !---------------------------------------------------------------------------
      ! mmm7_1  this call is setup to do a single grid point only:
      ! because of argument associatin we have to go through this crap:  HSJ
      !---------------------------------------------------------------------------
      IF(pertrb == izero)mmm_unstable = .FALSE.
      npoints = 1
      rminlcl(:)       = mmm_rmin_ze(jm)        ;  rmajlcl(:)    = mmm_rmaj_ctr_ze(jm)
      elonglcl(:)      = mmm_elong_ze(jm)       ;  nelcl(:)      = mmm_ne_ze(jm)  
      nhlcl(:)         = mmm_nh_ze(jm)
      nzlcl(:)         = mmm_nz_ze(jm)          ;  nflcl(:)      = mmm_nf_ze(jm)
      zefflcl(:)       = mmm_zeff_ze(jm)        ;  telcl(:)      = mmm_te_ze(jm)  
      qlcl(:)          = mmm_q_ze(jm)           ;  btorlcl(:)    = mmm_btor_ze(jm)     
      zimplcl(:)       = mmm_zimp_ze(jm)        ;  aimplcl       = mmm_aimp_ze(jm)  
      aimasslcl(:)     = mmm_aimass_ze(jm)      ;  wexbslcl(:)   = mmm_wexbs_ze(jm)
      gnelcl(:)        = mmm_gne_ze(jm)         ;  gnilcl(:)     = mmm_gni_ze(jm) 
      gnzlcl(:)        = mmm_gnz_ze(jm)         ;  gtelcl(:)     = mmm_gte_ze(jm)  
      gtilcl(:)        = mmm_gti_ze(jm)         ;  gnhlcl(:)     = mmm_gnh_ze(jm)   
      gqlcl(:)         = mmm_gq_ze(jm)          ;  ahydlcl       = mmm_ahyd_ze(jm)
      vtorlcl(:)       = mmm_vtor_ze(jm)        ;  vpollcl(:)    = mmm_vpol_ze(jm)   
      vparlcl(:)       = mmm_vpar_ze(jm)        ;  gvtorlcl(:)   = mmm_gvtor_ze(jm)  
      gvpollcl(:)      = mmm_gvpol_ze(jm)       ;  gvparlcl      = mmm_gvpar_ze(jm)  
      gelonglcl        = mmm_gelong_ze(jm)      ;  tilcl(:)      = mmm_ti_ze(jm) 
      !OUTPUTS:
      xtilcl(:)        = zeroc                  ;  xdilcl(:)     = zeroc     
      xdzlcl(:)        = zeroc                  ;  xtelcl(:)     = zeroc 
      xvtlcl(:)        = zeroc                  ;  xvplcl(:)     = zeroc 
      xtiW20lcl(:)     = zeroc                  ;  xdiW20lcl(:)  = zeroc 
      xteW20lcl(:)     = zeroc                  ;  xtiDBMlcl(:)  = zeroc
      xdiDBMlcl(:)     = zeroc                  ;  xteETGlcl(:)  =  zeroc
      xteDBMlcl(:)     = zeroc                  
      gammaDBMlcl(:)   = zeroc                                            
      omegaDBMlcl(:)   = zeroc
      gammaW20lcl(:,:) = zeroc
      omegaW20lcl(:,:) = zeroc
      vconvlcl(:,:)    = zeroc  
      vfluxlcl(:,:)    = zeroc  
      IF(.NOT.  testing_mmm)THEN                                                  
         CALL mmm7_1 (                                                                      &  
              rmin = rminlcl,rmaj = rmajlcl,rmaj0 = mmm_rmaj0, elong = elonglcl,   &
              ne = nelcl, nh=nhlcl,nz=nzlcl,nf=nflcl,                              &
              zeff=zefflcl,te=telcl,ti=tilcl,q=qlcl,                               &
              btor=btorlcl,zimp=zimplcl,  aimp=aimplcl,   ahyd=ahydlcl,            &
              aimass=aimasslcl, wexbs=wexbslcl,                                    &
              gne=gnelcl,   gni=gnilcl,    gnh=gnhlcl,                             &
              gnz=gnzlcl,   gte=gtelcl,    gti=gtilcl,                             &
              gq=gqlcl,                                                            &
              vtor=vtorlcl,    vpol=vpollcl,    vpar=vparlcl,                      &
              gvtor=gvtorlcl,  gvpol=gvpollcl,  gvpar=gvparlcl,                    &   
              gelong = gelonglcl,                                                  & ! outputs
              xti=xtilcl, xdi= xdilcl, xte=xtelcl,xdz=xdzlcl,                      &
              xvt=xvtlcl,    xvp=xvplcl,                                           &
              xtiW20=xtiW20lcl, xdiW20=xdiW20lcl, xteW20 = xteW20lcl,              &
              xtiDBM=xtiDBMlcl,                                                    &
              xdiDBM=xdiDBMlcl, xteDBM=xteDBMlcl, xteETG=xteETGlcl,                &
              gammaW20=gammaW20lcl, omegaW20=omegaW20lcl,                          &
              gammaDBM=gammaDBMlcl,                                                &
              omegaDBM=omegaDBMlcl, vconv=vconvlcl,  vflux=vfluxlcl,               &
                                !.. Feature controls
              npoints=npoints, lprint=lprint, nprout=nprout, nerr=nerr,            &
              cmodel= mmm_cmodel, cswitch=mmm_cswitch,                             & 
              lswitch=mmm_lswitch )

         !effective ion thermal diffusivity, M^2/sec
         mmm_xti_ze(jm)          = xtilcl(1)                    ;  mmm_xdi_ze(jm)     = xdilcl(1)
         mmm_xdz_ze(jm)          = xdzlcl(1)                    ;  mmm_xte_ze(jm)     = xtelcl(1)
         mmm_xvt_ze(jm)          = xvtlcl(1)                    ;  mmm_xvp_ze(jm)     = xvplcl(1) 
         mmm_xtiW20_ze(jm)       = xtiW20lcl(1)                 ;  mmm_xdiW20_ze(jm)  = xdiW20lcl(1)  
         mmm_xteW20_ze(jm)       = xteW20lcl(1)                 ;  mmm_xtiDBM_ze(jm)  = xtiDBMlcl(1)
         mmm_xdiDBM_ze(jm)       = xdiDBMlcl(1)                 ;  mmm_xteETG_ze(jm)  = xteETGlcl(1)
         mmm_xteDBM_ze(jm)       = xteDBMlcl(1)                 
         mmm_gammaDBM_ze(jm)     = gammaDBMlcl(1)                                             
         mmm_omegaDBM_ze(jm)     = omegaDBMlcl(1) 

         mmm_gammaW20_ze(:,jm) = gammaW20lcl(:,1)
         mmm_omegaW20_ze(:,jm) = omegaW20lcl(:,1)
         mmm_vconv_ze(:,jm)    = vconvlcl(:,1)
         mmm_vflux_ze(:,jm)    = vfluxlcl(:,1)


 
         !---------------------------------------------------------------------------------------
         ! -- apparently the energy vflux output needs to be multiplied by density and temperature
         ! -- WhaT appears to be put out by mmm7_1 is JUST VFLUX = XT(E,I)*(-rDT/DR)(1/(R*T))
         ! --------------------------------------------------------------------------------------
         mmm_vflux_ze(1,jm) =  mmm_vflux_ze(1,jm)*mmm_ti_ze(jm)*mmm_nh_ze(jm)*joupkev ! W/m^2)
         mmm_vflux_ze(3,jm) =  mmm_vflux_ze(3,jm)*mmm_te_ze(jm)*mmm_ne_ze(jm)*joupkev ! W/m^2)


      ELSE    ! testing_mmm set some simple values
         !effective ion thermal diffusivity, M^2/sec
         mmm_xti_ze(jm)          = 1._DP                        ;  mmm_xdi_ze(jm)     = 1._DP
         mmm_xdz_ze(jm)          = zeroc                        ;  mmm_xte_ze(jm)     = 1._DP
         mmm_xvt_ze(jm)          = zeroc                        ;  mmm_xvp_ze(jm)     = zeroc 
         mmm_xtiW20_ze(jm)       = 1._DP                        ;  mmm_xdiW20_ze(jm)  = 1._DP 
         mmm_xteW20_ze(jm)       = 1._DP                        ;  mmm_xtiDBM_ze(jm)  = zeroc
         mmm_xdiDBM_ze(jm)       = zeroc                        ;  mmm_xteETG_ze(jm)  = zeroc
         mmm_xteDBM_ze(jm)       = zeroc                
         mmm_gammaDBM_ze(jm)     = zeroc                                            
         mmm_omegaDBM_ze(jm)     = zeroc

         mmm_gammaW20_ze(:,jm) = zeroc
         mmm_omegaW20_ze(:,jm) = zeroc
         mmm_vconv_ze(:,jm)    = zeroc
         !mmm_vflux_ze(:,jm)   = zeroc
         dtidrho = (mmm_ti_ze(jm+1)- mmm_ti_ze(jm) )/( r(jm+1)-r(jm) )
         IF(jm .gt. 1)THEN
            dtidrho =0.5_DP*(mmm_ti_ze(jm+1)-mmm_ti_ze(jm-1))/( r(jm+1)-r(jm) )
            dtidrho = (mmm_ti_ze(jm+1)- mmm_ti_ze(jm-1) )/( r(jm+1)-r(jm-1) )
         else
            dtidrho =0.0_DP

         endif

         ! flux,ion thermal ,w/m^2
         mmm_vflux_ze(1,jm)       =  mmm_xti_ze(jm)*joupkev*1.e20* mmm_gti_ze(jm)*mmm_ti_ze(jm)/mmm_rmaj_ctr_ze(jm)
         mmm_vflux_ze(1,jm)       = -mmm_xti_ze(jm)*joupkev*1.0e20*dtidrho

         ! flux,hydrogenic,particle
         mmm_vflux_ze(2,jm)       =  mmm_xdi_ze(jm)*mmm_gnh_ze(jm)*mmm_nh_ze(jm)/mmm_rmaj_ctr_ze(jm)


         ! flux,electron  thermal ,w/m^2
         mmm_vflux_ze(3,jm)       =  mmm_xte_ze(jm)*joupkev*1.5e20*mmm_gte_ze(jm)*mmm_te_ze(jm)/mmm_rmaj_ctr_ze(jm)  ! 8888899999
         !mmm_vflux_ze(3,jm)       = - mmm_xte_ze(jm)*joupkev*1.5e20*dtedrho



         ! flux,impurities,particle
         mmm_vflux_ze(4,jm)       =  mmm_xdz_ze(jm)*mmm_gnz_ze(jm)*mmm_nz_ze(jm)/mmm_rmaj_ctr_ze(jm)

      ENDIF


      !------------------------------------------------------------------------------
      ! --  check for unstable growth rates. This info is used to determine if
      ! --  any variables need to be perturbed on this grid point, at this time step.
      !------------------------------------------------------------------------------
      IF(pertrb == izero)THEN
         IF(mmm_gammaDBM_ze(jm) .Gt. zeroc)mmm_unstable = .TRUE.
         DO  j=1,4
            IF(mmm_gammaW20_ze(j,jm) .GT. zeroc)mmm_unstable = .TRUE.
         ENDDO
      ENDIF




      ! restore the dependent variable to its base value
      ! the perturbations are recorded in mmm_val_pert
      IF(select_solver == 'nwt_pred_core')THEN
         IF(pertrb == 1) mmm_nh_ze(jm) = mmm_val_base_ze( 1,jm)
         IF(pertrb == 2) mmm_te_ze(jm) = mmm_val_base_ze( 2,jm)
         IF(pertrb == 3) mmm_ti_ze(jm) = mmm_val_base_ze( 3,jm)

         IF(pertrb == 5)THEN
            mmm_angrot_ze(jm) = mmm_val_base_ze(5,jm)
            mmm_vtor_ze(jm)   = dischg%rmajor*mmm_val_base_ze(5,jm)
            ! mmm_vexb_ze(jm)  = -mmm_theta_ze(jm)*mmm_vtor_ze(jm) - mmm_vdia_ze(jm)    &
            !                            + 2._DP* mmm_vpol_ze(jm)/hcap(jm) ! m/sec (1,nj)  ( Vexb =- Er/Bt )


         ENDIF

!         IF(pertrb == -1) mmm_gnh_ze(jm) = -denhgrdm*mmm_rmaj_ctr_ze(jm)/mmm_nh_ze(jm)
!         IF(pertrb == -2) mmm_gte_ze(jm) = -dtegrdm*mmm_rmaj_ctr_ze(jm)/mmm_te_ze(jm)
!         IF(pertrb == -3) mmm_gti_ze(jm) = -dtigrdm*mmm_rmaj_ctr_ze(jm)/mmm_ti_ze(jm)

         IF(pertrb == -2) mmm_gte_ze(jm) = dtegrdm_sv
         IF(pertrb == -3) mmm_gti_ze(jm) = dtigrdm_sv


!        IF(pertrb == -1) mmm_gnh_ze(jm) = mmm_val_base_ze( -1,jm)
!        IF(pertrb == -2) mmm_gte_ze(jm) = mmm_val_base_ze( -2,jm)
!        IF(pertrb == -3) mmm_gti_ze(jm) = mmm_val_base_ze( -3,jm)

         IF(pertrb == -5)THEN
            mmm_gangrot_ze(jm)           =  dwgrdm*mmm_rmaj_ctr_ze(jm)/mmm_angrot_ze(jm)
            mmm_gvtor_ze(jm)             =  dvtgrdm*mmm_rmaj_ctr_ze(jm)/mmm_angrot_ze(jm)
         ENDIF

      ENDIF

      mmm_loaded_ze = .TRUE.



      IF(nerr .ne. 0)THEN ! nerr is global to this routine
         lerrno = 57
         CALL terminate(lerrno,nlog)
      ENDIF

      RETURN

    END SUBROUTINE mmm_single_grid_point






    SUBROUTINE mmm_load_statefile_vectors(nj_local)
    !-----------------------------------------------------------------------------------
    ! -- load state file output arrays
    !-----------------------------------------------------------------------------------
      USE Plasma_properties,                                  ONLY : diffuse
      USE vector_class,                                       ONLY : new_Vector,      &
                                                                     delete_Vector_nf

      USE Grid_class,                                         ONLY : r

      USE solcon_gcnmp,                                       ONLY : use_mmm_flux

      USE error_handler,                                      ONLY : iomaxerr, lerrno,terminate

      USE io_gcnmp,                                           ONLY : nlog

      USE common_constants,                                   ONLY : zeroc,joupkev
      USE MPI_data,                                           ONLY : master,myid  ! for debug
       REAL(DP),DIMENSION(nj_local) :: work
       REAL(DP),DIMENSION(nj_local-1) ::rold_zc

      INTEGER(i4B)k, nj_local,oknf,ngW20,jj,njnew,nold,j,njoldin, nsp,kit,ksp


      njnew = nj_local ; nold = nj_local ! added
      IF(.NOT. mmm_loaded_ze)THEN
         ! if using _zct values then load  _ze values for output
         njnew = nj_local ; nold = nj_local
         DO j=1,nj_local-1
            rold_zc(j) = 0.5_DP*(r(j+1) +  r(j))
         ENDDO
         CALL mmm_zct2ze(nold,rold_zc,r,njnew)  
      ENDIF


      oknf = delete_Vector_nf(diffuse%mmm_xti)  ! non fatal delete
      !IF(select_solver == 'nwt_pred_cor')CALL zct2ze(mmm_xti_zct,mmm_xti_ze)
 
      diffuse%mmm_xti           = new_Vector(nj_local,mmm_xti_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xdi)
      diffuse%mmm_xdi           = new_Vector(nj_local,mmm_xdi_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xte)
      diffuse%mmm_xte           = new_Vector(nj_local,mmm_xte_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xdz)
      diffuse%mmm_xdz           = new_Vector(nj_local,mmm_xdz_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xvt)
      diffuse%mmm_xvt           = new_Vector(nj_local,mmm_xvt_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xvp)
      diffuse%mmm_xvp           = new_Vector(nj_local,mmm_xvp_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xtiW20)
      diffuse%mmm_xtiW20        = new_Vector(nj_local,mmm_xtiW20_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xdiW20)
      diffuse%mmm_xdiW20        = new_Vector(nj_local,mmm_xdiW20_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xteW20)
      diffuse%mmm_xteW20        = new_Vector(nj_local,mmm_xteW20_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xtiDBM)
      diffuse%mmm_xtiDBM        = new_Vector(nj_local,mmm_xtiDBM_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xdiDBM)
      diffuse%mmm_xdiDBM        = new_Vector(nj_local,mmm_xdiDBM_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xteDBM)
      diffuse%mmm_xteDBM        = new_Vector(nj_local,mmm_xteDBM_ze)

      oknf = delete_Vector_nf(diffuse%mmm_xteETG)
      diffuse%mmm_xteETG        = new_Vector(nj_local,mmm_xteETG_ze)

      oknf = delete_Vector_nf(diffuse%mmm_gammaDBM)
      diffuse%mmm_gammaDBM      = new_Vector(nj_local,mmm_gammaDBM_ze)

      oknf = delete_Vector_nf(diffuse%mmm_omegaDBM)
      diffuse%mmm_omegaDBM      = new_Vector(nj_local,mmm_omegaDBM_ze)
 
      ngW20 = 4
      IF(ASSOCIATED(diffuse%mmm_gammaW20))DEALLOCATE(diffuse%mmm_gammaW20)
      ALLOCATE(diffuse%mmm_gammaW20(ngW20))     
      IF(ASSOCIATED(diffuse%mmm_omegaW20))DEALLOCATE(diffuse%mmm_omegaW20)
      ALLOCATE(diffuse%mmm_omegaW20(ngW20))    
      IF(ASSOCIATED(diffuse%mmm_vflux))DEALLOCATE(diffuse%mmm_vflux)
      ALLOCATE(diffuse%mmm_vflux(ngW20))    
      IF(ASSOCIATED(diffuse%mmm_vconv))DEALLOCATE(diffuse%mmm_vconv)
      ALLOCATE(diffuse%mmm_vconv(ngW20+2))    
 
        DO jj =1,ngW20+2
           IF(jj .LE. ngW20)THEN
              work(:)                  = mmm_omegaW20_ze(jj,:)
              diffuse%mmm_omegaW20(jj) = new_Vector(nj_local,work)
              work(:)                  = mmm_gammaW20_ze(jj,:)
              diffuse%mmm_gammaW20(jj) = new_Vector(nj_local,work)
              work(:)                  = mmm_vflux_ze(jj,:)
              diffuse%mmm_vflux(jj)    = new_Vector(nj_local,work)
              ! If the flux was based on diffusivity type assumption we have to
              ! use the following output for Consistency:
              IF(jj == 1 .AND. use_mmm_flux(3) == 0)THEN ! ion thermal flux
                 k = SIZE(mmm_e_flux_zc,1)
                 njoldin = k+1
                 IF( k == nj_local-1)THEN
                    work(1) =  zeroc
                    DO j = 2,nj_local-1
                       work(j) =  joupkev*0.5_DP*(mmm_e_flux_zc(j-1,2,0)+ mmm_e_flux_zc(j,2,0))
                    ENDDO
                    work(nj_local) =  joupkev*mmm_e_flux_zc(k,2,0)
                 ELSE
                    nsp = 2      ! ion thermal flux
                    ksp = SIZE(mmm_e_flux_zc,2)
                    kit = lbound(mmm_e_flux_zc,3)
                    CALL mmm_flux_processing(njoldin,njnew,mmm_e_flux_zc,k,ksp,kit,nsp,work)
                    !lerrno = 62
                    !CALL terminate(lerrno,nlog)
                 ENDIF
                 diffuse%mmm_vflux(jj)    = new_Vector(nj_local,work)
              ENDIF
 
              IF(jj == 2 .AND. use_mmm_flux(1) == 0)THEN ! hydrogenic  particle flux
                 k = SIZE(mmm_p_flux_zc,1)
                 njoldin = k+1
                 IF( k == nj_local-1)THEN
                    work(1) = zeroc
                    DO j=2,nj_local-1
                       work(j) =  0.5_DP*(mmm_p_flux_zc(j-1,2,0) + mmm_p_flux_zc(j,2,0))
                    ENDDO
                    work(nj_local) = mmm_p_flux_zc(k,2,0)
                 ELSE
                    nsp = 2      ! ion particle flux
                    ksp = SIZE(mmm_p_flux_zc,2)
                    kit = SIZE(mmm_p_flux_zc,3)
                    CALL mmm_flux_processing(njoldin,njnew,mmm_p_flux_zc,k,ksp,kit,nsp,work)
                    !lerrno = 62
                    !CALL terminate(lerrno,nlog)
                 ENDIF
                 diffuse%mmm_vflux(jj)    = new_Vector(nj_local,work)
              ENDIF
       if(myid == master) write(6,FMT='("mmm_load_statefile_vectors,line3215")')
              IF(jj == 3 .AND. use_mmm_flux(2) == 0)THEN ! elect thermal
                 k = SIZE(mmm_e_flux_zc,1)
                 njoldin = k+1
                 IF( k == nj_local-1)THEN
                    work(1) = zeroc
                    DO j=2,nj_local-1
                       work(j) =  joupkev*0.5_DP*(mmm_e_flux_zc(j-1,3,0) + mmm_e_flux_zc(j,3,0))
                    ENDDO
                    work(nj_local) = joupkev*mmm_e_flux_zc(k,3,0)
       if(myid == master) write(6,FMT='("mmm_load_statefile_vectors,line3225")')
                 ELSE
                    nsp = 3      ! elct thermal flux
                    ksp = SIZE(mmm_e_flux_zc,2)
                    kit = SIZE(mmm_e_flux_zc,3)
       if(myid == master) write(6,FMT='("mmm_load_statefile_vectors,line3230")')
                    CALL mmm_flux_processing(njoldin,njnew,mmm_e_flux_zc,k,ksp,kit,nsp,work)
                    !lerrno = 62
                    !CALL terminate(lerrno,nlog)

                 ENDIF
                    diffuse%mmm_vflux(jj)    = new_Vector(nj_local,work)
 
              ENDIF
 
              IF(jj == 4 .AND. use_mmm_flux(1) == 0)THEN ! impurity particle flux
                 k = SIZE(mmm_p_flux_zc,1)
                 njoldin = k+1
                 IF( k == nj_local-1)THEN
                    DO j=2,nj_local-1
                       work(j) = 0.5_DP*( mmm_p_flux_zc(j-1,4,0) + mmm_p_flux_zc(j,4,0))       ! impurity ion flux m^2/sec
                    ENDDO
                    work(nj_local) = mmm_p_flux_zc(k,4,0)
                 ELSE
                    nsp = 4      ! imp partcl
                    ksp = SIZE(mmm_p_flux_zc,2)
                    kit = SIZE(mmm_p_flux_zc,3)
                    CALL mmm_flux_processing(njoldin,njnew,mmm_p_flux_zc,k,ksp,kit,nsp,work)
                    !lerrno = 62
                    !CALL terminate(lerrno,nlog)
                 ENDIF
                 diffuse%mmm_vflux(jj)    = new_Vector(nj_local,work)
              ENDIF
           ENDIF
       if(myid == master) write(6,FMT='("mmm_load_statefile_vectors,line3259")')

           work(:)                     = mmm_vconv_ze(jj,:)
           diffuse%mmm_vconv(jj)       = new_Vector(nj_local,work)
        ENDDO
      RETURN 

    END SUBROUTINE mmm_load_statefile_vectors




   SUBROUTINE mmm_zct2ze(nold,rold_zc,r,njnew)
!-------------------------------------------------------------------------------------
! -- zone center(_zct) to zone edge (_ze) conversion
!-------------------------------------------------------------------------------------

         IMPLICIT NONE
         REAL(DP), ALLOCATABLE,DIMEnSION(:)   :: work_temp_zc,work_temp_ze
         REAL(DP), ALLOCATABLE,DIMENSION(:,:) :: work_mmm_ze,work_mmm_zc
         INTEGER(I4B) j,nold,njnew
         REAL(DP) rold_zc(nold),r(njnew)

#ifndef ONETWO
         INTERFACE
            SUBROUTINE resize_extrap (f_ze,f_zc,nold_ze,rold_zc,rnew_ze,nnew_ze)

              USE NRTYPE ,         ONLY : I4B,DP 

              IMPLICIT  NONE
              REAL(DP),DIMENSION(:),ALLOCATABLE :: f_zc,f_ze
              REAL(DP),DIMENSION(nold_ze)       :: rold_ze
              REAL(DP),DIMENSION(nold_ze-1)     :: rold_zc
              REAL(DP),DIMENSION(nnew_ze)       :: rnew_ze
              REAL(DP) fl,fr,factor
              INTEGER(I4B) j,nnew_ze,nold_ze
            END SUBROUTINE resize_extrap
        END INTERFACE  



                CALL resize_extrap(mmm_gammaDBM_ze,mmm_gammaDBM_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_omegaDBM_ze,mmm_omegaDBM_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xdi_ze,mmm_xdi_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xti_ze,mmm_xti_zct,nold,rold_zc,r,njnew)

                CALL resize_extrap(mmm_xte_ze,mmm_xte_zct,nold,rold_zc,r,njnew)

                CALL resize_extrap(mmm_xdz_ze,mmm_xdz_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xvt_ze,mmm_xvt_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xvp_ze,mmm_xvp_zct,nold,rold_zc,r,njnew)
 
                CALL resize_extrap(mmm_xtiW20_ze,mmm_xtiW20_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xteW20_ze,mmm_xteW20_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xdiW20_ze,mmm_xdiW20_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xtiDBM_ze,mmm_xtiDBM_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xteDBM_ze,mmm_xteDBM_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xdiDBM_ze,mmm_xdiDBM_zct,nold,rold_zc,r,njnew)
                CALL resize_extrap(mmm_xteETG_ze,mmm_xTeETG_zct,nold,rold_zc,r,njnew)

             ! process mmm_vconv:
                IF(ALLOCATED(work_mmm_ze))DEALLOCATE(work_mmm_ze)
                ALLOCATE(work_mmm_ze(SIZE(mmm_vconv_zct,1),njnew))
                IF(ALLOCATED(work_mmm_zc))DEALLOCATE(work_mmm_zc)
                ALLOCATE(work_mmm_zc(SIZE(mmm_vconv_zct,1),njnew-1))
                DO j=1,SIZE(mmm_vconv_zct,1)
                   IF(ALLOCATED(work_temp_zc))DEALLOCATE(work_temp_zc)  ! needed because size changes in resize_extrap
                   ALLOCATE(work_temp_zc(SIZE(mmm_vconv_zct,2)))
                   IF(ALLOCATED(work_temp_ze))DEALLOCATE(work_temp_ze)
                   ALLOCATE(work_temp_ze(SIZE(mmm_vconv_zct,2)+1))
                   work_temp_zc(:) =  mmm_vconv_zct(j,:)

                   CALL resize_extrap(work_temp_ze,work_temp_zc,nold,rold_zc,r,njnew)
                   work_mmm_zc(j,:) = work_temp_zc(:)
                   work_mmm_ze(j,:) = work_temp_ze(:)
                ENDDO

                IF(ALLOCATED(mmm_vconv_zct))DEALLOCATE(mmm_vconv_zct)
                ALLOCATE(mmm_vconv_zct(SIZE(work_mmm_zc,1),SIZE(work_mmm_zc,2)))
                mmm_vconv_zct(:,:) = work_mmm_zc(:,:) 

                IF(ALLOCATED(mmm_vconv_ze))DEALLOCATE(mmm_vconv_ze)
                ALLOCATE(mmm_vconv_ze(SIZE(work_mmm_ze,1),SIZE(work_mmm_ze,2)))
                mmm_vconv_ze(:,:)  = work_mmm_ze(:,:)

                DEALLOCATE(work_mmm_ze,work_mmm_zc)
                DEALLOCATE(work_temp_zc)


             ! process mmm_omegaW20:
                IF(ALLOCATED(work_mmm_ze))DEALLOCATE(work_mmm_ze)
                ALLOCATE(work_mmm_ze(SIZE(mmm_omegaW20_zct,1),njnew))
                IF(ALLOCATED(work_mmm_zc))DEALLOCATE(work_mmm_zc)
                ALLOCATE(work_mmm_zc(SIZE(mmm_omegaW20_zct,1),njnew-1))
                DO j=1,SIZE(mmm_omegaW20_zct,1)
                   IF(ALLOCATED(work_temp_zc))DEALLOCATE(work_temp_zc) ! needed because size changes in resize_extrap
                   ALLOCATE(work_temp_zc(SIZE(mmm_omegaW20_zct,2)))
                   IF(ALLOCATED(work_temp_ze))DEALLOCATE(work_temp_ze)
                   ALLOCATE(work_temp_ze(SIZE(mmm_omegaW20_zct,2)+1))
                   work_temp_zc(:) =  mmm_omegaW20_zct(j,:)

                   CALL resize_extrap(work_temp_ze,work_temp_zc,nold,rold_zc,r,njnew)
                   work_mmm_zc(j,:) = work_temp_zc(:)
                   work_mmm_ze(j,:) = work_temp_ze(:)
                ENDDO

                IF(ALLOCATED(mmm_omegaW20_zct))DEALLOCATE(mmm_omegaW20_zct)
                ALLOCATE(mmm_omegaW20_zct(SIZE(work_mmm_zc,1),SIZE(work_mmm_zc,2)))
                mmm_omegaW20_zct(:,:) = work_mmm_zc(:,:) 

                IF(ALLOCATED(mmm_omegaW20_ze))DEALLOCATE(mmm_omegaW20_ze)
                ALLOCATE(mmm_omegaW20_ze(SIZE(work_mmm_ze,1),SIZE(work_mmm_ze,2)))
                mmm_omegaW20_ze(:,:)  = work_mmm_ze(:,:)
 
                DEALLOCATE(work_mmm_ze,work_mmm_zc)
                DEALLOCATE(work_temp_zc)

             ! process mmm_gammaW20:
                IF(ALLOCATED(work_mmm_ze))DEALLOCATE(work_mmm_ze)
                ALLOCATE(work_mmm_ze(SIZE(mmm_gammaW20_zct,1),njnew))
                IF(ALLOCATED(work_mmm_zc))DEALLOCATE(work_mmm_zc)
                ALLOCATE(work_mmm_zc(SIZE(mmm_gammaW20_zct,1),njnew-1))
                DO j=1,SIZE(mmm_gammaW20_zct,1)
                   IF(ALLOCATED(work_temp_zc))DEALLOCATE(work_temp_zc) ! needed because size changes in resize_extrap
                   ALLOCATE(work_temp_zc(SIZE(mmm_gammaW20_zct,2)))
                   IF(ALLOCATED(work_temp_ze))DEALLOCATE(work_temp_ze)
                   ALLOCATE(work_temp_ze(SIZE(mmm_gammaW20_zct,2)+1))
                   work_temp_zc(:) =  mmm_gammaW20_zct(j,:)

                   CALL resize_extrap(work_temp_ze,work_temp_zc,nold,rold_zc,r,njnew)
                   work_mmm_zc(j,:) = work_temp_zc(:)
                   work_mmm_ze(j,:) = work_temp_ze(:)
                ENDDO

                IF(ALLOCATED(mmm_gammaW20_zct))DEALLOCATE(mmm_gammaW20_zct)
                ALLOCATE(mmm_gammaW20_zct(SIZE(work_mmm_zc,1),SIZE(work_mmm_zc,2)))
                mmm_gammaW20_zct(:,:) = work_mmm_zc(:,:) 

                IF(ALLOCATED(mmm_gammaW20_ze))DEALLOCATE(mmm_gammaW20_ze)
                ALLOCATE(mmm_gammaW20_ze(SIZE(work_mmm_ze,1),SIZE(work_mmm_ze,2)))
                mmm_gammaW20_ze(:,:)  = work_mmm_ze(:,:)
 
                DEALLOCATE(work_mmm_ze,work_mmm_zc)
                DEALLOCATE(work_temp_zc)


            ! process mmm_vflux:

                IF(ALLOCATED(work_mmm_ze))DEALLOCATE(work_mmm_ze)
                ALLOCATE(work_mmm_ze(SIZE(mmm_vflux_zct,1),njnew))
                IF(ALLOCATED(work_mmm_zc))DEALLOCATE(work_mmm_zc)
                ALLOCATE(work_mmm_zc(SIZE(mmm_vflux_zct,1),njnew-1))
                DO j=1,SIZE(mmm_vflux_zct,1)
                   IF(ALLOCATED(work_temp_zc))DEALLOCATE(work_temp_zc) ! needed because size changes in resize_extrap
                   ALLOCATE(work_temp_zc(SIZE(mmm_vflux_zct,2)))
                   IF(ALLOCATED(work_temp_ze))DEALLOCATE(work_temp_ze)
                   ALLOCATE(work_temp_ze(SIZE(mmm_vflux_zct,2)+1))
                   work_temp_zc(:) =  mmm_vflux_zct(j,:)

                   CALL resize_extrap(work_temp_ze,work_temp_zc,nold,rold_zc,r,njnew)
                   work_mmm_zc(j,:) = work_temp_zc(:)
                   work_mmm_ze(j,:) = work_temp_ze(:)
                ENDDO

                IF(ALLOCATED(mmm_vflux_zct))DEALLOCATE(mmm_vflux_zct)
                ALLOCATE(mmm_vflux_zct(SIZE(work_mmm_zc,1),SIZE(work_mmm_zc,2)))
                mmm_vflux_zct(:,:) = work_mmm_zc(:,:) 

                IF(ALLOCATED(mmm_vflux_ze))DEALLOCATE(mmm_vflux_ze)
                ALLOCATE(mmm_vflux_ze(SIZE(work_mmm_ze,1),SIZE(work_mmm_ze,2)))
                mmm_vflux_ze(:,:)  = work_mmm_ze(:,:)
 
                DEALLOCATE(work_mmm_ze,work_mmm_zc)
                DEALLOCATE(work_temp_zc,work_temp_ze)


                mmm_loaded_ze = .TRUE.

#endif
        RETURN

   END SUBROUTINE mmm_zct2ze


 END MODULE mmm_in
