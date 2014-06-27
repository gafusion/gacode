
  SUBROUTINE gcnmp_namelist_defaults
! ----------------------------------------------------------------------
! --- set some defaults before reading namelist:
! ------------------------------------------------------------HSJ-------

    USE nrtype,          ONLY : DP,I4B,I2B

    USE bc_values_gcnmp, ONLY : u_vloop_bc,axis_bc_type,robin_bc, &
                                conserve_bc_flux,flux_mult_bc

    USE common_constants,ONLY : izero,zeroc

    USE curden_terms,    ONLY : rbsaxis

    USE error_handler 

    USE MPI_data,        ONLY : parallel_model

    USE solcon_gcnmp,    ONLY  : time0,time_max,dt,dtmin,dtmax,     &
                                dtstart,wdelt,max_steps,vloopvb,    &
                                use_constd,use_consts,use_glf23,    &
                                use_constconv,use_avg_glf_chi,      &
                                use_avg_chi,set_chie_chii,theta,    &
                                nsteps_plot,plot_time_incr,         &
                                save_incr_plot_file,iters_freeze,   &
                                nsteps_restart,bwcycle,max_bwcycle, &
                                max_bandwidth, min_bandwidth,       &
                                cycle_bandwidth, cgd,cgexp,         &
                                crit_grad, save_incr_restart_file,  &
                                ts_smfactor,use_glf23_flux,         &
                                single_density_simulation,disspc,   &
                                use_hyperbolic_den_eq,comp_var,     &
                                reuse_tglf,delta_u_crit,tglfvb,     &
                                outputvb,conve_mult,convi_mult,     &
                                use_stab_flux,asym_time,jac_skip,   &
                                asym_dt_max,select_solver,          &
                                use_flow_eq,set_parabolic_profiles, &
                                angrcple,it_imp_steps,prt_turbm,    &
                                core_iters_max,ambip_diffusion,     &
                                tol_corrector,rconsts_min,          &
                                rconsts_max,use_mmm_flux,use_mmm,   &
                                use_qrad_input,cyclotron_reflection,&
                                use_mmm_pinch,use_mmm_vpol_input,   &
                                use_mmm_vpar_input,use_mmm_er_input,&
                                elc_eng_fixed_transpt,              &
                                ion_eng_fixed_transpt,              &
                                ion_den_fixed_transpt,              &
                                ion_torrot_fixed_transpt,           &
                                curden_fixed_transpt,               &
                                use_input_confinement,              &
                                frz_tglf,frz_glf,frz_mmm,           &
                                use_frz_chi,                        &
                                fixed_profile_multiplier,           &
                                time_deriv_in_flux

    USE  neo_tport,            ONLY : jneo,neocl_mult,use_forcebal_chie, &
                                      use_forcebal_chii

    USE  io_gcnmp,             ONLY : rpc,iterdb_12_input_file,fixt_filename

    USE source_terms_gcnmp,    ONLY : beam_th_mult, stfuse_mult, sbfuse_mult, &
                                      qrad_mult,qfuse_mult_e,qfuse_mult_i,    &
                                      qrf_mult_e,qrf_mult_i,qbeam_mult_e,     &
                                      qbeam_mult_i

    USE fusion_gcnmp,          ONLY : differential_burnup,internal_thermal_fusion, &
                                      alpha_slow

    USE grid_class,            ONLY : nj,nj_out,nj_max,nj_start,allow_regrid, &
                                      use_compact_schemes

    USE tglfin,                ONLY : tglf_v_mult,tglf_d_mult,tglf_version_no

    USE glf23_gcnmp,           ONLY : glf_v_mult,glf_d_mult,dv_method

    USE mmm_in,                ONLY : mmm_v_mult,mmm_d_mult,mmm_scfctr, &
                                      use_mmm_zct

    USE turb_tport,            ONLY : wturbd,wturbv,gmax_crit

    USE Plasma_properties,     ONLY : mhd_dat

    USE ions_gcnmp,            ONLY : nion

    USE dep_var,               ONLY : dp4

    IMPLICIT NONE
 

    lerrno = izero ;  dtstart = 1.E-5_DP ; wdelt = 1.0_DP ;max_steps = 1
    dtmax = 0.010_DP  ; dtmin = 1.0E-06_DP
    u_vloop_bc = izero ; vloopvb =izero ;  rbsaxis = zeroc
    time0    = -1._DP  ! a value <0.0 means that tGCNMs given in iterdb file
                       ! will be used for the start time of the analysis
    time_max = -1._DP  ! a value <0.0 means that tGCNMf given in iterdb file
                       ! will be used for the end  time of the analysis
    neocl_mult(:,:)     = 1._DP
    ts_smfactor         = zeroc
    nsteps_plot         = izero
    nsteps_restart      = izero
    save_incr_plot_file = izero
    save_incr_restart_file  = izero
    set_parabolic_profiles  = izero       ! turn off writting parabolic profiles
    plot_time_incr          = zeroc
    use_hyperbolic_den_eq   = izero
    use_stab_flux       = .FALSE.
    use_qrad_input      = izero           ! do not use qrad input from statefile
    use_forcebal_chie   = .FALSE.
    use_forcebal_chii   = .FALSE.
    asym_time           = HUGE(time_max)
    asym_dt_max         = 0.0001_DP
    jneo                = izero
    jac_skip            = 1_I4B
    use_constd (:)      = izero
    use_consts (:)      = izero
    use_glf23(:)        = izero 
    use_mmm(:)          = izero
    use_constconv(:)    = izero 
    use_glf23_flux(:)   = izero
    use_mmm_flux(:)     = izero
    use_mmm_vpol_input  = izero
    use_mmm_vpar_input  = izero
    use_mmm_er_input    = izero
    mmm_scfctr          = 0.05 ! used to control perturbations in nwt_pred_core
    !use_mmm_zct         = .TRUE.  set in mmm_input_defaults.f90
    use_avg_glf_chi     = izero 
    use_avg_chi         = izero
    set_chie_chii       = zeroc
    axis_bc_type        = izero
    use_compact_schemes = izero
    rpc                 = izero
    prt_turbm           = 1000
    iters_freeze        = 10000000
    nj_out              = nj
    nj_start            = nj	
    allow_regrid        = .TRUE.

    rconsts_min(:)      = -1._DP  ! const source (if selected)
    rconsts_max(:)      =  2._DP  ! is on from 0<= r <=1.

    nj_max              = 257   ! must be 2^n+1 for integer n
    theta               = 1.0_DP
    conve_mult          = 1.0_DP
    convi_mult          = 1.0_DP
    angrcple            = 1.0_DP
    tol_corrector       = 1.e-4_DP
    parallel_model      = izero !parallel_model selects type of parallel 
                                !computationm to be done
                                !at present only parallel_model = 0 is working.
                                !parallel_model  = 0 standard domain 
                                !decompostion over the grid(default)
!                       = 1     above  plus independent columns of jacobian model
!                       = 2     OpenMP model
!                       = 3     ....?? not yet defined but will deal with parallel
!                               solution of systems of equations beyond what is
!                               done above.


     bwcycle         = izero
     max_bwcycle     = 3_I4B
     max_bandwidth   = 5            !NOTE that this is addition to the
     min_bandwidth   = izero        !standard tridiagoanl system bandwidth
     cycle_bandwidth = .FALSE. 
     cgd(:) = 0.0_DP ; cgexp(:) = 0.0_DP ; crit_grad(:)  = 0.0_DP
     comp_var(:)     = izero
     iterdb_12_input_file =' '
     single_density_simulation = .FALSE.
     beam_th_mult    = 1._DP
     stfuse_mult     = 1._DP
     sbfuse_mult     = 1._DP
     qfuse_mult_e    = 1._DP
     qfuse_mult_i    = 1._DP
     qrf_mult_e      = 1._DP
     qrf_mult_i      = 1._DP
     qrad_mult       = 1._DP
     cyclotron_reflection = 1._DP  ! means  all cyclotron radiation is reabosrbed
     qbeam_mult_e    = 1._DP
     qbeam_mult_i    = 1._DP
     differential_burnup = .FALSE.
     internal_thermal_fusion = .TRUE.
     alpha_slow = 'starting'        ! means slowing down of alphas starts at current time.
                                    ! use 'asymptotic' to assume slowig down is in
                                    ! in steady state 
     disspc(:) = zeroc
     reuse_tglf = 1                 ! this means tglf will be called on every time step
                                    ! set to number > 1 to skip some tglf calls
     delta_u_crit    = 0.01_DP
     outputvb        = izero
     tglfvb          = izero
     tglf_v_mult     = 1._DP   ; tglf_d_mult = 1._DP ; tglf_version_no = 1.94    
     glf_v_mult      = 1._DP   ; glf_d_mult  = 1._DP
     mmm_v_mult      = 1._DP   ; mmm_d_mult  = 1._DP
     use_mmm_pinch = izero
     select_solver   = 'none'
     dv_method       = izero   ! do not change (there is no dv_method)
     use_flow_eq     = -1
     it_imp_steps    = izero  ! if nwt_pred_core is used then turn off iterative improvement of turbulent models
     core_iters_max  = izero   ! allow this many iterative improvement steps . Note that 
                           ! we can carry out iterative improvement without invoking further changes in the
                           ! turbulent confinement 

    wturbd(:,:) = 1._dp
    wturbv(:,:) = 1._dp
    gmax_crit   = 0.05     ! max allowed relative change in gradient
    
    ambip_diffusion = .TRUE.

    mhd_dat%rhon_betan_ped = -1.   ! -1 means not used
    mhd_dat%betan_ped      = -1.



   !---------------------------------------------------------------------
   ! For running equations with fixed diffusivities  we use the following 
   ! switches. Switches will be effective only if corresponding profile 
   ! is run in simulation mode.
   ! The INPUT value is a number ,delta_t,
   ! which measures the elapsed time since the begining of the run
   ! The corresponding chi or flux  becomes frozen when 
   !       time = time0 + delta_t
   ! For time0 < time < time0+delta_t the profile is evolved with whatever
   ! transport model is in effect.
   ! For time0 + delta_t < time <= timmax the corresponding profile is fixed
   ! at the value in effect at time = time0 + delta_t 
   ! The default values are negative numbers so that the frozen profiles are not
   ! used.
   ! To run with the diffusivity frozen at the inital value,
   ! taken from the statefile input, set delta_t = 0.0
   ! 
   ! The ion_den_fixed transport APPLIES TO ALL NPRIM ION SPECIES
   ! SIMULTANEOUSLY, individual selection is not implemented.
   ! It is the users responsibility to align the density input properly
   ! For example if the original statefile used to generate the 
   ! fluxes was run with species 'dt' but the current run
   ! with fixed transport has in it species 'd' or 't' or 'd' and 't'
   ! then the constructed chis  will be WRONG!! The code has no way of
   ! knowing under what circumstances the original statefile was created.
   ! Of course normally you wouldnt mix and match statefiles in this way
   ! but you were warned.
   !---------------------------------------------------------------HSJ------
        elc_eng_fixed_transpt              = .FALSE.
        ion_eng_fixed_transpt              = .FALSE.
        ion_den_fixed_transpt              = .FALSE.
        ion_torrot_fixed_transpt           = .FALSE.
        curden_fixed_transpt               = .FALSE.
        use_input_confinement              = .FALSE.
        frz_mmm                            = .FALSE.
        frz_glf                            = .FALSE.
        frz_tglf                           = .FALSE.
        use_frz_chi                        = .TRUE.



        dbg_print                          =  .FALSE.

        fixt_filename = 'NONE'
       
        robin_bc                           = .FALSE.

        conserve_bc_flux                   = .FALSE.

        time_deriv_in_flux                 = .FALSE.

        flux_mult_bc(:)                    = 1._DP

    CALL set_forcebal_default_switches   ! turns off call to forcebal

    RETURN


  END SUBROUTINE gcnmp_namelist_defaults
