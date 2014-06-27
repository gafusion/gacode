  MODULE GCNMP_namelist

    USE nrtype,   ONLY : DP,I4B,I2B

    USE error_handler

    USE Plasma_properties ,  ONLY : mhd_dat,pellet,equil_type_size

    USE solcon_gcnmp,        ONLY : time0,time,time_max,dt,dtmin,                &
                               dtmax,use_glf23, dtstart,max_steps,itte,itti,     &
                               itxj,itw,itenp,iteni,freeze_type,tGCNMf,tGCNMs,   &
                               switch_method,tot_iters_max,non_lin_method,       &
                               fdigits,gradtol,steptol,bandwidth,jacobian_type,  &
                               fvectol,gradtol,ssqrmin,maxfev,steady_state,      &
                               conv_skip,wrt_nwt,outputvb,theta,                 &
                               iters_freeze,random_pert,jac_skip,tglfvb,         &
                               asym_time,asym_dt_max,create_plot_file,           &
                               use_constd,constd_values,wdelt,consts_values,     &
                               use_consts,use_avg_glf_chi,vloopvb,use_constconv, &
                               constconv_values,nmlt0,set_chie_chii,use_avg_chi, &
                               nsteps_plot,plot_time_incr,save_incr_plot_file,   &
                               nsteps_restart,create_restart_file,prt_turbm,     &
                               restart_time_incr,save_incr_restart_file,         &
                               max_bwcycle,max_bandwidth,min_bandwidth,          &
                               cycle_bandwidth,eq_split, cgd,cgexp,crit_grad,    &
                               ts_smfactor,single_density_simulation,            &
                               use_hyperbolic_den_eq,comp_var,use_glf23_flux,    &
                               itran_max,use_tglf,select_solver,disspc,          &
                               reuse_tglf,delta_u_crit,conve_mult,convi_mult,    &
                               use_stab_flux,td_glf,ss_glf,td_tglf,ss_tglf,      &
                               td_neo,ss_neo,glf_flx_use,use_flow_eq,glf_on,     &
                               tglf_on,set_parabolic_profiles,angrcple,          &
                               it_imp_steps,core_iters_max,ambip_diffusion,      &
                               use_paleo,tol_corrector,rconsts_min,rconsts_max,  &
                               td_mmm,ss_mmm,use_mmm,use_mmm_flux,mmm_on,        &
                               mmm_flx_use,use_qrad_input,use_forcebal,          &
                               cyclotron_reflection,use_mmm_pinch,               &
                               use_mmm_vpol_input,use_mmm_vpar_input,            &
                               use_mmm_er_input,elc_eng_fixed_transpt,           &
                               ion_eng_fixed_transpt,                            &
                               ion_den_fixed_transpt,                            &
                               ion_torrot_fixed_transpt,                         &
                               curden_fixed_transpt,                             &
                               use_input_confinement,frz_tglf,frz_glf,frz_mmm,   &
                               use_frz_chi,fixed_profile_multiplier,             &
                               time_deriv_in_flux

    USE MPI_data,        ONLY : myid,master,parallel_model

    USE neo_tport,       ONLY : neocl_mult,jneo,use_forcebal_chie,use_forcebal_chii

    USE common_constants,ONLY : izero,zeroc

    USE grid_class,      ONLY : nj,nj_out,allow_regrid,reset_cap_parms,       &
                                set_cap1,nj_max,nj_start,use_compact_schemes, &
                                fit_grid,frac_core_points

    USE io_gcnmp,        ONLY : namelist_filename,nlog,ncrt,switch_iterdb_output,&
                                rpc,iterdb_12_input_file,fixt_filename,maxchr

    USE ions_gcnmp,      ONLY : nprim,nimp,nion,name_size,ni_sc

    USE curden_terms,    ONLY : jboot,rbsaxis,irfc,ibcur,currf,curbeam

    USE bc_values_gcnmp, ONLY : u_vloop_bc,mult_den,mult_flux,gammat_bc,    &
                                axis_bc_type,fix_edge_ni_bc,fix_edge_te_bc, &
                                fix_edge_ti_bc,fix_edge_rot_bc,robin_bc,    &
                                conserve_bc_flux,flux_mult_bc 

    USE glf23_gcnmp,     ONLY : exbmult_glf,x_alpha_glf,                    &
                                jroot_glf,limp_glf,ibtflag_glf,irotstab,    &
                                dv_method,spline_gsc,dv_delt,               &
                                itte_dv,itti_dv,itangrot_dv,itenp_dv,       &
                                itene_dv,jeigen_glf,i_delay,use_mask,       &
                                t_delay_glf23,glf_v_mult,glf_d_mult

    USE glf_tport_gcnmp, ONLY : glf_input_defaults

    USE source_terms_gcnmp, ONLY: beam_th_mult, stfuse_mult, sbfuse_mult,      &
                                  qrad_mult,qfuse_mult_e,qfuse_mult_i,         &
                                  qrf_mult_e,qrf_mult_i,qbeam_mult_e,          &
                                  qbeam_mult_i

    USE fusion_gcnmp,    ONLY : differential_burnup,internal_thermal_fusion,   &
                                alpha_slow

    USE tglfin,          ONLY :                                                & 
                            adi_elec_tg,adiabatic_elec_tg,                     & !aliases
                            find_width_tg,aiwt_tg,alpha_n_tg,alpha_t_tg,       &
                            nb_max_tg,nbasis_max_tg,                           & ! aliases  
                            nb_min_tg,nbasis_min_tg,                           & ! aliases  
                            nx_tg,nxgrid_tg,                                   & ! aliases  
                            nspecies_tg,new_eikonal_tg,ibranch_tg,             &
                            alpha_kx_e_tg,alpha_kx_p_tg,alpha_kx_n_tg,         &
                            alpha_kx_t_tg,                                     &
                            debye_fac_tg,debye_factor_tg,                      &  !aliases
                            nmodes_tg, iflux_tg,drmindx_tg,drmajdx_tg,         &
                            xky0_tg,ky_tg,                                     &  !aliases 
                            width_max_tg, width_min_tg, nwidth_tg,             &
                            park_tg, ghat_tg, gchat_tg, alpha_e_tg,alpha_dia,  &
                            gamma_e_tg, alpha_p_tg, gamma_p_tg, igeo_tg,       &
                            theta_trap_tg,theta_trapped_tg,                    &  !aliases 
                            cnorm_tg, theta0_tg,rmin_tg,                       &
                            rmaj_tg,use_bisection_tg,mass_tg,mass_ion_tg,      &
                            mass_imp_tg,zs_tg, taus_tg,j_surface_tg,           &
                            as_tg, rlns_tg, rlts_tg,q_tg, xnuei_tg, wd_zero_tg,&
                            betae_tg, shat_tg, alpha_tg,xwell_tg, kappa_tg,    &
                            s_kappa_tg, delta_tg, s_delta_tg, zeff_tg,         &
                            debye_tg, use_bper_tg, use_bpar_tg,q_prime_tg,     & 
                            p_prime_tg, filter_tg,ns_tg,kx0_tg,                &
                            Linsker_tg,Linsker_factor_tg,                      & ! alias
                            gradB_tg,gradB_factor_tg,                          & ! alias
                            nky_tg, b_model_tg, ft_model_tg, xnuei_fac_tg,     &
                            shift_tg,aimpwt_tg,vpar_tg,                        &
                            etg_fac_tg,etg_factor_tg,                          & ! alias
                            xnu_fac_tg,xnu_factor_tg,                          & ! alias
                            ialphastab,xalpha,test_tglf,                       &
                            tglf_v_mult,tglf_d_mult,kygrid_model_tg,           &
                            tglf_grid_point_parallel,tglf_wave_no_parallel,    &
                            xnu_model_tg,vexb_mach_tg,vpar_mach_tg ,           &
                            alpha_quench_tg,sat_rule_tg,alpha_mhd_tg,          &
                            rlne_tg,rlni_tg,rlnimp_tg,rlti_tg,xnu_tg,x_psi_tg, &
                            zion_tg,zimp_tg,tglf_flux_mult,tglf_implicit_grad, &
                            tglf_version_no,dzmajdx_tg,use_TM_tg,xnue_tg,      &
                            zmaj_tg,use_mhd_rule_tg,damp_psi_tg,damp_sig_tg,   &
                            vexb_shear_tg, vpar_shear_tg,vns_shear_tg,         &
                            vts_shear_tg, nfourier_tg,fourier_tg,vexb_tg,      &
                            vpar_shear_model_tg,vpar_model_tg,sign_Bt_tg,      &
                            sign_It_tg,zeta_tg,s_zeta_tg,tglf_rmin_sa_in,      &
                            tglf_q_sa_in,tglf_shat_sa_in,tglf_alpha_sa_in,     &
                            tglf_xwell_sa_in,tglf_theta0_sa_in,                &
                            tglf_b_model_sa_in,tglf_ft_model_sa_in,            &
                            sign_It_tg_state,sign_Bt_tg_state

    USE mmm_in,          ONLY :                                                &
                                mmm_v_mult,mmm_d_mult,mmm_cmodel,mmm_cW20,     &
                                mmm_cDBM,mmm_cETG,mmm_lW20, mmm_lDBM,          &
                                mmm_lETG,mmm_cswitch,mmm_lswitch,              &
                                BADREAL,BADINT,mmm_scfctr,use_mmm_zct

   USE FORCEBAL_DATA_MOD,ONLY :  l_banana,       & ! forcebal_data_mod.f90
                                 l_pfirsch,      &
                                 l_classical,    &
                                 l_potato,       &
                                 l_squeeze,      &
                                 l_reduce_out,   &
                                 k_order

   USE modmmm7_1,        ONLY : MAXNOPT
                            

    USE turb_tport,      ONLY : wturbd,wturbv,gmax_crit

    USE  vector_class,   ONLY : zero_vector

    IMPLICIT NONE 

    INTEGER(I2B)  ionml,ionml_temp
    INTEGER(I2B), EXTERNAL :: get_next_io_unit 
    INTEGER(I4B)  j,ienct,run_mhd
    INTEGER, PARAMETER :: nion_max = 10
    REAL(DP) fix_edge_nion_bc(nion_max)                          ! this is reqired because allocatable 
                                                                 ! arrays are not allowe in namelists
    REAL(DP) irfc_mult,ibcur_mult,nictr,niavg,nistd,nirho, &
             rhon_betan_ped,betan_ped,nispline
    INTEGER(I4B),DIMENSION(5) :: itenid,itenpd                   ! namelist input
    REAL(DP),DIMENSION(5)     :: mul_den,mul_flux,gamma_bc       ! namelist input
    REAL(DP), DIMENSION(50)   :: rhon_spline,ni_knot
    REAL(DP) pellet_rmin,pellet_dep_width,pellet_np,pellet_freq
!    CHARACTER(LEN = LEN(namelist_filename)+2) name_temp
    CHARACTER(LEN = maxchr+2) name_temp
    CHARACTER(LEN =132) jbootstrap
    CHARACTER(len = name_size),PUBLIC ::  pellet_name
    LOGICAL pellet_inject



    NAMELIST /run_data/                                                                       &
      allow_regrid    ,      alpha_slow      ,      ambip_diffusion ,      angrcple        ,  &
      asym_dt_max     ,      asym_time       ,      axis_bc_type    ,      bandwidth       ,  &
      beam_th_mult    ,      betan_ped       ,      cgd             ,      cgexp           ,  &
      comp_var        ,      conserve_bc_flux,                                                &     
      constconv_values,      constd_values   ,      consts_values   ,                         &
      conve_mult      ,      convi_mult      ,      conv_skip       ,      core_iters_max  ,  &
      create_plot_file,      create_restart_file,   crit_grad       ,      cycle_bandwidth ,  &
      curden_fixed_transpt,                                                                   &
      cyclotron_reflection,  dbg_print, differential_burnup,   disspc          ,      dtmax,  &
      dtmin           ,      dtstart         ,      dv_delt         ,      dv_method       ,  &
      elc_eng_fixed_transpt,                                                                  &
      eq_split        ,      exbmult_glf     ,      fdigits         ,      fit_grid        ,  &
      fix_edge_nion_bc,      fix_edge_rot_bc ,      fix_edge_te_bc  ,      fix_edge_ti_bc  ,  &
      fixed_profile_multiplier, fixt_filename,      flux_mult_bc,                             &
      frac_core_points,      freeze_type     ,      frz_tglf        ,      frz_glf         ,  &
      frz_mmm,               fvectol         ,      gamma_bc        ,                         &
      gradtol         ,      ibcur_mult      ,      ibtflag_glf ,   internal_thermal_fusion,  &
      ion_den_fixed_transpt, ion_eng_fixed_transpt,  ion_torrot_fixed_transpt,                &
      irfc_mult       ,      irotstab        ,      itangrot_dv     ,      itene_dv        ,  &
      itenid          ,      itenpd          ,      itenp_dv        ,iterdb_12_input_file  ,  &
      iters_freeze    ,      itte            ,      itte_dv         ,      itti            ,  &
      itti_dv         ,      itw             ,      itxj            ,      it_imp_steps    ,  &
      i_delay         ,      jacobian_type   ,      jac_skip        ,      jbootstrap      ,  &
      jeigen_glf      ,      jneo            ,      jroot_glf       ,      limp_glf        ,  &
      maxfev          ,      max_bandwidth   ,      max_bwcycle     ,      max_steps       ,  &
      min_bandwidth   ,      mul_den         ,      mul_flux        ,      neocl_mult      ,  &
      niavg           ,      nictr           ,      nirho           ,      nispline        ,  &
      nistd           ,      ni_knot         ,      nj_max          ,      nj_out          ,  &
      nj_start        ,      non_lin_method  ,      nsteps_plot     ,      nsteps_restart  ,  &
      outputvb        ,      parallel_model  ,      pellet_dep_width,      pellet_freq     ,  &
      pellet_inject   ,      pellet_name     ,      pellet_np       ,      pellet_rmin     ,  &
      plot_time_incr  ,      prt_turbm       ,      qbeam_mult_e    ,      qbeam_mult_i    ,  &
      qfuse_mult_e    ,      qfuse_mult_i    ,      qrad_mult       ,      qrf_mult_e      ,  &
      qrf_mult_i      ,      random_pert     ,      rbsaxis         ,      rconsts_max     ,  &
      rconsts_min     ,      restart_time_incr,     rhon_betan_ped  ,      rhon_spline     ,  &
      robin_bc        ,                                                                       &
      rpc             ,      run_mhd         ,   save_incr_plot_file,  save_incr_restart_file,&
      sbfuse_mult     ,      select_solver   ,      set_cap1        ,      set_chie_chii   ,  &
      set_parabolic_profiles,single_density_simulation, spline_gsc  ,      ssqrmin         ,  &
      steady_state    ,      steptol         ,      stfuse_mult     ,  switch_iterdb_output,  &
      switch_method   ,      tglfvb          ,      theta           ,      time0           ,  &
      time_deriv_in_flux,                                                                     &
      time_max        ,      tol_corrector   ,      tot_iters_max   ,      ts_smfactor     ,  &
      t_delay_glf23   ,      use_avg_chi     ,      use_avg_glf_chi ,   use_compact_schemes,  &
      use_constconv   ,      use_constd      ,      use_consts      ,      use_flow_eq     ,  &
      use_glf23       ,      use_glf23_flux  , use_hyperbolic_den_eq,      use_mask        ,  &
      use_input_confinement, use_forcebal_chie,use_forcebal_chii, use_frz_chi,                &
      use_paleo       ,      use_qrad_input  ,      use_stab_flux   ,      u_vloop_bc      ,  &
      vloopvb         ,      wdelt           ,      wrt_nwt         ,      x_alpha_glf     ,  &


!  tglf namelist input:
                        adi_elec_tg , adiabatic_elec_tg,                         & !aliases
                        aiwt_tg,aimpwt_tg,alpha_mhd_tg,                          &
                        alpha_dia   , alpha_e_tg  ,alpha_quench_tg, alpha_n_tg,  & 
                        alpha_p_tg, alpha_t_tg,alpha_kx_e_tg,alpha_kx_p_tg,      &
                        alpha_kx_n_tg,alpha_kx_t_tg,                             &
                        alpha_tg    , as_tg       , betae_tg    , b_model_tg  ,  &
                        cnorm_tg    , delta_u_crit,                              &
                        debye_fac_tg,debye_factor_tg,                            & ! aliases 
                        debye_tg    , delta_tg,    drmindx_tg,dzmajdx_tg,        &
                        etg_fac_tg  ,etg_factor_tg,                              & ! aliases 
                        filter_tg   , find_width_tg,fourier_tg,                  &
                        ft_model_tg , gamma_e_tg  , gamma_p_tg  , gchat_tg    ,  &
                        ghat_tg     , glf_d_mult  , glf_v_mult  , gmax_crit   ,  &
                        gradB_tg    , ialphastab  , ibranch_tg  , iflux_tg    ,  &
                        igeo_tg     , kappa_tg    , kygrid_model_tg,             &
                        Linsker_tg  , Linsker_factor_tg,                         &  
                        mass_tg     , mass_ion_tg, mass_imp_tg,                  &
                        nb_max_tg   , nb_min_tg   , new_eikonal_tg  , nky_tg  ,  &
                        nfourier_tg,kx0_tg,                                      &
                        nmodes_tg   , nspecies_tg , ns_tg, nwidth_tg   , nx_tg,  &
                        park_tg     , p_prime_tg  , q_prime_tg  , q_tg        ,  &
                        reuse_tglf  , rlne_tg,rlni_tg,rlnimp_tg,rlns_tg,         &
                        rlti_tg,      rlts_tg, rmaj_tg     , rmin_tg,            &
                        sat_rule_tg,  shat_tg     , shift_tg , s_delta_tg,       &
                        s_kappa_tg  ,s_zeta_tg, taus_tg     , test_tglf   ,      &
                        tglf_d_mult ,tglf_v_mult , tglf_wave_no_parallel,        &
                        tglf_grid_point_parallel  ,tglf_flux_mult,               &
                        tglf_implicit_grad,theta0_tg, tglf_version_no,           &
                        tglf_rmin_sa_in,                                         &
                        tglf_q_sa_in,tglf_shat_sa_in,tglf_alpha_sa_in,           &
                        tglf_xwell_sa_in,tglf_theta0_sa_in,                      &
                        tglf_b_model_sa_in,tglf_ft_model_sa_in,                  & 
                        theta_trap_tg,sign_Bt_tg, sign_It_tg,                    &
                        theta_trapped_tg,use_bisection_tg,                       &
                        use_bpar_tg , use_bper_tg , use_mhd_rule_tg,             &
                        use_tglf    ,                                            &
                        vexb_mach_tg, vexb_shear_tg, vns_shear_tg,vts_shear_tg,  &
                        vexb_tg, vpar_shear_tg,vpar_tg,vpar_shear_model_tg,      &
                        vpar_model_tg,                                           &
                        vpar_mach_tg,wd_zero_tg    ,                             &
                        width_max_tg, width_min_tg, wturbd      , wturbv      ,  &
                        xalpha      , xky0_tg,xnu_tg, xnuei_fac_tg, xnuei_tg  ,  &
                        xnue_tg,xnu_model_tg, x_psi_tg,xwell_tg, zeff_tg,        &
                        zeta_tg,zion_tg,zimp_tg,zs_tg,use_Tm_tg,zmaj_tg,         &

! multimode namelist input gcnmp driver related:
                        mmm_cmodel,mmm_cW20,mmm_cDBM, mmm_cETG,                  &
                        mmm_lW20, mmm_lDBM, mmm_lETG,                            &
                        use_mmm_pinch,mmm_scfctr,use_mmm_zct,                    &
                        use_mmm,use_mmm_flux,mmm_v_mult,mmm_d_mult,              &
                        use_mmm_vpol_input,use_mmm_vpar_input,                   &
                        use_mmm_er_input,                                        &

! forcebal/nclass namelist input:                                    
                        use_forcebal,l_banana,l_pfirsch,l_classical,             &
                        l_potato,l_squeeze,l_reduce_out,k_order

     CONTAINS  

        SUBROUTINE read_namelist
! ------------------------------------------------------------------------------
! SET SOME RUN DIRECTIVE DEFAULTS AND read namelist "run_data" 
! NOTE: This routine is called only by the master process
! ---------------------------------------------------------------------HSJ-05/19/05

#ifndef GCNMP
     USE string_util,             ONLY : to_upper_case1
#endif

    INTEGER(I4B)  glf_flx_ctr,mmm_flx_ctr
    CHARACTER(len = equil_type_size)  equil_type


    INTERFACE
#ifdef GCNMP
      SUBROUTINE to_upper_case(string)
      USE nrtype,            ONLY : I4B
      IMPLICIT NONE
      INTEGER(I4B) l
      CHARACTER*(*), INTENT (INOUT) :: string
      END SUBROUTINE to_upper_case
#endif
    END INTERFACE


    !set local namelist defaults

        nictr = -1._DP ;  niavg = -1._DP  ;   nistd = -1._DP ; nirho =-1. ; nispline = -1
        rhon_betan_ped = -1. ; betan_ped = -1._DP
        rhon_spline(:) = -1._DP ; ni_knot(:) = -1._DP
        switch_iterdb_output = 0
        run_mhd = 0  ; equil_type ='none' ! valid 'fixbd_var_grid'
        set_cap1 = izero  ; eq_split = izero
        itenpd(:) = izero  ; itenid(:) = izero 
        IF(SIZE(itenid) .LT. SIZE(iteni))THEN

        ENDIF
        IF(SIZE(itenpd) .LT. SIZE(itenp))THEN
           lerrno =iomaxerr + 5_I4B
           CALL terminate(lerrno,nlog)
        ENDIF

        IF((SIZE(mul_den) .LT. SIZE(mult_den)) .OR.      &
                     (SIZE(mul_flux) .LT. SIZE(mult_flux)) )THEN
           lerrno =iomaxerr + 9_I4B
           CALL terminate(lerrno,nlog)
        ENDIF

        mul_flux(:) = zeroc ; mul_den(:) = 1._DP ; gamma_bc(:) = zeroc

        !ibcur and irfc are set in read_iterdb
        ibcur_mult = FLOAT(ibcur)   ; irfc_mult = FLOAT(irfc)
        
        pellet_inject                 = .FALSE.
        pellet_name(1:name_size)      = ' '
        pellet_rmin                   = zeroc
        pellet_dep_width              = zeroc
        pellet_np                     = zeroc
        pellet_freq                   = zeroc
        frac_core_points              = 1.0_DP
        fit_grid                      = .FALSE.
  

    !---------------------------------------------------------------------------
    ! set global namelist defaults:
    ! note that fix_edge_ni_bc(),fix_edge_te_bc,fix_edge_ti_bc,fix_edge_rot_bc
    ! have default values that were read in the statefile. Hence no new 
    ! defaults are given here. If the values from the statefile 
    ! are changed we need to redistribute them to all processors.
    ! We need to handle fix_edge_nion_bc here because of allocatable
    ! array issue with namelists. Oher fix edge parms are scalars so
    ! no special treateent is required.
    !---------------------------------------------------------------------------
        IF(nion .le. nion_max)THEN
              fix_edge_nion_bc(1:nion) = fix_edge_ni_bc(1:nion)
        ELSE
22        FORMAT(" ERROR nio_max in sub read_namelist is too small ")
          WRITE(ncrt, 22)
          lerrno = 50_I4B 
          CALL terminate(lerrno,nlog)
        ENDIF

        CALL gcnmp_namelist_defaults
        CALL glf_input_defaults
        CALL tglf_input_defaults
        CALL mmm_input_defaults

        ionml      = get_next_io_unit()
        ionml_temp = get_next_io_unit()
 

    name_temp = ADJUSTL(namelist_filename(1:LEN_TRIM(namelist_filename)))
#ifdef GCNMP
    CALL to_upper_case(name_temp)
#else
    CALL to_upper_case1(name_temp) 
#endif
#ifndef DARWIN
    IF(name_temp(1:LEN_TRIM(name_temp)) .EQ.                &
         namelist_filename(1:LEN_TRIM(name_temp)))          &
         name_temp = name_temp(1:LEN_TRIM(name_temp))//'_1' 
#else
      name_temp = name_temp(1:LEN_TRIM(name_temp))//'_1' 
#endif

    !------------------------------------------------------------------------------
    ! process namelist input file
    !------------------------------------------------------------------------------
    write (*,'(a)') name_temp,namelist_filename
    OPEN (unit = ionml, file = namelist_filename, status = 'OLD')
    OPEN (unit = ionml_temp,file =name_temp,status ='UNKNOWN')
    CALL strip_comments (ionml, ionml_temp)
    CLOSE(unit = ionml,status = 'KEEP')
    REWIND(UNIT=ionml_temp)
    READ(ionml_temp,nml=run_data,END = 2) !  better to let it
                                          !  fail than redirect with ERR = 
    CLOSE(unit = ionml_temp,status = 'DELETE')
 

    !Overide some user inputs if necessary for use_input_confinement scheme
    !logic for fixed input confinement . 
    ! check user input:
    IF(ion_torrot_fixed_transpt)THEN
      WRITE(ncrt, FMT='("option ion_torrot_fixed_transpt not implemented" )')
      lerrno = 66
      CALL terminate(lerrno,nlog)
      ion_torrot_fixed_transpt      = .FALSE.   ! true not implemented  
    ENDIF  
    IF(ion_den_fixed_transpt)THEN
      WRITE(ncrt, FMT='("option ion_den_fixed_transpt not implemented" )')
      lerrno = 66
      CALL terminate(lerrno,nlog)
      ion_den_fixed_transpt   = .FALSE.   ! true not implemented
    ENDIF
    IF(curden_fixed_transpt)THEN
          WRITE(ncrt, FMT='("option  curden_fixed_transpt not implemented" )')
          lerrno = 66
          CALL terminate(lerrno,nlog)
          ! curden_fixed_transpt          = .FALSE.   ! true not implemented
    ENDIF
    IF(elc_eng_fixed_transpt)    use_input_confinement     = .TRUE.
    IF(ion_eng_fixed_transpt)    use_input_confinement     = .TRUE.
    IF(ion_torrot_fixed_transpt) use_input_confinement     = .TRUE.
    IF(ion_den_fixed_transpt )   use_input_confinement     = .TRUE.
    IF(curden_fixed_transpt)     use_input_confinement     = .TRUE.

 

    IF(use_input_confinement)THEN
       set_chie_chii = zeroc
       allow_regrid  = .FALSE.
!       select_solver = "full_newton"  ! 8888889999999 see time_advance.f90 for use of diffusive_solve
       use_flow_eq = 0
       !neocl_mult(:,:) = zeroc ! no neoclassical transport
       !jneo = 20 ! use this to indicate no neoclassical transport
       IF( .NOT. frz_glf)THEN
          use_glf23(:)   = izero
          use_glf23_flux = izero
       ENDIF
       IF(.NOT. frz_tglf)THEN
          use_tglf       = izero 
       ENDIF
       IF(.NOT. frz_mmm)THEN
          use_mmm        = izero
          use_mmm_flux   = izero
       ENDIF
       steady_state   = 1  ! may allow 0 later
       use_forcebal   = izero
       use_constd(:)  = izero
       use_consts(:)  = izero
       use_hyperbolic_den_eq = izero
       use_avg_chi    = izero
    ENDIF


     IF(.NOT. allow_regrid) nj_out = nj
    !fix_edge_nion_bc(1:nion) either contains a copy of fix_edge_ni_bc(1:nion)
    ! or new values if they were set in the namelist. So reset fix_edge_ni_bc(1:nion)
    fix_edge_ni_bc(1:nion) = fix_edge_nion_bc(1:nion)


    dt = dtstart
    nmlt0 = .FALSE.

    ni_sc(1)  = nictr ; ni_sc(2) = niavg ; ni_sc(3)= nistd ; ni_sc(4) = nirho
    mhd_dat%ni_spline = nispline
    mhd_dat%rhon_betan_ped = rhon_betan_ped
    mhd_dat%betan_ped      = betan_ped
    mhd_dat%run_mhd        = run_mhd
    mhd_dat%equil_type     = equil_type
    IF(mhd_dat%ni_spline .GT. 1)THEN
       mhd_dat%ni_knot   = zero_vector(mhd_dat%ni_spline )
       mhd_dat%rhon_spline = zero_vector(mhd_dat%ni_spline )
       DO j=1,mhd_dat%ni_spline
          mhd_dat%ni_knot%data(j)     = ni_knot(j)
          mhd_dat%rhon_spline%data(j)   = rhon_spline(j)
       ENDDO
    ELSE
       mhd_dat%ni_knot   = zero_vector(1)
       mhd_dat%rhon_spline = zero_vector(1)
    ENDIF
    IF(time0 .LT. zeroc)THEN
        time0 = time          !set initial time to iterdb file time
        WRITE(nlog,FMT='(" start time =" ,f12.6," taken from iterdb file")')time0
    ELSE
       WRITE(nlog,FMT='(" start time =" ,f12.6," taken from namelist file")')time0
       nmlt0 = .TRUE.
    ENDIF
    WRITE(nlog,FMT ='(" max number of time steps allowed = ",i12,/)')max_steps
    tGCNMs = time0

    IF(time_max .LT. zeroc .AND. steady_state == 1 )THEN
       time_max = tGCNMf !set final time to iterdb file time
       WRITE(nlog,FMT='(" end  time = ",f12.6," taken from iterdb file")')time_max
    ELSEIF(time_max .LT. zeroc .AND. steady_state == 0 )THEN
       WRITE(nlog,FMT='(" run into steady state")')
    ELSE
      IF(steady_state == 1)THEN
          WRITE(nlog,FMT='(" end time = ",f16.6," taken from namelist file")')time_max
      ELSE
          WRITE(nlog,FMT='(" run into steady state")')
      ENDIF
    ENDIF
    IF(steady_state == 0 )eq_split = izero !no time stepping so no equation splitting

    IF(time_max .LT. time0) lerrno = iomaxerr + 1_I4B
    IF(lerrno .GT. 0) CALL terminate(lerrno,nlog)
    tGCNMf = time_max
 
    jboot = -1
    IF( jbootstrap == 'Sauter')THEN
       jboot = 110                               ! used internally in gcnm
       !note jboot = 111 option (repalces eta with Sauter eta) is not implemented
    ELSEIF( jbootstrap == 'Sauter_zeff')THEN
       jboot = 112 ! use zeff instead of z of primary ion
    ELSE
       WRITE(ncrt,FMT='("jbootstrap = ",a," not recognized")')jbootstrap
       lerrno = iomaxerr + 8_I4B 
       CALL terminate(lerrno,nlog)
    ENDIF

    IF(ABS(neocl_mult(4,4) ) .LT. TINY(neocl_mult(4,4)))THEN
       WRITE(ncrt,15)
       WRITE(nlog,15)
15     FORMAT(2x,'Caution: neocl_mult(4,4) = 0.0 causes resitivity and related',/ &
           2x,'quantities to be zero')
    ENDIF

    itenp(1:nprim) = itenpd(1:nprim) !because allocatable arrays 
                                     !are not allowed in namelists
    iteni(1:nimp)  = itenid(1:nimp)
    IF(set_cap1 .GT. 0 ) CALL reset_cap_parms(set_cap1)

    IF(ABS(irfc_mult -FLOAT( irfc)) .GT. 1.e-5 .AND. irfc .NE. 0)THEN
       IF(ABS(irfc_mult) .GT. 1.e-8)THEN
          currf(:) = currf(:)*irfc_mult
       ELSE
          currf(:) = 0.0_DP
          irfc = 0
       ENDIF

    ENDIF

    DO j =1,nion
       IF(mul_den(j) == 0.0 .AND. mul_flux(j) == 0.0)THEN
          WRITE(nlog, 20)j
20        FORMAT(" ERROR ambiguous input of mul_den and  mul_flux ", &
                 " For species k= " ,i5,/," cant both be zero")
          WRITE(ncrt, 20)j
          lerrno = iomaxerr + 21_I4B 
          CALL terminate(lerrno,nlog)
       ENDIF
    ENDDO

    mult_den(1:nion) = mul_den(1:nion) ; mult_flux(1:nion) =  mul_flux(1:nion)

    gammat_bc(1:nion) = gamma_bc(1:nion)
      
    jac_skip = MAX(jac_skip,1_I4B)

    IF(ABS(ibcur_mult -FLOAT( ibcur)) .GT. 1.e-5 .AND. ibcur .NE. 0) THEN
       IF(ABS(ibcur_mult) .GT. 1.e-8)THEN
          curbeam(:) = curbeam(:)*ibcur_mult
       ELSE
          curbeam(:) = 0.0_DP
          ibcur = 0
       ENDIF
    ENDIF

    IF(steady_state == 0)use_mask =0      ! make sure multipliers of
                                          ! transport coefficients,eg cm_xke,etc
                                          ! remain at  1.0 for steady state runs.

    IF(jacobian_type .EQ. 1 .AND. parallel_model == 1)THEN
       lerrno = iomaxerr + 103
       CALL terminate(lerrno,nlog)
    ENDIF


! ------------------------------------------------------------------
!   check logic for density_simulation. Turn it off if 
!   densities are not evolved:
! ------------------------------------------------------------------
    ienct = izero
    DO j=1,nprim
       IF( itenp(j) .NE. 0) ienct = ienct+1
    ENDDO
    DO j=1,nimp
       IF(iteni(j) .NE. 0)  ienct = ienct +1
    ENDDO


! ------------------------------------------------------------------
!   check logic for analysis/simulation mode  Turn confinement models
!   off if corresponding profile is not evolved:
! ------------------------------------------------------------------
    IF(use_input_confinement)THEN
       use_tglf   = 0 
       use_glf23  = 0   
       use_mmm    =0
    ENDIF

    ! eventually turbulent transpor in one channel and fixed transport in
    ! another may be allowed. Then the follwoing logic needs to be activated
    ! after eliminating the above zeroing of all channels
    IF(ienct == izero .OR. (ienct ==1 .AND. ion_den_fixed_transpt ))THEN
       single_density_simulation  = .FALSE.
       use_tglf(1)   = 0                         ! no density simulations so these are off
       use_glf23(1)  = 0                         ! (even if on in namelist)
       use_mmm(1)    = 0
    ENDIF

    IF(itte == 0 .OR. (itte ==1 .AND. elc_eng_fixed_transpt) )THEN
       use_tglf(2)   = 0                         ! no turbulent electron temperature simulation 
                                                 ! so these are off (even if on in namelist)
                                                 ! fixed transport overides turbulent!!!
       use_glf23(2)  = 0                         
       use_mmm(2)    = 0
    ENDIF

    IF(itti == 0 )THEN
       use_tglf(3)   = 0                         ! no ion temperature  simulation so these are off
       use_glf23(3)  = 0                         ! (even if on in namelist)
       use_mmm(3)    = 0
    ENDIF
    IF(itw == 0)THEN
       use_tglf(5)   = 0                         ! no toroidal rotation  simulation so these are off
       use_glf23(5)  = 0                         ! (even if on in namelist)
       use_mmm(5)    = 0
    ENDIF

!---------------------------------------------------------------------------------------------
! -- sort  out turbulent transport specifications
!---------------------------------------------------------------------------------------------

    td_glf       = .FALSE.                             ! time dependent glf23
    ss_glf       = .FALSE.                             ! steady state glf23

    td_tglf      = .FALSE.                             ! time dependent glf23
    ss_tglf      = .FALSE.                             ! steady state glf23

    td_neo       = .FALSE.                             ! time dependent neo
    ss_neo       = .FALSE.                             ! steady state neo

    td_mmm       = .FALSE.
    ss_mmm       = .FALSE.

    tglf_on      = .FALSE.
    glf_on       = .FALSE.
    mmm_on       = .FALSE.

    glf_flx_use  = .FALSE.
    mmm_flx_use  = .FALSE.





! -------------------------------------------------------------------------------------------
! -- check logic for tglf model usage. tglf takes precedence over glf so
! -- it is checked first. NOTE include_tglf is set in SUBROUTINE set_itran
! -- use_tglf and use_glf select confinement models but do not
! -- control analysis/simulation (which is done with switches itenp,iteni,itte,itti,itxj,itw)
! --------------------------------------------------------------------------------------------



    DO j=1,itran_max
       IF(use_tglf(j)  .NE.  izero )THEN        ! tglf is on for this  dep variable
          tglf_on = .TRUE.
          IF( j == 4)THEN
             lerrno = 11
             CALL terminate(lerrno,nlog)
          ENDIF
          IF( j == 5)THEN                       ! tglf v1.93 allows this
             !lerrno = 12 
             !CALL terminate(lerrno,nlog)
          ENDIF
       ENDIF
       IF(j .NE. 4)THEN                         ! 4 ==> poloidal B field
          if(use_glf23(j) .NE. izero)glf_on = .TRUE.
          IF(use_mmm(j)   .NE. izero)mmm_on = .TRUE.
       ENDIF
    ENDDO
    IF( (glf_on .AND. tglf_on)  .OR. (glf_on .AND. mmm_on)    &
         .OR.  (tglf_on .AND. mmm_on))THEN   ! only one turbulent  model can  be on
       lerrno = 23
       CALL terminate(lerrno,nlog)
    ENDIF
    IF(tglf_on)THEN                             ! tglf takes precedence
       use_glf23(:)      = izero                ! make sure glf23,mmm are off for this run.
       use_glf23_flux(:) = izero
       use_mmm(:)        = izero
       use_mmm_flux(:)   = izero
       glf_on = .FALSE.            
       mmm_on = .FALSE.
       if(use_tglf(1) > 0 .AND. itenpd(1)  == 0) lerrno = 15
       if(use_tglf(2) > 0 .AND. itte       == 0) lerrno = 16
       if(use_tglf(3) > 0 .AND. itti       == 0) lerrno = 17
       if(use_tglf(5) > 0 .AND. itw        == 0) lerrno = 18
       IF(lerrno .NE. izero)THEN
          CALL terminate(lerrno,nlog)
       ENDIF
    ENDIF


! -------------------------------------------------------------------
! -- check logic for glf23 flux model usage
! -- NOTE include_glf is set in SUBROUTINE set_itran
! -------------------------------------------------------------------
     lerrno = izero
     IF(glf_on)THEN
          if(use_glf23(1) > 0 .AND. itenpd(1) == 0)lerrno = 15
          if(use_glf23(2) > 0 .AND. itte      == 0)lerrno = 16
          if(use_glf23(3) > 0 .AND. itti      == 0)lerrno = 17
          if(use_glf23(5) > 0 .AND. itw       == 0)lerrno = 18
          IF(lerrno .NE. izero)THEN
             CALL terminate(lerrno,nlog)
          ENDIF
     ENDIF
    glf_flx_ctr = izero
    DO j=1,itran_max
       if(use_glf23_flux(j) .NE. izero) glf_flx_ctr = glf_flx_ctr + 1
       IF(use_glf23(j)  == izero )THEN                ! glf23 will not be used
          IF(use_glf23_flux(j) .NE. izero)THEN        ! glf model and flux go together
             lerrno = j
          ENDIF
       ENDIF
    ENDDO
    IF(lerrno .NE. izero)THEN
        IF(lerrno > 0) lerrno = 8_I4B
       ! IF(lerrno < 0) lerrno = 9_I4B removed 2/29/2012 HSJ
        CALL terminate(lerrno,nlog)
    ENDIF
! -------------------------------------------------------------------
! -- check logic for mmm  flux model usage
! -- NOTE include_glf is set in SUBROUTINE set_itran
! -------------------------------------------------------------------
     lerrno = izero
     IF(mmm_on)THEN
          if(use_mmm(1) > izero .AND. itenpd(1) == 0)lerrno = 15
          if(use_mmm(2) > izero .AND. itte      == 0)lerrno = 16
          if(use_mmm(3) > izero .AND. itti      == 0)lerrno = 17
          if(use_mmm(5) > izero .AND. itw       == 0)lerrno = 18
          IF(lerrno .NE. izero)THEN
             CALL terminate(lerrno,nlog)
          ENDIF
     ENDIF
    mmm_flx_ctr = izero
    DO j=1,itran_max
       if(use_mmm_flux(j) .NE. izero) mmm_flx_ctr = mmm_flx_ctr + 1
       IF(use_mmm(j)  == izero )THEN                ! mmm will not be used
          IF(use_mmm_flux(j) .NE. izero)THEN        ! mmm  model and flux go together
             lerrno = j
          ENDIF
       ELSE
          use_forcebal = 1   ! turn on forceball if any use_mmm  .ne. 0
       ENDIF
    ENDDO
    IF(use_forcebal .NE. 1)THEN
       use_forcebal_chie = .FALSE.
       use_forcebal_chii = .FALSE.
    ENDIF
   
    IF(lerrno .NE. izero)THEN
        lerrno = 9_I4B                              ! new hsj 2/29/2012
        CALL terminate(lerrno,nlog)
    ENDIF


!----------------------------------------------------------------------------------
! -- for turbulent transport the solver may be  forced by the confinment model
! -- for both time dependent and time independent cases:
!-----------------------------------------------------------------------------------   

    IF(glf_on .AND. glf_flx_ctr .GT. izero .AND. steady_state ==0  )select_solver='full_newton'
    IF(glf_on .AND. glf_flx_ctr  ==  izero .AND. steady_state ==0  )select_solver='full_newton'

    IF(mmm_on .AND. mmm_flx_ctr .GT. izero .AND. steady_state ==0  )select_solver='full_newton'
    IF(mmm_on .AND. mmm_flx_ctr  ==  izero .AND. steady_state ==0  )select_solver='full_newton'
    IF(tglf_on .AND. steady_state ==1  )THEN
       IF(select_solver .ne. 'maccormack')   select_solver='nwt_pred_cor'
       IF(select_solver .ne. 'nwt_pred_cor') select_solver='maccormack'
    ENDIF

    IF(tglf_on .AND. steady_state ==0  )select_solver='full_newton'

    IF(select_solver == 'none')THEN
       ! no turbulent confinement selected. Solver selection is not
       ! forced by confinement model but user has not specified which
       ! solver to use
          lerrno = 14
          CALL terminate(lerrno,nlog)
    ENDIF

    IF(select_solver == 'full_newton' .AND. .NOT. glf_on .AND. .NOT. tglf_on .AND. .NOT. mmm_on)THEN
          IF(use_flow_eq == -1)THEN
             ! Newton method for time dependent or time independent  non turbulent cases
             ! use of diffusive or flow  equations must be explicitely
             ! given in this case
             !lerrno = 39
             !CALL terminate(lerrno,nlog)
          ENDIF
    ELSEIF(select_solver == 'maccormack' .AND. .NOT. glf_on .AND. .NOT. tglf_on .AND. .not. mmm_on)THEN
       IF(steady_state == 0)THEN
          !macormack  method for time independent non turbulent cases does not exist
          lerrno = 38
          CALL terminate(lerrno,nlog)
       ELSE
          !macormack  method for time dependent non turbulent cases
          use_flow_eq = 1
!          CALL terminate(lerrno,nlog)
       ENDIF
    ENDIF
!-------------------------------------------------------------------------------   
! -- if solver is maccormack and glf23 is selected make sure flux option is set.
! -- (since tglf does not have a diffusive  solver we dont need to check for it).
! -- Maccormack only supports time dependent solutions!
!-------------------------------------------------------------------------------   

    IF( select_solver == 'maccormack')THEN
       reuse_tglf = MAX(reuse_tglf,1)            ! reuse_tglf= 0 causes crash in modulo function
       IF(steady_state .NE. 1)                               lerrno = 24
       IF(glf_on )THEN
          td_glf = .TRUE.     
       ELSEIF(tglf_on)THEN
          td_tglf = .TRUE.  
       ELSEIF(mmm_on)THEN
          td_mmm  = .TRUE.
       ENDIF
       td_neo = .TRUE.
       IF(lerrno .NE. izero) CALL terminate(lerrno,nlog)
    ENDIF


    IF( select_solver == 'full_newton')THEN
       IF(steady_state == 1)THEN              ! time dependent case
                                              ! here only diffusive glf
                                              ! or diffusive mmm
                                              ! or neoclassical  is valid:
          IF( glf_flx_ctr .NE. izero .OR. mmm_flx_ctr .NE. izero  .OR. tglf_on)THEN
             lerrno = 36
             CALL terminate(lerrno,nlog)
          ENDIF
          IF(glf_on)td_glf = .TRUE. 
          IF(mmm_on)td_mmm = .TRUE. 
          td_neo = .TRUE.
       ELSE                          
         ! time indepenent case
         ! here we have diffusive glf or flow based glf
         ! or flow based tglf
         ! or diffuse mmm ! or flow based mmm
          ss_neo = .TRUE.
          IF(glf_on)THEN
              ss_glf = .TRUE.                  
          ELSEIF(tglf_on)THEN
             ss_tglf = .TRUE.
          ELSEIF(mmm_on)THEN
             ss_mmm = .TRUE.             
          ENDIF
       ENDIF

    ENDIF

    IF(select_solver == 'nwt_pred_cor')THEN      ! only time dependent case
       IF(glf_on)td_glf   = .TRUE.
       td_neo             = .TRUE.
       IF(tglf_on)td_tglf = .TRUE.
       IF(mmm_on)td_mmm   = .TRUE.
    ENDIF

    IF(td_tglf)THEN
       IF(tglf_wave_no_parallel .AND. .NOT. tglf_grid_point_parallel)THEN
          ! tglf_wave_no_parallel can be on only if tglf_grid_point_parallel
          ! is on
          lerrno = 58
          CALL terminate(lerrno,nlog)
       ENDIF
    ENDIF

    IF(glf_flx_ctr > 0)glf_flx_use = .TRUE.
    IF(mmm_flx_ctr > 0)mmm_flx_use = .TRUE.

    IF((jneo .lt. 0 .OR. jneo .GT. 7 ) .AND. jneo .NE. 20)THEN
       lerrno = 37
       CALL terminate(lerrno,nlog)
    ENDIF





!---------------------------------------------------------------------------
! -- check multimode input
! -- Assign user specified parameters if they  are set. Otherwise use
! -- the default values set in mmm_input_defaults
! -- note Weiland,DRBM,ETG  weights are on by default,
! -- see mmm_input_defaults.f90  mmm_cmodel = (1.,1.,1.)
! -- namelist can only turn these off
!---------------------------------------------------------------------------
    IF(mmm_on)THEN

       CALL mmm_check_input( mmm_cmodel(1), "mmm_cmodel" )
       CALL mmm_check_input( mmm_cmodel(2), "mmm_cmodel" )
       CALL mmm_check_input( mmm_cmodel(3), "mmm_cmodel" )


       IF ( ABS( mmm_cW20(1) - BADREAL ) > 1E-6_DP ) mmm_cswitch(1,1)= mmm_cW20(1)
       IF ( ABS( mmm_cW20(2) - BADREAL ) > 1E-6_DP ) mmm_cswitch(2,1)= mmm_cW20(2)
       IF ( ABS( mmm_cW20(3) - BADREAL ) > 1E-6_DP ) mmm_cswitch(3,1)= mmm_cW20(3)
       IF ( ABS( mmm_cW20(4) - BADREAL ) > 1E-6_DP ) mmm_cswitch(4,1)= mmm_cW20(4)
       IF ( ABS( mmm_cW20(5) - BADREAL ) > 1E-6_DP ) mmm_cswitch(5,1)= mmm_cW20(5)
       IF ( ABS( mmm_cW20(6) - BADREAL ) > 1E-6_DP ) mmm_cswitch(6,1)= mmm_cW20(6)


       IF ( ABS( mmm_cDBM(1) - BADREAL ) > 1E-6_DP ) mmm_cswitch(1,2)= mmm_cDBM(1)


       IF ( ABS( mmm_cETG(1) - BADREAL ) > 1E-6_DP ) mmm_cswitch(1,3)= mmm_cETG(1)
       IF ( ABS( mmm_cETG(2) - BADREAL ) > 1E-6_DP ) mmm_cswitch(2,3)= mmm_cETG(2)


       IF ( mmm_lW20(1) /= BADINT ) mmm_lswitch(1:MAXNOPT,1)= mmm_lW20

       IF ( mmm_lDBM(1) /= BADINT ) mmm_lswitch(1:MAXNOPT,2)= mmm_lDBM

       IF ( mmm_lETG(1) /= BADINT ) mmm_lswitch(1:MAXNOPT,3)= mmm_lETG


    ENDIF




! ---------------------------------------------------------------------
! -- Load pellet info into type pellet_params  from local namelist data:
! -- if this option is used (eg pellet_inject = T)
! -- then the pellet source in iterdb file will be
! -- ignored !!! 
! ----------------------------------------------------------------------

       pellet%inject                 = pellet_inject
       pellet%np                     = pellet_np
       pellet%name(1:name_size)      = pellet_name(1:name_size)
       pellet%dep_width              = pellet_dep_width
       pellet%rmin                   = pellet_rmin
       pellet%freq                   = pellet_freq

!-------------------------------------------------------------------------
! -- take care of some tglf aliases, v1.94 has new names for some things
! -- make sure we pick them up in a compatible manner here
!--------------------------------------------------------------------------
    IF(adi_elec_tg .EQV. .TRUE.)THEN  ! assumes both are defaulted to false
       adiabatic_elec_tg = .TRUE.
    ELSEIF(adiabatic_elec_tg .EQV. .TRUE.)THEN
       adi_elec_tg = .TRUE.
    ENDIF

    IF(nbasis_max_tg .NE. 4)THEN       ! assumes both are defaulted to 4
       !This means user set it in namelist so use it
       nb_max_tg  = nbasis_max_tg 
    ELSEIF(nb_max_tg .NE.  4)THEN
       nbasis_max_tg = nb_max_tg
    ENDIF

    IF(nbasis_min_tg .NE. 1)THEN       ! assumes both are defaulted to 1
       !This means user set it in namelist so use it
       nb_min_tg  = nbasis_min_tg 
    ELSEIF(nb_min_tg .NE.  1)THEN
       nbasis_min_tg = nb_min_tg
    ENDIF

    IF(nxgrid_tg .NE. 16)THEN       ! assumes both are defaulted to 16
       !This means user set it in namelist so use it
       nx_tg  = nxgrid_tg 
    ELSEIF(nx_tg .NE.  16)THEN
       nxgrid_tg = nx_tg
    ENDIF

    IF(ABS(debye_factor_tg -1._DP) .GT. 1.E-10_DP )THEN       ! assumes both are defaulted to 16
       !This means user set it in namelist so use it
       debye_fac_tg = debye_factor_tg
    ELSEIF(ABS(debye_fac_tg -1._DP) .GT. 1.E-10_DP )THEN
       debye_factor_tg = debye_fac_tg
    ENDIF

    IF(ABS(ky_tg - 0.3_DP) .GT. 1.E-10_DP)THEN
       xky0_tg = ky_tg
    ELSEIF(ABS(xky0_tg - 0.3_DP) .GT. 1.E-10_DP)THEN
       ky_tg = xky0_tg
    ENDIF

   IF(ABS(theta_trapped_tg - 0.7_DP) .GT. 1.E-10_DP)THEN
       theta_trap_tg = theta_trapped_tg 
    ELSEIF(ABS(theta_trap_tg  - 0.7_DP) .GT. 1.E-10_DP)THEN
       theta_trapped_tg = theta_trap_tg 
    ENDIF

   IF(ABS(Linsker_factor_tg - zeroc ) .GT. 1.E-10_DP)THEN
       Linsker_tg =  Linsker_factor_tg
    ELSEIF(ABS(Linsker_tg  - zeroc ) .GT. 1.E-10_DP)THEN
       Linsker_factor_tg = Linsker_tg
    ENDIF

   IF(ABS(gradB_factor_tg - zeroc ) .GT. 1.E-10_DP)THEN
       gradB_tg =  gradB_factor_tg
    ELSEIF(ABS(gradB_tg  - zeroc ) .GT. 1.E-10_DP)THEN
       gradB_factor_tg = gradB_tg
    ENDIF

   IF(ABS(etg_factor_tg - 4.0_DP ) .GT. 1.E-10_DP)THEN
       etg_fac_tg =  etg_factor_tg
    ELSEIF(ABS(etg_fac_tg  - 4.0_DP ) .GT. 1.E-10_DP)THEN
       etg_factor_tg = etg_fac_tg
    ENDIF


   IF(ABS(xnu_factor_tg - 1.0_DP ) .GT. 1.E-10_DP)THEN
       xnu_fac_tg =  xnu_factor_tg
    ELSEIF(ABS(xnu_fac_tg  - 1.0_DP ) .GT. 1.E-10_DP)THEN
       xnu_factor_tg = xnu_fac_tg
    ENDIF


   IF(nspecies_tg .NE. 2)THEN
       ns_tg = nspecies_tg
    ELSEIF(ns_tg .NE. 2)THEN
       nspecies_tg = ns_tg
    ENDIF
    IF(nspecies_tg .LT. 2)THEN
       lerrno = iomaxerr + 4_I4B
       CALL terminate(lerrno,nlog)
    ENDIF

   IF(ABS(tglf_rmin_sa_in -0.5_DP).GT. 1.E-10_DP)THEN
      rmin_tg  = tglf_rmin_sa_in
   ELSEIF(ABS(rmin_tg -0.5_DP).GT. 1.E-10_DP)THEN
      tglf_rmin_sa_in = rmin_tg
   ENDIF

   IF(ABS(tglf_q_sa_in -2.0_DP).GT. 1.E-10_DP)THEN
      q_tg  = tglf_q_sa_in
   ELSEIF(ABS(q_tg -2.0_DP).GT. 1.E-10_DP)THEN
      tglf_q_sa_in = q_tg
   ENDIF

   IF(ABS(tglf_shat_sa_in -1.0_DP).GT. 1.E-10_DP)THEN
      shat_tg  = tglf_shat_sa_in
   ELSEIF(ABS(shat_tg -1.0_DP).GT. 1.E-10_DP)THEN
      tglf_shat_sa_in = shat_tg
   ENDIF

   IF(ABS(tglf_alpha_sa_in).GT. 1.E-10_DP)THEN
      alpha_tg  = tglf_alpha_sa_in
   ELSEIF(ABS(alpha_tg).GT. 1.E-10_DP)THEN
      tglf_alpha_sa_in = alpha_tg
   ENDIF

   IF(ABS(tglf_xwell_sa_in).GT. 1.E-10_DP)THEN
      xwell_tg  = tglf_xwell_sa_in
   ELSEIF(ABS(xwell_tg).GT. 1.E-10_DP)THEN
      tglf_xwell_sa_in = xwell_tg
   ENDIF

   IF(ABS(tglf_theta0_sa_in) .NE. 1)THEN
      theta0_tg  = tglf_theta0_sa_in
   ELSEIF(ABS(theta0_tg).GT. 1.E-10_DP)THEN
      tglf_theta0_sa_in = theta0_tg
   ENDIF

   IF(tglf_b_model_sa_in  .NE. 1)THEN
      b_model_tg  = tglf_b_model_sa_in
   ELSEIF(b_model_tg .NE. 1)THEN
      tglf_b_model_sa_in = b_model_tg
   ENDIF


   !here is where the sign on bt and ip is set for the run
   ! do not change these values further. 
   ! They will be distributed in distribute namelist
   ! It is assumed that the state file was read
   ! If set in namelist then state file values are NOT used !
   IF( sign_Bt_tg .LT. -1.e5 )THEN
     !  sign_Bt_tg has the default value of -1.e6 so ign_Bt_tg_state,
     !  ign_It_tg_state are available. 
     !  it was not set in namelist. This means we take the
     ! statefile value:
      sign_Bt_tg = sign_Bt_tg_state
   ENDIF
   IF( sign_It_tg .LT. -1.e5 )THEN 
     !  sign_It_tg has the default value of -1.e6 so 
     !  it was not set in namelist. This means we take the
     !  statefile value:
      sign_It_tg = sign_It_tg_state
   ENDIF
  RETURN
1   lerrno = 4_I4B
    CALL terminate(lerrno,nlog)
  RETURN
2   lerrno = iomaxerr + 2_I4B
    CALL terminate(lerrno,nlog)
  RETURN
3   lerrno = iomaxerr + 3_I4B
    CALL terminate(lerrno,nlog)


  RETURN


  END SUBROUTINE read_namelist


  END MODULE GCNMP_namelist
