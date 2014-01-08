     MODULE SOLCON_GCNMP

      USE nrtype,                             ONLY : DP,I4B,I2B
      USE common_constants,                   ONLY : izero
      IMPLICIT NONE
!
      SAVE 
      INTEGER(I4B),PARAMETER,PUBLIC :: itran_max = 5
      INTEGER(I4B),PARAMETER,PUBLIC :: std_file_size = 256
      INTEGER(I4B),PARAMETER,PUBLIC :: std_char_size = 8
      INTEGER(I4B), PUBLIC :: prtcl_transpt,enrgy_transpt,ic_step,ip_step
      CHARACTER*8,PUBLIC   ::    istep
      REAL(DP),PUBLIC      :: q0_max,q0_radius,q0_mult,set_chie_chii
      REAL(DP),PUBLIC      :: cur_seed_amp,cur_seed_ctr,cur_seed_width
      REAL(DP),PUBLIC      :: time0, time, eqtime,timnew, time_max, dt,dtstart
      REAL(DP),PUBLIC      :: dtmin, dtmax, dtevmin, dtmult,relmin, relmax, relit
      REAL(DP),PUBLIC      :: theta, dvres, timzen, erneumax,ts_smfactor
      REAL(DP),PUBLIC      :: eqtim0, relaxsol, timav, dtsumi,cximult,dbgsvx_rcon
      REAL(DP),PUBLIC      :: timecrit, dtmaxcrit,cparam,dcparam_in,tol_corrector 
      REAL(DP),PUBLIC      :: cpu_time,ssqr_global,time_p,wdelt,delta_u_crit
      REAL(DP),PUBLIC      :: time_tol, tGCNMf,tGCNMs,asym_time,asym_dt_max
      REAL(DP),PUBLIC      :: conve_mult,convi_mult,te_mult_conve,ti_mult_convi
      REAL(DP),PUBLIC      :: cyclotron_reflection,fixed_profile_multiplier
      REAL(DP),PUBLIC,DIMENSION(:,:),ALLOCATABLE :: u,usave,udisp,dudr,dudt,dudr_save,dgudu, &
                                                    u_prev_it
                                                    ! see also dudtsv in source_terms_gcnmp
      REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE :: xscale,u_scale,edge_flux,edge_flux_save
      REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE :: cm_xke,cm_xki,cm_xkw,Mu,Ms,ssqr_local
      REAL(DP),PUBLIC,DIMENSION(itran_max) :: constd_values, consts_values, &
                                      constconv_values

      REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE ::     t_slave_compt_tot,               &
                                                      t_slave_compt_alltot             
                                               ! this array accumulates the entries
                                               ! in t_slave_compt_tot for each processor
                                               ! OVER THE ENTIRE RUN. (t_slave_compt_tot
                                               ! just gives the processor time for a single time step.)

      REAL(DP),PUBLIC,DIMENSION(:),ALLOCATABLE :: msg_monitr
!      REAL(DP),PARAMETER :: angrcple = 1.0_DP  
      REAL(DP)                  angrcple

      REAL(DP), DIMENSION(itran_max) :: cgd,cgexp,crit_grad,disspc, &
                                        rconsts_min,rconsts_max
!
      INTEGER(I4B) steady_state,eq_split,use_hyperbolic_den_eq
      INTEGER(I4B),DIMENSION(itran_max) :: use_constd,use_glf23,use_consts, &
                                   use_constconv,comp_var,use_glf23_flux,   &
                                   use_tglf, use_paleo,use_mmm,use_mmm_flux
      INTEGER(I4B),PUBLIC :: include_constd,include_constconv,include_consts, &
                             include_glf,include_neocl ,use_avg_chi,          &
                             increase_dt,include_tglf,reuse_tglf,             &
                             set_parabolic_profiles,ssqr_inc,core_iters_max,  &
                             it_imp_steps,include_mmm,use_qrad_input,         &
                             use_mmm_vpol_input,use_mmm_vpar_input,           &
                             use_forcebal,use_mmm_pinch,use_mmm_er_input
      INTEGER(I4B),PUBLIC,DIMENSION(:),ALLOCATABLE :: itran,itran_global,itenp,iteni, &
                                                      freeze_ni
      INTEGER(I4B),PUBLIC :: itte,itti,itxj,itw,ntot,nkt,nequations,nu,nkt_global
      INTEGER(I4B),PUBLIC :: max_steps, nsteps, itmax, irfsw1, isecrem0
      INTEGER(I4B),PUBLIC :: istop, icallneu, icallnub, icallrf,icalln,icall_tglf_ss_fdjac
      INTEGER(I4B),PUBLIC :: ilimdt, idterr, ifreya_old, use_avg_glf_chi
      INTEGER(I4B),PUBLIC :: ineucg, ifreya, freeze_type,switch_method,use_flow_eq
      INTEGER(I4B),PUBLIC :: glf_debug,external_beam_cur,conv_mntr,ihalve_nl
      LOGICAL,     PUBLIC ::  te_range_check,write_profiles,nmlt0,last_chance
      LOGICAL,     PUBLIC :: ipturbmax_quit,first_plot,first_restart_file
      LOGICAL,     PUBLIC :: cycle_bandwidth,single_density_simulation
      LOGICAL,     PUBLIC :: fd_switch,fd_pred,fd_corct,refactor_called,use_stab_flux
      LOGICAL,     PUBLIC :: td_glf,ss_glf,td_tglf,ss_tglf, td_neo,ss_neo,glf_flx_use, &
                             glf_on,tglf_on,norm_known,reset_source,tglf_ss_jac,       &
                             smooth_mhd_input,ambip_diffusion,td_mmm,ss_mmm,           &
                             mmm_on,mmm_flx_use

      LOGICAL,     PUBLIC :: use_fixed_elc_eng_transpt,use_fixed_ion_eng_transpt,  &
                             use_fixed_ion_den_transpt, use_fixed_torrot_transpt,  &
                             use_fixed_curden_transpt,use_input_confinement,       &
                             elc_eng_fixed_transpt,                                &
                             ion_eng_fixed_transpt,                                &
                             ion_den_fixed_transpt,                                &
                             ion_torrot_fixed_transpt,                             &
                             curden_fixed_transpt,                                 &
                             time_deriv_in_flux

      LOGICAL,     PUBLIC :: frz_glf,frz_tglf,frz_mmm,use_frz_chi

      INTEGER(I4B),PARAMETER,PUBLIC :: iangrot = 1
      INTEGER(I4B),PARAMETER,PRIVATE :: nbp0_ic = 2
      INTEGER(I4B),PARAMETER,Public  :: solv_size = 16
      CHARACTER(len=8) ,PUBLIC ::  bp0_icv(nbp0_ic),bp0_ic
      CHARACTER(len=solv_size) ,PUBLIC :: select_solver
      CHARACTER(len=std_char_size) ,PUBLIC, DIMENSION(:),ALLOCATABLE   :: nameu
      CHARACTER(len=std_char_size), PUBLIC, DIMENSION(:),POINTER       :: nameu_hist,nameg_hist
      CHARACTER(len=std_char_size), DIMENSION(:),ALLOCATABLE  :: var_type
      CHARACTER(len=std_char_size), DIMENSION(:),ALLOCATABLE  :: nameu_gp
      CHARACTER(len=1)dbgsvx_equed
      CHARACTER*10,DIMENSION(:),ALLOCATABLE  :: run_mode
      REAL(DP) ,   DIMENSION(:),ALLOCATABLE  :: eqtn_scale
      INTEGER(I4b),DIMENSION(:),ALLOCATABLE  :: grid_point
      INTEGER(I4b),DIMENSION(:),ALLOCATABLE  :: u_k_ind
      INTEGER(I4b),DIMENSION(:),ALLOCATABLE  :: map_k_2_pert
      INTEGER(I4b),DIMENSION(:),ALLOCATABLE  :: n_slave_grid_tot
      INTEGER(I4b)  freeze_xte,freeze_xti,freeze_xrbp,freeze_xwr
      INTEGER(I4b)  freeze_xnue,freeze_alpha_exp,freeze_xni
      INTEGER(I4b)  tot_iters,jac_skip,pick_eq_set
      INTEGER(I4b)   random_pert,conv_skip,dbgsvx_info
      INTEGER(I4b)  wrt_nwt,tipts,iters_freeze,iteration_method
      INTEGER(I4b)  tot_iters_max, non_lin_method,bandwidth
      INTEGER(I4b)  max_bandwidth,min_bandwidth,max_bwcycle
      INTEGER(I4b)  bwcycle,cycle_index,prev_cycle,lin_crit_grad
      INTEGER(I4b)  get_typf,jacobian_type,maxfev,outputvb,vloopvb,tglfvb,glfvb 
      INTEGER(I4b)  info_conv,fdigits,prt_turbm
      INTEGER(I4b)  create_plot_file,create_restart_file
      INTEGER(I4b)  save_incr_plot_file, nsteps_plot
      INTEGER(I4b)  save_incr_restart_file, nsteps_restart

      REAL(DP) plot_time_incr,restart_time_incr

      REAL(DP), DIMENSION(:),ALLOCATABLE ::typx,typf,sf,sx,             &
                            resid,xdep,fvp,dudtb,gamma_surfc,stable
      REAL(DP), DIMENSION(:),POINTER :: time_hist,ssqr_time_hist,       &
                                        te0_hist,ti0_hist,rbp0_hist,    &
                                        angrot0_hist,ni0_hist,dt_hist,  &
                                        gmax_hist,umax_hist,rg_hist,    &
                                        ru_hist
      REAL(DP)     ssqr,gradmax,ssqr_hist,gradtol,stepsize_factor,      &
                   dbar_stab,fvectol,steptol ,ssqrmin 

      LOGICAL xvset,jac_known,first_step,last_step,rescale,             &
              write_forcebal_summary




      DATA first_plot,first_restart_file /.TRUE.,.TRUE./
!      DATA select_solver /'full_newton'/
      DATA xvset,jac_known /.FALSE., .FALSE./ 
      DATA tipts /0/    !total iters per time step counter
      DATA ic_step,ip_step /0,0/    !counts number of predictor/corrector steps
      DATA time_tol,iteration_method /1.e-9,0/
      DATA bp0_icv /'eqdsk','analytic'/
      DATA icalln,cycle_index,prev_cycle,lin_crit_grad,icall_tglf_ss_fdjac   /0,1,1,0,0/
      DATA tglf_ss_jac /.FALSE./


  END MODULE SOLCON_GCNMP
