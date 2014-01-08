      MODULE glf23_gcnmp
      USE nrtype,                      ONLY : DP,I4B
      IMPLICIT NONE
!
!
!    jmaxm is used both as a dimension in glf23 and as the size
!    of the working arrays. For 1d arrays this is not a
!    problem but for the single 2d array, egamma_d, that is passed to
!    glf this requires special treatment. Hence this include file was turned
!    into a module instead. HSJ 03/25/02
       INTEGER(I4B),PARAMETER :: itport_glf    = 5
       INTEGER(I4B),PARAMETER :: nspecies_glf   = 3             !  electrons,prim ion,fast ions
       INTEGER(I4B)  navg_chi, glf23_iglf, reset_mask,         &
              init_egamma_exp,freeze_alpha_exp
       INTEGER,SAVE :: first_step
       INTEGER(I4B), SAVE, DIMENSION(:), ALLOCATABLE:: itport_pt
       REAL(DP) iglf_chie,iglf_chii,iglf_chiv,iglf_d
       REAL(DP), SAVE,DIMENSION(:,:), ALLOCATABLE :: egamma_d ,dcoefp_glf23
       REAL(DP), SAVE,DIMENSION(:), ALLOCATABLE ::             &
              te_glf,ti_glf,ne_glf, ni_glf,                    &
              q_glf,  shat_exp,  elong_glf,  zeff_exp,         &
              rho_glf,  grho1_glf,  grho2_glf,  gradte_glf,    &
              zpte_glf,  gradti_glf,  zpti_glf,                &
              gradne_glf,  zpne_glf,  gradni_glf,  tau_glf,    &
              zpni_glf,  betae_glf,  betai_glf,  beta_glf,     &
              alpha_exp,  rmin_glf,  rmaj_glf,  drhogf,        &
              ns_glf,  bteff_glf,alpha_exp_save,               &
              csda_glf,  rhosda_glf,  drhodr,   drhodrrrho,    &
              angrotp_exp,  egamma_exp,  gamma_p_exp,          &
              vphim_glf,  vparm_glf,  vperm_glf,  ve_glf,      &
              fcglf,  ak1glf,  ak2glf,  aneoglf,vpolm_glf,     &
              diff_m,  chie_m,  chii_m,                        &
              etaphi_m,  etapar_m,  etaper_m,                  &
              exch_m , egamma_m ,                              &
              gamma_p_m ,anrate_m , anrate2_m ,                &
              anfreq_m , anfreq2_m ,vp_eff_m,                  &
              ve_eff_m,vi_eff_m,vrot_eff_m,                    &
              diff_mdv,  chie_mdv,  chii_mdv,                  &
              etaphi_mdv,  etapar_mdv,  etaper_mdv,            &
              exch_mdv , egamma_mdv,egamma_exp_initial,        &
               gamma_p_mdv ,anrate_mdv , anrate2_mdv ,         &
              anfreq_mdv , anfreq2_mdv,zpti_save,              &
              chii_save,chi_e_avg,chi_i_avg,exchmdv,chii_ms,   &
              chie_ms,chie_avg,chii_avg,etaphi_avg,            &
              etaphi_ms,d_avg,xkang_avg,xchie_sv,xchii_sv,     &
              xkangwlsv,xke_glf23,xki_glf23,xkw_glf23,diff_ms, &
              diff_avg,niflux_gf_m, nimpflux_gf_m,tiflux_gf_m, &
              teflux_gf_m,vparflux_gf_m,vperflux_gf_m,         &
              vphiflux_gf_m,etgflux_gf_m,d_sv,                 &
              niflux_gf_m_loc,nimpflux_gf_m_loc,               &
              tiflux_gf_m_loc,teflux_gf_m_loc,                 &
              vparflux_gf_m_loc,vperflux_gf_m_loc,             &
              vphiflux_gf_m_loc,etgflux_gf_m_loc,gamma_net_i,  &
              gamma_net_e,gamma_net_i_loc,gamma_net_e_loc,     &

              glf_te_fg,glf_ti_fg,glf_ne_fg,glf_ni_fg,         &
              glf_ns_fg,glf_angrot_fg,glf_exch_m,              &
              glf_egamma_exp,glf_gamma_p_exp,glf_vphim,        &
              glf_vparm,glf_vperm,glf_zeff,glf_rho,glf_grho1,  &
              glf_grho2,glf_rmin,glf_rmaj,glf_q,glf_shat,      &
              glf_alpha,glf_elong,glf_diff_m,glf_chie_m,       &
              glf_chii_m,glf_etaphi_m,glf_etapar_m,            &
              glf_etaper_m,glf_egamma_m,                       &
              glf_gamma_p_m,glf_anrate,glf_anrate2,glf_anfreq, &
              glf_anfreq2,glf_ni_flux,glf_nimp_flux,           &
              glf_ti_flux,glf_te_flux,glf_vparflux_gf_m,       &
              glf_vperflux_gf_m,glf_vphiflux_gf_m,             &
              glf_etgflux_gf_m,glf_gamma_net_i,glf_gamma_net_e,&
              glf_egamma_zc,glf_gamma_p_zc,dqrm_drho,          &
              dvexb_base_drho,vexb_base,vexb_term1,diamag_term,&
              vperm_term1,vparm_term1,ddiamag_drho,            &
              dvexb_term1_drho,dvpar_term1_drho,dvphim_drho,   &
              egamma_exp_zc,gamma_p_exp_zc,dangrot_drho,rdrho
              





     REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:)              ::  &
          glf_send_packed, glf_receive_packed,glf_typ_val,     &
          glf_etg_output,glf_gamma_net_i_output,               &
          glf_gamma_net_e_output,glf_anfreq_output,            &
          glf_anfreq2_output,glf_anrate_output,                &
          glf_anrate2_output
     REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:)            ::  &
          glf_p_output,glf_e_output,glf_val_pert,              &
          glf_val_base,glf_egamma_d,glf_m_output,              &
          glf_eff_chi_output
     REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:,:)          ::  &
             glf_p_flux,glf_e_flux,glf_m_flux,glf_vp,          &
             glf_vm,glf_d

     REAL(DP), SAVE     ::                                     &
                          nimpflux_gf,tiflux_gf,teflux_gf,     &     
                          vparflux_gf,vperflux_gf,vphiflux_gf, &
                          etgflux_gf,glf_v_mult,glf_d_mult,    &
                          niflux_gf,glf_corot

!
      REAL(DP)  chie_dv(-1:1),chii_dv(-1:1),chitorot_dv(-1:1), &
              zptij_p(-1:1),zptej_p(-1:1)
      INTEGER(I4B),PARAMETER       ::               pdelay =10 
      INTEGER(I4B) glf_max_pack
      LOGICAL recall_glf,use_glf_pert,glf_unstable

!     INPUT switches:
      REAL(DP) exbmult_glf,x_alpha_glf,dv_delt,t_delay_glf23
      INTEGER(I4B) jroot_glf,limp_glf,ibtflag_glf,irotstab,    &
                   dv_method,spline_gsc,                       &
                   itte_dv,itti_dv,itangrot_dv,itenp_dv,       &
                   itene_dv,jeigen_glf,i_delay,use_mask
!     pdelay is sized according to i_delay
!     the value 10 is hardwired in glf23 for dimension of egamma_d
 
        END MODULE  glf23_gcnmp
