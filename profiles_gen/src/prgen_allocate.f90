subroutine prgen_allocate(ftype)

  use prgen_globals

  implicit none

  character(len=*), intent(in) :: ftype

  allocate(rho(nx))
  allocate(dpsi(nx))
  allocate(rmin(nx))
  allocate(rmaj(nx))
  allocate(q(nx))
  allocate(p_tot(nx))
  allocate(fpol(nx))
  allocate(omega0(nx))
  allocate(vpolc_exp(nx))
  allocate(vtorc_exp(nx))
  allocate(te_kev(nx))
  allocate(ti_kev(nx))
  allocate(ne_e19m3(nx))
  allocate(ni_e19m3(nx))
  allocate(nz_e19m3(nx))
  allocate(zeff(nx))

  allocate(jtot(nx))
  allocate(johm(nx))
  allocate(jbs(nx))
  allocate(jnb(nx))
  allocate(jrf(nx))

  allocate(kappa(nx)) ; kappa = 1.0
  allocate(delta(nx)) ; delta = 0.0
  allocate(zmag(nx))  ; zmag = 0.0
  allocate(zeta(nx))  ; zeta = 0.0
  allocate(shape_sin(6,nx))   ; shape_sin = 0.0
  allocate(shape_cos(0:6,nx)) ; shape_cos = 0.0 

  allocate(qspow_e(nx))
  allocate(qspow_i(nx))
  allocate(qpow_e(nx))
  allocate(qpow_i(nx))
  allocate(qohm(nx))

  if (ftype == 'iterdb') then

     allocate(onetwo_ti(nx))
     allocate(onetwo_te(nx))
     allocate(onetwo_rho_grid(nx))
     allocate(onetwo_ene(nx))
     allocate(onetwo_enalp(nx))
     allocate(onetwo_angrot(nx))
     allocate(onetwo_psi(nx))

     allocate(onetwo_enion(nx,onetwo_nion))
     allocate(onetwo_enbeam(nx,onetwo_nbion))
     allocate(onetwo_enn(nx,onetwo_nbion))
     allocate(onetwo_ennw(nx,onetwo_nbion))
     allocate(onetwo_ennv(nx,onetwo_nbion))
     allocate(onetwo_pressb(nx,onetwo_nbion))

     allocate(onetwo_sion(nx,onetwo_nion))
     allocate(onetwo_srecom(nx,onetwo_nion))
     allocate(onetwo_scx(nx,onetwo_nion))
     allocate(onetwo_sbcx(nx,onetwo_nion))
     allocate(onetwo_s(nx,onetwo_nion))
     allocate(onetwo_dudt(nx,onetwo_nion))

     allocate(onetwo_rho_mhd_gridnpsi(onetwo_npsi))
     allocate(onetwo_rmajavnpsi(onetwo_npsi))
     allocate(onetwo_rminavnpsi(onetwo_npsi))
     allocate(onetwo_elongxnpsi(onetwo_npsi))
     allocate(onetwo_psivolpnpsi(onetwo_npsi))
     allocate(onetwo_triangnpsi_u(onetwo_npsi))
     allocate(onetwo_triangnpsi_l(onetwo_npsi))

     allocate(onetwo_volume(nx))
     allocate(onetwo_hcap(nx)) 
     allocate(onetwo_qbeame(nx)) 
     allocate(onetwo_qrfe(nx)) 
     allocate(onetwo_qrad(nx)) 
     allocate(onetwo_qione(nx)) 
     allocate(onetwo_qioni(nx)) 
     allocate(onetwo_dpedt(nx)) 
     allocate(onetwo_qfuse(nx)) 
     allocate(onetwo_qbeami(nx)) 
     allocate(onetwo_qrfi(nx)) 
     allocate(onetwo_qcx(nx))
     allocate(onetwo_qsync(nx)) 
     allocate(onetwo_dpidt(nx)) 
     allocate(onetwo_qfusi(nx)) 
     allocate(onetwo_qdelt(nx)) 
     allocate(sion_d(nx)) 
     allocate(sbcx_d(nx)) 
     allocate(onetwo_sbeam(nx)) 
     allocate(onetwo_sbeame(nx)) 
     allocate(onetwo_sscxl(nx)) 
     allocate(onetwo_storqueb(nx))

     allocate(onetwo_enion_vec(n_ion_max,nx))
     allocate(onetwo_Tion_vec(n_ion_max,nx))

  else if (ftype == 'pfile') then

     allocate(peqdsk_psi(nx))
     allocate(peqdsk_ne(nx))
     allocate(peqdsk_te(nx))
     allocate(peqdsk_ni(nx))
     allocate(peqdsk_ti(nx))
     allocate(peqdsk_nz(3,nx))
     allocate(peqdsk_nb(nx))
     allocate(peqdsk_pb(nx))
     allocate(peqdsk_omegat(nx))
     allocate(peqdsk_omgeb(nx))
     allocate(peqdsk_pow_i(nx))
     allocate(peqdsk_pow_e(nx))

  else if (ftype == 'swim') then

     allocate(plst_alla_name(plst_dp1_nspec_alla))
     allocate(plst_all_name(plst_dp1_nspec_all))
     allocate(plst_q_all(plst_dp1_nspec_all))
     allocate(plst_m_all(plst_dp1_nspec_all))
     allocate(plst_ns(nx,plst_dp1_nspec_all))
     allocate(plst_ts(nx,plst_dp1_nspec_all))
     allocate(plst_nb(nx,plst_dim_nspec_beam))
     allocate(plst_nmini(nx))
     allocate(plst_nfusi(nx))
     allocate(plst_tb(nx,plst_dim_nspec_beam))
     allocate(plst_tmini(nx))
     allocate(plst_tfusi(nx))
     allocate(plst_epar(nx))
     allocate(plst_eperp(nx))
     allocate(plst_bepar(nx,plst_dim_nspec_beam))
     allocate(plst_beperp(nx,plst_dim_nspec_beam))
     allocate(plst_vol(nx))
     allocate(plst_surf(nx))
     allocate(plst_phit(nx))
     allocate(plst_psipol(nx))
     allocate(plst_elong(nx))
     allocate(plst_triang(nx))
     allocate(plst_triangu(nx))
     allocate(plst_iota(nx))
     allocate(plst_r_midp_in(nx))
     allocate(plst_r_midp_out(nx))
     allocate(plst_z_midp(nx))
     allocate(plst_epot(nx))
     allocate(plst_omegat(nx))

     allocate(plst_pe_trans(nx))
     allocate(plst_pi_trans(nx))
     allocate(plst_peech(nx))
     allocate(plst_picrf_totals(nx,plst_dp1_nspec_all))
     allocate(plst_pmini(nx))
     allocate(plst_pminth(nx))
     allocate(plst_picth(nx))
     allocate(plst_pmine(nx))
     allocate(plst_pohme(nx))
     allocate(plst_qie(nx))
     allocate(plst_pbe(nx))
     allocate(plst_pbi(nx))
     allocate(plst_pbth(nx))
     allocate(plst_pfusi(nx))
     allocate(plst_pfusth(nx))
     allocate(plst_pfuse(nx))
     allocate(plst_prad(nx))
     allocate(plst_prad_br(nx))
     allocate(plst_prad_cy(nx))
     allocate(plst_prad_li(nx))
     allocate(plst_tq_trans(nx))
     allocate(plst_sn_trans(nx))

     allocate(plst_curr_ohmic(nx))
     allocate(plst_curr_bootstrap(nx))
     allocate(plst_curt(nx))

  endif

end subroutine prgen_allocate
