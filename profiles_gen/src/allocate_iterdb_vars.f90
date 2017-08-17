subroutine allocate_iterdb_vars

  use prgen_globals

  implicit none

  allocate(onetwo_ti(nx))
  allocate(onetwo_te(nx))
  allocate(onetwo_rho_grid(nx))
  allocate(onetwo_ene(nx))
  allocate(onetwo_enalp(nx))
  allocate(onetwo_zeff(nx))
  allocate(onetwo_angrot(nx))
  allocate(onetwo_psi(nx))

  allocate(onetwo_enion(nx,onetwo_nion))
  allocate(onetwo_enbeam(nx,onetwo_nbion))
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
  allocate(onetwo_qohm(nx)) 
  allocate(onetwo_qrad(nx)) 
  allocate(onetwo_qione(nx)) 
  allocate(onetwo_qioni(nx)) 
  allocate(onetwo_dpedt(nx)) 
  allocate(onetwo_qfuse(nx)) 
  allocate(onetwo_qbeami(nx)) 
  allocate(onetwo_qrfi(nx)) 
  allocate(onetwo_qcx(nx)) 
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

end subroutine allocate_iterdb_vars
