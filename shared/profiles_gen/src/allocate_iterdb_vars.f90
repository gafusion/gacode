subroutine allocate_iterdb_vars

  use prgen_globals

  implicit none

  allocate(onetwo_ti(onetwo_nj))
  allocate(onetwo_te(onetwo_nj))
  allocate(onetwo_rho_grid(onetwo_nj))
  allocate(onetwo_ene(onetwo_nj))
  allocate(onetwo_enalp(onetwo_nj))
  allocate(onetwo_zeff(onetwo_nj))
  allocate(onetwo_angrot(onetwo_nj))
  allocate(onetwo_psi(onetwo_nj))

  allocate(onetwo_enion(onetwo_nj,onetwo_nion))
  allocate(onetwo_enbeam(onetwo_nj,onetwo_nbion))
  allocate(onetwo_pressb(onetwo_nj,onetwo_nbion))
  allocate(onetwo_press(onetwo_nj))

  allocate(onetwo_sion(onetwo_nj,onetwo_nion))
  allocate(onetwo_srecom(onetwo_nj,onetwo_nion))
  allocate(onetwo_scx(onetwo_nj,onetwo_nion))
  allocate(onetwo_sbcx(onetwo_nj,onetwo_nion))
  allocate(onetwo_s(onetwo_nj,onetwo_nion))
  allocate(onetwo_dudt(onetwo_nj,onetwo_nion))

  allocate(onetwo_rho_mhd_gridnpsi(onetwo_npsi))
  allocate(onetwo_rmajavnpsi(onetwo_npsi))
  allocate(onetwo_rminavnpsi(onetwo_npsi))
  allocate(onetwo_elongxnpsi(onetwo_npsi))
  allocate(onetwo_psivolpnpsi(onetwo_npsi))
  allocate(onetwo_triangnpsi_u(onetwo_npsi))
  allocate(onetwo_triangnpsi_l(onetwo_npsi))

  allocate(onetwo_volume(onetwo_nj))
  allocate(onetwo_hcap(onetwo_nj)) 
  allocate(onetwo_qbeame(onetwo_nj)) 
  allocate(onetwo_qrfe(onetwo_nj)) 
  allocate(onetwo_qohm(onetwo_nj)) 
  allocate(onetwo_qrad(onetwo_nj)) 
  allocate(onetwo_qione(onetwo_nj)) 
  allocate(onetwo_qioni(onetwo_nj)) 
  allocate(onetwo_dpedt(onetwo_nj)) 
  allocate(onetwo_qfuse(onetwo_nj)) 
  allocate(onetwo_qbeami(onetwo_nj)) 
  allocate(onetwo_qrfi(onetwo_nj)) 
  allocate(onetwo_qcx(onetwo_nj)) 
  allocate(onetwo_dpidt(onetwo_nj)) 
  allocate(onetwo_qfusi(onetwo_nj)) 
  allocate(onetwo_qdelt(onetwo_nj)) 
  allocate(sion_d(onetwo_nj)) 
  allocate(sbcx_d(onetwo_nj)) 
  allocate(onetwo_sbeam(onetwo_nj)) 
  allocate(onetwo_storqueb(onetwo_nj))
  
  
  allocate(onetwo_enion_vec(onetwo_nion+onetwo_nbion+1,onetwo_nj))
  allocate(onetwo_Tion_vec(onetwo_nion+onetwo_nbion+1,onetwo_nj))

end subroutine allocate_iterdb_vars
