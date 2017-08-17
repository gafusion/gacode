subroutine allocate_ufile_vars

  use prgen_globals

  implicit none

  allocate(ufile_ne(nx))
  allocate(ufile_ni(nx,n_ion_max))
  ufile_ni = 0.0
  allocate(ufile_te(nx))
  allocate(ufile_ti(nx,n_ion_max))
  ufile_ti = 0.0
  allocate(ufile_zeff(nx))
  allocate(ufile_pres(nx))
  allocate(ufile_vrot(nx))
  allocate(ufile_vrotm(nx))
  allocate(ufile_volume(nx))
  allocate(ufile_qnbii(nx))
  allocate(ufile_qnbie(nx))
  allocate(ufile_qicrhi(nx))
  allocate(ufile_qicrhe(nx))
  allocate(ufile_qei(nx))
  allocate(ufile_qrad(nx))
  allocate(ufile_qeche(nx))
  allocate(ufile_qechi(nx))
  allocate(ufile_qohm(nx))
  allocate(ufile_qwalli(nx))
  allocate(ufile_qwalle(nx))

  allocate(ufile_qfuse(nx))
  allocate(ufile_qfusi(nx))
  allocate(ufile_qlhe(nx))
  allocate(ufile_qlhi(nx))
  allocate(ufile_snbie(nx))
  allocate(ufile_snbii(nx))
  allocate(ufile_swall(nx))
  allocate(ufile_torq(nx))

end subroutine allocate_ufile_vars
