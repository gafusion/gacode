subroutine allocate_corsica_vars

  use prgen_globals

  implicit none

  allocate(corsica_rho(corsica_nvals))
  allocate(corsica_r_a(corsica_nvals))
  allocate(corsica_psin(corsica_nvals))
  allocate(corsica_vl(corsica_nvals))
  allocate(corsica_te(corsica_nvals))
  allocate(corsica_ti(corsica_nvals))
  allocate(corsica_ne(corsica_nvals))
  allocate(corsica_ndt(corsica_nvals))
  allocate(corsica_nz(corsica_nvals))
  allocate(corsica_nalpha(corsica_nvals))
  allocate(corsica_zeff(corsica_nvals))
  allocate(corsica_q(corsica_nvals))
  allocate(corsica_j(corsica_nvals))
  allocate(corsica_jbs(corsica_nvals))

end subroutine allocate_corsica_vars
