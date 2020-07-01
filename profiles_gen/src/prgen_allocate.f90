subroutine prgen_allocate

  use prgen_globals

  implicit none

  allocate(rho(nx))
  allocate(dpsi(nx))
  allocate(rmin(nx))
  allocate(rmaj(nx))
  allocate(q(nx))
  allocate(p_tot(nx))
  allocate(omega0(nx))
  allocate(vpolc_exp(nx))
  allocate(vtorc_exp(nx))

  allocate(johm(nx))
  allocate(jbs(nx))
  allocate(jnb(nx))
  allocate(jrf(nx))

  allocate(kappa(nx)) ; kappa = 1.0
  allocate(delta(nx)) ; delta = 0.0
  allocate(zmag(nx))  ; zmag = 0.0
  allocate(zeta(nx))  ; zeta = 0.0
  allocate(shape_sin3(nx)) ; shape_sin3 = 0.0
  allocate(shape_cos0(nx)) ; shape_cos0 = 0.0
  allocate(shape_cos1(nx)) ; shape_cos1 = 0.0
  allocate(shape_cos2(nx)) ; shape_cos2 = 0.0
  allocate(shape_cos3(nx)) ; shape_cos3 = 0.0

  allocate(qspow_e(nx))
  allocate(qspow_i(nx))
  allocate(qpow_e(nx))
  allocate(qpow_i(nx))

end subroutine prgen_allocate
