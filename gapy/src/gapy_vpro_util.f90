subroutine vpro_compute_derived

  use vpro
  use util

  implicit none

  double precision :: rhod(nexp)

  rhod(:) = arho_exp*rho(:)

  ! b_unit
  call util_bound_deriv(bunit,rhod**2,rmin**2,nexp)
  bunit(:) = bt_exp*bunit

  ! s
  call util_bound_deriv(s,q,rmin,nexp)
  s(:) = (rmin(:)/q(:))*s(:)

end subroutine vpro_compute_derived

