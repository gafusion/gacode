subroutine allocate_peqdsk_vars

  use prgen_globals

  implicit none

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

end subroutine allocate_peqdsk_vars
