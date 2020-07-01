subroutine prgen_swap

  use expro
  
  implicit none

  ! Perform reordering (may just be identity)

  call swapc(expro_name)
  call swapc(expro_type)
  call swap1(expro_z,expro_n_ion)
  call swap1(expro_mass,expro_n_ion)
  call swap2(expro_ni,expro_n_ion,expro_n_exp)
  call swap2(expro_ti,expro_n_ion,expro_n_exp)
  call swap2(expro_vpol,expro_n_ion,expro_n_exp)
  call swap2(expro_vtor,expro_n_ion,expro_n_exp)

end subroutine prgen_swap

subroutine swap1(f,nion)

  use prgen_globals

  implicit none

  integer, intent(in) :: nion
  real, intent(inout), dimension(nion) :: f

  integer :: i,ip
  real, dimension(nion) :: ftmp

  ftmp = f
  do i=1,nion
     ip   = reorder_vec(i)
     f(i) = ftmp(ip)
  enddo

end subroutine swap1

subroutine swap2(f,nion,nexp)

  use prgen_globals

  implicit none

  integer, intent(in) :: nion,nexp
  real, intent(inout), dimension(nion,nexp) :: f

  integer :: i,ip
  real, dimension(nion,nexp) :: ftmp

  ftmp = f
  do i=1,nion
     ip     = reorder_vec(i)
     f(i,:) = ftmp(ip,:)
  enddo

end subroutine swap2

subroutine swapc(f)

  use prgen_globals

  implicit none

  character(len=10), intent(inout), dimension(20) :: f

  integer :: i,ip
  character(len=10), dimension(20) :: ftmp
  
  ftmp = f
  do i=1,10
     ip   = reorder_vec(i)
     f(i) = ftmp(ip)
  enddo

end subroutine swapc
