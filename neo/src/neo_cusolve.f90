subroutine SOLVE_cuda(nnz_coo, icoo, jcoo, vcoo, m_asize, n, rhs, sol)
  use matconv
  implicit none

  integer, intent (in) :: nnz_coo
  integer, dimension(*), intent (in) :: icoo
  integer, dimension(*), intent (in) :: jcoo
  real, dimension(*), intent (in) :: vcoo
  integer, intent (in) :: m_asize
  integer, intent (in) :: n
  real, dimension(*), intent (in) :: rhs
  real, dimension(*), intent (inout) :: sol

  integer :: nnz_csr
  integer, allocatable, target :: icsr(:), jcsr(:)
  real, allocatable, target :: vcsr(:)

  call coo_to_csr( n, nnz_coo, icoo, jcoo, vcoo, nnz_csr, icsr, jcsr, vcsr )
  call cusolve_sparse(n, nnz_csr, vcsr, icsr, jcsr, rhs, sol)

  if (allocated(icsr)) deallocate(icsr)
  if (allocated(jcsr)) deallocate(jcsr)
  if (allocated(vcsr)) deallocate(vcsr)

end subroutine SOLVE_cuda
