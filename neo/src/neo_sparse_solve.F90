module neo_sparse_solve

  implicit none

  public :: SOLVE_sparse

contains

  subroutine SOLVE_sparse(m_size, m_iindx, m_jindx, m, m_asize, m_rows, rhs, sol)
    use neo_globals
    use neo_umfpack
    implicit none
    integer, intent (in) :: m_size                !number of matrix values in COO format
    integer, dimension(:), intent (in) :: m_iindx !matrix row-indices of COO format
    integer, dimension(:), intent (in) :: m_jindx !matrix col-indices of COO format
    real, dimension(:), intent (in) :: m          !matrix values of COO format
    integer, intent (in) :: m_asize               !total size of array m
    integer, intent (in) :: m_rows                !number of rows in matrix
    real, dimension(:), intent (in) :: rhs        !right-hand-side vector
    real, dimension(:), intent (inout) :: sol     !solution vector

    if (use_petsc == 1) then
#ifdef NEO_HAVE_PETSC
      call SOLVE_petsc(m_size, m_iindx, m_jindx, m, m_asize, m_rows, rhs, sol)
      if(error_status == 0) return
#endif
      call neo_error('ERROR: (NEO) PetSc package not used')
      error_status = 0 ! try with other solvers
    endif

    if (use_slu == 1) then
#ifdef NEO_HAVE_SUPERLU
      call SOLVE_superlu(m_size, m_iindx, m_jindx, m, m_asize, m_rows, rhs, sol)
      if(error_status == 0) return
#endif
      call neo_error('ERROR: (NEO) SuperLU package not used')
      error_status = 0 ! try with other solvers
    endif

    !UMFPACK is the default option
    call SOLVE_umfpack(m_size, m_iindx, m_jindx, m, m_asize, m_rows, rhs, sol)

  end subroutine SOLVE_sparse

end module neo_sparse_solve
