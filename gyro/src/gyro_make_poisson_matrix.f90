!-----------------------------------------------------
! gyro_make_poisson_matrix.f90
!
! PURPOSE:
!  Define sparse form of Poisson matrix used for 
!  EXPLICIT field solve.
!-----------------------------------------------------

subroutine gyro_make_poisson_matrix

  use gyro_globals
  use math_constants
  use gyro_poisson_private

  !------------------------
  implicit none
  !------------------------

  !----------------------------------------------------------------
  ! Begin by doing a hand-count of the nonzero elements in 
  ! the sparse Maxwell matrix:
  !
  ! ELECTROSTATIC

  ! Row dimension of sparse Poisson matrix
  n_poisson_row = n_x*n_blend

  ! radial band width
  n_gyro = 2*m_gyro-i_gyro+1

  if (boundary_method == 1) then

     ! nonzero elements in n=0 PP matrix:
     n_zero = (n_x-1)*n_gyro*n_blend**2+n_x*n_blend

     ! nonzero elements in n>0 PP matrix:
     n_fini = n_x*n_gyro*n_blend**2

  else

     ! nonzero elements in n=0 nonperiodic PP matrix:
     n_zero = (n_x*n_gyro-m_gyro*(m_gyro+1))*n_blend**2 ! MPP

     ! same number of nonzeros in n>0 nonperiodic matrices:
     n_fini = n_zero

  endif
  !----------------------------------------------------------------


  if (n_1(in_1) == 0) then
     n_poisson = n_zero
  else
     n_poisson = n_fini
  endif

  lindx(1)  = 2*n_poisson
  lvalue(1) = n_poisson

  allocate(m_poisson(lvalue(1)))
  allocate(indx_poisson(lindx(1)))

  if (n_1(in_1) == 0 .and.  boundary_method == 1) then
     n_x_max  = n_x-1
  else
     n_x_max  = n_x
  endif

  allocate(ap_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))

  if (electron_method /= 2) then

     ! Adiabatic electrons or ions: (2-R) phi

     call gyro_blend_poisson(1)

  else

     ! Kinetic electrons: (1-R) phi

     call gyro_blend_poisson(0)

  endif

  call gyro_sparse_solve_umfpack(n_poisson,n_poisson_row,1,0)

  !--------------------------------------------------------
  ! These are large matrices and deallocation is important:
  !
  deallocate(ap_mm)
  !--------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_make_poisson_matrix done]'
  endif

end subroutine gyro_make_poisson_matrix
