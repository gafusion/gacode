!-----------------------------------------------------
! gyro_make_poissonaperp_matrix.f90
!
! PURPOSE:
!  Define sparse form of Poisson-AmperePerp (PB) 
!  matrix and factorize using UMFPACK.
!-----------------------------------------------------

subroutine gyro_make_poissonaperp_matrix

  use gyro_globals
  use gyro_poissonaperp_private

  !------------------------
  implicit none
  !------------------------

  !----------------------------------------------------------------
  ! Begin by doing a hand-count of the nonzero elements in 
  ! the sparse PB matrix:
  !
  ! Row dimension of sparse PB submatrix
  n_poissonaperp_row = 2*n_x*n_blend

  ! radial bandwidth
  n_gyro = 2*m_gyro-i_gyro+1

  if (boundary_method == 1) then

     ! nonzero elements in n=0 PB matrix:
     n_zero = 4*(n_x-1)*n_gyro*n_blend**2+2*n_x*n_blend 

     ! nonzero elements in n>0 PB matrix:
     n_fini = 4*n_x*n_gyro*n_blend**2

  else

     ! nonzero elements in nonperiodic PB matrices:
     n_zero = 4*(n_x*n_gyro-m_gyro*(m_gyro+1))*n_blend**2 ! MPP,MPB,MBP,MBB

     n_fini = n_zero

  endif

  !----------------------------------------------------------------

  if (n_1(in_1) == 0) then
     n_poissonaperp = n_zero
  else
     n_poissonaperp = n_fini
  endif

  lindx(4)  = 2*n_poissonaperp
  lvalue(4) = n_poissonaperp
 
  allocate(m_poissonaperp(lvalue(4)))
  allocate(indx_poissonaperp(lindx(4)))

  if (n_1(in_1) == 0 .and. boundary_method == 1) then
     n_x_max  = n_x-1
  else
     n_x_max  = n_x
  endif

  ! Poisson MPP matrix
  allocate(ap_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
  if (electron_method /= 2) then

     ! Adiabatic electrons or ions: (2-R) phi

     call gyro_blend_poisson(1)

  else

     ! Kinetic electrons: (1-R) phi

     call gyro_blend_poisson(0)

  endif

  ! Ampere Perp MBB matrix, and the coupled matrices MBP and MPB
  allocate(ab_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
  allocate(abp_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
  call gyro_blend_ampereperp
  call gyro_make_jperp_elec

  call gyro_sparse_solve_umfpack(n_poissonaperp,n_poissonaperp_row,4,0)

  !---------------------------------------------
  ! These are large matrices and deallocation is
  ! important:
  !
  deallocate(ap_mm)
  deallocate(ab_mm)
  deallocate(abp_mm)
  !---------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_make_poissonaperp_matrix done]'
  endif

end subroutine gyro_make_poissonaperp_matrix
