module tgyro_iteration_variables

  !------------------------------------------------------
  ! Shared/standard variables
  !
  integer :: correct_flag
  !
  integer :: i
  integer :: p
  integer :: pp
  integer :: ip
  integer :: i_tran_loop
  real :: dx
  !
  ! To be allocated in tgyro_iteration_driver: 
  !
  integer, dimension(:), allocatable :: ipiv
  real, dimension(:,:), allocatable :: jf
  real, dimension(:,:), allocatable :: jg
  real, dimension(:,:), allocatable :: jfg
  real, dimension(:), allocatable :: x_vec
  real, dimension(:), allocatable :: f_vec
  real, dimension(:), allocatable :: g_vec
  real, dimension(:), allocatable :: x_vec0
  real, dimension(:), allocatable :: f_vec0
  real, dimension(:), allocatable :: g_vec0
  real, dimension(:), allocatable :: b
  !------------------------------------------------------

  !------------------------------------------------------
  ! Blocked serial/parallel variables:
  !
  integer :: i_worker
  integer :: search_index=1 
  integer :: search_max=4
  real, dimension(6,4) :: search
  ! Idea for searching is to try
  ! full Newton -> 1.0
  ! tiny step   -> 0.1
  ! larger      -> 0.75
  ! smaller     -> 0.5, etc.
  data search/&
        1.0, 0.1, 0.75, 0.5, 0.25, 0.5,&
       -1.0,-0.1,-0.75,-0.5,-0.25,-0.5,&
        2.0, 0.2, 1.50, 1.1, 0.60, 1.1,&
       -2.0,-0.2,-1.50,-1.1,-0.60,-1.1/
  !
  ! To be allocated in tgyro_iteration_serial:
  !
  real, dimension(:), allocatable :: res1
  real, dimension(:), allocatable :: x_vec1
  real, dimension(:), allocatable :: f_vec1
  real, dimension(:), allocatable :: g_vec1
  !
  ! To be further allocated in tgyro_iteration_parallel:
  !
  real, dimension(:,:), allocatable :: res_vec
  real, dimension(:,:), allocatable :: x_vec_vec
  real, dimension(:,:), allocatable :: f_vec_vec
  real, dimension(:,:), allocatable :: g_vec_vec
  !------------------------------------------------------

end module tgyro_iteration_variables
