subroutine tgyro_residual(f,g,res,n,method)

  use tgyro_iteration_variables

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: method
  real, intent(in), dimension(n) :: f
  real, intent(in), dimension(n) :: g
  real, intent(inout), dimension(n) :: res
  
  select case (method)

  case (2)

     ! ABSOLUTE VALUE NORM
     res = abs(f-g)

  case (3)

     ! SQUARE RESIDUAL
     res = (f-g)**2

  case default

     print '(a)','ERROR: (tgyro_residual) Invalid choice of residual function.'
     stop

  end select
  
end subroutine tgyro_residual
