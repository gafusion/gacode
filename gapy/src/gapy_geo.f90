subroutine geov(t,tc,n)

  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: t
  double precision, intent(inout), dimension(n) :: tc

!f2py intent(in) n
!f2py intent(in) t
!f2py intent(out) tc

  tc = sin(t)
  
end subroutine geov
