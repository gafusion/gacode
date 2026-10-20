!---------------------------------------------------------------
! neo_compute_fcoll(m0,lambda,f,fb)
!
! INPUT
!  m0     (integer)
!  lambda (real)
!
! OUTPUT
!  f (-m0:m0,-m0:m0)  (real)
!  fb(-m0:m0,-m0:m0)  (real)
!
!  Where f is the normalized field-particle collision integral: 
!
!                        F(m,n,lambda)
!
!  and fb is the complement:  
!
!                        FB(m,n,lambda)
!
! NOTES:
!  We use the internal variables g,gb in terms of which
!
!               G(m,n,lambda) = GB(n,m,1/lambda)
!
!      G(m,n,lambda) + GB(m,n,lambda) = beta[(m+1)/2,(n+1)/2]
!
!  Results should be good to machine precision in all cases
!----------------------------------------------------------------

subroutine neo_compute_fcoll(m0,lambda_in,f,fb)

  use iso_fortran_env, only: qp => real128
  implicit none

  ! arguments
  integer, intent(in) :: m0
  real(qp), intent(in) :: lambda_in
  real(qp), dimension(-m0:m0,-m0:m0), intent(inout) :: f,fb

  ! internal variables
  real(qp), dimension(:,:), allocatable :: g,gb
  real(qp), dimension(:), allocatable :: gamma
  real(qp), dimension(:,:), allocatable :: beta
  real(qp) :: lambda
  real(qp) :: pi
  real(qp) :: x
  real(qp) :: r
  real(qp) :: c
  integer :: i
  integer :: m
  integer :: n

  real(qp), parameter :: lambda_large = 10.0_qp

  pi = 4.0_qp*atan(1.0_qp)

  lambda = lambda_in

  allocate(g(-m0:m0,-m0:m0))
  allocate(gb(-m0:m0,-m0:m0))
  allocate(gamma(m0+m0+2))
  allocate(beta(0:m0,0:m0))

  !-------------------------------------------------
  ! Let GAMMA be the true gamma function; then
  !
  ! gamma(i) = GAMMA(i/2), or
  !
  ! gamma(1) = GAMMA(1/2)
  ! gamma(2) = GAMMA(1)
  ! gamma(3) = GAMMA(3/2)
  ! gamma(4) = GAMMA(2)
  ! gamma(5) = GAMMA(5/2)
  !
  gamma(2) = 1.0
  do i=1,m0
     gamma(2*(i+1)) = i*gamma(2*i)
  enddo
  !
  gamma(1) = sqrt(pi)
  do i=1,m0
     gamma(2*i+1) = (i-0.5_qp)*gamma(2*i-1)
  enddo
  !-------------------------------------------------

  !-------------------------------------------------
  ! Also need the beta function
  !
  do m=0,m0
     do n=0,m0
        beta(m,n) = gamma(m+1)*gamma(n+1)/gamma(m+n+2)
     enddo
  enddo
  !-------------------------------------------------

  g(:,:)  = 0.0
  gb(:,:) = 0.0

  !-------------------------------------------------------
  ! STAGE A: Compute GB(m,n,lambda)
  !-------------------------------------------------------

  r = 1.0_qp/(1.0_qp+lambda)

  do m=0,m0,2
     do n=1,m0,2

        !------------------------
        ! CASE 1: m even, n odd
        !------------------------
       
        x = 0.0
        do i=0,(n-1)/2
           x = x+(lambda/(1.0_qp+lambda))**i/(1.0_qp+lambda)**((m+1)/2.0_qp)*&
                gamma(2*i+m+1)/gamma(m+n+2)*&
                gamma(n+1)/gamma(2*i+2)
        enddo
        gb(m,n) = x

     enddo
  enddo

  do m=1,m0,2
     do n=0,m0,2

        !------------------------
        ! CASE 2: m odd, n even
        !------------------------

        if (lambda < lambda_large) then

           ! lambda is small enough to use original sum

           x = 0.0
           do i=0,(m-1)/2
              x = x+r**i*gamma(2*i+1)/(gamma(2*i+2)*gamma(1))
           enddo

           gb(m,n) = beta(m,n)*(1.0_qp-x*sqrt(lambda/(1.0_qp+lambda)))

        else

           ! Use complementary large-lambda sum (16 digits of precision)

           ! c(mbar)
           c = gamma(m)/(gamma(m+1)*gamma(1))*r**((m-1)/2)

           x = 0.0
           do i=(m+1)/2,(m+1)/2+int(16*log(10.0_qp)/log(lambda))
              c = c*(i-0.5_qp)/(i)*r
              x = x+c
           enddo

           gb(m,n) = beta(m,n)*x*sqrt(lambda/(1.0_qp+lambda))

        endif

        x = 0.0
        do i=0,n/2-1
           x = x+(lambda/(1+lambda))**(i+0.5_qp)/(1+lambda)**(0.5_qp+m/2.0_qp)*&
                gamma(2*i+2+m)/gamma(m+n+2)*&
                gamma(n+1)/gamma(2*i+3)
        enddo
        gb(m,n) = gb(m,n)+x

     enddo
  enddo

  !-------------------------------------------------------
  ! STAGE B: Compute G(n,m,lambda) = GB(m,n,1/lambda)
  !-------------------------------------------------------

  lambda = 1/lambda_in
  r = 1.0_qp/(1.0_qp+lambda)

  do m=0,m0,2
     do n=1,m0,2

        !------------------------
        ! CASE 1: m even, n odd
        !------------------------

        x = 0.0
        do i=0,(n-1)/2
           x = x+(lambda/(1.0_qp+lambda))**i/(1.0_qp+lambda)**((m+1)/2.0_qp)*&
                gamma(2*i+m+1)/gamma(m+n+2)*&
                gamma(n+1)/gamma(2*i+2)
        enddo

        g(n,m) = x

     enddo
  enddo

  do m=1,m0,2
     do n=0,m0,2

        !------------------------
        ! CASE 2: m odd, n even
        !------------------------

        if (lambda < lambda_large) then

           ! lambda is small enough to use original sum

           x = 0.0
           do i=0,(m-1)/2
              x = x+r**i*gamma(2*i+1)/(gamma(2*i+2)*gamma(1))
           enddo

           g(n,m) = beta(m,n)*(1.0_qp-x*sqrt(lambda/(1.0_qp+lambda)))

        else

           ! Use complementary large-lambda sum (16 digits of precision)

           ! c(mbar)
           c = gamma(m)/(gamma(m+1)*gamma(1))*r**((m-1)/2)

           x = 0.0
           do i=(m+1)/2,(m+1)/2+int(16*log(10.0_qp)/log(lambda))
              c = c*(i-0.5_qp)/i*r
              x = x+c
           enddo

           g(n,m) = beta(m,n)*x*sqrt(lambda/(1.0_qp+lambda))

        endif

        x = 0.0
        do i=0,n/2-1
           x = x+(lambda/(1+lambda))**(i+0.5_qp)/(1+lambda)**(0.5_qp+m/2.0_qp)*&
                gamma(2*i+2+m)/gamma(m+n+2)*&
                gamma(n+1)/gamma(2*i+3)
        enddo
        g(n,m) = g(n,m)+x

     enddo
  enddo

  !-------------------------------------------------------
  ! STAGE C: Compute special (negative) elements
  !-------------------------------------------------------

  lambda = lambda_in

  !--------------------------------
  ! Case 1: Special elements of g
  !--------------------------------

  if (lambda < 1/lambda_large) then

     ! Use asymptotic series for small lambda

     do n=2,m0
        m = 1-n
        c = 1.0_qp/(n+1)
        g(m,n) = c
        do i=1,16
           c = -c*(1+0.5_qp/i)*(2*i+n-1.0_qp)/(2*i+n+1.0_qp)*lambda
           g(m,n) = g(m,n)+c
        enddo
        g(m,n) = 2*g(m,n)*lambda**(0.5_qp*(n+1))

        do m=3-n,-1,2
           g(m,n) = 2.0_qp/(m+n)*(0.5_qp*(m-1)*g(m-2,n)+lambda**((n+1)/2.0_qp)/(1+lambda)**((m+n)/2.0_qp))
        enddo
     enddo

  else

     g(-1,0) = log((sqrt(1+lambda)+sqrt(lambda))/(sqrt(1.0_qp+lambda)-sqrt(lambda)))
     !g(-1,1) = log(1.0+lambda)

     m = -1
     do n=2,m0
        g(m,n) = 2.0_qp/(m+n)*(&
             0.5_qp*(n-1)*g(m,n-2)-lambda**((n-1)/2.0_qp)/(1+lambda)**((m+n)/2.0_qp))
     enddo

     do n=0,m0
        do m=0,-m0+2,-1
           if (m+n >= 1) then
              g(m-2,n) = 2.0_qp/(m-1.0_qp)*(&
                   0.5_qp*(m+n)*g(m,n)-lambda**((n+1)/2.0_qp)/(1+lambda)**((m+n)/2.0_qp))
           endif
        enddo
     enddo

  endif

  !--------------------------------
  ! Case 2: Special elements of gb
  !--------------------------------

  if (lambda > lambda_large) then

     ! Use small-lambda asymptotic series for g to get gb:
     ! 
     !           gb(n,m,lambda) = g(m,n,1/lambda)

     lambda = 1/lambda_in

     do n=2,m0

        m = 1-n
        c = 1.0_qp/(n+1)
        gb(n,m) = c
        do i=1,16
           c = -c*(1+0.5_qp/i)*(2*i+n-1.0_qp)/(2*i+n+1.0_qp)*lambda
           gb(n,m) = gb(n,m)+c
        enddo
        gb(n,m) = 2*gb(n,m)*lambda**(0.5_qp*(n+1))

        do m=3-n,-1,2
           gb(n,m) = 2.0_qp/(m+n)*(0.5_qp*(m-1)*gb(n,m-2)+lambda**((n+1)/2.0_qp)/(1+lambda)**((m+n)/2.0_qp))
        enddo

     enddo

     lambda = lambda_in

  else

     gb(0,-1) = log((sqrt(1.0_qp+lambda)+1.0_qp)/(sqrt(1.0_qp+lambda)-1.0_qp))
     !gb(1,-1) = log(1.0+1.0/lambda)

     n = -1
     do m=2,m0
        gb(m,n) = 2.0_qp/(m+n)*(&
             0.5_qp*(m-1)*gb(m-2,n)-lambda**((n+1)/2.0_qp)/(1+lambda)**((m+n)/2.0_qp))
     enddo

     do m=0,m0
        do n=0,-m0+2,-1

           if (m+n >= 1) then

              gb(m,n-2) = 2/(n-1.0_qp)*(&
                   0.5_qp*(m+n)*gb(m,n)-lambda**((n-1)/2.0_qp)/(1+lambda)**((m+n)/2.0_qp))

           endif

        enddo
     enddo

  endif

  !-------------------------------------
  ! CONCLUSION: Assign f,fb given g,gb:
  !-------------------------------------

  f(:,:)  = 0.0
  fb(:,:) = 0.0

  do m=-m0,m0
     do n=-m0,m0
        if (m+n >= -1) then
           f(m,n)  =  g(m,n)*0.25_qp*gamma(m+n+2)/lambda**((n+1)/2.0_qp)
           fb(m,n) = gb(m,n)*0.25_qp*gamma(m+n+2)/lambda**((n+1)/2.0_qp)
        endif
     enddo
  enddo

  deallocate(g)
  deallocate(gb)
  deallocate(gamma)
  deallocate(beta)

end subroutine neo_compute_fcoll

