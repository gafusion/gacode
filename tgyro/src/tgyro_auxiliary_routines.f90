!-----------------------------------------------------------
! sigv.f90
!
! PURPOSE:
!  Compute D-T fusion reactivity <sigma*v> in cm^3/s.
!----------------------------------------------------------

real function sigv(ti,type)

  implicit none

  real, intent(in) :: ti
  character(len=*), intent(in) :: type

  real :: c1,c2,c3,c4,c5,c6,c7
  real :: r0,theta,xi,bg,er

  select case (trim(type))

  case ('hively')

     ! L.M. Hively, Nucl. Fusion 17 (1977) 873.

     ! Table 1: (S5)
     r0  = 0.2935
     c1 = -21.377692
     c2 = -25.204054
     c3 = -7.1013427e-2
     c4 = 1.9375451e-4
     c5 = 4.9246592e-6
     c6 = -3.9836572e-8

     ! Eq. (5)
     sigv = exp(c1/ti**r0+c2+ti*(c3+ti*(c4+ti*(c5+ti*c6))))

  case ('bosch')

     ! H.-S. Bosch and G.M. Hale, Nucl. Fusion 32 (1992) 611.

     ! Table VII:
     c1 = 1.17302e-9
     c2 = 1.51361e-2
     c3 = 7.51886e-2
     c4 = 4.60643e-3
     c5 = 1.3500e-2
     c6 = -1.06750e-4
     c7 = 1.36600e-5
     bg = 34.3827
     er = 1.124656e6

     ! Eq. (12) 
     r0    = ti*(c2+ti*(c4+ti*c6))/(1.0+ti*(c3+ti*(c5+ti*c7)))
     theta = ti/(1.0-r0)
     xi    = (bg**2/(4.0*theta))**(1.0/3.0)

     sigv = c1*theta*sqrt(xi/(er*ti**3))*exp(-3.0*xi)

  case default

     print '(a)','ERROR: (sigv) Unknown reactivity type.'
     sigv = 0.0
     stop

  end select

end function sigv

!--------------------------------------------------------
! sivukhin.f90
!
! PURPOSE:
!  Compute a low-accuracy but fast approximation to the 
!  ion-alpha heating fraction.
!
!               x
!            1  /     dy
!    F(x) = --- | -----------
!            x  /  1+y^(3/2)
!               0
!          
!  Here, F is the fraction of the alpha energy transferred
!  to ions (at common temperature Ti) by collisions, and
!
!               x = E_alpha/E_crit  
!
!  Details are given in Stix, Plasma Phys. 14 (1972) 367.
!  The function F is derived from Sivukhin's energy loss
!  equation and so that is the rationale for the name.
!---------------------------------------------------------

real function sivukhin(x)

  use math_constants

  implicit none

  real, intent(in) :: x
  integer :: i
  real :: f,dy,yi

  integer, parameter :: n=12


  if (x > 0.1) then

     if (x > 4.0) then

        ! Large-x asymptotic formula

        f = (2*pi/3)/sin(2*pi/3)-2.0/sqrt(x)+0.5/(x*x)

     else

        ! Numerical integration

        dy = x/(n-1)
        f  = 0.0
        do i=1,n
           yi = (i-1)*dy
           if (i == 1 .or. i == n) then
              f = f+0.5/(1.0+yi**1.5) 
           else
              f = f+1.0/(1.0+yi**1.5) 
           endif
        enddo
        f = f*dy

     endif

     sivukhin = f/x

  else

     ! Small-x asymptotic series

     sivukhin = 1.0-0.4*x**1.5

  endif

end function sivukhin

subroutine rad_alpha(ne,ni,te,ti,s_alpha_he,s_alpha_i,s_alpha_e,frac_ai,e_cross,n_alpha,t_alpha,n,nion)

  use tgyro_globals, only : &
       pi,&
       me,&
       mi,&
       malpha,&
       therm_flag,&
       zi_vec,&
       e_alpha,&
       k,&
       dt_flag,&
       tgyro_input_fusion_scale

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: nion
  real, intent(in) :: ni(nion,n)
  real, intent(in) :: ne(n)
  real, intent(in) :: ti(nion,n)
  real, intent(in) :: te(n)
  real, intent(inout) :: s_alpha_he(n)
  real, intent(inout) :: s_alpha_i(n)
  real, intent(inout) :: s_alpha_e(n)
  real, intent(inout) :: frac_ai(n)
  real, intent(inout) :: e_cross(n)
  real, intent(inout) :: n_alpha(n)
  real, intent(inout) :: t_alpha(n)

  real, external :: sivukhin
  real, external :: sigv

  real :: x_a
  real :: a,i2,i4,taus
  real :: n_d,n_t
  real :: s_alpha
  real, dimension(:), allocatable :: c_a

  integer :: i

  allocate(c_a(n))

  ! Alpha heating coefficients [Stix, Plasma Phys. 14 (1972) 367] 
  ! See in particular Eqs. 15 and 17.
  c_a(:) = 0.0
  do i=1,nion
     if (therm_flag(i) == 1) then
        c_a(:) = c_a(:)+(ni(i,:)/ne(:))*zi_vec(i)**2/(mi(i)/malpha)
     endif
  enddo

  n_d = 0.0
  n_t = 0.0

  do i=1,n
     
     if (dt_flag == 1) then
        ! D and T given by ion 1 and ion 2 (order doesn't matter)
        n_d = ni(1,i)
        n_t = ni(2,i)
     else if (dt_flag == 2) then
        ! D and T given by half of ion 1
        n_d = ni(1,i) / 2
        n_t = ni(1,i) / 2
     endif

     ! e_cross = (1/2) malpha v_cross^2 [agrees with Estrada 2006]
     e_cross(i) = k*te(i)*(4*sqrt(me/malpha)/(3*sqrt(pi)*c_a(i)))**(-2.0/3.0)
     ! x_a = (v_alpha/v_cross)^2
     x_a = e_alpha/e_cross(i)
     frac_ai(i) = sivukhin(x_a)

     ! Alpha particle source and power 
     !  - Can use 'hively' or 'bosch' formulae.
     !  - sigv in cm^3/s
     s_alpha_he(i) = n_d*n_t*sigv(ti(1,i)/1e3,'bosch') * tgyro_input_fusion_scale

     s_alpha      = s_alpha_he(i)*e_alpha
     s_alpha_i(i) = s_alpha*frac_ai(i)
     s_alpha_e(i) = s_alpha*(1-frac_ai(i))

     ! JC 2022:
     ! Adding alpha-particle density/temperature calculation [Estrada 2006]
     ! Eqs (9),(10), where a=v_cross/v_alpha
     a = 1/sqrt(x_a)
     i2 = (1/3.0)*log((1+a**3)/a**3)
     i4 = 0.5-a**2*((1/6.0)*log((1-a+a**2)/(1+a)**2)+ &
          1/sqrt(3.0)*(atan((2-a)/(a*sqrt(3.0)))+pi/6))

     taus = 0.0

     ! Eqs (7),(15)
     !n_alpha(i) = s_alpha_he(i)*taus*i2
     !t_alpha(i) = 2*i4/(3*i2)*e_alpha

     n_alpha(i) = 1.0 ; t_alpha(i) = 1.0
     
  enddo

  deallocate(c_a)

end subroutine rad_alpha

subroutine collision_rates(ne,ni,te,ti,nui,nue,nu_exch,taus,n,nion)

  use tgyro_globals, only : &
       pi,&
       me,&
       mi,&
       malpha,&
       e, &
       k, &
       therm_flag,&
       zi_vec

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: nion
  real, intent(in) :: ni(nion,n)
  real, intent(in) :: ne(n)
  real, intent(in) :: ti(nion,n)
  real, intent(in) :: te(n)
  real, intent(inout) :: nui(nion,n)
  real, intent(inout) :: nue(n)
  real, intent(inout) :: nu_exch(n)
  real, intent(inout) :: taus(n)

  real, dimension(:), allocatable :: loglam
  real :: c_exch
  real :: z_a
  
  integer :: i

  allocate(loglam(n))
  
  ! Coulomb logarithm
  loglam(:) = 24.0-log(sqrt(ne(:))/te(:))

  ! 1/tau_ii (Belli 2008) in 1/s
  do i=1,nion
     nui(i,:) = sqrt(2.0)*pi*ni(i,:)*(zi_vec(i)*e)**4*loglam(:) &
          /(sqrt(mi(i))*(k*ti(i,:))**1.5)
  enddo

  ! 1/tau_ee (Belli 2008) in 1/s
  nue(:) = sqrt(2.0)*pi*ne(:)*e**4*loglam(:)/(sqrt(me)*(k*te(:))**1.5)

  ! NOTE: 
  ! c_exch = 1.8e-19 is the formulary exch. coefficient
  c_exch = 2.0*(4.0/3)*sqrt(2.0*pi)*e**4/k**1.5

  ! nu_exch in 1/s
  nu_exch(:) = 0.0
  do i=1,nion
     if (therm_flag(i) == 1) then
        nu_exch(:) = nu_exch(:)+c_exch*sqrt(me*mi(i))*zi_vec(i)**2 &
             *ni(i,:)*loglam(:)/(me*ti(i,:)+mi(i)*te(:))**1.5
     endif
  enddo

  ! JC: 2022
  ! Alpha slowing-down time in terms of tau_ee = 1/nue [need to double-check)
  z_a = 2.0
  taus(:) = (3/4.0)*sqrt(pi)*(malpha/me)*1/(z_a**2*nue(:))
  
  deallocate(loglam)

end subroutine collision_rates
