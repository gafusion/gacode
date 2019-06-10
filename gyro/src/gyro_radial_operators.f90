!-------------------------------------------------------
! gyro_radial_operators.f90
!
! PURPOSE:
!  This routine calculates complex forward and inverse
!  FT coefficients and finite-element coefficients,
!  derivative operators, and is the top-level driver
!  for field-matrix and gyroaverage construction.
!-----------------------------------------------------

subroutine gyro_radial_operators

  use gyro_globals
  use math_constants

  !--------------------------------------
  implicit none  
  !
  integer :: p
  integer :: n_node
  integer :: i_lump
  integer :: j_lump
  !
  real :: t0
  real :: denom
  complex :: z0
  real, dimension(:), allocatable :: r_node
  real, dimension(-m_dx:m_dx-i_dx) :: w_temp
  real, dimension(-m_dx:m_dx) :: s_temp
  !
  real, external :: BLEND_f3
  !--------------------------------------

  !---------------------------------------------------
  ! COMPLEX FOURIER COEFFICIENTS:
  !
  do i=1,n_x
     p = i-n_x/2-1
     do ip=1,n_x
        m = ip-n_x/2-1

        ! c's are the matrix elements for m -> p

        cr(i,ip) = exp(i_c*m*p*(pi_2/n_x))

     enddo
  enddo

  ! ci's are the matrix elements for the inverse 
  ! transformation, namely p -> m.

  cri = conjg(cr)/n_x
  !---------------------------------------------------

  !---------------------------------------------------
  ! SOURCE EXPANSION COEFFICIENTS:
  !
  if (source_method > 1) then

     n_node = n_lump+3

     allocate(r_node(n_node))

     ! Node points are equally-spaced in r.

     do i=1,n_node
        r_node(i) = r(1)+(r(n_x)-r(1))*(i-1)/(n_node-1.0)
     enddo

     do i_lump=1,n_lump
        do i=1,n_x
           t0 = (r(i)-r_node(1))/(r_node(2)-r_node(1))
           b_src(i,i_lump) = BLEND_f3(i_lump-1,t0)
        enddo
     enddo

     m_src(:,:) = 0.0
     do i_lump=1,n_lump
        do j_lump=1,n_lump
           do i=1,n_x
              m_src(i_lump,j_lump) = m_src(i_lump,j_lump)+ & 
                   b_src(i,i_lump)*b_src(i,j_lump)
           enddo
        enddo
     enddo

     call DGETRF(n_lump,n_lump,m_src,n_lump,src_piv,info)

     deallocate(r_node)

  endif
  !
  ! Set source amplitude coefficient to zero if n /= 0.
  if (n_1(in_1) /= 0) nu_source = 0.0
  !---------------------------------------------------

  !---------------------------------------------------
  ! Radial index mapping
  !
  if (boundary_method == 1) then

     do i=1,n_x
        i_cyc(i+n_x) = i
        i_cyc(i) = i
        i_cyc(i-n_x) = i
     enddo

     i_loop(:) = i_cyc(:)

  else

     ! These forms will be used in gyro_blend_*:
  
     do i=1,n_x
        i_cyc(i+n_x) = n_x
        i_cyc(i) = i
        i_cyc(i-n_x) = 1
     enddo

     do i=1,n_x
        i_loop(i+n_x) = i+n_x
        i_loop(i) = i
        i_loop(i-n_x) = i-n_x
     enddo

     ! Arrays for efficient implementation of buffer damping
     explicit_damp_vec(:,:) = 0.0
     do i=1,n_x
        if (n_1(in_1) == 0) then
           if (i <= n_explicit_damp .or. i > n_x-n_explicit_damp) then
              explicit_damp_vec(:,i)      = explicit_damp
              explicit_damp_vec(indx_e,i) = explicit_damp
           endif
        endif
     enddo

  endif
  !---------------------------------------------------

  ! Define delta-function (for use with d/dr)

  w_d0(:) = (0.0,0.0)
  w_d0(0) = (1.0,0.0)

  ! Define delta-function (for use with << >>)

  w_g0(:) = (0.0,0.0)
  w_g0(0) = (1.0,0.0)

  ! Variable width delta-function

  w_gd0(:) = (0.0,0.0)
  w_gd0(0) = (1.0,0.0)

  !-----------------------------------
  ! Use PolyDiff to get 1st-derivative
  ! using 2*m_dx+1 points:
  !
  if (i_dx == 0) then

     ! Finite-difference derivative and dissipation

     call polydiff(2*m_dx,m_dx,w_temp,denom)
     w_d1 = w_temp/d_x

     call pascal(2*m_dx+1,s_temp)
     s_d1 = abs(radial_upwind)*s_temp*abs(w_temp(m_dx)/d_x)

     call poly2diff(m_dx,w_temp)
     w_d2 = w_temp/(d_x*d_x)

  else

     ! Pseudo-spectral derivative (FD dissipation)

     w_d1(:) = (0.0,0.0)

     do m=-m_dx,m_dx-1
        do p=-m_dx,m_dx-1
           z0 = exp(i_c*p*m*(pi_2/n_x))
           w_d1(m) = w_d1(m)+p*z0
        enddo
     enddo

     w_d1(:) = -(1.0/n_x)*(pi_2*i_c/x_length)*w_d1(:)

     w_d2(:) = (0.0,0.0)

     do m=-m_dx,m_dx-1
        do p=-m_dx,m_dx-1
           z0 = exp(i_c*p*m*(pi_2/n_x))
           w_d2(m) = w_d2(m)+p*p*z0
        enddo
     enddo

     ! Note the sign change
     w_d2(:) = (1.0/n_x)*(pi_2*i_c/x_length)**2*w_d2(:)

     if (m_dx == 2) then

        s_d1(-2) = -2.0
        s_d1(-1) = 4.0
        s_d1(0)  = -6.0
        s_d1(1)  = 4.0 

        ! Usual 3rd order upwind normalization
        s_d1(-2:1) = abs(radial_upwind)*s_d1(-2:1)/(12.0*d_x)

     else

        w_temp = 0.0
        s_temp = 0.0

        call polydiff(2*m_dx-2,m_dx,w_temp(-m_dx+1:m_dx-1),denom)
        call pascal(2*m_dx-1,s_temp(-m_dx+1:m_dx-1))

        s_d1 = abs(radial_upwind)*s_temp*abs(w_temp(m_dx-1)/d_x)

     endif

  endif
  !-----------------------------------

  !-----------------------------------
  ! Find tau-space gyroaverage:
  !
  call gyro_bessel_stencils
  !-----------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_radial_operators done]'
  endif

end subroutine gyro_radial_operators
