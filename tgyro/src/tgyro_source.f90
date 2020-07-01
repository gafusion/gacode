!------------------------------------------------------------
! tgyro_source.f90
!
! PURPOSE:
!  Calculation of primitive source elements, and volume 
!  integration into required form.
!------------------------------------------------------------

subroutine tgyro_source 

  use tgyro_globals

  implicit none

  integer :: i,i_ion
  real, external :: sigv
  real, external :: dtrate_dv
  real :: n_d,n_t
  real :: s_alpha
  real :: g,phi,wpe,wce
  real, parameter :: r_coeff=0.8

  !-------------------------------------------------------
  ! Source terms (erg/cm^3/s):
  !
  do i=1,n_r

     if (loc_scenario > 2) then
        !-------------------------------------------------------
        ! Alpha power
        !  - sigv in cm^3/s
        if (tgyro_dt_method == 1) then
           ! Assume D and T given by ion 1 and ion 2 
           ! (order doesn't matter)
           n_d = ni(1,i)
           n_t = ni(2,i)
        else
           ! Assume ion 1 is DT hybrid.
           n_d = 0.5*ni(1,i)
           n_t = 0.5*ni(1,i)
        endif
        ! Alpha particle source and power 
        ! - Can use 'hively' or 'bosch' formulae.
        sn_alpha(i) = n_d*n_t*sigv(ti(1,i)/1e3,'bosch') * tgyro_input_fusion_scale

        s_alpha      = sn_alpha(i)*e_alpha
        s_alpha_i(i) = s_alpha*frac_ai(i)
        s_alpha_e(i) = s_alpha*frac_ae(i)
     else
        sn_alpha(i)  = 0.0
        s_alpha_i(i) = 0.0
        s_alpha_e(i) = 0.0
     endif
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Bremsstrahlung radiation
     ! - From NRL formulary 
     ! - 1 W/cm^3 = 1e7 erg/cm^3/s

     s_brem(i) = 1e7*1.69e-32*ne(i)**2*sqrt(te(i))*z_eff(i)
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Synchrotron radiation
     ! - Trubnikov, JETP Lett. 16 (1972) 25.
     wpe = sqrt(4*pi*ne(i)*e**2/me)
     wce = e*abs(b_ref(i))/(me*c)
     g   = k*te(i)/(me*c**2)
     phi = 60*g**1.5*sqrt((1.0-r_coeff)*(1+1/aspect_rat/sqrt(g))/(r_min*wpe**2/c/wce))

     s_sync(i) = me/(3*pi*c)*g*(wpe*wce)**2*phi
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Classical electron-ion energy exchange
     ! - Positive as defined on RHS of ion equation
     ! - Multiply formulary expression by (3/2)ne:

     s_exch(i) = 1.5*nu_exch(i)*ne(i)*k*(te(i)-ti(1,i))
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Anomalous electron-ion energy exchange
     ! - Positive as defined on RHS of ion equation
     ! - Skip exchange with fast ions

     s_expwd(i) = 0.0
     do i_ion=1,loc_n_ion
        if (therm_flag(i_ion) == 1) then
           s_expwd(i) = s_expwd(i)+expwd_i_tur(i_ion,i)*s_gb(i)
        endif
     enddo
     !-------------------------------------------------------

  enddo
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Powers in units of erg/s

  ! Get integrated alpha-power
  call tgyro_volume_int(s_alpha_i,p_i_fus)
  call tgyro_volume_int(s_alpha_e,p_e_fus)
  call tgyro_volume_int(sn_alpha,f_he_fus)

  ! Get integrated collisional exchange power
  call tgyro_volume_int(s_exch,p_exch)
  
  ! Get integrated anomalous exchange power
  call tgyro_volume_int(s_expwd,p_expwd)

  ! Get integrated Bremsstrahlung power
  call tgyro_volume_int(s_brem,p_brem)

  ! Get integrated Synchrotron power
  call tgyro_volume_int(s_sync,p_sync)
  !-------------------------------------------------------

  !-------------------------------------------------------
  select case (loc_scenario)

  case (1)

     ! Experimental, static exchange

     p_i(:) = p_i_in(:) ! Total ion input power (including exchange)
     p_e(:) = p_e_in(:) ! Total electon input power (including exchange)

  case (2)

     ! Experimental, dynamic exchange

     p_i(:) = p_i_in(:) &              ! Total ion input power 
          +(p_exch(:)-p_exch_in(:)) &  ! Consistent e-i exchange
          +p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

     p_e(:) = p_e_in(:) &              ! Total electron input power 
          -(p_exch(:)-p_exch_in(:)) &  ! Consistent e-i exchange
          -p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange
     
  case (3)

     ! Reactor with consistent alpha power, exchange and radiation.

     p_i(:) = &
          +p_i_fus(:) &                ! Fusion power to ions
          +p_i_aux_in(:) &             ! Auxiliary ion heating [fixed]
          +p_exch(:) &                 ! Collisional exchange
          +p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

     p_e(:) = &
          +p_e_fus(:) &                ! Fusion power to electrons
          +p_e_aux_in(:) &             ! Auxiliary electron heating [fixed]
          -p_exch(:)   &               ! Collisional exchange
          -p_brem(:) &                 ! Bremsstrahlung radiation
          -p_sync(:) &                 ! Synchrotron radiation
          -p_line_in(:) &              ! Line radiation [fixed] 
          -p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

  end select
  !-------------------------------------------------------

  !------------------------------------------------
  ! Target energy fluxes in erg/s/cm^2 
  eflux_i_target(1) = 0.0 
  eflux_e_target(1) = 0.0

  eflux_i_target(2:n_r) = p_i(2:n_r)/volp(2:n_r) 
  eflux_e_target(2:n_r) = p_e(2:n_r)/volp(2:n_r)
  !------------------------------------------------

  !------------------------------------------------
  ! Target particle fluxes in 1/s/cm^2
  select case (loc_pflux_method)

  case (1)

     pflux_e_target(:) = 0.0

  case (2)

     pflux_e_target(1) = 0.0
     pflux_e_target(2:n_r) = f_b_in(2:n_r)/volp(2:n_r)

  case (3)

     pflux_e_target(1) = 0.0
     pflux_e_target(2:n_r) = (f_b_in(2:n_r) + f_w_in(2:n_r)) &
          /volp(2:n_r) 

  end select
  !------------------------------------------------

  !------------------------------------------------
  ! Target angular momentum fluxes in erg/cm^2
  ! 
  mflux_target(1) = 0.0 
  mflux_target(2:n_r) = mf_in(2:n_r)/volp(2:n_r)
  !------------------------------------------------

  !------------------------------------------------
  ! Target He ash flux in 1/s/cm^2
  !
  pflux_he_target(1) = 0.0
  pflux_he_target(2:n_r) = f_he_fus(2:n_r)/volp(2:n_r)
  !------------------------------------------------

  !------------------------------------------------
  ! Target fluxes in GB units
  eflux_i_target(:)  = eflux_i_target(:)/q_gb(:)
  eflux_e_target(:)  = eflux_e_target(:)/q_gb(:)
  pflux_e_target(:)  = pflux_e_target(:)/gamma_gb(:)
  mflux_target(:)    = mflux_target(:)/pi_gb(:)
  pflux_he_target(:) = pflux_he_target(:)/gamma_gb(:)
  !------------------------------------------------

end subroutine tgyro_source

