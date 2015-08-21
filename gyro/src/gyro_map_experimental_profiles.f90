!-----------------------------------------------------------
! gyro_map_experimental_profiles.f90
!
! PURPOSE:
!  Generate radial profile functions on the simulation 
!  (slice) r-grid.
!
! NOTES:
! "_p" -> experimental (coarse) grid.
! "_s" -> simulation (fine) grid.
!
!  Nonuniform grid technique:
!  -------------------------
!
!  r_e -> equally spaced grid.
!  r   -> physical grid (may be nonuniform).    
!-----------------------------------------------------------

subroutine gyro_map_experimental_profiles

  use gyro_globals
  use gyro_profile_exp

  !-------------------------------------------------
  implicit none
  ! 
  real :: inter_fac
  !
  integer :: j1
  integer :: j2
  !-------------------------------------------------

  !-------------------------------------------------
  ! Grid mapping:
  !
  select case (nonuniform_grid_flag)

  case (0)

     ! UNIFORM

     r_s(:) = r_e(:)

     dr_eodr(:) = 1.0

  case (1)

     ! NONUNIFORM - changed endpoints

     r_s(ir_norm) = r_e(ir_norm)

     !   shoot outward r_s(i) to r_s(i+1)

     j1 = 0
     j2 = 0

     do i=ir_norm,n_x-1
        do i_exp=1,n_grid_exp-1
           if (r_s(i) >= r_p(i_exp) .and. r_s(i) < r_p(i_exp+1)) then
              j1 = i_exp
              j2 = i_exp+1
           endif
        enddo
        if (j1 == 0) then
           call catch_error("Error in map_experimental_profiles")
        endif

        inter_fac = (r_s(i)-r_p(j1))/(r_p(j2)-r_p(j1))

        rhosda_s(i) = rhosda_p(j1)+(rhosda_p(j2)-rhosda_p(j1))*inter_fac

        r_s(i+1) = r_s(i)+(r_e(i+1)-r_e(i))*(rhosda_s(i)/rhosda_s(ir_norm))

     enddo

     !   shoot inward r_s(i) to r_s(i-1)

     do i=ir_norm,2,-1
        do i_exp=1,n_grid_exp-1
           if (r_s(i) >= r_p(i_exp) .and. r_s(i) < r_p(i_exp+1)) then
              j1 = i_exp
              j2 = i_exp+1
           endif
        enddo
        inter_fac = (r_s(i)-r_p(j1))/(r_p(j2)-r_p(j1))

        rhosda_s(i) = rhosda_p(j1)+(rhosda_p(j2)-rhosda_p(j1))*inter_fac

        r_s(i-1) = r_s(i)+(r_e(i-1)-r_e(i))*(rhosda_s(i)/rhosda_s(ir_norm))

     enddo

     ! r_s unequal spaced grid found

     ! now reset r to r_s and compute dr_eodr

     r(:) = r_s(:)

     do i=1,n_x-1
        dr_eodr(i)=(r_e(i+1)-r_e(i))/(r(i+1)-r(i))
     enddo
     dr_eodr(n_x) = dr_eodr(n_x-1)

  end select
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get the GYRO slice 
  ! profiles (_s) from experimental (_p) profiles.
  ! 
  call cub_spline(r_p,rhogrid_exp,n_grid_exp,r_s,rhogrid_s,n_x)

  call cub_spline(r_p,rmaj_exp/rmin_exp(n_grid_exp),n_grid_exp,r_s,rmaj_s,n_x)
  call cub_spline(r_p,drmaj_p,n_grid_exp,r_s,drmaj_s,n_x)
  call cub_spline(r_p,q_p,n_grid_exp,r_s,q_s,n_x)
  call cub_spline(r_p,shat_p,n_grid_exp,r_s,shat_s,n_x)
  call cub_spline(r_p,kappa_exp,n_grid_exp,r_s,kappa_s,n_x)
  call cub_spline(r_p,s_kappa_p,n_grid_exp,r_s,s_kappa_s,n_x)
  call cub_spline(r_p,delta_exp,n_grid_exp,r_s,delta_s,n_x)
  call cub_spline(r_p,s_delta_p,n_grid_exp,r_s,s_delta_s,n_x)
  call cub_spline(r_p,zeta_exp,n_grid_exp,r_s,zeta_s,n_x)
  call cub_spline(r_p,s_zeta_p,n_grid_exp,r_s,s_zeta_s,n_x)
  call cub_spline(r_p,zmag_exp/rmin_exp(n_grid_exp),n_grid_exp,r_s,zmag_s,n_x)
  call cub_spline(r_p,dzmag_p,n_grid_exp,r_s,dzmag_s,n_x)

  do is=1,n_spec
     call cub_spline(r_p,den_exp(is,:),n_grid_exp,r_s,den_s(is,:),n_x)
     call cub_spline(r_p,tem_exp(is,:),n_grid_exp,r_s,tem_s(is,:),n_x)
  enddo
  call cub_spline(r_p,ptot_exp,n_grid_exp,r_s,ptot_s,n_x)

  call cub_spline(r_p,z_eff_exp,n_grid_exp,r_s,z_eff_s,n_x)
  call cub_spline(r_p,b_unit_p,n_grid_exp,r_s,b_unit_s,n_x)
  call cub_spline(r_p,rhosda_p,n_grid_exp,r_s,rhosda_s,n_x)
  call cub_spline(r_p,csda_p,n_grid_exp,r_s,csda_s,n_x)

  do is=1,n_spec
     call cub_spline(r_p,dlntdr_p(is,:),n_grid_exp,r_s,dlntdr_s(is,:),n_x)
     call cub_spline(r_p,dlnndr_p(is,:),n_grid_exp,r_s,dlnndr_s(is,:),n_x)
  enddo
  call cub_spline(r_p,dlnptotdr_p,n_grid_exp,r_s,dlnptotdr_s,n_x)

  call cub_spline(r_p,gamma_e_exp,n_grid_exp,r_s,gamma_e_s,n_x)
  call cub_spline(r_p,gamma_p_exp,n_grid_exp,r_s,gamma_p_s,n_x)
  call cub_spline(r_p,mach_exp,n_grid_exp,r_s,mach_s,n_x)
  call cub_spline(r_p,w0_exp,n_grid_exp,r_s,w0_s,n_x)
  call cub_spline(r_p,w0p_exp,n_grid_exp,r_s,w0p_s,n_x)

  call cub_spline(r_p,beta_unit_p,n_grid_exp,r_s,beta_unit_s,n_x)
  call cub_spline(r_p,beta_unit_ptot_p,n_grid_exp,r_s,beta_unit_ptot_s,n_x)

  if (num_equil_flag == 1) then
     do j1=1,8
        do j2=0,n_fourier_geo
           call cub_spline(r_p,geo_p(j1,j2,:),n_grid_exp,r_s,a_fourier_geo_s(j1,j2,:),n_x)
        enddo
     enddo
  else
     a_fourier_geo_s(:,:,:) = 0.0
  endif
  !------------------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[map_exp_profiles done]'
  endif

end subroutine gyro_map_experimental_profiles
