!--------------------------------------------------------------
! gyro_read_input.f90
!
! PURPOSE:
!  Complete read of input parameters.
!
! NOTES:
!  Input parameters are automatically broadcast 
!  to entire GYRO_COMM_WORLD. 
!
!   =============================================
!   HOW TO ADD INPUT VARIABLES 
!
!   Adding a new input variable is now simple.
!   Just append a line analogous to one of those 
!   below.  The variable will be read from the 
!   file 'input.dat' and broadcast to all members 
!   of GYRO_COMM_WORLD.
!   =============================================
!---------------------------------------------------------------

subroutine gyro_read_input

  use gyro_globals

  !-------------------------
  implicit none
  !
  real :: x
  !-------------------------

  if (i_proc == 0) open(unit=1,file=trim(path)//'input.gyro.gen',status='old')

  !--------------------------------------------------------
  ! BEGIN reading data:
  !
  call readbc_int(n_x)
  call readbc_int(n_theta_section)
  call readbc_int(n_pass)
  call readbc_int(n_trap)
  call readbc_int(n_energy)
  call readbc_int(n_blend)
  call readbc_int(blend_fit_order)
  call readbc_int(n_theta_plot)
  call readbc_int(n0)
  call readbc_int(n_n)
  call readbc_int(d_n)
  call readbc_int(n_ref)
  call readbc_real(q0)
  call readbc_real(s0)
  call readbc_real(delta0)
  call readbc_real(zeta0)
  call readbc_real(s_delta0)
  call readbc_real(s_zeta0)
  call readbc_real(kappa0)
  call readbc_real(s_kappa0)
  call readbc_real(drmaj0)
  call readbc_real(r_maj)
  call readbc_real(r0)
  call readbc_real(box_multiplier)
  call readbc_real(rho_star)
  call readbc_real(dt)
  call readbc_real(time_max)
  call readbc_int(time_skip)
  call readbc_int(boundary_method)
  call readbc_int(m_dx)
  call readbc_int(m_gyro)
  call readbc_int(n_explicit_damp)
  call readbc_real(explicit_damp)
  call readbc_real(explicit_damp_elec)
  call readbc_int(debug_flag)
  call readbc_int(nonlinear_flag)
  call readbc_real(amp_n)
  call readbc_real(amp_0)
  call readbc_real(radial_upwind)
  call readbc_int(electron_method)
  call readbc_int(radial_profile_method)
  call readbc_int(plot_u_flag)
  call readbc_int(plot_epar_flag)
  call readbc_int(plot_n_flag)
  call readbc_int(plot_e_flag)
  call readbc_int(plot_v_flag)
  call readbc_real(z_eff)
  call readbc_real(betae_unit)
  call readbc_real(ampere_scale)
  call readbc_int(n_field)
  call readbc_int(source_flag)
  call readbc_real(nu_source)
  call readbc_int(verbose_flag)
  call readbc_int(nonuniform_grid_flag)
  call readbc_real(s_grid)
  call readbc_real(gamma_e)
  call readbc_real(nu_ei)
  call readbc_real(nu_ei_scale)
  call readbc_real(nu_ii_scale)
  call readbc_real(nu_i_krook)
  call readbc_real(plot_filter)
  call readbc_int(n_source)
  call readbc_int(flat_profile_flag)
  call readbc_int(density_method)
  call readbc_int(integrator_method)
  call readbc_real(energy_max)

  ! Species: (1-3) ions and electrons.

  mu_vec = 0.0
  call readbc_real(x) ; mu_vec(1) = x
  call readbc_real(x) ; mu_vec(2) = x
  call readbc_real(x) ; mu_vec(3) = x
  call readbc_real(x) ; mu_vec(4) = x
  call readbc_real(x) ; mu_vec(5) = x
  call readbc_real(x) ; mu_vec(0) = x

  dlnndr_vec = 0.0
  call readbc_real(x) ; dlnndr_vec(1) = x
  call readbc_real(x) ; dlnndr_vec(2) = x
  call readbc_real(x) ; dlnndr_vec(3) = x
  call readbc_real(x) ; dlnndr_vec(4) = x
  call readbc_real(x) ; dlnndr_vec(5) = x
  call readbc_real(x) ; dlnndr_vec(0) = x

  dlntdr_vec = 0.0
  call readbc_real(x) ; dlntdr_vec(1) = x
  call readbc_real(x) ; dlntdr_vec(2) = x
  call readbc_real(x) ; dlntdr_vec(3) = x
  call readbc_real(x) ; dlntdr_vec(4) = x
  call readbc_real(x) ; dlntdr_vec(5) = x
  call readbc_real(x) ; dlntdr_vec(0) = x

  n_vec = 0.0
  call readbc_real(x) ; n_vec(1) = x
  call readbc_real(x) ; n_vec(2) = x
  call readbc_real(x) ; n_vec(3) = x
  call readbc_real(x) ; n_vec(4) = x
  call readbc_real(x) ; n_vec(5) = x
  n_vec(0) = 1.0  

  t_vec = 0.0
  call readbc_real(x) ; t_vec(1) = x
  call readbc_real(x) ; t_vec(2) = x
  call readbc_real(x) ; t_vec(3) = x
  call readbc_real(x) ; t_vec(4) = x
  call readbc_real(x) ; t_vec(5) = x
  t_vec(0) = 1.0

  eps_dlnndr_vec = 0.0
  call readbc_real(x) ; eps_dlnndr_vec(1) = x
  call readbc_real(x) ; eps_dlnndr_vec(2) = x
  call readbc_real(x) ; eps_dlnndr_vec(3) = x
  call readbc_real(x) ; eps_dlnndr_vec(4) = x
  call readbc_real(x) ; eps_dlnndr_vec(5) = x
  call readbc_real(x) ; eps_dlnndr_vec(0) = x

  eps_dlntdr_vec = 0.0
  call readbc_real(x) ; eps_dlntdr_vec(1) = x
  call readbc_real(x) ; eps_dlntdr_vec(2) = x
  call readbc_real(x) ; eps_dlntdr_vec(3) = x
  call readbc_real(x) ; eps_dlntdr_vec(4) = x
  call readbc_real(x) ; eps_dlntdr_vec(5) = x
  call readbc_real(x) ; eps_dlntdr_vec(0) = x

  z_vec = 0.0
  call readbc_real(x) ; z_vec(1) = x
  call readbc_real(x) ; z_vec(2) = x
  call readbc_real(x) ; z_vec(3) = x
  call readbc_real(x) ; z_vec(4) = x
  call readbc_real(x) ; z_vec(5) = x
  z_vec(0) = -1.0

  orbit_upwind_vec = 0.0
  call readbc_real(x) ; orbit_upwind_vec(1) = x
  call readbc_real(x) ; orbit_upwind_vec(2) = x
  call readbc_real(x) ; orbit_upwind_vec(3) = x
  call readbc_real(x) ; orbit_upwind_vec(4) = x
  call readbc_real(x) ; orbit_upwind_vec(5) = x
  call readbc_real(x) ; orbit_upwind_vec(0) = x

  ! More profile parameters
  call readbc_real(pgamma0)
  call readbc_real(pgamma0_scale)
  call readbc_real(mach0)
  call readbc_real(mach0_scale)
  call readbc_int(lindiff_method)
  call readbc_int(trapdiff_flag)
  call readbc_int(restart_new_flag)
  call readbc_int(restart_data_skip)
  call readbc_int(kill_i_parallel_flag)
  call readbc_int(kill_i_drift_flag)
  call readbc_int(kill_e_drift_flag)
  call readbc_int(kill_coll_flag)
  call readbc_real(doppler_scale)
  call readbc_int(nl_method)
  call readbc_int(kill_gyro_b_flag)
  call readbc_int(velocity_output_flag)
  call readbc_real(q_scale)
  call readbc_int(dist_print)
  call readbc_int(nint_ORB_s)
  call readbc_int(nint_ORB_do)
  call readbc_int(udsymmetry_flag)
  call readbc_int(gyro_method)
  call readbc_int(sparse_method)
  call readbc_int(n_mumps_max)
  call readbc_int(n_study)
  call readbc_real(amp_study)
  call readbc_real(lambda_debye_scale)
  call readbc_real(lambda_debye)
  call readbc_int(n_x_offset)
  call readbc_int(n_theta_mult)
  call readbc_int(silent_flag)
  call readbc_int(nonlinear_transfer_flag)
  call readbc_real(l_x)
  call readbc_real(l_y)
  call readbc_int(entropy_flag)
  call readbc_int(ord_rbf)
  call readbc_int(num_equil_flag)
  call readbc_real(zmag0)
  call readbc_real(dzmag0)
  call readbc_int(output_flag)
  call readbc_real(ipccw)
  call readbc_real(btccw)
  call readbc_int(geo_gradbcurv_flag)
  call readbc_int(geo_fastionbeta_flag)
  call readbc_real(geo_betaprime_scale)
  call readbc_int(poisson_z_eff_flag)
  call readbc_int(z_eff_method)
  call readbc_int(truncation_method)
  call readbc_real(fluxaverage_window)
  call readbc_int(huge_restart)

  ! GK eigenvalue solver inputs
  call readbc_int(gkeigen_proc_mult)
  call readbc_int(gkeigen_method)
  call readbc_int(gkeigen_matrixonly)
  call readbc_int(gkeigen_mwrite_flag)
  call readbc_int(gkeigen_kspace_dim)
  call readbc_int(gkeigen_n_values)
  call readbc_int(gkeigen_iter)
  call readbc_real(gkeigen_tol)
  call readbc_real(gkeigen_omega_target)
  call readbc_real(gkeigen_gamma_target)
  call readbc_int(linsolve_method)

  ! Field eigenvalue solver inputs
  call readbc_int(fieldeigen_root_method)
  call readbc_real(fieldeigen_wr)
  call readbc_real(fieldeigen_wi)
  call readbc_real(fieldeigen_tol)

  call readbc_int(reintegrate_flag)
  call readbc_int(ic_method)
  call readbc_int(zf_test_flag)
  call readbc_int(lock_ti_flag)
  ! DONE reading data.
  !--------------------------------------------------------

  if (i_proc == 0) close(1)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_read_input done]'
  endif

end subroutine gyro_read_input

!------------------------------------------------------------
! Service routines: 
!
! (1) read and broadcast an integer:
!
subroutine readbc_int(p)

  use mpi
  use gyro_globals

  implicit none
  integer, intent(inout) :: p

  if (i_proc == 0) read(1,*) p

  call MPI_BCAST(p,1,MPI_INTEGER,0,GYRO_COMM_WORLD,i_err)

end subroutine readbc_int
!
! (2) read and broadcast a real:
!
subroutine readbc_real(x)
  
  use mpi
  use gyro_globals

  implicit none
  real, intent(inout) :: x

  if (i_proc == 0) read(1,*) x

  call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,i_err)

end subroutine readbc_real
!------------------------------------------------------------
