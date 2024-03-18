!------------------------------------------------------------
! qlgyro_run.f90
!
! PURPOSE:
!  Driver for QLGYRO
!  Runs GYRO for all the KYs specified
!   
!
! 
!------------------------------------------------------------

subroutine qlgyro_run_gyro

  use mpi
  use qlgyro_globals
  use gyro_interface
  use tglf_interface
  
  implicit none
  integer :: i_ky, i_run, i_colour
  integer :: ios

  character(len=5) :: kystr
  real :: ky, timestep, maxtime, freq_tol

  integer, dimension(:), allocatable :: ky_run
  integer :: run_status

  logical :: file_exists=.false.
  
  ! Map TGYRO parameters to GYRO
  
  gyrotest_flag=0

  ! Reads in GYRO input and sets up GYRO variables
  call gyro_read_input
  call gyro_read_input_extra
  call map_global2interface
  
  gyro_silent_flag_in=1

  ! Maps GYRO variables to QL-GYRO
  call qlgyro_gyro_map

  ! Maps QL-GYRO variable to TGLF
  call qlgyro_tglf_map

  call tglf_dump_local()

  call qlgyro_ky_spectrum

  n_ky = tglf_nky_in

  timestep = gyro_time_step_in
  maxtime = gyro_time_max_in
  freq_tol = abs(maxtime - int(maxtime))
  
  tglf_nmodes_in = n_modes

  ! Ensure GYRO prints out wavefunction and flux data
  gyro_plot_u_flag_in = 1
  gyro_lindiff_method_in = 3

  if (gyro_theta_plot_in .eq. 1) then
     gyro_theta_plot_in = 64
  end if
  
  if (i_tran == 1) then

     ! Eigenvalues
     allocate(tglf_eigenvalue_spectrum_out(2, tglf_nky_in, tglf_nmodes_in))
     
     ! Fields
     allocate(tglf_field_spectrum_out(4, tglf_nky_in, tglf_nmodes_in))
     
     ! Fluxes
     allocate(tglf_flux_spectrum_out(5, n_species, 3, tglf_nky_in, tglf_nmodes_in))
     !allocate(gyro_fluxes_out(loc_n_ion+1, gyro_n_field_in, 4, gyro_radial_grid_in, tglf_nky_in, n_inst))

     ! Array to check if simulation started for each ky
     allocate(ky_run(tglf_nky_in), ky_color(tglf_nky_in))

     ! Grad r geometry factor
     allocate(sat_geo_spectrum(tglf_nky_in))
     
     ! Initial status of each sim
     ky_run = 0

     ! Initial proc no. for ky run
     ky_color = -1
     
     inquire(file=trim(path)//'out.qlgyro.status', exist=file_exists)

     ! Initialise status file with all zeros if doesnt exist
     call init_qlgyro_status
     
  end if

  call MPI_BCAST(ky_run, tglf_nky_in, MPI_INTEGER, 0, QLGYRO_COMM_WORLD, ierr)
  
  tglf_eigenvalue_spectrum_out = 0.0
  tglf_flux_spectrum_out = 0.0
  tglf_field_spectrum_out = 0.0

  call MPI_barrier(QLGYRO_COMM_WORLD, ierr)

  do i_ky=1,tglf_nky_in

     ! Initial assignment of jobs to cores
     if (i_ky .le. n_parallel .and. i_ky .ne. color+1) then
        cycle
     end if

     call get_qlgyro_status(i_ky, run_status)
  
     ! If run is being done then skip to next ky
     if (run_status .ne. 0) then
        cycle
     end if

     run_status = 1
     call read_qlgyro_status(ky_run, ky_color)
     call write_qlgyro_status(ky_run, i_ky, run_status, color)

     if (adjoint .eq. 0) then
        
        print*, ' '
        write(*,21) 'Group ',  color, ' running ky = ', tglf_ky_spectrum_out(i_ky), ky_color(i_ky)
     end if

     ! Ensure linear settings for GYRO
     gyro_toroidal_grid_in = 1
     gyro_nonlinear_flag_in = 0
     gyro_box_multiplier_in = -1
     
     ky = tglf_ky_spectrum_out(i_ky)
     gyro_l_y_in = ky

     ! Scales timestep by 1/ly if larger than 1.0
     if (ky .le. 1.0) then
        gyro_time_step_in = timestep
        gyro_time_max_in  = maxtime
     else
        gyro_time_step_in = timestep * 1.0/ky
        if (maxtime .gt. 0.0) then
           gyro_time_max_in = int(maxtime * 1.0/ky) + freq_tol !- 1.0
        else
           gyro_time_max_in = int(maxtime * 1.0/ky) - freq_tol - 1.0
        end if
     end if

     ! Creates directory for each ky and changes into it for that run
     write(kystr, '(F5.2)') ky
     runpath = 'KY_'//adjustl(kystr)
     
     call system("mkdir -p "//runpath)

     inquire(file=trim(path)//'input.gacode', exist=file_exists)
     if (file_exists .eqv. .true.) then        
        call system("cp input.gacode  "//runpath)
     end if
     
     gyro_path_in = trim(runpath)//'/'
     
     gyro_restart_method=1
     call MPI_barrier(qlgyro_comm, ierr)
     
     do i_run=1,n_runs

        ! Call GYRO
        call gyro_run(gyrotest_flag,gyro_restart_method,transport_method)

        call qlgyro_gyro_cleanup
        
        ! Check if converged
        if (abs(gyro_omega_error_out) .gt. freq_tol) then

           gyro_restart_method=1
           
           if ( i_run == n_runs) then
              run_status = 2              
           end if
        else
           run_status = 3           
           exit
        end if           
     end do
     
     call read_qlgyro_status(ky_run, ky_color)
     call write_qlgyro_status(ky_run, i_ky, run_status, color)

     ! Set run status
     if (adjoint == 0 ) then
        if (run_status .eq. 3) then
           write(*,22) 'GYRO run converged at ky =', ky
        else if (run_status .eq. 2) then
           write(*,22) 'GYRO not converged for ky = ', ky
        end if
     end if

     ! GYRO eigenvalues
     tglf_eigenvalue_spectrum_out(2,i_ky,1) = real(gyro_omega_out)
     tglf_eigenvalue_spectrum_out(1,i_ky,1) = aimag(gyro_omega_out)

     ! GYRO has electrons at last index, TGLF has them in first so load them seperately
     tglf_flux_spectrum_out(1:4, 1, :gyro_n_field_in, i_ky, 1) = gyro_gbflux_out(1:4, n_species, :gyro_n_field_in) * gyro_l_y_in
     
     tglf_flux_spectrum_out(1:4, 2:, :gyro_n_field_in, i_ky, 1) = gyro_gbflux_out(1:4, 1:n_species-1, :gyro_n_field_in) * gyro_l_y_in
     
     !if (adjoint .eq. 0) print*, i_ky, tglf_flux_spectrum_out(2, :, 1, i_ky, 1)
     ! If run hasn't converged, sent eigenvalue/growth to 0
     if (run_status .eq. 2) then
        tglf_eigenvalue_spectrum_out(2,i_ky,1) = 0.0
        tglf_eigenvalue_spectrum_out(1,i_ky,1) = 0.0
        
        ! GYRO has electrons at last index, TGLF has them in first so load them seperately
        tglf_flux_spectrum_out(1:4, :, :gyro_n_field_in, i_ky, 1) = 0.0
     end if

     call gyro_sat_geo_factor(sat_geo_spectrum(i_ky))

     call MPI_BARRIER(qlgyro_comm, ierr)
  end do

  call MPI_BARRIER(QLGYRO_COMM_WORLD, ierr)

  if (i_proc_global == 0) then
     print*, 'Ran all kys'
     open(unit=1,file=trim(runfile),position='append')
     write(1, *) 'INFO: Ran all kys'
     write(1, *) 'B_unit = ', bunit
     close(1)
  endif

  
21 format(A6, I3, A24, F5.2, I3)
22 format(A30, F5.2)
 
end subroutine qlgyro_run_gyro


subroutine gyro_sat_geo_factor(sat_geo_phi_average)
  
  use tglf_interface
  use tglf_max_dimensions
  use tglf_sgrid
  use qlgyro_globals
  use gyro_interface
  
  implicit none

  real :: drmindx, kx0

  integer :: it, tglf_n_theta, jt, t1, t2

  real :: sat_geo_phi_average, qrat_geo_phi, phi_norm
  real :: qrat_1, qrat_2, gyro_theta

  real, dimension(gyro_theta_plot_in*gyro_radial_grid_in) :: gyro_qrat_geo, gyro_int_factor

  real, dimension(:), allocatable :: dR, dZ, dtheta, dl_dtheta, int_factor
  real :: int_factor_1, int_factor_2, d_theta
  
  drmindx = 1.0
  kx0 = 0.0
  
  CALL put_Miller_geometry(tglf_rmin_loc_in, &
       tglf_rmaj_loc_in, &
       tglf_zmaj_loc_in, &
       drmindx, &
       tglf_drmajdx_loc_in, &
       tglf_dzmajdx_loc_in, &
       tglf_kappa_loc_in, &
       tglf_s_kappa_loc_in, &
       tglf_delta_loc_in, &
       tglf_s_delta_loc_in, &
       tglf_zeta_loc_in, &
       tglf_s_zeta_loc_in, &
       tglf_q_loc_in, &
       tglf_q_prime_loc_in, &
       tglf_p_prime_loc_in, &
       kx0)
  
  call tglf_setup_geometry

  ! t_s contains data on theta
  ! Now goes from 2*pi -> 0
  t_s = t_s + 2*pi
  
  tglf_n_theta = ms+1

  allocate(dR(0:ms), dZ(0:ms), dtheta(0:ms), int_factor(0:ms), &
       dl_dtheta(0:ms))

  ! Find integration factor, then interpolating. Not sure if correct order there
  ! Shift arrays forward/back to find dl/dtheta

  dZ = cshift(Z, -1) - cshift(Z, 1)
  dR = cshift(R, -1) - cshift(R, 1)  
  dtheta = cshift(t_s, -1) - cshift(t_s, 1)
  
  dtheta = modulo(dtheta, 2*pi)

  dl_dtheta = sqrt( (dR/dtheta)**2 + (dZ/dtheta)**2)

  ! Integration factor = dl/dtheta * B/Bp
  int_factor = b_geo * dl_dtheta / Bp

  ! Create qrat_geo and int_factor on GYRO thetab grid by interpolating
  do it=1, gyro_theta_plot_in*gyro_radial_grid_in
     gyro_theta = modulo(gyro_thetab_out(it), 2*pi)

     ! Shift it to positive if negative
     if (gyro_theta .lt. 0.0) gyro_theta = 2*pi - gyro_theta
     
     do jt=0, tglf_n_theta-1
        if ( t_s(jt) .le. gyro_theta) then
           ! Index of t_s above and below index of interest
           t1 = jt-1
           t2 = jt
           exit
        end if
     end do
        
     ! TGLF uses qrat_geo**2
     qrat_1 = qrat_geo(t1)**2
     qrat_2 = qrat_geo(t2)**2

     ! Integration measure interpolated
     int_factor_1 = int_factor(t1)
     int_factor_2 = int_factor(t2)
     
     ! Linearly extrapolate and multiply by |Phi|
     gyro_qrat_geo(it) = (qrat_1 + (qrat_2 - qrat_1)*(gyro_theta - t_s(t1)) / &
          (t_s(t2) - t_s(t1))) * abs(gyro_wavefunction_out(1, it))**2
     
     gyro_int_factor(it) = int_factor_1 + (int_factor_2 -int_factor_1) * &
          (gyro_theta - t_s(t1)) / (t_s(t2) - t_s(t1))
     
  end do

  qrat_geo_phi = 0.0
  
  ! Trapezoid rule
  do it=1, gyro_theta_plot_in*gyro_radial_grid_in -1

     ! dtheta
     d_theta = gyro_thetab_out(it+1) - gyro_thetab_out(it)

     ! Wavefunction average using composite trapezoidal rule
     qrat_geo_phi = qrat_geo_phi +  d_theta * gyro_int_factor(it) * &
          (gyro_qrat_geo(it+1) + gyro_qrat_geo(it))/2
     
     phi_norm = phi_norm + d_theta * gyro_int_factor(it) * &
          (abs(gyro_wavefunction_out(1, it+1))**2 + abs(gyro_wavefunction_out(1, it))**2)/2

  end do

  ! SAT_geo is 1/qrat_geo
  sat_geo_phi_average = 1.0 / (qrat_geo_phi / phi_norm)
  
end subroutine gyro_sat_geo_factor


