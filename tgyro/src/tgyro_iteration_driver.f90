!-----------------------------------------------------------
! tgyro_iteration_driver.f90
!
! PURPOSE:
!  Main driver for local solver.
!
!  Quantities to evolve:
!   Ti: dlntidr
!   Te: dlntedr
!   ne: dlnnedr
!   Er: f_rot = [a/c_s(0)]*r_maj(0)*gamma_p/r_maj
!----------------------------------------------------------

subroutine tgyro_iteration_driver

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables

  implicit none
  integer :: i_ion

  n_r   = n_inst+1
  p_max = n_evolve*(n_r-1)

  flux_counter = 0

  allocate(ipiv(p_max))
  allocate(jf(p_max,p_max))
  allocate(jg(p_max,p_max))
  allocate(jfg(p_max,p_max))
  allocate(x_vec(p_max)) 
  allocate(f_vec(p_max))
  allocate(g_vec(p_max))
  allocate(x_vec0(p_max)) 
  allocate(f_vec0(p_max))
  allocate(g_vec0(p_max))
  allocate(b(p_max))

  call tgyro_allocate_globals
 
  !---------------------------------------
  ! Error monitoring variables
  !
  error_flag = 0
  error_msg = 'INFO: clean exit from TGYRO'
  !
  b_flag(:) = ' ' 
  !---------------------------------------

  ! Generate ALL radial profiles.
  call tgyro_init_profiles

  ! Step-length for Jacobian
  dx = loc_dx/r_min

  !----------------------------------------------
  ! Choose flux_method based on path information.
  !
  if (tgyro_noturb_flag == 1) then

     flux_method = 0

  else if (lcode == 'ifs') then

     ! IFS-PPPL
     flux_method = 1

  else if (lcode == 'tglf') then

     ! TGLF
     flux_method = 2

  else if (lcode == 'etg') then
     
     ! ETG critical gradient
     flux_method = 5

  else if (lcode == 'mmm') then
     
     ! Multi-Mode model gradient
     flux_method = 6

  else if (lcode == 'qlgyro') then

     ! QLGYRO
     flux_method = 8

  endif

  !---------------------------------------------
  allocate(flux_method_vec(n_inst))
  call MPI_ALLGATHER(flux_method,1,MPI_INTEGER,flux_method_vec,1,MPI_INTEGER,gyro_adj,ierr)
  call tgyro_write_input

  if (loc_restart_flag == 0) then
     ! Initialize relaxation parameters to starting value.
     relax(:) = 1.0 ; res(:) = 0.0
     ! Create, but do not write to, datafiles.
     call tgyro_write_data(0)
  endif
  call tgyro_mpi_info('INFO: (tgyro_iteration_driver) Survived initialization, starting iterations')

  correct_flag = 0

  mask = 0
  ! Mapping function from radius/field to p
  p = 0
  do i=2,n_r
     ip = 0
     if (loc_ti_feedback_flag == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = dlntidr(1,i)
        mask(p,1) = 1 
     endif
     if (loc_te_feedback_flag == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = dlntedr(i)
        mask(p,2) = 1
     endif
     if (loc_er_feedback_flag == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = f_rot(i)
        mask(p,3) = 1
     endif
     if (evo_e(0) == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = dlnnedr(i)
        mask(p,4) = 1
     endif
     do i_ion=1,loc_n_ion
        if (evo_e(i_ion) >= 1) then
           p  = p+1
           ip = ip+1
           pmap(i,ip) = p
           x_vec(p) = dlnnidr(i_ion,i)
           mask(p,4+i_ion) = 1
        endif
     enddo
  enddo

  ! Make some resets if we are in test mode
  if (gyrotest_flag == 1) then
     tgyro_iteration_method = 1
     tgyro_relax_iterations = 0
  endif
 
  select case (tgyro_iteration_method) 

  case (1) 

     call tgyro_iteration_zero
     call tgyro_iteration_standard

  case (2) 

     call tgyro_iteration_zero
     call tgyro_iteration_diagonal

  case (4) 

     call tgyro_iteration_serial

  case (5) 

     call tgyro_iteration_parallel

  case (6) 

     call tgyro_iteration_zero
     call tgyro_iteration_simplerelax

  end select
    
end subroutine tgyro_iteration_driver
