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

  use tgyro_globals
  use tgyro_iteration_variables
  use tgyro_ped
  use mpi

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

  else if (lcode == "ifs") then

     ! IFS-PPPL
     flux_method = 1

  else if (lcode == "tglf") then

     ! TGLF
     flux_method = 2

  else if (lcode == "glf23") then

     ! GLF23
     flux_method = 3

  else

     ! GYRO
     flux_method = 4

     ! Reset step-length for Jacobian
     dx = loc_dx_gyro/r_min

  endif
  !---------------------------------------------
  allocate(flux_method_vec(n_inst))
  call MPI_ALLGATHER(flux_method,1,MPI_INTEGER,flux_method_vec,1,MPI_INTEGER,gyro_adj,ierr)

  call tgyro_write_input

  ! NOTE: See gyro/src/gyro_globals.f90 for definition of transport_method

  ! Standard transport calculation
  if (tgyro_gyro_restart_flag == 0) then
     transport_method = 2
  else
     transport_method = 3
  endif

  if (loc_restart_flag == 0) then
     ! Create, but do not write to, datafiles.
     call tgyro_write_data(0)
     ! Initialize relaxation parameters to starting value.
     relax(:) = 1.0
  endif
  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,'(t2,a)') 'INFO: (TGYRO) Survived initialization and starting iterations'
     close(1)
  endif

  correct_flag = 0

  ! Mapping function from radius/field to p
  p = 0
  do i=2,n_r
     ip = 0
     if (loc_ti_feedback_flag == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = dlntidr(1,i)
     endif
     if (loc_te_feedback_flag == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = dlntedr(i)
     endif
     if (loc_er_feedback_flag == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = f_rot(i)
     endif
     if (evo_e(0) == 1) then
        p  = p+1
        ip = ip+1
        pmap(i,ip) = p
        x_vec(p) = dlnnedr(i)
     endif
     do i_ion=1,loc_n_ion
        if (evo_e(i_ion) == 1) then
           p  = p+1
           ip = ip+1
           pmap(i,ip) = p
           x_vec(p) = dlnnidr(i_ion,i)
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

     call tgyro_iteration_standard

  case (4) 

     call tgyro_iteration_serial

  case (5) 

     call tgyro_iteration_parallel

  case (6) 

     call tgyro_iteration_simplerelax

  end select

end subroutine tgyro_iteration_driver
