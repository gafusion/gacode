!-----------------------------------------------------------------
! qlgyro_run.f90
!
! PURPOSE:
!  Manage call to QLGYRO simulation for both standalone and
!  TGYRO usage.
!-----------------------------------------------------------------

subroutine qlgyro_run(lpath_in, qlgyro_comm_in, i_tran_in)

  use mpi
  use qlgyro_globals
  use qlgyro_cgyro_interface

  !-----------------------------------------------------------------
  implicit none
  !
  !-----------------------------------------------------------------

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: lpath_in
  integer, intent(in) :: qlgyro_comm_in
  integer, intent(in) :: i_tran_in
  character(len=13) :: format_string

  !-----------------------------------------------------------------
  
  !-----------------------------------------------------------------
  ! Query MPI for dimensions
  !
  QLGYRO_COMM_WORLD = qlgyro_comm_in
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(QLGYRO_COMM_WORLD,i_proc_global,ierr)
  call MPI_COMM_SIZE(QLGYRO_COMM_WORLD,n_proc_global,ierr)
  n_inst = n_proc_global
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Transport solver operation
  transport_method = 1
  tag_removal = 1
  i_tran = i_tran_in
  !

  path = lpath_in
  cgyro_path_in = path

  if (transport_method .eq. 1) then
     if (i_tran .lt. 10) then
        format_string = "(A10, I1, A1)"
     else if (i_tran .lt. 100) then
        format_string = "(A10, I2, A1)"
     else if (i_tran .lt. 1000) then
        format_string = "(A10, I3, A1)"
     else
        write(*, *) "Too many iterations... over 1000..."
        stop
     end if

     write(iter_path, format_string) "ITERATION_", i_tran, "/"

     ! Adjust path directory if running TGYRO in parallel
     if (cgyro_iteration_method_in .eq. 5) then
        write(worker_str, *) cgyro_worker_index_in
        iter_path = trim(iter_path)//trim(adjustl(worker_str))//"/"
     end if
     if (i_proc_global .eq. 0) then
        call system("mkdir -p "//trim(path)//trim(iter_path))
     end if
  else
     iter_path = "./"
  end if

  call qlgyro_allocate_globals

  call qlgyro_read_input

  call qlgyro_comm_setup

  CGYRO_COMM_WORLD = qlgyro_comm
  call cgyro_init(path, CGYRO_COMM_WORLD)

  if (cgyro_n_toroidal_in .gt. 1) then
     call qlgyro_run_cgyro_fluxtube
  else 
     call qlgyro_run_cgyro_balloon
  end if

  ! Synchronise all cores
  call qlgyro_comm_sync

  ! Apply saturation rule
  if (sat_rule .gt. 0) then
     call qlgyro_sat1
  else if (sat_rule .eq. -1) then
     call qlgyro_sat_mg
  else
     open(unit=1,file=trim(runfile),position='append')
     write(1,*) '-----------------------------------------'
     write(1,*) 'Sat rule ', sat_rule, 'not a valid option'
     write(1,*) '-----------------------------------------'
     close(1)
     stop
  end if

  ! Sum flux over ky spectrum
  call qlgyro_sum_fluxes

  if (i_proc_global .eq. 0) then
     call write_qlgyro_ky_spectrum
     call write_qlgyro_QL_weight_spectrum
     call write_qlgyro_eigenvalue_spectrum
     call write_qlgyro_flux_spectrum
     call write_qlgyro_field_spectrum
     call write_qlgyro_gbflux
     call write_qlgyro_units
     call write_qlgyro_sat_geo_spectrum
     call write_qlgyro_kxrms_spectrum

     open(unit=1,file=trim(runfile),position='append')
     write(1,*) 'Completed QLGYRO run'
     close(1)

  end if

  ! Need to deallocate array which get re-allocated in next call
  call tglf_deallocate

end subroutine qlgyro_run
