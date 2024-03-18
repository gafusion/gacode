  !--------------------------------------------------------
  ! qlgyro_comm_sync.f90
  !
  ! PURPOSE
  ! Synchronises (broadcast) outputs on all color groups
  ! -------------------------------------------------------


subroutine qlgyro_comm_sync

  use mpi
  use qlgyro_globals
  use tglf_interface
  
  implicit none

  integer :: i_kypx0, core, i_theta
  integer, dimension(n_kypx0) :: ky_run
  
  ! Sync all proc
  call MPI_BARRIER(QLGYRO_COMM_WORLD, ierr)
  
  ! Get core for each 
  call read_qlgyro_status(ky_run, ky_color)
  
  do i_kypx0 = 1, n_kypx0

     ! Core to bcast data
     core = ky_color(i_kypx0) * procs
     
     call MPI_BCAST(qlgyro_eigenvalue_spectrum_out(1:2, i_ky(i_kypx0), i_px0(i_kypx0)) , &
          2, &
          MPI_DOUBLE, &
          core, &
          QLGYRO_COMM_WORLD, &
          ierr)

     call MPI_BCAST(qlgyro_flux_spectrum_out(1:5, 1:n_species, 1:3,  i_ky(i_kypx0), i_px0(i_kypx0)), &
          (5* n_species* 3),  &
          MPI_DOUBLE, &
          core, &
          QLGYRO_COMM_WORLD, &
          ierr)

     call MPI_BCAST(sat_geo_spectrum(i_ky(i_kypx0)), &
          1,  &
          MPI_DOUBLE, &
          core, &
          QLGYRO_COMM_WORLD, &
          ierr)

     call MPI_BCAST(qlgyro_field_spectrum_out(1:3, 1:n_thetab, i_ky(i_kypx0), i_px0(i_kypx0)), &
          (6*n_thetab),  &
          MPI_DOUBLE, &
          core, &
          QLGYRO_COMM_WORLD, &
          ierr)

     call MPI_BCAST(qlgyro_k_perp(1:n_thetab, i_ky(i_kypx0), i_px0(i_kypx0)), &
          n_thetab,  &
          MPI_DOUBLE, &
          core, &
          QLGYRO_COMM_WORLD, &
          ierr)
          
     call MPI_BCAST(qlgyro_jacobian(1:n_thetab, i_ky(i_kypx0), i_px0(i_kypx0)), &
          n_thetab,  &
          MPI_DOUBLE, &
          core, &
          QLGYRO_COMM_WORLD, &
          ierr)

  end do

end subroutine qlgyro_comm_sync
