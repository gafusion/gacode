!-----------------------------------------------------------
! gyro_cleanup.f90
!
! PURPOSE:
!  Assign output to interface, Deallocate and clean up.
!-----------------------------------------------------------

subroutine qlgyro_gyro_cleanup

  use mpi
  use qlgyro_globals
  use gyro_globals
  use gyro_interface

  implicit none

  integer :: l, i_field, n_ball, io, pos, n_flux
  complex :: omega_data(2)
  complex :: ftemp(gyro_theta_plot_in*gyro_radial_grid_in)
  real(4) :: flux_temp(n_kinetic, gyro_n_field_in, p_moment, n_x), flux_loc_temp(n_kinetic, gyro_n_field_in, p_moment)
  character(len=24), dimension(3) :: field_files

  io = 21
  
  ! Write eigenvalues to interface from final line in out.gyro.freq file
  open(unit=io, file=trim(runpath)//'/out.gyro.freq', position="append")
  backspace(unit=io)
  read(io, 20) omega_data
  close(io)
  
20 format(4(es11.4,1x))

  gyro_omega_out = omega_data(1)
  gyro_omega_error_out = omega_data(2)

  ! Fluxes - reorder output for QLGYRO
  if(.not.allocated(gyro_gbflux_out)) allocate(gyro_gbflux_out(4, n_species, gyro_n_field_in))

  n_flux = n_kinetic * gyro_n_field_in * 4 * n_x

  ! Read out.gyro.gbflux_i file and skip to the final time step
  pos = (data_step * n_flux* BYTE) +1
  
  open(unit=io,file=trim(runpath)//'/bin.gyro.gbflux_i',status='old',access='stream', form='unformatted')
  read(io, pos=pos) flux_temp
  close(io)

  ! Calculate average in radial positions of interest
  flux_loc_temp = 0.0
  do l=1+n_explicit_damp,n_x-n_explicit_damp
     flux_loc_temp(:,:,:) = flux_loc_temp(:,:,:)+flux_temp(:,:,:,l)
  enddo ! l
  
  flux_loc_temp(:,:,:) = flux_loc_temp(:,:,:)/(n_x-2*n_explicit_damp)

  ! Re order output for QL-GYRO
  do l=1,p_moment
     gyro_gbflux_out(l, :, :) = real(flux_loc_temp(:, :, l))
  end do

  
  ! Wavefunction read from file
  field_files = (/'/out.gyro.balloon_phi  ', '/out.gyro.balloon_a    ', '/out.gyro.balloon_aperp'/)

  n_ball =  gyro_theta_plot_in*gyro_radial_grid_in
  
  if(.not.allocated(gyro_wavefunction_out)) allocate(gyro_wavefunction_out(gyro_n_field_in, n_ball))
  if(.not.allocated(gyro_thetab_out)) allocate(gyro_thetab_out(n_ball))
  
  do i_field=1,gyro_n_field_in

     ! Start at end and work back to final timestep
     open(unit=io,file=trim(runpath)//trim(field_files(i_field)), status='old', position='append')

     do l=1, n_ball
        backspace(io)
     end do
     
     do l=1, n_ball
        read(io, fmtstr2) ftemp(l)
     end do
     close(io)
     
     gyro_wavefunction_out(i_field, :) = ftemp

  end do

  ! Theta ballooning grid
  do l=1, n_ball
     gyro_thetab_out(l) = -(1.0+gyro_radial_grid_in)+(2.0*gyro_radial_grid_in*(l-1.0)/n_ball)
  end do
  
  gyro_thetab_out = gyro_thetab_out * pi

end subroutine qlgyro_gyro_cleanup

