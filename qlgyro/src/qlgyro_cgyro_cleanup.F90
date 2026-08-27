!-----------------------------------------------------------
! qlgyro_cgyro_cleanup.f90
!
! PURPOSE:
!  Assign output to interface, Deallocate and clean up.
!-----------------------------------------------------------

subroutine qlgyro_cgyro_cleanup

  use mpi
  use parallel_lib
  use timer_lib
  use cgyro_globals
  use cgyro_field_mod, only : field
  use cgyro_flux_mod, only : gflux
  use qlgyro_cgyro_interface
  use qlgyro_globals

  implicit none

  integer :: l, i_field, ir, it, pos, start, end
  complex :: a_norm
  complex :: ftemp(cgyro_n_theta_in, cgyro_n_radial_in)

  ! No need for cleanup in test mode
  if (test_flag == 1) return

  ! Manage exit message
  cgyro_error_status_out = error_status
  cgyro_signal_out = signal
  
  if (signal == 1) then
     cgyro_error_message_out = 'Linear converged'
  else
     cgyro_error_message_out = 'Linear terminated at max time'
  endif

  ! Write eigenvalues to interface
  cgyro_omega_out = freq
  ! EAB
  !cgyro_omega_error_out = freq_err
  cgyro_omega_error_out = 0.0

  bgs2 = b_gs2

  ! Convert to right units now
  cgyro_k_perp_out = k_perp * rho / bgs2
  
  do ir=1, cgyro_n_radial_in
     start = (ir - 1)* cgyro_n_theta_in + 1
     end =  ir * cgyro_n_theta_in
     cgyro_jacobian_out(start:end) = g_theta / bmag
  end do

  ! Fluxes - reorder output for QLGYRO
  ! Re-normalise 
  do i_tor=1,cgyro_n_toroidal_in
     do l=1,3
        cgyro_gbflux_out(l, :, :, i_tor) = real(gflux(0, :, l, :, i_tor))
     end do
  end do

  ! Wavefunction
  do i_tor=1,cgyro_n_toroidal_in
     do i_field=1,cgyro_n_field_in
        do ir=1,cgyro_n_radial_in
           do it=1,cgyro_n_theta_in
              ftemp(it,ir) = field(i_field,ic_c(ir,it),i_tor)
           enddo
        enddo

        if (i_field == 1) then
           it = maxloc(abs(ftemp(:,cgyro_n_radial_in/2+1)),dim=1)
           a_norm = ftemp(it,cgyro_n_radial_in/2+1)
        endif
        call extended_ang(ftemp)
     
        ftemp = ftemp / a_norm
        
        cgyro_wavefunction_out(i_tor, i_field, :) = reshape(ftemp, (/cgyro_n_radial_in * cgyro_n_theta_in/))
  end do
enddo

  if (auto_box_size .eq. 0) then
     ! Get ballooning theta from out.cgyro.grids
     open(unit=io, file=trim(runpath)//'out.cgyro.grids', status='old')
     
     ! Skip to thetab line in file
     pos = (11 + cgyro_n_radial_in + cgyro_n_theta_in + cgyro_n_energy_in + cgyro_n_xi_in)
     do l=1,pos
        read(io, *)
     end do
    
     read(io, '(1pe13.6)') cgyro_thetab_out
     close(io)
  else
     cgyro_thetab_out = 0.0
  end if
  ! Free up commuicators
  call MPI_COMM_FREE(NEW_COMM_1,i_err)
  call MPI_COMM_FREE(NEW_COMM_2,i_err)

  ! If in TGYRO remove tag
  if (tag_removal .eq. 1) then
     if (adjoint .eq. 0) then
        open(unit=io, iostat=i_err, file=trim(runpath)//runfile_restart_tag, status='old')
        if (i_err == 0) close(io, status='delete')
        if (restart_mode .eq. 1) then
           open(unit=io, iostat=i_err, file=trim(runpath)//runfile_restart, status='old')
           if (i_err == 0) close(io, status='delete')
        end if
     end if
  end if
  
  ! Reset timer info
  timer_cpu_maxindx = 0
  timer_cpu_tag = ''
  timer_cpu=0.0
  timer_cpu_in=0.0

  ! Deallocate arrays
  !call cgyro_deallocate_arrays
  
end subroutine qlgyro_cgyro_cleanup


subroutine qlgyro_cgyro_deallocate_arrays

  !use cgyro_globals
  !use parallel_lib
  implicit none

  call cgyro_cleanup
  
end subroutine qlgyro_cgyro_deallocate_arrays

