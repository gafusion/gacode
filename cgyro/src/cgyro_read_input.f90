subroutine cgyro_read_input

  use mpi
  use cgyro_globals

  implicit none

  integer :: is
  real :: x
  character (len=1) :: cdummy

  if (i_proc == 0) open(unit=1,file=trim(path)//'input.cgyro.gen',status='old')

  call cgyro_readbc_int(n_energy)
  call cgyro_readbc_int(n_xi)
  call cgyro_readbc_int(n_theta)
  call cgyro_readbc_int(n_radial)
  call cgyro_readbc_int(n_toroidal)
  call cgyro_readbc_int(n_field)
  call cgyro_readbc_real(e_max)
  call cgyro_readbc_int(e_method)
  call cgyro_readbc_real(delta_t)
  call cgyro_readbc_real(max_time)
  call cgyro_readbc_int(print_step)
  call cgyro_readbc_int(restart_step)
  call cgyro_readbc_real(freq_tol)
  call cgyro_readbc_int(restart_mode)
  call cgyro_readbc_real(up_radial)
  call cgyro_readbc_real(up_theta)
  call cgyro_readbc_real(up_alpha)
  call cgyro_readbc_int(nup_radial)
  call cgyro_readbc_int(nup_theta)
  call cgyro_readbc_int(nup_alpha)
  call cgyro_readbc_int(n_wave)
  call cgyro_readbc_int(constant_stream_flag)
  call cgyro_readbc_real(ky)
  call cgyro_readbc_int(box_size)
  call cgyro_readbc_real(ipccw)
  call cgyro_readbc_real(btccw)
  call cgyro_readbc_int(silent_flag)
  call cgyro_readbc_int(profile_model)
  call cgyro_readbc_int(equilibrium_model)
  call cgyro_readbc_int(collision_model)
  call cgyro_readbc_int(collision_mom_restore)
  call cgyro_readbc_int(collision_ene_restore)
  call cgyro_readbc_int(collision_ene_diffusion)
  call cgyro_readbc_int(collision_kperp)
  call cgyro_readbc_int(collision_field_model)
  call cgyro_readbc_int(collision_ion_model)
  call cgyro_readbc_real(collision_ele_scale)
  call cgyro_readbc_real(z_eff)
  call cgyro_readbc_int(z_eff_method)
  call cgyro_readbc_int(zf_test_flag)
  call cgyro_readbc_int(nonlinear_flag)
  call cgyro_readbc_int(nonlinear_method)
  call cgyro_readbc_real(te_ade)
  call cgyro_readbc_real(ne_ade)
  call cgyro_readbc_real(masse_ade)
  call cgyro_readbc_real(dlntdre_ade)
  call cgyro_readbc_real(dlnndre_ade)
  call cgyro_readbc_real(lambda_star)
  call cgyro_readbc_int(test_flag)
  call cgyro_readbc_int(h_print_flag)
  call cgyro_readbc_int(moment_print_flag)
  call cgyro_readbc_int(kxkyflux_print_flag)
  call cgyro_readbc_real(amp0)
  call cgyro_readbc_real(amp)
  call cgyro_readbc_real(gamma_e)
  call cgyro_readbc_real(gamma_p)
  call cgyro_readbc_real(mach)
  call cgyro_readbc_int(rotation_model)
  call cgyro_readbc_real(error_tol)
  call cgyro_readbc_real(adapt_tol)
  call cgyro_readbc_int(mpi_rank_order)
  call cgyro_readbc_int(hiprec_flag)
  call cgyro_readbc_int(udsymmetry_flag)
  call cgyro_readbc_int(shear_method)
  call cgyro_readbc_int(n_global)
  call cgyro_readbc_real(nu_global)
  call cgyro_readbc_int(psym_flag)
  call cgyro_readbc_int(profile_shear_flag)
  call cgyro_readbc_int(theta_plot)
  call cgyro_readbc_int(mpiio_small_stripe_factor)
  call cgyro_readbc_int(mpiio_stripe_factor)
  call cgyro_readbc_int(use_bin)

  call cgyro_readbc_real(rmin)
  call cgyro_readbc_real(rmaj)
  call cgyro_readbc_real(q)
  call cgyro_readbc_real(s)
  call cgyro_readbc_real(shift)    
  call cgyro_readbc_real(kappa)   
  call cgyro_readbc_real(s_kappa)  
  call cgyro_readbc_real(delta)       
  call cgyro_readbc_real(s_delta)
  call cgyro_readbc_real(zeta)      
  call cgyro_readbc_real(s_zeta)
  call cgyro_readbc_real(zmag)       
  call cgyro_readbc_real(dzmag)
  call cgyro_readbc_real(betae_unit)
  call cgyro_readbc_int(subroutine_flag)
  call cgyro_readbc_int(n_species)

  call cgyro_readbc_real(nu_ee)

  do is=1,6
     call cgyro_readbc_real(x)   ; z(is) = x
     call cgyro_readbc_real(x)   ; mass(is) = x
     call cgyro_readbc_real(x)   ; dens(is) = x
     call cgyro_readbc_real(x)   ; temp(is) = x
     call cgyro_readbc_real(x)   ; dlnndr(is) = x
     call cgyro_readbc_real(x)   ; dlntdr(is) = x
     call cgyro_readbc_real(x)   ; sdlnndr(is) = x
     call cgyro_readbc_real(x)   ; sdlntdr(is) = x
  enddo

  call cgyro_readbc_int(quasineutral_flag)
  call cgyro_readbc_real(lambda_star_scale)
  call cgyro_readbc_real(gamma_e_scale)
  call cgyro_readbc_real(gamma_p_scale)
  call cgyro_readbc_real(mach_scale)
  call cgyro_readbc_real(q_scale)
  call cgyro_readbc_real(s_scale)
  call cgyro_readbc_real(shift_scale)
  call cgyro_readbc_real(kappa_scale)
  call cgyro_readbc_real(delta_scale)
  call cgyro_readbc_real(zeta_scale)
  call cgyro_readbc_real(s_kappa_scale)
  call cgyro_readbc_real(s_delta_scale)
  call cgyro_readbc_real(s_zeta_scale)
  call cgyro_readbc_real(beta_star_scale)
  call cgyro_readbc_real(betae_unit_scale)
  call cgyro_readbc_real(nu_ee_scale)
  do is=1,6
     call cgyro_readbc_real(x)   ; dlnndr_scale(is) = x
     call cgyro_readbc_real(x)   ; dlntdr_scale(is) = x
  enddo

  if (i_proc == 0) close(1)

  ! GEO fourier coefficients
  geo_ny_in = 0
  geo_yin_in(:,:) = 0.0
  if (subroutine_flag == 0 .and. equilibrium_model == 3 & 
       .and. profile_model == 1) then
     if (i_proc == 0) then
        open(unit=1,file=trim(path)//'input.geo',status='old')
        ! header skip
        do
           read(1,'(a)') cdummy
           if (cdummy /= '#') exit
        enddo
        backspace 1
        ! n_fourier
        read(1,*) geo_ny_in
        ! fourier coefficients
        read(1,*) geo_yin_in(:,0:geo_ny_in)
        close(1)
     endif
     call MPI_BCAST(geo_ny_in,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)
     call MPI_BCAST(geo_yin_in,size(geo_yin_in),MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)
  endif

end subroutine cgyro_read_input

!------------------------------------------------------------
! Service routines: 
!
! (1) read and broadcast an integer:
!
subroutine cgyro_readbc_int(p)

  use mpi
  use cgyro_globals

  implicit none
  integer, intent(inout) :: p

  if (i_proc == 0) read(1,*) p

  call MPI_BCAST(p,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_readbc_int
!
! (2) read and broadcast a real:
!
subroutine cgyro_readbc_real(x)
  
  use mpi
  use cgyro_globals

  implicit none
  real, intent(inout) :: x

  if (i_proc == 0) read(1,*) x

  call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_readbc_real
!------------------------------------------------------------
