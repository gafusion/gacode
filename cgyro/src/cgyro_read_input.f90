subroutine cgyro_read_input

  use mpi
  use cgyro_globals

  implicit none

  integer :: is
  character (len=1) :: cdummy

  if (i_proc == 0) open(unit=1,file=trim(path)//'input.cgyro.gen',status='old')

  call cgyro_readbc_int(n_energy)
  call cgyro_readbc_int(n_xi)
  call cgyro_readbc_int(n_theta)
  call cgyro_readbc_int(n_radial)
  call cgyro_readbc_int(n_toroidal)
  call cgyro_readbc_int(n_field)
  call cgyro_readbc_real(e_max)
  call cgyro_readbc_real(alpha_poly)
  call cgyro_readbc_int(e_method)
  call cgyro_readbc_int(e_fix)
  call cgyro_readbc_int(delta_t_method)
  call cgyro_readbc_real(delta_t)
  call cgyro_readbc_real(error_tol)
  call cgyro_readbc_real(max_time)
  call cgyro_readbc_int(print_step)
  call cgyro_readbc_int(restart_step)
  call cgyro_readbc_real(freq_tol)
  call cgyro_readbc_real(up_radial)
  call cgyro_readbc_real(up_theta)
  call cgyro_readbc_real(up_alpha)
  call cgyro_readbc_int(nup_radial)
  call cgyro_readbc_int(nup_theta)
  call cgyro_readbc_int(nup_alpha)
  call cgyro_readbc_int(n_wave)
  call cgyro_readbc_int(constant_stream_flag)
  call cgyro_readbc_int(explicit_trap_flag)
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
  call cgyro_readbc_int(collision_precision_mode)
  call cgyro_readbc_real(z_eff)
  call cgyro_readbc_int(z_eff_method)
  call cgyro_readbc_int(zf_test_mode)
  call cgyro_readbc_int(nonlinear_flag)
  call cgyro_readbc_int(nonlinear_method)
  call cgyro_readbc_int(ae_flag)
  call cgyro_readbc_real(temp_ae)
  call cgyro_readbc_real(dens_ae)
  call cgyro_readbc_real(mass_ae)
  call cgyro_readbc_real(dlntdr_ae)
  call cgyro_readbc_real(dlnndr_ae)
  call cgyro_readbc_real(lambda_star)
  call cgyro_readbc_int(h_print_flag)
  call cgyro_readbc_int(moment_print_flag)
  call cgyro_readbc_int(gflux_print_flag)
  call cgyro_readbc_int(field_print_flag)
  call cgyro_readbc_real(amp0)
  call cgyro_readbc_real(amp)
  call cgyro_readbc_real(gamma_e)
  call cgyro_readbc_real(gamma_p)
  call cgyro_readbc_real(mach)
  call cgyro_readbc_int(rotation_model)
  call cgyro_readbc_int(mpi_rank_order)
  call cgyro_readbc_int(velocity_order)
  call cgyro_readbc_int(hiprec_flag)
  call cgyro_readbc_int(udsymmetry_flag)
  call cgyro_readbc_int(shear_method)
  call cgyro_readbc_int(n_global)
  call cgyro_readbc_real(nu_global)
  call cgyro_readbc_int(profile_shear_flag)
  call cgyro_readbc_int(theta_plot)
  call cgyro_readbc_int(gpu_bigmem_flag)
  call cgyro_readbc_int(upwind_single_flag)
  call cgyro_readbc_real(px0)
  call cgyro_readbc_int(stream_term)
  call cgyro_readbc_real(stream_factor)
  call cgyro_readbc_int(exch_flag)

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
  do is=3,n_shape
     call cgyro_readbc_real(shape_sin(is))
     call cgyro_readbc_real(shape_s_sin(is))
  enddo
  do is=0,n_shape
     call cgyro_readbc_real(shape_cos(is))
     call cgyro_readbc_real(shape_s_cos(is))
  enddo
  call cgyro_readbc_real(betae_unit)
  call cgyro_readbc_int(n_species)
  call cgyro_readbc_real(nu_ee)

  ! vectors
  is = size(z)
  call cgyro_readbc_realv(z,is)   
  call cgyro_readbc_realv(mass,is) 
  call cgyro_readbc_realv(dens,is)   
  call cgyro_readbc_realv(temp,is)
  call cgyro_readbc_realv(dlnndr,is)
  call cgyro_readbc_realv(dlntdr,is)
  call cgyro_readbc_realv(sdlnndr,is)
  call cgyro_readbc_realv(sdlntdr,is)
  call cgyro_readbc_realv(dlnndr_scale,is)
  call cgyro_readbc_realv(dlntdr_scale,is)

  call cgyro_readbc_int(quasineutral_flag)
  call cgyro_readbc_real(lambda_star_scale)
  call cgyro_readbc_real(gamma_e_scale)
  call cgyro_readbc_real(gamma_p_scale)
  call cgyro_readbc_real(mach_scale)
  call cgyro_readbc_real(beta_star_scale)
  call cgyro_readbc_real(betae_unit_scale)
  call cgyro_readbc_real(nu_ee_scale)

  if (i_proc == 0) close(1)

end subroutine cgyro_read_input

!------------------------------------------------------------
! Service routines: 
!
! (1) read and broadcast an integer:
!
subroutine cgyro_readbc_int(p)

  use mpi
  use cgyro_globals, only : i_proc,i_err,CGYRO_COMM_WORLD

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
  use cgyro_globals, only : i_proc,i_err,CGYRO_COMM_WORLD

  implicit none
  real, intent(inout) :: x

  if (i_proc == 0) read(1,*) x

  call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_readbc_real
!
! (2) read and broadcast a vector of reals
!
subroutine cgyro_readbc_realv(x,n)

  use mpi
  use cgyro_globals, only : i_proc,i_err,CGYRO_COMM_WORLD

  implicit none
  integer, intent(in) :: n
  real, intent(inout) :: x(n)
  integer :: i
  
  if (i_proc == 0) then
     do i=1,n
        read(1,*) x(i)
     enddo
  endif

  call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_readbc_realv
!------------------------------------------------------------
