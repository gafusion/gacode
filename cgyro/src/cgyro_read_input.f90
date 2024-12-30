subroutine cgyro_read_input

  use mpi
  use cgyro_globals

  implicit none

  integer :: is

  if (i_proc == 0) open(unit=1,file=trim(path)//'input.cgyro.gen',status='old')

  call cgyro_readbc_int(n_energy,'N_ENERGY')
  call cgyro_readbc_int(n_xi,'N_XI')
  call cgyro_readbc_int(n_theta,'N_THETA')
  call cgyro_readbc_int(n_radial,'N_RADIAL')
  call cgyro_readbc_int(n_toroidal,'N_TOROIDAL')
  call cgyro_readbc_int(n_field,'N_FIELD')
  call cgyro_readbc_real(e_max)
  call cgyro_readbc_real(alpha_poly)
  call cgyro_readbc_int(e_fix,'E_FIX')
  call cgyro_readbc_int(delta_t_method,'DELTA_T_METHOD')
  call cgyro_readbc_real(delta_t)
  call cgyro_readbc_real(error_tol)
  call cgyro_readbc_real(max_time)
  call cgyro_readbc_int(print_step,'PRINT_STEP')
  call cgyro_readbc_int(restart_step,'RESTART_STEP')
  call cgyro_readbc_int(restart_preservation_mode,'RESTART_PRESERVATION_MODE')
  call cgyro_readbc_int(mpiio_stripe_factor,'MPIIO_STRIPE_FACTOR')
  call cgyro_readbc_int(mpiio_small_stripe_factor,'MPIIO_SMALL_STRIPE_FACTOR')
  call cgyro_readbc_real(freq_tol)
  call cgyro_readbc_real(up_radial)
  call cgyro_readbc_real(up_theta)
  call cgyro_readbc_real(up_alpha)
  call cgyro_readbc_int(nup_radial,'NUP_RADIAL')
  call cgyro_readbc_int(nup_theta,'NUP_THETA')
  call cgyro_readbc_int(nup_alpha,'NUP_ALPHA')
  call cgyro_readbc_int(n_wave,'N_WAVE')
  call cgyro_readbc_int(constant_stream_flag,'CONSTANT_STREAM_FLAG')
  call cgyro_readbc_int(explicit_trap_flag,'EXPLICIT_TRAP_FLAG')
  call cgyro_readbc_real(ky)
  call cgyro_readbc_int(box_size,'BOX_SIZE')
  call cgyro_readbc_real(ipccw)
  call cgyro_readbc_real(btccw)
  call cgyro_readbc_int(silent_flag,'SILENT_FLAG')
  call cgyro_readbc_int(profile_model,'PROFILE_MODEL')
  call cgyro_readbc_int(equilibrium_model,'EQUILIBRIUM_MODEL')
  call cgyro_readbc_int(collision_model,'COLLISION_MODEL')
  call cgyro_readbc_int(collision_mom_restore,'COLLISION_MOM_RESTORE')
  call cgyro_readbc_int(collision_ene_restore,'COLLISION_ENE_RESTORE')
  call cgyro_readbc_int(collision_ene_diffusion,'COLLISION_ENE_DIFFUSION')
  call cgyro_readbc_int(collision_kperp,'COLLISION_KPERP')
  call cgyro_readbc_int(collision_field_model,'COLLISION_FIELD_MODEL')
  call cgyro_readbc_int(collision_ion_model,'COLLISION_ION_MODEL')
  call cgyro_readbc_int(collision_precision_mode,'COLLISION_PRECISION_MODE')
  call cgyro_readbc_int(collision_test_mode,'COLLISION_TEST_MODE')
  call cgyro_readbc_int(collision_field_max_l,'COLLISION_FIELD_MAX_L')
  call cgyro_readbc_int(collision_test_max_l,'COLLISION_TEST_MAX_L')
  call cgyro_readbc_real(z_eff)
  call cgyro_readbc_int(z_eff_method,'Z_EFF_METHOD')
  call cgyro_readbc_int(zf_test_mode,'ZF_TEST_MODE')
  call cgyro_readbc_int(nonlinear_flag,'NONLINEAR_FLAG')
  call cgyro_readbc_int(ae_flag,'AE_FLAG')
  call cgyro_readbc_real(temp_ae)
  call cgyro_readbc_real(dens_ae)
  call cgyro_readbc_real(mass_ae)
  call cgyro_readbc_real(dlntdr_ae)
  call cgyro_readbc_real(dlnndr_ae)
  call cgyro_readbc_real(lambda_star)
  call cgyro_readbc_int(h_print_flag,'H_PRINT_FLAG')
  call cgyro_readbc_int(moment_print_flag,'MOMENT_PRINT_FLAG')
  call cgyro_readbc_int(gflux_print_flag,'GFLUX_PRINT_FLAG')
  call cgyro_readbc_int(field_print_flag,'FIELD_PRINT_FLAG')
  call cgyro_readbc_int(triad_print_flag,'TRIAD_PRINT_FLAG')
  call cgyro_readbc_real(amp0)
  call cgyro_readbc_real(amp)
  call cgyro_readbc_real(gamma_e)
  call cgyro_readbc_real(gamma_p)
  call cgyro_readbc_real(mach)
  call cgyro_readbc_int(rotation_model,'ROTATION_MODEL')
  call cgyro_readbc_int(nt_loc,'TOROIDALS_PER_PROC')
  call cgyro_readbc_int(mpi_rank_order,'MPI_RANK_ORDER')
  call cgyro_readbc_int(velocity_order,'VELOCITY_ORDER')
  call cgyro_readbc_int(hiprec_flag,'HIPREC_FLAG')
  call cgyro_readbc_int(udsymmetry_flag,'UDSYMMETRY_FLAG')
  call cgyro_readbc_int(shear_method,'SHEAR_METHOD')
  call cgyro_readbc_int(global_flag,'GLOBAL_FLAG')
  call cgyro_readbc_int(n_global,'N_GLOBAL')
  call cgyro_readbc_real(nu_global)
  call cgyro_readbc_int(theta_plot,'THETA_PLOT')
  call cgyro_readbc_int(gpu_bigmem_flag,'GPU_BIGMEM_FLAG')
  call cgyro_readbc_int(upwind_single_flag,'UPWIND_SINGLE_FLAG')
  call cgyro_readbc_real(px0)
  call cgyro_readbc_int(stream_term,'STREAM_TERM')
  call cgyro_readbc_real(stream_factor)
  call cgyro_readbc_int(exch_flag,'EXCH_FLAG')
  call cgyro_readbc_real(res_weight_power)

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
  call cgyro_readbc_int(n_species,'N_SPECIES')
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
  call cgyro_readbc_realv(sbeta,is)
  call cgyro_readbc_realv(dlnndr_scale,is)
  call cgyro_readbc_realv(dlntdr_scale,is)

  call cgyro_readbc_int(quasineutral_flag,'QUASINEUTRAL_FLAG')
  call cgyro_readbc_real(lambda_star_scale)
  call cgyro_readbc_real(gamma_e_scale)
  call cgyro_readbc_real(gamma_p_scale)
  call cgyro_readbc_real(mach_scale)
  call cgyro_readbc_real(beta_star_scale)
  call cgyro_readbc_real(betae_unit_scale)
  call cgyro_readbc_real(nu_ee_scale)
  call cgyro_readbc_real(zf_scale)

  if (i_proc == 0) close(1)

end subroutine cgyro_read_input

!------------------------------------------------------------
! Service routines: 
!
! (1) read and broadcast an integer:
!
subroutine cgyro_readbc_int(p,label)

  use mpi
  use cgyro_globals, only : i_proc,i_err,CGYRO_COMM_WORLD

  implicit none
  integer, intent(inout) :: p
  character (len=*), intent(in) :: label
  character (len=40) :: actual_label

  if (i_proc == 0) then
       read(1,*) p,actual_label
       if (label /= actual_label) then
          write(*,*) 'ERROR: Invalid label found in input file', actual_label, ' != ', label
          call MPI_ABORT(CGYRO_COMM_WORLD,1,i_err)
          STOP 'Invalid label found in input file'
       endif
  endif

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
