!-----------------------------------------------------------
! tgyro_iteration_standard.f90
!
! PURPOSE:
!  Control of original and related iteration schemes.
!  For complementary methods, see tgyro_iteration_pppl.f90
!----------------------------------------------------------

subroutine tgyro_iteration_standard

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables

  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,'(t2,a)') 'INFO: (TGYRO) Beginning standard iterations'
     close(1)
  endif

  !---------------------------------------------------
  ! ZEROTH ITERATION
  !
  ! One pass to get fluxes.
  !
  ! We *do* want the iteration number to increase in 
  ! the case where we restart without iterating (this
  ! can be used to compare transport models).
  !  
  if (tgyro_relax_iterations == 0 .and. loc_restart_flag == 1) i_tran=i_tran+1
  !
  call tgyro_target_vector(x_vec,g_vec)
  if (loc_restart_flag == 0 .or. tgyro_relax_iterations == 0) then
     ! Need to determine initial fluxes
     gyro_restart_method = 1
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     gyro_restart_method = 2
  else
     ! Initial fluxes already computed
     call tgyro_flux_set(f_vec)
     ! GYRO restart data available
     gyro_restart_method = 2
  endif
  res0 = 0.0
  if (gyrotest_flag == 0) call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)

  if (loc_restart_flag == 0 .or. tgyro_relax_iterations == 0) then
     call tgyro_write_data(1)
  endif
  !----------------------------------------------------

  do i_tran_loop=1,tgyro_relax_iterations

     i_tran = i_tran+1

     ! Initialize gradients

     x_vec0 = x_vec
     g_vec0 = g_vec
     f_vec0 = f_vec
     res0   = res

     !----------------------------------------------
     ! Build dQ^T/dz (dense matrix)
     !
     do p=1,p_max

        x_vec(:) = x_vec0(:)
        x_vec(p) = x_vec0(p)+dx

        call tgyro_target_vector(x_vec,g_vec)
        jg(:,p) = (g_vec(:)-g_vec0(:))/dx

     enddo

     ! Reset gradients
     x_vec = x_vec0

     ! Reset profiles to be consistent with gradient.
     call tgyro_profile_set(x_vec,0.0,0)
     call tgyro_profile_functions 
     !----------------------------------------------

     !----------------------------------------------
     ! Build dQ/dz (block diagonal matrix)
     !
     ! (p  ,p) (p  ,p+1) (p  ,p+2)
     ! (p+1,p) (p+1,p+1) (p+1,p+2)
     ! (p+2,p) (p+2,p+1) (p+2,p+2)

     jf(:,:) = 0.0

     do ip=0,n_evolve-1
        call tgyro_flux_vector(x_vec,f_vec,dx,evolve_indx(ip+1))
        do p=1,p_max,n_evolve
           do pp=0,n_evolve-1
              jf(p+pp,p+ip) = (f_vec(p+pp)-f_vec0(p+pp))/dx
           enddo
        enddo
     enddo
     !
     !----------------------------------------------

     !----------------------------------------------
     ! Total Jacobian: (dQ/dz-dQ^T/dz)
     !
     jfg(:,:) = jf(:,:)-jg(:,:)
     !----------------------------------------------

     !----------------------------------------------------
     ! Compute target.  Relaxation is added to move less
     ! aggressively to target solution, f0=g0.
     !
     b(:) = -(f_vec0(:)-g_vec0(:))*relax(:)
     !----------------------------------------------------

     ! LAPACK matrix factorization into L/U components
     call DGETRF(p_max,p_max,jfg,p_max,ipiv,ierr) 

     ! LAPACK matrix solve (jfg)x=b (pass jfg and b, return x in b).
     call DGETRS('N',p_max,1,jfg,p_max,ipiv,b,p_max,ierr)

     if (ierr < 0) then
        call tgyro_catch_error('ERROR: DGETRS failed in tgyro_iteration_standard')
     endif

     ! Check to see if step length exceeds maximum 
     do p=1,p_max
        if (abs(b(p)) > loc_dx_max/r_min) then
           b(p) = sign(loc_dx_max/r_min,b(p))
           b_flag(p) = '*'
        endif
     enddo

     !----------------------------------------------------
     ! Update gradient using Newton-step.
     !
     x_vec(:) = x_vec0(:)+b(:)
     !----------------------------------------------------

     !-----------------------------------------------------
     ! Correction step:
     !  strategy to cope with an increasing residual.
     ! 
     call tgyro_target_vector(x_vec,g_vec)
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     !
     ! Compute initial residual
     call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
     call tgyro_write_intermediate(0,res)
     !
     do p=1,p_max
        if (res0(p) < res(p) .and. loc_relax > 1.0) then

           correct_flag = 1

           ! Correct solution vector and try relaxation
           x_vec(p) = x_vec0(p)
           relax(p) = relax(p)/loc_relax

           ! If relaxation gets too small, try large value.
           if (relax(p) < 1/loc_relax**3) then
              relax(p) = 0.75*loc_relax
              x_vec(p) = x_vec0(p)+2*b(p)
           endif
        else

           ! Reset relaxation
           relax(p) = 1.0

        endif
     enddo
     !
     if (correct_flag == 1) then

        ! Recompute solution
        call tgyro_target_vector(x_vec,g_vec)
        gyro_restart_method = 1
        call tgyro_flux_vector(x_vec,f_vec,0.0,0)
        gyro_restart_method = 2

        ! Recompute residual
        call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
        call tgyro_write_intermediate(1,res)

        correct_flag = 0

     endif
     !----------------------------------------------------- 

     call tgyro_write_data(1)

  enddo

end subroutine tgyro_iteration_standard
