!-----------------------------------------------------------
! tgyro_iteration_serial.f90
!
! PURPOSE:
!  Main driver for serial "blocked" solver.  This is the
!  same algorithm as in tgyro_iteration_parallel.
!----------------------------------------------------------

subroutine tgyro_iteration_serial

  use tgyro_globals
  use tgyro_iteration_variables

  implicit none

  allocate(res1(p_max))
  allocate(x_vec1(p_max))
  allocate(f_vec1(p_max))
  allocate(g_vec1(p_max))

  do i_tran_loop=1,tgyro_relax_iterations

     i_tran = i_tran+1

     x_vec0 = x_vec
     call tgyro_target_vector(x_vec,g_vec)
     g_vec0 = g_vec

     !-----------------------------------------------------------
     ! BEGIN TARGET (SOURCE) JACOBIAN: dQ^T/dz
     !-----------------------------------------------------------

     do p=1,p_max

        x_vec(:) = x_vec0(:)
        x_vec(p) = x_vec0(p)+dx

        call tgyro_target_vector(x_vec,g_vec)
        jg(:,p) = (g_vec(:)-g_vec0(:))/dx

     enddo

     !-----------------------------------------------------------
     ! END TARGET (SOURCE) JACOBIAN
     !-----------------------------------------------------------

     !-----------------------------------------------------------
     ! BEGIN FLUX JACOBIAN: dQ/dz
     !-----------------------------------------------------------

     ! Reset gradients
     x_vec = x_vec0
     
     call tgyro_target_vector(x_vec,g_vec)
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     
     f_vec0 = f_vec

     call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
     call tgyro_write_intermediate(0,res)

     !-------------------------------------------------------------
     ! If this is the first iteration, write data for the ZEROTH 
     ! iteration.
     !
     if (i_tran_loop == 1 .and. loc_restart_flag == 0) then
        i_tran = 0
        call tgyro_write_data(1)
        i_tran = 1
     endif
     !-------------------------------------------------------------

     !  Build dQ/dz (block diagonal matrix)
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

     !-----------------------------------------------------------
     ! END FLUX JACOBIAN
     !-----------------------------------------------------------

     !----------------------------------------------
     ! Total Jacobian: (dQ/dz-dQ^T/dz)
     !
     jfg(:,:) = jf(:,:)-jg(:,:)
     !----------------------------------------------

     !----------------------------------------------------
     ! Compute actual-target. Later a search vector will
     ! be applied to the solultions so as to move less
     ! than a full Newron step.
     !
     b(:) = -(f_vec0(:)-g_vec0(:))
     !----------------------------------------------------

     ! LAPACK matrix factorization into L/U components
     call DGETRF(p_max,p_max,jfg,p_max,ipiv,ierr) 

     ! LAPACK matrix solve (jfg)x=b (pass jfg and b, return x in b).
     call DGETRS('N',p_max,1,jfg,p_max,ipiv,b,p_max,ierr)

     if (ierr < 0) then
        call tgyro_catch_error('ERROR: (tgyro_iteration_serial) DGETRS failed.')
     endif

     ! Check to see if step length exceeds maximum 
     do p=1,p_max
        if (abs(b(p)) > loc_dx_max/r_min) then
           b(p) = sign(loc_dx_max/r_min,b(p))
           b_flag(p) = '*'
        endif
     enddo

     x_vec = x_vec0

     correct_flag = 0
     
     ! Each worker (serial) gets a different test vector
     ! NOTE: this is done to match the parallel algorithm, where extra
     !       parallel workers are assigned to do parallel Jacobian
     do i_worker=1,n_evolve+1

        ! Update gradient using search vector
        x_vec1(:) = x_vec0(:)+b(:)*search(i_worker,search_index)
        call tgyro_target_vector(x_vec1,g_vec1)
        call tgyro_flux_vector(x_vec1,f_vec1,0.0,0)
        call tgyro_residual(f_vec1,g_vec1,res1,p_max,loc_residual_method)
        call tgyro_write_intermediate(i_worker,res1)

        ! Test to see if GLOBAL residual is reduced
        if (sum(res1) < sum(res)) then
           res   = res1
           x_vec = x_vec1
           f_vec = f_vec1
           call tgyro_target_vector(x_vec,g_vec)
           correct_flag = 1
           relax(:) = search(i_worker,search_index)
        endif
        
     enddo

     if (correct_flag == 1) then
        ! The correction step reduced the residual, so reset the index.
        search_index = 1
     else
        ! Residual not reduced, so try different search length 
        search_index = search_index+1
        if (search_index > search_max) then
           error_flag = 1
           error_msg  = 'ERROR: (tgyro_iteration_serial) convergence failure'
        else
           relax(:) = 0.0
        endif
     endif

     ! Write current profiles, gradients, fluxes, targets.
     call tgyro_write_data(1)

     if (error_flag == 1) exit

  enddo


end subroutine tgyro_iteration_serial
