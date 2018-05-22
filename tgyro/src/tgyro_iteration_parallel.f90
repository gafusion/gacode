!-----------------------------------------------------------
! tgyro_iteration_parallel.f90
!
! PURPOSE:
!  Main driver for parallel "blocked" solver.  This is the
!  same algorithm as in tgyro_iteration_serial.
!----------------------------------------------------------

subroutine tgyro_iteration_parallel

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables

  implicit none

  allocate(res1(p_max))
  allocate(x_vec1(p_max)) 
  allocate(f_vec1(p_max))
  allocate(g_vec1(p_max))

  allocate(res_vec(p_max,n_worker))
  allocate(x_vec_vec(p_max,n_worker)) 
  allocate(f_vec_vec(p_max,n_worker))
  allocate(g_vec_vec(p_max,n_worker))

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

     x_vec = x_vec0
     call tgyro_target_vector(x_vec,g_vec)

     ! Worker index tells which function to evaluate 
     call tgyro_flux_vector(x_vec,f_vec,dx,worker_index)

     if (worker == 0) f_vec0 = f_vec
     call MPI_BCAST(f_vec0,size(f_vec0),MPI_DOUBLE_PRECISION,0,gyro_rad,ierr)

     if (worker == 0) call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
     call MPI_BCAST(res,size(res),MPI_DOUBLE_PRECISION,0,gyro_rad,ierr)

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

     !----------------------------------------------
     ! Block diagonal matrix
     !
     ! (p  ,p) (p  ,p+1) (p  ,p+2)
     ! (p+1,p) (p+1,p+1) (p+1,p+2)
     ! (p+2,p) (p+2,p+1) (p+2,p+2)

     jf(:,:) = 0.0

     call MPI_ALLGATHER(f_vec,size(f_vec),MPI_DOUBLE_PRECISION,&
          f_vec_vec,size(f_vec),MPI_DOUBLE_PRECISION,gyro_rad,ierr)

     do i_worker=2,n_worker 
        ip = i_worker-2
        do p=1,p_max,n_evolve
           do pp=0,n_evolve-1
              jf(p+pp,p+ip) = (f_vec_vec(p+pp,i_worker)-f_vec0(p+pp))/dx
           enddo
        enddo
     enddo
     !
     !----------------------------------------------

     !-----------------------------------------------------------
     ! END FLUX JACOBIAN
     !-----------------------------------------------------------

     !----------------------------------------------
     ! Total Jacobian: (dQ/dz-dQ^T/dz)
     !
     jfg(:,:) = jf(:,:)-jg(:,:)
     !----------------------------------------------

     !----------------------------------------------------
     ! Compute actual-target: Relaxation is added to move 
     ! less aggressively to target solution, f0=g0.
     !
     b(:) = -(f_vec0(:)-g_vec0(:))
     !----------------------------------------------------

     ! LAPACK matrix factorization into L/U components
     call DGETRF(p_max,p_max,jfg,p_max,ipiv,ierr) 

     ! LAPACK matrix solve (jfg)x=b (pass jfg and b, return x in b).
     call DGETRS('N',p_max,1,jfg,p_max,ipiv,b,p_max,ierr)

     if (ierr < 0) then
        call tgyro_catch_error('ERROR: DGETRS failed in tgyro_iteration_parallel')
     endif

     ! Check to see if step length exceeds maximum 
     do p=1,p_max
        if (abs(b(p)) > loc_dx_max/r_min) then
           b(p) = sign(loc_dx_max/r_min,b(p))
           b_flag(p) = '*'
        endif
     enddo

     x_vec = x_vec0

     ! Each worker (in parallel) gets a different search vector
     ! NOTE: this is done since we have extra workers already
     !       assigned to do parallel Jacobian
     
     x_vec1(:) = x_vec0(:)+b(:)*search(worker+1,search_index)
     call tgyro_target_vector(x_vec1,g_vec1)
     call tgyro_flux_vector(x_vec1,f_vec1,0.0,0)
     call tgyro_residual(f_vec1,g_vec1,res1,p_max,loc_residual_method)

     call MPI_ALLGATHER(res1,size(res1),MPI_DOUBLE_PRECISION,&
          res_vec,size(res1),MPI_DOUBLE_PRECISION,gyro_rad,ierr)
     call MPI_ALLGATHER(f_vec1,size(f_vec1),MPI_DOUBLE_PRECISION,&
          f_vec_vec,size(f_vec1),MPI_DOUBLE_PRECISION,gyro_rad,ierr)
     call MPI_ALLGATHER(g_vec1,size(g_vec1),MPI_DOUBLE_PRECISION,&
          g_vec_vec,size(g_vec1),MPI_DOUBLE_PRECISION,gyro_rad,ierr)
     call MPI_ALLGATHER(x_vec1,size(x_vec1),MPI_DOUBLE_PRECISION,&
          x_vec_vec,size(x_vec1),MPI_DOUBLE_PRECISION,gyro_rad,ierr)

     correct_flag = 0
     do i_worker=1,n_worker
        res1   = res_vec(:,i_worker)
        x_vec1 = x_vec_vec(:,i_worker)
        f_vec1 = f_vec_vec(:,i_worker)
        g_vec1 = g_vec_vec(:,i_worker)
        call tgyro_write_intermediate(i_worker,res1)
        if (sum(res1) < sum(res)) then
           res = res1
           x_vec = x_vec1
           f_vec = f_vec1
           call tgyro_target_vector(x_vec,g_vec)
           correct_flag = 1
           relax(:) = search(i_worker,search_index)
        endif
     enddo

     if (correct_flag == 1) then
        search_index = 1
     else
        search_index = search_index+1
        if (search_index > search_max) then
           error_flag = 1
           error_msg  = 'ERROR: convergence failure'
        else
           relax(:) = 0.0
        endif
     endif

     ! Write current profiles, gradients, fluxes, targets.
     call tgyro_write_data(1)

     if (error_flag == 1) exit

  enddo

end subroutine tgyro_iteration_parallel
