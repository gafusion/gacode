!-----------------------------------------------------------
! tgyro_iteration_simplerelax.f90
!
! PURPOSE:
! Trial simple relaxation approach to look for more robust 
! solution method.
!----------------------------------------------------------

subroutine tgyro_iteration_simplerelax

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables

  real :: simpledz
  
  quasifix = 0
  res0 = 0.0
  
  do i_tran_loop=1,tgyro_relax_iterations

     i_tran = i_tran+1

     ! Initialize gradients

     x_vec0 = x_vec
     g_vec0 = g_vec
     f_vec0 = f_vec
     res0   = res

     !----------------------------------------------
     ! Relax based on (flux-target)/max(target,1)
     ! dz/z = -loc_relax*(Q_tot - Qtarget)/max(Qtarget,Qtot)

     if (i_tran_loop == 1) then
        relax = loc_relax
     endif

     ! Construct adjustment vector
     do p=1,p_max
        simpledz = relax(p)*(f_vec(p)-g_vec(p))/sqrt(f_vec(p)**2+g_vec(p)**2)                   
        if (abs(simpledz) > loc_dx_max) then
           simpledz = loc_dx_max*simpledz/abs(simpledz)
        endif

        if (mask(p,3) == 1) then
           ! evolving rotation, need to take a bit of extra care
           x_vec(p) = x_vec(p)*(1.0-simpledz*sign(1.0,g_vec(p)))
        else
           ! Method updated 09/2023 to work for x_vec(p) < 0
           x_vec(p) = x_vec(p)-simpledz*abs(x_vec(p))
        endif

     enddo

     ! Update targets/profiles
     call tgyro_target_vector(x_vec,g_vec)
     ! Update fluxes
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     
     ! Compute initial residual
     call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
     call tgyro_write_intermediate(0,res)

     ! Output results
     call tgyro_write_data(1)

  enddo
  
end subroutine tgyro_iteration_simplerelax
