subroutine tgyro_profile_reintegrate

  use tgyro_globals
  use tgyro_ped
  use tgyro_iteration_variables

  implicit none

  integer :: i_ion
  integer :: i_star  
  real :: w

  if (tgyro_ped_model > 1) then

     ! Map data over r(n_r) < r < a

     call tgyro_pedestal_map(dlnnedr(n_r),zn_top,n_top(1)*n_frac,nn_vec(:,2)*n_frac,&
          i_star,exp_ne)
     call tgyro_pedestal_map(dlntedr(n_r),zt_top,t_top(1)*t_frac,t_vec(:)*t_frac,&
          i_star,exp_te)

     ! Map ion densities
     ! NOTE: Assumption is that ion profiles are 
     !
     !    ni(j,r) = ne(r)*ratio(j) for r > r(n_r)
     !
     ! where ratio(j) is the pivot density ratio at t=0 (see tgyro_pedestal).
     !
     do i_ion=1,loc_n_ion
        if (therm_flag(i_ion) == 1) then
           w = n_ratio(i_ion)*n_frac
           call tgyro_pedestal_map(dlnnidr(i_ion,n_r),zn_top,w*n_top(1),w*nn_vec(:,2),&
                i_star,exp_ni(i_ion,:))
           w = t_ratio(i_ion)*t_frac
           call tgyro_pedestal_map(dlntidr(i_ion,n_r),zt_top,w*t_top(1),w*t_vec(:),&
                i_star,exp_ti(i_ion,:))
        endif
     enddo

  endif

  ! Map data inside r < r(n_r)

  call tgyro_expro_map(r,dlnnedr,n_r,ne(n_r),rmin_exp,exp_ne,n_exp,'log')
  call tgyro_expro_map(r,dlntedr,n_r,te(n_r),rmin_exp,exp_te,n_exp,'log')
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        call tgyro_expro_map(r,dlnnidr(i_ion,:),n_r,ni(i_ion,n_r),rmin_exp,exp_ni(i_ion,:),n_exp,'log')
        call tgyro_expro_map(r,dlntidr(i_ion,:),n_r,ti(i_ion,n_r),rmin_exp,exp_ti(i_ion,:),n_exp,'log')
     endif
  enddo
  call tgyro_expro_map(r,w0p,n_r,w0(n_r),rmin_exp,exp_w0,n_exp,'lin')

  ! Compute pressure: ptot_exp
  call tgyro_pressure

end subroutine tgyro_profile_reintegrate
