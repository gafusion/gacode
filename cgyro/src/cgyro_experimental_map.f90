!-----------------------------------------------------------
! cgyro_experimental_map.f90 
!
! PURPOSE:
!  Generate radial profile functions on the simulation 
!  (slice) r-grid.
!-----------------------------------------------------------

subroutine cgyro_experimental_map

  use cgyro_globals
  use cgyro_experimental_globals

  implicit none
  
  integer :: i, j

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_exp) ones.
  ! 
  call cub_spline(rmin_exp,rmaj_exp,n_grid_exp,rmin,rmaj,1)
  call cub_spline(rmin_exp,q_exp,n_grid_exp,rmin,q,1)
  call cub_spline(rmin_exp,s_exp,n_grid_exp,rmin,s,1)
  call cub_spline(rmin_exp,shift_exp,n_grid_exp,rmin,shift,1)
  call cub_spline(rmin_exp,kappa_exp,n_grid_exp,rmin,kappa,1)
  call cub_spline(rmin_exp,s_kappa_exp,n_grid_exp,rmin,s_kappa,1)
  call cub_spline(rmin_exp,delta_exp,n_grid_exp,rmin,delta,1)
  call cub_spline(rmin_exp,s_delta_exp,n_grid_exp,rmin,s_delta,1)
  call cub_spline(rmin_exp,zeta_exp,n_grid_exp,rmin,zeta,1)
  call cub_spline(rmin_exp,s_zeta_exp,n_grid_exp,rmin,s_zeta,1)
  call cub_spline(rmin_exp,zmag_exp,n_grid_exp,rmin,zmag,1)
  call cub_spline(rmin_exp,dzmag_exp,n_grid_exp,rmin,dzmag,1)
  call cub_spline(rmin_exp,gamma_e_exp,n_grid_exp,rmin,gamma_e,1)
  call cub_spline(rmin_exp,gamma_p_exp,n_grid_exp,rmin,gamma_p,1)
  call cub_spline(rmin_exp,mach_exp,n_grid_exp,rmin,mach,1)
  call cub_spline(rmin_exp,rhos_exp,n_grid_exp,rmin,rhos,1)
  call cub_spline(rmin_exp,z_eff_exp,n_grid_exp,rmin,z_eff,1)
  call cub_spline(rmin_exp,b_unit_exp,n_grid_exp,rmin,b_unit,1)

  call cub_spline(rmin_exp,te_ade_exp,n_grid_exp,rmin,te_ade,1)
  call cub_spline(rmin_exp,ne_ade_exp,n_grid_exp,rmin,ne_ade,1)
  call cub_spline(rmin_exp,dlntdre_ade_exp,n_grid_exp,rmin,dlntdre_ade,1)
  call cub_spline(rmin_exp,dlnndre_ade_exp,n_grid_exp,rmin,dlnndre_ade,1)
 
  do i=1,n_species
     ! Note: maping is only done for n_species (not n_species_exp)
     call cub_spline(rmin_exp,dens_exp(i,:),n_grid_exp,rmin,dens(i),1)
     call cub_spline(rmin_exp,temp_exp(i,:),n_grid_exp,rmin,temp(i),1)
     call cub_spline(rmin_exp,dlntdr_exp(i,:),n_grid_exp,rmin,dlntdr(i),1)
     call cub_spline(rmin_exp,dlnndr_exp(i,:),n_grid_exp,rmin,dlnndr(i),1)
     call cub_spline(rmin_exp,sdlntdr_exp(i,:),n_grid_exp,rmin,sdlntdr(i),1)
     call cub_spline(rmin_exp,sdlnndr_exp(i,:),n_grid_exp,rmin,sdlnndr(i),1)
  enddo

  if (geo_numeq_flag == 1) then
     do i=1,8
        do j=0,geo_ny
           call cub_spline(rmin_exp,geo_yin_exp(i,j,:),n_grid_exp,rmin, &
                geo_yin(i,j),1)
        enddo
     enddo
  else
     geo_yin(:,:) = 0.0
  endif

end subroutine cgyro_experimental_map
