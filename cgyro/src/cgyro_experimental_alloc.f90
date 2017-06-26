subroutine cgyro_experimental_alloc(flag)

  use cgyro_globals
  use cgyro_experimental_globals

  implicit none
  integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
  
  if (flag == 1) then

     allocate(rmin_exp(n_grid_exp))
     allocate(rmaj_exp(n_grid_exp))
     allocate(q_exp(n_grid_exp))
     allocate(s_exp(n_grid_exp))
     allocate(shift_exp(n_grid_exp))
     allocate(kappa_exp(n_grid_exp))
     allocate(s_kappa_exp(n_grid_exp))
     allocate(delta_exp(n_grid_exp))
     allocate(s_delta_exp(n_grid_exp))
     allocate(zeta_exp(n_grid_exp))
     allocate(s_zeta_exp(n_grid_exp))
     allocate(zmag_exp(n_grid_exp))
     allocate(dzmag_exp(n_grid_exp))

     allocate(te_ade_exp(n_grid_exp))
     allocate(ne_ade_exp(n_grid_exp))
     allocate(dlntdre_ade_exp(n_grid_exp))
     allocate(dlnndre_ade_exp(n_grid_exp))

     allocate(temp_exp(n_species_exp,n_grid_exp))
     allocate(dens_exp(n_species_exp,n_grid_exp))
     allocate(dlntdr_exp(n_species_exp,n_grid_exp))
     allocate(dlnndr_exp(n_species_exp,n_grid_exp))
     allocate(sdlntdr_exp(n_species_exp,n_grid_exp))
     allocate(sdlnndr_exp(n_species_exp,n_grid_exp))

     allocate(z_eff_exp(n_grid_exp))
     allocate(b_unit_exp(n_grid_exp))
     allocate(gamma_e_exp(n_grid_exp))
     allocate(gamma_p_exp(n_grid_exp))
     allocate(mach_exp(n_grid_exp))
     allocate(rhos_exp(n_grid_exp))

     allocate(geo_yin_exp(8,0:geo_ny,n_grid_exp))
     geo_yin_exp(:,:,:)=0.0

  else

     if(allocated(rmin_exp))        deallocate(rmin_exp)
     if(allocated(rmaj_exp))        deallocate(rmaj_exp)
     if(allocated(q_exp))           deallocate(q_exp)
     if(allocated(s_exp))           deallocate(s_exp)
     if(allocated(shift_exp))       deallocate(shift_exp)
     if(allocated(kappa_exp))       deallocate(kappa_exp)
     if(allocated(s_kappa_exp))     deallocate(s_kappa_exp)
     if(allocated(delta_exp))       deallocate(delta_exp)
     if(allocated(s_delta_exp))     deallocate(s_delta_exp)
     if(allocated(zeta_exp))        deallocate(zeta_exp)
     if(allocated(s_zeta_exp))      deallocate(s_zeta_exp)
     if(allocated(zmag_exp))        deallocate(zmag_exp)
     if(allocated(dzmag_exp))      deallocate(dzmag_exp)

     if(allocated(te_ade_exp))      deallocate(te_ade_exp)
     if(allocated(ne_ade_exp))      deallocate(ne_ade_exp)
     if(allocated(dlntdre_ade_exp)) deallocate(dlntdre_ade_exp)
     if(allocated(dlnndre_ade_exp)) deallocate(dlnndre_ade_exp)

     if(allocated(temp_exp))        deallocate(temp_exp)
     if(allocated(dens_exp))        deallocate(dens_exp)
     if(allocated(dlntdr_exp))      deallocate(dlntdr_exp)
     if(allocated(dlnndr_exp))      deallocate(dlnndr_exp)

     if(allocated(z_eff_exp))       deallocate(z_eff_exp)
     if(allocated(b_unit_exp))      deallocate(b_unit_exp)
     if(allocated(gamma_e_exp))     deallocate(gamma_e_exp)
     if(allocated(gamma_p_exp))     deallocate(gamma_p_exp)
     if(allocated(mach_exp))        deallocate(mach_exp)
     if(allocated(rhos_exp))        deallocate(rhos_exp)

     if(allocated(geo_yin_exp))     deallocate(geo_yin_exp)

  endif


end subroutine cgyro_experimental_alloc
