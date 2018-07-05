subroutine EXPRO_locsim_alloc(flag)

  use EXPRO_locsim_interface
  use EXPRO_interface

  implicit none
  integer, intent (in) :: flag  

  ! flag=0: deallocate
  ! flag=1: allocate

  if (flag == 1) then

     allocate(rmin_exp(EXPRO_n_exp))

     allocate(temp_exp(n_species_exp,EXPRO_n_exp))
     allocate(dens_exp(n_species_exp,EXPRO_n_exp))
     allocate(dlntdr_exp(n_species_exp,EXPRO_n_exp))
     allocate(dlnndr_exp(n_species_exp,EXPRO_n_exp))
     allocate(sdlntdr_exp(n_species_exp,EXPRO_n_exp))
     allocate(sdlnndr_exp(n_species_exp,EXPRO_n_exp))

     allocate(gamma_e_exp(EXPRO_n_exp))
     allocate(gamma_p_exp(EXPRO_n_exp))
     allocate(mach_exp(EXPRO_n_exp))

  else

     if(allocated(rmin_exp))        deallocate(rmin_exp)

     if(allocated(temp_exp))        deallocate(temp_exp)
     if(allocated(dens_exp))        deallocate(dens_exp)
     if(allocated(dlntdr_exp))      deallocate(dlntdr_exp)
     if(allocated(dlnndr_exp))      deallocate(dlnndr_exp)
     if(allocated(sdlntdr_exp))     deallocate(sdlntdr_exp)
     if(allocated(sdlnndr_exp))     deallocate(sdlnndr_exp)
     
     if(allocated(gamma_e_exp))     deallocate(gamma_e_exp)
     if(allocated(gamma_p_exp))     deallocate(gamma_p_exp)
     if(allocated(mach_exp))        deallocate(mach_exp)
     
     if(allocated(geo_yin_exp))     deallocate(geo_yin_exp)

  endif

end subroutine EXPRO_locsim_alloc
