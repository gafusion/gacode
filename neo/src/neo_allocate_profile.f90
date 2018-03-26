module neo_allocate_profile

  implicit none

  public :: PROFILE_SIM_alloc

  logical, private :: initialized_sim = .false.

contains

  !  Allocate/Deallocate simulation-grid profile functions.
  subroutine PROFILE_SIM_alloc(flag)

    use neo_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if(flag == 1) then
       if(initialized_sim) return

       allocate(r(n_radial))
       allocate(rmaj(n_radial))
       allocate(dphi0dr(n_radial))
       allocate(epar0(n_radial))
       allocate(q(n_radial))
       allocate(rho(n_radial))   
       allocate(shear(n_radial))
       allocate(shift(n_radial))
       allocate(kappa(n_radial))
       allocate(s_kappa(n_radial))
       allocate(delta(n_radial))
       allocate(s_delta(n_radial))
       allocate(zeta(n_radial))
       allocate(s_zeta(n_radial))
       allocate(zmag(n_radial))
       allocate(s_zmag(n_radial))
       allocate(beta_star(n_radial))
       allocate(te_ade(n_radial))
       allocate(ne_ade(n_radial))
       allocate(dlntdre_ade(n_radial))
       allocate(dlnndre_ade(n_radial))
       allocate(dens(n_species,n_radial))
       allocate(temp(n_species,n_radial))
       allocate(vth(n_species,n_radial))
       allocate(dlnndr(n_species,n_radial))
       allocate(dlntdr(n_species,n_radial))
       allocate(nu(n_species,n_radial))
       allocate(z(n_species))
       allocate(mass(n_species))
       allocate(aniso_model(n_species))
       allocate(vth_para(n_species,n_radial))
       allocate(temp_para(n_species,n_radial))
       allocate(dlntdr_para(n_species,n_radial))
       allocate(temp_perp(n_species,n_radial))
       allocate(dlntdr_perp(n_species,n_radial))
       allocate(b_unit(n_radial))
       allocate(dens_norm(n_radial))
       allocate(temp_norm(n_radial))
       allocate(vth_norm(n_radial))
       allocate(omega_rot(n_radial))
       allocate(omega_rot_deriv(n_radial))
       allocate(psiN_polflux(n_radial))
       allocate(rhoN_torflux(n_radial))
       
       geo_numeq_flag = 0
       geo_ny = 0
       allocate(geo_yin(8,0:geo_ny,n_radial))
       geo_yin(:,:,:) = 0.0

       initialized_sim = .true.

    else

       if(.NOT. initialized_sim) return

       deallocate(r)
       deallocate(rmaj)
       deallocate(dphi0dr)
       deallocate(epar0)
       deallocate(q)
       deallocate(rho)   
       deallocate(shear)
       deallocate(shift)
       deallocate(kappa)
       deallocate(s_kappa)
       deallocate(delta)
       deallocate(s_delta)
       deallocate(zeta)
       deallocate(s_zeta)
       deallocate(zmag)
       deallocate(s_zmag)
       deallocate(beta_star)
       deallocate(te_ade)
       deallocate(ne_ade)
       deallocate(dlntdre_ade)
       deallocate(dlnndre_ade)
       deallocate(dens)
       deallocate(temp)
       deallocate(vth)
       deallocate(dlnndr)
       deallocate(dlntdr)
       deallocate(nu)
       deallocate(z)
       deallocate(mass)
       deallocate(aniso_model)
       deallocate(vth_para)
       deallocate(temp_para)
       deallocate(dlntdr_para)
       deallocate(temp_perp)
       deallocate(dlntdr_perp)
       deallocate(b_unit)
       deallocate(dens_norm)
       deallocate(temp_norm)
       deallocate(vth_norm)
       deallocate(omega_rot)
       deallocate(omega_rot_deriv)
       deallocate(psiN_polflux)
       deallocate(rhoN_torflux)
       
       deallocate(geo_yin)

       initialized_sim = .false.

    endif

  end subroutine PROFILE_SIM_alloc

end module neo_allocate_profile
