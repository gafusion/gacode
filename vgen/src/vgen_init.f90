subroutine vgen_init

  use mpi
  use vgen_globals
  use neo_interface
  use expro

  implicit none
  integer :: j

  call neo_init_serial(path)
  call neo_read_input()
  call map_global2interface()

  neo_n_radial_in = 1
  neo_profile_model_in = 1
  
  ! Set qn option for EXPRO

  EXPRO_ctrl_quasineutral_flag = 1
  EXPRO_ctrl_n_ion = neo_n_species_in + neo_ae_flag_in - 1
  n_ions = EXPRO_ctrl_n_ion
  
  ! Set equilibrium option for EXPRO

  if (neo_equilibrium_model_in == 3) then
     EXPRO_ctrl_numeq_flag = 1
  else
     EXPRO_ctrl_numeq_flag = 0
  endif

  call expro_read('input.gacode')
   
  if (EXPRO_error == 1) then
     if (i_proc == 0) then
        print '(a)', 'ERROR: (VGEN) Negative ion density'
     endif
     call MPI_finalize(i_err)
     stop
  endif

  ! Set mass and charge from EXPRO
  do j=1, EXPRO_ctrl_n_ion
     neo_z_in(j)    = expro_z(j)
     neo_mass_in(j) = expro_mass(j)/2.0
  enddo
  if(neo_ae_flag_in == 0) then
     neo_z_in(neo_n_species_in)    = expro_ze
     neo_mass_in(neo_n_species_in) = expro_masse/2.0
  endif
  
  ! Set sign of btccw and ipccw from sign of b and q from EXPRO
  neo_btccw_in = -EXPRO_signb
  neo_ipccw_in = -EXPRO_signb*EXPRO_signq

  ! Check ion velocity ix index
  
  if (erspecies_indx > EXPRO_ctrl_n_ion) then
     if (i_proc == 0) then
        print '(a)', 'ERROR: (VGEN) Invalid species index'
     endif
     call MPI_finalize(i_err)
     stop
  endif
  
  ! Set nn option for neoclassical solution
  if (nn_flag == 1) then
     neo_sim_model_in = 4
     ! Presently only computes jpar
     if(er_method /= 4) then
        if (i_proc == 0) then
           print '(a)','ERROR: (VGEN) NEO NN requires er_method=4'
        endif
        call MPI_finalize(i_err)
        stop
     endif
  endif
  
  ! Storage 
  allocate(vtor_measured(EXPRO_n_exp))
  allocate(pflux_sum(EXPRO_n_exp))
  allocate(jbs_neo(EXPRO_n_exp))
  allocate(jbs_sauter(EXPRO_n_exp))
  allocate(jsigma_neo(EXPRO_n_exp))
  allocate(jsigma_sauter(EXPRO_n_exp))
  allocate(jtor_neo(EXPRO_n_exp))
  allocate(jtor_sauter(EXPRO_n_exp))
  
  pflux_sum(:)     = 0.0
  jbs_neo(:)       = 0.0
  jbs_sauter(:)    = 0.0
  jsigma_neo(:)    = 0.0
  jsigma_sauter(:) = 0.0
  jtor_neo(:)      = 0.0
  jtor_sauter(:)   = 0.0

  do j=1,EXPRO_n_ion
     if (j == erspecies_indx) then
        vtor_measured(:) = EXPRO_vtor(j,:)
     endif
     EXPRO_vpol(j,:) = 0.0
     EXPRO_vtor(j,:) = 0.0
  enddo
  
end subroutine vgen_init
