  subroutine vgen_init
    use mpi
    use vgen_globals
    use neo_interface
    use EXPRO_interface

    implicit none
    integer :: num_ele
    integer :: indx_ele
    integer :: j
    
    call neo_init_serial(path)
    call neo_read_input()
    call map_global2interface()
    
    neo_n_radial_in = 1
    neo_profile_model_in = 1
    !neo_sim_model_in = 1
    
    zfac(:) = 0
    do j=1,neo_n_species_in
       zfac(j)=neo_z_in(j)
    enddo
    
    
    ! Species checks
    num_ele = 0
    indx_ele = 0
    do j=1,neo_n_species_in
       if (zfac(j) < 0.0) then
          num_ele  = num_ele + 1
          indx_ele = j
       endif
    enddo
    
    if (num_ele == 0) then
       n_ions = neo_n_species_in
    else if (num_ele == 1) then
       n_ions = neo_n_species_in - 1
       if (indx_ele /= neo_n_species_in) then
          if (i_proc == 0) then
             print '(a)','ERROR: (VGEN) Electron species must be n_species'
          endif
          call MPI_finalize(i_err)
          stop
       endif
    else
       if (i_proc == 0) then
          print '(a)', 'ERROR: (VGEN) Only one electron species allowed'
       endif
       call MPI_finalize(i_err)
       stop
    endif
    
    if (n_ions < 1) then
       if (i_proc == 0) then
          print '(a)', 'ERROR: (VGEN) There must be at least one ion species'
       endif
       call MPI_finalize(i_err)
       stop
    endif
    
    if (erspecies_indx > n_ions) then
       if (i_proc == 0) then
          print '(a)', 'ERROR: (VGEN) Invalid species index'
       endif
       call MPI_finalize(i_err)
       stop
    endif

    ! Read experimental profiles using EXPRO library.

  call EXPRO_palloc(MPI_COMM_WORLD,path,1)

  ! Reset species 1 density for quasineutrality
 
  EXPRO_ctrl_quasineutral_flag = 1
  EXPRO_ctrl_z(:) = 0.0
  EXPRO_ctrl_n_ion = 0
  do j=1,neo_n_species_in
     if (zfac(j) > 0.0) then
        EXPRO_ctrl_z(j) = 1.0*zfac(j)
        EXPRO_ctrl_n_ion = EXPRO_ctrl_n_ion+1
     endif
  enddo
  
  ! Set equilibrium option for EXPRO
  
  if (neo_equilibrium_model_in == 3) then
     EXPRO_ctrl_numeq_flag = 1
  else
     EXPRO_ctrl_numeq_flag = 0
  endif

  call EXPRO_pread

  ! Write the derived quantities to input.profiles.extra

  if (i_proc == 0) call EXPRO_write_derived(1,'input.profiles.extra')

  ! Set sign of btccw and ipccw from sign of b and q from EXPRO
  neo_btccw_in = -EXPRO_signb
  neo_ipccw_in = -EXPRO_signb*EXPRO_signq
 
  ! Storage 
  allocate(vtor_measured(EXPRO_n_exp))
  allocate(jbs_neo(EXPRO_n_exp))
  allocate(jbs_sauter(EXPRO_n_exp))
  allocate(jbs_koh(EXPRO_n_exp))
  allocate(jbs_nclass(EXPRO_n_exp))
  allocate(jtor_neo(EXPRO_n_exp))
  allocate(jtor_sauter(EXPRO_n_exp))
  allocate(pflux_sum(EXPRO_n_exp))
  jbs_neo(:)     = 0.0
  jbs_sauter(:)  = 0.0
  jbs_koh(:)     = 0.0
  jbs_nclass(:)  = 0.0
  jtor_neo(:)    = 0.0
  jtor_sauter(:) = 0.0
  pflux_sum(:)   = 0.0
  
  do j=1,10
     if (j == erspecies_indx) then
        vtor_measured(:) = EXPRO_vtor(j,:)
     endif
     EXPRO_vpol(j,:) = 0.0
     EXPRO_vtor(j,:) = 0.0
  enddo
  
end subroutine vgen_init
