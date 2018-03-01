module neo_neural

  implicit none

  public :: NEURAL_alloc, NEURAL_do, NEURAL_write
  real :: jpar_nn_neo, jtor_nn_neo
  real, dimension(:), allocatable :: pflux_nn_neo, eflux_nn_neo
  real, dimension(:), allocatable :: kbig_upar_nn_neo, vpol_th0_nn_neo
  character(len=80),private :: runfile_nn = 'out.neo.transport_nn'
  integer, parameter, private :: io=1
  logical, private :: initialized = .false.

contains

  subroutine NEURAL_alloc(flag)
    use neo_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if(flag == 1) then
       if(initialized) return

       allocate(pflux_nn_neo(n_species))
       allocate(eflux_nn_neo(n_species))
       allocate(kbig_upar_nn_neo(n_species))
       allocate(vpol_th0_nn_neo(n_species))
       
       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io,file=trim(path)//runfile_nn,status='replace')
          close(io)
       endif

       initialized = .true.

    else
       if(.NOT. initialized) return
       deallocate(pflux_nn_neo)
       deallocate(eflux_nn_neo)
       deallocate(kbig_upar_nn_neo)
       deallocate(vpol_th0_nn_neo)
       initialized = .false.

    endif

  end subroutine NEURAL_alloc
       
  subroutine NEURAL_do(ir)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    
    call compute_nn_flow(ir)
    if(error_status > 0) return
    
  end subroutine NEURAL_do

  ! The NN presently assumes 2 ion species: (D + C) and kinetic electrons.
  ! Strong rotation effects are not included.
  subroutine compute_nn_flow(ir)
    use mpi
    use neo_globals
    use neo_equilibrium
    implicit none
    
    integer, intent (in) :: ir
    integer :: is, k, ierr
    integer :: is_i1, is_i2
    real    :: d_max
    real(4), dimension(18) :: nn_out
    real(4), dimension(6)  :: nn_in
    real, dimension(6) :: xmin = (/ 0.05,0.323,1.0,-2.0,0.6,1.0 /)
    real, dimension(6) :: xmax = (/ 0.35,0.766,10.0,1.0,0.99,3.0 /)
    real, dimension(6) :: C_ln, C_ke, C_ki1, C_ki2
    real :: ke, ki1, ki2
    character(len=218) :: root
    character(len=255) :: data
    
    include 'brainfuse_lib.inc'
    
    ! The nn assumes 2 ion species + electrons
    if(n_species /= 3) then
       call neo_error('ERROR: (NEO) NN requires 2 ion species + electrons')
       return
    endif
    if(adiabatic_ele_model == 1) then
       call neo_error('ERROR: (NEO) NN requires 2 ion species + electrons')
       return
    endif
    
    ! Identify the main ion species
    is_i1 = -1
    d_max = -1.0
    do is=1, n_species
       if(Z(is) > 0 .and. dens(is,ir) > d_max) then
          is_i1 = is
          d_max  = dens(is,ir) 
       endif
    enddo
    if(abs(Z(is_i1) - 1.0) > epsilon(0.)) then
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'WARNING: (NEO) NN assumes Z(ion1)=1.0'
          close(io_neoout)
       endif
    endif
    
    ! Identify the secondary ion species
    do is=1, n_species
       if(is /= is_ele .and. is /= is_i1) then
          is_i2 = is
       endif
    enddo
    if(abs(Z(is_i2) - 6.0) > epsilon(0.)) then
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'WARNING: (NEO) NN assumes Z(ion2)=6.0'
          close(io_neoout)
       endif
    endif

    if(rotation_model == 2) then
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'WARNING: (NEO) NN assumes weak rotation'
          close(io_neoout)
       endif
    endif
    
    ! Set the input parameters for the NN
    ! epsilon=r/R
    nn_in(1) = r(ir)/rmaj(ir)
    ! f_trap
    nn_in(2) = geo_param(ir,2)
    ! |q|
    nn_in(3) = abs(q(ir))                        
    ! log(nuee/cs/R)
    nn_in(4) = log10(nu(is_ele,ir)*rmaj(ir)/sqrt(temp(is_ele,ir)))
    ! n_i1/ne
    nn_in(5) = dens(is_i1,ir)/dens(is_ele,ir)
    ! T_i1/Te
    nn_in(6) = temp(is_i1,ir)/temp(is_ele,ir)      

    ! Re-scale nn_in if out-of-range of training data (except for f_trap)
    do k=1,6
       if(k /= 2) then
          if(nn_in(k) > xmax(k))  then
             nn_in(k) = xmax(k)
          endif
          if(nn_in(k) < xmin(k))  then
             nn_in(k) = xmin(k)
          endif
       endif
    enddo
    
    call get_environment_variable('GACODE_ROOT',root)

    ! Get coefficients computed by the NN 
    ! (1) Cne, (2) Cni1, (3) Cni2, (4) Cte, (5) Cti1, (6) Cti2
    C_ln(1) = dlnndr(is_ele,ir)-1.5*dlntdr(is_ele,ir)
    C_ln(2) = dlnndr(is_i1,ir)-1.5*dlntdr(is_i1,ir)
    C_ln(3) = dlnndr(is_i2,ir)-1.5*dlntdr(is_i2,ir)
    C_ln(4) = dlntdr(is_ele,ir)
    C_ln(5) = dlntdr(is_i1,ir)
    C_ln(6) = dlntdr(is_i2,ir)
    
    ! flow K
    data = trim(root)//'/../neural/neonn/flow_nn/'
    ierr=load_anns(2, trim(data)//char(0),'brainfuse'//char(0))
    if(ierr == 0) then
       call neo_error('ERROR: (NEO) Neural network loading failed.')
       return
    endif
    ierr=load_anns_inputs(nn_in)
    ierr=run_anns()
    ierr=get_anns_avg_array(nn_out)
    if(ierr > 0) then
       call neo_error('ERROR: (NEO) Flow neural network failed.')
       return
    endif
    do k=1,6
       C_ke(k)   = nn_out(k)
       C_ki1(k)  = nn_out(k+6)
       C_ki2(k)  = nn_out(k+12)
    enddo
    ! Reconstruct the flow coefficients from NEO NN: K_a B_unit/(v_norm n_norm)
    ! The NN output is <B^2/Bunit^2> K_a B_unit/(c_s n_a)
    ! ~ - (I/psip) rho_s sum_b (Cnb dlnnb/dr + Ctb dlntb/dr)
    kbig_upar_nn_neo(:) = 0.0
    do k=1,6
       kbig_upar_nn_neo(is_ele)  = kbig_upar_nn_neo(is_ele) &
            + geo_param(ir,1) * rho(ir) / geo_param(ir,3) * temp(is_ele,ir) &
            * dens(is_ele,ir) * C_ke(k)*C_ln(k)
       kbig_upar_nn_neo(is_i1) = kbig_upar_nn_neo(is_i1) &
            + geo_param(ir,1) * rho(ir) / geo_param(ir,3) * temp(is_ele,ir) &
            * dens(is_i1,ir) * C_ki1(k)*C_ln(k)
       kbig_upar_nn_neo(is_i2) = kbig_upar_nn_neo(is_i2) &
            + geo_param(ir,1) * rho(ir) / geo_param(ir,3) * temp(is_ele,ir) &
            * dens(is_i2,ir) * C_ki2(k)*C_ln(k)
    enddo
    vpol_th0_nn_neo(:) = kbig_upar_nn_neo(:)/dens(:,ir)*geo_param(ir,4)
    
    ! Reconstruct jpar from NEO NN 
    jpar_nn_neo = 0.0
    do is=1,n_species
       ! diamagnetic component and neoclassical flow component 
       jpar_nn_neo = jpar_nn_neo + geo_param(ir,1) * rho(ir) &
            * dens(is,ir)*temp(is,ir)* (dlnndr(is,ir) + dlntdr(is,ir)) &
            + geo_param(ir,3) * Z(is)*kbig_upar_nn_neo(is)
    enddo
    
    ! Toroidal component of Bootstrap current = sum <Z*n*utor/R>/<1/R>   
    jtor_nn_neo = jpar_nn_neo*Btor2_avg/Bmag2_avg
    do is=1, n_species
       jtor_nn_neo = jtor_nn_neo + rho(ir)*geo_param(ir,1) &
            * dens(is,ir)*temp(is,ir) &
            * (dlnndr(is,ir) + dlntdr(is,ir))*(1.0-Btor2_avg/Bmag2_avg)
    enddo
    jtor_nn_neo = jtor_nn_neo/(Btor_th0*bigR_th0*bigRinv_avg)
    
    ! Particle Flux
    data = trim(root)//'/../neural/neonn/pflux_nn/'
    ierr=load_anns(3, trim(data)//char(0),'brainfuse'//char(0))
    if(ierr == 0) then
       call neo_error('ERROR: (NEO) Neural network loading failed.')
       return
    endif
    ierr=load_anns_inputs(nn_in)
    ierr=run_anns()
    ierr=get_anns_avg_array(nn_out)
    if(ierr > 0) then
       call neo_error('ERROR: (NEO) Particle flux neural network failed.')
       return
    endif
    do k=1,6
       C_ke(k)   = nn_out(k)
       C_ki1(k)  = nn_out(k+6)
       C_ki2(k)  = nn_out(k+12)
    enddo
    pflux_nn_neo(:) = 0.0
    do k=1,6
       pflux_nn_neo(is_ele) = pflux_nn_neo(is_ele) &
            + (geo_param(ir,1) * rho(ir))**2 / geo_param(ir,3) &
            * nu(is_ele,ir) * C_ke(k)*C_ln(k)
       pflux_nn_neo(is_i1) = pflux_nn_neo(is_i1) &
            + (geo_param(ir,1) * rho(ir))**2 / geo_param(ir,3) &
            * nu(is_ele,ir) * C_ki1(k)*C_ln(k)
       pflux_nn_neo(is_i2) = pflux_nn_neo(is_i2) &
            + (geo_param(ir,1) * rho(ir))**2 / geo_param(ir,3) &
            * nu(is_ele,ir) * C_ki2(k)*C_ln(k)
    enddo

    ! Energy Flux
    data = trim(root)//'/../neural/neonn/eflux_nn/'
    ierr=load_anns(4, trim(data)//char(0),'brainfuse'//char(0))
    if(ierr == 0) then
       call neo_error('ERROR: (NEO) Neural network loading failed.')
       return
    endif
    ierr=load_anns_inputs(nn_in)
    ierr=run_anns()
    ierr=get_anns_avg_array(nn_out)
    if(ierr > 0) then
       call neo_error('ERROR: (NEO) Energy flux neural network failed.')
       return
    endif
    do k=1,6
       C_ke(k)   = nn_out(k)
       C_ki1(k)  = nn_out(k+6)
       C_ki2(k)  = nn_out(k+12)
    enddo
    eflux_nn_neo(:) = 0.0
    do k=1,6
       eflux_nn_neo(is_ele) = eflux_nn_neo(is_ele) &
            + (geo_param(ir,1) * rho(ir))**2 / geo_param(ir,3) &
            * nu(is_ele,ir) * C_ke(k)*C_ln(k)
       eflux_nn_neo(is_i1) = eflux_nn_neo(is_i1) &
            + (geo_param(ir,1) * rho(ir))**2 / geo_param(ir,3) &
            * nu(is_ele,ir) * C_ki1(k)*C_ln(k)
       eflux_nn_neo(is_i2) = eflux_nn_neo(is_i2) &
            + (geo_param(ir,1) * rho(ir))**2 / geo_param(ir,3) &
            * nu(is_ele,ir) * C_ki2(k)*C_ln(k)
    enddo
    
  end subroutine compute_nn_flow
  
  subroutine NEURAL_write(ir)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    integer :: is
    
    if(silent_flag == 0 .and. i_proc == 0) then
       open(io,file=trim(path)//runfile_nn,status='old',position='append')
       write(io,'(e16.8)',advance='no') r(ir)
       write(io,'(e16.8)',advance='no') jpar_nn_neo
       do is=1,n_species
          write(io,'(e16.8)',advance='no') pflux_nn_neo(is)
          write(io,'(e16.8)',advance='no') eflux_nn_neo(is)
          write(io,'(e16.8)',advance='no') kbig_upar_nn_neo(is)
          write(io,'(e16.8)',advance='no') vpol_th0_nn_neo(is)
       enddo
    endif
    
  end subroutine NEURAL_write

end module neo_neural
