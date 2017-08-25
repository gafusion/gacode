module neo_neural

  implicit none

  public :: NEURAL_alloc, NEURAL_do, NEURAL_write
  real :: jpar_nn_neo, jtor_nn_neo, jpar_nn_sau, jtor_nn_sau
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

       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io,file=trim(path)//runfile_nn,status='replace')
          close(io)
       endif

       initialized = .true.

    else
       if(.NOT. initialized) return

       initialized = .false.

    endif

  end subroutine NEURAL_alloc
       
  subroutine NEURAL_do(ir)
    implicit none
    integer, intent (in) :: ir
    
    call compute_nn_jpar(ir)

  end subroutine NEURAL_do

  ! The NN for the bootstrap current presently assumes 2 ion species
  ! (D + C) and kinetic electrons.  Strong rotation effects are not included.
  subroutine compute_nn_jpar(ir)
    use mpi
    use neo_globals
    use neo_equilibrium
    implicit none
    
    integer, intent (in) :: ir
    integer :: ierr
    integer :: is
    integer :: is_ele, is_i1, is_i2
    real    :: d_max
    character(len=1000) :: jbsnn_model
    real(4), dimension(5)  :: nn_in
    real(4), dimension(12) :: nn_out
    real :: CTi2_neo, CTi1_neo, CTe_neo, CNi2_neo, CNi1_neo, CNe_neo
    real :: CTi2_sau, CTi1_sau, CTe_sau, CNi2_sau, CNi1_sau, CNe_sau
    
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
    
    ! Identify the electron species
    do is=1, n_species
       if(Z(is) < 0.0) then
          is_ele = is
          exit
       endif
    enddo
    
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
    nn_in(4) = r(ir)/rmaj(ir)                      ! r/R
    nn_in(3) = abs(q(ir))                          ! q
    nn_in(1) = dens(is_i1,ir)/dens(is_ele,ir)      ! ni1/ne
    nn_in(5) = temp(is_i1,ir)/temp(is_ele,ir)      ! Ti1/Te
    ! log(nuee/cs/R)
    nn_in(2) = log10(nu(is_ele,ir)*rmaj(ir)/sqrt(temp(is_ele,ir)))  
    
    ! Run the NN
    call get_environment_variable('JBSNN_MODEL_DIR',jbsnn_model)
    ierr=load_anns(0, trim(jbsnn_model)//char(0),'brainfuse'//char(0))
    ierr=load_anns_inputs(nn_in)
    ierr=run_anns()
    ierr=get_anns_avg_array(nn_out)
    
    ! Get coeffcients computed by the NN
    CTi2_neo = nn_out(1)
    CTi1_neo = nn_out(2)
    CTe_neo  = nn_out(3)
    CNi2_neo = nn_out(4)
    CNi1_neo = nn_out(5)
    CNe_neo  = nn_out(6)
    CTi2_sau = nn_out(7)
    CTi1_sau = nn_out(8)
    CTe_sau  = nn_out(9)
    Cni2_sau = nn_out(10)
    Cni1_sau = nn_out(11)
    Cne_sau  = nn_out(12)
    
    ! Reconstruct jpar from NEO NN and Sauter NN
    
    jpar_nn_neo = I_div_psip * rho(ir) * temp(is_ele,ir) &
         * (abs(Z(is_ele))*dens(is_ele,ir) &
         * (CTe_neo*dlntdr(is_ele,ir) + Cne_neo*dlnndr(is_ele,ir)) &
         + abs(Z(is_i1))*dens(is_i1,ir) &
         * (CTi1_neo*dlntdr(is_i1,ir) + Cni1_neo*dlnndr(is_i1,ir)) &
         + abs(Z(is_i2))*dens(is_i2,ir) &
         * (CTi2_neo*dlntdr(is_i2,ir) + Cni2_neo*dlnndr(is_i2,ir)))
    
    jpar_nn_sau = I_div_psip * rho(ir) * temp(is_ele,ir) &
         * (abs(Z(is_ele))*dens(is_ele,ir) &
         * (CTe_sau*dlntdr(is_ele,ir) + Cne_sau*dlnndr(is_ele,ir)) &
         + abs(Z(is_i1))*dens(is_i1,ir) &
         * (CTi1_sau*dlntdr(is_i1,ir) + Cni1_sau*dlnndr(is_i1,ir)) &
         + abs(Z(is_i2))*dens(is_i2,ir) &
         * (CTi2_sau*dlntdr(is_i2,ir) + Cni2_sau*dlnndr(is_i2,ir)))
    
    ! Toroidal component of Bootstrap current = sum <Z*n*utor/R>/<1/R>
    
    jtor_nn_neo = jpar_nn_neo*Btor2_avg/Bmag2_avg
    jtor_nn_sau = jpar_nn_sau*Btor2_avg/Bmag2_avg
    do is=1, n_species
       jtor_nn_neo = jtor_nn_neo + rho(ir)*I_div_psip*dens(is,ir)*temp(is,ir) &
            * (dlnndr(is,ir) + dlntdr(is,ir))*(1.0-Btor2_avg/Bmag2_avg)
       jtor_nn_sau = jtor_nn_sau + rho(ir)*I_div_psip*dens(is,ir)*temp(is,ir) &
            * (dlnndr(is,ir) + dlntdr(is,ir))*(1.0-Btor2_avg/Bmag2_avg)
    enddo
    jtor_nn_neo = jtor_nn_neo/(Btor_th0*bigR_th0*bigRinv_avg)
    jtor_nn_sau = jtor_nn_sau/(Btor_th0*bigR_th0*bigRinv_avg)
    
  end subroutine compute_nn_jpar
  
  subroutine NEURAL_write(ir)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    
    if(silent_flag == 0 .and. i_proc == 0) then
       open(io,file=trim(path)//runfile_nn,status='old',position='append')
       write(io,'(e16.8)',advance='no') r(ir)
       write(io,'(e16.8)',advance='no') jpar_nn_neo
       write(io,'(e16.8)',advance='no') jpar_nn_sau
    endif
    
  end subroutine NEURAL_write

end module neo_neural
