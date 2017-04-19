subroutine tgyro_write_input

  use mpi
  use tgyro_globals

  implicit none

  integer :: i_ion,is
  character(len=7) :: ttext
  character(len=1) :: itag

  !----------------------------------------------------------------
  ! Trap miscellaneous errors
  !
  ! - Bogus matching radius
  !
  if (i_bc < 0) then
     error_flag = 1
     error_msg = 'ERROR: (TGYRO) Problem with matching radius.'
  endif
  !
  ! - Advanced iteration methods cannot do zero iterations:
  !
  if (tgyro_iteration_method >= 4 .and. tgyro_relax_iterations == 0) then
     error_flag = 1
     error_msg = 'ERROR: (TGYRO) TGYRO_ITERATION_METHOD > 4 requires TGYRO_RELAX_ITERATIONS > 0.'
  endif
  !
  ! - Need rotation physics to evolve Er
  !  
  if (tgyro_rotation_flag == 0 .and. loc_er_feedback_flag == 1) then
     error_flag = 1
     error_msg = 'ERROR: (TGYRO) TGYRO_ROTATION_FLAG must be 1 to evolve Er.'
  endif
  !
  ! - Need to set at least one profile to evolve, even for zero iterations, just to 
  !   define a residual error.
  !  
  if (loc_er_feedback_flag == 0 .and. &
       loc_te_feedback_flag == 0 .and. &
       loc_ti_feedback_flag == 0 .and. &
       maxval(evo_e) <= 0) then
     error_flag = 1
     error_msg = 'ERROR: (TGYRO) Must set one profile to evolve, even if running for 0 iterations.'
  endif
  if (maxval(evo_e(0:loc_n_ion)) == -1) then
     error_flag = 1
     error_msg = 'ERROR: (TGYRO) All species densities cannot be simultaneously floated.'
  endif


  !----------------------------------------------------------------

  if (i_proc_global == 0) then

     open(unit=1,file=trim(runfile),position='append')

     write(1,*) 'Job control'
     write(1,*) 

     !--------------------------------------------------------
     select case (loc_restart_flag)

     case (0)

        write(1,10) 'LOC_RESTART_FLAG','New simulation'

     case (1)

        write(1,10) 'LOC_RESTART_FLAG','Restarting old simulation'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_RESTART_FLAG'

     end select
     !--------------------------------------------------------

     write(1,*)
     write(1,*) 'Iteration control'
     write(1,*) 

     select case (tgyro_iteration_method)

     case (1)

        write(1,10) 'TGYRO_ITERATION_METHOD','Local Residual Minimization (original scheme)'

     case (4)

        write(1,10) 'TGYRO_ITERATION_METHOD','Block method (serial)'

     case (5)

        write(1,10) 'TGYRO_ITERATION_METHOD','Block method (parallel)'

     case (6)

        write(1,10) 'TGYRO_ITERATION_METHOD','Simple Relax (*IN DEVELOPMENT*)'

     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_ITERATION_METHOD'

     end select

     !--------------------------------------------------------
     select case (loc_residual_method)

     case (1)

        write(1,10) 'LOC_RESIDUAL_METHOD','(f-g)^2/(f^2+g^2)'

     case (2)

        write(1,10) 'LOC_RESIDUAL_METHOD','abs(f-g)'

     case (3)

        write(1,10) 'LOC_RESIDUAL_METHOD','(f-g)^2'

     case (4)

        write(1,10) 'LOC_RESIDUAL_METHOD','(f-g)^2/MAX(1,(f^2+g^2))'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_RESIDUAL_METHOD'

     end select
     !--------------------------------------------------------

     write(1,20) 'LOC_DX (Jacobian dx)',loc_dx
     if (maxval(flux_method_vec) == 4) then
        write(1,20) 'LOC_DX_GYRO (GYRO Jacobian dx)',loc_dx_gyro
     endif
     write(1,20) 'LOC_DX_MAX (maximum dx)',loc_dx_max
     write(1,20) 'LOC_RELAX (conv. relaxation)',loc_relax

     write(1,*)
     write(1,*) 'Scenario control'
     write(1,*) 

     !--------------------------------------------------------
     select case (TGYRO_EXPWD_FLAG)

     case (0)

        write(1,10) 'TGYRO_EXPWD_FLAG','Turbulent exchange OFF'

     case (1)

        write(1,10) 'TGYRO_EXPWD_FLAG','Turbulent exchange ON'

     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_EXPWD_FLAG'

     end select
     !--------------------------------------------------------


     !--------------------------------------------------------
     select case (loc_scenario)

     case (1)

        write(1,10) 'LOC_SCENARIO','Experimental, static exchange'

     case (2)

        write(1,10) 'LOC_SCENARIO','Experimental, dynamic exchange'

     case (3)

        write(1,10) 'LOC_SCENARIO','Reactor with input aux. power'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_SCENARIO'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     if (loc_scenario > 2) then

        select case (tgyro_dt_method)

        case (1)

           write(1,10) 'TGYRO_DT_METHOD','Reaction cross section <n1*n2> (use with separate D and T)'
           if (loc_n_ion == 1) then
              error_flag = 1
              error_msg = 'ERROR: (tgyro) Need LOC_N_ION > 1 for D-T reaction cross-section'
           endif

        case (2)

           write(1,10) 'TGYRO_DT_METHOD','Reaction cross section <n1*n1>/4 (use with single main ion)'

        case default

           error_flag = 1
           error_msg = 'Error: TGYRO_DT_METHOD'

        end select

     endif
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_neo_method)

     case (1)

        if (loc_chang_hinton == 1) then
           write(1,10) 'LOC_NEO_METHOD','Hinton-Hazeltine theory with Chang-Hinton Qi (ignores impurities)'
        else
           write(1,10) 'LOC_NEO_METHOD','Hinton-Hazeltine theory (ignores impurities)'
        endif

     case (2)

        write(1,10) 'LOC_NEO_METHOD','NEO code'
        if (loc_n_ion >= 4) then
           write(1,10) 'INFO: (tgyro) Using reduced energy resolution to cope with so many ions.'
        endif

     case (3)

        write(1,10) 'LOC_NEO_METHOD','Hirshman-Sigmar theory for all species.'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_NEO_METHOD'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (tgyro_tglf_revision)

     case(0)

        write(1,10) 'TGYRO_TGLF_REVISION', 'DEFAULTS + OVERWRITES'

     case (1)

        write(1,10) 'TGYRO_TGLF_REVISION','APS07'

     case (2)

        write(1,10) 'TGYRO_TGLF_REVISION','TGLF-09'

     case (3)

        write(1,10) 'TGYRO_TGLF_REVISION','in development'

     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_TGLF_REVISION'

     end select
     !--------------------------------------------------------

     select case(tgyro_glf23_revision)
     case (1)

        write(1,10) 'TGYRO_GLF23_REVISION','Original GLF23'

     case (2)

        write(1,10) 'TGYRO_GLF23_REVISION','retuned GLF23 v1.61'

     case(3)

        write(1,10) 'TGYRO_GLF23_REVISION','renormed GLF23'
     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_GLF23_REVISION'

     end select
     !--------------------------------------------------------

     !---------------------------------------------------------------------------------------------------

     write(1,*) 
     write(1,*) 'Profile evolution settings'
     write(1,*) 

     select case (loc_ti_feedback_flag)

     case (0)
        write(1,10) 'LOC_TI_FEEDBACK_FLAG','Ti evolution OFF'

     case (1)
        write(1,10) 'LOC_TI_FEEDBACK_FLAG','Ti evolution ON'

     case default
        error_flag = 1
        error_msg = 'Error: LOC_TI_FEEDBACK_FLAG'

     end select

     select case (loc_te_feedback_flag)

     case (0)
        write(1,10) 'LOC_TE_FEEDBACK_FLAG','Te evolution OFF'

     case (1)
        write(1,10) 'LOC_TE_FEEDBACK_FLAG','Te evolution ON'

     case default
        error_flag = 1
        error_msg = 'Error: LOC_TE_FEEDBACK_FLAG'

     end select

     select case (loc_er_feedback_flag)

     case (0)
        write(1,10) 'LOC_ER_FEEDBACK_FLAG','Er evolution OFF'

     case (1)
        write(1,10) 'LOC_ER_FEEDBACK_FLAG','Er evolution ON'

     case default
        error_flag = 1
        error_msg = 'Error: LOC_ER_FEEDBACK_FLAG'

     end select

     select case (evo_e(0))

     case (-1)
        write(1,10) 'TGYRO_DEN_METHOD0','ne profile set using quasineutrality'
     case (0)
        write(1,10) 'TGYRO_DEN_METHOD0','ne profile fixed'
     case (1)
        write(1,10) 'TGYRO_DEN_METHOD0','ne evolution ON'
     case default
        error_flag = 1
        error_msg = 'Error: TGYRO_EVO_E(0)'
     end select

     do is=1,loc_n_ion
        itag = trim(ion_tag(is))
        select case (evo_e(is))

        case (-1)
           write(1,10) 'TGYRO_DEN_METHOD'//itag,'ni'//itag//' profile set using quasineutrality'
        case (0)
           write(1,10) 'TGYRO_DEN_METHOD'//itag,'ni'//itag//' profile fixed'
        case (1)
           write(1,10) 'TGYRO_DEN_METHOD'//itag,'ni'//itag//' evolution ON'
        case default
           error_flag = 1
           error_msg = 'Error: TGYRO_DEN_METHOD'//itag//': invalid value.'
        end select
     enddo

     !---------------------------------------------------------------------------------------------------

     write(1,*)
     write(1,*) 'Physics parameters'
     write(1,*) 

     write(1,20) 'LOC_RMIN',r(2)/r_min
     write(1,20) 'LOC_RMAX',r(n_r)/r_min
     write(1,20) 'Pivot radius',r(i_bc)/r_min
     write(1,20) 'LOC_NU_SCALE',loc_nu_scale
     write(1,20) 'LOC_BETAE_SCALE',loc_betae_scale
     if (loc_betae_scale > 0.0) then
        write(1,10) '-> Fluctuations','Electromagnetic'
     else
        write(1,10) '-> Fluctuations','Electrostatic'
     endif
     if (tgyro_ptot_flag == 0) then
        write(1,10) '-> Pressure','Using pressure (beta) computed by summing over included species.'
     else
        write(1,10) '-> Pressure','Using total pressure from GFILE.'
     endif

     !--------------------------------------------------------
     select case (loc_zeff_flag)

     case (0)

        write(1,10) 'LOC_ZEFF_FLAG','Z_eff=1'

     case (1)

        write(1,10) 'LOC_ZEFF_FLAG','Z_eff from data'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_ZEFF_FLAG'

     end select
     !--------------------------------------------------------
     !--------------------------------------------------------
     select case (loc_num_equil_flag)

     case (0)

        write(1,10) 'LOC_NUM_EQUIL_FLAG','Using Miller model equilibrium in TGLF/GYRO/NEO.'

     case (1)

        write(1,10) 'LOC_NUM_EQUIL_FLAG','Using general Fourier equilibrium in TGLF/GYRO/NEO.'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_NUM_EQUIL_FLAG'

     end select
     !--------------------------------------------------------
     !--------------------------------------------------------
     select case (loc_lock_profile_flag)

     case (0)

        write(1,10) 'LOC_LOCK_PROFILE_FLAG','Initial profiles from gradients'

     case (1)

        write(1,10) 'LOC_LOCK_PROFILE_FLAG','Initial profiles from data'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_LOCK_PROFILE_FLAG'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_evolve_grad_only_flag)

     case (0)

        write(1,10) 'LOC_EVOLVE_GRAD_ONLY_FLAG','Self-consistently evolve profiles with gradients'

     case (1)

        write(1,10) 'LOC_EVOLVE_GRAD_ONLY_FLAG','Evolve gradients but *NOT* corresponding profiles'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_LOCK_PROFILE_FLAG'

     end select

     !--------------------------------------------------------

     write(1,*)
     write(1,*) 'Pedestal parameters'
     write(1,*) 

     select case (tgyro_ped_model)

     case (1)

        write(1,10) 'TGYRO_PED_MODEL','No pedestal.  Use fixed pivot.'

     case (2)

        write(1,10) 'TGYRO_PED_MODEL','Dynamic EPED1_NN pedestal.'

     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_PED_MODEL'

     end select

     write(1,*)
     write(1,*) 'Ion parameters'
     write(1,*) 

     do i_ion=1,loc_n_ion
        if (therm_flag(i_ion) == 1) then
           ttext = 'thermal'
           if (i_ion == i_ash) ttext = 'ash'
        else
           ttext = 'fast'
        endif
        if (i_ion == 1) then
           write(1,40) 'ion 1 [mass,charge,type]',mi_vec(i_ion),zi_vec(i_ion),ttext
        else
           write(1,40) 'ion '//trim(ion_tag(i_ion))//' [mass,charge,type]',mi_vec(i_ion),zi_vec(i_ion),ttext
        endif
     enddo
     write(1,*) 
     write(1,10) 'INFO: (tgyro)','GyroBohm factors defined by ion 1 mass.'
     write(1,*)
     !--------------------------------------------------------

     !--------------------------------------------------------
     write(1,*)
     write(1,*) 'Rotation and field orientation'
     write(1,*) 

     if (tgyro_rotation_flag == 1) then
        write(1,10) 'TGYRO_ROTATION_FLAG','Rotation effects ON'
     else
        write(1,10) 'TGYRO_ROTATION_FLAG','Rotation effects OFF'
     endif
     write(1,30) 'IPCCW',-signb*signq
     write(1,30) 'BTCCW',-signb
     !--------------------------------------------------------

     !--------------------------------------------------------
     write(1,*)
     write(1,*) 'Input profile rescaling factors'
     write(1,*) 
     write(1,20) 'TGYRO_INPUT_DEN_SCALE',tgyro_input_den_scale
     write(1,20) 'TGYRO_INPUT_TE_SCALE',tgyro_input_te_scale
     write(1,20) 'TGYRO_INPUT_TI_SCALE',tgyro_input_ti_scale
     write(1,20) 'TGYRO_INPUT_W0_SCALE',tgyro_input_w0_scale
     write(1,20) 'TGYRO_INPUT_PAUX_SCALE',tgyro_input_paux_scale

     write(1,*)'-----------------------------------------------------------------'

     close(1)

  endif

  !--------------------------------------------------------------
  ! Trap GYRO input data error and stop 
  call MPI_BCAST(error_flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(error_msg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !
  if (error_flag == 1) then
     call tgyro_catch_error(trim(error_msg))
  endif
  !--------------------------------------------------------------

10 format(t2,a,t33,':',t35,a)
20 format(t2,a,t33,':',t35,f8.4,1x,a)
30 format(t2,a,t33,':',t35,i4)
40 format(t2,a,t33,':',t35,f6.2,1x,f6.2,3x,a)

end subroutine tgyro_write_input

