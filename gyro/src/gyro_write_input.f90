!------------------------------------------------------
! gyro_write_stdout.f90 [caller BigScience]
!
! PURPOSE:
!  This subroutine summarizes most input parameters
!  and derived quantities of interest.
!------------------------------------------------------

subroutine gyro_write_input

  use mpi
  use gyro_globals
  use gyro_profile_exp
  use gyro_pointers
  use math_constants
  use geo
  use ompdata

  !--------------------------------------------------
  implicit none
  !
  integer :: i_ion
  !
  real :: co_r_loc,co_r
  real :: co_tau
  real :: neutral
  real :: c_min
  real :: c_max
  real :: num_rho
  real :: num_rho_l
  real :: num_rho_r
  real :: num_rho_min
  real :: x_s
  real :: den_total
  real :: v_alfven
  real :: rhos_abs
  !
  real, dimension(n_ion) :: band
  character (len=16), dimension(n_spec) :: n_tag
  character (len=16), dimension(n_spec) :: dn_tag
  character (len=16), dimension(n_spec) :: t_tag
  character (len=16), dimension(n_spec) :: dt_tag
  character (len=16), dimension(n_spec) :: mu_tag
  character (len=20), dimension(n_spec) :: emax_tag
  character (len=16), dimension(n_spec) :: z_tag
  character (len=16), dimension(n_spec) :: nu_tag
  character (len=2), dimension(n_ion) :: e_tag
  !--------------------------------------------------

  rhos_abs = abs(rhos_norm)

  !--------------------------------------------------
  ! Define some character tags for printing of local 
  ! parameters:
  !
  e_tag(1) = ''
  do i_ion=2,n_ion
     e_tag(i_ion) = '_'//achar(i_ion-1+iachar("1"))  
  enddo

  n_tag(1) = '.NI_OVER_NE'
  do i_ion=2,n_ion
     n_tag(i_ion) = '.NI_OVER_NE'//e_tag(i_ion)
  enddo
  n_tag(n_spec) = '(null)'

  dn_tag(1) = '.DLNNDR'
  do i_ion=2,n_ion
     dn_tag(i_ion) = '.DLNNDR'//e_tag(i_ion)
  enddo
  dn_tag(n_spec) = '.DLNNDR_ELECTRON'

  t_tag(1) = '.TI_OVER_TE'
  do i_ion=2,n_ion
     t_tag(i_ion) = '.TI_OVER_TE'//e_tag(i_ion)
  enddo
  t_tag(n_spec) = '(null)'

  dt_tag(1) = '.DLNTDR'
  do i_ion=2,n_ion
     dt_tag(i_ion) = '.DLNTDR'//e_tag(i_ion)
  enddo
  dt_tag(n_spec) = '.DLNTDR_ELECTRON'

  mu_tag(1) = '.MU'
  do i_ion=2,n_ion
     mu_tag(i_ion) = '.MU'//e_tag(i_ion)
  enddo
  mu_tag(n_spec) = '.MU_ELECTRON'

  emax_tag(1) = '.ENERGY_MAX'
  do i_ion=2,n_ion
     emax_tag(i_ion) = '.ENERGY_MAX'//e_tag(i_ion)
  enddo
  emax_tag(indx_e) = '.ENERGY_MAX_ELECTRON'

  z_tag(1) = '.Z'
  do i_ion=2,n_ion
     z_tag(i_ion) = '.Z'//e_tag(i_ion)
  enddo
  z_tag(n_spec) = '(Z_ELECTRON)'

  nu_tag(1) = '(NU_II)'
  do i_ion=2,n_ion
     nu_tag(i_ion) = '(NU_II'//e_tag(i_ion)//')'
  enddo
  nu_tag(n_spec) = '(NU_EI)'
  !--------------------------------------------------

  if (i_proc == 0 .and. output_flag == 1 .and. gkeigen_j_set == 0) then
     open(unit=1,file=trim(runfile),status='old',position='append')

     write(1,*) '---------------- PARLLELISM DIMENSIONS --------------------'
     write(1,10) 'MPI tasks',n_proc
     write(1,10) 'OpenMP threads',n_omp
     write(1,*) '----------------- GRID DIMENSIONS -------------------------'
     write(1,10) 'n_n',n_n
     write(1,10) 'n_x',n_x
     write(1,10) 'n_stack',n_stack
     write(1,10) 'n_blend',n_blend
     write(1,10) 'n_pass',n_pass
     write(1,10) 'n_trap',n_trap
     write(1,10) 'n_energy',n_energy
     write(1,10) 'n_field',n_field
     write(1,*) '--'
     write(1,10) 'n_spec',n_spec
     write(1,10) 'n_ion',n_ion
     write(1,10) 'n_kinetic',n_kinetic
     write(1,10) 'n_gk',n_gk
     write(1,10) 'indx_e',indx_e
     write(1,*) '--'
     write(1,30) 'FIELD POINTS:',n_n*n_x*n_blend*n_field
     write(1,30) 'DIST. POINTS:',n_n*n_x*n_stack*(n_pass+n_trap)*n_energy*n_kinetic
     write(1,*) '--'
     write(1,20) 'energy_max',energy_max
     write(1,20) 'dt',dt

     write(1,*) '--------------------- LOCAL PARAMETERS --------------------'
     write(1,*) ' NOTE: use abs(SAFETY_FACTOR) as input'
     write(1,20) '# RADIUS [INPUT]',r0
     write(1,20) '.RADIUS',r_norm
     write(1,20) '.ASPECT_RATIO',rmaj_s(ir_norm)
     write(1,20) '.SHIFT',drmaj_s(ir_norm)
     write(1,20) '.ZMAG',zmag_s(ir_norm)
     write(1,20) '.DZMAG',dzmag_s(ir_norm)
     write(1,20) '.KAPPA',kappa_s(ir_norm)
     write(1,20) '.S_KAPPA',s_kappa_s(ir_norm)
     write(1,20) '.DELTA',delta_s(ir_norm)
     write(1,20) '.S_DELTA',s_delta_s(ir_norm)
     write(1,20) '.ZETA',zeta_s(ir_norm)
     write(1,20) '.S_ZETA',s_zeta_s(ir_norm)
     write(1,20) '.SAFETY_FACTOR',abs(q_norm)
     write(1,20) '.SHEAR',shat_norm 
     write(1,20) '.RHO_STAR',rhos_norm
     write(1,20) '.Z_EFF',z_eff_s(ir_norm)
     write(1,20) '.MACH',mach_s(ir_norm)
     write(1,20) '.PGAMMA',gamma_p_s(ir_norm)
     write(1,20) '.GAMMA_E',gamma_e_s(ir_norm)
!     write(1,20) '.LAMBDA_DEBYE',lambda_debye
     write(1,20) '.NU_EI',nu_s(indx_e,ir_norm)
     write(1,20) '.IPCCW',ipccw
     write(1,20) '.BTCCW',btccw
     i = ir_norm
     do is=1,n_spec-1
        if (electron_method == 3) then
           i_ion = n_spec-is+1
           write(1,'(a)') ' # Adiabatic Ion '
        else
           i_ion = is
           write(1,'(a,i2)') ' # Ion ',i_ion
        endif
        write(1,20) n_tag(is),den_s(i_ion,i)
        write(1,20) t_tag(is),tem_s(i_ion,i)
        write(1,20) dn_tag(is),dlnndr_s(i_ion,i)
        write(1,20) dt_tag(is),dlntdr_s(i_ion,i)
        write(1,20) z_tag(is),z(i_ion)
        write(1,20) nu_tag(is),nu_s(i_ion,ir_norm)
        write(1,20) mu_tag(is),mu(i_ion) 
     enddo

     if (electron_method > 1) then
        write(1,*) '# Electrons'
        write(1,20) dn_tag(n_spec),dlnndr_s(indx_e,ir_norm)
        write(1,20) dt_tag(n_spec),dlntdr_s(indx_e,ir_norm)
        write(1,20) z_tag(n_spec),z(indx_e)
        write(1,20) '.BETAE_UNIT',betae_unit_norm
        write(1,20) mu_tag(n_spec),mu(indx_e)
     endif

     !--------------------------------------------
     ! Check for charge balance:
     !
     neutral = 0.0
     do i=1,n_x
        neutral = neutral+sum(z(:)*den_s(:,i))/n_x
     enddo
     !--------------------------------------------

     write(1,*) '------------- LOCAL PARAMETERS (diagnostic) ---------------'
     write(1,20) 'n_i*z_i - n_e: ',neutral
     write(1,20) 'r/R0',r(ir_norm)/rmaj_s(ir_norm)
     if (radial_profile_method == 3) write(1,20) 'rho_norm',rhogrid_s(ir_norm)
     write(1,20) 'b_unit',b_unit_norm
     write(1,20) 'beta_unit_norm',beta_unit_s(ir_norm)
     i = ir_norm
     if (electron_method /= 3) then
        do is=1,n_spec-1
           write(1,20) 'betai_unit'//e_tag(is),beta_unit_s(i)* &
                pr_s(is,i)/sum(pr_s(1:n_spec,i))
        enddo
     else
        do is=1,n_spec-1
           write(1,20) 'betai_unit'//e_tag(is),beta_unit_s(i)* &
                pr_s(n_spec-is+1,i)/sum(pr_s(1:n_spec,i))
        enddo
     endif
     write(1,20) 'betae_unit_norm',betae_unit_norm
     write(1,20) 'beta_*',beta_star_s(ir_norm)
     write(1,20) 'alpha_MHD',q_norm**2*rmaj_s(ir_norm)*beta_star_s(ir_norm)
     write(1,20) 'lambda_*',lambda_debye/rhos_norm
     write(1,20) 'volume (1/a^3)',GEO_volume
     write(1,20) 'omega_00 (c_s/a)',w0_s(ir_norm)
     write(1,*) '* Note that f = f_sim exp(i n omega0[r0] t)'
     if (betae_unit_norm > 0.0) then
        den_total = sum(den_s(1:n_ion,ir_norm)/den_s(indx_e,ir_norm)/mu(1:n_ion)**2)
        v_alfven  = sqrt(2.0/(betae_unit_norm*den_total))

        write(1,*) '----------- ALFVEN WAVE PARAMETERS (diagnostic) -----------'
        write(1,20) '(v_A/c_s)',v_alfven
        write(1,20) 'Omega_TAE',v_alfven/(2*q_norm*rmaj_s(ir_norm))
        write(1,20) 'Omega_A',v_alfven/(q_norm*rmaj_s(ir_norm))
     endif
     write(1,*) '------------------- UPWIND PARAMETERS ---------------------'
     write(1,20) 'radial_upwind',radial_upwind
     do i_ion=1,n_spec-1
        write(1,20) 'orbit_upwind'//e_tag(i_ion),orbit_upwind_vec(i_ion)
     enddo
     if (electron_method == 2 .or. electron_method == 4) then
        write(1,20) 'orbit_upwind (elec)',orbit_upwind_vec(0)
     endif
     if (source_method > 1) then
        write(1,*) '----------- SOURCE PARAMETERS -----------------'
        write(1,20) 'nu_source',nu_source
        write(1,10) 'n_source',n_source
     endif
     write(1,*) '--------------- RADIAL DOMAIN PARAMETERS ------------------'
     write(1,20) 'box_multiplier',box_multiplier
     write(1,20) 'L/a',x_length
     write(1,*) '--'
     write(1,40) 'explicit_damp',n_explicit_damp,explicit_damp
     write(1,*) '--'
     write(1,20) 'LEFT : r_a',r(1)
     write(1,20) 'LEFT : r_a_physical',r(1+n_explicit_damp) 
     write(1,20) 'NORM : r(ir_norm)',r(ir_norm)
     write(1,20) 'RIGHT: r_b_physical',r(n_x-n_explicit_damp)
     write(1,20) 'RIGHT: r_b',r(n_x)

  endif

  if (gyrotest_flag == 0 .and. output_flag == 1) then

     co_tau = maxval(abs(v_theta(:,:,:,:)))*dt/d_tau(1)*&
          maxval(mu(1:n_kinetic))
     co_r_loc = maxval(abs(omega_dr(:,:,:,:)))*dt/d_x

     call MPI_ALLREDUCE(co_r_loc,&
          co_r,&
          1,&
          MPI_DOUBLE_PRECISION,&
          MPI_MAX,&
          GYRO_COMM_WORLD,&
          i_err)

     if (collision_flag == 1) then
        call MPI_ALLREDUCE(condition_number,&
             c_max,&
             1,&
             MPI_DOUBLE_PRECISION,&
             MPI_MAX,&
             GYRO_COMM_WORLD,&
             i_err)
     endif

     if (i_proc == 0) then
        write(1,*) '----------- PARALLELIZATION PARAMETERS --------'
        write(1,30) '(nek) per subgroup:',n_nek_loc_1
        write(1,30) '(ine) per subgroup:',n_ine_loc_1
        write(1,*) '----------- TIME STEPPING PARAMETERS ----------'
        write(1,50) 'time_skip',time_skip
        write(1,50) 'restart_data_skip',restart_data_skip
        write(1,*) '----------- STABILITY PARAMETERS ----------------'
        write(1,20) 'd/dtau Courant',co_tau
        write(1,20) '  d/dr Courant',co_r 
        if (collision_flag == 1) then
           write(1,20) 'Log(RBF Cond. num.)',log10(c_max)
        endif
     endif

  endif

  if (i_proc == 0 .and. output_flag == 1) then
     write(1,*) '------------- CENTRAL WAVENUMBERS SIMULATED ---------------'
     write(1,*) '    (k_y = nq/r, rho = rho_sD_unit)'
  endif

  c_min = q_s(1)/r_s(1)*r_s(ir_norm)/q_s(ir_norm)*&
       sqrt(tem_s(n_spec,1)/tem_s(n_spec,ir_norm))*b_unit_s(ir_norm)/b_unit_s(1)

  c_max = q_s(n_x)/r_s(n_x)*r_s(ir_norm)/q_s(ir_norm)*&
       sqrt(tem_s(n_spec,n_x)/tem_s(n_spec,ir_norm))*b_unit_s(ir_norm)/b_unit_s(n_x)

  if (gyrotest_flag == 0) then

     call collect_real(krho_i(in_1,ir_norm),krho_collect)

     if (output_flag == 1) then

        do in=1,n_n

           if (i_proc == 0) then

              if (krho_collect(in)*c_min >= 10.0 .or. &
                   krho_collect(in)*c_max >= 10.0) then
                 write(1,'(t2,a,i5,3x,a,f6.2,1x,2(a,f6.2,a))') &
                      'n = ',n(in), 'k_y rho = ', &
                      krho_collect(in),&
                      '[',krho_collect(in)*c_min,']',&
                      '[',krho_collect(in)*c_max,']'
              else
                 write(1,'(t2,a,i5,3x,a,f6.3,1x,2(a,f6.3,a))') &
                      'n = ',n(in), 'k_y rho = ', &
                      krho_collect(in),&
                      '[',krho_collect(in)*c_min,']',&
                      '[',krho_collect(in)*c_max,']'
              endif

           endif

        enddo ! in

     endif

  else

     i = ir_norm
     krho_collect(:) = n(:)*q_s(i)/r_s(i)*rhos_norm/b_unit_s(i)

     if (output_flag == 1) then

        do in=1,n_n
           if (krho_collect(in)*c_min >= 10.0 .or. &
                krho_collect(in)*c_max >= 10.0) then
              write(1,'(t2,a,i5,3x,a,f6.2,1x,2(a,f6.2,a))') &
                   'n = ',n(in), 'k_y rho = ', &
                   krho_collect(in),&
                   '[',krho_collect(in)*c_min,']',&
                   '[',krho_collect(in)*c_max,']'
           else
              write(1,'(t2,a,i5,3x,a,f6.3,1x,2(a,f6.3,a))') &
                   'n = ',n(in), 'k_y rho = ', &
                   krho_collect(in),&
                   '[',krho_collect(in)*c_min,']',&
                   '[',krho_collect(in)*c_max,']'
           endif
        enddo ! in

     endif

  endif


  if (i_proc == 0 .and. output_flag == 1) then
     write(1,*) ' '
     write(1,'(t2,a,t27,f10.6)') 'min resolved k_x*rho_sD : ', &
          2*pi/(x_length/rhos_abs)
     write(1,'(t2,a,t27,f10.6)') 'max resolved k_x*rho_sD : ', &
          2*pi/(x_length/rhos_abs)*n_x/4.0

     !-----------------------------------------------------------
     ! Compute central box size
     !
     write(1,*) '--------------- CENTRAL BOX SIZE -----------------'
     write(1,25) 'abs(Lx/rho_sD)',x_length/abs(rhos_abs)
     if (krho_collect(1) == 0.0) then
        if (n_n > 1) then
           write(1,25) 'abs(Ly/rho_sD)',abs(2*pi/krho_collect(2))
        endif
     else
        write(1,25) 'abs(Ly/rho_sD)',abs(2*pi/krho_collect(1))
     endif
     write(1,25) '=> abs(dx/rho_sD)',x_length/rhos_abs/n_x
     write(1,*) '--------------------------------------------------'
     !
     ! Check for sufficient gyro bandwidth
     !
     write(1,*) '------------ GYRO-OPERATOR BANDWIDTH -------------'
     if (m_gyro == n_x/2) then
        write(1,*) '=> EXACT BESSEL.'
     else
        i = ir_norm
        do i_ion=1,n_ion
           x_s = sqrt(tem_s(i_ion,i))/mu(i_ion)/abs(z(i_ion))
           num_rho = x_length/(rhos_abs*x_s)
           band(i_ion) = (2*m_gyro+1.0)*num_rho/n_x
           write(1,25) 'Ion '//achar(i_ion-1+iachar("1")),band(i_ion)
        enddo
        if (minval(band) < 10.0 .and. nonlinear_flag == 1) then
           write(1,*) 'WARNING: RADIAL_GYRO_BAND should be > 10.0'
        endif
     endif

     write(1,*) '--------------------------------------------------'
     !
     ! Check for sufficient boundary layer width
     !
     if (boundary_method == 2) then
        write(1,*) '---------------- BUFFER WIDTH --------------------'
        num_rho_min = 8.0
        do i_ion=1,n_ion
           i = n_explicit_damp/2+1
           x_s = sqrt(tem_s(i_ion,i))/mu(i_ion)/abs(z(i_ion))
           num_rho_l = (r(1+n_explicit_damp)-r(1))/(rhos_abs*x_s)
           i = n_x-n_explicit_damp/2
           x_s = sqrt(tem_s(i_ion,i))/mu(i_ion)/abs(z(i_ion))
           num_rho_r = (r(n_x)-r(n_x-n_explicit_damp))/(rhos_abs*x_s)
           write(1,25) 'Ion '//achar(i_ion-1+iachar("1")),num_rho_l,num_rho_r
           if (num_rho_l < num_rho_min) num_rho_min = num_rho_l
           if (num_rho_r < num_rho_min) num_rho_min = num_rho_r
        enddo
        if (num_rho_min < 8.0) then
           write(1,*) 'WARNING: All buffer widths should be > 8.0'
        endif
        write(1,*) '--------------------------------------------------'
     endif

     write(1,20) 'Startup time',startup_time
     write(1,*) '--------------------------------------------------'
     ! File list
     write(1,*) 'PLEASE SEE: '
     write(1,*) ' - out.gyro.units for normalizing parameters'
     write(1,*) ' - out.gyro.memory for memory usage'
     write(1,*) ' - out.gyro.efficiency for parallelization efficiency'
     write(1,*) ' - out.gyro.phase_space for velocity-space nodes and weights'
     close(1)

  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_write_input done]'
  endif

10 format(t2,a,t23,': ',i4) 
20 format(t2,a,t23,': ',20(f10.6,1x)) 
25 format(t2,a,t23,': ',20(f12.6,1x)) 
30 format(t2,a,t23,': ',i9) 
40 format(t2,a,t23,': ',i3,' (',f8.4,')') 
50 format(t2,a,t23,': ',i5) 

end subroutine gyro_write_input

