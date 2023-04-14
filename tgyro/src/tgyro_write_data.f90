!-----------------------------------------------------------
! tgyro_write_data.f90
!
! PURPOSE:
!  Manage I/O for everything but iteration diagnostics.
!
! NOTES:
!  Files are opened and closed at each iteration.
!----------------------------------------------------------

subroutine tgyro_write_data(i_print)

  use mpi
  use tgyro_globals
  use tgyro_ped
  use tgyro_iteration_variables, only : i_tran_loop
  use expro

  implicit none

  integer :: i,i_exp
  integer :: ip
  integer :: i_print
  integer :: is,i_ion
  integer :: p
  real, dimension(2:n_r,loc_n_ion+4) :: res2,relax2
  character(len=2) :: itag
  character(len=6) :: ntag,ttag
  character(len=50) :: date_str,time_str
  logical :: converged
  real, dimension(:,:), allocatable :: dum2d

  ! Renormalize residuals so the error estimates are comparable

  select case (loc_residual_method) 

  case (2)

     ! ABSOLUTE VALUE NORM
     res_norm = res

  case (3)

     ! SQUARE RESIDUAL
     res_norm = sqrt(res)

  end select

  ! Convergence status
  converged = (sum(res_norm)/size(res_norm) < tgyro_residual_tol)
  
  !====================================================
  ! input.gacode
  !====================================================

  if (tgyro_write_profiles_flag /= 0 .and. i_print > 0) then 

     call expro_read('input.gacode',MPI_COMM_WORLD)

     call tgyro_profile_reintegrate

     ! First, hose the original profiles 
     expro_ni = 0.0
     expro_ti = 0.0

     expro_ptot = ptot_exp
     expro_ne   = exp_ne*1e-13
     expro_te   = exp_te*1e-3

     expro_ni(1:loc_n_ion,:) = exp_ni(1:loc_n_ion,:)*1e-13
     expro_ti(1:loc_n_ion,:) = exp_ti(1:loc_n_ion,:)*1e-3
     expro_w0   = exp_w0
     expro_ptot = ptot_exp ! already in Pa
     do i_exp=1,expro_n_exp
        expro_z_eff(i_exp) = sum(exp_ni(1:loc_n_ion,i_exp)*zi_vec(1:loc_n_ion)**2)/&
             exp_ne(i_exp)
     enddo

     call tgyro_mpi_info('rad_alpha')

     ! Alpha power density
     call rad_alpha(exp_ne,&
          exp_ni(1:loc_n_ion,:),& ! 1/cm^3
          exp_te,               & ! eV
          exp_ti(1:loc_n_ion,:),& ! eV
          expro_qei,&             ! dummy arg
          expro_qfusi,&           ! used (out) erg/cm^3/s
          expro_qfuse,&           ! used (out) erg/cm^3/s
          expro_qbrem,&           ! dummy arg
          expro_qsync,&           ! dummy arg
          expro_qei,&             ! dummy arg
          expro_qline,&           ! dummy arg
          expro_n_exp,&
          loc_n_ion)

     call tgyro_mpi_info('out')

     allocate(dum2d(loc_n_ion,expro_n_exp))
     
     ! Collisional exchange power density
     call collision_rates(&
          exp_ne,&       ! 1/cm^3
          exp_ni,&       ! 1/cm^3
          exp_te,&       ! eV
          exp_ti,&       ! eV
          dum2d,&        ! dummy arg (nui)
          expro_qsync,&  ! dummy arg (nue)
          exp_nu_exch, & ! 1/s
          expro_qline, & ! dummy arg (taus)
          expro_n_exp,&
          loc_n_ion)

     deallocate(dum2d)

     expro_qei = 1.5*exp_nu_exch(:)*exp_ne(:)*k*(exp_te(:)-exp_ti(1,:))

     ! Synchrotron radiation
     call rad_sync(aspect_rat,r_min,1e4*expro_bt0,exp_ne,exp_te,expro_qsync,expro_n_exp)

     ! Bremsstrahlung and line radiation (Post 1977) 
     call rad_ion_adas(&
          exp_te,&                ! eV
          exp_ne,&                ! 1/cm^3
          exp_ni(1:loc_n_ion,:),& ! 1/cm^3   
          zi_vec(1:loc_n_ion),&
          ion_name,&
          expro_qbrem,&
          expro_qline,&
          loc_n_ion,&
          expro_n_exp)

     ! Convert erg/cm^3/s -> W/cm^3 = MW/m^3
     expro_qbrem = expro_qbrem*1e-7
     expro_qsync = expro_qsync*1e-7
     expro_qline = expro_qline*1e-7
     expro_qfuse = expro_qfuse*1e-7
     expro_qfusi = expro_qfusi*1e-7
     expro_qei   = expro_qei*1e-7

     if (i_proc_global == 0) then

        call date_and_time(DATE=date_str,TIME=time_str)
        expro_head_tgyro = '#     tgyro : '//trim(date_str)//' '//trim(time_str)

        if (tgyro_write_profiles_flag == -1) then
           ! Output for each iteration
           write(ntag,'(i0)') i_tran
           call expro_write('input.gacode.'//trim(ntag))
        endif

        if (i_tran_loop == tgyro_relax_iterations .or. converged) then
           ! Output for last iteration
           call expro_write('input.gacode.new')
        endif

     endif

  endif

  if (i_proc_global == 0) then

     !------------------------------------------------------------------------------------------
     ! Initialization:
     !
     if (i_print == 0) then

        open(unit=1,file='out.tgyro.geometry.1',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.geometry.2',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.nu_rho',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.power_e',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.power_i',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.alpha',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.gyrobohm',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.residual',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.evo_ti',status='replace')
        close(1)
        open(unit=1,file='out.tgyro.evo_te',status='replace')
        close(1)
        open(unit=1,file='out.tgyro.evo_er',status='replace')
        close(1)
        open(unit=1,file='out.tgyro.evo_ne',status='replace')
        close(1)
        do i_ion=1,loc_n_ion
           open(unit=1,file='out.tgyro.evo_n'//trim(ion_tag(i_ion)),status='replace')
           close(1)
        enddo

        open(unit=1,file='out.tgyro.flux_e',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.profile',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.profile_e',status='replace')
        close(1)

        do i_ion=1,loc_n_ion
           open(unit=1,file='out.tgyro.flux_i'//trim(ion_tag(i_ion)),status='replace')
           close(1)
           open(unit=1,file='out.tgyro.profile_i'//trim(ion_tag(i_ion)),status='replace')
           close(1)
        enddo

        if (tgyro_ped_model > 1) then
           open(unit=1,file='out.tgyro.ped',status='replace')
           close(1)
        endif

        open(unit=1,file='out.tgyro.control',status='replace')
        write(1,*) n_r
        write(1,*) n_evolve
        write(1,*) 0
        close(1)

        ! Special case: iteration status file (see loc_write_intermediate)
        open(unit=1,file='out.tgyro.iterate',status='replace')
        close(1)

        open(unit=1,file='out.tgyro.prec',status='replace')
        close(1)

        return

     endif
     !------------------------------------------------------------------------------------------

     !------------------------------------------------------------------------------------------

     !====================================================
     ! Geometry 1 [** constant **]
     !====================================================

     open(unit=1,file='out.tgyro.geometry.1',status='old',position='append')

     write(1,20) 'r/a','rho','q','s','kappa','s_kappa','delta','s_delta','shift','rmaj/a','b_unit'
     write(1,20) '','','','','','','','','','','(T)'
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             rho(i),&
             q(i),&
             s(i),&
             kappa(i),&
             s_kappa(i),&
             delta(i),&
             s_delta(i),&
             shift(i),&
             r_maj(i)/r_min,&
             b_unit(i)/1e4
     enddo

     close(1)

     !====================================================
     ! Geometry 2 [** constant **]
     !====================================================

     open(unit=1,file='out.tgyro.geometry.2',status='old',position='append')

     write(1,20) 'r/a','zmag/a','dzmag','zeta','s_zeta','volume','d(vol)/dr','<|grad_r|>','rmin'
     write(1,20) '','','','','','(m^3)','(m^2)','','(cm)'
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             zmag(i)/r_min,&
             dzmag(i),&
             zeta(i),&
             s_zeta(i),&
             vol(i)/1e6,&
             volp(i)/1e4,&
             ave_grad_r(i),&
             r(i)

     enddo

     close(1)

     !====================================================
     ! Collision rates and gyroradii
     !====================================================

     open(unit=1,file='out.tgyro.nu_rho',status='old',position='append')

     write(1,20) 'r/a','(a/cs)/t_ii','(a/cs)/t_ee','nue_star','(a/cs)nu_exch','rho_i/a','rho_s/a','frac_ae'
     write(1,20) '','','','','','','',''
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             nui(1,i)*r_min/c_s(i),&
             nue(i)*r_min/c_s(i),&
             nue_star(i),&
             nu_exch(i)*r_min/c_s(i),&
             rho_i(i)/r_min,&
             rho_s(i)/r_min,&
             frac_ae(i)
     enddo

     close(1)

     !====================================================
     ! Electron powers (erg/s = 1e-7 W = 1e-7*1e-6 MW)
     !====================================================

     open(unit=1,file='out.tgyro.power_e',status='old',position='append')

     write(1,20) 'r/a','p_e_fus','p_e_aux','p_e_ohm','p_brem','p_sync','p_line','p_exch[-]','p_expwd[-]','p_e_tot'
     write(1,20) '','(MW)','(MW)','(MW)','(MW)','(MW)','(MW)','(MW)','(MW)','(MW)'
     do i=1,n_r
        write(1,10) &
             r(i)/r_min,&
             +p_e_fus(i)*1e-7*1e-6,&
             +p_e_aux_in(i)*1e-7*1e-6,&
             +p_e_ohmic_in(i)*1e-7*1e-6,&
             -p_brem(i)*1e-7*1e-6,&
             -p_sync(i)*1e-7*1e-6,&
             -p_line(i)*1e-7*1e-6, &
             -p_exch(i)*1e-7*1e-6,&
             -p_expwd(i)*1e-7*1e-6,&
             +p_e(i)*1e-7*1e-6
     enddo

     close(1) 

     !====================================================
     ! Ion powers (erg/s = 1e-7 W = 1e-7*1e-6 MW)
     !====================================================

     open(unit=1,file='out.tgyro.power_i',status='old',position='append')

     write(1,20) 'r/a','p_i_fus','p_i_aux','p_exch','p_expwd','p_i_tot'
     write(1,20) '','(MW)','(MW)','(MW)','(MW)','(MW)'
     do i=1,n_r
        write(1,10) &
             r(i)/r_min,&
             p_i_fus(i)*1e-7*1e-6,&
             p_i_aux_in(i)*1e-7*1e-6,&
             p_exch(i)*1e-7*1e-6,&
             p_expwd(i)*1e-7*1e-6,&
             p_i(i)*1e-7*1e-6 
     enddo

     close(1) 

     !====================================================
     ! Alpha parameters (out.tgyro.alpha) (erg/s = 1e-7 W = 1e-7*1e-6 MW)
     !====================================================

     open(unit=1,file='out.tgyro.alpha',status='old',position='append')

     write(1,20) 'r/a','s_alpha','s_alpha_i','s_alpha_e','frac_ai','frac_ae','E_alpha/E_c'
     write(1,20) '','(MW/cm^3)','(MW/cm^3)','(MW/cm^3)','','',''
     do i=1,n_r
        write(1,10) &
             r(i)/r_min,&
             (s_alpha_i(i)+s_alpha_e(i))*1e-7*1e-6,&
             s_alpha_i(i)*1e-7*1e-6,&
             s_alpha_e(i)*1e-7*1e-6,&
             frac_ai(i),&
             frac_ae(i),&
             e_alpha/e_cross(i)
     enddo

     close(1) 

     !====================================================
     ! gyroBohm factors in physical units
     !====================================================

     open(unit=1,file='out.tgyro.gyrobohm',status='old',position='append')

     write(1,20) 'r/a','Chi_GB','Q_GB','Gamma_GB','Pi_GB','S_GB','c_s'
     write(1,20) '','m^2/s','MW/m^2','10^19/m^2/s','J/m^2','MW/m^3','m/s'
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             chi_gb(i)*1e-4,&
             q_gb(i)*1e-7*1e-6/1e-4,&
             gamma_gb(i)*1e-19/1e-4,&
             pi_gb(i)*1e-7/1e-4,&
             s_gb(i)*1e-7*1e-6/1e-6,&
             c_s(i)/100.0
     enddo

     close(1)

     !====================================================
     ! Electron particle and energy fluxes (out.tgyro.flux_e)
     !====================================================

     open(unit=1,file='out.tgyro.flux_e',status='old',position='append')

     write(1,20) 'r/a','pflux_e_neo','pflux_e_tur','eflux_e_neo','eflux_e_tur',&
          'mflux_e_neo','mflux_e_tur','expwd_e_tur'
     write(1,20) '','(GB)','(GB)','(GB)','(GB)','(GB)','(GB)','(GB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             pflux_e_neo(i),&
             pflux_e_tur(i),&
             eflux_e_neo(i),&
             eflux_e_tur(i),&
             mflux_e_neo(i),&
             mflux_e_tur(i),&
             expwd_e_tur(i)
     enddo

     close(1)

     !====================================================
     ! Electron temperature and density profiles
     !====================================================

     open(unit=1,file='out.tgyro.profile_e',status='old',position='append')

     write(1,20) 'r/a','ne','a/Lne', 'te','a/Lte','betae_unit'
     write(1,20) '','(1/cm^3)','','(keV)','',''
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             ne(i),&
             r_min*dlnnedr(i),&
             te(i)/1e3,&
             r_min*dlntedr(i),&
             betae_unit(i)
     enddo

     close(1)

     do i_ion=1,loc_n_ion

        itag = 'i'//trim(ion_tag(i_ion))     

        !====================================================
        ! Ion particle and energy fluxes
        !====================================================

        open(unit=1,file='out.tgyro.flux_'//itag,status='old',position='append')

        ttag = itag//'_tur'     
        ntag = itag//'_neo'

        write(1,20) 'r/a','pflux_'//ntag,'pflux_'//ttag,'eflux_'//ntag,&
             'eflux_'//ttag,'mflux_'//ntag,'mflux_'//ttag,'expwd_'//ttag
        write(1,20) '','(GB)','(GB)','(GB)','(GB)','(GB)','(GB)','(GB)'
        do i=1,n_r
           write(1,10) r(i)/r_min,&
                pflux_i_neo(i_ion,i),&
                pflux_i_tur(i_ion,i),&
                eflux_i_neo(i_ion,i),&
                eflux_i_tur(i_ion,i),&
                mflux_i_neo(i_ion,i),&
                mflux_i_tur(i_ion,i),&
                expwd_i_tur(i_ion,i)
        enddo

        close(1)

        !====================================================
        ! Ion profiles
        !====================================================

        open(unit=1,file='out.tgyro.profile_'//itag,status='old',position='append')

        write(1,20) 'r/a','n'//itag,'a/Ln'//itag,'t'//itag,'a/Lt'//itag,'beta'//itag//'_unit'
        write(1,20) '','(1/cm^3)','','(keV)','',''
        do i=1,n_r
           write(1,10) r(i)/r_min,&
                ni(i_ion,i),&
                r_min*dlnnidr(i_ion,i),&
                ti(i_ion,i)/1e3,&
                r_min*dlntidr(i_ion,i),&
                beta_unit(i)*ni(i_ion,i)*k*ti(i_ion,i)/pr(i)
        enddo

        close(1)

     enddo ! i_ion

     !====================================================
     ! Ti evolution
     !====================================================

     open(unit=1,file='out.tgyro.evo_ti',status='old',position='append')

     write(1,20) 'r/a','eflux_i_tot','eflux_i_target'
     write(1,20) '','(GB)','(GB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,eflux_i_tot(i),eflux_i_target(i)
     enddo
     close(1)

     !====================================================
     ! Te evolution
     !====================================================

     open(unit=1,file='out.tgyro.evo_te',status='old',position='append')

     write(1,20) 'r/a','eflux_e_tot','eflux_e_target'
     write(1,20) '','(GB)','(GB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,eflux_e_tot(i),eflux_e_target(i)
     enddo
     close(1)

     !====================================================
     ! Er evolution
     !====================================================

     open(unit=1,file='out.tgyro.evo_er',status='old',position='append')

     write(1,20) 'r/a','mflux_tot','mflux_target'
     write(1,20) '','(GB)','(GB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,mflux_tot(i),mflux_target(i)
     enddo
     close(1)

     !====================================================
     ! Density evolution
     !====================================================

     open(unit=1,file='out.tgyro.evo_ne',status='old',position='append')

     write(1,20) 'r/a','pflux_e_tot','pflux_e_target'
     write(1,20) '','(GB)','(GB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,pflux_e_tot(i),pflux_e_target(i)*evo_c(0)
     enddo
     close(1)

     do i_ion=1,loc_n_ion
        itag = 'i'//trim(ion_tag(i_ion))     
        open(unit=1,file='out.tgyro.evo_n'//trim(ion_tag(i_ion)),status='old',position='append')
        write(1,20) 'r/a','pflux_'//itag//'_tot','pflux_'//itag//'_target'
        write(1,20) '','(GB)','(GB)'
        do i=1,n_r
           if (evo_e(i_ion) == 2) then
              write(1,10) r(i)/r_min,pflux_i_tot(i_ion,i),pflux_he_target(i)
           else
              write(1,10) r(i)/r_min,pflux_i_tot(i_ion,i),pflux_e_target(i)*evo_c(i_ion)
           endif
        enddo
        close(1)
     enddo

     !====================================================
     ! Additional profiles
     !====================================================

     open(unit=1,file='out.tgyro.profile',status='old',position='append')

     write(1,20) 'r/a','beta_unit','a*beta_*','a*gamma_p/cs','a*gamma_e/cs','a*f_rot','M=wR/cs'
     write(1,20) '','','','','','',''
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             beta_unit(i),&
             r_min*beta_unit(i)*dlnpdr(i),&
             r_min/c_s(i)*gamma_p(i),&
             r_min/c_s(i)*gamma_p(i)*r(i)/(q(i)*r_maj(i)),&
             r_min*f_rot(i),&
             u00(i)/c_s(i)
     enddo

     close(1)

     ! If we are in test mode, there is no point in going further
     if (gyrotest_flag == 1) return

     !====================================================
     ! Residuals
     !====================================================

     open(unit=1,file='out.tgyro.residual',status='old',position='append')

     write(1,50,advance='no') 'r/a','E(eflux_i)','R1','E(eflux_e)','R2','E(mflux)','R3','E(pflux_e)','R4','E(pflux_i1)','R5'
     do i=2,loc_n_ion
        write(1,55,advance='no') 'E(pflux_i',i,')','R',i+4
     enddo
     write(1,'(a)')

     if (tgyro_relax_iterations == 0) then
        write(1,30) 'ITERATION*: ',i_tran,sum(res_norm)/size(res_norm),flux_counter*n_worker*n_inst
     else 
        write(1,30) 'ITERATION : ',i_tran,sum(res_norm)/size(res_norm),flux_counter*n_worker*n_inst
     endif

     res2(:,:)   = 0.0
     relax2(:,:) = 0.0 

     p = 0
     write(1,40) (0.0,ip=0,2*(loc_n_ion+4))
     do i=2,n_r
        if (loc_ti_feedback_flag == 1) then
           p  = p+1
           res2(i,1) = res_norm(p)
           relax2(i,1) = relax(p)
        endif
        if (loc_te_feedback_flag == 1) then
           p  = p+1
           res2(i,2) = res_norm(p)
           relax2(i,2) = relax(p)
        endif
        if (loc_er_feedback_flag == 1) then
           p  = p+1
           res2(i,3) = res_norm(p)
           relax2(i,3) = relax(p)
        endif
        do is=0,loc_n_ion
           if (evo_e(is) >= 1) then
              p  = p+1
              res2(i,4+is) = res_norm(p)
              relax2(i,4+is) = relax(p)
           endif
        enddo
        write(1,40) &
             r(i)/r_min,(res2(i,ip),relax2(i,ip),ip=1,loc_n_ion+4)
     enddo

     close(1)

     !====================================================
     ! Precision (just the total residual)
     !====================================================

     open(unit=1,file='out.tgyro.prec',status='old',position='append')
     write(1,*) sum(res_norm)/size(res_norm)
     close(1)

     !====================================================
     ! Pedestal
     !====================================================

     if (tgyro_ped_model > 1) then
        open(unit=1,file='out.tgyro.ped',status='old',position='append')
        write(1,20) 'r_*/a','r_top/a','psi_top','n_top','t_top','p_top','zn_top','zt_top','betan'
        write(1,20) '[-]','[-]','[-]','[1/cm^3]','[keV]','[Pa]','[1/cm]','[1/cm]','%'
        write(1,10) r(n_r)/r_min,r_top(1)/r_min,psi_top(1),n_top,t_top/1e3,p_top,zn_top,zt_top,betan_in
        close(1)
     endif

     !====================================================
     ! Control (control.out)
     !====================================================

     open(unit=1,file='out.tgyro.control',status='old',position='append')
     backspace(1)
     write(1,*) i_tran
     close(1)
     !-------------------------------------------------------------------------------------------

     !-------------------------------------------------------------------------------------------
     ! Write progress to screen
     if (i_tran < 10) then
        print '(a,i1,a,1pe10.3,a)', 'INFO: (TGYRO) Finished iteration ',i_tran,' [',sum(res_norm)/size(res),']'
     else if (i_tran < 100) then
        print '(a,i2,a,1pe10.3,a)', 'INFO: (TGYRO) Finished iteration ',i_tran,' [',sum(res_norm)/size(res),']'
     else
        print '(a,i3,a,1pe10.3,a)', 'INFO: (TGYRO) Finished iteration ',i_tran,' [',sum(res_norm)/size(res),']'
     endif
     !-------------------------------------------------------------------------------------------

  endif

  ! Exit if residual less than tolerance
  if (converged .and. i_tran_loop > 1) then
     call MPI_FINALIZE(ierr)
     stop
  endif

  ! Data
10 format(t1,11(1pe13.6,2x))
  ! Text headers
20 format(t2,a,t17,a,t32,a,t47,a,t62,a,t77,a,t92,a,t107,a,t122,a,t137,a,t152,a)
  ! Residual header
30 format(t2,a,i3,1pe12.5,2x,'[',i6,']')
  ! Residuals
40 format(t2,f8.6,14(1x,2(1pe10.3,1x)))
50 format(t2,a,t12,a,t26,a,t35,a,t49,a,t58,a,t72,a,t81,a,t95,a,t104,a,t118,a)
55 format(7x,a,i0,a,3x,a,i0)

end subroutine tgyro_write_data
