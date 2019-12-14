!
      PROGRAM tglf_driver_mpi
!
      USE tglf_pkg
      USE tglf_tg
      USE tglf_mpi
!
      IMPLICIT NONE
! local variables
      INTEGER :: i,j,ierr
      REAL :: cmult,a_unit=100.0,B_unit=10000.0
!
  ! initialize MPI
      call MPI_INIT(ierr)

      iCommTglf = MPI_COMM_WORLD
      call MPI_COMM_RANK(iCommTglf,iProcTglf,ierr)
      call MPI_COMM_SIZE(iCommTglf,nProcTglf,ierr)
      write(*,*)"ierr=",ierr
      write(*,*)"iProcTglf=",iProcTglf
      write(*,*)"nProcTglf=",nProcTglf

      OPEN (unit=3,file='./tglfin',status='old')
      READ(3,nml=tglfin)
      CLOSE(3)
!
      CALL put_signs(sign_Bt_tg,sign_It_tg)
      CALL put_rare_switches(theta_trapped_tg,park_tg,ghat_tg,gchat_tg, &
       wd_zero_tg,Linsker_factor_tg,gradB_factor_tg,filter_tg,damp_psi_tg,damp_sig_tg)
      CALL put_switches(iflux_tg,use_bper_tg,use_bpar_tg,use_mhd_rule_tg,use_bisection_tg, &
       ibranch_tg,nmodes_tg,nbasis_max_tg,nbasis_min_tg,nxgrid_tg,nky_tg)
!
      CALL put_species(ns_tg,zs_tg,mass_tg) 
!
      CALL put_model_parameters(adiabatic_elec_tg,alpha_e_tg,alpha_p_tg,alpha_mach_tg    &
      ,alpha_quench_tg,alpha_zf_tg,xnu_factor_tg,debye_factor_tg                       &
      ,etg_factor_tg,sat_rule_tg,kygrid_model_tg,xnu_model_tg              &
      ,vpar_model_tg,vpar_shear_model_tg,write_wavefunction_flag_tg)
!      
      CALL put_kys(ky_tg)
!
      CALL put_gaussian_width(width_max_tg,width_min_tg,nwidth_tg    &
      ,find_width_tg)
!
      CALL put_gradients(rlns_tg,rlts_tg,vpar_shear_tg,vexb_shear_tg)
!
      CALL put_averages(taus_tg,as_tg,vpar_tg,vexb_tg,betae_tg,xnue_tg,zeff_tg,debye_tg)
!
      if(igeo_tg.eq.0)then
        CALL put_s_alpha_geometry(rmin_tg,rmaj_tg,q_tg,shat_tg,alpha_tg, &
         xwell_tg,theta0_tg,b_model_tg,ft_model_tg)
      elseif(igeo_tg.eq.1)then
!        q_prime_tg = shat_tg*(q_tg/rmin_tg)**2
        CALL put_Miller_geometry(rmin_tg,rmaj_tg,zmaj_tg,drmindx_tg,drmajdx_tg,dzmajdx_tg, &
         kappa_tg,s_kappa_tg,delta_tg,s_delta_tg,zeta_tg,s_zeta_tg,q_tg,q_prime_tg,p_prime_tg,kx0_tg)
      elseif(igeo_tg.eq.2)then
        CALL put_Fourier_geometry(q_tg,q_prime_tg,p_prime_tg,nfourier_tg,fourier_tg)
      elseif(igeo_tg.eq.3)then
         OPEN(unit=10,file='test132010.surf',status='old')
         READ(10,*)n_elite_tg,n_surface_tg
         n_elite_tg = n_elite_tg-1
         if(j_surface_tg.gt.n_surface_tg)j_surface_tg=n_surface_tg
         do j=1,j_surface_tg
           do i=0,n_elite_tg
             READ(10,'(1p,3d26.16)')R_elite_tg(i),Z_elite_tg(i),Bp_elite_tg(i)
           enddo
         enddo
         CLOSE(10)
         a_unit = 1.0
         B_unit = 1.0
         do i=0,n_elite_tg
           R_elite_tg(i) = R_elite_tg(i)/a_unit
           Z_elite_tg(i) = Z_elite_tg(i)/a_unit
           Bp_elite_tg(i) = Bp_elite_tg(i)/B_unit
         enddo
!         n_elite_tg=300
!         dt = 2.0*pi/REAL(n_elite_tg)
!         do i=0,n_elite_tg
!           t = REAL(i)*dt
!           R_elite_tg(i) = rmaj_tg + rmin_tg*COS(t)
!           Z_elite_tg(i) = rmin_tg*SIN(t)
!           Bp_elite_tg(i) = rmin_tg/(q_tg*R_elite_tg(i))
!         enddo
         CALL put_ELITE_geometry(n_elite_tg,q_tg,q_prime_tg,p_prime_tg,R_elite_tg,Z_elite_tg,Bp_elite_tg)
      else
        write(*,*)"igeo_tg invalid",igeo_tg
!        stop
      endif
!
      call tglf_setup_geometry

      CALL put_eikonal(new_eikonal_tg)

      CALL tglf_TM_mpi

      if(iProcTglf.eq.0)then
         write(*,*)"R2_ave = ",get_R2_ave()
         write(*,*)"B2_ave =",get_B2_ave()*drmindx_tg**2
         write(*,*)"Rbt_ave = ",get_Rbt_ave()*drmindx_tg
         write(*,*)"B_ave=",get_B_ave()
         write(*,*)"Bt_ave=",get_Bt_ave()
         do j=1,ns_tg
         cmult = 0.7967
         if(j.eq.2)cmult=1.207
!         cmult = 1.059
!         if(j.eq.2)cmult=1.385
          cmult = 1.0
          write(*,*)"total flux for species ",j
          write(*,*)"electrostatic"
!          write(*,*)"particle_flux=",get_particle_flux(j,1)/rlns_tg(j)          
!          write(*,*)"energy_flux=",cmult*get_energy_flux(j,1)/rlts_tg(j)
!          write(*,*)"particle_flux=",get_particle_flux(j,1)*drmindx_tg**2         
!          write(*,*)"energy_flux=",cmult*get_energy_flux(j,1)*drmindx_tg**2
          write(*,*)"particle_flux=",get_particle_flux(j,1)        
          write(*,*)"energy_flux=",get_energy_flux(j,1)
          write(*,*)"stress_par=",get_stress_par(j,1)
          write(*,*)"stress_tor=",get_stress_tor(j,1)
!          write(*,*)"n_bar_sum=",get_n_bar_sum(j)
!          write(*,*)"t_bar_sum=",get_t_bar_sum(j)
!          if(gamma_e_tg.ne.0.0)write(*,*)"stress_par=",get_stress_par(j,1)/gamma_e_tg
!          if(gamma_e_tg.ne.0.0)write(*,*)"stress_par=",get_stress_par(j,1)*drmindx_tg**2
!           write(*,*)"stress_tor=",get_stress_tor(j,1)
!          if(gamma_p_tg.ne.0.0)then
!           write(*,*)"stress_par=",cmult*get_stress_par(j,1)/gamma_p_tg
!           write(*,*)"stress_per=",cmult*(get_stress_par(j,1)-get_stress_tor(j,1)/rmaj_tg)/gamma_p_tg
!           write(*,*)"stress_tor=",cmult*get_stress_tor(j,1)
!           write(*,*)"Prandtl = ",get_stress_tor(j,1)*rlts_tg(2)/(get_energy_flux(j,1)*gamma_p_tg*rmaj_tg)
!          endif
          write(*,*)"exchange=",get_exchange(j,1)
         if(use_bper_tg)then
          write(*,*)"B_per"
          write(*,*)"particle_flux=",get_particle_flux(j,2)          
          write(*,*)"energy_flux=",get_energy_flux(j,2)
          write(*,*)"stress_par=",get_stress_par(j,2)
          write(*,*)"stress_tor=",get_stress_tor(j,2)
          write(*,*)"exchange=",get_exchange(j,2)
         endif
         if(use_bpar_tg)then
          write(*,*)"B_par"
          write(*,*)"particle_flux=",get_particle_flux(j,3)          
          write(*,*)"energy_flux=",get_energy_flux(j,3)
          write(*,*)"stress_par=",get_stress_par(j,3)
          write(*,*)"stress_tor=",get_stress_tor(j,3)
          write(*,*)"exchange=",get_exchange(j,3)
         endif
          write(*,*)"q_low=",get_q_low(j)*drmindx_tg**2
          write(*,*)"q_high=",(get_energy_flux(j,1)+get_energy_flux(j,2) - get_q_low(j))*drmindx_tg**2
!          write(*,*)"n_bar_sum=",get_n_bar_sum(j)
!          write(*,*)"t_bar_sum=",get_t_bar_sum(j)
        enddo
!
        CALL write_tglf_input
      endif 
!        
      new_eikonal_tg=.FALSE.
      CALL put_eikonal(new_eikonal_tg)
      CALL tglf_TM_mpi
     if(iProcTglf.eq.0)then
       write(*,*)"second call to TGLF_TM put_eikonal test"
       do j=1,ns_tg
          write(*,*)"total flux for species ",j
          write(*,*)"particle_flux=",get_particle_flux(j,1)          
          write(*,*)"energy_flux=",get_energy_flux(j,1)
          write(*,*)"stress_par=",get_stress_par(j,1)
          write(*,*)"stress_tor=",get_stress_tor(j,1)
          write(*,*)"n_bar_sum=",get_n_bar_sum(j)
          write(*,*)"t_bar_sum=",get_t_bar_sum(j)
       enddo  
     endif        
!
      STOP
!  
      END   !program tglf_driver_mpi
