!-----------------------------------------------------------
! gyro_write_initdata.f90
!
! PURPOSE:
!  Consolidated code to write all initial data.  
!  file1: profile data
!  file2: units
!  file3: geometry data
!
!-----------------------------------------------------------

subroutine gyro_write_initdata(datafile1,datafile2,datafile3,io)

  use gyro_globals
  use gyro_profile_exp
  use math_constants
  use geo

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile1
  character (len=*), intent(in) :: datafile2
  character (len=*), intent(in) :: datafile3
  !
  real :: kt
  real :: theta
  real :: dr
  integer :: n_wedge
  real, dimension(:), allocatable :: ttmp
  !---------------------------------------------------

  if (output_flag == 0 .or. i_proc /= 0) return

  open(unit=io,file=datafile1,status='replace')

  ! Essential scalars
  write(io,20) n_x
  write(io,20) n_theta_section
  write(io,20) n_pass
  write(io,20) n_trap
  write(io,20) n_energy
  write(io,20) n_theta_plot
  write(io,20) n0
  write(io,20) n_n
  write(io,20) d_n
  write(io,20) n_explicit_damp
  write(io,20) nonlinear_flag
  write(io,20) electron_method
  write(io,20) n_field
  write(io,20) n_ion
  write(io,20) n_kinetic
  write(io,20) n_spec
  write(io,20) n_grid_exp
  write(io,20) boundary_method

  ! Basic profile data (** CORRECT IN VUGYRO **)
  write(io,10) r(:)
  write(io,10) q(:)
  write(io,10) r_s(:)
  write(io,10) q_s(:)
  write(io,10) dlntdr_s(:,:)
  write(io,10) dlnndr_s(:,:)
  write(io,10) tem_s(:,:)
  write(io,10) den_s(:,:)
  write(io,10) rmaj_s(:)/r_s(:)
  write(io,10) delta_s(:)
  write(io,10) zeta_s(:)
  write(io,10) kappa_s(:)
  write(io,10) drmaj_s(:)
  write(io,10) shat_s(:)
  write(io,10) s_delta_s(:)
  write(io,10) s_zeta_s(:)
  write(io,10) s_kappa_s(:)
  write(io,10) zmag_s(:)
  write(io,10) dzmag_s(:)
  write(io,10) beta_unit_s(:)
  write(io,10) gamma_e_s(:)
  write(io,10) gamma_p_s(:) 
  write(io,10) mach_s(:) 
  write(io,10) b_unit_s(:)
  write(io,10) dr_eodr(:)
  write(io,10) z_eff_s(:)
  write(io,10) nu_s(n_spec,:)
  write(io,10) w0_s(:)

  write(io,10) box_multiplier

  write(io,10) lambda(ir_norm,:)
  write(io,10) energy(:)
  write(io,10) lambda_tp(ir_norm)
  write(io,10) krho_collect(:)
  write(io,10) rhos_norm
  write(io,10) z(:)

  ! Recent additions
  write(io,20) n_theta_plot*n_theta_mult

  ! Number of flux moments
  write(io,20) n_moment

  close(io)

  !================================================================

  open(unit=io,file=datafile2,status='replace')

  ! kT in MJ (note the conversion 1.6022e-22 MJ/keV)
  kt = 1.6022e-22*tem_norm

  write(io,26) 2.0*kg_proton,'m_ref (kg)'
  write(io,26) b_unit_norm,'b_unit (Tesla)'
  write(io,26) a_meters,'a (m)'
  write(io,26) csda_norm,'csD/a (1/s)'
  write(io,26) csda_norm*a_meters,'csD (m/s)'
  write(io,26) tem_norm,'Te (keV)'
  write(io,26) den_norm,'ne (10^19/m^3)'
  write(io,26) rhos_norm*a_meters,'rho_sD (m)'
  write(io,26) csda_norm*(rhos_norm*a_meters)**2,&
       'chi_gBD (m^2/s)'

  write(io,26) 1e19*den_norm*(csda_norm*a_meters)*rhos_norm**2/0.624e22,&
       'Gamma_gBD (0.624e22/m^2/s) = (MW/keV/m^2)'

  write(io,26) 1e19*den_norm*(csda_norm*a_meters)*kt*rhos_norm**2,&
       'Q_gBD (MJ/m^2/s) = (MW/m^2)'

  write(io,26) 1e19*den_norm*a_meters*kt*rhos_norm**2*1e6,&
       'Pi_gBD (J/m^2) = (Nm/m^2)'

  write(io,26) 1e19*den_norm*csda_norm*kt*rhos_norm**2,&
       'S_gBD (MJ/m^3/s) = (MW/m^3)'

  close(io)

  !================================================================

  open(unit=io,file=datafile3,status='replace')

  n_wedge = n_theta_plot*n_theta_mult
  allocate(ttmp(n_wedge))
  
  do i=1,n_x

     if (flat_profile_flag == 0) then

        ! All profiles are global and radial variation 
        ! is consistent

        GEO_rmin_in      = r_s(i)
        GEO_rmaj_in      = rmaj_s(i)
        GEO_drmaj_in     = drmaj_s(i)
        GEO_zmag_in      = zmag_s(i)
        GEO_dzmag_in     = dzmag_s(i)
        GEO_q_in         = q_s(i)
        GEO_s_in         = shat_s(i)
        GEO_kappa_in     = kappa_s(i)
        GEO_s_kappa_in   = s_kappa_s(i)
        GEO_delta_in     = delta_s(i)
        GEO_s_delta_in   = s_delta_s(i)
        GEO_zeta_in      = zeta_s(i)
        GEO_s_zeta_in    = s_zeta_s(i)
        GEO_beta_star_in = beta_star_s(i)

     else

        ! Profiles are flat and so some parameters need
        ! to be linearly extrapolated.

        dr = r(i)-r(ir_norm)

        GEO_rmin_in  = r(i)
        GEO_rmaj_in  = rmaj_s(ir_norm)+drmaj_s(ir_norm)*dr
        GEO_drmaj_in = drmaj_s(ir_norm)
        GEO_zmag_in  = zmag_s(ir_norm)+dzmag_s(ir_norm)*dr
        GEO_dzmag_in = dzmag_s(ir_norm)
        GEO_q_in     = q(i)
        GEO_s_in     = shat_s(ir_norm)
        GEO_kappa_in = kappa_s(ir_norm)+&
             kappa_s(ir_norm)*s_kappa_s(ir_norm)/r(ir_norm)*dr
        GEO_s_kappa_in = s_kappa_s(ir_norm)
        GEO_delta_in   = delta_s(ir_norm)+s_delta_s(ir_norm)/r(ir_norm)*dr
        GEO_s_delta_in = s_delta_s(ir_norm)
        GEO_zeta_in    = zeta_s(ir_norm)+s_zeta_s(ir_norm)/r(ir_norm)*dr
        GEO_s_zeta_in  = s_zeta_s(ir_norm)
        GEO_beta_star_in = beta_star_s(ir_norm)

     endif

     GEO_fourier_in(:,:) = a_fourier_geo_s(:,0:n_fourier_geo,i)

     do j=1,n_wedge
        ttmp(j) = -pi+(j-1)*pi_2/n_wedge
     enddo
     if (n_wedge == 1) ttmp = 0.0

     call geo_interp(n_wedge,ttmp,.true.)
     
     do j=1,n_wedge

        write(io,10) GEO_nu(j)
        write(io,10) GEO_gsin(j)
        write(io,10) GEO_gcos1(j)
        write(io,10) GEO_gcos2(j)
        write(io,10) GEO_usin(j)
        write(io,10) GEO_ucos(j)
        write(io,10) GEO_b(j)
        write(io,10) GEO_g_theta(j)
        write(io,10) GEO_grad_r(j)
        write(io,10) GEO_gq(j)
        write(io,10) GEO_captheta(j)
        write(io,10) GEO_theta_nc(j)

     enddo ! j
  enddo ! i

  close(io)

  deallocate(ttmp)
  
  if (debug_flag == 1 ) then
     print *,'[gyro_write_initdata done]'
  endif

10 format(1pe15.8)
20 format(i4)
26 format(t2,1pe15.8,2x,a) 

end subroutine gyro_write_initdata
