!--------------------------------------------------------------
! prgen_map_plasmastate.f90
!
! PURPOSE:
!  Map native plasmastate data onto input.profiles standard.
!
! NOTES:
!  See prgen_map_iterdb.f90 for analogous routine for 
!  iterbd data.
!--------------------------------------------------------------

subroutine prgen_map_plasmastate

  use prgen_globals
  use EXPRO_interface

  implicit none

  integer :: i,j
  integer :: ip,ix
  real, dimension(nx) :: dphidpsi
  real, dimension(:), allocatable :: f1_therm,f2_therm
  real, dimension(:), allocatable :: f1_lump,f2_lump,f3_lump
  real, dimension(:), allocatable :: f1_fast,f2_fast,f3_fast,f4_fast
  real :: z_eff_lump
  real :: m_eff_lump

  !--------------------------------------------------------------------
  ! Calculate integrated powers from input sources
  !
  pow_e(1)     = 0.0
  pow_i(1)     = 0.0
  pow_ei(1)    = 0.0
  flow_mom(1)  = 0.0
  flow_beam(1) = 0.0

  pow_e_fus(1) = 0.0 
  pow_i_fus(1) = 0.0
  pow_e_sync(1) = 0.0
  pow_e_brem(1) = 0.0
  pow_e_line(1) = 0.0

  pow_e_ohm(1) = 0.0
  pow_e_nb(1)  = 0.0
  pow_i_nb(1)  = 0.0
  pow_e_rf(1)  = 0.0
  pow_i_rf(1)  = 0.0

  do i=2,nx

     ! Total powers to electrons and ions "per zone"
     ! Integrated power is thus a partial sum.
     ! Factor of 1e-6 converts plasmastate (W) to input.profiles (MW).

     pow_e(i) = pow_e(i-1)+1e-6*plst_pe_trans(i-1)
     pow_i(i) = pow_i(i-1)+1e-6*plst_pi_trans(i-1)

     ! Collisional exchange
     !
     ! plst_qie : power from ions to electrons
     ! pow_ei   : power from electrons to ions
     ! 
     ! Thus, we need negative sign here:
     pow_ei(i) = pow_ei(i-1)-1e-6*plst_qie(i-1)

     ! Reactor power breakdown

     pow_e_fus(i) = pow_e_fus(i-1)+1e-6*plst_pfuse(i-1)
     pow_i_fus(i) = pow_i_fus(i-1)+1e-6*(plst_pfusi(i-1)+plst_pfusth(i-1))

     pow_e_ohm(i) = pow_e_ohm(i-1)+1e-6*plst_pohme(i-1)
     pow_e_nb(i)  = pow_e_nb(i-1)+1e-6*plst_pbe(i-1)
     pow_i_nb(i)  = pow_i_nb(i-1)+1e-6*(plst_pbi(i-1)+plst_pbth(i-1))
     pow_e_rf(i)  = pow_e_rf(i-1)+1e-6*(plst_peech(i-1)+plst_pmine(i-1))
     pow_i_rf(i)  = pow_i_rf(i-1)+1e-6*(plst_pmini(i-1)+plst_pminth(i-1)+plst_picth(i-1))

     pow_e_sync(i) = pow_e_sync(i-1)+1e-6*plst_prad_cy(i-1)
     pow_e_brem(i) = pow_e_brem(i-1)+1e-6*plst_prad_br(i-1)
     pow_e_line(i) = pow_e_line(i-1)+1e-6*plst_prad_li(i-1)

     ! Momentum source
     !
     ! tq_trans already in Nm.
     ! COORDINATES: -ipccw accounts for plasmastate toroidal angle convention
     flow_mom(i) = flow_mom(i-1)+plst_tq_trans(i-1)*(-ipccw)

     ! Particle source
     !
     ! MW/keV = 0.624e22/s

     flow_beam(i) = flow_beam(i-1)+plst_sn_trans(i-1)/0.624e22

  enddo

  ! Check radiated powers
  if (minval(pow_e_sync) < 0.0) then
     print '(a)','WARNING: (prgen) Found negative sync radiation.'
  endif
  if (minval(pow_e_brem) < 0.0) then
     print '(a)','WARNING: (prgen) Found negative brem radiation.'
  endif
  if (minval(pow_e_line) < 0.0) then
     print '(a)','WARNING: (prgen) Found negative line radiation.'
  endif

  ! Manage auxiliary powers
  if (true_aux_flag == 1) then
     pow_e_aux(:) = pow_e_ohm+pow_e_nb+pow_e_rf
     pow_i_aux(:) =          +pow_i_nb+pow_i_rf
     print '(a)','INFO: (prgen) Setting aux. power as ohmic+NB+RF.'
  else
     pow_e_aux(:) = pow_e-(pow_e_fus-pow_ei-pow_e_sync-pow_e_brem-pow_e_line)
     pow_i_aux(:) = pow_i-(pow_i_fus+pow_ei)
     print '(a)','INFO: (prgen) Setting aux. power as total-fus-rad.'
  endif
  !
  pow_e_err = abs(1.0-(pow_e_fus(nx)+pow_e_aux(nx)-pow_ei(nx)- &
       pow_e_sync(nx)-pow_e_brem(nx)-pow_e_line(nx))/pow_e(nx))
  pow_i_err = abs(1.0-(pow_i_fus(nx)+pow_i_aux(nx)+pow_ei(nx))/pow_i(nx))
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! COORDINATES: set sign of poloidal flux 
  dpsi(:) = abs(dpsi(:))*(-ipccw)
  ! 
  ! Convert potential to Omega
  !
  ! omega0 = -c d(Phi)/dpsi
  !
  ! phi [statvolt] = (10/3) phi (kV)
  ! psi [Maxwell]  = 10^8 psi [Weber]
  ! c [cm/s] = 2.9979e10
  !
  ! NOTE: dpsi = plst_psipol-plst_psipol(1)
  !
  call bound_deriv(dphidpsi,plst_epot,dpsi,nx)
  !
  omega0 = -2.9979e10*dphidpsi*(10.0/3.0)/1e8
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Lump main ions and/or fast ions
  !
  allocate(f1_therm(nx))
  allocate(f1_fast(nx))
  allocate(f1_lump(nx))
  allocate(f2_therm(nx))
  allocate(f2_fast(nx))
  allocate(f2_lump(nx))
  allocate(f3_fast(nx))
  allocate(f3_lump(nx))
  allocate(f4_fast(nx))
  !
  f1_therm(:) = 0.0
  f2_therm(:) = 0.0
  do i=2,plst_dp1_nspec_th
     f1_therm(:) = f1_therm(:)+plst_ns(:,i)*plst_q_all(i)/1.6022e-19  
     f2_therm(:) = f2_therm(:)+plst_ns(:,i)*(plst_q_all(i)/1.6022e-19)**2   
  enddo

  print '(a)','INFO: (prgen) Found these ion species:'
  do i=2,plst_dp1_nspec_all
     print '(t6,i2,1x,a)', i-1,trim(plst_all_name(i))
  enddo

  f1_fast(:) = 0.0
  f2_fast(:) = 0.0
  f3_fast(:) = 0.0
  f4_fast(:) = 0.0
  do i=plst_dp1_nspec_th+1,ntop
     f1_fast(:) = f1_fast(:)+plst_ns(:,i)*plst_q_all(i)/1.6022e-19 
     f2_fast(:) = f2_fast(:)+plst_ns(:,i)*(plst_q_all(i)/1.6022e-19)**2
     f3_fast(:) = f3_fast(:)+plst_ns(:,i)*plst_m_all(i)*plst_q_all(i)/1.6022e-19
     f4_fast(:) = f4_fast(:)+plst_ns(:,i)*plst_ts(:,i)*plst_q_all(i)/1.6022e-19
  enddo

  ! Lump main ions
  if (n_lump > 1) then

     f1_lump(:) = 0.0
     f2_lump(:) = 0.0
     f3_lump(:) = 0.0

     ! Add 1 to account for electrons at index 1
     lump_vec(:) = lump_vec(:)+1
     do j=1,n_lump
        i = lump_vec(j)
        ! f1 = sum_i ni*Zi
        f1_lump(:) = f1_lump(:)+plst_ns(:,i)*plst_q_all(i)/1.6022e-19
        ! f2 = sum_i ni*Zi^2
        f2_lump(:) = f2_lump(:)+plst_ns(:,i)*(plst_q_all(i)/1.6022e-19)**2
        ! f3 = sum_i ni*mi*Zi
        f3_lump(:) = f3_lump(:)+plst_ns(:,i)*plst_m_all(i)*plst_q_all(i)/1.6022e-19
     enddo

     z_eff_lump = nint(sum(f2_lump(:)/f1_lump(:))/nx)

     ! Charge-density-weighted effective mass: sum_i (ni*mi*Zi)/(ni*Zi)
     m_eff_lump = sum(f3_lump(:)/f1_lump(:))/nx

     ! Replace first lumped species with lumped density
     plst_ns(:,lump_vec(1))     = f1_lump(:)/z_eff_lump 
     plst_q_all(lump_vec(1))    = z_eff_lump*1.6022e-19  
     plst_m_all(lump_vec(1))    = m_eff_lump  
     plst_all_name(lump_vec(1)) = '[lumped]'

     ! Remove others and restack 
     do j=2,n_lump
        ix = lump_vec(j)-j+2
        do i=ix+1,plst_dp1_nspec_all
           plst_ns(:,i-1)  = plst_ns(:,i)
           plst_ts(:,i-1)  = plst_ts(:,i)
           plst_q_all(i-1) = plst_q_all(i)
           plst_m_all(i-1) = plst_m_all(i)
           plst_all_name(i-1) = plst_all_name(i)
        enddo
        plst_dp1_nspec_all = plst_dp1_nspec_all-1
        plst_dp1_nspec_th = plst_dp1_nspec_th-1
     enddo

  endif

  ! Recompute thermal fractions after lumping
  f1_therm(:) = 0.0
  f2_therm(:) = 0.0
  do i=2,plst_dp1_nspec_th
     f1_therm(:) = f1_therm(:)+plst_ns(:,i)*(plst_q_all(i)/1.6022e-19) 
     f2_therm(:) = f2_therm(:)+plst_ns(:,i)*(plst_q_all(i)/1.6022e-19)**2  
  enddo

  ! Lump fast ions
  if (lump_fast_flag == 1) then

     z_eff_lump = nint(sum(f2_fast(:)/f1_fast(:))/nx)
     ! Charge-density-weighted effective mass: sum_i (ni*mi*Zi)/(ni*Zi)
     m_eff_lump = sum(f3_fast(:)/f1_fast(:))/nx

     ! Replace first lumped species with lumped density
     ix = plst_dp1_nspec_th+1 
     plst_ns(:,ix)     = f1_fast(:)/z_eff_lump 
     ! Charge-density-weighted effective temp: sum_i (ni*Ti*Zi)/(ni*Zi)
     plst_ts(:,ix)     = f4_fast(:)/f1_fast(:)
     plst_q_all(ix)    = z_eff_lump*1.6022e-19  
     plst_m_all(ix)    = m_eff_lump
     plst_all_name(ix) = '[fast]'

     plst_dp1_nspec_all = ix

  endif

  plst_zeff(:) = (f2_therm(:)+f2_fast(:))/(f1_therm(:)+f1_fast(:))

  ! Compute the quasineutrality error with max 5 ions:

  quasi_err = 0.0
  ix = min(plst_dp1_nspec_th+1,n_ion_max+1)
  do i=1,nx
     quasi_err = quasi_err+sum(plst_ns(i,2:ix)*plst_q_all(2:ix))
  enddo
  quasi_err = abs(quasi_err/sum(-plst_ns(:,1)*plst_q_all(1))-1.0)

  deallocate(f1_therm,f2_therm)
  deallocate(f1_lump,f2_lump,f3_lump)
  deallocate(f1_fast,f2_fast,f3_fast,f4_fast)
  !-------------------------------------------------------------------------------

  !---------------------------------------------------------
  ! Map profile data into EXPRO interface variables
  !
  EXPRO_n_exp = nx
  call EXPRO_alloc('./',1)
  !
  EXPRO_rho  = plst_rho
  EXPRO_rmin = rmin(:)
  EXPRO_rmaj = rmaj(:)
  ! COORDINATES: set sign of q
  EXPRO_q = abs(q(:))*ipccw*btccw
  EXPRO_kappa = kappa(:)
  EXPRO_delta = delta(:)
  EXPRO_te = plst_ts(:,1)
  EXPRO_ne = plst_ns(:,1)*1e-19
  EXPRO_z_eff = plst_zeff(:)
  EXPRO_w0 = omega0(:) 
  EXPRO_flow_mom = flow_mom(:)
  EXPRO_pow_e = pow_e(:)
  EXPRO_pow_i = pow_i(:)
  EXPRO_pow_ei = pow_ei(:)
  EXPRO_zeta = zeta(:)
  EXPRO_flow_beam = flow_beam(:)
  EXPRO_flow_wall = 0.0
  EXPRO_zmag = zmag(:)
  EXPRO_ptot = p_tot ! total pressure, thermal + fast ion
  ! COORDINATES: This poloidal flux has correct sign (see above).
  EXPRO_polflux = dpsi(:)

  EXPRO_ni = 0.0
  EXPRO_ti = 0.0

  ! ni,ti,vphi
  do i=1,plst_dp1_nspec_all-1
     ip = reorder_vec(i)
     EXPRO_ni(i,:) = plst_ns(:,ip+1)*1e-19
     EXPRO_ti(i,:) = plst_ts(:,ip+1)
     if (trim(plst_all_name(ip+1)) == 'C') then
        ! COORDINATES: -ipccw accounts for plasmastate toroidal angle convention
        EXPRO_vtor(i,:) = -ipccw*plst_omegat(:)*(rmaj(:)+rmin(:))
     endif
  enddo

  !---------------------------------------------------
  ! Read the cer file and overlay
  !
  if (cer_file /= "null") then
     rho = plst_rho
     allocate(vpolc_exp(nx))
     allocate(vtorc_exp(nx))
     call prgen_read_cer
     EXPRO_w0 = omega0(:)
     do i=1,plst_dp1_nspec_all
        if (trim(plst_all_name(ip+1)) == 'C') then
           EXPRO_vtor(i,:) = vtorc_exp(:)
           EXPRO_vpol(i,:) = vpolc_exp(:)
        endif
     enddo
  endif
  !---------------------------------------------------

  ! Additional powers (fusion and radiation)
  EXPRO_pow_e_fus = pow_e_fus(:)
  EXPRO_pow_i_fus  = pow_i_fus(:)
  EXPRO_pow_e_sync = pow_e_sync(:)
  EXPRO_pow_e_brem = pow_e_brem(:)
  EXPRO_pow_e_line = pow_e_line(:)

  ! Additional powers (external heating)
  EXPRO_pow_e_aux = pow_e_aux(:)
  EXPRO_pow_i_aux = pow_i_aux(:)
  !---------------------------------------------------------

  ! Ion reordering diagnostics

  print '(a)','INFO: (prgen) Created these species:'
  do i=1,n_ion_max
     ip = reorder_vec(i)
     if (ip >= plst_dp1_nspec_all) then
        print '(t6,i2,1x,3(a))',i,'[null]'
     else
        if (ip+1 > plst_dp1_nspec_th) then
           ion_type(i) = type_fast
        else
           ion_type(i) = type_therm
        endif
        print '(t6,i2,1x,a,1x,a)',i,trim(plst_all_name(ip+1)),ion_type(i)
        ion_mass(i) = plst_m_all(ip+1)/1.66e-27
        ion_z(i)    = nint(plst_q_all(ip+1)/1.6022e-19)
        call prgen_ion_name(nint(ion_mass(i)),ion_z(i),ion_name(i))     
     endif
  enddo

end subroutine prgen_map_plasmastate
