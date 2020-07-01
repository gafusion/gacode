!--------------------------------------------------------------
! prgen_map_plasmastate.f90
!
! PURPOSE:
!  Map native plasmastate data onto input.profiles standard.
!--------------------------------------------------------------

subroutine prgen_map_plasmastate

  use prgen_globals
  use expro

  implicit none

  integer :: i,j
  integer :: ix
  real :: dvol
  real, dimension(nx) :: dphidpsi
  real, dimension(:), allocatable :: f1_therm,f2_therm
  real, dimension(:), allocatable :: f1_lump,f2_lump,f3_lump
  real, dimension(:), allocatable :: f1_fast,f2_fast,f3_fast,f4_fast
  real :: z_eff_lump
  real :: m_eff_lump

  !----------------------------------------------------------------------
  ! shot number
  !
  expro_shot = plst_shot_number

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

  print '(a)','INFO: (prgen_map_plasmastate) Found these ion species:'
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

  ! Compute the quasineutrality error with max 9 ions:

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
  ! Map profile data into expro interface variables
  !
  expro_n_exp = nx
  expro_n_ion = plst_dp1_nspec_all-1
  call expro_init(1)
  !
  expro_rho  = rho
  expro_rmin = rmin
  expro_rmaj = rmaj
  expro_te = plst_ts(:,1)
  expro_ne = plst_ns(:,1)*1e-19
  expro_z_eff = plst_zeff(:)
  expro_w0 = omega0(:) 
  expro_ptot = p_tot ! total pressure, thermal + fast ion

  expro_ni = 0.0
  expro_ti = 0.0

  ! ni,ti,vphi
  do i=1,plst_dp1_nspec_all-1
     expro_ni(i,:) = plst_ns(:,i+1)*1e-19
     expro_ti(i,:) = plst_ts(:,i+1)
     if (trim(plst_all_name(i+1)) == 'C') then
        ! COORDINATES: -ipccw accounts for plasmastate toroidal angle convention
        expro_vtor(i,:) = -ipccw*plst_omegat(:)*(rmaj(:)+rmin(:))
     endif
  enddo

  !---------------------------------------------------
  ! Read the cer file and overlay
  !
  if (file_cer /= "null") then
     call prgen_read_cer
     expro_w0 = omega0(:)
     do i=1,plst_dp1_nspec_all
        if (trim(plst_all_name(i+1)) == 'C') then
           expro_vtor(i,:) = vtorc_exp(:)
           expro_vpol(i,:) = vpolc_exp(:)
        endif
     enddo
  endif
  !---------------------------------------------------

  ! Ion identification

  do i=1,expro_n_ion
     expro_mass(i) = plst_m_all(i+1)/1.66e-27
     expro_z(i)    = nint(plst_q_all(i+1)/1.6022e-19)
     call prgen_ion_name(nint(expro_mass(i)),expro_z(i),expro_name(i))     
     if (i+1 > plst_dp1_nspec_th) then
        expro_type(i) = type_fast
     else
        expro_type(i) = type_therm
     endif
  enddo

  !--------------------------------------------------------------------
  ! Convert SWIM integrated powers (W) to densities (W/m^3)
  do i=2,nx
     
     dvol = plst_vol(i)-plst_vol(i-1)

     ! Total powers to electrons and ions "per zone"
     ! Integrated power is thus a partial sum.
     ! Factor of 1e-6 converts plasmastate (W) to input.gacode (MW).

     ! plasmastate supplies totals which may not equal sum of components
     qpow_e(i) = 1e-6*plst_pe_trans(i-1)/dvol
     qpow_i(i) = 1e-6*plst_pi_trans(i-1)/dvol

     ! Collisional exchange
     expro_qei(i) = -1e-6*plst_qie(i-1)/dvol

     ! Fusion power
     expro_qfuse(i) = 1e-6*plst_pfuse(i-1)/dvol
     expro_qfusi(i) = 1e-6*(plst_pfusi(i-1)+plst_pfusth(i-1))/dvol

     ! Heating
     expro_qohme(i)  = 1e-6*plst_pohme(i-1)/dvol
     expro_qbeame(i) = 1e-6*plst_pbe(i-1)/dvol
     expro_qbeami(i) = 1e-6*(plst_pbi(i-1)+plst_pbth(i-1))/dvol
     expro_qrfe(i)   = 1e-6*(plst_peech(i-1)+plst_pmine(i-1))/dvol
     expro_qrfi(i)   = 1e-6*(plst_pmini(i-1)+plst_pminth(i-1)+plst_picth(i-1))/dvol

     ! Radiation
     expro_qsync(i) = 1e-6*plst_prad_cy(i-1)/dvol
     expro_qbrem(i) = 1e-6*plst_prad_br(i-1)/dvol
     expro_qline(i) = 1e-6*plst_prad_li(i-1)/dvol

     ! Momentum source (tq_trans already in Nm)
     !   -ipccw accounts for plasmastate toroidal angle convention
     expro_qmom(i) = -ipccw*plst_tq_trans(i-1)/dvol

     ! Particle source
     expro_qpar_beam(i) = plst_sn_trans(i-1)/dvol
     expro_qpar_wall(i) = 0.0

  enddo
  expro_qei(1)    = expro_qei(2)
  expro_qfuse(1)  = expro_qfuse(2)
  expro_qfusi(1)  = expro_qfusi(2)
  expro_qohme(1)  = expro_qohme(2)
  expro_qbeame(1) = expro_qbeame(2)
  expro_qbeami(1) = expro_qbeami(2)
  expro_qrfe(1)   = expro_qrfe(2)
  expro_qrfi(1)   = expro_qrfi(2)
  expro_qsync(1)  = expro_qsync(2)
  expro_qbrem(1)  = expro_qbrem(2)
  expro_qline(1)  = expro_qline(2)
  expro_qmom(1)   = expro_qmom(2)
  expro_qpar_beam(1)   = expro_qpar_beam(2)
  expro_qpar_wall(1)   = expro_qpar_wall(2)
  qpow_e(1)       = qpow_e(2)
  qpow_i(1)       = qpow_i(2)
  
  qspow_e = expro_qohme+expro_qbeame+expro_qrfe+expro_qfuse-expro_qei &
       -expro_qsync-expro_qbrem-expro_qline
  qspow_i =             expro_qbeami+expro_qrfi+expro_qfusi+expro_qei
  
  ! Manage auxiliary powers
  if (true_aux_flag == 1) then
     print '(a)','INFO: (prgen_map_plasmastate) Setting aux. power as ohmic+NB+RF.'
  else
     ! DEFAULT:
     ! Put missing auxiliary power into cx, giving correct total powers
     ! This assumes that total powers are correct.
     expro_qione = qpow_e-qspow_e
     expro_qioni = qpow_i-qspow_i
 
     print '(a)','INFO: (prgen_map_plasmastate) Setting aux. power as total-fus-rad'
     print '(a)','INFO: (prgen_map_plasmastate) Missing aux. power stored in qione,qioni'
  endif
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Convert SWIM currents (MA) to current densities (MA/m^2)
  expro_johm=0 ; expro_jbs = 0 ; expro_jbstor = 0
  do i=2,nx
     dvol = (plst_vol(i)-plst_vol(i-1))/(2*pi*rmaj(i))
     expro_johm(i) = plst_curr_ohmic(i)/dvol*1e-6
     expro_jbs(i) = plst_curr_bootstrap(i)/dvol*1e-6
     expro_jbstor(i) = (plst_curt(i)-plst_curt(i-1))/dvol*1e-6
  enddo
  expro_johm(1) = expro_johm(2)
  expro_jbs(1) = expro_jbs(2)
  expro_jbstor(1) = expro_jbstor(2)
  !--------------------------------------------------------------------

end subroutine prgen_map_plasmastate
