       MODULE tglf_tg
!
      IMPLICIT NONE
      SAVE
!
      INTEGER,PARAMETER :: nsm=3, nfm=3
!
! external controls
!
      LOGICAL :: change_defaults_tg=.TRUE.
      INTEGER :: nfields_tg=1 
      INTEGER :: igeo_tg=0
!
! user specified normalizations for TGLF 
!
      REAL :: a0 ! unit of length
      REAL :: T0 ! unit of temperature
      REAL :: N0 ! unit of density
      REAL :: B0 ! unit of magnetic field
      REAL :: m0 ! unit of mass
!
! input data for the tglf model
!
      LOGICAL :: iflux_tg=.TRUE.         ! if .TRUE. compute quasilinear fluxes 
      LOGICAL :: use_bper_tg=.FALSE.     ! include perpendicular magnetic field fluctuations
      LOGICAL :: use_bpar_tg=.FALSE.     ! include parallel magnetic field fluctuations
      LOGICAL :: use_bisection_tg=.FALSE.! use bisection search method to find width that maximizes growthrate
      LOGICAL :: find_width_tg=.TRUE.    ! find the width that maximizes the growthrate
      LOGICAL :: adi_elec_tg=.FALSE.     ! if .TRUE. use adiabatic electron approximation
      LOGICAL :: new_eikonal_tg=.TRUE.   ! if .TRUE. compute the eikonal, if .FALSE. re-use eikonal from previous call
      INTEGER :: nspecies_tg=2           ! number of species electrons + ions limited to (max=3)
      INTEGER :: ibranch_tg=-1           ! if -1 order linear mode spectrum by growthrate, 
      INTEGER :: nmodes_tg=2             ! number of ustable modes to use in computing fluxes (max=4)
      INTEGER :: nky_tg=15               ! number of modes for high-k sector
      INTEGER :: nb_max_tg=4             ! maximum number of poloidal hermite basis functions
      INTEGER :: nb_min_tg=1             ! minumum number of poloidal hermite basis functions
      INTEGER :: nx_tg=16                ! number of nodes for the guass-hermite quadrature
      INTEGER :: nwidth_tg=21            ! maximum number of widths to use in search for maximum growthrate
      INTEGER :: b_model_tg=0            ! do 1 or don't 0 use B(theta) factor in the s-alpha geometry k_per**2
      INTEGER :: ft_model_tg=1           ! trapped fraction model for s-alpha geometry
      INTEGER :: kygrid_model_tg=1       ! select version of ky-grid to use
      INTEGER :: xnu_model_tg            ! select version of trapped-passing boundary ei-collision model to use
      REAL :: xky0_tg=0.3                ! k_theta*rhos for a single k linear stability call to tglf
      REAL :: width_max_tg=1.65          ! maximum width to use for the search for the maximum growthrate
      REAL :: width_min_tg=0.3           ! minimum width to use for the search for the maximum growthrate
      REAL :: park_tg=1.0                ! multiplies the parallel gradient operator
      REAL :: ghat_tg=1.0                ! multiplies the toridal drift terms (not closure terms)
      REAL :: gchat_tg=1.0               ! multiplies the toriodal drift closure terms
      REAL :: betae_tg=0.0               ! electron beta
      REAL :: xnuei_tg=0.0               ! electron-ion collision frequency 
      REAL :: rmin_tg=0.5                ! minor radius of local flux surface
      REAL :: rmaj_tg=3.0                ! major radius of local flux surface
      REAL :: zeff_tg=1.0                ! effective charge of ions 
      REAL :: debye_tg=0.0               ! ratio of Debye length to the gyro-radius scale rhos0
      REAL :: taus_tg(nsm)               ! temperature T/T0 for each species 
      REAL :: as_tg(nsm)                 ! density N/N0 for each species
      REAL :: rlns_tg(nsm)               ! -a0*dLog(N)/drmin for each species 
      REAL :: rlts_tg(nsm)               ! -a0*dLog(T)/drmin for each species
      REAL :: zs_tg(nsm)                 ! charge/e for each species
      REAL :: mass_tg(nsm)               ! mass/m0 for each species
      REAL :: alpha_p_tg=0.0             ! multiplies the parallel velocity shear
      REAL :: alpha_e_tg=0.0             ! multiplies the ExB velocity shear
      REAL :: gamma_e_tg=0.0             ! ExB velocity shear
      REAL :: gamma_p_tg=0.0             ! parallel velocity shear
      REAL :: vexb_mach_tg=0.0           ! ExB velocity mach number
      REAL :: vpar_mach_tg=0.0           ! parallel velocity mach number
      REAL :: wd_zero_tg=0.1             ! cutoff for curvature drift eigenvalues ABS(wd)>=wd_zero
      REAL :: Linsker_tg=0.0             ! multiplies Linker terms (parallel gradient of FLR terms)
      REAL :: gradB_tg=0.0               ! multiplies parallel gradient of Log(B) terms
      REAL :: theta_trap_tg=0.7          ! parameter for trapped fraction model
      REAL :: xnuei_fac_tg=1.0           ! multiplies the trapped-passing boundary terms in the electron-ion collisions
      REAL :: debye_fac_tg=1.0           ! multiplies the debye length
      REAL :: etg_fac_tg = 4.0           ! exponent parameter for ETG saturation rule
      REAL :: filter_tg = 2.0            ! filter to remove spurious high-frequency magnetic fluctuations
      REAL :: x_psi_tg=0.0               ! shift polorization current due to magnetic fluctuations
! Shifted cicle geometry inputs
      REAL :: theta0_tg=0.0              ! kx/(shat ky) radial mode number
      REAL :: shat_tg=1.0                ! magnetic shear
      REAL :: alpha_tg=0.0               ! normalized pressure gradient
      REAL :: xwell_tg=0.0               ! magnetic well parameter
      REAL :: q_tg=2.0                   ! safety factor
! Miller geometry inputs
      REAL :: delta_tg=0.0               ! triangularity of flux surface 
      REAL :: kappa_tg=1.0               ! elongation of flux surface
      REAL :: shift_tg=0.0               ! dRmajor/drminor
      REAL :: s_delta_tg=0.0             ! shear in triagularity
      REAL :: s_kappa_tg=0.0             ! shear in elongation
      REAL :: q_prime_tg=0.0             ! safety factor derivative wrt poloidal flux
      REAL :: p_prime_tg=0.0             ! pressure derivative wrt poloidal flux
!
! output results 
!
      REAL :: v_bar_sum_tg               ! sum over the TM spectrum of mode intensities (if zero no unstable modes were found)  
      REAL :: particle_flux_tg(nsm,nfm)  ! total particle flux for each species and field
                                         ! (nfm = electric potential, perpendicular magnetic field, parallel magnetic field)
      REAL :: energy_flux_tg(nsm,nfm)    ! total energy flux for each species and field
!
      END module tglf_tg
!
!--------------------------------------------------------
!
      SUBROUTINE tglf_startup
!
! This subroutine calls the put_ routines that set global parameters
! and switches for the TGLF transport model. It only needs to be called 
! if the default values are to be changed by reading in the file tglfin.
!
      USE tglf_tg
      USE tglf_pkg
      IMPLICIT NONE
! 
! namelist ....................................................

      NAMELIST /tglfin/ adi_elec_tg, find_width_tg, new_eikonal_tg, &
        nb_max_tg, nb_min_tg, nx_tg, ibranch_tg, nspecies_tg, &
        nmodes_tg, iflux_tg, xky0_tg, width_max_tg, width_min_tg, &
        nwidth_tg, park_tg, ghat_tg, gchat_tg, alpha_e_tg, gamma_e_tg, &
        alpha_p_tg, gamma_p_tg, igeo_tg, theta_trap_tg, &
        theta0_tg,taus_tg,as_tg,rlns_tg,rlts_tg,mass_tg,zs_tg, &
        rmin_tg, rmaj_tg, use_bisection_tg, vexb_mach_tg, vpar_mach_tg, &
        q_tg, xnuei_tg, wd_zero_tg, betae_tg, shat_tg, alpha_tg, &
        xwell_tg, kappa_tg, s_kappa_tg, delta_tg, s_delta_tg, shift_tg,  &
        zeff_tg, debye_tg, use_bper_tg, use_bpar_tg,q_prime_tg, &
        p_prime_tg, filter_tg, Linsker_tg, gradB_tg, x_psi_tg, &
        b_model_tg, ft_model_tg, xnuei_fac_tg, debye_fac_tg, &
        nky_tg,etg_fac_tg
!
      OPEN (unit=3,file='tglfin',status='old')
      READ(3,nml=tglfin)
      CLOSE(3)
!
! sequence of calls to set tglf switches and model parameters for all future calls
! all of the put_ calls are optional if defaults are to be used
!
      CALL put_species(nspecies_tg,zs_tg,mass_tg) 
!      
      CALL put_kys(xky0_tg)
!
      CALL put_gaussian_width(width_max_tg,width_min_tg,nwidth_tg    &
      ,find_width_tg)
!
      CALL put_switches(iflux_tg,use_bper_tg,use_bpar_tg,use_bisection_tg, &
       ibranch_tg,nmodes_tg,nb_max_tg,nb_min_tg,nx_tg,nky_tg,  &
       kygrid_model_tg,xnu_model_tg)
!
      CALL put_model_parameters(adiabatic_elec_tg,alpha_e_tg,alpha_p_tg    &
      ,alpha_quench_tg,xnu_factor_tg,debye_factor_tg                       &
      ,etg_factor_tg,sat_rule_tg,kygrid_model_tg,xnu_model_tg              &
      ,vpar_model_tg,vpar_shear_model_tg,write_wavefunction_flag_tg)
!
      RETURN
      END SUBROUTINE tglf_startup
!
!----------------------------------------------------------
!
      SUBROUTINE tglf_TM_interface
!
! This is a transport code interface for the TGLF transport model.
! It is assumed that all of the input variables have been set before calling 
! this interface and are stored in the module tglf_tg
! The TGLF dimensionless particle and energy fluxes are output and stored in tglf_tg
!
!
      USE tglf_tg
      USE tglf_pkg
      IMPLICIT NONE
      INTEGER :: j,k
!
      CALL put_averages(taus_tg,as_tg,vpar_tg,vexb_tg,betae_tg,xnue_tg,zeff_tg,debye_tg)
!
!
      CALL put_gradients(rlns_tg,rlts_tg,vpar_shear_tg,vexb_shear_tg)
!
!
      if(igeo_tg.eq.0)then
        CALL put_s_alpha_geometry(rmin_tg,rmaj_tg,q_tg,shat_tg,alpha_tg, &
         xwell_tg,theta0_tg,b_model_tg,ft_model_tg)
      elseif(igeo_tg.eq.1)then
        CALL put_Miller_geometry(rmin_tg,rmaj_tg,q_tg,q_prime_tg,p_prime_tg, &
          shift_tg,kappa_tg,s_kappa_tg,delta_tg,s_delta_tg)
      endif
      CALL put_eikonal(new_eikonal_tg)
      CALL tglf_TM
      new_eikonal_tg=.TRUE. ! reset to default
!   
      v_bar_sum_tg = get_v_bar_sum()  ! sum of mode amplitudes to check if there were any unstable modes
!
! collect normalized fluxes
!
      nfields_tg=1
      if(use_bper_tg)nfields_tg=nfields_tg+1
      if(use_bpar_tg)nfields_tg=nfields_tg+1
!
      do j=1,nspecies_tg
        do k=1,nfields_tg
         particle_flux_tg(j,k) = get_particle_flux(j,k)          
         energy_flux_tg(j,k) = get_energy_flux(j,k)
        enddo
      enddo  
!
      RETURN
!  
      END   !SUBROUTINE tglf_TM_interface
!
      PROGRAM tglf_TM_driver
!***********************************************************************************
! This driver is an example of how a transport code might use the TGLF transport model
! It contains transport code specific coding that will need to be adapted by the installer.
! For this example the units of the transport code are assumed to be CGS except temperature 
! which is in ev like in the NRL formulary.
! The module tglf_tg provides the storage for the input and output for the TGLF transport model. 
! This example driver computes the TGLF inputs from the local plasma profiles and magnetic geometry
! and defines the normalization used for dimensionless TGLF inputs and outputs.
! The installer can choose these normalizations to suit their preference.
!***********************************************************************************
      USE tglf_tg
      USE tglf_pkg
      IMPLICIT NONE
      INTEGER, PARAMETER ::  ngrid=51
      INTEGER :: i,j,k,j0,bt_flag
      REAL :: Bunit    ! local magnetic field normalizations
      REAL :: e0,c0,k0,mp,pi     ! natural constants
      REAL :: cs0,rhos0,omega0                ! plasma units of velocity, gyroradius, gyrofrequency
      REAL :: rhostar2                        ! gyrobohm reduction factor (rhos/a0)**2 
      REAL :: rho(ngrid)                      ! radial coordinate (cm): toroidal magnetic flux = B0 Pi rho^2  
      REAL :: ne(ngrid),ni(ngrid),nz(ngrid)   ! electron, ion, impurity densities (1/cm^3) on x-grid
      REAL :: te(ngrid),ti(ngrid),tz(ngrid)   ! electron, ion, impurity temperatures (ev) on x-grid
      REAL :: ptot(ngrid)                     
      REAL :: vexb(ngrid),vpar(ngrid)         ! exb and parallel velocity (cm/sec) in x-grid
      REAL :: rmajor(ngrid),rminor(ngrid),elongation(ngrid),triangularity(ngrid),q(ngrid) ! geometry
      REAL :: lnlamda,taue,drhodr,dr,cexb     ! local variables                                   
!
! test case setup
!
! cgs natural constants
      k0 = 1.6022E-12        ! energy of 1ev (erg/ev)
      e0 = 4.8032E-10        ! elementary charge in statcoulombs
      c0 = 2.9979E+10        ! speed of light in cm/sec
      mp = 1.6726E-24        ! proton mass in gm
      pi = 3.141592653589793 
! user defined units for the transport code (cgs + ev in this case)
      a0 = 100.0   ! length unit in cm
      T0 = 1000.0  ! temperature unit in ev  
      N0 = 1.0E13  ! density unit in 1/cm^3 
!
! test case setup 
! profiles and geometry arrays would normally be filled with real data outside of this routine
!
      do i=1,ngrid
        rho(i) = a0*REAL(i-1)/REAL(ngrid-1)
        ne(i) = 5.0*N0*EXP(-rho(i)/a0)
        ni(i) = 5.0*N0*EXP(-rho(i)/a0)
        nz(i) = (ne(i)-ni(i))/6.0  ! charge balance for carbon impurity
        te(i) = 5.0*T0*EXP(-3.0*rho(i)/a0)
        ti(i) = 5.0*T0*EXP(-3.0*rho(i)/a0)
        tz(i) = 5.0*T0*EXP(-3.0*rho(i)/a0)
        ptot(i) = ne(i)*te(i)+(ni(i)+nz(i))*ti(i)    ! should actually use the total pressure from the MHD equilibrium
! The following use Miller geometry definitions
        elongation(i) = 1.0          ! elongation of the flux surface
        triangularity(i) = 0.0       ! triangularity of the flux surface
        q(i) = 1.1 + 3.0*rho(i)/a0   ! safety factor
        rmajor(i) = 300.0            ! major radius in cm
        rminor(i) = rho(i)           ! minor radius in cm 
!note  vexb = c R(0) q/(rho B0) dphi/drho = c R(0) q/(rminor Bunit) dphi/dr where phi = electostatic potential
        vexb(i) = 2.0E7*(1.0-0.8*rho(i)/a0)          ! cm/sec
        vpar(i) = vexb(i)*rho(i)/(q(i)*rmajor(1))    ! cm/sec  
        vpar(i) = 0.0  ! not active yet
      enddo
!
! end test case setup
!
! choose TGLF switches and model parameters using external tglfin file read by tglf_startup
! note: tglf_startup will overwrite all of the *_tg inputs parameters so it has to be called before setting
! the inputs below
!
      if(change_defaults_tg)CALL tglf_startup  ! only needs to be called once 
      change_defaults_tg=.FALSE.
!
! here it is assumed that j0 is the index of the flux surface labled by rho(j0)
! it could be set outside of this routine
! The flux will be computed on the half grid point 1/2(rho(j0)+rho(j0+1)) in this example
! finite differences are used for derivatives but TGLF is local to a flux
! surface so any numerical derivative scheme can be used as required by the transport code solution method.
!
      j0=25 ! for example
!
! The scales N0, T0, a0, m0 used to normalize the TGLF inputs are arbitrary. But B0 is specific.
! For this example we use the GYRO conventions
      N0 = (ne(j0+1)+ne(j0))/2.0   ! density scale used by GYRO
      T0 = (te(j0+1)+te(j0))/2.0   ! temperature scale used by GYRO
      a0 = rminor(ngrid)           ! length scale used by GYRO
      m0 = 2.0*mp                  ! main ion mass in gm (deuterium)
      B0 = 2.0E4    ! magnetic field defined so that toroidal magnetic flux = B0 pi rho^2 in Gauss.
!
!local magnetic field unit
!
      dr = (rminor(j0+1)-rminor(j0))/a0    ! gradients are taken with respect to the minor radius even for s-alpha geometry
      drhodr = (rho(j0+1)-rho(j0))/(rminor(j0+1)-rminor(j0))
      Bunit = B0*drhodr*(rho(j0+1)+rho(j0))/(rminor(j0+1)+rminor(j0))  ! Miller geometry magnetic field unit
! note for s-alpha geometry (igeo_tg=0) you can use Bunit = B0 instead here by setting bt_flag=0.
! This was the default for GLF23 while still using rminor for derivatives. 
! You could also use rho for derivatives for s-alpha geometry but this has not been the convention.
      bt_flag = 1
      if(igeo_tg.eq.0.and.bt_flag.eq.0)bunit=B0
!
! derived units for the plasma
!
      cs0 = SQRT(k0*T0/m0)       ! thermal velocity unit cm/sec
      omega0 = e0*Bunit/(m0*c0)  ! gyrofrequency unit 1/sec
      rhos0 = cs0/omega0         ! gyroradius unit cm
!
!
! local field averages
!
      as_tg(1) = (ne(j0)+ne(j0+1))/(2.0*N0)
      as_tg(2) = (ni(j0)+ni(j0+1))/(2.0*N0)
      as_tg(3) = (nz(j0)+nz(j0+1))/(2.0*N0)
      taus_tg(1) = (te(j0)+te(j0+1))/(2.0*T0)
      taus_tg(2) = (ti(j0)+ti(j0+1))/(2.0*T0)
      taus_tg(3) = (tz(j0)+tz(j0+1))/(2.0*T0)
      zeff_tg = (as_tg(2)*zs_tg(2)**2+as_tg(3)*zs_tg(3)**2)/(as_tg(1)*zs_tg(1)**2)  ! need to include all ions, fast ions too
      debye_tg = SQRT(k0*T0/(4.0*pi*N0*e0**2))     !Debye length unit in cm
      debye_tg = debye_tg/rhos0   ! TGLF normalized Debye length unit
!  lnlamda and taue from NRL formulary 
!  note that the local electron temperature and density need to be used for xnuei_tg
      lnlamda = 24.0 -0.5*LOG(as_tg(1)*N0)+LOG(taus_tg(1)*T0)       
      taue = (3.44E5)*((taus_tg(1)*T0)**1.5)/(as_tg(1)*N0*lnlamda)  !  sec
!  xnuei = 3/4 (Pi**0.5)/taue
      xnuei_tg = 0.75*SQRT(pi)/taue          ! 1/sec 
      xnuei_tg = xnuei_tg*a0/cs0             ! normalized electron-ion collision rate
      betae_tg = (8.0*pi*k0*N0*T0/Bunit**2)  ! electron beta unit
!
! local field derivatives
!
      rlns_tg(1) = -(ne(j0+1)-ne(j0))/(dr*as_tg(1)*N0)
      rlns_tg(2) = -(ni(j0+1)-ni(j0))/(dr*as_tg(2)*N0)
      if(as_tg(3).eq.0.0)then
        rlns_tg(3) = 0.0
      else
        rlns_tg(3) = -(nz(j0+1)-nz(j0))/(dr*as_tg(3)*N0)
      endif
      rlts_tg(1) = -(te(j0+1)-te(j0))/(dr*taus_tg(1)*T0)
      rlts_tg(2) = -(ti(j0+1)-ti(j0))/(dr*taus_tg(2)*T0)
      rlts_tg(3) = -(te(j0+1)-te(j0))/(dr*taus_tg(3)*T0)
!
      cexb =(rminor(j0+1)+rminor(j0))/(rmajor(1)*(q(j0+1)+q(j0)))      ! r/(q R(0))
      gamma_e_tg = cexb*(vexb(j0+1)-vexb(j0))/(dr*cs0) ! Waltz-Miller definition      
      gamma_p_tg = (vpar(j0+1)-vpar(j0))/(dr*cs0)
!
! local magnetic geometry
!
      rmin_tg = (rminor(j0+1)+rminor(j0))/(2.0*a0)
      rmaj_tg = (rmajor(j0+1)+rmajor(j0))/(2.0*a0)
      q_tg = (q(j0+1)+q(j0))/2.0
!
      if(igeo_tg.eq.0)then  ! s-alpha geometry
        shat_tg = (rmin_tg/q_tg)*(q(j0+1)-q(j0))/dr        ! r/q dq/dr 
        alpha_tg = -(8.0*pi*k0/Bunit**2)*q_tg**2*rmaj_tg*(ptot(j0+1)-ptot(j0))/dr
        xwell_tg = 0.0  ! magnetic well parameter (not used)
        theta0_tg = 0.0 ! radial wavenumber normalized to ky*shat (the most unstable modes is typically at kx=0)
     elseif(igeo_tg.eq.1)then ! Miller geometry
       q_prime_tg = (q_tg/rmin_tg)*(q(j0+1)-q(j0))/dr            
       p_prime_tg = (k0/Bunit**2)*(q_tg/rmin_tg)*(ptot(j0+1)-ptot(j0))/dr 
       shift_tg = (rmajor(j0+1)-rmajor(j0))/(dr*a0)
       kappa_tg = (elongation(j0+1)+elongation(j0))/2.0
       s_kappa_tg = rmin_tg*(elongation(j0+1)-elongation(j0))/dr
       delta_tg = (triangularity(j0+1)+triangularity(j0))/2.0
       s_delta_tg = (rmin_tg/SQRT(1.0-delta_tg**2))*(triangularity(j0+1)-triangularity(j0))/dr
     else  ! more options to come later
        write(*,*)"igeo_tg invalid",igeo_tg
        stop
     endif
!
! 
! control logic to save computation 
!
!     new_eikonal_tg = .TRUE.  will cause TGLF to compute the eikonal and save it.
!     new_eikonal_tg = .FALSE. will cause TGLF to skip the eikonal calulation and reuse the saved result. It also
!     skips the gaussian width search using the previously saved widths. This option 
!     can only be used if a call to TGLF has been made with new_eikonal=.TRUE.
!     new_eikonal_tg=.FALSE. is useful if this call to TGLF is a perturbation of a previous call with 
!     new_eikonal_tg=.TRUE.
!     For new_eikonal_tg=.FALSE. is is best to check that v_bar_sum_tg > 0 otherwise the first call found no
!     unstable modes and the pertubation call should just be skipped and the fluxes set to zero. 
!     it is ok to make several perturbative calls in a row (e.g. to vary temperature and density gradients) 
!     but you have to ! reset new_eikonal=.FASLE. since it is reset to .TRUE. by tglf_TM_interface to prevent mistakes.
!
      new_eikonal_tg = .TRUE. 
!
! call the TGLF model to compute fluxes
!
     CALL tglf_TM_interface
!
! debug 
     if(igeo_tg.eq.0)write(*,*)"shat_tg=",shat_tg,"alpha_tg=",alpha_tg
     if(igeo_tg.eq.1)then
      write(*,*)"q_prime_tg=",q_prime_tg,"p_prime_tg=",p_prime_tg
      write(*,*)"shift_tg=",shift_tg,"kappa_tg=",kappa_tg
      write(*,*)"s_kappa_tg=",s_kappa_tg,"delta_tg=",delta_tg
      write(*,*)"s_delta_tg=",s_delta_tg
     endif
     write(*,*)"q_tg =",q_tg
     write(*,*)"rmin_tg=",rmin_tg,"rmaj_tg=",rmaj_tg
     write(*,*)"mass_tg(1)=",(mass_tg(j),j=1,nspecies_tg)
     write(*,*)"as_tg(1)=",(as_tg(j),j=1,nspecies_tg)
     write(*,*)"taus_tg(1)=",(taus_tg(j),j=1,nspecies_tg)
     write(*,*)"rlns_tg=",(rlns_tg(j),j=1,nspecies_tg)
     write(*,*)"rlts_tg=",(rlts_tg(j),j=1,nspecies_tg)
     write(*,*)"betae_tg=", betae_tg," debye_tg=",debye_tg
     write(*,*)"xnuei_tg=",xnuei_tg
     write(*,*)"gamma_e_tg=",gamma_e_tg,"gamma_p_tg=",gamma_p_tg
     CALL put_kys(xky0_tg/sqrt(taus_tg(1)*mass_tg(2)))
     CALL put_gaussian_width(width_max_tg,width_min_tg,nwidth_tg &
      ,find_width_tg)
     CALL tglf
     write(*,*)"gamma =",get_growthrate(1)*cs0/a0
     write(*,*) 'freq_tg(1)  = ',get_frequency(1)*cs0/a0
     rhostar2=(rhos0/a0)**2
     write(*,*) 'v_bar(1) = ',get_v_bar(1)*rhostar2
     write(*,*) 'QL_particle_flux_tg(1,1,1) = ',get_QL_phi(1)*get_QL_particle_flux(1,1,1)*N0*cs0
     write(*,*) 'QL_energy_flux_tg(1,1,1) = ',get_QL_phi(1)*get_QL_energy_flux(1,1,1)*N0*T0*cs0
!
!
!     restore the units used in the transport code (un-normalize)
!
! note that an extra factor of drhodr is needed for transport codes using the rho-grid so that 
! the divergence of the flux density Qrho is 1/(dV/drho)d/drho(Qrho dV/drho) with V=volume of the flux surface. 
! Because of the use of r=rminor for derivatives, TGLF naturally gives flux densities Qr on the r-grid where the divergence is
! 1/(dV/dr)d/dr(Qr dV/dr) so Qrho = Qr drhodr.
!
      rhostar2=(rhos0/a0)**2
      do i=1,nspecies_tg
        do j=1,nfields_tg
          particle_flux_tg(i,j) = drhodr*N0*cs0*rhostar2*particle_flux_tg(i,j)
          energy_flux_tg(i,j) =drhodr* N0*T0*cs0*rhostar2*energy_flux_tg(i,j)
          write(*,*)"species = ",i," field = ",j
          write(*,*)"particle_flux = ",particle_flux_tg(i,j)," cm^-2 sec^-1"
          write(*,*)"energy_flux =",energy_flux_tg(i,j)," ev cm^-2 sec -1"
        enddo
      enddo    
!    
      STOP
!
      END  !PROGRAM tglf_TM_driver


