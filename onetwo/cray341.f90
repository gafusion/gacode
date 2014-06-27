MODULE PELLET_MOD
!-------------------------------------------------------------------------------
!PELLET_MOD is an F90 module of routines that calculates the ablation of solid
!  hydrogenic or impurity pellets in a plasma
!
!References:
!
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routine:
!
!  PELLET   - solves for pellet ablation along a trajectory in a plasma
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE PRL_MOD

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  PELLET_KUT,		& !hydrogenic ablation rate for Kuteev model
  PELLET_KUT_IM,	 & !impurity ablation rate for Kuteev model
  PELLET_MAC,		& !hydrogenic ablation rate for Macaulay model
  PELLET_NGS,		& !hydrogenic  ablation rate for NGS model
  PELLET_NGS_ABL,	& !auxiliary ablation routine for NGS model
  PELLET_NGS_QE,	 & !electron heat flux at pellet surface for NGS model
  PELLET_NGS_QF,	 & !fast ion heat flux at pellet surface for NGS model
  PELLET_PARKS,	  & !hydrogenic ablation rate for Parks model
  PELLET_PARKSQ,	 & !hydrogenic ablation rate for Parks model arbitrary Q
  PELLET_PARKS_IM,	 & !impurity ablation rate for Parks model
  PELLET_RK4		 !4th order Runge Kutta routine for time stepping

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
INTEGER, PRIVATE :: &
  k_pel_pl,		& !ablation model [hydrogenic 0-5, impurity 10-11]
				 !=0 ref NGS MODEL, e distr, ellipt shield
				 !=1 Milora NGS model, single e energy, spher shield
				 !=2 NGPS model, e distr, 1mm neutral layer
				 !=3 Macaulay hydrogenic model
				 !=4 Kuteev hydrogenic model
				 !=5 Parks hydrogenic model
				 !=6 Parks hydrogenic model arbitrary Q
				 !=10 Parks impurity model
				 !=11 Kuteev impurity model
  neg_pl			 !number of electron energy groups in plasma [-]

REAL(KIND=rspec), PRIVATE :: &
  amup_pl,		 & !atomic mass number of pellet atoms [-]
				 !used for either hydrogenic or impurity pellets
  denm_pl,		 & !molecular density in pellet
  ellipt_pl,		 & !ionized cloud ellipticity factor [-]
  fionc_pl,		& !ionized cloud thickness factor [-]
  fqe_pl,		  & !electron heat flux attenuation factor [-]
  qeo_pl,		  & !incident electron heat flux from plasma [keV/m**2/s]
  rcl_pl,		  & !neutral cloud radius [m]
  rhosrp0_pl,		& !pellet mass density * initial pellet radius [kg/m**2]
  rpel_pl,		 & !initial pellet radius [m]
  vpel_pl,		 & !pellet velocity [m/s]
  xrhoro_pl		  !last solution to normalized cloud thickness [-]

!Fast ions
LOGICAL, PRIVATE :: &
  l_fast_pl		  !option to include fast ions [logical]

INTEGER, PRIVATE :: &
  nf_pl,		   & !number of fast ion species (H,D,T,He3,He4 only) [-]
  nefmax_pl		  !max no. of energy groups for any fast ions [-]

INTEGER, PRIVATE, ALLOCATABLE :: &
  izf_pl(:),		 & !charge of fast ions (1 or 2 only) [-]
  nef_pl(:)		  !number of energy groups for fast ions [-]

REAL(KIND=rspec), PRIVATE, ALLOCATABLE :: &
  amuf_pl(:),		& !atomic mass number of fast ions [-]
  vcf_pl(:),		 & !local critical velocity for fast ions [m/s]
  ef_pl(:,:),		& !energy groups for fast ions [keV]
  denf_pl(:,:),	  & !local fast ion density by energy group [/m**3]
  efo_pl(:,:),	   & !local average fast ion energy in energy group [keV]
  qfo_pl(:,:)		!local energy flux in energy group [keV/m**2/s]

!Physical and conversion constants
REAL(KIND=rspec), PRIVATE, PARAMETER :: &
  z_coulomb=1.6022e-19,	& !coulomb charge [C]
  z_eion=0.0326,		 & !hydrogen ionization energy [keV]
  z_electronmass=9.1095e-31, & !electron mass [kg]
  z_epsilon0=8.8542e-12,	 & !permittivity of free space [F/m]
  z_evap=1.0e-5,		 & !solid hydrogen evaporation energy [keV]
  z_gam=1.4,			 & !ratio of hydrogen specific heats [-]
  z_j7kv=1.6022e-16,	   & !conversion constant [J/keV]
  z_navogadro=6.0221e23,	 & !Avogadro's number [-]
  z_pi=3.141592654,		& !pi [-]
  z_protonmass=1.6726e-27,   & !proton mass [kg]
  z_tolc=0.1,			& !tolerance for convergence of cloud solution [-]
  z_tolr=0.1,			& !tolerance for fractional radius of pellet remaining [-]
  z_tolt=0.2,			& !tolerance for fractional change in te per step [-]
  z_rclmin=0.001		   !minimum cloud radius [m]

!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE PELLET(kprint, nprt, ndat, &
			k_pel,amupel,rpel,vpel,n_r,dvol_r,den0_r,te0_r,n_p,map_p,s_p,&
			pden_r,iflag,message, &
			NF,NE_F,IZ_F,AMU_F,E_EF,VC_RF,DEN_REF, &
			PEL_IONS,T_P,RPEL1_P,SRC_P,DEN0_P,DEN1_P,TE0_P,TE1_P, &
			R0,A0,BT0,NCSOL,K_PRL,NPRLCLD,IPRLCLD,FPELPRL,PRLINJANG,PRLQ_R,PRLDEP)
!-------------------------------------------------------------------------------
!PELLET calculates the ablation profile for solid pellets injected into a plasma
!
!References:
!  W.A.Houlberg, S.L.Milora, S.E.Attenberger, Nucl Fusion 28 (1988) 595
!  S.L.Milora, ORNL/TM-8616 (1983)
!  W.A.Houlberg, M.A.Iskra, H.C.Howe, S.E.Attenberger, ORNL/TM-6549 (1979)
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!  L.R. Baylor, Added linkages to PRL model routine, Parks,Baylor, submitted 2004 
!-------------------------------------------------------------------------------

!Declaration of input variables
integer, intent(in) :: &
	kprint,		& 	! kprint-option for detail of printout.
				! = 0 errors in input and solution.
				! = 1 above plus input values and options.
				! = 2 above plus radial input and calculated parameters.
				! = 3 above plus parameters along pellet path.
	nprt,		&	! file ID for output
	ndat			! peldat ID for output
INTEGER, INTENT(IN) :: &
  k_pel,		   & !ablation model [hydrogenic 0-5, impurity 10-11]
				 !=0 ref NGS MODEL, e distr, ellipt shield
				 !=1 Milora NGS model, single e energy, spher shield
				 !=2 NGPS model, e distr, 1mm neutral layer
				 !=3 Macaulay hydrogenic model
				 !=4 Kuteev hydrogenic model
				 !=5 Parks hydrogenic model
				 !=6 Parks hydrogenic model arbitrary Q
				 !=10 Parks impurity model
				 !=11 Kuteev impurity model
				 !=else failure
  n_r,			 & !number of radial plasma cells [-]
  n_p			  !number of segments along pellet path [-]

INTEGER, INTENT(IN) :: &
  map_p(:)		   !plasma cell for each path segment [-]

REAL(KIND=rspec), INTENT(IN) :: &
  amupel,		  & !atomic mass number of pellet atoms [-]
				 !used for either hydrogenic or impurity pellets
  rpel,			& !initial pellet spherical radius [m]
  vpel			 !pellet velocity [m/s]

REAL(KIND=rspec), INTENT(IN) :: &
  dvol_r(:),		 & !volume of plasma cell i [m**3]
  den0_r(:),		 & !initial electron density in radial cells [/m**3]
  te0_r(:),		& !initial electron temperature in radial cells [keV]
  s_p(:)			 !distance to the beginning of path segment [m]
				 !s_p(1) is the entry point and can be non-zero

!Declaration of optional input variables
INTEGER, INTENT(IN), OPTIONAL :: &
  NF,			& !number of fast ion species (H,D,T,He3,He4) [-]
  NE_F(:),		 & !number of energy groups for fast ions [-]
  IZ_F(:),		 & !charge of fast ions (1,2) [-]
  NCSOL,		   & !Number of SOL grid points in n_r   [-]
  K_PRL,		   & !flag for PRL model calculation
  NPRLCLD,		 & !Number of PRL cloudlets from pellet [-]
  IPRLCLD			!ID number of PRL cloudlet for diag output [-]

REAL(KIND=rspec), INTENT(IN), OPTIONAL :: &
  AMU_F(:),		& !atomic mass number of fast ions [-]
  E_EF(:,:),		 & !energy group boundaries for fast ions [keV]
				 !E_EF(1,NF)= the birth energy
				 !E_EF(NE_F+1,NF)= energy where becomes thermal
				 !E_EF(j,NF) > E_EF(j+1,NF)
  VC_RF(:,:),		& !critical velocity profile for fast ions [m/s]
  DEN_REF(:,:,:),	& !density profile of fast ion per energy group [1/m**3]
  R0,			& !Major radius [m]
  A0,			& !Minor radius [m]
  BT0,			 & !Magnetic field on axis used by PRL [T]
  FPELPRL,		 & !PRL fraction of pellet mass to drift [-]
  PRLINJANG,		 & !PRL injection angle (LFS=0, HFS=pi) [radians]
  PRLQ_R(:)		  !q profile for PRL model [-]


!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message			!warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag			!error and warning flag [-]
				 !=-1 warning
				 !=0 no warnings or errors
				 !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  pden_r(:)		  !increase in plasma density in radial cells [/m**3]

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  PEL_IONS		   !number of ions in pellet [-]

REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  T_P(:),		  & !time of pellet entry to path cells [s]
  RPEL1_P(:),		& !pellet effective radius at exit from path cells [m]
  SRC_P(:),		& !source rate in path cells [/s]
  DEN0_P(:),		 & !electron density at entrance to path cells [keV]
  DEN1_P(:),		 & !electron density at exit from path cells [keV]
  TE0_P(:),		& !electron temperature at entrance to path cells [keV]
  TE1_P(:),		& !electron temperature at exit from path cells [keV]
  PRLDEP(:)		  !resulting density deposition from PRL [/m**3]

!-------------------------------------------------------------------------------
!Declartation of local variables
LOGICAL :: &
  l_inout,l_inside

INTEGER :: &
  i,l, ii, idiag

REAL(KIND=rspec) :: &
  dennew,denold,dt,rp,rpold,src,srcp,t,tenew,teold


!
! Local variables for PRL extension
!
REAL(KIND=rspec) :: &
  fi,nfloat,neasum,rho_p,rpellet,rmaj,amin, &
  pnum,fnumcld,fcld,fnumpart,nump,fncp,vpelprl

REAL(KIND=rspec), ALLOCATABLE :: &
  nea(:),		 & !density profile after cloudlet [cm^-3]
  tea(:),		 & !temperature profile after cloudlet [eV]
  neb(:),		 & !density profile before cloudlet [cm^-3]
  teb(:),		 & !temperature profile before cloudlet [eV]
  rgrid(:)		!rho grid for PRL [-]

INTEGER :: icld, ngrid

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Set internal data from external input
k_pel_pl=k_pel
amup_pl=amupel
rpel_pl=rpel
vpel_pl=vpel
!Hydrogenic molecular density
denm_pl=-8.6857e26*amup_pl**2+6.3023e27*amup_pl+2.1200e28

!Internal data differing from default
IF(k_pel_pl == 0) THEN

  !hydrogenic NGS model: e distribution, elliptical neutral shield
  neg_pl=10
  ellipt_pl=15.0
  fionc_pl=0
  rcl_pl=rpel_pl+z_rclmin
	  
ELSEIF(k_pel_pl == 1) THEN

  !Milora hydrogenic NGS model: single e energy, spherical neutral shield
  neg_pl=1
  ellipt_pl=1
  fionc_pl=0

ELSEIF(k_pel_pl == 2) THEN

  !hydrogenic NGPS model: e distribution, 1mm neutral layer thickness
  neg_pl=10
  ellipt_pl=1
  fionc_pl=1
  rcl_pl=rpel_pl+z_rclmin

ELSEIF(k_pel_pl == 3) THEN

  !Macaulay hydrogenic NGS
  neg_pl=1
  ellipt_pl=1
  fionc_pl=0

ELSEIF(k_pel_pl == 4) THEN

  !Kuteev hydrogenic NGS
  neg_pl=1
  ellipt_pl=1
  fionc_pl=0

ELSEIF(k_pel_pl == 5) THEN

  !Parks hydrogenic NGS
  neg_pl=1
  ellipt_pl=1
  fionc_pl=0

ELSEIF(k_pel_pl == 6) THEN

  !Parks hydrogenic NGS-Q
  neg_pl=1
  ellipt_pl=1
  fionc_pl=0

ELSEIF(k_pel_pl == 10) THEN

  !Parks imopurity model
  neg_pl=1
  ellipt_pl=1
  fionc_pl=0

ELSEIF(k_pel_pl == 11) THEN

  !Kuteev impurity model
  neg_pl=1
  ellipt_pl=1
  fionc_pl=0

ELSE

  iflag=1
  message='PELLET/ERROR(1):invalid model choice'
  GOTO 9999

ENDIF
!Fast ions
IF(PRESENT(NF) .AND. &
   NF > 0 .AND. &
   PRESENT(NE_F) .AND. &
   PRESENT(IZ_F) .AND. &
   PRESENT(AMU_F) .AND. &
   PRESENT(E_EF) .AND. &
   PRESENT(VC_RF) .AND. &
   PRESENT(DEN_REF)) THEN

  l_fast_pl=.TRUE.
  nf_pl=NF
  nefmax_pl=0

  DO i=1,nf_pl

	IF(NE_F(i) > nefmax_pl) nefmax_pl=NE_F(i)

  ENDDO

  ALLOCATE(nef_pl(1:nf_pl), &
		 izf_pl(1:nf_pl), &
		 amuf_pl(1:nf_pl), &
		 ef_pl(1:nefmax_pl+1,1:nf_pl), &
		 vcf_pl(1:nf_pl), &
		 denf_pl(1:nefmax_pl,1:nf_pl), &
		 efo_pl(1:nefmax_pl,1:nf_pl), &
		 qfo_pl(1:nefmax_pl,1:nf_pl))

  nef_pl(1:nf_pl)=NE_F(1:nf_pl)
  izf_pl(1:nf_pl)=IZ_F(1:nf_pl)
  amuf_pl(1:nf_pl)=AMU_F(1:nf_pl)
  ef_pl(1:nefmax_pl+1,1:nf_pl)=E_EF(1:nefmax_pl+1,1:nf_pl)
  vcf_pl(:)=0
  denf_pl(:,:)=0
  efo_pl(:,:)=0
  qfo_pl(:,:)=0

ELSE

  l_fast_pl=.FALSE.

ENDIF
!Zero output arrays
pden_r(:)=0
IF(PRESENT(T_P)) T_P(:)=0
IF(PRESENT(RPEL1_P)) RPEL1_P(:)=0.
IF(PRESENT(SRC_P)) SRC_P(:)=0
IF(PRESENT(DEN0_P)) DEN0_P(:)=0
IF(PRESENT(DEN1_P)) DEN1_P(:)=0
IF(PRESENT(TE0_P)) TE0_P(:)=0
IF(PRESENT(TE1_P)) TE1_P(:)=0

!Set private pellet parameters for normalizations
rhosrp0_pl=(2*amup_pl*z_protonmass*denm_pl)*rpel_pl

!Initialize parameters at beginning of pellet trajectory
t=0
rp=rpel_pl
xrhoro_pl=0


IF(PRESENT(NPRLCLD)) THEN
   ngrid = N_R - NCSOL
   ALLOCATE(neb(1:ngrid), &
		 teb(1:ngrid), &
		 nea(1:ngrid), &
		 tea(1:ngrid), &
		 rgrid(1:ngrid))

   icld = 0  ! PRL cloudlet counter
   fcld = 0.0
ENDIF

!Optional output
IF(PRESENT(PEL_IONS)) PEL_IONS=4*z_pi/3*rpel_pl**3*(2*denm_pl)

!-------------------------------------------------------------------------------
!Follow the pellet path and determine the ablation rate
!-------------------------------------------------------------------------------
!Flag to indicate whether pellet has entered plasma
l_inside=.FALSE.

!Flag to indicate whether pellet has entered and exited plasma
l_inout=.FALSE.
print*, "n_p = ", n_p
DO l=1,n_p
	print*, "l = ", l
  IF(rp > 1.0e-6*rpel_pl) THEN

	i=map_p(l)

	IF(((i <= 0) .OR. (i > n_r)) .AND. &
	 (.NOT. l_inside)) THEN

	!Pellet has not entered plasma yet
	!Increment parameters for path outside plasma
	dt=(s_p(l+1)-s_p(l))/vpel_pl
	t=t+dt
	IF(PRESENT(T_P)) T_P(l)=t

	ELSEIF((i > 0 .AND. i <= n_r) .AND. &
		 (.NOT. l_inout)) THEN

	!Pellet is inside plasma and never exited
	l_inside=.TRUE.
	denold=den0_r(i)+pden_r(i)
	teold=(den0_r(i)*te0_r(i)-pden_r(i)*z_eion/1.5)/denold
	IF(teold < z_eion) teold=z_eion
	IF(l_fast_pl) THEN

	  vcf_pl(1:nf_pl)=VC_RF(i,1:nf_pl)
	  denf_pl(1:nefmax_pl,1:nf_pl)=DEN_REF(i,1:nefmax_pl,1:nf_pl)

	ENDIF

	dt=(s_p(l+1)-s_p(l))/vpel_pl
	rpold=rp
	print*, "te into rk4 = ", teold
	CALL PELLET_RK4(dvol_r(i),dt, &
				t,rp,teold,denold, &
				srcp,iflag,message)
	pden_r(i)=pden_r(i)+srcp
	dennew=den0_r(i)+pden_r(i)
	tenew=(den0_r(i)*te0_r(i)-pden_r(i)*z_eion/1.5)/dennew
	src = srcp * dvol_r(i)
	IF(tenew < z_eion) tenew=z_eion
!	if (kprint .gt. 2) then
!			write(ndat, '(5(E16.8))') s_p(l) * 1e2, src, srcp * 1e-6, qeo_pl * 1e-1, fqe_pl
!	endif
		
	!Set optional output parameters along path
	IF(PRESENT(T_P)) T_P(l)=t
	IF(PRESENT(RPEL1_P)) RPEL1_P(l)=rp
	IF(PRESENT(SRC_P)) SRC_P(l)=srcp*dvol_r(i)/dt
	IF(PRESENT(DEN0_P)) DEN0_P(l)=denold
	IF(PRESENT(DEN1_P)) DEN1_P(l)=dennew
	IF(PRESENT(TE0_P)) TE0_P(l)=teold
	IF(PRESENT(TE1_P)) TE1_P(l)=tenew

	IF (K_PRL .gt.0 .and. PRESENT(NPRLCLD) .and. rp > 0.0) THEN   ! Call PRL drift model
!
! Check to see if cloudlet mass has been ablated away since last cloudlet
!
	  fnumcld = NPRLCLD

	  !Check for remaining pellet size to do drift calc
	  if (rp**3 .lt. ((1.0-(fcld/fnumcld))*rpel**3)) then

		fcld = fcld + 1.0  ! Increment cloud counters
		icld = icld + 1
		! Do PRL calculation
		rgrid(1) = 0.0
		do ii=1,ngrid !Profile generation
		fi = float(ii)
		nfloat = float(ngrid)
		rgrid(ii) = fi/nfloat
		teb(ii) = te0_r(ii)*1000.0   ! convert to eV
		neb(ii) = den0_r(ii)/10**6   ! convert to cm^-3
		nea(ii) = 0.0
		tea(ii) = 0.0
		enddo !Profile generation
		rgrid(ngrid) = 1.0
		fi = l
		fncp = float(ngrid)
		!rho_p = 1.0-fi/nfloat
		!rho_p = lc(l)/fncp
		rho_p = map_p(l)/nfloat
		if (rho_p .gt. 1.0) rho_p = 1.0
		rpellet = rp*100 ! cm from m
		rmaj = R0*100	! cm from m
		amin = A0*100	! cm from m
		pnum = PEL_IONS/nprlcld
		vpelprl = vpel*-100.0 ! cm/s from m/s

		  if (iprlcld == icld) then	! diagnostic output flag
			idiag = 1
		  else
			idiag = 0
		  endif
		call prl(icld,idiag,BT0,rmaj,amin,PRLINJANG,rpellet,  &
				 rho_p, pnum, vpelprl, ngrid, teb, neb,  &
				 PRLQ_R, tea, nea, rgrid, nump)

		!Save PRL deposition profile and adjust Te profile
		neasum=0

		do i=1,ngrid
		neasum = neasum+nea(i)
		enddo

		fnumpart = (4.0/3.0)*z_pi*(rpel**3)*2*denm_pl/nprlcld !Num pellet particles

		if(neasum > 0.0) then

		do i=1,ngrid  
		  !nea(i) = (fnumpart/neasum)*nea(i)/dvol(i)
		  nea(i) = nea(i)/dvol_r(i) !Convert to density
		enddo

		endif

		do i=1,ngrid  ! Determine dep profile

		PRLDEP(i) = PRLDEP(i) + nea(i)
		! Adiabatic approximation
		!te(i) = te(i)*den(i)/(den(i) + nea(i))
		!den(i) = den(i) + nea(i)

		enddo

	  endif

	endif

	ELSEIF(((i <= 0) .OR. (i > n_r)) .AND. l_inside) THEN

	!Pellet has been in the plasma but is now outside
	l_inout=.TRUE.

	ENDIF !present
  ENDIF  !rp>1e-6

ENDDO

!
! PRL deposition profile
!
if (PRESENT(NPRLCLD) .and. fcld > 1 ) then ! Done with PRL - determine prldep
   do i=1,ngrid  ! Determine dep profile
		PRLDEP(i) = FPELPRL*PRLDEP(i)*fnumcld/fcld + &
					   (1-fpelprl)*pden_r(i)
   enddo

   write(*,*) ' Fueling efficiency: ',neasum/(fnumpart*nprlcld)

!	  if (kprl.gt.0) then  ! New output of PRL dep calculation
!		 write(nerr,*) ' '
!		 write(nerr,*) ' PRL output for num clouds = ', icld
!		 do i=1,nc
!		 dennew=den(i)+prldep(i)
!		   write(nerr,1230) i,den(i),dennew,dvol(i),prldep(i)
!		 enddo
!	  endif
endif

if (kprint .gt. 2) then
	write(ndat, *) '\n  init. Te (eV)   init. de (/cc)  dp inc. (/cc)'
	do i = 1, n_r
		write(ndat, '(3(e16.8))') te0_r(i) * 1e3, den0_r(i) * 1e-6, pden_r(i) * 1e-6
	end do
endif

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET

SUBROUTINE PELLET_KUT(rp,te,den, &
				rdot)
!-------------------------------------------------------------------------------
!PELLET_KUT determines the pellet ablation rate using a fit generated by Kuteev
!
!References:
!  B.V.Kuteev, Nucl Fusion 35 (1995) 431
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den			  !electron density [/m**3]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  rdot			 !rate of change in pellet radius [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec), PARAMETER :: &
  temin=1.0e-3

REAL(KIND=rspec) :: &
  dndt,dinf,rpcm,tinf

!-------------------------------------------------------------------------------
!Check pellet size and minimum electron temperature
!-------------------------------------------------------------------------------
IF((rp < (z_tolr*rpel_pl)) .OR. &
   (te < temin)) THEN

  rdot=0
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Calculate ablation rate
!-------------------------------------------------------------------------------
!Convert to Kuteev's units of cm and eV
tinf=te*1.0e3
rpcm=rp*1.0e2
dinf=den*1.0e-6

!Compute the ablation rate in atoms/s, Eq 72
dndt=3.46e14*tinf**(1.72)*dinf**(0.453)*rpcm**(1.443)*amup_pl**(-0.283)

!Compute dr/dt in m/s
rdot=-dndt/((2*denm_pl)*4*z_pi*rp**2)

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_KUT

SUBROUTINE PELLET_KUT_IM(rp,te,den, &
				 rdot)
!-------------------------------------------------------------------------------
!PELLET_KUT_IM determines the pellet impurity ablation rate using a fit
!  generated by Kuteev, et al
!
!References:
!  B.V.Kuteev, V.Yu.Sergeev, S.Sudo, Nucl Fusion 35 (1995) 1167
!  B.V.Kuteev, Y.Yu.Sergeev, A.Yu.Kostrukov, V.A.Segal, O.A.Bakhareva,
!	P.B.Parks, Bull Am Phys Soc 44 (1999) 115
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The Kr coefficient given in the APS poster is used (0.94e15) instead of the
!	earlier published value (2.15e15)
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den			  !electron density [/m**3]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  rdot			 !rate of change in pellet radius [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec), PARAMETER :: &
  temin=1.0e-3

REAL(KIND=rspec) :: &
  coeff,		   & !
  dena,			& !atomic density [atoms/m**3]
  dndt,			& !ablation rate [atoms/s]
  dinf,			& !e density at infinity [/cm**3]
  rhomass,		 & !solid mass density [g/m**3]
  rpcm,			& !pellet radius [cm]
  tinf			 !e temp at infinity [eV]

!-------------------------------------------------------------------------------
!Check pellet size and minimum electron temperature
!-------------------------------------------------------------------------------
IF((rp < (z_tolr*rpel_pl)) .OR. &
   (te < temin)) THEN

  rdot=0
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Calculate ablation rate
!-------------------------------------------------------------------------------
!Convert to Kuteev's units of cm and eV
tinf=te*1.0e3
rpcm=rp*1.0e2
dinf=den*1.0e-6

!Determine coefficient to use (function of impurity)
IF(amup_pl >=6.0 .AND. &
   amup_pl <=7.0) THEN

  !Lithium (nominally 6.9)
  coeff=1.04e15
  rhomass=0.534e6

ELSEIF(amup_pl >= 8.5 .AND. &
	 amup_pl <= 9.5) THEN

  !Beryllium (nominally 9.0)
  coeff=0.91e15
  rhomass=1.848e6

ELSEIF(amup_pl >= 10.0 .AND. &
	 amup_pl <= 11.0) THEN

  !Boron (nominally 10.8)
  coeff=0.67e15
  rhomass=2.340e6

ELSEIF(amup_pl >= 11.5 .AND. &
	 amup_pl <= 12.5) THEN

  !Carbon (nominally 12.0)
  coeff=0.52e15
  rhomass=2.000e6

ELSEIF(amup_pl >= 20.0 .AND. &
	 amup_pl <= 21.0) THEN

  !Neon (nominally 20.2)
  coeff=2.32e15
  rhomass=2.205e6

ELSEIF(amup_pl >= 39.0 .AND. &
	 amup_pl <= 40.0) THEN

  !Argon (nominally 39.9)
  coeff=1.38e15
  rhomass=1.400e6

ELSEIF(amup_pl >= 47.0 .AND. &
	 amup_pl <= 48.0) THEN

  !Titanium (nominally 47.9)
  coeff=0.49e15
  rhomass=4.505e6

ELSEIF(amup_pl >= 83.0 .AND. &
	 amup_pl <= 84.0) THEN

  !Krypton (nominally 83.8)
  coeff=0.94e15
  rhomass=2.155e6

ELSEIF(amup_pl >= 95.0 .AND. &
	 amup_pl <= 96.0) THEN

  !Molybdenum (nominally 95.9)
  coeff=0.35e15
  rhomass=10.220e6

ELSEIF(amup_pl >= 130.5 .AND. &
	 amup_pl <= 131.5) THEN

  !Xenon (nominally 131.0)
  coeff=0.72e15
  rhomass=3.520e6

ELSEIF(amup_pl >= 183.5 .AND. &
	 amup_pl <= 184.5) THEN

  !Tungsten (nominally 184.0)
  coeff=0.26e15
  rhomass=19.350e6

ELSE

  coeff=0
  rhomass=1.0e6

ENDIF

dena=rhomass*z_navogadro/amup_pl
	
!Compute the ablation rate, in atoms/s
dndt=coeff*tinf**(1.64)*dinf**(0.333)*rpcm**(1.333)*amup_pl**(-0.333)

!Compute dr/dt in m/s
rdot=-dndt/(dena*4*z_pi*rp**2)

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_KUT_IM

SUBROUTINE PELLET_MAC(rp,te,den, &
				rdot)
!-------------------------------------------------------------------------------
!PELLET_MAC determines the pellet ablation rate using a fit to a 2D simulation
!  by Macaulay
!
!References:
!  A.K.Macaulay, Nucl Fusion 34 (1994) 43
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den			  !electron density [/m**3]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  rdot			 !rate of change in pellet radius [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec), PARAMETER :: &
  temin=1.0e-3

REAL(KIND=rspec) :: &
  corr,dinf,dnstar,fdens,fmass,g2dgs,rpcm,testar,tinf

!-------------------------------------------------------------------------------
!Check pellet size and minimum electron temperature
!-------------------------------------------------------------------------------
IF((rp < (z_tolr*rpel_pl)) .OR. &
   (te < temin)) THEN

  rdot=0
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Calculate ablation rate
!-------------------------------------------------------------------------------
!Convert from deuterium to arbitrary hydrogenic species
fdens=denm_pl/(-8.6857e26_rspec*4.0_rspec  &
		   +6.3023e27_rspec*2.0_rspec  &
		   +2.1200e28_rspec)
fmass=(2/amup_pl)**(1.0/3.0)

!Convert to Macaulay's units of cm and eV
tinf=te*1.0e3
rpcm=rp*1.0e2
dinf=den*1.0e-6

!Compute the ablation rate, g2dgs, in atoms/s
dnstar=1.6e11*tinf**2/rpcm
testar=1.1e-7*((dinf*rpcm)**2/tinf)**(1.0/3.0)
corr=1.0+0.08*LOG(dnstar*SQRT(testar)/6.0e18)
g2dgs=9.0e15*corr*(1.0+0.09*LOG(250.0*testar))*tinf**(11.0/6.0) &
	 *(rpcm**4*dinf)**(1.0/3.0)/(LOG(0.147*tinf))**(2.0/3.0)

!Compute dr/dt in m/s
rdot=-g2dgs/((2*denm_pl)*4*z_pi*rp**2)

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_MAC

SUBROUTINE PELLET_NGS(l_newcell,rp,te,den,tfd, &
				rdot,iflag,message)
!-------------------------------------------------------------------------------
!PELLET_NGS determines the pellet ablation rate by iterating on the
!  dimensionless cloud thickness
!
!References:
!  W.A.Houlberg, S.L.Milora, S.E.Attenberger, Nucl Fusion 28 (1988) 595
!  S.L.Milora, ORNL/TM-8616 (1983)
!  Forsythe, Malcolm, Moler, Comp Meth for Math Computations, Prentice-Hall
!	(1977) 161
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The zeroin procedure from Forsythe, et al., is used to find the simultaneous
!	solution to two equations for dr/dt, which are non-linear functions of the
!	cloud thickness and plasma parameters
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den,			 & !electron density [/m**3]
  tfd			  !collisionless self-limiting parameter =[-]

!Declaration of in/out variables
LOGICAL, INTENT(INOUT) :: &
  l_newcell		  !flag to recalculate fast ion fluxes in plasma
				 !=TRUE new plasma cell
				 !=FALSE same plasma cell as previous call

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message			!warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag			!error and warning flag [-]
				 !=-1 warning
				 !=0 no warnings or errors
				 !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  rdot			 !rate of change in pellet radius [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec), PARAMETER :: &
  amin=1.0e-10,		  & !lower limit of norm cloud thickness range [-]
  bmax=0.1,			& !upper limit of norm cloud thickness range [-]
  temin=1.0e-3,		  & !cutoff electron temperature for ablation [keV]
  xrange=3.0			 !check xrhoro/xrange < b < xrhoro*xrange for soln

REAL(KIND=rspec) :: &
  a,b,			 & !normalized cloud thicknesses spanning soln
  c,d,difa,difb,difc,e,fm,p,q,r,rdot1a,rdot1b,rdot2a,rdot2b,rdot2c,s, &
  toltst,xrhor

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Check pellet size and minimum electron temperature
IF((rp < (z_tolr*rpel_pl)) .OR. &
   (te < temin)) THEN

  rdot=0
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Begin zeroin procedure
!-------------------------------------------------------------------------------
IF((xrhoro_pl < amin) .OR. &
   (xrhoro_pl > bmax)) THEN

!Use entire range
  a=amin
  b=bmax
  CALL PELLET_NGS_ABL(l_newcell,rp,te,den,a,tfd, &
				rdot1a,rdot2a,difa)
  l_newcell=.FALSE.
  CALL PELLET_NGS_ABL(l_newcell,rp,te,den,b,tfd, &
				rdot1b,rdot2b,difb)

  IF(ABS(SIGN(1.0_rspec,difa)+SIGN(1.0_rspec,difb)) > 1.0e-5) THEN

	iflag=1
	message='PELLET_NGS/ERROR(1):solution out of range'
	GOTO 9999

  ENDIF

ELSE

  !Evaluate at previous solution
  b=xrhoro_pl
  CALL PELLET_NGS_ABL(l_newcell,rp,te,den,b,tfd, &
				rdot1b,rdot2b,difb)
  l_newcell=.FALSE.

  !Check for convergence of zeroin procedure at previous solution
  IF(ABS(difb) < z_tolc) THEN

	xrhor=b
!	print*, "xrhor here: ", xrhor
	xrhoro_pl=xrhor
	rdot=rdot2b
	GOTO 9999

  ENDIF

  IF(ABS(difb) > (1.99999)) THEN

	!Start over using entire range
	a=amin
	b=bmax
	CALL PELLET_NGS_ABL(l_newcell,rp,te,den,a,tfd, &
				rdot1a,rdot2a,difa)
	l_newcell=.FALSE.
	CALL PELLET_NGS_ABL(l_newcell,rp,te,den,b,tfd, &
				rdot1b,rdot2b,difb)


	IF(ABS(SIGN(1.0_rspec,difa)+SIGN(1.0_rspec,difb)) > 1.0e-5_rspec) THEN

	iflag=1
	message='PELLET_NGS/ERROR(2):solution out of range'
	GOTO 9999

	ENDIF

	GOTO 10

  ENDIF

  !Evaluate at previous solution*xrange
  a=b
  difa=difb
  rdot2a=rdot2b
  b=xrhoro_pl*xrange
  IF(b > bmax) b=bmax
  CALL PELLET_NGS_ABL(l_newcell,rp,te,den,b,tfd, &
				rdot1b,rdot2b,difb)

  IF(ABS(SIGN(1.0_rspec,difa)+SIGN(1.0_rspec,difb)) < 1.0e-5) GOTO 10

  !Check which side to try next
  IF(ABS(difb) > ABS(difa)) THEN

	!Evaluate at previous solution/xrange
	b=a
	difb=difa
	rdot2b=rdot2a
	a=xrhoro_pl/xrange
	IF(a < amin) a=amin
	CALL PELLET_NGS_ABL(l_newcell,rp,te,den,a,tfd, &
				rdot1a,rdot2a,difa)

	IF(ABS(SIGN(1.0_rspec,difa)+SIGN(1.0_rspec,difb)) < 1.0e-5) GOTO 10

	!Evaluate at amin
	b=a
	difb=difa
	rdot2b=rdot2a
	a=amin
	CALL PELLET_NGS_ABL(l_newcell,rp,te,den,a,tfd, &
				rdot1a,rdot2a,difa)

	IF(ABS(SIGN(1.0_rspec,difa)+SIGN(1.0_rspec,difb)) > 1.0e-5) THEN

	iflag=1
	message='PELLET_NGS/ERROR(3):solution out of range'
	GOTO 9999

	ENDIF

  ELSE

	!Evaluate at bmax
	a=b
	difa=difb
	rdot2a=rdot2b
	b=bmax
	CALL PELLET_NGS_ABL(l_newcell,rp,te,den,b,tfd, &
				rdot1b,rdot2b,difb)

	IF(ABS(SIGN(1.0_rspec,difa)+SIGN(1.0_rspec,difb)) > 1.0e-5) THEN

	iflag=1
	message='PELLET_NGS/ERROR(4):solution out of range'
	GOTO 9999

	ENDIF

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Step through interval
!-------------------------------------------------------------------------------
!Initial conditions
10 c=a
   difc=difa
   rdot2c=rdot2a
   d=b-a
   e=d

!Identify closest solution
20 IF(ABS(difc) < ABS(difb)) THEN

	 !Reorder points a, b, and c.
	 a=b
	 b=c
	 c=a
	 difa=difb
	 rdot2a=rdot2b
	 difb=difc
	 rdot2b=rdot2c
	 difc=difa
	 rdot2c=rdot2a

   ENDIF

!Convergence test
   toltst=2.0e-6*ABS(b)
   fm=(c-b)/2

   IF((ABS(fm) <= toltst) .AND. (ABS(difb) >= z_tolc)) THEN

	 iflag=1
	 message='PELLET_NGS/ERROR(5):no convergence'
	 GOTO 9999

   ENDIF

   !Check for end of zeroin procedure
   IF((ABS(fm) <= toltst) .OR. &
	(ABS(difb) < z_tolc)) THEN

	 xrhor=b
!	 print*, "xrhor here: ", xrhor
	 xrhoro_pl=xrhor
	 rdot=rdot2b
	 GOTO 9999

   ENDIF

   !Check if bisection is necessary
   IF(ABS(e) < toltst) GOTO 30
   IF(ABS(difa) <= ABS(difb)) GOTO 30

   !Interpolate
   IF(a == c) THEN

	 !Linear interpolation
	 s=difb/difa
	 p=2*fm*s
	 q=1.0-s

   ELSE

	 !Inverse quadratic interpolation.
	 q=difa/difc
	 r=difb/difc
	 s=difb/difa
	 p=s*(2*fm*q*(q-r)-(b-a)*(r-1.0))
	 q=(q-1.0)*(r-1.0)*(s-1.0)

   ENDIF

   !Adjust signs
   IF(p > 0.0) q=-q
   p=ABS(p)

   !Check if interpolation is acceptable
   IF((2*p) >= (3*fm*q-ABS(toltst*q))) GOTO 30
   IF(p >= ABS(e*q/2)) GOTO 30
   e=d
   d=p/q
   GOTO 40

   !Bisection
30 d=fm
   e=d

  !Complete step
40 a=b
   difa=difb
   rdot2a=rdot2b
   IF(ABS(d) > toltst) b=b+d
   IF(ABS(d) <= toltst) b=b+SIGN(toltst,fm)
   CALL PELLET_NGS_ABL(l_newcell,rp,te,den,b,tfd, &
				 rdot1b,rdot2b,difb)
   IF((difb*(difc/ABS(difc))) > 0.0) GOTO 10
   GOTO 20

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_NGS

SUBROUTINE PELLET_NGS_ABL(l_newcell,rp,te,den,xrhor,tfd, &
				  rdot1,rdot2,difrdt)
!-------------------------------------------------------------------------------
!PELLET_NGS_ABL calculates the rates of recession of the pellet surface from two
!  equations using the cloud thickness as the independent parameter
!
!References:
!  W.A.Houlberg, S.L.Milora, S.E.Attenberger, Nucl Fusion 28 (1988) 595
!  S.L.Milora, ORNL/TM-8616 (1983)
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
LOGICAL, INTENT(IN) :: &
  l_newcell		  !logical to recalculate fast ion fluxes in plasma 
				 !=TRUE new plasma cell
				 !=FALSE same plasma cell as previous call

REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den,			 & !electron density [/m**3]
  xrhor,		   & !cloud thickness / rhosrp0 [-]
  tfd			  !collisionless self-limiting parameter =[-]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  rdot1,		   & !dr/dt from energy balance at pellet surface [m/s]
  rdot2,		   & !dr/dt from balance in cloud [m/s]
  difrdt			 !ratio of diff to avg of the dr/dt solutions [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i

REAL(KIND=rspec) :: &
  areae,areaf,areat,cloudi,cloudn,fqf,q,qfo, &
  qtc,			 & !total heat flux at neutral cloud surface [keV/m**2/s]
  qtp,			 & !total heat flux at pellet surface [keV/m**2/s]
  xmp,			 & !molecular mass of pellet species [kg]
  xr

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Molecular mass
xmp=2*amup_pl*z_protonmass

!Total effective pellet surface areas
areat=4*z_pi*rp**2

!Effective area parallel to field line for electrons
areae=areat/2

!Effective area perpendicular to field lines for fast ions
areaf=areat/2

!-------------------------------------------------------------------------------
!Heat fluxes to cloud and pellet surfaces
!-------------------------------------------------------------------------------
!Neutral and plasma cloud thicknesses
cloudn=ellipt_pl*xrhor*rhosrp0_pl/xmp
cloudi=fionc_pl*xrhor*rhosrp0_pl/xmp

!Electrons
CALL PELLET_NGS_QE(te,den,cloudn,cloudi,tfd)
qtp=qeo_pl*fqe_pl*areae/areat
qtc=qeo_pl*(1.0-fqe_pl)
!Fast ions
IF(l_fast_pl) THEN

  DO i=1,nf_pl

	CALL PELLET_NGS_QF(l_newcell,cloudn,cloudi,izf_pl(i),amuf_pl(i), &
				 vcf_pl(i),nef_pl(i),ef_pl(:,i),denf_pl(:,i), &
				 efo_pl(:,i),qfo_pl(:,i),qfo, &
				 fqf)

	qtp=qtp+qfo*fqf*areaf/areat
	qtc=qtc+qfo*(1.0-fqf)

  ENDDO

ENDIF

IF(qtc < 0.0) qtc=0
q=qtc/(rhosrp0_pl*xrhor)

!-------------------------------------------------------------------------------
!Pellet ablation rate
!-------------------------------------------------------------------------------
rdot1=-qtp/(z_evap*denm_pl)
xr=rp/rpel_pl
rdot2=-1.25*xrhor/xr*(z_j7kv*q*rp*(z_gam-1.0)/2)**(1.0/3.0)
difrdt=2*(rdot2-rdot1)/(rdot2+rdot1)
END SUBROUTINE PELLET_NGS_ABL

SUBROUTINE PELLET_NGS_QE(te,den,cloudn,cloudi,tfd)
!-------------------------------------------------------------------------------
!PELLET_NGS_QE calculates the electron heat flux incident on the pellet cloud
!  and the heat flux attenuation factor in the cloud for the NGS model
!
!References:
!  W.A.Houlberg, S.L.Milora, S.E.Attenberger, Nucl Fusion 28 (1988) 595
!  S.L.Milora, ORNL/TM-8616 (1983)
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!
!Internal:
!  qeo_pl		  -incident electron heat flux from plasma [keV/m**2/s]
!  fqe_pl=qep/qeo_pl   -electron heat flux attenuation factor [-]
!  es			-energy below which elastic scattering is used [keV]
!  qep			 -electron heat flux at pellet surface [keV/m**2/s]
!  alfe			-cross-section for elastic scattering
!  ce			-mean electron thermal speed in plasma [m/s]
!  cjeo			-random electron particle flux in plasma [keV/m**2/s]
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  te,			& !electron temperature in plasma-keV]
  den,			 & !electron density in plasma-/m**3]
  cloudn,		  & !neutral cloud thickness-molecule/m**2]
  cloudi,		  & !ionized cloud thickness [ 0.5*electrons/m**2]
  tfd			  !collisionless self-limiting parameter [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  jg

REAL(KIND=rspec), SAVE :: &
  ae(2)=	 (/ 2.35e21,  4.0e21/), &
  we1(1)=	(/ 1.00000/), &
  we5(5)=	(/ 3.58220,  2.39789,  1.76644,  1.25078,  0.50209/), &
  we10(10)=  (/ 4.26829,  3.14987,  2.60654,  2.22568,  1.91562, &
			1.64694,  1.39401,  1.14336,  0.87304,  0.38355/), &
  we20(20)=  (/ 4.91742,  3.82498,  3.33808,  2.98877,  2.72599, &
			2.50164,  2.31666,  2.14625,  1.99227,  1.84801, &
			1.71277,  1.58720,  1.46082,  1.33618,  1.21265, &
			1.08658,  0.95348,  0.80956,  0.63799,  0.29548/)

REAL(KIND=rspec), PARAMETER :: &
  alfe=1.8e-20, &
  c=2.0, &
  es=0.1

REAL(KIND=rspec) :: &
  c1,c2,ce,cjeo,eo,eog,eogp,eogp2,enlam,epg,ethr2,explim,qeo0,qeog, &
  qeogp,qep,qepg,qesg,qex,qfactr,sgm,tej,tfactr,xi,xi2

REAL(KIND=rspec) :: &
  we(20)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Number of electron energy groups and weightingg factors
IF(neg_pl > 10) THEN

  neg_pl=20
  we(1:neg_pl)=we20(1:neg_pl)

ELSEIF(neg_pl > 5) THEN

  neg_pl=10
  we(1:neg_pl)=we10(1:neg_pl)

ELSEIF(neg_pl > 1) THEN

  neg_pl=5
  we(1:neg_pl)=we5(1:neg_pl)

ELSE

  neg_pl=1
  we(1)=we1(1)

ENDIF

!Other physics parameters
eo=1.5*te
tej=te*z_j7kv
ce=SQRT((8/z_pi)*tej/z_electronmass)
cjeo=den*ce/4
qeo0=4*cjeo*eo/3
qep=0
c1=ae(1)/ae(2)
c2=2*cloudn/ae(2)
enlam=10.0
sgm=enlam*4*z_pi*(z_coulomb/(4*z_pi*z_epsilon0))**2
ethr2=2*cloudi*sgm*(z_coulomb/z_j7kv)**2
!write(*, *) eo, te, ce, cjeo, qeo0, c1, c2, z_j7kv
qeo_pl=0

!-------------------------------------------------------------------------------
!Calculate the heat flux incident on the pellet surface
!-------------------------------------------------------------------------------
!Consider an arbitraty number of energy groups
DO jg=1,neg_pl

  !Set group factors
  tfactr=tfd*SQRT(we(jg))

  IF(tfactr > 0.01) THEN

	qfactr=(1.0-EXP(-tfactr))/tfactr

  ELSE

	qfactr=1

  ENDIF

  qeog=qeo0*qfactr/neg_pl
  qeo_pl=qeo_pl+qeog
  eog=eo*we(jg)
  eogp2=eog**2-ethr2
  IF(eogp2 < 0.0) eogp2=0
  eogp=SQRT(eogp2)
  qeogp=qeog*eogp/eog
  qesg=qeog*es/eog

  !Set xi2 and incident flux
  IF(eogp <= es) THEN

	!Scattering regime for eogp < es
	xi2=0
	qex=qeogp

  ELSE

	!Scattering regime for eogp > es
	xi2=ae(1)*(eogp-es)+ae(2)*(eogp**2-es**2)/2
	qex=qesg

  ENDIF

  !Check whether epg > or < es
  xi=alfe*(cloudn-xi2)
  IF(xi >= 0.0) THEN

	!epg < es
	explim=85!LOG(3.4e38_rspec)
	IF(xi > explim) xi=explim
	qepg=qex*(c+1)/(c+EXP(xi))
!c	write(*, *) "1,", qex, c, xi
  ELSE

	!epg > es
	epg=-c1+SQRT((c1+eogp)**2-c2)
	qepg=qeog*epg/eog
!c	write(*, *) "2,", qeog, epg, eog 
  ENDIF
!write(*, *) "qepg: ", qepg
  !Add to energy flux at pellet surface
  qep=qep+qepg

ENDDO   

!Fraction of electron energy reaching pellet surface
fqe_pl=qep/qeo_pl
!write(*, *) "QE: ", qep, qeo_pl
END SUBROUTINE PELLET_NGS_QE

SUBROUTINE PELLET_NGS_QF(l_newcell,cloudn,cloudi,izf,af,vcf,nfg,efg,denfg, &
				 efgo,qfgo,qfo, &
				 fqf)
!-------------------------------------------------------------------------------
!PELLET_NGS_QF calculates the fast ion heat flux incident on the pellet cloud
!  and the heat flux attenuation factor in the cloud for either hydrogenic or
!  helium fast ions
!
!References:
!  W.A.Houlberg, S.L.Milora, S.E.Attenberger, Nucl Fusion 28 (1988) 595
!  S.L.Milora, ORNL/TM-8616 (1983)
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  Input/output variables are input when l_newcell=.FALSE
!  ah(3)		   -parameters for fit to fast H+ energy loss in H2 gas
!  ahe(2)		  -parameters for fit to fast He+ energy loss in H2 gas
!  qfp			 -fast ion heat flux to pellet surface [kev/m**2/s]
!-------------------------------------------------------------------------------

!Declaration of input variables
LOGICAL, INTENT(IN) :: &
  l_newcell		  !flag for entry to new cell [logical]
				 !(recalculate fast ion fluxes)

INTEGER :: &
  izf,			 & !fast ion charge number [-]
				 !=1 fast hydrogen ions
				 !=2 fast helium ions
  nfg			  !number of fast energy groups [-]

REAL(KIND=rspec), INTENT(IN) :: &
  af,			& !atomic mass number of fast ions [-]
  vcf,			 & !critical velocity in classical thermalization [m/s]
  cloudn,		  & !neutral cloud thickness [molecule/m**2]
  cloudi			 !ionized cloud thickness [0.5*electrons/m**2]

REAL(KIND=rspec), INTENT(IN) :: &
  efg(:),		  & !fast ion energy group boundaries, efg(j)>efg(j+1), [keV]
  denfg(:)		   !fast ion density in energy intervals [/m**3]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  qfo			  !fast ion heat flux to neutral cloud [keV/m**2/s]

REAL(KIND=rspec), INTENT(INOUT) :: &
  efgo(:),		 & !average energy in energy group [keV]
  qfgo(:)			!energy flux in energy group [keV/m**2/s]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  fqf			  !heat flux attenuation factor at pellet (=qfp/qfo) [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  jg

REAL(KIND=rspec), SAVE :: &
  ah(3)=	 (/ 1.0e-20,  30.93,	3.239/), &
  ahe(2)=	(/ 2.02e-22, 0.325/)

REAL(KIND=rspec) :: &
  arhmg,b,bb,c,c1,c2,cong,efcon,efgop,efgop2,efgot,egp,enlam,ethr2, &
  pemr,qfcon,qfo0,qfp,sgm,sq3,th1,th2,th3,tl1,tl2,tl3,vl,xlower, &
  xupper,xvh2,xvh3,xvl,xvl2,xvl3

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Check whether fast ion charge is appropriate 
IF((izf < 1) .OR. &
   (izf > 2)) GOTO 9999

!-------------------------------------------------------------------------------
!Set fast ion parameters if a new cell has been entered
!-------------------------------------------------------------------------------
IF(l_newcell) THEN

!Constants
  qfo0=0
  sq3=SQRT(3.0)

!Velocities and integral expressions at low end of group
  vl=SQRT(2*efg(1)*z_j7kv/(af*z_protonmass))
  xvl=vl/vcf
  xvl2=xvl**2
  xvl3=xvl*xvl2
  tl1=LOG(1.0+xvl3)
  tl2=LOG((1.0-xvl+xvl2)/(1.0+xvl)**2)
  tl3=ATAN((2*xvl-1.0)/sq3)

!Sum contribution over energy groups
  DO jg=1,nfg !Over fast ion energy groups

	!Shift expressions to boundaries of next lower energy group
	vl=SQRT(2*efg(jg+1)*z_j7kv/(af*z_protonmass))
	xvl=vl/vcf
	xvh2=xvl2
	xvl2=xvl**2
	xvh3=xvl3
	xvl3=xvl*xvl2
	th1=tl1
	tl1=LOG(1.0+xvl3)
	th2=tl2
	tl2=LOG((1.0-xvl+xvl2)/(1.0+xvl)**2)
	th3=tl3
	tl3=ATAN((2*xvl-1.0)/sq3)

	!Calculate mean group energy in plasma
	cong=3/LOG((1.0+xvh3)/(1.0+xvl3))
	efcon=af*z_protonmass*vcf**2/2/z_j7kv
	xupper=xvh2/2-(th2/2+sq3*th3)/3
	xlower=xvl2/2-(tl2/2+sq3*tl3)/3
	efgo(jg)=efcon*cong*(xupper-xlower)

	!Calculate mean group heat flux in plasma
	qfcon=denfg(jg)*efcon*(vcf/12)
	qfgo(jg)=qfcon*cong*((xvh3-th1)-(xvl3-tl1))
	qfo0=qfo0+qfgo(jg)

  ENDDO !Over fast ion energy groups

ENDIF

!-------------------------------------------------------------------------------
!Heat flux to pellet surface
!-------------------------------------------------------------------------------
!Initialize fluxes
qfo=0
qfp=0
fqf=0

!Check whether there are any fast ions present
IF(SUM(qfgo(1:nfg)) <= 0.0) GOTO 9999

!Set constants for all energy groups
pemr=af*z_protonmass/z_electronmass
enlam=10.0
sgm=enlam*4*z_pi*pemr*(z_coulomb/(4*z_pi*z_epsilon0))**2
ethr2=2*cloudi*sgm*(z_coulomb/z_j7kv)**2

!Set threshhold energy to zero for ions normal to ionized cloud
ethr2=0

!Determine fast ion species and get heat fluxes at pellet surface
IF(izf == 1) THEN

  !Hydrogenic ions
  b=ah(2)*SQRT(af)/ah(3)
  bb=b**2
  arhmg=cloudn*ah(1)/ah(3)

  !Calculate parameters incident on the pellet
  DO jg=1,nfg !Over fast ion energy groups

	egp=0
	efgop2=efgo(jg)**2-ethr2
	IF(efgop2 < 0.0) efgop2=0
	efgop=SQRT(efgop2)
	qfo=qfo+qfgo(jg)*efgop/efgo(jg)
	c=2*SQRT(efgop)/b+(efgop-arhmg)/bb
	IF(c > 0.0) egp=bb*(-1.0+SQRT(1.0+c))**2
	qfp=qfp+qfgo(jg)*egp/efgo(jg)

  ENDDO !Over fast ion energy groups

  fqf=qfp/qfo

ELSEIF(izf == 2) THEN

  !Helium ions
  c1=1.0-ahe(2)
  c2=cloudn*(4/af)**ahe(2)*ahe(1)*c1
  ethr2=4*ethr2

  !Calculate parameters incident on the pellet
  DO jg=1,nfg !Over fast ion energy groups

	egp=0
	efgop2=efgo(jg)**2-ethr2
	IF(efgop2 < 0.0) efgop2=0
	efgop=SQRT(efgop2)
	qfo=qfo+qfgo(jg)*efgop/efgo(jg)
	efgot=efgop**c1
	egp=MAX(efgot-c2,0.0_rspec)**(1/c1)
	qfp=qfp+qfgo(jg)*egp/efgo(jg)

  ENDDO !Over fast ion energy groups

	fqf=qfp/qfo

ELSE

  !Default calculation for izf > 2
  qfp=0
  fqf=0

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_NGS_QF




SUBROUTINE PELLET_PARKS(rp,te,den, &
				rdot)
!-------------------------------------------------------------------------------
!PELLET_PARKS determines the pellet ablation rate using a fit generated by Parks
!  based on IPADBASE results
!
!References:
!  P.B.Parks, M.N.Rosenbluth, Plasma Phys. (1998) 1380
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den			  !electron density [/m**3]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  rdot			 !rate of change in pellet radius [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec), PARAMETER :: &
  temin=1.0e-3

REAL(KIND=rspec) :: &
  dinf,rpcm,tinf,dndt,loglamH,Ihyd

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Check pellet size and minimum electron temperature
IF((rp < (z_tolr*rpel_pl)) .OR. &
   (te < temin)) THEN

  rdot=0
  GOTO 9999

ENDIF

!Convert to Park's units of cm and eV
tinf=te*1.0e3
rpcm=rp*1.0e2
dinf=den*1.0e-6

Ihyd = 7.514   ! Effective hydrogenic excitation energy
loglamH = log(2*tinf/Ihyd)

!-------------------------------------------------------------------------------
!Calculate ablation rate
!-------------------------------------------------------------------------------
!Compute the ablation rate in cm/s, Eqn 2
rdot=-8.2e15*(amup_pl**(-0.333)*dinf**(0.333)*rpcm**(-0.666)*tinf**(1.833))/ &
		 ((loglamH**0.666) *4*z_pi*(2*(denm_pl * 1e-6)))

!Compute dr/dt in m/s
rdot=rdot*0.01

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_PARKS



SUBROUTINE PELLET_PARKSQ(rp,te,den, &
				rdot)
!-------------------------------------------------------------------------------
!PELLET_PARKSQ determines the pellet ablation rate using a Q valid for all Te
!
!References:
!  P.B.Parks, MR.J. Turnbull, Plasma Phys. (1978) 1735
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den			  !electron density [/m**3]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  rdot			 !rate of change in pellet radius [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec), PARAMETER :: &
  temin=1.0e-3

REAL(KIND=rspec) :: &
  dinf,rpcm,tinf,dndt, &
  gam,lam,ehat,rhat,qhat,Q,Est,sig,Lst,Loss,A,lamst,Coef

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Check pellet size and minimum electron temperature
IF((rp < (z_tolr*rpel_pl)) .OR. &
   (te < temin)) THEN

  rdot=0
  GOTO 9999

ENDIF

!Convert to Park's units of cm and eV
tinf=te*1.0e3
rpcm=rp*1.0e2
dinf=den*1.0e-6

! Determination of coefficients
gam = 7/5.
lam = 0.961
ehat = 1.152
rhat = 0.645
qhat = 1.56
Q = 0.65

Est = 2*tinf/ehat
sig = (8.8e-13/ Est**1.71) - 1.62e-12/Est**1.932
A = 1./((Est/100)**0.823 + 1/((Est/60)**0.125) + 1/((Est/48)**1.94))
Lst = 8.62e-15*A
Loss = 2*Lst/Est
lamst = sig + loss
Coef = ((gam-1)**0.333 * lam)/(rhat**1.333 * qhat**0.333)

!-------------------------------------------------------------------------------
!Calculate ablation rate
!-------------------------------------------------------------------------------
!Compute the ablation rate in atoms/s, Eqn 2
dndt=5.0586e7*Coef*Q**0.333*(amup_pl**(-0.333)*dinf**(0.333)*rpcm**(1.333)* &
	 tinf**(0.5))/lamst**0.666

!Compute dr/dt in m/s
rdot=-dndt/((2*denm_pl)*4*z_pi*rp**2)


!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_PARKSQ




SUBROUTINE PELLET_PARKS_IM(rp,te,den, &
				   rdot)
!-------------------------------------------------------------------------------
!PELLET_PARKS_IM determines the pellet ablation rate for impurity pellets using
!  the model by Parks
!
!References:
!  P.B.Parks, Memo (Jul-1995)
!  B.V.Kuteev, Y.Yu.Sergeev, A.Yu.Kostrukov, V.A.Segal, O.A.Bakhareva,
!	P.B.Parks, Bull Am Phys Soc 44 (1999) 115
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  dena			-atomic density [atoms/m**3]
!  rhomass		 -solid mass density [g/m**3]
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  rp,			& !pellet radius [m]
  te,			& !electron temperature [keV]
  den			  !electron density [/m**3]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  rdot			 !rate of change in pellet radius [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec), PARAMETER :: &
  temin=1.0e-3

REAL(KIND=rspec) :: &
  coeff,dena,dinf,dndt,rhomass,rpcm,tinf

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Check pellet size and minimum electron temperature
IF((rp < (z_tolr*rpel_pl)) .OR. &
   (te < temin)) THEN

  rdot=0
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Calculate ablation rate
!-------------------------------------------------------------------------------
!Convert to Park's units of cm and eV
tinf=te*1.0e3
rpcm=rp*1.0e2
dinf=den*1.0e-6

!Determine coefficient to use (function of impurity)
IF(amup_pl >= 20.0 .AND. &
   amup_pl <= 21.0) THEN

  !Neon (nominally 20.2)
  coeff=2.343e15
  rhomass=2.205e6

ELSEIF(amup_pl >= 39.0 .AND. &
	 amup_pl <= 40.0) THEN

  !Argon (nominally 39.9)
  coeff=1.583e15
  rhomass=1.400e6

ELSEIF(amup_pl >= 83.0 .AND. &
	 amup_pl <= 84.0) THEN

  !Krypton (nominally 83.8)
  coeff=0.997e15
  rhomass=2.155e6

ELSEIF(amup_pl >= 130.5 .AND. &
	 amup_pl <= 131.5) THEN

  !Xenon (nominally 131.0)
  coeff=0.801e15
  rhomass=3.520e6

ELSE

  coeff=0
  rhomass=1.0e6

ENDIF

dena=rhomass*z_navogadro/amup_pl
	
!Compute the ablation rate, in atoms/s
dndt=coeff*tinf**(1.64)*dinf**(0.333)*rpcm**(1.333)*amup_pl**(-0.333)

!Compute dr/dt in m/s
rdot=-dndt/(dena*4*z_pi*rp**2)

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_PARKS_IM

SUBROUTINE PELLET_RK4(dvol,dt, &
				t,rp,te,den, &
				dden,iflag,message)
!-------------------------------------------------------------------------------
!PELLET_RK4 is a fourth order Runge-Kutta integration routine that updates the
!  pellet radius and plasma parameters as it passes through a cell
!
!References:
!  W.A.Houlberg, S.L.Milora, S.E.Attenberger, Nucl Fusion 28 (1988) 595
!  W.A.Houlberg, L.R.Baylor 6/2004
!  W.A.Houlberg, F90 free format 8/2004
!
!Comment:
!  l_newcell		-flag for entry to new cell [logical]
!				 -recalculate fast ion fluxes
!  nstep			-number of time steps to get through this cell [-]
!  fte			-estimated fractional change in te [-]
!  fdts			 -time step multiplier to keep d(te)/te < z_tolt [-]
!  s			  -function to evaluate density perturbation [/m**3]
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  dvol,			& !volume of this cell [m**3]
  dt			   !time pellet spends in this cell [s]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  t,			 & !time since injection at entry/exit [s]
  rp,			& !pellet radius at entry/exit [m]
  te,			& !electron temperature at entry/exit [keV]
  den			  !electron density at entry/exit [/m**3]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message			!warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag			!error and warning flag [-]
				 !=-1 warning
				 !=0 no warnings or errors
				 !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  dden			 !incremental electron density [/m**3]

!-------------------------------------------------------------------------------
!Declaration of local variables
LOGICAL :: &
  l_newcell

INTEGER, PARAMETER :: &
  mstep=500

INTEGER :: &
  i,istep,nstep

REAL(KIND=rspec), PARAMETER :: &
  adcon=25.0

REAL(KIND=rspec) :: &
  adfac,areac,areap,ce,ddens,dens,dnte,drav,drp3,dtold,dts,endot, &
  fdts,fioncn,fnmax,fntmax,fte,rdot,rdotav,rdotmx,rp0,rp3,rpmin, &
  rps,tcl,temin,tes,tfd,tnew,told

REAL(KIND=rspec) :: &
  dr(4)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Set parameters assuming a single step through the cell
l_newcell=.TRUE.
told=t
tnew=told+dt
dts=dt
rp0=rp
dden=0
temin=z_eion
nstep=0

!Set maximum delta(n)/n, pellet size at that max, and the max ablation rate
fntmax=(te-temin)/(z_eion*2/3+temin)
IF(fntmax <= 0.0) fntmax=0
fnmax=MIN(20.0_rspec,fntmax)
drp3=(3/(8*z_pi))*fnmax*dvol*(den/denm_pl)
rp3=rp**3

IF(drp3 < rp3) THEN

  rpmin=(rp3-drp3)**(1.0/3.0)

ELSE

  rpmin=0

ENDIF

rdotmx=-(rp0-rpmin)/dt

!-------------------------------------------------------------------------------
!Advance pellet through time interval in cell using up to mstep time steps
!-------------------------------------------------------------------------------
step_loop: DO istep=1,mstep  !Over time steps

  !Initialize intermediate runge-kutta parameters and step flag
  rps=rp
  nstep=nstep+1

  !Employ 4th order Runge-Kutta scheme for each step
  DO i=1,4 !Over RK steps

	!Adiabatic self-limiting ablation for large perturbations
	dden=2*denm_pl*4*z_pi/3/dvol*(rp0**3-rps**3)
	adfac=EXP(-(dden/(adcon*den))**2)
	dden=(1.0-adfac)*dden
	dens=den+dden
	tes=(den*te-2*dden*z_eion/3)/dens
	IF(tes < z_eion) tes=z_eion

	!Check whether to use collisionless self-limiting ablation
	IF(neg_pl == 1) THEN

	!Not applicable for single electron energy group
	tfd=0

	ELSE

	!Collisionless self-limiting ablation for multiple groups
	dtold=tnew-told
	ce=SQRT(8*tes*z_j7kv/z_pi/z_electronmass)
	tfd=dtold*(ce/4)*2*z_pi*rcl_pl**2/dvol

	ENDIF

	!Choose appropriate ablation model
	!Hydrogenic pellets
	IF(k_pel_pl >= 0 .AND. k_pel_pl <= 2) THEN

	!Neutral Gas Shielding model and variations
	print*, "te into ngs = ", tes
	CALL PELLET_NGS(l_newcell,rps,tes,dens,tfd, &
				rdot,iflag,message)

	IF(iflag > 0) THEN

	  message='PELLET_RK4/'//message
	  GOTO 9999

	ENDIF

	ELSEIF(k_pel_pl == 3) THEN

	!Macaulay hydrogenic model
	CALL PELLET_MAC(rps,tes,dens, &
				rdot)

	ELSEIF(k_pel_pl == 4) THEN

	!Kuteev hydrogenic model
	CALL PELLET_KUT(rps,tes,dens, &
				rdot)

	ELSEIF(k_pel_pl == 5) THEN

	!Parks hydrogenic model
	CALL PELLET_PARKS(rps,tes,dens, &
				rdot)

	ELSEIF(k_pel_pl == 6) THEN

	!Parks hydrogenic model for arbitrary heating coef Q
	CALL PELLET_PARKSQ(rps,tes,dens, &
				rdot)

	!Impurity pellets
	ELSEIF(k_pel_pl == 10) THEN

	!Parks impurity model
	CALL PELLET_PARKS_IM(rps,tes,dens, &
				   rdot)

	ELSEIF(k_pel_pl == 11) THEN

	!Kuteev impurity model
	CALL PELLET_KUT_IM(rps,tes,dens, &
				 rdot)

	ENDIF

	IF(ABS(rdot) > ABS(rdotmx)) rdot=rdotmx

	!Update neutral cloud size	 
	areap=4*z_pi*rps**2
	endot=-2*denm_pl*areap*rdot

	IF(k_pel_pl == 2) THEN

	rcl_pl=rps+z_rclmin
	tcl=2*rcl_pl/vpel_pl
	areac=2*z_pi*rcl_pl**2
	fioncn=endot*tcl/(areac*2*denm_pl*rpel_pl)/xrhoro_pl
	fionc_pl=(fionc_pl+fioncn)/2

	ENDIF

	dr(i)=dts*rdot

	IF(i <= 2) THEN

	IF((rp+dr(i)/2) < rpmin) dr(i)=2*(rpmin-rp)
	rps=rp+dr(i)/2

	ELSE

	IF((rp+dr(i)) < rpmin) dr(i)=rpmin-rp
	rps=rp+dr(i)

	ENDIF

	IF(i == 1) THEN

	!Check perturbation on background plasma
	IF(rps < z_tolr*rpel_pl/2) rps=z_tolr*rpel_pl/2
	ddens=2*denm_pl/dvol*4*z_pi/3 *(rp**3-rps**3)
	dnte=2*ddens/3*z_eion
	fte=2*(ddens+dnte/tes)/(dens+2*ddens)
	fdts=1
	IF(fte /= 0.0) fdts=z_tolt/fte

	IF(fdts < 1.0) THEN

	  !Reduce time step size and try again
	  dts=0.9*fdts*dts
	  CYCLE step_loop

	ENDIF

	ENDIF

  ENDDO !Over RK steps

  !Completed step
  drav=(dr(1)+2*(dr(2)+dr(3))+dr(4))/6
  rdotav=drav/dts

  IF((rp+drav) < (z_tolr*rpel_pl))THEN

	dts=-rp/rdotav
	drav=-rp

  ENDIF

  rp=rp+drav
  dden=2*denm_pl/dvol*4*z_pi/3*(rp0**3-rp**3)
  t=t+dts
  dts=tnew-t

  IF((dts <= 1.0e-6*tnew   ) .OR. &
	 (rp  <= 1.0e-6*rpel_pl)) GOTO 9999

ENDDO step_loop

!Exceeded maximum number of steps in cell
iflag=1
message='PELLET_RK4/ERROR(2):exceeded max iterations'

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_RK4

subroutine pelfil (ndat, nin, versid)
integer, intent(in) :: ndat, nin
integer :: i
!
! --- create pellet plot data file, and copy inone file text to it
!
      character*8  title(10)
      character *(*) versid
!
      call DESTROY ('peldat')
      open   (unit = ndat, file = 'peldat', status = 'UNKNOWN')
      write  (ndat, 8010) versid
  8010 format (' ******** ', a, ' ******** (33x65)' ///)
      rewind (unit = nin)
!
 2000 read   (nin , '(    9a8)', end=2020)  (title(i), i=1,9)
      write  (ndat, '(1x, 9a8)'          )  (title(i), i=1,9)
      go to 2000
!
2020	write  (ndat, '(''stop'')')
      rewind (unit = nin)
      return
!
end subroutine pelfil

SUBROUTINE propel

! ----------------------------------------------------------------------
!	 propel	   mar. 11, 1984	 Dave S., Bob S.
!
!	 ONETWO interface for Oak Ridge pellet fueling model.  the pellet
!	 model needs electron density and temperature profiles, and fast
!	 ion density profiles resolved into energy bins.  the model returns
!	 the increase in electron density in each plasma shell due to the
!	 complete ablation of a pellet.
!
!	 other ONETWO inputs via common /pelcom/:
!	 ipel   - pellet ion species index
!	 pelrad - initial pellet radius (cm)
!	 vpel   - pellet velocity (cm/sec)
!	 nbgpel - number of fast beam ion energy bins
!	 ipelet - 0 for no pellet; 1 to inject pellet
!	 npel   - pellet counter
!	 nupel  - I/O channel for pellet plot data (common /io/)
! ----------------------------------------------------------------------
!
!

	USE param
	USE fusion
	USE io 
	USE ions
	USE neut
	USE nub
	USE nub2
	USE solcon
	USE soln
	USE numbrs
	USE mesh
	USE sourc
	USE machin
	USE geom
	USE events
	USE pelcom

!	implicit integer (i-n), real*8 (a-h, o-z)
	integer, dimension(:), allocatable :: ne_f, iz_f
	real(KIND=rspec), dimension(:), allocatable :: amu_f
!	real(KIND=rspec) :: ne_f(nprim), iz_f(nprim), amu_f(nprim)
	real(KIND=rspec), dimension(:, :), allocatable :: e_ef, vc_rf
	real(KIND=rspec), dimension(:, :, :), allocatable :: den_ref 
	character(len = 63) :: message
!
!	include 'param.i'
!	include 'events.i'
!	include 'fusion.i'
!	include 'geom.i'
!	include 'io.i'
!	include 'ions.i'
!	include 'machin.i'
!	include 'mesh.i'
!	include 'neut.i'
!	include 'nub.i'
!	include 'nub2.i'
!	include 'numbrs.i'
!	include 'pelcom.i'
!	include 'solcon.i'
!	include 'soln.i'
!	include 'sourc.i'
!
!	 kj2 = 2*kj
!	 kbg = max. number of fast beam ion energy groups
!	 kag = max. number of fast alpha energy groups
!
	integer, parameter :: kj2 = kj * 2, kbg = 25, kag = 2
	real(KIND=rspec), dimension(1 : kj) :: dvol, elden, teev, pden, workp
	integer, dimension(1 : kj2) :: mapil
	real(KIND=rspec), dimension(1 : kj2) :: rchord
	real :: ebg(1 : kbg), vcb(1 : kj), denbig(1 : kj, 1 : kbg), denbm(1 : kbg)
	real :: eag(1 : kag), vca(1 : kj), denaig(1 : kj, 1 : kag)
	integer :: ibg, kprint, nprt, ndat, njm1, nseg, i, j, jb, nbg, jbm, jbe, nbge, jb1, jg, jbg, kpel, nag, iflag, k
	real, parameter :: pi = 3.14159265
	real :: debg, enbx, sum, rmass, rnfdot, engy, v0, vgy, tscat, svncx, cxr, taucx, rksum, atwalp, azf, eion, totnum, denion, totnew, totadd
!
!	 check array dimensions
!
	if (ipelet .eq. 0) return
	if (2 * nj .gt. kj2) call STOP('subroutine PROPEL: problem #1', 59)
	if (nbgpel .gt. kbg) call STOP ('subroutine PROPEL: problem #2', 60)
	ibg = nbgpel / 6
	if (ibg * 6 .ne. nbgpel) call STOP ('subroutine PROPEL: problem #3', 61)
!
!	 set output unit numbers
!
	kprint = 3
	nprt = nqik
	ndat = nupel
!
!	 increment pellet counter and set next pellet time
!
	npel = npel + 1
	write (ndat, 90) npel, time
	write (ncrt, 90) npel, time
90	format (' pellet number ', i3, ' at time = ', f7.3, ' seconds')
	timevent(8) = timpel(npel + 1)
	if (npel .eq. 10)  timevent(8) = 1000.0
!
!		ONETWO has nj grid points enclosing njm1 plasma shells.
!		the pellet is assumed to pass horizontally through
!		the plasma center, thus it sees nseg = 2*njm1 plasma segments.
!		initialize plasma shell volumes, and map pellet path
!		segments onto the plasma shells.
!
	njm1 = nj - 1
	nseg = 2 * njm1
	do j = 1, njm1
		dvol(j)  = 2.0 * pi ** 2 * rmajor * (r(j+1) ** 2 - r(j) ** 2)
	end do
	dvol(nj) = 0.0
	do i = 1, njm1
		mapil(nj - i) = i
		mapil(njm1 + i) = i
	end do
	do i = 1, nj
		rchord(i) = (r(nj) - r(nj + 1 - i)) / SQRT (kappa)
	end do
	do i = nj + 1, nseg + 1
		rchord(i) = (r(nj) + r(i - nj + 1)) / SQRT (kappa)
	end do
!
!		the pellet model uses the following strange grid:
!		te(1)  =  te(r=0)
!		te(i) = te(r = (r(i)+r(i+1))/2)
!		te(nj)  =   te(r=a)
!		convert ONETWO profiles to this grid.
!
	elden(1) = ene(1)
	teev(1) = te(1) * 1000.0
	do i = 2, njm1
		elden(i) = 0.5 * (ene(i) + ene(i + 1))
		teev(i) = 0.5 * (te(i) + te(i + 1)) * 1000.0
	enddo
	elden(nj) = ene(nj)
	teev(nj)  = te(nj) * 1000.0
!
!		initialize beam fast ion energy groups.
!		there are nbgpel energy intervals of width debg from
!		the full injection energy down to zero.  nbgpel should
!		be a multiple of 6 to allow a clean match with the
!		half and third beam energy components.  note that
!		the interval from debg to 0 is always left empty.
!
	ebg(1) = ebkev(1) * 1000.0
	debg = ebg(1) / nbgpel
	nbg	= nbgpel - 1
	do jg = 2, nbg
		ebg(jg)	= ebg(1) - (jg - 1) * debg
	end do
	ebg(nbg + 1) = debg
!
!	 at each point in space, sum the beam particle density
!	 from each energy component of each beamline.
!	 model the fast ion slowing down by an analytical expression
!	 derived by Callen and Rome, assuming steady state conditions.
!
	call zeroa(denbig, kj * kbg)
	call zeroa(denbm, kbg)
	do j = 1, nj
		vcb(j) = 1.385e6 * SQRT (14.8 * teev(j))
		do jbm = 1, nbeams
			do 210 jbe = 1, 3
				nbge = nbgpel / jbe - 1
				jb1  = nbg - nbge + 1
				enbx = 0.5 * (enb(j, jbe, jbm) + enb(j + 1, jbe, jbm))
				if (j .eq. 1 .or. j .eq. nj)  enbx = enb(j, jbe, jbm)
				sum  = 0.0
				if (enbx .eq. 0.0)  go to 210
				do jbg = jb1, nbg
					rmass = atw_beam * 1.6726e-24
					rnfdot = sbsav(j, jbe, jbm)
					engy = ebg(jbg)
!					e0 = ebg(1) / jbe
!					encrit = 0.5 * atw_beam * 1.6726e-24 * 1.6022e-12 * vcb(j) ** 2
					v0 = SQRT(2.0 * ebg(1) * 1.6022e-12 / (jbe * rmass))
					vgy	= SQRT (2.0 * engy * 1.6022e-12 / rmass)
					tscat = taus(j)
					svncx = (enn(j, 1) + enn(j, 2)) * cxr(engy * 0.001 / atw_beam)
					taucx = 1.0e30
					if (svncx .gt. 0.0)  taucx = 1.0 / svncx
					denbm(jbg) = (rnfdot * tscat / rmass) * vgy * (1.0 / (vgy ** 3 + vcb(j) ** 3)) * (((v0 ** 3 + vcb(j) ** 3) / (vgy ** 3+vcb(j) ** 3)) ** (tscat / taucx)) * debg * 1.6022e-12
					sum = sum + denbm(jbg)
				end do
				rksum = enbx / sum
				do jbg = jb1, nbg
					denbig(j, jbg) = denbig(j, jbg) + denbm(jbg) * rksum
				end do
210			continue
		end do
	end do
!
!	 neglect fast alpha particles for now
!
	nag	= 1
	atwalp = 4.0
	call zeroa (eag,kag)
	call zeroa (vca,kag)
	call zeroa (denaig,kag*kj)
!
!	 fling the pellet
!
	kpel = pelmod;
	iflag = -2;
	!message = ''
!PELLET(k_pel,amupel,rpel,vpel,n_r,dvol_r,den0_r,te0_r,n_p,map_p,s_p,pden_r,iflag,message,
!NF,NE_F,IZ_F,AMU_F,E_EF,VC_RF,DEN_REF,PEL_IONS,T_P,RPEL1_P,SRC_P,DEN0_P,DEN1_P,TE0_P,TE1_P,
!R0,A0,BT0,NCSOL,K_PRL,NPRLCLD,IPRLCLD,FPELPRL,PRLINJANG,PRLQ_R,PRLDEP)
	allocate(ne_f(nprim), iz_f(nprim), amu_f(nprim))
	if (nag .gt. nbg) then
		allocate(e_ef(nag, nprim), vc_rf(nag, nprim), den_ref(nag, nprim, 2))
	else
		allocate(e_ef(nbg, nprim), vc_rf(nbg, nprim), den_ref(nbg, nprim, 2))
	endif
	!if two specieses, then print message that He is excluded
		do i = 1,  nprim
			if (namep(i) .eq. 'he') then
				ne_f(i) = nag
				iz_f(i) = azf!(?)
				amu_f(i) = atwalp
				do j = 1, nag
					e_ef(j, i) = eag(j)
					vc_rf(j, i) = vca(j)
					do k = 1, 2
						den_ref(j, k, i) = denaig(j, k)
					end do
				end do			
			else 
				ne_f(i) = nbg
				iz_f(i) = azf
				amu_f(i) = atw_beam
				do j = 1, nbg
					e_ef(j, i) = ebg(j)
					vc_rf(j, i) = vcb(j)
					do k = 1, 2
						den_ref(j, k, i) = denbig(j, k)
					end do
				end do
			endif
		end do
		call pellet(kprint = kprint, &
			nprt = nprt, &
			ndat = ndat, &
			k_pel = kpel, &
			amupel = atw(ipel), &
			rpel = pelrad * 0.01, &
			vpel = vpel * 0.01, &
			n_r = nj, &
			dvol_r = dvol * 1e-6, &
			den0_r = elden * 1e6, &
			te0_r = teev * 1e-3, &
			n_p = nseg, &
			map_p = mapil, &
			s_p = rchord * 0.01, &
			pden_r = pden, &
			iflag = iflag, &
			message = message, &
			nf = nprim, &
			ne_f = ne_f, &
			iz_f = iz_f, &
			amu_f = amu_f, &
			e_ef = e_ef * 1e-3, &
			vc_rf = vc_rf * 1e-2, &
			den_ref = den_ref * 1e6, &
			k_prl = 0)
	!endif
!
!		convert pden to ONETWO grid
!
	pden = pden * 1e-6
	call copya (pden, workp, nj)
	pden(1)  = workp(1)
	pden(2)  = (workp(1) + 2.0*workp(2))/3.0
	do 250 i=3,njm1
 250  pden(i)  = 0.5 * (workp(i-1) + workp(i))
	pden(nj) = workp(nj)
!
!		calculate new plasma ene and te profiles.
!		note: 32.6 ev/pellet particle is dissipated in melting
!		and ionizing the pellet.  if fast beam ions are present,
!		it is assumed that they do all of the ablation, but the
!		ablation energy is not subtracted from the fast ion energy.
!
	eion = 0.0326
!
	do j=1,nj
	  if (enb(1,1,1) .ne. 0.0) then
	   te(j)  = te(j)*ene(j)/(ene(j)+pden(j))
	  end if
          if (enb(1,1,1) .eq. 0.0) then
	   te(j)  = (te(j)*ene(j)-eion*pden(j))/(ene(j)+pden(j))
	  end if
          ene(j) = ene(j) + pden(j)
	end do
!
!	 calculate new plasma ion density and temperature profiles
!
	do 310 j=1,nj
 310  workp(j) =  4.0 * pi**2 * rmajor*hcap(j)
	call trapv (r,en(1,ipel),workp,nj,totnum)
	do j=1,nj
	  denion = 0.0
	  do 320 i=1,nion
 320	denion	 = denion + en(j,i)
	  ti(j)	=  ti(j)*denion/(denion+pden(j))
	  en(j,ipel) =  en(j,ipel) + pden(j)
	end do
!
	call trapv (r,en(1,ipel),workp,nj,totnew)
!
!	 update the total number of ions for species ipel
!
	totadd = totnew - totnum
	if (ineut(ipel) .ne. 0)  snaddt(ipel) =  snaddt(ipel) + totadd
	return
!
	end subroutine propel

	real*8 function ssumpl (n, sx, incx)
!
!	implicit  integer (i-n), real*8 (a-h, o-z)
!
! ----------------------------------------------------------------------
! --- ssumpl adds up the elements of the array sx.
! --- last revision: 3/81 w.a.houlberg and s.e.attenberger ornl.
! --- other comments:
! --- use OMNILIB version for CRAY optimization
! ----------------------------------------------------------------------
!
	integer, intent(in) :: n
	integer, intent(in) :: incx
	real, intent(in) :: sx(:)

	real :: sum = 0.0
	integer :: i
	do i = 1, n, incx
	  sum = sum + sx(i)
	end do
	ssumpl = sum
	return
!
	end function


END MODULE PELLET_MOD
