!------------------------------------------------------------
! locpargen.f90
!
! PURPOSE:
!  Generate list of local parameters given input radius or 
!  input rho.
!------------------------------------------------------------

program locpargen

  use EXPRO_interface
  use GEO_interface

  implicit none

  integer :: j1,j2
  real :: r0
  real :: rho0
  real :: psi0
  real :: a
  real, dimension(1) :: x
  real, dimension(1) :: y
  real, dimension(5) :: z
  real, dimension(:), allocatable :: x_vec
  real, dimension(:,:,:), allocatable :: geo_p
  real :: ar, sf, shear

  open(unit=1,file='input.locpargen',status='old')
  read(1,*) r0
  read(1,*) rho0
  read(1,*) psi0 
  read(1,*) z(1)
  read(1,*) z(2)
  read(1,*) z(3)
  read(1,*) z(4)
  read(1,*) z(5)
  close(1)

  EXPRO_ctrl_n_ion = 5
  EXPRO_ctrl_quasineutral_flag = 0
  EXPRO_ctrl_z(1:5) = z(1:5)
  ! We don't need the numerical eq. flag set for this routine.
  EXPRO_ctrl_numeq_flag = 1

  call EXPRO_alloc('./',1) 
  call EXPRO_read


  print '(a)','INFO: (locpargen) Local input parameters:'
  print *

  allocate(x_vec(EXPRO_n_exp))

  ! Minor radius
  a = EXPRO_rmin(EXPRO_n_exp)

  if (r0 > 0.0) then

     ! Use local radius (r/a)

     x(1)  = r0
     x_vec = EXPRO_rmin/a

     ! RADIUS
     print 10,'RADIUS=',r0

  else if (rho0 > 0.0) then

     ! Use local rho

     x(1)  = rho0
     x_vec = EXPRO_rho

     ! RADIUS
     call cub_spline(x_vec,EXPRO_rmin/a,EXPRO_n_exp,x,y,1)
     print 10,'RADIUS=',y(1)

     r0 = y(1)

  else 

     ! Use local psi_N

     x(1)  = psi0*EXPRO_polflux(EXPRO_n_exp)
     x_vec = EXPRO_polflux

     ! RADIUS
     call cub_spline(x_vec,EXPRO_rmin/a,EXPRO_n_exp,x,y,1)
     print 10,'RADIUS=',y(1)

     r0 = y(1)

  endif

  GEO_rmin_in = r0 

  ! ASPECT_RATIO
  call cub_spline(x_vec,EXPRO_rmaj/a,EXPRO_n_exp,x,y,1)
  print 10,'ASPECT_RATIO=',y(1)
  ar = y(1)

  GEO_rmaj_in = y(1)

  ! SHIFT
  call cub_spline(x_vec,EXPRO_drmaj,EXPRO_n_exp,x,y,1)
  print 10,'SHIFT=',y(1)

  GEO_drmaj_in = y(1)

  ! ZMAG
  call cub_spline(x_vec,EXPRO_zmag/a,EXPRO_n_exp,x,y,1)
  print 10,'ZMAG=',y(1)

  GEO_zmag_in = y(1)

  ! DZMAG
  call cub_spline(x_vec,EXPRO_dzmag,EXPRO_n_exp,x,y,1)
  print 10,'DZMAG=',y(1)

  GEO_dzmag_in = y(1)

  ! SHEAR
  call cub_spline(x_vec,EXPRO_s,EXPRO_n_exp,x,y,1)
  print 10,'SHEAR=',y(1)
  shear = y(1)

  GEO_s_in = shear

  ! SAFETY_FACTOR
  call cub_spline(x_vec,EXPRO_q,EXPRO_n_exp,x,y,1)
  print 10,'SAFETY_FACTOR=',y(1)
  sf = y(1)

  GEO_q_in = y(1)

  ! KAPPA
  call cub_spline(x_vec,EXPRO_kappa,EXPRO_n_exp,x,y,1)
  print 10,'KAPPA=',y(1)

  GEO_kappa_in = y(1)

  ! S_KAPPA
  call cub_spline(x_vec,EXPRO_skappa,EXPRO_n_exp,x,y,1)
  print 10,'S_KAPPA=',y(1)

  GEO_s_kappa_in = y(1)

  ! DELTA
  call cub_spline(x_vec,EXPRO_delta,EXPRO_n_exp,x,y,1)
  print 10,'DELTA=',y(1)

  GEO_delta_in = y(1)

  ! S_DELTA
  call cub_spline(x_vec,EXPRO_sdelta,EXPRO_n_exp,x,y,1)
  print 10,'S_DELTA=',y(1)

  GEO_s_delta_in = y(1)

  ! ZETA
  call cub_spline(x_vec,EXPRO_zeta,EXPRO_n_exp,x,y,1)
  print 10,'ZETA=',y(1)

  GEO_zeta_in = y(1)

  ! S_ZETA
  call cub_spline(x_vec,EXPRO_szeta,EXPRO_n_exp,x,y,1)
  print 10,'S_ZETA=',y(1)

  GEO_s_zeta_in = y(1)

  ! TI_OVER_TE
  call cub_spline(x_vec,EXPRO_ti(1,:)/EXPRO_te,EXPRO_n_exp,x,y,1)
  print 10,'TI_OVER_TE=',y(1)

  ! TI_OVER_TE_2
  call cub_spline(x_vec,EXPRO_ti(2,:)/EXPRO_te,EXPRO_n_exp,x,y,1)
  print 10,'TI_OVER_TE_2=',y(1)

  ! NI_OVER_NE
  call cub_spline(x_vec,EXPRO_ni(1,:)/EXPRO_ne,EXPRO_n_exp,x,y,1)
  print 10,'NI_OVER_NE=',y(1)

  ! NI_OVER_NE_2
  call cub_spline(x_vec,EXPRO_ni(2,:)/EXPRO_ne,EXPRO_n_exp,x,y,1)
  print 10,'NI_OVER_NE_2=',y(1)

  ! DLNNDR
  call cub_spline(x_vec,a*EXPRO_dlnnidr(1,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNNDR=',y(1)

  ! DLNNDR_2
  call cub_spline(x_vec,a*EXPRO_dlnnidr(2,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNNDR_2=',y(1)

  ! DLNNDR_ELECTRON
  call cub_spline(x_vec,a*EXPRO_dlnnedr,EXPRO_n_exp,x,y,1)
  print 10,'DLNNDR_ELECTRON=',y(1)

  ! DLNTDR
  call cub_spline(x_vec,a*EXPRO_dlntidr(1,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNTDR=',y(1)

  ! DLNTDR_2
  call cub_spline(x_vec,a*EXPRO_dlntidr(2,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNTDR_2=',y(1)

  ! DLNTDR_ELECTRON
  call cub_spline(x_vec,a*EXPRO_dlntedr,EXPRO_n_exp,x,y,1)
  print 10,'DLNTDR_ELECTRON=',y(1)

  ! GAMMA_E
  call cub_spline(x_vec,EXPRO_gamma_e*a/EXPRO_cs,EXPRO_n_exp,x,y,1)
  print 10,'GAMMA_E=',y(1)

  ! RHO_STAR
  call cub_spline(x_vec,EXPRO_rhos/a,EXPRO_n_exp,x,y,1)
  print 10,'RHO_STAR=',y(1)

  !--------------------------------------------------------
  ! Model geometry output
  GEO_ntheta_in = 2001

  GEO_signb_in=1.0 
  GEO_beta_star_in=0.0

  GEO_model_in = 0
  call GEO_alloc(1)
  call GEO_do()
  call GEO_write('out.locpargen.geo',2)
  call GEO_alloc(0)
  !--------------------------------------------------------

  print *
  print '(a)','INFO: (locpargen) Additional quantities:'
  print *

  !---------------------------------------
  ! Some added physical quantities
  call cub_spline(x_vec,EXPRO_rho,EXPRO_n_exp,x,y,1)
  print 10,'rho         : ',y(1)
  call cub_spline(x_vec,EXPRO_rmin,EXPRO_n_exp,x,y,1)
  print 10,'rmin [m]    : ',y(1)
  call cub_spline(x_vec,EXPRO_polflux,EXPRO_n_exp,x,y,1)
  if (abs(EXPRO_polflux(EXPRO_n_exp)) > 0.0) then
     print 10,'psi_N       : ',y(1)/EXPRO_polflux(EXPRO_n_exp)
  else
     print '(a)','psi_N       : UNAVAILABLE'
  endif
  call cub_spline(x_vec,EXPRO_bunit,EXPRO_n_exp,x,y,1)
  print 10,'B_unit [T]  : ',y(1)
  call cub_spline(x_vec,EXPRO_cs,EXPRO_n_exp,x,y,1)
  print 10,'c_s   [m/s] : ',y(1)
  call cub_spline(x_vec,EXPRO_rhos,EXPRO_n_exp,x,y,1)
  print 10,'rhos   [m]  : ',y(1)
  call cub_spline(x_vec,EXPRO_vol,EXPRO_n_exp,x,y,1)
  print 10,'vol  [m^3]  : ',y(1)
  print 10,'vol/a^3     : ',y(1)/a**3
  call cub_spline(x_vec,sqrt(EXPRO_thetascale),EXPRO_n_exp,x,y,1)
  print 10,'thetascale[-]: ',y(1)
  !---------------------------------------

  !------------------------------------------------------------
  ! Create input.geo with local parameters for general geometry
  !
  if (EXPRO_nfourier > 0) then  

     allocate(geo_p(8,0:EXPRO_nfourier,EXPRO_n_exp))

     geo_p(1:4,:,:) = EXPRO_geo(:,:,:)/EXPRO_rmin(EXPRO_n_exp)
     geo_p(5:8,:,:) = EXPRO_dgeo(:,:,:)

     open(unit=1,file='input.geo',status='replace')
     write(1,'(a)') '# input.geo'
     write(1,'(a)') '#'
     write(1,'(a)') '# See https://fusion.gat.com/theory/input.geo for complete documentation.'
     write(1,'(a)') '#'
     write(1,'(a,f10.6)') '# NOTE: Derived from input.profiles.geo at r/a=',x
     write(1,'(a)') '# Lengths normalized to a' 
     write(1,'(a,f10.6)') '# ASPECT_RATIO=',ar
     write(1,'(a,f10.6)') '# SAFETY_FACTOR=',sf
     write(1,'(a,f10.6)') '# SHEAR=',shear
     write(1,'(a,i3)') '# BTCCW=',-EXPRO_signb
     write(1,'(a,i3)') '# IPCCW=',-EXPRO_signb*EXPRO_signq
     write(1,'(a)') '#'
     write(1,'(a)') '# File format:'
     write(1,'(a)') '#-------------------'
     write(1,'(a)') '# nfourier'
     write(1,'(a)') '# a[8,0:nfourier]'    
     write(1,'(a)') '#-------------------'
     write(1,'(a)') '#'
     write(1,*) EXPRO_nfourier

     GEO_model_in = 1
     GEO_nfourier_in = EXPRO_nfourier
     call GEO_alloc(1)

     do j2=0,EXPRO_nfourier
        do j1=1,8
           call cub_spline(x_vec,geo_p(j1,j2,:),EXPRO_n_exp,x,y,1)
           write(1,'(1pe20.13)') y(1)
           GEO_fourier_in(j1,j2) = y(1)
        enddo
     enddo

     call GEO_do()
     call GEO_write('out.locpargen.geo.2',2)
     call GEO_alloc(0)

     print *
     print '(a)','INFO: (locpargen) Wrote input.geo.'

     close(1)

  endif
  !------------------------------------------------------------

  call EXPRO_alloc('./',0) 

10 format(a,sp,1pe12.5)

end program locpargen
