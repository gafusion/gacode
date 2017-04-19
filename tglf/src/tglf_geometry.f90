SUBROUTINE get_xgrid_functions
  !***********************************************************
  !
  !***********************************************************
  USE tglf_global
  !
  new_width=.FALSE.
  ! write(*,*)"get_xgrid_functions"
  !
  if(igeo.eq.0)call xgrid_functions_sa
  if(igeo.ne.0)call xgrid_functions_geo
  !
END SUBROUTINE get_xgrid_functions
!
!
SUBROUTINE xgrid_functions_sa
  !***********************************************************
  !
  !***********************************************************
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_hermite
  USE tglf_xgrid
  !
  IMPLICIT NONE
  INTEGER :: i
  REAL :: thx,dthx,sn,cn,eps,Rx,Rx1,Rx2
  REAL :: kyi,wE,a0,vexb_shear_kx0
  !
  ! debug
  ! write(*,*)"shat_sa=",shat_sa,"alpha_sa=",alpha_sa
  ! write(*,*)"theta0_sa=",theta0_sa,"xwell_sa=",xwell_sa
  ! write(*,*)"rmin_sa=",rmin_sa,"rmaj_sa=",rmaj_sa
  ! write(*,*)"q_sa=",q_sa,"b_model_sa=",b_model_sa
  ! debug
  !
  ! set the units used by the tglf equations
  R_unit=rmaj_sa
  q_unit=q_sa
  !
  ! set global input values
  !
  rmin_input = rmin_sa
  Rmaj_input = rmaj_sa
  !
  ! fill the x-grid eikonal function arrays wdx and b0x
  eps = rmin_sa/rmaj_sa
  !
  ! generalized quench rule kx0 shift
  !
  !EPS2011 midplane_shear = shat_sa - alpha_sa
  !EPS2011 sign_kx0=1.0
  kx0 = 0.0
  kx0_e = 0.0
  if(alpha_quench_in.eq.0.0.and.gamma_reference_kx0(1).ne.0.0)then
     vexb_shear_kx0 = alpha_e_in*vexb_shear_s
     kyi = ky*vs(2)*mass(2)/ABS(zs(2))
     wE = MIN(kyi/0.3,1.0)*vexb_shear_kx0/gamma_reference_kx0(1)
     ! write(*,*)"wE=",wE
     kx0_e = -(0.36*vexb_shear_kx0/gamma_reference_kx0(1) + 0.38*wE*TANH((0.69*wE)**6))
     if(sat_rule_in.eq.1)kx0_e = -(0.53*vexb_shear_kx0/gamma_reference_kx0(1) + 0.25*wE*TANH((0.69*wE)**6))
     a0 = 1.3
     if(sat_rule_in.eq.1)a0=1.45
     if(ABS(kx0_e).gt.a0)kx0_e = a0*kx0_e/ABS(kx0_e)
!     a0 = alpha_e_in*2.0
!     if(alpha_e_in.ne.0.0)then
!        kx0_e = a0*TANH(kx0_e/a0)
!     else
!        kx0_e = 0.0
!     endif
     !
     ! EPS2011 kx0 = alpha_kx_e_in*0.19*TANH(vexb_shear_s*rmaj_sa/vs(2))*kyi*kyi/(kyi*kyi+0.001)/ky
     ! EPS2011 kx0 = kx0 -alpha_kx_p_in*sign_Bt_in*TANH((0.26/3.0)*vpar_shear_in(2)*rmaj_sa/vs(2)) &
     ! EPS2011 *(0.06*kyi*kyi/(kyi*kyi+0.001)+0.25*TANH((1.9*kyi)**3))/ky
     ! EPS2011 if(midplane_shear.gt.0.0)then
     ! EPS2011 kx0 = kx0 -alpha_kx_p_in*sign_Bt_in*TANH((0.26/3.0)*vpar_shear_in(2)*Rmaj_sa/vs(2)) &
     ! EPS2011 *midplane_shear*(1.215*EXP(-(5.4*kyi)**3)+1.274**TANH((1.7*kyi)**2))/(1.0 + kyi**2)
     ! EPS2011 else
     ! EPS2011 kx0 = kx0 -alpha_kx_p_in*sign_Bt_in*TANH((0.26/3.0)*vpar_shear_in(2)*Rmaj_sa/vs(2)) &
     ! EPS2011 *midplane_shear*(1.132**TANH((1.4*kyi)**2))/(1.0 + kyi**2)
     ! EPS2011 endif
     ! EPS2011 kx0 = kx0 + alpha_kx_n_in*0.19*TANH(vns_shear_in(2)*rmaj_sa/vs(2))*kyi*kyi/(kyi*kyi+0.001)/ky
     ! EPS2011 kx0 = kx0 - alpha_kx_t_in*0.19*TANH(vts_shear_in(2)*rmaj_sa/vs(2))*kyi*kyi/(kyi*kyi+0.001)/ky
     ! EPS2011 if(kx0.lt.0.0)sign_kx0=-1.0
     kx0 = sign_Bt_in*kx0_e ! this is here to cancel the sign_Bt_in factor in kxx below
  endif
 !
  do i=1,nx
     thx = width_in*x(i)
     sn = sin(thx)
     cn = cos(thx)
     Rx = 1.0 + eps*cn
     Bx(i) = 1.0/Rx
     kxx(i) = -shat_sa*(thx-theta0_sa) + alpha_sa*sn - kx0
     kxx(i) = sign_Bt_in*kxx(i)
     ! wdx(i) = -xwell_sa*MIN(1.0,alpha_sa)+ &
     ! cn+sn*(shat_sa*(thx-theta0_sa) - alpha_sa*sn)
     wdx(i) = -xwell_sa*MIN(1.0,alpha_sa)+ cn - sn*kxx(i)
     ! b0x(i) = 1.0+(shat_sa*(thx-theta0_sa) - alpha_sa*sn)**2
     b0x(i) = 1.0+(kxx(i))**2
     b2x(i) = 1.0
     wdpx(i) = 0.0   ! not used for s-alpha
     if(b_model_sa.eq.1)then
        ! put B dependence into k_per**2
        b2x(i) = Bx(i)**2
     endif
     if(b_model_sa.eq.2)then
        ! put 1/R(theta) factor into wd
        wdx(i) = wdx(i)/Rx
        ! put B dependence into k_per**2
        b2x(i) = Bx(i)**2
     endif
     ! debug
     ! write(*,*)i,thx,Bx(i),wdx(i),b0x(i)
     !
     ! momentum equation stress projection coefficients
     !
     cx_tor_par(i) = rmaj_sa*Rx*sign_Bt_in
     cx_tor_per(i) = -rmin_sa/q_sa
     cx_par_par(i) = Bx(i)
  enddo
  !
  ! compute the effective trapped fraction
  !
  call get_ft_sa
  !
  ! compute the flux surface averages
  !
  B2_ave_out = 0.0
  R2_ave_out = 0.0
  RBt_ave_out = rmaj_sa
  B_ave_out = 0.0
  Bt_ave_out = 0.0
  dthx = pi_2/REAL(2*nx)
  thx = 0.0
  do i=1,2*nx
     Rx1 = 1.0 + eps*COS(thx)
     thx = thx + dthx
     Rx2 = 1.0 + eps*COS(thx)
     R2_ave_out = R2_ave_out + dthx*(Rx1**2 + Rx2**2)/2.0
     B2_ave_out = B2_ave_out + dthx*(0.5/Rx1**2 + 0.5/Rx2**2)
     B_ave_out = B_ave_out + dthx*(0.5/Rx1 +0.5/Rx2)
     Bt_ave_out = Bt_ave_out + dthx*(0.5*Rx1 + 0.5*Rx2)
  enddo
  R2_ave_out = (R2_ave_out*rmaj_sa**2)/pi_2
  B2_ave_out = B2_ave_out/pi_2
  B_ave_out = B_ave_out/pi_2
  Bt_ave_out = Bt_ave_out/pi_2
  Grad_r_ave_out = 1.0
  !
  ! poloidal magnetic field at outboard midplane
  !
  Bp0_out = rmin_sa/(q_sa*(rmaj_sa+rmin_sa))
  !
END SUBROUTINE xgrid_functions_sa
!
SUBROUTINE xgrid_functions_geo
  !******************************************************************************!************************
  !
  ! PURPOSE: compute the geometric coefficients on the x-grid
  !
  !
  !******************************************************************************!************************
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_hermite
  USE tglf_xgrid
  USE tglf_sgrid
  IMPLICIT NONE
  !
  INTEGER :: i,m
  INTEGER :: m1,m2
  REAL :: y_x
  REAL :: Ly
  REAL :: thx
  REAL :: sign_theta,loops
  REAL :: dkxky1,dkxky2
  REAL :: wd1,wd2,wdp1,wdp2
  REAL :: b1,b2
  REAL :: y1,y2
  REAL :: kxx1,kxx2
  REAL :: cxtorper1,cxtorper2
  REAL :: B2x1,B2x2,R2x1,R2x2,norm_ave,dlp
  REAL :: kyi
  REAL :: wE,wd0,a0,vexb_shear_kx0
  !
  !
  ! find length along magnetic field y
  !
  y(0)=0.0
  ! pk_geo = 2 Bp/B = 2 ds/dy
  do m=1,ms
     ! y(m) = y(m-1)+0.5*(2.0/pk_geo(m)+2.0/pk_geo(m-1))*ds
     y(m) = y(m-1)+s_p(m)*ds*4.0/(pk_geo(m)+pk_geo(m-1))
     ! y(1)=y(0)+0.5*(2.0*ds/(0.25*(pk_geo(1)+pk_geo(ms-1))))
     ! do m=1,ms-1
     ! y(m+1)=y(m-1)+2.0*ds/(0.25*(pk_geo(m+1)+pk_geo(m-1)))
  enddo
  ! set the global units
  Ly=y(ms)
  ! R_unit = Rmaj_s*b_geo(0)/(qrat_geo(0)*(costheta_geo(0)+costheta_p_geo(0)))
  R_unit = Rmaj_s*b_geo(0)/(qrat_geo(0)*costheta_geo(0))
  q_unit = Ly/(pi_2*R_unit)
  ! midplane effective shear: reduces to s-alpha in shifted circle
  ! note: S_prime(0)=0.0, S_prime(ms)=-2 pi q_prime, y(0)=0.0, y(ms)=Ly
  ! midplane shear is average of left and right y-derivatives at midplane for general geometry
  midplane_shear = -(Ly/pi_2)*((rmin_s/q_s)**2) &
       *0.5*(S_prime(1)/y(1)+(S_prime(ms)-S_prime(ms-1))/(y(ms)-y(ms-1)))
  ! write(*,*)"midplane_shear = ",midplane_shear
  midplane_shear = midplane_shear + 0.11
  ! save f for output
  RBt_ave_out = f/B_unit
  !
  ! write(*,*)"qrat_geo(0)=",qrat_geo(0),"b_geo(0)=",b_geo(0)
  ! write(*,*)"costheta_geo(0)=",costheta_geo(0)
  ! write(*,*)"sintheta_geo(0)=",sintheta_geo(0)
  ! write(*,*)"Ly = ",Ly,"R_unit = ",R_unit,"Rmaj_s =",Rmaj_s
  ! write(*,*)"q_unit=",q_unit
  ! open(2,file='y_s.dat')
  ! do m=0,ms
  ! write(2,*)m,pi_2*y(m)/Ly
  ! enddo
  ! close(2)
  !
  ! generalized quench rule kx0 shift
  !
  kx0 = kx0_loc/ky ! note that kx0 is kx/ky
  kx0_e=0.0
  kx0_p=0.0
  !EPS2011 sign_kx0=1.0
  if(alpha_quench_in.eq.0.0.and.gamma_reference_kx0(1).ne.0.0)then
     vexb_shear_kx0 = alpha_e_in*vexb_shear_s
     kyi = ky*vs(2)*mass(2)/ABS(zs(2))
     wE=0.0
     wd0 = ABS(ky/Rmaj_s)
     kx0_factor = ABS(b_geo(0)/qrat_geo(0)**2)
     kx0_factor = 1.0+0.40*(kx0_factor-1.0)**2
     ! write(*,*)"kx0_factor=",kx0_factor
     wE = kx0_factor*MIN(kyi/0.3,1.0)*vexb_shear_kx0/gamma_reference_kx0(1)
     ! write(*,*)"wE=",wE
     ! kx0_factor = alpha_kx_p_in
     !! kx0_e = -(0.36*vexb_shear_s/gamma_reference_kx0(1) + 0.29*wE*TANH((0.71*wE)**6))
     kx0_e = -(0.36*vexb_shear_kx0/gamma_reference_kx0(1) + 0.38*wE*TANH((0.69*wE)**6))
     if(sat_rule_in.eq.1)kx0_e = -(0.53*vexb_shear_kx0/gamma_reference_kx0(1) + 0.25*wE*TANH((0.69*wE)**6))
!     a0 = alpha_e_in*2.0
!     if(alpha_e_in.ne.0.0)then
!        kx0_e = a0*TANH(kx0_e/a0)
!     else
!        kx0_e = 0.0
!     endif
     a0 = 1.3
     if(sat_rule_in.eq.1)a0=1.45
     if(ABS(kx0_e).gt.a0)kx0_e = a0*kx0_e/ABS(kx0_e)
     !
     !APS2010 kx0 = alpha_kx_e_in*0.19*TANH(vexb_shear_s*Rmaj_s/vs(2))*kyi*kyi/(kyi*kyi+0.001)/ky
     !EPS2011 kx0 = kx0 -alpha_kx_p_in*sign_Bt_in*TANH((0.26/3.0)*vpar_shear_in(2)*Rmaj_s/vs(2)) &
     !EPS2011 *(0.06*kyi*kyi/(kyi*kyi+0.001)+0.25*TANH((1.9*kyi)**3))/ky
     !EPS2011 if(midplane_shear.gt.0.0)then
     !EPS2011 kx0_p = -alpha_kx_p_in*sign_Bt_in*TANH((0.26/3.0)*vpar_shear_in(2)*Rmaj_s/vs(2)) &
     !EPS2011 *midplane_shear*(1.43*EXP(-(5.4*kyi)**3)+1.50*TANH((1.7*kyi)**2))/(1.0 + kyi**2)
     !EPS2011 else
     !EPS2011 kx0_p = -alpha_kx_p_in*sign_Bt_in*TANH((0.26/3.0)*vpar_shear_in(2)*Rmaj_s/vs(2)) &
     !EPS2011 *midplane_shear*(2.13*TANH((1.4*kyi)**2))/(1.0 + kyi**2)
     !EPS2011 endif
     ! kx0 = kx0 + alpha_kx_n_in*0.19*TANH(vns_shear_in(2)*Rmaj_s/vs(2))*kyi*kyi/(kyi*kyi+0.001)/ky
     ! kx0 = kx0 - alpha_kx_t_in*0.19*TANH(vts_shear_in(2)*Rmaj_s/vs(2))*kyi*kyi/(kyi*kyi+0.001)/ky
     !
     !ESP2011 if(kx0.lt.0.0)sign_kx0=-1.0
     !EPS2011 kx0 = kx0_e + kx0_p
     ! write(*,*)ky,"kx0_e",kx0_e,"kx0_p=",kx0_p,"kx0=",kx0
     kx0 = sign_Bt_in*kx0_e ! this is here to cancel the sign_Bt_in factor in kxx below
  endif
  ! kx0 = alpha_kx_e_in
  ! write(*,*)"1/qrat=",1.0/qrat_geo(0)
  !
  !*************************************************************
  ! begin calculation of wdx and b0x
  !*************************************************************
  ! compute drift wdx and pependicular wavenumber squared b0x
  ! at the Hermite nodes x(i)
  ! thx is the ballooning angle = 2 pi y/Ly
  ! x is the argument of the Hermite basis functions = thx/width_in
  !
  do i=1,nx
     thx = width_in*x(i)
     sign_theta=1.0
     if(thx.lt.0.0) sign_theta=-1.0
     !
     ! use quasi-periodic property of S_prime to evaluate it when thx > 2pi
     ! S_prime(0)=0
     ! S_prime(t+2pi*m) = m*S_prime(2pi) + S_prime(t)
     !
     loops = REAL(INT(ABS(thx/pi_2)))
     y_x = Ly*(ABS(thx) - loops*pi_2)/pi_2
     if(thx.lt.0.0)then
        !
        ! not always up/down symmetric so can't just use symmetry
        ! for negative theta = -(t+m*2pi) use
        ! S_prime(-t-m*2pi) = -(m+1)*S_prime(2pi)+S_prime(2pi-t)
        !
        y_x=Ly-y_x
        loops=loops+1.0
     endif
     do m=1,ms
        if(y(m).ge.y_x) exit
     enddo
     ! write(*,*)"exit at",m,ms,y_x,y(ms)
     if(m.gt.ms)m=ms
     m1=m-1
     m2=m
     !
     ! dkxky is the offset for S_prime due to the number of loops
     !
     dkxky1 = sign_theta*loops*S_prime(ms)
     y1=y(m1)
     dkxky2 = sign_theta*loops*S_prime(ms)
     y2=y(m2)
     ! write(*,*)"check interpolation",m1,m2
     ! write(*,*)"y=",y1,y2
     ! write(*,*)S_prime(m1)+dkxky1,S_prime(m2)+dkxky2
     !
     ! note that costheta_geo and sintheta_geo are periodic so
     ! we can use f(-t-m*2pi) = f(2pi-t) if f(0)=f(2pi) and 0<t<2pi
     !
     ! intepolate kxx
     kxx1 = (kx_factor(m1)*(S_prime(m1)+dkxky1)-kx0*b_geo(m1)/qrat_geo(m1)**2)*qrat_geo(m1)/b_geo(m1)
     kxx2 = (kx_factor(m2)*(S_prime(m2)+dkxky2)-kx0*b_geo(m2)/qrat_geo(m2)**2)*qrat_geo(m2)/b_geo(m2)
     kxx(i) = kxx1 +(kxx2-kxx1)*(y_x-y1)/(y2-y1)
     kxx(i) = sign_Bt_in*kxx(i)
     ! interpolate wdx
     wd1 = (qrat_geo(m1)/b_geo(m1))*(costheta_geo(m1) &
          +(kx_factor(m1)*(S_prime(m1)+dkxky1)-kx0*b_geo(m1)/qrat_geo(m1)**2)*sintheta_geo(m1))
     wd2 = (qrat_geo(m2)/b_geo(m2))*(costheta_geo(m2) &
          +(kx_factor(m2)*(S_prime(m2)+dkxky2)-kx0*b_geo(m2)/qrat_geo(m2)**2)*sintheta_geo(m2))
     ! write(*,*)"wd1,,wd2=",wd1,wd2
     wdx(i) = wd1 +(wd2-wd1)*(y_x-y1)/(y2-y1)
     wdx(i) = (R_unit/Rmaj_s)*wdx(i)
     ! write(*,*)i,"wdx = ",wdx(i),y_x,x(i),y1,y2
     ! interpolate wdpx
     wdp1 = (qrat_geo(m1)/b_geo(m1))*costheta_p_geo(m1)
     wdp2 = (qrat_geo(m2)/b_geo(m2))*costheta_p_geo(m2)
     wdpx(i) = wdp1 + (wdp1-wdp2)*(y_x-y1)/(y2-y1)
     wdpx(i) = (R_unit/Rmaj_s)*wdpx(i)
     ! interpolate b2x = b_geo**2
     b1 = b_geo(m1)**2
     b2 = b_geo(m2)**2
     b2x(i) =  b1 +(b2-b1)*(y_x-y1)/(y2-y1)
     ! interpolate b0x
     b1 = (1.0+(kx_factor(m1)*(S_prime(m1)+dkxky1)-kx0*b_geo(m1)/qrat_geo(m1)**2)**2) &
          *qrat_geo(m1)**2
     b2 = (1.0+(kx_factor(m2)*(S_prime(m2)+dkxky2)-kx0*b_geo(m2)/qrat_geo(m2)**2)**2) &
          *qrat_geo(m2)**2
     ! write(*,*)"b1,b2,b3,b4=",b1,b2,b3,b4
     b0x(i) =  b1 +(b2-b1)*(y_x-y1)/(y2-y1)
     if(b0x(i).lt.0.0)then
        write(*,*)"interpolation error b0x < 0",i,b0x(i),b1,b2
        b0x(i)=(b1+b2)/2.0
     endif
     !
     ! interpolate viscous stress projection coefficients
     !
     cxtorper1 = -R(m1)*Bp(m1)/b_geo(m1)
     cxtorper2 = -R(m2)*Bp(m2)/b_geo(m2)
     cx_tor_par(i) = f/b_geo(m1) + (f/b_geo(m2)-f/b_geo(m1))*(y_x-y1)/(y2-y1)
     cx_tor_par(i) = sign_Bt_in*cx_tor_par(i)
     cx_tor_per(i) = cxtorper1 + (cxtorper2-cxtorper1)*(y_x-y1)/(y2-y1)
     cx_par_par(i) = b_geo(m1) + (b_geo(m2)-b_geo(m1))*(y_x-y1)/(y2-y1)
     !
     ! debug
     ! write(*,*)"b0x = ",b0x(i)
     ! write(*,*)" y_x =",y_x," m = ",m
     ! write(*,*)i," loops = ",loops
  enddo
  ! write(*,*)"ave_shear =",ave_shear,midplane_shear
  !
  ! compute flux surface averages
  !
  B2_ave_out = 0.0
  R2_ave_out = 0.0
  B_ave_out = 0.0
  Bt_ave_out = 0.0
  Grad_r_ave_out = 0.0
  norm_ave=0.0
  do i=1,ms
     dlp = s_p(i)*ds*(0.5/Bp(i)+0.5/Bp(i-1))
     norm_ave = norm_ave + dlp
     B2x1 = b_geo(i-1)**2
     B2x2 = b_geo(i)**2
     B2_ave_out = B2_ave_out + dlp*(B2x1+B2x2)/2.0
     R2x1 = R(i-1)**2
     R2x2 = R(i)**2
     R2_ave_out = R2_ave_out + dlp*(R2x1+R2x2)/2.0
     B_ave_out = B_ave_out + dlp*(b_geo(i-1)+b_geo(i))/2.0
     Bt_ave_out = Bt_ave_out + dlp*(f/b_geo(i-1)+f/b_geo(i))/(2.0*Rmaj_s)
     Grad_r_ave_out = Grad_r_ave_out + dlp*0.5*((R(i-1)*Bp(i-1))+(R(i)*Bp(i)))*(q_s/rmin_s)
  enddo
  R2_ave_out = R2_ave_out/norm_ave
  B2_ave_out = B2_ave_out/norm_ave
  B2_ave_out = B2_ave_out/B_unit**2
  B_ave_out = B_ave_out/norm_ave
  Bt_ave_out = Bt_ave_out/norm_ave
  Grad_r_ave_out = Grad_r_ave_out/norm_ave
  !
  ! poloidal magnetic field on outboard midplane
  !
  Bp0_out = Bp(0)/B_unit
  !
  ! write(*,*)"R2_ave_out=",R2_ave_out
  ! write(*,*)"B2_ave_out=",B2_ave_out
  ! write(*,*)"B_ave_out=",B_ave_out
  ! write(*,*)"Bt_ave_out=",Bt_ave_out
  !
  ! do m=0,ms
  ! write(*,*)m,s_prime(ms-m),s_prime(ms)-s_prime(m)
  ! enddo
  !*************************************************************
  ! end of calculation of wdx and b0x
  !*************************************************************
  !
  ! compute the effective trapped fraction
  call get_ft_geo
  ! write(*,*)"ft_geo=",ft
  !
END SUBROUTINE xgrid_functions_geo
!
!
SUBROUTINE get_ft_sa
  !
  ! shifted circle version of get_ft
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_hermite
  !
  IMPLICIT NONE
  INTEGER :: i
  REAL :: norm,ww
  REAL :: eps,theta_max
  REAL :: cn,thx,ftx,Bmax,Bmin
  REAL :: Rx,Bx
  !
  ! compute pitch angle at bounce average boundary
  !
  eps = rmin_sa/rmaj_sa
  ! vshear_eff = 2.0*(ky/R_unit)*ABS(vpar_in)*sqrt_two*width_in
  ! theta_eff = width_in/MAX(1.0-vshear_eff,0.1)
  ! theta_max = theta_trapped_in*theta_eff*pi/sqrt_two
  theta_max = theta_trapped_in*width_in*pi/sqrt_two
  if(theta_max.gt.pi)theta_max=pi
  Bmax = 1.0/(1.0 + eps*COS(theta_max))
  ! write(*,*)"theta_max = ",theta_max," Bmax = ",Bmax
  Bmin = 1.0/(1.0 + eps)
  modB_min = Bmin
  modB_test = 0.5*(Bmax + Bmin)/Bmin
  !
  if(ft_model_sa.eq.0)then
     ! Gauss-Chebyshev integration of trappped fraction
     ! over the bounce angle weighted by gaussian wavefunction
     ft = 0.0
     norm = 0.0
     do i=1,nx
        ww = COS(pi*REAL(2*i-1)/REAL(2*nx))
        thx = ww*theta_max
        cn = COS(thx)
        Rx = (1.0 + eps*cn)
        Bx = 1.0/Rx
        ftx = 1.0 - Bx/Bmax
        ftx = SQRT(ftx)
        ww = SQRT(1.0 - ww**2)*EXP(-(thx/width_in)**2)
        ft = ft + ww*ftx
        norm = norm + ww
     enddo
     ft = ft/norm
     ! write(*,*)"ft = ",ft
  endif
  !
  if(ft_model_sa.eq.1)ft = SQRT(1.0-Bmin/Bmax)
  !
  if(ft_model_sa.eq.2)then
     ! Gauss-Chebyshev integration of trappped fraction
     ! over the bounce angle
     ft = 0.0
     norm = 0.0
     do i=1,nx
        ww = COS(pi*REAL(2*i-1)/REAL(2*nx))
        thx = ww*theta_max
        cn = COS(thx)
        Rx = (1.0 + eps*cn)
        Bx = 1.0/Rx
        ftx = 1.0 - Bx/Bmax
        ftx = SQRT(ftx)
        ww = SQRT(1.0 - ww**2)
        ft = ft + ww*ftx
        norm = norm + ww
     enddo
     ft = ft/norm
  endif
  !
  if(ft_model_sa.eq.3)then
     ! hermite integration averaged ft with Gaussian envelope of wavefunction
     ft = 0.0
     do i=1,nx
        ww = wx(i)*h(1,i)*h(1,i)
        thx = x(i)*width_in + theta0_sa
        Rx = 1.0 + eps*COS(thx)
        Bx = 1.0/Rx
        ftx = SQRT(MAX(1.0 - Bx/Bmax,0.0))
        ! write(*,*)i,"thx = ",thx,"ftx=",ftx
        ft = ft + ww*ftx
     enddo
     ! write(*,*)"ft = ",ft,"ft0=",SQRT(1.0-Bmin/Bmax)
  endif
  ! write(*,*)"ft = ",ft,"ft_model_sa =",ft_model_sa
  !
END SUBROUTINE get_ft_sa
!
SUBROUTINE get_ft_geo
  !
  ! general geometry version of get_ft
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_sgrid
  !
  IMPLICIT NONE
  !
  INTEGER,PARAMETER :: nb_grid=25
  INTEGER :: i,j,m,m_max,m_min,j_max
  INTEGER :: pm(2,0:nb_grid),qm
  REAL :: Bmax,Bmin,By(0:nb_grid),delta_y(0:nb_grid)
  REAL :: Ly
  REAL :: B_bounce,kpar
  REAL :: db,test1,test2,bounce_y
  !
  !*************************************************************
  ! begin trapped fraction model
  !*************************************************************
  ! find global Bmax and Bmin
  !
  ! test case for debug
  ! do m=0,ms/2
  ! b_geo(m)= 1.0 +2.0*(REAL(ABS(m))/REAL(ms/2))**2
  ! b_geo(ms-m)=b_geo(m) +0.01*REAL(ABS(ms/2-m)*m)/REAL(ms/2)
  ! b_geo(m)=b_geo(m) -0.01*REAL(ABS(ms/2-m)*m)/REAL(ms/2)
  ! enddo
  !
  Ly=y(ms)
  Bmax = b_geo(0)
  Bmin = b_geo(0)
  m_max = 0
  m_min = 0
  do m=1,ms
     if(b_geo(m).gt.Bmax)then
        Bmax = b_geo(m)
        m_max = m
     endif
     if(b_geo(m).lt.Bmin)then
        Bmin = b_geo(m)
        m_min = m
     endif
  enddo
  ! write(*,*)"global Bmax = ",Bmax,"at m =",m_max
  ! write(*,*)"global Bmin = ",Bmin,"at m =",m_min
  ! write(*,*)"Bmax/Bmin = ",Bmax/Bmin
  !
  ! make a table of Bmin=< By <= Bmax
  !
  By(0)=Bmin
  dB = (Bmax - Bmin)/REAL(nb_grid)
  do i=1,nb_grid
     By(i) = By(i-1) + db
     ! write(*,*)i,"By=",By(i)
  enddo

  !
  ! find pairs of m's at the same B starting at Bmin and taking the farthest pair
  !
  pm(1,0)=m_min
  pm(2,0)=m_min
  qm=0
  do i=1,nb_grid-1
     ! go clockwise
     j_max = m_max - m_min
     if(j_max.lt.0)j_max=m_max+ms-m_min
     do m=1,j_max
        ! find the farthest gridpoint where b_geo=By
        j=m_min+m
        if(j.gt.ms)j=j-ms
        test1=b_geo(j-1)-By(i)
        test2=b_geo(j)-By(i)
        if(test1*test2.le.0.0)then
           if(ABS(test1).lt.ABS(test2))then
              qm=j-1
           else
              qm=j
           endif
           ! write(*,*)i,j,"qm =",qm,"test1 =",test1,"test2=",test2
        endif
     enddo
     pm(1,i)=qm
     ! go counterclockwise
     j_max = m_max - m_min
     if(j_max.lt.0)then
        j_max=ABS(j_max)
     else
        j_max=ms-j_max
     endif
     qm=0
     do m=1,j_max
        ! find the farthest gridpoint where b_geo=By
        j=m_min-m
        if(j.lt.0)j=j+ms
        test1=b_geo(j+1)-By(i)
        test2=b_geo(j)-By(i)
        if(test1*test2.le.0.0)then
           if(ABS(test1).lt.ABS(test2))then
              qm=j+1
           else
              qm=j
           endif
           ! write(*,*)i,j,"qm =",qm,"test1 =",test1,"test2=",test2
        endif
     enddo
     pm(2,i)=qm
  enddo
  pm(1,nb_grid)=m_max
  pm(2,nb_grid)=m_max
  do i=0,nb_grid
     if(pm(1,i).gt.ms.or.pm(2,i).gt.ms)then
        write(*,*)"error in get_ft_geo: pm out of bounds",pm(1,i),pm(2,i),ms
     endif
     ! write(*,*)i,"pm(1,i) =",pm(1,i)," pm(2,i)=",pm(2,i)
  enddo
  !
  ! now pm contains the pairs of m's with the same value of B
  !
  ! make a table of distances along field lines between the turning points for lookup
  !
  delta_y(0)=0.0
  do i=1,nb_grid
     if(y(pm(1,i)).gt.y(pm(2,i)))then
        delta_y(i) = y(pm(1,i))-y(pm(2,i))
     else
        delta_y(i)=Ly+y(pm(1,i))-y(pm(2,i))
     endif
     ! write(*,*)i,"delta_y =",delta_y(i)
  enddo
  !
  ! compute trapped fraction
  !
  kpar= pi_2/(Ly*sqrt_two*width_in)
  bounce_y = MIN(Ly,pi*theta_trapped_in/kpar)
  ! kpar = pi_2/(Ly*sqrt_two*width_in*theta_trapped_in) &
  ! +xnu_factor_in*xnuei_in*(Ly/pi)*SQRT(mass_in(1))
  ! bounce_y = MIN(Ly,pi/kpar)
  ! write(*,*)"bounce_y =",bounce_y,"kpar =",kpar,"Ly=",Ly
  ! write(*,*)"pi*theta_trapped/kpar =",pi*theta_trapped_in/kpar
  B_bounce = Bmax
  if(bounce_y.lt.Ly)then
     do i=1,nb_grid
        if(delta_y(i).gt.bounce_y)exit
     enddo
     B_bounce = By(i-1)+(By(i)-By(i-1))* &
          (bounce_y-delta_y(i-1))/(delta_y(i)-delta_y(i-1))
     ! write(*,*)i,"B_bounce =",B_bounce,Bmax
  endif
  ft = SQRT(1.0 - Bmin/B_bounce)
  modB_min = ABS(Bmin)
  modB_test = 0.5*(Bmax + Bmin)/Bmin
  ! write(*,*)"ft = ",ft
  !*************************************************************
  ! end of trapped fraction model
  !*************************************************************
END SUBROUTINE get_ft_geo
!
!*************************************************************
!
!---------------------------------------------------------------
! mercier_luc.f [various callers]
!
! PURPOSE: Compute ballooning mode eikonal form of poloidally varying terms
! in the TGLF equations using the mercier-luc local equilibrium expansion
!
! DEFINITIONS:
! Length is in units of rmin at the last closed
! flux surface (a_unit).
! NOTES:
! To follow what's going on in this routine, it is
! necessary to have the following papers at hand:
!
! (1) R.L. Miller, M.S. Chu, J.M. Greene, Y.R. Lin-liu
! and R.E. Waltz, 'Noncircular, finite aspect ratio,
! local equilibrium model', Phys. Plasmas 5 (1998) 973.
!
! (2) R.E. Waltz and R.L. Miller, 'Ion temperature gradient
! turbulence simulations and plasma flux surface shape',
! Phys. Plasmas 6 (1999) 4265.
!
! **See MILLER_inputs for control variables**
!
! REVISIONS:
! 13 Aug 02: jc
! Documentation after MILLER de-spaghetti.
! 18 Nov 03: jc
! Removed reset of s_delta and s_kappa when they
! are input as zero.
! 24 Nov 03: jc
! Entire code has been rewritten for greater
! efficiency and clarity. Documentation greatly
! improved as well.
! 13 Jan 10: gms
! added vprime and vpp and interchange stability calculation
!---------------------------------------------------------------

SUBROUTINE mercier_luc

  !-------------------------------------------
  ! the following must be defined from a previous call to one of the
  ! geometry routines miller_geo, fourier_geo,ELITE_geo and stored in tglf_sgrid:
  ! ms ! the number of points in the s-grid (flux surface contour)
  ! ds ! the arc length differential on a flux surface
  ! R(ms) ! the major radius on the s-grid
  ! Z(ms) ! the vertical coordinate on the s-grid
  ! Bp(ms) ! the poloidal magnetic field on the s-grid normalized to B_unit
  ! q_s = local flux surface safety factor
  ! q_prime_s = dq/dpsi
  ! p_prime_s = dp/dpsi
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_sgrid
  !
  IMPLICIT NONE
  INTEGER :: m,m1,m2,m3,m4
  !
  REAL :: sin_u(0:ms),r_curv(0:ms)
  REAL :: psi_x(0:ms)
  REAL :: Bt(0:ms),B(0:ms)
  REAL :: d_0(0:ms),d_p(0:ms),d_ffp(0:ms)
  REAL :: delta_s, ds2
  !
  REAL :: R_s, Z_s
  REAL :: R_ss, Z_ss
  REAL :: error_check
  REAL :: dq1,dq2
  REAL :: d0_s1,d0_s2
  REAL :: dp_s1,dp_s2
  REAL :: dffp_s1,dffp_s2
  REAL :: vprime, vpp, dvpp1,dvpp2
  REAL :: ave_M1,ave_M2,ave_M3,ave_M4
  REAL :: p_prime_M,q_prime_M,H
  !
  !-----------------------
  !
  ! compute the first and second derivatives of R,Z on the s-grid
  ! and the local radius of curvature.
  ! Note that for the Mercier-Luc coordinate dR/ds = cos(u), dZ/ds = -sin(u)
  ! so (dR/ds)**2+(dZ/ds)**2 = 1, error_check compute the error in this relation
  ! to make sure that the input flux surface coordinates R(s), Z(s) are ok.
  !
  delta_s = 12.0*ds
  ds2 = 12.0*ds**2
  ! write(*,*)"ms = ",ms," ds = ",ds
  ! do m=0,ms
  ! write(*,*)m,"R = ",R(m)," Z = ",Z(m)
  ! enddo
  ! note that the point 0 and ms are the same so m+1->1 and m-1->ms-1 at m=0
  error_check=0.0
  do m=0,ms
     m1 = MOD(ms+m-2,ms)
     m2 = MOD(ms+m-1,ms)
     m3 = MOD(m+1,ms)
     m4 = MOD(m+2,ms)
     R_s = (R(m1)-8.0*R(m2)+8.0*R(m3)-R(m4))/delta_s
     Z_s = (Z(m1)-8.0*Z(m2)+8.0*Z(m3)-Z(m4))/delta_s
     s_p(m) = sqrt(R_s**2 + Z_s**2)
     R_ss = (-R(m1)+16.0*R(m2)-30.0*R(m)+16.0*R(m3)-R(m4))/ds2
     Z_ss = (-Z(m1)+16.0*Z(m2)-30.0*Z(m)+16.0*Z(m3)-Z(m4))/ds2
     r_curv(m) = (s_p(m)**3)/(R_s*Z_ss - Z_s*R_ss)
     sin_u(m) = -Z_s/s_p(m)
     ! write(*,*)m,m1,m2,m3,m4
     ! write(*,*)m," r_curv =",r_curv(m),"sin_u =",sin_u(m)
     ! write(*,*)"error=",m,ABS(s_p(m) -1.0)
     ! error_check = error_check + (s_p(m)**2 -1.0)**2
  enddo
  ! open(3,file='curv.dat',status='replace')
  ! do m=0,ms
  ! write(3,*)r_curv(m),sin_u(m),s_p(m)
  ! enddo
  ! close(3)
  !
  error_check = SQRT(error_check/real(ms))
  if(error_check .gt. 1.00) &
       write(*,*)"error in s-grid derivative = ",error_check
  !---------------------------------------------------------------
  ! Compute f=R*Bt such that the eikonal S which solves
  ! B*Grad(S)=0 has the correct quasi-periodicity S(s+Ls)=S(s)-2*pi*q_s
  !
  ! Ls
  ! / ds
  ! 1/f = | ---------------
  ! / R**2 Bp 2pi q
  ! 0
  !
  ! f -> R for a circular flux-surface
  !
  ! Define psi_x
  do m=0,ms
     psi_x(m)=R(m)*Bp(m)
  enddo
  !
  ! First compute 2 pi q/f:
  !
  f = 0.0
  do m=1,ms
     f = f &
          +0.5*ds*(s_p(m-1)/(R(m-1)*psi_x(m-1)) + s_p(m)/(R(m)*psi_x(m)))
     ! fixed error that changes results slightly for TGLF_1.93 compared to previous.
     ! +ds*s_p(m)*2.0/(R(m-1)*psi_x(m-1)+R(m)*psi_x(m))
  enddo
  !
  f = pi_2*q_s/f
  ! write(*,*)"f = ",f,q_s
  ! write(*,*)"ds=",ds
  !
  !-----------------------------------------------------------
  !-----------------------------------------------------------
  ! Compute toroidal and total fields:
  !
  ! 2 2 2
  ! B = B + B
  ! p t
  !
  ! B = f/R
  ! t
  !
  do m=0,ms
     Bt(m) = f/R(m)
     B(m) = SQRT(Bt(m)**2 + Bp(m)**2)
     ! write(*,*)m,"Bt =",Bt(m)," B =",B(m)
  enddo
  !----------------------------------------------------------
  !---------------------------------------------------------------
  ! Compute Miller's D , D and D needed for kx.
  ! 0 p ff'
  !
  d_0(0) = 0.0
  d_p(0) = 0.0
  d_ffp(0) = 0.0
  !
  dq1 = ds*s_p(0)*f/(R(0)*psi_x(0)**2)
  d0_s1 = -dq1*(2.0/r_curv(0)+2.0*sin_u(0)/R(0))
  dp_s1 = dq1*4.0*pi*R(0)/Bp(0)
  dffp_s1 = dq1*(R(0)/Bp(0))*(B(0)/f)**2
  !
  do m=1,ms
     dq2 = ds*s_p(m)*f/(R(m)*psi_x(m)**2)
     d0_s2 = -dq2*(2.0/r_curv(m)+2.0*sin_u(m)/R(m))
     dp_s2 = dq2*4.0*pi*R(m)/Bp(m)
     dffp_s2 = dq2*(R(m)/Bp(m))*(B(m)/f)**2
     !
     d_0(m) = d_0(m-1)+0.5*(d0_s1+d0_s2)
     d_p(m) = d_p(m-1)+0.5*(dp_s1+dp_s2)
     d_ffp(m) = d_ffp(m-1)+0.5*(dffp_s1+dffp_s2)
     !
     d0_s1 = d0_s2
     dp_s1 = dp_s2
     dffp_s1 = dffp_s2
     !
  enddo
  !
  ! d_p(0)=0.0
  ! do m=1,ms
  ! d_p(m)=d_p(m-1)+d_0(ms-(m-1))-d_0(ms-m)
  ! write(*,*)t_s(m),d_0(m),d_p(m),d_ffp(m)
  ! write(*,*)m,d_p(m),d_0(ms-m)-d_0(ms),d_0(m)
  ! enddo
  !
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Begin computing geometric quantities required for solution
  ! of gyrokinetic equation:
  !
  ! - b_geo replaces bmaj(j)=b_theta(j)/b_unit
  !
  ! - pk_geo is close to pk=2*rmin/(rmaj*q), the coefficient
  ! of d/dtheta
  !
  ! - qrat_geo -> 1 in a circle
  !
  ! Note that for the physical quantity
  !
  ! k_theta = nq/r
  !
  ! we use
  !
  ! kyrhos_s = n*q_s/rmin_s*rhos_unit_s
  !
  ! which is exactly the same as for the circle.
  !
  ! Also, "omega_star" remains unchanged from circle with
  ! logarithmic density gradients along minor axis.
  !
  ! - "ky*rhos" in Bessel function is kyrhos_s*qrat_geo(j)/b_geo(j)
  !
  ! - "kx*rhos" is kxoky_geo(j)*ky*rhos
  !
  do m=0,ms
     b_geo(m) = B(m)
     pk_geo(m) = 2.0*Bp(m)/B(m)
     qrat_geo(m) = (rmin_s/R(m))*(B(m)/Bp(m))/q_s
  enddo
  !---------------------------------------------------------------
  !
  ! Determine ff_prime from:
  !
  ! 2 pi q_prime = d_0(ms)
  ! +d_p(ms)*p_prime
  ! +d_ffp(ms)*ff_prime
  !
  ff_prime = (pi_2*q_prime_s-d_0(ms)-d_p(ms)*p_prime_s) &
       /d_ffp(ms)
  ! write(*,*)"ff_prime=",ff_prime,"f=",f
  ! write(*,*)"d_0,d_p,d_ffp",d_0(ms),d_p(ms),d_ffp(ms)
  !---------------------------------------------------------------

  !--------------------------------------------------------------
  ! Compute [[kx/ky]] (herein, kxoky_geo) from Waltz-Miller [2]
  ! paper.
  ! 2
  ! (R B_p) S1
  ! kxoky_geo = ---------- ------
  ! B R B_p
  !
  ! S1
  ! ------ = -(d_0(theta)+d_p(theta)*p_prime+d_ffp(theta)*ff_prime)
  ! R B_p
  !
  do m=0,ms

     S_prime(m) = -(d_0(m)+d_p(m)*p_prime_s+d_ffp(m)*ff_prime)
     kx_factor(m) = (psi_x(m)**2)/B(m)
     kxoky_geo(m) = S_prime(m)*kx_factor(m)
     ! write(*,*)"check s_prime",S_prime(m)

  enddo
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  ! Compute drift coefficients:
  !
  !
  ! p_prime_zero forces grad-B-curvature to zero to compensates
  ! for b_par =0
  !
  p_prime_zero_s = -1.0
  if(use_mhd_rule_in)p_prime_zero_s = 0.0
  !
  ! write(*,*)"debug p_prime_zero",p_prime_zero_s

  do m=0,ms

     epsl_geo(m) = 2.0/rmaj_s*qrat_geo(m)/b_geo(m)

     ! Waltz/Miller [[cos_p]] 

     costheta_p_geo(m) = p_prime_zero_s*rmaj_s*(Bp(m)/B(m)**2)*(4.0*pi*R(m)*p_prime_s)

     ! Waltz/Miller [[cos]] 

     costheta_geo(m) = &
          -rmaj_s*(Bp(m)/B(m)**2)*(Bp(m)/r_curv(m) &
          -(f**2/(Bp(m)*R(m)**3))*sin_u(m)) 

  enddo

  !----------------------------------------------------------
  ! Functions which require theta-derivatives:
  !
  do m=0,ms

     ! Waltz/Miller [[sin]]
     m1=MOD(ms+m-2,ms)
     m2=MOD(ms+m-1,ms)
     m3=MOD(m+1,ms)
     m4=MOD(m+2,ms)
     sintheta_geo(m) = -rmaj_s*(f/(R(m)*B(m)**2))* &
          (B(m1)-8.0*B(m2)+8.0*B(m3)-B(m4))/(delta_s*s_p(m))
     ! write(*,*)m,m1,m2,m3,m4,"sintheta_geo=",sintheta_geo(m)

  enddo
  !
  ! compute vprime and vpp
  !
  vprime = 0.0
  vpp = 0.0
  dvpp1 = (s_p(0)*R(0)/psi_x(0)**3) &
       *(4.0*pi*p_prime_s*R(0)**2 + ff_prime - 2.0*psi_x(0)/r_curv(0))
  do m=1,ms
     dvpp2 = (s_p(m)*R(m)/psi_x(m)**3) &
          *(4.0*pi*p_prime_s*R(m)**2 + ff_prime - 2.0*psi_x(m)/r_curv(m))
     vprime = vprime + 0.5*ds*(s_p(m-1)/Bp(m-1) + s_p(m)/Bp(m))
     vpp = vpp + 0.5*ds*(dvpp1 + dvpp2)
     dvpp1 = dvpp2
  enddo
  vprime = pi_2*vprime
  vpp = pi_2*vpp
  ! write(*,*)"vprime = ",vprime,"vpp = ",vpp
  !
  ave_M1 = 0.0
  ave_M2 = 0.0
  ave_M3 = 0.0
  ave_M4 = 0.0
  do m=1,ms
     ave_M1 = ave_M1 + 0.5*ds*(s_p(m-1)/Bp(m-1)**3 + s_p(m)/Bp(m)**3)
     ave_M2 = ave_M2 + 0.5*ds*(s_p(m-1)*R(m-1)/psi_x(m-1)**3 + s_p(m)*R(m)/psi_x(m)**3)
     ave_M3 = ave_M3 + 0.5*ds*(s_p(m-1)*(B(m-1)/psi_x(m-1))**2/Bp(m-1) + s_p(m)*(B(m)/psi_x(m))**2/Bp(m))
     ave_M4 = ave_M4 + 0.5*ds*(s_p(m-1)*B(m-1)**2/Bp(m-1) + s_p(m)*B(m)**2/Bp(m))
  enddo
  ! write(*,*)"M1=",ave_M1,"M2=",ave_M2,"M3=",ave_M3,"M4=",ave_M4
  p_prime_M = 4.0*pi*p_prime_s
  q_prime_M = pi_2*q_prime_s
  DM_out = 0.25 + (p_prime_M/MAX(q_prime_M**2,1e-12))*((vpp/pi_2 - p_prime_M*ave_M1)*ave_M3 &
       + (f**2*p_prime_M*ave_M2 - q_prime_M*f)*ave_M2)
  H = (f*p_prime_M*q_prime_M/MAX(q_prime_M**2,1e-12))*ave_M3*(ave_M2/ave_M3 - vprime/(pi_2*ave_M4))
  ! write(*,*)"H = ",H
  DR_out = DM_out - (0.5 - H)**2
  !
  !
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !
  ! call mercier_write
  ! open(2,file='functions_geo.dat',status='replace')
  ! do m=ms/2,0,-1
  ! write(2,*)t_s(m),costheta_geo(m),sintheta_geo(m),kxoky_geo(m)
  ! enddo
  ! do m=ms-1,ms/2,-1
  ! write(2,*)t_s(m)+pi_2,costheta_geo(m),sintheta_geo(m) &
  ! ,kx_factor(m)*(S_prime(m)-S_prime(ms))
  ! enddo
  close(2)
  !
END SUBROUTINE mercier_luc
!
!
!---------------------------------------------------------------
! mercier_write.f
!
! PURPOSE:
! Output of some internal MILLER variables, as well
! as selected values of required external variables
! like costheta, etc.
!
! REVISIONS:
! 30 june 2005: gms
! Based on MILLER_write.f90 from GYRO
!---------------------------------------------------------------

SUBROUTINE mercier_write
  !
  USE tglf_global
  USE tglf_sgrid
  !
  IMPLICIT NONE
  INTEGER :: j,jj
  character(len=7) ang(0:4)
  character(len=11) names(0:12)
  REAL :: theta_write(0:4)
  !---------------------------------------------------

  open(unit=3,file='mercier.out',status='replace')
  write(3,10) 'rmin',rmin_loc
  write(3,10) 'rmaj0',rmaj_loc
  write(3,10) 'kappa',kappa_loc
  write(3,10) 's_kappa',s_kappa_loc
  write(3,10) 'delta',delta_loc
  write(3,10) 's_delta',s_delta_loc
  write(3,10) 'drmindx',drmindx_loc
  write(3,10) 'drmajdx',drmajdx_loc
  write(3,10) 'q',q_loc
  write(3,10) 'f',f
  write(3,10) 'ff_prime',ff_prime
  write(3,10) 'q_prime',q_prime_s
  write(3,10) 'p_prime',p_prime_s

  ang(0) = '0.00*pi'
  ang(1) = '0.50*pi'
  ang(2) = '1.00*pi'
  ang(3) = '1.50*pi'
  ang(4) = '2.00*pi'
  names(0) = 'theta '
  names(1) = 't_s '
  names(2) = 'b_geo '
  names(3) = 'grad_r '
  names(4) = 'R '
  names(5) = 'Z '
  names(6) = 'pk_geo '
  names(7) = 'qrat_geo '
  names(8) = 'cos '
  names(9) = 'cos_p '
  names(10) = 'sin '
  names(11) = 'epsl '
  names(12) = 'kxoky_geo '
  !

  theta_write(0) = 0.0
  theta_write(1) = 0.5*pi
  theta_write(2) = 1.0*pi
  theta_write(3) = 1.5*pi
  theta_write(4) = 2.0*pi
  write(3,*)
  write(3,30)names(0),names(1),(names(j),j=2,7)
  write(3,*)
  !
  do jj=0,4
     do j=ms,2,-1
        if(t_s(j).gt.theta_write(jj)-pi_2)exit
     enddo
     if(ABS(t_s(j-1)+pi_2-theta_write(jj)) &
          .lt.ABS(t_s(j)+pi_2-theta_write(jj)))j=j-1
     !
     write(3,20)ang(jj),t_s(j),b_geo(j),q_s*R(j)*Bp(j)/rmin_loc, &
          R(j),Z(j),pk_geo(j),qrat_geo(j)
     !
  enddo
  !
  write(3,*)
  write(3,30)names(0),names(1),(names(j),j=8,12)
  write(3,*)
  !
  do jj=0,4
     do j=ms,2,-1
        if(t_s(j).gt.theta_write(jj)-pi_2)exit
     enddo
     if(ABS(t_s(j-1)+pi_2-theta_write(jj)) &
          .lt.ABS(t_s(j)+pi_2-theta_write(jj)))j=j-1
     !
     write(3,20) ang(jj),t_s(j),costheta_geo(j),costheta_p_geo(j), &
          sintheta_geo(j),epsl_geo(j),kxoky_geo(j), &
          (S_prime(j)-S_prime(ms))*kx_factor(j)
     !
  enddo
  close(3)

10 format(t2,a,':',t20,f10.6)
20 format(t2,a,1x,10(es10.3,1x))
30 format(t2,10a)
  !
END SUBROUTINE mercier_write
!
!
!---------------------------------------------------------------
! miller_geo.f [various callers]
!
! PURPOSE:
! Core routine for calculation of MILLER shaped flux surface
! quantities R(theta), Z(theta), grad_r(theta)
!
! DEFINITIONS:
! Length is in units of rmin at the last closed
! flux surface (a_unit).
!
! R(r,theta) = R0(r) + r cos[theta + x_delta*sin(theta)]
! Z(r,theta) = kappa r sin(theta)
!
!
! NOTES:
! To follow what's going on in this routine, it is
! necessary to have the following papers at hand:
!
! (1) R.L. Miller, M.S. Chu, J.M. Greene, Y.R. Lin-liu
! and R.E. Waltz, 'Noncircular, finite aspect ratio,
! local equilibrium model', Phys. Plasmas 5 (1998) 973.
!
! (2) R.E. Waltz and R.L. Miller, 'Ion temperature gradient
! turbulence simulations and plasma flux surface shape',
! Phys. Plasmas 6 (1999) 4265.
!
! **See MILLER_inputs for control variables**
!
! REVISIONS:
! 13 Aug 02: jc
! Documentation after MILLER de-spaghetti.
! 18 Nov 03: jc
! Removed reset of s_delta and s_kappa when they
! are input as zero.
! 24 Nov 03: jc
! Entire code has been rewritten for greater
! efficiency and clarity. Documentation greatly
! improved as well.
! 24 June 05: gms
! produced this version which only computes R,Z,Bp for input
! into mercier_luc.f which completes the calculation of the Waltz-Miller
! functions [[sin]],[[cos]], etc.
! 15 June 2010: updated to GYRO conventions with squarness (zeta_loc) and
! squarness shear s_zeta_loc = rmin*d(zeta)/dx. Included elevation Zmaj_loc.
! Also changed definition of s_delta = rmin*d(delta)/dx from Waltz-Miller convention to GYRO's.
!---------------------------------------------------------------

SUBROUTINE miller_geo
  !
  USE tglf_global
  USE tglf_sgrid
  !
  IMPLICIT NONE
  !
  !-------------------------------------------
  !
  INTEGER,PARAMETER :: mts=5
  INTEGER :: m
  !-----------------------------------------------
  !
  REAL :: theta, x_delta
  REAL :: dtheta
  REAL :: arg_r,darg_r,arg_z,darg_z
  REAL :: R_t,Z_t
  REAL :: R_r, Z_r
  REAL :: l_t, grad_r, det
  REAL :: scale_max, l_t1, arclength
  !
  !-------------------------------------------
  ! set global input values
  !
  rmin_input = rmin_loc
  Rmaj_input = rmaj_loc

  ! write(*,*)"miller_geo"
  x_delta = ASIN(delta_loc)
  ! write(*,*)"pi = ",pi," x_delta = ",x_delta
  !
  ! set the flux surface constants needed for mercier-luc
  !
  Rmaj_s = rmaj_loc
  Zmaj_s = zmaj_loc
  q_s = q_loc
  if(rmin_loc.lt.0.00001)rmin_loc=0.00001
  rmin_s = rmin_loc
  p_prime_s = p_prime_loc
  q_prime_s = q_prime_loc
  !
  !--------------------------------------------------------------
  !
  ! compute the arclength around the flux surface
  !
  theta = 0.0
  arg_r = theta+x_delta*sin(theta)
  darg_r = 1.0+x_delta*cos(theta)
  arg_z = theta + zeta_loc*sin(2.0*theta)
  darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)
  r_t = -rmin_loc*sin(arg_r)*darg_r
  z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z
  l_t = SQRT(r_t**2+z_t**2)
  ! scale dtheta by l_t to keep mts points in each ds interval of size pi_2/ms
  dtheta = pi_2/(REAL(mts*ms)*l_t)
  !
  l_t1 = l_t
  scale_max=l_t
  arclength = 0.0
  do while(theta.lt.pi_2)
     theta = theta + dtheta
     if(theta.gt.pi_2)then
        theta=theta-dtheta
        dtheta=pi_2-theta
        theta = pi_2
     endif
     ! write(*,*)"theta = ",theta,"dtheta=",dtheta
     arg_r = theta+x_delta*sin(theta)
     ! d(arg_r)/dtheta
     darg_r = 1.0+x_delta*cos(theta)
     ! dR/dtheta
     r_t = -rmin_loc*sin(arg_r)*darg_r
     !
     arg_z = theta + zeta_loc*sin(2.0*theta)
     ! d(arg_z)/dtheta
     darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)
     ! dZ/dtheta
     z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z
     ! dl/dtheta
     l_t = SQRT(r_t**2+z_t**2)
     ! arclength along flux surface in poloidal direction
     arclength = arclength + 0.50*(l_t + l_t1)*dtheta
     ! save maximum expansion scale for later
     if(l_t.gt.scale_max) scale_max = l_t
     l_t1 = l_t
  enddo
  Ls = arclength
  !
  ! debug
  ! write(*,*)"arclength = ", arclength
  ! write(*,*)"scale_max =",scale_max
  !
  ! Find the theta points which map to an equally spaced s-grid of ms points along the arclength
  ! going clockwise from the outboard midplane around the flux surface
  ! by searching for the theta where dR**2 + dZ**2 >= ds**2 for a centered difference df=f(m+1)-f(m-1).
  ! This keeps the finite difference error of dR/ds, dZ/ds on the s-grid small
  !
  ds = arclength/REAL(ms)
  ! write(*,*)"ds=",ds
  t_s(0)=0.0
  t_s(ms)=-pi_2
  ! make a first guess based on theta=0.0
  theta=0.0
  arg_r = theta+x_delta*sin(theta)
  darg_r = 1.0+x_delta*cos(theta)
  arg_z = theta + zeta_loc*sin(2.0*theta)
  darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)
  r_t = -rmin_loc*sin(arg_r)*darg_r
  z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z
  l_t = SQRT(r_t**2+z_t**2)
  dtheta = -ds/l_t
  theta=dtheta
  l_t1=l_t
  !
  do m=1,ms/2
     arg_r = theta+x_delta*sin(theta)
     darg_r = 1.0+x_delta*cos(theta)
     arg_z = theta + zeta_loc*sin(2.0*theta)
     darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)
     r_t = -rmin_loc*sin(arg_r)*darg_r
     z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z
     l_t = SQRT(r_t**2+z_t**2)
     dtheta = -ds/(0.5*(l_t+l_t1))
     t_s(m)=t_s(m-1)+dtheta
     theta = t_s(m) +dtheta
     l_t1=l_t
  enddo
  ! distribute endpoint error over interior points
  dtheta = (t_s(ms/2)-(-pi))/REAL(ms/2)
  ! write(*,*)"enpoint error =",dtheta
  ! dtheta=0.0
  ! t_s(ms/2)=-pi
  do m=1,ms/2
     t_s(m) = t_s(m)-REAL(m)*dtheta
     t_s(ms-m)=-pi_2 - t_s(m)
  enddo
  ! write(*,*)"t_s(ms/2)+pi=",t_s(ms/2)+pi
  !
  ! open(2,file='t_s.dat',status='replace')
  ! do m=0,ms
  ! write(2,*)m,t_s(m)
  ! enddo
  ! close(2)
  !
  !--------------------------------------------------------------
  !
  !
  !---------------------------------------------------------------
  ! all equilibrium functions satisfy
  !
  ! f(0) = f(l_theta)
  !
  ! Loop to compute most geometrical quantities needed for Mercie-Luc expansion
  ! R, Z, R*Bp on flux surface s-grid
  !
  ! NOTES:
  ! If grad_r_theta diverges because denominator goes
  ! through zero, magnetic field lines are intersecting
  ! and the magnetic surfaces are not nested.
  !
  do m=0,ms

     theta = t_s(m)
     arg_r = theta + x_delta*sin(theta)
     darg_r = 1.0 + x_delta*cos(theta)
     arg_z = theta + zeta_loc*sin(2.0*theta)
     darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)

     ! R(theta)
     ! Z(theta)

     R(m) = rmaj_loc + rmin_loc*cos(arg_r)
     Z(m) = Zmaj_loc + kappa_loc*rmin_loc*sin(arg_z)

     ! dR/dtheta
     ! dZ/dtheta

     R_t = -rmin_loc*sin(arg_r)*darg_r
     Z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z

     ! dl/dtheta

     l_t = SQRT(R_t**2+Z_t**2)
     ! dR/dr
     ! dZ/dr
     R_r = drmajdx_loc + drmindx_loc*cos(arg_r) &
          -sin(arg_r)*s_delta_loc*sin(theta)/sqrt(1.0 - delta_loc**2)
     Z_r = dzmajdx_loc + kappa_loc*sin(arg_z)*(drmindx_loc +s_kappa_loc) &
          +kappa_loc*cos(arg_z)*s_zeta_loc*sin(2.0*theta)
     ! Jacobian
     det = R_r*z_t - R_t*Z_r
     ! grad_r
     grad_r = ABS(l_t/det)

     if(m.eq.0)then
        B_unit = 1.0/grad_r ! B_unit choosen to make bx(0)=ky**2 i.e. qrat_geo(0)/b_geo(0)=1.0
        if(drmindx_loc.eq.1.0)B_unit=1.0 ! Waltz-Miller convention
        ! write(*,*)"B_unit = ",B_unit
     endif
     ! Bp = (r/q) B_unit*Abs(grad_r)/R


     Bp(m) = (rmin_s/(q_s*R(m)))*grad_r*B_unit
     !
  enddo
  !
  !
  ! write(*,*)"rmin=",rmin_loc,"rmaj=",rmaj_loc,"q_loc=",q_loc
  ! write(*,*)shift_loc,kappa_loc,delta_loc,s_kappa_loc,s_delta_loc
  ! open(3,file='Bp.dat',status='replace')
  ! open(2,file='RZ.dat',status='replace')
  ! do m=0,ms
  ! write(3,*)m,Bp(m)
  ! write(2,*)m,R(m),Z(m)
  ! enddo
  ! close(3)
  ! close(2)
  !
  ! Note the definitions:
  !
  ! q_prime -> dq/dpsi = dq/dr dr/dpsi
  ! p_prime -> dp/dpsi = dp/dr dr/dpsi
  !
  ! R B_p
  ! and dpsi/dr = -------- = b_unit (r/q)
  ! |grad r|
  !
  ! So, we can write:
  ! 2
  ! q_prime = (q/r) s/b_unit
  ! p_prime = (q/r) (1/b_unit) dp/dr
  ! = (q/r) (1/b_unit) p * dlnpdr
  ! rescale p_prime and q_prime to keep normalized eikonal invariant
  p_prime_s = p_prime_s*B_unit ! B_unit**2/B_unit
  q_prime_s = q_prime_s/B_unit !
  ! write(*,*)"p_prime_s =",p_prime_s,"q_prime_s =",q_prime_s
  !
END SUBROUTINE miller_geo
!
!---------------------------------------------------------------
!
SUBROUTINE fourier_geo
  ! maps the fourier representation of R,Z onto the contour s-grid
  ! and computes Bp,t_s, Ls used by mercier-luc
  !
  USE tglf_global
  USE tglf_sgrid
  !
  IMPLICIT NONE
  !
  !-------------------------------------------
  !
  INTEGER,PARAMETER :: mts=5
  INTEGER :: n, m
  !-----------------------------------------------
  !
  REAL :: theta, nr
  REAL :: dtheta
  REAL :: arg
  REAL :: R_t,Z_t
  REAL :: R_r, Z_r
  REAL :: l_t, grad_r, det
  REAL :: scale_max, l_t1, arclength
  !
  !--------------------------------------------------------------
  !
  ! set the flux surface constants needed for mercier-luc
  !
  Rmaj_s = fourier_in(1,0)/2.0
  Zmaj_s = fourier_in(3,0)/2.0
  q_s = q_fourier_in
  rmin_s = fourier_in(1,1)
  ! write(*,*)"Rmaj_s=",Rmaj_s,"Zmaj_s=",Zmaj_s,"q_s=",q_s,"rmin_s=",rmin_s
  rmin_s=MAX(rmin_s,0.00001)
  p_prime_s = p_prime_fourier_in
  q_prime_s = q_prime_fourier_in
  !
  ! set global input values
  !
  rmin_input = rmin_s
  Rmaj_input = Rmaj_s
  !
  ! compute the arclength around the flux surface
  !
  theta = 0.0
  R_t = 0.0
  Z_t = 0.0
  ! write(*,*)"nfourier_in=",nfourier_in
  do n=1,nfourier_in
     nr = REAL(n)
     arg = nr*theta
     R_t = R_t -nr*fourier_in(1,n)*SIN(arg) + nr*fourier_in(2,n)*COS(arg)
     Z_t = Z_t -nr*fourier_in(3,n)*SIN(arg) + nr*fourier_in(4,n)*COS(arg)
  enddo
  l_t = SQRT(r_t**2+z_t**2)
  ! scale dtheta by l_t to keep mts points in each ds interval of size pi_2/ms
  dtheta = pi_2/(REAL(mts*ms)*l_t)
  !
  l_t1 = l_t
  scale_max=l_t
  arclength = 0.0
  do while(theta.lt.pi_2)
     theta = theta + dtheta
     if(theta.gt.pi_2)then
        theta=theta-dtheta
        dtheta=pi_2-theta
        theta = pi_2
     endif
     R_t = 0.0
     Z_t = 0.0
     do n=1,nfourier_in
        nr = REAL(n)
        arg = nr*theta
        R_t = R_t -nr*fourier_in(1,n)*SIN(arg) + nr*fourier_in(2,n)*COS(arg)
        Z_t = Z_t -nr*fourier_in(3,n)*SIN(arg) + nr*fourier_in(4,n)*COS(arg)
     enddo
     ! dl/dtheta
     l_t = SQRT(r_t**2+z_t**2)
     ! arclength along flux surface in poloidal direction
     arclength = arclength + 0.50*(l_t + l_t1)*dtheta
     ! save maximum expansion scale for later
     if(l_t.gt.scale_max) scale_max = l_t
     l_t1 = l_t
  enddo
  Ls = arclength
  !
  ! debug
  ! write(*,*)"arclength = ", arclength
  ! write(*,*)"scale_max =",scale_max
  !
  ! Find the theta points which map to an equally spaced s-grid of ms points along the arclength
  ! going clockwise from the outboard midplane around the flux surface
  ! by searching for the theta where dR**2 + dZ**2 >= ds**2 for a centered difference df=f(m+1)-f(m-1).
  ! This keeps the finite difference error of dR/ds, dZ/ds on the s-grid small
  !
  ds = arclength/REAL(ms)
  ! write(*,*)"ds=",ds
  t_s(0)=0.0
  t_s(ms)=-pi_2
  ! make a first guess based on theta=0.0
  theta=0.0
  R_t = 0.0
  Z_t = 0.0
  do n=1,nfourier_in
     nr = REAL(n)
     arg = nr*theta
     R_t = R_t -nr*fourier_in(1,n)*SIN(arg) + nr*fourier_in(2,n)*COS(arg)
     Z_t = Z_t -nr*fourier_in(3,n)*SIN(arg) + nr*fourier_in(4,n)*COS(arg)
  enddo
  l_t = SQRT(r_t**2+z_t**2)
  dtheta = -ds/l_t
  theta=dtheta
  l_t1=l_t
  !
  do m=1,ms/2
     R_t = 0.0
     Z_t = 0.0
     do n=1,nfourier_in
        nr = REAL(n)
        arg = nr*theta
        R_t = R_t -nr*fourier_in(1,n)*SIN(arg) + nr*fourier_in(2,n)*COS(arg)
        Z_t = Z_t -nr*fourier_in(3,n)*SIN(arg) + nr*fourier_in(4,n)*COS(arg)
     enddo
     l_t = SQRT(r_t**2+z_t**2)
     dtheta = -ds/(0.5*(l_t+l_t1))
     t_s(m)=t_s(m-1)+dtheta
     theta = t_s(m) +dtheta
     l_t1=l_t
  enddo
  ! distribute endpoint error over interior points
  dtheta = (t_s(ms/2)-(-pi))/REAL(ms/2)
  ! write(*,*)"enpoint error =",dtheta
  ! dtheta=0.0
  ! t_s(ms/2)=-pi
  do m=1,ms/2
     t_s(m) = t_s(m)-REAL(m)*dtheta
     t_s(ms-m)=-pi_2 - t_s(m)
  enddo
  ! write(*,*)"t_s(ms/2)+pi=",t_s(ms/2)+pi
  !
  ! open(2,file='t_s.dat',status='replace')
  ! do m=0,ms
  ! write(2,*)m,t_s(m)
  ! enddo
  ! close(2)
  !
  !---------------------------------------------------------------
  ! all equilibrium functions satisfy
  !
  ! f(0) = f(l_theta)
  !
  ! Loop to compute most geometrical quantities needed for Mercie-Luc expansion
  ! R, Z, R*Bp on flux surface s-grid
  !
  ! NOTES:
  ! If grad_r_theta diverges because denominator goes
  ! through zero, magnetic field lines are intersecting
  ! and the magnetic surfaces are not nested.
  !
  !--------------------------------------------------------------
  !
  ! compute R,Z,BP on the arc-length s-grid
  !
  do m=0,ms
     theta = t_s(m)
     R(m) = fourier_in(1,0)/2.0
     Z(m) = fourier_in(3,0)/2.0
     R_t = 0.0
     Z_t = 0.0
     R_r = fourier_in(5,0)/2.0
     Z_r = fourier_in(7,0)/2.0
     do n=1,nfourier_in
        ! R(theta)
        ! Z(theta)
        nr=REAL(n)
        arg = theta*nr
        R(m) = R(m) + fourier_in(1,n)*COS(arg) + fourier_in(2,n)*SIN(arg)
        Z(m) = Z(m) + fourier_in(3,n)*COS(arg) + fourier_in(4,n)*SIN(arg)

        ! dR/dtheta
        ! dZ/dtheta
        ! note these are derivatives wrt the fourier theta variable
        ! the factors dtheta/ds needed to take s-derivatives drop out of the grad_r calculation

        R_t = R_t -nr*fourier_in(1,n)*SIN(arg) + nr*fourier_in(2,n)*COS(arg)
        Z_t = Z_t -nr*fourier_in(3,n)*SIN(arg) + nr*fourier_in(4,n)*COS(arg)

        ! dR/dr
        ! dZ/dr
        R_r = R_r + fourier_in(5,n)*COS(arg) + fourier_in(6,n)*SIN(arg)
        Z_r = Z_r + fourier_in(7,n)*COS(arg) + fourier_in(8,n)*SIN(arg)
     enddo
     ! dl/dtheta

     l_t = SQRT(R_t**2+Z_t**2)
     ! Jacobian
     det = R_r*z_t - R_t*Z_r
     ! grad_r
     grad_r = ABS(l_t/det)
     !
     B_unit = 1.0
     ! if(m.eq.0)then
     ! B_unit = 1.0/grad_r ! B_unit choosen to make bx(0)=ky**2 i.e. qrat_geo(0)/b_geo(0)=1.0
     ! if(drmindx_loc.eq.1.0)B_unit=1.0 ! Waltz-Miller convention
     ! write(*,*)"B_unit = ",B_unit
     ! endif
     ! Bp = (r/q) B_unit*Abs(grad_r)/R


     Bp(m) = (rmin_s/(q_s*R(m)))*grad_r*B_unit
     !
  enddo

  !
  ! open(3,file='Bp.dat',status='replace')
  ! open(2,file='RZ.dat',status='replace')
  ! do m=0,ms
  ! write(3,*)m,Bp(m)
  ! write(2,*)m,R(m),Z(m)
  ! enddo
  ! close(3)
  ! close(2)
  !
  ! Note the definitions:
  !
  ! q_prime -> dq/dpsi = dq/dr dr/dpsi
  ! p_prime -> dp/dpsi = dp/dr dr/dpsi
  !
  ! R B_p
  ! and dpsi/dr = -------- = b_unit (r/q)
  ! |grad r|
  !
  ! So, we can write:
  ! 2
  ! q_prime = (q/r) s/b_unit
  ! p_prime = (q/r) (1/b_unit) dp/dr
  ! = (q/r) (1/b_unit) p * dlnpdr
  ! rescale p_prime and q_prime to keep normalized eikonal invariant
  p_prime_s = p_prime_s*B_unit ! B_unit**2/B_unit
  q_prime_s = q_prime_s/B_unit !
  ! write(*,*)"p_prime_s =",p_prime_s,"q_prime_s =",q_prime_s
  !
END SUBROUTINE fourier_geo
!
!--------------------------------------------------------------------
!
SUBROUTINE ELITE_geo
  !
  ! interpolates the input R_ELITE,Z_ELITE,Bp_ELITE onto the s-grid
  ! used by mercier_luc and sets the arclength Ls and differential ds
  !
  USE tglf_global
  USE tglf_sgrid
  !
  IMPLICIT NONE
  !
  INTEGER :: i,j,k,imax,im,ip
  REAL :: arclength,drde,dzde,de
  REAL :: Rmax,Zmax,Bpmax,emax,Zmin,Rmin
  REAL :: e_length,s_length
  REAL :: da
  REAL :: a,b,c
  REAL :: area,R0,Z0
  !
  ! compute the arclength, area and centroid cooridinates of the flux surface
  !
  Rmax = R_ELITE(0)
  Rmin = R_ELITE(0)
  Zmax = Z_ELITE(0)
  Zmin = Z_ELITE(0)
  imax=0
  drde = (R_ELITE(1)-R_ELITE(n_ELITE))/2.0
  dzde = (Z_ELITE(1)-Z_ELITE(n_ELITE))/2.0
  da = SQRT(drde**2 + dzde**2)
  arclength = da
  area = ABS(drde)*ABS(Z_ELITE(0))
  ! integral from 0.0 to ABS(Z_ELITE(i)) at each R_ELITE assuming Z=0 is inside of the flux surface
  ! this makes it unecessary to find the Z reflected across Z=0 a the same R
  Z0 = ABS(drde)*ABS(Z_ELITE(0))*Z_ELITE(0)/2.0
  R0 = ABS(drde)*ABS(Z_ELITE(0))*(R_ELITE(1)+R_ELITE(n_ELITE))/2.0
  do i=1,n_ELITE-1
     if(Z_ELITE(i).gt.Zmax)Zmax = Z_ELITE(i)
     if(Z_ELITE(i).lt.Zmin)Zmin = Z_ELITE(i)
     if(R_ELITE(i).lt.Rmin)Rmin = R_ELITE(i)
     if(R_ELITE(i).gt.Rmax)then
        Rmax = R_ELITE(i)
        imax = i
     endif
     drde = (R_ELITE(i+1)-R_ELITE(i-1))/2.0
     dzde = (Z_ELITE(i+1)-Z_ELITE(i-1))/2.0
     da = SQRT(drde**2 + dzde**2)
     arclength = arclength + da
     area = area + ABS(drde)*ABS(Z_ELITE(i))
     Z0 = Z0 + ABS(drde)*ABS(Z_ELITE(i))*Z_ELITE(i)/2.0
     R0 = R0 + ABS(drde)*ABS(Z_ELITE(i))*(R_ELITE(i+1)+R_ELITE(i-1))/2.0
  enddo
  write(*,*)"imax = ",imax,"Rmax = ",Rmax,"Rmin = ",Rmin
  write(*,*)"Zmax = ",Zmax,"Zmin = ",Zmin
  de = arclength/REAL(n_ELITE)
  ds = arclength/REAL(ms)
  Z0 = Z0/area
  R0 = R0/area
  write(*,*)"area = ",area,"volume = ",pi_2*R0*area
  write(*,*)"Z0 = ",Z0,"R0 = ",R0
  write(*,*)"arclength = ",arclength
  !
  ! set global input values
  !
  rmin_input = (Rmax-Rmin)/2.0
  Rmaj_input = (Rmax+Rmin)/2.0

  !
  ! find the major and minor radius at Z0
  !
  ! signdZ=1.0
  ! if(Z_ELITE(0).lt.Z)signZ=-1.0
  !
  ! i1 = 0
  ! do i=0,n_ELITE
  ! if(Z_ELITE
  !
  ! enddo
  !
  ! find maximum of quardratic fit through R(e): R(e) = a + b e + c e^2
  !
  im = imax-1
  if(im.lt.0)im=n_ELITE-1
  ip = imax + 1
  a = R_ELITE(im)
  b = (2.0*R_ELITE(imax)-0.5*R_ELITE(ip)-1.5*R_ELITE(im))/de
  c = (0.5*R_ELITE(ip)+0.5*R_ELITE(im)-R_ELITE(imax))/de**2
  ! dR/de=0 at emax=-b/(2 c)
  emax= -b/(2.0*c)
  write(*,*)"emax/de=",emax/de
  Rmax = a + b*emax + c*emax**2
  ! interpolate Z,Bp onto the point emax
  a = Z_ELITE(im)
  b = (2.0*Z_ELITE(imax)-0.5*Z_ELITE(ip)-1.5*Z_ELITE(im))/de
  c = (0.5*Z_ELITE(ip)+0.5*Z_ELITE(im)-Z_ELITE(imax))/de**2
  Zmax = a + b*emax + c*emax**2
  a = Bp_ELITE(im)
  b = (2.0*Bp_ELITE(imax)-0.5*Bp_ELITE(ip)-1.5*Bp_ELITE(im))/de
  c = (0.5*Bp_ELITE(ip)+0.5*Bp_ELITE(im)-Bp_ELITE(imax))/de**2
  Bpmax = a + b*emax + c*emax**2
  !
  ! compute B_unit = d psi/dr *(q/r)
  !
  rmin_s = arclength/pi_2
  B_unit = Rmax*Bpmax*q_ELITE/rmin_s ! this choice makes bx(0) = ky^2 i.e. qrat_geo(0)/b_geo(0)=1.0
  B_unit = 1.0
  !
  ! interpolate R,Z,Bp onto the s-grid
  !
  if(emax/de.lt.1.0)then
     imax=im
  else
     emax = emax - de
  endif
  e_length = -emax
  s_length = ds
  R(0) = Rmax
  Z(0) = 0.0
  Bp(0) = Bpmax/B_unit
  k = 0
  do i=1,n_ELITE-1
     j = imax + i
     if(j .gt. n_ELITE)then
        ! wrap around the flux surface
        j = imax + i - n_ELITE
     endif
     e_length = e_length + de
     if(e_length.gt.s_length)then
        k = k+1
        R(k) = R_ELITE(j-1)+(R_ELITE(j)-R_ELITE(j-1))*(e_length-s_length)/de
        Z(k) = Z_ELITE(j-1)+(Z_ELITE(j)-Z_ELITE(j-1))*(e_length-s_length)/de - Zmax
        Bp(k) = (Bp_ELITE(j-1)+(Bp_ELITE(j)-Bp_ELITE(j-1))*(e_length-s_length)/de)/B_unit
        s_length = s_length + ds
     endif
  enddo
  R(ms) = R(0)
  Z(ms) = Z(0)
  Bp(ms) = Bp(0)
  write(*,*)"k = ",k
  !
  ! set remaining variables for tglf_sgrid input
  !
  Ls = arclength
  Rmaj_s = Rmax
  p_prime_s = p_prime_ELITE
  q_s = q_ELITE
  q_prime_s = q_prime_ELITE
  ! debug
  write(*,*)"Rmax = ",Rmax
  write(*,*)"Zmax = ",Zmax
  write(*,*)"Bpmax = ",Bpmax
  write(*,*)"B_unit = ",B_unit
  write(*,*)"Ls = ",Ls
  write(*,*)"ds = ",ds
  write(*,*)"Rmaj_s =",Rmaj_s
  write(*,*)"rmin_s = ",rmin_s
  write(*,*)"p_prime_s = ",p_prime_s
  write(*,*)"q_s = ",q_s
  write(*,*)"q_prime_s = ",q_prime_s
  ! open(3,file='Bp.dat',status='replace')
  ! open(2,file='RZ.dat',status='replace')
  ! do i=0,ms
  ! write(3,*)i,Bp(i)
  ! write(2,*)i,R(i),Z(i)
  ! enddo
  ! close(3)
  ! close(2)
  !
END SUBROUTINE ELITE_geo
