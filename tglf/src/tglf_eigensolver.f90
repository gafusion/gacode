      SUBROUTINE tglf_eigensolver
!*********************************************
!
!*********************************************
      USE tglf_dimensions
      USE tglf_global
      USE tglf_closure
      USE tglf_species
      USE tglf_eigen
      USE tglf_coeff
      USE tglf_nuei_coeff
      USE tglf_sgrid
!
      IMPLICIT NONE
!
      CHARACTER(1) :: rightvectors
      INTEGER :: is,js
      INTEGER ::  j1, j2, j, i
      INTEGER :: ifail,ib, jb, ia, ja, ia0, ja0
      INTEGER :: info,lwork
      INTEGER :: xnu_model
      INTEGER :: idum=0
      REAL :: ft2,ft3,ft4,ft5
      REAL :: Linsker,am,bm
      REAL :: w_s, w_dh, w_dg, w_d1, modw_d1,w_d0, w_cd
      REAl :: wd_psi,modwd_psi
      REAL :: E_i,M_i,N_j,J_j
      REAL :: k_par0,k_par1,modk_par0,modk_par1
      REAL :: h10n,h10p1,h10p3,h10r13,h10r33
      REAL :: hnb0,hp1b0,hp3b0,hr11b0,hr13b0,hr33b0
      REAL :: hw113b0,hw133b0,hw333b0
      REAL :: hnbp,hp1bp,hp3bp,hr11bp,hr13bp,hr33bp
      REAL :: hw113bp,hw133bp,hw333bp
      REAL :: c_tor_par_hp1,c_tor_par_hr11
      REAL :: c_tor_par_hr13
      REAL :: g10n,g10p1,g10p3,g10r13,g10r33
      REAL :: gnb0,gp1b0,gp3b0,gr11b0,gr13b0,gr33b0
      REAL :: gw113b0,gw133b0,gw333b0
      REAL :: gnbp,gp1bp,gp3bp,gr11bp,gr13bp,gr33bp
      REAL :: gw113bp,gw133bp,gw333bp
      REAL :: c_tor_par_gp1,c_tor_par_gr11
      REAL :: c_tor_par_gr13
      REAL :: betae_psi,betae_sig
      REAL :: max_freq,test
      REAL :: gradB1
      REAL :: c35
      REAL :: xnuei,d_ab,d_ij,d_ee,d_ab_psi,d_11,d_1
      REAL :: xnuion
      REAL :: xnu_bndry,xnu_hat,xnu_a,xnu_b,xnu_c
      REAL :: xnu_phi_b
      REAL :: xnu_n_b
      REAL :: xnu_p1_b
      REAL :: xnu_p3_b
      REAL :: xnu_u_b
      REAL :: xnu_q1_b
      REAL :: xnu_q3_b
      REAL :: xnu_p1_1
      REAL :: xnu_u_u_1,xnu_u_q3_1
      REAL :: xnu_q1_u_1,xnu_q1_q1_1,xnu_q1_q3_1
      REAL :: xnu_q3_u_1,xnu_q3_q3_1
      REAL :: c01,c02,c03,c04,c05
      REAL :: c06,c07,c08,c09,c010
      REAL :: cb1,cb2,cb3,cb4,cb5,cb6,cb7,cb8
      REAL :: an,ap3,ap1,bn,bp3,bp1
      REAL :: cnuei,kparvthe
      REAL :: nuei_p1_p1_1,nuei_p1_p3_1
      REAL :: nuei_u_u_1,nuei_u_q3_1
      REAL :: nuei_q1_u_1,nuei_q1_q1_1,nuei_q1_q3_1
      REAL :: nuei_q3_u_1,nuei_q3_q3_1
      REAL :: nuei_p1_p1_t,nuei_p1_p3_t
      REAL :: nuei_u_u_t,nuei_u_q3_t
      REAL :: nuei_q1_u_t,nuei_q1_q1_t,nuei_q1_q3_t
      REAL :: nuei_q3_u_t,nuei_q3_q3_t
      REAL :: nuei_n_n,nuei_n_p1,nuei_n_p3
      REAL :: nuei_p3_n,nuei_p3_p1,nuei_p3_p3
      REAL :: nuei_p1_n,nuei_p1_p3,nuei_p1_p1
      REAL :: nuei_u_u,nuei_u_q1,nuei_u_q3
      REAL :: nuei_q1_u,nuei_q1_q1,nuei_q1_q3
      REAL :: nuei_q3_u,nuei_q3_q1,nuei_q3_q3
      REAL :: k1,k2,k3,k4,k5
      REAL :: ki,ks0,charge_tot,gradne,gradne_s
      REAL :: beta2,bs
      REAL :: damp_psi,damp_sig
      COMPLEX :: phi_A,phi_B,phi_AU,phi_BU
      COMPLEX :: psi_A,psi_B,psi_AN,psi_BN  
      COMPLEX :: sig_A,sig_B
      COMPLEX :: k_par,k_par_psi
      REAL,ALLOCATABLE,DIMENSION(:) :: rwork
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: at,bt
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: vleft,vright
      COMPLEX,ALLOCATABLE,DIMENSION(:) :: work
      COMPLEX,ALLOCATABLE,DIMENSION(:) :: zomega
!
      ifail = 0
      lwork = 8*iur
      ALLOCATE(rwork(lwork))
      ALLOCATE(at(iur,iur))
      ALLOCATE(bt(iur,iur))
      ALLOCATE(vleft(iur,iur))
      ALLOCATE(vright(iur,iur))
      lwork = 33*iur
      ALLOCATE(work(lwork))
      ALLOCATE(zomega(iur))
!      write(*,*)"eigensolver allocation done"
!
      c35 = 3.0/5.0
      ft2 = ft*ft
      ft3 = ft*ft2
      ft4 = ft*ft3
      ft5 = ft*ft4
!      write(*,*)"R_unit = ",R_unit
!      write(*,*)"B_unit = ",B_unit
!      write(*,*)"q_unit = ",q_unit
!
       ky = ky_in
       k_par0 = park_in/(R_unit*q_unit*width_in)
       w_d0 = ky/R_unit
       w_cd = -gchat_in*w_d0
       w_s = -ky/B_unit
       wd_psi = w_cd
       modwd_psi = ABS(wd_psi)
       betae_psi = 0.0
       damp_psi = 0.0
       if(use_bper_in.eqv..FALSE.)then
           hnb0 = 0.0
           hp1b0 = 0.0
           hp3b0 = 0.0
           hr11b0 = 0.0
           hr13b0 = 0.0
           hr33b0 = 0.0
           hw113b0 = 0.0
           hw133b0 = 0.0
           hw333b0 = 0.0
           kpar_hp1b0 = 0.0
           kpar_hr11b0 = 0.0
           kpar_hr13b0 = 0.0
           wdhp1b0 = 0.0
           wdhr11b0 = 0.0
           wdhr13b0 = 0.0
           gnb0 = 0.0
           gp1b0 = 0.0
           gp3b0 = 0.0
           gr11b0 = 0.0
           gr13b0 = 0.0
           gr33b0 = 0.0
           gw113b0 = 0.0
           gw133b0 = 0.0
           gw333b0 = 0.0
           kpar_gp1b0 = 0.0
           kpar_gr11b0 = 0.0
           kpar_gr13b0 = 0.0
           wdgp1b0 = 0.0
           wdgr11b0 = 0.0
           wdgr13b0 = 0.0
       endif
       if(vpar_model_in.ne.0)then
           hnbp = 0.0
           hp1bp = 0.0
           hp3bp = 0.0
           hr11bp = 0.0
           hr13bp = 0.0
           hr33bp = 0.0
           hw113bp = 0.0
           hw133bp = 0.0
           hw333bp = 0.0
           wdhp1bp = 0.0
           wdhr11bp = 0.0
           wdhr13bp = 0.0
           kpar_hnbp = 0.0
           kpar_hp1bp = 0.0
           kpar_hp3bp = 0.0
           kpar_hr11bp = 0.0
           kpar_hr13bp = 0.0
           gnbp = 0.0
           gp1bp = 0.0
           gp3bp = 0.0
           gr11bp = 0.0
           gr13bp = 0.0
           gr33bp = 0.0
           gw113bp = 0.0
           gw133bp = 0.0
           gw333bp = 0.0
           wdgp1bp = 0.0
           wdgr11bp = 0.0
           wdgr13bp = 0.0
           kpar_gnbp = 0.0
           kpar_gp1bp = 0.0
           kpar_gp3bp = 0.0
           kpar_gr11bp = 0.0
           kpar_gr13bp = 0.0
       endif
       if(use_bpar_in.eqv..FALSE.)then
           h10n = 0.0
           h10p1 = 0.0
           h10p3 = 0.0
           h10r13 = 0.0
           h10r33 = 0.0
           g10n = 0.0
           g10p1 = 0.0
           g10p3 = 0.0
           g10r13 = 0.0
           g10r33 = 0.0
       endif
       if(use_bper_in)then
         betae_psi = 0.5*betae_in/(ky*ky)
!       write(*,*)"betae_psi = ",betae_psi
         damp_psi = damp_psi_in/MAX(betae_psi*vs(1)*vs(1),0.001)
       endif
       max_freq=2.0*ABS(ave_wdh(1,1))/R_unit
!       if(filter_in.gt.0.0)then
         do is=ns0,ns
           test = ABS(as(is)*zs(is)*(ave_hp3p0(is,1,1)*rlns(is)  &
                +1.5*(ave_hr13p0(is,1,1)-ave_hp3p0(is,1,1))*rlts(is)))
           max_freq = MAX(max_freq,test)
         enddo 
         max_freq = filter_in*ABS(ky)*max_freq
!         max_freq = 2.0*ABS(ky)*max_freq/ABS(as(1)*zs(1))
!         write(*,*)"max_freq = ",max_freq,filter_in
!       endif
       betae_sig = 0.0
       damp_sig = 0.0
       if(use_bpar_in)then
          betae_sig = 0.5*betae_in/(ky*ky)
          damp_sig = damp_sig_in/MAX(betae_sig*vs(1)*vs(1),0.001)
       endif
!       write(*,*)"betae_psi=",betae_psi,"betae_sig=",betae_sig
!       write(*,*)"damp_psi=",damp_psi,"damp_sig=",damp_sig
!       write(*,*)"k_par0=",k_par0,"w_d0=",w_d0,"w_cd=",w_cd,"w_s=",w_s
!
       Linsker = 0.5*Linsker_factor_in
       if(nbasis.eq.1)Linsker=0.0
!       write(*,*)"Linsker = ",Linsker
       am = 1.0
       bm = 0.0
!       if(iflagin_gf(5).eq.1)bm = 1.0
       if(Linsker.ne.0.0)then
         am = am/2.0
         bm = bm/2.0
       endif
!
!  GLF toroidal closure coefficients
!
       call get_v
       if(nroot.gt.6)then
        call get_u(ft)
         u2_r = u2_r*ft2
         u2_i = u2_i*ft2
         u3_r = u3_r/ft2
         u3_i = u3_i/ft2
         u5_r = u5_r*ft2
         u5_i = u5_i*ft2
         u7_r = u7_r*ft2
         u7_i = u7_i*ft2
         u9_r = u9_r/ft2
         u9_i = u9_i/ft2
         ub2_r = ub2_r*ft2
         ub2_i = ub2_i*ft2
         ub3_r = ub3_r/ft2
         ub3_i = ub3_i/ft2
         ub5_r = ub5_r*ft2
         ub5_i = ub5_i*ft2
         ub7_r = ub7_r*ft2
         ub7_i = ub7_i*ft2
         ub9_r = ub9_r/ft2
         ub9_i = ub9_i/ft2
       endif
!
! GLF parallel closure coefficients
!
       bpar_HP = 3.0+(32.0 - 9.0*pi)/(3.0*pi - 8.0)
       bper_DH = 1.0
! note sqrt(2) included in definition of d_per,d_par
       dper_DH = sqrt_two*sqrt_pi/2.0
       dpar_HP = 2.0*sqrt_two*sqrt_pi/(3.0*pi - 8.0)
       b1 = bpar_HP
       d1 = dpar_HP
       b3 = bper_DH
       d3 = dper_DH
       b33 = (b1 - b3)/(3.0)
       d33 = (d1 - d3)/(3.0)
!
!       if(vpar_shear_model_in.eq.1)then  
! include R(theta)/R0 factor like gyro convetions. Note that sign_Bt_in is in ave_c_tor_par
         do is=1,ns
           vpar_shear_s(is)=sign_It_in*alpha_p_in*vpar_shear_in(is)*ave_c_tor_par(1,1)/Rmaj_input
           if(vpar_model_in.eq.0)vpar_s(is) = alpha_mach_in*sign_It_in*vpar_in(is)*ave_c_tor_par(1,1)/Rmaj_input
         enddo
!       endif
!
!  compute electron-ion collsion model coefficients
!
       xnu_model = xnu_model_in
! different trapped/passing boundary layer models
! xnu_model = 0  version 1.80 large xnu limit = 0.0
! xnu_model = 1  version 1.81 used for APS07 , large xnu limit = adiabatic
! xnu_model = 2  version 1.85 large xnu_limit = circulating response 
! xnu_model = 3  retuned trapped boundary term to fit CGYRO with Lorentz operator 2/8/2017
       k1=0.0
       k2=0.0
       k3=0.0
       k4=0.0
       k5=0.0
       if(xnu_model.le.1)then
         k1 = 0.601248 + 1.12838*zeff_in
         k2 = 0.797885 + 1.12838*zeff_in
         k3 = 1.795240 + 2.25676*zeff_in
       endif
       xnu_p1_1 = -(4.0/5.0)*k2
       xnu_u_u_1 = -(2.0/3.0)*k1
       xnu_u_q3_1 = k1 - (2.0/5.0)*k2
       xnu_q1_u_1 = -(4.0/5.0)*k2
       xnu_q1_q1_1 = -(16.0/35.0)*k3
       xnu_q1_q3_1 = (6.0/5.0)*k2 + (12.0/35.0)*k3
       xnu_q3_u_1 = -(4.0/9.0)*k2    
       xnu_q3_q3_1 = (2.0/3.0)*k2 - (4.0/15.0)*k3
!       xnu_p1_1 = -(8.0/5.0)*k2
!       xnu_u_u_1 = -(4.0/3.0)*k1
!       xnu_u_q3_1 = 2.0*k1 - (4.0/5.0)*k2
!       xnu_q1_u_1 = -(8.0/5.0)*k2
!       xnu_q1_q1_1 = -(32.0/35.0)*k3
!       xnu_q1_q3_1 = (12.0/5.0)*k2 + (24.0/35.0)*k3
!       xnu_q3_u_1 = -(8.0/9.0)*k2
!       xnu_q3_q3_1 = (4.0/3.0)*k2 - (8.0/15.0)*k3
       if(park_in.eq.0.0)then
         xnu_u_u_1 = 0.0
         xnu_u_q3_1 =0.0
         xnu_q1_u_1=0.0
         xnu_q1_q1_1=0.0
         xnu_q1_q3_1=0.0
         xnu_q3_u_1=0.0
         xnu_q3_q3_1=0.0
       endif
       if(nroot.gt.6)then
!         xnu_bndry=(1.0 -ft2)* &
!          (1.96*ky/R_unit+0.60*xnue_in)/ &
!          (0.75*ky/R_unit + xnue_in)
         xnu_hat = xnue_in/(ky*taus(1)/R_unit)
!         xnu_bndry = (1.0 - ft2)*(((9.34*xnu_hat)**4)/(1+(9.34*xnu_hat)**4))* &
!                      (3.21 + 1.64*xnu_hat)/(1.0 + 2.40*xnu_hat)
!
! model for effective ion wavenumber averaged with ion charge densities zs*as
! ki = k_theta*sqrt(sum_s(rho_s**2*as*zs))
!
         ki=0.0
         charge_tot=0.0
         do is=2,ns
           ki = ki + taus(is)*mass(is)*as(is)*zs(is)
           charge_tot = charge_tot + as(is)*zs(is)
         enddo
         ki = SQRT(ki/charge_tot)*ky_in
         ks0 = ky*SQRT(taus(1)*mass(2))
!
!         write(*,*)"ki=",ki,SQRT(taus(2))*ky
         xnu_phi_b=0.0
         if(xnu_model.eq.1)xnu_phi_b=1.0
         if(xnu_phi_b.eq.0.0)then
! boundary collision model without phi terms
           gradne = rlns(1)*R_unit
           gradne_s = MAX(gradne+10.8,1.8)
           xnu_c = gradne_s*(1.5*(1.0-TANH((gradne_s/12.6)**2))+0.13)
!           write(*,*)"xnu_hat = ",xnu_hat,"xnu_c =",xnu_c
           xnu_a = ki*(MAX(0.36+0.10*gradne,0.0) + xnu_c*(ki/ks0)*(1.0-TANH(ki/0.55)))
           xnu_b = 3.1/(1.0+(2.1*ki + 8.0*ki*ki)*xnu_hat)
         else
! boundary collsion model with phi terms
           xnu_a = 0.41 + ((ki/ks0)**1.7)*0.70*(1.0 + 1.4*ki/0.38)/(1.0 + (ki/0.38)**4)
           xnu_b = 3.79/(1.0 + 4.63*ki*xnu_hat)
         endif
         xnu_bndry = (1.0 - ft2)*xnu_a*xnu_b
!         write(*,*)"xnu_bndry =",xnu_a*xnu_b
!         xnu_bndry = (1.0 - ft2)
         xnu_n_b = xnu_factor_in*xnu_bndry*xnu_p1_1
         xnu_p3_b =  xnu_n_b
         xnu_p1_b = xnu_n_b
         xnu_u_b = xnu_factor_in*xnu_bndry*xnu_q1_q1_1 
         xnu_q3_b = xnu_u_b
         xnu_q1_b = xnu_u_b
          if(park_in.eq.0.0)then
           xnu_u_b = 0.0
           xnu_q1_b=0.0
           xnu_q3_b =0.0
        endif
      endif
! debug
!      write(*,*)"check xnu"
!      write(*,*)xnu_p1_1,xnu_u_u_1 
!      write(*,*)xnu_u_q3_1,xnu_q3_u_1
!      write(*,*)xnu_q3_q3_1,xnu_q1_u_1
!      write(*,*)xnu_q1_q3_1,xnu_q1_q1_1
!      write(*,*)xnu_n_b,xnu_p3_b
!      write(*,*)xnu_p1_b,xnu_u_b
!      write(*,*)xnu_q3_b,xnu_q1_b
!
      cnuei = 0.0
      if(xnu_model.eq.2.or.xnu_model.eq.3)cnuei = xnue_in
!      kparvthe = ABS(k_par0)*vs(1)*xnu_factor_in/sqrt_two
      kparvthe = ABS(k_par0)*vs(1)/sqrt_two
      kparvthe=MAX(kparvthe,1.0E-10)
!      write(*,*)"cnuei=",cnuei,"kparvthe=",kparvthe,"ft=",ft
!
! full velocity space terms 
!
      k1 = cnuei*(nuei_c1(1) +zeff_in*nuei_c1(2))
      k2 = cnuei*(nuei_c1(3) +zeff_in*nuei_c1(4))
      k3 = cnuei*(nuei_c1(5) +zeff_in*nuei_c1(6))
      k4 = cnuei*(nuei_c1(7) +zeff_in*nuei_c1(8))
      k5 = cnuei*(nuei_c1(9) +zeff_in*nuei_c1(10))
      c01 = (4.0/5.0)*k4
      c02 = (2.0/5.0)*k2 - k1
      c03 = (2.0/3.0)*k1
      c04 = (4.0/15.0)*k3-(4.0/3.0)*k2+(5.0/3.0)*k1
      c05 = (16.0/35.0)*k5
      c06 = 0.0
      c07 = 0.0
      c08 = 0.0
      c09 = 0.0
      c010 = 0.0

!      write(*,*)"nuei_c1"
!      do i=1,nc1
!       write(*,*)i,nuei_c1(i)
!      enddo
!      write(*,*)"c0=",c01,c02,c03,c04,c05,c06,c07,c08,c09,c010
      nuei_p1_p1_1 = c01
      nuei_p1_p3_1 = -nuei_p1_p1_1
      nuei_u_q3_1 = c02
      nuei_u_u_1 = c03 -(5.0/3.0)*nuei_u_q3_1
      nuei_q3_q3_1 = c04
      nuei_q3_u_1 = (10.0/9.0)*c02 -(5.0/3.0)*nuei_q3_q3_1
      nuei_q3_u_1 = nuei_q3_u_1 + (5.0/3.0)*nuei_u_u_1
      nuei_q3_q3_1 = nuei_q3_q3_1 + (5.0/3.0)*nuei_u_q3_1
      nuei_q1_q1_1 = c05
      nuei_q1_q3_1 = -(9.0/5.0)*nuei_q1_q1_1
      nuei_q1_u_1 = (9.0/5.0)*nuei_q3_u_1
      nuei_q1_q3_1 = nuei_q1_q3_1 + (9.0/5.0)*nuei_q3_q3_1
! debug compare with xnu_*_1 terms
!      write(*,*)"check nuei_1"
!      write(*,*)nuei_p1_p1_1,nuei_u_u_1+(5./3.)*nuei_u_q3_1
!      write(*,*)nuei_u_q3_1,nuei_q3_u_1+(5./3.)*nuei_q3_q3_1
!      write(*,*)nuei_q3_q3_1,nuei_q1_u_1+3.0*nuei_q3_q3_1
!      write(*,*)nuei_q1_q3_1,nuei_q1_q1_1
!
      nuei_p1_p1_t = c01
      nuei_p1_p3_t = -ft2*nuei_p1_p1_t
      nuei_u_q3_t = c02
      nuei_u_u_t = c03 -(5.0/3.0)*nuei_u_q3_t
      nuei_q3_q3_t = c04
      nuei_q3_u_t = (10.0/9.0)*c02 -(5.0/3.0)*nuei_q3_q3_t
      nuei_q3_u_t = nuei_q3_u_t + (5.0/3.0)*nuei_u_u_t
      nuei_q3_q3_t = nuei_q3_q3_t + (5.0/3.0)*nuei_u_q3_t
      nuei_q1_q1_t = c05
      nuei_q1_q3_t = -(9.0/5.0)*ft2*nuei_q1_q1_t
      nuei_q1_u_t = (9.0/5.0)*ft2*nuei_q3_u_t
      nuei_q1_q3_t = nuei_q1_q3_t + (9.0/5.0)*ft2*nuei_q3_q3_t
! debug compare with xnu_*_1 terms
!      write(*,*)"check nuei_1"
!      write(*,*)nuei_p1_p1_t,nuei_u_u_t+(5./3.)*nuei_u_q3_t
!      write(*,*)nuei_u_q3_t,nuei_q3_u_t+(5./3.)*nuei_q3_q3_t
!      write(*,*)nuei_q3_q3_t,nuei_q1_u_t+3.0*nuei_q3_q3_t
!      write(*,*)nuei_q1_q3_t,nuei_q1_q1_t
!
       an = 0.75 
       ap3 = 1.25 
       ap1 = 2.25 
       bn = ft
       bp3 = ft
       bp1 = ft3
      cb3=0.0
      cb5 =0.0
!recalibrated 8/20/14      cb1 = 0.114*SQRT(kparvthe*cnuei*(1.0 + 0.82*zeff_in))
      cb1 = 0.163*SQRT(kparvthe*cnuei*(1.0 + 0.82*zeff_in))
      if(xnu_model_in.eq.3)cb1 = 0.50*(kparvthe**0.34)*(cnuei*(1.0 + 0.82*zeff_in))**0.66
      cb1 = cb1*xnu_factor_in
      cb2 = cb1
      cb4 = cb1
! even trapped region terms
      nuei_n_n = (1.0 -ft2)*cb1
!      nuei_n_n = (1.0 -ft2)*cnuei*cb1
      nuei_n_p3 = (1.0 - ft2)*cnuei*cb3
!      nuei_n_p1 = cb4*(1.0 -ft2)/MAX(ft4,1.D-4)
      nuei_n_p1 = 0.0
      nuei_n_p3 = nuei_n_p3 - ft2*nuei_n_p1
      nuei_n_n = nuei_n_n -nuei_n_p3 - ft2*nuei_n_p1
      nuei_p3_n = (2.0/3.0)*(1.0 - ft2)*cnuei*cb3
      nuei_p3_p3 = (1.0 -ft2)*cb2
!      nuei_p3_p3 = (1.0 -ft2)*cnuei*cb2
!      nuei_p3_p1 = (1.0 -ft2)*cb5/MAX(ft4,1.D-4)
      nuei_p3_p1 = 0.0
      nuei_p3_p3 = nuei_p3_p3 - ft2*nuei_p3_p1
      nuei_p3_n = nuei_p3_n -nuei_p3_p3 -ft2*nuei_p3_p1
      nuei_p3_n = nuei_p3_n +nuei_n_n
      nuei_p3_p3 = nuei_p3_p3 +nuei_n_p3
      nuei_p3_p1 = nuei_p3_p1 +nuei_n_p1
      nuei_p1_n = 0.0
      nuei_p1_p3 = 0.0
      nuei_p1_p1 = (1.0 - ft2)*cb4
!      nuei_p1_p1 = (1.0 - ft2)*cnuei*cb4
!      nuei_p1_p1=0.0
      nuei_p1_p3 = nuei_p1_p3 -ft2*nuei_p1_p1
      nuei_p1_n = nuei_p1_n -nuei_p1_p3 -ft2*nuei_p1_p1
      nuei_p1_n = nuei_p1_n + ft2*nuei_p3_n
      nuei_p1_p3 = nuei_p1_p3 +ft2*nuei_p3_p3
      nuei_p1_p1 = nuei_p1_p1 +ft2*nuei_p3_p1      
! odd trapped region terms
      cb5=0.0
      cb6=0.0
      cb7=0.0
      cb8=0.0
      nuei_u_u = (1.0 -ft2)*cnuei*cb5
      nuei_u_q3 = (1.0 -ft2)*cnuei*cb7
!      nuei_u_q1 = cb5*(1.0 -ft2)/MAX(ft4,1.D-4)
      nuei_u_q1 = 0.0
      nuei_u_q3 = nuei_u_q3 - (9.0/5.0)*ft2*nuei_u_q1
      nuei_u_u = nuei_u_u -(5.0/3.0)*nuei_u_q3 -3.0*ft2*nuei_u_q1
      nuei_q3_u = (10.0/9.0)*(1.0 -ft2)*cnuei*cb7
      nuei_q3_q3 = (1.0 - ft2)*cnuei*cb6
!      nuei_q3_q1 = cb6*(1.0 -ft2)/MAX(ft4,1.D-4)
      nuei_q3_q1 = 0.0
      nuei_q3_q3 = nuei_q3_q3 -(9.0/5.0)*ft2*nuei_q3_q1
      nuei_q3_u = nuei_q3_u - (5.0/3.0)*nuei_q3_q3 -3.0*ft2*nuei_q3_q1
      nuei_q3_u = nuei_q3_u + (5.0/3.0)*nuei_u_u
      nuei_q3_q3 = nuei_q3_q3 +(5.0/3.0)*nuei_u_q3
      nuei_q3_q1 = nuei_q3_q1 +(5.0/3.0)*nuei_u_q1
      nuei_q1_u = 0.0
      nuei_q1_q3 = 0.0
      nuei_q1_q1 = (1.0 -ft2)*cnuei*cb8
      nuei_q1_q3 = nuei_q1_q3 - (9.0/5.0)*ft2*nuei_q1_q1
      nuei_q1_u = nuei_q1_u - (5.0/3.0)*nuei_q1_q3 - 3.0*ft2*nuei_q1_q1
      nuei_q1_u = nuei_q1_u +(9.0/5.0)*ft2*nuei_q3_u
      nuei_q1_q3 = nuei_q1_q3 +(9.0/5.0)*ft2*nuei_q3_q3
      nuei_q1_q1 = nuei_q1_q1 + (9.0/5.0)*ft2*nuei_q3_q1

! debug
!      write(*,*)"check nuei",ft
!      write(*,*)nuei_n_n,nuei_n_p3,nuei_n_p1
!      write(*,*)nuei_p3_n,nuei_p3_p3,nuei_p3_p1
!      write(*,*)nuei_p1_n,nuei_p1_p3,nuei_p1_p1
!      write(*,*)nuei_u_u,nuei_u_q3,nuei_u_q1
!      write(*,*)nuei_q3_u,nuei_q3_q3,nuei_q3_q1
!      write(*,*)nuei_q1_u,nuei_q1_q3,nuei_q1_q1
!
! done with electron-ion collision model
!
!
! start of loop over species is,js for amat
!
      do is = ns0,ns
!
      do js = ns0,ns
!
! start of loop over basis ib,jb for amat
!
       do ib = 1,nbasis
       do jb = 1,nbasis
! trick to fool compiler into not messing with these loops
         if(idum.eq.1)write(*,*)"fooled you"
!        if(js.eq.ns0)then
!         write(*,*)"is=",is,"js=",js,"ib=",ib,"jb=",jb
!
! collision  model
!
        xnuei = 0.0
        if(is.eq.1)then
          xnuei = xnue_in
        endif
        xnuion = 0.0
        d_ab=0.0
        d_ij=0.0
        d_ee=0.0
        d_ab_psi=0.0
        d_11 = 0.0
        d_1 = 0.0
        if(is.eq.1)d_1=1.0
        if(is.eq.1.and.js.eq.1)d_11=1.0
        if(is.eq.js)d_ij=1.0
        if(ib.eq.jb.and.is.eq.js)d_ab=1.0
        if(is.eq.1.and.d_ab.eq.1.0)d_ee=1.0
        if(ib.eq.jb)d_ab_psi=1.0
        b1 = bpar_HP
        d1 = dpar_HP
        b3 = bper_DH
        d3 = dper_DH
        bs = 20.0*SQRT(taus(is)*mass(is))*ky/ABS(zs(is))
!        b3 = MIN(bs,1.0)*bper_DH
!        d3 = MIN(bs,1.0)*dper_DH
        b33 = (b1 - b3)/(3.0)
        d33 = (d1 - d3)/(3.0)
!
         hn = ave_hnp0(is,ib,jb)
         hp1 = ave_hp1p0(is,ib,jb)
         hp3 = ave_hp3p0(is,ib,jb)
         hr11 = ave_hr11p0(is,ib,jb)
         hr13 = ave_hr13p0(is,ib,jb)
         hr33 = ave_hr33p0(is,ib,jb)
         c_tor_par_hp1 = ave_c_tor_par_hp1p0(is,ib,jb)/Rmaj_input
         c_tor_par_hr11 = ave_c_tor_par_hr11p0(is,ib,jb)/Rmaj_input
         c_tor_par_hr13 = ave_c_tor_par_hr13p0(is,ib,jb)/Rmaj_input
!        write(*,*)is,ib,jb
!        write(*,*)hn,hp1,hp3
!        write(*,*)hr11,hr13,hr33
         if(use_bper_in)then
           hnb0 = ave_hnb0(is,ib,jb)
           hp1b0 = ave_hp1b0(is,ib,jb)
           hp3b0 = ave_hp3b0(is,ib,jb)
           hr11b0 = ave_hr11b0(is,ib,jb)
           hr13b0 = ave_hr13b0(is,ib,jb)
           hr33b0 = ave_hr33b0(is,ib,jb)
           hw113b0 = ave_hw113b0(is,ib,jb)
           hw133b0 = ave_hw133b0(is,ib,jb)
           hw333b0 = ave_hw333b0(is,ib,jb)
           if(vpar_model_in.eq.0)then
            hnbp = ave_hnbp(is,ib,jb)
            hp1bp = ave_hp1bp(is,ib,jb)
            hp3bp = ave_hp3bp(is,ib,jb)
            hr11bp = ave_hr11bp(is,ib,jb)
            hr13bp = ave_hr13bp(is,ib,jb)
            hr33bp = ave_hr33bp(is,ib,jb)
            hw113bp = ave_hw113bp(is,ib,jb)
            hw133bp = ave_hw133bp(is,ib,jb)
            hw333bp = ave_hw333bp(is,ib,jb)
           endif
         endif
         if(use_bpar_in)then
           h10n = 1.5*(hnb0-hp3b0)
           h10p1 = 2.5*hp1b0 - 1.5*hr13b0
           h10p3 = 2.5*hp3b0 - 1.5*hr33b0
           h10r13 = 3.5*hr13b0 - 1.5*hw133b0
           h10r33 = 3.5*hr33b0 - 1.5*hw333b0
         endif
!           write(*,*)"check hb0",is,js,ib,jb
!           write(*,*)h10n,h10p1,h10p3
!           write(*,*)h10r13,h10r33
!           write(*,*)hnb0,hp1b0,hp3b0
!           write(*,*)hr11b0,hr13b0,hr33b0
!           write(*,*)hw113b0,hw133b0,hw333b0
         hu1 = ave_hu1(is,ib,jb)
         hu3 = ave_hu3(is,ib,jb)
         ht1 = ave_ht1(is,ib,jb)
         ht3 = ave_ht3(is,ib,jb)
!         write(*,*)hu1,hu3,ht1,ht3
         wdhp1p0 = ave_wdhp1p0(is,ib,jb)
         wdhr11p0 = ave_wdhr11p0(is,ib,jb)
         wdhr13p0 = ave_wdhr13p0(is,ib,jb)
         if(use_bper_in)then
           wdhp1b0 = ave_wdhp1b0(is,ib,jb)
           wdhr11b0 = ave_wdhr11b0(is,ib,jb)
           wdhr13b0 = ave_wdhr13b0(is,ib,jb)
           if(vpar_model_in.eq.0)then
             wdhp1bp = ave_wdhp1bp(is,ib,jb)
             wdhr11bp = ave_wdhr11bp(is,ib,jb)
             wdhr13bp = ave_wdhr13bp(is,ib,jb)
           endif
         endif
         wdhu1 = ave_wdhu1(is,ib,jb)
         wdhu3 = ave_wdhu3(is,ib,jb)
         wdhu3ht1 = ave_wdhu3ht1(is,ib,jb)
         wdhu3ht3 = ave_wdhu3ht3(is,ib,jb)
         wdhu33 = ave_wdhu33(is,ib,jb)
         wdhu33ht1 = ave_wdhu33ht1(is,ib,jb)
         wdhu33ht3 = ave_wdhu33ht3(is,ib,jb)
         modwdhu1 = ave_modwdhu1(is,ib,jb)
         modwdhu3 = ave_modwdhu3(is,ib,jb)
         modwdhu3ht1 = ave_modwdhu3ht1(is,ib,jb)
         modwdhu3ht3 = ave_modwdhu3ht3(is,ib,jb)
         modwdhu33 = ave_modwdhu33(is,ib,jb)
         modwdhu33ht1 = ave_modwdhu33ht1(is,ib,jb)
         modwdhu33ht3 = ave_modwdhu33ht3(is,ib,jb)
         hv1r = (v1_r-vb1_r)*ave_modwdh(ib,jb)+vb1_r*modwdhu3*c35
         hv2r = (v2_r-vb2_r)*ave_modwdh(ib,jb)+vb2_r*modwdhu3*c35
         hv3r = (v3_r-vb3_r)*ave_modwdh(ib,jb)+vb3_r*modwdhu33*c35
         hv4r = (v4_r-vb4_r)*ave_modwdh(ib,jb)+vb4_r*modwdhu33*c35
         hv5r = (v5_r-vb5_r)*ave_modwdh(ib,jb)+vb5_r*modwdhu3*c35
         hv6r = (v6_r-vb6_r)*ave_modwdh(ib,jb)+vb6_r*modwdhu3*c35
         hv7r = (v7_r-vb7_r)*ave_modwdh(ib,jb)+vb7_r*modwdhu3*c35
         hv8r = (v8_r-vb8_r)*ave_modwdh(ib,jb)+vb8_r*modwdhu33*c35
         hv9r = (v9_r-vb9_r)*ave_modwdh(ib,jb)+vb9_r*modwdhu33*c35
         hv10r = (v10_r-vb10_r)*ave_modwdh(ib,jb) &
          +vb10_r*modwdhu33*c35
         hv1rht1 = (v1_r-vb1_r)*ave_modwdht1(is,ib,jb) &
          +vb1_r*modwdhu3ht1*c35
         hv2rht3 = (v2_r-vb2_r)*ave_modwdht3(is,ib,jb) &
          +vb2_r*modwdhu3ht3*c35
         hv3rht1 = (v3_r-vb3_r)*ave_modwdht1(is,ib,jb) &
          +vb3_r*modwdhu33ht1*c35
         hv4rht3 = (v4_r-vb4_r)*ave_modwdht3(is,ib,jb) &
          +vb4_r*modwdhu33ht3*c35
         hv6rhu1 = (v6_r*modwdhu1)
         hv7rhu3 = (v7_r*modwdhu3)
         hv9rhu1 = (v9_r*modwdhu1)
         hv10rhu3= (v10_r*modwdhu3)
         hv1i = (v1_i-vb1_i)*ave_wdh(ib,jb)+vb1_i*wdhu3*c35
         hv2i = (v2_i-vb2_i)*ave_wdh(ib,jb)+vb2_i*wdhu3*c35
         hv3i = (v3_i-vb3_i)*ave_wdh(ib,jb)+vb3_i*wdhu33*c35
         hv4i = (v4_i-vb4_i)*ave_wdh(ib,jb)+vb4_i*wdhu33*c35
         hv5i = (v5_i-vb5_i)*ave_wdh(ib,jb)+vb5_i*wdhu3*c35
         hv6i = (v6_i-vb6_i)*ave_wdh(ib,jb)+vb6_i*wdhu3*c35
         hv7i = (v7_i-vb7_i)*ave_wdh(ib,jb)+vb7_i*wdhu3*c35
         hv8i = (v8_i-vb8_i)*ave_wdh(ib,jb)+vb8_i*wdhu33*c35
         hv9i = (v9_i-vb9_i)*ave_wdh(ib,jb)+vb9_i*wdhu33*c35
         hv10i = (v10_i-vb10_i)*ave_wdh(ib,jb)+vb10_i*wdhu33*c35
         hv1iht1 = (v1_i-vb1_i)*ave_wdht1(is,ib,jb) &
          +vb1_i*wdhu3ht1*c35
         hv2iht3 = (v2_i-vb2_i)*ave_wdht3(is,ib,jb) &
          +vb2_i*wdhu3ht3*c35
         hv3iht1 = (v3_i-vb3_i)*ave_wdht1(is,ib,jb) &
          +vb3_i*wdhu33ht1*c35
         hv4iht3 = (v4_i-vb4_i)*ave_wdht3(is,ib,jb) &
          +vb4_i*wdhu33ht3*c35
         hv6ihu1 = (v6_i*wdhu1)
         hv7ihu3 = (v7_i*wdhu3)
         hv9ihu1 = (v9_i*wdhu1)
         hv10ihu3= (v10_i*wdhu3)

! debug
!         write(*,*)"hv1r = ",hv1r
!         write(*,*)"hv2r = ",hv2r
!         write(*,*)"hv3r = ",hv3r
!         write(*,*)"hv4r = ",hv4r
!         write(*,*)"hv5r = ",hv5r
!         write(*,*)"hv6r = ",hv6r
!         write(*,*)"hv7r = ",hv7r
!         write(*,*)"hv8r = ",hv8r
!         write(*,*)"hv9r = ",hv9r
!         write(*,*)"hv10r = ",hv10r
!         write(*,*)"hv1i = ",hv1i
!         write(*,*)"hv2i = ",hv2i
!         write(*,*)"hv3i = ",hv3i
!         write(*,*)"hv4i = ",hv4i
!         write(*,*)"hv5i = ",hv5i
!         write(*,*)"hv6i = ",hv6i
!         write(*,*)"hv7i = ",hv7i
!         write(*,*)"hv8i = ",hv8i
!         write(*,*)"hv9i = ",hv9i
!         write(*,*)"hv10i = ",hv10i
!
!          write(*,*)is,ib,jb
!          write(*,*)"ave_kpar_eff =",ave_kpar_eff(is,ib,jb)
!          write(*,*)"ave_modkpar_eff =",ave_modkpar_eff(is,ib,jb)
!
          kpar_hnp0 = k_par0*ave_kparhnp0(is,ib,jb)
          kpar_hp1p0 = k_par0*ave_kparhp1p0(is,ib,jb)
          kpar_hp3p0 = k_par0*ave_kparhp3p0(is,ib,jb)
          if(use_bper_in)then
            kpar_hp1b0 = k_par0*ave_kparhp1b0(is,ib,jb)
            kpar_hr11b0 = k_par0*ave_kparhr11b0(is,ib,jb)
            kpar_hr13b0 = k_par0*ave_kparhr13b0(is,ib,jb)
            if(vpar_model_in.eq.0)then
              kpar_hnbp = k_par0*ave_kparhnbp(is,ib,jb)
              kpar_hp1bp = k_par0*ave_kparhp1bp(is,ib,jb)
              kpar_hp3bp = k_par0*ave_kparhp3bp(is,ib,jb)
              kpar_hr11bp = k_par0*ave_kparhr11bp(is,ib,jb)
              kpar_hr13bp = k_par0*ave_kparhr13bp(is,ib,jb)
            endif
          endif
          kpar_hu1 = ave_kparhu1(is,ib,jb)
          kpar_hu3 = ave_kparhu3(is,ib,jb)
          kpar_hb1 = b1*ave_kpar_eff(is,ib,jb)
          kpar_hb3 = b3*ave_kpar_eff(is,ib,jb)
          kpar_hb33 = b33*ave_kpar_eff(is,ib,jb)
          kpar_hb1ht1 = b1*ave_kparht1(is,ib,jb)
          kpar_hb3ht3 = b3*ave_kparht3(is,ib,jb)
          kpar_hb33ht1 = b33*ave_kparht1(is,ib,jb)
          modkpar_hd1 = d1*ave_modkpar_eff(is,ib,jb)
          modkpar_hd3 = d3*ave_modkpar_eff(is,ib,jb)
          modkpar_hd33 = d33*ave_modkpar_eff(is,ib,jb)
          modkpar_hd1hu1 = d1*ave_modkparhu1(is,ib,jb)
          modkpar_hd3hu3 = d3*ave_modkparhu3(is,ib,jb)
          modkpar_hd33hu1 = d33*ave_modkparhu1(is,ib,jb)
!          grad_hu1 = ave_gradhu1(is,ib,jb)
!          grad_hu3 = ave_gradhu3(is,ib,jb)
          grad_hu1 = 0.0
          grad_hu3 = 0.0
          dhr13 = -b3*(ave_kparht3(is,ib,jb)-ave_kparht1(is,ib,jb)/3.0 &
              -ave_kparhu3(is,ib,jb) + ave_kparhu1(is,ib,jb)/3.0)
!
         if(Linsker.eq.0.0)then
           gradhp1=0.0
           gradhr11=0.0
           gradhr13=0.0
           gradhp1p1=0.0
           gradhr11p1=0.0
           gradhr13p1=0.0
         else
          gradhp1 = Linsker*ave_gradhp1p0(is,ib,jb)
          gradhr11 = Linsker*ave_gradhr11p0(is,ib,jb)
          gradhr13 = Linsker*ave_gradhr13p0(is,ib,jb)
          gradhp1p1 = Linsker*ave_gradhp1p1(is,ib,jb)
          gradhr11p1 = Linsker*ave_gradhr11p1(is,ib,jb)
          gradhr13p1 = Linsker*ave_gradhr13p1(is,ib,jb) 
        endif    
!
        if(nroot.gt.6)then
         gn = ave_gnp0(is,ib,jb)
         gp1 = ave_gp1p0(is,ib,jb)
         gp3 = ave_gp3p0(is,ib,jb)
         gr11 = ave_gr11p0(is,ib,jb)
         gr13 = ave_gr13p0(is,ib,jb)
         gr33 = ave_gr33p0(is,ib,jb)
         c_tor_par_gp1 = ave_c_tor_par_gp1p0(is,ib,jb)/Rmaj_input
         c_tor_par_gr11 = ave_c_tor_par_gr11p0(is,ib,jb)/Rmaj_input
         c_tor_par_gr13 = ave_c_tor_par_gr13p0(is,ib,jb)/Rmaj_input
!         c_tor_par_gp1 = B2_ave_out*gp1
!         c_tor_par_gr11 = B2_ave_out*gr11
!         c_tor_par_gr13 = B2_ave_out*gr13
!        write(*,*)is,js,ib,jb
!        write(*,*)gn,gp1,gp3
!        write(*,*)gr11,gr13,gr33
         if(use_bper_in)then
           gnb0 = ave_gnb0(is,ib,jb)
           gp1b0 = ave_gp1b0(is,ib,jb)
           gp3b0 = ave_gp3b0(is,ib,jb)
           gr11b0 = ave_gr11b0(is,ib,jb)
           gr13b0 = ave_gr13b0(is,ib,jb)
           gr33b0 = ave_gr33b0(is,ib,jb)
           gw113b0 = ave_gw113b0(is,ib,jb)
           gw133b0 = ave_gw133b0(is,ib,jb)
           gw333b0 = ave_gw333b0(is,ib,jb)
           if(vpar_model_in.eq.0)then
            gnbp = ave_gnbp(is,ib,jb)
            gp1bp = ave_gp1bp(is,ib,jb)
            gp3bp = ave_gp3bp(is,ib,jb)
            gr11bp = ave_gr11bp(is,ib,jb)
            gr13bp = ave_gr13bp(is,ib,jb)
            gr33bp = ave_gr33bp(is,ib,jb)
            gw113bp = ave_gw113bp(is,ib,jb)
            gw133bp = ave_gw133bp(is,ib,jb)
            gw333bp = ave_gw333bp(is,ib,jb)
           endif
         endif
         if(use_bpar_in)then
           g10n = 1.5*(gnb0-gp3b0)
           g10p1 = 2.5*gp1b0 - 1.5*gr13b0
           g10p3 = 2.5*gp3b0 - 1.5*gr33b0
           g10r13 = 3.5*gr13b0 - 1.5*gw133b0
           g10r33 = 3.5*gr33b0 - 1.5*gw333b0
         endif
!           write(*,*)"check gb0",ft,is,js,ib,jb
!           write(*,*)g10n,g10p1,g10p3
!           write(*,*)g10r13,g10r33
!           write(*,*)gnb0,gp1b0,gp3b0
!           write(*,*)gr11b0,gr13b0,gr33b0
!           write(*,*)gw113b0,gw133b0,gw333b0
         gu1 = ave_gu1(is,ib,jb)
         gu3 = ave_gu3(is,ib,jb)
         gt1 = ave_gt1(is,ib,jb)
         gt3 = ave_gt3(is,ib,jb)
!         write(*,*)gu1,gu3,gt1,gt3
!         gu1 = hu1*ft2
!         gu3 = hu3
!         gt1 = ht1*ft2
!         gt3 = ht3
         wdgp1p0 = ave_wdgp1p0(is,ib,jb)
         wdgr11p0 = ave_wdgr11p0(is,ib,jb)
         wdgr13p0 = ave_wdgr13p0(is,ib,jb)
         if(use_bper_in)then
           wdgp1b0 = ave_wdgp1b0(is,ib,jb)
           wdgr11b0 = ave_wdgr11b0(is,ib,jb)
           wdgr13b0 = ave_wdgr13b0(is,ib,jb)
           if(vpar_model_in.eq.0)then
            wdgp1bp = ave_wdgp1bp(is,ib,jb)
            wdgr11bp = ave_wdgr11bp(is,ib,jb)
            wdgr13bp = ave_wdgr13bp(is,ib,jb)
           endif
         endif
         wdgu1 = ave_wdgu1(is,ib,jb)
         wdgu3 = ave_wdgu3(is,ib,jb)
         wdgu33 = ave_wdgu33(is,ib,jb)
         wdgu3gt1 = ave_wdgu3gt1(is,ib,jb)
         wdgu3gt3 = ave_wdgu3gt3(is,ib,jb)
         wdgu33gt1 = ave_wdgu33gt1(is,ib,jb)
         wdgu33gt3 = ave_wdgu33gt3(is,ib,jb)
         modwdgu1 = ave_modwdgu1(is,ib,jb)
         modwdgu3 = ave_modwdgu3(is,ib,jb)
         modwdgu33 = ave_modwdgu33(is,ib,jb)
         modwdgu3gt1 = ave_modwdgu3gt1(is,ib,jb)
         modwdgu3gt3 = ave_modwdgu3gt3(is,ib,jb)
         modwdgu33gt1 = ave_modwdgu33gt1(is,ib,jb)
         modwdgu33gt3 = ave_modwdgu33gt3(is,ib,jb)
         gu1r = (u1_r-ub1_r)*ave_modwdg(ib,jb)+ub1_r*modwdgu3*c35
         gu2r = (u2_r-ub2_r)*ave_modwdg(ib,jb)+ub2_r*modwdgu3*c35
         gu3r = (u3_r-ub3_r)*ave_modwdg(ib,jb)+ub3_r*modwdgu33*c35
         gu4r = (u4_r-ub4_r)*ave_modwdg(ib,jb)+ub4_r*modwdgu33*c35
         gu5r = (u5_r-ub5_r)*ave_modwdg(ib,jb)+ub5_r*modwdgu3*c35
         gu6r = (u6_r-ub6_r)*ave_modwdg(ib,jb)+ub6_r*modwdgu3*c35
         gu7r = (u7_r-ub7_r)*ave_modwdg(ib,jb)+ub7_r*modwdgu3*c35
         gu8r = (u8_r-ub8_r)*ave_modwdg(ib,jb)+ub8_r*modwdgu33*c35
         gu9r = (u9_r-ub9_r)*ave_modwdg(ib,jb)+ub9_r*modwdgu33*c35
         gu10r = (u10_r-ub10_r)*ave_modwdg(ib,jb) &
                 +ub10_r*modwdgu33*c35
         gu1rgt1 = (u1_r-ub1_r)*ave_modwdgt1(is,ib,jb) &
                 +ub1_r*modwdgu3gt1*c35
         gu2rgt3 = (u2_r-ub2_r)*ave_modwdgt3(is,ib,jb) &
                 +ub2_r*modwdgu3gt3*c35
         gu3rgt1 = (u3_r-ub3_r)*ave_modwdgt1(is,ib,jb) &
                 +ub3_r*modwdgu33gt1*c35
         gu4rgt3 = (u4_r-ub4_r)*ave_modwdgt3(is,ib,jb) &
                 +ub4_r*modwdgu33gt3*c35
         gu6rgu1 = (u6_r*modwdgu1)
         gu7rgu3 = (u7_r*modwdgu3)
         gu9rgu1 = (u9_r*modwdgu1)
         gu10rgu3= (u10_r*modwdgu3)
         gu1i = (u1_i-ub1_i)*ave_wdg(ib,jb)+ub1_i*wdgu3*c35
         gu2i = (u2_i-ub2_i)*ave_wdg(ib,jb)+ub2_i*wdgu3*c35
         gu3i = (u3_i-ub3_i)*ave_wdg(ib,jb)+ub3_i*wdgu33*c35
         gu4i = (u4_i-ub4_i)*ave_wdg(ib,jb)+ub4_i*wdgu33*c35
         gu5i = (u5_i-ub5_i)*ave_wdg(ib,jb)+ub5_i*wdgu3*c35
         gu6i = (u6_i-ub6_i)*ave_wdg(ib,jb)+ub6_i*wdgu3*c35
         gu7i = (u7_i-ub7_i)*ave_wdg(ib,jb)+ub7_i*wdgu3*c35
         gu8i = (u8_i-ub8_i)*ave_wdg(ib,jb)+ub8_i*wdgu33*c35
         gu9i = (u9_i-ub9_i)*ave_wdg(ib,jb)+ub9_i*wdgu33*c35
         gu10i = (u10_i-ub10_i)*ave_wdg(ib,jb) &
                 +ub10_i*wdgu33*c35
         gu1igt1 = (u1_i-ub1_i)*ave_wdgt1(is,ib,jb) &
                 +ub1_i*wdgu3gt1*c35
         gu2igt3 = (u2_i-ub2_i)*ave_wdgt3(is,ib,jb) &
                 +ub2_i*wdgu3gt3*c35
         gu3igt1 = (u3_i-ub3_i)*ave_wdgt1(is,ib,jb) &
                 +ub3_i*wdgu33gt1*c35
         gu4igt3 = (u4_i-ub4_i)*ave_wdgt3(is,ib,jb) &
                 +ub4_i*wdgu33gt3*c35
         gu6igu1 = (u6_i*wdgu1)
         gu7igu3 = (u7_i*wdgu3)
         gu9igu1 = (u9_i*wdgu1)
         gu10igu3= (u10_i*wdgu3)
! debug
!         write(*,*)"c35*wdgu3-wd",c35*wdgu3-ave_wdg(ib,jb)
!         write(*,*)"c35*wdgu33-wd",c35*wdgu33-ave_wdg(ib,jb)
!         write(*,*)"c35*modwdgu3-modwd",c35*modwdgu3-ave_modwd(ib,jb)
!         write(*,*)"c35*modwdgu33-modwd",c35*modwdgu33-ave_modwd(ib,jb)
!         write(*,*)"ave_modwd(1,1)=",ave_modwd(1,1)
!         write(*,*)"c35*gu3-1",c35*ave_gu3(is,ib,jb)-1.0
!         write(*,*)"c35*gu33-1",c35*ave_gu33(is,ib,jb)-1.0
!         write(*,*)"gu1r = ",gu1r
!         write(*,*)"gu2r = ",gu2r
!         write(*,*)"gu3r = ",gu3r
!         write(*,*)"gu4r = ",gu4r
!         write(*,*)"gu5r = ",gu5r
!         write(*,*)"gu6r = ",gu6r
!         write(*,*)"gu7r = ",gu7r
!         write(*,*)"gu8r = ",gu8r
!         write(*,*)"gu9r = ",gu9r
!         write(*,*)"gu10r = ",gu10r
!         write(*,*)"gu1i = ",gu1i
!         write(*,*)"gu2i = ",gu2i
!         write(*,*)"gu3i = ",gu3i
!         write(*,*)"gu4i = ",gu4i
!         write(*,*)"gu5i = ",gu5i
!         write(*,*)"gu6i = ",gu6i
!         write(*,*)"gu7i = ",gu7i
!         write(*,*)"gu8i = ",gu8i
!         write(*,*)"gu9i = ",gu9i
!         write(*,*)"gu10i = ",gu10i
!
          kpar_gnp0 = k_par0*ave_kpargnp0(is,ib,jb)
          kpar_gp1p0 = k_par0*ave_kpargp1p0(is,ib,jb)
          kpar_gp3p0 = k_par0*ave_kpargp3p0(is,ib,jb)
          if(use_bper_in)then
            kpar_gp1b0 = k_par0*ave_kpargp1p0(is,ib,jb)
            kpar_gr11b0 = k_par0*ave_kpargr11b0(is,ib,jb)
            kpar_gr13b0 = k_par0*ave_kpargr13b0(is,ib,jb)
            if(vpar_model_in.eq.0)then
              kpar_gnbp = k_par0*ave_kpargnbp(is,ib,jb)
              kpar_gp1bp = k_par0*ave_kpargp1bp(is,ib,jb)
              kpar_gp3bp = k_par0*ave_kpargp3bp(is,ib,jb)
              kpar_gr11bp = k_par0*ave_kpargr11bp(is,ib,jb)
              kpar_gr13bp = k_par0*ave_kpargr13bp(is,ib,jb)
            endif
          endif
          kpar_gu1 = ave_kpargu1(is,ib,jb)
          kpar_gu3 = ave_kpargu3(is,ib,jb)
          kpar_gb1 = ft2*b1*ave_kpar_eff(is,ib,jb)
          kpar_gb3 = ft2*b3*ave_kpar_eff(is,ib,jb)
          kpar_gb33 = b33*ave_kpar_eff(is,ib,jb)
          kpar_gb1gt1 = ft2*b1*ave_kpargt1(is,ib,jb)
          kpar_gb3gt3 = ft2*b3*ave_kpargt3(is,ib,jb)
          kpar_gb33gt1 = b33*ave_kpargt1(is,ib,jb)
          modkpar_gd1 = ft*d1*ave_modkpar_eff(is,ib,jb)
          modkpar_gd3 = ft*d3*ave_modkpar_eff(is,ib,jb)
          modkpar_gd33 = (d33/ft)*ave_modkpar_eff(is,ib,jb)
          modkpar_gd1gu1 = ft*d1*ave_modkpargu1(is,ib,jb)
          modkpar_gd3gu3 = ft*d3*ave_modkpargu3(is,ib,jb)
          modkpar_gd33gu1 = (d33/ft)*ave_modkpargu1(is,ib,jb)
!          grad_gu1 = ave_gradgu1(is,ib,jb)
!          grad_gu3 = ave_gradgu3(is,ib,jb)
          grad_gu1 = 0.0
          grad_gu3 = 0.0
          dgr13 = -b3*(ft2*ave_kpargt3(is,ib,jb)-ave_kpargt1(is,ib,jb)/3.0 &
              -ft2*ave_kpargu3(is,ib,jb) + ave_kpargu1(is,ib,jb)/3.0)
!
         if(nbasis.eq.1.or.Linsker.eq.0.0)then
           gradgp1=0.0
           gradgr11=0.0
           gradgr13=0.0
           gradgp1p1=0.0
           gradgr11p1=0.0
           gradgr13p1=0.0
         else
          gradgp1 = Linsker*ave_gradgp1p0(is,ib,jb)
          gradgr11 = Linsker*ave_gradgr11p0(is,ib,jb)
          gradgr13 = Linsker*ave_gradgr13p0(is,ib,jb)
          gradgp1p1 = Linsker*ave_gradgp1p1(is,ib,jb)
          gradgr11p1 = Linsker*ave_gradgr11p1(is,ib,jb)
          gradgr13p1 = Linsker*ave_gradgr13p1(is,ib,jb) 
         endif    
!
        endif  ! nroot>6
!        write(*,*)is,js,"dhr13 = ",dhr13
!        write(*,*)is,js,"dgr13 = ",dgr13
!
        w_d1 = -ghat_in*w_d0
        k_par1 = k_par0
        if(js.ne.is)then
          w_d1 = 0.0
          k_par1 = 0.0
        endif
        modw_d1 = ABS(w_d1)
        modk_par1 = ABS(k_par1)
        modk_par0 = ABS(k_par0)
        w_dh= w_d1*ave_wdh(ib,jb)
        w_dg= w_d1*ave_wdg(ib,jb)
        k_par = k_par1*ave_kpar_eff(is,ib,jb)
        k_par_psi = k_par0*ave_kpar_eff(is,ib,jb)
        gradB1=k_par1*gradB_factor_in
!        write(*,*)is,js,ib,jb,"k_par=",k_par
!        write(*,*)"gradB=",ib,jb,gradB1
        if(nbasis.eq.1.or.gradB1.eq.0.0)then
         gradB=0.0
         gradBhp1=0.0
         gradBhp3=0.0
         gradBhr11=0.0
         gradBhr13=0.0
         gradBhr33=0.0
         gradBhu1=0.0
         gradBhu3=0.0
         gradBhu33=0.0
         if(nroot.gt.6)then
          gradBgp1=0.0
          gradBgp3=0.0
          gradBgr11=0.0
          gradBgr13=0.0
          gradBgr33=0.0
          gradBgu1=0.0
          gradBgu3=0.0
          gradBgu33=0.0
         endif
        else
         gradB = gradB1*ave_gradB(ib,jb)
         gradBhp1=gradB1*ave_gradBhp1(is,ib,jb)
         gradBhp3=gradB1*ave_gradBhp3(is,ib,jb)
         gradBhr11=gradB1*ave_gradBhr11(is,ib,jb)
         gradBhr13=gradB1*ave_gradBhr13(is,ib,jb)
         gradBhr33=gradB1*ave_gradBhr33(is,ib,jb)
         gradBhu1=gradB1*ave_gradBhu1(is,ib,jb)
         gradBhu3=gradB1*ave_gradBhu3(is,ib,jb)
         gradBhu33=gradB1*ave_gradBhu33(is,ib,jb)
         if(nroot.gt.6)then
          gradBgp1=gradB1*ave_gradBgp1(is,ib,jb)
          gradBgp3=gradB1*ave_gradBgp3(is,ib,jb)
          gradBgr11=gradB1*ave_gradBgr11(is,ib,jb)
          gradBgr13=gradB1*ave_gradBgr13(is,ib,jb)
          gradBgr33=gradB1*ave_gradBgr33(is,ib,jb)
          gradBgu1=gradB1*ave_gradBgu1(is,ib,jb)
          gradBgu3=gradB1*ave_gradBgu3(is,ib,jb)
          gradBgu33=gradB1*ave_gradBgu33(is,ib,jb)
         endif
        endif
!
! debug
!        write(*,*)"eigensolver",is,ib,jb
!        write(*,*)kpar_hb1
!        write(*,*)kpar_hb3
!        write(*,*)kpar_hb1ht1
!        write(*,*)kpar_hb3ht3
!        write(*,*)kpar_hb33
!        write(*,*)kpar_hb33ht1
!        write(*,*)modkpar_hd1
!        write(*,*)modkpar_hd3
!        write(*,*)modkpar_hd1hu1
!        write(*,*)modkpar_hd3hu3
!        write(*,*)modkpar_hd33
!        write(*,*)modkpar_hd33hu1
!        write(*,*)kpar_gb1
!        write(*,*)kpar_gb3
!        write(*,*)kpar_gb1gt1
!        write(*,*)kpar_gb3gt3
!        write(*,*)kpar_gb33
!        write(*,*)kpar_gb33gt1
!        write(*,*)modkpar_gd1
!        write(*,*)modkpar_gd3
!        write(*,*)modkpar_gd1gu1
!        write(*,*)modkpar_gd3gu3
!        write(*,*)modkpar_gd33
!        write(*,*)modkpar_gd33gu1
!         write(*,*)kpar_hp3
!         write(*,*)kpar_hr11
!         write(*,*)kpar_hr13
!         write(*,*)kpar_hu1
!         write(*,*)kpar_hu3
!         write(*,*)ave_kparht1(is,ib,jb)
!         write(*,*)ave_kparht3(is,ib,jb)
!
! 
           M_i = zs(is)*vs(is)/taus(is)
           J_j = zs(js)*as(js)*vs(js)
           E_i = zs(is)/taus(is)
           N_j = zs(js)*as(js)
!
!        xnu_therm=0.0
!        write(*,*)is,js,ib,jb,xnuei,d_ab
!        write(*,*)is,js,ib,jb,gradB,mod_kpar
!        write(*,*)is,js,ib,jb,w_s,w_cd,w_dh,modw_dh
!         write(*,*)is,js,ib,jb,d_ee
!
! matrix in order n, u, p1, p3, q1,q3
!
! n_u equ #1
!
      ia0 = (is-ns0)*nroot*nbasis
      ja0 = (js-ns0)*nroot*nbasis
!      write(*,*)"ia0 = ",ia0,"ja0 = ",ja0
!
      ia = ib + ia0
!
!  untrapped terms
!
      ja = jb + ja0
!
      phi_A = N_j*xi*w_s*(rlns(is)*hn + rlts(is)*1.5*(hp3-hn)) 
      if(vpar_model_in.eq.0)then
        phi_A = phi_A + N_j*E_i*kpar_hnp0*vpar_s(is)
      endif
      phi_B = -hn*E_i*N_j
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
        sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*h10n + rlts(is)*1.5*(h10p3-h10n)))
        sig_B = betae_sig*h10n*as(js)*taus(js)*zs(is)*zs(is) &
         /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*hp1b0
       if(vpar_model_in.eq.0)then
         psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdhp1b0
         psi_B = betae_psi*M_i*J_j*vpar_s(is)*hp1b0/vs(is)
         phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*hnbp + rlts(is)*1.5*(hp3bp-hnbp)) &
          + E_i*kpar_hnbp*vpar_s(is)) 
         phi_BU = -betae_psi*U0*E_i*J_j*vpar_s(is)*hp1bp
         psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdhp1bp+w_s*vpar_shear_s(is)*hp1bp)
         psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*hp1bp/vs(is)
       endif
      endif
!
      amat(ia,ja) = phi_A + psi_AN  
      bmat(ia,ja) = d_ab + phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU -k_par*vs(is) + am*gradB*vs(is)
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A -0.5*xi*w_dh*taus(is)/zs(is)
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A -1.5*xi*w_dh*taus(is)/zs(is)
      bmat(ia,ja) = 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0

      if(nroot.gt.6)then
!
!  n_u ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A)
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A)
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
!  n_u trapped particle terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A
      bmat(ia,ja) = 1.5*sig_B
!
      endif !  nroot.gt.6
!
! u_par_u equ #2
!
      ia = nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*vpar_shear_s(is)*hp1/vs(is)
      phi_B = 0.0 
      if(vpar_model_in.eq.0)then
        phi_A = phi_A  +N_j*xi*w_cd*wdhp1p0*vpar_s(is)/vs(is)  &
          + d_1*(nuei_u_u_1+nuei_u_q3_1*5.0/3.0)*hp1*E_i*N_j*vpar_s(is)/vs(is)
        phi_B = -E_i*N_j*hp1*vpar_s(is)/vs(is)
      endif
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*hp1b0+1.5*rlts(is)*(hr13b0-hp1b0))  
       psi_B = betae_psi*M_i*J_j*hp1b0
       psi_A = psi_A - damp_psi*psi_B
       if(vpar_model_in.eq.0)then
        psi_A = psi_A  -betae_psi*J_j*M_i*kpar_hp1b0*vpar_s(is)
        phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear_s(is)*hp1bp +w_cd*wdhp1bp*vpar_s(is))/vs(is) &
          + d_1*(nuei_u_u_1+nuei_u_q3_1*5.0/3.0)*hp1*E_i*betae_psi*U0*J_j*vpar_s(is)/vs(is)
        phi_BU = -betae_psi*U0*E_i*J_j*hp1bp*vpar_s(is)/vs(is)
        psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*hp1bp+1.5*rlts(is)*(hr13bp-hp1bp))  &
          + M_i*kpar_hp1bp*vpar_s(is))
        psi_BN = -betae_psi*U0*M_i*N_j*hp1bp
       endif
      endif
!
! u_par_u  untrapped terms
!
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) =  psi_A +phi_AU  &
       -d_ee*nuei_u_u_1  &
       +xnuei*(d_ab*xnu_u_u_1 - d_ij*hu3*xnu_u_q3_1)
!      +resist(is,js)
      bmat(ia,ja) = d_ab + psi_B +phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) =  -(k_par - k_par1*gradhp1p1)*vs(is) &
       + (am+bm*0.5)*gradB*vs(is) 
      bmat(ia,ja) = 0.0
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) =  - bm*1.5*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = -0.5*xi*w_dh*taus(is)/zs(is)
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = -1.5*xi*w_dh*taus(is)/zs(is) &
       -d_ee*nuei_u_q3_1 &
       +xnuei*d_ab*xnu_u_q3_1
      bmat(ia,ja) = 0.0

      if(nroot.gt.6)then
!
!  u_par_u ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 0.0 
      bmat(ia,ja) = 0.0
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! u_par_u trapped particle terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      endif  !  nroot.gt.6
!
! p_par_u equ #3
!
      ia = 2*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*hp1 + rlts(is)*1.5*(hr13-hp1))  
      if(vpar_model_in.eq.0)then
        phi_A = phi_A +N_j*E_i*kpar_hp1p0*vpar_s(is)
      endif
      phi_B = -hp1*E_i*N_j
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
        sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*h10p1 + rlts(is)*1.5*(h10r13-h10p1)))
        sig_B = betae_sig*h10p1*as(js)*taus(js)*zs(is)*zs(is) &
        /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*hr11b0
       if(vpar_model_in.eq.0)then
        psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdhr11b0
        psi_B = betae_psi*M_i*J_j*vpar_s(is)*hr11b0/vs(is)
        phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*hp1bp + rlts(is)*1.5*(hr13bp-hp1bp)) &
         +E_i*kpar_hp1bp*vpar_s(is))
        phi_BU = -betae_psi*U0*hp1bp*E_i*J_j
        psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdhr11bp+w_s*vpar_shear_s(is)*hr11bp)
        psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*hr11bp/vs(is)
       endif
      endif
!
!  p_par_u  untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN  &
       +2.0*taus(is)*(modw_d1*hv1rht1/ABS(zs(is)) +w_d1*xi*hv1iht1/zs(is)) &
       +2.0*taus(is)*(modw_d1*hv2rht3/ABS(zs(is)) +w_d1*xi*hv2iht3/zs(is)) &
      -xnuei*d_ij*xnu_p1_1*(ht1-ht3)
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = k_par1*grad_hu1*vs(is) + psi_A + phi_AU
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) =  -0.5*sig_A   & 
       -xi*w_d1*(taus(is)/zs(is))*(0.5*wdhu1+1.5*wdhu3) &
       -2.0*taus(is)*(modw_d1*hv1r/ABS(zs(is)) +w_d1*xi*hv1i/zs(is)) &
       -d_ee*nuei_p1_p1_1  &
       +xnuei*d_ab*xnu_p1_1
      bmat(ia,ja) = d_ab -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A -2.0*taus(is)* &
           (modw_d1*hv2r/ABS(zs(is)) +w_d1*xi*hv2i/zs(is)) &
       -d_ee*nuei_p1_p3_1  &
       -xnuei*d_ab*xnu_p1_1
      bmat(ia,ja) = 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = -k_par*vs(is) + (am +bm)*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = -bm*3.0*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
      if(nroot.gt.6)then
!
!   p_par_u ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A)
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A)
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
!  p_par_u trapped particle terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A
      bmat(ia,ja) = 1.5*sig_B
!
      endif   ! nroot >6
!
! p_tot_u equ #4
!
      ia = 3*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*hp3+rlts(is)*1.5*(hr33-hp3))  
      if(vpar_model_in.eq.0)then
        phi_A = phi_A + N_j*E_i*kpar_hp3p0*vpar_s(is)
      endif
      phi_B = -hp3*E_i*N_j
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
         sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*h10p3 + rlts(is)*1.5*(h10r33-h10p3)))
         sig_B = betae_sig*h10p3*as(js)*taus(js)*zs(is)*zs(is) &
         /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*hr13b0
       if(vpar_model_in.eq.0)then
         psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdhr13b0
         psi_B = betae_psi*M_i*J_j*vpar_s(is)*hr13b0/vs(is)
         phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*hp3bp+rlts(is)*1.5*(hr33bp-hp3bp))    &
          + E_i*kpar_hp3bp*vpar_s(is))
         phi_BU = -betae_psi*U0*hp3bp*E_i*J_j
         psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdhr13bp +w_s*vpar_shear_s(is)*hr13bp)
         psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*hr13bp/vs(is)
        endif
      endif
!
!   p_tot_u untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN  &
       +2.0*taus(is)*(modw_d1*hv3rht1/ABS(zs(is)) +w_d1*xi*hv3iht1/zs(is)) &
       +2.0*taus(is)*(modw_d1*hv4rht3/ABS(zs(is)) +w_d1*xi*hv4iht3/zs(is))
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = k_par1*grad_hu3*vs(is) + psi_A + phi_AU
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A  &
       -xi*w_d1*(taus(is)/zs(is))*0.5*wdhu3 &
       -2.0*taus(is)*(modw_d1*hv3r/ABS(zs(is)) +w_d1*xi*hv3i/zs(is))
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) =  1.5*sig_A  &
       -xi*w_d1*(taus(is)/zs(is))*1.5*wdhu33 &
       -2.0*taus(is)*(modw_d1*hv4r/ABS(zs(is)) +w_d1*xi*hv4i/zs(is))
      bmat(ia,ja) = d_ab + 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = -k_par*vs(is) + am*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
      if(nroot.gt.6)then
!
!   p_tot_u ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A)
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A)
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
!   p_tot_u trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A
      bmat(ia,ja) = 1.5*sig_B
!
      endif  ! nroot>6
!
! q_par_u equ #5
!
      ia = 4*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*hr11*vpar_shear_s(is)/vs(is) 
      phi_B = 0.0
      if(vpar_model_in.eq.0)then
        phi_A = phi_A +N_j*xi*w_cd*wdhr11p0*vpar_s(is)/vs(is) &
         + d_1*(nuei_q1_u_1 + (5.0/3.0)*nuei_q1_q3_1+3.0*nuei_q1_q1_1)*hr11*E_i*N_j*vpar_s(is)/vs(is)
        phi_B = -E_i*N_j*hr11*vpar_s(is)/vs(is)
      endif
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*hr11b0+1.5*rlts(is)*(hw113b0-hr11b0)) 
       psi_B =betae_psi*M_i*J_j*hr11b0
       psi_A = psi_A - damp_psi*psi_B
       if(vpar_model_in.eq.0)then
        psi_A = psi_A  -betae_psi*J_j*M_i*kpar_hr11b0*vpar_s(is)
        phi_AU = betae_psi*U0*J_j*xi*(w_s*hr11bp*vpar_shear_s(is) +w_cd*wdhr11bp*vpar_s(is))/vs(is) &
         + d_1*(nuei_q1_u_1 + (5.0/3.0)*nuei_q1_q3_1+3.0*nuei_q1_q1_1)*hr11bp*E_i*betae_psi*U0*J_j*vpar_s(is)/vs(is)
        phi_BU = -betae_psi*U0*E_i*J_j*hr11bp*vpar_s(is)/vs(is)
        psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*hr11bp+1.5*rlts(is)*(hw113bp-hr11bp)) &
         +M_i*kpar_hr11bp*vpar_s(is))
        psi_BN =-betae_psi*U0*M_i*N_j*hr11bp
       endif
      endif
!
!  q_par_u untrapped terms
!
      ja =          jb + ja0
      amat(ia,ja) = k_par1*kpar_hb1ht1*vs(is) + phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja =   nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU +modk_par1*modkpar_hd1hu1*vs(is) &
       - taus(is)*(modw_d1*hv5r/ABS(zs(is)) + w_d1*xi*hv5i/zs(is)) &
      -d_ee*nuei_q1_u_1  &
      +xnuei*(d_ab*xnu_q1_u_1 - d_ij*hu1*xnu_q1_q1_1 &
       - d_ij*hu3*xnu_q1_q3_1)
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -k_par1*(kpar_hu1 -grad_hu1 - gradhr11p1)*vs(is) &
       - k_par1*kpar_hb1*vs(is) + (am+bm*1.5)*gradBhu1*vs(is)      &
       -d_11*k_par*vs(is)*c06
      bmat(ia,ja) = 0.0
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = -bm*4.5*gradBhu3*vs(is)    &
        +d_11*k_par*vs(is)*c06 - d_11*k_par*vs(is)*c08
      bmat(ia,ja) = 0.0
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) =   &
        -modk_par1*modkpar_hd1*vs(is) &
        - taus(is)*(modw_d1*hv6r/ABS(zs(is)) + w_d1*xi*hv6i/zs(is)) &
       -d_ee*nuei_q1_q1_1  &
       +xnuei*d_ab*xnu_q1_q1_1
      bmat(ia,ja) = d_ab
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = - taus(is)*(modw_d1*hv7r/ABS(zs(is)) + w_d1*xi*hv7i/zs(is)) &
       -d_ee*nuei_q1_q3_1   &
       +xnuei*d_ab*xnu_q1_q3_1
      bmat(ia,ja) = 0.0
!
      if(nroot.gt.6)then
!
!  q_par_u ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
!  q_par_u trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      endif   ! nroot >6
!
! q_tot_u equ #6
!
      ia = 5*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*vpar_shear_s(is)*hr13/vs(is) 
      phi_B = 0.0
      if(vpar_model_in.eq.0)then
        phi_A = phi_A +N_j*xi*w_cd*wdhr13p0*vpar_s(is)/vs(is) &
         + d_1*(nuei_q3_u_1 + (5.0/3.0)*nuei_q3_q3_1)*hr13*E_i*N_j*vpar_s(is)/vs(is)
        phi_B = -E_i*N_j*hr13*vpar_s(is)/vs(is)
      endif
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*hr13b0+1.5*rlts(is)*(hw133b0-hr13b0)) 
       psi_B = hr13b0*betae_psi*M_i*J_j
       psi_A = psi_A - damp_psi*psi_B
       if(vpar_model_in.eq.0)then
        psi_A = psi_A -betae_psi*J_j*M_i*kpar_hr13b0*vpar_s(is)
        phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear_s(is)*hr13bp +w_cd*wdhr13bp*vpar_s(is))/vs(is) &
        + d_1*(nuei_q3_u_1 + (5.0/3.0)*nuei_q3_q3_1)*hr13bp*E_i*betae_psi*U0*J_j*vpar_s(is)/vs(is)
        phi_BU = -betae_psi*U0*E_i*J_j*hr13bp*vpar_s(is)/vs(is)
        psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*hr13bp+1.5*rlts(is)*(hw133bp-hr13bp)) &
         + M_i*kpar_hr13bp*vpar_s(is))
        psi_BN = -betae_psi*U0*hr13bp*M_i*N_j
       endif
      endif
!
!  q_tot_u untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN  &
       +k_par1*(kpar_hb3ht3 -dhr13+ kpar_hb33ht1)*vs(is)
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) =  psi_A  + phi_AU + &
       modk_par1*(modkpar_hd3hu3 + modkpar_hd33hu1)*vs(is) &
        - taus(is)*(modw_d1*hv8r/ABS(zs(is)) + w_d1*xi*hv8i/zs(is)) &
       -d_ee*nuei_q3_u_1  &
       +xnuei*(d_ab*xnu_q3_u_1 - d_ij*hu3*xnu_q3_q3_1)
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = - k_par1*(kpar_hu3 -grad_hu3 - gradhr13p1)*vs(is) & 
       - k_par1*kpar_hb33*vs(is) + (am+bm*0.5)*gradBhu3*vs(is)      &
       -d_11*k_par*vs(is)*c07
      bmat(ia,ja) = 0.0
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = -k_par1*kpar_hb3*vs(is) - bm*1.5*gradBhu33*vs(is) &
       +d_11*k_par*vs(is)*c07 
      bmat(ia,ja) = 0.0
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = -modk_par1*modkpar_hd33*vs(is) &
       - taus(is)*(modw_d1*hv9r/ABS(zs(is)) + w_d1*xi*hv9i/zs(is))
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) =   &
       -modk_par1*modkpar_hd3*vs(is) &
       - taus(is)*(modw_d1*hv10r/ABS(zs(is)) + w_d1*xi*hv10i/zs(is)) &
       -d_ee*nuei_q3_q3_1  &
       +xnuei*d_ab*xnu_q3_q3_1
      bmat(ia,ja) = d_ab
!
      if(nroot.gt.6)then
!
!  q_tot_u ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
!  q_tot_u trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0

      endif    ! nroot>6
!
      if(nroot.gt.6)then
!
! n_g equ #7
!
      ia = 6*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*gn + rlts(is)*1.5*(gp3-gn))  
      if(vpar_model_in.eq.0)then
        phi_A = phi_A + N_j*E_i*kpar_gnp0*vpar_s(is)
      endif
      phi_B = -E_i*N_j*gn
      phi_A = phi_A +xnu_phi_b*xnuei*xnu_n_b*phi_B
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
         sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*g10n + rlts(is)*1.5*(g10p3-g10n)))
         sig_B = betae_sig*g10n*as(js)*taus(js)*zs(is)*zs(is) &
         /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*gp1b0
       if(vpar_model_in.eq.0)then
        psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdgp1b0
        psi_B = betae_psi*M_i*J_j*vpar_s(is)*gp1b0/vs(is)
        phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*gnbp + rlts(is)*1.5*(gp3bp-gnbp))  &
          + E_i*kpar_gnbp*vpar_s(is))
        phi_BU = -betae_psi*U0*E_i*J_j*gnbp
        psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdgp1bp+w_s*vpar_shear_s(is)*gp1bp)
        psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*gp1bp/vs(is)
       endif
      endif
!
!  n_g untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_n_n*bn
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A +d_ee*nuei_n_p1*bp1  &
       -d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A +d_ee*nuei_n_p3*bp3  &
       +d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
      bmat(ia,ja) = 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
!  n_g ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A +psi_AN)  &
       -d_ee*nuei_n_n  &
       +xnuei*d_ab*xnu_n_b 
      bmat(ia,ja) = d_ab - 1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -k_par*vs(is) + am*gradB*vs(is) -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A) &
        -0.5*xi*w_dg*taus(is)/zs(is)   &
        -d_ee*nuei_n_p1
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A) &
        -1.5*xi*w_dg*taus(is)/zs(is)  &
        -d_ee*nuei_n_p3
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! n_g trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A 
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A  
      bmat(ia,ja) = 1.5*sig_B
!
! u_par_g equ #8
!
      ia = 7*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*vpar_shear_s(is)*gp1/vs(is)
      if(vpar_model_in.eq.0)then
        phi_A = phi_A  + N_j*xi*w_cd*wdgp1p0*vpar_s(is)/vs(is)
        phi_B = -E_i*N_j*gp1*vpar_s(is)/vs(is)
      endif
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*gp1b0+1.5*rlts(is)*(gr13b0-gp1b0))
       psi_B =betae_psi*M_i*J_j*gp1b0
       psi_A = psi_A - damp_psi*psi_B
       if(vpar_model_in.eq.0)then
        psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gp1b0*vpar_s(is)
        phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear_s(is)*gp1bp +w_cd*wdgp1bp*vpar_s(is))/vs(is)
        phi_BU = -betae_psi*U0*E_i*J_j*gp1bp*vpar_s(is)/vs(is)
        psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*gp1bp+1.5*rlts(is)*(gr13bp-gp1bp)) &
         + M_i*kpar_gp1bp*vpar_s(is)) 
        psi_BN = -betae_psi*U0*M_i*N_j*gp1bp
       endif
      endif
!
! u_par_g untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU + d_ee*ft3*nuei_u_u 
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = +d_ee*ft5*nuei_u_q1  &
            -d_ee*(1.0 -ft2)*(ft3*nuei_u_u*1.25 +ft3*nuei_u_q3*35.0/12.0 +ft5*nuei_u_q1*6.25) 
                   
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = d_ee*ft3*nuei_u_q3  &
        +d_ee*(1.0-ft2)*(ft3*nuei_u_u*2.25 +ft3*nuei_u_q3*5.25 +ft5*nuei_u_q1*11.25)  
                  
      bmat(ia,ja) = 0.0
!
! u_par_g ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU)  &
       -d_ee*nuei_u_u_t -d_ee*nuei_u_u  &
       +xnuei*(d_ab*xnu_u_u_1 - d_ij*gu3*xnu_u_q3_1) &
       +xnuei*d_ab*xnu_u_b 
!    &   +resist(is,js)
      bmat(ia,ja) = d_ab -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -(k_par - k_par1*gradgp1p1)*vs(is) &
       +(am + bm*0.5)*gradB*vs(is) 
      bmat(ia,ja) = 0.0
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 0.0 &
       -bm*1.5*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = -0.5*xi*w_dg*taus(is)/zs(is)  &
        -d_ee*nuei_u_q1
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = -1.5*xi*w_dg*taus(is)/zs(is) &
       -d_ee*nuei_u_q3_t -d_ee*nuei_u_q3  &
       +xnuei*d_ab*xnu_u_q3_1
      bmat(ia,ja) = 0.0
!
! u_par_g trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! p_par_g equ #9
!
      ia = 8*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*gp1 + rlts(is)*1.5*(gr13-gp1))  
      if(vpar_model_in.eq.0)then
        phi_A = phi_A  +N_j*E_i*kpar_gp1p0*vpar_s(is)
      endif
      phi_B = -E_i*N_j*gp1
      phi_A = phi_A +xnu_phi_b*xnuei*xnu_p1_b*phi_B
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
        sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*g10p1 + rlts(is)*1.5*(g10r13-g10p1)))
        sig_B = betae_sig*g10p1*as(js)*taus(js)*zs(is)*zs(is) &
        /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*gr11b0
       if(vpar_model_in.eq.0)then
        psi_A = psi_A  -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdgr11b0
        psi_B = betae_psi*M_i*J_j*vpar_s(is)*gr11b0/vs(is)
        phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*gp1bp + rlts(is)*1.5*(gr13bp-gp1bp)) &
          + E_i*kpar_gp1bp*vpar_s(is))
        phi_BU = -betae_psi*U0*E_i*J_j*gp1bp
        psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdgr11bp+w_s*vpar_shear_s(is)*gr11bp)
        psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*gr11bp/vs(is)
       endif
      endif
!
! p_par_g untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p1_n*bn
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p1_p1*bp1  &
       -d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A +d_ee*nuei_p1_p3*bp3    &
       +d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
      bmat(ia,ja) = 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! p_par_g ghost terms
!
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN) &
        +2.0*taus(is)*(modw_d1*gu1rgt1/ABS(zs(is)) +w_d1*xi*gu1igt1/zs(is)) &
        +2.0*taus(is)*(modw_d1*gu2rgt3/ABS(zs(is)) +w_d1*xi*gu2igt3/zs(is)) &
        -d_ee*nuei_p1_n  &
        -xnuei*d_ij*xnu_p1_1*(gt1-ft2*gt3) 
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = k_par1*grad_gu1*vs(is) -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A)  &
        -xi*w_d1*(taus(is)/zs(is))*(0.5*wdgu1+1.5*wdgu3) &
        -2.0*taus(is)*(modw_d1*gu1r/ABS(zs(is)) +w_d1*xi*gu1i/zs(is)) & 
       -d_ee*nuei_p1_p1_t -d_ee*nuei_p1_p1  &
       +xnuei*d_ab*(xnu_p1_1 + xnu_p1_b)
      bmat(ia,ja) = d_ab -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A) &
       -2.0*taus(is)*(modw_d1*gu2r/ABS(zs(is)) +w_d1*xi*gu2i/zs(is)) &
       -d_ee*nuei_p1_p3_t -d_ee*nuei_p1_p3  &
       -xnuei*d_ab*xnu_p1_1*ft2
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = -k_par*vs(is) +(am +bm)*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = -bm*3.0*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
! p_par_g trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A 
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A 
      bmat(ia,ja) = 1.5*sig_B
!
! p_tot_g equ #10
!
      ia = 9*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*gp3+rlts(is)*1.5*(gr33-gp3))
      if(vpar_model_in.eq.0)then 
        phi_A = phi_A + N_j*E_i*kpar_gp3p0*vpar_s(is) 
      endif
      phi_B = -gp3*E_i*N_j
      phi_A = phi_A +xnu_phi_b*xnuei*xnu_p3_b*phi_B
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
        sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*g10p3 + rlts(is)*1.5*(g10r33-g10p3)))
        sig_B = betae_sig*g10p3*as(js)*taus(js)*zs(is)*zs(is) &
        /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*gr13b0
       if(vpar_model_in.eq.0)then
        psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdgr13b0
        psi_B = betae_psi*M_i*J_j*vpar_s(is)*gr13b0/vs(is)
        phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*gp3bp+rlts(is)*1.5*(gr33bp-gp3bp))  &
         + E_i*kpar_gp3bp*vpar_s(is)) 
        phi_BU = -betae_psi*U0*gp3bp*E_i*J_j
        psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdgr13bp+w_s*vpar_shear_s(is)*gr13bp)
        psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*gr13bp/vs(is)
       endif
      endif
!
!  p_tot_g untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p3_n*bn
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p3_p1*bp1  &
       -d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A +d_ee*nuei_p3_p3*bp3    &
       +d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
      bmat(ia,ja) = 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
!  p_tot_g ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN) &
       +2.0*taus(is)*(modw_d1*gu3rgt1/ABS(zs(is)) +w_d1*xi*gu3igt1/zs(is)) &
       +2.0*taus(is)*(modw_d1*gu4rgt3/ABS(zs(is)) +w_d1*xi*gu4igt3/zs(is)) &
       -d_ee*nuei_p3_n
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = k_par1*grad_gu3*vs(is) -1.0*(psi_A + phi_AU)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A) &
        -xi*w_d1*(taus(is)/zs(is))*0.5*wdgu3 & 
        -2.0*taus(is)*(modw_d1*gu3r/ABS(zs(is)) +w_d1*xi*gu3i/zs(is))  &
        -d_ee*nuei_p3_p1
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A)  &
       -xi*w_d1*(taus(is)/zs(is))*1.5*wdgu33 &
       -2.0*taus(is)*(modw_d1*gu4r/ABS(zs(is)) +w_d1*xi*gu4i/zs(is)) &
       -d_ee*nuei_p3_p3  &
       +xnuei*d_ab*xnu_p3_b
      bmat(ia,ja) = d_ab + 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = -k_par*vs(is) + am*gradB*vs(is)
      bmat(ia,ja) = 0.0
!
!  p_tot_g trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A 
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A 
      bmat(ia,ja) = 1.5*sig_B
!
!  q_par_g equ #11
!
      ia = 10*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*vpar_shear_s(is)*gr11/vs(is)
      if(vpar_model_in.eq.0)then
        phi_A = phi_A + N_j*xi*w_cd*wdgr11p0*vpar_s(is)/vs(is)
        phi_B = -E_i*N_j*gr11*vpar_s(is)/vs(is)
      endif
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*gr11b0+1.5*rlts(is)*(gw113b0-gr11b0))
       psi_B = gr11b0*betae_psi*M_i*J_j
       psi_A = psi_A - damp_psi*psi_B
       if(vpar_model_in.eq.0)then
        psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gr11b0*vpar_s(is)
        phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear_s(is)*gr11bp +w_cd*wdgr11bp*vpar_s(is))/vs(is)
        phi_BU = -betae_psi*U0*E_i*J_j*gr11bp*vpar_s(is)/vs(is)
        psi_AN = -betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*gr11bp+1.5*rlts(is)*(gw113bp-gr11bp)) &
         +M_i*kpar_gr11bp*vpar_s(is))
        psi_BN = -gr11bp*betae_psi*U0*M_i*N_j
       endif
      endif
!
! q_par_g untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU +d_ee*ft3*nuei_q1_u
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = +d_ee*ft5*nuei_q1_q1 &
        -d_ee*(1.0 -ft2)*(ft3*nuei_q1_u*1.25 +ft3*nuei_q1_q3*35.0/12.0+ft5*nuei_q1_q1*6.25)
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = d_ee*ft3*nuei_q1_q3  &
       +d_ee*(1.0 -ft2)*(ft3*nuei_q1_u*2.25 +ft3*nuei_q1_q3*5.25 +ft5*nuei_q1_q1*11.25)
      bmat(ia,ja) = 0.0
!
! q_par_g ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN) + k_par1*kpar_gb1gt1*vs(is) 
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU) +modk_par1*modkpar_gd1gu1*vs(is)  &
        - taus(is)*(modw_d1*gu5r/ABS(zs(is)) + w_d1*xi*gu5i/zs(is)) &
        -d_ee*nuei_q1_u_t -d_ee*nuei_q1_u &
        +xnuei*(d_ab*xnu_q1_u_1*ft2 - d_ij*gu1*xnu_q1_q1_1 &
        - d_ij*gu3*xnu_q1_q3_1*ft2)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) =  -k_par1*(kpar_gu1 -grad_gu1 - gradgr11p1)*vs(is) &
        - k_par1*kpar_gb1*vs(is) + (am+bm*1.5)*gradBgu1*vs(is)
      bmat(ia,ja) = 0.0
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = -bm*4.5*gradBgu3*vs(is) 
      bmat(ia,ja) = 0.0
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) =    &
        -modk_par1*modkpar_gd1*vs(is) &
        - taus(is)*(modw_d1*gu6r/ABS(zs(is)) + w_d1*xi*gu6i/zs(is)) &
        -d_ee*nuei_q1_q1_t -d_ee*nuei_q1_q1  &
        +xnuei*d_ab*(xnu_q1_q1_1 + xnu_q1_b)
      bmat(ia,ja) = d_ab
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = - taus(is)*(modw_d1*gu7r/ABS(zs(is)) + w_d1*xi*gu7i/zs(is)) &
        -d_ee*nuei_q1_q3_t -d_ee*nuei_q1_q3  &
        +xnuei*d_ab*xnu_q1_q3_1*ft2
      bmat(ia,ja) = 0.0
!
! q_par_g trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! q_tot_g equ #12
!
      ia = 11*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*vpar_shear_s(is)*gr13/vs(is)
      if(vpar_model_in.eq.0)then
        phi_A = phi_A + N_j*xi*w_cd*wdgr13p0*vpar_s(is)/vs(is)
        phi_B = -E_i*N_j*gr13*vpar_s(is)/vs(is)
      endif
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*gr13b0+1.5*rlts(is)*(gw133b0-gr13b0))
       psi_B = gr13b0*betae_psi*M_i*J_j
       psi_A = psi_A - damp_psi*psi_B
       if(vpar_model_in.eq.0)then
        psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gr13b0*vpar_s(is)
        phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear_s(is)*gr13bp +w_cd*wdgr13bp*vpar_s(is))/vs(is)
        phi_BU = -betae_psi*U0*E_i*J_j*gr13bp*vpar_s(is)/vs(is)
        psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*gr13bp+1.5*rlts(is)*(gw133bp-gr13bp)) &
         +M_i*kpar_gr13bp*vpar_s(is)) 
        psi_BN = -gr13bp*betae_psi*U0*M_i*N_j
       endif
      endif
!
! q_tot_g untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = psi_A + phi_AU + d_ee*ft3*nuei_q3_u
      bmat(ia,ja) = psi_B + phi_BU
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = d_ee*ft5*nuei_q3_q1   &
       -d_ee*(1.0 -ft2)*(ft3*nuei_q3_u*1.25 +ft3*nuei_q3_q3*35.0/12.0 +ft5*nuei_q3_q1*6.25)
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = d_ee*ft3*nuei_q3_q3  &
       +d_ee*(1.0 -ft2)*(ft3*nuei_q3_u*2.25 + ft3*nuei_q3_q3*5.25 + ft5*nuei_q3_q1*11.25)
      bmat(ia,ja) = 0.0
!
! q_tot_g ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN) + k_par1*(kpar_gb3gt3 -dgr13 + kpar_gb33gt1)*vs(is)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(psi_A + phi_AU) +  &
       modk_par1*(modkpar_gd3gu3 + modkpar_gd33gu1)*vs(is) &
        - taus(is)*(modw_d1*gu8r/ABS(zs(is)) + w_d1*xi*gu8i/zs(is)) &
       -d_ee*nuei_q3_u_t -d_ee*nuei_q3_u  &
       +xnuei*(d_ab*xnu_q3_u_1 - d_ij*gu3*xnu_q3_q3_1)
      bmat(ia,ja) = -1.0*(psi_B + phi_BU)
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = - k_par1*(kpar_gu3 -grad_gu3 - gradgr13p1)*vs(is) &
       - k_par1*kpar_gb33*vs(is) + (am+bm*0.5)*gradBgu3*vs(is)
      bmat(ia,ja) = 0.0
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = &
        - k_par1*kpar_gb3*vs(is) - bm*1.5*gradBgu33*vs(is) 
      bmat(ia,ja) = 0.0
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = -modk_par1*modkpar_gd33*vs(is) &
       - taus(is)*(modw_d1*gu9r/ABS(zs(is)) + w_d1*xi*gu9i/zs(is))  &
       -d_ee*nuei_q3_q1
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) =   &
        -modk_par1*modkpar_gd3*vs(is) &
        - taus(is)*(modw_d1*gu10r/ABS(zs(is)) +w_d1* xi*gu10i/zs(is)) &
        -d_ee*nuei_q3_q3_t -d_ee*nuei_q3_q3  &
        +xnuei*d_ab*(xnu_q3_q3_1 + xnu_q3_b)
      bmat(ia,ja) = d_ab
!
! q_tot_g trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! n_t equ #13
!
      ia = 12*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*gn + rlts(is)*1.5*(gp3-gn))       
      phi_B = -gn*E_i*N_j
      phi_A = phi_A +xnu_phi_b*xnuei*xnu_n_b*phi_B
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
        sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*g10n + rlts(is)*1.5*(g10p3-g10n)))
        sig_B = betae_sig*g10n*as(js)*taus(js)*zs(is)*zs(is) &
        /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*gp1b0
       if(vpar_model_in.eq.0)then
        psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdgp1b0
        psi_B = betae_psi*M_i*J_j*vpar_s(is)*gp1b0/vs(is)
        phi_AU = betae_psi*U0*J_j*xi*w_s*(rlns(is)*gnbp + rlts(is)*1.5*(gp3bp-gnbp))
        phi_BU = -betae_psi*U0*gnbp*E_i*J_j
        psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdgp1bp+w_s*vpar_shear_s(is)*gp1bp)
        psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*gp1bp/vs(is)
       endif
      endif
!
! n_t untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_n_n*bn
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A +d_ee*nuei_n_p1*bp1  &
       -d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A +d_ee*nuei_n_p3*bp3    &
       +d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
      bmat(ia,ja) = 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! n_t ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN)
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A) 
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A) 
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! n_t trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN   &
        -d_ee*nuei_n_n  &
        +xnuei*d_ab*xnu_n_b 
      bmat(ia,ja) = d_ab + phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A -0.5*xi*w_dg*taus(is)/zs(is)  &
        -d_ee*nuei_n_p1
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A -1.5*xi*w_dg*taus(is)/zs(is)  &
       -d_ee*nuei_n_p3
      bmat(ia,ja) = 1.5*sig_B
!
! p_par_t equ #14
!
      ia = 13*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*gp1 + rlts(is)*1.5*(gr13-gp1))     
      phi_B = -gp1*E_i*N_j
      phi_A = phi_A +xnu_phi_b*xnuei*xnu_p1_b*phi_B
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
        sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*g10p1 + rlts(is)*1.5*(g10r13-g10p1)))
        sig_B = betae_sig*g10p1*as(js)*taus(js)*zs(is)*zs(is) &
        /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*gr11b0
       if(vpar_model_in.eq.0)then
        psi_A = psi_A  -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdgr11b0
        psi_B = betae_psi*M_i*J_j*vpar_s(is)*gr11b0/vs(is)
        phi_AU = betae_psi*U0*J_j*xi*w_s*(rlns(is)*gp1bp + rlts(is)*1.5*(gr13bp-gp1bp))
        phi_BU = -betae_psi*U0*gp1bp*E_i*J_j
        psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdgr11bp+ w_s*vpar_shear_s(is)*gr11bp)
        psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*gr11bp/vs(is) 
       endif
      endif
!
! p_par_t  untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p1_n*bn
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p1_p1*bp1  &
       -d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A
      bmat(ia,ja) = 1.5*sig_B +d_ee*nuei_p1_p3*bp3   &
       +d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! p_par_t  ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN) 
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A)
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A)
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! p_par_t  trapped terms
!
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A +psi_AN &
        +2.0*taus(is)*(modw_d1*gu1rgt1/ABS(zs(is)) +w_d1*xi*gu1igt1/zs(is)) &
        +2.0*taus(is)*(modw_d1*gu2rgt3/ABS(zs(is)) +w_d1*xi*gu2igt3/zs(is)) &
        -d_ee*nuei_p1_n  &
        -xnuei*d_ij*xnu_p1_1*(gt1-ft2*gt3) 
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) =  -0.5*sig_A   &
        -xi*w_d1*(taus(is)/zs(is))*(0.5*wdgu1+1.5*wdgu3) &
        -2.0*taus(is)*(modw_d1*gu1r/ABS(zs(is)) +w_d1*xi*gu1i/zs(is)) &
        -d_ee*nuei_p1_p1_t -d_ee*nuei_p1_p1  &
        +xnuei*d_ab*(xnu_p1_1 + xnu_p1_b)
      bmat(ia,ja) = d_ab -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A  &   
       -2.0*taus(is)*(modw_d1*gu2r/ABS(zs(is)) +w_d1*xi*gu2i/zs(is)) &
       -d_ee*nuei_p1_p3_t -d_ee*nuei_p1_p3  &
       -xnuei*d_ab*xnu_p1_1*ft2
      bmat(ia,ja) = 1.5*sig_B
!
! p_tot_t equ #15
!
      ia = 14*nbasis+ib + ia0
!
      phi_A = N_j*xi*w_s*(rlns(is)*gp3+rlts(is)*1.5*(gr33-gp3)) 
      phi_B = -gp3*E_i*N_j
      phi_A = phi_A +xnu_phi_b*xnuei*xnu_p3_b*phi_B
      sig_A = 0.0
      sig_B = 0.0
      psi_A = 0.0
      psi_B = 0.0
      phi_AU = 0.0
      phi_BU = 0.0
      psi_AN = 0.0
      psi_BN = 0.0
      if(use_bpar_in)then
        sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
        (xi*w_s*(rlns(is)*g10p3 + rlts(is)*1.5*(g10r33-g10p3)))
        sig_B = betae_sig*g10p3*as(js)*taus(js)*zs(is)*zs(is) &
        /(taus(is)*mass(is))
        sig_A = sig_A - damp_sig*sig_B
      endif
      if(use_bper_in)then
       psi_A = -betae_psi*J_j*xi*w_s*vpar_shear_s(is)*gr13b0
       if(vpar_model_in.eq.0)then
        psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar_s(is)*wdgr13b0
        psi_B = betae_psi*M_i*J_j*vpar_s(is)*gr13b0/vs(is)
        phi_AU = betae_psi*U0*J_j*xi*w_s*(rlns(is)*gp3bp+rlts(is)*1.5*(gr33bp-gp3bp)) 
        phi_BU = -betae_psi*U0*gp3bp*E_i*J_j
        psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar_s(is)*wdgr13bp+w_s*vpar_shear_s(is)*gr13bp) 
        psi_BN = -betae_psi*U0*M_i*N_j*vpar_s(is)*gr13bp/vs(is)
       endif
      endif
!
! p_tot_t untrapped terms
!
      ja = jb + ja0
      amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p3_n*bn
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 2*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p3_p1*bp1  &
       -d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 3*nbasis+jb + ja0
      amat(ia,ja) = 1.5*sig_A +d_ee*nuei_p3_p3*bp3    &
       +d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
      bmat(ia,ja) = 1.5*sig_B
!
      ja = 4*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 5*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! p_tot_t ghost terms
!
      ja = 6*nbasis+jb + ja0
      amat(ia,ja) = -1.0*(phi_A + psi_AN) 
      bmat(ia,ja) = -1.0*(phi_B + psi_BN)
!
      ja = 7*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 8*nbasis+jb + ja0
      amat(ia,ja) = -0.5*(-1.0*sig_A) 
      bmat(ia,ja) = -0.5*(-1.0*sig_B)
!
      ja = 9*nbasis+jb + ja0
      amat(ia,ja) = 1.5*(-1.0*sig_A) 
      bmat(ia,ja) = 1.5*(-1.0*sig_B)
!
      ja = 10*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
      ja = 11*nbasis+jb + ja0
      amat(ia,ja) = 0.0
      bmat(ia,ja) = 0.0
!
! p_tot_t trapped terms
!
      ja = 12*nbasis+jb + ja0
      amat(ia,ja) = phi_A + psi_AN &
       +2.0*taus(is)*(modw_d1*gu3rgt1/ABS(zs(is)) +w_d1*xi*gu3igt1/zs(is)) &
       +2.0*taus(is)*(modw_d1*gu4rgt3/ABS(zs(is)) +w_d1*xi*gu4igt3/zs(is)) &
       -d_ee*nuei_p3_n 
      bmat(ia,ja) = phi_B + psi_BN
!
      ja = 13*nbasis+jb + ja0
      amat(ia,ja) = -0.5*sig_A -xi*w_d1*(taus(is)/zs(is))*0.5*wdgu3   &
        -2.0*taus(is)*(modw_d1*gu3r/ABS(zs(is)) +w_d1*xi*gu3i/zs(is)) &
        -d_ee*nuei_p3_p1
      bmat(ia,ja) = -0.5*sig_B
!
      ja = 14*nbasis+jb + ja0
      amat(ia,ja) =    &
       +1.5*sig_A -xi*w_d1*(taus(is)/zs(is))*1.5*wdgu33 &
       -2.0*taus(is)*(modw_d1*gu4r/ABS(zs(is)) +w_d1*xi*gu4i/zs(is)) &
       -d_ee*nuei_p3_p3  &
       +xnuei*d_ab*xnu_p3_b
      bmat(ia,ja) = d_ab + 1.5*sig_B
!
      endif ! nroot .gt. 6
!
      enddo ! end of js loop 
      enddo ! end of ib loop
      enddo ! end of jb loop
      enddo ! end of is loop
!
!
!..find the eigenvalues and eigenvectors
!
!
      lwork=33*iur
      rightvectors="N"
!      if(iflux_in)rightvectors="V"
      do i=1,iur
      do j=1,iur
        at(i,j) = amat(i,j)
        bt(i,j) = bmat(i,j)
!        write(*,*)i,j,at(i,j),bt(i,j)
      enddo
      enddo
!      call system_clock(cpucount1)
      call zggev("N",rightvectors,iur,at,iur,bt,iur,    &
                  alpha,beta,vleft,iur,vright,iur,work,lwork,rwork,info)
!      call system_clock(cpucount2,cpurate)
!      write(*,*)"cputime for zggev =",REAL(cpucount2-cpucount1)/REAL(cpurate)
!      write(*,*)"jmax = ",jmax,alpha(jmax)/beta(jmax)
!      write(*,*)"work(1)",work(1)

      if (info /= 0) then
         call tglf_error(1,"ERROR: ZGGEV failed in tglf_eigensolver.f90")
      endif

!      cputime2=MPI_WTIME()
      do j1=1,iur
        beta2=REAL(CONJG(beta(j1))*beta(j1))
        if(beta2.ne.0.0)then  ! zomega = -xi*(frequency+xi*growthrate)
          zomega(j1)=alpha(j1)*CONJG(beta(j1))/beta2
        else
          zomega(j1)=0.0
        endif
!        write(*,*)j1,"zomega=",zomega(j1)
        rr(j1) = REAL(zomega(j1))
        ri(j1) = AIMAG(zomega(j1))
! filter out numerical instabilities that sometimes occur 
! with high mode frequency
        if(filter_in.gt.0.0)then
          if(rr(j1).gt.0.0.and.ABS(ri(j1)).gt.max_freq)then
            rr(j1)=-rr(j1)
!            write(*,*)"filtered mode ",j1,rr(j1),ri(j1)
!            write(*,*)"beta = ",beta(j1),"alpha =",alpha(j1)
!            write(*,*)(vright(j1,j2),j2=1,iur)
          endif
        endif
        do j2=1,iur
          if(iflux_in)then
            vr(j1,j2) = REAL(vright(j1,j2))
            vi(j1,j2) = AIMAG(vright(j1,j2))
          else
            vr(j1,j2) = 0.0
            vi(j1,j2) = 0.0
          endif
        enddo
      enddo
!
      if(ALLOCATED(rwork))DEALLOCATE(rwork)
      if(ALLOCATED(at))DEALLOCATE(at)
      if(ALLOCATED(bt))DEALLOCATE(bt)
      if(ALLOCATED(vleft))DEALLOCATE(vleft)
      if(ALLOCATED(vright))DEALLOCATE(vright)
      if(ALLOCATED(work))DEALLOCATE(work)
      if(ALLOCATED(zomega))DEALLOCATE(zomega)
!
      END SUBROUTINE tglf_eigensolver
!
!
      SUBROUTINE get_v
!*********************************************************
!     compute the toroidal drift closure coefficients for
!     trapped fraction ft=1
!*********************************************************
      USE tglf_dimensions
      USE tglf_closure
!
      IMPLICIT NONE
      REAL :: v(20),vb(20)
!
!    ft = 1.0
      data v / &
       1.826847832031249,-0.861310124999999, &
       0.6575781562499995,2.2825501875,0.4241429609375004, &
       0.915902718750001,0.1407968007812502, &
       0.6294216650390628,-47.27620127343747, &
       -12.47474156250001,1.301349406249996, &
       10.06786334436035,20.47157314453126, &
       7.759854870117187,-22.84431474609372, &
       -3.688952213541671,-1.045452443847656, &
       1.344488212402342,15.16781977294922, &
       4.613708160644531 /
!
      data vb / &
       1.0720427495888059,-1.0352701554051484, &
       -1.1587645266235862,-1.539881808285635, &
       0.7536017884881802,1.3101251808957506, &
       0.12907305691341753,0.41497311224424527, &
       -34.27184170047316,27.21047057820418, &
       3.1070862741557397,-16.217231447363456, &
       28.24681814543993,5.152204009143581, &
       -13.451163138468832,5.097927633723974, &
       -0.08081708900093537,-8.719713606189327, &
       9.245808116510783,13.203793959198736 /
!
        v1_r = v(1)
        v1_i = v(2)
        v2_r = v(3)
        v2_i = v(4)
        v3_r = v(5)
        v3_i = v(6)
        v4_r = v(7)
        v4_i = v(8)
        v5_r = v(9)
        v5_i = v(10)
        v6_r = v(11)
        v6_i = v(12)
        v7_r = v(13)
        v7_i = v(14)
        v8_r = v(15)
        v8_i = v(16)
        v9_r = v(17)
        v9_i = v(18)
        v10_r = v(19)
        v10_i = v(20)
!
        vb1_r = vb(1)
        vb1_i = vb(2)
        vb2_r = vb(3)
        vb2_i = vb(4)
        vb3_r = vb(5)
        vb3_i = vb(6)
        vb4_r = vb(7)
        vb4_i = vb(8)
        vb5_r = vb(9)
        vb5_i = vb(10)
        vb6_r = vb(11)
        vb6_i = vb(12)
        vb7_r = vb(13)
        vb7_i = vb(14)
        vb8_r = vb(15)
        vb8_i = vb(16)
        vb9_r = vb(17)
        vb9_i = vb(18)
        vb10_r = vb(19)
        vb10_i = vb(20)
!
!  debug
!      write(*,*)" v1_r = ",v1_r," v1_i = ",v1_i
!      write(*,*)" v2_r = ",v2_r," v2_i = ",v2_i
!      write(*,*)" v3_r = ",v3_r," v3_i = ",v3_i
!      write(*,*)" v4_r = ",v4_r," v4_i = ",v4_i
!      write(*,*)" v5_r = ",v5_r," v5_i = ",v5_i
!      write(*,*)" v6_r = ",v6_r," v6_i = ",v6_i
!      write(*,*)" v7_r = ",v7_r," v7_i = ",v7_i
!      write(*,*)" v8_r = ",v8_r," v8_i = ",v8_i
!      write(*,*)" v9_r = ",v9_r," v9_i = ",v9_i
!      write(*,*)" v10_r = ",v10_r," v10_i = ",v10_i
!      write(*,*)"debug get_v"
!      do i=1,20
!       write(*,*)i,v(i)
!      enddo
!
      END SUBROUTINE get_v
!
      SUBROUTINE get_u(fx)
!*********************************************************
!     compute the toroidal drift closure coefficients for
!     trapped fraction fx
!*********************************************************
      USE tglf_closure
!
      IMPLICIT NONE
!
      INTEGER :: j1,j2,i
      REAL :: fx,df,fm(21),v(20),vb(20)
      REAL :: vm(21,20),vbm(21,20)
!
      DATA fm / &
       0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, &
       0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, &
       0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 /
!
!    fm = 0.0  actually fit at ft=0.01
!
      data vm(1,1:20)  / &
       1.575574876228914,-0.6066701129968165, &
       0.6040857655895493,1.987597754698735, &
       0.4509018695196247,0.946277522991777, &
       0.1894230720895441,0.7011692753158356, &
       -48.87702963095776,-21.86758647252671, &
       1.498549195246206,8.79874389216027, &
       20.32547237814752,9.39260392021951, &
       -25.01389867263115,-5.586906841157134, &
       -1.201891743762065,1.685043677575747, &
       15.39835508978011,3.977625916718097 /
!
      data vbm(1,1:20)  / &
       2.100242896319225,0.41667690846001193, &
       -0.23584358880505907,2.0264076479472095, &
       -0.027047472681656068,1.9243266978211786, &
       -0.5910878004088855,0.4585457567254142, &
       -79.60792795010252,-26.402851545261953, &
       -1.0808479486740747,19.235421162998012, &
       25.048572608410915,-34.60901938698052, &
       -41.884208839717886,-7.355935633790796, &
       -2.1073425532059207,5.311681937985481, &
       17.975800034278148,-10.334403370288097 /
!
!    fm = 0.05
!
      data vm(2,1:20)  / &
       1.584786447454745,-0.6005752129838008, &
       0.6028112671584963,1.98765210456097, &
       0.4477119220748415,0.953518302602162, &
       0.1880421376412414,0.6992027675321516, &
       -48.86098519124041,-21.89932338258506, &
       1.491696713468709,8.84295868558821, &
       20.33016280877602,9.37862756702879, &
       -24.98559467865923,-5.516911032433131, &
       -1.198040231611221,1.672188726738941, &
       15.40076145875806,3.979856813017581 /
!
      data vbm(2,1:20)  / &
       2.121457471029516,0.38518780537092984, &
       -0.34081443468939987,1.9803641807448933, &
       0.017338123513881434,1.9057456774658739, &
       -0.6123675735911711,0.40651219567855923, &
       -78.66396042500197,-25.547026168613534, &
       -0.6229671173913928,18.821351431504898, &
       23.85578343658188,-34.573500517308155, &
       -41.730003436393375,-6.8110515127692715, &
       -2.032397881331805,4.828801761804959, &
       18.059324409672904,-9.870490325012534 /
!
!    fm = 0.1
!
       data vm(3,1:20)  / &
       1.595418415489227,-0.5964744511322662, &
       0.6019894105218661,1.989859605060335, &
       0.444599724006794,0.960239982479518, &
       0.1863648539556405,0.696829193093178, &
       -48.85797928040578,-21.86785630615234, &
       1.481310183083416,8.89130515738148, &
       20.34018889309613,9.35582274907792, &
       -24.95167599410477,-5.438059174404275, &
       -1.192188727174912,1.658711694223377, &
       15.39909222122237,3.984030512951928 /
!
       data vbm(3,1:20)  / &
       2.1418890848782683,0.37391202098800846, &
       -0.4733533815130547,1.8897731789776775, &
       0.06991178836242341,1.9049122783441002, &
       -0.6236716217351208,0.33211780692692583, &
       -76.93296863080818,-23.775733986611037, &
       -0.25079191521391675,18.06271730470707, &
       22.719793749125607,-34.65689366769612, &
       -41.28107178077724,-6.038166234724489, &
       -1.9345231030850365,3.8941949691975553, &
       18.187918678454132,-8.812937790189816 /
!
!    fm = 0.15
!
      data vm(4,1:20)  / &
       1.607008019221601,-0.59461395979712, &
       0.6017517435582062,1.994720762935243, &
       0.441374757274443,0.965945095701004, &
       0.1843143567369409,0.6943771736985549, &
       -48.86771017056885,-21.77126756732583, &
       1.468984485137806,8.93991595036156, &
       20.35577378239828,9.33081906985168, &
       -24.91335710987816,-5.352420447248309, &
       -1.184636668413774,1.645751401933153, &
       15.39442522231617,3.989804820318835 /
!
      data vbm(4,1:20)  / &
       2.1569540608973403,0.36115860765518937, &
       -0.6037670682564563,1.789981699244811, &
       0.12116427792447931,1.9107022285542576, &
       -0.6271207860584389,0.2674056416480938, &
       -75.16655459776027,-21.575076212895837, &
       0.11236196918187102,17.329983922005297, &
       21.83414619412646,-35.066597547635695, &
       -40.67100668056874,-5.391219852432574, &
       -1.8397200302749055,3.150643178962453, &
       18.193568634338387,-7.9395835947655815 /
!
!    fm = 0.2
!
      data vm(5,1:20)  / &
       1.619763658028577,-0.5948050641585538, &
       0.6020433583098876,2.001998535933301, &
       0.4378853583252869,0.970311497439482, &
       0.1820248257258578,0.6915676799985609, &
       -48.89215624869315,-21.61359805805468, &
       1.454080163462317,8.99105256180684, &
       20.37554741544917,9.30117158542314, &
       -24.86964094415601,-5.257134877591955, &
       -1.175635707528013,1.632487441471199, &
       15.38653661711697,3.997104219626161 /
!
      data vbm(5,1:20)  / &
       2.148636487619103,0.3458545440796659, &
       -0.702054730530762,1.6846886581127625, &
       0.15288867876905898,1.9089796726777277, &
       -0.6130587509583316,0.21155509624058244, &
       -73.19041343501488,-19.297921478439903, &
       0.3108214915128239,16.528358533147603, &
       21.116195545576957,-35.78224239554675, &
       -40.027132725246105,-5.03499402515285, &
       -1.7483678120930428,2.5164881621105573, &
       18.120837538360504,-7.088913923897836 /
!
!    fm = 0.25
!
      data vm(6,1:20) / &
       1.636176554086759,-0.599485028334832, &
       0.6035244686514099,2.012388177728817, &
       0.4337176032317322,0.973475292138934, &
       0.1791083677917222,0.6884454722124726, &
       -48.94401587488863,-21.38496273382624, &
       1.436516504636102,9.04845368990215, &
       20.40049020229807,9.26173684650638, &
       -24.81680828589113,-5.15721385907243, &
       -1.164160090387029,1.619581402172636, &
       15.37373768022986,4.006870967609708 /
!
      data vbm(6,1:20) / &
       2.1546290496634555,0.33440129956941367, &
       -0.7731880292188994,1.5527084406569154, &
       0.18082634981556542,1.9193109637246648, &
       -0.5806554547601828,0.162235503251981, &
       -71.47501312013225,-17.245684967327854, &
       0.35002420215408037,15.80148999344889, &
       20.84613227338142,-36.35483098353724, &
       -39.23267113476709,-4.870905399506814, &
       -1.6401198987738763,2.016416796562974, &
       17.870648459921636,-6.352073408510565 /
!
!    fm = 0.3
!
      data vm(7,1:20)  / &
       1.646052871314648,-0.605082037178735, &
       0.6055113026131101,2.023326788177401, &
       0.4301553789993933,0.973779598263391, &
       0.1771513984371125,0.6857453498972519, &
       -48.98036849212899,-21.07947502930006, &
       1.423703176051638,9.08208578883894, &
       20.42410557436844,9.22713508991918, &
       -24.76775651809943,-5.073501091069786, &
       -1.154437561548364,1.607649627592847, &
       15.36167636402217,4.017306549075079 / 
!
      data vbm(7,1:20)  / &
       2.1513012555337845,0.348153357177939, &
       -0.8257247675545437,1.401677671547645, &
       0.20337562189294361,1.945892154602142, &
       -0.5433669011675593,0.12537519571250227, &
       -69.49442209055175,-14.492172241451986, &
       0.4132517144676316,15.007054068687692, &
       20.704516614813326,-37.25387168378119, &
       -38.15479808875955,-4.8340458002796804, &
       -1.4911028115450349,1.3514857885811122, &
       17.443287906219293,-5.139217968050646 /
!
!    fm = 0.35
!
      data vm(8,1:20)  / &
       1.654948217986326,-0.6097696412961999, &
       0.6078191158184839,2.03685276356042, &
       0.4263185123879022,0.972563893396645, &
       0.1751589652078732,0.683151509011474, &
       -49.01828106889333,-20.75996617488966, &
       1.409498078231957,9.1227103584037,20.44582926796566, &
       9.18980152123914,-24.71677816313796, &
       -4.994371517341897,-1.145204351464681, &
       1.596920319198234,15.35007386678302, &
       4.028384606743625 /
!
      data vbm(8,1:20)  / &
       2.1513684857990096,0.3702774338504988, &
       -0.8324362851031883,1.2605015031903262, &
       0.21918428871663032,1.9726461682591983, &
       -0.498273178512208,0.11115438196927782, &
       -67.50308119529046,-11.820695139846723, &
       0.39545618609342703,14.432461687744366, &
       20.64194571685896,-38.932850877889294, &
       -37.11556234315118,-5.189528502715698, &
       -1.2394869588902089,0.7234934628378644, &
       16.86156394994619,-3.846720035966003 /
!
!    fm = 0.4
!
      data vm(9,1:20)  / &
       1.669556840339296,-0.6163861614985355, &
       0.6109430394067002,2.052800053034924, &
       0.4205037341885763,0.97268547908153, &
       0.1719945345142707,0.679863400135331, &
       -49.05804491389187,-20.37694331845068, &
       1.39161366819924,9.17777702052686, &
       20.4780011276474,9.14521858065844, &
       -24.63901128377355,-4.884470921605778, &
       -1.134286840623679,1.580131422830657, &
       15.33596597620725,4.047976178953485 /
!
      data vbm(9,1:20)  / &
       2.1605903453628823,0.3700298029682928, &
       -0.8141109030844055,1.1573003975400917, &
       0.22878472493280338,1.9863959681736238, &
       -0.4572676977192338,0.11237931651934113, &
       -66.08196782656529,-9.307047697662053, &
       0.3803809291109026,13.832745596066626, &
       20.816375232092817,-39.81795422793366, &
       -36.08658933485706,-5.311976024040643, &
       -0.9776583270154982,0.24895580859852995, &
       16.181286484837152,-2.7199294308067667 /
!
!    fm = 0.45
!
      data vm(10,1:20)  / &
       1.685489987880985,-0.6245492149090265, &
       0.6147884692362493,2.071982074584791, &
       0.4140025996164738,0.971243789082111, &
       0.1683159113731278,0.6767387703441327, &
       -49.10216013588896,-19.93233925504785, &
       1.375069858958643,9.23927595107621, &
       20.50633752178151,9.08136523137532, &
       -24.55536955622268,-4.770719088343094, &
       -1.121343210981297,1.562118246797276, &
       15.32015210825958,4.069787697394207 /
!
      data vbm(10,1:20)  / &
       2.1730855874909913,0.35178115552541456, &
       -0.7373358116921315,1.0588292749680648, &
       0.23692533123572343,1.9878713414347633, &
       -0.39936043468928534,0.129581224006158, &
       -64.40345283700135,-5.6474803990666915, &
       0.4464564895667831,13.392063022241036, &
       20.833709529476486,-41.82009108881107, &
       -34.78015959024863,-5.509723441949259, &
       -0.5123995424609449,-0.5871599259397859, &
       15.15811380312606,-0.7640251210131144 /
!
!    fm =0.5
!
      data vm(11,1:20)  / &
       1.701119018866823,-0.6358858676431006, &
       0.6190005747094675,2.091782511975988, &
       0.4079031724903291,0.968248270993724, &
       0.1646800838958645,0.6731004782085041, &
       -49.12902757284288,-19.42220295554368, &
       1.361537584770168,9.29948649350977, &
       20.53210290482125,9.0091512533601, &
       -24.45417369583586,-4.656134527700456, &
       -1.110739123412469,1.542219994292789, &
       15.30417165363005,4.095882951668317 /
!
      data vbm(11,1:20)  / &
       2.1798977680160343,0.31842602898883743, &
       -0.6113895619337811,0.9822163960742745, &
       0.2367958335142667,1.9733640324150512, &
       -0.33059638633218674,0.157378137551125, &
       -62.65355852659036,-2.3394699250483972, &
       0.8298447761463912,13.007330572646818, &
       19.936564143039803,-43.83657346835537, &
       -33.841069900509694,-5.779941344277579, &
       -0.04942125216637194,-1.595543277010143, &
       14.488041866787265,1.4254200018898384 /
!
!    fm = 0.55
!
      data vm(12,1:20)  / &
       1.716892973056783,-0.6477529032534867, &
       0.6241498106473077,2.112078264044541, &
       0.4011386584722641,0.965401091991952, &
       0.161465903259113,0.6694004716950336, &
       -49.12825994378126,-18.91389210521596, &
       1.344586958673598,9.37344886344791, &
       20.55336881718631,8.93308048981402, &
       -24.34235102083391,-4.55367679970704, &
       -1.101014925186506,1.524121056743954, &
       15.28852480402593,4.12592484192357 /
!
      data vbm(12,1:20)  / &
       2.173648528496596,0.27785866403912574, &
       -0.46247319359590766,0.9416659079146991, &
       0.22191374311651668,1.94473989212783, &
       -0.25990281944354193,0.18961221391701988, &
       -61.01938317104004,0.08263174360865744, &
       1.2650072807109733,12.793047034813716, &
       18.35779387020247,-45.03275605095713, &
       -33.309164184931554,-6.057448189457001, &
       0.38610353254962426,-2.489147078018906, &
       14.00043303499043,3.636275515024897 /
!
!    fm = 0.6
!
      data vm(13,1:20)  / &
       1.731618494393875,-0.6612253701707162, &
       0.6295055256270565,2.132570937900934, &
       0.3946273078920458,0.961435171907832, &
       0.1583969621180752,0.6653382696813155, &
       -49.09996992204884,-18.30745757310675, &
       1.337051030914882,9.43225306615724, &
       20.56807177474404,8.84312928412704, &
       -24.22206858066527,-4.454290444171457, &
       -1.090751845081662,1.504376120166766, &
       15.27277475505978,4.161361435396874 /
!
      data vbm(13,1:20)  / &
       2.165258153153121,0.22195400023092907, &
       -0.2901337475507507,0.9532539006444243, &
       0.19800467822129636,1.9011326617000779, &
       -0.19012642241663424,0.2356712050549581, &
       -60.57641811358452,1.4496797124330127, &
       2.200251819913351,12.191110932520557, &
       16.043516600571373,-43.986697403333636, &
       -33.41358162751744,-5.51982578292761, &
       1.2375113222744485,-3.4191580742017926, &
       13.035784948780702,5.5431029192454755 /
!
!    fm = 0.65
!
      data vm(14,1:20) / &
       1.744154605621783,-0.67497804789661, &
       0.6354178905303506,2.153262444200677, &
       0.3891788046272644,0.957217432595458, &
       0.1554626054403879,0.6616475171248536, &
       -49.02975299551767,-17.72046710040578, &
       1.329076571485967,9.4954572031658, &
       20.57560678699514,8.76150355762339, &
       -24.10238049004819,-4.362273732623926, &
       -1.080767411765934,1.485144657126245, &
       15.26088199740946,4.198260206745224 /
!
      data vbm(14,1:20) / &
       2.1049984608851573,0.11645015751885165, &
       -0.10395333126146412,1.0078968720804031, &
       0.1646608550696791,1.8043731514533325, &
       -0.10610473521681932,0.2956978733437321, &
       -61.12512767055042,2.754735795597491, &
       3.077275272606009,9.927614765896221, &
       15.162660776107955,-39.203830127748546, &
       -33.29198784379133,-4.424152252648941, &
       2.9324912850104514,-4.546752758247106, &
       10.362309180270689,7.851420565503595 /
!
!    fm = 0.7
!
      data vm(15,1:20) / &
       1.762022616214332,-0.6919511224009725, &
       0.6418710767955955,2.176317811011077, &
       0.3842066305375349,0.953908562270084, &
       0.1513639533899996,0.6580961117405976, &
       -48.92907568841273,-17.07884298927553, &
       1.328765142155775,9.57927586700208, &
       20.57098635060782,8.65565771204697, &
       -23.94605298593593,-4.258451114175736, &
       -1.072234560754618,1.461005254873122, &
       15.24914497697904,4.248041948327185 /
!
      data vbm(15,1:20) / &
       2.0074767391530557,-0.13300715407205188, &
       0.08421706263334022,1.0713833305728522, &
       0.16978619081169366,1.6320856047728134, &
       -0.014336307912426814,0.37965026966983934, &
       -63.53838090113402,4.839663833914924, &
       3.882749378977139,4.753876209534907, &
       18.697405235967835,-29.532698337703987, &
       -31.31194564874733,-2.9361386230989397, &
       4.553784170968756,-5.314482967771423, &
       6.338268974035085,8.999014230913438 /
!
!    fm = 0.75
!
      data vm(16,1:20) / &
       1.776960187792968,-0.7095598478437503, &
       0.6474974757597657,2.19861746036987, &
       0.3807881827185059,0.94951127725586, &
       0.1477442199999997,0.6542037461914064, &
       -48.78291749420161,-16.38933840107416, &
       1.335337506445313,9.65318304972656, &
       20.55476735449216,8.53470494046875, &
       -23.78040297187498,-4.167801433007817, &
       -1.065760895937499,1.436230282499998, &
       15.23978196152344,4.303312619787586 /
!
      data vbm(16,1:20) / &
       1.9413441419131594,-0.5049472835755147, &
       0.26497351630208316,1.1432622465887607, &
       0.25367726103644006,1.4401699918313549, &
       0.05201400426822691,0.4891893699410206, &
       -64.85745175492326,6.498146802728398, &
       4.337550191674632,0.029548206534631305, &
       23.420265770436774,-18.656063859250136, &
       -26.80857718315157,0.040037426996491376, &
       4.714428315830718,-6.1077393595000355, &
       3.5447028425963043,8.998433691340423 /
!
!    fm = 0.8
!
      data vm(17,1:20) / &
       1.79322326703252,-0.7293955903237455, &
       0.6524793684524797,2.221090476405281, &
       0.3792445396644045,0.945789732336901, &
       0.1438597673821279,0.6501996334180939, &
       -48.56753770817832,-15.65933567752927, &
       1.347472416019409,9.74498344553637, &
       20.52131557344067,8.38753043128974, &
       -23.58176112140332,-4.071796164464213, &
       -1.063067273056987,1.41023887588886, &
       15.23338175107462,4.365806299534065 /
!
      data vbm(17,1:20) / &
       1.9302452318300933,-0.852951492526195, &
       0.3959214391718825,1.2179370094826645, &
       0.35429785060954677,1.299423631777687, &
       0.06411587583140571,0.5953019409078614, &
       -63.591094720220845,9.52458307472122, &
       3.9486119177739134,-3.6935258168374183, &
       28.18323197405138,-10.84654875537802, &
       -21.36141608219227,4.003742699637443, &
       4.165609291655199,-7.788310799129874, &
       2.515757872673094,10.054115856246224 /
!
!    fm = 0.85
!
      data vm(18,1:20) / &
       1.809790739201659,-0.753398392226074, &
       0.6553053728729949,2.244545248168946, &
       0.3820171485000003,0.940944877000977, &
       0.1392198314375001,0.6460151436641313, &
       -48.28193564849856,-14.8908941046982, &
       1.365923048447265,9.85216421640625, &
       20.46679076367182,8.21070356616212, &
       -23.35185989431153,-3.978412261040036, &
       -1.064765803271483,1.383540860839843, &
       15.22942607636106,4.440690919453123 /
!
      data vbm(18,1:20) / &
       1.824690945158358,-1.0224701209263978, &
       0.3691794928586602,1.22527764787106, &
       0.37864068383857724,1.1842894816200285, &
       0.02062242953488158,0.6492874977441266, &
       -59.24214375751758,13.925388922392322, &
       3.2016455900875123,-6.969968706290668, &
       31.507070361877837,-6.474531032349542, &
       -16.31038380476514,5.564670756021075, &
       3.5126476169498666,-9.092970157298982, &
       2.458883888335885,12.36013932503142 /
!
!    fm = 0.9
!
      data vm(19,1:20) / &
       1.822111588391397,-0.783258322876798, &
       0.6572177250761032,2.262749100604942, &
       0.3897142779239253,0.935242513662237, &
       0.1367331502602441,0.6411619118830746, &
       -47.93253809288446,-14.04953522268507, &
       1.377555833953782,9.93688437453102, &
       20.42566526134788,8.03091107166026, &
       -23.12418382285787,-3.88725752681089, &
       -1.064935205386257,1.360150435943781, &
       15.22069324776645,4.516858461693895 /
!
      data vbm(19,1:20) / &
       1.559564910391775,-1.0823458305029623, &
       0.1375482462215561,1.1075017578135444, &
       0.41608866355887975,1.051766857566625, &
       -0.055887342912950375,0.6527962777369825, &
       -51.070813584066585,19.993379644496837, &
       1.8400262011996587,-10.822932773743421, &
       34.218411655749954,-4.068833327477977, &
       -11.992929268209483,5.2946439162902745, &
       2.327701879113646,-9.316567784118169, &
       3.6669828464070724,13.602189229004901 /
!
!    fm = 0.95
!
      data vm(20,1:20) / &
       1.829432893074816,-0.821264250926512, &
       0.6566431623090833,2.27642809366214, &
       0.4051807883283995,0.926469588593079, &
       0.1367587925338443,0.6353655200204377, &
       -47.55020933452914,-13.22722187190671, &
       1.359783044628289,10.02544116333976, &
       20.41217259599056,7.860814446748875, &
       -22.92923357820408,-3.788525103368127, &
       -1.060362392568311,1.346813822506071, &
       15.19983614838365,4.58234936299015 /
!
      data vbm(20,1:20) / &
       1.4440415836960916,-1.1118801469051225, &
       -0.24131271266938503,0.5436925664278842, &
       0.5438178906177233,1.0815083368294336, &
       -0.06582725902585534,0.5923741177286641, &
       -47.931312608227685,26.76489912248597, &
       2.4371207962908974,-12.736608147977805, &
       33.42304757059201,-3.184062077650802, &
       -12.163214267960935,5.573309385568713, &
       1.9040506168618645,-9.506701820528814, &
       4.908946246863524,14.948551545336201 /
!
!    fm = 1.0
!
      data vm(21,1:20) / &
       1.826847832031249,-0.861310124999999, &
       0.6575781562499995,2.2825501875, &
       0.4241429609375004,0.915902718750001, &
       0.1407968007812502,0.6294216650390628, &
       -47.27620127343747,-12.47474156250001, &
       1.301349406249996,10.06786334436035, &
       20.47157314453126,7.759854870117187, &
       -22.84431474609372,-3.688952213541671, &
       -1.045452443847656,1.344488212402342, &
       15.16781977294922,4.613708160644531 /
!
      data vbm(21,1:20) / &
       1.0720427495888059,-1.0352701554051484, &
       -1.1587645266235862,-1.539881808285635, &
       0.7536017884881802,1.3101251808957506, &
       0.12907305691341753,0.41497311224424527, &
       -34.27184170047316,27.21047057820418, &
       3.1070862741557397,-16.217231447363456, &
       28.24681814543993,5.152204009143581, &
       -13.451163138468832,5.097927633723974, &
       -0.08081708900093537,-8.719713606189327, &
       9.245808116510783,13.203793959198736 /
!
!    interpolate to fx
!
      j1 = INT(20*fx)+1
      if(j1.eq.21)j1=20
!      write(*,*)"fx = ",fx," j1 = ",j1,fm(j1),fm(j1+1)
      j2 = j1+1
      df = (fx-fm(j1))/(fm(j2)-fm(j1))
!      write(*,*)"df = ",df
      do i=1,20
        v(i) = vm(j1,i) + (vm(j2,i)-vm(j1,i))*df  
        vb(i) = vbm(j1,i) + (vbm(j2,i)-vbm(j1,i))*df    
      enddo
!
        u1_r = v(1)
        u1_i = v(2)
        u2_r = v(3)
        u2_i = v(4)
        u3_r = v(5)
        u3_i = v(6)
        u4_r = v(7)
        u4_i = v(8)
        u5_r = v(9)
        u5_i = v(10)
        u6_r = v(11)
        u6_i = v(12)
        u7_r = v(13)
        u7_i = v(14)
        u8_r = v(15)
        u8_i = v(16)
        u9_r = v(17)
        u9_i = v(18) 
        u10_r = v(19)
        u10_i = v(20)
!
        ub1_r = vb(1)
        ub1_i = vb(2)
        ub2_r = vb(3)
        ub2_i = vb(4)
        ub3_r = vb(5)
        ub3_i = vb(6)
        ub4_r = vb(7)
        ub4_i = vb(8)
        ub5_r = vb(9)
        ub5_i = vb(10)
        ub6_r = vb(11)
        ub6_i = vb(12)
        ub7_r = vb(13)
        ub7_i = vb(14)
        ub8_r = vb(15)
        ub8_i = vb(16)
        ub9_r = vb(17)
        ub9_i = vb(18) 
        ub10_r = vb(19)
        ub10_i = vb(20)
!
! debug
!      write(*,*)"debug get_u","ft = ",fx
!      do i=1,20
!       write(*,*)i,v(i)
!      enddo
!
      END SUBROUTINE get_u
!
!--------------------------------------------------------------
!

