!
!--------------------------------------------------------------
!
      SUBROUTINE get_multiscale_spectrum
!
!***********************************************************************
!  questions  should be addressed to
!  Gary Staebler 858-455-3466  or email: gary.staebler@gat.com
!***********************************************************************
!
!  TGLF multiscale saturation rule: sat_rule_in = 1 option
!  April 8, 2016: G.M. Staebler, J. Candy, N. T. Howard, and C. Holland
!                  Physics of Plasmas, 23 (2016) 062518
!  June 22, 2017: Retuned after coding error to Laplacian terms in Ampere 
!                 and Poisson equations was fixed
!  TGLF multiscale saturation rule: sat_rule_in = 2 option
!  Aug. 6, 2020 based on papers submitted to PPCF & NF
!
!***********************************************************************
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      USE tglf_xgrid
      IMPLICIT NONE
      !
      LOGICAL :: USE_MIX=.TRUE.
      LOGICAL :: first_pass = .TRUE.
      LOGICAL :: USE_SUB1=.FALSE.
      INTEGER :: i,is,k,l,j,j1,j2,jmax1,jmax2,m
      
      INTEGER :: expsub=2, exp_ax, jmax_mix
      REAL :: test,testmax1,testmax2
      REAL :: gammamax1,kymax1,gammamax2,kymax2,ky0,ky1,ky2
!      REAL :: f0,f1,f2,a,b,c,x0,x02,dky,xmax
      REAL :: gamma0,gamma1,gammaeff,dgamma
      REAL :: cnorm, phinorm, czf, cz1, cz2, kyetg
      REAL :: kycut
      REAL :: cky,sqcky,delta,ax,ay,kx
      REAL :: mix1,mix2,mixnorm,gamma_ave
      REAL :: vzf,dvzf,vzf1,vzf2
      REAL :: vzf_mix, kymax_mix
      REAL :: etg_streamer
      REAL :: kxzf, sat_geo_factor
      REAL :: b0,b1,b2,b3
      REAL :: d1,d2,Gq,kx_width,dlnpdr,ptot
      REAL :: measure
      ! SAT3
      REAL :: kmax, gmax, kmin, aoverb, coverb, k0, kP, c_1
      REAL :: kT, doversig0, eoversig0, foversig0, vzf_out_fp, scal
      REAL :: x_ITG, x_TEM, Y_ITG, Y_TEM, x, Y, sig_ratio, Fky, YT
      REAL,DIMENSION(nkym) :: sum_W_i, abs_W_ratio, dummy_interp
      REAL,DIMENSION(nkym) :: gamma_fp=0.0
      REAL,DIMENSION(nmodes_in) :: Ys, xs, QLA_P, QLA_E, QLA_O, YTs
      !
      REAL,DIMENSION(nkym) :: gamma_net=0.0
      REAL,DIMENSION(nkym) :: gamma=0.0
      REAL,DIMENSION(nkym) :: gamma_kymix=0.0
      REAL,PARAMETER :: small=1.0E-10
      !
      
      if(jmax_out.eq.0)then
         first_pass = .TRUE.
      else
        first_pass = .FALSE.
      endif
      if(first_pass)then  ! first pass for spectral shift model or only pass for quench rule
        gamma_net(:) = eigenvalue_spectrum_out(1,:,1) 
        CALL get_zonal_mixing(nky,ky_spectrum,gamma_net,vzf_mix,kymax_mix,jmax_mix) 
        vzf_out = vzf_mix 
        kymax_out = kymax_mix 
        jmax_out = jmax_mix
        gamma_net(:) = eigenvalue_spectrum_out(1,:,1)
!        write(*,*)"FIRST PASS: vzf_out = ",vzf_out,"  jmax_out = ",jmax_out
!        write(*,*)"FIRST PASS: kymax_out = ",kymax_out,"  gammamax_out = ",kymax_out*vzf_out
      endif
      ! model fit parameters
      ! need to set alpha_zf_in = 1.0
      ! Miller geometry values igeo=1
      if(rlnp_cutoff_in.gt.0.0)then
        if(beta_loc.eq.0.0)then
          dlnpdr = 0.0
          ptot = 0.0
 !        do is=1,nstotal_in   ! include all species even non-kinetic ones like fast ions
          do is=1,ns
            ptot = ptot + as(is)*taus(is)
            dlnpdr = dlnpdr + as(is)*taus(is)*(rlns(is)+rlts(is))
          enddo
          dlnpdr = rmaj_input*dlnpdr/MAX(ptot,0.01)
        else
          dlnpdr = -p_prime_loc*(8.0*pi/beta_loc)*(rmin_input/q_loc)*rmaj_input
        endif
        if(dlnpdr .ge. rlnp_cutoff_in)dlnpdr = rlnp_cutoff_in
        if(dlnpdr .lt. 4.0)dlnpdr = 4.0
      else
         dlnpdr = 12.0
      endif
!         write(*,*)"dlnpdr = ",dlnpdr
!
      czf = ABS(alpha_zf_in)
!
! coefficients for SAT_RULE = 1
      if(sat_rule_in.eq.1) then
         kyetg=1.28
         cnorm=14.21
        if(USE_SUB1)then
          cnorm=12.12
          expsub=1
        endif
        cz1=0.48*czf
        cz2=1.0*czf
        cnorm = 14.29
        measure = SQRT(taus(1)*mass(2))
!
!        if(USE_X3)then
!         cnorm = 12.94  ! note this is normed to GASTD CGYRO units
!         cz1 = 0.0
!         cz2=1.4*czf
!         etg_streamer = 1.0
!         kyetg=etg_streamer*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
!         cky=3.0
!         sqcky=SQRT(cky)
!        endif
        if(USE_MIX)then
!original        kyetg=1.9*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
          cky=3.0
          sqcky=SQRT(cky)
          if(USE_SUB1)cnorm=12.12
          cz1=0.48*czf
!original         cz2=0.92*czf
! retuned June 22,2017
          cz2 = 1.0*czf
          etg_streamer=1.05
          if(alpha_quench_in .ne. 0.0)etg_streamer=2.1
!          kyetg=etg_streamer*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
          kyetg=etg_streamer/rho_ion  ! fixed to ion gyroradius
        endif
        if(igeo.eq.0)then ! s-alpha
           cnorm=14.63
           cz1=0.90*czf
           cz2=1.0*czf
        endif
     endif
! coefficents for SAT_RULE = 2
     if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
       ! SAT2 fit for CGYRO linear modes NF 2021 paper
!        b0 = 0.72
        b0 = 0.76
        b1 = 1.22
!        b2 = 24.52  ! note this is b2**2 in PPCF paper 2020
!        b2 = 11.21
!        b2 = 8.44
!        b2 = 2.13
        b2 = 3.74
        if(nmodes_in.gt.1)b2 = 3.55
!        b3 = 0.88
!        b3 = 2.4
        b3 = 1.0
!        d1 = 0.5475*((Bt0_out/B_geo0_out)**4)/grad_r0_out**2
        d1 = (Bt0_out/B_geo0_out)**4    ! PPCF paper 2020
        d1 = d1/grad_r0_out
        Gq = B_geo0_out/grad_r0_out
        d2 = b3/Gq**2
        cnorm = b2*(12.0/dlnpdr)
        kyetg = 1000.0   ! does not impact SAT2
        cky=3.0
        sqcky=SQRT(cky)
        kycut = b0*kymax_out
        cz1=0.0
        cz2 = 1.05*czf
        measure = 1.0/kymax_out     ! changed from alpha_i on Aug 8, 2021
      endif
!      write(*,*)"Bt0 = ",Bt0_out," B0 = ",B_geo0_out," gradr0 = ",grad_r0_out
!      write(*,*)"G1 = ",sat_geo1_out,"  G2 = ",sat_geo2_out,"  sat_geo0 = ",sat_geo0_out
!      write(*,*)"d1 = ",d1,"  d2 = ",d2
      if(sat_rule_in.eq.3)then
       kmax = kymax_out
       gmax = vzf_out * kymax_out
       kmin = 0.685 * kmax
       aoverb = - 1.0 / (2 * kmin)
       coverb = - 0.751 * kmax
       kT = 1.0 / rho_ion ! SAT3 used up to ky rho_av = 1.0, then SAT2
       k0 = 0.6 * kmin
       kP = 2.0 * kmin
       c_1 = - 2.42
       x_ITG = 0.8
       x_TEM = 1.0
       Y_ITG = 3.3 * (gmax ** 2) / (kmax ** 5)
       Y_TEM = 12.7 * (gmax ** 2) / (kmax ** 4)
       scal = 0.82 ! Q(SAT3 GA D) / (2 * QLA(ITG,Q) * Q(SAT2 GA D))
       Ys = 0
       xs = 0
       do k=1,nmodes_in
        sum_W_i = 0
        do is = 2, ns ! sum over ion species, requires electrons to be species 1
         sum_W_i = sum_W_i + QL_flux_spectrum_out(2,is,1,:,k)
        end do
        ! check for singularities in weight ratio near kmax
        i = 1
        do while (ky_spectrum(i) < kmax)
	     i = i + 1
	    end do
	    if(sum_W_i(i).eq.0.0 .OR. sum_W_i(i-1).eq.0.0)then
	     x = 0.5
	    else
	     abs_W_ratio = abs(QL_flux_spectrum_out(2,1,1,:,k) / sum_W_i)
         x = linear_interpolation(ky_spectrum, abs_W_ratio, kmax)
	    end if
	    xs(k) = x
        Y = mode_transition_function(x, Y_ITG, Y_TEM)  
        Ys(k) = Y
       enddo
      endif 
! coefficients for spectral shift model for ExB shear
      ax=0.0
      ay=0.0
      exp_ax = 1
      if(alpha_quench_in.eq.0.0)then
      !spectral shift model parameters
        ax = 1.15
        ay = 0.56
        exp_ax = 4
       if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
         ax = 1.21
         ay = 1.0
         exp_ax = 2
       endif
      endif
!    write(*,*)" ax= ",ax," ay= ",ay," kycut = ",kycut
!  apply the spectral shift temporal suppression factor
!
     if(first_pass .eqv. .FALSE.)then  ! second pass for spectral shift model
        do i=1,nky
          kx = spectral_shift_out(i)
          if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
            ky0=ky_spectrum(i)
            if(ky0.lt.kycut)then
              kx_width = kycut/grad_r0_out
            else
              kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
            endif
            kx = kx*ky0/kx_width
          endif
          gamma_net(i) = eigenvalue_spectrum_out(1,i,1)/(1.0 + ABS(ax*kx)**exp_ax)
!         write(*,*)ky0,kx,(ax*kx)**4
!          write(*,*)i,"gamma_net = ",gamma_net(i)
        enddo
!     write(*,*)"gammamax_out = ",vzf_out*kymax_out, gamma_net(jmax_out),"  vexb_shear = ",vexb_shear_s!!
! find the maximum of gamma_net/ky
        if(sat_rule_in.eq.1)then
          CALL get_zonal_mixing(nky,ky_spectrum,gamma_net,vzf_mix,kymax_mix,jmax_mix)
          vzf_out = vzf_mix
          kymax_out = kymax_mix
          jmax_out = jmax_mix
        else
          vzf_out_fp = vzf_out ! used for recreation of first pass for SAT3
          vzf_out = vzf_out*gamma_net(jmax_out)/MAX(eigenvalue_spectrum_out(1,jmax_out,1),small)
        endif
!        write(*,*)"2nd PASS: vzf_out = ",vzf_out,"  jmax_out = ",jmax_out
!        write(*,*)"2nd PASS: kymax_out = ",kymax_out,"  gammamax_out = ",kymax_out*vzf_out
!        write(*,*)"2nd PASS: gamma_net(jmax) = ",gamma_net(jmax_out)
      endif   ! second pass complete


! compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
      ! note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
      gammamax1= vzf_out*kymax_out
      kymax1 = kymax_out
      jmax1 = jmax_out
      vzf1 = vzf_out
!      vzf2 = gammamax2/kymax2
!      dvzf = MAX(vzf2-vzf1,0.0)
!      dgamma= dvzf*kymax1
!      gamma1=gammamax1
!      write(*,*)"dvzf = ",dvzf," vzf1 = ",vzf1," vzf2 = ",vzf2
!      write(*,*)"gammamax1 = ",gammamax1
!      write(*,*)"kymax1 = ",kymax1,"  kycut = ",kycut
      do j=1,nky
! include zonal flow effects on growthrate model
!          gamma=0.0
          gamma0 = gamma_net(j)
          ky0=ky_spectrum(j)
          if(sat_rule_in.eq.1)then
            if(ky0.lt.kymax1)then
              gamma(j) = MAX(gamma0  - cz1*(kymax1 - ky0)*vzf1,0.0)
            else
              gamma(j) = cz2*gammamax1  +  Max(gamma0 - cz2*vzf1*ky0,0.0)          
            endif
          elseif(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
! note that the gamma mixing is for gamma_model not G*gamma_model
            if(ky0.lt.kymax1)then
              gamma(j) = gamma0
            else
              gamma(j) = gammamax1 +Max(gamma0 - cz2*vzf1*ky0,0.0)
            endif
          endif
          gamma_kymix(j) = gamma(j)     ! used if USE_MIX=F 
      enddo
    if(USE_MIX)then
      !mix over ky > kymax with integration weight = sqcky*ky0**2/(ky0**2 + cky*(ky-ky0)**2)
      do j=jmax1+2,nky
        gamma_ave = 0.0
        ky0 = ky_spectrum(j)
        mixnorm = ky0*(ATAN(sqcky*(ky_spectrum(nky)/ky0-1.0))-  &
                  ATAN(sqcky*(ky_spectrum(jmax1+1)/ky0-1.0)))
        do i=jmax1+1,nky-1
          ky1 = ky_spectrum(i)
          ky2 = ky_spectrum(i+1)
          mix1 = ky0*(ATAN(sqcky*(ky2/ky0-1.0))- ATAN(sqcky*(ky1/ky0-1.0)))
          delta = (gamma(i+1)-gamma(i))/(ky2-ky1)
          mix2 = ky0*mix1 + (ky0*ky0/(2.0*sqcky))*(LOG(cky*(ky2-ky0)**2+ky0**2)- &
                 LOG(cky*(ky1-ky0)**2+ky0**2))
          gamma_ave = gamma_ave + (gamma(i)-ky1*delta)*mix1 + delta*mix2
        enddo  
        gamma_kymix(j) = gamma_ave/mixnorm
!        write(*,*)j,ky0,gamma(j),gamma_kymix(j)
      enddo  
    endif
    
    ! recreate first pass for ExB rule in SAT3
    if(sat_rule_in.eq.3)then
     if(first_pass .eqv. .FALSE.)then
      gamma_fp=0.0
      do j=1,nky
       gamma0 = eigenvalue_spectrum_out(1,j,1)
       ky0=ky_spectrum(j)
       if(ky0.lt.kymax1)then
        gamma(j) = gamma0
       else
        gamma(j) = (gammamax1 * (vzf_out_fp/vzf_out)) + Max(gamma0 - cz2*vzf_out_fp*ky0,0.0)
       endif
       gamma_fp(j) = gamma(j)
      enddo
      if(USE_MIX)then
       do j=jmax1+2,nky
        gamma_ave = 0.0
        ky0 = ky_spectrum(j)
        mixnorm = ky0*(ATAN(sqcky*(ky_spectrum(nky)/ky0-1.0))-  &
                  ATAN(sqcky*(ky_spectrum(jmax1+1)/ky0-1.0)))
        do i=jmax1+1,nky-1
         ky1 = ky_spectrum(i)
         ky2 = ky_spectrum(i+1)
         mix1 = ky0*(ATAN(sqcky*(ky2/ky0-1.0))- ATAN(sqcky*(ky1/ky0-1.0)))
         delta = (gamma(i+1)-gamma(i))/(ky2-ky1)
         mix2 = ky0*mix1 + (ky0*ky0/(2.0*sqcky))*(LOG(cky*(ky2-ky0)**2+ky0**2)- &
                 LOG(cky*(ky1-ky0)**2+ky0**2))
         gamma_ave = gamma_ave + (gamma(i)-ky1*delta)*mix1 + delta*mix2
        enddo  
        gamma_fp(j) = gamma_ave/mixnorm
       enddo  
      endif 
     else
      gamma_fp = gamma_kymix
     end if
    end if
    
    ! generate SAT2 potential for SAT3 to connect to for electron scale
    if(sat_rule_in.eq.3)then 
     if(ky_spectrum(nky) >= kT)then
      dummy_interp=0.0
      k = 1
      do while (ky_spectrum(k) < kT)
       k=k+1
      end do
      do l=k-1,k
       gamma0 = eigenvalue_spectrum_out(1,l,1)
	   ky0 = ky_spectrum(l)
	   kx = spectral_shift_out(l)
	   if(ky0.lt. kycut)then
	    kx_width = kycut/grad_r0_out
	    sat_geo_factor = SAT_geo0_out*d1*SAT_geo1_out
	   else
	    kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
	    sat_geo_factor = SAT_geo0_out*(d1*SAT_geo1_out*kycut +  &
					 (ky0 - kycut)*d2*SAT_geo2_out)/ky0
	   endif
	   kx = kx*ky0/kx_width
	   gammaeff = 0.0
	   if(gamma0.gt.small)gammaeff = gamma_fp(l)
	   ! potentials without multimode and ExB effects, added later
	   dummy_interp(l) = scal*measure*cnorm*(gammaeff/(kx_width*ky0))**2
	   if(units_in.ne.'GYRO')dummy_interp(l) = sat_geo_factor*dummy_interp(l)
	  end do
	  YT = linear_interpolation(ky_spectrum, dummy_interp, kT)
	  do l=1,nmodes_in
	   YTs(l) = YT
	  end do
	 else
	  if(aoverb*(kP**2)+kP+coverb-((kP-kT)*(2*aoverb*kP+1))==0)then
	   YTs=0.0
	  else
	   do l=1,nmodes_in
	    YTs(l) = Ys(l)*(((aoverb*(k0**2)+k0+coverb)/(aoverb*(kP**2)+ &
	                  kP+coverb-((kP-kT)*(2*aoverb*kP+1))))**abs(c_1))
	   end do
	  end if
	 end if
	end if
  
    
    
! intensity model  
       do j=1,nky 
        gamma0 = eigenvalue_spectrum_out(1,j,1) 
        ky0 = ky_spectrum(j) 
        kx = spectral_shift_out(j) 
        if(sat_rule_in.eq.1)then
          sat_geo_factor = SAT_geo0_out
          kx_width = ky0
        endif
        if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
          if(ky0.lt. kycut)then
            kx_width = kycut/grad_r0_out
            sat_geo_factor = SAT_geo0_out*d1*SAT_geo1_out
           else
             kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
             sat_geo_factor = SAT_geo0_out*(d1*SAT_geo1_out*kycut +  &
                              (ky0 - kycut)*d2*SAT_geo2_out)/ky0
           endif
           kx = kx*ky0/kx_width
    !       write(*,*)"sat2 ",B_geo0_out,grad_r0_out,SAT_geo1_out,SAT_geo2_out,sat_geo_factor
        endif
        if(sat_rule_in.eq.1 .OR. sat_rule_in.eq.2)then
         do i=1,nmodes_in
           gammaeff = 0.0
           if(gamma0.gt.small)gammaeff = &
                gamma_kymix(j)*(eigenvalue_spectrum_out(1,j,i)/gamma0)**expsub
           if(ky0.gt.kyetg)gammaeff = gammaeff*SQRT(ky0/kyetg)
           field_spectrum_out(2,j,i) = measure*cnorm*((gammaeff/(kx_width*ky0))/(1.0+ay*kx**2))**2
           if(units_in.ne.'GYRO')field_spectrum_out(2,j,i) = sat_geo_factor*field_spectrum_out(2,j,i)
         enddo       
        elseif(sat_rule_in.eq.3)then ! SAT3
         if(gamma_fp(j)==0)then
          Fky=0.0
         else
		  Fky =  ((gamma_kymix(j) / gamma_fp(j))**2) / ((1.0 + ay*(kx**2))**2)
		 end if
         do i=1,nmodes_in
          field_spectrum_out(2,j,i) = 0.0  
          if(gamma0.gt.small) then
           if (ky0 <= kP) then ! initial quadratic
			sig_ratio = (aoverb * (ky0 ** 2) + ky0 + coverb) / (aoverb * (k0 ** 2) + k0 + coverb)	
			field_spectrum_out(2,j,i) = Ys(i) * (sig_ratio ** c_1) * Fky *(eigenvalue_spectrum_out(1,j,i)/gamma0)**(2 * expsub)
		   else if (ky0 <= kT) then ! connecting quadratic
		    if(YTs(i)==0.0 .OR. kP==kT)then
		     field_spectrum_out(2,j,i) = 0.0
		    else
			 doversig0 = ((Ys(i) / YTs(i))**(1.0/abs(c_1)))-((aoverb*(kP**2)+kP+coverb-((kP-kT)*(2*aoverb*kP+1)))/(aoverb*(k0**2)+k0+coverb))
			 doversig0 = doversig0 * (1.0/((kP-kT)**2))
			 eoversig0 = - 2 * doversig0 * kP + ((2 * aoverb * kP + 1)/(aoverb * (k0 ** 2) + k0 + coverb))
			 foversig0 = ((Ys(i) / YTs(i))**(1.0/abs(c_1))) - eoversig0 * kT - doversig0 * (kT ** 2)
			 sig_ratio = doversig0 * (ky0 ** 2) + eoversig0 * ky0 + foversig0
			 field_spectrum_out(2,j,i) = Ys(i) * (sig_ratio ** c_1) * Fky *(eigenvalue_spectrum_out(1,j,i)/gamma0)**(2 * expsub)
			end if
		   else ! SAT2 for electron scale
			gammaeff = gamma_kymix(j)*(eigenvalue_spectrum_out(1,j,i)/gamma0)**expsub
			if(ky0.gt.kyetg)gammaeff = gammaeff*SQRT(ky0/kyetg)
			field_spectrum_out(2,j,i) = scal*measure*cnorm*((gammaeff/(kx_width*ky0))/(1.0+ay*kx**2))**2
			if(units_in.ne.'GYRO')field_spectrum_out(2,j,i) = sat_geo_factor*field_spectrum_out(2,j,i)
		   end if
          end if
         enddo 
        endif
       enddo
       
	  !SAT3 QLA here
	  QLA_P=0.0
	  QLA_E=0.0
	  if(sat_rule_in.eq.3)then
	   do k=1,nmodes_in
	    ! factor of 2 included for real symmetry
	    QLA_P(k) = 2 * mode_transition_function(xs(k), 1.1, 0.6)
	    QLA_E(k) = 2 * mode_transition_function(xs(k), 0.75, 0.6)
	   enddo
	   QLA_O = 2 * 0.8
	  else
	   QLA_P = 1.0
	   QLA_E = 1.0
	   QLA_O = 1.0 
	  endif
	   
     ! recompute the intensity and flux spectra
      do j=1,nky
         do i=1,nmodes_in
            phinorm=field_spectrum_out(2,j,i) 
            field_spectrum_out(1,j,i) = phinorm
            field_spectrum_out(3,j,i) = QL_field_spectrum_out(3,j,i)*phinorm
            field_spectrum_out(4,j,i) = QL_field_spectrum_out(4,j,i)*phinorm
            do is=1,ns
               intensity_spectrum_out(1,is,j,i) = QL_intensity_spectrum_out(1,is,j,i)*phinorm
               intensity_spectrum_out(2,is,j,i) = QL_intensity_spectrum_out(2,is,j,i)*phinorm
               intensity_spectrum_out(3,is,j,i) = QL_intensity_spectrum_out(3,is,j,i)*phinorm
               intensity_spectrum_out(4,is,j,i) = QL_intensity_spectrum_out(4,is,j,i)*phinorm
               do k=1,3
                  flux_spectrum_out(1,is,k,j,i) = QL_flux_spectrum_out(1,is,k,j,i)*phinorm*QLA_P(i)
                  flux_spectrum_out(2,is,k,j,i) = QL_flux_spectrum_out(2,is,k,j,i)*phinorm*QLA_E(i)
                  flux_spectrum_out(3,is,k,j,i) = QL_flux_spectrum_out(3,is,k,j,i)*phinorm*QLA_O(i)
                  flux_spectrum_out(4,is,k,j,i) = QL_flux_spectrum_out(4,is,k,j,i)*phinorm*QLA_O(i)
                  flux_spectrum_out(5,is,k,j,i) = QL_flux_spectrum_out(5,is,k,j,i)*phinorm*QLA_O(i)
              enddo
           enddo
        enddo
      enddo   
      ! SAT3 functions
      contains 
		! mode transition function
		real function mode_transition_function (x, y1, y2)
		implicit none
		
			real :: x, y1, y2, y

			if (x < x_ITG) then
				y = y1
			else if (x > x_TEM) then
				y = y2
			else
				y = y1 * ((x_TEM - x) / (x_TEM - x_ITG)) + y2 * ((x - x_ITG) / (x_TEM - x_ITG))
			end if

			mode_transition_function = y
			
		end function mode_transition_function
		
		! linear interpolation
		real function linear_interpolation (x, y, x0)
		implicit none
		
			real, dimension(NKY) :: x, y
			real :: x0, y0
			
			i = 1
			do while (x(i) < x0)
				i = i + 1
			end do
			y0 = ((y(i) - y(i-1)) * x0 + (x(i) * y(i - 1) - x(i-1) * y(i))) / (x(i) - x(i-1)) ! y = m x0 + c

			linear_interpolation = y0
		end function linear_interpolation 
      
      END SUBROUTINE get_multiscale_spectrum
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE get_zonal_mixing(nmix,ky_mix,gamma_mix,vzf_mix,kymax_mix,jmax_mix)
!
!  finds the maximum of gamma/ky spectrum vzf_out and kymax_out
!
    USE tglf_dimensions
    USE tglf_global
    USE tglf_species
    IMPLICIT NONE
    LOGICAL :: use_kymin = .false.
    INTEGER :: nmix,i,k,j,j1,j2,jmax1,jmin
    REAL :: test,testmax,peakmax
    REAL :: kx, kyhigh, kycut, kymin
    REAL :: gammamax1,kymax1,testmax1,ky0,ky1,ky2
    REAL :: f0,f1,f2,a,b,c,x1,deltaky,xmax,xmin
    REAL :: vzf1, vzf2
    REAL :: kymax_mix, vzf_mix
    INTEGER :: jmax_mix, down
    REAL, DIMENSION(nkym) :: gamma_mix, ky_mix
!  initialize output of subroutine
    vzf_mix = 0.0
    kymax_mix = 0.0
    jmax_mix = 1
    xmin = 0.0
    if(alpha_zf_in.lt.0.0)use_kymin = .true. 
!
! find the local maximum of gamma_mix/ky_mix with the largest gamma_mix/ky_mix^2
!
      gammamax1= gamma_mix(1)
      kymax1 = ky_mix(1)
      testmax = 0.0
      peakmax=0.0
      jmax_mix=1
      kycut=0.8/rho_ion
      kymin = 0.0
      jmin = 0
      if(use_kymin)kymin = 0.173*sqrt(2.0)/rho_ion
      if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
        kycut = grad_r0_out*kycut
        kymin = grad_r0_out*kymin
      endif
!      write(*,*)" kycut = ",kycut," kymin = ",kymin
      ! find the low and high ky peaks of gamma/ky
      do j=1,nky-1
       if(ky_mix(j).lt.kymin)jmin = j
       if((ky_mix(j+1).ge.kymin).and.(ky_mix(j).le.kycut))then
         j1=j
         kymax1 = ky_mix(j)
         testmax1 = gamma_mix(j)/kymax1
!         write(*,*)"j=",j,"ky = ",ky0," gamma_net = ",gamma_net(j)
         if(testmax1.gt.testmax)then
           testmax = testmax1
           jmax_mix = j
          endif
         endif
      enddo
      if(testmax.eq.0.0)jmax_mix=j1
! no unstable modes in range set kymax index to end of range
! this is cut of at j1 since a maximum may not exist in the low-k range
      kymax1 = ky_mix(jmax_mix)
      gammamax1 = gamma_mix(jmax_mix)
      if(kymax1.lt.kymin)then
          kymax1 = kymin
!interpolate to find a more accurate low-k maximum gamma/ky
         gammamax1 = gamma_mix(1)   &
          +(gamma_mix(2)-gamma_mix(1))*(kymin-ky_mix(1))/(ky_mix(2)-ky_mix(1))
      endif
!        write(*,*)" jmax_mix = ",jmax_mix,"  gammamax1 = ",gammamax1," kymax1 = ",kymax1
      if(jmax_mix .gt. 1 .and. jmax_mix .lt. j1)then
 !        write(*,*)"refining low-k maximum"
! determine kymax1 and gammamax1 bounded by the tree points f0,f1,f2
! use a quadratic fit: f = a + b x + c x^2  to f = gamma/ky centered at jmax1
         jmax1 = jmax_mix
         f0 =  gamma_mix(jmax1-1)/ky_mix(jmax1-1)
         f1 =  gamma_mix(jmax1)/ky_mix(jmax1)
         f2 =  gamma_mix(jmax1+1)/ky_mix(jmax1+1)
         deltaky = ky_mix(jmax1+1)-ky_mix(jmax1-1)
         x1 = (ky_mix(jmax1)-ky_mix(jmax1-1))/deltaky
         a = f0
         b = (f1 - f0*(1-x1*x1)-f2*x1*x1)/(x1-x1*x1)
         c = f2 - f0 - b
!         write(*,*)"f0 = ",f0,"  f1 = ",f1,"  f2 = ",f2,"  x1 = ",x1
         if(f0 .ge.f1)then
! if f0>f1 then f1 is not a local maximum
             kymax1 = ky_mix(jmax1-1)
             gammamax1 = f0*kymax1
             if(kymax1.lt.kymin)then
!interpolate to find the value of gammamax1 at kymin
               kymax1 = kymin
               xmin = (kymin - ky_mix(jmax1-1))/deltaky
               gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
             endif
        endif
        if(f0.lt.f1  )then
!        if f0<f1 then f1>f2 due to the maximum search
! use the quadratic fit to refine the local maximum:
             xmax = -b/(2.0*c)
             xmin = 0.0
             if(ky_mix(jmax1-1).lt.kymin)then
               xmin = (kymin - ky_mix(jmax1-1))/deltaky
             endif
             if(xmax .ge. 1.0)then
! if xmax >= 1  use f2 as the maximum
               kymax1 = ky_mix(jmax1+1)
               gammamax1 = f2*kymax1
             elseif(xmax.le.xmin)then
               if(xmin .gt. 0.0)then
                 kymax1 = kymin
! use the quadratic fit to determine gammamax1 at kymin
                 gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
               elseif(xmax .le. 0.0)then
! if xmax<=0 use f0 as the maximum
                 kymax1 = ky_mix(jmax1-1)
                 gammamax1 = f0*kymax1
               endif
             else
! the conditions f0<f1<f2 and xmin<xmax<1 are satisfied
! use the quadratic fit to determine gammamax1 and kymax1
               kymax1 = ky_mix(jmax1-1)+deltaky*xmax
               gammamax1 = (a+b*xmax+c*xmax*xmax)*kymax1
             endif !xmax tests
           endif !f0 < f1
       endif  ! jmax_mix > 1
       vzf_mix = gammamax1/kymax1
       kymax_mix = kymax1
!      jmax_mix = jmax1
!      write(*,*)"get_zonal_mxing"
!      write(*,*)"xmax = ",xmax, "  xmin = ",xmin
!      write(*,*)"gammamax1 = ",gammamax1," kymax1 = ",kymax1, "  kymin = ",kymin
 
 END SUBROUTINE get_zonal_mixing
