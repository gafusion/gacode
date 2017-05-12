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
      LOGICAL :: USE_TTF=.FALSE.
      INTEGER :: i,is,k,j,j1,j2,jmax1,jmax2
      REAL :: test,testmax1,testmax2
      REAL :: gammamax1,kymax1,gammamax2,kymax2,ky0,ky1,ky2
      REAL :: f0,f1,f2,a,b,c,x0,x02,dky,xmax
      REAL :: gamma0,gamma1,gammaeff,dgamma
      REAL :: cnorm, phinorm, czf, cz1, cz2, kyetg
      REAL :: kyhigh, kycut
      REAL :: cky,sqcky,delta,ax,ay,kx
      REAL :: mix1,mix2,mixnorm,gamma_ave
      REAL :: vzf,dvzf,vzf1,vzf2,vzf3,vzf4
      REAL :: bz1,bz2
      REAL,DIMENSION(nkym) :: gamma_net=0.0
      REAL,DIMENSION(nkym) :: gamma=0.0
      REAL,DIMENSION(nkym) :: gamma_mix=0.0
      REAL,PARAMETER :: small=1.0E-10
      !
      ! model fit parameters
      ! need to set alpha_zf_in = 1.0
      ! Miller geometry values igeo=1
      if(xnu_model_in.eq.3)USE_TTF=.TRUE.
      czf = alpha_zf_in
      bz1=0.0
      bz2=0.0
      kyetg=1.28
      cnorm=14.21
      cz1=0.48*czf
      cz2=1.0*czf
      ax=0.0
      ay=0.0
      if(alpha_quench_in.eq.0.0)then
      !spectral shift model parameters
        ax = 1.15
        ay = 0.56
      endif
      !   write(*,*)" ax= ",ax," ay= ",ay
      do i=1,nky
         kx=spectral_shift_out(i)
         gamma_net(i) = eigenvalue_spectrum_out(1,i,1)/(1.0 + (ax*kx)**4)
      !   write(*,*)i,"gamma_net = ",gamma_net(i)
      enddo
      if(USE_MIX)then
        kyetg=1.9*zs(2)/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
        cky=3.0
        sqcky=SQRT(cky)
        cnorm = 14.29
        cz1=0.48*czf
        cz2=0.92*czf  
        if(USE_TTF)then
           bz1=1.0
           bz2=0.18
           cz1=0.48*czf
           cz2=1.6*czf
!           cz2=1.35*czf*(1.563824/q_in)
!           cz1=0.48*czf*((3.0*2.0/0.5)*(rmin_input/(Rmaj_input*q_in)))**2
!           cz2=1.35*czf*((3.098143*1.563824/0.600049)*(rmin_input/(Rmaj_input*q_in)))
           kyetg=0.8*0.04/SQRT(taus(1)*mass(1))  ! fixed streamer size to electron gyroradius
           cky=3.0
           sqcky=SQRT(cky)
        endif  
      endif   
      if(igeo.eq.0)then ! s-alpha 
       cnorm=14.63
       cz1=0.90*czf
       cz2=1.0*czf
      endif
      !
      ! renormalize the fluxes and intensities to the phi-norm from the v-norm
      do j=1,nky
      !  write(*,*)"spectal_shift_out(",j,") = ",spectral_shift_out(j)
         do i=1,nmodes_in
            phinorm=1.0
            if(ABS(field_spectrum_out(2,j,i)).gt.small)phinorm=field_spectrum_out(2,j,i)
            do is=1,ns
               intensity_spectrum_out(1,is,j,i) = intensity_spectrum_out(1,is,j,i)/phinorm
               intensity_spectrum_out(2,is,j,i) = intensity_spectrum_out(2,is,j,i)/phinorm
               do k=1,3
                  flux_spectrum_out(1,is,k,j,i) = flux_spectrum_out(1,is,k,j,i)/phinorm
                  flux_spectrum_out(2,is,k,j,i) = flux_spectrum_out(2,is,k,j,i)/phinorm
                  flux_spectrum_out(3,is,k,j,i) = flux_spectrum_out(3,is,k,j,i)/phinorm
                  flux_spectrum_out(4,is,k,j,i) = flux_spectrum_out(4,is,k,j,i)/phinorm
                  flux_spectrum_out(5,is,k,j,i) = flux_spectrum_out(5,is,k,j,i)/phinorm
              enddo
           enddo
        enddo
      enddo
      ! find the maximum of gamma/ky 
      gammamax1= gamma_net(1)
      kymax1 = ky_spectrum(1)
      testmax1 = gammamax1/kymax1
      testmax2 = 0.0
      jmax1=1
      jmax2=0
      kycut=0.8/SQRT(taus_in(2)*mass_in(2))
      kyhigh=0.15/SQRT(taus_in(1)*mass_in(1))
!      write(*,*)" kycut = ",kycut," kyhigh = ",kyhigh
      j1=1
      j2=1
      ! find the low and high ky peaks of gamma/ky
      do j=2,nky
         ky0 = ky_spectrum(j)
         if(ky0 .le. kycut)j1=j1+1
         if(ky0 .lt. kyhigh)j2=j2+1
!         write(*,*)"j=",j,"ky = ",ky0," gamma_net = ",gamma_net(j)
         test = gamma_net(j)/ky0
         if(ky0 .le. kycut)then
            if(test .gt. testmax1)then
              testmax1=test
              jmax1=j
            endif 
         endif 
         if(ky0 .gt. kycut)then
           if(test .gt. testmax2)then    
             testmax2 = test
             jmax2=j
           endif
         endif        
      enddo
!      write(*,*)"j1 = ",j1," j2 = ",j2
      ! handle exceptions
      if(j1.eq.nky)then  ! the maximum ky in the ky-spectrum is less than kycut
         j1=nky-1        ! note that j2=nky in this case     
      endif
      if(jmax2.eq.0)jmax2=j2    ! there was no high-k peak set kymax2 to kyhigh or the highest ky 
!      if(jmax1.lt.nky)then
!        test=gamma_net(jmax1+1)/ky_spectrum(jmax1+1)
!        if(testmax1.le.test)then
          ! there is no low-k peak 
!          write(*,*)" gammamax1/kymax1 = ",gammamax1/kymax1, test
!          jmax1=jmax2  
!        endif
!      endif
      gammamax2 = gamma_net(jmax2)
      kymax2 = ky_spectrum(jmax2)
      gammamax1 = gamma_net(jmax1)
      kymax1 = ky_spectrum(jmax1)
 !     write(*,*)" jmax1 = ",jmax1," jmax2= ",jmax2
 !     write(*,*)" g/k 1 = ",gammamax1/kymax1," g/k 2 = ",gammamax2/kymax2
      !interpolate to find a more accurate low-k maximum gamma/ky 
      ! this is cut of at j1 since a maximum may not exist in the low-k range
      if(jmax1.gt.1.and.jmax1.lt.j1)then
!      write(*,*)"refining low-k maximum"
         f0 =  gamma_net(jmax1-1)/ky_spectrum(jmax1-1)
         f1 =  gamma_net(jmax1)/ky_spectrum(jmax1)
         f2 =  gamma_net(jmax1+1)/ky_spectrum(jmax1+1)
         dky = (ky_spectrum(jmax1+1)-ky_spectrum(jmax1-1))
         x0 = (ky_spectrum(jmax1)-ky_spectrum(jmax1-1))/dky
         a = f0
         x02 = x0*x0
         b = (f1 - f0*(1-x02)-f2*x02)/(x0-x02)
         c = f2 - f0 - b
         xmax = -b/(2.0*c)
         if(xmax .ge. 1.0)then
           kymax1 = ky_spectrum(jmax1+1)
           gammamax1 = f2*kymax1
         elseif(xmax.lt.0.0)then
           kymax1 = ky_spectrum(jmax1-1)
           gammamax1 = f0*kymax1
         else
           kymax1 = ky_spectrum(jmax1-1)+dky*xmax
           gammamax1 = (a+b*xmax+c*xmax*xmax)*kymax1
         endif     
      endif
!      write(*,*)"gammamax1 = ",gammamax1," kymax1 = ",kymax1
!      write(*,*)"gammamax2 = ",gammamax2," kymax2 = ",kymax2
      ! compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
      ! note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
      vzf1 = gammamax1/kymax1
      vzf2 = gammamax2/kymax2
      dvzf = MAX(vzf2-vzf1,0.0)
      dgamma= dvzf*kymax1
      gamma1=gammamax1
!      write(*,*)"dvzf = ",dvzf," vzf1 = ",vzf1," vzf2 = ",vzf2
!      write(*,*)"gamma1 = ",gamma1
      vzf3 = MAX(vzf1 - bz1*dvzf,0.0)
      vzf4 = MAX(vzf1 - bz2*dvzf,0.0)
      do j=1,nky
! include zonal flow effects on growthrate model
!          gamma=0.0
          gamma0 = gamma_net(j)
          ky0=ky_spectrum(j)
          if(USE_TTF)then
            if(ky0.lt.kymax1)then
              gamma(j) = MAX(gamma0  - cz1*(kymax1 - ky0)*vzf3,0.0)
            else
              gamma(j) = MAX(gamma0  - cz2*vzf4*ky0 ,gamma1)          
            endif 
          else
            if(ky0.lt.kymax1)then
              gamma(j) = MAX(gamma0  - cz1*(kymax1 - ky0)*vzf1,0.0)
            else
              gamma(j) = cz2*gammamax1  +  Max(gamma0 - cz2*vzf1*ky0,0.0)          
            endif 
          endif
          gamma_mix(j) = gamma(j)
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
        gamma_mix(j) = gamma_ave/mixnorm  
!        write(*,*)j,ky0,gamma(j),gamma_mix(j)
      enddo  
    endif      
! intensity model
      do j=1,nky
        gamma0 = eigenvalue_spectrum_out(1,j,1)
        ky0 = ky_spectrum(j)
        kx = spectral_shift_out(j)
        do i=1,nmodes_in
          gammaeff = 0.0
          if(gamma0.gt.small)gammaeff = gamma_mix(j)*(eigenvalue_spectrum_out(1,j,i)/gamma0)**2
          if(USE_TTF)then
            if(ky0.gt.kyetg)gammaeff = gammaeff*(ky0/kyetg)
          else
            if(ky0.gt.kyetg)gammaeff = gammaeff*SQRT(ky0/kyetg)
          endif
          field_spectrum_out(2,j,i) = (cnorm*gammaeff*gammaeff/ky0**4)/(1.0+ay*kx**2)**2
        enddo
     enddo
     ! recompute the intensity and flux spectra
      do j=1,nky
         do i=1,nmodes_in
            do is=1,ns
               intensity_spectrum_out(1,is,j,i) = intensity_spectrum_out(1,is,j,i)*field_spectrum_out(2,j,i)
               intensity_spectrum_out(2,is,j,i) = intensity_spectrum_out(2,is,j,i)*field_spectrum_out(2,j,i)
               do k=1,3
                  flux_spectrum_out(1,is,k,j,i) = flux_spectrum_out(1,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(2,is,k,j,i) = flux_spectrum_out(2,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(3,is,k,j,i) = flux_spectrum_out(3,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(4,is,k,j,i) = flux_spectrum_out(4,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(5,is,k,j,i) = flux_spectrum_out(5,is,k,j,i)*field_spectrum_out(2,j,i)
              enddo
           enddo
        enddo
      enddo    

      END SUBROUTINE get_multiscale_spectrum
      
