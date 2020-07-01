!-------------------------------------------------------------------------------
!
      SUBROUTINE get_ky_spectrum
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      IMPLICIT NONE
!
      integer :: nk_zones,nky0,nky1,nky2 
      INTEGER :: spectrum_type=0
      INTEGER :: i
      INTEGER :: ntest
      REAL :: ky_min=0.05
      REAL :: ky_max=0.7
      REAL :: ky0,ky1,lnky,dky0
      REAL :: ky_cut
      REAL :: ky_factor

!
!  spectrum_type = 0 for linear GYRO spectrum
!  spectrum_type = 1 for APS07 spectrum
!  spectrum_type = 2 for IAEA08 spectrum
!  spectrum_type = 3 for spectrum with ky_min = ky_cut
!
      new_kyspectrum=.FALSE.
!
!
      spectrum_type = kygrid_model_in
!
      nk_zones=3
!
      nky = nky_in
!
      if(units_in.eq.'GYRO')then
        ky_factor = 1.0
      else
        ky_factor = grad_r0_out
      endif
      ky0 = ky_min
      ky1 = ky_max
      if(nk_zones.ge.2)ky1 = ky_max/SQRT(mass(1))
      if(spectrum_type.eq.0)then
!        dky0 = (ky_max-ky_min)/REAL(nky-1)
        ky1=ky_in
        dky0=ky1/REAL(nky)
        do i=1,nky
!          ky_spectrum(i) = ky0 + REAL(i-1)*dky0
          ky_spectrum(i) = REAL(i)*dky0
          dky_spectrum(i) = dky0
        enddo
!        if(use_highk)then
!          ky0=1.0
!          dky0 = (ky1-ky0)/REAL(nky-1)
!          do i=1,nky
!            ky_spectrum(nky+i) = ky0 + REAL(i-1)*dky0
!            dky_spectrum(nky+i) = dky0
!          enddo
!          nky=2*nky
!        endif
      endif
      if(spectrum_type.eq.1)then   ! APS07 spectrum
        nky=9
        ky_max = 0.9*ky_factor*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))  !k_theta*rho_ion = 0.9
        dky0 = ky_max/REAL(nky)
        do i=1,nky
          ky_spectrum(i) = REAL(i)*dky0
          dky_spectrum(i) = dky0
        enddo
        ky0 = ky_max+dky0
        ky1 = 0.4*ky_factor*ABS(zs(1))/SQRT(taus_in(1)*mass_in(1))  !k_theta*rho_e = 0.4
        dky0 = LOG(ky1/ky0)/REAL(nky_in-1)
        lnky = LOG(ky0)
        if(nky_in.gt.0)then
          do i=nky+1,nky+nky_in     
            ky_spectrum(i) = EXP(lnky)
            dky_spectrum(i) = ky_spectrum(i)*dky0
            lnky = lnky + dky0
          enddo
        endif
        nky = nky + nky_in
      endif
      if(spectrum_type.eq.2)then  ! IAEA08 spectrum
        nky1=8
        dky0=0.05*ky_factor*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))
        do i=1,nky1
          ky_spectrum(i) = REAL(i)*dky0
          dky_spectrum(i) = dky0
        enddo
!        dky0 = ky_max/REAL(nky)
        dky0=0.2/SQRT(taus_in(2)*mass_in(2))
        ky0 = ky_spectrum(nky1)
        nky2=7
        nky = nky1+nky2
        do i=nky1+1,nky
          ky_spectrum(i) = ky0 + REAL(i-nky1)*dky0
          dky_spectrum(i) = dky0
        enddo
!        ky_max = 0.9*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))  !k_theta*rho_ion = 0.9
!        ky0 = ky_max+dky0
        ky0 = ky_spectrum(nky) + dky0
        ky1 = 0.4*ky_factor*ABS(zs(1))/SQRT(taus_in(1)*mass_in(1))  !k_theta*rho_e = 0.4
        dky0 = LOG(ky1/ky0)/REAL(nky_in-1)
        lnky = LOG(ky0)
        if(nky_in.gt.0)then
          do i=nky+1,nky+nky_in     
            ky_spectrum(i) = EXP(lnky)
            dky_spectrum(i) = ky_spectrum(i)*dky0
            lnky = lnky + dky0
          enddo
          nky = nky + nky_in
        endif
      endif
      if(spectrum_type.eq.3)then   ! ky_min=ky_in spectrum similar to APS07
!      the value of ky_in must be set externally e.g. ky_in = rhos*q/r
!        ky_max = 6.0*ky_min  ! n=6
!        nky0=6
        ky_max = 1.0*ky_factor*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))  !k_theta*rho_ion = 1.0
        nky0 = INT(ky_max/ky_in)-1
        ky_min = ky_in
        dky0 = ky_min
        ky_spectrum(1) = ky_min
        dky_spectrum(1) = ky_min
        do i=2,nky0
          ky_spectrum(i) = ky_spectrum(i-1) + dky0
          dky_spectrum(i) = dky0
        enddo
        if(ky_spectrum(nky0).lt.ky_max)then
          nky1 = 1
          ky_min = ky_spectrum(nky0)
          dky0 = (ky_max-ky_min)/REAL(nky1)
          do i=1,nky1
            ky_spectrum(nky0+i) = ky_spectrum(nky0+i-1) + dky0
            dky_spectrum(nky0+i) = dky0
          enddo 
        else
          nky1=0
          ky_max = ky_spectrum(nky0)
        endif       
        ky0 = ky_max
        ky1 = 0.4*ky_factor*ABS(zs(1))/SQRT(taus_in(1)*mass_in(1))  !k_theta*rho_e = 0.4
        dky0 = LOG(ky1/ky0)/REAL(nky_in)
        lnky = LOG(ky0+dky0)
        do i=1,nky_in     
          ky_spectrum(nky0+nky1+i) = EXP(lnky)
          dky_spectrum(nky0+nky1+i) = ky_spectrum(nky0+nky1+i)*dky0
          lnky = lnky + dky0
        enddo
        nky = nky0 + nky1 + nky_in
      endif
      if(spectrum_type.eq.4)then   ! APS07 spectrum with ky_min=0.05
        nky=11
        ky_min = 0.05*ky_factor*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))
        ky_spectrum(1) = ky_min
        dky_spectrum(1) = ky_min
        ky_spectrum(2) = 2.0*ky_min
        dky_spectrum(2) = ky_min
        ky_spectrum(3) = 3.0*ky_min
        dky_spectrum(3) = ky_min
        ky_spectrum(4) = 4.0*ky_min
        dky_spectrum(4) = ky_min
        ky_min = ky_spectrum(4)
        ky_max = 1.0*ky_factor*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))
        dky0 = 0.1*ky_factor*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))
        do i=5,nky
          ky_spectrum(i) = ky_min + REAL(i-4)*dky0
          dky_spectrum(i) = dky0
        enddo
        ky0 = 1.0*ky_factor*ABS(zs(2))/SQRT(taus_in(2)*mass_in(2))
        ky1 = 0.4*ky_factor*ABS(zs(1))/SQRT(taus_in(1)*mass_in(1))  !k_theta*rho_e = 0.4
        dky0 = LOG(ky1/ky0)/REAL(nky_in-1)
        lnky = LOG(ky0)
        if(nky_in.gt.0)then
          dky_spectrum(nky)=ky0-ky_spectrum(nky)
          do i=nky+1,nky+nky_in     
            ky_spectrum(i) = EXP(lnky)
            dky_spectrum(i) = ky_spectrum(i)*dky0
            lnky = lnky + dky0
          enddo
        endif
        nky = nky + nky_in
      endif
!
      if(spectrum_type.eq.5)then
        nky=nky_in
        ky1=ky_in
        dky0=ky1/REAL(nky)
        ntest=INT(0.5/dky0)+1
        nky=11 + nky_in - ntest
        do i=1,10
          ky_spectrum(i) = REAL(i)*0.05
          dky_spectrum(i) = 0.05
        enddo
        ky_spectrum(11) = REAL(ntest)*dky0
        dky_spectrum(i) = ky_spectrum(11)-0.5
        do i=12,nky
          ky_spectrum(i) = ky_spectrum(i-1) + dky0
          dky_spectrum(i) = dky0
        enddo
      endif


! debug
!      write(*,*)"ky_min=",ky_min,"ky_max=",ky_max
!      write(*,*)"nky = ",nky,"ky0 = ",ky0," ky1 = ",ky1
!      do i=1,nky
!      write(*,*)i,"ky=",ky_spectrum(i),"dky=",dky_spectrum(i)
!      enddo
!
      END SUBROUTINE get_ky_spectrum

