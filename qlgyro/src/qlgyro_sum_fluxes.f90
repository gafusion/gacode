SUBROUTINE qlgyro_sum_fluxes
!
!  computes spectral integrals of field, intensity and fluxes.
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      USE tglf_interface
      USE qlgyro_globals
      IMPLICIT NONE
!
      INTEGER :: i,j,is,imax
      REAL :: dky, qlgyro_ky_in
      REAL :: phi_bar0,phi_bar1
      REAl :: v_bar0,v_bar1
      REAL :: dky0,dky1,ky0,ky1
      REAL :: pflux0(nsm,3),eflux0(nsm,3)
      REAL :: stress_par0(nsm,3),stress_tor0(nsm,3)
      REAL :: exch0(nsm,3)
      REAL :: nsum0(nsm),tsum0(nsm)
      REAL :: pflux1(nsm,3),eflux1(nsm,3)
      REAL :: stress_par1(nsm,3),stress_tor1(nsm,3)
      REAL :: exch1(nsm,3)
      REAL :: nsum1(nsm),tsum1(nsm)
!
!      CALL tglf_startup  
!
! initialize fluxes
!
      do is=1,tglf_ns_in
        do j=1,3 
          pflux0(is,j) = 0.0
          eflux0(is,j) = 0.0
          stress_par0(is,j) = 0.0
          stress_tor0(is,j) = 0.0
          exch0(is,j) = 0.0
        enddo
        n_bar_sum_out(is) = 0.0
        t_bar_sum_out(is) = 0.0
        q_low_out(is) = 0.0
        nsum0(is) = 0.0
        tsum0(is) = 0.0
      enddo
      phi_bar_sum_out = 0.0
      v_bar_sum_out = 0.0
      v_bar0 = 0.0
      phi_bar0 = 0.0     
!
! compute the flux spectrum
!
!      CALL get_bilinear_spectrum_mpi
!
! sum over ky spectrum
!
      iflux_in=.TRUE. 
      dky0=0.0
      ky0=0.0 
      do i=1,nky
        qlgyro_ky_in = tglf_ky_spectrum_out(i)
        dky = tglf_dky_spectrum_out(i)
        ky1=qlgyro_ky_in
        if(i.eq.1)then
          dky1=ky1
        else
          dky = LOG(ky1/ky0)/(ky1-ky0)
          dky1 = ky1*(1.0 - ky0*dky)
          dky0 = ky0*(ky1*dky - 1.0)
        endif
! normalize the ky integral to make it independent of the 
! choice of temperature and mass scales 
        dky0 = dky0*SQRT(tglf_taus_in(1)*tglf_mass_in(2))
        dky1 = dky1*SQRT(tglf_taus_in(1)*tglf_mass_in(2))
!
! compute the field integrals
!
        v_bar1 = 0.0
        phi_bar1 = 0.0
        do imax = 1,tglf_nmodes_in
          v_bar1 = v_bar1 + field_spectrum_out(1,i,imax)
          phi_bar1 = phi_bar1 + field_spectrum_out(2,i,imax)
        enddo
        phi_bar_sum_out = phi_bar_sum_out + dky0*phi_bar0 + dky1*phi_bar1
        v_bar_sum_out = v_bar_sum_out + dky0*v_bar0 + dky1*v_bar1
        phi_bar0 = phi_bar1
        v_bar0 = v_bar1
!
! compute the intensity integrals
!
!        do is=1,ns
!          nsum1(is) = 0.0
!          tsum1(is) = 0.0
!          do imax = 1,tglf_nmodes_in
!            nsum1(is) = nsum1(is) + intensity_spectrum_out(1,is,i,imax)
!            tsum1(is) = tsum1(is) + intensity_spectrum_out(2,is,i,imax)
!          enddo
!          n_bar_sum_out(is) = n_bar_sum_out(is) &
!              + dky0*nsum0(is) + dky1*nsum1(is)
!          t_bar_sum_out(is) = t_bar_sum_out(is) &
!              + dky0*tsum0(is) + dky1*tsum1(is)
!          nsum0(is) = nsum1(is)
!          tsum0(is) = tsum1(is)
!        enddo
!
! compute the flux integrals
!
        do is=1,tglf_ns_in
          do j=1,3
            pflux1(is,j) = 0.0
            eflux1(is,j) = 0.0
            stress_tor1(is,j) = 0.0
            stress_par1(is,j) = 0.0
            exch1(is,j) = 0.0
            do imax = 1,tglf_nmodes_in
              pflux1(is,j) = pflux1(is,j) + tglf_flux_spectrum_out(1,is,j,i,imax)
              eflux1(is,j) = eflux1(is,j) + tglf_flux_spectrum_out(2,is,j,i,imax)
              stress_tor1(is,j) = stress_tor1(is,j) + &
                 tglf_flux_spectrum_out(3,is,j,i,imax)
              stress_par1(is,j) = stress_par1(is,j) + &
                 tglf_flux_spectrum_out(4,is,j,i,imax)
              exch1(is,j) = exch1(is,j) + tglf_flux_spectrum_out(5,is,j,i,imax)
           enddo !imax
            tglf_particle_flux_out(is,j) = tglf_particle_flux_out(is,j) &
              + dky0*pflux0(is,j) + dky1*pflux1(is,j)
            tglf_energy_flux_out(is,j) = tglf_energy_flux_out(is,j) &
              + dky0*eflux0(is,j) + dky1*eflux1(is,j)
            tglf_stress_tor_out(is,j) = tglf_stress_tor_out(is,j) &
              + dky0*stress_tor0(is,j) + dky1*stress_tor1(is,j)
            tglf_stress_par_out(is,j) = tglf_stress_par_out(is,j) &
              + dky0*stress_par0(is,j) + dky1*stress_par1(is,j)
            tglf_exchange_out(is,j) = tglf_exchange_out(is,j) &
              + dky0*exch0(is,j) + dky1*exch1(is,j)
!            write(*,*)is,j,i
!            write(*,*)"ky0=",ky0,"ky1=",ky1
!            write(*,*)"pflux0=",pflux0,"pflux1=",pflux1
!            write(*,*)"eflux0=",eflux0,"eflux1=",eflux1
!            write(*,*)dky0*pflux0+dky1*pflux1
!            write(*,*)dky0*eflux0+dky1*eflux1
!            write(*,*)"tglf_stress_tor_out=",tglf_stress_tor_out(is,1)
             pflux0(is,j) = pflux1(is,j)
             eflux0(is,j) = eflux1(is,j)
             stress_par0(is,j) = stress_par1(is,j)
             stress_tor0(is,j) = stress_tor1(is,j)
             exch0(is,j) = exch1(is,j)
           enddo  ! j
           if(qlgyro_ky_in*SQRT(tglf_taus_in(2)*tglf_mass_in(2)).le.1.0)then
             q_low_out(is) = tglf_energy_flux_out(is,1)+tglf_energy_flux_out(is,2)
           endif
         enddo  ! is 
!
        ky0 = ky1
      enddo  ! i
!
!      CALL tglf_shutdown
      !

      tglf_elec_pflux_out  = sum(tglf_particle_flux_out(1, :))
      tglf_elec_eflux_out  = sum(tglf_energy_flux_out(1, :))
      tglf_elec_mflux_out  = sum(tglf_stress_tor_out(1, :))
      tglf_elec_expwd_out  = sum(tglf_exchange_out(1, :))

      do i=2,tglf_ns_in
         tglf_ion_pflux_out(i-1)  = sum(tglf_particle_flux_out(i, :))
         tglf_ion_eflux_out(i-1)  = sum(tglf_energy_flux_out(i, :))
         tglf_ion_mflux_out(i-1)  = sum(tglf_stress_tor_out(i, :))
         tglf_ion_expwd_out(i-1)  = sum(tglf_exchange_out(i, :))
      end do

      if (i_proc_global == 0) then
         open(unit=1,file=trim(runfile),position='append')
         write(1,*) 'QLGYRO Heat     fluxes (Qe, Qi): ', tglf_elec_eflux_out, tglf_ion_eflux_out(1)
         write(1,*) 'QLGYRO Particle fluxes (Qe, Qi): ', tglf_elec_pflux_out, tglf_ion_pflux_out(1)
         close(1)
      end if
21 format(A35, F8.3, F8.3)
    END SUBROUTINE qlgyro_sum_fluxes

