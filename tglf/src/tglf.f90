!---------------------------------------------------------
! tglf_run.f90
!
! PURPOSE:
!  Manage standalone TGLF run by calling
!
!    call tglf_read_input()
!    call tglf_run() 
!
!  and then adding extra I/O. 
!---------------------------------------------------------

program tglf

  use tglf_pkg
  use tglf_interface
  use tglf_global

  implicit none

  integer :: i
  integer :: n
  character (len=4) :: tag(5)=(/'ion1','ion2','ion3','ion4','ion5'/)
  real :: prec

  call tglf_read_input()
  call tglf_run()

  if(units_in.eq.'GENE')then
     print 30,'GENE reference units used'
     print 30,'Conversion to TGLF units:'
     print 30,'Bunit/Bref = ',1.0/Bref_out
     print 30,'Te/Tref = ',taus_in(1)
     print 30,'mi/mref = ',mass_in(2)
     print 30,'a/Lref = ',1.0
     print 30,'cs/cref = ',SQRT(taus_in(1)/mass_in(2))
     print 30,'rhos/rhoref = ',SQRT(mass_in(2)*taus_in(1))*Bref_out
  endif
  if(alpha_zf_in.lt.0.0)print 30,' kx_geo0_out = ',kx_geo0_out, &
      ' SAT_geo0_out = ',SAT_geo0_out

! write interchange stability criteria with ELITE conventions
  Print 30,'  D(R) = ',-interchange_DR,'  D(I) = ',0.25-interchange_DM

  if (tglf_use_transport_model_in) then

     ! Output to screen

     print 20,'Gam/Gam_GB','    Q/Q_GB','Q_low/Q_GB','  Pi/Pi_GB', '    S/S_GB'
     print 10,'elec',&
          tglf_elec_pflux_out,&
          tglf_elec_eflux_out,&
          tglf_elec_eflux_low_out,&
          tglf_elec_mflux_out,&
          tglf_elec_expwd_out

     prec = abs(tglf_elec_pflux_out)+&
          abs(tglf_elec_eflux_out)+&
          abs(tglf_elec_eflux_low_out)+&
          abs(tglf_elec_mflux_out)

     do i=1,tglf_ns_in-1
        print 10,tag(i),&
             tglf_ion_pflux_out(i),&
             tglf_ion_eflux_out(i),&
             tglf_ion_eflux_low_out(i),&
             tglf_ion_mflux_out(i),&
             tglf_ion_expwd_out(i)

        prec = prec+&
             abs(tglf_ion_pflux_out(i))+&
             abs(tglf_ion_eflux_out(i))+&
             abs(tglf_ion_eflux_low_out(i))+&
             abs(tglf_ion_mflux_out(i))
     enddo

     ! Output to file

     n = tglf_ns_in-1
     open(unit=1,file='out.tglf.gbflux',status='replace')
     write(1,'(32(1pe11.4,1x))') tglf_elec_pflux_out,tglf_ion_pflux_out(1:n),&
          tglf_elec_eflux_out,tglf_ion_eflux_out(1:n),&
          tglf_elec_mflux_out,tglf_ion_mflux_out(1:n),&
          tglf_elec_expwd_out,tglf_ion_expwd_out(1:n)
     close(1)

     open(unit=1,file='out.tglf.grid',status='replace')
     write(1,'(i2)') tglf_ns_in,tglf_nxgrid_in
     close(1)

     ! write ky spectrum to file out.tglf.ky_spectrum
     CALL write_tglf_ky_spectrum

     ! write flux spectrum summed over nmodes to file out.tglf.sum_flux_spectrum
     ! this can be compared directly with the flux spectrum from CGYRO 

     CALL write_tglf_sum_flux_spectrum

     ! write density fluctuation amplitude spectrum to file out.tglf.density_spectrum
     CALL write_tglf_density_spectrum

     ! write temperature fluctuation amplitude spectrum to file out.tglf.temperature_spectrum
     CALL write_tglf_temperature_spectrum

     ! write intensity fluctuation amplitude spectrum per mode to file out.tglf.intensity_spectrum
     CALL write_tglf_intensity_spectrum

     ! write field fluctuation amplitude per mode spectrum to file out.tglf.field_spectrum
     CALL write_tglf_field_spectrum

     ! write eigenvalue spectrum to file out.tglf.eigenvalue_spectrum
     CALL write_tglf_eigenvalue_spectrum

     ! write ne-te crossphase spectrum to file out.tglf.nete_crossphase_spectrum
     CALL write_tglf_nete_crossphase_spectrum

     ! write ns-ts crossphase spectrum to file out.tglf.nsts_crossphase_spectrum
     CALL write_tglf_nsts_crossphase_spectrum

     ! write QL flux (weight) spectrum per mode to file out.tglf.QL_weight_spectrum
     CALL write_tglf_QL_flux_spectrum

     ! write saturation geometry facotor spectrum per mode to file out.tglf.
     CALL write_tglf_sat_geo_spectrum

     ! write kx/ky-spectral shift spectrum to file out.tglf.spectral_shift
     CALL write_tglf_spectral_shift

     ! write intensity fluctuation amplitude spectrum per mode to file out.tglf.scalar_saturation_parameters
     CALL write_tglf_scalar_saturation_parameters

  else

     print 10,'     ky:',tglf_ky_in
     print 10,'     ft:',get_ft()
     print 10,'Guassian width = ',get_gaussian_width()
     ! Collect linear eigenvalues
     do i=1,tglf_nmodes_in
        print 10,'(wr,wi):',tglf_eigenvalue_out(i)
     enddo

     prec = sum(abs(tglf_eigenvalue_out))

    ! write single ky-eigenmode wavefunction to file
    CALL write_wavefunction_out('out.tglf.wavefunction')

  endif

  open(unit=1,file=trim(tglf_path_in)//'out.tglf.prec')
  write(1,*) prec
  close(1)

10 format(a,10(1x,1pe11.4))
20 format(t7,a,t19,a,t31,a,t43,a,t55,a)
30 format(a,1pe11.4,a,1pe11.4)

end program tglf
