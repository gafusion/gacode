program tglf

  use tglf_mpi
  use tglf_pkg
  use tglf_interface

  implicit none

  integer :: i
  integer :: n
  character (len=4) :: tag(5)=(/'ion1','ion2','ion3','ion4','ion5'/)
  real(8) :: prec
  integer :: ierr

  ! initialize MPI
  call MPI_INIT(ierr)

  iCommTglf = MPI_COMM_WORLD
  call MPI_COMM_RANK(iCommTglf,iProcTglf,ierr)
  call MPI_COMM_SIZE(iCommTglf,nProcTglf,ierr)

  call tglf_read_input()
  call tglf_run_mpi()

  if (iProcTglf == iProc0Tglf) then

! write interchange stability criteria with ELITE conventions
  Print 30,'  D(R) = ',-interchange_DR,'  D(I) = ',0.25-interchange_DM

     if (tglf_use_transport_model_in) then

        print 20,'Gam/Gam_GB',' Q/Q_GB','Q_low/Q_GB',' Pi/Pi_GB', ' S/S_GB'
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

     ! write flux spectrum to file out.tglf.flux_spectrum

     CALL write_tglf_flux_spectrum

     ! write density fluctuation amplitude spectrum to file out.tglf.density_spectrum
     CALL write_tglf_density_spectrum

     ! write temperature fluctuation amplitude spectrum to file out.tglf.temperature_spectrum
     CALL write_tglf_temperature_spectrum

     ! write intensity fluctuation amplitude spectrum to file out.tglf.intensity_spectrum
     CALL write_tglf_intensity_spectrum

     ! write field fluctuation amplitude spectrum to file out.tglf.field_spectrum
     CALL write_tglf_field_spectrum

     ! write eigenvalue spectrum to file out.tglf.eigenvalue_spectrum
     CALL write_tglf_eigenvalue_spectrum

     ! write ne-te crossphase spectrum to file out.tglf.nete_crossphase_spectrum
     CALL write_tglf_nete_crossphase_spectrum

     else

        print 10,'     ky:',tglf_ky_in
        print 10,'Guassian width = ',get_gaussian_width()
        ! Collect linear eigenvalues
        do i=1,tglf_nmodes_in
           print 10,'(wr,wi):',tglf_eigenvalue_out(i)
        enddo

        prec = sum(abs(tglf_eigenvalue_out))

     endif

     open(unit=1,file='out.tglf.grid',status='replace')
     write(1,'(i2)') tglf_ns_in,tglf_nxgrid_in
     close(1)

     open(unit=1,file=trim(tglf_path_in)//'out.tglf.prec')
     write(1,*) prec
     close(1)

  endif

10 format(a,10(1x,1pe11.4))
20 format(t7,a,t19,a,t31,a,t43,a,t55,a)
30 format(a,1pe11.4,a,1pe11.4)


end program tglf
