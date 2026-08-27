!
SUBROUTINE write_qlgyro_flux_spectrum
  !
  use qlgyro_globals
  use tglf_interface
  !   
  IMPLICIT NONE
  CHARACTER(24) :: fluxfile="out.qlgyro.flux_spectrum"
  INTEGER :: i,j,is,imax,jmax
  real :: ky
  REAL :: dky
  REAL :: dky0,dky1,ky0,ky1
  REAL :: pflux0(n_species,3),eflux0(n_species,3)
  REAL :: stress_par0(n_species,3),stress_tor0(n_species,3)
  REAL :: exch0(n_species,3)
  REAL :: pflux1(n_species,3),eflux1(n_species,3)
  REAL :: stress_par1(n_species,3),stress_tor1(n_species,3)
  REAL :: exch1(n_species,3)
  REAL :: pflux_out,eflux_out,mflux_tor_out,mflux_par_out,exch_out
  CHARACTER (117) :: full_fluxfile

  full_fluxfile = trim(path)//trim(iter_path)//fluxfile

  !
  OPEN(unit=33,file=full_fluxfile,status='replace')
  !
  ! initialize fluxes
  !
!  write(*,*)"ns0,ns,nky,nmodes",ns0,ns,nky,nmodes_in
!
  do is=1,n_species
     do j=1,3 
        pflux0(is,j) = 0.0
        eflux0(is,j) = 0.0
        stress_par0(is,j) = 0.0
        stress_tor0(is,j) = 0.0
        exch0(is,j) = 0.0
     enddo
  enddo
  ! set number of field
  jmax = 1
  if(tglf_use_bper_in)jmax=2
  if(tglf_use_Bpar_in)jmax=3
  ! loop over species and fields
  do is=1,n_species
     do j=1,jmax
        write(33,*)"species = ",is,"field =",j
        write(33,*)" ky,particle flux,energy flux,toroidal stress,parallel stress,exchange"

        !
        ! loop over ky spectrum
        !
        tglf_iflux_in=.TRUE. 
        dky0=0.0
        ky0=0.0 
        do i=1,n_ky
           ky = tglf_ky_spectrum_out(i)
           dky = tglf_dky_spectrum_out(i)
           ky1=ky
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
           ! compute the fluxes in the same way as tglf_TM.f90 
           !
           pflux1(is,j) = 0.0
           eflux1(is,j) = 0.0
           stress_tor1(is,j) = 0.0
           stress_par1(is,j) = 0.0
           exch1(is,j) = 0.0
           do imax = 1,n_modes
              pflux1(is,j) = pflux1(is,j) + tglf_flux_spectrum_out(1,is,j,i,imax)
              eflux1(is,j) = eflux1(is,j) + tglf_flux_spectrum_out(2,is,j,i,imax)
              stress_tor1(is,j) = stress_tor1(is,j) + &
                   tglf_flux_spectrum_out(3,is,j,i,imax)
              !stress_par1(is,j) = stress_par1(is,j) + &
              !     flux_spectrum_out(4,is,j,i,imax)
              exch1(is,j) = exch1(is,j) + tglf_flux_spectrum_out(4,is,j,i,imax)
           enddo !imax
           pflux_out = dky0*pflux0(is,j) + dky1*pflux1(is,j)
           eflux_out = dky0*eflux0(is,j) + dky1*eflux1(is,j)
           mflux_tor_out = dky0*stress_tor0(is,j) + dky1*stress_tor1(is,j)
           mflux_par_out = 0.0 !dky0*stress_par0(is,j) + dky1*stress_par1(is,j)
           exch_out = dky0*exch0(is,j) + dky1*exch1(is,j)
           write(33,*)ky,pflux_out,eflux_out,mflux_tor_out,mflux_par_out,exch_out
           pflux0(is,j) = pflux1(is,j)
           eflux0(is,j) = eflux1(is,j)
           stress_par0(is,j) = stress_par1(is,j)
           stress_tor0(is,j) = stress_tor1(is,j)
           exch0(is,j) = exch1(is,j)
           ky0 = ky1
        enddo  ! i
     enddo  ! j
  enddo  ! is 
  !
  CLOSE(33)
  !    
END SUBROUTINE write_qlgyro_flux_spectrum
!-----------------------------------------------------------------


SUBROUTINE write_qlgyro_field_spectrum
  !
  USE qlgyro_globals
  USE tglf_interface

  IMPLICIT NONE
  CHARACTER(29) :: fluxfile="out.qlgyro.field_spectrum"
  INTEGER :: i
  REAL :: phi,a_par,b_par
  CHARACTER (122) :: full_fluxfile

  full_fluxfile = trim(path)//trim(iter_path)//fluxfile

  OPEN(unit=33,file=full_fluxfile,status='replace')
!
  write(33,*)"gyro-bohm normalized field fluctuation intensity spectra per mode"
  write(33,*)"potential,A_par,B_par"
  write(33,*)"index limits: nky"
  
  write(33,*)n_ky

  if(tglf_use_bper_in)then
     write(33,*)"a_par_yes"
  else
     write(33,*)"a_par_no"
  endif
  if(tglf_use_bpar_in) then
     write(33,*)"b_par_yes"
  else
     write(33,*)"b_par_no"
  endif

  do i=1,n_ky
    phi = 0.0
    a_par=0.0
    b_par=0.0
    phi = phi + tglf_field_spectrum_out(2,i,1)
    if(tglf_use_bper_in)a_par = a_par + tglf_field_spectrum_out(3,i,1)
    if(tglf_use_bpar_in)b_par = b_par + tglf_field_spectrum_out(4,i,1)
    write(33,*)phi,a_par,b_par
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_qlgyro_field_spectrum

SUBROUTINE write_qlgyro_ky_spectrum
!
  USE qlgyro_globals
  use tglf_interface
!
  IMPLICIT NONE
  CHARACTER(22) :: fluxfile="out.qlgyro.ky_spectrum"
  INTEGER :: i
  CHARACTER (115) :: full_fluxfile

  full_fluxfile = trim(path)//trim(iter_path)//fluxfile

  OPEN(unit=33,file=full_fluxfile,status='replace')
!
  write(33,*)"index limits: nky"
  write(33,*)tglf_nky_in
!
  do i=1,tglf_nky_in
    write(33,*)tglf_ky_spectrum_out(i)
  enddo
  CLOSE(33)
!
 END SUBROUTINE write_qlgyro_ky_spectrum


SUBROUTINE write_qlgyro_eigenvalue_spectrum
  !
  USE qlgyro_globals
  use tglf_interface
  
  IMPLICIT NONE
  CHARACTER(30) :: fluxfile="out.qlgyro.eigenvalue_spectrum"
  INTEGER :: i,n
  CHARACTER (123) :: full_fluxfile

  full_fluxfile = trim(path)//trim(iter_path)//fluxfile

  !
  OPEN(unit=33,file=full_fluxfile,status='replace')
!
  write(33,*)"gyro-bohm normalized eigenvalue spectra"
  write(33,*)"(gamma(n),freq(n),n=1,nmodes_in)"
  do i=1,n_ky
    write(33,*)(tglf_eigenvalue_spectrum_out(1,i,n),tglf_eigenvalue_spectrum_out(2,i,n),n=1,tglf_nmodes_in)
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_qlgyro_eigenvalue_spectrum

 SUBROUTINE write_qlgyro_QL_weight_spectrum
   !

   use qlgyro_globals
   use tglf_interface
  IMPLICIT NONE
  CHARACTER(35) :: fluxfile="out.qlgyro.QL_weight_spectrum"
  INTEGER :: i,j,k,is
  REAL :: phinorm
  REAL,PARAMETER :: small=1.0E-10
  CHARACTER (128) :: full_fluxfile

  full_fluxfile = trim(path)//trim(iter_path)//fluxfile

  OPEN(unit=33,file=full_fluxfile,status='replace')
!
  write(33,*)"QL weights per mode:"
  write(33,*)"type: 1=particle,2=energy,3=toroidal stress,4=parallel stress,5=exchange"
  write(33,*)"index limits: nky,nmodes,ns,field,type"
  write(33,*)n_ky,n_modes,n_species,3,5
!
      ! renormalize the fluxes and intensities to the phi-norm from the v-norm
      do j=1,n_ky
         do i=1,n_modes
            phinorm=1.0
            if(ABS(tglf_field_spectrum_out(2,j,i)).gt.small)phinorm=tglf_field_spectrum_out(2,j,i)
            do is=1,n_species
               do k=1,3
                  write(33,*)tglf_flux_spectrum_out(1,is,k,j,i)/phinorm
                  write(33,*)tglf_flux_spectrum_out(2,is,k,j,i)/phinorm
                  write(33,*)tglf_flux_spectrum_out(3,is,k,j,i)/phinorm
                  write(33,*)tglf_flux_spectrum_out(4,is,k,j,i)/phinorm
                  write(33,*)tglf_flux_spectrum_out(5,is,k,j,i)/phinorm
              enddo
           enddo
        enddo
      enddo
  CLOSE(33)
!
 END SUBROUTINE write_qlgyro_QL_weight_spectrum

 subroutine write_qlgyro_gbflux

   use qlgyro_globals
   use tglf_interface
   implicit none

   CHARACTER (17) :: fluxfile="out.qlgyro.gbflux"
   CHARACTER (110) :: full_fluxfile

   full_fluxfile = trim(path)//trim(iter_path)//fluxfile

   open(unit=1,file=full_fluxfile,status='replace')
   write(1,'(32(1pe11.4,1x))') tglf_elec_pflux_out,tglf_ion_pflux_out(1:n_ion),&
        tglf_elec_eflux_out,tglf_ion_eflux_out(1:n_ion),&
        tglf_elec_mflux_out,tglf_ion_mflux_out(1:n_ion),&
        tglf_elec_expwd_out,tglf_ion_expwd_out(1:n_ion)
   close(1)

 end subroutine write_qlgyro_gbflux

 subroutine write_qlgyro_units

   use mpi
   use qlgyro_globals
   implicit none

   integer :: io=21
   character(16) :: datafile='out.qlgyro.units'
   real :: kt
   CHARACTER (109) :: full_datafile

   full_datafile = trim(path)//trim(iter_path)//datafile

   open(unit=io,file=full_datafile,status='replace')
   
   ! kT in MJ (note the conversion 1.6022e-22 MJ/keV)
   kt = 1.6022e-22*temp_norm
   
   write(io,26) 2.0*kg_proton,'m_ref (kg)'
   write(io,26) bunit,'b_unit (Tesla)'
   write(io,26) a_meters,'a (m)'
   write(io,26) csda_norm,'csD/a (1/s)'
   write(io,26) csda_norm*a_meters,'csD (m/s)'
   write(io,26) temp_norm,'Te (keV)'
   write(io,26) dens_norm,'ne (10^19/m^3)'
   write(io,26) rho_star_norm*a_meters,'rho_sD (m)'
   write(io,26) csda_norm*(rho_star_norm*a_meters)**2,&
        'chi_gBD (m^2/s)'

   write(io,26) 1e19*dens_norm*(csda_norm*a_meters)*rho_star_norm**2/0.624e22,&
        'Gamma_gBD (0.624e22/m^2/s) = (MW/keV/m^2)'
   
   write(io,26) 1e19*dens_norm*(csda_norm*a_meters)*kt*rho_star_norm**2,&
        'Q_gBD (MJ/m^2/s) = (MW/m^2)'

   write(io,26) 1e19*dens_norm*a_meters*kt*rho_star_norm**2*1e6,&
        'Pi_gBD (J/m^2) = (Nm/m^2)'
   
   write(io,26) 1e19*dens_norm*csda_norm*kt*rho_star_norm**2,&
        'S_gBD (MJ/m^3/s) = (MW/m^3)'
   
   close(io)

26 format(t2,1pe15.8,2x,a) 

 end subroutine write_qlgyro_units

  SUBROUTINE write_qlgyro_sat_geo_spectrum
!
  USE qlgyro_globals
 
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.qlgyro.sat_geo_spectrum"
  INTEGER :: i
  CHARACTER (126) :: full_fluxfile

  full_fluxfile = trim(path)//trim(iter_path)//fluxfile
  !
  OPEN(unit=33,file=full_fluxfile,status='replace')
!
  write(33,*)"saturation model geometry factor 1/(<phi|(B/(B_unit grad_r))**2\phi>/<phi|phi>) per mode:"
  write(33,*)"index limits: nky"
  write(33,*)n_ky
!
  do i=1,n_ky
    write(33,*)(sat_geo_spectrum(i))
  enddo
  CLOSE(33)
!
 END SUBROUTINE write_qlgyro_sat_geo_spectrum
!-----------------------------------------------------------------
!

   SUBROUTINE write_qlgyro_kxrms_spectrum
!
  USE qlgyro_globals
 
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.qlgyro.kxrms_spectrum"
  INTEGER :: i
  CHARACTER (126) :: full_fluxfile

  full_fluxfile = trim(path)//trim(iter_path)//fluxfile
  !
  OPEN(unit=33,file=full_fluxfile,status='replace')
!
  write(33,*)"kx rms model used"
  write(33,*)"index limits: nky"
  write(33,*)n_ky
!
  do i=1,n_ky
    write(33,*)(kxrms_spectrum(i))
  enddo
  CLOSE(33)
!
 END SUBROUTINE write_qlgyro_kxrms_spectrum
!-----------------------------------------------------------------
!

 subroutine init_qlgyro_status
   ! Reads in status file for a given ky
   
   use mpi
   use qlgyro_globals
   use tglf_interface
   implicit none
   
   integer :: file_handle, i, colour
   character(len=55) :: outstr
   CHARACTER(18) :: statusfile="out.qlgyro.status"
   CHARACTER(98) :: full_statusfile

   full_statusfile = trim(path)//statusfile

   call MPI_FILE_DELETE(full_statusfile, MPI_INFO_NULL, ierr)

   call MPI_FILE_OPEN(QLGYRO_COMM_WORLD, full_statusfile, &
        MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, file_handle, ierr)
   
   do i=1,n_kypx0

      colour = -1
      if (i .le. n_parallel) colour = i-1
      write(outstr, 44) "KY = ", tglf_ky_spectrum_out(i_ky(i)), ", PX0 = ",&
      px0_spectrum(i_px0(i)), ": fresh       ", 0, colour, NEW_LINE(' ')
      
      if (i_proc_global .eq. 0) then
         call MPI_FILE_WRITE(file_handle, TRIM(outstr), LEN(TRIM(outstr)), &
              MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      end if
   end do
   
   statstr = LEN(TRIM(outstr))
   call MPI_file_close(file_handle, ierr)
   
   
44 format(A5, F7.3, A8, F7.3, A14, 2I4, A1)
   
 end subroutine init_qlgyro_status

 
 subroutine write_qlgyro_status(runs, colors, i_kypx0, status, color_in)
   ! Writes out status file for each ky value
   ! ky_run :: integer(n_kypx0) status for each ky run
   ! i_kypx0   :: indice where to overwrite value
   ! status :: overwrite value
   ! color  :: color of cores that ran this ky

   
   use qlgyro_globals
   use tglf_interface
   use mpi
   implicit none
   
   CHARACTER(18) :: statusfile="out.qlgyro.status"
   integer, dimension(n_kypx0), intent(inout) :: runs
   integer, dimension(n_kypx0), intent(inout) :: colors
   integer, intent(in) :: i_kypx0, status, color_in
   integer :: file_handle
   integer(kind=MPI_OFFSET_KIND) :: offset
   CHARACTER(statstr) :: statusline
   CHARACTER(98) :: full_statusfile
   CHARACTER(14) :: ky_status

   full_statusfile = trim(path)//statusfile

   offset = (i_kypx0 - 1) * statstr

   runs(i_kypx0) = status
   colors(i_kypx0) = color_in

   select case(runs(i_kypx0))

   case (0)
      ky_status = ": fresh       "
   case (1)
      ky_status = ": running     "
   case (2)
      ky_status = ": unconverged "
   case (3)
      ky_status = ": converged   "
   end select

   write(statusline, 44) "KY = ", tglf_ky_spectrum_out(i_ky(i_kypx0)), ", PX0 = ",&
   px0_spectrum(i_px0(i_kypx0)), ky_status, runs(i_kypx0), colors(i_kypx0), NEW_LINE(' ')

   call MPI_FILE_OPEN(qlgyro_comm, full_statusfile, &
        MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)

   call MPI_FILE_WRITE_AT(file_handle, offset, statusline, &
           statstr, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, ierr)

   call MPI_file_close(file_handle, ierr)
   
44 format(A5, F7.3, A8, F7.3, A14, 2I4, A1)
    
 end subroutine write_qlgyro_status

 
 subroutine get_qlgyro_status(i_kypx0, status, loop_cycle)
   
   use mpi
   use qlgyro_globals
   use tglf_interface
   implicit none
      
   integer :: file_handle
   integer(kind=MPI_OFFSET_KIND) :: offset
   integer, intent(in) :: i_kypx0
   integer, intent(out) :: status
   logical, intent(out) :: loop_cycle
   CHARACTER(18) :: statusfile="out.qlgyro.status"
   CHARACTER(1) :: outstr
   CHARACTER(98) :: full_statusfile

   full_statusfile = trim(path)//statusfile

   offset = (i_kypx0-1) * statstr + 44
   call MPI_FILE_OPEN(qlgyro_comm, full_statusfile, &
        MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, ierr)

   ! Location of run status
   call MPI_FILE_READ_AT(file_handle, offset, outstr, &
        1, MPI_CHARACTER, &
        MPI_STATUS_IGNORE, ierr)

   read(outstr, *) status
   call MPI_file_close(file_handle, ierr)

   loop_cycle = .true.
   if (status .eq. 0) loop_cycle = .false.

 end subroutine get_qlgyro_status

 
 subroutine read_qlgyro_status(runs, cores)
   ! Reads in status file for a given ky
   ! i_kypx0 :: indice for ky run
   ! status :: integer stating run status
   
   use qlgyro_globals
   use tglf_interface
   use mpi
   implicit none

   integer, dimension(n_kypx0), intent(out) :: runs, cores
   integer :: file_handle
   CHARACTER(18) :: statusfile="out.qlgyro.status"
   CHARACTER(1) :: status_str
   CHARACTER(3) :: core_str
   INTEGER :: i_kypx0
   INTEGER(kind=MPI_OFFSET_KIND) :: offset1, offset2
   CHARACTER(98) :: full_statusfile

   full_statusfile = trim(path)//statusfile

   call MPI_FILE_OPEN(qlgyro_comm, full_statusfile, &
        MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, ierr)

   do i_kypx0=1,n_kypx0
      
      ! Location of run status
      offset1 = (i_kypx0-1) * statstr + 44
      call MPI_FILE_READ_AT(file_handle, offset1, status_str, &
           1, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, ierr)
      
      read(status_str, *) runs(i_kypx0)

      ! Location of run color
      offset2 = offset1 + 2
      call MPI_FILE_READ_AT(file_handle, offset2, core_str, &
           3, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, ierr)
      
      read(core_str, *) cores(i_kypx0)
      
   end do
   
   call MPI_file_close(file_handle, ierr)

end subroutine read_qlgyro_status

 ! SUBROUTINE write_tglf_density_spectrum
!   !
!   USE tglf_dimensions
!   USE tglf_global
!   USE tglf_species
!   USE tglf_kyspectrum
!   !   
!   IMPLICIT NONE
!   CHARACTER(25) :: fluxfile="out.tglf.density_spectrum"
!   INTEGER :: i,is,n
!   REAL :: density(nsm)
!   !
!   if(new_start)then
!      write(*,*)"error: tglf_TM must be called before write_tglf_density_spectrum"
!      write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
!   endif
!   !
!   OPEN(unit=33,file=fluxfile,status='replace')
! !
!   write(33,*)"gyro-bohm normalized density fluctuation amplitude spectra"
!   write(33,*)"ky,(density(is),is=1,ns_in)"
!   do i=1,nky
!     do is=1,ns
!       density(is) = 0.0
!       do n = 1,nmodes_in
!         density(is) = density(is) + intensity_spectrum_out(1,is,i,n)
!       enddo
!       density(is) = SQRT(density(is))
!     enddo
!     write(33,*)ky_spectrum(i),(density(is),is=1,ns)
!   enddo
! !
!   CLOSE(33)
! !
! END SUBROUTINE write_tglf_density_spectrum
! !-----------------------------------------------------------------

! SUBROUTINE write_tglf_temperature_spectrum
!   !
!   USE tglf_dimensions
!   USE tglf_global
!   USE tglf_species
!   USE tglf_kyspectrum
!   !   
!   IMPLICIT NONE
!   CHARACTER(29) :: fluxfile="out.tglf.temperature_spectrum"
!   INTEGER :: i,is,n
!   REAL :: temp(nsm)
!   !
!   if(new_start)then
!      write(*,*)"error: tglf_TM must be called before write_tglf_temperature_spectrum"
!      write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
!   endif
!   !
!   OPEN(unit=33,file=fluxfile,status='replace')
! !
!   write(33,*)"gyro-bohm normalized temperature fluctuation amplitude spectra"
!   write(33,*)"ky,(temperature(is),is=1,ns_in)"
!   do i=1,nky
!     do is=1,ns
!       temp(is) = 0.0
!       do n = 1,nmodes_in
!         temp(is) = temp(is) + intensity_spectrum_out(2,is,i,n)
!       enddo
!       temp(is) = SQRT(temp(is))
!     enddo
!     write(33,*)ky_spectrum(i),(temp(is),is=1,ns)
!   enddo
! !
!   CLOSE(33)
! !
! END SUBROUTINE write_tglf_temperature_spectrum
! !-----------------------------------------------------------------

! SUBROUTINE write_tglf_intensity_spectrum
!   !
!   USE tglf_dimensions
!   USE tglf_global
!   USE tglf_species
!   USE tglf_kyspectrum
!   !   
!   IMPLICIT NONE
!   CHARACTER(27) :: fluxfile="out.tglf.intensity_spectrum"
!   INTEGER :: i,is,n
!   REAL :: den,tem, u_par, q_tot
!   !
!   if(new_start)then
!      write(*,*)"error: tglf_TM must be called before write_tglf_intensity_spectrum"
!      write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
!   endif
!   !
!   OPEN(unit=33,file=fluxfile,status='replace')
! !
!   write(33,*)"gyro-bohm normalized intensity fluctuation amplitude spectra"
!   write(33,*)"ky,density,temperature,parallel velocity, parallel energy"
!   do is=ns0,ns
!     write(33,*)"species ",is
!     do i=1,nky
!       den = 0.0
!       tem = 0.0
!       u_par = 0.0
!       q_tot = 0.0
!       do n = 1,nmodes_in
!         den = den + intensity_spectrum_out(1,is,i,n)
!         tem = tem + intensity_spectrum_out(2,is,i,n)
!         u_par = u_par + intensity_spectrum_out(3,is,i,n)
!         q_tot = q_tot + intensity_spectrum_out(4,is,i,n)
!       enddo
!       write(33,*)ky_spectrum(i),den,tem,u_par,q_tot
!     enddo
!   enddo
! !
!   CLOSE(33)
! !
! END SUBROUTINE write_tglf_intensity_spectrum

! !-----------------------------------------------------------------

! SUBROUTINE write_tglf_intensity_spectrum_per_mode
!   !
!   USE tglf_dimensions
!   USE tglf_global
!   USE tglf_species
!   USE tglf_kyspectrum
!   !   
!   IMPLICIT NONE
!   CHARACTER(36) :: fluxfile="out.tglf.intensity_spectrum_per_mode"
!   INTEGER :: i,is,n
!   REAL :: den,tem, u_par, q_tot
!   !
!   if(new_start)then
!      write(*,*)"error: tglf_TM must be called before write_tglf_intensity_spectrum_per_mode"
!      write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
!   endif
!   !
!   OPEN(unit=33,file=fluxfile,status='replace')
! !
!   write(33,*)"gyro-bohm normalized intensity fluctuation amplitude spectra per mode"
!   write(33,*)"density,temperature,parallel velocity,parallel energy"
!   write(33,*)"index limits: ns,nky,nmodes"
!   write(33,*)ns,nky,nmodes_in
!   do is=ns0,ns
!     do i=1,nky
!       den = 0.0
!       tem = 0.0
!       u_par = 0.0
!       q_tot = 0.0
!       do n = 1,nmodes_in
!         den = intensity_spectrum_out(1,is,i,n)
!         tem = intensity_spectrum_out(2,is,i,n)
!         u_par = intensity_spectrum_out(3,is,i,n)
!         q_tot = intensity_spectrum_out(4,is,i,n)
!         write(33,*) den
!         write(33,*) tem
!         write(33,*) u_par
!         write(33,*) q_tot
!       enddo      
!     enddo
!   enddo
! !
!   CLOSE(33)
! !
! END SUBROUTINE write_tglf_intensity_spectrum_per_mode

! !-----------------------------------------------------------------

! !-----------------------------------------------------------------

! SUBROUTINE write_tglf_field_spectrum_per_mode
!   !
!   USE tglf_dimensions
!   USE tglf_global
!   USE tglf_species
!   USE tglf_kyspectrum
!   !   
!   IMPLICIT NONE
!   CHARACTER(36) :: fluxfile="out.tglf.field_spectrum_per_mode"
!   INTEGER :: i,n
!   REAL :: phi,a_par,b_par
!   !
!   if(new_start)then
!      write(*,*)"error: tglf_TM must be called before write_tglf_field_spectrum_per_mode"
!      write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
!   endif
!   !
!   OPEN(unit=33,file=fluxfile,status='replace')
! !
!   write(33,*)"gyro-bohm normalized field fluctuation intensity spectra per mode"
!   write(33,*)"potential,A_par,B_par"
!   write(33,*)"index limits: nky,nmodes"
!   write(33,*)nky,nmodes_in
!   if(use_bper_in)then
!      write(33,*)"a_par_yes"
!   else
!      write(33,*)"a_par_no"
!   endif
!   if(use_bpar_in) then
!      write(33,*)"b_par_yes"
!   else
!      write(33,*)"b_par_no"
!   endif
!   do i=1,nky
!     phi = 0.0
!     a_par=0.0
!     b_par=0.0
!     do n = 1,nmodes_in
!       phi = field_spectrum_out(2,i,n)
!       if(use_bper_in)a_par = field_spectrum_out(3,i,n)
!       if(use_bpar_in)b_par = field_spectrum_out(4,i,n)
!       write(33,*)phi
!       if(use_bper_in)write(33,*)a_par
!       if(use_bpar_in)write(33,*)b_par
!     enddo
!   enddo
! !
!   CLOSE(33)
! !
! END SUBROUTINE write_tglf_field_spectrum_per_mode
! !-----------------------------------------------------------------


