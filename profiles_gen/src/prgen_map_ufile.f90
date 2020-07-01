!------------------------------------------------------------
! prgen_map_ufile.f90
!
! PURPOSE:
!  Map UFILE data onto input.profiles standard.  
!  
! NOTES:
!  - See UFILE documentation at
!    http://tokamak-profiledb.ccfe.ac.uk/DOCS/PR08MAN/pdbman.html
!------------------------------------------------------------

subroutine prgen_map_ufile

  use prgen_globals
  use expro

  implicit none

  integer :: i

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  expro_n_exp = nx
  expro_n_ion = ufile_nion
  call expro_init(1)
  !
  expro_rho(:)   = rho(:)
  expro_rmin(:)  = rmin(:)
  expro_rmaj(:)  = rmaj(:)
  expro_te(:)    = ufile_te(:)*1e-3
  expro_ne(:)    = ufile_ne(:)*1e-19
  expro_z_eff(:) = ufile_zeff(:)
  expro_w0(:)    = 0.0
  expro_ptot(:)  = ufile_pres(:)

  print '(a,i2)','INFO: (prgen) Number of ions: ',ufile_nion

  !-----------------------------------------------------------------
  ! Construct ion densities and temperatures
  !
  do i=1,ufile_nion
     expro_ni(i,:) = ufile_ni(:,i)*1e-19
     expro_ti(i,:) = ufile_ti(:,i)*1e-3
  enddo

  ! vphi
  expro_vtor(:,:) = 0.0

  ! Insert carbon toroidal velocity
  do i=1,ufile_nion
     if (nint(ufile_m(i)) == 12) then
        print '(a)', 'INFO: (prgen) Assuming VROT is the carbon toroidal rotation.'
        expro_vtor(i,:) = -ufile_vrot(:)*(rmaj(:)+rmin(:))
     endif
  enddo

  ! vpol
  expro_vpol(:,:) = 0.0

  ! Heating powers
  
  expro_qohme = ufile_qohm
  
  expro_qbeami = ufile_qnbii
  expro_qrfi   = ufile_qicrhi+ufile_qechi+ufile_qlhi

  expro_qbeame = ufile_qnbie
  expro_qrfe   = ufile_qicrhe+ufile_qeche+ufile_qlhe

  expro_qfuse = ufile_qfuse
  expro_qfusi = ufile_qfusi

  expro_qei = ufile_qei
  expro_qbeami = -ufile_qwalli
  
  print '(a)','INFO: (prgen) i-power: Setting expro_qbeami = -QWALLI'
  print '(a)','INFO: (prgen) i-power: (QNBII+QICHRI+QLHI+QECHI)+QEI-QWALLI'

  expro_qbeame = -ufile_qwalle
 
  ! JC: need radiated power breakdown
  !  ? = ufile_qrad(:)
  print '(a)','INFO: (prgen) e-power: Setting expro_qbeame = -QWALLE'
  print '(a)','INFO: (prgen) e-power: Setting QRAD=0'
  print '(a)','INFO: (prgen) e-power: (QNBIE+QICHRE+QLHE+QECHE+QOHM)-QE-QRAD-QWALLE'

  ! No particle or momentum source
  expro_qpar_beam = 0.0
  expro_qpar_wall = 0.0
  expro_qmom = 0.0

end subroutine prgen_map_ufile
