!-------------------------------------------------------------------
! prgen.f90
!
! PURPOSE:
!  Generate input.profiles from prexisting data.  This executable is
!  managed and called by gacode/shared/bin/profiles_gen.
!
!  Autodetected profile formats:
!   GACODE      (input.profiles)
!   ITERDB      (text iterdb)
!   ITERDBNC    (netCDF iterdb)
!   SWIM        (plasmastate)
!   PFILE       (peqdsk)
!   GENF        (General Fusion)
!   MANUAL      (Auto-import of text files)
!
!  Autodetected geometry formats:
!   GFILE       (geqdsk equilibrium data)
!--------------------------------------------------------------------

program prgen

  use prgen_globals

  implicit none
  
  integer :: i

  ! Default reordering vector
  do i=1,n_ion_max
     reorder_vec(i) = i
  enddo

  !--------------------------------------------------------------------
  ! Parse the config file
  !
  open(unit=1,file='.prgenconfig',status='old')
  read(1,'(a)') date
  read(1,'(a)') file_state
  read(1,'(a)') raw_data_type
  read(1,'(a)') file_g
  read(1,'(a)') file_cer
  read(1,'(a)') file_ti
  read(1,*) efit_method
  read(1,*) verbose_flag
  read(1,*) ipccw
  read(1,*) btccw
  read(1,*) nfourier
  read(1,*) lump_fast_flag
  read(1,*) true_aux_flag
  read(1,*) reorder_vec(1:10)
  read(1,*) n_lump
  allocate(lump_vec(n_lump))
  read(1,*) lump_vec(:)
  close(1)
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Read the iterdb file and define standard variables.
  !
  ! Note that nx will be the experimental vector length in ALL cases:
  !
  select case (trim(raw_data_type))

  case('null') ! gfile only
     print '(a)','INFO: (prgen) Merging gfile data into otherwise empty input.gacode.'
     format_type = 0
     ! Nothing to read

  case ('ITERDB') ! Old text format
     print '(a)','INFO: (prgen) Assuming old iterdb text format.'

     format_type = 1
     call prgen_read_iterdb 

  case ('ITERDBNC') ! New NetCDF format
     print '(a)','INFO: (prgen) Assuming iterdb NetCDF format.'

     format_type = 1
     call prgen_read_iterdb_nc

  case ('SWIM') ! Plasmastate format
     print '(a)','INFO: (prgen) Assuming SWIM (plasmastate) format.'

     format_type = 2
     call prgen_read_plasmastate

  case ('PFILE') ! peqdsk format
     print '(a)','INFO: (prgen) Assuming PFILE (peqdsk) format.'

     format_type = 3
     call prgen_read_peqdsk

  case ('GENF') ! General Fusion format
     print '(a)','INFO: (prgen) Assuming GENF (General Fusion) format.'

     format_type = 4
     call prgen_read_genf

  case ('manual') ! MANUAL import
     print '(a)','INFO: (prgen) Assuming MANUAL file import.'

     format_type = 5
     call prgen_read_manual

  case('GACODE') ! Note (we may or may not have gmerge_flag == 1)
     print '(a)','INFO: (prgen) Assuming input.gacode (GACODE) format.'

     format_type = 7
     call prgen_read_inputgacode

  case('LEGACY') ! Note (we may or may not have gmerge_flag == 1)
     print '(a)','INFO: (prgen) Assuming input.profiles (LEGACY GACODE) format.'

     format_type = 8
     call prgen_read_inputprofiles

  case default ! No case matched

     print '(a)','ERROR: (prgen) Unmatched options.'
     stop

  end select
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Contour the EFIT data for "better" geometry analysis and calculation
  ! of shape parameters.  At this point, the shape coefficients have
  ! already been computed and we are just reading the output.
  !
  select case (efit_method)
  case (0)
     ! Use geometry data contained in profile data 
     print '(a)','WARNING: (prgen) Using original geometry data. Better to use gfile.'
  case (1)
     ! Use EFIT mapper
     call prgen_geometry
  end select
  !--------------------------------------------------------------------

  !-----------------------------------------------------
  ! Set ipccw and btccw to standard DIII-D configuration
  ! if not defined at input
  if (ipccw == 0) ipccw = 1
  if (btccw == 0) btccw = -1
  !------------------------------------------------------

  select case (format_type)

  case (0)
     call prgen_map_null
  case (1)
     call prgen_map_iterdb
  case (2) 
     call prgen_map_plasmastate
  case (3) 
     call prgen_map_peqdsk
  case (4)
     call prgen_map_genf
  case (5)
     call prgen_map_manual
  case (7,8)
     call prgen_map_inputgacode
  end select

  call prgen_swap
  call prgen_write

end program prgen
