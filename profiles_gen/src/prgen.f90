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
!   CORSICA     (Corsica) 
!   UFILE       (ITPA profile database format)
!
!  Autodetected geometry formats:
!   GFILE       (geqdsk equilibrium data)
!   DSKGATO_OLD (old-type dskgato flux-surface data)
!   DSKGATO_NEW (new-type dskgato flux-surface data)
!--------------------------------------------------------------------

program prgen

  use prgen_globals

  implicit none

  !--------------------------------------------------------------------
  ! Parse the config file
  !
  open(unit=1,file='.prgenconfig',status='old')
  read(1,'(a)') date
  read(1,'(a)') file_state
  read(1,'(a)') raw_data_type
  read(1,'(a)') file_g
  read(1,'(a)') file_cer
  read(1,*) efit_method
  read(1,*) verbose_flag
  read(1,*) ipccw
  read(1,*) btccw
  read(1,*) nfourier
  read(1,*) n_null
  read(1,*) lump_fast_flag
  read(1,*) true_aux_flag
  read(1,*) reorder_vec(:)
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
  if (trim(raw_data_type) == 'GACODE') then

     ! Note (we may or may not have gmerge_flag == 1)
     print '(a)','INFO: (prgen) Assuming input.gacode (GACODE) format.'

     call prgen_read_inputgacode

     format_type = 7

  else if (trim(raw_data_type) == 'LEGACY') then

     ! Note (we may or may not have gmerge_flag == 1)
     print '(a)','INFO: (prgen) Assuming input.profiles (LEGACY GACODE) format.'

     call prgen_read_inputprofiles

     format_type = 8

  else if (trim(raw_data_type) == 'null') then

     ! Pure gfile parsing

     format_type = 0

     ! Minimal processing required to merge gfile data into otherwise empty input.gacode 
     call prgen_read_null

  else if (trim(raw_data_type) == 'ITERDBNC') then

     ! New NetCDF format
     print '(a)','INFO: (prgen) Assuming iterdb NetCDF format.'

     format_type = 1

     call prgen_read_iterdb_nc

  else if (trim(raw_data_type) == 'SWIM') then

     ! Plasmastate format
     print '(a)','INFO: (prgen) Assuming SWIM (plasmastate) format.'

     format_type = 2

     call prgen_read_plasmastate

  else if (trim(raw_data_type) == 'PFILE') then

     ! peqdsk format
     print '(a)','INFO: (prgen) Assuming PFILE (peqdsk) format.'

     format_type = 3

     call prgen_read_peqdsk

  else if (trim(raw_data_type) == 'CORSICA') then

     ! corsica format
     print '(a)','INFO: (prgen) Assuming CORSICA format.'

     format_type = 5

     call prgen_read_corsica

  else if (trim(raw_data_type) == 'UFILE') then

     ! UFILE format
     print '(a)','INFO: (prgen) Assuming UFILE format.'

     format_type = 6

     call prgen_read_ufile

  else if (trim(raw_data_type) == 'GENF') then

     ! General Fusion format
     print '(a)','INFO: (prgen) Assuming GENF (General Fusion) format.'

     format_type = 4

     call prgen_read_genf

  else if (trim(raw_data_type) == 'ITERDB') then

     ! Old text format
     print '(a)','INFO: (prgen) Assuming old iterdb text format.'

     format_type = 1

     call prgen_read_iterdb 

  else

     ! No case matched

     print '(a)','ERROR: (prgen) Unmatched options.'
     stop

  endif
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Contour the EFIT data for "better" geometry analysis and calculation
  ! of shape parameters.  At this point, the shape coefficients have
  ! already been computed and we are just reading the output.
  !
  select case (efit_method)
  case (0)
     ! Use geometry data contained in profile data 
     print '(a)','INFO: (prgen) Using original geometry data.'
  case (1)
     ! Use OMFIT-EFIT mapper
     call prgen_read_omfit
  end select
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! High-resolution geometry
  !
  if (efit_method == 1) then
     if (format_type == 1 .or. format_type == 2) then
        print '(a,2(f10.8,a))', &
             'INFO: (prgen) dpsi = ',dpsi_data,' (statefile) ',dpsi_efit,' (new)'
     endif
  endif
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
     call prgen_map_corsica
  case (6)
     call prgen_map_ufile
  case (7,8)
     call prgen_map_inputgacode
  end select

  call prgen_swap
  call prgen_write

  ! Successful completion
  open(unit=1,file='success',status='replace')
  write(1,*) 1
  close(1)

end program prgen
