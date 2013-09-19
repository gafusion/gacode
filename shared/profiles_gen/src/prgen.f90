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
  use EXPRO_interface

  !--------------------------------------------------
  implicit none
  !--------------------------------------------------

  !--------------------------------------------------
  ! Parse the config file
  !
  open(unit=1,file='.config',status='old')
  read(1,'(a)') date
  read(1,'(a)') raw_data_file
  read(1,'(a)') raw_data_type
  read(1,'(a)') cer_file
  read(1,*) efit_method
  read(1,*) nogatoq_flag
  read(1,*) verbose_flag
  read(1,*) gmerge_flag
  read(1,*) ipccw
  read(1,*) btccw
  read(1,*) nfourier
  read(1,*) reorder_vec(:)
  close(1)
  !--------------------------------------------------

  n_indx  = size(EXPRO_tag) 

  !------------------------------------------------------------------
  ! Read the iterdb file and define standard variables.
  !
  ! Note that nx will be the experimental vector length in ALL cases:
  !
  if (gmerge_flag == 1) then

     call prgen_read_inputprofiles

     format_type = 7

  else if (trim(raw_data_type) == 'null') then

     ! Pure gfile parsing

     format_type = 0

     ! Minimal processing required to merge gfile data into otherwise
     ! empty input.profiles output file.
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

  else if (trim(raw_data_type) == 'ITERDB') then

     ! Old text format
     print '(a)','INFO: (prgen) Assuming old iterdb text format.'

     format_type = 1

     call prgen_read_iterdb 

  endif
  !------------------------------------------------------------------

  !---------------------------------------------------
  ! Read the GATO file for "better" geometry.  At this
  ! point, GATO has already run and we are just reading 
  ! the output.
  !
  select case (efit_method)
  case (1)
     ! Use geometry data contained in profile data 
  case (2)
     ! Use GATO-EFIT mapper
     call prgen_read_gato
     print '(a)','INFO: (prgen) Wrote input.profiles.geo.'
  case (3)
     ! Use OMFIT-EFIT mapper
     call prgen_read_omfit
     print '(a)','INFO: (prgen) Wrote input.profiles.geo'
  case (4,5)
     ! Use DSKGATO data
     call prgen_read_dskgato
     print '(a)','INFO: (prgen) Wrote input.profiles.geo'
  end select
  !---------------------------------------------------

  select case (format_type)

  case (0)
     call prgen_map_null
  case (1)
     call prgen_map_iterdb
  case (2) 
     call prgen_map_plasmastate
  case (3) 
     call prgen_map_peqdsk
  case (5) 
     call prgen_map_corsica
  case (6)
     call prgen_map_ufile
  case (7)
     call prgen_map_inputprofiles

  end select

  call prgen_write

  ! Successful completion
  open(unit=1,file='success',status='replace')
  write(1,*) 1
  close(1)

end program prgen
