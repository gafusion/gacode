!--------------------------------------------------------------
! prgen_read_iterdb.f90
!
! PURPOSE:
!  Extract iterdb data from native TEXT file.
!--------------------------------------------------------------

subroutine prgen_read_iterdb

  use prgen_globals

  implicit none

  character (len=100) :: t
  integer :: i
  real :: x
  real, dimension(:), allocatable :: xv
  real, dimension(:), allocatable :: xvv

  !----------------------------------------------------
  ! Read the iterdb file
  !
  open(unit=1,file=raw_data_file,status='old')
  read(1,*) t

  read(1,*) t ; read(1,*) onetwo_ishot
  read(1,*) t ; read(1,*) onetwo_nj
  read(1,*) t ; read(1,*) onetwo_nion
  read(1,*) t ; read(1,*) onetwo_nprim
  read(1,*) t ; read(1,*) onetwo_nimp
  read(1,*) t ; read(1,*) onetwo_nneu
  read(1,*) t ; read(1,*) i
  read(1,*) t ; read(1,*) onetwo_namep(1:onetwo_nprim)
  read(1,*) t ; read(1,*) onetwo_namei(1:onetwo_nimp)
  onetwo_nameb(1) = onetwo_namep(i)
  onetwo_nbion = 1 ! Assume 1 for now
  read(1,*) t ; read(1,*) t
  read(1,*) t ; read(1,*) onetwo_time
  read(1,*) t ; read(1,*) onetwo_Rgeom
  read(1,*) t ; read(1,*) onetwo_Rmag
  read(1,*) t ; read(1,*) onetwo_R0
  read(1,*) t ; read(1,*) onetwo_kappa
  read(1,*) t ; read(1,*) onetwo_delta
  read(1,*) t ; read(1,*) x
  read(1,*) t ; read(1,*) onetwo_volo
  read(1,*) t ; read(1,*) onetwo_cxareao
  read(1,*) t ; read(1,*) onetwo_Btor
  read(1,*) t ; read(1,*) x
  read(1,*) t ; read(1,*) x
  read(1,*) t ; read(1,*) x
  read(1,*) t ; read(1,*) x
  read(1,*) t ; read(1,*) x ! Te0
  read(1,*) t ; read(1,*) x ! Ti0

  nx = onetwo_nj

  call allocate_internals
  call allocate_iterdb_vars

  allocate(xv(onetwo_nj))

  read(1,*) t ; read(1,*) onetwo_psi ! psi on rho grid
  read(1,*) t ; read(1,*) onetwo_rho_grid
  read(1,*) t ; read(1,*) xv ! fcap
  read(1,*) t ; read(1,*) xv ! gcap
  read(1,*) t ; read(1,*) onetwo_hcap
  read(1,*) t ; read(1,*) onetwo_te
  read(1,*) t ; read(1,*) onetwo_ti
  read(1,*) t ; read(1,*) q
  read(1,*) t ; read(1,*) onetwo_ene


  do i=1,onetwo_nion
     read(1,*) t ; read(1,*) onetwo_enion(:,i)
  enddo


  sbcx_d(:) = 0.0
  sion_d(:) = 0.0
  do i=1,onetwo_nprim
     read(1,*) t ; read(1,*) onetwo_sion(:,i)
     read(1,*) t ; read(1,*) onetwo_srecom(:,i)
     read(1,*) t ; read(1,*) onetwo_scx(:,i)
     read(1,*) t ; read(1,*) onetwo_sbcx(:,i)
     read(1,*) t ; read(1,*) onetwo_s(:,i)
     read(1,*) t ; read(1,*) onetwo_dudt(:,i)
     sbcx_d(:) = sbcx_d(:)+onetwo_sbcx(:,i)
     sion_d(:) = sion_d(:)+onetwo_sion(:,i)
  enddo

  read(1,*) t ; read(1,*) onetwo_enbeam(:,1) ! fast ion density


  do i=1,onetwo_nprim
     read(1,*) t ; read(1,*) xv ! neutral density
  enddo

  do i=1,onetwo_nprim
     read(1,*) t ; read(1,*) xv ! neutral density from wall source
  enddo

  do i=1,onetwo_nprim
     read(1,*) t ; read(1,*) xv ! neutral density from volume source
  enddo

  do i=1,onetwo_nprim
     read(1,*) t ; read(1,*) xv ! volume source of neutrals
  enddo

  read(1,*) t ; read(1,*) xv ! sbion, beam electron source
  read(1,*) t ; read(1,*) onetwo_sbeam ! sbion, beam thermal ion source
  read(1,*) t ; read(1,*) xv ! total current density
  read(1,*) t ; read(1,*) xv ! ohmic current density
  read(1,*) t ; read(1,*) xv ! bootstrap current density
  read(1,*) t ; read(1,*) xv ! beam-driven current density
  read(1,*) t ; read(1,*) xv ! RF current density
  read(1,*) t ; read(1,*) xv ! rho*bp0*fcap*gcap*hcap, tesla*meters
  read(1,*) t ; read(1,*) onetwo_zeff   ! 31
  read(1,*) t ; read(1,*) onetwo_angrot ! 32
  read(1,*) t ; read(1,*) xv ! electron thermal diffusivity
  read(1,*) t ; read(1,*) xv ! ion thermal diffusivity
  read(1,*) t ; read(1,*) xv ! ion nc thermal diffusivity
  read(1,*) t ; read(1,*) onetwo_dpedt ! 36
  read(1,*) t ; read(1,*) onetwo_dpidt ! 37
  read(1,*) t ; read(1,*) xv ! electron conduction
  read(1,*) t ; read(1,*) xv ! ion conduction
  read(1,*) t ; read(1,*) xv ! electron convection
  read(1,*) t ; read(1,*) xv ! ion convection
  read(1,*) t ; read(1,*) onetwo_qbeame  ! 42
  read(1,*) t ; read(1,*) onetwo_qdelt   ! 43
  read(1,*) t ; read(1,*) onetwo_qbeami  ! 44
  read(1,*) t ; read(1,*) onetwo_qrfe    ! 45
  read(1,*) t ; read(1,*) onetwo_qrfi    ! 46
  read(1,*) t ; read(1,*) onetwo_qione   ! 47
  read(1,*) t ; read(1,*) onetwo_qioni   ! 48
  read(1,*) t ; read(1,*) onetwo_qcx     ! 49
  read(1,*) t ; read(1,*) xv ! 2d electron heating
  read(1,*) t ; read(1,*) xv ! 2d ion heating
  read(1,*) t ; read(1,*) onetwo_qfuse  ! 52
  read(1,*) t ; read(1,*) onetwo_qfusi  ! 53
  read(1,*) t ; read(1,*) xv ! beam fusion electron heating
  read(1,*) t ; read(1,*) xv ! beam fusion ion heating
  read(1,*) t ; read(1,*) xv ! qmag electron heating
  read(1,*) t ; read(1,*) xv ! sawtooth electron heating
  read(1,*) t ; read(1,*) xv ! sawtooth ion heating
  read(1,*) t ; read(1,*) onetwo_qrad ! 59
  read(1,*) t ; read(1,*) onetwo_qohm ! 60
  read(1,*) t ; read(1,*) rmaj ! 61
  read(1,*) t ; read(1,*) rmin ! 62
  read(1,*) t ; read(1,*) onetwo_volume ! 63
  read(1,*) t ; read(1,*) kappa ! 64
  read(1,*) t ; read(1,*) delta ! 65
  read(1,*) t ; read(1,*) xv ! indentation of each flux surface
  read(1,*) t
  read(1,*) t ; read(1,*) xv ! surface area each flux surface
  read(1,*) t ; read(1,*) xv ! cross-sectional area each surface
  read(1,*) t ; read(1,*) xv ! flux surface average absolute grad rho
  read(1,*) t ; read(1,*) xv ! flux surface  grad_rho_sq
  read(1,*) t ; read(1,*) onetwo_nb !number points in plasma boundary

  allocate(xvv(onetwo_nb))

  read(1,*) t ; read(1,*) xvv !plasma boundary r
  read(1,*) t ; read(1,*) xvv !plasma boundary z

  ! Torque density may be missing on iterdb file
  read(1,'(a)',iostat=i) t
  if (i == 0) then
     print '(3(a))', 'INFO: (prgen) Assuming "', trim(t), '" is beam torque density.'
     read(1,*) onetwo_storqueb !torque density nt-m/m**3
  else
     onetwo_storqueb(:) = 0.0
  endif

  ! Beam pressure may be missing on iterdb file
  read(1,'(a)',iostat=i) t
  if (i == 0) then
     print '(3(a))', 'INFO: (prgen) Assuming "', trim(t), '" is beam pressure.'
     read(1,*) onetwo_pressb(:,1)
  else
     onetwo_pressb(:,:) = 0.0
     onetwo_nbion = 0 ! Need beam pressure to get effective beam temp
  endif

  ! Total pressure may be missing on iterdb file
  read(1,'(a)',iostat=i) t
  if (i == 0) then
     print '(3(a))', 'INFO: (prgen) Assuming "', trim(t), '" is total pressure.'
    read(1,*) onetwo_press(:)
  else
    onetwo_press(:) = 0.0
  endif

  dpsi(:) = onetwo_psi(:)-onetwo_psi(1)

  ! No squareness
  zeta(:) = 0.0

  ! No elevation
  zmag(:) = 0.0

end subroutine prgen_read_iterdb
