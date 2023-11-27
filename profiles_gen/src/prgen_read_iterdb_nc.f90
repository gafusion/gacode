!--------------------------------------------------------------
! prgen_read_iterdb_nc.f90
!
! PURPOSE:
!  Extract iterdb data from native NETCDF file.
!--------------------------------------------------------------

subroutine prgen_read_iterdb_nc

  use prgen_globals
  use netcdf

  implicit none

  integer :: i

  ! NetCDF variables
  integer :: ncid
  integer :: varid
  integer :: err

  real, dimension(:), allocatable :: work

  ! Open the file (NF90_NOWRITE means read-only)
  err = nf90_open(file_state,NF90_NOWRITE,ncid)

  err = nf90_inq_varid(ncid,trim('shot'),varid)
  err = nf90_get_var(ncid,varid,onetwo_ishot)

  err = nf90_inq_varid(ncid,trim('nj'),varid)
  err = nf90_get_var(ncid,varid,nx)

  err = nf90_inq_dimid(ncid,trim('dim_npsi'),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=onetwo_npsi)

  err = nf90_inq_varid(ncid,trim('nion'),varid)
  err = nf90_get_var(ncid,varid,onetwo_nion)

  err = nf90_inq_varid(ncid,trim('nprim'),varid)
  err = nf90_get_var(ncid,varid,onetwo_nprim)

  err = nf90_inq_varid(ncid,trim('nimp'),varid)
  err = nf90_get_var(ncid,varid,onetwo_nimp)

  err = nf90_inq_varid(ncid,trim('nneu'),varid)
  err = nf90_get_var(ncid,varid,onetwo_nneu)

  err = nf90_inq_varid(ncid,trim('nbion'),varid)
  err = nf90_get_var(ncid,varid,onetwo_nbion)

  err = nf90_inq_varid(ncid,trim('namep'),varid)
  err = nf90_get_var(ncid,varid,onetwo_namep(1:onetwo_nprim))

  err = nf90_inq_varid(ncid,trim('namei'),varid)
  err = nf90_get_var(ncid,varid,onetwo_namei(1:onetwo_nimp))

  err = nf90_inq_varid(ncid,trim('nameb'),varid)
  err = nf90_get_var(ncid,varid,onetwo_nameb(1:onetwo_nbion))

  err = nf90_inq_varid(ncid,trim('time'),varid)
  err = nf90_get_var(ncid,varid,onetwo_time)

  err = nf90_inq_varid(ncid,trim('rgeom'),varid)
  err = nf90_get_var(ncid,varid,onetwo_Rgeom)

  err = nf90_inq_varid(ncid,trim('rma'),varid)
  err = nf90_get_var(ncid,varid,onetwo_Rmag)

  err = nf90_inq_varid(ncid,trim('rmajor'),varid)
  err = nf90_get_var(ncid,varid,onetwo_R0)

  err = nf90_inq_varid(ncid,trim('kappa'),varid)
  err = nf90_get_var(ncid,varid,onetwo_kappa)

  err = nf90_inq_varid(ncid,trim('deltao'),varid)
  err = nf90_get_var(ncid,varid,onetwo_delta)

  err = nf90_inq_varid(ncid,trim('volume'),varid)
  err = nf90_get_var(ncid,varid,onetwo_volo)

  err = nf90_inq_varid(ncid,trim('areao'),varid)
  err = nf90_get_var(ncid,varid,onetwo_cxareao)

  err = nf90_inq_varid(ncid,trim('btor'),varid)
  err = nf90_get_var(ncid,varid,bcentr)

  err = nf90_inq_varid(ncid,trim('tot_cur'),varid)
  err = nf90_get_var(ncid,varid,current)

  err = nf90_inq_varid(ncid,trim('Ipsign'),varid)
  if (err == 0) then
     err = nf90_get_var(ncid,varid,onetwo_ipccw)
     if (ipccw == 0) then
        ipccw = onetwo_ipccw
     endif
  endif
  
  call prgen_allocate('iterdb')

  allocate(work(onetwo_npsi))

  err = nf90_inq_varid(ncid,trim('rho_grid'),varid)
  err = nf90_get_var(ncid,varid,onetwo_rho_grid)
  ! Will use this for normalized rho-grid
  do i=1,nx
     rho(i) = (i-1)/(nx-1.0)
  enddo

  err = nf90_inq_varid(ncid,trim('hcap'),varid)
  err = nf90_get_var(ncid,varid,onetwo_hcap)

  err = nf90_inq_varid(ncid,trim('Te'),varid)
  err = nf90_get_var(ncid,varid,onetwo_te)

  err = nf90_inq_varid(ncid,trim('Ti'),varid)
  err = nf90_get_var(ncid,varid,onetwo_ti)

  err = nf90_inq_varid(ncid,trim('psir_grid'),varid)
  err = nf90_get_var(ncid,varid,onetwo_psi)

  err = nf90_inq_varid(ncid,trim('q_value'),varid)
  err = nf90_get_var(ncid,varid,q)

  err = nf90_inq_varid(ncid,trim('curohm'),varid)
  err = nf90_get_var(ncid,varid,johm)

  err = nf90_inq_varid(ncid,trim('curboot'),varid)
  err = nf90_get_var(ncid,varid,jbs)

  err = nf90_inq_varid(ncid,trim('currf'),varid)
  err = nf90_get_var(ncid,varid,jrf)

  err = nf90_inq_varid(ncid,trim('curbeam'),varid)
  err = nf90_get_var(ncid,varid,jnb)

  err = nf90_inq_varid(ncid,trim('ene'),varid)
  err = nf90_get_var(ncid,varid,onetwo_ene)

  err = nf90_inq_varid(ncid,trim('enion'),varid)
  err = nf90_get_var(ncid,varid,onetwo_enion(:,1:onetwo_nion))

  err = nf90_inq_varid(ncid,trim('enbeam'),varid)
  err = nf90_get_var(ncid,varid,onetwo_enbeam(:,1:onetwo_nbion))
  
  err = nf90_inq_varid(ncid,trim('enalp'),varid)
  err = nf90_get_var(ncid,varid,onetwo_enalp)

  ! Total plasma pressure
  err = nf90_inq_varid(ncid,trim('press'),varid)
  err = nf90_get_var(ncid,varid,p_tot)

  err = nf90_inq_varid(ncid,trim('pressb'),varid)
  err = nf90_get_var(ncid,varid,onetwo_pressb(:,1:onetwo_nbion))

!==============================================================
!    volumetric_electron_heating_terms = SortedDict(list(zip(
!        ['qohm', 'qdelt', 'qrad', 'qione', 'qbeame', 'qrfe', 'qfuse'],
!        [[1, 7], [-1, 11], [-1, 10], [-1, 602], [1, 2], [1, 3], [1, 6]]  # [sign, index_id]
!    )))
!--------------------------------------------------------------
  err = nf90_inq_varid(ncid,trim('qohm'),varid)
  err = nf90_get_var(ncid,varid,qohm) ; qohm = 1e-6*qohm ! W/m^3 -> MW/m^3

! QUESTION: qdelt flipped sign after ONETWO v5.8
!           version info is contained in netcdf global 'title' variable
!           The sign flip refers to which version?
  err = nf90_inq_varid(ncid,trim('qdelt'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qdelt)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qdelt(:) = -1*onetwo_qdelt(:)

  err = nf90_inq_varid(ncid,trim('qrad'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qrad)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qrad(:) = -1*onetwo_qrad(:)

  err = nf90_inq_varid(ncid,trim('qione'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qione)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qione(:) = -1*onetwo_qione(:)

  err = nf90_inq_varid(ncid,trim('qbeame'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qbeame)

  err = nf90_inq_varid(ncid,trim('qrfe'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qrfe)

  err = nf90_inq_varid(ncid,trim('qfuse'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qfuse)

!==============================================================
!    volumetric_ion_heating_terms = SortedDict(list(zip(
!        ['qdelt', 'qioni', 'qcx', 'qbeami', 'qrfi', 'qfusi'],
!        [[1, 11], [1, 602], [-1, 305], [1, 2], [1, 5], [1, 6]]
!    )))
!--------------------------------------------------------------
!  err = nf90_inq_varid(ncid,trim('qdelt'),varid)
!  err = nf90_get_var(ncid,varid,onetwo_qdelt)
!  !Correct for difference in sign convention with ASCII format iterdb
!  onetwo_qdelt(:) = -1*onetwo_qdelt(:)

  err = nf90_inq_varid(ncid,trim('qioni'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qioni)

  err = nf90_inq_varid(ncid,trim('qcx'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qcx)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qcx(:) = -1*onetwo_qcx(:)

  err = nf90_inq_varid(ncid,trim('qbeami'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qbeami)

  err = nf90_inq_varid(ncid,trim('qrfi'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qrfi)

  err = nf90_inq_varid(ncid,trim('qfusi'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qfusi)

!==============================================================
!    volumetric_electron_particles_terms = SortedDict(list(zip(
!        ['sion_thermal_e', 'sbeame', 'srecom_e', 'sion_imp_e', 's2d_e', 'ssaw_e', 'spellet'],
!        [[1, 601], [1, 2], [1, 602], [1, 0], [1, 0], [1, 0], [1, 14]]
!    )))
!--------------------------------------------------------------
! NOTE: sum over sion is equivalent to sion_thermal_e
! Source of electrons due to electron impact ionization of thermal neutrals
!  err = nf90_inq_varid(ncid,trim('sion_thermal_e'),varid)
!  err = nf90_get_var(ncid,varid,sion_d)
  err = nf90_inq_varid(ncid,trim('sion'),varid)
  err = nf90_get_var(ncid,varid,onetwo_sion(:,1:onetwo_nion))
  sion_d(:) = 0.0
  do i=1,onetwo_nprim
     sion_d(:) = sion_d(:)+onetwo_sion(:,i)
  enddo

! QUESTION: why should this be negative?
! Source of fast ions (and thermal neutrals) due to charge exchange of beam neutrals with thermal ions
!  err = nf90_inq_varid(ncid,trim('sbcx'),varid)
!  err = nf90_get_var(ncid,varid,onetwo_sbcx(:,1:onetwo_nion))
!  sbcx_d(:) = 0.0
!  do i=1,onetwo_nprim
!     sbcx_d(:) = sbcx_d(:)+onetwo_sbcx(:,i)
!  enddo

!  sbeame : beam electron source,#/(meter^3*second)
  err = nf90_inq_varid(ncid,trim('sbeame'),varid)
  err = nf90_get_var(ncid,varid,onetwo_sbeame)

!  sbeam : beam ion  source,#/(meter^3*second), species: d
  err = nf90_inq_varid(ncid,trim('sbeam'),varid)
  err = nf90_get_var(ncid,varid,onetwo_sbeam)

! Source of electrons due to recombination
!  err = nf90_inq_varid(ncid,trim('srecom_e'),varid)
!  err = nf90_get_var(ncid,varid,onetwo_srecom_e)

!  err = nf90_inq_varid(ncid,trim('sion_imp_e'),varid)
!  err = nf90_get_var(ncid,varid,onetwo_sion_imp_e)

!  err = nf90_inq_varid(ncid,trim('s2d_e'),varid)
!  err = nf90_get_var(ncid,varid,onetwo_s2d_e)

!  err = nf90_inq_varid(ncid,trim('ssaw_e'),varid)
!  err = nf90_get_var(ncid,varid,onetwo_ssaw_e)

!  err = nf90_inq_varid(ncid,trim('spellet'),varid)
!  err = nf90_get_var(ncid,varid,onetwo_spellet)

!==============================================================
!    volumetric_momentum_terms = SortedDict(list(zip(
!        ['storqueb', ],  # this is the total rather than individual components
!        [[1, 1]]
!    )))
!--------------------------------------------------------------
  err = nf90_inq_varid(ncid,trim('storqueb'),varid)
  err = nf90_get_var(ncid,varid,onetwo_storqueb)

!==============================================================

  err = nf90_inq_varid(ncid,trim('zeff'),varid)
  err = nf90_get_var(ncid,varid,zeff)

  err = nf90_inq_varid(ncid,trim('angrot'),varid)
  err = nf90_get_var(ncid,varid,onetwo_angrot)

  err = nf90_inq_varid(ncid,trim('dpedt'),varid)
  err = nf90_get_var(ncid,varid,onetwo_dpedt)

  err = nf90_inq_varid(ncid,trim('dpidt'),varid)
  err = nf90_get_var(ncid,varid,onetwo_dpidt)

  err = nf90_inq_varid(ncid,trim('rho_mhd_gridnpsi'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_rho_mhd_gridnpsi,onetwo_npsi)

  ! shrink flux
  onetwo_rho_mhd_gridnpsi = onetwo_rho_mhd_gridnpsi &
       * onetwo_rho_grid(nx) &
       / onetwo_rho_mhd_gridnpsi(onetwo_npsi)

  onetwo_rho_mhd_gridnpsi(onetwo_npsi) = onetwo_rho_grid(nx)

  err = nf90_inq_varid(ncid,trim('rmajavnpsi'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_rmajavnpsi,onetwo_npsi)

  err = nf90_inq_varid(ncid,trim('rminavnpsi'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_rminavnpsi,onetwo_npsi)

  err = nf90_inq_varid(ncid,trim('psivolpnpsi'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_psivolpnpsi,onetwo_npsi)

  err = nf90_inq_varid(ncid,trim('elongxnpsi'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_elongxnpsi,onetwo_npsi)

  err = nf90_inq_varid(ncid,trim('triangnpsi_u'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_triangnpsi_u,onetwo_npsi)

  err = nf90_inq_varid(ncid,trim('triangnpsi_l'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_triangnpsi_l,onetwo_npsi)

  ! New sscxl added by Kinsey
  err = nf90_inq_varid(ncid,trim('sscxl'),varid)
  err = nf90_get_var(ncid,varid,onetwo_sscxl)
  
  call cub_spline(onetwo_rho_mhd_gridnpsi,onetwo_rmajavnpsi,onetwo_npsi,&
       onetwo_rho_grid,rmaj,nx)

  call cub_spline(onetwo_rho_mhd_gridnpsi,onetwo_rminavnpsi,onetwo_npsi,&
       onetwo_rho_grid,rmin,nx)

  call cub_spline(onetwo_rho_mhd_gridnpsi,onetwo_elongxnpsi,onetwo_npsi,&
       onetwo_rho_grid,kappa,nx)
 
  work = 0.5*(onetwo_triangnpsi_u+onetwo_triangnpsi_l)

  call cub_spline(onetwo_rho_mhd_gridnpsi,work,onetwo_npsi,&
       onetwo_rho_grid,delta,nx)

  err = nf90_close(ncid)

  dpsi(:) = onetwo_psi(:)-onetwo_psi(1)

  ! Compute torflux(a) [will be overwritten by gfile]
  torfluxa = 0.5*bcentr*onetwo_rho_grid(nx)**2

end subroutine prgen_read_iterdb_nc

subroutine reorder(x,xt,n)

  integer, intent(in) :: n
  real, dimension(n) :: x
  real, dimension(n) :: xt
  integer :: i

  do i=1,n
     xt(i) = x(n-i+1)
  enddo

end subroutine reorder
