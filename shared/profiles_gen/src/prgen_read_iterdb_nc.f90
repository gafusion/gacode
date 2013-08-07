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
  err = nf90_open(raw_data_file,NF90_NOWRITE,ncid)

  err = nf90_inq_varid(ncid,trim('shot'),varid)
  err = nf90_get_var(ncid,varid,onetwo_ishot)

  err = nf90_inq_varid(ncid,trim('nj'),varid)
  err = nf90_get_var(ncid,varid,onetwo_nj)

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
  err = nf90_get_var(ncid,varid,onetwo_Btor)

  nx = onetwo_nj

  call allocate_internals
  call allocate_iterdb_vars

  allocate(work(onetwo_npsi))

  err = nf90_inq_varid(ncid,trim('rho_grid'),varid)
  err = nf90_get_var(ncid,varid,onetwo_rho_grid)

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
  err = nf90_get_var(ncid,varid,onetwo_press)

  err = nf90_inq_varid(ncid,trim('pressb'),varid)
  err = nf90_get_var(ncid,varid,onetwo_pressb(:,1:onetwo_nbion))

  err = nf90_inq_varid(ncid,trim('sbcx'),varid)
  err = nf90_get_var(ncid,varid,onetwo_sbcx(:,1:onetwo_nion))

  err = nf90_inq_varid(ncid,trim('sion'),varid)
  err = nf90_get_var(ncid,varid,onetwo_sion(:,1:onetwo_nion))

  sbcx_d(:) = 0.0
  sion_d(:) = 0.0
  do i=1,onetwo_nprim
     sbcx_d(:) = sbcx_d(:)+onetwo_sbcx(:,i)
     sion_d(:) = sion_d(:)+onetwo_sion(:,i)
  enddo

  err = nf90_inq_varid(ncid,trim('sbeam'),varid)
  err = nf90_get_var(ncid,varid,onetwo_sbeam)

  err = nf90_inq_varid(ncid,trim('zeff'),varid)
  err = nf90_get_var(ncid,varid,onetwo_zeff)

  err = nf90_inq_varid(ncid,trim('angrot'),varid)
  err = nf90_get_var(ncid,varid,onetwo_angrot)

  err = nf90_inq_varid(ncid,trim('dpedt'),varid)
  err = nf90_get_var(ncid,varid,onetwo_dpedt)

  err = nf90_inq_varid(ncid,trim('dpidt'),varid)
  err = nf90_get_var(ncid,varid,onetwo_dpidt)

  err = nf90_inq_varid(ncid,trim('qbeame'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qbeame)

  err = nf90_inq_varid(ncid,trim('qdelt'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qdelt)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qdelt(:) = -1*onetwo_qdelt(:)

  err = nf90_inq_varid(ncid,trim('qbeami'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qbeami)

  err = nf90_inq_varid(ncid,trim('qrfe'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qrfe)

  err = nf90_inq_varid(ncid,trim('qrfi'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qrfi)

  err = nf90_inq_varid(ncid,trim('qione'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qione)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qione(:) = -1*onetwo_qione(:)

  err = nf90_inq_varid(ncid,trim('qioni'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qioni)

  err = nf90_inq_varid(ncid,trim('qcx'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qcx)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qcx(:) = -1*onetwo_qcx(:)

  err = nf90_inq_varid(ncid,trim('qfuse'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qfuse)

  err = nf90_inq_varid(ncid,trim('qfusi'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qfusi)

  err = nf90_inq_varid(ncid,trim('qrad'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qrad)
  !Correct for difference in sign convention with ASCII format iterdb
  onetwo_qrad(:) = -1*onetwo_qrad(:)

  err = nf90_inq_varid(ncid,trim('qohm'),varid)
  err = nf90_get_var(ncid,varid,onetwo_qohm)

  err = nf90_inq_varid(ncid,trim('storqueb'),varid)
  err = nf90_get_var(ncid,varid,onetwo_storqueb)

  err = nf90_inq_varid(ncid,trim('rho_mhd_gridnpsi'),varid)
  err = nf90_get_var(ncid,varid,work)

  call reorder(work,onetwo_rho_mhd_gridnpsi,onetwo_npsi)

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

  call cub_spline(onetwo_rho_mhd_gridnpsi,onetwo_rmajavnpsi,onetwo_npsi,&
       onetwo_rho_grid,rmaj,onetwo_nj)

  call cub_spline(onetwo_rho_mhd_gridnpsi,onetwo_rminavnpsi,onetwo_npsi,&
       onetwo_rho_grid,rmin,onetwo_nj)

  call cub_spline(onetwo_rho_mhd_gridnpsi,onetwo_elongxnpsi,onetwo_npsi,&
       onetwo_rho_grid,kappa,onetwo_nj)

  work = 0.5*(onetwo_triangnpsi_u+onetwo_triangnpsi_l)

  call cub_spline(onetwo_rho_mhd_gridnpsi,work,onetwo_npsi,&
       onetwo_rho_grid,delta,onetwo_nj)

  err = nf90_close(ncid)

  dpsi(:) = onetwo_psi(:)-onetwo_psi(1)

  ! No squareness 
  zeta(:) = 0.0

  ! No elevation 
  zmag(:) = 0.0

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
