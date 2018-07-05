!--------------------------------------------------------------
! EXPRO_pread.f90
!
! PURPOSE:
!  Read and broadcast experimental profiles (parallel).
!
! NOTES:
!  Variable naming corresponds to input.profiles 
!  documentation at 
!
!    http://fusion.gat.com/theory/input.profiles
!--------------------------------------------------------------

subroutine EXPRO_pread

  use mpi
  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer :: ierr
  integer :: i_proc

  call MPI_COMM_RANK(comm,i_proc,ierr)

  if (i_proc == 0) then
     call EXPRO_read_driver
  endif

  call MPI_BCAST(EXPRO_shot,1,MPI_INTEGER,0,comm,ierr)
  call MPI_BCAST(EXPRO_n_ion,1,MPI_INTEGER,0,comm,ierr)
  call MPI_BCAST(EXPRO_n_exp,1,MPI_INTEGER,0,comm,ierr)
  call MPI_BCAST(EXPRO_b_ref,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_ip_exp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_rvbv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_arho,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  ! 1-5
  call MPI_BCAST(EXPRO_rho,size(EXPRO_rho),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_rmin,size(EXPRO_rho),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_rmaj,size(EXPRO_rmaj),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_q,size(EXPRO_q),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_kappa,size(EXPRO_kappa),MPI_DOUBLE_PRECISION,0,comm,ierr)

  ! 6-10
  call MPI_BCAST(EXPRO_delta,size(EXPRO_delta),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_te,size(EXPRO_te),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_ne,size(EXPRO_ne),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_z_eff,size(EXPRO_z_eff),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_w0,size(EXPRO_w0),MPI_DOUBLE_PRECISION,0,comm,ierr)

  ! 11-15
  call MPI_BCAST(EXPRO_flow_mom,size(EXPRO_flow_mom),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_e,size(EXPRO_pow_e),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_i,size(EXPRO_pow_i),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_ei,size(EXPRO_pow_ei),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_zeta,size(EXPRO_zeta),MPI_DOUBLE_PRECISION,0,comm,ierr)

  ! 16-20
  call MPI_BCAST(EXPRO_flow_beam,size(EXPRO_flow_beam),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_flow_wall,size(EXPRO_flow_wall),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_zmag,size(EXPRO_zmag),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_ptot,size(EXPRO_ptot),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_polflux,size(EXPRO_polflux),MPI_DOUBLE_PRECISION,0,comm,ierr)
  ! dummy
 
  ! 21-25
  call MPI_BCAST(EXPRO_ni,size(EXPRO_ni),MPI_DOUBLE_PRECISION,0,comm,ierr)
 
  ! 26-30
  call MPI_BCAST(EXPRO_ti,size(EXPRO_ti),MPI_DOUBLE_PRECISION,0,comm,ierr)

  ! 31-35
  call MPI_BCAST(EXPRO_vtor,size(EXPRO_vtor),MPI_DOUBLE_PRECISION,0,comm,ierr)

  ! 36-40
  call MPI_BCAST(EXPRO_vpol,size(EXPRO_vpol),MPI_DOUBLE_PRECISION,0,comm,ierr)

  ! 41-45
  call MPI_BCAST(EXPRO_pow_e_fus,size(EXPRO_pow_e_fus),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_i_fus,size(EXPRO_pow_i_fus),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_e_sync,size(EXPRO_pow_e_sync),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_e_brem,size(EXPRO_pow_e_brem),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_e_line,size(EXPRO_pow_e_line),MPI_DOUBLE_PRECISION,0,comm,ierr)

  ! 46-50
  call MPI_BCAST(EXPRO_pow_e_aux,size(EXPRO_pow_e_aux),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_pow_i_aux,size(EXPRO_pow_i_aux),MPI_DOUBLE_PRECISION,0,comm,ierr)

  call MPI_BCAST(EXPRO_sbeame,size(EXPRO_sbeame),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_sbcx,size(EXPRO_sbcx),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call MPI_BCAST(EXPRO_sscxl,size(EXPRO_sscxl),MPI_DOUBLE_PRECISION,0,comm,ierr)

  if (EXPRO_nfourier > 0) then

     call MPI_BCAST(EXPRO_geo,&
          size(EXPRO_geo),&
          MPI_DOUBLE_PRECISION,&
          0,&
          comm,&
          ierr)

     call MPI_BCAST(EXPRO_dgeo,&
          size(EXPRO_dgeo),&
          MPI_DOUBLE_PRECISION,&
          0,&
          comm,&
          ierr)

  endif

  ! Now all processes have these profiles, so no further 
  ! broadcasts are required.

  call EXPRO_compute_derived

end subroutine EXPRO_pread

