!-----------------------------------------------------------------
! cgyro_nl_comm.F90
!
! PURPOSE:
!  Nonlinear communication routines
!-----------------------------------------------------------------

module cgyro_nl_comm

  implicit none

contains

!
! Comm is a transposea
! Reminder: nc ~= n_radial*n_theta
! First half of the transpose is done locally
!  from (theta,radial,nv_loc,nt_loc) -> (radial, nt_loc, theta, nv_loc)
! Then AlltoAll finishes the transpose
!  from (radial, nt_loc, theta, nv_loc_1, nv_loc_2) x toroidal_procs -> (radial, nt_loc, theta, nv_loc_1 , toroidal_procs) x nv_loc_2
! Implies nv_loc_2 == toroidal_procs
!

! NOTE: call cgyro_nl_fftw_comm1/2_async before cgyro_nl_fftw
subroutine cgyro_nl_fftw_comm1_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,itor
  integer :: iexch

  call timer_lib_in('nl_mem')

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(iexch) &
!$acc&         present(ic_c,h_x,fpack) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial) default(none)
#else
!$omp parallel do collapse(4) private(iexch)
#endif
  do it=1,n_theta
   do iv_loc_m=1,nv_loc
    do itor=nt1,nt2
     do ir=1,n_radial
       iexch = iv_loc_m + (it-1)*nv_loc
       fpack(ir,itor-nt1+1,iexch) = h_x(ic_c(ir,it),iv_loc_m,itor)
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(fpack)
#endif
    do iexch=nv_loc*n_theta+1,nsplit*n_toroidal_procs
      fpack(1:n_radial,1:nt_loc,iexch) = (0.0,0.0)
    enddo
  endif

  call parallel_slib_f_nc_async(fpack,f_nl,f_req)

  call timer_lib_out('nl_mem')

end subroutine cgyro_nl_fftw_comm1_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm1_test
  use parallel_lib
  use cgyro_globals

  implicit none

  call parallel_slib_test(f_req)

end subroutine cgyro_nl_fftw_comm1_test

subroutine cgyro_nl_fftw_comm1_r(ij)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  integer :: ir,it,iv_loc_m,ic_loc_m,itor
  integer :: iexch
  complex :: my_psi
  real :: psi_mul

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc(f_nl,fpack)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  psi_mul = (q*rho/rmin)*(2*pi/length)

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch,ic_loc_m,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(iexch,ic_loc_m,my_psi) &
!$acc&         present(ic_c,px,rhs,fpack) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial) copyin(ij) default(none)
#else
!$omp parallel do collapse(4) private(iexch,ic_loc_m,my_psi)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           iexch = iv_loc_m + (it-1)*nv_loc
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              my_psi = fpack(ir,itor-nt1+1,iexch)
           endif           
           if (itor == 0) then
              my_psi = my_psi*zf_scale
           endif
           
           ! RHS -> -[f,g] = [f,g]_{r,-alpha}
           rhs(ic_loc_m,iv_loc_m,itor,ij) = rhs(ic_loc_m,iv_loc_m,itor,ij)+psi_mul*my_psi
        enddo
      enddo
    enddo
  enddo

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw_comm1_r


!
! Comm2 is a transpose
! Reminder: nc ~= n_radial*n_theta
! First half of the transpose is done locally with sub-sampling
!  from (n_field,n_theta,n_radial,nt_loc) -> (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! Then AlltoAll finishes the transpose
!  (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_proc)xn_toroidal_proc -> (n_field,n_radial,n_jtheta,nt_loc,n_toroida_procl)xn_toroidal_proc
! 

subroutine cgyro_nl_fftw_comm2_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,it_loc,itm,itl,itf
  integer :: itor,mytor
  integer :: iltheta_min
  complex :: gval

  call timer_lib_in('nl_mem')

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(5) &
!$omp&         private(itor,it,iltheta_min,mytor,gval)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(5) independent &
!$acc&         private(itor,it,iltheta_min,mytor,gval) &
!$acc&         present(field,gpack) &
!$acc&         present(n_toroidal_procs,nt_loc,n_jtheta,nv_loc,nt1) &
!$acc&         present(n_theta,n_radial,n_field,nsplit) &
!$acc&         default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(it_loc,itor,mytor,it,ir,iltheta_min,gval)
#endif
  do itm=1,n_toroidal_procs
   do itl=1,nt_loc
    do it_loc=1,n_jtheta
     do ir=1,n_radial
      do itf=1,n_field
       iltheta_min = 1+((itm-1)*nsplit)/nv_loc
       it = it_loc+iltheta_min-1
       itor = itl+(itm-1)*nt_loc
       gval = (0.0,0.0)
       if (it <= n_theta) then
         mytor = nt1+itl-1
         ! ic_c(ir,it) = (ir-1)*n_theta+it
         gval = field(itf,(ir-1)*n_theta+it,mytor)
       endif
       ! else just padding
       gpack(itf,ir,it_loc,itor) = gval
      enddo
     enddo
    enddo
   enddo
  enddo

  call parallel_slib_f_fd_async(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)

  call timer_lib_out('nl_mem')

end subroutine cgyro_nl_fftw_comm2_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm2_test
  use parallel_lib
  use cgyro_globals

  implicit none

  call parallel_slib_test(g_req)

end subroutine cgyro_nl_fftw_comm2_test

end module cgyro_nl_comm

