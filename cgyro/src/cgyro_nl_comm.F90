!-----------------------------------------------------------------
! cgyro_nl_comm.F90
!
! PURPOSE:
!  Nonlinear communication routines
!-----------------------------------------------------------------

module cgyro_nl_comm

  implicit none

contains

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm_test
  use parallel_lib
  use cgyro_globals

  implicit none

  if (fA_req_valid) call parallel_slib_test(fA_req)
  if (g_req_valid) call parallel_slib_test(g_req)
  if (fB_req_valid) call parallel_slib_test(fB_req)

end subroutine cgyro_nl_fftw_comm_test

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
subroutine cgyro_nl_fftw_comm1_f64_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex :: h_loc

  call timer_lib_in('nl_mem')

  if (nsplitB > 0) then

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA,fpackB) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA,nsplitB) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#endif
  do it=1,n_theta
   do iv_loc_m=1,nv_loc
    do itor=nt1,nt2
     do ir=1,n_radial
       h_loc = h_x(ic_c(ir,it),iv_loc_m,itor)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
       else
          iexch_base = 1+itor0*nsplitB
          fpackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA)) = h_loc
       endif
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA,fpackB,nsplit,nsplitA,nsplitB)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
       else
          iexch_base = 1+itor0*nsplitB
          fpackB(1:n_radial,1:nt_loc,iexch_base+(isplit0-nsplitA)) = (0.0,0.0)
       endif
    enddo
  endif

  else ! nsplitB==0

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#endif
  do it=1,n_theta
   do iv_loc_m=1,nv_loc
    do itor=nt1,nt2
     do ir=1,n_radial
       h_loc = h_x(ic_c(ir,it),iv_loc_m,itor)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA,nsplit,nsplitA)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
    enddo
  endif

  endif ! if nspliB>0

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! split the comm in two, so we can start working on first as soon as it is ready
  call parallel_slib_f_nc_async(nsplitA,fpackA,fA_nl,fA_req)
  fA_req_valid = .TRUE.
  ! send only the first half, use comm3 to send the other half
  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm1_f64_async

subroutine cgyro_nl_fftw_comm1_f32_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex(KIND=REAL32) :: h_loc

  call timer_lib_in('nl_mem')

  if (nsplitB > 0) then

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA32,fpackB32) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA,nsplitB) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#endif
  do it=1,n_theta
   do iv_loc_m=1,nv_loc
    do itor=nt1,nt2
     do ir=1,n_radial
       h_loc = h_x(ic_c(ir,it),iv_loc_m,itor)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA32(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
       else
          iexch_base = 1+itor0*nsplitB
          fpackB32(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA)) = h_loc
       endif
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA32,fpackB32,nsplit,nsplitA,nsplitB)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA32(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
       else
          iexch_base = 1+itor0*nsplitB
          fpackB32(1:n_radial,1:nt_loc,iexch_base+(isplit0-nsplitA)) = (0.0,0.0)
       endif
    enddo
  endif

  else ! nsplitB==0

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA32) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#endif
  do it=1,n_theta
   do iv_loc_m=1,nv_loc
    do itor=nt1,nt2
     do ir=1,n_radial
       h_loc = h_x(ic_c(ir,it),iv_loc_m,itor)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA32(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA32,nsplit,nsplitA)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA32(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
    enddo
  endif

  endif ! if nspliB>0

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! split the comm in two, so we can start working on first as soon as it is ready
  call parallel_slib_f_nc32_async(nsplitA,fpackA32,fA_nl32,fA_req)
  fA_req_valid = .TRUE.
  ! send only the first half, use comm3 to send the other half
  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm1_f32_async

subroutine cgyro_nl_fftw_comm1_async

  use cgyro_globals

  implicit none
  !-----------------------------------

  if (nl_single_flag > 1) then
    call cgyro_nl_fftw_comm1_f32_async
  else
    call cgyro_nl_fftw_comm1_f64_async
  endif
end subroutine cgyro_nl_fftw_comm1_async

subroutine cgyro_nl_fftw_comm3_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  if (nsplitB > 0) then
    call timer_lib_in('nl_comm')

    if (nl_single_flag > 1) then
      call parallel_slib_f_nc32_async(nsplitB,fpackB32,fB_nl32,fB_req)
    else
      call parallel_slib_f_nc_async(nsplitB,fpackB,fB_nl,fB_req)
    endif
    fB_req_valid = .TRUE.

    call timer_lib_out('nl_comm')
  endif
end subroutine cgyro_nl_fftw_comm3_async

subroutine cgyro_nl_fftw_comm1_r64(ij)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  integer :: ir,it,iv_loc_m,ic_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex :: my_psi
  real :: psi_mul

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc_wait(nsplitA,fA_nl,fpackA,fA_req)
  fA_req_valid = .FALSE.
  if (nsplitB > 0) then
    ! no major compute to overlap
    call parallel_slib_r_nc_wait(nsplitB,fB_nl,fpackB,fB_req)
    fB_req_valid = .FALSE.
  endif
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  psi_mul = (q*rho/rmin)*(2*pi/length)

  if (nsplitB > 0) then

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(ic_c,px,rhs,fpackA,fpackB) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA,nsplitB) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              if (isplit0 < nsplitA) then
                 iexch_base = 1+itor0*nsplitA
                 my_psi = fpackA(ir,itor-nt1+1,iexch_base+isplit0)
              else
                 iexch_base = 1+itor0*nsplitB
                 my_psi = fpackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))
              endif
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

  else ! nsplitB==0

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(ic_c,px,rhs,fpackA) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              iexch_base = 1+itor0*nsplitA
              my_psi = fpackA(ir,itor-nt1+1,iexch_base+isplit0)
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

  endif ! if nsplitB>0

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw_comm1_r64


subroutine cgyro_nl_fftw_comm1_r32(ij)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  integer :: ir,it,iv_loc_m,ic_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex(KIND=REAL32) :: my_psi
  real(KIND=REAL32) :: psi_mul

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc32_wait(nsplitA,fA_nl32,fpackA32,fA_req)
  fA_req_valid = .FALSE.
  if (nsplitB > 0) then
    ! no major compute to overlap
    call parallel_slib_r_nc32_wait(nsplitB,fB_nl32,fpackB32,fB_req)
    fB_req_valid = .FALSE.
  endif
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  psi_mul = (q*rho/rmin)*(2*pi/length)

  if (nsplitB > 0) then

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(ic_c,px,rhs,fpackA32,fpackB32) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA,nsplitB) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              if (isplit0 < nsplitA) then
                 iexch_base = 1+itor0*nsplitA
                 my_psi = fpackA32(ir,itor-nt1+1,iexch_base+isplit0)
              else
                 iexch_base = 1+itor0*nsplitB
                 my_psi = fpackB32(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))
              endif
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

  else ! nsplitB==0

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(ic_c,px,rhs,fpackA32) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              iexch_base = 1+itor0*nsplitA
              my_psi = fpackA32(ir,itor-nt1+1,iexch_base+isplit0)
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

  endif ! if nsplitB>0

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw_comm1_r32

subroutine cgyro_nl_fftw_comm1_r(ij)
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  if (nl_single_flag .EQ. 0) then
    call cgyro_nl_fftw_comm1_r64(ij)
  else
    call cgyro_nl_fftw_comm1_r32(ij)
  endif

end subroutine cgyro_nl_fftw_comm1_r

subroutine cgyro_nl_fftw_comm1_r_triad(ij)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------
  
  integer :: is,ix,ie
  integer :: id,itd,itd_class,jr0(0:2),itorbox,jc
  integer :: ir,it,iv_loc_m,ic_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex :: my_psi
  real :: psi_mul

  real :: dv,dvr,dvp,rval,rval2
  complex :: cprod,cprod2,thfac

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector independent &
!$acc&         present(triad_loc_old,triad_loc) &
!$acc&         present(nt1,nt2,n_radial,n_species) default(none)
#else
!$omp parallel do private(ir,is)
#endif
  do itor=nt1,nt2
      do ir=1,n_radial
        do is=1,n_species
          triad_loc_old(is,ir,itor,3) = triad_loc(is,ir,itor,3)
          triad_loc_old(is,ir,itor,4) = triad_loc(is,ir,itor,4)
          triad_loc(is,ir,itor,:) = 0.0
        enddo
      enddo
  enddo

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc_wait(nsplitA,fA_nl,fpackA,fA_req)
  call parallel_slib_r_nc_wait(nsplitA,eA_nl,epackA,eA_req)
  fA_req_valid = .FALSE.
  if (nsplitB > 0) then
    ! no major compute to overlap
    call parallel_slib_r_nc_wait(nsplitB,fB_nl,fpackB,fB_req)
    call parallel_slib_r_nc_wait(nsplitB,eB_nl,epackB,eB_req)
    fB_req_valid = .FALSE.
  endif
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  psi_mul = (q*rho/rmin)*(2*pi/length)

  if (nsplitB > 0) then

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(id,itorbox,jr0,jc,itd,itd_class,thfac) &
!$omp&         private(ic_loc_m,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         private(id,itorbox,jr0,jc,itd,itd_class,thfac) &
!$acc&         present(h_x,g_x,cap_h_c,cap_h_ct,field,dvjvec_c,jvec_c) &
!$acc&         present(ic_c,is_v,ix_v,ie_v,w_exi,w_theta,dens2_rot,z,temp) &
!$acc&         present(omega_stream,vel,xi,thfac_itor,cderiv,uderiv) &
!$acc&         present(px,rhs,fpackA,fpackB,epackA,epackB,diss_r,triad_loc) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA,nsplitB) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,is,ix,ie,dv,dvr,dvp,rval,rval2,cprod,cprod2) &
!$omp&         private(id,itorbox,jr0,jc,itd,itd_class,thfac)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do ir=1,n_radial
        do it=1,n_theta
           itorbox = itor*box_size*sign_qs
           jr0(0) = n_theta*modulo(ir-itorbox-1,n_radial)
           jr0(1) = n_theta*(ir-1)
           jr0(2) = n_theta*modulo(ir+itorbox-1,n_radial)

           ic_loc_m = ic_c(ir,it)

           is = is_v(iv_loc_m +nv1 -1 )
           ix = ix_v(iv_loc_m +nv1 -1 )
           ie = ie_v(iv_loc_m +nv1 -1 )
           dv = w_exi(ie,ix)
           dvr  = w_theta(it)*dens2_rot(it,is)*dv
           dvp = w_theta(it)*dv*(ir-1-nx0/2)**2

           ! Density moment
           cprod = w_theta(it)*cap_h_c(ic_loc_m,iv_loc_m,itor)*dvjvec_c(1,ic_loc_m,iv_loc_m,itor)/z(is)
           cprod = -(dvr*z(is)/temp(is)*field(1,ic_loc_m,itor)-cprod)
           cprod2= ( jvec_c(1,ic_loc_m,iv_loc_m,itor)*z(is)/temp(is) )*conjg(field(1,ic_loc_m,itor) )

           iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
           itor0 = iexch0/nsplit
           isplit0 = modulo(iexch0,nsplit)
           if (isplit0 < nsplitA) then
              iexch_base = 1+itor0*nsplitA
              my_psi = fpackA(ir,itor-nt1+1,iexch_base+isplit0)

              ! 1. Triad energy transfer (all)
              triad_loc(is,ir,itor,1) = triad_loc(is,ir,itor,1) &
               + fpackA(ir,itor-nt1+1,iexch_base+isplit0)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr*psi_mul
              ! 2. Triad energy transfer ( {NZF-NZF} coupling )
              if (itor == 0) then
                ! New : Diag. direct ZF production [A. Ishizawa PRL 2019 ]
                triad_loc(is,ir,itor,2) = triad_loc(is,ir,itor,2)  &
               + epackA(ir,itor-nt1+1,iexch_base+isplit0)*cprod2*dvr*psi_mul
              else
              triad_loc(is,ir,itor,2) = triad_loc(is,ir,itor,2) &
               + epackA(ir,itor-nt1+1,iexch_base+isplit0)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr*psi_mul
              endif
           else
              iexch_base = 1+itor0*nsplitB
              my_psi = fpackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))

              ! 1. Triad energy transfer (all)
              triad_loc(is,ir,itor,1) = triad_loc(is,ir,itor,1) &
               + fpackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr*psi_mul
              ! 2. Triad energy transfer ( {NZF-NZF} coupling )
              if (itor == 0) then
                ! New : Diag. direct ZF production [A. Ishizawa PRL 2019 ]
                triad_loc(is,ir,itor,2) = triad_loc(is,ir,itor,2)  &
               + epackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))*cprod2*dvr*psi_mul
              else
                triad_loc(is,ir,itor,2) = triad_loc(is,ir,itor,2)  &
               + epackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr*psi_mul
              endif
           endif

           ! 3. Entropy 
           triad_loc(is,ir,itor,3) = triad_loc(is,ir,itor,3) &
                + cap_h_c(ic_loc_m,iv_loc_m,itor)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr &
                - field(1,ic_loc_m,itor)*conjg(field(1,ic_loc_m,itor))*(z(is)/temp(is))**2*dvr &
                - 2.0*cprod*conjg(field(1,ic_loc_m,itor))*(z(is)/temp(is))
           ! 4. Field potential
           triad_loc(is,ir,itor,4) = triad_loc(is,ir,itor,4)  & 
                + sum( field(:,ic_loc_m,itor)*conjg(field(:,ic_loc_m,itor)) )*dvp
           ! 5. Diss. (radial)
           triad_loc(is,ir,itor,5) = triad_loc(is,ir,itor,5)  &  
                + diss_r(ic_loc_m,iv_loc_m,itor)*h_x(ic_loc_m,iv_loc_m,itor)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr
           ! 6. Diss. (theta )
           rval = omega_stream(it,is,itor)*vel(ie)*xi(ix)
           rval2 = abs(omega_stream(it,is,itor))
           cprod = 0.0 
           cprod2= 0.0

          !icd_c(ic, id, itor)     = ic_c(jr,modulo(it+id-1,n_theta)+1)
          !jc = icd_c(ic, id, itor)
          !dtheta(ic, id, itor)    := cderiv(id)*thfac
          !dtheta_up(ic, id, itor) := uderiv(id)*thfac*up_theta
          itd = n_theta+it-nup_theta
          itd_class = 0
          jc = jr0(itd_class)+itd
          thfac = thfac_itor(itd_class,itor)

           do id=-nup_theta,nup_theta
              if (itd > n_theta) then
                ! move to next itd_class of compute
                itd = itd - n_theta
                itd_class = itd_class + 1
                jc = jr0(itd_class)+itd
                thfac = thfac_itor(itd_class,itor)
              endif

              cprod2 = cprod2 - rval* thfac*cderiv(id) *cap_h_c(jc,iv_loc_m,itor)
              cprod = cprod - rval2* uderiv(id)*up_theta *g_x(jc,iv_loc_m,itor)
              itd = itd + 1
              jc = jc + 1
           enddo 

           triad_loc(is,ir,itor,6) = triad_loc(is,ir,itor,6) + cprod*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr
           ! 7. Diss. (Coll. = Implicit advance - theta_streaming )
           triad_loc(is,ir,itor,7) = triad_loc(is,ir,itor,7) &
                 + ( cap_h_ct(iv_loc_m,itor,ic_loc_m)/delta_t + cprod2*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor)) )*dvr


           if ( (itor == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
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

  else ! nsplitB==0

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(id,itorbox,jr0,jc,itd,itd_class,thfac) &
!$omp&         private(ic_loc_m,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         private(id,itorbox,jr0,jc,itd,itd_class,thfac) &
!$acc&         present(h_x,g_x,cap_h_c,cap_h_ct,field,dvjvec_c,jvec_c) &
!$acc&         present(ic_c,is_v,ix_v,ie_v,w_exi,w_theta,dens2_rot,z,temp) &
!$acc&         present(omega_stream,vel,xi,thfac_itor,cderiv,uderiv) &
!$acc&         present(px,rhs,fpackA,epackA,diss_r,triad_loc) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,is,ix,ie,dv,dvr,dvp,rval,rval2,cprod,cprod2) &
!$omp&         private(id,itorbox,jr0,jc,itd,itd_class,thfac)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do ir=1,n_radial
        do it=1,n_theta
           itorbox = itor*box_size*sign_qs
           jr0(0) = n_theta*modulo(ir-itorbox-1,n_radial)
           jr0(1) = n_theta*(ir-1)
           jr0(2) = n_theta*modulo(ir+itorbox-1,n_radial)

           ic_loc_m = ic_c(ir,it)

           is = is_v(iv_loc_m +nv1 -1 )
           ix = ix_v(iv_loc_m +nv1 -1 )
           ie = ie_v(iv_loc_m +nv1 -1 )
           dv = w_exi(ie,ix)
           dvr  = w_theta(it)*dens2_rot(it,is)*dv
           dvp = w_theta(it)*dv*(ir-1-nx0/2)**2

           ! Density moment
           cprod = w_theta(it)*cap_h_c(ic_loc_m,iv_loc_m,itor)*dvjvec_c(1,ic_loc_m,iv_loc_m,itor)/z(is)
           cprod = -(dvr*z(is)/temp(is)*field(1,ic_loc_m,itor)-cprod)
           cprod2= ( jvec_c(1,ic_loc_m,iv_loc_m,itor)*z(is)/temp(is) )*conjg(field(1,ic_loc_m,itor) )


           iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
           itor0 = iexch0/nsplit
           isplit0 = modulo(iexch0,nsplit)
           iexch_base = 1+itor0*nsplitA
           my_psi = fpackA(ir,itor-nt1+1,iexch_base+isplit0)
     

           ! 1. Triad energy transfer (all)
           triad_loc(is,ir,itor,1) = triad_loc(is,ir,itor,1) &
               + fpackA(ir,itor-nt1+1,iexch_base+isplit0)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr*psi_mul
           ! 2. Triad energy transfer ( {NZF-NZF} coupling )
           if (itor == 0) then
             ! New : Diag. direct ZF production [A. Ishizawa PRL 2019 ]
             triad_loc(is,ir,itor,2) = triad_loc(is,ir,itor,2)  &
               + epackA(ir,itor-nt1+1,iexch_base+isplit0)*cprod2*dvr*psi_mul
           else
           triad_loc(is,ir,itor,2) = triad_loc(is,ir,itor,2) &
               + epackA(ir,itor-nt1+1,iexch_base+isplit0)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr*psi_mul
           endif
           
           ! 3. Entropy 
           triad_loc(is,ir,itor,3) = triad_loc(is,ir,itor,3) &
                + cap_h_c(ic_loc_m,iv_loc_m,itor)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr &
                - field(1,ic_loc_m,itor)*conjg(field(1,ic_loc_m,itor))*(z(is)/temp(is))**2*dvr &
                - 2.0*cprod*conjg(field(1,ic_loc_m,itor))*(z(is)/temp(is))
           ! 4. Field potential
           triad_loc(is,ir,itor,4) = triad_loc(is,ir,itor,4)  & 
                + sum( field(:,ic_loc_m,itor)*conjg(field(:,ic_loc_m,itor)) )*dvp
           ! 5. Diss. (radial)
           triad_loc(is,ir,itor,5) = triad_loc(is,ir,itor,5)  &  
                + diss_r(ic_loc_m,iv_loc_m,itor)*h_x(ic_loc_m,iv_loc_m,itor)*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr
           ! 6. Diss. (theta )
           rval = omega_stream(it,is,itor)*vel(ie)*xi(ix)
           rval2 = abs(omega_stream(it,is,itor))
           cprod = 0.0 
           cprod2= 0.0

          !icd_c(ic, id, itor)     = ic_c(jr,modulo(it+id-1,n_theta)+1)
          !jc = icd_c(ic, id, itor)
          !dtheta(ic, id, itor)    := cderiv(id)*thfac
          !dtheta_up(ic, id, itor) := uderiv(id)*thfac*up_theta
          itd = n_theta+it-nup_theta
          itd_class = 0
          jc = jr0(itd_class)+itd
          thfac = thfac_itor(itd_class,itor)

           do id=-nup_theta,nup_theta
              if (itd > n_theta) then
                ! move to next itd_class of compute
                itd = itd - n_theta
                itd_class = itd_class + 1
                jc = jr0(itd_class)+itd
                thfac = thfac_itor(itd_class,itor)
              endif

              cprod2 = cprod2 - rval* thfac*cderiv(id) *cap_h_c(jc,iv_loc_m,itor)
              cprod = cprod - rval2* uderiv(id)*up_theta *g_x(jc,iv_loc_m,itor)
              itd = itd + 1
              jc = jc + 1
           enddo 

           triad_loc(is,ir,itor,6) = triad_loc(is,ir,itor,6) + cprod*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor))*dvr
           ! 7. Diss. (Coll. = Implicit advance - theta_streaming )
           triad_loc(is,ir,itor,7) = triad_loc(is,ir,itor,7) &
                 + ( cap_h_ct(iv_loc_m,itor,ic_loc_m)/delta_t + cprod2*conjg(cap_h_c(ic_loc_m,iv_loc_m,itor)) )*dvr


           if ( (itor == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
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

  endif ! if nsplitB>0

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw_comm1_r_triad

!
! Comm2 is a transpose
! Reminder: nc ~= n_radial*n_theta
! First half of the transpose is done locally with sub-sampling
!  from (n_field,n_theta,n_radial,nt_loc) -> (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! Then AlltoAll finishes the transpose
!  (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_proc)xn_toroidal_proc -> (n_field,n_radial,n_jtheta,nt_loc,n_toroida_procl)xn_toroidal_proc
! 

subroutine cgyro_nl_fftw_comm2_f64_async
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

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd_async(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  g_req_valid = .TRUE.

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm2_f64_async

subroutine cgyro_nl_fftw_comm2_f32_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,it_loc,itm,itl,itf
  integer :: itor,mytor
  integer :: iltheta_min
  complex(KIND=REAL32) :: gval

  call timer_lib_in('nl_mem')

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(5) &
!$omp&         private(itor,it,iltheta_min,mytor,gval)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(5) independent &
!$acc&         private(itor,it,iltheta_min,mytor,gval) &
!$acc&         present(field,gpack32) &
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
       gpack32(itf,ir,it_loc,itor) = gval
      enddo
     enddo
    enddo
   enddo
  enddo

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd32_async(n_field,n_radial,n_jtheta,gpack32,g_nl32,g_req)
  g_req_valid = .TRUE.

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm2_f32_async

subroutine cgyro_nl_fftw_comm2_async

  use cgyro_globals

  implicit none
  !-----------------------------------

  if (nl_single_flag > 1) then
    call cgyro_nl_fftw_comm2_f32_async
  else
    call cgyro_nl_fftw_comm2_f64_async
  endif
end subroutine cgyro_nl_fftw_comm2_async

end module cgyro_nl_comm

