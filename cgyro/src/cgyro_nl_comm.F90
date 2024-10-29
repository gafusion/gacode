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
subroutine cgyro_nl_fftw_comm1_async
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

end subroutine cgyro_nl_fftw_comm1_async

subroutine cgyro_nl_fftw_comm3_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  if (nsplitB > 0) then
    call timer_lib_in('nl_comm')

    call parallel_slib_f_nc_async(nsplitB,fpackB,fB_nl,fB_req)
    fB_req_valid = .TRUE.

    call timer_lib_out('nl_comm')
  endif
end subroutine cgyro_nl_fftw_comm3_async

subroutine cgyro_nl_fftw_comm1_r(ij)
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

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd_async(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  g_req_valid = .TRUE.

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm2_async


!
! Comm2 is a transpose
! Reminder: nc ~= n_radial*n_theta
! First half of the transpose is done locally with sub-sampling
!  from (n_field,n_theta,n_radial,nt_loc) -> (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! Then AlltoAll finishes the transpose
!  (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_proc)xn_toroidal_proc -> (n_field,n_radial,n_jtheta,nt_loc,n_toroida_procl)xn_toroidal_proc
!

subroutine cgyro_nl_fftw_comm2_async_stress(itf)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: itf
  integer :: ir,it,it_loc,itm,itl
  integer :: itor,mytor
  integer :: iltheta_min
  complex :: gval

  call timer_lib_in('nl_mem')

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(itor,it,iltheta_min,mytor,gval)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(4) independent &
!$acc&         private(itor,it,iltheta_min,mytor,gval) &
!$acc&         present(field,gpack) &
!$acc&         present(n_toroidal_procs,nt_loc,n_jtheta,nv_loc,nt1) &
!$acc&         present(n_theta,n_radial,n_field,nsplit) &
!$acc&         copyin(itf) &
!$acc&         default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(it_loc,itor,mytor,it,ir,iltheta_min,gval)
#endif
  do itm=1,n_toroidal_procs
   do itl=1,nt_loc
    do it_loc=1,n_jtheta
     do ir=1,n_radial
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

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd_async(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  g_req_valid = .TRUE.

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm2_async_stress


subroutine cgyro_nl_fftw_comm1_r_stress(itf)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: itf
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
!$acc&         present(ic_c,px,stress,fpackA,fpackB) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA,nsplitB) copyin(itf) default(none)
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
           stress(ic_loc_m,iv_loc_m,itor,itf) = stress(ic_loc_m,iv_loc_m,itor,itf)+psi_mul*my_psi
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
!$acc&         present(ic_c,px,stress,fpackA) copyin(psi_mul,zf_scale) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA) copyin(itf) default(none)
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
           stress(ic_loc_m,iv_loc_m,itor,itf) = stress(ic_loc_m,iv_loc_m,itor,itf)+psi_mul*my_psi
        enddo
      enddo
    enddo
  enddo

  endif ! if nsplitB>0

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw_comm1_r_stress


end module cgyro_nl_comm

