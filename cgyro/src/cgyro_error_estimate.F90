!------------------------------------------------------------------------------
! cgyro_error_estimate.f90
!
! PURPOSE:
!  Compute two time-integration error estimates based on 
!  (1) difference between computed field and quadratic extrapolation, 
!  (2) 3rd-order estimate of collisionless error based on RK algebra
!------------------------------------------------------------------------------

subroutine cgyro_error_estimate

  use mpi
  use cgyro_globals
  use cgyro_io
  use timer_lib

  implicit none

  integer :: i_f,itor
  real, dimension(2) :: norm_loc,norm
  real, dimension(2) :: pair_loc,pair
  real, dimension(2) :: error_loc

  real :: norm_loc_s,error_loc_s,h_s,r_s

#if (!defined(OMPGPU)) && defined(_OPENACC)
  ! launch Estimate of collisionless error via 3rd-order linear estimate async ahead of time on GPU/ACC
  ! CPU-only and OMPGPU code will work on it later
  ! NOTE: If I have multiple itor, sum them all together
  h_s=0.0
  r_s=0.0
!$acc parallel loop collapse(3) independent gang vector &
!$acc&         present(h_x,rhs(:,:,:,1)) reduction(+:h_s,r_s) async(2)
  do itor=nt1,nt2
   do iv_loc=1,nv_loc
     do ic=1,nc
        h_s = h_s + abs(h_x(ic,iv_loc,itor))
        r_s = r_s + abs(rhs(ic,iv_loc,itor,1))
     enddo
   enddo
  enddo
#endif

  call timer_lib_in('field')

  norm_loc_s = 0.0
  error_loc_s = 0.0

  ! field_olds are always only in system memory... too expensive to keep in GPU memory
  ! assuming field was already synched to system memory
!$omp parallel do collapse(3) reduction(+:norm_loc_s,error_loc_s)
  do itor=nt1,nt2
   do ic=1,nc
     do i_f=1,n_field

        ! 1. Estimate of total (field) error via quadratic interpolation

        field_loc(i_f,ic,itor) = 3*field_old(i_f,ic,itor) - &
                3*field_old2(i_f,ic,itor) + &
                field_old3(i_f,ic,itor)
        field_dot(i_f,ic,itor) = (3*field(i_f,ic,itor) - &
                4*field_old(i_f,ic,itor) + &
                field_old2(i_f,ic,itor) )/(2*delta_t)

        ! Define norm and error for each mode number n
        norm_loc_s  = norm_loc_s  + abs(field(i_f,ic,itor))
        error_loc_s = error_loc_s + abs(field(i_f,ic,itor)-field_loc(i_f,ic,itor))

        ! save old values for next iteration
        field_old3(i_f,ic,itor) = field_old2(i_f,ic,itor)
        field_old2(i_f,ic,itor) = field_old(i_f,ic,itor)
        field_old(i_f,ic,itor)  = field(i_f,ic,itor)
     enddo
   enddo
  enddo

  norm_loc(1)  = norm_loc_s
  error_loc(1) = error_loc_s

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc) &
!$acc&         present(cap_h_c_dot,cap_h_c,cap_h_c_old,cap_h_c_old2) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) copyin(delta_t) default(none)
#else
!$omp parallel do collapse(3) private(iv_loc)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        cap_h_c_dot(ic,iv_loc,itor) = (3*cap_h_c(ic,iv_loc,itor) - &
                4*cap_h_c_old(ic,iv_loc,itor) + &
                cap_h_c_old2(ic,iv_loc,itor) )/(2*delta_t)
        cap_h_c_old2(ic,iv_loc,itor) = cap_h_c_old(ic,iv_loc,itor)
        cap_h_c_old(ic,iv_loc,itor) = cap_h_c(ic,iv_loc,itor)
     enddo
   enddo
  enddo

  call timer_lib_out('field')

  ! 2. Estimate of collisionless error via 3rd-order linear estimate

  call timer_lib_in('str')

#if (!defined(OMPGPU)) && defined(_OPENACC)
  ! wait for the async GPU compute to be completed
!$acc wait(2)
#else
  ! NOTE: If I have multiple itor, sum them all together
  h_s=0.0
  r_s=0.0
#if defined(OMPGPU)
  ! no async for OMPG{U for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&    reduction(+:h_s,r_s)
#else
!$omp parallel do collapse(3) reduction(+:h_s,r_s)
#endif
  do itor=nt1,nt2
   do iv_loc=1,nv_loc
     do ic=1,nc
        h_s = h_s + abs(h_x(ic,iv_loc,itor))
        r_s = r_s + abs(rhs(ic,iv_loc,itor,1))
     enddo
   enddo
  enddo

#endif

  pair_loc(1) = h_s
  pair_loc(2) = r_s

  call timer_lib_out('str')

  call timer_lib_in('str_comm')

  ! sum over velocity space
  call MPI_ALLREDUCE(pair_loc,&
       pair,&
       2,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  norm_loc(2) = pair(1)
  error_loc(2) = pair(2)

  ! Get sum of all errors
  call MPI_ALLREDUCE(error_loc, &
       integration_error, &
       2, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  ! Get sum of all norms
  call MPI_ALLREDUCE(norm_loc, &
       norm, &
       2, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call timer_lib_out('str_comm')

  ! 1: total (field) error
  ! 2: collisionless error
  integration_error = integration_error/norm

  ! Trigger code shutdown on large collisionless error
  if (integration_error(2) > 5e2*error_tol .and. i_time > 2) then
     call cgyro_error('Integration error exceeded limit.')
  endif

end subroutine cgyro_error_estimate
