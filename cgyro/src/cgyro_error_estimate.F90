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
  use cgyro_field_mod, only : cgyro_field_e_compute, field_cur, field_dot
  use timer_lib

  implicit none

  integer :: i_f,itor
  real, dimension(2) :: norm_loc,norm
  real, dimension(2) :: pair_loc,pair
  real, dimension(2) :: error_loc

  real :: norm_loc_s,error_loc_s,h_s,r_s

  ! epar
  complex :: wderiv,thfac
  integer :: ir,it
  ! ir loop specific
  integer :: itorbox
  integer :: jr0(0:2)   ! n_theta*(pre-compute jr-1)
  ! it loop specific
  integer :: id
  integer :: itd   ! precompute modulo(it+id-1,n_theta)+1, use for iteration
  integer :: itd_class
  integer :: jc

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

  call cgyro_field_e_compute(delta_t, norm_loc_s,error_loc_s)

  call timer_lib_in('field')

  norm_loc(1)  = norm_loc_s
  error_loc(1) = error_loc_s

  if (nonlinear_flag == 0) then

     epar(:,nt1:nt2) = 0.0
     if (n_field >= 2) epar(:,nt1:nt2) = -field_dot(2,:,nt1:nt2)
!$omp parallel do collapse(2) private(itorbox,jr0,wderiv,ic,itd,itd_class,jc,thfac,id)
     do itor=nt1,nt2
        do ir=1,n_radial
           itorbox = itor*box_size*sign_qs

           jr0(0) = n_theta*modulo(ir-itorbox-1,n_radial)
           jr0(1) = n_theta*(ir-1)
           jr0(2) = n_theta*modulo(ir+itorbox-1,n_radial)

           do it=1,n_theta

              wderiv = 0.0
              ic = (ir-1)*n_theta + it

              itd = n_theta+it-nup_theta
              itd_class = 0
              jc = jr0(itd_class)+itd
              thfac = thfac_itor(itd_class,itor)
              do id=-nup_theta,nup_theta
                 if (itd > n_theta) then
                    itd = itd - n_theta
                    itd_class = itd_class + 1
                    jc = jr0(itd_class)+itd
                    thfac = thfac_itor(itd_class,itor)
                 endif
                 wderiv = wderiv + thfac*cderiv(id)*field_cur(1,jc,itor)
                 itd = itd + 1
                 jc = jc + 1
              enddo

              ! E_par = - (1/c) (d/dt) A_par - b dot grad phi
              epar(ic,itor) = epar(ic,itor) - wderiv/(q*rmaj*g_theta(it))

           enddo
        enddo
     enddo

  endif

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
  ! no async for OMPGPU for now
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
