!-----------------------------------------------------------------
! cgyro_upwind.f90
!
! PURPOSE:
!  Compute corrected distribution used in conservative advection scheme:
!
!                  /
!              J0  | dv J0 |vp| g 
!                  /
!  |vp| g -  ----------------------
!                   /
!                   | dv J0^2 
!                   /
!
!-----------------------------------------------------------------

subroutine cgyro_upwind

  use cgyro_globals
  use mpi
  use timer_lib

  implicit none

  integer :: is,ie,ix
  complex, dimension(nc,n_species) :: res_loc 
  complex, dimension(nc,n_species) :: res

  call timer_lib_in('str')
  res_loc(:,:) = (0.0,0.0)

!$omp parallel private(iv_loc,ic,is)
!$omp do reduction(+:res_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        res_loc(ic,is) = res_loc(ic,is)+upfac1(ic,iv_loc)*g_x(ic,iv_loc)
     enddo
  enddo
!$omp end do
!$omp end parallel

  call timer_lib_out('str')

  call timer_lib_in('str_comm')
  call MPI_ALLREDUCE(res_loc(:,:),&
       res(:,:),&
       size(res(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
  call timer_lib_out('str_comm')

  call timer_lib_in('str')

!$omp parallel do private(iv_loc,is,ix,ie,ic)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        g_x(ic,iv_loc) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc) &
             -upfac2(ic,iv_loc)*res(ic,is)
     enddo
  enddo

  call timer_lib_out('str')

end subroutine cgyro_upwind
