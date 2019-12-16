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
  complex, dimension(nc,n_species,2) :: res_loc 
  complex, dimension(nc,n_species,2) :: res
  complex :: res_loc_one, res_loc_two

  call timer_lib_in('str')

!$acc data create(res_loc,res) present(g_x)

!$acc parallel loop collapse(2) gang &
!$acc&         private(res_loc_one,iv) present(upfac1,is_v) default(none) 
  do is=1,n_species
     do ic=1,nc
       res_loc_one = (0.0,0.0)
       res_loc_two = (0.0,0.0)

!$acc loop vector private(iv_loc) reduction(+:res_loc_one)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_one = res_loc_one+upfac1(ic,iv_loc,1)*g_x(ic,iv_loc)
          endif
       enddo
       res_loc(ic,is,1) = res_loc_one

!$acc loop vector private(iv_loc) reduction(+:res_loc_two)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_two = res_loc_two+upfac1(ic,iv_loc,2)*g_x(ic,iv_loc)
          endif
       enddo
       res_loc(ic,is,2) = res_loc_two
    enddo

  enddo

  call timer_lib_out('str')

  call timer_lib_in('str_comm')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(res_loc)
#else
!$acc host_data use_device(res_loc,res)
#endif

#ifdef SUMMIT
  call MPI_ALLREDUCE(res_loc(:,:,:),&
       res(:,:,:),&
       2*size(res(:,:,:)),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#else
  call MPI_ALLREDUCE(res_loc(:,:,:),&
       res(:,:,:),&
       size(res(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#endif

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(res)
#else
!$acc end host_data
#endif

  call timer_lib_out('str_comm')

  call timer_lib_in('str')

!$acc parallel loop collapse(2) independent &
!$acc&         present(is_v,ix_v,ie_v,xi,vel,upfac2) &
!$acc&         private(iv_loc,is,ix,ie) default(none)
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        g_x(ic,iv_loc) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc) &
             -upfac2(ic,iv_loc,1)*res(ic,is,1) &
             -upfac2(ic,iv_loc,2)*res(ic,is,2) 
     enddo
  enddo

!$acc end data

  call timer_lib_out('str')

end subroutine cgyro_upwind
