subroutine cgyro_field_coefficients

  use mpi
  use cgyro_globals

  implicit none

  integer :: ir,it,is,ie,ix
  real, dimension(nc) :: sum_loc
  real, dimension(:), allocatable :: pb11,pb12,pb21,pb22
 
  !-------------------------------------------------------------------------
  ! Field equation prefactors, sums.
  !
  sum_loc(:) = 0.0
!$omp parallel private(iv_loc,ic)
!$omp do reduction(+:sum_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        it = it_c(ic)
        sum_loc(ic) = sum_loc(ic)+vfac(iv_loc)*dens_rot(it,is) &
             *(1.0-jvec_c(1,ic,iv_loc)**2) 
     enddo
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLREDUCE(sum_loc,&
       sum_den_x,&
       size(sum_den_x),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (ae_flag == 1) then
     do ic=1,nc
        it = it_c(ic)
        sum_den_x(ic) = sum_den_x(ic)+dens_ele*dens_ele_rot(it)/temp_ele
     enddo
  endif
  !----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Field-solve coefficients (i.e., final numerical factors).
  !
  do ic=1,nc
     ir = ir_c(ic) 
     it = it_c(ic)
     if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
        fcoef(:,ic) = 0.0
     else
        fcoef(1,ic) = 1.0/(k_perp(ic)**2*lambda_debye**2*dens_ele/temp_ele &
             + sum_den_h(it))
        if (n_field > 1) fcoef(2,ic) = 1.0/(-2.0*k_perp(ic)**2* &
             rho**2/betae_unit*dens_ele*temp_ele)
        if (n_field > 2) fcoef(3,ic) = -betae_unit/(2.0*dens_ele*temp_ele)
     endif
  enddo

  if (n_field > 1) then

     sum_loc(:) = 0.0
!$omp parallel private(iv_loc,ic)
!$omp do reduction(+:sum_loc)
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        do ic=1,nc
           it = it_c(ic)
           sum_loc(ic) = sum_loc(ic)+vfac(iv_loc)*dens_rot(it,is) &
                *jvec_c(2,ic,iv_loc)**2 
        enddo
     enddo
!$omp end do
!$omp end parallel

     call MPI_ALLREDUCE(sum_loc,&
          sum_cur_x,&
          size(sum_cur_x),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

  endif

  if (n_field == 1 .or. n_field == 2) then
     do ic=1,nc
        if (k_perp(ic) > 0.0) then
           gcoef(1,ic) = 1.0/(k_perp(ic)**2*lambda_debye**2*&
                dens_ele/temp_ele+sum_den_x(ic))
        endif
     enddo
  endif

  if (n_field > 1) then
     do ic=1,nc
        if (k_perp(ic) > 0.0) then
           gcoef(2,ic) = 1.0/(-2.0*k_perp(ic)**2*&
                rho**2/betae_unit*dens_ele*temp_ele-sum_cur_x(ic))
        endif
     enddo
  endif

  if (n_field > 2) then
     allocate(pb11(nc))
     allocate(pb12(nc))
     allocate(pb21(nc))
     allocate(pb22(nc))

     do ic=1,nc
        pb11(ic) = k_perp(ic)**2*lambda_debye**2* &
             dens_ele/temp_ele+sum_den_x(ic)
     enddo

     sum_loc(:)  = 0.0
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        do ic=1,nc
           ir = ir_c(ic) 
           it = it_c(ic)
           sum_loc(ic) = sum_loc(ic)-w_xi(ix)*w_e(ie)*dens(is) &
                *z(is)*jvec_c(1,ic,iv_loc)*jvec_c(3,ic,iv_loc) &
                *z(is)/temp(is)*dens_rot(it,is)
        enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          pb12,&
          size(pb12),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     pb21(:) = pb12(:)*betae_unit/(-2*dens_ele*temp_ele)

     sum_loc(:)  = 0.0
     iv_loc = 0
     do iv=nv1,nv2
        iv_loc = iv_loc+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        do ic=1,nc
           ir = ir_c(ic) 
           it = it_c(ic)
           sum_loc(ic) = sum_loc(ic)+w_xi(ix)*w_e(ie)*dens(is) &
                *temp(is)*jvec_c(3,ic,iv_loc)**2 &
                *(z(is)/temp(is))**2 * dens_rot(it,is)
        enddo
     enddo
     call MPI_ALLREDUCE(sum_loc,&
          pb22,&
          size(pb22),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     pb22(:) = 1.0-pb22(:)*betae_unit/(-2*dens_ele*temp_ele) 

     ! Determinant
     do ic=1,nc
        if (k_perp(ic) > 0.0) then
           sum_loc(ic) = pb11(ic)*pb22(ic)-pb12(ic)*pb21(ic)
        else
           sum_loc(ic) = 1.0
        endif
     enddo

     pb11 = pb11/sum_loc
     pb12 = pb12/sum_loc
     pb21 = pb21/sum_loc
     pb22 = pb22/sum_loc

     gcoef(3,:) = pb11
     gcoef(1,:) = pb22
     gcoef(4,:) = -pb12
     gcoef(5,:) = -pb21

     deallocate(pb11)
     deallocate(pb12)
     deallocate(pb21)
     deallocate(pb22)
  endif

  ! Set selected zeros
  do ic=1,nc
     ir = ir_c(ic) 
     if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
        gcoef(:,ic) = 0.0
     endif
  enddo

  ! Arrays to speed up velocity integrals in Maxwell equations
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ie = ie_v(iv)
     ix = ix_v(iv)
     do ic=1,nc
        it = it_c(ic)
        dvjvec_c(:,ic,iv_loc) = dens_rot(it,is)*w_e(ie)*w_xi(ix)*z(is)*dens(is)* &
             jvec_c(:,ic,iv_loc)
     enddo
  enddo
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
     it = it_c(ic)
     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        dvjvec_v(:,ic_loc,iv) = dens_rot(it,is)*w_e(ie)*w_xi(ix)*z(is)*dens(is)* &
             jvec_v(:,ic_loc,iv)
     enddo
  enddo
  !-------------------------------------------------------------------------

end subroutine cgyro_field_coefficients
