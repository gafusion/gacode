subroutine cgyro_field_coefficients

  use mpi
  use cgyro_globals

  implicit none

  integer :: ir,it,is,ie,ix,itor
  real :: sum_one
  real, dimension(:,:), allocatable :: sum_loc
  real, dimension(:,:), allocatable :: pb11,pb12,pb21,pb22
 
  !-------------------------------------------------------------------------
  ! Field equation prefactors, sums.
  !
  allocate(sum_loc(nc,nt1:nt2))
!$omp parallel do collapse(2) private(iv,iv_loc,is,it,sum_one) shared(sum_loc)
  do itor=nt1,nt2
    do ic=1,nc
      it = it_c(ic)
      sum_one = 0.0
      do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        sum_one = sum_one + vfac(iv_loc)*dens_rot(it,is) &
             *(1.0-jvec_c(1,ic,iv_loc,itor)**2) 
      enddo
      sum_loc(ic,itor) = sum_one
    enddo
  enddo

  call MPI_ALLREDUCE(sum_loc,&
       sum_den_x,&
       size(sum_den_x),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (ae_flag == 1) then
    do itor=nt1,nt2
     do ic=1,nc
        it = it_c(ic)
        sum_den_x(ic,itor) = sum_den_x(ic,itor)+dens_ele*dens_ele_rot(it)/temp_ele
     enddo
    enddo
  endif
  !----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Field-solve coefficients (i.e., final numerical factors).
  !
  do itor=nt1,nt2
   do ic=1,nc
     ir = ir_c(ic) 
     it = it_c(ic)
     if (itor == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
        fcoef(:,ic,itor) = 0.0
     else
        fcoef(1,ic,itor) = 1.0/(k_perp(ic,itor)**2*lambda_debye**2*dens_ele/temp_ele &
             + sum_den_h(it))
        if (n_field > 1) fcoef(2,ic,itor) = 1.0/(-2.0*k_perp(ic,itor)**2* &
             rho**2/betae_unit*dens_ele*temp_ele)
        if (n_field > 2) fcoef(3,ic,itor) = -betae_unit/(2.0*dens_ele*temp_ele)
     endif
   enddo
  enddo

  if (n_field > 1) then

!$omp parallel do collapse(2) private(iv,iv_loc,is,it,sum_one) shared(sum_loc)
     do itor=nt1,nt2
       do ic=1,nc
         it = it_c(ic)
         sum_one = 0.0
         do iv=nv1,nv2
           iv_loc = iv-nv1+1
           is = is_v(iv)
           sum_one = sum_one + vfac(iv_loc)*dens_rot(it,is) &
                *jvec_c(2,ic,iv_loc,itor)**2 
         enddo
         sum_loc(ic,itor) = sum_one
       enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          sum_cur_x,&
          size(sum_cur_x),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

  endif

  if (n_field == 1 .or. n_field == 2) then
    do itor=nt1,nt2
     do ic=1,nc
        if (k_perp(ic,itor) > 0.0) then
           gcoef(1,ic,itor) = 1.0/(k_perp(ic,itor)**2*lambda_debye**2*&
                dens_ele/temp_ele+sum_den_x(ic,itor))
        endif
     enddo
    enddo
  endif

  if (n_field > 1) then
    do itor=nt1,nt2
     do ic=1,nc
        if (k_perp(ic,itor) > 0.0) then
           gcoef(2,ic,itor) = 1.0/(-2.0*k_perp(ic,itor)**2*&
                rho**2/betae_unit*dens_ele*temp_ele-sum_cur_x(ic,itor))
        endif
     enddo
    enddo
  endif

  if (n_field > 2) then
     allocate(pb11(nc,nt1:nt2))
     allocate(pb12(nc,nt1:nt2))
     allocate(pb21(nc,nt1:nt2))
     allocate(pb22(nc,nt1:nt2))

     do itor=nt1,nt2
      do ic=1,nc
        pb11(ic,itor) = k_perp(ic,itor)**2*lambda_debye**2* &
             dens_ele/temp_ele+sum_den_x(ic,itor)
      enddo
     enddo

!$omp parallel do collapse(2) shared(sum_loc) &
!$omp&         private(iv,iv_loc,ir,it,sum_one,is,ix,ie)
     do itor=nt1,nt2
       do ic=1,nc
         ir = ir_c(ic) 
         it = it_c(ic)
         sum_one = 0.0
         do iv=nv1,nv2
           iv_loc = iv-nv1+1
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           sum_one = sum_one - w_exi(ie,ix)*dens2_rot(it,is) &
                *z(is)*jvec_c(1,ic,iv_loc,itor)*jvec_c(3,ic,iv_loc,itor) &
                *z(is)/temp(is)
         enddo
         sum_loc(ic,itor) = sum_one
       enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          pb12,&
          size(pb12),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

!$omp parallel do collapse(2) shared(sum_loc) &
!$omp&         private(iv,iv_loc,ir,it,sum_one,is,ix,ie)
     do itor=nt1,nt2
       do ic=1,nc
         ir = ir_c(ic) 
         it = it_c(ic)
         sum_one = 0.0
         do iv=nv1,nv2
           iv_loc = iv-nv1+1
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           sum_one = sum_one + w_exi(ie,ix)*dens2_rot(it,is) &
                *temp(is)*jvec_c(3,ic,iv_loc,itor)**2 &
                *(z(is)/temp(is))**2
         enddo
         sum_loc(ic,itor) = sum_one
       enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          pb22,&
          size(pb22),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     ! Determinant
!$omp parallel do collapse(2) shared(sum_loc) private(sum_one)
     do itor=nt1,nt2
       do ic=1,nc
         pb21(ic,itor) = pb12(ic,itor)*betae_unit/(-2*dens_ele*temp_ele)
         pb22(ic,itor) = 1.0-pb22(ic,itor)*betae_unit/(-2*dens_ele*temp_ele) 

         if (k_perp(ic,itor) > 0.0) then
           sum_one = pb11(ic,itor)*pb22(ic,itor)-pb12(ic,itor)*pb21(ic,itor)
         else
           sum_one = 1.0
         endif

         gcoef(3,ic,itor) = pb11(ic,itor)/sum_one
         gcoef(1,ic,itor) = pb22(ic,itor)/sum_one
         gcoef(4,ic,itor) = -pb12(ic,itor)/sum_one
         gcoef(5,ic,itor) = -pb21(ic,itor)/sum_one
       enddo
     enddo

     deallocate(pb11)
     deallocate(pb12)
     deallocate(pb21)
     deallocate(pb22)
  endif

  deallocate(sum_loc)

!$omp parallel do private(iv,iv_loc,is,ie,ix,ic,it,ic_loc,ir) shared(gcoef,dvjvec_c,dvjvec_v)
  do itor=nt1,nt2
   ! Set selected zeros
   do ic=1,nc
     ir = ir_c(ic) 
     if (itor == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
        gcoef(:,ic,itor) = 0.0
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
        dvjvec_c(:,ic,iv_loc,itor) = dens2_rot(it,is)*w_exi(ie,ix)*z(is)* &
             jvec_c(:,ic,iv_loc,itor)
     enddo
   enddo
   do ic=nc1,nc2
     ic_loc = ic-nc1+1
     it = it_c(ic)
     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        dvjvec_v(:,ic_loc,itor,iv) = dens2_rot(it,is)*w_exi(ie,ix)*z(is)* &
             jvec_v(:,ic_loc,itor,iv)
     enddo
   enddo
  enddo
  !-------------------------------------------------------------------------

#if defined(OMPGPU)
!$omp target update to(fcoef,gcoef,dvjvec_c,dvjvec_v)
#elif defined(_OPENACC)
!$acc update device(fcoef,gcoef,dvjvec_c,dvjvec_v)
#endif

end subroutine cgyro_field_coefficients
