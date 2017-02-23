!-----------------------------------------------------------------
! cgyro_field.f90
!
! PURPOSE:
!  Perform field solves for data distributed in both the velocity 
!  and configuration indices:
!  
!  - cgyro_field_v() (configuration distributed)
!  - cgyro_field_c() (velocity distributed)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Velocity (configuration-distributed) field solve
!-----------------------------------------------------------------
subroutine cgyro_field_v

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is,ie,ix,it
  complex :: fac

  
  call timer_lib_in('field_H')

  field_loc(:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of H

!$omp parallel do private(ic_loc,iv,is,ix,ie,fac,it)
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
     it = it_c(ic)
     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        fac = w_e(ie)*w_xi(ix)*z(is)*dens(is)*dens_rot(it,is)*cap_h_v(ic_loc,iv)
        field_loc(:,ic) = field_loc(:,ic)+fac*jvec_v(:,ic_loc,iv) 
     enddo
  enddo

  call MPI_ALLREDUCE(field_loc(:,:),&
       field(:,:),&
       size(field(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! Poisson LHS factors
  if (n == 0 .and. ae_flag == 1) then
     call cgyro_field_ae('v')
  else
     field(:,:) = fcoef(:,:)*field(:,:)
  endif

  call timer_lib_out('field_H')

end subroutine cgyro_field_v


!-----------------------------------------------------------------
! Configuration (velocity-distributed) field solve
!-----------------------------------------------------------------
subroutine cgyro_field_c

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is,it
  complex, dimension(nc) :: tmp
  
  call timer_lib_in('field_h')

  field_loc(:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

!$omp parallel private(ic)
!$omp do reduction(+:field_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        it = it_c(ic)
        field_loc(:,ic) = field_loc(:,ic)+(dvfac(iv_loc)*jvec_c(:,ic,iv_loc))&
             *dens_rot(it,is)*h_x(ic,iv_loc)
     enddo
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLREDUCE(field_loc(:,:),&
       field(:,:),&
       size(field(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (n_field > 2) then
     field(3,:) = field(3,:)*fcoef(3,:)
  endif

  ! Poisson LHS factors
 if (n == 0 .and. ae_flag == 1) then
    call cgyro_field_ae('c')
  else
     if (n_field > 2) then
!$omp workshare 
        tmp(:) = field(1,:)
        field(1,:) = gcoef(1,:)*field(1,:)+gcoef(4,:)*field(3,:)
        field(2,:) = gcoef(2,:)*field(2,:)
        field(3,:) = gcoef(3,:)*field(3,:)+gcoef(5,:)*tmp(:)
!$omp end workshare
     else
!$omp workshare
        field(:,:) = gcoef(:,:)*field(:,:)
!$omp end workshare
     endif
  endif

!$omp parallel do private(iv,iv_loc,ic,is)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        chi(ic,iv_loc) = sum(jxvec_c(:,ic,iv_loc)*field(:,ic))
        cap_h_c(ic,iv_loc) = h_x(ic,iv_loc)+psi(ic,iv_loc)*z(is)/temp(is)
     enddo
  enddo

  call timer_lib_out('field_h')

end subroutine cgyro_field_c


!-----------------------------------------------------------------
! Adiabatic electron field solves for n=0
!-----------------------------------------------------------------
subroutine cgyro_field_ae(space)

  use cgyro_globals
  
  implicit none

  character(len=1), intent(in) :: space
  integer :: ir,i,j
  real, dimension(n_theta) :: pvec_inr, pvec_ini

  if (space == 'c') then
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
           do i=1,n_theta
              pvec_outr(i) = 0
              pvec_outi(i) = 0
           enddo
           do i=1,n_theta
              pvec_inr(i) =  real(field(1,ic_c(ir,i)))
              pvec_ini(i) = aimag(field(1,ic_c(ir,i)))
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_outr(i) = pvec_outr(i)+xzf(ir,i,j)*pvec_inr(j)
                 pvec_outi(i) = pvec_outi(i)+xzf(ir,i,j)*pvec_ini(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i)) = cmplx(pvec_outr(i),pvec_outi(i))
           enddo
        endif
     enddo
  else
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
           do i=1,n_theta
              pvec_inr(i) =  real(field(1,ic_c(ir,i)))
              pvec_ini(i) = aimag(field(1,ic_c(ir,i)))
           enddo
           do i=1,n_theta
              pvec_outr(i) = 0
              pvec_outi(i) = 0
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_outr(i) = pvec_outr(i)+hzf(ir,i,j)*pvec_inr(j)
                 pvec_outi(i) = pvec_outi(i)+hzf(ir,i,j)*pvec_ini(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i)) = cmplx(pvec_outr(i),pvec_outi(i))
           enddo
        endif
     enddo
  endif

end subroutine cgyro_field_ae
