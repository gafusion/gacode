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

  
  call timer_lib_in('field')

  field_loc(:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of H

!$omp parallel do private(ic_loc,iv)
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
     do iv=1,nv
        field_loc(:,ic) = field_loc(:,ic)+dvjvec_v(:,ic_loc,iv)*cap_h_v(ic_loc,iv)
     enddo
  enddo

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call MPI_ALLREDUCE(field_loc(:,:),&
       field(:,:),&
       size(field(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

   call timer_lib_out('field_com')
  
  call timer_lib_in('field')

  ! Poisson LHS factors
  if (n == 0 .and. ae_flag == 1) then
     call cgyro_field_ae('v')
  else
     field(:,:) = fcoef(:,:)*field(:,:)
  endif

  call timer_lib_out('field')

end subroutine cgyro_field_v


!-----------------------------------------------------------------
! Configuration (velocity-distributed) field solve
!-----------------------------------------------------------------
subroutine cgyro_field_c

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is
  
  complex, dimension(nc) :: tmp
  
  call timer_lib_in('field')

  field_loc(:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

!$omp parallel private(iv_loc,ic)
!$omp do reduction(+:field_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        field_loc(:,ic) = field_loc(:,ic)+dvjvec_c(:,ic,iv_loc)*h_x(ic,iv_loc)
     enddo
  enddo
!$omp end do
!$omp end parallel

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call MPI_ALLREDUCE(field_loc(:,:),&
       field(:,:),&
       size(field(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  call timer_lib_out('field_com')

  if (n_field > 2) then
     field(3,:) = field(3,:)*fcoef(3,:)
  endif

  call timer_lib_in('field')

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

!$omp parallel do private(iv_loc,is,ic)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        psi(ic,iv_loc) = sum( jvec_c(:,ic,iv_loc)*field(:,ic))
        chi(ic,iv_loc) = sum(jxvec_c(:,ic,iv_loc)*field(:,ic))
        cap_h_c(ic,iv_loc) = h_x(ic,iv_loc)+psi(ic,iv_loc)*z(is)/temp(is)
     enddo
  enddo

  call timer_lib_out('field')

end subroutine cgyro_field_c

!-----------------------------------------------------------------
! Adiabatic electron field solves for n=0
!-----------------------------------------------------------------
subroutine cgyro_field_ae(space)

  use cgyro_globals
  use timer_lib  

  implicit none

  character(len=1), intent(in) :: space
  integer :: ir,i,j
  complex, dimension(n_theta) :: pvec_in,pvec_out

  call timer_lib_in('field')

  if (space == 'c') then
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
           do i=1,n_theta
              pvec_out(i) = 0.0
              pvec_in(i)  = field(1,ic_c(ir,i))
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_out(i) = pvec_out(i)+xzf(ir,i,j)*pvec_in(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i)) = pvec_out(i)
           enddo
        endif
     enddo
  else
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
           do i=1,n_theta
              pvec_out(i) = 0.0
              pvec_in(i)  =  field(1,ic_c(ir,i))
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_out(i) = pvec_out(i)+hzf(ir,i,j)*pvec_in(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i)) = pvec_out(i)
           enddo
        endif
     enddo
  endif

  call timer_lib_out('field')

end subroutine cgyro_field_ae
