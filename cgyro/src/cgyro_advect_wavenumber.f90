!---------------------------------------------------------
! cgyro_advect_wavenumber.f90
!
! PURPOSE:
!  Manage shearing by wavenumber advection.
!---------------------------------------------------------

subroutine cgyro_advect_wavenumber(ij)

  use cgyro_globals
  use timer_lib

  implicit none

  integer, intent(in) :: ij
  integer :: ir,l,ll,j
  complex, dimension(:,:),allocatable :: h0
  complex, dimension(nv_loc) :: dh

  call timer_lib_in('shear')

  ! Wavenumber advection ExB shear
  if (shear_method == 2) then
     allocate(h0(nv_loc,1-2*n_wave:n_radial+2*n_wave))
     ! Zero wavenumber (probably causes boundary reflection)
!$omp parallel do private(j,ir,iv_loc,h0,dh,l,ll,ic)
     do j=1,n_theta
        h0 = 0.0
        do ir=1,n_radial
           h0(:,ir) = omega_eb*h_x(ic_c(ir,j),:)
        enddo
        do ir=1,n_radial
           dh(:) = 0.0
           do l=1,n_wave
              ll = 2*l-1
              dh(:) = dh(:)+c_wave(l)*(h0(:,ir+ll)-h0(:,ir-ll))
           enddo
           ic = ic_c(ir,j)
           rhs(ic,:,ij) = rhs(ic,:,ij)+dh(:)
        enddo
     enddo
     deallocate(h0)
  endif

  !-------------------------------------------------------------------------

  ! Wavenumber advection profile shear
  if (profile_shear_flag == 1) then
     allocate(h0(nv_loc,1-2*n_wave:n_radial+2*n_wave))
     ! Zero wavenumber (probably causes boundary reflection)
!$omp parallel do private(j,ir,iv_loc,h0,dh,l,ll,ic)
     do j=1,n_theta
        h0 = 0.0
        do ir=1,n_radial
           do iv_loc=1,nv_loc
              h0(iv_loc,ir) = sum(omega_ss(:,ic_c(ir,j),iv_loc)*field(:,ic_c(ir,j)))
           enddo
        enddo
        do ir=1,n_radial
           dh(:) = 0.0
           do l=1,n_wave
              ll = 2*l-1
              dh(:) = dh(:)+c_wave(l)*(h0(:,ir+ll)-h0(:,ir-ll))
           enddo
           ic = ic_c(ir,j)
           rhs(ic,:,ij) = rhs(ic,:,ij)+dh(:)
        enddo
     enddo
     deallocate(h0)
  endif

  call timer_lib_out('shear')

end subroutine cgyro_advect_wavenumber
