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
  integer :: irp,irm
  complex, dimension(:,:),allocatable :: h0
  complex, dimension(nv_loc) :: dh
  real :: scale

  ! Wavenumber advection ExB shear
  if (shear_method == 2) then
     call timer_lib_in('shear')
     allocate(h0(nv_loc,1-2*n_wave:n_radial+2*n_wave))
     ! Zero wavenumber (probably causes boundary reflection)
!$omp parallel do private(j,ir,h0,dh,l,ll,ic)
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

     ! Zonal damping (to reduce box-size correlation)
     if (n == 0) then
        irp = +1+1+n_radial/2
        irm = -1+1+n_radial/2
        scale = nu_global*abs(gamma_e)/2.0
!$omp parallel do private(j,h0,dh,ic)
        do j=1,n_theta

           ! p = +1
           h0(:,irp) = h_x(ic_c(irp,j),:)

           ! p = -1
           h0(:,irm) = h_x(ic_c(irm,j),:)

           ! p = +1
           dh(:) = -scale*(h0(:,irp)-h0(:,irm))
           ic = ic_c(irp,j)
           rhs(ic,:,ij) = rhs(ic,:,ij)+dh(:)

           ! p = -1
           dh(:) = -scale*(h0(:,irm)-h0(:,irp))
           ic = ic_c(irm,j)
           rhs(ic,:,ij) = rhs(ic,:,ij)+dh(:)

        enddo
     endif

     deallocate(h0)
     call timer_lib_out('shear')
  endif

  !-------------------------------------------------------------------------

  ! Wavenumber advection profile shear
  if (profile_shear_flag == 1) then
     call timer_lib_in('shear')
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
     call timer_lib_out('shear')
  endif

end subroutine cgyro_advect_wavenumber
