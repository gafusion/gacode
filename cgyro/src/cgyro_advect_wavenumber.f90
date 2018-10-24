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
  integer :: ir,l,ll,j,icc,in
  integer :: irp,irm,irpc,irmc
  complex, dimension(:,:),allocatable :: he
  complex :: dh
  real :: scale

  if (nonlinear_flag == 0) return
  
  if (profile_shear_flag == 1 .or. shear_method == 2) then
     call timer_lib_in('shear')
     allocate(he(n_theta,1-2*n_wave:n_radial+2*n_wave))

!$omp parallel do private(j,ir,in,icc,l,ll,he) &
!$omp&            private(scale,irm,irp,irmc,irpc,dh)
     do in=1,nv_loc
        he(:,1-2*n_wave:0) =0.0
        he(:,n_radial+1:n_radial+2*n_wave) =0.0

        ! Wavenumber advection ExB shear

        if (shear_method == 2) then

           do ir=1,n_radial
              icc = (ir-1)*n_theta
              do j=1,n_theta
                 he(j,ir) = omega_eb*h_x(icc+j,in)
              enddo
           enddo

           do ir=1,n_radial
              icc = (ir-1)*n_theta
              do l=1,n_wave
                 ll = 2*l-1
                 do j=1,n_theta
                    rhs(icc+j,in,ij) = rhs(icc+j,in,ij)+c_wave(l)*(he(j,ir+ll)-he(j,ir-ll))
                 enddo
              enddo
           enddo
        endif

        ! Wavenumber advection profile shear

        if (profile_shear_flag == 1) then

           do ir=1,n_radial
              icc = (ir-1)*n_theta
              do j=1,n_theta
                 he(j,ir) = sum(omega_ss(:,icc+j,in)*field(:,icc+j))
              enddo
           enddo

           do ir=1,n_radial
              icc = (ir-1)*n_theta
              do l=1,n_wave
                 ll = 2*l-1
                 do j=1,n_theta
                    rhs(icc+j,in,ij) = rhs(icc+j,in,ij)+c_wave(l)*(he(j,ir+ll)-he(j,ir-ll))
                 enddo
              enddo
           enddo
        endif

        ! Zonal damping (to reduce box-size correlation)

        if (n == 0) then
           scale = nu_global*abs(gamma_e)/2.0

           irm = -1+1+n_radial/2
           irp = irm+2
           irmc = (irm-1)*n_theta
           irpc = irmc + 2*n_theta

           do j=1,n_theta
              dh = scale*(h_x(irpc+j,in)-h_x(irmc+j,in))
              ! p = +1
              rhs(irpc+j,in,ij) = rhs(irpc+j,in,ij)-dh
              ! p = -1
              rhs(irmc+j,in,ij) = rhs(irmc+j,in,ij)+dh
           enddo
        endif

     enddo ! in=1,nv_loc

     deallocate(he)
     call timer_lib_out('shear')

  endif

end subroutine cgyro_advect_wavenumber
