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

  if (profile_shear_flag == 1 .or. shear_method == 2) then
     call timer_lib_in('shear')
     allocate(he(n_theta,1-2*n_wave:n_radial+2*n_wave))

!$acc parallel loop gang private(in,ir,l,icc,ll,he) &
!$acc&                   present(rhs(:,:,ij),omega_ss,field,h_x,c_wave) &
!$acc&                   private(scale,irm,irp,irmc,irpc) &
!$acc&                   vector_length(n_theta)
     do in=1,nv_loc
       he(:,1-2*n_wave:0) =0.0
       he(:,n_radial+1:n_radial+2*n_wave) =0.0

       ! Wavenumber advection ExB shear

       if (shear_method == 2) then
         !ic_c(ir,j) = j + (ir-1)*n_theta

!$acc loop seq
         do ir=1,n_radial
           !icc =ic_c(ir,1)-1
           icc = (ir-1)*n_theta
!$acc loop vector private(j)
           do j=1,n_theta
             he(j,ir) = omega_eb*h_x(icc+j,in)
           enddo
         enddo

!$acc loop seq
         do ir=1,n_radial
           !icc =ic_c(ir,1)-1
           icc = (ir-1)*n_theta
!$acc loop seq
           do l=1,n_wave
              ll = 2*l-1
!$acc loop vector private(j)
              do j=1,n_theta
                rhs(icc+j,in,ij) = rhs(icc+j,in,ij)+c_wave(l)*(he(j,ir+ll)-he(j,ir-ll))
              enddo
           enddo
         enddo
       endif

       ! Wavenumber advection profile shear

       if (profile_shear_flag == 1) then
         !ic_c(ir,j) = j + (ir-1)*n_theta
 
         do ir=1,n_radial
           !icc =ic_c(ir,1)-1
           icc = (ir-1)*n_theta
!$acc loop vector private(j)
           do j=1,n_theta
              he(j,ir) = sum(omega_ss(:,icc+j,in)*field(:,icc+j))
           enddo
         enddo

!$acc loop seq
         do ir=1,n_radial
           !icc =ic_c(ir,1)-1
           icc = (ir-1)*n_theta
!$acc loop seq
           do l=1,n_wave
             ll = 2*l-1
!$acc loop vector private(j)
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
         !irmc = ic_c(irm,j)-1
         !irpc = ic_c(irp,j)-1
         irmc = (irm-1)*n_theta
         irpc = irmc + 2*n_theta

!$acc loop vector private(j,dh)
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
