!---------------------------------------------------------
! cgyro_advect_wavenumber.gpu.f90
!
! PURPOSE:
!  Manage shearing by wavenumber advection (GPU version).
!---------------------------------------------------------

subroutine cgyro_advect_wavenumber(ij)

  use cgyro_globals
  use timer_lib

  implicit none

  integer, intent(in) :: ij
  integer :: ir,l,ll,j,icc,in
  complex, dimension(:,:),allocatable :: he

  if (nonlinear_flag == 0) return

  if (source_flag == 1) then
     call timer_lib_in('shear')
     allocate(he(n_theta,1-2*n_wave:n_radial+2*n_wave))

!$acc parallel loop gang private(in,ir,l,icc,ll,he) &
!$acc&                   present(rhs(:,:,ij),omega_ss,field,h_x,c_wave) &
!$acc&                   vector_length(n_theta)
     do in=1,nv_loc
       he(:,1-2*n_wave:0) =0.0
       he(:,n_radial+1:n_radial+2*n_wave) =0.0

       ! Wavenumber advection ExB shear

       if (shear_method == 2) then

!$acc loop seq
         do ir=1,n_radial
           icc = (ir-1)*n_theta
!$acc loop vector private(j)
           do j=1,n_theta
             he(j,ir) = omega_eb*h_x(icc+j,in)
           enddo
         enddo

!$acc loop seq
         do ir=1,n_radial
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

!$acc loop seq 
         do ir=1,n_radial
           icc = (ir-1)*n_theta
!$acc loop vector private(j)
           do j=1,n_theta
              he(j,ir) = sum(omega_ss(:,icc+j,in)*field(:,icc+j))
           enddo
         enddo

!$acc loop seq
         do ir=1,n_radial
           icc = (ir-1)*n_theta
!$acc loop seq
           do l=1,n_wave
             ll = 2*l-1
!$acc loop vector private(j)
             do j=1,n_theta
                rhs(icc+j,in,ij) = rhs(icc+j,in,ij)-c_wave(l)*(he(j,ir+ll)-he(j,ir-ll))
              enddo
           enddo
         enddo
       endif

    enddo

    deallocate(he)
    call timer_lib_out('shear')

  endif 

end subroutine cgyro_advect_wavenumber
