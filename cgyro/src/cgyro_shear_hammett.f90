!-------------------------------------------------------------------
! cgyro_shear_hammett.f90
!
! PURPOSE:
!  Spectral shear algorithm using the Hammett discrete shift method.
!  Wavenumbers are shifted to the left or right at selected time 
!  intervals depending upon sign of omega_eb.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!--------------------------------------------------------------------

subroutine cgyro_shear_hammett

  use cgyro_globals
  use timer_lib

  implicit none

  integer :: ir
  complex, dimension(n_theta,nv_loc) :: a1

  
  gtime = gtime+omega_eb*delta_t

  ! Forward shearing
  if (gtime > 0.5) then

     call timer_lib_in('shear')
     gtime = gtime-1.0

     a1 = h_x(ic_c(1,:),:)

     do ir=2,n_radial
        h_x(ic_c(ir-1,:),:) = h_x(ic_c(ir,:),:)
     enddo

     h_x(ic_c(n_radial,:),:) = 0.0
     call timer_lib_out('shear')

     call cgyro_field_c

  endif

  ! Backward shearing
  if (gtime < -0.5) then

     call timer_lib_in('shear')
     gtime = gtime+1.0

     a1 = h_x(ic_c(n_radial,:),:) 

     do ir=n_radial-1,1,-1
        h_x(ic_c(ir+1,:),:) = h_x(ic_c(ir,:),:)
     enddo

     h_x(ic_c(1,:),:) = 0.0
     call timer_lib_out('shear')

     call cgyro_field_c

  endif

end subroutine cgyro_shear_hammett
