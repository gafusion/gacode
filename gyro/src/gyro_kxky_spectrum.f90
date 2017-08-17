!-------------------------------------------------------
! get_field_spectrum.f90
!
! PURPOSE:
!  This routine calculates the instantaneous potential
!  spectrum, averaged over theta.
!
! NOTES:
!
! Let <> = theta average
!
! kxkyspec = < |phi_np|^2 >
!
!             L
!          1  /
!  NOTE:  --- | dr [<|phi_n|^2>] = sum_p < |phi_np|^2 >
!          L  /
!             0
!
! Also note that 
!
! sum_n sum_p |phi_np|^2 = Volume average of phi^2
!-------------------------------------------------------

subroutine gyro_kxky_spectrum

  use gyro_globals
  use math_constants
  use ompdata

  !-------------------------------------
  implicit none
  !
  integer :: p
  real :: temp
  complex :: f_bar(n_theta_int,n_x)
  !-------------------------------------

  f_bar = (0.0,0.0)

  ! Compute FT of phi.
!$omp parallel private(i,p,j)
  do i=1,n_x  
     do p = ibeg, iend
        do j=1,n_theta_int
           f_bar(j,p) = f_bar(j,p)+phi(j,i,1)*cRi(p,i)
        enddo ! j
     enddo ! p
  enddo ! i
!$omp end parallel

  ! Find spectral sum:

  ! Previous version of GYRO eliminated certain spectral 
  ! components.  Now, we just do the full (periodic) spectral 
  ! decomposition to avoid confusion.

  do p=1,n_x

     temp = 0.0
     do j=1,n_theta_int
        temp = temp+abs(f_bar(j,p))**2
     enddo
     kxkyspec(p) = temp/n_theta_int

  enddo ! p

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_field_spectrum done]' 
  endif

end subroutine gyro_kxky_spectrum
