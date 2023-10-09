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

  integer :: ir,itor
  integer :: n_changes
  integer, dimension(nt1:nt2) :: change_sign
  integer :: my_sign

  
  n_changes = 0

  call timer_lib_in('shear')
  do itor=nt1,nt2
   gtime(itor) = gtime(itor)+omega_eb_base*itor*delta_t

   ! Forward shearing
   if (gtime(itor) > 0.5) then

     gtime(itor) = gtime(itor)-1.0
     n_changes = n_changes + 1
     change_sign(itor) = 1

   ! Backward shearing
   else if (gtime(itor) < -0.5) then

     gtime(itor) = gtime(itor)+1.0
     n_changes = n_changes + 1
     change_sign(itor) = -1

   else
     change_sign(itor) = 0
   endif
  enddo

  if (n_changes > 0) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(2) &
!$omp&  private(ir,my_sign) map(to:change_sign)
#elif defined(_OPENACC)
!$acc parallel loop collapse(2) independent gang vector_length(n_theta) &
!$acc&              private(ir,my_sign) present(h_x,ic_c) copyin(change_sign)
#else
!$omp parallel do collapse(2) private(ir,my_sign)
#endif
    do itor=nt1,nt2
     do iv_loc=1,nv_loc

      my_sign = change_sign(itor)
      if (my_sign == 1) then
       ! Forward shearing
#if (!defined(OMPGPU)) && defined(_OPENACC)
!acc loop seq
#endif
       do ir=2,n_radial
         h_x(ic_c(ir-1,:),iv_loc,itor) = h_x(ic_c(ir,:),iv_loc,itor)
       enddo

       h_x(ic_c(n_radial,:),iv_loc,itor) = 0.0
      else if (my_sign == -1) then
       ! Backward shearing
#if (!defined(OMPGPU)) && defined(_OPENACC)
!acc loop seq
#endif
       do ir=n_radial-1,1,-1
         h_x(ic_c(ir+1,:),iv_loc,itor) = h_x(ic_c(ir,:),iv_loc,itor)
       enddo

       h_x(ic_c(1,:),iv_loc,itor) = 0.0
      endif
     enddo
    enddo
  endif

  call timer_lib_out('shear')

  if (n_changes > 0) then
     call cgyro_field_c
  endif

end subroutine cgyro_shear_hammett
