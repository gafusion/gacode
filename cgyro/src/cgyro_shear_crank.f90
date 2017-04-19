!-------------------------------------------------------------------
! cgyro_shear_crank.f90
!
! PURPOSE:
!  Crank-Nicholson (trapezoidal) advection solve.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!--------------------------------------------------------------------

subroutine cgyro_shear_crank

  use cgyro_globals

  implicit none

  integer :: ir,j
  real :: s0
  complex, dimension(n_radial-1) :: al,au
  complex, dimension(n_radial) :: ad,h0,b

  s0 = omega_eb*delta_t/2

  do iv_loc=1,nv_loc
     do j=1,n_theta
        h0(:) = h_x(ic_c(:,j),iv_loc)
        b(1) = h0(1)-s0*h0(2)
        do ir=2,n_radial-1
           b(ir) = h0(ir)+s0*(h0(ir-1)-h0(ir+1)) 
        enddo
        b(n_radial) = h0(n_radial)+s0*h0(n_radial-1)
        al(:) = -s0
        ad(:) = 1.0
        au(:) = s0
        call ZGTSV(n_radial,1,al,ad,au,b,n_radial,info)
        h_x(ic_c(:,j),iv_loc) = b
     enddo
  enddo

  call cgyro_field_c

end subroutine cgyro_shear_crank
