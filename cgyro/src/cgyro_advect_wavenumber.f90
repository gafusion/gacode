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
  integer :: ir,ip,j
  complex, dimension(:),allocatable :: h0
  complex :: dh,uh
  complex, dimension(n_radial-1) :: al,au
  complex, dimension(n_radial) :: ad
  complex, dimension(n_radial,n_theta) :: b

  call timer_lib_in('shear')

  ! Wavenumber advection ExB shear
  if (shear_method == 2) then
     allocate(h0(1-nup_wave:n_radial+nup_wave))
!$omp parallel do private(j,h0,ir,dh,uh,ip,ic)
     do iv_loc=1,nv_loc
        do j=1,n_theta
           h0(1:n_radial) = h_x(ic_c(:,j),iv_loc)
           h0(1-nup_wave:0) = h0(1)
           h0(n_radial+1) = h0(n_radial)
           do ir=1,n_radial
              dh = 0.0
              uh = 0.0
              do ip=-nup_wave,nup_wave
                 dh = dh+der_wave(ip)*h0(ir+ip)
                 uh = uh+dis_wave(ip)*h0(ir+ip)
              enddo
              ic = ic_c(ir,j)
              rhs(ic,iv_loc,ij) = rhs(ic,iv_loc,ij)+&
                   omega_eb*dh-abs(omega_eb)*uh*up_wave
           enddo
        enddo
     enddo
     deallocate(h0)
  endif

  if (shear_method == 5) then
     allocate(h0(-1:n_radial+2))
     do iv_loc=1,nv_loc

        h0 = 0.0
        select case (nup_wave)
        case (1)

           ! Lele 4th order compact
           do j=1,n_theta
              h0(1:n_radial) = h_x(ic_c(:,j),iv_loc)
              do ir=1,n_radial
                 b(ir,j) = 11.0*(h0(ir+1)-h0(ir-1))+0.5*(h0(ir+2)-h0(ir-2))
              enddo
           enddo
           al(:) = 5.0
           ad(:) = 14.0
           au(:) = 5.0
       
        case (2)
           
           ! 4th order compact
           do j=1,n_theta
              h0(1:n_radial) = h_x(ic_c(:,j),iv_loc)
              do ir=1,n_radial
                 b(ir,j) = 3*(h0(ir+1)-h0(ir-1)) 
              enddo
           enddo
           al(:) = 1.0
           ad(:) = 4.0
           au(:) = 1.0

        case (3)
         ! 6th order compact
          do j=1,n_theta
              h0(1:n_radial) = h_x(ic_c(:,j),iv_loc)
              do ir=1,n_radial
                 b(ir,j) = 7/3.0*(h0(ir+1)-h0(ir-1))+1/12.0*(h0(ir+2)-h0(ir-2))
              enddo
           enddo
           al(:) = 1.0
           ad(:) = 3.0
           au(:) = 1.0
        end select

        call ZGTSV(n_radial,n_theta,al,ad,au,b,n_radial,info)
        do j=1,n_theta
           rhs(ic_c(:,j),iv_loc,ij) = rhs(ic_c(:,j),iv_loc,ij)+omega_eb*b(:,j)
        enddo
     enddo
     deallocate(h0)
  endif

  !-------------------------------------------------------------------------

  ! Wavenumber advection profile shear
  if (profile_shear_flag == 1) then
!$omp parallel do private(j,h0,ir,dh,ip)
     do iv_loc=1,nv_loc
        h0 = 0.0
        do j=1,n_theta
           do ir=1,n_radial
              h0(ir) = sum(omega_ss(:,ic_c(ir,j),iv_loc)*field(:,ic_c(ir,j)))
           enddo
           do ir=1,n_radial
              dh = 0.0
              do ip=-nup_wave,nup_wave
                 dh = dh+der_wave(ip)*h0(ir+ip)
              enddo
              ic = ic_c(ir,j)
              rhs(ic,iv_loc,ij) = rhs(ic,iv_loc,ij)+dh
           enddo
        enddo
     enddo
  endif

  call timer_lib_out('shear')

end subroutine cgyro_advect_wavenumber
