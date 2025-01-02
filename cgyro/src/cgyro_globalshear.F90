!-----------------------------------------------------------------
! cgyro_globalshear.F90
!
! PURPOSE:
!  Manage shearing by wavenumber advection (new Dudkovskaia terms)
!-----------------------------------------------------------------

subroutine cgyro_globalshear(ij)

  use cgyro_globals
  use timer_lib

  implicit none

  integer, intent(in) :: ij
  integer :: ir,l,ll,j,iccj,ivc,itor,llnt
  complex :: rl,h1,h2

  if (nonlinear_flag == 0 .or. source_flag == 0) return

  call timer_lib_in('shear')
  
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp private(ivc,ir,l,iccj,j,ll,rl,llnt,h1,h2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector &
!$acc private(ivc,ir,l,iccj,j,ll,rl,llnt,h1,h2) &
!$acc present(rhs(:,:,:,ij),omega_ss,omega_sbeta,field,cap_h_c,h_x,c_wave)
#else
!$omp parallel do collapse(3) private(ivc,ir,l,iccj,j,ll,rl,llnt,h1,h2) 
#endif
  do itor=nt1,nt2
     do ivc=1,nv_loc
        do ir=1,n_radial
           do j=1,n_theta

              iccj = (ir-1)*n_theta+j
              rl   = rhs(iccj,ivc,itor,ij)
              
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq
#endif
              do l=1,n_wave

                 ll = (2*l-1)
                 llnt = ll*n_theta

                 if ( (ir+ll) <= n_radial ) then
                    ! ExB shear
                    h1 = omega_eb_base*itor*h_x(iccj+llnt,ivc,itor)
                    ! beta_star shear
                    h1 = h1+omega_sbeta(iccj+llnt,ivc,itor)*( &
                         (1-sbeta_h)*cap_h_c(iccj+llnt,ivc,itor)+sbeta_h*h_x(iccj+llnt,ivc,itor))
                    ! omega_star shear
                    h1 = h1+sum(omega_ss(:,iccj+llnt,ivc,itor)*field(:,iccj+llnt,itor))
                 else
                    h1 = 0.0
                 endif
                
                 if ( (ir-ll) >= 1 ) then
                    ! ExB shear
                    h2 = omega_eb_base*itor*h_x(iccj-llnt,ivc,itor)
                    ! beta_star shear
                    h2 = h2+omega_sbeta(iccj-llnt,ivc,itor)*( &
                         (1-sbeta_h)*cap_h_c(iccj-llnt,ivc,itor)+sbeta_h*h_x(iccj-llnt,ivc,itor))
                    ! omega_star shear
                    h2 = h2+sum(omega_ss(:,iccj-llnt,ivc,itor)*field(:,iccj-llnt,itor))
                 else
                    h2 = 0.0
                 endif
                 
                 rl = rl+c_wave(l)*(h1-h2)

              enddo

              rhs(iccj,ivc,itor,ij) = rl

           enddo
        enddo
     enddo
  enddo

  call timer_lib_out('shear')

end subroutine cgyro_globalshear
