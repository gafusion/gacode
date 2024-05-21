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
  integer :: ir,l,ll,j,iccj,ivc,itor,llnt
  complex :: rl,he1,he2

  if (nonlinear_flag == 0) return

  if (source_flag == 1) then
     call timer_lib_in('shear')

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(ivc,ir,l,iccj,j,ll,rl,llnt,he1,he2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector &
!$acc&         private(ivc,ir,l,iccj,j,ll,rl,llnt,he1,he2) &
!$acc&         present(rhs(:,:,:,ij),omega_ss,field,h_x,c_wave)
#else
!$omp parallel do collapse(4) private(ivc,ir,l,iccj,j,ll,rl,llnt,he1,he2)
#endif
     do itor=nt1,nt2
      do ivc=1,nv_loc
       do ir=1,n_radial
         do j=1,n_theta
           iccj = (ir-1)*n_theta+j

           ! Wavenumber advection ExB shear
           if (shear_method == 2) then
                 rl = 0.0
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq
#endif
                 do l=1,n_wave
                    ll = (2*l-1)
                    llnt = ll*n_theta
                    ! was he(j,ir+ll)
                    if ( (ir+ll) <= n_radial ) then
                       he1 = h_x(iccj+llnt,ivc,itor)
                    else
                       he1 = 0.0
                    endif
                    ! was he(j,ir-ll)
                    if ( (ir-ll) >= 1 ) then
                       he2 = h_x(iccj-llnt,ivc,itor)
                    else
                       he2 = 0.0
                    endif
                    ! Sign throughout paper is incorrect (or gamma -> - gamma)
                    ! Thus sign below has been checked and is correct
                    rl = rl+c_wave(l)*(he1-he2)
                 enddo
                 rhs(iccj,ivc,itor,ij) = rhs(iccj,ivc,itor,ij) + omega_eb_base*itor*rl
           endif

           ! Wavenumber advection profile shear
           if (profile_shear_flag == 1) then
                 iccj = (ir-1)*n_theta+j
                 rl = rhs(iccj,ivc,itor,ij)
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq
#endif
                 do l=1,n_wave
                    ll = 2*l-1
                    llnt = ll*n_theta
                    ! was he(j,ir+ll)
                    if ( (ir+ll) <= n_radial ) then
                       he1 = sum(omega_ss(:,iccj+llnt,ivc,itor)*field(:,iccj+llnt,itor))
                    else
                       he1 = 0.0
                    endif
                    ! was he(j,ir-ll)
                    if ( (ir-ll) >= 1 ) then
                       he2 = sum(omega_ss(:,iccj-llnt,ivc,itor)*field(:,iccj-llnt,itor))
                    else
                       he2 = 0.0
                    endif
                    ! Note opposite sign to ExB shear
                    rl = rl-c_wave(l)*(he1-he2)
                 enddo
                 rhs(iccj,ivc,itor,ij) = rl
           endif
         enddo
       enddo
      enddo
     enddo

     call timer_lib_out('shear')

  endif

end subroutine cgyro_advect_wavenumber
