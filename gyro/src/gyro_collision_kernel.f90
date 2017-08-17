!---------------------------------------------------
! gyro_collision_kernel.f90
!
! PURPOSE:
!  Take a collision step using the pitch-angle 
!  scattering operator with the RBF method.  
!---------------------------------------------------

subroutine gyro_collision_kernel(ic)

  use gyro_globals
  use gyro_collision_private
  use gyro_pointers
  use math_constants

  !--------------------------------------------------------------
  implicit none
  !
  integer, intent(in) :: ic
  integer :: p,pp
  complex, dimension(n_rbf) :: fc,fcp,cphase
  !---------------------------------------------------------------

  p_ine_loc = 0

  do p_ine = 1+i_proc_1,n_ine_1,n_proc_1

     p_ine_loc = p_ine_loc+1

     i  = ine_i(p_ine)
     ie = ine_e(p_ine)

!$omp parallel do private(p,m) schedule(static)
     do k=1,n_lambda
        do m=1,n_stack
           p = m + (k - 1)*n_stack
           cphase(p) = exp(i_c*angp(i)*theta_t(i,k,m))
           fc(p) = h_c(m,k,p_ine_loc)*cphase(p)
        enddo
     enddo

     !----------------------------------------------- 
     ! This is essentially the full collision advance
     !
!$omp parallel do private(pp) schedule(static)
     do p=1,n_rbf
        fcp(p) = 0.0
        do pp=1,n_rbf
           fcp(p) = fcp(p)+d_rbf(pp,p,p_ine_loc,ic)*fc(pp)
        enddo
     enddo
     !----------------------------------------------- 

!$omp parallel do private(p,m) schedule(static)
     do k=1,n_lambda
        do m=1,n_stack
           p = m + (k - 1)*n_stack
           h_c(m,k,p_ine_loc) = fcp(p)/cphase(p)
        enddo
     enddo

  enddo ! p_ine_loc

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_collision_kernel done]'
  endif

end subroutine gyro_collision_kernel
